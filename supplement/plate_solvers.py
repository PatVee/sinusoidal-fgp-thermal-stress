"""Local linear thermoelastic benchmark: full 3D Q2 solids and sinusoidal plates.

Units internally: length / a, stress / E0, E0 = 70 GPa. Temperature in K.
This is NOT an implementation of Eringen nonlocality or a full Maxwell solver.
The optional magnetic operator is an explicitly idealised transverse restoring
energy (cB/2)*integral((dw/dx)^2 dV); its dimensionless strength is etaB.
"""
from dataclasses import dataclass, asdict
from functools import lru_cache
import time
import numpy as np
from numpy.polynomial.legendre import leggauss, Legendre
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve

E0 = 70e9
ALPHA0 = 23.4e-6


@dataclass(frozen=True)
class Case:
    slenderness: float = 10.0
    aspect: float = 1.0  # b/a
    gradation: float = 1.0
    porosity: float = 0.1  # mean true void fraction (assumed effective laws)
    pattern: str = 'uniform'
    deltaT: float = 20.0
    bottomT: float = 0.0  # rise above stress-free reference
    pressure: float = 0.0  # p/E0; positive in +z
    support: str = 'SS'
    etaB: float = 0.0
    E_exponent: float = 2.0
    k_exponent: float = 1.0
    E_scale: float = 1.0
    alpha_scale: float = 1.0
    nonlocal_length: float = 0.0
    material: str = 'FGM'

    @property
    def h(self): return 1.0 / self.slenderness

    @property
    def cB(self): return self.etaB * self.h**2 / (12*(1-0.3**2))

    def validate(self):
        if self.nonlocal_length != 0:
            raise NotImplementedError('Finite nonlocal length is not implemented; no stiffness rescaling is used.')
        if self.support not in ['SS','CCCC']:
            raise ValueError('Supports: SS diaphragm or CCCC only.')
        if not (0 <= self.porosity < 0.45): raise ValueError('Mean porosity outside benchmark domain.')
        if self.gradation < 0: raise ValueError('Negative gradation exponent.')


def props(z, c):
    z=np.asarray(z); t=np.clip(z/c.h+0.5,0,1)
    v=t**c.gradation
    if c.material=='Al': v=np.zeros_like(t)
    if c.material=='SiC': v=np.ones_like(t)
    if c.pattern=='uniform': phi=np.zeros_like(t)+c.porosity
    elif c.pattern=='centre': phi=c.porosity*(1+np.cos(2*np.pi*(t-.5)))
    elif c.pattern=='faces': phi=c.porosity*(1-np.cos(2*np.pi*(t-.5)))
    else: raise ValueError(c.pattern)
    E=(1+(427/70-1)*v)*(1-phi)**c.E_exponent*c.E_scale
    nu=.3+(.17-.3)*v
    al=(23.4e-6+(4.3e-6-23.4e-6)*v)*c.alpha_scale
    k=(233+(65-233)*v)*(1-phi)**c.k_exponent
    return E,nu,al,k,phi


@lru_cache(maxsize=512)
def thermal(c):
    z=np.linspace(-c.h/2,c.h/2,4097)
    k=props(z,c)[3]
    resistance=cumulative_trapezoid(1/k,z,initial=0)
    T=c.bottomT+c.deltaT*resistance/resistance[-1]
    return CubicSpline(z,T), -c.deltaT/resistance[-1]


def C3(E,nu):
    E,nu=np.broadcast_arrays(E,nu)
    lam=E*nu/((1+nu)*(1-2*nu)); mu=E/(2*(1+nu))
    C=np.zeros(E.shape+(6,6))
    C[...,:3,:3]=lam[...,None,None]
    for i in range(3): C[...,i,i]+=2*mu
    for i in range(3,6): C[...,i,i]=mu
    return C


def Cplate(E,nu):
    # strain order xx, yy, xy (engineering), xz, yz
    E,nu=np.broadcast_arrays(E,nu); C=np.zeros(E.shape+(5,5))
    Q=E/(1-nu**2)
    C[...,0,0]=C[...,1,1]=Q
    C[...,0,1]=C[...,1,0]=nu*Q
    for j in [2,3,4]: C[...,j,j]=E/(2*(1+nu))
    return C


def lagrange2(x):
    return np.array([x*(x-1)/2,1-x*x,x*(x+1)/2]), np.array([x-.5,-2*x,x+.5])


@lru_cache(maxsize=128)
def solid_mesh(nx,ny,nz,b,h):
    gx=np.linspace(0,1,2*nx+1); gy=np.linspace(0,b,2*ny+1); gz=np.linspace(-h/2,h/2,2*nz+1)
    nodes=np.stack(np.meshgrid(gx,gy,gz,indexing='ij'),axis=-1).reshape(-1,3)
    ids=np.arange(len(nodes)).reshape(len(gx),len(gy),len(gz))
    elems=np.array([ids[2*i:2*i+3,2*j:2*j+3,2*k:2*k+3].ravel()
                    for i in range(nx) for j in range(ny) for k in range(nz)])
    dofs=(3*elems[:,:,None]+np.arange(3)).reshape(len(elems),81)
    return nodes,elems,dofs


def shape3(r,s,t,dx,dy,dz):
    Lx,Dx=lagrange2(r); Ly,Dy=lagrange2(s); Lz,Dz=lagrange2(t)
    N=np.einsum('i,j,k->ijk',Lx,Ly,Lz).ravel()
    G=np.stack([np.einsum('i,j,k->ijk',Dx,Ly,Lz).ravel()*2/dx,
                np.einsum('i,j,k->ijk',Lx,Dy,Lz).ravel()*2/dy,
                np.einsum('i,j,k->ijk',Lx,Ly,Dz).ravel()*2/dz],axis=1)
    B=np.zeros((6,81))
    B[0,0::3]=G[:,0]; B[1,1::3]=G[:,1]; B[2,2::3]=G[:,2]
    B[3,0::3]=G[:,1]; B[3,1::3]=G[:,0]
    B[4,0::3]=G[:,2]; B[4,2::3]=G[:,0]
    B[5,1::3]=G[:,2]; B[5,2::3]=G[:,1]
    return N,G,B


def solid(c,nx=8,ny=None,nz=4,patch=None):
    """27-node tensor-product quadratic solid elements with full integration."""
    c.validate(); start=time.perf_counter(); ny=ny or nx
    nodes,elems,dofs=solid_mesh(nx,ny,nz,c.aspect,c.h)
    ne=len(elems); dx=1/nx; dy=c.aspect/ny; dz=c.h/nz
    centres=nodes[elems[:,13]]
    Ke=np.zeros((ne,81,81)); Fe=np.zeros((ne,81)); g,w=leggauss(3)
    Tfun,_=thermal(c)
    # Matrices only vary through thickness; assemble nz prototype elements.
    protoK=np.zeros((nz,81,81)); protoF=np.zeros((nz,81));
    zc=np.linspace(-c.h/2+dz/2,c.h/2-dz/2,nz)
    for i,r in enumerate(g):
      for j,s in enumerate(g):
       for k,t in enumerate(g):
        N,G,B=shape3(r,s,t,dx,dy,dz)
        z=zc+t*dz/2; E,nu,al,_,_=props(z,c); C=C3(E,nu)
        fac=w[i]*w[j]*w[k]*dx*dy*dz/8
        et=np.zeros((nz,6)); et[:,:3]=(al*Tfun(z))[:,None]
        protoK+=np.einsum('ia,kij,jb->kab',B,C,B,optimize=True)*fac
        protoF+=np.einsum('ia,kij,kj->ka',B,C,et,optimize=True)*fac
        if c.cB:
            bm=np.zeros(81); bm[2::3]=G[:,0]
            protoK+=c.cB*np.outer(bm,bm)[None,:,:]*fac
    ez=np.tile(np.arange(nz),nx*ny); Ke[:]=protoK[ez]; Fe[:]=protoF[ez]
    if c.pressure:
      top=np.where(ez==nz-1)[0]
      for i,r in enumerate(g):
       for j,s in enumerate(g):
        N,_,_=shape3(r,s,1,dx,dy,dz)
        Fe[top,2::3]+=c.pressure*N[None,:]*w[i]*w[j]*dx*dy/4
    rows=np.repeat(dofs,81,axis=1).ravel(); cols=np.tile(dofs,(1,81)).ravel()
    K=coo_matrix((Ke.ravel(),(rows,cols)),shape=(3*len(nodes),)*2).tocsr()
    F=np.zeros(3*len(nodes)); np.add.at(F,dofs.ravel(),Fe.ravel())
    x,y,z=nodes.T; onx=np.isclose(x,0)|np.isclose(x,1); ony=np.isclose(y,0)|np.isclose(y,c.aspect)
    fixed=np.zeros((len(nodes),3),dtype=bool)
    prescribed=np.zeros((len(nodes),3))
    if patch is not None:
        boundary=onx|ony|np.isclose(z,-c.h/2)|np.isclose(z,c.h/2)
        fixed[boundary,:]=True
        prescribed[:]=patch(nodes)
    elif c.support=='SS':
        fixed[onx,1:3]=True; fixed[ony,0]=True; fixed[ony,2]=True
    else: fixed[onx|ony,:]=True
    fixed=fixed.ravel(); free=np.where(~fixed)[0]
    u=prescribed.ravel().copy()
    rhs=F[free]-K[free][:,fixed]@u[fixed]
    u[free]=spsolve(K[free][:,free],rhs)
    residual=np.linalg.norm((K@u-F)[free])/max(np.linalg.norm(rhs),1e-30)
    return dict(case=c,kind='Q2_3D',u=u,nodes=nodes,elems=elems,dofs=dofs,
                mesh=(nx,ny,nz),seconds=time.perf_counter()-start,residual=residual,
                energy=float(.5*u@(K@u)),ndof=len(free))


def eval_solid(result,xyz):
    c=result['case']; nx,ny,nz=result['mesh']; dx=1/nx; dy=c.aspect/ny; dz=c.h/nz
    xyz=np.asarray(xyz); values=[]; stresses=[]
    for x,y,z in xyz:
        ix=min(nx-1,max(0,int(x/dx))); iy=min(ny-1,max(0,int(y/dy))); iz=min(nz-1,max(0,int((z+c.h/2)/dz)))
        r=2*(x-ix*dx)/dx-1; s=2*(y-iy*dy)/dy-1; t=2*(z+c.h/2-iz*dz)/dz-1
        N,G,B=shape3(r,s,t,dx,dy,dz)
        e=(ix*ny+iy)*nz+iz; ue=result['u'][result['dofs'][e]]
        E,nu,al,_,_=props(z,c); eps=B@ue; eps[:3]-=al*thermal(c)[0](z)
        values.append(N@ue.reshape(27,3)); stresses.append(C3(E,nu)@eps)
    return np.array(values),np.array(stresses)


def basis_1d(x,axis,component,support,p):
    """Symmetry-adapted Ritz basis, derivatives through second order."""
    x=np.asarray(x)
    odd=(component in ['u','rx'] and axis=='x') or (component in ['v','ry'] and axis=='y')
    vals=[]; d1=[]; d2=[]
    for i in range(p):
        if support=='SS':
            n=(2*i+1)*np.pi
            if odd: v=np.cos(n*x); d=-n*np.sin(n*x); dd=-n*n*v
            else: v=np.sin(n*x); d=n*np.cos(n*x); dd=-n*n*v
            scale=1.0
        else:
            # Legendre polynomials on [0,1], with essential boundary bubbles.
            degree=2*i+(1 if odd else 0)
            P=Legendre.basis(degree).convert(kind=np.polynomial.Polynomial)
            s=np.polynomial.Polynomial([-1,2]); poly=P(s)
            bubble=np.polynomial.Polynomial([0,1,-1])
            poly=poly*bubble**(2 if component=='w' else 1)
            v=poly(x); d=poly.deriv()(x); dd=poly.deriv(2)(x)
            scale=max(np.max(np.abs(poly(np.linspace(0,1,1001)))),1e-12)
        vals.append(v/scale); d1.append(d/scale); d2.append(dd/scale)
    return np.array(vals).T,np.array(d1).T,np.array(d2).T


def plate_basis(x,y,c,p):
    x=np.asarray(x); y=np.asarray(y)/c.aspect; out={}
    for comp in ['u','v','w','rx','ry']:
        bx=basis_1d(x,'x',comp,c.support,p)
        by=basis_1d(y,'y',comp,c.support,p)
        combos={}
        for name,i,j in [('v',0,0),('x',1,0),('y',0,1),('xx',2,0),('yy',0,2),('xy',1,1)]:
            combos[name]=np.einsum('qi,qj->qij',bx[i],by[j]).reshape(len(x),p*p)/c.aspect**j
        out[comp]=combos
    return out


def plate_B(b,z,c,p):
    n=len(b['w']['v']); m=p*p; B=np.zeros((n,5,5*m)); f=c.h/np.pi*np.sin(np.pi*z/c.h); fp=np.cos(np.pi*z/c.h)
    sl=lambda i:slice(i*m,(i+1)*m)
    B[:,0,sl(0)]=b['u']['x']; B[:,1,sl(1)]=b['v']['y']
    B[:,2,sl(0)]=b['u']['y']; B[:,2,sl(1)]=b['v']['x']
    B[:,0,sl(2)]=-z*b['w']['xx']; B[:,1,sl(2)]=-z*b['w']['yy']; B[:,2,sl(2)]=-2*z*b['w']['xy']
    B[:,0,sl(3)]=f*b['rx']['x']; B[:,1,sl(4)]=f*b['ry']['y']
    B[:,2,sl(3)]=f*b['rx']['y']; B[:,2,sl(4)]=f*b['ry']['x']
    B[:,3,sl(3)]=fp*b['rx']['v']; B[:,4,sl(4)]=fp*b['ry']['v']
    return B


def plate_ss(c,p=41):
    """Analytically orthogonal Navier modes for the symmetric SS benchmark.

    p odd harmonics per axis: 1,3,...,2p-1. Five unknowns per mode.
    Exact in-plane integration prevents numerical quadrature aliasing.
    """
    c.validate(); start=time.perf_counter()
    ax,by=np.meshgrid((2*np.arange(p)+1)*np.pi,
                      (2*np.arange(p)+1)*np.pi/c.aspect,indexing='ij')
    ax=ax.ravel(); by=by.ravel(); n=p*p
    K=np.zeros((n,5,5)); F=np.zeros((n,5)); area=c.aspect/4
    integral=4/(ax*by)
    gz,wz=leggauss(64); Tf,_=thermal(c)
    for z,weight in zip(gz*c.h/2,wz*c.h/2):
        E,nu,al,_,_=props(z,c); C=Cplate(E,nu)
        f=c.h/np.pi*np.sin(np.pi*z/c.h); fp=np.cos(np.pi*z/c.h)
        B=np.zeros((n,5,5))
        B[:,0,0]=-ax; B[:,0,2]=z*ax**2; B[:,0,3]=-f*ax
        B[:,1,1]=-by; B[:,1,2]=z*by**2; B[:,1,4]=-f*by
        B[:,2,0]=by; B[:,2,1]=ax; B[:,2,2]=-2*z*ax*by
        B[:,2,3]=f*by; B[:,2,4]=f*ax
        B[:,3,3]=fp; B[:,4,4]=fp
        K+=np.einsum('nia,ij,njb->nab',B,C,B,optimize=True)*area*weight
        eth=np.array([al*Tf(z),al*Tf(z),0,0,0])
        F+=np.einsum('nia,i->na',B,C@eth)*integral[:,None]*weight
    F[:,2]+=c.pressure*integral
    K[:,2,2]+=c.cB*c.h*ax**2*area
    d=np.sqrt(np.diagonal(K,axis1=1,axis2=2))
    Ks=K/d[:,:,None]/d[:,None,:]
    q=np.linalg.solve(Ks,(F/d)[...,None])[...,0]/d
    residual=np.linalg.norm(np.einsum('nij,nj->ni',K,q)-F)/max(np.linalg.norm(F),1e-30)
    return dict(case=c,kind='SPT',u=q.T.ravel(),p=p,
                seconds=time.perf_counter()-start,residual=residual,ndof=q.size)


def plate(c,p=5):
    if c.support=='SS': return plate_ss(c,p)
    c.validate(); start=time.perf_counter(); ng=max(20,4*p+4)
    g,w=leggauss(ng); x=(g+1)/2; y=x*c.aspect
    xx,yy=np.meshgrid(x,y,indexing='ij'); weights=(w[:,None]*w[None,:]*c.aspect/4).ravel()
    b=plate_basis(xx.ravel(),yy.ravel(),c,p); m=p*p; K=np.zeros((5*m,5*m)); F=np.zeros(5*m)
    gz,wz=leggauss(32); Tfun,_=thermal(c)
    for zz,ww in zip(gz*c.h/2,wz*c.h/2):
        E,nu,al,_,_=props(zz,c); C=Cplate(E,nu); B=plate_B(b,zz,c,p)
        K+=np.einsum('qia,ij,qjb,q->ab',B,C,B,weights,optimize=True)*ww
        et=np.array([al*Tfun(zz),al*Tfun(zz),0,0,0])
        F+=np.einsum('qia,i,q->a',B,C@et,weights,optimize=True)*ww
    F[2*m:3*m]+=c.pressure*np.einsum('qa,q->a',b['w']['v'],weights)
    if c.cB: K[2*m:3*m,2*m:3*m]+=c.cB*c.h*np.einsum('qa,qb,q->ab',b['w']['x'],b['w']['x'],weights)
    d=np.sqrt(np.maximum(np.diag(K),1e-30)); Ks=K/d[:,None]/d[None,:]
    u=solve(Ks,F/d,assume_a='pos')/d
    residual=np.linalg.norm(K@u-F)/max(np.linalg.norm(F),1e-30)
    return dict(case=c,kind='SPT',u=u,p=p,seconds=time.perf_counter()-start,residual=residual,ndof=len(u))


def eval_plate_ss(result,xyz):
    c=result['case']; p=result['p']; xyz=np.asarray(xyz)
    ax,by=np.meshgrid((2*np.arange(p)+1)*np.pi,
                      (2*np.arange(p)+1)*np.pi/c.aspect,indexing='ij')
    ax=ax.ravel(); by=by.ravel(); U,V,W,X,Y=result['u'].reshape(5,-1)
    x,y,z=xyz.T
    sx=np.sin(x[:,None]*ax); cx=np.cos(x[:,None]*ax)
    sy=np.sin(y[:,None]*by); cy=np.cos(y[:,None]*by)
    ss=sx*sy; cs=cx*sy; sc=sx*cy; cc=cx*cy
    f=c.h/np.pi*np.sin(np.pi*z/c.h); fp=np.cos(np.pi*z/c.h)
    u=np.c_[cs@U-z*(cs@(ax*W))+f*(cs@X),
            sc@V-z*(sc@(by*W))+f*(sc@Y),ss@W]
    eps=np.c_[ss@(-ax*U)+z*(ss@(ax*ax*W))-f*(ss@(ax*X)),
              ss@(-by*V)+z*(ss@(by*by*W))-f*(ss@(by*Y)),
              cc@(by*U+ax*V)-2*z*(cc@(ax*by*W))+f*(cc@(by*X+ax*Y)),
              fp*(cs@X),fp*(sc@Y)]
    tf=ss@(16/(c.aspect*ax*by))
    E,nu,al,_,_=props(z,c); eps[:,:2]-=(al*thermal(c)[0](z)*tf)[:,None]
    stress=np.zeros((len(xyz),6))
    stress[:,[0,1,3,4,5]]=np.einsum('qij,qj->qi',Cplate(E,nu),eps)
    return u,stress


def eval_plate(result,xyz):
    if result['case'].support=='SS': return eval_plate_ss(result,xyz)
    c=result['case']; p=result['p']; m=p*p; xyz=np.asarray(xyz); q=result['u']; b=plate_basis(xyz[:,0],xyz[:,1],c,p)
    dis=np.zeros((len(xyz),3)); stress=np.zeros((len(xyz),6))
    for k,(x,y,z) in enumerate(xyz):
        f=c.h/np.pi*np.sin(np.pi*z/c.h)
        dis[k,0]=b['u']['v'][k]@q[:m]-z*b['w']['x'][k]@q[2*m:3*m]+f*b['rx']['v'][k]@q[3*m:4*m]
        dis[k,1]=b['v']['v'][k]@q[m:2*m]-z*b['w']['y'][k]@q[2*m:3*m]+f*b['ry']['v'][k]@q[4*m:]
        dis[k,2]=b['w']['v'][k]@q[2*m:3*m]
        bb={comp:{key:v[k:k+1] for key,v in d.items()} for comp,d in b.items()}
        E,nu,al,_,_=props(z,c); eps=plate_B(bb,z,c,p)[0]@q; eps[:2]-=al*thermal(c)[0](z)
        s=Cplate(E,nu)@eps; stress[k,[0,1,3,4,5]]=s
    return dis,stress


def compare(c,nx=10,nz=32,p=101,npath=81,reference='modal'):
    if reference=='modal':
        from modal_reference import modal_solid,eval_modal
        a=modal_solid(c,p=p,nz=nz); eva=eval_modal
    else:
        a=solid(c,nx=nx,nz=nz); eva=eval_solid
    b=plate(c,p=p)
    z=np.linspace(-c.h/2,c.h/2,npath)
    centre=np.c_[np.full(npath,.5),np.full(npath,c.aspect/2),z]
    shear=np.c_[np.full(npath,.25),np.full(npath,c.aspect/2),z]
    edge_x=min(.25,c.h/2)
    edge=np.c_[np.full(npath,edge_x),np.full(npath,c.aspect/2),z]
    ua,sa=eva(a,centre); ub,sb=eval_plate(b,centre)
    _,ta=eva(a,shear); _,tb=eval_plate(b,shear)
    _,ea=eva(a,edge); _,eb=eval_plate(b,edge)
    mid=npath//2; scale=max(np.max(np.abs(sa[:,:2])),1e-12)
    norm=lambda v:np.sqrt(np.trapezoid(v*v,z)/c.h)
    sxpeak=max(np.max(np.abs(sa[:,0])),1e-15)
    row=asdict(c)|dict(reference=a['kind'],nx=nx if reference!='modal' else 0,nz=nz,p=p,w3=ua[mid,2],wp=ub[mid,2],
        w_error=100*abs(ua[mid,2]-ub[mid,2])/max(abs(ua[mid,2]),1e-15),
        sx3_max=np.max(sa[:,0])*E0/1e6,sxp_max=np.max(sb[:,0])*E0/1e6,
        sx3_abs=np.max(np.abs(sa[:,0]))*E0/1e6,sxp_abs=np.max(np.abs(sb[:,0]))*E0/1e6,
        sxpeak_error=100*abs(np.max(np.abs(sa[:,0]))-np.max(np.abs(sb[:,0])))/sxpeak,
        sx_error=100*norm(sa[:,0]-sb[:,0])/max(norm(sa[:,0]),1e-15),
        edge_x=edge_x,edge_sx_error=100*norm(ea[:,0]-eb[:,0])/max(norm(ea[:,0]),1e-15),
        edge_sx3_abs=np.max(np.abs(ea[:,0]))*E0/1e6,edge_sxp_abs=np.max(np.abs(eb[:,0]))*E0/1e6,
        edge_sz3_abs=np.max(np.abs(ea[:,2]))*E0/1e6,
        txz3_max=np.max(np.abs(ta[:,4]))*E0/1e6,txzp_max=np.max(np.abs(tb[:,4]))*E0/1e6,
        txz_error=100*norm(ta[:,4]-tb[:,4])/max(norm(ta[:,4]),1e-15),
        txz_diff_MPa=norm(ta[:,4]-tb[:,4])*E0/1e6,
        txz_error_inplane_scale=100*norm(ta[:,4]-tb[:,4])/scale,
        sz_error_scale=100*norm(sa[:,2]-sb[:,2])/scale,
        sx3_peak_zh=z[np.argmax(sa[:,0])]/c.h,sxp_peak_zh=z[np.argmax(sb[:,0])]/c.h,
        solid_seconds=a['seconds'],plate_seconds=b['seconds'],solid_residual=a['residual'],plate_residual=b['residual'],
        w_over_h=abs(ua[mid,2])/c.h)
    profiles=dict(z_h=z/c.h,u3=ua,up=ub,s3=sa,sp=sb,t3=ta,tp=tb,e3=ea,ep=eb)
    return row,profiles,a,b
