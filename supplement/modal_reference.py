"""Three-dimensional elasticity, Fourier in x,y and quadratic FE in z.

No plate kinematic assumptions: U(z), V(z), W(z) are independent fields.
Only the symmetric, uniformly heated/loaded SS diaphragm case is supported.
The same full C3 material law as the independent 27-node solid is used.
"""
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve_banded
from plate_solvers import C3, props, thermal, lagrange2


def modal_solid(c,p=101,nz=24):
    c.validate()
    if c.support!='SS': raise ValueError('Modal reference requires SS diaphragms.')
    start=time.perf_counter(); nd=3*(2*nz+1); dz=c.h/nz
    # Coefficients of 1, alpha, beta, alpha², beta², alpha*beta.
    bands=np.zeros((6,17,nd)); loads=np.zeros((3,nd)); g,w=leggauss(4)
    for ez in range(nz):
        dofs=np.arange(6*ez,6*ez+9)
        ke=np.zeros((6,9,9)); fe=np.zeros((3,9))
        for xi,wi in zip(g,w):
            L,D=lagrange2(xi); D=D*2/dz
            z=-c.h/2+(ez+(xi+1)/2)*dz
            E,nu,al,_,_=props(z,c); C=C3(E,nu)
            B0=np.zeros((6,9)); Ba=B0.copy(); Bb=B0.copy()
            B0[2,2::3]=D; B0[4,0::3]=D; B0[5,1::3]=D
            Ba[0,0::3]=-L; Ba[3,1::3]=L; Ba[4,2::3]=L
            Bb[1,1::3]=-L; Bb[3,0::3]=L; Bb[5,2::3]=L
            weight=wi*dz/2
            ke+=np.array([B0.T@C@B0,
                B0.T@C@Ba+Ba.T@C@B0,B0.T@C@Bb+Bb.T@C@B0,
                Ba.T@C@Ba,Bb.T@C@Bb,Ba.T@C@Bb+Bb.T@C@Ba])*weight
            mag=np.zeros(9); mag[2::3]=L
            ke[3]+=c.cB*np.outer(mag,mag)*weight
            eth=np.r_[np.full(3,al*thermal(c)[0](z)),np.zeros(3)]
            fe+=np.array([B0.T@C@eth,Ba.T@C@eth,Bb.T@C@eth])*weight
        for i,di in enumerate(dofs):
            loads[:,di]+=fe[:,i]
            for j,dj in enumerate(dofs): bands[:,8+di-dj,dj]+=ke[:,i,j]
    ax,by=np.meshgrid((2*np.arange(p)+1)*np.pi,
                      (2*np.arange(p)+1)*np.pi/c.aspect,indexing='ij')
    ax=ax.ravel(); by=by.ravel(); q=np.zeros((p*p,nd)); worst=0.0
    for k,(a,b) in enumerate(zip(ax,by)):
        ab=np.einsum('i,ijk->jk',[1,a,b,a*a,b*b,a*b],bands)
        f=(loads[0]+a*loads[1]+b*loads[2])*16/(c.aspect*a*b)
        f[-1]+=c.pressure*16/(c.aspect*a*b)
        q[k]=solve_banded((8,8),ab,f,check_finite=False)
        if k in [0,p*p//2,p*p-1]:
            residual=-f.copy()
            for row in range(17):
                offset=row-8
                j=np.arange(max(0,-offset),min(nd,nd-offset)); i=j+offset
                residual[i]+=ab[row,j]*q[k,j]
            worst=max(worst,np.linalg.norm(residual)/max(np.linalg.norm(f),1e-30))
    return dict(case=c,kind='3D_Fourier_Q2',u=q.reshape(p*p,2*nz+1,3),
                ax=ax,by=by,p=p,nz=nz,ndof=q.size,
                seconds=time.perf_counter()-start,residual=worst)


def eval_modal(result,xyz):
    c=result['case']; nz=result['nz']; dz=c.h/nz
    ax=result['ax']; by=result['by']; dis=[]; stresses=[]
    for x,y,z in np.asarray(xyz):
        ez=min(nz-1,max(0,int((z+c.h/2)/dz)))
        xi=2*(z+c.h/2-ez*dz)/dz-1; L,D=lagrange2(xi); D*=2/dz
        local=result['u'][:,2*ez:2*ez+3,:]
        U,V,W=np.einsum('i,kij->jk',L,local)
        Uz,Vz,Wz=np.einsum('i,kij->jk',D,local)
        sx=np.sin(ax*x); cx=np.cos(ax*x); sy=np.sin(by*y); cy=np.cos(by*y)
        ss=sx*sy; cs=cx*sy; sc=sx*cy; cc=cx*cy
        dis.append([cs@U,sc@V,ss@W])
        eps=np.array([ss@(-ax*U),ss@(-by*V),ss@Wz,
                      cc@(by*U+ax*V),cs@(Uz+ax*W),sc@(Vz+by*W)])
        # Reconstruct the modal stress, using the same Fourier projection of
        # the uniform thermal eigenstrain as in the weak-form right-hand side.
        # This avoids adding an unprojected thermal tail to truncated stresses.
        tf=ss@(16/(c.aspect*ax*by))
        E,nu,al,_,_=props(z,c); eps[:3]-=al*thermal(c)[0](z)*tf
        stresses.append(C3(E,nu)@eps)
    return np.array(dis),np.array(stresses)
