"""Independent quarter-domain Q2 solids, with refinement near both supports.

Symmetry: u=0 at x=a/2, v=0 at y=b/2; other symmetry tractions natural.
First check centre displacement against the saved full-domain uniform mesh,
then resolve the support region using graded element widths. Direct stress
recovery at an element interface can select a different one-sided derivative.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,time
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent; S=ROOT/'supplement'; sys.path.insert(0,str(S))
from plate_solvers import *
from modal_reference import modal_solid,eval_modal

def graded_solid(c,n=10,nz=8,power=2):
    start=time.perf_counter()
    xe=.5*np.linspace(0,1,n+1)**power; ye=c.aspect*xe; ze=np.linspace(-c.h/2,c.h/2,nz+1)
    def midpoints(e):
        out=np.empty(2*len(e)-1); out[::2]=e; out[1::2]=(e[:-1]+e[1:])/2; return out
    gx,gy,gz=map(midpoints,[xe,ye,ze]); nodes=np.stack(np.meshgrid(gx,gy,gz,indexing='ij'),axis=-1).reshape(-1,3)
    ids=np.arange(len(nodes)).reshape(len(gx),len(gy),len(gz))
    elems=np.array([ids[2*i:2*i+3,2*j:2*j+3,2*k:2*k+3].ravel() for i in range(n) for j in range(n) for k in range(nz)])
    dofs=(3*elems[:,:,None]+np.arange(3)).reshape(len(elems),81)
    Ke=np.zeros((len(elems),81,81)); Fe=np.zeros((len(elems),81)); g,w=leggauss(3)
    dz=c.h/nz; zc=(ze[:-1]+ze[1:])/2
    for ix in range(n):
      for iy in range(n):
        dx=xe[ix+1]-xe[ix]; dy=ye[iy+1]-ye[iy]
        pk=np.zeros((nz,81,81)); pf=np.zeros((nz,81))
        for i,r in enumerate(g):
         for j,s in enumerate(g):
          for k,t in enumerate(g):
            N,G,B=shape3(r,s,t,dx,dy,dz); z=zc+t*dz/2
            E,nu,al,_,_=props(z,c); C=C3(E,nu); fac=w[i]*w[j]*w[k]*dx*dy*dz/8
            et=np.zeros((nz,6)); et[:,:3]=(al*thermal(c)[0](z))[:,None]
            pk+=np.einsum('ia,kij,jb->kab',B,C,B,optimize=True)*fac
            pf+=np.einsum('ia,kij,kj->ka',B,C,et,optimize=True)*fac
        sl=slice((ix*n+iy)*nz,(ix*n+iy+1)*nz); Ke[sl]=pk; Fe[sl]=pf
    rows=np.repeat(dofs,81,axis=1).ravel(); cols=np.tile(dofs,(1,81)).ravel()
    K=coo_matrix((Ke.ravel(),(rows,cols)),shape=(3*len(nodes),)*2).tocsr()
    F=np.zeros(3*len(nodes)); np.add.at(F,dofs.ravel(),Fe.ravel())
    x,y,z=nodes.T; fixed=np.zeros((len(nodes),3),bool)
    fixed[np.isclose(x,0),1:3]=True
    fixed[np.isclose(y,0),0]=True; fixed[np.isclose(y,0),2]=True
    fixed[np.isclose(x,.5),0]=True; fixed[np.isclose(y,c.aspect/2),1]=True
    free=np.where(~fixed.ravel())[0]; u=np.zeros(3*len(nodes)); u[free]=spsolve(K[free][:,free],F[free])
    return dict(case=c,u=u,dofs=dofs,xe=xe,ye=ye,ze=ze,n=n,nz=nz,
                residual=np.linalg.norm((K@u-F)[free])/np.linalg.norm(F[free]),seconds=time.perf_counter()-start)

def evaluate(r,points):
    c=r['case']; us=[]; stresses=[]
    for xyz in points:
        index=[]; xi=[]; delta=[]
        for x,edges in zip(xyz,[r['xe'],r['ye'],r['ze']]):
            i=min(len(edges)-2,max(0,np.searchsorted(edges,x,side='right')-1)); dx=edges[i+1]-edges[i]
            index.append(i); delta.append(dx); xi.append(2*(x-edges[i])/dx-1)
        N,G,B=shape3(*xi,*delta); ix,iy,iz=index; e=(ix*r['n']+iy)*r['nz']+iz
        ue=r['u'][r['dofs'][e]]; E,nu,al,_,_=props(xyz[2],c); eps=B@ue; eps[:3]-=al*thermal(c)[0](xyz[2])
        us.append(N@ue.reshape(27,3)); stresses.append(C3(E,nu)@eps)
    return np.array(us),np.array(stresses)

def main():
    c=Case(); z=np.linspace(-c.h/2,c.h/2,81)
    xyz=np.concatenate([np.c_[np.full(81,x),np.full(81,.5),z] for x in [.5,.05]])
    modal=modal_solid(c,p=241,nz=48); um,sm=eval_modal(modal,xyz)
    rows=[]; norm=lambda v:np.sqrt(np.trapezoid(v*v,z))
    for n,nz,power in [(6,8,1),(10,8,2),(14,10,2)]:
        a=graded_solid(c,n=n,nz=nz,power=power); u,s=evaluate(a,xyz)
        rec=dict(quarter_n=n,nz=nz,grading_power=power,seconds=a['seconds'],residual=a['residual'],
            w_difference_percent=100*abs(u[40,2]/um[40,2]-1),
            centre_stress_difference_percent=100*norm(s[:81,0]-sm[:81,0])/norm(sm[:81,0]),
            edge_stress_difference_percent=100*norm(s[81:,0]-sm[81:,0])/norm(sm[81:,0]),
            top_edge_sigma_zz_MPa=s[-1,2]*70000)
        if power==1:
            old=np.load(S/'results/paper_thermal_solid_n12.npz')['solid_stress']
            rec['quarter_full_stress_difference_percent']=100*np.linalg.norm(old-s)/np.linalg.norm(old)
        rows.append(rec); print(rec,flush=True)
        pd.DataFrame(rows).to_csv(S/'tables/paper_graded_solid_check.csv',index=False)
        np.savez_compressed(S/'results'/f'paper_graded_n{n}_nz{nz}.npz',z_h=z/c.h,solid_stress=s,modal_stress=sm)

if __name__=='__main__':main()
