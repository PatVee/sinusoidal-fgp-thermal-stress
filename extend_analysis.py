"""Additional thermal verification and edge-distance analysis for the paper."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import sys,json,time
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent; S=ROOT/'supplement'; sys.path.insert(0,str(S))
from plate_solvers import Case,solid,eval_solid,plate,eval_plate,E0
from modal_reference import modal_solid,eval_modal

def norm(v,z): return np.sqrt(np.trapezoid(v*v,z))

def main():
    start=time.perf_counter(); rows=[]; cross=[]
    for slender in [5,10,20,40]:
        c=Case(slenderness=slender)
        a=modal_solid(c,p=241,nz=48); b=plate(c,p=241)
        z=np.linspace(-c.h/2,c.h/2,81)
        distances=sorted(set([d for d in [.25,.5,.75,1,1.5,2,2.5,3,4,5,7.5,10,15,20] if d<=slender/2]+[slender/2]))
        for d in distances:
            xyz=np.c_[np.full(81,d*c.h),np.full(81,.5),z]
            ua,sa=eval_modal(a,xyz); ub,sb=eval_plate(b,xyz)
            row=dict(slenderness=slender,d_over_h=d,x_over_a=d*c.h,
                 stress_error=100*norm(sa[:,0]-sb[:,0],z)/norm(sa[:,0],z),
                 ref_stress_norm_MPa=norm(sa[:,0],z)/np.sqrt(c.h)*70000,
                 ref_peak_MPa=np.max(np.abs(sa[:,0]))*70000,
                 plate_peak_MPa=np.max(np.abs(sb[:,0]))*70000)
            rows.append(row)
            np.savez_compressed(S/'results'/f'paper_edge_s{slender}_d{d:g}.npz',z_h=z/c.h,s3=sa,sp=sb,u3=ua,up=ub)
        print('DISTANCE',slender,rows[-1],flush=True)
        pd.DataFrame(rows).to_csv(S/'tables/paper_edge_distance.csv',index=False)
        if slender==10:
            allpts=np.concatenate([np.c_[np.full(81,x),np.full(81,.5),z] for x in [.5,.05]])
            um,sm=eval_modal(a,allpts)
            for nx in [12,16]:
                sol=solid(c,nx=nx,nz=8); uf,sf=eval_solid(sol,allpts)
                r=dict(nx=nx,nz=8,solid_seconds=sol['seconds'],modal_p=241,modal_nz=48,
                    w_difference_percent=100*abs(uf[40,2]/um[40,2]-1),
                    centre_stress_difference_percent=100*norm(sf[:81,0]-sm[:81,0],z)/norm(sm[:81,0],z),
                    edge_stress_difference_percent=100*norm(sf[81:,0]-sm[81:,0],z)/norm(sm[81:,0],z),
                    solid_residual=sol['residual'])
                cross.append(r); print('THERMAL CROSSCHECK',r,flush=True)
                pd.DataFrame(cross).to_csv(S/'tables/paper_thermal_crosscheck.csv',index=False)
                np.savez_compressed(S/'results'/f'paper_thermal_solid_n{nx}.npz',z_h=z/c.h,solid_stress=sf,modal_stress=sm)
    meta=dict(elapsed_seconds=time.perf_counter()-start,case_count=4,solid_checks=2,paths=len(rows))
    (S/'results/paper_extension_metadata.json').write_text(json.dumps(meta,indent=2))
    print('DONE',meta,flush=True)

if __name__=='__main__': main()
