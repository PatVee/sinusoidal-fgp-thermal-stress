"""Run the reproducible, explicitly illustrative local plate study.

Recommended: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run_study.py
All outputs are relative to this script. See README.md before interpreting them.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import argparse, json, platform, time
from dataclasses import asdict, replace
import numpy as np
import pandas as pd
import scipy
from scipy.stats import qmc
from plate_solvers import *
from modal_reference import modal_solid, eval_modal

ROOT=Path(__file__).resolve().parent
for folder in ['tables','results','figures']: (ROOT/folder).mkdir(exist_ok=True)


def safe_json(obj):
    if isinstance(obj,np.generic): return obj.item()
    if isinstance(obj,np.ndarray): return obj.tolist()
    raise TypeError(type(obj).__name__)


def verify():
    records=[]
    def check(name,value,limit,unit):
        records.append(dict(check=name,value=float(value),limit=limit,unit=unit,passed=bool(value<=limit)))
        print('VERIFY',name,value,'PASS' if value<=limit else 'FAIL',flush=True)
    c=Case(material='Al',porosity=0,deltaT=0)
    A=np.array([[.001,.0002,0],[.0001,-.0004,0],[0,0,.0003]])
    r=solid(c,nx=2,nz=2,patch=lambda x:x@A.T)
    exact=r['nodes']@A.T
    check('3D affine displacement patch',np.max(np.abs(r['u'].reshape(-1,3)-exact))/np.max(np.abs(exact)),1e-10,'relative')
    c=replace(c,bottomT=20)
    r=solid(c,nx=2,nz=2,patch=lambda x:ALPHA0*20*x)
    _,s=eval_solid(r,[[.2,.3,-.02],[.6,.4,0],[.7,.8,.02]])
    check('Free thermal expansion stress',np.max(np.abs(s))*E0/1e6,1e-6,'MPa')
    c=Case(pattern='centre'); z=np.linspace(-c.h/2,c.h/2,501)
    T,flux=thermal(c); q=-props(z,c)[3]*T.derivative()(z)
    check('Steady heat flux constancy',np.max(np.abs(q-flux))/abs(flux),2e-5,'relative')
    c=Case(material='Al',porosity=0,deltaT=0,pressure=1e-9,slenderness=100)
    m,n=np.meshgrid(2*np.arange(201)+1,2*np.arange(201)+1,indexing='ij')
    D=c.h**3/(12*(1-.3**2))
    wkl=np.sum(16*c.pressure*np.sin(m*np.pi/2)*np.sin(n*np.pi/2)/
               (np.pi**6*m*n*D*(m*m+n*n)**2))
    for name,r,ev in [('SPT',plate(c,p=61),eval_plate),('3D modal',modal_solid(c,p=61,nz=16),eval_modal)]:
        u,_=ev(r,[[.5,.5,0]])
        check(name+' thin-pressure Navier limit',abs(u[0,2]/wkl-1),.002,'relative')
    c=Case(material='Al',porosity=0,deltaT=0,pressure=1e-7,slenderness=20)
    a=solid(c,nx=10,nz=6); b=modal_solid(c,p=121,nz=32)
    points=np.c_[np.full(81,.5),np.full(81,.5),np.linspace(-c.h/2,c.h/2,81)]
    ua,sa=eval_solid(a,points); ub,sb=eval_modal(b,points)
    check('Independent solid/modal pressure displacement',abs(ua[40,2]/ub[40,2]-1),.005,'relative')
    check('Independent solid/modal pressure normal stress',np.linalg.norm(sa[:,0]-sb[:,0])/np.linalg.norm(sb[:,0]),.02,'relative L2')
    pd.DataFrame(records).to_csv(ROOT/'tables/verification.csv',index=False)
    if not all(r['passed'] for r in records): raise RuntimeError('A verification gate failed. Review before running study.')
    return records


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--quick',action='store_true')
    args=parser.parse_args(); start=time.perf_counter()
    p,nz=(61,16) if args.quick else (121,32)
    base=Case(); rows=[]; cache={}
    checks=verify()
    def execute(tag,c,group,p_override=None,nz_override=None,nx=10):
        pp=p_override or (p if c.support=='SS' else 5)
        zz=nz_override or (nz if c.support=='SS' else 8)
        ref='modal' if c.support=='SS' else 'solid'
        key=(c,pp,zz,nx,ref)
        if key in cache: row,pr=cache[key]
        else:
            row,pr,*_=compare(c,nx=nx,nz=zz,p=pp,reference=ref)
            cache[key]=(row,pr)
        row=row|dict(case_id=tag,group=group)
        row['screen_pass']=bool(row['w_error']<=5 and row['sx_error']<=10 and row['edge_sx_error']<=10)
        rows.append(row)
        np.savez_compressed(ROOT/'results'/f'{tag}.npz',**pr)
        pd.DataFrame(rows).to_csv(ROOT/'tables/all_cases.csv',index=False)
        print(f"RUN {len(rows):02d} {tag}: w error={row['w_error']:.3f}%, centre stress={row['sx_error']:.3f}%, edge stress={row['edge_sx_error']:.3f}%",flush=True)
        return row

    execute('baseline',base,'rq1')
    execute('thick_high_gradation',replace(base,slenderness=5,gradation=5,porosity=.2),'rq1')
    execute('pressure_diagnostic',Case(material='Al',porosity=0,deltaT=0,pressure=1e-7,slenderness=20),'diagnostic')
    mode_levels=[31,61,121] if args.quick else [31,61,121,241]
    for pp in mode_levels:
        execute(f'conv_modes_{pp}',base,'convergence_modes',p_override=pp,nz_override=32)
    for zz in ([8,16,32] if args.quick else [8,16,32,64]):
        execute(f'conv_thickness_{zz}',base,'convergence_thickness',p_override=p,nz_override=zz)
    if not args.quick:
        execute('conv_thick_fine',replace(base,slenderness=5,gradation=5,porosity=.2),'convergence_extreme',p_override=241,nz_override=64)
        execute('conv_thin_fine',replace(base,slenderness=40,gradation=.5,porosity=.2),'convergence_extreme',p_override=241,nz_override=64)
    for phi in [0,.2]:
      for grad in [.5,1,5]:
       for slender in [5,10,20,40]:
        execute(f'map_s{slender}_g{grad}_phi{phi}',replace(base,slenderness=slender,gradation=grad,porosity=phi),'rq2')
    for eta in [0,10,30,100,300,1000]:
        execute(f'magnetic_eta{eta}',replace(base,etaB=eta),'rq3')
    bounds=dict(gradation=(.5,2.0),porosity=(.05,.2),E_exponent=(1.5,2.5),
                k_exponent=(.5,1.5),E_scale=(.85,1.15),alpha_scale=(.85,1.15))
    for name,(lo,hi) in bounds.items():
        for label,value in [('low',lo),('high',hi)]:
            execute(f'oat_{name}_{label}',replace(base,**{name:value}),'rq4_oat')
    for pattern in ['uniform','centre','faces']:
        execute(f'pattern_{pattern}',replace(base,pattern=pattern),'rq4_pattern')
    execute('restraint_CCCC',replace(base,support='CCCC'),'rq4_restraint',nx=10,nz_override=8)
    if not args.quick:
        execute('restraint_CCCC_fine',replace(base,support='CCCC'),'convergence_restraint',nx=12,nz_override=8,p_override=6)
    design=qmc.LatinHypercube(d=len(bounds),seed=20260912).random(8 if args.quick else 16)
    sampled=qmc.scale(design,[v[0] for v in bounds.values()],[v[1] for v in bounds.values()])
    for i,values in enumerate(sampled):
        execute(f'bounded_joint_{i:02d}',replace(base,**dict(zip(bounds,values))),'rq4_joint')
    df=pd.DataFrame(rows)
    for name,group in df.groupby('group'):
        group.to_csv(ROOT/'tables'/f'{name}.csv',index=False)
    b=df[df.case_id=='baseline'].iloc[0]
    last=df[df.case_id=='magnetic_eta1000'].iloc[0]
    extremes=df[df.group=='rq2']; joint=df[df.group=='rq4_joint']
    summary=dict(baseline={k:b[k] for k in ['w3','wp','w_over_h','w_error','sx3_abs','sxp_abs','sx_error','edge_sx_error']},
                 magnetic_eta1000=dict(deflection_change_percent=100*(abs(last.w3)/abs(b.w3)-1),
                     centre_peak_stress_change_percent=100*(last.sx3_abs/b.sx3_abs-1)),
                 rq2=dict(cases=len(extremes),screen_passes=int(extremes.screen_pass.sum()),
                     w_error_range=[extremes.w_error.min(),extremes.w_error.max()],
                     centre_stress_error_range=[extremes.sx_error.min(),extremes.sx_error.max()],
                     edge_stress_error_range=[extremes.edge_sx_error.min(),extremes.edge_sx_error.max()]),
                 bounded_design=dict(samples=len(joint),w_over_h_range=[joint.w_over_h.min(),joint.w_over_h.max()],
                     peak_stress_MPa_range=[joint.sx3_abs.min(),joint.sx3_abs.max()]))
    metadata=dict(run_utc=pd.Timestamp.now(tz='UTC').isoformat(),python=platform.python_version(),
                  numpy=np.__version__,scipy=scipy.__version__,quick=args.quick,main_p=p,main_nz=nz,
                  elapsed_seconds=time.perf_counter()-start,row_count=len(rows),unique_comparisons=len(cache),
                  baseline=asdict(base),bounds=bounds,random_seed=20260912,
                  interpretation='Deterministic sensitivity design; no inferred probabilities or posterior.',
                  verification=checks,summary=summary)
    (ROOT/'results/run_metadata.json').write_text(json.dumps(metadata,indent=2,default=safe_json))
    print(json.dumps(summary,indent=2,default=safe_json),flush=True)
    print('DONE',metadata['elapsed_seconds'],'seconds',flush=True)


if __name__=='__main__': main()
