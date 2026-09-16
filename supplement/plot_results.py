"""Regenerate every figure from saved CSV/NPZ outputs; no solver reruns."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT=Path(__file__).resolve().parent
F=ROOT/'figures'; F.mkdir(exist_ok=True)
D=pd.read_csv(ROOT/'tables/all_cases.csv')
META=json.loads((ROOT/'results/run_metadata.json').read_text())
B=D[D.case_id=='baseline'].iloc[0]
NAVY='#16324f'; BLUE='#247ba0'; ORANGE='#db6d28'; GREEN='#3a806d'; GRAY='#697785'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
    'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,
    'axes.labelcolor':NAVY,'text.color':NAVY,'axes.titleweight':'bold',
    'grid.alpha':.2,'figure.facecolor':'white','savefig.facecolor':'white'})

def save(fig,name):
    fig.savefig(F/f'{name}.png',dpi=190,bbox_inches='tight')
    fig.savefig(F/f'{name}.svg',bbox_inches='tight')
    plt.close(fig)

def profile(tag): return np.load(ROOT/'results'/f'{tag}.npz')

def stress_panel(ax,pr,key,title):
    ax.plot(pr[key+'3'][:,0]*70000,pr['z_h'],color=BLUE,lw=2.2,label='3D elasticity')
    ax.plot(pr['sp' if key=='s' else 'ep'][:,0]*70000,pr['z_h'],color=ORANGE,lw=2,ls='--',label='Sinusoidal plate')
    ax.set(xlabel=r'$\sigma_{xx}$ (MPa)',ylabel=r'Thickness coordinate $z/h$',title=title)
    ax.grid(True); ax.axvline(0,color=GRAY,lw=.7)


def overview():
    fig,axs=plt.subplots(2,2,figsize=(11.7,8.0),layout='constrained')
    pr=profile('baseline')
    stress_panel(axs[0,0],pr,'s','Centre: stress profiles nearly coincide')
    stress_panel(axs[0,1],pr,'e',f'Near edge: profile error {B.edge_sx_error:.1f}%')
    axs[0,0].legend(frameon=False,fontsize=9)
    cases=D[(D.group=='rq2')&(D.porosity==.2)]
    for g,col in zip([.5,1,5],[BLUE,GREEN,ORANGE]):
        t=cases[cases.gradation==g].sort_values('slenderness')
        axs[1,0].plot(t.slenderness,t.w_error,'o-',color=col,label=f'Grading exponent {g:g}')
    axs[1,0].axhline(5,color=GRAY,ls=':',label='5% screening level')
    axs[1,0].set(xlabel='Slenderness a/h',ylabel='Centre deflection error (%)',title='Thickness affects deflection accuracy')
    axs[1,0].legend(frameon=False,fontsize=8); axs[1,0].grid(True)
    t=D[D.group=='rq3'].sort_values('etaB'); t0=t.iloc[0]
    axs[1,1].plot(t.etaB,np.abs(t.w3/t0.w3),'o-',color=BLUE,label='Centre deflection / baseline')
    axs[1,1].plot(t.etaB,t.sx3_abs/t0.sx3_abs,'s-',color=ORANGE,label='Centre peak stress / baseline')
    axs[1,1].set(xscale='symlog',xlabel=r'Idealised restoring strength $\eta_B$',ylabel='Ratio to zero-field operator',title='Restoring action trades deflection for stress')
    axs[1,1].legend(frameon=False,fontsize=8); axs[1,1].grid(True)
    fig.suptitle('Executed Python benchmark: local, linear thermoelasticity\nIllustrative inputs; magnetic panel is a restoring-operator sensitivity',fontsize=15)
    save(fig,'00_overview')


def profiles():
    fig,axs=plt.subplots(1,3,figsize=(12.5,4.8),layout='constrained')
    pr=profile('baseline')
    stress_panel(axs[0],pr,'s',r'Centre: $x/a=y/b=0.5$')
    stress_panel(axs[1],pr,'e',r'Near edge: $x=h/2,\ y/b=0.5$')
    axs[0].legend(frameon=False,fontsize=9)
    axs[2].plot(pr['s3'][:,2]*70000,pr['z_h'],color=BLUE,label='3D, centre')
    axs[2].plot(pr['e3'][:,2]*70000,pr['z_h'],color=GREEN,label='3D, near edge')
    axs[2].axvline(0,color=ORANGE,ls='--',label='Plate: identically zero')
    axs[2].set(xlabel=r'$\sigma_{zz}$ (MPa)',ylabel=r'$z/h$',title='Transverse normal stress')
    axs[2].legend(frameon=False,fontsize=8); axs[2].grid(True)
    fig.suptitle('RQ1 | Stress accuracy depends on where it is measured\nBaseline: a/h = 10, grading exponent = 1, mean porosity = 0.10, top rise = 20 K',fontsize=14)
    save(fig,'01_stress_profiles')


def applicability():
    fig,axs=plt.subplots(2,3,figsize=(12.8,7.4),layout='constrained')
    specs=[('w_error','Centre deflection',15),('sx_error','Centre normal stress',1),('edge_sx_error','Near-edge normal stress',60)]
    for i,phi in enumerate([0,.2]):
      data=D[(D.group=='rq2')&(D.porosity==phi)]
      for j,(key,title,vmax) in enumerate(specs):
        tab=data.pivot(index='gradation',columns='slenderness',values=key).sort_index()
        vals=tab.to_numpy(); im=axs[i,j].imshow(vals,cmap='YlOrRd',vmin=0,vmax=max(vmax,vals.max()),aspect='auto')
        axs[i,j].set(xticks=np.arange(4),xticklabels=['5','10','20','40'],yticks=np.arange(3),yticklabels=['0.5','1','5'],xlabel='Slenderness a/h',ylabel='Grading exponent',title=f'{title}\nMean porosity = {phi:g}')
        for (r,c),v in np.ndenumerate(vals):
            axs[i,j].text(c,r,f'{v:.2f}' if j==1 else f'{v:.1f}',ha='center',va='center',fontsize=10,color='white' if v>.64*max(vmax,vals.max()) else NAVY)
        fig.colorbar(im,ax=axs[i,j],shrink=.78,label='Error (%)')
    fig.suptitle('RQ2 | Error maps for 24 thermal cases\nNear-edge path stays at x = h/2; these are sampled errors, not a validated design envelope',fontsize=14)
    save(fig,'02_error_maps')


def magnetic():
    fig,axs=plt.subplots(1,3,figsize=(12.8,4.3),layout='constrained')
    t=D[D.group=='rq3'].sort_values('etaB')
    axs[0].plot(t.etaB,np.abs(t.w3)*t.slenderness,'o-',color=BLUE,label='3D')
    axs[0].plot(t.etaB,np.abs(t.wp)*t.slenderness,'s--',color=ORANGE,label='Plate')
    axs[0].set(ylabel=r'Centre displacement $|w|/h$',title='Deflection'); axs[0].legend(frameon=False)
    axs[1].plot(t.etaB,t.sx3_abs,'o-',color=ORANGE,label='Centre path')
    axs[1].plot(t.etaB,t.edge_sx3_abs,'s-',color=GREEN,label='Near-edge path')
    axs[1].set(ylabel=r'Peak $|\sigma_{xx}|$ along path (MPa)',title='3D stress'); axs[1].legend(frameon=False)
    axs[2].plot(t.etaB,t.w_error,'o-',color=BLUE,label='Deflection')
    axs[2].plot(t.etaB,t.edge_sx_error,'s-',color=ORANGE,label='Near-edge stress')
    axs[2].set(ylabel='Plate / 3D difference (%)',title='Approximation error'); axs[2].legend(frameon=False,fontsize=9)
    for ax in axs:
        ax.set_xscale('symlog',linthresh=10); ax.set_xlabel(r'Restoring parameter $\eta_B$'); ax.grid(True)
    fig.suptitle('RQ3 | Conditional result for an idealised transverse restoring energy\nNo conversion to tesla, full Maxwell coupling, buckling, or dynamic stability is claimed',fontsize=14)
    save(fig,'03_restoring_operator')


def sensitivity():
    specs=[('gradation','Grading exponent'),('porosity','Mean porosity'),('E_exponent','Porosity–modulus exponent'),('k_exponent','Porosity–conductivity exponent'),('E_scale','Modulus scale'),('alpha_scale','Expansion scale')]
    fig,axs=plt.subplots(1,2,figsize=(12.2,5.0),layout='constrained')
    out=[]
    for ax,key,title in zip(axs,['w3','sx3_abs'],['Centre displacement','Centre peak normal stress']):
      for i,(name,label) in enumerate(specs):
        lo=D[D.case_id==f'oat_{name}_low'].iloc[0]; hi=D[D.case_id==f'oat_{name}_high'].iloc[0]
        a=100*(lo[key]/B[key]-1); b=100*(hi[key]/B[key]-1)
        ax.plot([a,b],[i,i],color='#bac4ce',lw=4,zorder=1)
        ax.scatter(a,i,color=BLUE,s=40,label='Low input' if i==0 else None,zorder=2)
        ax.scatter(b,i,color=ORANGE,s=40,label='High input' if i==0 else None,zorder=2)
        out.append(dict(parameter=name,response=key,low_input=lo[name],high_input=hi[name],low_change_percent=a,high_change_percent=b))
      ax.set(yticks=range(len(specs)),yticklabels=[s[1] for s in specs],xlabel='Change from baseline (%)',title=title)
      ax.invert_yaxis(); ax.axvline(0,color=GRAY,lw=1); ax.grid(axis='x',alpha=.2); ax.legend(frameon=False,fontsize=9)
    pd.DataFrame(out).to_csv(ROOT/'tables/sensitivity_changes.csv',index=False)
    fig.suptitle('RQ4 | One-at-a-time bounded sensitivities\nAssumed ranges for exploration; no posterior or confidence interval',fontsize=14)
    save(fig,'04_sensitivity')


def joint_design():
    t=D[D.group=='rq4_joint']; fig,ax=plt.subplots(figsize=(7.8,5.4),layout='constrained')
    im=ax.scatter(t.w_over_h,t.sx3_abs,c=t.porosity,cmap='viridis',s=75,edgecolor='white',lw=.6)
    ax.scatter(B.w_over_h,B.sx3_abs,marker='*',s=200,color=ORANGE,label='Baseline',zorder=3)
    ax.set(xlabel=r'3D centre displacement $|w|/h$',ylabel=r'Centre-path peak $|\sigma_{xx}|$ (MPa)',title=f'RQ4 | Joint bounded design: {len(t)} deterministic samples')
    ax.grid(True); ax.legend(frameon=False); fig.colorbar(im,ax=ax,label='Mean porosity')
    save(fig,'05_joint_design')


def convergence():
    fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained'); records=[]
    for ax,group,key,title in [(axs[0],'convergence_modes','p','Fourier truncation'),(axs[1],'convergence_thickness','nz','Through-thickness elements')]:
        t=D[D.group==group].sort_values(key); fine=t.iloc[-1]; pref=profile(fine.case_id)
        wdiff=[]; sdiff=[]; ediff=[]; levels=[]
        for _,r in t.iloc[:-1].iterrows():
            pp=profile(r.case_id); z=pp['z_h']
            norm=lambda x:np.sqrt(np.trapezoid(x*x,z))
            wd=100*abs(r.w3/fine.w3-1)
            sd=100*norm(pp['s3'][:,0]-pref['s3'][:,0])/norm(pref['s3'][:,0])
            ed=100*norm(pp['e3'][:,0]-pref['e3'][:,0])/norm(pref['e3'][:,0])
            wdiff.append(wd); sdiff.append(sd); ediff.append(ed); levels.append(r[key])
            records.append(dict(group=group,level=r[key],reference_level=fine[key],w_change_percent=wd,centre_stress_profile_change_percent=sd,edge_stress_profile_change_percent=ed))
        for val,label,col in [(wdiff,'Deflection',BLUE),(sdiff,'Centre stress',GREEN),(ediff,'Near-edge stress',ORANGE)]:
            ax.semilogy(levels,val,'o-',color=col,label=label)
        ax.set(xlabel='Odd harmonics per axis' if key=='p' else 'Quadratic elements in z',ylabel='3D change from finest run (%)',title=title)
        ax.grid(True,which='both'); ax.legend(frameon=False,fontsize=9)
    pd.DataFrame(records).to_csv(ROOT/'tables/convergence_deltas.csv',index=False)
    fig.suptitle('Numerical convergence | Baseline local thermal case\nComparison of numerical resolutions; not physical validation',fontsize=14)
    save(fig,'06_convergence')


if __name__=='__main__':
    overview(); profiles(); applicability(); magnetic(); sensitivity(); joint_design(); convergence()
    print('Saved seven figures as PNG and SVG.')
