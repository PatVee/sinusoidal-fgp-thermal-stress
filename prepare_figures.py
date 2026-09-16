"""Editable vector problem/method figures and plots from executed data only."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,FancyBboxPatch,FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap

R=Path(__file__).resolve().parent; S=R/'supplement'; F=R/'figures'
D=pd.read_csv(S/'tables/all_cases.csv'); B=D[D.case_id=='baseline'].iloc[0]
N='#17324D'; C='#237C9E'; O='#D16A32'; G='#368270'; P='#885A9D'; GR='#788590'; LIGHT='#EDF2F5'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold',
 'text.color':N,'axes.labelcolor':N,'grid.alpha':.18,'pdf.fonttype':42,'ps.fonttype':42})
def save(fig,name):
    for ext in ['pdf','svg','png']:fig.savefig(F/f'{name}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def arrow(ax,a,b,color=GR,style='-|>',lw=1.2):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=11,lw=lw,color=color))
def box(ax,x,y,w,h,title,body,color=C):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.05,rounding_size=0.12',fc='white',ec=color,lw=1.4))
    ax.text(x+.14,y+h-.2,title,va='top',fontsize=10,fontweight='bold',color=color)
    ax.text(x+.14,y+h-.62,body,va='top',fontsize=8.7,linespacing=1.4,color=N)

def problem():
    fig=plt.figure(figsize=(7.8,5.6),layout='constrained'); gs=fig.add_gridspec(2,3,width_ratios=[1.2,1.2,1])
    ax=fig.add_subplot(gs[:,:2]); ax.set(xlim=(0,10),ylim=(0,6.5)); ax.axis('off')
    ax.text(.1,6.22,'(a) Plate and sampling paths',fontsize=11,fontweight='bold')
    def q(x,y,z):return np.array([1+5.8*x+1.6*y,2+1.7*y+.65*z])
    cmap=LinearSegmentedColormap.from_list('fgm',[C,'#CED9DD',O])
    for k in range(20):
        z0=k/20; z1=(k+1)/20
        ax.add_patch(Polygon([q(0,0,z0),q(1,0,z0),q(1,0,z1),q(0,0,z1)],fc=cmap((k+.5)/20),ec='none'))
        ax.add_patch(Polygon([q(1,0,z0),q(1,1,z0),q(1,1,z1),q(1,0,z1)],fc=cmap((k+.5)/20),ec='none',alpha=.8))
    ax.add_patch(Polygon([q(0,0,1),q(1,0,1),q(1,1,1),q(0,1,1)],fc='#F4D6C2',ec=N,lw=1.1))
    for face in [[q(0,0,0),q(1,0,0),q(1,0,1),q(0,0,1)],[q(1,0,0),q(1,1,0),q(1,1,1),q(1,0,1)]]:
        ax.add_patch(Polygon(face,fill=False,ec=N,lw=1))
    for x in [.12,.35,.58,.81]:
        ax.add_patch(Polygon([q(x-.035,0,-.32),q(x+.035,0,-.32),q(x,0,-.03)],fc='#CDD5DA',ec=GR,lw=.8))
    for x,col,label in [(.5,G,'C'),(.05,P,'E')]:
        a=q(x,.5,-.1); b=q(x,.5,1.23); ax.plot([a[0],b[0]],[a[1],b[1]],ls='--',lw=2,color=col)
        ax.scatter(*q(x,.5,.5),s=30,color=col,zorder=10)
        ax.text(b[0]-.15,b[1]+.1,label,color=col,fontweight='bold',fontsize=12)
    ax.text(5.5,5.58,r'Hot face: $\theta^+=20\ \mathrm{K}$',color=O,ha='center',fontsize=11,fontweight='bold')
    for x in [.48,.65,.82]:
        t=q(x,.72,1); arrow(ax,t+[0,1.1],t+[0,.1],color=O)
    ax.text(3.9,1.39,r'Cold face: $\theta^-=0$',color=C,ha='center',fontsize=10)
    arrow(ax,(1,.96),(6.8,.96),style='<->'); ax.text(3.9,.75,r'$a$',ha='center',fontsize=12)
    arrow(ax,q(1.06,0,-.65),q(1.06,1,-.65),style='<->'); ax.text(8.3,2.15,r'$b$',fontsize=12)
    arrow(ax,q(1.1,1,0),q(1.1,1,1),style='<->'); ax.text(9.18,4.01,r'$h$',va='center',fontsize=12)
    origin=np.array([.58,4.65])
    for vec,lab in [([.7,0],'x'),([.35,.38],'y'),([0,.75],'z')]:
        arrow(ax,origin,origin+vec,color=N); ax.text(*(origin+np.array(vec)+[.04,.02]),f'${lab}$',fontsize=11)
    ax.text(.15,.38,r'$x=0,a:\ v=w=0,\ \sigma_{xx}=0$'+ '\n'+r'$y=0,b:\ u=w=0,\ \sigma_{yy}=0$',fontsize=8.8,va='center')
    ax.text(5.15,.38,'C: centre\nE: x = h/2, y = b/2',fontsize=8.6,va='center')
    z=np.linspace(-.5,.5,301)
    bx=fig.add_subplot(gs[0,2]);
    for g,col in zip([.5,1,5],[C,O,G]):bx.plot((z+.5)**g,z,lw=1.7,color=col,label=f'g = {g:g}')
    bx.set(xlabel=r'Ceramic fraction $V_c$',ylabel=r'$z/h$',title='(b) Grading'); bx.grid(True); bx.legend(frameon=False,fontsize=8)
    cx=fig.add_subplot(gs[1,2]);
    for val,lab,col,ls in [(.1+0*z,'Uniform',C,'-'),(.1*(1+np.cos(2*np.pi*z)),'Centre-rich',G,'--'),(.1*(1-np.cos(2*np.pi*z)),'Face-rich',O,'-.')]:cx.plot(val,z,lw=1.7,label=lab,color=col,ls=ls)
    cx.set(xlabel=r'Void fraction $\phi$',ylabel=r'$z/h$',title='(c) Pore profiles'); cx.grid(True); cx.legend(frameon=False,fontsize=8)
    save(fig,'fig01_problem')

def methodology():
    fig,ax=plt.subplots(figsize=(7.4,6.5)); ax.set(xlim=(0,12),ylim=(0,8.5)); ax.axis('off')
    box(ax,.2,7.42,11.6,.8,'COMMON PROBLEM DEFINITION','Geometry, effective properties, steady temperature and explicit supports',N)
    for x in [2,6,10]:arrow(ax,(x,7.35),(x,6.96))
    box(ax,.2,5.05,3.45,1.85,'3D solid elements','27-node solid elements\nFull and quarter domains\nUniform / graded meshes',C)
    box(ax,4.25,5.05,3.45,1.85,'Fourier / 3D reference','Fourier modes in x and y\nIndependent U(z), V(z), W(z)\nQuadratic z elements',G)
    box(ax,8.3,5.05,3.45,1.85,'Sinusoidal plate','Five displacement fields\nSinusoidal shear variation\nPlane-stress reduction',O)
    arrow(ax,(1.95,5),(1.95,4.45),color=C); arrow(ax,(6,5),(6,4.45),color=G); arrow(ax,(10,5),(10,4.45),color=O)
    box(ax,.2,3.13,3.45,1.22,'VERIFICATION','Patch / thin-plate limits\nThermal stress cross-check',C)
    box(ax,4.25,3.13,7.5,1.22,'MATCHED RESPONSES','Displacement and stress at common thickness paths\nCheck Fourier order, thickness and local mesh resolution',N)
    arrow(ax,(3.72,3.74),(4.14,3.74),color=C)
    for x in [1.6,4.5,7.5,10.4]:
        arrow(ax,(8,3.08),(x,2.57),color=GR)
    cards=[(.2,'RQ1 | Accuracy','Profiles and peaks\nCentre / support',C),
           (3.15,'RQ2 | Applicability','Thickness / grading\nSupport distance',G),
           (6.1,'RQ3 | Restoration','Positive stiffness\nStress / deflection',O),
           (9.05,'RQ4 | Robustness','Properties / pores\nRestraints / bounds',P)]
    for x,title,body,col in cards:box(ax,x,1.15,2.65,1.35,title,body,col)
    ax.text(6,.45,'Local-model verification; physical calibration and finite-size effects remain separate.',ha='center',fontsize=8.7,color=GR)
    fig.tight_layout(); save(fig,'fig02_methodology')

def verification():
    fig,axs=plt.subplots(1,3,figsize=(7.4,3.8),layout='constrained')
    a=np.load(S/'results/paper_thermal_solid_n16.npz'); b=np.load(S/'results/paper_graded_n14_nz10.npz')
    for ax,comp,title in [(axs[0],0,r'(a) Near-edge $\sigma_{xx}$'),(axs[1],2,r'(b) Near-edge $\sigma_{zz}$')]:
        ax.plot(b['modal_stress'][81:,comp]*70000,b['z_h'],color=G,lw=2,label='Fourier / 3D')
        ax.plot(a['solid_stress'][81:,comp]*70000,a['z_h'],color=GR,ls=':',lw=1.8,label='Uniform solid mesh')
        ax.plot(b['solid_stress'][81:,comp]*70000,b['z_h'],color=C,ls='--',lw=1.8,label='Graded solid mesh')
        ax.set(xlabel='Stress (MPa)',ylabel=r'$z/h$',title=title); ax.grid(True)
    axs[0].legend(frameon=False,fontsize=8)
    ax=axs[2]; e=.5*np.linspace(0,1,15)**2
    for v in e:ax.plot([v,v],[0,.5],color=C,lw=.55);ax.plot([0,.5],[v,v],color=C,lw=.55)
    ax.plot([.05,.05],[0,.5],color=P,ls='--',lw=1.2)
    ax.scatter(.05,.5,s=32,color=P,zorder=4)
    ax.set(xlim=(-.015,.515),ylim=(-.015,.515),aspect='equal',xlabel=r'$x/a$',ylabel=r'$y/b$',title='(c) Graded quarter mesh')
    ax.text(.25,.24,'14 × 14 in-plane\n10 through thickness',ha='center',va='center',fontsize=8,bbox=dict(fc='white',ec='none',alpha=.94))
    save(fig,'fig03_reference_check')

def baseline():
    fig,axs=plt.subplots(1,2,figsize=(7.4,3.8),layout='constrained')
    for ax,tag,title in [(axs[0],'paper_edge_s10_d5','(a) Centre'),(axs[1],'paper_edge_s10_d0.5','(b) Near edge, x = h/2')]:
        q=np.load(S/'results'/f'{tag}.npz')
        ax.plot(q['s3'][:,0]*70000,q['z_h'],color=G,lw=2,label='3D elasticity')
        ax.plot(q['sp'][:,0]*70000,q['z_h'],color=O,lw=1.8,ls='--',label='Sinusoidal plate')
        ax.set(xlabel=r'$\sigma_{xx}$ (MPa)',ylabel=r'$z/h$',title=title); ax.grid(True);ax.legend(frameon=False,fontsize=8)
    save(fig,'fig04_profiles')

def distance():
    d=pd.read_csv(S/'tables/paper_edge_distance.csv'); fig,ax=plt.subplots(figsize=(6.9,4.1),layout='constrained')
    for slender,col,mark in zip([5,10,20,40],[C,G,O,P],['o','s','^','D']):
        q=d[(d.slenderness==slender)&(d.d_over_h<=3)]
        ax.semilogy(q.d_over_h,q.stress_error,marker=mark,color=col,lw=1.6,markersize=4,label=f'a/h = {slender}')
    ax.axhline(10,ls=':',color=GR,label='10% comparison threshold')
    ax.axvspan(.75,1,color='#D9E6E1',alpha=.5)
    ax.set(xlabel=r'Distance from support $d/h$',ylabel=r'Normal-stress profile error $e_\sigma$ (%)',ylim=(.04,80),xlim=(.18,3.04))
    ax.grid(True,which='both');ax.legend(frameon=False,ncol=2,fontsize=8)
    save(fig,'fig06_edge_distance')

def error_map():
    fig,axs=plt.subplots(1,2,figsize=(7.4,3.6),layout='constrained'); data=D[(D.group=='rq2')&(D.porosity==.2)]
    for ax,key,title,limit in [(axs[0],'w_error','(a) Centre displacement',15),(axs[1],'edge_sx_error','(b) Near-edge stress profile',50)]:
        q=data.pivot(index='gradation',columns='slenderness',values=key).sort_index().to_numpy()
        im=ax.imshow(q,cmap='YlOrRd',vmin=0,vmax=limit,aspect='auto')
        for (i,j),v in np.ndenumerate(q):ax.text(j,i,f'{v:.2f}' if key=='w_error' else f'{v:.1f}',ha='center',va='center',color='white' if v>limit*.65 else N,fontsize=10)
        ax.set(xticks=range(4),xticklabels=[5,10,20,40],yticks=range(3),yticklabels=[.5,1,5],xlabel='Slenderness a/h',ylabel='Grading exponent g',title=title)
        fig.colorbar(im,ax=ax,shrink=.8,label='Error (%)')
    save(fig,'fig05_parameter_map')

def sensitivity():
    fig,axs=plt.subplots(1,3,figsize=(7.4,3.8),layout='constrained')
    names=['gradation','porosity','E_exponent','k_exponent','E_scale','alpha_scale']; labels=[r'$g$',r'$\bar\phi$',r'$r_E$',r'$r_k$',r'$s_E$',r'$s_\alpha$']
    q=pd.read_csv(S/'tables/sensitivity_changes.csv')
    for ax,key,title in [(axs[0],'w3','(a) Centre displacement'),(axs[1],'sx3_abs','(b) Centre peak stress')]:
        t=q[q.response==key].set_index('parameter').loc[names]
        for i,(_,r) in enumerate(t.iterrows()):
            ax.plot([r.low_change_percent,r.high_change_percent],[i,i],color='#CAD3DA',lw=3)
            ax.scatter(r.low_change_percent,i,s=25,color=C,label='Low input' if i==0 else None)
            ax.scatter(r.high_change_percent,i,s=25,color=O,label='High input' if i==0 else None)
        ax.set(yticks=range(6),yticklabels=labels,xlabel='Change (%)',title=title);ax.invert_yaxis();ax.grid(axis='x');ax.axvline(0,lw=.7,color=GR);ax.legend(frameon=False,fontsize=8)
    t=D[D.group=='rq4_joint'];ax=axs[2]
    im=ax.scatter(t.w_over_h*1e4,t.sx3_abs,c=t.porosity,cmap='viridis',s=32)
    ax.scatter(B.w_over_h*1e4,B.sx3_abs,s=100,marker='*',color=O,label='Baseline')
    ax.set(xlabel=r'$10^4|w|/h$',ylabel=r'Centre peak $|\sigma_{xx}|$ (MPa)',title='(c) Joint design');ax.grid(True);ax.legend(frameon=False,fontsize=8)
    fig.colorbar(im,ax=ax,shrink=.8,label='Mean porosity')
    save(fig,'fig07_sensitivity')

def pores():
    tags=['pattern_uniform','pattern_centre','pattern_faces','restraint_CCCC_fine']; labels=['Uniform\nSS','Centre-rich\nSS','Face-rich\nSS','Uniform\nCCCC']
    t=D.set_index('case_id').loc[tags];fig,axs=plt.subplots(1,2,figsize=(7.4,3.6),layout='constrained');x=np.arange(4)
    for ax,k1,k2,title in [(axs[0],'w3','wp','(a) Centre displacement'),(axs[1],'sx3_abs','sxp_abs','(b) Centre-path peak normal stress')]:
        fac=t.slenderness.to_numpy()*1e4 if k1=='w3' else 1
        ax.bar(x-.16,t[k1]*fac,width=.3,color=G,label='3D');ax.bar(x+.16,t[k2]*fac,width=.3,color=O,label='Plate')
        ax.set(xticks=x,xticklabels=labels,title=title,ylabel=r'$10^4w/h$' if k1=='w3' else r'Peak $|\sigma_{xx}|$ (MPa)');ax.axhline(0,color=GR,lw=.7);ax.grid(axis='y');ax.legend(frameon=False,fontsize=8)
    save(fig,'fig08_pores_restraint')

def restoration():
    t=D[D.group=='rq3'].sort_values('etaB');r0=t.iloc[0]
    fig,axs=plt.subplots(1,3,figsize=(7.4,3.5),layout='constrained')
    axs[0].plot(t.etaB,np.abs(t.w3/r0.w3),'o-',color=C,label='3D displacement')
    axs[0].plot(t.etaB,t.sx3_abs/r0.sx3_abs,'s-',color=O,label='3D centre stress')
    axs[0].set(ylabel='Ratio to unrestrained baseline',title='(a) Response trends');axs[0].legend(frameon=False,fontsize=8)
    axs[1].plot(t.etaB,t.sx3_abs,'o-',color=O,label='Centre');axs[1].plot(t.etaB,t.edge_sx3_abs,'s-',color=G,label='Near edge')
    axs[1].set(ylabel=r'Peak $|\sigma_{xx}|$ (MPa)',title='(b) Stress redistribution');axs[1].legend(frameon=False,fontsize=8)
    axs[2].plot(t.etaB,t.w_error,'o-',color=C,label='Displacement');axs[2].plot(t.etaB,t.edge_sx_error,'s-',color=O,label='Near-edge stress')
    axs[2].set(ylabel='Plate / 3D error (%)',title='(c) Approximation error');axs[2].legend(frameon=False,fontsize=8)
    for ax in axs:ax.set_xscale('symlog',linthresh=10);ax.set_xlabel(r'Restoring strength $\eta_R$');ax.grid(True)
    save(fig,'fig09_restoration')

if __name__=='__main__':
    F.mkdir(exist_ok=True)
    for fn in [problem,methodology,verification,baseline,distance,error_map,sensitivity,pores,restoration]:fn()
    print('Nine vector figures written as PDF, SVG and 300-dpi PNG.')
