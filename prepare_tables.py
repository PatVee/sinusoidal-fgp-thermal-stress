from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parent; S=R/'supplement'; T=R/'tables'; T.mkdir(exist_ok=True)
D=pd.read_csv(S/'tables/all_cases.csv'); B=D[D.case_id=='baseline'].iloc[0]
def table(name,cols,head,rows):
    body='\\begin{tabularx}{\\linewidth}{'+cols+'}\n\\toprule\n'+' & '.join(head)+r' \\'+ '\n\\midrule\n'
    body+='\n'.join(' & '.join(map(str,r))+r' \\' for r in rows)
    (T/f'{name}.tex').write_text(body+'\n\\bottomrule\n\\end{tabularx}\n')
def sci(v):
    m,e=f'{v:.3e}'.split('e');return '$'+m+r'\times10^{'+str(int(e))+'}$'
table('materials','Xrrl',['Property','Aluminium','Silicon carbide','Unit'],[
    [r'Young\textquotesingle s modulus, $E$',70,427,'GPa'],[r'Poisson ratio, $\nu$',.30,.17,'--'],
    [r'Thermal expansion, $10^6\alpha$',23.4,4.3,r'K$^{-1}$'],[r'Thermal conductivity, $k$',233,65,r'W m$^{-1}$ K$^{-1}$']])
v=pd.read_csv(S/'tables/verification.csv')
labels=['Affine displacement patch (relative error)','Free thermal expansion (stress, MPa)','Steady heat-flux constancy (relative error)','SPT / thin-plate pressure limit (relative error)','3D / thin-plate pressure limit (relative error)','Independent pressure displacement (relative error)','Independent pressure normal stress (relative error)']
table('verification','Xrrc',['Check','Measured','Tolerance','Result'],[[lab,sci(r.value),sci(r.limit),'Pass'] for lab,(_,r) in zip(labels,v.iterrows())])
c=pd.read_csv(S/'tables/paper_thermal_crosscheck.csv');g=pd.read_csv(S/'tables/paper_graded_solid_check.csv')
rows=[]
for _,r in c.iterrows():rows.append([f'Uniform, full ${int(r.nx)}^2\\times{int(r.nz)}$',f'{r.w_difference_percent:.4f}',f'{r.centre_stress_difference_percent:.3f}',f'{r.edge_stress_difference_percent:.3f}'])
for _,r in g[g.grading_power==2].iterrows():rows.append([f'Graded, quarter ${int(r.quarter_n)}^2\\times{int(r.nz)}$',f'{r.w_difference_percent:.4f}',f'{r.centre_stress_difference_percent:.3f}',f'{r.edge_stress_difference_percent:.3f}'])
table('thermal_check','Xrrr',['Solid mesh',r'$w_C$ error (\%)',r'$\sigma_C$ error (\%)',r'$\sigma_E$ error (\%)'],rows)
cv=pd.read_csv(S/'tables/convergence_deltas.csv')
table('convergence','Xrrr',['Refinement comparison',r'$w_C$ (\%)',r'$\sigma_C$ (\%)',r'$\sigma_E$ (\%)'],[
 [f"{'Modes' if r.group=='convergence_modes' else 'Thickness elements'}: {int(r.level)} to {int(r.reference_level)}",f'{r.w_change_percent:.6f}',f'{r.centre_stress_profile_change_percent:.3f}',f'{r.edge_stress_profile_change_percent:.3f}'] for _,r in cv.iterrows()])
table('baseline','Xrrr',['Quantity','3D','Plate',r'Difference (\%)'],[
 [r'Centre displacement, $10^5w_C/a$',f'{B.w3*1e5:.5f}',f'{B.wp*1e5:.5f}',f'{B.w_error:.3f}'],
 [r'Centre peak $|\sigma_{xx}|$ (MPa)',f'{B.sx3_abs:.3f}',f'{B.sxp_abs:.3f}',f'{B.sxpeak_error:.3f}'],
 [r'Edge-path peak $|\sigma_{xx}|$ (MPa)',f'{B.edge_sx3_abs:.3f}',f'{B.edge_sxp_abs:.3f}',f'{100*abs(B.edge_sxp_abs/B.edge_sx3_abs-1):.3f}'],
 [r'Centre profile, $e_{\sigma,C}$','--','--',f'{B.sx_error:.3f}'],
 [r'Edge profile, $e_{\sigma,E}$','--','--',f'{B.edge_sx_error:.3f}']])
q=D[(D.group=='rq2')&(D.porosity==.2)].sort_values(['gradation','slenderness'])
table('parameter_map','Xrrrr',['Grading $g$','$a/h$',r'$e_w$ (\%)',r'$e_{\sigma,C}$ (\%)',r'$e_{\sigma,E}$ (\%)'],[
 [f'{r.gradation:g}',f'{r.slenderness:g}',f'{r.w_error:.3f}',f'{r.sx_error:.3f}',f'{r.edge_sx_error:.3f}'] for _,r in q.iterrows()])
ed=pd.read_csv(S/'tables/paper_edge_distance.csv');rows=[]
for slender in [5,10,20,40]:
 q=ed[ed.slenderness==slender].set_index('d_over_h')
 rows.append([slender]+[f'{q.loc[d,"stress_error"]:.2f}' for d in [.25,.5,.75,1,1.5]])
table('distance','Xrrrrr',[r'$a/h$',r'$d/h=0.25$','0.50','0.75','1.00','1.50'],rows)
sens=pd.read_csv(S/'tables/sensitivity_changes.csv');rows=[]
for name,label in [('gradation','$g$'),('porosity',r'$\bar\phi$'),('E_exponent','$r_E$'),('k_exponent','$r_k$'),('E_scale','$s_E$'),('alpha_scale',r'$s_\alpha$')]:
 a=sens[(sens.parameter==name)&(sens.response=='w3')].iloc[0];b=sens[(sens.parameter==name)&(sens.response=='sx3_abs')].iloc[0]
 rows.append([label,f'{a.low_input:g}--{a.high_input:g}']+[f'{v:+.2f}' if abs(v)>.005 else '0.00' for v in [a.low_change_percent,a.high_change_percent,b.low_change_percent,b.high_change_percent]])
table('sensitivity','Xrrrrr',['Input','Range',r'$\Delta w_-$ (\%)',r'$\Delta w_+$ (\%)',r'$\Delta S_-$ (\%)',r'$\Delta S_+$ (\%)'],rows)
rows=[]
for tag,label in [('pattern_uniform','Uniform, SS'),('pattern_centre','Centre-rich, SS'),('pattern_faces','Face-rich, SS'),('restraint_CCCC_fine','Uniform, CCCC')]:
 r=D[D.case_id==tag].iloc[0];rows.append([label,f'{r.w3*r.slenderness*1e4:.4f}',f'{r.wp*r.slenderness*1e4:.4f}',f'{r.sx3_abs:.3f}',f'{r.sxp_abs:.3f}'])
table('pores','Xrrrr',['Profile and support',r'$10^4 w_{3D}/h$',r'$10^4 w_P/h$',r'$S_{3D}$ (MPa)',r'$S_P$ (MPa)'],rows)
q=D[D.group=='rq3'].sort_values('etaB')
table('restoration','Xrrrr',[r'$\eta_R$',r'$10^4|w_{3D}|/h$',r'$S_{3D,C}$ (MPa)',r'$S_{3D,E}$ (MPa)',r'$e_w$ (\%)'],[
 [int(r.etaB),f'{r.w_over_h*1e4:.4f}',f'{r.sx3_abs:.3f}',f'{r.edge_sx3_abs:.3f}',f'{r.w_error:.3f}'] for _,r in q.iterrows()])
print('Ten publication tables generated from saved outputs.')
