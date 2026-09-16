# Executed Python study: thermal porous FGM plates

The package contains executed deterministic analyses, source code, numerical verification, convergence runs, CSV tables, and seven figures in PNG and SVG. No Bayesian inference was performed.

**Run it**

Requires Python 3.10 or later with NumPy, SciPy, pandas and Matplotlib. The recorded run used Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0.

```bash
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run_study.py
python plot_results.py
```

On Windows, set the two thread environment variables using your shell's syntax, or simply run the Python commands. The full study took about three minutes in the execution environment; hardware and thread settings change runtime. `python run_study.py --quick` uses lower resolution and fewer joint samples; it overwrites the output files, so use a separate copy if preserving this full run.

`plot_results.py` reads saved results and never reruns the solvers. For a new simply supported case:

```python
from plate_solvers import Case, compare
case = Case(slenderness=10, gradation=1, porosity=0.10, deltaT=20)
row, profiles, reference_solution, plate_solution = compare(case, p=121, nz=32)
print(row['w_error'], row['edge_sx_error'])
```

For a clamped case use, for example, `compare(Case(support='CCCC'), nx=10, nz=8, p=5, reference='solid')`. The modal reference supports SS diaphragm boundaries only.

The optional `pressure` input is the top-face transverse traction divided by $E_0$, positive in the $+z$ direction.

**Files**

| File or directory | Purpose |
|---|---|
| `plate_solvers.py` | Material and thermal laws; independent 27-node 3D solid; sinusoidal plate; response comparisons |
| `modal_reference.py` | 3D elasticity with Fourier expansions in the plane and quadratic elements through thickness |
| `run_study.py` | Verification gates, convergence, parameter sweeps and saved outputs |
| `plot_results.py` | Separate plotting script |
| `RESULTS.md` | Findings, numerical qualifications and implications for the four questions |
| `tables/all_cases.csv` | 76 labelled rows representing 70 distinct model comparisons |
| `tables/verification.csv` | Seven passed checks and their tolerances |
| `tables/convergence_deltas.csv` | Numerical changes relative to finer calculations |
| `tables/rq*.csv` | Tables for each research question |
| `results/*.npz` | Through-thickness displacement and stress profiles for each labelled case |
| `results/run_metadata.json` | Inputs, bounds, seed, software versions, runtime and summaries |
| `results/execution.log` | Actual execution output |
| `figures/` | Seven figures, each exported as PNG and SVG |

**Selected benchmark inputs**

These inputs are analyst-selected; they are not asserted to be the original proposal's values.

| Input | Baseline |
|---|---|
| Geometry | Square plate, $a/h=10$, $b/a=1$ |
| Absolute size | Unspecified: calculations use lengths divided by $a$ |
| Constituent pair | Aluminium / silicon carbide |
| Aluminium | $E=70$ GPa, $\nu=0.30$, $\alpha=23.4\times10^{-6}$ K$^{-1}$, $k=233$ W m$^{-1}$ K$^{-1}$ |
| Silicon carbide | $E=427$ GPa, $\nu=0.17$, $\alpha=4.3\times10^{-6}$ K$^{-1}$, $k=65$ W m$^{-1}$ K$^{-1}$ |
| Grading | Ceramic fraction $V_c=t^g$, $t=z/h+1/2$, $g=1$ |
| Mean porosity | 0.10, uniform through thickness |
| Temperature rise | Bottom 0 K, top 20 K, measured from the stress-free temperature |
| Loading | Thermal only; no pressure or restoring operator in the baseline |
| Supports | SS diaphragms as defined below |
| Constitutive scope | Local, stationary, small-strain, linear thermoelasticity |

The nominal Al/SiC values are based on the material example in [Vel and Batra (2003)](https://doi.org/10.1016/S0020-7683(03)00361-5). The present homogenisation and porosity laws are assumptions for this benchmark; their paper is not being reproduced.

With $P_\mathrm{mix}=P_m+(P_c-P_m)V_c$, the effective laws are

$$E(z)=s_E E_\mathrm{mix}(z)[1-\phi(z)]^{r_E},\qquad k(z)=k_\mathrm{mix}(z)[1-\phi(z)]^{r_k},$$

$$\nu(z)=\nu_\mathrm{mix}(z),\qquad \alpha(z)=s_\alpha\alpha_\mathrm{mix}(z).$$

Baseline $r_E=2$, $r_k=1$, $s_E=s_\alpha=1$. Properties do not vary with temperature. Porosity alters these effective properties; individual pores and damage are not resolved.

The three pore profiles have the **same mean void fraction** $\bar\phi$:

$$\phi_\mathrm{uniform}=\bar\phi,\quad
\phi_\mathrm{centre}=\bar\phi[1+\cos(2\pi z/h)],\quad
\phi_\mathrm{faces}=\bar\phi[1-\cos(2\pi z/h)].$$

Steady conduction satisfies $d(kT_{,z})/dz=0$. The solution is obtained by integrating $1/k(z)$ between the prescribed face temperatures. There is no transient heat equation in this study.

**Mechanical models and boundary conditions**

The 3D reference uses the full isotropic elasticity tensor and thermal eigenstrain:

$$\boldsymbol\sigma=\mathsf C(z):[\boldsymbol\varepsilon-\alpha(z)T(z)\mathsf I].$$

For SS diaphragms, $v=w=0$ on the $x$ edges and $u=w=0$ on the $y$ edges, at every thickness coordinate. The remaining normal displacement is free and its traction is natural. These are specific diaphragm supports, not every possible interpretation of “simply supported”. For CCCC, all three displacements vanish on every lateral face. The top and bottom mechanical surfaces are traction-free except in the explicitly labelled pressure diagnostic.

The primary SS reference uses

$$u=\sum U_{mn}(z)\cos(\alpha_m x)\sin(\beta_n y),\quad
v=\sum V_{mn}(z)\sin(\alpha_m x)\cos(\beta_n y),\quad
w=\sum W_{mn}(z)\sin(\alpha_m x)\sin(\beta_n y).$$

Here $\alpha_m=m\pi/a$, $\beta_n=n\pi/b$, with positive odd $m,n$. Odd modes suffice for the symmetric geometry, material distribution and loads considered. **$U,V,W$ vary independently through thickness**; this is a discretisation of 3D elasticity, not a plate theory. Four-point quadrature is used in each quadratic thickness element. The main SS runs use 121 odd harmonics per axis and 32 quadratic thickness elements, with refinements to 241 harmonics and 64 elements. Fourier stresses use the same projected thermal eigenstrain as the modal weak form. Truncation is checked explicitly.

An independent 3D solver uses tensor-product 27-node quadratic solid elements and full $3\times3\times3$ integration. It supplies the independent checks and the clamped control cases. Clamped cases use $10\times10\times8$ and $12\times12\times8$ elements. Stress is evaluated directly from the finite element strain, without smoothing.

The selected five-field sinusoidal plate approximation is

$$u=u_0-z w_{0,x}+f(z)r_x,\quad
v=v_0-z w_{0,y}+f(z)r_y,\quad w=w_0,\qquad
f(z)=\frac{h}{\pi}\sin\frac{\pi z}{h}.$$

The reduced constitutive matrix assumes $\sigma_{zz}=0$. The shear factor $f'(z)$ vanishes at both faces. There is no additional shear correction factor. SS in-plane modal integrals are analytical; thickness integrals use 64-point Gauss quadrature. For CCCC, boundary-compatible polynomial Ritz functions enforce zero displacement and zero slope for $w_0$ at the edges.

This plate model omits transverse stretching. A thesis formulation that includes a separate normal-deformation field would be a different model and could change the results.

**Conditional magnetic-restoring experiment**

The proposal's exact Lorentz/Maxwell equations could not be recovered. RQ3 therefore uses the expressly assumed positive restoring energy

$$U_B=\frac{1}{2}\int_\Omega M_B(w_{,x})^2\,dV,\qquad
\frac{M_B}{E_0}=\frac{\eta_B(h/a)^2}{12(1-0.3^2)},\quad E_0=70\ \mathrm{GPa}.$$

The same term is added to both models, with $\eta_B=0,10,30,100,300,1000$. It tests the response to directional transverse restoring stiffness. **It is not a verified magnetic constitutive law, has no calibrated conversion to tesla, and includes no magnetic effect on material properties.** It also does not calculate magnetic damping, buckling, or dynamic stability. Results from this experiment are conditional on this operator.

**Response definitions and units**

Code lengths are divided by $a$; code stresses by $E_0=70$ GPa. `w3` and `wp` are $w/a$ at the centre of the midplane; multiply by a chosen physical $a$ to obtain displacement. `w_over_h` is $|w_{3D}|/h$. CSV stress peaks are in MPa; NPZ stresses are dimensionless and need multiplication by 70,000 for MPa.

Normal-stress profiles have 81 equally spaced thickness samples. The centre path is $(a/2,b/2,z)$; the near-edge path is $(h/2,b/2,z)$ in all executed thermal cases. The latter stays a fixed distance **in thickness units**, so it remains inside the edge region when the plate becomes thin. Shear is also sampled at $(a/4,b/2,z)$.

$$e_w=100\frac{|w_\mathrm{plate}-w_\mathrm{3D}|}{|w_\mathrm{3D}|},\qquad
e_\sigma=100\frac{\|\sigma_\mathrm{plate}-\sigma_\mathrm{3D}\|_{L^2(z)}}{\|\sigma_\mathrm{3D}\|_{L^2(z)}}.$$

`sx3_abs` is the maximum absolute $\sigma_{xx}$ on the **centre path**, not a global maximum, a von Mises stress, or a failure measure. `sx3_max` is the signed maximum on that path. `edge_sx3_abs` is the corresponding absolute peak on the near-edge path. `sx3_peak_zh` records the signed tensile-maximum location; it is not the absolute-peak location.

`txz_error` is retained for audit but should not be interpreted when thermal shear is nearly zero or unresolved. Use `txz_diff_MPa` and `txz_error_inplane_scale` with the convergence evidence. Likewise, a large percentage displacement error in a nearly motionless clamped case can be misleading.

`screen_pass` uses illustrative requirements of 5% deflection error, 10% centre normal-stress profile error and 10% near-edge profile error, jointly. These are exploratory thresholds, not a standard or a demonstrated safety limit.

**Deterministic robustness design**

Twelve one-at-a-time runs use the following input bounds. Sixteen Latin-hypercube points explore all six ranges jointly, with seed 20260912. The sampling is a space-filling design, **not a probability distribution for the inputs**. Observed minima and maxima are not guaranteed bounds over the entire parameter box.

| Parameter | Low | High |
|---|---:|---:|
| Grading exponent | 0.5 | 2.0 |
| Mean porosity | 0.05 | 0.20 |
| Porosity–modulus exponent | 1.5 | 2.5 |
| Porosity–conductivity exponent | 0.5 | 1.5 |
| Modulus scale | 0.85 | 1.15 |
| Expansion scale | 0.85 | 1.15 |

**Scope of evidence**

Seven numerical checks passed: affine solid patch, free thermal expansion, constant heat flux, thin-pressure limits for both models, and independent solid/modal displacement and stress comparisons. These verify selected numerical behaviours; they do not provide experimental validation. Residuals reported for the modal 3D solver sample its first, middle and final mode systems.

No finite nonlocal length, nonlinear geometry, moving plate, viscoelastic response, temperature-dependent properties, pore damage, buckling or dynamic stability is implemented. A nonzero `nonlocal_length` raises `NotImplementedError`; the code does not imitate nonlocality by rescaling elastic modulus. Absolute nanoscale behaviour cannot be inferred from this size-independent local model.

The supporting papers are useful for defining subsequent work, but are not reproduced here: [Gong et al. (2023)](https://www.amm.shu.edu.cn/EN/10.1007/s10483-023-2943-8), [Thi et al. (online 2025)](https://doi.org/10.1007/s10409-025-24824-x), and the related [Jalaei and Thai (2019)](https://doi.org/10.1016/j.compositesb.2019.107164). Their material systems, kinematics and couplings must be matched before using their numerical results as validation targets.

The next formulation-specific step requires the original proposal and supporting PDFs: extract their exact kinematics, effective properties, thermal and mechanical boundaries, nonlocal constitutive boundary conditions and electromagnetic law; replace the corresponding assumptions; then repeat the same verification and comparison workflow. A current literature comparison is also needed before claiming a research gap or novelty.
