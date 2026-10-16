# EM-4 Exercises

Conventions as in the notes: SI units, $k = 1/4\pi\varepsilon_0 = 8.988\times10^9$ N m²/C², $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$ with $\hat{\mathbf n}$ out of the dielectric, $\rho_b = -\nabla\cdot\mathbf P$, $\mathbf D = \varepsilon_0\mathbf E + \mathbf P$, $\varepsilon_r = 1 + \chi_e$. Equation numbers refer to the notes. Every number asked for below is printed or asserted by `project.py`; the solutions say where.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (The atom as a stretched ball.) Hydrogen has $\alpha/4\pi\varepsilon_0 = 0.667$ Å$^3$.
(a) Use (3.2) to find the radius $a$ of the uniform-ball model with this polarizability; compare with the Bohr radius $0.529$ Å.
(b) In $E = 10^6$ V/m find $p = \alpha E$ and the displacement $d = p/e$ of the nucleus from the cloud's centre.
(c) The restoring field at $d = a$ is $E_a = e/4\pi\varepsilon_0a^2$. Evaluate it and show $d/a = E/E_a$. When must $\mathbf p = \alpha\mathbf E$ fail?

**A2.** (Bound charge by hand.) A cube $0\le x,y,z\le L$ has $\mathbf P = kx\,\hat{\mathbf x}$ with $k = 10^{-4}$ C/m³, $L = 2$ cm. Find $\rho_b$, $\sigma_b$ on each of the six faces, and the total bound charge. Which faces carry none, and why?

**A3.** (Torque against thermal energy.) A water molecule has $p = 1.85$ D $= 6.17\times10^{-30}$ C m.
(a) The maximum torque in $E = 10^5$ V/m. (b) $2pE/k_BT$ at $300$ K. (c) The field at which $2pE/k_BT = 1$, compared with air's breakdown field $3$ MV/m. What does (b) say about the form of $P(E)$ for a polar liquid in laboratory fields?

**A4.** (A filled capacitor.) Plates of area $50$ cm², gap $1$ mm, filled with $\varepsilon_r = 2.1$ (Teflon), connected to $100$ V. Find $C$, $Q$, $\sigma_f$, $E$, $P$, and $\sigma_b$ on the dielectric face touching the positive plate. Verify $\sigma_f + \sigma_b = \varepsilon_0E$ and say what it means.

**A5.** (Refraction of field lines.) At an interface between $\varepsilon_1 = \varepsilon_0$ and $\varepsilon_2 = 5\varepsilon_0$ with no free charge, $\mathbf E$ in medium 1 makes $\theta_1 = 30°$ with the normal.
(a) Derive $\tan\theta_2/\tan\theta_1 = \varepsilon_2/\varepsilon_1$ from (3.16)–(3.17). (b) Find $\theta_2$ and $|\mathbf E_2|/|\mathbf E_1|$. (c) Toward or away from the normal? What does $\varepsilon_2\to\infty$ do to the lines on the vacuum side?

## B. Core (do all)

**B1.** (Radial polarization, two ways.) A sphere of radius $R$ has $\mathbf P = kr\,\hat{\mathbf r}$ ($k$ constant, C/m³) and no free charge.
(a) Find $\rho_b$, $\sigma_b$; show the total bound charge is zero.
(b) Find $\mathbf E$ everywhere by Gauss's law (EM-2 (3.5)) using the bound charge alone.
(c) Find $\mathbf D$ from (3.15) and symmetry, then $\mathbf E$ from $\mathbf D = \varepsilon_0\mathbf E + \mathbf P$. Do the routes agree?
(d) $R = 1$ cm, $k = 10^{-4}$ C/m³: $\rho_b$, $\sigma_b$, $E$ just inside $r = R$, and the jump of $E$ at $R$ against (3.16).

**B2.** (Total charge and total dipole moment.) Extend $\mathbf P$ of a bounded object by zero outside, and take $\rho_b = -\nabla\cdot\mathbf P$ in the distribution sense (the surface layer $\sigma_b$ is then included automatically).
(a) Prove $\int\rho_b\,d\tau = 0$ over all space without separating volume and surface terms.
(b) Prove $\int\mathbf r\,\rho_b\,d\tau = \int\mathbf P\,d\tau$, using $\nabla\cdot(x_i\mathbf P) = P_i + x_i\nabla\cdot\mathbf P$.
(c) Apply (b) to the uniformly polarized sphere. (d) What do (a) and (b) say about the field outside B1's sphere?
(e) In project Part 1, `P_blob_offset` has a net moment. Change its offset $0.5$ to $-1.0$, predict the dipole moment of the grid bound charge, then run.

**B3.** (A quadratic radial polarization, and a potential check.) $\mathbf P = kr^2\hat{\mathbf r}$ in a sphere of radius $R$, no free charge.
(a) $\rho_b$, $\sigma_b$, total. (b) $\mathbf E$ inside and outside. (c) $V(0)$ from $\mathbf E$. (d) $V(0)$ directly from $\frac1{4\pi\varepsilon_0}\left[\oint\sigma_b\,da/R + \int\rho_b\,d\tau/r\right]$; show agreement.
(e) $R = 1$ cm, $P(R) = 1\ \mu$C/m²: $V(0)$, $E$ just inside $R$, $\rho_b(R)$. Why is $E = 0$ outside although $\rho_b$ and $\sigma_b$ are not zero?

**B4.** (A point charge inside a dielectric sphere.) $q = 1$ nC sits at the centre of a linear dielectric sphere, $\varepsilon_r = 3$, $R = 5$ cm, vacuum outside.
(a) $\mathbf D$, $\mathbf E$, $\mathbf P$ everywhere. (b) The bound charge on $r = R$ and at the centre (use $\nabla\cdot(\hat{\mathbf r}/r^2) = 4\pi\delta^3$); show the total is zero and that the net charge at the centre is $q/\varepsilon_r$. (c) $\sigma_b(R)$, $E$ at $r = 2$ cm, and $E$ just inside and just outside $R$; check the jump against $\sigma_b/\varepsilon_0$. (d) Compare the centre's bound charge with (3.19).

**B5.** (Two layers.) A capacitor of plate area $20$ cm² is filled by layer 1 against the positive plate ($1$ mm, $\varepsilon_r = 2$) and layer 2 against the negative plate ($2$ mm, $\varepsilon_r = 5$), at $V = 200$ V.
(a) $D$, $E_1$, $E_2$; check $\sum E\,d = V$. (b) $C$ and $Q$. (c) $P$ in each layer and $\sigma_b$ on all four faces; check they sum to zero. (d) Show the interface charge equals $(\mathbf P_1 - \mathbf P_2)\cdot\hat{\mathbf n}_{12}$ and also $\varepsilon_0(E_2 - E_1)$. There is no free charge there and no bulk $\rho_b$ in either layer: what produces it?

**B6.** (Which energy?) The series half-filled capacitor of Example 4.3(a), $\varepsilon_r = 4$, slab of thickness $d/2$ against one plate, free charge density $D$.
(a) Per unit area, in units of $D^2d/\varepsilon_0$: $W = \frac12\int\mathbf D\cdot\mathbf E$, $\frac{\varepsilon_0}{2}\int E^2$, $\frac12\int\mathbf P\cdot\mathbf E$.
(b) Show $\frac12\int\mathbf P\cdot\mathbf E/W = (\varepsilon_r - 1)/\varepsilon_r(\varepsilon_r + 1) = 0.15$.
(c) Which of the first two equals $\frac12QV$? (d) Interpret $\frac12\mathbf P\cdot\mathbf E$ in terms of $\mathbf p = \alpha\mathbf E$ (3.24).

**B7.** (A cylinder in a uniform field; code.) A long dielectric cylinder, radius $R$, permittivity $\varepsilon_{\text{in}}$, in a medium $\varepsilon_{\text{out}}$ and a transverse uniform field $E_0\hat{\mathbf x}$.
(a) A cylinder polarized uniformly across its axis, $\mathbf P = P\hat{\mathbf x}$: show $\mathbf E_{\text{in}} = -\mathbf P/2\varepsilon_0$ with two offset cylinders of charge density $\pm\rho$. Compare with the sphere.
(b) Use EM-3 (3.24) with only $n = 1$: $V_{\text{in}} = -As\cos\phi$, $V_{\text{out}} = -E_0s\cos\phi + B\cos\phi/s$. Write both matching conditions and solve: $E_{\text{in}} = 2\varepsilon_{\text{out}}E_0/(\varepsilon_{\text{in}} + \varepsilon_{\text{out}})$.
(c) Rederive it for $\varepsilon_{\text{out}} = \varepsilon_0$ by self-consistency, using (a).
(d) In project Part 4, `cylinder_in_field(4.0, 1.0, L=L, n=n)` with `exact_bc=False` and `True`, for $(L, n) = (4R, 121), (8R, 241), (16R, 481)$ (the same $h = R/15$). Show (i) the staircase error is the same for all $L$ with exact boundary data; (ii) the difference between the two boundary conditions times $L^2$ is constant. Predict the truncation error at $L = 32R$ before computing.

**B8.** (The slab, both ways.) Plates $w = 10$ cm wide and $\ell = 10$ cm long, gap $d = 1$ cm; a slab of $\varepsilon_r = 4$ fills the gap to depth $x$.
(a) $C(x)$. (b) A battery holds $V = 1$ kV while the slab is pulled in from $x_1 = 2$ cm to $x_2 = 6$ cm: battery's work, change in field energy, mechanical work, force $F$.
(c) Repeat after disconnecting the battery at $x_1$: $W(x)$, the voltage and force at $x_2$, and the work done by the field. Why is it less than in (b)?
(d) Why is $F$ at $x_1$ the same in both?
(e) Integrate the force density $(\mathbf P\cdot\nabla)\mathbf E = \frac12\varepsilon_0\chi_e\nabla E^2$ over the slab to recover (3.25). Where along the slab is the force exerted?

## C. Challenge (do at least one)

**C1.** (Clausius–Mossotti: where it works and where it fails.)
(a) Liquid argon ($1.40$ g/cm³, molar mass $39.948$, $\alpha/4\pi\varepsilon_0 = 1.64$ Å$^3$): $N$, then $\varepsilon_r$ from (3.27); compare with the dilute $\chi_e = N\alpha/\varepsilon_0$ and the tabulated $1.5$.
(b) Water ($1.00$ g/cm³, $18.015$), electronic $\alpha/4\pi\varepsilon_0 = 1.45$ Å$^3$: $\varepsilon_r$ from (3.27) against $n^2 = 1.333^2$.
(c) Add Debye's $\alpha_{\text{or}} = p^2/3k_BT$ ($p = 6.17\times10^{-30}$ C m, $T = 293$ K). Compute $N\alpha/3\varepsilon_0$ and the $\varepsilon_r$ that (3.27) now returns. What went wrong?
(d) Onsager's equation (quoted, not derived): $\dfrac{(\varepsilon - n^2)(2\varepsilon + n^2)}{\varepsilon(n^2 + 2)^2} = \dfrac{Np^2}{9\varepsilon_0k_BT}$. Show it reduces to Debye's $\varepsilon - 1 = Np^2/3\varepsilon_0k_BT$ for a dilute gas, solve it for water, and find the factor $g$ by which $p^2$ would have to be multiplied to reach $80.1$.

**C2.** (A piezoelectric plate.) A quartz plate, thickness $t = 1$ mm, electrodes of area $1$ cm² on its faces, $\varepsilon_r = 4.5$, $d_{11} = 2.3$ pC/N, is stressed by $1$ MPa along its thickness; $P_{\text{pz}} = d_{11}\sigma$.
(a) $\sigma_b$, the charge $Q$ on each face, and the force applied. (b) The open-circuit voltage. (c) The charge through a wire if the electrodes are shorted. (d) The stored energy two ways. (e) The converse effect $S = d_{11}E$: the thickness change at $100$ V.

**C3.** (The averaging lemma, second proof.) Prove Lemma 3.6 from $\mathbf E = -\nabla V$.
(a) Show $\langle\mathbf E\rangle_{\text{ball}} = -\frac{1}{\mathcal V}\oint V\hat{\mathbf n}\,da$, taking care at the point charge.
(b) For $q$ at $\mathbf r'$ with $r' < R_s$, expand $V$ on the sphere with EM-3 (3.27) and derive (3.13).
(c) Repeat for $r' > R_s$. (d) Why does only $\ell = 1$ survive?
