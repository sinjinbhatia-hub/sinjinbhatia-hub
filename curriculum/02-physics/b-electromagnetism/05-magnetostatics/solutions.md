# EM-5 Solutions

Equation numbers refer to [notes.md](notes.md). Every numerical value below is printed by `project.py` (Part 9, or the Part named in the text); run it to reproduce them.

## A. Warm-up

**A1.** (a) $K = 5\ \text{MeV} = 5\times10^6\times1.602\times10^{-19} = 8.011\times10^{-13}$ J. Then
$$v = \sqrt{2K/m} = \sqrt{2(8.011\times10^{-13})/6.645\times10^{-27}} = 1.553\times10^7\ \text{m/s}.$$
By (3.2), with $q = 2e$,
$$r = \frac{mv}{qB} = \frac{(6.645\times10^{-27})(1.553\times10^7)}{2(1.602\times10^{-19})(0.8)} = 0.4025\ \text{m},\qquad f = \frac{qB}{2\pi m} = 6.140\ \text{MHz}\ (T = 0.1629\ \mu\text{s}).$$
The relativistic correction is $K/mc^2 = 1.3\times10^{-3}$, so the non-relativistic formulas are good to about $0.1\%$.

(b) The radius is $r = p/qB$, so the same $r$ and $B$ need $p/q$ equal: $p_p/e = p_\alpha/2e$, i.e. $p_p = p_\alpha/2$. Then
$$K_p = \frac{p_p^2}{2m_p} = \frac{p_\alpha^2}{8m_p} = \frac{m_\alpha}{4m_p}K_\alpha = \frac{6.645}{4(1.6726)}(5\ \text{MeV}) = 4.966\ \text{MeV}.$$
(c) The force $q\mathbf v\times\mathbf B$ is perpendicular to $\mathbf v$, so its power $\mathbf F\cdot\mathbf v$ is zero (Theorem 3.1); it turns $\mathbf v$ and cannot change $|\mathbf v|$.

**A2.** (a) $B = \dfrac{\mu_0I}{2\pi s} = \dfrac{2\times10^{-7}(30)}{0.05} = 1.200\times10^{-4}$ T $= 120.0\ \mu$T.

(b) The segment has half-length $0.10$ m, so at the midpoint-perpendicular $\sin\theta_2 = -\sin\theta_1 = 0.10/\sqrt{0.05^2 + 0.10^2} = 0.8944$. By (3.10),
$$B = \frac{\mu_0I}{4\pi s}(2\times0.8944) = \frac{10^{-7}(30)}{0.05}(1.7889) = 107.33\ \mu\text{T},$$
which is $0.8944$ of the long-wire value (the fraction is $\sin\theta$ itself, because $\mu_0I/4\pi s\cdot2\sin\theta = (\mu_0I/2\pi s)\sin\theta$).

(c) $F/L = \mu_0I_1I_2/2\pi d = 2\times10^{-7}(10)(25)/0.02 = 2.50\times10^{-3}$ N/m. Direction: put wire 1 on the $z$-axis with current $+z$; at wire 2 (at $x = d$) $\mathbf B_1 = B\hat{\mathbf y}$. Wire 2 carries $-z$, so by (3.6) the force per length is $(-I_2\hat{\mathbf z})\times B\hat{\mathbf y} = +I_2B\,\hat{\mathbf x}$: away from wire 1. **Repulsive** (antiparallel currents repel).

**A3.** (a) Centre: $B = \mu_0NI/2R = 4\pi\times10^{-7}(75)/0.08 = 1.1781$ mT. On axis, (3.13) with $I\to NI$ is $B(z) = B(0)(1 + z^2/R^2)^{-3/2}$. At $z = 3$ cm, $z/R = 0.75$, $1 + 0.5625 = 1.5625 = 1.25^2$, so the factor is $1.25^{-3} = 0.512$ and $B = 0.6032$ mT.
(b) $(1 + u^2)^{-3/2} = \tfrac12$ gives $1 + u^2 = 2^{2/3} = 1.5874$, $u = 0.7664$, so $z = 0.7664R = 3.066$ cm.

**A4.** Ampère on a circle of radius $s$ inside the doughnut encloses $N$ wires each carrying $I$ through the circle's disk: $B\,2\pi s = \mu_0NI$, so $B = \mu_0NI/2\pi s = 1.92\times10^{-4}/s$ (SI). At $s = 6, 8, 10$ cm: $3.20, 2.40, 1.92$ mT. In the hole ($s < 6$ cm) a circle encloses no wire, so $2\pi sB = 0$. Outside ($s > 10$ cm) the disk of the circle is pierced by the winding twice: $N$ turns of current go through one way (the hole side) and $N$ the other way (the outer side), so $I_{\rm enc} = NI - NI = 0$. Both cases rely on $\mathbf B = B(s)\hat{\boldsymbol\phi}$ (the symmetry assumption of §3.6).

**A5.** (a) Use $d\mathbf l'\times\hat{\boldsymbol{\mathscr r}}$ with $d\mathbf l' = \hat{\mathbf x}$. At $(0,0,1)$: $\hat{\boldsymbol{\mathscr r}} = \hat{\mathbf z}$ and $\hat{\mathbf x}\times\hat{\mathbf z} = -\hat{\mathbf y}$. At $(0,1,0)$: $\hat{\mathbf x}\times\hat{\mathbf y} = +\hat{\mathbf z}$. (Thumb along $+x$; the fingers curl from $+y$ toward $+z$ and then toward $-y$.) The project's integrator returns $-2\times10^{-7}\hat{\mathbf y}$ T and $+2\times10^{-7}\hat{\mathbf z}$ T for $1$ A.
(b) Clockwise from $+z$: the thumb points $-\hat{\mathbf z}$, so $\mathbf m = I\mathbf a$ points along $-\hat{\mathbf z}$, and on the axis $\mathbf B$ points along $\mathbf m$ ($-\hat{\mathbf z}$; the integrator gives $B_z = -6.29\ \mu$T at the centre of a $10$ cm, $1$ A loop). The loops **repel**: nearest sides carry antiparallel currents; equivalently $\mathbf m_2\cdot\mathbf B_1 < 0$ so $U = -\mathbf m_2\cdot\mathbf B_1 > 0$ falls as the loops separate. The force function of B8 gives $F_z > 0$ (away from the source) on a downward-pointing loop above a counterclockwise ring.
(c) $[\mu_0I/s]$ = (N/A$^2$)(A/m) = N/(A m), and from $F = ILB$, $[B]$ = N/(A m) = T. $[\mathbf A]$ = T m, so $[\oint\mathbf A\cdot d\mathbf l]$ = T m$^2$ = Wb. From $\varepsilon_0$ in C$^2$/(N m$^2$) = A$^2$s$^2$/(N m$^2$), $\mu_0\varepsilon_0$ has units s$^2$/m$^2$, so $1/\sqrt{\mu_0\varepsilon_0}$ is m/s; its value is $2.997925\times10^8$ m/s $= c$.

## B. Core

**B1.** (a) Orient $\mathcal S_1$ and $\mathcal S_2$ both by the right-hand rule around $\mathcal P$. The union $\mathcal S_1 - \mathcal S_2$ (with $\mathcal S_2$ reversed) is a closed surface bounding a volume $\mathcal V$; its outward normal agrees with the right-hand orientation of whichever surface lies outside and is opposite to the other's (in (b) it agrees with the hemisphere $\mathcal S_2$ and opposes the disk $\mathcal S_1$; if the surfaces cross, split the argument at the crossing). The overall sign does not matter, because the closed-surface integral is zero: By the divergence theorem (EM-1 (3.18)) and (3.8),
$$\int_{\mathcal S_1}\mathbf J\cdot d\mathbf a - \int_{\mathcal S_2}\mathbf J\cdot d\mathbf a = \pm\oint\mathbf J\cdot d\mathbf a = \pm\int_{\mathcal V}\nabla\cdot\mathbf J\,d\tau = 0.\ \blacksquare$$
(b) $I = J_0\pi a^2$. (i) On the flat disk $\hat{\mathbf n} = \hat{\mathbf z}$ and $\mathbf J\cdot\hat{\mathbf n} = J_0$ for $s < a$, $0$ for $a < s < b$: $I_{\rm enc} = J_0\pi a^2$. (ii) On the hemisphere $\hat{\mathbf n} = \hat{\mathbf r}$, $d\mathbf a = b^2\sin\theta\,d\theta\,d\phi$ and the wire occupies the cap $s = b\sin\theta < a$, i.e. $\theta < \theta_a = \arcsin(a/b)$, where $\mathbf J\cdot\hat{\mathbf n} = J_0\cos\theta$:
$$I_{\rm enc} = \int_0^{2\pi}d\phi\int_0^{\theta_a}J_0\cos\theta\,b^2\sin\theta\,d\theta = 2\pi J_0b^2\frac{\sin^2\theta_a}{2} = \pi J_0a^2.$$
The two agree, as (a) requires.
(c) Step (a) uses $\nabla\cdot\mathbf J = 0$. In a charging capacitor charge piles up on a plate, $\nabla\cdot\mathbf J = -\partial_t\rho\ne0$ there. A flat disk around the wire catches the wire current $I$; a bulging surface that passes between the plates catches no conduction current (the gap has $\mathbf J = 0$). The two surfaces share a boundary but differ in $I_{\rm enc}$ by $dQ/dt$, so (3.19) cannot hold for both. Ampère's law needs a repair term in time-dependent situations (EM-6's displacement current).

**B2.** (a) $I = \int_0^aJ_0\frac sa\,2\pi s\,ds = \frac{2\pi J_0}{a}\frac{a^3}{3} = \frac{2\pi}{3}J_0a^2$.
(b) $\mathbf B = B_\phi(s)\hat{\boldsymbol\phi}$ by the argument of §3.6. A circle of radius $s$ gives $2\pi sB_\phi = \mu_0I_{\rm enc}(s)$. For $s < a$, $I_{\rm enc} = 2\pi J_0s^3/3a$, so
$$B_\phi = \frac{\mu_0J_0s^2}{3a}\ (s < a),\qquad B_\phi = \frac{\mu_0J_0a^2}{3s} = \frac{\mu_0I}{2\pi s}\ (s > a).$$
Both equal $\mu_0J_0a/3$ at $s = a$.
(c) In cylindrical coordinates with $\mathbf B = B_\phi(s)\hat{\boldsymbol\phi}$ only, $(\nabla\times\mathbf B)_s = -\partial_zB_\phi = 0$ and $(\nabla\times\mathbf B)_\phi = \partial_zB_s - \partial_sB_z = 0$, and
$$(\nabla\times\mathbf B)_z = \frac1s\frac{d}{ds}(sB_\phi) = \frac1s\frac{d}{ds}\frac{\mu_0J_0s^3}{3a} = \frac{\mu_0J_0s}{a} = \mu_0J_z.\ \checkmark$$
(d) $J_0 = 3I/2\pi a^2 = 3(10)/(2\pi\times10^{-4}) = 4.775\times10^4$ A/m$^2$. $B(a) = \mu_0I/2\pi a = 200.0\ \mu$T. Inside, $B\propto s^2$, so $B(a/2) = B(a)/4 = 50.0\ \mu$T; uniform $\mathbf J$ would give $B\propto s$ and $100.0\ \mu$T. Half the current is enclosed where $(s/a)^3 = \tfrac12$: $s = 2^{-1/3}a = 0.7937a = 7.937$ mm (uniform $\mathbf J$: $a/\sqrt2 = 0.707a$; here more current sits at large $s$).

**B3.** Let $f(u) = (R^2 + u^2)^{-3/2}$ (even in $u$), so by (3.13) $B_z(z) = \tfrac12\mu_0IR^2[f(z - d/2) - f(z + d/2)]$.
(a) Replacing $z\to-z$ gives $f(z + d/2) - f(z - d/2) = -[\,f(z - d/2) - f(z + d/2)]$ by evenness of $f$: $B_z$ is odd, so $B_z(0) = 0$. Differentiating, $B_z'(0) = \tfrac12\mu_0IR^2[f'(-d/2) - f'(d/2)] = -\mu_0IR^2f'(d/2)$ ($f'$ is odd). With $f'(u) = -3u(R^2 + u^2)^{-5/2}$,
$$B_z'(0) = \tfrac32\mu_0IR^2\,d\,(R^2 + d^2/4)^{-5/2} = \frac{48\mu_0IR^2d}{(4R^2 + d^2)^{5/2}}.$$
(b) $f'' = 3(4u^2 - R^2)(R^2 + u^2)^{-7/2}$ (§3.4), and differentiating once more,
$$f''' = 3\left[8u(R^2 + u^2)^{-7/2} - 7u(4u^2 - R^2)(R^2 + u^2)^{-9/2}\right] = 15u(3R^2 - 4u^2)(R^2 + u^2)^{-9/2}.$$
$f'''$ is odd, so $B_z'''(0) = -\mu_0IR^2f'''(d/2)$, which vanishes (for $d\ne0$) iff $4(d/2)^2 = 3R^2$, i.e. $d = \sqrt3R$. (The project confirms this symbolically.)
(c) At $d = \sqrt3R$: $4R^2 + d^2 = 7R^2$, so $B_z'(0) = 48\mu_0IR^2\sqrt3R/(7R^2)^{5/2} = (48\sqrt3/7^{5/2})\mu_0I/R^2 = 0.6413\,\mu_0I/R^2$.
(d) For an axisymmetric field, $\nabla\cdot\mathbf B = \frac1s\partial_s(sB_s) + \partial_zB_z = 0$. Near the axis $\partial_zB_z\approx B_z'(0,z)$ independent of $s$, so $\partial_s(sB_s) = -sB_z'$, and integrating with $B_s(0) = 0$ (regularity) gives $sB_s = -s^2B_z'/2$: $B_s = -\tfrac12sB_z'$.
(e) $\mu_0I/R^2 = 1.2566\times10^{-6}(1000)/0.0025 = 0.5027$ T/m, so $B_z'(0) = 0.6413\times0.5027 = 0.3223$ T/m (more precisely $0.32235$; the integrator's central difference gives $0.3224$ T/m). Departure of $B_z$ from $B_z'(0)z$: $-0.20\%$ at $z = R/4$ and $-3.00\%$ at $R/2$; the ratio is $15\approx2^4$, so the leading correction is $\propto z^5$ (the $z^3$ term is gone). At $d = R$ (same opposite currents) the $z^3$ term survives and the departure at $R/4$ is $-6.6\%$, about $30$ times worse. The radial check: $B_s(0.1R) = -805.9\ \mu$T from the integrator, $-\tfrac12(0.1R)B_z' = -805.9\ \mu$T from (d).

**B4.** (a) $\mathbf A_1 = \tfrac12B_0(-y, x, 0)$: $(\nabla\times\mathbf A_1)_z = \partial_xA_y - \partial_yA_x = \tfrac12B_0 + \tfrac12B_0 = B_0$; the $x$- and $y$-components involve only $\partial_zA_y,\partial_zA_x,\partial_yA_z,\partial_xA_z$, all zero; $\nabla\cdot\mathbf A_1 = \partial_x(-\tfrac12B_0y) + \partial_y(\tfrac12B_0x) = 0$. $\mathbf A_2 = B_0(0, x, 0)$: $(\nabla\times\mathbf A_2)_z = \partial_x(B_0x) = B_0$, other components $0$; $\nabla\cdot\mathbf A_2 = \partial_y(B_0x) = 0$.
(b) $\mathbf A_2 - \mathbf A_1 = (\tfrac12B_0y, \tfrac12B_0x, 0) = \nabla\lambda$ with $\lambda = \tfrac12B_0xy$. $\nabla^2\lambda = \partial_x^2\lambda + \partial_y^2\lambda = 0$. So two different vector potentials, both in Coulomb gauge, describe the same field: the Coulomb gauge fixes $\mathbf A$ only up to $\nabla\lambda$ with $\nabla^2\lambda = 0$ (§3.7, step 3 solves Poisson for the *difference* between gauges). Boundary conditions (decay at infinity) fix it for localized sources; a uniform field extends to infinity, so they cannot.
(c) Take the square $(0,0)\to(a,0)\to(a,a)\to(0,a)\to(0,0)$, $a = 0.1$ m, $B_0 = 0.5$ T. For $\mathbf A_2$: the bottom ($A_x = 0$) and top ($A_x = 0$) contribute $0$; the left edge has $A_y = 0$; the right edge has $A_y = B_0a$ over length $a$: total $B_0a^2$. For $\mathbf A_1$: bottom: $A_x = -\tfrac12B_0y = 0$, so $0$; right edge ($x = a$, $dl = dy$): $A_y = \tfrac12B_0a$, giving $\tfrac12B_0a^2$; top ($y = a$, $x$ from $a$ to $0$): $A_x = -\tfrac12B_0a$, $\int A_xdx = (-\tfrac12B_0a)(0 - a) = \tfrac12B_0a^2$; left edge: $A_y = 0$. Total $B_0a^2$. Both give $5.000\times10^{-3}$ Wb $= B_0a^2$ = the flux (Theorem 3.6), as the gauge-invariance argument demands.
(d) $\mathbf A_1 = \tfrac12\mathbf B\times\mathbf r$ is azimuthal, $A_\phi = B_0s/2$, and respects the rotation symmetry; $\mathbf A_2$ depends on $x$ only and is better when the problem is translation-invariant in $y$ (the Landau gauge of QM-4).

**B5.** (a) With $\mathbf A = A_x(z)\hat{\mathbf x}$, $\nabla\times\mathbf A = (0, A_x'(z), 0)$. Matching $\mathbf B = -\tfrac12\mu_0K\,\mathrm{sgn}(z)\hat{\mathbf y}$ gives $A_x' = -\tfrac12\mu_0K\,\mathrm{sgn}(z)$, so
$$\mathbf A = -\tfrac12\mu_0K|z|\,\hat{\mathbf x} + \text{const}.$$
(The constant is a gauge choice, $\nabla(Cx)$; we take $0$.)
(b) $\nabla\cdot\mathbf A = \partial_xA_x = 0$. $\nabla^2A_x = A_x'' = -\tfrac12\mu_0K\cdot2\delta(z) = -\mu_0K\delta(z)$, and the sheet's volume current density is $J_x = K\delta(z)$, so $\nabla^2\mathbf A = -\mu_0\mathbf J$. $A_x(0^\pm) = 0$: continuous. $A_x'(0^+) - A_x'(0^-) = -\tfrac12\mu_0K - \tfrac12\mu_0K = -\mu_0K$, which is (3.29).
(c) $|A_x|$ grows linearly: the sheet is infinite, so the integral (3.24) $\int\mathbf J/\mathscr r\,d\tau'$ diverges, as for the infinite wire (log). EM-2 has the same structure for the infinite charged plane, $V = -\sigma|z|/2\varepsilon_0$, which also cannot vanish at infinity, so its constant is fixed by a convention instead.
(d) $|\mathbf B| = \mu_0K/2 = 4\pi\times10^{-7}(100)/2 = 62.83\ \mu$T on each side, in opposite directions; $\partial_zA_x$ jumps by $-\mu_0K = -125.66\ \mu$T.

**B6.** (a) $\mathbf M$ is uniform so $\mathbf J_b = \nabla\times\mathbf M = 0$. $\mathbf K_b = \mathbf M\times\hat{\mathbf n}$: on the curved side $\hat{\mathbf n} = \hat{\mathbf s}$ and $\hat{\mathbf z}\times\hat{\mathbf s} = \hat{\boldsymbol\phi}$, so $\mathbf K_b = M\hat{\boldsymbol\phi}$ ($10^6$ A/m); on the end faces $\hat{\mathbf n} = \pm\hat{\mathbf z}\parallel\mathbf M$ so $\mathbf K_b = 0$. The magnet is a solenoid with $nI = M$.
(b) (3.14) with $nI = M$, $L = 20$ mm, $R = 5$ mm. Centre ($z = 0$): the bracket is $2(L/2)/\sqrt{R^2 + L^2/4}$, so
$$B = \mu_0M\frac{L/2}{\sqrt{R^2 + L^2/4}} = 1.2566\ \text{T}\times\frac{0.010}{\sqrt{2.5\times10^{-5} + 10^{-4}}} = 1.2566\times0.8944 = 1.1240\ \text{T}.$$
End-face centre ($z = L/2$): the second term vanishes and the first is $L/\sqrt{R^2 + L^2}$:
$$B = \tfrac12\mu_0M\frac{L}{\sqrt{R^2 + L^2}} = 0.6283\times\frac{0.020}{0.020616} = 0.6096\ \text{T}.$$
(The ring-stack integrator gives $1.1243$ T and $0.6096$ T.)
(c) $\mathbf H = \mathbf B/\mu_0 - \mathbf M$. At the centre $B/\mu_0 = 8.944\times10^5$ A/m, so $H = -1.0557\times10^5$ A/m $= -0.106M$. At the end face, $B^\perp$ is continuous, $B/\mu_0 = 4.8507\times10^5$ A/m on both sides, and $\mathbf M$ is $M$ inside and $0$ outside: $H_{\rm inside} = 4.8507\times10^5 - 10^6 = -5.1493\times10^5$ A/m and $H_{\rm outside} = +4.8507\times10^5$ A/m. The jump is $4.8507\times10^5 + 5.1493\times10^5 = 10^6 = M$, as (3.42) requires.
(d) Rod ($H = 0$) $>$ this cylinder ($-0.106M$) $>$ sphere ($-0.333M$) $>$ thin disk ($-M$). $\mathbf H$ is the field of the fictitious surface charges $\pm M$ on the end faces (§3.15); the farther away those faces are from the point of interest (relative to its size), the weaker the opposing field. A rod puts the faces at infinity; a disk puts them at zero distance.

**B7.** (a) $\mathbf H$ obeys $\oint\mathbf H\cdot d\mathbf l = I_{f,\rm enc}$, so $H_\phi = Is/2\pi a^2$ for $s < a$ and $H_\phi = I/2\pi s$ for $s > a$. Inside, $\mathbf B = \mu_0\mu_r\mathbf H$ and $\mathbf M = \chi_m\mathbf H = (\mu_r - 1)\mathbf H$ (both along $\hat{\boldsymbol\phi}$); outside, $\mathbf M = 0$ and $\mathbf B = \mu_0\mathbf H$. At $s = a$: $H = 20/(2\pi\times0.005) = 636.62$ A/m, $B = 200\mu_0(636.62) = 0.1600$ T, $M = 199(636.62) = 1.2669\times10^5$ A/m.
(b) $\mathbf J_b = \nabla\times\mathbf M = \chi_m\nabla\times\mathbf H = \chi_m\mathbf J_f$ (3.43): $\chi_mI/\pi a^2 = 199(2.5465\times10^5) = 5.0675\times10^7$ A/m$^2$ along $+\hat{\mathbf z}$ — $199$ times the free current density. $\mathbf K_b = \mathbf M\times\hat{\mathbf s} = M(a)\,\hat{\boldsymbol\phi}\times\hat{\mathbf s} = -M(a)\hat{\mathbf z} = -1.2669\times10^5$ A/m. Totals: volume $\chi_mI = +3980$ A; surface $-M(a)\,2\pi a = -\chi_mH(a)\,2\pi a = -\chi_mI = -3980$ A. They cancel. They must: outside the wire the Ampère circulation of $\mathbf B$ equals $\mu_0$ times *all* current, and since $\mathbf M = 0$ there, $\mathbf B = \mu_0\mathbf H = \mu_0I/2\pi s$ depends on $I$ alone, so the net bound current is zero.
(c) At $s = 2a = 1$ cm: $B = \mu_0I/2\pi s = 0.400$ mT, the vacuum value; the permeability changes $B$ only inside the wire (at $s = a$, $0.1600$ T instead of $0.8$ mT).

**B8.** A reference implementation (it is the `small_loop_force_torque` of Part 9):

```python
def small_loop_force_torque(src, c, normal, side, I_s):
    n = np.asarray(normal, float); n /= np.linalg.norm(n)
    e1 = np.cross(n, [0, 1.0, 0]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)                          # e1 x e2 = n: counterclockwise about n
    h = side / 2; c = np.asarray(c, float)
    P = np.array([c - h*e1 - h*e2, c + h*e1 - h*e2, c + h*e1 + h*e2, c - h*e1 + h*e2])
    dl = np.roll(P, -1, axis=0) - P               # four edge vectors
    mid = P + dl / 2
    f = np.cross(dl, biot_savart(mid, *src))      # dl x B at the midpoints
    return I_s * f.sum(0), I_s * np.cross(mid - c, f).sum(0), I_s * side**2 * n
```
For the stated case $\mathbf m = (2.000\times10^{-6},\,0,\,3.464\times10^{-6})$ A m$^2$ and $\mathbf B(\mathbf c) = (9.837\times10^{-7},\,0,\,4.2136\times10^{-6})$ T. The loop gives
$$\mathbf F = (9.362\times10^{-12},\ 0,\ -2.5414\times10^{-10})\ \text{N},\qquad\nabla(\mathbf m\cdot\mathbf B) = (9.366\times10^{-12},\ 0,\ -2.5414\times10^{-10})\ \text{N},$$
a relative difference of $1.5\times10^{-5}$, and $\mathbf N = (0, -5.0185\times10^{-12}, 0)$ N m versus $\mathbf m\times\mathbf B = (0, -5.0196\times10^{-12}, 0)$ N m, a relative difference of $2.0\times10^{-4}$. The $y$-force and the $x$ and $z$ torques are zero by symmetry ($y = 0$ is a mirror plane of ring and loop, and $\mathbf m,\mathbf B$ lie in it). Also $U = -\mathbf m\cdot\mathbf B = -1.656\times10^{-11}$ J. The force points down (toward the ring: $\mathbf m$ has a component along $\mathbf B$, which grows toward the wire) and outward in $x$ (inside the ring $B_z$ grows toward the wire).

Halving the side from $4$ to $2$ to $1$ mm, the force errors are $5.9\times10^{-5}$, $1.5\times10^{-5}$, $3.5\times10^{-6}$ and the torque errors $8.1\times10^{-4}$, $2.0\times10^{-4}$, $5.1\times10^{-5}$: a factor $\approx4$ each time, i.e. $O(\text{side}^2)$. Reason: the midpoint rule is exact when $\mathbf B$ varies linearly along an edge, and the first neglected term is the curvature of $\mathbf B$, of relative size $(\text{side}/R)^2$. The identities (3.34) and (3.35) themselves drop terms of the same order (they are the first terms of a multipole expansion in side$/R$).

## C. Challenge

**C1.** (a) With $\mathbf B = \nabla\times\mathbf A$ taken in the distributional sense (the same sense in which (3.33) contains a delta function; $\mathbf A\sim1/r^2$ is integrable at the origin), the curl theorem (§3.13, step 3) gives $\int_{r<R}\mathbf B\,d\tau = \oint\hat{\mathbf r}\times\mathbf A\,da$. With $\mathbf A = \frac{\mu_0}{4\pi}\mathbf m\times\hat{\mathbf r}/R^2$ on the sphere and $da = R^2d\Omega$,
$$\hat{\mathbf r}\times(\mathbf m\times\hat{\mathbf r}) = \mathbf m - \hat{\mathbf r}(\hat{\mathbf r}\cdot\mathbf m),\qquad\oint[\mathbf m - \hat{\mathbf r}(\hat{\mathbf r}\cdot\mathbf m)]\,d\Omega = 4\pi\mathbf m - \tfrac{4\pi}{3}\mathbf m = \tfrac{8\pi}{3}\mathbf m,$$
using $\oint\hat r_i\hat r_j\,d\Omega = \tfrac{4\pi}{3}\delta_{ij}$. So $\int_{r<R}\mathbf B\,d\tau = \frac{\mu_0}{4\pi}\frac{8\pi}{3}\mathbf m = \tfrac23\mu_0\mathbf m$, for every $R$.
(b) The non-delta part is $\frac{\mu_0}{4\pi}[3(\mathbf m\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf m]/r^3$. On the shell $\epsilon < r < R$ the radial integral gives $\int r^2dr/r^3 = \ln(R/\epsilon)$ times the angular integral $\oint(3\hat r_i\hat r_j - \delta_{ij})\,d\Omega\,m_j = (4\pi - 4\pi)m_i = 0$. Zero for every $\epsilon$, so the principal-value integral over the ball vanishes.
(c) The total over the ball is $\tfrac23\mu_0\mathbf m$ by (a) and the non-delta part contributes $0$ by (b), so the delta term must integrate to $\tfrac23\mu_0\mathbf m$: it is $\tfrac{2\mu_0}{3}\mathbf m\,\delta^3(\mathbf r)$, which is the coefficient in (3.33).
(d) The uniformly magnetized sphere has uniform $\mathbf B_{\rm in} = \tfrac23\mu_0\mathbf M$, so $\int B_z\,d\tau = \tfrac23\mu_0M\cdot\tfrac{4\pi}{3}R^3 = \tfrac23\mu_0m$ with $m = \tfrac{4\pi}{3}R^3M$; the integral over any larger ball is the same, because outside the sphere the dipole field has zero angular average. The project gives $3.5103\times10^{-6}$ T m$^3$ from the ring model against $3.5092\times10^{-6}$ for $\tfrac23\mu_0m$. For the electric dipole, $\mathbf E = -\nabla V$ and $\int\mathbf E\,d\tau = -\oint V\hat{\mathbf n}\,da = -\frac{1}{4\pi\varepsilon_0}\oint(\mathbf p\cdot\hat{\mathbf r})\hat{\mathbf r}\,d\Omega = -\frac{\mathbf p}{3\varepsilon_0}$. The sign is opposite relative to the dipole because the magnetic surface integrand has an extra "$+\mathbf m$" term (the cross product $\hat{\mathbf r}\times(\mathbf m\times\hat{\mathbf r})$) that the scalar one lacks: $\mathbf B$ points *along* $\mathbf m$ inside the loop, $\mathbf E$ points *against* $\mathbf p$ inside the dipole.

**C2.** (a) The flux through core and gap is the same and the areas are equal, so one $B$ serves both. In the core $H_c = \dfrac{B_s}{\mu_0\mu_r}\,\mathrm{artanh}\dfrac{B}{B_s}$; in the gap $H = B/\mu_0$. Ampère for $\mathbf H$ around the circuit, $H_cl + Hg = NI$:
$$\frac{B_sl}{\mu_0\mu_r}\,\mathrm{artanh}\frac{B}{B_s} + \frac{Bg}{\mu_0} = NI. \tag{$*$}$$
(b) The left side increases monotonically from $0$ to $\infty$ as $B$ goes from $0$ to $B_s$, so a root exists and is unique; `brentq` finds it on $(0, B_s)$. Results (nonlinear / linear (3.46)):

| $NI$ (A) | 1000 | 3000 | 5000 | 10000 |
|---|---|---|---|---|
| nonlinear $B$ (T) | 0.2707 | 0.8023 | 1.2878 | 1.6000 |
| linear $B$ (T) | 0.2711 | 0.8132 | 1.3553 | 2.7106 |

(c) For $B = 1.0$ T: $B_sl/\mu_0\mu_r = 1.6(0.636)/(1.2566\times10^{-3}) = 809.8$ A, $\mathrm{artanh}(0.625) = 0.7332$, so the core needs $593.7$ A; the gap needs $Bg/\mu_0 = 3183.1$ A. Total $3776.8$ A, against $3689.2$ A in the linear model ($2.4\%$ more). For $1.5$ T the figures are $6165.0$ versus $5533.8$ A ($11\%$ more).
(d) Since $B_{\rm core} = B_s\tanh(\cdot) < B_s$ and the gap carries the same $B$, no $NI$ can push the gap field above $B_s = 1.6$ T in this model; at $NI = 10^4$ it is $1.59998$ T, while the linear formula would claim $2.71$ T. A $2$ T magnet needs either a tapered pole whose gap area $A_g$ is smaller than the core area $A_c$ (flux concentration: $B_g = B_cA_c/A_g$ with $B_c < B_s$) or a higher-$B_s$ material such as a cobalt–iron alloy.

**C3.** (a) Saturating at $H = +3$ gives $M = 0.997\,M_s$ (a few hysterons have thresholds above $3$ and keep their random starting state). Setting $H = 0$ leaves $M = 0.448\,M_s$, the remanence of Your turn 3.
(b) The values of $M$ at the turning points $H_k = 3(-0.9)^k$ are $-0.988$ ($k = 1$, $H = -2.7$), $+0.470$ ($k = 10$, $H = +1.05$), $+0.052$ ($k = 20$, $H = +0.365$), $+0.001$ ($k = 40$) and $+0.001$ ($k = 80$). (At $k = 10$ the field is still $+1.05$, so $M$ is large; these are not remanences.) After the final $H = 0$: $M = +0.0011$, about $400$ times smaller than $0.448$. Plotted against $H$ the turning points form nested minor loops that shrink toward the origin (`figures/demag.png`).
(c) With factor $0.5$ and $24$ steps the final $M = +0.0090$: still tiny, but $8$ times larger than with $0.9$.
(d) A hysteron flips up when $H$ exceeds its up-threshold $u$ and down when $H$ falls below its down-threshold $d$. In an alternating sequence of shrinking amplitude, the state of a hysteron is set by the *last* excursion that reached one of its thresholds. If $|u| < |d|$ the up-threshold keeps being reached by later, smaller positive excursions, so the hysteron ends up; if $|d| < |u|$ it ends down. Writing $u = h_u + h_c$, $d = h_u - h_c$, the condition $|u| < |d|$ is $h_u < 0$ (check the four sign cases, including both thresholds on one side of zero, which the final step to $H = 0$ settles the same way). So the final state is $-\mathrm{sign}(h_u)$: it depends on the hysteron's *bias*, not on the history, and with $h_u$ distributed symmetrically about zero the ensemble sums to $\approx0$. The project checks this: the prediction $-\overline{\mathrm{sign}(h_u)} = +0.0033$ (the finite-sample imbalance of $20{,}000$ hysterons), and the factor $0.99$ (800 steps) gives $+0.0035$. The factor $0.9$ gives $0.0011$ and $0.5$ gives $0.0090$: the deviations come from hysterons whose two thresholds differ in magnitude by less than one step, whose fate then depends on the parity of the step at which the sequence ends. Coarser steps catch more of them (that is the practical reason to decrease the amplitude slowly when demagnetizing a tool or a tape head).
