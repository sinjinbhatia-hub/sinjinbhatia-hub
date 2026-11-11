# CM-5 Solutions

Every number below is printed by `project.py` (Parts 9 and 10); every symbolic step marked "SymPy" is checked there.

## A. Warm-up

**A1.**
(a) $\mu/m_e = m_p/(m_p + m_e) = 1836.15267/1837.15267 = 0.99945568$. Ground state: $-13.605693\times0.99945568 = -13.5983$ eV. The finite proton mass raises every level by 5.4 parts in $10^4$.

(b) $\mu_D/m_e = 3670.48297/3671.48297 = 0.99972763$. Since $\lambda \propto 1/\mu$, $\lambda_D = \lambda_H\,\mu_H/\mu_D$ and
$$
\Delta\lambda = \lambda_H\Big(1 - \frac{\mu_H}{\mu_D}\Big) = 656.28\text{ nm}\times\Big(1 - \frac{0.99945568}{0.99972763}\Big) = 0.1785\text{ nm},
$$
toward the *blue* (deuterium's heavier nucleus gives a larger $\mu$, deeper levels, shorter wavelength). A shift of this size in the Balmer lines is how deuterium was identified spectroscopically in 1931–32 (Urey, Brickwedde and Murphy).

(c) Equal masses: $\mu = m_e/2$, so $E_1 = -13.605693/2 = -6.8028$ eV.

(d) $\mu/m_{\rm Moon} = m_\oplus/(m_\oplus + m_{\rm Moon}) = 5.9722/(5.9722 + 0.07342) = 0.98786$. By (3.2) each body sits at $(m_{\rm other}/M)$ times the separation from the barycenter, on the side away from the other body, so Earth's center is $(m_{\rm Moon}/M)\times384{,}400\text{ km} = 384{,}400\times0.07342/6.04562 = 4668$ km from it. That is less than 6371 km: the barycenter is inside the Earth, about 1700 km below the surface.

**A2.** Compare with (3.12): $c = 2$, $e = 0.5$, an ellipse. By (3.13): $r_p = c/(1 + e) = 1.3333$, $r_a = c/(1 - e) = 4.0000$, $a = c/(1 - e^2) = 2/0.75 = 2.6667$, $b = a\sqrt{1 - e^2} = 2.6667\times0.86603 = 2.3094$. From $c = \ell^2/\mu k$: $\ell = \sqrt2 = 1.4142$. From (3.15): $E = -k/2a = -0.1875$. Check (3.14): $1 + 2E\ell^2/\mu k^2 = 1 + 2(-0.1875)(2) = 0.25$, so $e = 0.5$ ✓. Period (3.16): $T = 2\pi a^{3/2} = 2\pi\times4.3546 = 27.36$. Speeds: $v_p = \ell/\mu r_p = 1.4142/1.3333 = 1.0607$, $v_a = 1.4142/4 = 0.3536$. Vis-viva at perihelion: $v_p^2 = 2/r_p - 1/a = 1.5 - 0.375 = 1.125$, $v_p = 1.0607$ ✓.

**A3.** $\sqrt{GM_\odot/a} = \sqrt{1.32712\times10^{20}/5.7909\times10^{10}} = 47.87$ km/s, and $\sqrt{(1 + e)/(1 - e)} = \sqrt{1.20563/0.79437} = 1.23196$. So $v_p = 47.87\times1.23196 = 58.98$ km/s and $v_a = 47.87/1.23196 = 38.86$ km/s. Ratio $v_p/v_a = 1.5177 = (1 + e)/(1 - e) = r_a/r_p$ ✓ (angular momentum: $r_pv_p = r_av_a$).

**A4.** By (3.16), $T \propto [G(M_\odot + m_J)]^{-1/2}$. Using $GM_\odot$ alone overestimates $T$ by the factor $\sqrt{1 + q}$, $q = GM_J/GM_\odot = 9.5459\times10^{-4}$:
$$
\Delta T = T\big(\sqrt{1 + q} - 1\big) = 4332.59\text{ d}\times4.772\times10^{-4} = 2.07\text{ d} \qquad(\approx \tfrac12qT;\ \tfrac12q = 4.773\times10^{-4}).
$$
Kepler, with periods good to perhaps a day, could not have seen it; any modern ephemeris must include it.

**A5.**
(a) $\int(R^2/4)\,d\Omega = (R^2/4)\cdot4\pi = \pi R^2$, the geometric shadow.

(b) With $d\Omega = 2\pi\sin\theta\,d\theta$ and $s = \sin(\theta/2)$: $ds = \tfrac12\cos(\theta/2)\,d\theta$, so $\sin\theta\,d\theta = 2s\cos(\theta/2)\,d\theta = 4s\,ds$. Then
$$
\sigma(\theta > \theta_0) = 2\pi\Big(\frac{k}{4E}\Big)^2\int_{s_0}^1\frac{4s\,ds}{s^4} = 8\pi\Big(\frac{k}{4E}\Big)^2\Big[-\frac{1}{2s^2}\Big]_{s_0}^1 = 4\pi\Big(\frac{k}{4E}\Big)^2\frac{1 - s_0^2}{s_0^2} = \pi\Big(\frac{k}{2E}\Big)^2\cot^2\frac{\theta_0}{2},
$$
which is $\pi b(\theta_0)^2$ by (3.26). The integral reproduces the geometric statement "everything with $b < b(\theta_0)$ scatters beyond $\theta_0$."

(c) As $\theta_0 \to 0$, $\cot^2(\theta_0/2) \approx 4/\theta_0^2 \to \infty$. A pure $1/r$ potential deflects every particle however far away it passes, so the total cross-section is infinite. In matter, the atomic electrons cancel the nuclear charge beyond about $10^{-10}$ m, so particles with larger $b$ are not deflected at all, which cuts off the small angles.

(d) $\cot^2(0.5°)/\cot^2(45°) = \cot^2(0.5°) = 13{,}131$. Almost all Rutherford scattering is at tiny angles; the large-angle events Geiger and Marsden counted are a minute fraction.

## B. Core

**B1.**
(a) $L = \tfrac12m_1|\dot{\mathbf r}_1|^2 + \tfrac12m_2|\dot{\mathbf r}_2|^2 - V(|\mathbf r_1 - \mathbf r_2|) + m_1\mathbf g\cdot\mathbf r_1 + m_2\mathbf g\cdot\mathbf r_2$. The kinetic and interaction terms become $\tfrac12M|\dot{\mathbf R}|^2 + \tfrac12\mu|\dot{\mathbf r}|^2 - V(r)$ exactly as in Theorem 3.1. The external term is $\mathbf g\cdot(m_1\mathbf r_1 + m_2\mathbf r_2) = M\mathbf g\cdot\mathbf R$. So
$$
L = \big[\tfrac12M|\dot{\mathbf R}|^2 + M\mathbf g\cdot\mathbf R\big] + \big[\tfrac12\mu|\dot{\mathbf r}|^2 - V(r)\big].
$$
The $\mathbf R$ equation is $M\ddot{\mathbf R} = M\mathbf g$: $\mathbf R(t) = \mathbf R_0 + \mathbf V_0t + \tfrac12\mathbf gt^2$.

(b) $\mathbf g$ appears only in the $\mathbf R$ bracket, so $\mu\ddot{\mathbf r} = -\nabla V(r)$, unchanged. A uniform field accelerates both bodies equally and cannot change their separation.

(c) The external potential energy is $m_1\Phi(\mathbf r_1) + m_2\Phi(\mathbf r_2)$ (it enters $L$ with a minus sign) with $\mathbf r_1 = \mathbf R + (m_2/M)\mathbf r$, $\mathbf r_2 = \mathbf R - (m_1/M)\mathbf r$. Taylor-expand $\Phi(\mathbf R + \boldsymbol\delta) = \Phi + \boldsymbol\delta\cdot\nabla\Phi + \tfrac12\boldsymbol\delta^\top H\boldsymbol\delta$ with $H_{ij} = \partial_i\partial_j\Phi(\mathbf R)$:
- zeroth order: $(m_1 + m_2)\Phi(\mathbf R) = M\Phi(\mathbf R)$;
- first order: $\big[m_1\tfrac{m_2}{M} - m_2\tfrac{m_1}{M}\big]\mathbf r\cdot\nabla\Phi = 0$;
- second order: $\tfrac12\big[m_1\tfrac{m_2^2}{M^2} + m_2\tfrac{m_1^2}{M^2}\big]\mathbf r^\top H\mathbf r = \tfrac12\mu\,\mathbf r^\top H(\mathbf R)\,\mathbf r$.

The first coupling is the second-order term $-\tfrac12\mu\,\mathbf r^\top H(\mathbf R)\mathbf r$ in $L$: it depends on $\mathbf R$ through the field's curvature and on $\mathbf r$ quadratically. This is the **tidal** coupling. For the Earth–Moon pair in the Sun's field it is the dominant perturbation of the Moon's orbit. The linear term vanished because $\mathbf R$ is the mass-weighted center, the same cancellation as in step 3 of Theorem 3.1.

**B2.**
(a) $u = (1 + e\cos\varphi)/c$ gives $u'' + u = 1/c$. By (3.10), $-\dfrac{\mu}{\ell^2u^2}F = \dfrac1c$, so $F = -\dfrac{\ell^2u^2}{\mu c} = -\dfrac{\ell^2/\mu c}{r^2}$: inverse square, with $k = \ell^2/\mu c$, consistent with $c = \ell^2/\mu k$.

(b) $u = \sec\varphi/(2R)$. Then $u' = \sec\varphi\tan\varphi/(2R)$ and
$$
u'' = \frac{\sec\varphi\tan^2\varphi + \sec^3\varphi}{2R} = \frac{\sec\varphi(2\sec^2\varphi - 1)}{2R}, \qquad u'' + u = \frac{\sec^3\varphi}{R}.
$$
With $\sec\varphi = 2Ru$: $u'' + u = 8R^2u^3$ (SymPy). By (3.10), $F = -\dfrac{\ell^2u^2}{\mu}(u'' + u) = -\dfrac{8\ell^2R^2}{\mu}u^5 = -\dfrac{8\ell^2R^2}{\mu r^5}$. An attractive inverse-fifth-power force.

(c) For any orbit, $\dot r = -(\ell/\mu)u'$ and $r\dot\varphi = \ell/\mu r = (\ell/\mu)u$, so $v^2 = (\ell/\mu)^2(u'^2 + u^2)$. Here $u'^2 + u^2 = \sec^2\varphi(\tan^2\varphi + 1)/(4R^2) = \sec^4\varphi/(4R^2) = 4R^2u^4$ (SymPy), so $\tfrac12\mu v^2 = 2\ell^2R^2/\mu r^4$. The potential with $F = -dV/dr$ and $V(\infty) = 0$ is $V = -2\ell^2R^2/\mu r^4$ (check: $-dV/dr = -8\ell^2R^2/\mu r^5$ ✓). So $E = 2\ell^2R^2/\mu r^4 - 2\ell^2R^2/\mu r^4 = 0$ at every point, and $v \propto 1/r^2 \to \infty$ at the center. Forces steeper than $1/r^3$ allow orbits that fall into the center (compare Theorem 3.9).

**B3.** $r_p = 6371 + 600 = 6971$ km, $r_a = 6371 + 39{,}700 = 46{,}071$ km. By (3.13): $a = (r_p + r_a)/2 = 26{,}521$ km and $e = (r_a - r_p)/(r_a + r_p) = 39{,}100/53{,}042 = 0.7372$. By (3.16) with $k/\mu = GM_\oplus$: $T = 2\pi\sqrt{a^3/GM_\oplus} = 11.940$ h, within 0.2% of half a sidereal day (11.967 h), so the ground track nearly repeats twice a day. Vis-viva:
$$
v_p = \sqrt{GM_\oplus\Big(\frac{2}{r_p} - \frac1a\Big)} = 9.966\text{ km/s}, \qquad v_a = \sqrt{GM_\oplus\Big(\frac{2}{r_a} - \frac1a\Big)} = 1.508\text{ km/s}.
$$
Check: $r_pv_p = 6971\times9.966$ and $r_av_a = 46{,}071\times1.508$ are both $6.948\times10^4$ km²/s ✓ (asserted equal to $10^{-12}$).

**B4.**
(a) The attractive force magnitude is $f = dV/dr = ke^{-r/\lambda}\big(1/r^2 + 1/\lambda r\big)$. By (3.20) the circular orbit is stable iff $\omega_r^2/\omega_\varphi^2 = 3 + r f'/f > 0$, and $rf'/f = d\ln f/d\ln r$. With $x = r/\lambda$:
$$
\ln f = \ln k - \frac{r}{\lambda} + \ln\big(r^{-2} + \lambda^{-1}r^{-1}\big) \;\Rightarrow\; r\frac{d\ln f}{dr} = -x + r\,\frac{-2r^{-3} - \lambda^{-1}r^{-2}}{r^{-2} + \lambda^{-1}r^{-1}} = -x - \frac{2 + x}{1 + x}
$$
(multiply numerator and denominator of the fraction by $r^2$). Then
$$
3 + \frac{rf'}{f} = 3 - x - \frac{2 + x}{1 + x} = \frac{(3 - x)(1 + x) - (2 + x)}{1 + x} = \frac{1 + x - x^2}{1 + x}.
$$
The denominator is positive, so stability $\Leftrightarrow x^2 - x - 1 < 0 \Leftrightarrow x < (1 + \sqrt5)/2 = 1.618034$ (SymPy factors the condition as $-(r^2 - r - 1)/(r + 1)$ in units $\lambda = 1$).

(b) As $\lambda \to \infty$, $x \to 0$ at any fixed $r_0$, the condition always holds, and $\omega_r^2/\omega_\varphi^2 \to 1$: the Kepler result.

(c) Screened Coulomb potentials (the Debye-screened field of a charge in a plasma or electrolyte) and Yukawa's potential for the exchange of a massive particle (the nuclear force) have this form. Inside about 1.6 screening lengths the attraction falls slowly enough to support stable circles; beyond that the exponential makes it fall faster than the centrifugal barrier, and a slightly perturbed circular orbit either escapes or spirals inward.

**B5.**
(a) $\mu\ddot x = -kx$, $\mu\ddot y = -ky$: two independent oscillators of the same frequency $\omega = \sqrt{k/\mu}$. So $\mathbf r(t) = \mathbf P\cos\omega t + \mathbf Q\sin\omega t$ for constant vectors $\mathbf P, \mathbf Q$. This is the image of the unit circle $(\cos\omega t, \sin\omega t)$ under the linear map with columns $\mathbf P, \mathbf Q$ (LA-1's column picture), which is an ellipse centered on the origin (a circle or a segment in degenerate cases). Rotating to its principal axes and shifting $t$ puts it in the form $x = A\cos\omega t$, $y = B\sin\omega t$.

(b) $E = \tfrac12\mu\omega^2(A^2\sin^2\omega t + B^2\cos^2\omega t) + \tfrac12k(A^2\cos^2\omega t + B^2\sin^2\omega t) = \tfrac12k(A^2 + B^2)$, using $\mu\omega^2 = k$. $\ell = \mu(x\dot y - y\dot x) = \mu(A\cos\omega t\cdot B\omega\cos\omega t + B\sin\omega t\cdot A\omega\sin\omega t) = \mu\omega AB$.

(c) $r^2 = A^2\cos^2\omega t + B^2\sin^2\omega t = \tfrac12(A^2 + B^2) + \tfrac12(A^2 - B^2)\cos2\omega t$. So $r$ oscillates at angular frequency $2\omega$ while the particle goes around once in $2\pi/\omega$: two radial cycles per revolution, maximum on the major axis and minimum on the minor axis, at every amplitude. This is $\omega_r/\omega_\varphi = \sqrt{1 + 3} = 2$ from (3.21), exact here and not only for small oscillations, as Bertrand's theorem requires.

(d) $E - \omega\ell = \tfrac12k(A^2 + B^2) - \mu\omega^2AB = \tfrac12k(A - B)^2 \ge 0$, with equality iff $A = B$, a circle. For fixed $\ell$ the circle is the lowest-energy orbit, i.e. it sits at the minimum of $V_{\rm eff}$. The Kepler analogue is $E \ge -\mu k^2/2\ell^2$ from (3.9), with equality for $e = 0$.

**B6.**
(a) At closest approach $\dot r = 0$, so by (3.8) with $V = +k/r$ and $\ell^2 = 2\mu Eb^2$: $E = Eb^2/r^2 + k/r$, i.e. $r^2 - (k/E)r - b^2 = 0$. The positive root is
$$
r_{\min} = \frac{k}{2E}\Big[1 + \sqrt{1 + (2Eb/k)^2}\Big] = \frac{k}{2E}(1 + e) = \frac{k}{2E}\Big[1 + \frac{1}{\sin(\theta/2)}\Big],
$$
using $e = \sqrt{1 + (2Eb/k)^2}$ and $\sin(\theta/2) = 1/e$ from Theorem 3.13. (Equivalently $r_{\min} = c/(e - 1)$ from step 3 of that theorem: $c = 2Eb^2/k = (k/2E)(e^2 - 1)$, and $(e^2 - 1)/(e - 1) = e + 1$.)

(b) Gold, 5.0 MeV: $k/2E = 22.75$ fm. At $60°$, $\sin30° = 0.5$, $r_{\min} = 22.75\times3 = 68.3$ fm. At $180°$, $r_{\min} = 22.75\times2 = 45.5$ fm $= k/E$. Both are far outside the 7 fm nucleus.

(c) Aluminium: $k = 2\times13\times1.43996 = 37.44$ MeV fm; contact distance $1.2(3 + 1.587) = 5.50$ fm. Head-on collisions reach contact when $k/E \le 5.50$ fm, i.e. $E \ge 37.44/5.50 = 6.80$ MeV. At 5 MeV the head-on $r_{\min}$ is already only 7.49 fm, so on light targets the backward-angle cross-section of higher-energy alphas departs from (3.27). Deviations of exactly this kind on light nuclei gave the first estimates of nuclear radii.

**B7.** All speeds from vis-viva with $GM_\odot$; $r_1 = 1$ AU, $r_2 = 1.524$ AU, $a_t = 1.262$ AU.
(a) Earth's circular speed is 29.78 km/s and the transfer ellipse's perihelion speed is 32.73 km/s, so $\Delta v_1 = 2.946$ km/s. At aphelion the transfer speed is 21.48 km/s and Mars's circular speed 24.13 km/s, so $\Delta v_2 = 2.650$ km/s. Transfer time: half the period of the transfer ellipse, $\pi\sqrt{a_t^3/GM_\odot} = 258.9$ d.

(b) The spacecraft arrives $180°$ from where it left. Mars's period is $T_M = 687.2$ d, so during the transfer Mars moves $360°\times258.9/687.2 = 135.6°$. To meet the spacecraft, Mars must start $180° - 135.6° = 44.4°$ ahead of Earth.

(c) The Earth–Mars angle repeats with the synodic period $1/(1/T_\oplus - 1/T_M) = 780$ d: a launch window about every 26 months, which is why Mars missions cluster in those years.

**B8.**
(a) At an apsis $\dot r = 0$, so $\mathbf p = p\,\hat{\boldsymbol\varphi}$ and $\mathbf p\times\mathbf L = p\ell\,\hat{\mathbf r}$ (as in Theorem 3.15). Then $\mathbf A = (p\ell - \mu k)\hat{\mathbf r} = (\ell^2/r - \mu k)\hat{\mathbf r}$, parallel to $\hat{\mathbf r}$. At perihelion $\ddot r \ge 0$, so by (3.7) $\ell^2/\mu r_p^3 \ge k/r_p^2 + \lambda/r_p^3$, giving $\ell^2/r_p \ge \mu k + \mu\lambda/r_p > \mu k$: $\mathbf A$ points along $+\hat{\mathbf r}$, at the particle. The code (Part 10) finds the angle between $\mathbf A$ and $\mathbf r$ at every perihelion below $10^{-6}$.

(b) Because $\mathbf A$ points at each perihelion, its direction at successive perihelia advances by the perihelion advance itself. Measured: $0.076163$ rad per orbit, equal to $2\pi(1/\beta - 1) = 0.076163$ from (3.22).

(c) Differentiate: $\dot{\mathbf A} = \dot{\mathbf p}\times\mathbf L - \mu k\,d\hat{\mathbf r}/dt$ with $\dot{\mathbf p} = (-k/r^2 - \lambda/r^3)\hat{\mathbf r}$. The $-k/r^2$ part cancels the last term exactly (CM-4 Example 4.6). What is left is
$$
\dot{\mathbf A} = -\frac{\lambda}{r^3}\,\hat{\mathbf r}\times\mathbf L = \frac{\lambda\ell}{r^3}\,\hat{\boldsymbol\varphi},
$$
using $\hat{\mathbf r}\times\hat{\mathbf z} = -\hat{\boldsymbol\varphi}$. This is perpendicular to $\mathbf A$ at the apsides (pure rotation there) but not elsewhere, so $|\mathbf A|$ changes within an orbit and returns to the same value each perihelion. Predict the extremes with $\mu = k = 1$, $\ell^2 = 0.84$, $\lambda = 0.02$, $r_p = 0.6$. Perihelion: $|\mathbf A|/\mu k = \ell^2/\mu kr_p - 1 = 0.84/0.6 - 1 = 0.4000$. Aphelion: by Theorem 3.11, $u = u_c(1 + e'\cos\beta\varphi)$ with $u_c = \mu k/\ell^2\beta^2 = 1/(\ell^2 - \mu\lambda) = 1/0.82 = 1.21951$; $u_p = 1/0.6$ gives $e' = 0.36667$, so $u_a = u_c(1 - e') = 0.77236$, $r_a = 1.29474$. There $\ell^2u_a = 0.6488 < \mu k$, so $\mathbf A$ points along $-\hat{\mathbf r}_a$ (still toward the perihelion side) with $|\mathbf A|/\mu k = 1 - 0.84\times0.77236 = 0.3512$. The code finds $|\mathbf A|/\mu k$ ranging over $[0.3512, 0.4000]$.

## C. Challenge

**C1.**
(a) Substitute $u = (GM/h^2)w$: $(GM/h^2)(w'' + w) = GM/h^2 + (3GM/c^2)(GM/h^2)^2w^2$. Divide by $GM/h^2$: $w'' + w = 1 + \delta w^2$, $\delta = 3(GM)^2/(h^2c^2)$.

(b) At order $\delta^0$: $w_0'' + w_0 = 1$, solved by $w_0 = 1 + e\cos\varphi$ (perihelion at $\varphi = 0$). At order $\delta$:
$$
w_1'' + w_1 = w_0^2 = 1 + 2e\cos\varphi + e^2\cos^2\varphi = \Big(1 + \frac{e^2}{2}\Big) + 2e\cos\varphi + \frac{e^2}{2}\cos2\varphi.
$$
The term $2e\cos\varphi$ is resonant: it drives the oscillator $w'' + w$ at its own natural frequency.

(c) Particular solutions term by term: the constant gives $1 + e^2/2$; $\tfrac{e^2}{2}\cos2\varphi$ gives $-\tfrac{e^2}{6}\cos2\varphi$ (since $(-4 + 1)C = e^2/2$); the resonant term gives $e\varphi\sin\varphi$, because $(\varphi\sin\varphi)'' = 2\cos\varphi - \varphi\sin\varphi$, so $(\varphi\sin\varphi)'' + \varphi\sin\varphi = 2\cos\varphi$. The first two pieces are bounded and $2\pi$-periodic: they reshape the orbit slightly but do not move the perihelion from orbit to orbit. The third grows. Keeping it,
$$
w \approx 1 + e\big[\cos\varphi + \delta\varphi\sin\varphi\big] \approx 1 + e\cos[(1 - \delta)\varphi],
$$
since $\cos[(1 - \delta)\varphi] = \cos\varphi\cos\delta\varphi + \sin\varphi\sin\delta\varphi \approx \cos\varphi + \delta\varphi\sin\varphi$ while $\delta\varphi \ll 1$. Maxima of $w$ (perihelia) are at $(1 - \delta)\varphi = 2\pi m$, i.e. $\varphi \approx 2\pi m(1 + \delta)$: an advance of $2\pi\delta$ per orbit. Finally, the unperturbed orbit has $c = \ell^2/\mu k = h^2/GM = a(1 - e^2)$, so $h^2 = GMa(1 - e^2)$ and
$$
2\pi\delta = \frac{6\pi(GM)^2}{h^2c^2} = \frac{6\pi GM}{c^2a(1 - e^2)},
$$
which is (3.23). For Mercury the project gets $5.019\times10^{-7}$ rad per orbit, $42.98''$ per century.

(d) The project finds $0.006299$ rad per orbit against $2\pi\delta = 0.006283$, a ratio $1.00252$. The first-order theory drops $O(\delta^2)$ terms. Estimate them from the near-circular orbit: the circular solution satisfies $w_c = 1 + \delta w_c^2$, so $w_c = 1 + \delta + O(\delta^2)$; linearizing $w = w_c + x$ gives $x'' + (1 - 2\delta w_c)x = 0$, frequency $\omega^2 = 1 - 2\delta - 2\delta^2$. The advance is $2\pi(\omega^{-1} - 1)$ with $\omega^{-1} = (1 - y)^{-1/2} = 1 + y/2 + 3y^2/8$, $y = 2\delta + 2\delta^2$: $\omega^{-1} = 1 + \delta + 2.5\delta^2$. So the advance is $2\pi\delta(1 + 2.5\delta)$, ratio $1.00250$ for $\delta = 10^{-3}$, matching the measured $1.00252$ (the remaining $2\times10^{-5}$ is shared between the $O(\delta^3)$ terms and the $e = 0.2$ correction).

**C2.**
(a) $(x + ae)^2/a^2 + y^2/b^2 = \cos^2\psi + \sin^2\psi = 1$: the ellipse (3.13), center at $(-ae, 0)$, focus at the origin. And
$$
x^2 + y^2 = a^2(\cos\psi - e)^2 + a^2(1 - e^2)\sin^2\psi = a^2\big[1 - 2e\cos\psi + e^2\cos^2\psi\big] = a^2(1 - e\cos\psi)^2,
$$
so $r = a(1 - e\cos\psi)$. At $\psi = 0$, $r = a(1 - e) = r_p$.

(b) The area swept from the focus is $dA = \tfrac12(x\,dy - y\,dx)$. With $dx = -a\sin\psi\,d\psi$, $dy = b\cos\psi\,d\psi$:
$$
x\,dy - y\,dx = ab(\cos\psi - e)\cos\psi\,d\psi + ab\sin^2\psi\,d\psi = ab(1 - e\cos\psi)\,d\psi.
$$
So $A(\psi) = \tfrac12ab(\psi - e\sin\psi)$ from perihelion. By (3.6) the area grows uniformly and reaches $\pi ab$ in time $T$: $A = \pi ab\,t/T$. Equate: $2\pi t/T = \psi - e\sin\psi$. (Given $t$, solve for $\psi$ numerically; the project uses Newton's method.)

(c) From (b), $dt = (T/2\pi)(1 - e\cos\psi)\,d\psi$. Hence
$$
\Big\langle\frac1r\Big\rangle = \frac1T\int_0^T\frac{dt}{r} = \frac{1}{2\pi}\int_0^{2\pi}\frac{(1 - e\cos\psi)\,d\psi}{a(1 - e\cos\psi)} = \frac1a,
$$
$$
\langle r\rangle = \frac{1}{2\pi}\int_0^{2\pi}a(1 - e\cos\psi)^2d\psi = \frac{a}{2\pi}\int_0^{2\pi}\big(1 - 2e\cos\psi + e^2\cos^2\psi\big)d\psi = a\Big(1 + \frac{e^2}{2}\Big).
$$
The first is the virial theorem result of CM-4 Exercise B6(d); the second shows the planet lingers far out. For the Molniya orbit the project's time-sampled average is $\langle r\rangle/a = 1.271696 = 1 + e^2/2$.

(d) $r > a \Leftrightarrow \cos\psi < 0 \Leftrightarrow \psi \in (\pi/2, 3\pi/2)$. The elapsed time is $(T/2\pi)\big[(\psi - e\sin\psi)\big]_{\pi/2}^{3\pi/2} = (T/2\pi)(\pi + 2e)$, a fraction $\tfrac12 + e/\pi = 0.7346$ for $e = 0.7372$. With apogee placed over high northern latitudes, the satellite spends about three-quarters of each 12-hour orbit high and slow above that hemisphere, where geostationary satellites sit low on the horizon.

**C3.**
(a) $F = -dV/dr = 2\alpha/r^3 = 2\alpha u^3$ (repulsive). By (3.10), $u'' + u = -\dfrac{\mu}{\ell^2u^2}2\alpha u^3 = -\dfrac{2\mu\alpha}{\ell^2}u$, i.e. $u'' + \gamma^2u = 0$ with $\gamma^2 = 1 + 2\mu\alpha/\ell^2$. Choose $\varphi = 0$ for the incoming direction ($u = 0$ there): $u = A\sin\gamma\varphi$. The particle is at infinity again at the next zero, $\varphi = \pi/\gamma$, so the position vector sweeps $\pi/\gamma$ in total. As in steps 4–5 of Theorem 3.13, the velocity at infinity is radial: incoming along polar angle $\pi$ (from $\varphi = 0$ toward the center), outgoing along $\pi/\gamma$. The deflection is $\theta = \pi - \pi/\gamma$. (Check: $\alpha \to 0$ gives $\gamma \to 1$, $\theta \to 0$.)

(b) With $\ell^2 = 2\mu Eb^2$: $\gamma^2 = 1 + \alpha/Eb^2$, so $1/\gamma = b/\sqrt{b^2 + \alpha/E}$ and $\theta = \pi\big(1 - b/\sqrt{b^2 + \alpha/E}\big)$. Let $x = 1 - \theta/\pi = b/\sqrt{b^2 + \alpha/E}$; then $b^2 = (\alpha/E)\,x^2/(1 - x^2)$. Write (3.24) as $\dfrac{d\sigma}{d\Omega} = \dfrac{1}{2\sin\theta}\Big|\dfrac{d(b^2)}{d\theta}\Big|$. Now $\dfrac{d(b^2)}{dx} = \dfrac{\alpha}{E}\cdot\dfrac{2x(1 - x^2) + 2x^3}{(1 - x^2)^2} = \dfrac{\alpha}{E}\cdot\dfrac{2x}{(1 - x^2)^2}$ and $|dx/d\theta| = 1/\pi$, so
$$
\frac{d\sigma}{d\Omega} = \frac{\alpha}{\pi E}\cdot\frac{x}{(1 - x^2)^2\sin\theta}.
$$
Substitute $x = (\pi - \theta)/\pi$ and $1 - x^2 = \theta(2\pi - \theta)/\pi^2$:
$$
\frac{d\sigma}{d\Omega} = \frac{\alpha}{\pi E}\cdot\frac{(\pi - \theta)/\pi}{\theta^2(2\pi - \theta)^2\pi^{-4}\sin\theta} = \frac{\alpha\pi^2}{E}\cdot\frac{\pi - \theta}{\theta^2(2\pi - \theta)^2\sin\theta}.
$$
At $90°$ this is $8\alpha/(9\pi E)$, i.e. $0.5659$ for $\alpha/E = 2$. The project integrates the trajectory with $b = 1$, $\alpha = 1$, $E = \tfrac12$ and gets $76.08°$, the formula's value, and checks the cross-section against a numerical $db/d\theta$.

(c) As $\theta \to 0$: $\pi - \theta \to \pi$, $2\pi - \theta \to 2\pi$, $\sin\theta \to \theta$, so $d\sigma/d\Omega \to \alpha\pi/(4E\theta^3)$.

*Impulse approximation.* For small deflection the path is nearly the straight line $(x, b)$, $x = vt$. For $V = \alpha r^{-s}$ the transverse force there is $F_\perp = -\partial V/\partial y = \alpha s\,b\,(x^2 + b^2)^{-(s+2)/2}$. The transverse momentum gained is
$$
\Delta p_\perp = \frac1v\int_{-\infty}^\infty F_\perp\,dx = \frac{\alpha s}{v}\,b^{-s}\,I_s, \qquad I_s = \int_{-\infty}^\infty\frac{d\xi}{(1 + \xi^2)^{(s+2)/2}}
$$
(substituting $x = b\xi$). So $\theta \approx \Delta p_\perp/\mu v = \alpha sI_s/(2Eb^s) \propto b^{-s}$, i.e. $b \propto \theta^{-1/s}$, and by (3.24) with $\sin\theta \approx \theta$,
$$
\frac{d\sigma}{d\Omega} \approx \frac{b}{\theta}\Big|\frac{db}{d\theta}\Big| = \frac1s\frac{b^2}{\theta^2} \propto \theta^{-2 - 2/s}.
$$
Checks: $s = 2$ gives $\theta^{-3}$, and with $I_2 = \pi/2$, $\theta \approx \pi\alpha/(2Eb^2)$, so $b^2 = \pi\alpha/(2E\theta)$ and $d\sigma/d\Omega \approx \tfrac12b^2/\theta^2 = \pi\alpha/(4E\theta^3)$, matching the exact limit. $s = 1$ with $\alpha = k$, $I_1 = 2$: $\theta \approx k/(Eb)$ and $d\sigma/d\Omega \approx b^2/\theta^2 = k^2/(E^2\theta^4)$, which is Rutherford's $(k/4E)^2\sin^{-4}(\theta/2) \approx (k/4E)^2\cdot16/\theta^4$ at small angles. The small-angle cross-section measures the power of the potential's tail.

**C4.**
(a) $\nabla\Phi_{\rm eff} = Gm_1\dfrac{\mathbf r - \mathbf r_1}{\rho_1^3} + Gm_2\dfrac{\mathbf r - \mathbf r_2}{\rho_2^3} - \Omega^2\mathbf r$, with $\rho_i = |\mathbf r - \mathbf r_i|$. Because the center of mass is at the origin, $m_1\mathbf r_1 + m_2\mathbf r_2 = \mathbf 0$, so $M\mathbf r = m_1(\mathbf r - \mathbf r_1) + m_2(\mathbf r - \mathbf r_2)$. With $\Omega^2 = GM/d^3$:
$$
\nabla\Phi_{\rm eff} = Gm_1(\mathbf r - \mathbf r_1)\Big(\frac{1}{\rho_1^3} - \frac{1}{d^3}\Big) + Gm_2(\mathbf r - \mathbf r_2)\Big(\frac{1}{\rho_2^3} - \frac{1}{d^3}\Big).
$$
At an equilateral point $\rho_1 = \rho_2 = d$, both brackets vanish: $\nabla\Phi_{\rm eff} = \mathbf 0$, for any masses. (The project evaluates it numerically for Sun–Jupiter: $2.2\times10^{-16}$.)

(b) If $\mathbf r$ is not on the line through $\mathbf r_1$ and $\mathbf r_2$, the vectors $\mathbf r - \mathbf r_1$ and $\mathbf r - \mathbf r_2$ are not parallel, hence linearly independent, so $\nabla\Phi_{\rm eff} = \mathbf 0$ forces both coefficients to vanish: $\rho_1 = \rho_2 = d$. Two circles of radius $d$ about $\mathbf r_1$ and $\mathbf r_2$ (which are $d$ apart) meet in exactly two points, the apexes of the two equilateral triangles: L4 and L5. On the line the two vectors are parallel, only one scalar equation remains, and it has three solutions, the collinear points L1, L2, L3 (Euler, 1767).

(c) $\tfrac12(1 - \sqrt{23/27}) = 0.03852$. Sun–Jupiter: $m_J/(M_\odot + m_J) = 9.537\times10^{-4} < 0.0385$, stable, and indeed thousands of Trojan asteroids librate about Jupiter's L4 and L5. Earth–Moon: $0.01214 < 0.0385$, stable in the restricted three-body model (the Sun's tidal pull, B1(c), complicates the real Earth–Moon L4 and L5). Stability at a *maximum* of $\Phi_{\rm eff}$ is possible only because the Coriolis force turns a drift away into a circulation around the point; that is what Routh's analysis, and the project's "Your turn" 1, show.
