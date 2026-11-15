# EM-7: Electromagnetic Waves: Wave Equation, Plane Waves, Polarization, Reflection and Refraction, Waveguides, Dispersion

**Time:** ~12 hours · **Prerequisites:** EM-6, F3 · **Unlocks:** EM-8 (radiation), QM-1 (photons, blackbody), QM-2 (wave packets)

> **Conventions in this module.** SI units; $c = 1/\sqrt{\mu_0\varepsilon_0}$, $Z_0 = \mu_0c = 376.73\ \Omega$. **Physics sign:** a monochromatic field is $\tilde{\mathbf E} = \tilde{\mathbf E}_0e^{i(\mathbf k\cdot\mathbf r - \omega t)}$ (tildes mark complex quantities) and the physical field is its real part. Engineers write $e^{j(\omega t - \mathbf k\cdot\mathbf r)}$, in effect $j = -i$, which conjugates every complex amplitude, impedance and permittivity (F1 §3.8; §8). Wave vector $\mathbf k = k\hat{\mathbf k}$, $k = 2\pi/\lambda$. Linear media: $\mathbf D = \varepsilon\mathbf E$, $\mathbf H = \mathbf B/\mu$, $\varepsilon = \varepsilon_r\varepsilon_0$, $\mu = \mu_r\mu_0$, index $n = \sqrt{\varepsilon_r\mu_r}$, speed $v = c/n$, impedance $Z = \sqrt{\mu/\varepsilon}$. Intensity $I = \langle S\rangle$. Polarization vector $\hat{\mathbf n}\perp\mathbf k$; Jones vector $(\tilde E_x, \tilde E_y)$ for $\mathbf k\parallel\hat{\mathbf z}$. **Interfaces:** boundary $z = 0$, medium 1 ($n_1$) in $z<0$, light incident from medium 1, angles from the normal, plane of incidence $xz$; $r = \tilde E_{0R}/\tilde E_{0I}$, $t = \tilde E_{0T}/\tilde E_{0I}$, power coefficients $R$, $T$. **p sign:** Griffiths' convention, the three p unit vectors coincide at normal incidence ($r_p = r_s$ there); Hecht and Born–Wolf flip $r_p$ (§8). TE$_{mn}$: $m$ half-wavelengths across the wide side $a$, $n$ across $b\le a$.

## 0. Where this sits

EM-6 derived the vacuum wave equation (EM-6 (3.24)), the interface conditions (EM-6 (3.23)), Poynting's theorem and field momentum (EM-6 (3.26), (3.31)), and ran a 1D FDTD solver. You also need F3's Fourier transform and transfer function (F3 (3.12), (3.23)) and CM-1 §3.8's damped driven oscillator. This module solves Maxwell's equations for waves and asks what matter does to them. EM-8 asks where waves come from (radiation); QM-1 quantizes them into photons; QM-2's free-particle wave packet is the Gaussian spreading computation of §3.14, line for line.

## 1. The problem that forced this idea

Maxwell's equations predicted transverse waves at $1/\sqrt{\mu_0\varepsilon_0}$, numerically the speed of light (EM-6 §1), and the prediction went untested for twenty years. Heinrich Hertz, at Karlsruhe in 1887–88, drove a spark gap at a few hundred megahertz, detected the waves with a second gap across the room, reflected them off a zinc sheet into standing waves, measured their wavelength, refracted them through a prism of pitch, and showed they were polarized. Light was electromagnetism.

Optics had been doing without Maxwell, with an awkward success. Malus found in 1808 that light reflected from glass is partly polarized; Brewster (1815) found the angle, $\tan\theta_B = n$, at which it is completely polarized. Augustin Fresnel, treating light as transverse vibrations of an elastic ether, derived in 1821–23 how much light reflects from glass at each angle and polarization. His formulas were right. His ether, a solid stiff enough to carry waves at $3\times10^8$ m/s yet offering no resistance to the planets, had no consistent mechanics, and he needed assumptions about its density inside glass that nobody could justify. Maxwell's theory recovers Fresnel's equations from two boundary conditions and no ether (§3.10).

The second puzzle was colour. Newton's prism (1666) bent violet more than red. Cauchy (1836) fitted $n(\lambda)\approx A + B/\lambda^2$; Sellmeier (1871) did better with poles at certain wavelengths, as if the glass contained oscillators. Lorentz, in the 1890s, supplied them: electrons bound in molecules, driven by the wave's field, responding differently at different frequencies.

So the problem: given Maxwell's equations, what do their solutions look like in empty space, at an interface, in a conductor, in a pipe; and why does glass have an index that depends on colour?

## 2. Intuition

**A wave is a shape that moves.** $f(z - vt)$ is the graph of $f$ slid right by $vt$. You know the 1D wave equation from PDEs and the string from F2 §1 and CM-6 Example 4.5. Everything here is that equation applied to the components of $\mathbf E$ and $\mathbf B$, plus the constraints Maxwell's equations add.

**A plane EM wave is a rigid pattern of crossed fields.** EM-6 §3.8 found that $\mathbf E = f(x - ct)\hat{\mathbf y}$ comes with $\mathbf B = f(x - ct)\hat{\mathbf z}/c$: $\mathbf E$, $\mathbf B$ and the direction of travel form a right-handed triad, the fields peak together, $E/B = c$. Guess: the energy is shared equally between $\mathbf E$ and $\mathbf B$ and flows along $\mathbf E\times\mathbf B$ at $c$; since momentum is energy over $c$ for anything moving at $c$, light absorbed by a black surface pushes with pressure $I/c$.

**Reflection is what the boundary forces.** Tie a light string to a heavy one and send a pulse along the light one. At the knot both strings move together, and the massless knot feels no net force, so the slopes match. One incoming wave cannot meet both conditions; a reflected and a transmitted wave appear to balance the books. A heavy second string acts almost like a wall, so the reflection comes back inverted. Glass is a "heavy string" for light (light is slower in it), so guess the reflection coefficient by putting $v = c/n$ into the string answer.

```
  light string (fast)         knot          heavy string (slow)
  ---->  incident             |
  <----  reflected, inverted  |  ---->  transmitted, smaller
```

At an angle, the conditions must hold at every point of the interface at every moment, so the crests of the incident, reflected and transmitted waves must meet on the interface. The slower wave's crests are closer together, so to keep step along the surface the transmitted wave tilts toward the normal. That is Snell's law, from geometry alone.

**Dispersion is driven electrons.** In glass each electron is a mass on a spring (CM-1 §3.8) driven by the passing $\mathbf E$. Below resonance it moves in phase with the drive, its dipole adds to the medium's response, and the wave slows: $n > 1$. Drive closer to the (ultraviolet) resonance and the response grows, so $n$ grows: blue bends more than red. At resonance the electrons absorb. Free electrons have no spring at all, always lag the drive by $180°$, and give $n < 1$; below a certain frequency no wave gets in.

**A pipe is a filter.** Between two metal plates a wave travels as two plane waves zig-zagging between the walls, arranged so that their tangential $\mathbf E$ cancels at each wall: a whole number of half-wavelengths fits across. Below the frequency where one half-wavelength fits, nothing propagates. Above it, energy moves down the pipe at $c\cos\theta < c$ while the pattern where crests hit the wall slides along at $c/\cos\theta > c$. Guess: $v_pv_g = c^2$.

## 3. Formal development

### 3.1 The 1D wave equation and d'Alembert's solution

$$\frac{\partial^2f}{\partial z^2} = \frac{1}{v^2}\frac{\partial^2f}{\partial t^2}. \tag{3.1}$$

▶ **Try it first.** Change variables to $\xi = z - vt$ and $\eta = z + vt$ and rewrite (3.1).

**Theorem 3.1 (d'Alembert, 1747).** Every $C^2$ solution of (3.1) on the whole line has the form
$$f(z, t) = g(z - vt) + h(z + vt). \tag{3.2}$$

*Derivation.*
1. By the chain rule, $\partial_z = \partial_\xi + \partial_\eta$ and $\partial_t = -v\partial_\xi + v\partial_\eta$.
2. So $\partial_z^2 = \partial_\xi^2 + 2\partial_\xi\partial_\eta + \partial_\eta^2$ and $v^{-2}\partial_t^2 = \partial_\xi^2 - 2\partial_\xi\partial_\eta + \partial_\eta^2$.
3. Subtracting, (3.1) becomes $4\,\partial_\xi\partial_\eta f = 0$.
4. $\partial_\eta(\partial_\xi f) = 0$ says $\partial_\xi f$ depends on $\xi$ only; integrate in $\xi$: $f = g(\xi) + h(\eta)$. ∎

■ **What just happened.** The wave equation says nothing more than "two shapes, one moving right and one moving left, at speed $v$." The initial shape and velocity fix $g$ and $h$. Every other result in this module is about *which* shapes Maxwell's equations and the boundaries allow.

### 3.2 Sinusoids, complex notation, time averages

The sinusoid $f = A\cos(kz - \omega t + \delta)$ has wavelength $\lambda = 2\pi/k$, period $2\pi/\omega$ and $\omega = vk$. We write it as $f = \operatorname{Re}\tilde f$ with $\tilde f = \tilde Ae^{i(kz - \omega t)}$, $\tilde A = Ae^{i\delta}$.

Why this is legitimate: if $L$ is linear with real coefficients (derivatives, real constants), $L(\operatorname{Re}\tilde f) = \operatorname{Re}(L\tilde f)$, so if $L\tilde f = 0$ then $L(\operatorname{Re}\tilde f) = 0$. Maxwell's equations in linear media are of this form, and linearity (superposition) is also why Fourier analysis (F3) decomposes any wave into sinusoids. Products break the rule: $\operatorname{Re}(\tilde f\tilde g)\ne\operatorname{Re}\tilde f\,\operatorname{Re}\tilde g$. For energy and the Poynting vector take real parts first, or use:

**Lemma 3.2 (time averages).** For complex amplitudes $\tilde A$, $\tilde B$,
$$\big\langle\operatorname{Re}(\tilde Ae^{-i\omega t})\,\operatorname{Re}(\tilde Be^{-i\omega t})\big\rangle = \tfrac12\operatorname{Re}(\tilde A\tilde B^*). \tag{3.3}$$

*Derivation.* 1. $\operatorname{Re}(\tilde Ae^{-i\omega t}) = \tfrac12(\tilde Ae^{-i\omega t} + \tilde A^*e^{i\omega t})$. 2. The product is $\tfrac14(\tilde A\tilde Be^{-2i\omega t} + \tilde A^*\tilde B^*e^{2i\omega t} + \tilde A\tilde B^* + \tilde A^*\tilde B)$. 3. The $e^{\pm2i\omega t}$ terms average to zero over a period, leaving $\tfrac14(\tilde A\tilde B^* + \text{c.c.}) = \tfrac12\operatorname{Re}(\tilde A\tilde B^*)$. ∎ (Project Part 1 checks this with SymPy.)

### 3.3 Two strings: the template for every interface

Strings of mass per length $\mu_1$ ($z<0$) and $\mu_2$ ($z>0$), joined by a massless knot at $z = 0$, under tension $T$: speeds $v_j = \sqrt{T/\mu_j}$. A wave $\tilde A_Ie^{i(k_1z - \omega t)}$ arrives from the left.

▶ **Try it first.** Write incident, reflected and transmitted waves with one common $\omega$. What two conditions hold at the knot? Solve for $\tilde A_R$, $\tilde A_T$.

*Derivation.*
1. Left: $\tilde A_Ie^{i(k_1z - \omega t)} + \tilde A_Re^{i(-k_1z - \omega t)}$. Right: $\tilde A_Te^{i(k_2z - \omega t)}$. One $\omega$ on both sides because the knot oscillates at one frequency; then $k_j = \omega/v_j$.
2. The string is continuous: $f(0^-, t) = f(0^+, t)$ gives $\tilde A_I + \tilde A_R = \tilde A_T$.
3. The knot is massless, so the transverse force on it, $T(\partial_zf|_{0^+} - \partial_zf|_{0^-})$, must vanish: $k_1(\tilde A_I - \tilde A_R) = k_2\tilde A_T$.
4. Solve the $2\times2$ system:
$$\tilde A_R = \frac{k_1 - k_2}{k_1 + k_2}\tilde A_I = \frac{v_2 - v_1}{v_2 + v_1}\tilde A_I,\qquad\tilde A_T = \frac{2k_1}{k_1 + k_2}\tilde A_I = \frac{2v_2}{v_1 + v_2}\tilde A_I. \tag{3.4}$$ ∎

■ **What just happened.** A slower second string ($v_2 < v_1$) reflects inverted; a fixed wall is $v_2\to0$ ($\tilde A_R = -\tilde A_I$). The recipe (one $\omega$ everywhere, as many conditions as unknown amplitudes, a small linear system) is the whole of §3.8–3.12 and §3.15.

### 3.4 Plane waves in vacuum

In vacuum every Cartesian component of $\mathbf E$ and $\mathbf B$ obeys the wave equation with speed $c$ (EM-6 (3.24)). Try $\tilde{\mathbf E} = \tilde{\mathbf E}_0e^{i(\mathbf k\cdot\mathbf r - \omega t)}$, $\tilde{\mathbf B} = \tilde{\mathbf B}_0e^{i(\mathbf k\cdot\mathbf r - \omega t)}$; on such fields $\nabla\to i\mathbf k$, $\partial_t\to-i\omega$. The wave equation follows from Maxwell's equations, but not every solution of it solves them.

▶ **Try it first.** Apply $\nabla\to i\mathbf k$, $\partial_t\to-i\omega$ to all four Maxwell equations in vacuum.

**Theorem 3.3 (Transverse, crossed, in phase).** A monochromatic plane wave solves Maxwell's equations in vacuum if and only if $\omega = ck$ and
$$\mathbf k\cdot\tilde{\mathbf E}_0 = 0,\qquad\mathbf k\cdot\tilde{\mathbf B}_0 = 0, \tag{3.5}$$
$$\tilde{\mathbf B}_0 = \frac1c\,\hat{\mathbf k}\times\tilde{\mathbf E}_0. \tag{3.6}$$

*Derivation.*
1. $\nabla\cdot\mathbf E = 0$ becomes $i\mathbf k\cdot\tilde{\mathbf E}_0 = 0$; $\nabla\cdot\mathbf B = 0$ becomes $i\mathbf k\cdot\tilde{\mathbf B}_0 = 0$. That is (3.5): both fields are transverse.
2. Faraday, $\nabla\times\mathbf E = -\partial_t\mathbf B$, becomes $i\mathbf k\times\tilde{\mathbf E}_0 = i\omega\tilde{\mathbf B}_0$, so $\tilde{\mathbf B}_0 = (k/\omega)\hat{\mathbf k}\times\tilde{\mathbf E}_0$.
3. Maxwell–Ampère, $\nabla\times\mathbf B = c^{-2}\partial_t\mathbf E$, becomes $i\mathbf k\times\tilde{\mathbf B}_0 = -i(\omega/c^2)\tilde{\mathbf E}_0$. Insert step 2 and use $\hat{\mathbf k}\times(\hat{\mathbf k}\times\tilde{\mathbf E}_0) = \hat{\mathbf k}(\hat{\mathbf k}\cdot\tilde{\mathbf E}_0) - \tilde{\mathbf E}_0 = -\tilde{\mathbf E}_0$ (by step 1): $-(k^2/\omega)\tilde{\mathbf E}_0 = -(\omega/c^2)\tilde{\mathbf E}_0$, i.e. $\omega = ck$. Then step 2 is (3.6). ∎

■ **What just happened.** $\mathbf E\perp\mathbf B\perp\mathbf k$, a right-handed triad; $\tilde{\mathbf B}_0$ is a real multiple of $\hat{\mathbf k}\times\tilde{\mathbf E}_0$, so the fields are in phase; $|\mathbf B| = |\mathbf E|/c$. The real field is $\mathbf E = E_0\hat{\mathbf n}\cos(\mathbf k\cdot\mathbf r - \omega t + \delta)$, $\mathbf B = (E_0/c)(\hat{\mathbf k}\times\hat{\mathbf n})\cos(\mathbf k\cdot\mathbf r - \omega t + \delta)$, with polarization vector $\hat{\mathbf n}\perp\hat{\mathbf k}$. Every solution of the vacuum equations is a superposition of these (F3: Fourier transform in $\mathbf r$).

### 3.5 Energy, intensity, momentum, radiation pressure

▶ **Try it first.** Insert the real fields above into $u = \tfrac12(\varepsilon_0E^2 + B^2/\mu_0)$ and $\mathbf S = \mathbf E\times\mathbf B/\mu_0$ (EM-6 (3.26)), then average over a period.

*Derivation.*
1. $B^2/\mu_0 = E^2/(c^2\mu_0) = \varepsilon_0E^2$: the magnetic energy equals the electric energy. So $u = \varepsilon_0E_0^2\cos^2(\cdot)$.
2. $\mathbf E\times\mathbf B = (E_0^2/c)\cos^2(\cdot)\,\hat{\mathbf n}\times(\hat{\mathbf k}\times\hat{\mathbf n}) = (E_0^2/c)\cos^2(\cdot)\,\hat{\mathbf k}$, so $\mathbf S = cu\,\hat{\mathbf k}$: energy moves at $c$ along $\hat{\mathbf k}$.
3. $\langle\cos^2\rangle = \tfrac12$ (Lemma 3.2 with $\tilde A = \tilde B = 1$):
$$\langle u\rangle = \tfrac12\varepsilon_0E_0^2,\qquad I = \langle S\rangle = \tfrac12c\varepsilon_0E_0^2. \tag{3.7}$$
4. The momentum density is $\mathbf g = \mathbf S/c^2$ (EM-6 (3.31)), so $\langle\mathbf g\rangle = (I/c^2)\hat{\mathbf k}$. In time $\Delta t$, the momentum arriving on area $A$ at normal incidence is that in a cylinder of length $c\Delta t$: $\Delta p = (I/c^2)Ac\,\Delta t$.
5. A black surface absorbs it; a mirror reverses it and so receives twice as much. Pressure $= \Delta p/(A\Delta t)$:
$$P_{\rm abs} = \frac Ic,\qquad P_{\rm mirror} = \frac{2I}c. \tag{3.8}$$ ∎

■ **What just happened.** $I/c = \langle u\rangle$: the pressure on an absorber equals the energy density, as for any beam of particles moving at $c$. Sunlight at Earth (the solar constant, 1361 W/m²) pushes with $4.54\ \mu$Pa on a black surface (Example 4.1). QM-1 reads $p = E/c$ photon by photon.

### 3.6 Polarization and the Jones calculus

Take $\mathbf k = k\hat{\mathbf z}$. By (3.5), $\tilde{\mathbf E}_0$ has only $x$ and $y$ components. The **Jones vector** is
$$\mathbf J = \begin{pmatrix}\tilde E_x\\\tilde E_y\end{pmatrix}\in\mathbb C^2,\qquad\mathbf E(z, t) = \operatorname{Re}\big[\mathbf Je^{i(kz - \omega t)}\big]. \tag{3.9}$$
Linear: components in phase, e.g. $(\cos\theta, \sin\theta)$. Circular: equal magnitudes, phases $\pm\pi/2$ apart; for $\mathbf J = \tfrac{E_0}{\sqrt2}(1, i)$, at $z = 0$, $\mathbf E = \tfrac{E_0}{\sqrt2}(\cos\omega t, \sin\omega t)$ turns from $\hat{\mathbf x}$ toward $\hat{\mathbf y}$, counterclockwise for an observer facing the oncoming beam. Anything else is elliptical.

Optical elements act linearly on $\mathbf J$: they are $2\times2$ complex matrices, composed right to left (LA-1).
$$P(\theta) = \mathbf u\mathbf u^\top,\ \mathbf u = \begin{pmatrix}\cos\theta\\\sin\theta\end{pmatrix};\qquad W(\Gamma) = \begin{pmatrix}1 & 0\\0 & e^{i\Gamma}\end{pmatrix};\qquad M(\theta) = R(\theta)\,M\,R(\theta)^\top. \tag{3.10}$$
$P(\theta)$, a polarizer with transmission axis at $\theta$ from $x$, is the orthogonal projection onto $\mathbf u$ (LA-5). $W(\Gamma)$ is a wave plate with fast axis along $x$: a birefringent slab of thickness $d$ whose slow component picks up the extra phase $\Gamma = (n_s - n_f)\omega d/c$, which with $e^{i(kz - \omega t)}$ appears as $e^{+i\Gamma}$. $\Gamma = \pi/2$ is a quarter-wave plate (QWP), $\Gamma = \pi$ a half-wave plate (HWP). $M(\theta)$ is any element $M$ rotated so its axis lies at $\theta$; $R$ is the rotation matrix.

**Malus's law.** Horizontal light through $P(\theta)$: $P(\theta)(1, 0)^\top = \cos\theta\,\mathbf u$, intensity $\cos^2\theta$. Unpolarized light through a polarizer keeps half its intensity.

**Stokes parameters and the Poincaré sphere.** With $S_0 = |\tilde E_x|^2 + |\tilde E_y|^2$,
$$S_1 = |\tilde E_x|^2 - |\tilde E_y|^2,\qquad S_2 = 2\operatorname{Re}(\tilde E_x^*\tilde E_y),\qquad S_3 = 2\operatorname{Im}(\tilde E_x^*\tilde E_y). \tag{3.11}$$
For a fully polarized beam $S_1^2 + S_2^2 + S_3^2 = S_0^2$, so $(S_1, S_2, S_3)/S_0$ lies on the unit sphere: horizontal at $(1, 0, 0)$, $45°$ linear at $(0, 1, 0)$, $\tfrac1{\sqrt2}(1, i)$ at the north pole. The overall phase of $\mathbf J$ drops out, so the sphere is exactly the set of polarization states.

**This is a qubit (QC-1).** A Jones vector normalized to $|\tilde E_x|^2 + |\tilde E_y|^2 = 1$, with overall phase ignored, is a qubit state $\alpha|H\rangle + \beta|V\rangle$, and the Poincaré sphere *is* the Bloch sphere, with axes relabelled (Bloch $z\leftrightarrow S_1$, $x\leftrightarrow S_2$, $y\leftrightarrow S_3$). Wave plates are unitary matrices, so they rotate the sphere: $W(\pi/2) = \operatorname{diag}(1, i)$ is the phase gate $S$, and a HWP with fast axis at $22.5°$ is exactly the Hadamard gate $\tfrac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}$ (project Part 2 asserts both). A polarizer is a projective measurement.

### 3.7 Waves in linear media

With no free charge or current, Maxwell's equations in a linear, homogeneous medium are the vacuum ones with $\mu_0\varepsilon_0\to\mu\varepsilon$ (EM-6 (3.22)). Everything in §3.4–3.5 carries over with $c\to v$:
$$v = \frac1{\sqrt{\mu\varepsilon}} = \frac cn,\quad n = \sqrt{\varepsilon_r\mu_r};\qquad B_0 = \frac{E_0}v;\qquad Z\equiv\frac{E_0}{H_0} = \sqrt{\frac\mu\varepsilon};\qquad I = \tfrac12\varepsilon vE_0^2 = \frac{E_0^2}{2Z}. \tag{3.12}$$
For non-magnetic media ($\mu = \mu_0$, almost every optical material) $n = \sqrt{\varepsilon_r}$, $Z = Z_0/n$, and $I = \tfrac12nc\varepsilon_0E_0^2$: the same $E_0$ carries $n$ times more power in glass. That factor of $n$ is why transmitted power is not $|t|^2$.

### 3.8 Reflection and transmission at normal incidence

▶ **Try it first.** Copy §3.3: incident $\tilde E_{0I}\hat{\mathbf x}$ travelling $+z$ in medium 1, reflected $\tilde E_{0R}\hat{\mathbf x}$ travelling $-z$, transmitted $\tilde E_{0T}\hat{\mathbf x}$. Get each $\mathbf B$ from (3.6) with $c\to v$, then impose the interface conditions.

*Derivation.*
1. By (3.6), $\tilde{\mathbf B}_I = (\tilde E_{0I}/v_1)\hat{\mathbf y}$, $\tilde{\mathbf B}_T = (\tilde E_{0T}/v_2)\hat{\mathbf y}$, and for the reflected wave $\hat{\mathbf k} = -\hat{\mathbf z}$, so $\tilde{\mathbf B}_R = -(\tilde E_{0R}/v_1)\hat{\mathbf y}$.
2. EM-6 (3.23) with no free charge or current: $D^\perp$, $B^\perp$, $\mathbf E^\parallel$, $\mathbf H^\parallel$ continuous. There are no $z$ components, so the first two say $0 = 0$.
3. $\mathbf E^\parallel$: $\tilde E_{0I} + \tilde E_{0R} = \tilde E_{0T}$.
4. $\mathbf H^\parallel$: $\frac1{\mu_1v_1}(\tilde E_{0I} - \tilde E_{0R}) = \frac1{\mu_2v_2}\tilde E_{0T}$, i.e. $\tilde E_{0I} - \tilde E_{0R} = \beta\tilde E_{0T}$ with $\beta = \frac{\mu_1v_1}{\mu_2v_2} = \frac{\mu_1n_2}{\mu_2n_1}$.
5. Add and subtract: $\tilde E_{0T} = \frac2{1+\beta}\tilde E_{0I}$, $\tilde E_{0R} = \frac{1-\beta}{1+\beta}\tilde E_{0I}$. For $\mu_1 = \mu_2$:
$$r = \frac{n_1 - n_2}{n_1 + n_2},\qquad t = \frac{2n_1}{n_1 + n_2}. \tag{3.13}$$ ∎

■ **What just happened.** This is (3.4) with $v_j = c/n_j$, as §2 guessed. Going into a denser medium ($n_2 > n_1$) gives $r < 0$: the reflected $\mathbf E$ is flipped, a phase shift of $\pi$. From dense to rare there is no flip. The project's FDTD sends a pulse into $n = 1.5$ and measures $r = -0.1999$, $t = 0.7998$.

**Theorem 3.4 (Energy conservation).** With $R = |r|^2$ and $T = (n_2/n_1)|t|^2$,
$$R + T = 1. \tag{3.14}$$
*Proof.* 1. $R = I_R/I_I$ and $T = I_T/I_I$. By (3.12) with $\mu = \mu_0$, $I\propto nE_0^2$, so $T = (n_2/n_1)|t|^2$. 2. $R + T = \frac{(n_1 - n_2)^2 + 4n_1n_2}{(n_1 + n_2)^2} = \frac{(n_1 + n_2)^2}{(n_1 + n_2)^2} = 1$. ∎

### 3.9 Oblique incidence: the laws of reflection and refraction

Incident, reflected and transmitted waves $\tilde{\mathbf E}_{0j}e^{i(\mathbf k_j\cdot\mathbf r - \omega_jt)}$, $j = I, R, T$, with $|\mathbf k_I| = |\mathbf k_R| = n_1\omega/c$, $|\mathbf k_T| = n_2\omega/c$.

▶ **Try it first.** Every boundary condition at $z = 0$ has the form $(\cdots)e^{i(\mathbf k_I\cdot\mathbf r - \omega_It)} + (\cdots)e^{i(\mathbf k_R\cdot\mathbf r - \omega_Rt)} = (\cdots)e^{i(\mathbf k_T\cdot\mathbf r - \omega_Tt)}$ for all $x$, $y$, $t$. What does that force on the exponents?

*Derivation.*
1. **Lemma.** If $\sum_ja_je^{iq_jx} = 0$ for all $x$ with distinct real $q_j$, then every $a_j = 0$. (Differentiate $0, 1, \dots, N-1$ times and set $x = 0$: the coefficient matrix $(iq_j)^m$ is Vandermonde, invertible for distinct $q_j$.)
2. Applied in $t$: all three frequencies are equal. Applied in $x$ and $y$ on $z = 0$: the tangential components of the three wave vectors are equal, $(\mathbf k_I)_x = (\mathbf k_R)_x = (\mathbf k_T)_x$ and likewise for $y$. (If not, the coefficients would all vanish, which allows only zero fields.)
3. Choose axes with $\mathbf k_I$ in the $xz$ plane. Then $\mathbf k_R$ and $\mathbf k_T$ have no $y$ component: **the three wave vectors lie in one plane, the plane of incidence.**
4. $x$ components: $k_I\sin\theta_I = k_R\sin\theta_R$ with $k_I = k_R$, so **$\theta_R = \theta_I$**; and $k_I\sin\theta_I = k_T\sin\theta_T$:
$$n_1\sin\theta_I = n_2\sin\theta_T\qquad\text{(Snell's law)}. \tag{3.15}$$ ∎

■ **What just happened.** Phase matching alone gives the laws of geometrical optics; the amplitudes will need the boundary conditions. F6 Example 4.3 got (3.15) from Fermat's least time, and the two derivations are one: Fermat says the phase $\int\mathbf k\cdot d\mathbf l$ is stationary along the ray, and stationarity at the interface is conservation of tangential $\mathbf k$ (the same stationary phase picks the classical path in QM-9's path integral).

### 3.10 The Fresnel equations

From now on $\mu_1 = \mu_2 = \mu_0$ unless stated, and the common phase factor $e^{i(k_xx - \omega t)}$ is dropped at $z = 0$. Define
$$\alpha = \frac{\cos\theta_T}{\cos\theta_I},\qquad\beta = \frac{\mu_1n_2}{\mu_2n_1}\ \Big(= \frac{n_2}{n_1}\Big).$$

**p polarization (TM):** $\mathbf E$ in the plane of incidence. With $\hat{\mathbf k}_I = (\sin\theta_I, 0, \cos\theta_I)$, $\hat{\mathbf k}_R = (\sin\theta_I, 0, -\cos\theta_I)$, $\hat{\mathbf k}_T = (\sin\theta_T, 0, \cos\theta_T)$, choose the field directions
$$\hat{\mathbf e}_I = (\cos\theta_I, 0, -\sin\theta_I),\quad\hat{\mathbf e}_R = (\cos\theta_I, 0, \sin\theta_I),\quad\hat{\mathbf e}_T = (\cos\theta_T, 0, -\sin\theta_T),$$
all equal to $\hat{\mathbf x}$ at normal incidence (Griffiths' Fig. 9.15).

▶ **Try it first.** Compute the three $\mathbf B$'s from (3.6), write all four conditions of EM-6 (3.23), and count how many are independent.

*Derivation.*
1. $\hat{\mathbf k}_I\times\hat{\mathbf e}_I = \hat{\mathbf y}$, $\hat{\mathbf k}_R\times\hat{\mathbf e}_R = -\hat{\mathbf y}$, $\hat{\mathbf k}_T\times\hat{\mathbf e}_T = \hat{\mathbf y}$ (expand the cross products). So $\tilde{\mathbf B}_I = \frac{\tilde E_{0I}}{v_1}\hat{\mathbf y}$, $\tilde{\mathbf B}_R = -\frac{\tilde E_{0R}}{v_1}\hat{\mathbf y}$, $\tilde{\mathbf B}_T = \frac{\tilde E_{0T}}{v_2}\hat{\mathbf y}$.
2. **$B^\perp$:** no $\mathbf B$ has a $z$ component. $0 = 0$.
3. **$\mathbf E^\parallel$** ($x$ component): $\cos\theta_I(\tilde E_{0I} + \tilde E_{0R}) = \cos\theta_T\tilde E_{0T}$, so $\tilde E_{0I} + \tilde E_{0R} = \alpha\tilde E_{0T}$.
4. **$\mathbf H^\parallel$** ($y$ component): $\frac1{\mu_1v_1}(\tilde E_{0I} - \tilde E_{0R}) = \frac1{\mu_2v_2}\tilde E_{0T}$, so $\tilde E_{0I} - \tilde E_{0R} = \beta\tilde E_{0T}$.
5. **$D^\perp$:** $\varepsilon_1\sin\theta_I(\tilde E_{0I} - \tilde E_{0R}) = \varepsilon_2\sin\theta_T\tilde E_{0T}$. By Snell, $\sin\theta_T/\sin\theta_I = n_1/n_2$, and $\varepsilon_j = n_j^2/(\mu_jc^2)$, so $\varepsilon_2\sin\theta_T/(\varepsilon_1\sin\theta_I) = \mu_1n_2/(\mu_2n_1) = \beta$. **This is step 4 again.**
6. Two independent equations (steps 3, 4), two unknowns. Add and subtract: $\tilde E_{0T} = \frac2{\alpha + \beta}\tilde E_{0I}$, $\tilde E_{0R} = \frac{\alpha - \beta}{\alpha + \beta}\tilde E_{0I}$. Multiply through by $n_1\cos\theta_I$:
$$r_p = \frac{\alpha - \beta}{\alpha + \beta} = \frac{n_1\cos\theta_T - n_2\cos\theta_I}{n_1\cos\theta_T + n_2\cos\theta_I},\qquad t_p = \frac2{\alpha + \beta} = \frac{2n_1\cos\theta_I}{n_1\cos\theta_T + n_2\cos\theta_I}. \tag{3.16}$$ ∎

■ **What just happened.** The normal conditions were redundant, and that is general: $(\nabla\times\mathbf E)_z = \partial_xE_y - \partial_yE_x$ involves only tangential derivatives of tangential components, so continuity of $\mathbf E^\parallel$ forces continuity of $\partial_tB_z$, hence of $B_z$; $(\nabla\times\mathbf H)_z = \partial_tD_z$ does the same for $D^\perp$. **For time-dependent fields impose the two tangential conditions; the normal ones follow.** At $\theta_I = 0$, (3.16) reduces to (3.13).

**s polarization (TE):** $\mathbf E = \tilde E_0\hat{\mathbf y}$ for all three waves. Outline: (3.6) gives $\mathbf B$'s in the plane of incidence with $x$ components $-\cos\theta_I\tilde E_{0I}/v_1$, $+\cos\theta_I\tilde E_{0R}/v_1$, $-\cos\theta_T\tilde E_{0T}/v_2$. $\mathbf E^\parallel$ gives $\tilde E_{0I} + \tilde E_{0R} = \tilde E_{0T}$; $H_x$ gives $\tilde E_{0I} - \tilde E_{0R} = \alpha\beta\tilde E_{0T}$. Hence
$$r_s = \frac{1 - \alpha\beta}{1 + \alpha\beta} = \frac{n_1\cos\theta_I - n_2\cos\theta_T}{n_1\cos\theta_I + n_2\cos\theta_T},\qquad t_s = \frac2{1 + \alpha\beta} = \frac{2n_1\cos\theta_I}{n_1\cos\theta_I + n_2\cos\theta_T}. \tag{3.17}$$
Exercise B2 fills in every step.

**Power.** A beam of intensity $I$ at angle $\theta$ delivers $I\cos\theta$ per unit area of interface. So $R = |r|^2$ and, with (3.12), $T = \frac{I_T\cos\theta_T}{I_I\cos\theta_I} = \frac{n_2\cos\theta_T}{n_1\cos\theta_I}|t|^2 = \alpha\beta|t|^2$. Then $R_p + T_p = \frac{(\alpha - \beta)^2 + 4\alpha\beta}{(\alpha + \beta)^2} = 1$ and $R_s + T_s = \frac{(1 - \alpha\beta)^2 + 4\alpha\beta}{(1 + \alpha\beta)^2} = 1$. The project checks $R + T = 1$ at 500 angles both ways and confirms all four interface conditions on the assembled fields to $10^{-12}$.

### 3.11 Brewster's angle; total internal reflection and the evanescent wave

▶ **Try it first.** Set $r_p = 0$ in (3.16) and eliminate $\theta_T$ with Snell's law.

**Theorem 3.5 (Brewster).** For $\mu_1 = \mu_2$ and $n_1\ne n_2$, $r_p = 0$ exactly when
$$\tan\theta_B = \frac{n_2}{n_1}. \tag{3.18}$$
*Derivation.*
1. $r_p = 0\iff n_1\cos\theta_T = n_2\cos\theta_I$. With Snell, $\cos\theta_T = (n_2/n_1)\cos\theta_I$ and $\sin\theta_T = (n_1/n_2)\sin\theta_I$.
2. Square and add: $1 = (n_2/n_1)^2\cos^2\theta_I + (n_1/n_2)^2\sin^2\theta_I$.
3. Put $u = \tan^2\theta_I$, so $\cos^2\theta_I = 1/(1+u)$, $\sin^2\theta_I = u/(1+u)$: $1 + u = (n_2/n_1)^2 + (n_1/n_2)^2u$, i.e. $u\big(1 - (n_1/n_2)^2\big) = (n_2/n_1)^2\big(1 - (n_1/n_2)^2\big)$.
4. Divide by $1 - (n_1/n_2)^2\ne0$: $\tan^2\theta_I = (n_2/n_1)^2$. ∎

■ **What just happened.** At $\theta_B$, $\sin\theta_T = (n_1/n_2)\sin\theta_B = \cos\theta_B$, so the reflected and transmitted rays are perpendicular. The reflected wave is radiated by the electrons of medium 2, oscillating along $\hat{\mathbf e}_T$, which now points along the reflected ray; a dipole radiates nothing along its own axis (EM-8). So unpolarized light reflects at $\theta_B$ as pure s light, $\mathbf E$ parallel to the surface. Glare off water or a road (Brewster $53°$ for water) has horizontal $\mathbf E$, and polarized sunglasses, with vertical transmission axes, block it. Air to glass, $\theta_B = 56.31°$.

**Total internal reflection.** If $n_1 > n_2$, Snell gives $\sin\theta_T = (n_1/n_2)\sin\theta_I > 1$ beyond the **critical angle** $\theta_c = \arcsin(n_2/n_1)$ ($41.81°$ glass to air).

▶ **Try it first.** Keep the formulas and let $\cos\theta_T = \sqrt{1 - \sin^2\theta_T}$ be imaginary. Which root makes physical sense with $e^{-i\omega t}$?

*Derivation.*
1. Tangential matching (§3.9) still holds: $k_{Tx} = k_1\sin\theta_I$. Then $k_{Tz}^2 = k_2^2 - k_{Tx}^2 = (\omega/c)^2(n_2^2 - n_1^2\sin^2\theta_I) < 0$.
2. So $k_{Tz} = \pm i\kappa$, $\kappa = (\omega/c)\sqrt{n_1^2\sin^2\theta_I - n_2^2}$, and the transmitted field is $\propto e^{\mp\kappa z}e^{i(k_{Tx}x - \omega t)}$. Only $k_{Tz} = +i\kappa$ (field decays into medium 2) is allowed; that is $\cos\theta_T = +i\sqrt{\sin^2\theta_T - 1}$.
3. The field falls by $1/e$ in
$$d = \frac1\kappa = \frac{\lambda_0}{2\pi\sqrt{n_1^2\sin^2\theta_I - n_2^2}}, \tag{3.19}$$
with $\lambda_0$ the vacuum wavelength.
4. In (3.17), $r_s = (a - ib)/(a + ib)$ with $a = n_1\cos\theta_I$, $b = n_2|\cos\theta_T|$ real: a complex number over its conjugate, so $|r_s| = 1$. The same holds for $r_p$. ∎

■ **What just happened.** Beyond $\theta_c$ all the energy reflects, with phase shifts that differ between s and p (Exercise C2 turns the difference into Fresnel's glass quarter-wave plate). The field in medium 2 is an **evanescent wave**: it runs along the surface, decays away from it, and carries no average energy across it (Exercise B4). Glass to air at $45°$: $d = 0.450\lambda_0$, 285 nm at 633 nm. With the engineering sign the physical root is $-j\kappa$ (§8).

**Fibres and frustrated TIR.** An optical fibre is a glass core of index slightly above its cladding; light that meets the core wall beyond the critical angle bounces down the core by TIR with no loss per bounce, so a long link is limited only by absorption and scattering in the glass (about 0.2 dB/km near 1.55 µm). Bring a second glass block within a fraction of $d$ and the evanescent field reaches it and propagates again: **frustrated TIR**, the optical twin of tunnelling (QM-3).

### 3.12 Conductors: complex wave number and skin depth

In an ohmic conductor $\mathbf J_f = \sigma\mathbf E$ (EM-6 (3.1)); free charge decays as $e^{-\sigma t/\varepsilon}$ (continuity plus Gauss), so set $\rho_f = 0$. Then $\nabla\cdot\mathbf E = 0$, $\nabla\cdot\mathbf B = 0$, $\nabla\times\mathbf E = -\partial_t\mathbf B$, $\nabla\times\mathbf B = \mu\sigma\mathbf E + \mu\varepsilon\partial_t\mathbf E$.

▶ **Try it first.** Take the curl of Faraday as in EM-6 §3.8. Insert $e^{i(\tilde kz - \omega t)}$.

*Derivation.*
1. $\nabla\times(\nabla\times\mathbf E) = -\nabla^2\mathbf E = -\partial_t(\mu\sigma\mathbf E + \mu\varepsilon\partial_t\mathbf E)$, so $\nabla^2\mathbf E = \mu\sigma\partial_t\mathbf E + \mu\varepsilon\partial_t^2\mathbf E$. The first-derivative term is a diffusion term: Joule loss.
2. With $\partial_z\to i\tilde k$, $\partial_t\to-i\omega$:
$$\tilde k^2 = \mu\varepsilon\omega^2 + i\mu\sigma\omega. \tag{3.20}$$
3. Write $\tilde k = k + i\kappa$: the field is $e^{-\kappa z}e^{i(kz - \omega t)}$. Real and imaginary parts of (3.20): $k^2 - \kappa^2 = \mu\varepsilon\omega^2$, $2k\kappa = \mu\sigma\omega$. Squaring and adding, $(k^2 + \kappa^2)^2 = (\mu\varepsilon\omega^2)^2 + (\mu\sigma\omega)^2$, which with the first equation gives $\kappa = \omega\sqrt{\tfrac{\varepsilon\mu}2}\big[\sqrt{1 + (\sigma/\varepsilon\omega)^2} - 1\big]^{1/2}$ (and the same with $+1$ for $k$).
4. **Good conductor**, $\sigma\gg\varepsilon\omega$: drop the first term of (3.20). $\tilde k^2 = i\mu\sigma\omega$ and $\sqrt i = e^{i\pi/4} = (1 + i)/\sqrt2$, so $k = \kappa = \sqrt{\mu\sigma\omega/2}$ and the **skin depth** is
$$\delta = \frac1\kappa = \sqrt{\frac2{\mu\sigma\omega}}. \tag{3.21}$$ ∎

**Phase lag.** Faraday with complex $\tilde k$ gives $\tilde B_0 = (\tilde k/\omega)\tilde E_0$, and $\tilde k = Ke^{i\phi}$ with $K = \sqrt{k^2 + \kappa^2}$, $\phi = \arctan(\kappa/k)$:
$$B_0 = \frac K\omega E_0,\qquad\delta_B = \delta_E + \phi. \tag{3.22}$$
With $e^{-i\omega t}$, $\cos(kz - \omega t + \delta_E + \phi)$ peaks a time $\phi/\omega$ after $\cos(kz - \omega t + \delta_E)$: $\mathbf B$ lags $\mathbf E$ by $\phi$, which is $45°$ in a good conductor. And $B_0/E_0 = \sqrt{\mu\sigma/\omega}\gg1/c$: inside a metal the field is almost all magnetic.

**Why metals shine.** At a vacuum–metal face, §3.8 goes through with $v_2 = \omega/\tilde k$: $\tilde\beta = c\tilde k/\omega$, $r = (1 - \tilde\beta)/(1 + \tilde\beta)$. For copper at 1 GHz $|\tilde\beta|\approx3\times10^4$, so $r\approx-1$ and $R = 0.99991$ (Example 4.5). At optical frequencies the static $\sigma$ fails (electrons no longer collide within a cycle) and the free-electron plasma of §3.13 takes over, with the same conclusion below the plasma frequency.

### 3.13 The Lorentz oscillator model of dispersion

**Model.** A bound electron (charge $q$, mass $m$) is a damped oscillator with natural frequency $\omega_0$ and damping rate $\gamma$, driven by the wave's field (in a dilute medium the local field is the applied field; the magnetic force is smaller by $v/c$):
$$m(\ddot x + \gamma\dot x + \omega_0^2x) = q\tilde E_0e^{-i\omega t}. \tag{3.23}$$
This is CM-1 §3.8's oscillator, and its steady response is F3 (3.23)'s transfer function with the time sign flipped ($+ic\omega$ there becomes $-i\gamma\omega$ here, the $e^{\pm i\omega t}$ swap).

▶ **Try it first.** Solve (3.23) for the steady state. Turn the dipole moment $qx$ into a polarization $P$ and then into $\varepsilon$.

*Derivation.*
1. Try $\tilde x = \tilde x_0e^{-i\omega t}$: $(-\omega^2 - i\gamma\omega + \omega_0^2)\tilde x_0 = (q/m)\tilde E_0$, so $\tilde x_0 = \frac{q/m}{\omega_0^2 - \omega^2 - i\gamma\omega}\tilde E_0$.
2. Dipole moment $\tilde p = q\tilde x$. With $N$ molecules per volume, each having $f_j$ electrons of frequency $\omega_j$ and damping $\gamma_j$ ($f_j$ is the **oscillator strength**), $\tilde{\mathbf P} = N\sum_jf_jq\tilde x_j = \varepsilon_0\tilde\chi_e\tilde{\mathbf E}$ (EM-4 (3.5), (3.18)).
3. $\tilde\varepsilon_r = 1 + \tilde\chi_e$:
$$\tilde\varepsilon_r(\omega) = 1 + \frac{Nq^2}{m\varepsilon_0}\sum_j\frac{f_j}{\omega_j^2 - \omega^2 - i\gamma_j\omega}. \tag{3.24}$$
4. The wave equation in the medium is (3.20) with $\sigma = 0$ and $\tilde\varepsilon$: $\tilde k = \sqrt{\tilde\varepsilon_r}\,\omega/c$. Write $\tilde k = k + i\kappa$; intensity falls as $e^{-\alpha z}$ with **absorption coefficient** $\alpha = 2\kappa$.
5. For a gas, $\tilde\chi_e\ll1$ and $\sqrt{1 + \tilde\chi_e}\approx1 + \tilde\chi_e/2$. Separating real and imaginary parts:
$$n = 1 + \frac{Nq^2}{2m\varepsilon_0}\sum_j\frac{f_j(\omega_j^2 - \omega^2)}{(\omega_j^2 - \omega^2)^2 + \gamma_j^2\omega^2},\qquad\alpha = \frac{Nq^2\omega^2}{m\varepsilon_0c}\sum_j\frac{f_j\gamma_j}{(\omega_j^2 - \omega^2)^2 + \gamma_j^2\omega^2}. \tag{3.25}$$ ∎

■ **What just happened.** Away from resonances $n$ rises with $\omega$: **normal dispersion**. Within about $\gamma_j$ of a resonance $n$ falls with $\omega$ (**anomalous dispersion**), exactly where $\alpha$ peaks: anomalous dispersion and absorption come together (`lorentz_model.png`).

**Sellmeier and Cauchy.** With $\gamma_j = 0$ and $\lambda = 2\pi c/\omega$, (3.24) becomes $n^2 = 1 + \sum_jB_j\lambda^2/(\lambda^2 - \lambda_j^2)$ with $B_j = Nq^2f_j/(m\varepsilon_0\omega_j^2)$: Sellmeier's formula is the lossless Lorentz model, and glass catalogues still use it. Far below a single UV resonance ($\lambda\gg\lambda_1$), $\lambda^2/(\lambda^2 - \lambda_1^2)\approx1 + \lambda_1^2/\lambda^2$, so $n^2\approx n_0^2 + B_1\lambda_1^2/\lambda^2$ with $n_0^2 = 1 + B_1$, and to first order
$$n(\lambda)\approx A + \frac B{\lambda^2},\qquad A = n_0,\quad B = \frac{B_1\lambda_1^2}{2n_0}. \tag{3.26}$$
That is Cauchy's formula. The project's toy glass ($B_1 = 1.1$ at $\lambda_1 = 0.09\ \mu$m, plus an infrared resonance) has $n = 1.4571$ at 589 nm; a Cauchy fit over 0.4–0.7 µm is good to $1.4\times10^{-4}$, with $B = 0.00350\ \mu\text{m}^2$ against (3.26)'s $0.00307$, the gap coming from the next term of the expansion and the infrared resonance's tail.

**Why the sky is blue.** Far below resonance, step 1 gives $\tilde x_0\approx q\tilde E_0/(m\omega_0^2)$: an induced dipole of amplitude $p_0 = q^2E_0/(m\omega_0^2)$, independent of $\omega$. An oscillating dipole radiates a field proportional to its acceleration $\ddot p = -\omega^2p$, so the radiated power goes as $|\ddot p|^2$: $\langle P\rangle = \mu_0p_0^2\omega^4/(12\pi c)$ (Larmor's formula; EM-8 derives it). With $p_0$ fixed, scattered power $\propto\omega^4$ (Rayleigh, 1871): 450 nm blue is scattered $(700/450)^4 = 5.86$ times more than 700 nm red. The sky is scattered light, so it is blue; at sunset the direct beam crosses so much air that its blue is scattered out, so the sun is red. ("Blue reflects off the ocean" and "blue has more energy" are wrong.) Violet scatters more still, but sunlight has less of it and your eye is less sensitive to it.

**Plasmas.** For free electrons ($\omega_j = 0$, $f = 1$) and negligible damping, (3.24) is
$$\varepsilon_r = 1 - \frac{\omega_p^2}{\omega^2},\qquad\omega_p^2 = \frac{Ne^2}{m\varepsilon_0},\qquad\text{so}\qquad n^2 = 1 - \frac{\omega_p^2}{\omega^2},\quad\omega^2 = \omega_p^2 + c^2k^2. \tag{3.27}$$
Below the **plasma frequency** $\omega_p$, $k$ is imaginary: the wave is evanescent and, by step 4 of §3.11, totally reflected. The ionosphere has $N\sim10^{12}$ m⁻³, $f_p = 8.98$ MHz: shortwave radio bounces off it and travels beyond the horizon, while FM and GPS pass through. Copper's conduction electrons ($N = 8.47\times10^{28}$ m⁻³) give $\omega_p = 1.64\times10^{16}$ rad/s, $\lambda_p = 115$ nm: metals reflect visible light and turn transparent in the far ultraviolet. Above $\omega_p$, $n < 1$ and $v_p > c$, harmlessly (§3.14).

### 3.14 Group velocity and the spreading of wave packets

**Two waves.** $\cos(k_1z - \omega_1t) + \cos(k_2z - \omega_2t) = 2\cos\big(\tfrac{\Delta k}2z - \tfrac{\Delta\omega}2t\big)\cos(\bar kz - \bar\omega t)$: a carrier moving at $\bar\omega/\bar k$ inside an envelope moving at $\Delta\omega/\Delta k$. Shrink the spacing and the envelope speed becomes a derivative.

A pulse given at $z = 0$ as a time signal $E(0, t) = \int\hat A(\omega)e^{-i\omega t}d\omega$ propagates through a medium with wave number $k(\omega)$ as
$$E(z, t) = \int\hat A(\omega)\,e^{i(k(\omega)z - \omega t)}\,d\omega. \tag{3.28}$$
(Each frequency travels independently because the equations are linear and time-invariant: F3's "pure tones are eigenfunctions." The project does exactly this with an FFT.)

▶ **Try it first.** Let $\hat A$ be peaked at $\omega_c$, expand $k(\omega) = k_0 + k'q + \tfrac12k''q^2$ with $q = \omega - \omega_c$, and keep only the linear term.

*Derivation (first order).*
1. $E(z, t) = e^{i(k_0z - \omega_ct)}\int\hat A(\omega_c + q)e^{-iq(t - k'z)}dq = e^{i(k_0z - \omega_ct)}F(t - k'z)$, where $F(s) = \int\hat A(\omega_c + q)e^{-iqs}dq$ is the envelope.
2. The envelope is rigid and arrives at depth $z$ at time $k'z$. The **group velocity** is
$$v_g = \frac1{dk/d\omega} = \frac{d\omega}{dk},\qquad\text{and with }k = n(\omega)\omega/c:\quad v_g = \frac c{n + \omega\,dn/d\omega}\equiv\frac c{n_g}. \tag{3.29}$$ ∎

The carrier crests move at the **phase velocity** $v_p = \omega/k = c/n$. In normal dispersion $dn/d\omega > 0$, so $v_g < v_p$. The project's Lorentz medium at the carrier frequency has $n = 1.5774$, $v_p = 0.634c$, $v_g = 0.569c$, and the FDTD pulse peak travels at $0.570c$.

**Second order: Gaussian spreading.** ▶ **Try it first.** With $\hat A(\omega) = e^{-\tau_0^2q^2}$, the $q$ integral in step 1 picks up a factor $e^{ik''zq^2/2}$. Do the Gaussian integral with a complex coefficient.

*Derivation.*
1. At $z = 0$: $\int e^{-\tau_0^2q^2 - iqs}dq = \sqrt{\pi}/\tau_0\,e^{-s^2/(4\tau_0^2)}$ (the Gaussian transform, F3 (3.12) in other variables), so the intensity $|E|^2\propto e^{-s^2/(2\tau_0^2)}$ has standard deviation $\tau_0$.
2. At depth $z$, with $s = t - k'z$: $\int\exp\big(-aq^2 - iqs\big)dq$ with $a = \tau_0^2 - ik''z/2$.
3. The Gaussian integral holds for complex $a$ with $\operatorname{Re}a > 0$ (complete the square and shift the contour): the result is $\sqrt{\pi/a}\,e^{-s^2/(4a)}$.
4. $|e^{-s^2/(4a)}|^2 = \exp\big(-\tfrac{s^2}2\operatorname{Re}\tfrac1a\big)$ and $\operatorname{Re}\tfrac1a = \tau_0^2/|a|^2$. So the intensity is Gaussian in $s$ with standard deviation $|a|/\tau_0$:
$$\tau(z) = \tau_0\sqrt{1 + \Big(\frac{k''z}{2\tau_0^2}\Big)^2}. \tag{3.30}$$ ∎

■ **What just happened.** The pulse stays Gaussian, its centre moves at $v_g$, and it widens at a rate set by the **group-velocity dispersion** $k''$; short pulses (broad spectra) spread fastest. $\operatorname{Im}(1/a)$ adds a phase quadratic in $s$, so the instantaneous frequency sweeps across the pulse: it becomes **chirped**, red in front when $k'' > 0$. Run it backwards: a pulse launched with the opposite chirp (spectral phase $e^{ib_0q^2/2}$, $b_0 < 0$) is compressed until $b_0 + k''z = 0$ and spreads afterwards, $\tau = \tau_0\sqrt{1 + ((b_0 + k''z)/2\tau_0^2)^2}$. The project launches such a chirp into a Lorentz medium and measures widths 75.5, 60.4, 113.4 (FDTD) against 75.2, 60.0, 109.0 from this formula at depths 5, 1000, 3000 (grid units).

**QM-2 preview.** A free particle has $\psi(x, t) = \int\phi(k)e^{i(kx - \omega(k)t)}dk$ with $\omega = \hbar k^2/2m$ exactly. Swap the roles of $(z, t)$ and $(t, x)$: the same four steps, with $k''z$ replaced by $\omega''t = \hbar t/m$, give $\sigma(t) = \sigma_0\sqrt{1 + (\hbar t/2m\sigma_0^2)^2}$, exact because $\omega(k)$ is exactly quadratic (Exercise C3).

**Faster than light?** In an absorption band (3.25) makes $dn/d\omega$ large and negative, so $n_g = n + \omega\,dn/d\omega$ can drop below 1 or below 0: $v_g > c$ or negative (Exercise C1: for a dilute gas, $n_g(\omega_0) = 1 - Nq^2f/(m\varepsilon_0\gamma^2)$). It has been measured, and it does not violate relativity. First, (3.29) assumed $k(\omega)$ linear across the pulse spectrum; in an absorption band $\alpha$ varies strongly across it, the pulse is reshaped, and the output peak is built from the front of the input pulse, which was already there. Second, information rides on a front (a switch-on), which contains arbitrarily high frequencies, where every material has $\varepsilon\to\varepsilon_0$ ((3.24) as $\omega\to\infty$): fronts move at $c$ (Sommerfeld and Brillouin, 1914). The deep reason is causality: the polarization cannot respond before the field arrives, which constrains $\tilde\varepsilon(\omega)$ (the Kramers–Kronig relations, Jackson §7.10).

### 3.15 Guided waves

**Parallel plates (the picture).** Perfectly conducting plates at $y = 0$ and $y = b$. Subtract two plane waves with $\mathbf k = (0, \pm k_y, k_z)$ and $\mathbf E\parallel\hat{\mathbf x}$: the sum is $\propto\sin(k_yy)e^{i(k_zz - \omega t)}\hat{\mathbf x}$. $\mathbf E^\parallel = 0$ on both plates forces $k_y = n\pi/b$, and $k_y^2 + k_z^2 = \omega^2/c^2$ gives $\omega^2 = c^2k_z^2 + (n\pi c/b)^2$: §2's zig-zag made precise. Plates also carry a uniform mode with $\mathbf E\parallel\hat{\mathbf y}$ at $c$ for every frequency.

**Rectangular waveguide.** A hollow, perfectly conducting pipe, $0\le x\le a$, $0\le y\le b$, $b\le a$, vacuum inside. Look for $\tilde{\mathbf E} = \tilde{\mathbf E}_0(x, y)e^{i(kz - \omega t)}$, $\tilde{\mathbf B} = \tilde{\mathbf B}_0(x, y)e^{i(kz - \omega t)}$. At the walls $\mathbf E^\parallel = 0$ and $B^\perp = 0$ (fields vanish inside a perfect conductor; surface charge and current take up the jumps in EM-6 (3.23)).

▶ **Try it first.** Write the six components of the two curl equations with $\partial_z\to ik$, $\partial_t\to-i\omega$. Assume $E_z = 0$ (TE) and express $E_x$, $E_y$, $B_x$, $B_y$ through $B_z$.

*Derivation (TE modes).*
1. The six components:
$\partial_yE_z - ikE_y = i\omega B_x$, $\ ikE_x - \partial_xE_z = i\omega B_y$, $\ \partial_xE_y - \partial_yE_x = i\omega B_z$;
$\partial_yB_z - ikB_y = -\tfrac{i\omega}{c^2}E_x$, $\ ikB_x - \partial_xB_z = -\tfrac{i\omega}{c^2}E_y$, $\ \partial_xB_y - \partial_yB_x = -\tfrac{i\omega}{c^2}E_z$.
2. With $E_z = 0$, the second gives $B_y = (k/\omega)E_x$. Put it in the fourth: $\partial_yB_z = i(k^2/\omega - \omega/c^2)E_x$. With the first, $B_x = -(k/\omega)E_y$; in the fifth: $\partial_xB_z = (i/\omega)(\omega^2/c^2 - k^2)E_y$. So, writing $\gamma^2\equiv\omega^2/c^2 - k^2$,
$$E_x = \frac{i\omega}{\gamma^2}\partial_yB_z,\quad E_y = -\frac{i\omega}{\gamma^2}\partial_xB_z,\quad B_x = \frac{ik}{\gamma^2}\partial_xB_z,\quad B_y = \frac{ik}{\gamma^2}\partial_yB_z. \tag{3.31}$$
3. $B_z$ obeys the wave equation (EM-6 (3.24)), i.e. $(\partial_x^2 + \partial_y^2 + \gamma^2)B_z = 0$.
4. Boundary conditions: $E_y = 0$ at $x = 0, a$ means $\partial_xB_z = 0$ there; $E_x = 0$ at $y = 0, b$ means $\partial_yB_z = 0$ there. ($B^\perp = 0$ then holds automatically by (3.31).)
5. Separate, $B_z = X(x)Y(y)$: $X''/X + Y''/Y = -\gamma^2$, so $X'' = -k_x^2X$, $Y'' = -k_y^2Y$, $k_x^2 + k_y^2 = \gamma^2$. $X'(0) = 0$ picks $X = \cos k_xx$; $X'(a) = 0$ needs $\sin k_xa = 0$, $k_x = m\pi/a$. Likewise $k_y = n\pi/b$.
6. Result, for integers $m, n\ge0$ not both zero:
$$B_z = B_0\cos\frac{m\pi x}a\cos\frac{n\pi y}b,\qquad\omega^2 = c^2k^2 + \omega_{mn}^2,\qquad\omega_{mn} = c\pi\sqrt{\frac{m^2}{a^2} + \frac{n^2}{b^2}}. \tag{3.32}$$ ∎

(Project Part 8 checks with SymPy that (3.31)–(3.32) satisfy all of Maxwell's equations and the wall conditions.) TM modes ($B_z = 0$) follow the same steps with $E_z\propto\sin(m\pi x/a)\sin(n\pi y/b)$, which needs $m, n\ge1$.

■ **What just happened.** Below the **cutoff frequency** $\omega_{mn}$, $k$ is imaginary and the mode dies exponentially: the guide is a high-pass filter, mode by mode. The lowest is TE$_{10}$, $\omega_{10} = c\pi/a$: half a wavelength fits across the wide side. From (3.32),
$$v_p = \frac\omega k = \frac c{\sqrt{1 - (\omega_{mn}/\omega)^2}} > c,\qquad v_g = \frac{d\omega}{dk} = \frac{c^2k}\omega = c\sqrt{1 - (\omega_{mn}/\omega)^2} < c,\qquad v_pv_g = c^2. \tag{3.33}$$
The plasma (3.27) has the same dispersion relation with $\omega_p$ in place of $\omega_{mn}$. Energy travels at $v_g$; the faster phase pattern carries nothing. Example 4.6 computes both for an X-band guide.

**Theorem 3.6 (No TEM mode in a hollow guide).** A hollow guide whose cross-section is bounded by a single closed conducting wall cannot carry a wave with $E_z = B_z = 0$.

*Proof.*
1. With $E_z = B_z = 0$, the third component equation gives $\partial_xE_y - \partial_yE_x = 0$, and Gauss gives $\partial_xE_x + \partial_yE_y = 0$ (no $\partial_zE_z$).
2. The cross-section is simply connected, so the curl-free $\mathbf E_0$ is a gradient (EM-1 Theorem 3.16): $\mathbf E_0 = -\nabla\varphi$. Then step 1's divergence says $\nabla^2\varphi = 0$ in the cross-section.
3. $\mathbf E^\parallel = 0$ on the wall: the derivative of $\varphi$ along the wall vanishes, so $\varphi$ is constant on the (connected) wall.
4. A harmonic function equal to a constant on the whole boundary is that constant inside (EM-3 Theorem 3.3, uniqueness, in 2D). So $\mathbf E_0 = 0$, and then the second and first component equations give $B_x = B_y = 0$. ∎

**The coaxial line has one.** Two conductors, radii $a < b$, held at different potentials: step 3 now allows a different constant on each, and by symmetry the harmonic $\varphi$ depends on $s$ alone: $\varphi\propto\ln s$. So $\mathbf E_0 = (A/s)\hat{\mathbf s}$. With $E_z = B_z = 0$, the component equations give $B_y = (k/\omega)E_x$ (second) and $B_y = (\omega/kc^2)E_x$ (fourth); both hold only if $k = \omega/c$. Hence
$$\mathbf E = \frac As\cos(kz - \omega t)\,\hat{\mathbf s},\qquad\mathbf B = \frac A{cs}\cos(kz - \omega t)\,\hat{\boldsymbol\phi},\qquad k = \frac\omega c, \tag{3.34}$$
a TEM wave at $c$ ($c/n$ with a dielectric filling) with no cutoff, which is why coax carries everything from DC upward. Exercise B7 gets its impedance, $(Z_0/2\pi n)\ln(b/a)$: $50\ \Omega$ for polyethylene with $b/a = 3.5$.

**Optical fibre.** A fibre is a dielectric waveguide: the field oscillates in the core and is evanescent in the cladding, and the modes obey matching conditions like (3.16)–(3.17) at the core boundary instead of $\mathbf E^\parallel = 0$. A thin enough core carries only the lowest mode, which has no cutoff.

**Where waves come from.** Everything here is source-free. Waves are launched by accelerating charges: the Rayleigh dipole, the current in a half-wave antenna. EM-8 solves Maxwell's equations with sources (the inhomogeneous wave equations EM-6 (3.37), by retarded potentials) and derives the $\omega^4$ dipole power used above.

## 4. Worked examples

### Example 4.1: A 1 kW/m² plane wave; the pressure of sunlight (project Part 1)

From (3.7), $E_0 = \sqrt{2I/(c\varepsilon_0)} = \sqrt{2000/(2.998\times10^8\times8.854\times10^{-12})} = \sqrt{7.535\times10^5} = 868.0$ V/m. Then $B_0 = E_0/c = 2.895\times10^{-6}$ T, small next to Earth's 25–65 µT. Energy density $\langle u\rangle = I/c = 3.336\times10^{-6}$ J/m³. By (3.8), on a black surface $P = 3.336\ \mu$Pa, on a mirror $6.671\ \mu$Pa.

Sunlight above the atmosphere, $I = 1361$ W/m²: $E_0 = 1012.6$ V/m, $P_{\rm abs} = 4.540\ \mu$Pa, $P_{\rm mirror} = 9.080\ \mu$Pa. A $100\times100$ m mirror sail facing the Sun feels $9.080\times10^{-6}\times10^4 = 90.8$ mN, continuously and with no fuel: the case for solar sails (JAXA's IKAROS flew one in 2010).

### Example 4.2: A quarter-wave plate turns linear into circular (project Part 2)

Light polarized at $45°$, $\mathbf J = \tfrac1{\sqrt2}(1, 1)$, meets a QWP with fast axis along $x$, $W(\pi/2) = \operatorname{diag}(1, i)$. Output $\tfrac1{\sqrt2}(1, i)$. Stokes (3.11): $S_1 = \tfrac12 - \tfrac12 = 0$, $S_2 = 2\operatorname{Re}(\tfrac12\cdot i) = 0$, $S_3 = 2\operatorname{Im}(\tfrac i2) = 1$: the north pole. Physical field at $z = 0$: $E_x = \tfrac{E_0}{\sqrt2}\cos\omega t$, $E_y = \tfrac{E_0}{\sqrt2}\operatorname{Re}(ie^{-i\omega t}) = \tfrac{E_0}{\sqrt2}\sin\omega t$, a circle traversed from $x$ toward $y$.

Test: through a polarizer at any $\theta$, $|\mathbf u\cdot\mathbf J|^2 = \tfrac12|\cos\theta + i\sin\theta|^2 = \tfrac12$; the $45°$ input would give $\cos^2(\theta - 45°)$.

Turn it around: horizontal light, QWP with fast axis at $45°$. $R(45°)W(\pi/2)R(45°)^\top = \tfrac12\begin{pmatrix}1+i & 1-i\\1-i & 1+i\end{pmatrix}$, so $(1, 0)\mapsto\tfrac12(1+i, 1-i) = e^{i\pi/4}\tfrac1{\sqrt2}(1, -i)$: circular, the other way round, the south pole. Rotating the plate through $180°$ drags the output along a figure-eight, $(\cos^22\theta, \sin2\theta\cos2\theta, -\sin2\theta)$, on the Poincaré sphere (`poincare.png`). In QC-1's language, the $S$ gate takes $|+\rangle$ to $|{+i}\rangle$.

### Example 4.3: Normal incidence on glass (project Parts 3–4)

Air to glass, $n_1 = 1$, $n_2 = 1.5$. By (3.13), $r = (1 - 1.5)/2.5 = -0.2$ and $t = 2/2.5 = 0.8$. $R = 0.04$ and $T = 1.5\times0.64 = 0.96$; $R + T = 1$. Note $|t|^2 = 0.64$ is not the transmitted power. The reflected wave is inverted.

A window has two faces: $0.96^2 = 0.9216$, or, counting the bounces between faces incoherently, $T^2(1 + R^2 + R^4 + \cdots) = T^2/(1 - R^2) = 0.9231$. A lens with 12 uncoated surfaces would transmit $0.96^{12} = 0.613$ and fill the image with ghosts, which is why lenses are coated (Exercise B8). The project's FDTD measures $r = -0.1999$, $t = 0.7998$, $R = 0.0400$, $T = 0.9600$.

### Example 4.4: Fresnel coefficients at 45° (project Part 3)

Air to glass, $\theta_I = 45°$: $\sin\theta_T = 0.70711/1.5 = 0.47140$, $\theta_T = 28.126°$, $\cos\theta_T = 0.88192$. The four products: $n_1\cos\theta_I = 0.70711$, $n_2\cos\theta_T = 1.32288$, $n_2\cos\theta_I = 1.06066$, $n_1\cos\theta_T = 0.88192$.

s: $r_s = (0.70711 - 1.32288)/2.02999 = -0.30334$, $t_s = 1.41421/2.02999 = 0.69666$.
p: $r_p = (0.88192 - 1.06066)/1.94258 = -0.09201$, $t_p = 1.41421/1.94258 = 0.72801$.

$\alpha\beta = (0.88192/0.70711)(1.5) = 1.87083$. $R_s = 0.09201$, $T_s = 1.87083\times0.69666^2 = 0.90799$; $R_p = 0.008466$, $T_p = 1.87083\times0.72801^2 = 0.99153$. Both pairs sum to 1. Unpolarized light reflects $(R_s + R_p)/2 = 0.05024$, and the reflected light has degree of polarization $(R_s - R_p)/(R_s + R_p) = 0.8315$: 83% polarized at $45°$, 100% at $\theta_B = 56.31°$. (With Hecht's sign convention $r_p = +0.09201$; §8.)

### Example 4.5: The skin depth of copper at 60 Hz and 1 GHz (project Part 7)

$\sigma = 5.96\times10^7$ S/m, $\mu = \mu_0$. Check the good-conductor condition at 1 GHz: $\sigma/(\varepsilon_0\omega) = 5.96\times10^7/(8.854\times10^{-12}\times6.283\times10^9) = 1.07\times10^9\gg1$, more so at 60 Hz. By (3.21):

60 Hz: $\mu_0\sigma\omega = 4\pi\times10^{-7}\times5.96\times10^7\times377.0 = 2.824\times10^4$, $\delta = \sqrt{2/2.824\times10^4} = 8.42$ mm.
1 GHz: $\mu_0\sigma\omega = 4.706\times10^{11}$, $\delta = 2.062\ \mu$m.

The project computes $\kappa$ from the exact (3.20) and agrees with (3.21) to $10^{-6}$; the phase of $\tilde k$ is $45.0000°$. At 1 GHz, $cB_0/E_0 = c|\tilde k|/\omega = 3.27\times10^4$ and the reflectance is $R = 0.999914$. Consequences: mains conductors much thicker than $2\delta$ waste their centres; microwave parts need only a few micrometres of silver or gold plating, because that is where the current flows.

### Example 4.6: An X-band waveguide at 10 GHz (project Part 8)

$a = 2.3$ cm, $b = 1.0$ cm. Cutoffs $f_{mn} = \tfrac c2\sqrt{(m/a)^2 + (n/b)^2}$: TE$_{10}$ 6.517 GHz, TE$_{20}$ 13.03, TE$_{01}$ 14.99, TE$_{11}$ (and TM$_{11}$) 16.35 GHz. At 10 GHz only TE$_{10}$ propagates: the guide is single-mode from 6.52 to 13.03 GHz, which is why the aspect ratio is about 2:1.

$\sqrt{1 - (6.517/10)^2} = \sqrt{1 - 0.42473} = 0.75846$. By (3.33): $v_p = c/0.75846 = 1.3185c = 3.953\times10^8$ m/s, $v_g = 0.7585c = 2.274\times10^8$ m/s, and $v_pv_g = c^2$. The guide wavelength is $\lambda_g = 2\pi/k = \lambda_0/0.75846 = 3.953$ cm (2.998 cm in free space), $k = 158.96$ rad/m, and a pulse takes 4.398 ns per metre instead of 3.336.

## 5. Why it matters and where it connects

- **EM-8.** The Rayleigh dipole and the antenna teaser become retarded potentials and the Larmor formula; Brewster's dipole argument becomes the $\sin^2\theta$ radiation pattern.
- **QM-1, QM-2, QM-3.** Radiation pressure (3.8) is photon momentum $\hbar k$; blackbody mode counting uses (3.32)-type modes in a 3D box; (3.30) is the free particle's spreading; frustrated TIR is tunnelling.
- **QC-1.** The Jones vector is a qubit and the Poincaré sphere the Bloch sphere; wave plates are $S$ and Hadamard gates (§3.6).
- **F3.** (3.28) is the Fourier solution of a linear time-invariant system: $e^{ik(\omega)z}$ is the transfer function of a slab of glass and $z\,dk/d\omega$ its group delay.
- **CM-6.** The loaded string's $\omega(\kappa)$ (CM-6 Example 4.5) is a dispersive medium with $v_g < v_p$ and a band edge.
- **Engineering.** Waveguides and coax in radar and satellite links, optical coatings (Exercise B8), fibre links whose bit rate is limited by (3.30), skin-depth heating and shielding, ionospheric radio and GPS delay corrections (Exercise B5).

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills: a plane wave's fields and intensity, a laser pointer's field and push, polarizers and wave plates, water's $R$, $\theta_B$, $\theta_c$, and seawater's skin depth. B1–B8 are the module: B1 (energy and impedance at a string junction) sets up the analogy; B2 (the s-polarization Fresnel derivation in full) and B4 (the evanescent wave carries no energy) are essential, as are B5 (plasma dispersion and the ionosphere), B6 (waveguide modes), B7 (coax impedance) and B8 (the anti-reflection coating, derived and checked in `project.py`). C1 (superluminal group velocity), C2 (Fresnel's rhomb) and C3 (QM-2's spreading packet) go further.

## 7. Coding project

`project.py` runs in about 4 seconds; figures go to `figures/`.

1. **Part 1, plane waves:** SymPy checks Theorem 3.3's wave against all four Maxwell equations, Lemma 3.2, and (3.4); Example 4.1.
2. **Part 2, Jones calculus:** Malus, three polarizers, Example 4.2, QWP $= S$, HWP $=$ Hadamard; the Poincaré sphere with rotating-plate trajectories (`poincare.png`).
3. **Part 3, Fresnel:** (3.16)–(3.17) with complex $\cos\theta_T$; $R + T = 1$ at every angle; all four interface conditions checked directly; Brewster zero, $|r| = 1$ beyond $\theta_c$, (3.19) (`fresnel.png`).
4. **Part 4, FDTD with matter:** EM-6's Yee leapfrog in grid units plus an **auxiliary differential equation** for the polarization, $\ddot P + \gamma\dot P + \omega_0^2P = \omega_p^2E$, coupled through $\partial_tD = -\partial_xB$. One function covers dielectrics, Lorentz media and plasmas; a glass half-space reproduces Example 4.3.
5. **Part 5, dispersion:** a chirped pulse in a Lorentz medium, widths from FDTD, from FFT propagation of (3.28), and from (3.30); the measured $v_g$ against (3.29) (`dispersive_fdtd.png`).
6. **Parts 6–8:** the two-resonance Lorentz toy and Cauchy fit (`lorentz_model.png`), plasma frequencies, Example 4.5, SymPy verification of the TE fields, mode plots, dispersion curves, $v_p$ and $v_g$ (`waveguide.png`), Example 4.6, coax.
7. **Your turn** (reference solutions for 1 and 3): the quarter-wave anti-reflection coating in FDTD, frustrated TIR, a plasma mirror, the opposite chirp. Part 9 prints every number in `solutions.md`.

**Build it yourself:** before reading `fdtd`, add the polarization to your own EM-6 leapfrog: derive the centred-difference update for $P^{n+1}$ from (3.23) (the factor $1 + \gamma\Delta t/2$ comes from averaging $\dot P$), subtract $P^{n+1} - P^n$ in the $E$ update, and check that a static field gives $\varepsilon(0) = \varepsilon_\infty + \omega_p^2/\omega_0^2$. Then write `fresnel` and find the Brewster zero.

## 8. Where you'll struggle

- **The sign of $i$.** You will mix $e^{-i\omega t}$ and $e^{j\omega t}$ within a page, especially with engineering figures. The tells: engineering permittivity is $\varepsilon' - j\varepsilon''$ (physics $\varepsilon' + i\varepsilon''$), the TIR root is $-j\kappa$ (physics $+i\kappa$), a slow axis gets $e^{-j\Gamma}$ (physics $e^{+i\Gamma}$). NumPy's FFT synthesizes with $e^{+i\omega t}$, so the project's propagator conjugates (`apply_filter`). Fix: write the convention at the top of the page and run one sanity check per calculation: does the evanescent wave decay? Is $\operatorname{Im}\tilde\varepsilon > 0$ for an absorber?
- **Four boundary conditions, two equations.** The Fresnel derivation looks over-determined until you see that the normal conditions repeat the tangential ones (§3.10). Impose $\mathbf E^\parallel$ and $\mathbf H^\parallel$; use $D^\perp$ as a check. Using $B^\perp$ for p light gives $0 = 0$ and the feeling that an equation went missing.
- **Fresnel signs and handedness.** Griffiths (and this module) make $r_p = r_s$ at normal incidence; Hecht and Born–Wolf make $r_p = -r_s$ there, keeping a right-handed $(\hat{\mathbf e}_p, \hat{\mathbf e}_s, \hat{\mathbf k})$ frame for every beam. Mixing them corrupts phase differences (C2 uses $-r_p$ for this reason). "Right circular" is clockwise to the receiver in Hecht and the opposite in IEEE usage; describe the rotation in words before naming it.
- **$T\ne|t|^2$.** Transmitted power carries $n_2\cos\theta_T/(n_1\cos\theta_I)$. If your $R + T\ne1$, this is the reason nine times in ten.
- **Why the sky is blue.** The right answer has two parts: the dipole is driven far below resonance, so its amplitude does not depend on $\omega$, and radiated power goes as $\omega^4$. If you can't say both, you don't have it.
- **$v_g > c$.** It feels like it must break relativity. It doesn't (end of §3.14). Do C1 and compare the absorption length with any distance over which a "superluminal" advance could be seen.
- **Waveguide indexing.** Some books put $a$ on the short side or call the dominant mode TE$_{01}$; TM$_{mn}$ needs both indices nonzero; a circular guide's dominant mode is TE$_{11}$. Ask where the half-wavelengths fit.
- **FDTD with dispersion.** The polarization equation adds a time scale: you need $\omega_0\Delta t\ll1$ for accuracy (the centred oscillator update is unstable for $\omega_0\Delta t > 2$) on top of EM-6's Courant condition. Mur's boundary absorbs only in vacuum at normal incidence, so end the grid in vacuum or make it long. Grid dispersion accumulates: at depth 3000 cells the project's FDTD width is 2% above the FFT width, and it took $S = 0.9$ instead of $0.5$ to get that close.

## 9. Check yourself

1. Derive d'Alembert's solution. Why may you use complex fields in Maxwell's equations but not in $\mathbf S$?
2. Show from Maxwell's equations that a vacuum plane wave is transverse with $\tilde{\mathbf B}_0 = \hat{\mathbf k}\times\tilde{\mathbf E}_0/c$.
3. Derive $I = \tfrac12c\varepsilon_0E_0^2$ and the pressure on a mirror.
4. Write the Jones matrices of a polarizer at $\theta$ and a QWP with fast axis along $x$. What does the QWP do to $45°$ light?
5. Derive $r$, $t$ at normal incidence and show $R + T = 1$. When does the reflected $\mathbf E$ flip?
6. Derive Snell's law from phase matching.
7. For p polarization, which two boundary conditions are independent, and why are the other two automatic?
8. Derive $\tan\theta_B = n_2/n_1$ and explain the polarization of glare.
9. Derive the evanescent decay length. Which root of $\cos\theta_T$, and why?
10. Derive $\delta = \sqrt{2/\mu\sigma\omega}$ and estimate it for copper at 60 Hz.
11. Derive $\tilde\varepsilon_r(\omega)$ for bound electrons and $n^2 = 1 - \omega_p^2/\omega^2$ for free ones. Where is dispersion anomalous?
12. Derive $v_g = d\omega/dk$, the Gaussian width $\tau(z)$, and the TE$_{mn}$ cutoff; prove a hollow guide has no TEM mode.

## 10. Pair with

- **Griffiths, *Introduction to Electrodynamics*, ch. 9.** The same path (strings, plane waves, Fresnel for p, conductors, Lorentz dispersion, waveguides, coax) with more problems; §9.3.3's figure is the p geometry used here.
- **Hecht, *Optics*, ch. 2–4 and 8.** The optics side done properly: Fresnel's equations with many plots, TIR and Fresnel's rhomb, wave plates and polarization by every mechanism. Watch for his $r_p$ sign.
- **Feynman, *Lectures on Physics* Vol. I, ch. 28–33.** Refraction derived from the radiation of driven electrons (ch. 31), which explains *why* a slower wave appears; the best physical account of the index of refraction.
- **MIT 8.07 lecture notes, lectures 20–24.** A compact second pass at a more mathematical level, including dispersion and causality.
- **MIT 8.03 (Vibrations and Waves), OCW materials.** The wave fundamentals (strings, impedance matching, dispersion, group velocity) with demonstrations, if §3.1–3.3 felt fast.
