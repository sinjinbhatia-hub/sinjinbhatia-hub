# QM-1 Solutions

Every number below is printed by `project.py` (Parts 2, 6 and 8). $hc = 1239.84$ eV nm.

## A. Warm-up

**A1.** $E = hc/\lambda$. At 400 nm: $1239.84/400 = 3.100$ eV. At 700 nm: $1239.84/700 = 1.771$ eV. At 1 MHz: $E = hf = (6.626\times10^{-34})(10^6)/(1.602\times10^{-19}) = 4.136\times10^{-9}$ eV.

Only the 400 nm photon exceeds $W = 2.3$ eV. Threshold: $\lambda_0 = hc/W = 1239.84/2.3 = 539.1$ nm.

The radio transmitter: each electron absorbs one photon at a time (3.3). One radio photon delivers $4\times10^{-9}$ eV, about $10^{-9}$ of what is needed. The number of photons sets how many absorption events happen, not how much energy each one delivers. Two photons can be absorbed by one electron at once, but at these intensities and energies the probability is negligible (it needs about $5\times10^8$ simultaneous photons).

**A2.** $E_n = -13.606/n^2$ eV. $E_2 = -13.606/4 = -3.401$ eV; $E_3 = -13.606/9 = -1.512$ eV.

Lyman-$\alpha$: $1/\lambda = R_\infty(1 - \frac14) = \frac34R_\infty$, $\lambda = 4/(3\times1.097373\times10^7)$ m $= 121.50$ nm.

H$\beta$: $1/\lambda = R_\infty(\frac14 - \frac1{16}) = \frac{3}{16}R_\infty$, $\lambda = 16/(3R_\infty) = 486.01$ nm.

Ionization from $n = 2$: $0 - E_2 = 3.401$ eV (a 364.5 nm photon, the Balmer series limit).

**A3.** (a) $\hat S_x|{-x}\rangle = \frac\hbar2\begin{pmatrix}0&1\\1&0\end{pmatrix}\frac1{\sqrt2}\begin{pmatrix}1\\-1\end{pmatrix} = \frac\hbar2\frac1{\sqrt2}\begin{pmatrix}-1\\1\end{pmatrix} = -\frac\hbar2|{-x}\rangle$. Eigenvalue $-\hbar/2$.

(b) $\langle{+x}|{\downarrow}\rangle = \frac1{\sqrt2}(1, 1)(0, 1)^\top = \frac1{\sqrt2}$, so $P(S_x = +\hbar/2) = \frac12$.

(c) $|\chi\rangle = \frac1{\sqrt2}(1, i)^\top$ is $|{+y}\rangle$ by (3.14).
- $S_x$: $\langle{\pm x}|\chi\rangle = \frac12(1\pm i)$, $|\cdot|^2 = \frac14(1 + 1) = \frac12$ each.
- $S_y$: $\langle{+y}|\chi\rangle = \frac12(\bar1\cdot1 + \bar i\cdot i) = \frac12(1 + (-i)(i)) = \frac12(1 + 1) = 1$, so $P(+) = 1$; $\langle{-y}|\chi\rangle = \frac12(1 + \overline{(-i)}\,i) = \frac12(1 + i\cdot i) = 0$, so $P(-) = 0$. The conjugate on the bra is essential: without it you would get $\frac12(1 + i^2) = 0$ for $P(+)$, which is wrong.
- $S_z$: $|\langle{\uparrow}|\chi\rangle|^2 = \frac12$, $|\langle{\downarrow}|\chi\rangle|^2 = |i/\sqrt2|^2 = \frac12$.

**A4.** (a) For any observable with eigenvectors $|a_i\rangle$: $|\langle a_i|e^{i\alpha}\psi\rangle|^2 = |e^{i\alpha}|^2|\langle a_i|\psi\rangle|^2 = |\langle a_i|\psi\rangle|^2$, using linearity in the second slot and $|e^{i\alpha}| = 1$. After the outcome $a_i$, (3.10) gives $|a_i\rangle\langle a_i|e^{i\alpha}\psi\rangle/|\langle a_i|e^{i\alpha}\psi\rangle| = e^{i\alpha}\times(\text{the state from }|\psi\rangle)$, which differs only by the same global phase. Since every future probability is computed by the same two operations (and unitary evolution is linear, so the phase rides along), no experiment can detect $\alpha$.

(b) $|{+x}\rangle = \frac1{\sqrt2}(1, 1)^\top$ and $|{-x}\rangle = \frac1{\sqrt2}(1, -1)^\top$ give identical $S_z$ probabilities ($\frac12$, $\frac12$) but an $S_x$ measurement distinguishes them perfectly: $P(S_x = +\hbar/2) = 1$ for the first and 0 for the second. A relative phase of $e^{i\pi} = -1$ between components is not a global phase: you cannot remove it by multiplying the whole vector by one number.

**A5.** $k_BT = (1.381\times10^{-23})(300)/(1.602\times10^{-19}) = 0.02585$ eV. $p = \sqrt{2m_nk_BT}$ with $m_n = 1.675\times10^{-27}$ kg: $p = \sqrt{2(1.675\times10^{-27})(4.142\times10^{-21})} = 3.725\times10^{-24}$ kg m/s, so $\lambda = h/p = 0.1779$ nm. That sits inside the 0.1–0.3 nm range of interatomic spacings, so crystals diffract thermal neutrons strongly (the Bragg condition $2d\sin\theta = \lambda$ has solutions at sizeable angles). Neutrons are neutral and scatter off nuclei, so they see light atoms (hydrogen) and magnetic order that X-rays miss.

## B. Core

**B1.** (a) The electron's kinetic energy is $hc/\lambda - hc/\lambda'$, largest when $\lambda'$ is largest. By (3.7), $\lambda'$ is largest when $1 - \cos\theta$ is largest: $\theta = 180°$ (back-scattering), where $\lambda' = \lambda + 2h/m_ec$.

(b) $\lambda = 1239.84/(10^5)$ nm $= 12.398$ pm. $\lambda' = 12.398 + 2(2.4263) = 17.251$ pm. $E' = 1239.84/0.017251$ eV $= 71.87$ keV. Electron kinetic energy: $100 - 71.87 = 28.13$ keV, the "Compton edge" of a 100 keV gamma-ray spectrum.

(c) In energies, (3.7) at $\theta = 180°$ is $\frac{1}{E'} = \frac1E + \frac{2}{m_ec^2}$ (divide (3.7) by $hc$). As $E\to\infty$, $1/E'\to2/(m_ec^2)$, so $E'\to m_ec^2/2 = 255.5$ keV. However energetic the incoming photon, the back-scattered one carries at most about 256 keV.

**B2.** (a) Let the electron be at distance $r_e$ and the nucleus at $r_N$ from the center of mass, with $r = r_e + r_N$ the separation and $mr_e = Mr_N$, so $r_e = Mr/(m + M)$, $r_N = mr/(m + M)$. Both circle at angular velocity $\Omega$.
1. Newton for the electron: $m\Omega^2r_e = \kappa/r^2$, i.e. $\mu\Omega^2r = \kappa/r^2$ with $\mu = mM/(m + M)$ (substitute $r_e$).
2. Total angular momentum: $L = m\Omega r_e^2 + M\Omega r_N^2 = \Omega\frac{mM^2 + Mm^2}{(m + M)^2}r^2 = \mu\Omega r^2$. Set $L = n\hbar$.
3. Steps 1 and 2 are those of Theorem 3.1 with $m\to\mu$ and $v\to\Omega r$: from 2, $\Omega = n\hbar/(\mu r^2)$; into 1, $\mu n^2\hbar^2/(\mu^2r^3) = \kappa/r^2$, so $r = n^2\hbar^2/(\mu\kappa)$.
4. Energy: kinetic $\frac12m\Omega^2r_e^2 + \frac12M\Omega^2r_N^2 = \frac12\mu\Omega^2r^2 = \frac12\kappa/r$ by step 1, so $E = -\kappa/(2r) = -\mu\kappa^2/(2\hbar^2n^2)$. Hence $E_R\to E_R\,\mu/m$ and $R_M = R_\infty\mu/m = R_\infty/(1 + m/M)$.

(b) $m/m_p = 1/1836.15$, $m/m_d = 1/3670.48$. $R_H = 1.097373\times10^7/(1 + 1/1836.15) = 1.096776\times10^7$ m$^{-1}$; $R_D = 1.097074\times10^7$ m$^{-1}$. H$\alpha$: $36/(5R_H) = 656.470$ nm; D$\alpha$: $36/(5R_D) = 656.291$ nm. Shift: $0.179$ nm, deuterium on the short-wavelength side. Urey, Brickwedde and Murphy saw faint companions displaced by this amount after concentrating the heavy isotope by evaporating liquid hydrogen.

**B3.** (a) In the box $|p| = \sqrt{2mE}$ is constant. One cycle goes $0\to L$ with $p = +|p|$ and back with $p = -|p|$, so $\oint p\,dq = |p|L + (-|p|)(-L) = 2|p|L$. Setting $2|p|L = nh$: $|p| = nh/(2L)$ and $E_n = p^2/2m = n^2h^2/(8mL^2)$.

(b) The orbit $H = E$ is an ellipse in the $(q, p)$ plane with semi-axes $q_0 = \sqrt{2E/(m\omega^2)}$ and $p_0 = \sqrt{2mE}$. $\oint p\,dq$ is its area, $\pi q_0p_0 = \pi\sqrt{2E/(m\omega^2)}\sqrt{2mE} = 2\pi E/\omega$ (`project.py` Part 2 checks this integral in SymPy). Setting it to $nh$: $E_n = nh\omega/(2\pi) = n\hbar\omega$.

(c) The box is exactly right, provided $n\ge1$ ($n = 0$ would mean $p = 0$, a particle at rest; the rule does not say to exclude it, the wave picture does: a standing wave that vanishes at both walls must have at least one half-wavelength in the box). The oscillator has the right spacing $\hbar\omega$ but every level is too low by $\frac12\hbar\omega$, and the rule allows $E = 0$, a particle sitting motionless at the bottom of the well, which QM-2's uncertainty principle forbids. The box works because its walls reflect the wave with a fixed phase that happens to compensate exactly; the smooth turning points of the oscillator contribute a phase loss of $\pi/2$ each, which is the $\frac12$ that WKB (QM-8) restores: $\oint p\,dq = (n + \frac12)h$.

(d) The rule's ground state ($n = 0$) has $E = 0$; the correct one is $\frac12\hbar\omega$. The ratio is 0: the rule misses the zero-point energy entirely. (If you insist on starting at $n = 1$, the rule's lowest level $\hbar\omega$ is twice the true ground state.)

**B4.** (a) *Projector identity.* $P = \frac12(I + \hat{\mathbf n}'\cdot\boldsymbol\sigma)$ satisfies $P|{\pm n'}\rangle = \frac12(1\pm1)|{\pm n'}\rangle$, since $|{\pm n'}\rangle$ are eigenvectors of $\hat{\mathbf n}'\cdot\boldsymbol\sigma$ with eigenvalues $\pm1$ (LA-8 Proposition 3.22). So $P$ fixes $|{+n'}\rangle$ and kills $|{-n'}\rangle$; $|{+n'}\rangle\langle{+n'}|$ does the same; two matrices that agree on a basis are equal. *Probability.* $|\langle{+n'}|\psi\rangle|^2 = \langle\psi|{+n'}\rangle\langle{+n'}|\psi\rangle = \langle\psi|P|\psi\rangle = \frac12(1 + \hat{\mathbf n}'\cdot\langle\psi|\boldsymbol\sigma|\psi\rangle) = \frac12(1 + \hat{\mathbf n}'\cdot\hat{\mathbf n})$, using $\langle\boldsymbol\sigma\rangle = \hat{\mathbf n}$ from Theorem 3.5. With $\gamma$ the angle between the axes this is $\frac12(1 + \cos\gamma) = \cos^2(\gamma/2)$.

(b) Two states with Bloch vectors $\hat{\mathbf n}$, $\hat{\mathbf n}'$ are orthogonal iff $|\langle{+n'}|\psi\rangle|^2 = 0$ iff $\hat{\mathbf n}\cdot\hat{\mathbf n}' = -1$ iff $\hat{\mathbf n}' = -\hat{\mathbf n}$.

(c) $a = 1/\sqrt5$ is already real and positive, so $\cos(\theta/2) = 1/\sqrt5$, $e^{i\phi}\sin(\theta/2) = 2i/\sqrt5$, giving $\sin(\theta/2) = 2/\sqrt5$ and $\phi = 90°$. $\theta = 2\arccos(1/\sqrt5) = 126.87°$. $\cos\theta = 2\cos^2(\theta/2) - 1 = \frac25 - 1 = -0.6$; $\sin\theta = 2\cdot\frac1{\sqrt5}\cdot\frac2{\sqrt5} = 0.8$. So $\hat{\mathbf n} = (0.8\cos90°, 0.8\sin90°, -0.6) = (0, 0.8, -0.6)$: in the $yz$-plane, below the equator.

**B5.** Let $\theta = 60°$, so $\theta/2 = 30°$ and $|{\pm n}\rangle = (\cos30°, \sin30°)^\top$, $(\sin30°, -\cos30°)^\top$ (Proposition 3.22 with $\phi = 0$).

(a) $P(+n) = |\langle{+n}|{\uparrow}\rangle|^2 = \cos^230° = 0.75$. Afterwards the state is $|{+n}\rangle$.

(b) Through $+n$: $P(+n)\cdot|\langle{\downarrow}|{+n}\rangle|^2 = 0.75\times\sin^230° = 0.75\times0.25 = 0.1875$. Through $-n$: $P(-n) = \sin^230° = 0.25$, then $|\langle{\downarrow}|{-n}\rangle|^2 = \cos^230° = 0.75$, product $0.1875$. Either port, with path recorded: the two routes are distinguishable, so probabilities add: $0.375$.

(c) With the beams recombined and no record, the routes are indistinguishable and amplitudes add (rule ii): the amplitude to end in $|{\downarrow}\rangle$ is $\langle{\downarrow}|{+n}\rangle\langle{+n}|{\uparrow}\rangle + \langle{\downarrow}|{-n}\rangle\langle{-n}|{\uparrow}\rangle = \sin30°\cos30° + (-\cos30°)\sin30° = 0$. Equivalently the analyzer acts as $\sum_k|{k}\rangle\langle k| = I$. Probability 0. In (b) the detectors make the routes distinguishable, so probabilities add (rule iii) and the interference term, which is exactly $-0.375$, is gone.

(d) `project.py` Part 3, seed 0: 1874/10000 = 0.1874 through $+n$ then $z-$, and 3752/10000 = 0.3752 through either. Both match within 4 binomial standard deviations ($\sigma\approx39$ and $48$).

**B6.** (a) $\sigma_y\sigma_z = \begin{pmatrix}0&-i\\i&0\end{pmatrix}\begin{pmatrix}1&0\\0&-1\end{pmatrix} = \begin{pmatrix}0&i\\i&0\end{pmatrix} = i\sigma_x$; $\sigma_z\sigma_y = \begin{pmatrix}1&0\\0&-1\end{pmatrix}\begin{pmatrix}0&-i\\i&0\end{pmatrix} = \begin{pmatrix}0&-i\\-i&0\end{pmatrix} = -i\sigma_x$. So $[\sigma_y, \sigma_z] = 2i\sigma_x$ and $[\hat S_y, \hat S_z] = \frac{\hbar^2}4\cdot2i\sigma_x = i\hbar\hat S_x$.

$\sigma_z\sigma_x = \begin{pmatrix}1&0\\0&-1\end{pmatrix}\begin{pmatrix}0&1\\1&0\end{pmatrix} = \begin{pmatrix}0&1\\-1&0\end{pmatrix} = i\sigma_y$ (since $i\sigma_y = \begin{pmatrix}0&1\\-1&0\end{pmatrix}$); $\sigma_x\sigma_z = \begin{pmatrix}0&-1\\1&0\end{pmatrix} = -i\sigma_y$. So $[\hat S_z, \hat S_x] = i\hbar\hat S_y$.

(b) $\hat S_x^2 = \frac{\hbar^2}4\sigma_x^2 = \frac{\hbar^2}4I$, so $\langle\hat S_x^2\rangle = \hbar^2/4$ in every state. By Theorem 3.5, $\langle\hat S_x\rangle = \frac\hbar2n_x$. So $(\Delta S_x)^2 = \frac{\hbar^2}4 - \frac{\hbar^2}4n_x^2 = \frac{\hbar^2}4(1 - n_x^2)$, and the same for $y$.

(c) $\frac14|\langle[\hat S_x, \hat S_y]\rangle|^2 = \frac14|i\hbar\langle\hat S_z\rangle|^2 = \frac14\hbar^2\frac{\hbar^2}{4}n_z^2 = \frac{\hbar^4}{16}n_z^2$. The left side is $\frac{\hbar^4}{16}(1 - n_x^2)(1 - n_y^2)$. Using $n_z^2 = 1 - n_x^2 - n_y^2$, the inequality becomes $1 - n_x^2 - n_y^2 + n_x^2n_y^2\ge1 - n_x^2 - n_y^2$, i.e. $n_x^2n_y^2\ge0$. True always; equality iff $n_x = 0$ or $n_y = 0$: Bloch vectors in the $yz$-plane or the $xz$-plane. $|{+x}\rangle$ ($n_y = 0$) gives equality: $\Delta S_x = 0$ and $\langle\hat S_z\rangle = 0$, so both sides are 0. "Minimum uncertainty" here means only the bound is attained, and in this state the bound itself is zero, which shows that Robertson's bound depends on the state and can be uninformative. `project.py` Part 6 checks both closed forms on 2000 random states.

**B7.** (a) $\hat H = \omega_0\hat S_z = \frac{\hbar\omega_0}2\sigma_z$ is diagonal, so $U(t) = e^{-i\hat Ht/\hbar} = \operatorname{diag}(e^{-i\omega_0t/2}, e^{i\omega_0t/2})$ (LA-7 §3.9: the exponential of a diagonal matrix is entrywise). $|\psi(t)\rangle = \frac1{\sqrt2}(e^{-i\omega_0t/2}, e^{i\omega_0t/2})^\top$. Then $\langle{+x}|\psi(t)\rangle = \frac12(e^{-i\omega_0t/2} + e^{i\omega_0t/2}) = \cos(\omega_0t/2)$ and $P(S_x = +\hbar/2) = \cos^2(\omega_0t/2)$.

(b) With $a = e^{-i\omega_0t/2}/\sqrt2$, $b = e^{i\omega_0t/2}/\sqrt2$: $\bar ab = \frac12e^{i\omega_0t}$. By Theorem 3.5, step 4: $\langle\hat S_x\rangle = \frac\hbar2\cdot2\operatorname{Re}(\bar ab) = \frac\hbar2\cos\omega_0t$, $\langle\hat S_y\rangle = \frac\hbar2\sin\omega_0t$, $\langle\hat S_z\rangle = \frac\hbar2(|a|^2 - |b|^2) = 0$. The mean spin rotates counterclockwise (seen from $+z$) in the $xy$-plane at angular frequency $\omega_0$: Larmor precession. Note the half-angle: the state returns to $|{+x}\rangle$ up to a sign after $\omega_0t = 2\pi$, and returns exactly (including sign) only after $4\pi$.

(c) $U^\dagger U = \operatorname{diag}(e^{i\omega_0t/2}e^{-i\omega_0t/2}, e^{-i\omega_0t/2}e^{i\omega_0t/2}) = I$.

(d) $\omega_0 = g\mu_BB/\hbar = 2.0023\times9.274\times10^{-24}\times1/(1.0546\times10^{-34}) = 1.7609\times10^{11}$ rad/s; period $2\pi/\omega_0 = 35.68$ ps (a 28 GHz precession, the frequency of electron spin resonance at 1 T).

**B8.** (a) $\bar\phi_1\phi_2 = A^2e^{ikx}e^{ikx} = A^2e^{2ikx}$. So $2\operatorname{Re}[\bar\phi_1\phi_2\langle t_1|t_2\rangle] = 2A^2r\operatorname{Re}e^{i(2kx + \delta)} = 2A^2r\cos(2kx + \delta)$, and with $|\phi_1|^2 = |\phi_2|^2 = A^2$, (3.16) gives $P = 2A^2[1 + r\cos(2kx + \delta)]$. Over a fringe (with $A$ slowly varying), $P_{\max} = 2A^2(1 + r)$ and $P_{\min} = 2A^2(1 - r)$, so the visibility is $\frac{2r}{2} = r$. The phase $\delta$ shifts the whole fringe pattern sideways by $\delta/(2k)$ without changing its contrast.

(b) The overlap is $\langle t_1|t_2\rangle = (1, 0)\cdot(\cos\beta, e^{i\delta}\sin\beta)^\top = \cos\beta$, so the visibility is $|\cos\beta|$: 1, 0.866, 0.5, 0. A reference implementation, added after the existing loop in Part 5:

```python
for beta in (0, np.pi / 6, np.pi / 3, np.pi / 2):
    t1, t2 = np.array([1, 0]), np.array([np.cos(beta), np.exp(0.4j) * np.sin(beta)])
    ov = np.vdot(t1, t2)                                 # <t1|t2>, conjugating t1
    p1, p2 = amps(xs)
    joint = p1[:, None] * t1[None, :] + p2[:, None] * t2[None, :]
    P = (np.abs(joint)**2).sum(1)                        # sum over tag outcomes, Parseval
    assert np.allclose(P, pattern(xs, ov))
    assert abs(visibility(xs, P) - abs(np.cos(beta))) < 1e-3
```

(c) The overlap $\langle t_1|t_2\rangle$ is a property of the tag system, not of the electron's trajectory: the derivation of (3.16) never changes $\phi_1$ or $\phi_2$. Interference disappears because the two routes end in orthogonal joint states, so their amplitudes can no longer add; any disturbance of the electron is beside the point, and schemes with negligible momentum kick still destroy the fringes when they record the path.

**B9.** $|\psi\rangle = \frac1{\sqrt5}(1, 2i)^\top$, so $a = 1/\sqrt5$, $b = 2i/\sqrt5$, $\bar ab = 2i/5$.

(a) $\langle\hat S_x\rangle = \frac\hbar2\cdot2\operatorname{Re}(\bar ab) = 0$. $\langle\hat S_y\rangle = \frac\hbar2\cdot2\operatorname{Im}(\bar ab) = \frac\hbar2\cdot\frac45 = 0.4\hbar$. $\langle\hat S_z\rangle = \frac\hbar2(|a|^2 - |b|^2) = \frac\hbar2(\frac15 - \frac45) = -0.3\hbar$. (Bloch vector $(0, 0.8, -0.6)$, as in B4(c).)

(b) $\langle{+y}|\psi\rangle = \frac1{\sqrt2}(1, -i)\cdot\frac1{\sqrt5}(1, 2i)^\top = \frac{1 + 2}{\sqrt{10}} = \frac3{\sqrt{10}}$, $p_+ = 0.9$. $\langle{-y}|\psi\rangle = \frac1{\sqrt{10}}(1 + i\cdot2i) = -\frac1{\sqrt{10}}$, $p_- = 0.1$. $\sum a_ip_i = \frac\hbar2(0.9 - 0.1) = 0.4\hbar$. ✓

(c) $\overline{\langle\psi|\hat A|\psi\rangle} = \langle\hat A\psi|\psi\rangle$ (conjugate symmetry of the inner product) $= \langle\psi|\hat A^\dagger\psi\rangle$ (definition of the adjoint, LA-8 §3.3) $= \langle\psi|\hat A\psi\rangle$ (Hermitian). A number equal to its conjugate is real.

## C. Challenge

**C1.** (a) $u_\lambda\propto\lambda^{-5}(e^{hc/\lambda k_BT} - 1)^{-1}$. In terms of $x = hc/(\lambda k_BT)$, $\lambda\propto1/x$, so $u\propto g(x) = x^5/(e^x - 1)$ and maximizing over $\lambda$ is maximizing $g$ over $x$. $g'(x) = \frac{5x^4(e^x - 1) - x^5e^x}{(e^x - 1)^2} = 0$ gives $5(e^x - 1) = xe^x$, i.e. $x = 5(1 - e^{-x})$. Iterating $x\leftarrow5(1 - e^{-x})$ from $x = 5$: $4.966, 4.9651, \dots$; Brent's method in `project.py` gives $x = 4.965114$. Then $b = \lambda_{\max}T = hc/(k_Bx) = 2.8978\times10^{-3}$ m K $= 2.8978$ mm K.

(b) $\lambda_{\max} = 2.8978\times10^{-3}/5772$ m $= 502.0$ nm (blue-green). At 5000 K, 579.6 nm.

(c) $\int_0^\infty u_\lambda d\lambda$: substitute $\lambda = hc/(xk_BT)$, $d\lambda = -\frac{hc}{k_BT}\frac{dx}{x^2}$. Then $\frac{8\pi hc}{\lambda^5} = 8\pi hc\frac{x^5k_B^5T^5}{h^5c^5}$ and the integral is $8\pi hc\frac{k_B^5T^5}{h^5c^5}\frac{hc}{k_BT}\int_0^\infty\frac{x^3\,dx}{e^x - 1} = \frac{8\pi k_B^4T^4}{h^3c^3}\cdot\frac{\pi^4}{15}$. So $a = 8\pi^5k_B^4/(15h^3c^3) = 7.566\times10^{-16}$ J m$^{-3}$ K$^{-4}$ (`project.py` Part 1 checks the integral numerically to 0.2%). With (3.1), $\int\lambda^{-4}d\lambda$ diverges as $\lambda^{-3}$ at the short-wavelength end: the ultraviolet catastrophe.

**C2.** (a) *Repeatability:* an atom with $m_z > 0$ passes SGz$+$ and, since $\mathbf m$ is unchanged, gives $+$ again. *Experiment 2:* given $m_z > 0$, $\mathbf m$ is uniform on the upper hemisphere, which is symmetric under $m_x\to-m_x$, so $P(m_x > 0) = \frac12$.

(b) Given $m_z > 0$, we need the fraction of the upper hemisphere with $\mathbf m\cdot\hat{\mathbf n} > 0$. The set $\{\mathbf m\cdot\hat{\mathbf z} > 0\}\cap\{\mathbf m\cdot\hat{\mathbf n} > 0\}$ is a spherical lune bounded by two great circles whose planes meet at angle $\pi - \theta$. A lune of dihedral angle $\alpha$ is the fraction $\alpha/2\pi$ of the unit sphere, so its area is $4\pi\cdot\alpha/2\pi = 2\alpha$. The two hemispheres' normals are at angle $\theta$, so their intersection is a lune of angle $\pi - \theta$ (check: $\theta = 0$ gives the whole hemisphere, angle $\pi$), with area $2(\pi - \theta)$. The upper hemisphere has area $2\pi$. So $P(+n\mid z+) = (\pi - \theta)/\pi = 1 - \theta/\pi$.

Comparison (`project.py` Part 8, $2\times10^5$ samples):

| $\theta$ | model $1 - \theta/\pi$ | Monte Carlo | quantum $\cos^2(\theta/2)$ |
|---|---|---|---|
| 60° | 0.6667 | 0.6670 | 0.7500 |
| 90° | 0.5000 | 0.4990 | 0.5000 |
| 120° | 0.3333 | 0.3337 | 0.2500 |

The model agrees with quantum mechanics only at $0°$, $90°$ and $180°$, the angles in experiments 1–2. At $60°$ it is off by 8 percentage points, easily measurable.

(c) The code (as in Part 8):

```python
m = rng.normal(size=(200000, 3)); m /= np.linalg.norm(m, axis=1)[:, None]   # uniform on sphere
mz = m[m[:, 2] > 0]                                                         # passed SGz+
for th in (np.pi / 3, np.pi / 2, 2 * np.pi / 3):
    n = np.array([np.sin(th), 0, np.cos(th)])
    print(th, (mz @ n > 0).mean(), 1 - th / np.pi, np.cos(th / 2)**2)
```

Normalizing Gaussian vectors gives the uniform distribution on the sphere because the 3D Gaussian is rotation-invariant.

(d) With $\mathbf m$ unchanged, experiment 3 gives SGz$+$, SGx$+$, SGz: always $+$, because $m_z > 0$ is still true. The data say 50/50. To rescue the model, let each SGx re-randomize $\mathbf m$ uniformly over the hemisphere $\mathbf m\cdot\hat{\mathbf x} > 0$ selected by the outcome. Then experiment 3 gives 50/50. But (b) still fails, so you would also have to change the response rule (for example, outcome probabilities depending on $\mathbf m$ in a tuned way). Such single-particle models can be built (Bell 1966 gave one for spin-1/2), so one particle alone cannot refute hidden variables. What cannot be built, Bell showed in 1964, is a *local* hidden-variable model reproducing the correlations of two entangled spins measured far apart; QM-7 derives the inequality and the experimental violation.

**C3.** (a) From Theorem 3.1, $v_n = n\hbar/(mr_n)$ and $r_n = n^2a_0$, so $f_{\rm orb} = \frac{v_n}{2\pi r_n} = \frac{n\hbar}{2\pi mr_n^2} = \frac{\hbar}{2\pi ma_0^2n^3}$. Since $E_R = \kappa/(2a_0)$ and $a_0 = \hbar^2/(m\kappa)$, $E_R = \hbar^2/(2ma_0^2)$, so $\hbar/(ma_0^2) = 2E_R/\hbar$ and $f_{\rm orb} = \frac{2E_R}{2\pi\hbar n^3} = \frac{2E_R}{hn^3}$.

(b) $f_{n\to n-1} = \frac{E_R}{h}\Big(\frac1{(n-1)^2} - \frac1{n^2}\Big) = \frac{E_R}{h}\frac{2n - 1}{n^2(n - 1)^2}$. Ratio: $\frac{f_{n\to n-1}}{f_{\rm orb}} = \frac{(2n - 1)n}{2(n - 1)^2} = \frac{2n^2 - n}{2n^2 - 4n + 2}\to1$ as $n\to\infty$; expanding, the ratio is $1 + \frac{3}{2n} + O(n^{-2})$. At $n = 10$: $\frac{190}{162} = 1.17284$. At $n = 100$: $\frac{19900}{19602} = 1.01520$.

(c) Classically (EM-8 §3.8–3.9), an electron orbiting at frequency $f_{\rm orb}$ radiates at $f_{\rm orb}$ (and harmonics, for elliptical orbits). For large orbits, the atom is a macroscopic object and classical electrodynamics must be right. So a quantum theory must make the emitted frequency for $n\to n - 1$ approach the orbital frequency for large $n$ (and $n\to n - k$ approach the $k$th harmonic). Bohr used exactly this demand in 1913 to fix the constant in his model, and later made it a general principle, the **correspondence principle**. The ratio $1 + 3/(2n)$ shows how slowly the classical limit is approached.
