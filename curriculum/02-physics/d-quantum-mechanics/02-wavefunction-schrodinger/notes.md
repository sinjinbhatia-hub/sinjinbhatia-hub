# QM-2: The Wavefunction and the Schrödinger Equation: Probability, Operators, Expectation Values, Uncertainty, Momentum Space

**Time:** ~12 hours · **Prerequisites:** QM-1, F3 · **Unlocks:** QM-3

> **Conventions in this module.** $\hbar = h/2\pi$ is written explicitly everywhere. $\Psi(x,t)$ is the wavefunction of one particle of mass $m$ on the line; $|\Psi(x,t)|^2\,dx$ is the probability of finding it in $[x, x+dx]$ at time $t$, and $\int_{-\infty}^\infty|\Psi|^2dx = 1$. $V(x)$ is the potential energy, always **real**. The $L^2$ inner product is $\langle f, g\rangle = \int\overline{f(x)}\,g(x)\,dx$, conjugate on the **first** slot (LA-8, F3). Operators carry hats: $\hat x f = xf$, $\hat p = -i\hbar\,\partial/\partial x$, $\hat H = \hat p^2/2m + V(\hat x)$. $\langle Q\rangle$ is an expectation value, $\sigma_Q^2 = \langle Q^2\rangle - \langle Q\rangle^2$. Time dependence uses the **physics sign** $e^{-i\omega t}$, $e^{-iEt/\hbar}$ (F1 §3.8 warning): a plane wave is $e^{i(kx-\omega t)}$. Momentum space uses the symmetric convention $\Phi(p) = (2\pi\hbar)^{-1/2}\int\Psi(x)e^{-ipx/\hbar}dx$; in $k = p/\hbar$, $\phi(k) = (2\pi)^{-1/2}\int\Psi e^{-ikx}dx$. Relative to F3's $\hat f(\xi) = \int f e^{-2\pi ix\xi}dx$: $\Phi(p) = (2\pi\hbar)^{-1/2}\,\hat\Psi\big(p/2\pi\hbar\big)$, i.e. $p = 2\pi\hbar\,\xi$ (§3.5). Simulations in `project.py` use $\hbar = m = 1$.

## 0. Where this sits

QM-1 told you which experiments broke classical physics and previewed the postulates as "a state is a unit vector, an observable is a Hermitian operator." This module makes that concrete for one particle on a line: the state is a complex function $\Psi(x,t)$, its evolution is the Schrödinger equation, and the momentum operator, the uncertainty principle, and the energy eigenvalue problem all come out of it. You need F3's Fourier transform, Plancherel (F3 (3.15)) and uncertainty theorem (F3 Theorem 3.11), LA-8's adjoint and Hermitian operators, and CM-8's Hamilton–Jacobi equation (CM-8 (3.28)).

QM-3 solves $\hat H\psi = E\psi$ (§3.9 here) for wells, oscillators and barriers. QM-4 rewrites everything here in Dirac notation, where §3.6's "same state in two bases" becomes the whole theory. QM-8's WKB method starts from §3.11.

## 1. The problem that forced this idea

By 1925 de Broglie's relation $\lambda = h/p$ (QM-1) had a wave for every particle and no equation for the wave. Peter Debye, after Schrödinger presented de Broglie's thesis at a Zürich colloquium in late 1925, remarked that talking about waves without a wave equation was childish; that is Felix Bloch's recollection, fifty years later (*Physics Today*, December 1976). Within weeks Schrödinger had one.

He did not get it by guessing. He got it from the analogy you met in CM-8 §3.10. In optics, the eikonal equation of geometric optics is the short-wavelength limit of the wave equation: rays are what waves do when the wavelength is small compared with everything else. Hamilton had noticed in the 1830s that the Hamilton–Jacobi equation for a particle has the same form as the eikonal equation, with trajectories as rays and surfaces of constant action $S$ as wavefronts. Nobody needed a "wave mechanics" behind classical mechanics because classical mechanics worked. In 1926 it no longer did. Schrödinger's argument, set out in the second of his four papers *Quantisierung als Eigenwertproblem* (Annalen der Physik, 1926), was: classical mechanics is the geometric-optics limit of some wave theory, as ray optics is the limit of wave optics, so find the wave equation whose short-wavelength limit is Hamilton–Jacobi. The first paper wrote the time-independent equation and got the hydrogen spectrum; the fourth wrote the time-dependent one. §3.11 runs the argument backwards: put $\Psi = \sqrt\rho\,e^{iS/\hbar}$ into the equation and Hamilton–Jacobi falls out, plus one extra term that vanishes as $\hbar\to0$.

What the wave *meant* was unsettled. Schrödinger hoped $|\Psi|^2$ was a smeared-out charge density. That fails: a free wave packet spreads without limit (§3.8), and an electron hitting a screen arrives at one point. Max Born, analyzing scattering in 1926 (*Zur Quantenmechanik der Stoßvorgänge*), proposed that $\Psi$ gives the *probability* of an outcome. In the paper's famous footnote, added in proof, he corrected himself to say the probability is proportional to the *square* of the amplitude. That footnote is the statistical interpretation of §3.2.

Two more pieces arrived in 1927. Paul Ehrenfest showed that expectation values obey Newton's law in the form $d\langle p\rangle/dt = \langle -V'\rangle$ (§3.4). Heisenberg published the uncertainty principle with a heuristic estimate, and Earle Kennard proved the sharp form $\sigma_x\sigma_p \ge \hbar/2$ the same year. You already proved it in F3; §3.5 adds the $\hbar$.

The equation itself is a postulate. Nobody derives it from something more basic; §3.1 *motivates* it and §8 says so again.

## 2. Intuition

**A complex field whose modulus squared is a probability.** $\Psi(x,t)$ is a complex number at each point. Its phase is invisible to a position measurement, which only sees $|\Psi|^2$. But the phase is not decoration: how the phase varies in space is momentum, and how it varies in time is energy. Two wavefunctions with the same $|\Psi|^2$ and different phase patterns move differently.

**A wave packet is a sum of plane waves (F3).** A plane wave $e^{ikx}$ has wavelength $2\pi/k$, so by de Broglie it has momentum $p = \hbar k$, and a definite one. Its $|\Psi|^2 = 1$ everywhere: perfectly known momentum, completely unknown position. To localize a particle you add plane waves with a spread of $k$ so that they cancel away from one region. F3 §2 told you the rule: narrow in $x$ needs wide in $k$, and the product of the widths has a floor. Multiply by $\hbar$ to convert $k$ to $p$. You should be able to guess the theorem now: $\sigma_x\sigma_p \ge \hbar\times(\text{F3's floor in } x\text{-}k\text{ units}) = \hbar/2$, with equality for a Gaussian.

**Guess the momentum operator.** You want an operation that, acting on $e^{ikx}$, returns $\hbar k$ times $e^{ikx}$. Differentiation does it: $\partial_x e^{ikx} = ik\,e^{ikx}$, so $-i\hbar\,\partial_x e^{ikx} = \hbar k\,e^{ikx}$. Guess: $\hat p = -i\hbar\,\partial_x$. For a packet, which is a superposition of many $k$'s, $\Psi^*\hat p\Psi$ integrated should average $\hbar k$ over the mixture. §3.3 derives this guess from the equation instead of assuming it.

**Guess the equation.** The same trick in time: $i\hbar\,\partial_t e^{-i\omega t} = \hbar\omega\,e^{-i\omega t} = E\,e^{-i\omega t}$ by Planck–Einstein $E = \hbar\omega$. So "energy" is $i\hbar\partial_t$ and "kinetic energy" is $\hat p^2/2m = -(\hbar^2/2m)\partial_x^2$. Demanding $E = p^2/2m + V$ as an operator statement gives
$$
i\hbar\,\partial_t\Psi = \Big(-\frac{\hbar^2}{2m}\partial_x^2 + V\Big)\Psi .
$$
The equation is first order in time (unlike the classical wave equation), complex, and linear. First order in time means $\Psi(x,0)$ alone determines the future, which is what a state should do.

**Expectation values are weighted averages.** If $|\Psi|^2$ is a probability density, the average position over many identically prepared particles is $\int x|\Psi|^2dx$. For momentum you cannot weight by $|\Psi|^2$, which knows nothing about $p$; the momentum information is in the phase, and $\hat p$ extracts it.

**The picture to keep.** Position space and momentum space are two descriptions of one state, related by a Fourier transform. A measurement of $x$ samples $|\Psi(x)|^2$; a measurement of $p$ samples $|\Phi(p)|^2$. The free Schrödinger equation leaves $|\Phi(p)|^2$ alone (nothing pushes the particle) and lets the different-$p$ components run apart in $x$ at different speeds $p/m$, so $|\Psi|^2$ spreads.

## 3. Formal development

### 3.1 The Schrödinger equation

**Postulate (Schrödinger equation).** For a particle of mass $m$ in a real potential $V(x)$,
$$
i\hbar\frac{\partial\Psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2\Psi}{\partial x^2} + V(x)\Psi \equiv \hat H\Psi. \tag{3.1}
$$

This is not derived. Here is the motivation, and its limits.

▶ **Try it first.** Substitute $\Psi = e^{i(kx-\omega t)}$ into (3.1) with $V$ constant. What relation between $\omega$ and $k$ does it impose, and what is it in terms of $E = \hbar\omega$, $p = \hbar k$?

**Motivation.**

1. $\partial_t e^{i(kx-\omega t)} = -i\omega\,e^{i(kx-\omega t)}$, so the left side is $\hbar\omega\,\Psi$.
2. $\partial_x^2 e^{i(kx-\omega t)} = -k^2e^{i(kx-\omega t)}$, so the right side is $(\hbar^2k^2/2m + V)\Psi$.
3. Equality for all $x, t$ requires
$$
\hbar\omega = \frac{\hbar^2k^2}{2m} + V, \qquad\text{i.e.}\qquad E = \frac{p^2}{2m} + V. \tag{3.2}
$$
4. Conversely: if you want an equation that is linear (so superpositions of solutions are solutions, which interference experiments demand), first order in $t$ (so $\Psi(x,0)$ fixes the future), and has plane waves obeying $E = p^2/2m + V$ when $V$ is constant, then (3.1) is the simplest choice. The second-order wave equation $\partial_t^2\Psi = c^2\partial_x^2\Psi$ gives $\omega = ck$, the photon's relation, not the electron's.

■ **What just happened.** The Schrödinger equation is the dispersion relation $E = p^2/2m + V$ with $E\to i\hbar\partial_t$ and $p\to -i\hbar\partial_x$. For non-constant $V$ we keep the same form; that is the leap, and experiment (the hydrogen spectrum, QM-6) justified it.

### 3.2 The statistical interpretation and conservation of probability

**Born rule.** $\rho(x,t) = |\Psi(x,t)|^2$ is the probability density for the particle's position:
$$
P(a\le x\le b,\ t) = \int_a^b|\Psi(x,t)|^2dx, \qquad \int_{-\infty}^\infty|\Psi(x,t)|^2dx = 1. \tag{3.3}
$$
Because (3.1) is linear, $c\Psi$ solves it whenever $\Psi$ does, so any square-integrable solution can be **normalized** by choosing $|c|$. Multiplying by a constant phase $e^{i\theta}$ changes nothing physical (Exercise A3). But normalization is imposed at one time. If the equation let $\int|\Psi|^2$ drift, the theory would be inconsistent.

**Theorem 3.1 (Probability is conserved locally).** Let $\Psi$ solve (3.1) with real $V$. Define the **probability current**
$$
J(x,t) = \frac{\hbar}{2mi}\Big(\Psi^*\frac{\partial\Psi}{\partial x} - \Psi\frac{\partial\Psi^*}{\partial x}\Big) = \frac{\hbar}{m}\operatorname{Im}\Big(\Psi^*\frac{\partial\Psi}{\partial x}\Big). \tag{3.4}
$$
Then
$$
\frac{\partial\rho}{\partial t} + \frac{\partial J}{\partial x} = 0, \tag{3.5}
$$
and if $\Psi\to0$ as $|x|\to\infty$ (so $J\to0$), $\frac{d}{dt}\int|\Psi|^2dx = 0$.

▶ **Try it first.** Write $\partial_t|\Psi|^2 = \Psi^*\partial_t\Psi + \Psi\,\partial_t\Psi^*$. You need $\partial_t\Psi^*$: take the complex conjugate of (3.1). Then look for a total $x$-derivative.

**Derivation.**

1. From (3.1): $\partial_t\Psi = \frac{i\hbar}{2m}\partial_x^2\Psi - \frac{i}{\hbar}V\Psi$.
2. Conjugate, using $V^* = V$: $\partial_t\Psi^* = -\frac{i\hbar}{2m}\partial_x^2\Psi^* + \frac{i}{\hbar}V\Psi^*$. This is the only place reality of $V$ enters.
3. $\partial_t\rho = \Psi^*\partial_t\Psi + \Psi\partial_t\Psi^* = \frac{i\hbar}{2m}\big(\Psi^*\partial_x^2\Psi - \Psi\partial_x^2\Psi^*\big) - \frac{i}{\hbar}V|\Psi|^2 + \frac{i}{\hbar}V|\Psi|^2.$ The potential terms cancel.
4. The bracket is a derivative: $\partial_x(\Psi^*\partial_x\Psi - \Psi\partial_x\Psi^*) = \partial_x\Psi^*\partial_x\Psi + \Psi^*\partial_x^2\Psi - \partial_x\Psi\partial_x\Psi^* - \Psi\partial_x^2\Psi^* = \Psi^*\partial_x^2\Psi - \Psi\partial_x^2\Psi^*$.
5. So $\partial_t\rho = \frac{i\hbar}{2m}\partial_x(\Psi^*\partial_x\Psi - \Psi\partial_x\Psi^*) = -\partial_x J$, since $\frac{i\hbar}{2m} = -\frac{\hbar}{2mi}$. That is (3.5).
6. Integrate (3.5) over $\mathbb R$: $\frac{d}{dt}\int\rho\,dx = -[J]_{-\infty}^{\infty} = 0$. ∎

■ **What just happened.** Probability flows like charge in EM-6's continuity equation $\partial_t\rho_e + \nabla\cdot\mathbf J_e = 0$: it can move but not appear or disappear. If $V$ had an imaginary part $-i\Gamma$, step 3 would leave $-\frac{2\Gamma}{\hbar}\rho$: probability would leak away. Physicists use exactly that as a model of absorption or decay, and it is why physical potentials are real. The project checks (3.5) on the grid to $5\times10^{-6}$.

For a plane wave $Ae^{i(kx-\omega t)}$, (3.4) gives $J = |A|^2\hbar k/m = \rho\,v$: density times velocity $p/m$. Example 4.2 does this and the standing wave.

### 3.3 Expectation values and the momentum operator, derived

The **expectation value** of position is the mean of the distribution (3.3):
$$
\langle x\rangle = \int_{-\infty}^\infty x\,|\Psi(x,t)|^2dx. \tag{3.6}
$$
It is the average over many identically prepared particles, each measured once. It is not the average of repeated measurements on one particle: the first measurement changes the state (QM-4).

**Theorem 3.2 (The momentum operator).** If $\Psi$ solves (3.1) and $\Psi$, $\partial_x\Psi$ vanish fast enough at infinity that boundary terms below are zero, then
$$
m\frac{d\langle x\rangle}{dt} = \int\Psi^*\Big(-i\hbar\frac{\partial}{\partial x}\Big)\Psi\,dx. \tag{3.7}
$$

▶ **Try it first.** $\frac{d}{dt}\langle x\rangle = \int x\,\partial_t\rho\,dx$. Use (3.5) and integrate by parts once. Then integrate by parts again on half of what you get.

**Derivation.**

1. $\frac{d\langle x\rangle}{dt} = \int x\,\partial_t\rho\,dx = -\int x\,\partial_xJ\,dx$ by (3.5).
2. By parts: $-\int x\,\partial_xJ\,dx = -[xJ]_{-\infty}^\infty + \int J\,dx = \int J\,dx$. The boundary term needs $xJ\to0$; for a normalizable state it does.
3. Insert (3.4): $\int J\,dx = \frac{\hbar}{2mi}\int(\Psi^*\partial_x\Psi - \Psi\partial_x\Psi^*)dx$.
4. By parts on the second term: $\int\Psi\,\partial_x\Psi^*dx = [\Psi\Psi^*]_{-\infty}^\infty - \int\Psi^*\partial_x\Psi\,dx = -\int\Psi^*\partial_x\Psi\,dx$.
5. So $\int J\,dx = \frac{\hbar}{2mi}\cdot2\int\Psi^*\partial_x\Psi\,dx = \frac{1}{m}\int\Psi^*(-i\hbar\partial_x)\Psi\,dx$. Multiply by $m$. ∎

■ **What just happened.** We never said what momentum is. We defined it as $m\,d\langle x\rangle/dt$, the classical definition applied to the mean, and the equation forced the form
$$
\langle p\rangle = \int\Psi^*\,\hat p\,\Psi\,dx, \qquad \hat p = -i\hbar\frac{\partial}{\partial x}. \tag{3.8}
$$
That is the guess from §2, now a consequence. The project's split-operator solver satisfies (3.7) to $10^{-13}$, because its kinetic step moves $\langle x\rangle$ by exactly $\langle p\rangle\,dt/m$.

**The general rule.** For any classical quantity $Q(x,p)$, the quantum expectation value is
$$
\langle Q\rangle = \int\Psi^*\,Q(\hat x, \hat p)\,\Psi\,dx, \qquad \hat x = x,\ \hat p = -i\hbar\,\partial_x. \tag{3.9}
$$
For example $\langle T\rangle = -\frac{\hbar^2}{2m}\int\Psi^*\partial_x^2\Psi\,dx$ and $\langle H\rangle = \int\Psi^*\hat H\Psi\,dx$. When $Q$ contains a product like $xp$, the order matters because $\hat x\hat p \ne \hat p\hat x$; the fix is in Exercise A5. Rule (3.9) is a postulate with a strong track record, not a theorem; Theorem 3.2 proves it only for $Q = p$, and §3.6 proves it for any $Q(p)$.

### 3.4 Hermitian operators and Ehrenfest's theorem

**Definition 3.3.** An operator $\hat Q$ on square-integrable functions is **Hermitian** if
$$
\langle f, \hat Qg\rangle = \langle\hat Qf, g\rangle \quad\text{for all (nice) } f, g. \tag{3.10}
$$
This is LA-8's $A^\dagger = A$ in function space, with $\langle\cdot,\cdot\rangle$ the $L^2$ inner product. "Nice" means the functions vanish at infinity and are differentiable enough; QM-4 is careful about domains, we are not.

**Theorem 3.4.** $\hat x$, $\hat p$ and $\hat H$ (with $V$ real) are Hermitian. Consequently every expectation value $\langle Q\rangle = \langle\Psi, \hat Q\Psi\rangle$ of a Hermitian $\hat Q$ is real.

▶ **Try it first.** For $\hat p$, write $\langle f, \hat pg\rangle = \int f^*(-i\hbar g')dx$ and move the derivative onto $f$ by parts. Watch the $i$.

**Derivation.**

1. $\hat x$: $\langle f, xg\rangle = \int f^*xg = \int(xf)^*g = \langle xf, g\rangle$ because $x$ is real.
2. $\hat p$: $\langle f, \hat pg\rangle = -i\hbar\int f^*g'\,dx = -i\hbar[f^*g]_{-\infty}^\infty + i\hbar\int(f')^*g\,dx = \int(-i\hbar f')^*g\,dx = \langle\hat pf, g\rangle$. The boundary term vanishes; $(-i\hbar)^* = i\hbar$ is what makes the sign come out. Without the $i$, $\partial_x$ alone is **anti**-Hermitian.
3. $\hat p^2$: $\langle f, \hat p\hat pg\rangle = \langle\hat pf, \hat pg\rangle = \langle\hat p\hat pf, g\rangle$, using step 2 twice. $V$ real is Hermitian by step 1's argument. Sums of Hermitian operators are Hermitian, so $\hat H$ is.
4. Reality: $\langle Q\rangle^* = \langle\Psi, \hat Q\Psi\rangle^* = \langle\hat Q\Psi, \Psi\rangle = \langle\Psi, \hat Q\Psi\rangle = \langle Q\rangle$, using $\langle u, v\rangle^* = \langle v, u\rangle$ and (3.10). ∎

■ **What just happened.** Observables must have real averages, and Hermiticity guarantees it (LA-8 Theorem 3.12 is the finite-dimensional version: Hermitian matrices have real eigenvalues and orthogonal eigenvectors). Every proof used integration by parts with a boundary term dropped. That drop is legitimate only for normalizable states; §8 returns to this.

**Time evolution of any expectation value.** Write (3.1) as $\partial_t\Psi = -\frac{i}{\hbar}\hat H\Psi$. For an operator $\hat Q$ with no explicit time dependence,

▶ **Try it first.** Differentiate $\langle\Psi, \hat Q\Psi\rangle$ with the product rule, substitute for $\partial_t\Psi$, and use Hermiticity of $\hat H$ to move it to the right. Remember that the inner product is conjugate-linear in the first slot.

1. $\frac{d}{dt}\langle Q\rangle = \langle\partial_t\Psi, \hat Q\Psi\rangle + \langle\Psi, \hat Q\partial_t\Psi\rangle$.
2. $\langle -\tfrac{i}{\hbar}\hat H\Psi, \hat Q\Psi\rangle = +\tfrac{i}{\hbar}\langle\hat H\Psi, \hat Q\Psi\rangle = \tfrac{i}{\hbar}\langle\Psi, \hat H\hat Q\Psi\rangle$ (conjugate-linearity, then (3.10)).
3. $\langle\Psi, \hat Q(-\tfrac{i}{\hbar}\hat H\Psi)\rangle = -\tfrac{i}{\hbar}\langle\Psi, \hat Q\hat H\Psi\rangle$.
4. Adding:
$$
\frac{d\langle Q\rangle}{dt} = \frac{i}{\hbar}\big\langle[\hat H, \hat Q]\big\rangle, \qquad [\hat A, \hat B] = \hat A\hat B - \hat B\hat A. \tag{3.11}
$$ ∎

**Theorem 3.5 (Ehrenfest).** For $\hat H = \hat p^2/2m + V(x)$,
$$
\frac{d\langle x\rangle}{dt} = \frac{\langle p\rangle}{m}, \qquad \frac{d\langle p\rangle}{dt} = \Big\langle -\frac{\partial V}{\partial x}\Big\rangle. \tag{3.12}
$$

**Derivation.** Compute two commutators acting on a test function $f$.

1. $[\hat H, \hat x]f = -\frac{\hbar^2}{2m}\big((xf)'' - xf''\big) = -\frac{\hbar^2}{2m}\cdot2f' = -\frac{i\hbar}{m}(-i\hbar f')$, so $[\hat H, \hat x] = -\frac{i\hbar}{m}\hat p$. ($V$ commutes with $x$.) Then (3.11): $\frac{d\langle x\rangle}{dt} = \frac{i}{\hbar}\cdot\big(-\frac{i\hbar}{m}\big)\langle p\rangle = \frac{\langle p\rangle}{m}$. This is Theorem 3.2 again by a second route.
2. $[\hat H, \hat p]f = V(-i\hbar f') - (-i\hbar)(Vf)' = -i\hbar Vf' + i\hbar V'f + i\hbar Vf' = i\hbar V'f$. ($\hat p^2$ commutes with $\hat p$.) Then $\frac{d\langle p\rangle}{dt} = \frac{i}{\hbar}\cdot i\hbar\langle V'\rangle = -\langle V'\rangle$. ∎

■ **What just happened.** Newton's law holds for averages, but with $\langle V'(x)\rangle$, not $V'(\langle x\rangle)$. They agree when $V'$ is linear in $x$ (free particle, uniform force, harmonic oscillator), so in those cases $\langle x\rangle(t)$ follows a classical trajectory exactly. For $V = \frac14x^4$ they differ, and the project shows $\langle x\rangle$ drifting off the classical orbit by up to $1.98$ (in units where the amplitude is $2$) within $t = 10$. Exercise B1 derives (3.12) the long way, by integration by parts; doing it once is worth the hour.

### 3.5 The uncertainty principle

The standard deviations are $\sigma_x^2 = \langle x^2\rangle - \langle x\rangle^2$ and $\sigma_p^2 = \langle p^2\rangle - \langle p\rangle^2$, with $\langle p^2\rangle = \langle\Psi, \hat p^2\Psi\rangle = \|\hat p\Psi\|^2 = \hbar^2\int|\partial_x\Psi|^2dx$ (Hermiticity, step 3 of Theorem 3.4).

**The de Broglie argument.** A packet of width $\Delta x$ containing waves near $k_0$ must be built from wavenumbers spread by $\Delta k$ such that the extreme components drift one full cycle out of phase across the packet: $\Delta k\,\Delta x \sim 2\pi$. Hence $\Delta p\,\Delta x \sim h$. This gets the scale and misses the constant; the theorem gets the constant.

**Definition 3.6 (Momentum-space wavefunction).**
$$
\Phi(p,t) = \frac{1}{\sqrt{2\pi\hbar}}\int_{-\infty}^\infty\Psi(x,t)\,e^{-ipx/\hbar}\,dx, \qquad \Psi(x,t) = \frac{1}{\sqrt{2\pi\hbar}}\int_{-\infty}^\infty\Phi(p,t)\,e^{ipx/\hbar}\,dp. \tag{3.13}
$$

**Conversion from F3.** With $\xi = p/(2\pi\hbar)$, F3's transform is $\hat\Psi(\xi) = \int\Psi e^{-2\pi ix\xi}dx = \int\Psi e^{-ipx/\hbar}dx$, so
$$
\Phi(p) = \frac{1}{\sqrt{2\pi\hbar}}\,\hat\Psi\Big(\frac{p}{2\pi\hbar}\Big), \qquad |\Phi(p)|^2dp = |\hat\Psi(\xi)|^2d\xi. \tag{3.14}
$$
The second identity follows from $dp = 2\pi\hbar\,d\xi$. It says the probability measures agree: momentum $p$ and spatial frequency $\xi$ are the same random variable in different units. The sanity check is F3's: $\Phi(0)\sqrt{2\pi\hbar} = \int\Psi\,dx$.

**Theorem 3.7 (Heisenberg–Kennard).** If $\|\Psi\| = 1$ and $x\Psi$, $\partial_x\Psi$ are square integrable, then
$$
\sigma_x\,\sigma_p \ge \frac{\hbar}{2}, \tag{3.15}
$$
with equality exactly for $\Psi = C\exp\!\big(-\frac{(x-x_0)^2}{4a^2} + \frac{ip_0x}{\hbar}\big)$, $a>0$: a Gaussian times a plane wave.

▶ **Try it first.** You need two facts: $|\Phi(p)|^2$ is the momentum distribution (so $\sigma_p$ is its standard deviation, §3.6), and F3 Theorem 3.11. Convert units.

**Derivation.**

1. *$\sigma_p$ is a spread of $|\Phi|^2$.* §3.6 (Theorem 3.8) shows $\langle p^n\rangle = \int p^n|\Phi(p)|^2dp$. So $\sigma_p$ is the standard deviation of the density $|\Phi|^2$.
2. *Change variables.* By (3.14), the density of $\xi = p/2\pi\hbar$ is $|\hat\Psi(\xi)|^2$, so $\sigma_p = 2\pi\hbar\,\sigma_\xi$ where $\sigma_\xi$ is the spread of $|\hat\Psi|^2$.
3. *F3.* F3 Theorem 3.11 (with the centering remark before it, which removes $\langle x\rangle$ and $\langle p\rangle$ by translation and modulation) gives $\sigma_x\sigma_\xi \ge 1/4\pi$.
4. *Multiply.* $\sigma_x\sigma_p = 2\pi\hbar\,\sigma_x\sigma_\xi \ge 2\pi\hbar/4\pi = \hbar/2$. Equality in F3 holds iff the centered function is $Ce^{-\pi a'x^2}$; undoing the translation (by $x_0$) and modulation (by $e^{ip_0x/\hbar}$) gives the stated family. ∎

■ **What just happened.** No new physics was used beyond "momentum is $\hbar$ times the Fourier variable of position," which is de Broglie plus (3.13). The uncertainty principle is a theorem about Fourier pairs. It is not about measurement disturbing the particle; it says no state *has* sharp $x$ and sharp $p$. QM-4 derives the general form $\sigma_A\sigma_B \ge \frac12|\langle[\hat A,\hat B]\rangle|$ from Cauchy–Schwarz with abstract operators; with $[\hat x, \hat p] = i\hbar$ it reproduces (3.15).

**The minimum-uncertainty Gaussian.** With $x_0 = p_0 = 0$ and normalization,
$$
\Psi(x) = (2\pi a^2)^{-1/4}e^{-x^2/4a^2}, \qquad \sigma_x = a, \qquad \sigma_p = \frac{\hbar}{2a}. \tag{3.16}
$$
Example 4.4 computes its $\Phi(p)$ and confirms $\sigma_x\sigma_p = \hbar/2$; the project gets $0.500000$. A Gaussian with a **chirp** (phase quadratic in $x$) is *not* minimal, even though its $|\Psi|^2$ is the same Gaussian (Exercise B5).

### 3.6 Momentum space

**Why $|\Phi|^2$ is the momentum distribution.** (3.13) writes $\Psi$ as a superposition of plane waves $e^{ipx/\hbar}/\sqrt{2\pi\hbar}$, each of definite momentum $p$, with amplitude $\Phi(p)$. By analogy with a finite superposition $\sum_nc_n\psi_n$ of eigenstates, where $|c_n|^2$ are the outcome probabilities (§3.9, and QM-4's measurement postulate), $|\Phi(p)|^2dp$ is the probability that a momentum measurement lands in $[p, p+dp]$. This is a motivation; the next theorem shows it is consistent with (3.8).

**Theorem 3.8 (Momentum space).** For nice $\Psi$:

(a) (Plancherel) $\displaystyle\int|\Phi(p)|^2dp = \int|\Psi(x)|^2dx$; more generally $\langle\Phi_1, \Phi_2\rangle = \langle\Psi_1, \Psi_2\rangle$. $\quad(3.17)$

(b) The transform of $\hat p\Psi$ is $p\,\Phi(p)$, so $\langle p^n\rangle = \int p^n|\Phi(p)|^2dp$. $\quad(3.18)$

(c) The transform of $x\Psi$ is $i\hbar\,\partial_p\Phi$, so in momentum space $\hat x = i\hbar\,\partial/\partial p$ and $\langle x\rangle = \int\Phi^*(i\hbar\partial_p)\Phi\,dp$. $\quad(3.19)$

▶ **Try it first.** For (b), transform $-i\hbar\Psi'$ and integrate by parts. For (c), differentiate (3.13) under the integral sign with respect to $p$.

**Derivation.**

1. (a) By (3.14) and F3 (3.15): $\int|\Phi(p)|^2dp = \int|\hat\Psi(\xi)|^2d\xi = \int|\Psi|^2dx$. The polarized version follows the same way from F3's $\langle f, g\rangle = \langle\hat f, \hat g\rangle$.
2. (b) $\frac{1}{\sqrt{2\pi\hbar}}\int(-i\hbar\Psi')e^{-ipx/\hbar}dx = \frac{-i\hbar}{\sqrt{2\pi\hbar}}\Big([\Psi e^{-ipx/\hbar}] - \int\Psi\cdot\big(-\tfrac{ip}{\hbar}\big)e^{-ipx/\hbar}dx\Big) = \frac{-i\hbar}{\sqrt{2\pi\hbar}}\cdot\frac{ip}{\hbar}\int\Psi e^{-ipx/\hbar}dx = p\,\Phi(p)$. Repeating, $\hat p^n\Psi \mapsto p^n\Phi$. Then $\langle p^n\rangle = \langle\Psi, \hat p^n\Psi\rangle = \langle\Phi, p^n\Phi\rangle = \int p^n|\Phi|^2dp$ by (a).
3. (c) $\partial_p\Phi = \frac{1}{\sqrt{2\pi\hbar}}\int\Psi\cdot\big(-\tfrac{ix}{\hbar}\big)e^{-ipx/\hbar}dx$, so $i\hbar\,\partial_p\Phi = \frac{1}{\sqrt{2\pi\hbar}}\int x\Psi\,e^{-ipx/\hbar}dx$, the transform of $x\Psi$. Then $\langle x\rangle = \langle\Psi, x\Psi\rangle = \langle\Phi, i\hbar\partial_p\Phi\rangle$ by (a). ∎

■ **What just happened.** The two representations are mirror images: in $x$-space, $\hat x$ multiplies and $\hat p = -i\hbar\partial_x$ differentiates; in $p$-space, $\hat p$ multiplies and $\hat x = +i\hbar\partial_p$ differentiates. The sign flip comes from the sign in the exponent. In LA-4 language, $\Psi(x)$ and $\Phi(p)$ are coordinates of **one vector** in two bases, position eigenfunctions and momentum eigenfunctions, and the Fourier transform is the change-of-basis matrix. Plancherel says it is unitary (LA-8 Theorem 3.20: unitary maps preserve inner products). A unitary change of basis that diagonalizes $\hat p$ is exactly what LA-8's spectral theorem promises for a Hermitian operator; QM-4 states this in Dirac notation as $\langle p|\Psi\rangle = \Phi(p)$, $\langle x|\Psi\rangle = \Psi(x)$, $\langle x|p\rangle = e^{ipx/\hbar}/\sqrt{2\pi\hbar}$. Momentum space will feel redundant until then. The project uses it as the reason the split-operator method works: $\hat p^2/2m$ is diagonal in the FFT basis.

### 3.7 The free particle and wave packets

With $V = 0$, try a separable plane wave.

▶ **Try it first.** Is $e^{i(kx - \hbar k^2t/2m)}$ a solution of (3.1) with $V = 0$? Can it be normalized?

1. By §3.1 step 3 with $V = 0$, $\Psi_k = e^{i(kx-\omega(k)t)}$ solves (3.1) iff
$$
\omega(k) = \frac{\hbar k^2}{2m}. \tag{3.20}
$$
2. $|\Psi_k|^2 = 1$ for all $x$, so $\int|\Psi_k|^2dx = \infty$. A plane wave is not a physical state; it is a building block.
3. Linearity: any superposition solves (3.1). The general free solution is
$$
\Psi(x,t) = \frac{1}{\sqrt{2\pi}}\int_{-\infty}^\infty\phi(k)\,e^{i(kx - \hbar k^2t/2m)}\,dk, \qquad \phi(k) = \frac{1}{\sqrt{2\pi}}\int\Psi(x,0)e^{-ikx}dx. \tag{3.21}
$$
At $t = 0$ this is Fourier inversion, so it matches any initial condition. In momentum language: $\Phi(p,t) = \Phi(p,0)\,e^{-ip^2t/2m\hbar}$, so **$|\Phi(p,t)|^2$ never changes** for a free particle. The project asserts this to $10^{-9}$.

**Phase and group velocity.** A single component's crests move at the **phase velocity**; the envelope of a packet concentrated near $k_0$ moves at the **group velocity** (EM-7 (3.29) derived $v_g = d\omega/dk$ for light in a medium).

▶ **Try it first.** Expand $\omega(k)$ to first order about $k_0$ in (3.21) and show the envelope is a rigid function of $x - \omega'(k_0)t$.

1. $v_p = \omega/k = \hbar k/2m$.
2. Expand $\omega(k) \approx \omega_0 + \omega_0'(k-k_0)$, $\omega_0' = d\omega/dk|_{k_0} = \hbar k_0/m$. Then
$\Psi \approx e^{i(k_0x - \omega_0t)}\frac{1}{\sqrt{2\pi}}\int\phi(k)e^{i(k-k_0)(x - \omega_0't)}dk = e^{i(k_0x-\omega_0t)}\,F(x - \omega_0't)$,
a carrier wave times an envelope $F$ that translates rigidly at $v_g = \hbar k_0/m$.
3. So
$$
v_g = \frac{d\omega}{dk} = \frac{\hbar k_0}{m} = \frac{p_0}{m}, \qquad v_p = \frac{\hbar k_0}{2m} = \frac{v_g}{2}. \tag{3.22}
$$

■ **What just happened.** The particle moves at the group velocity, which is the classical velocity $p/m$. The phase velocity is half of it and has no direct physical meaning (the overall phase of $\Psi$ is not observable; shifting $V$ by a constant changes $\omega$ and $v_p$ and nothing measurable). The dispersion $\omega \propto k^2$ is nonlinear, so the second-order term $\frac12\omega''(k-k_0)^2t$ dropped in step 2 eventually matters: the packet spreads. The project's figure `phase_group.png` shows crests sliding backward through the envelope; the measured $v_p$ at the packet center is $1.005$ for $k_0 = 2$ ($\hbar = m = 1$), the $0.5\%$ excess being the spreading correction.

### 3.8 The Gaussian wave packet, exactly

Take the minimum-uncertainty state (3.16), centered at $x = 0$ and given mean momentum $\hbar k_0$:
$$
\Psi(x,0) = (2\pi a^2)^{-1/4}\exp\Big(-\frac{x^2}{4a^2} + ik_0x\Big). \tag{3.23}
$$

▶ **Try it first.** Compute $\phi(k)$ from (3.21) (it is F3 (3.12) after completing the square), put it back into (3.21), and do the $k$-integral with $\int e^{-\alpha k^2 + bk}dk = \sqrt{\pi/\alpha}\,e^{b^2/4\alpha}$, valid for $\operatorname{Re}\alpha > 0$.

**Derivation.** Write $\beta = \hbar t/2ma^2$ (dimensionless time).

1. *Momentum amplitude.* $\phi(k) = \frac{(2\pi a^2)^{-1/4}}{\sqrt{2\pi}}\int e^{-x^2/4a^2 - i(k-k_0)x}dx$. The Gaussian integral with $\alpha = 1/4a^2$, $b = -i(k-k_0)$ gives $\sqrt{4\pi a^2}\,e^{-a^2(k-k_0)^2}$. Collecting constants:
$$
\phi(k) = \Big(\frac{2a^2}{\pi}\Big)^{1/4}e^{-a^2(k-k_0)^2}. \tag{3.24}
$$
2. *Exponent of the time integral.* In (3.21): $-a^2(k-k_0)^2 + ikx - i\frac{\hbar t}{2m}k^2 = -\alpha k^2 + bk - a^2k_0^2$ with $\alpha = a^2 + \frac{i\hbar t}{2m} = a^2(1+i\beta)$ and $b = 2a^2k_0 + ix$. $\operatorname{Re}\alpha = a^2 > 0$, so the formula applies.
3. *Integrate.* $\Psi(x,t) = \frac{1}{\sqrt{2\pi}}\big(\frac{2a^2}{\pi}\big)^{1/4}\sqrt{\frac{\pi}{\alpha}}\exp\big(\frac{b^2}{4\alpha} - a^2k_0^2\big)$. The prefactor simplifies to $(2\pi a^2)^{-1/4}(1+i\beta)^{-1/2}$.
4. *Complete the square (the error-prone step).* Put $v = \hbar k_0/m$. Then
$\frac{b^2}{4\alpha} - a^2k_0^2 = \frac{(2a^2k_0+ix)^2 - 4\alpha a^2k_0^2}{4\alpha} = \frac{-x^2 + 4ia^2k_0x - 2ia^2\hbar k_0^2t/m}{4\alpha}.$
We claim this equals $-\frac{(x-vt)^2}{4\alpha} + ik_0x - \frac{i\hbar k_0^2t}{2m}$. Multiply the claim by $4\alpha$ and expand: $-x^2 + 2xvt - v^2t^2 + 4ia^2k_0x - \frac{2ia^2\hbar k_0^2t}{m} + \frac{2i\hbar t}{m}ik_0x - \frac{2i\hbar t}{m}\cdot\frac{i\hbar k_0^2t}{2m}$. The terms $2xvt$ and $-\frac{2\hbar tk_0x}{m}$ cancel; $-v^2t^2$ and $+\frac{\hbar^2k_0^2t^2}{m^2}$ cancel. What remains is the left side. So
$$
\Psi(x,t) = \frac{(2\pi a^2)^{-1/4}}{\sqrt{1+i\beta}}\exp\Big(-\frac{(x - vt)^2}{4a^2(1+i\beta)} + ik_0x - \frac{i\hbar k_0^2t}{2m}\Big), \qquad \beta = \frac{\hbar t}{2ma^2}. \tag{3.25}
$$
5. *Density.* $|1+i\beta|^{-1} = (1+\beta^2)^{-1/2}$, and $\operatorname{Re}\frac{1}{1+i\beta} = \frac{1}{1+\beta^2}$. So
$$
|\Psi(x,t)|^2 = \frac{1}{\sqrt{2\pi}\,\sigma(t)}\exp\Big(-\frac{(x-vt)^2}{2\sigma(t)^2}\Big), \qquad \sigma(t) = a\sqrt{1 + \Big(\frac{\hbar t}{2ma^2}\Big)^2}. \tag{3.26}
$$ ∎

■ **What just happened.** The density stays Gaussian, its center moves at $v = \hbar k_0/m$ (group velocity), and its width grows. For $t \ll \tau$ the width is nearly constant; for $t\gg\tau$ it grows linearly, $\sigma \approx \hbar t/2ma = \sigma_p t/m$, which is what a cloud of classical particles with velocity spread $\sigma_p/m$ would do. The crossover time is the **spreading time**
$$
\tau = \frac{2ma^2}{\hbar}. \tag{3.27}
$$
A narrower packet spreads faster ($\tau\propto a^2$) because it contains a wider range of velocities. $\sigma_p = \hbar/2a$ is constant (§3.7), so $\sigma_x\sigma_p = \frac{\hbar}{2}\sqrt{1+\beta^2}$ grows: the evolved state is a chirped Gaussian, no longer minimal. SymPy in the project verifies that (3.25) solves (3.1); the split-operator solution matches (3.25) to $10^{-9}$ and $\sigma(20) = 10.04988$ for $a = 1$.

### 3.9 The time-independent Schrödinger equation and stationary states

▶ **Try it first.** Put $\Psi(x,t) = \psi(x)\varphi(t)$ into (3.1), divide by $\psi\varphi$, and argue that each side is a constant.

**Derivation.**

1. Substituting: $i\hbar\psi\varphi' = \varphi\big(-\frac{\hbar^2}{2m}\psi'' + V\psi\big)$. Divide by $\psi\varphi$: $i\hbar\frac{\varphi'}{\varphi} = \frac{1}{\psi}\big(-\frac{\hbar^2}{2m}\psi'' + V\psi\big)$.
2. The left side depends only on $t$, the right only on $x$. Both equal a constant; call it $E$.
3. Time part: $\varphi' = -\frac{iE}{\hbar}\varphi$, so $\varphi = e^{-iEt/\hbar}$ (constant absorbed into $\psi$).
4. Space part:
$$
-\frac{\hbar^2}{2m}\psi'' + V\psi = E\psi, \qquad\text{i.e.}\qquad \hat H\psi = E\psi. \tag{3.28}
$$
This is the **time-independent Schrödinger equation**: an eigenvalue problem for $\hat H$. ∎

**Theorem 3.9 (Stationary states).** Let $\Psi = \psi(x)e^{-iEt/\hbar}$ with $\hat H\psi = E\psi$, $\|\psi\| = 1$.

(a) $|\Psi|^2 = |\psi|^2$ and every $\langle Q(x,p)\rangle$ are independent of $t$.

(b) $E$ is real, $\langle H\rangle = E$ and $\sigma_H = 0$: the energy is definite.

(c) If $\{\psi_n\}$ are orthonormal eigenfunctions with energies $E_n$, then
$$
\Psi(x,t) = \sum_nc_n\psi_n(x)e^{-iE_nt/\hbar}, \qquad c_n = \langle\psi_n, \Psi(\cdot,0)\rangle \tag{3.29}
$$
solves (3.1) with $\Psi(x,0) = \sum c_n\psi_n$, and
$$
\sum_n|c_n|^2 = 1, \qquad \langle H\rangle = \sum_n|c_n|^2E_n. \tag{3.30}
$$

▶ **Try it first.** For (b), compute $\hat H^2\psi$. For (c), expand $\langle\Psi, \hat H\Psi\rangle$ as a double sum and use orthonormality.

**Derivation.**

1. (a) $|e^{-iEt/\hbar}| = 1$ for real $E$, and in $\int\Psi^*Q(x, -i\hbar\partial_x)\Psi\,dx$ the factors $e^{+iEt/\hbar}$ and $e^{-iEt/\hbar}$ cancel because $Q$ does not differentiate in $t$.
2. (b) $E = E\langle\psi,\psi\rangle = \langle\psi, \hat H\psi\rangle$, real by Theorem 3.4. So $\langle H\rangle = E$. $\hat H^2\psi = \hat H(E\psi) = E^2\psi$, so $\langle H^2\rangle = E^2$ and $\sigma_H^2 = E^2 - E^2 = 0$.
3. (c) Each term solves (3.1): $i\hbar\partial_t(\psi_ne^{-iE_nt/\hbar}) = E_n\psi_ne^{-iE_nt/\hbar} = \hat H(\psi_ne^{-iE_nt/\hbar})$. By linearity the sum does. At $t = 0$ it is $\sum c_n\psi_n$; taking $\langle\psi_m, \cdot\rangle$ and using $\langle\psi_m, \psi_n\rangle = \delta_{mn}$ gives $c_m = \langle\psi_m, \Psi(\cdot,0)\rangle$.
4. Norm: $\langle\Psi,\Psi\rangle = \sum_{m,n}c_m^*c_ne^{i(E_m-E_n)t/\hbar}\delta_{mn} = \sum|c_n|^2$.
5. Energy: $\hat H\Psi = \sum c_nE_n\psi_ne^{-iE_nt/\hbar}$, so $\langle\Psi, \hat H\Psi\rangle = \sum_{m,n}c_m^*c_nE_n\,e^{i(E_m-E_n)t/\hbar}\delta_{mn} = \sum|c_n|^2E_n$. ∎

■ **What just happened.** Stationary states are not static (the phase rotates) but nothing observable changes. All the dynamics of a general state is in the *relative* phases $e^{-i(E_m-E_n)t/\hbar}$ between components. The $|c_n|^2$ are constant in time and sum to 1; QM-4's measurement postulate says they are the probabilities of measuring $E_n$, and (3.30) is then just "expected value = $\sum$ probability × outcome."

**Two-state superposition.** Take real $\psi_1, \psi_2$ and $\Psi = \frac{1}{\sqrt2}\big(\psi_1e^{-iE_1t/\hbar} + \psi_2e^{-iE_2t/\hbar}\big)$. Then
$$
|\Psi|^2 = \tfrac12\psi_1^2 + \tfrac12\psi_2^2 + \psi_1\psi_2\cos(\omega_{21}t), \qquad \omega_{21} = \frac{E_2 - E_1}{\hbar}, \tag{3.31}
$$
because $\frac12\psi_1\psi_2\big(e^{i(E_2-E_1)t/\hbar} + e^{-i(E_2-E_1)t/\hbar}\big) = \psi_1\psi_2\cos\omega_{21}t$. The density sloshes at the **Bohr frequency** $\omega_{21}$. A charged particle whose charge distribution oscillates at $\omega_{21}$ radiates at that frequency (EM-8's dipole radiation), which is why spectral lines sit at $\Delta E/h$; QM-8 makes this a calculation. Example 4.5 puts numbers in.

### 3.10 Theorems about energy eigenfunctions in one dimension

Throughout, $V$ is real and $\psi$ is a normalizable solution of (3.28).

**Theorem 3.10.** (a) The eigenfunctions can be taken real. (b) $E > V_{\min} = \inf_xV(x)$. (c) Eigenfunctions with different energies are orthogonal. (d) (Stated) For the potentials in QM-3, the eigenfunctions, together with any continuum states, are complete: every normalizable $\Psi(x,0)$ is a sum (3.29) (plus an integral over the continuum).

▶ **Try it first.** For (a), conjugate (3.28). For (b), write $E = \langle\psi, \hat H\psi\rangle$ and look at the sign of each piece.

**Derivation.**

1. (a) Conjugating (3.28) with $V$, $E$ real: $\psi^*$ solves it with the same $E$. Linearity: $\psi_R = \frac12(\psi + \psi^*)$ and $\psi_I = \frac{1}{2i}(\psi - \psi^*)$ are real solutions, and $\psi = \psi_R + i\psi_I$. So every eigenfunction is a combination of real ones, and the eigenspace has a real basis. (In 1D bound states are nondegenerate, QM-3, so $\psi$ itself is real up to a constant phase.)
2. (b) $E = \langle\psi, \hat H\psi\rangle = \frac{1}{2m}\langle\hat p\psi, \hat p\psi\rangle + \int V|\psi|^2dx = \frac{\hbar^2}{2m}\int|\psi'|^2dx + \int V|\psi|^2dx$. The first term is $\ge0$; it is $>0$ because $\psi' \equiv 0$ would make $\psi$ constant, not normalizable. The second is $\ge V_{\min}\int|\psi|^2 = V_{\min}$. So $E > V_{\min}$.
3. (c) $E_m\langle\psi_n, \psi_m\rangle = \langle\psi_n, \hat H\psi_m\rangle = \langle\hat H\psi_n, \psi_m\rangle = E_n\langle\psi_n, \psi_m\rangle$ (Hermiticity, $E_n$ real). So $(E_m - E_n)\langle\psi_n, \psi_m\rangle = 0$. ∎
4. (d) is not proved here. (3.28) is a Sturm–Liouville problem (F2's Sturm–Liouville section, $p = 1$, $q = V$, weight $1$); on a finite interval with regular boundary conditions completeness is the Sturm–Liouville theorem, and on the whole line it is part of the spectral theorem for unbounded self-adjoint operators (Reed and Simon, vol. I, ch. VIII).

■ **What just happened.** Every 1D eigenfunction fact used either "$V$ is real" or "$\hat H$ is Hermitian." (b) is a classical statement in disguise: kinetic energy is positive, and confinement makes it strictly positive. It is why there is no state at the bottom of a well: zero-point energy. (c) is LA-8's spectral theorem in function space, the reason (3.29) has easy coefficients.

### 3.11 The Hamilton–Jacobi connection

This is the derivation CM-8 §3.10 promised. CM-8 wrote $\psi = Ae^{iS/\hbar}$; here $A = \sqrt\rho$ to make the probability meaning explicit. (Erwin Madelung, 1927, read the result as hydrodynamics; David Bohm, 1952, built a hidden-variable theory on it. We use it as mathematics.)

▶ **Try it first.** Put $\Psi = R\,e^{iS/\hbar}$ with $R = \sqrt\rho \ge 0$ and $S$ real into (3.1). Divide by $e^{iS/\hbar}$ and separate real and imaginary parts.

**Derivation** (subscripts are partial derivatives).

1. $\Psi_t = \big(R_t + \frac{i}{\hbar}RS_t\big)e^{iS/\hbar}$ and $\Psi_{xx} = \big(R_{xx} + \frac{2i}{\hbar}R_xS_x + \frac{i}{\hbar}RS_{xx} - \frac{1}{\hbar^2}RS_x^2\big)e^{iS/\hbar}$.
2. (3.1) divided by $e^{iS/\hbar}$: $i\hbar R_t - RS_t = -\frac{\hbar^2}{2m}R_{xx} - \frac{i\hbar}{2m}(2R_xS_x + RS_{xx}) + \frac{1}{2m}RS_x^2 + VR$.
3. Real part, divided by $-R$:
$$
\frac{\partial S}{\partial t} + \frac{1}{2m}\Big(\frac{\partial S}{\partial x}\Big)^2 + V + Q = 0, \qquad Q = -\frac{\hbar^2}{2m}\frac{\partial_x^2\sqrt\rho}{\sqrt\rho}. \tag{3.32}
$$
4. Imaginary part, times $2R/\hbar$: $2RR_t + \frac1m(2RR_xS_x + R^2S_{xx}) = 0$, i.e.
$$
\frac{\partial\rho}{\partial t} + \frac{\partial}{\partial x}\Big(\rho\,\frac{S_x}{m}\Big) = 0. \tag{3.33}
$$ ∎

■ **What just happened.** The Schrödinger equation is exactly two real equations. (3.33) is the continuity equation (3.5): indeed $J = \frac{\hbar}{m}\operatorname{Im}(\Psi^*\Psi_x) = \rho S_x/m$, so probability flows with velocity $S_x/m$, and $S_x$ plays the role of momentum, as $\partial S/\partial q = p$ did in CM-8. (3.32) is CM-8's Hamilton–Jacobi equation (3.28) for $H = p^2/2m + V$ with one extra term, the **quantum potential** $Q$. It is the only place $\hbar$ appears. When $\rho$ varies slowly on the scale of the wavelength $h/|S_x|$, $Q$ is negligible and the phase of $\Psi$ obeys classical Hamilton–Jacobi: that is the $\hbar\to0$ limit, the geometric optics of Schrödinger's 1926 argument, and the starting point of QM-8's WKB. $Q$ depends on the *shape* of $\rho$, not its size, which is why it does not fade for a dilute packet. For the free Gaussian, $Q$ pushes outward and causes the spreading of §3.8 (Exercise C2). The project computes $S = \hbar\arg\Psi$ and $\rho$ from a numerical free packet and finds (3.32) satisfied to relative error $1.5\times10^{-8}$; dropping $Q$ leaves a residual of $5\%$ of $S_t$.

### 3.12 Units and when quantum mechanics matters

$\hbar = 1.054572\times10^{-34}$ J s. For an electron, the natural length built from $\hbar$, $m_e$, $c$ is the Compton wavelength
$$
\lambda_C = \frac{h}{m_ec} = 2.4263\times10^{-12}\ \text{m}, \qquad \frac{\hbar}{m_ec} = 3.8616\times10^{-13}\ \text{m}. \tag{3.34}
$$
Localizing an electron below $\lambda_C$ requires momentum spread $\gtrsim m_ec$, i.e. relativistic energies that create pairs; the one-particle Schrödinger equation is meaningless there. Atoms are $10^{-10}$ m, far above $\lambda_C$, which is why non-relativistic QM works for chemistry.

When is a system "quantum"? Compare the de Broglie wavelength of a typical particle with the length scales of the problem. At temperature $T$, the **thermal de Broglie wavelength** is
$$
\lambda_{\rm th} = \frac{h}{\sqrt{2\pi mk_BT}}. \tag{3.35}
$$
SM-5's criterion: a gas is classical when $n\lambda_{\rm th}^3 \ll 1$, i.e. the interparticle spacing $n^{-1/3}$ is much larger than $\lambda_{\rm th}$. Numbers from the project: N$_2$ at 300 K has $\lambda_{\rm th} = 1.90\times10^{-11}$ m against an air spacing of $3.4\times10^{-9}$ m, so $n\lambda_{\rm th}^3 = 1.7\times10^{-7}$: classical. A conduction electron at 300 K has $\lambda_{\rm th} = 4.3\times10^{-9}$ m against a spacing of $2.3\times10^{-10}$ m in copper: deeply quantum, which is why metals need Fermi–Dirac statistics (SM-5).

## 4. Worked examples

### Example 4.1: A parabolic wavefunction

$\Psi(x,0) = A(a^2 - x^2)$ for $|x|\le a$, zero outside. Normalize, then find $\langle x\rangle$, $\langle x^2\rangle$, $\sigma_x$, $\langle p\rangle$, $\langle p^2\rangle$, $\sigma_p$, and check (3.15).

*Normalize.* $\int_{-a}^a(a^2-x^2)^2dx = 2\int_0^a(a^4 - 2a^2x^2 + x^4)dx = 2\big(a^5 - \frac23a^5 + \frac15a^5\big) = \frac{16a^5}{15}$. So $A = \sqrt{15/16a^5}$ (taken real and positive).

*Position.* $\langle x\rangle = 0$: the integrand $x(a^2-x^2)^2$ is odd. $\langle x^2\rangle = A^2\cdot2\int_0^ax^2(a^2-x^2)^2dx = A^2\cdot2\big(\frac{a^7}{3} - \frac{2a^7}{5} + \frac{a^7}{7}\big) = \frac{15}{16a^5}\cdot\frac{16a^7}{105} = \frac{a^2}{7}$. So $\sigma_x = a/\sqrt7 = 0.3780a$.

*Momentum.* $\langle p\rangle = -i\hbar A^2\int_{-a}^a(a^2-x^2)(-2x)dx = 0$ (odd integrand). It had to be: for a real $\Psi$, $\langle p\rangle = -i\hbar\int\Psi\Psi'dx = -\frac{i\hbar}{2}[\Psi^2] = 0$, and a real number that is $-i$ times something real must vanish. For $\langle p^2\rangle$ use $\hbar^2\int|\Psi'|^2dx$ (§3.5): $\Psi' = -2Ax$ inside, so $\langle p^2\rangle = 4A^2\hbar^2\int_{-a}^ax^2dx = 4\cdot\frac{15}{16a^5}\cdot\frac{2a^3}{3}\hbar^2 = \frac{5\hbar^2}{2a^2}$.

(Using $-\hbar^2\int\Psi\Psi''$ instead: $\Psi'' = -2A$ inside gives $2A^2\hbar^2\int(a^2-x^2) = 2\cdot\frac{15}{16a^5}\cdot\frac{4a^3}{3}\hbar^2 = \frac{5\hbar^2}{2a^2}$. Same answer. $\Psi'$ jumps at $x = \pm a$, so $\Psi''$ has delta functions there, but they multiply $\Psi(\pm a) = 0$.)

*Check.* $\sigma_p = \hbar\sqrt{5/2}/a$, so $\sigma_x\sigma_p = \hbar\sqrt{5/14} = 0.5976\hbar > 0.5\hbar$. The parabola is $20\%$ above the Gaussian floor. The project confirms all of this in SymPy.

### Example 4.2: Probability current of a plane wave and a standing wave

*Plane wave.* $\Psi = Ae^{i(kx-\omega t)}$: $\partial_x\Psi = ik\Psi$, so $\Psi^*\partial_x\Psi = ik|A|^2$ and $J = \frac{\hbar}{m}\operatorname{Im}(ik|A|^2) = |A|^2\frac{\hbar k}{m}$. With $\rho = |A|^2$, $J = \rho v$ where $v = \hbar k/m$ is the classical velocity. For an electron with $k = 10^{10}$ m$^{-1}$, $v = \hbar k/m_e = 1.158\times10^6$ m/s. A beam of $n$ electrons per m$^3$ described by this wave carries $nv$ electrons per m$^2$ per second.

*Standing wave.* $\Psi = A\cos(kx)e^{-i\omega t}$: $\Psi^*\partial_x\Psi = |A|^2\cos(kx)(-k\sin kx)$, which is real. Its imaginary part is zero, so $J = 0$. Equivalently $\Psi = \frac{A}{2}(e^{ikx} + e^{-ikx})e^{-i\omega t}$ carries equal currents both ways (Exercise B3 shows the cross terms drop out). Any $\Psi$ that is a real function times a phase depending only on $t$ has $J = 0$; every stationary bound state in QM-3 is of this type (Theorem 3.10(a)).

### Example 4.3: How fast does a localized particle spread?

Use (3.26)–(3.27) with $a = \sigma_x(0)$.

*Electron, $a = 10^{-10}$ m.* $\tau = 2m_ea^2/\hbar = 2(9.109\times10^{-31})(10^{-20})/(1.0546\times10^{-34}) = 1.728\times10^{-16}$ s. An electron localized to an atom-sized region, if freed, doubles its width in about $\sqrt3\,\tau\approx3\times10^{-16}$ s. After 1 s, $\sigma \approx \hbar t/2m_ea = 579$ km. That is not a paradox: it means a free electron localized to 1 Å has a velocity spread $\sigma_p/m_e = \hbar/2m_ea = 5.8\times10^5$ m/s.

*Dust grain, $m = 1\ \mu$g $= 10^{-9}$ kg.* Localized to $a = 1\ \mu$m: $\tau = 2(10^{-9})(10^{-12})/(1.0546\times10^{-34}) = 1.90\times10^{13}$ s $= 6.0\times10^5$ years. Even localized to 1 Å, $\tau = 1.90\times10^5$ s $= 2.2$ days, and in practice collisions with air molecules localize it again far sooner (QM-9's decoherence). Spreading is a quantum effect you never see in a grain because $\tau\propto m$.

### Example 4.4: The Gaussian in momentum space

For (3.16), $\Psi(x) = (2\pi a^2)^{-1/4}e^{-x^2/4a^2}$, find $\Phi(p)$ and $\sigma_x\sigma_p$.

By (3.13), $\Phi(p) = \frac{(2\pi a^2)^{-1/4}}{\sqrt{2\pi\hbar}}\int e^{-x^2/4a^2 - ipx/\hbar}dx$. The Gaussian integral with $\alpha = 1/4a^2$, $b = -ip/\hbar$ gives $\sqrt{4\pi a^2}\,e^{-a^2p^2/\hbar^2}$. So
$$
\Phi(p) = \Big(\frac{2a^2}{\pi\hbar^2}\Big)^{1/4}e^{-a^2p^2/\hbar^2}.
$$
Check the constant: $(2\pi a^2)^{-1/4}(4\pi a^2)^{1/2}(2\pi\hbar)^{-1/2} = 2^{-1/4+1-1/2}\pi^{-1/4+1/2-1/2}a^{1/2}\hbar^{-1/2} = (2a^2/\pi\hbar^2)^{1/4}$. Normalization: $\int|\Phi|^2dp = \sqrt{2a^2/\pi\hbar^2}\sqrt{\pi\hbar^2/2a^2} = 1$, as Plancherel demands. $|\Phi|^2 \propto e^{-2a^2p^2/\hbar^2} = e^{-p^2/2\sigma_p^2}$ with $\sigma_p = \hbar/2a$. So $\sigma_x\sigma_p = a\cdot\hbar/2a = \hbar/2$ exactly. The project computes $\Phi$ from the FFT, matches this formula to $10^{-10}$, and gets $\sigma_p = 0.500000$ for $a = \hbar = 1$. In $k$-units this is (3.24) with $k_0 = 0$, and $\Phi(p) = \hbar^{-1/2}\phi(p/\hbar)$.

### Example 4.5: Two states of an electron in a box

Take a particle in an infinite well $0 < x < L$. QM-3 derives $\psi_n = \sqrt{2/L}\sin(n\pi x/L)$, $E_n = n^2\pi^2\hbar^2/2mL^2$; check it now: $-\frac{\hbar^2}{2m}\psi_n'' = \frac{\hbar^2}{2m}\frac{n^2\pi^2}{L^2}\psi_n$, and $\psi_n(0) = \psi_n(L) = 0$. Let
$$
\Psi(x,t) = \tfrac{1}{\sqrt2}\big(\psi_1e^{-iE_1t/\hbar} + \psi_2e^{-iE_2t/\hbar}\big). \tag{4.5}
$$
*Density.* By (3.31), $|\Psi|^2 = \frac{1}{L}\big(\sin^2\frac{\pi x}{L} + \sin^2\frac{2\pi x}{L}\big) + \frac{2}{L}\sin\frac{\pi x}{L}\sin\frac{2\pi x}{L}\cos\omega_{21}t$, with $\omega_{21} = 3\pi^2\hbar/2mL^2$. At $t = 0$ the cross term is positive where $\sin(2\pi x/L) > 0$, the left half: the particle is on the left. Half a period later it is on the right.

*Energy.* $c_1 = c_2 = 1/\sqrt2$, so (3.30) gives $\langle H\rangle = \frac12(E_1 + E_2) = \frac52E_1$, and a measurement gives $E_1$ or $E_2$, each with probability $\frac12$; never $\frac52E_1$.

*Numbers for an electron in $L = 1$ nm.* $E_1 = \pi^2\hbar^2/2m_eL^2 = 6.02\times10^{-20}$ J $= 0.3760$ eV. $E_2 - E_1 = 3E_1 = 1.128$ eV. The slosh frequency is $f = (E_2 - E_1)/h = 2.728\times10^{14}$ Hz, period $3.67\times10^{-15}$ s, and light emitted at that frequency has $\lambda = c/f = 1.10\ \mu$m (near infrared). $\langle H\rangle = 0.940$ eV. In units $\hbar = m = L = 1$ the slosh period is $2\pi/\omega_{21} = 4/3\pi = 0.4244$ and $\langle x\rangle$ swings by $\pm16/9\pi^2 = \pm0.1801$ about $\frac12$ (Exercise B4).

### Example 4.6: Ehrenfest in a harmonic well

Units $\hbar = m = \omega = 1$, $V = \frac12x^2$. Start from the ground-state Gaussian displaced to $x_0 = 2$: $\Psi(x,0) = \pi^{-1/4}e^{-(x-2)^2/2}$ (this is (3.16) with $a = 1/\sqrt2$, so $\sigma_x = 0.7071$, $\sigma_p = 0.7071$).

*Ehrenfest.* (3.12) gives $\frac{d\langle x\rangle}{dt} = \langle p\rangle$ and $\frac{d\langle p\rangle}{dt} = -\langle x\rangle$, because $V' = x$ is linear, so $\langle V'\rangle = \langle x\rangle$. The initial values are $\langle x\rangle = 2$ and $\langle p\rangle = 0$ (real $\Psi$). So $\langle x\rangle = 2\cos t$, $\langle p\rangle = -2\sin t$: the classical orbit, exactly. At $t = \pi/2$, $\langle x\rangle = 0$ and $\langle p\rangle = -2$.

*Energy.* $\langle T\rangle = \frac12\langle p^2\rangle = \frac12\sigma_p^2 = \frac14$; $\langle V\rangle = \frac12\langle x^2\rangle = \frac12(x_0^2 + \sigma_x^2) = \frac12(4 + \frac12) = \frac94$. So $\langle H\rangle = \frac52$, conserved.

*What the project shows.* The split-operator run matches $2\cos t$ to $1.6\times10^{-5}$, holds $\langle H\rangle = 2.5$ to a relative $5\times10^{-6}$, and keeps $\sigma_x\sigma_p = 0.500000$ at every step: this displaced Gaussian (a coherent state, QM-3) stays minimum-uncertainty forever. Switch to $V = \frac14x^4$ from the same start: now $\langle -V'\rangle = -\langle x^3\rangle = -(x_0^3 + 3x_0\sigma_x^2) = -11$ at $t = 0$, while $-V'(\langle x\rangle) = -8$. The mean is pushed harder than a classical particle at the mean would be, and by $t = 10$ the two disagree by $1.98$.

## 5. Why it matters and where it connects

- **QM-3** solves (3.28) for the infinite well (Example 4.5's $\psi_n$), the oscillator (whose ground state is (3.16) with $a^2 = \hbar/2m\omega$), the delta and finite wells, and the step and barrier that the project previews. Its transmission and reflection coefficients are ratios of the currents (3.4) (Exercise B3).
- **QM-4** turns §3.6 into the general framework: $\Psi(x) = \langle x|\Psi\rangle$, $\Phi(p) = \langle p|\Psi\rangle$, Theorem 3.4 into "observables are self-adjoint," (3.11) into the Heisenberg equation, and Theorem 3.7 into the Robertson relation $\sigma_A\sigma_B \ge \frac12|\langle[\hat A,\hat B]\rangle|$.
- **QM-8 and QM-9.** WKB is (3.32) with $Q$ treated perturbatively. The Bohr frequency (3.31) becomes the resonance condition of Rabi oscillations and Fermi's golden rule. QM-9's Ehrenfest-based classical limit and decoherence explain why Example 4.3's dust grain never spreads in practice.
- **CM-8.** (3.32) is CM-8 (3.28) plus $Q$; $S_x = p$ in both. Poisson brackets $\{x, H\}$, $\{p, H\}$ are the classical versions of the commutators in Theorem 3.5.
- **EM-6, EM-7.** (3.5) is the continuity equation; (3.22) is EM-7's group velocity for the dispersion relation $\omega = \hbar k^2/2m$. EM-7 C3 is this module's spreading packet from the optics side.
- **F3 and QC-5.** The split-operator method alternates between position and momentum with the FFT; QC-5's quantum Fourier transform is the same change of basis on $2^n$ amplitudes, and QC-9's Trotterization is the splitting formula (7.1).
- **SM-5.** The criterion $n\lambda_{\rm th}^3\sim1$ of §3.12 decides when Bose–Einstein and Fermi–Dirac statistics replace Maxwell–Boltzmann.

## 6. Exercises

See [exercises.md](exercises.md). **A1–A5** are drills: normalizing an exponential, plane-wave velocities, why a global phase is invisible, the scales of §3.12 for liquid helium, and why $\hat x\hat p$ is not Hermitian. **B1–B8** are the module: B1 (Ehrenfest by integration by parts), B2 (the exponential in momentum space), B3 (current of counter-propagating waves), B4 ($\langle x\rangle$ and $\langle p\rangle$ in the box superposition), B5 (the chirped Gaussian), B6 (quartic Ehrenfest, code), B7 (the virial theorem), B8 (order of accuracy of the split-operator method, code). Essential: B1, B2, B4, B5. **C1–C3** go further: the uniform force in momentum space, the quantum potential as the cause of spreading, and the box ground state's uncertainty product.

## 7. Coding project

`project.py` (about 4 s) builds a split-operator solver and checks every number in §3–§4 and in the solutions. The scheme: over a step $dt$,
$$
e^{-i\hat H\,dt/\hbar} = e^{-iV\,dt/2\hbar}\,e^{-i\hat p^2dt/2m\hbar}\,e^{-iV\,dt/2\hbar} + O(dt^3). \tag{7.1}
$$
$V$ is diagonal on the position grid; $\hat p^2/2m$ is diagonal on the FFT grid (§3.6: $\hat p\mapsto p$); the FFT carries $\Psi$ between them. Each factor multiplies by phases of modulus 1, so the norm is conserved to rounding. The error is $O(dt^3)$ per step because the symmetric arrangement cancels the $O(dt^2)$ commutator term $[\hat T, V]dt^2/2\hbar^2$ of the Baker–Campbell–Hausdorff expansion; Exercise B8 measures the resulting second-order global error. The markdown cell in Part 3 repeats this.

The parts: (1) the scales of §3.12; (2) SymPy checks of Example 4.1, of (3.25) solving (3.1) and having width (3.26), and of the Madelung split (3.32)–(3.33); (3) the solver; (4) the free Gaussian: width against (3.26) to $10^{-8}$, $|\Phi|^2$ fixed, norm drift $2\times10^{-13}$, $\sigma_x\sigma_p$ growing as $\frac12\sqrt{1+(t/2)^2}$, and Example 4.4; (5) phase vs group velocity (`phase_group.png`, a grid of $\operatorname{Re}\Psi$ and $|\Psi|$ snapshots with the two speeds marked); (6) the current of Example 4.2 and the continuity equation on the grid; (7) Ehrenfest for harmonic and quartic potentials, energy conservation, and the coherent state's $\sigma_x\sigma_p = \hbar/2$ (`ehrenfest.png`); (8) a packet with $\langle E\rangle = 2.005 > V_0 = 1.5$ hitting a step: $13.2\%$ reflects, against $13.1\%$ from averaging the plane-wave $R(k) = ((k-q)/(k+q))^2$ over $|\phi(k)|^2$, and $R(k_0) = 1/9$ (QM-3 derives $R$); (9) the Hamilton–Jacobi check of §3.11; (10) the SI numbers of Examples 4.3 and 4.5; (11) the exercise numbers.

**Build it yourself first:** write `evolve(grid, psi, V, dt, n)` from (7.1) using `np.fft`, and test it on the free Gaussian against (3.25) before reading the file's version. Then do the "Your turn" items: a uniform force (Ehrenfest says $\langle x\rangle = x_0 + Ft^2/2m$ with the free width), the box slosh period measured with the solver (high walls), and the free equation solved directly in momentum space.

## 8. Where you'll struggle

- **The equation is not derived, and you will want it to be.** Every attempt to "derive" (3.1) assumes something equivalent to it. §3.1 is a motivation. Accept it as a postulate that has passed every non-relativistic test, and spend your energy on what follows from it.
- **The sign $e^{-iEt/\hbar}$ and the Fourier conventions will cost you an hour.** F1 used $e^{+i\omega t}$, F3 used $e^{-2\pi ix\xi}$ with no prefactor, this module uses $(2\pi\hbar)^{-1/2}e^{-ipx/\hbar}$. Fix: before any Fourier computation, write (3.14) on the page and check $\int|\Phi|^2 = 1$ for your answer. If it comes out $2\pi$ or $1/\hbar$, a prefactor is wrong.
- **Integration by parts drops boundary terms, and only normalizable states earn that.** Theorems 3.1, 3.2, 3.4 all threw away $[\ldots]_{-\infty}^\infty$. For a plane wave the boundary terms do not vanish, and "$\hat p$ is Hermitian" fails as stated. When you see a derivation produce nonsense, look first at a dropped boundary term.
- **"The wavefunction is the particle" is a misreading.** $|\Psi|^2$ spreading to 579 km (Example 4.3) does not mean the electron is 579 km wide. Each detection finds it at one point; $\Psi$ encodes the distribution of outcomes.
- **The Gaussian-packet algebra is long, and completing the square is where errors live.** Step 4 of §3.8 is where everyone loses a factor of 2 or an $i$. Fix: do it with $k_0 = 0$ first (much shorter), then get $k_0 \ne 0$ by noting the $k_0$ packet is the $k_0 = 0$ one boosted (multiply by $e^{ik_0x - i\hbar k_0^2t/2m}$ and shift $x\to x - vt$). Check your result by putting it into (3.1), as the project does.
- **Momentum space feels redundant.** Until QM-4 it partly is. Its payoffs here: (3.18) makes $\langle p^n\rangle$ an ordinary average, the free particle is trivial in it, and the split-operator method lives on it.
- **The split-operator time step must resolve the fastest phase.** The kinetic factor $e^{-ip^2dt/2m\hbar}$ is exact per mode, but the splitting error grows with $p_{\max}^2dt/m\hbar$ and $V_{\max}dt/\hbar$. If $V_{\max}dt/\hbar$ approaches $\pi$, the potential phase aliases and the result is garbage while the norm stays perfect. A conserved norm is not evidence of accuracy (LA-13's stability-vs-accuracy distinction). Check energy and halve $dt$.

## 9. Check yourself

1. Write the Schrödinger equation and say which part is the dispersion relation $E = p^2/2m + V$.
2. Where in the proof of $\frac{d}{dt}\int|\Psi|^2 = 0$ is $V$ real used?
3. Write the probability current and evaluate it for $e^{ikx}$ and $\cos kx$.
4. How is $\hat p = -i\hbar\partial_x$ obtained from the definition $\langle p\rangle = m\,d\langle x\rangle/dt$?
5. State Ehrenfest's theorem. For which potentials does $\langle x\rangle(t)$ follow a classical trajectory exactly?
6. Why must an observable's operator be Hermitian? Show $\hat p$ is.
7. Convert F3's $\sigma_x\sigma_\xi \ge 1/4\pi$ into $\sigma_x\sigma_p\ge\hbar/2$.
8. What is $\hat x$ in momentum space, and where does its sign come from?
9. Give $v_p$ and $v_g$ for a free electron with wavenumber $k$. Which one is the particle's velocity?
10. Write $\sigma(t)$ for a free Gaussian packet and the spreading time. Why does a narrower packet spread faster?
11. Prove that a stationary state has $\sigma_H = 0$, and that $E > V_{\min}$.
12. Put $\Psi = \sqrt\rho\,e^{iS/\hbar}$ into the Schrödinger equation: what are the two real equations, and which term disappears as $\hbar\to0$?

## 10. Pair with

- **Griffiths & Schroeter, *Introduction to Quantum Mechanics* (3rd ed.), ch. 1 and §2.1, §2.4.** The standard route through §3.2–3.5 and §3.9 with many more short problems; its §1.5–1.6 derivation of $\hat p$ is the one used here.
- **MIT 8.04 (Zwiebach), lectures 1–6, and Zwiebach, *Mastering Quantum Mechanics*, ch. 3–5.** More careful than Griffiths about momentum space, wave packets, and the free-particle Gaussian; good on why the equation must be complex and first order.
- **Shankar, *Principles of Quantum Mechanics*, ch. 4 (optional now).** The postulates stated in Hilbert-space form; read it after QM-4 if this module's "two bases" remark left you wanting the formal version.
- **Feynman Lectures on Physics, vol. III, ch. 16.** Feynman's own heuristic route to the equation and to the momentum amplitude, with a physicist's view of why $|\Phi(p)|^2$ is the momentum distribution.
