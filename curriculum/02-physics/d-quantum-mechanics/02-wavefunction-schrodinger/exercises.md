# QM-2 Exercises

Conventions as in the notes: $\hat p = -i\hbar\,\partial_x$, $\langle f, g\rangle = \int\bar f g\,dx$, time dependence $e^{-iEt/\hbar}$, $\Phi(p) = (2\pi\hbar)^{-1/2}\int\Psi e^{-ipx/\hbar}dx$, $J = \frac{\hbar}{m}\operatorname{Im}(\Psi^*\partial_x\Psi)$. Where a problem says "units $\hbar = m = 1$" it matches `project.py`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** $\Psi(x,t) = Ae^{-\lambda|x|}e^{-i\omega t}$ with $\lambda > 0$.
(a) Normalize it.
(b) Find $\langle x\rangle$, $\langle x^2\rangle$ and $\sigma_x$.
(c) What is the probability of finding the particle more than one standard deviation from $\langle x\rangle$? Give a number.

**A2.** (a) Show that $e^{i(kx-\omega t)}$ solves the free Schrödinger equation if and only if $\omega = \hbar k^2/2m$.
(b) For an electron with de Broglie wavelength $1$ Å, compute $k$, $v_g$, $v_p$ and the kinetic energy in eV.

**A3.** (a) Show that $\Psi$ and $e^{i\theta}\Psi$ ($\theta$ a real constant) give the same $|\Psi|^2$ and the same $\langle Q\rangle$ for every $Q(x,p)$.
(b) Now let $\theta = \theta(x)$. Show that $|\Psi|^2$ is unchanged but $\langle p\rangle$ shifts by $\hbar\langle\theta'(x)\rangle$. What does $\theta = k_0x$ do to a state?

**A4.** (a) Compute $\hbar/m_ec$ and explain in one sentence why the Schrödinger equation is unreliable for an electron localized more tightly than this.
(b) Liquid helium-4 at $4.2$ K has density $125$ kg/m³. Compute $\lambda_{\rm th}$ (3.35), the interatomic spacing $n^{-1/3}$, and $n\lambda_{\rm th}^3$. Is the liquid classical?

**A5.** (a) Show that $\hat x\hat p$ is not Hermitian by computing its adjoint: show $(\hat x\hat p)^\dagger = \hat p\hat x = \hat x\hat p - i\hbar$.
(b) For any real normalizable $\Psi$, compute $\langle\Psi, \hat x\hat p\Psi\rangle$. Why is the answer a red flag for an "observable"?
(c) Show that $\frac12(\hat x\hat p + \hat p\hat x)$ is Hermitian.

## B. Core (do all)

**B1.** (Ehrenfest the long way.) Without using commutators, prove $\frac{d\langle p\rangle}{dt} = \langle -\partial_xV\rangle$: differentiate $\langle p\rangle = -i\hbar\int\Psi^*\partial_x\Psi\,dx$ under the integral, substitute the Schrödinger equation and its conjugate, and integrate by parts. Say which boundary terms you drop and why.

**B2.** (The two-sided exponential in momentum space.) For $\Psi(x) = \sqrt\lambda\,e^{-\lambda|x|}$:
(a) Compute $\Phi(p)$ and check $\int|\Phi|^2dp = 1$.
(b) Compute $\langle p^2\rangle$ two ways: $\hbar^2\int|\Psi'|^2dx$ and $\int p^2|\Phi|^2dp$.
(c) Find $\sigma_x\sigma_p$ and compare with F3 Exercise B2's $\sigma_x\sigma_\xi = \sqrt2/4\pi$ through (3.14).
(d) $\Psi$ has a kink at $0$. Why is $\langle p^2\rangle$ still finite, while $\langle p^4\rangle$ is infinite? Use the decay of $\Phi$.

**B3.** (Current of counter-propagating waves.) Let $\Psi = Ae^{ikx} + Be^{-ikx}$ with $k > 0$.
(a) Show that $J = \frac{\hbar k}{m}(|A|^2 - |B|^2)$: the interference terms in $|\Psi|^2$ carry no current.
(b) A wave $Ae^{ikx} + Be^{-ikx}$ for $x < 0$ meets a step and becomes $Ce^{iqx}$ for $x > 0$. Using (3.5) for a stationary state, explain why the reflection and transmission probabilities are $R = |B|^2/|A|^2$ and $T = \frac{q}{k}\frac{|C|^2}{|A|^2}$, and why $R + T = 1$.

**B4.** (The box superposition.) Units $\hbar = m = L = 1$ unless stated. For Example 4.5's state (4.5):
(a) Compute $x_{12} = \int_0^L\psi_1x\psi_2\,dx$ and show $\langle x\rangle(t) = \frac L2 - \frac{16L}{9\pi^2}\cos\omega_{21}t$.
(b) Compute $\langle p\rangle(t)$ two ways: as $m\,d\langle x\rangle/dt$, and directly from $p_{12} = \int\psi_1(-i\hbar\partial_x)\psi_2\,dx$. Show they agree.
(c) What are the slosh period and the amplitude of $\langle x\rangle$ in units $\hbar = m = L = 1$?

**B5.** (The chirped Gaussian.) Let $\Psi(x,0) = (2\pi a^2)^{-1/4}\exp\big(-\frac{x^2}{4a^2} + icx^2\big)$, $c$ real.
(a) Show $\sigma_x = a$ and $\sigma_p^2 = \hbar^2\big(\frac{1}{4a^2} + 4c^2a^2\big)$, so $\sigma_x\sigma_p = \frac\hbar2\sqrt{1 + 16c^2a^4}$.
(b) Show that the free packet (3.25) at time $t$ (with $k_0 = 0$) is a chirped Gaussian, and find its width and chirp in terms of $\beta = \hbar t/2ma^2$. Check that (a) reproduces $\sigma_x\sigma_p = \frac\hbar2\sqrt{1+\beta^2}$.
(c) Units $\hbar = m = 1$, $a = 1$, $c = -\frac14$. Show this state is a minimum-uncertainty Gaussian evolved for a *negative* time, so it first narrows. Find when its width is smallest, and the minimum width.

**B6.** (Ehrenfest in a quartic well; code.) $V = \frac14x^4$, units $\hbar = m = 1$, $\Psi(x,0)$ a real Gaussian of width $a$ centered at $x_0 = 2$.
(a) Show that $\langle -V'\rangle = -(x_0^3 + 3x_0a^2)$ at $t = 0$, versus $-V'(\langle x\rangle) = -x_0^3$.
(b) Estimate, from (a), how long $\langle x\rangle$ stays within $0.1$ of the classical trajectory starting at rest from $x_0$, as a function of $a$.
(c) Extend Part 7 of `project.py` to measure that time for $a = 0.3$, $1/\sqrt2$ and $1$, and compare with (b). Does a narrower packet stay classical longer indefinitely? What limits it?

**B7.** (Virial theorem.) (a) Use (3.11) to show $\frac{d}{dt}\langle xp\rangle = 2\langle T\rangle - \langle x\,V'(x)\rangle$, where $\langle xp\rangle$ means the expectation of $\hat x\hat p$.
(b) Conclude that in a stationary state $2\langle T\rangle = \langle xV'\rangle$. For $V = \alpha|x|^n$ show $2\langle T\rangle = n\langle V\rangle$; for the oscillator, $\langle T\rangle = \langle V\rangle$.
(c) Evaluate $\frac{d}{dt}\langle xp\rangle$ at $t = 0$ for Example 4.6's state and check it against the classical $\frac{d}{dt}(xp)$ along $x = 2\cos t$.

**B8.** (Order of the split-operator method; code.) (a) Using $e^{A}e^{B} = e^{A + B + \frac12[A,B] + \cdots}$, show that the unsymmetric splitting $e^{-iV\,dt/\hbar}e^{-iT\,dt/\hbar}$ has a local error of order $dt^2$, and that the symmetric (7.1) cancels it.
(b) In `project.py`, evolve the displaced Gaussian of Example 4.6 in $V = \frac14x^4$ to $t = 2$ with $dt = 0.01, 0.005, 0.0025$, compare each with a reference run at $dt = 0.000625$, and report the ratios of successive errors. What order does this show?
(c) For a grid with $dx = 0.047$, what is $p_{\max}$? Explain why a perfectly conserved norm does not tell you $dt$ is small enough.

## C. Challenge (do at least one)

**C1.** (Uniform force in momentum space.) For $V = -Fx$ the momentum-space equation is $i\hbar\partial_t\Phi = \frac{p^2}{2m}\Phi - Fi\hbar\,\partial_p\Phi$ (use (3.19)).
(a) Solve it by characteristics: show
$\Phi(p,t) = \Phi_0(p - Ft)\exp\!\Big(-\frac{i}{6m\hbar F}\big[p^3 - (p-Ft)^3\big]\Big)$.
(b) Deduce that $|\Phi|^2$ translates rigidly, $\langle p\rangle = p_0 + Ft$ and $\sigma_p$ is constant. Is $\sigma_x(t)$ the same as for a free particle?
(c) Check (a) against the split-operator solution for a Gaussian with $a = 1$, $F = 0.1$, $t = 6$ ($\hbar = m = 1$). Why is the split-operator method exact here up to a global phase?

**C2.** (The quantum potential spreads the packet.) For the free Gaussian (3.26), $\sqrt\rho\propto\exp(-(x-vt)^2/4\sigma^2)$.
(a) Compute the quantum potential $Q$ (3.32) and the force $-\partial_xQ$. Which way does it push?
(b) The velocity field is $S_x/m$, with $J = \rho S_x/m$. Show that the curves $x(t) = vt + y_0\,\sigma(t)/a$ satisfy $\dot x = S_x/m$ for each $y_0$, and that $m\ddot x = -\partial_xQ$ along them.
(c) Explain in two sentences why the free packet spreads even though $V = 0$, in this language.

**C3.** (Uncertainty in a box.) (a) For the box ground state $\psi_1 = \sqrt{2/L}\sin(\pi x/L)$, compute $\sigma_x$, $\sigma_p$ and $\sigma_x\sigma_p$.
(b) Can any state confined to $0 < x < L$ reach $\sigma_x\sigma_p = \hbar/2$? Prove your answer using Theorem 3.7's equality case.
