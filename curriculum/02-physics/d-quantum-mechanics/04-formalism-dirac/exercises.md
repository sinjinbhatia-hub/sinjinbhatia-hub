# QM-4 Exercises

Conventions as in the notes: $\langle\phi|\psi\rangle$ is conjugate-linear in the first slot, $\langle x|p\rangle = (2\pi\hbar)^{-1/2}e^{ipx/\hbar}$, $\hat U(t) = e^{-i\hat Ht/\hbar}$, Schrödinger picture unless stated. Where a problem says "units $\hbar = m = \omega = 1$" it matches `project.py`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** In $\mathbb C^2$ let $|\phi\rangle = \frac{1}{\sqrt2}(i, 1)^\top$ and $|\psi\rangle = \frac{1}{\sqrt5}(1, 2i)^\top$.
(a) Compute $\langle\phi|\psi\rangle$ and $\langle\psi|\phi\rangle$, and check they are complex conjugates.
(b) Write the matrix of $\hat M = |\psi\rangle\langle\phi|$. Compute $\hat M|\phi\rangle$ and $\operatorname{Tr}\hat M$, and explain why $\operatorname{Tr}(|\psi\rangle\langle\phi|) = \langle\phi|\psi\rangle$ in general.
(c) Show $\hat M^\dagger = |\phi\rangle\langle\psi|$ by computing both matrices.

**A2.** (a) From the definition (3.9), prove $(\hat A\hat B)^\dagger = \hat B^\dagger\hat A^\dagger$ and $\langle\phi|\hat A|\psi\rangle^* = \langle\psi|\hat A^\dagger|\phi\rangle$.
(b) For Hermitian $\hat A$, $\hat B$, which of $\hat A\hat B$, $[\hat A, \hat B]$, $i[\hat A, \hat B]$, $\hat A\hat B + \hat B\hat A$ are Hermitian? Prove each answer.
(c) Show that $\langle\psi|\hat A^\dagger\hat A|\psi\rangle\ge0$ for every $\hat A$ and $|\psi\rangle$.

**A3.** (Projectors.) Let $\{|n\rangle\}$ be an orthonormal basis of $\mathbb C^N$ and $\hat P_n = |n\rangle\langle n|$.
(a) Show $\hat P_n\hat P_m = \delta_{nm}\hat P_n$ and $\hat P_n^\dagger = \hat P_n$.
(b) Show that a Hermitian $\hat P$ with $\hat P^2 = \hat P$ has eigenvalues only $0$ and $1$.
(c) Use (3.3) to show $\sum_n\operatorname{Tr}\hat P_n = N$, and use (3.10) to write $\operatorname{Tr}\hat Q = \sum_nq_n$.

**A4.** (Commutators against Poisson brackets.) Using only $[\hat x, \hat p] = i\hbar$ and (3.14):
(a) Compute $[\hat x, \hat p^2]$ and $[\hat x^2, \hat p]$.
(b) Prove by induction $[\hat x, \hat p^n] = i\hbar n\hat p^{n-1}$, and conclude $[\hat x, f(\hat p)] = i\hbar f'(\hat p)$ for a polynomial $f$.
(c) Compute the Poisson bracket $\{x, p^n\}$ (CM-8 (3.12)) and check the dictionary (3.15).

**A5.** (A degenerate measurement.) Let $\hat Q = \begin{pmatrix}2&i&0\\-i&2&0\\0&0&3\end{pmatrix}$ and $|\psi\rangle = \frac{1}{\sqrt3}(1, 1, 1)^\top$.
(a) Find the eigenvalues and an orthonormal eigenbasis. Which eigenvalue is degenerate?
(b) Find the probability of each outcome and the post-measurement state for each, using the Lüders rule (3.12).
(c) Compute $\langle Q\rangle$ two ways and $\sigma_Q$.
(d) Show that a different orthonormal basis of the degenerate eigenspace gives the same probabilities and the same post-measurement state.

## B. Core (do all)

**B1.** (Compatibility, degenerate case.) Let $\hat A = \operatorname{diag}(1, 1, 2)$ and $\hat B = \begin{pmatrix}0&1&0\\1&0&0\\0&0&5\end{pmatrix}$.
(a) Show $[\hat A, \hat B] = 0$.
(b) Show that the standard basis is an eigenbasis of $\hat A$ but not of $\hat B$. Which step of Theorem 3.6(b) fails?
(c) Construct a common orthonormal eigenbasis following the proof of Theorem 3.6(c), and list the joint eigenvalues $(a, b)$. Is $\{\hat A, \hat B\}$ a complete set of commuting observables on $\mathbb C^3$?

**B2.** (The Schrödinger–Robertson inequality.)
(a) Keep the real part in step 3 of the proof of Theorem 3.7 to show
$$\sigma_A^2\sigma_B^2\ \ge\ \Big(\tfrac12\langle\{\hat A, \hat B\}\rangle - \langle A\rangle\langle B\rangle\Big)^2 + \Big(\tfrac{1}{2i}\langle[\hat A, \hat B]\rangle\Big)^2,$$
where $\{\hat A, \hat B\} = \hat A\hat B + \hat B\hat A$. The first term is the squared covariance.
(b) For the chirped Gaussian $\psi(x) = (2\pi a^2)^{-1/4}\exp(-x^2/4a^2 + icx^2)$ ($c$ real), QM-2 Exercise B5 found $\sigma_x = a$ and $\sigma_x\sigma_p = \frac\hbar2\sqrt{1 + 16c^2a^4}$. Compute the covariance $\frac12\langle\{\hat x, \hat p\}\rangle - \langle x\rangle\langle p\rangle$ and show the chirped Gaussian saturates the inequality of (a) for every $c$, though it violates equality in (3.18) for $c\ne0$.
(c) Evaluate both sides for $a = 1$, $c = -\frac14$ (units $\hbar = 1$).

**B3.** (Saturating states.) Let $\hat A$, $\hat B$ be Hermitian and $\lambda$ real.
(a) Suppose $(\hat A - i\lambda\hat B)|\psi\rangle = \mu|\psi\rangle$ with $|\psi\rangle$ normalized. Show $\langle A\rangle = \operatorname{Re}\mu$ and $\langle B\rangle = -\operatorname{Im}\mu/\lambda$, and hence that $(\hat A - \langle A\rangle)|\psi\rangle = i\lambda(\hat B - \langle B\rangle)|\psi\rangle$.
(b) Conclude that $|\psi\rangle$ gives equality in (3.17) (and in B2(a)).
(c) For $\hat A = \hat x$, $\hat B = \hat p$, $\lambda = -1/m\omega$, show $\hat A - i\lambda\hat B$ is proportional to the oscillator's lowering operator $\hat a$. Why would $\lambda = +1/m\omega$ give nothing? What are the minimum-uncertainty states with $\sigma_x = \sqrt{\hbar/2m\omega}$?

**B4.** (Spreading from the Heisenberg picture.) For the free particle, $\hat x_H(t) = \hat x + \hat pt/m$ (§3.11).
(a) Show $\sigma_x^2(t) = \sigma_x^2(0) + \frac{2t}{m}\operatorname{cov}(x, p) + \frac{t^2}{m^2}\sigma_p^2$, with $\operatorname{cov}(x, p) = \frac12\langle\{\hat x, \hat p\}\rangle - \langle x\rangle\langle p\rangle$ at $t = 0$.
(b) For the minimum-uncertainty Gaussian of width $a$, recover QM-2's $\sigma_x(t) = a\sqrt{1 + (\hbar t/2ma^2)^2}$. Give $\sigma_x(2)$ for $\hbar = m = a = 1$.
(c) Compute $[\hat x_H(t), \hat x_H(0)]$ and use (3.17) to show $\sigma_x(t)\sigma_x(0)\ge\hbar t/2m$. Interpret: why can't a free particle's position be sharp at two different times?

**B5.** (Neutrino oscillations.) Two flavor states are related to two mass states by $|\nu_e\rangle = \cos\theta|1\rangle + \sin\theta|2\rangle$, $|\nu_\mu\rangle = -\sin\theta|1\rangle + \cos\theta|2\rangle$, with $\hat H|i\rangle = E_i|i\rangle$ and $E_i\approx pc + m_i^2c^4/2pc$.
(a) Write $\hat H$ in the flavor basis in the form (3.28) and identify $\mathbf h$.
(b) Starting in $|\nu_\mu\rangle$, derive $P_{\nu_\mu\to\nu_e}(t) = \sin^22\theta\sin^2(\Delta m^2c^4t/4\hbar E)$, with $\Delta m^2 = m_2^2 - m_1^2$ and $E\approx pc$, both directly from (3.23) and from (3.30).
(c) With $t = L/c$, show the phase is $1.267\,\Delta m^2[\text{eV}^2]\,L[\text{km}]/E[\text{GeV}]$, and find the first maximum $L$ for $\Delta m^2 = 2.5\times10^{-3}$ eV$^2$, $E = 1$ GeV.

**B6.** (Translations, parity, Bloch.)
(a) Using $[\hat x, \hat p] = i\hbar$ and (3.14), show $[\hat x, \hat T(a)] = a\hat T(a)$ for $\hat T(a) = e^{-i\hat pa/\hbar}$ (expand the exponential), and hence $\hat T(a)^\dagger\hat x\hat T(a) = \hat x + a$. Conclude that $\hat T(a)|x\rangle = |x + a\rangle$.
(b) Show $\hat T(a)$ commutes with $\hat H = \hat p^2/2m + V(\hat x)$ if and only if $V(x + a) = V(x)$ for all $x$ (use $\hat T^\dagger V(\hat x)\hat T = V(\hat x + a)$, from (a) and (3.16)).
(c) For parity $\hat\Pi$, show $\hat\Pi^2 = \hat 1$, that $\hat\Pi$ is Hermitian and unitary, and that $\hat\Pi\hat p\hat\Pi = -\hat p$. Then show that if $V$ is even, every nondegenerate eigenstate of $\hat H$ has definite parity.
(d) Explain why Theorem 3.15 does not give a degeneracy for the double well of `project.py` Part 8 but does for the ring.

**B7.** (Energy–time bound in a two-state system; code.) Let $\hat H = \varepsilon\sigma_z + \Delta\sigma_x$, start in $|{\uparrow}\rangle$, units $\hbar = 1$.
(a) Show $\sigma_H = |\Delta|$ in this state and that it is constant in time.
(b) Using (3.21), show $|d\langle\sigma_z\rangle/dt|\le2\sigma_H\sigma_{\sigma_z}$. In `project.py`, check this along the trajectory for $\varepsilon = 0.3$, $\Delta = 0.4$ on $t\in[0, 10]$.
(c) The Mandelstam–Tamm "speed limit" says no state reaches an orthogonal state before $t_\perp = \pi\hbar/2\sigma_H$. Show that for $\varepsilon = 0$ the state reaches $|{\downarrow}\rangle$ exactly at $t_\perp$, and give $t_\perp$ for $\Delta = 0.4$. Why does $\varepsilon\ne0$ never reach $|{\downarrow}\rangle$?

**B8.** ($\hat x$ in the momentum representation.)
(a) Using (3.6), show $\langle p|\hat x|\psi\rangle = i\hbar\,\partial_p\Phi(p)$. (Write $x\langle p|x\rangle = i\hbar\,\partial_p\langle p|x\rangle$ and insert (3.4).)
(b) Hence write the time-independent Schrödinger equation for $V = -Fx$ (uniform force) in momentum space and solve it.
(c) Show directly in the momentum representation that $[\hat x, \hat p] = i\hbar$.

## C. Challenge (do at least one)

**C1.** (Interaction picture and the rotating-wave approximation; code.) A spin in a static field along $z$ plus a transverse oscillating field: $\hat H = \frac{\hbar\omega_0}{2}\sigma_z + \hbar\Omega\cos(\omega t)\,\sigma_x$, $\Omega\ll\omega_0$.
(a) With $\hat H_0 = \frac{\hbar\omega_0}{2}\sigma_z$, compute $\hat V_I(t)$ from (3.27). Use $e^{i\theta\sigma_z/2}\sigma_\pm e^{-i\theta\sigma_z/2} = e^{\pm i\theta}\sigma_\pm$ with $\sigma_\pm = (\sigma_x\pm i\sigma_y)/2$.
(b) Drop the terms oscillating at $\omega_0 + \omega$ (the rotating-wave approximation). Show the remaining $\hat V_I$ is time independent in a frame rotating at the detuning $\delta = \omega_0 - \omega$, and that the flip probability is $\frac{\Omega^2}{\Omega^2 + \delta^2}\sin^2\big(\frac12\sqrt{\Omega^2 + \delta^2}\,t\big)$. Match it to (3.30).
(c) In `project.py`, integrate the full equation for $\omega_0 = \omega = 1$, $\Omega = 0.05$ from $|{\downarrow}\rangle$ and report the flip probability at $t = \pi/\Omega$. Why is it not exactly 1?

**C2.** (Why the canonical commutator needs infinite dimensions.)
(a) Prove that no $N\times N$ matrices satisfy $XP - PX = i\hbar I$.
(b) For the truncated ladder matrices of Example 4.2, prove $[a, a^\dagger] = I - N|N - 1\rangle\langle N - 1|$ and hence $[X, P] = i\hbar(I - N|N - 1\rangle\langle N - 1|)$.
(c) Prove that if $[\hat x, \hat p] = i\hbar$ held for bounded operators, then $[\hat x, \hat p^n] = i\hbar n\hat p^{n-1}$ would give $n\hbar\|\hat p^{n-1}\|\le2\|\hat x\|\|\hat p\|\|\hat p^{n-1}\|$ for all $n$, a contradiction. (So at least one of $\hat x$, $\hat p$ is unbounded: the domain issues of §3.1 are unavoidable.)

**C3.** (A mixed state.) A source emits $|{\uparrow}\rangle$ or $|{+x}\rangle = (|{\uparrow}\rangle + |{\downarrow}\rangle)/\sqrt2$ with equal probability.
(a) Write $\hat\rho$ (3.32) as a $2\times2$ matrix and check $\operatorname{Tr}\hat\rho = 1$.
(b) Write $\hat\rho = \frac12(I + \mathbf r\cdot\boldsymbol\sigma)$ and find $\mathbf r$. Show $\operatorname{Tr}\hat\rho^2 = \frac12(1 + |\mathbf r|^2)$ and evaluate it.
(c) Show no ket $|\chi\rangle$ has $|\chi\rangle\langle\chi| = \hat\rho$, and find a *different* ensemble of pure states with the same $\hat\rho$. Why does this mean no measurement can tell the two sources apart?
