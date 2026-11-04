# QM-3 Exercises

Conventions as in the notes: infinite well on $[0, a]$; finite well on $-a < x < a$ with $z_0 = a\sqrt{2mV_0}/\hbar$; $a_\pm = (2\hbar m\omega)^{-1/2}(\mp i\hat p + m\omega\hat x)$; $T$ and $R$ defined by currents (3.18). Where a problem says "$\hbar = m = 1$", lengths and energies are in those units.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** An electron is confined to an infinite well of width $a = 1$ nm. Compute $E_1$ in eV, the energy of the $2\to1$ transition, and the wavelength of the emitted photon. Which part of the spectrum is it in?

**A2.** Verify $\langle\psi_1|\psi_3\rangle = 0$ and $\langle\psi_2|\psi_2\rangle = 1$ for the infinite well by doing the integrals directly (not by quoting (3.4)).

**A3.** Apply $a_+$ to $\psi_0$ from (3.11) as a differential operator and show the result is $\psi_1 = (m\omega/\pi\hbar)^{1/4}\sqrt{2m\omega/\hbar}\,x\,e^{-m\omega x^2/2\hbar}$, normalized. Then apply $a_-$ to $\psi_1$ and confirm you get $\psi_0$ back with coefficient $\sqrt1$.

**A4.** A particle hits a step of height $V_0$ with $E = \frac43V_0$. Compute $k'/k$, $R$ and $T$, check $R + T = 1$, and compute the naive $|F/A|^2$.

**A5.** For the delta-well bound state (3.20), compute the probability of finding the particle at $|x| > 1/\kappa$. Sketch $\psi$ and mark where $\psi'$ is discontinuous.

## B. Core (do all)

**B1.** (Ladder matrix elements.) (a) Using (3.13) and (3.14), show $\langle\psi_m|\hat x|\psi_n\rangle = \sqrt{\hbar/2m\omega}\,\big(\sqrt n\,\delta_{m,n-1} + \sqrt{n+1}\,\delta_{m,n+1}\big)$. (b) This "selection rule" says a dipole ($\propto x$) can only change $n$ by $\pm1$. Explain in one sentence why a diatomic molecule's vibrational spectrum (in the oscillator approximation) has a single line. (c) Compute $\langle\psi_0|\hat x^4|\psi_0\rangle$ by ladders, and check against the Gaussian integral.

**B2.** (A discontinuous initial state.) In the infinite well, $\Psi(x, 0) = \sqrt{2/a}$ for $0 < x < a/2$ and $0$ for $a/2 < x < a$. (a) Check normalization. (b) Find $c_n$. (c) Compute $P(E_1)$ and $P(E_2)$. (d) Show that $\sum|c_n|^2E_n$ diverges, and explain physically why this initial state has infinite energy while every finite truncation of its expansion has finite energy.

**B3.** (Sloshing oscillator.) $\Psi(x, 0) = (\psi_0 + \psi_1)/\sqrt2$. Find $\langle x\rangle(t)$ and $\langle p\rangle(t)$ using ladders, and verify Ehrenfest's theorem $d\langle p\rangle/dt = -m\omega^2\langle x\rangle$ (QM-2 (3.12)). What is $\langle H\rangle$?

**B4.** (Where nondegeneracy fails.) (a) A particle on a ring of circumference $L$ ($\psi(x + L) = \psi(x)$, $V = 0$) has eigenstates $e^{\pm 2\pi inx/L}$. They are degenerate. Find the step of the proof of Theorem 3.2 that fails. (b) Same question for the free-particle scattering states $e^{\pm ikx}$ on the line. (c) Explain why Theorem 3.4's step 3 needs nondegeneracy, and show with $e^{\pm ikx}$ (for $V = 0$, which is even) that an eigenfunction of an even potential need not have definite parity.

**B5.** (A real finite well.) An electron sits in a finite well of width $2a = 1$ nm and depth $V_0 = 10$ eV. (a) Compute $z_0$ and the number of bound states. (b) Estimate $E_1$ from the infinite-well approximation and say whether the true value is higher or lower. (c) For a very shallow well ($z_0\ll1$) show from (3.22) that $\kappa a \approx z_0^2$, so the binding energy is $|E| \approx z_0^4\hbar^2/2ma^2$.

**B6.** (Above the barrier.) (a) Starting from $M_{22}$ in §3.8 and the substitution $\kappa\to il$, $l = \sqrt{2m(E - V_0)}/\hbar$, derive $T = [1 + V_0^2\sin^2(la)/4E(E - V_0)]^{-1}$ for $E > V_0$. (b) At which energies is $T = 1$? (c) Show that $T\to1$ as $E\to\infty$ and find $T$ in the limit $E\to V_0^+$; compare with the limit $E\to V_0^-$ of (3.26).

**B7.** (Current conservation.) (a) For any real $V(x)$ that is constant (possibly different) on the far left and far right, prove $R + T = 1$ from the $x$-independence of $J$ in a stationary state. (b) Verify $R + T = 1$ explicitly for the delta potential (3.21). (c) Show that the transfer matrix of any region has $\det = k_0/k_N$ and explain how this encodes (a).

**B8.** (Code: the step in `transfer()`.) (a) Use `transfer(E, [0], [0, V0])` from `project.py` to compute $R$ and $T$ for a step over $E/V_0 \in (1, 5]$, and plot them against (3.23). (b) Plot the naive $|F/A|^2$ on the same axes and find where it exceeds 1. (c) Compute $R$ for a step **down**, $V = 0 \to -3E$. Explain why a particle reflects off a drop.

**B9.** (Parity.) (a) Show that under $x\to-x$, $a_+ \to -a_+$ (as an operator acting on functions of $x$), and conclude that $\psi_n$ has parity $(-1)^n$. (b) For any even $V$, prove $\langle x\rangle = 0$ in every bound eigenstate. (c) Is $\langle x\rangle = 0$ for the superposition in B3? Why does that not contradict (b)?

## C. Challenge (do at least one)

**C1.** (The two-delta double well.) $V = -\alpha[\delta(x - b) + \delta(x + b)]$, units $\hbar = m = 1$. (a) Show the even bound state satisfies $\kappa(1 + \tanh\kappa b) = 2\alpha$ and the odd one $\kappa(1 + \coth\kappa b) = 2\alpha$. (b) Show the odd state exists only if $\alpha b > \frac12$. (c) For $\alpha b\gg1$, show the splitting is $\Delta \approx 2\alpha^2e^{-2\alpha b} = 4|E_b|e^{-\kappa_0d}$, with $d = 2b$, $\kappa_0 = \alpha$, $E_b$ the single-delta energy. (d) For $\alpha = 1$ and $b = 1, 2, 3$, compare the exact $\Delta$ (solve numerically) with (c).

**C2.** (Code: the triangular well.) $V(x) = Fx$ for $x > 0$, $V = \infty$ for $x < 0$ (a particle on a floor under gravity, or an electron in a MOSFET inversion layer). (a) Show that with $x = \ell(y + \epsilon)$, $\ell = (\hbar^2/2mF)^{1/3}$, the equation becomes Airy's equation $u'' = yu$, and that the boundary conditions give $E_n = -a_n(\hbar^2F^2/2m)^{1/3}$ with $a_n$ the zeros of $\operatorname{Ai}$. (b) With $\hbar = m = F = 1$, compute $E_1, E_2, E_3$ using `scipy.special.ai_zeros` and check them with `fd_eigs`.

**C3.** (The Gamow integral.) (a) Evaluate $\gamma = \frac1\hbar\int_{r_1}^{r_2}\sqrt{2m(\frac{2Zk_ee^2}{r} - E)}\,dr$ exactly, with $r_2 = 2Zk_ee^2/E$, by substituting $r = r_2\sin^2u$, to get (3.30). (b) Expand for $r_1\ll r_2$ to get $\gamma \approx \frac{\pi}{2}\frac{\sqrt{2mE}}{\hbar}r_2 - \frac{2\sqrt{2mE}}{\hbar}\sqrt{r_1r_2}$, and show the first term is $\propto Z/\sqrt E$. (c) Use this to explain, in two sentences, why a factor of 2 in $E$ gives tens of powers of 10 in $t_{1/2}$.
