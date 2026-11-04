# QM-3: One-Dimensional Problems: Infinite Well, Harmonic Oscillator, Free Particle, Delta and Finite Wells, Scattering and Tunneling

**Time:** ~16 hours · **Prerequisites:** QM-2 · **Unlocks:** QM-4, SM-5 (box states; the oscillator), CM-6 (quantized modes)

> **Conventions in this module.** $\hbar$ is explicit. One particle of mass $m$ on the line; $\hat H = -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x)$ and we solve the time-independent equation $\hat H\psi = E\psi$ (QM-2 (3.28)). Inner product $\langle f, g\rangle = \int\bar f g\,dx$, conjugate on the first slot; $\langle\psi_m|\psi_n\rangle$ means the same thing. The **infinite well** sits on $[0, a]$ with $\psi_n = \sqrt{2/a}\sin(n\pi x/a)$, $n = 1, 2, \dots$ (Griffiths). The **finite well** is $V = -V_0$ on $-a < x < a$ (width $2a$, Griffiths), with $z_0 = a\sqrt{2mV_0}/\hbar$. **Oscillator** ladder operators are $a_\pm = (2\hbar m\omega)^{-1/2}(\mp i\hat p + m\omega\hat x)$, with $[a_-, a_+] = 1$ (Griffiths; Shankar writes $a = a_-$, $a^\dagger = a_+$). States $\psi_n$ start at $n = 0$ for the oscillator and $n = 1$ for the box. Wavenumbers: $k = \sqrt{2mE}/\hbar$ where $E > V$, $\kappa = \sqrt{2m(V_0 - E)}/\hbar$ where $E < V$. Scattering: incident amplitude $A$ from the left, reflected $B$, transmitted $F$; $R = |B/A|^2$ and $T = (k_{\rm right}/k_{\rm left})|F/A|^2$, which is $|F/A|^2$ only when the speeds on the two sides are equal (§3.4). For the Dirac comb (§3.12), $K$ is the Bloch wavenumber and $k$ is still $\sqrt{2mE}/\hbar$. `project.py` uses $\hbar = m = 1$ except where SI units are stated.

## 0. Where this sits

You need QM-2: the Schrödinger equation, the current $J$ (QM-2 (3.4)), stationary states and the expansion $\Psi = \sum c_n\psi_ne^{-iE_nt/\hbar}$ (QM-2 (3.29)), and the 1D eigenfunction theorems (QM-2 Theorem 3.10). This module solves the exactly solvable 1D problems and builds three tools: the **ladder algebra** (QM-5's angular momentum, QC's photons), **matching** organized as a **transfer matrix** (LA-4's composition), and the **finite-difference Hamiltonian** (symmetric tridiagonal, LA-8). SM-5 needs the box and the oscillator; CM-6's modes become quantized oscillators.

## 1. The problem that forced this idea

By 1925 the Bohr–Sommerfeld rule (QM-1 §3.3) produced discrete energies for a few periodic orbits but had no mechanism, and it got the oscillator wrong: $E_n = n\hbar\omega$, while molecular band spectra (Mulliken's 1924 isotope shifts) pointed to half-integers. Schrödinger's first 1926 paper solved hydrogen; his second solved the oscillator and found $E_n = (n + \frac12)\hbar\omega$ with no extra rule. Discreteness came from one requirement: the wavefunction must be normalizable. A confined wave can only take certain shapes.

The particle in a box is not a toy. Sommerfeld's 1927–28 theory of metals treated conduction electrons as fermions in a box and fixed the heat-capacity puzzle of classical electron theory (SM-5). Semiconductor quantum dots are boxes a few nanometres across whose color is set by $E_n \propto 1/a^2$; the $\pi$ electrons of a conjugated dye behave roughly like particles in a box the length of the chain, which is why longer chains absorb redder light.

The second problem was alpha decay. Uranium-238 emits alpha particles of about 4.2 MeV, yet Rutherford's scattering experiments showed the Coulomb barrier around the nucleus is higher than that. Classically the alpha cannot get out. Worse, half-lives range from microseconds to billions of years while energies change by only a factor of about two (Geiger and Nuttall's empirical correlation, 1911). In 1928 Gamow, and independently Gurney and Condon, explained both with barrier penetration: in a region where $E < V$ the wave decays exponentially but does not vanish, and the leaked fraction depends exponentially on $\sqrt{V - E}$. It was the first quantum calculation of a decay rate. Tunneling then became engineering: Fowler and Nordheim's field emission (1928), Esaki's tunnel diode (1957), and Binnig and Rohrer's scanning tunneling microscope (IBM Zürich, 1981), whose current changes by about a factor of ten per ångström of gap.

## 2. Intuition

**A confined wave must fit.** A string fixed at both ends vibrates only at frequencies whose half-wavelengths fit an integer number of times into its length (F2, CM-6). Here wavelength is momentum ($\lambda = h/p$) and frequency is energy. Fit $n$ half-wavelengths into a box of length $a$: $\lambda_n = 2a/n$, $p_n = nh/2a$, $E_n = n^2h^2/8ma^2$. That is the whole infinite-well answer. Squeezing a particle costs energy $\propto 1/a^2$: the uncertainty principle as a spectrum.

**Curvature is set by $E - V$.** Write the equation as $\psi'' = \frac{2m}{\hbar^2}(V - E)\psi$. Where $E > V$ (classically allowed), $\psi''$ has the opposite sign to $\psi$: it curves back toward the axis and oscillates, faster where $E - V$ is larger. Where $E < V$ (forbidden), $\psi''$ has the same sign as $\psi$: it grows or decays exponentially at rate $\kappa = \sqrt{2m(V - E)}/\hbar$. A wavefunction can leak into a forbidden region; it just decays there.

Now guess the module.

- **Bound states.** If $E < V$ on both sides, a solution decaying on the left will generically blow up on the right. Only special $E$ give decay on both sides: quantization is a shooting problem, and each higher state turns over once more inside the well, so the $n$-th state has $n - 1$ nodes.
- **The oscillator.** Classical oscillator frequency does not depend on amplitude, so guess evenly spaced levels. The lowest cannot be at the bottom: a Gaussian of width $\sigma$ costs kinetic energy $\hbar^2/8m\sigma^2$ and potential energy $\frac12m\omega^2\sigma^2$; minimizing gives $\sigma^2 = \hbar/2m\omega$ and $E = \frac12\hbar\omega$. Guess $E_n = (n + \frac12)\hbar\omega$. Evenly spaced levels suggest operators that move up or down one rung; factoring $p^2 + m^2\omega^2x^2 = (m\omega x - ip)(m\omega x + ip)$, exact for numbers and off by a commutator for operators, builds them.
- **Tunneling.** Inside a barrier of width $a$ the wave decays like $e^{-\kappa x}$, so the probability reaching the far side is down by about $e^{-2\kappa a}$. That is the thick-barrier law (3.27) up to a prefactor of order 1.
- **Reflection without a wall.** Waves reflect from any abrupt change of medium (a light rope tied to a heavy one, light hitting glass). Quantum particles reflect from a step they classically clear.
- **Two wells.** Two identical wells coupled by tunneling behave like two coupled pendulums (CM-6): the common level splits into an even and an odd combination, by an amount set by the tunneling amplitude $\sim e^{-\kappa b}$.

## 3. Formal development

### 3.1 The infinite square well

$V = 0$ for $0 < x < a$ and $V = \infty$ elsewhere. Outside, $\psi = 0$ (any nonzero $\psi$ would make $\int V|\psi|^2 = \infty$). Inside, $-\frac{\hbar^2}{2m}\psi'' = E\psi$:
$$\psi'' = -k^2\psi, \qquad k = \frac{\sqrt{2mE}}{\hbar}. \tag{3.1}$$
($E > 0$ by QM-2 Theorem 3.10(b), since $V_{\min} = 0$.) We require $\psi$ continuous, so $\psi(0) = \psi(a) = 0$. We do **not** require $\psi'$ continuous at the walls: an infinite jump in $V$ allows a kink (the delta function of §3.5 shows how a finite kink comes from an infinite spike).

▶ **Try it first.** Solve (3.1) with $\psi(0) = \psi(a) = 0$ and normalize.

1. General solution: $\psi = A\sin kx + B\cos kx$.
2. $\psi(0) = B = 0$.
3. $\psi(a) = A\sin ka = 0$. $A = 0$ gives nothing, so $\sin ka = 0$: $k = n\pi/a$. $n = 0$ gives $\psi \equiv 0$; negative $n$ gives $-\psi_{|n|}$, the same state. So $n = 1, 2, 3, \dots$
4. Energies: $E_n = \hbar^2k_n^2/2m$,
$$k_n = \frac{n\pi}{a}, \qquad E_n = \frac{n^2\pi^2\hbar^2}{2ma^2}. \tag{3.2}$$
5. Normalize: $\int_0^a|A|^2\sin^2(n\pi x/a)\,dx = |A|^2a/2 = 1$, so $A = \sqrt{2/a}$ (phase chosen real):
$$\psi_n(x) = \sqrt{\frac2a}\,\sin\frac{n\pi x}{a}. \tag{3.3}$$

■ **What just happened.** The string's fitting condition from §2, with $h^2/8 = \pi^2\hbar^2/2$. Quantization came from the second boundary condition: one free constant ($E$) was tuned to satisfy it.

**Four properties.** (i) About the center $x = a/2$, the $\psi_n$ alternate even ($n$ odd) and odd ($n$ even). (ii) $\psi_n$ has $n - 1$ interior nodes. (iii) **Orthonormality**:
$$\langle\psi_m|\psi_n\rangle = \delta_{mn}. \tag{3.4}$$
*Proof.* For $m \ne n$, $\frac2a\int_0^a\sin\frac{m\pi x}{a}\sin\frac{n\pi x}{a}dx = \frac1a\int_0^a\big[\cos\frac{(m-n)\pi x}{a} - \cos\frac{(m+n)\pi x}{a}\big]dx = \frac1\pi\Big[\frac{\sin(m-n)\pi}{m-n} - \frac{\sin(m+n)\pi}{m+n}\Big] = 0$. For $m = n$ it is step 5. ∎ (This is also QM-2 Theorem 3.10(c).) (iv) **Completeness** (stated): any $f$ on $[0, a]$ with $\int|f|^2 < \infty$ is $\sum c_n\psi_n$ in the $L^2$ sense. This is the Fourier sine series of F2; the infinite well's eigenfunctions *are* F2's sine basis, with $\sqrt{2/a}$ for normalization.

**The time-dependent solution.** By QM-2 (3.29),
$$\Psi(x, t) = \sum_{n=1}^\infty c_n\sqrt{\frac2a}\sin\frac{n\pi x}{a}\,e^{-iE_nt/\hbar}, \qquad c_n = \langle\psi_n|\Psi(\cdot, 0)\rangle = \sqrt{\frac2a}\int_0^a\sin\frac{n\pi x}{a}\Psi(x, 0)\,dx. \tag{3.5}$$
The formula for $c_n$ is "Fourier's trick": take the inner product of both sides with $\psi_m$ and use (3.4), the same move as LA-5's coefficients in an orthonormal basis. $\sum|c_n|^2 = 1$ and $\langle H\rangle = \sum|c_n|^2E_n$ (QM-2 (3.30)). Example 4.1 does a full case.

**Sloshing and revivals.** A two-state superposition $(\psi_1 + \psi_2)/\sqrt2$ sloshes left and right at $\omega_{21} = (E_2 - E_1)/\hbar = 3\pi^2\hbar/2ma^2$ (QM-2 Example 4.5). A general state does something more surprising. Every phase is $E_nt/\hbar = n^2\pi^2\hbar t/2ma^2$. At
$$T_{\rm rev} = \frac{4ma^2}{\pi\hbar} \tag{3.6}$$
this equals $2\pi n^2$, a multiple of $2\pi$ for every $n$, so $\Psi(x, T_{\rm rev}) = \Psi(x, 0)$ exactly: any state of the box, however scrambled it looks in between, reassembles perfectly. The reason is that the spectrum is commensurate ($E_n/E_1 = n^2$ are integers). States containing only odd $n$ (such as Example 4.1's) already revive at $T_{\rm rev}/2$, where every phase is $\pi n^2 \equiv \pi$, a global sign; for $(\psi_1 + \psi_2)/\sqrt2$ the overlap at $T_{\rm rev}/2$ is $0$ (the project checks both).

### 3.2 The harmonic oscillator: the algebraic method

$V = \frac12m\omega^2x^2$, so
$$\hat H = \frac{1}{2m}\big[\hat p^2 + (m\omega\hat x)^2\big].$$

**The canonical commutator.** Everything rests on one identity.

▶ **Try it first.** Apply $\hat x\hat p - \hat p\hat x$ to an arbitrary test function $f(x)$.

1. $\hat x\hat pf = x(-i\hbar f')$.
2. $\hat p\hat xf = -i\hbar(xf)' = -i\hbar f - i\hbar xf'$.
3. Subtract: $(\hat x\hat p - \hat p\hat x)f = i\hbar f$ for every $f$. An operator identity is a statement about its action on all test functions, so
$$[\hat x, \hat p] = i\hbar. \tag{3.7}$$

■ **What just happened.** The product rule, nothing more. The test function matters: operators are only defined by what they do to functions, and $[\hat x, \hat p]$ is the operator "multiply by $i\hbar$."

**Factoring the Hamiltonian.** Define
$$a_\pm = \frac{1}{\sqrt{2\hbar m\omega}}\big(\mp i\hat p + m\omega\hat x\big). \tag{3.8}$$

▶ **Try it first.** Multiply out $a_-a_+$ and use (3.7) to eliminate $\hat x\hat p - \hat p\hat x$.

1. $a_-a_+ = \frac{1}{2\hbar m\omega}(i\hat p + m\omega\hat x)(-i\hat p + m\omega\hat x) = \frac{1}{2\hbar m\omega}\big[\hat p^2 + m^2\omega^2\hat x^2 - im\omega(\hat x\hat p - \hat p\hat x)\big]$.
2. The bracket's last term is $-im\omega(i\hbar) = m\omega\hbar$. So $a_-a_+ = \frac{1}{\hbar\omega}\hat H + \frac12$.
3. The same expansion of $a_+a_-$ has $+im\omega(\hat x\hat p - \hat p\hat x) = -m\omega\hbar$: $a_+a_- = \frac{1}{\hbar\omega}\hat H - \frac12$.
4. Hence
$$\hat H = \hbar\omega\Big(a_+a_- + \frac12\Big) = \hbar\omega\Big(a_-a_+ - \frac12\Big), \qquad [a_-, a_+] = 1. \tag{3.9}$$

■ **What just happened.** For numbers, $u^2 + v^2 = (v - iu)(v + iu)$ exactly. For operators the factorization is off by the commutator, and the leftover $\frac12$ is the zero-point energy.

**The ladder property.**

▶ **Try it first.** If $\hat H\psi = E\psi$, compute $\hat H(a_+\psi)$ using (3.9) and $a_-a_+ = a_+a_- + 1$.

1. $\hat H(a_+\psi) = \hbar\omega(a_+a_- + \frac12)a_+\psi = \hbar\omega\,a_+(a_-a_+ + \frac12)\psi = \hbar\omega\,a_+(a_+a_- + 1 + \frac12)\psi = a_+(\hat H + \hbar\omega)\psi = (E + \hbar\omega)(a_+\psi)$.
2. $\hat H(a_-\psi) = \hbar\omega(a_-a_+ - \frac12)a_-\psi = \hbar\omega\,a_-(a_+a_- - \frac12)\psi = a_-(\hat H - \hbar\omega)\psi = (E - \hbar\omega)(a_-\psi)$.

So
$$\hat H(a_\pm\psi) = (E \pm \hbar\omega)(a_\pm\psi). \tag{3.10}$$

■ **What just happened.** $a_+$ maps an eigenfunction to another eigenfunction one quantum higher (or to zero); $a_-$ one quantum lower. No differential equation was solved.

**The ground state.** Lowering repeatedly would produce energies below $V_{\min} = 0$, forbidden by QM-2 Theorem 3.10(b). So the ladder must stop: there is a state with $a_-\psi_0 = 0$ (zero is not an eigenfunction, so the chain ends there).

▶ **Try it first.** Write $a_-\psi_0 = 0$ as a first-order ODE and solve it.

1. $\frac{1}{\sqrt{2\hbar m\omega}}\big(\hbar\frac{d}{dx} + m\omega x\big)\psi_0 = 0$, i.e. $\psi_0'/\psi_0 = -m\omega x/\hbar$.
2. Integrate: $\psi_0 = A e^{-m\omega x^2/2\hbar}$. Normalize with $\int e^{-m\omega x^2/\hbar}dx = \sqrt{\pi\hbar/m\omega}$:
$$\psi_0(x) = \Big(\frac{m\omega}{\pi\hbar}\Big)^{1/4}e^{-m\omega x^2/2\hbar}, \qquad E_0 = \frac12\hbar\omega. \tag{3.11}$$
3. $E_0$: by (3.9), $\hat H\psi_0 = \hbar\omega(a_+a_-\psi_0 + \frac12\psi_0) = \frac12\hbar\omega\psi_0$.
4. A first-order linear ODE has a one-dimensional solution space, so $\psi_0$ is the *only* bottom rung. Every eigenstate lowers in finitely many steps to $\psi_0$; so every energy is $E_0 + n\hbar\omega$ and every eigenstate is $\propto a_+^n\psi_0$:
$$E_n = \Big(n + \frac12\Big)\hbar\omega, \qquad n = 0, 1, 2, \dots \tag{3.12}$$

■ **What just happened.** The Gaussian is QM-2 (3.16) with $\sigma_x^2 = \hbar/2m\omega$, the minimum-uncertainty state. The whole spectrum came from the commutator plus "energy is bounded below."

**Normalization constants.** We need $\langle f|a_\pm g\rangle = \langle a_\mp f|g\rangle$: the adjoint of $a_+$ is $a_-$. This holds because $\hat x$ and $\hat p$ are Hermitian (QM-2 Theorem 3.4) and conjugating the coefficient $-i$ gives $+i$: $\langle f|(-i\hat p + m\omega\hat x)g\rangle = \langle(i\hat p + m\omega\hat x)f|g\rangle$.

▶ **Try it first.** Write $a_+\psi_n = c_n\psi_{n+1}$ and compute $|c_n|^2 = \langle a_+\psi_n|a_+\psi_n\rangle$ by moving one $a_+$ across.

1. $|c_n|^2 = \langle a_+\psi_n|a_+\psi_n\rangle = \langle\psi_n|a_-a_+\psi_n\rangle$.
2. By (3.9), $a_-a_+ = \hat H/\hbar\omega + \frac12$, which acts on $\psi_n$ as $n + 1$. So $|c_n|^2 = n + 1$.
3. For the lowering constant $d_n$ in $a_-\psi_n = d_n\psi_{n-1}$: $|d_n|^2 = \langle\psi_n|a_+a_-\psi_n\rangle = n$, since $a_+a_- = \hat H/\hbar\omega - \frac12$ acts as $n$.
4. Choosing real positive phases:
$$a_+\psi_n = \sqrt{n+1}\,\psi_{n+1}, \qquad a_-\psi_n = \sqrt n\,\psi_{n-1}, \qquad \psi_n = \frac{1}{\sqrt{n!}}(a_+)^n\psi_0. \tag{3.13}$$

■ **What just happened.** $\hat N = a_+a_-$ is the **number operator**, $\hat N\psi_n = n\psi_n$, $\hat H = \hbar\omega(\hat N + \frac12)$. The $\sqrt{n+1}$ and $\sqrt n$ are forced by norms. In the basis $\{\psi_n\}$, $a_-$ is the matrix with $\sqrt1, \sqrt2, \sqrt3, \dots$ immediately above the diagonal and $a_+$ its transpose; the project builds them and checks $[a_-, a_+] = 1$.

**Orthogonality from the ladder.** $\langle\psi_m|a_+a_-\psi_n\rangle = n\langle\psi_m|\psi_n\rangle$; moving both operators across, $\langle a_+a_-\psi_m|\psi_n\rangle = m\langle\psi_m|\psi_n\rangle$. Subtract: $(n - m)\langle\psi_m|\psi_n\rangle = 0$, so $\langle\psi_m|\psi_n\rangle = 0$ for $m \ne n$. ∎ (QM-2 Theorem 3.10(c) again, in ladder language.)

**Expectation values.** Adding and subtracting the two lines of (3.8):
$$\hat x = \sqrt{\frac{\hbar}{2m\omega}}(a_+ + a_-), \qquad \hat p = i\sqrt{\frac{\hbar m\omega}{2}}(a_+ - a_-). \tag{3.14}$$
$\hat x\psi_n$ and $\hat p\psi_n$ are combinations of $\psi_{n\pm1}$, which are orthogonal to $\psi_n$: $\langle x\rangle = \langle p\rangle = 0$ in every eigenstate. Second moments need two ladder steps that return to $n$; Example 4.2 does them and finds $\langle T\rangle = \langle V\rangle = \frac12E_n$ (the oscillator's virial theorem, CM-4 Theorem 3.9 with $V\propto x^2$) and $\sigma_x\sigma_p = (n + \frac12)\hbar$.

### 3.3 The harmonic oscillator: the analytic method

Now solve the differential equation directly. It is longer, but the argument that quantizes the energy is one you have seen before and will see again.

▶ **Try it first.** Substitute $\xi = \sqrt{m\omega/\hbar}\,x$ and $K = 2E/\hbar\omega$, and find the equation for $h$ in $\psi = h(\xi)e^{-\xi^2/2}$.

1. The equation becomes $\frac{d^2\psi}{d\xi^2} = (\xi^2 - K)\psi$. For large $\xi$, $\psi \approx e^{\pm\xi^2/2}$; only the minus sign is normalizable, which motivates the substitution.
2. With $\psi = he^{-\xi^2/2}$: $\psi'' = (h'' - 2\xi h' + (\xi^2 - 1)h)e^{-\xi^2/2}$. So
$$h'' - 2\xi h' + (K - 1)h = 0 \qquad\text{(Hermite's equation).} \tag{3.15}$$
3. Power series $h = \sum_ja_j\xi^j$. The coefficient of $\xi^j$: $(j+2)(j+1)a_{j+2} - 2ja_j + (K - 1)a_j = 0$:
$$a_{j+2} = \frac{2j + 1 - K}{(j+1)(j+2)}\,a_j. \tag{3.16}$$
Even and odd coefficients decouple: $a_0$ seeds an even series, $a_1$ an odd one.
4. **Termination.** If the series does not stop, then for large $j$, $a_{j+2}/a_j \approx 2/j$. The series for $e^{\xi^2} = \sum\xi^{2k}/k!$ has ratio $\frac{1/(k+1)!}{1/k!} = \frac{1}{k+1} \approx \frac2j$ with $j = 2k$. So $h$ grows like $e^{\xi^2}$ and $\psi = he^{-\xi^2/2}$ like $e^{+\xi^2/2}$: not normalizable. The series must terminate: for some $n$, $2n + 1 - K = 0$, and the other-parity series must be zero ($a_1 = 0$ for even $n$, $a_0 = 0$ for odd $n$).
5. $K = 2n + 1$ is $E_n = (n + \frac12)\hbar\omega$, (3.12) again.

■ **What just happened.** This is the Frobenius argument of EM-3 §3.5: Legendre's equation has a series that diverges at $x = \pm1$ unless it terminates, which forces $\ell(\ell+1)$ and the Legendre polynomials. Here, divergence at $\xi\to\infty$ forces $K = 2n+1$. Same structure, same quantization mechanism: a physical requirement (finite, normalizable) chops an infinite series into a polynomial. Notice it next time; it happens again for hydrogen (QM-6).

**Hermite polynomials.** The terminating $h$ with leading coefficient $2^n$ is $H_n(\xi)$: $H_0 = 1$, $H_1 = 2\xi$, $H_2 = 4\xi^2 - 2$, $H_3 = 8\xi^3 - 12\xi$. In $\xi$, $a_\pm = \frac{1}{\sqrt2}(\xi \mp \frac{d}{d\xi})$. Applying $a_+$ to $H_ne^{-\xi^2/2}$ gives $\frac{1}{\sqrt2}(2\xi H_n - H_n')e^{-\xi^2/2}$, and applying $a_-$ gives $\frac{1}{\sqrt2}H_n'e^{-\xi^2/2}$. Matching to (3.13) with the normalization below gives $H_{n+1} = 2\xi H_n - H_n'$ and $H_n' = 2nH_{n-1}$, which combine into the recurrence
$$H_{n+1}(\xi) = 2\xi H_n(\xi) - 2nH_{n-1}(\xi), \qquad \psi_n(x) = \Big(\frac{m\omega}{\pi\hbar}\Big)^{1/4}\frac{1}{\sqrt{2^nn!}}H_n(\xi)e^{-\xi^2/2}. \tag{3.17}$$
The project builds $H_n$ by (3.17) and checks against `scipy.special.hermite`.

**What the states look like.** $\psi_n$ has $n$ nodes and parity $(-1)^n$ and decays beyond the turning points $\pm A_n$, $A_n = \sqrt{(2n+1)\hbar/m\omega}$. The classical density at energy $E_n$ is $P_{\rm cl}(x) = 1/(\pi\sqrt{A_n^2 - x^2})$ (time spent $\propto 1/|v|$), peaked at the turning points; $|\psi_0|^2$ is peaked at the center, the opposite. For large $n$, $|\psi_n|^2$ oscillates about $P_{\rm cl}$ and its local average approaches it (CM-8 §3.10's WKB result $|\psi|^2\propto1/p$, checked there for $n = 20$).

**Why it matters.** Every smooth minimum is an oscillator to second order, and every normal mode of CM-6 quantizes to a ladder with spacing $\hbar\omega$: phonons and photons (SM-5). The zero-point energy is real; it is why helium stays liquid at absolute zero under its own vapor pressure.

### 3.4 The free particle, and bound vs scattering states

For $V = 0$ everywhere, $\psi = Ae^{ikx} + Be^{-ikx}$ with $E = \hbar^2k^2/2m$ for any $k$; QM-2 §3.7–3.8 built normalizable packets from these and found $v_g = \hbar k/m$. Plane waves are not normalizable, but they idealize a long packet with sharp momentum (project Part 7).

**Definitions.** For a potential with limits $V(\pm\infty)$: a **bound state** has $E < V(-\infty)$ and $E < V(+\infty)$; it decays on both sides, is normalizable, and occurs only at discrete $E$. A **scattering state** has $E$ above at least one of the limits; it oscillates at infinity, is not normalizable, and exists for a continuum of $E$. For $V\to0$ at $\pm\infty$: $E < 0$ bound, $E > 0$ scattering.

**Reflection and transmission from currents.** Send a wave in from the left. Far left, $\psi = Ae^{ik_Lx} + Be^{-ik_Lx}$; far right, $\psi = Fe^{ik_Rx}$. The current of $Ae^{ikx}$ is $J = \hbar k|A|^2/m$ (QM-2 Example 4.2): probability density $|A|^2$ moving at speed $\hbar k/m$.

▶ **Try it first.** Show $J$ is independent of $x$ in any stationary state, then compute $J$ on both sides.

1. In a stationary state $\rho = |\psi|^2$ does not depend on $t$, so the continuity equation $\partial_t\rho + \partial_xJ = 0$ (QM-2 (3.5)) gives $\partial_xJ = 0$.
2. Left: $J = \frac{\hbar}{m}\operatorname{Im}(\psi^*\psi')$ with $\psi = Ae^{ik_Lx} + Be^{-ik_Lx}$. The cross terms are real (QM-2 Exercise B3), so $J_L = \frac{\hbar k_L}{m}(|A|^2 - |B|^2)$.
3. Right: $J_R = \frac{\hbar k_R}{m}|F|^2$.
4. $J_L = J_R$. Divide by the incident current $\hbar k_L|A|^2/m$:
$$R + T = 1, \qquad R = \frac{|B|^2}{|A|^2}, \qquad T = \frac{k_R}{k_L}\frac{|F|^2}{|A|^2}. \tag{3.18}$$

■ **What just happened.** $R$ and $T$ are ratios of **fluxes** (probability per unit time), not of densities. When the speeds differ, a transmitted wave of larger amplitude can carry less flux. The project's packet calculation confirms that these are the probabilities a real packet splits into.

### 3.5 The delta-function potential

$V(x) = -\alpha\delta(x)$, $\alpha > 0$: an infinitely narrow, infinitely deep well of finite "strength" (F3 treats $\delta$ as the functional $\int\delta f = f(0)$). Away from $0$, $V = 0$.

▶ **Try it first.** Integrate the Schrödinger equation from $-\epsilon$ to $\epsilon$ and let $\epsilon\to0$.

1. $-\frac{\hbar^2}{2m}\big[\psi'(\epsilon) - \psi'(-\epsilon)\big] + \int_{-\epsilon}^\epsilon V\psi\,dx = E\int_{-\epsilon}^\epsilon\psi\,dx$.
2. The right side goes to $0$ ($\psi$ is bounded). For a finite jump in $V$, the $V$ integral goes to $0$ too, so $\psi'$ is continuous; that is why finite wells require $\psi'$ continuity. For the delta, $\int V\psi = -\alpha\psi(0)$. So
$$\Delta\psi'(0) \equiv \psi'(0^+) - \psi'(0^-) = -\frac{2m\alpha}{\hbar^2}\psi(0). \tag{3.19}$$
$\psi$ itself stays continuous (a jump in $\psi$ would make $\psi''$ contain $\delta'$, which nothing in the equation can match).

**The bound state.** $E < 0$, $\kappa = \sqrt{-2mE}/\hbar$. Away from $0$, $\psi = $ combinations of $e^{\pm\kappa x}$; normalizability keeps $e^{\kappa x}$ for $x < 0$ and $e^{-\kappa x}$ for $x > 0$; continuity makes the two coefficients equal: $\psi = Be^{-\kappa|x|}$.

1. $\Delta\psi' = -\kappa B - \kappa B = -2\kappa B$. With (3.19), $-2\kappa B = -\frac{2m\alpha}{\hbar^2}B$, so $\kappa = m\alpha/\hbar^2$.
2. Energy and normalization ($\int B^2e^{-2\kappa|x|}dx = B^2/\kappa$):
$$E_b = -\frac{m\alpha^2}{2\hbar^2}, \qquad \psi(x) = \frac{\sqrt{m\alpha}}{\hbar}e^{-m\alpha|x|/\hbar^2}. \tag{3.20}$$
There is exactly one bound state, whatever $\alpha$.

**Scattering.** $E > 0$. Left: $Ae^{ikx} + Be^{-ikx}$; right: $Fe^{ikx}$. Let $\beta = m\alpha/\hbar^2k$.

▶ **Try it first.** Two conditions (continuity and (3.19)), two unknowns ($B$, $F$ in units of $A$).

1. Continuity: $A + B = F$.
2. Jump: $ikF - ik(A - B) = -\frac{2m\alpha}{\hbar^2}F$, i.e. $F - A + B = 2i\beta F$.
3. Substitute $B = F - A$: $2F - 2A = 2i\beta F$, so $F = \dfrac{A}{1 - i\beta}$ and $B = F - A = \dfrac{i\beta A}{1 - i\beta}$.
4. Same $k$ on both sides, so $T = |F/A|^2$:
$$R = \frac{\beta^2}{1+\beta^2}, \qquad T = \frac{1}{1+\beta^2} = \frac{1}{1 + m\alpha^2/2\hbar^2E}. \tag{3.21}$$

■ **What just happened.** Low energies reflect, high energies transmit, crossover at $E\sim|E_b|$. The barrier $+\alpha\delta(x)$ has $\beta \to -\beta$ and identical $R$ and $T$: they depend on $\alpha^2$. Classically a barrier reflects everything below its top and a well reflects nothing; quantum mechanically the delta well reflects exactly as much as the delta barrier. (The bound state is visible in the scattering solution: the transmission amplitude $F/A = 1/(1 - i\beta)$ has a pole at $\beta = -i$, i.e. $k = i m\alpha/\hbar^2 = i\kappa$. Bound states are poles of $t(E)$ continued to $E < 0$, a theme of QM-8's scattering theory.)

### 3.6 The finite square well

$V = -V_0$ for $|x| < a$, $0$ outside. Bound states have $-V_0 < E < 0$. Let $l = \sqrt{2m(E + V_0)}/\hbar$ inside and $\kappa = \sqrt{-2mE}/\hbar$ outside. $V$ is even, so (§3.11) we can take the solutions even or odd.

▶ **Try it first.** For the even solution, $\psi = D\cos lx$ inside and $Fe^{-\kappa x}$ for $x > a$. Impose continuity of $\psi$ and $\psi'$ at $x = a$ and eliminate $D/F$.

1. Continuity: $Fe^{-\kappa a} = D\cos la$. Derivative: $-\kappa Fe^{-\kappa a} = -lD\sin la$.
2. Divide: $\kappa = l\tan la$.
3. Dimensionless variables: $z = la$, $z_0 = a\sqrt{2mV_0}/\hbar$. Since $\kappa^2 + l^2 = 2mV_0/\hbar^2$, $\kappa a = \sqrt{z_0^2 - z^2}$. So
$$\text{even:}\quad \tan z = \sqrt{(z_0/z)^2 - 1}, \qquad\qquad \text{odd:}\quad -\cot z = \sqrt{(z_0/z)^2 - 1}. \tag{3.22}$$
4. Odd: $\psi = D\sin lx$ inside gives $-\kappa = l\cot la$, the second equation.
5. Energies: $E + V_0 = \hbar^2l^2/2m = \hbar^2z^2/2ma^2$, i.e. $E/V_0 = z^2/z_0^2 - 1$.

■ **What just happened.** Quantization again from matching: the outside decay rate must equal the inside "log-derivative." (3.22) is **transcendental**: there is no closed-form solution, and none will appear if you try harder. Solve it graphically or numerically.

**The graphical solution.** Plot $\sqrt{(z_0/z)^2 - 1}$ (decreasing from $+\infty$ at $z \to 0$ to $0$ at $z = z_0$) against $\tan z$ (positive and increasing on each $[n\pi, n\pi + \pi/2)$) and $-\cot z$ (positive and increasing on each $[n\pi + \pi/2, (n+1)\pi)$). Each branch that starts below $z_0$ crosses the curve exactly once. Even and odd states alternate, and the number of bound states is $\lfloor 2z_0/\pi\rfloor + 1$ (when $2z_0/\pi$ is not an integer). The project's `finite_well.png` is this picture for $z_0 = 8$; Example 4.5 gives the numbers.

**Wide and deep.** For large $z_0$, the curve is large over most of $(0, z_0)$, so the crossings sit slightly below the poles of $\tan$ and $-\cot$, at $z \approx n\pi/2$. Then $l\cdot2a \approx n\pi$ and $E_n + V_0 \approx n^2\pi^2\hbar^2/2m(2a)^2$: the infinite well of width $2a$, measured from the bottom. The true $z_n$ are a little smaller, because the wavefunction leaks out and "sees" a wider box.

**Shallow: there is always a bound state.**

**Theorem 3.1.** For every $z_0 > 0$ the finite well has at least one bound state (even). If $z_0 < \pi/2$, it has exactly one.

*Proof.* Let $g(z) = \tan z - \sqrt{z_0^2/z^2 - 1}$ on $(0, z_*)$, $z_* = \min(z_0, \pi/2)$. As $z\to0^+$, $g\to-\infty$. If $z_0 \le \pi/2$, $g(z_0) = \tan z_0 > 0$; if $z_0 > \pi/2$, $g\to+\infty$ as $z\to\pi/2^-$. $g$ is continuous on $(0, z_*)$, so by the intermediate value theorem it has a root: an even bound state. If $z_0 < \pi/2$: $g$ is strictly increasing on $(0, z_0)$ ($\tan$ increases, the square root decreases), so the root is unique; and an odd state needs $-\cot z > 0$, i.e. $z \in (\pi/2, \pi)$ modulo $\pi$, impossible for $z < z_0 < \pi/2$. ∎

For small $z_0$ the binding is tiny, $\kappa a \approx z_0^2$ (Exercise B5), and the state spreads far outside, but it exists. In 1D every attractive well binds (QM-8, variationally); in 3D a spherical well needs a minimum depth (QM-6).

Scattering off the finite well, with its transmission resonances, comes out of the barrier calculation in one line; it is in §3.8.

### 3.7 The potential step

$V = 0$ for $x < 0$, $V_0 > 0$ for $x > 0$.

**$E > V_0$: partial reflection.** $k = \sqrt{2mE}/\hbar$ on the left, $k' = \sqrt{2m(E - V_0)}/\hbar$ on the right. $\psi = Ae^{ikx} + Be^{-ikx}$ ($x < 0$), $Fe^{ik'x}$ ($x > 0$).

▶ **Try it first.** Match $\psi$ and $\psi'$ at $0$, find $B/A$ and $F/A$, and then compute $R$ and $T$ using (3.18). Then compute $|F/A|^2$ and compare.

1. $A + B = F$ and $ik(A - B) = ik'F$.
2. Add $k\times$(first) to the second divided by $i$: $2kA = (k + k')F$, so $F/A = 2k/(k + k')$ and $B/A = F/A - 1 = (k - k')/(k + k')$.
3. By (3.18):
$$R = \Big(\frac{k - k'}{k + k'}\Big)^2, \qquad T = \frac{k'}{k}\Big(\frac{2k}{k + k'}\Big)^2 = \frac{4kk'}{(k + k')^2}. \tag{3.23}$$
4. Check: $(k - k')^2 + 4kk' = (k + k')^2$, so $R + T = 1$.
5. But $|F/A|^2 = 4k^2/(k + k')^2 > 1$ whenever $k' < k$.

■ **What just happened.** The transmitted wave is *taller* than the incident one ($F/A > 1$), yet it carries less flux, because it moves slower ($k' < k$). Probability piles up where the particle is slow, exactly as for the classical density $1/v$. $T = |F/A|^2$ would say more than 100% transmits. The current, not the amplitude, carries the flux. A particle with $E > V_0$ classically always clears the step; quantum mechanically it reflects with probability $R > 0$, and also reflects from a step *down* ($V_0 < 0$, $k' > k$), since $R$ depends only on $(k - k')^2$.

**$E < V_0$: total reflection with penetration.** On the right $\psi = Fe^{-\kappa x}$. Matching replaces $ik'$ by $-\kappa$: $B/A = (k - i\kappa)/(k + i\kappa)$, a ratio of complex conjugates, so $|B/A| = 1$ and $R = 1$. The wave penetrates a distance $\sim1/\kappa$ but carries no current there ($Fe^{-\kappa x}$ is a real function times a constant, $J = 0$). The reflected wave picks up the phase $-2\arctan(\kappa/k)$.

### 3.8 The rectangular barrier and tunneling

$V = V_0$ for $0 < x < a$, $0$ elsewhere, $0 < E < V_0$. Left: $Ae^{ikx} + Be^{-ikx}$; inside: $Ce^{\kappa x} + De^{-\kappa x}$; right: $Fe^{ikx}$. Four matching conditions (ψ and ψ' at $0$ and at $a$), four unknowns ($B, C, D, F$ in units of $A$). You can eliminate by brute force (Griffiths Problem 2.33; the project does it with SymPy). That algebra is where errors live, so we organize it as matrices.

**The propagator of $(\psi, \psi')$.** In a region where $V$ is constant, $\psi'' = -q^2\psi$ with $q^2 = 2m(E - V)/\hbar^2$ ($q$ real if $E > V$, $q = i\kappa$ if $E < V$). Given $\psi$ and $\psi'$ at the left edge, the solution at distance $d$ is $\psi(d) = \psi(0)\cos qd + \psi'(0)\sin(qd)/q$, and differentiating, $\psi'(d) = -q\psi(0)\sin qd + \psi'(0)\cos qd$. So
$$\begin{pmatrix}\psi\\ \psi'\end{pmatrix}_{\!d} = P_q(d)\begin{pmatrix}\psi\\ \psi'\end{pmatrix}_{\!0}, \qquad P_q(d) = \begin{pmatrix}\cos qd & \sin(qd)/q\\ -q\sin qd & \cos qd\end{pmatrix}, \qquad \det P_q = 1. \tag{3.24}$$
For $q = i\kappa$: $\cos(i\kappa d) = \cosh\kappa d$, $\sin(i\kappa d)/(i\kappa) = \sinh(\kappa d)/\kappa$, $-q\sin qd = \kappa\sinh\kappa d$. Because $\psi$ and $\psi'$ are continuous at each interface, crossing several regions means multiplying their $P$'s (rightmost first, LA-1/LA-4's composition).

**From $(\psi, \psi')$ to amplitudes.** In a region with real wavenumber $k$, $\psi = Pe^{ikx} + Qe^{-ikx}$ gives
$$\begin{pmatrix}\psi\\ \psi'\end{pmatrix} = W_k(x)\begin{pmatrix}P\\ Q\end{pmatrix}, \qquad W_k(x) = \begin{pmatrix}e^{ikx} & e^{-ikx}\\ ike^{ikx} & -ike^{-ikx}\end{pmatrix}, \qquad \det W_k = -2ik. \tag{3.25}$$

▶ **Try it first.** For the barrier, write the amplitudes on the right (referenced to $x = a$, i.e. $\psi = \tilde Fe^{ik(x - a)}$) as $M$ times $(A, B)$, with $M = W_k(0)^{-1}P_{i\kappa}(a)W_k(0)$, and compute $M_{22}$.

1. $(\tilde F, 0)^\top = M(A, B)^\top$. The second row gives $0 = M_{21}A + M_{22}B$, so $B/A = -M_{21}/M_{22}$. The first row then gives $\tilde F/A = (M_{11}M_{22} - M_{12}M_{21})/M_{22} = \det M/M_{22}$, and $\det M = \det P = 1$. So $\tilde F/A = 1/M_{22}$ and $T = 1/|M_{22}|^2$.
2. $W_k(0)^{-1} = \begin{pmatrix}\frac12 & \frac{1}{2ik}\\ \frac12 & -\frac{1}{2ik}\end{pmatrix}$. With $P = \begin{pmatrix}p_{11} & p_{12}\\ p_{21} & p_{22}\end{pmatrix}$, the $(2,2)$ entry of $W^{-1}PW$ is $\frac12(p_{11} + p_{22}) - \frac{ik}{2}p_{12} - \frac{1}{2ik}p_{21}$.
3. Insert $p_{11} = p_{22} = \cosh\kappa a$, $p_{12} = \sinh(\kappa a)/\kappa$, $p_{21} = \kappa\sinh\kappa a$:
$$M_{22} = \cosh\kappa a + i\,\frac{\kappa^2 - k^2}{2k\kappa}\sinh\kappa a.$$
(The SymPy elimination in the project gives $F/A = e^{-ika}/M_{22}$, the same thing with the phase reference moved to $x = 0$.)
4. $|M_{22}|^2 = \cosh^2\kappa a + \frac{(\kappa^2 - k^2)^2}{4k^2\kappa^2}\sinh^2\kappa a = 1 + \Big[1 + \frac{(\kappa^2 - k^2)^2}{4k^2\kappa^2}\Big]\sinh^2\kappa a = 1 + \frac{(\kappa^2 + k^2)^2}{4k^2\kappa^2}\sinh^2\kappa a$.
5. $\kappa^2 + k^2 = 2mV_0/\hbar^2$ and $k^2\kappa^2 = 4m^2E(V_0 - E)/\hbar^4$, so the coefficient is $V_0^2/4E(V_0 - E)$:
$$T = \Big[1 + \frac{V_0^2\sinh^2(\kappa a)}{4E(V_0 - E)}\Big]^{-1}. \tag{3.26}$$

■ **What just happened.** The barrier is "propagate $(\psi, \psi')$ through a forbidden region," a $2\times 2$ hyperbolic rotation, sandwiched between two changes of variables to amplitudes. Classically $T = 0$. Quantum mechanically $T > 0$ for any finite barrier: **tunneling**.

**The thick barrier.** For $\kappa a \gg 1$, $\sinh^2\kappa a \approx e^{2\kappa a}/4$ and the 1 is negligible:
$$T \approx \frac{16E(V_0 - E)}{V_0^2}\,e^{-2\kappa a}. \tag{3.27}$$
The prefactor is at most 4; the exponential does everything. $\ln T$ is linear in $a$ with slope $-2\kappa$ (the project fits it: $-1.99988$ against $-2$). Double the width and you square the (small) transmission; Example 4.6 makes this concrete.

**Above the barrier and the finite well's resonances.** For $E > V_0$, replace $\kappa$ by $il$ with $l = \sqrt{2m(E - V_0)}/\hbar$: $\sinh(\kappa a) \to i\sin(la)$ and $V_0 - E \to -(E - V_0)$, so $T = [1 + V_0^2\sin^2(la)/4E(E - V_0)]^{-1}$. Now $T = 1$ whenever $la = n\pi$: a whole number of half-wavelengths fits over the barrier and the reflections from the two edges cancel. For the **finite well** of depth $V_0$ and width $w$ ($= 2a$ in §3.6), replace $V_0$ by $-V_0$:
$$T = \Big[1 + \frac{V_0^2\sin^2(lw)}{4E(E + V_0)}\Big]^{-1}, \qquad l = \frac{\sqrt{2m(E + V_0)}}{\hbar}, \qquad T = 1 \text{ at } E + V_0 = \frac{n^2\pi^2\hbar^2}{2mw^2}. \tag{3.28}$$
The well becomes perfectly transparent at the energies of the infinite well of the same width (measured from its bottom). This is the 1D model of the **Ramsauer–Townsend effect** (Ramsauer 1921; Townsend and Bailey 1922): electrons of about 1 eV cross argon and the heavier noble gases almost unscattered. The real effect is 3D, but the cancellation mechanism is this one.

**WKB (stated).** For a smooth barrier $V(x)$, QM-8 derives
$$T \approx \exp\Big(-2\int_{x_1}^{x_2}\kappa(x)\,dx\Big), \qquad \kappa(x) = \frac{\sqrt{2m(V(x) - E)}}{\hbar}, \tag{3.29}$$
with $x_1, x_2$ the classical turning points. For a rectangle it reproduces the exponent of (3.27). It is (3.27) applied to thin slices, multiplied together, with the prefactors dropped.

**Applications with numbers.**

*Alpha decay (Gamow 1928).* An alpha particle of energy $E$ sits inside a nucleus of radius $r_1$; outside, the daughter (charge $Ze$) repels it with $V(r) = 2Zk_ee^2/r$, $k_e = 1/4\pi\epsilon_0$. The barrier extends from $r_1$ to the turning point $r_2 = 2Zk_ee^2/E$. With (3.29), $T = e^{-2\gamma}$, and the integral can be done exactly (Exercise C3):
$$\gamma = \frac1\hbar\int_{r_1}^{r_2}\sqrt{2m\big(V(r) - E\big)}\,dr = \frac{\sqrt{2mE}}{\hbar}\,r_2\Big[\arccos\sqrt{\tfrac{r_1}{r_2}} - \sqrt{\tfrac{r_1}{r_2}\big(1 - \tfrac{r_1}{r_2}\big)}\Big], \qquad t_{1/2} = \ln2\cdot\frac{2r_1}{v}\,e^{2\gamma}. \tag{3.30}$$
Lifetime model: the alpha hits the wall every $2r_1/v$, $v = \sqrt{2E/m}$, and escapes with probability $e^{-2\gamma}$ per hit. For $r_1 \ll r_2$, $\gamma \approx \frac{\pi}{2}\frac{\sqrt{2mE}}{\hbar}r_2 \propto Z/\sqrt E$: $\ln t_{1/2}$ is linear in $Z/\sqrt E$, the form of Geiger and Nuttall's empirical rule. Example 4.8 puts in numbers.

*The STM.* An electron at the Fermi level of a metal sees a barrier of height about the work function $\phi \approx 4.5$ eV across the vacuum gap $d$. $\kappa = \sqrt{2m_e\phi}/\hbar = 1.087\times10^{10}$ m$^{-1}$, and the current $\propto T \propto e^{-2\kappa d}$ drops by $e^{2\kappa\cdot1\,\text{Å}} = 8.79$ per ångström. Hold the current fixed with feedback and the tip height traces the surface: 0.01 Å of height changes the current by a factor $1.022$, a 2% signal, which is why sub-ångström vertical resolution is routine.

*Field emission (Fowler–Nordheim 1928).* A field $F$ turns the surface barrier into a triangle, $V - E = \phi - eFx$, of length $\phi/eF$. (3.29) gives $2\int_0^{\phi/eF}\frac{\sqrt{2m(\phi - eFx)}}{\hbar}dx = \frac43\frac{\sqrt{2m}\,\phi^{3/2}}{\hbar eF}$. For $\phi = 4.5$ eV and $F = 5$ V/nm, $T\sim e^{-13.04} = 2.2\times10^{-6}$. Halving the field doubles the exponent: $e^{-26.08} = 4.7\times10^{-12}$.

### 3.9 The transfer-matrix method for any piecewise-constant potential

The barrier calculation generalizes with no new ideas. Take regions $j = 0, 1, \dots, N$ with constant $V_j$, interfaces $x_1 < \dots < x_N$, and allow a delta $\alpha_j\delta(x - x_j)$ at each interface.

▶ **Try it first.** What matrix represents crossing a delta $\alpha\,\delta(x - x_j)$, acting on $(\psi, \psi')$? (Use (3.19) with $-\alpha \to +\alpha$.)

1. Inside region $j$ of width $d_j = x_{j+1} - x_j$: multiply by $P_{q_j}(d_j)$ from (3.24).
2. Across a delta: $\psi$ continuous, $\psi'$ jumps by $\frac{2m\alpha_j}{\hbar^2}\psi$: multiply by $D_j = \begin{pmatrix}1 & 0\\ 2m\alpha_j/\hbar^2 & 1\end{pmatrix}$, $\det D_j = 1$.
3. Total in $(\psi, \psi')$ form: $\Phi = D_NP_{N-1}D_{N-1}\cdots P_1D_1$, $\det\Phi = 1$.
4. Convert to amplitudes with (3.25): $M = W_{k_N}(x_N)^{-1}\,\Phi\,W_{k_0}(x_1)$. Its determinant is $\det W_{k_0}/\det W_{k_N} = k_0/k_N$.
5. With no wave incoming from the right, as in the barrier: $r = B/A = -M_{21}/M_{22}$, transmitted amplitude $\det M/M_{22}$, and by (3.18)
$$T = \frac{k_N}{k_0}\Big|\frac{\det M}{M_{22}}\Big|^2 = \frac{k_0/k_N}{|M_{22}|^2}, \qquad R = \Big|\frac{M_{21}}{M_{22}}\Big|^2. \tag{3.31}$$
(If $k_N$ is imaginary, $E$ is below the right-hand potential and $T = 0$.)

■ **What just happened.** Any piecewise-constant potential (and any smooth one, approximated by thin slices) is a product of $2\times 2$ unimodular matrices. Unimodularity is current conservation in matrix form. `transfer()` uses the equivalent form $S_j = W_{j+1}(0)^{-1}D_jW_j(d_j)$, phases referenced to each region's left edge so $e^{\kappa x}$ never overflows.

The project checks `transfer()` against (3.26), (3.23), (3.21) and (3.28), then does what closed forms cannot: the **double barrier**, two barriers with a well between, reaches $T = 1$ at the well's quasi-bound energy although each barrier alone transmits 1.6% there (resonant tunneling).

### 3.10 The double well and tunneling splitting

Two identical wells separated by a barrier of width $b$ and height above the levels set by $\kappa$. Let $\phi_L$, $\phi_R$ be the ground state of each well alone, with energy $E_0$. When $b$ is large, each tail reaches into the other well only as $e^{-\kappa b}$.

Because $V$ is even, eigenstates are even or odd (§3.11). The two lowest are close to
$$\psi_\pm \approx \frac{1}{\sqrt2}(\phi_L \pm \phi_R), \qquad E_\pm = E_0 \mp \frac{\Delta}{2}, \qquad \Delta \propto e^{-\kappa b}. \tag{3.32}$$
The even state $\psi_+$ has no node and is lower (node theorem, §3.11); the odd state has a node in the barrier, where the even state has a small bump. The splitting is the tunneling coupling: the energy shift from $\phi_L$'s tail overlapping well $R$ is $\propto\phi_L(\text{at well } R)\sim e^{-\kappa b}$. QM-8 makes this precise; Exercise C1 does the two-delta well exactly and finds $\Delta \approx 4|E_b|e^{-\kappa d}$ for separation $d$. The project's finite-difference double well gives $d\ln\Delta/db = -3.924$ against $-\kappa = -3.924$.

**Dynamics.** Start in $\phi_L = (\psi_+ + \psi_-)/\sqrt2$: $\Psi(t) \propto \psi_+ + e^{-i\Delta t/\hbar}\psi_-$ reaches $\phi_R$ at $t = \pi\hbar/\Delta$, exponentially slowly in $b$.

**Ammonia.** In NH$_3$, the nitrogen can sit above or below the plane of the three hydrogens: a double well. The inversion splitting is $\Delta = h\nu$ with $\nu = 23.87$ GHz, $\Delta = 9.87\times10^{-5}$ eV, and the first maser (Gordon, Zeiger and Townes, 1954) ran on this transition. QM-4 treats exactly this as a two-level system, $\hat H = E_0\mathbb 1 - \frac{\Delta}{2}\sigma_x$ in the $\{\phi_L, \phi_R\}$ basis; QC-9's superconducting flux qubits are double wells of this kind.

### 3.11 General theorems in one dimension

**Theorem 3.2 (No degeneracy).** Two bound states of a 1D Hamiltonian with the same energy are proportional.

▶ **Try it first.** Multiply the equation for $\psi_1$ by $\psi_2$, the equation for $\psi_2$ by $\psi_1$, and subtract.

*Proof.* 1. Both satisfy $\psi'' = \frac{2m}{\hbar^2}(V - E)\psi$. So $\psi_1''\psi_2 - \psi_2''\psi_1 = 0$.
2. The left side is $\frac{d}{dx}(\psi_1'\psi_2 - \psi_2'\psi_1)$, so the Wronskian $W = \psi_1'\psi_2 - \psi_2'\psi_1$ is constant.
3. Bound states vanish (with their derivatives) at $\pm\infty$, so $W = 0$ everywhere.
4. On any interval where $\psi_1, \psi_2 \ne 0$: $\psi_1'/\psi_1 = \psi_2'/\psi_2$, so $\ln\psi_1 = \ln\psi_2 + c$ and $\psi_1 = e^c\psi_2$. Two solutions of a second-order linear ODE that are proportional on an interval agree everywhere (uniqueness of the initial-value problem at any point of the interval), so $\psi_1 = \text{const}\cdot\psi_2$ on the whole line. ∎

■ **What just happened.** Step 3 is the only place "bound" was used. Scattering states *are* degenerate ($e^{ikx}$ and $e^{-ikx}$ both have energy $\hbar^2k^2/2m$; their Wronskian is $-2ik \ne 0$), and so are states on a ring, where periodic boundary conditions replace decay (Exercise B4). In 3D, degeneracy is the rule (QM-5, QM-6).

**Theorem 3.3 (Node theorem, stated).** Order the bound states by energy, $E_1 < E_2 < \cdots$. Then $\psi_n$ has exactly $n - 1$ nodes. This is Sturm's oscillation theorem (1836) for Sturm–Liouville problems; Courant and Hilbert, *Methods of Mathematical Physics*, vol. I, ch. VI, prove it. The project checks it for the box. The ground state has no nodes, so it can be chosen positive.

**Theorem 3.4 (Parity).** If $V(-x) = V(x)$ and $\psi(x)$ solves $\hat H\psi = E\psi$, so does $\psi(-x)$. Bound states can then be chosen with definite parity, and in fact each one is either even or odd.

*Proof.* 1. Let $\chi(x) = \psi(-x)$. Then $\chi''(x) = \psi''(-x)$, and $-\frac{\hbar^2}{2m}\chi''(x) + V(x)\chi(x) = -\frac{\hbar^2}{2m}\psi''(-x) + V(-x)\psi(-x) = E\psi(-x) = E\chi(x)$.
2. Then $\psi_\pm = \psi(x) \pm \psi(-x)$ are solutions (or zero) with definite parity, and $\psi = \frac12(\psi_+ + \psi_-)$, so the eigenspace has a basis of definite-parity functions.
3. For a bound state, Theorem 3.2 says $\chi = c\psi$. Applying $x\to -x$ twice: $\psi = c^2\psi$, so $c = \pm1$. ∎

Consequences: in any bound eigenstate of an even potential, $|\psi|^2$ is even and $\langle x\rangle = 0$. For the oscillator, $a_+$ is odd under $x \to -x$ (both $\hat x$ and $\hat p$ flip sign), so $\psi_n = a_+^n\psi_0/\sqrt{n!}$ has parity $(-1)^n$.

### 3.12 Periodic potentials in one page

A crystal is a periodic potential, $V(x + a) = V(x)$. **Bloch's theorem (stated):** the eigenstates can be chosen to satisfy
$$\psi(x + a) = e^{iKa}\psi(x) \tag{3.33}$$
for some real **Bloch wavenumber** $K$, defined modulo $2\pi/a$. (Why: translation by $a$ commutes with $\hat H$, so they share eigenstates, and a translation that preserves the norm can only multiply by a phase. QM-4's symmetry discussion makes this precise; Ashcroft and Mermin ch. 8 prove it.)

**The Dirac comb.** $V(x) = \alpha\sum_j\delta(x - ja)$, $\alpha > 0$ (repulsive, as in Griffiths §5.3). In $0 < x < a$, $V = 0$ and
$$\psi = A\sin kx + B\cos kx, \qquad k = \sqrt{2mE}/\hbar. \tag{3.34}$$
By (3.33), in $-a < x < 0$, $\psi(x) = e^{-iKa}\big[A\sin k(x + a) + B\cos k(x + a)\big]$.

▶ **Try it first.** Impose continuity of $\psi$ at $x = 0$ and the jump condition $\Delta\psi' = (2m\alpha/\hbar^2)\psi(0)$; eliminate $A$ and $B$.

1. Continuity: $B = e^{-iKa}(A\sin ka + B\cos ka)$, so $A\sin ka = (e^{iKa} - \cos ka)B$.
2. Jump: $kA - e^{-iKa}k(A\cos ka - B\sin ka) = \frac{2m\alpha}{\hbar^2}B$.
3. Multiply step 2 by $\sin ka$ and use step 1 for $A\sin ka$: $k(e^{iKa} - \cos ka)B - e^{-iKa}k\big[(e^{iKa} - \cos ka)\cos ka - \sin^2ka\big]B = \frac{2m\alpha}{\hbar^2}\sin(ka)\,B$.
4. The bracket is $e^{iKa}\cos ka - 1$, so the left side is $k\big[e^{iKa} - \cos ka - \cos ka + e^{-iKa}\big]B = 2k(\cos Ka - \cos ka)B$.
5. Cancel $2kB$:
$$\cos Ka = \cos ka + \frac{m\alpha}{\hbar^2k}\sin ka = \cos ka + \beta\,\frac{\sin ka}{ka}, \qquad \beta = \frac{m\alpha a}{\hbar^2}. \tag{3.35}$$

■ **What just happened.** For real $K$, the left side lies in $[-1, 1]$. The right side $f(E)$ does not always: wherever $|f(E)| > 1$, there is no real $K$, no Bloch state, and no allowed energy. The spectrum breaks into **bands** (where $|f| \le 1$) separated by **gaps**.

**The bands.** At $E\to0$, $f\to1 + \beta > 1$: a gap at the bottom. Each band ends where $ka = n\pi$, $f = (-1)^n$; for $\beta = 3$, $a = 1$ (project units) the first band runs from $E = 1.954$ to $4.936 = \pi^2/2$, the second from $9.463$ to $19.74 = 4\pi^2/2$. As $E$ grows, $\beta\sin(ka)/ka\to0$ and the gaps shrink. A finite comb of $N$ deltas transmits in windows matching the bands and is opaque in the gaps, where $T$ decays like $e^{-2N\operatorname{Im}(K)a}$; `transmission.png` shows $N = 12$.

**Conductors and insulators.** Each band holds $2N$ electrons for $N$ cells (two spins, $N$ values of $K$). Filled bands cannot respond to a field without lifting electrons across a gap (insulator; semiconductor if the gap is small); a partly filled band can (metal). SM-5's free electron gas is the $\beta\to0$ limit.

## 4. Worked examples

### Example 4.1: The box with $\Psi(x, 0) = Ax(a - x)$

*Normalize.* $\int_0^ax^2(a - x)^2dx = a^5/30$, so $A = \sqrt{30/a^5}$.

*Coefficients.* By (3.5), $c_n = \sqrt{\frac2a}\sqrt{\frac{30}{a^5}}\int_0^ax(a - x)\sin\frac{n\pi x}{a}dx$. Substitute $u = x/a$: the integral is $a^3\int_0^1u(1 - u)\sin(n\pi u)du$. Integrating by parts twice (the boundary terms vanish because $u(1-u)$ is zero at both ends, and the second-derivative term is $-2$):
$$\int_0^1u(1 - u)\sin(n\pi u)\,du = \frac{2\big(1 - (-1)^n\big)}{(n\pi)^3}.$$
So $c_n = \frac{\sqrt{60}}{a^3}\cdot a^3\cdot\frac{4}{(n\pi)^3} = \frac{8\sqrt{15}}{(n\pi)^3}$ for odd $n$, and $c_n = 0$ for even $n$ (the initial state is even about $a/2$; the even-$n$ states are odd about $a/2$).

*Probabilities.* $c_1 = 8\sqrt{15}/\pi^3 = 0.999277$, so $P(E_1) = |c_1|^2 = 960/\pi^6 = 0.998555$. $P(E_3) = 960/(729\pi^6) = 0.001370$. The parabola is almost exactly the ground state.

*Energy.* $\sum|c_n|^2E_n = \sum_{n\text{ odd}}\frac{960}{n^6\pi^6}\cdot\frac{n^2\pi^2\hbar^2}{2ma^2} = \frac{480\hbar^2}{\pi^4ma^2}\sum_{n\text{ odd}}\frac{1}{n^4} = \frac{480\hbar^2}{\pi^4ma^2}\cdot\frac{\pi^4}{96} = \frac{5\hbar^2}{ma^2}$. Directly: $\langle H\rangle = -\frac{\hbar^2}{2m}\int\Psi\Psi''dx = -\frac{\hbar^2}{2m}\cdot\frac{30}{a^5}\int_0^ax(a - x)(-2)\,dx = \frac{\hbar^2}{m}\cdot\frac{30}{a^5}\cdot\frac{a^3}{6} = \frac{5\hbar^2}{ma^2}$. Compare $E_1 = \pi^2\hbar^2/2ma^2 = 4.9348\,\hbar^2/ma^2$: slightly above $E_1$, as it must be. The project sums 2000 odd terms and gets $5.000000$.

### Example 4.2: $\langle x^2\rangle$, $\langle p^2\rangle$, and $\sigma_x\sigma_p$ in state $\psi_n$ by ladders

Square (3.14): $\hat x^2 = \frac{\hbar}{2m\omega}(a_+^2 + a_+a_- + a_-a_+ + a_-^2)$. In $\langle\psi_n|\cdot|\psi_n\rangle$, $a_\pm^2\psi_n\propto\psi_{n\pm2}$ is orthogonal to $\psi_n$. By (3.13), $a_+a_-\psi_n = n\psi_n$ and $a_-a_+\psi_n = (n+1)\psi_n$. So
$$\langle x^2\rangle_n = \frac{\hbar}{2m\omega}(2n + 1) = \Big(n + \frac12\Big)\frac{\hbar}{m\omega}.$$
With the signs of (3.14): $\hat p^2 = -\frac{\hbar m\omega}{2}(a_+^2 - a_+a_- - a_-a_+ + a_-^2)$, giving $\langle p^2\rangle_n = \frac{\hbar m\omega}{2}(2n+1) = (n + \frac12)\hbar m\omega$.

*Virial.* $\langle V\rangle = \frac12m\omega^2\langle x^2\rangle = \frac12(n + \frac12)\hbar\omega$ and $\langle T\rangle = \langle p^2\rangle/2m = \frac12(n + \frac12)\hbar\omega$. Equal, each half of $E_n$.

*Uncertainty.* $\langle x\rangle = \langle p\rangle = 0$ (§3.2), so $\sigma_x\sigma_p = \sqrt{\langle x^2\rangle\langle p^2\rangle} = (n + \frac12)\hbar$: $\frac12\hbar, \frac32\hbar, \frac52\hbar, \frac72\hbar$ for $n = 0, 1, 2, 3$. Only the ground state saturates QM-2's bound. The project computes these from the $60\times60$ ladder matrices.

### Example 4.3: Termination of the Hermite series at $n = 2$

Take $K = 5$ ($E = \frac52\hbar\omega$) and the even series, $a_1 = 0$. From (3.16): $a_2 = \frac{0 + 1 - 5}{1\cdot2}a_0 = -2a_0$; $a_4 = \frac{2\cdot2 + 1 - 5}{3\cdot4}a_2 = 0$; then $a_6 = a_8 = \dots = 0$. So $h = a_0(1 - 2\xi^2)$. Check in (3.15) with $K - 1 = 4$: $h'' = -4a_0$, $-2\xi h' = 8a_0\xi^2$, $4h = 4a_0 - 8a_0\xi^2$; sum $0$. With the Hermite normalization, $h = -\frac{a_0}{2}H_2$, $H_2 = 4\xi^2 - 2$, and
$$\psi_2 = \Big(\frac{m\omega}{\pi\hbar}\Big)^{1/4}\frac{1}{\sqrt8}(4\xi^2 - 2)e^{-\xi^2/2}.$$
*Cross-check.* With $a_+ = (\xi - d/d\xi)/\sqrt2$ and $\psi_1 \propto \sqrt2\,\xi e^{-\xi^2/2}$, $a_+\psi_1 \propto (2\xi^2 - 1)e^{-\xi^2/2} = \sqrt2\psi_2$ (same prefactor), as (3.13) requires. If instead $K = 5$ but $a_1 \ne 0$, the odd series has $a_{j+2}/a_j = (2j - 4)/((j+1)(j+2))$, never zero for odd $j$, so it diverges like $e^{\xi^2}$; that is why $a_1 = 0$ is forced.

### Example 4.4: The delta well: bound state, and transmission at $E = |E_b|$

Bound state (3.20): $\kappa = m\alpha/\hbar^2$, $E_b = -m\alpha^2/2\hbar^2$, $\psi = \sqrt\kappa\,e^{-\kappa|x|}$. Now scatter at $E = |E_b| = m\alpha^2/2\hbar^2$. Then $k = \sqrt{2mE}/\hbar = m\alpha/\hbar^2$, so $\beta = m\alpha/\hbar^2k = 1$. By (3.21), $T = R = \frac12$; equivalently $T = 1/(1 + m\alpha^2/2\hbar^2E) = 1/(1 + 1)$. Half the particles reflect from a well whose only bound state has the same $|E|$. The amplitudes: $F/A = 1/(1 - i)$, $|F/A| = 1/\sqrt2$, phase $+45°$; $B/A = i/(1 - i)$. Flip the sign of $\alpha$ and $T$ is still $\frac12$ (project: $0.5000$ both).

### Example 4.5: The finite well with $z_0 = 8$

$\lfloor 2\cdot8/\pi\rfloor + 1 = \lfloor5.09\rfloor + 1 = 6$ bound states. Solving (3.22) numerically (the project uses Brent's method on $z\sin z - \sqrt{z_0^2 - z^2}\cos z$ and $z\cos z + \sqrt{z_0^2 - z^2}\sin z$, which have the same roots as (3.22) but no poles):

| state | parity | $z_n$ | $n\pi/2$ (infinite well) | $E_n/V_0 = z_n^2/z_0^2 - 1$ |
|---|---|---|---|---|
| 1 | even | 1.3955 | 1.5708 | −0.9696 |
| 2 | odd | 2.7859 | 3.1416 | −0.8787 |
| 3 | even | 4.1648 | 4.7124 | −0.7290 |
| 4 | odd | 5.5214 | 6.2832 | −0.5237 |
| 5 | even | 6.8307 | 7.8540 | −0.2710 |
| 6 | odd | 7.9573 | — | −0.0106 |

The infinite-well estimate $z \approx n\pi/2$ is 12.6% high for $n = 1$ and worse higher up, where states leak further. The sixth state, near $z_0$, has $\kappa a = \sqrt{64 - 7.9573^2} = 0.83$: its decay length $a/0.83$ is longer than the half-width. With $\hbar = m = a = 1$, $V_0 = 32$, the energies are $E_n = 32\,(E_n/V_0) = -31.026, -28.119, -23.327, -16.757, -8.671, -0.341$, and the finite-difference solver reproduces them to a relative $10^{-4}$ in $E + V_0$.

### Example 4.6: An electron tunneling through a 5 eV-deep forbidden region

$E = 5$ eV, $V_0 = 10$ eV, $a = 0.5$ nm. Then $\kappa = \sqrt{2m_e(5\text{ eV})}/\hbar = 1.1456\times10^{10}$ m$^{-1}$ and $\kappa a = 5.728$. The prefactor $V_0^2/4E(V_0 - E) = 100/100 = 1$, so (3.26) is $T = 1/(1 + \sinh^2\kappa a) = 1/\cosh^2\kappa a$:
$$T(0.5\text{ nm}) = 4.2353\times10^{-5}.$$
The thick-barrier formula (3.27): $16\cdot\frac12\cdot\frac12 e^{-11.456} = 4e^{-11.456} = 4.2354\times10^{-5}$, indistinguishable. Now double the width to 1 nm: $\kappa a = 11.456$, $T = 4.4846\times10^{-10}$, smaller by a factor $9.44\times10^4 = e^{11.456} = e^{2\kappa\,(0.5\text{ nm})}$: the extra 0.5 nm costs one more factor $e^{-2\kappa(0.5\text{ nm})}$, and the prefactor 4 is common to both. Doubling the width squares $T/4$.

### Example 4.7: The step with $E = 2V_0$

$k' = \sqrt{2m(E - V_0)}/\hbar = k/\sqrt2$. With $r = k'/k = 1/\sqrt2$, (3.23) gives
$$R = \Big(\frac{1 - r}{1 + r}\Big)^2 = \Big(\frac{\sqrt2 - 1}{\sqrt2 + 1}\Big)^2 = (\sqrt2 - 1)^4 = 17 - 12\sqrt2 = 0.0294,$$
$$T = \frac{4r}{(1 + r)^2} = \frac{4\sqrt2}{(\sqrt2 + 1)^2} = \frac{4\sqrt2}{3 + 2\sqrt2} = 4\sqrt2(3 - 2\sqrt2) = 12\sqrt2 - 16 = 0.9706.$$
$R + T = 1$. The naive $|F/A|^2 = 4/(1 + r)^2 = 8/(3 + 2\sqrt2) = 24 - 16\sqrt2 = 1.3726$: "137% transmission," the error the current definition prevents.

### Example 4.8: Alpha decay (nuclear physics)

*Alpha decay.* Use (3.30) with the alpha mass $m = 6.645\times10^{-27}$ kg and the touching radius $r_1 = 1.2\text{ fm}\,(A_d^{1/3} + 4^{1/3})$ of daughter (mass number $A_d$) and alpha. This is a model choice, and the absolute half-lives move by a factor of several if you change the 1.2 fm.

- $^{238}$U → $^{234}$Th: $E = 4.27$ MeV, $Z = 90$. $\gamma = 44.4$, $t_{1/2} = 9.87\times10^9$ yr (measured $4.47\times10^9$ yr).
- $^{212}$Po → $^{208}$Pb: $E = 8.95$ MeV, $Z = 82$. $\gamma = 16.3$, $t_{1/2} = 8.85\times10^{-8}$ s (measured $2.99\times10^{-7}$ s).

The energies differ by a factor of $2.1$; the predicted half-lives by $10^{24.5}$, the measured ones by $10^{23.7}$. A 1D tunneling estimate with no nuclear structure gets both within a factor of about 3, because nearly all the variation is in $e^{2\gamma}$. The project's `gamow.png` plots $\log_{10}t_{1/2}$ against $E$.

## 5. Why it matters and where it connects

- **QM-4**: $a_\pm$ become $a$, $a^\dagger$, and (3.13) reads $a^\dagger|n\rangle = \sqrt{n+1}|n+1\rangle$; the double well becomes the two-level system $E_0 - \frac\Delta2\sigma_x$; (3.33) becomes "translations commute with $\hat H$."
- **QM-5** reuses the ladder trick for $L_\pm$; the norm computation of (3.13) gives $\sqrt{\ell(\ell+1) - m(m\pm1)}$. **QM-6** quantizes hydrogen by the asymptotic-factor-plus-termination argument of §3.3.
- **QM-8**: WKB derives (3.29); perturbation theory starts from these spectra; the variational method generalizes Theorem 3.1.
- **SM-5**: (3.2) in 3D is the Fermi gas; the $\hbar\omega$ ladder is the phonon (Debye) and the photon (blackbody); §3.12 refines the free-electron metal. **CM-6**: each normal mode becomes an independent ladder.
- **LA-8, LA-13**: the FD Hamiltonian is real symmetric tridiagonal (real spectrum, orthogonal eigenvectors) with $O(h^2)$ discretization error.
- **EM-3, F2, EM-7**: §3.3's termination is EM-3 §3.5's Legendre argument; (3.5) is a Fourier sine series; the step's $R$ and $T$ are normal-incidence Fresnel coefficients with $k$ as refractive index, and the double barrier is a Fabry–Pérot etalon.

## 6. Exercises

See [exercises.md](exercises.md). **A1–A5** are drills (a 1 nm box, orthogonality, $a_+\psi_0$, a step, the delta well). **B1–B9** are the module: ladder matrix elements (B1), a discontinuous state with infinite $\langle H\rangle$ (B2), oscillator sloshing (B3), where nondegeneracy fails (B4), a real finite well (B5), the barrier above its top (B6), $R + T = 1$ in general (B7), the step in `transfer()` (B8, code), parity (B9). Essential: B1, B2, B6, B7, B8. **C1–C3** go further: the two-delta double well's splitting, the triangular well and Airy functions (code), and the Gamow integral.

## 7. Coding project

`project.py` (about 4 s) checks every number in §3–§4 and in the solutions. (1) **Finite-difference eigensolver** `fd_eigs`: the tridiagonal $H_{jj} = \hbar^2/mh^2 + V_j$, $H_{j,j\pm1} = -\hbar^2/2mh^2$, solved with `scipy.linalg.eigh_tridiagonal`; box and oscillator spectra and the node theorem. (2) Example 4.1 and the revival (3.6). (3) **Oscillator**: Hermite recurrence (3.17) vs `scipy.special.hermite`, Example 4.3 in SymPy, $60\times60$ ladder matrices ($[a_-, a_+] = 1$ except a last diagonal entry $-(N-1)$, a truncation artifact you should be able to explain), Example 4.2, and a coherent state $|\alpha = 2\rangle$ that oscillates rigidly with $\sigma_x = 1/\sqrt2$ forever (`oscillator.png`; the state of a laser mode, previewing QC-9). (4) **Finite well**: Example 4.5 by root-finding and by FD; Theorem 3.1; the double-well splitting (`finite_well.png`). (5) **Transfer matrices** `transfer()` implementing (3.31): barrier, step, delta, well resonances, $\ln T$ slope, double barrier, Dirac comb (`transmission.png`). (6) SI numbers and Gamow (`gamow.png`). (7) **A packet hits a barrier** with QM-2's split-operator solver: transmitted probability $0.3483$ against the stationary $T(k)$ averaged over $|\phi(k)|^2$, $0.3482$ ($T(k_0) = 0.3470$): stationary scattering is the long-packet limit (`packet_barrier.png`). (8) Numbers for the solutions.

**Build it yourself first:** write `fd_eigs` (ten lines; `np.linalg.eigh` on a dense matrix is fine) and a single-barrier transfer matrix from (3.24)–(3.25), and check against (3.26) before reading the file's version. Then do the "Your turn" list: the half-oscillator (predict first which states survive), the triangular well with `scipy.special.airy`, the two-delta splitting, Kronig–Penney barriers, and a narrow packet.

## 8. Where you'll struggle

- **The ladder algebra feels like a trick.** It is the most important technique in the module. Redo §3.2 from a blank page until $\hat H = \hbar\omega(a_+a_- + \frac12)$, the ladder property and the normalization move ($|c|^2 = \langle\psi|a_\mp a_\pm\psi\rangle$, then (3.9)) take ten minutes without notes. QM-5's $L_\pm$ and all second quantization use that move.
- **Which order, $a_+a_-$ or $a_-a_+$?** You will put the $\frac12$ on the wrong side. Check with $\psi_0$: $a_-\psi_0 = 0$, so $\hat H\psi_0 = \frac12\hbar\omega\psi_0$ needs $\hat H = \hbar\omega(a_+a_- + \frac12)$.
- **The finite-well equation has no closed-form solution.** Students waste hours trying to solve (3.22) by algebra. Draw the graph, count crossings, use a root-finder. Analytically you can count states, take limits, and bound the roots.
- **"$T = |F/A|^2$" is wrong when the speeds differ.** It is correct only with equal speeds on both sides; for a step it gives 137% (Example 4.7). Define $T$ by currents, (3.18).
- **The barrier algebra is where arithmetic errors live.** A sign error in one exponent is invisible until $T > 1$; that is why §3.8 uses matrices. Check $T\to1$ as $V_0\to0$, $T\to0$ as $a\to\infty$, $R + T = 1$, and `transfer()`.
- **Forgetting that a delta makes $\psi'$ jump** (and that a finite step does not). Integrate across the point every time; (3.19) is two lines.
- **The Frobenius termination argument is the same argument you saw for Legendre (EM-3).** If "why must it terminate?" feels arbitrary, redo step 4 of §3.3: a non-terminating series behaves like $e^{\xi^2}$. Hydrogen (QM-6) needs exactly this.
- **Bloch's theorem is stated, not proved here.** It is QM-4's "commuting operators share eigenstates" applied to translations, and will feel unmotivated until then.
- **The double-well splitting is easy to state and hard to compute well.** The exponential $e^{-\kappa b}$ is robust; the prefactor is not. C1 gives one exact case; QM-8 the general one.

## 9. Check yourself

1. Derive $E_n$ for the infinite well in three lines. Why is $\psi'$ allowed to jump at the walls?
2. What is the revival time of the infinite well, and why do all states revive?
3. Starting from $\hat p = -i\hbar\,d/dx$, prove $[\hat x, \hat p] = i\hbar$.
4. Write $\hat H$ in terms of $a_\pm$ and prove $\hat H(a_+\psi) = (E + \hbar\omega)a_+\psi$.
5. Why must there be a state with $a_-\psi_0 = 0$? Solve for it.
6. Derive $a_+\psi_n = \sqrt{n+1}\psi_{n+1}$.
7. Why must the Hermite series terminate, and what does termination quantize?
8. Derive the jump condition for $\psi'$ at a delta and the delta well's bound-state energy.
9. Sketch the graphical solution of the finite well. How many bound states for $z_0 = 3$?
10. Why is $T \ne |F/A|^2$ for a step? Write $T$ in terms of currents.
11. Give the thick-barrier tunneling formula. If $\kappa a = 10$, by what factor does $T$ change when $a$ increases 10%?
12. State Bloch's theorem and the Dirac comb's band condition. Where are the band edges?

## 10. Pair with

- **Griffiths & Schroeter, *Introduction to Quantum Mechanics* (3rd ed.), ch. 2 and §5.3.** The standard treatment in these conventions; its delta, barrier and finite-well problems are the practice set.
- **MIT 8.04 (Zwiebach, 2016), lectures 7–16** (OCW). Better than these notes on the node theorem and the shooting picture of quantization.
- **Zwiebach, *Mastering Quantum Mechanics*, ch. 6–9.** More careful 1D theorems and a fuller treatment of phase shifts and time delays.
- **Shankar, *Principles of Quantum Mechanics*, ch. 5, 7.** The oscillator in the energy basis with $a$, $a^\dagger$ from the start, the version QM-4 uses.
- **Townsend, *A Modern Approach to Quantum Mechanics*, ch. 6–7.** Operator-language oscillator, and the double well and ammonia in more depth.
