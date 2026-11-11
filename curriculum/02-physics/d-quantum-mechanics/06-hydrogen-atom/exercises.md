# QM-6 Exercises

Conventions as in the notes: $\psi_{n\ell m} = R_{n\ell}(r)Y_\ell^m(\theta,\phi)$, $u = rR$, $\int_0^\infty u^2dr = 1$, $n = n_r + \ell + 1$, $a = a_0 = 0.052918$ nm, $E_n = -13.606\ \text{eV}/n^2$ (use $\mu$ when the problem says so), modern Laguerre polynomials (SciPy `genlaguerre`), Condon–Shortley $Y_\ell^m$. Atomic units ($\hbar = m_e = e^2/4\pi\varepsilon_0 = 1$) are allowed wherever SI is not asked for.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** List every state $|n\ell m\rangle$ with $n = 4$ in a table with columns $\ell$, spectroscopic name, number of radial nodes, allowed $m$. Count them with and without spin and check against (3.14).

**A2.** (a) Show that $R_{21} = \frac{1}{\sqrt{24}}a^{-3/2}\frac ra e^{-r/2a}$ is normalized. (b) Find the most probable radius for the $2p$ state. (c) Without integrating, explain why the most probable radius for $3d$ is $9a$.

**A3.** Compute the vacuum wavelengths, using $R_{\rm H} = 1.09678\times10^7$ m$^{-1}$, of Lyman-$\alpha$ ($2\to1$), H$\beta$ ($4\to2$), Paschen-$\alpha$ ($4\to3$), and the Balmer series limit. Which of these can your eye see?

**A4.** For $\text{He}^+$, use (3.21) to give the ground-state energy (with the reduced mass; $m_\alpha = 7294.3\,m_e$), the Bohr radius, the ionization energy, and the wavelength of the $2\to1$ line. Why is ionizing $\text{He}^+$ harder than ionizing hydrogen by a factor of 4 and not 2?

**A5.** (a) Compute $|\psi_{100}(0)|^2$ in m$^{-3}$. (b) Show that $\psi_{n\ell m}(0) = 0$ for every $\ell\ge1$ and explain in one sentence why only s-states can be found at the nucleus. (c) Is $|\psi|^2$ or $r^2|R|^2$ largest at $r = 0$ for $1s$?

**A6.** (a) Show that $d_{xz}\equiv(Y_2^{-1} - Y_2^1)/\sqrt2 = \sqrt{15/4\pi}\,xz/r^2$, using $Y_2^{\pm1} = \mp\sqrt{15/8\pi}\sin\theta\cos\theta\,e^{\pm i\phi}$. (b) If you measure $\hat L_z$ in the state $R_{32}d_{xz}$, what values can you get and with what probabilities? (c) Is $R_{32}d_{xz}$ an energy eigenstate?

## B. Core (do all)

**B1.** (The machinery once more.) (a) Starting from (3.8), substitute $u = \rho^{\ell+1}e^{-\rho}v$ and derive (3.10) yourself, showing every term. (b) Run the recursion (3.11) for $n = 3$, $\ell = 0$ and show $v = c_0(1 - 2\rho + \frac23\rho^2)$. (c) Convert to $r$ ($\rho = r/3a$), normalize, and check against $R_{30}$ in (3.16). (d) Find the two radial nodes in units of $a$.

**B2.** (a) Compute $\langle r\rangle$ and $\langle r^2\rangle$ in the $2p$ state by direct integration and check against (3.17). (b) Find $\sigma_r$ and $\sigma_r/\langle r\rangle$. (c) Compute $\langle r\rangle$ for $2s$ from (3.17) and explain in words why it is larger than $2p$'s even though the $2s$ has probability near the nucleus.

**B3.** (Virial theorem.) (a) Prove $2\langle T\rangle = \langle\mathbf r\cdot\nabla V\rangle$ in any normalizable stationary state (repeat the proof of Theorem 3.2 with every commutator written out). (b) For the 3D oscillator $V = \frac12m\omega^2r^2$, find $\langle T\rangle$ and $\langle V\rangle$ in a state of energy $E$. (c) For hydrogen in the state $n = 3$, give $\langle T\rangle$, $\langle V\rangle$, and $\langle1/r\rangle$ in eV and $a^{-1}$, without any integral.

**B4.** (Infinite spherical well.) (a) Solve the $\ell = 0$ case completely: $u$, the normalized $\psi_{n00}$, $E_{n0}$. (b) Show $j_1(x) = \sin x/x^2 - \cos x/x$ satisfies the $\ell = 1$ radial equation (3.2) with $V = 0$, and that its zeros satisfy $\tan x = x$. (c) Using the zero $\beta_{11} = 4.4934$, give $E_{11}/E_{10}$, and for an electron in a well of radius $a = 1$ nm give $E_{10}$, $E_{11}$ and $E_{20}$ in eV. (d) Explain why the levels depend on $\ell$ at fixed "number of nodes plus $\ell$," unlike hydrogen.

**B5.** (Muonic hydrogen.) A muon ($m_\mu = 206.768\,m_e$) binds to a proton ($m_p = 1836.15\,m_e$). (a) Compute the reduced mass, $E_1$, and the Bohr radius. (b) Treat the proton as a sphere of radius $r_p = 0.84$ fm and compute $P(r < r_p)$ in the $1s$ state for ordinary and muonic hydrogen, using the small-$x$ expansion of Example 4.6. (c) Show the ratio is close to $(\mu_\mu/\mu_e)^3$ and explain why this makes muonic hydrogen the sensitive probe of the proton radius.

**B6.** (Isotope shift.) In 1932 Harold Urey and collaborators found deuterium by observing faint companion lines displaced from the hydrogen Balmer lines. (a) Using $R_M = R_\infty/(1 + m_e/M)$ with $m_d = 3670.48\,m_e$, compute the Balmer-$\alpha$ wavelengths for H and D. (b) Show that to first order in $m_e/M$, $\lambda_{\rm H} - \lambda_{\rm D}\approx\lambda_\infty\,m_e(1/m_p - 1/m_d)$ and evaluate it. (c) On which side of the H line does the D line sit?

**B7.** (Correspondence.) For the circular state $\ell = n - 1$: (a) derive $\langle r\rangle$ and $\sigma_r$ from $P(r)\propto r^{2n}e^{-2r/na}$ using the Gamma integral, (b) evaluate $\sigma_r/\langle r\rangle$ for $n = 1, 10, 100$, (c) compare $\langle L^2\rangle$ with Bohr's $L^2 = n^2\hbar^2$ and compute the relative difference at $n = 100$. (d) Modify Part 5 of `project.py` to solve $\ell = 19$ on a suitable grid and assert that the peak of $|u|^2$ is at $400a$.

## C. Challenge (do at least one)

**C1.** ($SO(4)$ counting.) Accept from §3.12 that on a bound level the operators $\hat{\mathbf J}_{1,2} = \frac12(\hat{\mathbf L}\pm\hat{\mathbf M})$ are two commuting angular momenta with $\hat J_1^2 = \hat J_2^2$, and that $\hat{\mathbf L} = \hat{\mathbf J}_1 + \hat{\mathbf J}_2$. (a) If $j_1 = j_2 = j$, show using QM-5's addition rule (3.42) that the level contains $\ell = 0, 1, \ldots, 2j$, and that $n = 2j + 1$ reproduces the allowed $\ell$ of §3.6. (b) Count the dimension two ways. (c) Explain why half-integer $j$ is allowed here although $\ell$ must be an integer.

**C2.** (Feynman–Hellmann.) The Feynman–Hellmann theorem says $\partial E_n/\partial\lambda = \langle\partial\hat H/\partial\lambda\rangle$ for a normalized eigenstate of $\hat H(\lambda)$. (a) Prove it in two lines. (b) Treat $e^2$ as the parameter in the hydrogen radial Hamiltonian and derive $\langle1/r\rangle = 1/(n^2a)$ for every $n, \ell$. (c) Treat $\ell$ as a continuous parameter in the radial Hamiltonian (3.3), remembering that $E$ depends on $\ell$ through $n = n_r + \ell + 1$, and derive $\langle1/r^2\rangle = 1/[n^3(\ell + \frac12)a^2]$. (d) Check (c) for $1s$ by direct integration.

**C3.** (Why potassium's 19th electron is $4s$.) Model the potential seen by potassium's outer electron by $V(r) = -\frac{1}{r}\big[1 + (Z-1)e^{-r/d}\big]$ (atomic units, $Z = 19$): the full nuclear charge near $r = 0$, charge 1 far away. (a) Check both limits. (b) With `radial_solve`, compute the fourth $\ell = 0$ eigenvalue (the $4s$; $1s$, $2s$, $3s$ are occupied core states) and the lowest $\ell = 2$ eigenvalue (the $3d$) for $d = 0.5\,a_0$, in eV. Compare with potassium's measured $4s$ and $3d$ binding energies, $4.34$ and $1.67$ eV. (c) Compare with the unscreened hydrogenic values and explain the reversal using the inner humps of §3.7. (d) Find the $d$ at which the model's ordering flips, and say what physical claim the model is (and is not) making.
