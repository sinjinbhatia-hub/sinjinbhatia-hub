# QM-5 Exercises

Conventions as in the notes: $[\hat L_i, \hat L_j] = i\hbar\varepsilon_{ijk}\hat L_k$, $\hat J_\pm = \hat J_x \pm i\hat J_y$, matrices in the basis $m = j, j-1, \ldots, -j$, Condon–Shortley phases, $\hat{\mathbf S} = \frac\hbar2\boldsymbol\sigma$, $|\hat{\mathbf n}\rangle = \cos\frac\theta2|{\uparrow}\rangle + e^{i\phi}\sin\frac\theta2|{\downarrow}\rangle$, $\hat H = -\gamma\mathbf B\cdot\hat{\mathbf S}$.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Using only $[\hat x_i, \hat p_j] = i\hbar\delta_{ij}$, compute $[\hat L_z, \hat x]$, $[\hat L_z, \hat y]$, $[\hat L_z, \hat z]$ and $[\hat L_z, \hat p_x]$. Then show $[\hat L_z, \hat x^2 + \hat y^2] = 0$ and say in one sentence why that had to happen.

**A2.** For $\ell = 1$, write $\hat L_+$, $\hat L_-$, $\hat L_x$, $\hat L_y$, $\hat L_z$ as $3\times3$ matrices using (3.14). Verify $[\hat L_x, \hat L_y] = i\hbar\hat L_z$ by multiplying matrices, and find the eigenvalues of $\hat L_x$ without diagonalizing (one sentence).

**A3.** (a) Compute $\sigma_x\sigma_y\sigma_z$. (b) Show $\operatorname{tr}(\sigma_a\sigma_b\sigma_c) = 2i\varepsilon_{abc}$. (c) For $\mathbf a = (1, 2, 2)$, compute $(\mathbf a\cdot\boldsymbol\sigma)^2$ from (3.26) and deduce the eigenvalues of $\mathbf a\cdot\boldsymbol\sigma$.

**A4.** A spin-$\frac12$ particle is in $|{+y}\rangle = (|{\uparrow}\rangle + i|{\downarrow}\rangle)/\sqrt2$. (a) Find its $(\theta, \phi)$ and Bloch vector. (b) Give the probabilities of $\pm\hbar/2$ for $\hat S_x$, $\hat S_y$, $\hat S_z$. (c) Which state is orthogonal to it, and where is that on the sphere?

**A5.** (a) Compute the Larmor frequency $f = |\gamma|B/2\pi$ of a proton at 3 T and of an electron at 1 T, using $\gamma_p = 2.6752\times10^8$ and $|\gamma_e| = 1.7609\times10^{11}$ rad s$^{-1}$T$^{-1}$. (b) Which way does each precess about $\mathbf B$? (c) Why does electron spin resonance use microwaves and NMR radio waves?

**A6.** Use (3.42) to decompose (a) $2\otimes1$ and (b) $\frac12\otimes\frac12\otimes\frac12$ into multiplets, and check the dimension count in each case.

## B. Core (do all)

**B1.** (Central potentials.) (a) Show $[\hat L_i, \hat r^2] = 0$ and $[\hat L_i, \hat p^2] = 0$, where $\hat r^2 = \hat x_j\hat x_j$ and $\hat p^2 = \hat p_j\hat p_j$. (b) Conclude $[\hat{\mathbf L}, \hat H] = 0$ for $\hat H = \hat p^2/2m + V(\hat r)$, and say what this implies for the energy levels.

**B2.** (Uncertainty in $|\ell, m\rangle$.) (a) Show $\langle\hat L_x\rangle = \langle\hat L_y\rangle = 0$ in $|\ell, m\rangle$. (b) Show $\langle\hat L_x^2\rangle = \langle\hat L_y^2\rangle = \frac{\hbar^2}{2}\big(\ell(\ell+1) - m^2\big)$. (c) Check the uncertainty relation $\sigma_{L_x}\sigma_{L_y}\ge\frac\hbar2|\langle\hat L_z\rangle|$ and find when it is an equality. (d) In the "vector model," $\mathbf L$ lies on a cone about $z$ of half-angle $\alpha$ with $\cos\alpha = m/\sqrt{\ell(\ell+1)}$. Compute the smallest possible $\alpha$ for $\ell = 2$. Why can't $\alpha$ be zero?

**B3.** Derive $Y_2^1$ from $Y_2^2 = \sqrt{\frac{15}{32\pi}}\sin^2\theta\,e^{2i\phi}$ by applying $\hat L_-$ (3.18) and dividing by the coefficient from (3.14). Check against (3.22), and check $\hat L_+Y_2^2 = 0$.

**B4.** (A general spin state.) $|\chi\rangle = A(3i, 4)^\top$.
(a) Normalize it and compute $\langle\hat S_x\rangle$, $\langle\hat S_y\rangle$, $\langle\hat S_z\rangle$.
(b) Find its Bloch vector and $(\theta,\phi)$. Find the probabilities of $+\hbar/2$ for $\hat S_x$ and for $\hat S_y$.
(c) Prove that $|\hat{\mathbf n}\rangle\langle\hat{\mathbf n}| = \frac12(I + \hat{\mathbf n}\cdot\boldsymbol\sigma)$ and hence $|\langle\hat{\mathbf n}'|\hat{\mathbf n}\rangle|^2 = \frac12(1 + \hat{\mathbf n}\cdot\hat{\mathbf n}')$.

**B5.** (Gates as rotations.) (a) Show $R_x(\pi) = -i\sigma_x$ and $R_y(\pi/2)|{\uparrow}\rangle = |{+x}\rangle$. (b) Show the Hadamard $H = (\sigma_x + \sigma_z)/\sqrt2$ equals $iR_{\hat{\mathbf n}}(\pi)$ with $\hat{\mathbf n} = (\hat{\mathbf x} + \hat{\mathbf z})/\sqrt2$. (c) Use Theorem 3.4 to show $H\sigma_zH = \sigma_x$ without multiplying matrices.

**B6.** (Rabi line shape.) (a) At detuning $\Delta = \omega_1$, find the maximum flip probability and the first time it is reached. (b) Show that the time average of (3.35) is $\frac12\,\omega_1^2/(\omega_1^2 + \Delta^2)$ and find its full width at half maximum in $\Delta$. (c) For protons with $B_1 = 10\ \mu$T, give $\omega_1/2\pi$, the $\pi$-pulse length, and the detuning in Hz at which the maximum flip probability drops to $\frac12$.

**B7.** (Spin-1 Stern–Gerlach.) (a) Derive $\hat S_x$ for spin 1 from (3.14) and find its eigenvectors. (b) Atoms in $|1, 1\rangle$ enter an SG$x$ magnet. Give the three beam probabilities. (c) The $m_x = 0$ beam alone enters an SG$z$ magnet. Give the three probabilities and explain in words why the middle beam is empty. (d) Run Part 8 of `project.py` and compare.

**B8.** (The singlet.) (a) Show that each of $\hat S_x$, $\hat S_y$, $\hat S_z$ (total) annihilates $|0,0\rangle = \frac{1}{\sqrt2}(|{\uparrow\downarrow}\rangle - |{\downarrow\uparrow}\rangle)$. (b) Show $\boldsymbol\sigma_1\cdot\boldsymbol\sigma_2 = 2\hat S^2/\hbar^2 - 3$ and evaluate it on the singlet and on the triplet. (c) Write the singlet in the $x$ basis $|{\pm x}\rangle$ and comment.

## C. Challenge (do at least one)

**C1.** (Why orbital $\ell$ can't be $\frac12$.) Suppose $f(\theta,\phi) = \sqrt{\sin\theta}\,e^{i\phi/2}$ were "$Y_{1/2}^{1/2}$." (a) Show $\hat L_+f = 0$ and $\hat L_zf = \frac\hbar2f$. (b) Compute $g = \hat L_-f$; it should be the bottom rung $m = -\frac12$. (c) Show $\hat L_-g\ne0$, so the ladder fails to terminate, contradicting Theorem 3.1. (d) What else is wrong with $f$ as a wavefunction?

**C2.** (Clebsch–Gordan by code.) Write a function that, given $j_1$, $j_2$, builds the total $\hat J$ with `np.kron`, starts from the top state, lowers with $\hat J_-$, and for each smaller $j$ finds the top state as the vector in the $m = j$ subspace orthogonal to all previously found states, with the Condon–Shortley sign. Finish Example 4.7 by hand (the three $m < 0$ coupled states of $1\otimes\frac12$) and check them with your code; then check $\frac32\otimes1$ against `sympy.physics.quantum.cg.CG`. Compare with `multiplets_by_lowering` in `project.py`.

**C3.** (The 2-to-1 map.) For $U\in SU(2)$ define $R_{ab} = \frac12\operatorname{tr}(\sigma_aU\sigma_bU^\dagger)$. (a) Show $U\sigma_bU^\dagger = R_{ab}\sigma_a$ and that $R$ is real and orthogonal. (b) Show $R(-U) = R(U)$, and that $R(U) = I$ only if $U = \pm I$. (c) Compute $R$ for $U = R_y(\pi/2)$ and identify it.
