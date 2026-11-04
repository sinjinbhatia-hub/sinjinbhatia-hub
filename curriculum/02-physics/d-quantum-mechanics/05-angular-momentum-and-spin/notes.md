# QM-5: Angular Momentum and Spin: the Commutator Algebra, Ladder Operators, Spherical Harmonics, Spin-1/2, Pauli Matrices, Stern–Gerlach, the Bloch Sphere

**Time:** ~14 hours · **Prerequisites:** QM-4 · **Unlocks:** QM-6, QM-7, QC-1

> **Conventions in this module.** $\hbar$ explicit (`project.py` sets $\hbar = 1$ outside its SI part). $\varepsilon_{ijk}$ is Levi-Civita ($\varepsilon_{xyz} = +1$), repeated indices summed. $\hat{\mathbf L} = \hat{\mathbf r}\times\hat{\mathbf p}$, $[\hat x_i, \hat p_j] = i\hbar\delta_{ij}$. A generic angular momentum $\hat{\mathbf J}$ is any Hermitian triple with $[\hat J_i, \hat J_j] = i\hbar\varepsilon_{ijk}\hat J_k$; $\hat{\mathbf L}$ and $\hat{\mathbf S}$ are cases. $\hat J_\pm = \hat J_x \pm i\hat J_y$. $\hat J^2|j,m\rangle = \hbar^2j(j+1)|j,m\rangle$, $\hat J_z|j,m\rangle = \hbar m|j,m\rangle$; matrices list $m$ from $+j$ down to $-j$. $\theta$ is polar from $+z$, $\phi$ azimuthal. $Y_\ell^m$ carry the **Condon–Shortley phase** $(-1)^m$ for $m > 0$, as in Griffiths, Jackson, SciPy's `sph_harm_y` and EM-3 §3.8 (no sign difference from EM-3). Spin-$\frac12$: $\hat{\mathbf S} = \frac\hbar2\boldsymbol\sigma$, $|{\uparrow}\rangle = (1,0)^\top$, $|{\downarrow}\rangle = (0,1)^\top$ (QC-1's $|0\rangle$, $|1\rangle$). Bloch parametrization $|\hat{\mathbf n}\rangle = \cos\frac\theta2|{\uparrow}\rangle + e^{i\phi}\sin\frac\theta2|{\downarrow}\rangle$, with **half** angles. $\boldsymbol\mu = \gamma\hat{\mathbf S}$, $\hat H = -\boldsymbol\mu\cdot\mathbf B$; $\gamma > 0$ for the proton, $< 0$ for the electron.

## 0. Where this sits

QM-4 gave you the formalism: Hilbert space, Hermitian observables, compatible observables sharing an eigenbasis, Heisenberg's $d\hat A/dt = \frac i\hbar[\hat H, \hat A]$. QM-1 §3.8 wrote down the spin-$\frac12$ matrices because they fit the Stern–Gerlach data. This module derives them, and EM-3 Theorem 3.6's spherical harmonics, from one commutation relation. You need QM-4, QM-3's oscillator ladder, EM-3's Legendre functions, and LA-12's Kronecker product (§3.14) and Appendix A on $SO(3)$ and $SU(2)$.

QM-6 separates hydrogen as $R(r)Y_\ell^m$ using §3.3. QM-7 starts from §3.10's singlet and §3.11's spin-statistics statement. QM-8 does spin–orbit coupling and the Zeeman effect in §3.10's coupled basis. QC-1 is §3.5–3.8 renamed.

## 1. The problem that forced this idea

In 1896 Pieter Zeeman saw sodium's spectral lines broaden in a magnetic field. Lorentz explained the simplest case at once: an orbiting electron's magnetic moment shifts its energy according to the orbit's orientation, and one line becomes three. Most lines split into four, six, or more. By the early 1920s Landé could fit these "anomalous" Zeeman patterns empirically, but only with half-integer quantum numbers, which made no sense for an orbit.

In 1922 Otto Stern and Walther Gerlach sent silver atoms through an inhomogeneous magnet and saw the beam split in two (QM-1 §1.4). An orbital angular momentum $\ell$ has $2\ell + 1$ orientations, always an odd number. They read the result as Bohr–Sommerfeld space quantization; nobody yet knew that silver's outer electron has $\ell = 0$ and that they had measured something else.

In 1925 Wolfgang Pauli, working out how electrons fill atomic shells, concluded that the electron needs a fourth quantum number, a "two-valuedness not describable classically." Ralph Kronig suggested to Pauli that this was the electron rotating, and Pauli talked him out of publishing. Later that year George Uhlenbeck and Samuel Goudsmit in Leiden proposed the same thing: an intrinsic angular momentum $\hbar/2$. Lorentz, whom they consulted, pointed out that a spinning electron of the classical radius would need a surface moving far faster than light (§3.4 does the calculation). Ehrenfest, their adviser, had already submitted the paper; the name "spin" stuck although the picture is wrong. Pauli wrote down the $2\times2$ matrices in 1927; in 1928 Dirac found that a relativistic wave equation for the electron forces spin-$\frac12$ and $g = 2$ on you.

Underneath is a problem about rotations, which in three dimensions do not commute (LA-1 Example 4.3, CM-7). Angular momentum generates rotations (CM-4: it is the Noether charge of rotational symmetry). So the components of angular momentum cannot commute either, and the way they fail to commute is one equation, $[\hat J_x, \hat J_y] = i\hbar\hat J_z$. This module shows that this equation alone, with no wavefunctions, fixes every possible angular momentum spectrum, including half-integer ones that no orbital motion can produce. Those are spin.

## 2. Intuition

**Angular momentum generates rotations.** In CM-4 and CM-8, $L_z$ generates rotations about $z$ through the Poisson bracket, and $\{L_x, L_y\} = L_z$. LA-12 Appendix A.3 says the same about the generators of $SO(3)$. Dirac's rule (QM-1 §3.11) replaces $\{\cdot,\cdot\}$ by $\frac{1}{i\hbar}[\cdot,\cdot]$, so expect $[\hat L_x, \hat L_y] = i\hbar\hat L_z$.

**You can know $L^2$ and one component, not two components.** Observables share an eigenbasis exactly when they commute (QM-4). The components don't, so no state has definite $L_x$ and $L_z$ (unless both vanish). The length squared $L^2$ is rotation-invariant, so it should commute with every component. Label states by $L^2$ and $L_z$.

**Ladders, again.** QM-3's oscillator had $[\hat H, \hat a_\pm] = \pm\hbar\omega\,\hat a_\pm$, so $\hat a_+$ climbed the energy ladder. Define $\hat L_\pm = \hat L_x \pm i\hat L_y$ and expect $[\hat L_z, \hat L_\pm] = \pm\hbar\hat L_\pm$: $\hat L_+$ raises $L_z$ by $\hbar$ at fixed $L^2$. The new feature is that this ladder cannot go on forever, because $L_z^2 \le L^2$. So it has a top rung and a bottom rung spaced by whole steps of $\hbar$. If the top is $m_{\rm top}\hbar$ and, by $z\to-z$ symmetry, the bottom is $-m_{\rm top}\hbar$, then $2m_{\rm top}$ is a whole number. Guess the theorem: $m$ comes in strings $-j, \ldots, j$ with $j$ an integer **or a half-integer**, $2j + 1$ states per string.

**Why half-integers can't be orbits.** An orbital wavefunction $\psi(\theta,\phi)$ with definite $L_z = \hbar m$ goes like $e^{im\phi}$ (since $\hat L_z = -i\hbar\partial_\phi$). Going once around the axis must return the same function, so $m$ is an integer; $e^{i\phi/2}$ changes sign after $2\pi$. The algebra allows $j = \frac12$, functions on the sphere refuse it, so whatever has $j = \frac12$ is not a function of angles. That is spin.

**Two states make a qubit.** $j = \frac12$ has two states: $\mathbb C^2$, like a polarization (EM-7) or a qubit (QC-1). A unit vector in $\mathbb C^2$ up to phase has two real parameters, so the states fill a sphere, the Bloch sphere, and "spin up along $\hat{\mathbf n}$" is the point $\hat{\mathbf n}$. Up and down along an axis are opposite points, $180°$ apart, though they are orthogonal vectors. That factor of 2 is why half angles appear everywhere below, and why a $360°$ rotation multiplies a spin-$\frac12$ state by $-1$.

**A magnetic field turns the sphere.** A moment $\boldsymbol\mu = \gamma\mathbf S$ in a field feels the torque $\boldsymbol\mu\times\mathbf B$ (EM-5), so $d\mathbf S/dt = \gamma\mathbf S\times\mathbf B$: gyroscopic precession at $\gamma B$. Expect the Bloch vector to do exactly this, and a weak field rotating at that frequency to tip it over (magnetic resonance).

## 3. Formal development

### 3.1 Orbital angular momentum and its commutators

$$
\hat L_x = \hat y\hat p_z - \hat z\hat p_y, \qquad \hat L_y = \hat z\hat p_x - \hat x\hat p_z, \qquad \hat L_z = \hat x\hat p_y - \hat y\hat p_x, \qquad\text{i.e. } \hat L_i = \varepsilon_{ijk}\hat x_j\hat p_k. \tag{3.1}
$$
No ordering problem arises: each product pairs a coordinate with a momentum in a different direction, and those commute. So each $\hat L_i$ is Hermitian (product of commuting Hermitian operators).

▶ **Try it first.** Compute $[\hat L_z, \hat x]$ and $[\hat L_z, \hat p_x]$ using only $[\hat x_i, \hat p_j] = i\hbar\delta_{ij}$. Guess the general pattern.

**Derivation of the commutators.** We use the product rule $[\hat A, \hat B\hat C] = [\hat A, \hat B]\hat C + \hat B[\hat A, \hat C]$ and $[\hat x_i, \hat x_j] = [\hat p_i, \hat p_j] = 0$.

1. $[\hat L_i, \hat x_l] = \varepsilon_{ijk}[\hat x_j\hat p_k, \hat x_l] = \varepsilon_{ijk}\hat x_j[\hat p_k, \hat x_l] = -i\hbar\varepsilon_{ijl}\hat x_j = i\hbar\varepsilon_{ilj}\hat x_j$.
2. $[\hat L_i, \hat p_l] = \varepsilon_{ijk}[\hat x_j, \hat p_l]\hat p_k = i\hbar\varepsilon_{ilk}\hat p_k$.
   Together:
$$
[\hat L_i, \hat x_j] = i\hbar\varepsilon_{ijk}\hat x_k, \qquad [\hat L_i, \hat p_j] = i\hbar\varepsilon_{ijk}\hat p_k. \tag{3.2}
$$
   Any vector operator $\hat{\mathbf V}$ built from $\hat{\mathbf r}$ and $\hat{\mathbf p}$ obeys the same rule, $[\hat L_i, \hat V_j] = i\hbar\varepsilon_{ijk}\hat V_k$, because the product rule passes it through products. That is what "transforms as a vector under rotations" means for operators.
3. Apply this to $\hat L_y = \hat z\hat p_x - \hat x\hat p_z$: $[\hat L_x, \hat L_y] = [\hat L_x, \hat z]\hat p_x + \hat z[\hat L_x, \hat p_x] - [\hat L_x, \hat x]\hat p_z - \hat x[\hat L_x, \hat p_z]$. By (3.2), $[\hat L_x, \hat z] = i\hbar\varepsilon_{xzy}\hat y = -i\hbar\hat y$, $[\hat L_x, \hat p_x] = [\hat L_x, \hat x] = 0$, $[\hat L_x, \hat p_z] = -i\hbar\hat p_y$. So $[\hat L_x, \hat L_y] = -i\hbar\hat y\hat p_x + i\hbar\hat x\hat p_y = i\hbar\hat L_z$.
$$
[\hat L_x, \hat L_y] = i\hbar\hat L_z, \quad [\hat L_y, \hat L_z] = i\hbar\hat L_x, \quad [\hat L_z, \hat L_x] = i\hbar\hat L_y; \qquad [\hat L_i, \hat L_j] = i\hbar\varepsilon_{ijk}\hat L_k. \tag{3.3}
$$
   The others follow by the cyclic relabeling $x\to y\to z\to x$, which maps (3.1) to itself. Example 4.1 redoes the first by brute force.

Define
$$
\hat L^2 = \hat L_x^2 + \hat L_y^2 + \hat L_z^2. \tag{3.4}
$$

4. $[\hat L^2, \hat L_z] = [\hat L_x^2, \hat L_z] + [\hat L_y^2, \hat L_z] + 0$. By the product rule, $[\hat L_x^2, \hat L_z] = \hat L_x[\hat L_x, \hat L_z] + [\hat L_x, \hat L_z]\hat L_x = -i\hbar(\hat L_x\hat L_y + \hat L_y\hat L_x)$, and $[\hat L_y^2, \hat L_z] = \hat L_y(i\hbar\hat L_x) + (i\hbar\hat L_x)\hat L_y = +i\hbar(\hat L_y\hat L_x + \hat L_x\hat L_y)$. They cancel. By cyclic symmetry the same holds for $x$ and $y$:
$$
[\hat L^2, \hat L_i] = 0. \tag{3.5}
$$

■ **What just happened.** Only $[\hat x_i, \hat p_j] = i\hbar\delta_{ij}$ went in. Out came (3.3), the quantum version of CM-4's $\{L_i, L_j\} = \varepsilon_{ijk}L_k$, and (3.5): $\hat L^2$ commutes with every generator (a **Casimir**). So $\hat L^2$ and $\hat L_z$ have a common eigenbasis (QM-4); $\hat L_x$ and $\hat L_z$ do not. Everything below uses only (3.3) and (3.5), so we write $\hat{\mathbf J}$ for any triple obeying them.

### 3.2 The algebraic derivation of the spectrum

This is the derivation in quantum mechanics that you should be able to reproduce from memory. Expect to need three attempts.

**Theorem 3.1 (Angular momentum spectrum).** Let $\hat J_x, \hat J_y, \hat J_z$ be Hermitian with $[\hat J_i, \hat J_j] = i\hbar\varepsilon_{ijk}\hat J_k$, acting on a space where $\hat J^2$ and $\hat J_z$ have a common eigenvector. Then the eigenvalues of $\hat J^2$ are $\hbar^2j(j+1)$ with $j\in\{0, \frac12, 1, \frac32, \ldots\}$; for each such $j$ the eigenvalues of $\hat J_z$ are $\hbar m$ with $m = -j, -j+1, \ldots, j$ ($2j + 1$ values); and with suitable phases
$$
\hat J_\pm|j, m\rangle = \hbar\sqrt{j(j+1) - m(m\pm1)}\;|j, m\pm1\rangle. \tag{3.14}
$$

▶ **Try it first.** Define $\hat J_\pm = \hat J_x \pm i\hat J_y$. Compute $[\hat J_z, \hat J_\pm]$ and then $\hat J_z(\hat J_+|\psi\rangle)$ for an eigenvector $|\psi\rangle$ of $\hat J_z$. Then express $\hat J_-\hat J_+$ in terms of $\hat J^2$ and $\hat J_z$.

**Derivation.**

1. *Ladder operators.* Let
$$
\hat J_\pm = \hat J_x \pm i\hat J_y, \qquad \hat J_\pm^\dagger = \hat J_\mp. \tag{3.6}
$$
   They are not Hermitian; they are each other's adjoints because $\hat J_x$, $\hat J_y$ are Hermitian.
2. *Their commutators.* $[\hat J_z, \hat J_\pm] = [\hat J_z, \hat J_x] \pm i[\hat J_z, \hat J_y] = i\hbar\hat J_y \pm i(-i\hbar\hat J_x) = \pm\hbar(\hat J_x \pm i\hat J_y)$, and $[\hat J^2, \hat J_\pm] = 0$ by (3.5):
$$
[\hat J_z, \hat J_\pm] = \pm\hbar\hat J_\pm, \qquad [\hat J^2, \hat J_\pm] = 0. \tag{3.7}
$$
3. *The products.* $\hat J_\mp\hat J_\pm = (\hat J_x \mp i\hat J_y)(\hat J_x \pm i\hat J_y) = \hat J_x^2 + \hat J_y^2 \pm i[\hat J_x, \hat J_y] = \hat J_x^2 + \hat J_y^2 \mp \hbar\hat J_z$. So
$$
\hat J_\mp\hat J_\pm = \hat J^2 - \hat J_z^2 \mp \hbar\hat J_z. \tag{3.8}
$$
4. *Raising and lowering.* Let $\hat J^2|\lambda, m\rangle = \hbar^2\lambda|\lambda, m\rangle$ and $\hat J_z|\lambda, m\rangle = \hbar m|\lambda, m\rangle$, with $\lambda$, $m$ real (Hermitian operators). By (3.7),
$$
\hat J_z\big(\hat J_\pm|\lambda,m\rangle\big) = \big(\hat J_\pm\hat J_z \pm \hbar\hat J_\pm\big)|\lambda,m\rangle = \hbar(m\pm1)\,\hat J_\pm|\lambda,m\rangle, \qquad \hat J^2\hat J_\pm|\lambda,m\rangle = \hbar^2\lambda\,\hat J_\pm|\lambda,m\rangle. \tag{3.9}
$$
   So $\hat J_\pm|\lambda,m\rangle$ is either zero or an eigenvector with the same $\lambda$ and $m$ shifted by $\pm1$.
5. *The bound.* $\hat J_x^2 + \hat J_y^2$ is a sum of squares of Hermitian operators, so its expectation value is $\|\hat J_x\psi\|^2 + \|\hat J_y\psi\|^2 \ge 0$. In $|\lambda, m\rangle$ it equals $\hbar^2(\lambda - m^2)$. Hence
$$
m^2 \le \lambda. \tag{3.10}
$$
6. *Termination.* Starting at any $m$ and applying $\hat J_+$ repeatedly gives eigenvectors with $m + 1, m + 2, \ldots$ unless one of them is zero. By (3.10) the sequence cannot go on, so there is a largest value $m_t$ with $\hat J_+|\lambda, m_t\rangle = 0$. Likewise there is a smallest $m_b$ with $\hat J_-|\lambda, m_b\rangle = 0$, reached from $m_t$ in an integer number $N\ge0$ of lowering steps: $m_t - m_b = N$.
7. *Top and bottom rungs.* Apply $\hat J_-$ to $\hat J_+|\lambda, m_t\rangle = 0$ and use (3.8): $0 = (\hat J^2 - \hat J_z^2 - \hbar\hat J_z)|\lambda, m_t\rangle = \hbar^2(\lambda - m_t^2 - m_t)|\lambda, m_t\rangle$. Apply $\hat J_+$ to $\hat J_-|\lambda, m_b\rangle = 0$: $0 = \hbar^2(\lambda - m_b^2 + m_b)|\lambda,m_b\rangle$. So
$$
\lambda = m_t(m_t + 1) = m_b(m_b - 1). \tag{3.11}
$$
8. *Solve.* $m_t^2 + m_t - m_b^2 + m_b = (m_t + m_b)(m_t - m_b + 1) = 0$. The second factor is $N + 1 \ge 1$, so $m_b = -m_t$. Write $j \equiv m_t$. Then $N = 2j$ is a non-negative integer, so $j \in \{0, \frac12, 1, \ldots\}$, and
$$
\hat J^2 \to \hbar^2j(j+1), \qquad \hat J_z \to \hbar m, \quad m = -j, -j+1, \ldots, j-1, j. \tag{3.12}
$$
9. *Matrix elements.* By (3.9), $\hat J_+|j,m\rangle = c_m|j, m+1\rangle$ for some number $c_m$. Its squared norm is, by (3.6) and (3.8),
$$
|c_m|^2 = \langle j,m|\hat J_-\hat J_+|j,m\rangle = \hbar^2\big(j(j+1) - m^2 - m\big) = \hbar^2\big(j(j+1) - m(m+1)\big). \tag{3.13}
$$
   Choosing every $c_m$ real and positive (this fixes the relative phases within a multiplet; it is the Condon–Shortley convention) gives the $+$ sign of (3.14). For $\hat J_-$ use $\langle j,m|\hat J_+\hat J_-|j,m\rangle = \hbar^2(j(j+1) - m(m-1))$. ∎

■ **What just happened.** The algebra (3.3) and positivity (3.10) did everything: the ladder stops at both ends, the stopping conditions force a symmetric string, and its length $2j$ is a whole number, so $j$ can be half-integer. Nothing said "orbital." (3.14) is the module's most-used formula; the project's `jmat(j)` builds every angular momentum matrix from it.

**Degeneracy.** If $\hat H$ commutes with $\hat{\mathbf J}$ (rotational symmetry), it commutes with $\hat J_\pm$, so $\hat H\hat J_-|j,m\rangle = E\hat J_-|j,m\rangle$: all $2j + 1$ states of a multiplet share one energy. That is CM-4's "symmetry implies conservation" in the stronger quantum form "symmetry implies degeneracy" (QM-6's multiplets). **Comparison with QM-3.** Same method as the oscillator: shift operators plus positivity. The oscillator's $[\hat a_-, \hat a_+] = 1$ gives a bottom and no top; (3.8) bounds $m$ on both sides, so every multiplet is finite.

### 3.3 Spherical harmonics

Now make $\hat{\mathbf L}$ a differential operator, $\hat{\mathbf p} = -i\hbar\nabla$.

▶ **Try it first.** In spherical coordinates $\nabla = \hat{\mathbf r}\,\partial_r + \hat{\boldsymbol\theta}\,\frac1r\partial_\theta + \hat{\boldsymbol\phi}\,\frac{1}{r\sin\theta}\partial_\phi$ (EM-1). Compute $\mathbf r\times\nabla$ with $\mathbf r = r\hat{\mathbf r}$, using $\hat{\mathbf r}\times\hat{\boldsymbol\theta} = \hat{\boldsymbol\phi}$ and $\hat{\mathbf r}\times\hat{\boldsymbol\phi} = -\hat{\boldsymbol\theta}$. Which variable disappears?

**Derivation of the angular operators.**

1. $\mathbf r\times\nabla = r\hat{\mathbf r}\times\big(\hat{\mathbf r}\partial_r + \hat{\boldsymbol\theta}\frac1r\partial_\theta + \hat{\boldsymbol\phi}\frac{1}{r\sin\theta}\partial_\phi\big) = \hat{\boldsymbol\phi}\,\partial_\theta - \hat{\boldsymbol\theta}\,\frac{1}{\sin\theta}\partial_\phi$. The radial derivative and $r$ drop out:
$$
\hat{\mathbf L} = -i\hbar\Big(\hat{\boldsymbol\phi}\,\partial_\theta - \hat{\boldsymbol\theta}\,\frac{1}{\sin\theta}\partial_\phi\Big). \tag{3.15}
$$
2. Insert the Cartesian components $\hat{\boldsymbol\theta} = (\cos\theta\cos\phi, \cos\theta\sin\phi, -\sin\theta)$ and $\hat{\boldsymbol\phi} = (-\sin\phi, \cos\phi, 0)$:
$$
\hat L_x = i\hbar\big(\sin\phi\,\partial_\theta + \cot\theta\cos\phi\,\partial_\phi\big), \qquad \hat L_y = i\hbar\big(-\cos\phi\,\partial_\theta + \cot\theta\sin\phi\,\partial_\phi\big), \tag{3.16}
$$
$$
\hat L_z = -i\hbar\,\frac{\partial}{\partial\phi}. \tag{3.17}
$$
   (For $\hat L_z$: $\hat{\boldsymbol\phi}_z = 0$ and $-\hat{\boldsymbol\theta}_z/\sin\theta = 1$.) Rotations about $z$ shift $\phi$, and $\hat L_z$ generates them.
3. Combine: $\hat L_\pm = \hat L_x \pm i\hat L_y = i\hbar\big[(\sin\phi \mp i\cos\phi)\partial_\theta + \cot\theta(\cos\phi \pm i\sin\phi)\partial_\phi\big]$. Since $\sin\phi \mp i\cos\phi = \mp ie^{\pm i\phi}$,
$$
\hat L_\pm = \pm\hbar\,e^{\pm i\phi}\Big(\frac{\partial}{\partial\theta} \pm i\cot\theta\,\frac{\partial}{\partial\phi}\Big). \tag{3.18}
$$
4. Use (3.8) in the form $\hat L^2 = \hat L_-\hat L_+ + \hat L_z^2 + \hbar\hat L_z$. Acting on $f(\theta,\phi)$: $\hat L_-\hat L_+f = -\hbar^2e^{-i\phi}(\partial_\theta - i\cot\theta\partial_\phi)\big[e^{i\phi}(f_\theta + i\cot\theta f_\phi)\big]$. The $\partial_\theta$ gives $e^{i\phi}(f_{\theta\theta} - i\csc^2\theta f_\phi + i\cot\theta f_{\theta\phi})$; the $-i\cot\theta\partial_\phi$ gives $-i\cot\theta e^{i\phi}(if_\theta + f_{\theta\phi} - \cot\theta f_\phi + i\cot\theta f_{\phi\phi})$. Adding, the $f_{\theta\phi}$ terms cancel and
   $\hat L_-\hat L_+f = -\hbar^2\big(f_{\theta\theta} + \cot\theta f_\theta + \cot^2\theta f_{\phi\phi} + i(\cot^2\theta - \csc^2\theta)f_\phi\big) = -\hbar^2\big(f_{\theta\theta} + \cot\theta f_\theta + \cot^2\theta f_{\phi\phi} - if_\phi\big)$.
   Adding $\hat L_z^2f = -\hbar^2f_{\phi\phi}$ and $\hbar\hat L_zf = -i\hbar^2f_\phi$, the $f_\phi$ terms cancel and $\cot^2\theta + 1 = \csc^2\theta$:
$$
\hat L^2 = -\hbar^2\left[\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\Big(\sin\theta\frac{\partial}{\partial\theta}\Big) + \frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\phi^2}\right]. \tag{3.19}
$$
   The project checks this identity in SymPy on an arbitrary $f(\theta,\phi)$.

■ **What just happened.** The bracket in (3.19) is the angular part of the Laplacian (EM-3 (3.26)), so
$$
\nabla^2 = \frac{1}{r^2}\frac{\partial}{\partial r}\Big(r^2\frac{\partial}{\partial r}\Big) - \frac{\hat L^2}{\hbar^2r^2}. \tag{3.20}
$$
Kinetic energy splits into radial and centrifugal parts, $\hat L^2/2mr^2$ being CM-5's effective potential $L^2/2mr^2$. QM-6 starts here.

**The eigenfunctions.** Seek $Y(\theta,\phi)$ with $\hat L^2Y = \hbar^2\ell(\ell+1)Y$ and $\hat L_zY = \hbar mY$.

1. *The $\phi$ equation.* $-i\hbar\partial_\phi Y = \hbar mY$ gives $Y = \Theta(\theta)e^{im\phi}$. Single-valuedness, $Y(\theta, \phi + 2\pi) = Y(\theta,\phi)$, forces $e^{2\pi im} = 1$: **$m$ is an integer.** Since $m$ runs from $-\ell$ to $\ell$ in integer steps, **$\ell$ is an integer.** The half-integer multiplets of Theorem 3.1 have no realization as functions on the sphere. Whatever carries $j = \frac12$ is a new degree of freedom (§3.4). Exercise C1 shows that forcing $m = \frac12$ anyway gives a ladder that fails to terminate.
2. *The $\theta$ equation.* Putting $\Theta e^{im\phi}$ into (3.19) and writing $x = \cos\theta$ gives $\frac{d}{dx}\big[(1-x^2)\Theta'\big] + \big[\ell(\ell+1) - \frac{m^2}{1-x^2}\big]\Theta = 0$, the associated Legendre equation of EM-3 §3.8. EM-3's Frobenius argument gives regular solutions $P_\ell^m(\cos\theta)$ only for integer $\ell\ge|m|$. The ladder gives them more directly.
3. *The top rung by a first-order equation.* $\hat L_+Y_\ell^\ell = 0$ with $Y_\ell^\ell = \Theta e^{i\ell\phi}$ is, by (3.18), $\Theta' + i\cot\theta(i\ell)\Theta = \Theta' - \ell\cot\theta\,\Theta = 0$. So $d\ln\Theta/d\theta = \ell\cot\theta$, $\Theta = c\sin^\ell\theta$. Normalizing, $\oint|Y|^2d\Omega = 2\pi|c|^2\int_0^\pi\sin^{2\ell+1}\theta\,d\theta = 2\pi|c|^2\frac{2^{2\ell+1}(\ell!)^2}{(2\ell+1)!} = 1$, and the Condon–Shortley sign is $(-1)^\ell$:
$$
Y_\ell^\ell(\theta,\phi) = \frac{(-1)^\ell}{2^\ell\ell!}\sqrt{\frac{(2\ell+1)!}{4\pi}}\,\sin^\ell\theta\,e^{i\ell\phi}. \tag{3.21}
$$
4. *Every other rung by lowering.* $Y_\ell^{m-1} = \hat L_-Y_\ell^m/\big(\hbar\sqrt{\ell(\ell+1) - m(m-1)}\big)$ from (3.14). Example 4.3 gets $Y_1^0$ from $Y_1^1$; Exercise B3 gets $Y_2^1$ from $Y_2^2$. The result for $\ell\le2$ (the same table as EM-3 §3.8):
$$
\begin{aligned}
&Y_0^0 = \tfrac{1}{\sqrt{4\pi}}, \qquad Y_1^0 = \sqrt{\tfrac{3}{4\pi}}\cos\theta, \qquad Y_1^{\pm1} = \mp\sqrt{\tfrac{3}{8\pi}}\sin\theta\,e^{\pm i\phi},\\
&Y_2^0 = \sqrt{\tfrac{5}{16\pi}}(3\cos^2\theta - 1), \quad Y_2^{\pm1} = \mp\sqrt{\tfrac{15}{8\pi}}\sin\theta\cos\theta\,e^{\pm i\phi}, \quad Y_2^{\pm2} = \sqrt{\tfrac{15}{32\pi}}\sin^2\theta\,e^{\pm2i\phi}.
\end{aligned} \tag{3.22}
$$
   Negative $m$ satisfy $Y_\ell^{-m} = (-1)^m(Y_\ell^m)^*$.
5. *Orthonormality and completeness.* $\hat L^2$ and $\hat L_z$ are Hermitian on $L^2(S^2, d\Omega)$, and $Y_\ell^m$ with different $(\ell, m)$ differ in at least one eigenvalue, so they are orthogonal (QM-4, LA-8 Theorem 3.12): $\oint Y_\ell^{m*}Y_{\ell'}^{m'}d\Omega = \delta_{\ell\ell'}\delta_{mm'}$. Completeness (Sturm–Liouville in $\theta$ plus F2's Fourier series in $\phi$) we state, as EM-3 Theorem 3.6 did.
6. *Parity.* The inversion $\mathbf r\to-\mathbf r$ sends $(\theta,\phi)\to(\pi-\theta, \phi+\pi)$. It commutes with $\hat{\mathbf L}$ (both $\mathbf r$ and $\mathbf p$ flip, so $\mathbf r\times\mathbf p$ doesn't), so it maps each multiplet to itself, and since it commutes with $\hat L_-$ it acts on every member as on the top. On (3.21): $\sin^\ell(\pi - \theta) = \sin^\ell\theta$ and $e^{i\ell(\phi+\pi)} = (-1)^\ell e^{i\ell\phi}$. So
$$
Y_\ell^m(\pi - \theta, \phi + \pi) = (-1)^\ell\,Y_\ell^m(\theta,\phi). \tag{3.23}
$$

■ **What just happened.** Single-valuedness killed the half-integers; the top rung solved a first-order ODE; lowering gave the rest, with no series solution. Parity $(-1)^\ell$ is behind QM-9's dipole selection rule $\Delta\ell = \pm1$ ($\mathbf r$ is odd).

**Pictures and multipoles.** `spherical_harmonics.png` plots $|Y_\ell^m|^2$ as a surface whose radius in direction $(\theta,\phi)$ is the angular probability density; all are symmetric about $z$ because $|e^{im\phi}| = 1$. Chemistry's $p_x$, $p_y$ dumbbells are the real combinations of $Y_1^{\pm1}$, eigenfunctions of $\hat L^2$ but not of $\hat L_z$. EM-3's multipole expansion is a $Y_\ell^m$ expansion via the addition theorem (Jackson §3.6): the dipole moment is a charge distribution's $\ell = 1$ content, the quadrupole its $\ell = 2$ content.

### 3.4 Spin and the Pauli matrices

**Postulate (spin).** Each particle species carries an intrinsic angular momentum $\hat{\mathbf S}$, acting on an internal space $\mathbb C^{2s+1}$, with $[\hat S_i, \hat S_j] = i\hbar\varepsilon_{ijk}\hat S_k$ and a fixed $s$: ($\frac12$ for electrons, protons, neutrons, quarks; $1$ for photons; $0$ for pions). $\hat{\mathbf S}$ commutes with $\hat{\mathbf r}$ and $\hat{\mathbf p}$; an electron's state lives in $L^2(\mathbb R^3)\otimes\mathbb C^2$ (LA-12 §3.11); the generator of rotations of everything is $\hat{\mathbf J} = \hat{\mathbf L} + \hat{\mathbf S}$.

▶ **Try it first.** Use (3.14) with $j = \frac12$ to write $\hat S_+$ and $\hat S_-$ as $2\times2$ matrices in the basis $(|{\uparrow}\rangle, |{\downarrow}\rangle) = (|\frac12,\frac12\rangle, |\frac12,-\frac12\rangle)$. Then recover $\hat S_x$ and $\hat S_y$.

**Derivation of the spin-$\frac12$ matrices.**

1. $\hat S_z = \frac\hbar2\operatorname{diag}(1, -1)$ from (3.12).
2. (3.14) with $j = \frac12$, $m = -\frac12$: $\hat S_+|{\downarrow}\rangle = \hbar\sqrt{\frac34 + \frac14}|{\uparrow}\rangle = \hbar|{\uparrow}\rangle$, and $\hat S_+|{\uparrow}\rangle = 0$. So $\hat S_+ = \hbar\begin{pmatrix}0&1\\0&0\end{pmatrix}$, $\hat S_- = \hat S_+^\dagger = \hbar\begin{pmatrix}0&0\\1&0\end{pmatrix}$.
3. $\hat S_x = \frac12(\hat S_+ + \hat S_-)$ and $\hat S_y = \frac{1}{2i}(\hat S_+ - \hat S_-)$:
$$
\hat S_x = \frac\hbar2\begin{pmatrix}0&1\\1&0\end{pmatrix}, \qquad \hat S_y = \frac\hbar2\begin{pmatrix}0&-i\\i&0\end{pmatrix}, \qquad \hat S_z = \frac\hbar2\begin{pmatrix}1&0\\0&-1\end{pmatrix}. \tag{3.24}
$$
   So $\hat{\mathbf S} = \frac\hbar2\boldsymbol\sigma$, with the **Pauli matrices** $\sigma_x$, $\sigma_y$, $\sigma_z$ the matrices in the brackets. QM-1 (3.13) guessed them from the data; here they are forced by the algebra and the phase convention.

**Theorem 3.2 (Pauli algebra).** For $a, b\in\{x,y,z\}$,
$$
\sigma_a\sigma_b = \delta_{ab}I + i\varepsilon_{abc}\sigma_c, \tag{3.25}
$$
and consequently, for real vectors $\mathbf a$, $\mathbf b$,
$$
(\mathbf a\cdot\boldsymbol\sigma)(\mathbf b\cdot\boldsymbol\sigma) = (\mathbf a\cdot\mathbf b)\,I + i(\mathbf a\times\mathbf b)\cdot\boldsymbol\sigma. \tag{3.26}
$$

**Proof.**

1. *Squares.* $\sigma_x^2 = \begin{pmatrix}0&1\\1&0\end{pmatrix}^2 = I$, $\sigma_y^2 = \begin{pmatrix}0\cdot0 + (-i)(i) & 0\\ 0 & i(-i)\end{pmatrix} = I$, $\sigma_z^2 = I$.
2. *One product.* $\sigma_x\sigma_y = \begin{pmatrix}0&1\\1&0\end{pmatrix}\begin{pmatrix}0&-i\\i&0\end{pmatrix} = \begin{pmatrix}i&0\\0&-i\end{pmatrix} = i\sigma_z$.
3. *Anticommutation.* $(\sigma_x\sigma_y)^\dagger = \sigma_y\sigma_x$ (both Hermitian), and $(i\sigma_z)^\dagger = -i\sigma_z$. So $\sigma_y\sigma_x = -i\sigma_z = -\sigma_x\sigma_y$. The cyclic relabeling gives $\sigma_y\sigma_z = i\sigma_x$, $\sigma_z\sigma_x = i\sigma_y$ (check one directly: $\sigma_z\sigma_x = \begin{pmatrix}0&1\\-1&0\end{pmatrix} = i\sigma_y$ ✓). Steps 1–3 are (3.25).
4. *Vector form.* $(\mathbf a\cdot\boldsymbol\sigma)(\mathbf b\cdot\boldsymbol\sigma) = a_ab_b\sigma_a\sigma_b = a_ab_b(\delta_{ab} + i\varepsilon_{abc}\sigma_c) = \mathbf a\cdot\mathbf b + i(\mathbf a\times\mathbf b)_c\sigma_c$. ∎

Consequences: each $\sigma_a$ is Hermitian and squares to $I$, hence unitary with eigenvalues $\pm1$ (LA-8 Proposition 3.22); $[\sigma_a, \sigma_b] = 2i\varepsilon_{abc}\sigma_c$ is (3.3) again; $\{\sigma_a, \sigma_b\} = 2\delta_{ab}I$; $(\hat{\mathbf n}\cdot\boldsymbol\sigma)^2 = I$; and $\operatorname{tr}(\sigma_a\sigma_b) = 2\delta_{ab}$, so $\{I, \boldsymbol\sigma\}$ is an orthogonal basis of $2\times2$ matrices.

**Spin along any axis.** Measuring spin along the unit vector $\hat{\mathbf n} = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$ is measuring $\hat{\mathbf n}\cdot\hat{\mathbf S} = \frac\hbar2\hat{\mathbf n}\cdot\boldsymbol\sigma$, where
$$
\hat{\mathbf n}\cdot\boldsymbol\sigma = \begin{pmatrix}\cos\theta & e^{-i\phi}\sin\theta\\ e^{i\phi}\sin\theta & -\cos\theta\end{pmatrix}. \tag{3.27}
$$
Its square is $I$, so every spin component has outcomes $\pm\hbar/2$: two Stern–Gerlach spots for every magnet orientation.

▶ **Try it first.** Find the $+1$ eigenvector of (3.27). Use $1 - \cos\theta = 2\sin^2\frac\theta2$ and $\sin\theta = 2\sin\frac\theta2\cos\frac\theta2$.

1. $(\hat{\mathbf n}\cdot\boldsymbol\sigma - I)(a, b)^\top = 0$, top row: $(\cos\theta - 1)a + e^{-i\phi}\sin\theta\,b = 0$, so $b/a = e^{i\phi}(1 - \cos\theta)/\sin\theta = e^{i\phi}\,\frac{2\sin^2(\theta/2)}{2\sin(\theta/2)\cos(\theta/2)} = e^{i\phi}\tan\frac\theta2$.
2. Normalize with $a$ real and positive:
$$
|\hat{\mathbf n}, +\rangle = \cos\tfrac\theta2\,|{\uparrow}\rangle + e^{i\phi}\sin\tfrac\theta2\,|{\downarrow}\rangle, \tag{3.28}
$$
$$
|\hat{\mathbf n}, -\rangle = \sin\tfrac\theta2\,|{\uparrow}\rangle - e^{i\phi}\cos\tfrac\theta2\,|{\downarrow}\rangle. \tag{3.29}
$$
   (3.29) is (3.28) with $\theta\to\pi - \theta$, $\phi\to\phi+\pi$: the eigenvector for the opposite direction $-\hat{\mathbf n}$, as it must be. Special cases: $\hat{\mathbf x}$ ($\theta = \frac\pi2$, $\phi = 0$): $|{\pm x}\rangle = (|{\uparrow}\rangle \pm |{\downarrow}\rangle)/\sqrt2$. $\hat{\mathbf y}$ ($\phi = \frac\pi2$): $|{\pm y}\rangle = (|{\uparrow}\rangle \pm i|{\downarrow}\rangle)/\sqrt2$.

■ **What just happened.** The eigenvectors involve $\theta/2$, not $\theta$: $P(S_z = +\frac\hbar2) = \cos^2\frac\theta2 = \frac12(1 + \cos\theta)$. This is LA-8 Proposition 3.22, now derived; Example 4.4 redoes LA-8's $\theta = 60°$, $\phi = 90°$ case.

**QM-1's Stern–Gerlach sequence, one line each.** (i) SG$z$ twice: $|\langle{\uparrow}|{\uparrow}\rangle|^2 = 1$. (ii) SG$z{+}$ then SG$x$: $|\langle{\pm x}|{\uparrow}\rangle|^2 = \frac12$. (iii) SG$z{+}$, SG$x{+}$, SG$z$: $|\langle{\downarrow}|{+x}\rangle|^2 = \frac12$, so the $z$-information was erased. (iv) SG$z{+}$ then SG at angle $\theta$ in the $xz$-plane: $P_+ = \cos^2\frac\theta2$. Part 8 of the project runs the cascade for any $j$ ($2j + 1$ beams).

**Spin is not the electron spinning.** Model the electron as a uniform solid sphere of mass $m_e$ and the classical electron radius $r_e = e^2/4\pi\varepsilon_0m_ec^2 = 2.818\times10^{-15}$ m (any smaller radius makes it worse). Angular momentum $\frac\hbar2 = I\omega = \frac25m_er_e^2\omega$ gives an equatorial speed $v = \omega r_e = \frac54\frac{\hbar}{m_er_e} = \frac{5}{4\alpha}c \approx 171c$, using $\hbar/m_ec = r_e/\alpha$ with $\alpha\approx1/137$. That is Lorentz's objection, and experiments put the electron's size below $10^{-18}$ m, which makes it far worse. More basically, a spinning body has orbital angular momentum, which §3.3 makes an integer multiple of $\hbar$. Spin obeys (3.3), adds to $\hat{\mathbf L}$ to make the conserved total, and exerts torques, but it is not motion of anything.

### 3.5 The Bloch sphere

**Theorem 3.3 (Bloch sphere).** Every unit vector $|\chi\rangle\in\mathbb C^2$ equals $e^{i\alpha}|\hat{\mathbf n},+\rangle$ for a unique $\hat{\mathbf n}$ on the unit sphere and some phase $\alpha$. Its spin expectation is
$$
\langle\chi|\hat{\mathbf S}|\chi\rangle = \frac\hbar2\,\hat{\mathbf n}, \qquad\text{i.e.}\qquad \langle\boldsymbol\sigma\rangle = \hat{\mathbf n}. \tag{3.30}
$$
Two states are orthogonal if and only if their Bloch vectors are antipodal.

▶ **Try it first.** Write $|\chi\rangle = (a, b)$ with $|a|^2 + |b|^2 = 1$. Strip a phase to make $a$ real and non-negative. Why can you then write $a = \cos\frac\theta2$ with $\theta\in[0,\pi]$?

**Proof.**

1. *Every state is a point.* Multiply by $e^{-i\arg a}$ (if $a = 0$ the state is $|{\downarrow}\rangle$, $\theta = \pi$). Now $a\in[0, 1]$, so $a = \cos\frac\theta2$ for a unique $\theta\in[0,\pi]$, and $|b| = \sin\frac\theta2$, so $b = e^{i\phi}\sin\frac\theta2$, with $\phi$ unique except at the poles. That is (3.28). The half angle lets $\theta\in[0,\pi]$ cover $a\in[0,1]$ once.
2. *Expectation value.* $\langle\sigma_z\rangle = |a|^2 - |b|^2 = \cos\theta$; $\langle\sigma_x\rangle = 2\operatorname{Re}(\bar ab) = 2\cos\frac\theta2\sin\frac\theta2\cos\phi = \sin\theta\cos\phi$; $\langle\sigma_y\rangle = 2\operatorname{Im}(\bar ab) = \sin\theta\sin\phi$. That is $\hat{\mathbf n}$; the global phase cancels.
3. *Antipodes.* Up to phase, the only unit vector orthogonal to $|\hat{\mathbf n},+\rangle$ is $|\hat{\mathbf n},-\rangle$ (3.29), whose Bloch vector is $-\hat{\mathbf n}$. Quantitatively $|\langle\hat{\mathbf n}'|\hat{\mathbf n}\rangle|^2 = \frac12(1 + \hat{\mathbf n}\cdot\hat{\mathbf n}')$ (Exercise B4(c)). ∎

■ **What just happened.** The states of a two-level system (rays in $\mathbb C^2$) form a sphere, and a state's point on it is the direction of its spin expectation. Orthogonal means opposite, not perpendicular. Any two-level system (qubit, two atomic levels, polarization) has this picture. Mixed states (QM-7) fill the ball: $\rho = \frac12(I + \mathbf r\cdot\boldsymbol\sigma)$, $|\mathbf r|\le1$.

### 3.6 Spin in a magnetic field: Larmor precession

A spin carries a magnetic moment $\hat{\boldsymbol\mu} = \gamma\hat{\mathbf S}$ with **gyromagnetic ratio** $\gamma$: $\gamma_e = -g\,e/2m_e$ ($g\approx2.0023$) for the electron, $\gamma_p/2\pi = 42.577$ MHz/T for the proton. With the energy $-\boldsymbol\mu\cdot\mathbf B$ (EM-5),
$$
\hat H = -\gamma\,\mathbf B\cdot\hat{\mathbf S}. \tag{3.31}
$$

▶ **Try it first.** For uniform constant $\mathbf B$, compute $\frac{d}{dt}\hat S_x = \frac i\hbar[\hat H, \hat S_x]$ using (3.3) with $\mathbf B = B\hat{\mathbf z}$.

**Derivation 1 (Heisenberg).**

1. $\frac{d\hat S_i}{dt} = \frac i\hbar[\hat H, \hat S_i] = -\frac{i\gamma}{\hbar}B_j[\hat S_j, \hat S_i] = -\frac{i\gamma}{\hbar}B_j\,i\hbar\varepsilon_{jik}\hat S_k = \gamma\varepsilon_{jik}B_j\hat S_k = \gamma\varepsilon_{ikj}\hat S_kB_j$.
2. That is
$$
\frac{d\hat{\mathbf S}}{dt} = \gamma\,\hat{\mathbf S}\times\mathbf B, \tag{3.32}
$$
   the classical torque equation, exact for the operators (any $s$) and so for their expectation values.
3. With $\mathbf B = B\hat{\mathbf z}$ and $\omega_0 \equiv \gamma B$: $\dot S_x = \omega_0S_y$, $\dot S_y = -\omega_0S_x$, $\dot S_z = 0$. So $\langle\hat S_x\rangle + i\langle\hat S_y\rangle$ obeys $\dot u = -i\omega_0u$, and for an initial Bloch vector at polar angle $\theta_0$ and azimuth $\phi_0$,
$$
\langle\hat{\mathbf S}\rangle(t) = \frac\hbar2\big(\sin\theta_0\cos(\phi_0 - \omega_0t),\ \sin\theta_0\sin(\phi_0 - \omega_0t),\ \cos\theta_0\big). \tag{3.33}
$$

**Derivation 2 (Schrödinger).**

1. With $\mathbf B = B\hat{\mathbf z}$, $\hat H = -\frac{\hbar\omega_0}{2}\sigma_z$ is diagonal with energies $E_\uparrow = -\hbar\omega_0/2$, $E_\downarrow = +\hbar\omega_0/2$ (for $\gamma > 0$ the moment along $\mathbf B$, spin up, is lower).
2. Each component picks up $e^{-iEt/\hbar}$: $|\chi(t)\rangle = e^{i\omega_0t/2}\cos\frac{\theta_0}{2}|{\uparrow}\rangle + e^{-i\omega_0t/2}e^{i\phi_0}\sin\frac{\theta_0}{2}|{\downarrow}\rangle = e^{i\omega_0t/2}\big(\cos\frac{\theta_0}{2}|{\uparrow}\rangle + e^{i(\phi_0 - \omega_0t)}\sin\frac{\theta_0}{2}|{\downarrow}\rangle\big)$.
3. Up to the global phase, this is (3.28) with $\phi(t) = \phi_0 - \omega_0t$ and $\theta$ fixed. By (3.30) the Bloch vector is (3.33).

■ **What just happened.** Both pictures give **Larmor precession**: the Bloch vector sweeps a cone about $\mathbf B$ at $\omega_0 = \gamma B$, clockwise seen from above for $\gamma > 0$, at fixed angle to the field (energy is conserved). The frequency is the splitting over $\hbar$, $\omega_0 = (E_\downarrow - E_\uparrow)/\hbar$. The project compares `solve_ivp` on (3.32), `expm` on the state, and (3.33) (`larmor_bloch.png`).

### 3.7 Magnetic resonance and the Rabi formula

Add to the static field $B_0\hat{\mathbf z}$ a weak field $B_1$ rotating in the $xy$-plane at angular frequency $\omega$, in the same sense as the precession:
$$
\mathbf B(t) = B_0\hat{\mathbf z} + B_1(\cos\omega t\,\hat{\mathbf x} - \sin\omega t\,\hat{\mathbf y}), \qquad \hat H(t) = -\omega_0\hat S_z - \omega_1(\cos\omega t\,\hat S_x - \sin\omega t\,\hat S_y), \tag{3.34}
$$
with $\omega_0 = \gamma B_0$, $\omega_1 = \gamma B_1$ (the **Rabi frequency**), and detuning $\Delta = \omega_0 - \omega$. Start in $|{\uparrow}\rangle$. What is the probability $P_{\uparrow\to\downarrow}(t)$?

▶ **Try it first.** Classically, go to a frame rotating with the field. In that frame $B_1$ is static. What happened to $B_0$? (In a frame rotating at $\omega$, precession at $\omega_0$ looks like precession at $\omega_0 - \omega$.)

**Derivation.**

1. *The time dependence is a rotation.* Let $\hat U(t) = e^{-i\omega t\hat S_z/\hbar}$, a rotation by $\omega t$ about $z$. By §3.8 ((3.37) applied to $\hat{\mathbf S}$), $\hat U^\dagger\hat S_x\hat U = \cos\omega t\,\hat S_x - \sin\omega t\,\hat S_y$. (Directly: $\hat U = \operatorname{diag}(e^{-i\omega t/2}, e^{i\omega t/2})$, $\hat U^\dagger\sigma_x\hat U = \begin{pmatrix}0&e^{i\omega t}\\e^{-i\omega t}&0\end{pmatrix}$.) $\hat S_z$ commutes with $\hat U$. So $\hat H(t) = \hat U^\dagger\hat K\hat U$ with the constant $\hat K = -\omega_0\hat S_z - \omega_1\hat S_x$.
2. *Rotating frame.* Define $|\tilde\chi\rangle = \hat U(t)|\chi\rangle$. Then $i\hbar\partial_t|\tilde\chi\rangle = i\hbar\dot{\hat U}|\chi\rangle + \hat U\,i\hbar\partial_t|\chi\rangle = \omega\hat S_z\hat U|\chi\rangle + \hat U\hat U^\dagger\hat K\hat U|\chi\rangle = (\omega\hat S_z + \hat K)|\tilde\chi\rangle$.
3. *Effective Hamiltonian.* $\hat H_{\rm eff} = -(\omega_0 - \omega)\hat S_z - \omega_1\hat S_x = -\Delta\hat S_z - \omega_1\hat S_x$, time-independent: a static effective field along $(\omega_1, 0, \Delta)/\gamma$. The rotating frame removed $\omega$ of the $\omega_0$; at resonance ($\Delta = 0$) the static field disappears entirely.
4. *Solve.* $\hat H_{\rm eff} = -\frac{\hbar\Omega}{2}\hat{\mathbf m}\cdot\boldsymbol\sigma$ with $\Omega = \sqrt{\omega_1^2 + \Delta^2}$ and $\hat{\mathbf m} = (\omega_1, 0, \Delta)/\Omega$. By (3.36) below, $e^{-i\hat H_{\rm eff}t/\hbar} = \cos\frac{\Omega t}{2}I + i\sin\frac{\Omega t}{2}\,\hat{\mathbf m}\cdot\boldsymbol\sigma$. Acting on $|{\uparrow}\rangle$, the $|{\downarrow}\rangle$ amplitude is $i\sin\frac{\Omega t}{2}\,m_x = i\frac{\omega_1}{\Omega}\sin\frac{\Omega t}2$.
5. *Back to the lab.* $\hat U$ is diagonal, so it changes only phases: $|\langle{\downarrow}|\chi\rangle| = |\langle{\downarrow}|\tilde\chi\rangle|$. Hence the **Rabi formula**
$$
P_{\uparrow\to\downarrow}(t) = \frac{\omega_1^2}{\omega_1^2 + \Delta^2}\,\sin^2\!\Big(\frac{\sqrt{\omega_1^2 + \Delta^2}\;t}{2}\Big). \tag{3.35}
$$

■ **What just happened.** In the rotating frame the spin precesses about a fixed tilted field. At resonance the field is in the equatorial plane and the spin swings from pole to pole: full inversion at $t = \pi/\omega_1$ (a **$\pi$-pulse**), an equal superposition at $\pi/2\omega_1$ (a **$\pi/2$-pulse**). Off resonance it never reaches the south pole; the maximum $\omega_1^2/(\omega_1^2 + \Delta^2)$ is a Lorentzian in $\Delta$ of full width $2\omega_1$ (F3 Example 3.C). The project integrates the lab-frame equation with the full (3.34) and matches (3.35) to $10^{-7}$ (`rabi.png`).

This is how NMR and MRI work. Protons in $B_0 = 1.5$ T or $3$ T precess at $63.9$ or $127.7$ MHz. A coil applies a microtesla $B_1$ at that frequency for about a millisecond, a $\pi/2$-pulse; the tipped magnetization then precesses and induces a voltage in the coil, and field gradients that make $\omega_0$ depend on position encode the image. It is also QC-1's single-qubit gate as a physical operation: a resonant pulse of length $t$ applies $R_{\hat{\mathbf x}}(\omega_1t)$ in the rotating frame, the drive's phase picks the axis, and superconducting and trapped-ion qubits are driven exactly this way.

### 3.8 The rotation operator and $SU(2)\to SO(3)$

A rotation by $\theta$ about $\hat{\mathbf n}$ acts on states as $e^{-i\theta\hat{\mathbf n}\cdot\hat{\mathbf J}/\hbar}$ (QM-4: $\hat{\mathbf J}$ generates rotations as $\hat p$ generates translations). For spin-$\frac12$, $\hat{\mathbf J}/\hbar = \boldsymbol\sigma/2$.

▶ **Try it first.** Expand $e^{-i\alpha\,\hat{\mathbf n}\cdot\boldsymbol\sigma}$ as a power series and use $(\hat{\mathbf n}\cdot\boldsymbol\sigma)^2 = I$. Sort into even and odd powers.

**Derivation.**

1. Let $A = \hat{\mathbf n}\cdot\boldsymbol\sigma$, so $A^2 = I$ by (3.26), $A^{2k} = I$, $A^{2k+1} = A$.
2. $e^{-i\alpha A} = \sum_k\frac{(-i\alpha)^{2k}}{(2k)!}I + \sum_k\frac{(-i\alpha)^{2k+1}}{(2k+1)!}A = \sum_k\frac{(-1)^k\alpha^{2k}}{(2k)!}I - i\sum_k\frac{(-1)^k\alpha^{2k+1}}{(2k+1)!}A$.
3. These are the cosine and sine series (F1's Euler formula). With $\alpha = \theta/2$:
$$
\hat U_{\hat{\mathbf n}}(\theta) \equiv e^{-i\theta\,\hat{\mathbf n}\cdot\boldsymbol\sigma/2} = \cos\frac\theta2\,I - i\sin\frac\theta2\,\hat{\mathbf n}\cdot\boldsymbol\sigma. \tag{3.36}
$$
   It is unitary with determinant $\cos^2\frac\theta2 + \sin^2\frac\theta2 = 1$: an element of $SU(2)$ (LA-12 (A.7)).

**Theorem 3.4 (Spin-$\frac12$ rotations rotate the Bloch vector).** For every real vector $\mathbf v$,
$$
\hat U_{\hat{\mathbf n}}(\theta)\,(\mathbf v\cdot\boldsymbol\sigma)\,\hat U_{\hat{\mathbf n}}(\theta)^\dagger = (R\mathbf v)\cdot\boldsymbol\sigma, \qquad R\mathbf v = \cos\theta\,\mathbf v + \sin\theta\,\hat{\mathbf n}\times\mathbf v + (1-\cos\theta)(\hat{\mathbf n}\cdot\mathbf v)\hat{\mathbf n}, \tag{3.37}
$$
where $R$ is the rotation by $\theta$ about $\hat{\mathbf n}$ (Rodrigues, LA-12 (A.4)). Consequently the Bloch vector of $\hat U|\chi\rangle$ is $R$ times the Bloch vector of $|\chi\rangle$.

**Proof.** Write $c = \cos\frac\theta2$, $s = \sin\frac\theta2$, $N = \hat{\mathbf n}\cdot\boldsymbol\sigma$, $V = \mathbf v\cdot\boldsymbol\sigma$.

1. $(c - isN)V(c + isN) = c^2V + ics(VN - NV) + s^2NVN$.
2. By (3.26), $VN - NV = 2i(\mathbf v\times\hat{\mathbf n})\cdot\boldsymbol\sigma$, so $ics(VN - NV) = -2cs(\mathbf v\times\hat{\mathbf n})\cdot\boldsymbol\sigma = 2cs(\hat{\mathbf n}\times\mathbf v)\cdot\boldsymbol\sigma$.
3. $NVN$: by (3.26), $NV = (\hat{\mathbf n}\cdot\mathbf v) + i(\hat{\mathbf n}\times\mathbf v)\cdot\boldsymbol\sigma$, and then $NVN = (\hat{\mathbf n}\cdot\mathbf v)N + i\big[(\hat{\mathbf n}\times\mathbf v)\cdot\hat{\mathbf n} + i((\hat{\mathbf n}\times\mathbf v)\times\hat{\mathbf n})\cdot\boldsymbol\sigma\big] = (\hat{\mathbf n}\cdot\mathbf v)N - \big(\mathbf v - (\hat{\mathbf n}\cdot\mathbf v)\hat{\mathbf n}\big)\cdot\boldsymbol\sigma$, using $(\hat{\mathbf n}\times\mathbf v)\times\hat{\mathbf n} = \mathbf v - (\hat{\mathbf n}\cdot\mathbf v)\hat{\mathbf n}$. So $NVN = \big(2(\hat{\mathbf n}\cdot\mathbf v)\hat{\mathbf n} - \mathbf v\big)\cdot\boldsymbol\sigma$: the reflection of $\mathbf v$ through the axis.
4. Collect: $(c^2 - s^2)\mathbf v + 2cs\,\hat{\mathbf n}\times\mathbf v + 2s^2(\hat{\mathbf n}\cdot\mathbf v)\hat{\mathbf n}$. With $c^2 - s^2 = \cos\theta$, $2cs = \sin\theta$, $2s^2 = 1 - \cos\theta$, this is $R\mathbf v$.
5. For the Bloch vector: $\langle\hat U\chi|\sigma_a|\hat U\chi\rangle = \langle\chi|\hat U^\dagger\sigma_a\hat U|\chi\rangle$. By (3.37) with $\hat U^\dagger = \hat U_{\hat{\mathbf n}}(-\theta)$, $\hat U^\dagger\sigma_a\hat U = (R^{-1}\mathbf e_a)\cdot\boldsymbol\sigma$, so the new Bloch vector has components $(R^{-1}\mathbf e_a)\cdot\mathbf r = \mathbf e_a\cdot R\mathbf r$ (as $R^{-1} = R^\top$). That is $R\mathbf r$. ∎

■ **What just happened.** $e^{-i\theta\hat{\mathbf n}\cdot\boldsymbol\sigma/2}$ turns the Bloch sphere rigidly by $\theta$: the half angle becomes a full angle because $\hat U$ appears twice in (3.37). This is LA-12 Appendix A.6's 2-to-1 map $SU(2)\to SO(3)$ for every axis; $\pm\hat U$ give the same $R$. At $\theta = 2\pi$, (3.36) gives $\hat U = -I$: the Bloch vector returns but the state changes sign. Only $4\pi$ gives $+I$.

**Is the $-1$ observable?** A global phase is not; a relative phase is. In 1975 Rauch and collaborators (at the Institut Laue–Langevin) and Werner, Colella, Overhauser and Eagen (University of Missouri) split a neutron beam in a silicon-crystal interferometer, precessed the spin on one path by a controllable angle in a magnetic field, and recombined. The interference oscillated with period $4\pi$ in the precession angle, not $2\pi$.

**The QC-1 vocabulary.** The single-qubit rotation gates are (3.36) about the coordinate axes:
$$
R_x(\theta) = \begin{pmatrix}\cos\frac\theta2 & -i\sin\frac\theta2\\ -i\sin\frac\theta2 & \cos\frac\theta2\end{pmatrix}, \quad R_y(\theta) = \begin{pmatrix}\cos\frac\theta2 & -\sin\frac\theta2\\ \sin\frac\theta2 & \cos\frac\theta2\end{pmatrix}, \quad R_z(\theta) = \begin{pmatrix}e^{-i\theta/2} & 0\\ 0 & e^{i\theta/2}\end{pmatrix}. \tag{3.38}
$$
Every $2\times2$ unitary is $e^{i\alpha}R_z(\beta)R_y(\gamma)R_z(\delta)$ (Euler angles, CM-7; QC-1 proves it). Familiar gates are rotations up to phase: $X = \sigma_x = iR_x(\pi)$, $Z = iR_z(\pi)$, and the Hadamard $H = (\sigma_x + \sigma_z)/\sqrt2 = iR_{\hat{\mathbf n}}(\pi)$ with $\hat{\mathbf n} = (\hat{\mathbf x} + \hat{\mathbf z})/\sqrt2$ (Exercise B5).

### 3.9 Higher spins and the photon

(3.14) builds the matrices for any $j$. For spin 1 (basis $m = 1, 0, -1$), $\hat S_+|1,0\rangle = \sqrt2\hbar|1,1\rangle$ and $\hat S_+|1,-1\rangle = \sqrt2\hbar|1,0\rangle$, so
$$
\hat S_x = \frac{\hbar}{\sqrt2}\begin{pmatrix}0&1&0\\1&0&1\\0&1&0\end{pmatrix}, \qquad \hat S_y = \frac{\hbar}{\sqrt2}\begin{pmatrix}0&-i&0\\i&0&-i\\0&i&0\end{pmatrix}, \qquad \hat S_z = \hbar\begin{pmatrix}1&0&0\\0&0&0\\0&0&-1\end{pmatrix}. \tag{3.39}
$$
$\hat S_x^2 = \frac{\hbar^2}{2}\begin{pmatrix}1&0&1\\0&2&0\\1&0&1\end{pmatrix}$ is not a multiple of $I$, so there is no product rule like (3.25). A spin-1 state has four real parameters after normalization and phase, so it is not a point on a sphere: besides $\langle\hat{\mathbf S}\rangle$ it has a quadrupole part $\frac12\langle\hat S_i\hat S_j + \hat S_j\hat S_i\rangle - \frac23\hbar^2\delta_{ij}$ ("Your turn" item 1). A Stern–Gerlach magnet splits spin 1 into three beams (Exercise B7).

The photon has spin 1, but a massless particle has no rest frame and only its helicity (angular momentum along the motion) is defined. Only $m = \pm1$ occur; $m = 0$ would be a longitudinal wave, which $\nabla\cdot\mathbf E = 0$ forbids in vacuum (EM-7). The two helicities are EM-7's two circular polarizations. We state this; deriving it needs relativistic field theory.

### 3.10 Addition of angular momenta

Two subsystems with angular momenta $\hat{\mathbf J}_1$ (on $V_1$, dimension $2j_1 + 1$) and $\hat{\mathbf J}_2$ (on $V_2$) combine on $V_1\otimes V_2$ (LA-12 §3.11). The total is
$$
\hat{\mathbf J} = \hat{\mathbf J}_1\otimes I + I\otimes\hat{\mathbf J}_2, \tag{3.40}
$$
usually written $\hat{\mathbf J}_1 + \hat{\mathbf J}_2$; in the project it is `np.kron(J1, I2) + np.kron(I1, J2)` (LA-12 §3.14). Operators on different factors commute, so $\hat{\mathbf J}$ obeys (3.3) and Theorem 3.1 applies. The product basis $|j_1m_1\rangle|j_2m_2\rangle$ diagonalizes $\hat J_z$ ($m = m_1 + m_2$) but not $\hat J^2$; the **coupled basis** $|j, m\rangle$ diagonalizes $\hat J_1^2$, $\hat J_2^2$, $\hat J^2$, $\hat J_z$. The coefficients $\langle j_1m_1;j_2m_2|jm\rangle$ relating them are the **Clebsch–Gordan coefficients**.

**Two spin-$\frac12$'s.** Four product states $|{\uparrow\uparrow}\rangle, |{\uparrow\downarrow}\rangle, |{\downarrow\uparrow}\rangle, |{\downarrow\downarrow}\rangle$ with $m = 1, 0, 0, -1$.

▶ **Try it first.** $|{\uparrow\uparrow}\rangle$ has the largest possible $m = 1$, so it must be the top of a $j = 1$ multiplet. Apply $\hat S_- = \hat S_{1-} + \hat S_{2-}$ to it. Then find the $m = 0$ state orthogonal to what you got.

**Derivation.**

1. *Top.* $m = 1$ occurs once, so $|{\uparrow\uparrow}\rangle$ is an $\hat S^2$ eigenvector ($\hat S^2$ commutes with $\hat S_z$); $\hat S_+$ kills it (no $m = 2$), so it is $|1, 1\rangle$.
2. *Lower.* Left side by (3.14) for $j = 1$: $\hat S_-|1,1\rangle = \hbar\sqrt2|1,0\rangle$. Right side, with each $\hat S_{i-}|{\uparrow}\rangle = \hbar|{\downarrow}\rangle$: $(\hat S_{1-} + \hat S_{2-})|{\uparrow\uparrow}\rangle = \hbar(|{\downarrow\uparrow}\rangle + |{\uparrow\downarrow}\rangle)$. So $|1,0\rangle = \frac{1}{\sqrt2}(|{\uparrow\downarrow}\rangle + |{\downarrow\uparrow}\rangle)$.
3. *Lower again.* $\hat S_-|1,0\rangle = \hbar\sqrt2|1,-1\rangle$ and $(\hat S_{1-} + \hat S_{2-})\frac{1}{\sqrt2}(|{\uparrow\downarrow}\rangle + |{\downarrow\uparrow}\rangle) = \frac{\hbar}{\sqrt2}(2|{\downarrow\downarrow}\rangle)$. So $|1,-1\rangle = |{\downarrow\downarrow}\rangle$.
4. *The remaining state.* The $m = 0$ subspace is two-dimensional; one direction is $|1,0\rangle$. The orthogonal direction $\frac{1}{\sqrt2}(|{\uparrow\downarrow}\rangle - |{\downarrow\uparrow}\rangle)$ is an $\hat S_z$ eigenvector with $m = 0$. $\hat S_+$ on it gives $\frac{\hbar}{\sqrt2}(|{\uparrow\uparrow}\rangle - |{\uparrow\uparrow}\rangle) = 0$, so it is a top rung with $m = 0$: $j = 0$.

$$
\text{triplet } (s=1):\ |{\uparrow\uparrow}\rangle,\ \tfrac{1}{\sqrt2}(|{\uparrow\downarrow}\rangle + |{\downarrow\uparrow}\rangle),\ |{\downarrow\downarrow}\rangle; \qquad \text{singlet } (s = 0):\ \tfrac{1}{\sqrt2}(|{\uparrow\downarrow}\rangle - |{\downarrow\uparrow}\rangle). \tag{3.41}
$$

■ **What just happened.** $\frac12\otimes\frac12 = 1\oplus0$, $4 = 3 + 1$. The triplet is symmetric under swapping the spins, the singlet antisymmetric. Every component of $\hat{\mathbf S}$ annihilates the singlet (Exercise B8), so it is rotation-invariant: the same form along every axis. Its coefficient matrix $\frac{1}{\sqrt2}\begin{pmatrix}0&1\\-1&0\end{pmatrix}$ has rank 2, so it is not a product state (LA-12 §3.12): it is the curriculum's first **entangled** state. Measure one spin along any axis and the other is certainly opposite; QM-7 builds EPR and Bell on it. The project diagonalizes the $4\times4$ $\hat S^2$: eigenvalues $0, 2\hbar^2, 2\hbar^2, 2\hbar^2$.

**The general rule.** The top state $|j_1j_1\rangle|j_2j_2\rangle$ starts $j = j_1 + j_2$. The $m = j_1 + j_2 - 1$ subspace is two-dimensional: lowering uses one direction, the other starts $j_1 + j_2 - 1$. Each step down in $m$ adds one product state, hence one new multiplet, until $m = |j_1 - j_2|$, after which the count per $m$ stops growing. Hence
$$
j_1\otimes j_2 = (j_1 + j_2)\oplus(j_1 + j_2 - 1)\oplus\cdots\oplus|j_1 - j_2|, \qquad (2j_1 + 1)(2j_2 + 1) = \sum_{j = |j_1 - j_2|}^{j_1 + j_2}(2j+1). \tag{3.42}
$$
Check the count for $j_1\ge j_2$: $\sum_{k=0}^{2j_2}\big(2(j_1 - j_2 + k) + 1\big) = (2j_2 + 1)(2j_1 - 2j_2 + 1) + 2j_2(2j_2+1) = (2j_2+1)(2j_1 + 1)$. ✓ Every Clebsch–Gordan table (Griffiths Table 4.8, the Particle Data Group's) is the output of this procedure with the sign rule $\langle j_1j_1; j_2(j - j_1)|jj\rangle > 0$. To read one: columns are $(j, m)$, rows $(m_1, m_2)$, and the square root is implied (the PDG's $-1/3$ means $-\sqrt{1/3}$).

### 3.11 Spin–orbit coupling, spin and statistics, the Zeeman effect (previews)

**Spin–orbit coupling (QM-8).** In the electron's rest frame the nucleus circles it and makes a magnetic field $\propto\mathbf L$, adding $\xi(r)\hat{\mathbf L}\cdot\hat{\mathbf S}$ to $\hat H$. Since $\hat{\mathbf L}\cdot\hat{\mathbf S} = \frac12(\hat J^2 - \hat L^2 - \hat S^2)$, this commutes with $\hat J^2$, $\hat J_z$, $\hat L^2$, $\hat S^2$ but not with $\hat L_z$ or $\hat S_z$: the coupled basis is the right one. A $p$ electron splits into $j = \frac32$ and $j = \frac12$ (Example 4.7) with $\langle\hat{\mathbf L}\cdot\hat{\mathbf S}\rangle = \frac{\hbar^2}{2}\big(j(j+1) - 2 - \frac34\big) = \frac{\hbar^2}{2}$ and $-\hbar^2$. That is the sodium D doublet (589.0 and 589.6 nm) and hydrogen's fine structure.

**Spin and statistics (QM-7).** Relativistic quantum field theory proves (Pauli 1940) that identical integer-spin particles (bosons) have exchange-symmetric wavefunctions and half-integer-spin particles (fermions) antisymmetric ones. Antisymmetry gives the exclusion principle and the periodic table. QM-7 takes this as a postulate.

**The Zeeman effect (QM-8).** In a weak field $B\hat{\mathbf z}$ a level with total angular momentum $j$ splits into $2j + 1$ equally spaced sublevels, $\Delta E = g_J\mu_BBm$, with $\mu_B = e\hbar/2m_e$ and the Landé factor $g_J$ from projecting $\hat{\mathbf L} + 2\hat{\mathbf S}$ onto $\hat{\mathbf J}$. Spin's $g = 2$ against the orbit's $g = 1$ is why the patterns of §1 were "anomalous," and half-integer $j$ gives an even number of sublevels.

## 4. Worked examples

### Example 4.1: $[\hat L_x, \hat L_y]$ by brute force

Expand both products from (3.1), keeping operator order:
$\hat L_x\hat L_y = (\hat y\hat p_z - \hat z\hat p_y)(\hat z\hat p_x - \hat x\hat p_z) = \hat y\hat p_z\hat z\hat p_x - \hat y\hat p_z\hat x\hat p_z - \hat z\hat p_y\hat z\hat p_x + \hat z\hat p_y\hat x\hat p_z$,
$\hat L_y\hat L_x = (\hat z\hat p_x - \hat x\hat p_z)(\hat y\hat p_z - \hat z\hat p_y) = \hat z\hat p_x\hat y\hat p_z - \hat z\hat p_x\hat z\hat p_y - \hat x\hat p_z\hat y\hat p_z + \hat x\hat p_z\hat z\hat p_y$.

Compare term by term.

- Terms 2 and 7: $\hat y\hat p_z\hat x\hat p_z = \hat x\hat y\hat p_z^2 = \hat x\hat p_z\hat y\hat p_z$ (all factors commute: $\hat x, \hat y$ with $\hat p_z$). They cancel in the difference.
- Terms 3 and 6: $\hat z\hat p_y\hat z\hat p_x = \hat z^2\hat p_x\hat p_y = \hat z\hat p_x\hat z\hat p_y$ ($\hat p_x, \hat p_y$ commute with $\hat z$). Cancel.
- Terms 1 and 5: $\hat y\hat p_z\hat z\hat p_x - \hat z\hat p_x\hat y\hat p_z = \hat y\hat p_x(\hat p_z\hat z - \hat z\hat p_z) = \hat y\hat p_x(-i\hbar)$.
- Terms 4 and 8: $\hat z\hat p_y\hat x\hat p_z - \hat x\hat p_z\hat z\hat p_y = \hat x\hat p_y(\hat z\hat p_z - \hat p_z\hat z) = \hat x\hat p_y(i\hbar)$.

Total: $[\hat L_x, \hat L_y] = i\hbar(\hat x\hat p_y - \hat y\hat p_x) = i\hbar\hat L_z$. Only the pairs where $\hat z$ meets $\hat p_z$ survive. The project checks this in SymPy on an arbitrary $g(x,y,z)$.

### Example 4.2: the termination step, with numbers

Suppose you have found an eigenvector of $\hat L^2$ with eigenvalue $\lambda\hbar^2 = 2\hbar^2$ and of $\hat L_z$ with $m = 1$. Can you raise it?

$\|\hat L_+|\lambda, 1\rangle\|^2 = \langle\lambda,1|\hat L_-\hat L_+|\lambda,1\rangle = \hbar^2(\lambda - m^2 - m) = \hbar^2(2 - 1 - 1) = 0$. So $\hat L_+|\lambda, 1\rangle = 0$: $m = 1$ is the top rung, and $\lambda = m_t(m_t + 1) = 2$ is consistent with $\ell = 1$.

Suppose instead someone claims $\lambda = 2$ and $m = 2$. Then (3.10) fails: $\langle\hat L_x^2 + \hat L_y^2\rangle = \hbar^2(2 - 4) = -2\hbar^2 < 0$. Equivalently, (3.13) would give $\|\hat L_+|\lambda,2\rangle\|^2 = \hbar^2(2 - 6) = -4\hbar^2$, a negative squared length. The ladder stops because continuing would need vectors of negative length. At the bottom: $\|\hat L_-|\lambda,-1\rangle\|^2 = \hbar^2(2 - 1 - 1) = 0$ and $\lambda = m_b(m_b - 1) = 2$ ✓, so $m_t = -m_b = 1$: three states.

### Example 4.3: $Y_1^0$ from $Y_1^1$

Start from (3.21) with $\ell = 1$: $Y_1^1 = -\frac{1}{2}\sqrt{\frac{3!}{4\pi}}\sin\theta e^{i\phi} = -\sqrt{\frac{3}{8\pi}}\sin\theta\,e^{i\phi}$. Check the top-rung condition: by (3.18), $\hat L_+Y_1^1 \propto e^{i\phi}(\partial_\theta + i\cot\theta\,\partial_\phi)(\sin\theta e^{i\phi}) = e^{2i\phi}(\cos\theta - \cot\theta\sin\theta) = 0$. ✓

Lower: $\hat L_- = -\hbar e^{-i\phi}(\partial_\theta - i\cot\theta\,\partial_\phi)$.

1. $\partial_\theta(\sin\theta e^{i\phi}) = \cos\theta\,e^{i\phi}$.
2. $-i\cot\theta\,\partial_\phi(\sin\theta e^{i\phi}) = -i\cot\theta\cdot i\sin\theta e^{i\phi} = \cos\theta\,e^{i\phi}$.
3. Sum: $2\cos\theta\,e^{i\phi}$. So $\hat L_-Y_1^1 = -\hbar e^{-i\phi}\cdot\big(-\sqrt{\tfrac{3}{8\pi}}\big)\cdot2\cos\theta\,e^{i\phi} = 2\hbar\sqrt{\frac{3}{8\pi}}\cos\theta$.
4. (3.14): $\hat L_-|1,1\rangle = \hbar\sqrt{1\cdot2 - 1\cdot0}\,|1,0\rangle = \sqrt2\hbar|1,0\rangle$. Therefore $Y_1^0 = \frac{2}{\sqrt2}\sqrt{\frac{3}{8\pi}}\cos\theta = \sqrt{\frac{3}{4\pi}}\cos\theta$, matching (3.22). ✓

Check normalization: $2\pi\cdot\frac{3}{4\pi}\int_{-1}^1x^2dx = 1$ ✓. The Condon–Shortley minus in $Y_1^1$ and the minus in $\hat L_-$ combine to make $Y_1^0$ positive; drop either and the signs disagree with (3.14). The project checks $\hat L_-Y_1^1 = \sqrt2\,Y_1^0$ by finite differences on SciPy's harmonics.

### Example 4.4: spin along $\theta = 60°$, $\phi = 90°$

$\hat{\mathbf n} = (\sin60°\cos90°, \sin60°\sin90°, \cos60°) = (0, \frac{\sqrt3}{2}, \frac12)$: in the $yz$-plane, $60°$ from $z$. By (3.27), $\hat{\mathbf n}\cdot\boldsymbol\sigma = \begin{pmatrix}\frac12 & -i\frac{\sqrt3}{2}\\ i\frac{\sqrt3}2 & -\frac12\end{pmatrix}$.

By (3.28) with $\theta/2 = 30°$ and $e^{i\phi} = i$: $|\hat{\mathbf n}, +\rangle = \cos30°|{\uparrow}\rangle + i\sin30°|{\downarrow}\rangle = (\frac{\sqrt3}2, \frac i2)^\top = (0.8660, 0.5i)$. Check: the top entry of $(\hat{\mathbf n}\cdot\boldsymbol\sigma)|\hat{\mathbf n},+\rangle$ is $\frac12\cdot\frac{\sqrt3}2 + (-i\frac{\sqrt3}2)(\frac i2) = \frac{\sqrt3}4 + \frac{\sqrt3}4 = \frac{\sqrt3}2$ ✓; the bottom is $i\frac{\sqrt3}{2}\cdot\frac{\sqrt3}2 - \frac12\cdot\frac i2 = \frac{3i}4 - \frac i4 = \frac i2$ ✓. By (3.29), $|\hat{\mathbf n}, -\rangle = (\sin30°, -i\cos30°) = (\frac12, -i\frac{\sqrt3}2)$.

Measure $\hat S_z$ on $|\hat{\mathbf n},+\rangle$: $P(+\frac\hbar2) = \cos^2 30° = \frac34$, $P(-\frac\hbar2) = \frac14$; on $|\hat{\mathbf n},-\rangle$ the reverse. Check with (3.30): $\langle\hat S_z\rangle = \frac\hbar2(\frac34 - \frac14) = \frac\hbar2\cos60°$ ✓, and $\langle\sigma_y\rangle = 2\operatorname{Im}(\bar ab) = 2\cdot\frac{\sqrt3}2\cdot\frac12 = \frac{\sqrt3}2 = n_y$ ✓. The azimuth enters only as the phase $i$, invisible to $\hat S_z$. With $\theta$ in place of $\theta/2$ you would get $P(+) = \cos^260° = \frac14$, wrong. This is LA-8 Example 4.4 with a physical meaning: a Stern–Gerlach magnet tilted $60°$ from $z$ toward $y$.

### Example 4.5: a proton in 1 T, and MRI

$\gamma_p = 2.675222\times10^8$ rad s$^{-1}$T$^{-1}$ (CODATA). At $B = 1$ T: $f_0 = \gamma_pB/2\pi = 42.577$ MHz, period $23.49$ ns. Clinical scanners: $1.5$ T gives $63.87$ MHz, $3$ T gives $127.73$ MHz, radio frequencies.

Start a proton in $|{+x}\rangle$ ($\theta_0 = 90°$, $\phi_0 = 0$). By (3.33), $\langle\hat S_x\rangle(t) = \frac\hbar2\cos\omega_0t$, $\langle\hat S_y\rangle(t) = -\frac\hbar2\sin\omega_0t$, $\langle\hat S_z\rangle = 0$. At $t = 5$ ns: $\omega_0t = 1.3376$ rad, so $\langle\hat S_x\rangle = 0.1155\hbar$, $\langle\hat S_y\rangle = -0.4865\hbar$. The spin has turned $76.6°$ clockwise (seen from $+z$) in 5 ns.

Why MRI needs many protons: $\hbar\omega_0 = 1.76\times10^{-7}$ eV at 1 T against $k_BT = 0.027$ eV at body temperature, so the excess population in $|{\uparrow}\rangle$ is $\hbar\omega_0/2k_BT\approx3.3\times10^{-6}$. Water's $6.7\times10^{22}$ protons per cm³ make up for it.

### Example 4.6: Rabi flips at resonance and at 10% detuning

Units of $\omega_1$. **At resonance**, $\Delta = 0$: (3.35) gives $P = \sin^2(\omega_1t/2)$. A $\pi$-pulse, $t_\pi = \pi/\omega_1$, gives $P = \sin^2(\pi/2) = 1.000000$: complete inversion.

**At $\Delta = 0.1\,\omega_1$**, same pulse length: $\Omega = \omega_1\sqrt{1.01} = 1.004988\,\omega_1$, so $\Omega t_\pi/2 = 1.578631$ rad, $\sin^2 = 0.999939$, and $P = \frac{1}{1.01}\times0.999939 = 0.990038$. The best you can do at this detuning is $P_{\max} = 1/1.01 = 0.990099$, reached at $\omega_1t = \pi/\sqrt{1.01} = 3.1260$. A detuning of 10% of the Rabi frequency costs you 1% of the flip; the error is second order in $\Delta/\omega_1$, which is why resonant pulses are forgiving. At $\Delta = \omega_1$, $P_{\max} = \frac12$ (Exercise B6).

For MRI: $B_1 = 10\ \mu$T gives $\omega_1/2\pi = 425.8$ Hz and $t_\pi = 1.174$ ms. The project integrates the lab-frame equation for $\Delta/\omega_1 \in\{0, 0.1, 0.5, 1, 2\}$ (`rabi.png`).

### Example 4.7: $\ell = 1$ with spin $\frac12$, by lowering; the counting

An electron in a $p$ orbital: $j_1 = \ell = 1$, $j_2 = s = \frac12$. Product states: $3\times2 = 6$. (3.42) predicts $1\otimes\frac12 = \frac32\oplus\frac12$, and $4 + 2 = 6$ ✓. Write $|m_\ell\rangle|m_s\rangle$ as $|m_\ell, {\uparrow}\rangle$ etc.

1. *Top.* $|\frac32,\frac32\rangle = |1, {\uparrow}\rangle$ (the only state with $m = \frac32$).
2. *Lower.* Left: $\hat J_-|\frac32,\frac32\rangle = \hbar\sqrt{\frac{15}4 - \frac32\cdot\frac12}\,|\frac32,\frac12\rangle = \hbar\sqrt3\,|\frac32,\frac12\rangle$. Right: $\hat L_-|1\rangle = \hbar\sqrt{2}|0\rangle$ and $\hat S_-|{\uparrow}\rangle = \hbar|{\downarrow}\rangle$, so $(\hat L_- + \hat S_-)|1,{\uparrow}\rangle = \hbar\sqrt2|0,{\uparrow}\rangle + \hbar|1,{\downarrow}\rangle$. Hence
$$|\tfrac32,\tfrac12\rangle = \sqrt{\tfrac23}\,|0,{\uparrow}\rangle + \sqrt{\tfrac13}\,|1,{\downarrow}\rangle .$$
3. *Orthogonal complement.* The $m = \frac12$ subspace is spanned by $|0,{\uparrow}\rangle$ and $|1,{\downarrow}\rangle$. The orthogonal unit vector, with the Condon–Shortley sign (coefficient of the $m_1 = j_1 = 1$ state positive), is
$$|\tfrac12,\tfrac12\rangle = \sqrt{\tfrac23}\,|1,{\downarrow}\rangle - \sqrt{\tfrac13}\,|0,{\uparrow}\rangle .$$
   Check that it is a top rung: $\hat L_+|1\rangle = 0$, $\hat S_+|{\downarrow}\rangle = \hbar|{\uparrow}\rangle$, $\hat L_+|0\rangle = \hbar\sqrt2|1\rangle$, so $\hat J_+|\frac12,\frac12\rangle = \hbar\big(\sqrt{\frac23} - \sqrt{\frac13}\sqrt2\big)|1,{\uparrow}\rangle = 0$ ✓. A top rung with $m = \frac12$: $j = \frac12$.
4. *Finish by lowering* (Exercise C2 completes the rest): $|\frac32,-\frac12\rangle = \sqrt{\frac13}|{-1},{\uparrow}\rangle + \sqrt{\frac23}|0,{\downarrow}\rangle$, $|\frac32,-\frac32\rangle = |{-1},{\downarrow}\rangle$, $|\frac12,-\frac12\rangle = \sqrt{\frac13}|0,{\downarrow}\rangle - \sqrt{\frac23}|{-1},{\uparrow}\rangle$.

Read off: $\langle1,0;\frac12,\frac12|\frac32,\frac12\rangle = \sqrt{2/3} = 0.8165$, $\langle1,1;\frac12,-\frac12|\frac12,\frac12\rangle = \sqrt{2/3}$, $\langle1,0;\frac12,\frac12|\frac12,\frac12\rangle = -\sqrt{1/3} = -0.5774$, matching Griffiths Table 4.8 and SymPy's `CG` (checked for every entry, with dimensions $\{4, 2\}$ and $J^2$ eigenvalues $\frac{15}4\hbar^2$, $\frac34\hbar^2$). In $|\frac12,\frac12\rangle$ a measurement of $L_z$ gives $m_\ell = 1$ with probability $\frac23$. These are hydrogen's $2P_{3/2}$ and $2P_{1/2}$ levels (§3.11).

## 5. Why it matters and where it connects

- **QM-6.** Hydrogen is $R_{n\ell}(r)Y_\ell^m$; the radial equation carries the centrifugal term of (3.20); the $2\ell + 1$ degeneracy is §3.2's symmetry-degeneracy remark, and the extra degeneracy in $\ell$ is the Laplace–Runge–Lenz vector (CM-4 Example 4.6).
- **QM-7.** The singlet (3.41) is the first entangled state; the Bloch ball becomes the density matrix; §3.11's spin-statistics statement fixes identical-particle symmetry.
- **QM-8, QM-9.** Spin–orbit and Zeeman shifts are first-order perturbation theory in Example 4.7's coupled basis. Selection rules come from parity (3.23) and (3.42) applied to the $\ell = 1$ dipole operator.
- **QC-1.** Qubit $=$ spin-$\frac12$: Theorem 3.3 is the Bloch sphere, (3.36)–(3.38) are the gates, the $-1$ at $2\pi$ is why gates are "up to global phase," and Rabi pulses (3.35) are how gates are applied.
- **CM-7, LA-12.** Rigid rotations are $SO(3)$; spin-$\frac12$ lives on its double cover $SU(2)$; Theorem 3.4 completes LA-12 Appendix A.6 for every axis.
- **EM-3, EM-7, SM-4.** EM-3's $Y_\ell^m$ are §3.3's; EM-7's polarizations are photon helicities; SM-4's spin-$\frac12$ paramagnet is §3.6's two levels in a heat bath, and Example 4.5's $3\times10^{-6}$ polarization is its high-temperature limit.

## 6. Exercises

See `exercises.md`. A (6): commutators, $\ell = 1$ matrices, Pauli products, spin states, gyromagnetic numbers, counting. B (8): central potentials, uncertainty in $|\ell, m\rangle$, $Y_2^1$ by lowering, a Bloch vector, gates as rotations, the Rabi line shape, spin-1 Stern–Gerlach, the singlet. C (3): half-integer $\ell$ on the sphere, a general Clebsch–Gordan code, the $SU(2)\to SO(3)$ kernel. B2, B3, B4, B8 are essential.

## 7. Coding project

`project.py` (about 2 s). Part 1: `jmat(j)` builds $\hat J_{x,y,z,\pm}$ from (3.14) and asserts (3.3), (3.5), (3.7), (3.8), $\hat J^2 = j(j+1)$ for $j = \frac12, 1, \frac32, 2$. Part 2: SciPy's `sph_harm_y` against finite-difference (3.17)–(3.19) for $\ell\le3$, parity, the table (3.22), Example 4.3, a Gauss–Legendre Gram matrix, SymPy checks of (3.19) and Example 4.1, and `spherical_harmonics.png`. Part 3: the Pauli algebra, Example 4.4, Theorem 3.3. Part 4: (3.36) against `expm`, Theorem 3.4 for random axes, the $2\pi$ sign. Part 5: Larmor in SI (Example 4.5) and Schrödinger vs Heisenberg vs (3.33), with `larmor_bloch.png`. Part 6: the lab-frame equation with (3.34) integrated and compared with (3.35) at five detunings (`rabi.png`). Part 7: `np.kron` addition, (3.41), and `multiplets_by_lowering` checked against SymPy's `CG`. Part 8: QM-1's Stern–Gerlach simulator for any $j$: spin 1 gives beams $\frac14, \frac12, \frac14$, the $m_x = 0$ beam through SG$z$ gives $\frac12, 0, \frac12$, spin $\frac32$ at $60°$ gives $\frac{27}{64}, \frac{27}{64}, \frac9{64}, \frac1{64}$.

**Build it yourself first:** write `jmat(j)` from (3.14) and assert the commutators, and write the lowering-plus-orthogonal-complement procedure for general $(j_1, j_2)$ (Exercise C2), before reading the file's versions. Then the "Your turn" items: spin-1 Larmor precession with its quadrupole components, $1\otimes1$ Clebsch–Gordan coefficients by hand checked against SymPy, and the resonance line shape.

## 8. Where you'll struggle

- **The ladder derivation takes three attempts.** Everyone gets lost between "the ladder must stop" and "so $m_b = -m_t$." Memorize the skeleton, not the algebra: (i) $[\hat J_z, \hat J_\pm] = \pm\hbar\hat J_\pm$, so $\hat J_\pm$ shift $m$; (ii) $m^2\le\lambda$, so there is a top and a bottom; (iii) $\hat J_\mp\hat J_\pm = \hat J^2 - \hat J_z^2\mp\hbar\hat J_z$ at each end gives $\lambda = m_t(m_t+1) = m_b(m_b-1)$; (iv) factor. Write it on a blank sheet, check against §3.2, repeat a day later.
- **The half-integer puzzle gets fudged.** Books say "for orbital angular momentum $\ell$ is an integer" without saying why. The algebra allows $\frac12$; single-valuedness of $e^{im\phi}$ forbids it on the sphere (§3.3 step 1, Exercise C1). And spin is not the electron spinning: the surface would move at $171c$ (§3.4).
- **Condon–Shortley signs will make numerical checks fail by a sign.** A $Y_\ell^m$ without the $(-1)^m$ makes $\hat L_-Y_\ell^m$ disagree in sign with (3.14), and your Clebsch–Gordan coefficients disagree with tables. Before trusting a source, evaluate its $Y_1^1$ at $\theta = \phi = \pi/2$: Condon–Shortley gives $-i\sqrt{3/8\pi}$. SciPy's `sph_harm_y` and Griffiths agree with these notes; SciPy's older `sph_harm` took $(m, \ell, \phi_{\rm azimuth}, \theta_{\rm polar})$, an argument swap that has broken many scripts.
- **The $\theta/2$.** Everyone writes $\cos\theta|{\uparrow}\rangle + e^{i\phi}\sin\theta|{\downarrow}\rangle$ once. Test: at $\theta = \pi/2$ that is $|{\downarrow}\rangle$ up to phase, not an equatorial state. The same half angle is in (3.36).
- **The wrong square root.** $\hat J_\pm$ carries $\sqrt{j(j+1) - m(m\pm1)} = \sqrt{(j\mp m)(j\pm m+1)}$. Check that it vanishes at $m = j$ for $\hat J_+$ and $m = -j$ for $\hat J_-$.
- **Addition is bookkeeping; only practice makes it fast.** Do $\frac12\otimes\frac12$, $1\otimes\frac12$ and $1\otimes1$ completely by hand once each and check against a table. Usual errors: forgetting that $\hat J_-$ acts on both factors; dropping the $\sqrt{\cdot}$ on the left side of (3.14); choosing the sign of the orthogonal complement arbitrarily instead of by Condon–Shortley.
- **Precession signs.** For $\gamma > 0$ the spin precesses clockwise about $\mathbf B$ seen from above; electrons go the other way. A rotating-frame calculation that gives nonsense usually has the drive rotating the wrong way.

## 9. Check yourself

1. Derive $[\hat L_z, \hat x] = i\hbar\hat y$ from $[\hat x_i, \hat p_j] = i\hbar\delta_{ij}$.
2. Write $\hat J_-\hat J_+$ in terms of $\hat J^2$ and $\hat J_z$, and use it to get the top-rung condition $\lambda = j(j+1)$.
3. Why is $m_b = -m_t$, and why does that allow half-integer $j$?
4. Why must orbital $m$ be an integer? What would $m = \frac12$ do to $\psi$ on a trip around the $z$-axis?
5. Write $Y_\ell^\ell$ up to normalization from $\hat L_+Y = 0$, and state the parity of $Y_\ell^m$.
6. State (3.26) and use it to show $(\hat{\mathbf n}\cdot\boldsymbol\sigma)^2 = I$.
7. Write $|\hat{\mathbf n}, +\rangle$ and prove $\langle\boldsymbol\sigma\rangle = \hat{\mathbf n}$. Why are orthogonal spin states antipodal?
8. Give the Larmor frequency of a proton at 3 T and state the Rabi formula; what is a $\pi$-pulse?
9. Why does $e^{-i2\pi\hat{\mathbf n}\cdot\boldsymbol\sigma/2} = -I$, and how was the sign observed?
10. Build the triplet and singlet from $|{\uparrow\uparrow}\rangle$, and check $(2j_1+1)(2j_2+1) = \sum(2j+1)$ for $1\otimes\frac12$.

## 10. Pair with

- **Griffiths & Schroeter (3rd ed.), §4.3–4.4.** The standard route through §3.1–3.10 with many short problems; Table 4.8 is the Clebsch–Gordan table of Example 4.7.
- **Townsend, *A Modern Approach to Quantum Mechanics*, ch. 1–3 and 9.** Spin and Stern–Gerlach first; rotation operators and $SU(2)$ more carefully than Griffiths, with the neutron-interferometry discussion.
- **Shankar, *Principles of Quantum Mechanics*, ch. 12 and 14.** Rotations generated by $\hat{\mathbf L}$ from the translation analogy; magnetic resonance worked at length.
- **Sakurai & Napolitano, *Modern Quantum Mechanics*, ch. 3.** The graduate version: $d^j_{m'm}(\beta)$, full addition theory, Wigner–Eckart. Read after this module.
- **MIT 8.05 (Zwiebach), angular momentum lectures.** Clean algebraic derivations of §3.2 and §3.10 and of spin–orbit coupling.
- **LA-12 Appendix A.** The mathematics of §3.8. Reread it before that section.

