# QM-7: Composite Systems and Entanglement: Tensor Products, Identical Particles, Bell States, EPR, Bell's Inequality, Density Matrices

**Time:** ~12 hours · **Prerequisites:** QM-5, LA-12 · **Unlocks:** QC-1, QC-2, QC-7, SM-5

> **Conventions in this module.** $\hbar$ explicit in the notes; `project.py` sets $\hbar = 1$. A two-part system has state space $\mathcal H_A\otimes\mathcal H_B$ and $|ab\rangle := |a\rangle\otimes|b\rangle$ with the **left factor being system $A$** (QC convention: qubit 0 leftmost; LA-12 §3.14). Pair indices are ordered row-major, $(i,j)\mapsto i\,d_B + j$, so `np.kron(A, B)` is $\hat A\otimes\hat B$ and `psi.reshape(dA, dB)` is the coefficient matrix $C$. Spin-$\frac12$: $|0\rangle = |{\uparrow}\rangle = (1,0)^\top$, $|1\rangle = |{\downarrow}\rangle$, $\hat{\mathbf S} = \frac\hbar2\boldsymbol\sigma$ (QM-5). Bell states $|\Phi^\pm\rangle = (|00\rangle \pm |11\rangle)/\sqrt2$, $|\Psi^\pm\rangle = (|01\rangle \pm |10\rangle)/\sqrt2$; $|\Psi^-\rangle$ is QM-5's singlet (3.41). Density operator $\hat\rho$: Hermitian, $\operatorname{Tr}\hat\rho = 1$, $\hat\rho\ge0$. Partial trace $\operatorname{Tr}_B$. Purity $\operatorname{Tr}\hat\rho^2$. Von Neumann entropy $S = -\operatorname{Tr}\hat\rho\ln\hat\rho$ in **nats** (QC-7 uses bits: divide by $\ln 2$). CHSH combination $S_{\rm CHSH} = E(a,b) - E(a,b') + E(a',b) + E(a',b')$. Bosons: wavefunction symmetric under exchange; fermions: antisymmetric.

## 0. Where this sits

You need QM-4's postulates (states, Born rule, projection, unitary evolution), QM-5's spin-$\frac12$, Pauli matrices and singlet (3.41), and LA-12's tensor products: the coefficient matrix (§3.12), the Schmidt decomposition (Theorem 3.16) and Kronecker products (§3.14). LA-12 did the mathematics; this module does the physics that made people care.

QC-1 to QC-3 run on §3.1–3.11: qubits are spins, two-qubit gates are Kronecker matrices, and QC-3 plays CHSH as a game. QC-7 makes the density operator of §3.4–3.7 the main object (channels, Kraus operators, decoherence, §3.15 here as preview). SM-5 turns §3.12–3.14's symmetrization into Bose–Einstein and Fermi–Dirac statistics. QM-8 computes helium's energies better than §3.14 does.

## 1. The problem that forced this idea

For two particles, Schrödinger's equation lives on configuration space: $\psi(\mathbf r_1, \mathbf r_2)$ is one function of six coordinates, not two functions of three. Heisenberg and Dirac (1926) saw that for identical particles it must be symmetric or antisymmetric, and that antisymmetry is Pauli's 1925 exclusion rule ("no two electrons share all four quantum numbers"), which he had extracted from atomic spectra. It explains the periodic table and why matter takes up space. Why some particles are one kind and some the other waited for Pauli's 1940 spin-statistics theorem, whose proof needs relativity.

In 1927 Landau and von Neumann introduced the density matrix, forced by two problems: real preparations are imperfect, so you have a probability distribution over states; and a part of a larger system generally has no state vector of its own. If two particles have interacted, the pair has a wavefunction and neither particle does.

In 1935 Einstein, Podolsky and Rosen (EPR) used exactly such a pair to argue that quantum mechanics is incomplete. After the particles separate, measuring one tells you the other's position or momentum, your choice, without touching it; so, they argued, the far particle must have had both values all along, which the theory cannot express. Schrödinger, the same year, named the phenomenon *Verschränkung*, entanglement, "not one but rather *the* characteristic trait of quantum mechanics." Bohr replied that the experimental arrangement is part of what a quantity means; physicists accepted his authority more than his argument, and the question went quiet. Bohm (1951) rewrote EPR with spins, the version used here, and in 1952 published a hidden-variable theory that reproduces quantum mechanics but is explicitly nonlocal.

John Bell (1964) turned the philosophy into arithmetic. Any theory in which outcomes are fixed in advance by local hidden variables obeys an inequality on correlations; quantum mechanics violates it. Clauser, Horne, Shimony and Holt (1969) gave the form used in experiments (CHSH). Freedman and Clauser (1972) and Aspect's group (1982) measured violations; in 2015 three groups (Delft, Vienna, NIST) closed the main loopholes at once. Aspect, Clauser and Zeilinger shared the 2022 Nobel Prize. Local realism is ruled out by experiment, not by taste.

The forcing problem: two particles' state space is $\mathcal H_A\otimes\mathcal H_B$ (LA-12 §3.11), almost none of its vectors are products, and the correlations in the others have no classical description. Handling them needs mathematics (Schmidt, partial trace, density operators) and physics (no-signaling, Bell).

## 2. Intuition

**A two-particle state is a matrix.** The amplitudes of $|\psi\rangle = \sum_{ij}C_{ij}|i\rangle|j\rangle$ form a $d_A\times d_B$ matrix $C$ (LA-12 §3.12). A product $|a\rangle|b\rangle$ has $C = \mathbf a\mathbf b^\top$, rank 1; everything else is **entangled**. For two qubits, rank 1 means $\det C = 0$: $|00\rangle + |01\rangle + |10\rangle + |11\rangle$ (det 0) is $|{+}\rangle|{+}\rangle$ in disguise; $|00\rangle + |11\rangle$ (det 1) is a product of nothing. Guess the test: the number of nonzero singular values, the Schmidt rank.

**The whole can be pure while the parts are random.** In $|\Phi^+\rangle = (|00\rangle + |11\rangle)/\sqrt2$, measure $A$ in $z$: 50/50, and $B$ then certainly matches. Since $|\Phi^+\rangle = (|{++}\rangle + |{--}\rangle)/\sqrt2$ (expand it), the same holds in $x$. Someone holding only $B$ sees a fair coin in every basis, and no state vector does that (each $|\hat{\mathbf n}\rangle$ is certain along $\hat{\mathbf n}$). Guess: a subsystem is described by a positive trace-1 matrix (here $\frac12 I$, the center of the Bloch ball), and its mixedness measures its entanglement with the rest.

**Collapse is real for bookkeeping, useless for signaling.** When $A$ gets "0", $B$'s description jumps to $|0\rangle$, but Bob doesn't know $A$ got "0". Averaged over $A$'s outcomes he sees the same fair coin whether $A$ measured $z$, $x$, or nothing. Guess the theorem: nothing $A$ does locally changes any probability at $B$; correlations appear only when the two lists are compared, via an ordinary signal.

**A mixture is a list of states with odds, not a superposition.** "50% $|0\rangle$, 50% $|1\rangle$" and "50% $|{+}\rangle$, 50% $|{-}\rangle$" are different recipes with identical statistics for every measurement. The superposition $|{+}\rangle$ is different: it is certain in $x$. The density matrix forgets the recipe, correctly, because the recipe is not physically accessible.

**Bell's inequality is about pre-written answers.** Suppose each particle carries a card listing its answer ($\pm1$) to each question. Alice is asked $a$ or $a'$, Bob $b$ or $b'$. For one pair of cards, $A(B - B') + A'(B + B') = \pm2$, because one of $B\mp B'$ is zero and the other $\pm2$. So its average over many pairs lies in $[-2, 2]$. Quantum mechanics predicts $2\sqrt2$. No cards written in advance can do that, unless Alice's card depends on Bob's question.

**Identical particles have no labels.** Swapping two identical electrons changes nothing measurable, and swapping twice is the identity, so the swap multiplies $\psi$ by $\pm1$. With $-1$, both in one state gives $\phi(x_1)\phi(x_2) - \phi(x_2)\phi(x_1) = 0$: exclusion in one line. Guess for positions: antisymmetric $\psi$ vanishes at $x_1 = x_2$, so fermions sit farther apart than distinguishable particles in the same orbitals, bosons closer, with no force involved.

## 3. Formal development

### 3.1 The composite-system postulate and Kronecker products

**Postulate (composite systems; QM-1's sixth postulate).** If system $A$ has state space $\mathcal H_A$ and $B$ has $\mathcal H_B$, the combined system has state space $\mathcal H_A\otimes\mathcal H_B$. If $\{|a_i\rangle\}$ and $\{|b_j\rangle\}$ are orthonormal bases, $\{|a_i\rangle\otimes|b_j\rangle\}$ is an orthonormal basis (LA-12 Theorem 3.14), so $\dim = d_Ad_B$. Ten qubits: $2^{10} = 1024$ dimensions, not $20$; this exponential growth is why QC exists.

An observable of $A$ alone, $\hat A$, acts on the composite as $\hat A\otimes\hat 1$; one of $B$ alone as $\hat 1\otimes\hat B$. On products, $(\hat A\otimes\hat B)(|a\rangle\otimes|b\rangle) = \hat A|a\rangle\otimes\hat B|b\rangle$, extended linearly (LA-12 (3.10)).

▶ **Try it first.** Show $[\hat A\otimes\hat 1, \hat 1\otimes\hat B] = 0$ using LA-12's identity $(A\otimes B)(C\otimes D) = AC\otimes BD$.

**Theorem 3.1 (Local operators commute).** For any operators $\hat A$ on $\mathcal H_A$, $\hat B$ on $\mathcal H_B$, $[\hat A\otimes\hat1, \hat1\otimes\hat B] = 0$.

*Proof.*
1. $(\hat A\otimes\hat1)(\hat1\otimes\hat B) = \hat A\hat1\otimes\hat1\hat B = \hat A\otimes\hat B$ (LA-12 Theorem 3.17(a)).
2. $(\hat1\otimes\hat B)(\hat A\otimes\hat1) = \hat1\hat A\otimes\hat B\hat1 = \hat A\otimes\hat B$.
3. Subtract. ∎

■ **What just happened.** Measurements on different subsystems are compatible (QM-4): they share an eigenbasis, their order doesn't matter, and you can measure both. This holds *however* entangled the state is; entanglement is a property of states, not operators.

**Matrices.** With row-major pair ordering, the matrix of $\hat A\otimes\hat B$ is the Kronecker product, LA-12 (3.12):
$$
\hat A\otimes\hat B = \begin{pmatrix}A_{00}B & A_{01}B\\ A_{10}B & A_{11}B\end{pmatrix}\quad(\text{for } d_A = 2), \qquad |ij\rangle \leftrightarrow \mathbf e_{i d_B + j}. \tag{3.1}
$$
So $|00\rangle, |01\rangle, |10\rangle, |11\rangle$ are $\mathbf e_0, \ldots, \mathbf e_3$, and $\sigma_z\otimes\hat1 = \operatorname{diag}(1,1,-1,-1)$ acts on the left qubit. Example 4.1 does a pair in full. **The ordering is a convention, and getting it backwards produces transposed blocks with no error message.** `np.kron(sz, I2)` is $\sigma_z$ on qubit 0; `np.kron(I2, sz)` is on qubit 1.

**Wavefunctions.** For two particles on a line, the position basis is $|x_1\rangle\otimes|x_2\rangle$ and the components $\psi(x_1, x_2) = (\langle x_1|\otimes\langle x_2|)|\psi\rangle$ are the two-particle wavefunction, with $|\psi(x_1, x_2)|^2dx_1dx_2$ the probability of finding particle 1 in $dx_1$ and particle 2 in $dx_2$. A product state is $\psi(x_1, x_2) = \phi(x_1)\chi(x_2)$. The coefficient matrix becomes the function of two variables; LA-12's rank becomes "number of terms in the expansion $\sum_k\sigma_k\phi_k(x_1)\chi_k(x_2)$." An electron with spin has state space $L^2(\mathbb R^3)\otimes\mathbb C^2$, with "wavefunction" $\psi(\mathbf r, s)$, $s\in\{\uparrow, \downarrow\}$.

### 3.2 Product states, entangled states, and the Schmidt decomposition

**Definition.** $|\psi\rangle\in\mathcal H_A\otimes\mathcal H_B$ is a **product state** if $|\psi\rangle = |\alpha\rangle\otimes|\beta\rangle$ for some $|\alpha\rangle, |\beta\rangle$; otherwise it is **entangled**.

With $|\psi\rangle = \sum_{ij}C_{ij}|i\rangle|j\rangle$, LA-12 Theorem 3.15 says: product $\iff\operatorname{rank}C = 1$. LA-12 Theorem 3.16 (the SVD $C = U\Sigma V^\dagger$ substituted into the expansion) gives

**Theorem 3.2 (Schmidt decomposition).** Every normalized $|\psi\rangle\in\mathcal H_A\otimes\mathcal H_B$ can be written
$$
|\psi\rangle = \sum_{k=1}^r\lambda_k|u_k\rangle\otimes|v_k\rangle,\qquad \lambda_k > 0,\quad\sum_k\lambda_k^2 = 1, \tag{3.2}
$$
with $\{|u_k\rangle\}$ orthonormal in $\mathcal H_A$, $\{|v_k\rangle\}$ orthonormal in $\mathcal H_B$. The $\lambda_k$ are the singular values of $C$; $|u_k\rangle$ has components $U_{ik}$, $|v_k\rangle$ has components $\overline{V_{jk}}$; the **Schmidt rank** $r = \operatorname{rank}C\le\min(d_A, d_B)$. The state is entangled iff $r\ge2$.

The derivation is LA-12 §3.13 steps 1–4; reread it now rather than have it copied here. What is new is the reading: in the Schmidt basis the two subsystems are perfectly correlated, "$A$ in $|u_k\rangle$ exactly when $B$ is in $|v_k\rangle$," with probability $\lambda_k^2$.

**The Bell states.** Their coefficient matrices are
$$
C_{\Phi^\pm} = \tfrac{1}{\sqrt2}\begin{pmatrix}1&0\\0&\pm1\end{pmatrix},\qquad C_{\Psi^\pm} = \tfrac{1}{\sqrt2}\begin{pmatrix}0&1\\\pm1&0\end{pmatrix}. \tag{3.3}
$$
Each is $\frac1{\sqrt2}$ times a unitary, so both singular values are $\frac1{\sqrt2}$: rank 2, $\lambda_1 = \lambda_2$. This equal-weight case is **maximally entangled** (§3.5 makes "maximal" precise). The four are orthonormal (Exercise A4), so they form a basis of $\mathbb C^2\otimes\mathbb C^2$, the **Bell basis**; and each is a local Pauli acting on $|\Phi^+\rangle$: $|\Phi^-\rangle = (\hat1\otimes\sigma_z)|\Phi^+\rangle$, $|\Psi^+\rangle = (\hat1\otimes\sigma_x)|\Phi^+\rangle$, $|\Psi^-\rangle = (\hat 1\otimes\sigma_x\sigma_z)|\Phi^+\rangle$. QC-4's teleportation and superdense coding use exactly this. The Schmidt form is not unique when $\lambda_1 = \lambda_2$: $|\Phi^+\rangle = (|00\rangle + |11\rangle)/\sqrt2 = (|{++}\rangle + |{--}\rangle)/\sqrt2$, which is §2's "same correlation in every basis."

### 3.3 Measuring one subsystem

Alice measures an observable $\hat A = \sum_a a\hat P_a$ on her half ($\hat P_a$ the eigenprojector for eigenvalue $a$). On the composite, she is measuring $\hat A\otimes\hat1 = \sum_a a\,(\hat P_a\otimes\hat1)$, whose eigenprojectors are $\hat P_a\otimes\hat1$.

▶ **Try it first.** Apply QM-4's Born rule and projection postulate to $\hat A\otimes\hat1$. Then, for $\hat P_a = |a\rangle\langle a|$, write the post-measurement state as $|a\rangle\otimes(\text{something in }\mathcal H_B)$.

**Derivation.**
1. *Probability.* Born rule for the composite observable:
$$
p(a) = \langle\psi|(\hat P_a\otimes\hat1)|\psi\rangle. \tag{3.4}
$$
2. *Post-measurement state.* Projection postulate:
$$
|\psi\rangle\ \longmapsto\ |\psi_a\rangle = \frac{(\hat P_a\otimes\hat1)|\psi\rangle}{\sqrt{p(a)}}. \tag{3.5}
$$
3. *Bob's conditional state.* For nondegenerate $\hat P_a = |a\rangle\langle a|$, $(\hat P_a\otimes\hat1)\sum_{ij}C_{ij}|i\rangle|j\rangle = |a\rangle\otimes\sum_j\big(\sum_i\langle a|i\rangle C_{ij}\big)|j\rangle$. So after outcome $a$ the state is the product $|a\rangle\otimes|\beta_a\rangle$ with
$$
|\beta_a\rangle = \frac{1}{\sqrt{p(a)}}\sum_j\Big(\sum_i\overline{\alpha_i}\,C_{ij}\Big)|j\rangle,\qquad \alpha_i = \langle i|a\rangle,\qquad p(a) = \|C^\top\bar{\boldsymbol\alpha}\|^2. \tag{3.6}
$$

■ **What just happened.** Measuring $A$ destroys the entanglement and leaves Bob in $|\beta_a\rangle$, which depends on Alice's basis. For $|\Phi^+\rangle$, $|\beta_a\rangle = \overline{|a\rangle}$: Alice's $|{+}\rangle$ leaves Bob in $|{+}\rangle$, her $|{+y}\rangle$ leaves him in $|{-y}\rangle$. Can Bob tell anything happened? §3.8.

### 3.4 The density operator: mixtures

A source emits $|\psi_1\rangle$ with probability $p_1$, $|\psi_2\rangle$ with $p_2$, and so on ($p_i\ge0$, $\sum p_i = 1$; the $|\psi_i\rangle$ normalized, not necessarily orthogonal). This is an **ensemble** $\{p_i, |\psi_i\rangle\}$.

▶ **Try it first.** Write $\langle\hat Q\rangle$ for this source as an average of averages, then use $\langle\psi|\hat Q|\psi\rangle = \operatorname{Tr}(\hat Q|\psi\rangle\langle\psi|)$ (cyclicity of the trace, LA-7) to pull the sum inside.

**Derivation.**
1. Classical probability over which state, quantum expectation within each: $\langle\hat Q\rangle = \sum_ip_i\langle\psi_i|\hat Q|\psi_i\rangle$.
2. For any $|\psi\rangle$, $\langle\psi|\hat Q|\psi\rangle = \sum_n\langle n|\hat Q|\psi\rangle\langle\psi|n\rangle$ (insert $\sum_n|n\rangle\langle n|$) $= \operatorname{Tr}(\hat Q|\psi\rangle\langle\psi|)$.
3. Linearity of the trace:
$$
\langle\hat Q\rangle = \operatorname{Tr}(\hat\rho\hat Q),\qquad \hat\rho := \sum_ip_i|\psi_i\rangle\langle\psi_i|. \tag{3.7}
$$
4. Outcome probabilities are expectations of projectors: $p(q) = \operatorname{Tr}(\hat\rho\hat P_q)$. After outcome $q$, each branch is projected by (3.5), so $\hat\rho\mapsto\hat P_q\hat\rho\hat P_q/p(q)$.

■ **What just happened.** Every prediction depends on the ensemble only through $\hat\rho$. That makes $\hat\rho$, the **density operator**, the state. A pure state $|\psi\rangle$ is the special case $\hat\rho = |\psi\rangle\langle\psi|$, which also removes the global-phase redundancy of kets.

**Theorem 3.3 (Properties).** Any $\hat\rho$ of the form (3.7) is (i) Hermitian, (ii) of trace 1, (iii) positive semidefinite: $\langle\phi|\hat\rho|\phi\rangle\ge0$ for all $|\phi\rangle$.

*Proof.* (i) Each $|\psi_i\rangle\langle\psi_i|$ is Hermitian and the $p_i$ are real. (ii) $\operatorname{Tr}|\psi_i\rangle\langle\psi_i| = \langle\psi_i|\psi_i\rangle = 1$, so $\operatorname{Tr}\hat\rho = \sum p_i = 1$. (iii) $\langle\phi|\hat\rho|\phi\rangle = \sum_ip_i|\langle\phi|\psi_i\rangle|^2\ge0$. ∎

Conversely any operator with (i)–(iii) is a density operator: diagonalize it (spectral theorem, LA-8), $\hat\rho = \sum_k\mu_k|k\rangle\langle k|$; (iii) makes $\mu_k\ge0$ and (ii) makes them sum to 1, so it is the ensemble $\{\mu_k, |k\rangle\}$. Its diagonal entries in any basis are probabilities; the off-diagonal entries $\rho_{mn}$ are the **coherences**.

**Theorem 3.4 (Purity).** For a density operator, $\operatorname{Tr}\hat\rho^2\le1$, and the following are equivalent: (a) $\hat\rho = |\psi\rangle\langle\psi|$ for some $|\psi\rangle$; (b) $\hat\rho^2 = \hat\rho$; (c) $\operatorname{Tr}\hat\rho^2 = 1$.

*Proof.* Write $\hat\rho = \sum_k\mu_k|k\rangle\langle k|$ with $\mu_k\in[0,1]$, $\sum\mu_k = 1$. Then $\hat\rho^2 = \sum\mu_k^2|k\rangle\langle k|$ and $\operatorname{Tr}\hat\rho^2 = \sum\mu_k^2\le\sum\mu_k = 1$, since $\mu_k^2\le\mu_k$ on $[0,1]$.
(a)⇒(b): $(|\psi\rangle\langle\psi|)^2 = |\psi\rangle\langle\psi|\psi\rangle\langle\psi| = |\psi\rangle\langle\psi|$.
(b)⇒(c): $\operatorname{Tr}\hat\rho^2 = \operatorname{Tr}\hat\rho = 1$.
(c)⇒(a): $\sum_k(\mu_k - \mu_k^2) = 1 - 1 = 0$, a sum of nonnegative terms, so each $\mu_k\in\{0,1\}$; with $\sum\mu_k = 1$ exactly one is 1. So $\hat\rho = |k\rangle\langle k|$. ∎

■ **What just happened.** Purity $\operatorname{Tr}\hat\rho^2$ runs from 1 (pure) down to $1/d$ (the **maximally mixed** state $\hat1/d$, all $\mu_k = 1/d$). A single number tells you whether there is any classical uncertainty left in the description.

### 3.5 The reduced density operator and the partial trace

Now the second motivation. The pair is in a pure entangled state $|\psi\rangle$ and Alice only ever looks at $A$. She needs an object $\hat\rho_A$ on $\mathcal H_A$ that gives all her predictions.

▶ **Try it first.** Demand $\langle\psi|(\hat A\otimes\hat1)|\psi\rangle = \operatorname{Tr}(\hat\rho_A\hat A)$ for every $\hat A$. Expand the left side in components with $C$, and read off $\hat\rho_A$.

**Derivation.**
1. $\langle\psi|(\hat A\otimes\hat1)|\psi\rangle = \sum_{kl}\sum_{ij}\overline{C_{kl}}C_{ij}\langle k|\hat A|i\rangle\langle l|j\rangle = \sum_{ik}A_{ki}\sum_j C_{ij}\overline{C_{kj}}$.
2. $\sum_jC_{ij}\overline{C_{kj}} = (CC^\dagger)_{ik}$, so the left side is $\sum_{ik}A_{ki}(CC^\dagger)_{ik} = \operatorname{Tr}(\hat A\,CC^\dagger)$.
3. This holds for every $\hat A$, and the trace pairing $(X, Y)\mapsto\operatorname{Tr}(XY)$ is nondegenerate, so $\hat\rho_A$ is unique:
$$
\hat\rho_A = CC^\dagger,\qquad \hat\rho_B = C^\top\overline C\quad(\text{same argument with }\hat1\otimes\hat B). \tag{3.8}
$$
4. *The general operation.* For any operator $\hat\rho$ on $\mathcal H_A\otimes\mathcal H_B$ with matrix elements $\rho_{ij,kl} = \langle ij|\hat\rho|kl\rangle$, define the **partial trace**
$$
(\operatorname{Tr}_B\hat\rho)_{ik} = \sum_j\rho_{ij,kj},\qquad(\operatorname{Tr}_A\hat\rho)_{jl} = \sum_i\rho_{ij,il}. \tag{3.9}
$$
For $\hat\rho = |\psi\rangle\langle\psi|$, $\rho_{ij,kl} = C_{ij}\overline{C_{kl}}$, and (3.9) gives back (3.8). The same step 1 for a mixed $\hat\rho$ shows $\operatorname{Tr}[\hat\rho(\hat A\otimes\hat1)] = \operatorname{Tr}[(\operatorname{Tr}_B\hat\rho)\hat A]$ in general. On products, $\operatorname{Tr}_B(\hat X\otimes\hat Y) = \hat X\operatorname{Tr}\hat Y$.

■ **What just happened.** The **reduced density operator** $\hat\rho_A = \operatorname{Tr}_B\hat\rho$ is forced by asking for Alice's predictions. In code: reshape the $d_Ad_B\times d_Ad_B$ matrix to a four-index array $\rho[i,j,k,l]$ and sum the two $B$ indices together, `np.einsum('ijkj->ik', rho.reshape(dA, dB, dA, dB))` (LA-12's einsum). For a pure state it is even simpler: `C = psi.reshape(dA, dB); rhoA = C @ C.conj().T`. It is a legitimate density operator: Hermitian ($(CC^\dagger)^\dagger = CC^\dagger$), trace $\|C\|_F^2 = 1$, and positive, $\langle\phi|CC^\dagger|\phi\rangle = \|C^\dagger\phi\|^2\ge0$.

**Theorem 3.5 (Bell pairs have maximally mixed halves).** For each Bell state, $\hat\rho_A = \hat\rho_B = \frac12\hat1$.

*Proof.* By (3.3) $C = W/\sqrt2$ with $W$ unitary (and real), so $CC^\dagger = WW^\dagger/2 = \hat1/2$ and $C^\top\overline C = W^\top\overline W/2 = \overline{W^\dagger W}/2 = \hat1/2$. ∎

The whole is pure and each part has the smallest possible purity, $\frac12$. Classically, a certain joint distribution has certain marginals; quantum mechanically you can know everything about the pair and nothing about either member.

**Theorem 3.6 (Entanglement entropy).** For a pure $|\psi\rangle$ with Schmidt coefficients $\lambda_k$,
$$
S(\hat\rho_A) = S(\hat\rho_B) = -\sum_{k}\lambda_k^2\ln\lambda_k^2. \tag{3.10}
$$
It is $0$ iff $|\psi\rangle$ is a product, and $\le\ln d$ with $d = \min(d_A, d_B)$, with equality iff all $d$ Schmidt coefficients equal $1/\sqrt d$.

*Proof.*
1. SVD: $C = U\Sigma V^\dagger$, so $\hat\rho_A = CC^\dagger = U\Sigma^2U^\dagger$; its eigenvalues are $\lambda_k^2$ (and zeros). Then $S = -\operatorname{Tr}\hat\rho_A\ln\hat\rho_A = -\sum_k\lambda_k^2\ln\lambda_k^2$ (with $0\ln0 = 0$).
2. $\hat\rho_B = C^\top\overline C = \overline{C^\dagger C} = \overline V\Sigma^2V^\top$: the same nonzero eigenvalues. So $S(\hat\rho_A) = S(\hat\rho_B)$, even when $d_A\neq d_B$.
3. With $p_k = \lambda_k^2$ a probability vector on $r$ entries, $S = -\sum p_k\ln p_k = 0$ iff one $p_k = 1$, i.e. $r = 1$, a product (Theorem 3.2).
4. Upper bound: $\ln$ is concave, so by Jensen $\sum_kp_k\ln(1/p_k)\le\ln\sum_kp_k(1/p_k) = \ln r\le\ln d$, with equality iff all $1/p_k$ are equal and $r = d$. ∎

■ **What just happened.** For pure global states, the mixedness of a part measures the entanglement: $0$ for products, $\ln 2$ (one bit) for a Bell pair. **Warning:** for mixed global states $S(\hat\rho_A)$ is no longer an entanglement measure ($\frac14\hat1\otimes\hat1$ has $S(\hat\rho_A) = \ln2$ and no entanglement). QC-7 handles that case.

### 3.6 The Bloch ball

▶ **Try it first.** Expand a $2\times2$ density matrix in the basis $\{\hat1, \sigma_x, \sigma_y, \sigma_z\}$ of Hermitian matrices. Which coefficient does trace 1 fix? What does positivity say about the rest, given QM-5's $(\mathbf r\cdot\boldsymbol\sigma)^2 = |\mathbf r|^2\hat1$?

**Derivation.**
1. Any $2\times2$ Hermitian matrix is a real combination $c_0\hat1 + \mathbf c\cdot\boldsymbol\sigma$ (four real parameters, four independent Hermitian matrices). $\operatorname{Tr}\sigma_i = 0$ and $\operatorname{Tr}\hat1 = 2$, so trace 1 forces $c_0 = \frac12$. Write
$$
\hat\rho = \tfrac12(\hat1 + \mathbf r\cdot\boldsymbol\sigma),\qquad \mathbf r\in\mathbb R^3,\qquad r_i = \operatorname{Tr}(\hat\rho\sigma_i) = \langle\sigma_i\rangle. \tag{3.11}
$$
(The last equality uses $\operatorname{Tr}(\sigma_i\sigma_j) = 2\delta_{ij}$.)
2. By QM-5's (3.26), $(\mathbf r\cdot\boldsymbol\sigma)^2 = |\mathbf r|^2\hat1$ and $\operatorname{Tr}(\mathbf r\cdot\boldsymbol\sigma) = 0$, so $\mathbf r\cdot\boldsymbol\sigma$ has eigenvalues $\pm|\mathbf r|$, and $\hat\rho$ has eigenvalues $\frac12(1\pm|\mathbf r|)$.
3. Positivity $\iff\frac12(1 - |\mathbf r|)\ge0\iff|\mathbf r|\le1$.
4. Purity: $\operatorname{Tr}\hat\rho^2 = \frac14\operatorname{Tr}(\hat1 + 2\mathbf r\cdot\boldsymbol\sigma + |\mathbf r|^2\hat1) = \frac12(1 + |\mathbf r|^2)$. By Theorem 3.4, pure iff $|\mathbf r| = 1$.

■ **What just happened.** Qubit states are the points of the solid unit ball. The surface is QM-5's Bloch sphere of pure states (Theorem 3.3 there: $|\hat{\mathbf n}\rangle\langle\hat{\mathbf n}| = \frac12(\hat1 + \hat{\mathbf n}\cdot\boldsymbol\sigma)$); the interior holds mixtures, with the length of the **Bloch vector** $\mathbf r$ measuring how much is known; the center is $\frac12\hat1$. A mixture $\sum p_i\hat\rho_i$ has Bloch vector $\sum p_i\mathbf r_i$: mixing is taking centers of mass. Half of a Bell pair sits at the center.

### 3.7 The von Neumann equation

▶ **Try it first.** Differentiate $\hat\rho(t) = \sum_ip_i|\psi_i(t)\rangle\langle\psi_i(t)|$ using $i\hbar\,\partial_t|\psi\rangle = \hat H|\psi\rangle$ and its adjoint.

**Derivation.**
1. $\partial_t|\psi_i\rangle = -\frac i\hbar\hat H|\psi_i\rangle$; taking the adjoint ($\hat H$ Hermitian), $\partial_t\langle\psi_i| = \frac i\hbar\langle\psi_i|\hat H$.
2. Product rule: $\partial_t(|\psi_i\rangle\langle\psi_i|) = -\frac i\hbar\big(\hat H|\psi_i\rangle\langle\psi_i| - |\psi_i\rangle\langle\psi_i|\hat H\big)$.
3. Sum with the time-independent weights $p_i$:
$$
i\hbar\frac{d\hat\rho}{dt} = [\hat H, \hat\rho],\qquad \hat\rho(t) = \hat U(t)\hat\rho(0)\hat U(t)^\dagger,\quad \hat U = e^{-i\hat Ht/\hbar}. \tag{3.12}
$$

■ **What just happened.** This is CM-8's Liouville equation $\partial_t\rho = \{H, \rho\}$ with Dirac's replacement $\{\cdot,\cdot\}\to\frac{1}{i\hbar}[\cdot,\cdot]$ (QM-1 §3.11): the dictionary row "phase-space density ↔ density operator." Note the sign is opposite to Heisenberg's $d\hat A/dt = \frac i\hbar[\hat H, \hat A]$ (QM-4); states and observables move in opposite directions so that $\operatorname{Tr}(\hat\rho\hat A)$ comes out the same in both pictures. Unitary evolution preserves eigenvalues, hence purity and entropy: closed systems never get more mixed. A *subsystem* of an interacting pair obeys no such equation with $\hat H_A$ alone; its purity can fall. That is §3.15, and QC-7's master equations.

### 3.8 The no-signaling theorem

Alice either does nothing, or measures in one basis, or another, or applies any unitary to her qubit. Bob, far away, measures his. Can his statistics depend on her choice?

▶ **Try it first.** Bob doesn't know Alice's outcome, so his state is the average of the branches $(\hat P_a\otimes\hat1)\hat\rho(\hat P_a\otimes\hat1)$. Take $\operatorname{Tr}_A$. Move $\hat P_a$ around inside a trace over $A$.

**Lemma 3.7.** For operators $\hat X, \hat Y$ on $\mathcal H_A$ and any $\hat\rho$ on $\mathcal H_A\otimes\mathcal H_B$: $\operatorname{Tr}_A[(\hat X\otimes\hat1)\hat\rho(\hat Y\otimes\hat1)] = \operatorname{Tr}_A[(\hat Y\hat X\otimes\hat1)\hat\rho]$.

*Proof.* In components, $[(\hat X\otimes\hat1)\hat\rho(\hat Y\otimes\hat1)]_{ij,kl} = \sum_{m,n}X_{im}\rho_{mj,nl}Y_{nk}$. Set $k = i$ and sum over $i$ (3.9): $\sum_{i,m,n}Y_{ni}X_{im}\rho_{mj,nl} = \sum_{m,n}(YX)_{nm}\rho_{mj,nl}$, which is $\operatorname{Tr}_A[(\hat Y\hat X\otimes\hat1)\hat\rho]_{jl}$. ∎

**Theorem 3.8 (No-signaling).** Let Alice perform a projective measurement $\{\hat P_a\}$ on $A$ (outcome unknown to Bob), or a unitary $\hat U_A$. Then Bob's reduced state, and hence the probability of every outcome of every measurement Bob can make, is unchanged.

*Proof.*
1. Measurement: Bob's state after Alice's measurement, averaged over her outcomes with weights $p(a)$, is $\hat\rho_B' = \sum_ap(a)\operatorname{Tr}_A\big[(\hat P_a\otimes\hat1)\hat\rho(\hat P_a\otimes\hat1)/p(a)\big]$.
2. Lemma 3.7 with $\hat X = \hat Y = \hat P_a$ and $\hat P_a^2 = \hat P_a$: $\hat\rho_B' = \sum_a\operatorname{Tr}_A[(\hat P_a\otimes\hat1)\hat\rho] = \operatorname{Tr}_A[(\textstyle\sum_a\hat P_a\otimes\hat1)\hat\rho] = \operatorname{Tr}_A\hat\rho = \hat\rho_B$, using completeness $\sum_a\hat P_a = \hat1$.
3. Unitary: Lemma 3.7 with $\hat X = \hat U_A$, $\hat Y = \hat U_A^\dagger$ gives $\operatorname{Tr}_A[(\hat U_A^\dagger\hat U_A\otimes\hat1)\hat\rho] = \hat\rho_B$.
4. Bob's probabilities are $\operatorname{Tr}(\hat\rho_B'\hat Q_b)$ (§3.5), the same as before. ∎

■ **What just happened.** Collapse changes Bob's description *conditional on Alice's result*; averaged over results, which is all Bob has, nothing changes. Alice's choice prepares different *ensembles* at Bob (for $|\Phi^+\rangle$: $\{|0\rangle, |1\rangle\}$ if she measures $z$, $\{|{+}\rangle, |{-}\rangle\}$ if $x$) with the same $\hat\rho_B$. The correlations show only when the lists are compared over a classical channel. Completeness plus cyclicity do the work; a nonlinear quantum mechanics, or a way for Bob to identify ensembles, would allow signaling.

### 3.9 What a mixed state is not

**Not a superposition.** $|{+}\rangle = (|0\rangle + |1\rangle)/\sqrt2$ and the 50/50 mixture of $|0\rangle$ and $|1\rangle$ give the same $z$ statistics. Their density matrices are
$$
|{+}\rangle\langle{+}| = \tfrac12\begin{pmatrix}1&1\\1&1\end{pmatrix}\ (\mathbf r = \hat{\mathbf x}),\qquad \tfrac12\hat 1\ (\mathbf r = 0).
$$
The coherences (off-diagonal $\frac12$) are the difference: measure $\sigma_x$ and the superposition gives $+1$ every time, the mixture 50/50. Interference lives in the coherences.

**Not a unique recipe.** $\frac12|0\rangle\langle0| + \frac12|1\rangle\langle1| = \frac12|{+}\rangle\langle{+}| + \frac12|{-}\rangle\langle{-}| = \frac12\hat1$ (Example 4.5). Every observable gives the same statistics for both recipes, because every prediction is $\operatorname{Tr}(\hat\rho\hat Q)$. In Bloch language: any set of points with center of mass at $\mathbf r$ is a recipe for $\hat\rho(\mathbf r)$. The general statement (Hughston–Jozsa–Wootters; Nielsen & Chuang Theorem 2.6) is that two ensembles $\{\sqrt{p_i}|\psi_i\rangle\}$ and $\{\sqrt{q_j}|\phi_j\rangle\}$ give the same $\hat\rho$ iff they are related by a unitary matrix. The question "which mixture is it really?" has no experimental answer; if it had, Alice could use §3.8's two ensembles to signal. So $\hat\rho$ is the complete description, and the ensemble ambiguity is required for consistency with relativity.

**Not always ignorance.** Half of a Bell pair has $\hat\rho_B = \frac12\hat1$, but Bob's qubit has no state vector of its own to be ignorant of. This *improper* mixture is mathematically identical to a coin-flip mixture; §3.15 shows this is how the classical world arises.

### 3.10 EPR, in Bohm's version

Take a spin singlet $|\Psi^-\rangle = (|{\uparrow\downarrow}\rangle - |{\downarrow\uparrow}\rangle)/\sqrt2$ and separate the two spins. QM-5 showed the singlet is rotation-invariant, so it has the same form along every axis: $|\Psi^-\rangle = (|{+\hat{\mathbf n}}, {-\hat{\mathbf n}}\rangle - |{-\hat{\mathbf n}}, {+\hat{\mathbf n}}\rangle)/\sqrt2$ up to phase. So measuring $A$ along any $\hat{\mathbf n}$ predicts $B$'s result along $\hat{\mathbf n}$ with certainty.

**The EPR argument.** (1) *Reality criterion:* "If, without in any way disturbing a system, we can predict with certainty the value of a physical quantity, then there exists an element of physical reality corresponding to this physical quantity." (2) *Locality:* measuring $A$ does not disturb $B$. (3) Alice can predict $B$'s $S_z$ with certainty by measuring her $S_z$, or $B$'s $S_x$ by measuring her $S_x$. By (1) and (2), $B$'s $S_z$ and $S_x$ are both elements of reality, whatever Alice does. (4) Quantum mechanics has no state with definite $S_z$ and $S_x$. So it is incomplete.

**Bohr's response** (1935), honestly summarized: the "without in any way disturbing" clause is ambiguous; there is no mechanical disturbance, but the choice of measurement on $A$ defines which predictions about $B$ make sense at all, so the two predictions can't be combined. The reply denies EPR's premise about what "physical reality" means rather than refuting the logic. Most physicists sided with Bohr without being able to say what experiment would decide. Bell supplied one.

**The hidden-variable program.** Supplement $|\psi\rangle$ by variables $\lambda$ that fix the outcomes, so quantum probabilities are ignorance of $\lambda$. Von Neumann's 1932 "impossibility proof" assumed an unjustified linearity (Grete Hermann 1935, Bell 1966); Bohm's 1952 theory showed hidden variables are possible, at the price of nonlocality. Bell asked whether that price is always required.

### 3.11 Bell's inequality (CHSH)

**Setup.** Alice chooses setting $a$ or $a'$, Bob $b$ or $b'$; each gets an outcome $\pm1$. The **correlation** is $E(a,b) = \langle\text{product of outcomes}\rangle$. **Local realism** means three things: (i) a variable $\lambda$ with distribution $\rho(\lambda)$ determines outcomes $A(a,\lambda), B(b,\lambda)\in\{\pm1\}$ (realism); (ii) $A$ does not depend on $b$ nor $B$ on $a$ (locality); (iii) $\rho(\lambda)$ does not depend on the settings (free choice, "measurement independence"). Then $E(a,b) = \int\rho(\lambda)A(a,\lambda)B(b,\lambda)\,d\lambda$. (Stochastic local models reduce to this by enlarging $\lambda$.)

▶ **Try it first.** For fixed $\lambda$, with $A = A(a,\lambda)$ etc., find the possible values of $A(B - B') + A'(B + B')$.

**Theorem 3.9 (CHSH inequality).** Under local realism, $|S_{\rm CHSH}|\le2$, where
$$
S_{\rm CHSH} = E(a,b) - E(a,b') + E(a',b) + E(a',b'). \tag{3.13}
$$

*Proof.*
1. For each $\lambda$, $B, B'\in\{\pm1\}$, so either $B = B'$ (then $B - B' = 0$, $B + B' = \pm2$) or $B = -B'$ (then $B + B' = 0$, $B - B' = \pm2$).
2. So $s(\lambda) := A(B - B') + A'(B + B') = \pm2$, because $A, A' = \pm1$. Equivalently: of the 16 assignments of $(A, A', B, B')$, all give $s = \pm2$ (the project enumerates them).
3. By (ii), the same $A(a,\lambda)$ appears in $E(a,b)$ and $E(a,b')$, and by (iii) all four correlations average over the same $\rho(\lambda)$. So $S_{\rm CHSH} = \int\rho(\lambda)s(\lambda)\,d\lambda$, an average of numbers in $\{-2, 2\}$: $|S_{\rm CHSH}|\le2$. ∎

**The quantum prediction.** Spin measurement along unit vector $\mathbf a$ is the observable $\mathbf a\cdot\boldsymbol\sigma$, eigenvalues $\pm1$.

▶ **Try it first.** Use $(\sigma_i\otimes\hat1 + \hat1\otimes\sigma_i)|\Psi^-\rangle = 0$ (the singlet has total spin 0; QM-5 §3.10) to move Bob's operator onto Alice's side.

**Derivation of $E(\mathbf a, \mathbf b) = -\mathbf a\cdot\mathbf b$.**
1. $E(\mathbf a,\mathbf b) = \langle\Psi^-|(\mathbf a\cdot\boldsymbol\sigma)\otimes(\mathbf b\cdot\boldsymbol\sigma)|\Psi^-\rangle = \sum_{ij}a_ib_j\langle\Psi^-|\sigma_i\otimes\sigma_j|\Psi^-\rangle$.
2. Total spin zero: $(\hat1\otimes\sigma_j)|\Psi^-\rangle = -(\sigma_j\otimes\hat1)|\Psi^-\rangle$. (Check for $j = z$: $\sigma_z$ on the second factor maps $|01\rangle - |10\rangle\mapsto -|01\rangle - |10\rangle$, and on the first factor maps it to $|01\rangle + |10\rangle$.)
3. So $(\sigma_i\otimes\sigma_j)|\Psi^-\rangle = (\sigma_i\otimes\hat1)(\hat1\otimes\sigma_j)|\Psi^-\rangle = -(\sigma_i\sigma_j\otimes\hat1)|\Psi^-\rangle$.
4. $\langle\Psi^-|(\sigma_i\sigma_j\otimes\hat1)|\Psi^-\rangle = \operatorname{Tr}(\hat\rho_A\sigma_i\sigma_j) = \frac12\operatorname{Tr}(\sigma_i\sigma_j) = \delta_{ij}$ by §3.5 and Theorem 3.5.
5. Hence $\langle\sigma_i\otimes\sigma_j\rangle = -\delta_{ij}$ and
$$
E(\mathbf a,\mathbf b) = -\mathbf a\cdot\mathbf b = -\cos\theta_{ab}. \tag{3.14}
$$

Example 4.6 redoes this with explicit $4\times4$ matrices. Each side alone is a fair coin (Theorem 3.5), whatever the settings: no-signaling. Only the correlation depends on the relative angle.

**Theorem 3.10 (The quantum maximum; Tsirelson's bound).** With (3.14), $\max|S_{\rm CHSH}| = 2\sqrt2$, attained with $\mathbf a\perp\mathbf a'$, $\mathbf b\perp\mathbf b'$ and each at $45°$ to its neighbours.

*Proof.*
1. $|S_{\rm CHSH}| = |\mathbf a\cdot(\mathbf b - \mathbf b') + \mathbf a'\cdot(\mathbf b + \mathbf b')|$.
2. For fixed $\mathbf b, \mathbf b'$, the maximum over unit $\mathbf a, \mathbf a'$ is reached with $\mathbf a\parallel(\mathbf b - \mathbf b')$, $\mathbf a'\parallel(\mathbf b + \mathbf b')$: $|\mathbf b - \mathbf b'| + |\mathbf b + \mathbf b'|$.
3. With $\gamma$ the angle between $\mathbf b, \mathbf b'$: $|\mathbf b\mp\mathbf b'| = \sqrt{2\mp2\cos\gamma}$, i.e. $2\sin\frac\gamma2$ and $2\cos\frac\gamma2$. Their sum $2(\sin\frac\gamma2 + \cos\frac\gamma2) = 2\sqrt2\sin(\frac\gamma2 + \frac\pi4)$ is maximal at $\gamma = 90°$: $2\sqrt2$.
4. Then $\mathbf a\parallel\mathbf b - \mathbf b'$ and $\mathbf a'\parallel\mathbf b + \mathbf b'$ are perpendicular and each at $45°$ to $\mathbf b$ and $\mathbf b'$. ∎

That no state and no $\pm1$ observables beat $2\sqrt2$ is Tsirelson's theorem (1980); Exercise C1 proves it from $\hat{\mathcal C}^2 = 4\cdot\hat1 - [\hat A, \hat A']\otimes[\hat B, \hat B']$. Ranking: local realism $\le2$, quantum $\le2\sqrt2\approx2.83$, no-signaling alone $\le4$ (Popescu–Rohrlich boxes). Why nature stops at $2\sqrt2$ is open.

**Photons.** Aspect-type experiments use polarization-entangled photons in $|\Phi^+\rangle = (|HH\rangle + |VV\rangle)/\sqrt2$. Polarization is a two-state system whose Bloch sphere (the Poincaré sphere) doubles angles: $H$ and $V$ are $90°$ apart in space but antipodal on the sphere. So a polarizer at angle $\chi$ is a Bloch axis at $2\chi$, and $E(\chi_a, \chi_b) = \cos2(\chi_a - \chi_b)$. The optimal polarizer angles are $a = 0°$, $a' = 45°$, $b = 22.5°$, $b' = 67.5°$, half the spin angles. Example 4.6 does both.

**Experiments, with numbers.**
- *Freedman and Clauser (1972)*, calcium cascade photons, fixed polarizers: violation by about 6 standard deviations.
- *Aspect, Grangier, Roger (1982)*: $S_{\rm CHSH} = 2.697\pm0.015$, against the quantum prediction $2.70\pm0.05$ for their setup's efficiency and angles, more than 40 standard deviations above 2. Aspect, Dalibard and Roger (1982) switched the analyzer settings during the photons' flight, addressing locality.
- *2015 loophole-free tests.* Earlier tests left the **locality loophole** (a light-speed signal could carry the settings) or the **detection loophole** (most photons undetected, the detected subset possibly unrepresentative). Hensen et al. (Delft) entangled electron spins in two diamonds 1.3 km apart: $S_{\rm CHSH} = 2.42\pm0.20$ from 245 trials, $p\approx0.039$. Giustina et al. (Vienna) and Shalm et al. (NIST) used high-efficiency photon detectors and an Eberhard-type inequality, with $p$-values about $3.7\times10^{-31}$ and $2.3\times10^{-7}$.

**What is ruled out, and what isn't.** Ruled out: every theory satisfying (i)–(iii), the class EPR hoped for. Not ruled out: *nonlocal* hidden variables (Bohm's theory reproduces every prediction, its nonlocality hidden by no-signaling); *superdeterminism* (violating (iii): settings and particles correlated from the past, which would make every randomized experiment suspect); *many-worlds* (no single outcome, so (i) is not even formulated); and views that deny (i) while keeping a modified locality (QBism, relational QM). This module's position is only what the experiments say: **no local hidden-variable theory in the sense of (i)–(iii) describes nature.** Which remaining option is right is a real, undecided question.

**Mermin's device** (*Physics Today* 1985, required reading) packages Bell's argument as two detectors with three switch positions and red/green lights: Example 4.7.

**GHZ.** Three particles make the contradiction all-or-nothing instead of statistical (Greenberger, Horne, Zeilinger 1989; Mermin 1990). In $|{\rm GHZ}\rangle = (|000\rangle + |111\rangle)/\sqrt2$, the observables $X_1X_2X_3$, $X_1Y_2Y_3$, $Y_1X_2Y_3$, $Y_1Y_2X_3$ commute and the state is an eigenstate of all four, with eigenvalues $+1, -1, -1, -1$ (the project asserts it; $\sigma_x = X$, $\sigma_y = Y$). A local-realist assignment $x_i, y_i = \pm1$ must satisfy $x_1y_2y_3 = y_1x_2y_3 = y_1y_2x_3 = -1$; multiplying these, each $y_i^2 = 1$ drops out and $x_1x_2x_3 = -1$, contradicting $X_1X_2X_3 = +1$. One run of each kind refutes it. QC-3 recasts CHSH and GHZ as nonlocal games.

### 3.12 Identical particles

**The exchange operator.** On two-particle states define $\hat P_{12}$ by $\hat P_{12}(|\phi\rangle\otimes|\chi\rangle) = |\chi\rangle\otimes|\phi\rangle$, so $(\hat P_{12}\psi)(x_1, x_2) = \psi(x_2, x_1)$. Then $\hat P_{12}^2 = \hat 1$, so its eigenvalues satisfy $\lambda^2 = 1$: $\lambda = \pm1$. It is Hermitian and unitary. For identical particles every observable treats the two the same, $\hat P_{12}\hat Q\hat P_{12} = \hat Q$, so $[\hat Q, \hat P_{12}] = 0$; in particular $[\hat H, \hat P_{12}] = 0$ and exchange symmetry is conserved.

**Symmetrization postulate.** States of identical particles are eigenstates of every $\hat P_{ij}$: **bosons** with eigenvalue $+1$ (symmetric), **fermions** with $-1$ (antisymmetric). Why it's a postulate: indistinguishability alone requires only $|\psi(x_2, x_1)|^2 = |\psi(x_1, x_2)|^2$, i.e. $\psi(x_2, x_1) = e^{i\theta}\psi(x_1, x_2)$. In three dimensions, exchanging twice is continuously deformable to doing nothing, which forces $e^{2i\theta} = 1$, so $\pm1$. In two dimensions it is not, and **anyons** with other phases exist as quasiparticles (the fractional quantum Hall effect).

**Spin-statistics theorem (stated).** Particles of integer spin are bosons; particles of half-integer spin are fermions. Electrons, protons, neutrons, quarks: fermions. Photons, gluons, the Higgs, ${}^4$He atoms: bosons. Composites count constituents: an odd number of fermions makes a fermion. The proof (Fierz 1939, Pauli 1940) needs relativistic quantum field theory: with the wrong statistics, either the energy is unbounded below or field operators at spacelike separation fail to commute (causality fails). Nonrelativistic quantum mechanics cannot derive it; we use it as a law.

▶ **Try it first.** Build a normalized antisymmetric state from two orthonormal one-particle states $\phi_a, \phi_b$. Then set $b = a$.

**Theorem 3.11 (Pauli exclusion).** Two identical fermions cannot occupy the same one-particle state.

*Proof.*
1. The antisymmetric combination is $\psi_-(x_1, x_2) = \frac{1}{\sqrt2}[\phi_a(x_1)\phi_b(x_2) - \phi_b(x_1)\phi_a(x_2)]$; for $a\neq b$, orthonormality makes its norm 1. Check $\hat P_{12}\psi_- = -\psi_-$.
2. For $b = a$ it is $\phi_a(x_1)\phi_a(x_2) - \phi_a(x_1)\phi_a(x_2) = 0$, not a state.
3. More generally, any antisymmetric $\psi$ obeys $\psi(x, x) = -\psi(x, x)$, so $\psi(x, x) = 0$: two fermions with the same spin are never at the same point. ∎

**Slater determinants.** For $N$ fermions in orthonormal orbitals $\phi_1, \ldots, \phi_N$ (here "orbital" includes spin),
$$
\psi(x_1, \ldots, x_N) = \frac{1}{\sqrt{N!}}\det\begin{pmatrix}\phi_1(x_1) & \phi_2(x_1) & \cdots & \phi_N(x_1)\\ \vdots & & & \vdots\\ \phi_1(x_N) & \phi_2(x_N) & \cdots & \phi_N(x_N)\end{pmatrix}. \tag{3.15}
$$
LA-6 does all the work: swapping particles $i$ and $j$ swaps two rows, which flips the sign of the determinant (antisymmetry); two equal orbitals make two equal columns, so the determinant is zero (Pauli); the $N!$ terms of the Leibniz expansion are mutually orthogonal products, so $1/\sqrt{N!}$ normalizes. Bosons use the **permanent** (same expansion, all signs $+$), plus a normalization correction when orbitals repeat.

### 3.13 Exchange forces

Two noninteracting particles on a line occupy orthonormal orbitals $\phi_a, \phi_b$. Three cases:
$$
\psi_D = \phi_a(x_1)\phi_b(x_2),\qquad \psi_\pm = \tfrac{1}{\sqrt2}[\phi_a(x_1)\phi_b(x_2)\pm\phi_b(x_1)\phi_a(x_2)] \tag{3.16}
$$
(distinguishable; bosons $+$; fermions $-$, for spatial wavefunctions with spins in the same state). Notation: $\langle x^n\rangle_a = \int x^n|\phi_a|^2dx$, $\langle x\rangle_{ab} = \int x\,\overline{\phi_a}\phi_b\,dx$.

▶ **Try it first.** Compute $\langle x_1^2\rangle$ and $\langle x_1x_2\rangle$ for $\psi_\pm$. Which cross terms vanish by orthogonality and which survive?

**Derivation.**
1. *Distinguishable:* $\langle(x_1 - x_2)^2\rangle_D = \langle x^2\rangle_a + \langle x^2\rangle_b - 2\langle x\rangle_a\langle x\rangle_b$.
2. *$\langle x_1^2\rangle_\pm$:* $\frac12\big[\langle x^2\rangle_a + \langle x^2\rangle_b\pm2\operatorname{Re}\big(\int\overline{\phi_a}x_1^2\phi_b\,dx_1\int\overline{\phi_b}\phi_a\,dx_2\big)\big]$. The second integral is $\langle\phi_b|\phi_a\rangle = 0$. So $\langle x_1^2\rangle_\pm = \frac12(\langle x^2\rangle_a + \langle x^2\rangle_b)$, and the same for $\langle x_2^2\rangle_\pm$.
3. *$\langle x_1x_2\rangle_\pm$:* the direct terms give $\frac12(\langle x\rangle_a\langle x\rangle_b + \langle x\rangle_b\langle x\rangle_a)$; the cross terms give $\pm\frac12(\langle x\rangle_{ab}\langle x\rangle_{ba} + \langle x\rangle_{ba}\langle x\rangle_{ab}) = \pm|\langle x\rangle_{ab}|^2$. So $\langle x_1x_2\rangle_\pm = \langle x\rangle_a\langle x\rangle_b\pm|\langle x\rangle_{ab}|^2$.
4. Combine:
$$
\langle(x_1 - x_2)^2\rangle_\pm = \langle(x_1 - x_2)^2\rangle_D\mp2|\langle x\rangle_{ab}|^2. \tag{3.17}
$$

■ **What just happened.** Bosons sit closer together and fermions farther apart than distinguishable particles in the same orbitals, by an amount set by the overlap integral $\langle x\rangle_{ab}$, which is nonzero only where the two orbitals overlap. Nothing pushes or pulls: there is no term in $\hat H$. The "exchange force" is a correlation built into the allowed states. Its energy consequences are real once there *is* an interaction: a repulsive $V(|x_1 - x_2|)$ costs fermions less, because they avoid each other. Example 4.8 puts in numbers for the infinite well. Spin matters: two electrons in the spin singlet have a symmetric spatial part and huddle like bosons (next section).

### 3.14 Spin and space together: helium, H₂, the periodic table

For electrons the *total* state, space ⊗ spin, is antisymmetric. With two electrons, $\psi(\mathbf r_1, \mathbf r_2)\chi(s_1, s_2)$ antisymmetric means either symmetric space with antisymmetric spin, or antisymmetric space with symmetric spin. QM-5 (3.41): the triplet is symmetric and the singlet antisymmetric. So
$$
\text{spin singlet}\ \Longleftrightarrow\ \text{symmetric }\psi(\mathbf r_1, \mathbf r_2),\qquad \text{spin triplet}\ \Longleftrightarrow\ \text{antisymmetric }\psi(\mathbf r_1, \mathbf r_2). \tag{3.18}
$$

**Helium.** $\hat H = \hat h_1 + \hat h_2 + V_{ee}$, with $\hat h$ hydrogenic with $Z = 2$ and $V_{ee} = e^2/(4\pi\varepsilon_0|\mathbf r_1 - \mathbf r_2|)$. Ignoring $V_{ee}$, the orbitals are QM-6's $\phi_{n\ell m}$ with energies $-Z^2(13.6\ \text{eV})/n^2$.
- *Ground state $1s^2$:* the only spatial state is $\phi_{1s}(\mathbf r_1)\phi_{1s}(\mathbf r_2)$, symmetric, so by (3.18) the spins must be in the singlet. There is no triplet $1s^2$.
- *Excited $1s\,2s$:* spatial states $\psi_\pm$ of (3.16) with $a = 1s$, $b = 2s$. Singlet with $\psi_+$ (**parahelium**), triplet with $\psi_-$ (**orthohelium**).

▶ **Try it first.** Compute $\langle\psi_\pm|V_{ee}|\psi_\pm\rangle$, writing it as a "direct" term plus or minus an "exchange" term.

**Derivation.**
1. Expand: $\langle\psi_\pm|V_{ee}|\psi_\pm\rangle = \frac12[\langle ab|V|ab\rangle + \langle ba|V|ba\rangle\pm\langle ab|V|ba\rangle\pm\langle ba|V|ab\rangle]$, where $|ab\rangle = \phi_a(\mathbf r_1)\phi_b(\mathbf r_2)$.
2. $V_{ee}$ is symmetric under $1\leftrightarrow2$, so $\langle ba|V|ba\rangle = \langle ab|V|ab\rangle =: J$ and $\langle ba|V|ab\rangle = \langle ab|V|ba\rangle =: K$ (real for real orbitals).
3. So $E_\pm = E_a + E_b + J\pm K$, with
$$
J = \iint\frac{e^2|\phi_a(\mathbf r_1)|^2|\phi_b(\mathbf r_2)|^2}{4\pi\varepsilon_0r_{12}},\qquad K = \iint\frac{e^2\,\phi_a(\mathbf r_1)\phi_b(\mathbf r_1)\,\phi_a(\mathbf r_2)\phi_b(\mathbf r_2)}{4\pi\varepsilon_0r_{12}}. \tag{3.19}
$$
4. $J$ is the classical repulsion of two charge clouds. $K$ is the Coulomb self-energy of the "overlap charge" $\phi_a\phi_b$, which is positive (the Coulomb kernel is a positive operator). So **orthohelium (triplet) lies below parahelium (singlet)** by $2K$.

**Numbers in the crude model** (unscreened hydrogenic orbitals, $Z = 2$; project Part 9 computes the integrals on a radial grid): $E_{1s} + E_{2s} = -2.5$ Ha $= -68.03$ eV, $J = \frac{34}{81}$ Ha $= 11.42$ eV, $K = \frac{32}{729}$ Ha $= 1.194$ eV. Orthohelium $-57.80$ eV, parahelium $-55.41$ eV, splitting $2K = 2.39$ eV. Measured: the $2^3S$ and $2^1S$ levels lie $19.82$ and $20.62$ eV above the ground state, a splitting of $0.80$ eV. Right sign, factor 3 too big (screening shrinks the overlap; QM-8 does better). The point: the energy depends on the spin state although $\hat H$ contains no spin. Pauli ties spin to spatial symmetry and Coulomb does the rest. This is the origin of Hund's first rule and of ferromagnetism: with $\hat{\mathbf S}_1\cdot\hat{\mathbf S}_2 = \frac{\hbar^2}{2}(s(s+1) - \frac32)$ (QM-5), $E_\pm$ is $E_a + E_b + J - \frac K2 - \frac{2K}{\hbar^2}\hat{\mathbf S}_1\cdot\hat{\mathbf S}_2$, Heisenberg's spin–spin coupling.

**The covalent bond in H₂.** Two protons $A, B$, two electrons. Heitler and London (1927) took $\psi_\pm\propto\phi_A(\mathbf r_1)\phi_B(\mathbf r_2)\pm\phi_B(\mathbf r_1)\phi_A(\mathbf r_2)$ with $\phi_{A,B}$ hydrogen $1s$ orbitals on each proton (now overlapping, so the normalization is $[2(1\pm\langle\phi_A|\phi_B\rangle^2)]^{-1/2}$). The energy has (3.19)'s direct-plus-or-minus-exchange structure, now including electron–proton attractions, and here the exchange term is *negative* at bonding distances: the symmetric $\psi_+$, which puts both electrons in the region between the protons where both attract them, is lower and binds; $\psi_-$ vanishes on the midplane and does not. By (3.18) $\psi_+$ is the spin singlet: **a covalent bond is a pair of electrons with opposite spins**, the chemists' Lewis pair. (Measured: $4.75$ eV deep at $0.74$ Å; QM-8 does the variational calculation.)

**The periodic table.** QM-6's orbitals $(n, \ell, m)$ each take two electrons ($m_s = \pm\frac12$), so subshells hold $2(2\ell + 1) = 2, 6, 10$; screening sets the filling order $1s, 2s, 2p, 3s, 3p, 4s, 3d, \ldots$ and the row lengths $2, 8, 8, 18, 18, 32$. Without antisymmetry every electron would sit in $1s$ and chemistry would not exist; Dyson and Lenard (1967) proved that bulk matter is stable because of the exclusion principle.

**Thermal consequence (SM-5).** For noninteracting identical particles, a many-body state is specified by occupation numbers $n_k$ of one-particle states; antisymmetry restricts fermions to $n_k\in\{0, 1\}$, bosons to $n_k\in\{0, 1, 2, \ldots\}$. Summing the grand-canonical Boltzmann factor over these allowed values gives $\langle n_k\rangle = 1/(e^{\beta(\varepsilon_k - \mu)} + 1)$ (Fermi–Dirac) and $1/(e^{\beta(\varepsilon_k - \mu)} - 1)$ (Bose–Einstein). White dwarfs, metals, Bose–Einstein condensates and blackbody radiation are this section at finite temperature; SM-5 derives them.

### 3.15 Decoherence: the classical world as a reduced description

A system $S$ in superposition $\alpha|0\rangle + \beta|1\rangle$ interacts with an environment $E$ that starts in $|e\rangle$. If the interaction records which state $S$ is in, $|0\rangle|e\rangle\mapsto|0\rangle|e_0\rangle$, $|1\rangle|e\rangle\mapsto|1\rangle|e_1\rangle$ (a CNOT-like unitary when $E$ is a qubit), linearity gives
$$
(\alpha|0\rangle + \beta|1\rangle)|e\rangle\ \longmapsto\ \alpha|0\rangle|e_0\rangle + \beta|1\rangle|e_1\rangle. \tag{3.20}
$$

▶ **Try it first.** Trace out $E$ with (3.9). What multiplies the off-diagonal element?

**Derivation (which-path).** QM-1's two-slit experiment: the particle passes slit $L$ or $R$, $|0\rangle = |L\rangle$, $|1\rangle = |R\rangle$, $\alpha = \beta = \frac{1}{\sqrt2}$, and $|e_0\rangle, |e_1\rangle$ normalized states of a which-path detector.
1. $\hat\rho_S = \operatorname{Tr}_E|\Psi\rangle\langle\Psi| = |\alpha|^2|L\rangle\langle L| + |\beta|^2|R\rangle\langle R| + \alpha\bar\beta\,\langle e_1|e_0\rangle\,|L\rangle\langle R| + \bar\alpha\beta\,\langle e_0|e_1\rangle\,|R\rangle\langle L|$. (Use $\operatorname{Tr}_E(|x\rangle\langle y|\otimes|e_i\rangle\langle e_j|) = \langle e_j|e_i\rangle|x\rangle\langle y|$.)
2. At screen position $x$, slit $L$ contributes amplitude $\psi_L(x) = \langle x|L\rangle$ and $R$ contributes $\psi_R(x)$. The detection probability is $P(x) = \langle x|\hat\rho_S|x\rangle$ (on the system, whatever the environment does): $P(x) = \frac12\big(|\psi_L|^2 + |\psi_R|^2 + 2\operatorname{Re}[\psi_L\overline{\psi_R}\langle e_1|e_0\rangle]\big)$.
3. Take equal envelopes $|\psi_L| = |\psi_R| = A$ and relative phase $\varphi(x)$, $\psi_L\overline{\psi_R} = A^2e^{i\varphi}$, and write $\langle e_1|e_0\rangle = |\gamma|e^{i\delta}$: $P = A^2[1 + |\gamma|\cos(\varphi + \delta)]$.
4. The **visibility** $(P_{\max} - P_{\min})/(P_{\max} + P_{\min})$ is
$$
\mathcal V = |\langle e_0|e_1\rangle|. \tag{3.21}
$$

■ **What just happened.** If the environment records the path perfectly ($\langle e_0|e_1\rangle = 0$), the fringes vanish and the screen shows the sum of two one-slit patterns, as if the particle went through one slit or the other at random. If the record is imperfect, the fringes survive with reduced contrast. Nobody needs to *read* the detector: the entanglement alone kills interference at $S$, and no operation on $S$ alone can restore it. A macroscopic object is recorded constantly by air molecules and photons, each multiplying the coherence by some $\langle e_0^{(k)}|e_1^{(k)}\rangle$ of modulus below 1, so coherences decay exponentially in the number of records (project Part 10). The reduced state becomes diagonal in the basis the environment monitors (position, for objects): an improper mixture (§3.9) indistinguishable from a classical distribution. Decoherence explains why big things don't interfere; it does not say why one outcome occurs, which is where §3.11's interpretations part. QM-9 and QC-7 make it quantitative ($T_2$, master equations).

## 4. Worked examples

### Example 4.1: $\sigma_z\otimes\hat1$, $\hat1\otimes\sigma_x$, and their commutator

By (3.1), with $\sigma_z = \operatorname{diag}(1, -1)$ and $\sigma_x = \begin{pmatrix}0&1\\1&0\end{pmatrix}$:
$$
\sigma_z\otimes\hat1 = \begin{pmatrix}1\cdot\hat1 & 0\\ 0 & -1\cdot\hat1\end{pmatrix} = \begin{pmatrix}1&0&0&0\\0&1&0&0\\0&0&-1&0\\0&0&0&-1\end{pmatrix},\qquad \hat1\otimes\sigma_x = \begin{pmatrix}\sigma_x&0\\0&\sigma_x\end{pmatrix} = \begin{pmatrix}0&1&0&0\\1&0&0&0\\0&0&0&1\\0&0&1&0\end{pmatrix}.
$$
Left-multiplying by $\sigma_z\otimes\hat1$ negates rows 3–4 of $\hat1\otimes\sigma_x$; right-multiplying negates columns 3–4. Both give
$$
\sigma_z\otimes\sigma_x = \begin{pmatrix}0&1&0&0\\1&0&0&0\\0&0&0&-1\\0&0&-1&0\end{pmatrix},
$$
so the commutator is zero (Theorem 3.1). Check on a basis vector: $\hat1\otimes\sigma_x$ flips qubit 1, $|10\rangle = \mathbf e_2\mapsto|11\rangle = \mathbf e_3$ (column 2 is $\mathbf e_3$). The trap: $\sigma_x\otimes\hat1$ is a different matrix (it swaps $\mathbf e_0\leftrightarrow\mathbf e_2$), and same-qubit operators need not commute: $[\sigma_z\otimes\hat1, \sigma_x\otimes\hat1] = 2i\sigma_y\otimes\hat1$. The project asserts all of these.

### Example 4.2: Schmidt decomposition of $(|00\rangle + |01\rangle + |10\rangle)/\sqrt3$

1. Coefficient matrix: $C = \frac{1}{\sqrt3}\begin{pmatrix}1&1\\1&0\end{pmatrix}$, $\det C = -\frac13\neq0$, so entangled.
2. $C$ is real symmetric, so its SVD comes from its eigendecomposition. $M = \begin{pmatrix}1&1\\1&0\end{pmatrix}$ has $\mu^2 - \mu - 1 = 0$, $\mu = \varphi = \frac{1+\sqrt5}{2} = 1.618$ and $\mu = -1/\varphi = -0.618$, with eigenvectors $\mathbf v_1 = (\varphi, 1)/\sqrt{\varphi^2 + 1} = (0.8507, 0.5257)$ and $\mathbf v_2 = (1, -\varphi)/\sqrt{1 + \varphi^2} = (0.5257, -0.8507)$.
3. $C = \frac{\varphi}{\sqrt3}\mathbf v_1\mathbf v_1^\top - \frac{1}{\varphi\sqrt3}\mathbf v_2\mathbf v_2^\top$. Singular values must be positive, so absorb the sign into the $B$ vector: $\lambda_1 = \varphi/\sqrt3 = 0.9342$, $|u_1\rangle = |v_1\rangle = 0.8507|0\rangle + 0.5257|1\rangle$; $\lambda_2 = 1/(\varphi\sqrt3) = 0.3568$, $|u_2\rangle = 0.5257|0\rangle - 0.8507|1\rangle$, $|v_2\rangle = -|u_2\rangle$.
$$
|\psi\rangle = 0.9342\,|u_1\rangle|u_1\rangle + 0.3568\,|u_2\rangle(-|u_2\rangle).
$$
4. Check: $\lambda_1^2 + \lambda_2^2 = (\varphi^2 + \varphi^{-2})/3 = 3/3 = 1$, since $\varphi^2 + \varphi^{-2} = (\varphi + 1) + (2 - \varphi) = 3$. The squares are $\lambda_{1,2}^2 = (3\pm\sqrt5)/6 = 0.8727, 0.1273$.
5. Entanglement entropy (3.10): $S = -0.8727\ln0.8727 - 0.1273\ln0.1273 = 0.3813$ nats $= 0.5500$ bits, a bit over half a Bell pair's worth.
6. Reduced state: $\hat\rho_A = CC^\dagger = \frac13\begin{pmatrix}2&1\\1&1\end{pmatrix}$, trace 1, eigenvalues $(3\pm\sqrt5)/6$ as promised by Theorem 3.6. Bloch vector $\mathbf r = (\frac23, 0, \frac13)$, $|\mathbf r| = \sqrt5/3 = 0.745$, and $\frac12(1 + |\mathbf r|) = 0.8727$. ✓

### Example 4.3: Reduced states of $|\Phi^+\rangle$ and of $\cos\alpha|00\rangle + \sin\alpha|11\rangle$

$|\Phi^+\rangle$: $C = \hat1/\sqrt2$, $\hat\rho_A = \hat1/2$, purity $\frac12$, $S = \ln2 = 0.6931$. By (3.9) directly: the nonzero entries of $\rho$ are $\rho_{00,00} = \rho_{00,11} = \rho_{11,00} = \rho_{11,11} = \frac12$, so $(\rho_A)_{00} = \rho_{00,00} + \rho_{01,01} = \frac12$, $(\rho_A)_{01} = \rho_{00,10} + \rho_{01,11} = 0$, $(\rho_A)_{11} = \frac12$. The coherence $\rho_{00,11}$ that carries the entanglement pairs $j = 0$ with $l = 1$, and the trace keeps only $j = l$.

$|\psi_\alpha\rangle = \cos\alpha|00\rangle + \sin\alpha|11\rangle$ ($0\le\alpha\le\frac\pi2$): already in Schmidt form, $\lambda = (\cos\alpha, \sin\alpha)$. $\hat\rho_A = \operatorname{diag}(\cos^2\alpha, \sin^2\alpha)$, Bloch vector $(0, 0, \cos2\alpha)$, and
$$
S(\alpha) = -\cos^2\alpha\ln\cos^2\alpha - \sin^2\alpha\ln\sin^2\alpha.
$$
$S(0) = 0$ (product $|00\rangle$). $S(\pi/6)$: $\cos^2 = \frac34$, $\sin^2 = \frac14$, $S = -\frac34\ln\frac34 - \frac14\ln\frac14 = 0.2158 + 0.3466 = 0.5623$ nats. $S(\pi/4) = \ln2$. With $p = \cos^2\alpha$, $dS/dp = \ln\frac{1-p}{p}$, zero only at $p = \frac12$, and $S$ is concave in $p$, so the maximum is at $\alpha = \pi/4$. The project plots $S(\alpha)$ and asserts the maximum location.

### Example 4.4: No-signaling, checked in two bases

Take $|\psi\rangle = \cos\alpha|00\rangle + \sin\alpha|11\rangle$ with $\alpha = \pi/6$. Before anything, $\hat\rho_B = \operatorname{diag}(\frac34, \frac14)$.

*Alice measures $z$.* Outcome 0 with $p = \cos^2\alpha = \frac34$, Bob left in $|0\rangle$; outcome 1 with $p = \frac14$, Bob in $|1\rangle$. Bob's average: $\frac34|0\rangle\langle0| + \frac14|1\rangle\langle1| = \operatorname{diag}(\frac34, \frac14)$.

*Alice measures $x$.* By (3.6) with $|a\rangle = |{\pm}\rangle$, $\bar{\boldsymbol\alpha} = (1, \pm1)/\sqrt2$: $C^\top\bar{\boldsymbol\alpha} = \frac{1}{\sqrt2}(\cos\alpha, \pm\sin\alpha)$, so $p(\pm) = \frac12$ and Bob is left in $|\beta_\pm\rangle = \cos\alpha|0\rangle\pm\sin\alpha|1\rangle$. Each has coherences $\pm\cos\alpha\sin\alpha = \pm0.4330$; the average is
$$
\tfrac12|\beta_+\rangle\langle\beta_+| + \tfrac12|\beta_-\rangle\langle\beta_-| = \begin{pmatrix}\cos^2\alpha & 0\\0 & \sin^2\alpha\end{pmatrix} = \operatorname{diag}(\tfrac34, \tfrac14).
$$
The coherences cancel. The two ensembles $\{\frac34: |0\rangle, \frac14: |1\rangle\}$ and $\{\frac12: |\beta_+\rangle, \frac12: |\beta_-\rangle\}$ are different (the second has Bloch vectors $(\pm\sin2\alpha, 0, \cos2\alpha) = (\pm0.866, 0, 0.5)$, the first $(0,0,\pm1)$) with the same center of mass $(0, 0, 0.5)$. Bob's probability of "0" is $\frac34$ either way. The project repeats this for random states and random bases.

### Example 4.5: Two mixtures, one $\hat\rho = \frac12\hat1$

$\frac12|0\rangle\langle0| + \frac12|1\rangle\langle1| = \frac12\begin{pmatrix}1&0\\0&0\end{pmatrix} + \frac12\begin{pmatrix}0&0\\0&1\end{pmatrix} = \frac12\hat1$.
$\frac12|{+}\rangle\langle{+}| + \frac12|{-}\rangle\langle{-}| = \frac14\begin{pmatrix}1&1\\1&1\end{pmatrix} + \frac14\begin{pmatrix}1&-1\\-1&1\end{pmatrix} = \frac12\hat1$.
Any set of Bloch vectors with center of mass zero works too (three at $120°$ on the equator: Exercise B4). For every $\hat Q$, $\langle\hat Q\rangle = \frac12\operatorname{Tr}\hat Q$ in each. By §3.8 the two above are what Bob holds after Alice measures half of $|\Phi^+\rangle$ in $z$ or in $x$ without telling him.

### Example 4.6: The singlet correlation and the CHSH angles

*Explicit computation.* Let $\mathbf a = \hat{\mathbf z}$ and $\mathbf b = (\sin\theta, 0, \cos\theta)$, so $\mathbf b\cdot\boldsymbol\sigma = \begin{pmatrix}\cos\theta & \sin\theta\\ \sin\theta & -\cos\theta\end{pmatrix}$. Apply $\sigma_z\otimes(\mathbf b\cdot\boldsymbol\sigma)$ to $|\Psi^-\rangle = (|01\rangle - |10\rangle)/\sqrt2$:
- $|01\rangle\mapsto(+1)|0\rangle\otimes(\sin\theta|0\rangle - \cos\theta|1\rangle)$,
- $|10\rangle\mapsto(-1)|1\rangle\otimes(\cos\theta|0\rangle + \sin\theta|1\rangle)$.

So the image is $\frac{1}{\sqrt2}[\sin\theta|00\rangle - \cos\theta|01\rangle + \cos\theta|10\rangle + \sin\theta|11\rangle]$. Its inner product with $\frac{1}{\sqrt2}(\langle01| - \langle10|)$ is $\frac12(-\cos\theta - \cos\theta) = -\cos\theta$. This is (3.14) with $\mathbf a\cdot\mathbf b = \cos\theta$. Outcome probabilities: the same-sign probability is $\frac12(1 + E) = \sin^2\frac\theta2$ and the opposite-sign probability $\cos^2\frac\theta2$.

*CHSH angles for spins* (all in the $xz$-plane, angle from $z$): $a = 0°$, $a' = 90°$, $b = 45°$, $b' = 135°$. Then $E(a,b) = -\cos45° = -0.7071$, $E(a,b') = -\cos135° = +0.7071$, $E(a',b) = -\cos45° = -0.7071$, $E(a',b') = -\cos45° = -0.7071$, and
$$
S_{\rm CHSH} = -0.7071 - 0.7071 - 0.7071 - 0.7071 = -2\sqrt2 = -2.8284.
$$
$|S| = 2\sqrt2 > 2$. (The sign is the singlet's anticorrelation; flip one of Bob's outcome labels and it becomes $+2\sqrt2$.)

*For photons*, $E = \cos2(\chi_a - \chi_b)$ with polarizer angles $\chi_a = 0°$, $\chi_{a'} = 45°$, $\chi_b = 22.5°$, $\chi_{b'} = 67.5°$: $E(a,b) = \cos45° = 0.7071$, $E(a,b') = \cos135° = -0.7071$, $E(a',b) = \cos45° = 0.7071$, $E(a',b') = \cos(-45°) = 0.7071$, so $S = 0.7071 + 0.7071 + 0.7071 + 0.7071 = 2\sqrt2$. Same geometry, half the angles. The project simulates $10^5$ singlet pairs with the spin angles, drawing outcomes by the Born rule, and estimates $S$ with its standard error.

### Example 4.7: Mermin's device

Two detectors, each with a switch (1, 2, 3) and a light that flashes red or green, and a source between them; switch settings are chosen at random and independently for each run. The data: (a) when the switches agree, the lights always flash the same color; (b) over all runs, the lights agree exactly half the time.

*Instruction sets.* With no wires between detectors, (a) suggests each particle carries the same instructions "color for 1, 2, 3". Take RRG. The nine equally likely setting pairs: the three with equal settings agree. Of the six unequal pairs, (1,2) and (2,1) agree (R, R), and (1,3), (3,1), (2,3), (3,2) disagree. Total agreement $5/9$. RRR or GGG agree $9/9$. Every set gives at least $5/9 > 1/2$. So (b) is impossible for instruction sets, whatever their distribution.

*Quantum version.* Use the singlet, axes for settings 1, 2, 3 in a plane at $0°$, $120°$, $240°$, and let Bob's detector flash red for spin down (so "same color" means opposite spins). Equal settings: opposite spins with certainty, so (a) holds. Unequal settings, $\theta = 120°$: opposite-spin probability $\cos^2\frac\theta2 = \cos^260° = \frac14$. Total: $\frac39\cdot1 + \frac69\cdot\frac14 = \frac13 + \frac16 = \frac12$. That is (b). The project checks the $5/9$ bound over all eight sets and the quantum $\frac12$ by simulation.

### Example 4.8: Exchange forces in the infinite well, $n = 1, 2$

Well of width $a$, $\phi_n = \sqrt{2/a}\sin(n\pi x/a)$. Needed integrals (Griffiths Problem 2.4 and Example 5.1): $\langle x\rangle_n = \frac a2$, $\langle x^2\rangle_n = a^2\big(\frac13 - \frac{1}{2n^2\pi^2}\big)$, and
$$
\langle x\rangle_{12} = \frac2a\int_0^a x\sin\frac{\pi x}{a}\sin\frac{2\pi x}{a}dx = -\frac{16a}{9\pi^2} = -0.18013\,a.
$$
(Use $\sin u\sin2u = \frac12(\cos u - \cos3u)$ and $\int_0^\pi u\cos(ku)\,du = ((-1)^k - 1)/k^2$: $\frac{2a}{\pi^2}\cdot\frac12\big(-2 + \frac29\big) = -\frac{16a}{9\pi^2}$.)

Numbers: $\langle x^2\rangle_1 = 0.28267a^2$, $\langle x^2\rangle_2 = 0.32067a^2$, $2\langle x\rangle_1\langle x\rangle_2 = 0.5a^2$, $2|\langle x\rangle_{12}|^2 = 0.06489a^2$.
$$
\langle(x_1 - x_2)^2\rangle_D = 0.10334a^2,\qquad \langle(x_1 - x_2)^2\rangle_+ = 0.03845a^2,\qquad \langle(x_1 - x_2)^2\rangle_- = 0.16823a^2.
$$
RMS separations $0.321a$, $0.196a$, $0.410a$. The exchange correction is more than half the distinguishable value, in both directions: bosons are found 39% closer, fermions 28% farther, with no force. The project computes all three by 2D quadrature on a grid and plots the three densities $|\psi(x_1, x_2)|^2$: the fermion density vanishes along the diagonal $x_1 = x_2$, the boson density piles up on it.

## 5. Why it matters and where it connects

- **QC-1 to QC-3.** Qubits are §3.1's spins, registers are tensor products, entangling gates are Kronecker matrices with the qubit-0-leftmost convention. QC-3's CHSH game is Theorem 3.9 ($\cos^2\frac\pi8\approx0.854$ winning probability is $\frac12 + \frac{2\sqrt2}{8}$) and QC-4's teleportation uses the Bell basis and §3.3's conditional states.
- **QC-7.** Everything in §3.4–3.9 is QC-7's starting point; §3.15 is the simplest **quantum channel** (dephasing), and $\operatorname{Tr}_E[\hat U(\hat\rho\otimes|e\rangle\langle e|)\hat U^\dagger]$ is the Stinespring form of a Kraus map. Entropies switch to bits.
- **SM-5.** §3.14's occupation-number argument gives Fermi–Dirac and Bose–Einstein. The von Neumann entropy is SM-3's Gibbs entropy $-\sum p\ln p$ applied to $\hat\rho$'s eigenvalues; a canonical ensemble is $\hat\rho = e^{-\beta\hat H}/Z$.
- **LA-12, LA-9, LA-6.** Schmidt = SVD; partial trace = einsum; Slater = determinant. Linear algebra you have already done is doing the physics.
- **QM-8.** Helium's $J$ and $K$ reappear as first-order perturbation theory, then the variational method fixes the screening; H₂ is the textbook variational molecule.
- **CM-8.** (3.12) is CM-8's Liouville equation quantized; the Bloch ball plays the role of the space of phase-space distributions.

## 6. Exercises

See `exercises.md`. A (6): Kronecker ordering, product-or-entangled, a Bloch vector, the Bell basis, Pauli counting, mixture vs superposition. B (9): local operators and partial traces, a Schmidt decomposition with entropy, no-signaling for unitaries, ensemble ambiguity, $|\Phi^+\rangle$'s correlations and CHSH, Werner states, oscillator exchange forces, helium's spin–space pairing, decoherence by many environment qubits. C (3): Tsirelson's operator bound, purification and steering, local strategies by enumeration (code). B2, B3, B5, B7, B9 are essential.

## 7. Coding project

`project.py` (a few seconds). Parts 1–3: Kronecker toolkit (Example 4.1), Bell states (rank 2, $\hat\rho_{A,B} = \frac12\hat1$, $S = \ln2$), Example 4.2, `einsum` partial trace against $CC^\dagger$, $S(\alpha)$ plotted. Part 4: random states in the Bloch ball, purity $\frac12(1 + |\mathbf r|^2)$. Part 5: the von Neumann equation by `solve_ivp`. Part 6: no-signaling for random states and measurements; Examples 4.4–4.5. Part 7: CHSH with $10^5$ Born-rule pairs, the 16-assignment enumeration and a local model with $|S|\le2$, Mermin, GHZ. Part 8: two particles in a box, (3.17) by 2D quadrature, densities plotted, a Slater determinant. Part 9: helium's $J$, $K$. Part 10: decoherence and visibility $|\langle e_0|e_1\rangle|$.

**Build it yourself first:** `ptrace(rho, dA, dB, keep)` by reshape and `einsum`, and the CHSH Monte Carlo (draw Alice's outcome with probability $\frac12$, then Bob's conditional on it from (3.6)), before reading the file's versions. Then the "Your turn" items: a random $3\times4$ state's Schmidt decomposition and equal reduced spectra, a Werner state's entanglement against its mixing parameter, and a three-fermion density.

## 8. Where you'll struggle

- **"Measuring A collapses B instantly."** True as a statement about the description conditional on A's result; false as a statement about anything B can detect. Most confusion about entanglement and relativity comes from skipping Theorem 3.8. The fix: do Example 4.4 by hand once, then state in words which object changed (B's conditional state) and which didn't (B's unconditional $\hat\rho_B$).
- **The partial trace is an einsum and you'll make it hard.** People write loops with four indices and get the pairing wrong. Reshape to `(dA, dB, dA, dB)`, then `'ijkj->ik'` keeps $A$ and `'ijil->jl'` keeps $B$. Test on $\hat X\otimes\hat Y$, whose partial traces are $\hat X\operatorname{Tr}\hat Y$ and $\hat Y\operatorname{Tr}\hat X$; a pairing bug shows up immediately.
- **The tensor ordering convention will transpose your matrices silently.** `np.kron(A, B)` puts $A$ on the left factor. If you reshape a state as `(dB, dA)` by mistake, $\hat\rho_A$ and $\hat\rho_B$ swap, and for asymmetric states you get the wrong answer with no error. QC-2 inherits this: a CNOT with control on qubit 1 is a different matrix from control on qubit 0. Always test with a product $|0\rangle|1\rangle$ whose reduced states you know.
- **Mixed state vs superposition is the module's most common conceptual error.** $(|0\rangle + |1\rangle)/\sqrt2$ is not "50% 0 and 50% 1." Ask which measurement distinguishes them ($\sigma_x$) and look at the off-diagonal entries. If you catch yourself saying "the particle is in both states," write the density matrix.
- **The Bell counting proof is simple and you still won't believe the conclusion.** Everyone tries to invent a clever local model that beats 2. Write your own in the project's Part 7 and watch it fail; Theorem 3.9's step 2 is why. Then read Mermin.
- **Spin-statistics is stated, not derived.** Nonrelativistic QM cannot prove it; don't hunt for a proof in Griffiths. Weinberg's *Quantum Theory of Fields* vol. 1 ch. 5 has one.
- **Exchange forces aren't forces.** There is no term in $\hat H$; (3.17) is pure correlation from the symmetry of the state. The energy splittings in §3.14 come from the ordinary Coulomb term evaluated in differently-shaped states. The name misleads everyone once.
- **Getting the helium pairing backwards.** Triplet = symmetric spin = antisymmetric space = electrons apart = lower Coulomb energy. Write the chain each time.
- **Interpretation questions are real.** Don't dismiss them and don't think one is proven. The experiments rule out (i)–(iii) together and nothing else; keep the three in front of you when reading popular accounts.

## 9. Check yourself

1. What is the dimension of the state space of 3 qubits and one spin-1? Why not the sum?
2. Why do $\hat A\otimes\hat1$ and $\hat1\otimes\hat B$ commute, and what does that mean physically?
3. State the rank criterion for a product state. Is $|00\rangle - |01\rangle + |10\rangle - |11\rangle$ entangled?
4. Derive $\hat\rho_A = CC^\dagger$ from $\langle\hat A\otimes\hat1\rangle = \operatorname{Tr}(\hat\rho_A\hat A)$.
5. Prove $\operatorname{Tr}\hat\rho^2 = 1$ iff pure.
6. Give the eigenvalues and purity of $\frac12(\hat1 + \mathbf r\cdot\boldsymbol\sigma)$ in terms of $|\mathbf r|$.
7. State the von Neumann equation and explain its sign relative to Heisenberg's.
8. Prove no-signaling for a projective measurement in three lines.
9. Write two different ensembles with $\hat\rho = \operatorname{diag}(\frac34, \frac14)$.
10. Derive $|S_{\rm CHSH}|\le2$ for local hidden variables, and name the three assumptions.
11. Compute $E(\mathbf a,\mathbf b)$ for the singlet and the angles that give $2\sqrt2$.
12. Why does antisymmetry imply Pauli exclusion, and why does the helium triplet lie lower?
13. Why does a which-path record with $\langle e_0|e_1\rangle = 0.6$ reduce the fringe visibility to $0.6$?

## 10. Pair with

- **Griffiths & Schroeter (3rd ed.), ch. 5 (identical particles) and ch. 12 (afterword).** The standard treatment of §3.12–3.14 with Example 5.1's exchange integrals; ch. 12 is a clean EPR/Bell/no-clone summary.
- **Townsend, *A Modern Approach to Quantum Mechanics*, ch. 5.** Two-spin systems, the singlet and Bell's inequality done carefully with spins first; better than Griffiths on EPR.
- **Nielsen & Chuang, §2.3–2.6.** The density operator, partial trace, Schmidt decomposition, purification and HJW theorem with full proofs. Read now; QC-1 assumes it.
- **Sakurai & Napolitano, §3.4.** The density operator with ensemble averages and the connection to quantum statistical mechanics.
- **Mermin, "Is the moon there when nobody looks? Reality and the quantum theory," *Physics Today* 38(4), 38 (1985).** Required. The clearest account of Bell's theorem ever written; Example 4.7 is its device.
- **Preskill, Ph219 lecture notes, ch. 2 (free online).** Density operators, Bloch ball, Schmidt, HJW, and the beginnings of open systems, at QC level.
