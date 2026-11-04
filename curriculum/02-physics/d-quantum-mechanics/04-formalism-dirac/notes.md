# QM-4: Formalism: Hilbert Space, Dirac Notation, Hermitian Operators, Eigenbases, the Generalized Statistical Interpretation, the Uncertainty Principle Derived, Time Evolution, Pictures

**Time:** ~16 hours · **Prerequisites:** QM-3, LA-8, LA-12 · **Unlocks:** QM-5, QM-7, QC-1

> **Conventions in this module.** $\hbar = h/2\pi$ is explicit; `project.py` uses $\hbar = m = \omega = 1$. Kets $|\psi\rangle$ are vectors; bras $\langle\phi|$ are linear functionals, the Riesz partner of $|\phi\rangle$ (LA-12 Theorem 3.9). The inner product $\langle\phi|\psi\rangle$ is **conjugate-linear in the first slot** and linear in the second (LA-8's physics convention; Axler's book is linear in the first slot, so its formulas carry conjugates in different places). $\hat A^\dagger$ is the adjoint, defined by $\langle\phi|\hat A\psi\rangle = \langle\hat A^\dagger\phi|\psi\rangle$. $[\hat A, \hat B] = \hat A\hat B - \hat B\hat A$. Resolutions of the identity: $\sum_n|n\rangle\langle n| = \hat 1$ and $\int|x\rangle\langle x|\,dx = \hat 1$ with $\langle x|x'\rangle = \delta(x - x')$. Momentum eigenfunctions are $\langle x|p\rangle = (2\pi\hbar)^{-1/2}e^{ipx/\hbar}$, matching QM-2's $\Phi(p) = (2\pi\hbar)^{-1/2}\int\Psi e^{-ipx/\hbar}dx$. Time dependence uses the physics sign $e^{-iEt/\hbar}$. The **Schrödinger picture** (states move, operators fixed) is the default; the propagator is $\hat U(t) = e^{-i\hat Ht/\hbar}$. Pauli matrices $\sigma_x, \sigma_y, \sigma_z$ are as in QM-1 (3.13), with $\boldsymbol\sigma = (\sigma_x, \sigma_y, \sigma_z)$.

## 0. Where this sits

You need QM-2 ($\hat p = -i\hbar\partial_x$, Ehrenfest (3.11), momentum space), QM-3 (box, ladder operators, Fourier trick), LA-8 (adjoints; spectral theorem, Theorem 3.12; Cauchy–Schwarz, Theorem 3.4), LA-12 §3.6 (Riesz, bras), F3's delta function, and CM-8's Poisson brackets (3.12)–(3.14) and Dirac's rule (3.19).

This module turns QM-2/3 into linear algebra. QM-5 is commutator algebra built on §3.7; QM-7 needs outer products and the density operator (§3.13); QC-1's qubit is §3.12's two-state system.

## 1. The problem that forced this idea

In July 1925 Werner Heisenberg sent off a paper that replaced electron orbits with arrays of numbers, one for each pair of states, and found that his rule for multiplying the arrays was not commutative. Max Born recognized them as matrices, and Born, Heisenberg and Pascual Jordan (the "three-man paper," 1926) built matrix mechanics with $\hat x\hat p - \hat p\hat x = i\hbar$ as its central law. Dirac reached the same law independently from the Poisson bracket (CM-8 §1(d)). Months later Schrödinger's wave mechanics appeared: a differential equation, no matrices anywhere, and the same hydrogen spectrum. In 1926 Schrödinger showed why: the integrals $\int\psi_m^*\hat Q\psi_n\,dx$ between his eigenfunctions are Heisenberg's matrices.

That was a computation, not a framework. John von Neumann (papers from 1927; *Mathematical Foundations of Quantum Mechanics*, 1932) identified the state space as a Hilbert space and observables as self-adjoint operators on it. Matrix mechanics is the theory written in the energy eigenbasis; wave mechanics is the same theory in the "position basis." Von Neumann saw that the position basis is a problem: the eigenfunctions of $\hat x$ are delta functions and those of $\hat p$ are plane waves, neither square-integrable. He called the delta function a mathematical fiction and built his theory without it.

Dirac's *Principles of Quantum Mechanics* (1930) went the other way, with abstract vectors and delta-normalized continuous bases; the bra-ket notation came in a 1939 paper and has been standard since the *Principles*' third edition (1947). It is designed so that the spectral theorem does the work: drop $\hat 1 = \sum_n|n\rangle\langle n|$ between any two symbols and every change of basis, probability and matrix element falls out. The mathematics that makes the continuous bases honest (rigged Hilbert spaces, Gelfand and collaborators, 1950s–60s) came a generation after physicists were already using them correctly.

The forcing problem for you is the same: wave mechanics works but hides its structure. QM-3's eigenfunction expansions, Fourier transforms and ladder operators were one operation in different bases. If you did LA-8 and LA-12, you have done the hard part: QM-2/3 was finite-dimensional linear algebra with an infinite index.

## 2. Intuition

Take $\mathbf v\in\mathbb C^3$ with components $v_i = \mathbf e_i^\dagger\mathbf v$, and let the index become continuous: $\psi(x)$ for $x\in\mathbb R$. **A wavefunction is a column vector whose rows are labeled by $x$.** The dot product $\mathbf u^\dagger\mathbf v = \sum_i\bar u_iv_i$ becomes $\int\overline{\phi(x)}\psi(x)\,dx$. An operator is a matrix with two continuous indices, $Q(x, x')$, and acting on $\psi$ means $\int Q(x, x')\psi(x')\,dx'$. Multiplication by $V(x)$ is a diagonal matrix; differentiation is a matrix concentrated next to the diagonal (on QM-2's grid, literally the central-difference matrix).

With that one picture, translate QM-2 and QM-3 line by line:

| Wave mechanics | Linear algebra |
|---|---|
| $\int\phi^*\psi\,dx$ | $\mathbf u^\dagger\mathbf v$ |
| $\hat Q$ Hermitian | $A^\dagger = A$ (symmetric with conjugates) |
| $\hat H\psi_n = E_n\psi_n$ | eigenvectors of a Hermitian matrix |
| eigenfunctions orthogonal, complete | spectral theorem: orthonormal eigenbasis |
| Fourier trick $c_n = \int\psi_n^*\psi\,dx$ | coordinates in an orthonormal basis, $c_n = \mathbf q_n^\dagger\mathbf v$ (LA-8 Theorem 3.6) |
| $\Phi(p)$, the Fourier transform | the same vector in another basis; the Fourier kernel is the change-of-basis matrix |
| "measure energy, get $E_n$ with probability $|c_n|^2$" | squared length of the projection onto the $n$-th eigenvector |

Two guesses you should be able to make before §3:

1. If $\hat Q$ has eigenvectors $|n\rangle$ and $|\psi\rangle = \sum c_n|n\rangle$, then $\langle\psi|\hat Q|\psi\rangle = \sum_n q_n|c_n|^2$. In matrix terms, a Hermitian quadratic form diagonalized.
2. Two observables are both sharp only in a common eigenvector. Commuting matrices share an eigenbasis (LA-8); for non-commuting ones the commutator should bound how sharp both can be. That bound is the uncertainty principle, and it will come from Cauchy–Schwarz.

Measurement in this picture: an observable is a set of orthogonal directions (its eigenvectors) with a number on each, and measuring snaps the state's unit arrow onto one direction with probability equal to the squared cosine of the angle. QM-1's spin measurements were this in $\mathbb C^2$.

## 3. Formal development

### 3.1 Hilbert space, and what is swept under the rug

The state space of one particle on a line is $L^2(\mathbb R)$: complex functions with $\int|\psi|^2dx < \infty$, with the inner product
$$
\langle\phi|\psi\rangle = \int_{-\infty}^\infty\overline{\phi(x)}\,\psi(x)\,dx, \qquad \langle\phi|\psi\rangle = \overline{\langle\psi|\phi\rangle}, \qquad \langle\psi|\psi\rangle\ge0 \text{ with equality iff } \psi = 0. \tag{3.1}
$$
It is a vector space (LA-3: a sum of square-integrable functions is square-integrable, by $|f + g|^2\le2|f|^2 + 2|g|^2$) and an inner-product space (LA-8). A **Hilbert space** is an inner-product space that is **complete**: every sequence $\psi_k$ with $\|\psi_k - \psi_l\|\to0$ converges to something in the space. Completeness is what lets an infinite expansion $\sum c_n\psi_n$ with $\sum|c_n|^2 < \infty$ define a vector. $L^2$ is complete (the Riesz–Fischer theorem); the continuous functions with the same inner product are not, because a sequence of smooth functions can converge in norm to a step.

The honest page: three things differ from LA-8's $\mathbb C^n$, and this course treats all three informally.

1. **Elements are equivalence classes.** Functions differing on a set of measure zero are the same vector, so $\psi$ "at a point" is not defined for a general $L^2$ element.
2. **Unbounded operators have domains.** $\hat x\psi = x\psi$ and $\hat p\psi = -i\hbar\psi'$ are not defined on all of $L^2$: $x\cdot(1 + x^2)^{-1/2}$ is not square-integrable even though $(1 + x^2)^{-1/2}$ is. Hermiticity proofs drop boundary terms (QM-2 Theorem 3.4), legitimate only on a suitable domain. "Hermitian" (symmetric) and "self-adjoint" (same domain as the adjoint) differ here, and the spectral theorem needs self-adjoint; we say the first and mean the second.
3. **Some eigenvectors are not in the space.** $\hat x$'s eigenfunction with eigenvalue $x_0$ is $\delta(x - x_0)$; $\hat p$'s is $e^{ip_0x/\hbar}$. Neither has finite norm. They live in a bigger space of distributions (F3), and the structure $\Phi\subset L^2\subset\Phi'$ (nice test functions, Hilbert space, distributions) is called a **rigged Hilbert space** or Gelfand triple. In it, the spectral theorem for $\hat x$ and $\hat p$ holds with delta-normalized eigenvectors.

We proceed as physicists do: use $|x\rangle$ and $|p\rangle$ freely, delta-normalized, and apply the results to normalizable states.

### 3.2 Kets, bras, and outer products

▶ **Try it first.** For $\mathbf v = (1, i)^\top$, is the functional $\mathbf u\mapsto\langle\mathbf v, \mathbf u\rangle$ the row $(1, i)$ or $(1, -i)$?

A **ket** $|\psi\rangle$ is an element of the Hilbert space $\mathcal H$. A **bra** $\langle\phi|$ is a continuous linear functional on $\mathcal H$: it eats a ket and returns a number, linearly. The Riesz representation theorem (LA-12 Theorem 3.9 in finite dimensions; it holds verbatim for continuous functionals on any Hilbert space) says every bra is of the form $|\psi\rangle\mapsto\langle\phi|\psi\rangle$ for exactly one ket $|\phi\rangle$. So bras and kets are in one-to-one correspondence, and the correspondence is conjugate-linear: the bra of $c|\phi\rangle$ is $\bar c\langle\phi|$. In coordinates a ket is a column and its bra is the conjugate-transposed row, which answers the prompt: $(1, -i)$.

The **outer product** $|\alpha\rangle\langle\beta|$ is the operator
$$
\big(|\alpha\rangle\langle\beta|\big)|\psi\rangle = |\alpha\rangle\,\langle\beta|\psi\rangle, \tag{3.2}
$$
a column times a row: a rank-one map whose range is the line through $|\alpha\rangle$ (LA-4). Its adjoint is $|\beta\rangle\langle\alpha|$. If $\langle n|n\rangle = 1$, then $\hat P_n = |n\rangle\langle n|$ is the orthogonal projection onto the line through $|n\rangle$: $\hat P_n^2 = |n\rangle\langle n|n\rangle\langle n| = \hat P_n$ and $\hat P_n^\dagger = \hat P_n$ (LA-5, LA-8).

### 3.3 Resolutions of the identity, and what a wavefunction is

▶ **Try it first.** Let $\{|n\rangle\}$ be an orthonormal basis. Write $|\psi\rangle = \sum_nc_n|n\rangle$, take the inner product with $\langle m|$, and solve for $c_m$. Then put the answer back into the expansion.

**Derivation (discrete).**

1. Orthonormal means $\langle m|n\rangle = \delta_{mn}$. Basis (with completeness of $\mathcal H$) means every $|\psi\rangle$ is $\sum_nc_n|n\rangle$ with $\sum|c_n|^2 < \infty$.
2. Apply $\langle m|$ and use linearity in the second slot: $\langle m|\psi\rangle = \sum_nc_n\langle m|n\rangle = c_m$.
3. Substitute: $|\psi\rangle = \sum_n|n\rangle\langle n|\psi\rangle = \big(\sum_n|n\rangle\langle n|\big)|\psi\rangle$ for every $|\psi\rangle$. So
$$
\sum_n|n\rangle\langle n| = \hat 1. \tag{3.3}
$$ ∎

**Derivation (continuous).** Postulate kets $|x\rangle$, one for each $x\in\mathbb R$, that are eigenvectors of position, $\hat x|x\rangle = x|x\rangle$, with the continuous version of orthonormality, $\langle x|x'\rangle = \delta(x - x')$, and suppose they span: every $|\psi\rangle = \int c(x')|x'\rangle\,dx'$.

1. Apply $\langle x|$: $\langle x|\psi\rangle = \int c(x')\langle x|x'\rangle\,dx' = \int c(x')\delta(x - x')\,dx' = c(x)$. The delta does what $\delta_{mn}$ did.
2. Substitute: $|\psi\rangle = \int|x'\rangle\langle x'|\psi\rangle\,dx'$ for all $|\psi\rangle$, so
$$
\langle x|x'\rangle = \delta(x - x'), \qquad \int|x\rangle\langle x|\,dx = \hat 1. \tag{3.4}
$$ ∎

The coefficients have names. In the energy basis of QM-3 they are $c_n = \langle n|\psi\rangle$. In the position basis they are the wavefunction:
$$
\boxed{\psi(x) = \langle x|\psi\rangle}, \qquad \psi_n = \langle n|\psi\rangle. \tag{3.5}
$$
**The wavefunction is the list of coefficients of the state in the position basis.** It is not the state; it is the state's coordinates in one particular basis, the way $(3, 2)$ is not an arrow but the arrow's coordinates. Inserting (3.4) into $\langle\phi|\psi\rangle$ recovers (3.1): $\langle\phi|\psi\rangle = \int\langle\phi|x\rangle\langle x|\psi\rangle\,dx = \int\overline{\phi(x)}\psi(x)\,dx$, since $\langle\phi|x\rangle = \overline{\langle x|\phi\rangle}$.

■ **What just happened.** "Insert a complete set of states" is the only move in Dirac notation. (3.3) and (3.4) are the spectral theorem's statement that an orthonormal eigenbasis exists, written so that it can be dropped between any two symbols.

### 3.4 The momentum basis, and the Fourier transform as a change of basis

In the position basis, $\hat p$ acts as QM-2 derived: $\langle x|\hat p|\psi\rangle = -i\hbar\,\partial_x\langle x|\psi\rangle$. We want $\langle x|p\rangle$, the position-basis components of the momentum eigenket $|p\rangle$.

▶ **Try it first.** Write $\hat p|p\rangle = p|p\rangle$ in components using the rule just stated, solve the ODE, and fix the constant by requiring $\langle p|p'\rangle = \delta(p - p')$. You need F3's $\int e^{ikx}dx = 2\pi\delta(k)$.

**Derivation.**

1. Take $\langle x|$ of $\hat p|p\rangle = p|p\rangle$: $-i\hbar\,\partial_x\langle x|p\rangle = p\langle x|p\rangle$. The solution is $\langle x|p\rangle = N_pe^{ipx/\hbar}$.
2. Insert (3.4): $\langle p|p'\rangle = \int\langle p|x\rangle\langle x|p'\rangle\,dx = \bar N_pN_{p'}\int e^{i(p' - p)x/\hbar}dx$.
3. F3's $\hat\delta = 1$ inverted gives $\int e^{ikx}dx = 2\pi\delta(k)$. With $k = (p' - p)/\hbar$ and $\delta(k) = \hbar\,\delta(p' - p)$: $\langle p|p'\rangle = |N_p|^2\,2\pi\hbar\,\delta(p - p')$.
4. Dirac orthonormality $\langle p|p'\rangle = \delta(p - p')$ requires $|N_p| = (2\pi\hbar)^{-1/2}$; choose it real and positive:
$$
\langle x|p\rangle = \frac{e^{ipx/\hbar}}{\sqrt{2\pi\hbar}}, \qquad \langle p|p'\rangle = \delta(p - p'), \qquad \int|p\rangle\langle p|\,dp = \hat 1. \tag{3.6}
$$
5. The momentum-space wavefunction is the coefficient in the momentum basis, and one insertion of (3.4) computes it:
$$
\Phi(p) = \langle p|\psi\rangle = \int\langle p|x\rangle\langle x|\psi\rangle\,dx = \frac{1}{\sqrt{2\pi\hbar}}\int e^{-ipx/\hbar}\psi(x)\,dx. \tag{3.7}
$$ ∎

■ **What just happened.** QM-2's "the Fourier transform of $\Psi$ is the momentum amplitude" is a one-line change of basis (LA-4). The "matrix" $\langle p|x\rangle$ is unitary, and Plancherel is unitarity. On a periodic $N$-point grid it is the DFT matrix $F_{kj} = e^{-ikx_j}/\sqrt N$; `project.py` Part 3 checks $F^\dagger F = I$, Parseval, and that the central-difference $\hat p$ is diagonal in that basis with eigenvalues $\hbar\sin(k\,\Delta x)/\Delta x$.

**Matrix elements.** Inserting (3.3) twice into $\hat Q|\psi\rangle$,
$$
\langle m|\hat Q|\psi\rangle = \sum_n\langle m|\hat Q|n\rangle\langle n|\psi\rangle, \qquad Q_{mn} = \langle m|\hat Q|n\rangle. \tag{3.8}
$$
The state becomes a column $(\psi_n)$ and the operator a matrix $(Q_{mn})$. Heisenberg's 1925 arrays were exactly these, in the energy eigenbasis. In the position basis $\langle x|\hat x|x'\rangle = x\,\delta(x - x')$ and $\langle x|\hat p|x'\rangle = -i\hbar\,\delta'(x - x')$: §2's picture.

### 3.5 Observables: Hermitian operators and their theorems

The adjoint is defined by
$$
\langle\phi|\hat A\psi\rangle = \langle\hat A^\dagger\phi|\psi\rangle \quad\text{for all }\phi, \psi, \tag{3.9}
$$
which in an orthonormal basis makes $(A^\dagger)_{mn} = \overline{A_{nm}}$ (LA-8 Theorem 3.8). Read (3.9) carefully: $\hat A$ acting to the right on the ket equals $\hat A^\dagger$ acting on the bra's ket. Equivalently $\langle\phi|\hat A|\psi\rangle = \overline{\langle\psi|\hat A^\dagger|\phi\rangle}$. Properties: $(\hat A\hat B)^\dagger = \hat B^\dagger\hat A^\dagger$, $(c\hat A)^\dagger = \bar c\hat A^\dagger$, $(|\alpha\rangle\langle\beta|)^\dagger = |\beta\rangle\langle\alpha|$. An **observable** is a Hermitian operator, $\hat Q^\dagger = \hat Q$.

**Theorem 3.1.** $\hat x$, $\hat p$, and $\hat H = \hat p^2/2m + V(\hat x)$ with real $V$ are Hermitian.

*Proof.* (QM-2 Theorem 3.4 in Dirac form.) $\langle\phi|\hat x|\psi\rangle = \int\bar\phi\,x\,\psi\,dx = \overline{\langle\psi|\hat x|\phi\rangle}$ since $x$ is real. For $\hat p$, integrate by parts once and drop the boundary term (on the domain of §3.1): $\int\bar\phi(-i\hbar\psi')dx = \int\overline{(-i\hbar\phi')}\psi\,dx$. Products: $(\hat p^2)^\dagger = \hat p^\dagger\hat p^\dagger = \hat p^2$; $V(\hat x)$ is multiplication by a real function. Sums of Hermitian operators are Hermitian. ∎

The next three are LA-8 Theorems 3.9, 3.10 and 3.12, restated with their quantum meaning.

**Theorem 3.2 (Real eigenvalues).** If $\hat Q^\dagger = \hat Q$ and $\hat Q|q\rangle = q|q\rangle$ with $|q\rangle\ne0$, then $q\in\mathbb R$. *Proof.* $q\langle q|q\rangle = \langle q|\hat Qq\rangle = \langle\hat Qq|q\rangle = \bar q\langle q|q\rangle$, and $\langle q|q\rangle > 0$. ∎ Measurement outcomes are real numbers.

**Theorem 3.3 (Orthogonality).** Eigenvectors of a Hermitian $\hat Q$ with distinct eigenvalues are orthogonal. *Proof.* $q'\langle q|q'\rangle = \langle q|\hat Qq'\rangle = \langle\hat Qq|q'\rangle = q\langle q|q'\rangle$ (using $q$ real), so $(q' - q)\langle q|q'\rangle = 0$. ∎ Distinct outcomes are perfectly distinguishable states.

**Theorem 3.4 (Completeness; a postulate in infinite dimensions).** For an observable $\hat Q$, the eigenvectors (discrete, delta-normalized continuous, or both) can be chosen orthonormal and complete:
$$
\hat Q = \sum_nq_n|n\rangle\langle n| + \int q\,|q\rangle\langle q|\,dq, \qquad \sum_n|n\rangle\langle n| + \int|q\rangle\langle q|\,dq = \hat 1. \tag{3.10}
$$
In $\mathbb C^N$ this is LA-8 Theorem 3.12, proved. For self-adjoint operators on a Hilbert space it is von Neumann's spectral theorem, whose proof needs measure theory (Reed and Simon vol. I, ch. VII). For QM-3's box and oscillator the discrete part alone is complete (Sturm–Liouville); for $\hat x$ and $\hat p$ only the continuous part exists. Physics postulates (3.10).

**Degeneracy.** If several orthonormal eigenvectors share one eigenvalue $q$, the eigenspace $E_q$ has dimension $>1$ and its basis is not unique. Any orthonormal basis of $E_q$ (Gram–Schmidt, LA-5) works, and the projector $\hat P_q = \sum_{i}|q, i\rangle\langle q, i|$ is the same for every choice.

**Theorem 3.5 (Determinate states are eigenstates).** For Hermitian $\hat Q$ and a normalized $|\psi\rangle$, $\sigma_Q = 0$ if and only if $\hat Q|\psi\rangle = \langle Q\rangle|\psi\rangle$.

▶ **Try it first.** Write $\sigma_Q^2 = \langle\psi|(\hat Q - \langle Q\rangle)^2|\psi\rangle$ and move one factor onto the bra.

*Proof.* Let $\hat D = \hat Q - \langle Q\rangle$, Hermitian since $\langle Q\rangle$ is real. Then $\sigma_Q^2 = \langle\psi|\hat D\hat D\psi\rangle = \langle\hat D\psi|\hat D\psi\rangle = \|\hat D\psi\|^2$, which is zero iff $\hat D|\psi\rangle = 0$. ∎

■ **What just happened.** "A state with a definite value of $Q$" and "an eigenvector of $\hat Q$" are the same thing, and $\sigma_Q = \|(\hat Q - \langle Q\rangle)\psi\|$ is the distance by which $\hat Q|\psi\rangle$ misses being parallel to $|\psi\rangle$.

### 3.6 The generalized statistical interpretation

With (3.10) in hand, QM-1's postulates can be stated in full.

**Postulate (Born rule, Lüders collapse).** Measuring the observable $\hat Q$ on the normalized state $|\psi\rangle$:

- (discrete) gives eigenvalue $q$ with probability
$$
P(q) = \langle\psi|\hat P_q|\psi\rangle = \sum_i|\langle q, i|\psi\rangle|^2, \qquad\text{nondegenerate: } P(q_n) = |\langle n|\psi\rangle|^2 = |c_n|^2; \tag{3.11}
$$
- (continuous) gives a value in $[q, q + dq]$ with probability $|\langle q|\psi\rangle|^2dq$;
- leaves the system, immediately after outcome $q$, in the normalized projection
$$
|\psi\rangle\ \longmapsto\ \frac{\hat P_q|\psi\rangle}{\|\hat P_q|\psi\rangle\|}. \tag{3.12}
$$

The degenerate form of (3.12) is due to Gerhart Lüders (1951); von Neumann's 1932 version was ambiguous there. Position and momentum are the special cases: $|\langle x|\psi\rangle|^2dx = |\psi(x)|^2dx$ is QM-2's Born rule and $|\langle p|\psi\rangle|^2dp = |\Phi(p)|^2dp$ is QM-2's momentum distribution. Nothing in QM-2 is lost; all of it is now one rule applied in different bases.

**Consistency checks.** (i) The probabilities sum to one: $\sum_n|c_n|^2 = \sum_n\langle\psi|n\rangle\langle n|\psi\rangle = \langle\psi|\hat 1|\psi\rangle = 1$, by (3.3). That is Parseval (LA-8 Theorem 3.6). (ii) The expectation value:

▶ **Try it first.** Show that the average of the outcomes, $\sum q_nP(q_n)$, equals $\langle\psi|\hat Q|\psi\rangle$. Insert (3.3) to the right of $\hat Q$.

1. $\langle\psi|\hat Q|\psi\rangle = \langle\psi|\hat Q\sum_n|n\rangle\langle n|\psi\rangle = \sum_n\langle\psi|\hat Q|n\rangle\,c_n$.
2. $\hat Q|n\rangle = q_n|n\rangle$, so $\langle\psi|\hat Q|n\rangle = q_n\langle\psi|n\rangle = q_n\bar c_n$.
3. Hence
$$
\langle Q\rangle = \sum_nq_n|c_n|^2 = \langle\psi|\hat Q|\psi\rangle. \tag{3.13}
$$ ∎

(iii) Repeatability: measuring $\hat Q$ again right after outcome $q$ gives $q$ with probability one, because the collapsed state lies in $E_q$ and $\hat P_q$ acts on it as the identity. `project.py` Part 7 simulates $10^4$ measurements of a three-level observable, compares the frequencies with (3.11), and checks repeatability.

■ **What just happened.** Probabilities are squared lengths of projections, the mean is the quadratic form, and collapse is projection followed by renormalization. All three are LA-5/LA-8 geometry. That nature does this is physics, not geometry; QM-9's decoherence explains when (3.12) is a good effective description.

### 3.7 Commutators, compatibility, and functions of operators

The commutator $[\hat A, \hat B] = \hat A\hat B - \hat B\hat A$ is bilinear and antisymmetric. Two identities do most of the work:
$$
[\hat A, \hat B\hat C] = [\hat A, \hat B]\hat C + \hat B[\hat A, \hat C], \qquad [\hat A, [\hat B, \hat C]] + [\hat B, [\hat C, \hat A]] + [\hat C, [\hat A, \hat B]] = 0. \tag{3.14}
$$

▶ **Try it first.** Prove the first by adding and subtracting $\hat B\hat A\hat C$.

**Derivation.**

1. $[\hat A, \hat B\hat C] = \hat A\hat B\hat C - \hat B\hat C\hat A = (\hat A\hat B\hat C - \hat B\hat A\hat C) + (\hat B\hat A\hat C - \hat B\hat C\hat A) = [\hat A, \hat B]\hat C + \hat B[\hat A, \hat C]$.
2. Jacobi: $[\hat A, [\hat B, \hat C]] = \hat A\hat B\hat C - \hat A\hat C\hat B - \hat B\hat C\hat A + \hat C\hat B\hat A$. The other two terms are its cyclic relabelings $A\to B\to C\to A$: $\hat B\hat C\hat A - \hat B\hat A\hat C - \hat C\hat A\hat B + \hat A\hat C\hat B$ and $\hat C\hat A\hat B - \hat C\hat B\hat A - \hat A\hat B\hat C + \hat B\hat A\hat C$. Each of the six orderings $ABC, ACB, BCA, BAC, CAB, CBA$ appears once with $+$ and once with $-$. Sum: $0$. ∎

Compare CM-8 (3.12)–(3.13): the Poisson bracket is bilinear, antisymmetric, obeys Leibniz $\{F, GK\} = \{F, G\}K + G\{F, K\}$ and Jacobi. Every one of these holds for
$$
\{A, B\}\ \longrightarrow\ \frac{1}{i\hbar}[\hat A, \hat B], \qquad \{x, p\} = 1\ \longrightarrow\ [\hat x, \hat p] = i\hbar, \tag{3.15}
$$
term for term (divide (3.14) by $i\hbar$ or $(i\hbar)^2$; the factors cancel in each identity). One place the dictionary breaks: ordering. Classically $xp = px$; quantum mechanically you must choose, and the Groenewold–van Hove theorem (1946) says no assignment of operators to all polynomials in $x, p$ can preserve (3.15) exactly. Exercise A4 shows it works for $[\hat x, \hat p^n] = i\hbar n\hat p^{n-1}$, matching $\{x, p^n\} = np^{n-1}$.

**Theorem 3.6 (Compatible observables).** Let $\hat A$, $\hat B$ be Hermitian on a finite-dimensional space (or have discrete complete eigenbases). (a) If they have a common orthonormal eigenbasis, they commute. (b) If $[\hat A, \hat B] = 0$ and $\hat A$ is nondegenerate, every eigenvector of $\hat A$ is an eigenvector of $\hat B$. (c) If $[\hat A, \hat B] = 0$ in general, a common orthonormal eigenbasis exists.

▶ **Try it first.** For (b), apply $\hat A$ to $\hat B|a\rangle$.

*Proof.* (a) If $\hat A|n\rangle = a_n|n\rangle$, $\hat B|n\rangle = b_n|n\rangle$ for a basis, then $\hat A\hat B|n\rangle = a_nb_n|n\rangle = \hat B\hat A|n\rangle$, and operators that agree on a basis are equal.

(b) Let $\hat A|a\rangle = a|a\rangle$. Then $\hat A(\hat B|a\rangle) = \hat B\hat A|a\rangle = a(\hat B|a\rangle)$, so $\hat B|a\rangle$ is in the $a$-eigenspace. That eigenspace is one-dimensional, so $\hat B|a\rangle = b|a\rangle$.

(c) The argument of (b) shows $\hat B$ maps each eigenspace $E_a$ of $\hat A$ into itself. Restricted to $E_a$, $\hat B$ is still Hermitian, so by LA-8 Theorem 3.12 it has an orthonormal eigenbasis of $E_a$. Each of those vectors is an eigenvector of $\hat A$ (eigenvalue $a$) and of $\hat B$. Collect them over all $a$. ∎

■ **What just happened.** Commuting is the algebraic test for "can be sharp together." In (c), $a$ alone may not label a state, but $(a, b)$ might; if not, add a third commuting observable. A set whose joint eigenvalues label a basis uniquely is a **complete set of commuting observables** (hydrogen's $(n, \ell, m)$ in QM-6).

**Functions of an operator.** For Hermitian $\hat Q = \sum_nq_n|n\rangle\langle n|$ and any function $f$ on its spectrum,
$$
f(\hat Q) = \sum_nf(q_n)|n\rangle\langle n|. \tag{3.16}
$$
It agrees with the direct definition for polynomials, because $\hat Q^k = \sum_nq_n^k|n\rangle\langle n|$ (cross terms vanish by orthonormality), hence with LA-7's power series for $\exp$.

### 3.8 The uncertainty principle, derived in general

**Theorem 3.7 (Robertson, 1929).** For Hermitian $\hat A$, $\hat B$ and any normalized $|\psi\rangle$ in their domains,
$$
\sigma_A^2\,\sigma_B^2\ \ge\ \left(\frac{1}{2i}\langle[\hat A, \hat B]\rangle\right)^2. \tag{3.17}
$$
(The right side is a nonnegative real number: $[\hat A, \hat B]^\dagger = -[\hat A, \hat B]$, so $\langle[\hat A, \hat B]\rangle$ is purely imaginary.)

▶ **Try it first.** Define $|f\rangle = (\hat A - \langle A\rangle)|\psi\rangle$ and $|g\rangle = (\hat B - \langle B\rangle)|\psi\rangle$. You know $\sigma_A^2 = \langle f|f\rangle$ from Theorem 3.5. Apply Cauchy–Schwarz to $\langle f|g\rangle$ and throw away its real part.

**Derivation.**

1. By the proof of Theorem 3.5, $\sigma_A^2 = \langle f|f\rangle$ and $\sigma_B^2 = \langle g|g\rangle$.
2. Cauchy–Schwarz (LA-8 Theorem 3.4): $\langle f|f\rangle\langle g|g\rangle\ge|\langle f|g\rangle|^2$.
3. For any complex $z$, $|z|^2 = (\operatorname{Re}z)^2 + (\operatorname{Im}z)^2\ge(\operatorname{Im}z)^2 = \left(\frac{z - \bar z}{2i}\right)^2$.
4. Compute $\langle f|g\rangle$, using Hermiticity of $\hat A - \langle A\rangle$ to move it to the right: $\langle f|g\rangle = \langle\psi|(\hat A - \langle A\rangle)(\hat B - \langle B\rangle)|\psi\rangle = \langle\hat A\hat B\rangle - \langle A\rangle\langle B\rangle$. Likewise $\overline{\langle f|g\rangle} = \langle g|f\rangle = \langle\hat B\hat A\rangle - \langle A\rangle\langle B\rangle$.
5. So $\langle f|g\rangle - \overline{\langle f|g\rangle} = \langle\hat A\hat B - \hat B\hat A\rangle = \langle[\hat A, \hat B]\rangle$. Combining 1–5 gives (3.17). ∎

■ **What just happened.** The uncertainty principle is Cauchy–Schwarz plus "drop the real part." This is the proof the LA checkpoint previewed. No measurement disturbance entered: (3.17) bounds the spreads of two separate ensembles of measurements on identical states. Keeping the real part gives the stronger Schrödinger–Robertson relation (Exercise B2).

**Corollary (Heisenberg; Kennard 1927).** With $[\hat x, \hat p] = i\hbar$, the right side of (3.17) is $(\hbar/2)^2$:
$$
\sigma_x\sigma_p\ \ge\ \frac\hbar2. \tag{3.18}
$$
QM-2 got this from F3's Fourier uncertainty theorem. Here it needs only the commutator.

**Theorem 3.8 (Minimum-uncertainty states are Gaussians).** Equality in (3.18) holds if and only if $\psi(x) = N\exp\!\big(-\frac{(x - \langle x\rangle)^2}{4\sigma^2} + \frac{i\langle p\rangle x}{\hbar}\big)$ for some $\sigma > 0$.

▶ **Try it first.** Equality needs equality in both step 2 and step 3. What does each say about $|g\rangle$ and $|f\rangle$?

**Derivation.**

1. Equality in Cauchy–Schwarz holds iff $|g\rangle = c|f\rangle$ for some $c\in\mathbb C$ (or $|f\rangle = 0$, which would make $\sigma_x = 0$, impossible for a normalizable state).
2. Then $\langle f|g\rangle = c\sigma_x^2$. Equality in step 3 needs $\operatorname{Re}\langle f|g\rangle = 0$, so $c = ia$ with $a$ real.
3. In the position basis, $|g\rangle = ia|f\rangle$ reads
$$
\big(-i\hbar\partial_x - \langle p\rangle\big)\psi = ia\,(x - \langle x\rangle)\psi. \tag{3.19}
$$
4. Divide by $-i\hbar\psi$: $\psi'/\psi = i\langle p\rangle/\hbar - (a/\hbar)(x - \langle x\rangle)$. Integrate: $\psi = N\exp\!\big(-\frac{a}{2\hbar}(x - \langle x\rangle)^2 + \frac{i\langle p\rangle x}{\hbar}\big)$, absorbing constant phases into $N$.
5. Normalizable iff $a > 0$. Matching $|\psi|^2\propto e^{-(x - \langle x\rangle)^2/2\sigma^2}$ gives $a = \hbar/2\sigma^2$. Conversely, every such Gaussian satisfies (3.19), hence equality. ∎

■ **What just happened.** QM-2 found that Gaussians are optimal by computing; here they are forced. The same argument for general $\hat A$, $\hat B$ says the saturating states solve $(\hat A - \langle A\rangle)\psi = i\lambda(\hat B - \langle B\rangle)\psi$, i.e. they are eigenvectors of the non-Hermitian operator $\hat A - i\lambda\hat B$ (Exercise B3). For $\hat x$ and $\hat p$ that operator is proportional to the oscillator's $\hat a$, so the minimum-uncertainty states are QM-3's coherent states. `project.py` Part 5 checks (3.17) on 2000 random states for random $4\times4$ Hermitian $A$, $B$, finds the floor ratio $1.000000$ by BFGS, and shows an eigenvector of $A - 0.7iB$ saturates to 12 digits.

### 3.9 The generalized Ehrenfest theorem and energy–time uncertainty

**Theorem 3.9 (Generalized Ehrenfest).** If $i\hbar\frac{d}{dt}|\psi\rangle = \hat H|\psi\rangle$ with $\hat H$ Hermitian, then for any $\hat Q$ (possibly depending on $t$ explicitly),
$$
\frac{d}{dt}\langle Q\rangle = \frac{i}{\hbar}\big\langle[\hat H, \hat Q]\big\rangle + \Big\langle\frac{\partial\hat Q}{\partial t}\Big\rangle. \tag{3.20}
$$

*Derivation.* 1. Product rule: $\frac{d}{dt}\langle\psi|\hat Q|\psi\rangle = \langle\dot\psi|\hat Q\psi\rangle + \langle\psi|\partial_t\hat Q|\psi\rangle + \langle\psi|\hat Q\dot\psi\rangle$. 2. $|\dot\psi\rangle = -\frac i\hbar\hat H|\psi\rangle$; in the first slot the constant conjugates: $\langle\dot\psi|\hat Q\psi\rangle = \frac{i}{\hbar}\langle\hat H\psi|\hat Q\psi\rangle = \frac i\hbar\langle\psi|\hat H\hat Q|\psi\rangle$ (by (3.9) with $\hat H^\dagger = \hat H$). 3. The third term is $-\frac i\hbar\langle\psi|\hat Q\hat H|\psi\rangle$. Add. ∎ This is QM-2 (3.11) with the explicit-time term, and with $\hat Q = \hat x, \hat p$ it gives QM-2's Ehrenfest equations again. `project.py` Part 6 checks (3.20) by finite differences of the exact evolution for $\hat x$, $\hat x^2$ (oscillator) and $\sigma_x$, $\sigma_z$ (two-state).

**Theorem 3.10 (Mandelstam–Tamm, 1945).** For time-independent $\hat Q$, put $\hat A = \hat H$, $\hat B = \hat Q$ in (3.17) and use (3.20): $\sigma_H\sigma_Q\ge\frac12|\langle[\hat H, \hat Q]\rangle| = \frac\hbar2\big|\frac{d\langle Q\rangle}{dt}\big|$. Define $\Delta t_Q = \sigma_Q/|d\langle Q\rangle/dt|$, the time for $\langle Q\rangle$ to change by one standard deviation. Then
$$
\Delta E\,\Delta t_Q\ \ge\ \frac\hbar2, \qquad \Delta E = \sigma_H. \tag{3.21}
$$

Read (3.21) for what it is. Time is a parameter, not an operator (Pauli, 1933: a $\hat t$ with $[\hat H, \hat t] = i\hbar$ would make $\hat H$ unbounded below). $\Delta t_Q$ depends on the state and on the chosen $\hat Q$; it is not the spread of any measured distribution. The inequality says a state with sharp energy changes slowly in every respect; a stationary state does not change at all. Popular uses ("virtual particles borrow $\Delta E$ for time $\hbar/\Delta E$") are heuristics that do not follow from (3.21). Exercise B7 shows a two-state system reaching an orthogonal state in exactly the minimum time $\pi\hbar/2\sigma_H$.

### 3.10 Time evolution: the propagator

For time-independent $\hat H$, the Schrödinger equation $i\hbar\frac{d}{dt}|\psi(t)\rangle = \hat H|\psi(t)\rangle$ is a linear ODE with constant coefficients, and LA-7's matrix exponential solves it:
$$
|\psi(t)\rangle = \hat U(t)|\psi(0)\rangle, \qquad \hat U(t) = e^{-i\hat Ht/\hbar} = \sum_ne^{-iE_nt/\hbar}|n\rangle\langle n|, \tag{3.22}
$$
the second form by (3.16) with $\hat H|n\rangle = E_n|n\rangle$. Applying $\hat U$ to $|\psi(0)\rangle = \sum c_n|n\rangle$ gives QM-2/3's general solution
$$
|\psi(t)\rangle = \sum_nc_ne^{-iE_nt/\hbar}|n\rangle, \qquad c_n = \langle n|\psi(0)\rangle. \tag{3.23}
$$
Stationary states are the eigenvectors of $\hat U(t)$: they acquire only a phase, and every expectation value in them is constant.

**Theorem 3.11 ($\hat U$ is unitary).** If $\hat H^\dagger = \hat H$ then $\hat U(t)^\dagger\hat U(t) = \hat 1$, so $\langle\psi(t)|\phi(t)\rangle = \langle\psi(0)|\phi(0)\rangle$; in particular the norm is conserved.

*Proof.* From (3.22), $\hat U^\dagger = \sum_n\overline{e^{-iE_nt/\hbar}}|n\rangle\langle n| = \sum_ne^{+iE_nt/\hbar}|n\rangle\langle n|$, using $E_n$ real (Theorem 3.2). Then $\hat U^\dagger\hat U = \sum_{n,m}e^{iE_nt/\hbar}e^{-iE_mt/\hbar}|n\rangle\langle n|m\rangle\langle m| = \sum_n|n\rangle\langle n| = \hat 1$. ∎ Without Hermiticity the phases would not have modulus one and probability would leak.

**Generators.** Conversely, every continuous one-parameter family of unitaries $\hat U(s)$ with $\hat U(s + s') = \hat U(s)\hat U(s')$ is $e^{-is\hat G}$ for a unique Hermitian $\hat G$: Stone's theorem (1932; in finite dimensions $\hat G = i\,\frac{d\hat U}{ds}\big|_0$, and LA-12 Appendix A does the case of rotations, where the generators are the angular momenta). Near the identity, $\hat U(\varepsilon)\approx\hat 1 - i\varepsilon\hat G$, and unitarity to first order, $(\hat 1 + i\varepsilon\hat G^\dagger)(\hat 1 - i\varepsilon\hat G) = \hat 1 + i\varepsilon(\hat G^\dagger - \hat G) + O(\varepsilon^2)$, forces $\hat G^\dagger = \hat G$. Time evolution is generated by $\hat H/\hbar$; differentiating $|\psi(t + \varepsilon)\rangle = (\hat 1 - i\varepsilon\hat H/\hbar)|\psi(t)\rangle$ is the Schrödinger equation.

**Time-dependent $\hat H(t)$.** $\hat U$ is still unitary, but if $[\hat H(t), \hat H(t')]\ne0$ it is the time-ordered exponential $\mathcal T\exp(-\frac i\hbar\int_0^t\hat H\,dt')$, the Dyson series $\hat 1 - \frac i\hbar\int_0^t\hat H(t_1)dt_1 + (-\frac i\hbar)^2\int_0^tdt_1\int_0^{t_1}dt_2\,\hat H(t_1)\hat H(t_2) + \cdots$, with later times to the left. QM-8's time-dependent perturbation theory is this series truncated.

### 3.11 The Heisenberg picture

Experiments see expectation values, and $\langle\psi(t)|\hat A|\psi(t)\rangle = \langle\psi(0)|\hat U^\dagger\hat A\hat U|\psi(0)\rangle$. Group the $\hat U$'s with the operator instead of the state:
$$
\hat A_H(t) = \hat U(t)^\dagger\hat A\,\hat U(t), \qquad |\psi_H\rangle = |\psi(0)\rangle. \tag{3.24}
$$

**Theorem 3.12 (Heisenberg picture).** (a) $\langle\psi_H|\hat A_H(t)|\psi_H\rangle = \langle\psi(t)|\hat A|\psi(t)\rangle$ for every state and operator: the two pictures give identical predictions. (b) For time-independent $\hat H$,
$$
\frac{d\hat A_H}{dt} = \frac i\hbar[\hat H, \hat A_H] + \Big(\frac{\partial\hat A}{\partial t}\Big)_H. \tag{3.25}
$$
(c) Commutators are preserved: if $[\hat A, \hat B] = \hat C$ then $[\hat A_H(t), \hat B_H(t)] = \hat C_H(t)$; operator products map to products.

▶ **Try it first.** Differentiate $\hat U^\dagger\hat A\hat U$ with $\dot{\hat U} = -\frac i\hbar\hat H\hat U$ and $\dot{\hat U}^\dagger = \frac i\hbar\hat U^\dagger\hat H$.

*Derivation.* (a) is (3.24). (b): 1. Product rule: $\dot{\hat A}_H = \dot{\hat U}^\dagger\hat A\hat U + \hat U^\dagger(\partial_t\hat A)\hat U + \hat U^\dagger\hat A\dot{\hat U}$. 2. Substitute: $= \frac i\hbar\hat U^\dagger\hat H\hat A\hat U - \frac i\hbar\hat U^\dagger\hat A\hat H\hat U + (\partial_t\hat A)_H$. 3. $\hat H$ commutes with $\hat U = e^{-i\hat Ht/\hbar}$, so $\hat U^\dagger\hat H\hat A\hat U = \hat H\hat U^\dagger\hat A\hat U = \hat H\hat A_H$, and $\hat U^\dagger\hat A\hat H\hat U = \hat A_H\hat H$. That is (3.25). (c): $\hat U^\dagger\hat A\hat B\hat U = \hat U^\dagger\hat A\hat U\hat U^\dagger\hat B\hat U$ since $\hat U\hat U^\dagger = \hat 1$. ∎

■ **What just happened.** Compare CM-8 (3.14), $\dot F = \{F, H\} + \partial_tF$. With (3.15), $\{F, H\}\to\frac{1}{i\hbar}[\hat F, \hat H] = \frac i\hbar[\hat H, \hat F]$: (3.25) is Hamilton's equation of motion for an arbitrary observable, with hats. This is Heisenberg's 1925 formulation, in which operators move and states sit still; CM-8's dictionary is now complete. Taking the expectation of (3.25) in the fixed $|\psi_H\rangle$ gives (3.20) back.

**The free particle.** $\hat H = \hat p^2/2m$. By (3.14), $[\hat H, \hat x] = \frac{1}{2m}(\hat p[\hat p, \hat x] + [\hat p, \hat x]\hat p) = -\frac{i\hbar}{m}\hat p$ and $[\hat H, \hat p] = 0$, so $\hat x_H(t) = \hat x + \hat pt/m$: classical free motion with operators (Exercise B4 gets packet spreading from it).

**The oscillator.** $\hat H = \hat p^2/2m + \frac12m\omega^2\hat x^2$. Also $[\hat H, \hat p] = \frac12m\omega^2[\hat x^2, \hat p] = \frac12m\omega^2(2i\hbar\hat x) = i\hbar m\omega^2\hat x$. Then (3.25) gives $\dot{\hat x}_H = \hat p_H/m$ and $\dot{\hat p}_H = -m\omega^2\hat x_H$: Hamilton's equations for the oscillator, a linear ODE system for two operators. With initial values $\hat x$, $\hat p$ at $t = 0$,
$$
\hat x_H(t) = \hat x\cos\omega t + \frac{\hat p}{m\omega}\sin\omega t, \qquad \hat p_H(t) = \hat p\cos\omega t - m\omega\,\hat x\sin\omega t. \tag{3.26}
$$
Check: differentiate and substitute; and $[\hat x_H(t), \hat p_H(t)] = (\cos^2\omega t + \sin^2\omega t)\,i\hbar = i\hbar$ as (c) requires. In ladder form, $\hat a_H(t) = \hat ae^{-i\omega t}$. `project.py` Part 2 checks (3.26) as a matrix identity, $\hat U^\dagger\hat X\hat U = \hat X\cos t + \hat P\sin t$, in the truncated $N = 60$ number basis.

**The interaction picture (one paragraph; QM-8 uses it).** Split $\hat H = \hat H_0 + \hat V(t)$ with $\hat H_0$ solvable. Move by $\hat H_0$ only: $|\psi_I(t)\rangle = e^{i\hat H_0t/\hbar}|\psi(t)\rangle$ and $\hat V_I(t) = e^{i\hat H_0t/\hbar}\hat V(t)e^{-i\hat H_0t/\hbar}$. Differentiating, $i\hbar\frac{d}{dt}|\psi_I\rangle = -\hat H_0|\psi_I\rangle + e^{i\hat H_0t/\hbar}(\hat H_0 + \hat V)|\psi\rangle$, and the $\hat H_0$ terms cancel because $\hat H_0$ commutes with its own exponential:
$$
i\hbar\frac{d}{dt}|\psi_I(t)\rangle = \hat V_I(t)|\psi_I(t)\rangle. \tag{3.27}
$$
States move only because of $\hat V$, slowly if $\hat V$ is weak: the basis of perturbation theory. Exercise C1 uses it for a driven spin.

### 3.12 Two-state systems: the complete worked instance

The smallest nontrivial Hilbert space is $\mathbb C^2$. Any $2\times2$ Hermitian matrix is (LA-8 Proposition 3.22: $\{I, \sigma_x, \sigma_y, \sigma_z\}$ is a real basis of the Hermitian $2\times2$ matrices)
$$
\hat H = E_0\hat 1 + \mathbf h\cdot\boldsymbol\sigma = \begin{pmatrix}E_0 + h_z & h_x - ih_y\\ h_x + ih_y & E_0 - h_z\end{pmatrix}, \qquad E_\pm = E_0\pm|\mathbf h|, \tag{3.28}
$$
with $E_0$ and $\mathbf h\in\mathbb R^3$ in energy units. The eigenvalues: $(\mathbf h\cdot\boldsymbol\sigma)^2 = |\mathbf h|^2\hat 1$ (expand using $\sigma_j\sigma_k + \sigma_k\sigma_j = 2\delta_{jk}$), so the eigenvalues of $\mathbf h\cdot\boldsymbol\sigma$ are $\pm|\mathbf h|$. The eigenvectors are the spin states along $\pm\hat{\mathbf h}$: if $\hat{\mathbf h}$ has polar angles $(\theta, \phi)$, $|{+}\rangle = (\cos\frac\theta2, e^{i\phi}\sin\frac\theta2)^\top$ (LA-8 (3.15)). $E_0$ adds only a global phase; drop it.

A state has a Bloch vector $\mathbf n = \langle\boldsymbol\sigma\rangle$, $|\mathbf n| = 1$ (QM-1).

**Theorem 3.13 (Precession).** Under (3.28),
$$
\frac{d\mathbf n}{dt} = \frac2\hbar\,\mathbf h\times\mathbf n: \tag{3.29}
$$
the Bloch vector rotates rigidly about $\hat{\mathbf h}$ at angular frequency $2|\mathbf h|/\hbar = (E_+ - E_-)/\hbar$.

▶ **Try it first.** Use (3.20) with $\hat Q = \sigma_j$ and $[\sigma_k, \sigma_j] = 2i\epsilon_{kjl}\sigma_l$ (QM-1 Proposition 3.7).

*Derivation.* 1. $\frac{d}{dt}\langle\sigma_j\rangle = \frac i\hbar\sum_kh_k\langle[\sigma_k, \sigma_j]\rangle = \frac i\hbar\sum_{k,l}h_k\,2i\epsilon_{kjl}\langle\sigma_l\rangle = -\frac2\hbar\sum_{k,l}\epsilon_{kjl}h_kn_l$. 2. $\epsilon_{kjl} = -\epsilon_{jkl}$, so this is $\frac2\hbar\sum_{k,l}\epsilon_{jkl}h_kn_l = \frac2\hbar(\mathbf h\times\mathbf n)_j$. 3. $\dot{\mathbf n} = \boldsymbol\Omega\times\mathbf n$ with constant $\boldsymbol\Omega = 2\mathbf h/\hbar$ is rigid rotation about $\boldsymbol\Omega$ at rate $|\boldsymbol\Omega|$ (CM-7's angular velocity). ∎

■ **What just happened.** All two-level dynamics is a gyroscope. Larmor precession, ammonia, neutrino flavors, MRI, and QC-1's single-qubit gate $e^{-i\theta\hat{\mathbf n}\cdot\boldsymbol\sigma/2}$ (a rotation of the Bloch sphere by $\theta$ about $\hat{\mathbf n}$) are this one theorem. QM-8's Rabi oscillations are §3.12 in a rotating frame (Exercise C1).

**The flip probability.** Start in $|{\uparrow}\rangle$ ($\mathbf n = \hat{\mathbf z}$) with $\mathbf h = (\Delta, 0, \varepsilon)$, i.e. $\hat H = \varepsilon\sigma_z + \Delta\sigma_x$. Rotating $\hat{\mathbf z}$ about $\hat{\mathbf h}$ by angle $\alpha = 2|\mathbf h|t/\hbar$: the component along $\hat{\mathbf h}$, $\cos\chi$ with $\sin\chi = \Delta/|\mathbf h|$, is fixed; the perpendicular part, of length $\sin\chi$, turns. So $n_z(t) = \cos^2\chi + \sin^2\chi\cos\alpha$. The probability of $|{\downarrow}\rangle$ is $(1 - n_z)/2 = \sin^2\chi\,\frac{1 - \cos\alpha}{2}$:
$$
P_{\uparrow\to\downarrow}(t) = \frac{\Delta^2}{\varepsilon^2 + \Delta^2}\,\sin^2\!\Big(\frac{\sqrt{\varepsilon^2 + \Delta^2}\,t}{\hbar}\Big). \tag{3.30}
$$
This is the **Rabi formula with detuning**: complete flipping only when $\varepsilon = 0$ (the two basis states have equal energy), and a Lorentzian $\Delta^2/(\varepsilon^2 + \Delta^2)$ maximum otherwise.

**Three instances.**

*Spin in a magnetic field.* $\hat H = -\gamma\hat{\mathbf S}\cdot\mathbf B = -\frac{\gamma\hbar}{2}\mathbf B\cdot\boldsymbol\sigma$, so $\mathbf h = -\gamma\hbar\mathbf B/2$ and (3.29) says the spin precesses at $2|\mathbf h|/\hbar = |\gamma|B$, the Larmor frequency. For the electron $|\gamma|/2\pi = 28.025$ GHz/T; for the proton $42.577$ MHz/T, so protons in a 1.5 T MRI magnet precess at about 64 MHz.

*The ammonia molecule.* NH$_3$ is a pyramid, and the nitrogen can sit above or below the plane of hydrogens: two states $|L\rangle$, $|R\rangle$ of equal energy $E_0$, coupled by tunneling through the barrier: $\hat H = E_0\hat 1 - A\sigma_x$. The eigenstates are $(|L\rangle\pm|R\rangle)/\sqrt2$ split by $2A$; the measured inversion line is $\nu = 23.87$ GHz so $2A = h\nu = 9.87\times10^{-5}$ eV. A molecule prepared in $|L\rangle$ is in $|R\rangle$ after $1/(2\nu) = 20.95$ ps, by (3.30) with $\varepsilon = 0$. The two eigenstates respond oppositely to an electric field, which is how Gordon, Zeiger and Townes separated them in 1954 to build the first maser (Feynman vol. III ch. 8–9).

*Neutrino oscillations.* Flavor states ($\nu_e$, $\nu_\mu$) are superpositions of mass states with energies $E_i\approx pc + m_i^2c^4/2pc$. In the flavor basis this is (3.28) with $|\mathbf h| = \Delta m^2c^4/4E$ and mixing angle $\theta$ (Exercise B5); (3.30) becomes $P_{\nu_\mu\to\nu_e} = \sin^22\theta\sin^2(\Delta m^2c^4L/4\hbar cE)$. With $\Delta m^2\approx2.5\times10^{-3}$ eV$^2$ and $E = 1$ GeV the first maximum is at $L\approx496$ km. Super-Kamiokande's 1998 atmospheric-neutrino data fit this formula: neutrinos have mass.

### 3.13 Symmetries and conservation (CM-4 completed)

A **symmetry** is a unitary $\hat U_s$ with $[\hat U_s, \hat H] = 0$, so $\hat U_s\hat U(t) = \hat U(t)\hat U_s$. For a continuous family $\hat U_s(\theta) = e^{-i\theta\hat G}$ (Stone, §3.10), differentiating $[\hat U_s(\theta), \hat H] = 0$ at $\theta = 0$ gives $[\hat G, \hat H] = 0$.

**Theorem 3.14 (Noether, quantum version).** If $[\hat G, \hat H] = 0$ and $\hat G$ has no explicit time dependence, then $\frac{d}{dt}\langle G\rangle = 0$ in every state; indeed $\hat G_H(t) = \hat G$, and the whole probability distribution of $G$ is constant. *Proof.* (3.20) or (3.25) with $[\hat H, \hat G] = 0$; for the distribution, every function $f(\hat G)$ also commutes with $\hat H$ (use (3.16) and Theorem 3.6). ∎ Compare CM-4 Theorem 3.3; here the conserved charge is literally the symmetry's generator.

**Translations are generated by $\hat p$.** Define $\hat T(a) = e^{-i\hat pa/\hbar}$.

▶ **Try it first.** Expand the exponential in powers of $a$ and act on $\psi(x)$ with $\hat p = -i\hbar\partial_x$. Recognize the result.

1. $\hat T(a)\psi = \sum_{n\ge0}\frac{1}{n!}\big(-\frac{ia}{\hbar}\big)^n(-i\hbar)^n\partial_x^n\psi = \sum_n\frac{(-a)^n}{n!}\psi^{(n)}(x)$, since $(-i)(-i) = -1$.
2. That is the Taylor series of $\psi$ about $x$ evaluated at $x - a$:
$$
\big(e^{-i\hat pa/\hbar}\psi\big)(x) = \psi(x - a). \tag{3.31}
$$ ∎

(For analytic $\psi$; for general $\psi$, the momentum basis gives the same result without series: $\hat T(a)$ multiplies $\Phi(p)$ by $e^{-ipa/\hbar}$, which by (3.7) is a shift.) $[\hat T(a), \hat H] = 0$ for all $a$ iff $[\hat p, \hat H] = 0$ iff $V$ is constant; then momentum is conserved. Time translation is generated by $\hat H$ ((3.22)) and rotations by $\hat{\mathbf L}$ (QM-5: $e^{-i\phi\hat L_z/\hbar}$). **Parity** $\hat\Pi\psi(x) = \psi(-x)$ is a discrete symmetry: $\hat\Pi^2 = \hat 1$, eigenvalues $\pm1$, Hermitian and unitary, with $\hat\Pi\hat x\hat\Pi = -\hat x$ and $\hat\Pi\hat p\hat\Pi = -\hat p$, so $[\hat\Pi, \hat H] = 0$ iff $V$ is even.

**Theorem 3.15 (Symmetry implies degeneracy).** If $[\hat U_s, \hat H] = 0$ and $\hat H|n\rangle = E_n|n\rangle$, then $\hat U_s|n\rangle$ is also an eigenvector with eigenvalue $E_n$. If two symmetries $\hat U_1$, $\hat U_2$ commute with $\hat H$ but not with each other, some energy level is degenerate.

*Proof.* $\hat H\hat U_s|n\rangle = \hat U_s\hat H|n\rangle = E_n\hat U_s|n\rangle$. For the second claim: if every level were nondegenerate, each $|n\rangle$ would be an eigenvector of both $\hat U_1$ and $\hat U_2$ (Theorem 3.6(b), whose proof works for unitary $\hat U_s$ unchanged), so $\hat U_1$ and $\hat U_2$ would be diagonal in the same basis and commute. ∎

This is the CM-4 promise kept: hydrogen's $(2\ell + 1)$-fold $m$-degeneracy is rotational symmetry (QM-5/6), its extra $\ell$-degeneracy the Runge–Lenz symmetry of CM-5. In one dimension bound states are nondegenerate (QM-3), so parity only labels them. `project.py` Part 8 shows the labels on a double well and a genuine degeneracy on a ring, where translation and parity commute with $\hat H$ but not with each other, forcing $E(k) = E(-k)$.

**Theorem 3.16 (Bloch).** If $V(x + a) = V(x)$, the eigenstates of $\hat H$ can be chosen with $\psi(x + a) = e^{ika}\psi(x)$, i.e. $\psi = e^{ikx}u(x)$ with $u$ periodic. *Proof.* $[\hat T(a), \hat H] = 0$, and $\hat T(a)$ is unitary, so the argument of Theorem 3.6(c) gives common eigenvectors. A unitary's eigenvalues have modulus one, so $\hat T(a)\psi = e^{-ika}\psi$ for some real $k$, i.e. $\psi(x - a) = e^{-ika}\psi(x)$; replace $x$ by $x + a$. Set $u = e^{-ikx}\psi$; then $u(x + a) = e^{-ikx}e^{-ika}e^{ika}\psi(x) = u(x)$. ∎ QM-3 stated this for the Dirac comb; here it is two lines.

**The density operator (preview; QM-7 does it).** For a pure state, $\hat\rho = |\psi\rangle\langle\psi|$ and $\langle Q\rangle = \langle\psi|\hat Q|\psi\rangle = \operatorname{Tr}(\hat\rho\hat Q)$ (insert (3.3) and use $\operatorname{Tr}(|\psi\rangle\langle\phi|) = \langle\phi|\psi\rangle$). A source that emits $|\psi_k\rangle$ with probability $w_k$ is described by
$$
\hat\rho = \sum_kw_k|\psi_k\rangle\langle\psi_k|, \qquad \langle Q\rangle = \operatorname{Tr}(\hat\rho\hat Q), \qquad \operatorname{Tr}\hat\rho^2\le1, \tag{3.32}
$$
with equality only for a pure state. Mixed states are not kets; half of an entangled pair is one (QM-7), and QC-7/8's noise models are written in them. Exercise C3 computes one.

## 4. Worked examples

### Example 4.1: The box in Dirac notation

QM-3's infinite well on $[0, L]$ (QM-3 calls the width $a$; here $a$ is taken by the ladder operator) has $\langle x|n\rangle = \sqrt{2/L}\sin(n\pi x/L)$ and $E_n = n^2\pi^2\hbar^2/2mL^2$. Take $|\psi\rangle$ with $\langle x|\psi\rangle = \sqrt{30/L^5}\,x(L - x)$.

*Components.* By (3.5) and (3.4), $c_n = \langle n|\psi\rangle = \int_0^L\langle n|x\rangle\langle x|\psi\rangle\,dx = \sqrt{60}/L^3\int_0^L\sin(n\pi x/L)\,x(L - x)\,dx$. Two integrations by parts give $\int_0^L x(L - x)\sin(n\pi x/L)dx = \frac{2L^3}{(n\pi)^3}(1 - (-1)^n)$, so $c_n = \frac{8\sqrt{15}}{(n\pi)^3}$ for odd $n$ and $0$ for even $n$. (The even states are odd about $L/2$ and $\psi$ is even: a parity selection rule, §3.13.)

*Probabilities* (3.11): $|c_1|^2 = 960/\pi^6 = 0.998555$, $|c_3|^2 = 1.370\times10^{-3}$, $|c_5|^2 = 6.39\times10^{-5}$. Parseval: $\sum_{\text{odd}}960/(n\pi)^6 = \frac{960}{\pi^6}\cdot\frac{\pi^6}{960} = 1$, using $\sum_{\text{odd}}n^{-6} = \pi^6/960$.

*Energy two ways* (3.13). From the components: $\langle H\rangle = \sum|c_n|^2E_n = \frac{960}{\pi^6}\frac{\pi^2\hbar^2}{2mL^2}\sum_{\text{odd}}n^{-4} = \frac{480\hbar^2}{\pi^4mL^2}\cdot\frac{\pi^4}{96} = \frac{5\hbar^2}{mL^2}$. Directly: $\langle\psi|\hat H|\psi\rangle = -\frac{\hbar^2}{2m}\int_0^L\psi\psi''dx = -\frac{\hbar^2}{2m}\cdot\frac{30}{L^5}\int_0^Lx(L - x)(-2)\,dx = \frac{30\hbar^2}{mL^5}\cdot\frac{L^3}{6} = \frac{5\hbar^2}{mL^2}$. In units of $E_1$: $\langle H\rangle/E_1 = 10/\pi^2 = 1.01321$. The parabola is 99.86% ground state with energy 1.3% above $E_1$ (QM-8's variational principle exploits this). `project.py` Part 1 does both computations in SymPy and sums the series.

### Example 4.2: Heisenberg's matrices for the oscillator

With $\hat a|n\rangle = \sqrt n|n - 1\rangle$, $\hat a^\dagger|n\rangle = \sqrt{n + 1}|n + 1\rangle$ (QM-3's $a_-$ and $a_+$; $a_+ = a_-^\dagger$ is now a statement about adjoints) and $\hat x = \sqrt{\hbar/2m\omega}(\hat a + \hat a^\dagger)$, $\hat p = i\sqrt{m\hbar\omega/2}(\hat a^\dagger - \hat a)$, the matrix elements (3.8) are
$$
\langle m|\hat x|n\rangle = \sqrt{\tfrac{\hbar}{2m\omega}}\big(\sqrt n\,\delta_{m,n-1} + \sqrt{n + 1}\,\delta_{m,n+1}\big), \qquad \langle m|\hat p|n\rangle = i\sqrt{\tfrac{m\hbar\omega}{2}}\big(\sqrt{n + 1}\,\delta_{m,n+1} - \sqrt n\,\delta_{m,n-1}\big).
$$
Truncated to $n = 0, 1, 2, 3$ (units $\hbar = m = \omega = 1$):
$$
X = \frac{1}{\sqrt2}\begin{pmatrix}0&1&0&0\\1&0&\sqrt2&0\\0&\sqrt2&0&\sqrt3\\0&0&\sqrt3&0\end{pmatrix}, \qquad P = \frac{i}{\sqrt2}\begin{pmatrix}0&-1&0&0\\1&0&-\sqrt2&0\\0&\sqrt2&0&-\sqrt3\\0&0&\sqrt3&0\end{pmatrix}.
$$
Both are Hermitian; these are Heisenberg's 1925 matrices.

*The commutator.* $[\hat x, \hat p] = \frac{i\hbar}{2}\big([\hat a + \hat a^\dagger, \hat a^\dagger - \hat a]\big) = \frac{i\hbar}{2}\big([\hat a, \hat a^\dagger] - [\hat a^\dagger, \hat a]\big) = i\hbar[\hat a, \hat a^\dagger]$. For the truncated $N\times N$ matrices, $\hat a\hat a^\dagger = \operatorname{diag}(1, 2, \ldots, N - 1, 0)$ and $\hat a^\dagger\hat a = \operatorname{diag}(0, 1, \ldots, N - 1)$, so
$$
[X, P] = i\hbar\operatorname{diag}(1, 1, \ldots, 1, -(N - 1)).
$$
For $N = 4$: $i\hbar\operatorname{diag}(1, 1, 1, -3)$. For $N = 60$ the last entry is $-59i\hbar$. This is no accident: $\operatorname{Tr}[X, P] = 0$ for finite matrices while $\operatorname{Tr}(i\hbar I_N)\ne0$, so **no finite matrices satisfy $[X, P] = i\hbar I$**, and truncation dumps the error into the edge. Low-lying physics is fine. Similarly $\frac12(P^2 + X^2)$ built from the truncated matrices equals $\operatorname{diag}(n + \frac12)$ except in the last entry ($29.5$ instead of $59.5$ for $N = 60$), which is why the project builds $\hat H$ as $\operatorname{diag}(n + \frac12)$ directly.

### Example 4.3: A three-level observable

Let $\hat Q = \begin{pmatrix}0&-i&0\\ i&0&-i\\ 0&i&0\end{pmatrix}$ (it is $\sqrt2\hat S_y/\hbar$ for spin 1; QM-5) and $|\psi\rangle = \frac{1}{\sqrt2}(1, i, 0)^\top$.

*Eigenvalues.* Expand $\det(\hat Q - q)$ along the first row: $-q\big[(-q)^2 - (-i)(i)\big] - (-i)\big[i(-q) - 0\big] = -q(q^2 - 1) - q\cdot(-1)\cdot(-1)\cdot(-1)$, i.e. $-q^3 + q + q = -q^3 + 2q$. Roots $q = -\sqrt2, 0, \sqrt2$.

*Eigenvectors.* Solving $(\hat Q - q)\mathbf v = 0$: $|0\rangle = \frac{1}{\sqrt2}(1, 0, 1)$, $|{\pm\sqrt2}\rangle = \frac12(\mp i, \sqrt2, \pm i)$. Check $q = \sqrt2$ on $(-i, \sqrt2, i)$ row by row: $-i\cdot\sqrt2 = \sqrt2\cdot(-i)$ ✓; $i(-i) - i(i) = 2 = \sqrt2\cdot\sqrt2$ ✓; $i\sqrt2 = \sqrt2\cdot i$ ✓.

*Probabilities* (3.11), remembering to conjugate the bra:
- $\langle0|\psi\rangle = \frac12(1\cdot1 + 0 + 0) = \frac12$, $P(0) = \frac14$.
- $\langle\sqrt2|\psi\rangle = \frac{1}{2\sqrt2}\big(\overline{-i}\cdot1 + \sqrt2\cdot i\big) = \frac{i(1 + \sqrt2)}{2\sqrt2}$, $P(\sqrt2) = \frac{(1 + \sqrt2)^2}{8} = \frac{3 + 2\sqrt2}{8} = 0.72855$.
- $\langle-\sqrt2|\psi\rangle = \frac{1}{2\sqrt2}(\overline{i} + \sqrt2 i) = \frac{i(\sqrt2 - 1)}{2\sqrt2}$, $P(-\sqrt2) = \frac{3 - 2\sqrt2}{8} = 0.02145$.

Sum: $\frac14 + \frac68 = 1$ ✓.

*$\langle Q\rangle$ two ways.* From (3.13): $\sqrt2(0.72855 - 0.02145) = \sqrt2\cdot\frac{4\sqrt2}{8} = 1$. Directly: $\hat Q|\psi\rangle = \frac{1}{\sqrt2}(-i\cdot i, i, i\cdot i) = \frac{1}{\sqrt2}(1, i, -1)$, and $\langle\psi|\hat Q\psi\rangle = \frac12(1\cdot1 + \bar i\cdot i + 0) = \frac12(1 + 1) = 1$ ✓. Also $\langle Q^2\rangle = 2(P(\sqrt2) + P(-\sqrt2)) = \frac32$, so $\sigma_Q^2 = \frac12$.

*Post-measurement state.* If the outcome is $\sqrt2$, (3.12) leaves $|\sqrt2\rangle = \frac12(-i, \sqrt2, i)$ (up to the irrelevant phase $i$ of the projection). Measuring $\hat Q$ again gives $\sqrt2$ with certainty. `project.py` Part 1 checks every number and Part 7 samples $10^4$ outcomes: frequencies $(0.0222, 0.2422, 0.7356)$ against $(0.0214, 0.25, 0.7286)$, all within four standard errors.

### Example 4.4: The spin uncertainty product, derived

For spin-1/2, $[\hat S_x, \hat S_y] = i\hbar\hat S_z$ (QM-1 Proposition 3.7), so (3.17) reads $\sigma_{S_x}\sigma_{S_y}\ge\frac\hbar2|\langle S_z\rangle|$. A state with Bloch vector $\mathbf n$ has $\langle\hat S_j\rangle = \frac\hbar2n_j$ and $\langle\hat S_j^2\rangle = \frac{\hbar^2}{4}$ (since $\sigma_j^2 = I$), so $\sigma_{S_j}^2 = \frac{\hbar^2}{4}(1 - n_j^2)$. The inequality becomes $(1 - n_x^2)(1 - n_y^2)\ge n_z^2 = 1 - n_x^2 - n_y^2$, i.e. $n_x^2n_y^2\ge0$, QM-1 Exercise B6's result, now derived from the general theorem.

Numbers: $\theta = \pi/3$, $\phi = \pi/4$, so $\mathbf n = (0.61237, 0.61237, 0.5)$ and $1 - n_x^2 = 1 - 0.375 = 0.625$. Then $\sigma_{S_x}\sigma_{S_y} = \frac{\hbar^2}{4}(0.625) = 0.15625\hbar^2$ against the bound $\frac\hbar2\cdot\frac\hbar2(0.5) = 0.125\hbar^2$. Equality needs $n_x = 0$ or $n_y = 0$.

### Example 4.5: The two-state system $\hat H = \varepsilon\sigma_z + \Delta\sigma_x$

Units $\hbar = 1$, $\varepsilon = 0.3$, $\Delta = 0.4$, so $\mathbf h = (0.4, 0, 0.3)$, $|\mathbf h| = 0.5$.

*Spectrum* (3.28): $E_\pm = \pm0.5$, with $|{+}\rangle$ along $\hat{\mathbf h} = (0.8, 0, 0.6)$, i.e. $\theta = 53.13°$: $|{+}\rangle = (\cos26.57°, \sin26.57°) = (2, 1)/\sqrt5$. Check: $\hat H(2, 1)^\top = (0.6 + 0.4, 0.8 - 0.3) = 0.5(2, 1)$ ✓.

*Precession* (3.29): starting from $\mathbf n = \hat{\mathbf z}$, the Bloch vector circles $\hat{\mathbf h}$ at angular frequency $2|\mathbf h| = 1$, period $2\pi = 6.2832$, on a cone of half-angle $\chi = 53.13°$.

*Flip probability* (3.30): $P(t) = 0.64\sin^2(0.5t)$. Maximum $0.64$ first at $t = \pi$; at $t = 1$, $P = 0.64\sin^2(0.5) = 0.64\times0.229849 = 0.147103$. Scanning $\varepsilon$, the maximum traces the Lorentzian $\Delta^2/(\varepsilon^2 + \Delta^2)$ of half-width $\Delta$, the resonance curve of spectroscopy. `project.py` Part 4 checks all of this and plots the precession on the Bloch sphere (`two_state.png`).

### Example 4.6: The oscillator in the Heisenberg picture, coherent initial state

Units $\hbar = m = \omega = 1$. A coherent state is an eigenvector of the (non-Hermitian) lowering operator, $\hat a|\alpha\rangle = \alpha|\alpha\rangle$ with $\alpha\in\mathbb C$; in the number basis $\langle n|\alpha\rangle = e^{-|\alpha|^2/2}\alpha^n/\sqrt{n!}$ (check by applying $\hat a$; QM-3's project evolves one). With with $\hat x = (\hat a + \hat a^\dagger)/\sqrt2$, $\langle x\rangle_0 = \sqrt2\operatorname{Re}\alpha$, $\langle p\rangle_0 = \sqrt2\operatorname{Im}\alpha$. Take $\alpha = 1 + i$: $\langle x\rangle_0 = \langle p\rangle_0 = \sqrt2 = 1.41421$.

By (3.26), with no need to evolve the state, $\langle x(t)\rangle = \langle x\rangle_0\cos t + \langle p\rangle_0\sin t = \sqrt2(\cos t + \sin t) = 2\cos(t - \pi/4)$: amplitude $2 = \sqrt2|\alpha|$, phase $\arg\alpha$. The width: $\hat x_H(t) - \langle x(t)\rangle = (\hat x - \langle x\rangle_0)\cos t + (\hat p - \langle p\rangle_0)\sin t$, so $\sigma_x^2(t) = \sigma_x^2\cos^2t + \sigma_p^2\sin^2t + 2\operatorname{cov}(x, p)\sin t\cos t$ with $\operatorname{cov} = \frac12\langle\{\hat x - \langle x\rangle, \hat p - \langle p\rangle\}\rangle$. For a coherent state $\sigma_x^2 = \sigma_p^2 = \frac12$ and $\operatorname{cov} = 0$, so $\sigma_x(t) = 1/\sqrt2 = 0.70711$ forever: the packet oscillates like a classical particle without spreading. `project.py` Part 2 builds $|\alpha\rangle$ in the $N = 60$ basis, checks $\hat a|\alpha\rangle = \alpha|\alpha\rangle$ there, and checks $\langle x(t)\rangle$ and $\sigma_x$ against these formulas.

### Example 4.7: The translation operator on a Gaussian

Let $\psi(x) = (2\pi s^2)^{-1/4}e^{-x^2/4s^2}$ (so $\sigma_x = s$ and $\sigma_p = \hbar/2s$). By (3.31), $\hat T(a)\psi = \psi(x - a)$: the same Gaussian centered at $a$, and $\langle x\rangle$ goes from $0$ to $a$. In the Heisenberg-like language, $\hat T^\dagger\hat x\hat T = \hat x + a$ (Exercise B6).

The overlap $\langle\psi|\hat T(a)|\psi\rangle$ can be computed two ways. In position space: $\int\psi(x)\psi(x - a)dx = (2\pi s^2)^{-1/2}\int e^{-[x^2 + (x - a)^2]/4s^2}dx$. Complete the square, $x^2 + (x - a)^2 = 2(x - a/2)^2 + a^2/2$, so the integral is $(2\pi s^2)^{-1/2}e^{-a^2/8s^2}\sqrt{2\pi s^2} = e^{-a^2/8s^2}$. In momentum space: $\langle\psi|e^{-i\hat pa/\hbar}|\psi\rangle = \int|\Phi(p)|^2e^{-ipa/\hbar}dp$, the characteristic function of a Gaussian with variance $\sigma_p^2 = \hbar^2/4s^2$, which is $e^{-a^2\sigma_p^2/2\hbar^2} = e^{-a^2/8s^2}$ ✓. For $s = 1$, $a = 2$: $e^{-1/2} = 0.60653$.

`project.py` Part 8 builds $\hat T = \exp(-i\hat Pa/\hbar)$ by `expm` on a 256-point grid with the spectral momentum $\hat P = F^\dagger\operatorname{diag}(\hbar k)F$: the shift is exact to $10^{-10}$ and the overlap is $0.60653$. With the central-difference $\hat p$ (eigenvalues $\hbar\sin(k\Delta x)/\Delta x$) the shift errs by $\sim10^{-3}$, because Fourier components travel different distances.

## 5. Why it matters and where it connects

- **QM-5** derives angular momentum spectra from $[\hat L_x, \hat L_y] = i\hbar\hat L_z$ with (3.14), Theorem 3.6 (common eigenstates of $\hat L^2$, $\hat L_z$) and ladder operators. **QM-6:** hydrogen's $(n, \ell, m)$ are a complete set of commuting observables; its degeneracies are Theorem 3.15 for rotations and Runge–Lenz (CM-5 Theorem 3.10).
- **QM-7** builds $\mathcal H_A\otimes\mathcal H_B$ on LA-12 §3.12 with outer products (3.2) and the density operator (3.32).
- **QM-8.** The interaction picture (3.27) and Dyson series are time-dependent perturbation theory; Rabi oscillation is (3.30) in a rotating frame; the variational principle is LA-8's Rayleigh quotient (Theorem 3.17), with Example 4.1 as a first instance.
- **QC-1, QC-2, QC-5, QC-9.** A qubit is §3.12's $\mathbb C^2$; a gate is $e^{-i\theta\hat{\mathbf n}\cdot\boldsymbol\sigma/2}$; measurement is (3.11)–(3.12); the quantum Fourier transform is Part 3's DFT matrix on $2^n$ amplitudes; Hamiltonian simulation approximates (3.22).
- **CM-8, CM-4.** (3.25) is CM-8 (3.14) under (3.15); Theorem 3.14 is CM-4 Theorem 3.3 with the generator promoted to an operator.

## 6. Exercises

See [exercises.md](exercises.md). **A1–A5** drill bra-ket arithmetic, adjoints, projectors, $[\hat x, \hat p^n]$, and a degenerate measurement. **B1–B8** are the module; essential: B1 (compatibility, degenerate case), B2 (Schrödinger–Robertson), B4 (spreading from the Heisenberg picture), B6 (translations, parity, Bloch). B7 is a coding exercise. **C1–C3**: the driven spin and the rotating-wave approximation (code), why no finite matrices satisfy $[X, P] = i\hbar$, a mixed state.

## 7. Coding project

`project.py` (about 3 s) is a 100-line library plus applications. `State` holds $\langle n|\psi\rangle$; `Operator` holds $\langle m|\hat Q|n\rangle$ and implements `dag` (3.9), `expect` (3.13), `sigma` (Theorem 3.5), `born` ((3.11): eigenvalues grouped by degeneracy, probabilities, projectors), `measure` (sample, then collapse by (3.12)) and `U` ((3.22) by `scipy.linalg.expm`). Parts: (1) Examples 4.1, 4.3; (2) the $N = 60$ oscillator, its $-59i$ commutator edge, (3.26) as a matrix identity, Example 4.6; (3) the DFT as change of basis; (4) the two-state system against (3.29)–(3.30) (`two_state.png`); (5) Theorem 3.7 on random matrices, the BFGS floor, Example 4.4; (6) (3.20) by finite differences; (7) $10^4$ measurements (`measurements.png`); (8) translation, parity, the ring degeneracy; (9) every number in the solutions.

**Build it yourself first:** write `born` and `measure` for a possibly degenerate observable. The trap: `eigh` returns an arbitrary basis of each eigenspace, so group eigenvalues within a tolerance and build $\hat P_q$ from the whole group. Then do "Your turn": the interaction picture for a driven spin, $\operatorname{Tr}(\hat\rho\hat Q)$ and purity, and Bloch's theorem on an 8-cell comb (reference checks follow each prompt).

## 8. Where you'll struggle

- **If LA-8 and LA-12 are solid, this module is a relief; if not, it is where QM becomes incomprehensible.** Test: can you prove closed-book that Hermitian matrices have real eigenvalues and orthogonal eigenvectors, and say why $\langle\phi|$ is the conjugate transpose of $|\phi\rangle$? If not, redo LA-8 §3.3–3.4 and LA-12 §3.6 first.
- **The conjugate-linear slot.** $\langle c\phi|\psi\rangle = \bar c\langle\phi|\psi\rangle$. Forgetting it causes most wrong $i$'s. Fix: conjugate every number the moment you turn a ket into a bra, and check that probabilities sum to 1.
- **The direction of the adjoint.** $\langle\phi|\hat A|\psi\rangle = \langle\hat A^\dagger\phi|\psi\rangle$; writing $\langle\hat A\phi|\psi\rangle$ is wrong unless $\hat A$ is Hermitian. Fix: in coordinates it is $\boldsymbol\phi^\dagger A\boldsymbol\psi = (A^\dagger\boldsymbol\phi)^\dagger\boldsymbol\psi$.
- **The continuous spectrum is swept under the rug** by every undergraduate course, this one included: $|x\rangle$ and $|p\rangle$ are not in $L^2$. The rug is called a rigged Hilbert space (§3.1). When a calculation produces $\delta(0)$, that is the rug showing; go back to normalizable packets.
- **The energy–time relation is not like the others.** $\Delta t$ in (3.21) is not a standard deviation and there is no time operator. Whenever you use it, say which $\hat Q$ defines $\Delta t$.
- **Heisenberg-picture signs.** $\hat A_H = \hat U^\dagger\hat A\hat U$ (not $\hat U\hat A\hat U^\dagger$), and $\dot{\hat A}_H = +\frac i\hbar[\hat H, \hat A_H]$ while $\frac{d}{dt}|\psi\rangle = -\frac i\hbar\hat H|\psi\rangle$. Check on the oscillator: $\dot{\hat x}_H = +\hat p_H/m$.
- **"Complete set of commuting observables"** is jargon until QM-5 and QM-6 need two or three quantum numbers per state. Tolerate it.
- **Symmetry implies degeneracy is one paragraph here and organizes all of atomic physics.** Do not skim Theorem 3.15: hydrogen's $m$ and $\ell$ degeneracies and the $\pm k$ of Bloch bands are all this theorem.
- **Truncated bases lie at the edge.** Example 4.2's $-59i\hbar$ is not a bug. Increase $N$ until the answer you want stops changing.

## 9. Check yourself

1. Why is $\psi(x) = \langle x|\psi\rangle$, and what is $\Phi(p)$ in the same notation?
2. Derive $\langle x|p\rangle$ up to normalization, and say where the $(2\pi\hbar)^{-1/2}$ comes from.
3. Write $\hat P = |n\rangle\langle n|$ as a matrix for $|n\rangle = (1, i)/\sqrt2$ and check $\hat P^2 = \hat P$.
4. Prove that a Hermitian operator has real eigenvalues, in one line.
5. State the Born rule and the collapse rule for a degenerate eigenvalue.
6. Show $\sum_nq_n|c_n|^2 = \langle\psi|\hat Q|\psi\rangle$.
7. Prove $[\hat A, \hat B\hat C] = [\hat A, \hat B]\hat C + \hat B[\hat A, \hat C]$, and state its Poisson-bracket twin.
8. Where in the proof of Robertson's inequality is Cauchy–Schwarz used, and what is thrown away? Why must the minimum-uncertainty wavefunction be a Gaussian?
9. Prove $e^{-i\hat Ht/\hbar}$ is unitary for Hermitian $\hat H$.
10. Write the Heisenberg equation and solve it for $\hat x_H(t)$ in the oscillator.
11. What is the precession frequency of the Bloch vector under $\hat H = \mathbf h\cdot\boldsymbol\sigma$, and what is the maximum flip probability for $\varepsilon\sigma_z + \Delta\sigma_x$?
12. Show $e^{-i\hat pa/\hbar}\psi(x) = \psi(x - a)$, and explain why $[\hat p, \hat H] = 0$ means momentum is conserved.

## 10. Pair with

- **Griffiths & Schroeter, ch. 3.** The same material more gently, with more short problems; its §3.5 on energy–time matches §3.9.
- **Shankar, ch. 1 and 4.** The best formal treatment at this level; read ch. 1 if §3.1–3.4 moved too fast.
- **Townsend, ch. 1–4.** Spin-first Dirac notation from Stern–Gerlach before any wavefunctions; read in parallel, the best preparation for QC-1.
- **MIT 8.05 (Zwiebach), lectures 1–8; *Mastering Quantum Mechanics*, ch. 13–16.** Careful on bras as functionals, both uncertainty relations, and spin precession.
- **Dirac, *Principles*, ch. 1–3.** The original, short and clear on why states form a vector space.
- **Sakurai & Napolitano, ch. 1–2 (ambitious).** Translations, pictures and propagators at graduate level.
