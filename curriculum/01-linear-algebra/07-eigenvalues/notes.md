# LA-7: Eigenvalues and Eigenvectors: Diagonalization, Powers, the Matrix Exponential and ODE Systems, Complex Eigenvalues, Markov Chains

**Time:** ~14 hours · **Prerequisites:** LA-6 (determinants), F1 (complex numbers and $e^{i\theta}$) · **Unlocks:** LA-8 (spectral theorem), LA-11 (Jordan form), CM-6 (normal modes)

> **Conventions in this module.** $A$ is an $n\times n$ matrix, real unless said otherwise. An eigenvalue is $\lambda$ (lambda), an eigenvector is $\mathbf{x}$, and an eigenvector is **never the zero vector**. Eigenvalues of a real matrix may be complex, so when needed we work in $\mathbb{C}^n$. The characteristic polynomial is $p_A(\lambda) = \det(A - \lambda I)$ (Strang's convention; Axler uses $\det(\lambda I - A)$, which differs by the factor $(-1)^n$ and has the same roots). $X$ is the matrix whose columns are eigenvectors and $\Lambda = \operatorname{diag}(\lambda_1, \dots, \lambda_n)$. A **Markov matrix** is column-stochastic: entries $\ge 0$, every column sums to $1$; $a_{ij}$ is the probability of moving *from* state $j$ *to* state $i$; probability vectors are columns and the chain is $\mathbf{u}_{k+1} = A\mathbf{u}_k$. (Probability texts, including F5, use row vectors and the transpose. Translate.) The complex inner product is $\langle \mathbf{u}, \mathbf{v}\rangle = \bar{\mathbf{u}}^\top\mathbf{v}$, conjugating the **first** slot (the physics convention of LA-8 and QM-4), so $\langle\mathbf{v},\mathbf{v}\rangle = \sum_i |v_i|^2 = \|\mathbf{v}\|^2$. $i = \sqrt{-1}$; $\omega$ is an angular frequency; for a $2\times 2$ matrix $T = \operatorname{tr}A$ and $D = \det A$. Time derivatives are $\dot{\mathbf{u}}$ or $d\mathbf{u}/dt$. The exponential of a matrix is defined by its power series (3.12), never entrywise.

## 0. Where this sits

You need: null spaces and independence (LA-3), the idea that $M^{-1}AM$ is the same linear map in a different basis (LA-4), and from LA-6 the facts that $\det A = 0$ exactly when $A$ is singular (Theorem 3.7 there), $\det(AB) = \det A\det B$ (Theorem 3.10 there), $\det A^\top = \det A$ (Theorem 3.11), the Leibniz formula (Theorem 3.9, eq. 3.9), and cofactor expansion (Theorem 3.12, eq. 3.12). From F1 you need $e^{i\theta} = \cos\theta + i\sin\theta$, conjugates, and that a degree-$n$ polynomial has $n$ complex roots counted with multiplicity.

This module is the payoff of the first six. Everything in it rests on one question: *which vectors does $A$ merely stretch?* The answer makes powers $A^k$, the differential equation $\dot{\mathbf{u}} = A\mathbf{u}$, and Markov chains into collections of one-dimensional problems. LA-8 proves that symmetric matrices always have a full set of perpendicular such directions (we prove the easy half here). LA-11 handles the matrices that don't have enough of them. CM-6 is this module's §3.10 done with Lagrangians. QM-4 is this module with complex Hermitian matrices, where eigenvalues are the possible outcomes of a measurement.

## 1. The problem that forced this idea

Nobody set out to find eigenvalues. They kept showing up as the roots of a polynomial that appeared whenever a linear system was pushed through time.

**Vibration.** In the 1740s and 1750s d'Alembert, Daniel Bernoulli, and Euler argued about the vibrating string. Lagrange, in his early work (1759–1765), attacked the discretized version: $n$ equal masses on a weightless string, each pulled by its neighbors. The equations of motion are linear, $M\ddot{\mathbf{x}} = -K\mathbf{x}$, and Lagrange found that there are exactly $n$ special "shapes" in which every mass oscillates at a single common frequency, and that every motion is a sum of these. Finding the frequencies meant finding the roots of a degree-$n$ polynomial. Those shapes are the eigenvectors of $M^{-1}K$ and the squared frequencies are its eigenvalues; §3.10 derives this in a page. The two-mass version is Example 4.3, and CM-6 is the full story.

**The stability of the solar system.** Lagrange and Laplace (1770s–1780s) asked whether the slow, cumulative ("secular") perturbations the planets exert on one another would eventually destroy the system. Linearizing the long-term drift of the eccentricities gives two vectors of variables, $\mathbf{h}$ and $\mathbf{k}$ (one entry per planet), with $\dot{\mathbf{h}} = B\mathbf{k}$ and $\dot{\mathbf{k}} = -B\mathbf{h}$, where $B$ is a constant matrix that becomes symmetric once each planet's variables are rescaled by a factor built from its mass and orbit size. As a first-order system $\dot{\mathbf{u}} = A\mathbf{u}$ for $\mathbf{u} = (\mathbf{h}, \mathbf{k})$, the matrix is $A = \begin{pmatrix}0 & B\\ -B & 0\end{pmatrix}$, which is *not* symmetric. Its eigenvalues are $\pm ig$, where $g$ runs over the roots of $\det(B - gI) = 0$, the **secular equation**. If every $g$ is real, the eigenvalues of $A$ are purely imaginary and every solution is a sum of bounded oscillations $e^{\pm igt}$ (the boundary case of Theorem 3.13 below); a non-real $g$ would give an eigenvalue of $A$ with positive real part, and exponential growth. Differentiating once gives $\ddot{\mathbf{h}} = -B^2\mathbf{h}$, a second-order system of the kind in §3.10. Laplace gave an argument specific to the planets: a sum weighted by their (positive) masses stays constant, which keeps the eccentricities bounded. Cauchy (1829), in a paper on this secular equation, proved the general algebraic fact: a real symmetric matrix has only real eigenvalues. That is Theorem 3.17 below, proved the way Cauchy's successors would, with the complex inner product. Cauchy's 1840 name for the polynomial, *équation caractéristique*, is the one that stuck in French and English.

**Rigid bodies.** Euler (1765) found that any solid body has three perpendicular axes about which it can spin without wobbling. Those are the eigenvectors of the inertia tensor (CM-7); the symmetry of the tensor is why they are perpendicular (Theorem 3.17 again).

**Iteration.** Andrey Markov (1906–1913) studied sequences of dependent random variables where the next state depends only on the current one. His 1913 example counted vowels and consonants in 20,000 letters of Pushkin's *Eugene Onegin*: a two-state chain. The long-run fraction of vowels is an eigenvector of the transition matrix for the eigenvalue $1$. Perron (1907) and Frobenius (1908–1912) proved the general facts about positive matrices that make these chains converge. Today the same eigenvector, for a matrix with billions of states, is PageRank.

**The name.** Hilbert, working on integral equations around 1904, called the special values *Eigenwerte*: "own values," the values proper to the operator. English resisted the hybrid "eigenvalue" for decades ("characteristic value," "proper value," and Sylvester's "latent root" all circulated) and then gave in.

The common thread: in a direction $A$ doesn't rotate, $A$ acts like a *number*. If you can find $n$ independent such directions, you have replaced one $n$-dimensional problem with $n$ one-dimensional ones, and one-dimensional linear problems are things you solved in your first calculus course: $u_{k+1} = \lambda u_k$ gives $\lambda^k u_0$, and $\dot u = \lambda u$ gives $e^{\lambda t}u_0$.

## 2. Intuition

### Watch a matrix act on the unit circle

Take $A = \begin{pmatrix} 2 & 1 \\ 1 & 2\end{pmatrix}$ and feed it unit vectors from every direction. Each output is rotated by some angle and stretched by some factor:

| input $\mathbf{v}$ (unit) | output $A\mathbf{v}$ | rotation | stretch |
|---|---|---|---|
| $(1, 0)$ | $(2, 1)$ | $+26.6°$ | $2.24$ |
| $(0, 1)$ | $(1, 2)$ | $-26.6°$ | $2.24$ |
| $(2, 1)/\sqrt5$ | $(5, 4)/\sqrt5$ | $+12.1°$ | $2.86$ |
| $(1, 2)/\sqrt5$ | $(4, 5)/\sqrt5$ | $-12.1°$ | $2.86$ |
| $(1, 1)/\sqrt2$ | $(3, 3)/\sqrt2$ | $0°$ | $3$ |
| $(1, -1)/\sqrt2$ | $(1, -1)/\sqrt2$ | $0°$ | $1$ |

Most directions get turned. Two do not: $(1, 1)$ is tripled and $(1, -1)$ is left alone. Those two are the **eigenvectors**, and the factors $3$ and $1$ are the **eigenvalues**. Every vector between them gets rotated *toward* $(1, 1)$, because its $(1,1)$ component is tripled while its $(1,-1)$ component isn't.

```
                 x2
                  |        /  (1,1): A stretches by 3, no rotation
                  |      /
                  |    /
                  |  /
   ---------------+---------------- x1
                  |  \
                  |    \
                  |      \
                  |        \  (1,-1): A stretches by 1, no rotation
```

The project draws the full picture, a fan of unit vectors and where they land, when you run it (Part 0 saves `unit_vectors_and_images.png` into `figures/`).

### Why these directions make dynamics one-dimensional

Now iterate: $\mathbf{u}_{k+1} = A\mathbf{u}_k$. Pick any start, say $\mathbf{u}_0 = (1, 0)$. Write it in the eigenvector basis: $(1, 0) = \tfrac12(1, 1) + \tfrac12(1, -1)$. Apply $A$ once: the first piece triples, the second stays. Apply it $k$ times:

$$
\mathbf{u}_k = A^k\mathbf{u}_0 = \tfrac12\, 3^k (1, 1) + \tfrac12\, 1^k (1, -1) = \tfrac12\begin{pmatrix} 3^k + 1 \\ 3^k - 1\end{pmatrix}.
$$

You never multiplied a matrix. Each component along an eigenvector evolves on its own, as a number times a power. After a few steps $\mathbf{u}_k$ is essentially along $(1, 1)$: the eigenvector with the largest eigenvalue takes over.

The same thing happens for a differential equation $\dot{\mathbf{u}} = A\mathbf{u}$. Along $(1, 1)$ the equation reads $\dot c = 3c$, so that component grows like $e^{3t}$; along $(1, -1)$ it reads $\dot c = c$, so $e^{t}$. The solution is $c_1 e^{3t}(1,1) + c_2 e^{t}(1, -1)$. Compare: $\lambda^k$ for iteration, $e^{\lambda t}$ for flow. Same eigenvectors, same eigenvalues, two different functions of them.

### When it fails, and what the failure looks like

Not every matrix cooperates. A rotation by $90°$ turns *every* real vector, so there is no real eigenvector at all. A shear $\begin{pmatrix}1 & 1\\ 0 & 1\end{pmatrix}$ fixes the $x$-axis and tilts everything else, so there is only *one* eigenvector direction, and one direction isn't enough to decompose $\mathbb{R}^2$. The first failure is cured by allowing complex numbers (§3.11): the rotation has eigenvalues $e^{\pm i\theta}$, and complex eigenvalues are what oscillation looks like in this language. The second failure is real and permanent (§3.14); its cure is LA-11.

### What you should now be able to guess

- Eigenvectors for *different* eigenvalues point in genuinely different directions (independent), so $n$ distinct eigenvalues give a basis.
- In that basis $A$ is diagonal: $A = X\Lambda X^{-1}$.
- $A^k = X\Lambda^k X^{-1}$ and $e^{At} = X e^{\Lambda t} X^{-1}$: any function of $A$ is that function applied to the eigenvalues, in the eigenvector basis.
- Long-run behavior is set by the eigenvalue of largest modulus (iteration) or largest real part (flow).
- Symmetric matrices, like the example above, have perpendicular eigenvectors.

Section 3 proves each of these.

## 3. Formal development

### 3.1 Definition and eigenspaces

**Definition 3.1 (Eigenvalue, eigenvector).** Let $A$ be $n \times n$. A scalar $\lambda \in \mathbb{C}$ is an **eigenvalue** of $A$ if there is a **nonzero** vector $\mathbf{x} \in \mathbb{C}^n$ with
$$
A\mathbf{x} = \lambda\mathbf{x}. \tag{3.1}
$$
Any such $\mathbf{x}$ is an **eigenvector** for $\lambda$.

Two remarks on the definition. First, $\mathbf{x} = \mathbf{0}$ is excluded because $A\mathbf{0} = \lambda\mathbf{0}$ for every $\lambda$; allowing it would make every number an eigenvalue. Second, $\lambda = 0$ is allowed: $A\mathbf{x} = 0\cdot\mathbf{x}$ says $\mathbf{x}$ is in the null space, so **$0$ is an eigenvalue exactly when $A$ is singular**.

**Theorem 3.2 (Eigenspaces are subspaces).** For a fixed $\lambda$, the set $E_\lambda = \{\mathbf{x} : A\mathbf{x} = \lambda\mathbf{x}\}$ is a subspace, namely $E_\lambda = N(A - \lambda I)$. Its nonzero members are the eigenvectors for $\lambda$.

*Proof.* $A\mathbf{x} = \lambda\mathbf{x} \iff A\mathbf{x} - \lambda I\mathbf{x} = \mathbf{0} \iff (A - \lambda I)\mathbf{x} = \mathbf{0}$. So $E_\lambda$ is the null space of $A - \lambda I$, and null spaces are subspaces (LA-3). ∎

Consequence: eigenvectors are never unique. If $\mathbf{x}$ works, so does $c\mathbf{x}$ for every $c \ne 0$ (and, for complex $\mathbf{x}$, every complex $c$). "The eigenvector" means "any nonzero vector in $E_\lambda$." Keep this in mind when you compare with NumPy in §8.

### 3.2 Finding eigenvalues: the characteristic equation

▶ **Try it first.** $\lambda$ is an eigenvalue exactly when $A - \lambda I$ has a nonzero vector in its null space. Turn that into a statement about a determinant using LA-6.

**Derivation.**

1. $\lambda$ is an eigenvalue $\iff$ there is $\mathbf{x} \ne \mathbf{0}$ with $(A - \lambda I)\mathbf{x} = \mathbf{0}$ $\iff$ $N(A - \lambda I) \ne \{\mathbf{0}\}$.
2. A square matrix has a nontrivial null space exactly when it is singular (LA-2: fewer than $n$ pivots; LA-3: columns dependent).
3. A square matrix is singular exactly when its determinant is zero (LA-6, Theorem 3.7). So
$$
\lambda \text{ is an eigenvalue of } A \iff \det(A - \lambda I) = 0. \tag{3.2}
$$
∎

■ **What just happened.** We traded a question about vectors ("is there an $\mathbf{x}$?") for a question about one number ("is this determinant zero?"). The vectors come back in a second step: once you know $\lambda$, solve $(A - \lambda I)\mathbf{x} = \mathbf{0}$ by elimination (LA-2) to get a basis of $E_\lambda$.

**Definition 3.3 (Characteristic polynomial).** $p_A(\lambda) = \det(A - \lambda I)$.

*Why it is a polynomial of degree $n$, and what its top coefficient is.* By the Leibniz formula (LA-6, Theorem 3.9, eq. 3.9), $\det(A - \lambda I)$ is a sum over permutations $\sigma$ of $\pm$ products of $n$ entries of $A - \lambda I$, one from each row and column. The off-diagonal entries are the constants $a_{ij}$; the diagonal entries are $a_{ii} - \lambda$. The identity permutation contributes $\prod_i (a_{ii} - \lambda)$, which has degree $n$ with leading term $(-\lambda)^n = (-1)^n\lambda^n$. Any other permutation moves at least two indices, so it picks at most $n - 2$ diagonal entries and has degree $\le n - 2$. Hence $p_A$ has degree exactly $n$ and leading coefficient $(-1)^n$; and **every term of degree $n-1$ comes from the identity permutation alone.** We use this in Theorem 3.4.

By the fundamental theorem of algebra (F1), $p_A$ has $n$ complex roots counted with multiplicity: every $n\times n$ matrix has $n$ eigenvalues in $\mathbb{C}$, some possibly repeated, some possibly non-real even when $A$ is real.

**The $2\times 2$ case, fully.** For $A = \begin{pmatrix} a & b \\ c & d\end{pmatrix}$,
$$
p_A(\lambda) = \det\begin{pmatrix} a - \lambda & b \\ c & d - \lambda\end{pmatrix} = (a-\lambda)(d-\lambda) - bc = \lambda^2 - (a + d)\lambda + (ad - bc) = \lambda^2 - T\lambda + D, \tag{3.3}
$$
with $T = \operatorname{tr}A$ and $D = \det A$. So
$$
\lambda_{1,2} = \frac{T \pm \sqrt{T^2 - 4D}}{2}. \tag{3.4}
$$
Real and distinct when $T^2 > 4D$, repeated when $T^2 = 4D$, a complex-conjugate pair when $T^2 < 4D$. You should be able to write (3.3) without thinking by the end of this module.

**The $3\times 3$ case.** Expand $\det(A - \lambda I)$ by cofactors along a row (LA-6, Theorem 3.12). You get a cubic. Find one root by testing the integer divisors of the constant term $\det A$ (if the matrix was designed by a human, there is one); divide it out; solve the remaining quadratic. Example 4.6 does this in full. Beyond $3\times 3$ this procedure is not how anyone computes eigenvalues; see §8.

### 3.3 Trace and determinant

**Theorem 3.4 (Trace and determinant from eigenvalues).** Let $\lambda_1, \dots, \lambda_n$ be the eigenvalues of $A$, counted with multiplicity, complex ones included. Then
$$
\lambda_1 + \cdots + \lambda_n = \operatorname{tr}A = \sum_i a_{ii}, \qquad \lambda_1 \lambda_2 \cdots \lambda_n = \det A.
$$

▶ **Try it first.** Write $p_A(\lambda)$ two ways: factored over its roots, and expanded from the determinant. Compare the coefficients of $\lambda^{n-1}$ and of $\lambda^0$.

**Derivation.**

1. *Root form.* $p_A$ has degree $n$, leading coefficient $(-1)^n$, and roots $\lambda_1, \dots, \lambda_n$, so $p_A(\lambda) = (-1)^n\prod_{i}(\lambda - \lambda_i) = \prod_i (\lambda_i - \lambda)$. Expand the product: taking $-\lambda$ from every factor gives $(-1)^n\lambda^n$; taking $-\lambda$ from all but one factor and $\lambda_i$ from that one gives $(-1)^{n-1}\lambda^{n-1}\sum_i\lambda_i$; setting $\lambda = 0$ gives the constant term $\prod_i\lambda_i$. So
$$
p_A(\lambda) = (-1)^n\lambda^n + (-1)^{n-1}\Big(\sum_i\lambda_i\Big)\lambda^{n-1} + \cdots + \prod_i\lambda_i.
$$
2. *Determinant form, coefficient of $\lambda^{n-1}$.* By the remark after Definition 3.3, only the identity permutation contributes degree $n-1$, via $\prod_i(a_{ii} - \lambda)$. Expanding that product exactly as in step 1 (with $a_{ii}$ in place of $\lambda_i$) gives a $\lambda^{n-1}$ coefficient of $(-1)^{n-1}\sum_i a_{ii} = (-1)^{n-1}\operatorname{tr}A$.
3. *Determinant form, constant term.* $p_A(0) = \det(A - 0\cdot I) = \det A$.
4. Two polynomials that are equal have equal coefficients. Comparing $\lambda^{n-1}$: $\sum\lambda_i = \operatorname{tr}A$. Comparing $\lambda^0$: $\prod\lambda_i = \det A$. ∎

■ **What just happened.** Two free checks on every eigenvalue computation you will ever do. For $2\times 2$ they are more than checks; (3.3) says they *are* the computation. For the $3\times 3$ in Example 4.6, "do the eigenvalues add to the trace?" catches most arithmetic slips.

### 3.4 Eigenvalues you can read off

**Theorem 3.5.**
(a) If $A$ is triangular (upper, lower, or diagonal), its eigenvalues are its diagonal entries.
(b) If $A\mathbf{x} = \lambda\mathbf{x}$, then for every integer $k \ge 1$, $A^k\mathbf{x} = \lambda^k\mathbf{x}$; if $A$ is invertible, $\lambda \ne 0$ and $A^{-1}\mathbf{x} = \lambda^{-1}\mathbf{x}$; and for any scalar $c$, $(A + cI)\mathbf{x} = (\lambda + c)\mathbf{x}$. Same eigenvector every time.
(c) $A$ and $A^\top$ have the same eigenvalues (but, in general, different eigenvectors).

*Proof.* (a) $A - \lambda I$ is triangular with diagonal entries $a_{ii} - \lambda$, so $\det(A - \lambda I) = \prod_i (a_{ii} - \lambda)$ (LA-6, Theorem 3.6), whose roots are the $a_{ii}$.

(b) Induction: $A^k\mathbf{x} = A(A^{k-1}\mathbf{x}) = A(\lambda^{k-1}\mathbf{x}) = \lambda^{k-1}A\mathbf{x} = \lambda^k\mathbf{x}$. Inverse: apply $A^{-1}$ to both sides of $A\mathbf{x} = \lambda\mathbf{x}$ to get $\mathbf{x} = \lambda A^{-1}\mathbf{x}$; if $\lambda$ were $0$ this would say $\mathbf{x} = \mathbf{0}$, impossible, so divide by $\lambda$. Shift: $(A + cI)\mathbf{x} = A\mathbf{x} + c\mathbf{x} = (\lambda + c)\mathbf{x}$.

(c) $\det(A^\top - \lambda I) = \det\big((A - \lambda I)^\top\big) = \det(A - \lambda I)$ by LA-6, Theorem 3.11. Same polynomial, same roots. ∎

**Warning.** The eigenvalues of $A + B$ are **not** $\lambda_i + \mu_j$ and the eigenvalues of $AB$ are **not** $\lambda_i\mu_j$ in general. Take $A = \begin{pmatrix}0&1\\0&0\end{pmatrix}$ and $B = \begin{pmatrix}0&0\\1&0\end{pmatrix}$. Both have eigenvalues $0, 0$. But $A + B = \begin{pmatrix}0&1\\1&0\end{pmatrix}$ has eigenvalues $\pm 1$, and $AB = \begin{pmatrix}1&0\\0&0\end{pmatrix}$ has eigenvalues $1, 0$. Theorem 3.5(b) works because $A^k$, $A^{-1}$, and $A + cI$ all *share eigenvectors* with $A$. Exercise B1(c) shows that a shared eigenbasis is sufficient. It is not necessary: $\begin{pmatrix}1&1\\0&2\end{pmatrix}$ and $\operatorname{diag}(2, 1)$ do not commute, so by B1(c) they have no common eigenbasis, yet their sum $\begin{pmatrix}3&1\\0&3\end{pmatrix}$ has eigenvalues $1 + 2$ and $2 + 1$ and their product $\begin{pmatrix}2&1\\0&2\end{pmatrix}$ has $1\cdot2$ and $2\cdot1$. Both are upper triangular, and for triangular matrices the diagonals add and multiply entry by entry.

### 3.5 Independence, and diagonalization

**Theorem 3.6 (Eigenvectors for distinct eigenvalues are independent).** If $\mathbf{x}_1, \dots, \mathbf{x}_k$ are eigenvectors of $A$ for eigenvalues $\lambda_1, \dots, \lambda_k$ that are **all different**, then $\mathbf{x}_1, \dots, \mathbf{x}_k$ are linearly independent.

▶ **Try it first.** Do $k = 2$. Suppose $c_1\mathbf{x}_1 + c_2\mathbf{x}_2 = \mathbf{0}$. Apply $A$ to this equation. Separately, multiply it by $\lambda_2$. Subtract.

**Proof (induction on $k$).**

1. $k = 1$: a single nonzero vector is independent.
2. Suppose the statement holds for $k - 1$ vectors, and suppose
$$
c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_k\mathbf{x}_k = \mathbf{0}. \tag{★}
$$
3. Apply $A$ to (★), using $A\mathbf{x}_i = \lambda_i\mathbf{x}_i$: $\;c_1\lambda_1\mathbf{x}_1 + \cdots + c_k\lambda_k\mathbf{x}_k = \mathbf{0}$.
4. Multiply (★) by the scalar $\lambda_k$: $\;c_1\lambda_k\mathbf{x}_1 + \cdots + c_k\lambda_k\mathbf{x}_k = \mathbf{0}$.
5. Subtract step 4 from step 3. The $\mathbf{x}_k$ terms cancel:
$$
c_1(\lambda_1 - \lambda_k)\mathbf{x}_1 + \cdots + c_{k-1}(\lambda_{k-1} - \lambda_k)\mathbf{x}_{k-1} = \mathbf{0}.
$$
6. By the induction hypothesis $\mathbf{x}_1, \dots, \mathbf{x}_{k-1}$ are independent, so every coefficient is zero: $c_i(\lambda_i - \lambda_k) = 0$ for $i < k$. Since $\lambda_i \ne \lambda_k$, each $c_i = 0$ for $i < k$.
7. Now (★) reads $c_k\mathbf{x}_k = \mathbf{0}$ with $\mathbf{x}_k \ne \mathbf{0}$, so $c_k = 0$. All coefficients vanish. ∎

**Corollary 3.7.** If an $n\times n$ matrix has $n$ distinct eigenvalues, it has $n$ independent eigenvectors, which form a basis of $\mathbb{C}^n$ (of $\mathbb{R}^n$ if the eigenvalues are real).

Now suppose $A$ has $n$ independent eigenvectors $\mathbf{x}_1, \dots, \mathbf{x}_n$ (distinct eigenvalues or not) with eigenvalues $\lambda_1, \dots, \lambda_n$. Put the eigenvectors in the columns of $X$ and the eigenvalues on the diagonal of $\Lambda$.

▶ **Try it first.** Compute $AX$ column by column (LA-1, eq. 3.3). Then compute $X\Lambda$ column by column. Compare.

**Derivation (diagonalization).**

1. The $j$-th column of $AX$ is $A\mathbf{x}_j = \lambda_j\mathbf{x}_j$.
2. The $j$-th column of $X\Lambda$ is $X(\Lambda\mathbf{e}_j) = X(\lambda_j\mathbf{e}_j) = \lambda_j X\mathbf{e}_j = \lambda_j\mathbf{x}_j$ (LA-1, Theorem 3.5: $X\mathbf{e}_j$ is column $j$).
3. Every column agrees, so
$$
AX = X\Lambda. \tag{3.5}
$$
4. $X$ is square with independent columns, hence invertible (LA-2/LA-3). Multiply (3.5) on the right by $X^{-1}$, or on the left:
$$
A = X\Lambda X^{-1}, \qquad X^{-1}AX = \Lambda. \tag{3.6}
$$
∎

**Theorem 3.8 (Diagonalizability).** $A$ can be written $A = X\Lambda X^{-1}$ with $\Lambda$ diagonal and $X$ invertible **if and only if** $A$ has $n$ linearly independent eigenvectors. In that case the columns of $X$ are those eigenvectors and the diagonal of $\Lambda$ holds the matching eigenvalues.

*Proof.* "If" is the derivation above. "Only if": from $A = X\Lambda X^{-1}$ we get $AX = X\Lambda$, and reading that column by column (steps 1–2 backwards) says $A\mathbf{x}_j = \lambda_j\mathbf{x}_j$; the $\mathbf{x}_j$ are independent because $X$ is invertible. ∎

■ **What just happened.** $X^{-1}AX = \Lambda$ is the change-of-basis formula of LA-4: $\Lambda$ is the matrix of the *same linear map* written in the eigenvector basis, and in that basis the map scales each basis vector, so its matrix is diagonal. Diagonalization is not a trick; it is choosing the right coordinates.

Two things that are **not** the same: *invertible* ("is $0$ an eigenvalue?") and *diagonalizable* ("are there enough eigenvectors?"). The four combinations all occur. $I$ is invertible and diagonalizable. The zero matrix is diagonalizable (it *is* diagonal) and singular. The shear $\begin{pmatrix}1&1\\0&1\end{pmatrix}$ is invertible and not diagonalizable (§3.14). $N = \begin{pmatrix}0&1\\0&0\end{pmatrix}$ is singular and not diagonalizable (eigenvalue $0$ twice, only the eigenvector direction $\mathbf{e}_1$).

### 3.6 Powers, dynamics, and power iteration

**Theorem 3.9 (Powers).** If $A = X\Lambda X^{-1}$ then for every $k \ge 0$
$$
A^k = X\Lambda^k X^{-1}, \qquad \Lambda^k = \operatorname{diag}(\lambda_1^k, \dots, \lambda_n^k). \tag{3.7}
$$
Consequently, if $\mathbf{u}_0 = X\mathbf{c}$ (that is, $\mathbf{c} = X^{-1}\mathbf{u}_0$ holds the coordinates of $\mathbf{u}_0$ in the eigenvector basis), then
$$
\mathbf{u}_k = A^k\mathbf{u}_0 = c_1\lambda_1^k\mathbf{x}_1 + c_2\lambda_2^k\mathbf{x}_2 + \cdots + c_n\lambda_n^k\mathbf{x}_n. \tag{3.8}
$$

*Proof.* $A^2 = X\Lambda X^{-1}X\Lambda X^{-1} = X\Lambda^2X^{-1}$ because the inner $X^{-1}X$ is $I$; induction gives $A^k$. A power of a diagonal matrix is diagonal with powered entries. Then $A^k\mathbf{u}_0 = X\Lambda^kX^{-1}\mathbf{u}_0 = X\Lambda^k\mathbf{c} = X(\lambda_1^kc_1, \dots, \lambda_n^kc_n)^\top = \sum_i c_i\lambda_i^k\mathbf{x}_i$. ∎

Equation (3.8) is the whole theory of linear recurrences. Each eigen-component grows, decays, or oscillates according to $|\lambda_i|$: $A^k\mathbf{u}_0 \to \mathbf{0}$ for every $\mathbf{u}_0$ exactly when all $|\lambda_i| < 1$; it stays bounded when all $|\lambda_i| \le 1$ (for diagonalizable $A$); it blows up along $\mathbf{x}_i$ when $|\lambda_i| > 1$ and $c_i \ne 0$.

**Power iteration.** Suppose one eigenvalue is strictly largest in modulus, $|\lambda_1| > |\lambda_2| \ge |\lambda_3| \ge \cdots$, and $c_1 \ne 0$.

▶ **Try it first.** Factor $\lambda_1^k$ out of (3.8). What does the bracket converge to, and how fast?

**Derivation.**

1. $A^k\mathbf{u}_0 = \lambda_1^k\Big[c_1\mathbf{x}_1 + \sum_{i\ge 2}c_i\big(\tfrac{\lambda_i}{\lambda_1}\big)^k\mathbf{x}_i\Big]$.
2. Each ratio satisfies $|\lambda_i/\lambda_1| \le |\lambda_2/\lambda_1| < 1$, so the bracket tends to $c_1\mathbf{x}_1$ with error
$$
\Big\|\sum_{i\ge2}c_i\big(\tfrac{\lambda_i}{\lambda_1}\big)^k\mathbf{x}_i\Big\| \le \Big(\sum_{i\ge2}|c_i|\,\|\mathbf{x}_i\|\Big)\Big|\frac{\lambda_2}{\lambda_1}\Big|^k = C\,\Big|\frac{\lambda_2}{\lambda_1}\Big|^k. \tag{3.9}
$$
3. In practice we cannot let $\lambda_1^k$ grow or shrink the vector, so after each multiplication we divide by the norm: $\mathbf{v}_{k+1} = A\mathbf{v}_k/\|A\mathbf{v}_k\|$. Normalizing doesn't change direction, so $\mathbf{v}_k \to \pm\mathbf{x}_1/\|\mathbf{x}_1\|$ (the sign alternates if $\lambda_1 < 0$), with the direction error shrinking like $|\lambda_2/\lambda_1|^k$.
4. Estimate the eigenvalue by the **Rayleigh quotient** $\lambda^{(k)} = \mathbf{v}_k^\top A\mathbf{v}_k / \mathbf{v}_k^\top\mathbf{v}_k$: if $\mathbf{v}_k$ were exactly $\mathbf{x}_1$ this would be exactly $\lambda_1$, and in general the error is also $O(|\lambda_2/\lambda_1|^k)$. For symmetric $A$ the error is $O(|\lambda_2/\lambda_1|^{2k})$: Exercise B8 measures both rates and shows why they differ. ∎

■ **What just happened.** Convergence is *linear* with ratio $|\lambda_2/\lambda_1|$: each step buys the same number of digits, and if the top two eigenvalues are close, that number is small. The fix is **inverse iteration with a shift**: by Theorem 3.5(b), $(A - \sigma I)^{-1}$ has eigenvalues $1/(\lambda_i - \sigma)$ with the same eigenvectors, so power iteration on it converges to the eigenvector whose $\lambda_i$ is *closest to $\sigma$*, at ratio $|\lambda_{\text{closest}} - \sigma|/|\lambda_{\text{next}} - \sigma|$, which is tiny when $\sigma$ is a good guess. Both are in the project, and B8(c) measures the second rate. LA-13 §3.8 builds the QR algorithm, which is what `np.linalg.eig` actually runs, on top of these two ideas (LA-13 Theorem 3.16: the QR algorithm is simultaneous power iteration).

### 3.7 Markov matrices

**Definition 3.10.** $A$ is a **Markov matrix** (column-stochastic) if every $a_{ij} \ge 0$ and every column sums to $1$: $\sum_i a_{ij} = 1$ for each $j$. Interpretation: $a_{ij}$ is the probability of being in state $i$ next step given state $j$ now. A **probability vector** has nonnegative entries summing to $1$; $u_i$ is the probability of being in state $i$. The chain is $\mathbf{u}_{k+1} = A\mathbf{u}_k$.

Writing $\mathbf{1} = (1, \dots, 1)^\top$, "every column sums to $1$" is the single equation $\mathbf{1}^\top A = \mathbf{1}^\top$. From it, $A$ sends probability vectors to probability vectors: entries of $A\mathbf{u}$ are sums of products of nonnegative numbers, and $\mathbf{1}^\top(A\mathbf{u}) = (\mathbf{1}^\top A)\mathbf{u} = \mathbf{1}^\top\mathbf{u} = 1$.

**Theorem 3.11.** For every Markov matrix $A$: (a) $\lambda = 1$ is an eigenvalue; (b) every eigenvalue satisfies $|\lambda| \le 1$.

▶ **Try it first (a).** $\mathbf{1}^\top A = \mathbf{1}^\top$ is an eigenvector statement, but about rows. Transpose it, and use Theorem 3.5(c).

**Proof of (a).**

1. Transposing $\mathbf{1}^\top A = \mathbf{1}^\top$ gives $A^\top\mathbf{1} = \mathbf{1}$: the vector $\mathbf{1}$ is an eigenvector of $A^\top$ with eigenvalue $1$. (One says $\mathbf{1}^\top$ is a *left eigenvector* of $A$.)
2. $A$ and $A^\top$ have the same eigenvalues (Theorem 3.5(c)).
3. Therefore $1$ is an eigenvalue of $A$. ∎

Its eigenvector is **not** $\mathbf{1}$ in general; you must solve $(A - I)\mathbf{x} = \mathbf{0}$ (Example 4.1).

▶ **Try it first (b).** Let $A^\top\mathbf{v} = \lambda\mathbf{v}$. Look at the entry of $\mathbf{v}$ with the largest modulus and apply the triangle inequality to that row of the equation.

**Proof of (b).**

1. Let $A^\top\mathbf{v} = \lambda\mathbf{v}$ with $\mathbf{v} \ne \mathbf{0}$, possibly complex (again, same eigenvalues as $A$). Choose $m$ with $|v_m| = \max_j|v_j| > 0$.
2. Row $m$ of $A^\top\mathbf{v} = \lambda\mathbf{v}$ reads $\lambda v_m = \sum_j (A^\top)_{mj}v_j = \sum_j a_{jm}v_j$. (The row of $A^\top$ is the column of $A$.)
3. Triangle inequality, then $|v_j| \le |v_m|$, then column $m$ sums to $1$:
$$
|\lambda|\,|v_m| = \Big|\sum_j a_{jm}v_j\Big| \le \sum_j a_{jm}|v_j| \le \sum_j a_{jm}|v_m| = |v_m|.
$$
4. Divide by $|v_m| > 0$: $|\lambda| \le 1$. ∎

■ **What just happened.** The proof used $A^\top$ because its *rows* sum to $1$, which is what you need when you take the modulus of a weighted average: a weighted average of numbers of modulus $\le |v_m|$ has modulus $\le |v_m|$. Nothing is lost by passing to the transpose because eigenvalues don't care.

**Steady state.** A probability vector $\mathbf{u}_\infty$ with $A\mathbf{u}_\infty = \mathbf{u}_\infty$ is a **steady state**: an eigenvector for $\lambda = 1$ scaled to sum to $1$. If $A$ is diagonalizable with $\lambda_1 = 1$ and $|\lambda_i| < 1$ for $i \ge 2$, then (3.8) gives $\mathbf{u}_k = c_1\mathbf{x}_1 + \sum_{i\ge2}c_i\lambda_i^k\mathbf{x}_i \to c_1\mathbf{x}_1$, and since every $\mathbf{u}_k$ sums to $1$, so does the limit: the chain converges to the steady state from **every** starting distribution, and the convergence rate is $|\lambda_2|^k$. Exercise B4 shows that the transient eigenvectors each sum to zero and that, when $\lambda = 1$ is a simple eigenvalue (as here, since $|\lambda_i| < 1$ for $i \ge 2$), $c_1$ does not depend on $\mathbf{u}_0$.

When are the other eigenvalues strictly inside the unit circle? Not always: $A = \begin{pmatrix}0&1\\1&0\end{pmatrix}$ is Markov with eigenvalues $1$ and $-1$; it swaps the two states forever and never settles. Positivity rules this out:

**Theorem 3.12 (Perron–Frobenius for Markov matrices; stated).** If $A$ is a Markov matrix with **all entries strictly positive** (or if some power $A^m$ has all entries strictly positive), then: $\lambda = 1$ has a one-dimensional eigenspace; its eigenvector can be chosen with all entries positive, giving a unique steady state $\mathbf{u}_\infty$; every other eigenvalue has $|\lambda| < 1$; and $A^k\mathbf{u}_0 \to \mathbf{u}_\infty$ for every probability vector $\mathbf{u}_0$.

We do not prove all of this here. Exercise C2 proves the two parts that follow from the argument above ($|\lambda| = 1 \Rightarrow \lambda = 1$, and the eigenspace is one-dimensional) and the convergence for diagonalizable $A$. The positivity of the eigenvector and the non-diagonalizable case are in Strang, *Introduction to Linear Algebra*, §10.3, and in full generality in Horn & Johnson, *Matrix Analysis*, ch. 8. F5 picks the probabilistic story up from here; ST-6 uses it to justify MCMC.

### 3.8 Differential equations: $\dot{\mathbf{u}} = A\mathbf{u}$

The system is
$$
\frac{d\mathbf{u}}{dt} = A\mathbf{u}, \qquad \mathbf{u}(0) = \mathbf{u}_0, \qquad \mathbf{u}(t) \in \mathbb{R}^n, \tag{3.10}
$$
with $A$ constant. For $n = 1$ this is $\dot u = au$ and the solution is $u_0e^{at}$.

▶ **Try it first.** Guess $\mathbf{u}(t) = e^{\lambda t}\mathbf{x}$ with a constant vector $\mathbf{x}$ and a constant $\lambda$. Substitute into (3.10). What must $\lambda$ and $\mathbf{x}$ be?

**Derivation.**

1. $\dfrac{d}{dt}\big(e^{\lambda t}\mathbf{x}\big) = \lambda e^{\lambda t}\mathbf{x}$, while $A\big(e^{\lambda t}\mathbf{x}\big) = e^{\lambda t}A\mathbf{x}$. These agree for all $t$ exactly when $A\mathbf{x} = \lambda\mathbf{x}$. So every eigenpair gives a solution $e^{\lambda t}\mathbf{x}$.
2. *Superposition.* If $\mathbf{u}_1(t)$ and $\mathbf{u}_2(t)$ solve $\dot{\mathbf{u}} = A\mathbf{u}$, so does $c_1\mathbf{u}_1 + c_2\mathbf{u}_2$: $\;(c_1\mathbf{u}_1 + c_2\mathbf{u}_2)' = c_1A\mathbf{u}_1 + c_2A\mathbf{u}_2 = A(c_1\mathbf{u}_1 + c_2\mathbf{u}_2)$. So for any constants,
$$
\mathbf{u}(t) = c_1e^{\lambda_1t}\mathbf{x}_1 + c_2e^{\lambda_2t}\mathbf{x}_2 + \cdots + c_ne^{\lambda_nt}\mathbf{x}_n \tag{3.11}
$$
solves the differential equation.
3. *Initial condition.* At $t = 0$, (3.11) reads $\mathbf{u}_0 = \sum_i c_i\mathbf{x}_i = X\mathbf{c}$. If $A$ has $n$ independent eigenvectors, $X$ is invertible and $\mathbf{c} = X^{-1}\mathbf{u}_0$. The same coordinates as in (3.8).
4. *Uniqueness.* The right side of (3.10) is linear in $\mathbf{u}$, hence Lipschitz, so the existence-and-uniqueness theorem from your ODE course (Picard–Lindelöf) says there is exactly one solution through $\mathbf{u}_0$. We have exhibited one; it is the one. ∎

■ **What just happened.** Compare (3.11) with (3.8): $\lambda^k$ became $e^{\lambda t}$ and everything else is identical. The eigenvectors are the same; only the one-dimensional problem along each changed from "multiply by $\lambda$ each step" to "grow at rate $\lambda$."

**Theorem 3.13 (Stability).** Let $A$ be diagonalizable. Every solution of (3.10) tends to $\mathbf{0}$ as $t\to\infty$ if and only if $\operatorname{Re}\lambda_i < 0$ for every eigenvalue. Every solution stays bounded if and only if $\operatorname{Re}\lambda_i \le 0$ for every $i$. If some $\operatorname{Re}\lambda_i > 0$, some solutions grow without bound.

*Proof.* $|e^{\lambda t}| = |e^{(\operatorname{Re}\lambda)t}e^{i(\operatorname{Im}\lambda)t}| = e^{(\operatorname{Re}\lambda)t}$ (F1). If every $\operatorname{Re}\lambda_i < 0$, each term of (3.11) has size $|c_i|e^{(\operatorname{Re}\lambda_i)t}\|\mathbf{x}_i\| \to 0$, so every solution decays; if every $\operatorname{Re}\lambda_i \le 0$, each term is bounded by $|c_i|\|\mathbf{x}_i\|$. Conversely, suppose some $\operatorname{Re}\lambda_j \ge 0$. We need a *real* solution that does not decay. If $\lambda_j$ is real, take $\mathbf{x}_j$ real (a null vector of the real matrix $A - \lambda_jI$); the solution starting at $\mathbf{u}_0 = \mathbf{x}_j$ is $e^{\lambda_jt}\mathbf{x}_j$ (all other $c_i = 0$ by independence), of size $e^{\lambda_jt}\|\mathbf{x}_j\|$. If $\lambda_j = a + ib$ with $b \ne 0$, then $\mathbf{x}_j = \mathbf{r} + i\mathbf{s}$ is complex and is not an admissible starting point in $\mathbb{R}^n$. Use instead the real solution $\mathbf{p}(t) = \operatorname{Re}\big(e^{\lambda_jt}\mathbf{x}_j\big) = e^{at}(\mathbf{r}\cos bt - \mathbf{s}\sin bt)$ from (3.19) in §3.11 below, which starts at $\mathbf{r}$. The vectors $\mathbf{r}, \mathbf{s}$ are independent (§3.11, step 4), so $\|\mathbf{r}\cos\theta - \mathbf{s}\sin\theta\|$ is never zero; it is a continuous function of $\theta$ with period $2\pi$, so it has a minimum $c > 0$ on $[0, 2\pi]$, and $\|\mathbf{p}(t)\| \ge c\,e^{at}$. In both cases the solution does not decay, and it grows without bound if $\operatorname{Re}\lambda_j > 0$. ∎

The diagonalizable hypothesis matters only on the boundary $\operatorname{Re}\lambda = 0$: a defective eigenvalue there produces $t\,e^{\lambda t}$ growth (§3.14).

### 3.9 The matrix exponential

**Definition 3.14.** For any square matrix $A$,
$$
e^{A} = \sum_{k=0}^\infty \frac{A^k}{k!} = I + A + \frac{A^2}{2} + \frac{A^3}{6} + \cdots \tag{3.12}
$$

*Why the series converges.* Use the Frobenius norm $\|A\|_F = \sqrt{\sum_{ij}a_{ij}^2}$ from LA-1, C3. That exercise proved $\|A\mathbf{x}\| \le \|A\|_F\|\mathbf{x}\|$; applying it to each column of $B$ gives $\|AB\|_F^2 = \sum_j\|A\mathbf{b}_j\|^2 \le \|A\|_F^2\sum_j\|\mathbf{b}_j\|^2 = \|A\|_F^2\|B\|_F^2$, so $\|A^k\|_F \le \|A\|_F^k$. Every entry of $A^k$ is bounded by $\|A^k\|_F \le \|A\|_F^k$, so each entry of the series is dominated by $\sum_k\|A\|_F^k/k! = e^{\|A\|_F} < \infty$ and converges absolutely. Absolute convergence is what licenses the rearrangements below.

**Theorem 3.15 (Properties of $e^{At}$).**
(a) If $A = X\Lambda X^{-1}$ then $e^{At} = Xe^{\Lambda t}X^{-1}$ with $e^{\Lambda t} = \operatorname{diag}(e^{\lambda_1t}, \dots, e^{\lambda_nt})$.
(b) $\dfrac{d}{dt}e^{At} = Ae^{At} = e^{At}A$ for **every** square $A$. Hence $\mathbf{u}(t) = e^{At}\mathbf{u}_0$ solves (3.10), diagonalizable or not.
(c) If $AB = BA$ then $e^{A+B} = e^Ae^B$. Without commutativity this fails.
(d) $e^{At}$ is invertible with inverse $e^{-At}$, and $e^{A(s+t)} = e^{As}e^{At}$.

▶ **Try it first (a).** Substitute $(At)^k = X(\Lambda t)^kX^{-1}$ (Theorem 3.9) into the series.

**Derivation of (a).**

1. $(At)^k = t^kA^k = t^kX\Lambda^kX^{-1} = X(\Lambda t)^kX^{-1}$.
2. For each partial sum, $\sum_{k=0}^N X(\Lambda t)^kX^{-1}/k! = X\Big(\sum_{k=0}^N(\Lambda t)^k/k!\Big)X^{-1}$; let $N\to\infty$.
3. $(\Lambda t)^k$ is diagonal with entries $(\lambda_it)^k$, so $\sum_k(\Lambda t)^k/k!$ is diagonal with entries $\sum_k(\lambda_it)^k/k! = e^{\lambda_it}$. ∎
$$
e^{At} = Xe^{\Lambda t}X^{-1}. \tag{3.13}
$$
Check against (3.11): $e^{At}\mathbf{u}_0 = Xe^{\Lambda t}X^{-1}\mathbf{u}_0 = Xe^{\Lambda t}\mathbf{c} = \sum_ic_ie^{\lambda_it}\mathbf{x}_i$. Same formula.

▶ **Try it first (b).** Differentiate the series term by term. Reindex.

**Derivation of (b).**

1. Each entry of $e^{At} = \sum_kA^kt^k/k!$ is a power series in $t$ that converges for every $t$ (by the bound above with $\|At\|_F = |t|\|A\|_F$). A power series with infinite radius of convergence may be differentiated term by term.
2. For $k \ge 1$: $\dfrac{d}{dt}\dfrac{A^kt^k}{k!} = \dfrac{A^k\,k\,t^{k-1}}{k!} = \dfrac{A^kt^{k-1}}{(k-1)!}$. The $k = 0$ term is the constant $I$ and differentiates to $0$.
3. Sum and factor one $A$ out on the left: $\displaystyle\sum_{k\ge1}\frac{A^kt^{k-1}}{(k-1)!} = A\sum_{k\ge1}\frac{A^{k-1}t^{k-1}}{(k-1)!} = A\sum_{m\ge0}\frac{A^mt^m}{m!} = Ae^{At}$. Factoring the $A$ out on the right instead gives $e^{At}A$. ∎

Then $\mathbf{u}(t) = e^{At}\mathbf{u}_0$ has $\dot{\mathbf{u}} = Ae^{At}\mathbf{u}_0 = A\mathbf{u}$ and $\mathbf{u}(0) = e^{0}\mathbf{u}_0 = I\mathbf{u}_0 = \mathbf{u}_0$. This is **the** solution of (3.10), and it needed no eigenvectors at all. When $A$ is defective, this is how you solve the system (Example 4.5).

▶ **Try it first (c).** Multiply the two series and collect terms by total degree $j + k = m$. Where do you need $AB = BA$?

**Derivation of (c).**

1. $e^Ae^B = \Big(\sum_j\dfrac{A^j}{j!}\Big)\Big(\sum_k\dfrac{B^k}{k!}\Big) = \sum_{j,k}\dfrac{A^jB^k}{j!\,k!}$, rearranged freely because both series converge absolutely.
2. Group by $m = j + k$: $\displaystyle\sum_{m\ge0}\sum_{j=0}^m\frac{A^jB^{m-j}}{j!\,(m-j)!} = \sum_{m\ge0}\frac{1}{m!}\sum_{j=0}^m\binom{m}{j}A^jB^{m-j}$.
3. The inner sum is the binomial expansion of $(A + B)^m$ **provided $AB = BA$**: expanding $(A+B)^m$ produces words like $ABAB$, and only commutativity lets you reorder each word into $A^jB^{m-j}$. So the sum is $\sum_m(A+B)^m/m! = e^{A+B}$. ∎

**Counterexample without commutativity.** $A = \begin{pmatrix}0&1\\0&0\end{pmatrix}$, $B = \begin{pmatrix}0&0\\1&0\end{pmatrix}$. Here $AB = \begin{pmatrix}1&0\\0&0\end{pmatrix} \ne \begin{pmatrix}0&0\\0&1\end{pmatrix} = BA$. Since $A^2 = B^2 = 0$, the series stop early: $e^A = I + A$, $e^B = I + B$, so
$$
e^Ae^B = (I + A)(I + B) = I + A + B + AB = \begin{pmatrix}2 & 1\\ 1 & 1\end{pmatrix}.
$$
But $S = A + B = \begin{pmatrix}0&1\\1&0\end{pmatrix}$ has $S^2 = I$, so even powers give $I$ and odd powers give $S$:
$$
e^{S} = \Big(\sum_{k\text{ even}}\tfrac{1}{k!}\Big)I + \Big(\sum_{k\text{ odd}}\tfrac1{k!}\Big)S = (\cosh 1)I + (\sinh 1)S = \begin{pmatrix}1.5431 & 1.1752\\ 1.1752 & 1.5431\end{pmatrix}.
$$
Not equal. The project checks both to machine precision. Note also $e^Ae^B \ne e^Be^A$ here.

*Proof of (d).* $A$ commutes with $-A$, so by (c) $e^{At}e^{-At} = e^{At - At} = e^{0} = I$. Likewise $As$ and $At$ commute, giving $e^{A(s+t)} = e^{As}e^{At}$. ∎

■ **What just happened.** The exponential of a matrix is defined by the one thing you know how to do with any matrix: multiply it by itself. Everything you expect from $e^x$ survives, except the law $e^{a+b} = e^ae^b$, which survives only for commuting matrices. QM-4's time-evolution operator $e^{-iHt/\hbar}$, QC-1's gates $e^{-i\theta\sigma/2}$, and LA-11's "functions of matrices" are all this definition.

### 3.10 Second-order systems: $M\ddot{\mathbf{x}} = -K\mathbf{x}$

Newton's law for coupled masses and springs reads
$$
M\ddot{\mathbf{x}} = -K\mathbf{x}, \tag{3.14}
$$
with $M$ the **mass matrix** (diagonal, positive, in our examples) and $K$ the symmetric **stiffness matrix**. $\mathbf{x}(t)$ holds the displacements from equilibrium. For one mass this is $m\ddot x = -kx$, solved by $e^{\pm i\omega t}$ with $\omega = \sqrt{k/m}$.

▶ **Try it first.** Try $\mathbf{x}(t) = \mathbf{v}e^{i\omega t}$ with a constant vector $\mathbf{v}$. Substitute. What equation must $\mathbf{v}$ and $\omega$ satisfy?

**Derivation.**

1. $\ddot{\mathbf{x}} = (i\omega)^2\mathbf{v}e^{i\omega t} = -\omega^2\mathbf{v}e^{i\omega t}$. Substituting: $-\omega^2M\mathbf{v}e^{i\omega t} = -K\mathbf{v}e^{i\omega t}$. Cancel $e^{i\omega t} \ne 0$ and the signs:
$$
K\mathbf{v} = \omega^2M\mathbf{v}. \tag{3.15}
$$
2. This is the **generalized eigenvalue problem**. When $M$ is invertible it is the ordinary problem $(M^{-1}K)\mathbf{v} = \omega^2\mathbf{v}$: $\omega^2$ is an eigenvalue of $M^{-1}K$ and $\mathbf{v}$ its eigenvector. A nonzero $\mathbf{v}$ exists exactly when $\det(K - \omega^2M) = 0$, a polynomial of degree $n$ in $\omega^2$.
3. For each solution $(\omega_j^2, \mathbf{v}_j)$ with $\omega_j^2 > 0$, both $\mathbf{v}_je^{i\omega_jt}$ and $\mathbf{v}_je^{-i\omega_jt}$ solve (3.14), so by superposition so do their real combinations $\mathbf{v}_j\cos\omega_jt$ and $\mathbf{v}_j\sin\omega_jt$.
4. With $n$ independent $\mathbf{v}_j$ the general real solution is
$$
\mathbf{x}(t) = \sum_{j=1}^n\mathbf{v}_j\big(a_j\cos\omega_jt + b_j\sin\omega_jt\big), \tag{3.16}
$$
whose $2n$ constants are fixed by $\mathbf{x}(0) = \sum_ja_j\mathbf{v}_j$ and $\dot{\mathbf{x}}(0) = \sum_jb_j\omega_j\mathbf{v}_j$. ∎

■ **What just happened.** Each $\mathbf{v}_j$ is a **normal mode**: a shape in which every mass oscillates at the single frequency $\omega_j$, all in step. Every motion of the system is a superposition of normal modes. That $\omega_j^2$ is real and positive follows from $K$ being symmetric positive definite and $M$ symmetric positive definite. Exercise B6(b) proves reality, and positivity is one more line from the same computation: B6(b) gives $\omega^2 = \bar{\mathbf{v}}^\top K\mathbf{v}\,/\,\bar{\mathbf{v}}^\top M\mathbf{v}$, and the lemma proved there for $M$ applies equally to $K$, so numerator and denominator are both positive. Exercise B6 also shows the modes are orthogonal in the $M$-weighted inner product, $\mathbf{v}_i^\top M\mathbf{v}_j = 0$, not in the ordinary one. CM-6 derives (3.14) from a Lagrangian for *any* system near a stable equilibrium and then does exactly this. Example 4.3 does the two-mass case.

### 3.11 Complex eigenvalues

**Rotation.** $R(\theta) = \begin{pmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{pmatrix}$ rotates the plane by $\theta$. For $\theta$ not a multiple of $\pi$, no real vector keeps its direction, so there can be no real eigenvector. There are complex ones.

▶ **Try it first.** Write (3.3) for $R(\theta)$. Solve. You should recognize Euler's formula.

**Derivation.**

1. $T = 2\cos\theta$, $D = \cos^2\theta + \sin^2\theta = 1$. So $\lambda^2 - 2\cos\theta\,\lambda + 1 = 0$.
2. $\lambda = \cos\theta \pm\sqrt{\cos^2\theta - 1} = \cos\theta\pm\sqrt{-\sin^2\theta} = \cos\theta \pm i\sin\theta$:
$$
\lambda = e^{\pm i\theta}. \tag{3.17}
$$
3. Eigenvector for $e^{i\theta}$: $R - e^{i\theta}I = \begin{pmatrix}-i\sin\theta & -\sin\theta\\ \sin\theta & -i\sin\theta\end{pmatrix}$. Row 1: $-i\sin\theta\,x_1 - \sin\theta\,x_2 = 0$, so $x_2 = -ix_1$: $\mathbf{x} = (1, -i)^\top$. For $e^{-i\theta}$: $\mathbf{x} = (1, i)^\top$. ∎

■ **What just happened.** $|\lambda| = 1$ because rotations preserve length; the *angle* of the eigenvalue is the rotation angle. The eigenvectors $(1, \mp i)$ don't depend on $\theta$: all rotations of the plane share eigenvectors, which is the matrix statement of the fact that they commute.

**Theorem 3.16 (Conjugate pairs).** If $A$ is real and $A\mathbf{x} = \lambda\mathbf{x}$, then $A\bar{\mathbf{x}} = \bar\lambda\bar{\mathbf{x}}$.

*Proof.* Conjugate both sides: $\overline{A\mathbf{x}} = \bar A\bar{\mathbf{x}} = A\bar{\mathbf{x}}$ since $A$ is real, and $\overline{\lambda\mathbf{x}} = \bar\lambda\bar{\mathbf{x}}$. ∎

So for a real matrix, complex eigenvalues come in conjugate pairs with conjugate eigenvectors. You compute one; the other is free.

**Real solutions from complex ones.** Suppose $A$ is real with eigenvalue $\lambda = a + ib$ ($b \ne 0$) and eigenvector $\mathbf{x} = \mathbf{r} + i\mathbf{s}$ ($\mathbf{r}, \mathbf{s}$ real vectors). Then $e^{\lambda t}\mathbf{x}$ solves $\dot{\mathbf{u}} = A\mathbf{u}$ but is complex. You want real solutions.

▶ **Try it first.** Write a complex solution as $\mathbf{u} = \mathbf{p}(t) + i\mathbf{q}(t)$ with $\mathbf{p}, \mathbf{q}$ real. Substitute into $\dot{\mathbf{u}} = A\mathbf{u}$ and compare real and imaginary parts.

**Derivation.**

1. $\dot{\mathbf{p}} + i\dot{\mathbf{q}} = A\mathbf{p} + iA\mathbf{q}$. Since $A$, $\mathbf{p}$, $\mathbf{q}$ are real, $A\mathbf{p}$ and $A\mathbf{q}$ are real, so equating real and imaginary parts gives $\dot{\mathbf{p}} = A\mathbf{p}$ and $\dot{\mathbf{q}} = A\mathbf{q}$:
$$
\text{if } \mathbf{u}(t) \text{ is a complex solution, then } \operatorname{Re}\mathbf{u}(t) \text{ and } \operatorname{Im}\mathbf{u}(t) \text{ are real solutions.} \tag{3.18}
$$
2. Expand $e^{\lambda t}\mathbf{x}$ with Euler: $e^{(a+ib)t}(\mathbf{r} + i\mathbf{s}) = e^{at}(\cos bt + i\sin bt)(\mathbf{r} + i\mathbf{s}) = e^{at}\big[(\mathbf{r}\cos bt - \mathbf{s}\sin bt) + i(\mathbf{r}\sin bt + \mathbf{s}\cos bt)\big]$.
3. Two real solutions:
$$
\mathbf{p}(t) = e^{at}(\mathbf{r}\cos bt - \mathbf{s}\sin bt), \qquad \mathbf{q}(t) = e^{at}(\mathbf{r}\sin bt + \mathbf{s}\cos bt). \tag{3.19}
$$
4. They are independent: $\mathbf{p}(0) = \mathbf{r}$ and $\mathbf{q}(0) = \mathbf{s}$, and $\mathbf{r}, \mathbf{s}$ are independent. (If $\mathbf{s} = c\mathbf{r}$ then $\mathbf{x} = (1 + ic)\mathbf{r}$, and dividing $A\mathbf{x} = \lambda\mathbf{x}$ by $1 + ic$ gives $A\mathbf{r} = \lambda\mathbf{r}$ with $\mathbf{r}$ real and nonzero; then $\lambda = (A\mathbf{r})_j/r_j$ for any nonzero component $r_j$ is a ratio of real numbers, contradicting $b \ne 0$. The case $\mathbf{r} = c\mathbf{s}$ is the same.) The conjugate pair $\bar\lambda, \bar{\mathbf{x}}$ gives $\mathbf{p}$ and $-\mathbf{q}$: nothing new. So one complex eigenpair supplies exactly the two real solutions it should. ∎

■ **What just happened.** $a = \operatorname{Re}\lambda$ controls growth or decay through $e^{at}$; $b = \operatorname{Im}\lambda$ controls the rate of rotation through $\cos bt$, $\sin bt$. Trajectories spiral in ($a < 0$), spiral out ($a > 0$), or trace closed ellipses ($a = 0$). This is **exactly what to do with complex eigenvalues**: compute one eigenvector, split $e^{\lambda t}\mathbf{x}$ into real and imaginary parts, and those are your two real solutions. Example 4.4 does it with numbers.

**Rotation with decay, in closed form.** The simplest matrix with eigenvalues $a \pm ib$ is
$$
A = \begin{pmatrix}a & -b\\ b & a\end{pmatrix} = aI + bJ, \qquad J = \begin{pmatrix}0 & -1\\ 1 & 0\end{pmatrix} = R(\tfrac\pi2), \qquad J^2 = -I.
$$
(By (3.3): $T = 2a$, $D = a^2 + b^2$, $\lambda = a \pm ib$.)

▶ **Try it first.** $aI$ and $bJ$ commute. Use Theorem 3.15(c), then compute $e^{J\theta}$ from the series using $J^2 = -I$.

**Derivation.**

1. $e^{At} = e^{aIt + bJt} = e^{at}\,e^{bJt}$ since $aI$ commutes with everything.
2. Powers of $J$ cycle: $J^0 = I$, $J^1 = J$, $J^2 = -I$, $J^3 = -J$, $J^4 = I$, $\dots$ So
$$
e^{J\theta} = \Big(1 - \tfrac{\theta^2}{2!} + \tfrac{\theta^4}{4!} - \cdots\Big)I + \Big(\theta - \tfrac{\theta^3}{3!} + \tfrac{\theta^5}{5!} - \cdots\Big)J = (\cos\theta)I + (\sin\theta)J = R(\theta).
$$
3. Therefore
$$
e^{At} = e^{at}R(bt) = e^{at}\begin{pmatrix}\cos bt & -\sin bt\\ \sin bt & \cos bt\end{pmatrix}. \tag{3.20}
$$
∎

■ **What just happened.** Solutions of $\dot{\mathbf{u}} = (aI + bJ)\mathbf{u}$ rotate at angular speed $b$ while scaling by $e^{at}$: a logarithmic spiral. Two independent derivations now agree on the eigenvalues of a rotation: (3.17) directly, and (3.20) plus Theorem 3.15(a) ($J$ has eigenvalues $\pm i$, so $e^{J\theta}$ has eigenvalues $e^{\pm i\theta}$). Every real $2\times2$ matrix with complex eigenvalues is similar to some $aI + bJ$ (§3.13, Exercise B5), so every such phase portrait is this spiral seen through a linear change of coordinates.

### 3.12 Symmetric matrices: real eigenvalues, orthogonal eigenvectors

**Theorem 3.17 (Real symmetric matrices).** Let $A$ be real with $A^\top = A$. Then (a) every eigenvalue of $A$ is real, and (b) eigenvectors belonging to **different** eigenvalues are orthogonal.

▶ **Try it first (a).** Let $A\mathbf{x} = \lambda\mathbf{x}$ with $\mathbf{x} \in \mathbb{C}^n$ possibly complex. Compute the number $\bar{\mathbf{x}}^\top A\mathbf{x}$ two ways: once using $A\mathbf{x} = \lambda\mathbf{x}$, once by first moving $A$ onto $\bar{\mathbf{x}}$.

**Proof of (a).**

1. $\bar{\mathbf{x}}^\top A\mathbf{x} = \bar{\mathbf{x}}^\top(\lambda\mathbf{x}) = \lambda\,\bar{\mathbf{x}}^\top\mathbf{x} = \lambda\|\mathbf{x}\|^2$, where $\|\mathbf{x}\|^2 = \sum_i|x_i|^2 > 0$.
2. Conjugate $A\mathbf{x} = \lambda\mathbf{x}$ (Theorem 3.16): $A\bar{\mathbf{x}} = \bar\lambda\bar{\mathbf{x}}$. Transpose: $\bar{\mathbf{x}}^\top A^\top = \bar\lambda\bar{\mathbf{x}}^\top$. Use $A^\top = A$: $\bar{\mathbf{x}}^\top A = \bar\lambda\bar{\mathbf{x}}^\top$. Multiply on the right by $\mathbf{x}$: $\bar{\mathbf{x}}^\top A\mathbf{x} = \bar\lambda\|\mathbf{x}\|^2$.
3. Subtract: $(\lambda - \bar\lambda)\|\mathbf{x}\|^2 = 0$, and $\|\mathbf{x}\|^2 \ne 0$, so $\lambda = \bar\lambda$: $\lambda$ is real. ∎

Since $\lambda$ is real, $A - \lambda I$ is a real singular matrix and has a real null vector: **the eigenvectors of a real symmetric matrix can always be taken real.**

**Proof of (b).** Let $A\mathbf{x}_1 = \lambda_1\mathbf{x}_1$, $A\mathbf{x}_2 = \lambda_2\mathbf{x}_2$, $\lambda_1 \ne \lambda_2$, both vectors real.

1. $\lambda_1\,\mathbf{x}_1^\top\mathbf{x}_2 = (\lambda_1\mathbf{x}_1)^\top\mathbf{x}_2 = (A\mathbf{x}_1)^\top\mathbf{x}_2 = \mathbf{x}_1^\top A^\top\mathbf{x}_2 = \mathbf{x}_1^\top A\mathbf{x}_2 = \mathbf{x}_1^\top(\lambda_2\mathbf{x}_2) = \lambda_2\,\mathbf{x}_1^\top\mathbf{x}_2$.
2. So $(\lambda_1 - \lambda_2)\,\mathbf{x}_1^\top\mathbf{x}_2 = 0$, and $\lambda_1 \ne \lambda_2$ forces $\mathbf{x}_1^\top\mathbf{x}_2 = 0$. ∎

■ **What just happened.** Both proofs have the same shape: move $A$ from one side of an inner product to the other using $A^\top = A$, and compare. This is Cauchy's 1829 theorem, the one that settled the secular equation. It is why the inertia tensor has perpendicular principal axes (CM-7), why the normal-mode frequencies in §3.10 are real (B6), and why a Hamiltonian's energy levels are real (QM-4, where "symmetric" becomes "Hermitian," $\bar A^\top = A$, and the same proof works with one change: in step 2, Theorem 3.16 needs $A$ real, so instead take the conjugate transpose of $A\mathbf{x} = \lambda\mathbf{x}$ to get $\bar{\mathbf{x}}^\top\bar A^\top = \bar\lambda\bar{\mathbf{x}}^\top$ and use $\bar A^\top = A$). **LA-8 (the spectral theorem)** completes the picture: a real symmetric matrix always has an *orthonormal basis* of eigenvectors, $A = Q\Lambda Q^\top$, *even when eigenvalues repeat*. What we proved here is the case of distinct eigenvalues; the repeated case needs one more idea.

### 3.13 Similar matrices

**Definition.** $B$ is **similar** to $A$ if $B = M^{-1}AM$ for some invertible $M$. (LA-4: $A$ and $B$ represent the same linear map in two bases, and $M$ is the change-of-basis matrix.)

**Theorem 3.18.** Similar matrices have the same characteristic polynomial, hence the same eigenvalues with the same multiplicities, the same trace, and the same determinant. If $A\mathbf{x} = \lambda\mathbf{x}$ then $B(M^{-1}\mathbf{x}) = \lambda(M^{-1}\mathbf{x})$.

*Proof.* $\det(B - \lambda I) = \det(M^{-1}AM - \lambda M^{-1}IM) = \det\big(M^{-1}(A - \lambda I)M\big) = \det(M^{-1})\det(A - \lambda I)\det(M) = \det(A - \lambda I)$, using the product rule (LA-6, Theorem 3.10) and $\det M^{-1} = 1/\det M$. For the eigenvector: $BM^{-1}\mathbf{x} = M^{-1}AMM^{-1}\mathbf{x} = M^{-1}A\mathbf{x} = \lambda M^{-1}\mathbf{x}$. ∎

■ **What just happened.** Eigenvalues belong to the linear map, not to the matrix; changing basis only re-coordinatizes the eigenvectors ($\mathbf{x} \mapsto M^{-1}\mathbf{x}$). Diagonalization (3.6) is the statement "$A$ is similar to a diagonal matrix." **The converse of Theorem 3.18 is false:** $\begin{pmatrix}1&1\\0&1\end{pmatrix}$ and $I$ have the same characteristic polynomial $(1-\lambda)^2$ but are not similar, since anything similar to $I$ is $M^{-1}IM = I$. The complete invariant for similarity is the Jordan form, LA-11.

### 3.14 Repeated eigenvalues: algebraic vs geometric multiplicity

**Definitions.** The **algebraic multiplicity** $\operatorname{AM}(\lambda)$ is the multiplicity of $\lambda$ as a root of $p_A$. The **geometric multiplicity** is $\operatorname{GM}(\lambda) = \dim E_\lambda = \dim N(A - \lambda I) = n - \operatorname{rank}(A - \lambda I)$: the number of independent eigenvectors for $\lambda$.

**Theorem 3.19.** For every eigenvalue, $1 \le \operatorname{GM}(\lambda) \le \operatorname{AM}(\lambda)$.

▶ **Try it first.** Pick a basis of $E_\lambda$, extend it to a basis of $\mathbb{C}^n$, and look at $M^{-1}AM$, where $M$ has that basis as its columns. What are its first $g$ columns? What does that force on its characteristic polynomial, which is also $p_A$ (Theorem 3.18)?

*Proof.*

1. $\operatorname{GM} \ge 1$ because $\lambda$ is an eigenvalue, so $E_\lambda$ contains a nonzero vector.
2. Let $g = \operatorname{GM}(\lambda)$. Choose a basis $\mathbf{x}_1, \dots, \mathbf{x}_g$ of $E_\lambda$ and extend it to a basis $\mathbf{x}_1, \dots, \mathbf{x}_n$ of $\mathbb{C}^n$ (LA-3). Let $M$ have these as columns; $M$ is invertible.
3. Let $B = M^{-1}AM$. For $j \le g$, column $j$ of $B$ is $B\mathbf{e}_j = M^{-1}AM\mathbf{e}_j = M^{-1}A\mathbf{x}_j = \lambda M^{-1}\mathbf{x}_j = \lambda\mathbf{e}_j$. So $B = \begin{pmatrix}\lambda I_g & C\\ 0 & D\end{pmatrix}$, block upper triangular with an $(n-g)\times(n-g)$ block $D$.
4. *Block-triangular determinants factor:* for $B = \begin{pmatrix}P & C\\ 0 & D\end{pmatrix}$ with $P$ of size $g\times g$, $\det B = \det P\det D$ (this is LA-6, Exercise B1(a), proved there from uniqueness). Reason, from the Leibniz formula (LA-6, Theorem 3.9, eq. 3.9), whose term for $\sigma$ is $\operatorname{sgn}(\sigma)\,b_{1\sigma(1)}\cdots b_{n\sigma(n)}$. The zero block is rows $g+1, \dots, n$ by columns $1, \dots, g$. If $\sigma$ sends some lower index $i > g$ to a column $\sigma(i) \le g$, the factor $b_{i\sigma(i)}$ is in the zero block and the term vanishes. The surviving $\sigma$ map $\{g+1, \dots, n\}$ into itself; being bijections, they then also map $\{1, \dots, g\}$ into itself. So each surviving $\sigma$ is a pair $(\sigma_1, \sigma_2)$, a permutation of each block. No pair of indices with one in each block is an inversion, so $\operatorname{sgn}\sigma = \operatorname{sgn}\sigma_1\operatorname{sgn}\sigma_2$, and the sum factors as $\det P\cdot\det D$.
5. Hence $p_B(\mu) = \det(B - \mu I) = \det\big((\lambda - \mu)I_g\big)\det(D - \mu I) = (\lambda - \mu)^g\,p_D(\mu)$. So $(\lambda - \mu)^g$ divides $p_B = p_A$ (Theorem 3.18): $\operatorname{AM}(\lambda) \ge g$. ∎

**Theorem 3.20.** $A$ is diagonalizable if and only if $\operatorname{GM}(\lambda) = \operatorname{AM}(\lambda)$ for every eigenvalue.

*Proof.* ($\Leftarrow$) Take a basis of each eigenspace; together these are $\sum_\lambda\operatorname{GM}(\lambda) = \sum_\lambda\operatorname{AM}(\lambda) = n$ vectors. They are independent: suppose a combination is zero and group it by eigenspace as $\mathbf{w}_1 + \cdots + \mathbf{w}_r = \mathbf{0}$ with $\mathbf{w}_i \in E_{\lambda_i}$. If any $\mathbf{w}_i \ne \mathbf{0}$, the nonzero ones are eigenvectors for distinct eigenvalues that are dependent, contradicting Theorem 3.6. So every $\mathbf{w}_i = \mathbf{0}$, and since each $\mathbf{w}_i$ is a combination of a basis of $E_{\lambda_i}$, all coefficients vanish. $n$ independent eigenvectors: diagonalizable by Theorem 3.8. ($\Rightarrow$) If $A = X\Lambda X^{-1}$, the columns of $X$ belonging to $\lambda$ are $\operatorname{AM}(\lambda)$ independent vectors in $E_\lambda$ (the number of times $\lambda$ appears on the diagonal of $\Lambda$ is its multiplicity as a root of $p_\Lambda = p_A$). So $\operatorname{GM} \ge \operatorname{AM}$, and with Theorem 3.19, equality. ∎

■ **What just happened.** The algebraic multiplicity is a ceiling on how many independent eigenvectors an eigenvalue can supply, and Theorem 3.19 says no eigenvalue ever exceeds it. The AMs add up to $n$, so you get a full basis of eigenvectors exactly when every eigenvalue reaches its ceiling. A simple eigenvalue (AM $= 1$) always does, since GM $\ge 1$; when every eigenvalue is simple, that is Corollary 3.7 again. Trouble can only come from repeated eigenvalues, and the test is one rank per repeated eigenvalue, $\operatorname{GM}(\lambda) = n - \operatorname{rank}(A - \lambda I)$.

**Defective matrices.** $A$ is **defective** if $\operatorname{GM}(\lambda) < \operatorname{AM}(\lambda)$ for some $\lambda$: it has too few eigenvectors and cannot be diagonalized. The standard example is the shear
$$
S = \begin{pmatrix}1 & 1\\ 0 & 1\end{pmatrix}: \quad p_S(\lambda) = (1-\lambda)^2,\ \operatorname{AM}(1) = 2;\qquad S - I = \begin{pmatrix}0&1\\0&0\end{pmatrix} \text{ has rank } 1,\ \operatorname{GM}(1) = 2 - 1 = 1.
$$
The only eigenvector direction is $\mathbf{e}_1$ (the $x$-axis, which the shear fixes). Every other vector gets tilted.

*What goes wrong.* Formula (3.8) needs a basis of eigenvectors; there isn't one. *What still works.* Direct computation: $S^k = \begin{pmatrix}1 & k\\ 0 & 1\end{pmatrix}$ (induction: $S^kS$ adds $1$ to the top-right entry). That is **polynomial** growth $k$ alongside the eigenvalue power $1^k$, something (3.8) can never produce. For the exponential, write $S = I + N$ with $N = \begin{pmatrix}0&1\\0&0\end{pmatrix}$, $N^2 = 0$, and $IN = NI$; then by Theorem 3.15(c)
$$
e^{St} = e^{It}e^{Nt} = e^t(I + Nt) = e^t\begin{pmatrix}1 & t\\ 0 & 1\end{pmatrix}.
$$
There is the $te^{\lambda t}$ you met in ODEs when the characteristic equation had a double root. **This is where it comes from.** Example 4.5 shows it in a critically damped oscillator. LA-11 (Theorem 3.10, the Jordan canonical form) proves that every matrix is $X(\Lambda + N)X^{-1}$ with $N$ nilpotent and commuting with $\Lambda$ (the Jordan form), so that $e^{At} = Xe^{\Lambda t}(I + Nt + N^2t^2/2 + \cdots)X^{-1}$ with finitely many terms, and supplies the missing directions as *generalized eigenvectors*.

Numerically, defectiveness is not an isolated accident: a matrix *near* a defective one has nearly parallel eigenvectors, so $X$ is nearly singular, $X^{-1}$ is huge, and computing with (3.6) in floating point loses digits: about $\log_{10}(1/\varepsilon)$ of them when the matrix is within $\varepsilon$ of a defective one, so all of them as $\varepsilon$ approaches machine precision (Exercise B7(d); LA-13).

### 3.15 The trace–determinant plane: classifying $2\times2$ phase portraits

Let $A$ be real $2\times2$ with $T = \operatorname{tr}A$, $D = \det A$, eigenvalues (3.4): $\lambda = \frac{T\pm\sqrt{T^2 - 4D}}{2}$. The phase portrait of $\dot{\mathbf{u}} = A\mathbf{u}$ is determined by $(T, D)$.

▶ **Try it first.** Draw axes $T$ (horizontal) and $D$ (vertical). Shade where the eigenvalues are real, where they are complex, where both are negative. Use $\lambda_1 + \lambda_2 = T$ and $\lambda_1\lambda_2 = D$.

**Derivation.**

1. *Discriminant.* $\Delta = T^2 - 4D$. $\Delta > 0$ (below the parabola $D = T^2/4$): two distinct real eigenvalues. $\Delta < 0$ (above it): complex pair $\lambda = \tfrac T2 \pm i\tfrac{\sqrt{-\Delta}}{2}$. $\Delta = 0$: a repeated real eigenvalue $T/2$.
2. *Real eigenvalues, $D < 0$.* Since $\lambda_1\lambda_2 = D < 0$, the eigenvalues have **opposite signs** (and $\Delta = T^2 - 4D > 0$ automatically). Solutions $c_1e^{\lambda_1t}\mathbf{x}_1 + c_2e^{\lambda_2t}\mathbf{x}_2$ decay along one eigenline and grow along the other: a **saddle**. Unstable for every $T$.
3. *Real eigenvalues, $D > 0$, $\Delta > 0$.* Same sign, and $\lambda_1 + \lambda_2 = T$ tells which. $T < 0$: both negative, every solution decays, a **stable node** (sink). $T > 0$: both positive, an **unstable node** (source). Trajectories leave or approach the origin tangent to the *slow* eigenvector (smaller $|\lambda|$), since the fast component dies first, and run parallel to the *fast* one far away.
4. *Complex eigenvalues, $\Delta < 0$* (which forces $D > T^2/4 \ge 0$). $\operatorname{Re}\lambda = T/2$. $T < 0$: **stable spiral**. $T > 0$: **unstable spiral**. $T = 0$: $\lambda = \pm i\sqrt D$, a **center**, closed ellipses, like the undamped oscillator.
5. *Boundaries.* On the parabola $\Delta = 0$: a repeated eigenvalue, a **star** if $\operatorname{GM} = 2$ (then $A = \tfrac T2 I$) or a **degenerate node** if defective (Example 4.5). On $D = 0$: a zero eigenvalue, a whole line of equilibria (the whole plane if $A = 0$). On $T = 0$, $D > 0$: centers, a knife edge between stable and unstable spirals.
6. *Stability summary.* All $\operatorname{Re}\lambda < 0$ $\iff$ $T < 0$ and $D > 0$. (Real case: two negatives have negative sum and positive product, and conversely. Complex case: $\operatorname{Re}\lambda = T/2 < 0$, and $D > 0$ holds automatically.) ∎

```
                              D
                              ^  <- positive D-axis (T = 0, D > 0): centers
          .                   |                   .
           . stable spiral    |  unstable spiral .
             .                |                .
  stable      .               |               .     unstable
  node          .             |             .       node
                  .           |           .
 D = T^2/4 -->      .         |         .
                       .      |      .
  ----------------------------+-----------------------------> T
                              |
           saddles (D < 0) below the T-axis, for every T
```

Above the $T$-axis, the left half ($T < 0$) is stable and the right half ($T > 0$) unstable; below it everything is a saddle, unstable. Spirals lie above the parabola $D = T^2/4$, nodes between the parabola and the $T$-axis, and centers on the positive $D$-axis. The project draws the same plane to scale (`figures/trace_det_plane.png`).

■ **What just happened.** Two numbers you can read off a $2\times2$ matrix without solving anything classify the entire long-run behavior of $\dot{\mathbf{u}} = A\mathbf{u}$. Dynamical systems courses linearize nonlinear systems at equilibria and then do precisely this. The project draws four representative portraits (a saddle, a stable node, a stable spiral, and a center on the boundary $T = 0$), with eigenvectors marked in the real cases.

## 4. Worked examples

### Example 4.1: Weather, a two-state Markov chain

Each day is sunny (S) or rainy (R). After a sunny day: $90\%$ sunny, $10\%$ rainy. After a rainy day: $50\%$ each. Columns are "from," rows are "to":
$$
A = \begin{pmatrix} 0.9 & 0.5\\ 0.1 & 0.5\end{pmatrix}, \qquad \mathbf{u}_k = \begin{pmatrix}P(\text{sunny on day }k)\\ P(\text{rainy on day }k)\end{pmatrix}.
$$
Columns sum to $1$ and entries are nonnegative: Markov.

*Eigenvalues.* $T = 1.4$, $D = 0.45 - 0.05 = 0.40$. By (3.3), $\lambda^2 - 1.4\lambda + 0.4 = 0$, so $\lambda = \frac{1.4\pm\sqrt{1.96 - 1.60}}{2} = \frac{1.4\pm0.6}{2}$: $\lambda_1 = 1$, $\lambda_2 = 0.4$. Theorem 3.11 predicted $1$ and $|\lambda_2| \le 1$. Check: $1 + 0.4 = 1.4 = T$, $1\times0.4 = D$. ✓

*Eigenvectors.* $\lambda_1 = 1$: $A - I = \begin{pmatrix}-0.1 & 0.5\\ 0.1 & -0.5\end{pmatrix}$; row 1 says $x_1 = 5x_2$, so $\mathbf{x}_1 = (5, 1)^\top$. Scaled to sum $1$: the steady state is $\mathbf{u}_\infty = (5/6, 1/6)^\top$. $\lambda_2 = 0.4$: $A - 0.4I = \begin{pmatrix}0.5 & 0.5\\ 0.1 & 0.1\end{pmatrix}$, so $\mathbf{x}_2 = (1, -1)^\top$. Note $\mathbf{x}_2$ sums to zero (B4(b) says it must).

*Start on a rainy day:* $\mathbf{u}_0 = (0, 1)^\top = c_1(5, 1)^\top + c_2(1, -1)^\top$. Rows: $5c_1 + c_2 = 0$ and $c_1 - c_2 = 1$. Adding: $6c_1 = 1$, so $c_1 = 1/6$, $c_2 = -5/6$. By (3.8),
$$
\mathbf{u}_k = \tfrac16\begin{pmatrix}5\\1\end{pmatrix} - \tfrac56(0.4)^k\begin{pmatrix}1\\-1\end{pmatrix} = \begin{pmatrix}\tfrac56 - \tfrac56(0.4)^k\\[2pt] \tfrac16 + \tfrac56(0.4)^k\end{pmatrix}.
$$
$k = 1$: $\tfrac56(0.4) = \tfrac13$, so $\mathbf{u}_1 = (\tfrac12, \tfrac12)^\top$. Direct check: $A(0, 1)^\top = $ second column $= (0.5, 0.5)^\top$. ✓
$k = 2$: $\tfrac56(0.16) = \tfrac{2}{15}$, so $\mathbf{u}_2 = (0.7, 0.3)^\top$. Direct: $A(0.5, 0.5)^\top = (0.45 + 0.25, 0.05 + 0.25)^\top = (0.7, 0.3)^\top$. ✓
$k\to\infty$: $\mathbf{u}_k \to (5/6, 1/6)^\top \approx (0.833, 0.167)^\top$. In the long run five days in six are sunny, whatever today is. The error shrinks by a factor $0.4$ per day; after ten days it is $(5/6)(0.4)^{10} \approx 10^{-4}$.

The project simulates this chain for 20,000 days and compares the empirical frequencies to $(5/6, 1/6)$.

### Example 4.2: Fibonacci in closed form

$F_{k+2} = F_{k+1} + F_k$, $F_0 = 0$, $F_1 = 1$. Package two consecutive terms: $\mathbf{u}_k = (F_{k+1}, F_k)^\top$. Then
$$
\mathbf{u}_{k+1} = \begin{pmatrix}F_{k+2}\\ F_{k+1}\end{pmatrix} = \begin{pmatrix}F_{k+1} + F_k\\ F_{k+1}\end{pmatrix} = \begin{pmatrix}1 & 1\\ 1 & 0\end{pmatrix}\mathbf{u}_k, \qquad \mathbf{u}_0 = \begin{pmatrix}1\\0\end{pmatrix}.
$$
*Eigenvalues.* $T = 1$, $D = -1$: $\lambda^2 - \lambda - 1 = 0$, $\lambda = \frac{1\pm\sqrt5}{2}$. Call them $\varphi = 1.6180\ldots$ (the golden ratio) and $\psi = -0.6180\ldots = 1 - \varphi$. $D < 0$: a saddle; one direction grows, one decays with alternating sign.

*Eigenvectors.* $A - \lambda I = \begin{pmatrix}1-\lambda & 1\\ 1 & -\lambda\end{pmatrix}$; row 2 gives $x_1 = \lambda x_2$, so $\mathbf{x} = (\lambda, 1)^\top$ for each eigenvalue: $\mathbf{x}_1 = (\varphi, 1)^\top$, $\mathbf{x}_2 = (\psi, 1)^\top$.

*Coordinates of $\mathbf{u}_0$.* $(1, 0)^\top = c_1(\varphi, 1)^\top + c_2(\psi, 1)^\top$. Row 2: $c_1 + c_2 = 0$. Row 1: $c_1\varphi + c_2\psi = c_1(\varphi - \psi) = c_1\sqrt5 = 1$. So $c_1 = 1/\sqrt5$, $c_2 = -1/\sqrt5$.

*Closed form.* By (3.8), $\mathbf{u}_k = \frac{1}{\sqrt5}\big[\varphi^k(\varphi, 1)^\top - \psi^k(\psi, 1)^\top\big]$. The second component is
$$
F_k = \frac{\varphi^k - \psi^k}{\sqrt5}.
$$
This is Binet's formula (1843), known to de Moivre a century earlier. Since $|\psi| < 1$, the second term is below $1/2$ in size for every $k \ge 0$, so $F_k$ is the integer nearest $\varphi^k/\sqrt5$, and $F_{k+1}/F_k \to \varphi$. Check: $k = 10$: $\varphi^{10}/\sqrt5 = 55.0036\ldots$, and $F_{10} = 55$. ✓ $F_{100} = 354{,}224{,}848{,}179{,}261{,}915{,}075$; the floating-point closed form gives $3.542248481792631\times10^{20}$, a relative error of about $3\times10^{-15}$: 14 to 15 correct significant digits, not the 16 a double can hold. Almost all of the error comes from rounding $\varphi$ itself to a double: a relative error $\delta$ in $\varphi$ becomes about $100\delta$ in $\varphi^{100}$. The project computes both.

### Example 4.3: Two masses, three springs (mechanics; CM-6 preview)

Wall — spring — $m_1$ — spring — $m_2$ — spring — wall, all three springs with stiffness $k = 1$, masses $m_1 = 1$ and $m_2 = 2$. Let $x_1, x_2$ be displacements from equilibrium (positive to the right).

*Equations of motion.* Force on $m_1$: the left spring pulls back with $-kx_1$; the middle spring is stretched by $x_2 - x_1$ and pulls $m_1$ forward with $+k(x_2 - x_1)$. Total: $-2x_1 + x_2$. Force on $m_2$: middle spring $-k(x_2 - x_1)$, right spring $-kx_2$. Total: $x_1 - 2x_2$. So
$$
\underbrace{\begin{pmatrix}1 & 0\\ 0 & 2\end{pmatrix}}_{M}\ddot{\mathbf{x}} = -\underbrace{\begin{pmatrix}2 & -1\\ -1 & 2\end{pmatrix}}_{K}\mathbf{x}.
$$
*Frequencies.* By (3.15), $K\mathbf{v} = \omega^2M\mathbf{v}$ needs $\det(K - \omega^2M) = 0$:
$$
\det\begin{pmatrix}2 - \omega^2 & -1\\ -1 & 2 - 2\omega^2\end{pmatrix} = (2 - \omega^2)(2 - 2\omega^2) - 1 = 2\omega^4 - 6\omega^2 + 3 = 0,
$$
$$
\omega^2 = \frac{6\pm\sqrt{36 - 24}}{4} = \frac{3\pm\sqrt3}{2}: \qquad \omega_1^2 = \frac{3-\sqrt3}{2} \approx 0.6340,\quad \omega_2^2 = \frac{3+\sqrt3}{2}\approx 2.3660,
$$
so $\omega_1 \approx 0.7962$ and $\omega_2 \approx 1.5382$. Cross-check with $M^{-1}K = \begin{pmatrix}2 & -1\\ -\tfrac12 & 1\end{pmatrix}$: $T = 3 = \omega_1^2 + \omega_2^2$ ✓, $D = 2 - \tfrac12 = \tfrac32 = \omega_1^2\omega_2^2 = \tfrac{9-3}{4}$ ✓.

*Mode shapes.* Row 1 of $(K - \omega^2M)\mathbf{v} = \mathbf{0}$: $(2 - \omega^2)v_1 - v_2 = 0$, so $\mathbf{v} = (1, 2 - \omega^2)^\top$.
Mode 1 (slow): $\mathbf{v}_1 = \big(1, \tfrac{1+\sqrt3}{2}\big)^\top \approx (1, 1.366)^\top$. Both masses move the same way; the heavy one swings farther.
Mode 2 (fast): $\mathbf{v}_2 = \big(1, \tfrac{1-\sqrt3}{2}\big)^\top \approx (1, -0.366)^\top$. The masses move oppositely, stretching the middle spring hard, hence the higher frequency.

*$M$-orthogonality (B6).* $\mathbf{v}_1^\top M\mathbf{v}_2 = 1\cdot1\cdot1 + 2\cdot\tfrac{1+\sqrt3}{2}\cdot\tfrac{1-\sqrt3}{2} = 1 + 2\cdot\tfrac{1 - 3}{4} = 1 - 1 = 0$. ✓ The ordinary dot product is $\mathbf{v}_1^\top\mathbf{v}_2 = 1 + \tfrac{1-3}{4} = \tfrac12 \ne 0$: the modes are *not* orthogonal in the usual sense, because $M^{-1}K$ is not symmetric. The right inner product for mechanics is the one weighted by mass.

*A specific motion.* Pull $m_1$ to $x_1 = 1$, hold $m_2$ at $0$, release from rest: $\mathbf{x}(0) = (1, 0)^\top$, $\dot{\mathbf{x}}(0) = \mathbf{0}$. In (3.16), zero initial velocity forces $b_1 = b_2 = 0$, and $a_1\mathbf{v}_1 + a_2\mathbf{v}_2 = (1, 0)^\top$:
$a_1 + a_2 = 1$ and $a_1\tfrac{1+\sqrt3}{2} + a_2\tfrac{1-\sqrt3}{2} = 0$. The second gives $a_1(1 + \sqrt3) = a_2(\sqrt3 - 1)$, so $a_1 = a_2\frac{\sqrt3-1}{\sqrt3+1} = a_2\frac{(\sqrt3-1)^2}{2} = a_2(2 - \sqrt3)$. Then $a_2(3 - \sqrt3) = 1$:
$$
a_2 = \frac{1}{3-\sqrt3} = \frac{3+\sqrt3}{6} \approx 0.7887, \qquad a_1 = \frac{3-\sqrt3}{6}\approx 0.2113.
$$
$$
\mathbf{x}(t) = 0.2113\begin{pmatrix}1\\1.366\end{pmatrix}\cos(0.7962\,t) + 0.7887\begin{pmatrix}1\\-0.366\end{pmatrix}\cos(1.5382\,t).
$$
The ratio $\omega_2/\omega_1 = \sqrt{(3+\sqrt3)/(3-\sqrt3)} = \sqrt{2+\sqrt3} \approx 1.932$ is irrational, so the motion never exactly repeats; energy sloshes back and forth between the masses. The project checks every number here against `scipy.linalg.eigh(K, M)` and the whole trajectory against `solve_ivp`, and draws both modes as frame sequences.

### Example 4.4: A spiral

$\dot{\mathbf{u}} = A\mathbf{u}$ with $A = \begin{pmatrix}-2 & -5\\ 1 & 0\end{pmatrix}$, $\mathbf{u}(0) = (2, 0)^\top$.

*Classify first (§3.15).* $T = -2$, $D = 5$, $\Delta = 4 - 20 = -16 < 0$, $T < 0$: a stable spiral.

*Eigenvalues.* $\lambda = \frac{-2\pm\sqrt{-16}}{2} = -1\pm2i$. Take $\lambda = -1 + 2i$, so $a = -1$, $b = 2$.

*One complex eigenvector.* $A - \lambda I = \begin{pmatrix}-1-2i & -5\\ 1 & 1-2i\end{pmatrix}$. Row 2: $x_1 + (1 - 2i)x_2 = 0$; with $x_2 = 1$, $x_1 = -1 + 2i$. So $\mathbf{x} = (-1 + 2i, 1)^\top = \mathbf{r} + i\mathbf{s}$ with $\mathbf{r} = (-1, 1)^\top$, $\mathbf{s} = (2, 0)^\top$. (Row 1 check: $(-1-2i)(-1+2i) - 5 = (1 + 4) - 5 = 0$. ✓)

*Two real solutions, from (3.19).*
$$
\mathbf{p}(t) = e^{-t}\big(\mathbf{r}\cos2t - \mathbf{s}\sin2t\big) = e^{-t}\begin{pmatrix}-\cos2t - 2\sin2t\\ \cos2t\end{pmatrix}, \qquad
\mathbf{q}(t) = e^{-t}\big(\mathbf{r}\sin2t + \mathbf{s}\cos2t\big) = e^{-t}\begin{pmatrix}2\cos2t - \sin2t\\ \sin2t\end{pmatrix}.
$$
*Match the initial condition.* $\mathbf{p}(0) = (-1, 1)^\top$, $\mathbf{q}(0) = (2, 0)^\top = \mathbf{u}(0)$. So $\mathbf{u} = \mathbf{q}$:
$$
\mathbf{u}(t) = e^{-t}\begin{pmatrix}2\cos2t - \sin2t\\ \sin2t\end{pmatrix}.
$$
*Verify.* $\dot{\mathbf{u}}(0)$: the first component's derivative at $0$ is $-(2) + (-2) = -4$, the second's is $2$. And $A\mathbf{u}(0) = A(2, 0)^\top = (-4, 2)^\top$. ✓ At $t = 0.7$: $\mathbf{u} \approx (-0.3206, 0.4894)^\top$, which the project confirms against `expm(A t) @ u0`.

*Picture.* Period $2\pi/b = \pi$; each turn shrinks the amplitude by $e^{-\pi}\approx0.043$. Direction of rotation: at $(2, 0)$ the velocity $(-4, 2)$ points up and left, so counterclockwise. General rule for $2\times2$: at $\mathbf{e}_1$ the velocity is the first column $(a_{11}, a_{21})$, so the sign of $a_{21}$ gives the sense of rotation ($a_{21} > 0$: counterclockwise). Here $a_{21} = 1$.

### Example 4.5: A defective matrix: critical damping

The oscillator $\ddot x + 2\dot x + x = 0$ ($m = 1$, damping $c = 2$, $k = 1$; $c^2 = 4mk$ is the critical value). Set $\mathbf{u} = (x, \dot x)^\top$:
$$
\dot{\mathbf{u}} = \begin{pmatrix}\dot x\\ -x - 2\dot x\end{pmatrix} = \underbrace{\begin{pmatrix}0 & 1\\ -1 & -2\end{pmatrix}}_{A}\mathbf{u}.
$$
$T = -2$, $D = 1$, $\Delta = 0$: we are on the parabola of §3.15. $\lambda = -1$, twice.

*Only one eigenvector.* $A + I = \begin{pmatrix}1 & 1\\ -1 & -1\end{pmatrix}$ has rank $1$, so $\operatorname{GM}(-1) = 1 < 2 = \operatorname{AM}(-1)$. The eigenvector is $(1, -1)^\top$. The one-parameter family $ce^{-t}(1, -1)^\top$ cannot match an arbitrary two-component initial condition. Formula (3.11) fails.

*The exponential doesn't care.* Write $A = -I + N$ with $N = A + I = \begin{pmatrix}1 & 1\\ -1 & -1\end{pmatrix}$. Then $N^2 = \begin{pmatrix}1-1 & 1-1\\ -1+1 & -1+1\end{pmatrix} = 0$, and $-I$ commutes with $N$, so by Theorem 3.15(c)
$$
e^{At} = e^{-t}e^{Nt} = e^{-t}(I + Nt) = e^{-t}\begin{pmatrix}1 + t & t\\ -t & 1 - t\end{pmatrix}.
$$
Check: $\frac{d}{dt}e^{At}\big|_{t=0} = -I + N = A$. ✓ Release from $x = 1$ at rest, $\mathbf{u}_0 = (1, 0)^\top$:
$$
\mathbf{u}(t) = e^{-t}\begin{pmatrix}1 + t\\ -t\end{pmatrix}, \qquad x(t) = (1 + t)e^{-t}.
$$
Check: $\dot x = e^{-t} - (1+t)e^{-t} = -te^{-t}$, which is the second component. ✓ Your ODE course's $(c_1 + c_2t)e^{-t}$ for a double root is this computation. The project compares $e^{At}$ at $t = 1.3$ with `scipy.linalg.expm`.

### Example 4.6: A $3\times3$ symmetric matrix, fully

$$
A = \begin{pmatrix}3 & 1 & -1\\ 1 & 3 & -1\\ -1 & -1 & 5\end{pmatrix}.
$$
Symmetric, so Theorem 3.17 promises real eigenvalues and orthogonal eigenvectors. Let's see them appear.

*Characteristic polynomial*, expanding $\det(A - \lambda I)$ along row 1 (LA-6, Theorem 3.12):
$$
p_A(\lambda) = (3-\lambda)\begin{vmatrix}3-\lambda & -1\\ -1 & 5-\lambda\end{vmatrix} - 1\begin{vmatrix}1 & -1\\ -1 & 5-\lambda\end{vmatrix} + (-1)\begin{vmatrix}1 & 3-\lambda\\ -1 & -1\end{vmatrix}.
$$
The minors: $(3-\lambda)(5-\lambda) - 1 = \lambda^2 - 8\lambda + 14$; $\;(5 - \lambda) - 1 = 4 - \lambda$; $\;-1 + (3 - \lambda) = 2 - \lambda$. So
$$
p_A(\lambda) = (3-\lambda)(\lambda^2 - 8\lambda + 14) - (4 - \lambda) - (2 - \lambda) = -\lambda^3 + 11\lambda^2 - 36\lambda + 36.
$$
(Expand: $(3-\lambda)(\lambda^2 - 8\lambda + 14) = -\lambda^3 + 11\lambda^2 - 38\lambda + 42$; then subtract $6 - 2\lambda$.) Checks from Theorem 3.4: the $\lambda^2$ coefficient $11 = 3 + 3 + 5 = \operatorname{tr}A$ ✓; the constant $36 = \det A$ (compute: $3(15-1) - 1(5-1) + (-1)(-1+3) = 42 - 4 - 2 = 36$ ✓).

*Roots.* Integer roots divide $36$. Try $\lambda = 2$: $-8 + 44 - 72 + 36 = 0$. ✓ Divide: $-\lambda^3 + 11\lambda^2 - 36\lambda + 36 = -(\lambda - 2)(\lambda^2 - 9\lambda + 18) = -(\lambda-2)(\lambda-3)(\lambda-6)$. Eigenvalues $2, 3, 6$. Sum $11$ ✓, product $36$ ✓.

*Eigenvectors*, by elimination on $A - \lambda I$:
- $\lambda = 2$: $\begin{pmatrix}1 & 1 & -1\\ 1 & 1 & -1\\ -1 & -1 & 3\end{pmatrix}$. Rows 1 and 2 agree; row 3 + row 1 gives $(0, 0, 2)$, so $x_3 = 0$, then $x_1 + x_2 = 0$: $\mathbf{x}_1 = (1, -1, 0)^\top$.
- $\lambda = 3$: $\begin{pmatrix}0 & 1 & -1\\ 1 & 0 & -1\\ -1 & -1 & 2\end{pmatrix}$. Row 1: $x_2 = x_3$; row 2: $x_1 = x_3$: $\mathbf{x}_2 = (1, 1, 1)^\top$.
- $\lambda = 6$: $\begin{pmatrix}-3 & 1 & -1\\ 1 & -3 & -1\\ -1 & -1 & -1\end{pmatrix}$. Row 1 $-$ row 2: $-4x_1 + 4x_2 = 0$, so $x_1 = x_2$; row 3: $x_3 = -x_1 - x_2 = -2x_1$: $\mathbf{x}_3 = (1, 1, -2)^\top$.

Verify one: $A(1, 1, -2)^\top = (3 + 1 + 2,\ 1 + 3 + 2,\ -1 - 1 - 10)^\top = (6, 6, -12)^\top = 6\,\mathbf{x}_3$. ✓

*Orthogonality.* $\mathbf{x}_1\cdot\mathbf{x}_2 = 1 - 1 + 0 = 0$; $\mathbf{x}_1\cdot\mathbf{x}_3 = 1 - 1 + 0 = 0$; $\mathbf{x}_2\cdot\mathbf{x}_3 = 1 + 1 - 2 = 0$. ✓ Three mutually perpendicular directions, as Theorem 3.17(b) demands. Normalize them to unit length and put them in $Q$: then $Q^\top Q = I$ and $A = Q\Lambda Q^\top$, which is LA-8's spectral theorem for this matrix.

*NumPy.* `np.linalg.eigh(A)` returns the eigenvalues as `[2., 3., 6.]` and eigenvectors $\pm(0.7071, -0.7071, 0)$, $(0.5774, 0.5774, 0.5774)$, $(-0.4082, -0.4082, 0.8165)$. The last is $-(1, 1, -2)/\sqrt6$. Not wrong: a different scale and sign for the same direction. See §8.

## 5. Why it matters and where it connects

- **LA-8 (spectral theorem).** Theorem 3.17 for distinct eigenvalues becomes $A = Q\Lambda Q^\top$ for every real symmetric $A$, repeated eigenvalues included. Positive definiteness is "all eigenvalues $> 0$." The Rayleigh quotient from §3.6 becomes a variational characterization of eigenvalues.
- **LA-9 (SVD).** The singular values of $A$ are the square roots of the eigenvalues of $A^\top A$; Example 4.6-style computations on $A^\top A$ are the first SVDs you compute by hand.
- **LA-11 (Jordan form).** What to do when §3.14 happens: generalized eigenvectors, $e^{At}$ for defective $A$, Cayley–Hamilton (C1 proves it for diagonalizable $A$), and functions of matrices in general.
- **LA-13 (numerical LA).** Power and inverse iteration (§3.6) are the seeds of the QR algorithm. The characteristic polynomial is never formed; §8 explains why. Conditioning of eigenvalues ties to how close $X$ is to singular (B7(d)).
- **CM-6 (normal modes).** Equation (3.15), $K\mathbf{v} = \omega^2M\mathbf{v}$, *is* the normal-mode problem. CM-6 derives $M$ and $K$ from the Lagrangian of any system near stable equilibrium and then uses this module's machinery, including the $M$-orthogonality of B6.
- **CM-7 (rigid bodies).** Principal axes are eigenvectors of the inertia tensor; Theorem 3.17 is why they are perpendicular.
- **F5 and ST-6 (Markov chains, MCMC).** Stationary distributions are the $\lambda = 1$ eigenvector of §3.7 (transposed to row convention). Mixing time is governed by $|\lambda_2|$, exactly as in Example 4.1. ST-7: an AR($p$) process is stationary iff the eigenvalues of its companion matrix lie inside the unit circle, which is (3.8).
- **QM-4 (formalism).** Observables are Hermitian matrices; their eigenvalues are the possible measurement outcomes (real, by the complex version of Theorem 3.17); the Schrödinger equation $i\hbar\dot\psi = H\psi$ is (3.10) with $A = -iH/\hbar$, solved by the matrix exponential $e^{-iHt/\hbar}$ of §3.9; eigenvectors of $H$ are stationary states. QC-1: single-qubit gates are $e^{-i\theta\,\sigma/2}$, computed exactly as in (3.20).
- **ML-2 (optimization).** Gradient descent on a quadratic is $\mathbf{u}_{k+1} = (I - \alpha H)\mathbf{u}_k$; its convergence rate is the largest $|1 - \alpha\lambda_i|$ over the eigenvalues of the Hessian, which is (3.8). ML-7/ML-8: exploding and vanishing gradients in recurrent networks are powers of a weight matrix, (3.8) again. ML-6: PCA is the eigendecomposition of a covariance matrix.
- **SM-6 (phase transitions).** The free energy per spin of the 1D Ising model, in the limit of infinitely many spins, is $-k_BT\log\lambda_1$ of a $2\times2$ transfer matrix: the dominant eigenvalue of §3.6.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 drill the $2\times2$ and triangular computations, a Markov chain, a rotation, and small exponentials. B1–B8 are the module: B2 (solve and classify an ODE system), B4 (Markov theory), B5 (complex eigenvalues to a real spiral), B6 (the generalized eigenproblem and $M$-orthogonality; the hardest), B7 (a defective matrix and what perturbing it does), and B8 (measure the convergence rates of power and inverse iteration with `project.py`, and explain the symmetric speed-up) are essential. C1 proves Cayley–Hamilton for diagonalizable matrices, C2 proves the Perron–Frobenius facts that follow from §3.7, and C3 finds the normal modes of an $n$-mass chain in closed form.

## 7. Coding project

`project.py` makes every object in this module compute:

1. **Unit vectors and their images** for the matrix of §2, with eigenvectors highlighted; then `power_iteration` and `inverse_iteration` (shifted), compared with `np.linalg.eig`, and a plot of the error against $k$ with the predicted slope $\log|\lambda_2/\lambda_1|$ from (3.9). The measured slopes are asserted: direction and nonsymmetric Rayleigh error at $\log|\lambda_2/\lambda_1|$, symmetric Rayleigh error at twice that, inverse iteration at $\log(|\lambda_{\text{closest}} - \sigma|/|\lambda_{\text{next}} - \sigma|)$.
2. **The weather chain** of Example 4.1: eigen-decomposition, the exact $\mathbf{u}_1, \mathbf{u}_2$, and a 20,000-step simulation whose empirical frequencies are compared with the $\lambda = 1$ eigenvector. Theorem 3.11 is checked on random Markov matrices.
3. **Fibonacci**: $A^k$ via (3.7), Binet's formula versus exact integer arithmetic, and where floating point gives out.
4. **$\dot{\mathbf{u}} = A\mathbf{u}$** solved by (3.11)/(3.13) and by `scipy.integrate.solve_ivp`, overlaid, for the spiral of Example 4.4 and the saddle of B2.
5. **The matrix exponential** by truncated series versus `scipy.linalg.expm`; $\frac{d}{dt}e^{At} = Ae^{At}$ by finite differences; the $e^Ae^B \ne e^{A+B}$ counterexample and the commuting case; Example 4.5's $e^{At} = e^{-t}(I + Nt)$ against `expm`; $e^{J\theta} = R(\theta)$.
6. **Phase portraits** of a saddle, a stable node, a stable spiral and a center, with eigenvectors drawn where they are real, plus the trace–determinant plane itself.
7. **Normal modes** of Example 4.3: frequencies and shapes from `scipy.linalg.eigh(K, M)`, $M$-orthogonality, the trajectory against `solve_ivp`, and each mode as a grid of frames.
8. **Identities**: trace and determinant on random matrices, eigenvalues of $A^2$, $A^{-1}$, $A + cI$, similar matrices, the shear's multiplicities, the nearly defective $A_\varepsilon$ of B7(d) (digits lost by $X\Lambda^kX^{-1}$ versus `matrix_power`), and Example 4.6's orthogonal eigenvectors.

**Build it yourself first:** before reading Part 1, write `power_iteration(A, k_max)` (multiply, normalize, repeat; return the Rayleigh quotient and the vector) and `inverse_iteration(A, shift, k_max)` (same loop with `np.linalg.solve(A - shift*I, v)` in place of the multiply). Before Part 4, write `solve_linear_ode(A, u0, t)` from (3.11): eigendecompose, find $\mathbf{c} = X^{-1}\mathbf{u}_0$, assemble $\sum c_ie^{\lambda_it}\mathbf{x}_i$, take the real part. Then compare.

**Your turn** (prompts at the bottom of the script): a three-state Markov chain of your own with its steady state and mixing rate; the normal modes of a three-mass chain against the closed form of C3; and $A^{100}$ computed by repeated multiplication and by $X\Lambda^{100}X^{-1}$, timed and compared.

## 8. Where you'll struggle

- **Complex eigenvalues.** They feel like a failure; they are an answer. Here is exactly what to do, every time: (1) compute **one** eigenvector $\mathbf{x}$ for **one** of the pair, $\lambda = a + ib$; (2) form $e^{\lambda t}\mathbf{x}$ and split it with Euler's formula into real and imaginary parts, (3.19); (3) those two real vector functions are your two independent real solutions, and the conjugate eigenvalue gives nothing new. Do not try to "make the eigenvector real"; it isn't. Do not add the conjugate solutions "to cancel the $i$"; taking the real part already did that. Example 4.4, then B5.
- **"My eigenvector is different from NumPy's."** Eigenvectors are defined up to a nonzero scalar (Theorem 3.2), and for complex eigenvalues up to a complex phase $e^{i\phi}$. NumPy returns unit-norm columns with an arbitrary sign (and phase). $(1, 1, -2)$ and $(-0.408, -0.408, 0.816)$ are the same eigenvector. Never compare eigenvectors entrywise. Either check $A\mathbf{v} \approx \lambda\mathbf{v}$ directly, or normalize both to the same convention (divide by the first nonzero entry), or check $|\mathbf{v}\cdot\mathbf{w}| = \|\mathbf{v}\|\|\mathbf{w}\|$. The project's `same_direction` helper does the last. Also: `np.linalg.eig` returns eigenvectors as **columns** of the second output, and the eigenvalue order is not sorted.
- **The characteristic polynomial is a $2\times2$ and $3\times3$ tool and nothing more.** For $n \ge 4$ you cannot factor it by hand, and in floating point you must not even form it: tiny errors in the coefficients of a polynomial produce large errors in its roots (Wilkinson's famous degree-20 example has roots that move by whole integers when one coefficient changes in its tenth digit). Real eigenvalue algorithms (LA-13: QR iteration, which `np.linalg.eig` uses) never compute $p_A$. For now: by hand up to $3\times3$, by `eig` beyond.
- **Defective matrices are real and you will meet them.** Whenever a damped system is critically damped, whenever an ODE's characteristic equation has a double root, the matrix is defective and the diagonalization formulas fail. Recognize it: a repeated eigenvalue with $\operatorname{rank}(A - \lambda I) > n - \operatorname{AM}(\lambda)$. Fix it: $e^{At}$ from the series, via $A = \lambda I + N$ with $N$ nilpotent (Example 4.5), or Jordan form (LA-11). Nearly defective matrices are worse: everything "works" and the answers silently lose digits, about half of them for a matrix $10^{-8}$ away from a defective one and all of them closer in (B7(d)).
- **Diagonalizable vs invertible.** Unrelated properties. Invertible asks whether $0$ is an eigenvalue. Diagonalizable asks whether there are $n$ independent eigenvectors. Keep the four examples from §3.5 in mind.
- **The Markov convention.** This module: columns sum to $1$, $\mathbf{u}_{k+1} = A\mathbf{u}_k$. Most probability books: rows sum to $1$, $\pi_{k+1} = \pi_kP$ with row vectors, $P = A^\top$. Both are correct; mixing them gives a matrix whose "steady state" is nonsense. Check which way the sums go before you compute anything.
- **Eigenvalues of sums and products.** In general the eigenvalues of $A + B$ are not sums $\lambda_i + \mu_j$, and those of $AB$ are not products. A common eigenbasis guarantees it (B1(c)); so does a common basis in which both are upper triangular (the second example in the §3.4 warning). Without some such structure, do not assume it. People get this wrong in ML-2 when reasoning about Hessians, and in QM-4 when adding non-commuting operators.
- **The hardest exercises.** B6 looks like bookkeeping and isn't: you have to see that the symmetric-matrix proof transfers to a *pair* $(K, M)$ only with the right inner product. C3's trig identity is a one-liner once you spot it and impossible until you do; the hint is the whole battle.
- **Losing the $\mathbf{x} \ne \mathbf{0}$ condition.** Every "proof" that every number is an eigenvalue, and every "eigenvector" that solves $(A - \lambda I)\mathbf{x} = \mathbf{0}$ trivially, comes from forgetting it.

## 9. Check yourself

1. Define eigenvalue and eigenvector. Why must $\mathbf{x} \ne \mathbf{0}$? Is $\lambda = 0$ allowed, and what does it mean?
2. Derive $\det(A - \lambda I) = 0$ from the definition in three lines.
3. Write the characteristic polynomial of a $2\times2$ matrix in terms of its trace and determinant. Why does this prove $\lambda_1 + \lambda_2 = T$ and $\lambda_1\lambda_2 = D$?
4. Prove that eigenvectors for two distinct eigenvalues are independent (the $k = 2$ case) in four lines.
5. Write $A^k$ and $e^{At}$ in terms of $X$ and $\Lambda$. Where in the derivation of $e^{At} = Xe^{\Lambda t}X^{-1}$ did you use $X^{-1}X = I$?
6. State and prove that $\lambda = 1$ is an eigenvalue of every Markov matrix. Why is $\mathbf{1}$ generally *not* its eigenvector?
7. $\dot{\mathbf{u}} = A\mathbf{u}$ with $A$ real $2\times2$, eigenvalues $-1 \pm 3i$, eigenvector $\mathbf{r} + i\mathbf{s}$ for $-1 + 3i$. Write the general real solution. What does the trajectory look like?
8. Under exactly what hypothesis is $e^{A+B} = e^Ae^B$? Give a $2\times2$ counterexample from memory.
9. What is the generalized eigenproblem for $M\ddot{\mathbf{x}} = -K\mathbf{x}$, and what physical object is an eigenvector?
10. Prove that a real symmetric matrix has real eigenvalues. Where did symmetry enter?
11. The matrix $\begin{pmatrix}2&1\\0&2\end{pmatrix}$: algebraic and geometric multiplicities of $2$? Diagonalizable? What is $e^{At}$?
12. A $2\times2$ real matrix has $T = -3$, $D = 5$. Classify the phase portrait without computing eigenvalues. Then $T = 1$, $D = -2$. Then $T = 0$, $D = 4$.

## 10. Pair with

- **Strang, *Introduction to Linear Algebra*, ch. 6 (§6.1–6.3).** The closest match to these notes in spirit and order; §6.2 on diagonalization and §6.3 on $e^{At}$ have the best collection of small worked examples anywhere. §6.4–6.5 (symmetric and positive definite) are LA-8. §10.3 has the Markov material with more probability.
- **MIT 18.06, lectures 21–25 (Strang, OCW).** Lecture 21 (eigenvalues), 22 (diagonalization and powers), 23 (differential equations and $e^{At}$), 24 (Markov matrices), 25 (symmetric matrices). Watch 23 after §3.9; Strang's blackboard derivation of $e^{At}$ is the one to copy.
- **Axler, *Linear Algebra Done Right*, ch. 5.** Eigenvalues *without determinants*: an operator on a complex space has an eigenvalue because $\mathbf{v}, T\mathbf{v}, \dots, T^n\mathbf{v}$ are dependent and the resulting polynomial factors. Read §5A–5B after §3.5 for a second, basis-free proof of Theorem 3.6 and for *invariant subspaces*, the concept this module keeps circling without naming. Axler delays Example 4.6-style computations and is better for *why* than for *how*.
- **3Blue1Brown, *Essence of Linear Algebra*, episode 14.** Eigenvectors as the directions that stay on their own span under the transformation; the clearest animation of §2. Fifteen minutes; watch before §3.
- **Strogatz, *Nonlinear Dynamics and Chaos*, §5.1–5.2.** The trace–determinant plane of §3.15 with pictures of every case and the "why" from a dynamics perspective. Optional, but the best treatment of what the four portraits look like.
