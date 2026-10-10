# LA-11: Jordan Form, Generalized Eigenvectors, the Minimal Polynomial, Cayley–Hamilton, and Functions of Matrices

**Time:** ~10 hours · **Prerequisites:** LA-7 (LA-10 is used in two proofs; see §0) · **Unlocks:** LA-12 (indirectly), CM-6 (the critically damped oscillator), QM-8 (degenerate perturbation theory intuition)

> **Conventions in this module.** Scalars are complex unless stated otherwise; $A$ is an $n\times n$ matrix, equivalently an operator on $V = \mathbb{C}^n$. $I$ is the identity, $\ker$ the null space, $\operatorname{im}$ the column space (range). The characteristic polynomial is the **monic** one, $p_A(x) = \det(xI - A) = (-1)^n\det(A - xI)$; its roots are the eigenvalues. $N_s$ is the $s\times s$ **shift**: ones on the superdiagonal, zeros elsewhere, so $N_s\mathbf{e}_1 = \mathbf{0}$ and $N_s\mathbf{e}_j = \mathbf{e}_{j-1}$. A **Jordan block** is $J_s(\lambda) = \lambda I_s + N_s$: $\lambda$ on the diagonal, ones on the **superdiagonal** (Axler, Horn & Johnson, Strang and SymPy all do this; some older books put the ones below the diagonal, which is the same thing with the basis reversed). A **chain** $\mathbf{v}_s \to \mathbf{v}_{s-1} \to \cdots \to \mathbf{v}_1 \to \mathbf{0}$ means $N\mathbf{v}_j = \mathbf{v}_{j-1}$ and $N\mathbf{v}_1 = \mathbf{0}$; as a basis it is always written **bottom first**, $(\mathbf{v}_1, \dots, \mathbf{v}_s)$, which is what puts the ones above the diagonal. For an eigenvalue $\lambda$: $d_k(\lambda) = \dim\ker(A - \lambda I)^k$, $a_\lambda$ is the algebraic multiplicity (multiplicity of $\lambda$ as a root of $p_A$), $g_\lambda = d_1(\lambda)$ the geometric multiplicity. Where an inner product appears (Proposition 3.12) it is $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^*\mathbf{v}$, conjugate-linear in the first slot, the physics convention of LA-8.

## 0. Where this sits

LA-7 ended with a gap. Diagonalization $A = PDP^{-1}$ needs $n$ independent eigenvectors, and you met matrices that don't have them: $\begin{pmatrix}1&1\\0&1\end{pmatrix}$ has eigenvalue $1$ twice and a single eigenvector. LA-7 called such matrices *defective* and moved on. This module closes the gap: every square matrix over $\mathbb{C}$ is similar to a matrix that is diagonal *plus a nilpotent shift*, the Jordan form, and that shift is exactly what produces the $t\,e^{\lambda t}$ solutions you have met in ODEs since your first course.

You need LA-7 cold: eigenvalues as roots of the characteristic polynomial, similarity $P^{-1}AP$, $e^{At}$ for diagonalizable $A$, and the fact that every complex matrix has at least one eigenvalue. You need rank–nullity (LA-4) and direct sums (LA-3). Schur's triangularization (LA-10) is used in two proofs; if you haven't done LA-10, all those proofs need is that every complex matrix is similar to an upper triangular one, and Exercise C1 proves that from scratch in half a page.

What it makes possible: CM-6 (the critically damped oscillator *is* a Jordan block; normal modes with repeated frequencies are *not*, and Proposition 3.12 says why), QM-8 (what a repeated eigenvalue does to eigenvectors under perturbation; the $\varepsilon$-demo of §3.11 is the intuition), ML-2 (heavy-ball momentum at its optimal setting is a Jordan block), ST-7 (companion matrices of AR processes are nonderogatory, Exercise B5), and LA-13 (the "never compute it in floating point" lesson, with the tools to say why).

## 1. The problem that forced this idea

Solve $y'' - 2y' + y = 0$. The characteristic equation $r^2 - 2r + 1 = (r-1)^2 = 0$ has one root, $r = 1$, so $e^t$ is a solution. A second-order equation needs two independent solutions. The other one is $t\,e^t$. You verified this in your ODE course by substitution and were told "for a repeated root, multiply by $t$." Euler had this rule in the 1740s. But *why*? Where does a polynomial factor $t$ come from, in a problem whose only ingredients are exponentials?

Turn it into a system, as LA-7 taught. With $x_1 = y$, $x_2 = y'$:
$$\mathbf{x}' = A\mathbf{x}, \qquad A = \begin{pmatrix} 0 & 1 \\ -1 & 2\end{pmatrix}.$$
LA-7's recipe: find eigenvalues, find eigenvectors, write $\mathbf{x}(t) = \sum_i c_i e^{\lambda_i t}\mathbf{v}_i$. The characteristic polynomial is $p_A(x) = x(x - 2) + 1 = (x-1)^2$, so $\lambda = 1$ twice. Eigenvectors: $A - I = \begin{pmatrix}-1&1\\-1&1\end{pmatrix}$ has rank $1$, so its kernel is a *line*, spanned by $(1,1)^\top$. One eigenvector. The recipe produces $c_1 e^t (1,1)^\top$ and stops. Half the solutions are missing. Diagonalization fails: there is no basis of eigenvectors, so there is no $P$ with $P^{-1}AP$ diagonal.

This is not a rare accident you can ignore. It happens in every critically damped oscillator (CM-6), in companion matrices of ODEs and recurrences with repeated roots (ST-7), in transition matrices with chains of transient states (F5), in momentum methods tuned to their fastest setting (ML-2), and in every matrix close enough to one of those that floating point can't tell the difference (§3.11).

The question "what is the simplest matrix similar to a given $A$, when diagonal is impossible?" was answered by Camille Jordan in his 1870 *Traité des substitutions et des équations algébriques*. Jordan was studying linear substitutions over finite fields for group theory, not ODEs, and the canonical form appears there as a tool. Weierstrass had obtained an equivalent result for complex matrices two years earlier, in his theory of elementary divisors; the name stuck to Jordan. Twelve years before either, Arthur Cayley's 1858 *Memoir on the Theory of Matrices* stated that every matrix satisfies its own characteristic equation, verified it for $2\times 2$ and $3\times 3$, and wrote that he had "not thought it necessary to undertake the labour of a formal proof of the theorem in the general case." Hamilton had the $3\times 3$ case in 1853 from quaternions. Frobenius gave the first general proof in 1878 and in the same paper introduced the minimal polynomial. Sylvester named nilpotent matrices in 1884.

The answer, in one sentence: **after a change of basis, every matrix is diagonal plus a shift, and the shift is what makes $t\,e^{\lambda t}$ appear.** The rest of this module makes each word precise, proves it, shows you how to compute with it, and tells you when you must not.

## 2. Intuition

### Stretch plus shift

Look at the $2\times 2$ Jordan block
$$J = \begin{pmatrix}\lambda & 1\\ 0 & \lambda\end{pmatrix} = \lambda I + N, \qquad N = \begin{pmatrix}0&1\\0&0\end{pmatrix}.$$
It does two things. The $\lambda I$ part stretches every vector by $\lambda$. The $N$ part is a **shift**: $N\mathbf{e}_2 = \mathbf{e}_1$ and $N\mathbf{e}_1 = \mathbf{0}$. It moves $\mathbf{e}_2$ onto $\mathbf{e}_1$ and kills $\mathbf{e}_1$. Apply it twice and everything is dead: $N^2 = 0$. An operator some power of which is zero is called **nilpotent**.

```
   e_2  --N-->  e_1  --N-->  0
```

So $J\mathbf{e}_1 = \lambda\mathbf{e}_1$ ($\mathbf{e}_1$ is an honest eigenvector) but $J\mathbf{e}_2 = \lambda\mathbf{e}_2 + \mathbf{e}_1$: $\mathbf{e}_2$ is stretched *and nudged* toward $\mathbf{e}_1$. It is "almost an eigenvector": $(J - \lambda I)\mathbf{e}_2 \ne \mathbf{0}$ but $(J-\lambda I)^2\mathbf{e}_2 = \mathbf{0}$. That is the whole content of the phrase *generalized eigenvector*.

A $3\times 3$ block is a longer chain,
```
   e_3  --N-->  e_2  --N-->  e_1  --N-->  0        N = J - λI,  N^3 = 0
```
and the general picture of a nilpotent operator is a **collection of chains** of various lengths, each ending in $\mathbf{0}$:
```
   v_3 --> v_2 --> v_1 --> 0          (a chain of length 3)
           w_2 --> w_1 --> 0          (a chain of length 2)
                   u_1 --> 0          (a chain of length 1: an eigenvector)
```
The bottom of each chain is a genuine eigenvector (it is in $\ker N$). The number of chains is $\dim\ker N$. A chain of length $s$ contributes one vector to $\ker N$, two to $\ker N^2$, and so on up to $s$: in the picture, $\ker N^k$ is spanned by the $k$ columns nearest the zero. So if you know $\dim\ker N^k$ for every $k$, you can count the chains of each length. That counting is the easy half of the module (Theorem 3.11). Proving that the chains always exist is the hard half (Theorem 3.6).

### Why a shift makes $t$ appear

Solve $\mathbf{y}' = J\mathbf{y}$ in chain coordinates for the $2\times 2$ block. The equations are
$$y_2' = \lambda y_2, \qquad y_1' = \lambda y_1 + y_2.$$
The top of the chain is on its own: $y_2 = c_2 e^{\lambda t}$. That solution is then *fed into* the equation for $y_1$ as a forcing term, and the forcing is at the natural rate of $y_1$, which is resonance: $y_1' - \lambda y_1 = c_2 e^{\lambda t}$ has the particular solution $c_2\, t\, e^{\lambda t}$. So $y_1 = (c_1 + c_2 t)e^{\lambda t}$. The $t$ is a resonance between a chain vector and the vector below it. A 3-chain feeds twice and gives $t^2 e^{\lambda t}/2$. You can guess the general formula now: $e^{Jt} = e^{\lambda t}(I + tN + \tfrac{t^2}{2}N^2 + \cdots)$, and the series stops because $N$ is nilpotent.

### What the theorem should say

Here is what you should be able to guess before §3. Take any $A$ on $\mathbb{C}^n$.

1. For each eigenvalue $\lambda$, instead of just $\ker(A - \lambda I)$ (the eigenvectors), take the bigger space $G_\lambda$ of vectors killed by *some* power of $A - \lambda I$. These are the generalized eigenvectors.
2. On $G_\lambda$, the operator $A - \lambda I$ is nilpotent, so it is a collection of chains.
3. The $G_\lambda$'s for different $\lambda$ fill the whole space with no overlap: $\mathbb{C}^n = G_{\lambda_1}\oplus\cdots\oplus G_{\lambda_m}$. (For a diagonalizable $A$ this is the eigenvector decomposition of LA-7; now it is true for every $A$.)
4. Put all the chains side by side as a basis. In that basis, $A$ is a block-diagonal matrix of Jordan blocks: **diagonal plus chains.**

Each step is where you'd expect a theorem to be, and each gets proved in §3. The thing to hold onto: an eigenvector is a chain of length 1; a Jordan block is a chain; diagonalizable means all chains have length 1.

### The one thing to be suspicious of

Diagonalizability is a *yes/no* property, and $\begin{pmatrix}\lambda&1\\ \varepsilon&\lambda\end{pmatrix}$ is diagonalizable for every $\varepsilon\ne0$ and not at $\varepsilon = 0$. So the Jordan form jumps discontinuously as the entries move. Floating point cannot see the difference between $\varepsilon = 0$ and $\varepsilon = 10^{-17}$. Keep that in mind throughout: this module is theory that explains ODEs, similarity, and matrix functions. It is not an algorithm you would run on data. §3.11 says what you run instead.

## 3. Formal development

### 3.1 Invariant subspaces

**Definition 3.1.** A subspace $U \subseteq V$ is **invariant** under $A$ if $A\mathbf{u} \in U$ for every $\mathbf{u}\in U$. Then $A$ restricts to an operator $A|_U$ on $U$.

Examples: $\ker A$ and $\operatorname{im} A$ are invariant (if $A\mathbf{u} = \mathbf{0}$ then $A(A\mathbf{u}) = \mathbf{0}$; if $\mathbf{u} = A\mathbf{w}$ then $A\mathbf{u} = A(A\mathbf{w})$). Every eigenspace is invariant. More generally, $\ker q(A)$ and $\operatorname{im} q(A)$ are invariant for every polynomial $q$, because $A$ commutes with $q(A)$: if $q(A)\mathbf{u} = \mathbf{0}$ then $q(A)(A\mathbf{u}) = A\,q(A)\mathbf{u} = \mathbf{0}$.

Why you care: if $V = U_1\oplus U_2$ with both $U_i$ invariant, and you pick a basis of $U_1$ followed by a basis of $U_2$, the matrix of $A$ in that basis is block diagonal, $\begin{pmatrix}A|_{U_1} & 0\\ 0 & A|_{U_2}\end{pmatrix}$ (LA-4: column $j$ of the matrix is the coordinates of $A$ applied to basis vector $j$). "Decompose $V$ into invariant pieces" is the same as "make $A$ block diagonal." The Jordan form is this with the best possible pieces.

### 3.2 Nilpotent operators

**Definition 3.2.** $N$ is **nilpotent** if $N^m = 0$ for some $m\ge1$. The smallest such $m$ is the **index** of $N$.

The shift $N_s$ is nilpotent of index $s$: $N_s^k\mathbf{e}_j = \mathbf{e}_{j-k}$ if $j > k$ and $\mathbf{0}$ otherwise, so $N_s^{s-1}\mathbf{e}_s = \mathbf{e}_1\neq\mathbf{0}$ but $N_s^s = 0$. Any strictly upper triangular matrix is nilpotent (Exercise A2). The matrix $\begin{pmatrix}1&1&-1\\1&1&-1\\2&2&-2\end{pmatrix} = \mathbf{u}\mathbf{v}^\top$ with $\mathbf{u} = (1,1,2)^\top$, $\mathbf{v} = (1,1,-1)^\top$, $\mathbf{v}^\top\mathbf{u} = 0$, is nilpotent of index 2 and looks nothing like a shift. A nilpotent operator has only the eigenvalue $0$: if $N\mathbf{v} = \lambda\mathbf{v}$ with $\mathbf{v}\ne\mathbf{0}$ then $\mathbf{0} = N^m\mathbf{v} = \lambda^m\mathbf{v}$ forces $\lambda = 0$.

**Lemma 3.3 (Kernels of powers).** For any operator $T$ on $V$ with $\dim V = n$:
(a) $\{\mathbf{0}\} = \ker T^0 \subseteq \ker T \subseteq \ker T^2 \subseteq \cdots$.
(b) If $\ker T^k = \ker T^{k+1}$ for some $k$, then $\ker T^{k} = \ker T^{k+j}$ for all $j \ge 0$: once the chain stops growing, it never grows again.
(c) $\ker T^n = \ker T^{n+1} = \cdots$.

▶ **Try it first.** For (b), take $\mathbf{v}\in\ker T^{k+2}$ and look at $T\mathbf{v}$.

**Derivation.**
1. (a) If $T^k\mathbf{v} = \mathbf{0}$ then $T^{k+1}\mathbf{v} = T(T^k\mathbf{v}) = \mathbf{0}$.
2. (b) Suppose $\ker T^k = \ker T^{k+1}$ and let $\mathbf{v}\in\ker T^{k+2}$. Then $T^{k+1}(T\mathbf{v}) = \mathbf{0}$, so $T\mathbf{v}\in\ker T^{k+1} = \ker T^k$, so $T^{k}(T\mathbf{v}) = T^{k+1}\mathbf{v} = \mathbf{0}$, so $\mathbf{v}\in\ker T^{k+1}$. Hence $\ker T^{k+2} \subseteq \ker T^{k+1}$, and with (a) they are equal. Repeat for $k+3, k+4, \dots$.
3. (c) The dimensions $\dim\ker T^k$ form a non-decreasing sequence of integers between $0$ and $n$. If the sequence strictly increased at each of the steps $0\to1\to\cdots\to n$ it would reach $n$ by step $n$ and could not grow further. If it ever fails to increase, (b) says it has stopped for good. Either way it is constant from $k = n$ on. ∎

■ **What just happened.** Kernels of powers grow, then freeze, and they freeze by step $n$ at the latest.

**Proposition 3.4.** If $N$ is nilpotent on an $n$-dimensional space with index $m$, then the kernels $\{\mathbf{0}\} \subsetneq \ker N \subsetneq \ker N^2 \subsetneq \cdots \subsetneq \ker N^m = V$ are **strictly** increasing up to $m$. Consequently $m\le n$ and $N^n = 0$.

*Proof.* If $\ker N^k = \ker N^{k+1}$ for some $k < m$, Lemma 3.3(b) gives $\ker N^k = \ker N^m = V$, so $N^k = 0$ with $k < m$, contradicting minimality of $m$. Each strict inclusion raises the dimension by at least $1$ starting from $0$, so $m \le \dim V = n$, and $N^n = N^{n-m}N^m = 0$. ∎

The sequence $d_k = \dim\ker N^k$, $k = 0, 1, \dots, m$, is the fingerprint of a nilpotent operator. Everything in §3.6 is read off it.

### 3.3 The chain theorem for nilpotent operators

**Definition 3.5.** A **chain** (Jordan chain) for $N$ of length $s$ is a list $\mathbf{v}_1, \dots, \mathbf{v}_s$ with $N\mathbf{v}_j = \mathbf{v}_{j-1}$ for $j \ge 2$ and $N\mathbf{v}_1 = \mathbf{0}$. Equivalently $\mathbf{v}_j = N^{s-j}\mathbf{v}_s$: the whole chain is generated by its top $\mathbf{v}_s$, and the condition on the top is $N^{s}\mathbf{v}_s = \mathbf{0} \ne N^{s-1}\mathbf{v}_s$.

In the basis $(\mathbf{v}_1, \dots, \mathbf{v}_s)$ (bottom first), the matrix of $N$ restricted to the span of the chain is the shift $N_s$: column $j$ records $N\mathbf{v}_j = \mathbf{v}_{j-1}$, a single $1$ in row $j-1$. This is why the ones sit on the superdiagonal.

**Theorem 3.6 (Jordan basis for a nilpotent operator).** Let $N$ be nilpotent on a finite-dimensional $V$. Then $V$ has a basis that is a union of chains for $N$. In that basis the matrix of $N$ is $N_{s_1}\oplus\cdots\oplus N_{s_r}$, a direct sum of shifts.

▶ **Try it first.** Induct on $\dim V$. The subspace $U = \operatorname{im}N$ is invariant and strictly smaller (why?), and $N|_U$ is nilpotent, so by induction $U$ has a chain basis. Each chain in $U$ has a top $\mathbf{u}$, and $\mathbf{u}\in\operatorname{im}N$ means $\mathbf{u} = N\mathbf{v}$ for some $\mathbf{v}$: the chain extends upward by one. What is still missing? Count dimensions with rank–nullity.

**Derivation.**
1. *Base case.* If $\dim V = 0$ there is nothing to do. If $N = 0$, any basis is a union of chains of length $1$. So assume $N\ne0$, $\dim V = n \ge 1$, and the theorem for all spaces of smaller dimension.
2. *The image is smaller.* $U = \operatorname{im}N$ is invariant (§3.1). Since $N$ is nilpotent, $\ker N \ne \{\mathbf{0}\}$: with $m$ the index, pick $\mathbf{w}$ with $N^{m-1}\mathbf{w}\ne\mathbf{0}$; then $N^{m-1}\mathbf{w}\in\ker N$. By rank–nullity (LA-4), $\dim U = n - \dim\ker N < n$. And $N|_U$ is nilpotent.
3. *Induction gives chains in $U$.* By the induction hypothesis, $U$ has a basis made of chains. Write them as $\mathbf{u}_{i,1}, \dots, \mathbf{u}_{i,s_i}$ for $i = 1,\dots,r$: chain $i$ has length $s_i$, bottom $\mathbf{u}_{i,1}\in\ker N$, top $\mathbf{u}_{i,s_i}$. In total $\sum_i s_i = \dim U$.
4. *Extend each chain upward by one.* Each top $\mathbf{u}_{i,s_i}$ lies in $U = \operatorname{im}N$, so choose $\mathbf{v}_i$ with $N\mathbf{v}_i = \mathbf{u}_{i,s_i}$. Now $\mathbf{u}_{i,1}, \dots, \mathbf{u}_{i,s_i}, \mathbf{v}_i$ is a chain of length $s_i + 1$.
5. *Fill out the kernel.* The bottoms $\mathbf{u}_{1,1}, \dots, \mathbf{u}_{r,1}$ lie in $\ker N$ and are independent (they are part of a basis of $U$). Extend them to a basis of $\ker N$ by adding vectors $\mathbf{w}_1, \dots, \mathbf{w}_p$, so $r + p = \dim\ker N$. Each $\mathbf{w}_k$ is a chain of length $1$.
6. *Count.* The proposed basis is $\{\mathbf{u}_{i,j}\} \cup \{\mathbf{v}_i\} \cup \{\mathbf{w}_k\}$, with $\dim U + r + p = \dim\operatorname{im}N + \dim\ker N = n$ vectors by rank–nullity. So it suffices to prove independence.
7. *Independence.* Suppose
$$\sum_{i,j} a_{ij}\mathbf{u}_{i,j} + \sum_i b_i\mathbf{v}_i + \sum_k c_k\mathbf{w}_k = \mathbf{0}.$$
Apply $N$. The $\mathbf{w}_k$ die. $N\mathbf{v}_i = \mathbf{u}_{i,s_i}$. $N\mathbf{u}_{i,j} = \mathbf{u}_{i,j-1}$ for $j\ge2$ and $\mathbf{0}$ for $j = 1$. So
$$\sum_{i}\sum_{j\ge2} a_{ij}\mathbf{u}_{i,j-1} + \sum_i b_i\mathbf{u}_{i,s_i} = \mathbf{0}.$$
This is a relation among the basis vectors of $U$. The top vector $\mathbf{u}_{i,s_i}$ of chain $i$ appears only in the second sum (the first sum only produces indices $j - 1 \le s_i - 1$), so every $b_i = 0$; then every $a_{ij}$ with $j\ge2$ is $0$.
8. What remains of the original relation is $\sum_i a_{i1}\mathbf{u}_{i,1} + \sum_k c_k\mathbf{w}_k = \mathbf{0}$, a relation among the basis of $\ker N$ chosen in step 5, so all $a_{i1} = 0$ and all $c_k = 0$. Independent, hence a basis. ∎

■ **What just happened.** Chains in the image pull back to longer chains in the whole space, and what the image misses is exactly the extra kernel directions, which are chains of length 1. Rank–nullity made the count come out. This is the hard half of the module. Read it twice, then do Example 4.3, which runs the same construction top-down as an algorithm.

### 3.4 Generalized eigenvectors and the decomposition theorem

**Definition 3.7.** For an eigenvalue $\lambda$ of $A$ on $V$ ($\dim V = n$), the **generalized eigenspace** is
$$G_\lambda = \ker(A - \lambda I)^n = \{\mathbf{v} : (A-\lambda I)^k\mathbf{v} = \mathbf{0}\text{ for some }k\ge1\}. \tag{3.2}$$
Its nonzero elements are **generalized eigenvectors**. The two descriptions agree by Lemma 3.3(c): if some power kills $\mathbf{v}$, the $n$-th does.

Three immediate facts. $G_\lambda \supseteq \ker(A - \lambda I) = E_\lambda$, the ordinary eigenspace. $G_\lambda$ is invariant under $A$ (§3.1: it is the kernel of a polynomial in $A$). And $(A - \lambda I)|_{G_\lambda}$ is nilpotent, by definition. So Theorem 3.6 applies to it, and that is where the Jordan blocks will come from. What remains is to show the $G_\lambda$'s fit together to fill $V$.

**Lemma 3.8.** If $\mu\ne\lambda$, then $A - \mu I$ is injective on $G_\lambda$. Consequently any polynomial $q(A)$ whose roots avoid $\lambda$ is injective on $G_\lambda$.

*Proof.* Let $\mathbf{v}\in G_\lambda$ with $(A - \mu I)\mathbf{v} = \mathbf{0}$, i.e. $A\mathbf{v} = \mu\mathbf{v}$. Then $(A-\lambda I)\mathbf{v} = (\mu - \lambda)\mathbf{v}$, and by induction $(A - \lambda I)^n\mathbf{v} = (\mu-\lambda)^n\mathbf{v}$. The left side is $\mathbf{0}$ because $\mathbf{v}\in G_\lambda$; since $\mu\ne\lambda$, $\mathbf{v} = \mathbf{0}$. For the consequence, factor $q(A) = c\prod_i (A - \mu_i I)$ over $\mathbb{C}$ with every $\mu_i \ne \lambda$; each factor is injective on $G_\lambda$ and maps $G_\lambda$ into itself (invariance), so the composition is injective on $G_\lambda$. ∎

**Theorem 3.9 (Generalized eigenspace decomposition).** Let $A$ act on $V = \mathbb{C}^n$ with distinct eigenvalues $\lambda_1, \dots, \lambda_m$ and algebraic multiplicities $a_1, \dots, a_m$, so $p_A(x) = \prod_i (x-\lambda_i)^{a_i}$ and $\sum_i a_i = n$. Then:
(a) each $G_{\lambda_i}$ is invariant under $A$;
(b) the $G_{\lambda_i}$ are independent: if $\mathbf{v}_1 + \cdots + \mathbf{v}_m = \mathbf{0}$ with $\mathbf{v}_i\in G_{\lambda_i}$, then every $\mathbf{v}_i = \mathbf{0}$;
(c) $\dim G_{\lambda_i} = a_i$;
and therefore $V = G_{\lambda_1}\oplus\cdots\oplus G_{\lambda_m}$.

▶ **Try it first.** For (b), hit the relation with the operator $Q_j = \prod_{i\ne j}(A - \lambda_i I)^n$, which kills every $G_{\lambda_i}$ with $i\ne j$. For (c), use an upper triangular form of $A$ (Schur, LA-10) with the $a_i$ copies of $\lambda_i$ placed first on the diagonal, and look at the kernel of $(T - \lambda_iI)^n$ block by block.

**Derivation.**
1. (a) was noted after Definition 3.7.
2. (b) Fix $j$ and apply $Q_j = \prod_{i\ne j}(A-\lambda_iI)^n$ to $\sum_i \mathbf{v}_i = \mathbf{0}$. For $i\ne j$, $Q_j$ contains the factor $(A-\lambda_iI)^n$, and the factors commute, so $Q_j\mathbf{v}_i = \mathbf{0}$. What survives is $Q_j\mathbf{v}_j = \mathbf{0}$. By Lemma 3.8, $Q_j$ is injective on $G_{\lambda_j}$ (its roots are the $\lambda_i$ with $i\neq j$, all different from $\lambda_j$). So $\mathbf{v}_j = \mathbf{0}$. Since $j$ was arbitrary, (b) holds.
3. (c) Schur's theorem (LA-10): $A = QTQ^*$ with $Q$ unitary and $T$ upper triangular. The construction starts from any eigenvector you like and recurses on the complement, so the eigenvalues can be placed on the diagonal of $T$ in any order you choose (Exercise C1 proves this directly, without inner products; unitarity is not needed here). Put the $a_i$ copies of $\lambda_i$ first. Then
$$T - \lambda_i I = \begin{pmatrix} N & B \\ 0 & R\end{pmatrix},$$
with $N$ an $a_i\times a_i$ strictly upper triangular block and $R$ upper triangular with diagonal entries $\lambda_j - \lambda_i \ne 0$, hence invertible (a triangular matrix is invertible iff its diagonal is nonzero, LA-6).
4. Powers of a block upper triangular matrix are block upper triangular with the diagonal blocks raised to the same power (multiply out two such matrices and look at the diagonal blocks): $(T - \lambda_iI)^n = \begin{pmatrix} N^n & X \\ 0 & R^n\end{pmatrix}$. $N$ is strictly upper triangular, so $N\mathbf{e}_k\in\operatorname{span}(\mathbf{e}_1,\dots,\mathbf{e}_{k-1})$ and hence $N^k\mathbf{e}_k = \mathbf{0}$ for every $k\le a_i$; so $N^{a_i} = 0$ and $N^n = 0$. $R^n$ is invertible.
5. Solve $\begin{pmatrix}0 & X\\0 & R^n\end{pmatrix}\begin{pmatrix}\mathbf{p}\\ \mathbf{q}\end{pmatrix} = \mathbf{0}$: the bottom block gives $R^n\mathbf{q} = \mathbf{0}$, so $\mathbf{q} = \mathbf{0}$, and then the top block is $X\mathbf{q} = \mathbf{0}$ automatically, so $\mathbf{p}$ is free. The kernel is $\{(\mathbf{p}, \mathbf{0})\}$, of dimension exactly $a_i$. Kernel dimension is a similarity invariant ($\ker(A-\lambda_iI)^n = Q\,\ker(T-\lambda_iI)^n$), so $\dim G_{\lambda_i} = a_i$.
6. By (b) the sum $G_{\lambda_1} + \cdots + G_{\lambda_m}$ is direct, so its dimension is $\sum_i a_i = n = \dim V$. It is all of $V$. ∎

■ **What just happened.** Independence came from a polynomial that kills all generalized eigenspaces but one and is invertible on that one. The dimension count came from the triangular form, where the generalized eigenspace for $\lambda_i$ is visibly the first $a_i$ coordinates. Together: $V$ splits into invariant pieces, one per eigenvalue, each of the "right" size, with $A - \lambda_i I$ nilpotent on the $i$-th piece.

A consequence worth isolating, because LA-7 only asserted it: **algebraic multiplicity equals $\dim G_\lambda$, and the geometric multiplicity satisfies $g_\lambda = \dim E_\lambda \le a_\lambda$** (since $E_\lambda\subseteq G_\lambda$). Equality for every $\lambda$ is exactly diagonalizability.

### 3.5 The Jordan canonical form

**Theorem 3.10 (Jordan canonical form).** Every $A\in\mathbb{C}^{n\times n}$ is similar to a block-diagonal matrix of Jordan blocks:
$$P^{-1}AP = J = J_{s_1}(\lambda_1)\oplus J_{s_2}(\lambda_2)\oplus\cdots\oplus J_{s_r}(\lambda_r), \tag{3.3}$$
where the $\lambda_j$ are eigenvalues of $A$ (repetitions allowed) and $s_1 + \cdots + s_r = n$. The multiset of pairs $(\lambda_j, s_j)$ is uniquely determined by $A$; only the order of the blocks can change. The columns of $P$ form a **Jordan basis**: a union of chains, each written bottom first.

**Derivation.**
1. *Existence.* By Theorem 3.9, $V = \bigoplus_\lambda G_\lambda$ with each $G_\lambda$ invariant. On $G_\lambda$ the operator $N_\lambda = (A-\lambda I)|_{G_\lambda}$ is nilpotent. Theorem 3.6 gives a basis of $G_\lambda$ made of chains for $N_\lambda$. In that basis $N_\lambda$ is $N_{s_1}\oplus\cdots$, so $A|_{G_\lambda} = \lambda I + N_\lambda$ is $J_{s_1}(\lambda)\oplus\cdots$. Concatenating the bases of all the $G_\lambda$ gives a basis of $V$ (direct sum) in which $A$ is block diagonal with these blocks. The change-of-basis matrix $P$ has these basis vectors as columns (LA-4), and $AP = PJ$ reads column by column as $A\mathbf{v}_1 = \lambda\mathbf{v}_1$, $A\mathbf{v}_j = \lambda\mathbf{v}_j + \mathbf{v}_{j-1}$: the chain relations.
2. *Uniqueness* is Theorem 3.11 below: the block sizes are computable from the numbers $\dim\ker(A-\lambda I)^k$, which depend only on $A$. ∎

Note what is unique and what is not. The *block structure* is unique. The *basis* $P$ is wildly non-unique (Example 4.3 finds a different $P$ from the one the matrix was built with). The diagonalizable case is $J$ diagonal, i.e. all $s_j = 1$, i.e. $P$ a basis of eigenvectors: LA-7's $A = PDP^{-1}$.

### 3.6 Reading the block sizes off kernel dimensions

**Theorem 3.11 (Block counting).** Let $d_k = \dim\ker(A - \lambda I)^k$ for $k = 0, 1, 2, \dots$ (so $d_0 = 0$ and $d_1 = g_\lambda$). Then
$$\#\{\text{Jordan blocks for }\lambda\text{ of size}\ge k\} = d_k - d_{k-1}, \qquad k \ge 1. \tag{3.4}$$
In particular the number of blocks for $\lambda$ is $g_\lambda = d_1$, the largest block has size equal to the smallest $k$ with $d_k = a_\lambda$, and the number of blocks of size exactly $k$ is $(d_k - d_{k-1}) - (d_{k+1} - d_k)$.

▶ **Try it first.** Compute $\dim\ker N_s^k$ for a single shift. Then ask what the blocks for *other* eigenvalues contribute to $\ker(J - \lambda I)^k$.

**Derivation.**
1. $\dim\ker$ is a similarity invariant, so compute with $J$ in place of $A$. $(J - \lambda I)^k$ is block diagonal with blocks $(J_{s}(\mu) - \lambda I)^k$, and the kernel of a block-diagonal matrix is the direct sum of the kernels of the blocks.
2. A block with $\mu\ne\lambda$: $J_s(\mu) - \lambda I$ is triangular with diagonal entries $\mu - \lambda \ne 0$, so invertible, so its $k$-th power is invertible. Contribution to the kernel: $0$.
3. A block with $\mu = \lambda$: $J_s(\lambda) - \lambda I = N_s$ and $N_s^k\mathbf{e}_j = \mathbf{e}_{j-k}$ for $j > k$, $\mathbf{0}$ for $j\le k$. So $\ker N_s^k = \operatorname{span}(\mathbf{e}_1, \dots, \mathbf{e}_{\min(s,k)})$ has dimension $\min(s, k)$.
4. Summing over the blocks for $\lambda$ with sizes $s_1, \dots, s_r$: $d_k = \sum_i \min(s_i, k)$.
5. Difference: $d_k - d_{k-1} = \sum_i [\min(s_i,k) - \min(s_i, k-1)]$. The bracket is $1$ if $s_i\ge k$ and $0$ if $s_i \le k-1$. So the difference counts blocks of size $\ge k$. ∎

■ **What just happened.** A chain of length $s$ contributes one new dimension to $\ker N^k$ for each $k\le s$ and then stops. The jumps $d_k - d_{k-1}$ therefore count how many chains are still "alive" at level $k$. This is the formula to learn cold: it needs only ranks, no eigenvectors, and it is how you *check* any Jordan form someone hands you.

**The algorithm for small matrices (by hand or in exact arithmetic).** For each eigenvalue $\lambda$:
1. Compute ranks of $(A-\lambda I), (A-\lambda I)^2, \dots$ until the rank stops dropping; $d_k = n - \operatorname{rank}$. The block sizes follow from (3.4).
2. To get the basis, work top-down. Let $m$ be the largest block size. Pick vectors in $\ker(A-\lambda I)^m$ that are independent modulo $\ker(A-\lambda I)^{m-1}$ (you need $d_m - d_{m-1}$ of them); each is the top of an $m$-chain; generate the chains by applying $A - \lambda I$ repeatedly. Drop to level $m-1$: you already hold the vectors $(A-\lambda I)\cdot(\text{tops})$ at that level; add new vectors from $\ker(A-\lambda I)^{m-1}$ until, together with $\ker(A-\lambda I)^{m-2}$, they span $\ker(A-\lambda I)^{m-1}$; each new one tops an $(m-1)$-chain. Continue down to level $1$, where the new vectors are eigenvectors that are chains of length $1$.
3. Write every chain bottom first, side by side, as the columns of $P$.

Theorem 3.6 guarantees the counts in step 2 always work out; `nilpotent_chains` in the project does exactly this and asserts the counts at every level. Example 4.3 runs it by hand.

**Proposition 3.12 (No blocks for symmetric, Hermitian, or normal matrices).** If $A$ is real symmetric (or complex Hermitian, or normal, $AA^* = A^*A$), then $\ker(A-\lambda I)^2 = \ker(A - \lambda I)$ for every $\lambda$, so every Jordan block has size $1$ and $A$ is diagonalizable.

*Proof.* Let $S = A - \lambda I$. If $A$ is normal so is $S$, and for normal $S$ and any $\mathbf{w}$, $\|S\mathbf{w}\|^2 = \langle S\mathbf{w}, S\mathbf{w}\rangle = \langle S^*S\mathbf{w}, \mathbf{w}\rangle = \langle SS^*\mathbf{w}, \mathbf{w}\rangle = \|S^*\mathbf{w}\|^2$, so $\ker S = \ker S^*$. Now suppose $S^2\mathbf{v} = \mathbf{0}$ and set $\mathbf{w} = S\mathbf{v}$. Then $S\mathbf{w} = \mathbf{0}$, so $\mathbf{w}\in\ker S = \ker S^*$, so $S^*S\mathbf{v} = S^*\mathbf{w} = \mathbf{0}$, so $\|S\mathbf{v}\|^2 = \langle S^*S\mathbf{v}, \mathbf{v}\rangle = 0$, so $S\mathbf{v} = \mathbf{0}$. Hence $d_2 = d_1$, and by Lemma 3.3(b) every $d_k = d_1$; by (3.4) no block has size $\ge 2$. (Symmetric and Hermitian matrices are normal, and for real symmetric $S$ the argument is one line: $\|S\mathbf{v}\|^2 = \mathbf{v}^\top S^\top S\mathbf{v} = \mathbf{v}^\top S^2\mathbf{v} = 0$.) ∎

Why it matters: CM-6's normal modes (symmetric $K$ and $M$) and QM's Hermitian Hamiltonians never have Jordan blocks, even with repeated eigenvalues; the ring of three masses in the project has frequency $\sqrt3$ twice, a two-dimensional eigenspace, and no chain. Jordan blocks live in *non-symmetric* problems: damping, friction, non-reciprocal couplings, companion matrices.

### 3.7 The minimal polynomial

**Definition 3.13.** The **minimal polynomial** $m_A$ of $A$ is the monic polynomial of least degree with $m_A(A) = 0$.

**Theorem 3.14 (Properties of the minimal polynomial).**
(a) $m_A$ exists and is unique, with degree $\le n^2$ (in fact $\le n$, once Cayley–Hamilton is proved in §3.8).
(b) If $q(A) = 0$ for a polynomial $q$, then $m_A$ divides $q$.
(c) The roots of $m_A$ are exactly the eigenvalues of $A$.

**Derivation.**
1. *(a) Existence.* The $n^2 + 1$ matrices $I, A, A^2, \dots, A^{n^2}$ live in the $n^2$-dimensional space of $n\times n$ matrices, so they are linearly dependent (LA-3): $\sum_{k=0}^{n^2} c_k A^k = 0$ with not all $c_k = 0$. That is a nonzero polynomial annihilating $A$. Among all nonzero annihilating polynomials pick one of least degree and divide by its leading coefficient.
2. *Uniqueness.* If $m$ and $\tilde m$ are both monic annihilators of least degree $d$, then $m - \tilde m$ annihilates $A$ and has degree $< d$, so it is the zero polynomial.
3. *(b)* Divide: $q = s\,m_A + r$ with $\deg r < \deg m_A$ (polynomial long division). Then $r(A) = q(A) - s(A)m_A(A) = 0$. By minimality of $\deg m_A$, $r = 0$.
4. *(c), eigenvalues are roots.* If $A\mathbf{v} = \lambda\mathbf{v}$ with $\mathbf{v}\ne\mathbf{0}$, then $A^k\mathbf{v} = \lambda^k\mathbf{v}$ and so $\mathbf{0} = m_A(A)\mathbf{v} = m_A(\lambda)\mathbf{v}$, forcing $m_A(\lambda) = 0$.
5. *(c), roots are eigenvalues.* If $m_A(\mu) = 0$, factor $m_A(x) = (x - \mu)q(x)$. Since $\deg q < \deg m_A$, $q(A)\ne0$; pick $\mathbf{w}$ with $\mathbf{u} = q(A)\mathbf{w} \ne \mathbf{0}$. Then $(A-\mu I)\mathbf{u} = m_A(A)\mathbf{w} = \mathbf{0}$: $\mathbf{u}$ is an eigenvector with eigenvalue $\mu$. ∎

**Theorem 3.15 (Minimal polynomial from the Jordan form).** Let $s_{\max}(\lambda)$ be the size of the largest Jordan block for $\lambda$. Then
$$m_A(x) = \prod_{\lambda}(x-\lambda)^{s_{\max}(\lambda)}. \tag{3.5}$$
Consequently $A$ is diagonalizable if and only if $m_A$ has distinct roots (no repeated factor).

▶ **Try it first.** When does a polynomial $q(x) = \prod_\mu (x - \mu)^{e_\mu}$ kill a single block $J_s(\lambda)$? Write $q(J_s(\lambda))$ as a product and ask which factor can be zero.

**Derivation.**
1. $q(A) = P\,q(J)\,P^{-1}$ (because $A^k = PJ^kP^{-1}$), so $q(A) = 0$ iff $q(J) = 0$ iff $q$ kills every block.
2. On one block $J_s(\lambda)$, $q(J_s(\lambda)) = \prod_\mu (J_s(\lambda) - \mu I)^{e_\mu}$. The factors with $\mu\ne\lambda$ are invertible (triangular with nonzero diagonal). The factor with $\mu = \lambda$ is $N_s^{e_\lambda}$. A product of invertible matrices with $N_s^{e_\lambda}$ is zero iff $N_s^{e_\lambda} = 0$ iff $e_\lambda \ge s$.
3. So $q(A) = 0$ iff $e_\lambda\ge s$ for every block $J_s(\lambda)$, i.e. iff $e_\lambda \ge s_{\max}(\lambda)$ for every eigenvalue $\lambda$. The annihilating polynomial of least degree takes $e_\lambda = s_{\max}(\lambda)$ exactly, which is (3.5).
4. Diagonalizable iff all blocks have size $1$ iff all $s_{\max}(\lambda) = 1$ iff every root of $m_A$ is simple. ∎

■ **What just happened.** The characteristic polynomial records $\dim G_\lambda$ (total size of the blocks for $\lambda$); the minimal polynomial records the *largest* block. They agree exactly when each eigenvalue has a single block. Example 4.5 shows a matrix where they differ. Exercise B4 proves "distinct roots $\Rightarrow$ diagonalizable" without the Jordan form, which is worth seeing because that proof is the one that survives in infinite dimensions.

### 3.8 Cayley–Hamilton

**Theorem 3.16 (Cayley–Hamilton).** Every matrix satisfies its own characteristic polynomial:
$$p_A(A) = 0. \tag{3.6}$$

**The fake proof, so you recognize it.** "$p_A(x) = \det(xI - A)$, so $p_A(A) = \det(AI - A) = \det(0) = 0$." This is wrong, and not in a fixable-by-care way. $\det(xI - A)$ is a *number* for each scalar $x$, computed by expanding a determinant whose entries are $x - a_{ii}$ and $-a_{ij}$. Substituting a matrix for $x$ inside that determinant makes no sense: $xI - A$ would have matrix entries, and "$AI - A = 0$" produces the scalar $\det 0 = 0$, whereas $p_A(A)$ is an $n\times n$ matrix. The argument confuses the polynomial $p_A$, which can be evaluated at a matrix, with the determinant function that happened to produce it. (There *is* a correct proof that starts from the adjugate identity $\operatorname{adj}(xI - A)\,(xI - A) = p_A(x)I$ and matches coefficients of matrix-valued polynomials; Exercise C3 walks through it. It is not the fake proof, and the difference is the whole point.)

▶ **Try it first.** Take $T$ upper triangular with diagonal $\lambda_1, \dots, \lambda_n$ and $V_k = \operatorname{span}(\mathbf{e}_1, \dots, \mathbf{e}_k)$. Show $(T - \lambda_k I)V_k \subseteq V_{k-1}$. Then chain these together.

**Derivation (via triangularization).**
1. By Schur (LA-10), $A = QTQ^*$ with $Q$ unitary and $T$ upper triangular (Exercise C1: any similarity to upper triangular form will do). The diagonal of $T$ lists the eigenvalues $\lambda_1, \dots, \lambda_n$ with multiplicity, so $p_A(x) = p_T(x) = \prod_{j=1}^n (x - \lambda_j)$ (determinant of a triangular matrix, LA-6). Since $p_A(A) = Q\,p_A(T)\,Q^*$, it suffices to show $p_A(T) = (T-\lambda_1I)(T-\lambda_2I)\cdots(T-\lambda_nI) = 0$.
2. Let $V_k = \operatorname{span}(\mathbf{e}_1, \dots, \mathbf{e}_k)$, $V_0 = \{\mathbf{0}\}$. Upper triangular means $T\mathbf{e}_k\in V_k$ (column $k$ has zeros below row $k$), so each $V_k$ is invariant under $T$ and under every $T - \lambda I$.
3. Moreover $(T - \lambda_kI)\mathbf{e}_k = T\mathbf{e}_k - t_{kk}\mathbf{e}_k\in V_{k-1}$: the diagonal entry cancels, leaving only the entries above it. Combined with step 2 (which handles $V_{k-1}\subseteq V_k$), $(T-\lambda_kI)V_k \subseteq V_{k-1}$.
4. Define $M_k = (T-\lambda_1I)(T-\lambda_2I)\cdots(T-\lambda_kI)$. Claim: $M_kV_k = \{\mathbf{0}\}$. For $k = 1$: $(T-\lambda_1I)\mathbf{e}_1 = (t_{11}-\lambda_1)\mathbf{e}_1 = \mathbf{0}$. Inductive step: $M_k = M_{k-1}(T-\lambda_kI)$; by step 3, $(T-\lambda_kI)V_k\subseteq V_{k-1}$, and by induction $M_{k-1}V_{k-1} = \{\mathbf{0}\}$. So $M_kV_k = \{\mathbf{0}\}$.
5. At $k = n$, $V_n = \mathbb{C}^n$, so $M_n = p_A(T) = 0$. ∎

■ **What just happened.** Each factor $T - \lambda_kI$ pushes the $k$-th basis vector down one level of the flag $V_0\subset V_1\subset\cdots\subset V_n$, and the earlier factors finish the job. The product kills a bigger subspace at each step. The project checks exactly this: $M_k$ restricted to the first $k$ columns is zero to rounding for the Schur form of a random matrix. The factors commute, so the order in the product is immaterial; the proof chose the convenient one.

**Using it.** Cayley–Hamilton says $A^n$ is a linear combination of $I, A, \dots, A^{n-1}$, with coefficients read off $p_A$. Two consequences you will use:

- *Inverse as a polynomial.* Write $p_A(x) = x^n + c_{n-1}x^{n-1} + \cdots + c_1x + c_0$, where $c_0 = p_A(0) = (-1)^n\det A$. Then $A(A^{n-1} + c_{n-1}A^{n-2} + \cdots + c_1I) = -c_0I$, so if $\det A\ne0$,
$$A^{-1} = -\frac{1}{c_0}\left(A^{n-1} + c_{n-1}A^{n-2} + \cdots + c_1I\right).$$
- *Powers reduced.* To compute $A^k$ for large $k$, divide $x^k$ by $p_A(x)$: $x^k = s(x)p_A(x) + r(x)$ with $\deg r < n$, and then $A^k = r(A)$. Every power, and in fact every function of $A$ (§3.9 and Exercise C4), is a polynomial in $A$ of degree $< n$. Example 4.4 does both. Corollary: $m_A$ divides $p_A$ (Theorem 3.14(b)), so $\deg m_A \le n$.

### 3.9 Functions of matrices

For a polynomial $f$, $f(A)$ means substitute $A$ for $x$. For a power series $f(x) = \sum_k c_kx^k$ converging on a disk containing all eigenvalues of $A$, $f(A) = \sum_k c_kA^k$ converges entrywise; this is the definition LA-7 used for $e^{A}$. The Jordan form turns both into a finite computation.

**Theorem 3.17 (A function of a Jordan block).** For $J = J_s(\lambda) = \lambda I + N$ (so $N = N_s$, $N^s = 0$) and $f$ a polynomial, or a power series converging on a disk around $\lambda$,
$$f(J) = \sum_{k=0}^{s-1}\frac{f^{(k)}(\lambda)}{k!}N^k = \begin{pmatrix} f(\lambda) & f'(\lambda) & \frac{f''(\lambda)}{2!} & \cdots & \frac{f^{(s-1)}(\lambda)}{(s-1)!}\\ & f(\lambda) & f'(\lambda) & \ddots & \vdots\\ && \ddots & \ddots & \frac{f''(\lambda)}{2!}\\ &&& f(\lambda) & f'(\lambda)\\ &&&& f(\lambda)\end{pmatrix}. \tag{3.7}$$

▶ **Try it first.** Taylor-expand a polynomial $f$ about $\lambda$: $f(x) = \sum_k \frac{f^{(k)}(\lambda)}{k!}(x-\lambda)^k$, exactly, with finitely many terms. Substitute $x = J$. What is $(J - \lambda I)^k$?

**Derivation.**
1. For a polynomial $f$ of degree $D$, Taylor's formula about $\lambda$ is an exact finite identity between polynomials: $f(x) = \sum_{k=0}^{D}\frac{f^{(k)}(\lambda)}{k!}(x-\lambda)^k$. A polynomial identity survives substitution of a matrix for $x$ (both sides are built from sums and products of $x$ and constants, and $J$ commutes with $I$), so $f(J) = \sum_{k=0}^D\frac{f^{(k)}(\lambda)}{k!}(J-\lambda I)^k$.
2. $(J-\lambda I)^k = N^k$, and $N^k = 0$ for $k\ge s$. The sum terminates at $k = s-1$, which is (3.7). The matrix picture: $N^k$ has ones on the $k$-th superdiagonal, so $f^{(k)}(\lambda)/k!$ sits on the $k$-th superdiagonal.
3. For a power series $f$, apply steps 1–2 to the partial sums $f_D$, which are polynomials: $f_D(J) = \sum_{k<s}\frac{f_D^{(k)}(\lambda)}{k!}N^k$. As $D\to\infty$, $f_D(J)\to f(J)$ by definition, and $f_D^{(k)}(\lambda)\to f^{(k)}(\lambda)$ because a power series may be differentiated term by term inside its disk of convergence. Finitely many terms, each converging, so (3.7) holds for $f$. ∎

■ **What just happened.** On an eigenvector, $f(A)$ acts as the *number* $f(\lambda)$. On a chain, $f(A)$ acts as the *Taylor polynomial* of $f$ at $\lambda$, with the shift $N$ playing the role of $(x - \lambda)$. Derivatives of $f$ enter only through Jordan blocks. This is why repeated eigenvalues bring $f'$ into the picture, and why $\frac{\partial}{\partial\lambda}e^{\lambda t} = te^{\lambda t}$ shows up in ODEs.

**The exponential.** Take $f(x) = e^{tx}$, so $f^{(k)}(\lambda) = t^ke^{\lambda t}$:
$$e^{tJ_2(\lambda)} = e^{\lambda t}\begin{pmatrix}1 & t\\0 & 1\end{pmatrix}, \tag{3.8}$$
$$e^{tJ_3(\lambda)} = e^{\lambda t}\begin{pmatrix}1 & t & t^2/2\\ 0 & 1 & t\\ 0 & 0 & 1\end{pmatrix}, \tag{3.9}$$
and in general $e^{tJ_s(\lambda)} = e^{\lambda t}\sum_{k<s}\frac{t^k}{k!}N^k$. Equivalently $e^{t(\lambda I + N)} = e^{\lambda tI}e^{tN}$ because $\lambda I$ and $N$ commute (LA-7: $e^{X+Y} = e^Xe^Y$ when $XY = YX$), and $e^{tN} = \sum_{k<s}t^kN^k/k!$ terminates.

**Any matrix.** With $A = PJP^{-1}$ and $J = \bigoplus_j J_{s_j}(\lambda_j)$,
$$f(A) = P\,f(J)\,P^{-1} = P\Big(\bigoplus_j f\big(J_{s_j}(\lambda_j)\big)\Big)P^{-1}, \tag{3.10}$$
since $A^k = PJ^kP^{-1}$ and $J^k$ is block diagonal with blocks $J_{s_j}(\lambda_j)^k$. The left side does not depend on $P$ (it is defined by the series), so any Jordan basis gives the same answer.

**The ODE, explained.** $\mathbf{x}' = A\mathbf{x}$ has solution $\mathbf{x}(t) = e^{At}\mathbf{x}_0 = Pe^{Jt}P^{-1}\mathbf{x}_0$. Each block contributes $e^{\lambda t}$ times polynomials in $t$ of degree $< s$. So: **the solutions of a constant-coefficient linear system are combinations of $t^ke^{\lambda t}$ with $0\le k < s_{\max}(\lambda)$, and nothing else.** The $t$ comes from the shift. Diagonalizable $A$: no shifts, no polynomials, pure exponentials, which is what LA-7 promised. Examples 4.1 and 4.6 do this for the system of §1 and for CM-6's critically damped oscillator.

### 3.10 Similarity: the complete answer

**Theorem 3.18 (Similarity classification).** Two matrices $A, B\in\mathbb{C}^{n\times n}$ are similar if and only if they have the same Jordan form (same multiset of blocks $(\lambda, s)$); equivalently, iff $\dim\ker(A-\lambda I)^k = \dim\ker(B-\lambda I)^k$ for every $\lambda\in\mathbb{C}$ and every $k\ge1$.

*Proof.* If $A$ and $B$ have the same Jordan form $J$, then $A = PJP^{-1}$ and $B = RJR^{-1}$, so $B = (RP^{-1})A(RP^{-1})^{-1}$. Conversely, if $B = SAS^{-1}$, then $(B-\lambda I)^k = S(A-\lambda I)^kS^{-1}$, so the kernel dimensions agree for every $\lambda$ and $k$, and by Theorem 3.11 those dimensions determine the block structure. ∎

So the Jordan form is a *complete* similarity invariant. The familiar invariants are not: trace, determinant, rank, characteristic polynomial, even the minimal polynomial together with the characteristic polynomial. The project exhibits two $7\times7$ nilpotent matrices, $J_3(0)\oplus J_3(0)\oplus J_1(0)$ and $J_3(0)\oplus J_2(0)\oplus J_2(0)$, with $p = x^7$, $m = x^3$, rank $4$ in both cases, that are not similar: $d_2 = 5$ versus $6$. Only the full sequence $d_k$ tells them apart.

Over $\mathbb{R}$: a real matrix with non-real eigenvalues has no real Jordan form in this sense; the *real Jordan form* replaces each conjugate pair of blocks by a block built from $2\times2$ rotation-scaling matrices $\begin{pmatrix}a & -b\\ b & a\end{pmatrix}$, as in LA-7's treatment of complex eigenvalues. Horn & Johnson §3.4 has it; we won't need it.

### 3.11 The numerical warning

Take
$$A_\varepsilon = \begin{pmatrix}\lambda & 1\\ \varepsilon & \lambda\end{pmatrix}.$$
Its characteristic polynomial is $(x-\lambda)^2 - \varepsilon$, eigenvalues $\lambda\pm\sqrt{\varepsilon}$, eigenvectors $(1, \pm\sqrt\varepsilon)^\top$. For every $\varepsilon\ne0$ the eigenvalues are distinct, so $A_\varepsilon$ is diagonalizable, with Jordan form $\operatorname{diag}(\lambda+\sqrt\varepsilon, \lambda-\sqrt\varepsilon)$. At $\varepsilon = 0$ the Jordan form is $J_2(\lambda)$. Three things go wrong at once as $\varepsilon\to0$:

1. The Jordan form jumps: it is a **discontinuous** function of the matrix entries. No amount of care fixes this; it is the definition.
2. The eigenvectors collapse onto each other. The angle between $(1,\sqrt\varepsilon)$ and $(1,-\sqrt\varepsilon)$ is $2\arctan\sqrt\varepsilon\approx2\sqrt\varepsilon$, and the condition number of the eigenvector matrix $V = \begin{pmatrix}1&1\\ \sqrt\varepsilon&-\sqrt\varepsilon\end{pmatrix}$ is exactly $1/\sqrt\varepsilon$ (Exercise C2). Diagonalizing $A_\varepsilon$ numerically means multiplying by $V^{-1}$, whose entries are of size $1/\sqrt\varepsilon$.
3. The eigenvalues are hypersensitive: a perturbation of size $\varepsilon$ in one entry moves them by $\sqrt\varepsilon$. At $\varepsilon = 10^{-16}$, machine precision, the eigenvalues move by $10^{-8}$, eight orders of magnitude more than the perturbation. For a block of size $s$ the exponent is $1/s$.

Since every floating-point matrix is an $\varepsilon$-perturbation of itself, a numerical routine cannot decide whether you meant $J_2(\lambda)$ or two nearly equal eigenvalues, and whichever it decides, the basis it returns is terribly conditioned. **Never compute a Jordan form in floating point.** This is not a limitation of current software. Golub and Wilkinson's 1976 paper "Ill-conditioned eigensystems and the computation of the Jordan canonical form" is the standard reference, and the conclusion has not changed.

What people do instead:

- For eigenvalues and invariant subspaces: the **Schur form** $A = QTQ^*$ (LA-10), computed by a backward-stable algorithm (LA-13), exists for every matrix, and its unitary $Q$ has condition number $1$. For real matrices, the real Schur form with $2\times2$ diagonal blocks for complex pairs.
- For $f(A)$: the **Schur–Parlett** algorithm computes $f(T)$ for the triangular Schur factor by a recurrence on the entries (Parlett's recurrence), grouping nearly equal eigenvalues into blocks and using a Taylor series on each block, then sets $f(A) = Qf(T)Q^*$. `scipy.linalg.funm` does this. For the exponential specifically, `scipy.linalg.expm` uses scaling and squaring with a Padé approximant and never looks at eigenvectors at all.
- In exact arithmetic (integers, rationals, symbols), the Jordan form is fine: `sympy.Matrix.jordan_form()` computes it exactly, which is how the project uses it.

The Jordan form is a theorem about what $A$ *is*: it tells you which $t^ke^{\lambda t}$ can appear, when two matrices are similar, what the minimal polynomial is. It is not a thing you compute from measured data.

