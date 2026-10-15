# LA-5: Orthogonality: Inner Products, Projections, Least Squares, Gram–Schmidt, QR

**Time:** ~12 hours · **Prerequisites:** LA-4 (and through it LA-1 to LA-3) · **Unlocks:** LA-8 (spectral theorem), F2 (Fourier series), ST-3 (regression)

> **Conventions in this module.** Vectors are columns, bold lowercase; $A$ is $m \times n$ with columns $\mathbf{a}_1, \dots, \mathbf{a}_n \in \mathbb{R}^m$ and rank $r$. The dot product is $\mathbf{u}\cdot\mathbf{v} = \mathbf{u}^\top\mathbf{v}$ and $\|\mathbf{v}\| = \sqrt{\mathbf{v}^\top\mathbf{v}}$ is the Euclidean norm; for a matrix, $\|M\|$ is the 2-norm $\max_{\mathbf{x}\ne\mathbf{0}}\|M\mathbf{x}\|/\|\mathbf{x}\|$ (used in §3.8), and $\varepsilon \approx 1.1\times10^{-16}$ is the unit roundoff of double precision; "orthogonal" means $\mathbf{u}^\top\mathbf{v} = 0$. $S^\perp$ (read "S perp") is the orthogonal complement of a subspace $S$. For a vector $\mathbf{b}$ and a subspace: $\mathbf{p}$ is the **projection** (a vector), $P$ is the **projection matrix** ($m\times m$), $\hat{\mathbf{x}}$ ("x-hat", $n$ entries) is the coefficient vector with $\mathbf{p} = A\hat{\mathbf{x}}$, and $\mathbf{e} = \mathbf{b} - \mathbf{p}$ is the error; in least squares the same vector is called the residual $\mathbf{r}$. A matrix $Q$ with $Q^\top Q = I$ has **orthonormal columns**; the phrase **orthogonal matrix** is reserved for the square case, where also $QQ^\top = I$. In $A = QR$ we normalize so that the diagonal of $R$ is positive; NumPy does not, see §8. For functions on $[-\pi, \pi]$ the inner product is $\langle f, g\rangle = \int_{-\pi}^{\pi} f(x)g(x)\,dx$ with no $1/\pi$ in front.

## 0. Where this sits

From LA-3 you have the four fundamental subspaces of a matrix, $C(A), N(A), C(A^\top), N(A^\top)$, with dimensions $r, n-r, r, m-r$. From LA-4 you have Rank–Nullity and the view of a matrix as a linear map. What you do not yet have is any notion of *angle* between those subspaces, or any way to answer "$A\mathbf{x} = \mathbf{b}$ has no solution; what is the best $\mathbf{x}$ anyway?" This module supplies both. The one idea is: **the closest point has a perpendicular error.** Everything else (normal equations, projection matrices, Gram–Schmidt, QR) is that sentence written in matrices.

It unlocks LA-8, where orthonormal bases of eigenvectors are the whole point; F2, where a Fourier series is literally a projection onto an orthogonal set of functions; and ST-3, where regression is this module with error bars.

## 1. The problem that forced this idea

On the first night of 1801 Giuseppe Piazzi, in Palermo, found a faint moving object in Taurus. He tracked it for about six weeks, through 11 February, and then it moved too close to the Sun to see. It was Ceres, the first asteroid. The question for the rest of the year was whether anyone could say where it would reappear.

An orbit around the Sun is fixed by six numbers (the orbital elements). Each night's observation gives two numbers (two angles on the sky). Piazzi had about two dozen observations, so there were many more equations than unknowns, and because every observation carries error, there is no set of six elements that satisfies all of them exactly. The equations are *inconsistent*. In the language of LA-3, $\mathbf{b}$ is not in the column space. The methods of the day determined an orbit from a handful of observations and had no systematic way to use the rest of the data.

Carl Friedrich Gauss, then twenty-four, computed an orbit from Piazzi's data, and on the last night of 1801 Franz Xaver von Zach recovered Ceres close to where Gauss said it would be. The orbit-determination method Gauss invented for this is its own story. What matters here is the rule he later described, in his 1809 account of the method, for making many inconsistent equations agree on a single answer: choose the unknowns so that the **sum of the squares of the errors is as small as possible**. Gauss said he had used this rule since 1795; how much of the 1801 Ceres computation leaned on it is something historians still argue about. Adrien-Marie Legendre published it first, in 1805, as an appendix to a book on comet orbits, under the name it still carries, *la méthode des moindres carrés*, the method of least squares. The priority dispute between them was bitter and never resolved.

Why squares? Other rules had been tried: Boscovich (1750s) minimized the sum of the *absolute* errors. That rule has no formula; you have to search. Squares have three advantages that this module makes precise. First, the sum of squares is a smooth function of the unknowns, so calculus finds its minimum and the equations it produces are *linear*. Second, those linear equations have a geometric meaning: the error vector is **perpendicular** to everything you can fit. That is the picture §2 builds and §3 proves. Third (Gauss, 1809): if the observation errors follow the bell curve, least squares is the most probable answer. That third reason belongs to ST-3; this module is about the first two.

One honesty note. Orbital elements enter the astronomical equations nonlinearly. Gauss linearized about a trial orbit and solved a *linear* least-squares problem for the correction, then repeated. The linear problem is this module; the outer loop is Gauss–Newton (ML-2).

## 2. Intuition

### Shadows

Stand a pole on flat ground with the Sun directly overhead. The shadow of the pole's tip is the point on the ground directly beneath it, and the line from tip to shadow is perpendicular to the ground. Among all points on the ground, the shadow is the one *closest* to the tip. Those two facts, perpendicular and closest, are the same fact. The shadow is the **projection** of the tip onto the ground.

With numbers: put the tip at $\mathbf{b} = (1, 2, 3)$ and let the ground be the $xy$-plane. The shadow is $\mathbf{p} = (1, 2, 0)$, and the error $\mathbf{e} = \mathbf{b} - \mathbf{p} = (0, 0, 3)$ points straight up, perpendicular to every floor vector $(x, y, 0)$. Any other floor point $(x, y, 0)$ is farther from the tip: its squared distance is $(x - 1)^2 + (y - 2)^2 + 9 \ge 9$, with equality only at the shadow.

In $\mathbb{R}^2$, project $\mathbf{b}$ onto the line through $\mathbf{a}$:

```
            b
            *
           /|
          / |  e = b - p      (perpendicular to the line)
         /  |
        /   |
   ----*----*-----------  line through a
       0    p = (a·b / a·a) a
```

You derived this in LA-1 (Example 4.5, and the Cauchy–Schwarz proof): $\mathbf{p} = t^*\mathbf{a}$ with $t^* = \mathbf{a}\cdot\mathbf{b}/\mathbf{a}\cdot\mathbf{a}$, and the error $\mathbf{e} = \mathbf{b} - \mathbf{p}$ satisfies $\mathbf{a}\cdot\mathbf{e} = 0$. One unknown ($t$), one perpendicularity condition.

Why is the perpendicular foot the closest point? Take any other point $\mathbf{v}$ on the line. Then $\mathbf{b} - \mathbf{v} = \mathbf{e} + (\mathbf{p} - \mathbf{v})$, and the two pieces are perpendicular (one is off the line at right angles, the other is along it). Pythagoras: $\|\mathbf{b} - \mathbf{v}\|^2 = \|\mathbf{e}\|^2 + \|\mathbf{p} - \mathbf{v}\|^2 \ge \|\mathbf{e}\|^2$. The right angle *is* the reason for the minimum. Keep that: it is the whole proof of Theorem 3.14.

### A plane, two conditions

Now the "ground" is a plane in $\mathbb{R}^3$ spanned by two vectors $\mathbf{a}_1, \mathbf{a}_2$, and $\mathbf{b}$ is a point above it. The shadow $\mathbf{p}$ is in the plane, so $\mathbf{p} = \hat{x}_1\mathbf{a}_1 + \hat{x}_2\mathbf{a}_2$: two unknowns. The error $\mathbf{e} = \mathbf{b} - \mathbf{p}$ must be perpendicular to the *whole plane*. A vector is perpendicular to a plane exactly when it is perpendicular to two vectors that span it. So: two unknowns, two conditions,

$$
\mathbf{a}_1\cdot(\mathbf{b} - \hat{x}_1\mathbf{a}_1 - \hat{x}_2\mathbf{a}_2) = 0, \qquad \mathbf{a}_2\cdot(\mathbf{b} - \hat{x}_1\mathbf{a}_1 - \hat{x}_2\mathbf{a}_2) = 0.
$$

Stack the two columns into $A = (\mathbf{a}_1\ \mathbf{a}_2)$. The two conditions say "each column of $A$, dotted with $\mathbf{b} - A\hat{\mathbf{x}}$, is zero," and the rows of $A^\top$ are the columns of $A$, so they are the single matrix statement $A^\top(\mathbf{b} - A\hat{\mathbf{x}}) = \mathbf{0}$. That is where the transpose comes from. It is not a trick; it is "perpendicular to every column," written once.

You should now be able to guess the main theorem: **the projection of $\mathbf{b}$ onto the column space of $A$ is $\mathbf{p} = A\hat{\mathbf{x}}$ where $A^\top A\hat{\mathbf{x}} = A^\top\mathbf{b}$.** And since $A\mathbf{x} = \mathbf{b}$ has no solution exactly when $\mathbf{b}$ is off the column space, the "best" $\mathbf{x}$ for an unsolvable system is the $\hat{\mathbf{x}}$ whose $A\hat{\mathbf{x}}$ is the shadow of $\mathbf{b}$. Least squares is projection.

### Why orthonormal bases are worth the trouble

Take a basis $\mathbf{q}_1, \mathbf{q}_2$ of the plane that is **orthonormal**: each has length 1 and they are perpendicular. To write $\mathbf{p} = c_1\mathbf{q}_1 + c_2\mathbf{q}_2$, dot both sides with $\mathbf{q}_1$: $\mathbf{q}_1\cdot\mathbf{p} = c_1(\mathbf{q}_1\cdot\mathbf{q}_1) + c_2(\mathbf{q}_1\cdot\mathbf{q}_2) = c_1$. The coordinate along $\mathbf{q}_1$ is a single dot product. No system to solve, no $A^\top A$ to invert. With a skewed basis, finding coordinates means solving $A\mathbf{x} = \mathbf{p}$; with an orthonormal one, every coordinate is independent of every other. That is why the second half of this module builds orthonormal bases (Gram–Schmidt) and packages them as $A = QR$: the hard part of projection, $(A^\top A)^{-1}$, becomes the identity.

How do you get an orthonormal basis from a skewed one? With shadows again. Take $\mathbf{a}_1 = (1, 0)$ and $\mathbf{a}_2 = (1, 1)$. Keep $\mathbf{a}_1$. Subtract from $\mathbf{a}_2$ its shadow on $\mathbf{a}_1$, which is $(1, 0)$; what remains, $(0, 1)$, is perpendicular to $\mathbf{a}_1$. Normalize, and with more vectors repeat: subtract every shadow on the vectors you already have. That is Gram–Schmidt. Now tilt $\mathbf{a}_2$ until it is nearly parallel to $\mathbf{a}_1$. Its shadow is then nearly all of it, and the perpendicular remainder is a tiny vector obtained by subtracting two nearly equal ones. That is where rounding error will bite (§3.8). The same nearly parallel columns make $A^\top A$ worse than $A$: if $A$ stretches one direction $10^4$ times more than another, $A^\top A$, which applies $A$ and then $A^\top$, stretches it $10^8$ times more. That squaring is the $\kappa(A)^2$ of §3.10.

## 3. Formal development

### 3.1 Orthogonal vectors

**Definition 3.1.** Vectors $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$ are **orthogonal**, written $\mathbf{u} \perp \mathbf{v}$, if $\mathbf{u}^\top\mathbf{v} = 0$. The zero vector is orthogonal to everything.

**Theorem 3.2 (Pythagoras).** $\mathbf{u} \perp \mathbf{v}$ if and only if
$$
\|\mathbf{u} + \mathbf{v}\|^2 = \|\mathbf{u}\|^2 + \|\mathbf{v}\|^2. \tag{3.1}
$$

*Proof.* $\|\mathbf{u}+\mathbf{v}\|^2 = (\mathbf{u}+\mathbf{v})^\top(\mathbf{u}+\mathbf{v}) = \|\mathbf{u}\|^2 + 2\,\mathbf{u}^\top\mathbf{v} + \|\mathbf{v}\|^2$. The cross term vanishes exactly when $\mathbf{u}^\top\mathbf{v} = 0$. ∎

**Theorem 3.3.** Nonzero vectors $\mathbf{v}_1, \dots, \mathbf{v}_k$ that are pairwise orthogonal are linearly independent.

*Proof.* Suppose $c_1\mathbf{v}_1 + \cdots + c_k\mathbf{v}_k = \mathbf{0}$. Dot both sides with $\mathbf{v}_i$: every term $c_j\mathbf{v}_i^\top\mathbf{v}_j$ with $j \ne i$ is zero, leaving $c_i\|\mathbf{v}_i\|^2 = 0$. Since $\mathbf{v}_i \ne \mathbf{0}$, $c_i = 0$. This holds for each $i$. ∎

The move in that proof, "dot with $\mathbf{v}_i$ to isolate $c_i$," is the engine of the whole module. In §3.7 it becomes "coefficients are dot products," and in F2 it becomes the formula for Fourier coefficients.

### 3.2 Orthogonal subspaces and the orthogonal complement

**Definition 3.4.** Subspaces $S, T \subseteq \mathbb{R}^n$ are **orthogonal** if $\mathbf{s}^\top\mathbf{t} = 0$ for every $\mathbf{s} \in S$ and $\mathbf{t} \in T$.

Two orthogonal subspaces meet only at $\mathbf{0}$: if $\mathbf{v}$ is in both, then $\mathbf{v}^\top\mathbf{v} = 0$, so $\mathbf{v} = \mathbf{0}$. Note that the floor and a wall of a room are *not* orthogonal subspaces: they share the line where they meet. Orthogonal subspaces are more restrictive than "at right angles" in everyday speech.

**Definition 3.5.** The **orthogonal complement** of a subspace $S \subseteq \mathbb{R}^n$ is
$$
S^\perp = \{\mathbf{v} \in \mathbb{R}^n : \mathbf{v}^\top\mathbf{s} = 0 \text{ for all } \mathbf{s} \in S\}. \tag{3.2}
$$

**Theorem 3.6.** $S^\perp$ is a subspace of $\mathbb{R}^n$.

*Proof.* $\mathbf{0} \in S^\perp$. If $\mathbf{v}, \mathbf{w} \in S^\perp$ and $c \in \mathbb{R}$, then for every $\mathbf{s} \in S$: $(\mathbf{v} + \mathbf{w})^\top\mathbf{s} = \mathbf{v}^\top\mathbf{s} + \mathbf{w}^\top\mathbf{s} = 0$ and $(c\mathbf{v})^\top\mathbf{s} = c\,\mathbf{v}^\top\mathbf{s} = 0$. Closed under both operations. ∎

**Theorem 3.7 (Dimension of the complement).** Let $S$ be a subspace of $\mathbb{R}^n$ with $\dim S = k$. Then
(a) $\dim S^\perp = n - k$;
(b) every $\mathbf{x} \in \mathbb{R}^n$ can be written in exactly one way as $\mathbf{x} = \mathbf{s} + \mathbf{t}$ with $\mathbf{s} \in S$, $\mathbf{t} \in S^\perp$;
(c) $(S^\perp)^\perp = S$.

▶ **Try it first.** Put a basis of $S$ into the *rows* of a $k \times n$ matrix $B$. What is $S^\perp$ in terms of $B$? Then use Rank–Nullity (LA-4).

**Derivation.**

1. Let $\mathbf{s}_1, \dots, \mathbf{s}_k$ be a basis of $S$ and let $B$ be the $k \times n$ matrix with these as rows. A vector $\mathbf{v}$ is orthogonal to every $\mathbf{s} \in S$ if and only if it is orthogonal to each basis vector (if $\mathbf{v}^\top\mathbf{s}_i = 0$ for all $i$ then $\mathbf{v}^\top\sum_i c_i\mathbf{s}_i = \sum_i c_i\mathbf{v}^\top\mathbf{s}_i = 0$; the converse is immediate). And "$\mathbf{v}^\top\mathbf{s}_i = 0$ for each $i$" says exactly $B\mathbf{v} = \mathbf{0}$ by the row form of matrix–vector multiplication (LA-1, eq. (3.2)). So $S^\perp = N(B)$.
2. $B$ has $k$ independent rows, so $\operatorname{rank} B = k$ (row rank equals column rank, LA-3). Rank–Nullity gives $\dim N(B) = n - k$. That is (a).
3. For (b), existence: $S \cap S^\perp = \{\mathbf{0}\}$ (shown after Definition 3.4), so a basis of $S$ together with a basis of $S^\perp$ is a set of $k + (n-k) = n$ vectors that is independent (a dependence $\sum a_i\mathbf{s}_i + \sum b_j\mathbf{t}_j = \mathbf{0}$ puts $\sum a_i\mathbf{s}_i = -\sum b_j\mathbf{t}_j$ in both subspaces, hence equal to $\mathbf{0}$, hence all coefficients vanish). $n$ independent vectors in $\mathbb{R}^n$ form a basis, so every $\mathbf{x}$ is a combination of them: $\mathbf{x} = \mathbf{s} + \mathbf{t}$.
4. Uniqueness: if $\mathbf{s} + \mathbf{t} = \mathbf{s}' + \mathbf{t}'$ then $\mathbf{s} - \mathbf{s}' = \mathbf{t}' - \mathbf{t}$ lies in $S \cap S^\perp = \{\mathbf{0}\}$.
5. For (c): every $\mathbf{s} \in S$ is orthogonal to every element of $S^\perp$ by definition, so $S \subseteq (S^\perp)^\perp$. By (a) applied twice, $\dim (S^\perp)^\perp = n - (n - k) = k = \dim S$. A subspace contained in another of the same dimension equals it. ∎

■ **What just happened.** The complement of a $k$-dimensional subspace is the null space of a rank-$k$ matrix, so its dimension is forced by Rank–Nullity. Part (b) is the statement that $\mathbb{R}^n = S \oplus S^\perp$: every vector splits into a piece in $S$ and a piece perpendicular to $S$, and the split is unique. The piece in $S$ is about to be called the projection.

### 3.3 The four subspaces, with right angles

**Lemma 3.8.** $\mathbf{v} \perp \operatorname{span}\{\mathbf{w}_1, \dots, \mathbf{w}_k\}$ if and only if $\mathbf{v}^\top\mathbf{w}_i = 0$ for each $i$.

This was step 1 of the last derivation. It is the reason transposes appear: "orthogonal to a subspace" reduces to finitely many dot products, one per spanning vector, and a stack of dot products is a matrix product.

**Theorem 3.9 (Fundamental Theorem of Linear Algebra, part 2).** For any $m \times n$ matrix $A$,
$$
N(A) = C(A^\top)^\perp \quad\text{in } \mathbb{R}^n, \qquad N(A^\top) = C(A)^\perp \quad\text{in } \mathbb{R}^m. \tag{3.3}
$$
Equivalently $C(A^\top) = N(A)^\perp$ and $C(A) = N(A^\top)^\perp$.

▶ **Try it first.** Write out what $A\mathbf{x} = \mathbf{0}$ says row by row. Each row of $A$ is a column of $A^\top$.

**Derivation.**

1. $\mathbf{x} \in N(A)$ means $A\mathbf{x} = \mathbf{0}$, i.e. by the row form, $\mathbf{r}_i^\top\mathbf{x} = 0$ for every row $\mathbf{r}_i^\top$ of $A$.
2. The rows of $A$ are the columns of $A^\top$, so step 1 says $\mathbf{x}$ is orthogonal to each column of $A^\top$.
3. By Lemma 3.8, that is the same as $\mathbf{x} \perp \operatorname{span}(\text{columns of } A^\top) = C(A^\top)$. So $\mathbf{x} \in N(A) \iff \mathbf{x} \in C(A^\top)^\perp$.
4. Apply steps 1–3 to the matrix $A^\top$ in place of $A$: $N(A^\top) = C((A^\top)^\top)^\perp = C(A)^\perp$.
5. The "equivalently" follows from Theorem 3.7(c): taking $\perp$ of both sides of $N(A) = C(A^\top)^\perp$ gives $N(A)^\perp = C(A^\top)$. ∎

■ **What just happened.** LA-3 told you the dimensions of the four subspaces add up: $r + (n-r) = n$ and $r + (m-r) = m$. Now you know *how* they fit: the row space and null space are perpendicular complements in $\mathbb{R}^n$, the column space and left null space are perpendicular complements in $\mathbb{R}^m$. The dimension count was never a coincidence; it is Theorem 3.7(a).

```
        R^n  (inputs)                                   R^m  (outputs)

      C(A^T)  row space  ─────────  A  ──────────►   C(A)  column space
       dim r                                           dim r
         │                                               │
         │  ⟂  (Theorem 3.9)                             │  ⟂  (Theorem 3.9)
         │                                               │
       N(A)   null space  ────────  A  ──────────►   0   N(A^T)  left null space
       dim n − r                                       dim m − r
```

Every $\mathbf{x} \in \mathbb{R}^n$ splits as $\mathbf{x} = \mathbf{x}_{\text{row}} + \mathbf{x}_{\text{null}}$ (Theorem 3.7(b)), and $A\mathbf{x} = A\mathbf{x}_{\text{row}}$: the matrix only sees the row-space part. Restricted to the row space, $A$ is a bijection onto the column space. One-to-one: if $\mathbf{x}, \mathbf{x}' \in C(A^\top)$ and $A\mathbf{x} = A\mathbf{x}'$, then $\mathbf{x} - \mathbf{x}'$ lies in $C(A^\top) \cap N(A)$, which is $\{\mathbf{0}\}$ because the two are orthogonal (Theorem 3.9 and the remark after Definition 3.4). Onto: every vector $A\mathbf{x}$ of $C(A)$ equals $A\mathbf{x}_{\text{row}}$, the image of a row-space vector. So every matrix of rank $r$ hides an invertible map between two $r$-dimensional spaces. The project checks the splitting and the one-to-one property.

**Corollary 3.10 (Fredholm alternative).** $A\mathbf{x} = \mathbf{b}$ has a solution if and only if $\mathbf{b}$ is orthogonal to every solution $\mathbf{y}$ of $A^\top\mathbf{y} = \mathbf{0}$.

*Proof.* Solvable iff $\mathbf{b} \in C(A)$ (LA-1, Theorem 3.6) iff $\mathbf{b} \in N(A^\top)^\perp$ by (3.3). ∎

### 3.4 Projection onto a line, then onto a subspace

**Onto a line.** For $\mathbf{a} \ne \mathbf{0}$ (LA-1 Example 4.5),
$$
\mathbf{p} = \frac{\mathbf{a}^\top\mathbf{b}}{\mathbf{a}^\top\mathbf{a}}\,\mathbf{a} = \underbrace{\frac{\mathbf{a}\mathbf{a}^\top}{\mathbf{a}^\top\mathbf{a}}}_{P}\,\mathbf{b}. \tag{3.4}
$$
Read the second form carefully: $\mathbf{a}\mathbf{a}^\top$ is an $m \times m$ matrix (column times row, rank 1), $\mathbf{a}^\top\mathbf{a}$ is a number. The matrix $P$ depends only on the line, not on $\mathbf{b}$. The error $\mathbf{e} = \mathbf{b} - \mathbf{p}$ satisfies $\mathbf{a}^\top\mathbf{e} = \mathbf{a}^\top\mathbf{b} - \frac{\mathbf{a}^\top\mathbf{b}}{\mathbf{a}^\top\mathbf{a}}\mathbf{a}^\top\mathbf{a} = 0$.

**Onto a subspace.** Let $A$ be $m \times n$ with *independent* columns, so $C(A)$ is an $n$-dimensional subspace of $\mathbb{R}^m$, and let $\mathbf{b} \in \mathbb{R}^m$. We want the point $\mathbf{p} \in C(A)$ with $\mathbf{b} - \mathbf{p} \perp C(A)$. (Theorem 3.7(b) says exactly one such $\mathbf{p}$ exists; now we compute it.)

▶ **Try it first.** Write $\mathbf{p} = A\hat{\mathbf{x}}$. Write "$\mathbf{b} - A\hat{\mathbf{x}}$ is orthogonal to each column $\mathbf{a}_j$" as $n$ equations. Stack them into one matrix equation and solve for $\hat{\mathbf{x}}$.

**Derivation.**

1. Since $\mathbf{p} \in C(A)$, it is a combination of the columns: $\mathbf{p} = A\hat{\mathbf{x}}$ for some $\hat{\mathbf{x}} \in \mathbb{R}^n$ (LA-1, eq. (3.1)).
2. The error $\mathbf{e} = \mathbf{b} - A\hat{\mathbf{x}}$ must be orthogonal to $C(A)$. By Lemma 3.8 it suffices that $\mathbf{a}_j^\top\mathbf{e} = 0$ for $j = 1, \dots, n$.
3. The $n$ numbers $\mathbf{a}_1^\top\mathbf{e}, \dots, \mathbf{a}_n^\top\mathbf{e}$ are the entries of $A^\top\mathbf{e}$, because the rows of $A^\top$ are $\mathbf{a}_1^\top, \dots, \mathbf{a}_n^\top$. So the $n$ conditions are one equation:
   $$A^\top(\mathbf{b} - A\hat{\mathbf{x}}) = \mathbf{0}. \tag{3.5}$$
4. Rearranged, these are the **normal equations**:
   $$A^\top A\,\hat{\mathbf{x}} = A^\top\mathbf{b}. \tag{3.6}$$
   $A^\top A$ is $n \times n$ and symmetric ($(A^\top A)^\top = A^\top A^{\top\top} = A^\top A$).
5. Theorem 3.11 below shows $A^\top A$ is invertible when the columns of $A$ are independent. Then
   $$\hat{\mathbf{x}} = (A^\top A)^{-1}A^\top\mathbf{b}, \tag{3.7}$$
   $$\mathbf{p} = A\hat{\mathbf{x}} = A(A^\top A)^{-1}A^\top\mathbf{b}, \qquad P = A(A^\top A)^{-1}A^\top. \tag{3.8}$$ ∎

■ **What just happened.** Three objects came out, and you must keep them apart. $\hat{\mathbf{x}} \in \mathbb{R}^n$ is the *recipe*: how much of each column. $\mathbf{p} \in \mathbb{R}^m$ is the *shadow*. $P$ is the $m \times m$ *machine* that casts shadows; it depends on $A$ only. When $n = 1$, (3.8) is (3.4): $(A^\top A)^{-1}$ is the number $1/\mathbf{a}^\top\mathbf{a}$. The transpose in (3.5) is Lemma 3.8 in matrix form: "perpendicular to every column of $A$" is "$A^\top$ times it is zero."

**Theorem 3.11.** $A^\top A$ is invertible if and only if the columns of $A$ are linearly independent. More precisely, $N(A^\top A) = N(A)$.

▶ **Try it first.** Suppose $A^\top A\mathbf{x} = \mathbf{0}$. Multiply on the left by $\mathbf{x}^\top$ and recognize a squared length.

**Derivation.**

1. If $A\mathbf{x} = \mathbf{0}$ then $A^\top A\mathbf{x} = A^\top\mathbf{0} = \mathbf{0}$. So $N(A) \subseteq N(A^\top A)$.
2. Conversely, if $A^\top A\mathbf{x} = \mathbf{0}$ then $\mathbf{x}^\top A^\top A\mathbf{x} = 0$. But $\mathbf{x}^\top A^\top A\mathbf{x} = (A\mathbf{x})^\top(A\mathbf{x}) = \|A\mathbf{x}\|^2$. A vector of length zero is zero: $A\mathbf{x} = \mathbf{0}$. So $N(A^\top A) \subseteq N(A)$.
3. Hence $N(A^\top A) = N(A)$. A square matrix is invertible iff its null space is $\{\mathbf{0}\}$ (LA-2), and $N(A) = \{\mathbf{0}\}$ iff the columns of $A$ are independent (LA-3). ∎

■ **What just happened.** The symmetric matrix $A^\top A$ (the *Gram matrix*; its entries are all the dot products $\mathbf{a}_i^\top\mathbf{a}_j$) inherits its invertibility from $A$. The identity $\mathbf{x}^\top A^\top A\mathbf{x} = \|A\mathbf{x}\|^2 \ge 0$ says $A^\top A$ is *positive semidefinite*, and positive definite when the columns are independent. That is the LA-8 definition arriving early; you will see $A^\top A$ again as the matrix whose eigenvalues are the squared singular values in LA-9.

If the columns of $A$ are dependent, $A^\top A$ is singular and (3.7) is meaningless, but the projection $\mathbf{p}$ still exists and is unique (Theorem 3.7(b)); only the recipe $\hat{\mathbf{x}}$ is non-unique. Throw away dependent columns, or wait for the pseudoinverse in LA-9.

### 3.5 What makes a matrix a projection

**Theorem 3.12.** Let $P = A(A^\top A)^{-1}A^\top$ with $A$ of independent columns. Then $P^2 = P$ and $P^\top = P$. Conversely, any matrix $P$ with $P^2 = P$ and $P^\top = P$ is the orthogonal projection onto its column space $C(P)$.

*Proof of the forward direction.* $P^2 = A(A^\top A)^{-1}\underbrace{A^\top A(A^\top A)^{-1}}_{I}A^\top = A(A^\top A)^{-1}A^\top = P$. For symmetry, use $(BCD)^\top = D^\top C^\top B^\top$ (LA-1 Theorem 3.8, twice) and $((A^\top A)^{-1})^\top = ((A^\top A)^\top)^{-1} = (A^\top A)^{-1}$: $P^\top = A^{\top\top}(A^\top A)^{-1}{}^\top A^\top = A(A^\top A)^{-1}A^\top = P$. ∎

▶ **Try it first** (converse). Given $P^2 = P = P^\top$ and any $\mathbf{b}$, split $\mathbf{b} = P\mathbf{b} + (I - P)\mathbf{b}$. Show the first piece is in $C(P)$ and the second is orthogonal to every vector $P\mathbf{w}$ in $C(P)$.

**Derivation.**

1. $P\mathbf{b} \in C(P)$ by definition of the column space.
2. For any $\mathbf{w}$: $(P\mathbf{w})^\top(I - P)\mathbf{b} = \mathbf{w}^\top P^\top(I - P)\mathbf{b} = \mathbf{w}^\top(P - P^2)\mathbf{b} = \mathbf{w}^\top\,0\,\mathbf{b} = 0$, using $P^\top = P$ then $P^2 = P$. So $(I - P)\mathbf{b} \perp C(P)$.
3. By Theorem 3.7(b) the decomposition of $\mathbf{b}$ into (piece in $C(P)$) + (piece in $C(P)^\perp$) is unique, and it is the definition of orthogonal projection. So $P\mathbf{b}$ is the projection of $\mathbf{b}$ onto $C(P)$. ∎

■ **What just happened.** Two algebraic properties, idempotent and symmetric, *characterize* orthogonal projections. You can test a matrix for being a projection without knowing what it projects onto. Each property has a meaning: $P^2 = P$ says projecting a shadow changes nothing; $P^\top = P$ says the projection is *orthogonal* (a non-symmetric idempotent is an oblique projection, a shadow cast by a slanted sun).

**Theorem 3.13.** If $P$ is the orthogonal projection onto $S$, then $I - P$ is the orthogonal projection onto $S^\perp$. For $P = A(A^\top A)^{-1}A^\top$: $I - P$ projects onto $C(A)^\perp = N(A^\top)$.

*Proof.* $(I-P)^2 = I - 2P + P^2 = I - P$ and $(I - P)^\top = I - P^\top = I - P$, so by Theorem 3.12 $I - P$ is the orthogonal projection onto its column space. That column space is $\{(I-P)\mathbf{b}\} = \{\mathbf{b} - \mathbf{p}\}$, the set of all errors, which is contained in $S^\perp$; and every $\mathbf{t} \in S^\perp$ has $P\mathbf{t} = \mathbf{0}$ (its projection onto $S$ is zero since $\mathbf{t} = \mathbf{0} + \mathbf{t}$ is its unique split), so $(I - P)\mathbf{t} = \mathbf{t}$ and $S^\perp \subseteq C(I - P)$. The last claim is (3.3). ∎

Also $P + (I - P) = I$ and $P(I - P) = P - P^2 = 0$: the two shadows reconstruct $\mathbf{b}$ and the second shadow of the first is nothing.

**Theorem 3.14 (Closest point).** Let $\mathbf{p}$ be the orthogonal projection of $\mathbf{b}$ onto a subspace $S$. Then for every $\mathbf{v} \in S$, $\|\mathbf{b} - \mathbf{v}\| \ge \|\mathbf{b} - \mathbf{p}\|$, with equality only if $\mathbf{v} = \mathbf{p}$.

*Proof.* $\mathbf{b} - \mathbf{v} = (\mathbf{b} - \mathbf{p}) + (\mathbf{p} - \mathbf{v})$. The first piece is in $S^\perp$, the second in $S$ (both $\mathbf{p}$ and $\mathbf{v}$ are). By Pythagoras (3.1), $\|\mathbf{b} - \mathbf{v}\|^2 = \|\mathbf{b} - \mathbf{p}\|^2 + \|\mathbf{p} - \mathbf{v}\|^2 \ge \|\mathbf{b} - \mathbf{p}\|^2$, with equality iff $\mathbf{p} = \mathbf{v}$. ∎

That is the §2 right-triangle argument, now in $m$ dimensions. "Perpendicular error" and "closest point" are the same thing.

### 3.6 Least squares

The system $A\mathbf{x} = \mathbf{b}$ with $m > n$ (more equations than unknowns) usually has no solution: $\mathbf{b} \notin C(A)$. The **least-squares** problem asks for the $\mathbf{x}$ that makes $A\mathbf{x}$ as close to $\mathbf{b}$ as possible:
$$
\text{minimize } f(\mathbf{x}) = \|A\mathbf{x} - \mathbf{b}\|^2 = \sum_{i=1}^m\big((A\mathbf{x})_i - b_i\big)^2. \tag{3.9}
$$

**Theorem 3.15 (Least squares).** $\hat{\mathbf{x}}$ minimizes (3.9) if and only if it satisfies the normal equations (3.6). If the columns of $A$ are independent, the minimizer is unique and given by (3.7); the minimum value is $\|\mathbf{e}\|^2 = \|\mathbf{b}\|^2 - \|\mathbf{p}\|^2$.

Two proofs, geometric and calculus. They must agree, and seeing why they agree is the point.

*Geometric proof.* As $\mathbf{x}$ ranges over $\mathbb{R}^n$, $A\mathbf{x}$ ranges over $C(A)$. By Theorem 3.14 the closest point of $C(A)$ to $\mathbf{b}$ is the projection $\mathbf{p}$, so $f$ is minimized exactly by those $\hat{\mathbf{x}}$ with $A\hat{\mathbf{x}} = \mathbf{p}$, i.e. by those with $\mathbf{b} - A\hat{\mathbf{x}} \perp C(A)$, i.e. (step 3 of §3.4) by the solutions of (3.5), which is (3.6). The minimum value is $\|\mathbf{b} - \mathbf{p}\|^2 = \|\mathbf{b}\|^2 - \|\mathbf{p}\|^2$ by Pythagoras, since $\mathbf{p} \perp \mathbf{e}$. ∎

▶ **Try it first** (calculus). Expand $f(\mathbf{x}) = (A\mathbf{x} - \mathbf{b})^\top(A\mathbf{x} - \mathbf{b})$ into three terms. Differentiate with respect to $x_k$. (For the middle term, write $\mathbf{b}^\top A\mathbf{x} = \sum_i (A^\top\mathbf{b})_i x_i$. Use summation letters other than $k$.)

**Derivation (calculus).**

1. Expand: $f(\mathbf{x}) = \mathbf{x}^\top A^\top A\mathbf{x} - 2\,\mathbf{b}^\top A\mathbf{x} + \mathbf{b}^\top\mathbf{b}$, using $(A\mathbf{x})^\top\mathbf{b} = \mathbf{b}^\top A\mathbf{x}$ (a number equals its transpose).
2. Write $G = A^\top A$ (symmetric, entries $g_{il}$) and $\mathbf{c} = A^\top\mathbf{b}$. Then $f = \sum_{i,l} g_{il}x_ix_l - 2\sum_i c_ix_i + \|\mathbf{b}\|^2$. The summation letters $i, l$ run over $1, \dots, n$; $k$ is reserved for the variable we differentiate by.
3. Fix $k$. The product rule gives $\partial(x_ix_l)/\partial x_k = \delta_{ik}x_l + x_i\delta_{lk}$, where $\delta_{ik}$ is $1$ if $i = k$ and $0$ otherwise. Each $\delta$ collapses a double sum to a single sum:
   $$\frac{\partial}{\partial x_k}\sum_{i,l} g_{il}x_ix_l = \sum_{i,l} g_{il}(\delta_{ik}x_l + x_i\delta_{lk}) = \sum_l g_{kl}x_l + \sum_i g_{ik}x_i.$$
   By symmetry $g_{ik} = g_{ki}$, so the two single sums are equal, and with $\partial(\sum_i c_ix_i)/\partial x_k = c_k$,
   $$\frac{\partial f}{\partial x_k} = 2\sum_l g_{kl}x_l - 2c_k = 2\,(G\mathbf{x} - \mathbf{c})_k.$$
   The diagonal term needs no special care: $g_{kk}x_k^2$ is the term $i = l = k$, and the product rule gives it $2g_{kk}x_k$, one copy in each single sum.
4. Collect the components into the gradient:
   $$\nabla f(\mathbf{x}) = 2A^\top A\mathbf{x} - 2A^\top\mathbf{b} = 2A^\top(A\mathbf{x} - \mathbf{b}). \tag{3.10}$$
5. Set $\nabla f = \mathbf{0}$: $A^\top(A\hat{\mathbf{x}} - \mathbf{b}) = \mathbf{0}$, which is (3.5), which is (3.6).
6. Is it a minimum? The Hessian $\partial^2 f/\partial x_j\partial x_k = 2g_{jk}$, i.e. $\nabla^2 f = 2A^\top A$, and $\mathbf{d}^\top(2A^\top A)\mathbf{d} = 2\|A\mathbf{d}\|^2 \ge 0$ for every direction $\mathbf{d}$. A quadratic with positive semidefinite Hessian is convex, so every stationary point is a global minimum. With independent columns the Hessian is positive definite and the minimum is unique. ∎

■ **What just happened.** Calculus says "set the gradient to zero"; geometry says "make the error perpendicular." Equation (3.10) shows they are the same sentence: the gradient *is* $2A^\top$ times the error, and it vanishes exactly when the error is orthogonal to every column. This is also your first matrix-calculus formula; LA-14 derives $\nabla(\mathbf{x}^\top G\mathbf{x}) = 2G\mathbf{x}$ for symmetric $G$ in general, and ML-2 minimizes (3.9) by gradient descent instead of solving (3.6) when $n$ is huge.

**Fitting a line.** Data $(t_1, b_1), \dots, (t_m, b_m)$; model $b = C + Dt$. Unknowns $\mathbf{x} = (C, D)^\top$. Each data point is one equation $C + Dt_i = b_i$:
$$
A = \begin{pmatrix} 1 & t_1 \\ 1 & t_2 \\ \vdots & \vdots \\ 1 & t_m \end{pmatrix}, \qquad
A^\top A = \begin{pmatrix} m & \sum t_i \\ \sum t_i & \sum t_i^2 \end{pmatrix}, \qquad
A^\top\mathbf{b} = \begin{pmatrix} \sum b_i \\ \sum t_ib_i \end{pmatrix}.
$$
$A$ is the **design matrix**: one row per observation, one column per unknown parameter. Its column space is the set of all outputs the model can produce; least squares projects the data onto it.

**Fitting a polynomial.** Model $b = c_0 + c_1t + \cdots + c_dt^d$. The design matrix has columns $\mathbf{1}, \mathbf{t}, \mathbf{t}^2, \dots, \mathbf{t}^d$ (entrywise powers). The model is nonlinear in $t$ but *linear in the unknowns* $c_k$, which is all least squares needs. The same goes for $b = c_1\sin t + c_2\cos t$, or any $b = \sum_k c_k\phi_k(t)$ with fixed basis functions $\phi_k$. Example 4.5 fits a quadratic to sensor readings.

### 3.7 Orthonormal bases and orthogonal matrices

**Definition 3.16.** Vectors $\mathbf{q}_1, \dots, \mathbf{q}_n \in \mathbb{R}^m$ are **orthonormal** if $\mathbf{q}_i^\top\mathbf{q}_j = 0$ for $i \ne j$ and $\mathbf{q}_i^\top\mathbf{q}_i = 1$. Put them as columns of an $m \times n$ matrix $Q$. Then
$$
Q^\top Q = I_n, \tag{3.11}
$$
because $(Q^\top Q)_{ij} = \mathbf{q}_i^\top\mathbf{q}_j$. If $Q$ is square ($m = n$), it is called an **orthogonal matrix**; then $Q^\top = Q^{-1}$ and also $QQ^\top = I$ (for square matrices a one-sided inverse is two-sided: LA-4, the corollary to Theorem 3.18; the proof of Theorem 3.17 below also derives $QQ^\top = I$ directly).

The name is unfortunate (the columns are ortho*normal*, not merely orthogonal), but it is universal. Rectangular $Q$ with orthonormal columns has $Q^\top Q = I_n$ but $QQ^\top \ne I_m$; see Theorem 3.17.

**Theorem 3.17 (Projection with an orthonormal basis).** If $Q$ has orthonormal columns, the projection onto $C(Q)$ is
$$
P = QQ^\top, \qquad \hat{\mathbf{x}} = Q^\top\mathbf{b}, \qquad \mathbf{p} = \sum_{i=1}^n(\mathbf{q}_i^\top\mathbf{b})\,\mathbf{q}_i. \tag{3.12}
$$
If $Q$ is square, then $QQ^\top = I$ and every $\mathbf{b} \in \mathbb{R}^n$ expands as
$$
\mathbf{b} = \sum_{i=1}^n(\mathbf{q}_i^\top\mathbf{b})\,\mathbf{q}_i, \qquad \|\mathbf{b}\|^2 = \sum_{i=1}^n(\mathbf{q}_i^\top\mathbf{b})^2 \quad\text{(Parseval)}. \tag{3.13}
$$

*Proof.* Substitute $A = Q$ in (3.7)–(3.8) and use (3.11): $(Q^\top Q)^{-1} = I$, so $\hat{\mathbf{x}} = Q^\top\mathbf{b}$ and $P = QQ^\top$. The $i$-th entry of $Q^\top\mathbf{b}$ is $\mathbf{q}_i^\top\mathbf{b}$, and $\mathbf{p} = Q\hat{\mathbf{x}} = \sum_i\hat{x}_i\mathbf{q}_i$. If $Q$ is square, $C(Q) = \mathbb{R}^n$, so the projection of $\mathbf{b}$ is $\mathbf{b}$ itself: $QQ^\top\mathbf{b} = \mathbf{b}$ for all $\mathbf{b}$, hence $QQ^\top = I$. Parseval: $\|\mathbf{b}\|^2 = \mathbf{b}^\top QQ^\top\mathbf{b} = (Q^\top\mathbf{b})^\top(Q^\top\mathbf{b}) = \|Q^\top\mathbf{b}\|^2$. ∎

**Bessel's inequality.** If $Q$ is $m\times n$ with $n < m$, the projection misses the part of $\mathbf{b}$ outside $C(Q)$, and Parseval becomes an inequality:
$$
\sum_{i=1}^n(\mathbf{q}_i^\top\mathbf{b})^2 = \|\mathbf{p}\|^2 \le \|\mathbf{b}\|^2, \quad\text{with equality if and only if } \mathbf{b} \in C(Q).
$$
*Proof.* $\|\mathbf{p}\|^2 = \|QQ^\top\mathbf{b}\|^2 = \mathbf{b}^\top Q(Q^\top Q)Q^\top\mathbf{b} = \|Q^\top\mathbf{b}\|^2$ by (3.11), and $\|\mathbf{b}\|^2 = \|\mathbf{p}\|^2 + \|\mathbf{e}\|^2$ by Pythagoras, since $\mathbf{e} \perp \mathbf{p}$ (Theorem 3.15). Equality holds iff $\mathbf{e} = \mathbf{0}$. ∎ A shadow is never longer than the thing that casts it. §3.12 and F2 apply this to Fourier coefficients.

This is the §2 promise kept: with an orthonormal basis, **coefficients are dot products** and the Gram matrix to invert is $I$. The sum $\sum_i(\mathbf{q}_i^\top\mathbf{b})\mathbf{q}_i$ in (3.13) is a Fourier series in disguise; §3.12 and F2 make the disguise explicit.

**Theorem 3.18 (Orthogonal matrices preserve geometry).** If $Q$ is $n\times n$ with $Q^\top Q = I$, then for all $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$:
(a) $(Q\mathbf{x})^\top(Q\mathbf{y}) = \mathbf{x}^\top\mathbf{y}$ (dot products are preserved);
(b) $\|Q\mathbf{x}\| = \|\mathbf{x}\|$ (lengths are preserved);
(c) the angle between $Q\mathbf{x}$ and $Q\mathbf{y}$ equals the angle between $\mathbf{x}$ and $\mathbf{y}$.

*Proof.* (a) $(Q\mathbf{x})^\top(Q\mathbf{y}) = \mathbf{x}^\top Q^\top Q\mathbf{y} = \mathbf{x}^\top I\mathbf{y} = \mathbf{x}^\top\mathbf{y}$. (b) Put $\mathbf{y} = \mathbf{x}$ in (a) and take square roots. (c) By LA-1's angle formula, $\cos\theta = \mathbf{x}^\top\mathbf{y}/(\|\mathbf{x}\|\|\mathbf{y}\|)$; numerator and denominator are unchanged by (a) and (b). ∎

Conversely, a square matrix that preserves all lengths is orthogonal (Exercise B7). So orthogonal matrices are exactly the **linear maps that preserve length**: rotations, reflections, and their products. (They are also exactly the distance-preserving maps $T$ of $\mathbb{R}^n$ with $T(\mathbf{0}) = \mathbf{0}$, but that needs one more step: such a $T$ is automatically linear. It preserves lengths and distances, so by polarization it preserves dot products; expanding $\|T(\mathbf{x} + \mathbf{y}) - T\mathbf{x} - T\mathbf{y}\|^2$ and $\|T(c\mathbf{x}) - cT\mathbf{x}\|^2$ into dot products gives $0$ for both.) Three families you will meet constantly:

- **Rotation** in $\mathbb{R}^2$ by $\theta$: $\begin{pmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{pmatrix}$. Columns are unit and perpendicular. Check.
- **Reflection** in $\mathbb{R}^2$ across the line at angle $\theta/2$: $\begin{pmatrix}\cos\theta & \sin\theta\\ \sin\theta & -\cos\theta\end{pmatrix}$. Symmetric, and its own inverse. (§3.11 generalizes this.)
- **Permutation** matrices (identity with rows reordered, LA-2): each column is a standard basis vector, and distinct columns are distinct basis vectors, so $Q^\top Q = I$.

Products and inverses of orthogonal matrices are orthogonal (Exercise B7), which is why they form a *group*, $O(n)$; LA-12 looks at its rotation subgroup $SO(3)$.

### 3.8 Gram–Schmidt

Given independent $\mathbf{a}_1, \dots, \mathbf{a}_n$, we want orthonormal $\mathbf{q}_1, \dots, \mathbf{q}_n$ with the *same span at every stage*: $\operatorname{span}\{\mathbf{q}_1, \dots, \mathbf{q}_j\} = \operatorname{span}\{\mathbf{a}_1, \dots, \mathbf{a}_j\}$ for each $j$.

▶ **Try it first.** $\mathbf{q}_1$ is forced: $\mathbf{a}_1/\|\mathbf{a}_1\|$. For $\mathbf{q}_2$, take $\mathbf{a}_2$ and remove the part along $\mathbf{q}_1$ (what formula from §3.4 removes it?). Normalize. For $\mathbf{q}_3$, remove the parts along $\mathbf{q}_1$ and $\mathbf{q}_2$. Why can you remove them one at a time?

**Derivation.**

1. $\mathbf{v}_1 = \mathbf{a}_1$, $\mathbf{q}_1 = \mathbf{v}_1/\|\mathbf{v}_1\|$.
2. Subtract from $\mathbf{a}_2$ its projection onto the line through $\mathbf{q}_1$, by (3.4) with $\mathbf{q}_1^\top\mathbf{q}_1 = 1$: $\mathbf{v}_2 = \mathbf{a}_2 - (\mathbf{q}_1^\top\mathbf{a}_2)\mathbf{q}_1$. Then $\mathbf{q}_1^\top\mathbf{v}_2 = \mathbf{q}_1^\top\mathbf{a}_2 - (\mathbf{q}_1^\top\mathbf{a}_2)(\mathbf{q}_1^\top\mathbf{q}_1) = 0$. And $\mathbf{v}_2 \ne \mathbf{0}$, else $\mathbf{a}_2$ would be a multiple of $\mathbf{q}_1$, hence of $\mathbf{a}_1$, contradicting independence. Set $\mathbf{q}_2 = \mathbf{v}_2/\|\mathbf{v}_2\|$.
3. In general, subtract from $\mathbf{a}_j$ its projection onto $\operatorname{span}\{\mathbf{q}_1, \dots, \mathbf{q}_{j-1}\}$. By (3.12) that projection is $\sum_{i<j}(\mathbf{q}_i^\top\mathbf{a}_j)\mathbf{q}_i$, a sum of one-dimensional projections *because the $\mathbf{q}_i$ are orthonormal* (that is what lets you remove the parts one at a time). So
   $$\mathbf{v}_j = \mathbf{a}_j - \sum_{i=1}^{j-1}(\mathbf{q}_i^\top\mathbf{a}_j)\,\mathbf{q}_i, \qquad \mathbf{q}_j = \frac{\mathbf{v}_j}{\|\mathbf{v}_j\|}. \tag{3.14}$$
4. $\mathbf{v}_j$ is the error of a projection onto $\operatorname{span}\{\mathbf{q}_1, \dots, \mathbf{q}_{j-1}\}$, hence orthogonal to each $\mathbf{q}_i$ with $i < j$ (Theorem 3.13). It is nonzero because $\mathbf{a}_j \notin \operatorname{span}\{\mathbf{a}_1, \dots, \mathbf{a}_{j-1}\} = \operatorname{span}\{\mathbf{q}_1, \dots, \mathbf{q}_{j-1}\}$.
5. Spans agree: $\mathbf{q}_j$ is a combination of $\mathbf{a}_j$ and $\mathbf{q}_1, \dots, \mathbf{q}_{j-1}$, so by induction of $\mathbf{a}_1, \dots, \mathbf{a}_j$; and $\mathbf{a}_j = \|\mathbf{v}_j\|\mathbf{q}_j + \sum_{i<j}(\mathbf{q}_i^\top\mathbf{a}_j)\mathbf{q}_i$ is a combination of $\mathbf{q}_1, \dots, \mathbf{q}_j$. ∎

**Theorem 3.19 (Gram–Schmidt).** For independent $\mathbf{a}_1, \dots, \mathbf{a}_n$, the vectors produced by (3.14) are orthonormal and $\operatorname{span}\{\mathbf{q}_1, \dots, \mathbf{q}_j\} = \operatorname{span}\{\mathbf{a}_1, \dots, \mathbf{a}_j\}$ for every $j$. In particular every subspace of $\mathbb{R}^m$ has an orthonormal basis.

■ **What just happened.** Gram–Schmidt is "project, subtract, normalize," repeated. Each new $\mathbf{q}_j$ is the normalized *error* of projecting $\mathbf{a}_j$ onto what you already have. The existence of orthonormal bases is what makes (3.12) usable for *any* subspace.

**Classical versus modified.** Formula (3.14) computes every coefficient $r_{ij} = \mathbf{q}_i^\top\mathbf{a}_j$ from the *original* $\mathbf{a}_j$. That is **classical Gram–Schmidt** (CGS). The **modified** version (MGS) subtracts the projections one at a time and computes each coefficient from the *current*, partly deflated vector:

$$
\mathbf{v} \leftarrow \mathbf{a}_j;\qquad \text{for } i = 1, \dots, j-1:\quad r_{ij} \leftarrow \mathbf{q}_i^\top\mathbf{v},\quad \mathbf{v} \leftarrow \mathbf{v} - r_{ij}\mathbf{q}_i;\qquad \mathbf{q}_j \leftarrow \mathbf{v}/\|\mathbf{v}\|.
$$

In exact arithmetic the two are identical: after subtracting the $\mathbf{q}_1$ component, $\mathbf{q}_2^\top\mathbf{v} = \mathbf{q}_2^\top\mathbf{a}_j - r_{1j}\mathbf{q}_2^\top\mathbf{q}_1 = \mathbf{q}_2^\top\mathbf{a}_j$ because $\mathbf{q}_2^\top\mathbf{q}_1 = 0$ (Exercise C4 does the general case). In floating point they are not identical, because $\mathbf{q}_2^\top\mathbf{q}_1$ is not exactly zero, and the discrepancy compounds. In CGS, an error in $\mathbf{q}_1$'s orthogonality is never seen by later steps; in MGS, each inner product is taken with a vector from which the earlier components have already been removed, so earlier errors get partially cancelled.

The measurable consequence needs three symbols. $\varepsilon \approx 1.1\times10^{-16}$ is the unit roundoff of double precision, the relative error of a single rounding (LA-13 makes this precise; NumPy's `np.finfo(float).eps` is $2\varepsilon$). $\kappa(A)$ is the condition number from LA-9. For a matrix $M$, $\|M\|$ is the 2-norm, the largest stretch $\max_{\mathbf{x}\ne\mathbf{0}}\|M\mathbf{x}\|/\|\mathbf{x}\|$ (LA-9 identifies it with the largest singular value). With $c$ a modest constant that depends only on the matrix size, the computed $Q$ satisfies
$$
\|Q^\top Q - I\| \lesssim c\,\kappa(A)^2\,\varepsilon \ \ \text{(classical, provided } c\,\kappa(A)^2\varepsilon < 1\text{)}, \qquad \|Q^\top Q - I\| \lesssim c\,\kappa(A)\,\varepsilon \ \ \text{(modified, provided } c\,\kappa(A)\varepsilon < 1\text{)}.
$$
The modified bound is Björck's (1967); the classical one is from Giraud, Langou, Rozložník and van den Eshof (*Numerische Mathematik*, 2005). LA-13 states both again. For the proof of the modified bound see Björck's paper or Higham, *Accuracy and Stability of Numerical Algorithms*, §19.8; for the classical bound, which is more recent than Higham's book, see the 2005 paper. These are upper bounds, and they say nothing past their thresholds. Past the threshold, orthogonality is lost entirely and $\|Q^\top Q - I\|$ is of order $1$.

On the $10\times10$ Hilbert matrix ($\kappa \approx 1.6\times10^{13}$) the project gets $\|Q^\top Q - I\| \approx 3.0$ for classical, versus $\approx 1.8\times10^{-4}$ for modified, versus $\approx 10^{-15}$ for Householder (§3.11). Classical is far past its threshold there ($\kappa^2\varepsilon \approx 3\times10^{10}$), which is why the loss is $3$ and not $10^{10}$: the columns are not orthogonal *at all*. Modified is inside its threshold ($\kappa\varepsilon \approx 2\times10^{-3}$) and lands below the bound. The project checks both bounds at every size from 2 to 12 (§7). Never use classical Gram–Schmidt for numbers you care about. §8 says more.

### 3.9 QR factorization

Step 5 of the Gram–Schmidt derivation says each $\mathbf{a}_j$ is a combination of $\mathbf{q}_1, \dots, \mathbf{q}_j$ *only*, with known coefficients:
$$
\mathbf{a}_j = r_{1j}\mathbf{q}_1 + r_{2j}\mathbf{q}_2 + \cdots + r_{jj}\mathbf{q}_j, \qquad r_{ij} = \mathbf{q}_i^\top\mathbf{a}_j \ (i < j), \quad r_{jj} = \|\mathbf{v}_j\|.
$$

▶ **Try it first.** Write these $n$ equations as one matrix equation $A = QR$. What shape is $R$, and which of its entries are forced to be zero?

**Derivation.**

1. The $j$-th column of $QR$ is $Q$ times the $j$-th column of $R$ (LA-1, eq. (3.3)), which is $\sum_i r_{ij}\mathbf{q}_i$.
2. Matching with the display above requires $r_{ij} = 0$ for $i > j$: $\mathbf{a}_j$ has no component along $\mathbf{q}_{j+1}, \dots, \mathbf{q}_n$. Those $\mathbf{q}$'s were built *after* $\mathbf{a}_j$, by subtracting projections, and are orthogonal to $\operatorname{span}\{\mathbf{a}_1, \dots, \mathbf{a}_j\}$. So $R$ is **upper triangular**.
3. Hence
   $$A = QR, \qquad Q\ (m\times n,\ Q^\top Q = I_n), \qquad R\ (n \times n,\ \text{upper triangular},\ r_{jj} = \|\mathbf{v}_j\| > 0). \tag{3.15}$$
   Also $R = Q^\top A$ (multiply (3.15) by $Q^\top$), which is a second way to read off $r_{ij} = \mathbf{q}_i^\top\mathbf{a}_j$. ∎

■ **What just happened.** QR is Gram–Schmidt with the bookkeeping kept. $Q$ holds the orthonormal basis; $R$ records how the original columns are built from it. $R$ is triangular for a *reason*: the $j$-th original vector lives in the span of the first $j$ new ones. This is the second factorization in the course that is a genuine computational tool (LU in LA-2 was the first; the SVD in LA-9 completes the core trio).

**Theorem 3.20 (Existence and uniqueness).** Every $m \times n$ matrix $A$ with independent columns has a factorization (3.15), and it is unique if the diagonal of $R$ is required to be positive.

*Proof.* Existence is Theorem 3.19 plus the derivation above. Uniqueness: suppose $A = Q_1R_1 = Q_2R_2$ with both $R$'s having positive diagonals. Then $Q_2^\top Q_1 R_1 = Q_2^\top Q_2R_2 = R_2$, so $Q_2^\top Q_1 = R_2R_1^{-1}$. The right side is upper triangular with positive diagonal (products and inverses of upper triangular matrices are upper triangular, LA-2, and the diagonal of the product is the product of the diagonals). The left side, call it $M$, satisfies $M^\top M = Q_1^\top Q_2Q_2^\top Q_1$. Now $Q_2Q_2^\top$ is the projection onto $C(Q_2) = C(A) = C(Q_1)$ (Theorem 3.17), so $Q_2Q_2^\top Q_1 = Q_1$ and $M^\top M = Q_1^\top Q_1 = I$: $M$ is orthogonal. An orthogonal upper triangular matrix is diagonal, because its inverse $M^\top$ is lower triangular and also upper triangular (the inverse of an upper triangular matrix is upper triangular); a diagonal orthogonal matrix has entries $\pm1$; with positive diagonal, $M = I$. So $Q_1 = Q_2$ and $R_1 = R_2$. ∎

The "$\pm1$" in that proof is the **sign ambiguity**: without the positivity convention, you may flip the sign of any $\mathbf{q}_i$ together with row $i$ of $R$ and still have a valid QR. NumPy's `np.linalg.qr` does not impose positive diagonals. §8.

**Full versus reduced.** (3.15) is the *reduced* QR (Trefethen & Bau's term; $Q$ is $m\times n$). Extend $\mathbf{q}_1, \dots, \mathbf{q}_n$ to an orthonormal basis of all of $\mathbb{R}^m$ (Gram–Schmidt on $\mathbf{a}_1, \dots, \mathbf{a}_n, \mathbf{e}_1, \dots, \mathbf{e}_m$, skipping dependent vectors) to get a square orthogonal $Q_{\text{full}}$ and an $m\times n$ $R_{\text{full}}$ whose bottom $m - n$ rows are zero. The extra columns of $Q_{\text{full}}$ are an orthonormal basis of $C(A)^\perp = N(A^\top)$. NumPy's default is reduced.

### 3.10 Least squares via QR

▶ **Try it first.** Substitute $A = QR$ into the normal equations (3.6) and simplify using (3.11). Cancel what can be cancelled.

**Derivation.**

1. $A^\top A = (QR)^\top(QR) = R^\top Q^\top QR = R^\top R$, and $A^\top\mathbf{b} = R^\top Q^\top\mathbf{b}$.
2. So (3.6) reads $R^\top R\hat{\mathbf{x}} = R^\top Q^\top\mathbf{b}$.
3. $R$ is triangular with nonzero diagonal, hence invertible, hence so is $R^\top$. Multiply by $(R^\top)^{-1}$:
   $$R\hat{\mathbf{x}} = Q^\top\mathbf{b}. \tag{3.16}$$
4. Solve (3.16) by back substitution (LA-2), $n^2$ operations. The projection is $\mathbf{p} = Q(Q^\top\mathbf{b})$ and the residual is $\mathbf{e} = \mathbf{b} - QQ^\top\mathbf{b}$. ∎

■ **What just happened.** The normal equations contain $A^\top A$; with $A = QR$, the $Q$'s cancel and $A^\top A$ collapses to $R^\top R$, with $R^\top$ then cancelling from both sides. What remains, (3.16), is the normal equations for the *orthonormal* basis: $\hat{\mathbf{x}}$ is read off by expressing $Q^\top\mathbf{b}$ (the coordinates of the shadow, Theorem 3.17) back in terms of the original columns via $R$.

**Why QR and not the normal equations.** Both give the same $\hat{\mathbf{x}}$ in exact arithmetic. In floating point, forming $A^\top A$ is the mistake. The sensitivity of a linear system to rounding is measured by the condition number $\kappa$ (defined in LA-9), and
$$
\kappa(A^\top A) = \kappa(A)^2.
$$
This is stated here and proved in LA-13 (Theorem 3.12; it is one line once you have the SVD: the singular values of $A^\top A$ are the squares of those of $A$). With double precision ($\varepsilon \approx 10^{-16}$), solving a system loses about $\log_{10}\kappa$ digits. A design matrix with $\kappa(A) = 10^8$, which is routine for polynomial fits on badly scaled data, has $\kappa(A^\top A) = 10^{16}$: the normal equations are singular to working precision and the answer is noise. QR works with $A$ itself. On a quadratic fit with $\kappa(A) \approx 6\times10^7$ the project measures relative errors of $4\times10^{-5}$ for the normal equations, $3\times10^{-11}$ for QR by modified Gram–Schmidt (run on $(A\ \ \mathbf{b})$, see the next paragraph), and $7\times10^{-10}$ for Householder QR: the normal equations lose about eleven of sixteen digits, QR five to seven.

Two honest caveats. First, the benefit requires the factorization to be computed *stably*. Householder QR (§3.11, LA-13) is; modified Gram–Schmidt gets most of the way *if* $Q^\top\mathbf{b}$ is obtained by running MGS on the augmented matrix $(A\ \ \mathbf{b})$, so that $\mathbf{b}$ is orthogonalized against the $\mathbf{q}$'s by the same process (the project does this); modified Gram–Schmidt followed by a separately computed $Q^\top\mathbf{b}$ is no better than the normal equations on hard problems (on the fit above it gives $2\times10^{-5}$). Second, when $m$ is enormous and $n$ small, $A^\top A$ is small and cheap, and people use it anyway after scaling the columns; that is a conditioning decision, not a theorem.

### 3.11 Householder reflections (a preview of LA-13)

For $\mathbf{v} \ne \mathbf{0}$ in $\mathbb{R}^m$, the **Householder reflection** is
$$
H = I - 2\frac{\mathbf{v}\mathbf{v}^\top}{\mathbf{v}^\top\mathbf{v}}. \tag{3.17}
$$
Compare with (3.4): $\mathbf{v}\mathbf{v}^\top/\mathbf{v}^\top\mathbf{v}$ is the projection $P_{\mathbf{v}}$ onto the line through $\mathbf{v}$, so $H = I - 2P_{\mathbf{v}}$: go to the foot of the perpendicular and then the same distance again. That is a reflection across the hyperplane $\mathbf{v}^\perp$.

**Theorem 3.21.** $H$ is symmetric and orthogonal; $H\mathbf{v} = -\mathbf{v}$; $H\mathbf{u} = \mathbf{u}$ for every $\mathbf{u} \perp \mathbf{v}$; and $H^2 = I$.

▶ **Try it first.** Transpose (3.17). Then compute $H^\top H$ using $P_{\mathbf{v}}^2 = P_{\mathbf{v}}$.

**Derivation.**

1. $P_{\mathbf{v}} = \mathbf{v}\mathbf{v}^\top/\mathbf{v}^\top\mathbf{v}$ is symmetric and idempotent (Theorem 3.12 with $A = \mathbf{v}$), so $H^\top = I - 2P_{\mathbf{v}}^\top = I - 2P_{\mathbf{v}} = H$.
2. $H^\top H = H^2 = (I - 2P_{\mathbf{v}})^2 = I - 4P_{\mathbf{v}} + 4P_{\mathbf{v}}^2 = I - 4P_{\mathbf{v}} + 4P_{\mathbf{v}} = I$. So $H$ is orthogonal and is its own inverse.
3. $H\mathbf{v} = \mathbf{v} - 2\mathbf{v}(\mathbf{v}^\top\mathbf{v})/(\mathbf{v}^\top\mathbf{v}) = \mathbf{v} - 2\mathbf{v} = -\mathbf{v}$. If $\mathbf{u}^\top\mathbf{v} = 0$ then $H\mathbf{u} = \mathbf{u} - 2\mathbf{v}(\mathbf{v}^\top\mathbf{u})/(\mathbf{v}^\top\mathbf{v}) = \mathbf{u}$. ∎

■ **What just happened.** $H$ flips the component along $\mathbf{v}$ and fixes everything perpendicular to it: a mirror. Being orthogonal, it preserves lengths (Theorem 3.18). The trick that makes it useful: given $\mathbf{x} \ne \mathbf{0}$, let $s = \operatorname{sign}(x_1)$ (take $s = 1$ if $x_1 = 0$), so that $s\,x_1 = |x_1|$, and choose
$$
\mathbf{v} = \mathbf{x} + s\|\mathbf{x}\|\,\mathbf{e}_1.
$$
Then $\mathbf{v}^\top\mathbf{v} = \|\mathbf{x}\|^2 + 2s\|\mathbf{x}\|x_1 + \|\mathbf{x}\|^2 = 2\|\mathbf{x}\|(\|\mathbf{x}\| + |x_1|)$ and $\mathbf{v}^\top\mathbf{x} = \|\mathbf{x}\|^2 + s\|\mathbf{x}\|x_1 = \|\mathbf{x}\|(\|\mathbf{x}\| + |x_1|)$, so $2\mathbf{v}^\top\mathbf{x}/\mathbf{v}^\top\mathbf{v} = 1$ and
$$
H\mathbf{x} = \mathbf{x} - \mathbf{v} = -s\|\mathbf{x}\|\,\mathbf{e}_1.
$$
One reflection sends any vector onto the first axis, zeroing all entries below the first.

Why that sign? The first entry of $\mathbf{v}$ is $x_1 + s\|\mathbf{x}\| = s(|x_1| + \|\mathbf{x}\|)$, a sum of two numbers of the same sign: no cancellation, and $\mathbf{v} \ne \mathbf{0}$ whenever $\mathbf{x} \ne \mathbf{0}$. The other mirror, $\mathbf{v} = \mathbf{x} - s\|\mathbf{x}\|\mathbf{e}_1$, also works in exact arithmetic (it sends $\mathbf{x}$ to $+s\|\mathbf{x}\|\mathbf{e}_1$), but its first entry is $s(|x_1| - \|\mathbf{x}\|)$, a difference of nearly equal numbers whenever $\mathbf{x}$ is close to the $\mathbf{e}_1$ axis. When $\mathbf{x}$ lies on the $\mathbf{e}_1$ axis that entry is exactly $0$, so $\mathbf{v} = \mathbf{0}$ and (3.17) divides by zero. (For $x_1 > 0$ this other mirror is the textbook-looking $\mathbf{v} = \mathbf{x} - \|\mathbf{x}\|\mathbf{e}_1$, which is why it is worth warning about.) The project shows the failure: for $\mathbf{x} = (1, 10^{-9}, 0)$ the computed $\|\mathbf{x}\|$ rounds to $1$, the wrong sign gives $\mathbf{v} = (0, 10^{-9}, 0)$, and $H\mathbf{x} = (1, -10^{-9}, 0)$: the entry that should have become zero only changed sign. The right sign returns $(-1, 0, 0)$. The sign rule is also why `np.linalg.qr` hands back negative diagonal entries in $R$ (§8): $H\mathbf{x} = -s\|\mathbf{x}\|\mathbf{e}_1$ makes the diagonal entry negative whenever the leading entry of the column being reduced is positive.

Apply this to the first column of $A$, then to the second column below the diagonal, and so on: after $n$ reflections, $H_n\cdots H_1A = R$ is upper triangular and $Q = H_1\cdots H_n$ is orthogonal *by construction*. A reflection (3.17) is orthogonal for every nonzero $\mathbf{v}$, so each step applies a map that is orthogonal to working precision however ill-conditioned $A$ is, and the sign rule computes $\mathbf{v}$ accurately, so each step really does zero the entries it targets. That is why Householder QR does not lose orthogonality, and it is LA-13's starting point.

### 3.12 Orthogonal functions (a preview of F2)

Replace $\mathbb{R}^m$ by continuous functions on $[-\pi, \pi]$ and the dot product $\sum_i u_iv_i$ by the integral $\langle f, g\rangle = \int_{-\pi}^{\pi}f(x)g(x)\,dx$ (a sum over a continuum of indices). It is symmetric, linear in each slot, and positive definite: $\langle f, f\rangle = \int f^2 \ge 0$, with equality only for $f = 0$ (if $f(x_0) \ne 0$, continuity keeps $f^2$ positive on a small interval around $x_0$, and the integral is positive). Those are the three properties of the dot product from LA-1 that this module's proofs use; positive definiteness is what makes Theorem 3.3 and step 2 of Theorem 3.11 work. So every theorem about projecting onto the span of *finitely many* functions transfers word for word: the normal equations of §3.4, the closest-point and least-squares Theorems 3.14 and 3.15, Theorem 3.17 with Bessel's inequality, and Gram–Schmidt (Theorem 3.19). Theorem 3.7's statements about complements ($\dim S^\perp = n - k$, $(S^\perp)^\perp = S$, and $\mathbb{R}^n = S \oplus S^\perp$ for every subspace) use the finite dimension of $\mathbb{R}^n$ and do not transfer as stated; LA-8 sets up general inner product spaces, and F2 treats the infinite-dimensional case. (For piecewise continuous functions such as the square wave below, $\int f^2 = 0$ only forces $f = 0$ away from finitely many points. F2 handles this by treating two functions that differ only at finitely many points as the same.)

**Claim.** For positive integers $m \ne n$,
$$
\int_{-\pi}^{\pi}\sin(mx)\sin(nx)\,dx = 0, \qquad \int_{-\pi}^{\pi}\sin^2(nx)\,dx = \pi. \tag{3.18}
$$

▶ **Try it first.** Use $\sin\alpha\sin\beta = \tfrac12[\cos(\alpha - \beta) - \cos(\alpha + \beta)]$.

**Derivation.**

1. $\sin(mx)\sin(nx) = \tfrac12\cos((m-n)x) - \tfrac12\cos((m+n)x)$.
2. For a nonzero integer $k$, $\int_{-\pi}^{\pi}\cos(kx)\,dx = \big[\sin(kx)/k\big]_{-\pi}^{\pi} = 0$, since $\sin(\pm k\pi) = 0$.
3. If $m \ne n$ both $m - n$ and $m + n$ are nonzero integers, so both integrals vanish. If $m = n$, the first term is $\tfrac12\cos 0 = \tfrac12$, integrating to $\pi$, and the second still vanishes. ∎

■ **What just happened.** $\sin x, \sin 2x, \sin 3x, \dots$ are an orthogonal set, each of "length" $\sqrt\pi$. Dividing by $\sqrt\pi$ makes them orthonormal, and (3.12) says the projection of $f$ onto the span of $\sin x, \dots, \sin Nx$ has coefficients $\langle f, \sin nx\rangle/\pi$. Those are the Fourier sine coefficients $b_n = \frac1\pi\int_{-\pi}^{\pi}f(x)\sin(nx)\,dx$. A Fourier series is an orthogonal projection, and the coefficients are dot products.

The coefficients against the *normalized* sines are $\langle f, \sin nx/\sqrt\pi\rangle = \sqrt\pi\,b_n$, so Bessel's inequality (after Theorem 3.17) gives
$$
\pi\sum_{n=1}^{N} b_n^2 \le \int_{-\pi}^{\pi} f^2\,dx \quad\text{for every } N.
$$
Do not upgrade this to Parseval's equality. (3.13) needed a square $Q$, an orthonormal basis of the *whole* space, and the sines alone are not a basis for all functions: the constant $f \equiv 1$ has every $b_n = 0$ (its product with $\sin nx$ is odd) but $\int f^2 = 2\pi$. For *odd* $f$ the sines are enough, and equality holds in the limit $N \to \infty$ (F2 Theorem 3.8). For general $f$, F2 adds $1$ and $\cos nx$ and gets the full Parseval identity, and asks when the projection converges to $f$ pointwise. The project projects the square wave $\operatorname{sign}(x)$, which is odd, onto the first few sines, and you can watch it happen. It also checks the equality: $b_n = 4/(n\pi)$ for odd $n$ and $0$ for even $n$, so $\pi\sum b_n^2 = \frac{16}{\pi}\big(1 + \frac19 + \frac1{25} + \cdots\big) = \frac{16}{\pi}\cdot\frac{\pi^2}{8} = 2\pi = \int_{-\pi}^{\pi}\operatorname{sign}(x)^2\,dx$.

## 4. Worked examples

### Example 4.1: Projecting onto a plane in $\mathbb{R}^3$

Project $\mathbf{b} = (1, 2, 3)^\top$ onto the plane spanned by $\mathbf{a}_1 = (1, 0, 1)^\top$ and $\mathbf{a}_2 = (1, 1, 0)^\top$.

$A = \begin{pmatrix}1 & 1\\0 & 1\\1 & 0\end{pmatrix}$. The Gram matrix and right side of (3.6):
$$
A^\top A = \begin{pmatrix}\mathbf{a}_1^\top\mathbf{a}_1 & \mathbf{a}_1^\top\mathbf{a}_2\\ \mathbf{a}_2^\top\mathbf{a}_1 & \mathbf{a}_2^\top\mathbf{a}_2\end{pmatrix} = \begin{pmatrix}2 & 1\\1 & 2\end{pmatrix}, \qquad A^\top\mathbf{b} = \begin{pmatrix}1 + 0 + 3\\ 1 + 2 + 0\end{pmatrix} = \begin{pmatrix}4\\3\end{pmatrix}.
$$
Solve $\begin{pmatrix}2 & 1\\1 & 2\end{pmatrix}\hat{\mathbf{x}} = \begin{pmatrix}4\\3\end{pmatrix}$: the inverse is $\frac13\begin{pmatrix}2 & -1\\-1 & 2\end{pmatrix}$, so $\hat{\mathbf{x}} = \frac13\begin{pmatrix}8 - 3\\-4 + 6\end{pmatrix} = \begin{pmatrix}5/3\\2/3\end{pmatrix}$.

Projection: $\mathbf{p} = \tfrac53\mathbf{a}_1 + \tfrac23\mathbf{a}_2 = \big(\tfrac53 + \tfrac23,\ \tfrac23,\ \tfrac53\big)^\top = \big(\tfrac73, \tfrac23, \tfrac53\big)^\top$.

Error: $\mathbf{e} = \mathbf{b} - \mathbf{p} = \big(1 - \tfrac73,\ 2 - \tfrac23,\ 3 - \tfrac53\big)^\top = \big(-\tfrac43, \tfrac43, \tfrac43\big)^\top$.

Check perpendicularity (3.5): $\mathbf{a}_1^\top\mathbf{e} = -\tfrac43 + 0 + \tfrac43 = 0$ ✓, $\mathbf{a}_2^\top\mathbf{e} = -\tfrac43 + \tfrac43 + 0 = 0$ ✓. Indeed $\mathbf{e} = \tfrac43(-1, 1, 1)^\top$, and $(-1, 1, 1)^\top$ is the normal to the plane ($= \mathbf{a}_1 \times \mathbf{a}_2$). Pythagoras: $\|\mathbf{b}\|^2 = 14$, $\|\mathbf{p}\|^2 = \tfrac{49 + 4 + 25}{9} = \tfrac{26}{3}$, $\|\mathbf{e}\|^2 = \tfrac{3\cdot16}{9} = \tfrac{16}{3}$, and $\tfrac{26}{3} + \tfrac{16}{3} = 14$ ✓.

The projection matrix (3.8): $A(A^\top A)^{-1} = \tfrac13\begin{pmatrix}1 & 1\\0 & 1\\1 & 0\end{pmatrix}\begin{pmatrix}2 & -1\\-1 & 2\end{pmatrix} = \tfrac13\begin{pmatrix}1 & 1\\-1 & 2\\2 & -1\end{pmatrix}$, then times $A^\top = \begin{pmatrix}1 & 0 & 1\\1 & 1 & 0\end{pmatrix}$:
$$
P = \frac13\begin{pmatrix}2 & 1 & 1\\1 & 2 & -1\\1 & -1 & 2\end{pmatrix}.
$$
Symmetric ✓. $P\mathbf{b} = \tfrac13(2 + 2 + 3,\ 1 + 4 - 3,\ 1 - 2 + 6)^\top = \tfrac13(7, 2, 5)^\top = \mathbf{p}$ ✓. And $I - P = \tfrac13\begin{pmatrix}1 & -1 & -1\\-1 & 1 & 1\\-1 & 1 & 1\end{pmatrix} = \mathbf{n}\mathbf{n}^\top/\mathbf{n}^\top\mathbf{n}$ with $\mathbf{n} = (-1, 1, 1)^\top$: the projection onto the normal line, as Theorem 3.13 says. The project checks $P^2 = P$ and draws the picture.

### Example 4.2: The least-squares line

Fit $b = C + Dt$ to the points $(0, 1), (1, 2), (2, 2), (3, 4)$. No line passes through all four (the slopes between consecutive points are $1, 0, 2$).

$$
A = \begin{pmatrix}1 & 0\\1 & 1\\1 & 2\\1 & 3\end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix}1\\2\\2\\4\end{pmatrix}, \quad
A^\top A = \begin{pmatrix}4 & 6\\6 & 14\end{pmatrix}, \quad A^\top\mathbf{b} = \begin{pmatrix}1 + 2 + 2 + 4\\0 + 2 + 4 + 12\end{pmatrix} = \begin{pmatrix}9\\18\end{pmatrix}.
$$
Normal equations: $4C + 6D = 9$ and $6C + 14D = 18$. Multiply the first by 3 and the second by 2: $12C + 18D = 27$, $12C + 28D = 36$; subtract: $10D = 9$, so $D = 0.9$; then $4C = 9 - 5.4 = 3.6$, $C = 0.9$. The line is $b = 0.9 + 0.9t$.

Fitted values $\mathbf{p} = A\hat{\mathbf{x}} = (0.9, 1.8, 2.7, 3.6)^\top$. Residual $\mathbf{r} = \mathbf{b} - \mathbf{p} = (0.1, 0.2, -0.7, 0.4)^\top$.

Check (3.5): $\mathbf{r}$ against column 1: $0.1 + 0.2 - 0.7 + 0.4 = 0$ ✓ (the residuals of any fit with a constant term sum to zero, Exercise B5). Against column 2: $0 + 0.2 - 1.4 + 1.2 = 0$ ✓. Minimum value $\|\mathbf{r}\|^2 = 0.01 + 0.04 + 0.49 + 0.16 = 0.70$; and $\|\mathbf{b}\|^2 - \|\mathbf{p}\|^2 = 25 - 24.3 = 0.7$ ✓ (Theorem 3.15).

### Example 4.3: Gram–Schmidt in $\mathbb{R}^4$

Orthonormalize $\mathbf{a}_1 = (1, 1, 1, 1)^\top$, $\mathbf{a}_2 = (1, 3, 3, 1)^\top$, $\mathbf{a}_3 = (1, 1, 1, 4)^\top$.

*Step 1.* $\|\mathbf{a}_1\| = 2$, so $\mathbf{q}_1 = \tfrac12(1, 1, 1, 1)^\top$.

*Step 2.* $r_{12} = \mathbf{q}_1^\top\mathbf{a}_2 = \tfrac12(1 + 3 + 3 + 1) = 4$. $\mathbf{v}_2 = \mathbf{a}_2 - 4\mathbf{q}_1 = (1, 3, 3, 1)^\top - (2, 2, 2, 2)^\top = (-1, 1, 1, -1)^\top$. $r_{22} = \|\mathbf{v}_2\| = 2$, $\mathbf{q}_2 = \tfrac12(-1, 1, 1, -1)^\top$. Check $\mathbf{q}_1^\top\mathbf{q}_2 = \tfrac14(-1 + 1 + 1 - 1) = 0$ ✓.

*Step 3.* $r_{13} = \mathbf{q}_1^\top\mathbf{a}_3 = \tfrac12(1 + 1 + 1 + 4) = \tfrac72$. $r_{23} = \mathbf{q}_2^\top\mathbf{a}_3 = \tfrac12(-1 + 1 + 1 - 4) = -\tfrac32$.
$$
\mathbf{v}_3 = \mathbf{a}_3 - \tfrac72\mathbf{q}_1 + \tfrac32\mathbf{q}_2 = \begin{pmatrix}1\\1\\1\\4\end{pmatrix} - \begin{pmatrix}7/4\\7/4\\7/4\\7/4\end{pmatrix} + \begin{pmatrix}-3/4\\3/4\\3/4\\-3/4\end{pmatrix} = \begin{pmatrix}-3/2\\0\\0\\3/2\end{pmatrix}.
$$
$r_{33} = \|\mathbf{v}_3\| = \sqrt{9/4 + 9/4} = 3/\sqrt2$, $\mathbf{q}_3 = \tfrac{1}{\sqrt2}(-1, 0, 0, 1)^\top$. Checks: $\mathbf{q}_1^\top\mathbf{q}_3 = 0$ ✓, $\mathbf{q}_2^\top\mathbf{q}_3 = \tfrac{1}{2\sqrt2}(1 + 0 + 0 - 1) = 0$ ✓.

Result:
$$
Q = \begin{pmatrix}1/2 & -1/2 & -1/\sqrt2\\ 1/2 & 1/2 & 0\\ 1/2 & 1/2 & 0\\ 1/2 & -1/2 & 1/\sqrt2\end{pmatrix}, \qquad
R = \begin{pmatrix}2 & 4 & 7/2\\ 0 & 2 & -3/2\\ 0 & 0 & 3/\sqrt2\end{pmatrix}.
$$
Verify column 3 of $QR$: $\tfrac72\mathbf{q}_1 - \tfrac32\mathbf{q}_2 + \tfrac{3}{\sqrt2}\mathbf{q}_3 = (\tfrac74, \tfrac74, \tfrac74, \tfrac74) + (\tfrac34, -\tfrac34, -\tfrac34, \tfrac34) + (-\tfrac32, 0, 0, \tfrac32) = (1, 1, 1, 4) = \mathbf{a}_3$ ✓.

### Example 4.4: QR of a $3\times2$, and least squares via $R\hat{\mathbf{x}} = Q^\top\mathbf{b}$

**(a)** Take the $A$ of Example 4.1. $\mathbf{q}_1 = \mathbf{a}_1/\sqrt2 = \tfrac{1}{\sqrt2}(1, 0, 1)^\top$, $r_{11} = \sqrt2$. $r_{12} = \mathbf{q}_1^\top\mathbf{a}_2 = \tfrac{1}{\sqrt2}$. $\mathbf{v}_2 = \mathbf{a}_2 - \tfrac{1}{\sqrt2}\mathbf{q}_1 = (1, 1, 0)^\top - \tfrac12(1, 0, 1)^\top = (\tfrac12, 1, -\tfrac12)^\top$, $r_{22} = \sqrt{\tfrac14 + 1 + \tfrac14} = \sqrt{3/2} = \tfrac{\sqrt6}{2}$, $\mathbf{q}_2 = \tfrac{1}{\sqrt6}(1, 2, -1)^\top$.
$$
Q = \begin{pmatrix}1/\sqrt2 & 1/\sqrt6\\ 0 & 2/\sqrt6\\ 1/\sqrt2 & -1/\sqrt6\end{pmatrix}, \qquad R = \begin{pmatrix}\sqrt2 & 1/\sqrt2\\ 0 & \sqrt6/2\end{pmatrix}.
$$
Now the same $\mathbf{b} = (1, 2, 3)^\top$: $Q^\top\mathbf{b} = \big(\tfrac{1 + 3}{\sqrt2},\ \tfrac{1 + 4 - 3}{\sqrt6}\big)^\top = \big(2\sqrt2,\ \tfrac{2}{\sqrt6}\big)^\top$. Back substitution in (3.16): second row, $\tfrac{\sqrt6}{2}\hat{x}_2 = \tfrac{2}{\sqrt6}$, so $\hat{x}_2 = \tfrac{4}{6} = \tfrac23$. First row, $\sqrt2\hat{x}_1 + \tfrac{1}{\sqrt2}\cdot\tfrac23 = 2\sqrt2$, so $\sqrt2\hat{x}_1 = 2\sqrt2 - \tfrac{\sqrt2}{3} = \tfrac{5\sqrt2}{3}$, $\hat{x}_1 = \tfrac53$. Same $\hat{\mathbf{x}}$ as Example 4.1 ✓. The projection via (3.12): $\mathbf{p} = 2\sqrt2\,\mathbf{q}_1 + \tfrac{2}{\sqrt6}\mathbf{q}_2 = (2, 0, 2)^\top + \tfrac13(1, 2, -1)^\top = (\tfrac73, \tfrac23, \tfrac53)^\top$ ✓, and $QQ^\top$ multiplies out to the same $P$ as before.

**(b)** The line of Example 4.2. $\mathbf{a}_1 = (1, 1, 1, 1)^\top$, $\mathbf{q}_1 = \tfrac12\mathbf{a}_1$, $r_{11} = 2$. $r_{12} = \mathbf{q}_1^\top(0, 1, 2, 3)^\top = 3$. $\mathbf{v}_2 = (0, 1, 2, 3)^\top - \tfrac32(1, 1, 1, 1)^\top = (-\tfrac32, -\tfrac12, \tfrac12, \tfrac32)^\top$, $r_{22} = \sqrt{\tfrac{9 + 1 + 1 + 9}{4}} = \sqrt5$, $\mathbf{q}_2 = \tfrac{1}{2\sqrt5}(-3, -1, 1, 3)^\top$. So $R = \begin{pmatrix}2 & 3\\0 & \sqrt5\end{pmatrix}$ and $Q^\top\mathbf{b} = \big(\tfrac92,\ \tfrac{-3 - 2 + 2 + 12}{2\sqrt5}\big)^\top = \big(\tfrac92, \tfrac{9}{2\sqrt5}\big)^\top$. Back substitution: $\sqrt5 D = \tfrac{9}{2\sqrt5} \Rightarrow D = \tfrac{9}{10}$; $2C + 3(0.9) = 4.5 \Rightarrow C = 0.9$. Same line ✓, and notice the second column of $Q$ is the *centered* time $t - \bar t$, normalized; Exercise B5 explains why that is not a coincidence.

### Example 4.5: Calibrating a pressure transducer (an ME example)

A transducer's output voltage $V$ is read at six known pressures. Write $u$ for the applied pressure (the letters $P$ and $\mathbf{p}$ are taken by the projection): $u = 0, 1, 2, 3, 4, 5$ bar.

| $u$ (bar) | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| $V$ (V) | 0.45 | 1.39 | 2.22 | 3.04 | 3.95 | 5.05 |

The datasheet says the response is nearly linear with a small quadratic term, so fit $V = c_0 + c_1u + c_2u^2$. The unknowns are $\hat{\mathbf{x}} = (c_0, c_1, c_2)^\top$, the data vector is $\mathbf{V}$ (it plays the role of $\mathbf{b}$), and the design matrix has columns $\mathbf{1}, \mathbf{u}, \mathbf{u}^2$. With $\sum u = 15$, $\sum u^2 = 55$, $\sum u^3 = 225$, $\sum u^4 = 979$:
$$
A^\top A = \begin{pmatrix}6 & 15 & 55\\15 & 55 & 225\\55 & 225 & 979\end{pmatrix}, \qquad
A^\top\mathbf{V} = \begin{pmatrix}\sum V\\ \sum uV\\ \sum u^2V\end{pmatrix} = \begin{pmatrix}16.10\\56.00\\227.08\end{pmatrix}.
$$
(Check the middle entry: $0 + 1.39 + 4.44 + 9.12 + 15.80 + 25.25 = 56.00$.) Eliminate (LA-2): row 2 minus $2.5\times$ row 1 gives $17.5c_1 + 87.5c_2 = 15.75$, i.e. $c_1 + 5c_2 = 0.9$. Row 3 minus $\tfrac{11}{3}\times$ row 2 gives $\tfrac{70}{3}c_1 + 154c_2 = 21.7467$, i.e. $70c_1 + 462c_2 = 65.24$. Substitute $c_1 = 0.9 - 5c_2$: $63 - 350c_2 + 462c_2 = 65.24$, so $112c_2 = 2.24$, $c_2 = 0.02$; then $c_1 = 0.8$; then $6c_0 = 16.1 - 15(0.8) - 55(0.02) = 3$, $c_0 = 0.5$.

Calibration: $V = 0.5 + 0.8u + 0.02u^2$ (volts, bar). Residuals $\mathbf{r} = \mathbf{V} - A\hat{\mathbf{x}} = (-0.05, 0.07, 0.04, -0.04, -0.07, 0.05)^\top$ V, $\|\mathbf{r}\|^2 = 0.018$, RMS $0.055$ V. Check (3.5): $\sum r_i = 0$, $\sum u_ir_i = 0 + 0.07 + 0.08 - 0.12 - 0.28 + 0.25 = 0$, $\sum u_i^2r_i = 0 + 0.07 + 0.16 - 0.36 - 1.12 + 1.25 = 0$ ✓.

To *use* the calibration you invert it: given a reading $V$, solve $0.02u^2 + 0.8u + (0.5 - V) = 0$ for $u$. Two engineering remarks. These readings were built so the arithmetic is exact (the residual is a multiple of $(-5, 7, 4, -4, -7, 5)$, which is orthogonal to all three columns by construction); real readings are not this tidy, and Part 4 of the project also fits a noisy synthetic quadratic. And if the pressures had been in pascals ($u \sim 10^5$), the columns $\mathbf{1}, \mathbf{u}, \mathbf{u}^2$ would range over $1, 10^5, 10^{10}$ and $A^\top A$ would be hopelessly ill-conditioned: *center and scale* the input ($u \to (u - \bar u)/s$ with $\bar u$ the mean and $s$ a typical spread) before fitting. That is §3.10's $\kappa(A)^2$ in practice.

## 5. Why it matters and where it connects

- **LA-8.** The spectral theorem says a symmetric matrix has an *orthonormal* basis of eigenvectors, $A = Q\Lambda Q^\top$ with $Q$ orthogonal (Theorem 3.18's matrices). The identity $\mathbf{x}^\top A^\top A\mathbf{x} = \|A\mathbf{x}\|^2$ from Theorem 3.11 is the definition of positive semidefinite. Inner product spaces generalize §3.12.
- **LA-9.** The SVD $A = U\Sigma V^\top$ has orthogonal $U, V$; the pseudoinverse $A^+$ equals $(A^\top A)^{-1}A^\top$ when the columns are independent, so (3.7) is $\hat{\mathbf{x}} = A^+\mathbf{b}$, and it still gives the minimum-norm least-squares solution when they are not. $\kappa(A)$ is defined there; $\kappa(A^\top A) = \kappa(A)^2$ follows in one line.
- **LA-13.** Householder QR (§3.11) is the stable way to compute (3.15), and LA-13 explains why (the full backward-stability proof is in Trefethen & Bau, lecture 16). LA-13 Theorem 3.12 proves $\kappa(A^\top A) = \kappa(A)^2$, the root of the normal equations' trouble, and LA-13 states the Gram–Schmidt bounds of §3.8 again (proofs: Björck, 1967, and Higham, *Accuracy and Stability of Numerical Algorithms*, §19.8, for modified; Giraud, Langou, Rozložník and van den Eshof, 2005, for classical). QR iteration for eigenvalues is $A \to RQ$ repeated.
- **LA-14.** Equation (3.10), $\nabla\|A\mathbf{x} - \mathbf{b}\|^2 = 2A^\top(A\mathbf{x} - \mathbf{b})$, is the first entry in the matrix-calculus table; the Hessian $2A^\top A$ is the second.
- **F2.** A Fourier series is the projection (3.12) onto $\{1, \cos nx, \sin nx\}$ with the inner product of §3.12; the coefficient formulas are "dot with $\mathbf{q}_i$"; Bessel's inequality (after Theorem 3.17) bounds the coefficients (F2 Theorem 3.6), and Parseval is (3.13) once F2 takes the trigonometric system to be complete (F2 Theorem 3.8). EM-3's Legendre polynomials are Gram–Schmidt applied to $1, x, x^2, \dots$ on $[-1, 1]$ (Exercise C2).
- **ST-3.** Regression is Example 4.2 with a probability model on the residuals. The "hat matrix" of statistics is $P$ of (3.8); fitted values are $P\mathbf{y}$, residuals are $(I - P)\mathbf{y}$, and their orthogonality (Theorem 3.13) is why the sums of squares decompose (ANOVA). ST-3 also supplies §1's postponed third reason for squares (Gauss, 1809: with normally distributed errors, least squares is the maximum-likelihood estimate) and its distribution-free successor, the Gauss–Markov theorem (with errors of mean zero, equal variance and no correlation, least squares is the best linear unbiased estimator). Weighted least squares (Exercise B8) is the first step toward generalized least squares.
- **ML-1 to ML-3.** Linear regression is the first learning algorithm; (3.10) is what gradient descent follows when $A^\top A$ is too big to form; Gauss–Newton (nonlinear least squares, §1's outer loop) linearizes and solves (3.6) repeatedly.
- **CM-6.** Normal modes are orthogonal with respect to the *mass-weighted* inner product $\mathbf{u}^\top M\mathbf{v}$; projecting an initial condition onto them is (3.12) with that inner product.
- **QM-4, QC-1.** An orthonormal eigenbasis $\{|\phi_i\rangle\}$ expands any state as $|\psi\rangle = \sum_i\langle\phi_i|\psi\rangle|\phi_i\rangle$, which is (3.13) with complex entries; the probabilities $|\langle\phi_i|\psi\rangle|^2$ sum to 1 by Parseval. Quantum gates are unitary, the complex version of Theorem 3.18: they preserve lengths, so probabilities stay normalized.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: orthogonality checks, a complement, a line projection and its matrix, a three-point least-squares line, an orthogonal matrix, a two-vector Gram–Schmidt. B1–B8 are the module: B2 (split a vector into row-space and null-space parts), B4 (the projection onto a plane two ways), B6 (the gradient derivation, entry by entry), and B8 (weighted least squares in `project.py`) are essential. C1–C4 go further: uniqueness of QR, Legendre polynomials by Gram–Schmidt, the log-decrement fit for a damped oscillator, and why modified Gram–Schmidt beats classical.

## 7. Coding project

`project.py` builds everything in §3 from scratch and checks it against NumPy and SciPy:

1. **Four subspaces at right angles.** Builds a rank-2 $6\times4$ matrix, gets orthonormal bases of $C(A), N(A^\top), C(A^\top), N(A)$, and asserts Theorem 3.9 (all cross dot products zero), Theorem 3.7 (dimensions add up; every $\mathbf{x}$ splits into row-space and null-space parts), that $A$ restricted to the row space is one-to-one (§3.3), and the Fredholm alternative.
2. **Projection.** `project_onto_subspace(A, b)` solves (3.6) and returns $\mathbf{p}, \hat{\mathbf{x}}, \mathbf{e}, P$. Example 4.1 is reproduced to the digit; $P^2 = P$, $P^\top = P$, $I - P = \mathbf{n}\mathbf{n}^\top/\mathbf{n}^\top\mathbf{n}$, and the closest-point property are asserted. A 3D figure shows the plane, $\mathbf{b}$, $\mathbf{p}$, and $\mathbf{e}$.
3. **Least squares.** `least_squares_normal` reproduces Example 4.2 and matches `np.linalg.lstsq`; the gradient formula (3.10) is checked against central finite differences of (3.9) at random points, and symbolically with SymPy (gradient and Hessian of $\|A\mathbf{x} - \mathbf{b}\|^2$ for a symbolic $3\times2$ system).
4. **Polynomial fits.** Example 4.5's calibration, then a noisy quadratic with a residual plot, including what the residuals look like when you fit the wrong model (a line): a visible pattern instead of noise.
5. **Gram–Schmidt.** `gram_schmidt_classical` and `gram_schmidt_modified` reproduce Example 4.3. On Hilbert matrices $H_{ij} = 1/(i + j - 1)$ of size 2 to 12, $\|Q^\top Q - I\|$ (2-norm, as in §3.8 and LA-13) is computed for both and for Householder; at $n = 10$ the three are about $3.0$, $1.8\times10^{-4}$, and $10^{-15}$. At every size the script asserts the §3.8 bounds: modified below $10\,\kappa\varepsilon$, classical below $10\,\kappa^2\varepsilon$ wherever $\kappa^2\varepsilon < 1$.
6. **QR.** `qr_decompose` (modified Gram–Schmidt, positive diagonal) and `least_squares_qr` (which runs MGS on $(A\ \ \mathbf{b})$ so that $Q^\top\mathbf{b}$ and $\|\mathbf{e}\|$ fall out of the factorization) reproduce Example 4.4 and match `np.linalg.qr` up to column signs. The identity $\kappa(A^\top A) = \kappa(A)^2$ is checked numerically, then an ill-conditioned quadratic fit ($t \in [100, 105]$) shows relative errors of $4\times10^{-5}$ for the normal equations, $2\times10^{-5}$ for naive MGS (asserted to be within a factor of ten of the normal equations), $3\times10^{-11}$ for augmented MGS, and $7\times10^{-10}$ for Householder.
7. **Orthogonal matrices and Householder.** Rotations, reflections, permutations, and random orthogonal matrices preserve lengths and dot products; `householder(v)` is symmetric and orthogonal, and with the sign rule of §3.11 (`householder_vector(x)`) sends $\mathbf{x}$ to $-\operatorname{sign}(x_1)\|\mathbf{x}\|\mathbf{e}_1$; the opposite sign is shown failing on $\mathbf{x} = (1, 10^{-9}, 0)$; a 2D figure reflects a shape across the mirror line $\mathbf{v}^\perp$.
8. **Orthogonal functions.** SymPy verifies (3.18) for $m, n \le 4$ and computes the sine coefficients of a square wave ($4/n\pi$ for odd $n$), checks Bessel's inequality $\pi\sum_{n\le N}b_n^2 \le 2\pi$ and, since the square wave is odd, the limit $2\pi$ (exactly, with SymPy's $\sum_{\text{odd}}1/n^2 = \pi^2/8$); a figure shows the projection onto the first 1, 3, 7 sines.

**Build it yourself first:** `project_onto_subspace`, `gram_schmidt_classical`, and `least_squares_qr` (the naive version, $R\hat{\mathbf{x}} = Q^\top\mathbf{b}$ with your own $Q, R$). Then compare. The "Your turn" cells add weighted least squares (with Exercise B8's weighted calibration checked), the damped-oscillator envelope, centering, and a second pass of classical Gram–Schmidt.

## 8. Where you'll struggle

- **"Where did the transpose come from?"** Everyone asks this at (3.5), and the usual answer ("multiply both sides by $A^\top$") is not an explanation, it is a recipe. The explanation is Lemma 3.8: the error must be orthogonal to the *whole* column space, which is the same as orthogonal to each column, which is $n$ dot products $\mathbf{a}_j^\top\mathbf{e} = 0$, and the matrix whose rows are $\mathbf{a}_j^\top$ is $A^\top$. The transpose is "every column, dotted with," in one symbol. If you can say that sentence out loud, you own the normal equations. If you memorized "$A^\top$ on both sides," you will forget it in six weeks.
- **Three objects, one letter.** $\mathbf{p}$ (the shadow, in $\mathbb{R}^m$), $P$ (the matrix that casts shadows, $m\times m$, depends only on the subspace), and $\hat{\mathbf{x}}$ (the coefficients, in $\mathbb{R}^n$, *not* in the subspace at all). Students write $P = (A^\top A)^{-1}A^\top$ (that is the matrix giving $\hat{\mathbf{x}}$, $n\times m$) or report $\hat{\mathbf{x}}$ when asked for the projection. Write the shapes above every symbol until it stops happening. Also: $A\hat{\mathbf{x}} = \mathbf{p}$, never $\mathbf{b}$. Least squares does not solve $A\mathbf{x} = \mathbf{b}$; it solves the nearest solvable system.
- **$Q$ must have *orthonormal* columns, not merely orthogonal ones.** The formula $P = QQ^\top$ silently uses $\mathbf{q}_i^\top\mathbf{q}_i = 1$. If your columns are orthogonal but have lengths $\ell_i$, the projection is $\sum_i(\mathbf{q}_i^\top\mathbf{b}/\ell_i^2)\mathbf{q}_i$, which is (3.8) with the diagonal Gram matrix $\operatorname{diag}(\ell_i^2)$ inverted. Forgetting to normalize in Gram–Schmidt is the most common bug in the project; the symptom is $Q^\top Q$ diagonal but not $I$.
- **Sign ambiguity in QR.** Theorem 3.20 makes QR unique *only after* you fix the signs of the diagonal of $R$. `np.linalg.qr` returns whatever signs its Householder reflections produce: by the sign rule of §3.11, a diagonal entry comes out negative whenever the leading entry of the column being reduced is positive (for a design matrix whose first column is all ones, $r_{11} = -\sqrt m$). Your Gram–Schmidt $Q$ will differ from NumPy's by $\pm1$ per column, and `np.allclose(Q, Q_np)` will fail while both are right. Compare $Q D$ with $D = \operatorname{diag}(\operatorname{sign}(\operatorname{diag} R_{\text{np}}))$, or compare $QR$ products, or compare $QQ^\top$ (which is sign-invariant). The project shows the fix.
- **You will be tempted to use classical Gram–Schmidt anyway.** It is shorter, it vectorizes (all the $r_{ij}$ for column $j$ are one product $Q_{:,1:j}^\top\mathbf{a}_j$), and on the $4\times3$ textbook example it gives the same answer to fourteen digits. Then you orthogonalize 50 columns of real data and $\|Q^\top Q - I\| = 0.3$ with no warning. Modified Gram–Schmidt costs the same number of flops. Use it, or use Householder (`np.linalg.qr`). And know that even MGS is not the end: for least squares, a separately computed $Q^\top\mathbf{b}$ throws the accuracy away again (project, Part 6); orthogonalize $\mathbf{b}$ along with the columns, or use Householder.
- **Never form the inverse in (3.7).** $(A^\top A)^{-1}$ is a theorem, not an algorithm. Solve (3.6) with `np.linalg.solve`, or better, use (3.16). And if you can avoid forming $A^\top A$ at all, do: $\kappa^2$ is a steep price.
- **Dependent columns.** If $A^\top A$ comes out singular, your columns are dependent (Theorem 3.11), typically because you included both a constant column and a set of indicator columns that sum to it, or $t$ and $2t$. The projection still exists; the recipe does not. Remove the redundancy (or wait for LA-9's pseudoinverse).
- **Exercise B6 is harder than it looks** if you have never differentiated a double sum. Do it for $n = 2$ with the four terms written out before doing it in general. Keeping the summation letters apart from the differentiation index, and applying the product rule to $x_ix_l$ (step 3 of the calculus derivation), is the entire difficulty.

## 9. Check yourself

1. Define $S^\perp$. If $S$ is a plane through the origin in $\mathbb{R}^5$, what is $\dim S^\perp$, and why?
2. State the two orthogonality relations among the four fundamental subspaces, and say in which space ($\mathbb{R}^n$ or $\mathbb{R}^m$) each lives.
3. Derive the normal equations from "the error is perpendicular to every column" in three lines. Where exactly does $A^\top$ enter?
4. What are the shapes of $\hat{\mathbf{x}}$, $\mathbf{p}$, $P$, and $A^\top A$ for $A$ of size $m\times n$? Which one depends on $\mathbf{b}$?
5. Why is $A^\top A$ invertible when the columns of $A$ are independent? (The proof is one identity: name it.)
6. Give two properties that together characterize an orthogonal projection matrix, and say what each means geometrically.
7. What does the gradient of $\|A\mathbf{x} - \mathbf{b}\|^2$ equal, and why does setting it to zero give the same equations as the geometric argument?
8. If $Q$ has orthonormal columns, what is the projection matrix onto $C(Q)$, and what are the coefficients of the projection of $\mathbf{b}$?
9. Prove in two lines that an orthogonal matrix preserves lengths. Then show that $H = I - 2\mathbf{v}\mathbf{v}^\top/\mathbf{v}^\top\mathbf{v}$ satisfies $H^2 = I$. What does $H$ do to $\mathbf{v}$?
10. Describe Gram–Schmidt in one sentence. Why is $R$ in $A = QR$ upper triangular?
11. Derive $R\hat{\mathbf{x}} = Q^\top\mathbf{b}$ from the normal equations. Why is it numerically preferable, in terms of $\kappa$?
12. Why is $\int_{-\pi}^{\pi}\sin 2x\sin 3x\,dx = 0$, and what does that have to do with this module?

## 10. Pair with

- **Strang, *Introduction to Linear Algebra*, ch. 4** (§4.1 orthogonality of the four subspaces, §4.2 projections, §4.3 least squares, §4.4 orthonormal bases and Gram–Schmidt). Strang's pictures of the big picture with right angles and his "three ways to see $A^\top A\hat{\mathbf{x}} = A^\top\mathbf{b}$" are the model for §3.3–3.6; his worked line fit is the classic.
- **MIT 18.06, lectures 14–17 (Strang).** 14: orthogonal vectors and subspaces. 15: projections onto subspaces. 16: projection matrices and least squares. 17: orthogonal matrices and Gram–Schmidt. Watch 15 and 16 after §3.4–3.6; seeing him derive $P$ at the board, slowly, is worth the hour.
- **Trefethen & Bau, *Numerical Linear Algebra*, lectures 6–9.** 6: projectors (the cleanest treatment of $P^2 = P$, including oblique ones). 7: QR factorization, full and reduced, with the "why is $R$ triangular" argument. 8: classical and modified Gram–Schmidt, with the operation count. 9: the classical-versus-modified and loss-of-orthogonality experiments that the project reproduces. These are the standard reference for §3.8–3.10, written for people who will compute.
- **Optional:** Trefethen & Bau lecture 10 (Householder triangularization) and lecture 11 (least squares, three algorithms compared) if §3.10–3.11 left you wanting the full story before LA-13.
