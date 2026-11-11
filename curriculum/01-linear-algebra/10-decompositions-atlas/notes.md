# LA-10: Decompositions Atlas: LU, Cholesky, QR, Schur, Spectral, SVD, Polar; Which to Use When

**Time:** ~6 hours · **Prerequisites:** LA-9 (and through it LA-2, LA-5, LA-7, LA-8) · **Unlocks:** LA-13

> **Conventions in this module.** Matrices are uppercase, vectors are bold lowercase columns. $A^*$ is the conjugate transpose $\bar A^{\top}$ (the physics modules write $A^\dagger$); for real $A$, $A^* = A^\top$. $Q$ is **unitary** if $Q^*Q = I$ and **orthogonal** if it is real and unitary. $A$ is **symmetric positive definite (SPD)** if it is real, $A = A^\top$, and $\mathbf{x}^\top A\mathbf{x} > 0$ for every $\mathbf{x} \ne \mathbf{0}$; **Hermitian positive definite** is the complex version. $L$ is lower triangular, $U$ and $R$ upper triangular; a *unit* triangular matrix has ones on the diagonal. The **Cholesky factor** is the lower triangular $L$ with *positive* diagonal in $A = LL^\top$ (SciPy returns the upper factor $L^\top$ unless you pass `lower=True`). Singular values are ordered $\sigma_1 \ge \sigma_2 \ge \cdots \ge 0$. In §3.7, where $U$ is the polar factor, the left singular vectors are written $W$ so $A = W\Sigma V^\top$; the polar decomposition is $A = UP$ here, while LA-9 §3.11 wrote the same factorization $A = QS$. Point sets for Procrustes are the *columns* of $X$ (source) and $Y$ (target). A **flop** is one floating-point add or multiply (Trefethen & Bau's convention); costs quote only the leading term in $n$. $A_k$ is the leading $k\times k$ principal submatrix of $A$.

## 0. Where this sits

By the end of LA-9 you have met six ways of writing a matrix as a product: $PA = LU$ (LA-2), $A = QR$ (LA-5), $A = X\Lambda X^{-1}$ (LA-7), $A = Q\Lambda Q^\top$ (LA-8), $A = U\Sigma V^\top$ (LA-9), and the pseudoinverse built from the last one. You met them one at a time, each as the climax of its own module. This module puts them side by side, adds the two you have not seen ($LDL^\top$ with Cholesky, and Schur), returns to the polar decomposition of LA-9 §3.11 for the problem it was made for (Procrustes), and answers the question none of the earlier modules could: *given a task, which one do you reach for, and why?*

It is a short module, organized as an atlas rather than a story. The new proofs are Cholesky (from $LDL^\top$, via the Schur complement) and Schur (by induction); Schur then gives the general spectral theorem for normal matrices in two lines. LA-13 takes the "how is it computed" column of the atlas and does it properly: floating point, conditioning, backward stability, Householder, the QR algorithm, Krylov methods. LA-11 handles the one factorization you should never compute numerically, the Jordan form. ST-3 (regression), CM-6 (normal modes), ML-6 (PCA) and QM-4 (observables) each live on one row of the table in §3.9.

## 1. The problem that forced this idea

Every numerical linear-algebra task is the same task: **factor $A$ into pieces that are easy, then do the easy thing to each piece.** Solving $A\mathbf{x} = \mathbf{b}$ is hard; solving $L\mathbf{y} = \mathbf{b}$ by forward substitution and then $U\mathbf{x} = \mathbf{y}$ by back substitution is easy. Computing $A^{1000}$ is hard; computing $X\Lambda^{1000}X^{-1}$ is one line. Finding the closest rank-$k$ matrix to $A$ looks impossible until you write $A = U\Sigma V^\top$ and keep $k$ terms (Eckart–Young, LA-9).

This "decompositional approach" is not how the subject started. Gauss (1809) eliminated unknowns from his least-squares normal equations one at a time; he was solving systems, not factoring matrices. Matrix forms of elimination appeared in the 1930s and 40s (Banachiewicz, Dwyer; LA-2), but the paper that made "elimination *is* a factorization $A = LU$" the standard view, and asked what rounding error does to it, was Alan Turing's, written in 1948 at the National Physical Laboratory, where he and James Wilkinson had been solving 18-unknown systems on desk calculators. That paper also coined the term *condition number*. The view spread through the 1950s and 60s. In 2000, when *Computing in Science & Engineering* named the ten algorithms with the greatest influence on twentieth-century science, one entry was not an algorithm at all but "the decompositional approach to matrix computation," which G. W. Stewart wrote up in an essay of that name. The point of the view is reuse: factor once, then solve for a thousand right-hand sides, or compute a determinant, or test positive definiteness, all from the same $L$ and $U$.

Why more than one factorization, then? Because "easy" depends on what you want to do next, and on what you already know about $A$.

- André-Louis Cholesky, a French artillery officer doing geodesic triangulation around 1910, needed to solve normal equations $A^\top A\mathbf{x} = A^\top\mathbf{b}$ by hand. The matrix is symmetric and positive definite. He saw that half of the elimination work is redundant and that the result can be written $A = LL^\top$, a matrix square root. He died in 1918; a colleague published the method in 1924. When a human is the computer, a factor of two in work is the difference between finishing the survey and not finishing it.
- Issai Schur proved in 1909 that every square complex matrix is unitarily similar to a triangular one. It looked like a lemma. Half a century later John Francis (1961; Vera Kublanovskaya independently, the same year) found an iteration, the QR algorithm, that converges to exactly this triangular form, and it became the way every eigenvalue in every library is computed. The Jordan form, which textbooks spend weeks on, is never computed numerically: an arbitrarily small perturbation changes it discontinuously, while the Schur form changes continuously.
- The SVD was found twice in the 1870s (Beltrami 1873, Jordan 1874) as a statement about bilinear forms and stayed a theorem rather than a tool for ninety years (Schmidt, Weyl and Eckart–Young used it on paper; LA-9 §1) until Golub and Kahan (1965) gave a stable algorithm. Only then did it become the everyday tool for rank, conditioning, least squares and data.
- The polar decomposition $A = UP$, a rotation times a stretch, is the matrix version of $z = e^{i\theta}r$. It answers a question geodesists, then crystallographers and roboticists, kept asking: *what rigid rotation best aligns these measured points with those reference points?* Schönemann (1966) solved this "orthogonal Procrustes problem" with the SVD; Kabsch (1976) made it standard in crystallography.

So the atlas is not a list of eight unrelated facts. It is eight answers to "what does *simple* mean for this task?" Students who learn the factorizations as isolated facts can compute each one and cannot choose between them. This module is about choosing.

## 2. Intuition

### Every factorization is a change of basis that makes $A$ simple

Write a generic factorization as
$$A = (\text{basis change}) \times (\text{simple matrix}) \times (\text{basis change}).$$
"Simple" always means **diagonal** or **triangular**. The real variety is in what kind of basis change you allow, and that is dictated by what you want to preserve:

| What you want to preserve | Allowed basis change | Name | What it can make $A$ |
|---|---|---|---|
| Solutions of $A\mathbf{x} = \mathbf{b}$ | Row operations on the left, $A \mapsto EA$ | (one-sided) | Triangular: $U = L^{-1}A$ |
| The quadratic form $\mathbf{x}^\top A\mathbf{x}$ (an energy) | $A \mapsto X^\top AX$, $X$ invertible | congruence | Diagonal $D$; only the signs of $D$ are invariant (Sylvester's law of inertia, LA-12) |
| Eigenvalues, powers, $e^{tA}$, dynamics | $A \mapsto X^{-1}AX$ | similarity | Diagonal $\Lambda$ when possible; triangular always |
| The above *and* lengths and angles | $A \mapsto Q^*AQ$, $Q$ unitary | unitary similarity | Triangular always (Schur); diagonal iff $A$ is normal |
| Lengths, with different bases for input and output | $A \mapsto U^*AV$, $U, V$ unitary | unitary equivalence | Diagonal always (SVD), even for rectangular $A$ |

Read down the table and the atlas falls out.

- **LU, $LDL^\top$, Cholesky** are the first row (row operations). Cholesky is also in the second row: $L^{-1}AL^{-\top} = I$ says an SPD matrix is congruent to the identity.
- **QR** is the first row with an orthogonal $Q$: $Q^\top A = R$. You pay twice LU and get an orthonormal basis of the column space for free.
- **Eigendecomposition** $A = X\Lambda X^{-1}$ is the third row. It does not always exist.
- **Schur** $A = QTQ^*$ is the fourth row. It always exists. It is the eigendecomposition's honest cousin: the same eigenvalues (on the diagonal of $T$), an orthonormal basis, and no promise of diagonality.
- **Spectral** $A = Q\Lambda Q^*$ is Schur when $T$ comes out diagonal. Guess now which matrices those are. (Hint: a diagonal $T$ commutes with $T^*$, and conjugating by $Q$ does not change whether a matrix commutes with its conjugate transpose.)
- **SVD** $A = U\Sigma V^\top$ is the fifth row.
- **Polar** $A = UP$ is the SVD regrouped: $W\Sigma V^\top = (WV^\top)(V\Sigma V^\top)$.

### The picture

Think of a linear map on $\mathbb{R}^2$. The SVD says: every such map is *rotate, stretch along the axes, rotate*. The polar decomposition says the same thing with the two rotations merged: *stretch along some orthogonal axes (the eigenvectors of $P$), then rotate*. The spectral theorem says: for a symmetric map there is a single rotation that does both jobs, so the map is *stretch along some orthogonal axes* and nothing else. Schur says: for any map at all, after one rotation the map is *triangular*: the first axis is mapped onto itself (it is an eigenvector), the second axis into the plane of the first two, and so on. Diagonal is "every axis is an eigenvector"; triangular is "every axis stays in the span of itself and the earlier ones."

```
   SVD:       rotate by V^T  ->  stretch by Σ (axes)      ->  rotate by W
   polar:     stretch by P (its own orthogonal axes)      ->  rotate by U
   spectral:  stretch by Λ along the orthonormal axes Q   (symmetric A only)
   Schur:     in the orthonormal basis Q, A is triangular (any square A)
```

### The guesses you should be able to make

1. Cholesky exists exactly when all the elimination pivots are positive, and that happens exactly when $A$ is SPD. (The pivots are the diagonal of $D$ in $A = LDL^\top$, and $\mathbf{x}^\top A\mathbf{x} = \mathbf{y}^\top D\mathbf{y}$ with $\mathbf{y} = L^\top\mathbf{x}$.)
2. Schur exists for every complex square matrix, because every complex square matrix has at least one eigenvector, and you can peel it off and repeat.
3. The matrices whose Schur form is diagonal are exactly those with $AA^* = A^*A$.
4. If $A$ is invertible its polar factors are forced: $P$ has to be $(A^*A)^{1/2}$, and then $U$ has to be $AP^{-1}$.
5. For "which factorization?", the deciding facts about $A$ are, in order: square? symmetric? positive definite? full rank? And about the task: solve once, solve many times, least squares, eigenvalues, or geometry?

## 3. Formal development

### 3.1 Three kinds of "the same matrix"

**Definition 3.1.** Let $A, B$ be $n\times n$. $B$ is **similar** to $A$ if $B = X^{-1}AX$ for some invertible $X$; **unitarily similar** if $X$ is unitary; **congruent** to $A$ if $B = X^*AX$ for some invertible $X$. For $m\times n$ matrices, $B$ is **unitarily equivalent** to $A$ if $B = U^*AV$ with $U, V$ unitary.

Similarity preserves eigenvalues: $\det(X^{-1}AX - \lambda I) = \det(X^{-1}(A - \lambda I)X) = \det(A - \lambda I)$ (LA-6, LA-7). Congruence preserves the quadratic form up to a change of variable, $\mathbf{y}^*B\mathbf{y} = (X\mathbf{y})^*A(X\mathbf{y})$, and so preserves positive definiteness and, more precisely, the numbers of positive, negative and zero eigenvalues. That last statement is **Sylvester's law of inertia**; we use it in this module (Example 4.1, Exercise B7) and LA-12 Theorem 3.12 proves it. Unitary equivalence preserves singular values, because $B^*B = V^*A^*AV$ is similar to $A^*A$. Unitary similarity is all three at once, which is why Schur is the factorization that respects everything.

### 3.2 Triangular factorizations: LU, $LDL^\top$, Cholesky

**Theorem 3.2 (LU with partial pivoting; LA-2).** Every $n\times n$ matrix $A$ can be written $PA = LU$ with $P$ a permutation, $L$ unit lower triangular with $|l_{ij}| \le 1$, and $U$ upper triangular. If $A$ is invertible, $U$ has nonzero diagonal. Cost: $\tfrac23 n^3$ flops. Each later solve of $A\mathbf{x} = \mathbf{b}$ costs $2n^2$ (two triangular substitutions).

The cost count is the one from LA-2 that you should redo now: step $k$ of elimination updates an $(n-k)\times(n-k)$ block at 2 flops per entry, and $\sum_{k=1}^{n}2(n-k)^2 \approx \tfrac23 n^3$.

Now specialize to symmetric $A$. If elimination *without* row swaps goes through (every pivot nonzero), then $A = LU$ with $L$ unit lower triangular, and symmetry forces more.

▶ **Try it first.** Suppose $A = A^\top = LU$ with $L$ unit lower triangular and $U$ upper triangular with nonzero diagonal. Let $D = \operatorname{diag}(u_{11}, \dots, u_{nn})$. Show that $U = DL^\top$.

**Theorem 3.3 ($LDL^\top$).** If $A$ is symmetric and elimination without pivoting produces nonzero pivots $d_1, \dots, d_n$, then $A = LDL^\top$ with $L$ unit lower triangular and $D = \operatorname{diag}(d_1, \dots, d_n)$. Cost: $\tfrac13 n^3$ flops.

**Derivation.**

1. Write $U = DM$ where $M = D^{-1}U$ is unit upper triangular (divide each row of $U$ by its pivot). So $A = LDM$.
2. Transpose: $A = A^\top = M^\top DL^\top$. Here $M^\top$ is unit lower triangular and $DL^\top$ is upper triangular with the same nonzero diagonal, so $A = M^\top(DL^\top)$ is *also* an LU factorization of $A$ with unit lower factor.
3. The LU factorization with unit lower factor is unique when it exists (LA-2: if $L_1U_1 = L_2U_2$ with $U_1$ invertible, then $L_2^{-1}L_1 = U_2U_1^{-1}$ is both unit lower triangular and upper triangular, hence $I$). So $M^\top = L$, i.e. $M = L^\top$, and $A = LDL^\top$.
4. Cost: the block that elimination updates at step $k$ is symmetric, so only its lower triangle, half of $(n-k)^2$ entries, needs computing. Half of $\tfrac23 n^3$ is $\tfrac13 n^3$. ∎

■ **What just happened.** Symmetry of $A$ forces symmetry of the factorization: the upper factor is the transpose of the lower one, up to the diagonal of pivots. Half the arithmetic was redundant.

The pivots $d_k$ are what connect elimination to positive definiteness.

▶ **Try it first.** With $A = LDL^\top$, substitute $\mathbf{y} = L^\top\mathbf{x}$ into $\mathbf{x}^\top A\mathbf{x}$. What do you get, and why does $\mathbf{y}$ range over all of $\mathbb{R}^n$ as $\mathbf{x}$ does?

**Theorem 3.4 (Positive pivots $\iff$ SPD).** Let $A = LDL^\top$ as in Theorem 3.3. Then $A$ is SPD if and only if every $d_k > 0$.

**Derivation.**

1. For any $\mathbf{x}$, with $\mathbf{y} = L^\top\mathbf{x}$,
$$\mathbf{x}^\top A\mathbf{x} = \mathbf{x}^\top LDL^\top\mathbf{x} = (L^\top\mathbf{x})^\top D(L^\top\mathbf{x}) = \mathbf{y}^\top D\mathbf{y} = \sum_{k=1}^n d_ky_k^2. \tag{3.1}$$
2. $L^\top$ is invertible (unit triangular, determinant 1), so $\mathbf{x} \mapsto \mathbf{y}$ is a bijection of $\mathbb{R}^n$, and $\mathbf{x} \ne \mathbf{0} \iff \mathbf{y} \ne \mathbf{0}$.
3. ($\Leftarrow$) If all $d_k > 0$, then for $\mathbf{y} \ne \mathbf{0}$ the sum in (3.1) is positive, so $A$ is SPD.
4. ($\Rightarrow$) If $A$ is SPD, take $\mathbf{y} = \mathbf{e}_k$, i.e. $\mathbf{x} = L^{-\top}\mathbf{e}_k \ne \mathbf{0}$. Then (3.1) reads $0 < \mathbf{x}^\top A\mathbf{x} = d_k$. ∎

■ **What just happened.** Congruence by $L^{-\top}$ turns $A$ into the diagonal matrix $D$, and a diagonal quadratic form is positive exactly when its coefficients are. This is also the fastest practical test for positive definiteness: run elimination (or Cholesky, below) and watch the pivots. No eigenvalues required.

There is a gap. Theorem 3.3 *assumed* that elimination without pivoting goes through. For SPD matrices it always does, and the cleanest way to see it is to do one step of elimination by hand and look at what is left.

▶ **Try it first.** Partition an SPD matrix as $A = \begin{pmatrix} a & \mathbf{b}^\top \\ \mathbf{b} & C\end{pmatrix}$ with $a$ a scalar. Show $a > 0$. Then multiply out $\begin{pmatrix} 1 & \mathbf{0}^\top \\ \mathbf{b}/a & I\end{pmatrix}\begin{pmatrix} a & \mathbf{0}^\top \\ \mathbf{0} & S\end{pmatrix}\begin{pmatrix} 1 & \mathbf{b}^\top/a \\ \mathbf{0} & I\end{pmatrix}$ and find the $S$ that makes it equal to $A$.

**Theorem 3.5 (Cholesky factorization).** Every SPD matrix $A$ can be written $A = LL^\top$ with $L$ lower triangular with positive diagonal. Elimination without pivoting never meets a zero or negative pivot. $L$ is unique (Exercise B1). Cost: $\tfrac13 n^3$ flops.

**Derivation** (induction on $n$).

1. $n = 1$: $A = (a)$ with $a > 0$; $L = (\sqrt a)$.
2. For $n > 1$, partition $A = \begin{pmatrix} a & \mathbf{b}^\top \\ \mathbf{b} & C\end{pmatrix}$ with $a \in \mathbb{R}$, $\mathbf{b} \in \mathbb{R}^{n-1}$, $C$ symmetric $(n-1)\times(n-1)$. Then $a = \mathbf{e}_1^\top A\mathbf{e}_1 > 0$: the first pivot is positive.
3. Multiply out the product in the prompt:
$$\begin{pmatrix} 1 & \mathbf{0}^\top \\ \mathbf{b}/a & I\end{pmatrix}\begin{pmatrix} a & \mathbf{0}^\top \\ \mathbf{0} & S\end{pmatrix}\begin{pmatrix} 1 & \mathbf{b}^\top/a \\ \mathbf{0} & I\end{pmatrix} = \begin{pmatrix} a & \mathbf{b}^\top \\ \mathbf{b} & \mathbf{b}\mathbf{b}^\top/a + S\end{pmatrix}.$$
This equals $A$ exactly when
$$S = C - \frac{\mathbf{b}\mathbf{b}^\top}{a}, \tag{3.2}$$
the **Schur complement** of $a$ in $A$: it is what elimination leaves in the bottom-right block after step 1.
4. $S$ is SPD. It is symmetric because $\mathbf{b}\mathbf{b}^\top$ is. For positivity, take $\mathbf{y} \in \mathbb{R}^{n-1}$, $\mathbf{y} \ne \mathbf{0}$, and set $\mathbf{x} = \begin{pmatrix} -\mathbf{b}^\top\mathbf{y}/a \\ \mathbf{y}\end{pmatrix} \ne \mathbf{0}$. Expanding $\mathbf{x}^\top A\mathbf{x}$ with $x_1 = -\mathbf{b}^\top\mathbf{y}/a$:
$$\mathbf{x}^\top A\mathbf{x} = ax_1^2 + 2x_1\mathbf{b}^\top\mathbf{y} + \mathbf{y}^\top C\mathbf{y} = \frac{(\mathbf{b}^\top\mathbf{y})^2}{a} - \frac{2(\mathbf{b}^\top\mathbf{y})^2}{a} + \mathbf{y}^\top C\mathbf{y} = \mathbf{y}^\top\Big(C - \frac{\mathbf{b}\mathbf{b}^\top}{a}\Big)\mathbf{y} = \mathbf{y}^\top S\mathbf{y}.$$
The left side is positive because $A$ is SPD. So $\mathbf{y}^\top S\mathbf{y} > 0$.
5. By the induction hypothesis $S = L_SL_S^\top$ with $L_S$ lower triangular with positive diagonal. Since $\operatorname{diag}(\sqrt a, L_S)\operatorname{diag}(\sqrt a, L_S)^\top = \operatorname{diag}(a, S)$, step 3 gives
$$A = \begin{pmatrix} 1 & \mathbf{0}^\top \\ \mathbf{b}/a & I\end{pmatrix}\begin{pmatrix} \sqrt a & \mathbf{0}^\top \\ \mathbf{0} & L_S\end{pmatrix}\begin{pmatrix} \sqrt a & \mathbf{0}^\top \\ \mathbf{0} & L_S\end{pmatrix}^{\!\top}\begin{pmatrix} 1 & \mathbf{b}^\top/a \\ \mathbf{0} & I\end{pmatrix} = LL^\top, \qquad L = \begin{pmatrix} \sqrt a & \mathbf{0}^\top \\ \mathbf{b}/\sqrt a & L_S\end{pmatrix},$$
and $L$ is lower triangular with positive diagonal.
6. Pivots: step 2 showed the first pivot is $a > 0$; the remaining pivots are the pivots of the SPD matrix $S$, positive by induction. Elimination never needs a row swap. ∎

■ **What just happened.** One step of elimination on an SPD matrix leaves a smaller SPD matrix. That is the entire reason Cholesky needs no pivoting, and it is why Cholesky is the definitive test for positive definiteness: if a pivot comes out $\le 0$, $A$ was not SPD (Theorem 3.4), and the step number tells you the size of the smallest leading block that already fails. Comparing with Theorem 3.3, $L_{\text{Chol}} = L_{\text{unit}}D^{1/2}$: $A = LDL^\top = (LD^{1/2})(LD^{1/2})^\top$ is the same factorization in two notations. One warning: an SPD matrix that is nearly singular can have a pivot that is positive in exact arithmetic but comes out slightly negative in floating point. The theorem is fine; your matrix is numerically singular. LA-13 quantifies this.

**The direct algorithm, by entry matching.** The induction is a proof; the algorithm is faster to state by equating entries of $A = LL^\top$. Since $L$ is lower triangular, $(LL^\top)_{ij} = \sum_{k=1}^{\min(i,j)}l_{ik}l_{jk}$. Take $j \le i$ and peel off the last term: $a_{ij} = \sum_{k<j}l_{ik}l_{jk} + l_{ij}l_{jj}$. Solve column by column, $j = 1, \dots, n$:
$$l_{jj} = \sqrt{\,a_{jj} - \sum_{k<j}l_{jk}^2\,}, \qquad l_{ij} = \frac{a_{ij} - \sum_{k<j}l_{ik}l_{jk}}{l_{jj}}\quad(i > j). \tag{3.3}$$
When we reach column $j$, every $l$ on the right side is already known (columns $< j$, and $l_{jj}$ before the $l_{ij}$). The quantity under the square root is the $j$-th pivot $d_j$ (because $l_{jj} = \sqrt{d_j}$), so (3.3) fails with a non-positive radicand exactly when $A$ is not SPD.

**Cost of (3.3).** Column $j$: $l_{jj}$ needs $j - 1$ multiplications and $j - 1$ subtractions; each of the $n - j$ entries $l_{ij}$ needs $2(j-1) + 1$ flops. Total
$$\sum_{j=1}^n 2(j-1)(n-j) \approx 2\int_0^n t(n-t)\,dt = 2\Big(\frac{n^3}{2} - \frac{n^3}{3}\Big) = \frac{n^3}{3}.$$
Half of LU, and no permutation bookkeeping, which matters more than the flop count for banded and sparse matrices: the pattern of nonzeros in $L$ is known before you start.

### 3.3 Orthogonal–triangular: QR

**Theorem 3.6 (QR; LA-5).** Every $m\times n$ matrix $A$ with $m \ge n$ can be written $A = QR$ with $Q$ ($m\times n$) having orthonormal columns and $R$ ($n\times n$) upper triangular (the *reduced* QR), or with $Q$ ($m\times m$) orthogonal and $R$ ($m\times n$) (the *full* QR). If $A$ has full column rank, the reduced factorization with positive diagonal in $R$ is unique. Cost by Householder reflections: $2mn^2 - \tfrac23 n^3$ flops ($\tfrac43 n^3$ for square $A$; forming $Q$ explicitly costs about as much again).

What it is for: least squares. $\min_{\mathbf{x}}\|A\mathbf{x} - \mathbf{b}\|$ is solved by $R\mathbf{x} = Q^\top\mathbf{b}$ (LA-5) without ever forming $A^\top A$. The reason not to form $A^\top A$ is one line:

**Theorem 3.7 (QR and Cholesky; conditioning).** If $A = QR$ is a reduced QR with $R$ having positive diagonal, then $A^\top A = R^\top R$ is the Cholesky factorization of $A^\top A$, and $\kappa_2(A^\top A) = \kappa_2(A)^2$.

*Proof.* $A^\top A = R^\top Q^\top QR = R^\top R$, and $R^\top$ is lower triangular with positive diagonal, so by uniqueness (Exercise B1) it is the Cholesky factor. For the condition number, write $A = U\Sigma V^\top$ (LA-9): $A^\top A = V\Sigma^\top\Sigma V^\top$ has singular values $\sigma_i^2$, so $\kappa_2(A^\top A) = \sigma_1^2/\sigma_n^2 = \kappa_2(A)^2$. ∎

A condition number of $10^8$ costs about 8 of your 16 digits in a solve (LA-9; LA-13 makes this precise); squaring it costs all 16. The normal equations are the one route to least squares you should refuse by default. How QR is computed: Householder reflections (LA-13), not Gram–Schmidt, which loses orthogonality in floating point (also LA-13).

### 3.4 Similarity to diagonal: eigendecomposition and the spectral theorem

**Eigendecomposition (LA-7).** If the $n\times n$ matrix $A$ has $n$ linearly independent eigenvectors, then $A = X\Lambda X^{-1}$ with $\Lambda = \operatorname{diag}(\lambda_1, \dots, \lambda_n)$ and the eigenvectors as the columns of $X$. It exists iff $A$ is diagonalizable (guaranteed when the eigenvalues are distinct; it fails for $\begin{pmatrix} 0 & 1 \\ 0 & 0\end{pmatrix}$). For: powers $A^k = X\Lambda^kX^{-1}$, the matrix exponential $e^{tA} = Xe^{t\Lambda}X^{-1}$ for ODE systems, Markov chains. Not for: numerical work when $X$ is ill-conditioned, because every error is amplified by $\kappa(X)$. Computed via Schur (§3.5), never by finding roots of the characteristic polynomial.

**Spectral theorem (LA-8).** If $A^* = A$ (Hermitian; real symmetric as a special case), then $A = Q\Lambda Q^*$ with $Q$ unitary and $\Lambda$ real diagonal: an orthonormal basis of eigenvectors with real eigenvalues. Here $X = Q$ has $\kappa(Q) = 1$, so the basis change is perfectly conditioned. For: quadratic forms, principal axes (the inertia tensor in CM-7), normal modes (CM-6), covariance matrices and PCA (LA-9, ML-6), observables (QM-4). Computed by Householder tridiagonalization followed by the implicit QR iteration or divide-and-conquer: about $\tfrac43 n^3$ flops for the eigenvalues, about $9n^3$ including eigenvectors. Jacobi rotations are an alternative, slower but more accurate for the small eigenvalues.

### 3.5 Schur: the factorization that always exists

**Theorem 3.8 (Schur, 1909).** Every $A \in \mathbb{C}^{n\times n}$ can be written
$$A = QTQ^*, \qquad Q\ \text{unitary},\quad T\ \text{upper triangular}. \tag{3.4}$$
The diagonal of $T$ consists of the eigenvalues of $A$, and they can be made to appear in any prescribed order.

▶ **Try it first.** Take a unit eigenvector $\mathbf{q}_1$ of $A$ and extend it to a unitary matrix $Q_1 = \begin{pmatrix} \mathbf{q}_1 & \hat Q_1\end{pmatrix}$. Compute $Q_1^*AQ_1$ in block form. What is its first column?

**Derivation** (induction on $n$).

1. $n = 1$: $A = (a)$, $Q = (1)$, $T = (a)$.
2. Let $n > 1$. The characteristic polynomial $\det(A - \lambda I)$ has degree $n \ge 1$ and so has a root $\lambda_1 \in \mathbb{C}$ (fundamental theorem of algebra, F1; this is the one place "complex" is used). Then $A - \lambda_1I$ is singular, so there is a unit vector $\mathbf{q}_1$ with $A\mathbf{q}_1 = \lambda_1\mathbf{q}_1$.
3. Extend $\mathbf{q}_1$ to an orthonormal basis $\mathbf{q}_1, \dots, \mathbf{q}_n$ of $\mathbb{C}^n$ (Gram–Schmidt, LA-5, applied to $\mathbf{q}_1$ followed by the standard basis, discarding dependent vectors). Let $Q_1 = \begin{pmatrix} \mathbf{q}_1 & \hat Q_1\end{pmatrix}$, where $\hat Q_1$ is $n\times(n-1)$ with $\hat Q_1^*\mathbf{q}_1 = \mathbf{0}$ and $\hat Q_1^*\hat Q_1 = I$.
4. The first column of $Q_1^*AQ_1$ is $Q_1^*A\mathbf{q}_1 = \lambda_1Q_1^*\mathbf{q}_1 = \lambda_1\mathbf{e}_1$, because the columns of $Q_1$ are orthonormal. So
$$Q_1^*AQ_1 = \begin{pmatrix} \lambda_1 & \mathbf{w}^* \\ \mathbf{0} & A_2\end{pmatrix}, \qquad \mathbf{w}^* = \mathbf{q}_1^*A\hat Q_1, \quad A_2 = \hat Q_1^*A\hat Q_1 \in \mathbb{C}^{(n-1)\times(n-1)}. \tag{3.5}$$
This is the **deflation** step: one eigenvalue has been split off and the rest of the problem is the smaller matrix $A_2$.
5. By the induction hypothesis, $A_2 = Q_2T_2Q_2^*$ with $Q_2$ unitary and $T_2$ upper triangular, both $(n-1)\times(n-1)$.
6. Set $Q = Q_1\begin{pmatrix} 1 & \mathbf{0}^* \\ \mathbf{0} & Q_2\end{pmatrix}$, a product of unitary matrices, hence unitary. Then
$$Q^*AQ = \begin{pmatrix} 1 & \mathbf{0}^* \\ \mathbf{0} & Q_2^*\end{pmatrix}\begin{pmatrix} \lambda_1 & \mathbf{w}^* \\ \mathbf{0} & A_2\end{pmatrix}\begin{pmatrix} 1 & \mathbf{0}^* \\ \mathbf{0} & Q_2\end{pmatrix} = \begin{pmatrix} \lambda_1 & \mathbf{w}^*Q_2 \\ \mathbf{0} & Q_2^*A_2Q_2\end{pmatrix} = \begin{pmatrix} \lambda_1 & \mathbf{w}^*Q_2 \\ \mathbf{0} & T_2\end{pmatrix} =: T,$$
which is upper triangular.
7. Eigenvalues: $T$ is similar to $A$, so they share eigenvalues, and the eigenvalues of a triangular matrix are its diagonal entries, since $\det(T - \lambda I) = \prod_i(t_{ii} - \lambda)$. Order: in step 2 we may choose *any* eigenvalue as $\lambda_1$, and recursively any order for the rest. ∎

■ **What just happened.** Three facts, used once each: a complex polynomial has a root (one eigenvector), Gram–Schmidt (complete it to a unitary basis), and $A\mathbf{q}_1 = \lambda_1\mathbf{q}_1$ (zeros under the diagonal in the first column). Everything else is block bookkeeping. Compare with the eigendecomposition: we never needed $n$ independent eigenvectors, only one at a time. The price is that $T$ keeps a strictly upper triangular part $\mathbf{w}^*Q_2$ recording how far $A$ is from having an orthonormal eigenbasis.

Two consequences you will use immediately.

**Corollary 3.8a (Unitary invariance of the Frobenius norm; departure from normality).** $\|A\|_F^2 = \|T\|_F^2 = \sum_i|\lambda_i|^2 + \sum_{i<j}|t_{ij}|^2$. So the strictly upper part of $T$ has a size $\sqrt{\|A\|_F^2 - \sum_i|\lambda_i|^2}$ that does not depend on which Schur form you compute, even though the individual entries $t_{ij}$ do (Example 4.2).

*Proof.* $\|QTQ^*\|_F^2 = \operatorname{tr}\big((QTQ^*)^*(QTQ^*)\big) = \operatorname{tr}(QT^*TQ^*) = \operatorname{tr}(T^*TQ^*Q) = \operatorname{tr}(T^*T) = \|T\|_F^2$, using $\operatorname{tr}(XY) = \operatorname{tr}(YX)$. Then split $T$ into its diagonal and its strictly upper part; they have no entries in common. ∎

**Corollary 3.8b (Eigenvectors from $T$).** If $\mathbf{y} \ne \mathbf{0}$ solves $(T - \lambda_kI)\mathbf{y} = \mathbf{0}$, then $\mathbf{x} = Q\mathbf{y}$ is an eigenvector of $A$ for $\lambda_k$. When the eigenvalues are distinct, $\mathbf{y}$ is found by back substitution: set $y_k = 1$ and $y_j = 0$ for $j > k$, then solve rows $k-1, \dots, 1$ upward. This is how `eig` gets eigenvectors: Schur first, then one triangular solve per eigenvalue.

*Proof.* $A\mathbf{x} = QTQ^*Q\mathbf{y} = QT\mathbf{y} = Q(\lambda_k\mathbf{y}) = \lambda_k\mathbf{x}$. ∎

**Normal matrices.** Now the question from §2: when does $T$ come out diagonal?

**Definition.** $A$ is **normal** if $AA^* = A^*A$. Hermitian ($A^* = A$), skew-Hermitian ($A^* = -A$) and unitary ($A^*A = I$) matrices are all normal. $\begin{pmatrix} 1 & 2 \\ 3 & 0\end{pmatrix}$ is not (Example 4.2).

▶ **Try it first.** Let $T$ be upper triangular with $TT^* = T^*T$. Write the $(1,1)$ entry of each side as a sum over the first *row* of $T$, respectively the first *column* of $T$. What does equality force?

**Lemma 3.9 (Normal and triangular implies diagonal).** If $T$ is upper triangular and normal, then $T$ is diagonal.

**Derivation.** We show, for $k = 1, 2, \dots, n$ in turn, that row $k$ of $T$ is zero off the diagonal.

1. The $(k,k)$ entry of $TT^*$ is $\sum_jt_{kj}\overline{t_{kj}} = \sum_{j\ge k}|t_{kj}|^2$: the squared length of **row** $k$ of $T$, which starts at the diagonal because $T$ is upper triangular.
2. The $(k,k)$ entry of $T^*T$ is $\sum_i\overline{t_{ik}}t_{ik} = \sum_{i\le k}|t_{ik}|^2$: the squared length of **column** $k$ of $T$, which ends at the diagonal.
3. $k = 1$: step 1 gives $|t_{11}|^2 + \sum_{j>1}|t_{1j}|^2$; step 2 gives $|t_{11}|^2$ alone, because column 1 of an upper triangular matrix has only its diagonal entry. Normality says they are equal, so $\sum_{j>1}|t_{1j}|^2 = 0$, so $t_{1j} = 0$ for all $j > 1$.
4. General $k$: suppose rows $1, \dots, k-1$ are already known to be diagonal. The entries above the diagonal in column $k$ are $t_{ik}$ for $i < k$; they sit in those rows, so they are zero, and step 2 collapses to $|t_{kk}|^2$. Equating with step 1: $\sum_{j>k}|t_{kj}|^2 = 0$, so row $k$ is zero off the diagonal.
5. After $k = n$, every off-diagonal entry is zero. ∎

■ **What just happened.** Row lengths and column lengths of $T$ agree when $T$ is normal; triangularity makes column 1 contain only its diagonal entry, which kills row 1's tail; then column 2 contains only its diagonal entry, and so on. The step people get wrong is step 4: you must use the zeros you already produced to *shorten the column sum* before comparing. Without that, comparing the $(k,k)$ entries gives $\sum_{j>k}|t_{kj}|^2 = \sum_{i<k}|t_{ik}|^2$, which is true and proves nothing on its own.

**Theorem 3.10 (Spectral theorem for normal matrices).** $A \in \mathbb{C}^{n\times n}$ is unitarily diagonalizable, $A = Q\Lambda Q^*$ with $Q$ unitary and $\Lambda$ diagonal, if and only if $A$ is normal.

*Proof.* ($\Rightarrow$) $AA^* = Q\Lambda Q^*Q\bar\Lambda Q^* = Q\Lambda\bar\Lambda Q^* = Q\bar\Lambda\Lambda Q^* = A^*A$, since diagonal matrices commute. ($\Leftarrow$) Take a Schur form $A = QTQ^*$. Then $TT^* = Q^*AQ\,Q^*A^*Q = Q^*AA^*Q = Q^*A^*AQ = T^*T$: unitary similarity preserves normality. By Lemma 3.9, $T$ is diagonal. ∎

**Corollaries.** With $T = \Lambda$ diagonal: Hermitian $A$ forces $\bar\lambda_i = \lambda_i$ (compare the diagonals of $T$ and $T^*$), so the eigenvalues are real, and LA-8's spectral theorem is recovered. Unitary $A$ forces $|\lambda_i| = 1$. Skew-Hermitian $A$ forces purely imaginary eigenvalues. For normal $A$, $\sigma_i = |\lambda_i|$ and $\|A\|_2 = \max_i|\lambda_i|$ (Exercise C1). For non-normal $A$ none of these hold, and the strictly upper part of $T$ measures the failure (Corollary 3.8a).

**Real Schur form.** If $A$ is real its eigenvalues come in conjugate pairs, so a real orthogonal $Q$ cannot in general make $T$ triangular: a real triangular matrix has real eigenvalues. What is true:

**Theorem 3.11 (Real Schur form; stated).** Every real $A \in \mathbb{R}^{n\times n}$ can be written $A = QTQ^\top$ with $Q$ real orthogonal and $T$ real *quasi*-upper-triangular: block upper triangular whose diagonal blocks are $1\times1$ (the real eigenvalues) or $2\times2$ with a pair of complex conjugate eigenvalues $a \pm ib$.

The proof is the same induction with one change: for a complex eigenvalue $a + ib$ with eigenvector $\mathbf{x} + i\mathbf{y}$, the real plane $\operatorname{span}\{\mathbf{x}, \mathbf{y}\}$ is invariant under $A$, and $A$ acts on it as $\begin{pmatrix} a & b \\ -b & a\end{pmatrix}$ in the basis $(\mathbf{x}, \mathbf{y})$ (Exercise B5 derives this). Deflate two dimensions at once using an orthonormal basis of that plane. The full proof is Horn & Johnson, *Matrix Analysis*, Theorem 2.3.4, or Golub & Van Loan, Theorem 7.4.1. This is what `scipy.linalg.schur` returns for real input unless you ask for `output='complex'`.

**How Schur is computed** (LA-13). Reduce $A$ to upper Hessenberg form (zeros below the first subdiagonal) by Householder reflections, $\tfrac{10}{3}n^3$ flops; then run the shifted QR iteration, which drives the subdiagonal to zero. About $10n^3$ flops for $T$ alone, $25n^3$ if $Q$ is accumulated. Every eigenvalue routine in LAPACK, hence in NumPy and SciPy, is this.

### 3.6 Two-sided orthogonal: SVD (recap)

**Theorem 3.12 (SVD; LA-9).** Every $A \in \mathbb{R}^{m\times n}$ can be written $A = U\Sigma V^\top$ with $U$ ($m\times m$) and $V$ ($n\times n$) orthogonal and $\Sigma$ ($m\times n$) diagonal with $\sigma_1 \ge \cdots \ge \sigma_r > 0 = \sigma_{r+1} = \cdots$, where $r = \operatorname{rank}A$. The $\sigma_i^2$ are the eigenvalues of $A^\top A$ (eigenvectors $V$) and of $AA^\top$ (eigenvectors $U$).

For: rank (count the $\sigma_i$ above a tolerance), null space and column space (the trailing $\mathbf{v}_i$ and the leading $\mathbf{u}_i$), condition number $\kappa_2 = \sigma_1/\sigma_n$, the pseudoinverse $A^+ = V\Sigma^+U^\top$ and rank-deficient least squares, the best rank-$k$ approximation $\sum_{i\le k}\sigma_i\mathbf{u}_i\mathbf{v}_i^\top$ (Eckart–Young), PCA. Relation to the others: for symmetric $A$, $\sigma_i = |\lambda_i|$ and $\mathbf{u}_i = \operatorname{sign}(\lambda_i)\mathbf{v}_i$, so the SVD is the spectral decomposition with the signs moved into $U$; for normal $A$ the same with $|\lambda_i|$ and phases. Computed by Golub–Kahan bidiagonalization and an implicit QR iteration on the bidiagonal (LA-13): about $4mn^2 - \tfrac43 n^3$ flops for $\Sigma$ alone, and roughly $20n^3$ for a square matrix with full $U$ and $V$. The most expensive factorization in the atlas, and the most informative.

### 3.7 Polar decomposition: rotation times stretch

**Theorem 3.13 (Polar decomposition).** Every $A \in \mathbb{C}^{n\times n}$ can be written $A = UP$ with $U$ unitary and $P$ Hermitian positive semidefinite. $P = (A^*A)^{1/2}$ is always unique. If $A$ is invertible, $U$ is unique too, and $P$ is positive definite.

▶ **Try it first.** Insert $V^*V = I$ into the middle of $A = W\Sigma V^*$ (SVD, with $W$ for the left factor so that $U$ is free) and regroup.

**Derivation.**

1. *Existence.* Let $A = W\Sigma V^*$ be an SVD. Then $A = W(V^*V)\Sigma V^* = (WV^*)(V\Sigma V^*) =: UP$. $U = WV^*$ is a product of unitaries, hence unitary. $P = V\Sigma V^*$ is Hermitian with eigenvalues $\sigma_i \ge 0$, hence positive semidefinite; it is positive definite iff all $\sigma_i > 0$ iff $A$ is invertible.
2. *$P$ is forced.* Suppose $A = UP$ with $U$ unitary and $P$ positive semidefinite. Then $A^*A = PU^*UP = P^2$. So $P$ must be a positive semidefinite square root of $A^*A$.
3. *The positive semidefinite square root is unique.* Let $B$ be positive semidefinite and let $P$ be positive semidefinite with $P^2 = B$. Diagonalize $P = VDV^*$ with $D = \operatorname{diag}(d_i)$, $d_i \ge 0$ (Theorem 3.10). Then $B = VD^2V^*$. For each eigenvalue $\mu$ of $B$, the $\mu$-eigenspace of $B$ is spanned by the columns $\mathbf{v}_i$ with $d_i^2 = \mu$, i.e. (as $d_i \ge 0$) with $d_i = \sqrt\mu$. So $P$ acts as multiplication by $\sqrt\mu$ on the $\mu$-eigenspace of $B$, for every $\mu$. Those eigenspaces span $\mathbb{C}^n$, so $B$ determines $P$ completely: $B^{1/2}$ is well defined.
4. *$U$ is forced when $A$ is invertible.* Then $A^*A$ is positive definite, so $P = (A^*A)^{1/2}$ is invertible, and $U = AP^{-1}$. ∎

■ **What just happened.** The SVD's two rotations $W$ and $V^*$ were regrouped into one rotation $U = WV^*$ and one symmetric stretch $P = V\Sigma V^*$ whose axes are the right singular vectors. The stretch is unique because the positive semidefinite square root is unique. When $A$ is singular, $U$ is *not* unique: it is free on the null space of $A$ (any choice of the last columns of $W$ works). In continuum mechanics this is $F = RU$, the deformation gradient as a rotation times the right stretch tensor (LA-9 §3.11 and Example 4.5, which wrote it $A = QS$); in attitude estimation it is "the nearest rotation to a noisy matrix." That last phrase is a theorem, and it is the orthogonal Procrustes problem.

**Theorem 3.14 (Nearest orthogonal matrix; orthogonal Procrustes).** Let $M \in \mathbb{R}^{n\times n}$ have SVD $M = W\Sigma V^\top$ with all $\sigma_i > 0$. Then
$$\arg\max_{R^\top R = I}\operatorname{tr}(R^\top M) = \arg\min_{R^\top R = I}\|M - R\|_F = WV^\top,$$
the orthogonal polar factor of $M$, and it is unique. Consequently, for point sets stored as the columns of $X, Y \in \mathbb{R}^{d\times k}$ and $M = YX^\top$, the orthogonal $R$ minimizing $\|RX - Y\|_F$ is $R = WV^\top$. (If some $\sigma_i = 0$, a minimizer still exists but need not be unique.)

▶ **Try it first.** Expand $\|RX - Y\|_F^2 = \operatorname{tr}\big((RX - Y)^\top(RX - Y)\big)$ and find the single term that depends on $R$.

**Derivation.**

1. Using $R^\top R = I$ and $\operatorname{tr}(XY) = \operatorname{tr}(YX)$,
$$\|RX - Y\|_F^2 = \operatorname{tr}(X^\top R^\top RX) - 2\operatorname{tr}(X^\top R^\top Y) + \operatorname{tr}(Y^\top Y) = \|X\|_F^2 + \|Y\|_F^2 - 2\operatorname{tr}(R^\top M), \qquad M = YX^\top. \tag{3.6}$$
Only the last term depends on $R$: minimizing the distance means **maximizing $\operatorname{tr}(R^\top M)$**. With $X = I$, (3.6) reads $\|M - R\|_F^2 = \|M\|_F^2 + n - 2\operatorname{tr}(R^\top M)$, the nearest-orthogonal-matrix problem.
2. Insert the SVD: $\operatorname{tr}(R^\top W\Sigma V^\top) = \operatorname{tr}(V^\top R^\top W\Sigma) = \operatorname{tr}(Z\Sigma)$ with $Z = V^\top R^\top W$, orthogonal as a product of orthogonal matrices. As $R$ ranges over all orthogonal matrices, so does $Z$.
3. $\operatorname{tr}(Z\Sigma) = \sum_iz_{ii}\sigma_i \le \sum_i|z_{ii}|\sigma_i \le \sum_i\sigma_i$, because every entry of an orthogonal matrix satisfies $|z_{ii}| \le 1$ (each column is a unit vector).
4. Equality requires $z_{ii} = 1$ for every $i$ (all $\sigma_i > 0$), and a unit column with a $1$ on the diagonal has zeros elsewhere, so $Z = I$. Then $V^\top R^\top W = I$, i.e. $R^\top = VW^\top$, $R = WV^\top$. ∎

■ **What just happened.** The distance from a point set to a rotated copy depends on $R$ only through $\operatorname{tr}(R^\top M)$, and the trace is maximized by lining the rotation up with the singular vectors of $M$. The minimum residual is $\|X\|_F^2 + \|Y\|_F^2 - 2\sum_i\sigma_i$, a number you can check before you compute $R$ (Example 4.3). Two practical refinements. Translation: center both point sets on their centroids first (Exercise C2 shows this is optimal). Reflections: if you need a *proper* rotation ($\det R = +1$) and $\det(WV^\top) = -1$, the best proper rotation is $W\operatorname{diag}(1, \dots, 1, -1)V^\top$ (Kabsch 1976, 1978; Umeyama 1991; we state this without proof). The project implements both.

### 3.8 Jordan form: in the atlas, not in the toolbox

$A = XJX^{-1}$ with $J$ block diagonal, each block $\lambda I + N$ where $N$ has ones on the superdiagonal and zeros elsewhere. It always exists over $\mathbb{C}$; LA-11 proves it. It is the complete answer to "what are all matrices up to similarity" and the right tool for proving things about $e^{tA}$ when $A$ is defective. It is never computed in floating point: the Jordan structure is a discontinuous function of the entries ($\begin{pmatrix} 0 & 1 \\ 0 & 0\end{pmatrix}$ has one $2\times2$ block; $\begin{pmatrix} 0 & 1 \\ \varepsilon & 0\end{pmatrix}$ has two $1\times1$ blocks for every $\varepsilon \ne 0$), so rounding changes the answer completely, and $\kappa(X)$ is unbounded. Schur is the numerical replacement: continuous, orthonormal, and sufficient for eigenvalues, invariant subspaces and matrix functions. Exercise C4 shows how Schur does a job you might have thought needed Jordan.

### 3.9 The decision table

Read the facts about $A$ off its shape and symmetry, test definiteness by attempting Cholesky, and then:

| Task | Facts about $A$ | Use | Why | Cost (flops) |
|---|---|---|---|---|
| Solve $A\mathbf{x} = \mathbf{b}$ once | square, nonsingular, no structure | PLU | cheapest general factorization; partial pivoting makes it stable in practice (LA-13) | $\tfrac23 n^3$ |
| Solve with many right-hand sides, same $A$ | as above | PLU once, then triangular solves | factor once, $2n^2$ per $\mathbf{b}$ | $\tfrac23 n^3 + 2n^2$ per $\mathbf{b}$ |
| Solve, $A$ SPD (stiffness, covariance, normal equations, Gram matrices) | SPD | Cholesky | half the flops; no pivoting, so band and sparsity structure survive; failure is a definiteness test | $\tfrac13 n^3$ |
| Solve, $A$ symmetric indefinite | symmetric, not SPD | $LDL^\top$ with symmetric pivoting (Bunch–Kaufman) | half the flops of LU; inertia for free | $\tfrac13 n^3$ |
| Least squares $\min\|A\mathbf{x} - \mathbf{b}\|$ | $m\times n$, $m \ge n$, full column rank | QR (Householder) | avoids the $\kappa^2$ of the normal equations (Theorem 3.7) | $2mn^2 - \tfrac23 n^3$ |
| Least squares, rank-deficient or nearly | small or zero $\sigma_i$ | SVD, truncated or pseudoinverse | only the SVD exposes which directions are undetermined | $\approx 4mn^2$ and up |
| Eigenvalues/vectors, symmetric | symmetric | spectral via tridiagonal QR algorithm or Jacobi | orthonormal eigenbasis, real eigenvalues, perfectly conditioned basis | $\tfrac43 n^3$ (values), $\approx 9n^3$ (vectors) |
| Eigenvalues/vectors, nonsymmetric | square | Schur, then eigenvectors from $T$ (Cor. 3.8b) | always exists, computed stably; eigenvectors only if needed | $\approx 10n^3$ ($T$), $25n^3$ (with $Q$) |
| Rank, null space, conditioning, low-rank approximation, PCA | any shape | SVD | Eckart–Young; $\kappa_2 = \sigma_1/\sigma_n$; rank is the count of non-negligible $\sigma_i$ | $\approx 20n^3$ (square, full) |
| Nearest orthogonal matrix; align point sets (Procrustes, registration) | square $M = YX^\top$ | polar via SVD, $U = WV^\top$ | Theorem 3.14 | one SVD of a $d\times d$ matrix |
| Matrix functions $e^{A}$, $A^{1/2}$, $\log A$ | symmetric: spectral; general: Schur | $f(A) = Qf(\Lambda)Q^\top$, or Schur–Parlett on $T$ | $f$ of a diagonal matrix is entrywise; of a triangular matrix, a recurrence | one spectral/Schur |
| Jordan form | — | never numerically | discontinuous in the entries (§3.8); use Schur | — |

Two rules bind the table together. First: **the cheaper the factorization, the more it assumes.** Cholesky ($\tfrac13 n^3$) assumes SPD; LU ($\tfrac23 n^3$) assumes square and nonsingular; QR ($\tfrac43 n^3$) assumes nothing but gives only the column space; SVD (around $20n^3$) assumes nothing and tells you everything. Second: **structure you fail to use is work you pay for and information you throw away.** Solving an SPD system by LU costs twice as much and does not tell you whether the matrix was SPD; solving least squares by the normal equations squares the condition number.

## 4. Worked examples

### Example 4.1: Cholesky by hand

$A = \begin{pmatrix} 4 & 2 & 2 \\ 2 & 5 & 3 \\ 2 & 3 & 6\end{pmatrix}$. Symmetric; is it SPD? Run (3.3) and find out.

Column 1: $l_{11} = \sqrt{a_{11}} = \sqrt4 = 2$. $l_{21} = a_{21}/l_{11} = 2/2 = 1$. $l_{31} = a_{31}/l_{11} = 2/2 = 1$.

Column 2: $l_{22} = \sqrt{a_{22} - l_{21}^2} = \sqrt{5 - 1} = 2$. $l_{32} = (a_{32} - l_{31}l_{21})/l_{22} = (3 - 1\cdot1)/2 = 1$.

Column 3: $l_{33} = \sqrt{a_{33} - l_{31}^2 - l_{32}^2} = \sqrt{6 - 1 - 1} = 2$.

$$L = \begin{pmatrix} 2 & 0 & 0 \\ 1 & 2 & 0 \\ 1 & 1 & 2\end{pmatrix}, \qquad LL^\top = \begin{pmatrix} 4 & 2 & 2 \\ 2 & 1 + 4 & 1 + 2 \\ 2 & 1 + 2 & 1 + 1 + 4\end{pmatrix} = A.\ ✓$$

All three radicands ($4, 4, 4$) were positive, so $A$ is SPD (Theorem 3.4), and they are the pivots: in $LDL^\top$ form $D = \operatorname{diag}(4, 4, 4)$ and $L_{\text{unit}} = L/2$; check $L_{\text{unit}}D^{1/2} = (L/2)(2I) = L$. Determinant: $\det A = (\det L)^2 = (2\cdot2\cdot2)^2 = 64$. (Cofactor expansion: $4(30 - 9) - 2(12 - 6) + 2(6 - 10) = 84 - 12 - 8 = 64$. ✓)

Solve $A\mathbf{x} = \mathbf{b}$ for $\mathbf{b} = (8, 10, 11)^\top$. Forward, $L\mathbf{y} = \mathbf{b}$: $y_1 = 8/2 = 4$; $y_2 = (10 - 1\cdot4)/2 = 3$; $y_3 = (11 - 1\cdot4 - 1\cdot3)/2 = 2$. Back, $L^\top\mathbf{x} = \mathbf{y}$: $x_3 = 2/2 = 1$; $x_2 = (3 - 1\cdot1)/2 = 1$; $x_1 = (4 - 1\cdot1 - 1\cdot1)/2 = 1$. So $\mathbf{x} = (1, 1, 1)^\top$, and $A\mathbf{x} = (8, 10, 11)^\top$. ✓

Now change one entry: $a_{33} = 1$ instead of $6$. Columns 1 and 2 are unchanged; column 3 gives radicand $1 - 1 - 1 = -1 < 0$. Cholesky fails at step 3, so the modified matrix is not SPD, and the pivots $(4, 4, -1)$ say it has exactly one negative eigenvalue (Sylvester's law of inertia, §3.1; proved in LA-12). The project confirms: eigenvalues $\approx 8.33, 2.45, -0.78$.

### Example 4.2: Schur form of a $2\times2$ by hand

$A = \begin{pmatrix} 1 & 2 \\ 3 & 0\end{pmatrix}$. Characteristic polynomial $\lambda^2 - \lambda - 6 = (\lambda - 3)(\lambda + 2)$: eigenvalues $3$ and $-2$. Not symmetric, so no spectral theorem. Not normal: $AA^\top = \begin{pmatrix} 5 & 3 \\ 3 & 9\end{pmatrix} \ne A^\top A = \begin{pmatrix} 10 & 2 \\ 2 & 4\end{pmatrix}$.

Follow the proof of Theorem 3.8 with $\lambda_1 = 3$. $(A - 3I)\mathbf{v} = \begin{pmatrix} -2 & 2 \\ 3 & -3\end{pmatrix}\mathbf{v} = \mathbf{0}$ gives $\mathbf{v} = (1, 1)^\top$; normalize: $\mathbf{q}_1 = \tfrac1{\sqrt2}(1, 1)^\top$. Complete to an orthonormal basis with $\mathbf{q}_2 = \tfrac1{\sqrt2}(-1, 1)^\top$:
$$Q = \frac1{\sqrt2}\begin{pmatrix} 1 & -1 \\ 1 & 1\end{pmatrix}\ (\text{a rotation by } 45°), \qquad T = Q^\top AQ.$$
Compute $AQ$ column by column: $A\mathbf{q}_1 = 3\mathbf{q}_1 = \tfrac1{\sqrt2}(3, 3)^\top$ and $A\mathbf{q}_2 = \tfrac1{\sqrt2}(-1 + 2,\ -3 + 0)^\top = \tfrac1{\sqrt2}(1, -3)^\top$. The entries of $T$ are the dot products $t_{ij} = \mathbf{q}_i^\top(A\mathbf{q}_j)$:
$$t_{11} = \mathbf{q}_1^\top(3\mathbf{q}_1) = 3, \quad t_{21} = \mathbf{q}_2^\top(3\mathbf{q}_1) = 0, \quad t_{12} = \tfrac12(1, 1)\cdot(1, -3) = -1, \quad t_{22} = \tfrac12(-1, 1)\cdot(1, -3) = -2,$$
$$T = \begin{pmatrix} 3 & -1 \\ 0 & -2\end{pmatrix}.$$
Upper triangular, eigenvalues on the diagonal, and $t_{21} = 0$ was automatic because $A\mathbf{q}_1 \parallel \mathbf{q}_1 \perp \mathbf{q}_2$: that is step 4 of the proof.

*Non-uniqueness, and what is invariant.* Had we chosen $\mathbf{q}_2 = \tfrac1{\sqrt2}(1, -1)^\top$, the same computation gives $T = \begin{pmatrix} 3 & 1 \\ 0 & -2\end{pmatrix}$. Had we started from $\lambda_1 = -2$, the diagonal would read $(-2, 3)$. What cannot change is $|t_{12}|$: by Corollary 3.8a, $t_{12}^2 = \|A\|_F^2 - (3^2 + (-2)^2) = (1 + 4 + 9 + 0) - 13 = 1$. The strictly upper part has size exactly $1$ in every Schur form of this $A$. SciPy's `schur` happens to return exactly our $Q$ and $T$ on the author's LAPACK build (another build may negate $\mathbf{q}_2$, and with it $t_{12}$); the project checks every number here, up to those signs.

*Eigenvector from $T$ (Corollary 3.8b).* For $\lambda = -2$, solve $(T + 2I)\mathbf{y} = \begin{pmatrix} 5 & -1 \\ 0 & 0\end{pmatrix}\mathbf{y} = \mathbf{0}$: set $y_2 = 1$, then $5y_1 - 1 = 0$, so $\mathbf{y} = (\tfrac15, 1)^\top \parallel (1, 5)^\top$. Then $\mathbf{x} = Q\mathbf{y} = \tfrac1{\sqrt2}(1 - 5,\ 1 + 5)^\top \parallel (-2, 3)^\top$. Check: $A(-2, 3)^\top = (-2 + 6,\ -6 + 0)^\top = (4, -6)^\top = -2\,(-2, 3)^\top$. ✓

### Example 4.3: Orthogonal Procrustes: align two triangles

A part drawing specifies three fiducial points; a scan measures them. Design points (columns of $X$) and measured points (columns of $Y$), already centered on their centroids:
$$X = \begin{pmatrix} -1 & 2 & -1 \\ -1 & -1 & 2\end{pmatrix}, \qquad Y = \begin{pmatrix} 3 & -1 & -2 \\ -3 & 2 & 1\end{pmatrix}.$$
(Check: every row sums to zero.) Find the rotation $R$ minimizing $\|RX - Y\|_F$.

*Step 1, the cross-covariance.* $M = YX^\top = \sum_i\mathbf{y}_i\mathbf{x}_i^\top$:
$$\begin{pmatrix} 3 \\ -3\end{pmatrix}(-1, -1) + \begin{pmatrix} -1 \\ 2\end{pmatrix}(2, -1) + \begin{pmatrix} -2 \\ 1\end{pmatrix}(-1, 2) = \begin{pmatrix} -3 & -3 \\ 3 & 3\end{pmatrix} + \begin{pmatrix} -2 & 1 \\ 4 & -2\end{pmatrix} + \begin{pmatrix} 2 & -4 \\ -1 & 2\end{pmatrix} = \begin{pmatrix} -3 & -6 \\ 6 & 3\end{pmatrix}.$$

*Step 2, the SVD of $M$.* $M^\top M = \begin{pmatrix} 9 + 36 & 18 + 18 \\ 18 + 18 & 36 + 9\end{pmatrix} = \begin{pmatrix} 45 & 36 \\ 36 & 45\end{pmatrix}$, eigenvalues $45 \pm 36 = 81, 9$ with eigenvectors $\mathbf{v}_1 = \tfrac1{\sqrt2}(1, 1)^\top$, $\mathbf{v}_2 = \tfrac1{\sqrt2}(1, -1)^\top$. So $\sigma_1 = 9$, $\sigma_2 = 3$, and $\mathbf{w}_i = M\mathbf{v}_i/\sigma_i$:
$$\mathbf{w}_1 = \frac1{9\sqrt2}\begin{pmatrix} -9 \\ 9\end{pmatrix} = \frac1{\sqrt2}\begin{pmatrix} -1 \\ 1\end{pmatrix}, \qquad \mathbf{w}_2 = \frac1{3\sqrt2}\begin{pmatrix} -3 + 6 \\ 6 - 3\end{pmatrix} = \frac1{\sqrt2}\begin{pmatrix} 1 \\ 1\end{pmatrix}.$$

*Step 3, the rotation.*
$$R = WV^\top = \frac12\begin{pmatrix} -1 & 1 \\ 1 & 1\end{pmatrix}\begin{pmatrix} 1 & 1 \\ 1 & -1\end{pmatrix} = \frac12\begin{pmatrix} 0 & -2 \\ 2 & 0\end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0\end{pmatrix},$$
a rotation by $+90°$, with $\det R = +1$, so no reflection fix is needed. $M$ itself is far from orthogonal; the SVD stripped out the stretch $V\Sigma V^\top = \begin{pmatrix} 6 & 3 \\ 3 & 6\end{pmatrix}$ and kept the rotation. (Check: $R\begin{pmatrix} 6 & 3 \\ 3 & 6\end{pmatrix} = \begin{pmatrix} -3 & -6 \\ 6 & 3\end{pmatrix} = M$. ✓)

*Step 4, the residual.* $RX = \begin{pmatrix} 1 & 1 & -2 \\ -1 & 2 & -1\end{pmatrix}$, so $Y - RX = \begin{pmatrix} 2 & -2 & 0 \\ -2 & 0 & 2\end{pmatrix}$ and $\|Y - RX\|_F^2 = 4 + 4 + 0 + 4 + 0 + 4 = 16$. Check against (3.6): $\|X\|_F^2 + \|Y\|_F^2 - 2(\sigma_1 + \sigma_2) = 12 + 28 - 24 = 16$. ✓ RMS error per point $\sqrt{16/3} \approx 2.31$. Before alignment, $\|Y - X\|_F^2 = 16 + 9 + 1 + 4 + 9 + 1 = 40$. The measured triangle is not a rotated copy of the design (if it were, the residual would be zero); the best rotation is still exactly $90°$.

### Example 4.4: Choosing the decomposition for three engineering tasks

**(a) A static stiffness solve, $K\mathbf{u} = \mathbf{f}$.** A finite-element stiffness matrix is symmetric, and for a properly constrained structure it is positive definite: $\tfrac12\mathbf{u}^\top K\mathbf{u}$ is the strain energy, positive for every nonzero displacement. It is also banded or sparse. **Cholesky.** Half the flops of LU, no pivoting so the band is preserved, and a failed pivot is a diagnosis: a zero pivot means a rigid-body mode you forgot to constrain. Concretely, two unit springs in a chain with the left end fixed give $K = \begin{pmatrix} 2 & -1 \\ -1 & 1\end{pmatrix}$, pivots $2$ and $1 - (-1)^2/2 = \tfrac12$: fine. Free the left end and $K = \begin{pmatrix} 1 & -1 \\ -1 & 1\end{pmatrix}$, pivots $1$ and $1 - 1 = 0$: Cholesky stops at step 2, which is exactly the statement that the whole chain can translate at no energy cost. If there are many load cases $\mathbf{f}_1, \mathbf{f}_2, \dots$, factor once and do two triangular solves ($2n^2$) per load.

**(b) Modal analysis, $K\boldsymbol\phi = \omega^2M\boldsymbol\phi$.** Two symmetric matrices, with $M$ SPD. The temptation is $M^{-1}K$, which is not symmetric and so throws away the spectral theorem. Instead, **Cholesky $M = LL^\top$, then spectral on $L^{-1}KL^{-\top}$**: substituting $\boldsymbol\phi = L^{-\top}\mathbf{y}$ gives $L^{-1}KL^{-\top}\mathbf{y} = \omega^2\mathbf{y}$, a *symmetric* standard eigenproblem (a congruence of $K$, so it keeps $K$'s inertia), with orthonormal $\mathbf{y}_i$ and therefore $M$-orthogonal modes: $\boldsymbol\phi_i^\top M\boldsymbol\phi_j = \mathbf{y}_i^\top L^{-1}LL^\top L^{-\top}\mathbf{y}_j = \mathbf{y}_i^\top\mathbf{y}_j = \delta_{ij}$. Numbers: $M = \operatorname{diag}(2, 1)$, $K = \begin{pmatrix} 3 & -1 \\ -1 & 1\end{pmatrix}$. Then $L = \operatorname{diag}(\sqrt2, 1)$ and $L^{-1}KL^{-\top} = \begin{pmatrix} 3/2 & -1/\sqrt2 \\ -1/\sqrt2 & 1\end{pmatrix}$, trace $\tfrac52$, determinant $\tfrac32 - \tfrac12 = 1$, so $\omega^2 = 2$ and $\omega^2 = \tfrac12$. Mode shapes from $(K - \omega^2M)\boldsymbol\phi = \mathbf{0}$: for $\omega^2 = \tfrac12$, $\begin{pmatrix} 2 & -1 \\ -1 & 1/2\end{pmatrix}\boldsymbol\phi = \mathbf{0}$ gives $\boldsymbol\phi_1 = (1, 2)^\top$; for $\omega^2 = 2$, $\begin{pmatrix} -1 & -1 \\ -1 & -1\end{pmatrix}\boldsymbol\phi = \mathbf{0}$ gives $\boldsymbol\phi_2 = (1, -1)^\top$. $M$-orthogonality: $\boldsymbol\phi_1^\top M\boldsymbol\phi_2 = 2\cdot1\cdot1 + 1\cdot2\cdot(-1) = 0$. ✓ This is CM-6's normal-mode theorem, and `scipy.linalg.eigh(K, M)` performs exactly this reduction.

**(c) A least-squares sensor calibration, $y_i \approx c_0 + c_1x_i + c_2x_i^2$.** The data matrix $A$ is $m\times3$ with $m \gg 3$: tall, no symmetry, no reason to be square. **Householder QR**, not the normal equations. If the $x_i$ span a narrow range, the columns $x$ and $x^2$ are nearly parallel, $\kappa(A)$ is large, and $\kappa(A^\top A) = \kappa(A)^2$ (Theorem 3.7) can pass $10^{14}$ while $\kappa(A)$ is a survivable $10^7$. If the fit is genuinely rank-deficient, or you want to know *which* combination of coefficients the data cannot determine, switch to the **SVD**: the right singular vector for the smallest $\sigma$ is that combination, and truncating it is the pseudoinverse solution. The project builds such an $A$ with $x \in [10, 10.1]$ and finds $\kappa(A) \approx 1.4\times10^7$, $\kappa(A^\top A) \approx 1.9\times10^{14}$.

## 5. Why it matters and where it connects

- **LA-13** is the "how it is computed" column of §3.9 done properly: floating-point error, backward stability of PLU and Householder QR, the QR algorithm converging to the Schur form, Golub–Kahan for the SVD, and Krylov methods for when $n$ is too large for any of this.
- **LA-11** proves the Jordan form that §3.8 declined to compute, and uses Schur's theorem directly in two of its proofs (Cayley–Hamilton among them: $p_A(T) = 0$ for triangular $T$ by a flag argument, so $p_A(A) = Qp_A(T)Q^* = 0$). Exercise C4 here gives a second proof of Cayley–Hamilton by continuity, and LA-11 §3.11's $\varepsilon$-demo is the same discontinuity as §3.8.
- **ST-3**: regression is least squares; the design matrix gets a QR, the covariance of the coefficient estimates is $\sigma^2(A^\top A)^{-1} = \sigma^2R^{-1}R^{-\top}$ from the same $R$, and ridge regression is the SVD with $\sigma_i \mapsto \sigma_i/(\sigma_i^2 + \lambda)$.
- **CM-6**: the normal-mode problem is Example 4.4(b); the Cholesky-then-spectral reduction is why normal modes are $M$-orthogonal.
- **LA-9 §3.11 (continuum mechanics)**: the polar decomposition $F = RU$ of the deformation gradient is Theorem 3.13 with $P$ renamed the right stretch tensor; LA-9 Example 4.5 works out a simple shear. No later module covers continuum kinematics.
- **CM-7**: the inertia tensor is symmetric, so its principal axes are the spectral theorem (§3.4); Euler's equations are written in that eigenbasis.
- **ML-4, ML-6**: PCA is the SVD row of the table; the Gram (kernel) matrix in ML-4 is SPD and is solved by Cholesky, whose diagonal also gives $\log\det K = 2\sum_i\log l_{ii}$ for the Gaussian-process likelihood.
- **F5**: sampling from $\mathcal N(\boldsymbol\mu, \Sigma)$ is $\boldsymbol\mu + L\mathbf{z}$ with $\Sigma = LL^\top$ and $\mathbf{z}$ standard normal: Cholesky as the square root of a covariance.
- **QM-4**: observables are Hermitian, so Theorem 3.10 gives their orthonormal eigenbasis with real eigenvalues; unitary time evolution has $|\lambda| = 1$, also Theorem 3.10. Normal operators are exactly those with a complete orthonormal set of eigenstates.
- **QC-1**: a single-qubit gate is a $2\times2$ unitary; Theorem 3.10 says it is a diagonal phase matrix in some orthonormal basis, which is why every gate is a rotation about some axis of the Bloch sphere.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: a $2\times2$ Cholesky, an indefinite $LDL^\top$, a flop budget, a second Schur form, a polar decomposition, and the decision table. B1–B8 are the module: B1 (Cholesky uniqueness), B2 (QR–Cholesky and $\kappa^2$), B5 (the $2\times2$ blocks of the real Schur form), B6 (polar decomposition and nearest orthogonal matrix) and B7 (extend `project.py` with $LDL^\top$ and inertia) are essential. C1–C4 go further: norms of normal matrices and transient growth, Procrustes with scale, Hadamard's inequality from Cholesky, and density of diagonalizable matrices from Schur.

## 7. Coding project

`project.py` builds the atlas as working code and checks every number in §4:

1. **Cholesky from (3.3)**, with a loud failure that names the pivot and its radicand when $A$ is not SPD. Checked against `scipy.linalg.cholesky` on Example 4.1 and on random SPD matrices; the failure diagnostics for the modified Example 4.1 matrix and for $\begin{pmatrix} 1 & 2 \\ 2 & 1\end{pmatrix}$ are shown, and the number of negative pivots is compared with the number of negative eigenvalues. LU without pivoting is run on random SPD matrices and never fails.
2. **Schur**: Example 4.2 by hand and by `scipy.linalg.schur`, reconstruction $A = QTQ^\top$, the invariant $|t_{12}| = 1$, the eigenvector from $T$, the real Schur form of a random $6\times6$. Then the normal-matrix theorem numerically: random $A = U\Lambda U^*$ with complex $\Lambda$ have Schur $T$ diagonal to $10^{-10}$, while random non-normal matrices have $\|\text{strict upper}(T)\|_F = \sqrt{\|A\|_F^2 - \sum|\lambda_i|^2}$. Figure: $|T|$ for a normal and a non-normal matrix.
3. **Polar and Procrustes**: `polar(A)` from the SVD with uniqueness checked ($P^2 = A^\top A$, $P$ SPD) and the nearest-orthogonal property tested against random orthogonal matrices; Example 4.3's triangles; a mirrored target that triggers the reflection fix; then registering a noisy scan of a 3D L-bracket to its CAD model, with a before/after figure and RMS error.
4. **Timing**: one SPD system of size $n = 200, \dots, 2000$ solved via LU, Cholesky, QR and SVD; measured times against the flop-count predictions on a log–log plot.
5. **`recommend(A, task)`**: inspects shape, symmetry, definiteness (by attempting Cholesky), rank and condition number, and prints which decomposition to use and why.
6. **Example 4.4 numbers**: stiffness pivots, modal frequencies and $M$-orthogonality via the Cholesky reduction and via `eigh(K, M)`, $\kappa(A^\top A) = \kappa(A)^2$ on the calibration matrix, the chain-of-springs pivots of Exercise B8, and the transient growth of Exercise C1.

**Build it yourself:** before reading the code, write `cholesky(A)` from (3.3) with the radicand check, and `procrustes(X, Y)` from Theorem 3.14 (center, $M = YX^\top$, SVD, $R = WV^\top$, determinant fix, $\mathbf{t} = \bar{\mathbf{y}} - R\bar{\mathbf{x}}$). Then compare. The "Your turn" cells ask for `sqrtm_spd` via the spectral decomposition, a proof-in-a-comment that unpivoted LU cannot fail on an SPD matrix, and Procrustes with a scale factor in 3D (Exercise C2).

## 8. Where you'll struggle

- **"I'll just call `np.linalg.solve`."** For a single dense nonsingular system, do. It costs you when (i) you solve with the same $A$ repeatedly and refactor each time: $k$ solves cost $k\cdot\tfrac23 n^3$ instead of $\tfrac23 n^3 + 2kn^2$, a factor of about $40$ at $n = 1000$, $k = 50$ (Exercise A3); (ii) $A$ is SPD and you pay double, lose the free definiteness test, and (in sparse codes) destroy the band structure; (iii) $A$ is ill-conditioned and `solve` returns garbage silently (it warns only when $A$ is numerically singular), whereas the SVD shows you the small $\sigma_i$ and lets you decide what to do about them; (iv) the problem is least squares and you form $A^\top A$ to make it square: that is Theorem 3.7 biting you.
- **Schur is the one nobody learns, and it is the key.** Courses go eigendecomposition → spectral theorem → SVD and skip Schur, which leaves the spectral theorem looking like a miracle and the eigendecomposition looking computable. Schur is what `eig` actually computes; the eigendecomposition is derived from it afterwards (Corollary 3.8b). Learn the proof of Theorem 3.8 until you can reproduce the block computation (3.5) from memory. It is the same deflation move as LA-8's proof of the spectral theorem, with "orthonormal eigenbasis" weakened to "orthonormal basis."
- **The diagonal-entry trick in Lemma 3.9.** $(TT^*)_{kk}$ is the squared length of *row* $k$ (sum over $j \ge k$); $(T^*T)_{kk}$ is the squared length of *column* $k$ (sum over $i \le k$). Three ways to get it wrong: swap row and column; forget to *shorten* the column sum using the zeros already proved in rows $1, \dots, k-1$; or compare the $(k,k)$ entries for all $k$ at once, which gives $n$ true equations and no conclusion. Do it in order, $k = 1, 2, \dots$, and after each step cross out the entries you have killed.
- **"SciPy's Schur form is different from mine, so I made a mistake."** Probably not. $Q$ and $T$ are not unique: the eigenvalue order and the signs (phases) of the columns of $Q$ are free, and the entries of the strictly upper part of $T$ change with them. What is invariant: the diagonal as a multiset, $\|T\|_F$, and $\sum_{i<j}|t_{ij}|^2$ (Corollary 3.8a). Check those.
- **Flop counts are not wall-clock times.** The timing figure will show Cholesky well below LU but not at exactly half, and the SVD at a different multiple than the flop count predicts. BLAS libraries are multithreaded and cache-blocked, and the factorizations block differently. The flop count tells you the exponent and the rough constant; only a measurement tells you the time.
- **Procrustes bookkeeping.** Three traps: forgetting to center both point sets (the rotation is then wrong, not only the translation); using $XY^\top$ in place of $YX^\top$, which hands you $R^\top$; and skipping the $\det = -1$ check, which quietly returns a reflection when the data are noisy and nearly planar. The residual identity (3.6) is your check: compute it before and after.
- **Cholesky "failing" on a matrix you know is SPD.** If $\kappa(A) \gtrsim 10^{16}$, a pivot that is positive in exact arithmetic can round to a negative number. The theorem is not wrong; your matrix is numerically singular. LA-13 shows how to shift or regularize.
- **The SVD versus the spectral decomposition of a symmetric matrix with negative eigenvalues.** They coincide only when $A$ is positive semidefinite. Otherwise $\sigma_i = |\lambda_i|$ and the sign moves into $\mathbf{u}_i = -\mathbf{v}_i$. If you take an SVD of a covariance-like matrix, every $\sigma_i$ looks fine, and the next step complains that the matrix is "not PSD," this is why.

## 9. Check yourself

1. Write each of LU, Cholesky, QR, spectral, Schur, SVD and polar as (basis change)(simple)(basis change), and name the relation each one uses (one-sided, congruence, similarity, unitary similarity, unitary equivalence).
2. Prove in three lines that $A = LDL^\top$ is SPD iff every $d_k > 0$.
3. Why does Cholesky need no pivoting? (One sentence about the Schur complement.)
4. Write the entry-matching formulas (3.3) for $l_{jj}$ and $l_{ij}$, and say which quantity under the square root is the pivot.
5. Why is $\tfrac13 n^3$ the cost of Cholesky and $\tfrac23 n^3$ the cost of LU?
6. State Schur's theorem with all hypotheses. Where exactly does the proof use that the field is $\mathbb{C}$?
7. Write the block form (3.5) of $Q_1^*AQ_1$ after one deflation step and explain why the first column is $\lambda_1\mathbf{e}_1$.
8. Prove that an upper triangular normal matrix is diagonal. Which entries of $TT^*$ and $T^*T$ do you compare first, and what does each sum run over?
9. Define a normal matrix; give three families of examples and one $2\times2$ non-example.
10. For invertible $A$, what are $U$ and $P$ in $A = UP$, and why is each forced?
11. For point sets $X, Y$ (columns), what matrix do you take the SVD of to find the best rotation, and what is the minimum residual in terms of $\|X\|_F$, $\|Y\|_F$ and the singular values?
12. You are handed (i) an SPD $5000\times5000$ system with 200 right-hand sides, (ii) a $10^6\times4$ least-squares problem, (iii) a $300\times300$ nonsymmetric matrix whose eigenvalues you need, (iv) a $50\times50$ matrix whose rank you doubt. Name the factorization for each and the cost in flops.

## 10. Pair with

- **Strang, *Linear Algebra and Learning from Data*, Part I** (I.1 for the factorizations as a group, I.8–I.9 for the SVD). Strang's framing of the factorizations as a family ($A = CR$, $LU$, $QR$, $Q\Lambda Q^\top$, $X\Lambda X^{-1}$, $U\Sigma V^\top$) is the ancestor of §2's table; he is better than these notes on the column-space view of each factorization and on why the SVD is the one for data.
- **Trefethen & Bau, *Numerical Linear Algebra*, lecture 5** ("More on the SVD"): the SVD as the master factorization, the cleanest proof that the $\sigma_i^2$ are eigenvalues of $A^\top A$, and the connection to norms. **Lecture 24** ("Eigenvalue problems"): why eigenvalue algorithms must be iterative, Schur as the target, and the symmetric/nonsymmetric split. Read lecture 24 before LA-13.
- **Golub & Van Loan, *Matrix Computations*, 4th ed., §4.2 (Cholesky), §7.1 and §7.4 (Schur and real Schur), §8.6 and Table 8.6.1 (SVD costs).** The source for every flop count quoted here, and for the algorithms LA-13 implements.
- **G. W. Stewart, "The decompositional approach to matrix computation," *Computing in Science & Engineering* 2(1), 50–59, 2000.** Ten pages on why §1's "factor, then do the easy thing" is the right way to think. Read it when you finish the module.
- **Horn & Johnson, *Matrix Analysis*, §2.3** (Schur's unitary triangularization, real Schur form) and **§2.5** (normal matrices): the full proofs of what this module stated.

## Appendix: the one-page cheat sheet

| Factorization | Form and shapes | Exists for | Unique? | Cost (flops) | Reach for it when | Not when |
|---|---|---|---|---|---|---|
| **PLU** | $PA = LU$; $P$ permutation, $L$ unit lower, $U$ upper; all $n\times n$ | every square $A$ | yes, given $P$ | $\tfrac23 n^3$ | general square solve; determinant; many right-hand sides (factor once) | $A$ is SPD (Cholesky); least squares; you need eigenvalues or rank |
| **$LDL^\top$** | $A = LDL^\top$; $L$ unit lower, $D$ diagonal | symmetric $A$ with nonzero pivots (any symmetric $A$ with symmetric pivoting and $2\times2$ blocks) | yes | $\tfrac13 n^3$ | symmetric indefinite solve; inertia from the signs of $D$ | $A$ not symmetric |
| **Cholesky** | $A = LL^\top$; $L$ lower with positive diagonal | SPD $A$, and only SPD | yes | $\tfrac13 n^3$ | SPD solve; definiteness test; $\log\det$; sampling $\mathcal N(\boldsymbol\mu, \Sigma)$; reducing $K\boldsymbol\phi = \omega^2M\boldsymbol\phi$ | $A$ not SPD (it will tell you); $A$ numerically singular |
| **QR** | $A = QR$; $A$ $m\times n$, $Q$ $m\times n$ orthonormal columns, $R$ $n\times n$ upper | every $A$ with $m \ge n$ | yes if full column rank and $\operatorname{diag}R > 0$ | $2mn^2 - \tfrac23 n^3$ | full-rank least squares; orthonormal basis of a column space | rank-deficient least squares (SVD); square solve (LU is half the cost) |
| **Eigen** | $A = X\Lambda X^{-1}$; $\Lambda$ diagonal | diagonalizable $A$ (e.g. distinct eigenvalues) | up to scaling and order | via Schur, $\approx 25n^3$ | $A^k$, $e^{tA}$, ODE systems, Markov chains, when $\kappa(X)$ is moderate | $A$ defective or $X$ ill-conditioned; $A$ symmetric (use spectral) |
| **Spectral** | $A = Q\Lambda Q^*$; $Q$ unitary, $\Lambda$ real diagonal | Hermitian $A$ (real symmetric) | up to order and bases of eigenspaces | $\tfrac43 n^3$ (values), $\approx 9n^3$ (with $Q$) | quadratic forms; principal axes; normal modes; PCA of a covariance; $f(A)$ for symmetric $A$ | $A$ not Hermitian |
| **Schur** | $A = QTQ^*$; $Q$ unitary, $T$ upper triangular (real form: $2\times2$ blocks) | every square complex $A$ | no: $Q$ phases and eigenvalue order are free; $\operatorname{diag}T$ and $\|\text{strict upper}\|_F$ invariant | $\approx 10n^3$ ($T$), $25n^3$ (with $Q$) | nonsymmetric eigenvalues; invariant subspaces; $f(A)$ for general $A$; proving things about all matrices | you need an actual eigenbasis (derive it from $T$ afterwards) |
| **Normal case** | $A = Q\Lambda Q^*$, $\Lambda$ complex diagonal | $AA^* = A^*A$ | as spectral | as Schur | unitary / Hermitian / skew-Hermitian: eigenvalues on the unit circle / real / imaginary | — |
| **SVD** | $A = U\Sigma V^\top$; $A$ $m\times n$, $U$ $m\times m$, $V$ $n\times n$ orthogonal, $\Sigma$ diagonal $\ge 0$ | every $A$ | $\Sigma$ yes; $U, V$ up to signs and rotations within equal $\sigma_i$ | $\approx 4mn^2$ ($\Sigma$ only); $\approx 20n^3$ square and full | rank; null and column space; $\kappa_2$; pseudoinverse; rank-deficient least squares; low-rank approximation; PCA | you only need to solve a well-conditioned square system (about $30\times$ the flops of LU) |
| **Polar** | $A = UP$; $U$ unitary, $P = (A^*A)^{1/2}$ PSD | every square $A$ | $P$ always; $U$ iff $A$ invertible | one SVD | nearest orthogonal matrix; Procrustes and registration; deformation $F = RU$; projecting onto $SO(n)$ | — |
| **Jordan** | $A = XJX^{-1}$; $J$ Jordan blocks | every square complex $A$ | up to block order | never compute | proofs about defective matrices (LA-11) | any floating-point computation (discontinuous; use Schur) |

**The question order.** Square? → Symmetric? → Positive definite (try Cholesky)? → Full rank (SVD if in doubt)? → Then the task: solve once (PLU) · solve often (PLU or Cholesky once, triangular solves after) · least squares (QR; SVD if rank-deficient) · eigenvalues (spectral if symmetric, Schur otherwise) · rank, conditioning, low rank (SVD) · nearest rotation (polar via SVD) · matrix function (spectral or Schur) · Jordan (never).
