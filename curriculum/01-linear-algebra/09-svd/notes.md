# LA-9: The Singular Value Decomposition and Low-Rank Approximation

**Time:** ~12 hours · **Prerequisites:** LA-8 (spectral theorem, positive semidefinite matrices, Rayleigh quotient) · **Unlocks:** LA-10, LA-14, ML-6, ST-3

> **Conventions in this module.** $A$ is a real $m \times n$ matrix of rank $r$. Its singular values are written in decreasing order, $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r > 0$, and we set $\sigma_{r+1} = \cdots = \sigma_{\min(m,n)} = 0$. In $A = U\Sigma V^\top$, $U$ is $m \times m$ orthogonal with columns $\mathbf{u}_i$ (left singular vectors), $V$ is $n \times n$ orthogonal with columns $\mathbf{v}_i$ (right singular vectors), and $\Sigma$ is $m \times n$ with $\sigma_i$ in position $(i,i)$ and zeros elsewhere. The *reduced* SVD keeps only the first $r$ columns: $U_r$ ($m \times r$), $\Sigma_r$ ($r \times r$), $V_r$ ($n \times r$). The sign of each pair $(\mathbf{u}_i, \mathbf{v}_i)$ is arbitrary but must be flipped *together*. $\|\mathbf{x}\|$ is the Euclidean norm of a vector; $\|A\|_2$ is the spectral (operator) norm and $\|A\|_F$ the Frobenius norm of a matrix; $\kappa$ means $\kappa_2 = \sigma_1/\sigma_n$ unless subscripted. A data matrix $X$ is $n \times p$: $n$ observations in rows, $p$ variables in columns, *centered* (column means zero); its sample covariance is $S = X^\top X/(n-1)$. The polar decomposition is written $A = QS$ with $Q$ orthogonal and $S$ symmetric positive semidefinite; continuum mechanics writes the same thing as $F = RU$, and *their* $U$ (the stretch) is *our* $S$, not our $U$. We flag this clash every time it appears.

## 0. Where this sits

You arrive with the spectral theorem from LA-8: a symmetric matrix has an orthonormal basis of eigenvectors, and in that basis it is diagonal. That theorem is about one basis and one square matrix. This module asks what survives when the matrix is rectangular, or square but not symmetric, so that eigenvectors are useless or absent. The answer is the singular value decomposition: use *two* orthonormal bases, one for the input space and one for the output space, and every matrix becomes diagonal. The proof is a half-page consequence of LA-8 applied to $A^\top A$.

The SVD is the capstone of the core linear algebra sequence (LA-1 through LA-9) and the tool you will reach for most often afterwards. LA-10 places it among the other decompositions. LA-13 explains how it is computed and why the obvious method is wrong. LA-14 differentiates functions of matrices whose answers come out in singular vectors. ST-3 uses the pseudoinverse and condition number for regression. ML-6 uses PCA, which is the SVD of a centered data matrix. CM-7 and the continuum-mechanics link in §3.11 use the polar decomposition, which is the SVD rearranged.

## 1. The problem that forced this idea

**Bilinear forms, 1873–1874.** The first appearance was not about data. Eugenio Beltrami (1873) and Camille Jordan (1874), working independently, asked a question about *bilinear forms*: an expression $f(\mathbf{x}, \mathbf{y}) = \mathbf{x}^\top A \mathbf{y}$ in two sets of variables. For a quadratic form $\mathbf{x}^\top A \mathbf{x}$ with $A$ symmetric, the spectral theorem (LA-8) says an orthogonal change of variables $\mathbf{x} = Q\boldsymbol{\xi}$ turns it into a sum of squares $\sum \lambda_i \xi_i^2$. Beltrami and Jordan asked the analogous question for the bilinear case, where $A$ is square but need not be symmetric: can *two* orthogonal substitutions, $\mathbf{x} = U\boldsymbol{\xi}$ and $\mathbf{y} = V\boldsymbol{\eta}$, reduce $f$ to the form $\sum_i \sigma_i\,\xi_i\eta_i$? Both proved they can, and both saw that the $\sigma_i^2$ are the eigenvalues of $A^\top A$ (equivalently of $AA^\top$). Sylvester rediscovered the result in 1889 and called the $\sigma_i$ "canonical multipliers." Erhard Schmidt (1907) carried it to integral equations, where the "matrix" is a function $K(s,t)$ and the decomposition is an infinite series; he also proved, in that setting, the approximation theorem we will meet as Eckart–Young. Hermann Weyl (1912) worked out how singular values move when the matrix is perturbed. The name "singular value" settled in during the 1930s. G. W. Stewart's 1993 survey *On the early history of the singular value decomposition* (SIAM Review) is the source for all of this and worth an hour.

**Psychometrics, 1936.** Carl Eckart and Gale Young published *The approximation of one matrix by another of lower rank* in the first volume of *Psychometrika*. The motivating data was a matrix of test scores: people down the rows, mental tests across the columns. Since Spearman (1904), psychologists had hoped such a matrix was "really" driven by a small number of underlying factors, which would make it close to a matrix of low rank. The practical question was: given the noisy observed matrix, what is the *nearest* matrix of rank $k$, and how far away is it? Eckart and Young proved that the answer is read off from the singular value decomposition: keep the $k$ largest singular values, discard the rest, and the error is $\sqrt{\sum_{i > k}\sigma_i^2}$. The theorem carries their names today (sometimes Eckart–Young–Mirsky, after Mirsky's 1960 extension to every norm that is invariant under orthogonal changes of basis). Every modern use of the SVD for compression, denoising, or dimension reduction is an application of their 1936 theorem.

**Computation, 1965.** For decades the SVD was a theorem, not a tool: forming $A^\top A$ and diagonalizing it squares the condition number and destroys the small singular values (the project shows this happening). Gene Golub and William Kahan (1965) found the way around it: reduce $A$ to bidiagonal form by orthogonal transformations applied directly to $A$, then compute the singular values of the bidiagonal matrix by an iteration that never forms $A^\top A$. Golub and Reinsch (1970) turned this into the algorithm that, with refinements, sits inside `numpy.linalg.svd` today. LA-13 derives it. Here you should only register the moral: *the derivation via $A^\top A$ is how you prove the SVD exists, and not how you compute it.*

**The modern problem.** You know from LA-7 and LA-8 what eigenvectors buy you: directions a square matrix merely scales. But a $100 \times 3$ matrix has no eigenvectors at all, since $A\mathbf{x} = \lambda\mathbf{x}$ does not even type-check when input and output live in different spaces. A square matrix that is not symmetric may have complex eigenvectors, non-orthogonal eigenvectors, or too few eigenvectors (LA-7, and the Jordan form in LA-11). So ask a weaker question that every matrix can answer. Instead of "which vectors does $A$ scale?" ask: **is there an orthonormal basis of the input space that $A$ maps to an orthogonal set in the output space?** If so, then in that input basis and the normalized output basis, $A$ is diagonal. The SVD says the answer is always yes, tells you how to find the bases, and the diagonal entries turn out to control everything you could want to know about $A$: its rank, its four fundamental subspaces, its norms, its sensitivity, its best low-rank approximations, and its (pseudo)inverse.

## 2. Intuition

### The unit circle goes to an ellipse

Take any $2 \times 2$ matrix and apply it to every unit vector. The result is an ellipse (possibly degenerate: a segment if the rank is 1, a point if $A = 0$). Here is the picture for $A = \begin{pmatrix} 3 & 0 \\ 4 & 5 \end{pmatrix}$, which we will decompose by hand in Example 4.1.

```
   input: unit circle                  output: the ellipse A(circle)
                                                  
          v1 = (1,1)/√2                              σ1 u1, length √45 ≈ 6.7
             ↗                                          ↗
      . ' ' ' ' .                              . - ' ' ' ' ' - .
    '             '                        . '                   ' .
   '   ↖ v2        '          A          .'      ↗ σ1 u1            '.
   '     (1,-1)/√2 '      -------->      '    ↙                       '
    '             '                       ' .   σ2 u2 (length √5 ≈ 2.2)
      ' . . . . '                            ' - . _ _ _ _ . - '
```

Three facts about this ellipse are the whole SVD:

1. **The ellipse has principal axes**, and they are perpendicular. The directions of the semi-axes are unit vectors $\mathbf{u}_1, \mathbf{u}_2$ in the output plane. These are the **left singular vectors**.
2. **The lengths of the semi-axes** are $\sigma_1 \ge \sigma_2 \ge 0$. These are the **singular values**. The longest semi-axis is the largest stretch $A$ applies to any unit vector, so $\sigma_1 = \max_{\|\mathbf{x}\|=1}\|A\mathbf{x}\|$. The shortest is the smallest stretch.
3. **The preimages of the axis endpoints** are unit vectors $\mathbf{v}_1, \mathbf{v}_2$ in the input plane, with $A\mathbf{v}_1 = \sigma_1\mathbf{u}_1$ and $A\mathbf{v}_2 = \sigma_2\mathbf{u}_2$. These are the **right singular vectors**. The remarkable part, which is not obvious from the picture, is that $\mathbf{v}_1 \perp \mathbf{v}_2$: the perpendicular axes of the output ellipse come from perpendicular directions on the input circle.

Put these together. $A$ sends the orthonormal pair $(\mathbf{v}_1, \mathbf{v}_2)$ to the orthogonal pair $(\sigma_1\mathbf{u}_1, \sigma_2\mathbf{u}_2)$. Any matrix that does this can be written as three simple steps:

- **Rotate** (or reflect) the input so that $\mathbf{v}_1, \mathbf{v}_2$ land on the coordinate axes. That is $V^\top$.
- **Stretch** along the axes by $\sigma_1$ and $\sigma_2$. That is the diagonal $\Sigma$.
- **Rotate** (or reflect) the axes to the directions $\mathbf{u}_1, \mathbf{u}_2$. That is $U$.

So $A = U\Sigma V^\top$: **rotate, stretch, rotate.** Every matrix, no matter how ugly, is a rotation, then an axis-aligned stretch, then another rotation. The stretch is the only part that changes lengths, and the two rotations say which directions get which stretch.

### Why $A^\top A$ is going to appear

Where are the axes of the ellipse? They are the directions on the unit circle where $\|A\mathbf{x}\|$ is largest and smallest. Now $\|A\mathbf{x}\|^2 = \mathbf{x}^\top A^\top A\,\mathbf{x}$, a quadratic form in $\mathbf{x}$ whose matrix $A^\top A$ is symmetric. LA-8 told you exactly where a quadratic form is largest and smallest on the unit sphere: at the eigenvectors of its matrix, with the eigenvalues as the extreme values (the Rayleigh quotient theorem). So the right singular vectors $\mathbf{v}_i$ *must* be eigenvectors of $A^\top A$, and $\sigma_i^2$ *must* be its eigenvalues. That is the whole derivation in §3.2, minus the bookkeeping. The orthogonality of the $\mathbf{v}_i$, which the picture does not make obvious, is now obvious: they are eigenvectors of a symmetric matrix.

### Higher dimensions and rectangular matrices

For an $m \times n$ matrix, the unit sphere in $\mathbb{R}^n$ goes to an ellipsoid in $\mathbb{R}^m$, but a possibly flattened one. If $A$ is $3 \times 2$ of rank 2, the unit circle in $\mathbb{R}^2$ lands as a flat ellipse lying in a plane in $\mathbb{R}^3$ (the column space). If $A$ is $2 \times 3$ of rank 2, the unit sphere in $\mathbb{R}^3$ is squashed onto a filled ellipse in $\mathbb{R}^2$, and one whole direction of the sphere (the null space) is crushed to the origin. In general the image is an $r$-dimensional ellipsoid with $r = \operatorname{rank}A$ semi-axes $\sigma_1\mathbf{u}_1, \dots, \sigma_r\mathbf{u}_r$, and $n - r$ input directions $\mathbf{v}_{r+1}, \dots, \mathbf{v}_n$ are sent to zero.

Now you can guess the theorem, and more. Guess (a): $A = U\Sigma V^\top$ exists for every real matrix. Guess (b): the $\mathbf{u}$'s that have nonzero $\sigma$'s span the column space, and the $\mathbf{v}$'s with zero $\sigma$ span the null space. Guess (c): the largest stretch $\|A\|_2$ is $\sigma_1$, and if $A$ is square and invertible, the worst-case amplification of errors in solving $A\mathbf{x} = \mathbf{b}$ is governed by how far the ellipse is from a circle, $\sigma_1/\sigma_n$. Guess (d): if you want to approximate $A$ by something simpler, drop the short axes first. All four guesses are correct, and §3 proves them.

## 3. Formal development

### 3.1 The statement, with shapes

**Definition 3.1 (Singular value decomposition).** Let $A$ be a real $m \times n$ matrix of rank $r$. A *singular value decomposition* of $A$ is a factorization

$$
A = U\Sigma V^\top, \tag{3.1}
$$

where $U$ is $m \times m$ orthogonal ($U^\top U = UU^\top = I_m$), $V$ is $n \times n$ orthogonal, and $\Sigma$ is $m \times n$ with entries $\Sigma_{ii} = \sigma_i$ for $i = 1, \dots, r$, where $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r > 0$, and all other entries zero. The $\sigma_i$ are the **singular values**, the columns $\mathbf{u}_1, \dots, \mathbf{u}_m$ of $U$ are the **left singular vectors**, and the columns $\mathbf{v}_1, \dots, \mathbf{v}_n$ of $V$ are the **right singular vectors**.

Since $\Sigma$ is "diagonal" but rectangular, its shape is worth drawing. For $m > n$ (tall) and $m < n$ (wide):

$$
\Sigma_{\,m > n} = \begin{pmatrix} \sigma_1 & & \\ & \ddots & \\ & & \sigma_n \\ 0 & \cdots & 0 \\ \vdots & & \vdots \end{pmatrix}, \qquad
\Sigma_{\,m < n} = \begin{pmatrix} \sigma_1 & & & 0 & \cdots \\ & \ddots & & \vdots & \\ & & \sigma_m & 0 & \cdots \end{pmatrix}.
$$

Multiplying out (3.1) column by column gives the statement you should actually keep in your head:

$$
A\mathbf{v}_i = \sigma_i\mathbf{u}_i \quad (i = 1, \dots, r), \qquad A\mathbf{v}_i = \mathbf{0} \quad (i = r+1, \dots, n). \tag{3.1$'$}
$$

*$A$ maps the orthonormal basis $\mathbf{v}_1, \dots, \mathbf{v}_n$ of $\mathbb{R}^n$ to the orthogonal set $\sigma_1\mathbf{u}_1, \dots, \sigma_r\mathbf{u}_r, \mathbf{0}, \dots, \mathbf{0}$.* Conversely, if you have orthonormal $\mathbf{v}$'s and $\mathbf{u}$'s satisfying (3.1$'$), you have an SVD: (3.1$'$) says $AV = U\Sigma$, and multiplying by $V^\top$ gives (3.1).

**The reduced SVD.** The zero rows and columns of $\Sigma$ multiply the columns $\mathbf{u}_{r+1}, \dots, \mathbf{u}_m$ and $\mathbf{v}_{r+1}, \dots, \mathbf{v}_n$ by nothing. Drop them:

$$
A = U_r\Sigma_rV_r^\top, \qquad U_r = \begin{pmatrix}\mathbf{u}_1 & \cdots & \mathbf{u}_r\end{pmatrix} \in \mathbb{R}^{m \times r}, \quad \Sigma_r = \operatorname{diag}(\sigma_1, \dots, \sigma_r), \quad V_r = \begin{pmatrix}\mathbf{v}_1 & \cdots & \mathbf{v}_r\end{pmatrix} \in \mathbb{R}^{n \times r}. \tag{3.2}
$$

Here $U_r^\top U_r = I_r$ and $V_r^\top V_r = I_r$ (orthonormal columns), but $U_rU_r^\top \ne I_m$ unless $r = m$: $U_r$ is not square and has no inverse. (3.1) and (3.2) are the same matrix $A$; (3.1) carries complete orthonormal bases of both spaces, (3.2) carries only the part that $A$ actually uses. Three shapes are in circulation and you need to recognize all of them:

| Name | $U$ | $\Sigma$ | $V$ | NumPy |
|---|---|---|---|---|
| Full | $m \times m$ | $m \times n$ | $n \times n$ | `np.linalg.svd(A)` returns `U, s, Vt` with `s` of length $\min(m,n)$ |
| Economy ("thin") | $m \times p$ | $p \times p$ | $n \times p$, $p = \min(m,n)$ | `np.linalg.svd(A, full_matrices=False)` |
| Reduced ("compact") | $m \times r$ | $r \times r$ | $n \times r$ | you slice to the nonzero `s` yourself |

NumPy returns $V^\top$ (called `Vt`), not $V$, and returns the singular values as a 1-D array, not as $\Sigma$. Economy and reduced coincide when $A$ has full rank $r = \min(m,n)$, which is the usual case for data; they differ for rank-deficient $A$, where the economy form carries singular values that are zero (or, numerically, tiny).

### 3.2 The derivation

**Theorem 3.2 (Existence).** Every real $m \times n$ matrix has a singular value decomposition.

▶ **Try it first.** $A^\top A$ is symmetric. Why is it positive semidefinite? Apply LA-8 to it and get an orthonormal eigenbasis $\mathbf{v}_i$ with eigenvalues $\lambda_i \ge 0$. Now compute $\|A\mathbf{v}_i\|^2$ and $(A\mathbf{v}_i)\cdot(A\mathbf{v}_j)$. That is the entire proof.

**Derivation.**

1. *$A^\top A$ is symmetric positive semidefinite.* It is $n \times n$ and $(A^\top A)^\top = A^\top (A^\top)^\top = A^\top A$. For any $\mathbf{x}$, $\mathbf{x}^\top A^\top A\mathbf{x} = (A\mathbf{x})^\top(A\mathbf{x}) = \|A\mathbf{x}\|^2 \ge 0$. By the spectral theorem (LA-8), there is an orthonormal basis $\mathbf{v}_1, \dots, \mathbf{v}_n$ of $\mathbb{R}^n$ and real numbers $\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_n$ with

$$
A^\top A\,\mathbf{v}_i = \lambda_i\mathbf{v}_i, \qquad \text{i.e.} \qquad A^\top A = V\Lambda V^\top, \tag{3.3}
$$

and $\lambda_i = \mathbf{v}_i^\top A^\top A\mathbf{v}_i = \|A\mathbf{v}_i\|^2 \ge 0$.

2. *Define the singular values.* Set $\sigma_i = \sqrt{\lambda_i} \ge 0$, so $\|A\mathbf{v}_i\| = \sigma_i$. Let $r$ be the number of strictly positive $\lambda_i$; then $\sigma_1 \ge \cdots \ge \sigma_r > 0$ and $\sigma_{r+1} = \cdots = \sigma_n = 0$.

3. *Define the left singular vectors and show they are orthonormal.* For $i \le r$, define

$$
\mathbf{u}_i = \frac{A\mathbf{v}_i}{\sigma_i}. \tag{3.4}
$$

Each $\mathbf{u}_i$ is a unit vector because $\|A\mathbf{v}_i\| = \sigma_i$. For the inner products,

$$
\mathbf{u}_i\cdot\mathbf{u}_j = \frac{(A\mathbf{v}_i)^\top(A\mathbf{v}_j)}{\sigma_i\sigma_j} = \frac{\mathbf{v}_i^\top (A^\top A\mathbf{v}_j)}{\sigma_i\sigma_j} = \frac{\mathbf{v}_i^\top(\lambda_j\mathbf{v}_j)}{\sigma_i\sigma_j} = \frac{\lambda_j}{\sigma_i\sigma_j}\,(\mathbf{v}_i\cdot\mathbf{v}_j) = \frac{\sigma_j^2}{\sigma_i\sigma_j}\,\delta_{ij} = \delta_{ij}.
$$

This is the step that carries the whole theorem. The images $A\mathbf{v}_i$ of orthonormal eigenvectors of $A^\top A$ are automatically orthogonal, because the inner product of two images is the quadratic form of $A^\top A$ evaluated on two *eigenvectors*, and eigenvectors of a symmetric matrix for different eigenvalues are orthogonal. Nothing about $A$ being square or symmetric was used.

4. *The defining relations.* For $i \le r$, (3.4) says $A\mathbf{v}_i = \sigma_i\mathbf{u}_i$. For $i > r$, $\|A\mathbf{v}_i\| = \sigma_i = 0$, so $A\mathbf{v}_i = \mathbf{0}$. This is (3.1$'$).

5. *Complete the output basis.* The vectors $\mathbf{u}_1, \dots, \mathbf{u}_r$ are orthonormal in $\mathbb{R}^m$, so $r \le m$. Extend them to an orthonormal basis $\mathbf{u}_1, \dots, \mathbf{u}_m$ of $\mathbb{R}^m$: append the standard basis vectors $\mathbf{e}_1, \dots, \mathbf{e}_m$ to the list and run Gram–Schmidt (LA-5), discarding the vectors that come out zero. Put the $\mathbf{u}_i$ into the columns of $U$, the $\mathbf{v}_i$ into the columns of $V$, and the $\sigma_i$ into the diagonal of the $m \times n$ matrix $\Sigma$.

6. *Assemble.* Column $j$ of $AV$ is $A\mathbf{v}_j$. Column $j$ of $U\Sigma$ is $U$ times column $j$ of $\Sigma$, which is $\sigma_j\mathbf{e}_j$ for $j \le r$ and $\mathbf{0}$ for $j > r$; so column $j$ of $U\Sigma$ is $\sigma_j\mathbf{u}_j$ for $j \le r$ and $\mathbf{0}$ otherwise. By step 4 the columns agree:

$$
AV = U\Sigma. \tag{3.5}
$$

Multiply on the right by $V^\top$ and use $VV^\top = I_n$: $A = U\Sigma V^\top$. ∎

7. *The number $r$ is the rank.* Multiplying by an invertible matrix does not change rank (LA-2), so $\operatorname{rank}A = \operatorname{rank}(AV) = \operatorname{rank}(U\Sigma) = \operatorname{rank}\Sigma = r$, since $\Sigma$ has exactly $r$ nonzero rows and they are independent. Theorem 3.3 below gives the same fact with more detail.

■ **What just happened.** We could not diagonalize $A$ directly, so we diagonalized the symmetric matrix $A^\top A$ instead, which LA-8 allows. The payoff came in step 3: the vectors $A\mathbf{v}_i$ turned out to be orthogonal with no further work, and normalizing them produced the output basis. The SVD is the spectral theorem for $A^\top A$ plus the one-line observation that eigenvectors of $A^\top A$ have orthogonal images. If you remember nothing else of this proof, remember step 3.

### 3.3 The other Gram matrix, and uniqueness

**$AA^\top$.** From (3.1), $AA^\top = U\Sigma V^\top V\Sigma^\top U^\top = U(\Sigma\Sigma^\top)U^\top$, and $\Sigma\Sigma^\top$ is the $m \times m$ diagonal matrix $\operatorname{diag}(\sigma_1^2, \dots, \sigma_r^2, 0, \dots, 0)$. So the left singular vectors are eigenvectors of $AA^\top$, with eigenvalues $\sigma_i^2$ for $i \le r$ and $0$ for $i > r$. The nonzero eigenvalues of $A^\top A$ ($n \times n$) and $AA^\top$ ($m \times m$) are the same numbers, $\sigma_1^2, \dots, \sigma_r^2$; the two matrices differ only in how many zeros they carry. When $m \ll n$ it is cheaper to start the derivation from $AA^\top$ (find the $\mathbf{u}_i$ first, then $\mathbf{v}_i = A^\top\mathbf{u}_i/\sigma_i$); Exercise B2 does this.

A trap: computing $V$ from an eigensolver applied to $A^\top A$ and $U$ *separately* from an eigensolver applied to $AA^\top$ gives vectors with unrelated signs, and $U\Sigma V^\top$ will then not equal $A$. The left vectors must come from the right vectors through (3.4), or vice versa.

**Uniqueness.** The singular values are unique: they are the square roots of the eigenvalues of $A^\top A$, which is determined by $A$. The singular vectors are not unique, and the freedom is exactly this:

- If $\sigma_i > 0$ is a *simple* singular value (it appears once), then $\mathbf{v}_i$ spans a one-dimensional eigenspace of $A^\top A$ and is determined up to sign. Then $\mathbf{u}_i = A\mathbf{v}_i/\sigma_i$ is determined by $\mathbf{v}_i$. The only freedom is the simultaneous flip $(\mathbf{u}_i, \mathbf{v}_i) \to (-\mathbf{u}_i, -\mathbf{v}_i)$.
- If $\sigma_i > 0$ is *repeated*, the span of the corresponding $\mathbf{v}$'s (an eigenspace of $A^\top A$) is determined, but any orthonormal basis of it may be used; the $\mathbf{u}$'s then follow from (3.4). The identity matrix is the extreme case: $I = QIQ^\top$ for every orthogonal $Q$.
- The vectors $\mathbf{v}_{r+1}, \dots, \mathbf{v}_n$ can be any orthonormal basis of $N(A)$, and $\mathbf{u}_{r+1}, \dots, \mathbf{u}_m$ any orthonormal basis of $N(A^\top)$ (Theorem 3.3).

So two correct programs can return different $U$ and $V$ for the same $A$ and both be right. Compare $U\Sigma V^\top$ to $A$, or compare subspaces, never the signs of individual columns. §8 returns to this.

### 3.4 The four fundamental subspaces

**Theorem 3.3 (The SVD exhibits the four subspaces).** Let $A = U\Sigma V^\top$ be an SVD of the $m \times n$ matrix $A$ of rank $r$. Then

- (a) $C(A) = \operatorname{span}\{\mathbf{u}_1, \dots, \mathbf{u}_r\}$ (column space),
- (b) $N(A^\top) = \operatorname{span}\{\mathbf{u}_{r+1}, \dots, \mathbf{u}_m\}$ (left null space),
- (c) $C(A^\top) = \operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_r\}$ (row space),
- (d) $N(A) = \operatorname{span}\{\mathbf{v}_{r+1}, \dots, \mathbf{v}_n\}$ (null space).

▶ **Try it first.** Write an arbitrary $\mathbf{x} \in \mathbb{R}^n$ in the basis of $\mathbf{v}$'s and apply (3.1$'$). Read off where $A\mathbf{x}$ lives and when it is zero. Then notice $A^\top = V\Sigma^\top U^\top$ is also an SVD.

**Proof.**

1. *Expansion.* Since $\mathbf{v}_1, \dots, \mathbf{v}_n$ is an orthonormal basis of $\mathbb{R}^n$, any $\mathbf{x}$ can be written $\mathbf{x} = \sum_{j=1}^n c_j\mathbf{v}_j$ with $c_j = \mathbf{v}_j\cdot\mathbf{x}$. By linearity and (3.1$'$),

$$
A\mathbf{x} = \sum_{j=1}^n c_j A\mathbf{v}_j = \sum_{j=1}^r c_j\sigma_j\mathbf{u}_j.
$$

2. *(a).* The right side is in $\operatorname{span}\{\mathbf{u}_1, \dots, \mathbf{u}_r\}$, so $C(A) \subseteq \operatorname{span}\{\mathbf{u}_1, \dots, \mathbf{u}_r\}$. Conversely, each $\mathbf{u}_i$ with $i \le r$ equals $A(\mathbf{v}_i/\sigma_i)$ and so lies in $C(A)$; since $C(A)$ is a subspace it contains their span. Equality.

3. *(d).* $A\mathbf{x} = \mathbf{0}$ iff $\sum_{j \le r} (c_j\sigma_j)\mathbf{u}_j = \mathbf{0}$. The $\mathbf{u}_1, \dots, \mathbf{u}_r$ are orthonormal, hence linearly independent, so this holds iff every coefficient $c_j\sigma_j$ is zero, iff $c_j = 0$ for all $j \le r$ (because $\sigma_j > 0$). That says $\mathbf{x} = \sum_{j > r}c_j\mathbf{v}_j$, i.e. $\mathbf{x} \in \operatorname{span}\{\mathbf{v}_{r+1}, \dots, \mathbf{v}_n\}$. Both directions of the "iff" were kept, so this is equality.

4. *(c) and (b).* Transpose (3.1): $A^\top = V\Sigma^\top U^\top$. Here $V$ and $U$ are orthogonal and $\Sigma^\top$ is $n \times m$ with the same $\sigma_1, \dots, \sigma_r$ on its diagonal, so this is an SVD of $A^\top$ with left singular vectors $\mathbf{v}_i$ and right singular vectors $\mathbf{u}_i$. Apply (a) to $A^\top$: $C(A^\top) = \operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_r\}$. Apply (d) to $A^\top$: $N(A^\top) = \operatorname{span}\{\mathbf{u}_{r+1}, \dots, \mathbf{u}_m\}$. ∎

■ **What just happened.** Everything you proved about the four subspaces in LA-3 and LA-5 falls out as bookkeeping. Rank–Nullity: $\dim N(A) = n - r$ because that is how many $\mathbf{v}$'s are left over. Row rank equals column rank: both are $r$. $C(A) \perp N(A^\top)$ and together they fill $\mathbb{R}^m$: because $\mathbf{u}_1, \dots, \mathbf{u}_m$ is an orthonormal basis split into two pieces. Same for $C(A^\top) \perp N(A)$ in $\mathbb{R}^n$.

**The thesis of this module.** *The SVD is the four-subspace picture with orthonormal bases chosen so that $A$ is diagonal.* In LA-3 the four subspaces were sets. In LA-5 they came in two orthogonal pairs. The SVD picks a specific orthonormal basis inside each: $\mathbf{v}_1, \dots, \mathbf{v}_r$ for the row space, $\mathbf{v}_{r+1}, \dots, \mathbf{v}_n$ for the null space, $\mathbf{u}_1, \dots, \mathbf{u}_r$ for the column space, $\mathbf{u}_{r+1}, \dots, \mathbf{u}_m$ for the left null space. And these bases are matched: $A$ sends the $i$-th row-space basis vector to $\sigma_i$ times the $i$-th column-space basis vector, and kills the null space. In the language of LA-4, $\Sigma = U^\top AV$ is the matrix of the linear map $\mathbf{x} \mapsto A\mathbf{x}$ with respect to the input basis $V$ and the output basis $U$. Every linear map between finite-dimensional inner product spaces is diagonal in some pair of orthonormal bases. The rest of this module is a list of things you can read off a diagonal matrix.

### 3.5 $A$ as a sum of rank-one pieces, ordered by importance

Recall from LA-1 that a product $XY$ can be computed as a sum of (column of $X$)(row of $Y$) outer products. Apply this to $(U\Sigma)V^\top$. Column $i$ of $U\Sigma$ is $\sigma_i\mathbf{u}_i$ (step 6 above), and row $i$ of $V^\top$ is $\mathbf{v}_i^\top$, so

$$
A = \sum_{i=1}^{r}\sigma_i\,\mathbf{u}_i\mathbf{v}_i^\top. \tag{3.6}
$$

Each term $\mathbf{u}_i\mathbf{v}_i^\top$ is an $m \times n$ matrix of rank one (every column is a multiple of $\mathbf{u}_i$) and of unit Frobenius norm, since $\|\mathbf{u}_i\mathbf{v}_i^\top\|_F^2 = \sum_{a,b}u_{ia}^2v_{ib}^2 = \|\mathbf{u}_i\|^2\|\mathbf{v}_i\|^2 = 1$. The terms are mutually orthogonal in the Frobenius inner product $\langle X, Y\rangle_F = \operatorname{tr}(X^\top Y)$:

$$
\langle\mathbf{u}_i\mathbf{v}_i^\top, \mathbf{u}_j\mathbf{v}_j^\top\rangle_F = \operatorname{tr}(\mathbf{v}_i\mathbf{u}_i^\top\mathbf{u}_j\mathbf{v}_j^\top) = (\mathbf{u}_i\cdot\mathbf{u}_j)\operatorname{tr}(\mathbf{v}_i\mathbf{v}_j^\top) = (\mathbf{u}_i\cdot\mathbf{u}_j)(\mathbf{v}_i\cdot\mathbf{v}_j) = \delta_{ij}.
$$

So (3.6) writes $A$ as a combination of $r$ orthonormal rank-one matrices with coefficients $\sigma_1 \ge \sigma_2 \ge \cdots$. The first term is the single rank-one matrix that captures the most of $A$; the second is the best rank-one correction to the first; and so on. The truncation

$$
A_k = \sum_{i=1}^{k}\sigma_i\,\mathbf{u}_i\mathbf{v}_i^\top = U_k\Sigma_kV_k^\top \qquad (k \le r)
$$

has rank exactly $k$, and §3.7 proves it is the best rank-$k$ approximation of $A$. For the record, storing $A_k$ costs $k(m + n + 1)$ numbers instead of $mn$; for a $1000 \times 1000$ matrix and $k = 20$ that is $40\,020$ numbers instead of a million.

### 3.6 Matrix norms from singular values

**Definition 3.4.** The **spectral norm** (operator norm, 2-norm) of $A$ is $\|A\|_2 = \max_{\|\mathbf{x}\| = 1}\|A\mathbf{x}\|$, the largest factor by which $A$ stretches any vector. The **Frobenius norm** is $\|A\|_F = \sqrt{\sum_{i,j}a_{ij}^2} = \sqrt{\operatorname{tr}(A^\top A)}$, the Euclidean length of $A$ viewed as a list of $mn$ numbers. (LA-1 Exercise C3 introduced both and asked whether the smallest constant $K$ in $\|A\mathbf{x}\| \le K\|\mathbf{x}\|$ is $\|A\|_F$. It is not. Here is what it is.)

**Theorem 3.5.** $\|A\|_2 = \sigma_1$. If $m \ge n$, also $\min_{\|\mathbf{x}\|=1}\|A\mathbf{x}\| = \sigma_n$.

▶ **Try it first.** Expand a unit vector $\mathbf{x}$ in the $\mathbf{v}$-basis, apply (3.1$'$), and compute $\|A\mathbf{x}\|^2$ using the orthonormality of the $\mathbf{u}$'s.

**Derivation.**

1. Let $\|\mathbf{x}\| = 1$ and write $\mathbf{x} = \sum_{j=1}^n c_j\mathbf{v}_j$. Because the $\mathbf{v}_j$ are orthonormal, $\|\mathbf{x}\|^2 = \sum_j c_j^2 = 1$.
2. $A\mathbf{x} = \sum_{j \le r}c_j\sigma_j\mathbf{u}_j$, and because the $\mathbf{u}_j$ are orthonormal, $\|A\mathbf{x}\|^2 = \sum_{j \le r}\sigma_j^2c_j^2$.
3. Since $\sigma_j \le \sigma_1$ for every $j$: $\|A\mathbf{x}\|^2 \le \sigma_1^2\sum_j c_j^2 = \sigma_1^2$. Equality holds at $\mathbf{x} = \mathbf{v}_1$ (all $c_j = 0$ except $c_1 = 1$), where $\|A\mathbf{v}_1\| = \sigma_1$. So the maximum is $\sigma_1$.
4. For the minimum when $m \ge n$: all $n$ of the $\sigma_j$ (including zeros, if any) appear, $\|A\mathbf{x}\|^2 = \sum_{j=1}^n\sigma_j^2c_j^2 \ge \sigma_n^2\sum_j c_j^2 = \sigma_n^2$, with equality at $\mathbf{x} = \mathbf{v}_n$. ∎

■ **What just happened.** In the $\mathbf{v}$-basis, $\|A\mathbf{x}\|^2$ is a weighted average of the $\sigma_j^2$ with weights $c_j^2$ summing to one. An average is at most the largest entry and at least the smallest. The spectral norm is the longest semi-axis of the ellipsoid, and (for $m \ge n$) the smallest singular value is the shortest.

$$
\|A\|_2 = \sigma_1, \qquad \min_{\|\mathbf{x}\|=1}\|A\mathbf{x}\| = \sigma_n \ (m \ge n). \tag{3.7}
$$

**Theorem 3.6.** $\|A\|_F^2 = \sigma_1^2 + \cdots + \sigma_r^2$.

**Derivation.** We need the trace identity $\operatorname{tr}(XY) = \operatorname{tr}(YX)$ for $X$ ($p \times q$) and $Y$ ($q \times p$): both sides equal $\sum_{i=1}^p\sum_{j=1}^q x_{ij}y_{ji}$. Then

$$
\|A\|_F^2 = \operatorname{tr}(A^\top A) = \operatorname{tr}(V\Sigma^\top\Sigma V^\top) = \operatorname{tr}(\Sigma^\top\Sigma\,V^\top V) = \operatorname{tr}(\Sigma^\top\Sigma) = \sum_{i=1}^r\sigma_i^2,
$$

moving the leading $V$ to the back with the trace identity and using $V^\top V = I$. ∎

$$
\|A\|_F^2 = \sum_{i=1}^{r}\sigma_i^2. \tag{3.8}
$$

Two consequences. First, $\sigma_1 \le \|A\|_F \le \sqrt{r}\,\sigma_1$, with $\|A\|_2 = \|A\|_F$ exactly when $A$ has rank one; this settles LA-1 C3. Second, both norms are **orthogonally invariant**: for orthogonal $Q$ ($m \times m$) and $W$ ($n \times n$), $QAW = (QU)\Sigma(W^\top V)^\top$ is an SVD of $QAW$ with the same $\Sigma$, so $\|QAW\|_2 = \|A\|_2$ and $\|QAW\|_F = \|A\|_F$. Rotating the input or output space changes neither the biggest stretch nor the total energy. This invariance is what makes the Eckart–Young proof work.

### 3.7 Eckart–Young: the best low-rank approximation

**Theorem 3.7 (Eckart–Young, Frobenius norm).** Let $A$ have SVD $A = \sum_{i=1}^r\sigma_i\mathbf{u}_i\mathbf{v}_i^\top$ and let $k < r$. For every $m \times n$ matrix $B$ with $\operatorname{rank}B \le k$,

$$
\|A - B\|_F^2 \ \ge\ \sum_{i=k+1}^{r}\sigma_i^2, \tag{3.9}
$$

with equality for $B = A_k = \sum_{i=1}^k\sigma_i\mathbf{u}_i\mathbf{v}_i^\top$. So $A_k$ is a best rank-$k$ approximation to $A$ in the Frobenius norm, and the error is $\|A - A_k\|_F = \sqrt{\sigma_{k+1}^2 + \cdots + \sigma_r^2}$.

The statement is easy to remember; the proof is not, and most courses skip it. We will not. The strategy: any rank-$k$ $B$ lives in a $k$-dimensional column space; project $A$ onto that space and onto its complement; the complement part of the error is unavoidable; then show no $k$-dimensional subspace can capture more of $A$'s energy than $\mathbf{u}_1, \dots, \mathbf{u}_k$ do.

▶ **Try it first.** Show the easy half: $\|A - A_k\|_F^2 = \sum_{i > k}\sigma_i^2$. Then convince yourself that for *any* orthogonal projector $P$ of rank $k$, $\|PA\|_F^2 = \sum_j\|A^\top\mathbf{q}_j\|^2$ where $\mathbf{q}_1, \dots, \mathbf{q}_k$ is an orthonormal basis of the range of $P$.

**Proof.**

1. *$A_k$ achieves the bound.* $A - A_k = \sum_{i > k}\sigma_i\mathbf{u}_i\mathbf{v}_i^\top$ is itself written in SVD form, with singular values $\sigma_{k+1}, \dots, \sigma_r$. By (3.8), $\|A - A_k\|_F^2 = \sum_{i > k}\sigma_i^2$.

2. *Set up the projector.* Let $B$ have rank $p \le k$. Let $\mathbf{q}_1, \dots, \mathbf{q}_p$ be an orthonormal basis of $C(B)$ (Gram–Schmidt, LA-5) and $P = \sum_{j=1}^p\mathbf{q}_j\mathbf{q}_j^\top$ the orthogonal projector onto $C(B)$. Then $P^2 = P = P^\top$, $PB = B$ (every column of $B$ is in $C(B)$) and $(I - P)B = 0$.

3. *Split the error by Pythagoras.* For any matrix $M$, the columns of $PM$ lie in $C(P)$ and the columns of $(I - P)M$ lie in $C(P)^\perp$, so column by column $\|M\mathbf{e}_j\|^2 = \|PM\mathbf{e}_j\|^2 + \|(I-P)M\mathbf{e}_j\|^2$; summing over columns, $\|M\|_F^2 = \|PM\|_F^2 + \|(I-P)M\|_F^2$. Apply this to $M = A - B$:

$$
\|A - B\|_F^2 = \|P(A-B)\|_F^2 + \|(I-P)(A-B)\|_F^2 \ \ge\ \|(I-P)(A-B)\|_F^2 = \|(I-P)A\|_F^2,
$$

using $(I-P)B = 0$ in the last step. Whatever $B$ is, the error is at least the part of $A$ that sticks out of $C(B)$.

4. *Rewrite the unavoidable part.* Apply the Pythagoras split to $M = A$: $\|(I-P)A\|_F^2 = \|A\|_F^2 - \|PA\|_F^2$. So far: $\|A - B\|_F^2 \ge \|A\|_F^2 - \|PA\|_F^2$, and it remains to show $\|PA\|_F^2 \le \sum_{i \le k}\sigma_i^2$: *no $k$-dimensional subspace captures more of $A$ than the top $k$ left singular vectors.*

5. *Express $\|PA\|_F^2$ through the $\mathbf{q}_j$.* Using the trace identity and $P^\top P = P$,
$\|PA\|_F^2 = \operatorname{tr}(A^\top P^\top PA) = \operatorname{tr}(A^\top PA) = \sum_{j=1}^p\operatorname{tr}(A^\top\mathbf{q}_j\mathbf{q}_j^\top A) = \sum_{j=1}^p\|A^\top\mathbf{q}_j\|^2.$

6. *Expand each $\mathbf{q}_j$ in the $\mathbf{u}$-basis.* $\mathbf{q}_j = \sum_{i=1}^m(\mathbf{u}_i\cdot\mathbf{q}_j)\mathbf{u}_i$. From the SVD of $A^\top$ (§3.4, step 4), $A^\top\mathbf{u}_i = \sigma_i\mathbf{v}_i$ for $i \le r$ and $A^\top\mathbf{u}_i = \mathbf{0}$ for $i > r$. So $A^\top\mathbf{q}_j = \sum_{i \le r}(\mathbf{u}_i\cdot\mathbf{q}_j)\sigma_i\mathbf{v}_i$ and, by orthonormality of the $\mathbf{v}_i$, $\|A^\top\mathbf{q}_j\|^2 = \sum_{i \le r}\sigma_i^2(\mathbf{u}_i\cdot\mathbf{q}_j)^2$. Summing over $j$:

$$
\|PA\|_F^2 = \sum_{i=1}^r\sigma_i^2c_i, \qquad c_i := \sum_{j=1}^p(\mathbf{u}_i\cdot\mathbf{q}_j)^2 = \|P\mathbf{u}_i\|^2.
$$

7. *Constraints on the weights $c_i$.* Each $c_i = \|P\mathbf{u}_i\|^2$ satisfies $0 \le c_i \le \|\mathbf{u}_i\|^2 = 1$, because an orthogonal projection never lengthens a vector. And their total is
$\sum_{i=1}^m c_i = \sum_{j=1}^p\sum_{i=1}^m(\mathbf{u}_i\cdot\mathbf{q}_j)^2 = \sum_{j=1}^p\|\mathbf{q}_j\|^2 = p \le k,$
using that the $\mathbf{u}_i$ are an orthonormal basis of $\mathbb{R}^m$ (Parseval). So $\sum_{i \le r}c_i \le k$ as well.

8. *The weighted sum is maximized by loading the top $k$ slots.* We have weights $0 \le c_i \le 1$ with $\sum_{i\le r}c_i \le k$, and we want to bound $\sum_{i \le r}\sigma_i^2c_i$ where $\sigma_1^2 \ge \sigma_2^2 \ge \cdots$. Split the sum at $k$:

$$
\sum_{i > k}\sigma_i^2c_i \ \le\ \sigma_k^2\sum_{i > k}c_i \ \le\ \sigma_k^2\sum_{i \le k}(1 - c_i),
$$

where the first step uses $\sigma_i^2 \le \sigma_k^2$ for $i > k$ and the second uses $\sum_{i > k}c_i \le k - \sum_{i \le k}c_i = \sum_{i \le k}(1 - c_i)$. Therefore

$$
\sum_{i \le r}\sigma_i^2c_i \ \le\ \sum_{i \le k}\bigl[\sigma_i^2c_i + \sigma_k^2(1 - c_i)\bigr] \ \le\ \sum_{i \le k}\bigl[\sigma_i^2c_i + \sigma_i^2(1 - c_i)\bigr] = \sum_{i \le k}\sigma_i^2,
$$

where the middle step uses $\sigma_k^2 \le \sigma_i^2$ for $i \le k$ and $1 - c_i \ge 0$. In words: a budget of $k$ units of weight, at most one unit per slot, buys the most when spent entirely on the $k$ largest $\sigma_i^2$.

9. *Combine.* $\|A - B\|_F^2 \ \ge\ \|A\|_F^2 - \|PA\|_F^2 \ \ge\ \sum_{i \le r}\sigma_i^2 - \sum_{i \le k}\sigma_i^2 = \sum_{i > k}\sigma_i^2$, using (3.8). With step 1, this is (3.9). ∎

■ **What just happened.** Two ideas did all the work. (i) Orthogonal invariance of the Frobenius norm lets you split any error into "inside $C(B)$" and "outside $C(B)$" pieces and throw the first away. (ii) The energy of $A$ captured by a $k$-dimensional subspace is a weighted average of the $\sigma_i^2$ with total weight $k$ and no weight exceeding 1, and such an average is maximized by the top $k$. Step 8 is a one-line linear program; the rest is bookkeeping with orthonormal bases. If $\sigma_k > \sigma_{k+1}$, the minimizer is unique (Exercise B4 looks at what happens when it is not).

**Theorem 3.8 (Eckart–Young, spectral norm).** For every $B$ with $\operatorname{rank}B \le k < r$, $\|A - B\|_2 \ge \sigma_{k+1}$, with equality for $B = A_k$.

**Proof.** This one is short enough to do in full. $N(B)$ has dimension $n - \operatorname{rank}B \ge n - k$ (Rank–Nullity). $W = \operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_{k+1}\}$ has dimension $k + 1$. For subspaces of $\mathbb{R}^n$, $\dim(N(B) \cap W) \ge \dim N(B) + \dim W - n \ge (n - k) + (k + 1) - n = 1$ (LA-3: $\dim(S + T) = \dim S + \dim T - \dim(S \cap T) \le n$). So there is a unit vector $\mathbf{w} \in N(B) \cap W$. Write $\mathbf{w} = \sum_{i \le k+1}c_i\mathbf{v}_i$ with $\sum c_i^2 = 1$. Then $(A - B)\mathbf{w} = A\mathbf{w} = \sum_{i \le k+1}c_i\sigma_i\mathbf{u}_i$, so $\|(A-B)\mathbf{w}\|^2 = \sum_{i \le k+1}c_i^2\sigma_i^2 \ge \sigma_{k+1}^2$. Hence $\|A - B\|_2 \ge \|(A-B)\mathbf{w}\| \ge \sigma_{k+1}$. For $B = A_k$, $A - A_k$ has largest singular value $\sigma_{k+1}$, so by (3.7) $\|A - A_k\|_2 = \sigma_{k+1}$. ∎

Both theorems say the same thing: to approximate $A$ with rank $k$, keep the $k$ longest axes of the ellipsoid. The Frobenius version measures the total discarded energy, the spectral version the single largest discarded stretch. Mirsky (1960) showed the same $A_k$ wins in every orthogonally invariant norm.

### 3.8 The condition number

**Definition.** For an invertible $n \times n$ matrix $A$, the (2-norm) **condition number** is

$$
\kappa_2(A) = \|A\|_2\,\|A^{-1}\|_2 = \frac{\sigma_1}{\sigma_n}.
$$

The second equality: $A^{-1} = (U\Sigma V^\top)^{-1} = V\Sigma^{-1}U^\top$, and $\Sigma^{-1} = \operatorname{diag}(1/\sigma_1, \dots, 1/\sigma_n)$. Reversing the order of the columns of $V$ and $U$ makes this a genuine SVD of $A^{-1}$ with singular values $1/\sigma_n \ge \cdots \ge 1/\sigma_1$, so $\|A^{-1}\|_2 = 1/\sigma_n$ by (3.7). Geometrically, $\kappa_2$ is the ratio of the longest to the shortest semi-axis of the ellipsoid $A(\text{unit sphere})$: $\kappa_2 = 1$ means the ellipsoid is a sphere ($A$ is a multiple of an orthogonal matrix, Exercise B7), and large $\kappa_2$ means a needle. For a singular matrix $\sigma_n = 0$ and we set $\kappa_2 = \infty$.

**Theorem 3.9 (Error amplification in $A\mathbf{x} = \mathbf{b}$).** Let $A$ be invertible, $A\mathbf{x} = \mathbf{b} \ne \mathbf{0}$, and $A(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} + \delta\mathbf{b}$. Then

$$
\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \ \le\ \kappa_2(A)\,\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}, \tag{3.10}
$$

and there are choices of $\mathbf{b}$ and $\delta\mathbf{b}$ for which equality holds.

▶ **Try it first.** Subtract the two equations to isolate $\delta\mathbf{x}$. Bound $\|\delta\mathbf{x}\|$ above using $\|A^{-1}\|_2$. Bound $\|\mathbf{x}\|$ *below* using $\|A\|_2$. Divide.

**Derivation.**

1. Subtracting, $A\,\delta\mathbf{x} = \delta\mathbf{b}$, so $\delta\mathbf{x} = A^{-1}\delta\mathbf{b}$ and by the definition of the spectral norm $\|\delta\mathbf{x}\| \le \|A^{-1}\|_2\|\delta\mathbf{b}\| = \|\delta\mathbf{b}\|/\sigma_n$.
2. $\|\mathbf{b}\| = \|A\mathbf{x}\| \le \|A\|_2\|\mathbf{x}\| = \sigma_1\|\mathbf{x}\|$, so $1/\|\mathbf{x}\| \le \sigma_1/\|\mathbf{b}\|$.
3. Multiply the two inequalities: $\|\delta\mathbf{x}\|/\|\mathbf{x}\| \le (\sigma_1/\sigma_n)\,\|\delta\mathbf{b}\|/\|\mathbf{b}\|$.
4. *Sharpness.* Take $\mathbf{b} = \mathbf{u}_1$, so $\mathbf{x} = A^{-1}\mathbf{u}_1 = \mathbf{v}_1/\sigma_1$, and $\delta\mathbf{b} = \varepsilon\mathbf{u}_n$, so $\delta\mathbf{x} = \varepsilon\mathbf{v}_n/\sigma_n$. Then $\|\delta\mathbf{x}\|/\|\mathbf{x}\| = \varepsilon\sigma_1/\sigma_n$ while $\|\delta\mathbf{b}\|/\|\mathbf{b}\| = \varepsilon$. Equality. ∎

■ **What just happened.** The worst case is a right-hand side along the long axis of the ellipse (cheap to reach: a short $\mathbf{x}$ does it) perturbed along the short axis (expensive to reach: a long $\delta\mathbf{x}$ is needed). The condition number is the ratio of the two costs. The same argument with $A$ perturbed instead of $\mathbf{b}$ gives $\|\delta\mathbf{x}\|/\|\mathbf{x} + \delta\mathbf{x}\| \le \kappa_2(A)\,\|\delta A\|_2/\|A\|_2$: from $(A + \delta A)(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b}$ subtract $A\mathbf{x} = \mathbf{b}$ to get $A\,\delta\mathbf{x} = -\delta A(\mathbf{x} + \delta\mathbf{x})$, then $\|\delta\mathbf{x}\| \le \|A^{-1}\|_2\|\delta A\|_2\|\mathbf{x} + \delta\mathbf{x}\|$ and multiply top and bottom by $\|A\|_2$.

**The rule of thumb.** Double-precision floating point represents numbers with relative error about $10^{-16}$, so the $\mathbf{b}$ and $A$ your computer actually works with are perturbed by about that much. By (3.10), the computed $\mathbf{x}$ can then be wrong by a relative error of about $\kappa_2 \cdot 10^{-16}$. Starting from 16 significant digits, **you lose about $\log_{10}\kappa_2$ of them.** With $\kappa_2 = 10^4$ you keep 12 digits; with $\kappa_2 = 10^{12}$ you keep 4; with $\kappa_2 \ge 10^{16}$ the answer is noise. This is a property of the *problem* $(A, \mathbf{b})$, not of the algorithm: no method can recover information that the rounding of the inputs destroyed. (A bad algorithm can lose more; LA-13 separates *conditioning* of the problem from *stability* of the algorithm, and shows that forming $A^\top A$ squares $\kappa$, which is why the §3.2 derivation is not how the SVD is computed.) Two more cautions. The condition number depends on the norm: $\kappa_1 = \|A\|_1\|A^{-1}\|_1$ and $\kappa_\infty$ are different numbers from $\kappa_2$ (Exercise B7 shows a matrix with $\kappa_2 = \phi^2 \approx 2.62$ and $\kappa_1 = \kappa_\infty = 4$); they agree to within factors of $n$, which is what matters for the digit count. And $\kappa$ is not the determinant: $\varepsilon I$ has determinant $\varepsilon^n$ and $\kappa_2 = 1$.

### 3.9 The pseudoinverse

If $A$ is square and invertible, $A^{-1} = V\Sigma^{-1}U^\top$: invert the stretches and swap the rotations. If $A$ is rectangular or singular, some $\sigma_i$ are zero and cannot be inverted. Invert the ones you can and leave the rest at zero.

**Definition (Moore–Penrose pseudoinverse).** Let $A = U\Sigma V^\top$ have rank $r$. Define $\Sigma^+$ as the $n \times m$ matrix with $(\Sigma^+)_{ii} = 1/\sigma_i$ for $i \le r$ and zeros elsewhere, and

$$
A^+ = V\Sigma^+U^\top = \sum_{i=1}^r\frac{1}{\sigma_i}\mathbf{v}_i\mathbf{u}_i^\top \qquad (n \times m). \tag{3.11}
$$

$A^+$ sends $\mathbf{u}_i \mapsto \mathbf{v}_i/\sigma_i$ for $i \le r$ and kills $\mathbf{u}_{r+1}, \dots, \mathbf{u}_m$ (the left null space). It does not depend on the choice of singular vectors, because it is characterized by Theorem 3.10. (Moore defined it in 1920 and Penrose, independently, in 1955 by four algebraic equations; Exercise B5 verifies them.)

**Theorem 3.10 (Minimum-norm least squares).** For every $\mathbf{b} \in \mathbb{R}^m$, the vector $\mathbf{x}^+ = A^+\mathbf{b}$ minimizes $\|A\mathbf{x} - \mathbf{b}\|$ over all $\mathbf{x} \in \mathbb{R}^n$. The set of all minimizers is $\mathbf{x}^+ + N(A)$, and among them $\mathbf{x}^+$ is the unique one of smallest norm. Moreover $\mathbf{x}^+ \in C(A^\top)$.

▶ **Try it first.** Expand $\mathbf{x}$ in the $\mathbf{v}$-basis and $\mathbf{b}$ in the $\mathbf{u}$-basis. Write $\|A\mathbf{x} - \mathbf{b}\|^2$ as a sum of squares, one per basis direction. Which coefficients of $\mathbf{x}$ affect the residual? Which affect only $\|\mathbf{x}\|$?

**Derivation.**

1. Write $\mathbf{x} = \sum_{j=1}^nc_j\mathbf{v}_j$ and $\mathbf{b} = \sum_{i=1}^md_i\mathbf{u}_i$ with $d_i = \mathbf{u}_i\cdot\mathbf{b}$. By (3.1$'$), $A\mathbf{x} = \sum_{j \le r}\sigma_jc_j\mathbf{u}_j$, so

$$
A\mathbf{x} - \mathbf{b} = \sum_{j \le r}(\sigma_jc_j - d_j)\mathbf{u}_j - \sum_{j > r}d_j\mathbf{u}_j, \qquad \|A\mathbf{x} - \mathbf{b}\|^2 = \sum_{j \le r}(\sigma_jc_j - d_j)^2 + \sum_{j > r}d_j^2 .
$$

2. The second sum does not depend on $\mathbf{x}$ at all: it is the part of $\mathbf{b}$ in $N(A^\top) = C(A)^\perp$, which no choice of $\mathbf{x}$ can reach. The first sum is a sum of squares that is minimized, with value zero, exactly when $c_j = d_j/\sigma_j$ for every $j \le r$. The coefficients $c_{r+1}, \dots, c_n$ do not appear. So the minimizers are precisely

$$
\mathbf{x} = \sum_{j \le r}\frac{d_j}{\sigma_j}\mathbf{v}_j + \sum_{j > r}c_j\mathbf{v}_j, \qquad c_{r+1}, \dots, c_n \text{ arbitrary},
$$

and the second sum ranges over exactly $N(A) = \operatorname{span}\{\mathbf{v}_{r+1}, \dots, \mathbf{v}_n\}$ (Theorem 3.3(d)). This is the set $\mathbf{x}^+ + N(A)$ with $\mathbf{x}^+ := \sum_{j \le r}(d_j/\sigma_j)\mathbf{v}_j$.

3. By orthonormality, $\|\mathbf{x}\|^2 = \sum_{j \le r}(d_j/\sigma_j)^2 + \sum_{j > r}c_j^2$, which is smallest, uniquely, when all $c_j = 0$ for $j > r$. That minimizer is $\mathbf{x}^+$, and it lies in $\operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_r\} = C(A^\top)$.

4. Finally $\mathbf{x}^+ = \sum_{j \le r}\frac{\mathbf{u}_j\cdot\mathbf{b}}{\sigma_j}\mathbf{v}_j = \Bigl(\sum_{j \le r}\frac{1}{\sigma_j}\mathbf{v}_j\mathbf{u}_j^\top\Bigr)\mathbf{b} = A^+\mathbf{b}$. ∎

■ **What just happened.** In the singular bases the least-squares problem is $n$ decoupled scalar problems $\sigma_jc_j \approx d_j$ plus $m - r$ unreachable components of $\mathbf{b}$. The pseudoinverse solves each solvable scalar problem exactly and sets every free coefficient to zero. The minimum residual is $\sqrt{\sum_{j > r}d_j^2}$, the distance from $\mathbf{b}$ to $C(A)$.

**Theorem 3.11 (Special cases).** (a) If $\operatorname{rank}A = n$ (independent columns), $A^+ = (A^\top A)^{-1}A^\top$. (b) If $\operatorname{rank}A = m$ (independent rows), $A^+ = A^\top(AA^\top)^{-1}$. (c) If $A$ is square and invertible, $A^+ = A^{-1}$. (d) In general $AA^+ = U_rU_r^\top$ is the orthogonal projector onto $C(A)$ and $A^+A = V_rV_r^\top$ is the orthogonal projector onto $C(A^\top)$.

**Proof.** (a) With $r = n$, $\Sigma$ is $m \times n$ with the invertible $\Sigma_n = \operatorname{diag}(\sigma_1, \dots, \sigma_n)$ on top of $m - n$ zero rows, and $\Sigma^\top\Sigma = \Sigma_n^2$ ($n \times n$, invertible). Then
$(A^\top A)^{-1}A^\top = (V\Sigma_n^2V^\top)^{-1}V\Sigma^\top U^\top = V\Sigma_n^{-2}V^\top V\Sigma^\top U^\top = V(\Sigma_n^{-2}\Sigma^\top)U^\top,$
and $\Sigma_n^{-2}\Sigma^\top$ is the $n \times m$ matrix $\begin{pmatrix}\Sigma_n^{-1} & 0\end{pmatrix} = \Sigma^+$. (b) is (a) applied to $A^\top$, using $(A^\top)^+ = (A^+)^\top$ (transpose (3.11)). (c) With $r = m = n$, $\Sigma^+ = \Sigma^{-1}$. (d) $AA^+ = U\Sigma V^\top V\Sigma^+U^\top = U(\Sigma\Sigma^+)U^\top$, and $\Sigma\Sigma^+$ is $m \times m$ with ones in the first $r$ diagonal positions and zeros elsewhere, so $U(\Sigma\Sigma^+)U^\top = \sum_{i \le r}\mathbf{u}_i\mathbf{u}_i^\top = U_rU_r^\top$, the projector onto $\operatorname{span}\{\mathbf{u}_1, \dots, \mathbf{u}_r\} = C(A)$. Likewise $A^+A = V(\Sigma^+\Sigma)V^\top = V_rV_r^\top$. ∎

Part (a) is the least-squares formula of LA-5, $\hat{\mathbf{x}} = (A^\top A)^{-1}A^\top\mathbf{b}$, and part (d) says $AA^+ = A(A^\top A)^{-1}A^\top$ is the LA-5 projection matrix $P$ onto the column space. The pseudoinverse is what least squares becomes when the columns are allowed to be dependent. Numerically, computing $A^+\mathbf{b}$ through the SVD (what `np.linalg.lstsq` and `np.linalg.pinv` do) is preferred over forming $A^\top A$, because $\kappa_2(A^\top A) = \kappa_2(A)^2$ (Exercise B7). In ST-3 you will see *ridge regression* replace $1/\sigma_i$ in (3.11) by $\sigma_i/(\sigma_i^2 + \lambda)$: it is a pseudoinverse that also damps the small singular values instead of inverting them.

### 3.10 Principal component analysis

**Setup.** You have $n$ observations of $p$ variables: a raw data matrix $\tilde X$ of size $n \times p$, one observation per row. Let $\boldsymbol{\mu} = \tilde X^\top\mathbf{1}/n$ be the vector of column means. **Center** the data:

$$
X = \tilde X - \mathbf{1}\boldsymbol{\mu}^\top, \qquad \text{so that every column of } X \text{ sums to zero.}
$$

The **sample covariance matrix** is

$$
S = \frac{X^\top X}{n-1} \qquad (p \times p). \tag{3.12}
$$

Its $(j,k)$ entry is $\frac{1}{n-1}\sum_{i}X_{ij}X_{ik}$, the sample covariance of variables $j$ and $k$; the diagonal entries are the sample variances. $S$ is symmetric positive semidefinite because it is a Gram matrix (step 1 of §3.2). (F5 and ST-1 discuss the $n - 1$. Nothing below depends on it.)

**Variance along a direction.** Pick a unit vector $\mathbf{w} \in \mathbb{R}^p$. The projection of observation $i$ (row $i$ of $X$) onto the line through $\mathbf{w}$ has coordinate $(X\mathbf{w})_i$; the $n$ numbers $\mathbf{t} = X\mathbf{w}$ are the **scores** along $\mathbf{w}$. Their mean is $\mathbf{1}^\top X\mathbf{w}/n = 0$ because the columns of $X$ sum to zero, so their sample variance is

$$
\frac{\|X\mathbf{w}\|^2}{n-1} = \frac{\mathbf{w}^\top X^\top X\mathbf{w}}{n-1} = \mathbf{w}^\top S\mathbf{w}.
$$

PCA asks: which direction has the most variance? Then, perpendicular to it, which next?

**Theorem 3.12 (PCA).** Let the centered $n \times p$ data matrix have SVD $X = U\Sigma V^\top$ with singular values $\sigma_1 \ge \cdots \ge \sigma_p \ge 0$ (padding with zeros if $\operatorname{rank}X < p$). Then:

- (a) The eigenvectors of $S$ are the right singular vectors $\mathbf{v}_1, \dots, \mathbf{v}_p$, with eigenvalues $\lambda_i = \sigma_i^2/(n-1)$.
- (b) Over all unit $\mathbf{w}$, the variance $\mathbf{w}^\top S\mathbf{w}$ of the scores is maximized at $\mathbf{w} = \pm\mathbf{v}_1$, with maximum value $\sigma_1^2/(n-1)$.
- (c) Over all unit $\mathbf{w}$ orthogonal to $\mathbf{v}_1, \dots, \mathbf{v}_{k-1}$, the variance is maximized at $\mathbf{w} = \pm\mathbf{v}_k$, with value $\sigma_k^2/(n-1)$.
- (d) The total variance $\operatorname{tr}S = \sum_j S_{jj}$ equals $\sum_i\sigma_i^2/(n-1) = \|X\|_F^2/(n-1)$, so the **fraction of variance explained** by the first $k$ components is $\sum_{i \le k}\sigma_i^2\big/\sum_{i \le p}\sigma_i^2$.

The $\mathbf{v}_i$ are the **principal components** (also *principal axes* or *loadings*), and the $\sigma_i\mathbf{u}_i = X\mathbf{v}_i$ are the **score vectors**.

▶ **Try it first.** Part (a) is one line from $S = X^\top X/(n-1)$ and (3.3). For (b), recognize $\mathbf{w}^\top S\mathbf{w}$ on unit vectors as the Rayleigh quotient of $S$ and quote LA-8, or redo the LA-8 argument by expanding $\mathbf{w}$ in the eigenbasis.

**Derivation.**

1. *(a).* $S = \dfrac{X^\top X}{n-1} = \dfrac{V\Sigma^\top\Sigma V^\top}{n-1} = V\,\dfrac{\Sigma^\top\Sigma}{n-1}\,V^\top$, and $\Sigma^\top\Sigma = \operatorname{diag}(\sigma_1^2, \dots, \sigma_p^2)$. This is an orthogonal diagonalization of $S$: eigenvectors $\mathbf{v}_i$, eigenvalues $\sigma_i^2/(n-1)$.

2. *(b).* The function $R(\mathbf{w}) = \mathbf{w}^\top S\mathbf{w}/\mathbf{w}^\top\mathbf{w}$ is the Rayleigh quotient of the symmetric matrix $S$. LA-8 proved: $\max_{\mathbf{w} \ne 0}R(\mathbf{w}) = \lambda_{\max}(S)$, attained exactly at the eigenvectors for $\lambda_{\max}$. Here is the argument again in this notation, since it is short. Write $\mathbf{w} = \sum_ic_i\mathbf{v}_i$ with $\sum_ic_i^2 = 1$. Then $S\mathbf{w} = \sum_i\lambda_ic_i\mathbf{v}_i$ and $\mathbf{w}^\top S\mathbf{w} = \sum_i\lambda_ic_i^2 \le \lambda_1\sum_ic_i^2 = \lambda_1$, with equality iff all the weight $c_i^2$ sits on indices with $\lambda_i = \lambda_1$. If $\lambda_1 > \lambda_2$ that means $\mathbf{w} = \pm\mathbf{v}_1$. The maximum variance is $\lambda_1 = \sigma_1^2/(n-1)$.

3. *(c).* The constraint $\mathbf{w} \perp \mathbf{v}_1, \dots, \mathbf{v}_{k-1}$ means $c_1 = \cdots = c_{k-1} = 0$, so $\mathbf{w}^\top S\mathbf{w} = \sum_{i \ge k}\lambda_ic_i^2 \le \lambda_k$, with equality at $\mathbf{w} = \pm\mathbf{v}_k$.

4. *(d).* The trace of a matrix is the sum of its eigenvalues (LA-7), so $\operatorname{tr}S = \sum_i\lambda_i = \sum_i\sigma_i^2/(n-1)$, which is $\|X\|_F^2/(n-1)$ by (3.8). ∎

■ **What just happened.** PCA is the spectral theorem applied to the covariance matrix, and the SVD of the data matrix hands you that spectral decomposition without forming the covariance matrix. The eigenvalues of $S$ are variances; the singular values of $X$ are their square roots times $\sqrt{n-1}$. Everything is the Rayleigh quotient of LA-8 wearing statistical clothes.

**The two-sided view: scores and loadings.** Write (3.6) for the data matrix:

$$
X = \sum_{i=1}^p\sigma_i\mathbf{u}_i\mathbf{v}_i^\top, \qquad XV = U\Sigma .
$$

The right singular vector $\mathbf{v}_i \in \mathbb{R}^p$ is a direction in *variable space*: its entries (the loadings) say how much each original variable contributes to component $i$. The vector $\sigma_i\mathbf{u}_i = X\mathbf{v}_i \in \mathbb{R}^n$ lives in *observation space*: its entries are the scores, the coordinate of each observation along component $i$. Row $i$ of $X$ (one observation) is $\sum_k(\sigma_ku_{ik})\mathbf{v}_k^\top$: the observation reconstructed from its scores. Keeping $k$ components, the $n \times k$ matrix of scores $U_k\Sigma_k = XV_k$ is the reduced-dimension data, and $X_k = U_k\Sigma_kV_k^\top$ is its reconstruction in the original coordinates. By Eckart–Young, $X_k$ is the best rank-$k$ approximation to the centered data in the Frobenius norm, which is to say: *the $k$-dimensional subspace through the centroid that minimizes the sum of squared perpendicular distances to the data points is spanned by $\mathbf{v}_1, \dots, \mathbf{v}_k$.* This is Pearson's 1901 formulation (lines and planes of closest fit); Hotelling's 1933 formulation is the variance-maximizing one in Theorem 3.12. Exercise B6 shows they agree. Two practical notes. PCA is not scale-invariant: measuring one variable in millimeters instead of meters multiplies its row and column of $S$ by $10^3$ and reorients every component; standardizing each column to unit variance (PCA on the correlation matrix) is a common choice, and it is a choice. And centering is not optional: §8 explains what goes wrong without it.

### 3.11 The polar decomposition and the deformation gradient

**Theorem 3.13 (Polar decomposition).** Every square $n \times n$ matrix $A$ can be written

$$
A = QS, \qquad Q \text{ orthogonal}, \quad S \text{ symmetric positive semidefinite}, \tag{3.13}
$$

with $S = \sqrt{A^\top A}$ (the unique positive semidefinite square root) and $Q = UV^\top$. If $A$ is invertible, $Q$ is unique as well. There is also a left version $A = S'Q$ with the same $Q$ and $S' = \sqrt{AA^\top} = QSQ^\top$.

▶ **Try it first.** Insert $I = V^\top V$ into the middle of $U\Sigma V^\top$ and regroup.

**Derivation.**

1. $A = U\Sigma V^\top = U(V^\top V)\Sigma V^\top = (UV^\top)(V\Sigma V^\top)$. Set $Q = UV^\top$ and $S = V\Sigma V^\top$.
2. $Q$ is a product of orthogonal matrices, so $Q^\top Q = VU^\top UV^\top = I$: orthogonal.
3. $S^\top = V\Sigma^\top V^\top = S$ (square $\Sigma$ is diagonal), and $S = V\Sigma V^\top$ is an orthogonal diagonalization with eigenvalues $\sigma_i \ge 0$: positive semidefinite. Its eigenvectors are the right singular vectors $\mathbf{v}_i$.
4. $S^2 = V\Sigma V^\top V\Sigma V^\top = V\Sigma^2V^\top = A^\top A$, so $S$ is a PSD square root of $A^\top A$. LA-8 showed a PSD matrix has exactly one PSD square root; so $S$ is determined by $A$ even though $U$ and $V$ are not.
5. If $A$ is invertible then so is $S$ (all $\sigma_i > 0$) and $Q = AS^{-1}$ is determined.
6. Left version: insert $U^\top U$ instead: $A = (U\Sigma U^\top)(UV^\top) = S'Q$ with $S' = U\Sigma U^\top = (UV^\top)(V\Sigma V^\top)(VU^\top) = QSQ^\top$, and $S'^2 = AA^\top$. ∎

■ **What just happened.** "Rotate, stretch, rotate" regrouped as "stretch (along the $\mathbf{v}_i$ axes, in place), then rotate" or as "rotate, then stretch (along the $\mathbf{u}_i$ axes)." The stretch $S$ is the symmetric part of $A$ in the multiplicative sense, the way $\sqrt{z\bar z}$ is the modulus of a complex number and $z/|z|$ its phase: $A = QS$ is the matrix version of $z = e^{i\theta}|z|$. Note $\det Q = \pm 1$ with the sign of $\det A$, since $\det A = \det Q\det S$ and $\det S = \prod\sigma_i \ge 0$. Exercise C1 shows $Q$ is also the *nearest* orthogonal matrix to $A$.

**The deformation gradient.** This is where the ME student should sit up. In continuum mechanics a deformation carries each material point $\mathbf{X}$ of a body to a spatial point $\mathbf{x}(\mathbf{X})$, and the **deformation gradient** is the Jacobian $F = \partial\mathbf{x}/\partial\mathbf{X}$ ($3 \times 3$, or $2 \times 2$ in plane problems). A short material fiber $d\mathbf{X}$ becomes $d\mathbf{x} = F\,d\mathbf{X}$. Physically, $F$ does two distinct things: it changes the lengths of fibers and the angles between them (strain), and it rotates the neighborhood rigidly (which produces no stress). The polar decomposition separates them. Mechanics texts write

$$
F = RU = VR \qquad \text{(mechanics notation)},
$$

with $R$ the **rotation tensor**, $U = \sqrt{F^\top F}$ the **right stretch tensor** and $V = \sqrt{FF^\top}$ the **left stretch tensor**. In our notation $R = Q = U_{\text{svd}}V_{\text{svd}}^\top$, $U_{\text{mech}} = S = V_{\text{svd}}\Sigma V_{\text{svd}}^\top$, and $V_{\text{mech}} = S' = U_{\text{svd}}\Sigma U_{\text{svd}}^\top$. Read that twice: *the mechanics $U$ is our $S$, and the mechanics $V$ is our $S'$; neither is our $U$ or $V$.* With that translation, the SVD of $F$ says exactly:

- The **principal stretches** $\lambda_1 \ge \lambda_2 \ge \lambda_3 > 0$ of the deformation are the singular values of $F$. (Here $\det F > 0$ for a physical deformation, so $F$ is invertible and $\det R = +1$: a proper rotation.)
- The **principal material directions** (eigenvectors of $U_{\text{mech}}$, the fibers that are stretched without being sheared relative to each other) are the right singular vectors $\mathbf{v}_i$ of $F$.
- The **principal spatial directions** (eigenvectors of $V_{\text{mech}}$, where those fibers end up) are the left singular vectors $\mathbf{u}_i = R\mathbf{v}_i$.
- The **right Cauchy–Green tensor** $C = F^\top F = U_{\text{mech}}^2$ has eigenvalues $\lambda_i^2 = \sigma_i^2$, and the **Green–Lagrange strain** $E = \tfrac12(C - I)$ has eigenvalues $\tfrac12(\sigma_i^2 - 1)$, which vanish exactly when $F$ is a pure rotation (all $\sigma_i = 1$).

So "$F = RU$: a deformation is a stretch followed by a rotation" is a theorem of this module, and the recipe for computing $R$ and $U$ by hand is: diagonalize $C = F^\top F$. Example 4.5 does it for simple shear.

### 3.12 Relation to the eigendecomposition

**Symmetric positive (semi)definite.** If $A = Q\Lambda Q^\top$ with all $\lambda_i \ge 0$ (LA-8), this is already an SVD: $U = V = Q$ and $\sigma_i = \lambda_i$ (after sorting). Covariance matrices, Gram matrices, inertia tensors, stiffness matrices: for all of these "SVD" and "eigendecomposition" are the same computation.

**Symmetric indefinite.** If $A = Q\Lambda Q^\top$ with some $\lambda_i < 0$, write $\lambda_i = s_i|\lambda_i|$ with $s_i = \pm 1$ and let $D = \operatorname{diag}(s_i)$. Then $A = (QD)\,|\Lambda|\,Q^\top$ with $QD$ orthogonal. So $\sigma_i = |\lambda_i|$, $V = Q$, and $U = QD$: the left singular vector for a negative eigenvalue is *minus* the eigenvector. The signs are absorbed into $U$, and the order of the singular values (by $|\lambda|$) may differ from the order of the eigenvalues. Example: $A = \operatorname{diag}(3, -2)$ has $\Sigma = \operatorname{diag}(3, 2)$, $V = I$, $U = \operatorname{diag}(1, -1)$ (Exercise A1).

**Nonsymmetric square.** No simple relation. The matrix $\begin{pmatrix}0 & 1\\0 & 0\end{pmatrix}$ has both eigenvalues zero and singular values $1, 0$. What survives: $\prod_i\sigma_i = |\det A| = \prod_i|\lambda_i|$ (take determinants in (3.1)), and $\sigma_n \le |\lambda| \le \sigma_1$ for every eigenvalue $\lambda$ (if $A\mathbf{x} = \lambda\mathbf{x}$ with $\|\mathbf{x}\| = 1$ then $|\lambda| = \|A\mathbf{x}\|$, and (3.7) bounds this). For *normal* matrices (LA-8: $A^\top A = AA^\top$, which includes orthogonal and skew-symmetric matrices) one has $\sigma_i = |\lambda_i|$ again. For everything else, eigenvalues describe what $A$ does when iterated (LA-7: powers, $e^{tA}$, stability) and singular values describe what $A$ does *once* (stretch, approximation, conditioning). Different questions, different numbers.

### 3.13 Why low-rank structure is everywhere

The SVD would be a curiosity if real matrices had slowly decaying singular values. They usually do not, and it is worth knowing why.

- **Images.** A photograph as a matrix of pixel intensities has strongly correlated rows: neighboring rows look alike, large regions are smooth or textured in a repeating way. Correlated rows means nearly dependent rows means small trailing singular values. Fine detail and noise live in the small ones. The project compresses a $256 \times 256$ image at ranks 1, 5, 20, 50.
- **Covariance matrices.** If $p$ observed variables are all driven by $k \ll p$ latent factors plus independent noise, then the population covariance is (rank $k$) + (diagonal), and the sample covariance has $k$ large eigenvalues and a floor of small ones. This was Spearman's 1904 observation about test scores and the reason Eckart and Young wrote their paper. The scree plot (eigenvalues in decreasing order) often shows an elbow at $k$.
- **Co-occurrence and ratings.** A term–document matrix (how often word $i$ appears in document $j$) or a user–item matrix (how user $i$ rated item $j$) is well approximated by a rank-$k$ matrix when a few topics or tastes explain most of the data. Latent semantic analysis (Deerwester et al., 1990) is the SVD of the term–document matrix; ML-6 and ML-9 build on this view.
- **Smooth functions.** If $a_{ij} = f(x_i, y_j)$ for a smooth kernel $f$, the matrix is close to low rank because a smooth function of two variables is well approximated by a short sum $\sum_k g_k(x)h_k(y)$ (think Taylor or Fourier). This is Schmidt's 1907 setting and the reason integral operators, Green's functions (EM-3), and discretized PDE solution operators have fast-decaying singular values.
- **Linear physical systems.** A response matrix (inputs to outputs) of a system with a few dominant modes has a few dominant singular values; CM-6's normal modes are the eigenvector version of the same statement.

What is *not* low rank: the identity, random matrices (all singular values comparable; the project's "noise floor" shows this), and any matrix whose rows carry genuinely independent information. Adding noise to a low-rank matrix produces a spectrum with $k$ large singular values and a flat floor of roughly equal small ones. Eckart–Young then says: truncating at $k$ removes most of the noise and little of the signal. That is why the SVD is, in practice, a denoising tool as much as a compression tool. (The precise shape of the noise floor is random-matrix theory; F5 touches the Gaussian case.)

## 4. Worked examples

### Example 4.1: The SVD of a $2 \times 2$ matrix by hand

$$
A = \begin{pmatrix} 3 & 0 \\ 4 & 5 \end{pmatrix}.
$$

*Step 1: $A^\top A$ and its eigenvalues.*
$A^\top A = \begin{pmatrix}3 & 4\\0 & 5\end{pmatrix}\begin{pmatrix}3 & 0\\4 & 5\end{pmatrix} = \begin{pmatrix}9 + 16 & 0 + 20\\0 + 20 & 0 + 25\end{pmatrix} = \begin{pmatrix}25 & 20\\20 & 25\end{pmatrix}.$
Trace $50$, determinant $625 - 400 = 225$. Eigenvalues solve $\lambda^2 - 50\lambda + 225 = 0$: $\lambda = 25 \pm \sqrt{625 - 225} = 25 \pm 20$, so $\lambda_1 = 45$, $\lambda_2 = 5$. (Check: $\det(A^\top A) = (\det A)^2 = 15^2 = 225 = 45 \cdot 5$. ✓)

*Step 2: singular values.* $\sigma_1 = \sqrt{45} = 3\sqrt5 \approx 6.708$, $\sigma_2 = \sqrt5 \approx 2.236$. Checks: $\sigma_1\sigma_2 = \sqrt{225} = 15 = |\det A|$ ✓; $\sigma_1^2 + \sigma_2^2 = 50 = 9 + 0 + 16 + 25 = \|A\|_F^2$ ✓ (Theorem 3.6).

*Step 3: right singular vectors.* For $\lambda_1 = 45$: $(A^\top A - 45I)\mathbf{v} = \begin{pmatrix}-20 & 20\\20 & -20\end{pmatrix}\mathbf{v} = \mathbf{0}$ gives $v_1 = v_2$, so $\mathbf{v}_1 = \tfrac{1}{\sqrt2}(1, 1)^\top$. For $\lambda_2 = 5$: $\begin{pmatrix}20 & 20\\20 & 20\end{pmatrix}\mathbf{v} = \mathbf{0}$ gives $\mathbf{v}_2 = \tfrac{1}{\sqrt2}(1, -1)^\top$. Orthogonal, as the spectral theorem promised.

*Step 4: left singular vectors from (3.4).*
$A\mathbf{v}_1 = \tfrac{1}{\sqrt2}\begin{pmatrix}3\\4 + 5\end{pmatrix} = \tfrac{1}{\sqrt2}\begin{pmatrix}3\\9\end{pmatrix}$, and $\|A\mathbf{v}_1\| = \tfrac{1}{\sqrt2}\sqrt{90} = \sqrt{45} = \sigma_1$ ✓. So $\mathbf{u}_1 = \dfrac{A\mathbf{v}_1}{\sigma_1} = \dfrac{(3, 9)^\top}{\sqrt2\cdot3\sqrt5} = \dfrac{(1, 3)^\top}{\sqrt{10}}$.
$A\mathbf{v}_2 = \tfrac{1}{\sqrt2}\begin{pmatrix}3\\4 - 5\end{pmatrix} = \tfrac{1}{\sqrt2}\begin{pmatrix}3\\-1\end{pmatrix}$, $\|A\mathbf{v}_2\| = \sqrt{10/2} = \sqrt5 = \sigma_2$ ✓. So $\mathbf{u}_2 = \dfrac{(3, -1)^\top}{\sqrt2\cdot\sqrt5} = \dfrac{(3, -1)^\top}{\sqrt{10}}$.
Check orthogonality (the step 3 miracle of §3.2): $\mathbf{u}_1\cdot\mathbf{u}_2 = (3 - 3)/10 = 0$ ✓.

*Step 5: assemble.*
$$
A = U\Sigma V^\top = \frac{1}{\sqrt{10}}\begin{pmatrix}1 & 3\\3 & -1\end{pmatrix}\begin{pmatrix}3\sqrt5 & 0\\0 & \sqrt5\end{pmatrix}\frac{1}{\sqrt2}\begin{pmatrix}1 & 1\\1 & -1\end{pmatrix}.
$$
Multiply back: $U\Sigma = \tfrac{\sqrt5}{\sqrt{10}}\begin{pmatrix}3 & 3\\9 & -1\end{pmatrix} = \tfrac{1}{\sqrt2}\begin{pmatrix}3 & 3\\9 & -1\end{pmatrix}$, then $(U\Sigma)V^\top = \tfrac12\begin{pmatrix}3 + 3 & 3 - 3\\9 - 1 & 9 + 1\end{pmatrix} = \begin{pmatrix}3 & 0\\4 & 5\end{pmatrix}$ ✓.

*Geometry.* The unit circle maps to an ellipse with semi-axes of length $6.708$ along $(1, 3)/\sqrt{10}$ (about $71.6°$ from the $x$-axis) and $2.236$ along $(3, -1)/\sqrt{10}$. The input directions that land on the axes are the diagonals $(1, 1)/\sqrt2$ and $(1, -1)/\sqrt2$. $\|A\|_2 = 6.708$, $\|A\|_F = \sqrt{50} = 7.071$, $\kappa_2 = \sqrt{45}/\sqrt5 = 3$. The rank-one approximation is $A_1 = \sigma_1\mathbf{u}_1\mathbf{v}_1^\top = \tfrac{3\sqrt5}{\sqrt{20}}\begin{pmatrix}1 & 1\\3 & 3\end{pmatrix} = \begin{pmatrix}1.5 & 1.5\\4.5 & 4.5\end{pmatrix}$, and $\|A - A_1\|_F^2 = 1.5^2 + 1.5^2 + 0.5^2 + 0.5^2 = 5 = \sigma_2^2$, as Eckart–Young says.

*What NumPy returns.* `np.linalg.svd(A)` gives `U = [[-0.316, -0.949], [-0.949, 0.316]]` and `Vt = [[-0.707, -0.707], [-0.707, 0.707]]`: it has chosen $-\mathbf{u}_1, -\mathbf{v}_1$ and $-\mathbf{u}_2, -\mathbf{v}_2$. Both pairs are flipped together, so $U\Sigma V^\top$ is still $A$. Neither answer is wrong (§3.3).

### Example 4.2: A rank-one $3 \times 2$ matrix and its four subspaces

$$
A = \begin{pmatrix}1 & 2\\2 & 4\\3 & 6\end{pmatrix} = \begin{pmatrix}1\\2\\3\end{pmatrix}\begin{pmatrix}1 & 2\end{pmatrix}.
$$

Rank 1 by inspection (column 2 is twice column 1). So $r = 1$ and we expect one nonzero singular value.

*$A^\top A$.* $A^\top A = \begin{pmatrix}1 + 4 + 9 & 2 + 8 + 18\\2 + 8 + 18 & 4 + 16 + 36\end{pmatrix} = \begin{pmatrix}14 & 28\\28 & 56\end{pmatrix} = 14\begin{pmatrix}1 & 2\\2 & 4\end{pmatrix}$. Eigenvalues: trace $70$, determinant $0$, so $\lambda_1 = 70$, $\lambda_2 = 0$. $\sigma_1 = \sqrt{70} \approx 8.367$, $\sigma_2 = 0$. Check: $\|A\|_F^2 = 1 + 4 + 9 + 4 + 16 + 36 = 70$ ✓.

*Right singular vectors.* For $\lambda_1 = 70$: $\begin{pmatrix}-56 & 28\\28 & -14\end{pmatrix}\mathbf{v} = \mathbf{0}$ gives $v_2 = 2v_1$: $\mathbf{v}_1 = (1, 2)^\top/\sqrt5$. For $\lambda_2 = 0$: $\mathbf{v}_2 \perp \mathbf{v}_1$: $\mathbf{v}_2 = (2, -1)^\top/\sqrt5$. Check $A\mathbf{v}_2 = \tfrac{1}{\sqrt5}(2 - 2, 4 - 4, 6 - 6)^\top = \mathbf{0}$ ✓.

*Left singular vector.* $A\mathbf{v}_1 = \tfrac{1}{\sqrt5}(1 + 4, 2 + 8, 3 + 12)^\top = \tfrac{5}{\sqrt5}(1, 2, 3)^\top = \sqrt5\,(1, 2, 3)^\top$, with norm $\sqrt5\sqrt{14} = \sqrt{70} = \sigma_1$ ✓. So $\mathbf{u}_1 = (1, 2, 3)^\top/\sqrt{14}$.

*Complete $U$.* We need $\mathbf{u}_2, \mathbf{u}_3$ orthonormal and orthogonal to $(1, 2, 3)$. Take $\mathbf{u}_2 = (2, -1, 0)^\top/\sqrt5$ (dot with $(1,2,3)$: $2 - 2 = 0$ ✓). For $\mathbf{u}_3$, run Gram–Schmidt on $\mathbf{e}_3$ against $\mathbf{u}_1, \mathbf{u}_2$, or (in $\mathbb{R}^3$ only) take the cross product $(1, 2, 3) \times (2, -1, 0) = (2\cdot 0 - 3\cdot(-1),\ 3\cdot 2 - 1\cdot 0,\ 1\cdot(-1) - 2\cdot 2) = (3, 6, -5)$, norm $\sqrt{9 + 36 + 25} = \sqrt{70}$. So $\mathbf{u}_3 = (3, 6, -5)^\top/\sqrt{70}$.

*Full and reduced forms.*
$$
A = \underbrace{\begin{pmatrix}\frac{1}{\sqrt{14}} & \frac{2}{\sqrt5} & \frac{3}{\sqrt{70}}\\ \frac{2}{\sqrt{14}} & \frac{-1}{\sqrt5} & \frac{6}{\sqrt{70}}\\ \frac{3}{\sqrt{14}} & 0 & \frac{-5}{\sqrt{70}}\end{pmatrix}}_{U\ (3 \times 3)}\underbrace{\begin{pmatrix}\sqrt{70} & 0\\0 & 0\\0 & 0\end{pmatrix}}_{\Sigma\ (3 \times 2)}\underbrace{\frac{1}{\sqrt5}\begin{pmatrix}1 & 2\\2 & -1\end{pmatrix}}_{V^\top\ (2 \times 2)}
= \underbrace{\frac{1}{\sqrt{14}}\begin{pmatrix}1\\2\\3\end{pmatrix}}_{U_1}\underbrace{\sqrt{70}}_{\Sigma_1}\underbrace{\frac{1}{\sqrt5}\begin{pmatrix}1 & 2\end{pmatrix}}_{V_1^\top}.
$$
The reduced form multiplies out to $\tfrac{\sqrt{70}}{\sqrt{70}}(1, 2, 3)^\top(1, 2) = A$ ✓.

*The four subspaces, read off (Theorem 3.3).* $C(A) = \operatorname{span}\{\mathbf{u}_1\}$: the line through $(1, 2, 3)$ in $\mathbb{R}^3$. $N(A^\top) = \operatorname{span}\{\mathbf{u}_2, \mathbf{u}_3\}$: the plane $x + 2y + 3z = 0$ (both $\mathbf{u}_2, \mathbf{u}_3$ satisfy it, and so does any combination). $C(A^\top) = \operatorname{span}\{\mathbf{v}_1\}$: the line through $(1, 2)$ in $\mathbb{R}^2$. $N(A) = \operatorname{span}\{\mathbf{v}_2\}$: the line through $(2, -1)$. Dimensions $1 + 2 = 3 = m$ and $1 + 1 = 2 = n$ ✓. The pseudoinverse is $A^+ = \tfrac{1}{\sigma_1}\mathbf{v}_1\mathbf{u}_1^\top = \tfrac{1}{\sqrt{70}}\cdot\tfrac{1}{\sqrt5\sqrt{14}}\begin{pmatrix}1\\2\end{pmatrix}\begin{pmatrix}1 & 2 & 3\end{pmatrix} = \tfrac{1}{70}\begin{pmatrix}1 & 2 & 3\\2 & 4 & 6\end{pmatrix} = A^\top/70$; for a rank-one matrix, $A^+ = A^\top/\|A\|_F^2$ always (Exercise A5).

### Example 4.3: PCA by hand on four points in the plane

Data (four observations of two variables), one per row:
$$
\tilde X = \begin{pmatrix}5 & 4\\4 & 1\\1 & 0\\2 & 3\end{pmatrix}.
$$

*Center.* Column means: $\mu_1 = (5 + 4 + 1 + 2)/4 = 3$, $\mu_2 = (4 + 1 + 0 + 3)/4 = 2$. Subtract:
$$
X = \begin{pmatrix}2 & 2\\1 & -1\\-2 & -2\\-1 & 1\end{pmatrix}.
$$
Each column sums to zero ✓.

*Covariance.* $X^\top X = \begin{pmatrix}4 + 1 + 4 + 1 & 4 - 1 + 4 - 1\\ \cdot & 4 + 1 + 4 + 1\end{pmatrix} = \begin{pmatrix}10 & 6\\6 & 10\end{pmatrix}$, so $S = X^\top X/3 = \begin{pmatrix}10/3 & 2\\2 & 10/3\end{pmatrix}$. Each variable has sample variance $10/3$ and their covariance is $2$ (correlation $0.6$).

*Eigenvectors.* $X^\top X$ has trace $20$ and determinant $64$: eigenvalues $10 \pm 6 = 16, 4$. So $\sigma_1 = 4$, $\sigma_2 = 2$, and the eigenvalues of $S$ are $\lambda_1 = 16/3 \approx 5.333$, $\lambda_2 = 4/3 \approx 1.333$. Eigenvectors: for $16$, $(10 - 16)v_1 + 6v_2 = 0$ gives $\mathbf{v}_1 = (1, 1)^\top/\sqrt2$; for $4$, $\mathbf{v}_2 = (1, -1)^\top/\sqrt2$. The first principal component is the diagonal direction: the two variables tend to be high together.

*Variance explained.* Total variance $\operatorname{tr}S = 20/3$; component 1 explains $\lambda_1/\operatorname{tr}S = 16/20 = 80\%$, component 2 the remaining $20\%$. Check Theorem 3.12(b): the variance along $\mathbf{v}_1$ is $\mathbf{v}_1^\top S\mathbf{v}_1 = \tfrac12(10/3 + 2 + 2 + 10/3) = 16/3$ ✓, and along the original $x$-axis it is only $S_{11} = 10/3$.

*Scores.* $X\mathbf{v}_1 = \tfrac{1}{\sqrt2}(2 + 2,\ 1 - 1,\ -2 - 2,\ -1 + 1)^\top = (2\sqrt2, 0, -2\sqrt2, 0)^\top$ and $X\mathbf{v}_2 = \tfrac{1}{\sqrt2}(0, 2, 0, -2)^\top = (0, \sqrt2, 0, -\sqrt2)^\top$. Sum of squares: $8 + 8 = 16 = \sigma_1^2$ and $2 + 2 = 4 = \sigma_2^2$ ✓, so the score vectors are $\sigma_1\mathbf{u}_1$ and $\sigma_2\mathbf{u}_2$ with $\mathbf{u}_1 = (1, 0, -1, 0)^\top/\sqrt2$ and $\mathbf{u}_2 = (0, 1, 0, -1)^\top/\sqrt2$. Observations 1 and 3 are extreme along the first component; 2 and 4 along the second.

*Rank-one reconstruction.* $X_1 = \sigma_1\mathbf{u}_1\mathbf{v}_1^\top = 4\cdot\tfrac{1}{\sqrt2}(1, 0, -1, 0)^\top\cdot\tfrac{1}{\sqrt2}(1, 1) = \begin{pmatrix}2 & 2\\0 & 0\\-2 & -2\\0 & 0\end{pmatrix}$: observations 1 and 3 are reproduced exactly (they lie on the first axis) and 2 and 4 collapse to the centroid. Error $\|X - X_1\|_F^2 = 1 + 1 + 1 + 1 = 4 = \sigma_2^2$ ✓. In original coordinates add back $(3, 2)$ to each row.

*What centering did.* The uncentered $\tilde X^\top\tilde X = \begin{pmatrix}46 & 30\\30 & 26\end{pmatrix}$ has eigenvalues $36 \pm \sqrt{1000} \approx 67.6, 4.4$ and top eigenvector at about $36°$ from the $x$-axis: close to the direction of the mean $(3, 2)$ (at $33.7°$), not the $45°$ of the centered cloud. Without centering, the first "component" mostly reports where the data sits, not how it varies. §8.

### Example 4.4: Condition number and error amplification

$$
A = \begin{pmatrix}1 & 1\\1 & 1.0001\end{pmatrix}, \qquad \mathbf{b} = \begin{pmatrix}2\\2.0001\end{pmatrix}, \qquad \text{solution } \mathbf{x} = \begin{pmatrix}1\\1\end{pmatrix}.
$$

*Singular values.* $A$ is symmetric positive definite (trace $2.0001 > 0$, det $= 1.0001 - 1 = 10^{-4} > 0$), so by §3.12 its singular values are its eigenvalues: $\lambda = \tfrac12\bigl(2.0001 \pm \sqrt{2.0001^2 - 4\cdot10^{-4}}\bigr) = \tfrac12\bigl(2.0001 \pm \sqrt{4.00000001}\bigr) = \tfrac12(2.0001 \pm 2.0000000025)$. So $\sigma_1 = 2.00005000125$, $\sigma_2 = 0.00004999875 \approx 5.0\times10^{-5}$, and

$$
\kappa_2(A) = \frac{\sigma_1}{\sigma_2} \approx 40\,002.
$$

(Quick estimate: $\sigma_1 \approx 2$, $\sigma_1\sigma_2 = \det A = 10^{-4}$, so $\sigma_2 \approx 5\times10^{-5}$ and $\kappa \approx 4\times10^4$.) Rule of thumb: expect to lose $\log_{10}(4\times10^4) \approx 4.6$ digits.

*A concrete perturbation.* Change $\mathbf{b}$ by $\delta\mathbf{b} = (10^{-4}, -10^{-4})^\top$: a relative change of $\|\delta\mathbf{b}\|/\|\mathbf{b}\| = 1.414\times10^{-4}/2.8285 = 5.0\times10^{-5}$, i.e. $0.005\%$. The inverse is $A^{-1} = \dfrac{1}{10^{-4}}\begin{pmatrix}1.0001 & -1\\-1 & 1\end{pmatrix} = \begin{pmatrix}10001 & -10000\\-10000 & 10000\end{pmatrix}$, so

$$
\delta\mathbf{x} = A^{-1}\delta\mathbf{b} = 10^{-4}\begin{pmatrix}10001 + 10000\\-10000 - 10000\end{pmatrix} = \begin{pmatrix}2.0001\\-2\end{pmatrix}.
$$

The solution moves from $(1, 1)$ to $(3.0001, -1)$. Check: $A(3.0001, -1)^\top = (3.0001 - 1,\ 3.0001 - 1.0001)^\top = (2.0001, 2.0000)^\top = \mathbf{b} + \delta\mathbf{b}$ ✓. Relative change in $\mathbf{x}$: $\|\delta\mathbf{x}\|/\|\mathbf{x}\| = 2.8285/1.4142 = 2.00005$, i.e. $200\%$. Amplification factor: $2.00005/(5.0\times10^{-5}) = 40\,002 = \kappa_2(A)$, to the digits shown. The bound (3.10) is attained because $\delta\mathbf{b} \propto (1, -1)$ is (to four decimals) the direction $\mathbf{u}_2$ of the short axis and $\mathbf{b} \propto (1, 1)$ is (to four decimals) the long axis $\mathbf{u}_1$: this is step 4 of Theorem 3.9.

*Geometry.* The ellipse $A(\text{unit circle})$ has semi-axes $2.00005$ and $0.00005$: a needle $40\,000$ times longer than it is wide, pointing along $(1, 1)$. Any right-hand side with a component across the needle requires an enormous $\mathbf{x}$ to produce.

### Example 4.5: Polar decomposition of a simple shear (a mechanics example)

The deformation gradient of simple shear with shear parameter $\gamma = 1$ is

$$
F = \begin{pmatrix}1 & 1\\0 & 1\end{pmatrix}: \qquad x_1 = X_1 + X_2, \quad x_2 = X_2.
$$

Horizontal fibers are unchanged, vertical fibers tilt by $45°$. We want $F = RU_{\text{mech}}$: how much is rigid rotation, how much is stretch, and along which fibers? ($\det F = 1$: area is preserved, as it must be for shear.)

*Right Cauchy–Green tensor.* $C = F^\top F = \begin{pmatrix}1 & 0\\1 & 1\end{pmatrix}\begin{pmatrix}1 & 1\\0 & 1\end{pmatrix} = \begin{pmatrix}1 & 1\\1 & 2\end{pmatrix}$. Trace $3$, determinant $1$: eigenvalues $\lambda = \tfrac{3 \pm \sqrt5}{2}$. With the golden ratio $\phi = \tfrac{1 + \sqrt5}{2} \approx 1.618$, note $\phi^2 = \phi + 1 = \tfrac{3 + \sqrt5}{2}$ and $\phi^{-2} = \tfrac{3 - \sqrt5}{2}$. So the **principal stretches** are

$$
\sigma_1 = \phi \approx 1.618, \qquad \sigma_2 = 1/\phi \approx 0.618, \qquad \sigma_1\sigma_2 = 1 = \det F \ ✓.
$$

*Principal material directions (right singular vectors).* For $\lambda_1 = \phi^2$: $(1 - \phi^2)v_1 + v_2 = 0$, so $v_2 = (\phi^2 - 1)v_1 = \phi v_1$ (using $\phi^2 - 1 = \phi$). $\mathbf{v}_1 = (1, \phi)^\top/\sqrt{1 + \phi^2}$, at angle $\arctan\phi \approx 58.28°$ from the $X_1$ axis. For $\lambda_2 = \phi^{-2}$: $\mathbf{v}_2 \perp \mathbf{v}_1$, $\mathbf{v}_2 = (\phi, -1)^\top/\sqrt{1 + \phi^2}$, at $-31.72°$. Note $1 + \phi^2 = \phi + 2 \approx 3.618$. The fiber initially at $58.28°$ is stretched by $\phi$ and the fiber at $-31.72°$ is compressed by $1/\phi$; those two fibers stay perpendicular to each other during the deformation. Every other fiber is sheared relative to its neighbors.

*Where they go (left singular vectors).* $\mathbf{u}_1 = F\mathbf{v}_1/\sigma_1 = \dfrac{(1 + \phi, \phi)^\top}{\phi\sqrt{\phi + 2}} = \dfrac{(\phi, 1)^\top}{\sqrt{\phi + 2}}$ using $(1 + \phi)/\phi = \phi^2/\phi = \phi$; angle $\arctan(1/\phi) \approx 31.72°$. $\mathbf{u}_2 = F\mathbf{v}_2/\sigma_2 = \phi\,\dfrac{(\phi - 1, -1)^\top}{\sqrt{\phi + 2}} = \dfrac{(\phi^2 - \phi, -\phi)^\top}{\sqrt{\phi + 2}} = \dfrac{(1, -\phi)^\top}{\sqrt{\phi + 2}}$. Check $\mathbf{u}_1\cdot\mathbf{u}_2 = (\phi - \phi)/(\phi + 2) = 0$ ✓. The stretched fiber rotated from $58.28°$ to $31.72°$: clockwise by $26.57°$.

*Rotation.* $R = Q = \mathbf{u}_1\mathbf{v}_1^\top + \mathbf{u}_2\mathbf{v}_2^\top = \dfrac{1}{\phi + 2}\left[\begin{pmatrix}\phi\\1\end{pmatrix}\begin{pmatrix}1 & \phi\end{pmatrix} + \begin{pmatrix}1\\-\phi\end{pmatrix}\begin{pmatrix}\phi & -1\end{pmatrix}\right] = \dfrac{1}{\phi + 2}\begin{pmatrix}2\phi & \phi^2 - 1\\1 - \phi^2 & 2\phi\end{pmatrix} = \dfrac{\phi}{\phi + 2}\begin{pmatrix}2 & 1\\-1 & 2\end{pmatrix}.$
Since $\sqrt5\,\phi = \tfrac{\sqrt5 + 5}{2} = \phi + 2$, the prefactor is $1/\sqrt5$:

$$
R = \frac{1}{\sqrt5}\begin{pmatrix}2 & 1\\-1 & 2\end{pmatrix} = \begin{pmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{pmatrix}, \qquad \tan\theta = -\tfrac12, \quad \theta = -26.57°.
$$

A clockwise rotation by $\arctan(1/2) = 26.57°$. (The general result for simple shear is $\tan\theta = -\gamma/2$.)

*Stretch.* $U_{\text{mech}} = S = R^\top F = \dfrac{1}{\sqrt5}\begin{pmatrix}2 & -1\\1 & 2\end{pmatrix}\begin{pmatrix}1 & 1\\0 & 1\end{pmatrix} = \dfrac{1}{\sqrt5}\begin{pmatrix}2 & 1\\1 & 3\end{pmatrix}$. Symmetric ✓; trace $5/\sqrt5 = \sqrt5 = \phi + 1/\phi$ ✓; determinant $(6 - 1)/5 = 1$ ✓. Its eigenvectors are $\mathbf{v}_1, \mathbf{v}_2$ with eigenvalues $\phi, 1/\phi$. And $U_{\text{mech}}^2 = \tfrac15\begin{pmatrix}5 & 5\\5 & 10\end{pmatrix} = C$ ✓. The left stretch is $V_{\text{mech}} = RU_{\text{mech}}R^\top = \tfrac{1}{\sqrt5}\begin{pmatrix}3 & 1\\1 & 2\end{pmatrix}$, whose square is $FF^\top = \begin{pmatrix}2 & 1\\1 & 1\end{pmatrix}$ ✓.

*Reading.* Simple shear by $\gamma = 1$ equals: stretch by $1.618$ along the fiber at $58.3°$ and by $0.618$ along the fiber at $-31.7°$ (no rotation of those fibers yet), then rotate the whole neighborhood clockwise by $26.6°$. The Green–Lagrange strain $E = \tfrac12(C - I) = \tfrac12\begin{pmatrix}0 & 1\\1 & 1\end{pmatrix}$ has eigenvalues $\tfrac12(\phi^2 - 1) = \phi/2 \approx 0.809$ and $\tfrac12(\phi^{-2} - 1) = -\phi^{-1}/2 \approx -0.309$ in those same two directions.

## 5. Why it matters and where it connects

- **LA-10 (decompositions atlas).** The SVD is the most expensive and most informative factorization in the atlas; LA-10 compares it with QR and the eigendecomposition and says when each is the right tool. The polar decomposition (3.13) appears there as "SVD regrouped."
- **LA-13 (numerical linear algebra).** Theorem 3.9 is the *conditioning* half of LA-13's conditioning-versus-stability distinction. LA-13 derives the Golub–Kahan bidiagonalization that computes the SVD without forming $A^\top A$, explains why `svd_via_eig` in the project loses the small singular values, and uses $\sigma_n$ as the distance to the nearest singular matrix (Exercise C3). Randomized SVD (LA-13) is Eckart–Young plus a random projection.
- **LA-14 (matrix calculus).** The gradient of $\|A - B\|_F^2$, of $\operatorname{tr}(Q^\top A)$ over orthogonal $Q$ (Procrustes, Exercise C1), and of $\mathbf{w}^\top S\mathbf{w}$ on the sphere are all computed there; their stationary points are the singular vectors found here.
- **ST-3 (linear regression).** Least squares is $\hat{\boldsymbol{\beta}} = X^+\mathbf{y}$ (Theorem 3.10); *multicollinearity* means a small $\sigma_p$, so $\kappa_2(X)$ is large and $\hat{\boldsymbol{\beta}}$ is unstable (Theorem 3.9); *ridge regression* replaces $1/\sigma_i$ by $\sigma_i/(\sigma_i^2 + \lambda)$ in (3.11); the "hat matrix" is $XX^+ = U_rU_r^\top$ (Theorem 3.11(d)). ST-3 lists LA-9 as a prerequisite for exactly these four facts.
- **ML-6 (unsupervised learning).** PCA via SVD is Theorem 3.12 verbatim. Spectral clustering uses the eigenvectors of a graph Laplacian, a symmetric PSD matrix where SVD and eigendecomposition coincide (§3.12). Manifold methods start from the failure of PCA on curved data.
- **ML-9 (transformers).** Attention weight matrices and embedding matrices are analyzed and compressed by their singular values; low-rank adaptation of large models is Eckart–Young applied to the update of a weight matrix.
- **F5 (multivariate Gaussian).** The covariance $\Sigma_{\text{cov}} = Q\Lambda Q^\top$ of a Gaussian is symmetric PSD; its level sets are ellipsoids whose axes are the principal components, with semi-axes $\propto\sqrt{\lambda_i}$. Sampling from $\mathcal N(\mathbf{0}, \Sigma_{\text{cov}})$ is $Q\Lambda^{1/2}\mathbf{z}$: rotate-stretch applied to white noise. The project's 3-D cloud is exactly this.
- **CM-6 / CM-7 (mechanics).** The deformation gradient $F = RU$ of §3.11 is the polar decomposition; the inertia tensor in CM-7 is symmetric PD, so its principal axes are its singular vectors; the stiffness and mass matrices in CM-6 are symmetric PD and the normal-mode problem is their generalized eigenproblem.
- **QM-7 / QC-2 (entanglement).** A two-particle state $|\psi\rangle = \sum_{ij}c_{ij}|i\rangle|j\rangle$ is a matrix $C = (c_{ij})$; its SVD $C = U\Sigma V^\top$ is the *Schmidt decomposition* $|\psi\rangle = \sum_k\sigma_k|u_k\rangle|v_k\rangle$. The state is entangled iff $\operatorname{rank}C > 1$, and the entanglement entropy is $-\sum_k\sigma_k^2\log\sigma_k^2$. The Schmidt decomposition is the SVD, with no changes.
- **F3 / EM-3.** The DFT matrix is unitary (all singular values equal), which is why the Fourier transform neither amplifies nor attenuates. Green's functions discretize to matrices of rapidly decaying singular values (§3.13).

## 6. Exercises

See [exercises.md](exercises.md). **A1–A6** are drills: SVDs you can read off, the rank-one approximation and inverse of Example 4.1, norms and $\kappa$, the shapes of the three forms, the pseudoinverse of a rank-one matrix, and what the singular values do under transpose, scaling, orthogonal multiplication and inversion. **B1–B8** are the module: B2 (full SVD of a $2 \times 3$ matrix by hand through $AA^\top$), B4 (uniqueness and non-uniqueness in Eckart–Young), B5 (the Penrose equations and the normal equations), B6 (Pearson's and Hotelling's PCA agree; what uncentered PCA measures), and B7 (condition numbers in three norms, $\kappa_2(A^\top A) = \kappa_2(A)^2$) are essential; B1 is the norm inequalities you will use constantly. B3 and B8 extend `project.py`. **C1–C3**: the Procrustes problem, Courant–Fischer and Weyl's inequality with a second proof of Eckart–Young, and the distance to singularity with a Hilbert-matrix experiment.

## 7. Coding project

`project.py` builds and tests everything in §3 numerically. In order:

1. **`svd_via_eig(A)`.** The §3.2 derivation as code: `eigh` of $A^\top A$, sort, $\sigma_i = \sqrt{\lambda_i}$, $\mathbf{u}_i = A\mathbf{v}_i/\sigma_i$, complete $U$ by Gram–Schmidt. It is checked against `np.linalg.svd` on random matrices (reconstruction, singular values, and the four subspaces as projectors, since signs will differ), then shown failing on a matrix with singular values $1, 10^{-5}, 10^{-9}$: the small ones come out wrong because $\lambda_3 = 10^{-18}$ is below the rounding noise of $A^\top A$. This is why LA-13 exists.
2. **Worked examples.** Every number in Examples 4.1–4.5 is asserted.
3. **Image compression.** A $256 \times 256$ synthetic image (gradient background, disc, rectangle, diagonal band, sinusoidal texture) is approximated at ranks 1, 5, 20, 50 in a figure grid; the error $\|A - A_k\|_F$ is plotted against $k$ and asserted equal to the Eckart–Young prediction $\sqrt{\sum_{i > k}\sigma_i^2}$ for every $k$. Noise is then added and the scree plot develops a floor.
4. **PCA on a 3-D Gaussian cloud** with a known covariance: the cloud with principal axes drawn with length proportional to the standard deviation $\sigma_i/\sqrt{n-1}$ along each, a scree plot, and the 2-D projection onto the first two components. PCA via `eigh` of the covariance is asserted equal to PCA via SVD of the centered data (eigenvalues to $\sigma_i^2/(n-1)$, eigenvectors up to sign).
5. **Condition number.** Matrices $A = U\operatorname{diag}(1, 1/\kappa)V^\top$ for $\kappa$ from $1$ to $10^{12}$; the worst-case relative error from a fixed relative perturbation of $\mathbf{b}$ is plotted against $\kappa$ (it is a straight line of slope 1 on log–log axes) alongside the actual rounding error of `np.linalg.solve`, and the bound (3.10) is asserted.
6. **Pseudoinverse** on a rank-deficient least-squares problem: $A^+\mathbf{b}$ has the minimal residual, is orthogonal to $N(A)$, and has smaller norm than every other minimizer $A^+\mathbf{b} + \mathbf{z}$ tried. $A^+ = (A^\top A)^{-1}A^\top$ is asserted for full column rank.
7. **Polar decomposition** of random matrices: $Q^\top Q = I$, $S = S^\top \succeq 0$, $QS = A$, and $S = \sqrt{A^\top A}$ via `scipy.linalg.sqrtm`.
8. **$\|A\|_2 = \sigma_1$** by sampling $20\,000$ random unit vectors: the sample maximum never exceeds $\sigma_1$ and gets within a few percent of it; $\mathbf{v}_1$ attains it exactly.

**Build it yourself first:** write `svd_via_eig(A)` (Part 1) and `pca(X)` (Part 4, returning components, explained variance, and scores) before reading the provided versions. The "Your turn" cell asks you to compress your own synthetic pattern, to generate a rank-$k$-plus-noise dataset and recover $k$ from the scree plot, and to solve the orthogonal Procrustes problem with the SVD (Exercise C1 proves why it works).

## 8. Where you'll struggle

- **The derivation everyone sees and nobody remembers.** Students can recite "$\sigma_i^2$ are the eigenvalues of $A^\top A$" and cannot say why the $\mathbf{u}_i$ are orthonormal. The entire theorem is step 3 of §3.2: $(A\mathbf{v}_i)\cdot(A\mathbf{v}_j) = \mathbf{v}_i^\top A^\top A\mathbf{v}_j = \lambda_j\,\mathbf{v}_i\cdot\mathbf{v}_j$. Write that line out until you can produce it cold. If you can, you can reconstruct the whole proof in two minutes; if you cannot, you have memorized a recipe.
- **"NumPy's SVD is wrong."** It is not. The signs of $(\mathbf{u}_i, \mathbf{v}_i)$ pairs are arbitrary (§3.3), and the basis for a repeated singular value is arbitrary. Your hand computation and NumPy's output will disagree in signs more often than not (Example 4.1). Check $U\Sigma V^\top = A$, or check that $\mathbf{u}_i\mathbf{u}_i^\top$ and $\mathbf{v}_i\mathbf{v}_i^\top$ agree, or check the subspace projectors $U_rU_r^\top$. Never check columns entrywise. The related bug: computing $U$ from `eigh(A @ A.T)` and $V$ from `eigh(A.T @ A)` separately and wondering why $U\Sigma V^\top \ne A$.
- **Shapes.** Full ($m \times m$, $m \times n$, $n \times n$), economy ($m \times p$, $p \times p$, $n \times p$ with $p = \min(m,n)$), reduced ($m \times r$, $r \times r$, $n \times r$). NumPy returns `Vt`, not `V`, and `s` as a vector, not $\Sigma$. Reconstructing with `U @ np.diag(s) @ Vt` fails for non-square $A$ in the full form because `np.diag(s)` is $p \times p$ while $\Sigma$ must be $m \times n$. Use `full_matrices=False` for reconstruction, and write the shapes above each product until you stop making this mistake.
- **Eckart–Young's proof is harder than its statement.** The statement is one line; the proof is nine steps. Most courses state it and move on. Do not. Steps 3–4 (split the error by a projector, use orthogonal invariance) and step 8 (a weight budget of $k$ with at most $1$ per slot is best spent on the top $k$) are reusable arguments; you will see the second one again as Ky Fan's maximum principle and in the proof of Weyl's inequality (Exercise C2). Work through the proof with $m = n = 3$, $k = 1$ and concrete $\mathbf{q}_1$ if the general case will not stick.
- **PCA without centering.** If you forget to subtract the column means, $X^\top X/(n-1)$ is not the covariance and the first "principal component" points roughly toward the mean of the data (Example 4.3 shows $36°$ instead of $45°$). Every PCA pipeline must center; many should also standardize. The second PCA failure: comparing eigenvalues of $S$ to singular values of $X$ and forgetting the $\sqrt{n-1}$.
- **$\kappa$ depends on the norm, and is not the determinant.** $\kappa_2 = \sigma_1/\sigma_n$ is the one with the clean geometric meaning; $\kappa_1$ and $\kappa_\infty$ differ from it by factors up to $n$ (Exercise B7). A tiny determinant does not mean ill-conditioned ($\varepsilon I$) and a determinant of $1$ does not mean well-conditioned (Example 4.5's shear has $\det F = 1$ and $\kappa_2 = \phi^2$; with $\gamma = 100$ it has $\det = 1$ and $\kappa_2 \approx 10^4$).
- **Reading the stretch out of the polar decomposition.** In $F = RU$ the mechanics $U$ is a symmetric stretch (our $S = V\Sigma V^\top$), not our orthogonal $U$. Its eigenvectors are the *right* singular vectors of $F$, and its eigenvalues are the singular values. Students who mix the two notations get the principal directions rotated by $R$ (they report the $\mathbf{u}_i$ instead of the $\mathbf{v}_i$). Say "stretch tensor" and "singular-vector matrix" out loud instead of "$U$".
- **Harder than they look.** B2 (the $2 \times 3$ SVD by hand) is routine but long; budget 40 minutes and check $\|A\|_F^2 = \sum\sigma_i^2$ at the end. B4(c) is an equality-case analysis of a nine-step proof; keep the proof open beside you. B6(a) needs Eckart–Young, not just the Rayleigh quotient. C2 is a real proof; the Courant–Fischer step is where most attempts stall.

## 9. Check yourself

1. State the SVD with all shapes for a $5 \times 3$ matrix of rank 2, in full and reduced form. Which columns of $U$ and $V$ span each of the four subspaces?
2. Why are the $\mathbf{u}_i = A\mathbf{v}_i/\sigma_i$ orthonormal? One line.
3. Why do $A^\top A$ and $AA^\top$ have the same nonzero eigenvalues? What are they?
4. In what sense is the SVD "unique"? What is the freedom when $\sigma_2 = \sigma_3$?
5. Prove $\|A\|_2 = \sigma_1$ in three lines. Prove $\|A\|_F^2 = \sum\sigma_i^2$ in one.
6. State Eckart–Young for both norms. What is $\|A - A_k\|_F$? Sketch the two ideas of the Frobenius proof.
7. Derive $\|\delta\mathbf{x}\|/\|\mathbf{x}\| \le \kappa_2\|\delta\mathbf{b}\|/\|\mathbf{b}\|$. For which $\mathbf{b}, \delta\mathbf{b}$ is it an equality? How many digits do you lose when $\kappa_2 = 10^7$?
8. Write $A^+$ in terms of the SVD. Why is $A^+\mathbf{b}$ the minimum-norm least-squares solution? What is $AA^+$?
9. Given centered data $X$ with SVD $U\Sigma V^\top$, what are the principal components, the variances along them, and the scores? Why must $X$ be centered?
10. Derive $A = QS$ from the SVD. In the deformation gradient $F = RU$, which singular vectors of $F$ are the principal material directions, and what are the principal stretches?
11. If $A$ is symmetric with eigenvalues $3, -5, 1$, what are its singular values, and how do $U$ and $V$ relate to the eigenvectors?
12. Give a $2 \times 2$ matrix whose eigenvalues are both $0$ but whose largest singular value is $10^6$.

## 10. Pair with

- **Strang, *Introduction to Linear Algebra* (5th ed.), ch. 7.** §7.1 opens with image compression, §7.2 is the derivation and the four subspaces, §7.3 is PCA, §7.4 the geometry. Strang's version of the "rotate, stretch, rotate" picture and his insistence on $A\mathbf{v}_i = \sigma_i\mathbf{u}_i$ as the thing to remember are the model for §2 and §3.1 here.
- **Strang, *Linear Algebra and Learning from Data*, §I.8–I.9.** §I.8 is the SVD with a proof of Eckart–Young in the spectral norm and a discussion of which norms the theorem holds in; §I.9 is PCA and the "closest rank-$k$ matrix" with real data. Better than these notes on the data-science motivation and on the geometry of $\mathbf{v}_1$ as the direction of maximum variance.
- **MIT 18.065 (Strang), lectures 4–7**, on OCW. Lecture 4 (eigenvalues) and 5 (positive definite matrices) are the LA-8 material in Strang's voice; **lecture 6 is the SVD** and lecture 7 is Eckart–Young. There is no 3Blue1Brown episode on the SVD (the *Essence of Linear Algebra* series stops at eigenvectors and change of basis), so lecture 6 is the video to watch for the picture; Strang draws the ellipse live and does the $A^\top A$ derivation on the board.
- **Trefethen & Bau, *Numerical Linear Algebra*, lectures 4–5.** Lecture 4 is the cleanest written statement of the geometric picture (hyperellipse) and of existence and uniqueness; lecture 5 proves Eckart–Young in the spectral norm by the dimension-counting argument of Theorem 3.8. Short and exact; read after you have done the derivation here.
- **G. W. Stewart, "On the early history of the singular value decomposition," *SIAM Review* 35 (1993).** The source for §1. Ten pages, no prerequisites beyond this module.
