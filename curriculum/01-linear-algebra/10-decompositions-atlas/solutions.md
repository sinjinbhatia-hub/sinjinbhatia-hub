# LA-10 Solutions

## A. Warm-up

**A1.** From (3.3): $l_{11} = \sqrt9 = 3$; $l_{21} = a_{21}/l_{11} = 3/3 = 1$; $l_{22} = \sqrt{a_{22} - l_{21}^2} = \sqrt{5 - 1} = 2$. So $L = \begin{pmatrix} 3 & 0 \\ 1 & 2\end{pmatrix}$, and $LL^\top = \begin{pmatrix} 9 & 3 \\ 3 & 1 + 4\end{pmatrix} = A$. ✓

Forward, $L\mathbf{y} = (12, 8)^\top$: $y_1 = 12/3 = 4$; $y_2 = (8 - 1\cdot4)/2 = 2$. Back, $L^\top\mathbf{x} = \mathbf{y}$: $x_2 = 2/2 = 1$; $x_1 = (4 - 1\cdot1)/3 = 1$. So $\mathbf{x} = (1, 1)^\top$. Check: $A\mathbf{x} = (9 + 3, 3 + 5)^\top = (12, 8)^\top$. ✓

**A2.** Cholesky: $l_{11} = 1$, $l_{21} = 2/1 = 2$, and the radicand for $l_{22}$ is $a_{22} - l_{21}^2 = 1 - 4 = -3 < 0$. It fails at step 2: the leading $2\times2$ block (all of $B$) is not positive definite.

$LDL^\top$ (same recurrence without square roots): $d_1 = 1$; $l_{21} = a_{21}/d_1 = 2$; $d_2 = a_{22} - l_{21}^2d_1 = 1 - 4 = -3$. So
$$L = \begin{pmatrix} 1 & 0 \\ 2 & 1\end{pmatrix}, \quad D = \operatorname{diag}(1, -3), \quad LDL^\top = \begin{pmatrix} 1 & 0 \\ 2 & 1\end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & -3\end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 2 & 4 - 3\end{pmatrix} = B. ✓$$
$B$ is congruent to $D$, so by Sylvester's law of inertia (notes §3.1; LA-12 Theorem 3.12) it has one positive and one negative eigenvalue. Check: $\lambda^2 - 2\lambda - 3 = (\lambda - 3)(\lambda + 1)$, eigenvalues $3$ and $-1$. ✓

**A3.** $n = 1000$, so $\tfrac23 n^3 = 6.67\times10^8$, $\tfrac13 n^3 = 3.33\times10^8$, and one pair of triangular solves is $2n^2 = 2\times10^6$.

(a) $50\times6.67\times10^8 = 3.33\times10^{10}$ flops.
(b) $6.67\times10^8 + 50\times2\times10^6 = 6.67\times10^8 + 1.0\times10^8 = 7.67\times10^8$. Ratio (a)/(b) $\approx 43$.
(c) $3.33\times10^8 + 1.0\times10^8 = 4.33\times10^8$. Ratio (a)/(c) $\approx 77$.

For 1000 right-hand sides: (a) $6.67\times10^{11}$; (b) $6.67\times10^8 + 2\times10^9 = 2.67\times10^9$; ratio $\approx 250$. The triangular solves cost more than the factorization once $2kn^2 > \tfrac23 n^3$, i.e. $k > n/3 \approx 333$ right-hand sides. Beyond that point the factorization is a small fraction of the work, and what matters is the $O(n^2)$ per solve, which no factorization can beat for a dense matrix.

**A4.** Eigenvector for $3$: $(A - 3I)\mathbf{v} = \begin{pmatrix} -1 & 1 \\ 0 & 0\end{pmatrix}\mathbf{v} = \mathbf{0}$ gives $\mathbf{v} = (1, 1)^\top$, so $\mathbf{q}_1 = \tfrac1{\sqrt2}(1, 1)^\top$; complete with $\mathbf{q}_2 = \tfrac1{\sqrt2}(1, -1)^\top$. Then $A\mathbf{q}_1 = 3\mathbf{q}_1$ and $A\mathbf{q}_2 = \tfrac1{\sqrt2}(2 - 1, -3)^\top = \tfrac1{\sqrt2}(1, -3)^\top$. Entries $t_{ij} = \mathbf{q}_i^\top A\mathbf{q}_j$:
$$t_{11} = 3, \quad t_{21} = 0, \quad t_{12} = \tfrac12(1, 1)\cdot(1, -3) = -1, \quad t_{22} = \tfrac12(1, -1)\cdot(1, -3) = 2, \qquad T = \begin{pmatrix} 3 & -1 \\ 0 & 2\end{pmatrix}.$$
$\|A\|_F^2 = 4 + 1 + 9 = 14$ and $\|T\|_F^2 = 9 + 1 + 4 = 14$. ✓ Corollary 3.8a predicts $t_{12}^2 = 14 - (9 + 4) = 1$, and indeed $|t_{12}| = 1$, the same as in the original $T = A$ where $t_{12} = 1$. (Choosing $\mathbf{q}_2 = \tfrac1{\sqrt2}(-1, 1)^\top$ instead gives $t_{12} = +1$; the sign is not invariant, the magnitude is.)

**A5.** $A^\top A = \begin{pmatrix} 0 & 3 \\ -2 & 0\end{pmatrix}\begin{pmatrix} 0 & -2 \\ 3 & 0\end{pmatrix} = \begin{pmatrix} 9 & 0 \\ 0 & 4\end{pmatrix}$, already diagonal, so $P = \operatorname{diag}(3, 2)$ (the unique positive semidefinite square root, Theorem 3.13 step 3). $U = AP^{-1} = \begin{pmatrix} 0 & -2 \\ 3 & 0\end{pmatrix}\operatorname{diag}(\tfrac13, \tfrac12) = \begin{pmatrix} 0 & -1 \\ 1 & 0\end{pmatrix}$. $U^\top U = I$ ✓, $\det U = 1$: a rotation by $+90°$. So $A$ = stretch $x$ by $3$ and $y$ by $2$, then rotate by $90°$.

Left version: $AA^\top = \begin{pmatrix} 0 & -2 \\ 3 & 0\end{pmatrix}\begin{pmatrix} 0 & 3 \\ -2 & 0\end{pmatrix} = \begin{pmatrix} 4 & 0 \\ 0 & 9\end{pmatrix}$, so $P' = \operatorname{diag}(2, 3)$, and $P'U = \begin{pmatrix} 2 & 0 \\ 0 & 3\end{pmatrix}\begin{pmatrix} 0 & -1 \\ 1 & 0\end{pmatrix} = \begin{pmatrix} 0 & -2 \\ 3 & 0\end{pmatrix} = A$. ✓ This is rotate first, then stretch $x$ by $2$ and $y$ by $3$. The relation is $P' = UPU^\top$: the stretch axes have been carried along by the rotation ($U\operatorname{diag}(3, 2)U^\top = \operatorname{diag}(2, 3)$, since $U$ swaps the axes). In general $A = UP = (UPU^\top)U$.

**A6.** (a) **PLU.** $\det A = \det P^{-1}\det L\det U = (\pm1)\prod_iu_{ii}$; $\tfrac23 n^3$ flops, and nothing cheaper exists for a dense matrix. (b) **Cholesky** $\Sigma = LL^\top$. If $\mathbf{z} \sim \mathcal N(\mathbf{0}, I)$ then $L\mathbf{z}$ has covariance $LL^\top = \Sigma$ (F5): one factorization of $\tfrac13 50^3$ flops, then $10^4$ matrix–vector products. (c) **Cholesky of $M$, then spectral** on $L^{-1}KL^{-\top}$ (Example 4.4b, `eigh(K, M)`): symmetric-definite generalized eigenproblem, $M$-orthogonal modes. (d) **SVD** of the centered $10^5\times3$ data matrix: the plane's normal is the right singular vector of the smallest $\sigma$ (total least squares; this is PCA, LA-9). If instead you fit $z = ax + by + c$ by ordinary least squares, **QR** of the $10^5\times3$ matrix $[\mathbf{x}, \mathbf{y}, \mathbf{1}]$. (e) **SVD**: look for a gap after $\sigma_3$; rank is not a question LU or eigenvalues can answer reliably. (f) **Polar via SVD** (Procrustes): center both point sets, $M = YX^\top$ ($3\times3$), $R = W\operatorname{diag}(1, 1, \det(WV^\top))V^\top$.

## B. Core

**B1.** From $L_1L_1^\top = L_2L_2^\top$, multiply on the left by $L_2^{-1}$ and on the right by $L_1^{-\top}$ (both exist: triangular with nonzero diagonal):
$$L_2^{-1}L_1 = L_2^\top L_1^{-\top}.$$
The left side is a product of lower triangular matrices (the inverse of a lower triangular matrix is lower triangular, by back substitution column by column), hence lower triangular. The right side is $L_2^\top$ (upper) times $L_1^{-\top} = (L_1^{-1})^\top$ (upper), hence upper triangular. A matrix that is both has zeros strictly below *and* strictly above the diagonal: it is diagonal. Call it $\Delta$. Then $L_1 = L_2\Delta$.

Substitute into $A$: $L_2L_2^\top = L_1L_1^\top = L_2\Delta\Delta^\top L_2^\top = L_2\Delta^2L_2^\top$. Cancel the invertible $L_2$ on the left and $L_2^\top$ on the right: $\Delta^2 = I$, so each $\delta_{ii} = \pm1$. Now the diagonal of $L_1 = L_2\Delta$ is $(L_1)_{ii} = (L_2)_{ii}\delta_{ii}$. Both $(L_1)_{ii}$ and $(L_2)_{ii}$ are positive, so $\delta_{ii} = +1$. Hence $\Delta = I$ and $L_1 = L_2$. ∎

Positivity was used only in the last sentence. Without it, $L$ is unique up to the signs of its columns ($L \mapsto L\Delta$ with $\Delta = \operatorname{diag}(\pm1)$), which is exactly the freedom in $A^\top A = R^\top R$ that B2 fixes.

**B2.** (a) $A^\top A = (QR)^\top(QR) = R^\top Q^\top QR = R^\top R$ because $Q^\top Q = I_n$ (orthonormal columns). $A$ has full column rank, so $R$ is invertible and its diagonal is nonzero. Flipping the sign of row $i$ of $R$ and column $i$ of $Q$ leaves $QR$ unchanged, so we may take $\operatorname{diag}R > 0$. Then $R^\top$ is lower triangular with positive diagonal, $A^\top A = R^\top R$ has the Cholesky form, and by B1 it *is* the Cholesky factorization. (So the normal-equations route and the QR route compute the same triangular matrix, one of them from a squared matrix.)

(b) Let $A = U\Sigma V^\top$ with $\sigma_1 \ge \cdots \ge \sigma_n > 0$. Then $A^\top A = V\Sigma^\top U^\top U\Sigma V^\top = V(\Sigma^\top\Sigma)V^\top$, an orthogonal diagonalization with eigenvalues $\sigma_i^2$. For an SPD matrix the singular values are the eigenvalues, so $\kappa_2(A^\top A) = \sigma_1^2/\sigma_n^2 = \kappa_2(A)^2$. ∎

(c) QR route: $R$ has the same singular values as $A$ (because $A = QR$ with $Q^\top Q = I$ gives $R^\top R = A^\top A$, so $\sigma_i(R)^2 = \sigma_i(A)^2$), hence $\kappa_2(R) = 10^7$; the triangular solve loses about 7 digits and leaves about 9. Normal equations: $\kappa_2(R^\top R) = 10^{14}$; the solve loses about 14 digits and leaves about 2. (Forming $A^\top A$ in floating point also rounds, and the loss there alone can be fatal.) Honest caveat for LA-13: the least-squares *problem* can itself be sensitive to perturbations at a level involving $\kappa^2$ when the residual is large (Trefethen & Bau, lecture 18); Householder QR is backward stable and delivers whatever accuracy the problem allows, whereas the normal equations lose digits even for well-posed problems.

**B3.** (a) Using $\|B\|_F^2 = \operatorname{tr}(B^*B)$ and $\operatorname{tr}(XY) = \operatorname{tr}(YX)$:
$$\|A\|_F^2 = \operatorname{tr}\big((QTQ^*)^*(QTQ^*)\big) = \operatorname{tr}(QT^*Q^*QTQ^*) = \operatorname{tr}(QT^*TQ^*) = \operatorname{tr}(T^*TQ^*Q) = \operatorname{tr}(T^*T) = \|T\|_F^2.$$
$T$ is upper triangular, so $\|T\|_F^2 = \sum_{i\le j}|t_{ij}|^2 = \sum_i|t_{ii}|^2 + \sum_{i<j}|t_{ij}|^2$, and $t_{ii} = \lambda_i$ (Theorem 3.8). Hence $\sum_{i<j}|t_{ij}|^2 = \|A\|_F^2 - \sum_i|\lambda_i|^2$, a quantity depending only on $A$.

$A$ normal $\iff$ $T$ normal (unitary similarity preserves $TT^* = T^*T$, as in the proof of Theorem 3.10) $\iff$ $T$ diagonal (Lemma 3.9 one way; a diagonal matrix commutes with its conjugate transpose the other way) $\iff$ $\sum_{i<j}|t_{ij}|^2 = 0$ $\iff$ $\|A\|_F^2 = \sum_i|\lambda_i|^2$. ∎

(b) $\begin{pmatrix} 1 & 2 \\ 3 & 0\end{pmatrix}$: $\|A\|_F^2 = 1 + 4 + 9 + 0 = 14$, eigenvalues $3, -2$, $\sum|\lambda_i|^2 = 13$, so $\sum_{i<j}|t_{ij}|^2 = 1$ (matching Example 4.2). $\begin{pmatrix} 0 & 1 \\ 0 & 0\end{pmatrix}$: $\|A\|_F^2 = 1$, both eigenvalues $0$, so the strictly upper part has squared size $1$: this matrix is already its own Schur form and is as far from normal as a matrix of Frobenius norm 1 can be.

(c) Any nilpotent $2\times2$ matrix of Frobenius norm $10$: $\begin{pmatrix} 0 & 10 \\ 0 & 0\end{pmatrix}$, or $\begin{pmatrix} 5 & -5 \\ 5 & -5\end{pmatrix}$ (trace $0$, determinant $-25 + 25 = 0$, so both eigenvalues are $0$; $\|A\|_F^2 = 4\cdot25 = 100$). By (a), its Schur form has $|t_{12}| = 10$.

**B4.** Write $A = Q\Lambda Q^*$ with $Q$ unitary and $\Lambda = \operatorname{diag}(\lambda_i)$.

(a) All $\lambda_i$ real means $\Lambda^* = \bar\Lambda = \Lambda$, so $A^* = Q\Lambda^*Q^* = Q\Lambda Q^* = A$. ∎

(b) $|\lambda_i| = 1$ means $\Lambda^*\Lambda = \operatorname{diag}(|\lambda_i|^2) = I$, so $A^*A = Q\Lambda^*Q^*Q\Lambda Q^* = Q\Lambda^*\Lambda Q^* = QQ^* = I$. ∎

(c) $\lambda_i \in \{0, 1\}$ gives $\Lambda^2 = \Lambda$ and $\Lambda$ real. So $A^2 = Q\Lambda^2Q^* = Q\Lambda Q^* = A$, and $A^* = A$ by (a). An idempotent Hermitian matrix is the orthogonal projection onto its column space, here the span of the columns of $Q$ with $\lambda_i = 1$. ∎

(d) Directly from normality: $\|A\mathbf{x}\|^2 = \mathbf{x}^*A^*A\mathbf{x} = \mathbf{x}^*AA^*\mathbf{x} = \|A^*\mathbf{x}\|^2$. Via $\Lambda$, with $\mathbf{z} = Q^*\mathbf{x}$: $\|A\mathbf{x}\|^2 = \|Q\Lambda\mathbf{z}\|^2 = \|\Lambda\mathbf{z}\|^2 = \sum_i|\lambda_i|^2|z_i|^2 = \|\bar\Lambda\mathbf{z}\|^2 = \|Q\bar\Lambda Q^*\mathbf{x}\|^2 = \|A^*\mathbf{x}\|^2$. ∎ (Conversely, $\|A\mathbf{x}\| = \|A^*\mathbf{x}\|$ for all $\mathbf{x}$ implies normality, by polarization; this is sometimes taken as the definition.)

**B5.** (a) $A$, $\mathbf{x}$, $\mathbf{y}$ are real, so $A(\mathbf{x} + i\mathbf{y}) = A\mathbf{x} + iA\mathbf{y}$ has real part $A\mathbf{x}$ and imaginary part $A\mathbf{y}$. The right side $(a + ib)(\mathbf{x} + i\mathbf{y}) = (a\mathbf{x} - b\mathbf{y}) + i(b\mathbf{x} + a\mathbf{y})$. Equate: $A\mathbf{x} = a\mathbf{x} - b\mathbf{y}$ and $A\mathbf{y} = b\mathbf{x} + a\mathbf{y}$. ∎

(b) Both $A\mathbf{x}$ and $A\mathbf{y}$ are combinations of $\mathbf{x}$ and $\mathbf{y}$, so $A$ maps $W = \operatorname{span}\{\mathbf{x}, \mathbf{y}\}$ into itself. The matrix of $A|_W$ in the basis $(\mathbf{x}, \mathbf{y})$ has as columns the coordinates of $A\mathbf{x}$ and $A\mathbf{y}$: $A\mathbf{x} = a\mathbf{x} + (-b)\mathbf{y}$ gives the column $(a, -b)^\top$, and $A\mathbf{y} = b\mathbf{x} + a\mathbf{y}$ gives $(b, a)^\top$. So the matrix is $\begin{pmatrix} a & b \\ -b & a\end{pmatrix}$, with characteristic polynomial $(\lambda - a)^2 + b^2$ and eigenvalues $a \pm ib$, as it must be.

(c) Suppose $\mathbf{x}, \mathbf{y}$ are dependent. If $\mathbf{x} \ne \mathbf{0}$, then $\mathbf{y} = c\mathbf{x}$ for a real $c$, and the eigenvector is $\mathbf{v} = (1 + ic)\mathbf{x}$. From $A\mathbf{v} = \lambda\mathbf{v}$, dividing by the nonzero scalar $1 + ic$: $A\mathbf{x} = \lambda\mathbf{x}$. The left side is a real vector and $\mathbf{x} \ne \mathbf{0}$ is real, so $\lambda$ is real, contradicting $b \ne 0$. If $\mathbf{x} = \mathbf{0}$, then $\mathbf{v} = i\mathbf{y}$ with $\mathbf{y} \ne \mathbf{0}$ real and $A\mathbf{y} = \lambda\mathbf{y}$, the same contradiction. So $\mathbf{x}, \mathbf{y}$ are independent and $W$ is a genuine plane. ∎

(d) In the proof of Theorem 3.8, replace "pick a unit eigenvector $\mathbf{q}_1$" by: if $\lambda_1$ is real, proceed as before with a real eigenvector; if $\lambda_1 = a + ib$ is not real, take an orthonormal basis $\mathbf{q}_1, \mathbf{q}_2$ of the invariant plane $W$ (Gram–Schmidt on $\mathbf{x}, \mathbf{y}$) and complete it to a real orthogonal $Q_1 = \begin{pmatrix} \mathbf{q}_1 & \mathbf{q}_2 & \hat Q_1\end{pmatrix}$. Because $A\mathbf{q}_1, A\mathbf{q}_2 \in W$, the first two columns of $Q_1^\top AQ_1$ are zero below row 2, so $Q_1^\top AQ_1 = \begin{pmatrix} B & * \\ 0 & A_2\end{pmatrix}$ with $B$ a real $2\times2$ block similar to $\begin{pmatrix} a & b \\ -b & a\end{pmatrix}$ (eigenvalues $a \pm ib$) and $A_2$ real $(n-2)\times(n-2)$. Deflate by two and recurse on $A_2$. The result is a real orthogonal $Q$ and a real quasi-triangular $T$. ∎

**B6.** (a) $A^\top A = \begin{pmatrix} 3 & 4 \\ 0 & 5\end{pmatrix}\begin{pmatrix} 3 & 0 \\ 4 & 5\end{pmatrix} = \begin{pmatrix} 25 & 20 \\ 20 & 25\end{pmatrix}$. Eigenvalues $25 \pm 20 = 45, 5$ with eigenvectors $\mathbf{v}_1 = \tfrac1{\sqrt2}(1, 1)^\top$, $\mathbf{v}_2 = \tfrac1{\sqrt2}(1, -1)^\top$. Singular values $\sigma_1 = \sqrt{45} = 3\sqrt5$, $\sigma_2 = \sqrt5$. Then
$$P = V\operatorname{diag}(3\sqrt5, \sqrt5)V^\top = \frac12\begin{pmatrix} 3\sqrt5 + \sqrt5 & 3\sqrt5 - \sqrt5 \\ 3\sqrt5 - \sqrt5 & 3\sqrt5 + \sqrt5\end{pmatrix} = \sqrt5\begin{pmatrix} 2 & 1 \\ 1 & 2\end{pmatrix}.$$
Check: $P^2 = 5\begin{pmatrix} 5 & 4 \\ 4 & 5\end{pmatrix} = \begin{pmatrix} 25 & 20 \\ 20 & 25\end{pmatrix}$. ✓

(b) $P^{-1} = \dfrac{1}{\sqrt5}\cdot\dfrac13\begin{pmatrix} 2 & -1 \\ -1 & 2\end{pmatrix}$ (determinant of $\begin{pmatrix} 2 & 1 \\ 1 & 2\end{pmatrix}$ is $3$). So
$$U = AP^{-1} = \frac{1}{3\sqrt5}\begin{pmatrix} 3 & 0 \\ 4 & 5\end{pmatrix}\begin{pmatrix} 2 & -1 \\ -1 & 2\end{pmatrix} = \frac{1}{3\sqrt5}\begin{pmatrix} 6 & -3 \\ 8 - 5 & -4 + 10\end{pmatrix} = \frac1{\sqrt5}\begin{pmatrix} 2 & -1 \\ 1 & 2\end{pmatrix}.$$
$U^\top U = \tfrac15\begin{pmatrix} 4 + 1 & -2 + 2 \\ -2 + 2 & 1 + 4\end{pmatrix} = I$. ✓ $UP = \begin{pmatrix} 2 & -1 \\ 1 & 2\end{pmatrix}\begin{pmatrix} 2 & 1 \\ 1 & 2\end{pmatrix} = \begin{pmatrix} 4 - 1 & 2 - 2 \\ 2 + 2 & 1 + 4\end{pmatrix} = \begin{pmatrix} 3 & 0 \\ 4 & 5\end{pmatrix} = A$. ✓

(c) By (3.6) with $X = I$: for orthogonal $Q$, $\|A - Q\|_F^2 = \|A\|_F^2 + 2 - 2\operatorname{tr}(Q^\top A)$, and Theorem 3.14 says $\operatorname{tr}(Q^\top A) \le \sigma_1 + \sigma_2$ with equality iff $Q = WV^\top = U$. So $U$ minimizes the distance. Its value: $\|A - U\|_F^2 = \|U(P - I)\|_F^2 = \|P - I\|_F^2$ (Frobenius norm is unitarily invariant) $= \|V(\Sigma - I)V^\top\|_F^2 = \sum_i(\sigma_i - 1)^2 = (3\sqrt5 - 1)^2 + (\sqrt5 - 1)^2 = (46 - 6\sqrt5) + (6 - 2\sqrt5) = 52 - 8\sqrt5 \approx 34.11$. Cross-check with (3.6): $\|A\|_F^2 + 2 - 2(\sigma_1 + \sigma_2) = 50 + 2 - 8\sqrt5$. ✓

(d) $U = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{pmatrix}$ with $\cos\theta = 2/\sqrt5$, $\sin\theta = 1/\sqrt5$: $\theta = \arctan\tfrac12 \approx 26.57°$. So $A$ is "stretch by $3\sqrt5$ along $(1, 1)$ and by $\sqrt5$ along $(1, -1)$, then rotate by $26.57°$."

**B7.** Reference implementation (also in `project.py`, Part 1):

```python
def ldlt(A):
    """A = L diag(d) L^T, L unit lower triangular (Theorem 3.3), no pivoting."""
    A = np.asarray(A, float); n = A.shape[0]
    L = np.eye(n); d = np.zeros(n)
    for j in range(n):
        d[j] = A[j, j] - (L[j, :j] ** 2) @ d[:j]
        if d[j] == 0:
            raise ZeroDivisionError(f"pivot {j+1} is exactly zero; symmetric pivoting needed")
        L[j + 1:, j] = (A[j + 1:, j] - (L[j + 1:, :j] * d[:j]) @ L[j, :j]) / d[j]
    return L, d
```

Derivation of the recurrence: $(LDL^\top)_{ij} = \sum_{k\le\min(i,j)}l_{ik}d_kl_{jk}$. For $j = i$: $a_{jj} = \sum_{k<j}l_{jk}^2d_k + d_j$ (using $l_{jj} = 1$). For $i > j$: $a_{ij} = \sum_{k<j}l_{ik}d_kl_{jk} + l_{ij}d_j$. Solve for $d_j$ then $l_{ij}$.

(a) On the Example 4.1 matrix: $d_1 = 4$; $l_{21} = l_{31} = 2/4 = \tfrac12$; $d_2 = 5 - (\tfrac12)^2\cdot4 = 4$; $l_{32} = (3 - \tfrac12\cdot4\cdot\tfrac12)/4 = \tfrac12$; $d_3 = 6 - \tfrac14\cdot4 - \tfrac14\cdot4 = 4$. So $d = (4, 4, 4)$ and $L = \begin{pmatrix} 1 & 0 & 0 \\ 1/2 & 1 & 0 \\ 1/2 & 1/2 & 1\end{pmatrix} = L_{\text{Chol}}/2$. ✓ (Since $L_{\text{Chol}} = LD^{1/2} = L\cdot2I$.)

(b) $\begin{pmatrix} 1 & 2 \\ 2 & 1\end{pmatrix}$: $d = (1, -3)$, one negative; eigenvalues $3, -1$, one negative. ✓ Modified Example 4.1 ($a_{33} = 1$): columns 1–2 unchanged, $d_3 = 1 - 1 - 1 = -1$, so $d = (4, 4, -1)$, one negative; eigenvalues $\approx 8.33, 2.45, -0.78$, one negative. ✓ This is Sylvester's law of inertia (LA-12 Theorem 3.12): $A = LDL^\top$ is a congruence, and congruence preserves the counts of positive, negative and zero eigenvalues.

(c) By Theorem 3.4, `cholesky` needs $l_{jj} = \sqrt{d_j}$ to be real and positive, so it fails exactly when some $d_j \le 0$. `ldlt` never takes a square root; it only divides by $d_j$, so it survives negative pivots and fails only on an exactly zero pivot. That is why $LDL^\top$ is the factorization for symmetric *indefinite* matrices (with pivoting, below).

(d) $\begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix}$ is invertible ($\det = -1$) but $d_1 = a_{11} = 0$, so unpivoted `ldlt` divides by zero. Symmetric pivoting swaps row $i$ with row $j$ *and* column $i$ with column $j$ (the matrix $\Pi A\Pi^\top$, still symmetric), bringing a different diagonal entry to the pivot position. Here it cannot help: every diagonal entry is $0$. That is why the Bunch–Kaufman algorithm allows $2\times2$ *blocks* on the diagonal of $D$; for this matrix $D = A$ itself is one $2\times2$ block and $L = I$. With $1\times1$ pivots only, $LDL^\top$ exists iff no leading principal minor vanishes after the permutation, which cannot always be arranged.

**B8.** (a) $K_3 = \begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1\end{pmatrix}$. $d_1 = 2$; $l_{21} = -1/2$, $l_{31} = 0$; $d_2 = 2 - l_{21}^2d_1 = 2 - \tfrac14\cdot2 = \tfrac32$; $l_{32} = (a_{32} - l_{31}d_1l_{21})/d_2 = (-1 - 0)/\tfrac32 = -\tfrac23$; $d_3 = 1 - l_{31}^2d_1 - l_{32}^2d_2 = 1 - 0 - \tfrac49\cdot\tfrac32 = 1 - \tfrac23 = \tfrac13$. Pivots $(2, \tfrac32, \tfrac13)$, all positive, so $K_3$ is SPD, and $\det K_3 = 2\cdot\tfrac32\cdot\tfrac13 = 1$. (Direct check: $2(2 - 1) - (-1)(-1 - 0) = 2 - 1 = 1$. ✓)

(b) Eliminating the first variable of a tridiagonal matrix with first pivot $d$ and first off-diagonal $-1$ changes only the $(2,2)$ entry of what remains: the Schur complement (3.2) is $C - \mathbf{b}\mathbf{b}^\top/d$ with $\mathbf{b} = (-1, 0, \dots, 0)^\top$, and $\mathbf{b}\mathbf{b}^\top/d$ has a single nonzero entry $1/d$ in position $(1,1)$ of $C$. So after $k - 1$ steps the remaining matrix is again tridiagonal with $-1$ off the diagonal, diagonal entries $(d_k, 2, \dots, 2, 1)$, and the recurrence is $d_{k+1} = 2 - 1/d_k$ for $k + 1 < n$ and $d_n = 1 - 1/d_{n-1}$. Induction: $d_1 = 2 = 2/1$; if $d_k = (k+1)/k$ then $d_{k+1} = 2 - k/(k+1) = (k+2)/(k+1)$. ✓ Finally $d_n = 1 - (n-1)/n = 1/n$. The determinant is the product of the pivots: $\prod_{k<n}\frac{k+1}{k}\cdot\frac1n = \frac{n}{1}\cdot\frac1n = 1$ (the product telescopes). ∎ (For $n = 3$: $2\cdot\tfrac32\cdot\tfrac13 = 1$. ✓)

(c) Free–free chain: diagonal $(1, 2, \dots, 2, 1)$. Same recurrence with $d_1 = 1$: $d_2 = 2 - 1/1 = 1$, $d_3 = 2 - 1 = 1$, and so on: $d_k = 1$ for all $k < n$. Then $d_n = 1 - 1/d_{n-1} = 1 - 1 = 0$. The radicand at the last Cholesky step is exactly $0$: $K$ is positive *semi*definite and singular. Null vector: every row of $K$ sums to zero, so $K\mathbf{1} = \mathbf{0}$. Mechanically, $\mathbf{u} = \mathbf{1}$ displaces every mass by the same amount, stretches no spring, and stores no energy ($\tfrac12\mathbf{1}^\top K\mathbf{1} = 0$): the rigid-body translation mode. Cholesky found it at the last step because the first $n - 1$ leading blocks are the stiffness matrices of *constrained* sub-chains and are positive definite; only the full free structure is singular.

(d) $K$ is SPD (if constrained) with bandwidth $1$ here and small bandwidth $w$ for any mesh after good ordering. Cholesky without pivoting preserves the band: in (3.3), $l_{ij}$ for $i > j + w$ is a combination of entries $a_{ij} = 0$ and products of already-zero $l$'s, so it is zero, and $L$ has the same bandwidth as $K$; the cost drops from $\tfrac13 n^3$ to about $nw^2$, and the storage from $n^2$ to $nw$. For $n = 10^5$ that is the difference between $3\times10^{14}$ flops and (say, $w = 300$) $10^{10}$. Partial pivoting in LU chooses rows by magnitude, not by structure: a row swap can move entries far from the diagonal and create fill-in that destroys the band (for a tridiagonal matrix the damage is mild, bandwidth $2$ in $U$; for a 2D or 3D mesh it is severe and unpredictable). Cholesky needs no pivoting (Theorem 3.5), so the nonzero pattern of $L$ is known before any arithmetic, which is what sparse direct solvers exploit. Twice fewer flops is a bonus; the predictable structure is the reason.

## C. Challenge

**C1.** (a) $A$ normal gives $A = Q\Lambda Q^*$ (Theorem 3.10). Then $A^*A = Q\bar\Lambda\Lambda Q^* = Q\operatorname{diag}(|\lambda_i|^2)Q^*$, so the eigenvalues of $A^*A$ are $|\lambda_i|^2$ and the singular values of $A$ are $|\lambda_i|$ (in decreasing order). Explicitly, writing $\lambda_i = |\lambda_i|e^{i\phi_i}$ and $D_\phi = \operatorname{diag}(e^{i\phi_i})$ (unitary), $A = (QD_\phi)\operatorname{diag}(|\lambda_i|)Q^*$ is an SVD. Hence $\|A\|_2 = \sigma_1 = \max_i|\lambda_i|$. ∎

(b) For any eigenpair $A\mathbf{v} = \lambda\mathbf{v}$ with $\|\mathbf{v}\| = 1$: $|\lambda| = \|A\mathbf{v}\| \le \|A\|_2\|\mathbf{v}\| = \|A\|_2$. ∎ Example: $A = \begin{pmatrix} 1 & 10^6 \\ 0 & 1\end{pmatrix}$ has both eigenvalues $1$, but $\|A\|_2 \ge \|A\mathbf{e}_2\| = \sqrt{10^{12} + 1} > 10^6$. (Or $\begin{pmatrix} 0 & 10^6 \\ 0 & 0\end{pmatrix}$, where the ratio is infinite.)

(c) With $a = 0.9$, $b = 1000$ and $N = \begin{pmatrix} 0 & 1 \\ 0 & 0\end{pmatrix}$, $A = aI + bN$ and $N^2 = 0$, so by the binomial theorem (valid because $I$ and $N$ commute) $A^k = a^kI + ka^{k-1}bN = \begin{pmatrix} a^k & ka^{k-1}b \\ 0 & a^k\end{pmatrix}$. The $(1,2)$ entry $f(k) = 1000\,k\,0.9^{k-1}$ dominates the norm. It grows while $f(k+1)/f(k) = 0.9(k+1)/k > 1$, i.e. while $k < 9$, and then decays; $f(9) = f(10) = 1000\cdot9\cdot0.9^8 \approx 3874$. So $\|A^k\|_2$ rises from about $1000$ at $k = 1$ to about $3.9\times10^3$ at $k = 9$–$10$ before decaying like $k\,0.9^k \to 0$ (the project prints the exact values). For the iteration $\mathbf{x}_{k+1} = A\mathbf{x}_k$, "spectral radius $< 1$" guarantees eventual decay but says nothing about the transient: a perturbation can be amplified by thousands before it dies, which in a control loop or a time-stepping scheme is the difference between stable and useless. For $\dot{\mathbf{x}} = B\mathbf{x}$ the same happens with $e^{tB}$: $\operatorname{Re}\lambda_i < 0$ guarantees $\|e^{tB}\| \to 0$, but non-normal $B$ (non-symmetric damping, gyroscopic terms, shear flows) can have $\|e^{tB}\| \gg 1$ for intermediate $t$. For normal $A$ none of this is possible: by (a) applied to $A^k$ (also normal, with eigenvalues $\lambda_i^k$), $\|A^k\|_2 = (\max_i|\lambda_i|)^k$ exactly, monotone from the start, and $\|e^{tB}\|_2 = e^{t\max_i\operatorname{Re}\lambda_i}$.

**C2.** (a) $f(s, R, \mathbf{t}) = \sum_{i=1}^k\|sR\mathbf{x}_i + \mathbf{t} - \mathbf{y}_i\|^2$. Setting the gradient with respect to $\mathbf{t}$ to zero: $2\sum_i(sR\mathbf{x}_i + \mathbf{t} - \mathbf{y}_i) = \mathbf{0}$, so $\mathbf{t}^* = \bar{\mathbf{y}} - sR\bar{\mathbf{x}}$ with $\bar{\mathbf{x}} = \tfrac1k\sum_i\mathbf{x}_i$, $\bar{\mathbf{y}} = \tfrac1k\sum_i\mathbf{y}_i$ ($f$ is a convex quadratic in $\mathbf{t}$, so this is the minimum). Substituting, $f = \sum_i\|sR(\mathbf{x}_i - \bar{\mathbf{x}}) - (\mathbf{y}_i - \bar{\mathbf{y}})\|^2$: the same problem for the centered sets, with no translation. So assume $X, Y$ centered from now on.

(b) Expanding as in (3.6) with $sR$ in place of $R$ and using $R^\top R = I$: $f(s, R) = s^2\|X\|_F^2 - 2s\operatorname{tr}(R^\top M) + \|Y\|_F^2$, $M = YX^\top$. For fixed $R$ this is a quadratic in $s$ with positive leading coefficient, minimized at $s^*(R) = \operatorname{tr}(R^\top M)/\|X\|_F^2$.

(c) Substituting $s^*(R)$: $f(s^*(R), R) = \|Y\|_F^2 - \operatorname{tr}(R^\top M)^2/\|X\|_F^2$. Minimizing over $R$ means maximizing $\operatorname{tr}(R^\top M)^2$. By Theorem 3.14, $\operatorname{tr}(R^\top M)$ ranges over $[-\sum_i\sigma_i, \sum_i\sigma_i]$ (the lower end at $R = -WV^\top$), so the maximum of the square is $(\sum_i\sigma_i)^2$, attained at $R = \pm WV^\top$; the constraint $s > 0$ selects $R = WV^\top$ with $s^* = (\sum_i\sigma_i)/\|X\|_F^2 > 0$. Minimum residual: $\|Y\|_F^2 - (\sum_i\sigma_i)^2/\|X\|_F^2$. ∎ Note that $R$ does not depend on the scale at all: replacing $Y$ by $cY$ multiplies $M$ by $c$, which leaves $W$ and $V$ unchanged for $c > 0$. The rotation estimate of plain Procrustes is already scale-invariant; only the residual and $s$ change.

(d) Code (the "Your turn" cell in `project.py` gives the setup):

```python
def procrustes_scaled(X, Y):
    xbar = X.mean(1, keepdims=True); ybar = Y.mean(1, keepdims=True)
    Xc, Yc = X - xbar, Y - ybar
    W, s, Vt = np.linalg.svd(Yc @ Xc.T)
    D = np.eye(X.shape[0]); D[-1, -1] = np.sign(np.linalg.det(W @ Vt))
    R = W @ D @ Vt
    scale = (s * np.diag(D)).sum() / np.sum(Xc**2)
    t = (ybar - scale * R @ xbar).ravel()
    return scale, R, t
```
On the bracket cloud of Part 3 with $Y = 1.7R_{\text{true}}X + \mathbf{t} + \text{noise}$ ($\sigma = 0.02$), this returns $s \approx 1.700$ (three digits), with the rotation error below $0.5°$ as before. (With the determinant fix, the trace that enters $s^*$ is $\sum_i\sigma_id_{ii}$, which is what the code computes.)

**C3.** (a) $A = LL^\top$ gives $a_{jj} = (LL^\top)_{jj} = \sum_kl_{jk}l_{jk} = \sum_{k\le j}l_{jk}^2$ (the sum stops at $k = j$ because $L$ is lower triangular). Every term is $\ge 0$ and one of them is $l_{jj}^2$, so $a_{jj} \ge l_{jj}^2$, with equality iff $l_{jk} = 0$ for all $k < j$.

(b) $\det A = \det L\det L^\top = (\det L)^2 = \prod_jl_{jj}^2 \le \prod_ja_{jj}$ by (a), since all factors are positive. Equality iff equality holds in (a) for every $j$, i.e. iff every off-diagonal entry of $L$ is zero, i.e. $L$ is diagonal, i.e. $A = LL^\top$ is diagonal. ∎

(c) If $B$ is singular, $\det B = 0$ and there is nothing to prove. Otherwise $B^\top B$ is SPD (for $\mathbf{x} \ne \mathbf{0}$, $\mathbf{x}^\top B^\top B\mathbf{x} = \|B\mathbf{x}\|^2 > 0$), with $(B^\top B)_{jj} = \mathbf{b}_j^\top\mathbf{b}_j = \|\mathbf{b}_j\|^2$ and $\det(B^\top B) = (\det B)^2$. By (b), $(\det B)^2 \le \prod_j\|\mathbf{b}_j\|^2$; take square roots. ∎ Geometry: $|\det B|$ is the volume of the parallelepiped spanned by the columns $\mathbf{b}_j$ (LA-6). Hadamard's inequality says that volume is at most the product of the edge lengths, with equality iff the edges are mutually orthogonal (a rectangular box), because equality in (b) means $B^\top B$ is diagonal, i.e. $\mathbf{b}_i^\top\mathbf{b}_j = 0$ for $i \ne j$. Tilting a box never increases its volume.

**C4.** (a) The diagonal entries of $T_\varepsilon$ are $t_{ii} + i\varepsilon$. Two of them coincide, $t_{ii} + i\varepsilon = t_{jj} + j\varepsilon$ with $i < j$, iff $\varepsilon = (t_{ii} - t_{jj})/(j - i)$: at most one bad $\varepsilon$ per pair, so at most $\binom n2$ bad values in total. For every other $\varepsilon > 0$, $T_\varepsilon$ is upper triangular with $n$ distinct diagonal entries, hence has $n$ distinct eigenvalues, hence is diagonalizable (LA-7: eigenvectors for distinct eigenvalues are independent). $A_\varepsilon = QT_\varepsilon Q^*$ is similar to $T_\varepsilon$, so it is diagonalizable too. And $A_\varepsilon - A = Q(\varepsilon\operatorname{diag}(1, \dots, n))Q^*$ has 2-norm $\varepsilon n$ (unitary similarity preserves the 2-norm, and the norm of a diagonal matrix is its largest entry). ∎

(b) Let $p_\varepsilon(\lambda) = \det(A_\varepsilon - \lambda I)$. For good $\varepsilon$, $A_\varepsilon$ is diagonalizable, so $p_\varepsilon(A_\varepsilon) = 0$ by the one-liner. Now let $\varepsilon \to 0$ through good values (possible: only finitely many are bad). The coefficients of $p_\varepsilon$ are sums of products of entries of $A_\varepsilon$ (a determinant is a polynomial in the entries), hence continuous in those entries, and $A_\varepsilon \to A$ entrywise, so each coefficient of $p_\varepsilon$ tends to the corresponding coefficient of $p_A$. Matrix addition and multiplication are continuous, so $p_\varepsilon(A_\varepsilon) = \sum_mc_m(\varepsilon)A_\varepsilon^m \to \sum_mc_mA^m = p_A(A)$. A limit of zero matrices is the zero matrix: $p_A(A) = 0$. ∎ (LA-11 proves Cayley–Hamilton a different way, directly from the triangular Schur factor: $(T - \lambda_kI)$ maps $\operatorname{span}(\mathbf{e}_1, \dots, \mathbf{e}_k)$ into $\operatorname{span}(\mathbf{e}_1, \dots, \mathbf{e}_{k-1})$, so the product of all $n$ factors is zero. That proof needs no limit; this one needs only Schur plus continuity.)

(c) For every good $\varepsilon$, $A_\varepsilon$ is diagonalizable, so its Jordan form is diagonal (no block larger than $1\times1$). If the Jordan form were continuous, the Jordan form of $A$ would also be diagonal, i.e. every matrix would be diagonalizable, which is false ($\begin{pmatrix} 0 & 1 \\ 0 & 0\end{pmatrix}$). So the Jordan form is *discontinuous* at every defective matrix: the Jordan forms of $A_\varepsilon$ do not converge to the Jordan form of $A$, even though $A_\varepsilon \to A$. In floating point every matrix *is* an $A_\varepsilon$ for some rounding-sized $\varepsilon$, so a computed Jordan form would report "diagonalizable with nearly equal eigenvalues" for every defective matrix, which is the correct answer to the wrong question. Continuous quantities (eigenvalues as a multiset, $p_A$, invariant subspaces for well-separated eigenvalue clusters, the Schur form up to its non-uniqueness) can be computed; the Jordan block structure cannot. This is the content of §3.8 and the reason LA-11 treats Jordan as a tool for proofs, not computation.
