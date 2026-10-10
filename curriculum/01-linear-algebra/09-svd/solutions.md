# LA-9 Solutions

Signs of singular-vector pairs are arbitrary; if your $(\mathbf{u}_i, \mathbf{v}_i)$ are the negatives of the ones below, you are right.

## A. Warm-up

**A1.** (a) $A = \operatorname{diag}(3, -2)$. $A^\top A = \operatorname{diag}(9, 4)$, so $\sigma_1 = 3$, $\sigma_2 = 2$, $V = I$. Then $\mathbf{u}_1 = A\mathbf{e}_1/3 = \mathbf{e}_1$ and $\mathbf{u}_2 = A\mathbf{e}_2/2 = -\mathbf{e}_2$: $U = \operatorname{diag}(1, -1)$, $\Sigma = \operatorname{diag}(3, 2)$, $V^\top = I$. Check: $\operatorname{diag}(1,-1)\operatorname{diag}(3,2) = \operatorname{diag}(3,-2)$ ✓. The negative eigenvalue became a positive singular value with its sign absorbed into $U$ (§3.12).

(b) $A = \begin{pmatrix}0 & 2\\1 & 0\end{pmatrix}$: $A\mathbf{e}_1 = \mathbf{e}_2$ (length 1), $A\mathbf{e}_2 = 2\mathbf{e}_1$ (length 2). So $A$ maps the orthonormal basis $(\mathbf{e}_2, \mathbf{e}_1)$ to the orthogonal set $(2\mathbf{e}_1, \mathbf{e}_2)$, which is (3.1$'$) with $\mathbf{v}_1 = \mathbf{e}_2, \sigma_1 = 2, \mathbf{u}_1 = \mathbf{e}_1$ and $\mathbf{v}_2 = \mathbf{e}_1, \sigma_2 = 1, \mathbf{u}_2 = \mathbf{e}_2$. $U = I$, $\Sigma = \operatorname{diag}(2, 1)$, $V = V^\top = \begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}$. Check: $I\operatorname{diag}(2,1)\begin{pmatrix}0&1\\1&0\end{pmatrix} = \begin{pmatrix}0 & 2\\1 & 0\end{pmatrix}$ ✓.

(c) $2I_3$: $A^\top A = 4I$, every unit vector is an eigenvector, $\sigma_1 = \sigma_2 = \sigma_3 = 2$. Choice 1: $U = V = I$. Choice 2: $U = V = Q$ for any orthogonal $Q$, e.g. the permutation $\begin{pmatrix}0&1&0\\1&0&0\\0&0&1\end{pmatrix}$ or a rotation by $30°$ about the $z$-axis: $Q(2I)Q^\top = 2QQ^\top = 2I$ ✓. §3.3 permits this because the singular value $2$ is repeated three times, so its "eigenspace" of $A^\top A$ is all of $\mathbb{R}^3$ and any orthonormal basis of it is a valid $V$; $U$ then follows from $\mathbf{u}_i = A\mathbf{v}_i/2 = \mathbf{v}_i$. For $I_2$: $I = Q_\theta I Q_\theta^\top$ for every rotation $Q_\theta = \begin{pmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{pmatrix}$, one SVD per $\theta \in [0, 2\pi)$: infinitely many.

**A2.** (a) From Example 4.1, $\sigma_1 = 3\sqrt5$, $\mathbf{u}_1 = (1,3)^\top/\sqrt{10}$, $\mathbf{v}_1 = (1,1)^\top/\sqrt2$. $A_1 = 3\sqrt5\cdot\frac{1}{\sqrt{20}}\begin{pmatrix}1\\3\end{pmatrix}\begin{pmatrix}1 & 1\end{pmatrix} = \frac32\begin{pmatrix}1 & 1\\3 & 3\end{pmatrix} = \begin{pmatrix}1.5 & 1.5\\4.5 & 4.5\end{pmatrix}$. Then $A - A_1 = \begin{pmatrix}1.5 & -1.5\\-0.5 & 0.5\end{pmatrix}$, which equals $\sigma_2\mathbf{u}_2\mathbf{v}_2^\top = \sqrt5\cdot\frac{1}{\sqrt{20}}\begin{pmatrix}3\\-1\end{pmatrix}\begin{pmatrix}1 & -1\end{pmatrix} = \frac12\begin{pmatrix}3 & -3\\-1 & 1\end{pmatrix}$ ✓. $\|A - A_1\|_F^2 = 2.25 + 2.25 + 0.25 + 0.25 = 5$, so $\|A - A_1\|_F = \sqrt5 = \sigma_2$ ✓. $A - A_1$ is rank one with the single singular value $\sqrt5$, so $\|A - A_1\|_2 = \sqrt5 = \sigma_2$ ✓.

(b) $A^{-1} = V\Sigma^{-1}U^\top$. First $\Sigma^{-1}U^\top = \operatorname{diag}\!\left(\tfrac{1}{3\sqrt5}, \tfrac{1}{\sqrt5}\right)\frac{1}{\sqrt{10}}\begin{pmatrix}1 & 3\\3 & -1\end{pmatrix} = \frac{1}{\sqrt{50}}\begin{pmatrix}1/3 & 1\\3 & -1\end{pmatrix}$. Then $V(\cdot) = \frac{1}{\sqrt2}\begin{pmatrix}1 & 1\\1 & -1\end{pmatrix}\frac{1}{\sqrt{50}}\begin{pmatrix}1/3 & 1\\3 & -1\end{pmatrix} = \frac{1}{10}\begin{pmatrix}1/3 + 3 & 1 - 1\\1/3 - 3 & 1 + 1\end{pmatrix} = \begin{pmatrix}1/3 & 0\\-4/15 & 1/5\end{pmatrix}$. Formula: $A^{-1} = \frac{1}{15}\begin{pmatrix}5 & 0\\-4 & 3\end{pmatrix} = \begin{pmatrix}1/3 & 0\\-4/15 & 1/5\end{pmatrix}$ ✓.

**A3.** (a) $\|A\|_2 = \sigma_1 = \sqrt{45} = 3\sqrt5 \approx 6.708$; $\|A\|_F = \sqrt{45 + 5} = \sqrt{50} = 5\sqrt2 \approx 7.071$; $\kappa_2 = \sqrt{45}/\sqrt5 = 3$.
(b) $\operatorname{diag}(1, 10^{-3})$: $\kappa_2 = 10^3$. $1000\operatorname{diag}(1, 10^{-3}) = \operatorname{diag}(1000, 1)$: $\kappa_2 = 10^3$ (scaling cancels in the ratio). $\operatorname{diag}(10^{-8}, 10^{-8}) = 10^{-8}I$: $\kappa_2 = 1$.
(c) Digits lost $\approx \log_{10}\kappa_2$: about $3$, $3$, and $0$. The smallest determinant is $10^{-16}$ for $10^{-8}I$, which is the *best* conditioned of the three. The determinant measures volume scaling, not the shape of the ellipse, and says nothing about digits lost.

**A4.** (a) Full: $U$ is $5 \times 5$, $\Sigma$ is $5 \times 3$, $V$ is $3 \times 3$. Economy ($p = \min(5,3) = 3$): $U$ is $5 \times 3$, $\Sigma$ is $3 \times 3$, $V$ is $3 \times 3$. Reduced ($r = 2$): $U_r$ is $5 \times 2$, $\Sigma_r$ is $2 \times 2$, $V_r$ is $3 \times 2$.
(b) $C(A) = \operatorname{span}\{\mathbf{u}_1, \mathbf{u}_2\}$, $N(A^\top) = \operatorname{span}\{\mathbf{u}_3, \mathbf{u}_4, \mathbf{u}_5\}$, $C(A^\top) = \operatorname{span}\{\mathbf{v}_1, \mathbf{v}_2\}$, $N(A) = \operatorname{span}\{\mathbf{v}_3\}$ (Theorem 3.3).
(c) `U` is $5 \times 5$, `s` has length $3$, `Vt` is $3 \times 3$. One entry, `s[2]`, is zero up to rounding (about $10^{-16}$ times `s[0]`). `np.diag(s)` is $3 \times 3$, and `U @ np.diag(s)` tries to multiply $5 \times 5$ by $3 \times 3$: shape error. With `full_matrices=False`, `U` is $5 \times 3$ and $(5 \times 3)(3 \times 3)(3 \times 3)$ is fine. To reconstruct from the full form build `Sigma = np.zeros((5, 3)); Sigma[:3, :3] = np.diag(s)`.

**A5.** (a) $A^\top A = \sigma^2\mathbf{v}\mathbf{u}^\top\mathbf{u}\mathbf{v}^\top = \sigma^2\mathbf{v}\mathbf{v}^\top$ has the single nonzero eigenvalue $\sigma^2$ with eigenvector $\mathbf{v}$, so the only singular value is $\sigma$, with $\mathbf{v}_1 = \mathbf{v}$ and $\mathbf{u}_1 = A\mathbf{v}/\sigma = \mathbf{u}$. By (3.11), $A^+ = \sigma^{-1}\mathbf{v}\mathbf{u}^\top$. Since $A^\top = \sigma\mathbf{v}\mathbf{u}^\top$ and $\|A\|_F^2 = \sigma^2$ (Theorem 3.6), $A^\top/\|A\|_F^2 = \sigma^{-1}\mathbf{v}\mathbf{u}^\top = A^+$ ✓.
(b) $A = \begin{pmatrix}1\\2\end{pmatrix}\begin{pmatrix}1 & 2\end{pmatrix} = 5\,\mathbf{u}\mathbf{u}^\top$ with $\mathbf{u} = (1,2)^\top/\sqrt5$: $\sigma = 5$. $A^+ = A^\top/25 = \frac{1}{25}\begin{pmatrix}1 & 2\\2 & 4\end{pmatrix}$. $AA^+ = \frac{1}{25}\begin{pmatrix}1&2\\2&4\end{pmatrix}^2 = \frac{1}{25}\begin{pmatrix}5 & 10\\10 & 20\end{pmatrix} = \frac15\begin{pmatrix}1&2\\2&4\end{pmatrix} = \mathbf{u}\mathbf{u}^\top$, the orthogonal projector onto the line through $(1, 2)$ (Theorem 3.11(d)). Then $AA^+A = \mathbf{u}\mathbf{u}^\top\cdot5\mathbf{u}\mathbf{u}^\top = 5\mathbf{u}\mathbf{u}^\top = A$ ✓.

**A6.** (a) $A^\top = V\Sigma^\top U^\top$ is an SVD with the same diagonal: same $\sigma_i$. (b) $cA = U(c\Sigma)V^\top$; if $c < 0$ write it as $(-U)(|c|\Sigma)V^\top$: singular values $|c|\sigma_i$. (c) $QAW = (QU)\Sigma(W^\top V)^\top$ with $QU$ and $W^\top V$ orthogonal: same $\sigma_i$. (d) $A^{-1} = V\Sigma^{-1}U^\top$: singular values $1/\sigma_i$ (in the reversed order $1/\sigma_n \ge \cdots \ge 1/\sigma_1$). (e) $A^\top A = V(\Sigma^\top\Sigma)V^\top$ is symmetric PSD with eigenvalues $\sigma_i^2$, and for a PSD matrix the singular values are the eigenvalues (§3.12): $\sigma_i^2$. (f) $\begin{pmatrix}A\\0\end{pmatrix} = \begin{pmatrix}U & 0\\0 & I\end{pmatrix}\begin{pmatrix}\Sigma\\0\end{pmatrix}V^\top$; the block matrix is orthogonal and the stacked $\Sigma$ is diagonal with the same entries: same $\sigma_i$.

## B. Core

**B1.** (a) $a_{ij} = \mathbf{e}_i^\top A\mathbf{e}_j$, so by Cauchy–Schwarz and the definition of the spectral norm, $|a_{ij}| \le \|\mathbf{e}_i\|\,\|A\mathbf{e}_j\| \le \|A\|_2\|\mathbf{e}_j\| = \sigma_1$. ∎

(b) For any $\mathbf{x}$: $\|AB\mathbf{x}\| \le \|A\|_2\|B\mathbf{x}\| \le \|A\|_2\|B\|_2\|\mathbf{x}\|$ (definition of the spectral norm, twice). Taking the maximum over unit $\mathbf{x}$ gives $\|AB\|_2 \le \|A\|_2\|B\|_2$. For Frobenius, with $\mathbf{b}_j$ the columns of $B$: $\|AB\|_F^2 = \sum_j\|A\mathbf{b}_j\|^2 \le \sum_j\|A\|_2^2\|\mathbf{b}_j\|^2 = \|A\|_2^2\|B\|_F^2$. ∎

(c) $\kappa_2(AB) = \|AB\|_2\|(AB)^{-1}\|_2 = \|AB\|_2\|B^{-1}A^{-1}\|_2 \le \|A\|_2\|B\|_2\,\|B^{-1}\|_2\|A^{-1}\|_2 = \kappa_2(A)\kappa_2(B)$ by (b). ∎

(d) $A = Q\Lambda Q^\top$ symmetric. By §3.12 the singular values are $|\lambda_i|$, so $\|A\|_2 = \sigma_1 = \max_i|\lambda_i|$. Nonsymmetric example: $A = \begin{pmatrix}0 & 1\\0 & 0\end{pmatrix}$ has both eigenvalues $0$ but $\|A\|_2 = 1$ ($A^\top A = \operatorname{diag}(0, 1)$). Or the shear $\begin{pmatrix}1 & 1\\0 & 1\end{pmatrix}$: eigenvalues $1, 1$, but $\|A\|_2 = \phi \approx 1.618$ (Example 4.5).

**B2.** $A = \begin{pmatrix}1 & 0 & 1\\0 & 1 & 1\end{pmatrix}$.

*Gram matrix.* $AA^\top = \begin{pmatrix}1 + 0 + 1 & 0 + 0 + 1\\0 + 0 + 1 & 0 + 1 + 1\end{pmatrix} = \begin{pmatrix}2 & 1\\1 & 2\end{pmatrix}$. Trace $4$, determinant $3$: eigenvalues $3$ and $1$. So $\sigma_1 = \sqrt3$, $\sigma_2 = 1$. Check (3.8): $\|A\|_F^2 = 1 + 1 + 1 + 1 = 4 = 3 + 1$ ✓.

*Left singular vectors* (eigenvectors of $AA^\top$): for $3$: $\begin{pmatrix}-1 & 1\\1 & -1\end{pmatrix}\mathbf{u} = 0$, $\mathbf{u}_1 = (1, 1)^\top/\sqrt2$. For $1$: $\mathbf{u}_2 = (1, -1)^\top/\sqrt2$.

*Right singular vectors* via $\mathbf{v}_i = A^\top\mathbf{u}_i/\sigma_i$ (the §3.3 recipe with the roles of $U$ and $V$ swapped): $A^\top = \begin{pmatrix}1 & 0\\0 & 1\\1 & 1\end{pmatrix}$. $A^\top\mathbf{u}_1 = (1, 1, 2)^\top/\sqrt2$, norm $\sqrt{6/2} = \sqrt3 = \sigma_1$ ✓, so $\mathbf{v}_1 = (1, 1, 2)^\top/\sqrt6$. $A^\top\mathbf{u}_2 = (1, -1, 0)^\top/\sqrt2$, norm $1 = \sigma_2$ ✓, so $\mathbf{v}_2 = (1, -1, 0)^\top/\sqrt2$.

*Null space.* $A\mathbf{x} = 0$: $x_1 + x_3 = 0$, $x_2 + x_3 = 0$, so $\mathbf{x} \propto (1, 1, -1)$: $\mathbf{v}_3 = (1, 1, -1)^\top/\sqrt3$. Orthogonality checks: $\mathbf{v}_1\cdot\mathbf{v}_3 = (1 + 1 - 2)/\sqrt{18} = 0$, $\mathbf{v}_2\cdot\mathbf{v}_3 = (1 - 1 + 0) = 0$, $\mathbf{v}_1\cdot\mathbf{v}_2 = (1 - 1 + 0) = 0$ ✓. $A\mathbf{v}_3 = (1 - 1, 1 - 1)^\top/\sqrt3 = \mathbf{0}$ ✓.

*Full SVD.*
$$
A = \underbrace{\frac{1}{\sqrt2}\begin{pmatrix}1 & 1\\1 & -1\end{pmatrix}}_{U\ (2\times2)}\underbrace{\begin{pmatrix}\sqrt3 & 0 & 0\\0 & 1 & 0\end{pmatrix}}_{\Sigma\ (2\times3)}\underbrace{\begin{pmatrix}1/\sqrt6 & 1/\sqrt6 & 2/\sqrt6\\1/\sqrt2 & -1/\sqrt2 & 0\\1/\sqrt3 & 1/\sqrt3 & -1/\sqrt3\end{pmatrix}}_{V^\top\ (3\times3)}.
$$
Check: $U\Sigma = \frac{1}{\sqrt2}\begin{pmatrix}\sqrt3 & 1 & 0\\\sqrt3 & -1 & 0\end{pmatrix}$; times $V^\top$: row 1 $= \frac{1}{\sqrt2}\left(\sqrt3\cdot\frac{(1,1,2)}{\sqrt6} + \frac{(1,-1,0)}{\sqrt2}\right) = \frac12(1,1,2) + \frac12(1,-1,0) = (1, 0, 1)$ ✓; row 2 $= \frac12(1,1,2) - \frac12(1,-1,0) = (0, 1, 1)$ ✓.

*Four subspaces.* $r = 2 = m$. $C(A) = \operatorname{span}\{\mathbf{u}_1, \mathbf{u}_2\} = \mathbb{R}^2$. $N(A^\top) = \{\mathbf{0}\}$. $C(A^\top) = \operatorname{span}\{\mathbf{v}_1, \mathbf{v}_2\}$ = the plane $x + y - z = 0$ (the orthogonal complement of $\mathbf{v}_3$). $N(A) = \operatorname{span}\{(1, 1, -1)\}$.

*Rank-one approximation.* $A_1 = \sqrt3\cdot\frac{(1,1)^\top}{\sqrt2}\cdot\frac{(1,1,2)}{\sqrt6} = \frac{\sqrt3}{\sqrt{12}}\begin{pmatrix}1 & 1 & 2\\1 & 1 & 2\end{pmatrix} = \frac12\begin{pmatrix}1 & 1 & 2\\1 & 1 & 2\end{pmatrix}$. $A - A_1 = \frac12\begin{pmatrix}1 & -1 & 0\\-1 & 1 & 0\end{pmatrix} = \sigma_2\mathbf{u}_2\mathbf{v}_2^\top$ ✓. $\|A - A_1\|_F^2 = 4\cdot\frac14 = 1$, so $\|A - A_1\|_F = 1 = \sigma_2$ and $\|A - A_1\|_2 = \sigma_2 = 1$ (rank one, single singular value $1$).

**B3.** (a) The wide case through $AA^\top$ is the same derivation with the roles of $U$ and $V$ exchanged, which is the same as applying the tall-case routine to $A^\top$: if $A^\top = U'\Sigma V'^\top$ then $A = V'\Sigma U'^\top$.

```python
def svd_via_eig_any(A):
    """Economy SVD for any shape via the smaller Gram matrix (Exercise B3)."""
    A = np.asarray(A, float)
    m, n = A.shape
    if m >= n:
        return svd_via_eig(A)                       # eigh(A^T A), n x n
    # m < n: diagonalize A A^T (m x m), then v_i = A^T u_i / σ_i  (§3.3)
    lam, U = np.linalg.eigh(A @ A.T)
    order = np.argsort(lam)[::-1]; lam, U = lam[order], U[:, order]
    s = np.sqrt(np.clip(lam, 0, None))
    tol = max(m, n) * np.sqrt(EPS) * s[0]
    r = int(np.sum(s > tol))
    V_r = (A.T @ U[:, :r]) / s[:r]                  # n x r
    V = complete_orthonormal(V_r, n)[:, :m]         # n x m
    return U, s, V.T

for shape in [(3, 7), (2, 50)]:
    A = rng.standard_normal(shape)
    U, s, Vt = svd_via_eig_any(A)
    U2, s2, Vt2 = np.linalg.svd(A, full_matrices=False)
    assert np.allclose(s, s2) and np.allclose(U @ np.diag(s) @ Vt, A)
    assert np.allclose(U @ U.T, U2 @ U2.T) and np.allclose(Vt.T @ Vt, Vt2.T @ Vt2)
```
(Equivalently: `U, s, Vt = svd_via_eig(A.T); return Vt.T, s, U.T`.)

(b) Code:
```python
sig = 10.0 ** -np.arange(0, 13, 2)                  # 1, 1e-2, ..., 1e-12
A = random_orthogonal(7) @ np.diag(sig) @ random_orthogonal(7).T
_, s_e, _ = svd_via_eig(A, tol=0.0)
s_l = np.linalg.svd(A, compute_uv=False)
re_e, re_l = np.abs(s_e - sig) / sig, np.abs(s_l - sig) / sig
plt.loglog(sig, re_e, "o-", label="via eigh(A^T A)"); plt.loglog(sig, re_l, "s-", label="np.linalg.svd")
plt.loglog(sig, EPS / (2 * sig**2), "k--", label="ε/(2σ²)"); plt.loglog(sig, EPS / sig, "k:", label="ε/σ")
plt.xlabel("true σ_i"); plt.ylabel("relative error"); plt.legend()
```
Typical output (one run): relative errors via `eigh` of $10^{-16}, 2\times10^{-14}, 8\times10^{-10}, 2\times10^{-6}, 0.26, 36, 1$ for $\sigma_i = 1, 10^{-2}, \dots, 10^{-12}$; via `np.linalg.svd`: $10^{-16}, 4\times10^{-15}, 4\times10^{-14}, 8\times10^{-12}, 2\times10^{-9}, 3\times10^{-7}, 1.4\times10^{-5}$.

Explanation. The computed eigenvalue is $\hat\lambda_i = \sigma_i^2 + \delta$ with $|\delta| \approx \varepsilon\|A^\top A\|_2 = \varepsilon\sigma_1^2$. Then $\hat\sigma_i = \sqrt{\sigma_i^2 + \delta} \approx \sigma_i + \delta/(2\sigma_i)$ when $|\delta| \ll \sigma_i^2$, so the relative error is $|\delta|/(2\sigma_i^2) \approx \varepsilon\sigma_1^2/(2\sigma_i^2)$: on log–log axes against $\sigma_i$ that is a line of **slope $-2$**. Once $\sigma_i^2 < \varepsilon\sigma_1^2$, i.e. $\sigma_i < \sqrt\varepsilon\,\sigma_1 \approx 10^{-8}$, the eigenvalue is swamped by noise: $\hat\lambda_i$ is a random number of size $\varepsilon\sigma_1^2$ (possibly negative, clipped to $0$) and $\hat\sigma_i \approx 10^{-8}$ regardless of the true value, so the relative error is $O(1)$ or worse (the $36$ above is $\hat\sigma = 3.7\times10^{-9}$ reported for a true $10^{-10}$). Golub–Kahan works on $A$ itself and gets absolute error $\approx\varepsilon\sigma_1$ in $\sigma_i$: relative error $\varepsilon\sigma_1/\sigma_i$, **slope $-1$**, correct to about $16 - \log_{10}(\sigma_1/\sigma_i)$ digits all the way down to $\sigma_i \approx 10^{-16}$. Squaring the matrix squares the condition number (B7(b)); that is the whole story.

**B4.** (a) $B_t = \operatorname{diag}(2 - t, 0)$ has one nonzero entry $2 - t \ge 1$ for $t \in [0, 1]$, so rank $1$. $A - B_t = \operatorname{diag}(t, 1)$ has singular values $t$ and $1$, so $\|A - B_t\|_2 = \max(t, 1) = 1 = \sigma_2$. Every $B_t$ with $0 \le t \le 1$ is a best rank-1 approximation in the spectral norm; $A_1 = B_0 = \operatorname{diag}(2, 0)$ is only one of them.

(b) $\|A - B_t\|_F^2 = t^2 + 1$, minimized only at $t = 0$, where it equals $\sigma_2^2 = 1$. In the Frobenius norm the minimizer is unique.

(c) Suppose $\operatorname{rank}B \le k$ and $\|A - B\|_F^2 = \sum_{i > k}\sigma_i^2$. Then every inequality in the proof of Theorem 3.7 is an equality.

*Step 3 equality* says $\|P(A - B)\|_F^2 = 0$, i.e. $PA = PB = B$.

*Step 9 equality* says $\|PA\|_F^2 = \sum_{i \le k}\sigma_i^2$, so the two inequalities of step 8 are equalities. The first was $\sum_{i > k}\sigma_i^2c_i \le \sigma_k^2\sum_{i > k}c_i$; since $\sigma_i \le \sigma_{k+1} < \sigma_k$ for every $i > k$, each term satisfies $\sigma_i^2c_i \le \sigma_{k+1}^2c_i$, which is strictly less than $\sigma_k^2c_i$ unless $c_i = 0$. So $c_i = 0$ for all $i > k$. The last was $\sum_{i \le k}[\sigma_i^2c_i + \sigma_i^2(1 - c_i)] = \sum_{i \le k}\sigma_i^2$ reached from $\sum_{i\le k}[\sigma_i^2 c_i + \sigma_k^2(1-c_i)]$; with the total now $\sum_{i \le k}\sigma_i^2c_i = \sum_{i\le k}\sigma_i^2$ and each $0 \le c_i \le 1$, $\sigma_i > 0$, we must have $c_i = 1$ for all $i \le k$.

*Geometry of $c_i = 1$.* $c_i = \|P\mathbf{u}_i\|^2 = 1 = \|\mathbf{u}_i\|^2$, and $\|\mathbf{u}_i\|^2 = \|P\mathbf{u}_i\|^2 + \|(I - P)\mathbf{u}_i\|^2$, so $(I - P)\mathbf{u}_i = 0$: $\mathbf{u}_i \in C(P)$ for $i = 1, \dots, k$. Thus $\operatorname{span}\{\mathbf{u}_1, \dots, \mathbf{u}_k\} \subseteq C(P)$; but $\dim C(P) = p \le k$, so $C(P) = \operatorname{span}\{\mathbf{u}_1, \dots, \mathbf{u}_k\}$ and $P = U_kU_k^\top$.

*Conclusion.* $B = PA = U_kU_k^\top\sum_i\sigma_i\mathbf{u}_i\mathbf{v}_i^\top = \sum_{i \le k}\sigma_i\mathbf{u}_i\mathbf{v}_i^\top = A_k$, because $U_kU_k^\top\mathbf{u}_i = \mathbf{u}_i$ for $i \le k$ and $= \mathbf{0}$ for $i > k$. ∎ (In (a), $\sigma_1 = 2 > 1 = \sigma_2$ and the Frobenius minimizer was indeed unique, while the spectral norm never had uniqueness.)

**B5.** Write $A = U\Sigma V^\top$, $A^+ = V\Sigma^+U^\top$. Note $\Sigma\Sigma^+$ is the $m \times m$ diagonal matrix $D_m$ with $r$ ones then zeros, and $\Sigma^+\Sigma$ is the $n \times n$ analogue $D_n$; also $D_m\Sigma = \Sigma = \Sigma D_n$ and $D_n\Sigma^+ = \Sigma^+ = \Sigma^+D_m$ (the ones sit exactly where the nonzero $\sigma_i$ are).

(a) $AA^+A = U\Sigma V^\top V\Sigma^+U^\top U\Sigma V^\top = U(\Sigma\Sigma^+\Sigma)V^\top = U(D_m\Sigma)V^\top = U\Sigma V^\top = A$. $A^+AA^+ = V(\Sigma^+\Sigma\Sigma^+)U^\top = V(D_n\Sigma^+)U^\top = A^+$. $AA^+ = UD_mU^\top$ and $(UD_mU^\top)^\top = UD_m^\top U^\top = UD_mU^\top$ since $D_m$ is diagonal. Likewise $A^+A = VD_nV^\top$ is symmetric. ∎

(b) $A^+ = V\Sigma^+U^\top = \sum_{i \le r}\sigma_i^{-1}\mathbf{v}_i\mathbf{u}_i^\top$ is itself in SVD form: nonzero singular values $\sigma_i^{-1}$, left vectors $\mathbf{v}_i$, right vectors $\mathbf{u}_i$ (reorder to make them decreasing; the formula (3.11) does not care about order). Applying (3.11) to it inverts each $\sigma_i^{-1}$ and swaps the vectors back: $(A^+)^+ = \sum_{i \le r}\sigma_i\mathbf{u}_i\mathbf{v}_i^\top = A$. For the transpose, $A^\top = V\Sigma^\top U^\top$ is an SVD of $A^\top$, so $(A^\top)^+ = U(\Sigma^\top)^+V^\top$. Both $(\Sigma^\top)^+$ and $(\Sigma^+)^\top$ are the $m \times n$ matrix with $1/\sigma_i$ in the first $r$ diagonal positions, so $(A^\top)^+ = U(\Sigma^+)^\top V^\top = (V\Sigma^+U^\top)^\top = (A^+)^\top$. ∎

(c) $A = (1\ \ 1)$, $B = (1\ \ 0)^\top$. $AB = (1)$, so $(AB)^+ = (1)$. By A5, $B^+ = B^\top/\|B\|^2 = (1\ \ 0)$ and $A^+ = A^\top/\|A\|^2 = (1/2,\ 1/2)^\top$. Then $B^+A^+ = (1\ \ 0)(1/2,\ 1/2)^\top = 1/2 \ne 1$.

(d) Let $\mathbf{x}$ minimize $f(\mathbf{x}) = \|A\mathbf{x} - \mathbf{b}\|^2 = \mathbf{x}^\top A^\top A\mathbf{x} - 2\mathbf{b}^\top A\mathbf{x} + \mathbf{b}^\top\mathbf{b}$. The gradient $2A^\top A\mathbf{x} - 2A^\top\mathbf{b}$ must vanish: $A^\top A\mathbf{x} = A^\top\mathbf{b}$ (LA-5; equivalently, the residual $A\mathbf{x} - \mathbf{b}$ must be orthogonal to $C(A)$, i.e. $A^\top(A\mathbf{x} - \mathbf{b}) = 0$). In the singular bases, with $\mathbf{x} = \sum_jc_j\mathbf{v}_j$ and $\mathbf{b} = \sum_id_i\mathbf{u}_i$: $A^\top A\mathbf{x} = \sum_{j \le r}\sigma_j^2c_j\mathbf{v}_j$ and $A^\top\mathbf{b} = \sum_{j \le r}\sigma_jd_j\mathbf{v}_j$. Comparing coefficients: $\sigma_j^2c_j = \sigma_jd_j$, i.e. $c_j = d_j/\sigma_j$ for $j \le r$, and no condition on $c_{r+1}, \dots, c_n$. So the solution set of the normal equations is $\mathbf{x}^+ + N(A)$, exactly the minimizer set of Theorem 3.10 (as it must be, since both describe the same minimizers). The unique one in $C(A^\top) = \operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_r\}$ has $c_j = 0$ for $j > r$: it is $\mathbf{x}^+ = \sum_{j \le r}(d_j/\sigma_j)\mathbf{v}_j = A^+\mathbf{b}$. ∎

**B6.** (a) Since $P$ is symmetric, row $i$ of $XP$ is $\mathbf{x}_i^\top P = (P\mathbf{x}_i)^\top$, so row $i$ of $X - XP$ is $(\mathbf{x}_i - P\mathbf{x}_i)^\top$ and $\|X - XP\|_F^2 = \sum_i\|\mathbf{x}_i - P\mathbf{x}_i\|^2$. Also $\operatorname{rank}(XP) \le \operatorname{rank}P = k$, because every row of $XP$ is $(P\mathbf{x}_i)^\top$ and lies in the $k$-dimensional $C(P)$. By Theorem 3.7, $\|X - XP\|_F^2 \ge \sum_{i > k}\sigma_i^2$ for every such $P$. For $P = V_kV_k^\top$: $XV_kV_k^\top = \sum_i\sigma_i\mathbf{u}_i\mathbf{v}_i^\top V_kV_k^\top$, and $\mathbf{v}_i^\top V_kV_k^\top = \mathbf{v}_i^\top$ for $i \le k$ (since $\mathbf{v}_i \in C(V_k)$) and $= \mathbf{0}$ for $i > k$ (since $\mathbf{v}_i \perp C(V_k)$). So $XV_kV_k^\top = \sum_{i \le k}\sigma_i\mathbf{u}_i\mathbf{v}_i^\top = X_k$ and $\|X - X_k\|_F^2 = \sum_{i > k}\sigma_i^2$: the bound is attained. The best-fitting $k$-plane through the centroid is spanned by $\mathbf{v}_1, \dots, \mathbf{v}_k$. ∎

(b) Let $Q$ ($p \times k$) have orthonormal columns spanning $C(P)$, so $P = QQ^\top$. The data projected into that subspace, in its own coordinates, is $XQ$ ($n \times k$, still centered), with covariance $Q^\top SQ$ and total variance $\operatorname{tr}(Q^\top SQ) = \operatorname{tr}(Q^\top X^\top XQ)/(n-1) = \|XQ\|_F^2/(n-1)$. And $\|XQ\|_F^2 = \operatorname{tr}(Q^\top X^\top XQ) = \operatorname{tr}(X^\top XQQ^\top) = \operatorname{tr}(X^\top XP) = \operatorname{tr}(P X^\top X P) = \|XP\|_F^2$, using the trace identity and $P^2 = P = P^\top$. So $\sum_i\|P\mathbf{x}_i\|^2 = \|XP\|_F^2 = (n-1)\times(\text{captured variance})$. Pythagoras for each row (the projection $P\mathbf{x}_i$ is orthogonal to the residual $\mathbf{x}_i - P\mathbf{x}_i$): $\|\mathbf{x}_i\|^2 = \|P\mathbf{x}_i\|^2 + \|\mathbf{x}_i - P\mathbf{x}_i\|^2$; summing, $\|X\|_F^2 = \|XP\|_F^2 + \|X - XP\|_F^2$. The left side does not depend on $P$, so maximizing the captured variance $\|XP\|_F^2$ and minimizing the residual $\|X - XP\|_F^2$ are the same problem with the same solution $P = V_kV_k^\top$, and the maximum captured variance is $\sum_{i \le k}\sigma_i^2/(n-1) = \lambda_1 + \cdots + \lambda_k$. ∎

(c) Points $(0, 10), (1, 10), (2, 10)$; mean $(1, 10)$; centered rows $(-1, 0), (0, 0), (1, 0)$. $X^\top X = \begin{pmatrix}2 & 0\\0 & 0\end{pmatrix}$: first principal component $(1, 0)$ (horizontal), variance $2/2 = 1$; second component has variance $0$. Uncentered: $\tilde X^\top\tilde X = \begin{pmatrix}0 + 1 + 4 & 0 + 10 + 20\\30 & 300\end{pmatrix} = \begin{pmatrix}5 & 30\\30 & 300\end{pmatrix}$, trace $305$, determinant $600$, eigenvalues $\frac{305 \pm \sqrt{305^2 - 2400}}{2} = \frac{305 \pm 301.04}{2} = 303.02, 1.98$. Top eigenvector: $(5 - 303.02)w_1 + 30w_2 = 0$, so $w_2 = 9.93\,w_1$: direction $\approx(0.100, 0.995)$, at $84.3°$, essentially the direction of the mean $(1, 10)$ (which is at $\arctan 10 = 84.3°$) and nearly **perpendicular** to the direction in which the data actually varies. The uncentered "first component" measures where the cloud sits relative to the origin, not how it is spread.

**B7.** (a) $\kappa_2 = \sigma_1/\sigma_n \ge 1$ because $\sigma_1 \ge \sigma_n$. Equality iff $\sigma_1 = \sigma_n$ iff all singular values equal some $c > 0$ iff $\Sigma = cI$ iff $A = U(cI)V^\top = c\,UV^\top = cQ$ with $Q = UV^\top$ orthogonal. Conversely, if $A = cQ$ then $A^\top A = c^2I$, all $\sigma_i = |c|$, $\kappa_2 = 1$. ∎

(b) The singular values of $A^\top A$ are $\sigma_i^2$ (A6(e)), so $\kappa_2(A^\top A) = \sigma_1^2/\sigma_n^2 = \kappa_2(A)^2$. ∎ Consequence: solving the normal equations $A^\top A\mathbf{x} = A^\top\mathbf{b}$ exposes you to $\kappa_2(A)^2$, so you lose about $2\log_{10}\kappa_2(A)$ digits instead of $\log_{10}\kappa_2(A)$. With $\kappa_2(A) = 10^8$ (quite ordinary for polynomial-fitting design matrices), QR or SVD keeps about 8 digits while the normal equations keep none. This is why `np.linalg.lstsq` uses the SVD and ST-3 will warn you off the normal equations.

(c) $F = \begin{pmatrix}1 & 1\\0 & 1\end{pmatrix}$, $F^{-1} = \begin{pmatrix}1 & -1\\0 & 1\end{pmatrix}$. From Example 4.5, $\sigma_1 = \phi$, $\sigma_2 = 1/\phi$, so $\kappa_2 = \phi^2 = \frac{3 + \sqrt5}{2} \approx 2.618$. Column sums of $|F|$: $1, 2$; of $|F^{-1}|$: $1, 2$; so $\kappa_1 = 2\cdot2 = 4$. Row sums of $|F|$: $2, 1$; of $|F^{-1}|$: $2, 1$; so $\kappa_\infty = 2\cdot2 = 4$. Three norms, two different condition numbers ($2.618$ vs $4$); all within the factor $n = 2$ of each other that the general theory guarantees.

(d) $cA$ has singular values $|c|\sigma_i$ (A6(b)), and the ratio $\sigma_1/\sigma_n$ is unchanged. $\det = 1$ with $\kappa_2 \ge 10^4$: $\operatorname{diag}(100, 0.01)$ ($\kappa_2 = 10^4$ exactly), or the shear $\begin{pmatrix}1 & 100\\0 & 1\end{pmatrix}$: $\sigma_1\sigma_2 = |\det| = 1$ and $\sigma_1^2 + \sigma_2^2 = \|F\|_F^2 = 10002$, so $\sigma_1 \approx 100.01$, $\sigma_2 \approx 0.0099990$, $\kappa_2 \approx 1.0002\times10^4$. $\det = 10^{-20}$ with $\kappa_2 = 1$: $10^{-10}I_2$.

**B8.** (a) Add after the Frobenius loop:
```python
ks2 = np.arange(0, 256, 5)
spec = np.array([np.linalg.norm(img - rank_k(U, s, Vt, k), 2) for k in ks2])
assert np.allclose(spec, s[ks2]), "Theorem 3.8: ||A - A_k||_2 = σ_{k+1}"
ax0.semilogy(ks2[1:], spec[1:], "C0.", label="‖A − A_k‖_2 = σ_{k+1}")
```
(Each spectral norm costs an SVD, hence the stride of 5.) (b):
```python
tau = 0.05
noisy = img + tau * rng.standard_normal(img.shape)
Un, sn, Vtn = np.linalg.svd(noisy, full_matrices=False)
k_range = np.arange(1, 150)
err_clean = np.array([np.linalg.norm(rank_k(Un, sn, Vtn, k) - img) for k in k_range])
k_opt = k_range[np.argmin(err_clean)]
print(k_opt, err_clean.min(), np.linalg.norm(noisy - img))
```
With the seed in `project.py` the optimum lands at $k \approx 24$–$25$ with $\|\tilde A_k - A\|_F \approx 8.7$, against $\|\tilde A - A\|_F \approx 12.8$ for the untruncated noisy image ($= \tau\sqrt{N^2} = 0.05\cdot256$): truncation removed about a third of the noise energy's square root while keeping almost all of the signal.

Why the optimum is finite. Split $\tilde A_k - A = (\tilde A_k - \tilde A) + (\tilde A - A) = -(\text{discarded part of }\tilde A) + E$ where $E$ is the noise. Roughly, $\|\tilde A_k - A\|_F^2 \approx (\text{discarded signal}) + (\text{kept noise}) = \sum_{i > k}\sigma_i^2 + k\cdot\tau^2N$, because the noise matrix spreads its total energy $\tau^2N^2$ evenly over $N$ singular directions (each noise singular value is of order $\tau\sqrt N$, the flat floor in the scree plot), so keeping $k$ directions keeps about $k/N$ of the noise energy, while the signal's energy is concentrated in the first few directions. The first term decreases with $k$ and the second grows linearly; the minimum is where the marginal signal $\sigma_{k+1}^2$ falls to the marginal noise $\tau^2N = 0.64$, i.e. $\sigma_{k+1} \approx 0.8$, which the clean spectrum reaches around $k \approx 40$. The measured optimum is a bit smaller because the noise also perturbs the *directions* $\tilde{\mathbf{u}}_i, \tilde{\mathbf{v}}_i$ of the weak components, so those components carry less usable signal than their $\sigma_i$ suggests. Truncation is a denoiser exactly when the signal is low rank and the noise is not (§3.13).

## C. Challenge

**C1.** (a) $\|A - Q\|_F^2 = \operatorname{tr}\bigl((A - Q)^\top(A - Q)\bigr) = \operatorname{tr}(A^\top A) - 2\operatorname{tr}(Q^\top A) + \operatorname{tr}(Q^\top Q) = \|A\|_F^2 - 2\operatorname{tr}(Q^\top A) + n$. Only the middle term depends on $Q$, so minimizing the distance is maximizing $\operatorname{tr}(Q^\top A)$. With $A = U\Sigma V^\top$ and the trace identity, $\operatorname{tr}(Q^\top A) = \operatorname{tr}(Q^\top U\Sigma V^\top) = \operatorname{tr}(V^\top Q^\top U\,\Sigma) = \operatorname{tr}(Z\Sigma)$ where $Z = V^\top Q^\top U$ is a product of orthogonal matrices, hence orthogonal. Since $\Sigma$ is diagonal, $\operatorname{tr}(Z\Sigma) = \sum_iz_{ii}\sigma_i$. Each row of an orthogonal matrix is a unit vector, so $|z_{ii}| \le 1$, and therefore $\operatorname{tr}(Z\Sigma) \le \sum_i\sigma_i$, with equality iff $z_{ii} = 1$ for every $i$ with $\sigma_i > 0$. If $A$ is invertible, all $\sigma_i > 0$, so all $z_{ii} = 1$; a unit row with a $1$ in position $i$ has zeros elsewhere, so $Z = I$, i.e. $V^\top Q^\top U = I$, $Q^\top = VU^\top$, $Q = UV^\top$: the polar factor, and it is the unique maximizer. (If $A$ is singular, $Q = UV^\top$ still attains the maximum, with $Z = I$, but any $Z$ that is the identity on the first $r$ coordinates also does, so the nearest orthogonal matrix is not unique.) The minimum distance is $\|A\|_F^2 - 2\sum\sigma_i + n = \sum_i(\sigma_i - 1)^2$: the distance from $A$ to the orthogonal group is the distance from its singular values to $(1, \dots, 1)$. ∎

(b) $\|A - BQ\|_F^2 = \|A\|_F^2 - 2\operatorname{tr}(Q^\top B^\top A) + \|BQ\|_F^2$, and $\|BQ\|_F = \|B\|_F$ by orthogonal invariance. So we maximize $\operatorname{tr}(Q^\top M)$ with $M = B^\top A$, which by the argument of (a) applied to $M = U_M\Sigma_MV_M^\top$ gives $Q = U_MV_M^\top$, the polar factor of $B^\top A$. ∎

(c) In the "Your turn" cell:
```python
def procrustes(A, B):
    """Orthogonal Q minimizing ||A - B Q||_F: the polar factor of B^T A (Exercise C1)."""
    Um, _, Vmt = np.linalg.svd(B.T @ A)
    return Um @ Vmt

def rot(t):
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

B = rng.standard_normal((50, 2)); t0 = 0.7
A = B @ rot(t0) + 0.01 * rng.standard_normal(B.shape)
Q = procrustes(A, B)
assert np.allclose(Q.T @ Q, np.eye(2)) and np.allclose(Q, rot(t0), atol=0.01)
angles = np.linspace(-np.pi, np.pi, 7201)
brute = np.array([np.linalg.norm(A - B @ rot(t)) for t in angles])
assert abs(angles[brute.argmin()] - t0) < 0.002
assert np.linalg.norm(A - B @ Q) <= brute.min() + 1e-12
```
The brute force only searches rotations ($\det = +1$); the SVD answer may be a reflection when the data call for one, which is why the general statement in (a)–(b) is about orthogonal $Q$. (Forcing a rotation: flip the sign of the last column of $U_M$ if $\det(U_MV_M^\top) < 0$; this is the Kabsch algorithm used to align molecules and point clouds.)

**C2.** (a) *Max–min.* Let $S$ be any $i$-dimensional subspace and $T_i = \operatorname{span}\{\mathbf{q}_i, \dots, \mathbf{q}_n\}$, of dimension $n - i + 1$. Since $\dim(S \cap T_i) \ge \dim S + \dim T_i - n = 1$ (LA-3), there is a unit $\mathbf{w} \in S \cap T_i$; writing $\mathbf{w} = \sum_{j \ge i}c_j\mathbf{q}_j$ with $\sum c_j^2 = 1$ gives $\mathbf{w}^\top M\mathbf{w} = \sum_{j \ge i}\lambda_jc_j^2 \le \lambda_i$. Hence $\min_{\mathbf{w} \in S}\mathbf{w}^\top M\mathbf{w} \le \lambda_i$ for every $S$, so the max over $S$ is $\le \lambda_i$. For $S = \operatorname{span}\{\mathbf{q}_1, \dots, \mathbf{q}_i\}$, every unit $\mathbf{w} = \sum_{j \le i}c_j\mathbf{q}_j$ has $\mathbf{w}^\top M\mathbf{w} = \sum_{j \le i}\lambda_jc_j^2 \ge \lambda_i$, with equality at $\mathbf{q}_i$; so the min over this $S$ is exactly $\lambda_i$, and the max–min equals $\lambda_i$. *Min–max.* Let $T$ be any subspace of dimension $n - i + 1$ and $S_i = \operatorname{span}\{\mathbf{q}_1, \dots, \mathbf{q}_i\}$; $\dim(T \cap S_i) \ge 1$, and a unit $\mathbf{w}$ in it has $\mathbf{w}^\top M\mathbf{w} \ge \lambda_i$, so $\max_{\mathbf{w} \in T} \ge \lambda_i$ for every $T$; for $T = T_i$ the max is exactly $\lambda_i$ (attained at $\mathbf{q}_i$). ∎

(b) Apply (a) to the symmetric matrix $M = A^\top A$: $\mathbf{w}^\top M\mathbf{w} = \|A\mathbf{w}\|^2$, the eigenvalues are $\sigma_i^2$ (padded with zeros), and the eigenvectors are the $\mathbf{v}_i$. Taking square roots (a monotone operation that commutes with max and min) gives both formulas for $\sigma_i(A)$. ∎

(c) Let $T_A = \operatorname{span}\{\mathbf{v}_i(A), \dots, \mathbf{v}_n(A)\}$ and $T_B = \operatorname{span}\{\mathbf{v}_j(B), \dots, \mathbf{v}_n(B)\}$, of dimensions $n - i + 1$ and $n - j + 1$. For unit $\mathbf{w} \in T_A$, $\|A\mathbf{w}\|^2 = \sum_{l \ge i}\sigma_l(A)^2c_l^2 \le \sigma_i(A)^2$; the same expansion in the right singular vectors of $B$ gives $\|B\mathbf{w}\| \le \sigma_j(B)$ for unit $\mathbf{w} \in T_B$. Now $\dim(T_A \cap T_B) \ge (n - i + 1) + (n - j + 1) - n = n - (i + j - 1) + 1 \ge 1$. Choose a subspace $T \subseteq T_A \cap T_B$ of dimension exactly $n - (i + j - 1) + 1$. For every unit $\mathbf{w} \in T$, by the triangle inequality $\|(A + B)\mathbf{w}\| \le \|A\mathbf{w}\| + \|B\mathbf{w}\| \le \sigma_i(A) + \sigma_j(B)$. The min–max formula of (b) with index $i + j - 1$ says $\sigma_{i+j-1}(A + B)$ is the minimum over all subspaces $T'$ of that dimension of $\max_{\mathbf{w} \in T'}\|(A + B)\mathbf{w}\|$, which is at most the value for our $T$: $\sigma_{i+j-1}(A + B) \le \sigma_i(A) + \sigma_j(B)$. ∎

(d) If $\operatorname{rank}B \le k$ then $\sigma_{k+1}(B) = 0$. Apply Weyl to $(A - B) + B = A$ with $j = k + 1$: $\sigma_{i+k}(A) \le \sigma_i(A - B) + \sigma_{k+1}(B) = \sigma_i(A - B)$. Then $\|A - B\|_F^2 = \sum_i\sigma_i(A - B)^2 \ge \sum_i\sigma_{i+k}(A)^2 = \sum_{l > k}\sigma_l(A)^2$, which is Theorem 3.7; and $\|A - B\|_2 = \sigma_1(A - B) \ge \sigma_{k+1}(A)$, which is Theorem 3.8. Equality at $B = A_k$ in both as before. ∎ This proof is shorter than the one in §3.7 once Weyl's inequality is available, and Weyl's inequality is the tool for everything about how singular values move under perturbation (LA-13).

**C3.** (a) "$\le$": take $E = -\sigma_n\mathbf{u}_n\mathbf{v}_n^\top$. Then $\|E\|_2 = \sigma_n$ and $A + E = \sum_{i < n}\sigma_i\mathbf{u}_i\mathbf{v}_i^\top$ has rank $n - 1$: singular. "$\ge$": if $A + E$ is singular then $\operatorname{rank}(A + E) \le n - 1$, so by Theorem 3.8 with $k = n - 1$ and $B = A + E$, $\|A - B\|_2 = \|E\|_2 \ge \sigma_n$. Hence $\min\{\|E\|_2 : A + E \text{ singular}\} = \sigma_n(A)$, and dividing by $\|A\|_2 = \sigma_1$, the relative distance to the nearest singular matrix is $\sigma_n/\sigma_1 = 1/\kappa_2(A)$. ∎ A matrix with $\kappa_2 = 10^{12}$ is within relative distance $10^{-12}$ of a singular one; rounding its entries to 16 digits leaves it indistinguishable from one at 4 digits past that.

(b) Code:
```python
from scipy.linalg import hilbert
for n in range(2, 15):
    H = hilbert(n); kappa = np.linalg.cond(H)
    x = np.linalg.solve(H, H @ np.ones(n))
    digits = -np.log10(max(np.abs(x - 1).max(), 1e-17))
    print(f"n={n:2d}  log10 κ={np.log10(kappa):5.2f}  correct digits={digits:5.1f}  16−log10 κ={16-np.log10(kappa):5.1f}")
```
Output (this platform; your last digit may differ):

| $n$ | $\log_{10}\kappa_2$ | correct digits | $16 - \log_{10}\kappa_2$ |
|---|---|---|---|
| 2 | 1.29 | 15.2 | 14.7 |
| 4 | 4.19 | 13.3 | 11.8 |
| 6 | 7.17 | 10.2 | 8.8 |
| 8 | 10.18 | 6.9 | 5.8 |
| 10 | 13.20 | 4.6 | 2.8 |
| 11 | 14.72 | 2.6 | 1.3 |
| 12 | 16.21 | 0.6 | −0.2 |
| 13 | 17.67 | −1.7 | — |
| 14 | 17.40 | −0.7 | — |

The rule holds to within about a digit: correct digits $\approx 16.5 - \log_{10}\kappa_2$ down the whole table (the bound (3.10) is a worst case over $\delta\mathbf{b}$, and the actual rounding error is a little kinder). Once $\kappa_2 > 10^{16}$ ($n \ge 12$) the computed solution has no correct digits at all, with no exception raised and at most a `LinAlgWarning`; worse, the *computed* $\kappa_2$ stops growing ($17.67$, then $17.40$) because the computed $\sigma_n$ is itself just rounding noise of size $\varepsilon\sigma_1$, so you can no longer even measure how bad the problem is. A least-squares fit of $\log_{10}\kappa_2(H_n)$ against $n$ over $n = 2, \dots, 10$ gives a slope of $1.49$ decades per unit $n$, matching the classical asymptotic $\kappa_2(H_n) \sim (1 + \sqrt2)^{4n}/\sqrt n$, whose slope is $4\log_{10}(1 + \sqrt2) = 1.53$. Fitting a polynomial of degree $11$ to data by the monomial basis produces exactly this matrix as $A^\top A$ (B7(b)), which is why nobody does that.
