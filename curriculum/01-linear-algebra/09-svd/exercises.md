# LA-9 Exercises

Conventions as in the notes: $A$ is $m \times n$ of rank $r$ with SVD $A = U\Sigma V^\top$, singular values $\sigma_1 \ge \cdots \ge \sigma_r > 0$; $\|A\|_2 = \sigma_1$, $\|A\|_F^2 = \sum\sigma_i^2$, $\kappa_2 = \sigma_1/\sigma_n$; $A_k = \sum_{i \le k}\sigma_i\mathbf{u}_i\mathbf{v}_i^\top$; $A^+ = V\Sigma^+U^\top$. Signs of singular-vector pairs are arbitrary; any answer that flips a $(\mathbf{u}_i, \mathbf{v}_i)$ pair together is correct.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Write down an SVD (with $\sigma$'s in decreasing order) of each matrix by inspection, no $A^\top A$ needed:
(a) $A = \begin{pmatrix}3 & 0\\0 & -2\end{pmatrix}$; (b) $A = \begin{pmatrix}0 & 2\\1 & 0\end{pmatrix}$; (c) $A = 2I_3$.
For (c), exhibit two different valid choices of $(U, V)$, and explain what in §3.3 permits this. Then show $I_2$ has infinitely many SVDs.

**A2.** For $A = \begin{pmatrix}3 & 0\\4 & 5\end{pmatrix}$ (Example 4.1), (a) compute $A_1 = \sigma_1\mathbf{u}_1\mathbf{v}_1^\top$ and verify $\|A - A_1\|_F = \sigma_2$ and $\|A - A_1\|_2 = \sigma_2$; (b) compute $A^{-1} = V\Sigma^{-1}U^\top$ and check it against the $2 \times 2$ inverse formula.

**A3.** (a) For the matrix of A2 give $\|A\|_2$, $\|A\|_F$, $\kappa_2(A)$. (b) Give $\kappa_2$ for $\operatorname{diag}(1, 10^{-3})$, for $1000\operatorname{diag}(1, 10^{-3})$, and for $\operatorname{diag}(10^{-8}, 10^{-8})$. (c) For each matrix in (b), how many significant digits do you expect to lose solving $A\mathbf{x} = \mathbf{b}$ in double precision? Which of the three has the smallest determinant, and does that matter?

**A4.** $A$ is $5 \times 3$ of rank 2. (a) Give the shapes of $U, \Sigma, V$ in the full, economy, and reduced SVDs. (b) Which columns of $U$ and $V$ span $C(A)$, $N(A^\top)$, $C(A^\top)$, $N(A)$? (c) What shapes does `U, s, Vt = np.linalg.svd(A)` return, how many entries of `s` are (numerically) zero, and why does `U @ np.diag(s) @ Vt` raise a shape error while `np.linalg.svd(A, full_matrices=False)` does not?

**A5.** Let $A = \sigma\mathbf{u}\mathbf{v}^\top$ with $\|\mathbf{u}\| = \|\mathbf{v}\| = 1$, $\sigma > 0$ (any rank-one matrix). (a) Show $A^+ = \sigma^{-1}\mathbf{v}\mathbf{u}^\top = A^\top/\|A\|_F^2$. (b) For $A = \begin{pmatrix}1 & 2\\2 & 4\end{pmatrix}$ compute $A^+$, verify $AA^+A = A$, and identify $AA^+$ as a projection onto which line.

**A6.** If $A$ has singular values $\sigma_1, \dots, \sigma_r$, what are the nonzero singular values of (a) $A^\top$, (b) $cA$ for real $c \ne 0$, (c) $QAW$ for orthogonal $Q, W$, (d) $A^{-1}$ (when $A$ is invertible), (e) $A^\top A$, (f) $\begin{pmatrix}A\\0\end{pmatrix}$ (zero rows appended)? One line of justification each, from the SVD.

## B. Core (do all)

**B1.** Prove, using only the SVD and the definition of the spectral norm:
(a) $\sigma_1(A) \ge \max_{i,j}|a_{ij}|$. (Hint: $a_{ij} = \mathbf{e}_i^\top A\mathbf{e}_j$.)
(b) $\|AB\|_2 \le \|A\|_2\|B\|_2$ and $\|AB\|_F \le \|A\|_2\|B\|_F$.
(c) $\kappa_2(AB) \le \kappa_2(A)\,\kappa_2(B)$ for invertible $A, B$.
(d) If $A$ is symmetric with eigenvalues $\lambda_i$, then $\|A\|_2 = \max_i|\lambda_i|$. Give a nonsymmetric $2 \times 2$ matrix for which $\|A\|_2 > \max_i|\lambda_i|$.

**B2.** Compute a full SVD of $A = \begin{pmatrix}1 & 0 & 1\\0 & 1 & 1\end{pmatrix}$ by hand, going through the smaller Gram matrix $AA^\top$ ($2 \times 2$): find $\sigma_1, \sigma_2$ and $\mathbf{u}_1, \mathbf{u}_2$ from $AA^\top$, then $\mathbf{v}_i = A^\top\mathbf{u}_i/\sigma_i$, then complete $\mathbf{v}_3$ from $N(A)$. Verify $\|A\|_F^2 = \sigma_1^2 + \sigma_2^2$ and $A\mathbf{v}_3 = \mathbf{0}$. Write the four fundamental subspaces. Compute $A_1$, $\|A - A_1\|_F$, and $\|A - A_1\|_2$.

**B3.** (Code.) In `project.py`, Part 1, `svd_via_eig(A)` handles $m \ge n$ via $A^\top A$ and refuses $m < n$. (a) Extend it to the wide case by working with $AA^\top$ (as B2 did by hand), and assert agreement with `np.linalg.svd` on random $3 \times 7$ and $2 \times 50$ matrices: compare singular values, the reconstruction $U\Sigma V^\top$, and the projectors $U_rU_r^\top$ and $V_rV_r^\top$. (b) For $A = U\operatorname{diag}(1, 10^{-2}, 10^{-4}, \dots, 10^{-12})V^\top$ with random orthogonal $U, V$, plot the relative error of each computed $\sigma_i$ from `svd_via_eig` and from `np.linalg.svd` against $\sigma_i$ on log–log axes. Explain the slope you see using "the computed eigenvalue $\lambda_i$ of $A^\top A$ carries an absolute error of order $\varepsilon\,\sigma_1^2$ with $\varepsilon \approx 10^{-16}$."

**B4.** (Uniqueness in Eckart–Young.) Let $A = \operatorname{diag}(2, 1)$.
(a) Show that $B_t = \operatorname{diag}(2 - t, 0)$ has rank 1 and $\|A - B_t\|_2 = 1 = \sigma_2$ for every $t \in [0, 1]$. So the best rank-1 approximation in the spectral norm is not unique.
(b) Compute $\|A - B_t\|_F$ and show it is minimized only at $t = 0$.
(c) Prove in general: if $\sigma_k > \sigma_{k+1}$, then $A_k$ is the *unique* rank-$\le k$ minimizer of $\|A - B\|_F$. (Follow the equality cases through the proof of Theorem 3.7: equality in step 3 forces $B = PA$; equality in step 8 forces $c_i = 1$ for $i \le k$ and $c_i = 0$ for $i > k$; conclude $C(P) = \operatorname{span}\{\mathbf{u}_1, \dots, \mathbf{u}_k\}$ and hence $B = A_k$.)
(d) Now let $\sigma_k = \sigma_{k+1}$: take $A = I_2$ and $k = 1$. Show that every $B = \mathbf{u}\mathbf{u}^\top$ with $\|\mathbf{u}\| = 1$ is a best rank-1 approximation in the Frobenius norm, so there are infinitely many. Which step of (c) breaks?

**B5.** (The pseudoinverse.) From $A^+ = V\Sigma^+U^\top$:
(a) Prove the four Penrose equations $AA^+A = A$, $A^+AA^+ = A^+$, $(AA^+)^\top = AA^+$, $(A^+A)^\top = A^+A$.
(b) Prove $(A^+)^+ = A$ and $(A^\top)^+ = (A^+)^\top$.
(c) Show by a counterexample with $A = \begin{pmatrix}1 & 1\end{pmatrix}$ and $B = \begin{pmatrix}1\\0\end{pmatrix}$ that $(AB)^+ \ne B^+A^+$ in general.
(d) Show that every minimizer of $\|A\mathbf{x} - \mathbf{b}\|$ satisfies the normal equations $A^\top A\mathbf{x} = A^\top\mathbf{b}$ (LA-5), and that $A^+\mathbf{b}$ is the unique solution of the normal equations that lies in $C(A^\top)$.

**B6.** (PCA two ways, and what centering does.) Let $X$ be the centered $n \times p$ data matrix with rows $\mathbf{x}_i^\top$ and SVD $U\Sigma V^\top$.
(a) (Pearson, 1901.) Show that for any rank-$k$ orthogonal projector $P$ on $\mathbb{R}^p$, $\sum_i\|\mathbf{x}_i - P\mathbf{x}_i\|^2 = \|X - XP\|_F^2$, that $XP$ has rank $\le k$, and hence by Eckart–Young that the minimum over $P$ is $\sum_{i > k}\sigma_i^2$, attained at $P = V_kV_k^\top$ (check that $XV_kV_k^\top = X_k$).
(b) (Hotelling, 1933.) Show that $\sum_i\|P\mathbf{x}_i\|^2 = \|XP\|_F^2$ is the total variance (times $n - 1$) captured by the subspace $C(P)$, and that (a) and (b) are the same optimization problem because $\|\mathbf{x}_i\|^2 = \|P\mathbf{x}_i\|^2 + \|\mathbf{x}_i - P\mathbf{x}_i\|^2$.
(c) Take the three points $(0, 10), (1, 10), (2, 10)$. Find the first principal component with centering and without centering (i.e. the top eigenvector of $\tilde X^\top\tilde X$). Explain in one sentence what the uncentered "component" is measuring.

**B7.** (Condition numbers.)
(a) Prove $\kappa_2(A) \ge 1$, with equality iff $A = cQ$ for some scalar $c \ne 0$ and orthogonal $Q$.
(b) Prove $\kappa_2(A^\top A) = \kappa_2(A)^2$ for $A$ with independent columns (for non-square $A$, define $\kappa_2 = \sigma_1/\sigma_n$ with $n$ the number of columns). What does this say about solving least squares through the normal equations $A^\top A\mathbf{x} = A^\top\mathbf{b}$ rather than through QR or the SVD, in terms of digits lost?
(c) For $F = \begin{pmatrix}1 & 1\\0 & 1\end{pmatrix}$ compute $\kappa_2$, $\kappa_1 = \|F\|_1\|F^{-1}\|_1$, and $\kappa_\infty = \|F\|_\infty\|F^{-1}\|_\infty$, where $\|A\|_1$ is the largest column sum of $|a_{ij}|$ and $\|A\|_\infty$ the largest row sum. (You may use $\kappa_2 = \sigma_1/\sigma_2$ with the singular values from Example 4.5.)
(d) Show $\kappa_2(cA) = \kappa_2(A)$ for $c \ne 0$. Give a $2 \times 2$ matrix with $\det = 1$ and $\kappa_2 \ge 10^4$, and a matrix with $\det = 10^{-20}$ and $\kappa_2 = 1$.

**B8.** (Code.) Extend Part 3 of `project.py`: (a) add the spectral-norm error $\|A - A_k\|_2$ to the error-versus-rank plot and assert it equals $\sigma_{k+1}$ for every $k$. (b) Part 3 adds Gaussian noise of standard deviation $\tau = 0.05$ (the image has values in $[0, 1]$) to get $\tilde A$ (`noisy`) and plots the singular values of $A$ and $\tilde A$ on one log-scale scree plot. Compute $\|\tilde A_k - A\|_F$ (distance from the truncated *noisy* image to the *clean* one) as a function of $k$, plot it, and find the $k$ that minimizes it; Part 3's last asserts check your answer. Explain in a paragraph why this optimal $k$ is finite: what does truncation remove, and what does it keep, when the singular values of the noise form a flat floor? Then predict the optimal $k$ from the *clean* singular values and the noise level $\tau\sqrt N$ alone.

## C. Challenge (do at least one)

**C1.** (Orthogonal Procrustes.) (a) Prove that the orthogonal matrix nearest to a square $A$ in the Frobenius norm is the polar factor $Q = UV^\top$: expand $\|A - Q\|_F^2 = \|A\|_F^2 - 2\operatorname{tr}(Q^\top A) + n$, write $\operatorname{tr}(Q^\top A) = \operatorname{tr}(Z\Sigma)$ with $Z = V^\top Q^\top U$ orthogonal, and show $\operatorname{tr}(Z\Sigma) \le \sum_i\sigma_i$ with equality iff $Z$ fixes every $\mathbf{e}_i$ with $\sigma_i > 0$. (b) Deduce the solution of $\min_{Q \text{ orthogonal}}\|A - BQ\|_F$ for given $A, B$ of the same shape (hint: expand and reduce to maximizing $\operatorname{tr}(Q^\top B^\top A)$; the answer is the polar factor of $B^\top A$). (c) Implement (b) in the "Your turn" cell of `project.py` and verify it on $A = BQ_0 + \text{noise}$ for a random rotation $Q_0$, and against a brute-force search over rotation angles in two dimensions.

**C2.** (Courant–Fischer, Weyl, and a second proof of Eckart–Young.) Let $M$ be symmetric $n \times n$ with eigenvalues $\lambda_1 \ge \cdots \ge \lambda_n$ and orthonormal eigenvectors $\mathbf{q}_i$.
(a) Prove $\lambda_i = \max_{\dim S = i}\ \min_{\mathbf{w} \in S, \|\mathbf{w}\| = 1}\mathbf{w}^\top M\mathbf{w} = \min_{\dim T = n - i + 1}\ \max_{\mathbf{w} \in T, \|\mathbf{w}\| = 1}\mathbf{w}^\top M\mathbf{w}$. (Hint: for the first, any $i$-dimensional $S$ meets $\operatorname{span}\{\mathbf{q}_i, \dots, \mathbf{q}_n\}$ nontrivially by dimension counting; $S = \operatorname{span}\{\mathbf{q}_1, \dots, \mathbf{q}_i\}$ attains it.)
(b) Deduce for any $m \times n$ matrix $A$: $\sigma_i(A) = \max_{\dim S = i}\min_{\mathbf{w} \in S, \|\mathbf{w}\|=1}\|A\mathbf{w}\| = \min_{\dim T = n - i + 1}\max_{\mathbf{w} \in T, \|\mathbf{w}\|=1}\|A\mathbf{w}\|$.
(c) (Weyl.) Prove $\sigma_{i + j - 1}(A + B) \le \sigma_i(A) + \sigma_j(B)$ whenever $i + j - 1 \le n$. (Hint: intersect $T_A = \operatorname{span}\{\mathbf{v}_i(A), \dots, \mathbf{v}_n(A)\}$ with $T_B = \operatorname{span}\{\mathbf{v}_j(B), \dots, \mathbf{v}_n(B)\}$; on the intersection $\|(A + B)\mathbf{w}\| \le \sigma_i(A) + \sigma_j(B)$; use the min–max form of (b).)
(d) Deduce that if $\operatorname{rank}B \le k$ then $\sigma_i(A - B) \ge \sigma_{i + k}(A)$ for all $i$, and from this both Theorem 3.7 and Theorem 3.8 in two lines each.

**C3.** (Distance to singularity, and the Hilbert matrix.) (a) For square invertible $A$ prove $\sigma_n(A) = \min\{\|E\|_2 : A + E \text{ singular}\}$, so $1/\kappa_2(A) = \min\|E\|_2/\|A\|_2$ is the relative distance to the nearest singular matrix. (Use Theorem 3.8 with $k = n - 1$ for "$\ge$" and an explicit $E$ for "$\le$".) (b) (Code.) The Hilbert matrix $H_n$ has entries $1/(i + j - 1)$, $i, j = 1, \dots, n$ (`scipy.linalg.hilbert`). For $n = 2, \dots, 14$ compute $\kappa_2(H_n)$, solve $H_n\mathbf{x} = H_n\mathbf{1}$ with `np.linalg.solve`, and tabulate the number of correct digits $-\log_{10}\|\mathbf{x} - \mathbf{1}\|_\infty$ next to $\log_{10}\kappa_2$. How well does "digits lost $\approx \log_{10}\kappa$" hold? What happens once $\kappa_2 > 10^{16}$? Fit $\log_{10}\kappa_2(H_n)$ against $n$ and report the growth rate per unit $n$.
