# LA-13 Exercises

Conventions as in the notes: $\varepsilon_{\text{mach}} = 2^{-53} \approx 1.11\times10^{-16}$ is the unit roundoff, $\|\cdot\|$ is the 2-norm unless subscripted, $\kappa(A) = \|A\|\,\|A^{-1}\|$, Householder vectors use the sign rule (3.13), and Example 4.2's matrix is $A = \begin{pmatrix}0.50005 & 0.49995\\0.49995 & 0.50005\end{pmatrix}$ with $b = (1,1)^\top$ and $x = (1,1)^\top$. Every number in the solutions is checked in Part 7 of `project.py`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Spacing of doubles.)
(a) What is the distance between consecutive doubles in $[2^{10}, 2^{11})$? In $[2^{-10}, 2^{-9})$?
(b) What is the smallest positive integer that is not a double? What does `float(2**53 + 1)` return, and why that value rather than the other neighbour?
(c) `0.1 + 0.2 == 0.3` is `False` in Python. Use `decimal.Decimal(0.1)` and friends to see the stored values, and explain why the computed sum lands exactly one spacing above the stored $0.3$.

**A2.** (Whose fault is the cancellation?) Let $f(x) = 1 - \cos x$ and $x = 10^{-8}$.
(a) Evaluate `1 - np.cos(1e-8)` and explain the result.
(b) Show $f(x) = 2\sin^2(x/2)$ and evaluate that instead.
(c) Compute the condition number $\kappa(f; x) = |x f'(x)/f(x)|$. Is the *problem* ill-conditioned? Which single step of the naive algorithm is, and what is its condition number by (3.4)?

**A3.** (Condition number of a nearly singular matrix.) Let $A = \begin{pmatrix}1 & 1\\1 & 1.0001\end{pmatrix}$.
(a) Compute $A^{-1}$ exactly.
(b) Compute $\kappa_\infty(A) = \|A\|_\infty\|A^{-1}\|_\infty$ (max absolute row sum).
(c) Compute $\kappa_2(A)$ from the eigenvalues ($A$ is symmetric positive definite).
(d) About how many correct digits does a backward stable solver deliver for $Ax = b$?

**A4.** (One reflection.) For $x = (-3, 4, 0)^\top$, compute the Householder vector $v$ by (3.13), $v^\top v$, the matrix $H$, and $Hx$. Across which plane does $H$ reflect? What would the other sign of $v$ have produced?

**A5.** (One sweep each.) For Example 4.5's system $\begin{pmatrix}4&1\\1&3\end{pmatrix}x = \begin{pmatrix}1\\2\end{pmatrix}$:
(a) Do one Jacobi step and one Gauss–Seidel step from $x_0 = 0$.
(b) Write both iteration matrices $B = I - M^{-1}A$ and find their spectral radii.
(c) How many decimal digits per step does each gain asymptotically? How are the two radii related?

**A6.** (Shifts.) A symmetric matrix has eigenvalues $10$, $9$ and $1$.
(a) How many steps of power iteration reduce the eigenvector error by a factor $10^{-6}$?
(b) How many steps of inverse iteration with shift $\mu = 9.8$?
(c) With $\mu = 9.4$ instead, which eigenvector does inverse iteration find, and how many steps does it need?

## B. Core (do all)

**B1.** (The residual is the backward error.) Let $A$ be invertible, $\hat x \ne 0$ any approximate solution of $Ax = b$, and $r = b - A\hat x$.
(a) Show that $\delta A = r\hat x^\top/\|\hat x\|^2$ satisfies $(A + \delta A)\hat x = b$ and $\|\delta A\| = \|r\|/\|\hat x\|$.
(b) Show that no smaller perturbation works: if $(A + E)\hat x = b$ then $\|E\| \ge \|r\|/\|\hat x\|$. Conclude that the backward error of $\hat x$ is exactly $\|r\|/(\|A\|\|\hat x\|)$.
(c) Take Example 4.2 and $\hat x = (2, 0)^\top$, the solution of the worst-case perturbed problem. Compute $r$, the backward error, and the forward error $\|\hat x - x\|/\|\hat x\|$. Show that forward error $= \kappa\times$ backward error *exactly*, and say which theorem's bound is attained.

**B2.** (Both perturbations at once; distance to singularity.)
(a) Let $Ax = b$ and $(A + \delta A)(x + \delta x) = b + \delta b$ with $\|A^{-1}\|\|\delta A\| < 1$. Prove
$$\frac{\|\delta x\|}{\|x\|} \le \frac{\kappa(A)}{1 - \kappa(A)\|\delta A\|/\|A\|}\left(\frac{\|\delta A\|}{\|A\|} + \frac{\|\delta b\|}{\|b\|}\right).$$
(b) State the first-order version and compare it with Theorems 3.5 and 3.6.
(c) Prove that $\min\{\|\delta A\| : A + \delta A \text{ singular}\} = \sigma_n(A)$, so the relative distance from $A$ to the nearest singular matrix is $1/\kappa_2(A)$. Find that nearest singular matrix for Example 4.2's $A$.

**B3.** (Partial pivoting's worst case is real.) Let $W_n$ have $1$ on the diagonal, $-1$ everywhere below it, $1$ in the last column, and $0$ elsewhere.
(a) For $n = 4$, run Gaussian elimination with partial pivoting by hand. Show that no row swaps are needed and find $L$, $U$ and the growth factor $\rho$.
(b) Prove that for general $n$ the last column of $U$ is $(1, 2, 4, \dots, 2^{n-1})^\top$, so $\rho = 2^{n-1}$.
(c) In NumPy, take $n = 60$, a random $x$, $b = W_{60}x$, and solve with `np.linalg.solve` (LAPACK, partial pivoting) and with `householder_qr`. Report $\kappa(W_{60})$ and the forward and backward errors. Reconcile what you see with Theorem 3.9 and (3.9).

**B4.** (Läuchli's matrix: when the normal equations break.) Let $\delta = 10^{-8}$,
$$A = \begin{pmatrix}1 & 1\\ \delta & 0\\ 0 & \delta\end{pmatrix}, \qquad b = \begin{pmatrix}2\\ \delta\\ \delta\end{pmatrix},$$
so that $x = (1, 1)^\top$ solves $Ax = b$ exactly.
(a) Compute $A^\top A$ exactly, then in double precision. What happens when you try to solve the normal equations?
(b) Find the singular values of $A$ and $\kappa_2(A)$. What accuracy does (3.9) predict for a backward stable least-squares solver?
(c) Solve by Householder QR and explain the difference with Theorem 3.12.

**B5.** (Ritz values move outward.) Run Lanczos (3.33) on symmetric $A$ and let $\theta^{(k)}_1 \ge \cdots \ge \theta^{(k)}_k$ be the eigenvalues of $T_k$.
(a) Show that $T_k$ is the leading $k\times k$ principal submatrix of $T_{k+1}$.
(b) Prove $\theta^{(k)}_1 \le \theta^{(k+1)}_1 \le \lambda_{\max}(A)$ and $\theta^{(k)}_k \ge \theta^{(k+1)}_{k+1} \ge \lambda_{\min}(A)$.
(c) For SPD $A$, conclude that $\theta^{(k)}_1/\theta^{(k)}_k$ is a lower bound for $\kappa(A)$ that never decreases with $k$. Check it against the project's Lanczos run on the $N = 100$ Laplacian.

**B6.** (A perfect shift.)
(a) Let $A = \begin{pmatrix}3 & 1\\1 & 1\end{pmatrix}$ (Example 4.4) and $\mu = 2 - \sqrt2$. Compute the QR factorization of $A - \mu I$ by hand (positive $r_{11}$), show $r_{22} = 0$, and show that $A_1 = RQ + \mu I = \operatorname{diag}(2 + \sqrt2, 2 - \sqrt2)$.
(b) Prove in general: if $H$ is upper Hessenberg with nonzero subdiagonal and $\mu$ is an eigenvalue, then one QR step with shift $\mu$ produces a matrix whose last row is $(0, \dots, 0, \mu)$.
(c) Show that $\begin{pmatrix}0&1\\1&0\end{pmatrix}$ is a fixed point of the unshifted QR algorithm, and also of the algorithm with the Rayleigh shift $\mu = a_{nn}$. Compute its Wilkinson shift (3.19) and say what one step with it does.

**B7.** (Steepest descent zigzags.) Use Example 4.5's system and $x_0 = 0$.
(a) Run two steps of steepest descent (§3.10(ii)) in exact fractions.
(b) Show that $r_2$ is a multiple of $r_0$. Deduce that the iteration repeats itself, scaled, forever, and never reaches $x_*$.
(c) Compute $\|e_k\|_A^2$ for $k = 0, 1, 2$, the per-step reduction factor, and compare it with the bound $(\kappa - 1)/(\kappa + 1)$.

**B8.** (Code: Hessenberg reduction for nonsymmetric matrices.)
(a) Copy `householder_tridiagonalize` from `project.py` and turn it into `householder_hessenberg(A)` for nonsymmetric $A$. Which two lines must change, and why? Check $QHQ^\top = A$, $Q^\top Q = I$, and that $H$ is upper Hessenberg.
(b) With `rng = np.random.default_rng(7)`, build $A = X\operatorname{diag}(6, 5, 4, 3, 2, 1)X^{-1}$ for a random $6\times6$ $X$, reduce it, and run 150 unshifted QR steps on $H$. Verify that every iterate stays Hessenberg (Theorem 3.17) and that the diagonal converges to $6, 5, \dots, 1$.
(c) Explain why one QR step on an $n\times n$ Hessenberg matrix costs $O(n^2)$ flops, against $O(n^3)$ for a full matrix.

## C. Challenge (do at least one)

**C1.** (Jacobi can diverge on an SPD matrix; Gauss–Seidel cannot.) Let $A = (1 - a)I + a\mathbf 1\mathbf 1^\top$ be $3\times3$, where $\mathbf 1 = (1,1,1)^\top$.
(a) Find the eigenvalues of $A$ and the values of $a$ for which $A$ is SPD.
(b) Find the eigenvalues of the Jacobi iteration matrix and the values of $a$ for which Jacobi converges.
(c) Prove the Ostrowski–Reich theorem for this splitting: if $A$ is SPD then Gauss–Seidel converges. (Write $A = L + D + L^\top$ with $L$ strictly lower triangular, $M = D + L$, $B = I - M^{-1}A$. For an eigenpair $Bx = \lambda x$, set $y = M^{-1}Ax$ and show $(1 - |\lambda|^2)\,x^*Ax = |1 - \lambda|^2\,x^*Dx$.)
(d) For $a = 0.8$, compute $\rho(B_J)$ and $\rho(B_{GS})$ numerically.

**C2.** (The Chebyshev facts behind (3.30).)
(a) Show that $T_k(\cos\theta) = \cos k\theta$ defines a polynomial of degree $k$, using $T_{k+1}(x) = 2xT_k(x) - T_{k-1}(x)$.
(b) Show that for $x \ge 1$, $T_k(x) = \cosh(k\operatorname{arccosh}x) = \tfrac12\big[(x + \sqrt{x^2-1})^k + (x - \sqrt{x^2-1})^k\big]$.
(c) Prove the minimax property used in Theorem 3.24: among polynomials of degree at most $k$ with $p(0) = 1$, the maximum of $|p|$ on $[\lambda_{\min}, \lambda_{\max}]$ (an interval not containing $0$) is minimized by $T_k(\ell(\lambda))/T_k(\gamma)$.
(d) For Example 4.6 ($\kappa = 178.1$), show $\operatorname{arccosh}\gamma = \ln\frac{\sqrt\kappa + 1}{\sqrt\kappa - 1}$, find the smallest $k$ with $1/T_k(\gamma) \le 10^{-8}$, and compare with the $k$ from the bound $2\rho^k \le 10^{-8}$. Why do they agree?

**C3.** (Code: why the power iterations must be re-orthonormalized.) Modify `randomized_svd` so that step 3 of Algorithm 3.4 forms $Y = (AA^\top)^qA\Omega$ directly, with a single QR at the end. On the project's $400\times300$ test matrix with $\sigma_j = 1/j$, $k = 10$, $p = 10$, compare the mean of $\|A - QQ^\top A\|_2/\sigma_{11}$ over 5 trials with and without re-orthonormalization for $q = 0, 4, 8, 12$. Explain the pattern quantitatively using $\varepsilon_{\text{mach}}$, and predict how many directions survive for $q = 8$ and $q = 12$.
