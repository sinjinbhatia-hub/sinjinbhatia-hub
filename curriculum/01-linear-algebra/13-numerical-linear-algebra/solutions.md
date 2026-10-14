# LA-13 Solutions

Every number below is asserted or printed in Part 7 of `project.py` (and Examples 4.2, 4.4, 4.5 in Parts 2–5).

## A. Warm-up

**A1.** (a) A double in $[2^e, 2^{e+1})$ is $(1.b_1\cdots b_{52})_2\times2^e$, so consecutive doubles differ by one unit in the 52nd bit: $2^{e-52}$. For $e = 10$ the spacing is $2^{-42} = 2.274\times10^{-13}$; for $e = -10$ it is $2^{-62} = 2.168\times10^{-19}$. Same relative spacing, wildly different absolute spacing: this is why rounding error is relative (3.1).

(b) Every integer up to $2^{53}$ fits in a 53-bit significand, so all of them are doubles. In $[2^{53}, 2^{54})$ the spacing is $2^{53-52} = 2$, so $2^{53} + 1$ is the first integer that is not a double. It lies exactly halfway between $2^{53}$ and $2^{53} + 2$. Ties go to the neighbour whose last significand bit is $0$: $2^{53} = (1.00\cdots0)_2\times2^{53}$ is even and $2^{53} + 2 = (1.00\cdots01)_2\times2^{53}$ is odd. So `float(2**53 + 1)` returns $2^{53} = 9007199254740992$.

(c) None of $0.1$, $0.2$, $0.3$ is a double (their binary expansions repeat). The stored values are
$$\mathrm{fl}(0.1) = 0.1000000000000000055511\ldots,\quad \mathrm{fl}(0.2) = 0.2000000000000000111022\ldots,\quad \mathrm{fl}(0.3) = 0.2999999999999999888977\ldots$$
The exact sum of the first two lies $2.7756\times10^{-17} = 2^{-55}$ above $\mathrm{fl}(0.3)$. In $[0.25, 0.5)$ the spacing is $2^{-54} = 5.551\times10^{-17}$, so the exact sum is exactly halfway between $\mathrm{fl}(0.3)$ and the next double, $0.3000000000000000444089\ldots$. The tie is broken toward the even significand, which is the upper one. Hence $\mathrm{fl}(0.1 + 0.2) - \mathrm{fl}(0.3) = 2^{-54} = 5.551\times10^{-17}$: one spacing. The lesson is not "floating point is broken". Each operation obeyed (3.2); the inputs were already rounded. Never test floats for equality; test $|a - b| \le \text{tol}\cdot|b|$.

**A2.** (a) $\cos(10^{-8}) = 1 - 5\times10^{-17} + O(10^{-33})$. The doubles just below $1$ are spaced $2^{-53} = 1.11\times10^{-16}$ apart, and $5\times10^{-17}$ is less than half that spacing, so $\mathrm{fl}(\cos 10^{-8}) = 1$ and the computed $1 - \cos x$ is $0.0$. Relative error $100\%$.

(b) From $\cos x = 1 - 2\sin^2(x/2)$, $f(x) = 2\sin^2(x/2)$. Evaluating it: $\sin(5\times10^{-9}) = 5\times10^{-9}$ to double precision, squared and doubled gives $5.0000\times10^{-17}$, correct to full precision. No subtraction, no cancellation.

(c) $f'(x) = \sin x$, so $\kappa(f; x) = \dfrac{x\sin x}{1 - \cos x} \approx \dfrac{x\cdot x}{x^2/2} = 2$; the project gets $2.0000$. The *problem* is perfectly conditioned: a relative change $\delta$ in $x$ changes $f$ by $2\delta$. The ill-conditioned step is the subtraction $1 - \mathrm{fl}(\cos x)$, with
$$\kappa_{\text{sub}} = \frac{|1| + |\cos x|}{|1 - \cos x|} \approx \frac{2}{5\times10^{-17}} = 4\times10^{16}$$
by (3.4). The rounding of $\cos x$ (relative error up to $1.1\times10^{-16}$) is amplified by $4\times10^{16}$ into an order-one error. This is an *unstable algorithm for a well-conditioned problem*: the opposite of Example 4.2, and the algorithm's fault. The cure, as in (3.5), is to rewrite the formula so the dangerous subtraction never happens.

**A3.** (a) $\det A = 1.0001 - 1 = 10^{-4}$, so
$$A^{-1} = \frac{1}{10^{-4}}\begin{pmatrix}1.0001 & -1\\ -1 & 1\end{pmatrix} = \begin{pmatrix}10001 & -10000\\ -10000 & 10000\end{pmatrix}.$$
(b) Row sums: $\|A\|_\infty = \max(2, 2.0001) = 2.0001$ and $\|A^{-1}\|_\infty = \max(20001, 20000) = 20001$. So $\kappa_\infty = 2.0001\times20001 = 40004.0001$.

(c) Trace $2.0001$, determinant $10^{-4}$, so $\lambda = \tfrac12\big(2.0001 \pm\sqrt{2.0001^2 - 4\times10^{-4}}\big) = \tfrac12\big(2.0001 \pm \sqrt{4.00000001}\big)$. The large one is $\lambda_1 = 2.000050$; the small one is best computed as $\det/\lambda_1 = 10^{-4}/2.00005 = 4.999875\times10^{-5}$ (the stable-quadratic trick of (3.5) again). For SPD $A$, $\sigma_i = \lambda_i$ and $\kappa_2 = \lambda_1/\lambda_2 = \lambda_1^2/\det A = 40002.0$.

(d) $16 - \log_{10}(40002) \approx 16 - 4.6 = 11.4$ correct digits. The two norms give nearly the same $\kappa$; they always agree within a factor depending only on $n$, so the choice of norm never changes the number of digits by more than about one.

**A4.** $\|x\| = 5$ and $x_1 = -3 < 0$, so $\operatorname{sign}(x_1) = -1$ and (3.13) gives
$$v = x - 5\mathbf e_1 = (-8, 4, 0)^\top, \qquad v^\top v = 64 + 16 = 80, \qquad v^\top x = 24 + 16 = 40 = \tfrac12v^\top v.$$
So $Hx = x - 2v\frac{v^\top x}{v^\top v} = x - v = (5, 0, 0)^\top = -\operatorname{sign}(x_1)\|x\|\mathbf e_1$. The matrix:
$$H = I - \frac{2}{80}vv^\top = I - \frac{1}{40}\begin{pmatrix}64 & -32 & 0\\ -32 & 16 & 0\\ 0 & 0 & 0\end{pmatrix} = \begin{pmatrix}-0.6 & 0.8 & 0\\ 0.8 & 0.6 & 0\\ 0 & 0 & 1\end{pmatrix}.$$
$H$ reflects across the plane $v^\perp = \{y : -8y_1 + 4y_2 = 0\} = \{y_2 = 2y_1\}$, spanned by $(1, 2, 0)^\top$ and $\mathbf e_3$; it acts as a reflection of the $(y_1, y_2)$-plane and leaves $y_3$ alone. The other sign, $v' = x + 5\mathbf e_1 = (2, 4, 0)^\top$, gives $H'x = (-5, 0, 0)^\top$. Its first entry $v'_1 = x_1 + \|x\| = -3 + 5$ subtracts two magnitudes, with $\kappa_{\text{sub}} = 8/2 = 4$: harmless here, but catastrophic when $x \approx -\|x\|\mathbf e_1$, which is exactly the case in which (3.13) picks the other sign.

**A5.** (a) Jacobi ($M = D = \operatorname{diag}(4, 3)$): $x_1 = D^{-1}b = (1/4, 2/3)$. Gauss–Seidel: first component $x^{(1)}_1 = 1/4$ as before; the second uses it immediately, $x^{(1)}_2 = (2 - 1\cdot\tfrac14)/3 = 7/12$. So $x_1^{GS} = (1/4, 7/12)$.

(b) $B_J = I - D^{-1}A = \begin{pmatrix}0 & -1/4\\ -1/3 & 0\end{pmatrix}$, with $\lambda^2 = 1/12$, so $\rho_J = 1/\sqrt{12} = 0.2887$. For Gauss–Seidel, $M = D + L = \begin{pmatrix}4 & 0\\ 1 & 3\end{pmatrix}$, $M^{-1} = \begin{pmatrix}1/4 & 0\\ -1/12 & 1/3\end{pmatrix}$, and $B_{GS} = M^{-1}N$ with $N = M - A = \begin{pmatrix}0 & -1\\ 0 & 0\end{pmatrix}$:
$$B_{GS} = \begin{pmatrix}0 & -1/4\\ 0 & 1/12\end{pmatrix}, \qquad \rho_{GS} = 1/12 = 0.0833.$$
(c) Digits per step $= -\log_{10}\rho$: $0.54$ for Jacobi, $1.08$ for Gauss–Seidel. And $\rho_{GS} = \rho_J^2$, the same relation as for the model problem in §3.9: Gauss–Seidel is exactly twice as fast.

**A6.** (a) Power iteration converges at rate $|\lambda_2/\lambda_1| = 0.9$. We need $0.9^k \le 10^{-6}$, i.e. $k \ge \ln(10^{-6})/\ln 0.9 = 131.1$: $132$ steps.

(b) With $\mu = 9.8$ the eigenvalues of $(A - \mu I)^{-1}$ are $1/(10 - 9.8) = 5$, $1/(9 - 9.8) = -1.25$, $1/(1 - 9.8) = -0.114$. The rate is $1.25/5 = 0.25$, and $0.25^k \le 10^{-6}$ needs $k \ge 9.97$: $10$ steps, each a solve with a factored $A - \mu I$.

(c) With $\mu = 9.4$ the distances are $0.6$ (to $10$), $0.4$ (to $9$), $8.4$ (to $1$). The closest is $9$, so inverse iteration converges to the eigenvector of $\lambda = 9$, at rate $0.4/0.6 = 2/3$: $\ln(10^{-6})/\ln(2/3) = 34.07$, so $35$ steps. A shift finds the eigenvalue *nearest to it*, not the one you hoped for, and a shift near the midpoint between two eigenvalues is slow.

## B. Core

**B1.** (a) $(A + \delta A)\hat x = A\hat x + r\dfrac{\hat x^\top\hat x}{\|\hat x\|^2} = A\hat x + r = b$. For the norm: a rank-one matrix $uw^\top$ satisfies $\|uw^\top z\| = \|u\|\,|w^\top z| \le \|u\|\|w\|\|z\|$ by Cauchy–Schwarz, with equality at $z = w$, so $\|uw^\top\| = \|u\|\|w\|$. Hence $\|\delta A\| = \|r\|\|\hat x\|/\|\hat x\|^2 = \|r\|/\|\hat x\|$.

(b) If $(A + E)\hat x = b$ then $E\hat x = b - A\hat x = r$, so $\|r\| = \|E\hat x\| \le \|E\|\|\hat x\|$ and $\|E\| \ge \|r\|/\|\hat x\|$. Part (a) attains this, so the smallest relative perturbation of $A$ that makes $\hat x$ exact is
$$\frac{\|\delta A\|}{\|A\|} = \frac{\|r\|}{\|A\|\,\|\hat x\|},$$
which is the backward error (Definition 3.2 with the data perturbation placed in $A$; Higham §7.1 gives the version that also perturbs $b$).

(c) $A\hat x = (2\times0.50005,\ 2\times0.49995) = (1.0001, 0.9999)$, so $r = (-10^{-4}, 10^{-4})$ and $\|r\| = \sqrt2\times10^{-4}$. With $\|A\| = 1$ and $\|\hat x\| = 2$, the backward error is $\sqrt2\times10^{-4}/2 = 7.0711\times10^{-5}$. The forward error is $\|\hat x - x\|/\|\hat x\| = \|(1, -1)\|/2 = \sqrt2/2 = 0.7071$, and $\kappa\times$ backward error $= 10^4\times7.0711\times10^{-5} = 0.7071$. Equal.

This is Theorem 3.6's bound (3.8) attained: $\hat x - x = -A^{-1}\delta A\,\hat x = -A^{-1}r$. Both inequalities in its proof are tight because $r \propto (-1, 1)$ is the direction $A^{-1}$ stretches by the full $10^4$, and $\|\delta A\,\hat x\| = \|\delta A\|\|\hat x\|$ since $\delta A$ is rank one along $\hat x$. The lesson: a residual that looks small ($7\times10^{-5}$ relative) certifies only a small backward error. The forward error can be $\kappa$ times larger.

**B2.** (a) Expand $(A + \delta A)(x + \delta x) = b + \delta b$ and subtract $Ax = b$:
1. $A\,\delta x = \delta b - \delta A\,x - \delta A\,\delta x$, so $\delta x = A^{-1}(\delta b - \delta A\,x - \delta A\,\delta x)$.
2. Take norms: $\|\delta x\| \le \|A^{-1}\|\big(\|\delta b\| + \|\delta A\|\|x\| + \|\delta A\|\|\delta x\|\big)$.
3. Move the $\delta x$ term left: $\|\delta x\|\big(1 - \|A^{-1}\|\|\delta A\|\big) \le \|A^{-1}\|\big(\|\delta b\| + \|\delta A\|\|x\|\big)$. The bracket on the left is positive by hypothesis.
4. Divide by $\|x\|$ and use $\|b\| \le \|A\|\|x\|$, i.e. $\|\delta b\|/\|x\| \le \|A\|\,\|\delta b\|/\|b\|$:
$$\frac{\|\delta x\|}{\|x\|} \le \frac{\|A^{-1}\|\|A\|}{1 - \|A^{-1}\|\|\delta A\|}\left(\frac{\|\delta b\|}{\|b\|} + \frac{\|\delta A\|}{\|A\|}\right).$$
5. Write $\|A^{-1}\|\|\delta A\| = \kappa(A)\|\delta A\|/\|A\|$ in the denominator. ∎

(b) To first order the denominator is $1$: $\|\delta x\|/\|x\| \lesssim \kappa(A)\big(\|\delta A\|/\|A\| + \|\delta b\|/\|b\|\big)$. With $\delta A = 0$ this is exactly Theorem 3.5. With $\delta b = 0$ it is Theorem 3.6, measured relative to $x$ instead of $\tilde x$ (the two agree to first order). Relative errors in $A$ and $b$ simply add, each amplified by $\kappa$. The denominator says the bound breaks down when $\kappa\|\delta A\|/\|A\|$ approaches $1$, which part (c) explains: that is when $A + \delta A$ can be singular.

(c) *Lower bound.* If $A + \delta A$ is singular, pick a unit $z$ with $(A + \delta A)z = 0$. Then $Az = -\delta A\,z$, so $\sigma_n = \min_{\|y\| = 1}\|Ay\| \le \|Az\| = \|\delta A\,z\| \le \|\delta A\|$ (LA-9 Theorem 3.5 for the first equality). *Attained.* With the SVD $A = \sum_i\sigma_iu_iv_i^\top$, take $\delta A = -\sigma_nu_nv_n^\top$, of norm $\sigma_n$. Then $(A + \delta A)v_n = \sigma_nu_n - \sigma_nu_n = 0$, so $A + \delta A$ is singular. Hence the minimum is $\sigma_n$, and relative to $\|A\| = \sigma_1$ it is $\sigma_n/\sigma_1 = 1/\kappa_2(A)$.

For Example 4.2's $A$ (SPD, so $u_2 = v_2 = (1, -1)/\sqrt2$ and $\sigma_2 = 10^{-4}$):
$$\delta A = -10^{-4}\cdot\frac12\begin{pmatrix}1 & -1\\ -1 & 1\end{pmatrix}, \qquad A + \delta A = \begin{pmatrix}0.5 & 0.5\\ 0.5 & 0.5\end{pmatrix},$$
which is visibly singular; the project confirms $\sigma_{\min}(A + \delta A) = 1.7\times10^{-17}$, zero to rounding. A condition number of $10^4$ means: a $0.01\%$ change in the matrix makes it singular.

**B3.** (a)
$$W_4 = \begin{pmatrix}1 & 0 & 0 & 1\\ -1 & 1 & 0 & 1\\ -1 & -1 & 1 & 1\\ -1 & -1 & -1 & 1\end{pmatrix}.$$
*Column 1.* All four candidates have $|{\pm1}| = 1$, a tie; partial pivoting keeps the first, so no swap. Multipliers $l_{21} = l_{31} = l_{41} = -1$; adding row 1 to rows 2–4 gives rows $(0, 1, 0, 2)$, $(0, -1, 1, 2)$, $(0, -1, -1, 2)$.
*Column 2.* Pivot $1$ ties with the $-1$s below; no swap. $l_{32} = l_{42} = -1$; adding row 2 gives $(0, 0, 1, 4)$ and $(0, 0, -1, 4)$.
*Column 3.* $l_{43} = -1$; adding row 3 gives $(0, 0, 0, 8)$.
$$L = \begin{pmatrix}1 & & & \\ -1 & 1 & & \\ -1 & -1 & 1 & \\ -1 & -1 & -1 & 1\end{pmatrix}, \qquad U = \begin{pmatrix}1 & 0 & 0 & 1\\ 0 & 1 & 0 & 2\\ 0 & 0 & 1 & 4\\ 0 & 0 & 0 & 8\end{pmatrix}, \qquad P = I,$$
and $\rho = \max|u_{ij}|/\max|a_{ij}| = 8/1 = 8 = 2^3$.

(b) Claim: before step $k$, row $k$ is $(0, \dots, 0, 1, 0, \dots, 0, 2^{k-1})$ (the $1$ in position $k$), and each row $i > k$ has zeros in columns $1..k-1$, $-1$ in columns $k..i-1$, and $2^{k-1}$ in column $n$; if $i < n$ it also has $1$ in column $i$ and zeros in columns $i+1..n-1$. (For $i = n$ the diagonal entry *is* the last-column entry.) For $k = 1$ this is $W_n$ itself. Step $k$: column $k$ has $1$ on the diagonal and $-1$ below, a tie, so no swap; every multiplier is $-1$, so row $i \leftarrow$ row $i$ + row $k$. Row $k$ has zeros in columns $k+1..n-1$, so those entries of row $i$ are unchanged; column $k$ becomes $0$; the last column becomes $2^{k-1} + 2^{k-1} = 2^k$. That is the claim for $k+1$. At the end $u_{kn} = 2^{k-1}$, so the last column of $U$ is $(1, 2, \dots, 2^{n-1})^\top$ and $\rho = 2^{n-1}$. ∎

(c) The project ($n = 60$): $\kappa(W_{60}) = 26.8$, a perfectly well-conditioned matrix that should give 14–15 digits. LAPACK's partial pivoting returns forward error $0.26$ and backward error $1.4\times10^{-2}$. Householder QR returns forward error $2.6\times10^{-15}$. Reconciliation: Theorem 3.9 bounds the backward error by $O(\rho\varepsilon)$, and here $\rho\varepsilon = 2^{59}\cdot2^{-53} = 64$: no guarantee at all. Entries of size up to $2^{59}$ in $U$ swamp the order-one data, whose low-order bits are rounded away. The measured backward error $1.4\times10^{-2}$ is about $10^{14}\varepsilon$, and by (3.9) the forward error is at most of order $\kappa\times1.4\times10^{-2} \approx 0.4$, consistent with $0.26$. Householder QR applies only orthogonal transformations, has no growth factor, and is backward stable for every matrix (Theorem 3.11). Moral: "GEPP is backward stable" means "in practice"; this matrix is the counterexample.

**B4.** (a) Exactly,
$$A^\top A = \begin{pmatrix}1 + \delta^2 & 1\\ 1 & 1 + \delta^2\end{pmatrix}, \qquad A^\top b = (2 + \delta^2)\begin{pmatrix}1\\1\end{pmatrix}.$$
With $\delta^2 = 10^{-16}$, below half the spacing $2.2\times10^{-16}$ of doubles just above $1$, $\mathrm{fl}(1 + \delta^2) = 1$. The computed $A^\top A$ is $\begin{pmatrix}1 & 1\\1 & 1\end{pmatrix}$, exactly singular, and `np.linalg.solve` raises `LinAlgError`. All information about $\delta$ was destroyed when the product was formed, before any solver ran.

(b) The eigenvalues of $A^\top A$ are $(1 + \delta^2) \pm 1$, i.e. $2 + \delta^2$ (eigenvector $(1,1)$) and $\delta^2$ (eigenvector $(1,-1)$). So $\sigma_1 = \sqrt{2 + \delta^2}$, $\sigma_2 = \delta$, and $\kappa_2(A) = \sqrt{2 + \delta^2}/\delta = 1.414\times10^8$. Since the residual is zero, a backward stable least-squares solver should deliver forward error about $\kappa\varepsilon \approx 1.6\times10^{-8}$: eight digits.

(c) Householder QR returns $x = (1, 1)$ to full precision (relative error below $10^{-14}$, better than the worst case). Theorem 3.12 is the explanation: $\kappa(A^\top A) = \kappa(A)^2 = 2\times10^{16} > 1/\varepsilon$, so the normal equations ask for more than 16 digits, while QR works with $A$ itself and its $\kappa = 1.4\times10^8$.

**B5.** (a) The Lanczos coefficients $\alpha_1, \dots, \alpha_k$ and $\beta_1, \dots, \beta_{k-1}$ are computed at steps $1..k$ and never revised; step $k+1$ appends $\beta_k$ and $\alpha_{k+1}$. Equivalently $(T_{k+1})_{ij} = q_i^\top Aq_j$ for $i, j \le k+1$, and restricting to $i, j \le k$ gives $T_k$. So $T_k$ is the leading principal $k\times k$ submatrix of $T_{k+1}$.

(b) By LA-8 Theorem 3.17, $\theta^{(k)}_1 = \max_{y\in\mathbb R^k, y\ne0}y^\top T_ky/y^\top y$. For any such $y$, the padded vector $\tilde y = (y, 0) \in \mathbb R^{k+1}$ has $\tilde y^\top T_{k+1}\tilde y = y^\top T_ky$ by (a), and $\tilde y^\top\tilde y = y^\top y$. So the maximum over $\mathbb R^{k+1}$ is at least the maximum over $\mathbb R^k$: $\theta^{(k)}_1 \le \theta^{(k+1)}_1$. With $x = Q_{k+1}\tilde y$, $\theta^{(k+1)}_1 = \max_{x\in\mathcal K_{k+1}}x^\top Ax/x^\top x \le \max_{x\in\mathbb R^n}x^\top Ax/x^\top x = \lambda_{\max}(A)$. The statements about the smallest Ritz value are the same argument with min in place of max. (Cauchy interlacing, LA-8 Exercise C1, gives the full interlacing pattern.)

(c) For SPD $A$ all Ritz values are positive, and $\lambda_{\min} \le \theta^{(k)}_k \le \theta^{(k)}_1 \le \lambda_{\max}$, so $\theta^{(k)}_1/\theta^{(k)}_k \le \kappa(A)$. By (b) the numerator never decreases and the denominator never increases, so the estimate never decreases. The project checks the monotonicity at every $k$ from $1$ to $200$ on the $N = 100$ Laplacian: the estimate is $574$ at $k = 30$ and $4134$ at $k = 200$, equal to the exact $\kappa$. Thirty matrix–vector products already tell you $\kappa \ge 574$. That is a lower bound on $\kappa$, so it says how much the guarantee (3.30) *cannot* promise: with $\kappa \ge 574$, $\rho \ge 0.920$, and $2\rho^k \le 10^{-8}$ needs $k \ge 229$. (It does not prove CG will be that slow; a clustered spectrum can beat (3.30). On this Laplacian, whose spectrum fills the interval, CG in fact needs about $300$ iterations: $297$ in the project's $N = 100$ run of Part 5.)

**B6.** (a) $A - \mu I = \begin{pmatrix}1 + \sqrt2 & 1\\ 1 & \sqrt2 - 1\end{pmatrix}$, with determinant $(\sqrt2 + 1)(\sqrt2 - 1) - 1 = 0$. QR by Gram–Schmidt:
- $r_{11} = \sqrt{(1 + \sqrt2)^2 + 1} = \sqrt{4 + 2\sqrt2} = 2.6131$, $q_1 = (1 + \sqrt2, 1)/r_{11}$.
- $r_{12} = q_1^\top(1, \sqrt2 - 1) = \big((1 + \sqrt2) + (\sqrt2 - 1)\big)/r_{11} = 2\sqrt2/r_{11} = 1.0824$.
- $r_{22} = |\det(A - \mu I)|/r_{11} = 0$: the second column is a multiple of the first. Take $q_2 = (-1, 1 + \sqrt2)/r_{11}$, the unit vector orthogonal to $q_1$.

Now $RQ$. Its second row is $0$ because $R$'s is. Its first row is $(r_{11}q_{11} + r_{12}q_{21},\ r_{11}q_{12} + r_{12}q_{22})$. Using $r_{12}/r_{11} = 2\sqrt2/(4 + 2\sqrt2) = \sqrt2/(2 + \sqrt2) = \sqrt2 - 1$:
$$(RQ)_{11} = (1 + \sqrt2) + (\sqrt2 - 1) = 2\sqrt2, \qquad (RQ)_{12} = -1 + (\sqrt2 - 1)(1 + \sqrt2) = -1 + 1 = 0.$$
So $A_1 = RQ + \mu I = \begin{pmatrix}2\sqrt2 + 2 - \sqrt2 & 0\\ 0 & 2 - \sqrt2\end{pmatrix} = \operatorname{diag}(2 + \sqrt2,\ 2 - \sqrt2)$: both eigenvalues in one step, against eleven unshifted steps to reach $10^{-8}$ (Example 4.4).

(b) Let $B = H - \mu I$, Hessenberg with nonzero subdiagonal entries $h_{j+1,j}$. Rows $2..n$ of columns $1..n-1$ of $B$ form an upper triangular $(n-1)\times(n-1)$ matrix with diagonal $h_{21}, \dots, h_{n,n-1} \ne 0$, so the first $n - 1$ columns of $B$ are linearly independent. In $B = QR$, $|r_{jj}|$ is the distance from column $j$ to the span of the earlier columns, so $r_{jj} \ne 0$ for $j \le n - 1$. Since $\mu$ is an eigenvalue, $0 = |\det B| = \prod_j|r_{jj}|$, which forces $r_{nn} = 0$: the last row of $R$ is zero. Row $n$ of $RQ$ is $\sum_lR_{nl}Q_{l,:} = 0$, so the last row of $RQ + \mu I$ is $(0, \dots, 0, \mu)$. ∎ In floating point the shift is only nearly an eigenvalue, so $r_{nn}$ is tiny rather than zero, and the last subdiagonal entry becomes tiny: this is why Wilkinson-shifted QR deflates after about two steps per eigenvalue (§3.8(c)).

(c) The columns of $A = \begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}$ are $\mathbf e_2, \mathbf e_1$, already orthonormal, so $Q = A$ and $R = I$, and $RQ = A$: a fixed point. The Rayleigh shift is $\mu = a_{22} = 0$, which gives the same step, so it stalls too; $0$ is equidistant from the eigenvalues $\pm1$, and inverse iteration with that shift has ratio $1$. The Wilkinson shift (3.19) has $a = c = 0$, $b = 1$, $\delta = 0$, $\operatorname{sign}(0) = 1$:
$$\mu = 0 - \frac{1\cdot1^2}{0 + \sqrt{0 + 1}} = -1,$$
an exact eigenvalue. By (b), one step deflates: $A + I = \begin{pmatrix}1 & 1\\1 & 1\end{pmatrix}$ has $r_{22} = 0$, and the project confirms $A_1 = \operatorname{diag}(1, -1)$.

**B7.** (a) From $x_0 = 0$: $r_0 = (1, 2)$, $Ar_0 = (6, 7)$, $\alpha_0 = 5/20 = 1/4$, $x_1 = (1/4, 1/2)$, $r_1 = (-1/2, 1/4)$, the same first step as CG in Example 4.5. Next, $Ar_1 = (-2 + 1/4,\ -1/2 + 3/4) = (-7/4, 1/4)$, $r_1^\top r_1 = 5/16$, $r_1^\top Ar_1 = 7/8 + 1/16 = 15/16$, so $\alpha_1 = 1/3$ and
$$x_2 = \Big(\frac14 - \frac16,\ \frac12 + \frac1{12}\Big) = \Big(\frac1{12}, \frac7{12}\Big), \qquad r_2 = r_1 - \tfrac13Ar_1 = \Big(\frac1{12}, \frac16\Big).$$

(b) $r_2 = r_0/12$. Since $r = Ae$, also $e_2 = e_0/12$. The step length $\alpha = r^\top r/r^\top Ar$ does not change when $r$ is scaled, so from $x_2$ the method repeats its first two steps exactly, scaled by $1/12$: $e_3 = e_1/12$, $e_4 = e_0/144$, and so on. Because $e_0 \ne 0$, no $e_k$ is ever zero. In two dimensions steepest descent bounces between two fixed directions forever. CG, which replaced the second direction $r_1$ by the conjugate direction $p_1$, finished at step 2.

(c) $Ae_0 = r_0 = b$, so $\|e_0\|_A^2 = e_0^\top b = \tfrac1{11} + \tfrac{14}{11} = \tfrac{15}{11}$. Next, $e_1 = x_* - x_1 = (-7/44, 6/44)$ and $Ae_1 = r_1$, so $\|e_1\|_A^2 = \tfrac{7}{88} + \tfrac{6}{176} = \tfrac{20}{176} = \tfrac5{44}$. Finally $\|e_2\|_A^2 = \|e_0\|_A^2/144 = \tfrac{5}{528}$. The ratios $\|e_1\|_A/\|e_0\|_A = \sqrt{(5/44)/(15/11)} = 1/\sqrt{12}$ and $\|e_2\|_A/\|e_1\|_A = \sqrt{44/528} = 1/\sqrt{12}$ are both $0.2887$. The eigenvalues of $A$ are $(7 \pm \sqrt5)/2$, so $\kappa = (7 + \sqrt5)/(7 - \sqrt5) = 1.9387$ and the Kantorovich factor is
$$\frac{\kappa - 1}{\kappa + 1} = \frac{2\sqrt5}{14} = \frac{\sqrt5}{7} = 0.3194.$$
The observed $0.2887$ is below the bound, which is attained only from the worst starting vector.

**B8.** (a) In `householder_tridiagonalize` the reflector loop is already correct for nonsymmetric $A$. The left reflector acts on rows $k+1..n$ and needs only columns $k..n$, because columns before $k$ are zero below row $k+1$ after earlier steps; the right reflector acts on columns $k+1..n$ of every row. The two lines that assume symmetry are at the end. `T = np.triu(np.tril(T, 1), -1)` cuts everything above the first superdiagonal, which is zero only for symmetric $A$: replace it by `np.triu(H, -1)`, which zeros only the rounding-level entries below the subdiagonal. And `0.5 * (T + T.T)` must go. The reference `householder_hessenberg` in Part 7 does exactly this.

(b) The project's run: $\|QHQ^\top - A\|/\|A\| = 2.0\times10^{-16}$, $Q^\top Q = I$, and the largest entry below the subdiagonal over all 150 QR iterates is below $10^{-12}$ relative (exactly $0$ in the reference run), so Hessenberg form is preserved as Theorem 3.17 says. The diagonal of $A_{150}$ matches $(6, 5, 4, 3, 2, 1)$ to $2.7\times10^{-13}$, consistent with the slowest ratio $\lambda_2/\lambda_1 = 5/6$: $(5/6)^{150} = 1.3\times10^{-12}$. (For nonsymmetric $A$ there is no Rayleigh-quotient squaring, so the diagonal error tracks the subdiagonal entries rather than their squares.)

(c) A Hessenberg matrix has $n - 1$ subdiagonal entries, so its QR factorization takes $n - 1$ Givens rotations (3.18). Each touches two rows of length at most $n$, about $6n$ flops, so forming $R$ costs $O(n^2)$; applying the same rotations to the columns to form $RQ$ costs $O(n^2)$ more. A full Householder QR costs $\tfrac43n^3$. With about two steps per eigenvalue, that is the difference between $O(n^3)$ and $O(n^4)$ for the whole eigenvalue computation.

## C. Challenge

**C1.** (a) $\mathbf 1\mathbf 1^\top$ has eigenvalue $3$ on $\mathbf 1$ and $0$ on the plane $\mathbf 1^\perp$. So $A = (1 - a)I + a\mathbf 1\mathbf 1^\top$ has eigenvalues $1 + 2a$ (once) and $1 - a$ (twice). $A$ is symmetric, so it is SPD iff both are positive: $-\tfrac12 < a < 1$.

(b) $D = I$, so $B_J = I - A = a(I - \mathbf 1\mathbf 1^\top)$, with eigenvalues $a(1 - 3) = -2a$ and $a$ (twice). So $\rho(B_J) = 2|a|$, and Jacobi converges iff $|a| < \tfrac12$. For $\tfrac12 \le a < 1$ the matrix is SPD and Jacobi diverges. At $a = 0.8$: eigenvalues of $A$ are $0.2, 0.2, 2.6$, and $\rho(B_J) = 1.6$.

(c) Write $A = L + D + L^\top$, $M = D + L$, $B = I - M^{-1}A$. Let $Bx = \lambda x$ with $x \ne 0$ (complex in general), and set $y = M^{-1}Ax$, so $Bx = x - y$ and $Ax = My$.
1. $\lambda \ne 1$: otherwise $y = 0$, so $Ax = 0$ and $x = 0$. Hence $y = (1 - \lambda)x \ne 0$.
2. Expand, using $A^\top = A$: $(x - y)^*A(x - y) = x^*Ax - y^*Ax - x^*Ay + y^*Ay$.
3. $y^*Ax = y^*My$, and $x^*Ay = (Ax)^*y = (My)^*y = y^*M^\top y$. So $(Bx)^*A(Bx) = x^*Ax - y^*(M + M^\top - A)y$.
4. $M + M^\top - A = D + L + D + L^\top - (L + D + L^\top) = D$.
5. The left side is $|\lambda|^2x^*Ax$ and $y = (1 - \lambda)x$, so $(1 - |\lambda|^2)\,x^*Ax = |1 - \lambda|^2\,x^*Dx$.
6. $x^*Ax > 0$ because $A$ is SPD; $D$ has positive diagonal $d_{ii} = \mathbf e_i^\top A\mathbf e_i > 0$, so $x^*Dx > 0$; and $|1 - \lambda| > 0$ by step 1. So $1 - |\lambda|^2 > 0$, i.e. $|\lambda| < 1$ for every eigenvalue, and Gauss–Seidel converges by Theorem 3.19. ∎

(d) The project: $\rho(B_J) = 1.6000$ (diverges), $\rho(B_{GS}) = 0.7155$ (converges). Jacobi's convergence needs more than SPD (diagonal dominance is enough); Gauss–Seidel's needs only SPD.

**C2.** (a) From $\cos((k+1)\theta) + \cos((k-1)\theta) = 2\cos\theta\cos k\theta$, the functions $T_k(\cos\theta) = \cos k\theta$ satisfy $T_{k+1}(x) = 2xT_k(x) - T_{k-1}(x)$ with $T_0 = 1$, $T_1 = x$. By induction each $T_k$ is a polynomial of degree $k$ with leading coefficient $2^{k-1}$ ($k \ge 1$), and $|T_k(x)| \le 1$ on $[-1, 1]$.

(b) For $x \ge 1$ put $x = \cosh t$ with $t \ge 0$. The identity $\cosh((k+1)t) + \cosh((k-1)t) = 2\cosh t\cosh kt$ shows that $\cosh(kt)$ satisfies the same recurrence with the same $T_0 = 1$, $T_1 = x$, so $T_k(\cosh t) = \cosh kt$. Since $\sinh t = \sqrt{x^2 - 1}$, $e^{\pm t} = \cosh t \pm \sinh t = x \pm \sqrt{x^2 - 1}$, and $\cosh kt = \tfrac12(e^{kt} + e^{-kt})$ is the stated formula.

(c) Let $p^*(\lambda) = T_k(\ell(\lambda))/T_k(\gamma)$, with $T_k(\gamma) > 0$ by (b). At the $k + 1$ points $\lambda_j = \ell^{-1}(\cos(j\pi/k))$, $j = 0..k$, all in $[\lambda_{\min}, \lambda_{\max}]$, $p^*$ takes the values $(-1)^j/T_k(\gamma)$, the extremes of $\pm\max|p^*|$ with alternating signs. Suppose $q$ has degree at most $k$, $q(0) = 1$, and $\max|q| < 1/T_k(\gamma)$ on the interval. At each $\lambda_j$, $|q(\lambda_j)| < |p^*(\lambda_j)|$, so $p^* - q$ has the sign $(-1)^j$. Between consecutive $\lambda_j$ it changes sign, giving $k$ zeros inside the interval. Also $(p^* - q)(0) = 0$, and $0$ is outside the interval: $k + 1$ distinct zeros of a polynomial of degree at most $k$. So $p^* - q \equiv 0$, contradicting $|q(\lambda_0)| < |p^*(\lambda_0)|$. Hence no admissible polynomial is smaller than $p^*$ on the interval. ∎

(d) $\operatorname{arccosh}\gamma = \ln(\gamma + \sqrt{\gamma^2 - 1})$, and the sketch of Theorem 3.24 showed $\gamma + \sqrt{\gamma^2 - 1} = (\sqrt\kappa + 1)/(\sqrt\kappa - 1)$. So $\operatorname{arccosh}\gamma = \ln\frac{\sqrt\kappa + 1}{\sqrt\kappa - 1} = -\ln\rho = 0.15016$ for $\kappa = 178.1$. Then $1/T_k(\gamma) = 1/\cosh(0.15016k) \le 10^{-8}$ requires $0.15016k \ge \operatorname{arccosh}(10^8) = 19.114$, i.e. $k \ge 127.3$: $k = 128$. The bound $2\rho^k = 2e^{-0.15016k} \le 10^{-8}$ requires $0.15016k \ge \ln(2\times10^8) = 19.114$: also $k = 128$. They agree because $\operatorname{arccosh}y = \ln(y + \sqrt{y^2 - 1}) = \ln(2y) - \tfrac{1}{4y^2} + O(y^{-4})$, and at $y = 10^8$ the correction is $2.5\times10^{-17}$. Dropping the $e^{-kt}$ term in step (b) costs nothing. The gap between (3.30) and practice in Example 4.6 ($63$ measured iterations on the residual, against the guaranteed $145$ for the residual and $128$ for the $A$-norm error) comes from replacing the actual eigenvalues by the whole interval, not from this step.

**C3.** The project's results, mean of $\|A - QQ^\top A\|_2/\sigma_{11}$ over 5 trials (without / with re-orthonormalization): $q = 0$: $1.31$ / $1.31$; $q = 4$: $0.57$ / $0.57$; $q = 8$: $1.07$ / $0.54$; $q = 12$: $2.14$ / $0.53$. (These differ slightly from the table in §3.13, which averages 10 different trials.)

Explanation. The component of a column of $(AA^\top)^qA\Omega$ along $u_j$ is $\sigma_j^{2q+1}$ times a Gaussian. Each floating-point product adds rounding errors of relative size $\varepsilon$ compared with the *largest* component, which is of order $\sigma_1^{2q+1}$, and those errors point in every direction. So direction $j$ survives only while $(\sigma_1/\sigma_j)^{2q+1} < 1/\varepsilon \approx 9\times10^{15}$; below that it is rounding noise and the final QR cannot recover it. Here $\sigma_j = 1/j$, so direction $j$ survives iff $j < (1/\varepsilon)^{1/(2q+1)}$.

- $q = 4$: even the last of the $k + p = 20$ directions has $(\sigma_1/\sigma_{20})^9 = 20^9 = 5\times10^{11} < 9\times10^{15}$. Nothing is lost, and the two versions agree.
- $q = 8$: $20^{17} = 1\times10^{22}$, far above $1/\varepsilon$. The threshold is $(1/\varepsilon)^{1/17} = 8.68$: about 8 or 9 directions survive, so $Q$ behaves like a rank-9 basis and the error should be near $\sigma_9$ to $\sigma_{10}$, a ratio of $1.10$ to $1.22$ against $\sigma_{11}$. Measured: $1.07$, close to the lower end.
- $q = 12$: threshold $(1/\varepsilon)^{1/25} = 4.35$, about 4 directions, error near $\sigma_5$, ratio $11/5 = 2.2$. Measured: $2.14$.

More power iterations made the unorthonormalized version *worse*, exactly as power iteration on many vectors without re-orthonormalization collapses them onto the dominant eigenvector (the reason the QR algorithm of §3.8 orthonormalizes at every step). Re-orthonormalizing after each product with $A$ or $A^\top$ limits the amplification within one product to $\sigma_1/\sigma_{20} = 20$, far below $1/\varepsilon$, and the error keeps improving with $q$.
