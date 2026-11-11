# LA-14 Solutions

## A. Warm-up

**A1.** We need $\tfrac12\mathbf{x}^\top A\mathbf{x} = 3x_1^2 - 2x_1x_2 + x_2^2$ with $A$ symmetric. The diagonal entries are twice the squared-term coefficients and the off-diagonal entries are the cross coefficient (since $\tfrac12(a_{12} + a_{21})x_1x_2 = a_{12}x_1x_2$):
$$A = \begin{pmatrix}6 & -2\\-2 & 2\end{pmatrix}, \qquad \mathbf{b} = \begin{pmatrix}5\\0\end{pmatrix}.$$
Check: $\tfrac12(6x_1^2 - 4x_1x_2 + 2x_2^2) = 3x_1^2 - 2x_1x_2 + x_2^2$. ✓

By (3.22) and (3.9), $\nabla f = A\mathbf{x} + \mathbf{b}$ and $H = A$. At $(1, 2)^\top$: $A\mathbf{x} = (6 - 4,\; -2 + 4)^\top = (2, 2)^\top$, so $\nabla f = (7, 2)^\top$. Direct partials: $\partial f/\partial x_1 = 6x_1 - 2x_2 + 5 = 6 - 4 + 5 = 7$; $\partial f/\partial x_2 = -2x_1 + 2x_2 = -2 + 4 = 2$. ✓

Eigenvalues of $H$: trace $8$, determinant $12 - 4 = 8$, so $\lambda = 4 \pm \sqrt{16 - 8} = 4 \pm 2\sqrt2 \approx 6.83, 1.17$. Both positive: a bowl (strictly convex). Its minimizer solves $A\mathbf{x} = -\mathbf{b}$: $A^{-1} = \tfrac18\begin{pmatrix}2 & 2\\2 & 6\end{pmatrix}$, $\mathbf{x}^* = -\tfrac18(10, 10)^\top = (-1.25, -1.25)^\top$. Check: $6(-1.25) - 2(-1.25) + 5 = 0$ and $-2(-1.25) + 2(-1.25) = 0$. ✓

**A2.** $J(x, y) = \begin{pmatrix}\partial_x(x^2y) & \partial_y(x^2y)\\\partial_x(\sin x + y) & \partial_y(\sin x + y)\end{pmatrix} = \begin{pmatrix}2xy & x^2\\\cos x & 1\end{pmatrix}$. At $(1, 2)$: $J = \begin{pmatrix}4 & 1\\0.5403 & 1\end{pmatrix}$.

$\mathbf{F}(1, 2) = (2,\; \sin 1 + 2)^\top = (2, 2.8415)^\top$. $J\mathbf{h} = (0.4 - 0.1,\; 0.05403 - 0.1)^\top = (0.3, -0.0460)^\top$. Linear approximation: $(2.3, 2.7955)^\top$.

Exact: $\mathbf{F}(1.1, 1.9) = (1.21 \times 1.9,\; \sin 1.1 + 1.9)^\top = (2.299, 0.8912 + 1.9)^\top = (2.299, 2.7912)^\top$.

Error: $(-0.0010, -0.0043)^\top$, of size about $0.0044$. Compare $\|\mathbf{h}\| = 0.1414$ and $\|\mathbf{h}\|^2 = 0.02$: the error is $0.22\|\mathbf{h}\|^2$, far smaller than $\|\mathbf{h}\|$, as (3.1) requires. In fact it is the second-order Taylor term (3.22′): for $F_2 = \sin x + y$ the Hessian is $\operatorname{diag}(-\sin x, 0)$, so $\tfrac12\mathbf{h}^\top H\mathbf{h} = \tfrac12(-0.8415)(0.01) = -0.0042$ ✓; for $F_1 = x^2y$ the Hessian is $\begin{pmatrix}2y & 2x\\2x & 0\end{pmatrix} = \begin{pmatrix}4 & 2\\2 & 0\end{pmatrix}$ and $\tfrac12\mathbf{h}^\top H\mathbf{h} = \tfrac12(0.04 - 0.04) = 0$, so the $-0.001$ is the cubic term $h_1^2h_2 = (0.01)(-0.1)$ exactly.

**A3.** (a) $\operatorname{tr}(AB) = \sum_i(AB)_{ii} = \sum_i\sum_j A_{ij}B_{ji}$. Swap the order of summation: $= \sum_j\sum_i B_{ji}A_{ij} = \sum_j(BA)_{jj} = \operatorname{tr}(BA)$. ∎ Note $AB$ is $m \times m$ and $BA$ is $n \times n$; the traces agree anyway.

(b) $(\mathbf{v}\mathbf{u}^\top)_{ii} = v_iu_i$, so $\operatorname{tr}(\mathbf{v}\mathbf{u}^\top) = \sum_iu_iv_i = \mathbf{u}^\top\mathbf{v}$. ∎

(c) $\big(A^\top(B \odot C)\big)_{jj} = \sum_i(A^\top)_{ji}(B \odot C)_{ij} = \sum_iA_{ij}B_{ij}C_{ij}$, so the trace is $\sum_{ij}A_{ij}B_{ij}C_{ij}$. Likewise $\big((A \odot B)^\top C\big)_{jj} = \sum_iA_{ij}B_{ij}C_{ij}$. Same triple sum. ∎

Cyclic: apply (a) with the two factors $AB$ and $C$: $\operatorname{tr}((AB)C) = \operatorname{tr}(C(AB))$. Counterexample to swapping non-cyclically: $A = \begin{pmatrix}0 & 1\\0 & 0\end{pmatrix}$, $B = \begin{pmatrix}0 & 0\\1 & 0\end{pmatrix}$, $C = \begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}$. $AB = \begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}$, $ABC = \begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}$, trace $1$. $AC = \begin{pmatrix}0 & 0\\0 & 0\end{pmatrix}$, so $\operatorname{tr}(ACB) = 0$. $1 \ne 0$.

**A4.** $e^1 + e^2 + e^3 = 2.7183 + 7.3891 + 20.0855 = 30.1929$, so $\mathbf{s} = (0.0900, 0.2447, 0.6652)^\top$.
$$\mathbf{s}\mathbf{s}^\top = \begin{pmatrix}0.0081 & 0.0220 & 0.0599\\0.0220 & 0.0599 & 0.1628\\0.0599 & 0.1628 & 0.4425\end{pmatrix}, \qquad J = \operatorname{diag}(\mathbf{s}) - \mathbf{s}\mathbf{s}^\top = \begin{pmatrix}0.0819 & -0.0220 & -0.0599\\-0.0220 & 0.1848 & -0.1628\\-0.0599 & -0.1628 & 0.2227\end{pmatrix}.$$
Row sums: $0.0819 - 0.0220 - 0.0599 = 0.0000$; $-0.0220 + 0.1848 - 0.1628 = 0.0000$; $-0.0599 - 0.1628 + 0.2227 = 0.0000$ (to rounding). $J\mathbf{1}$ is the vector of row sums, so $J\mathbf{1} = \mathbf{0}$: algebraically, $J\mathbf{1} = \mathbf{s} - \mathbf{s}(\mathbf{s}^\top\mathbf{1}) = \mathbf{s} - \mathbf{s} = \mathbf{0}$ because $\sum_ks_k = 1$.

Cross-entropy with $\mathbf{y} = \mathbf{e}_1$: by (3.21), $\nabla_{\mathbf{z}}L = \mathbf{s} - \mathbf{e}_1 = (-0.9100, 0.2447, 0.6652)^\top$. The correct logit is pushed up (negative gradient) by $0.91$, the wrong ones pushed down in proportion to their probabilities.

**A5.** The central difference of $e^x$ at $0$ is $\sinh(h)/h = 1 + h^2/6 + h^4/120 + \cdots$, so the truncation error is $h^2/6$ exactly to leading order ($f''' = 1$). The rounding estimate is $\epsilon|f|/h = 2.2 \times 10^{-16}/h$. Typical computed errors:

| $h$ | error observed | $h^2/6$ (truncation) | $\epsilon/h$ (rounding) | dominant |
|---|---|---|---|---|
| $10^{-1}$ | $1.67 \times 10^{-3}$ | $1.67 \times 10^{-3}$ | $2 \times 10^{-15}$ | truncation |
| $10^{-3}$ | $1.67 \times 10^{-7}$ | $1.67 \times 10^{-7}$ | $2 \times 10^{-13}$ | truncation |
| $10^{-6}$ | $\sim 3 \times 10^{-11}$ | $1.7 \times 10^{-13}$ | $2 \times 10^{-10}$ | rounding |
| $10^{-10}$ | $\sim 8 \times 10^{-8}$ | $1.7 \times 10^{-21}$ | $2 \times 10^{-6}$ | rounding |

The first two rows match the truncation formula to three digits. In the last two the error is rounding noise: its exact value (and even its sign) depends on the machine, but its size is within the $\epsilon/h$ estimate. $h = 10^{-6}$ is closest to $\epsilon^{1/3} \approx 6 \times 10^{-6}$ and is indeed the most accurate of the four, by three orders of magnitude over $h = 10^{-10}$. "Smaller $h$ is more accurate" is false below $\epsilon^{1/3}$.

## B. Core

**B1.** (a) Let $F(A, B) = AB$ and perturb both arguments: $(A + dA)(B + dB) - AB = (dA)B + A(dB) + (dA)(dB)$. The first two terms are linear in the increment $(dA, dB)$. For the third, entry $(i, k)$ of $(dA)(dB)$ is $\sum_j dA_{ij}dB_{jk}$, the dot product of row $i$ of $dA$ with column $k$ of $dB$, so by Cauchy–Schwarz (LA-1, Theorem 3.11) $|((dA)(dB))_{ik}|^2 \le \|\text{row}_i(dA)\|^2\|\text{col}_k(dB)\|^2$. Summing over $i$ and $k$: $\|(dA)(dB)\|_F^2 \le \big(\sum_i\|\text{row}_i(dA)\|^2\big)\big(\sum_k\|\text{col}_k(dB)\|^2\big) = \|dA\|_F^2\|dB\|_F^2$. With $\rho = \sqrt{\|dA\|_F^2 + \|dB\|_F^2}$ the size of the increment, $\|dA\|_F\|dB\|_F \le \tfrac12\rho^2 = o(\rho)$. So $F(A + dA, B + dB) = F(A, B) + [(dA)B + A(dB)] + o(\rho)$, and by Definition 3.7 the bracket is $dF$. ∎

(b) $AA^{-1} = I$ holds for all invertible $A$ near a given one, so both sides have the same differential. The right side is constant: $d(I) = 0$. The left, by (3.5): $(dA)A^{-1} + A\,d(A^{-1}) = 0$. Multiply by $A^{-1}$ on the left: $d(A^{-1}) = -A^{-1}(dA)A^{-1}$. ∎

(c) $d(A^2) = d(A \cdot A) = (dA)A + A(dA)$. $d(A^3) = d(A^2 \cdot A) = d(A^2)A + A^2\,dA = (dA)A^2 + A(dA)A + A^2(dA)$. With $A = \begin{pmatrix}0 & 1\\0 & 0\end{pmatrix}$, $dA = \begin{pmatrix}0 & 0\\1 & 0\end{pmatrix}$: $(dA)A = \begin{pmatrix}0 & 0\\0 & 1\end{pmatrix}$, $A(dA) = \begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}$, so $d(A^2) = I$. But $2A\,dA = \begin{pmatrix}2 & 0\\0 & 0\end{pmatrix} \ne I$. The scalar rule $d(a^2) = 2a\,da$ needs $a$ and $da$ to commute.

(d) $d\operatorname{tr}(A^2) = \operatorname{tr}\big((dA)A + A(dA)\big) = 2\operatorname{tr}(A\,dA) = 2\operatorname{tr}\big((A^\top)^\top dA\big)$, so $\nabla_A\operatorname{tr}(A^2) = 2A^\top$. (Even though $d(A^2) \ne 2A\,dA$, the trace makes the two terms equal by the cyclic rule.) $d\operatorname{tr}(A^3) = \operatorname{tr}\big((dA)A^2 + A(dA)A + A^2(dA)\big) = 3\operatorname{tr}(A^2\,dA)$ (each term cycles to $\operatorname{tr}(A^2\,dA)$), so $\nabla_A\operatorname{tr}(A^3) = 3(A^2)^\top = 3(A^\top)^2$.

**B2.** (a) With $\mathbf{r} = A\mathbf{x} - \mathbf{b}$, $d\mathbf{r} = A\,d\mathbf{x}$, and by (3.11)–(3.12):
$$df = 2\mathbf{r}^\top A\,d\mathbf{x} + 2\lambda\mathbf{x}^\top d\mathbf{x} = \big(2A^\top(A\mathbf{x} - \mathbf{b}) + 2\lambda\mathbf{x}\big)^\top d\mathbf{x}.$$
So $\nabla f = 2A^\top(A\mathbf{x} - \mathbf{b}) + 2\lambda\mathbf{x}$, and $\nabla f = \mathbf{0}$ iff $A^\top A\mathbf{x} + \lambda\mathbf{x} = A^\top\mathbf{b}$, i.e. $(A^\top A + \lambda I)\mathbf{x} = A^\top\mathbf{b}$.

(b) $\nabla f$ is linear in $\mathbf{x}$: $d(\nabla f) = 2(A^\top A + \lambda I)\,d\mathbf{x}$, so $H = 2(A^\top A + \lambda I)$. For $\mathbf{v} \ne \mathbf{0}$: $\mathbf{v}^\top H\mathbf{v} = 2(\|A\mathbf{v}\|^2 + \lambda\|\mathbf{v}\|^2) \ge 2\lambda\|\mathbf{v}\|^2 > 0$. So $H \succ 0$ regardless of whether $A^\top A$ is singular ($A^\top A \succeq 0$ always, by LA-8, since $\mathbf{v}^\top A^\top A\mathbf{v} = \|A\mathbf{v}\|^2$). A quadratic with $H \succ 0$ satisfies $f(\mathbf{x}) = f(\mathbf{x}^*) + \tfrac12(\mathbf{x} - \mathbf{x}^*)^\top H(\mathbf{x} - \mathbf{x}^*)$ at its stationary point $\mathbf{x}^*$ (expand; the linear term vanishes because $\nabla f(\mathbf{x}^*) = \mathbf{0}$), so $f(\mathbf{x}) > f(\mathbf{x}^*)$ for every $\mathbf{x} \ne \mathbf{x}^*$: the minimizer is unique, and $A^\top A + \lambda I$ is invertible, so it exists.

(c) $A^\top A = \begin{pmatrix}2 & 2\\2 & 2\end{pmatrix}$, $A^\top A + I = \begin{pmatrix}3 & 2\\2 & 3\end{pmatrix}$, $A^\top\mathbf{b} = (2, 2)^\top$. Determinant $9 - 4 = 5$, inverse $\tfrac15\begin{pmatrix}3 & -2\\-2 & 3\end{pmatrix}$, so $\mathbf{x} = \tfrac15(6 - 4,\; -4 + 6)^\top = (0.4, 0.4)^\top$.

For general $\lambda$: $A^\top A = 2\mathbf{1}\mathbf{1}^\top$ and $A^\top\mathbf{b} = 2\cdot\mathbf{1}$. Try $\mathbf{x} = c\mathbf{1}$: $(2\mathbf{1}\mathbf{1}^\top + \lambda I)c\mathbf{1} = (4c + \lambda c)\mathbf{1} = 2\cdot\mathbf{1}$ gives $c = 2/(4 + \lambda)$. ($\lambda = 1$: $c = 0.4$ ✓.) As $\lambda \to 0$, $c \to \tfrac12$: $\mathbf{x} \to (0.5, 0.5)^\top$.

Pseudoinverse: $A = \mathbf{1}\mathbf{1}^\top$ has rank one with SVD $\sigma_1 = 2$, $\mathbf{u}_1 = \mathbf{v}_1 = \mathbf{1}/\sqrt2$, so $A^+ = \mathbf{v}_1\mathbf{u}_1^\top/\sigma_1 = \mathbf{1}\mathbf{1}^\top/4$ (LA-9) and $A^+\mathbf{b} = \mathbf{1}(\mathbf{1}^\top\mathbf{b})/4 = \tfrac24\mathbf{1} = (0.5, 0.5)^\top$. ✓ The unregularized problem has infinitely many least-squares solutions ($x_1 + x_2 = 1$: any of them projects $\mathbf{b}$ onto $\operatorname{span}\{\mathbf{1}\}$), and ridge picks out the minimum-norm one in the limit, which is what $A^+$ returns.

**B3.** Let $F = AX - B$ ($m \times n$). By (3.5), $dF = (dA)X + A(dX)$. By (3.18) and rule (h), $d\|F\|_F^2 = 2\operatorname{tr}(F^\top dF)$.

*With respect to $A$* ($dX = 0$): $2\operatorname{tr}(F^\top(dA)X) = 2\operatorname{tr}(XF^\top dA)$ (cyclic) $= 2\operatorname{tr}\big((FX^\top)^\top dA\big)$ since $(FX^\top)^\top = XF^\top$. So $\nabla_A\|AX - B\|_F^2 = 2(AX - B)X^\top$. Shape: $(m \times n)(n \times p) = m \times p$, the shape of $A$. ✓

*With respect to $X$* ($dA = 0$): $2\operatorname{tr}(F^\top A\,dX) = 2\operatorname{tr}\big((A^\top F)^\top dX\big)$. So $\nabla_X\|AX - B\|_F^2 = 2A^\top(AX - B)$. Shape: $(p \times m)(m \times n) = p \times n$. ✓

Setting the second to zero: $A^\top AX = A^\top B$. Column $j$ of this equation is $A^\top A\mathbf{x}_j = A^\top\mathbf{b}_j$, the LA-5 normal equations for $\min\|A\mathbf{x}_j - \mathbf{b}_j\|^2$; this is right because $\|AX - B\|_F^2 = \sum_j\|A\mathbf{x}_j - \mathbf{b}_j\|^2$ and the columns of $X$ do not interact, so minimizing the sum is minimizing each term.

Code to add to the `identities` dictionary in Part 2 (the surrounding random matrices `A`, `P`, `M` have compatible shapes: use `M` ($4 \times 5$) as $A$, `P` ($5 \times 4$) as $X$, and a fresh $4 \times 4$ target):
```python
Bt = rng.standard_normal((4, 4))
identities["(B3) grad_A ||A X - B||_F^2 = 2(AX-B)X^T"] = (
    lambda W: np.sum((W @ P - Bt) ** 2), lambda W: 2 * (W @ P - Bt) @ P.T, M)
identities["(B3) grad_X ||A X - B||_F^2 = 2A^T(AX-B)"] = (
    lambda W: np.sum((M @ W - Bt) ** 2), lambda W: 2 * M.T @ (M @ W - Bt), P)
```
Both pass `grad_check` with relative error around $10^{-11}$.

**B4.** (a) $g = \tfrac12\mathbf{x}^\top\begin{pmatrix}2 & 4\\4 & 2\end{pmatrix}\mathbf{x}$, so by (3.22) $\nabla g = (2x + 4y,\; 4x + 2y)^\top$ and $H = \begin{pmatrix}2 & 4\\4 & 2\end{pmatrix}$. Eigenvalues: trace $4$, determinant $-12$, $\lambda = 2 \pm 4 = 6, -2$. Eigenvectors: for $6$, $(H - 6I)\mathbf{v} = 0$ gives $-4v_1 + 4v_2 = 0$, $\mathbf{q}_1 = (1, 1)^\top/\sqrt2$; for $-2$, $\mathbf{q}_2 = (1, -1)^\top/\sqrt2$. Mixed signs: the origin is a saddle (Theorem 3.14). Fastest increase along $\mathbf{q}_1$: $g(t, t) = t^2 + 4t^2 + t^2 = 6t^2 = \lambda_1t^2$. Decrease along $\mathbf{q}_2$: $g(t, -t) = t^2 - 4t^2 + t^2 = -2t^2 = \lambda_2t^2$.

(b) $\nabla g(1, 1) = (6, 6)^\top$. $H^{-1} = \frac{1}{4 - 16}\begin{pmatrix}2 & -4\\-4 & 2\end{pmatrix} = -\frac1{12}\begin{pmatrix}2 & -4\\-4 & 2\end{pmatrix}$. $H^{-1}\nabla g = -\frac1{12}(12 - 24,\; -24 + 12)^\top = -\frac1{12}(-12, -12)^\top = (1, 1)^\top$. Newton step (3.23): $\mathbf{x}_1 = (1, 1)^\top - (1, 1)^\top = (0, 0)^\top$. It lands exactly on the saddle (one step, because $g$ is quadratic). Newton's method solves $\nabla g = \mathbf{0}$; it does not know or care whether the solution is a minimum. When $H$ is indefinite, the Newton direction need not even be a descent direction. This is why ML-2's practical second-order methods modify $H$ (add $\lambda I$, use $J_r^\top J_r \succeq 0$ as in Gauss–Newton, or use a trust region).

(c) $\nabla h = (3x^2 - 3y^2,\; -6xy)^\top$, which is $\mathbf{0}$ at the origin. $H = \begin{pmatrix}6x & -6y\\-6y & -6x\end{pmatrix}$, which is the zero matrix at the origin, so Theorem 3.14 is silent. On the $x$-axis, $h(t, 0) = t^3$, which is positive for $t > 0$ and negative for $t < 0$: every neighborhood of the origin contains points where $h$ is above and below $h(0, 0) = 0$, so the origin is neither a minimum nor a maximum. (It is a "monkey saddle": three valleys and three ridges, at $60°$ spacing, since $h = \operatorname{Re}(x + iy)^3$.)

**B5.** (a) A $p \times q$ times $q \times r$ product costs $2pqr$ flops.
(i) Right to left: the $k - 1$ matrices $J_1, \dots, J_{k-1}$ are $d \times d$, so the $k - 2$ products among them each cost $2d^3$; the final product $J_k \cdot (d \times d)$ costs $2d^2$. Total $2(k - 2)d^3 + 2d^2 \sim 2kd^3$.
(ii) Reverse mode (left to right): $J_kJ_{k-1}$ is $(1 \times d)(d \times d)$, $2d^2$; the result is $1 \times d$ again, so each of the $k - 1$ products costs $2d^2$. Total $2(k - 1)d^2 \sim 2kd^2$.
(iii) Forward mode: one JVP $J_k(\cdots(J_1\mathbf{v}))$ is $k - 1$ products of type $(d \times d)(d \times 1)$ plus one $(1 \times d)(d \times 1)$: $2(k - 1)d^2 + 2d$. The gradient needs $d$ of them ($\mathbf{v} = \mathbf{e}_1, \dots, \mathbf{e}_d$): $2(k - 1)d^3 + 2d^2 \sim 2kd^3$.

(b) $d = 1000$, $k = 10$: (i) $2 \cdot 8 \cdot 10^9 + 2 \cdot 10^6 \approx 1.6 \times 10^{10}$; (ii) $2 \cdot 9 \cdot 10^6 = 1.8 \times 10^7$; (iii) $2 \cdot 9 \cdot 10^9 \approx 1.8 \times 10^{10}$. Reverse is a thousand times cheaper than either alternative.

(c) To evaluate $J_\ell$ you need $\mathbf{x}_{\ell-1}$, so reverse mode stores $\mathbf{x}_0, \dots, \mathbf{x}_{k-1}$: $kd$ numbers (and in a real network also the pre-activations, another $kd$). Forward mode carries only the current $(\mathbf{x}_\ell, J_\ell\cdots J_1\mathbf{v})$: $2d$ numbers.

(d) For $f : \mathbb{R} \to \mathbb{R}^d$, $J = J_k\cdots J_1$ with $J_1 \in \mathbb{R}^{d \times 1}$ a column and $J_k \in \mathbb{R}^{d \times d}$. One JVP with the scalar direction $v = 1$ costs $k - 1$ products of type $(d \times d)(d \times 1)$, $2(k - 1)d^2$, and yields the entire Jacobian (which is a single column). Reverse mode would need one VJP per output, $d$ of them, each $2(k - 1)d^2$: $2(k - 1)d^3$. Forward wins by a factor $d$. General rule: the Jacobian of $f : \mathbb{R}^n \to \mathbb{R}^m$ costs $n$ JVPs or $m$ VJPs, each about the cost of one evaluation of $f$; use forward when $n \le m$, reverse when $m < n$. A scalar loss is $m = 1$, the extreme case for reverse.

**B6.** (a) $W_1 \in \mathbb{R}^{h_1 \times d}$, $\mathbf{b}_1, \mathbf{z}_1, \mathbf{a}_1 \in \mathbb{R}^{h_1}$; $W_2 \in \mathbb{R}^{h_2 \times h_1}$, $\mathbf{b}_2, \mathbf{z}_2, \mathbf{a}_2 \in \mathbb{R}^{h_2}$; $W_3 \in \mathbb{R}^{m \times h_2}$, $\mathbf{b}_3, \mathbf{z}_3, \mathbf{y} \in \mathbb{R}^m$.

(b) Every line is (3.5), (3.14), (3.19) and the identification rule, exactly as in §3.7.
1. $dL = (\mathbf{z}_3 - \mathbf{y})^\top d\mathbf{z}_3$. Define $\boldsymbol{\delta}_3 = \mathbf{z}_3 - \mathbf{y} \in \mathbb{R}^m$.
2. $d\mathbf{z}_3 = (dW_3)\mathbf{a}_2 + W_3\,d\mathbf{a}_2 + d\mathbf{b}_3$. So $\nabla_{W_3}L = \boldsymbol{\delta}_3\mathbf{a}_2^\top \in \mathbb{R}^{m \times h_2}$, $\nabla_{\mathbf{b}_3}L = \boldsymbol{\delta}_3 \in \mathbb{R}^m$, and $\nabla_{\mathbf{a}_2}L = W_3^\top\boldsymbol{\delta}_3 \in \mathbb{R}^{h_2}$.
3. $d\mathbf{a}_2 = \sigma'(\mathbf{z}_2) \odot d\mathbf{z}_2$, so $\boldsymbol{\delta}_2 := \nabla_{\mathbf{z}_2}L = (W_3^\top\boldsymbol{\delta}_3) \odot \sigma'(\mathbf{z}_2) \in \mathbb{R}^{h_2}$.
4. $d\mathbf{z}_2 = (dW_2)\mathbf{a}_1 + W_2\,d\mathbf{a}_1 + d\mathbf{b}_2$. So $\nabla_{W_2}L = \boldsymbol{\delta}_2\mathbf{a}_1^\top \in \mathbb{R}^{h_2 \times h_1}$, $\nabla_{\mathbf{b}_2}L = \boldsymbol{\delta}_2 \in \mathbb{R}^{h_2}$, $\nabla_{\mathbf{a}_1}L = W_2^\top\boldsymbol{\delta}_2 \in \mathbb{R}^{h_1}$.
5. $d\mathbf{a}_1 = \sigma'(\mathbf{z}_1) \odot d\mathbf{z}_1$, so $\boldsymbol{\delta}_1 := \nabla_{\mathbf{z}_1}L = (W_2^\top\boldsymbol{\delta}_2) \odot \sigma'(\mathbf{z}_1) \in \mathbb{R}^{h_1}$.
6. $d\mathbf{z}_1 = (dW_1)\mathbf{x} + d\mathbf{b}_1$. So $\nabla_{W_1}L = \boldsymbol{\delta}_1\mathbf{x}^\top \in \mathbb{R}^{h_1 \times d}$, $\nabla_{\mathbf{b}_1}L = \boldsymbol{\delta}_1 \in \mathbb{R}^{h_1}$.

The pattern is the recursion stated after (3.33), applied twice.

(c) Only step 1 changes. With $\mathbf{s} = \operatorname{softmax}(\mathbf{z}_3)$ and $L = -\mathbf{y}^\top\log\mathbf{s}$, (3.21) gives $dL = (\mathbf{s} - \mathbf{y})^\top d\mathbf{z}_3$, so $\boldsymbol{\delta}_3 = \mathbf{s} - \mathbf{y}$. Steps 2–6 used $d\mathbf{z}_3$ only through $\boldsymbol{\delta}_3$, so every other formula is unchanged. (Here $\mathbf{y}$ is a probability vector, typically one-hot.)

**B7.** (a) $\partial_ar_i = e^{-bt_i}$, $\partial_br_i = -at_ie^{-bt_i}$. Second derivatives: $\partial_{aa}r_i = 0$, $\partial_{ab}r_i = -t_ie^{-bt_i}$, $\partial_{bb}r_i = at_i^2e^{-bt_i}$. So
$$\nabla^2r_i = \begin{pmatrix}0 & -t_ie^{-bt_i}\\-t_ie^{-bt_i} & at_i^2e^{-bt_i}\end{pmatrix}.$$

(b) The off-diagonal of $\sum_ir_i\nabla^2r_i$ is $-\sum_ir_it_ie^{-bt_i}$. The second entry of $J_r^\top\mathbf{r}$ is $\sum_i(\partial_br_i)r_i = -a\sum_ir_it_ie^{-bt_i}$. Hence off-diagonal $= (J_r^\top\mathbf{r})_2/a$. At a stationary point $\nabla f = J_r^\top\mathbf{r} = \mathbf{0}$, so the off-diagonal of the dropped term is exactly zero there. (The $(2, 2)$ entry $a\sum_ir_it_i^2e^{-bt_i}$ has no such reason to vanish and in Example 4.5 equals $-0.046$ at the solution.)

(c) At $(a, b) = (1.5, 0.5)$ with $t = (0, 1, 2)$, $y = (2.0, 1.2, 0.8)$: $e^{-bt} = (1, 0.6065, 0.3679)$, $\mathbf{r} = (1.5 - 2,\; 0.9098 - 1.2,\; 0.5518 - 0.8) = (-0.5, -0.2902, -0.2482)$. $t_ie^{-bt_i} = (0, 0.6065, 0.7358)$, $t_i^2e^{-bt_i} = (0, 0.6065, 1.4715)$.

Off-diagonal: $-[(-0.2902)(0.6065) + (-0.2482)(0.7358)] = -[-0.1760 - 0.1826] = 0.3586$. Entry $(2,2)$: $1.5[(-0.2902)(0.6065) + (-0.2482)(1.4715)] = 1.5[-0.1760 - 0.3652] = -0.8118$.
$$\sum_ir_i\nabla^2r_i = \begin{pmatrix}0 & 0.3586\\0.3586 & -0.8118\end{pmatrix}, \qquad J_r^\top J_r = \begin{pmatrix}1.5032 & -0.9578\\-0.9578 & 2.0457\end{pmatrix}, \qquad H = \begin{pmatrix}1.5032 & -0.5992\\-0.5992 & 1.2339\end{pmatrix}.$$
(Check: $(J_r^\top\mathbf{r})_2/a = 0.5379/1.5 = 0.3586$. ✓) The dropped term is about 40% of the $(2, 2)$ entry and a third of the off-diagonal: the approximation is poor at the start. The method still converged because (i) $J_r^\top J_r \succ 0$, so the Gauss–Newton step is a descent direction: $\nabla f^\top\mathbf{h} = -\mathbf{r}^\top J_r(J_r^\top J_r)^{-1}J_r^\top\mathbf{r} < 0$; (ii) the undamped step happened to be short enough to reduce $f$ (from $0.198$ to $0.00136$), which is luck, not theorem; and (iii) after that step the residuals are small, the dropped term is small, and Gauss–Newton converges fast, though linearly: each step is about $0.017$ times the previous one, the spectral radius of $(J_r^\top J_r)^{-1}\sum_ir_i\nabla^2r_i$ at the solution (Example 4.5). Only a zero-residual fit would make it quadratic like Newton. Far from the solution, undamped Gauss–Newton can overshoot (try the same problem from $(1, 1)$: the first step sends $b$ negative). Levenberg–Marquardt's $\lambda I$ (ML-2) is the fix.

**B8.** (a) $\ell(\boldsymbol{\mu}, \Sigma) = \sum_{i=1}^n\ell_i = -\tfrac{np}{2}\log 2\pi - \tfrac n2\log\det\Sigma - \tfrac12\sum_i(\mathbf{x}_i - \boldsymbol{\mu})^\top\Sigma^{-1}(\mathbf{x}_i - \boldsymbol{\mu})$. Gradients add, so with $\mathbf{u}_i = \mathbf{x}_i - \boldsymbol{\mu}$ and $S_{\boldsymbol{\mu}} = \sum_i\mathbf{u}_i\mathbf{u}_i^\top$, (3.38)–(3.39) give
$$\nabla_{\boldsymbol{\mu}}\ell = \Sigma^{-1}\sum_i\mathbf{u}_i = \Sigma^{-1}\Big(\sum_i\mathbf{x}_i - n\boldsymbol{\mu}\Big), \qquad \nabla_\Sigma\ell = \tfrac12\big(\Sigma^{-1}S_{\boldsymbol{\mu}}\Sigma^{-1} - n\Sigma^{-1}\big).$$

(b) $\nabla_{\boldsymbol{\mu}}\ell = \mathbf{0}$: multiply by $\Sigma$ to get $\sum_i\mathbf{x}_i = n\boldsymbol{\mu}$, so $\hat{\boldsymbol{\mu}} = \bar{\mathbf{x}}$, for every $\Sigma$. $\nabla_\Sigma\ell = 0$: $\Sigma^{-1}S\Sigma^{-1} = n\Sigma^{-1}$; multiply by $\Sigma$ on the left and on the right: $S = n\Sigma$, so $\hat\Sigma = S/n$ with $S$ evaluated at $\boldsymbol{\mu} = \bar{\mathbf{x}}$ (substituting the first result into the second).

(c) $\nabla_{\boldsymbol{\mu}}\ell = \Sigma^{-1}(\sum_i\mathbf{x}_i - n\boldsymbol{\mu})$ is affine in $\boldsymbol{\mu}$; $d(\nabla_{\boldsymbol{\mu}}\ell) = -n\Sigma^{-1}d\boldsymbol{\mu}$, so $H_{\boldsymbol{\mu}} = -n\Sigma^{-1}$. $\Sigma \succ 0$ has eigenvalues $\lambda_i > 0$, so $\Sigma^{-1}$ has eigenvalues $1/\lambda_i > 0$ (LA-8) and $\Sigma^{-1} \succ 0$; thus $H_{\boldsymbol{\mu}} \prec 0$. $\ell$ is a concave quadratic in $\boldsymbol{\mu}$ with a negative-definite Hessian, so its stationary point $\bar{\mathbf{x}}$ is the unique global maximum in $\boldsymbol{\mu}$ (same argument as B2(b) with the sign flipped). ($-H_{\boldsymbol{\mu}} = n\Sigma^{-1}$ is the Fisher information for $\boldsymbol{\mu}$, ST-1.)

(d) $\bar{\mathbf{x}} = \big(\tfrac{1 + 3 + 2}{3}, \tfrac{0 + 2 + 4}{3}\big) = (2, 2)^\top$. Deviations: $(-1, -2)^\top$, $(1, 0)^\top$, $(0, 2)^\top$. Outer products: $\begin{pmatrix}1 & 2\\2 & 4\end{pmatrix} + \begin{pmatrix}1 & 0\\0 & 0\end{pmatrix} + \begin{pmatrix}0 & 0\\0 & 4\end{pmatrix} = \begin{pmatrix}2 & 2\\2 & 8\end{pmatrix} = S$. $\hat\Sigma = S/3 = \begin{pmatrix}2/3 & 2/3\\2/3 & 8/3\end{pmatrix}$.

Verification: at $\Sigma = \hat\Sigma$ we have $S = 3\hat\Sigma$, so $\nabla_\Sigma\ell = \tfrac12\big(\hat\Sigma^{-1}(3\hat\Sigma)\hat\Sigma^{-1} - 3\hat\Sigma^{-1}\big) = \tfrac12(3\hat\Sigma^{-1} - 3\hat\Sigma^{-1}) = 0$. ✓ Explicitly, $\det\hat\Sigma = \tfrac{16}{9} - \tfrac49 = \tfrac43$ and $\hat\Sigma^{-1} = \tfrac34\begin{pmatrix}8/3 & -2/3\\-2/3 & 2/3\end{pmatrix} = \begin{pmatrix}2 & -0.5\\-0.5 & 0.5\end{pmatrix}$; $\hat\Sigma^{-1}S\hat\Sigma^{-1} = 3\hat\Sigma^{-1}\hat\Sigma\hat\Sigma^{-1} = 3\hat\Sigma^{-1}$, same conclusion.

## C. Challenge

**C1.** (a) By (3.6), $d(X^{-1}) = -X^{-1}(dX)X^{-1} = A\,(dX)\,B$ with $A = -X^{-1}$, $B = X^{-1}$. By the `vec` identity, $\operatorname{vec}(d(X^{-1})) = (B^\top \otimes A)\operatorname{vec}(dX) = \big(X^{-\top} \otimes (-X^{-1})\big)\operatorname{vec}(dX) = -(X^{-\top} \otimes X^{-1})\operatorname{vec}(dX)$. By Theorem 3.9 this matrix is the Jacobian. ∎

For $X = \begin{pmatrix}2 & 1\\1 & 3\end{pmatrix}$: $X^{-1} = \tfrac15\begin{pmatrix}3 & -1\\-1 & 2\end{pmatrix} = \begin{pmatrix}0.6 & -0.2\\-0.2 & 0.4\end{pmatrix}$, symmetric, so $X^{-\top} = X^{-1}$. The Kronecker product $X^{-\top} \otimes X^{-1}$ is the block matrix $\begin{pmatrix}0.6X^{-1} & -0.2X^{-1}\\-0.2X^{-1} & 0.4X^{-1}\end{pmatrix}$, so
$$-(X^{-\top} \otimes X^{-1}) = \begin{pmatrix}-0.36 & 0.12 & 0.12 & -0.04\\0.12 & -0.24 & -0.04 & 0.08\\0.12 & -0.04 & -0.24 & 0.08\\-0.04 & 0.08 & 0.08 & -0.16\end{pmatrix}.$$
Check of the $(1, 1)$ entry: with column-major `vec`, both the first output and the first input are the $(1,1)$ entries, so this is $\partial(X^{-1})_{11}/\partial x_{11}$. $(X^{-1})_{11} = x_{22}/(x_{11}x_{22} - x_{12}x_{21})$, so $\partial/\partial x_{11} = -x_{22}\cdot x_{22}/\det^2 = -9/25 = -0.36$. ✓

(b) $d\operatorname{tr}(AX^{-1}) = \operatorname{tr}(A\,d(X^{-1})) = -\operatorname{tr}(AX^{-1}(dX)X^{-1}) = -\operatorname{tr}(X^{-1}AX^{-1}\,dX)$ (cyclic). Write this as $\operatorname{tr}(G^\top dX)$ with $G^\top = -X^{-1}AX^{-1}$, so $G = -(X^{-1}AX^{-1})^\top = -X^{-\top}A^\top X^{-\top}$. ∎ With $A = I$: $\nabla_X\operatorname{tr}(X^{-1}) = -X^{-\top}X^{-\top} = -(X^{-1}X^{-1})^\top = -(X^{-2})^\top$. ∎

(c) Add to Part 2 of `project.py`:
```python
X4 = np.array([[2., 1.], [1., 3.]]); Xi4 = np.linalg.inv(X4)
vec_inv = lambda v: np.linalg.inv(v.reshape(2, 2, order="F")).ravel(order="F")
assert rel_error(num_jac(vec_inv, X4.ravel(order="F")), -np.kron(Xi4.T, Xi4)) < 1e-7
assert rel_error(num_grad(lambda W: np.trace(np.linalg.inv(W)), Xmat), -np.linalg.inv(Xmat @ Xmat).T) < 1e-7
```

**C2.** Reference implementation (append to `project.py`, after the `Tensor` class):
```python
def jvp(F, x, v, h=None):
    """J_F(x) v by a central difference along v (one directional derivative)."""
    h = EPS ** (1 / 3) if h is None else h
    return (np.asarray(F(x + h * v)) - np.asarray(F(x - h * v))) / (2 * h)

dN = 1000
Wbig = rng.standard_normal((dN, dN)) / np.sqrt(dN); xbig = rng.standard_normal(dN)
fbig = lambda v: np.sum(np.tanh(Wbig @ v))
t0 = time.perf_counter()
g_fwd = np.array([jvp(fbig, xbig, e) for e in np.eye(dN)])      # (i) 1000 JVPs
t_fwd = time.perf_counter() - t0
t0 = time.perf_counter()
_, (g_rev,) = tensor_grads(lambda xt: (Tensor(Wbig) @ xt).tanh().sum(), [xbig[:, None]])  # (ii) one reverse pass
t_rev = time.perf_counter() - t0
g_closed = Wbig.T @ (1 - np.tanh(Wbig @ xbig) ** 2)                # (iii) eq. (3.19) with a = 1
assert rel_error(g_fwd, g_closed) < 1e-7 and rel_error(g_rev.ravel(), g_closed) < 1e-10
print(f"forward {t_fwd*1e3:.0f} ms, reverse {t_rev*1e3:.1f} ms, ratio {t_fwd/t_rev:.0f}")
```
Counting: each JVP is two evaluations of $f$, each a matrix–vector product ($2d^2$ flops) plus $d$ `tanh`s, so (i) costs about $4d^3 = 4 \times 10^9$ flops. The reverse pass is one forward evaluation ($2d^2$), a VJP through `sum` and `tanh` ($O(d)$), and a VJP through the matmul (the rule (3.24): $W^\top\mathbf{g}$, $2d^2$): about $4d^2 = 4 \times 10^6$ flops. The flop ratio is $d = 1000$, the factor of §3.6. The wall clock shows that factor less cleanly. One run on the author's machine: (i) 0.25 s, (ii) 0.01 s through `Tensor`, (iii) 0.4 ms for the closed form, so (i)/(iii) $\approx 600$ and (i)/(ii) $\approx 25$. Two things are going on. The $4 \times 10^9$ flops of (i) are 2000 BLAS matrix–vector products, which run close to memory bandwidth, so they take a quarter of a second, not the tens of seconds a Python loop would. And the `Tensor` pass does whole-matrix work that the flop count ignores: it copies $W$ into a `Tensor`, allocates a $d \times d$ zero gradient for it, and in the matmul backward computes the unrequested $\nabla_W L = \mathbf{g}\mathbf{x}^\top$ (a $d \times d$ outer product) next to the $W^\top\mathbf{g}$ we want. Each of those is $O(d^2)$, the same order as evaluating $f$ once, so the reverse pass is still a small constant times the cost of $f$, which is the cheap gradient principle; the constant is larger than in the hand-written (iii), and because it is dominated by memory allocation, (ii) varies a lot between machines (we have seen $0.01$ to $0.08$ s). Scale $d$ to $2000$ and the same run gave (i) 1.7 s, (ii) 0.35 s, (iii) 2.4 ms: (i) grew by about $7\times$ ($d^3$ predicts $8\times$), (iii) by about $6\times$ ($d^2$ predicts $4\times$), and (ii) by far more, because its $d \times d$ allocations are now 32 MB each. The ratio that tracks §3.6 is (i)/(iii), about 700, near $d$. (The closed form (iii) is the reverse pass written by hand: $\nabla f = W^\top\sigma'(W\mathbf{x})$.)

Dual numbers: represent each intermediate as a pair $(u, \dot u)$ where $\dot u$ is its directional derivative along $\mathbf{v}$, and propagate $(u, \dot u) + (w, \dot w) = (u + w, \dot u + \dot w)$, $(u, \dot u)(w, \dot w) = (uw, \dot uw + u\dot w)$, $\tanh(u, \dot u) = (\tanh u, (1 - \tanh^2u)\dot u)$, and so on; the seed is $(\mathbf{x}, \mathbf{v})$. Every operation does a constant amount of extra work, so one pass costs a small multiple of $f$, and the full gradient still needs $d$ passes: the same $O(d \cdot \text{cost}(f))$ count as (i). The difference is that the tangent rules are exact derivative formulas, so there is no step $h$, no truncation error and no cancellation: the result is correct to rounding, unlike the $10^{-11}$ of central differences.

**C3.** (a) Let $\phi(A) = \log\det A$ and $g(A) := \operatorname{tr}(A^{-1}H)$ with the symmetric increment $H$ held fixed, so $d\phi = g(A)$ along $H$. Differentiating $g$ along the same increment $dA = H$, by (3.6): $dg = \operatorname{tr}\big(d(A^{-1})H\big) = -\operatorname{tr}(A^{-1}HA^{-1}H)$. So the second differential is $d^2\phi[H, H] = -\operatorname{tr}(A^{-1}HA^{-1}H)$. ∎

(b) For $A \succ 0$ the spectral theorem (LA-8) gives $A = Q\Lambda Q^\top$ with $\Lambda \succ 0$, and $A^{-1/2} := Q\Lambda^{-1/2}Q^\top$ is symmetric with $A^{-1/2}A^{-1/2} = A^{-1}$. Then
$$\operatorname{tr}(A^{-1}HA^{-1}H) = \operatorname{tr}(A^{-1/2}A^{-1/2}HA^{-1/2}A^{-1/2}H) = \operatorname{tr}\big(A^{-1/2}HA^{-1/2}\cdot A^{-1/2}HA^{-1/2}\big) = \operatorname{tr}(M^2),$$
moving the leading $A^{-1/2}$ to the end by the cyclic rule, with $M = A^{-1/2}HA^{-1/2}$. $M$ is symmetric ($M^\top = A^{-1/2}H^\top A^{-1/2} = M$), so $\operatorname{tr}(M^2) = \operatorname{tr}(M^\top M) = \|M\|_F^2 \ge 0$, with equality iff $M = 0$ iff $H = 0$ ($A^{-1/2}$ is invertible). Hence $d^2\phi[H, H] = -\|A^{-1/2}HA^{-1/2}\|_F^2 \le 0$ for every symmetric direction, which is the second-order condition for concavity (Theorem 3.14's reasoning, applied along every line in the cone). The positive-definite cone is convex (a convex combination of PD matrices is PD), so "concave on the cone" makes sense, and the inequality is strict: $\log\det$ is strictly concave there. ∎

(c) With $\boldsymbol{\mu} = \bar{\mathbf{x}}$ and $P = \Sigma^{-1}$: $\log\det\Sigma = -\log\det P$ and $\sum_i\mathbf{u}_i^\top\Sigma^{-1}\mathbf{u}_i = \sum_i\operatorname{tr}(P\mathbf{u}_i\mathbf{u}_i^\top) = \operatorname{tr}(PS)$, so
$$\ell(P) = \tfrac n2\log\det P - \tfrac12\operatorname{tr}(PS) + \text{const}.$$
The first term is concave in $P$ by (b) (positive multiple of a concave function); the second is linear in $P$; a concave plus a linear function is concave. So $\ell$ is strictly concave on the PD cone, and any stationary point is its unique global maximum. Stationary point: by (3.17) and (3.15), $\nabla_P\ell = \tfrac n2P^{-\top} - \tfrac12S^\top = \tfrac n2\Sigma - \tfrac12S$ (both symmetric), which vanishes iff $\Sigma = S/n$. Hence $\hat P = (S/n)^{-1}$ is the global maximizer in $P$. The map $\Sigma \mapsto \Sigma^{-1}$ is a bijection of the PD cone onto itself, and $\ell(\Sigma) = \ell(P(\Sigma))$ takes the same values at corresponding points, so the maximum of $\ell$ over $\Sigma$ is attained exactly at $\hat\Sigma = \hat P^{-1} = S/n$. (Concavity itself does not transfer: in $p = 1$, $\ell(\sigma^2) = -\tfrac n2\log\sigma^2 - \tfrac{S}{2\sigma^2}$ has second derivative $\tfrac{n}{2\sigma^4} - \tfrac{S}{\sigma^6}$, which is positive for $\sigma^2 > 2S/n$: not concave in $\sigma^2$, though concave in $1/\sigma^2$.)

(d) (3.39) is the gradient over all $p^2$ entries. If $\Sigma$ is parametrized by the $p(p + 1)/2$ entries $\theta_{ij}$, $i \le j$, with $\Sigma_{ij} = \Sigma_{ji} = \theta_{ij}$, then the chain rule gives $\partial\ell/\partial\theta_{ij} = \partial\ell/\partial\Sigma_{ij} + \partial\ell/\partial\Sigma_{ji} = 2G_{ij}$ for $i < j$ and $= G_{ii}$ on the diagonal, where $G$ is (3.39): the off-diagonal derivatives double. The location of the maximum does not change because $G$ is symmetric and $G = 0$ iff $2G_{ij} = 0$ for all $i < j$ and $G_{ii} = 0$: the two gradients vanish at the same points. A gradient-descent *path* would differ between the two parametrizations (different directions), but the fixed point is the same.
