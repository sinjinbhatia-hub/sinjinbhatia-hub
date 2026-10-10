# LA-1 Solutions

## A. Warm-up

**A1.** Column form: $A\mathbf{x} = 2\begin{pmatrix}1\\-2\\4\end{pmatrix} + (-1)\begin{pmatrix}3\\0\\1\end{pmatrix} = \begin{pmatrix}2\\-4\\8\end{pmatrix} + \begin{pmatrix}-3\\0\\-1\end{pmatrix} = \begin{pmatrix}-1\\-4\\7\end{pmatrix}$.

Row form: row 1: $1\cdot 2 + 3\cdot(-1) = -1$; row 2: $-2\cdot 2 + 0\cdot(-1) = -4$; row 3: $4\cdot 2 + 1\cdot(-1) = 7$. Same vector $(-1, -4, 7)^\top$. ✓

**A2.** $A = \begin{pmatrix}3 & 1\\1 & -2\end{pmatrix}$, $\mathbf{b} = \begin{pmatrix}7\\0\end{pmatrix}$. Column form: $x\begin{pmatrix}3\\1\end{pmatrix} + y\begin{pmatrix}1\\-2\end{pmatrix} = \begin{pmatrix}7\\0\end{pmatrix}$.

Solve: from the second equation $x = 2y$. Substitute: $6y + y = 7$, so $y = 1$, $x = 2$. Verify: $2\begin{pmatrix}3\\1\end{pmatrix} + 1\begin{pmatrix}1\\-2\end{pmatrix} = \begin{pmatrix}6\\2\end{pmatrix} + \begin{pmatrix}1\\-2\end{pmatrix} = \begin{pmatrix}7\\0\end{pmatrix}$. ✓

**A3.** $\mathbf{u}\cdot\mathbf{v} = 2 - 2 + 0 = 0$. $\|\mathbf{u}\| = \sqrt{1+4+4} = 3$. $\|\mathbf{v}\| = \sqrt{4+1+0} = \sqrt{5}$. $\cos\theta = 0/(3\sqrt5) = 0$, so $\theta = 90°$. The vectors are orthogonal.

**A4.** Using (3.3), the columns of $AB$ are $A\mathbf{b}_1$ and $A\mathbf{b}_2$:
$A\begin{pmatrix}0\\1\end{pmatrix} = \begin{pmatrix}2\\1\end{pmatrix}$, $A\begin{pmatrix}1\\0\end{pmatrix} = \begin{pmatrix}1\\0\end{pmatrix}$. So $AB = \begin{pmatrix}2 & 1\\1 & 0\end{pmatrix}$.

Columns of $BA$ are $B\mathbf{a}_1$ and $B\mathbf{a}_2$: $B\begin{pmatrix}1\\0\end{pmatrix} = \begin{pmatrix}0\\1\end{pmatrix}$, $B\begin{pmatrix}2\\1\end{pmatrix} = \begin{pmatrix}1\\2\end{pmatrix}$. So $BA = \begin{pmatrix}0 & 1\\1 & 2\end{pmatrix}$.

Not equal. Geometrically: $A$ is a **shear** (it fixes $\mathbf{e}_1$ and sends $\mathbf{e}_2 \mapsto (2, 1)^\top$; horizontal lines slide right in proportion to their height). $B$ is the **reflection** across the line $y = x$ (it swaps $\mathbf{e}_1$ and $\mathbf{e}_2$). $AB$ = reflect, then shear. $BA$ = shear, then reflect. Reflecting a sheared square is not the same as shearing a reflected square.

**A5.** $\|\mathbf{v}\|_1 = 7$, $\|\mathbf{v}\|_2 = 5$, $\|\mathbf{v}\|_\infty = 4$. Largest is $\ell^1$. Conjecture: $\|\mathbf{v}\|_\infty \le \|\mathbf{v}\|_2 \le \|\mathbf{v}\|_1$ always. (True; proved in C2.)

## B. Core

**B1.** Many answers. One set:

(a) $A = \begin{pmatrix}1 & 0\\0 & 1\end{pmatrix}$, $\mathbf{b} = \begin{pmatrix}1\\2\end{pmatrix}$. Row picture: lines $x = 1$ and $y = 2$ cross at $(1, 2)$. Column picture: columns $\mathbf{e}_1, \mathbf{e}_2$ span the plane; unique recipe $1\mathbf{e}_1 + 2\mathbf{e}_2$.

(b) $A = \begin{pmatrix}1 & 1\\1 & 1\end{pmatrix}$, $\mathbf{b} = \begin{pmatrix}1\\2\end{pmatrix}$. Row picture: $x + y = 1$ and $x + y = 2$, parallel lines, never meet. Column picture: both columns are $(1,1)^\top$; their span is the line $y = x$; $\mathbf{b} = (1, 2)^\top$ is not on it.

(c) Same $A$, $\mathbf{b} = \begin{pmatrix}1\\1\end{pmatrix}$. Row picture: the same line twice; every point on it is a solution. Column picture: $\mathbf{b}$ is on the span line; $x(1,1)^\top + y(1,1)^\top = (1,1)^\top$ holds for any $x + y = 1$.

**B2.** From (3.2), $(A(\mathbf{x} + \mathbf{y}))_i = \sum_j a_{ij}(x_j + y_j) = \sum_j a_{ij}x_j + \sum_j a_{ij}y_j = (A\mathbf{x})_i + (A\mathbf{y})_i$. Since this holds for every entry $i$, $A(\mathbf{x}+\mathbf{y}) = A\mathbf{x} + A\mathbf{y}$. Likewise $(A(c\mathbf{x}))_i = \sum_j a_{ij}(cx_j) = c\sum_j a_{ij}x_j = c(A\mathbf{x})_i$.

Why the same result: (3.1) and (3.2) are two ways of writing the *same* vector (proved in §3.2), so any property proved of one is a property of the other.

**B3.** $A\mathbf{0} = \sum_j 0\cdot\mathbf{a}_j = \mathbf{0}$. If $A\mathbf{x} = \mathbf{b}$ and $A\mathbf{y} = \mathbf{b}$, then by linearity $A(\mathbf{x} - \mathbf{y}) = A\mathbf{x} - A\mathbf{y} = \mathbf{b} - \mathbf{b} = \mathbf{0}$. If $\mathbf{x} \ne \mathbf{y}$ then $\mathbf{z} = \mathbf{x} - \mathbf{y} \ne \mathbf{0}$ and $A\mathbf{z} = \mathbf{0}$. So non-uniqueness of solutions to $A\mathbf{x} = \mathbf{b}$ is equivalent to the existence of a nonzero solution to the *homogeneous* system $A\mathbf{z} = \mathbf{0}$, and that depends only on $A$, not on $\mathbf{b}$.

**B4.** (a) $\mathbf{v}_1 + \mathbf{v}_2 = (1, 1, 2)^\top = \mathbf{v}_3$. So $\mathbf{v}_3 = 1\cdot\mathbf{v}_1 + 1\cdot\mathbf{v}_2$.

(b) Since $\mathbf{v}_3$ is a combination of the other two, $\operatorname{span}\{\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3\} = \operatorname{span}\{\mathbf{v}_1, \mathbf{v}_2\}$. Those two are not parallel, so the span is a plane through the origin. A normal vector $\mathbf{n}$ must satisfy $\mathbf{n}\cdot\mathbf{v}_1 = 0$ and $\mathbf{n}\cdot\mathbf{v}_2 = 0$: $n_1 + n_3 = 0$ and $n_2 + n_3 = 0$. Take $n_3 = -1$: $\mathbf{n} = (1, 1, -1)^\top$. The plane is $x + y - z = 0$. (Check: every $\mathbf{v}_i$ satisfies it.)

(c) $(1,1,1)^\top$: $1 + 1 - 1 = 1 \ne 0$. Not in the span. $(2, 3, 5)^\top$: $2 + 3 - 5 = 0$. In the span; explicitly $2\mathbf{v}_1 + 3\mathbf{v}_2 = (2, 3, 5)^\top$.

**B5.** Let $\mathbf{p} = \sum_i a_i\mathbf{v}_i$ and $\mathbf{q} = \sum_i b_i\mathbf{v}_i$. Then $\mathbf{p} + \mathbf{q} = \sum_i (a_i + b_i)\mathbf{v}_i$, a linear combination with coefficients $a_i + b_i$, hence in the span. And $c\mathbf{p} = \sum_i (ca_i)\mathbf{v}_i$, in the span with coefficients $ca_i$. Both steps use the distributive rules from Definition 3.1. ∎

**B6.** Expand both:
$\|\mathbf{u}+\mathbf{v}\|^2 = \|\mathbf{u}\|^2 + 2\,\mathbf{u}\cdot\mathbf{v} + \|\mathbf{v}\|^2$,
$\|\mathbf{u}-\mathbf{v}\|^2 = \|\mathbf{u}\|^2 - 2\,\mathbf{u}\cdot\mathbf{v} + \|\mathbf{v}\|^2$.
Add: the cross terms cancel, leaving $2\|\mathbf{u}\|^2 + 2\|\mathbf{v}\|^2$. ∎

Geometry: $\mathbf{u}$ and $\mathbf{v}$ are two adjacent sides of a parallelogram; $\mathbf{u}+\mathbf{v}$ and $\mathbf{u}-\mathbf{v}$ are its diagonals. The identity says: sum of squared diagonals = sum of squared sides (all four).

Counterexample for $\ell^1$: $\mathbf{u} = (1, 0)^\top$, $\mathbf{v} = (0, 1)^\top$. Then $\|\mathbf{u}+\mathbf{v}\|_1 = 2$, $\|\mathbf{u}-\mathbf{v}\|_1 = 2$, left side $= 4 + 4 = 8$. Right side $= 2\cdot 1 + 2\cdot 1 = 4$. $8 \ne 4$. So the $\ell^1$ norm does not come from any inner product.

**B7.** (a) $f(t) = (\mathbf{u} - t\mathbf{v})\cdot(\mathbf{u} - t\mathbf{v}) = \|\mathbf{u}\|^2 - 2t\,\mathbf{u}\cdot\mathbf{v} + t^2\|\mathbf{v}\|^2$. $f'(t) = -2\,\mathbf{u}\cdot\mathbf{v} + 2t\|\mathbf{v}\|^2 = 0 \Rightarrow t^* = \dfrac{\mathbf{u}\cdot\mathbf{v}}{\|\mathbf{v}\|^2}$. $f'' = 2\|\mathbf{v}\|^2 > 0$, so it's a minimum.

(b) $(\mathbf{u} - t^*\mathbf{v})\cdot\mathbf{v} = \mathbf{u}\cdot\mathbf{v} - t^*\|\mathbf{v}\|^2 = \mathbf{u}\cdot\mathbf{v} - \mathbf{u}\cdot\mathbf{v} = 0$. ∎

(c) Write $\mathbf{u} = t^*\mathbf{v} + (\mathbf{u} - t^*\mathbf{v})$. The two pieces are orthogonal by (b), so $\|\mathbf{u}\|^2 = \|t^*\mathbf{v}\|^2 + 2\,(t^*\mathbf{v})\cdot(\mathbf{u} - t^*\mathbf{v}) + \|\mathbf{u} - t^*\mathbf{v}\|^2$ and the middle term is $0$. ∎

(d) $\mathbf{u} = (2, 3, 6)^\top$, $\mathbf{v} = (1,1,1)^\top$. $\mathbf{u}\cdot\mathbf{v} = 11$, $\|\mathbf{v}\|^2 = 3$, $t^* = 11/3$. Projection $\mathbf{p} = (11/3, 11/3, 11/3)^\top$. Residual $\mathbf{r} = \mathbf{u} - \mathbf{p} = (-5/3, -2/3, 7/3)^\top$. Check (b): $\mathbf{r}\cdot\mathbf{v} = (-5 - 2 + 7)/3 = 0$. ✓ Check (c): $\|\mathbf{u}\|^2 = 4 + 9 + 36 = 49$. $\|\mathbf{p}\|^2 = 3\cdot(121/9) = 121/3$. $\|\mathbf{r}\|^2 = (25 + 4 + 49)/9 = 78/9 = 26/3$. Sum: $121/3 + 26/3 = 147/3 = 49$. ✓

**B8.** Reference implementation (also in `project.py` under "Your turn", so you can compare):

```python
def matmat_columns(A, B):
    """AB whose k-th column is A @ B[:, k], eq. (3.3)."""
    A = np.asarray(A, float); B = np.asarray(B, float)
    m, n = A.shape; n2, p = B.shape
    assert n == n2, "inner dimensions must match"
    C = np.zeros((m, p))
    for k in range(p):
        C[:, k] = matvec_columns(A, B[:, k])
    return C

A = rng.standard_normal((4, 3)); B = rng.standard_normal((3, 5))
assert np.allclose(matmat_columns(A, B), A @ B)
```

Why $(AB)\mathbf{x} = A(B\mathbf{x})$ is automatic: by (3.1), $(AB)\mathbf{x} = \sum_k x_k\,(\text{col}_k \text{ of } AB) = \sum_k x_k\,(A\mathbf{b}_k)$. By linearity of $A$, that equals $A\left(\sum_k x_k\mathbf{b}_k\right) = A(B\mathbf{x})$. The definition (3.3) was *chosen* to make this one-line argument work; "row times column" is the entry-level consequence, not the reason.

## C. Challenge

**C1.** Symmetric: $\sum_i w_i u_i v_i = \sum_i w_i v_i u_i$. Linear in the first slot: $\langle c\mathbf{u} + d\mathbf{u}', \mathbf{v}\rangle_w = \sum_i w_i(cu_i + du_i')v_i = c\langle\mathbf{u},\mathbf{v}\rangle_w + d\langle\mathbf{u}',\mathbf{v}\rangle_w$; by symmetry also linear in the second. Positive definite: $\langle\mathbf{v},\mathbf{v}\rangle_w = \sum_i w_i v_i^2 \ge 0$ since every $w_i > 0$, with equality only if every $v_i = 0$.

Cauchy–Schwarz: write $\langle\cdot,\cdot\rangle$ for $\langle\cdot,\cdot\rangle_w$ and $\|\mathbf{v}\|_w^2 = \langle\mathbf{v},\mathbf{v}\rangle$. If $\mathbf{v} = \mathbf{0}$, done. Else for all real $t$, by positive definiteness and bilinearity:
$$0 \le \langle\mathbf{u} - t\mathbf{v}, \mathbf{u} - t\mathbf{v}\rangle = \|\mathbf{u}\|_w^2 - 2t\langle\mathbf{u},\mathbf{v}\rangle + t^2\|\mathbf{v}\|_w^2.$$
A quadratic in $t$ with positive leading coefficient that is never negative has discriminant $\le 0$: $4\langle\mathbf{u},\mathbf{v}\rangle^2 - 4\|\mathbf{u}\|_w^2\|\mathbf{v}\|_w^2 \le 0$. ∎ Nothing about $\mathbb{R}^n$ or geometry was used: only the three properties. This is exactly the LA-8 proof.

Where: weighted least squares (ST-3) minimizes $\sum_i w_i (y_i - \hat y_i)^2$, which is $\|\mathbf{y} - \hat{\mathbf{y}}\|_w^2$. The Fisher information (ST-1) defines a weighted inner product on parameter space. The Mahalanobis distance (F5, ML-6) is the case where the weights become a full positive-definite matrix.

**C2.** Let $M = \|\mathbf{v}\|_\infty = \max_i|v_i|$.

(i) $\|\mathbf{v}\|_\infty \le \|\mathbf{v}\|_2$: $M^2 = v_k^2$ for the maximizing $k$, and $v_k^2 \le \sum_i v_i^2$. Equality iff all other $v_i = 0$, e.g. $\mathbf{e}_1$.

(ii) $\|\mathbf{v}\|_2 \le \|\mathbf{v}\|_1$: $\|\mathbf{v}\|_1^2 = \left(\sum_i|v_i|\right)^2 = \sum_i v_i^2 + \sum_{i\ne j}|v_i||v_j| \ge \sum_i v_i^2$. Equality iff all cross terms vanish iff at most one $v_i \ne 0$, e.g. $\mathbf{e}_1$.

(iii) $\|\mathbf{v}\|_1 \le \sqrt n\|\mathbf{v}\|_2$: Cauchy–Schwarz with $\mathbf{a} = (|v_1|, \dots, |v_n|)^\top$ and $\mathbf{1} = (1, \dots, 1)^\top$: $\|\mathbf{v}\|_1 = \mathbf{a}\cdot\mathbf{1} \le \|\mathbf{a}\|\|\mathbf{1}\| = \|\mathbf{v}\|_2\sqrt n$. Equality iff $\mathbf{a} \parallel \mathbf{1}$ iff all $|v_i|$ equal, e.g. $\mathbf{1}$.

(iv) $\sqrt n\|\mathbf{v}\|_2 \le n\|\mathbf{v}\|_\infty$: $\|\mathbf{v}\|_2^2 = \sum_i v_i^2 \le nM^2$, so $\|\mathbf{v}\|_2 \le \sqrt n M$, multiply by $\sqrt n$. Equality iff all $|v_i| = M$, e.g. $\mathbf{1}$.

**C3.** (a) $\sum_{i,j}a_{ij}^2 = \sum_j\left(\sum_i a_{ij}^2\right) = \sum_j\|\mathbf{a}_j\|^2$. And $(A^\top A)_{jj} = \sum_i (A^\top)_{ji}a_{ij} = \sum_i a_{ij}^2 = \|\mathbf{a}_j\|^2$, so $\operatorname{tr}(A^\top A) = \sum_j\|\mathbf{a}_j\|^2$.

(b) $(A\mathbf{x})_i = \mathbf{r}_i\cdot\mathbf{x}$, so by Cauchy–Schwarz $|(A\mathbf{x})_i| \le \|\mathbf{r}_i\|\|\mathbf{x}\|$. Then $\|A\mathbf{x}\|^2 = \sum_i (A\mathbf{x})_i^2 \le \sum_i\|\mathbf{r}_i\|^2\|\mathbf{x}\|^2 = \|A\|_F^2\|\mathbf{x}\|^2$ (sum of squared row lengths is also $\|A\|_F^2$). ∎

(c) Strict: $A = I_2$, $\mathbf{x} = \mathbf{e}_1$: $\|A\mathbf{x}\| = 1 < \sqrt2 = \|I\|_F$. Equality: $A = \begin{pmatrix}1 & 1\\0 & 0\end{pmatrix}$, $\mathbf{x} = (1, 1)^\top/\sqrt2$: $A\mathbf{x} = (\sqrt2, 0)^\top$, $\|A\mathbf{x}\| = \sqrt2 = \|A\|_F\cdot 1$. In general $\|A\|_F$ is **not** the smallest constant: for $A = I_n$ the smallest $K$ is $1$ while $\|I_n\|_F = \sqrt n$. The smallest $K$ is the operator norm $\|A\|_2 = \max_{\|\mathbf{x}\|=1}\|A\mathbf{x}\|$, and LA-9 shows $\|A\|_2 = \sigma_{\max}(A) \le \|A\|_F = \sqrt{\sum_k\sigma_k^2}$, with equality iff $A$ has rank $\le 1$.
