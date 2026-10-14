# LA-5 Solutions

## A. Warm-up

**A1.** (a) $\mathbf{u}^\top\mathbf{v} = 2 - 4 + 2 = 0$, so $\mathbf{u} \perp \mathbf{v}$. $\mathbf{u} + \mathbf{v} = (3, 0, 3)^\top$, $\|\mathbf{u} + \mathbf{v}\|^2 = 18$; $\|\mathbf{u}\|^2 = 1 + 4 + 4 = 9$, $\|\mathbf{v}\|^2 = 4 + 4 + 1 = 9$; $9 + 9 = 18$ ✓.

(b) We need $w_1 - 2w_2 + 2w_3 = 0$ and $2w_1 + 2w_2 + w_3 = 0$. Subtract twice the first from the second: $6w_2 - 3w_3 = 0$, so $w_3 = 2w_2$; then $w_1 = 2w_2 - 2w_3 = -2w_2$. Take $w_2 = 1$: $\mathbf{w} = (-2, 1, 2)^\top$. (The cross product $\mathbf{u}\times\mathbf{v} = (-6, 3, 6)^\top$ is a multiple.) Check: $\mathbf{u}^\top\mathbf{w} = -2 - 2 + 4 = 0$, $\mathbf{v}^\top\mathbf{w} = -4 + 2 + 2 = 0$ ✓.

(c) Theorem 3.3: three pairwise orthogonal nonzero vectors are linearly independent, and three independent vectors in $\mathbb{R}^3$ are a basis. No elimination needed.

**A2.** $\mathbf{x} \in S^\perp$ iff $(1, 2, 2)\mathbf{x} = 0$, i.e. $S^\perp = N(B)$ with $B = \begin{pmatrix}1 & 2 & 2\end{pmatrix}$, the plane $x_1 + 2x_2 + 2x_3 = 0$. Free variables $x_2, x_3$: $x_1 = -2x_2 - 2x_3$, basis $(-2, 1, 0)^\top, (-2, 0, 1)^\top$. $\dim S^\perp = 2$, and $1 + 2 = 3$ ✓ (Theorem 3.7(a) with $k = 1$).

$(S^\perp)^\perp$: vectors $\mathbf{n}$ with $-2n_1 + n_2 = 0$ and $-2n_1 + n_3 = 0$, i.e. $n_2 = n_3 = 2n_1$: $\mathbf{n} = n_1(1, 2, 2)^\top$. That is $S$ ✓ (Theorem 3.7(c)).

**A3.** $\mathbf{a}^\top\mathbf{b} = 3 + 0 + 6 = 9$, $\mathbf{a}^\top\mathbf{a} = 1 + 4 + 4 = 9$, so $\hat x = 1$ and $\mathbf{p} = \mathbf{a} = (1, 2, 2)^\top$. $\mathbf{e} = \mathbf{b} - \mathbf{p} = (2, -2, 1)^\top$; $\mathbf{a}^\top\mathbf{e} = 2 - 4 + 2 = 0$ ✓.

$P = \mathbf{a}\mathbf{a}^\top/9 = \frac19\begin{pmatrix}1 & 2 & 2\\2 & 4 & 4\\2 & 4 & 4\end{pmatrix}$. $P^2 = \frac{1}{81}\mathbf{a}(\mathbf{a}^\top\mathbf{a})\mathbf{a}^\top = \frac{9}{81}\mathbf{a}\mathbf{a}^\top = P$ ✓. (Also $\operatorname{tr}P = (1 + 4 + 4)/9 = 1$, the dimension of the line; B4(d).)

Pythagoras: $\|\mathbf{b}\|^2 = 9 + 0 + 9 = 18$, $\|\mathbf{p}\|^2 = 9$, $\|\mathbf{e}\|^2 = 4 + 4 + 1 = 9$; $9 + 9 = 18$ ✓.

**A4.** $A = \begin{pmatrix}1 & 0\\1 & 1\\1 & 2\end{pmatrix}$, $\mathbf{b} = (0, 1, 3)^\top$. $A^\top A = \begin{pmatrix}3 & 3\\3 & 5\end{pmatrix}$, $A^\top\mathbf{b} = \begin{pmatrix}0 + 1 + 3\\0 + 1 + 6\end{pmatrix} = \begin{pmatrix}4\\7\end{pmatrix}$.

$3C + 3D = 4$, $3C + 5D = 7$. Subtract: $2D = 3$, $D = \tfrac32$. Then $3C = 4 - \tfrac92 = -\tfrac12$, $C = -\tfrac16$. Line: $b = -\tfrac16 + \tfrac32t$.

Fitted values $\mathbf{p} = (-\tfrac16, \tfrac43, \tfrac{17}{6})^\top$. Residual $\mathbf{r} = \mathbf{b} - \mathbf{p} = (\tfrac16, -\tfrac13, \tfrac16)^\top$. Column 1: $\tfrac16 - \tfrac26 + \tfrac16 = 0$ ✓. Column 2: $0 - \tfrac13 + \tfrac26 = 0$ ✓. $\|\mathbf{r}\|^2 = \tfrac{1 + 4 + 1}{36} = \tfrac16$.

**A5.** (a) Columns $\mathbf{q}_1 = \tfrac13(1, -2, 2)^\top$, $\mathbf{q}_2 = \tfrac13(2, 2, 1)^\top$, $\mathbf{q}_3 = \tfrac13(-2, 1, 2)^\top$. Dot products (times 9): $\mathbf{q}_1\cdot\mathbf{q}_1$: $1 + 4 + 4 = 9$; $\mathbf{q}_2\cdot\mathbf{q}_2$: $4 + 4 + 1 = 9$; $\mathbf{q}_3\cdot\mathbf{q}_3$: $4 + 1 + 4 = 9$; $\mathbf{q}_1\cdot\mathbf{q}_2$: $2 - 4 + 2 = 0$; $\mathbf{q}_1\cdot\mathbf{q}_3$: $-2 - 2 + 4 = 0$; $\mathbf{q}_2\cdot\mathbf{q}_3$: $-4 + 2 + 2 = 0$. So $Q^\top Q = I$ and $Q^{-1} = Q^\top = \frac13\begin{pmatrix}1 & -2 & 2\\2 & 2 & 1\\-2 & 1 & 2\end{pmatrix}$.

(b) $Q\mathbf{x} = \tfrac13\begin{pmatrix}1 + 4 - 6\\-2 + 4 + 3\\2 + 2 + 6\end{pmatrix} = \tfrac13(-1, 5, 10)^\top$. $\|Q\mathbf{x}\|^2 = \tfrac{1 + 25 + 100}{9} = 14 = 1 + 4 + 9 = \|\mathbf{x}\|^2$ ✓ (Theorem 3.18).

(c) Coefficients $\mathbf{q}_i^\top\mathbf{x}$: $\mathbf{q}_1^\top\mathbf{x} = \tfrac13(1 - 4 + 6) = 1$, $\mathbf{q}_2^\top\mathbf{x} = \tfrac13(2 + 4 + 3) = 3$, $\mathbf{q}_3^\top\mathbf{x} = \tfrac13(-2 + 2 + 6) = 2$. So $\mathbf{x} = \mathbf{q}_1 + 3\mathbf{q}_2 + 2\mathbf{q}_3$; check the first entry: $\tfrac13(1 + 6 - 4) = 1$ ✓. Parseval: $1^2 + 3^2 + 2^2 = 14 = \|\mathbf{x}\|^2$ ✓.

**A6.** $\|\mathbf{a}_1\| = 3$, $\mathbf{q}_1 = \tfrac13(1, 2, 2)^\top$, $r_{11} = 3$. $r_{12} = \mathbf{q}_1^\top\mathbf{a}_2 = \tfrac13(0 + 6 + 0) = 2$. $\mathbf{v}_2 = \mathbf{a}_2 - 2\mathbf{q}_1 = (0, 3, 0)^\top - (\tfrac23, \tfrac43, \tfrac43)^\top = (-\tfrac23, \tfrac53, -\tfrac43)^\top$. $r_{22} = \|\mathbf{v}_2\| = \tfrac13\sqrt{4 + 25 + 16} = \tfrac{\sqrt{45}}{3} = \sqrt5$. $\mathbf{q}_2 = \mathbf{v}_2/\sqrt5 = \tfrac{1}{3\sqrt5}(-2, 5, -4)^\top$. Check $\mathbf{q}_1^\top\mathbf{q}_2 = \tfrac{1}{9\sqrt5}(-2 + 10 - 8) = 0$ ✓.

$R = \begin{pmatrix}3 & 2\\0 & \sqrt5\end{pmatrix}$. Verify: $r_{12}\mathbf{q}_1 + r_{22}\mathbf{q}_2 = \tfrac23(1, 2, 2)^\top + \tfrac13(-2, 5, -4)^\top = \tfrac13(2 - 2,\ 4 + 5,\ 4 - 4)^\top = (0, 3, 0)^\top = \mathbf{a}_2$ ✓.

## B. Core

**B1.** (a) If $\mathbf{v} \in S$ and $\mathbf{v} \in S^\perp$, then $\mathbf{v}$ is orthogonal to itself: $\mathbf{v}^\top\mathbf{v} = \|\mathbf{v}\|^2 = 0$, so $\mathbf{v} = \mathbf{0}$. ∎

(b) Let $\mathbf{v} \in T^\perp$, so $\mathbf{v}^\top\mathbf{t} = 0$ for all $\mathbf{t} \in T$. Every $\mathbf{s} \in S$ is also in $T$, so $\mathbf{v}^\top\mathbf{s} = 0$ for all $\mathbf{s} \in S$, i.e. $\mathbf{v} \in S^\perp$. ∎ (Bigger subspace, smaller complement.)

(c) ($\subseteq$) If $\mathbf{v} \perp S + T$, then in particular $\mathbf{v} \perp \mathbf{s} + \mathbf{0}$ for all $\mathbf{s} \in S$ and $\mathbf{v} \perp \mathbf{0} + \mathbf{t}$ for all $\mathbf{t} \in T$, so $\mathbf{v} \in S^\perp \cap T^\perp$. ($\supseteq$) If $\mathbf{v} \in S^\perp \cap T^\perp$, then $\mathbf{v}^\top(\mathbf{s} + \mathbf{t}) = \mathbf{v}^\top\mathbf{s} + \mathbf{v}^\top\mathbf{t} = 0 + 0$ for all $\mathbf{s}, \mathbf{t}$, so $\mathbf{v} \perp S + T$. ∎

**B2.** (a) $A\mathbf{x} = \mathbf{0}$: $x_1 + x_2 + x_3 + x_4 = 0$ and $x_1 + 2x_2 + 3x_3 + 4x_4 = 0$. Subtracting, $x_2 + 2x_3 + 3x_4 = 0$, so $x_2 = -2x_3 - 3x_4$ and then $x_1 = -x_2 - x_3 - x_4 = x_3 + 2x_4$. Free $x_3, x_4$: basis $\mathbf{n}_1 = (1, -2, 1, 0)^\top$, $\mathbf{n}_2 = (2, -3, 0, 1)^\top$. Against row 1 $(1,1,1,1)$: $1 - 2 + 1 + 0 = 0$, $2 - 3 + 0 + 1 = 0$ ✓. Against row 2 $(1,2,3,4)$: $1 - 4 + 3 + 0 = 0$, $2 - 6 + 0 + 4 = 0$ ✓. Dimensions: $\operatorname{rank}A = 2$, $\dim N(A) = 2$, $2 + 2 = 4$ ✓.

(b) $C(A^\top)$ has basis the two rows, i.e. the columns of $A^\top$ (a $4\times2$ matrix with independent columns). Normal equations (3.6) with $A^\top$ in place of $A$: $(A^\top)^\top A^\top\,\mathbf{y} = (A^\top)^\top\mathbf{x}$, i.e. $AA^\top\mathbf{y} = A\mathbf{x}$. $AA^\top = \begin{pmatrix}4 & 10\\10 & 30\end{pmatrix}$, $A\mathbf{x} = (1, 1)^\top$. Determinant $120 - 100 = 20$, inverse $\frac{1}{20}\begin{pmatrix}30 & -10\\-10 & 4\end{pmatrix}$, so $\mathbf{y} = \frac{1}{20}(20, -6)^\top = (1, -0.3)^\top$.

$\mathbf{x}_{\text{row}} = A^\top\mathbf{y} = 1\cdot(1, 1, 1, 1)^\top - 0.3\,(1, 2, 3, 4)^\top = (0.7, 0.4, 0.1, -0.2)^\top$.
$\mathbf{x}_{\text{null}} = \mathbf{x} - \mathbf{x}_{\text{row}} = (0.3, -0.4, -0.1, 0.2)^\top$.

Check $\mathbf{x}_{\text{null}} \in N(A)$: $0.3 - 0.4 - 0.1 + 0.2 = 0$ and $0.3 - 0.8 - 0.3 + 0.8 = 0$ ✓; in the basis, $\mathbf{x}_{\text{null}} = -0.1\,\mathbf{n}_1 + 0.2\,\mathbf{n}_2$ (fourth entry forces $0.2$, third forces $-0.1$; first entry $-0.1 + 0.4 = 0.3$ ✓). Check $A\mathbf{x}_{\text{row}} = (0.7 + 0.4 + 0.1 - 0.2,\ 0.7 + 0.8 + 0.3 - 0.8)^\top = (1, 1)^\top = A\mathbf{x}$ ✓. The matrix sees only the row-space part.

**B3.** By Theorem 3.11, $N(A^\top A) = N(A)$. Both $A$ ($m\times n$) and $A^\top A$ ($n\times n$) have $n$ columns, so Rank–Nullity gives $\operatorname{rank}(A^\top A) = n - \dim N(A^\top A) = n - \dim N(A) = \operatorname{rank}(A)$.

$C(A^\top A) \subseteq C(A^\top)$: every vector $A^\top A\mathbf{x} = A^\top(A\mathbf{x})$ is $A^\top$ times something. Both subspaces have dimension $\operatorname{rank}(A)$ (row rank equals column rank, so $\dim C(A^\top) = \operatorname{rank}A$). A subspace contained in another of the same dimension equals it. ∎

If $A^\top A = 0$ then $N(A^\top A) = \mathbb{R}^n = N(A)$, so $A\mathbf{x} = \mathbf{0}$ for every $\mathbf{x}$, including each $\mathbf{e}_j$; thus every column of $A$ is zero and $A = 0$. (Alternatively: the diagonal of $A^\top A$ is $\|\mathbf{a}_j\|^2$.)

For $AA^\top$ apply everything to $A^\top$: $N(AA^\top) = N(A^\top)$, $\operatorname{rank}(AA^\top) = \operatorname{rank}(A)$, $C(AA^\top) = C(A)$.

**B4.** (a) $A = \begin{pmatrix}1 & 0\\-1 & 1\\0 & -1\end{pmatrix}$, $A^\top A = \begin{pmatrix}2 & -1\\-1 & 2\end{pmatrix}$, $(A^\top A)^{-1} = \tfrac13\begin{pmatrix}2 & 1\\1 & 2\end{pmatrix}$. $A(A^\top A)^{-1} = \tfrac13\begin{pmatrix}2 & 1\\-2 + 1 & -1 + 2\\-1 & -2\end{pmatrix} = \tfrac13\begin{pmatrix}2 & 1\\-1 & 1\\-1 & -2\end{pmatrix}$. Times $A^\top = \begin{pmatrix}1 & -1 & 0\\0 & 1 & -1\end{pmatrix}$: row 1: $(2, -2 + 1, -1) = (2, -1, -1)$; row 2: $(-1, 1 + 1, -1) = (-1, 2, -1)$; row 3: $(-1, 1 - 2, 2) = (-1, -1, 2)$. So
$$P = \frac13\begin{pmatrix}2 & -1 & -1\\-1 & 2 & -1\\-1 & -1 & 2\end{pmatrix}.$$

(b) $\mathbf{n}\mathbf{n}^\top/\mathbf{n}^\top\mathbf{n} = \tfrac13\begin{pmatrix}1 & 1 & 1\\1 & 1 & 1\\1 & 1 & 1\end{pmatrix}$, and $I - $ that $= \tfrac13\begin{pmatrix}2 & -1 & -1\\-1 & 2 & -1\\-1 & -1 & 2\end{pmatrix}$, the same ✓. They must agree because the orthogonal projection onto a subspace is unique (Theorem 3.7(b)): it does not matter which basis of $S$, or which description of $S$ ($=\mathbf{n}^\perp$), you compute it from.

(c) $P\mathbf{b} = \tfrac13(2 - 2 - 3,\ -1 + 4 - 3,\ -1 - 2 + 6)^\top = (-1, 0, 1)^\top = \mathbf{p}$. $\mathbf{e} = \mathbf{b} - \mathbf{p} = (2, 2, 2)^\top = 2\mathbf{n}$ ✓, and indeed $\mathbf{n}^\top\mathbf{b}/\mathbf{n}^\top\mathbf{n} = 6/3 = 2$. Check $\mathbf{p} \in S$: $-1 + 0 + 1 = 0$ ✓.

(d) Write $P = \tfrac13M$. $M^2$, first row: $(4 + 1 + 1,\ -2 - 2 + 1,\ -2 + 1 - 2) = (6, -3, -3) = 3\,(2, -1, -1)$; by the symmetry of $M$ (each row is a cyclic shift of the others) every row of $M^2$ is $3$ times the corresponding row of $M$, so $M^2 = 3M$ and $P^2 = \tfrac19\cdot3M = P$ ✓. $\operatorname{tr}P = \tfrac13(2 + 2 + 2) = 2 = \dim S$.

Trace lemma (this is LA-4 Theorem 3.10; we re-prove it): $\operatorname{tr}(XY) = \sum_i(XY)_{ii} = \sum_i\sum_j x_{ij}y_{ji} = \sum_j\sum_i y_{ji}x_{ij} = \sum_j(YX)_{jj} = \operatorname{tr}(YX)$. ∎ Now with $X = A$ ($m\times n$) and $Y = (A^\top A)^{-1}A^\top$ ($n\times m$): $\operatorname{tr}\big(A(A^\top A)^{-1}A^\top\big) = \operatorname{tr}\big((A^\top A)^{-1}A^\top A\big) = \operatorname{tr}(I_n) = n$. ∎ So the trace of a projection matrix is the dimension of the subspace it projects onto. LA-4 Exercise C3(b) already explains why: in a basis made of a basis of $C(P)$ followed by a basis of $N(P)$, a projection has matrix $\operatorname{diag}(1, \dots, 1, 0, \dots, 0)$ with $\operatorname{rank}P$ ones, and the trace does not depend on the basis. In LA-7's language (Theorem 3.4, trace equals the sum of the eigenvalues), the eigenvalues of $P$ are $1$, $n$ times, and $0$, $m - n$ times.

**B5.** (a) $A = \mathbf{1}$ (the $m\times1$ column of ones). $A^\top A = \mathbf{1}^\top\mathbf{1} = m$, $A^\top\mathbf{b} = \sum_ib_i$. Normal equation: $mC = \sum_ib_i$, $C = \bar b$. ∎ The mean is the orthogonal projection of the data onto the line of constant vectors.

(b) The normal equations (3.5) say $A^\top\mathbf{r} = \mathbf{0}$, one entry per column. The entry for the column of ones is $\mathbf{1}^\top\mathbf{r} = \sum_ir_i = 0$. ∎

(c) $\bar t = 1.5$, $\mathbf{t}' = (-1.5, -0.5, 0.5, 1.5)^\top$. $A' = (\mathbf{1}\ \ \mathbf{t}')$, $A'^\top A' = \begin{pmatrix}4 & \sum t_i'\\ \sum t_i' & \sum t_i'^2\end{pmatrix} = \begin{pmatrix}4 & 0\\0 & 5\end{pmatrix}$, diagonal because $\sum t_i' = \sum t_i - m\bar t = 0$: centering makes the two columns orthogonal. $A'^\top\mathbf{b} = (9,\ -1.5 - 1 + 1 + 6)^\top = (9, 4.5)^\top$. So $C' = 9/4 = 2.25$, $D' = 4.5/5 = 0.9$. The fitted line is $b = 2.25 + 0.9(t - 1.5) = 0.9 + 0.9t$: the same line as Example 4.2.

Why it must be: $\operatorname{span}\{\mathbf{1}, \mathbf{t}'\} = \operatorname{span}\{\mathbf{1}, \mathbf{t}\}$ (each spanning set is a combination of the other), so both fits project $\mathbf{b}$ onto the *same* subspace, and the projection is unique (Theorem 3.7(b) / Theorem 3.14). Only the coordinates $\hat{\mathbf{x}}$ change, because the basis changed.

(d) From (c), with a diagonal Gram matrix, $D = \sum_it_i'b_i\big/\sum_it_i'^2$. Since $\sum_it_i' = 0$, $\sum_it_i'\bar b = 0$, so the numerator equals $\sum_it_i'(b_i - \bar b) = \sum_i(t_i - \bar t)(b_i - \bar b)$. ∎ That is the sample covariance of $t$ and $b$ over the sample variance of $t$. And $C' = \bar b$ by (a) applied to the orthogonal column: the line passes through $(\bar t, \bar b)$.

**B6.** (a) $f = \mathbf{x}^\top G\mathbf{x} - 2\mathbf{c}^\top\mathbf{x} + \|\mathbf{b}\|^2 = g_{11}x_1^2 + 2g_{12}x_1x_2 + g_{22}x_2^2 - 2c_1x_1 - 2c_2x_2 + \|\mathbf{b}\|^2$, using $g_{21} = g_{12}$. Then
$$\frac{\partial f}{\partial x_1} = 2g_{11}x_1 + 2g_{12}x_2 - 2c_1 = 2(G\mathbf{x} - \mathbf{c})_1, \qquad \frac{\partial f}{\partial x_2} = 2g_{12}x_1 + 2g_{22}x_2 - 2c_2 = 2(G\mathbf{x} - \mathbf{c})_2. $$
So $\nabla f = 2(A^\top A\mathbf{x} - A^\top\mathbf{b}) = 2A^\top(A\mathbf{x} - \mathbf{b})$, which is (3.10). The cross term $2g_{12}x_1x_2$ contributes to *both* partials, which is the bookkeeping that the general proof's "two sums" handles.

(b) $f(\mathbf{x} + h\mathbf{d}) = \|A\mathbf{x} - \mathbf{b} + hA\mathbf{d}\|^2 = \|A\mathbf{x} - \mathbf{b}\|^2 + 2h\,(A\mathbf{d})^\top(A\mathbf{x} - \mathbf{b}) + h^2\|A\mathbf{d}\|^2$. The coefficient of $h$ is the directional derivative $\nabla f(\mathbf{x})^\top\mathbf{d}$, and $(A\mathbf{d})^\top(A\mathbf{x} - \mathbf{b}) = \mathbf{d}^\top A^\top(A\mathbf{x} - \mathbf{b})$. Since this holds for every $\mathbf{d}$, $\nabla f = 2A^\top(A\mathbf{x} - \mathbf{b})$. ∎ The $h^2$ coefficient, $\|A\mathbf{d}\|^2 = \tfrac12\mathbf{d}^\top(2A^\top A)\mathbf{d} \ge 0$, is the Hessian form.

(c) Yes, a minimizer exists. The projection $\mathbf{p}$ of $\mathbf{b}$ onto $C(A)$ exists regardless of independence (Theorem 3.7(b)), and $\mathbf{p} \in C(A)$ means $\mathbf{p} = A\hat{\mathbf{x}}$ for some $\hat{\mathbf{x}}$; by Theorem 3.14 this $\hat{\mathbf{x}}$ minimizes $f$. It is not unique: if $A\mathbf{z} = \mathbf{0}$ then $A(\hat{\mathbf{x}} + \mathbf{z}) = \mathbf{p}$ too. Conversely any minimizer $\mathbf{x}^*$ has $A\mathbf{x}^* = \mathbf{p}$ (the closest point is unique), so $A(\mathbf{x}^* - \hat{\mathbf{x}}) = \mathbf{0}$. The set of minimizers is the affine subspace $\hat{\mathbf{x}} + N(A)$, of dimension $n - \operatorname{rank}A$. The normal equations $A^\top A\mathbf{x} = A^\top\mathbf{b}$ are still consistent (by B3, $A^\top\mathbf{b} \in C(A^\top) = C(A^\top A)$), they have many solutions. LA-9's pseudoinverse picks the one of smallest norm, which is the one in $C(A^\top)$.

**B7.** (a) $(Q_1Q_2)^\top(Q_1Q_2) = Q_2^\top Q_1^\top Q_1Q_2 = Q_2^\top IQ_2 = I$. For the inverse, $Q^{-1} = Q^\top$, and $(Q^\top)^\top Q^\top = QQ^\top = I$ (square orthogonal matrices satisfy both $Q^\top Q = I$ and $QQ^\top = I$). ∎

(b) Let $Q = \begin{pmatrix}a & b\\c & d\end{pmatrix}$ with orthonormal columns. $a^2 + c^2 = 1$, so $(a, c) = (\cos\theta, \sin\theta)$ for some $\theta$. The second column is a unit vector orthogonal to $(\cos\theta, \sin\theta)$; the unit vectors orthogonal to it are exactly $\pm(-\sin\theta, \cos\theta)$. The plus sign gives the rotation matrix, the minus sign gives $\begin{pmatrix}\cos\theta & \sin\theta\\ \sin\theta & -\cos\theta\end{pmatrix}$. Apply the latter to $\mathbf{u} = (\cos\tfrac\theta2, \sin\tfrac\theta2)^\top$: first entry $\cos\theta\cos\tfrac\theta2 + \sin\theta\sin\tfrac\theta2 = \cos(\theta - \tfrac\theta2) = \cos\tfrac\theta2$; second entry $\sin\theta\cos\tfrac\theta2 - \cos\theta\sin\tfrac\theta2 = \sin(\theta - \tfrac\theta2) = \sin\tfrac\theta2$. So $Q\mathbf{u} = \mathbf{u}$: the line at angle $\theta/2$ is fixed pointwise. Being symmetric and orthogonal it is its own inverse, and it sends the perpendicular direction $(-\sin\tfrac\theta2, \cos\tfrac\theta2)^\top$ to its negative (same computation with $\theta/2 + \pi/2$). That is the reflection across the line at angle $\theta/2$. ∎

(c) Polarization: $\|\mathbf{x} + \mathbf{y}\|^2 - \|\mathbf{x} - \mathbf{y}\|^2 = (\|\mathbf{x}\|^2 + 2\mathbf{x}^\top\mathbf{y} + \|\mathbf{y}\|^2) - (\|\mathbf{x}\|^2 - 2\mathbf{x}^\top\mathbf{y} + \|\mathbf{y}\|^2) = 4\mathbf{x}^\top\mathbf{y}$. If $Q$ preserves all lengths then $(Q\mathbf{x})^\top(Q\mathbf{y}) = \tfrac14(\|Q(\mathbf{x} + \mathbf{y})\|^2 - \|Q(\mathbf{x} - \mathbf{y})\|^2) = \tfrac14(\|\mathbf{x} + \mathbf{y}\|^2 - \|\mathbf{x} - \mathbf{y}\|^2) = \mathbf{x}^\top\mathbf{y}$. So $\mathbf{x}^\top(Q^\top Q)\mathbf{y} = \mathbf{x}^\top\mathbf{y}$ for all $\mathbf{x}, \mathbf{y}$; with $\mathbf{x} = \mathbf{e}_i$, $\mathbf{y} = \mathbf{e}_j$ this says $(Q^\top Q)_{ij} = \delta_{ij}$, i.e. $Q^\top Q = I$. ∎

**B8.** (a) $\sum_iw_i\big((A\mathbf{x})_i - b_i\big)^2 = \sum_i\big(\sqrt{w_i}(A\mathbf{x} - \mathbf{b})_i\big)^2 = \|W^{1/2}(A\mathbf{x} - \mathbf{b})\|^2 = \|(W^{1/2}A)\mathbf{x} - W^{1/2}\mathbf{b}\|^2$ with $W^{1/2} = \operatorname{diag}(\sqrt{w_i})$, since a diagonal matrix scales entries. This is (3.9) for $\tilde A = W^{1/2}A$, $\tilde{\mathbf{b}} = W^{1/2}\mathbf{b}$. Normal equations: $\tilde A^\top\tilde A\hat{\mathbf{x}} = \tilde A^\top\tilde{\mathbf{b}}$, i.e. $A^\top W^{1/2}W^{1/2}A\hat{\mathbf{x}} = A^\top W^{1/2}W^{1/2}\mathbf{b}$:
$$A^\top WA\,\hat{\mathbf{x}} = A^\top W\mathbf{b}.$$
Geometrically it is still a projection, but with the inner product $\mathbf{u}^\top W\mathbf{v}$ (LA-1, Exercise C1) in place of the dot product.

(b) If row $i$ appears $w_i$ times in the data, its residual appears $w_i$ times in the sum of squares: $\sum_iw_ir_i^2$. The two problems have identical objective functions, so identical minimizers. (Equivalently, stacking $w_i$ copies of row $i$ gives a design matrix $\check A$ with $\check A^\top\check A = A^\top WA$ and $\check A^\top\check{\mathbf{b}} = A^\top W\mathbf{b}$.) A weight is a count of how many times you would be willing to have observed the point.

(c) Reference implementation (also in `project.py` under "Your turn"):

```python
def least_squares_weighted(A, b, w):
    """argmin Σ w_i (Ax − b)_i² = argmin ‖W^{1/2}(Ax − b)‖²: scale row i by √w_i, then ordinary least squares."""
    s = np.sqrt(np.asarray(w, float))
    return least_squares_qr(s[:, None] * np.asarray(A, float), s * np.asarray(b, float))[0]

assert np.allclose(least_squares_weighted(A2, y2, np.ones(4)), xh2)        # weights 1: Example 4.2 again
xw = least_squares_weighted(A2, y2, [1, 1, 1e6, 1])
assert abs((A2 @ xw)[2] - y2[2]) < 1e-5                                    # the line passes through (2, 2)
c5w = least_squares_weighted(A5, V5, [1, 1, 1, 1, 1, 10])                  # Example 4.5, last reading weighted 10
assert np.allclose(c5w, [2377 / 4700, 7367 / 9400, 233 / 9400])
```

With the huge weight, the fit is $C + 2D = 2$ to five decimals: the third residual is forced to (nearly) zero and the other three points share what freedom is left (one parameter). For Example 4.5 with $w = (1, 1, 1, 1, 1, 10)$, solving $A^\top WA\,\hat{\mathbf{x}} = A^\top W\mathbf{V}$ exactly gives $(c_0, c_1, c_2) = \big(\tfrac{2377}{4700}, \tfrac{7367}{9400}, \tfrac{233}{9400}\big) \approx (0.5057, 0.7837, 0.0248)$, up from $(0.5, 0.8, 0.02)$ in $c_0$ and $c_2$ and down in $c_1$. The last reading (which had residual $+0.05$) is pulled toward the curve: the residuals become $(-0.0557, 0.0757, 0.0477, -0.0400, -0.0872, 0.0060)$. The sixth shrinks from $+0.050$ to $+0.006$; the first, second, third and fifth grow in size (the fifth from $-0.070$ to $-0.087$); the fourth stays at $-0.040$, because the change in the fitted quadratic happens to vanish at $u = 3$ (it is $\tfrac{1}{9400}(54 - 153u + 45u^2)$, and $54 - 459 + 405 = 0$). The weighted residual is no longer orthogonal to the columns of $A$; it is orthogonal to them in the $W$-inner product, $A^\top W\mathbf{r} = \mathbf{0}$. The last cell of `project.py` asserts these values against the exact fractions and prints them.

## C. Challenge

**C1.** (a) Let $M$ be upper triangular and invertible (so every $m_{ii} \ne 0$). Column $j$ of $M^{-1}$ is the solution $\mathbf{x}$ of $M\mathbf{x} = \mathbf{e}_j$. Back substitution from the bottom: row $n$ reads $m_{nn}x_n = (\mathbf{e}_j)_n = 0$ when $j < n$, so $x_n = 0$. Row $n - 1$ reads $m_{n-1,n-1}x_{n-1} + m_{n-1,n}x_n = 0$ when $j < n - 1$, so $x_{n-1} = 0$. Continue upward: every $x_i$ with $i > j$ is $0$. So column $j$ of $M^{-1}$ has zeros below position $j$: $M^{-1}$ is upper triangular. ∎

(b) Let $M$ be upper triangular with $M^\top M = I$. Then $M^{-1} = M^\top$ is lower triangular, and by (a) also upper triangular, hence diagonal; so $M = (M^\top)^\top$ is diagonal too. A diagonal matrix with $M^\top M = I$ has $m_{ii}^2 = 1$, so $m_{ii} = \pm1$. ∎

(c) Suppose $A = Q_1R_1 = Q_2R_2$ (both reduced, no sign convention). As in the notes, $M = Q_2^\top Q_1 = R_2R_1^{-1}$ is upper triangular and orthogonal, so by (b) $M = D$ is diagonal with entries $\pm1$. Then $R_2 = DR_1$ and $Q_1 = Q_2D$ (from $Q_2^\top Q_1 = D$ and $Q_2Q_2^\top Q_1 = Q_1$, since $Q_2Q_2^\top$ projects onto $C(A) \ni$ columns of $Q_1$). So any two QR factorizations differ by flipping the signs of some columns of $Q$ together with the same rows of $R$. There are $2^n$ choices of $D$, hence exactly $2^n$ reduced QR factorizations, and exactly one with positive diagonal. ∎

(d)
```python
Q, R = np.linalg.qr(A)
D = np.diag(np.sign(np.diag(R)))     # entries ±1, so D @ D = I
Q, R = Q @ D, D @ R
```
It is valid because $(QD)(DR) = Q(D^2)R = QR = A$, $QD$ still has orthonormal columns ($(QD)^\top(QD) = DQ^\top QD = D^2 = I$), and $DR$ is still upper triangular (scaling rows by $\pm1$) with diagonal $|r_{ii}| > 0$. For Example 4.2, $D = \operatorname{diag}(-1, -1)$ turns NumPy's $R = \begin{pmatrix}-2 & -3\\0 & -\sqrt5\end{pmatrix}$ into $\begin{pmatrix}2 & 3\\0 & \sqrt5\end{pmatrix}$, the $R$ of Example 4.4(b).

**C2.** Inner products needed: $\langle x^a, x^b\rangle = \int_{-1}^1x^{a+b}dx = \frac{2}{a + b + 1}$ if $a + b$ is even, and $0$ if $a + b$ is odd (odd integrand on a symmetric interval). That parity rule is why "half the inner products were zero": every even-degree monomial is orthogonal to every odd-degree one, so Gram–Schmidt never mixes the two families.

$p_0 = 1$, $\langle p_0, p_0\rangle = 2$.
$p_1 = x - \frac{\langle x, 1\rangle}{2}\cdot1 = x$ (since $\langle x, 1\rangle = 0$), $\langle p_1, p_1\rangle = \tfrac23$.
$p_2 = x^2 - \frac{\langle x^2, 1\rangle}{\langle 1, 1\rangle}\cdot1 - \frac{\langle x^2, x\rangle}{\langle x, x\rangle}\,x = x^2 - \frac{2/3}{2} - 0 = x^2 - \tfrac13$. Norm: $\langle p_2, p_2\rangle = \int(x^4 - \tfrac23x^2 + \tfrac19) = \tfrac25 - \tfrac49 + \tfrac29 = \tfrac25 - \tfrac29 = \tfrac{8}{45}$.
$p_3 = x^3 - \frac{\langle x^3, 1\rangle}{2} - \frac{\langle x^3, x\rangle}{2/3}\,x - \frac{\langle x^3, p_2\rangle}{8/45}\,p_2 = x^3 - 0 - \frac{2/5}{2/3}x - 0 = x^3 - \tfrac35x$ (both zero inner products by parity).

Rescale so the value at $1$ is $1$: $p_2(1) = \tfrac23$, so $\tfrac32p_2 = \tfrac32x^2 - \tfrac12 = \tfrac12(3x^2 - 1) = P_2$ ✓. $p_3(1) = \tfrac25$, so $\tfrac52p_3 = \tfrac52x^3 - \tfrac32x = \tfrac12(5x^3 - 3x) = P_3$ ✓. (And $\|P_2\|^2 = \tfrac94\cdot\tfrac{8}{45} = \tfrac25 = \tfrac{2}{2\cdot2 + 1}$, the general $\|P_n\|^2 = \tfrac{2}{2n + 1}$.)

SymPy for degree 4:
```python
import sympy as sp
x = sp.symbols("x")
ip = lambda f, g: sp.integrate(f * g, (x, -1, 1))
ps = []
for k in range(5):
    v = x**k
    for p in ps:
        v -= ip(x**k, p) / ip(p, p) * p
    ps.append(sp.expand(v))
p4 = ps[4]                                   # x**4 - 6*x**2/7 + 3/35
P4 = sp.expand(p4 / p4.subs(x, 1))           # 35*x**4/8 - 15*x**2/4 + 3/8
assert sp.expand(P4 - (35*x**4 - 30*x**2 + 3)/8) == 0
```
The output is $p_4 = x^4 - \tfrac67x^2 + \tfrac{3}{35}$, and $p_4(1) = \tfrac{8}{35}$, so $\tfrac{35}{8}p_4 = \tfrac18(35x^4 - 30x^2 + 3) = P_4$ ✓. Gram–Schmidt on monomials with this inner product *is* the Legendre family; with the weight $e^{-x^2}$ on $\mathbb{R}$ it gives Hermite polynomials (QM-3's oscillator), and with $1/\sqrt{1 - x^2}$ it gives Chebyshev.

**C3.** (a) $y'(t) = -Ae^{-\gamma t}(\gamma\cos\omega t + \omega\sin\omega t)$ (take $A > 0$). Write $\varphi = \arctan(\gamma/\omega)$, so $\sin\varphi = \gamma/\sqrt{\omega^2 + \gamma^2}$ and $\cos\varphi = \omega/\sqrt{\omega^2 + \gamma^2}$. The angle-addition formula gives $\gamma\cos\omega t + \omega\sin\omega t = \sqrt{\omega^2 + \gamma^2}\,\sin(\omega t + \varphi)$, so
$$y'(t) = -A\sqrt{\omega^2 + \gamma^2}\,e^{-\gamma t}\sin(\omega t + \varphi).$$
Setting $y' = 0$: $\sin(\omega t + \varphi) = 0$, equivalently $\tan\omega t = -\gamma/\omega$. The solutions are $\omega t_k = k\pi - \varphi$, spaced by $\pi$. Which are maxima? The sign of $y'$ is the sign of $-\sin(\omega t + \varphi)$. As $\omega t + \varphi$ increases through $k\pi$ with $k$ even, $\sin$ goes from negative to positive, so $y'$ goes from positive to negative: a local maximum. With $k$ odd, $\sin$ goes from positive to negative, $y'$ from negative to positive: a local minimum. So the maxima are the critical points with even $k$. At each, $\cos\omega t_k = \cos(k\pi - \varphi) = \cos\varphi = \dfrac{\omega}{\sqrt{\omega^2 + \gamma^2}}$ (using $\cos k\pi = 1$ and $\sin k\pi = 0$ for even $k$), the same number at every peak. So $y_k = y(t_k) = A'e^{-\gamma t_k}$ with $A' = A\omega/\sqrt{\omega^2 + \gamma^2}$, and $\ln y_k = \ln A' - \gamma t_k$ exactly: a straight line in $t_k$ with slope $-\gamma$. Design matrix columns $\mathbf{1}$ and $\mathbf{t}_{\text{peaks}}$; unknowns $(\ln A', -\gamma)$. (The classical log decrement $\delta = \ln(y_k/y_{k+1}) = \gamma\cdot2\pi/\omega$ is the same slope estimated from two peaks; least squares uses all of them.)

(b)
```python
A0, gamma, omega = 2.0, 0.15, 2 * np.pi
T = 2 * np.pi / omega
t = np.linspace(0, 20, 4001)
for noise in (0.01, 0.1):
    rng_c3 = np.random.default_rng(1)
    y = A0 * np.exp(-gamma * t) * np.cos(omega * t) + noise * rng_c3.standard_normal(t.size)
    tp, yp = [], []
    for k in range(int(t[-1] / T)):                  # the largest sample in each period
        win = (t >= k * T) & (t < (k + 1) * T)
        i = np.argmax(y[win]); tp.append(t[win][i]); yp.append(y[win][i])
    tp, yp = np.array(tp), np.array(yp)
    Ad = np.c_[np.ones(tp.size), tp]
    lnA, slope = least_squares_qr(Ad, np.log(yp))[0]
    excess = np.mean(yp - A0 * np.exp(-gamma * tp) * np.cos(omega * tp))
    print(f"noise {noise}: gamma_hat = {-slope:.3f} (true {gamma}); mean excess of picked peaks = {excess:.3f}")
```
With seed 1: $\hat\gamma = 0.145$ at $\sigma = 0.01$ and $\hat\gamma = 0.101$ at $\sigma = 0.1$. Both are low, the second badly.

The bias: near a maximum the signal is flat over tens of samples, so the largest of them is the one where the *noise* happened to be largest. That is biased upward: the picked value exceeds the true peak by about one to two $\sigma$ (the printed mean excess is $0.013$ at $\sigma = 0.01$ and $0.21$ at $\sigma = 0.1$). The excess is roughly the same for every peak, but in log space it lifts a peak of height $0.1$ by far more than a peak of height $2$. The late points are raised, the fitted line flattens, and $\hat\gamma$ comes out too small. Weighting by $y_k^2$ (B8) attacks variance, not bias, and barely helps: $0.148$ and $0.111$.

(c)
```python
    idx = np.round(np.arange(20) * T / (t[1] - t[0])).astype(int)   # samples at t_k = kT, where cos = 1
    tk, yk = t[idx], y[idx]
    keep = yk > 3 * noise                                            # drop peaks buried in the noise
    tk, yk = tk[keep], yk[keep]
    Ak = np.c_[np.ones(tk.size), tk]
    lnA_u, slope_u = least_squares_qr(Ak, np.log(yk))[0]
    lnA_w, slope_w = least_squares_weighted(Ak, np.log(yk), yk ** 2)
    print(f"noise {noise}: kept {keep.sum()} peaks; unweighted gamma_hat = {-slope_u:.3f}, "
          f"weighted gamma_hat = {-slope_w:.3f}, weighted intercept = {lnA_w:.3f} (ln A = {np.log(A0):.3f})")
```
(Indented to sit inside the loop of (b).) With seed 1: at $\sigma = 0.01$ all 20 samples are kept and $\hat\gamma = 0.152$ unweighted, $0.151$ weighted, intercept $0.695$ against $\ln A = 0.693$. At $\sigma = 0.1$ only 13 samples survive the $3\sigma$ cut; unweighted $\hat\gamma = 0.133$, weighted $\hat\gamma = 0.145$ with intercept $0.696$. (At $t_k = kT$ the intercept is $\ln A$, not $\ln A'$, because these are the cosine maxima rather than the true peaks; the difference is $\ln A' - \ln A = -0.0003$ here, invisible.)

Why $w_k = y_k^2$: if the sample is $y_k + \delta_k$ with $\delta_k$ of fixed standard deviation $\sigma$, then $\ln(y_k + \delta_k) = \ln y_k + \ln(1 + \delta_k/y_k) \approx \ln y_k + \delta_k/y_k$. The error in the fitted quantity has standard deviation $\sigma/y_k$: small peaks are noisy in log space. The right weight is the inverse variance (ST-3, Gauss–Markov with unequal variances), $w_k \propto y_k^2/\sigma^2$, i.e. $w_k = y_k^2$ up to a constant that does not change the minimizer. Scaling row $k$ by $\sqrt{w_k} = y_k$ is exactly B8.

**Bias or bad luck?** One seed cannot tell a biased estimator from an unlucky draw, so repeat the whole experiment over 4000 seeds and look at the mean and spread of $\hat\gamma$. The code below also includes the remedy discussed after the table.
```python
idx = np.round(np.arange(20) * T / (t[1] - t[0])).astype(int)
tk_all = t[idx]

def log_fits(y, noise):
    """Unweighted and y²-weighted log-line fits to the samples at t_k = kT that survive the 3σ cut."""
    yk = y[idx]; keep = yk > 3 * noise
    Ak = np.c_[np.ones(keep.sum()), tk_all[keep]]
    g_u = -least_squares_qr(Ak, np.log(yk[keep]))[0][1]
    g_w = -least_squares_weighted(Ak, np.log(yk[keep]), yk[keep] ** 2)[1]
    return g_u, g_w

def fit_exp_gn(tk, yk, A=1.0, g=0.1, iters=30):
    """Gauss–Newton for min Σ (y_k − A e^{−g t_k})²: linearize in (A, g), solve a linear least-squares problem, repeat."""
    for _ in range(iters):
        f = A * np.exp(-g * tk)
        J = np.c_[np.exp(-g * tk), -tk * f]            # ∂f/∂A, ∂f/∂g
        dA, dg = least_squares_qr(J, yk - f)[0]
        A, g = A + dA, g + dg
    return g

for noise in (0.01, 0.1):
    est = []
    for seed in range(4000):
        rng_mc = np.random.default_rng(seed)
        y = A0 * np.exp(-gamma * t) * np.cos(omega * t) + noise * rng_mc.standard_normal(t.size)
        est.append([*log_fits(y, noise), fit_exp_gn(tk_all, y[idx])])
    est = np.array(est)
    print(noise, est.mean(axis=0).round(4), est.std(axis=0).round(4))   # unweighted, weighted, Gauss–Newton
```
Mean $\pm$ standard deviation of $\hat\gamma$ (true value $0.15$):

| $\sigma$ | log fit, unweighted | log fit, weighted $y_k^2$ | Gauss–Newton, no cut |
|---|---|---|---|
| $0.01$ | $0.1501 \pm 0.0019$ | $0.1499 \pm 0.0008$ | $0.1500 \pm 0.0008$ |
| $0.1$ | $0.131 \pm 0.014$ | $0.138 \pm 0.009$ | $0.150 \pm 0.008$ |

At $\sigma = 0.01$, sampling at known times has removed the peak-picking bias of (b) completely, and weighting cuts the spread by more than half. At $\sigma = 0.1$ a new bias appears: the log fits come out low by about $1.3$ standard deviations, unweighted and weighted alike. That is not thirteen points being noisy; it is systematic. Its source is the $3\sigma$ cut. The clean signal at $t_k = 11, 12, 13$ is $0.38, 0.33, 0.28$, right at the threshold $0.3$, so whether a late sample survives depends on its noise, and the survivors are mostly the ones whose noise was positive. That lifts the late log-values, flattens the line, and lowers $\hat\gamma$: the same mechanism as in (b), now produced by truncation instead of by picking maxima. Weighting reduces the damage, because it gives the late points little say, but cannot remove it.

The remedy is to stop truncating, and that means stopping taking logs, since the log of a sample at or below zero does not exist. Fit $y_k \approx Ae^{-\gamma t_k}$ to all twenty samples directly, by nonlinear least squares. Gauss–Newton (ML-2, and §1's outer loop) does it: linearize the model in $(A, \gamma)$ about the current guess, solve the linear least-squares problem for the correction with this module's QR, and repeat. The last column of the table shows it unbiased at both noise levels, with a spread no larger than the weighted log fit's.

**C4.** (a) Claim: in exact arithmetic, after subtracting the components along $\mathbf{q}_1, \dots, \mathbf{q}_{i-1}$, the current vector is $\mathbf{v}^{(i)} = \mathbf{a}_j - \sum_{l<i}r_{lj}\mathbf{q}_l$ with every $r_{lj} = \mathbf{q}_l^\top\mathbf{a}_j$. Base case $i = 1$: $\mathbf{v}^{(1)} = \mathbf{a}_j$ and $r_{1j} = \mathbf{q}_1^\top\mathbf{a}_j$ in both versions. Inductive step: assume it for $i$. Modified computes $r_{ij} = \mathbf{q}_i^\top\mathbf{v}^{(i)} = \mathbf{q}_i^\top\mathbf{a}_j - \sum_{l<i}r_{lj}\,\mathbf{q}_i^\top\mathbf{q}_l = \mathbf{q}_i^\top\mathbf{a}_j$, because $\mathbf{q}_i^\top\mathbf{q}_l = 0$ for $l < i$. So $r_{ij}$ agrees with classical, and $\mathbf{v}^{(i+1)} = \mathbf{v}^{(i)} - r_{ij}\mathbf{q}_i$ has the claimed form. ∎ The only thing the proof used is $\mathbf{q}_i^\top\mathbf{q}_l = 0$ *exactly*. In floating point that inner product is $\sim\varepsilon$, not $0$, and classical never gets to correct for it.

(b) With $\operatorname{fl}(1 + \varepsilon^2) = 1$:
- $\mathbf{a}_1 = (1, \varepsilon, 0, 0)^\top$, $\|\mathbf{a}_1\| = \sqrt{1 + \varepsilon^2} \to 1$, so $\mathbf{q}_1 = (1, \varepsilon, 0, 0)^\top$.
- $\mathbf{a}_2 = (1, 0, \varepsilon, 0)^\top$: $r_{12} = \mathbf{q}_1^\top\mathbf{a}_2 = 1$. $\mathbf{v}_2 = \mathbf{a}_2 - \mathbf{q}_1 = (0, -\varepsilon, \varepsilon, 0)^\top$, $\|\mathbf{v}_2\| = \varepsilon\sqrt2$, $\mathbf{q}_2 = \tfrac{1}{\sqrt2}(0, -1, 1, 0)^\top$. (Same in both versions: only one subtraction.)
- $\mathbf{a}_3 = (1, 0, 0, \varepsilon)^\top$. **Classical:** $r_{13} = \mathbf{q}_1^\top\mathbf{a}_3 = 1$, $r_{23} = \mathbf{q}_2^\top\mathbf{a}_3 = 0$ (computed from the *original* $\mathbf{a}_3$, whose entries 2 and 3 are zero). $\mathbf{v}_3 = \mathbf{a}_3 - \mathbf{q}_1 - 0 = (0, -\varepsilon, 0, \varepsilon)^\top$, $\mathbf{q}_3 = \tfrac{1}{\sqrt2}(0, -1, 0, 1)^\top$. Then $\mathbf{q}_2^\top\mathbf{q}_3 = \tfrac12(0 + 1 + 0 + 0) = \tfrac12$. Not orthogonal at all.
  **Modified:** $r_{13} = 1$, $\mathbf{v} = \mathbf{a}_3 - \mathbf{q}_1 = (0, -\varepsilon, 0, \varepsilon)^\top$. Now $r_{23} = \mathbf{q}_2^\top\mathbf{v} = \tfrac{1}{\sqrt2}(0 + \varepsilon + 0 + 0) = \varepsilon/\sqrt2$, and $\mathbf{v} \leftarrow \mathbf{v} - \tfrac{\varepsilon}{\sqrt2}\mathbf{q}_2 = (0, -\varepsilon, 0, \varepsilon)^\top - \tfrac\varepsilon2(0, -1, 1, 0)^\top = (0, -\tfrac\varepsilon2, -\tfrac\varepsilon2, \varepsilon)^\top$, $\|\mathbf{v}\| = \varepsilon\sqrt{3/2}$, $\mathbf{q}_3 = \tfrac{1}{\sqrt6}(0, -1, -1, 2)^\top$. Then $\mathbf{q}_2^\top\mathbf{q}_3 = \tfrac{1}{\sqrt{12}}(0 + 1 - 1 + 0) = 0$ ✓, and $\mathbf{q}_1^\top\mathbf{q}_3 = -\varepsilon/\sqrt6$, which is as good as the data allow.

The line where classical fails is $r_{23} = \mathbf{q}_2^\top\mathbf{a}_3 = 0$. In exact arithmetic $\mathbf{q}_2^\top\mathbf{a}_3 = \mathbf{q}_2^\top(\mathbf{a}_3 - r_{13}\mathbf{q}_1)$ because $\mathbf{q}_2 \perp \mathbf{q}_1$; but the computed $\mathbf{q}_1$ is not a unit vector to within $\varepsilon^2$ and the computed $\mathbf{q}_2$ is only orthogonal to it to within $\varepsilon$, so the two expressions differ by $\varepsilon/\sqrt2$, which is the *entire* size of the component that needed removing. Classical throws away the information that modified keeps. Numerically (`eps = 1e-8`, the matrix above): `gram_schmidt_classical` gives $\mathbf{q}_2^\top\mathbf{q}_3 = 0.5$ and `gram_schmidt_modified` gives $\mathbf{q}_2^\top\mathbf{q}_3 \approx 10^{-16}$.

(c) Suppose the first pass returns $Q_1$ with $Q_1^\top Q_1 = I + E$. The columns of $Q_1$ span $C(A)$ (each is a combination of the $\mathbf{a}$'s), so a second pass is Gram–Schmidt on a matrix with the *same* column space but a much better condition number than $A$'s $1.6\times10^{13}$ for $H_{10}$. How much better depends on $E$. The eigenvalues of $Q_1^\top Q_1$ are the squared singular values of $Q_1$, and they lie between $1 - \|E\|$ and $1 + \|E\|$ (2-norm), so when $\|E\| < 1$,
$$\kappa(Q_1)^2 \le \frac{1 + \|E\|}{1 - \|E\|}.$$
The second pass subtracts from each column its components along the earlier columns, which are now the residual non-orthogonality $E$ rather than the whole of $\mathbf{a}_j$. By the §3.8 estimate with $\kappa(Q_1)$ in place of $\kappa(A)$, the second pass leaves $\|Q_2^\top Q_2 - I\| \lesssim c\,\kappa(Q_1)^2\varepsilon$. So two passes reach $\varepsilon$ only if $\kappa(Q_1) = O(1)$, which by the inequality above means the first pass must have lost much less than $1$.

For $H_{10}$ the first pass loses everything ($\|Q_1^\top Q_1 - I\| = 3.0$), so the inequality says nothing and $Q_1$ is itself badly conditioned: $\kappa(Q_1) = 2.3\times10^7$. The estimate then predicts $\kappa(Q_1)^2\varepsilon \approx (2.3\times10^7)^2 \times 1.1\times10^{-16} \approx 0.06$, of order $10^{-1}$, after two passes; observed is $3.2\times10^{-5}$ (the bound is pessimistic, as it was for $H_{10}$ in §3.8). Now $\kappa(Q_2) \approx 1.00003$, and the third pass gives $4\times10^{-16}$. So $H_{10}$ needs three passes, not two:
```python
Q, _ = gram_schmidt_classical(sla.hilbert(10))
for k in range(1, 4):
    print(k, np.linalg.norm(Q.T @ Q - np.eye(10), 2))    # 3.0e+00, 3.2e-05, 4.2e-16
    Q, _ = gram_schmidt_classical(Q)
```
$H_{12}$ ($\kappa(Q_1) \approx 1.6\times10^{11}$) needs four: $5.0 \to 1.0 \to 3.5\times10^{-13} \to 5\times10^{-16}$.

The slogan "twice is enough" is about a different arrangement. Kahan's version (published in Parlett, *The Symmetric Eigenvalue Problem*, 1980) reorthogonalizes a single vector against a set that is already orthonormal, with a test for when a second pass is needed. For a whole factorization the corresponding result is for CGS2: orthogonalize *each column twice* against the $\mathbf{q}$'s already accepted, before moving on to the next column. Giraud, Langou, Rozložník and van den Eshof (*Numerische Mathematik*, 2005) prove that CGS2 gives $\|Q^\top Q - I\| = O(\varepsilon)$ whenever $c(m, n)\,\kappa(A)\,\varepsilon < 1$, a condition on $\kappa$, not $\kappa^2$. On $H_{10}$ it gives $2.6\times10^{-16}$ in one sweep:
```python
def cgs2(A):
    """Classical Gram–Schmidt with each column orthogonalized twice before moving on."""
    A = np.asarray(A, float); Q = np.zeros_like(A)
    for j in range(A.shape[1]):
        v = A[:, j].copy()
        for _ in range(2):
            v -= Q[:, :j] @ (Q[:, :j].T @ v)
        Q[:, j] = v / np.linalg.norm(v)
    return Q

Q = cgs2(sla.hilbert(10))
print(np.linalg.norm(Q.T @ Q - np.eye(10), 2))    # 2.6e-16
```
Re-running whole passes, as above, is weaker, because the second pass starts from a $Q_1$ whose conditioning the first pass has already wrecked. Each pass costs as much as the first, which is why Householder, which never loses orthogonality in the first place, is the production method (LA-13). Note also $A = Q_1R_1 = (Q_2R_2)R_1$, so the product $R_2R_1$ of the passes is the $R$ of $A$: reorthogonalization buys a good $Q$, not a different factorization.
