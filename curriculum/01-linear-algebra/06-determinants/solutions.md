# LA-6 Solutions

## A. Warm-up

**A1.** $A = \begin{pmatrix} 2 & 1 & 0 \\ 4 & 3 & 1 \\ -2 & 0 & 5 \end{pmatrix}$.

*Elimination.* Pivot $2$ in row 1, no swap needed. $R_2 \leftarrow R_2 - 2R_1 = (0, 1, 1)$. $R_3 \leftarrow R_3 + R_1 = (0, 1, 5)$. Pivot $1$ in row 2. $R_3 \leftarrow R_3 - R_2 = (0, 0, 4)$. Pivots $2, 1, 4$, zero swaps: $\det A = (+1)(2)(1)(4) = 8$.

*Row 1 cofactors.* $\det A = 2\,C_{11} + 1\,C_{12} + 0\,C_{13}$ with $C_{11} = +\begin{vmatrix} 3 & 1 \\ 0 & 5\end{vmatrix} = 15$ and $C_{12} = -\begin{vmatrix} 4 & 1 \\ -2 & 5\end{vmatrix} = -(20 + 2) = -22$. So $\det A = 30 - 22 = 8$. ✓

*Column 3 cofactors.* $\det A = 0\,C_{13} + 1\,C_{23} + 5\,C_{33}$ with $C_{23} = (-1)^{2+3}\begin{vmatrix} 2 & 1 \\ -2 & 0\end{vmatrix} = -(0 + 2) = -2$ and $C_{33} = (-1)^{3+3}\begin{vmatrix} 2 & 1 \\ 4 & 3\end{vmatrix} = 6 - 4 = 2$. So $\det A = -2 + 10 = 8$. ✓

**A2.** With $\det A = 5$, $n = 3$:
(a) One swap: $-5$. (b) (D3) on row 2: $3\cdot 5 = 15$. (c) Theorem 3.3: unchanged, $5$. (d) $2^3\cdot 5 = 40$ (Corollary 3.10′(d)). (e) $5$ (Theorem 3.11). (f) $1/5$ (Corollary 3.10′(a)). (g) $25$. (h) $(-1)^3\cdot 5 = -5$. (i) The swap gives $-5$; scaling a column by $-2$ multiplies by $-2$ (column linearity, via Theorem 3.11): $(-5)(-2) = 10$.

**A3.** $(3, 1, 2)$: inversions $(3,1), (3,2)$; $N = 2$, $\operatorname{sgn} = +1$.
$(2, 1, 4, 3)$: inversions $(2,1), (4,3)$; $N = 2$, $\operatorname{sgn} = +1$.
$(4, 3, 1, 2)$: inversions $(4,3), (4,1), (4,2), (3,1), (3,2)$; $N = 5$, $\operatorname{sgn} = -1$. Its term: $-a_{14}a_{23}a_{31}a_{42}$ (row $i$ uses column $\sigma(i)$).
The term $a_{13}a_{21}a_{34}a_{42}$ has row $1 \to$ column $3$, $2 \to 1$, $3 \to 4$, $4 \to 2$: $\sigma = (3, 1, 4, 2)$. Inversions: $(3,1), (3,2), (4,2)$; $N = 3$; sign $-1$, so the term enters as $-a_{13}a_{21}a_{34}a_{42}$.
Exactly half the permutations are even: swapping the first two entries is a bijection between even and odd permutations (Lemma 3.8), so $12$ have sign $+1$.

**A4.** $\det A = 3\cdot 2 - 1\cdot 4 = 2$. Cramer: $x_1 = \dfrac{1}{2}\begin{vmatrix} 5 & 1 \\ 8 & 2\end{vmatrix} = \dfrac{10 - 8}{2} = 1$, $x_2 = \dfrac{1}{2}\begin{vmatrix} 3 & 5 \\ 4 & 8\end{vmatrix} = \dfrac{24 - 20}{2} = 2$. Check: $3(1) + 1(2) = 5$ ✓, $4(1) + 2(2) = 8$ ✓.

Cofactors: $C_{11} = +2$, $C_{12} = -4$, $C_{21} = -1$, $C_{22} = +3$, so $C = \begin{pmatrix} 2 & -4 \\ -1 & 3\end{pmatrix}$ and $\operatorname{adj}A = C^\top = \begin{pmatrix} 2 & -1 \\ -4 & 3\end{pmatrix}$.
$A\operatorname{adj}A = \begin{pmatrix} 3 & 1 \\ 4 & 2\end{pmatrix}\begin{pmatrix} 2 & -1 \\ -4 & 3\end{pmatrix} = \begin{pmatrix} 6 - 4 & -3 + 3 \\ 8 - 8 & -4 + 6\end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2\end{pmatrix} = (\det A)I$. ✓
$A^{-1} = \tfrac12\begin{pmatrix} 2 & -1 \\ -4 & 3\end{pmatrix}$, and $A^{-1}\mathbf{b} = \tfrac12\begin{pmatrix} 10 - 8 \\ -20 + 24\end{pmatrix} = \begin{pmatrix} 1 \\ 2\end{pmatrix}$. ✓

**A5.** By (3.7), $\mathbf{u}\times\mathbf{v} = (u_2v_3 - u_3v_2,\ -(u_1v_3 - u_3v_1),\ u_1v_2 - u_2v_1) = (2\cdot 3 - 0\cdot 1,\ -(1\cdot 3 - 0\cdot 0),\ 1\cdot 1 - 2\cdot 0) = (6, -3, 1)^\top$.
Orthogonality: $(6, -3, 1)\cdot(1, 2, 0) = 6 - 6 + 0 = 0$ ✓; $(6, -3, 1)\cdot(0, 1, 3) = 0 - 3 + 3 = 0$ ✓.
Area $= \|\mathbf{u}\times\mathbf{v}\| = \sqrt{36 + 9 + 1} = \sqrt{46} \approx 6.782$.
Theorem 3.15(c): $\|\mathbf{u}\|^2\|\mathbf{v}\|^2 - (\mathbf{u}\cdot\mathbf{v})^2 = 5\cdot 10 - 2^2 = 46$ ✓.
Volume: $\mathbf{w}\cdot(\mathbf{u}\times\mathbf{v}) = 6 - 3 + 1 = 4$. Determinant with rows $\mathbf{w}, \mathbf{u}, \mathbf{v}$: $\begin{vmatrix} 1 & 1 & 1 \\ 1 & 2 & 0 \\ 0 & 1 & 3\end{vmatrix} = 1(6 - 0) - 1(3 - 0) + 1(1 - 0) = 4$ ✓. Volume $4$.

## B. Core

**B1.** (a) Fix $B$ and $D$ and define $f(A) = \det\begin{pmatrix} A & B \\ 0 & D\end{pmatrix}$ as a function of the $k\times k$ block $A$. Row $i \le k$ of $M$ is $(\mathbf{r}_i^\top, \mathbf{b}_i^\top)$, where $\mathbf{r}_i^\top$ is row $i$ of $A$ and $\mathbf{b}_i^\top$ is the (fixed) row $i$ of $B$.

*(D2) for $f$.* Swapping rows $i, j \le k$ of $A$ swaps rows $i, j$ of $M$ (the $B$ parts come along, but they are the same for $M$ and the swapped matrix only up to the swap, which is what we want: the new $M$ is the old $M$ with two rows swapped). So $f$ flips sign.

*(D3) for $f$.* If row $i$ of $A$ is $c\,\mathbf{r}^\top + \mathbf{r}'^\top$, then row $i$ of $M$ is $(c\,\mathbf{r}^\top + \mathbf{r}'^\top, \mathbf{b}_i^\top) = c(\mathbf{r}^\top, \mathbf{b}_i^\top) + (\mathbf{r}'^\top, \mathbf{b}_i^\top) - c(\mathbf{0}^\top, \mathbf{b}_i^\top)$... this does not split cleanly because $\mathbf{b}_i$ appears once, not $(1 + c)$ times. So we argue slightly differently: first kill $B$. Case 1: $D$ invertible. Using the rows of the $D$ block (which have zeros in the first $k$ columns), row operations of the form "add a combination of rows $k+1, \dots, k+m$ to row $i \le k$" can clear the $B$ block entirely (since the rows of $D$ span $\mathbb{R}^m$), and by Theorem 3.3 they do not change $\det M$. So $\det M = \det\begin{pmatrix} A & 0 \\ 0 & D\end{pmatrix} =: g(A)$. Now $g$ satisfies (D2) and (D3) in the rows of $A$ directly (row $i$ of the block matrix is $(\mathbf{r}_i^\top, \mathbf{0}^\top)$, which **is** linear in $\mathbf{r}_i$), so by the expansion argument of Theorem 3.9 steps 1–4 applied to the first $k$ rows, $g(A) = \det A\cdot g(I_k)$. Finally $g(I_k) = \det\begin{pmatrix} I_k & 0 \\ 0 & D\end{pmatrix}$, which by the same argument in the rows of $D$ equals $\det D\cdot\det I_{k+m} = \det D$. So $\det M = \det A\det D$.
Case 2: $D$ singular. Then the last $m$ rows of $M$, which are $(\mathbf{0}^\top, \mathbf{d}_i^\top)$, are linearly dependent (the rows of $D$ are, by Theorem 3.7 applied to $D$), so $\det M = 0$ by Theorem 3.5, and $\det A\det D = 0$ too. ∎

(b) Transpose: $\begin{pmatrix} A & 0 \\ C & D\end{pmatrix}^\top = \begin{pmatrix} A^\top & C^\top \\ 0 & D^\top\end{pmatrix}$, which by (a) has determinant $\det A^\top\det D^\top = \det A\det D$, and Theorem 3.11 says transposing changes nothing.

(c) Take $M = \begin{pmatrix} 0 & I_2 \\ I_2 & 0\end{pmatrix} = \begin{pmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0\end{pmatrix}$. This is the permutation matrix of $\sigma = (3, 4, 1, 2)$, with inversions $(3,1), (3,2), (4,1), (4,2)$, $N = 4$, so $\det M = +1$. (Or: swap rows $1 \leftrightarrow 3$ and $2 \leftrightarrow 4$ to reach $I$; two swaps, sign $+1$.) The false formula gives $\det 0\cdot\det 0 - \det I\cdot\det I = -1 \ne 1$.

**B2.** (a) Two matrices represent the same map $T$ in different bases exactly when they are related by $A' = S^{-1}AS$ with $S$ the change-of-basis matrix (LA-4). Corollary 3.10′(c) says $\det A' = \det A$. So "the determinant of (any matrix of) $T$" is a single well-defined number, and we call it $\det T$. (Corollary 3.16′ reaches the same number without ever picking a basis.)

(b) $Q^\top Q = I \Rightarrow \det(Q^\top)\det Q = 1 \Rightarrow (\det Q)^2 = 1$ by Theorems 3.10 and 3.11. So $\det Q = \pm 1$. For $R_\theta = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{pmatrix}$, $\det R_\theta = \cos^2\theta + \sin^2\theta = 1$.

(c) $\det P = \det(P^2) = (\det P)^2$, so $\det P \in \{0, 1\}$. If $\det P = 1$ then $P$ is invertible (Theorem 3.7); multiply $P^2 = P$ on the left by $P^{-1}$ to get $P = I$, contradicting $P \neq I$. So $\det P = 0$. (Geometrically: a projection onto a proper subspace flattens volume to zero.)

(d) Choose an orthonormal basis $\mathbf{u}, \mathbf{q}_2, \mathbf{q}_3$ of $\mathbb{R}^3$ (extend $\mathbf{u}$ by Gram–Schmidt, or any two orthonormal vectors in the plane $\mathbf{u}^\perp$). Then $H\mathbf{u} = \mathbf{u} - 2\mathbf{u}(\mathbf{u}^\top\mathbf{u}) = -\mathbf{u}$ and $H\mathbf{q}_i = \mathbf{q}_i - 2\mathbf{u}(\mathbf{u}^\top\mathbf{q}_i) = \mathbf{q}_i$. In this basis the matrix of $H$ is $\operatorname{diag}(-1, 1, 1)$, with determinant $-1$. By (a), $\det H = -1$ in every basis, including the standard one. (Consistent with (b): $H$ is orthogonal, so $\pm 1$ was the only possibility; the reflection takes the $-1$.)

**B3.** (a) *Multilinear:* $\mathbf{v}\times\mathbf{w}$ is linear in $\mathbf{v}$ and in $\mathbf{w}$ (each component in (3.7) is a $2\times 2$ determinant, linear in each row), and $\mathbf{u}\cdot\mathbf{z}$ is linear in $\mathbf{u}$ and in $\mathbf{z}$ (LA-1). Composing linear maps gives linearity in each of $\mathbf{u}, \mathbf{v}, \mathbf{w}$ separately. *Alternating:* if $\mathbf{v} = \mathbf{w}$ then $\mathbf{v}\times\mathbf{v} = \mathbf{0}$ (two equal rows in (3.7)); if $\mathbf{u} = \mathbf{v}$ then $\mathbf{v}\cdot(\mathbf{v}\times\mathbf{w}) = 0$ by Theorem 3.15(b); same if $\mathbf{u} = \mathbf{w}$. By Theorem 3.16, $f(\mathbf{u}, \mathbf{v}, \mathbf{w}) = f(\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3)\det(\mathbf{u}\ \mathbf{v}\ \mathbf{w})$, and $f(\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3) = \mathbf{e}_1\cdot(\mathbf{e}_2\times\mathbf{e}_3) = \mathbf{e}_1\cdot\mathbf{e}_1 = 1$. So $f = \det(\mathbf{u}\ \mathbf{v}\ \mathbf{w})$. The cyclic identity: $(\mathbf{v}, \mathbf{w}, \mathbf{u})$ is obtained from $(\mathbf{u}, \mathbf{v}, \mathbf{w})$ by two swaps (swap slots 1,2 to get $(\mathbf{v}, \mathbf{u}, \mathbf{w})$, then slots 2,3), so the determinant is multiplied by $(-1)^2 = 1$. Same for $(\mathbf{w}, \mathbf{u}, \mathbf{v})$.

(b) The dot product is bilinear (LA-1, §3.5) but $\mathbf{u}\cdot\mathbf{u} = \|\mathbf{u}\|^2 \ne 0$ for $\mathbf{u} \ne \mathbf{0}$, so it is not alternating. The function $g(\mathbf{u}, \mathbf{v}) = \|\mathbf{u}\times\mathbf{v}\|$ has $g(\mathbf{u}, \mathbf{u}) = 0$, but $g(-\mathbf{u}, \mathbf{v}) = \|-\mathbf{u}\times\mathbf{v}\| = g(\mathbf{u}, \mathbf{v}) \ne -g(\mathbf{u}, \mathbf{v})$ when $g \ne 0$, so it is not linear in the first slot. (It is $|\det|$ of a $2\times 2$ in the plane of the vectors, and absolute values break linearity.)

(c) By Theorem 3.16 with $n = 2$, $f(\mathbf{u}, \mathbf{v}) = c\det(\mathbf{u}\ \mathbf{v})$ with $c = f(\mathbf{e}_1, \mathbf{e}_2)$. From the data: $4 = c\begin{vmatrix} 1 & 1 \\ 1 & -1\end{vmatrix} = c(-1 - 1) = -2c$, so $c = -2$. Then $f\big((3,0)^\top, (0,2)^\top\big) = -2\begin{vmatrix} 3 & 0 \\ 0 & 2\end{vmatrix} = -2\cdot 6 = -12$.

**B4.** (a) By (3.11) with $x_i = i - 1$: $\det V = \prod_{1 \le i < j \le n}(j - i)$. Group by $j$: for fixed $j$, the factors $(j - i)$ for $i = 1, \dots, j-1$ are $j-1, j-2, \dots, 1$, with product $(j-1)!$. So $\det V = \prod_{j=2}^n (j-1)! = \prod_{j=1}^{n-1} j!$. For $n = 4$: $1!\cdot 2!\cdot 3! = 12$. Term by term: nodes $0, 1, 2, 3$, factors $(1-0)(2-0)(3-0)(2-1)(3-1)(3-2) = 1\cdot 2\cdot 3\cdot 1\cdot 2\cdot 1 = 12$. ✓

(b) Let $p(x) = c_0 + c_1x + \cdots + c_{n-1}x^{n-1}$ vanish at distinct $x_1, \dots, x_n$. Then $V\mathbf{c} = \mathbf{0}$ where $V$ is the Vandermonde matrix on these nodes and $\mathbf{c} = (c_0, \dots, c_{n-1})^\top$. Since the nodes are distinct, $\det V \neq 0$ by (3.11), so $V$ is invertible (Theorem 3.7) and $\mathbf{c} = V^{-1}\mathbf{0} = \mathbf{0}$. All coefficients vanish. ∎

(c) If $p$ and $q$ of degree $\le n - 1$ both pass through the $n$ points, $p - q$ has degree $\le n - 1$ and vanishes at the $n$ distinct $x_i$, so $p - q = 0$ by (b). (Existence: $V\mathbf{c} = \mathbf{y}$ is solvable because $V$ is invertible.)

(d) Cramer: $c_j = \det V_j/\det V$, with $V_j$ having column $j$ replaced by $\mathbf{y}$. As $x_2 \to x_1$, $\det V \to 0$ (rows 1 and 2 of $V$ become equal), continuously, since the determinant is a polynomial in the nodes. The numerators $\det V_j$ tend to the determinant of a matrix whose rows 1 and 2 agree **except** in column $j$, where they hold $y_1 \ne y_2$; such a matrix has rows 1 and 2 differing by $(y_1 - y_2)\mathbf{e}_j^\top$, and expanding shows $\det V_j \to (y_1 - y_2)\times(\text{a cofactor})$, which is nonzero for at least one $j$ (otherwise the limiting $V_j$ would all be singular, and one checks that the cofactor is a Vandermonde determinant on the remaining nodes, nonzero). So at least one $|c_j| \to \infty$. This is correct behavior: you are asking a polynomial to pass through two points with (nearly) the same $x$ and different $y$, which forces a (nearly) vertical segment and hence enormous coefficients; in the limit, no polynomial does it.

**B5.** (a) $\operatorname{vol}\Delta = \int_0^1\int_0^{1 - t_1}\int_0^{1 - t_1 - t_2}dt_3\,dt_2\,dt_1 = \int_0^1\int_0^{1 - t_1}(1 - t_1 - t_2)\,dt_2\,dt_1 = \int_0^1\frac{(1 - t_1)^2}{2}\,dt_1 = \frac{1}{2}\cdot\frac{1}{3} = \frac16$.

(b) Let $T$ have matrix $V = (\mathbf{v}_1\ \mathbf{v}_2\ \mathbf{v}_3)$. Then $T(t_1, t_2, t_3)^\top = t_1\mathbf{v}_1 + t_2\mathbf{v}_2 + t_3\mathbf{v}_3$ (column picture), so $T(\Delta) = \{t_1\mathbf{v}_1 + t_2\mathbf{v}_2 + t_3\mathbf{v}_3 : t_i \ge 0, \sum t_i \le 1\}$, which is the convex hull of $\mathbf{0}, \mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3$: the tetrahedron. By Theorem 3.17(c), $\operatorname{vol}T(\Delta) = |\det V|\operatorname{vol}\Delta = \tfrac16|\det V|$.

(c) Edge vectors from $P_0$: $(2, 0, 0), (0, 3, 0), (1, 1, 4)$. As columns, $V = \begin{pmatrix} 2 & 0 & 1 \\ 0 & 3 & 1 \\ 0 & 0 & 4\end{pmatrix}$, upper triangular, $\det V = 2\cdot 3\cdot 4 = 24$. Volume $= 24/6 = 4$. Check: the base triangle $P_0P_1P_2$ lies in the $xy$-plane with legs $2$ and $3$, area $3$; the height is the $z$-coordinate of $P_3$, which is $4$; $\tfrac13\cdot 3\cdot 4 = 4$. ✓

(d) Translate so $P_0$ is at the origin (translation preserves volume, (V2)): $\operatorname{vol} = \tfrac16|\det(P_1 - P_0\ \ P_2 - P_0\ \ P_3 - P_0)|$. For the $4\times 4$ form: subtract column 1 from columns 2, 3, 4 (no change in determinant, Theorem 3.3 for columns): the first row becomes $(1, 0, 0, 0)$ and the lower block becomes $(P_0\ \ P_1 - P_0\ \ P_2 - P_0\ \ P_3 - P_0)$. Expand along row 1: $\det = 1\cdot C_{11} = \det(P_1 - P_0\ \ P_2 - P_0\ \ P_3 - P_0)$. ✓

(e) Translating all vertices by $\mathbf{t}$ leaves every difference $P_i - P_0$ unchanged. Rotating by $R$ (orthogonal, $\det R = \pm 1$ by B2(b)) replaces $V$ by $RV$, and $|\det(RV)| = |\det R||\det V| = |\det V|$ (Theorem 3.10).

**B6.** (a) With $\mathbf{g}(r, \phi, \theta) = (r\sin\phi\cos\theta,\ r\sin\phi\sin\theta,\ r\cos\phi)$,
$$J = \begin{pmatrix} \sin\phi\cos\theta & r\cos\phi\cos\theta & -r\sin\phi\sin\theta \\ \sin\phi\sin\theta & r\cos\phi\sin\theta & r\sin\phi\cos\theta \\ \cos\phi & -r\sin\phi & 0\end{pmatrix}.$$
Expand along row 3: $\det J = \cos\phi\cdot C_{31} + (-r\sin\phi)\cdot C_{32} + 0\cdot C_{33}$.
$C_{31} = (+1)\begin{vmatrix} r\cos\phi\cos\theta & -r\sin\phi\sin\theta \\ r\cos\phi\sin\theta & r\sin\phi\cos\theta\end{vmatrix} = r^2\cos\phi\sin\phi\cos^2\theta + r^2\cos\phi\sin\phi\sin^2\theta = r^2\cos\phi\sin\phi$.
$C_{32} = (-1)\begin{vmatrix} \sin\phi\cos\theta & -r\sin\phi\sin\theta \\ \sin\phi\sin\theta & r\sin\phi\cos\theta\end{vmatrix} = -\left(r\sin^2\phi\cos^2\theta + r\sin^2\phi\sin^2\theta\right) = -r\sin^2\phi$.
So $\det J = r^2\cos^2\phi\sin\phi + (-r\sin\phi)(-r\sin^2\phi) = r^2\sin\phi(\cos^2\phi + \sin^2\phi) = r^2\sin\phi$.

(b) Column 1 is $\hat{\mathbf{r}} = (\sin\phi\cos\theta, \sin\phi\sin\theta, \cos\phi)$, length $\sqrt{\sin^2\phi + \cos^2\phi} = 1$. Column 2 is $r(\cos\phi\cos\theta, \cos\phi\sin\theta, -\sin\phi)$, length $r$. Column 3 is $r\sin\phi(-\sin\theta, \cos\theta, 0)$, length $r\sin\phi$. Dot products: $\mathbf{c}_1\cdot\mathbf{c}_2 = r(\sin\phi\cos\phi\cos^2\theta + \sin\phi\cos\phi\sin^2\theta - \cos\phi\sin\phi) = 0$; $\mathbf{c}_1\cdot\mathbf{c}_3 = r\sin^2\phi(-\cos\theta\sin\theta + \sin\theta\cos\theta) = 0$; $\mathbf{c}_2\cdot\mathbf{c}_3 = r^2\sin\phi\cos\phi(-\cos\theta\sin\theta + \sin\theta\cos\theta) = 0$. Three mutually orthogonal columns span a rectangular box, whose volume is the product of the edge lengths, $1\cdot r\cdot r\sin\phi$ (this is (V1) after a rotation that aligns the box with the axes; $|\det|$ is unchanged by the rotation, B5(e)). So $|\det J| = r^2\sin\phi$.

(c) For $0 < \phi < \pi$, $\sin\phi > 0$ and $r > 0$, so $\det J > 0$: the ordered triple $(\partial_r\mathbf{g}, \partial_\phi\mathbf{g}, \partial_\theta\mathbf{g})$ is right-handed. Ordering the variables $(r, \theta, \phi)$ swaps columns 2 and 3 and gives $\det J = -r^2\sin\phi$. For integration only $|\det J|$ enters (3.10), so the volume element is $r^2\sin\phi\,dr\,d\phi\,d\theta$ either way.

(d) $J = \begin{pmatrix} \cos\theta & -r\sin\theta & 0 \\ \sin\theta & r\cos\theta & 0 \\ 0 & 0 & 1\end{pmatrix}$ is block diagonal with blocks $\begin{pmatrix} \cos\theta & -r\sin\theta \\ \sin\theta & r\cos\theta\end{pmatrix}$ and $(1)$. By B1(a), $\det J = (r\cos^2\theta + r\sin^2\theta)\cdot 1 = r$.

**B7.** (a) From (3.6), $A\operatorname{adj}A = (\det A)I$. Take determinants: $\det A\cdot\det(\operatorname{adj}A) = \det\big((\det A)I\big) = (\det A)^n$ by Corollary 3.10′(d). Divide by $\det A \neq 0$: $\det(\operatorname{adj}A) = (\det A)^{n-1}$.

(b) Suppose $\det A = 0$ and, for contradiction, $\operatorname{adj}A$ is invertible. Then from $A\operatorname{adj}A = 0\cdot I = 0$, multiplying on the right by $(\operatorname{adj}A)^{-1}$ gives $A = 0$. But every cofactor of the zero matrix is the determinant of an $(n-1)\times(n-1)$ zero matrix, which is $0$ for $n \ge 2$, so $\operatorname{adj}(0) = 0$, which is not invertible. Contradiction. So $\operatorname{adj}A$ is singular, $\det(\operatorname{adj}A) = 0 = 0^{n-1}$.

(c) The minor of $A^\top$ obtained by deleting row $i$ and column $j$ is the transpose of the minor of $A$ obtained by deleting row $j$ and column $i$, so (Theorem 3.11) they have the same determinant: $(C^{A^\top})_{ij} = (-1)^{i+j}M^{A}_{ji} = C^{A}_{ji}$. Hence $C^{A^\top} = (C^A)^\top$ and $\operatorname{adj}(A^\top) = (C^{A^\top})^\top = C^A = (\operatorname{adj}A)^\top$. Apply (3.6) to $A^\top$: $A^\top\operatorname{adj}(A^\top) = (\det A^\top)I$, i.e., $A^\top(\operatorname{adj}A)^\top = (\det A)I$. Transpose both sides: $\operatorname{adj}(A)\,A = (\det A)I$. ∎

(d) ($\Leftarrow$) If $A$ is an integer matrix with $\det A = \pm 1$, every cofactor is a sum of products of integers, hence an integer, so $\operatorname{adj}A$ is an integer matrix and $A^{-1} = \operatorname{adj}A/\det A = \pm\operatorname{adj}A$ has integer entries. ($\Rightarrow$) If $A$ and $A^{-1}$ both have integer entries, then $\det A$ and $\det A^{-1}$ are both integers (the Leibniz formula is a sum of products of entries) and $\det A\cdot\det A^{-1} = 1$. The only integers whose product is $1$ are $\pm 1$. ∎

**B8.** Reference implementation (also in `project.py` under "Your turn"):

```python
def cofactor(A, i, j):
    """C_ij = (-1)^(i+j) * det(A with row i and column j deleted)."""
    minor = np.delete(np.delete(A, i, axis=0), j, axis=1)
    return (-1) ** (i + j) * det_elimination(minor)[0]

def adjugate(A):
    """adj A = C^T, so (adj A)[j, i] = C_ij.  Satisfies A @ adj A = det(A) I, eq. (3.6)."""
    A = np.asarray(A, float); n = A.shape[0]
    C = np.array([[cofactor(A, i, j) for j in range(n)] for i in range(n)])
    return C.T

for _ in range(20):
    A = rng.standard_normal((4, 4))
    assert np.allclose(A @ adjugate(A), np.linalg.det(A) * np.eye(4))

H = np.array([[1 / (i + j + 1) for j in range(4)] for i in range(4)])
H_exact_inv = np.array(sp.Matrix(4, 4, lambda i, j: sp.Rational(1, i + j + 1)).inv(), float)
err_adj = np.abs(adjugate(H) / det_elimination(H)[0] - H_exact_inv).max()
err_inv = np.abs(np.linalg.inv(H) - H_exact_inv).max()
print(f"kappa(H) = {np.linalg.cond(H):.3g}, adjugate error = {err_adj:.3g}, inv error = {err_inv:.3g}")
```

What you should see: $\kappa(H_4) \approx 1.55\times 10^4$ and both errors of order $10^{-11}$ to $10^{-12}$ against exact inverse entries as large as $6480$; relative to the entry size that is about $10^{-15}$–$10^{-16}$, i.e., $\kappa\cdot\varepsilon_{\text{machine}}$ divided by the entry scale. Neither method is reliably better: the adjugate formula does not escape the conditioning of the problem, because the error in $A^{-1}$ is governed by $\kappa(A)$ regardless of how you compute it (LA-13). What it does cost is $n^2$ determinants of size $n - 1$: by elimination that is $O(n^5)$ against $O(n^3)$ for solving $AX = I$ directly; by cofactors it is $O(n!)$. The adjugate is for proofs (B7) and symbolic $2\times 2$ and $3\times 3$ inverses, not for numbers.

## C. Challenge

**C1.** (a) By (3.3), $\det(I + tB) = \sum_\sigma\operatorname{sgn}(\sigma)\prod_{i=1}^n\big(\delta_{i\sigma(i)} + t\,b_{i\sigma(i)}\big)$, where $\delta_{ij}$ is $1$ if $i = j$ and $0$ otherwise. For $\sigma = \text{id}$ the product is $\prod_i(1 + tb_{ii}) = 1 + t\sum_i b_{ii} + O(t^2)$. For $\sigma \ne \text{id}$, at least two indices $i$ have $\sigma(i) \ne i$ (a permutation cannot move exactly one point), and for each such $i$ the factor is $0 + tb_{i\sigma(i)}$; so the product has a factor of $t^2$ and is $O(t^2)$. Summing: $\det(I + tB) = 1 + t\operatorname{tr}B + O(t^2)$, and differentiating at $t = 0$ gives $\operatorname{tr}B$.

(b) For small $h$, $A(t + h) = A(t) + hA'(t) + o(h) = A(t)\big(I + hA(t)^{-1}A'(t) + o(h)\big)$. Taking determinants with the product rule: $\det A(t + h) = \det A(t)\cdot\det\big(I + hA^{-1}A' + o(h)\big)$. The determinant is a polynomial in the entries, hence differentiable, so the $o(h)$ perturbation inside changes the determinant by $o(h)$, and by (a) with $B = A^{-1}A'$, $\det(I + hA^{-1}A') = 1 + h\operatorname{tr}(A^{-1}A') + O(h^2)$. Therefore $\det A(t + h) - \det A(t) = h\det A(t)\operatorname{tr}(A(t)^{-1}A'(t)) + o(h)$; divide by $h$ and let $h \to 0$. ∎

(c) $\operatorname{tr}(XY) = \sum_i(XY)_{ii} = \sum_i\sum_j x_{ij}y_{ji} = \sum_j\sum_i y_{ji}x_{ij} = \sum_j(YX)_{jj} = \operatorname{tr}(YX)$. ∎

(d) Let $y(t) = \det A(t)$. While $A(t)$ is invertible, (b) gives $y' = y\operatorname{tr}(A^{-1}A') = y\operatorname{tr}(A^{-1}BA) = y\operatorname{tr}(BAA^{-1}) = y\operatorname{tr}B$ by (c). With $y(0) = \det I = 1$, the scalar ODE $y' = (\operatorname{tr}B)y$ has the unique solution $y(t) = e^{t\operatorname{tr}B}$. This never vanishes, so $A(t)$ never becomes singular: the hypothesis "invertible" holds for all $t$ (if $A$ became singular at some first time $t^*$, $y$ would be continuous and nonzero on $[0, t^*)$ with limit $0$, contradicting $y = e^{t\operatorname{tr}B}$ on that interval). At $t = 1$: $\det e^{B} = e^{\operatorname{tr}B}$. ∎

(e) A linear flow $\dot{\mathbf{x}} = B\mathbf{x}$ moves a region $\Omega$ to $A(t)\Omega = e^{tB}\Omega$, whose volume is $e^{t\operatorname{tr}B}\operatorname{vol}\Omega$ by Theorem 3.17. If $\operatorname{tr}B = 0$, volume is preserved for all time. For a Hamiltonian system $\dot q = \partial H/\partial p$, $\dot p = -\partial H/\partial q$, the linearization has $\operatorname{tr}B = \partial^2H/\partial q\partial p - \partial^2H/\partial p\partial q = 0$, and Liouville's theorem (CM-8) is this computation applied locally and integrated along the flow.

**C2.** (a) $\det(A - \lambda I) = \sum_\sigma\operatorname{sgn}(\sigma)\prod_i\big(a_{i\sigma(i)} - \lambda\delta_{i\sigma(i)}\big)$. Each factor has degree $\le 1$ in $\lambda$, and degree exactly $1$ only when $\sigma(i) = i$. The identity permutation contributes $\prod_i(a_{ii} - \lambda) = (-\lambda)^n + (-\lambda)^{n-1}\sum_i a_{ii} + (\text{degree} \le n-2)$. Any other $\sigma$ fixes at most $n - 2$ indices (as in C1(a)), so contributes degree $\le n - 2$. Hence $p(\lambda) = (-1)^n\lambda^n + (-1)^{n-1}(\operatorname{tr}A)\lambda^{n-1} + \cdots$, degree exactly $n$. The constant term is $p(0) = \det A$.

(b) Expand $(-1)^n\prod_i(\lambda - \lambda_i) = (-1)^n\big(\lambda^n - (\sum_i\lambda_i)\lambda^{n-1} + \cdots + (-1)^n\prod_i\lambda_i\big)$. Matching the $\lambda^{n-1}$ coefficient with (a): $(-1)^{n+1}\sum_i\lambda_i = (-1)^{n-1}\operatorname{tr}A$, so $\sum_i\lambda_i = \operatorname{tr}A$. Matching constants: $(-1)^n(-1)^n\prod_i\lambda_i = \det A$, so $\prod_i\lambda_i = \det A$.

(c) For $A = \begin{pmatrix} a & b \\ c & d\end{pmatrix}$: $p(\lambda) = (a - \lambda)(d - \lambda) - bc = \lambda^2 - (a + d)\lambda + (ad - bc) = \lambda^2 - (\operatorname{tr}A)\lambda + \det A$. For $\begin{pmatrix} 2 & 1 \\ 1 & 2\end{pmatrix}$: $p(\lambda) = \lambda^2 - 4\lambda + 3 = (\lambda - 1)(\lambda - 3)$, roots $1$ and $3$. Sum $4 = \operatorname{tr}A$ ✓, product $3 = \det A$ ✓.

(d) $S^{-1}AS - \lambda I = S^{-1}AS - \lambda S^{-1}S = S^{-1}(A - \lambda I)S$, so $\det(S^{-1}AS - \lambda I) = \det(S^{-1})\det(A - \lambda I)\det S = \det(A - \lambda I)$ by Theorem 3.10 and Corollary 3.10′(a). Two matrices of the same linear map in different bases have the same characteristic polynomial, so its roots, the eigenvalues, belong to the map. Without this, "eigenvalue of $T$" would depend on the coordinates chosen to write $T$ down.

**C3.** (a) Write $A = QR$ (LA-5) with $Q$ orthogonal and $R$ upper triangular with diagonal entries $r_{jj}$ (which may be $0$ if $A$ is singular). Then $|\det A| = |\det Q||\det R| = 1\cdot\prod_j|r_{jj}|$ by B2(b) and Theorem 3.6. Column $j$ of $A = QR$ is $\mathbf{a}_j = \sum_{i \le j}r_{ij}\mathbf{q}_i$ (column picture, $R$ upper triangular), and since the $\mathbf{q}_i$ are orthonormal, $\|\mathbf{a}_j\|^2 = \sum_{i\le j}r_{ij}^2 \ge r_{jj}^2$. So $|r_{jj}| \le \|\mathbf{a}_j\|$ for every $j$, and multiplying, $|\det A| \le \prod_j\|\mathbf{a}_j\|$. Equality iff $r_{ij} = 0$ for all $i < j$ (or some $r_{jj} = 0$ with $\|\mathbf{a}_j\| = 0$, the zero-column case), i.e., $R$ diagonal, i.e., $\mathbf{a}_j = r_{jj}\mathbf{q}_j$ are mutually orthogonal. ∎

(b) $|\det A|$ is the volume of the parallelepiped with edges $\mathbf{a}_j$ (Theorem 3.17(a)); $\prod_j\|\mathbf{a}_j\|$ is the volume of a rectangular box with the same edge lengths. The inequality says: of all parallelepipeds with given edge lengths, the box has the largest volume, because tilting an edge away from perpendicular lowers its height above the face spanned by the others.

(c) If $k = n$, $\det(A^\top A) = \det A^\top\det A = (\det A)^2$. For $k = 2$, $G = \begin{pmatrix} \mathbf{a}_1\cdot\mathbf{a}_1 & \mathbf{a}_1\cdot\mathbf{a}_2 \\ \mathbf{a}_2\cdot\mathbf{a}_1 & \mathbf{a}_2\cdot\mathbf{a}_2\end{pmatrix}$, so $\det G = \|\mathbf{a}_1\|^2\|\mathbf{a}_2\|^2 - (\mathbf{a}_1\cdot\mathbf{a}_2)^2 = \|\mathbf{a}_1\times\mathbf{a}_2\|^2$ by Theorem 3.15(c), the squared area. In general, with the thin QR factorization $A = QR$ ($Q$ is $n\times k$ with orthonormal columns, $R$ is $k\times k$), $G = R^\top Q^\top QR = R^\top R$ and $\det G = (\det R)^2 = \prod_j r_{jj}^2$, where $|r_{jj}|$ is the height of $\mathbf{a}_j$ above the span of $\mathbf{a}_1, \dots, \mathbf{a}_{j-1}$ (Gram–Schmidt); the product of base times successive heights is the $k$-volume.

(d) Reference code (also in `project.py` under "Your turn"):

```python
# Hadamard: |det A| <= prod of column norms, equality for orthogonal columns.
for _ in range(1000):
    A = rng.standard_normal((5, 5))
    assert abs(np.linalg.det(A)) <= np.prod(np.linalg.norm(A, axis=0)) * (1 + 1e-12)
Q, _ = np.linalg.qr(rng.standard_normal((5, 5)))
A_orth = Q @ np.diag(rng.uniform(0.5, 2, 5))          # orthogonal columns, varied lengths
assert np.isclose(abs(np.linalg.det(A_orth)), np.prod(np.linalg.norm(A_orth, axis=0)))
# Gram determinant of a 3x2 matrix = squared area = |a1 x a2|^2.
A32 = rng.standard_normal((3, 2))
assert np.isclose(np.linalg.det(A32.T @ A32), np.linalg.norm(np.cross(A32[:, 0], A32[:, 1])) ** 2)
```
