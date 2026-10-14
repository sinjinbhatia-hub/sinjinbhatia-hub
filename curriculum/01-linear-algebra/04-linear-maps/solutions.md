# LA-4 Solutions

## A. Warm-up

**A1.** (a) **Linear.** $T(a(x,y) + b(x',y')) = T(ax + bx',\, ay + by') = \big((ax+bx') + 2(ay+by'),\; 3(ax+bx') - (ay+by'),\; ay+by'\big) = a(x+2y,\, 3x-y,\, y) + b(x'+2y',\, 3x'-y',\, y')$. Each output entry is a fixed linear combination of the inputs, which is all (3.1) asks.

(b) **Not linear.** $T(0, 0) = (1, 0) \neq (0, 0)$, but every linear map sends $\mathbf{0}$ to $\mathbf{0}$.

(c) **Linear.** $(ap + bq)(2) = a\,p(2) + b\,q(2)$ by the definition of adding and scaling polynomials.

(d) **Linear.** $x\,(ap + bq)' = x(ap' + bq') = a\,xp' + b\,xq'$, using linearity of differentiation and of multiplication by $x$.

(e) **Not linear.** $T(2A) = (2A)^\top(2A) = 4A^\top A = 4T(A) \neq 2T(A)$ whenever $A^\top A \neq 0$, e.g. $A = I$.

(f) **Linear.** $\operatorname{tr}(aA + bB) = a\operatorname{tr}A + b\operatorname{tr}B$ (sum of diagonal entries is linear), and multiplying a scalar by the fixed matrix $I$ is linear in the scalar.

**A2.** Columns are $T\mathbf{e}_1 = T(1, 0) = (1, 3, 0)$ and $T\mathbf{e}_2 = T(0, 1) = (2, -1, 1)$, so
$$[T] = \begin{pmatrix} 1 & 2 \\ 3 & -1 \\ 0 & 1\end{pmatrix}, \quad 3\times 2 \text{ (three output coordinates, two input coordinates).}$$
Directly: $T(1, -1) = (1 - 2,\; 3 + 1,\; -1) = (-1, 4, -1)$. Via (3.3): $[T](1, -1)^\top = 1\cdot(1, 3, 0) - 1\cdot(2, -1, 1) = (-1, 4, -1)$. ✓

Kernel: $T(x, y) = \mathbf{0}$ forces $y = 0$ (third entry) and then $x = 0$ (first entry). $\ker T = \{\mathbf{0}\}$, nullity $0$, $T$ injective. Image: the span of the two columns, which are not parallel, so a plane in $\mathbb{R}^3$; rank $2$. Rank–Nullity: $2 = 0 + 2$. ✓ Not surjective: the image is $2$-dimensional inside $\mathbb{R}^3$. Injective, not surjective, as Corollary 3.15 permits only because $2 \neq 3$.

**A3.** (a) Columns of $P$ are the new basis vectors: $P = \begin{pmatrix} 2 & 1 \\ 1 & 1\end{pmatrix}$, $\det P = 1$, $P^{-1} = \begin{pmatrix} 1 & -1 \\ -1 & 2\end{pmatrix}$. (Check: $PP^{-1} = \begin{pmatrix} 2-1 & -2+2 \\ 1-1 & -1+2\end{pmatrix} = I$.)

(b) Solving: $2a + b = 3$, $a + b = 2$, subtract: $a = 1$, then $b = 1$. So $[(3,2)]_{\mathcal{B}'} = (1, 1)^\top$. Matrix: new coordinates come from $P^{-1}$ (old in, new out): $P^{-1}(3, 2)^\top = (3 - 2,\; -3 + 4)^\top = (1, 1)^\top$. ✓ Check: $1\cdot(2,1) + 1\cdot(1,1) = (3, 2)$. ✓

(c) The vector with $\mathcal{B}'$-coordinates $(3, 2)$ is $3(2,1) + 2(1,1) = (8, 5)$, which is $P(3, 2)^\top = (6 + 2,\; 3 + 2)^\top$. Same matrix as (b), opposite question: $P$ converts new coordinates to old.

**A4.** Row reduce: row 2 is twice row 1, so the rref is $\begin{pmatrix} 1 & 2 & 3 \\ 0 & 0 & 0\end{pmatrix}$. Pivot column 1, free variables $x_2, x_3$. From $x_1 = -2x_2 - 3x_3$: kernel basis $(-2, 1, 0)^\top$ and $(-3, 0, 1)^\top$ (set one free variable to $1$ and the other to $0$). Check: $A(-2, 1, 0)^\top = (-2 + 2,\; -4 + 4)^\top = \mathbf{0}$, $A(-3, 0, 1)^\top = (-3 + 3,\; -6 + 6)^\top = \mathbf{0}$. ✓ Nullity $2$.

Image: all columns are multiples of $(1, 2)^\top$, so $\operatorname{im}A = \operatorname{span}\{(1, 2)^\top\}$, rank $1$. Rank–Nullity: $3 = 2 + 1$. ✓ Not injective (kernel is a plane), not surjective (image is a line in $\mathbb{R}^2$).

**A5.** (a)
$$R_\theta R_\varphi = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{pmatrix}\begin{pmatrix} \cos\varphi & -\sin\varphi \\ \sin\varphi & \cos\varphi\end{pmatrix} = \begin{pmatrix} \cos\theta\cos\varphi - \sin\theta\sin\varphi & -\cos\theta\sin\varphi - \sin\theta\cos\varphi \\ \sin\theta\cos\varphi + \cos\theta\sin\varphi & -\sin\theta\sin\varphi + \cos\theta\cos\varphi\end{pmatrix}.$$
By Theorem 3.7, $R_\theta R_\varphi$ is the matrix of "rotate by $\varphi$, then by $\theta$," which is the rotation by $\theta + \varphi$, whose matrix is $R_{\theta+\varphi}$ by §3.13. Comparing the $(1,1)$ and $(2,1)$ entries:
$$\cos(\theta + \varphi) = \cos\theta\cos\varphi - \sin\theta\sin\varphi, \qquad \sin(\theta+\varphi) = \sin\theta\cos\varphi + \cos\theta\sin\varphi.$$
(The other two entries repeat these.)

(b) With $c = \cos 2\varphi$, $s = \sin 2\varphi$: $F_\varphi^2 = \begin{pmatrix} c & s \\ s & -c\end{pmatrix}\begin{pmatrix} c & s \\ s & -c\end{pmatrix} = \begin{pmatrix} c^2 + s^2 & cs - sc \\ sc - cs & s^2 + c^2\end{pmatrix} = I$. Reflecting twice is the identity. $F_\varphi R_\theta$ rotates by $\theta$, then reflects. Multiply it out with $c = \cos 2\varphi$, $s = \sin 2\varphi$:
$$F_\varphi R_\theta = \begin{pmatrix} c & s \\ s & -c\end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{pmatrix} = \begin{pmatrix} c\cos\theta + s\sin\theta & -c\sin\theta + s\cos\theta \\ s\cos\theta - c\sin\theta & -s\sin\theta - c\cos\theta\end{pmatrix} = \begin{pmatrix} \cos(2\varphi - \theta) & \sin(2\varphi - \theta) \\ \sin(2\varphi - \theta) & -\cos(2\varphi - \theta)\end{pmatrix},$$
using the angle-subtraction formulas (part (a) with $-\theta$). That is $F_{\varphi - \theta/2}$: a reflection across the line at angle $\varphi - \theta/2$. Its determinant, by (3.8), is $-\cos^2(2\varphi-\theta) - \sin^2(2\varphi-\theta) = -1$: rotate-then-reflect reverses orientation.

**A6.** (a) $D(1) = 0$, $D(x) = 1$, $D(x^2) = 2x$:
$$[D] = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \\ 0 & 0 & 0\end{pmatrix}, \qquad [D]^2 = \begin{pmatrix} 0 & 0 & 2 \\ 0 & 0 & 0 \\ 0 & 0 & 0\end{pmatrix}, \qquad [D]^3 = 0.$$
$[D]^2(c_0, c_1, c_2)^\top = (2c_2, 0, 0)^\top$, which is zero iff $c_2 = 0$: the polynomials of degree $\le 1$, whose second derivative is zero. $D$ is nilpotent of index $3 = \dim\mathcal{P}_2$.

(b) $T(1) = x\cdot 0 = 0$, $T(x) = x\cdot 1 = x$, $T(x^2) = x\cdot 2x = 2x^2$. So $[T] = \operatorname{diag}(0, 1, 2)$. It is diagonal because each basis vector is sent to a multiple of itself: $T x^k = k\,x^k$. The monomials are *eigenvectors* of $T$ with eigenvalues $0, 1, 2$ (LA-7).

## B. Core

**B1.** (a) Suppose $(\mathbf{v}_i)$ is independent and $\sum_i c_i\,T\mathbf{v}_i = \mathbf{0}$. By linearity $T\big(\sum_i c_i\mathbf{v}_i\big) = \mathbf{0}$, so $\sum_i c_i\mathbf{v}_i \in \ker T = \{\mathbf{0}\}$ (injective, Theorem 3.13), so $\sum_i c_i\mathbf{v}_i = \mathbf{0}$, so all $c_i = 0$. Hence $(T\mathbf{v}_i)$ is independent. Conversely, suppose $(T\mathbf{v}_i)$ is independent and $\sum_i c_i\mathbf{v}_i = \mathbf{0}$. Apply $T$: $\sum_i c_i\,T\mathbf{v}_i = \mathbf{0}$, so all $c_i = 0$. (This direction uses only linearity, not injectivity.)

(b) Suppose $(\mathbf{v}_i)$ spans $V$ and let $\mathbf{w} \in W$. Since $T$ is surjective, $\mathbf{w} = T\mathbf{v}$ for some $\mathbf{v} = \sum_i c_i\mathbf{v}_i$, so $\mathbf{w} = \sum_i c_i\,T\mathbf{v}_i$: the $T\mathbf{v}_i$ span $W$. Conversely, suppose $(T\mathbf{v}_i)$ spans $W$ and let $\mathbf{v} \in V$. Then $T\mathbf{v} = \sum_i c_i\,T\mathbf{v}_i = T\big(\sum_i c_i\mathbf{v}_i\big)$ for some $c_i$, and injectivity gives $\mathbf{v} = \sum_i c_i\mathbf{v}_i$. So the $\mathbf{v}_i$ span $V$.

(c) $\kappa_{\mathcal{B}}$ is an isomorphism (Theorem 3.4), so by (a) polynomials are independent iff their coordinate vectors are. Coordinate vectors of $1 + x,\; 1 - x,\; 1 + x^2$ in $(1, x, x^2)$: $(1, 1, 0)$, $(1, -1, 0)$, $(1, 0, 1)$. Suppose $a(1,1,0) + b(1,-1,0) + c(1,0,1) = \mathbf{0}$: third entry gives $c = 0$; second gives $a = b$; first gives $a + b = 0$, so $a = b = 0$. Independent. (Three independent vectors in the $3$-dimensional $\mathcal{P}_2$: they form a basis.)

**B2.** (a) $J(1) = x$, $J(x) = x^2/2$, $J(x^2) = x^3/3$; $D(1) = 0$, $D(x) = 1$, $D(x^2) = 2x$, $D(x^3) = 3x^2$:
$$[J] = \begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & \tfrac12 & 0 \\ 0 & 0 & \tfrac13\end{pmatrix} \ (4\times 3), \qquad [D] = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3\end{pmatrix} \ (3\times 4).$$

(b) $[D][J]$: column $j$ is $[D]$ times column $j$ of $[J]$: $[D](0,1,0,0)^\top = (1,0,0)^\top$, $[D](0,0,\tfrac12,0)^\top = (0,1,0)^\top$, $[D](0,0,0,\tfrac13)^\top = (0,0,1)^\top$. So $[D][J] = I_3$: differentiating $\int_0^x p$ gives back $p$ (the fundamental theorem of calculus) for every $p \in \mathcal{P}_2$.

$[J][D]$: columns $[J](0,0,0)^\top = \mathbf{0}$, $[J](1,0,0)^\top = (0,1,0,0)^\top$, $[J](0,2,0)^\top = (0,0,1,0)^\top$, $[J](0,0,3)^\top = (0,0,0,1)^\top$. So $[J][D] = \operatorname{diag}(0, 1, 1, 1)$: $\int_0^x p'(t)\,dt = p(x) - p(0)$. It fixes every polynomial with $p(0) = 0$ (coordinates $(0, c_1, c_2, c_3)$) and kills the constants. $J\circ D \neq I$ because information ($p(0)$) was destroyed by $D$ and cannot be recovered.

(c) $\ker J$: $Jp = 0$ means $(0, c_0, c_1/2, c_2/3) = \mathbf{0}$, so $p = 0$. $\ker J = \{0\}$, nullity $0$. $\operatorname{im}J$ = span of the columns $= \{$polynomials in $\mathcal{P}_3$ with zero constant term$\}$, rank $3$. Rank–Nullity for $J$: $\dim\mathcal{P}_2 = 3 = 0 + 3$. ✓
$\ker D$ = constants, nullity $1$. $\operatorname{im}D$ = span of $(1,0,0), (0,2,0), (0,0,3)$ $= \mathcal{P}_2$, rank $3$. Rank–Nullity for $D$: $\dim\mathcal{P}_3 = 4 = 1 + 3$. ✓

(d) A two-sided inverse would make $T: \mathcal{P}_3 \to \mathcal{P}_2$ an isomorphism, forcing $\dim\mathcal{P}_3 = \dim\mathcal{P}_2$ by Theorem 3.17(c), i.e. $4 = 3$. Impossible. More concretely: any $T: \mathcal{P}_3 \to \mathcal{P}_2$ has $\operatorname{rank}T \le 3$, so $\dim\ker T = 4 - \operatorname{rank}T \ge 1$ by (3.7), so $T$ is never injective. $J$ has the left inverse $D$ (that is $[D][J] = I_3$). A right inverse $R$ would satisfy $J\circ R = I_{\mathcal{P}_3}$, making $J$ surjective; but $\operatorname{rank}J \le \dim\mathcal{P}_2 = 3 < 4 = \dim\mathcal{P}_3$, so $J$ is not surjective. The asymmetry is dimensional: a map into a bigger space can be injective (left-invertible) but never surjective (right-invertible).

**B3.** (a) $P = \begin{pmatrix} 2 & 1 \\ 1 & 1\end{pmatrix}$, $P^{-1} = \begin{pmatrix} 1 & -1 \\ -1 & 2\end{pmatrix}$ (A3). First $AP$: columns $A(2,1)^\top = (4 + 1,\; 2 + 2)^\top = (5, 4)^\top$ and $A(1,1)^\top = (3, 3)^\top$:
$$AP = \begin{pmatrix} 5 & 3 \\ 4 & 3\end{pmatrix}.$$
Then
$$P^{-1}AP = \begin{pmatrix} 1 & -1 \\ -1 & 2\end{pmatrix}\begin{pmatrix} 5 & 3 \\ 4 & 3\end{pmatrix} = \begin{pmatrix} 5 - 4 & 3 - 3 \\ -5 + 8 & -3 + 6\end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 3 & 3\end{pmatrix}.$$

(b) $T\mathbf{b}'_1 = A(2,1)^\top = (5, 4)^\top$. In $\mathcal{B}'$-coordinates: $a(2,1) + b(1,1) = (5,4)$ gives $2a + b = 5$, $a + b = 4$, so $a = 1$, $b = 3$: column $(1, 3)^\top$. ✓ $T\mathbf{b}'_2 = A(1,1)^\top = (3, 3)^\top = 3\cdot(1, 1) = 0\cdot\mathbf{b}'_1 + 3\cdot\mathbf{b}'_2$: column $(0, 3)^\top$. ✓ Matches (a).

(c) $AP^{-1}$: columns $A(1,-1)^\top = (2 - 1,\; 1 - 2)^\top = (1, -1)^\top$ and $A(-1, 2)^\top = (-2 + 2,\; -1 + 4)^\top = (0, 3)^\top$, so $AP^{-1} = \begin{pmatrix} 1 & 0 \\ -1 & 3\end{pmatrix}$. Then
$$PAP^{-1} = \begin{pmatrix} 2 & 1 \\ 1 & 1\end{pmatrix}\begin{pmatrix} 1 & 0 \\ -1 & 3\end{pmatrix} = \begin{pmatrix} 2 - 1 & 0 + 3 \\ 1 - 1 & 0 + 3\end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 0 & 3\end{pmatrix}.$$
Trace $4$ and determinant $3$, exactly the same as $\begin{pmatrix} 1 & 0 \\ 3 & 3\end{pmatrix}$ (trace $4$, det $3$) and as $A$ (trace $4$, det $3$). Invariants cannot catch a direction error, because $PAP^{-1}$ is *also* similar to $A$: it is the matrix of $T$ in the basis whose vectors (in old coordinates) are the columns of $P^{-1}$, namely $(1, -1)$ and $(-1, 2)$. Wrong basis, not wrong mathematics. The only reliable test is to apply $P$ to $\mathbf{e}_1$ and see whether you get the first new basis vector.

(d) $A(1,1)^\top = (3, 3)^\top = 3(1,1)^\top$ and $A(1,-1)^\top = (1, -1)^\top$. So with $\mathcal{B}'' = \big((1,1),(1,-1)\big)$, $Q = \begin{pmatrix} 1 & 1 \\ 1 & -1\end{pmatrix}$, $Q^{-1} = \tfrac12 Q$: $AQ = \begin{pmatrix} 3 & 1 \\ 3 & -1\end{pmatrix}$, and $Q^{-1}AQ = \tfrac12\begin{pmatrix} 1 & 1 \\ 1 & -1\end{pmatrix}\begin{pmatrix} 3 & 1 \\ 3 & -1\end{pmatrix} = \tfrac12\begin{pmatrix} 6 & 0 \\ 0 & 2\end{pmatrix} = \operatorname{diag}(3, 1)$. Trace $4$, det $3$ again. In this basis $T$ is "stretch by $3$ along $(1,1)$, leave $(1,-1)$ alone." LA-7 is about finding this basis systematically.

**B4.** (a) By Theorem 3.10 with $X' = XY$ and $Y' = Z$: $\operatorname{tr}((XY)Z) = \operatorname{tr}(Z(XY)) = \operatorname{tr}(ZXY)$. Example: $X = E_{11}$, $Y = E_{12}$, $Z = E_{21}$ (matrix units). $YZ = E_{12}E_{21} = E_{11}$, so $XYZ = E_{11}E_{11} = E_{11}$, trace $1$. $XZ = E_{11}E_{21} = 0$, so $YXZ = 0$, trace $0$. $1 \neq 0$. (And $ZXY = E_{21}E_{11}E_{12} = E_{21}E_{12} = E_{22}$, trace $1 = \operatorname{tr}(XYZ)$, as the cyclic rule promises.)

(b) $\operatorname{tr}I = 2 = \operatorname{tr}N$, $\det I = 1 = \det N$. But for any invertible $P$, $P^{-1}IP = P^{-1}P = I \neq N$. So $N$ is not similar to $I$: the only matrix similar to $I$ is $I$ itself. (Same for any scalar matrix $cI$.)

(c) Suppose $XY - YX = cI$ with $c \neq 0$ (the case $c = 1$ included). Take the trace of both sides. $\operatorname{tr}(XY - YX) = \operatorname{tr}(XY) - \operatorname{tr}(YX) = 0$ by linearity of the trace and Theorem 3.10, while $\operatorname{tr}(cI) = cn \neq 0$. Contradiction. The field does not matter: the proof of Theorem 3.10 only reorders a finite sum of products $x_{ik}y_{ki}$, which works the same for complex entries, so the argument covers $c = i\hbar$. So the canonical commutation relation needs infinite-dimensional operators (where the trace is not defined for all operators and the argument breaks).

(d) $\operatorname{tr}T := \operatorname{tr}[T]_{\mathcal{B}}$ for any basis $\mathcal{B}$; well defined by Theorem 3.10 since all such matrices are similar. $\operatorname{tr}D = 0$ (the matrix in Example 4.1 has zero diagonal in either basis). The projection of Example 4.3 has trace $\tfrac15 + \tfrac45 = 1$ (or, in the adapted basis, $1 + 0$). The transpose map: its matrix is the permutation matrix of B7(a), with diagonal entries $1, 0, 0, 1$, so the trace is $2$.

**B5.** (a) Suppose $AB = I$. *Step 1:* $B$ is injective. If $B\mathbf{x} = \mathbf{0}$, then $\mathbf{x} = I\mathbf{x} = (AB)\mathbf{x} = A(B\mathbf{x}) = A\mathbf{0} = \mathbf{0}$ (associativity, LA-1 §3.3). So $\ker B = \{\mathbf{0}\}$, and by **Theorem 3.13** the map $\mathbf{x} \mapsto B\mathbf{x}$ is injective. *Step 2:* by **Corollary 3.15** (same dimension $n$ in and out), it is also surjective, hence bijective, hence $B$ is invertible by Theorem 3.18 (a)$\iff$(b): there is $B^{-1}$ with $BB^{-1} = B^{-1}B = I$. *Step 3:* $A = AI = A(BB^{-1}) = (AB)B^{-1} = IB^{-1} = B^{-1}$ (associativity again). Therefore $BA = BB^{-1} = I$. ∎ The entire content is Step 2: without equal dimensions, injective does not give surjective, and the argument dies.

(b) $A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0\end{pmatrix}$, $B = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0\end{pmatrix}$. $AB = I_2$ (the first two columns of $B$ picked out). $BA = \operatorname{diag}(1, 1, 0) \neq I_3$. In general: $BA$ is $3\times 3$ and its image is contained in $\operatorname{im}B$, which has dimension at most $2$ ($B$ has $2$ columns), so $\operatorname{rank}(BA) \le 2 < 3 = \operatorname{rank}I_3$. Equivalently, by Rank–Nullity $A: \mathbb{R}^3 \to \mathbb{R}^2$ has a nonzero kernel ($3 = \dim\ker A + \operatorname{rank}A \le \dim\ker A + 2$), and any $\mathbf{z} \neq \mathbf{0}$ with $A\mathbf{z} = \mathbf{0}$ gives $BA\mathbf{z} = \mathbf{0} \neq \mathbf{z}$.

**B6.** (a) Let $T: V \to W$ and $S: W \to U$ be isomorphisms. $S\circ T$ is linear: $(S\circ T)(a\mathbf{u} + b\mathbf{v}) = S(aT\mathbf{u} + bT\mathbf{v}) = aST\mathbf{u} + bST\mathbf{v}$. It is bijective because compositions of bijections are bijections, with inverse $T^{-1}\circ S^{-1}$. For the inverse: $T^{-1}$ is linear by Theorem 3.17(a) and bijective (its inverse is $T$). Reflexivity of $\cong$ comes from $I_V$.

(b) $\Phi(c_0 + c_1x + c_2x^2 + c_3x^3) = \begin{pmatrix} c_0 & c_1 \\ c_2 & c_3\end{pmatrix}$, $\Phi^{-1}\begin{pmatrix} a & b \\ c & d\end{pmatrix} = a + bx + cx^2 + dx^3$. $\Phi$ is linear (coefficients add and scale entrywise), and $\Phi^{-1}$ undoes it in both orders, so $\Phi$ is bijective. Transport: $\Phi^{-1}$ gives $a + bx + cx^2 + dx^3$; $D$ gives $b + 2cx + 3dx^2$; $\Phi$ gives
$$(\Phi\circ D\circ\Phi^{-1})\begin{pmatrix} a & b \\ c & d\end{pmatrix} = \begin{pmatrix} b & 2c \\ 3d & 0\end{pmatrix}.$$
This operator on $\mathbb{R}^{2\times 2}$ is "differentiation in disguise": its fourth power is zero, its kernel is the span of $E_{11}$, its image is the matrices with zero $(2,2)$ entry.

(c) ($\Rightarrow$) If $T$ is an isomorphism and $(\mathbf{b}_j)$ is any basis of $V$, then by B1 the list $(T\mathbf{b}_j)$ is independent and spanning in $W$, i.e. a basis. ($\Leftarrow$) Suppose $(\mathbf{b}_1, \dots, \mathbf{b}_n)$ is a basis of $V$ and $(T\mathbf{b}_1, \dots, T\mathbf{b}_n)$ is a basis of $W$. Surjective: any $\mathbf{w} = \sum_j c_j T\mathbf{b}_j = T\big(\sum_j c_j\mathbf{b}_j\big)$. Injective: if $T\mathbf{v} = \mathbf{0}$ with $\mathbf{v} = \sum_j x_j\mathbf{b}_j$, then $\sum_j x_j T\mathbf{b}_j = \mathbf{0}$, so all $x_j = 0$ by independence of the $T\mathbf{b}_j$, so $\mathbf{v} = \mathbf{0}$. ∎ (Equivalently: the matrix of $T$ in these two bases is $I$.)

**B7.** (a) $T(E_{11}) = E_{11}$, $T(E_{12}) = E_{21}$, $T(E_{21}) = E_{12}$, $T(E_{22}) = E_{22}$. Columns are the coordinate vectors $(1,0,0,0)$, $(0,0,1,0)$, $(0,1,0,0)$, $(0,0,0,1)$:
$$[T] = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1\end{pmatrix},$$
a permutation matrix swapping slots $2$ and $3$. $[T]^2$: swapping twice is the identity, so $[T]^2 = I$ (multiply it out: each column of $[T]$ is sent back to the corresponding $\mathbf{e}_j$). Hence $T$ is invertible (its own inverse), $\ker T = \{0\}$, $\operatorname{im}T = \mathbb{R}^{2\times 2}$.

(b) $K(A) = A - T(A)$, so $[K] = I - [T]$ by Theorem 3.19 (the matrix construction is linear):
$$[K] = \begin{pmatrix} 0 & 0 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & 0 & 0\end{pmatrix}.$$
Kernel: $[K](a,b,c,d)^\top = (0,\; b - c,\; c - b,\; 0)^\top = \mathbf{0}$ iff $b = c$. So $\ker K = \{A : A = A^\top\}$, the symmetric matrices, with basis $E_{11}, E_{22}, E_{12} + E_{21}$; dimension $3$. Image: every output is $(b - c)(0, 1, -1, 0)^\top$, i.e. a multiple of $E_{12} - E_{21} = \begin{pmatrix} 0 & 1 \\ -1 & 0\end{pmatrix}$, the skew-symmetric matrices; dimension $1$. Rank–Nullity: $3 + 1 = 4$. ✓

(c) $S^2 = \tfrac14(I + T)^2 = \tfrac14(I + 2T + T^2) = \tfrac14(2I + 2T) = \tfrac12(I + T) = S$, using $T^2 = I$. Image: $S(A) = \tfrac12(A + A^\top)$ is always symmetric, and $S(A) = A$ for symmetric $A$, so $\operatorname{im}S$ = symmetric matrices (dimension $3$). Kernel: $A + A^\top = 0$, the skew-symmetric matrices (dimension $1$). $[S] = \tfrac12(I + [T])$ has diagonal entries $\tfrac12(1 + 1), \tfrac12(1 + 0), \tfrac12(1 + 0), \tfrac12(1 + 1) = 1, \tfrac12, \tfrac12, 1$, so $\operatorname{tr}[S] = 3 = \operatorname{rank}S$. C3(b) explains: a projection is $\operatorname{diag}(1,1,1,0)$ in an adapted basis, and trace is basis-independent.

**B8.** Reference implementation (also at the bottom of `project.py`, so you can compare):

```python
def is_linear(T, basis, trials=20, tol=1e-9):
    """Numerical test of (3.1) on random combinations of the basis vectors."""
    B = np.column_stack([np.asarray(b, float).ravel() for b in basis])
    for _ in range(trials):
        a, b = rng.standard_normal(2)
        u, v = B @ rng.standard_normal(B.shape[1]), B @ rng.standard_normal(B.shape[1])
        if not np.allclose(T(a * u + b * v), a * T(u) + b * T(v), atol=tol):
            return False
    return True


def kernel_basis(T, basis_in, basis_out):
    """Basis of ker T as ambient vectors: null space of [T], mapped back through basis_in."""
    M = matrix_of_map(T, basis_in, basis_out)
    B_in = np.column_stack([np.asarray(b, float).ravel() for b in basis_in])
    return [B_in @ c for c in null_space(M).T]          # c = coordinates in basis_in


def image_basis(T, basis_in, basis_out):
    """Basis of im T as ambient vectors: column space of [T], mapped back through basis_out."""
    M = matrix_of_map(T, basis_in, basis_out)
    B_out = np.column_stack([np.asarray(b, float).ravel() for b in basis_out])
    return [B_out @ c for c in orth(M).T]              # c = coordinates in basis_out
```

(a) `is_linear` returns `True` for `D` and for the transpose map and `False` for `square = lambda p: p**2`, which squares each coefficient (NumPy's `**` is entrywise; it is not the polynomial product). `matrix_of_map` returns a matrix for the non-linear map anyway, because it only ever evaluates $T$ on the four basis vectors and records the answers; nothing in Definition 3.5 checks linearity. For the coefficient-wise square, $T(\mathbf{e}_k) = \mathbf{e}_k$, so the matrix is $I$, which predicts $T(2\mathbf{e}_1) = 2\mathbf{e}_1$; the true value is $4\mathbf{e}_1$. The matrix is the matrix of the unique linear map agreeing with $T$ on the basis (Theorem 3.2), which is not $T$.

(b) For $D$ on $\mathcal{P}_3$: one kernel vector, proportional to $(1, 0, 0, 0)$ (the constants); three image vectors spanning the coefficient vectors with zero last entry ($\mathcal{P}_2$). For $A \mapsto A - A^\top$: three kernel vectors, each symmetric when reshaped; one image vector, skew-symmetric when reshaped. In both cases the asserts `np.allclose(T(k), 0)` for each kernel vector and `len(kernel) + len(image) == 4` pass.

(c) A null-space vector $\mathbf{c}$ of $[T]_{\mathcal{B},\mathcal{C}}$ is a $\mathcal{B}$-coordinate vector (it is an *input* to the matrix, (3.3)), so the vector it names is $\sum_j c_j\mathbf{b}_j$, built from `basis_in`. A column-space vector of $[T]_{\mathcal{B},\mathcal{C}}$ is an *output*, i.e. a $\mathcal{C}$-coordinate vector, and names $\sum_i c_i\mathbf{c}_i$, built from `basis_out`. Mixing them up is the same error as confusing $P$ with $P^{-1}$.

## C. Challenge

**C1.** (a) $\operatorname{im}(S\circ T) = S(\operatorname{im}T) \subseteq \operatorname{im}S$, so $\operatorname{rank}(S\circ T) \le \operatorname{rank}S$. For the other bound, let $S|_{\operatorname{im}T}$ be $S$ with its domain cut down to the subspace $\operatorname{im}T$. Then $\operatorname{im}(S\circ T) = \operatorname{im}(S|_{\operatorname{im}T})$, and Rank–Nullity applied to $S|_{\operatorname{im}T}$ says its rank is at most the dimension of its domain, which is $\operatorname{rank}T$. So $\operatorname{rank}(S\circ T) \le \operatorname{rank}T$. ∎

(b) Rank–Nullity for $S|_{\operatorname{im}T}: \operatorname{im}T \to U$:
$$\operatorname{rank}T = \dim\ker(S|_{\operatorname{im}T}) + \operatorname{rank}(S\circ T).$$
The kernel of the restriction is $\ker S \cap \operatorname{im}T \subseteq \ker S$, so $\dim\ker(S|_{\operatorname{im}T}) \le \dim\ker S$. Now Rank–Nullity for $S\circ T$ and for $T$:
$$\dim\ker(S\circ T) = \dim V - \operatorname{rank}(S\circ T) = \dim V - \operatorname{rank}T + \dim\ker(S|_{\operatorname{im}T}) = \dim\ker T + \dim\ker(S|_{\operatorname{im}T}) \le \dim\ker T + \dim\ker S. ∎$$

(c) From the displayed identity in (b), $\operatorname{rank}(S\circ T) = \operatorname{rank}T - \dim(\ker S\cap\operatorname{im}T) \ge \operatorname{rank}T - \dim\ker S = \operatorname{rank}T - (\dim W - \operatorname{rank}S)$, which is Sylvester's inequality. The displayed identity shows exactly when it is an equality: when $\ker S \subseteq \operatorname{im}T$ (so that $\dim(\ker S\cap\operatorname{im}T) = \dim\ker S$). Equality example: $S = T = I$ on $\mathbb{R}^n$, $n \ge n + n - n$; or $S = \begin{pmatrix} 1 & 0 \\ 0 & 0\end{pmatrix}$, $T = \begin{pmatrix} 0 & 0 \\ 0 & 1\end{pmatrix}$ on $\mathbb{R}^2$, where $\ker S = \operatorname{span}\mathbf{e}_2 = \operatorname{im}T$, $ST = 0$, and $0 \ge 1 + 1 - 2 = 0$. Strict example: $S = T = \begin{pmatrix} 1 & 0 \\ 0 & 0\end{pmatrix}$, where $\ker S = \operatorname{span}\mathbf{e}_2$ meets $\operatorname{im}T = \operatorname{span}\mathbf{e}_1$ only in $\mathbf{0}$, so $ST = S$ has rank $1 > 1 + 1 - 2 = 0$. ✓

**C2.** (a) If $N^j\mathbf{v} = \mathbf{0}$ then $N^{j+1}\mathbf{v} = N(N^j\mathbf{v}) = N\mathbf{0} = \mathbf{0}$. So $\ker N^j \subseteq \ker N^{j+1}$.

(b) Suppose $\ker N^j = \ker N^{j+1}$ and let $\mathbf{v} \in \ker N^{j+2}$. Then $N^{j+1}(N\mathbf{v}) = \mathbf{0}$, so $N\mathbf{v} \in \ker N^{j+1} = \ker N^j$, so $N^j(N\mathbf{v}) = N^{j+1}\mathbf{v} = \mathbf{0}$, so $\mathbf{v} \in \ker N^{j+1}$. With (a), $\ker N^{j+1} = \ker N^{j+2}$. By induction the chain is constant from index $j$ on.

(c) Write $\ker N^0 = \ker I = \{\mathbf{0}\}$. $\ker N^k = V$ (since $N^k = 0$) and $\ker N^{k-1} \neq V$ (since $N^{k-1} \neq 0$), so $\ker N^{k-1} \subsetneq \ker N^k$. By (b), if the chain were ever constant at some step $j < k - 1$, i.e. $\ker N^j = \ker N^{j+1}$, it would stay constant forever after, contradicting $\ker N^{k-1} \neq \ker N^k$. Hence every step before $k$ is strict:
$$\{\mathbf{0}\} = \ker N^0 \subsetneq \ker N \subsetneq \ker N^2 \subsetneq \cdots \subsetneq \ker N^k = V.$$
Each strict inclusion of subspaces raises the dimension by at least $1$ (LA-3: a proper subspace has strictly smaller dimension), and there are $k$ of them starting from dimension $0$, so $n = \dim V \ge k$. ∎

(d) On $\mathcal{P}_n$, $D^{n+1} = 0$ (every derivative past the $n$-th of a degree-$\le n$ polynomial vanishes) and $D^n x^n = n! \neq 0$. Index $n + 1 = \dim\mathcal{P}_n$, the maximum allowed by (c). On $\mathbb{R}^3$: $N = E_{12}$, i.e. $N(x, y, z) = (y, 0, 0)$. $N^2(x,y,z) = N(y, 0, 0) = (0,0,0)$, and $N \neq 0$. Index $2$. Here the chain is $\{\mathbf{0}\} \subsetneq \ker N = \operatorname{span}(\mathbf{e}_1, \mathbf{e}_3) \subsetneq \mathbb{R}^3$, with a jump of $2$ in the first step.

(e) Multiply out: $(I - N)(I + N + \cdots + N^{k-1}) = I + N + \cdots + N^{k-1} - N - N^2 - \cdots - N^k = I - N^k = I$. In the other order, $(I + N + \cdots + N^{k-1})(I - N) = I + N + \cdots + N^{k-1} - N - N^2 - \cdots - N^k = I - N^k = I$, since $N^jN = N^{j+1}$ term by term. Both products are $I$, so $I - N$ is invertible with the stated inverse. (This is the geometric series; it converges trivially because it terminates.)

**C3.** (a) Existence: $\mathbf{v} = (\mathbf{v} - P\mathbf{v}) + P\mathbf{v}$, where $P(\mathbf{v} - P\mathbf{v}) = P\mathbf{v} - P^2\mathbf{v} = \mathbf{0}$ puts the first piece in $\ker P$ and the second is in $\operatorname{im}P$ by definition. Uniqueness: suppose $\mathbf{v} = \mathbf{k} + \mathbf{w} = \mathbf{k}' + \mathbf{w}'$ with $\mathbf{k}, \mathbf{k}' \in \ker P$ and $\mathbf{w}, \mathbf{w}' \in \operatorname{im}P$. Apply $P$: $P\mathbf{v} = P\mathbf{w} = P\mathbf{w}'$. But $P$ fixes its image ($\mathbf{w} = P\mathbf{u}$ gives $P\mathbf{w} = P^2\mathbf{u} = P\mathbf{u} = \mathbf{w}$), so $\mathbf{w} = \mathbf{w}'$, and then $\mathbf{k} = \mathbf{v} - \mathbf{w} = \mathbf{k}'$. ∎ (In particular $\ker P\cap\operatorname{im}P = \{\mathbf{0}\}$ and $\dim\ker P + \dim\operatorname{im}P = \dim V$, which is also Rank–Nullity.)

(b) Let $(\mathbf{w}_1, \dots, \mathbf{w}_r)$ be a basis of $\operatorname{im}P$ and $(\mathbf{k}_1, \dots, \mathbf{k}_s)$ a basis of $\ker P$. By (a) the combined list spans $V$ (every vector is image-part plus kernel-part) and is independent (a relation $\sum a_i\mathbf{w}_i + \sum b_j\mathbf{k}_j = \mathbf{0}$ writes $\mathbf{0}$ as image-part plus kernel-part; uniqueness in (a) forces both parts to be $\mathbf{0}$, and then the $a_i$ and $b_j$ vanish by independence of each list). So it is a basis, $r + s = n$. In it: $P\mathbf{w}_i = \mathbf{w}_i$ (fixed, as shown in (a)), so column $i$ is $\mathbf{e}_i$; $P\mathbf{k}_j = \mathbf{0}$, so column $r + j$ is $\mathbf{0}$. Hence $[P] = \operatorname{diag}(1, \dots, 1, 0, \dots, 0)$ with $r = \operatorname{rank}P$ ones. Its trace is $r$, and by Theorem 3.10 the trace is the same in every basis: $\operatorname{tr}P = \operatorname{rank}P$. ∎

(c) Given $V = U\oplus K$ (every $\mathbf{v}$ uniquely $\mathbf{u} + \mathbf{k}$), define $P\mathbf{v} = \mathbf{u}$. Well defined by uniqueness. Linear: if $\mathbf{v} = \mathbf{u} + \mathbf{k}$ and $\mathbf{v}' = \mathbf{u}' + \mathbf{k}'$, then $a\mathbf{v} + b\mathbf{v}' = (a\mathbf{u} + b\mathbf{u}') + (a\mathbf{k} + b\mathbf{k}')$ is the unique decomposition, so $P(a\mathbf{v} + b\mathbf{v}') = a\mathbf{u} + b\mathbf{u}' = aP\mathbf{v} + bP\mathbf{v}'$. $P^2 = P$: $P\mathbf{u} = \mathbf{u}$ for $\mathbf{u} \in U$ (its decomposition is $\mathbf{u} + \mathbf{0}$), so $P(P\mathbf{v}) = P\mathbf{u} = \mathbf{u} = P\mathbf{v}$. Image $= U$, kernel $= K$ by construction. Uniqueness: any projection $P'$ with image $U$ and kernel $K$ satisfies $P'\mathbf{u} = \mathbf{u}$ on $U$ (fixes its image) and $P'\mathbf{k} = \mathbf{0}$ on $K$, so $P'(\mathbf{u} + \mathbf{k}) = \mathbf{u} = P(\mathbf{u} + \mathbf{k})$. ∎

(d) $(I - P)^2 = I - 2P + P^2 = I - 2P + P = I - P$. $(I - P)\mathbf{v} = \mathbf{v} - P\mathbf{v}$ is the kernel-part of $\mathbf{v}$ from (a), so $\operatorname{im}(I - P) = \ker P$ and $\ker(I - P) = \{\mathbf{v} : P\mathbf{v} = \mathbf{v}\} = \operatorname{im}P$. The two projections split $V$ the same way with the roles of the two pieces swapped.

(e) In the basis $\mathcal{B}' = \big((1, 2), (1, 0)\big)$ (image direction, then kernel direction), $[P]_{\mathcal{B}'} = \operatorname{diag}(1, 0)$ by (b). Convert to the standard basis with (3.6) run backwards: $A = Q\,\operatorname{diag}(1, 0)\,Q^{-1}$ where $Q = \begin{pmatrix} 1 & 1 \\ 2 & 0\end{pmatrix}$ has the $\mathcal{B}'$ vectors as columns. $\det Q = -2$, $Q^{-1} = -\tfrac12\begin{pmatrix} 0 & -1 \\ -2 & 1\end{pmatrix} = \begin{pmatrix} 0 & \tfrac12 \\ 1 & -\tfrac12\end{pmatrix}$. $Q\operatorname{diag}(1,0) = \begin{pmatrix} 1 & 0 \\ 2 & 0\end{pmatrix}$, and
$$A = \begin{pmatrix} 1 & 0 \\ 2 & 0\end{pmatrix}\begin{pmatrix} 0 & \tfrac12 \\ 1 & -\tfrac12\end{pmatrix} = \begin{pmatrix} 0 & \tfrac12 \\ 0 & 1\end{pmatrix}.$$
Checks: $A(1, 2)^\top = (1, 2)^\top$ ✓ (image fixed); $A(1, 0)^\top = \mathbf{0}$ ✓ (kernel killed); $A^2 = \begin{pmatrix} 0 & \tfrac12 \\ 0 & 1\end{pmatrix}\begin{pmatrix} 0 & \tfrac12 \\ 0 & 1\end{pmatrix} = \begin{pmatrix} 0 & \tfrac12 \\ 0 & 1\end{pmatrix} = A$ ✓; $\operatorname{tr}A = 1 = \operatorname{rank}A$ ✓. $A \neq A^\top$. Geometrically: each point slides *horizontally* until it hits the line $y = 2x$ (so $(x, y) \mapsto (y/2, y)$), whereas Example 4.3 slides each point *perpendicularly* onto the line. Both have image $y = 2x$ and satisfy $P^2 = P$; they differ in the kernel ($(1, 0)$ versus $(2, -1)$). LA-5 singles out the perpendicular one by the extra condition $P = P^\top$ (equivalently, kernel $\perp$ image), and proves it is the one that minimizes the distance moved. That is the whole of least squares.
