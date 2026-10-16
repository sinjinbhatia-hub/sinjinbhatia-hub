# LA-8 Solutions

Conventions as in the notes: complex inner products conjugate the **first** slot, $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^\dagger\mathbf{v}$; eigenvalues ordered $\lambda_1\ge\cdots\ge\lambda_n$.

## A. Warm-up

**A1.** (a) $A^\top B = \begin{pmatrix}1&0\\2&-1\end{pmatrix}\begin{pmatrix}3&0\\1&1\end{pmatrix} = \begin{pmatrix}3&0\\5&-1\end{pmatrix}$, trace $3 + (-1) = 2$. Entrywise: $1\cdot3 + 2\cdot0 + 0\cdot1 + (-1)\cdot1 = 2$. ✓

(b) $\|A\|_F = \sqrt{1 + 4 + 0 + 1} = \sqrt6$, $\|B\|_F = \sqrt{9 + 0 + 1 + 1} = \sqrt{11}$, $\cos\theta = 2/\sqrt{66}\approx0.246$.

(c) *Conjugate symmetry:* $\operatorname{tr}(B^\dagger A) = \operatorname{tr}\big((A^\dagger B)^\dagger\big) = \overline{\operatorname{tr}(A^\dagger B)}$, because the trace of $M^\dagger$ is the conjugate of the trace of $M$ (diagonal entries get conjugated). *Linear in the second slot:* $\operatorname{tr}(A^\dagger(cB + C)) = c\operatorname{tr}(A^\dagger B) + \operatorname{tr}(A^\dagger C)$ since matrix multiplication and trace are linear. *Positive definite:* $\operatorname{tr}(A^\dagger A) = \sum_j(A^\dagger A)_{jj} = \sum_j\sum_i\bar a_{ij}a_{ij} = \sum_{i,j}|a_{ij}|^2\ge0$, and it is $0$ only if every $a_{ij} = 0$. ∎

**A2.** (a) $\langle x, x^2\rangle = \int_0^1x^3\,dx = \tfrac14$. $\|x\|^2 = \int_0^1x^2 = \tfrac13$, $\|x^2\|^2 = \int_0^1x^4 = \tfrac15$. $\cos\theta = \dfrac{1/4}{\sqrt{1/3}\sqrt{1/5}} = \dfrac{\sqrt{15}}{4}\approx0.968$, so $\theta\approx14.5°$.

(b) $|\langle x, x^2\rangle| = 0.25\le\sqrt{1/15}\approx0.258$. ✓ (Close to equality because $x$ and $x^2$ are nearly parallel on $[0,1]$.)

(c) $\langle1, x - c\rangle = \int_0^1(x - c)\,dx = \tfrac12 - c = 0$, so $c = \tfrac12$. This is the first Gram–Schmidt step on $[0,1]$; $x - \tfrac12$ is the first *shifted* Legendre polynomial up to scale.

**A3.** (a) Trace $4$, determinant $3$: $\lambda^2 - 4\lambda + 3 = 0$, $\lambda_1 = 3$, $\lambda_2 = 1$. $(A - 3I)\mathbf{x} = 0$ gives $x_1 = x_2$; $(A - I)\mathbf{x} = 0$ gives $x_1 = -x_2$. So
$$Q = \frac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix},\qquad\Lambda = \begin{pmatrix}3&0\\0&1\end{pmatrix}.$$
Check $Q^\top Q = \tfrac12\begin{pmatrix}2&0\\0&2\end{pmatrix} = I$.

(b) $P_1 = \mathbf{q}_1\mathbf{q}_1^\top = \tfrac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$, $P_2 = \tfrac12\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$. $P_1 + P_2 = I$ ✓. $P_1P_2 = \tfrac14\begin{pmatrix}1-1&-1+1\\1-1&-1+1\end{pmatrix} = 0$ ✓. $3P_1 + P_2 = \tfrac12\begin{pmatrix}4&2\\2&4\end{pmatrix} = A$ ✓.

(c) $\sqrt A = \sqrt3P_1 + P_2 = \tfrac12\begin{pmatrix}\sqrt3+1&\sqrt3-1\\\sqrt3-1&\sqrt3+1\end{pmatrix}$. Square it directly: diagonal entry $\tfrac14\big((\sqrt3+1)^2 + (\sqrt3-1)^2\big) = \tfrac14(4 + 2\sqrt3 + 4 - 2\sqrt3) = 2$; off-diagonal $\tfrac14\cdot2(\sqrt3+1)(\sqrt3-1) = \tfrac14\cdot2\cdot2 = 1$. So $(\sqrt A)^2 = A$ ✓.
$A^{-1} = \tfrac13P_1 + P_2 = \begin{pmatrix}\tfrac16+\tfrac12&\tfrac16-\tfrac12\\\tfrac16-\tfrac12&\tfrac16+\tfrac12\end{pmatrix} = \begin{pmatrix}\tfrac23&-\tfrac13\\-\tfrac13&\tfrac23\end{pmatrix}$. Check $AA^{-1}$: $(1,1)$ entry $2\cdot\tfrac23 - \tfrac13 = 1$, $(1,2)$ entry $-\tfrac23 + \tfrac23 = 0$. ✓
$e^A = e^3P_1 + eP_2 = \tfrac12\begin{pmatrix}e^3+e&e^3-e\\e^3-e&e^3+e\end{pmatrix}$.

**A4.** (a) $A = \begin{pmatrix}2&-1\\-1&2\end{pmatrix}$, $\lambda = 3, 1$ (eigenvectors $(1,-1)/\sqrt2$ and $(1,1)/\sqrt2$). Positive definite. Ellipse with semi-axis $1/\sqrt3\approx0.577$ along $(1,-1)$ and semi-axis $1$ along $(1,1)$. Check: the point $(1,1)/\sqrt2$ gives $2\cdot\tfrac12 - 2\cdot\tfrac12 + 2\cdot\tfrac12 = 1$. ✓

(b) $A = \begin{pmatrix}1&2\\2&1\end{pmatrix}$, $\lambda = 3, -1$. Indefinite. Hyperbola: in eigen-coordinates $3y_1^2 - y_2^2 = 1$, opening along $\mathbf{q}_1 = (1,1)/\sqrt2$ with vertices at distance $1/\sqrt3$, asymptotes where $q = 0$: $x^2 + 4xy + y^2 = 0$, i.e. $y/x = -2\pm\sqrt3$.

(c) $A = \begin{pmatrix}1&1\\1&1\end{pmatrix}$, $\lambda = 2, 0$. Positive semidefinite. $q = (x+y)^2 = 1$ is the pair of parallel lines $x + y = \pm1$, unbounded in the direction $\mathbf{q}_2 = (1,-1)/\sqrt2$ of the zero eigenvalue.

(d) $A = \begin{pmatrix}-1&2\\2&-5\end{pmatrix}$, trace $-6$, det $5 - 4 = 1 > 0$: both eigenvalues negative, $\lambda = -3\pm2\sqrt2\approx-0.17, -5.83$. Negative definite. $q\le0$ everywhere, so $q = 1$ is **empty**; $q = -1$ is an ellipse with semi-axes $1/\sqrt{0.17}\approx2.41$ and $1/\sqrt{5.83}\approx0.41$.

**A5.** (a) $S^\top = S$ and $N^\top = N$: both symmetric. $S^\dagger = \begin{pmatrix}0&-i\\-i&0\end{pmatrix} = -S\ne S$; $N^\dagger = \begin{pmatrix}1&-i\\-i&-1\end{pmatrix}\ne N$: neither is Hermitian. $\sigma_y^\top = \begin{pmatrix}0&i\\-i&0\end{pmatrix} = -\sigma_y$, not symmetric; $\sigma_y^\dagger = \sigma_y$ (shown in Proposition 3.22), Hermitian.

(b) $\det(S - \lambda I) = \lambda^2 - i\cdot i = \lambda^2 + 1$, so $\lambda = \pm i$. Not real. Theorem 3.9 does not apply because $S$ is not self-adjoint.

(c) $N^2 = \begin{pmatrix}1\cdot1 + i\cdot i & 1\cdot i + i\cdot(-1)\\ i\cdot1 + (-1)\cdot i & i\cdot i + (-1)(-1)\end{pmatrix} = \begin{pmatrix}0&0\\0&0\end{pmatrix}$. If $N\mathbf{v} = \lambda\mathbf{v}$ then $0 = N^2\mathbf{v} = \lambda^2\mathbf{v}$, so every eigenvalue is $0$. Eigenvectors solve $N\mathbf{v} = 0$: $v_1 + iv_2 = 0$, a one-dimensional space spanned by $(-i, 1)^\top$ (equivalently $(1, i)^\top$). One independent eigenvector for a $2\times2$ matrix: $N$ is not diagonalizable at all, let alone by a unitary. Theorem 3.12 (an orthonormal eigenbasis exists) fails for complex symmetric matrices; it needs $A^\dagger = A$.

(d) $\det(\sigma_y - \lambda I) = \lambda^2 - (-i)(i) = \lambda^2 - 1$, so $\lambda = \pm1$: real, as Theorem 3.9 guarantees for the Hermitian $\sigma_y$.

**A6.** $\mathbf{x}^\top A\mathbf{x} = 3x_1^2 + 2x_1x_2 + 3x_2^2$.
$R(1,0) = 3/1 = 3$. $R(1,1) = (3 + 2 + 3)/2 = 4$. $R(1,-1) = (3 - 2 + 3)/2 = 2$. $R(2,1) = (12 + 4 + 3)/5 = 19/5 = 3.8$.
Eigenvalues: trace $6$, det $8$, $\lambda = 4, 2$ with eigenvectors $(1,1)/\sqrt2$ and $(1,-1)/\sqrt2$. All four values lie in $[2, 4]$. ✓ Equality at $(1,1)$ and $(1,-1)$ because they are eigenvectors (Theorem 3.17). The other two are weighted averages: $(1,0) = \tfrac1{\sqrt2}(\mathbf{q}_1 + \mathbf{q}_2)$ has weights $w_1 = w_2 = \tfrac12$, giving $\tfrac12\cdot4 + \tfrac12\cdot2 = 3$; $(2,1)$ has $\mathbf{y} = Q^\top\mathbf{x} = (3/\sqrt2, 1/\sqrt2)$, weights $(\tfrac9{10}, \tfrac1{10})$, giving $0.9\cdot4 + 0.1\cdot2 = 3.8$. ✓

## B. Core

**B1.** (a) *Symmetric:* $\langle\mathbf{v},\mathbf{u}\rangle_W = \mathbf{v}^\top W\mathbf{u} = (\mathbf{v}^\top W\mathbf{u})^\top = \mathbf{u}^\top W^\top\mathbf{v} = \mathbf{u}^\top W\mathbf{v}$. *Bilinear:* matrix products are linear in each factor. *Positive definite:* $\mathbf{u}^\top W\mathbf{u} > 0$ for $\mathbf{u}\ne\mathbf{0}$ is Theorem 3.15 (ii), which holds because $W$ is PD. ∎

(b) $\langle A\mathbf{u},\mathbf{v}\rangle_W = (A\mathbf{u})^\top W\mathbf{v} = \mathbf{u}^\top A^\top W\mathbf{v}$ and $\langle\mathbf{u}, A^\star\mathbf{v}\rangle_W = \mathbf{u}^\top WA^\star\mathbf{v}$. These agree for all $\mathbf{u},\mathbf{v}$ iff the matrices agree (take $\mathbf{u} = \mathbf{e}_i$, $\mathbf{v} = \mathbf{e}_j$ to compare $(i,j)$ entries): $A^\top W = WA^\star$, so
$$A^\star = W^{-1}A^\top W.$$
With $W = I$ this is $A^\top$. ✓

(c) Self-adjoint iff $W^{-1}A^\top W = A$ iff $A^\top W = WA$ iff $(WA)^\top = WA$, i.e. $WA$ is symmetric. ∎

(d) $A = M^{-1}K$ has $MA = K = K^\top$, symmetric, so by (c) $A$ is self-adjoint for $\langle\cdot,\cdot\rangle_M$. And $K\mathbf{v} = \lambda M\mathbf{v}\iff M^{-1}K\mathbf{v} = \lambda\mathbf{v}\iff A\mathbf{v} = \lambda\mathbf{v}$, so the generalized eigenvalues are the roots of $\det(K - \lambda M) = \det M\,\det(A - \lambda I)$.

*Reality.* On the real space $(\mathbb{R}^n, \langle\cdot,\cdot\rangle_M)$, Theorem 3.9 says nothing about complex roots, because an eigenvalue of an operator on a real space is real by definition. So extend to $\mathbb{C}^n$ with $\langle\mathbf{u},\mathbf{v}\rangle_M = \mathbf{u}^\dagger M\mathbf{v}$. This is an inner product by B4(d): $M$ is real symmetric, hence Hermitian, and for $\mathbf{u} = \mathbf{a} + i\mathbf{b}$ with $\mathbf{a}, \mathbf{b}$ real, $\mathbf{u}^\dagger M\mathbf{u} = \mathbf{a}^\top M\mathbf{a} + \mathbf{b}^\top M\mathbf{b} + i(\mathbf{a}^\top M\mathbf{b} - \mathbf{b}^\top M\mathbf{a}) = \mathbf{a}^\top M\mathbf{a} + \mathbf{b}^\top M\mathbf{b} > 0$ for $\mathbf{u}\ne\mathbf{0}$ (the imaginary part cancels because $M$ is symmetric). For this complex inner product, using $M^\dagger = M$ and $K^\dagger = K$ ($K$ is real symmetric),
$$\langle A\mathbf{u},\mathbf{v}\rangle_M = \mathbf{u}^\dagger A^\dagger M\mathbf{v} = \mathbf{u}^\dagger(MA)^\dagger\mathbf{v} = \mathbf{u}^\dagger K\mathbf{v} = \mathbf{u}^\dagger MA\mathbf{v} = \langle\mathbf{u}, A\mathbf{v}\rangle_M,$$
so $A$ is self-adjoint on $(\mathbb{C}^n, \langle\cdot,\cdot\rangle_M)$. Every complex root $\lambda$ of $\det(A - \lambda I)$ has a nonzero eigenvector in $\mathbb{C}^n$, and Theorem 3.9, which was proved for any self-adjoint operator on any inner product space, makes $\lambda$ real. (C3(a)–(b) reach the same conclusion by a different route.)

*Orthogonality.* For real $\lambda$, $A - \lambda I$ is a real singular matrix, so the eigenvectors can be taken real, and Theorem 3.10 for $\langle\cdot,\cdot\rangle_M$ gives $\langle\mathbf{v}_i,\mathbf{v}_j\rangle_M = \mathbf{v}_i^\top M\mathbf{v}_j = 0$ when $\lambda_i\ne\lambda_j$. In CM-6, $\lambda = \omega^2$ real is why normal modes oscillate rather than grow or decay, and $M$-orthogonality is why the modes decouple the kinetic energy $\tfrac12\dot{\mathbf{x}}^\top M\dot{\mathbf{x}}$.

**B2.** (a) $\langle S, K\rangle_F = \operatorname{tr}(S^\top K) = \operatorname{tr}(SK)$. Transposing inside the trace and using cyclicity: $\operatorname{tr}(SK) = \operatorname{tr}((SK)^\top) = \operatorname{tr}(K^\top S^\top) = \operatorname{tr}(-KS) = -\operatorname{tr}(KS) = -\operatorname{tr}(SK)$. So $\operatorname{tr}(SK) = 0$. ∎

(b) $A = \tfrac12(A + A^\top) + \tfrac12(A - A^\top)$; the first term is symmetric, the second antisymmetric. Uniqueness: if $A = S + K = S' + K'$ then $S - S' = K' - K$ is both symmetric and antisymmetric, so it equals its own negative, so it is $0$.

(c) A symmetric matrix is determined by its entries on and above the diagonal: $\dim\mathcal S = n(n+1)/2$. An antisymmetric matrix has zero diagonal and is determined by the entries strictly above it: $\dim\mathcal K = n(n-1)/2$. Sum: $n^2 = \dim\mathbb{R}^{n\times n}$. By (a), $\mathcal K\subseteq\mathcal S^\perp$. By LA-5, $\dim\mathcal S^\perp = n^2 - \dim\mathcal S = n(n-1)/2 = \dim\mathcal K$. A subspace of $\mathcal S^\perp$ with the same dimension is all of it: $\mathcal K = \mathcal S^\perp$. ∎

(d) For any symmetric $M$: $A - M = (A - S) + (S - M)$ with $A - S = K\in\mathcal S^\perp$ and $S - M\in\mathcal S$. These are orthogonal, so by Pythagoras $\|A - M\|_F^2 = \|A - S\|_F^2 + \|S - M\|_F^2\ge\|A - S\|_F^2$, with equality iff $M = S$. $S$ is the orthogonal projection of $A$ onto $\mathcal S$. ∎ (Checked numerically in the project's "Your turn" 1.)

**B3.** (a) $\nabla f = (3x^2 - 3y^2,\ -6xy)$. $H = \begin{pmatrix}f_{xx}&f_{xy}\\f_{xy}&f_{yy}\end{pmatrix} = \begin{pmatrix}6x&-6y\\-6y&-6x\end{pmatrix}$. Critical points: $xy = 0$ and $x^2 = y^2$. If $x = 0$ then $y = 0$; if $y = 0$ then $x = 0$. The only critical point is $(0, 0)$.

(b) $H$ has trace $0$ and determinant $-36x^2 - 36y^2$, so its eigenvalues are $\pm6\sqrt{x^2 + y^2}$: indefinite at every point other than the origin. At the critical point $H(0,0) = 0$: both eigenvalues are $0$, so $H$ is semidefinite and singular, the one case in which Theorem 3.19 gives no conclusion. (A singular $H$ with eigenvalues of both signs would still give a saddle.)

(c) Along the $x$-axis $f(t, 0) = t^3$, which is positive for $t > 0$ and negative for $t < 0$. So $f$ takes values above and below $f(0,0) = 0$ arbitrarily close to the origin: neither a maximum nor a minimum, hence a saddle. In polar coordinates $f = r^3\cos3\varphi$ (the real part of $(x + iy)^3$): three ridges and three valleys alternate around the origin, the "monkey saddle" (two legs and a tail).

(d) $\nabla h = (3x^2 - 3,\ 2y) = \mathbf{0}$ gives $x = \pm1$, $y = 0$. $H = \begin{pmatrix}6x&0\\0&2\end{pmatrix}$. At $(1, 0)$: $H = \operatorname{diag}(6, 2)$, eigenvalues $6, 2 > 0$, positive definite, strict local minimum, $h = -2$. At $(-1, 0)$: $H = \operatorname{diag}(-6, 2)$, eigenvalues $-6, 2$, indefinite, saddle, $h = 2$ ($h$ decreases along $x$, increases along $y$). The textbook rule agrees: $\det H = 12x$ is positive at $x = 1$ with $h_{xx} = 6 > 0$ (minimum) and negative at $x = -1$ (saddle).

**B4.** (a) Write $\mathbf{u} = \sum_iu_i\mathbf{e}_i$ and $\mathbf{v} = \sum_jv_j\mathbf{e}_j$. Bilinearity gives $\langle\mathbf{u},\mathbf{v}\rangle = \sum_{i,j}u_iv_j\langle\mathbf{e}_i,\mathbf{e}_j\rangle = \sum_{i,j}u_iG_{ij}v_j = \mathbf{u}^\top G\mathbf{v}$. ∎

(b) $G_{ij} = \langle\mathbf{e}_i,\mathbf{e}_j\rangle = \langle\mathbf{e}_j,\mathbf{e}_i\rangle = G_{ji}$: symmetric. $\mathbf{u}^\top G\mathbf{u} = \langle\mathbf{u},\mathbf{u}\rangle > 0$ for $\mathbf{u}\ne\mathbf{0}$: positive definite. ∎

(c) For symmetric PD $G$, $\mathbf{u}^\top G\mathbf{v}$ is symmetric ($G^\top = G$), bilinear, and $\mathbf{u}^\top G\mathbf{u} > 0$ for $\mathbf{u}\ne\mathbf{0}$ by Theorem 3.15 (ii). So inner products on $\mathbb{R}^n$ correspond exactly to symmetric positive definite matrices.

(d) *Complex version.* Let $\langle\cdot,\cdot\rangle$ be an inner product on $\mathbb{C}^n$ (Definition 3.2) and $G_{ij} = \langle\mathbf{e}_i,\mathbf{e}_j\rangle$. Then
$$\langle\mathbf{u},\mathbf{v}\rangle = \Big\langle\sum_iu_i\mathbf{e}_i, \sum_jv_j\mathbf{e}_j\Big\rangle = \sum_{i,j}\bar u_iv_jG_{ij} = \mathbf{u}^\dagger G\mathbf{v},$$
the conjugate coming from the **first** slot. $G$ is Hermitian: $G_{ji} = \langle\mathbf{e}_j,\mathbf{e}_i\rangle = \overline{\langle\mathbf{e}_i,\mathbf{e}_j\rangle} = \overline{G_{ij}}$. And $\mathbf{u}^\dagger G\mathbf{u} = \|\mathbf{u}\|^2 > 0$: positive definite. Conversely, let $G$ be Hermitian PD ($\mathbf{x}^\dagger G\mathbf{x} > 0$ for $\mathbf{x}\ne\mathbf{0}$, as defined after Theorem 3.14) and set $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^\dagger G\mathbf{v}$. It is linear in $\mathbf{v}$ because matrix multiplication is. It is conjugate symmetric: a $1\times1$ matrix's conjugate is its conjugate transpose, so $\overline{\langle\mathbf{v},\mathbf{u}\rangle} = (\mathbf{v}^\dagger G\mathbf{u})^\dagger = \mathbf{u}^\dagger G^\dagger\mathbf{v} = \mathbf{u}^\dagger G\mathbf{v} = \langle\mathbf{u},\mathbf{v}\rangle$. And $\langle\mathbf{u},\mathbf{u}\rangle = \mathbf{u}^\dagger G\mathbf{u} > 0$ for $\mathbf{u}\ne\mathbf{0}$. So it is an inner product in our convention, and inner products on $\mathbb{C}^n$ correspond exactly to Hermitian PD matrices. In Axler's convention the same computation produces $\mathbf{v}^\dagger G^\top\mathbf{u}$, i.e. the formula $\mathbf{u}^\top G\bar{\mathbf{v}}$, conjugating the second vector. ∎

**B5.** (a) Suppose $\mathbf{w} = A\mathbf{v} = \mathbf{0}$. Then $A^\top A\mathbf{v} = \mathbf{0}$, but $A^\top A\mathbf{v} = \lambda\mathbf{v}\ne\mathbf{0}$. Contradiction, so $\mathbf{w}\ne\mathbf{0}$. Then $AA^\top\mathbf{w} = A(A^\top A\mathbf{v}) = A(\lambda\mathbf{v}) = \lambda\mathbf{w}$. ∎

(b) Let $E = \{\mathbf{v} : A^\top A\mathbf{v} = \lambda\mathbf{v}\}$ and $F = \{\mathbf{w} : AA^\top\mathbf{w} = \lambda\mathbf{w}\}$. By (a), $\phi(\mathbf{v}) = A\mathbf{v}$ maps $E\to F$; it is linear, and injective because $A\mathbf{v} = \mathbf{0}$ with $\mathbf{v}\in E$ forces $\lambda\mathbf{v} = A^\top A\mathbf{v} = \mathbf{0}$, so $\mathbf{v} = \mathbf{0}$. By the same argument with $A$ and $A^\top$ swapped, $\psi(\mathbf{w}) = A^\top\mathbf{w}$ maps $F\to E$ injectively ($A^\top A(A^\top\mathbf{w}) = A^\top(\lambda\mathbf{w})$). Injective linear maps $E\to F$ and $F\to E$ give $\dim E\le\dim F\le\dim E$, so the eigenspaces have equal dimension. For the symmetric matrices $A^\top A$ and $AA^\top$ that is the same as equal multiplicity: writing $A^\top A = Q\Lambda Q^\top$ (Theorem 3.11), $\det(A^\top A - tI) = \prod_i(\lambda_i - t)$, so the multiplicity of $\lambda$ as a root counts the $i$ with $\lambda_i = \lambda$, and those $\mathbf{q}_i$ form a basis of $E$ (the Remark after Theorem 3.13, step 1). The same holds for $AA^\top$ and $F$. So the nonzero eigenvalues agree with multiplicity. ∎ (In fact $\psi\circ\phi(\mathbf{v}) = A^\top A\mathbf{v} = \lambda\mathbf{v}$, so $\psi/\lambda$ is the inverse of $\phi$.)

(c) If $A$ ($m\times n$, $m > n$) has independent columns, $A^\top A$ is $n\times n$ and positive definite (Theorem 3.16): $n$ positive eigenvalues, no zeros. $AA^\top$ is $m\times m$ with the same $n$ positive eigenvalues by (b), and rank $AA^\top = \operatorname{rank}A = n$, so its remaining $m - n$ eigenvalues are $0$. **$m - n$ more zeros.** The project's "Your turn" 3 checks a $5\times3$ example: three shared positive eigenvalues, two extra zeros, and the square roots of the shared eigenvalues equal `np.linalg.svd(A)`'s singular values.

**B6.** (a) $(e^{iA})^\dagger = \sum_i\overline{e^{i\lambda_i}}P_i^\dagger = \sum_ie^{-i\lambda_i}P_i$, using $P_i^\dagger = (\mathbf{q}_i\mathbf{q}_i^\dagger)^\dagger = P_i$ and $\overline{e^{i\lambda}} = e^{-i\lambda}$, **which requires $\lambda_i$ real** (for complex $\lambda$, $\overline{e^{i\lambda}} = e^{-i\bar\lambda}\ne e^{-i\lambda}$). Then
$$(e^{iA})^\dagger e^{iA} = \sum_{j,k}e^{-i\lambda_j}e^{i\lambda_k}P_jP_k = \sum_je^{-i\lambda_j}e^{i\lambda_j}P_j = \sum_jP_j = I,$$
using $P_jP_k = \delta_{jk}P_j$ and completeness. So $e^{iA}$ is unitary. ∎

(b) Let $H = \mathbf{n}\cdot\boldsymbol\sigma$, $H^2 = I$, and $P_\pm = \tfrac12(I\pm H)$. Then $P_\pm^2 = \tfrac14(I\pm2H + H^2) = \tfrac14(2I\pm2H) = P_\pm$; $P_+P_- = \tfrac14(I - H^2) = 0$; $P_+ + P_- = I$; each is Hermitian since $H$ is. And $HP_\pm = \tfrac12(H\pm H^2) = \tfrac12(H\pm I) = \pm P_\pm$, so the range of $P_\pm$ lies in the $\pm1$ eigenspace. Why these are *the* spectral projectors: write $E_\pm$ for the $\pm1$ eigenspaces. We have $\operatorname{range}(P_\pm)\subseteq E_\pm$, and since $P_+ + P_- = I$ every vector is $P_+\mathbf{x} + P_-\mathbf{x}$, so $\operatorname{range}(P_+) + \operatorname{range}(P_-) = \mathbb{C}^2$. Hence $\dim\operatorname{range}(P_+) + \dim\operatorname{range}(P_-)\ge 2\ge\dim E_+ + \dim E_-$ (the eigenspaces meet only in $\mathbf{0}$, Theorem 3.10), which forces $\operatorname{range}(P_\pm) = E_\pm$. An orthogonal projection is determined by its range (LA-5), so $P_\pm$ are the projections onto the eigenspaces, and $H = P_+ - P_-$ is (3.10) with equal eigenvalues grouped (the Remark after Theorem 3.13). ∎

(c) By (3.11) with $f(t) = e^{-i\theta t/2}$:
$$e^{-i\theta H/2} = e^{-i\theta/2}P_+ + e^{i\theta/2}P_- = \tfrac12\big(e^{-i\theta/2} + e^{i\theta/2}\big)I + \tfrac12\big(e^{-i\theta/2} - e^{i\theta/2}\big)H = \cos\tfrac\theta2\,I - i\sin\tfrac\theta2\,H.\ ∎$$
(The project checks this against `scipy.linalg.expm`.)

(d) At $\theta = 2\pi$: $\cos\pi\,I - i\sin\pi\,H = -I$. A full $2\pi$ rotation of the Bloch sphere returns every direction $\mathbf{n}$ to itself but multiplies every state by $-1$. For a single qubit this global phase is unobservable (QC-1); the sign is the fingerprint of $SU(2)$ double-covering $SO(3)$ (LA-12), and it is physically visible in interference between a rotated and an unrotated spin (neutron interferometry).

**B7.** (a) $\|Q^\top AQ\|_F^2 = \operatorname{tr}\big((Q^\top AQ)^\top(Q^\top AQ)\big) = \operatorname{tr}(Q^\top A^\top QQ^\top AQ) = \operatorname{tr}(Q^\top A^\top AQ) = \operatorname{tr}(A^\top AQQ^\top) = \operatorname{tr}(A^\top A) = \|A\|_F^2$, using $QQ^\top = I$ twice and the cyclic property $\operatorname{tr}(XY) = \operatorname{tr}(YX)$. ∎

(b) $G$ equals $I$ outside rows/columns $p, q$. For $i, j\notin\{p, q\}$: $(G^\top AG)_{ij} = \sum_{k,l}G_{ki}a_{kl}G_{lj} = a_{ij}$ because $G_{ki} = \delta_{ki}$ and $G_{lj} = \delta_{lj}$. So only rows and columns $p, q$ change.
For $i\notin\{p,q\}$: $(AG)_{ip} = ca_{ip} + sa_{iq}$ and $(AG)_{iq} = -sa_{ip} + ca_{iq}$ (columns $p,q$ of $G$ are $(c, s)$ and $(-s, c)$ in rows $p, q$), and left-multiplying by $G^\top$ leaves row $i$ alone. So $(a'_{ip}, a'_{iq})$ is the rotation by $-\theta$ of $(a_{ip}, a_{iq})$: $a'^2_{ip} + a'^2_{iq} = a_{ip}^2 + a_{iq}^2$. By symmetry the same holds for $(a'_{pi}, a'_{qi})$.
The $2\times2$ block: $\begin{pmatrix}a'_{pp}&a'_{pq}\\a'_{pq}&a'_{qq}\end{pmatrix} = R^\top\begin{pmatrix}a_{pp}&a_{pq}\\a_{pq}&a_{qq}\end{pmatrix}R$ with $R = \begin{pmatrix}c&-s\\s&c\end{pmatrix}$, an orthogonal similarity; by (a) its Frobenius norm is preserved: $a'^2_{pp} + 2a'^2_{pq} + a'^2_{qq} = a_{pp}^2 + 2a_{pq}^2 + a_{qq}^2$. With $a'_{pq} = 0$: $a'^2_{pp} + a'^2_{qq} = a_{pp}^2 + a_{qq}^2 + 2a_{pq}^2$.
Now $\|A'\|_F^2 = \|A\|_F^2$ by (a), and the sum of squared diagonal entries went up by exactly $2a_{pq}^2$ (only $a_{pp}, a_{qq}$ changed on the diagonal). Since $\operatorname{off}^2 = \|\cdot\|_F^2 - (\text{diagonal squares})$, $\operatorname{off}(A')^2 = \operatorname{off}(A)^2 - 2a_{pq}^2$. ∎

(c) There are $N = n(n-1)$ off-diagonal entries, and the largest of their squares is at least their average: $a_{pq}^2\ge\operatorname{off}(A)^2/N$. So $\operatorname{off}(A')^2\le\operatorname{off}(A)^2(1 - 2/N)$, and after $k$ rotations $\operatorname{off}(A_k)^2\le(1 - 2/N)^k\operatorname{off}(A_0)^2\to0$.

*The diagonal converges to the eigenvalues.* $\operatorname{off}(A_k)\to0$ only says $A_k$ approaches its own diagonal part $D_k = \operatorname{diag}(a^{(k)}_{11}, \dots, a^{(k)}_{nn})$. To conclude that the diagonal entries approach the eigenvalues you need a perturbation bound, and Weyl's inequality (the corollary to Theorem 3.18) is exactly that. Apply it with $X = A_k$ and $Y = D_k$. The eigenvalues of $A_k$ are those of $A$ (orthogonal similarity), the eigenvalues of $D_k$ are its diagonal entries, and $\|X - Y\|_F = \operatorname{off}(A_k)$. So if $d^{(k)}_{(1)}\ge\cdots\ge d^{(k)}_{(n)}$ is the diagonal of $A_k$ sorted in decreasing order,
$$|d^{(k)}_{(i)} - \lambda_i|\le\operatorname{off}(A_k)\to0\qquad\text{for every } i.$$
The sorted diagonal converges to the sorted eigenvalues, which is all the algorithm needs, since it sorts the diagonal at the end (`jacobi_eigh` does).

*Each entry converges too.* A rotation changes only $a_{pp}$ and $a_{qq}$. Column $p$ of $G$ is $c\,\mathbf{e}_p + s\,\mathbf{e}_q$, so $a'_{pp} = c^2a_{pp} + 2cs\,a_{pq} + s^2a_{qq}$. The zeroing condition from (b) is $(c^2 - s^2)a_{pq} = cs(a_{pp} - a_{qq})$; with $t = s/c = \tan\theta$,
$$a'_{pp} - a_{pp} = -s^2(a_{pp} - a_{qq}) + 2cs\,a_{pq} = -t(c^2 - s^2)a_{pq} + 2tc^2a_{pq} = t(c^2 + s^2)a_{pq} = t\,a_{pq},$$
and $a'_{qq} = a_{qq} - t\,a_{pq}$ because the trace of the $2\times2$ block is unchanged. The project picks $|\theta|\le\pi/4$, so $|t|\le1$ and each diagonal entry moves by at most $|a^{(k)}_{pq}|\le\operatorname{off}(A_k)\le(1 - 2/N)^{k/2}\operatorname{off}(A_0)$ at step $k$. These bounds have a finite (geometric) sum, so each diagonal entry is a Cauchy sequence and converges, and by the Weyl bound its limit is one of the $\lambda_i$. The project's "Your turn" 5 asserts both the Weyl bound and $|a'_{pp} - a_{pp}|\le|a_{pq}|$ after every rotation.

(In practice convergence becomes quadratic near the end: once the diagonal entries have separated, every $n(n-1)/2$ rotations roughly square $\operatorname{off}(A)$. The project's semilog plot shows the slope steepening over the last few dozen rotations.)

(d) Reference implementation (the same as "Your turn" 4 in `project.py`):
```python
def cyclic_jacobi_eigh(A, tol=1e-12, max_sweeps=50):
    A = np.array(A, float); n = A.shape[0]; Q = np.eye(n); rotations = 0
    for sweep in range(1, max_sweeps + 1):
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(A[p, q]) < tol:
                    continue
                apq = A[p, q]
                A_new, G = jacobi_rotate(A, p, q)
                assert_off_drop(A, A_new, apq)   # off(A')^2 = off(A)^2 - 2 a_pq^2, to rounding
                A, Q = A_new, Q @ G; rotations += 1
        if np.abs(A - np.diag(np.diag(A))).max() < tol:
            break
    w = np.diag(A); order = np.argsort(w)[::-1]
    return w[order], Q[:, order], rotations, sweep
```
`assert_off_drop` is the helper that `jacobi_eigh` already uses; its tolerance is $10^{-14}\,\|A\|_F\operatorname{off}(A)$, because rounding in the rotated rows and columns perturbs $\operatorname{off}(A)^2$ by about $\varepsilon\|A\|_F\operatorname{off}(A)$. On the project's `M8` the classical method needs **91** rotations and the cyclic method **124**, in 5 sweeps of $28$ pairs; both give the eigenvalues of `np.linalg.eigh` and an orthogonal $Q$ with $Q\Lambda Q^\top = $ `M8`.

*Why more rotations.* The cyclic method also rotates on small entries, and a rotation lowers $\operatorname{off}^2$ by only $2a_{pq}^2$, which is little when $a_{pq}$ is small.

*Why the identity still holds.* The proof of (b) never used that $a_{pq}$ was the largest entry: it holds for a rotation in any plane $(p, q)$. What does break is the convergence argument in (c), which needed $a_{pq}^2\ge\operatorname{off}(A)^2/N$. Convergence of the cyclic method is still true, but needs a different proof (Forsythe and Henrici, 1960).

*Why people use the cyclic version anyway.* Finding the largest of the $n(n-1)/2$ entries above the diagonal costs $O(n^2)$ comparisons per rotation, while the rotation itself changes only rows and columns $p$ and $q$, which is $O(n)$ work when coded that way. (Our `jacobi_rotate` forms $G^\top AG$ with full matrix products for clarity, at $O(n^3)$ cost, so timing this code would not show the difference.) For large $n$ the search dominates, and $124$ cheap rotations beat $91$ expensive ones.

**B8.** (a) $B^\top = (V^\top AV)^\top = V^\top A^\top V = V^\top AV = B$. ∎

(b) Let $\mathbf{x} = V\mathbf{z}$. *Orthogonal to $\mathbf{q}_1$:* $\mathbf{q}_1^\top\mathbf{x} = (\mathbf{q}_1^\top V)\mathbf{z} = \mathbf{0}^\top\mathbf{z} = 0$ since every column of $V$ is orthogonal to $\mathbf{q}_1$. *Norm:* $\|\mathbf{x}\|^2 = \mathbf{z}^\top V^\top V\mathbf{z} = \mathbf{z}^\top\mathbf{z}$. *Eigenvector:* first, $A\mathbf{x}\in\mathbf{q}_1^\perp$ because $\mathbf{q}_1^\top A\mathbf{x} = (A\mathbf{q}_1)^\top\mathbf{x} = \lambda_1\mathbf{q}_1^\top\mathbf{x} = 0$ (symmetry, then eigenvector). The columns of $V$ are an orthonormal basis of $\mathbf{q}_1^\perp$, so the orthogonal projection onto $\mathbf{q}_1^\perp$ is $VV^\top$ (LA-5), and it fixes any vector already in $\mathbf{q}_1^\perp$. Hence
$$A\mathbf{x} = VV^\top A\mathbf{x} = V(V^\top AV)\mathbf{z} = VB\mathbf{z} = \mu V\mathbf{z} = \mu\mathbf{x}.\ ∎$$

(c) $A\mathbf{e}_2 = (1, 1)^\top$, and $\mathbf{e}_1^\top(A\mathbf{e}_2) = 1\ne0$, so $A\mathbf{e}_2\notin\mathbf{e}_1^\perp = \operatorname{span}\{\mathbf{e}_2\}$: not invariant. Eigenvalues $2, 1$ (triangular). $\lambda = 2$: $(A - 2I) = \begin{pmatrix}0&1\\0&-1\end{pmatrix}$, eigenvector $\mathbf{e}_1$. $\lambda = 1$: $(A - I) = \begin{pmatrix}1&1\\0&0\end{pmatrix}$, eigenvector $(1, -1)^\top$. Their dot product is $1\ne0$. The line that breaks is step 7, $\mathbf{q}_1^\top(A\mathbf{x}) = (A^\top\mathbf{q}_1)^\top\mathbf{x} = (A\mathbf{q}_1)^\top\mathbf{x}$: the second equality needs $A^\top = A$, and here $A^\top\mathbf{e}_1 = (2, 1)^\top\ne A\mathbf{e}_1 = (2, 0)^\top$. Equivalently, in step 8 with $Q_1 = I$ the block matrix is $A$ itself, whose top-right entry is $1$, not $0$.

## C. Challenge

**C1.** Let $V$ be the first $n-1$ columns of $I_n$, so $B = V^\top AV$, $V^\top V = I_{n-1}$, and the range of $V$ is the hyperplane $W = \{\mathbf{x} : x_n = 0\}$. For $\mathbf{z}\in\mathbb{R}^{n-1}$, $\mathbf{z}^\top B\mathbf{z}/\mathbf{z}^\top\mathbf{z} = (V\mathbf{z})^\top A(V\mathbf{z})/\|V\mathbf{z}\|^2 = R_A(V\mathbf{z})$. The map $S'\mapsto VS'$ is a bijection from $k$-dimensional subspaces of $\mathbb{R}^{n-1}$ to $k$-dimensional subspaces of $W$. So Courant–Fischer for $B$ reads
$$\mu_k = \max_{\substack{S\subseteq W\\\dim S = k}}\ \min_{\mathbf{x}\in S\setminus0}R_A(\mathbf{x}).$$
*Upper bound.* The maximum over $k$-dimensional subspaces contained in $W$ is at most the maximum over all $k$-dimensional subspaces of $\mathbb{R}^n$, which is $\lambda_k$ by (3.14). So $\mu_k\le\lambda_k$.
*Lower bound.* Let $S$ be any $(k+1)$-dimensional subspace of $\mathbb{R}^n$. By LA-3's dimension formula (Theorem 3.35), $\dim(S\cap W)\ge(k+1) + (n-1) - n = k$. Choose $T\subseteq S\cap W$ with $\dim T = k$. Shrinking the set can only raise a minimum: $\min_{S}R_A\le\min_TR_A$. And $T$ is a $k$-dimensional subspace of $W$, so $\min_TR_A\le\mu_k$. Hence $\min_SR_A\le\mu_k$ for every such $S$; taking the maximum over $S$ gives $\lambda_{k+1}\le\mu_k$ by (3.14). ∎
*Check.* Example 4.2: $A$ has $\lambda = 4, 1, 1$; the leading block $B = \begin{pmatrix}2&1\\1&2\end{pmatrix}$ has $\mu = 3, 1$. Interlacing: $4\ge3\ge1\ge1\ge1$. ✓

**C2.** (a) $\|u+v\|^2 = \|u\|^2 + 2\langle u, v\rangle + \|v\|^2$ and $\|u-v\|^2 = \|u\|^2 - 2\langle u, v\rangle + \|v\|^2$. Subtract: $\|u+v\|^2 - \|u-v\|^2 = 4\langle u, v\rangle$. ∎

(b) For any scalar $\alpha$, using conjugate-linearity in the first slot and linearity in the second,
$$\|u + \alpha v\|^2 = \|u\|^2 + \langle u, \alpha v\rangle + \langle\alpha v, u\rangle + |\alpha|^2\|v\|^2 = \|u\|^2 + \alpha\langle u, v\rangle + \bar\alpha\overline{\langle u, v\rangle} + |\alpha|^2\|v\|^2 = \|u\|^2 + 2\operatorname{Re}\big(\alpha\langle u, v\rangle\big) + |\alpha|^2\|v\|^2.$$
$\alpha = 1$ minus $\alpha = -1$: $\|u+v\|^2 - \|u-v\|^2 = 4\operatorname{Re}\langle u, v\rangle$.
$\alpha = -i$ minus $\alpha = i$: $\|u - iv\|^2 - \|u + iv\|^2 = 2\operatorname{Re}(-i\langle u,v\rangle) - 2\operatorname{Re}(i\langle u,v\rangle) = -4\operatorname{Re}(i\langle u, v\rangle) = 4\operatorname{Im}\langle u, v\rangle$, since $\operatorname{Re}(iz) = -\operatorname{Im}z$.
So $\langle u, v\rangle = \operatorname{Re} + i\operatorname{Im} = \tfrac14(\|u+v\|^2 - \|u-v\|^2) + \tfrac i4(\|u-iv\|^2 - \|u+iv\|^2)$. ∎
In Axler's convention $\langle u, \alpha v\rangle = \bar\alpha\langle u, v\rangle$, so the cross term is $2\operatorname{Re}(\bar\alpha\langle u, v\rangle)$ and the imaginary part comes out with the opposite sign: $+\tfrac i4(\|u+iv\|^2 - \|u-iv\|^2)$. Same real part.

(c) Let $U$ be linear with $\|U\mathbf{x}\| = \|\mathbf{x}\|$ for all $\mathbf{x}$. By (b) applied to $U\mathbf{u}, U\mathbf{v}$ and linearity of $U$:
$$\langle U\mathbf{u}, U\mathbf{v}\rangle = \tfrac14\big(\|U(\mathbf{u}+\mathbf{v})\|^2 - \|U(\mathbf{u}-\mathbf{v})\|^2\big) + \tfrac i4\big(\|U(\mathbf{u}-i\mathbf{v})\|^2 - \|U(\mathbf{u}+i\mathbf{v})\|^2\big).$$
Each norm on the right equals the norm of the argument without $U$, so the right side is the polarization formula for $\langle\mathbf{u},\mathbf{v}\rangle$. Hence $\langle U\mathbf{u}, U\mathbf{v}\rangle = \langle\mathbf{u},\mathbf{v}\rangle$, which is (b) of Theorem 3.20, hence $U^\dagger U = I$. ∎

**C3.** (a) $M$ is SPD, so by Theorem 3.15 (iii) $M = B^\top B$ with $B$ invertible; set $L = B^\top$, so $M = LL^\top$ (Cholesky is the triangular choice). Substituting $\mathbf{v} = L^{-\top}\mathbf{w}$ into $K\mathbf{v} = \lambda M\mathbf{v} = \lambda LL^\top\mathbf{v}$: $KL^{-\top}\mathbf{w} = \lambda L\mathbf{w}$. Multiply by $L^{-1}$: $C\mathbf{w} = \lambda\mathbf{w}$ with $C = L^{-1}KL^{-\top}$. $C^\top = (L^{-\top})^\top K^\top L^{-\top} = L^{-1}KL^{-\top} = C$. ∎

(b) $C$ is real symmetric, so by Theorem 3.11 its eigenvalues are real and it has orthonormal eigenvectors $\mathbf{w}_i$. The $\mathbf{v}_i = L^{-\top}\mathbf{w}_i$ solve the generalized problem with the same $\lambda_i$, and
$$\mathbf{v}_i^\top M\mathbf{v}_j = \mathbf{w}_i^\top L^{-1}(LL^\top)L^{-\top}\mathbf{w}_j = \mathbf{w}_i^\top\mathbf{w}_j = \delta_{ij}.\ ∎$$

(c) $\det(K - \lambda M) = \det\begin{pmatrix}2-\lambda&-1\\-1&2-2\lambda\end{pmatrix} = (2-\lambda)(2-2\lambda) - 1 = 2\lambda^2 - 6\lambda + 3$. Roots $\lambda = \dfrac{6\pm\sqrt{36 - 24}}{4} = \dfrac{3\pm\sqrt3}{2}$: $\lambda_1\approx2.366$, $\lambda_2\approx0.634$. First row of $(K - \lambda M)\mathbf{v} = 0$: $(2-\lambda)v_1 - v_2 = 0$, so $\mathbf{v} = (1, 2-\lambda)^\top$.
$\mathbf{v}_1 = \big(1, \tfrac{1-\sqrt3}2\big)^\top\approx(1, -0.366)^\top$ (masses move oppositely: the fast mode).
$\mathbf{v}_2 = \big(1, \tfrac{1+\sqrt3}2\big)^\top\approx(1, 1.366)^\top$ (masses move together: the slow mode).
$M$-orthogonality: $\mathbf{v}_1^\top M\mathbf{v}_2 = 1\cdot1\cdot1 + 2\cdot\tfrac{1-\sqrt3}2\cdot\tfrac{1+\sqrt3}2 = 1 + 2\cdot\tfrac{1-3}4 = 1 - 1 = 0$. ✓ Note the ordinary dot product is $1 + \tfrac{1-3}4 = \tfrac12\ne0$: the eigenvectors are orthogonal in the $M$-geometry, not the Euclidean one.

(d) `np.linalg.eigh` accepts one matrix only. `scipy.linalg.eigh(K, M)` solves $K\mathbf{v} = \lambda M\mathbf{v}$ for symmetric $K$ and SPD $M$ by exactly the reduction in (a) (Cholesky of $M$, then a symmetric eigensolver), returning eigenvalues in ascending order and eigenvectors normalized so that $V^\top MV = I$. On this example it returns $\lambda\approx(0.634, 2.366)$ and columns proportional to $\mathbf{v}_2, \mathbf{v}_1$.

**C4.** (a) The characteristic polynomial $\det(A - \lambda I)$ has a root in $\mathbb{C}$ (fundamental theorem of algebra), hence an eigenvector; normalize it to $\mathbf{v}$. The columns $\mathbf{v}, \mathbf{w}$ are orthonormal so $U$ is unitary. First column of $T = U^\dagger AU$: $U^\dagger A\mathbf{v} = \lambda U^\dagger\mathbf{v} = \lambda(\mathbf{v}^\dagger\mathbf{v}, \mathbf{w}^\dagger\mathbf{v})^\top = (\lambda, 0)^\top$. So $T = \begin{pmatrix}\lambda&b\\0&\mu\end{pmatrix}$ for some $b, \mu$.

(b) $T^\dagger T = U^\dagger A^\dagger UU^\dagger AU = U^\dagger A^\dagger AU = U^\dagger AA^\dagger U = U^\dagger AUU^\dagger A^\dagger U = TT^\dagger$, using $UU^\dagger = I$ and $A^\dagger A = AA^\dagger$. ∎

(c) $(TT^\dagger)_{11} = |\lambda|^2 + |b|^2$ (row 1 of $T$ dotted with itself) and $(T^\dagger T)_{11} = |\lambda|^2$ (column 1 of $T$ dotted with itself). Normality says they are equal, so $|b|^2 = 0$, $b = 0$.

(d) $T = \operatorname{diag}(\lambda, \mu)$ and $A = UTU^\dagger$: unitarily diagonalizable. ∎ In general, Schur's theorem (LA-10) gives $A = UTU^\dagger$ with $T$ upper triangular for *any* square $A$; (b) shows $T$ is normal when $A$ is. Comparing the $(1,1)$ entries of $TT^\dagger$ and $T^\dagger T$ kills the whole first row off the diagonal ($\sum_j|t_{1j}|^2 = |t_{11}|^2$); then the $(2,2)$ entries kill the second row, and so on down: $T$ is diagonal.
