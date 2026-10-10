# LA-6 Exercises

Conventions as in the notes: $\det$ is defined on rows by (D1)–(D3); $M_{ij}$ is the minor, $C_{ij} = (-1)^{i+j}M_{ij}$ the cofactor, $\operatorname{adj}A = C^\top$; $\operatorname{sgn}(\sigma) = (-1)^{N(\sigma)}$ with $N$ the inversion count.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Let $A = \begin{pmatrix} 2 & 1 & 0 \\ 4 & 3 & 1 \\ -2 & 0 & 5 \end{pmatrix}$. Compute $\det A$ three ways: by elimination (3.2), recording any swaps; by cofactor expansion along row 1 (3.4); and by cofactor expansion down column 3. Confirm all three agree.

**A2.** A $3\times 3$ matrix $A$ has $\det A = 5$. Without knowing $A$, find the determinant of: (a) $A$ with rows 1 and 3 swapped; (b) $A$ with row 2 multiplied by $3$; (c) $A$ with $4\times$(row 1) added to row 2; (d) $2A$; (e) $A^\top$; (f) $A^{-1}$; (g) $A^2$; (h) $-A$; (i) $A$ with rows 1 and 2 swapped and then column 1 multiplied by $-2$.

**A3.** Compute $N(\sigma)$ and $\operatorname{sgn}(\sigma)$ for $\sigma = (3, 1, 2)$, $(2, 1, 4, 3)$, and $(4, 3, 1, 2)$. Write the term of the $4\times 4$ Leibniz formula (3.3) that corresponds to $(4, 3, 1, 2)$, with its sign. Then go the other way: which permutation produces the term $a_{13}a_{21}a_{34}a_{42}$, and what is its sign? How many of the $24$ permutations of $\{1,2,3,4\}$ have sign $+1$?

**A4.** Solve $\begin{pmatrix} 3 & 1 \\ 4 & 2 \end{pmatrix}\mathbf{x} = \begin{pmatrix} 5 \\ 8 \end{pmatrix}$ by Cramer's rule (3.5). Then compute the cofactor matrix $C$ and the adjugate, verify $A\operatorname{adj}A = (\det A)I$ by direct multiplication, and recover $\mathbf{x} = A^{-1}\mathbf{b}$ from (3.6).

**A5.** Let $\mathbf{u} = (1, 2, 0)^\top$ and $\mathbf{v} = (0, 1, 3)^\top$. Compute $\mathbf{u}\times\mathbf{v}$ by the cofactor expansion (3.7), check it is orthogonal to both $\mathbf{u}$ and $\mathbf{v}$, find the area of the parallelogram they span, and verify Theorem 3.15(c) numerically. Then compute the volume of the parallelepiped on $\mathbf{u}, \mathbf{v}, \mathbf{w} = (1,1,1)^\top$ two ways: as $|\mathbf{w}\cdot(\mathbf{u}\times\mathbf{v})|$ and as a $3\times 3$ determinant.

## B. Core (do all)

**B1.** (Block triangular.) Let $A$ be $k\times k$, $D$ be $m\times m$, $B$ be $k\times m$, and $M = \begin{pmatrix} A & B \\ 0 & D \end{pmatrix}$.
(a) Prove $\det M = \det A\cdot\det D$. Suggested route: with $B, D$ fixed, show that $A \mapsto \det M$ satisfies (D2) and (D3) as a function of the rows of $A$, so by the uniqueness argument of Theorem 3.9 it equals $\det A$ times its value at $A = I_k$; then evaluate that value.
(b) Deduce the same formula for $\begin{pmatrix} A & 0 \\ C & D\end{pmatrix}$.
(c) Show by a $4\times 4$ example with $2\times 2$ blocks that $\det\begin{pmatrix} A & B \\ C & D\end{pmatrix} = \det A\det D - \det B\det C$ is **false** in general.

**B2.** (Determinants of maps.)
(a) Explain, using Corollary 3.10′(c), why it makes sense to speak of the determinant of a linear map $T : \mathbb{R}^n \to \mathbb{R}^n$ without choosing a basis.
(b) Prove that an orthogonal matrix ($Q^\top Q = I$) has $\det Q = \pm 1$, and that a $2\times 2$ rotation matrix has $\det = +1$.
(c) Prove that a projection matrix ($P^2 = P$, $P \neq I$) has $\det P = 0$.
(d) The reflection across the plane orthogonal to a unit vector $\mathbf{u} \in \mathbb{R}^3$ is $H = I - 2\mathbf{u}\mathbf{u}^\top$. Prove $\det H = -1$ by choosing a convenient basis and using (a).

**B3.** (Recognizing alternating multilinear forms.)
(a) Show that $f(\mathbf{u}, \mathbf{v}, \mathbf{w}) = \mathbf{u}\cdot(\mathbf{v}\times\mathbf{w})$ is multilinear and alternating on $\mathbb{R}^3$. Use Theorem 3.16 to conclude $f(\mathbf{u}, \mathbf{v}, \mathbf{w}) = \det\begin{pmatrix}\mathbf{u} & \mathbf{v} & \mathbf{w}\end{pmatrix}$ **without expanding anything**, and deduce the cyclic identity $\mathbf{u}\cdot(\mathbf{v}\times\mathbf{w}) = \mathbf{v}\cdot(\mathbf{w}\times\mathbf{u}) = \mathbf{w}\cdot(\mathbf{u}\times\mathbf{v})$.
(b) Show that $(\mathbf{u}, \mathbf{v}) \mapsto \mathbf{u}\cdot\mathbf{v}$ is multilinear but not alternating, and that $(\mathbf{u}, \mathbf{v}) \mapsto \|\mathbf{u}\times\mathbf{v}\|$ vanishes whenever $\mathbf{u} = \mathbf{v}$ but is not multilinear.
(c) A function $f$ on pairs of vectors in $\mathbb{R}^2$ is bilinear and alternating, and $f\big((1,1)^\top, (1,-1)^\top\big) = 4$. Find $f\big((3,0)^\top, (0,2)^\top\big)$.

**B4.** (Vandermonde consequences.)
(a) For the nodes $x_i = i - 1$, $i = 1, \dots, n$ (that is $0, 1, \dots, n-1$), show that $\det V = \prod_{j=1}^{n-1} j!$. Evaluate for $n = 4$ and check it against (3.11) term by term.
(b) Use $\det V \neq 0$ to prove: a polynomial of degree $\le n-1$ that vanishes at $n$ distinct points is the zero polynomial.
(c) Conclude that the polynomial of degree $\le n-1$ through $n$ points with distinct $x$-coordinates is unique.
(d) Suppose $x_2 \to x_1$ with the other nodes fixed. Using Cramer's rule (3.5), describe what happens to the coefficients of the interpolating polynomial when $y_1 \neq y_2$. Why is this the right behavior?

**B5.** (Tetrahedra.)
(a) Compute the volume of the standard simplex $\Delta = \{(t_1, t_2, t_3) : t_i \ge 0,\ t_1 + t_2 + t_3 \le 1\}$ by a triple integral.
(b) Show that the tetrahedron with vertices $\mathbf{0}, \mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3$ is the image of $\Delta$ under the linear map with matrix $(\mathbf{v}_1\ \mathbf{v}_2\ \mathbf{v}_3)$, and conclude from Theorem 3.17 that its volume is $\tfrac16|\det(\mathbf{v}_1\ \mathbf{v}_2\ \mathbf{v}_3)|$.
(c) Find the volume of the tetrahedron with vertices $P_0 = (0,0,0)$, $P_1 = (2,0,0)$, $P_2 = (0,3,0)$, $P_3 = (1,1,4)$. Check it against $\tfrac13(\text{base area})(\text{height})$.
(d) Give the formula for general vertices $P_0, P_1, P_2, P_3$, and show it can also be written as $\tfrac16\left|\det\begin{pmatrix} 1 & 1 & 1 & 1 \\ P_0 & P_1 & P_2 & P_3\end{pmatrix}\right|$ (a $4\times 4$ determinant with the coordinates of $P_i$ as the lower three entries of column $i$).
(e) Explain why rotating or translating all four vertices leaves the volume unchanged, using properties from §3.4.

**B6.** (Spherical and cylindrical Jacobians.) Spherical coordinates (physics convention): $x = r\sin\phi\cos\theta$, $y = r\sin\phi\sin\theta$, $z = r\cos\phi$, with $\phi \in [0, \pi]$ the polar angle and $\theta \in [0, 2\pi)$ the azimuth.
(a) Write the $3\times 3$ Jacobian matrix of $(r, \phi, \theta) \mapsto (x, y, z)$ and compute its determinant by cofactor expansion along the third row. Simplify to $r^2\sin\phi$.
(b) Show that the three columns are mutually orthogonal, with lengths $1$, $r$, $r\sin\phi$, and explain why that alone gives $|\det J| = r^2\sin\phi$.
(c) What is the sign of $\det J$ for $0 < \phi < \pi$? What happens to the sign if you order the variables $(r, \theta, \phi)$ instead, and why doesn't it matter for integration?
(d) Compute the Jacobian determinant of cylindrical coordinates $(r, \theta, z) \mapsto (r\cos\theta, r\sin\theta, z)$ using B1.

**B7.** (Adjugate facts.) Let $A$ be $n\times n$ with $n \ge 2$.
(a) Prove $\det(\operatorname{adj}A) = (\det A)^{n-1}$ when $A$ is invertible.
(b) Prove the same identity when $A$ is singular (both sides are $0$). Hint: if $\operatorname{adj}A$ were invertible, what would (3.6) say about $A$?
(c) Prove that $\operatorname{adj}(A)\,A = (\det A)I$ as well (the notes proved $A\operatorname{adj}A$). Hint: apply (3.6) to $A^\top$ and show $\operatorname{adj}(A^\top) = (\operatorname{adj}A)^\top$.
(d) Prove: if $A$ has integer entries, then $A^{-1}$ has integer entries if and only if $\det A = \pm 1$.

**B8.** (Code.) In `project.py`, implement `adjugate(A)` using a `cofactor(A, i, j)` helper, and assert `A @ adjugate(A)` equals `det(A) * I` on random $4\times 4$ matrices. Then take the $4\times 4$ Hilbert matrix $H_{ij} = 1/(i + j - 1)$, compute its inverse as `adjugate(H) / det(H)` and as `np.linalg.inv(H)`, and compare both against the exact integer inverse (which SymPy gives you: `sympy.Matrix(4, 4, lambda i, j: sympy.Rational(1, i + j + 1)).inv()`). Report the maximum entrywise error of each, the condition number `np.linalg.cond(H)`, and write a short comment: does the adjugate formula buy accuracy? What does it cost for an $n\times n$ matrix?

## C. Challenge (do at least one)

**C1.** (Jacobi's formula and $\det e^{B} = e^{\operatorname{tr}B}$.)
(a) Using the Leibniz formula, show that $\det(I + tB) = 1 + t\operatorname{tr}B + O(t^2)$ for any $n\times n$ matrix $B$, where $\operatorname{tr}B = \sum_i b_{ii}$. Hence $\frac{d}{dt}\det(I + tB)\big|_{t=0} = \operatorname{tr}B$.
(b) Let $A(t)$ be a differentiable matrix function with $A(t)$ invertible. Prove $\frac{d}{dt}\det A(t) = \det A(t)\cdot\operatorname{tr}\big(A(t)^{-1}A'(t)\big)$.
(c) Prove $\operatorname{tr}(XY) = \operatorname{tr}(YX)$ for square $X, Y$.
(d) Let $A(t)$ solve the matrix ODE $A' = BA$, $A(0) = I$ (the solution is the matrix exponential $e^{tB}$, which LA-7 constructs; assume existence and uniqueness, as for any linear ODE). Show $\det A(t) = e^{t\operatorname{tr}B}$. Conclude $\det e^{B} = e^{\operatorname{tr}B}$, and that $A(t)$ is invertible for every $t$.
(e) Where this goes: in CM-8, Liouville's theorem says a Hamiltonian flow preserves phase-space volume. Explain in two sentences how (d) with $\operatorname{tr}B = 0$ is the linear version of that statement.

**C2.** (The characteristic polynomial: a preview of LA-7.) Let $p(\lambda) = \det(A - \lambda I)$.
(a) Using the Leibniz formula, show $p$ is a polynomial in $\lambda$ of degree exactly $n$ with leading coefficient $(-1)^n$, that the coefficient of $\lambda^{n-1}$ is $(-1)^{n-1}\operatorname{tr}A$, and that the constant term is $\det A$.
(b) If $p$ factors as $(-1)^n\prod_{i=1}^n(\lambda - \lambda_i)$ (the fundamental theorem of algebra, F1, guarantees this over $\mathbb{C}$), show $\sum_i\lambda_i = \operatorname{tr}A$ and $\prod_i\lambda_i = \det A$.
(c) Write $p$ for a general $2\times 2$ matrix, then find the roots for $A = \begin{pmatrix} 2 & 1 \\ 1 & 2\end{pmatrix}$ and check (b).
(d) Prove that $p$ is unchanged if $A$ is replaced by $S^{-1}AS$. Why does this matter for the definition of eigenvalues of a linear map?

**C3.** (Hadamard's inequality and the Gram determinant; needs the QR factorization from LA-5.)
(a) Prove $|\det A| \le \prod_{j=1}^n\|\mathbf{a}_j\|$ (product of the column lengths), with equality iff the columns are mutually orthogonal (or some column is zero). Route: $A = QR$ with $Q$ orthogonal and $R$ upper triangular; relate $|\det A|$ to the diagonal of $R$ and $\|\mathbf{a}_j\|$ to column $j$ of $R$.
(b) Interpret (a) as a statement about the volume of a parallelepiped versus the volume of a box with the same edge lengths.
(c) For an $n\times k$ matrix $A$ with $k \le n$, the Gram matrix is $G = A^\top A$ ($k\times k$). Show $\det G = (\det A)^2$ when $k = n$, and that for $k = 2$, $n = 3$, $\det G = \|\mathbf{a}_1\times\mathbf{a}_2\|^2$, the squared area of the parallelogram. (So $\sqrt{\det G}$ is the $k$-dimensional volume of the parallelepiped spanned by the columns, in any $\mathbb{R}^n$; this is how LA-12 and EM-1 measure area on surfaces.)
(d) (Code.) Add a cell to `project.py` that checks (a) on $1000$ random $5\times 5$ matrices, checks equality for a random matrix with orthogonal columns, and checks (c) for a random $3\times 2$ matrix against the cross product.
