# LA-8 Exercises

Conventions as in the notes: complex inner products are conjugate-linear in the **first** slot, $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^\dagger\mathbf{v}$; $A^\dagger$ is the conjugate transpose; eigenvalues of symmetric matrices are ordered $\lambda_1\ge\cdots\ge\lambda_n$; $\langle A, B\rangle_F = \operatorname{tr}(A^\dagger B)$.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Let $A = \begin{pmatrix}1&2\\0&-1\end{pmatrix}$ and $B = \begin{pmatrix}3&0\\1&1\end{pmatrix}$.
(a) Compute $\langle A, B\rangle_F$ twice: as $\operatorname{tr}(A^\top B)$ and as $\sum_{ij}a_{ij}b_{ij}$.
(b) Compute $\|A\|_F$, $\|B\|_F$, and the cosine of the angle between $A$ and $B$.
(c) Verify that (3.3), $\langle A, B\rangle_F = \operatorname{tr}(A^\dagger B)$ on complex $m\times n$ matrices, satisfies the three axioms of Definition 3.2.

**A2.** Use $\langle f, g\rangle = \int_0^1 f(x)g(x)\,dx$.
(a) Compute $\langle x, x^2\rangle$, $\|x\|$, $\|x^2\|$, and the angle between the functions $x$ and $x^2$.
(b) Check Cauchy–Schwarz for this pair.
(c) Find the constant $c$ such that $x - c$ is orthogonal to the constant function $1$.

**A3.** Let $A = \begin{pmatrix}2&1\\1&2\end{pmatrix}$.
(a) Find $Q$ and $\Lambda$ with $A = Q\Lambda Q^\top$, $Q$ orthogonal.
(b) Write the spectral decomposition $A = \lambda_1P_1 + \lambda_2P_2$ with explicit $2\times2$ projection matrices $P_i$, and check $P_1 + P_2 = I$, $P_1P_2 = 0$.
(c) Using (3.11), compute $\sqrt A$, $A^{-1}$, and $e^A$. Verify $(\sqrt A)^2 = A$ and $A^{-1}A = I$ by direct multiplication.

**A4.** Classify each quadratic form as positive definite, positive semidefinite, negative definite, or indefinite, and describe the level set $q(\mathbf{x}) = 1$ (ellipse, hyperbola, pair of lines, empty), with its axes.
(a) $q = 2x^2 - 2xy + 2y^2$. (b) $q = x^2 + 4xy + y^2$. (c) $q = x^2 + 2xy + y^2$. (d) $q = -x^2 + 4xy - 5y^2$.

**A5.** (Symmetric is not Hermitian.) Let $S = \begin{pmatrix}0&i\\i&0\end{pmatrix}$, $N = \begin{pmatrix}1&i\\i&-1\end{pmatrix}$, and $\sigma_y = \begin{pmatrix}0&-i\\i&0\end{pmatrix}$.
(a) Which of the three are symmetric ($A^\top = A$)? Which are Hermitian ($A^\dagger = A$)?
(b) Find the eigenvalues of $S$. Are they real?
(c) Compute $N^2$. What are the eigenvalues of $N$? How many independent eigenvectors does $N$ have? Which theorem of §3 does this show fails for complex symmetric matrices?
(d) Confirm the eigenvalues of $\sigma_y$ are real, even though $\sigma_y$ is not symmetric.

**A6.** Let $A = \begin{pmatrix}3&1\\1&3\end{pmatrix}$. Compute the Rayleigh quotient $R(\mathbf{x})$ at $\mathbf{x} = (1, 0)^\top$, $(1, 1)^\top$, $(1, -1)^\top$, $(2, 1)^\top$. Find the eigenvalues of $A$ and confirm every value lies between them. Which inputs achieve equality, and why?

## B. Core (do all)

**B1.** (The adjoint depends on the inner product.) Let $W$ be a real symmetric positive definite $n\times n$ matrix and define $\langle\mathbf{u},\mathbf{v}\rangle_W = \mathbf{u}^\top W\mathbf{v}$.
(a) Show that $\langle\cdot,\cdot\rangle_W$ is an inner product on $\mathbb{R}^n$.
(b) Find the adjoint $A^\star$ of a matrix $A$ with respect to $\langle\cdot,\cdot\rangle_W$, i.e. the matrix with $\langle A\mathbf{u},\mathbf{v}\rangle_W = \langle\mathbf{u}, A^\star\mathbf{v}\rangle_W$ for all $\mathbf{u},\mathbf{v}$. Check that it reduces to $A^\top$ when $W = I$.
(c) Show that $A$ is self-adjoint with respect to $\langle\cdot,\cdot\rangle_W$ if and only if $WA$ is symmetric.
(d) Let $K$ be symmetric and $M$ symmetric positive definite. Show that $A = M^{-1}K$ is self-adjoint with respect to $\langle\cdot,\cdot\rangle_M$, and conclude that the generalized eigenproblem $K\mathbf{v} = \lambda M\mathbf{v}$ has real eigenvalues and eigenvectors that are orthogonal in the $M$-inner product whenever the eigenvalues are distinct. (This is the normal-mode problem of CM-6.)

**B2.** (Symmetric and antisymmetric matrices are orthogonal complements.) In $\mathbb{R}^{n\times n}$ with the Frobenius inner product, let $\mathcal S$ be the symmetric matrices and $\mathcal K$ the antisymmetric ones ($K^\top = -K$).
(a) Prove $\langle S, K\rangle_F = 0$ for every $S\in\mathcal S$, $K\in\mathcal K$.
(b) Show every $A$ is uniquely $S + K$ with $S\in\mathcal S$, $K\in\mathcal K$, and give the formulas.
(c) Compute $\dim\mathcal S$ and $\dim\mathcal K$ and conclude $\mathcal K = \mathcal S^\perp$.
(d) Deduce that $\tfrac12(A + A^\top)$ is the symmetric matrix closest to $A$ in Frobenius norm. (Check numerically in the project's "Your turn.")

**B3.** (The second-derivative test, and when it fails.)
(a) For $f(x, y) = x^3 - 3xy^2$, compute the gradient and the Hessian matrix $H(x, y)$. Find all critical points.
(b) Find the eigenvalues of $H(x,y)$ at a general point. At which points is $H$ indefinite? What is $H$ at the critical point, and what does Theorem 3.19 conclude there?
(c) Classify the critical point anyway, by examining $f$ along the $x$-axis. (This surface is the "monkey saddle.")
(d) For $h(x, y) = x^3 - 3x + y^2$, find the critical points, the Hessian at each, its eigenvalues, and the classification.

**B4.** (Every inner product on $\mathbb{R}^n$ is a positive definite matrix.) Let $\langle\cdot,\cdot\rangle$ be any inner product on $\mathbb{R}^n$ and define $G_{ij} = \langle\mathbf{e}_i,\mathbf{e}_j\rangle$.
(a) Show $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^\top G\mathbf{v}$ for all $\mathbf{u},\mathbf{v}$.
(b) Show $G$ is symmetric positive definite.
(c) Conversely, show that for any symmetric positive definite $G$, $\mathbf{u}^\top G\mathbf{v}$ is an inner product (this is B1(a) again; give the one-line reason using Theorem 3.15).
(d) State and prove the complex version, being careful about which slot carries the conjugate.

**B5.** ($A^\top A$ and $AA^\top$ share nonzero eigenvalues; SVD preview.) Let $A$ be real $m\times n$.
(a) Show that if $A^\top A\mathbf{v} = \lambda\mathbf{v}$ with $\lambda\ne0$ and $\mathbf{v}\ne\mathbf{0}$, then $\mathbf{w} = A\mathbf{v}$ is nonzero and $AA^\top\mathbf{w} = \lambda\mathbf{w}$.
(b) Show the map $\mathbf{v}\mapsto A\mathbf{v}$ from the $\lambda$-eigenspace of $A^\top A$ to the $\lambda$-eigenspace of $AA^\top$ is injective, and that $\mathbf{w}\mapsto A^\top\mathbf{w}$ goes back. Conclude the nonzero eigenvalues agree *with multiplicity*.
(c) For $m > n$, how many more zero eigenvalues does $AA^\top$ have than $A^\top A$ when $A$ has full column rank? Verify numerically in the project's "Your turn."

**B6.** (Exponentials of Hermitian matrices are unitary.) Let $A = \sum_i\lambda_iP_i$ be Hermitian with spectral decomposition (3.10), $P_i = \mathbf{q}_i\mathbf{q}_i^\dagger$.
(a) Define $e^{iA} = \sum_ie^{i\lambda_i}P_i$ by (3.11). Prove $e^{iA}$ is unitary. Where is the reality of $\lambda_i$ used?
(b) For a unit vector $\mathbf{n}$, show the projectors onto the $\pm1$ eigenspaces of $\mathbf{n}\cdot\boldsymbol\sigma$ are $P_\pm = \tfrac12(I\pm\mathbf{n}\cdot\boldsymbol\sigma)$, using only $(\mathbf{n}\cdot\boldsymbol\sigma)^2 = I$.
(c) Prove $e^{-i\theta\,\mathbf{n}\cdot\boldsymbol\sigma/2} = \cos\tfrac\theta2\,I - i\sin\tfrac\theta2\,\mathbf{n}\cdot\boldsymbol\sigma$.
(d) What is this matrix at $\theta = 2\pi$? (QC-1 will tell you it is the rotation of the Bloch sphere by $\theta$ about $\mathbf{n}$; LA-12 explains the sign.)

**B7.** (Code: why Jacobi converges.) For a symmetric $A$ let $\operatorname{off}(A)^2 = \sum_{i\ne j}a_{ij}^2$.
(a) Prove $\|Q^\top AQ\|_F = \|A\|_F$ for any orthogonal $Q$. (Hint: $\|A\|_F^2 = \operatorname{tr}(A^\top A)$ and the trace is cyclic.)
(b) Let $G$ be the Jacobi rotation in the $(p, q)$ plane: $G = I$ except $G_{pp} = G_{qq} = \cos\theta$, $G_{pq} = -\sin\theta$, $G_{qp} = \sin\theta$ (the sign convention of `project.py`; the opposite convention flips the sign of $\theta$). With $\tan2\theta = 2a_{pq}/(a_{pp} - a_{qq})$ it zeroes $a'_{pq}$ in $A' = G^\top AG$. Show it changes only rows and columns $p, q$, that the $2\times2$ block $\begin{pmatrix}a_{pp}&a_{pq}\\a_{pq}&a_{qq}\end{pmatrix}$ undergoes an orthogonal similarity, and that for $i\ne p, q$ the pair $(a_{ip}, a_{iq})$ is rotated rigidly. Conclude
$$\operatorname{off}(A')^2 = \operatorname{off}(A)^2 - 2a_{pq}^2.$$
(c) Why does (b) imply the classical Jacobi method (always zero the largest off-diagonal entry) converges, in the sense $\operatorname{off}(A_k)\to0$? (Hint: the largest of the $n(n-1)$ off-diagonal squares is at least their average.) Then use Weyl's inequality (the corollary to Theorem 3.18) to show the sorted diagonal of $A_k$ converges to the sorted eigenvalues of $A$.
(d) (Code.) `project.py` Part 3 asserts the identity in (b) after every rotation of the classical method. Write `cyclic_jacobi_eigh(A)`, which sweeps through all pairs $p < q$ in row order and rotates whenever $|a_{pq}|$ is above the tolerance, with no search for the largest entry. Put the same assertion inside its loop, and compare its rotation count with the classical method on the project's random $8\times8$ matrix `M8`. Which needs more rotations, why does the identity in (b) still hold, and why is the cyclic version still the one people use?

**B8.** (The step people skip.) In the proof of Theorem 3.11, $\mathbf{q}_1$ is a unit eigenvector with $A\mathbf{q}_1 = \lambda_1\mathbf{q}_1$, $V$ is $n\times(n-1)$ with orthonormal columns spanning $\mathbf{q}_1^\perp$, and $B = V^\top AV$.
(a) Prove $B$ is symmetric.
(b) Prove that if $B\mathbf{z} = \mu\mathbf{z}$ then $\mathbf{x} = V\mathbf{z}$ satisfies $A\mathbf{x} = \mu\mathbf{x}$, $\mathbf{x}\perp\mathbf{q}_1$, and $\|\mathbf{x}\| = \|\mathbf{z}\|$. (You will need $VV^\top$ = the projection onto $\mathbf{q}_1^\perp$, from LA-5, and the invariance of $\mathbf{q}_1^\perp$.)
(c) Take $A = \begin{pmatrix}2&1\\0&1\end{pmatrix}$ (not symmetric). It has eigenvector $\mathbf{q}_1 = \mathbf{e}_1$. Compute $A\mathbf{e}_2$ and show $\mathbf{q}_1^\perp$ is not invariant. Find both eigenvectors of $A$ and show they are not orthogonal. Which line of the proof of Theorem 3.11 is the one that breaks?

## C. Challenge (do at least one)

**C1.** (Cauchy interlacing.) Let $A$ be symmetric $n\times n$ with eigenvalues $\lambda_1\ge\cdots\ge\lambda_n$, and let $B$ be the $(n-1)\times(n-1)$ matrix obtained by deleting the last row and column, with eigenvalues $\mu_1\ge\cdots\ge\mu_{n-1}$. Use Courant–Fischer (3.14) to prove
$$\lambda_k\ \ge\ \mu_k\ \ge\ \lambda_{k+1}\qquad(1\le k\le n-1).$$
Check it on Example 4.2's matrix and its leading $2\times2$ block. (Hint: $B = V^\top AV$ with $V$ the first $n-1$ columns of $I$, so $\mathbf{z}^\top B\mathbf{z}/\mathbf{z}^\top\mathbf{z} = R_A(V\mathbf{z})$; subspaces of $\mathbb{R}^{n-1}$ are subspaces of $\mathbb{R}^n$ inside the hyperplane $x_n = 0$.)

**C2.** (Polarization: the norm determines the inner product.)
(a) In a real inner product space prove $\langle u, v\rangle = \tfrac14\left(\|u+v\|^2 - \|u-v\|^2\right)$.
(b) In a complex inner product space with our convention (conjugate-linear first slot) prove
$$\langle u, v\rangle = \tfrac14\left(\|u+v\|^2 - \|u-v\|^2\right) + \tfrac i4\left(\|u - iv\|^2 - \|u + iv\|^2\right).$$
What changes in Axler's convention?
(c) Conclude: if a linear map $U$ satisfies $\|U\mathbf{x}\| = \|\mathbf{x}\|$ for all $\mathbf{x}$, then $\langle U\mathbf{u}, U\mathbf{v}\rangle = \langle\mathbf{u},\mathbf{v}\rangle$ for all $\mathbf{u},\mathbf{v}$, so $U$ is unitary. This completes Theorem 3.20.

**C3.** (The generalized eigenproblem; CM-6 preview.) Let $K$ be symmetric and $M$ symmetric positive definite, and consider $K\mathbf{v} = \lambda M\mathbf{v}$.
(a) By Theorem 3.15(iii) write $M = LL^\top$ with $L$ invertible. Substituting $\mathbf{w} = L^\top\mathbf{v}$, show the problem becomes $C\mathbf{w} = \lambda\mathbf{w}$ with $C = L^{-1}KL^{-\top}$ symmetric.
(b) Conclude that all generalized eigenvalues are real, and that the eigenvectors $\mathbf{v}_i = L^{-\top}\mathbf{w}_i$ can be chosen with $\mathbf{v}_i^\top M\mathbf{v}_j = \delta_{ij}$.
(c) Solve by hand for $M = \begin{pmatrix}1&0\\0&2\end{pmatrix}$, $K = \begin{pmatrix}2&-1\\-1&2\end{pmatrix}$ (two masses $1$ and $2$ coupled by unit springs): find both $\lambda$ and the eigenvectors, and verify $\mathbf{v}_1^\top M\mathbf{v}_2 = 0$ directly.
(d) Why does `np.linalg.eigh(K, M)` not exist but `scipy.linalg.eigh(K, M)` does, and what does it compute? Check your answer to (c) with it.

**C4.** (Normal ⇒ unitarily diagonalizable, $2\times2$ case.) Let $A$ be a $2\times2$ complex matrix with $AA^\dagger = A^\dagger A$.
(a) Take a unit eigenvector $\mathbf{v}$ of $A$ (why does one exist?) and a unit $\mathbf{w}\perp\mathbf{v}$. Let $U = \begin{pmatrix}\mathbf{v}&\mathbf{w}\end{pmatrix}$. Show $T = U^\dagger AU$ is upper triangular: $T = \begin{pmatrix}\lambda&b\\0&\mu\end{pmatrix}$.
(b) Show $T$ is normal.
(c) Compare the $(1,1)$ entries of $TT^\dagger$ and $T^\dagger T$ to conclude $b = 0$.
(d) Conclude $A = U\operatorname{diag}(\lambda,\mu)U^\dagger$. Explain in two sentences how this argument becomes the general proof once Schur's theorem ($A = UTU^\dagger$ with $T$ triangular, LA-10) is available.
