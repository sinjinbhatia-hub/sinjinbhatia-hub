# Checkpoint: Linear Algebra Core (after LA-9)

Closed book. Pen and paper; a calculator only for arithmetic. Budget: 3 hours. Then check against the answer key at the bottom. The pass line is honest: **every** Part A item correct, at least **8 of 10** in Part B, at least **2 of 4** in Part C. Below that, the module named next to the question is where to go back.

Why this exists: the advanced block (LA-10 to LA-14) and every other course in this curriculum assume these ideas are automatic, not merely familiar. Getting 60% here and moving on will cost you ten times the hours later.

---

## Part A — Must be automatic (10 items, 2 minutes each)

A1. (LA-1) Write $A\mathbf{x}$ for a $3 \times 2$ matrix as a combination of columns. What are the coefficients?

A2. (LA-2) What does it mean that a $4\times 4$ matrix has 3 pivots? Give two consequences.

A3. (LA-3) State the dimensions of the four fundamental subspaces of an $m \times n$ matrix of rank $r$. Which pairs are orthogonal complements, and in which space does each pair live?

A4. (LA-4) If $B = P^{-1}AP$, what do $A$ and $B$ have in common and why?

A5. (LA-5) Write the projection matrix onto the column space of $A$ (independent columns). What two properties characterize a projection matrix?

A6. (LA-6) Which of these change a determinant, and how: swapping two rows; adding 3 times row 1 to row 2; multiplying a row by 5; transposing?

A7. (LA-7) $A$ is $2\times2$ with eigenvalues $3$ and $-1$ and eigenvectors $\mathbf{x}_1, \mathbf{x}_2$. Write $A^{10}\mathbf{v}$ for $\mathbf{v} = 2\mathbf{x}_1 + 5\mathbf{x}_2$.

A8. (LA-8) State the spectral theorem for real symmetric matrices. What's the complex version?

A9. (LA-9) In $A = U\Sigma V^\top$, which columns of $U$ and $V$ span each of the four fundamental subspaces?

A10. (LA-9) What is the condition number $\kappa_2(A)$ in terms of singular values, and what does $\kappa = 10^6$ mean for solving $A\mathbf{x} = \mathbf{b}$?

## Part B — Compute and prove (10 items, 10–15 minutes each)

B1. (LA-2) Find $L$ and $U$ for $A = \begin{pmatrix} 2 & 1 & 1 \\ 4 & 3 & 3 \\ 8 & 7 & 9 \end{pmatrix}$ and solve $A\mathbf{x} = (4, 10, 24)^\top$ using them.

B2. (LA-3) Find bases for all four fundamental subspaces of $A = \begin{pmatrix} 1 & 2 & 0 & 1 \\ 2 & 4 & 1 & 3 \\ 1 & 2 & 1 & 2 \end{pmatrix}$. Confirm the dimensions add up correctly in each space.

B3. (LA-3) Prove: if $\mathbf{v}_1, \dots, \mathbf{v}_k$ are linearly independent and $\mathbf{w} \notin \operatorname{span}\{\mathbf{v}_i\}$, then $\mathbf{v}_1, \dots, \mathbf{v}_k, \mathbf{w}$ are linearly independent.

B4. (LA-4) Let $T: P_2 \to P_2$ be $T(p) = p' + p(0)$. Find its matrix in the basis $\{1, x, x^2\}$, its kernel, and its image. Verify rank–nullity.

B5. (LA-5) Find the least-squares line through $(0, 0), (1, 1), (2, 1), (3, 3)$. Compute the residual vector and verify it's orthogonal to both columns of the design matrix.

B6. (LA-5) Prove that if $A$ has independent columns then $A^\top A$ is invertible.

B7. (LA-6) Compute $\det\begin{pmatrix} 1 & 1 & 1 \\ a & b & c \\ a^2 & b^2 & c^2 \end{pmatrix}$ and factor it. When is it zero, and what does that say about polynomial interpolation?

B8. (LA-7) Solve $\dfrac{d\mathbf{u}}{dt} = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}\mathbf{u}$, $\mathbf{u}(0) = (1, 0)^\top$. Is the origin stable?

B9. (LA-8) Orthogonally diagonalize $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$. Sketch the level set $\mathbf{x}^\top A \mathbf{x} = 1$ with its principal axes and semi-axis lengths.

B10. (LA-9) Compute the SVD of $A = \begin{pmatrix} 3 & 0 \\ 4 & 5 \end{pmatrix}$ by hand (via $A^\top A$). Give $\|A\|_2$, $\|A\|_F$, and $\kappa_2(A)$.

## Part C — Understanding (4 items, 20 minutes each)

C1. (LA-5/LA-8) Prove Cauchy–Schwarz for an abstract inner product using only the axioms. Then explain in two sentences why the same proof gives the uncertainty principle in quantum mechanics (you may use: for Hermitian $A, B$ with $\langle\psi|A\psi\rangle = \langle\psi|B\psi\rangle = 0$, $\sigma_A^2 \sigma_B^2 = \langle A\psi|A\psi\rangle \langle B\psi|B\psi\rangle$).

C2. (LA-8) Prove that a real symmetric matrix has only real eigenvalues, and that eigenvectors for distinct eigenvalues are orthogonal. Then explain why a $2\times2$ rotation matrix (not symmetric) has neither property.

C3. (LA-9) State the Eckart–Young theorem. Explain, in the language of the four fundamental subspaces, what information the rank-$k$ truncation keeps and what it throws away. Why is this the right notion of "compression"?

C4. (LA-7/LA-9) The matrix $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ has one eigenvector but two singular vectors. Explain geometrically what each decomposition is telling you about this shear, and why the eigen picture fails here while the SVD picture doesn't.

---

## Answer key

### Part A

A1. $A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2$; the coefficients are the entries of $\mathbf{x}$.

A2. Rank 3. Consequences: $A$ is singular (not invertible); $N(A)$ has dimension 1, so $A\mathbf{x} = \mathbf{0}$ has nonzero solutions and $A\mathbf{x} = \mathbf{b}$ has either no solution or infinitely many; $\det A = 0$; $0$ is an eigenvalue; $C(A)$ is a 3-dimensional subspace of $\mathbb{R}^4$.

A3. $\dim C(A) = r$, $\dim N(A) = n - r$, $\dim C(A^\top) = r$, $\dim N(A^\top) = m - r$. In $\mathbb{R}^n$: $N(A) = C(A^\top)^\perp$. In $\mathbb{R}^m$: $N(A^\top) = C(A)^\perp$.

A4. Same eigenvalues, characteristic polynomial, determinant, trace, rank, and Jordan structure, because they represent the same linear map in two bases ($P$ is the change-of-basis matrix).

A5. $P = A(A^\top A)^{-1}A^\top$. $P^2 = P$ and $P^\top = P$.

A6. Swap: sign flips. Add a multiple of a row to another: unchanged. Multiply a row by 5: multiplied by 5. Transpose: unchanged.

A7. $A^{10}\mathbf{v} = 2\cdot 3^{10}\mathbf{x}_1 + 5\cdot(-1)^{10}\mathbf{x}_2 = 118098\,\mathbf{x}_1 + 5\mathbf{x}_2$.

A8. Every real symmetric $A$ equals $Q\Lambda Q^\top$ with $Q$ orthogonal and $\Lambda$ real diagonal: there is an orthonormal basis of eigenvectors with real eigenvalues. Complex: every Hermitian $A$ equals $U\Lambda U^\dagger$ with $U$ unitary, $\Lambda$ real.

A9. $C(A) = \operatorname{span}(\mathbf{u}_1..\mathbf{u}_r)$, $N(A^\top) = \operatorname{span}(\mathbf{u}_{r+1}..\mathbf{u}_m)$, $C(A^\top) = \operatorname{span}(\mathbf{v}_1..\mathbf{v}_r)$, $N(A) = \operatorname{span}(\mathbf{v}_{r+1}..\mathbf{v}_n)$.

A10. $\kappa_2 = \sigma_1/\sigma_n$. A relative perturbation of $10^{-k}$ in $\mathbf{b}$ can produce a relative change of up to $10^{6-k}$ in $\mathbf{x}$; you lose about 6 digits.

### Part B

B1. Multipliers $\ell_{21} = 2$, $\ell_{31} = 4$: rows become $(2,1,1), (0,1,1), (0,3,5)$. Then $\ell_{32} = 3$: $(0,0,2)$. So $U = \begin{pmatrix}2&1&1\\0&1&1\\0&0&2\end{pmatrix}$, $L = \begin{pmatrix}1&0&0\\2&1&0\\4&3&1\end{pmatrix}$. Forward: $L\mathbf{c} = \mathbf{b}$: $c = (4, 2, 2)^\top$. Back: $U\mathbf{x} = \mathbf{c}$: $x_3 = 1$, $x_2 = 1$, $x_1 = 1$. $\mathbf{x} = (1,1,1)^\top$.

B2. RREF: $R = \begin{pmatrix}1&2&0&1\\0&0&1&1\\0&0&0&0\end{pmatrix}$, pivots in columns 1, 3; $r = 2$. $C(A)$: columns 1 and 3 of $A$: $\{(1,2,1)^\top, (0,1,1)^\top\}$. $C(A^\top)$: nonzero rows of $R$: $\{(1,2,0,1)^\top, (0,0,1,1)^\top\}$. $N(A)$: free variables $x_2, x_4$: $\{(-2,1,0,0)^\top, (-1,0,-1,1)^\top\}$. $N(A^\top)$: solve $A^\top\mathbf{y} = \mathbf{0}$: row 3 of $A$ = row 2 − row 1, so $\mathbf{y} = (1, -1, 1)^\top$ up to scale... check: $A^\top(1,-1,1)^\top$: column 1: $1 - 2 + 1 = 0$ ✓; column 2: $2-4+2=0$ ✓; column 3: $0-1+1=0$ ✓; column 4: $1-3+2=0$ ✓. Dimensions: in $\mathbb{R}^4$, $2 + 2 = 4$ ✓; in $\mathbb{R}^3$, $2 + 1 = 3$ ✓.

B3. Suppose $c_1\mathbf{v}_1 + \cdots + c_k\mathbf{v}_k + d\mathbf{w} = \mathbf{0}$. If $d \neq 0$, then $\mathbf{w} = -\frac{1}{d}\sum c_i\mathbf{v}_i \in \operatorname{span}$, contradiction. So $d = 0$, and then $\sum c_i\mathbf{v}_i = \mathbf{0}$ forces all $c_i = 0$ by independence. ∎

B4. $T(1) = 0 + 1 = 1$; $T(x) = 1 + 0 = 1$; $T(x^2) = 2x + 0 = 2x$. Matrix (columns are coordinates of images): $\begin{pmatrix}1&1&0\\0&0&2\\0&0&0\end{pmatrix}$. Rank 2. Kernel: $T(a + bx + cx^2) = (a + b) + 2cx = 0 \Rightarrow c = 0, b = -a$: kernel is $\operatorname{span}\{1 - x\}$, dimension 1. Image: $\operatorname{span}\{1, x\}$, dimension 2. $1 + 2 = 3 = \dim P_2$ ✓.

B5. Design matrix $A = \begin{pmatrix}1&0\\1&1\\1&2\\1&3\end{pmatrix}$, $\mathbf{b} = (0,1,1,3)^\top$. $A^\top A = \begin{pmatrix}4&6\\6&14\end{pmatrix}$, $A^\top\mathbf{b} = (5, 12)^\top$. Solve: $\det = 56 - 36 = 20$; $\hat{\mathbf{x}} = \frac{1}{20}\begin{pmatrix}14&-6\\-6&4\end{pmatrix}(5,12)^\top = \frac{1}{20}(70 - 72, -30 + 48)^\top = (-0.1, 0.9)^\top$. Line: $y = -0.1 + 0.9x$. Fitted: $(-0.1, 0.8, 1.7, 2.6)$. Residual $\mathbf{r} = (0.1, 0.2, -0.7, 0.4)^\top$. $\mathbf{r}\cdot(1,1,1,1) = 0$ ✓; $\mathbf{r}\cdot(0,1,2,3) = 0.2 - 1.4 + 1.2 = 0$ ✓.

B6. Suppose $A^\top A\mathbf{x} = \mathbf{0}$. Then $\mathbf{x}^\top A^\top A\mathbf{x} = \|A\mathbf{x}\|^2 = 0$, so $A\mathbf{x} = \mathbf{0}$, so $\mathbf{x} = \mathbf{0}$ by independence of columns. Thus $N(A^\top A) = \{\mathbf{0}\}$ and the square matrix $A^\top A$ is invertible. ∎

B7. Vandermonde: $(b-a)(c-a)(c-b)$. (Subtract $a\times$row 2 from row 3 and $a\times$row 1 from row 2, expand along column 1, factor.) Zero iff two of $a, b, c$ coincide. Interpolation: a polynomial of degree $\le 2$ through 3 points with distinct $x$-values exists and is unique because this matrix (transposed) is invertible exactly then.

B8. $\det(A - \lambda I) = \lambda^2 + 3\lambda + 2 = (\lambda+1)(\lambda+2)$. $\lambda = -1$: eigenvector $(1, -1)^\top$. $\lambda = -2$: $(1, -2)^\top$. $\mathbf{u}(0) = c_1(1,-1) + c_2(1,-2) = (1,0)$: $c_1 + c_2 = 1$, $-c_1 - 2c_2 = 0 \Rightarrow c_2 = -1, c_1 = 2$. $\mathbf{u}(t) = 2e^{-t}(1,-1)^\top - e^{-2t}(1,-2)^\top$, i.e. $u_1 = 2e^{-t} - e^{-2t}$, $u_2 = -2e^{-t} + 2e^{-2t}$. Both eigenvalues negative: stable (a stable node). This is an overdamped oscillator $\ddot{x} + 3\dot{x} + 2x = 0$.

B9. Eigenvalues: $(2-\lambda)^2 - 1 = 0 \Rightarrow \lambda = 1, 3$. $\lambda = 1$: $(1,1)^\top/\sqrt2$. $\lambda = 3$: $(1,-1)^\top/\sqrt2$. $A = Q\Lambda Q^\top$ with $Q = \frac{1}{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}$, $\Lambda = \operatorname{diag}(1, 3)$. The level set $\mathbf{x}^\top A\mathbf{x} = 1$ is an ellipse with axes along $(1,1)$ (semi-axis $1/\sqrt{1} = 1$) and $(1,-1)$ (semi-axis $1/\sqrt3 \approx 0.577$). Long axis along $(1,1)$.

B10. $A^\top A = \begin{pmatrix}25&20\\20&25\end{pmatrix}$, eigenvalues $45, 5$; $\sigma_1 = \sqrt{45} = 3\sqrt5$, $\sigma_2 = \sqrt5$. $\mathbf{v}_1 = (1,1)^\top/\sqrt2$, $\mathbf{v}_2 = (1,-1)^\top/\sqrt2$. $\mathbf{u}_1 = A\mathbf{v}_1/\sigma_1 = (3, 9)^\top/(\sqrt2\cdot3\sqrt5) = (1,3)^\top/\sqrt{10}$. $\mathbf{u}_2 = A\mathbf{v}_2/\sigma_2 = (3,-1)^\top/(\sqrt2\cdot\sqrt5) = (3,-1)^\top/\sqrt{10}$. Check $\mathbf{u}_1\cdot\mathbf{u}_2 = 0$ ✓. $\|A\|_2 = 3\sqrt5 \approx 6.708$; $\|A\|_F = \sqrt{9+16+25} = \sqrt{50} = \sqrt{45 + 5}$ ✓; $\kappa_2 = 3\sqrt5/\sqrt5 = 3$.

### Part C

C1. For $\mathbf{v}\ne0$ and any scalar $t$: $0 \le \langle \mathbf{u} - t\mathbf{v}, \mathbf{u} - t\mathbf{v}\rangle = \|\mathbf{u}\|^2 - 2t\operatorname{Re}\langle\mathbf{u},\mathbf{v}\rangle + |t|^2\|\mathbf{v}\|^2$ (real case: drop Re). Choosing $t = \langle\mathbf{v},\mathbf{u}\rangle/\|\mathbf{v}\|^2$ gives $0 \le \|\mathbf{u}\|^2 - |\langle\mathbf{u},\mathbf{v}\rangle|^2/\|\mathbf{v}\|^2$, i.e. $|\langle\mathbf{u},\mathbf{v}\rangle| \le \|\mathbf{u}\|\|\mathbf{v}\|$. Uncertainty: apply this to $\mathbf{u} = A\psi$, $\mathbf{v} = B\psi$: $\sigma_A^2\sigma_B^2 \ge |\langle A\psi|B\psi\rangle|^2 \ge (\operatorname{Im}\langle A\psi|B\psi\rangle)^2 = \frac14|\langle\psi|[A,B]\psi\rangle|^2$. The uncertainty principle is Cauchy–Schwarz plus the observation that the commutator is the imaginary part of the inner product. (QM-4 does the last step.)

C2. Let $A\mathbf{x} = \lambda\mathbf{x}$ with $\mathbf{x} \ne 0$ possibly complex. Then $\bar{\mathbf{x}}^\top A\mathbf{x} = \lambda\|\mathbf{x}\|^2$. Also $\bar{\mathbf{x}}^\top A\mathbf{x} = (A\bar{\mathbf{x}})^\top\mathbf{x}$ (since $A = A^\top$) $= \overline{(A\mathbf{x})}^\top\mathbf{x}$ (since $A$ is real) $= \bar\lambda\|\mathbf{x}\|^2$. So $\lambda = \bar\lambda$. Orthogonality: if $A\mathbf{x} = \lambda\mathbf{x}$, $A\mathbf{y} = \mu\mathbf{y}$, $\lambda\ne\mu$: $\lambda\,\mathbf{x}^\top\mathbf{y} = (A\mathbf{x})^\top\mathbf{y} = \mathbf{x}^\top A\mathbf{y} = \mu\,\mathbf{x}^\top\mathbf{y}$, so $(\lambda-\mu)\mathbf{x}^\top\mathbf{y} = 0$, so $\mathbf{x}^\top\mathbf{y} = 0$. ∎ A rotation by $\theta \ne 0, \pi$ moves every real vector, so it has no real eigenvectors; its eigenvalues $e^{\pm i\theta}$ are complex. The proof above used $A = A^\top$ in the step $(A\mathbf{x})^\top\mathbf{y} = \mathbf{x}^\top A\mathbf{y}$, which fails for a rotation ($R^\top = R^{-1} \ne R$).

C3. Eckart–Young: among all matrices of rank $\le k$, $A_k = \sum_{i\le k}\sigma_i\mathbf{u}_i\mathbf{v}_i^\top$ minimizes $\|A - B\|$ in both the spectral and Frobenius norms, with errors $\sigma_{k+1}$ and $\sqrt{\sum_{i>k}\sigma_i^2}$. Truncation keeps the $k$ directions in $C(A^\top)$ (the $\mathbf{v}_i$) along which $A$ stretches the most, and the $k$ directions in $C(A)$ (the $\mathbf{u}_i$) they land on; it throws away the input directions with small stretch (treating them as if they were in $N(A)$) and the corresponding output directions. It's the right notion because the error is measured in the norm that counts total (or worst-case) stretch, and the discarded stretches are the smallest possible.

C4. The shear fixes the $x$-axis ($\lambda = 1$, eigenvector $\mathbf{e}_1$) and has no other invariant direction: it has one eigenvector because every other vector gets tilted. Eigen-analysis asks "which directions are only scaled?" and the honest answer is "one." The SVD asks a different question: "which orthonormal input basis is mapped to an orthonormal output basis, scaled?" That always has an answer: $A^\top A = \begin{pmatrix}1&1\\1&2\end{pmatrix}$ has eigenvalues $(3\pm\sqrt5)/2$, so $\sigma_1 \approx 1.618$, $\sigma_2 \approx 0.618$ (the golden ratio and its inverse). The unit circle maps to an ellipse with those semi-axes, tilted; the right singular vectors are the preimages of the ellipse's axes. The eigen picture fails because it demands input and output directions coincide; the SVD lets them differ by a rotation, which is exactly what a shear does.
