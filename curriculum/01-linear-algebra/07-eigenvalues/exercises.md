# LA-7 Exercises

Conventions as in the notes: $p_A(\lambda) = \det(A - \lambda I)$; Markov matrices are column-stochastic with $\mathbf{u}_{k+1} = A\mathbf{u}_k$; $T = \operatorname{tr}A$, $D = \det A$ for $2\times2$ matrices; eigenvectors are nonzero and only defined up to scale, so give any convenient representative.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** $A = \begin{pmatrix}4 & 1\\ 2 & 3\end{pmatrix}$. Find the eigenvalues from (3.3) and an eigenvector for each by elimination on $A - \lambda I$. Check $A\mathbf{x} = \lambda\mathbf{x}$ for both, and check Theorem 3.4 (sum and product).

**A2.** $U = \begin{pmatrix}2 & 7 & -3\\ 0 & -1 & 4\\ 0 & 0 & 5\end{pmatrix}$. Without computing a determinant, write down the eigenvalues of $U$, $U^2$, $U^{-1}$, $U + 3I$, and $U^\top$. Is $\mathbf{e}_1$ an eigenvector of each of these five matrices? Then find eigenvectors of $U$ for its other two eigenvalues.

**A3.** $A = \begin{pmatrix}0.8 & 0.3\\ 0.2 & 0.7\end{pmatrix}$. Verify it is Markov. Find both eigenvalues and the steady state. Starting from $\mathbf{u}_0 = (1, 0)^\top$, compute $\mathbf{u}_1$ and $\mathbf{u}_2$ directly, then write $\mathbf{u}_k$ in closed form via (3.8) and confirm it reproduces $\mathbf{u}_1, \mathbf{u}_2$.

**A4.** $R = \begin{pmatrix}0 & -1\\ 1 & 0\end{pmatrix}$ (rotation by $90°$). Find its eigenvalues and an eigenvector for each, and verify $R\mathbf{x} = \lambda\mathbf{x}$ with complex arithmetic. Then: $R^2 = -I$ has eigenvalue $-1$ with *every* vector as eigenvector. Reconcile this with Theorem 3.5(b).

**A5.** Compute $e^{At}$ from the definition (3.12) for
(a) $A = \begin{pmatrix}2 & 0\\ 0 & -3\end{pmatrix}$; (b) $N = \begin{pmatrix}0 & 1\\ 0 & 0\end{pmatrix}$; (c) $A = \begin{pmatrix}2 & 1\\ 0 & 2\end{pmatrix}$.
For (c), use Theorem 3.15(c) with $A = 2I + N$, and check $\frac{d}{dt}e^{At}\big|_{t=0} = A$.

## B. Core (do all)

**B1.** (Eigenvalues of related matrices.)
(a) Prove: if $A$ is invertible, then $\mathbf{x}$ is an eigenvector of $A$ with eigenvalue $\lambda$ if and only if $\mathbf{x}$ is an eigenvector of $A^{-1}$ with eigenvalue $1/\lambda$. Prove also that $0$ is an eigenvalue of $A$ if and only if $A$ is singular.
(b) Let $A = \begin{pmatrix}0 & 1\\ 0 & 0\end{pmatrix}$, $B = \begin{pmatrix}0 & 0\\ 1 & 0\end{pmatrix}$. Compute the eigenvalues of $A$, $B$, $A + B$, $AB$. Conclude that eigenvalues of sums and products are not sums and products of eigenvalues.
(c) Suppose $A$ and $B$ are **simultaneously diagonalizable**: $X^{-1}AX = \Lambda_A$ and $X^{-1}BX = \Lambda_B$ with the *same* $X$. Prove that the eigenvalues of $A + B$ are $\lambda_i + \mu_i$, the eigenvalues of $AB$ are $\lambda_i\mu_i$, and $AB = BA$.

**B2.** Solve $\dot{\mathbf{u}} = \begin{pmatrix}1 & 2\\ 2 & 1\end{pmatrix}\mathbf{u}$ with $\mathbf{u}(0) = (3, 1)^\top$ using (3.11). Verify $\dot{\mathbf{u}}(0) = A\mathbf{u}(0)$. Classify the phase portrait from $(T, D)$ before you compute eigenvalues, and say along which line the trajectory runs as $t\to+\infty$ and as $t\to-\infty$. Why are the two eigenvectors orthogonal?

**B3.** $A = \begin{pmatrix}2 & 1\\ 1 & 2\end{pmatrix}$ (the matrix of §2).
(a) Diagonalize: find $X$, $\Lambda$, $X^{-1}$.
(b) Derive the closed form $A^k = \tfrac12\begin{pmatrix}3^k + 1 & 3^k - 1\\ 3^k - 1 & 3^k + 1\end{pmatrix}$ from (3.7), and verify it for $k = 5$ by direct multiplication (you should get $122$ and $121$).
(c) Compute $\lim_{k\to\infty}A^k/3^k$. Show the limit is a projection matrix ($P^2 = P$, $P^\top = P$) and say onto what. Connect this to power iteration (§3.6) and to the projection formula of LA-1, Example 4.5.

**B4.** (Markov theory.) Let $A$ be an $n\times n$ Markov matrix.
(a) Prove that $A^k$ is Markov for every $k \ge 1$.
(b) Prove that if $A\mathbf{x} = \lambda\mathbf{x}$ with $\lambda \ne 1$, then $\mathbf{1}^\top\mathbf{x} = 0$: the entries of $\mathbf{x}$ sum to zero.
(c) Suppose $A$ is diagonalizable and $\lambda_1 = 1$ is a **simple** eigenvalue, so $\lambda_i \ne 1$ for $i \ge 2$. Using (b), explain why in $\mathbf{u}_k = c_1\mathbf{x}_1 + \sum_{i\ge2}c_i\lambda_i^k\mathbf{x}_i$ the coefficient $c_1$ is the same for every probability vector $\mathbf{u}_0$, and what the transient terms do to the probability mass. Then show the simplicity hypothesis is needed: for the Markov matrix $A = I$ ($2\times2$) with $\mathbf{x}_1 = \mathbf{e}_1$, $\mathbf{x}_2 = \mathbf{e}_2$, find $c_1$ for $\mathbf{u}_0 = (a, 1 - a)^\top$.
(d) For the weather matrix $A = \begin{pmatrix}0.9 & 0.5\\ 0.1 & 0.5\end{pmatrix}$ of Example 4.1, compute $A^k$ in closed form via (3.7) and find $\lim_{k\to\infty}A^k$. Interpret the columns of the limit.

**B5.** $A = \begin{pmatrix}1 & -1\\ 1 & 1\end{pmatrix}$, $\dot{\mathbf{u}} = A\mathbf{u}$, $\mathbf{u}(0) = (1, 0)^\top$.
(a) Classify from $(T, D)$. Find the eigenvalues and one complex eigenvector.
(b) Produce two real solutions via (3.19) and match the initial condition.
(c) Recognize $A = aI + bJ$ and recover the same answer from (3.20). Which way does the trajectory rotate, and does it spiral in or out?
(d) Show that any real $2\times2$ matrix $B$ with eigenvalue $a + ib$ ($b\ne0$) and eigenvector $\mathbf{r} + i\mathbf{s}$ satisfies $B\begin{pmatrix}\mathbf{r} & \mathbf{s}\end{pmatrix} = \begin{pmatrix}\mathbf{r} & \mathbf{s}\end{pmatrix}\begin{pmatrix}a & b\\ -b & a\end{pmatrix}$, so $B$ is similar to a rotation-with-scaling matrix. (Compare real and imaginary parts of $B(\mathbf{r} + i\mathbf{s}) = (a + ib)(\mathbf{r} + i\mathbf{s})$.)

**B6.** (Mechanics: the generalized eigenproblem.) Let $K$ be real symmetric and $M$ real symmetric positive definite, and suppose $K\mathbf{v}_1 = \omega_1^2M\mathbf{v}_1$, $K\mathbf{v}_2 = \omega_2^2M\mathbf{v}_2$ with $\omega_1^2 \ne \omega_2^2$.
(a) Prove $\mathbf{v}_1^\top M\mathbf{v}_2 = 0$ ($M$-orthogonality).
(b) Prove that every $\omega^2$ solving $K\mathbf{v} = \omega^2M\mathbf{v}$ (with $\mathbf{v}$ possibly complex) is real. You may use that $\bar{\mathbf{v}}^\top M\mathbf{v} > 0$ for $\mathbf{v}\ne\mathbf{0}$; prove this from positive definiteness of $M$ on real vectors by writing $\mathbf{v} = \mathbf{r} + i\mathbf{s}$.
(c) Verify (a) numerically for Example 4.3 and show that the ordinary dot product $\mathbf{v}_1^\top\mathbf{v}_2$ is *not* zero.
(d) Explain in one sentence why Theorem 3.17(b) does not apply directly to $M^{-1}K$.
(e) For diagonal $M = \operatorname{diag}(m_i)$ let $M^{1/2} = \operatorname{diag}(\sqrt{m_i})$. Substitute $\mathbf{v} = M^{-1/2}\mathbf{w}$ into $K\mathbf{v} = \omega^2M\mathbf{v}$ and show it becomes an ordinary eigenproblem for the *symmetric* matrix $\tilde K = M^{-1/2}KM^{-1/2}$. Show that $M$-orthogonality of the $\mathbf{v}$'s is ordinary orthogonality of the $\mathbf{w}$'s. Compute $\tilde K$ for Example 4.3 and check its trace and determinant against $\omega_1^2 + \omega_2^2$ and $\omega_1^2\omega_2^2$.

**B7.** (Defective, and nearly defective.) $A = \begin{pmatrix}3 & 1\\ 0 & 3\end{pmatrix}$, $B = \begin{pmatrix}3 & 0\\ 0 & 3\end{pmatrix}$.
(a) Show $A$ and $B$ have the same characteristic polynomial, that $A$ is not diagonalizable, and that $A$ is not similar to $B$.
(b) Compute $A^k$ and $e^{At}$ in closed form.
(c) Solve $\dot{\mathbf{u}} = A\mathbf{u}$, $\mathbf{u}(0) = (1, 2)^\top$, and verify by substitution.
(d) Let $A_\varepsilon = \begin{pmatrix}3 & 1\\ 0 & 3 + \varepsilon\end{pmatrix}$ with $\varepsilon \ne 0$ small. Find its eigenvalues and eigenvectors, the angle between the eigenvectors, and $X_\varepsilon^{-1}$. What happens as $\varepsilon\to0$, and why does this make $A_\varepsilon^k = X_\varepsilon\Lambda_\varepsilon^kX_\varepsilon^{-1}$ a bad way to compute $A_\varepsilon^k$ in floating point?

**B8.** (Code.) `project.py`, Part 1, implements `power_iteration` and `inverse_iteration`, and contains a reference `convergence_rate` that answers (a)–(c). Write your own version first, then compare.
(a) For $A = \begin{pmatrix}4 & 1\\ 2 & 3\end{pmatrix}$ starting from $\mathbf{v}_0 = (1, 0)^\top$, record the direction error $e_k = \|\mathbf{v}_k - \mathbf{x}_1\|$ (after flipping the sign of $\mathbf{v}_k$ to match $\mathbf{x}_1$, both unit length) and fit a straight line to $\log e_k$ against $k$ over the range where $e_k > 10^{-10}$. Compare the slope to $\log|\lambda_2/\lambda_1| = \log(2/5)$. Do the same for the symmetric $\begin{pmatrix}2 & 1\\ 1 & 2\end{pmatrix}$ (ratio $1/3$).
(b) Record also the Rayleigh-quotient error $|\lambda^{(k)} - \lambda_1|$. For the symmetric matrix you should find a slope close to $2\log|\lambda_2/\lambda_1|$, twice as steep. Explain this: write the unit iterate as $\mathbf{v}_k \propto \mathbf{x}_1 + \varepsilon_k\mathbf{x}_2$ with $\mathbf{x}_1\perp\mathbf{x}_2$ orthonormal, expand the Rayleigh quotient to second order in $\varepsilon_k$, and say what is different in the nonsymmetric case: expand it to first order with unit eigenvectors and $c = \mathbf{x}_1^\top\mathbf{x}_2 \ne 0$.
(c) Extend your function with a `shift` argument that runs inverse iteration (solve with $A - \sigma I$ in place of the multiply; the target is now the eigenvector whose $\lambda$ is closest to $\sigma$). For $A = \begin{pmatrix}4 & 1\\ 2 & 3\end{pmatrix}$ and $\sigma = 2.2$, predict the slope from §3.6 and check it. How many decimal digits does each step gain, compared with plain power iteration on the same matrix?

## C. Challenge (do at least one)

**C1.** (Cayley–Hamilton, diagonalizable case.) For a polynomial $p(\lambda) = \sum_jc_j\lambda^j$ define $p(A) = \sum_jc_jA^j$ (with $A^0 = I$). Prove: if $A$ is diagonalizable, then $p_A(A) = 0$, the matrix satisfies its own characteristic equation. (Hint: $p_A(A) = Xp_A(\Lambda)X^{-1}$.)
(a) Use it for $2\times2$: from $A^2 - TA + DI = 0$ derive $A^{-1} = (TI - A)/D$ when $D\ne0$, and check against the formula from LA-2.
(b) Verify $p_A(A) = 0$ numerically for a random $4\times4$ matrix (add a cell to `project.py`; `np.poly(A)` gives the coefficients of the monic version).
(c) Where does your proof fail for a defective $A$? (The theorem is still true; LA-11 proves it in general.)

**C2.** (Toward Perron–Frobenius.) Let $A$ be a Markov matrix with **all entries strictly positive**.
(a) Prove: if $\lambda$ is an eigenvalue with $|\lambda| = 1$, then $\lambda = 1$. (Follow the proof of Theorem 3.11(b) for $A^\top\mathbf{v} = \lambda\mathbf{v}$ and analyze when both inequalities are equalities. For the triangle inequality $|\sum_jz_j| \le \sum_j|z_j|$, equality holds iff all nonzero $z_j$ have the same argument.)
(b) Prove the eigenspace of $A^\top$ for $\lambda = 1$ is $\operatorname{span}\{\mathbf{1}\}$, and deduce that the eigenspace of $A$ for $\lambda = 1$ is one-dimensional. (Use $\operatorname{rank}(A - I) = \operatorname{rank}(A - I)^\top$, LA-3.)
(c) Conclude: if in addition $A$ is diagonalizable, then $A^k\mathbf{u}_0\to\mathbf{u}_\infty$ for every probability vector $\mathbf{u}_0$, where $\mathbf{u}_\infty$ is the unique steady state. Where did you use B4(b)?

**C3.** (The $n$-mass chain.) $n$ equal masses $m$ are connected in a line by $n + 1$ equal springs $k$, the two end springs attached to walls.
(a) Derive $M\ddot{\mathbf{x}} = -K\mathbf{x}$ with $M = mI$ and $K = k\,T_n$, where $T_n$ is the $n\times n$ matrix with $2$ on the diagonal, $-1$ on the two neighboring diagonals, and $0$ elsewhere.
(b) Let $\theta_j = \dfrac{j\pi}{n+1}$ for $j = 1, \dots, n$ and define $\mathbf{v}_j$ by $(\mathbf{v}_j)_i = \sin(i\theta_j)$, $i = 1, \dots, n$. Prove that $T_n\mathbf{v}_j = (2 - 2\cos\theta_j)\mathbf{v}_j$. (Hint: $\sin((i-1)\theta) + \sin((i+1)\theta) = 2\sin(i\theta)\cos\theta$; check that the boundary rows $i = 1$ and $i = n$ work because $\sin(0) = 0$ and $\sin((n+1)\theta_j) = 0$.) Explain why the $\mathbf{v}_j$ are nonzero, why the eigenvalues are distinct, and why the $\mathbf{v}_j$ are therefore orthogonal.
(c) Conclude $\omega_j = 2\sqrt{k/m}\,\sin\!\big(\tfrac{j\pi}{2(n+1)}\big)$. Check against the project's numbers for $n = 3$ (eigenvalues of $T_3$ are $2 - \sqrt2$, $2$, $2 + \sqrt2$).
(d) Let the spacing be $a$ and the total length $L = (n+1)a$. Show that for fixed $j$ and $n\to\infty$, $\omega_j \to \dfrac{j\pi}{L}\,a\sqrt{k/m}$. Interpret: these are the standing-wave frequencies of a string of length $L$ with wave speed $c = a\sqrt{k/m}$, which is where F2's Fourier series come from.
