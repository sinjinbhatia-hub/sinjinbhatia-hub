# LA-5 Exercises

Conventions as in the notes: $A$ is $m\times n$ with columns $\mathbf{a}_j$; $\hat{\mathbf{x}}$ solves the normal equations $A^\top A\hat{\mathbf{x}} = A^\top\mathbf{b}$; $\mathbf{p} = A\hat{\mathbf{x}}$ is the projection, $P$ the projection matrix, $\mathbf{e} = \mathbf{b} - \mathbf{p}$ the error (residual); $Q$ has orthonormal columns; in $A = QR$ the diagonal of $R$ is positive.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Let $\mathbf{u} = (1, -2, 2)^\top$ and $\mathbf{v} = (2, 2, 1)^\top$. (a) Show $\mathbf{u} \perp \mathbf{v}$ and verify Pythagoras (3.1) for $\mathbf{u} + \mathbf{v}$. (b) Find a nonzero $\mathbf{w}$ orthogonal to both. (c) Why are $\mathbf{u}, \mathbf{v}, \mathbf{w}$ automatically a basis of $\mathbb{R}^3$? Which theorem?

**A2.** Let $S = \operatorname{span}\{(1, 2, 2)^\top\} \subset \mathbb{R}^3$. Write $S^\perp$ as the null space of a matrix, find a basis of it, and check $\dim S + \dim S^\perp = 3$. Then compute $(S^\perp)^\perp$ directly and confirm it is $S$.

**A3.** Project $\mathbf{b} = (3, 0, 3)^\top$ onto the line through $\mathbf{a} = (1, 2, 2)^\top$ using (3.4). Give $\mathbf{p}$, $\mathbf{e}$, and the $3\times3$ matrix $P$. Verify $\mathbf{a}^\top\mathbf{e} = 0$, $P^2 = P$ (without multiplying out all nine entries: use $\mathbf{a}\mathbf{a}^\top\mathbf{a}\mathbf{a}^\top = (\mathbf{a}^\top\mathbf{a})\,\mathbf{a}\mathbf{a}^\top$), and $\|\mathbf{b}\|^2 = \|\mathbf{p}\|^2 + \|\mathbf{e}\|^2$.

**A4.** Find the least-squares line $b = C + Dt$ through $(0, 0), (1, 1), (2, 3)$. Write $A$, $A^\top A$, $A^\top\mathbf{b}$, solve, and give the residual vector. Check that the residual is orthogonal to both columns of $A$.

**A5.** Let $Q = \frac13\begin{pmatrix}1 & 2 & -2\\-2 & 2 & 1\\2 & 1 & 2\end{pmatrix}$ (its columns are $\mathbf{u}, \mathbf{v}, \mathbf{w}$ from A1, normalized). (a) Verify $Q^\top Q = I$ by computing the six distinct dot products. What is $Q^{-1}$? (b) For $\mathbf{x} = (1, 2, 3)^\top$ compute $Q\mathbf{x}$ and confirm $\|Q\mathbf{x}\| = \|\mathbf{x}\|$. (c) Expand $\mathbf{x}$ in the basis of columns of $Q$ using (3.13) and check Parseval.

**A6.** Run Gram–Schmidt on $\mathbf{a}_1 = (1, 2, 2)^\top$, $\mathbf{a}_2 = (0, 3, 0)^\top$. Give $\mathbf{q}_1, \mathbf{q}_2$, and $R$, and verify $\mathbf{a}_2 = r_{12}\mathbf{q}_1 + r_{22}\mathbf{q}_2$.

## B. Core (do all)

**B1.** For subspaces $S, T$ of $\mathbb{R}^n$, prove: (a) $S \cap S^\perp = \{\mathbf{0}\}$; (b) if $S \subseteq T$ then $T^\perp \subseteq S^\perp$; (c) $(S + T)^\perp = S^\perp \cap T^\perp$, where $S + T = \{\mathbf{s} + \mathbf{t} : \mathbf{s} \in S, \mathbf{t} \in T\}$.

**B2.** Let $A = \begin{pmatrix}1 & 1 & 1 & 1\\1 & 2 & 3 & 4\end{pmatrix}$.
(a) Find a basis of $N(A)$ and verify each basis vector is orthogonal to each row of $A$ (Theorem 3.9).
(b) Split $\mathbf{x} = (1, 0, 0, 0)^\top$ as $\mathbf{x}_{\text{row}} + \mathbf{x}_{\text{null}}$ with $\mathbf{x}_{\text{row}} \in C(A^\top)$ and $\mathbf{x}_{\text{null}} \in N(A)$. (Hint: $\mathbf{x}_{\text{row}}$ is the projection of $\mathbf{x}$ onto $C(A^\top)$; use (3.6) with $A^\top$ in the role of $A$.) Check $A\mathbf{x} = A\mathbf{x}_{\text{row}}$.

**B3.** Prove that $\operatorname{rank}(A^\top A) = \operatorname{rank}(A)$ and $C(A^\top A) = C(A^\top)$ for any $m\times n$ matrix $A$. (Use Theorem 3.11 and Rank–Nullity.) Deduce that $A^\top A = 0$ forces $A = 0$. State the corresponding facts for $AA^\top$.

**B4.** Let $S$ be the plane $x + y + z = 0$ in $\mathbb{R}^3$.
(a) Find the projection matrix $P$ onto $S$ from (3.8), using the basis $(1, -1, 0)^\top, (0, 1, -1)^\top$.
(b) Find it again as $I - \mathbf{n}\mathbf{n}^\top/\mathbf{n}^\top\mathbf{n}$ with $\mathbf{n} = (1, 1, 1)^\top$ (Theorem 3.13). They must agree; explain why in one sentence.
(c) Project $\mathbf{b} = (1, 2, 3)^\top$; give $\mathbf{p}$ and $\mathbf{e}$ and check $\mathbf{e} \parallel \mathbf{n}$.
(d) Verify $P^2 = P$ for this $P$ with a short computation, and show $\operatorname{tr}P = 2$. Then prove in general that $\operatorname{tr}\big(A(A^\top A)^{-1}A^\top\big) = n$ when $A$ is $m\times n$ with independent columns. (First prove $\operatorname{tr}(XY) = \operatorname{tr}(YX)$ for $X$ of size $p\times q$ and $Y$ of size $q\times p$.)

**B5.** (Constants, means, and centering.)
(a) Show that the least-squares fit of the one-parameter model $b_i \approx C$ (a column of ones) is the mean $C = \bar b$.
(b) Show that whenever the design matrix contains a column of ones, the residuals sum to zero.
(c) Refit Example 4.2 with the time column replaced by $t_i' = t_i - \bar t$. Show $A'^\top A'$ is diagonal, solve, and show the fitted line is the *same line* as before. Which theorem guarantees this?
(d) From (c), derive the statistics textbook formula $D = \sum_i(t_i - \bar t)(b_i - \bar b)\big/\sum_i(t_i - \bar t)^2$.

**B6.** (The gradient, by hand.) Let $f(\mathbf{x}) = \|A\mathbf{x} - \mathbf{b}\|^2$ with $A$ of size $m\times2$.
(a) Write $f$ out as an explicit polynomial in $x_1, x_2$ using $G = A^\top A$ (entries $g_{11}, g_{12}, g_{22}$) and $\mathbf{c} = A^\top\mathbf{b}$. Compute $\partial f/\partial x_1$ and $\partial f/\partial x_2$ and confirm they are the entries of $2(G\mathbf{x} - \mathbf{c})$.
(b) Give a second derivation of (3.10) using the directional derivative: expand $f(\mathbf{x} + h\mathbf{d})$ in powers of $h$.
(c) Suppose the columns of $A$ are dependent, so $A^\top A$ is singular. Does $f$ still have a minimizer? Is it unique? Describe the full set of minimizers.

**B7.** (Orthogonal matrices.)
(a) Prove that the product of two orthogonal matrices is orthogonal, and that the inverse of an orthogonal matrix is orthogonal.
(b) Show every $2\times2$ orthogonal matrix is either $\begin{pmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{pmatrix}$ or $\begin{pmatrix}\cos\theta & \sin\theta\\ \sin\theta & -\cos\theta\end{pmatrix}$ for some $\theta$, and that the second one fixes every vector on the line at angle $\theta/2$ (so it is the reflection across that line).
(c) Prove the converse of Theorem 3.18(b): if $Q$ is square and $\|Q\mathbf{x}\| = \|\mathbf{x}\|$ for every $\mathbf{x}$, then $Q^\top Q = I$. (Hint: the polarization identity $\mathbf{x}^\top\mathbf{y} = \tfrac14\big(\|\mathbf{x} + \mathbf{y}\|^2 - \|\mathbf{x} - \mathbf{y}\|^2\big)$.)

**B8.** (Code: weighted least squares.) Readings with different reliabilities should not count equally. Given weights $w_i > 0$, the weighted problem is to minimize $\sum_i w_i\big((A\mathbf{x})_i - b_i\big)^2$.
(a) Show this equals $\|W^{1/2}(A\mathbf{x} - \mathbf{b})\|^2$ with $W = \operatorname{diag}(w_i)$, hence is an ordinary least-squares problem for the scaled matrix $W^{1/2}A$ and vector $W^{1/2}\mathbf{b}$. Write down its normal equations.
(b) Explain why, for integer weights, this is the same as duplicating row $i$ of the data $w_i$ times.
(c) Implement `least_squares_weighted(A, b, w)` in `project.py` by scaling rows and calling `least_squares_qr`. Assert that all weights equal to 1 reproduces the unweighted answer for Example 4.2, and that a weight of $10^6$ on the third point forces the fitted line through $(2, 2)$. Then fit Example 4.5 with weights $w = (1, 1, 1, 1, 1, 10)$ and report how the coefficients move.

## C. Challenge (do at least one)

**C1.** (QR is unique up to signs.) (a) Prove that the inverse of an invertible upper triangular matrix is upper triangular. (b) Prove that an orthogonal upper triangular matrix is diagonal with entries $\pm1$. (c) Complete the proof of Theorem 3.20 and count the number of distinct reduced QR factorizations of an $m\times n$ matrix with independent columns when no sign convention is imposed. (d) `np.linalg.qr(A)` on Example 4.2's $A$ returns $R$ with diagonal $(-2, -\sqrt5)$. Write the two-line fix that converts NumPy's $(Q, R)$ into the positive-diagonal factorization and explain why it is still a valid QR.

**C2.** (Legendre polynomials are Gram–Schmidt.) On the space of polynomials with inner product $\langle f, g\rangle = \int_{-1}^{1}f(x)g(x)\,dx$, run Gram–Schmidt on $1, x, x^2, x^3$ *without* normalizing (so that each output has leading coefficient 1). Show the outputs are $1,\ x,\ x^2 - \tfrac13,\ x^3 - \tfrac35x$. Rescale each so that its value at $x = 1$ is $1$ and compare with the Legendre polynomials $P_2 = \tfrac12(3x^2 - 1)$, $P_3 = \tfrac12(5x^3 - 3x)$. Explain, using parity, why half the inner products were zero before you computed them. Then use SymPy to get the degree-4 polynomial and check it against $P_4 = \tfrac18(35x^4 - 30x^2 + 3)$. (These are the polynomials that EM-3 gets from separating Laplace's equation in spherical coordinates.)

**C3.** (The log-decrement method.) A lightly damped oscillator rings down as $y(t) = A e^{-\gamma t}\cos\omega t$.
(a) Show that the local maxima occur where $\tan\omega t = -\gamma/\omega$, and that at every such maximum $\cos\omega t$ has the same value $\omega/\sqrt{\omega^2 + \gamma^2}$. Conclude that the peak heights satisfy $\ln y_k = \ln A' - \gamma t_k$ *exactly*, with $A' = A\omega/\sqrt{\omega^2 + \gamma^2}$, so that $\gamma$ is the slope of a least-squares line through $(t_k, \ln y_k)$.
(b) Write code (in the "Your turn" cells of `project.py`) that generates $y$ with $A = 2$, $\gamma = 0.15$, $\omega = 2\pi$ on $t \in [0, 20]$ (4001 samples) with Gaussian noise of standard deviation $0.01$, picks the largest sample in each period, fits the line, and reports $\hat\gamma$. Repeat with noise $0.1$. You will find $\hat\gamma$ biased *low* at both noise levels. Explain the bias. (Hint: near a peak the signal is flat over many samples; the largest of two hundred noisy values is not the peak. Is the error positive or negative, and which peaks does it hurt most in log space?)
(c) Two fixes. First, the period is known, so sample $y$ at $t_k = kT$ with $T = 2\pi/\omega$, where $\cos\omega t_k = 1$ exactly: then $y_k = Ae^{-\gamma t_k}$ plus *unbiased* noise. Discard any $y_k$ below $3\sigma$. Second, using B8, weight each remaining point by $w_k = y_k^2$ and explain why that is the right weight (linearize $\ln(y_k + \delta)$). Report $\hat\gamma$ at both noise levels with and without the weights.

**C4.** (Why modified beats classical.) (a) Prove by induction on $i$ that in exact arithmetic, the coefficient $r_{ij}$ computed by modified Gram–Schmidt (an inner product of $\mathbf{q}_i$ with the *current* vector $\mathbf{a}_j - \sum_{l<i}r_{lj}\mathbf{q}_l$) equals $\mathbf{q}_i^\top\mathbf{a}_j$, the classical coefficient.
(b) Now compute by hand in floating point, with $\varepsilon = 10^{-8}$ and the single rounding rule $\operatorname{fl}(1 + \varepsilon^2) = 1$, both versions on the columns of
$$A = \begin{pmatrix}1 & 1 & 1\\ \varepsilon & 0 & 0\\ 0 & \varepsilon & 0\\ 0 & 0 & \varepsilon\end{pmatrix}.$$
Show that classical Gram–Schmidt returns $\mathbf{q}_2^\top\mathbf{q}_3 = \tfrac12$ while modified returns $\mathbf{q}_2^\top\mathbf{q}_3 = 0$. Point to the exact line where classical goes wrong. Confirm numerically.
(c) "Your turn" item 5 says that running classical Gram–Schmidt again on its own output (reorthogonalization) recovers orthogonality. Explain why a repeated pass works when the first did not, use the $\kappa^2\varepsilon$ estimate of §3.8 with $\kappa(Q_1)$ in place of $\kappa(A)$ to predict how many passes the $10\times10$ Hilbert matrix needs, and confirm numerically.
