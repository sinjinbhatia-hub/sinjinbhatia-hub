# LA-7 Solutions

## A. Warm-up

**A1.** $T = 7$, $D = 12 - 2 = 10$. $p_A(\lambda) = \lambda^2 - 7\lambda + 10 = (\lambda - 2)(\lambda - 5)$, so $\lambda = 2, 5$.

$\lambda = 2$: $A - 2I = \begin{pmatrix}2 & 1\\ 2 & 1\end{pmatrix}$. Row 1: $2x_1 + x_2 = 0$, so $\mathbf{x}_1 = (1, -2)^\top$. Check: $A\mathbf{x}_1 = (4 - 2, 2 - 6)^\top = (2, -4)^\top = 2\mathbf{x}_1$. ✓

$\lambda = 5$: $A - 5I = \begin{pmatrix}-1 & 1\\ 2 & -2\end{pmatrix}$. Row 1: $x_1 = x_2$, so $\mathbf{x}_2 = (1, 1)^\top$. Check: $A\mathbf{x}_2 = (5, 5)^\top = 5\mathbf{x}_2$. ✓

Theorem 3.4: $2 + 5 = 7 = T$ ✓, $2\cdot5 = 10 = D$ ✓.

**A2.** $U$ is upper triangular, so by Theorem 3.5(a) its eigenvalues are the diagonal entries $2, -1, 5$. By Theorem 3.5(b): $U^2$ has $4, 1, 25$; $U^{-1}$ has $\tfrac12, -1, \tfrac15$; $U + 3I$ has $5, 2, 8$. By Theorem 3.5(c), $U^\top$ has $2, -1, 5$.

$U\mathbf{e}_1 = $ first column $= (2, 0, 0)^\top = 2\mathbf{e}_1$, so $\mathbf{e}_1$ is an eigenvector of $U$ (eigenvalue $2$), and by Theorem 3.5(b) also of $U^2$ (eigenvalue $4$), $U^{-1}$ ($\tfrac12$), and $U + 3I$ ($5$). But $U^\top\mathbf{e}_1 = $ first column of $U^\top$ $= (2, 7, -3)^\top$, not a multiple of $\mathbf{e}_1$. Same eigenvalues, different eigenvectors.

Eigenvector for $-1$: $U + I = \begin{pmatrix}3 & 7 & -3\\ 0 & 0 & 4\\ 0 & 0 & 6\end{pmatrix}$. Rows 2, 3 force $x_3 = 0$; row 1: $3x_1 + 7x_2 = 0$. Take $x_2 = -3$, $x_1 = 7$: $(7, -3, 0)^\top$. Check: $U(7, -3, 0)^\top = (14 - 21, 3, 0)^\top = (-7, 3, 0)^\top = -1\cdot(7, -3, 0)^\top$. ✓

Eigenvector for $5$: $U - 5I = \begin{pmatrix}-3 & 7 & -3\\ 0 & -6 & 4\\ 0 & 0 & 0\end{pmatrix}$. Row 2: $x_2 = \tfrac23x_3$; take $x_3 = 3$, $x_2 = 2$. Row 1: $-3x_1 + 14 - 9 = 0$, $x_1 = \tfrac53$. Scale by $3$: $(5, 6, 9)^\top$. Check: $U(5, 6, 9)^\top = (10 + 42 - 27, -6 + 36, 45)^\top = (25, 30, 45)^\top = 5(5, 6, 9)^\top$. ✓

**A3.** Entries nonnegative, columns $0.8 + 0.2 = 1$ and $0.3 + 0.7 = 1$: Markov. $T = 1.5$, $D = 0.56 - 0.06 = 0.5$; $\lambda^2 - 1.5\lambda + 0.5 = (\lambda - 1)(\lambda - 0.5)$, so $\lambda = 1, 0.5$.

Steady state: $A - I = \begin{pmatrix}-0.2 & 0.3\\ 0.2 & -0.3\end{pmatrix}$, row 1: $2x_1 = 3x_2$, so $\mathbf{x}_1 = (3, 2)^\top$ and $\mathbf{u}_\infty = (0.6, 0.4)^\top$. Other eigenvector: $A - 0.5I = \begin{pmatrix}0.3 & 0.3\\ 0.2 & 0.2\end{pmatrix}$, $\mathbf{x}_2 = (1, -1)^\top$.

Directly: $\mathbf{u}_1 = A(1, 0)^\top = (0.8, 0.2)^\top$ (first column). $\mathbf{u}_2 = A(0.8, 0.2)^\top = (0.64 + 0.06, 0.16 + 0.14)^\top = (0.7, 0.3)^\top$.

Closed form: $(1, 0)^\top = c_1(3, 2)^\top + c_2(1, -1)^\top$: $3c_1 + c_2 = 1$, $2c_1 - c_2 = 0$. Adding: $5c_1 = 1$, $c_1 = 0.2$, $c_2 = 0.4$. So
$$\mathbf{u}_k = 0.2(3, 2)^\top + 0.4(0.5)^k(1, -1)^\top = \big(0.6 + 0.4\cdot0.5^k,\ 0.4 - 0.4\cdot0.5^k\big)^\top.$$
$k = 1$: $(0.6 + 0.2, 0.4 - 0.2) = (0.8, 0.2)$ ✓. $k = 2$: $(0.6 + 0.1, 0.4 - 0.1) = (0.7, 0.3)$ ✓. The transient halves every step.

**A4.** $T = 0$, $D = 1$: $\lambda^2 + 1 = 0$, $\lambda = \pm i$. (Also (3.17) with $\theta = \pi/2$: $e^{\pm i\pi/2} = \pm i$.)

$\lambda = i$: $R - iI = \begin{pmatrix}-i & -1\\ 1 & -i\end{pmatrix}$. Row 1: $-ix_1 - x_2 = 0$, so $x_2 = -ix_1$: $\mathbf{x} = (1, -i)^\top$. Check: $R(1, -i)^\top = (0\cdot1 + (-1)(-i),\ 1\cdot1 + 0)^\top = (i, 1)^\top$, and $i(1, -i)^\top = (i, -i^2)^\top = (i, 1)^\top$. ✓

$\lambda = -i$: $\mathbf{x} = (1, i)^\top$. Check: $R(1, i)^\top = (-i, 1)^\top$ and $-i(1, i)^\top = (-i, -i^2)^\top = (-i, 1)^\top$. ✓

Reconciliation: Theorem 3.5(b) says the eigenvectors $(1, \mp i)$ of $R$ are eigenvectors of $R^2$ with eigenvalues $(\pm i)^2 = -1$. That is true. It does not say these are the *only* eigenvectors of $R^2$. Here $R^2 = -I$ has eigenvalue $-1$ with geometric multiplicity $2$ (the whole space), and $(1, -i)$, $(1, i)$ are two independent vectors in that eigenspace, as they should be. Two distinct eigenvalues of $R$ merged into one eigenvalue of $R^2$ because $i^2 = (-i)^2$.

**A5.** (a) $A$ is diagonal, so $A^k = \operatorname{diag}(2^k, (-3)^k)$ and the series (3.12) acts entrywise on the diagonal: $e^{At} = \operatorname{diag}\big(\sum_k(2t)^k/k!,\ \sum_k(-3t)^k/k!\big) = \operatorname{diag}(e^{2t}, e^{-3t})$.

(b) $N^2 = \begin{pmatrix}0 & 1\\ 0 & 0\end{pmatrix}\begin{pmatrix}0 & 1\\ 0 & 0\end{pmatrix} = 0$, so all higher powers vanish and the series stops: $e^{Nt} = I + Nt = \begin{pmatrix}1 & t\\ 0 & 1\end{pmatrix}$.

(c) $A = 2I + N$. $2I$ commutes with every matrix, so Theorem 3.15(c) gives $e^{At} = e^{2It}e^{Nt} = e^{2t}I\cdot(I + Nt) = e^{2t}\begin{pmatrix}1 & t\\ 0 & 1\end{pmatrix}$. Derivative: $\frac{d}{dt}\big[e^{2t}(I + Nt)\big] = 2e^{2t}(I + Nt) + e^{2t}N$; at $t = 0$ this is $2I + N = A$. ✓ Note $A$ is defective ($\lambda = 2$ twice, $A - 2I = N$ has rank 1) and the $te^{2t}$ entry is the signature of that.

## B. Core

**B1.** (a) Suppose $A\mathbf{x} = \lambda\mathbf{x}$, $\mathbf{x}\ne\mathbf{0}$, $A$ invertible. Apply $A^{-1}$: $\mathbf{x} = \lambda A^{-1}\mathbf{x}$. If $\lambda = 0$ this says $\mathbf{x} = \mathbf{0}$, a contradiction, so $\lambda\ne0$ and $A^{-1}\mathbf{x} = \lambda^{-1}\mathbf{x}$. Conversely, if $A^{-1}\mathbf{x} = \mu\mathbf{x}$, the same argument with $A^{-1}$ in place of $A$ (its inverse is $A$) gives $\mu \ne 0$ and $A\mathbf{x} = \mu^{-1}\mathbf{x}$; with $\mu = 1/\lambda$ this is $A\mathbf{x} = \lambda\mathbf{x}$. So the eigenvectors coincide and the eigenvalues are reciprocals. ∎

$0$ is an eigenvalue $\iff$ there is $\mathbf{x}\ne\mathbf{0}$ with $A\mathbf{x} = \mathbf{0}$ $\iff$ $N(A)\ne\{\mathbf{0}\}$ $\iff$ $A$ is singular (LA-2/LA-3). Equivalently $p_A(0) = \det A = 0$. ∎

(b) $A$ is triangular with diagonal $0, 0$: eigenvalues $0, 0$. $B$ likewise: $0, 0$. $A + B = \begin{pmatrix}0 & 1\\ 1 & 0\end{pmatrix}$: $T = 0$, $D = -1$, $\lambda^2 - 1 = 0$, eigenvalues $\pm1$. $AB = \begin{pmatrix}1 & 0\\ 0 & 0\end{pmatrix}$: eigenvalues $1, 0$. Sums of eigenvalues of $A$ and $B$ can only be $0$; products can only be $0$. Neither $\pm1$ nor $1$ is obtained. ∎

(c) $A + B = X\Lambda_AX^{-1} + X\Lambda_BX^{-1} = X(\Lambda_A + \Lambda_B)X^{-1}$. This is a diagonalization of $A + B$, so by Theorem 3.8 its eigenvalues are the diagonal entries $\lambda_i + \mu_i$ (same eigenvectors, the columns of $X$). $AB = X\Lambda_AX^{-1}X\Lambda_BX^{-1} = X\Lambda_A\Lambda_BX^{-1}$, diagonal entries $\lambda_i\mu_i$. Finally $BA = X\Lambda_B\Lambda_AX^{-1} = X\Lambda_A\Lambda_BX^{-1} = AB$ because diagonal matrices commute. ∎ (The converse, that commuting diagonalizable matrices are simultaneously diagonalizable, is also true and is the mathematics behind "compatible observables" in QM-4.)

**B2.** $T = 2$, $D = 1 - 4 = -3 < 0$: a saddle, whatever $T$ is. Eigenvalues: $\lambda^2 - 2\lambda - 3 = (\lambda - 3)(\lambda + 1)$, so $3$ and $-1$.

$\lambda = 3$: $A - 3I = \begin{pmatrix}-2 & 2\\ 2 & -2\end{pmatrix}$, $\mathbf{x}_1 = (1, 1)^\top$. $\lambda = -1$: $A + I = \begin{pmatrix}2 & 2\\ 2 & 2\end{pmatrix}$, $\mathbf{x}_2 = (1, -1)^\top$.

$(3, 1)^\top = c_1(1, 1)^\top + c_2(1, -1)^\top$: $c_1 + c_2 = 3$, $c_1 - c_2 = 1$, so $c_1 = 2$, $c_2 = 1$. By (3.11):
$$\mathbf{u}(t) = 2e^{3t}\begin{pmatrix}1\\1\end{pmatrix} + e^{-t}\begin{pmatrix}1\\-1\end{pmatrix} = \begin{pmatrix}2e^{3t} + e^{-t}\\ 2e^{3t} - e^{-t}\end{pmatrix}.$$
Check: $\mathbf{u}(0) = (3, 1)^\top$ ✓. $\dot{\mathbf{u}}(0) = (6 - 1, 6 + 1)^\top = (5, 7)^\top$; $A\mathbf{u}(0) = (3 + 2, 6 + 1)^\top = (5, 7)^\top$ ✓.

As $t\to+\infty$ the $e^{3t}$ term dominates: the trajectory runs off along the line spanned by $(1, 1)$ (the unstable eigenline). As $t\to-\infty$ the $e^{-t}$ term dominates: the trajectory came in from far out along $(1, -1)$ (the stable eigenline). The eigenvectors are orthogonal, $(1, 1)\cdot(1, -1) = 0$, because $A$ is symmetric and the eigenvalues are distinct (Theorem 3.17(b)).

**B3.** (a) $T = 4$, $D = 3$, $\lambda^2 - 4\lambda + 3 = (\lambda - 1)(\lambda - 3)$. $\lambda = 3$: $A - 3I = \begin{pmatrix}-1 & 1\\ 1 & -1\end{pmatrix}$, $\mathbf{x}_1 = (1, 1)^\top$. $\lambda = 1$: $A - I = \begin{pmatrix}1 & 1\\ 1 & 1\end{pmatrix}$, $\mathbf{x}_2 = (1, -1)^\top$.
$$X = \begin{pmatrix}1 & 1\\ 1 & -1\end{pmatrix},\quad \Lambda = \begin{pmatrix}3 & 0\\ 0 & 1\end{pmatrix},\quad X^{-1} = \frac{1}{-2}\begin{pmatrix}-1 & -1\\ -1 & 1\end{pmatrix} = \frac12\begin{pmatrix}1 & 1\\ 1 & -1\end{pmatrix}.$$
(Check $XX^{-1} = \frac12\begin{pmatrix}2 & 0\\ 0 & 2\end{pmatrix} = I$ ✓.)

(b) By (3.7), $A^k = X\Lambda^kX^{-1} = \begin{pmatrix}1 & 1\\ 1 & -1\end{pmatrix}\begin{pmatrix}3^k & 0\\ 0 & 1\end{pmatrix}\cdot\frac12\begin{pmatrix}1 & 1\\ 1 & -1\end{pmatrix}$. First $X\Lambda^k = \begin{pmatrix}3^k & 1\\ 3^k & -1\end{pmatrix}$. Then
$$A^k = \frac12\begin{pmatrix}3^k & 1\\ 3^k & -1\end{pmatrix}\begin{pmatrix}1 & 1\\ 1 & -1\end{pmatrix} = \frac12\begin{pmatrix}3^k + 1 & 3^k - 1\\ 3^k - 1 & 3^k + 1\end{pmatrix}.$$
$k = 5$: $3^5 = 243$, so $A^5 = \frac12\begin{pmatrix}244 & 242\\ 242 & 244\end{pmatrix} = \begin{pmatrix}122 & 121\\ 121 & 122\end{pmatrix}$. Direct: $A^2 = \begin{pmatrix}5 & 4\\ 4 & 5\end{pmatrix}$, $A^4 = (A^2)^2 = \begin{pmatrix}41 & 40\\ 40 & 41\end{pmatrix}$, $A^5 = A^4A = \begin{pmatrix}82 + 40 & 41 + 80\\ 80 + 41 & 40 + 82\end{pmatrix} = \begin{pmatrix}122 & 121\\ 121 & 122\end{pmatrix}$. ✓

(c) $A^k/3^k = \frac12\begin{pmatrix}1 + 3^{-k} & 1 - 3^{-k}\\ 1 - 3^{-k} & 1 + 3^{-k}\end{pmatrix}\to P = \frac12\begin{pmatrix}1 & 1\\ 1 & 1\end{pmatrix}$. $P^2 = \frac14\begin{pmatrix}2 & 2\\ 2 & 2\end{pmatrix} = P$ ✓, $P^\top = P$ ✓. $P = \mathbf{x}_1\mathbf{x}_1^\top/\|\mathbf{x}_1\|^2$ with $\mathbf{x}_1 = (1, 1)^\top$: it is the orthogonal projection onto the line through $(1, 1)$, exactly LA-1's $t^*\mathbf{v} = (\mathbf{v}\cdot\mathbf{u}/\|\mathbf{v}\|^2)\mathbf{v}$ written as a matrix. Power iteration: $A^k\mathbf{u}_0/3^k \to P\mathbf{u}_0 = c_1\mathbf{x}_1$, the component of $\mathbf{u}_0$ along the dominant eigenvector, which is step 1 of the derivation in §3.6 with $\lambda_1 = 3$. The other eigen-component is killed by $(1/3)^k$.

**B4.** (a) Entries: $(A^k)_{ij} = \sum_{\ell}(A^{k-1})_{i\ell}a_{\ell j}$ is a sum of products of nonnegative numbers, so nonnegative (induction on $k$). Column sums: $\mathbf{1}^\top A^k = (\mathbf{1}^\top A)A^{k-1} = \mathbf{1}^\top A^{k-1} = \cdots = \mathbf{1}^\top$. So $A^k$ is Markov. ∎

(b) $\mathbf{1}^\top(A\mathbf{x}) = (\mathbf{1}^\top A)\mathbf{x} = \mathbf{1}^\top\mathbf{x}$. Also $\mathbf{1}^\top(A\mathbf{x}) = \mathbf{1}^\top(\lambda\mathbf{x}) = \lambda\mathbf{1}^\top\mathbf{x}$. So $(1 - \lambda)\mathbf{1}^\top\mathbf{x} = 0$, and $\lambda\ne1$ gives $\mathbf{1}^\top\mathbf{x} = 0$. ∎

(c) Apply $\mathbf{1}^\top$ to $\mathbf{u}_k = c_1\mathbf{x}_1 + \sum_{i\ge2}c_i\lambda_i^k\mathbf{x}_i$. Left side: $\mathbf{1}^\top\mathbf{u}_k = 1$ (probability vector, by (a)). Right side: $c_1\mathbf{1}^\top\mathbf{x}_1 + 0$, because (b) applies to every $\mathbf{x}_i$ with $i\ge2$. This is where simplicity is used: (b) needs $\lambda_i \ne 1$. So $c_1\mathbf{1}^\top\mathbf{x}_1 = 1$; in particular $\mathbf{1}^\top\mathbf{x}_1 \ne 0$, and $c_1 = 1/(\mathbf{1}^\top\mathbf{x}_1)$ depends only on the chosen scaling of $\mathbf{x}_1$, not on $\mathbf{u}_0$. The transient terms $c_i\lambda_i^k\mathbf{x}_i$ each have entries summing to zero: they move probability *between* states (some entries positive, some negative) without changing the total. Every starting distribution has the same steady-state component; only the transients differ.

Without simplicity this fails. $A = I$ is Markov, with $\lambda_1 = \lambda_2 = 1$ and eigenvectors $\mathbf{x}_1 = \mathbf{e}_1$, $\mathbf{x}_2 = \mathbf{e}_2$. Then $\mathbf{u}_0 = (a, 1 - a)^\top = a\,\mathbf{e}_1 + (1 - a)\,\mathbf{e}_2$ gives $c_1 = a$, which depends on $\mathbf{u}_0$. The step that breaks is (b): $\mathbf{x}_2$ has eigenvalue $1$, so nothing forces its entries to sum to zero (they sum to $1$). Physically, $A = I$ never moves any probability, so the chain remembers its start forever.

(d) $X = \begin{pmatrix}5 & 1\\ 1 & -1\end{pmatrix}$ (columns $\mathbf{x}_1 = (5, 1)^\top$ for $\lambda = 1$, $\mathbf{x}_2 = (1, -1)^\top$ for $0.4$), $\Lambda = \operatorname{diag}(1, 0.4)$, $\det X = -6$, $X^{-1} = \frac{1}{-6}\begin{pmatrix}-1 & -1\\ -1 & 5\end{pmatrix} = \frac16\begin{pmatrix}1 & 1\\ 1 & -5\end{pmatrix}$.
$$A^k = X\Lambda^kX^{-1} = \frac16\begin{pmatrix}5 & 0.4^k\\ 1 & -0.4^k\end{pmatrix}\begin{pmatrix}1 & 1\\ 1 & -5\end{pmatrix} = \frac16\begin{pmatrix}5 + 0.4^k & 5 - 5\cdot0.4^k\\ 1 - 0.4^k & 1 + 5\cdot0.4^k\end{pmatrix}.$$
Check $k = 1$: $\frac16\begin{pmatrix}5.4 & 3\\ 0.6 & 3\end{pmatrix} = \begin{pmatrix}0.9 & 0.5\\ 0.1 & 0.5\end{pmatrix}$ ✓. Limit: $A^\infty = \frac16\begin{pmatrix}5 & 5\\ 1 & 1\end{pmatrix} = \begin{pmatrix}\mathbf{u}_\infty & \mathbf{u}_\infty\end{pmatrix}$. Column $j$ of $A^k$ is $A^k\mathbf{e}_j$, the distribution after $k$ days given that day 0 was in state $j$ with certainty. Both columns converge to the same steady state: where you start stops mattering. And $A^\infty\mathbf{u}_0 = \mathbf{u}_\infty(\mathbf{1}^\top\mathbf{u}_0) = \mathbf{u}_\infty$ for every probability vector.

**B5.** (a) $T = 2$, $D = 1 + 1 = 2$, $\Delta = 4 - 8 = -4 < 0$, $T > 0$: unstable spiral. $\lambda = \frac{2\pm\sqrt{-4}}{2} = 1\pm i$. For $\lambda = 1 + i$: $A - \lambda I = \begin{pmatrix}-i & -1\\ 1 & -i\end{pmatrix}$; row 1: $-ix_1 - x_2 = 0$, $x_2 = -ix_1$: $\mathbf{x} = (1, -i)^\top$, so $\mathbf{r} = (1, 0)^\top$, $\mathbf{s} = (0, -1)^\top$. (Row 2 check: $1 + (-i)(-i) = 1 - 1 = 0$ ✓.)

(b) With $a = 1$, $b = 1$ in (3.19):
$\mathbf{p}(t) = e^t(\mathbf{r}\cos t - \mathbf{s}\sin t) = e^t(\cos t, \sin t)^\top$,
$\mathbf{q}(t) = e^t(\mathbf{r}\sin t + \mathbf{s}\cos t) = e^t(\sin t, -\cos t)^\top$.
$\mathbf{p}(0) = (1, 0)^\top = \mathbf{u}(0)$, so $\mathbf{u}(t) = e^t(\cos t, \sin t)^\top$. Check: $\dot{\mathbf{u}}(0) = (1, 1)^\top$ and $A(1, 0)^\top = (1, 1)^\top$ ✓.

(c) $A = I + J$ with $a = 1$, $b = 1$. By (3.20), $e^{At} = e^tR(t)$, and $e^{At}(1, 0)^\top = e^t(\cos t, \sin t)^\top$. Same. Rotation is counterclockwise ($b = a_{21} = 1 > 0$: at $(1, 0)$ the velocity $(1, 1)$ points up), and the radius $e^t$ grows: an outward counterclockwise spiral.

(d) $B(\mathbf{r} + i\mathbf{s}) = (a + ib)(\mathbf{r} + i\mathbf{s}) = (a\mathbf{r} - b\mathbf{s}) + i(b\mathbf{r} + a\mathbf{s})$. Real parts: $B\mathbf{r} = a\mathbf{r} - b\mathbf{s}$. Imaginary parts: $B\mathbf{s} = b\mathbf{r} + a\mathbf{s}$. Put these as columns: $B\begin{pmatrix}\mathbf{r} & \mathbf{s}\end{pmatrix} = \begin{pmatrix}a\mathbf{r} - b\mathbf{s} & b\mathbf{r} + a\mathbf{s}\end{pmatrix} = \begin{pmatrix}\mathbf{r} & \mathbf{s}\end{pmatrix}\begin{pmatrix}a & b\\ -b & a\end{pmatrix}$. With $P = \begin{pmatrix}\mathbf{r} & \mathbf{s}\end{pmatrix}$ invertible (step 4 of the derivation of (3.19): $\mathbf{r}, \mathbf{s}$ independent), $P^{-1}BP = \begin{pmatrix}a & b\\ -b & a\end{pmatrix}$, which is $aI - bJ$: rotation (clockwise by $bt$, in these coordinates) with scaling $e^{at}$. So every real $2\times2$ matrix with complex eigenvalues is similar to a rotation-with-scaling, and its phase portrait is a spiral or ellipse in skewed coordinates. ∎

**B6.** (a) $\omega_1^2\,\mathbf{v}_1^\top M\mathbf{v}_2 = (\omega_1^2M\mathbf{v}_1)^\top\mathbf{v}_2$ (using $M^\top = M$ to write $\mathbf{v}_1^\top M = (M\mathbf{v}_1)^\top$) $= (K\mathbf{v}_1)^\top\mathbf{v}_2 = \mathbf{v}_1^\top K^\top\mathbf{v}_2 = \mathbf{v}_1^\top K\mathbf{v}_2 = \mathbf{v}_1^\top(\omega_2^2M\mathbf{v}_2) = \omega_2^2\,\mathbf{v}_1^\top M\mathbf{v}_2$. So $(\omega_1^2 - \omega_2^2)\mathbf{v}_1^\top M\mathbf{v}_2 = 0$, and $\omega_1^2\ne\omega_2^2$ forces $\mathbf{v}_1^\top M\mathbf{v}_2 = 0$. ∎

(b) First the lemma. For $\mathbf{v} = \mathbf{r} + i\mathbf{s}$ with $\mathbf{r}, \mathbf{s}$ real, not both zero: $\bar{\mathbf{v}}^\top M\mathbf{v} = (\mathbf{r} - i\mathbf{s})^\top M(\mathbf{r} + i\mathbf{s}) = \mathbf{r}^\top M\mathbf{r} + i\,\mathbf{r}^\top M\mathbf{s} - i\,\mathbf{s}^\top M\mathbf{r} + \mathbf{s}^\top M\mathbf{s}$. Since $M$ is symmetric, $\mathbf{s}^\top M\mathbf{r} = (\mathbf{s}^\top M\mathbf{r})^\top = \mathbf{r}^\top M\mathbf{s}$, so the imaginary terms cancel and $\bar{\mathbf{v}}^\top M\mathbf{v} = \mathbf{r}^\top M\mathbf{r} + \mathbf{s}^\top M\mathbf{s} > 0$, because at least one of $\mathbf{r}, \mathbf{s}$ is nonzero and $M$ is positive definite on real vectors.

Now let $K\mathbf{v} = \omega^2M\mathbf{v}$, $\mathbf{v}\ne\mathbf{0}$. Multiply on the left by $\bar{\mathbf{v}}^\top$: $\bar{\mathbf{v}}^\top K\mathbf{v} = \omega^2\,\bar{\mathbf{v}}^\top M\mathbf{v}$. Conjugate the original equation ($K$, $M$ real): $K\bar{\mathbf{v}} = \overline{\omega^2}M\bar{\mathbf{v}}$; transpose ($K$, $M$ symmetric): $\bar{\mathbf{v}}^\top K = \overline{\omega^2}\,\bar{\mathbf{v}}^\top M$; multiply on the right by $\mathbf{v}$: $\bar{\mathbf{v}}^\top K\mathbf{v} = \overline{\omega^2}\,\bar{\mathbf{v}}^\top M\mathbf{v}$. Subtract: $(\omega^2 - \overline{\omega^2})\bar{\mathbf{v}}^\top M\mathbf{v} = 0$, and the second factor is positive, so $\omega^2 = \overline{\omega^2}$ is real. ∎ (This is Theorem 3.17(a) with $M$ playing the role of $I$.)

Positivity comes with it. The first equation gives $\omega^2 = \bar{\mathbf{v}}^\top K\mathbf{v}\,/\,\bar{\mathbf{v}}^\top M\mathbf{v}$. If $K$ is also positive definite (as for a spring system near a stable equilibrium), the lemma above applied to $K$ instead of $M$ makes the numerator positive, so $\omega^2 > 0$ and $\omega$ is a real frequency. This is the positivity claimed in notes §3.10.

(c) From Example 4.3, $\mathbf{v}_1 = (1, \tfrac{1+\sqrt3}{2})^\top$, $\mathbf{v}_2 = (1, \tfrac{1-\sqrt3}{2})^\top$, $M = \operatorname{diag}(1, 2)$. $\mathbf{v}_1^\top M\mathbf{v}_2 = 1\cdot1\cdot1 + 2\cdot\tfrac{(1+\sqrt3)(1-\sqrt3)}{4} = 1 + 2\cdot\tfrac{-2}{4} = 1 - 1 = 0$. ✓ Ordinary: $\mathbf{v}_1^\top\mathbf{v}_2 = 1 + \tfrac{-2}{4} = \tfrac12 \ne 0$.

(d) $M^{-1}K = \begin{pmatrix}2 & -1\\ -\tfrac12 & 1\end{pmatrix}$ is not symmetric, and the proof of Theorem 3.17(b) uses $A^\top = A$ in the step $\mathbf{x}_1^\top A^\top\mathbf{x}_2 = \mathbf{x}_1^\top A\mathbf{x}_2$.

(e) Substitute $\mathbf{v} = M^{-1/2}\mathbf{w}$: $KM^{-1/2}\mathbf{w} = \omega^2MM^{-1/2}\mathbf{w} = \omega^2M^{1/2}\mathbf{w}$. Multiply on the left by $M^{-1/2}$: $(M^{-1/2}KM^{-1/2})\mathbf{w} = \omega^2\mathbf{w}$. The matrix $\tilde K = M^{-1/2}KM^{-1/2}$ is symmetric: $\tilde K^\top = (M^{-1/2})^\top K^\top(M^{-1/2})^\top = M^{-1/2}KM^{-1/2}$ since $M^{-1/2}$ is diagonal. So $\omega^2$ are eigenvalues of a symmetric matrix (real, orthogonal eigenvectors $\mathbf{w}_j$ by Theorem 3.17), and $\mathbf{w}_1^\top\mathbf{w}_2 = (M^{1/2}\mathbf{v}_1)^\top(M^{1/2}\mathbf{v}_2) = \mathbf{v}_1^\top M\mathbf{v}_2$: ordinary orthogonality of the $\mathbf{w}$'s is $M$-orthogonality of the $\mathbf{v}$'s.

Example 4.3: $M^{-1/2} = \operatorname{diag}(1, 1/\sqrt2)$, so $\tilde K = \begin{pmatrix}2 & -1/\sqrt2\\ -1/\sqrt2 & 1\end{pmatrix}$. $\operatorname{tr}\tilde K = 3 = \omega_1^2 + \omega_2^2$ ✓; $\det\tilde K = 2 - \tfrac12 = \tfrac32 = \omega_1^2\omega_2^2 = \tfrac{(3-\sqrt3)(3+\sqrt3)}{4} = \tfrac64$ ✓. (This substitution is how `scipy.linalg.eigh(K, M)` solves the generalized problem, with a Cholesky factor of $M$ in place of $M^{1/2}$ for non-diagonal $M$; the Cholesky factorization is LA-10 §3.2.)

**B7.** (a) Both are upper triangular with diagonal $3, 3$: $p_A = p_B = (3 - \lambda)^2$. $A - 3I = \begin{pmatrix}0 & 1\\ 0 & 0\end{pmatrix}$ has rank $1$, so $\operatorname{GM}(3) = 2 - 1 = 1 < 2 = \operatorname{AM}(3)$: not diagonalizable (Theorem 3.20). If $A$ were similar to $B = 3I$, then $A = M(3I)M^{-1} = 3MM^{-1} = 3I$, false. So same characteristic polynomial, not similar.

(b) $A = 3I + N$ with $N = \begin{pmatrix}0 & 1\\ 0 & 0\end{pmatrix}$, $N^2 = 0$, and $3I$ commutes with $N$, so the binomial theorem applies: $A^k = \sum_j\binom kj(3I)^{k-j}N^j = 3^kI + k\,3^{k-1}N = \begin{pmatrix}3^k & k\,3^{k-1}\\ 0 & 3^k\end{pmatrix}$ (all terms with $j\ge2$ vanish). Check $k = 2$: $A^2 = \begin{pmatrix}9 & 6\\ 0 & 9\end{pmatrix}$, and directly $\begin{pmatrix}3 & 1\\ 0 & 3\end{pmatrix}^2 = \begin{pmatrix}9 & 3 + 3\\ 0 & 9\end{pmatrix}$ ✓. By Theorem 3.15(c), $e^{At} = e^{3t}e^{Nt} = e^{3t}(I + Nt) = e^{3t}\begin{pmatrix}1 & t\\ 0 & 1\end{pmatrix}$.

(c) $\mathbf{u}(t) = e^{At}\mathbf{u}_0 = e^{3t}\begin{pmatrix}1 & t\\ 0 & 1\end{pmatrix}\begin{pmatrix}1\\2\end{pmatrix} = e^{3t}\begin{pmatrix}1 + 2t\\ 2\end{pmatrix}$. Verify: $\dot{\mathbf{u}} = 3e^{3t}(1 + 2t, 2)^\top + e^{3t}(2, 0)^\top = e^{3t}(5 + 6t, 6)^\top$; $A\mathbf{u} = e^{3t}(3(1 + 2t) + 2, 6)^\top = e^{3t}(5 + 6t, 6)^\top$ ✓. $\mathbf{u}(0) = (1, 2)^\top$ ✓. Note the $te^{3t}$: with only one eigenvector, the second component of the initial condition cannot be matched by pure exponentials.

(d) $A_\varepsilon$ is triangular: eigenvalues $3$ and $3 + \varepsilon$, distinct, so diagonalizable. Eigenvectors: $A_\varepsilon - 3I = \begin{pmatrix}0 & 1\\ 0 & \varepsilon\end{pmatrix}$ gives $\mathbf{x}_1 = (1, 0)^\top$; $A_\varepsilon - (3 + \varepsilon)I = \begin{pmatrix}-\varepsilon & 1\\ 0 & 0\end{pmatrix}$ gives $x_2 = \varepsilon x_1$, $\mathbf{x}_2 = (1, \varepsilon)^\top$. Angle: $\cos\theta = \frac{1}{\sqrt{1 + \varepsilon^2}}$, so $\theta \approx \varepsilon$ radians $\to 0$: the eigenvectors become parallel. $X_\varepsilon = \begin{pmatrix}1 & 1\\ 0 & \varepsilon\end{pmatrix}$, $\det X_\varepsilon = \varepsilon$, $X_\varepsilon^{-1} = \frac1\varepsilon\begin{pmatrix}\varepsilon & -1\\ 0 & 1\end{pmatrix} = \begin{pmatrix}1 & -1/\varepsilon\\ 0 & 1/\varepsilon\end{pmatrix}$, with entries of size $1/\varepsilon\to\infty$.

Then $A_\varepsilon^k = X_\varepsilon\Lambda^kX_\varepsilon^{-1} = \begin{pmatrix}3^k & \big((3+\varepsilon)^k - 3^k\big)/\varepsilon\\ 0 & (3+\varepsilon)^k\end{pmatrix}$. The top-right entry is a difference quotient: two nearly equal large numbers, subtracted, divided by a small number. In floating point this is catastrophic cancellation. The two powers agree in about their first $\log_{10}(1/\varepsilon)$ digits, and the subtraction throws those away. With $\varepsilon = 10^{-8}$ roughly half of the 16 digits are wrong: the project measures relative errors of about $3\times10^{-9}$, $4\times10^{-10}$, $5\times10^{-9}$ for $k = 5, 10, 30$, against about $10^{-16}$ for `matrix_power`. As $\varepsilon$ approaches machine epsilon ($2\times10^{-16}$), no digits survive. The true value is close to $k\,3^{k-1}$ (the derivative of $x^k$ at $3$). The diagonalization exists but is numerically poor because $X_\varepsilon$ is nearly singular. The exponential has the same issue. This is why LA-13 computes with the Schur form, not with $X\Lambda X^{-1}$.

**B8.** Reference implementation, identical to `convergence_rate` in `project.py`, Part 1 (it returns the two error histories, their fitted slopes, and the predicted slope):

```python
def convergence_rate(A, v0, k_max=60, floor=1e-10, shift=None):
    A = np.asarray(A, float); n = A.shape[0]
    lam, X = np.linalg.eig(A)
    if shift is None:                                  # power iteration
        order = np.argsort(-np.abs(lam), kind="stable")
        step = lambda v: A @ v
    else:                                              # (c): inverse iteration
        order = np.argsort(np.abs(lam - shift), kind="stable")
        lu = lu_factor(A - shift * np.eye(n))
        step = lambda v: lu_solve(lu, v)
    lam, X = lam[order], X[:, order]
    x1 = X[:, 0].real; x1 = x1 / np.linalg.norm(x1); lam1 = lam[0].real
    v = np.asarray(v0, float); v = v / np.linalg.norm(v)
    dir_err, ray_err = [], []
    for _ in range(k_max):
        s = np.sign(v @ x1) or 1.0                     # align sign before comparing
        dir_err.append(np.linalg.norm(s * v - x1))
        ray_err.append(abs(v @ A @ v - lam1))
        w = step(v); v = w / np.linalg.norm(w)
    dir_err, ray_err = np.array(dir_err), np.array(ray_err)
    def slope(err):
        k = np.arange(len(err)); ok = (err > floor) & (err < 0.5)
        return np.polyfit(k[ok], np.log(err[ok]), 1)[0]
    if shift is None:
        pred = np.log(abs(lam[1] / lam[0]))
    else:
        pred = np.log(abs(lam[0] - shift) / abs(lam[1] - shift))
    return dir_err, ray_err, slope(dir_err), slope(ray_err), pred
```

(a) For $A = \begin{pmatrix}4 & 1\\ 2 & 3\end{pmatrix}$ the fitted direction slope is $\approx -0.917$, against $\log(2/5) = -0.916$, matching (3.9). For $\begin{pmatrix}2 & 1\\ 1 & 2\end{pmatrix}$ it is $\approx -1.098$, against $\log(1/3) = -1.099$. The project asserts both within a few percent.

(b) For the symmetric matrix the Rayleigh slope is $\approx -2.19$, against $2\log(1/3) = -2.197$. Why: let $\mathbf{x}_1, \mathbf{x}_2$ be orthonormal eigenvectors (possible since $A$ is symmetric, Theorem 3.17) and write the iterate as $\mathbf{v} = \mathbf{x}_1 + \varepsilon\mathbf{x}_2$ (not yet normalized). Then
$$\mathbf{v}^\top A\mathbf{v} = (\mathbf{x}_1 + \varepsilon\mathbf{x}_2)^\top(\lambda_1\mathbf{x}_1 + \varepsilon\lambda_2\mathbf{x}_2) = \lambda_1 + \varepsilon^2\lambda_2, \qquad \mathbf{v}^\top\mathbf{v} = 1 + \varepsilon^2,$$
because the cross terms $\mathbf{x}_1^\top\mathbf{x}_2$ vanish. So the Rayleigh quotient is $\frac{\lambda_1 + \varepsilon^2\lambda_2}{1 + \varepsilon^2} = \lambda_1 + \varepsilon^2(\lambda_2 - \lambda_1) + O(\varepsilon^4)$: the error is **second order** in $\varepsilon$. Since $\varepsilon_k \propto (\lambda_2/\lambda_1)^k$, the eigenvalue error goes like $(\lambda_2/\lambda_1)^{2k}$, twice the slope.

In the nonsymmetric case take $\mathbf{x}_1, \mathbf{x}_2$ of unit length with $c = \mathbf{x}_1^\top\mathbf{x}_2 \ne 0$. The cross terms no longer vanish, in the numerator *or* the denominator:
$$\mathbf{v}^\top A\mathbf{v} = \lambda_1 + \varepsilon(\lambda_1 + \lambda_2)c + \varepsilon^2\lambda_2, \qquad \mathbf{v}^\top\mathbf{v} = 1 + 2\varepsilon c + \varepsilon^2.$$
Dividing, with $1/(1 + 2\varepsilon c + \varepsilon^2) = 1 - 2\varepsilon c + O(\varepsilon^2)$:
$$\frac{\mathbf{v}^\top A\mathbf{v}}{\mathbf{v}^\top\mathbf{v}} = \lambda_1 + \varepsilon(\lambda_1 + \lambda_2)c - 2\varepsilon c\lambda_1 + O(\varepsilon^2) = \lambda_1 + \varepsilon\,c\,(\lambda_2 - \lambda_1) + O(\varepsilon^2).$$
The first-order term $\varepsilon c(\lambda_2 - \lambda_1)$ is nonzero whenever $c \ne 0$, so the eigenvalue error has the same slope as the direction error. The project measures $\approx -0.909$ for $\begin{pmatrix}4 & 1\\ 2 & 3\end{pmatrix}$, against $\log(2/5) = -0.916$, and asserts it.

(c) By Theorem 3.5(b), $(A - \sigma I)^{-1}$ has eigenvalues $1/(\lambda_i - \sigma)$ with the same eigenvectors. For $\sigma = 2.2$ these are $1/(2 - 2.2) = -5$ and $1/(5 - 2.2) \approx 0.357$. Power iteration on $(A - \sigma I)^{-1}$ therefore converges to the eigenvector for $\lambda = 2$, namely $(1, -2)^\top$, at ratio $0.357/5 = |2 - 2.2|/|5 - 2.2| = 0.2/2.8 = 1/14$. Predicted slope: $\log(1/14) = -2.639$; the project measures $\approx -2.637$ and asserts agreement within $5\%$. Each step gains $\log_{10}14 \approx 1.15$ decimal digits, against $\log_{10}(5/2) \approx 0.40$ for power iteration on the same matrix: about three times as fast, and a shift closer to $2$ would be faster still. The sign of $\mathbf{v}_k$ alternates (the dominant eigenvalue $-5$ is negative), which is why the code aligns signs before measuring.

## C. Challenge

**C1.** Let $A = X\Lambda X^{-1}$ and $p_A(\lambda) = \sum_jc_j\lambda^j$. By (3.7), $A^j = X\Lambda^jX^{-1}$, so
$$p_A(A) = \sum_jc_jX\Lambda^jX^{-1} = X\Big(\sum_jc_j\Lambda^j\Big)X^{-1} = X\,p_A(\Lambda)\,X^{-1}.$$
$p_A(\Lambda)$ is diagonal with $i$-th entry $\sum_jc_j\lambda_i^j = p_A(\lambda_i) = 0$, since each $\lambda_i$ is a root of $p_A$. So $p_A(\Lambda) = 0$ and $p_A(A) = 0$. ∎

(a) $p_A(\lambda) = \lambda^2 - T\lambda + D$, so $A^2 - TA + DI = 0$, i.e. $A(A - TI) = -DI$, i.e. $A\cdot\frac{TI - A}{D} = I$ when $D\ne0$. Hence $A^{-1} = (TI - A)/D$. For $A = \begin{pmatrix}a & b\\ c & d\end{pmatrix}$: $TI - A = \begin{pmatrix}d & -b\\ -c & a\end{pmatrix}$, and dividing by $D = ad - bc$ gives the LA-2 formula. ✓ (This holds for every $2\times2$ matrix, diagonalizable or not; the general proof is LA-11 Theorem 3.16 (Cayley–Hamilton), or via LA-6's adjugate identity.)

(b) Add to `project.py`:
```python
A = rng.standard_normal((4, 4))
coeffs = np.poly(A)                      # monic: lambda^4 + c3 lambda^3 + ... + c0
P = sum(c * np.linalg.matrix_power(A, 4 - j) for j, c in enumerate(coeffs))
assert np.allclose(P, 0, atol=1e-9 * np.abs(coeffs).max())
```
A random real $4\times4$ is diagonalizable with probability $1$ (repeated eigenvalues have measure zero), so the theorem as proved applies; the monic polynomial differs from $p_A$ by the sign $(-1)^4 = 1$ here, and in any case $p(A) = 0$ is unaffected by an overall sign.

(c) The proof needs $X^{-1}$, i.e. a basis of eigenvectors. For a defective $A$ there is none, so the identity $p_A(A) = Xp_A(\Lambda)X^{-1}$ has no meaning. The theorem survives: for the shear $S = \begin{pmatrix}1 & 1\\ 0 & 1\end{pmatrix}$, $p_S(S) = (S - I)^2 = N^2 = 0$ directly. In general LA-11 Theorem 3.16 proves it by triangularizing $A$ (the Schur form of LA-10 §3.5) and chaining the invariant subspaces, or from the adjugate identity $(A - \lambda I)\operatorname{adj}(A - \lambda I) = p_A(\lambda)I$ of LA-6, Theorem 3.14, or by continuity: defective matrices are limits of diagonalizable ones and $p_A(A)$ depends continuously on $A$.

**C2.** (a) Let $A^\top\mathbf{v} = \lambda\mathbf{v}$ with $|\lambda| = 1$, and let $m$ maximize $|v_j|$. From the proof of Theorem 3.11(b):
$$|v_m| = |\lambda||v_m| = \Big|\sum_ja_{jm}v_j\Big| \le \sum_ja_{jm}|v_j| \le \sum_ja_{jm}|v_m| = |v_m|.$$
The two ends are equal, so both inequalities are equalities.

*Second inequality is equality:* $\sum_ja_{jm}(|v_m| - |v_j|) = 0$ with every term $\ge0$ (as $a_{jm} > 0$ and $|v_j|\le|v_m|$). A sum of nonnegative terms is zero only if every term is zero, and $a_{jm} > 0$, so $|v_j| = |v_m|$ for **all** $j$.

*First inequality is equality:* the complex numbers $z_j = a_{jm}v_j$ all have the same argument. This is the equality case stated in the exercise; here is why it holds. F1 Theorem 3.22 is the two-term inequality, and its step 2 is the identity $|z + w|^2 = |z|^2 + 2\operatorname{Re}(z\bar w) + |w|^2$. For $n$ terms the same expansion gives $\big|\sum_jz_j\big|^2 = \sum_j|z_j|^2 + 2\sum_{j<l}\operatorname{Re}(z_j\bar z_l)$, while $\big(\sum_j|z_j|\big)^2 = \sum_j|z_j|^2 + 2\sum_{j<l}|z_j||z_l|$. Each $\operatorname{Re}(z_j\bar z_l) \le |z_j\bar z_l| = |z_j||z_l|$ (F1 step 3), so equality of the two squares (which is what $|\sum_jz_j| = \sum_j|z_j|$ says) forces $\operatorname{Re}(z_j\bar z_l) = |z_j||z_l|$ for every pair, i.e. $z_j\bar z_l$ is real and $\ge 0$. Writing $z_j = |z_j|e^{i\alpha_j}$, that says $\cos(\alpha_j - \alpha_l) = 1$ whenever both are nonzero: all nonzero $z_j$ have the same argument. Here every $z_j$ is nonzero, since $a_{jm} > 0$ and $|v_j| = |v_m| > 0$ by the previous paragraph. Since $a_{jm} > 0$, the $v_j$ all have the same argument $\theta$. Combined with equal moduli: $v_j = |v_m|e^{i\theta}$ for all $j$, i.e. $\mathbf{v} = c\mathbf{1}$ with $c\ne0$.

Then $A^\top\mathbf{v} = cA^\top\mathbf{1} = c\mathbf{1} = \mathbf{v}$, so $\lambda\mathbf{v} = \mathbf{v}$ and $\lambda = 1$. ∎

(b) The argument in (a) applies to every eigenvector of $A^\top$ with eigenvalue $1$ (which has $|\lambda| = 1$) and shows it is a multiple of $\mathbf{1}$. So $N(A^\top - I) = \operatorname{span}\{\mathbf{1}\}$, dimension $1$. Since $(A - I)^\top = A^\top - I$ and a matrix and its transpose have the same rank (LA-3), $\dim N(A - I) = n - \operatorname{rank}(A - I) = n - \operatorname{rank}(A^\top - I) = \dim N(A^\top - I) = 1$. The $\lambda = 1$ eigenspace of $A$ is a line; the steady state (the probability vector on it) is unique. ∎

(c) Order the eigenvalues with $\lambda_1 = 1$. By (b), $\operatorname{GM}(1) = 1$; since $A$ is diagonalizable, Theorem 3.20 gives $\operatorname{AM}(1) = 1$, so no other $\lambda_i$ equals $1$. By Theorem 3.11(b) every $|\lambda_i|\le1$, and by (a) any $|\lambda_i| = 1$ would force $\lambda_i = 1$; hence $|\lambda_i| < 1$ for $i\ge2$. From (3.8), $\mathbf{u}_k = c_1\mathbf{x}_1 + \sum_{i\ge2}c_i\lambda_i^k\mathbf{x}_i\to c_1\mathbf{x}_1$. By B4(b), each $\mathbf{x}_i$ with $i\ge2$ sums to zero, so $1 = \mathbf{1}^\top\mathbf{u}_k = c_1\mathbf{1}^\top\mathbf{x}_1$ for all $k$: the limit $c_1\mathbf{x}_1$ is the probability vector in $E_1$, which by (b) is the unique steady state $\mathbf{u}_\infty$, independent of $\mathbf{u}_0$. ∎ (B4(b) was used to identify the limit as $\mathbf{u}_\infty$ rather than some other multiple of $\mathbf{x}_1$. Without diagonalizability the conclusion still holds, but the transient terms carry polynomial factors, $k^j\lambda_i^k$, and the argument needs the Jordan form (LA-11 Theorem 3.10).)

**C3.** (a) Number the masses $1, \dots, n$ left to right, with displacements $x_i$ and the convention $x_0 = x_{n+1} = 0$ for the walls. The spring between masses $i - 1$ and $i$ is stretched by $x_i - x_{i-1}$ and pulls mass $i$ to the left with force $k(x_i - x_{i-1})$; the spring between $i$ and $i + 1$ is stretched by $x_{i+1} - x_i$ and pulls mass $i$ to the right with force $k(x_{i+1} - x_i)$. So
$$m\ddot x_i = -k(x_i - x_{i-1}) + k(x_{i+1} - x_i) = -k(2x_i - x_{i-1} - x_{i+1}),$$
and the wall convention makes the formula correct at $i = 1$ and $i = n$. In matrix form $mI\ddot{\mathbf{x}} = -kT_n\mathbf{x}$ with $(T_n)_{ii} = 2$, $(T_n)_{i,i\pm1} = -1$. (For $n = 2$ with $m_1 = m_2$ this is Example 4.3's $K$.)

(b) Row $i$ of $T_n\mathbf{v}_j$, for $1\le i\le n$, using the wall convention $\sin(0\cdot\theta_j) = 0$ for $i = 1$ and $\sin((n+1)\theta_j) = \sin(j\pi) = 0$ for $i = n$ (so the formula below is valid in the boundary rows too, where $T_n$ has only one $-1$):
$$(T_n\mathbf{v}_j)_i = -\sin((i-1)\theta_j) + 2\sin(i\theta_j) - \sin((i+1)\theta_j) = 2\sin(i\theta_j) - 2\sin(i\theta_j)\cos\theta_j = (2 - 2\cos\theta_j)\sin(i\theta_j),$$
by the identity $\sin((i-1)\theta) + \sin((i+1)\theta) = 2\sin(i\theta)\cos\theta$ (sum-to-product; expand both sides with the addition formula and the $\cos(i\theta)\sin\theta$ terms cancel). So $T_n\mathbf{v}_j = \mu_j\mathbf{v}_j$ with $\mu_j = 2 - 2\cos\theta_j = 4\sin^2(\theta_j/2)$. ∎

$\mathbf{v}_j\ne\mathbf{0}$: its first entry is $\sin\theta_j = \sin\frac{j\pi}{n+1}$, and $0 < \frac{j\pi}{n+1} < \pi$ for $1\le j\le n$, where sine is positive. Distinct eigenvalues: $\cos$ is strictly decreasing on $(0, \pi)$ and the $\theta_j$ are distinct points there, so the $\mu_j$ are distinct. $T_n$ is symmetric, so by Theorem 3.17(b) eigenvectors for distinct eigenvalues are orthogonal: $\mathbf{v}_j^\top\mathbf{v}_l = \sum_i\sin(i\theta_j)\sin(i\theta_l) = 0$ for $j\ne l$. (This is the discrete orthogonality of sines behind the discrete sine transform.)

(c) $K\mathbf{v} = \omega^2M\mathbf{v}$ with $K = kT_n$, $M = mI$ reads $kT_n\mathbf{v} = \omega^2m\mathbf{v}$, so $\omega_j^2 = \frac km\mu_j = \frac{4k}{m}\sin^2\frac{j\pi}{2(n+1)}$ and $\omega_j = 2\sqrt{k/m}\,\sin\frac{j\pi}{2(n+1)}$. For $n = 3$: $\theta_j = \pi/4, \pi/2, 3\pi/4$, so $\mu_j = 2 - 2\cos\theta_j = 2 - \sqrt2,\ 2,\ 2 + \sqrt2$, matching the project's `eigvalsh` output $0.5858, 2, 3.4142$. Mode shapes for $n = 3$: $\mathbf{v}_1 = (\sin\frac\pi4, \sin\frac\pi2, \sin\frac{3\pi}4) = (\tfrac{\sqrt2}2, 1, \tfrac{\sqrt2}2)$ (all in phase, middle largest), $\mathbf{v}_2 = (1, 0, -1)$ (middle mass still), $\mathbf{v}_3 = (\tfrac{\sqrt2}2, -1, \tfrac{\sqrt2}2)$ (alternating).

(d) With $L = (n+1)a$: $\frac{j\pi}{2(n+1)} = \frac{j\pi a}{2L}$. For fixed $j$ and $n\to\infty$ (so $a\to0$), $\sin\frac{j\pi a}{2L}\approx\frac{j\pi a}{2L}$, giving $\omega_j\approx2\sqrt{k/m}\cdot\frac{j\pi a}{2L} = \frac{j\pi}{L}\,a\sqrt{k/m}$. With $c = a\sqrt{k/m}$ this is $\omega_j = j\pi c/L$: the harmonics of a string of length $L$ fixed at both ends, $\omega_j = j\omega_1$, with mode shapes $\sin(j\pi x/L)$ (the continuum limit of $\sin(i\theta_j) = \sin(j\pi\,x_i/L)$ at $x_i = ia$). Expanding an initial shape in these modes is a Fourier sine series; F2 does the continuum version, and the wave equation $\ddot u = c^2u''$ is the $n\to\infty$ limit of (a). For finite $n$ the frequencies are *not* harmonic ($\sin$ is not linear): the chain is a dispersive medium. SM-5's Debye model of a solid keeps only the small-$k$ linear limit $\omega \approx ck$ of this dispersion relation, which is where its cutoff frequency comes from.
