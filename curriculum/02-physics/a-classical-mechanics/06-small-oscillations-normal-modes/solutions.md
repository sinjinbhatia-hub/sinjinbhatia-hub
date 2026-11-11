# CM-6 Solutions

Conventions as in the notes. $\lambda = \omega^2$ where convenient. Every number below is asserted in Part 9 of `project.py` (or in Parts 2–7 where the exercise reuses a worked example).

## A. Warm-up

**A1.** The springs give $V = \frac12k_1x_1^2 + \frac12k_2(x_2 - x_1)^2 + \frac12k_3x_2^2$, so
$$M = \begin{pmatrix}1 & 0\\0 & 2\end{pmatrix}, \qquad K = \begin{pmatrix}k_1 + k_2 & -k_2\\-k_2 & k_2 + k_3\end{pmatrix} = \begin{pmatrix}3 & -2\\-2 & 4\end{pmatrix}.$$
$\det(K - \lambda M) = (3 - \lambda)(4 - 2\lambda) - 4 = 2\lambda^2 - 10\lambda + 8 = 2(\lambda - 1)(\lambda - 4)$, so $\omega_1^2 = 1$ and $\omega_2^2 = 4$.

Shapes come from the first row of $(K - \lambda M)\mathbf{a} = 0$, which reads $(3 - \lambda)a_1 - 2a_2 = 0$. For $\lambda = 1$: $2a_1 = 2a_2$, so $\mathbf{a}_1 \propto (1, 1)^\top$. For $\lambda = 4$: $-a_1 = 2a_2$, so $\mathbf{a}_2 \propto (2, -1)^\top$.

$M$-norms: $(1,1)M(1,1)^\top = 1 + 2 = 3$ and $(2,-1)M(2,-1)^\top = 4 + 2 = 6$. So $\mathbf{a}_1 = (1, 1)^\top/\sqrt3$ and $\mathbf{a}_2 = (2, -1)^\top/\sqrt6$.

$M$-orthogonality: $(1,1)M(2,-1)^\top = 1\cdot2 + 2\cdot1\cdot(-1) = 0$. ✓ The ordinary dot product of the normalized shapes is $(2 - 1)/\sqrt{18} = 1/(3\sqrt2) = 0.2357 \ne 0$.

**A2.** With $m = \ell = 1$: $M = \begin{pmatrix}2 & 1\\1 & 1\end{pmatrix}$ and $K = g\operatorname{diag}(2, 1)$. Example 4.4 gives the squared norms $(1, \sqrt2)M(1, \sqrt2)^\top = 4 + 2\sqrt2 = 6.8284$ and $(1, -\sqrt2)M(1, -\sqrt2)^\top = 4 - 2\sqrt2 = 1.1716$. Dividing each shape by the square root of its norm (2.61313 and 1.08239):
$$A = \begin{pmatrix}0.38268 & 0.92388\\0.54120 & -1.30656\end{pmatrix}.$$
Check $A^\top MA$. First compute $MA = \begin{pmatrix}2(0.38268) + 0.54120 & 2(0.92388) - 1.30656\\0.38268 + 0.54120 & 0.92388 - 1.30656\end{pmatrix} = \begin{pmatrix}1.30656 & 0.54120\\0.92388 & -0.38268\end{pmatrix}$. Then $(A^\top MA)_{11} = 0.38268(1.30656) + 0.54120(0.92388) = 0.50000 + 0.50000 = 1$, $(A^\top MA)_{12} = 0.38268(0.54120) + 0.54120(-0.38268) = 0$, and $(A^\top MA)_{22} = 0.92388(0.54120) + (-1.30656)(-0.38268) = 0.50000 + 0.50000 = 1$. ✓

$A^\top KA$: by (3.11) its diagonal is $\omega_j^2 = (2 \mp \sqrt2)g = 5.7466$ and $33.4934$. Directly, $(A^\top KA)_{11} = g[2(0.38268)^2 + (0.54120)^2] = g(0.29289 + 0.29289) = 0.58579g = 5.7466$. The script checks the full matrix.

$AA^\top = M^{-1}$: from $A^\top MA = I$, the square matrix $A$ has left inverse $A^\top M$, so $A^{-1} = A^\top M$ and also $AA^\top M = I$. Multiply on the right by $M^{-1}$: $AA^\top = M^{-1}$. Numerically $AA^\top = \begin{pmatrix}0.14645 + 0.85355 & 0.20711 - 1.20711\\ \cdot & 0.29289 + 1.70711\end{pmatrix} = \begin{pmatrix}1 & -1\\-1 & 2\end{pmatrix}$, and $M^{-1} = \frac{1}{2 - 1}\begin{pmatrix}1 & -1\\-1 & 2\end{pmatrix}$. ✓ This is the completeness relation $\sum_j\mathbf{a}_j\mathbf{a}_j^\top = M^{-1}$, the $M$-version of $\sum_j\mathbf{e}_j\mathbf{e}_j^\top = I$.

**A3.** (a) $\mathbf{v} = (1, -1)^\top$: $\mathbf{v}^\top M\mathbf{v} = 1 - 2 - 2 + 1 = -2 < 0$. Positive definiteness fails (the eigenvalues of $M$ are $3$ and $-1$).

(b) $\mathbf{v}^\top M\mathbf{v} = v_1^2 + v_1v_2 + v_2^2 = (v_1 + \frac12v_2)^2 + \frac34v_2^2$, a sum of squares that vanishes only when $v_2 = 0$ and then $v_1 = 0$. So it is positive for $\mathbf{v} \ne 0$. Symmetry fails: $\langle\mathbf{e}_1, \mathbf{e}_2\rangle_M = M_{12} = 1$ but $\langle\mathbf{e}_2, \mathbf{e}_1\rangle_M = M_{21} = 0$.

(c) Derivation 1 gives $M_{ij} = \sum_\alpha m_\alpha\,\partial_i\mathbf{r}_\alpha\cdot\partial_j\mathbf{r}_\alpha$. That is symmetric, and $\dot q^\top M\dot q = \sum_\alpha m_\alpha|\dot{\mathbf{r}}_\alpha|^2 \ge 0$, so neither matrix can be a mass matrix. Matrix (a) would give some velocity negative kinetic energy. Matrix (b) is subtler: $T = \frac12\dot q^\top M\dot q$ sees only the symmetric part $\frac12(M + M^\top) = \begin{pmatrix}1 & \frac12\\\frac12 & 1\end{pmatrix}$, whose eigenvalues are $\frac12$ and $\frac32$. That part is a legitimate mass matrix. So if you write $T$ with a lopsided matrix, symmetrize it before you form $\langle\cdot,\cdot\rangle_M$ or call `eigh`.

**A4.** Since $M = I$, $\omega^2$ are the eigenvalues of $K$.

Equilibria: $\partial_xV = x(x^2 - 1) - 2(y - x)$ and $\partial_yV = 2(y - x)$. The second vanishes only on $y = x$; then the first is $x(x^2 - 1) = 0$. So the equilibria are $(0, 0)$, $(1, 1)$, $(-1, -1)$.

Hessian: $V_{xx} = 3x^2 - 1 + 2 = 3x^2 + 1$, $V_{xy} = -2$, $V_{yy} = 2$.

- At $(0, 0)$: $K = \begin{pmatrix}1 & -2\\-2 & 2\end{pmatrix}$, trace 3, determinant $-2$. The eigenvalues are $\frac{3 \pm\sqrt{17}}{2} = -0.5616$ and $3.5616$. $K$ is indefinite, so by Theorem 3.5 and the remark after it the equilibrium is unstable. Along the eigenvector of $-0.5616$, displacements grow like $e^{\gamma t}$ with $\gamma = \sqrt{0.5616} = 0.7494$. Along the other eigenvector the motion oscillates at $\sqrt{3.5616} = 1.8872$. It is a saddle.
- At $(\pm1, \pm1)$: $K = \begin{pmatrix}4 & -2\\-2 & 2\end{pmatrix}$, trace 6, determinant 4. The eigenvalues are $3 \mp\sqrt5 = 0.7639$ and $5.2361$, both positive. Both equilibria are stable minima, with $\omega = 0.8740$ and $2.2882$.

**A5.** Take $K/m = \omega_0^2T$ with $T = \operatorname{tridiag}(-1, 2, -1)$ of size $3\times3$.

- $T(1/\sqrt2, 1, 1/\sqrt2)^\top = (\sqrt2 - 1,\ 2 - \sqrt2,\ \sqrt2 - 1)^\top = (2 - \sqrt2)(1/\sqrt2, 1, 1/\sqrt2)^\top$, because $(2 - \sqrt2)/\sqrt2 = \sqrt2 - 1$. So $\omega^2/\omega_0^2 = 2 - \sqrt2$, with no interior node.
- $T(1, 0, -1)^\top = (2,\ -1 + 1,\ -2)^\top = 2(1, 0, -1)^\top$. So $\omega^2/\omega_0^2 = 2$, with one node, at the middle bead.
- $T(1/\sqrt2, -1, 1/\sqrt2)^\top = (\sqrt2 + 1,\ -2 - \sqrt2,\ \sqrt2 + 1)^\top = (2 + \sqrt2)(1/\sqrt2, -1, 1/\sqrt2)^\top$. So $\omega^2/\omega_0^2 = 2 + \sqrt2$, with two nodes (two sign changes).

(4.5) with $N = 3$ gives shapes $\sin(jk\pi/4)$: $(\frac{1}{\sqrt2}, 1, \frac{1}{\sqrt2})$, $(1, 0, -1)$, $(\frac{1}{\sqrt2}, -1, \frac{1}{\sqrt2})$. It gives $\omega_j^2/\omega_0^2 = 4\sin^2(j\pi/8) = 2 - 2\cos(j\pi/4) = 2 - \sqrt2,\ 2,\ 2 + \sqrt2$. Mode $j$ has $j - 1$ nodes. ✓

**A6.** $\omega_1 = \sqrt{g/\ell} = \sqrt{16.35} = 4.0435\ \mathrm{s^{-1}}$. $T_{\text{beat}} = 2\pi/(\omega_2 - \omega_1) = 12$ s gives $\Delta\omega = \pi/6 = 0.5236\ \mathrm{s^{-1}}$.

Exact: $\omega_2^2 - \omega_1^2 = 2k/m$, so $k = \frac m2(\omega_2 - \omega_1)(\omega_2 + \omega_1) = m\,\Delta\omega\,(\omega_1 + \Delta\omega/2) = 0.2(0.5236)(4.3053) = 0.4508$ N/m.

Weak coupling: $k \approx m\omega_1\Delta\omega = 0.2(4.0435)(0.5236) = 0.4234$ N/m.

The ratio is exactly $1 + \Delta\omega/(2\omega_1) = 1.065$, so the approximation is 6.1% low. The coupling is "weak" ($2k/m \approx 4.5 \ll g/\ell = 16.35$), but not weak enough for 1% work. A rubber band of this stiffness is easy to find.

## B. Core

**B1.** (a) Substitute $q = A\xi$ into $M\ddot q + C\dot q + Kq = F$ and multiply on the left by $A^\top$: $A^\top MA\,\ddot\xi + A^\top CA\,\dot\xi + A^\top KA\,\xi = A^\top F$. By (3.11), $A^\top MA = I$ and $A^\top KA = \Omega^2$. With $C = \alpha M + \beta K$, linearity gives $A^\top CA = \alpha A^\top MA + \beta A^\top KA = \alpha I + \beta\Omega^2$, a diagonal matrix. Row $j$ is $\ddot\xi_j + (\alpha + \beta\omega_j^2)\dot\xi_j + \omega_j^2\xi_j = f_j$. Define $\zeta_j$ by $2\zeta_j\omega_j = \alpha + \beta\omega_j^2$ to get (3.15). ∎

(b) From Example 4.2, $\mathbf{a}_1 = (0.45970, 0.62796)^\top$ and $\mathbf{a}_2 = (0.88807, -0.32506)^\top$. With $C = \operatorname{diag}(1, 0)$, $(A^\top CA)_{ij} = a_{1i}a_{1j}$, the product of the first components:
$$A^\top CA = \begin{pmatrix}0.45970^2 & 0.45970\cdot0.88807\\ \cdot & 0.88807^2\end{pmatrix} = \begin{pmatrix}0.21132 & 0.40825\\0.40825 & 0.78868\end{pmatrix} = \begin{pmatrix}\frac{3 - \sqrt3}{6} & \frac{1}{\sqrt6}\\\frac{1}{\sqrt6} & \frac{3 + \sqrt3}{6}\end{pmatrix}.$$
The exact forms use $a_{11}^2 = 1/(3 + \sqrt3)$ (Example 4.2's norm $4.73205 = 3 + \sqrt3$) and $a_{12}^2 = 1/(3 - \sqrt3)$. Not diagonal. And $\alpha M + \beta K = \begin{pmatrix}\alpha + 2\beta & -\beta\\-\beta & 2\alpha + 2\beta\end{pmatrix}$. Matching the off-diagonal $0$ forces $\beta = 0$. Matching the $(2,2)$ entry $0$ then forces $\alpha = 0$, but the $(1,1)$ entry must be 1. Contradiction.

(c) $\omega_1 = 0.79623$, $\omega_2 = 1.53819$. $\zeta_1 = \dfrac{0.1 + 0.05(0.63397)}{2(0.79623)} = \dfrac{0.13170}{1.59246} = 0.08270$. $\zeta_2 = \dfrac{0.1 + 0.05(2.36603)}{2(1.53819)} = \dfrac{0.21830}{3.07638} = 0.07096$.

(d) We need $\alpha + \beta\omega_j^2 = 2\zeta\omega_j$ for $j = 1, 2$. Subtracting the two equations gives $\beta(\omega_2^2 - \omega_1^2) = 2\zeta(\omega_2 - \omega_1)$, so
$$\beta = \frac{2\zeta}{\omega_1 + \omega_2}, \qquad \alpha = 2\zeta\omega_1 - \beta\omega_1^2 = \frac{2\zeta\,\omega_1\omega_2}{\omega_1 + \omega_2}.$$
Check: $\alpha + \beta\omega_j^2 = \frac{2\zeta}{\omega_1 + \omega_2}(\omega_1\omega_2 + \omega_j^2) = \frac{2\zeta\omega_j(\omega_{\text{other}} + \omega_j)}{\omega_1 + \omega_2} = 2\zeta\omega_j$. ✓ With $\zeta = 0.02$, $\omega_1 + \omega_2 = 2.33442$ and $\omega_1\omega_2 = \sqrt{\lambda_1\lambda_2} = \sqrt{3/2} = 1.22474$: $\beta = 0.017135$, $\alpha = 0.020986$.

(e) The off-diagonal entry $0.40825$ means that the damping force produced by mode 1's velocity pushes on mode 2, and vice versa. Energy leaks between modes through the dashpot, so a motion started in one mode does not stay in it. The free motions are then complex modes in which the two masses are neither in phase nor in antiphase, found from the $2n$-dimensional first-order system.

**B2.** (a) $\mathbf{r}_1 = \ell(\sin\theta_1, -\cos\theta_1)$ and $\mathbf{r}_2 = \mathbf{r}_1 + \ell(\sin\theta_2, -\cos\theta_2)$. As in Example 4.4, $T = \frac12m_1\ell^2\dot\theta_1^2 + \frac12m_2\ell^2[\dot\theta_1^2 + \dot\theta_2^2 + 2\dot\theta_1\dot\theta_2\cos(\theta_1 - \theta_2)]$ and $V = -(m_1 + m_2)g\ell\cos\theta_1 - m_2g\ell\cos\theta_2$. At second order $\cos(\theta_1 - \theta_2) \to 1$ in $T$ and $-\cos\theta \to -1 + \theta^2/2$ in $V$. So $M_{ij} = \partial^2T/\partial\dot\theta_i\partial\dot\theta_j = \ell^2\begin{pmatrix}m_1 + m_2 & m_2\\m_2 & m_2\end{pmatrix}$ and $K = g\ell\operatorname{diag}(m_1 + m_2, m_2)$.

(b) Put $\lambda = \omega^2\ell/g$ and factor out $g\ell$: $K - \omega^2M = g\ell\begin{pmatrix}(m_1 + m_2)(1 - \lambda) & -m_2\lambda\\-m_2\lambda & m_2(1 - \lambda)\end{pmatrix}$. The determinant vanishes when $(m_1 + m_2)m_2(1 - \lambda)^2 = m_2^2\lambda^2$. Divide by $(m_1 + m_2)m_2$: $(1 - \lambda)^2 = \mu\lambda^2$, so $1 - \lambda = \pm\sqrt\mu\,\lambda$ and $\lambda = 1/(1 \pm\sqrt\mu)$. The first row gives $(m_1 + m_2)(1 - \lambda)a_1 = m_2\lambda a_2$, so $a_2/a_1 = (1 - \lambda)/(\mu\lambda) = \pm\sqrt\mu\lambda/(\mu\lambda) = \pm1/\sqrt\mu$. The $+$ sign gives the smaller $\lambda$: the slow mode is in phase.

(c) Equal masses: $\mu = \frac12$, so $\lambda = \frac{1}{1 \pm 1/\sqrt2} = \frac{\sqrt2}{\sqrt2 \pm 1} = \sqrt2(\sqrt2 \mp 1) = 2 \mp\sqrt2$, and $a_2/a_1 = \pm\sqrt2$. This is (4.4). ✓ For $m_1 = 3$, $m_2 = 1$: $\mu = \frac14$, $\lambda = 1/(1 \pm \frac12) = \frac23, 2$, so $\omega^2 = \frac{2g}{3\ell}$ and $\frac{2g}{\ell}$, with shapes $(1, 2)^\top$ and $(1, -2)^\top$.

(d) $m_2 \to 0$ ($\mu \to 0$): both $\lambda \to 1$. The heavy upper pendulum swings at $\sqrt{g/\ell}$ as if alone, and the light lower pendulum is a pendulum of the same natural frequency driven at its pivot, so it is driven at resonance. That is why $a_2/a_1 = \pm1/\sqrt\mu \to \infty$: the small bob swings enormously. The splitting $\lambda \approx 1 \mp\sqrt\mu$ grows like the square root of the mass ratio. This is the principle of the tuned mass damper. $m_1 \to 0$ ($\mu \to 1$): $\lambda \to \frac12$ with $a_2/a_1 = 1$. The two strings form one straight pendulum of length $2\ell$, with $\omega^2 = g/(2\ell)$. ✓ The other $\lambda \to \infty$: a massless joint cannot sustain a kink, so that mode's frequency diverges.

**B3.** (a) $M^{-1} = L^{-\top}L^{-1}$, so $M^{-1}K = L^{-\top}L^{-1}KL^{-\top}L^\top = L^{-\top}\tilde KL^\top$. Write $S = L^\top$. Then $M^{-1}K = S^{-1}\tilde KS$: similar to the symmetric $\tilde K$. By the spectral theorem, $\tilde K = Q\Lambda Q^\top$ with $\Lambda$ real diagonal and $Q$ orthogonal. So $M^{-1}K = (S^{-1}Q)\Lambda(S^{-1}Q)^{-1}$ is diagonalizable with real eigenvalues $\Lambda$ and eigenvector matrix $S^{-1}Q = L^{-\top}Q$. Those eigenvectors are exactly the $\mathbf{a}_j$ of Theorem 3.4. Symmetry of $M^{-1}K$ was never needed, only similarity to a symmetric matrix.

(b) $M = \operatorname{diag}(1, 2)$, so $L = \operatorname{diag}(1, \sqrt2)$ and $\tilde K = L^{-1}KL^{-1} = \begin{pmatrix}2 & -1/\sqrt2\\-1/\sqrt2 & 1\end{pmatrix}$, symmetric. Its trace is 3 and its determinant $2 - \frac12 = \frac32$, so $\lambda^2 - 3\lambda + \frac32 = 0$ and $\lambda = (3 \mp\sqrt3)/2 = 0.63397, 2.36603$, the same as Example 4.2. Eigenvectors from the first row, $(2 - \lambda)y_1 = y_2/\sqrt2$, give $y_2/y_1 = \sqrt2(2 - \lambda)$. For $\lambda_1$ the ratio is $1.93185$, so $\mathbf{y}_1 = (1, 1.93185)/\sqrt{4.73205} = (0.45970, 0.88807)$. For $\lambda_2$ it is $-0.51764$, so $\mathbf{y}_2 = (0.88807, -0.45970)$. These are orthonormal in the ordinary sense: $0.45970\cdot0.88807 - 0.88807\cdot0.45970 = 0$. Map back with $\mathbf{a} = L^{-\top}\mathbf{y} = (y_1, y_2/\sqrt2)$: $\mathbf{a}_1 = (0.45970, 0.62796)$ and $\mathbf{a}_2 = (0.88807, -0.32506)$, matching Example 4.2. ✓

(c) The matrix of $M^{-1}K$ in the basis $\{\mathbf{a}_j\}$ is $A^{-1}(M^{-1}K)A = A^\top M\,M^{-1}KA = A^\top KA = \Omega^2$, using $A^{-1} = A^\top M$ (A2). `np.linalg.eig` on $M^{-1}K$ returns the right eigenvalues, because they are similarity invariants. But it normalizes each eigenvector to unit Euclidean length, so $\mathbf{v}_j^\top M\mathbf{v}_j \ne 1$, and (3.13) with those vectors gives wrong amplitudes. Inside a degenerate eigenspace it can also return vectors that are not $M$-orthogonal at all. If you must use it, rescale $\mathbf{v}_j \leftarrow \mathbf{v}_j/\sqrt{\mathbf{v}_j^\top M\mathbf{v}_j}$ and $M$-orthogonalize within each degenerate block. Better to call `eigh(K, M)`.

**B4.** (a) With $q = A\xi$: $\frac12\dot q^\top M\dot q = \frac12\dot\xi^\top(A^\top MA)\dot\xi = \frac12\sum_j\dot\xi_j^2$ and $\frac12q^\top Kq = \frac12\xi^\top(A^\top KA)\xi = \frac12\sum_j\omega_j^2\xi_j^2$. Add them: $E = \sum_jE_j$. Each $E_j$ is constant because $\dot E_j = \dot\xi_j(\ddot\xi_j + \omega_j^2\xi_j) = 0$ by (3.12).

(b) $m = \ell = 1$, so $M = I$, and $\xi_j(0) = \mathbf{a}_j^\top q_0 = \theta_0/\sqrt2$ for both modes, with $\dot\xi_j(0) = 0$. So $E_1 = \frac12(9.81)\frac{\theta_0^2}{2} = 2.4525\,\theta_0^2$ and $E_2 = \frac12(11.81)\frac{\theta_0^2}{2} = 2.9525\,\theta_0^2$ (joules, with $\theta_0$ in radians). Their sum is $5.405\,\theta_0^2 = \frac12K_{11}\theta_0^2 = \frac12(10.81)\theta_0^2$. ✓

(c) From (4.1), $\theta_1 = \theta_0\cos(\varepsilon t)\cos(\bar\omega t)$, which is $0$ at $\varepsilon t = \pi/2$. Differentiating, $\dot\theta_1 = -\theta_0[\varepsilon\sin(\varepsilon t)\cos(\bar\omega t) + \bar\omega\cos(\varepsilon t)\sin(\bar\omega t)]$. At $\varepsilon t = \pi/2$ this is $-\varepsilon\theta_0\cos(\bar\omega T_{\text{transfer}})$, so $|\dot\theta_1| \le \varepsilon\theta_0 = 0.152\theta_0$; the script finds $0.119\theta_0$. Pendulum 1's own energy $\frac12m\ell^2\dot\theta_1^2 + \frac12mg\ell\theta_1^2$ is then $0.13\%$ of $E$.

Reconciliation: $E_1$ and $E_2$ are energies of *modes*, not of pendulums. At $T_{\text{transfer}}$ both modes are still present with their original energies. Their contributions cancel at pendulum 1 (in-phase plus out-of-phase shapes with opposite signs) and add at pendulum 2. Split by pendulum, the energy contains cross terms $\xi_1\xi_2$ and $\dot\xi_1\dot\xi_2$ whose slowly varying part oscillates at $\omega_2 - \omega_1$, and the spring's energy belongs to neither pendulum. Energy flows between pendulums but never between modes.

**B5.** (a) *From Theorem 3.2.* For the string, $M = mI$ and the shapes $\mathbf{s}^{(j)}_k = \sin\frac{jk\pi}{N+1}$ are eigenvectors of $K$ with distinct $\omega_j^2$ ((4.5); $\sin$ is increasing on $(0, \pi/2)$). Theorem 3.2 gives $\mathbf{s}^{(j)\top}M\mathbf{s}^{(j')} = m\,\mathbf{s}^{(j)}\cdot\mathbf{s}^{(j')} = 0$ for $j \ne j'$. No trigonometry was used.

*Direct summation, both cases at once.* Let $\theta_p = p\pi/(N+1)$ and $S(p) = \sum_{k=1}^N\cos(k\theta_p)$. For $0 < p < 2(N+1)$, $\sin(\theta_p/2) \ne 0$. The telescoping identity $2\sin\frac\theta2\cos k\theta = \sin(k + \frac12)\theta - \sin(k - \frac12)\theta$ summed over $k = 1, \dots, N$ gives $S(p) = \dfrac{\sin(N + \frac12)\theta_p}{2\sin(\theta_p/2)} - \dfrac12$. Since $(N + 1)\theta_p = p\pi$, $\sin(N + \frac12)\theta_p = \sin(p\pi - \frac{\theta_p}{2}) = -\cos(p\pi)\sin\frac{\theta_p}{2} = -(-1)^p\sin\frac{\theta_p}{2}$. So
$$S(p) = -\tfrac12(-1)^p - \tfrac12 = \begin{cases}-1 & p \text{ even}\\0 & p \text{ odd}\end{cases}, \qquad S(0) = N.$$
Now $\sin a\sin b = \frac12[\cos(a - b) - \cos(a + b)]$, so $\sum_k\sin\frac{jk\pi}{N+1}\sin\frac{j'k\pi}{N+1} = \frac12[S(j - j') - S(j + j')]$. Since $S$ is even in $p$, use $|j - j'|$; both $|j - j'|$ and $j + j'$ lie in $[0, 2N]$. If $j \ne j'$, then $j - j'$ and $j + j'$ have the same parity, so $S(j - j') = S(j + j')$ and the sum is $0$. If $j = j'$, the sum is $\frac12[N - S(2j)] = \frac12[N + 1] = \frac{N+1}{2}$, since $2j$ is even and nonzero. ∎

(b) $\mathbf{s}^{(j)\top}M\mathbf{s}^{(j)} = m(N+1)/2$, so $\mathbf{a}_j = \sqrt{\dfrac{2}{m(N+1)}}\,\Big(\sin\dfrac{jk\pi}{N+1}\Big)_{k=1}^N$. This is what `string_closed_form` uses.

(c) Write $q_0 = \sum_jc_j\mathbf{s}^{(j)}$ and dot with $\mathbf{s}^{(j)}$: by (a), $c_j = \frac{2}{N+1}\,\mathbf{s}^{(j)}\cdot q_0 = \frac13\,\mathbf{s}^{(j)}\cdot q_0$. The values $\sin(k\pi/6)$ for $k = 1..5$ are $\frac12, \frac{\sqrt3}{2}, 1, \frac{\sqrt3}{2}, \frac12$.
- $j = 1$: $\frac13(\frac12 + \sqrt3 + 3 + \sqrt3 + \frac12) = \frac{4 + 2\sqrt3}{3} = 2.48803$.
- $j = 3$: $\sin(k\pi/2) = 1, 0, -1, 0, 1$, so $c_3 = \frac13(1 - 3 + 1) = -\frac13$.
- $j = 5$: $\sin(5k\pi/6) = \frac12, -\frac{\sqrt3}{2}, 1, -\frac{\sqrt3}{2}, \frac12$, so $c_5 = \frac13(\frac12 - \sqrt3 + 3 - \sqrt3 + \frac12) = \frac{4 - 2\sqrt3}{3} = 0.17863$.
- $j = 2, 4$: $\sin\frac{j(6 - k)\pi}{6} = \sin(j\pi - \frac{jk\pi}{6}) = -(-1)^j\sin\frac{jk\pi}{6}$. For even $j$ the shape is odd about the middle bead, while $q_0$ is even about it, so their dot product vanishes: $c_2 = c_4 = 0$.

Check: $\sum_jc_j\mathbf{s}^{(j)} = (1, 2, 3, 2, 1)$ (asserted).

(d) Released from rest, $\dot\xi_j(0) = 0$, so by (3.13) $y_k(t) = \sum_{j \in\{1,3,5\}}c_j\sin\frac{jk\pi}{6}\cos\omega_jt$ with $\omega_j = 2\omega_0\sin\frac{j\pi}{12}$: $\omega_1 = 0.5176\omega_0$, $\omega_3 = 1.4142\omega_0$, $\omega_5 = 1.9319\omega_0$. The ratios are irrational: $\omega_3/\omega_1 = \sin\frac{\pi}{4}/\sin\frac{\pi}{12} = 1 + \sqrt3$ and $\omega_5/\omega_1 = \cot\frac{\pi}{12} = 2 + \sqrt3$. (Non-integer alone would not be enough; a rational ratio still gives a periodic motion.) So there is no common period, and the shape never exactly repeats.

**B6.** The normalized modes are $\mathbf{a}_0 = (1, 1, 1)^\top/\sqrt{44}$ ($\omega = 0$), $\mathbf{a}_s = (1, 0, -1)^\top/\sqrt{32}$ ($\omega_s = 0.25$), and $\mathbf{a}_a = (1, -\frac83, 1)^\top/10.832$ ($\omega_a = 0.47871$), where $\sqrt{44}$, $\sqrt{32}$ and $10.832$ are $\sqrt{\mathbf{a}^\top M\mathbf{a}}$ for the unnormalized shapes.

(a) $q_0 = 0$, so every $\xi_j(0) = 0$. $Mv_0 = (16, 0, 0)^\top$, so $\dot\xi_j(0) = \mathbf{a}_j^\top Mv_0 = 16\,a_{1j}$: $\dot\xi_0 = 16/\sqrt{44} = 2.41209$, $\dot\xi_s = 16/\sqrt{32} = 2.82843$, $\dot\xi_a = 16/10.832 = 1.47710$.

(b) By (3.13), $q(t) = \mathbf{a}_0\dot\xi_0\,t + \mathbf{a}_s\dot\xi_s\frac{\sin\omega_st}{\omega_s} + \mathbf{a}_a\dot\xi_a\frac{\sin\omega_at}{\omega_a}$. The coefficient vectors are $\mathbf{a}_j\dot\xi_j = 16\,a_{1j}\mathbf{a}_j$:
- zero mode: $\frac{16}{44}(1, 1, 1) = \frac{4}{11}(1, 1, 1)$;
- symmetric stretch: $\frac{16}{32}(1, 0, -1) = \frac12(1, 0, -1)$;
- antisymmetric stretch: $\frac{16}{117.33}(1, -\frac83, 1) = \frac{3}{22}(1, -\frac83, 1) = (\frac{3}{22}, -\frac{4}{11}, \frac{3}{22})$, using $117.33 = 352/3$.

So
$$x_1 = \tfrac{4}{11}t + 2\sin(0.25t) + 0.28485\sin(0.47871t),\quad x_2 = \tfrac{4}{11}t - 0.75961\sin(0.47871t),\quad x_3 = \tfrac{4}{11}t - 2\sin(0.25t) + 0.28485\sin(0.47871t),$$
where $\frac{1/2}{0.25} = 2$, $\frac{3/22}{0.47871} = 0.28485$, $\frac{4/11}{0.47871} = 0.75961$. Check at $t = 0$: $\dot x_1 = \frac4{11} + \frac12 + \frac3{22} = 1$, $\dot x_2 = \frac4{11} - \frac4{11} = 0$, $\dot x_3 = \frac4{11} - \frac12 + \frac3{22} = 0$. ✓

(c) $x_{\text{cm}} = (16x_1 + 12x_2 + 16x_3)/44$. The symmetric-stretch terms give $16(2) - 16(2) = 0$. The antisymmetric terms give $\frac{3}{22}(16 + 16) - \frac{4}{11}(12) = \frac{48}{11} - \frac{48}{11} = 0$. So $x_{\text{cm}} = \frac{4}{11}t$, and the total momentum $P = 16$ over total mass $44$ is $\frac{4}{11}$. ✓ The only term that moves the center of mass is the zero mode.

(d) Let $\mathbf{a}_0 = c(1, \dots, 1)^\top$ with $c = (\sum_im_i)^{-1/2}$ so that $\mathbf{a}_0^\top M\mathbf{a}_0 = 1$. Then $\dot\xi_0 = \mathbf{a}_0^\top M\dot q = c\sum_im_i\dot q_i = cP$. Conservation: $\ddot\xi_0 = -0\cdot\xi_0 = 0$. Or directly, $\dot P = \mathbf{1}^\top M\ddot q = -\mathbf{1}^\top Kq = -(K\mathbf{1})^\top q = 0$, using $K = K^\top$ and $K\mathbf{1} = 0$. For a vibrational mode $j$, Theorem 3.2 gives $\mathbf{a}_0^\top M\mathbf{a}_j = c\sum_im_ia_{ij} = 0$, so the momentum it carries, $\sum_im_ia_{ij}\dot\xi_j$, is zero. This is CM-4's Noether pair (translation invariance, momentum) seen inside the modal expansion.

(e) `eigh` returns the zero eigenvalue as a tiny number of either sign; here it printed $-0$. `np.sqrt` of a tiny negative number is `nan`, and $0/0$ is `nan`. Either way the `nan` spreads through $q(t)$. `fit_initial_conditions` calls `frequencies`, which flags $|\omega^2| < 10^{-9}\max|\omega^2|$ as zero, sets $\omega = 0$ there, and uses $t$ in place of $\sin(\omega t)/\omega$ for those modes.

**B7.** (a) With $m = \ell = 1$ and $\Delta = \theta_1 - \theta_2$, $L = \frac12[2\dot\theta_1^2 + \dot\theta_2^2 + 2\dot\theta_1\dot\theta_2\cos\Delta] + g(2\cos\theta_1 + \cos\theta_2)$.
1. $\partial L/\partial\dot\theta_1 = 2\dot\theta_1 + \dot\theta_2\cos\Delta$. Its time derivative is $2\ddot\theta_1 + \ddot\theta_2\cos\Delta - \dot\theta_2(\dot\theta_1 - \dot\theta_2)\sin\Delta$. And $\partial L/\partial\theta_1 = -\dot\theta_1\dot\theta_2\sin\Delta - 2g\sin\theta_1$. Euler–Lagrange, after the $\dot\theta_1\dot\theta_2$ terms cancel: $2\ddot\theta_1 + \cos\Delta\,\ddot\theta_2 = -\dot\theta_2^2\sin\Delta - 2g\sin\theta_1$.
2. $\partial L/\partial\dot\theta_2 = \dot\theta_2 + \dot\theta_1\cos\Delta$. Its time derivative is $\ddot\theta_2 + \ddot\theta_1\cos\Delta - \dot\theta_1(\dot\theta_1 - \dot\theta_2)\sin\Delta$. And $\partial L/\partial\theta_2 = \dot\theta_1\dot\theta_2\sin\Delta - g\sin\theta_2$. So $\cos\Delta\,\ddot\theta_1 + \ddot\theta_2 = \dot\theta_1^2\sin\Delta - g\sin\theta_2$.

Linearizing ($\cos\Delta \to 1$, $\sin\theta \to \theta$, drop $\dot\theta^2$) recovers $M\ddot\theta = -K\theta$ of Example 4.4. ✓ The reference code (Part 9 of the script):

```python
def double_pendulum_rhs(t, y):
    t1, t2, w1, w2 = y; c, s = np.cos(t1 - t2), np.sin(t1 - t2)
    acc = np.linalg.solve([[2.0, c], [c, 1.0]], [-w2**2 * s - 2 * g * np.sin(t1), w1**2 * s - g * np.sin(t2)])
    return [w1, w2, acc[0], acc[1]]

modes44 = (np.array([2 - SQ2, 2 + SQ2]) * g, A44, M44)
for th0 in (0.1, 1.0):
    q_d, xi_d, _ = fit_initial_conditions([th0, 0.0], [0, 0.0], modes44)
    sol = solve_ivp(double_pendulum_rhs, (0, 10), [th0, 0, 0, 0], t_eval=t_d,
                    method="DOP853", rtol=1e-10, atol=1e-12)
    Xi, Xid = sol.y[:2].T @ M44 @ A44, sol.y[2:].T @ M44 @ A44
    E_t = 0.5 * (Xid**2 + modes44[0] * Xi**2)
```

(b) $\max_t|q_{\text{exact}} - q_{\text{modal}}|/\theta_0$ over 10 s is $0.171$ at $0.1$ rad and $1.434$ at $1.0$ rad.

(c) The relative swing of $E_1$ is $0.0085$ at $0.1$ rad and $0.852$ at $1.0$ rad.

(d) The unnormalized shapes are $\mathbf{b}_\pm = (1, \pm\sqrt2)^\top$, with $\mathbf{b}_\pm^\top M\mathbf{b}_\pm = 4 \pm 2\sqrt2 = 2(2 \pm\sqrt2)$. Here $Mq_0 = (2\theta_0, \theta_0)^\top$, so $\mathbf{b}_\pm^\top Mq_0 = (2 \pm\sqrt2)\theta_0$. Then $\xi_\pm(0) = \frac{(2 \pm\sqrt2)\theta_0}{\sqrt{2(2 \pm\sqrt2)}}$, $\xi_\pm^2 = \frac{(2 \pm\sqrt2)\theta_0^2}{2}$, and with $\omega_\pm^2 = (2 \mp\sqrt2)g$,
$$E_\pm = \tfrac12\omega_\pm^2\xi_\pm^2 = \tfrac12(2 \mp\sqrt2)(2 \pm\sqrt2)\frac{g\theta_0^2}{2} = \tfrac12g\theta_0^2,$$
because $(2 - \sqrt2)(2 + \sqrt2) = 2$. The total is $g\theta_0^2 = \frac12q_0^\top Kq_0 = \frac12(2g)\theta_0^2$. ✓ For $\theta_0 = 0.1$: $\xi = (0.13066, 0.054120)$ and $E_1 = E_2 = 0.04905$.

(e) Position error measures phase. In the exact system each mode's frequency shifts with amplitude, as for a single pendulum. A small frequency error $\delta\omega$ makes the phase drift by $\delta\omega\,t$, and after enough cycles (here about nine periods of the fast mode in 10 s) that is an $O(1)$ error in $q$, even though each mode keeps its energy. Energy exchange between modes is a different failure: the nonlinear terms couple the modes, and the decomposition itself stops describing the motion. At 0.1 rad the modes are still the right description, with slightly wrong frequencies. At 1.0 rad they are not: $E_1$ swings by 85%, and the motion is heading toward the chaotic regime of the double pendulum. Use the modal energies to decide whether the mode picture holds. Use the position error to decide whether linear frequencies are good enough for your time horizon.

**B8.** The modes are $\mathbf{a}_1 = (\frac12, \frac{1}{\sqrt2}, \frac12)$, $\mathbf{a}_2 = (\frac{1}{\sqrt2}, 0, -\frac{1}{\sqrt2})$, $\mathbf{a}_3 = (\frac12, -\frac{1}{\sqrt2}, \frac12)$, with $\omega^2 = 2 - \sqrt2,\ 2,\ 2 + \sqrt2$. This is the $N = 3$ string of A5 with $M = I$.

(a) At $\omega = 0$ the steady state solves $K\hat q = \mathbf{F}$, so $\hat q = K^{-1}\mathbf{F}$, and (3.16) gives $\hat q(0) = \sum_j\mathbf{a}_j\mathbf{a}_j^\top\mathbf{F}/\omega_j^2$. These agree for every $\mathbf{F}$, so $K^{-1} = \sum_j\mathbf{a}_j\mathbf{a}_j^\top/\omega_j^2$. Directly: $A^\top KA = \Omega^2$ gives $K = A^{-\top}\Omega^2A^{-1}$, hence $K^{-1} = A\Omega^{-2}A^\top$ (all $\omega_j > 0$ here), and $A\Omega^{-2}A^\top = \sum_j\mathbf{a}_j\mathbf{a}_j^\top/\omega_j^2$ column by column.

(b) Directly: $K^{-1} = \frac14\begin{pmatrix}3 & 2 & 1\\2 & 4 & 2\\1 & 2 & 3\end{pmatrix}$ (check $K\cdot$ first column: $2\cdot\frac34 - \frac12 = 1$, $-\frac34 + 1 - \frac14 = 0$, $-\frac12 + \frac12 = 0$). So $q = (\frac34, \frac12, \frac14)$. Modal: mode $j$ contributes $\mathbf{a}_j\,a_{1j}/\omega_j^2$:
- mode 1: $\frac{1/2}{2 - \sqrt2}\mathbf{a}_1 = 0.85355\,\mathbf{a}_1 = (0.42678, 0.60355, 0.42678)$;
- mode 2: $\frac{1/\sqrt2}{2}\mathbf{a}_2 = (0.25, 0, -0.25)$;
- mode 3: $\frac{1/2}{2 + \sqrt2}\mathbf{a}_3 = 0.14645\,\mathbf{a}_3 = (0.07322, -0.10355, 0.07322)$.

The sum is $(0.75, 0.5, 0.25)$. ✓ Mode 1 dominates, carrying 57% of $q_1$: the weight $1/\omega_j^2$ makes the softest mode the most compliant. This is why truncating to the lowest modes is accurate for slow loads.

(c) $a_{22} = 0$: mode 2 has a node at mass 2. A force on mass 2 has $\mathbf{a}_2^\top\mathbf{F} = a_{22} = 0$, so the $j = 2$ term of (3.16) is absent from every component of $\hat q$. Only two peaks appear (the script finds them at $0.765$ and $1.847$). For a force on mass 1, the response of mass 2 is $\hat q_2 = \sum_ja_{2j}a_{1j}/(\cdots)$. Its $j = 2$ term has the factor $a_{22} = 0$, so again two peaks. The two statements are the same fact, because the receptance matrix $\sum_j\mathbf{a}_j\mathbf{a}_j^\top/(\cdots)$ is symmetric (Maxwell reciprocity): a node cannot be excited, and it cannot see.

(d) $\zeta_j = (0.02 + 0.01\omega_j^2)/(2\omega_j) = 0.016892,\ 0.014142,\ 0.014651$. At $\omega = \omega_j$ the $j$-th denominator is $2i\zeta_j\omega_j^2$, so $|\hat q_1| \approx a_{1j}^2/(2\zeta_j\omega_j^2) = 12.632,\ 8.839,\ 2.499$. The full sum gives $12.649,\ 8.853,\ 2.557$. The single-mode estimate is good to 0.2% on the two big peaks and 2.3% on the smallest, where the off-resonance tails of the other two modes matter most. This is how modal testing reads $\zeta_j$ off a measured peak.

## C. Challenge

**C1.** (a) Label beads $(r, c)$, $r, c = 1..5$, with heights $y_{r,c}$. The frame is $y = 0$ wherever $r$ or $c$ is $0$ or $6$. In units $\tau/d = 1$, $V = \frac12\sum_r\sum_{c=0}^5(y_{r,c+1} - y_{r,c})^2 + \frac12\sum_c\sum_{r=0}^5(y_{r+1,c} - y_{r,c})^2$. For fixed $r$, the inner horizontal sum is the $N = 5$ string energy $\frac12\mathbf{y}_r^\top K_1\mathbf{y}_r$ of row $r$. Order the beads row by row, index $5(r - 1) + c$ with $c$ running fastest. Then the horizontal part is block diagonal with $K_1$ in each block, which is $I_5\otimes K_1$. The vertical part couples $(r, c)$ to $(r', c)$ with weight $(K_1)_{rr'}$, and $(K_1\otimes I_5)_{(r,c),(r',c')} = (K_1)_{rr'}\delta_{cc'}$ is exactly that. So $K = K_1\otimes I_5 + I_5\otimes K_1$.

(b) $(\mathbf{u}\otimes\mathbf{v})_{(r,c)} = u_rv_c$. Then $[(K_1\otimes I)(\mathbf{u}\otimes\mathbf{v})]_{(r,c)} = \sum_{r',c'}(K_1)_{rr'}\delta_{cc'}u_{r'}v_{c'} = (K_1\mathbf{u})_rv_c = s\,u_rv_c$, and $[(I\otimes K_1)(\mathbf{u}\otimes\mathbf{v})]_{(r,c)} = u_r(K_1\mathbf{v})_c = t\,u_rv_c$. Adding, $K(\mathbf{u}\otimes\mathbf{v}) = (s + t)\,\mathbf{u}\otimes\mathbf{v}$. With $M = I$ and $K_1$'s eigenpairs from (4.5) ($\omega_0 = 1$), $s_j = 4\sin^2\frac{j\pi}{12}$ for $j = 1..5$, i.e. $s = (2 - \sqrt3, 1, 2, 3, 2 + \sqrt3)$, and $\omega_{jl}^2 = s_j + s_l$. The shape $\sin\frac{jr\pi}{6}\sin\frac{lc\pi}{6}$ is the discrete drum mode, with $j - 1$ nodal lines one way and $l - 1$ the other.

(c) $(\mathbf{u}_j\otimes\mathbf{u}_l)\cdot(\mathbf{u}_{j'}\otimes\mathbf{u}_{l'}) = \sum_{r,c}u_{j,r}u_{l,c}u_{j',r}u_{l',c} = (\mathbf{u}_j\cdot\mathbf{u}_{j'})(\mathbf{u}_l\cdot\mathbf{u}_{l'}) = \delta_{jj'}\delta_{ll'}$ for orthonormal $\mathbf{u}$'s. That gives 25 orthonormal eigenvectors in $\mathbb{R}^{25}$: a basis, hence all the modes.

(d) There are 15 unordered pairs $\{j, l\}$. Pairs with $j \ne l$ give two modes, $(j, l)$ and $(l, j)$, at the same frequency. The reflection of the square across its diagonal swaps them, so these 2-fold degeneracies are symmetry-forced. Further coincidences come from $s_j + s_{6-j} = 4[\sin^2\frac{j\pi}{12} + \cos^2\frac{j\pi}{12}] = 4$, since $\sin\frac{(6 - j)\pi}{12} = \cos\frac{j\pi}{12}$. So $(1,5), (5,1), (2,4), (4,2), (3,3)$ all share $\omega^2 = 4$: a 5-fold degeneracy. Every other sum is distinct. The tally:
- 4 nondegenerate values: $(1,1), (2,2), (4,4), (5,5)$;
- 8 doubly degenerate values;
- 1 fivefold value.

That is $13$ distinct frequencies, and $4 + 16 + 5 = 25$ modes. ✓ The fivefold one is not explained by the square's symmetry group $D_4$, whose irreducible representations have dimension at most 2. It comes from a property of the chain, not of the square: flipping the sign on alternate beads ($u_k \to (-1)^ku_k$, the matrix $D = \operatorname{diag}((-1)^k)$) gives $DK_1D = 4I - K_1$, so $D$ maps the eigenvector for $s_j$ to one for $4 - s_j$. The chain spectrum is therefore symmetric about 2, $s_{6-j} = 4 - s_j$, and every pair $(j, 6 - j)$ lands on $\omega^2 = 4$. (The same $D\otimes D$ makes the whole lattice spectrum symmetric about 4.)

(e) The lowest six eigenvalues are $0.5359\ (1,1)$, $1.2679$ twice $\{(1,2), (2,1)\}$, $2\ (2,2)$, and $2.2679$ twice $\{(1,3), (3,1)\}$. So two degenerate pairs appear among them. `eigh` returns some rotation of each pair; you typically see diagonal-looking patterns instead of the product shapes. A 1% heavier bead breaks the symmetry and splits each pair; the shapes then lock to the perturbation.

**C2.** (a) Positive elongation means compression here; the sign does not matter for $V$. The front suspension spring is compressed by $z + a\theta - z_f$ (the body point above the front axle moves $z + a\theta$ for small nose-up $\theta$). The rear suspension is compressed by $z - b\theta - z_r$, and the tyres by $z_f$ and $z_r$. So
$$B = \begin{pmatrix}1 & a & -1 & 0\\1 & -b & 0 & -1\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{pmatrix}, \quad C = \operatorname{diag}(k_{sf}, k_{sr}, k_{tf}, k_{tr}), \quad M = \operatorname{diag}(1200, 1800, 90, 70),$$
$$K = B^\top CB = \begin{pmatrix}80000 & -3000 & -45000 & -35000\\-3000 & 133200 & -49500 & 52500\\-45000 & -49500 & 245000 & 0\\-35000 & 52500 & 0 & 235000\end{pmatrix}\ \text{(SI)}.$$
Here $K_{z\theta} = k_{sf}a - k_{sr}b = 49500 - 52500 = -3000$ couples bounce to pitch. It vanishes when $k_{sf}a = k_{sr}b$: a push at the center of mass then compresses the two springs so that their moments about the center of mass cancel.

(b) `normal_modes(K_car, M_car)` gives $f = \omega/2\pi = 1.170,\ 1.261,\ 8.323,\ 9.237$ Hz.

(c) Body modes, with the axles moving in between: each suspension spring is in series with its tyre, $k_{\text{eff},f} = \frac{45\cdot200}{245} = 36.735$ kN/m and $k_{\text{eff},r} = \frac{35\cdot200}{235} = 29.787$ kN/m. Bounce: $\omega^2 = (k_{\text{eff},f} + k_{\text{eff},r})/m_b = 55.435$, so $f = 1.185$ Hz. Pitch: $\omega^2 = (k_{\text{eff},f}a^2 + k_{\text{eff},r}b^2)/I = 61.928$, so $f = 1.252$ Hz. Wheel hop with the body fixed: front $\sqrt{245000/90}/2\pi = 8.304$ Hz, rear $\sqrt{235000/70}/2\pi = 9.222$ Hz. All four estimates are within 1.5% of the full result. The axles are 7–8 times faster than the body, so the two pairs barely interact. That separation of scales is why ride and wheel-hop can be designed separately.

(d) Scale each body mode so that its largest entry is 1. Mode 1 is $(z, \theta) = (1, 0.3155)$: the body point at $x$ ahead of the center of mass moves $z + x\theta$, which is zero at $x = -3.170$ m, behind the rear axle ($x = -1.5$). This is the **bounce** mode: the whole body rises, with the front moving $1.347$ and the rear $0.527$. Mode 2 is $(-0.4741, 1)$, with a node at $x = +0.474$ m, between the center of mass and the front axle. This is the **pitch** mode. The mixing is strong because the uncoupled bounce and pitch frequencies are close: $\omega^2 = 55.435$ and $61.928$, a gap of $6.49\ \mathrm{s^{-2}}$. Compare the mass-normalized coupling $(k_{\text{eff},f}a - k_{\text{eff},r}b)/\sqrt{m_bI} = -2.907\ \mathrm{s^{-2}}$, which is almost half the gap. For a symmetric $2\times2$ problem the mixing angle satisfies $\tan2\phi = 2|\text{coupling}|/\text{gap} = 0.90$, so $\phi \approx 21°$. Near-degenerate modes mix strongly under small couplings, the finite-gap version of §3.6, where an exact degeneracy lets any mixture through.

**C3.** (a) From $KA = MA\Omega^2$, $M^{-1}KA = A\Omega^2$, and by induction $(M^{-1}K)^lA = A\Omega^{2l}$. So $A^\top M(M^{-1}K)^lA = A^\top MA\,\Omega^{2l} = \Omega^{2l}$, and summing, $A^\top CA = \sum_lc_l\Omega^{2l} = p(\Omega^2)$, which is diagonal. For symmetry: $[M(M^{-1}K)^l]^\top = [(M^{-1}K)^l]^\top M = (KM^{-1})^lM$. Regrouping the product $KM^{-1}KM^{-1}\cdots KM^{-1}M$ as $M\cdot M^{-1}KM^{-1}K\cdots M^{-1}K$ gives $M(M^{-1}K)^l$. So each term, and $C$, is symmetric.

(b) We need $p(\omega_j^2) = 2\zeta_j\omega_j$ for $j = 1..n$. With $p(x) = \sum_{l=0}^{n-1}c_lx^l$ that is a linear system $V\mathbf{c} = \mathbf{y}$ whose Vandermonde matrix $V_{jl} = (\omega_j^2)^l$ is invertible iff the $\omega_j^2$ are distinct. So any $n$ target ratios can be met. Degree 1, $p(x) = c_0 + c_1x$, gives $C = c_0M + c_1K$: Rayleigh damping.

(c) The targets are $y_j = 2\zeta_j\omega_j = 2(0.02)(0.76537),\ 2(0.03)(1.41421),\ 2(0.05)(1.84776) = 0.030615,\ 0.084853,\ 0.184776$ at $x_j = 2 - \sqrt2, 2, 2 + \sqrt2$. The nodes are symmetric about 2, so write $p = d_0 + d_1u + d_2u^2$ with $u = x - 2 \in\{-\sqrt2, 0, \sqrt2\}$:
- $d_0 = y_2 = 0.084853$;
- $d_1 = (y_3 - y_1)/(2\sqrt2) = 0.054504$;
- $d_2 = (y_1 + y_3 - 2y_2)/4 = 0.011421$.

Expanding $u = x - 2$: $c_2 = d_2 = 0.011421$, $c_1 = d_1 - 4d_2 = 0.008819$, $c_0 = d_0 - 2d_1 + 4d_2 = 0.021529$. The script solves the Vandermonde system, gets the same coefficients, and checks that $A^\top CA = \operatorname{diag}(2\zeta_j\omega_j)$.

(d) First, the term $M(M^{-1}K)^2 = KM^{-1}K$ couples second neighbors, and higher powers fill $C$ in. In a finite-element model that destroys the sparsity the solver depends on. Second, $p(x)$ is controlled only at the fitted $\omega_j^2$. Beyond them it grows like $x^{n-1}$, or can turn negative, which gives the unmeasured high modes enormous or negative damping, and the Vandermonde system is badly conditioned when the $\omega_j^2$ are spread out. Damping is rarely known in more than two or three modes anyway.
