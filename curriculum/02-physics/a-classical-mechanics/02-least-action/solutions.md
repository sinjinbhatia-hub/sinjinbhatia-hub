# CM-2 Solutions

Conventions as in the notes. Every number below is printed or asserted in `project.py` (Parts 1–3 and 7).

## A. Warm-up

**A1.** (a) The true path is the straight line $x_{\text{cl}} = 1.5t$, so $\dot x = 1.5$ and $S = \frac12\cdot2\cdot1.5^2\cdot2 = 4.5$.

(b) $\dot x = 1.5t$, so $\int_0^2\dot x^2\,dt = 2.25\int_0^2t^2\,dt = 2.25\cdot\frac83 = 6$ and $S = \frac12\cdot2\cdot6 = 6$.

(c) $\eta = \frac34t^2 - 1.5t$, which vanishes at $t = 0$ and $t = 2$. $\dot\eta = 1.5(t - 1)$, so $\int_0^2\dot\eta^2\,dt = 2.25\int_0^2(t - 1)^2\,dt = 2.25\cdot\frac23 = 1.5$, and $\frac12 m\int\dot\eta^2 = \frac12\cdot2\cdot1.5 = 1.5$. Then $S_{\text{cl}} + 1.5 = 6$. ✓ The cross term is $\int_0^2m\dot x_{\text{cl}}\dot\eta\,dt = m\cdot1.5\,[\eta]_0^2 = 0$ because $\dot x_{\text{cl}}$ is constant and $\eta$ vanishes at both ends. That is the first variation vanishing on the true path (Lemma 3.1).

**A2.** The string fixes $y_2 = c - y$, so $\dot y_2 = -\dot y$. Then $T = \frac12 m_1\dot y^2 + \frac12 m_2\dot y^2 = \frac12(m_1 + m_2)\dot y^2$ and $V = m_1gy + m_2g(c - y) = (m_1 - m_2)gy + \text{const}$. Dropping the constant (Theorem 3.4 with $F = \text{const}\cdot t$),
$$L = \tfrac12(m_1 + m_2)\dot y^2 - (m_1 - m_2)gy.$$
$\partial L/\partial\dot y = (m_1 + m_2)\dot y$ and $\partial L/\partial y = -(m_1 - m_2)g$, so $(m_1 + m_2)\ddot y + (m_1 - m_2)g = 0$:
$$\ddot y = \frac{m_2 - m_1}{m_1 + m_2}\,g,$$
which is CM-1 (3.16). The tension is gone because $y$ parametrizes only configurations the string allows; the tension is the force that enforces the constraint, and it never enters $T$ or $V$. (Its total power is zero: with $T_s$ the tension, it delivers $T_s\dot y_1$ to one mass and $T_s\dot y_2 = -T_s\dot y_1$ to the other, so it does no net work.) CM-3 Example 4.2 recovers it with a multiplier. The project checks $\ddot y$ in SymPy.

**A3.** $L = T + V = \frac12 m\dot x^2 + \frac12 kx^2$. $\partial L/\partial\dot x = m\dot x$, $\partial L/\partial x = kx$, so $m\ddot x - kx = 0$, i.e. $\ddot x = (k/m)x$. With $x(0) = x_0$, $\dot x(0) = 0$: $x(t) = x_0\cosh(\sqrt{k/m}\,t)$. The mass runs away from equilibrium exponentially, as if the spring pushed outward with Hooke's law reversed. This is Theorem 3.3 with $\beta/\alpha = +1$.

$L = 3(T - V)$: $\partial L/\partial\dot x = 3m\dot x$, $\partial L/\partial x = -3kx$, so $\mathcal E = 3(m\ddot x + kx)$, which vanishes exactly when $m\ddot x + kx = 0$: the correct oscillator (Theorem 3.4(b) with $c = 3$).

**A4.** For each added term $G$ compute $\mathcal E[G] = \frac{d}{dt}\frac{\partial G}{\partial\dot q} - \frac{\partial G}{\partial q}$; the new equation is $m\ddot q + \mathcal E[G] = 0$.

(a) $G = 3q\dot q$: $\partial G/\partial\dot q = 3q$, whose time derivative is $3\dot q$; $\partial G/\partial q = 3\dot q$. $\mathcal E[G] = 0$: no change. $F = \frac32q^2$, since $\frac{d}{dt}\frac32q^2 = 3q\dot q$.

(b) $G = t\dot q$: $\partial G/\partial\dot q = t$, whose time derivative is $1$; $\partial G/\partial q = 0$. $\mathcal E[G] = 1$, so the equation becomes $m\ddot q = -1$: a constant force of size 1 in the $-q$ direction. The reason: $t\dot q = \frac{d}{dt}(tq) - q$, and $-q$ in $L$ is a potential $V = q$. The term is not $dF/dt$ for any $F(q, t)$.

(c) $G = \dot q\cos q$: $\partial G/\partial\dot q = \cos q$, whose time derivative is $-\sin q\,\dot q$; $\partial G/\partial q = -\dot q\sin q$. $\mathcal E[G] = 0$: no change. $F = \sin q$.

(d) $G = q + t\dot q$: $\partial G/\partial\dot q = t$, time derivative $1$; $\partial G/\partial q = 1$. $\mathcal E[G] = 0$: no change. $F = tq$, since $\frac{d}{dt}(tq) = q + t\dot q$.

The project's SymPy check returns $\{0, 1, 0, 0\}$ for (a)–(d).

**A5.** Differentiate slot by slot.

- $p_r = m\dot r$ (kg·m/s); $\mathcal F_r = mr\dot\theta^2 + mr\sin^2\theta\,\dot\phi^2 - V'(r)$.
- $p_\theta = mr^2\dot\theta$ (kg·m²/s); $\mathcal F_\theta = mr^2\sin\theta\cos\theta\,\dot\phi^2$.
- $p_\phi = mr^2\sin^2\theta\,\dot\phi$ (kg·m²/s); $\mathcal F_\phi = 0$.

$\phi$ is cyclic, so $p_\phi$ is conserved: it is $L_z$, the angular momentum about the $z$-axis. $p_\theta$ is not conserved, because $\mathcal F_\theta \ne 0$ whenever $\dot\phi \ne 0$, even though the total angular momentum is conserved: $|\mathbf L|^2 = m^2r^4(\dot\theta^2 + \sin^2\theta\,\dot\phi^2) = p_\theta^2 + p_\phi^2/\sin^2\theta$, and the two terms trade off as $\theta$ changes. Two of the three "momenta" are angular momenta, not mass times velocity.

## B. Core

**B1.** (a) $\partial L/\partial\dot q = f'(\dot q)$, whose time derivative is $f''(\dot q)\ddot q$; $\partial L/\partial q = 0$. So $\mathcal E = f''(\dot q)\ddot q = 0$, and since $f'' \ne 0$, $\ddot q = 0$.

(b) Now $\partial L/\partial q = -V'(q)$, so $\mathcal E = f''(\dot q)\ddot q + V'(q) = 0$. Suppose its solutions agree with Newton's for every $V$ and every initial condition. Take $V(q) = q$: Newton gives $\ddot q = -1/m$ along every solution, while E–L gives $f''(\dot q)\ddot q = -1$. Start a solution with arbitrary velocity $v_0$. At $t = 0$ both hold, so $f''(v_0)\cdot(-1/m) = -1$, i.e. $f''(v_0) = m$. Since $v_0$ was arbitrary, $f''(v) = m$ for all $v$, and integrating twice, $f(v) = \frac12 mv^2 + av + b$. Conversely, with this $f$, $f'' = m$ and E–L is Newton's law for every $V$.

(c) $a\dot q = \frac{d}{dt}(aq)$ and $b = \frac{d}{dt}(bt)$ are total derivatives of functions of $(q, t)$, so by Theorem 3.4 they change nothing. The kinetic term is fixed up to exactly the freedom Theorem 3.4 allows.

(d) Landau & Lifshitz require **Galilean relativity**: the equations of motion must be the same in every inertial frame. Under a boost $v \to v + u$ the Lagrangian may change by at most a total derivative. For $\frac12 mv^2$: $\frac12 m(v + u)^2 = \frac12 mv^2 + muv + \frac12 mu^2 = \frac12 mv^2 + \frac{d}{dt}\big(muq + \frac12 mu^2t\big)$, a total derivative. For $v^4$: $(v + u)^4 = v^4 + 4uv^3 + 6u^2v^2 + 4u^3v + u^4$, and $4uv^3$ is not $dF/dt$ for any $F(q, t)$, because $dF/dt = F_q\dot q + F_t$ is linear in $v$. So $v^4$ would give different physics in different inertial frames. (Homogeneity and isotropy of space and time are also used, to make $L$ depend on $|\mathbf v|^2$ alone.)

**B2.** (a) $\tilde L = \frac12 m(\dot r^2 + r^2\dot\theta^2) - \tilde V(r, \theta)$ with $\tilde V(r, \theta) = V(r\cos\theta, r\sin\theta)$. Then
$$\tilde{\mathcal E}_r = m\ddot r - mr\dot\theta^2 + \cos\theta\,V_x + \sin\theta\,V_y,\qquad \tilde{\mathcal E}_\theta = \frac{d}{dt}(mr^2\dot\theta) + \frac{\partial\tilde V}{\partial\theta} = m(r^2\ddot\theta + 2r\dot r\dot\theta) - r\sin\theta\,V_x + r\cos\theta\,V_y.$$

(b) $\dot x = \dot r\cos\theta - r\dot\theta\sin\theta$. Differentiating again,
$$\ddot x = \ddot r\cos\theta - 2\dot r\dot\theta\sin\theta - r\ddot\theta\sin\theta - r\dot\theta^2\cos\theta.$$
$\dot y = \dot r\sin\theta + r\dot\theta\cos\theta$, so
$$\ddot y = \ddot r\sin\theta + 2\dot r\dot\theta\cos\theta + r\ddot\theta\cos\theta - r\dot\theta^2\sin\theta.$$

(c) The Cartesian expressions are $\mathcal E_x = m\ddot x + V_x$ and $\mathcal E_y = m\ddot y + V_y$.

$\cos\theta\,\ddot x + \sin\theta\,\ddot y$: the $\ddot r$ terms give $\ddot r(\cos^2\theta + \sin^2\theta) = \ddot r$; the $\dot r\dot\theta$ terms give $2\dot r\dot\theta(-\sin\theta\cos\theta + \sin\theta\cos\theta) = 0$; the $\ddot\theta$ terms give $r\ddot\theta(-\sin\theta\cos\theta + \sin\theta\cos\theta) = 0$; the $\dot\theta^2$ terms give $-r\dot\theta^2$. So $\cos\theta\,\ddot x + \sin\theta\,\ddot y = \ddot r - r\dot\theta^2$, and
$$\cos\theta\,\mathcal E_x + \sin\theta\,\mathcal E_y = m(\ddot r - r\dot\theta^2) + \cos\theta\,V_x + \sin\theta\,V_y = \tilde{\mathcal E}_r.\ ✓$$

$-r\sin\theta\,\ddot x + r\cos\theta\,\ddot y$: the $\ddot r$ terms give $r\ddot r(-\sin\theta\cos\theta + \sin\theta\cos\theta) = 0$; the $\dot r\dot\theta$ terms give $2r\dot r\dot\theta(\sin^2\theta + \cos^2\theta) = 2r\dot r\dot\theta$; the $\ddot\theta$ terms give $r^2\ddot\theta(\sin^2\theta + \cos^2\theta) = r^2\ddot\theta$; the $\dot\theta^2$ terms give $r^2\dot\theta^2(\sin\theta\cos\theta - \sin\theta\cos\theta) = 0$. So
$$-r\sin\theta\,\mathcal E_x + r\cos\theta\,\mathcal E_y = m(r^2\ddot\theta + 2r\dot r\dot\theta) - r\sin\theta\,V_x + r\cos\theta\,V_y = \tilde{\mathcal E}_\theta.\ ✓$$

The coefficients $(\cos\theta, \sin\theta)$ and $(-r\sin\theta, r\cos\theta)$ are the columns of $J = \partial(x, y)/\partial(r, \theta)$, so this is $\tilde{\boldsymbol{\mathcal E}} = J^\top\boldsymbol{\mathcal E}$. The project repeats the check in SymPy.

**B3.** (a) Differentiate:
$\dot x = \dot X\cos\Omega t - \dot Y\sin\Omega t - \Omega(X\sin\Omega t + Y\cos\Omega t) = (\dot X - \Omega Y)\cos\Omega t - (\dot Y + \Omega X)\sin\Omega t$,
$\dot y = \dot X\sin\Omega t + \dot Y\cos\Omega t + \Omega(X\cos\Omega t - Y\sin\Omega t) = (\dot X - \Omega Y)\sin\Omega t + (\dot Y + \Omega X)\cos\Omega t$.
So $(\dot x, \dot y)$ is the vector $(\dot X - \Omega Y,\ \dot Y + \Omega X)$ rotated by the angle $\Omega t$. Rotation preserves length, so $\dot x^2 + \dot y^2 = (\dot X - \Omega Y)^2 + (\dot Y + \Omega X)^2$ and $\tilde L$ is as stated.

(b) $X$: $p_X = m(\dot X - \Omega Y)$, $\dot p_X = m(\ddot X - \Omega\dot Y)$, $\partial\tilde L/\partial X = m\Omega(\dot Y + \Omega X)$. E–L: $m(\ddot X - \Omega\dot Y) - m\Omega(\dot Y + \Omega X) = 0$, i.e.
$$m\ddot X = 2m\Omega\dot Y + m\Omega^2X.$$
$Y$: $p_Y = m(\dot Y + \Omega X)$, $\dot p_Y = m(\ddot Y + \Omega\dot X)$, $\partial\tilde L/\partial Y = -m\Omega(\dot X - \Omega Y)$. E–L: $m(\ddot Y + \Omega\dot X) + m\Omega(\dot X - \Omega Y) = 0$, i.e.
$$m\ddot Y = -2m\Omega\dot X + m\Omega^2Y.$$
The terms $m\Omega^2(X, Y)$ point away from the axis: the **centrifugal** force. The terms $2m\Omega(\dot Y, -\dot X)$ equal $-2m\,\boldsymbol\Omega\times\mathbf v_{\text{rot}}$ with $\boldsymbol\Omega = \Omega\hat{\mathbf z}$ and $\mathbf v_{\text{rot}} = (\dot X, \dot Y, 0)$, since $\hat{\mathbf z}\times(\dot X, \dot Y, 0) = (-\dot Y, \dot X, 0)$: the **Coriolis** force.

(c) $p_X = m(\dot X - \Omega Y)$. By (a), the inertial velocity is the rotation of $(\dot X - \Omega Y, \dot Y + \Omega X)$, so $(p_X, p_Y)$ is the inertial momentum $(m\dot x, m\dot y)$ resolved along the rotating axes. $m\dot X$ is the momentum *relative to* the rotating frame, which is a different vector. This is the rule $\tilde p_j = \sum_iJ_{ij}p_i$ of §3.6 with $J$ the rotation matrix.

(d) Theorem 3.5, which allows time-dependent point transformations. The map $(X, Y) \mapsto (x, y)$ is invertible for every $t$, so solutions of these equations are exactly the inertial solutions (straight lines at constant speed) written in rotating coordinates.

**B4.** (a) With $m = k = 1$, (3.11) is $\mathcal Q[\eta] = \frac\tau2\sum_nb_n^2\big[(n\pi/\tau)^2 - 1\big]$. The form is diagonal in the sine basis, so the descending directions are spanned by the $\sin(n\pi t/\tau)$ with negative brackets, i.e. with $n\pi/\tau < 1$, i.e. $n\pi < \tau$. Let $K$ be the number of such $n$. On their span, $\mathcal Q < 0$ for every nonzero $\eta$, so there are at least $K$ independent descending directions. There cannot be more: if $W$ were a subspace of dimension $K + 1$ on which $\mathcal Q < 0$, the $K$ linear conditions $b_1 = \dots = b_K = 0$ would leave a nonzero $\eta \in W$ with only $n > K$ modes, and for that $\eta$ every term of (3.11) is $\ge 0$, a contradiction. So the count is exactly $K$, and those $n$ are exactly the indices of the points $n\pi$ in $(0, \tau)$. These points are the conjugate points of $t = 0$: the solution $\sin t$ of $\ddot q = -q$ starts at 0 and vanishes again at each $n\pi$.

(b) $\tau = 2$: $\pi > 2$, so 0. $\tau = 4$: $\pi < 4 < 2\pi$, so 1. $\tau = 5$: $\pi < 5 < 2\pi \approx 6.28$, so 1. $\tau = 7$: $2\pi < 7 < 3\pi \approx 9.42$, so 2. The project's discrete Hessian has $\{0, 1, 1, 2\}$ negative eigenvalues for $\tau = 2, 4, 5, 7$. ✓

(c) The general solution is $q = A\cos t + B\sin t$; $q(0) = 0$ forces $A = 0$, and then $q(\pi) = B\sin\pi = 0$ for every $B$. So $q = B\sin t$ for any $B$: infinitely many true paths with the same endpoints. Their action is $\frac12B^2\int_0^\pi(\cos^2t - \sin^2t)\,dt = \frac12B^2\int_0^\pi\cos2t\,dt = \frac12B^2\big[\tfrac12\sin2t\big]_0^\pi = 0$ for every $B$. The difference between two members is $(B - B')\sin t$, which is the $n = 1$ mode $\sin(\pi t/\tau)$ with $\tau = \pi$; its bracket $(\pi/\pi)^2 - 1$ is zero, so $S$ is exactly flat along the family.

**B5.** (a) Integrate $\dot q^2$ by parts: $\int_0^\tau\dot q^2\,dt = [q\dot q]_0^\tau - \int_0^\tau q\ddot q\,dt$. Since $\ddot q = -\omega^2q$ on the true path,
$$S_{\text{cl}} = \frac m2\int_0^\tau(\dot q^2 - \omega^2q^2)\,dt = \frac m2\Big([q\dot q]_0^\tau - \int_0^\tau q(\ddot q + \omega^2q)\,dt\Big) = \frac m2\big(q_2\dot q(\tau) - q_1\dot q(0)\big).$$
From the given path, $\dot q(t) = \omega[-q_1\cos\omega(\tau - t) + q_2\cos\omega t]/\sin\omega\tau$, so
$\dot q(\tau) = \omega(q_2\cos\omega\tau - q_1)/\sin\omega\tau$ and $\dot q(0) = \omega(q_2 - q_1\cos\omega\tau)/\sin\omega\tau$. Substituting,
$$S_{\text{cl}} = \frac{m\omega}{2\sin\omega\tau}\Big[q_2^2\cos\omega\tau - q_1q_2 - q_1q_2 + q_1^2\cos\omega\tau\Big] = \frac{m\omega\big[(q_1^2 + q_2^2)\cos\omega\tau - 2q_1q_2\big]}{2\sin\omega\tau}.$$

(b) $\partial S_{\text{cl}}/\partial q_2 = m\omega(2q_2\cos\omega\tau - 2q_1)/(2\sin\omega\tau) = m\omega(q_2\cos\omega\tau - q_1)/\sin\omega\tau = m\dot q(\tau)$. ✓
$\partial S_{\text{cl}}/\partial q_1 = m\omega(q_1\cos\omega\tau - q_2)/\sin\omega\tau = -m\dot q(0)$. ✓

(c) Write $S_{\text{cl}} = \frac{m\omega}2\big[(q_1^2 + q_2^2)\cot\omega\tau - 2q_1q_2\csc\omega\tau\big]$. Then
$$\frac{\partial S_{\text{cl}}}{\partial\tau} = \frac{m\omega^2}2\Big[-\frac{q_1^2 + q_2^2}{\sin^2\omega\tau} + \frac{2q_1q_2\cos\omega\tau}{\sin^2\omega\tau}\Big] = -\frac{m\omega^2\,(q_1^2 + q_2^2 - 2q_1q_2\cos\omega\tau)}{2\sin^2\omega\tau}.$$
The energy, evaluated at $t = 0$: $E = \frac12 m\dot q(0)^2 + \frac12 m\omega^2q_1^2 = \frac{m\omega^2}{2\sin^2\omega\tau}\big[(q_2 - q_1\cos\omega\tau)^2 + q_1^2\sin^2\omega\tau\big] = \frac{m\omega^2}{2\sin^2\omega\tau}\big[q_2^2 - 2q_1q_2\cos\omega\tau + q_1^2\big]$. So $\partial S_{\text{cl}}/\partial\tau = -E$, as (3.13) predicts.

Numbers ($m = \omega = 1$, $q_1 = 0$, $q_2 = 1$, $\tau = 2$): $S_{\text{cl}} = \frac{\cos2}{2\sin2} = \frac12\cot2 = -0.22883$; $p(0) = 1/\sin2 = 1.09975$; $p(\tau) = \cos2/\sin2 = -0.45766$; $\partial S_{\text{cl}}/\partial\tau = -1/(2\sin^22) = -0.60473$, and $E = \frac12p(0)^2 = 0.60473$. The project checks all three derivatives by finite differences.

(d) As $\tau \to \pi/\omega$, $\sin\omega\tau \to 0$ and $\cos\omega\tau \to -1$, so the bracket tends to $-(q_1^2 + q_2^2) - 2q_1q_2 = -(q_1 + q_2)^2$. Unless $q_2 = -q_1$, $S_{\text{cl}}$ diverges (to $-\infty$ from below $\pi/\omega$, $+\infty$ from above), and so does the required launch momentum $p(0)$. At $\tau = \pi/\omega$ the boundary-value problem has no solution if $q_2 \ne -q_1$ and infinitely many if $q_2 = -q_1$ (B4(c)). This is the conjugate point where Theorem 3.6's $n = 1$ bracket changes sign: the principal function is singular exactly where the true path stops being a minimum.

**B6.** (a) Fields: $\mathbf E' = -\nabla(\Phi - \partial_t\chi) - \partial_t(\mathbf A + \nabla\chi) = \mathbf E + \nabla\partial_t\chi - \partial_t\nabla\chi = \mathbf E$, since mixed partials commute. $\mathbf B' = \nabla\times\mathbf A + \nabla\times\nabla\chi$, and the curl of a gradient vanishes component by component: $x$: $\partial_y\partial_z\chi - \partial_z\partial_y\chi = 0$; $y$: $\partial_z\partial_x\chi - \partial_x\partial_z\chi = 0$; $z$: $\partial_x\partial_y\chi - \partial_y\partial_x\chi = 0$. So $\mathbf B' = \mathbf B$.
Lagrangian: $L' - L = -q_e(\Phi' - \Phi) + q_e\dot{\mathbf r}\cdot(\mathbf A' - \mathbf A) = q_e\partial_t\chi + q_e\dot{\mathbf r}\cdot\nabla\chi = q_e\frac{d\chi}{dt}$, by the chain rule for $\chi(\mathbf r(t), t)$. By Theorem 3.4 the equations of motion are unchanged.

(b) $\mathbf A_L - \mathbf A_s = \big(0 + \frac B2y,\ Bx - \frac B2x,\ 0\big) = \frac B2(y, x, 0)$. With $\chi = \frac12Bxy$: $\nabla\chi = (\frac12By, \frac12Bx, 0)$. ✓

(c) Symmetric gauge, $L_s = \frac12 m|\mathbf v|^2 + \frac12q_eB(x\dot y - y\dot x)$: $p_x = m\dot x - \frac12q_eBy$, $p_y = m\dot y + \frac12q_eBx$.
Landau gauge, $L_L = \frac12 m|\mathbf v|^2 + q_eBx\dot y$: $p_x = m\dot x$, $p_y = m\dot y + q_eBx$.
In the Landau gauge $L_L$ contains neither $y$ nor $z$, so $p_y = m\dot y + q_eBx$ and $p_z = m\dot z$ are conserved. (Check: the E–L equation $m\ddot y = -q_eB\dot x$ says exactly $\frac{d}{dt}(m\dot y + q_eBx) = 0$.) $p_x = m\dot x$ is not conserved, because $\partial L_L/\partial x = q_eB\dot y \ne 0$. In the symmetric gauge both $x$ and $y$ appear, so neither $p_x$ nor $p_y$ is conserved; for example $\dot p_y = \partial L_s/\partial y = -\frac12q_eB\dot x$. What the symmetric gauge conserves instead is rotational: in polar coordinates $L_s = \frac12 m(\dot r^2 + r^2\dot\phi^2 + \dot z^2) + \frac12q_eBr^2\dot\phi$, $\phi$ is cyclic, and $p_\phi = mr^2\dot\phi + \frac12q_eBr^2$ is conserved. Between the two gauges, $\mathbf p_L - \mathbf p_s = q_e\nabla\chi$. The motion is gauge-independent: both Lagrangians give $m\ddot x = q_eB\dot y$, $m\ddot y = -q_eB\dot x$, $m\ddot z = 0$ (project, Part 4).

**B7.** (a) $\frac{dF}{dt} = \frac{\partial F}{\partial q}\dot q + \frac{\partial F}{\partial t} = 2q\dot q\sin t + q^2\cos t$.

(b) Prediction: $S' - S = F(q(2), 2) - F(q(0), 0) = 1^2\sin2 - 0 = \sin2 = 0.909297$, the same for every path with these endpoints. Code (the reference is in `project.py`, Part 7):

```python
Lp = lambda q, qd, t: 0.5 * qd**2 - 0.5 * q**2 + 2 * q * qd * np.sin(t) + q**2 * np.cos(t)
paths = [(lambda t: np.sin(t) / np.sin(2), lambda t: np.cos(t) / np.sin(2)),
         (lambda t: t / 2, lambda t: 0.5 + 0 * t),
         (lambda t: t**2 / 4, lambda t: t / 2)]
for qf, qdf in paths:
    S = action(qf, qdf, 0, 2, 1.0, lambda q: 0.5 * q**2)
    Sp = gl_integrate(lambda t: Lp(qf(t), qdf(t), t), 0, 2)
    assert np.isclose(Sp - S, np.sin(2))
```

All three differences equal $0.909297$.

(c) With the telescoping discretization, $S_h' = S_h + \sum_{k=0}^{N}\big[F(q_{k+1}, t_{k+1}) - F(q_k, t_k)\big] = S_h + F(q_{N+1}, t_{N+1}) - F(q_0, t_0)$, and both end values are fixed. So $S_h'$ and $S_h$ differ by a constant, have identical gradients, and have the same minimizer node for node: Theorem 3.4 survives discretization exactly. With a left-point rule, $h\big[2q_k\frac{q_{k+1} - q_k}h\sin t_k + q_k^2\cos t_k\big]$ differs from $F(q_{k+1}, t_{k+1}) - F(q_k, t_k)$ by terms of order $h^2$ that depend on the path, and $N \approx \tau/h$ of them add to a path-dependent error of order $h$. The added term is then no longer constant on paths, its gradient is a small spurious force, and the minimizer shifts by an amount that vanishes only as $h \to 0$.

## C. Challenge

**C1.** (a) $p = \dfrac{\partial L}{\partial y'} = \dfrac{1}{\sqrt{2gy}}\cdot\dfrac{y'}{\sqrt{1 + y'^2}}$, and
$$h = py' - L = \frac{y'^2}{\sqrt{2gy(1 + y'^2)}} - \frac{1 + y'^2}{\sqrt{2gy(1 + y'^2)}} = -\frac{1}{\sqrt{2gy(1 + y'^2)}}.$$

(b) $h$ is conserved. $L$ does not depend explicitly on the "time" $x$, and for such Lagrangians the energy function is constant along solutions: this is F6's Beltrami identity (F6 (3.11)), the same computation that gives energy conservation in mechanics. $p$ is not conserved, because $L$ depends on the "coordinate" $y$: the E–L equation says $dp/dx = \partial L/\partial y \ne 0$.

(c) With $y' = \tan\alpha$ and $|\alpha| < \pi/2$, $\sqrt{1 + y'^2} = 1/\cos\alpha$, so $h = -\cos\alpha/\sqrt{2gy} = -\cos\alpha/v$. With $\beta = \pi/2 - \alpha$, the angle from the vertical, $\cos\alpha = \sin\beta$, and conservation of $h$ says $\sin\beta/v = \text{const}$. This is **Snell's law** $n\sin\beta = \text{const}$ for a medium whose refractive index is proportional to $1/v$ (F6 Example 4.3). Johann Bernoulli's own 1697 solution used exactly this: he treated the bead as a light ray in a medium whose speed $\sqrt{2gy}$ increases with depth.

(d) $h = \text{const}$ gives $2gy(1 + y'^2) = 1/h^2$, i.e. $y(1 + y'^2) = 1/(2gh^2) =: 2a$. For the cycloid, $dx/d\theta = a(1 - \cos\theta)$ and $dy/d\theta = a\sin\theta$, so $y' = \sin\theta/(1 - \cos\theta)$ and
$$1 + y'^2 = \frac{(1 - \cos\theta)^2 + \sin^2\theta}{(1 - \cos\theta)^2} = \frac{2(1 - \cos\theta)}{(1 - \cos\theta)^2} = \frac{2}{1 - \cos\theta}.$$
Then $y(1 + y'^2) = a(1 - \cos\theta)\cdot\frac2{1 - \cos\theta} = 2a$. ✓

**C2.** (a) On the circle, $T = \frac12r^2\dot\theta^2 = \frac12$ and $V = -1$, so $L = T - V = \frac32$ and $S = \frac32\cdot2\pi = 3\pi = 9.424778$.

(b) At launch $r = 1$ and $v = 1$ whatever the direction, so $E = \frac12v^2 - 1/r = -\frac12$. With $E = -k/(2a)$ this gives $a = 1$ for every $\gamma$, hence $P = 2\pi a^{3/2} = 2\pi$: after $t = 2\pi$ the particle has completed exactly one orbit and is back at $(1, 0)$. The circle's velocity is perpendicular to $\mathbf r$; tilting it by $\gamma$ leaves a perpendicular component $\cos\gamma$, so $\ell = rv\cos\gamma = \cos\gamma$. From $\ell^2 = a(1 - e^2) = 1 - e^2$, $e^2 = 1 - \cos^2\gamma = \sin^2\gamma$, so $e = \sin\gamma$. For $\gamma = 0.3$, $e = 0.2955$; the project measures $e = (r_{\max} - r_{\min})/(r_{\max} + r_{\min})$ on the integrated orbit and gets the same value, with $a = \frac12(r_{\max} + r_{\min}) = 1$.

(c) Let $G = \mathbf r\cdot\mathbf p$ with $\mathbf p = m\dot{\mathbf r}$. Then $\dot G = \dot{\mathbf r}\cdot\mathbf p + \mathbf r\cdot\dot{\mathbf p} = m|\dot{\mathbf r}|^2 + \mathbf r\cdot\mathbf F = 2T + \mathbf r\cdot\big(-k\mathbf r/r^3\big) = 2T - k/r = 2T + V$. Average over one period $P$: $\langle\dot G\rangle = (G(P) - G(0))/P = 0$ because the motion is periodic. So $2\langle T\rangle + \langle V\rangle = 0$. Energy is constant, so $E = \langle T\rangle + \langle V\rangle = \frac12\langle V\rangle$, giving $\langle V\rangle = 2E = -k/a$ and $\langle T\rangle = k/(2a)$. Then
$$S = \int_0^P(T - V)\,dt = P\big(\langle T\rangle - \langle V\rangle\big) = P\cdot\frac{3k}{2a} = 2\pi\cdot\frac32 = 3\pi$$
for every orbit in (b), independent of $e$. The project integrates the $\gamma = 0.3$ ellipse, finds it back at $(1, 0)$ at $t = 2\pi$, and gets $S = 9.424778$, the same as the circle.

(d) With $s = \sin(t/2)$, $r = 1 + \varepsilon s$, $\dot r = \frac\varepsilon2\cos(t/2)$, $\dot\theta = 1$:
$\frac12\dot r^2 = \frac{\varepsilon^2}8\cos^2(t/2)$; $\frac12r^2\dot\theta^2 = \frac12 + \varepsilon s + \frac12\varepsilon^2s^2$; $-V = 1/r = 1 - \varepsilon s + \varepsilon^2s^2 + O(\varepsilon^3)$.
Sum: $L = \frac32 + \varepsilon^2\big[\frac18\cos^2(t/2) + \frac32\sin^2(t/2)\big] + O(\varepsilon^3)$. The $O(\varepsilon)$ terms cancel pointwise. On a general true path only their integral must vanish (the first variation); here they cancel at every instant because on the circle both $\partial L/\partial r = r\dot\theta^2 - 1/r^2$ and $p_r = \dot r$ are zero. Over $[0, 2\pi]$, $\int\cos^2(t/2)\,dt = \int\sin^2(t/2)\,dt = \pi$, so
$$S = 3\pi + \varepsilon^2\Big(\frac\pi8 + \frac{3\pi}2\Big) + O(\varepsilon^3) = 3\pi + \frac{13\pi}8\varepsilon^2 + O(\varepsilon^3),$$
with $13\pi/8 = 5.1051$; the project's second difference gives $5.1051$.

(e) Along the radial perturbation the action rises, but along the family of ellipses in (b), all true paths with the same endpoints and the same duration, the action is exactly constant. So the circle over one full period is **not a strict minimum**: there is a direction of zero second variation. A one-parameter family of true paths joining the same two events means $t = 2\pi$ is a **conjugate point** of $t = 0$ along the circle, the Kepler analogue of $\tau = \pi/\omega$ for the oscillator (B4(c)). By Jacobi's condition (F6 §3.11), for any interval extending past $t = 2\pi$ the circle is no longer a minimum at all.

**C3.** (a) $\partial L/\partial\dot q = e^{\gamma t}m\dot q$, with time derivative $e^{\gamma t}(m\ddot q + m\gamma\dot q)$; $\partial L/\partial q = -e^{\gamma t}kq$. E–L: $e^{\gamma t}(m\ddot q + m\gamma\dot q + kq) = 0$, and $e^{\gamma t} > 0$, so $m\ddot q + m\gamma\dot q + kq = 0$. (Project, Part 7.)

(b) $p = e^{\gamma t}m\dot q$. It is not the physical momentum $m\dot q$; it grows by the factor $e^{\gamma t}$.

(c) $h = p\dot q - L = e^{\gamma t}\big(m\dot q^2 - \frac12 m\dot q^2 + \frac12kq^2\big) = e^{\gamma t}\big(\frac12 m\dot q^2 + \frac12kq^2\big) = e^{\gamma t}E_{\text{mech}}$. Differentiate along a solution, using $\dot E_{\text{mech}} = \dot q(m\ddot q + kq) = -m\gamma\dot q^2$:
$$\frac{dh}{dt} = \gamma e^{\gamma t}E_{\text{mech}} - m\gamma e^{\gamma t}\dot q^2 = \gamma e^{\gamma t}\big(\tfrac12kq^2 - \tfrac12 m\dot q^2\big) = -\gamma L = -\frac{\partial L}{\partial t}.$$
So $h$ is not conserved (because $L$ depends explicitly on $t$), and it is not the mechanical energy, which is $e^{-\gamma t}h$.

(d) The E–L equation for $q_1$ is $m(a_1\ddot q_1 + \dot a_1\dot q_1) + b\,\partial V/\partial q_1 = 0$. Divide by $a_1$: $m\ddot q_1 + m(\dot a_1/a_1)\dot q_1 + (b/a_1)\partial V/\partial q_1 = 0$. For this to be $m\ddot q_1 + m\gamma_1\dot q_1 + \partial V/\partial q_1 = 0$ for all motions, with $\partial V/\partial q_1 \not\equiv 0$, we need $\dot a_1/a_1 = \gamma_1$ and $b = a_1$, so $b(t) = a_1(t) = C_1e^{\gamma_1t}$. For $q_2$ the E–L equation is $m(a_2\ddot q_2 + \dot a_2\dot q_2) + b\,\partial V/\partial q_2 = 0$; dividing by $a_2$ and matching with $m\ddot q_2 + m\gamma_2\dot q_2 + \partial V/\partial q_2 = 0$ (using $\partial V/\partial q_2 \not\equiv 0$) gives $\dot a_2/a_2 = \gamma_2$ and $b = a_2$, so $b(t) = C_2e^{\gamma_2t}$. Then $C_1e^{\gamma_1t} = C_2e^{\gamma_2t}$ for all $t$, which forces $\gamma_1 = \gamma_2$. So within this family there is no Lagrangian for coupled coordinates with different damping rates. The trick in (a) works because one overall factor $e^{\gamma t}$ can multiply the whole Lagrangian only when every coordinate is damped at the same rate. This is the sense of §3.11: there is no general Lagrangian for dissipation, and the special cases that exist produce a "momentum" and an "energy" that are not the physical ones. Rayleigh's dissipation function (CM-3 §3.9) is the general tool.
