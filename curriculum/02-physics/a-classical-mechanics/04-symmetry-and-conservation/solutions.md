# CM-4 Solutions

Conventions as in the notes and exercises: $p_i = \partial L/\partial\dot q_i$; $\delta L$ from (3.3); strict charge $Q = \sum_i p_iK_i$ (3.4); quasi-symmetry $\delta L = dF/dt$ with charge $Q = \sum_i p_iK_i - F$ (3.5); $h$ from (3.8). "On shell" means "using the Euler–Lagrange equations". Every symbolic claim below is re-derived, and every number printed, in the final section of `project.py` ("Reference checks for the exercises").

One fact is used over and over, so we prove it once.

**Lemma (when is $\delta L$ a total derivative?).** Let $G(q, \dot q, t)$ be given. There is an $F(q, t)$ with $G = dF/dt = \sum_j (\partial F/\partial q_j)\dot q_j + \partial F/\partial t$ only if $G$ is linear in the velocities, $G = \sum_j a_j(q, t)\dot q_j + b(q, t)$, with $a_j = \partial F/\partial q_j$ and $b = \partial F/\partial t$. In particular, if $G$ contains no velocities, then every $\partial F/\partial q_j = 0$, so $F = F(t)$ and $G = F'(t)$ must be a function of $t$ alone.

*Proof.* $dF/dt$ is affine in $\dot q$ with the stated coefficients, and two expressions that agree for all $(q, \dot q, t)$ have equal coefficients of each $\dot q_j$ and equal $\dot q$-free parts. ∎

## A. Warm-up

**A1.** (a) $\phi$ is cyclic: $p_\phi = m\ell^2\sin^2\theta\,\dot\phi$ is conserved. The symmetry is rotation about the vertical axis; gravity points along that axis, so it does not care about the azimuth. $\theta$ is not cyclic. (There is no explicit $t$ either, so $h = T + V$ is conserved by Theorems 3.6–3.7; that comes from time translation, not from a cyclic coordinate.)

(b) $x$ and $y$ are cyclic: $p_x = m\dot x$, $p_y = m\dot y$. Symmetries: horizontal translations. $z$ is not cyclic because gravity singles out a height. (Rotation about $\hat{\mathbf z}$ is also a symmetry, but no Cartesian coordinate is its angle; B1 finds its charge.)

(c) Nothing is cyclic, because $V$ contains $r$, $\theta$ and $z$. Try $K = (0, 1, a)$: rotate by $\varepsilon$ and advance by $a\varepsilon$, a screw motion. Then $z - a\theta \to (z + a\varepsilon) - a(\theta + \varepsilon) = z - a\theta$, so $V$ is unchanged; $K$ is constant, so $\dot K = 0$ and the kinetic term is untouched. Explicitly, with $V_2$ the derivative of $V$ in its second slot, $\partial L/\partial\theta = aV_2$ and $\partial L/\partial z = -V_2$, so
$$\delta L = 1\cdot aV_2 + a\cdot(-V_2) = 0.$$
Charge: $Q = p_\theta + a\,p_z = mr^2\dot\theta + ma\dot z$. Check by a coordinate change: with $\zeta = z - a\theta$, $L = \tfrac12 m\big(\dot r^2 + r^2\dot\theta^2 + (\dot\zeta + a\dot\theta)^2\big) - V(r, \zeta)$, $\theta$ is now cyclic, and $p_\theta = mr^2\dot\theta + ma(\dot\zeta + a\dot\theta) = mr^2\dot\theta + ma\dot z$. Same $Q$. A symmetry is cyclic in *some* coordinates; Noether's theorem saves you from having to find them.

(d) No coordinate is cyclic: $|\mathbf r|$ contains $x$, $y$ and $z$. Conserved anyway: $\mathbf L = \mathbf r\times m\dot{\mathbf r}$, all three components (rotations, §3.3b and Example 4.1); $h = \tfrac12 m|\dot{\mathbf r}|^2 - \kappa/r = E$ (no explicit $t$, Cartesian coordinates are fixed); and the Laplace–Runge–Lenz vector $\mathbf A$ (§3.8, no point symmetry behind it).

**A2.** (a) $p_r = m\dot r$, $p_\theta = mr^2\dot\theta$. $h = \dot r\,m\dot r + \dot\theta\,mr^2\dot\theta - L = m\dot r^2 + mr^2\dot\theta^2 - \tfrac12 m(\dot r^2 + r^2\dot\theta^2) + V = \tfrac12 m(\dot r^2 + r^2\dot\theta^2) + V(r)$. Conserved ($\partial L/\partial t = 0$) and equal to $T + V$ (polar coordinates do not depend on $t$, Theorem 3.7).

(b) $p = m\dot x$, $h = m\dot x^2 - L = \tfrac12 m\dot x^2 + \tfrac12 kx^2 - xF_0\cos\Omega t$. This is $T + V$ with the time-dependent potential $V = \tfrac12 kx^2 - xF_0\cos\Omega t$ (Theorem 3.7 allows $V(q, t)$). It is not conserved: by (3.9),
$$\frac{dh}{dt} = -\frac{\partial L}{\partial t} = -\frac{\partial}{\partial t}\big(xF_0\cos\Omega t\big) = xF_0\Omega\sin\Omega t.$$
That is the rate at which the drive changes the potential under the particle, $\partial V/\partial t$.

(c) $p_x = m\dot x$, $p_y = m\dot y + q_eBx$.
$$h = \dot x(m\dot x) + \dot y(m\dot y + q_eBx) - \tfrac12 m(\dot x^2 + \dot y^2) - q_eBx\dot y = \tfrac12 m(\dot x^2 + \dot y^2).$$
The magnetic term is linear in the velocities, so by Euler's theorem (3.10) with degree 1 it contributes itself to $\sum\dot q_i\,\partial L/\partial\dot q_i$ and then cancels against $-L$. $h$ is conserved (no explicit $t$) and equals the kinetic energy: the speed is constant, so the magnetic field does no work. $h$ is $T + V$ with $V = 0$; the velocity-dependent "potential" $U = -q_eBx\dot y$ of CM-3 does not appear in it.

**A3.** (a) $V = -\kappa(x^2 + y^2 + z^2)^{-1/2}$: $\partial V/\partial x = \kappa x/r^3$, and the same with $y$, $z$. So $x\,\partial_xV + y\,\partial_yV + z\,\partial_zV = \kappa(x^2 + y^2 + z^2)/r^3 = \kappa/r = (-1)\cdot V$. ✓

$V = \tfrac12 k(x^2 + y^2)$: $x(kx) + y(ky) = k(x^2 + y^2) = 2V$. ✓

$f = x^2y/z$: $x\,\partial_x f = 2x^2y/z$, $y\,\partial_y f = x^2y/z$, $z\,\partial_z f = z(-x^2y/z^2) = -x^2y/z$. Sum: $2x^2y/z = 2f$. ✓

(b) $f(\lambda x, \lambda y) = \lambda^2x^2 + \lambda y$. If this equalled $\lambda^k(x^2 + y)$ for all $\lambda > 0$, then at $(x, y) = (1, 0)$ we would need $\lambda^2 = \lambda^k$, so $k = 2$, and at $(0, 1)$ we would need $\lambda = \lambda^k$, so $k = 1$. Contradiction. Euler's identity: $x\cdot 2x + y\cdot 1 - kf = (2 - k)x^2 + (1 - k)y$. For $k = 1$ the defect is $x^2$, for $k = 2$ it is $-y$, and no $k$ kills both terms.

(c) $\partial T/\partial\dot q_i = \tfrac12\sum_k M_{ik}\dot q_k + \tfrac12\sum_j M_{ji}\dot q_j = \sum_k M_{ik}\dot q_k$, using $M_{ji} = M_{ij}$. Then $\sum_i\dot q_i\,\partial T/\partial\dot q_i = \sum_{i,k}M_{ik}\dot q_i\dot q_k = 2T$. ∎

**A4.** From Theorem 3.8 the exponents are: time $1 - k/2$, speed $k/2$, energy $k$, angular momentum $1 + k/2$.

| $k$ | time | speed | energy | ang. mom. |
|---|---|---|---|---|
| $-2$ | $2$ | $-1$ | $-2$ | $0$ |
| $-1$ | $3/2$ | $-1/2$ | $-1$ | $1/2$ |
| $1$ | $1/2$ | $1/2$ | $1$ | $3/2$ |
| $2$ | $0$ | $1$ | $2$ | $2$ |
| $4$ | $-1$ | $2$ | $4$ | $3$ |

(a) $V = mgz$ has $k = 1$, so period $\propto(\text{length})^{1/2}$. The string length is a constraint; scaling every length by $\alpha$ maps a motion on a string of length $\ell$ to a motion on a string of length $\alpha\ell$ with the same angular amplitude. So $T\propto\sqrt\ell$, as in $2\pi\sqrt{\ell/g}$.

(b) $k = -1$, $a'/a = 4$: period $\times 4^{3/2} = 8$, speed $\times 4^{-1/2} = 0.5$, energy $\times 4^{-1} = 0.25$, angular momentum $\times 4^{1/2} = 2$. Checks: $E = -\kappa/2a$ gives $E'/E = a/a' = 1/4$ ✓; $L^2 = m\kappa a(1 - e^2)$ at fixed $e$ gives $L'/L = \sqrt4 = 2$ ✓.

(c) $k = 4$: period $\propto A^{-1}$, so doubling the amplitude halves the period. (B6 gets $T(1) = 7.416299$, $T(2) = 3.708149$.)

**A5.** $\partial L/\partial x = 0$, $\partial L/\partial y = -V'(y)$, $\partial L/\partial\dot x = m\dot x$, $\partial L/\partial\dot y = m\dot y$.

(a) $K = (1, 0)$, $\dot K = 0$: $\delta L = 0$. Symmetry, $Q = m\dot x$.

(b) $K = (0, 1)$, $\dot K = 0$: $\delta L = -V'(y)$. Not a symmetry unless $V' \equiv 0$. By the Lemma, a velocity-free $\delta L$ is a total derivative only if it depends on $t$ alone, so we need $V'(y) = $ const: $V = cy$ (a uniform field). Then $F = -ct$ and $Q = m\dot y + ct$ (B1(b) with $c = mg$). For any other $V$, neither.

(c) $K = (-y, x)$, $\dot K = (-\dot y, \dot x)$: $\delta L = 0\cdot(-y) - V'(y)\,x + m\dot x(-\dot y) + m\dot y\,\dot x = -xV'(y)$. This depends on $x$, so by the Lemma it is not a total derivative unless $V' \equiv 0$. Neither.

(d) $K = (y, 0)$, $\dot K = (\dot y, 0)$: $\delta L = m\dot x\dot y$. Quadratic in velocities, so by the Lemma no $F(q, t)$ exists. Neither, for every $V$. (A shear changes the kinetic energy's shape; no total derivative can compensate.)

**A6.** (a) $K = \dot x$, $\dot K = \ddot x$. $\delta L = (-V'(x))\dot x + (m\dot x)\ddot x$. And $dL/dt = m\dot x\ddot x - V'(x)\dot x$ (no explicit $t$). Equal, so $F = L$.

(b) $Q = pK - F = m\dot x^2 - (\tfrac12 m\dot x^2 - V) = \tfrac12 m\dot x^2 + V = h$. ($F = L$ depends on $\dot x$, outside the $F(q, t)$ of Theorem 3.4; its proof only used $dQ/dt = \delta L - dF/dt$ along a path, which still holds.)

(c) With $V(x, t)$: $\delta L = -\partial_xV\,\dot x + m\dot x\ddot x$ and $dL/dt = m\dot x\ddot x - \partial_xV\,\dot x - \partial_tV$. So
$$\delta L - \frac{dL}{dt} = \frac{\partial V}{\partial t} = -\frac{\partial L}{\partial t}.$$
$\delta L$ is a total derivative iff $\partial_tV$ is. A function $G$ with $dG/dt = \partial_tV(x, t)$ cannot depend on $\dot x$ (its derivative would contain $\ddot x$), so $G = G(x, t)$, and by the Lemma $\partial_tV$ must be a function of $t$ alone: $V = V_0(x) + f(t)$. Then $\delta L = d(L + f)/dt$, $F = L + f$, and $Q = m\dot x^2 - L - f = \tfrac12 m\dot x^2 + V_0(x)$, conserved. Note $h = \tfrac12 m\dot x^2 + V_0 + f$ itself is *not* conserved ($dh/dt = f'$), but $f(t)$ is physically empty: $-f(t) = -\tfrac{d}{dt}\int f\,dt$ is a total derivative in $L$. For any other $V(x, t)$, $dh/dt = \partial_tV \ne 0$: whoever changes the potential does work on the particle.

## B. Core

**B1.** $\partial L/\partial x = \partial L/\partial y = 0$, $\partial L/\partial z = -mg$, $\partial L/\partial\dot{\mathbf r} = m\dot{\mathbf r}$.

(a) $K = (1, 0, 0)$ and $(0, 1, 0)$ give $\delta L = 0$; charges $m\dot x$, $m\dot y$. $K = (0, 0, 1)$ gives $\delta L = -mg \ne 0$: not strict.

(b) $\delta L = -mg = \tfrac{d}{dt}(-mgt)$, so $F = -mgt$ and $Q = m\dot z - (-mgt) = m\dot z + mgt$. On the solution $\dot z = \dot z_0 - gt$, $Q = m\dot z_0$: the momentum the particle had at $t = 0$. "$Q$ is constant" is $\dot p_z = -mg$ integrated once.

(c) About $\hat{\mathbf z}$: $K = (-y, x, 0)$, $\dot K = (-\dot y, \dot x, 0)$. $\delta L = 0 + 0 + 0 + m\dot x(-\dot y) + m\dot y\dot x + 0 = 0$. Symmetry, $Q = L_z = m(x\dot y - y\dot x)$.

About $\hat{\mathbf x}$: $K = (0, -z, y)$, $\dot K = (0, -\dot z, \dot y)$. $\delta L = (-mg)\,y + m\dot y(-\dot z) + m\dot z\,\dot y = -mgy$. This is velocity-free and depends on $y$, so by the Lemma it is not $dF/dt$ for any $F(q, t)$.

About $\hat{\mathbf y}$: $K = (z, 0, -x)$, $\dot K = (\dot z, 0, -\dot x)$. $\delta L = (-mg)(-x) + m\dot x\dot z - m\dot z\dot x = mgx$. Velocity-free, depends on $x$: not a total derivative. Gravity has a direction; tilting the experiment changes the physics.

(d) $K = (0, 0, t)$, $\dot K = (0, 0, 1)$. $\delta L = (-mg)t + m\dot z = \tfrac{d}{dt}\big(mz - \tfrac12 mgt^2\big)$. So $F = mz - \tfrac12 mgt^2$ and
$$Q = m\dot z\,t - mz + \tfrac12 mgt^2.$$
At $t = 0$, $Q = -mz_0$. Constancy gives $mz = m\dot z\,t + \tfrac12 mgt^2 + mz_0$. Substitute $\dot z = \dot z_0 - gt$ from (b): $z = z_0 + \dot z_0t - gt^2 + \tfrac12 gt^2 = z_0 + \dot z_0t - \tfrac12 gt^2$. The boost charge and the translation charge together *are* the trajectory. The other two boosts are ordinary: $K = (t, 0, 0)$ gives $\delta L = m\dot x = \tfrac{d}{dt}(mx)$ and $Q = m\dot xt - mx$; $K = (0, t, 0)$ gives $Q = m\dot yt - my$.

(e) The ten (§1): time translation, three translations, three rotations, three boosts. Surviving as symmetries or quasi-symmetries: time translation ($h = \tfrac12 m|\dot{\mathbf r}|^2 + mgz$, conserved since $\partial L/\partial t = 0$), all three translations, rotation about $\hat{\mathbf z}$, all three boosts. That is eight. Lost: rotations about $\hat{\mathbf x}$ and $\hat{\mathbf y}$.

They are not really lost. Put $Z = z + \tfrac12 gt^2$; then $\ddot Z = \ddot z + g = 0$, so in the coordinates $(x, y, Z)$ the particle is free and rotation about $\hat{\mathbf x}$ is $(y, Z)\to(y - \varepsilon Z, Z + \varepsilon y)$. In the original coordinates that is $K = (0, -z - \tfrac12 gt^2, y)$, $\dot K = (0, -\dot z - gt, \dot y)$. Then
$$\delta L = (-mg)\,y + m\dot y(-\dot z - gt) + m\dot z\,\dot y = -mg(y + t\dot y) = \frac{d}{dt}(-mgty),$$
so $F = -mgty$ and
$$Q = m\dot y\big(-z - \tfrac12 gt^2\big) + m\dot z\,y + mgty = m(y\dot z - z\dot y) + mg\big(ty - \tfrac12 t^2\dot y\big).$$
With $\dot Z = \dot z + gt$ this is $m(y\dot Z - Z\dot y)$: the $x$-component of angular momentum *measured in the freely falling frame*. About $\hat{\mathbf y}$: $K = (z + \tfrac12 gt^2, 0, -x)$, $\dot K = (\dot z + gt, 0, -\dot x)$, $\delta L = (-mg)(-x) + m\dot x(\dot z + gt) - m\dot z\dot x = mg(x + t\dot x) = \tfrac{d}{dt}(mgtx)$, $F = mgtx$, $Q = m(z\dot x - x\dot z) + mg(\tfrac12 t^2\dot x - tx) = m(Z\dot x - x\dot Z)$. All ten conserved quantities exist; two of them are explicitly time-dependent. This is the equivalence principle in Noether's language: uniform gravity is a free particle seen from an accelerating frame.

**B2.** (a) Invert: $\mathbf r_1 = \mathbf R + (m_2/M)\mathbf s$, $\mathbf r_2 = \mathbf R - (m_1/M)\mathbf s$ (check: $m_1\mathbf r_1 + m_2\mathbf r_2 = M\mathbf R$ and $\mathbf r_1 - \mathbf r_2 = \mathbf s$). Then
$$T = \tfrac12 m_1\big|\dot{\mathbf R} + \tfrac{m_2}{M}\dot{\mathbf s}\big|^2 + \tfrac12 m_2\big|\dot{\mathbf R} - \tfrac{m_1}{M}\dot{\mathbf s}\big|^2 = \tfrac12 M|\dot{\mathbf R}|^2 + \Big(\tfrac{m_1m_2}{M} - \tfrac{m_2m_1}{M}\Big)\dot{\mathbf R}\cdot\dot{\mathbf s} + \tfrac12\frac{m_1m_2^2 + m_2m_1^2}{M^2}|\dot{\mathbf s}|^2.$$
The cross term vanishes and the last coefficient is $m_1m_2(m_1 + m_2)/M^2 = \mu$. With $V = V(s)$: $L = \tfrac12 M|\dot{\mathbf R}|^2 + \tfrac12\mu|\dot{\mathbf s}|^2 - V(s)$.

(b) $X, Y, Z$ (the components of $\mathbf R$) are cyclic. Their momenta $M\dot{\mathbf R} = m_1\dot{\mathbf r}_1 + m_2\dot{\mathbf r}_2 = \mathbf P$: the translation charge of §3.3a. Translating both bodies is translating $\mathbf R$ with $\mathbf s$ fixed, which is why it became a coordinate.

(c) Expand with (a):
$$\mathbf L = m_1\big(\mathbf R + \tfrac{m_2}{M}\mathbf s\big)\times\big(\dot{\mathbf R} + \tfrac{m_2}{M}\dot{\mathbf s}\big) + m_2\big(\mathbf R - \tfrac{m_1}{M}\mathbf s\big)\times\big(\dot{\mathbf R} - \tfrac{m_1}{M}\dot{\mathbf s}\big).$$
The $\mathbf R\times\dot{\mathbf R}$ terms sum to $M\,\mathbf R\times\dot{\mathbf R}$; the mixed terms $\mathbf R\times\dot{\mathbf s}$ and $\mathbf s\times\dot{\mathbf R}$ have coefficients $m_1m_2/M - m_2m_1/M = 0$; the $\mathbf s\times\dot{\mathbf s}$ terms have coefficient $\mu$. So $\mathbf L = M\mathbf R\times\dot{\mathbf R} + \mu\,\mathbf s\times\dot{\mathbf s}$. Each piece is conserved: $\tfrac{d}{dt}(M\mathbf R\times\dot{\mathbf R}) = M\mathbf R\times\ddot{\mathbf R} = \mathbf 0$ because $\ddot{\mathbf R} = \mathbf 0$ by (b); $\tfrac{d}{dt}(\mu\,\mathbf s\times\dot{\mathbf s}) = \mathbf s\times\mu\ddot{\mathbf s} = \mathbf s\times(-V'(s)\,\mathbf s/s) = \mathbf 0$. Responsible symmetries: $K_{\mathbf R} = \hat{\mathbf n}\times\mathbf R$, $K_{\mathbf s} = \mathbf 0$ leaves $\tfrac12 M|\dot{\mathbf R}|^2$ invariant (a dot product, §3.3b step 1) and does not touch the rest, with charge $\hat{\mathbf n}\cdot M\mathbf R\times\dot{\mathbf R}$; $K_{\mathbf R} = \mathbf 0$, $K_{\mathbf s} = \hat{\mathbf n}\times\mathbf s$ leaves $\tfrac12\mu|\dot{\mathbf s}|^2 - V(s)$ invariant, with charge $\hat{\mathbf n}\cdot\mu\,\mathbf s\times\dot{\mathbf s}$. Rotating the whole system is doing both at once, and its charge is the sum.

(d) $K_{\mathbf R} = \mathbf vt$, $\dot K_{\mathbf R} = \mathbf v$, $K_{\mathbf s} = \mathbf 0$. $L$ does not contain $\mathbf R$, so $\delta L = M\dot{\mathbf R}\cdot\mathbf v = \tfrac{d}{dt}(M\mathbf v\cdot\mathbf R)$, $F = M\mathbf v\cdot\mathbf R$, and $Q = M\dot{\mathbf R}\cdot\mathbf vt - M\mathbf v\cdot\mathbf R = \mathbf v\cdot(\mathbf Pt - M\mathbf R)$. This is (3.7). In these coordinates the boost does not touch $\mathbf s$, so the potential term never appears and the computation is one line.

(e) $\sum\dot q_i\,\partial L/\partial\dot q_i = M|\dot{\mathbf R}|^2 + \mu|\dot{\mathbf s}|^2$, so $h = \tfrac12 M|\dot{\mathbf R}|^2 + \tfrac12\mu|\dot{\mathbf s}|^2 + V(s)$, conserved because $\partial L/\partial t = 0$. Since $\dot{\mathbf R}$ is constant, $\tfrac12 M|\dot{\mathbf R}|^2$ is constant, hence $E_{\text{rel}} = h - \tfrac12 M|\dot{\mathbf R}|^2 = \tfrac12\mu|\dot{\mathbf s}|^2 + V(s)$ is conserved by itself. $L$ splits as $L_{\text{cm}}(\dot{\mathbf R}) + L_{\text{rel}}(\mathbf s, \dot{\mathbf s})$ with no coupling, so $\mathbf s$ obeys $\mu\ddot{\mathbf s} = -\nabla V(s)$: one particle of mass $\mu$ in the central potential $V$, with conserved energy $E_{\text{rel}}$ and angular momentum $\mu\,\mathbf s\times\dot{\mathbf s}$. That is where CM-5 starts.

**B3.** $\partial L/\partial x = q_eB\dot y$, $\partial L/\partial y = 0$.

(a) $p_x = m\dot x$, $p_y = m\dot y + q_eBx$ (CM-3 (4.18) with $\mathbf A = (0, Bx, 0)$). Euler–Lagrange: $m\ddot x = q_eB\dot y$ and $\tfrac{d}{dt}(m\dot y + q_eBx) = 0$, i.e. $\ddot x = \omega_c\dot y$, $\ddot y = -\omega_c\dot x$. Let $w = \dot x + i\dot y$. Then $\dot w = \omega_c\dot y - i\omega_c\dot x = -i\omega_c w$, so $w = w_0e^{-i\omega_ct}$: the speed $v = |w|$ is constant and the velocity turns clockwise at angular frequency $\omega_c$. Integrating, $x + iy = z_c + \dfrac{i}{\omega_c}w(t)$ with $z_c$ a constant, a circle about $z_c$ of radius $|w|/\omega_c = mv/(q_eB)$.

(b) $p_y = m\dot y + q_eBx$ is conserved.

(c) $K = (1, 0)$: $\delta L = \partial L/\partial x = q_eB\dot y = \tfrac{d}{dt}(q_eBy)$. So $F = q_eBy$ and $Q_x = p_x - F = m\dot x - q_eBy$. In this gauge $p_x = m\dot x$, and $m\dot x$ is not conserved ($m\ddot x = q_eB\dot y$). The reason for the $F$: shifting $x$ shifts $\mathbf A = (0, Bx, 0)$ by $(0, B\varepsilon, 0) = \nabla(B\varepsilon y)$, a gauge transformation, which changes $L$ by $q_e\dot{\mathbf r}\cdot\nabla(B\varepsilon y) = \varepsilon\tfrac{d}{dt}(q_eBy)$. Space is homogeneous; the Landau gauge is not.

(d) From (a), $z_c = (x + iy) - \tfrac{i}{\omega_c}(\dot x + i\dot y)$, so $X_c = x + \dot y/\omega_c$ and $Y_c = y - \dot x/\omega_c$. Compare: $p_y/(q_eB) = x + m\dot y/(q_eB) = x + \dot y/\omega_c = X_c$ ✓ and $-Q_x/(q_eB) = y - m\dot x/(q_eB) = Y_c$ ✓. Distance: $(x - X_c)^2 + (y - Y_c)^2 = (\dot y^2 + \dot x^2)/\omega_c^2 = (mv/q_eB)^2 = \rho_L^2$.

(e) $K = (-y, x)$, $\dot K = (-\dot y, \dot x)$:
$$\delta L = q_eB\dot y(-y) + 0\cdot x + m\dot x(-\dot y) + (m\dot y + q_eBx)\dot x = q_eB(x\dot x - y\dot y) = \frac{d}{dt}\Big[\tfrac12 q_eB(x^2 - y^2)\Big].$$
$F = \tfrac12 q_eB(x^2 - y^2)$, and $Q_{\text{rot}} = -y\,m\dot x + x(m\dot y + q_eBx) - F = m(x\dot y - y\dot x) + \tfrac12 q_eB(x^2 + y^2)$.

In terms of the circle: let $\boldsymbol\rho = \mathbf r - \mathbf r_c$. By (d), $\rho_x = -\dot y/\omega_c$, $\rho_y = \dot x/\omega_c$, so $m\dot y = -q_eB\rho_x$ and $m\dot x = q_eB\rho_y$. Then $m(x\dot y - y\dot x) = -q_eB(x\rho_x + y\rho_y) = -q_eB\,\mathbf r\cdot\boldsymbol\rho$ and
$$Q_{\text{rot}} = \tfrac12 q_eB\big(r^2 - 2\mathbf r\cdot\boldsymbol\rho\big) = \tfrac12 q_eB\big(|\mathbf r - \boldsymbol\rho|^2 - \rho^2\big) = \tfrac12 q_eB\big(R_c^2 - \rho_L^2\big).$$
Not independent: $R_c^2 = (p_y^2 + Q_x^2)/(q_eB)^2$ and $\rho_L^2 = 2mh/(q_eB)^2$ with $h = \tfrac12 mv^2$, so $Q_{\text{rot}} = (p_y^2 + Q_x^2)/(2q_eB) - mh/(q_eB)$, a function of $p_y$, $Q_x$, $h$. (Phase space is four-dimensional; three independent constants is the most a non-stationary motion can have.)

(f) $p'_x = m\dot x - \tfrac12 q_eBy$, $p'_y = m\dot y + \tfrac12 q_eBx$. Rotation, $K = (-y, x)$: $\delta(\dot x^2 + \dot y^2) = 2\dot x(-\dot y) + 2\dot y\dot x = 0$ and $\delta(x\dot y - y\dot x) = (-y)\dot y + x\dot x - x\dot x - y(-\dot y) = 0$, so $\delta L' = 0$, strict; $Q = -y\,p'_x + x\,p'_y = m(x\dot y - y\dot x) + \tfrac12 q_eB(x^2 + y^2) = Q_{\text{rot}}$. $x$-translation: $\delta L' = \partial L'/\partial x = \tfrac12 q_eB\dot y$, $F = \tfrac12 q_eBy$, $Q = p'_x - F = m\dot x - q_eBy = Q_x$. $y$-translation: $\delta L' = \partial L'/\partial y = -\tfrac12 q_eB\dot x$, $F = -\tfrac12 q_eBx$, $Q = p'_y - F = m\dot y + q_eBx = p_y$. So in this gauge the rotation is strict and both translations are quasi; in the Landau gauge $y$-translation was strict and the other two quasi. The charges are identical. The two Lagrangians differ by a total derivative, $L - L' = \tfrac12 q_eB(x\dot y + y\dot x) = \tfrac{d}{dt}(\tfrac12 q_eBxy)$, so "strict versus quasi" is a property of the representative you wrote down, not of the physics. That is why the hypothesis of Noether's theorem must be invariance *up to a total derivative*.

**B4.** (a) $K = (-y, x)$, $\dot K = (-\dot y, \dot x)$: $\delta L = (-m\omega_1^2x)(-y) + (-m\omega_2^2y)(x) + m\dot x(-\dot y) + m\dot y\dot x = m(\omega_1^2 - \omega_2^2)xy$.

(b) $\dot L_z = m(\dot x\dot y + x\ddot y - \dot y\dot x - y\ddot x) = m(x\ddot y - y\ddot x)$. With $\ddot x = -\omega_1^2x$, $\ddot y = -\omega_2^2y$: $\dot L_z = m(\omega_1^2 - \omega_2^2)xy = \delta L$. This is steps 1–2 of the proof of Theorem 3.3: on shell, $dQ/dt = \delta L$ for any $K$. Step 3 (the hypothesis $\delta L = 0$) fails when $\omega_1 \ne \omega_2$, and $L_z$ changes at exactly the rate $\delta L$.

(c) Write $y = B(\cos\omega t\cos\phi + \sin\omega t\sin\phi)$. With $\cos\omega t = x/A$: $\sin\omega t\sin\phi = y/B - (x/A)\cos\phi$. Square, and add $\cos^2\omega t\sin^2\phi = (x/A)^2\sin^2\phi$:
$$\sin^2\phi = \frac{x^2}{A^2} + \frac{y^2}{B^2} - \frac{2xy\cos\phi}{AB}.$$
A quadratic form in $(x, y)$ with no linear terms, with determinant $(1 - \cos^2\phi)/(A^2B^2) > 0$ when $\sin\phi \ne 0$: an ellipse centered at the origin (a segment when $\sin\phi = 0$). The force center is the ellipse's *center*, not a focus as for Kepler. Angular momentum: $\dot x = -A\omega\sin\omega t$, $\dot y = -B\omega\sin(\omega t - \phi)$, so
$$x\dot y - y\dot x = AB\omega\big[\sin\omega t\cos(\omega t - \phi) - \cos\omega t\sin(\omega t - \phi)\big] = AB\omega\sin\phi,$$
and $L_z = m\omega AB\sin\phi$.

(d) $x = A\cos(\omega_1t - \phi_1)$, $y = B\cos(2\omega_1t - \phi_2)$: both repeat after $2\pi/\omega_1$, so the orbit closes (a Lissajous figure). $L$ is a sum of two uncoupled time-independent oscillators, so $E_x$ and $E_y$ are each conserved ($\dot E_x = m\dot x(\ddot x + \omega_1^2x) = 0$, and the same computation with $y$, $\omega_2$). $L_z$ is not: $\delta L = -3m\omega_1^2xy \ne 0$.

For $Q_3$: $\dot a_x = \ddot x + i\omega_1\dot x = -\omega_1^2x + i\omega_1\dot x = i\omega_1a_x$, and $\dot a_y = -\omega_2^2y + i\omega_2\dot y = i\omega_2a_y = 2i\omega_1a_y$. So
$$\frac{d}{dt}\big(a_x^2\bar a_y\big) = 2a_x\dot a_x\bar a_y + a_x^2\,\overline{\dot a_y} = 2i\omega_1a_x^2\bar a_y - 2i\omega_1a_x^2\bar a_y = 0.$$
Its real part: $a_x^2 = \dot x^2 - \omega_1^2x^2 + 2i\omega_1x\dot x$ and $\bar a_y = \dot y - 2i\omega_1y$, so $\operatorname{Re}(a_x^2\bar a_y) = (\dot x^2 - \omega_1^2x^2)\dot y + 4\omega_1^2x\dot x\,y = Q_3$. Conserved. (The imaginary part $2\omega_1[x\dot x\dot y - y(\dot x^2 - \omega_1^2x^2)]$ is a fourth constant, not independent: four dimensions allow at most three.)

Theorem 3.5 is about charges of the form $\sum p_iK_i(q, t) - F(q, t)$, linear in the momenta. $Q_3$ is cubic in the velocities, outside its scope: it neither promises a point symmetry behind $Q_3$ (there is none) nor rules out such conserved quantities. Like the Runge–Lenz vector, $Q_3$ comes from a transformation that moves velocities, and it exists only because $\omega_2/\omega_1$ is rational.

**B5.** (a) Position $(r\cos\omega t, r\sin\omega t, 0)$; velocity $\dot r(\cos\omega t, \sin\omega t, 0) + r\omega(-\sin\omega t, \cos\omega t, 0)$, so $T = \tfrac12 m(\dot r^2 + \omega^2r^2)$. The rod is horizontal, so $V = 0$ (the rod's vertical normal force cancels gravity). $L = \tfrac12 m(\dot r^2 + \omega^2r^2)$, $m\ddot r = m\omega^2r$, and $r(t) = C_1\cosh\omega t + C_2\sinh\omega t$.

(b) $h = \dot r\,m\dot r - L = \tfrac12 m(\dot r^2 - \omega^2r^2)$, conserved ($\partial L/\partial t = 0$). $E = T = \tfrac12 m(\dot r^2 + \omega^2r^2)$. $T_2 = \tfrac12 m\dot r^2$, $T_1 = 0$, $T_0 = \tfrac12 m\omega^2r^2$. Check (3.11): $T_2 - T_0 + V = h$ ✓, $E - h = m\omega^2r^2 = T_1 + 2T_0$ ✓.

(c) $r(0) = r_0$ gives $C_1 = r_0$; $\dot r(0) = C_2\omega = 0$ gives $C_2 = 0$. So $r = r_0\cosh\omega t$, $\dot r = r_0\omega\sinh\omega t$, and
$$h = \tfrac12 m r_0^2\omega^2(\sinh^2\omega t - \cosh^2\omega t) = -\tfrac12 m\omega^2r_0^2, \qquad E(t) = \tfrac12 m r_0^2\omega^2(\sinh^2\omega t + \cosh^2\omega t) = \tfrac12 m\omega^2r_0^2\cosh 2\omega t.$$
$E - h = \tfrac12 m\omega^2r_0^2(\cosh 2\omega t + 1) = m\omega^2r_0^2\cosh^2\omega t = m\omega^2r^2 = 2T_0$ ✓.

(d) $\tau = dL_z/dt = 2mr\dot r\omega$, power $\tau\omega = 2m\omega^2r\dot r$. Directly: $dE/dt = m\dot r\ddot r + m\omega^2r\dot r = m\dot r\omega^2r + m\omega^2r\dot r = 2m\omega^2r\dot r$ ✓. On the solution of (c): $2m\omega^2r_0^2\omega\cosh\omega t\sinh\omega t = m\omega^3r_0^2\sinh 2\omega t$, which is also $\tfrac{d}{dt}$ of $E(t)$ in (c).

(e) $\omega t = 2$; $\cosh 2 = 3.76220$, $\sinh 2 = 3.62686$, $\cosh 4 = 27.3082$, $\sinh 4 = 27.2899$.
- $r = 0.1\times 3.76220 = 0.37622$ m; $\dot r = 0.1\times 2\times 3.62686 = 0.72537$ m/s.
- $h = -\tfrac12(0.5)(4)(0.01) = -0.0100$ J; $E(0) = +0.0100$ J.
- $E(1\text{ s}) = 0.0100\times\cosh 4 = 0.27308$ J.
- Power $= m\omega^3r_0^2\sinh 4 = 0.5\times 8\times 0.01\times 27.2899 = 1.0916$ W.

(f) In the rotating frame, $h = \tfrac12 m\dot q^2 + V_{\text{eff}}$. For the hoop, $V_{\text{eff}} = -mgR\cos\theta - \tfrac12 mR^2\omega^2\sin^2\theta$ is bounded and, for $\omega^2 < g/R$, has a minimum at $\theta = 0$; the bead is trapped between turning points, and $E - h = mR^2\omega^2\sin^2\theta$ stays below $mR^2\omega^2$. For the rod, $V_{\text{eff}} = -\tfrac12 m\omega^2r^2$ is unbounded below, $r = 0$ is a maximum, and the bead runs away exponentially; $E - h = m\omega^2r^2$ grows without bound, all of it supplied by the motor. Gravity's restoring term is what the rod lacks.

**B6.** (a) Energy at the turning point: $\tfrac12\dot x^2 + V(x) = V(A)$, so $|\dot x| = \sqrt{2(V(A) - V(x))}$. The time from $0$ to $A$ is $\int_0^A dx/|\dot x|$, and since $V$ is even the period is four of these. Put $x = Au$: $V(A) - V(Au) = A^k(1 - u^k)/k$ and $dx = A\,du$, so
$$T(A) = 4\int_0^1\frac{A\,du}{\sqrt{2A^k(1 - u^k)/k}} = A^{1 - k/2}\cdot 4\int_0^1\frac{du}{\sqrt{2(1 - u^k)/k}} = A^{1 - k/2}\,T(1).$$

(b) $k = 4$: $T(1) = 4\sqrt2\int_0^1(1 - u^4)^{-1/2}du$. With $s = u^4$, $du = \tfrac14 s^{-3/4}ds$:
$$\int_0^1(1 - u^4)^{-1/2}du = \tfrac14\int_0^1 s^{-3/4}(1 - s)^{-1/2}ds = \tfrac14 B\big(\tfrac14, \tfrac12\big).$$
So $T(1) = \sqrt2\,B(\tfrac14, \tfrac12) = \sqrt2\,\Gamma(\tfrac14)\sqrt\pi/\Gamma(\tfrac34) = \sqrt{2\pi}\,\Gamma(\tfrac14)/\Gamma(\tfrac34)$. Numbers: $2.506628\times 3.625610/1.225417 = 7.416299$, matching the project. $T(2) = T(1)/2 = 3.708149$.

(c) $k = 4$: $2\langle T\rangle = 4\langle V\rangle$, so $\langle T\rangle = 2\langle V\rangle$ and $E = \langle T\rangle + \langle V\rangle = 3\langle V\rangle$: $\langle V\rangle = E/3$, $\langle T\rangle = 2E/3$. At $A = 1$, $E = V(1) = \tfrac14$: $\langle T\rangle = 1/6 = 0.166667$, $\langle V\rangle = 1/12 = 0.083333$. The project gets both.

(d) $k = -1$: $2\langle T\rangle = -\langle V\rangle$, so $E = \langle T\rangle + \langle V\rangle = \tfrac12\langle V\rangle$ and $\langle V\rangle = 2E = -\kappa/a$. Since $\langle V\rangle = -\kappa\langle 1/r\rangle$, $\langle 1/r\rangle = 1/a$ exactly, for every eccentricity (project: $1.000000$ at $a = 1$, $e = 0.5$). The time average of $r$ itself is *not* $a$: the project gets $1.125000 = a(1 + e^2/2)$. (Derivation in CM-5: with the eccentric anomaly $\psi$, $r = a(1 - e\cos\psi)$ and $dt \propto (1 - e\cos\psi)\,d\psi$, so $\langle r\rangle = \tfrac{a}{2\pi}\int_0^{2\pi}(1 - e\cos\psi)^2d\psi = a(1 + e^2/2)$.) The planet lingers near apocenter, which pulls $\langle r\rangle$ above $a$; the virial theorem only controls $\langle 1/r\rangle$.

**B7.** (a) $\mathbf p\cdot\mathbf L = \mathbf p\cdot(\mathbf r\times\mathbf p) = 0$, so by Lagrange's identity $|\mathbf p\times\mathbf L|^2 = p^2L^2 - (\mathbf p\cdot\mathbf L)^2 = p^2L^2$. By cyclic permutation of the triple product, $\mathbf r\cdot(\mathbf p\times\mathbf L) = \mathbf L\cdot(\mathbf r\times\mathbf p) = L^2$, so $\hat{\mathbf r}\cdot(\mathbf p\times\mathbf L) = L^2/r$. Now square (3.14):
$$A^2 = |\mathbf p\times\mathbf L|^2 - 2m\kappa\,\hat{\mathbf r}\cdot(\mathbf p\times\mathbf L) + m^2\kappa^2 = p^2L^2 - \frac{2m\kappa L^2}{r} + m^2\kappa^2 = 2mL^2\Big(\frac{p^2}{2m} - \frac{\kappa}{r}\Big) + m^2\kappa^2 = m^2\kappa^2 + 2mEL^2.$$

(b) Divide by $m^2\kappa^2$: $e = \sqrt{1 + 2EL^2/(m\kappa^2)}$. $E < 0$: $e < 1$, ellipse ($e = 0$, a circle, at the lowest energy allowed for given $L$, $E = -m\kappa^2/2L^2$). $E = 0$: $e = 1$, parabola. $E > 0$: $e > 1$, hyperbola. Project orbit: $e = \sqrt{1 + 2(-\tfrac12)(\tfrac34)} = \sqrt{1/4} = 0.5000$.

(c) At pericenter $\mathbf r\cdot\mathbf p = 0$, so BAC–CAB gives $\mathbf p\times\mathbf L = \mathbf p\times(\mathbf r\times\mathbf p) = \mathbf r\,p^2 - \mathbf p(\mathbf p\cdot\mathbf r) = p^2\mathbf r$, and with $p = L/r_p$ (momentum perpendicular to $\mathbf r$):
$$\mathbf A = \Big(\frac{L^2}{r_p} - m\kappa\Big)\hat{\mathbf r}_p.$$
Sign: the radial equation from the polar Lagrangian is $m\ddot r = L^2/(mr^3) - \kappa/r^2$. Pericenter is a minimum of $r$, so $\ddot r \ge 0$ there, i.e. $L^2/(mr_p^3) \ge \kappa/r_p^2$, i.e. $L^2/r_p \ge m\kappa$. The coefficient is $\ge 0$: $\mathbf A$ points from the force center to the pericenter (and vanishes for a circle). At apocenter the same algebra gives $\mathbf A = (L^2/r_a - m\kappa)\hat{\mathbf r}_a$, but now $r$ is a maximum, $\ddot r \le 0$, and the coefficient is $\le 0$: $\mathbf A$ points *opposite* to $\hat{\mathbf r}_a$, which is again toward the pericenter. One conserved vector, two consistent evaluations.

(d) Phase space of one particle is six-dimensional. With only $E$ and $\mathbf L$ (four independent numbers), an orbit lies on a $6 - 4 = 2$-dimensional surface: $\mathbf L$ fixes the plane, $E$ and $|\mathbf L|$ fix the radial turning points, and a generic orbit fills the annulus $r_{\min}\le r\le r_{\max}$ densely as a rosette. With $\mathbf A$ added, seven numbers minus two relations ($\mathbf A\cdot\mathbf L = 0$, part (a)) leave five, so the orbit is a curve: it closes. The project's $V = -\kappa/r - \epsilon/r^3$ is central but not $1/r$, so only $E$ and $\mathbf L$ survive: the rosette of `lrl_precession.png`.

**B8.** A compact version of the reference implementation (the script's final section):

```python
eps_list = np.array([0.0005, 0.001, 0.002, 0.004])
adv = []
for ep in eps_list:
    sol_e = integrate(lambda t_, s, ep=ep: kepler_rhs(t_, s, eps=ep), s0, 8 * T0,
                      n_out=800, events=pericenter)
    ps = sol_e.y_events[0][:, :3]                       # positions at pericenter passages
    ang = np.unwrap(np.arctan2(ps @ e2, ps @ e1))      # in-plane angle of each pericenter
    adv.append(np.mean(np.diff(ang)))                   # advance per orbit [rad]
adv_pred = 6 * np.pi * eps_list / Lmag**4               # first order, m = kappa = 1
rel_err = (np.array(adv) - adv_pred) / adv_pred
c_fit = np.sum(eps_list * rel_err) / np.sum(eps_list**2)   # least squares, line through 0
```

The line through the origin minimizes $\sum_i(r_i - c\,\epsilon_i)^2$; setting the derivative in $c$ to zero gives $c = \sum\epsilon_ir_i/\sum\epsilon_i^2$. Results ($L^4 = 0.5625$):

| $\epsilon$ | 0.0005 | 0.001 | 0.002 | 0.004 |
|---|---|---|---|---|
| measured [°] | 0.967 | 1.947 | 3.950 | 8.135 |
| first order [°] | 0.960 | 1.920 | 3.840 | 7.680 |
| relative error | 0.0070 | 0.0141 | 0.0287 | 0.0593 |

Slope $c \approx 14.7$: relative error $\approx 14.7\,\epsilon$. Linear extrapolation puts the 1% error at $\epsilon \approx 0.00068$ and 10% at $\epsilon \approx 0.0068$.

Why linear: the exact advance is a smooth function of $\epsilon$ vanishing at $\epsilon = 0$, so $\Delta\phi(\epsilon) = c_1\epsilon + c_2\epsilon^2 + O(\epsilon^3)$. First-order perturbation theory computes $c_1\epsilon$ exactly. The error is $c_2\epsilon^2 + \dots$, and relative to $c_1\epsilon$ it is $(c_2/c_1)\epsilon + O(\epsilon^2)$. The fit measures $c_2/c_1 \approx 14.7$; positive means the second-order term adds to the advance. Halving $\epsilon$ halves the relative error, which is what the table shows.

## C. Challenge

**C1.** Equations of motion: $\ddot x = -\omega^2x$, $\ddot y = -\omega^2y$.

(a) $\dot Q_1 = m(\dot x\ddot x - \dot y\ddot y) + m\omega^2(x\dot x - y\dot y) = m(-\omega^2x\dot x + \omega^2y\dot y) + m\omega^2(x\dot x - y\dot y) = 0$.
$\dot Q_2 = m(\ddot x\dot y + \dot x\ddot y + \omega^2\dot xy + \omega^2x\dot y) = m(-\omega^2x\dot y - \omega^2\dot xy + \omega^2\dot xy + \omega^2x\dot y) = 0$.
$\dot L_z = m(x\ddot y - y\ddot x) = m(-\omega^2xy + \omega^2yx) = 0$.

(b) $K = (\dot y, \dot x)$, $\dot K = (\ddot y, \ddot x)$:
$$\delta L = -m\omega^2x\dot y - m\omega^2y\dot x + m\dot x\ddot y + m\dot y\ddot x = \frac{d}{dt}\big[m\dot x\dot y - m\omega^2xy\big].$$
$F = m(\dot x\dot y - \omega^2xy)$, and $Q = m\dot x\dot y + m\dot y\dot x - F = m(\dot x\dot y + \omega^2xy) = Q_2$.

$K = (\dot x, -\dot y)$, $\dot K = (\ddot x, -\ddot y)$: $\delta L = -m\omega^2x\dot x + m\omega^2y\dot y + m\dot x\ddot x - m\dot y\ddot y = \tfrac{d}{dt}(L_x - L_y)$ with $L_x = \tfrac12 m\dot x^2 - \tfrac12 m\omega^2x^2$, $L_y = \tfrac12 m\dot y^2 - \tfrac12 m\omega^2y^2$. So $F = L_x - L_y$ and $Q = m\dot x^2 - m\dot y^2 - L_x + L_y = \tfrac12 m(\dot x^2 - \dot y^2) + \tfrac12 m\omega^2(x^2 - y^2) = Q_1$.

`total_derivative_potential` builds $F(q, t)$ by the Poincaré-lemma integral, which needs $\delta L$ to be linear in velocities with no accelerations. Velocity-dependent $K$ puts accelerations into $\dot K$, hence into $\delta L$, and $F$ depends on velocities; the function refuses, correctly, because it only searches the $F(q, t)$ class. The proof of Theorem 3.4 does not care: on shell $\tfrac{d}{dt}\sum p_iK_i = \sum(\dot p_iK_i + p_i\dot K_i) = \delta L$ for any $K$, and subtracting $dF/dt$ leaves zero whatever $F$ depends on.

(c) Let $\mathbf u = (\dot x, \omega x)$ and $\mathbf w = (\dot y, \omega y)$ in $\mathbb R^2$. Then $E = \tfrac12 m(|\mathbf u|^2 + |\mathbf w|^2)$, $Q_1 = \tfrac12 m(|\mathbf u|^2 - |\mathbf w|^2)$, $Q_2 = m\,\mathbf u\cdot\mathbf w$, and $u_1w_2 - u_2w_1 = \omega(\dot xy - x\dot y)$, so $\omega L_z = -m(u_1w_2 - u_2w_1)$. Hence
$$E^2 - Q_1^2 = m^2|\mathbf u|^2|\mathbf w|^2, \qquad Q_2^2 + \omega^2L_z^2 = m^2\big[(\mathbf u\cdot\mathbf w)^2 + (u_1w_2 - u_2w_1)^2\big] = m^2|\mathbf u|^2|\mathbf w|^2,$$
the last step being Lagrange's identity in $\mathbb R^2$ (expand both squares; the cross terms $\pm 2u_1w_1u_2w_2$ cancel). So $E^2 = Q_1^2 + Q_2^2 + \omega^2L_z^2$. ∎

(d) Four constants, one relation: three independent, in a four-dimensional phase space. Each orbit lies on a one-dimensional curve, so every bounded orbit closes (they are the centered ellipses of B4(c), all with period $2\pi/\omega$). Kepler is the same story one dimension up: five independent constants in six dimensions (B7(d)). These are the only two central potentials in which every bounded orbit closes (Bertrand's theorem, CM-5).

(e) $x = \cos t$, $y = 2\sin t$ solves $\ddot x = -x$, $\ddot y = -y$; $\dot x = -\sin t$, $\dot y = 2\cos t$.
- $E = \tfrac12(\sin^2t + 4\cos^2t) + \tfrac12(\cos^2t + 4\sin^2t) = \tfrac52 = 2.50$.
- $Q_1 = \tfrac12(\sin^2t - 4\cos^2t) + \tfrac12(\cos^2t - 4\sin^2t) = -\tfrac32 = -1.50$.
- $Q_2 = (-\sin t)(2\cos t) + (\cos t)(2\sin t) = 0$.
- $L_z = \cos t\cdot 2\cos t - 2\sin t\cdot(-\sin t) = 2$.
- $E^2 = 6.25 = 2.25 + 0 + 4$ ✓.

(f) $\dot a_x = \ddot x + i\omega\dot x = -\omega^2x + i\omega\dot x = i\omega(\dot x + i\omega x) = i\omega a_x$, and in the same way $\dot a_y = -\omega^2y + i\omega\dot y = i\omega a_y$. $|a_x|^2 = \dot x^2 + \omega^2x^2$, so $\tfrac12 m(|a_x|^2 + |a_y|^2) = E$ and $\tfrac12 m(|a_x|^2 - |a_y|^2) = Q_1$. And
$$m\,a_x\bar a_y = m(\dot x + i\omega x)(\dot y - i\omega y) = m\big[\dot x\dot y + \omega^2xy + i\omega(x\dot y - \dot xy)\big] = Q_2 + i\omega L_z.$$
The linear maps of $(a_x, a_y)\in\mathbb C^2$ that preserve $|a_x|^2 + |a_y|^2$ are the unitary matrices, $U(2)$. Every solution is $\mathbf a(t) = e^{i\omega t}\mathbf a(0)$, and multiplication by $e^{i\omega t}$ commutes with every matrix, so a unitary map of the initial data maps solutions to solutions of the same energy. The overall phase $e^{i\alpha}$ is time evolution itself (generated by $E$); what is left is $SU(2)$, three-dimensional, with the three charges $Q_1, Q_2, \omega L_z$. That is the hidden symmetry.

A rotation of the plane by $\alpha$ acts on positions and velocities alike, so $\mathbf u' = \cos\alpha\,\mathbf u - \sin\alpha\,\mathbf w$, $\mathbf w' = \sin\alpha\,\mathbf u + \cos\alpha\,\mathbf w$ (a real, hence unitary, matrix acting on $(a_x, a_y)$). Then $|\mathbf u'|^2 - |\mathbf w'|^2 = \cos2\alpha(|\mathbf u|^2 - |\mathbf w|^2) - 2\sin2\alpha\,\mathbf u\cdot\mathbf w$ and $\mathbf u'\cdot\mathbf w' = \tfrac12\sin2\alpha(|\mathbf u|^2 - |\mathbf w|^2) + \cos2\alpha\,\mathbf u\cdot\mathbf w$, and $u'_1w'_2 - u'_2w'_1 = u_1w_2 - u_2w_1$ (a determinant times $\det = 1$). So
$$Q_1' = \cos2\alpha\,Q_1 - \sin2\alpha\,Q_2, \qquad Q_2' = \sin2\alpha\,Q_1 + \cos2\alpha\,Q_2, \qquad L_z' = L_z.$$
The vector $(Q_1, Q_2, \omega L_z)$, of length $E$ by (c), turns by $2\alpha$ about its third axis. The angle doubles because $(Q_1, Q_2)$ encodes the orientation of the ellipse's axes, which returns to itself after a rotation by $\pi$.

**C2.** (a) $\int_0^R 4\pi r^3\,\dfrac{dP}{dr}dr = \big[4\pi r^3P\big]_0^R - \int_0^R 12\pi r^2P\,dr = 0 - 3\int P\,dV$, using $P(R) = 0$ and $r^3P \to 0$ at the center. The right side: $-\int_0^R 4\pi r^3\,\dfrac{Gm(r)\rho}{r^2}dr = -\int_0^R\dfrac{Gm(r)}{r}\,4\pi r^2\rho\,dr = -\int\dfrac{Gm\,dm}{r} = V_{\text{grav}}$. So $-3\int P\,dV = V_{\text{grav}}$. ∎

(b) $3\int P\,dV = 2\int u\,dV = 2U_{\text{th}}$, so $2U_{\text{th}} = -V_{\text{grav}}$, at every instant, for any static star. In the language of Theorem 3.9: $G = \sum_a\mathbf p_a\cdot\mathbf r_a = \tfrac12\,dI/dt$ with $I = \sum_a m_a r_a^2$, and a static star has constant (coarse-grained) $I$, so $dG/dt = \tfrac12\ddot I = 0$ and step 2 of that proof needs no time average.

(c) $m(r) = \tfrac43\pi r^3\rho$, so the equilibrium condition is $dP/dr = -\tfrac43\pi G\rho^2r$. The candidate $P = \tfrac{2\pi}{3}G\rho^2(R^2 - r^2)$ has $dP/dr = -\tfrac{4\pi}{3}G\rho^2r$ ✓ and $P(R) = 0$ ✓. Then
$$3\int P\,dV = 3\int_0^R\tfrac{2\pi}{3}G\rho^2(R^2 - r^2)\,4\pi r^2dr = 8\pi^2G\rho^2\Big(\frac{R^5}{3} - \frac{R^5}{5}\Big) = \frac{16\pi^2}{15}G\rho^2R^5,$$
and with $M = \tfrac43\pi R^3\rho$, $\tfrac35 GM^2/R = \tfrac35 G\cdot\tfrac{16\pi^2}{9}R^5\rho^2 = \tfrac{16\pi^2}{15}G\rho^2R^5$ ✓, equal to $-V_{\text{grav}}$ from Example 4.5(b). The project checks it numerically ($G = R = 1$, $\rho = 1.3$): both sides $17.79161$.

(d) $E = U_{\text{th}} + V_{\text{grav}} = U_{\text{th}} - 2U_{\text{th}} = -U_{\text{th}} = -\tfrac32 Nk_BT$, so $dE/dT = -\tfrac32 Nk_B < 0$.

*Star.* Radiating lowers $E$, so $T$ rises. With $k_BT = GM\bar m/(5R)$, $R\propto 1/T$: the star shrinks. Of the gravitational energy released, half is radiated and half heats the gas ($dV_{\text{grav}} = 2\,dE$, $dU_{\text{th}} = -dE$).

*Satellite.* Circular orbit: $v^2 = GM_E/r$ and, per kilogram, $E = \tfrac12 v^2 - GM_E/r = -GM_E/(2r) = -\tfrac12 v^2$, the virial relation $E = -T$. At $r = 6371 + 400 = 6771$ km: $v = \sqrt{3.986004\times10^{14}/6.771\times10^6} = 7672.6$ m/s. $dE/dr = GM_E/(2r^2)$, so descending 1 km costs $3.986004\times10^{14}\times10^3/(2\times(6.771\times10^6)^2) = 4347$ J/kg. $v\propto r^{-1/2}$ gives $dv = -\tfrac{v}{2r}dr$: descending 1 km gains $7672.6\times10^3/(2\times6.771\times10^6) = 0.567$ m/s. Drag removes energy and the satellite speeds up; gravity does twice the work drag does.

*Sun, redone.* $|V_{\text{grav}}| \approx 1.8\,GM^2/R$ is three times $0.6\,GM^2/R$, so the estimate triples: $T \approx 3\times 2.77\times10^6 = 8.32\times10^6$ K, about half the central $1.57\times10^7$ K.

(e) In the canonical ensemble (system in contact with a heat bath, SM-4) the heat capacity is $C = (\langle E^2\rangle - \langle E\rangle^2)/(k_BT^2) \ge 0$, so negative $C$ is forbidden there. In the microcanonical ensemble (isolated system at fixed $E$, $1/T = dS/dE$) nothing forbids it: $C < 0$ means $S(E)$ is locally convex. For ordinary systems with short-range forces such a region is unstable to splitting into phases and the ensembles agree; gravity is long-range, the energy is not additive over subsystems, and the ensembles genuinely differ. Put a self-gravitating gas against a bath and there is no stable equilibrium: if it is slightly hotter than the bath it loses heat, gets hotter, and loses more; if slightly colder it gains heat, cools, and expands. Lynden-Bell and Wood (1968) called the runaway the gravothermal catastrophe.

**C3.** (a) $k = -2$: $\beta = \alpha^{1 - k/2} = \alpha^2$, and angular momentum scales as $\alpha^{1 + k/2} = \alpha^0$. Suppose a circular orbit of radius $r_0$ exists with angular momentum $\ell$. Scaling by any $\alpha$ gives a circular orbit of radius $\alpha r_0$ with the *same* $\ell$. So for a given $\ell$ circular orbits come in a continuum of radii or not at all: no radius is selected, so none can sit at an isolated minimum of $V_{\text{eff}}$, and none is stable.

(b) $V_{\text{eff}} = \dfrac{\ell^2}{2mr^2} - \dfrac{c}{r^2} = \dfrac{\ell^2 - 2mc}{2mr^2}$.
- $\ell^2 < 2mc$: $V_{\text{eff}}$ is negative, increasing, $\to -\infty$ at $r = 0$. No inner turning point: an inward-moving particle (or any with $E < 0$) falls to the center.
- $\ell^2 = 2mc$: $V_{\text{eff}}\equiv 0$, $\ddot r = 0$. Circular at every radius if $\dot r = 0$ (the continuum of (a)); otherwise $r$ changes linearly in time, in to the center or out to infinity.
- $\ell^2 > 2mc$: a repulsive $1/r^2$ barrier, weaker than the free centrifugal one. The particle comes in to $r_{\min}$ and leaves to infinity. No bound orbits.

(c) $F = -dV/dr = -2c/r^3$, so $F(1/u) = -2cu^3$ and Binet gives $u'' + u = (m/\ell^2u^2)(2cu^3) = (2mc/\ell^2)u$:
$$u'' + \Big(1 - \frac{2mc}{\ell^2}\Big)u = 0.$$
Initial data: $\dot r = -\tfrac{\ell}{m}u'$ (from $r = 1/u$ and $\dot\theta = \ell u^2/m$), so starting at $r_0$ with $\dot r = 0$ means $u(0) = 1/r_0$, $u'(0) = 0$.
- $\ell^2 < 2mc$: $u'' = \lambda^2u$, $\lambda^2 = 2mc/\ell^2 - 1$. General solution $u = C_1\cosh\lambda\theta + C_2\sinh\lambda\theta$, which is finite at every finite $\theta$; $r\to 0$ ($u\to\infty$) needs $\theta\to\infty$. The particle reaches the center after infinitely many turns. With our data, $u = u_0\cosh\lambda\theta$.
- $\ell^2 = 2mc$: $u'' = 0$, $u = u_0 + u_0'\theta$: a circle if $u_0' = 0$, otherwise a hyperbolic spiral $r = 1/(u_0 + u_0'\theta)$.
- $\ell^2 > 2mc$: $u'' = -\lambda'^2u$, $\lambda'^2 = 1 - 2mc/\ell^2 \in (0, 1)$. With our data $u = u_0\cos\lambda'\theta$, which reaches $u = 0$ ($r = \infty$) at the finite angle $\theta_\infty = \pi/(2\lambda')$. Since $\lambda' < 1$, $\theta_\infty > \pi/2$: the attraction bends the path more than a free particle's.

(d) $m = c = 1$, $r(0) = 1$, $\dot r(0) = 0$, so $E = V_{\text{eff}}(1) = (\ell^2 - 2)/2$ and $\tfrac12\dot r^2 + (\ell^2 - 2)/(2r^2) = (\ell^2 - 2)/2$:
$$\dot r^2 = (\ell^2 - 2)\,\frac{r^2 - 1}{r^2}.$$
$\ell^2 < 2$: $r\le 1$ and $r\dot r = -\sqrt{2 - \ell^2}\sqrt{1 - r^2}$, i.e. $\tfrac{d}{dt}\sqrt{1 - r^2} = -r\dot r/\sqrt{1 - r^2} = \sqrt{2 - \ell^2}$. So $t = \sqrt{1 - r^2}/\sqrt{2 - \ell^2}$, and $r = 0$ is reached at the finite time $1/\sqrt{2 - \ell^2}$: infinitely many turns in finite time.
$\ell^2 > 2$: $r\ge 1$, $\tfrac{d}{dt}\sqrt{r^2 - 1} = \sqrt{\ell^2 - 2}$, $t = \sqrt{r^2 - 1}/\sqrt{\ell^2 - 2}$.

Numbers. $\ell^2 = 1.9$, $r = 0.05$: $t = \sqrt{0.9975}/\sqrt{0.1} = 3.1583$; $\lambda = \sqrt{2/1.9 - 1} = 0.229416$, and $u = \cosh\lambda\theta = 20$ at $\theta = \operatorname{arccosh}20/\lambda = 3.688254/0.229416 = 16.077$ rad $= 2.56$ turns. $\ell^2 = 2.1$, $r = 20$: $t = \sqrt{399}/\sqrt{0.1} = 63.166$; $\lambda' = \sqrt{1 - 2/2.1} = 0.218218$, $\cos\lambda'\theta = 0.05$ at $\theta = \arccos 0.05/\lambda' = 1.520775/0.218218 = 6.9691$ rad, and $\theta_\infty = \pi/(2\lambda') = 7.1983$ rad. The project's integrations agree with all four numbers to the printed digits (`inverse_square_orbits.png`).

(e) Dimensions: $[m] = \mathrm M$, $[c] = [V][r^2] = \mathrm{ML^4T^{-2}}$, $[\ell] = \mathrm{ML^2T^{-1}}$. For $m^ac^b\ell^d$ to be a length we need $a + b + d = 0$ (M), $4b + 2d = 1$ (L), $-2b - d = 0$ (T). The last gives $d = -2b$, and then the L equation reads $0 = 1$. For a time: $4b + 2d = 0$ and $-2b - d = 1$, which again contradict ($d = -2b$ makes the second read $0 = 1$). The only combination is the dimensionless $\ell^2/(mc)$, and it is exactly the parameter whose value, above or below 2, decides the motion in (b). This is (a) seen through units: with no length in the problem, nothing can pick the radius of a circular orbit or the size of a bound one. The scale is set by the initial conditions alone ($r_0$, or $E$). Compare Kepler: $\ell^2/(m\kappa)$ *is* a length (the semi-latus rectum), which is why Kepler orbits have definite sizes.

**C4.** (a) $p = me^{\gamma t}\dot x$, so $\dot p = me^{\gamma t}(\ddot x + \gamma\dot x)$; $\partial L/\partial x = -ke^{\gamma t}x$. Euler–Lagrange: $e^{\gamma t}(m\ddot x + \gamma m\dot x + kx) = 0$, and $e^{\gamma t}\ne 0$. ∎

(b) $h = \dot xp - L = e^{\gamma t}\big(\tfrac12 m\dot x^2 + \tfrac12 kx^2\big)$, and by (3.9) $dh/dt = -\partial L/\partial t = -\gamma L = -\gamma e^{\gamma t}\big(\tfrac12 m\dot x^2 - \tfrac12 kx^2\big)$. Not conserved; it oscillates in sign. (Direct check: $\tfrac{d}{dt}h = \gamma h + e^{\gamma t}\dot x(m\ddot x + kx) = \gamma h - \gamma me^{\gamma t}\dot x^2 = -\gamma L$.)

(c) $K = \dot x + \tfrac12\gamma x$, $\dot K = \ddot x + \tfrac12\gamma\dot x$.
$$\delta L = -ke^{\gamma t}x\big(\dot x + \tfrac12\gamma x\big) + me^{\gamma t}\dot x\big(\ddot x + \tfrac12\gamma\dot x\big) = e^{\gamma t}(m\dot x\ddot x - kx\dot x) + \gamma e^{\gamma t}\big(\tfrac12 m\dot x^2 - \tfrac12 kx^2\big).$$
And $dL/dt = \gamma L + e^{\gamma t}(m\dot x\ddot x - kx\dot x)$: the same. So $F = L$ and
$$Q = pK - L = me^{\gamma t}\dot x\big(\dot x + \tfrac12\gamma x\big) - e^{\gamma t}\big(\tfrac12 m\dot x^2 - \tfrac12 kx^2\big) = e^{\gamma t}\Big[\tfrac12 m\dot x^2 + \tfrac12 kx^2 + \tfrac12\gamma m\,x\dot x\Big].$$
Check: $\dot Q = \gamma Q + e^{\gamma t}\big[m\dot x\ddot x + kx\dot x + \tfrac12\gamma m(\dot x^2 + x\ddot x)\big]$. Substitute $m\ddot x = -\gamma m\dot x - kx$: $m\dot x\ddot x = -\gamma m\dot x^2 - kx\dot x$ and $\tfrac12\gamma mx\ddot x = -\tfrac12\gamma^2mx\dot x - \tfrac12\gamma kx^2$. The bracket becomes $-\tfrac12\gamma m\dot x^2 - \tfrac12\gamma^2mx\dot x - \tfrac12\gamma kx^2 = -\gamma e^{-\gamma t}Q$, so $\dot Q = 0$. ∎

Why the pieces fail separately: $K = \dot x$ gives $\delta L = dL/dt - \partial L/\partial t = dL/dt - \gamma L$, and $K = \tfrac12\gamma x$ gives $\delta L = \tfrac12\gamma\big(x\,\partial_xL + \dot x\,\partial_{\dot x}L\big) = \gamma L$ ($L$ is homogeneous of degree 2 in $(x, \dot x)$, Euler's theorem). $\gamma L$ is not a total derivative: a $G$ with $dG/dt = \gamma L$ cannot depend on $\dot x$ (no $\ddot x$ in $\gamma L$), and then $dG/dt$ is linear in $\dot x$ while $\gamma L$ has a $\dot x^2$ term. The sum cancels the $\mp\gamma L$. In words: the damped oscillator looks the same a moment later *if you also magnify the amplitude to undo the decay*.

(d) $x = e^{-\gamma t/2}X$, $\dot x = e^{-\gamma t/2}(\dot X - \tfrac12\gamma X)$, so $e^{\gamma t}x^2 = X^2$ and so on:
$$Q = \tfrac12 m\big(\dot X - \tfrac12\gamma X\big)^2 + \tfrac12 kX^2 + \tfrac12\gamma mX\big(\dot X - \tfrac12\gamma X\big) = \tfrac12 m\dot X^2 + \tfrac12\Big(k - \frac{m\gamma^2}{4}\Big)X^2 = \tfrac12 m\dot X^2 + \tfrac12 m\Omega^2X^2.$$
(The $X\dot X$ terms cancel; the $X^2$ terms are $\tfrac18 m\gamma^2 - \tfrac14 m\gamma^2 = -\tfrac18 m\gamma^2$.) And $X$ is an undamped oscillator: $\ddot x = e^{-\gamma t/2}(\ddot X - \gamma\dot X + \tfrac14\gamma^2X)$, and substituting into (a) gives $m\ddot X + (k - \tfrac14 m\gamma^2)X = 0$.

(e) $Q(0) = \tfrac12\cdot 0 + \tfrac12\cdot 1 + 0 = 0.5$; the project finds it constant to $1.5\times10^{-11}$ over $0\le t\le 20$.

Range of $h$. With $\omega_0^2 = k/m$, $h = e^{\gamma t}\tfrac{m}{2}(\dot x^2 + \omega_0^2x^2)$ and $Q = e^{\gamma t}\tfrac{m}{2}(\dot x^2 + \omega_0^2x^2 + \gamma x\dot x)$. Put $\dot x = \rho\cos\varphi$, $\omega_0x = \rho\sin\varphi$: then $\gamma x\dot x = \tfrac{\gamma}{2\omega_0}\rho^2\sin2\varphi$ and
$$\frac{h}{Q} = \frac{1}{1 + \frac{\gamma}{2\omega_0}\sin2\varphi} \in \Big[\frac{\omega_0}{\omega_0 + \gamma/2},\ \frac{\omega_0}{\omega_0 - \gamma/2}\Big].$$
The phase-plane angle $\varphi$ sweeps through all values every cycle ($\Omega\approx 1$, three cycles in 20 s), so both ends are reached: $h\in[0.5/1.1,\ 0.5/0.9] = [0.4545, 0.5556]$, which is what the integration finds. Note what this says: $h = e^{\gamma t}\times(\text{mechanical energy})$ stays within 11% of $Q$, because the mechanical energy decays like $e^{-\gamma t}$ up to a wobble.

$x(20)$: $X(0) = 1$, $\dot X(0) = \dot x(0) + \tfrac12\gamma x(0) = 0.1$, so $X = \cos\Omega t + (0.1/\Omega)\sin\Omega t$ with $\Omega = \sqrt{1 - 0.01} = 0.994987$, $0.1/\Omega = 0.100504$. Then $x(20) = e^{-2}\big[\cos(19.8997) + 0.100504\sin(19.8997)\big] = 0.135335\times[\cdots] = 0.079116$, matching the integration to $10^{-9}$.
