# CM-8 Solutions

Conventions as in the notes. Numbers marked "project" are printed or asserted by Part 13 of `project.py` (B2, B3, B6 in Parts 6 and 12).

## A. Warm-up

**A1.** (a) $f'(v) = v^3$ is increasing, so (3.2) gives $p = v^3$, $v^* = p^{1/3}$ (real cube root). Then
$$g(p) = p\cdot p^{1/3} - \tfrac14p^{4/3} = \tfrac34|p|^{4/3}, \qquad g'(p) = p^{1/3} = v^*. \checkmark$$
The exponents $4$ and $4/3$ satisfy $\frac14 + \frac34 = 1$: Legendre transforms of powers give conjugate (Hölder) exponents.

(b) $f'(v) = mc^2\cdot\frac{v/c^2}{\sqrt{1 - v^2/c^2}} = \gamma mv$ with $\gamma = (1 - v^2/c^2)^{-1/2}$, and $f'' = \gamma^3m > 0$, so $f$ is convex. Invert $p = \gamma mv$: $p^2(1 - v^2/c^2) = m^2v^2$ gives $v^* = \dfrac{pc^2}{\sqrt{p^2c^2 + m^2c^4}}$ and $\gamma = \dfrac{\sqrt{p^2c^2 + m^2c^4}}{mc^2}$. Write $\mathcal{E} = \sqrt{p^2c^2 + m^2c^4}$. Then
$$g = pv^* + \frac{mc^2}{\gamma} = \frac{p^2c^2}{\mathcal{E}} + \frac{m^2c^4}{\mathcal{E}} = \mathcal{E}, \qquad g'(p) = \frac{pc^2}{\mathcal{E}} = v^*. \checkmark$$
This $f$ is the free relativistic Lagrangian, so $g$ is $H = \sqrt{p^2c^2 + m^2c^4}$, rest energy included. Every real $p$ is allowed: $f'(v) \to \pm\infty$ as $v \to \pm c$, so the slopes of $f$ cover $\mathbb{R}$. The finite velocity interval maps onto the whole momentum line. (Project: symbolic check.)

**A2.** (a) $m_1$ descends $x$ while $m_2$ rises $x$: $T = \frac12(m_1 + m_2)\dot x^2$, $V = -m_1gx + m_2gx$. So $L = \frac12(m_1 + m_2)\dot x^2 + (m_1 - m_2)gx$, $p = (m_1 + m_2)\dot x$, and
$$H = \frac{p^2}{2(m_1 + m_2)} - (m_1 - m_2)gx, \qquad \dot x = \frac{p}{m_1 + m_2}, \qquad \dot p = (m_1 - m_2)g.$$
Hence $\ddot x = (m_1 - m_2)g/(m_1 + m_2)$.

(b) On the wire $\dot z = kx\dot x$, so $T = \frac12m(1 + k^2x^2)\dot x^2$ and $V = mgz = \frac12mgkx^2$. Then $p = \partial L/\partial\dot x = m(1 + k^2x^2)\dot x$, which is **not** $m\dot x$: the momentum conjugate to $x$ carries the vertical motion too. The constraint is fixed, so Theorem 3.3 gives
$$H = \frac{p^2}{2m(1 + k^2x^2)} + \frac12mgkx^2, \qquad \dot x = \frac{p}{m(1 + k^2x^2)}, \qquad \dot p = \frac{k^2xp^2}{m(1 + k^2x^2)^2} - mgkx.$$
(For $\dot p$: $\partial_x\big[p^2/2m(1 + k^2x^2)\big] = -k^2xp^2/m(1 + k^2x^2)^2$, and $\dot p = -\partial H/\partial x$.)

**A3.** By Leibniz, $\{q, p^2\} = p\{q, p\} + \{q, p\}p = 2p$, so $\{q^2, p^2\} = q\{q, p^2\} + \{q, p^2\}q = 4qp$.

$\{q, H\} = \frac{1}{2m}\{q, p^2\} + \{q, V(q)\} = \frac{p}{m} + 0$ (two functions of $q$ alone have zero bracket). $\{p, H\} = 0 + \{p, V\} = \frac{\partial p}{\partial q}\frac{\partial V}{\partial p} - \frac{\partial p}{\partial p}\frac{\partial V}{\partial q} = -V'(q)$. These are Hamilton's equations.

For the oscillator, $\{qp, \frac{p^2}{2m}\} = q\{p, \tfrac{p^2}{2m}\} + \{q, \tfrac{p^2}{2m}\}p = 0 + \frac{p^2}{m}$ and $\{qp, \frac12m\omega^2q^2\} = q\{p, \frac12m\omega^2q^2\} + 0 = -m\omega^2q^2$. So
$$\{qp, H\} = \frac{p^2}{m} - m\omega^2q^2 = 2(T - V).$$
By (3.14), $\frac{d}{dt}(qp) = 2(T - V)$. Average over one period $\tau$: the left side averages to $[qp]_0^\tau/\tau = 0$ because the motion is periodic. So $\langle T\rangle = \langle V\rangle = E/2$. (Project: symbolic check.)

**A4.** For one degree of freedom, canonical ⇔ $\det M = 1$ (§3.7).

(a) $M = \begin{pmatrix}\cos\alpha & -\sin\alpha\\ \sin\alpha & \cos\alpha\end{pmatrix}$, $\det = 1$. Canonical (it is the time-$(-\alpha)$ flow of the oscillator with $m = \omega = 1$, a rotation of the phase plane).
(b) $\det = 2\cdot\frac12 = 1$. Canonical.
(c) $M = \begin{pmatrix}2q & 0\\ -p/2q^2 & 1/2q\end{pmatrix}$, $\det = 1$. Canonical: the point transformation $f(q) = q^2$.
(d) $M = \begin{pmatrix}1 & 0\\ 2q & 1\end{pmatrix}$, $\det = 1$. Canonical: a shear in $p$, the same kind of map as a kick in (3.31), and the 1D version of the gauge change $\mathbf{p} \to \mathbf{p} + \nabla\chi$.
(e) $M = \begin{pmatrix}0 & 1\\ 1 & 0\end{pmatrix}$, $\det = -1$. **Not** canonical: $\dot Q = \dot p = -\partial H/\partial q = -\partial K/\partial P$, the wrong sign. Compare $Q = p$, $P = -q$, which is canonical.
(Project: all five determinants.)

**A5.** $V(x) = \frac13x^3 - x$. Fixed points need $p = 0$ and $V'(x) = x^2 - 1 = 0$: $x = \pm1$. The Hessian of $H$ is $\operatorname{diag}(V''(x), 1) = \operatorname{diag}(2x, 1)$.
- $(1, 0)$: $V'' = 2 > 0$, a minimum of $H$, **elliptic**, $V(1) = -\frac23$.
- $(-1, 0)$: $V'' = -2 < 0$, a saddle, **hyperbolic**, $V(-1) = \frac23$.

Separatrix energy $E_s = \frac23$. On $p = 0$: $\frac13x^3 - x = \frac23$, i.e. $x^3 - 3x - 2 = (x + 1)^2(x - 2) = 0$. The double root is the saddle; the separatrix loop crosses the axis again at $x = 2$. Sketch: closed orbits nested around $(1, 0)$, filling the teardrop bounded by the separatrix from $x = -1$ to $x = 2$, with a corner at the saddle. Motion is bounded for $-\frac23 \le E < \frac23$ with the particle inside the loop ($x > -1$). Everything else escapes to $x \to -\infty$, where $V \to -\infty$. (Project: roots and energies.)

**A6.** $E = \frac12m\omega^2A^2 = \frac12\cdot 2\cdot 9\cdot 0.25 = 2.25$ J. By (3.26), $J = 2\pi E/\omega = 1.5\pi = 4.7124$ J·s and $\nu = \omega/2\pi = 0.4775$ Hz. Bohr–Sommerfeld: $n = J/h = 4.7124/6.626\times10^{-34} = 7.112\times10^{33}$. Adjacent levels differ by $\Delta J = h$, a relative change of $1/n$, about $10^{-34}$. No measurement of a lab spring can resolve that, which is why classical mechanics works for springs. (Project prints all four.)

## B. Core

**B1.** (a) Index form, $i = 1, 2, 3$, with $\mathbf{A}$ and $\Phi$ functions of $(\mathbf{r}, t)$:

1. $\dot x_i = \partial H/\partial p_i = (p_i - q_eA_i)/m$, so the kinetic momentum is $m\dot x_i = p_i - q_eA_i$.
2. $\dot p_i = -\partial H/\partial x_i = \sum_j\frac{p_j - q_eA_j}{m}q_e\,\partial_iA_j - q_e\partial_i\Phi = q_e\sum_j\dot x_j\,\partial_iA_j - q_e\partial_i\Phi$.
3. Along the path, $\frac{d}{dt}A_i(\mathbf{r}(t), t) = \partial_tA_i + \sum_j\dot x_j\,\partial_jA_i$. So
$$m\ddot x_i = \dot p_i - q_e\frac{dA_i}{dt} = q_e\Big(-\partial_i\Phi - \partial_tA_i\Big) + q_e\sum_j\dot x_j\big(\partial_iA_j - \partial_jA_i\big).$$
4. The first bracket is $E_i$. For the second, $\mathbf{v}\times(\nabla\times\mathbf{A}) = \nabla(\mathbf{v}\cdot\mathbf{A}) - (\mathbf{v}\cdot\nabla)\mathbf{A}$ with $\mathbf{v}$ held fixed, whose $i$-th component is $\sum_jv_j\partial_iA_j - \sum_jv_j\partial_jA_i$. So $m\ddot{\mathbf{r}} = q_e(\mathbf{E} + \mathbf{v}\times\mathbf{B})$. ∎ (Project: symbolic, for arbitrary $\mathbf{A}$, $\Phi$.)

(b) Both potentials give $\nabla\times\mathbf{A} = B\hat{\mathbf{z}}$.
- Landau gauge $\mathbf{A} = (0, Bx, 0)$: $H = \frac{1}{2m}\big[p_x^2 + (p_y - q_eBx)^2 + p_z^2\big]$. $y$ is cyclic, so $p_y = m\dot y + q_eBx$ is conserved. $m\dot y$ is not: it changes as $x$ oscillates. The conserved $p_y$ fixes the $x$-coordinate of the orbit's center, $x_c = p_y/q_eB$.
- Symmetric gauge $\mathbf{A} = \frac{B}{2}(-y, x, 0)$: $\partial H/\partial y = \frac{p_x - q_eA_x}{m}\cdot\big(-q_e\partial_yA_x\big) = \dot x\,\frac{q_eB}{2}$, so $\dot p_y = -\frac{q_eB}{2}\dot x \ne 0$. Here $p_y$ is not conserved; the gauge is rotationally symmetric, so $p_\phi = xp_y - yp_x$ is instead.

$m\dot y$ is conserved in neither. Canonical momentum depends on the gauge, and its conservation reflects a symmetry of the chosen $\mathbf{A}$, not only of $\mathbf{B}$. The gauge-invariant (measurable) quantity is $m\mathbf{v} = \mathbf{p} - q_e\mathbf{A}$.

**B2.** (a) Explicit Euler: $M_{\text{EE}} = \begin{pmatrix}1 & h\\ -h & 1\end{pmatrix}$. Symplectic Euler: $p' = p - hq$, then $q' = q + hp' = (1 - h^2)q + hp$, so $M_{\text{SE}} = \begin{pmatrix}1 - h^2 & h\\ -h & 1\end{pmatrix}$.

(b) $\det M_{\text{EE}} = 1 + h^2$; $\det M_{\text{SE}} = (1 - h^2) + h^2 = 1$. After $N$ explicit steps the area is multiplied by $(1 + h^2)^N$; for $h = 0.1$, $N = 1000$: $(1.01)^{1000} = 20959.2$ (project).

(c) $M_{\text{EE}} = \sqrt{1 + h^2}\times(\text{rotation})$, so $q'^2 + p'^2 = (1 + h^2)(q^2 + p^2)$: each step multiplies $H$ by $1 + h^2$ and the orbit spirals out. The map is not canonical (it fails (3.20)), so it cannot be the flow of any Hamiltonian and nothing stops the growth. $M_{\text{SE}}$ has $\det = 1$ and $|\operatorname{tr}| = |2 - h^2| < 2$ for $h < 2$, so its eigenvalues lie on the unit circle; it is a canonical map that conserves $\tilde H = \frac12(p^2 + q^2) - \frac{h}{2}qp$ exactly (§3.11), and its orbits are the ellipses $\tilde H = $ const.

**B3.** Use $\{x, H\} = p_x$, $\{p_x, H\} = -x$, and the same for $y$.

(a) Write $H = H_x + H_y$ with $H_x = \frac12(p_x^2 + x^2)$, $H_y = \frac12(p_y^2 + y^2)$. These involve different pairs, so $\{H_x, H_y\} = 0$, and $\{A, H\} = \{H_x - H_y, H_x + H_y\} = 0$.
$\{B, H\} = p_x\{p_y, H\} + \{p_x, H\}p_y + x\{y, H\} + \{x, H\}y = -p_xy - xp_y + xp_y + p_xy = 0$.
$\{L, H\} = x\{p_y, H\} + \{x, H\}p_y - y\{p_x, H\} - \{y, H\}p_x = -xy + p_xp_y + yx - p_yp_x = 0$.

(b) Use $\{F, p_x\} = \partial F/\partial x$ and $\{F, x\} = -\partial F/\partial p_x$ (3.18).
- $\{H_x, B\} = \{H_x, p_x\}p_y + x\{H_x, y\} + \{H_x, x\}y = xp_y + 0 - p_xy = L$, and $\{H_y, B\} = p_x\{H_y, p_y\} + x\{H_y, y\} = p_xy - xp_y = -L$. So $\{A, B\} = L - (-L) = 2L$.
- $\{B, L\}$: $\{p_xp_y, xp_y\} = p_y\{p_x, x\}p_y = -p_y^2$; $\{p_xp_y, -yp_x\} = -p_x\{p_y, y\}p_x = p_x^2$; $\{xy, xp_y\} = x\cdot x\{y, p_y\} = x^2$; $\{xy, -yp_x\} = -y\cdot y\{x, p_x\} = -y^2$. Sum: $p_x^2 - p_y^2 + x^2 - y^2 = 2A$.
- $\{A, L\}$ is the rate of change of $A$ along the rotation generated by $L$: $\dot x = \{x, L\} = -y$, $\dot y = x$, $\dot p_x = \{p_x, L\} = -\partial L/\partial x = -p_y$, $\dot p_y = p_x$. Then $\frac{dA}{d\lambda} = p_x(-p_y) - p_y(p_x) + x(-y) - y(x) = -2B$. So $\{L, A\} = 2B$.

Poisson's theorem says the bracket of two conserved quantities is conserved. Here the brackets are $2L$, $2A$, $2B$, which (a) showed are conserved: consistent, and closed (an $\mathfrak{su}(2)$ algebra up to factors of 2). (Project: all brackets.)

(c) Let $\mathbf{u} = (x, p_x)$ and $\mathbf{w} = (y, p_y)$ in $\mathbb{R}^2$. Then $H = \frac12(|\mathbf{u}|^2 + |\mathbf{w}|^2)$, $A = \frac12(|\mathbf{u}|^2 - |\mathbf{w}|^2)$, $B = \mathbf{u}\cdot\mathbf{w}$, $L = u_1w_2 - u_2w_1$. Lagrange's identity in $\mathbb{R}^2$ gives $B^2 + L^2 = |\mathbf{u}|^2|\mathbf{w}|^2$, so
$$A^2 + B^2 + L^2 = \tfrac14(|\mathbf{u}|^2 - |\mathbf{w}|^2)^2 + |\mathbf{u}|^2|\mathbf{w}|^2 = \tfrac14(|\mathbf{u}|^2 + |\mathbf{w}|^2)^2 = H^2.$$
So among $H, A, B, L$ three are independent: one more than the two needed for integrability. The orbit lies in the intersection of three independent level sets in $\mathbb{R}^4$, a curve, so every orbit is closed (the familiar centered ellipses, both frequencies equal). The quantum counterpart is the $(n + 1)$-fold degeneracy of the 2D oscillator's levels.

**B4.** (a) $p_\theta = \partial L/\partial\dot\theta = mR^2\dot\theta$. Then
$$H = p_\theta\dot\theta - L = \frac{p_\theta^2}{2mR^2} - \frac12mR^2\Omega^2\sin^2\theta - mgR\cos\theta.$$
$T + V = \frac12mR^2(\dot\theta^2 + \Omega^2\sin^2\theta) - mgR\cos\theta$, so $H - (T + V) = -mR^2\Omega^2\sin^2\theta$. $L$ has no explicit $t$, so $H$ is conserved (§3.2). $T + V$ is not: in inertial coordinates the constraint moves (the hoop rotates), Theorem 3.3 does not apply, and the motor that holds $\Omega$ fixed does work on the bead whenever $\theta$ changes. $H$ is the energy in the rotating frame, with the centrifugal potential $-\frac12mR^2\Omega^2\sin^2\theta$.

(b) Fixed points: $p_\theta = 0$ and
$$\frac{\partial H}{\partial\theta} = -mR^2\Omega^2\sin\theta\cos\theta + mgR\sin\theta = mR\sin\theta\,(g - R\Omega^2\cos\theta) = 0.$$
So $\theta = 0$, $\theta = \pi$, and $\cos\theta^* = g/R\Omega^2$ (which exists only if $\Omega^2 > g/R$). The Hessian is $\operatorname{diag}(\partial^2_\theta H, 1/mR^2)$ with $\partial^2_\theta H = mR\cos\theta\,(g - R\Omega^2\cos\theta) + mR^2\Omega^2\sin^2\theta$.
- $\theta = 0$: $\partial^2_\theta H = mR(g - R\Omega^2)$. Elliptic for $\Omega^2 < g/R$, hyperbolic for $\Omega^2 > g/R$.
- $\theta = \pi$: $-mR(g + R\Omega^2) < 0$. Always hyperbolic.
- $\theta^*$: the first term vanishes, leaving $mR^2\Omega^2\sin^2\theta^* > 0$. Elliptic.

At $\Omega^2 = g/R$ the bottom loses stability and two new elliptic points split off at $\pm\theta^*$ (a pitchfork bifurcation). (Project: $H$, $H - (T + V)$, $\partial_\theta H$, $\partial^2_\theta H(0)$.)

**B5.** From the $F_1$ row of the table: $p = \partial F_1/\partial q = m\omega q\cot Q$ and $P = -\partial F_1/\partial Q = \frac12m\omega q^2\csc^2Q$. Solve the second for $q$: $q^2 = \frac{2P}{m\omega}\sin^2Q$, so $q = \sqrt{2P/m\omega}\,\sin Q$. Substitute into the first: $p = m\omega\sqrt{\frac{2P}{m\omega}}\sin Q\cdot\frac{\cos Q}{\sin Q} = \sqrt{2Pm\omega}\cos Q$. That is (3.21). $F_1$ has no $t$, so $K = H = \omega P$ (§3.7). (Project: symbolic.)

**B6.** (a) Free-flight trajectories are $q(t) = q_0 + p_0t/m$, $p(t) = p_0$, and $\rho$ is constant along them (§3.5), so $\rho(q, p, t) = \rho_0(q - pt/m, p)$. Check with (3.11): $\partial_t\rho = -\frac{p}{m}\partial_1\rho_0$ and $\{\rho, H\} = \partial_q\rho\cdot\frac{p}{m} = \frac{p}{m}\partial_1\rho_0$; they cancel.

(b) $z(t) = Mz_0$ is linear, so (taking zero means) $\Sigma(t) = \langle zz^\top\rangle = M\Sigma_0M^\top$. With $\Sigma_0 = \operatorname{diag}(\sigma_q^2, \sigma_p^2)$:
$$\Sigma(t) = \begin{pmatrix}\sigma_q^2 + \sigma_p^2t^2/m^2 & \sigma_p^2t/m\\ \sigma_p^2t/m & \sigma_p^2\end{pmatrix}, \qquad \det\Sigma(t) = \sigma_q^2\sigma_p^2 + \frac{\sigma_p^4t^2}{m^2} - \frac{\sigma_p^4t^2}{m^2} = \sigma_q^2\sigma_p^2.$$
Numbers: $\sigma_q^2(3) = 0.49 + 0.16\times 9 = 1.9300$, $\operatorname{cov}(q, p) = 0.16\times 3 = 0.4800$, $\det\Sigma = 0.49\times 0.16 = 0.0784$ at every $t$ (project).

(c) Liouville fixes the area of the $1\sigma$ ellipse, $\pi\sqrt{\det\Sigma}$, not its width. The ellipse shears: fast particles get ahead, $q$ and $p$ become correlated, and the ellipse gets long in $q$ and thin across. The QM-3 free-particle Gaussian wave packet spreads by the same formula, $\sigma_x^2(t) = \sigma_x^2 + \sigma_p^2t^2/m^2$.

**B7.** (a) $\frac{1}{2m}W'^2 = E$ gives $W = \sqrt{2mE}\,q$ (branch $p > 0$), so $S = \sqrt{2mE}\,q - Et$. Then $\beta = \partial S/\partial E = q\sqrt{m/2E} - t$, so $q = \sqrt{2E/m}\,(t + \beta)$: uniform motion at $v = \sqrt{2E/m}$, with $p = \partial S/\partial q = \sqrt{2mE} = mv$.

(b) $S_q = m(q - q_0)/t$ and $S_t = -m(q - q_0)^2/2t^2$, so $\frac{1}{2m}S_q^2 + S_t = \frac{m(q - q_0)^2}{2t^2} - \frac{m(q - q_0)^2}{2t^2} = 0$. The straight path from $(q_0, 0)$ to $(q, t)$ has velocity $(q - q_0)/t$, so its momentum is $m(q - q_0)/t = S_q$. Its Lagrangian $\frac12m(q - q_0)^2/t^2$ is constant, and times $t$ gives $m(q - q_0)^2/2t = S$. This is a different complete integral, with $q_0$ as the constant: $\beta = \partial S/\partial q_0 = -m(q - q_0)/t = -p$, conserved. $e^{iS/\hbar}$ is the phase of QM-9's free-particle propagator. (Project: both checks.)

**B8.** Code added to `project.py` (Part 13), reusing Part 5's integrators:
```python
def energy_errors(steps, periods):
    hh = T_pend / steps; n = periods * steps
    e_r = rk4_pendulum(th_amp, 0.0, hh, n, n)[-1]          # energy at the end
    e_v = verlet_pendulum(th_amp, 0.0, hh, n, 1)            # energy at every step
    return abs(e_r - E0) / abs(E0), np.max(np.abs(e_v - E0)) / abs(E0)
```
Results over $10^3$ periods (project): RK4 final error $2.639\times10^{-2}$ (50 steps) and $8.347\times10^{-4}$ (100 steps), ratio $31.6$; Verlet maximum error $1.815\times10^{-2}$ and $4.533\times10^{-3}$, ratio $4.00$.

(a) Verlet is second order and is, up to exponentially small terms, the exact flow of $\tilde H = H + O(h^2)$ (§3.11). It conserves $\tilde H$ up to those terms, so $|H - H_0| \le |H - \tilde H| + |\tilde H - H_0| = O(h^2)$ for exponentially long times (far beyond $10^4$ periods): halving $h$ divides the bound by 4. The same maximum appears at $10^4$ periods in Part 5; it does not grow.

(b) For the linear oscillator, one RK4 step multiplies the complex amplitude $q + ip/m\omega$ by $R(z)$ at $z = -i\omega h$, and $|R|$ is the same at $\pm i\omega h$. Write $y = \omega h$. Then $R(iy) = \big(1 - \frac{y^2}{2} + \frac{y^4}{24}\big) + i\big(y - \frac{y^3}{6}\big)$ and
$$|R(iy)|^2 = 1 - y^2 + \frac{y^4}{3} - \frac{y^6}{24} + \frac{y^8}{576} + y^2 - \frac{y^4}{3} + \frac{y^6}{36} = 1 - \frac{y^6}{72} + \frac{y^8}{576}.$$
The $y^2$ and $y^4$ terms cancel exactly. Energy is proportional to $|\text{amplitude}|^2$, so each step multiplies it by about $1 - (\omega h)^6/72$, and after time $t$ ($t/h$ steps) the relative loss is about $t\,\omega^6h^5/72 \propto h^5$. Halving $h$ divides it by $2^5 = 32$. The pendulum is nonlinear, but the measured $31.6$ shows the same scaling. The drift is still linear in $t$: RK4 trades accuracy per step for a slow leak that never stops. (Project: the $|R|^2$ identity symbolically.)

## C. Challenge

**C1.** (a) Energy conservation, $\frac12ml^2\dot\theta^2 - mgl\cos\theta = -mgl\cos\theta_0$, with $\cos\theta = 1 - 2\sin^2\frac{\theta}{2}$, gives $\dot\theta^2 = \frac{4g}{l}\big(k^2 - \sin^2\frac{\theta}{2}\big)$. A quarter period runs from $0$ to $\theta_0$:
$$T = 4\int_0^{\theta_0}\frac{d\theta}{\dot\theta} = 2\sqrt{\frac{l}{g}}\int_0^{\theta_0}\frac{d\theta}{\sqrt{k^2 - \sin^2(\theta/2)}}.$$
Substitute $\sin\frac{\theta}{2} = k\sin\varphi$: $\frac12\cos\frac{\theta}{2}\,d\theta = k\cos\varphi\,d\varphi$, so $d\theta = \frac{2k\cos\varphi\,d\varphi}{\sqrt{1 - k^2\sin^2\varphi}}$, and $\sqrt{k^2 - \sin^2\frac{\theta}{2}} = k\cos\varphi$. The limits become $0$ to $\pi/2$:
$$T = 2\sqrt{\frac{l}{g}}\int_0^{\pi/2}\frac{2k\cos\varphi}{k\cos\varphi\sqrt{1 - k^2\sin^2\varphi}}\,d\varphi = 4\sqrt{\frac{l}{g}}\,\mathrm{K}(k).$$

(b) $p = ml^2\dot\theta = 2ml\sqrt{gl}\,\sqrt{k^2 - \sin^2(\theta/2)} = 2ml\sqrt{gl}\,k\cos\varphi$. The loop is four copies of the quarter, so with $\Delta = \sqrt{1 - k^2\sin^2\varphi}$,
$$J = 4\int_0^{\theta_0}p\,d\theta = 16\,ml\sqrt{gl}\int_0^{\pi/2}\frac{k^2\cos^2\varphi}{\Delta}\,d\varphi.$$
Write $k^2\cos^2\varphi = k^2 - k^2\sin^2\varphi = (k^2 - 1) + \Delta^2$. Then the integral is $(k^2 - 1)\mathrm{K} + \mathrm{E}$, so $J = 16\,ml\sqrt{gl}\,[\mathrm{E} - (1 - k^2)\mathrm{K}]$.
Check $dJ/dE = T$. With the given derivatives,
$$\frac{d}{dk}\big[\mathrm{E} - (1 - k^2)\mathrm{K}\big] = \frac{\mathrm{E} - \mathrm{K}}{k} + 2k\mathrm{K} - \frac{\mathrm{E}}{k} + \frac{(1 - k^2)\mathrm{K}}{k} = -\frac{\mathrm{K}}{k} + 2k\mathrm{K} + \frac{\mathrm{K}}{k} - k\mathrm{K} = k\mathrm{K},$$
so $dJ/dk = 16\,ml\sqrt{gl}\,k\mathrm{K}$. And $E = -mgl\cos\theta_0 = mgl(2k^2 - 1)$, so $dE/dk = 4mglk$. Dividing: $\frac{dJ}{dE} = \frac{16\,ml\sqrt{gl}\,k\mathrm{K}}{4mglk} = 4\sqrt{l/g}\,\mathrm{K} = T$, as (3.25) requires. (Project: $dJ/dk = 16k\mathrm{K}$ with SymPy's elliptic functions, and $dJ/dE = T$ at five energies in Part 8.)

(c) $k'^2 = 1 - k^2 = \frac12(1 - E/mgl) = \frac{E_s - E}{2mgl}$. Then $\mathrm{K} \approx \ln\frac{4}{k'} = \frac12\ln\frac{16}{k'^2} = \frac12\ln\frac{32\,mgl}{E_s - E}$, and
$$T \approx 2\sqrt{\frac{l}{g}}\,\ln\frac{32\,mgl}{E_s - E}.$$
With $m = g = l = 1$ (project): $E = 0.99$: exact $16.1570$, asymptotic $16.1418$. $E = 0.999$: exact $20.7491$, asymptotic $20.7470$. Each factor of 10 closer to $E_s$ adds $2\ln 10 \approx 4.6$ time units: the period diverges, but only logarithmically, because the pendulum lingers near the hyperbolic point (Example 4.3).

**C2.** (a) $p_\phi = \ell$ is constant and $\phi$ advances by $2\pi$ per revolution, so $J_\phi = 2\pi\ell$. For $p_r$, with $E = -|E|$:
$$2m\Big(E + \frac{k}{r}\Big) - \frac{\ell^2}{r^2} = \frac{2m|E|}{r^2}\Big(-r^2 + \frac{k}{|E|}r - \frac{\ell^2}{2m|E|}\Big) = \frac{2m|E|}{r^2}(r - r_1)(r_2 - r),$$
with $r_1 + r_2 = k/|E|$ and $r_1r_2 = \ell^2/2m|E|$. So $J_r = 2\sqrt{2m|E|}\,I$ with $I = \int_{r_1}^{r_2}\frac{\sqrt{(r - r_1)(r_2 - r)}}{r}\,dr$. Put $a = \frac{r_1 + r_2}{2}$, $c = \frac{r_2 - r_1}{2}$, $r = a - c\cos u$, $u \in [0, \pi]$. Then $(r - r_1)(r_2 - r) = c^2\sin^2u$ and $dr = c\sin u\,du$, so $I = \int_0^\pi\frac{c^2\sin^2u}{a - c\cos u}du$. Since $c^2\sin^2u = (a - c\cos u)(a + c\cos u) - (a^2 - c^2)$,
$$I = \int_0^\pi(a + c\cos u)\,du - (a^2 - c^2)\int_0^\pi\frac{du}{a - c\cos u} = \pi a - (a^2 - c^2)\frac{\pi}{\sqrt{a^2 - c^2}} = \pi\Big(\frac{r_1 + r_2}{2} - \sqrt{r_1r_2}\Big),$$
using the standard integral $\int_0^\pi du/(a - c\cos u) = \pi/\sqrt{a^2 - c^2}$ for $a > c \ge 0$ (substitute $s = \tan\frac{u}{2}$). Therefore
$$J_r = 2\sqrt{2m|E|}\,\pi\Big(\frac{k}{2|E|} - \frac{\ell}{\sqrt{2m|E|}}\Big) = \pi k\sqrt{\frac{2m}{|E|}} - 2\pi\ell.$$

(b) $J_r + J_\phi = \pi k\sqrt{2m/|E|}$, so $E = -2\pi^2mk^2/(J_r + J_\phi)^2$. $E$ depends on the actions only through their sum, so
$$\nu_r = \frac{\partial E}{\partial J_r} = \frac{4\pi^2mk^2}{(J_r + J_\phi)^3} = \frac{\partial E}{\partial J_\phi} = \nu_\phi.$$
In one radial period $\phi$ advances by exactly $2\pi$, so the orbit returns to its starting point and closes. For a generic central potential $\nu_r \ne \nu_\phi$ and orbits precess (CM-5's rosettes). Substituting $J_r + J_\phi$: $\nu = \frac{2|E|^{3/2}}{\pi k\sqrt{2m}}$. Kepler's third law with $a = k/2|E|$ gives $T = 2\pi\sqrt{ma^3/k} = \frac{\pi k\sqrt{2m}}{2|E|^{3/2}} = 1/\nu$. ✓ This degeneracy is the classical face of hydrogen's energy depending only on $n$ in QM.

(c) Project, $m = k = 1$, $E = -0.3$, $\ell = 0.8$: $r \in [0.35857, 2.97476]$; $J_r$ numerical $= 3.08501 =$ closed form; radial period from the integrated orbit $13.5193 = 1/\nu$; $\phi$ advances $2.000000\,\pi$ per radial period.

**C3.** (a) $\theta' = \theta + p + K\sin\theta$, $p' = p + K\sin\theta$, so
$$M = \begin{pmatrix}1 + K\cos\theta & 1\\ K\cos\theta & 1\end{pmatrix}, \qquad \det M = 1 + K\cos\theta - K\cos\theta = 1.$$
Canonical. For $H = \frac12p^2 + K\cos\theta\sum_n\delta(t - n)$, the rotor turns freely between kicks ($\dot\theta = p$, $\dot p = 0$), and across a kick $\dot p = -\partial H/\partial\theta = K\sin\theta\,\delta(t - n)$ makes $p$ jump by $K\sin\theta$. Kick then drift for one unit of time is (3.31) with $h = 1$, $V = K\cos\theta$, $T = p^2/2$.

(b) On the torus, a fixed point needs $\sin\theta = 0$ and $p \equiv 0 \pmod{2\pi}$: $(0, 0)$ and $(\pi, 0)$. A determinant-1 matrix has eigenvalues $\lambda$, $1/\lambda$, which lie on the unit circle (elliptic) iff $|\operatorname{tr}M| < 2$. At $(0, 0)$: $\operatorname{tr} = 2 + K > 2$, hyperbolic for every $K > 0$. At $(\pi, 0)$: $\operatorname{tr} = 2 - K$, elliptic iff $0 < K < 4$. ($\theta = \pi$ is the minimum of $V = K\cos\theta$, the "bottom" of the kicked pendulum.)

(c) Project (`standard_map.png`), $\max_n|p_n - p_0|$ over 50 orbits and 2000 kicks: $2.45$ at $K = 0.5$, $4.13$ at $K = 0.97$, $28.71$ at $K = 1.5$ ($2\pi = 6.28$). At small $K$ most orbits lie on invariant circles that wrap around the cylinder, $p \approx f(\theta)$: KAM tori that survived. An invariant closed curve that wraps the cylinder splits it into two pieces the map cannot carry points between, so $p$ is trapped in a band. At $K = 0.97$ islands and chaotic layers have grown, but some wrapping circle survives and confinement holds; Greene (1979) located the breakup of the last one at $K_c = 0.9716$. At $K = 1.5$ none survive, and chaotic orbits diffuse in $p$ without bound. The blocking argument is two-dimensional: with more degrees of freedom tori do not separate phase space, and orbits can leak around them (Arnold diffusion).
