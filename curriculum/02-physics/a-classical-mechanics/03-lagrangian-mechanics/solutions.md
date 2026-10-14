# CM-3 Solutions

Conventions as in the notes: $L = T - V$, angles from the downward vertical, constraint force $= +\lambda\nabla f$, virtual displacements at frozen time, $g = 9.81\ \mathrm{m/s^2}$. Equation numbers (3.x), (4.x) refer to `notes.md`. Every number below is printed or asserted by `project.py`. The symbolic results of A3, A4, B1(c), B3, B4(c), B5, B6, B8(a), C1, C2 and C3 are checked with SymPy in its last cell.

## A. Warm-up

**A1.**

(a) Planar double pendulum: $2\times 2$ planar coordinates minus 2 rod constraints $= 2$. Coordinates $\theta_1, \theta_2$.

(b) Free rigid body: $6$. Center of mass $(X, Y, Z)$ plus three orientation angles (Euler angles, CM-7).

(c) Two masses on a rigid rod in 3D: $6 - 1 = 5$. Center of mass $(X, Y, Z)$ plus the rod's direction as a polar angle and an azimuth. Spinning the rod about its own axis moves no point mass, so it is not a degree of freedom.

(d) Bead on a fixed wire: $3 - 2 = 1$. Arc length $s$ along the wire.

(e) Three masses in a plane joined pairwise by rods: $6 - 3 = 3$. The triangle is rigid: center of mass $(X, Y)$ plus one orientation angle $\psi$.

(f) Disk rolling along a line, upright: center position $x$ and rotation angle $\phi$, tied by $\dot x = R\dot\phi$, which integrates to $x - R\phi = \text{const}$ (§3.2). One degree of freedom, coordinate $x$ (or $\phi$). Holonomic.

(g) Sphere rolling on a plane. A configuration needs the center's $(x, y)$ (its height is fixed at $R$) and three orientation angles: a 5-dimensional configuration space, every point of which is reachable by rolling (the argument of C1(a) adapts). But no-slip says the contact point is at rest, $\mathbf{v}_{\text{cm}} = \boldsymbol\omega\times R\hat{\mathbf{z}}$: two velocity conditions that do not integrate. So at each instant only $5 - 2 = 3$ velocities are independent (the three components of $\boldsymbol\omega$, which then fix $\mathbf{v}_{\text{cm}}$). "Degrees of freedom" can mean the number of coordinates (5) or the number of independent velocities (3). For holonomic systems these agree. For non-holonomic ones they do not, and that mismatch is exactly why the simple form of Lagrange's equations fails for them.

**A2.**

(a) Holonomic, scleronomic: $f = x^2 + y^2 + z^2 - R^2 = 0$.

(b) Not holonomic: an inequality. The particle moves freely inside and feels a force only on contact. While in contact, (a) applies. You handle it by switching regimes at each contact, which is the impact problem of §4.10.

(c) $x\dot x + y\dot y = \frac12\frac{d}{dt}(x^2 + y^2) = 0$, so $x^2 + y^2$ is constant: holonomic, $f = x^2 + y^2 - C = 0$ with $C$ set by the initial position. A velocity condition that is an exact time derivative is a holonomic constraint in disguise.

(d) Non-holonomic: this is the rolling disk of §3.2 with $\theta$ free; C1(a) proves it.

(e) With $\theta_0$ a constant, $\dot x\sin\theta_0 - \dot y\cos\theta_0 = \frac{d}{dt}(x\sin\theta_0 - y\cos\theta_0) = 0$, so $f = x\sin\theta_0 - y\cos\theta_0 - C = 0$: holonomic. The particle is confined to a fixed line of direction $(\cos\theta_0, \sin\theta_0)$. Compare (d) and C1(d): the same expression with $\theta$ a coordinate does not integrate. Whether a velocity constraint integrates depends on what is held fixed.

(f) Holonomic, rheonomic: $f = z - h(t) = 0$ (a particle on a floor that moves up and down on a schedule).

**A3.** The bob is at $\mathbf{r} = (\ell\sin\theta, -\ell\cos\theta)$, so $\partial\mathbf{r}/\partial\theta = \ell(\cos\theta, \sin\theta)$. Keep gravity in $V = -mg\ell\cos\theta$ and treat the push as non-conservative. From (3.6),
$$Q_\theta = F\hat{\mathbf{x}}\cdot\ell(\cos\theta, \sin\theta) = F\ell\cos\theta,$$
a torque: force $F$ times lever arm $\ell\cos\theta$ (the vertical distance from pivot to bob). With $L = \frac12 m\ell^2\dot\theta^2 + mg\ell\cos\theta$, (3.11) gives
$$m\ell^2\ddot\theta + mg\ell\sin\theta = F\ell\cos\theta.$$
Equilibrium ($\ddot\theta = 0$): $mg\sin\theta = F\cos\theta$, so $\tan\theta_{\text{eq}} = F/(mg)$. The bob hangs along the total applied force $(F, -mg)$.

Potential route: a constant force is $-\nabla V_F$ with $V_F = -Fx = -F\ell\sin\theta$. Then $L' = \frac12 m\ell^2\dot\theta^2 + mg\ell\cos\theta + F\ell\sin\theta$ and $\partial L'/\partial\theta = -mg\ell\sin\theta + F\ell\cos\theta$, which gives the same equation. It has to: step 6 of §3.4 showed $Q_i = -\partial V/\partial q_i$ for any conservative force.

**A4.** In terms of $(X, s)$: $X_{\text{wedge}} = X$, $x_{\text{block}} = X + s\cos\alpha$, $y_{\text{block}} = h - s\sin\alpha$. Each row of $J$ holds the derivatives of one Cartesian coordinate:
$$J = \begin{pmatrix}1 & 0\\ 1 & \cos\alpha\\ 0 & -\sin\alpha\end{pmatrix}, \qquad \mathcal{M}J = \begin{pmatrix}M & 0\\ m & m\cos\alpha\\ 0 & -m\sin\alpha\end{pmatrix},$$
$$J^\top\mathcal{M}J = \begin{pmatrix}M + m & m\cos\alpha\\ m\cos\alpha & m\cos^2\alpha + m\sin^2\alpha\end{pmatrix} = \begin{pmatrix}M + m & m\cos\alpha\\ m\cos\alpha & m\end{pmatrix},$$
the matrix of Example 4.6. Its trace is $M + 2m > 0$ and its determinant is $m(M + m) - m^2\cos^2\alpha = m(M + m\sin^2\alpha) > 0$. A symmetric $2\times2$ matrix has eigenvalues whose product is the determinant and whose sum is the trace. Positive product means the two eigenvalues have the same sign; positive sum means that sign is $+$. So the matrix is positive definite for every $\alpha$. This agrees with §3.6 step 4: the two columns of $J$ are independent for every $\alpha$ (the first has no $y$-component, the second has a nonzero $x_{\text{block}}$ or $y_{\text{block}}$ entry).

**A5.** $L = \frac12 m\dot x^2 - \frac12 kx^2$, $\mathcal{F} = \frac12 c\dot x^2$, $E = \frac12 m\dot x^2 + \frac12 kx^2$. Lagrange's equation with $Q = -\partial\mathcal{F}/\partial\dot x = -c\dot x$ is $m\ddot x + kx = -c\dot x$, the given equation. Along a solution,
$$\frac{dE}{dt} = m\dot x\ddot x + kx\dot x = \dot x(m\ddot x + kx) = \dot x(-c\dot x) = -c\dot x^2 = -2\left(\tfrac12 c\dot x^2\right) = -2\mathcal{F}.\ ✓$$
At $t = 0$: $E(0) = \frac12(1)(1)^2 + \frac12(4)(1)^2 = 0.5 + 2 = 2.50$ J and $dE/dt = -c\dot x(0)^2 = -0.5\times1 = -0.50$ W. The position $x(0)$ does not enter: only motion dissipates.

## B. Core

**B1.** (a) In the Lagrangian formalism the independent variables are $(q_1, \dots, q_n, \dot q_1, \dots, \dot q_n, t)$, and $\partial/\partial\dot q_i$ holds every $q_j$, every $\dot q_j$ with $j \ne i$, and $t$ fixed. Along a motion the chain rule gives
$$\dot{\mathbf{r}}(q, \dot q, t) = \sum_j\frac{\partial\mathbf{r}}{\partial q_j}(q, t)\,\dot q_j + \frac{\partial\mathbf{r}}{\partial t}(q, t),$$
and we take this as the definition of $\dot{\mathbf{r}}$ as a function of $(q, \dot q, t)$. The coefficients $\partial\mathbf{r}/\partial q_j$ and $\partial\mathbf{r}/\partial t$ depend on $(q, t)$ only, so $\partial/\partial\dot q_i$ treats them as constants, and $\partial\dot q_j/\partial\dot q_i = \delta_{ij}$. Hence $\partial\dot{\mathbf{r}}/\partial\dot q_i = \sum_j(\partial\mathbf{r}/\partial q_j)\delta_{ij} = \partial\mathbf{r}/\partial q_i$. ∎

(b) Let $\mathbf{g}_i(q, t) := \partial\mathbf{r}/\partial q_i$. Its total time derivative along a motion is
$$\frac{d\mathbf{g}_i}{dt} = \sum_j\frac{\partial\mathbf{g}_i}{\partial q_j}\dot q_j + \frac{\partial\mathbf{g}_i}{\partial t} = \sum_j\frac{\partial^2\mathbf{r}}{\partial q_j\,\partial q_i}\dot q_j + \frac{\partial^2\mathbf{r}}{\partial t\,\partial q_i}.$$
Differentiating the expression for $\dot{\mathbf{r}}$ in (a) with respect to $q_i$ at fixed $\dot q$ and $t$ gives
$$\frac{\partial\dot{\mathbf{r}}}{\partial q_i} = \sum_j\frac{\partial^2\mathbf{r}}{\partial q_i\,\partial q_j}\dot q_j + \frac{\partial^2\mathbf{r}}{\partial q_i\,\partial t}.$$
Because $\mathbf{r}$ is $C^2$, mixed partials commute (Schwarz's theorem), and the two lines agree term by term. ∎ The $\partial\mathbf{r}/\partial t$ term contributes the same mixed partial to both sides. That is what keeps the identity true for a moving constraint.

(c) $\mathbf{r} = R(\sin\theta\cos\omega t,\ \sin\theta\sin\omega t,\ -\cos\theta)$.

- $\dot{\mathbf{r}} = R\dot\theta(\cos\theta\cos\omega t, \cos\theta\sin\omega t, \sin\theta) + R\omega\sin\theta(-\sin\omega t, \cos\omega t, 0)$.
- $\partial\dot{\mathbf{r}}/\partial\dot\theta = R(\cos\theta\cos\omega t, \cos\theta\sin\omega t, \sin\theta)$.
- $\partial\mathbf{r}/\partial\theta = R(\cos\theta\cos\omega t, \cos\theta\sin\omega t, \sin\theta)$. Same vector: (3.7) ✓.
- $\partial\dot{\mathbf{r}}/\partial\theta$ (at fixed $\dot\theta$, $t$) $= R\dot\theta(-\sin\theta\cos\omega t, -\sin\theta\sin\omega t, \cos\theta) + R\omega\cos\theta(-\sin\omega t, \cos\omega t, 0)$.
- $\frac{d}{dt}(\partial\mathbf{r}/\partial\theta) = R\dot\theta(-\sin\theta\cos\omega t, -\sin\theta\sin\omega t, \cos\theta) + R\omega\cos\theta(-\sin\omega t, \cos\omega t, 0)$. The first term comes from $\theta(t)$, the second from the explicit $\omega t$. Same vector: (3.8) ✓.

If you forget the explicit time dependence when computing $d/dt$, you lose the second term and the identity appears to fail. That second term is the $\partial^2\mathbf{r}/\partial t\,\partial\theta$ of part (b).

**B2.** (a) $L = \frac12 m(\dot x^2 + \dot y^2) - mgy$ and $f = x^2 + y^2 - \ell^2$, so $\nabla f = (2x, 2y)$. Equations (3.18):
$$m\ddot x = 2\lambda x, \qquad m\ddot y + mg = 2\lambda y.$$
The constraint force is $\lambda\nabla f = 2\lambda(x, y)$, along the line from the pivot to the bob. For $\lambda > 0$ it points **outward**, away from the pivot. A string can only pull, so for a string $\lambda \le 0$. A rigid rod can also push; $\lambda > 0$ means the rod is in compression (this happens near the top of a slow swing over the top).

(b) Differentiate the constraint twice: $\dot f = 2(x\dot x + y\dot y) = 0$ and $\ddot f = 2(x\ddot x + y\ddot y + \dot x^2 + \dot y^2) = 0$. Substitute $\ddot x = 2\lambda x/m$ and $\ddot y = 2\lambda y/m - g$:
$$\frac{2\lambda}{m}(x^2 + y^2) - gy + \dot x^2 + \dot y^2 = 0 \quad\Longrightarrow\quad \lambda = \frac{m\,(gy - \dot x^2 - \dot y^2)}{2\ell^2},$$
using $x^2 + y^2 = \ell^2$. Polar form: $x = \ell\sin\theta$, $y = -\ell\cos\theta$, $\dot x^2 + \dot y^2 = \ell^2\dot\theta^2$, so
$$\lambda = \frac{m(-g\ell\cos\theta - \ell^2\dot\theta^2)}{2\ell^2} = -\frac{m(g\cos\theta + \ell\dot\theta^2)}{2\ell}.$$
The force is $2\lambda(x, y) = 2\lambda\ell\,\hat{\mathbf{r}} = -m(g\cos\theta + \ell\dot\theta^2)\,\hat{\mathbf{r}}$: inward, with magnitude
$$T = m(g\cos\theta + \ell\dot\theta^2) = -2\ell\lambda.$$
Newton check, radial direction: $T - mg\cos\theta = m\ell\dot\theta^2$ (centripetal). ✓

(c) $T = 1\times(9.81\times1 + 1\times2^2) = 13.81$ N.

(d) The reference implementation is `pendulum_dae_rhs_factory` in `project.py`:

```python
def pendulum_dae_rhs_factory(beta=0.0):
    def rhs(tt, y):
        xx, yy, xd, yd = y
        J = np.array([[2 * xx, 2 * yy]])                  # grad f
        f = xx**2 + yy**2 - ELL_P**2
        fd = 2 * (xx * xd + yy * yd)
        gamma = -2 * (xd**2 + yd**2) - 2 * beta * fd - beta**2 * f
        qdd, lam = dae_accelerations(M_pend, F_pend, J, [gamma])   # M = m I, F = (0, -m g)
        return np.concatenate([[xd, yd], qdd])
    return rhs
```

The constraint row's right side is $-\dot q^\top(\nabla^2 f)\dot q$ with $\nabla^2 f = 2I$, which is $-2(\dot x^2 + \dot y^2)$. Results, from $\theta_0 = 1$ rad at rest over 10 s with RK45 at `rtol=1e-6`: the recorded $\lambda(t)$ matches $-T/(2\ell)$ computed from the integrated angle to $0.1\%$ (the conversion to $\theta$ uses the slightly drifted $x, y$). The maximum constraint error $|f|$ is $8.76\times10^{-5}\ \mathrm{m^2}$ without stabilization and $8.69\times10^{-7}\ \mathrm{m^2}$ with Baumgarte $\beta = 5$, about 101 times smaller.

Why the plain version drifts: it enforces $\ddot f = 0$, so integration errors in $f$ and $\dot f$ are never corrected. Each step's error in $\dot f$ adds a linear-in-time error to $f$, and these accumulate. Baumgarte enforces $\ddot f + 2\beta\dot f + \beta^2 f = 0$ instead, a critically damped equation for $f$ whose solutions decay like $(1 + \beta t)e^{-\beta t}$. Errors are now pulled back toward $f = 0$ on the time scale $1/\beta = 0.2$ s. Choose $\beta$ fast compared to the motion but not so large that the equations become stiff.

**B3.** (a) The string does not slip on the rim, so the string speed equals the rim speed: $\dot x_1 = R\dot\phi$. With $R$ constant this is $\frac{d}{dt}(x_1 - R\phi) = 0$, which integrates to $x_1 - R\phi = \text{const}$. Choose the zero of $\phi$ to make the constant zero. Holonomic. Rolling is non-holonomic only when the rolling direction depends on the configuration (the disk's heading $\theta$ in §3.2). Here the string always meets the rim in the same way.

(b) Coordinates $(x_1, x_2, \phi)$ and $L = \frac12 m_1\dot x_1^2 + \frac12 m_2\dot x_2^2 + \frac12 I\dot\phi^2 + m_1 gx_1 + m_2 gx_2$. Gradients: $\nabla f_1 = (1, 1, 0)$, $\nabla f_2 = (1, 0, -R)$. Equations (3.18):
$$m_1\ddot x_1 - m_1 g = \lambda_1 + \lambda_2, \tag{B3.1}$$
$$m_2\ddot x_2 - m_2 g = \lambda_1, \tag{B3.2}$$
$$I\ddot\phi = -R\lambda_2. \tag{B3.3}$$

(c) Differentiate the constraints twice: $\ddot x_2 = -\ddot x_1$ and $\ddot\phi = \ddot x_1/R$. From (B3.3), $\lambda_2 = -I\ddot x_1/R^2$. From (B3.2), $\lambda_1 = -m_2\ddot x_1 - m_2 g$. Substitute both into (B3.1):
$$m_1\ddot x_1 - m_1 g = -m_2\ddot x_1 - m_2 g - \frac{I}{R^2}\ddot x_1 \quad\Longrightarrow\quad \ddot x_1 = \frac{(m_1 - m_2)\,g}{m_1 + m_2 + I/R^2}.$$
Then
$$\lambda_1 = -m_2(g + \ddot x_1) = -\frac{m_2 g\,(2m_1 + I/R^2)}{m_1 + m_2 + I/R^2}, \qquad \lambda_2 = -\frac{I}{R^2}\ddot x_1 = -\frac{I}{R^2}\,\frac{(m_1 - m_2)\,g}{m_1 + m_2 + I/R^2}.$$

(d) By (3.19), the generalized constraint force on $x_1$ is $\lambda_1\partial f_1/\partial x_1 + \lambda_2\partial f_2/\partial x_1 = \lambda_1 + \lambda_2$, on $x_2$ it is $\lambda_1$, and on $\phi$ it is $-R\lambda_2$ (a torque). The $x$'s point down and the strings pull up, so
$$T_1 = -(\lambda_1 + \lambda_2), \qquad T_2 = -\lambda_1.$$
Newton check: (B3.1) and (B3.2) rearrange to $T_1 = m_1(g - \ddot x_1)$ and $T_2 = m_2(g + \ddot x_1)$. ✓ Then
$$T_1 - T_2 = -\lambda_2 = \frac{I\ddot x_1}{R^2} = \frac{I\ddot\phi}{R}.$$
Read as $(T_1 - T_2)R = I\ddot\phi$, this is the pulley's rotational Newton law. The two string segments pull down on opposite sides of the rim, and the difference in their tensions is the net torque that spins the pulley up. The multipliers divide the labor: $\lambda_1$ enforces inextensibility (minus the tension on $m_2$'s side), and $\lambda_2$ enforces no-slip (minus the tension difference, which is the friction the rim exerts along the string).

(e) As $I\to0$: $\ddot x_1 \to (m_1 - m_2)g/(m_1 + m_2)$, $\lambda_2 \to 0$ (so $T_1 = T_2$: a massless pulley passes tension through unchanged), and $\lambda_1 \to -2m_1m_2g/(m_1 + m_2)$. That is (4.4). ✓ For $I > 0$, the pulley acts like an extra mass $I/R^2$ added to the denominator: both weights must also accelerate the rim.

**B4.** (a) The point at radius $r$ on the rotating wire is $\mathbf{r}(r, t) = (r\cos\omega t,\ r\sin\omega t,\ \frac12 kr^2)$. Then
$$\dot{\mathbf{r}} = (\dot r\cos\omega t - r\omega\sin\omega t,\ \dot r\sin\omega t + r\omega\cos\omega t,\ kr\dot r), \qquad |\dot{\mathbf{r}}|^2 = \dot r^2 + r^2\omega^2 + k^2r^2\dot r^2,$$
because the cross terms in the $x$ and $y$ components cancel. So $T = \frac12 m(1 + k^2r^2)\dot r^2 + \frac12 m\omega^2r^2$ and $V = mgz = \frac12 mgkr^2$, giving
$$L = \tfrac12 m(1 + k^2r^2)\dot r^2 + \tfrac12 m(\omega^2 - gk)r^2. ✓$$
Here $T_2 = \frac12 m(1 + k^2r^2)\dot r^2$, $T_1 = 0$, $T_0 = \frac12 m\omega^2r^2$.

(b) $\partial L/\partial\dot r = m(1 + k^2r^2)\dot r$, with time derivative $m(1 + k^2r^2)\ddot r + 2mk^2r\dot r^2$. Also $\partial L/\partial r = mk^2r\dot r^2 + m(\omega^2 - gk)r$. Subtracting,
$$m(1 + k^2r^2)\ddot r + mk^2r\dot r^2 = m(\omega^2 - gk)\,r, \qquad \ddot r = \frac{(\omega^2 - gk)\,r - k^2r\dot r^2}{1 + k^2r^2}.$$
The $\dot r^2$ term exists because the mass "matrix" $m(1 + k^2r^2)$ depends on $r$ (the fourth bullet of §8). Of the $2mk^2r\dot r^2$ from $d/dt$ and the $mk^2r\dot r^2$ from $\partial L/\partial r$, one survives. Physically, the wire gets steeper as $r$ grows. Moving out at constant $\dot r$ means the vertical speed $kr\dot r$ keeps increasing, which takes a vertical acceleration $k\dot r^2$. The component of that acceleration along the wire is the $k^2r\dot r^2$ term.

(c) Equilibria need $\dot r = 0$ and $\ddot r = 0$, so $(\omega^2 - gk)\,r = 0$.

- If $\omega^2 \ne gk$, the only equilibrium is $r = 0$. Linearizing (the $k^2$ terms are higher order) gives $\ddot r = (\omega^2 - gk)r$. For $\omega^2 < gk$ it is stable, with frequency $\sqrt{gk - \omega^2}$. For $\omega^2 > gk$ it is unstable, with growth rate $\sqrt{\omega^2 - gk}$.
- If $\omega^2 = gk$, every $r$ is an equilibrium: with $\dot r = 0$ the right side vanishes identically. The equilibria form a continuous line, and every one of them is neutral.

Compare the hoop. In both systems $V_{\text{eff}} = V - T_0$, and the bottom loses stability when the centrifugal term overcomes gravity's curvature there. On the hoop, $V_{\text{eff}} = -mgR\cos\theta - \frac12 mR^2\omega^2\sin^2\theta$, and the two pieces have different shapes in $\theta$. Past threshold they balance at a new angle $\cos\theta^* = g/(R\omega^2)$, and two stable branches appear: the pitchfork. Expanding to fourth order,
$$V_{\text{eff}}^{\text{hoop}} = \text{const} + \tfrac12 mR^2(\omega_0^2 - \omega^2)\theta^2 + mR^2\left(\frac{\omega^2}{6} - \frac{\omega_0^2}{24}\right)\theta^4 + \dots,$$
and at $\omega = \omega_0$ the quartic coefficient is $mR^2\omega_0^2/8 > 0$. The new minima sit in that positive quartic term. On the parabola, $V_{\text{eff}} = \frac12 m(gk - \omega^2)r^2$ is exactly quadratic: gravity's pull along the wire and the centrifugal push are both linear in $r$ with the same coefficient shape, so they cancel everywhere at once or nowhere. There is no quartic term and no new branch. Below threshold there is one stable equilibrium, at threshold a whole line of them, and above it none that is stable. This is a degenerate bifurcation, not a pitchfork.

The project's numerical version uses $k = 2\ \mathrm{m^{-1}}$, so $\omega_{\text{crit}} = \sqrt{gk} = 4.4294$ rad/s. At $\omega = \omega_{\text{crit}}/2$, a small oscillation ($r_0 = 0.02$ m) has period $1.6386$ s against $2\pi/\sqrt{gk - \omega^2} = 1.6379$ s. At $r_0 = 0.3$ m the period is $1.7767$ s, longer because the effective mass $m(1 + k^2r^2)$ grows away from the axis. At $\omega = \omega_{\text{crit}}$, a bead released at rest at $r = 0.3$ m stays there exactly. At $\omega = 1.5\,\omega_{\text{crit}}$ it reaches $r = 14.3$ m by $t = 6$ s.

(d) For $\omega^2 > gk$ the bead slides outward and up the wire, at first exponentially ($r \propto e^{\sqrt{\omega^2 - gk}\,t}$), and nothing stops it: it leaves the end of the wire. The conserved quantity is $h = T_2 - T_0 + V = \frac12 m(1 + k^2r^2)\dot r^2 + \frac12 m(gk - \omega^2)r^2$. As $r$ grows its second term goes to $-\infty$, so the first must grow to compensate. Meanwhile $T + V = h + 2T_0 = h + m\omega^2r^2$ grows without bound. The motor turning the wire is supplying that energy through the wire's normal force.

**B5.** With $\ell_1 = \ell_2 = \ell$, the mass matrix of (4.9) at $\theta_1 = \theta_2 = 0$ and the Hessian of $V$ there are
$$M_0 = \ell^2\begin{pmatrix}m_1 + m_2 & m_2\\ m_2 & m_2\end{pmatrix}, \qquad K = g\ell\begin{pmatrix}m_1 + m_2 & 0\\ 0 & m_2\end{pmatrix}.$$
(Linearizing (4.10) directly gives the same thing: $\cos\Delta\to1$, $\sin\theta_i\to\theta_i$, and the $\dot\theta^2\sin\Delta$ terms are second order.) Put $\Omega^2 = \mu g/\ell$, so $\Omega^2\ell^2 = \mu g\ell$ and
$$K - \Omega^2M_0 = g\ell\begin{pmatrix}(m_1 + m_2)(1 - \mu) & -m_2\mu\\ -m_2\mu & m_2(1 - \mu)\end{pmatrix}.$$
Setting the determinant to zero gives $(m_1 + m_2)m_2(1 - \mu)^2 - m_2^2\mu^2 = 0$, that is,
$$(1 - \mu)^2 = \rho\,\mu^2, \qquad \rho := \frac{m_2}{m_1 + m_2} \in (0, 1).$$
Taking square roots, $1 - \mu = \pm\sqrt\rho\,\mu$, so $\mu = 1/(1 \pm \sqrt\rho)$. The smaller root is $\mu_- = 1/(1 + \sqrt\rho)$ and the larger is $\mu_+ = 1/(1 - \sqrt\rho)$: the stated formula. ✓

Check $m_1 = m_2$: $\rho = \frac12$, and $\mu = 1/(1 \pm 1/\sqrt2) = \sqrt2/(\sqrt2 \pm 1) = \sqrt2(\sqrt2 \mp 1) = 2 \mp \sqrt2$. So $\mu_- = 0.5858$ and $\mu_+ = 3.4142$, as in Example 4.4. ✓

Mode shapes, from the first row: $(m_1 + m_2)(1 - \mu)a_1 = m_2\mu a_2$, so $a_2/a_1 = (1 - \mu)/(\rho\mu) = \pm\sqrt\rho/\rho = \pm1/\sqrt\rho$. The slow mode is in phase ($+$), the fast mode is in antiphase ($-$). For $m_1 = m_2$ this is $\pm\sqrt2$. ✓

Limits.
- $m_2\to0$: $\rho\to0$ and $\mu_\pm\to1$, so both frequencies tend to $\sqrt{g/\ell}$. The heavy upper bob swings as a simple pendulum, unaffected by the light one. The light lower bob is itself a pendulum of the same length hanging from that moving support, so it is driven exactly at resonance. The frequencies merge, and $a_2/a_1 = \pm1/\sqrt\rho\to\pm\infty$: in both modes the light bob swings far more than the heavy one.
- $m_1\to0$: $\rho\to1$, so $\mu_-\to\frac12$ and $\mu_+\to\infty$. The slow mode has $\Omega^2 = g/(2\ell)$ with $a_2/a_1\to1$: the chain moves straight, as a simple pendulum of length $2\ell$, because a massless middle joint cannot resist bending. The fast mode's frequency diverges. Bending the chain costs potential energy (gravity on $m_2$ restores it), but the joint has no inertia, so the kink relaxes infinitely fast. The system is left with one effective degree of freedom. This matches the straight-pendulum limit in Example 4.4, which is consistent only for $m_1 = 0$.

**B6.** Write $p := p_\phi$. From (4.16), $V_{\text{eff}}(\theta) = \dfrac{p^2}{2m\ell^2\sin^2\theta} - mg\ell\cos\theta$. Then
$$V_{\text{eff}}' = -\frac{p^2\cos\theta}{m\ell^2\sin^3\theta} + mg\ell\sin\theta, \qquad V_{\text{eff}}'' = \frac{p^2}{m\ell^2}\,\frac{\sin^2\theta + 3\cos^2\theta}{\sin^4\theta} + mg\ell\cos\theta,$$
using $\frac{d}{d\theta}\frac{\cos\theta}{\sin^3\theta} = \frac{-\sin\theta\cdot\sin^3\theta - \cos\theta\cdot3\sin^2\theta\cos\theta}{\sin^6\theta} = -\frac{\sin^2\theta + 3\cos^2\theta}{\sin^4\theta}$.

At the conical equilibrium, (4.15) gives $p = m\ell^2\sin^2\theta_0\,\dot\phi_0$, so $p^2/(m\ell^2\sin^4\theta_0) = m\ell^2\dot\phi_0^2$. The equilibrium condition $V_{\text{eff}}'(\theta_0) = 0$ gave $\dot\phi_0^2 = g/(\ell\cos\theta_0)$, i.e. $mg\ell\cos\theta_0 = m\ell^2\dot\phi_0^2\cos^2\theta_0$. Therefore
$$V_{\text{eff}}''(\theta_0) = m\ell^2\dot\phi_0^2(\sin^2\theta_0 + 3\cos^2\theta_0) + m\ell^2\dot\phi_0^2\cos^2\theta_0 = m\ell^2\dot\phi_0^2(1 + 3\cos^2\theta_0).$$
The $\theta$ equation is $m\ell^2\ddot\theta = -V_{\text{eff}}'(\theta)$. With $\theta = \theta_0 + \epsilon$ it linearizes to $m\ell^2\ddot\epsilon = -V_{\text{eff}}''(\theta_0)\,\epsilon$, so
$$\Omega^2 = \frac{V_{\text{eff}}''(\theta_0)}{m\ell^2} = \dot\phi_0^2(1 + 3\cos^2\theta_0) = \frac{g}{\ell\cos\theta_0}(1 + 3\cos^2\theta_0). \qquad ∎$$

Bounds: for $0 < \theta_0 < \pi/2$, $0 < \cos^2\theta_0 < 1$, so $1 < 1 + 3\cos^2\theta_0 < 4$ and $1 < \Omega/\dot\phi_0 < 2$, strictly.

Precession: seen from above, the bob's distance from the axis, $\ell\sin\theta$, oscillates at $\Omega$ while the azimuth advances at about $\dot\phi_0$. Between successive outermost points (one $\theta$-period $2\pi/\Omega$), the azimuth advances by $\Delta\phi = 2\pi\dot\phi_0/\Omega$, which lies strictly between $\pi$ and $2\pi$. A closed ellipse centered on the axis, which is the small-amplitude limit $\theta_0\to0$ where $\Omega\to2\dot\phi_0$ (the 2D isotropic oscillator), has its outermost points exactly $\pi$ apart. Here they are more than $\pi$ apart, so each successive outermost point is rotated forward by $\Delta\phi - \pi$: the near-ellipse precesses in the direction of motion. It closes only if $\Delta\phi/2\pi$ happens to be rational.

Numbers ($\ell = 1$ m, $\theta_0 = 40°$): $\cos\theta_0 = 0.76604$, $\dot\phi_0^2 = g/(\ell\cos\theta_0) = 12.8060$, $\dot\phi_0 = 3.5786$ rad/s. $1 + 3\cos^2\theta_0 = 2.76047$, so $\Omega = 3.5786\times\sqrt{2.76047} = 5.9456$ rad/s and $\Omega/\dot\phi_0 = 1.6615$. The $\theta$-period is $2\pi/\Omega = 1.0568$ s (the project measures $1.0578$ s on a perturbed conical orbit; the difference is the amplitude dependence). The azimuth advances $216.7°$ per $\theta$-period, so the outermost point precesses by $36.7°$ each time.

**B7.** (a) **Euler's theorem on homogeneous functions.** If $\mathcal{F}(q, s\dot q) = s^k\mathcal{F}(q, \dot q)$ for all $s > 0$, then $\sum_i\dot q_i\,\partial\mathcal{F}/\partial\dot q_i = k\mathcal{F}$. *Proof:* differentiate both sides with respect to $s$. The left side gives $\sum_i(\partial\mathcal{F}/\partial\dot q_i)(q, s\dot q)\,\dot q_i$ and the right gives $ks^{k-1}\mathcal{F}(q, \dot q)$. Set $s = 1$. ∎ Rayleigh's $\mathcal{F} = \frac12\sum_{ij}c_{ij}(q)\dot q_i\dot q_j$ has $k = 2$.

Now let $h = \sum_i\dot q_ip_i - L$ with $p_i = \partial L/\partial\dot q_i$. Then
$$\frac{dh}{dt} = \sum_i(\ddot q_ip_i + \dot q_i\dot p_i) - \sum_i\left(\frac{\partial L}{\partial q_i}\dot q_i + \frac{\partial L}{\partial\dot q_i}\ddot q_i\right) - \frac{\partial L}{\partial t}.$$
The $\ddot q_ip_i$ terms cancel, $\partial L/\partial t = 0$ by hypothesis, and Lagrange's equations (3.11) give $\dot p_i - \partial L/\partial q_i = Q_i^{\text{nc}} = -\partial\mathcal{F}/\partial\dot q_i$. So
$$\frac{dh}{dt} = \sum_i\dot q_iQ_i^{\text{nc}} = -\sum_i\dot q_i\frac{\partial\mathcal{F}}{\partial\dot q_i} = -2\mathcal{F}. \qquad ∎$$
For scleronomic constraints $h = T + V = E$ (§3.7), and $\mathcal{F}\ge0$ because $c_{ij}$ is a Gram matrix (§3.9 step 3). So the energy never increases.

(b) $\mathcal{F} = \frac12 c\left(|\dot{\mathbf{r}}_1|^2 + |\dot{\mathbf{r}}_2|^2\right)$. From Example 4.4, $|\dot{\mathbf{r}}_1|^2 = \ell_1^2\dot\theta_1^2$ and $|\dot{\mathbf{r}}_2|^2 = \ell_1^2\dot\theta_1^2 + \ell_2^2\dot\theta_2^2 + 2\ell_1\ell_2\dot\theta_1\dot\theta_2\cos\Delta$, so
$$\mathcal{F} = \tfrac12(2c)\ell_1^2\dot\theta_1^2 + \tfrac12 c\,\ell_2^2\dot\theta_2^2 + c\,\ell_1\ell_2\dot\theta_1\dot\theta_2\cos\Delta.$$
This is $T$ of (4.9) with $(m_1 + m_2, m_2)\to(2c, c)$: the same Gram structure with drag coefficients as the weights. The generalized forces are
$$Q_{\theta_1} = -\frac{\partial\mathcal{F}}{\partial\dot\theta_1} = -\left(2c\,\ell_1^2\dot\theta_1 + c\,\ell_1\ell_2\dot\theta_2\cos\Delta\right), \qquad Q_{\theta_2} = -\frac{\partial\mathcal{F}}{\partial\dot\theta_2} = -\left(c\,\ell_2^2\dot\theta_2 + c\,\ell_1\ell_2\dot\theta_1\cos\Delta\right).$$

(c) In `project.py`:

```python
F_dp = T_dp.subs({m1: c, m2: c})                                   # (m1+m2, m2) -> (2c, c)
Q_dp = [-sp.diff(F_dp, th1.diff(t)), -sp.diff(F_dp, th2.diff(t))]
dpd = lagrangian_to_eom(L_dp, [th1, th2], t, Q=Q_dp, params={**dp_params, c: 0.2})
sol = solve_ivp(dpd.rhs, (0, 60), [2.0, 2.0, 0.0, 0.0], t_eval=tt, rtol=1e-10, atol=1e-12)
```

Then evaluate $E = h$ (from `energy_function`) and $\mathcal{F}$ on the solution and compare `np.gradient(E, tt)` with $-2\mathcal{F}$. Results for $m_i = \ell_i = 1$ and $c = 0.2$ kg/s: $E(0) = V = -3g\cos2 = 12.247$ J. The resting value is $V_{\min} = -3g = -29.43$ J, so the pendulum starts $41.68$ J above rest. At $t = 60$ s, $E = -29.430$ J, only $2.85\times10^{-4}$ J above rest: essentially everything has been dissipated. The finite-difference $dE/dt$ matches $-2\mathcal{F}$ to $6.3\times10^{-5}$ of $\max|2\mathcal{F}|$. The plot shows the decay is not a clean exponential: $E$ falls in steps, fastest when the bobs are moving fastest, because $\mathcal{F}$ is quadratic in the velocities.

**B8.** Throughout, $\boldsymbol\omega = \omega\hat{\mathbf{z}}$ is constant and $V = 0$. With $\mathcal{R}(t)$ the rotation by $\omega t$ about $z$, $\mathbf{r}' = \mathcal{R}^\top\mathbf{r}$. Useful fact: $\hat{\mathbf{z}}\times(u, v, 0) = (-v, u, 0)$.

(a) The inertial position is fixed at $(r_0, 0, 0)$, so $\mathbf{r}' = (r_0\cos\omega t,\ -r_0\sin\omega t,\ 0)$. The rotating observer sees the particle circle backwards at rate $\omega$. Differentiating,
$$\dot{\mathbf{r}}' = r_0\omega(-\sin\omega t, -\cos\omega t, 0), \qquad \ddot{\mathbf{r}}' = -r_0\omega^2(\cos\omega t, -\sin\omega t, 0) = -\omega^2\mathbf{r}'.$$
Coriolis: $\boldsymbol\omega\times\dot{\mathbf{r}}' = r_0\omega^2(\cos\omega t, -\sin\omega t, 0) = \omega^2\mathbf{r}'$, so $-2m\boldsymbol\omega\times\dot{\mathbf{r}}' = -2m\omega^2\mathbf{r}'$, pointing inward.
Centrifugal: $\boldsymbol\omega\cdot\mathbf{r}' = 0$, so $\boldsymbol\omega\times(\boldsymbol\omega\times\mathbf{r}') = (\boldsymbol\omega\cdot\mathbf{r}')\boldsymbol\omega - \omega^2\mathbf{r}' = -\omega^2\mathbf{r}'$, and $-m\boldsymbol\omega\times(\boldsymbol\omega\times\mathbf{r}') = +m\omega^2\mathbf{r}'$, pointing outward.
With $V = 0$ and $\dot{\boldsymbol\omega} = 0$, (4.20) reads $m\ddot{\mathbf{r}}' = -2m\omega^2\mathbf{r}' + m\omega^2\mathbf{r}' = -m\omega^2\mathbf{r}'$, which matches $\ddot{\mathbf{r}}' = -\omega^2\mathbf{r}'$. ✓ The sum is the inward force $m\omega^2r_0$ that the apparent circular motion requires. Nothing real acts on the particle: the Coriolis force ($2m\omega^2r_0$ inward) beats the centrifugal force ($m\omega^2r_0$ outward), and the net is exactly centripetal.

(b) At rest in the rotating frame at $(r_0, 0, 0)$, the particle moves on a circle of radius $r_0$ at angular speed $\omega$ in the inertial frame. So the real force must be centripetal, $\mathbf{F} = -m\omega^2r_0\,\hat{\mathbf{r}}$, supplied by a string for example. In the rotating frame, $\ddot{\mathbf{r}}' = \dot{\mathbf{r}}' = \mathbf{0}$, the Coriolis term vanishes, and (4.20) with a real force $\mathbf{F}'$ added reads $\mathbf{0} = \mathbf{F}' + m\omega^2\mathbf{r}'$. So $\mathbf{F}' = -m\omega^2\mathbf{r}'$, the same force expressed in rotating components. ✓

(c) The Coriolis force $-2m\boldsymbol\omega\times\dot{\mathbf{r}}'$ is perpendicular to $\dot{\mathbf{r}}'$, so its power $\mathbf{F}\cdot\dot{\mathbf{r}}'$ is zero. It bends rotating-frame trajectories without changing rotating-frame speeds.

## C. Challenge

**C1.** (a) The configuration space is $\{(x, y, \theta, \phi)\}$. There are two motions the constraints allow. **Roll:** $\theta$ fixed, the contact point moves a signed distance $d$ along $(\cos\theta, \sin\theta)$, and $\phi$ changes by $d/R$. **Pivot:** $x, y, \phi$ fixed and $\theta$ changes. The contact point does not move, so $\dot x = \dot y = \dot\phi = 0$ satisfies both constraints.

*Claim:* from any $A = (x_0, y_0, \theta_0, \phi_0)$ you can reach any $B = (x_1, y_1, \theta_1, \phi_1)$.
1. Pivot to face $(x_1, y_1)$, roll the distance to it, and pivot to $\theta_1$. You are now at $(x_1, y_1, \theta_1, \phi')$ for some $\phi'$.
2. Let $\Delta = \phi_1 - \phi'$ and $d = R\Delta/4$. Roll $d$, pivot $+90°$, and repeat four times. The four displacements are $d$ times four unit vectors $90°$ apart, so they sum to zero. $\theta$ turns through $360°$ and is back where it started. $\phi$ changes by $4d/R = \Delta$. If $\Delta < 0$, take $d < 0$ and roll backwards.

Now suppose the constraints were holonomic, that is, equivalent to $f_k(x, y, \theta, \phi) = \text{const}$ for some non-constant functions $f_k$. Every allowed motion would keep each $f_k$ constant, because $df_k/dt = \nabla f_k\cdot\dot q = 0$ for every allowed $\dot q$. By the claim, every $B$ connects to $A$ through allowed motions, so $f_k(B) = f_k(A)$ for all $B$. Then $f_k$ is constant on the whole configuration space, a contradiction. ∎ Put geometrically: two holonomic constraints would confine the disk to a 2-dimensional surface, but the reachable set is all of the 4-dimensional space.

(b) Take $q = (x, y, \theta, \phi)$ with the given $T$ and $V = 0$. The constraints are $\mathbf{a}_1\cdot\dot q = \dot x - R\dot\phi\cos\theta = 0$ and $\mathbf{a}_2\cdot\dot q = \dot y - R\dot\phi\sin\theta = 0$. The Lagrange–d'Alembert equations, (3.18) with $\mathbf{a}_k$ in place of $\nabla f_k$, are
$$m\ddot x = \lambda_1, \quad m\ddot y = \lambda_2, \quad I_\theta\ddot\theta = 0, \quad I_\phi\ddot\phi = -R\cos\theta\,\lambda_1 - R\sin\theta\,\lambda_2,$$
together with the two constraints. That makes six equations for $x, y, \theta, \phi, \lambda_1, \lambda_2$. The derivation of §3.8 goes through unchanged: d'Alembert restricts $\delta q$ to $\mathbf{a}_k\cdot\delta q = 0$, and step 2 never used that $\mathbf{a}_k$ was a gradient.

(c) The third equation gives $\dot\theta = \Omega$, a constant. Differentiating the constraints gives $\ddot x = R\ddot\phi\cos\theta - R\dot\phi\dot\theta\sin\theta$ and $\ddot y = R\ddot\phi\sin\theta + R\dot\phi\dot\theta\cos\theta$. Put these, with $\lambda_1 = m\ddot x$ and $\lambda_2 = m\ddot y$, into the $\phi$ equation:
$$I_\phi\ddot\phi = -mR(\ddot x\cos\theta + \ddot y\sin\theta) = -mR\left[R\ddot\phi(\cos^2\theta + \sin^2\theta) + R\dot\phi\dot\theta(-\sin\theta\cos\theta + \sin\theta\cos\theta)\right] = -mR^2\ddot\phi.$$
So $(I_\phi + mR^2)\ddot\phi = 0$, which gives $\ddot\phi = 0$. The disk rolls at constant speed $v = R\dot\phi$ while its heading turns at the constant rate $\Omega$. The contact point traces a circle of radius $\rho = v/\Omega$, or a line if $\Omega = 0$.

The multipliers are $\lambda_1 = m\ddot x = -mR\dot\phi\dot\theta\sin\theta$ and $\lambda_2 = m\ddot y = mR\dot\phi\dot\theta\cos\theta$, that is,
$$(\lambda_1, \lambda_2) = mv\Omega\,(-\sin\theta, \cos\theta).$$
This is perpendicular to the heading $(\cos\theta, \sin\theta)$ with magnitude $mv\Omega = mv^2/\rho$: exactly the centripetal force that keeps the center of mass on its circle. Only the ground can supply it, so it is the **sideways static friction** at the contact point. Its generalized component on $\phi$ is $-R(\lambda_1\cos\theta + \lambda_2\sin\theta) = 0$, so it does not brake the rolling, consistent with $\ddot\phi = 0$. It also does no work, because the contact point is at rest, consistent with d'Alembert.

(d) Take $q = (x, y, \theta)$, $\mathbf{a} = (\sin\theta, -\cos\theta, 0)$, and $L = \frac12 m(\dot x^2 + \dot y^2) + \frac12 I\dot\theta^2 + mgx\sin\beta$. The Lagrange–d'Alembert equations are
$$m\ddot x - mg\sin\beta = \lambda\sin\theta, \qquad m\ddot y = -\lambda\cos\theta, \qquad I\ddot\theta = 0.$$
So $\ddot\theta = 0$ and the heading turns at its initial constant rate. The constraint force $\lambda(\sin\theta, -\cos\theta)$ is perpendicular to the blade and acts at the center of mass, so it exerts no torque about the vertical.

The substitution: put $\dot y = \dot x\tan\theta$ into $L$ to get $\tilde L = \frac12 m\dot x^2\sec^2\theta + \frac12 I\dot\theta^2 + mgx\sin\beta$. The Euler–Lagrange equation for $\theta$ is
$$\frac{d}{dt}(I\dot\theta) = \frac{\partial\tilde L}{\partial\theta} = m\dot x^2\sec^2\theta\tan\theta \quad\Longrightarrow\quad I\ddot\theta = m\dot x^2\sec^2\theta\tan\theta \ne 0$$
whenever the skate moves at an angle. This predicts a torque from nowhere: no force in the problem acts off the center of mass.

Which step is violated: §3.8 computes the brackets $E_i = \frac{d}{dt}\frac{\partial T}{\partial\dot q_i} - \frac{\partial T}{\partial q_i} - Q_i$ from the **unconstrained** $L(x, y, \theta, \dot x, \dot y, \dot\theta)$ and lets the constraint restrict only $\delta q$ (step 1). Substitution changes the function being differentiated: $\partial\tilde L/\partial\theta = \partial L/\partial\theta + (\partial L/\partial\dot y)(\partial\dot y/\partial\theta) = 0 + m\dot y\,\dot x\sec^2\theta$, a spurious term from differentiating the constraint relation itself. For a holonomic constraint the same substitution is legitimate. There $y = y(x, \theta)$ exists, $\mathbf{r}(x, \theta)$ is a genuine parametrization, and the identities (3.7)–(3.8) hold for it, so §3.4 applies to the reduced coordinates. For $\dot y = \dot x\tan\theta$ there is no such function: it would need $\partial y/\partial x = \tan\theta$ and $\partial y/\partial\theta = 0$, which contradict each other because $\partial(\tan\theta)/\partial\theta \ne 0$ (mixed partials would disagree). So the interchange identity (3.8), on which (3.9) and the whole of §3.4 rest, has no $\mathbf{r}(x, \theta)$ to apply to, and d'Alembert's principle does not produce the Euler–Lagrange equations of $\tilde L$. A velocity constraint restricts $\delta q$, not $L$.

**C2.** (a) The pivot is at $(0, a\cos\gamma t)$, so the bob is at $x = \ell\sin\theta$, $y = a\cos\gamma t - \ell\cos\theta$, with $\dot x = \ell\dot\theta\cos\theta$ and $\dot y = -a\gamma\sin\gamma t + \ell\dot\theta\sin\theta$. Then
$$T = \tfrac12 m\left(\ell^2\dot\theta^2 - 2a\gamma\ell\,\dot\theta\sin\gamma t\sin\theta + a^2\gamma^2\sin^2\gamma t\right), \qquad V = mgy = mga\cos\gamma t - mg\ell\cos\theta.$$
Drop the pure functions of $t$ ($\frac12 ma^2\gamma^2\sin^2\gamma t$ and $mga\cos\gamma t$):
$$L = \tfrac12 m\ell^2\dot\theta^2 - ma\gamma\ell\,\dot\theta\sin\gamma t\sin\theta + mg\ell\cos\theta.$$
$\partial L/\partial\dot\theta = m\ell^2\dot\theta - ma\gamma\ell\sin\gamma t\sin\theta$, with time derivative $m\ell^2\ddot\theta - ma\gamma^2\ell\cos\gamma t\sin\theta - ma\gamma\ell\,\dot\theta\sin\gamma t\cos\theta$. Also $\partial L/\partial\theta = -ma\gamma\ell\,\dot\theta\sin\gamma t\cos\theta - mg\ell\sin\theta$. Subtracting, the $\dot\theta\sin\gamma t\cos\theta$ terms cancel:
$$m\ell^2\ddot\theta - ma\gamma^2\ell\cos\gamma t\sin\theta + mg\ell\sin\theta = 0 \quad\Longrightarrow\quad \ddot\theta + \frac1\ell\left(g - a\gamma^2\cos\gamma t\right)\sin\theta = 0. ✓$$
In the pivot's frame the effective gravity is $g + \ddot y_p = g - a\gamma^2\cos\gamma t$.

(b) About $\theta = 0$, $\sin\theta\approx\theta$ gives $\ddot\theta + (\omega_0^2 - (a\gamma^2/\ell)\cos\gamma t)\,\theta = 0$ with $\omega_0^2 = g/\ell$. With $\tau = \gamma t/2$ this becomes $\theta'' + (\delta - 2q\cos2\tau)\theta = 0$, where $\delta = 4\omega_0^2/\gamma^2$ and $q = 2a/\ell$: the Mathieu equation in standard form. About $\theta = \pi + \epsilon$, $\sin(\pi + \epsilon) = -\sin\epsilon\approx-\epsilon$ gives
$$\ddot\epsilon - \left(\omega_0^2 - \frac{a\gamma^2}{\ell}\cos\gamma t\right)\epsilon = 0,$$
which is Mathieu with $\delta = -4\omega_0^2/\gamma^2 < 0$ and $q = -2a/\ell$ (the sign of $q$ does not matter: shifting $\tau$ by $\pi/2$ flips it). Negative $\delta$ means the undriven inverted pendulum is unstable. Kapitza's question is whether the drive can stabilize it.

(c) Write $\theta = \Theta + \xi$, where $\Theta$ varies on the slow time scale $1/\omega_0$ and $\xi$ is small and fast (time scale $1/\gamma$). To first order in $\xi$,
$$\ddot\Theta + \ddot\xi = -\frac1\ell\left(g - a\gamma^2\cos\gamma t\right)(\sin\Theta + \xi\cos\Theta).$$
*Fast part.* The largest fast term on the right is $(a\gamma^2/\ell)\cos\gamma t\sin\Theta$, large because $\gamma$ is large. Balance it against $\ddot\xi$, holding $\Theta$ fixed over one fast period: $\ddot\xi = (a\gamma^2/\ell)\cos\gamma t\sin\Theta$, so
$$\xi = -\frac{a}{\ell}\cos\gamma t\,\sin\Theta.$$
Its size is at most $a/\ell \ll 1$, which is consistent with the expansion. The terms dropped from the fast balance, $-(g/\ell)\xi\cos\Theta$ and $(a\gamma^2/\ell)\cos\gamma t\,\xi\cos\Theta$, are smaller than the kept one by factors $\omega_0^2/\gamma^2$ and $a/\ell$.

*Slow part.* Average over one fast period. $\langle\xi\rangle = \langle\ddot\xi\rangle = \langle\cos\gamma t\rangle = 0$, but the product of two fast oscillations does not average to zero:
$$\ddot\Theta = -\frac g\ell\sin\Theta + \frac{a\gamma^2}{\ell}\cos\Theta\,\langle\xi\cos\gamma t\rangle = -\frac g\ell\sin\Theta - \frac{a^2\gamma^2}{\ell^2}\sin\Theta\cos\Theta\,\langle\cos^2\gamma t\rangle = -\frac g\ell\sin\Theta - \frac{a^2\gamma^2}{2\ell^2}\sin\Theta\cos\Theta.$$
Multiply by $m\ell^2$: $m\ell^2\ddot\Theta = -dV_{\text{eff}}/d\Theta$ with
$$V_{\text{eff}}(\Theta) = -mg\ell\cos\Theta + \frac{ma^2\gamma^2}{4}\sin^2\Theta,$$
since $\frac{d}{d\Theta}\frac{ma^2\gamma^2}{4}\sin^2\Theta = \frac{ma^2\gamma^2}{2}\sin\Theta\cos\Theta$. ✓ Now $V_{\text{eff}}'(\pi) = 0$ and
$$V_{\text{eff}}''(\Theta) = mg\ell\cos\Theta + \frac{ma^2\gamma^2}{2}\cos2\Theta, \qquad V_{\text{eff}}''(\pi) = -mg\ell + \frac{ma^2\gamma^2}{2}.$$
This is positive, so $\Theta = \pi$ is a minimum, exactly when $a^2\gamma^2 > 2g\ell$. ∎ The slow oscillation frequency about the inverted position is $\Omega_{\text{eff}}^2 = V_{\text{eff}}''(\pi)/(m\ell^2) = a^2\gamma^2/(2\ell^2) - g/\ell$. (At $\Theta = 0$, $V_{\text{eff}}'' = mg\ell + ma^2\gamma^2/2 > 0$: the bottom stays stable and the shaking stiffens it.)

(d) The reference code is `kapitza_rhs` in `project.py`. It integrates the full equation from (a) with `max_step=0.005` s, so every drive period is resolved. With $\ell = 1$ m and $a = 0.05$ m, $2g\ell = 19.62$.
- $\gamma = 150$ rad/s: $a^2\gamma^2 = 56.25 > 19.62$, so the inverted position should be stable. The prediction is $\Omega_{\text{eff}}^2 = 56.25/2 - 9.81 = 18.315\ \mathrm{s^{-2}}$, $\Omega_{\text{eff}} = 4.280$ rad/s, slow period $1.47$ s. Integrated from $\theta = \pi + 0.1$ at rest, $\max|\theta - \pi| = 0.100$ rad over 10 s: the pendulum rocks slowly about the inverted position, with a small fast jitter of relative size $a/\ell$ on top.
- $\gamma = 60$ rad/s: $a^2\gamma^2 = 9.00 < 19.62$, and $\max|\theta - \pi| = 6.19$ rad: it falls over and swings around through the bottom.

Kapitza's criterion separates the two runs.

**C3.** (a) The bob is at $\mathbf{r} = \ell(t)(\sin\theta, -\cos\theta)$, so $\dot{\mathbf{r}} = \dot\ell(\sin\theta, -\cos\theta) + \ell\dot\theta(\cos\theta, \sin\theta)$. The two directions are orthogonal unit vectors, so $|\dot{\mathbf{r}}|^2 = \dot\ell^2 + \ell^2\dot\theta^2$, $T = \frac12 m\dot\ell^2 + \frac12 m\ell^2\dot\theta^2$, and $V = -mg\ell(t)\cos\theta$. Since $\ell(t)$ is prescribed, $\frac12 m\dot\ell^2$ is a pure function of $t$ with no $\theta$ or $\dot\theta$ in it. It contributes nothing to $\partial L/\partial\theta$ or $\partial L/\partial\dot\theta$ and can be dropped (step 5 of §3.10). So
$$L = \tfrac12 m\ell(t)^2\dot\theta^2 + mg\ell(t)\cos\theta.$$
$\partial L/\partial\dot\theta = m\ell^2\dot\theta$, with time derivative $m\ell^2\ddot\theta + 2m\ell\dot\ell\dot\theta$. Also $\partial L/\partial\theta = -mg\ell\sin\theta$. Dividing by $m\ell$:
$$\ell\ddot\theta + 2\dot\ell\dot\theta + g\sin\theta = 0. ✓$$
The $2\dot\ell\dot\theta$ term says that angular momentum $m\ell^2\dot\theta$ changes only through the gravitational torque. Shortening the string spins the bob up, as a skater does by pulling in their arms.

(b) $h = \dot\theta\,\partial L/\partial\dot\theta - L = \frac12 m\ell^2\dot\theta^2 - mg\ell\cos\theta$. By (3.16), $dh/dt = -\partial L/\partial t = -(m\ell\dot\ell\dot\theta^2 + mg\dot\ell\cos\theta) \ne 0$, because $L$ depends on $t$ explicitly through $\ell(t)$. $T + V = h + \frac12 m\dot\ell^2$ is not conserved either. Physically, the agent pulling the string does real work on the bob (tension times the radial speed $\dot\ell$). The constraint is rheonomic, so d'Alembert only promises that the constraint does no *virtual* work.

(c) For small angles, $\ell\ddot\theta = -2\dot\ell\dot\theta - g\theta$. Differentiate $E = \frac12 m\ell^2\dot\theta^2 + \frac12 mg\ell\theta^2$:
$$\frac{dE}{dt} = m\ell\dot\ell\dot\theta^2 + m\ell^2\dot\theta\ddot\theta + \tfrac12 mg\dot\ell\theta^2 + mg\ell\theta\dot\theta.$$
Substitute $m\ell^2\dot\theta\ddot\theta = m\ell\dot\theta(-2\dot\ell\dot\theta - g\theta) = -2m\ell\dot\ell\dot\theta^2 - mg\ell\theta\dot\theta$:
$$\frac{dE}{dt} = -m\ell\dot\ell\dot\theta^2 + \tfrac12 mg\dot\ell\theta^2 = \frac{\dot\ell}{\ell}\left(\tfrac12 mg\ell\theta^2 - m\ell^2\dot\theta^2\right) = \frac{\dot\ell}{\ell}(U - 2K),$$
where $K = \frac12 m\ell^2\dot\theta^2$ and $U = \frac12 mg\ell\theta^2$. This is exact for the linearized equation. Over one period, $\ell$ and $\dot\ell$ are nearly constant and the motion is nearly harmonic, so $\langle K\rangle = \langle U\rangle = E/2$, and
$$\left\langle\frac{dE}{dt}\right\rangle = \frac{\dot\ell}{\ell}\left(\frac E2 - E\right) = -\frac{\dot\ell}{2\ell}E.$$
So $d\ln E/dt = -\frac12\,d\ln\ell/dt$, which gives $E\propto\ell^{-1/2}$. Since $\omega = \sqrt{g/\ell}\propto\ell^{-1/2}$ as well, $E/\omega$ is constant to leading order. This is the adiabatic invariant; CM-8 identifies it as the action $\oint p\,dq/2\pi$. At a turning point $E = \frac12 mg\ell\,\theta_{\max}^2$, so $\theta_{\max}^2 = 2E/(mg\ell)\propto\ell^{-1/2}/\ell = \ell^{-3/2}$, and $\theta_{\max}\propto\ell^{-3/4}$. The approximation needs $\ell$ to change little within one period: $\dot\ell/(\ell\omega)\ll1$.

(d) The reference code is `adiabatic_rhs` in `project.py`: `[y[1], (-2*elld*y[1] - G*np.sin(y[0]))/ell_t]` with $\ell = 1 + 0.02t$, which is the full nonlinear equation from (a). Turning points are where $\dot\theta$ changes sign. Here $\omega$ runs from $3.132$ to $2.215$ rad/s, so the slowness parameter is $\dot\ell/(\ell\omega)\le0.0064$. Over the 42 turning points in 50 s, $E/\omega$ stays constant to $0.03\%$ (peak-to-peak over mean) and $\theta_{\max}\ell^{3/4}$ to $0.02\%$. Meanwhile $E$ itself falls to $0.708$ of its starting value, against the prediction $\ell^{-1/2} = 2^{-1/2} = 0.707$ as the length doubles. The invariant holds far better than the energy, which is the point: the energy is not conserved, and the adiabatic invariant replaces it as the quantity that survives slow changes.
