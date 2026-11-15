# CM-1 Exercises

Conventions as in the notes: SI units, $g = 9.81\ \mathrm{m/s^2}$, $y$ up, $V$ for potential energy with $\mathbf{F} = -\nabla V$, $\gamma = b/m$ and $\omega_0^2 = k/m$ for the oscillator, pendulum angle from the downward vertical.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Three particles: $m_1 = 2$ kg at $\mathbf{r}_1 = (1, 0, 0)$ m moving at $\mathbf{v}_1 = (0, 3, 0)$ m/s; $m_2 = 1$ kg at $(0, 2, 0)$ m moving at $(-1, 0, 0)$ m/s; $m_3 = 3$ kg at $(0, 0, -1)$ m moving at $(2, 1, 0)$ m/s. Compute $M$, $\mathbf{R}$, $\mathbf{P}$, and $\dot{\mathbf{R}}$. Verify $\mathbf{P} = M\dot{\mathbf{R}}$. Then compute the total angular momentum $\mathbf{L}$ about the origin.

**A2.** Decide whether each force is conservative, and find $V$ when it is:
(a) $\mathbf{F} = (y, x, 0)$; (b) $\mathbf{F} = (-y, x, 0)$; (c) $\mathbf{F} = (2xy, x^2 + z^2, 2yz)$; (d) $\mathbf{F} = -b\mathbf{v}$ (linear drag).

**A3.** A baseball ($m = 0.145$ kg) falls from rest with quadratic drag $c = 1.225\times 10^{-3}$ kg/m. Find the terminal speed $s_\infty$, the time constant $\tau = s_\infty/g$, and the time to reach 90% of terminal speed. Compare $s_\infty$ to the speed a dropped ball would have after 5 s with no drag.

**A4.** $V(x) = V_0\left(\dfrac{x}{a} + \dfrac{a}{x}\right)$ for $x > 0$, with $V_0, a > 0$. Find the equilibrium, show it is stable, and find the small-oscillation frequency using (3.8). Then confirm the result of (3.7) for a pure harmonic potential: show that $T_{\text{period}} = \sqrt{2m}\int_{-A}^{A}dx/\sqrt{E - \tfrac12 kx^2}$ with $E = \tfrac12 kA^2$ equals $2\pi\sqrt{m/k}$.

**A5.** Atwood machine with $m_1 = 0.8$ kg and $m_2 = 1.2$ kg. Compute the acceleration and the tension from (3.16). Check the limits $m_1 = m_2$ and $m_2 \gg m_1$ of (3.16) and say in one sentence why each makes sense. What is the force the pulley's axle must supply?

## B. Core (do all)

**B1.** (Splitting into CM and relative motion.) For a system of particles write $\mathbf{r}_i = \mathbf{R} + \mathbf{r}_i'$ where $\mathbf{R}$ is the center of mass. Show that $\sum_i m_i\mathbf{r}_i' = \mathbf{0}$ and hence:
(a) $\mathbf{L} = \mathbf{R}\times\mathbf{P} + \sum_i\mathbf{r}_i'\times m_i\dot{\mathbf{r}}_i'$ (angular momentum of the CM plus angular momentum about the CM);
(b) $T = \tfrac12 M\dot{\mathbf{R}}^2 + \sum_i\tfrac12 m_i\dot{\mathbf{r}}_i'^2$ (König's theorem).
Why do the cross terms vanish in both?

**B2.** (Curl-free is not enough.) Let $\mathbf{F} = \dfrac{(-y, x, 0)}{x^2 + y^2}$ on $\mathbb{R}^3$ minus the $z$-axis.
(a) Compute $\nabla\times\mathbf{F}$ and show it is $\mathbf{0}$ everywhere $\mathbf{F}$ is defined.
(b) Compute $\oint\mathbf{F}\cdot d\mathbf{r}$ counterclockwise around the unit circle in the $xy$-plane.
(c) Which hypothesis of Theorem 3.3 fails? Find a function $V$ with $\mathbf{F} = -\nabla V$ on the region $x > 0$, and explain why it cannot be extended to the whole punctured plane.

**B3.** (Double well.) $V(x) = -\tfrac12 x^2 + \tfrac14 x^4$ for a particle of mass $m = 1$.
(a) Find all equilibria and classify them (stable/unstable) via $V''$.
(b) Find the small-oscillation frequency about each stable equilibrium.
(c) What is the energy of the separatrix? Sketch the phase portrait, showing at least: one closed orbit inside each well, the separatrix (a figure-eight), and one orbit outside.
(d) A particle starts at rest at $x = 2$. Find its speed as it passes $x = 0$ and the turning point on the other side.

**B4.** (Period in a quartic well.) $V(x) = \alpha x^4$, $\alpha > 0$, amplitude $A$ (release from rest at $x = A$).
(a) Starting from (3.7), substitute $x = Au$ to show $T_{\text{period}} = \dfrac{4}{A}\sqrt{\dfrac{m}{2\alpha}}\displaystyle\int_0^1\frac{du}{\sqrt{1 - u^4}}$.
(b) The integral is $\dfrac{\Gamma(1/4)\Gamma(1/2)}{4\,\Gamma(3/4)} \approx 1.31103$ (you may quote this; the project checks it numerically). So the period is *inversely* proportional to amplitude. Explain in one sentence, from the shape of $V$, why a bigger swing is faster here, while for the pendulum it is slower.
(c) Verify the scaling numerically by extending `project.py` (see "Your turn").

**B5.** (Power at resonance.) For the steady state $x = A\cos(\omega t - \delta)$ of §3.8:
(a) Show the velocity amplitude $\omega A(\omega)$ is maximized exactly at $\omega = \omega_0$, for any damping.
(b) Show the time-averaged power delivered by the drive, $\langle F\dot{x}\rangle$, equals $\tfrac12 F_0\omega A\sin\delta$, and that this equals the time-averaged power dissipated by damping, $\langle b\dot{x}^2\rangle = \tfrac12 b\omega^2A^2$. (Use $\sin\delta = \gamma\omega/\sqrt{(\omega_0^2 - \omega^2)^2 + \gamma^2\omega^2}$.)
(c) Write $\langle P\rangle$ as a function of $\omega$ and show that, for $\gamma \ll \omega_0$, the two frequencies where it falls to half its peak are $\omega \approx \omega_0 \pm \gamma/2$, so the full width at half maximum is $\gamma$. (Approximate $\omega_0^2 - \omega^2 \approx 2\omega_0(\omega_0 - \omega)$ near the peak.) Compare with the project's measured FWHM for Example 4.4.

**B6.** (Rocket with gravity.) A rocket burns fuel at a constant rate $|\dot{m}| = \beta$ with exhaust speed $u$, starting from rest at mass $m_0$ in uniform gravity $g$.
(a) Show $v(t) = u\ln\dfrac{m_0}{m_0 - \beta t} - gt$.
(b) Show the rocket leaves the pad only if $u\beta > m_0 g$. What happens physically if this fails?
(c) For $m_0 = 1000$ kg, $u = 2500$ m/s, $\beta = 11.667$ kg/s, burn time 60 s, compute $v(60)$ and the height reached at burnout (integrate $v$; you will need $\int\ln(m_0 - \beta t)\,dt$). Compare $v(60)$ with Example 4.6.

**B7.** (Code: Coulomb friction.) A mass on a spring on a rough horizontal table feels friction of constant magnitude $\mu mg$ opposing the velocity: $m\ddot{x} = -kx - \mu mg\,\operatorname{sgn}(\dot{x})$.
(a) Over one half-cycle, from rest at $x = A_1 > 0$ to rest at $x = -A_2$, use the work–energy theorem (3.4) to show $A_1 - A_2 = 2\mu mg/k$. Conclude the amplitude decreases by $4\mu mg/k$ per full period, linearly, until the mass stops. (Contrast: linear damping gives exponential decay.)
(b) Add this system to `project.py` under "Your turn", integrate with RK4 using $m = 1$, $k = 100$, $\mu = 0.1$, $x_0 = 0.5$, and extract the successive turning-point amplitudes. Assert they decrease by $4\mu mg/k = 0.03924$ per period to within 1%. When does the mass stop for good? (Hint: it sticks when $k|x| \le \mu mg$ at a turning point.)

**B8.** (General bead on a wire.) A bead slides on a frictionless wire $y = f(x)$ under gravity. Following §3.10(b) with the normal direction $\propto(-f'(x), 1)$:
(a) Derive $\ddot{x} = -\dfrac{f'\left(g + f''\dot{x}^2\right)}{1 + f'^2}$ and $N = m\left(g + f''\dot{x}^2 + f'\ddot{x}\right)\sqrt{1 + f'^2}$.
(b) Check that $f = kx^2$ reproduces (3.14) and (3.15).
(c) Show that energy $E = \tfrac12 m(1 + f'^2)\dot{x}^2 + mgf(x)$ is conserved by (a) directly: differentiate $E$ and substitute $\ddot{x}$. (This is the calculation CM-3 will make unnecessary.)

## C. Challenge (do at least one)

**C1.** (Next term in the pendulum period.) Expand $(1 - k^2\sin^2\phi)^{-1/2}$ to order $k^4$ and use $\int_0^{\pi/2}\sin^4\phi\,d\phi = 3\pi/16$ to show $K(k) = \tfrac{\pi}{2}\left(1 + \tfrac14 k^2 + \tfrac{9}{64}k^4 + O(k^6)\right)$. Then expand $k^2 = \sin^2(\theta_0/2)$ to order $\theta_0^4$ and show
$$T_{\text{period}} = T_0\left(1 + \frac{\theta_0^2}{16} + \frac{11\,\theta_0^4}{3072} + O(\theta_0^6)\right).$$
Evaluate at $\theta_0 = 90°$ and compare with the exact $1.18034$ and the one-term value $1.15421$. The project's SymPy cell reproduces the series; make sure your hand calculation agrees with it.

**C2.** (Why symplectic Euler works.) For the system $\dot{x} = v$, $\dot{v} = a(x)$ with step $h$:
(a) Compute the Jacobian matrix of the explicit Euler map $(x, v)\mapsto(x + hv,\ v + ha(x))$ and show its determinant is $1 - h^2a'(x)$. For the harmonic oscillator this is $1 + h^2\omega^2$, recovering (3.17).
(b) Compute the Jacobian of the symplectic Euler map $(x, v)\mapsto(x + h(v + ha(x)),\ v + ha(x))$ and show its determinant is exactly $1$ for every $a(x)$.
(c) For the harmonic oscillator ($a = -\omega^2x$), show that symplectic Euler exactly conserves the *modified* energy $\tilde{E} = \tfrac12 v^2 + \tfrac12\omega^2x^2 - \tfrac12 h\omega^2xv$. Conclude that the true energy $E = \tfrac12 v^2 + \tfrac12\omega^2x^2$ stays within a bounded band of width $O(h)$ forever. (The project checks this identity with SymPy and numerically.)

**C3.** (Pendulum with a shaken pivot, Newtonian style.) The pivot of a simple pendulum moves horizontally as $X(t) = A\cos\Omega t$. The bob is at $(X + \ell\sin\theta, -\ell\cos\theta)$.
(a) Write Newton's equations for the bob in Cartesian components with the unknown tension $T$, as in §3.10(a).
(b) Eliminate $T$ to obtain $\ell\ddot{\theta} = -g\sin\theta + A\Omega^2\cos(\Omega t)\cos\theta$, and find $T$ in terms of $\theta, \dot{\theta}, t$.
(c) For small $\theta$ this is a driven oscillator. Identify $\omega_0$ and the effective drive, and write down the steady-state amplitude of $\theta$ from (3.10) with $\gamma = 0$. What happens as $\Omega \to \omega_0$, and what does the model leave out that would prevent it?
(d) One sentence: what will CM-3 need to write down to get (b)?
