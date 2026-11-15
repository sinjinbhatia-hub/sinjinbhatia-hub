# CM-4 Exercises

Conventions as in the notes: $L = T - V$, $p_i = \partial L/\partial\dot q_i$, the energy function is $h = \sum_i \dot q_i p_i - L$ (3.8), $\delta L$ is the first-order change (3.3) under $q_i \to q_i + \varepsilon K_i$, a *symmetry* has $\delta L = 0$ with charge $Q = \sum_i p_i K_i$ (3.4), and a *quasi-symmetry* has $\delta L = dF/dt$ with $F = F(q, t)$ and charge $Q = \sum_i p_i K_i - F$ (3.5). Italic $L$ is the Lagrangian, bold $\mathbf{L}$ angular momentum, $\kappa$ the Kepler constant, $q_e$ an electric charge. Every number you are asked for is printed or asserted by `project.py` in its final section, "Reference checks for the exercises". Do the exercise first, then run the script and compare.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Cyclic coordinates.) For each Lagrangian name the cyclic coordinates, write the conserved momenta, and say in words which physical symmetry each one is.
(a) The spherical pendulum (CM-3 Example 4.7): $L = \tfrac12 m\ell^2(\dot\theta^2 + \sin^2\theta\,\dot\phi^2) + mg\ell\cos\theta$.
(b) A particle in uniform gravity: $L = \tfrac12 m(\dot x^2 + \dot y^2 + \dot z^2) - mgz$.
(c) A particle in a helical potential, cylindrical coordinates: $L = \tfrac12 m(\dot r^2 + r^2\dot\theta^2 + \dot z^2) - V(r,\ z - a\theta)$, $a$ a fixed length. Nothing is cyclic. Find a generator $K = (K_r, K_\theta, K_z)$ for which (3.3) holds anyway, and write the charge (3.4).
(d) The Kepler problem in Cartesian coordinates, $L = \tfrac12 m|\dot{\mathbf r}|^2 + \kappa/|\mathbf r|$. Is anything cyclic? Is anything conserved?

**A2.** (The energy function, three times.) Compute $h$ from (3.8) for each Lagrangian; say whether it is conserved and whether it equals $T + V$.
(a) $L = \tfrac12 m(\dot r^2 + r^2\dot\theta^2) - V(r)$.
(b) The driven oscillator $L = \tfrac12 m\dot x^2 - \tfrac12 kx^2 + xF_0\cos\Omega t$. Also compute $dh/dt$ from (3.9).
(c) A charge $q_e$ in the field $\mathbf B = B\hat{\mathbf z}$, Landau gauge (CM-3 Example 4.8), moving in the plane: $L = \tfrac12 m(\dot x^2 + \dot y^2) + q_eBx\dot y$. What happened to the magnetic term, and what does that say about the work a magnetic field does?

**A3.** (Euler's theorem.) (a) Verify (3.10) by direct differentiation for $V = -\kappa/\sqrt{x^2 + y^2 + z^2}$ (degree $-1$), for $V = \tfrac12 k(x^2 + y^2)$ (degree $2$) and for $f = x^2y/z$ (degree $2$). (b) Show that $f(x, y) = x^2 + y$ is not homogeneous of any degree, and that (3.10) fails for it whatever $k$ you try. (c) Without using Euler's theorem, show directly that $T = \tfrac12\sum_{j,k}M_{jk}(q)\,\dot q_j\dot q_k$ with $M$ symmetric satisfies $\sum_i \dot q_i\,\partial T/\partial\dot q_i = 2T$.

**A4.** (Scaling exponents.) From Theorem 3.8, make a table of the exponents $\nu$ in $(\text{time}, \text{speed}, \text{energy}, \text{angular momentum}) \propto (\text{size})^{\nu}$ for $k = -2, -1, 1, 2, 4$. Then answer from the table alone: (a) how does a pendulum's period depend on its length (uniform gravity is $k = 1$)? (b) Two Kepler orbits of the same shape with $a' = 4a$: the ratios of periods, speeds, energies and angular momenta. Check two of them against $E = -\kappa/2a$ and $L^2 = m\kappa a(1 - e^2)$. (c) The amplitude of a quartic oscillator $V = \beta x^4$ is doubled; what happens to the period?

**A5.** (Checking $\delta L$.) $L = \tfrac12 m(\dot x^2 + \dot y^2) - V(y)$. Compute $\delta L$ from (3.3) for (a) $K = (1, 0)$, (b) $K = (0, 1)$, (c) $K = (-y, x)$, (d) $K = (y, 0)$ (a shear). Which are symmetries, which quasi-symmetries, which neither? Write the charge of every symmetry or quasi-symmetry. For (b), is there a $V$ for which it becomes a quasi-symmetry?

**A6.** (Time translation as a generator.) For $L = \tfrac12 m\dot x^2 - V(x)$ take $K = \dot x$, so $\dot K = \ddot x$. (a) Compute $\delta L$ from (3.3) and show it equals $dL/dt$, so that $F = L$. (b) Write the charge (3.5) and show it is $h$. (c) Repeat (a) with $V = V(x, t)$: what is $\delta L - dL/dt$ now, for which $V$ is it still a total derivative, and why is $h$ otherwise not conserved?

## B. Core (do all)

**B1.** (Uniform gravity: eight of the ten.) $L = \tfrac12 m(\dot x^2 + \dot y^2 + \dot z^2) - mgz$.
(a) Which of the three translations are strict symmetries? Write their charges.
(b) Show that $z \to z + \varepsilon$ is a quasi-symmetry: find $F$ and the charge $Q$. Interpret $Q$ using the explicit solution $\dot z(t) = \dot z_0 - gt$.
(c) Which rotations are symmetries? Show that the rotation about $\hat{\mathbf x}$, $K = (0, -z, y)$, is not even a quasi-symmetry: no $F(q, t)$ exists.
(d) The boost along $\hat{\mathbf z}$, $K = (0, 0, t)$. Find $F$ and $Q$, and show that "$Q$ is constant" is the same statement as $z(t) = z_0 + \dot z_0 t - \tfrac12 gt^2$.
(e) Count the Galilean generators (the ten of §1) that survive as symmetries or quasi-symmetries. Which two are lost? Then use the picture "uniform gravity is a free particle seen from a frame accelerating along $\hat{\mathbf z}$" to show that the two *conserved quantities* are not lost at all, only deformed: find a generator of the form $K = (0, -z - f(t), y)$ that is a quasi-symmetry, and its charge.

**B2.** (Two bodies, ten charges.) $L = \tfrac12 m_1|\dot{\mathbf r}_1|^2 + \tfrac12 m_2|\dot{\mathbf r}_2|^2 - V(|\mathbf r_1 - \mathbf r_2|)$, $M = m_1 + m_2$, $\mu = m_1m_2/M$.
(a) With $\mathbf R = (m_1\mathbf r_1 + m_2\mathbf r_2)/M$ and $\mathbf s = \mathbf r_1 - \mathbf r_2$, show $L = \tfrac12 M|\dot{\mathbf R}|^2 + \tfrac12\mu|\dot{\mathbf s}|^2 - V(s)$.
(b) Which coordinates are now cyclic? Which Noether charges of §3.3 are their momenta?
(c) Show $\mathbf L = \sum_a \mathbf r_a\times m_a\dot{\mathbf r}_a = M\,\mathbf R\times\dot{\mathbf R} + \mu\,\mathbf s\times\dot{\mathbf s}$, and that each term is conserved by itself. Which symmetry of which half of $L$ is responsible for each?
(d) The boost $K_{\mathbf R} = \mathbf v t$, $K_{\mathbf s} = \mathbf 0$: find $F$ and $Q$, and compare with §3.3c.
(e) Compute $h$ and show that $E_{\text{rel}} = \tfrac12\mu|\dot{\mathbf s}|^2 + V(s)$ is conserved by itself. Conclude that the relative motion is a one-body central-force problem for a particle of mass $\mu$ (the starting point of CM-5).

**B3.** (Charged particle in a magnetic field.) $L = \tfrac12 m(\dot x^2 + \dot y^2) + q_eBx\dot y$: CM-3 Example 4.8 with $\mathbf B = B\hat{\mathbf z}$, Landau gauge $\mathbf A = (0, Bx, 0)$, motion in the plane. Write $\omega_c = q_eB/m$.
(a) Write the canonical momenta $p_x, p_y$ (CM-3 eq. 4.18) and the Euler–Lagrange equations. Show the motion is uniform circular with angular frequency $\omega_c$ and radius $mv/(q_eB)$.
(b) $y$ is cyclic. Write the conserved $p_y$.
(c) $x$ is not cyclic, yet space is homogeneous in $x$. Show that $x \to x + \varepsilon$ is a quasi-symmetry: find $F$ and the charge $Q_x$. Why is $Q_x$ neither the kinetic momentum $m\dot x$ nor the canonical momentum $p_x$?
(d) Show that $X_c = p_y/(q_eB)$ and $Y_c = -Q_x/(q_eB)$ are the coordinates of the center of the circle, and that the distance from $(x, y)$ to $(X_c, Y_c)$ is the Larmor radius $mv/(q_eB)$. The two conserved momenta *are* the guiding center.
(e) The rotation $K = (-y, x)$: show it is a quasi-symmetry with $F = \tfrac12 q_eB(x^2 - y^2)$ and charge $Q_{\text{rot}} = m(x\dot y - y\dot x) + \tfrac12 q_eB(x^2 + y^2)$. Express $Q_{\text{rot}}$ through $R_c^2 = X_c^2 + Y_c^2$ and the Larmor radius $\rho_L$. Is it independent of $p_y$ and $Q_x$?
(f) In the symmetric gauge $\mathbf A = \tfrac12 B(-y, x, 0)$, $L' = \tfrac12 m(\dot x^2 + \dot y^2) + \tfrac12 q_eB(x\dot y - y\dot x)$. Which of the three transformations are now strict symmetries? Do the charges change? What does this say about "symmetry of $L$" versus "symmetry of the physics"?

**B4.** (Isotropic versus anisotropic oscillator.) $L = \tfrac12 m(\dot x^2 + \dot y^2) - \tfrac12 m(\omega_1^2x^2 + \omega_2^2y^2)$.
(a) Compute $\delta L$ for the rotation $K = (-y, x)$.
(b) Compute $dL_z/dt$ for $L_z = m(x\dot y - y\dot x)$ directly from the equations of motion and compare with (a). Which step of the proof of Theorem 3.3 is this?
(c) For $\omega_1 = \omega_2 = \omega$ the general solution is $x = A\cos\omega t$, $y = B\cos(\omega t - \phi)$. Show the orbit is an ellipse centered at the origin (not at a focus), and that $L_z = m\omega AB\sin\phi$.
(d) For $\omega_2 = 2\omega_1$: is the orbit closed? Which of $L_z$, $E_x = \tfrac12 m\dot x^2 + \tfrac12 m\omega_1^2x^2$ and $E_y$ are conserved? Show that $Q_3 = (\dot x^2 - \omega_1^2x^2)\,\dot y + 4\omega_1^2\,x\dot x\,y$ is also conserved (hint: $a_x = \dot x + i\omega_1 x$, $a_y = \dot y + i\omega_2 y$), and explain why Theorem 3.5 neither predicts nor forbids it.

**B5.** (Bead on a rotating rod: $h \ne E$ with the solution in closed form.) A bead of mass $m$ slides without friction along a horizontal rod that a motor turns at constant $\omega$ about a vertical axis through one end; $r$ is the bead's distance from the axis.
(a) Write $L$, the equation of motion, and its general solution.
(b) Compute $h$ and $E = T + V$. Which is conserved? Identify $T_2, T_1, T_0$ and check (3.11).
(c) The bead is released at $r_0$, at rest relative to the rod. Show $r(t) = r_0\cosh\omega t$ and compute $h$ and $E(t)$ explicitly. Check $E - h = 2T_0$.
(d) The bead's angular momentum about the axis is $L_z = mr^2\omega$. Compute the motor torque $\tau = dL_z/dt$ and its power $\tau\omega$, and show it equals $dE/dt$.
(e) Numbers: $m = 0.5$ kg, $\omega = 2$ rad/s, $r_0 = 0.1$ m. At $t = 1$ s give $r$, $\dot r$, $h$, $E(0)$, $E(1\ \text{s})$ and the motor power.
(f) The hoop of Example 4.3 has $E$ oscillating in a bounded range; here $E$ grows without bound. What in the two effective potentials makes the difference?

**B6.** (Mechanical similarity and the virial theorem, with numbers.) Take $m = 1$ and $V = |x|^k/k$.
(a) From energy conservation derive the period at amplitude $A$: $T(A) = 4\int_0^A dx\big/\sqrt{2\,(V(A) - V(x))}$. Substitute $x = Au$ and show $T(A) = A^{1 - k/2}\,T(1)$. This is Theorem 3.8 for this system, by hand.
(b) For $k = 4$ evaluate $T(1)$ in closed form. (Hint: $s = u^4$ turns the integral into a Beta function; $B(p, q) = \Gamma(p)\Gamma(q)/\Gamma(p + q)$, $\Gamma(\tfrac12) = \sqrt\pi$.) Compare with the project's $7.416299$ and give $T(2)$.
(c) Use the virial theorem for $k = 4$ to express $\langle T\rangle$ and $\langle V\rangle$ in terms of $E$; evaluate for $A = 1$ and compare with the project.
(d) For a bound Kepler orbit use the virial theorem and $E = -\kappa/2a$ to show that the time average of $1/r$ is exactly $1/a$. (The project also prints the time average of $r$ itself; is it $a$?)

**B7.** (The size and direction of the Laplace–Runge–Lenz vector.)
(a) Square (3.14) to prove $A^2 = m^2\kappa^2 + 2mEL^2$. You will need $|\mathbf p\times\mathbf L|^2 = p^2L^2$ and $\hat{\mathbf r}\cdot(\mathbf p\times\mathbf L) = L^2/r$; prove both.
(b) Hence $e = A/(m\kappa) = \sqrt{1 + 2EL^2/(m\kappa^2)}$. Classify $E < 0$, $E = 0$, $E > 0$ by $e$. Evaluate $e$ for the project's orbit: $m = \kappa = 1$, $E = -\tfrac12$, $L^2 = \tfrac34$.
(c) Show that $\mathbf A$ points from the force center to the pericenter: evaluate $\mathbf A$ at pericenter, where $\mathbf r\cdot\mathbf p = 0$, and fix the sign of the coefficient. What does the same computation give at apocenter?
(d) Count constants: $E, \mathbf L, \mathbf A$ are seven numbers with two relations. For a general central force only $E$ and $\mathbf L$ are conserved. What does each count say about the region of the orbital plane a bound orbit fills? Which count applies to the project's $V = -\kappa/r - \epsilon/r^3$?

**B8.** (Code.) Extend Part 3 of `project.py`: for $\epsilon \in \{0.0005, 0.001, 0.002, 0.004\}$ integrate the perturbed Kepler problem from the same initial state, measure the pericenter advance per orbit from the pericenter events, compare with the first-order prediction $6\pi\kappa\epsilon m^2/L^4$, and plot both against $\epsilon$. Then plot the relative error of the first-order formula against $\epsilon$ and fit a straight line through the origin. Report the slope and explain, from the structure of perturbation theory, why the error should be linear in $\epsilon$. (The reference implementation is in the script's final section; write yours first.)

## C. Challenge (do at least one)

**C1.** (The isotropic oscillator's hidden charges.) $L = \tfrac12 m(\dot x^2 + \dot y^2) - \tfrac12 m\omega^2(x^2 + y^2)$. Besides $E$ and $L_z = m(x\dot y - y\dot x)$ consider
$$Q_1 = \tfrac12 m(\dot x^2 - \dot y^2) + \tfrac12 m\omega^2(x^2 - y^2), \qquad Q_2 = m(\dot x\dot y + \omega^2 xy).$$
(a) Show by differentiating along solutions that $Q_1$, $Q_2$ and $L_z$ are conserved.
(b) Find generators. Show that $K = (\dot y, \dot x)$ gives $\delta L = dF/dt$ with $F = m(\dot x\dot y - \omega^2xy)$ and charge $Q_2$, and that $K = (\dot x, -\dot y)$ gives $Q_1$ (what is $F$?). Why does the project's `total_derivative_potential` refuse these cases even though the proof of Theorem 3.4 goes through?
(c) Prove $E^2 = Q_1^2 + Q_2^2 + \omega^2 L_z^2$.
(d) Count independent constants and conclude that every orbit closes. Compare with Kepler (B7d).
(e) Check everything on the solution $x = \cos t$, $y = 2\sin t$ ($m = \omega = 1$).
(f) With $a_x = \dot x + i\omega x$ and $a_y = \dot y + i\omega y$, show $\dot a_x = i\omega a_x$, $\dot a_y = i\omega a_y$, $E = \tfrac12 m(|a_x|^2 + |a_y|^2)$, $Q_1 = \tfrac12 m(|a_x|^2 - |a_y|^2)$ and $Q_2 + i\omega L_z = m\,a_x\bar a_y$. Which linear maps of $(a_x, a_y)$ preserve $E$? That group is the hidden symmetry. Say how an ordinary rotation of the plane by $\alpha$ acts on the vector $(Q_1, Q_2, \omega L_z)$.

**C2.** (The virial theorem, carefully: negative heat capacity.)
(a) A star in hydrostatic equilibrium obeys $dP/dr = -Gm(r)\rho(r)/r^2$, with $m(r)$ the mass inside radius $r$ and $P(R) = 0$ at the surface. Multiply by $4\pi r^3$, integrate from $0$ to $R$, and show $3\int P\,dV = -V_{\text{grav}}$, where $V_{\text{grav}} = -\int_0^R Gm(r)\,dm/r$.
(b) For a classical monatomic ideal gas $P = \tfrac23 u$, $u$ the thermal energy density. Conclude $2U_{\text{th}} = -V_{\text{grav}}$ exactly, for a static star, with no time average. This is hypothesis (iv) of Example 4.5(b) made honest.
(c) Verify (a) for the uniform-density sphere: show $P(r) = \tfrac{2\pi}{3}G\rho^2(R^2 - r^2)$ satisfies hydrostatic equilibrium, and compute both sides as $\tfrac35\,GM^2/R$.
(d) The total energy is $E = U_{\text{th}} + V_{\text{grav}} = -U_{\text{th}} = -\tfrac32 Nk_BT$. Show the heat capacity $dE/dT$ is negative and draw the two consequences: a star that radiates gets hotter and smaller (use $k_BT = GM\bar m/5R$ from Example 4.5(b)); a satellite in a circular orbit that loses energy to drag speeds up. Put numbers on the satellite: circular orbit $400$ km above Earth ($GM_E = 3.986\times10^{14}\ \text{m}^3\,\text{s}^{-2}$, $R_E = 6371$ km): the speed, the energy lost per kilogram per kilometer of descent, and the speed gained per kilometer. Also redo the Sun estimate of Example 4.5(b) with the realistic $|V_{\text{grav}}| \approx 1.8\,GM^2/R$.
(e) Why does negative heat capacity not contradict thermodynamics? Which ensemble forbids it, which allows it, and what happens to a self-gravitating gas put in contact with a heat bath?

**C3.** (Why $V = -c/r^2$ is special.) A particle of mass $m$ in $V = -c/r^2$, $c > 0$, with angular momentum $\ell$.
(a) $k = -2$: compute $\beta$ from (3.12) and the scaling exponent of the angular momentum. Explain why "the scaling preserves $\ell$" by itself rules out a stable circular orbit of definite radius.
(b) Write $V_{\text{eff}} = \ell^2/2mr^2 + V$ and describe the radial motion in the three cases $\ell^2 <,\ =,\ > 2mc$.
(c) Derive the orbit equation (Binet: $u = 1/r$, $u'' + u = -(m/\ell^2u^2)\,F(1/u)$ with $F = -dV/dr$) and solve it in the three cases. For $\ell^2 < 2mc$ show the particle reaches the center after infinitely many turns; for $\ell^2 > 2mc$ show it escapes to infinity at a finite angle $\theta_\infty$.
(d) Take $m = c = 1$ and start at $r = 1$ with $\dot r = 0$. Show the particle reaches radius $r$ at time $t = \sqrt{1 - r^2}/\sqrt{2 - \ell^2}$ when $\ell^2 < 2$ (so the fall time is finite) and at $t = \sqrt{r^2 - 1}/\sqrt{\ell^2 - 2}$ when $\ell^2 > 2$. For $\ell^2 = 1.9$ give the time and the angle swept when $r = 0.05$; for $\ell^2 = 2.1$ the time and angle at $r = 20$, and $\theta_\infty$. Compare with the project.
(e) Dimensional analysis: show that no length (and no time) can be built from $m$, $c$ and $\ell$ alone. Connect this to (a), and say what does set the scale.

**C4.** (A symmetry that is not a time translation: the damped oscillator.) The Caldirola–Kanai Lagrangian $L = e^{\gamma t}\big(\tfrac12 m\dot x^2 - \tfrac12 kx^2\big)$.
(a) Show the Euler–Lagrange equation is $m\ddot x + \gamma m\dot x + kx = 0$.
(b) Compute $h$ and $dh/dt$ from (3.9). Is $h$ conserved?
(c) Show that $K = \dot x + \tfrac12\gamma x$ (a time shift combined with a rescaling of $x$) gives $\delta L = dL/dt$, so $F = L$, and derive the charge $Q = e^{\gamma t}\big[\tfrac12 m\dot x^2 + \tfrac12 kx^2 + \tfrac12\gamma m\,x\dot x\big]$. Verify $dQ/dt = 0$ using (a). Explain why the time shift alone is not a quasi-symmetry but the combination is.
(d) Substitute $x = e^{-\gamma t/2}X$ and show $Q = \tfrac12 m\dot X^2 + \tfrac12 m\Omega^2X^2$ with $\Omega^2 = k/m - \gamma^2/4$: the conserved quantity is the energy of an undamped oscillator hiding inside the damped one.
(e) For $m = k = 1$, $\gamma = 0.2$, $x(0) = 1$, $\dot x(0) = 0$: what is $Q$? What range does $h$ cover for $0 \le t \le 20$ (closed form, using (d))? What is $x(20)$? Compare with the project.
