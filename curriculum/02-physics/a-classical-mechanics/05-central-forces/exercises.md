# CM-5 Exercises

Conventions as in the notes: $\mu$ is the reduced mass, $\ell = |\mathbf L|$, $u = 1/r$, primes are $d/d\varphi$, $V = -k/r$ ($k > 0$) for gravity and $V = +k/r$ for Coulomb repulsion, $\theta$ is the scattering angle. Constants: $GM_\odot = 1.32712440018\times10^{20}$ m³/s², $GM_\oplus = 3.986004418\times10^{14}$ m³/s², 1 AU $= 1.495978707\times10^{11}$ m, $e^2/4\pi\varepsilon_0 = 1.43996$ MeV fm. `project.py` Parts 9–10 print every number in the solutions; do the problems first.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Reduced masses.)
(a) Hydrogen: $m_p/m_e = 1836.15267$. Compute $\mu/m_e$ and the ground-state energy, $-13.605693$ eV $\times\,\mu/m_e$.
(b) Deuterium: $m_d/m_e = 3670.48297$. Hydrogen's Balmer-α line is at 656.28 nm, and wavelengths scale as $1/\mu$. How far, and in which direction, is deuterium's line shifted?
(c) Positronium (an electron bound to a positron): what is its ground-state energy?
(d) Earth–Moon: $m_\oplus = 5.9722\times10^{24}$ kg, $m_{\rm Moon} = 7.342\times10^{22}$ kg, separation 384,400 km. Find $\mu/m_{\rm Moon}$ and the barycenter's distance from Earth's center. Is it inside the Earth (radius 6371 km)?

**A2.** (Reading a conic.) In Kepler units ($\mu = k = 1$) an orbit is $r = 2/(1 + 0.5\cos\varphi)$. Find $e$, $c$, $r_p$, $r_a$, $a$, $b$, $\ell$, $E$, $T$, $v_p$, $v_a$. Check $E$ against (3.14) and $v_p$ against vis-viva (3.18).

**A3.** (Mercury.) $a = 5.7909\times10^{10}$ m, $e = 0.20563$. Compute $v_p$ and $v_a$ from (3.17), and check $v_p/v_a = r_a/r_p$.

**A4.** (Jupiter in Kepler III.) $GM_J = 1.26686534\times10^{17}$ m³/s², and Jupiter's sidereal period is 4332.59 d. If you computed the period from Jupiter's semi-major axis using $GM_\odot$ alone, how many days too long would your answer be?

**A5.** (Total cross-sections.)
(a) Integrate the hard-sphere result (3.25) over all solid angles.
(b) Integrate Rutherford's (3.27) over $\theta_0 < \theta \le \pi$ and show the result is $\pi b(\theta_0)^2$, as (3.28) claims.
(c) What happens as $\theta_0 \to 0$? What physical fact keeps real measurements finite?
(d) Compute $\sigma(\theta > 1°)/\sigma(\theta > 90°)$.

## B. Core (do all; these are the module)

**B1.** (Two bodies in a uniform field.) Two bodies interact through $V(|\mathbf r_1 - \mathbf r_2|)$ and both sit in a uniform gravitational field $\mathbf g$.
(a) Write $L$ in terms of $\mathbf R$ and $\mathbf r$ and find $\mathbf R(t)$.
(b) Show the relative motion is exactly unaffected by $\mathbf g$.
(c) Now let the external field be non-uniform, with potential $\Phi(\mathbf x)$ per unit mass. Expand the external potential energy to second order in $\mathbf r$ about $\mathbf R$. Where does coupling between $\mathbf R$ and $\mathbf r$ first appear, and what is it called?

**B2.** (Binet backward: force from orbit.)
(a) Show that if every orbit has the form $r = c/(1 + e\cos\varphi)$, the force is $F = -(\ell^2/\mu c)/r^2$.
(b) A particle moves on the circle $r = 2R\cos\varphi$, which passes through the force center. Find $F(r)$.
(c) Show that the orbit in (b) has energy exactly zero, so the particle reaches the center with infinite speed.

**B3.** (Molniya orbit.) A satellite has perigee altitude 600 km and apogee altitude 39,700 km ($R_\oplus = 6371$ km). Find $a$, $e$, $T$ (compare with half a sidereal day, 11.967 h), and $v_p$, $v_a$ by vis-viva. Check angular momentum.

**B4.** (Stability in a Yukawa potential.) $V = -(k/r)e^{-r/\lambda}$ with $k, \lambda > 0$.
(a) Using (3.20), show that circular orbits are stable if and only if $r_0 < \tfrac{1 + \sqrt5}{2}\lambda$.
(b) Check the limit $\lambda \to \infty$.
(c) Where does a potential like this occur, and what does (a) say physically?

**B5.** (Isotropic oscillator, $n = 1$.) $F = -kr$.
(a) Solve the equations of motion in Cartesian coordinates and show every orbit is an ellipse centered on the origin.
(b) With axes along the ellipse, $x = A\cos\omega t$, $y = B\sin\omega t$. Find $E$ and $\ell$ in terms of $A$, $B$.
(c) Show that $r$ goes through two full oscillations per revolution at every amplitude, consistent with (3.21).
(d) Show $E \ge \omega\ell$ with equality only for circles. What is the Kepler analogue of this inequality?

**B6.** (Closest approach and nuclear size.)
(a) For repulsive Coulomb scattering, show $r_{\min} = (k/2E)[1 + 1/\sin(\theta/2)]$.
(b) For 5.0 MeV alpha particles on gold, evaluate $r_{\min}$ at $\theta = 60°$ and $180°$.
(c) For aluminium ($Z = 13$, $A = 27$), take the contact distance as $1.2(27^{1/3} + 4^{1/3})$ fm. Above what alpha energy do head-on collisions reach contact, so that (3.27) should fail at large angles?

**B7.** (Hohmann to Mars.) Treat Earth's and Mars's orbits as coplanar circles of radius 1 AU and 1.524 AU.
(a) Find $\Delta v_1$, $\Delta v_2$ and the transfer time.
(b) At launch, how far ahead of Earth (in heliocentric angle) must Mars be?
(c) How often do launch windows recur?

**B8.** (Code: the LRL vector of a precessing orbit.) In `project.py` Part 5 the force is $-k/r^2 - \lambda/r^3$. Compute $\mathbf A = \mathbf p\times\mathbf L - \mu k\hat{\mathbf r}$ along the integrated orbit.
(a) Show (by hand and in code) that at every perihelion $\mathbf A$ points exactly at the particle.
(b) Measure the angle $\mathbf A$ advances between successive perihelia and compare with (3.22).
(c) Is $|\mathbf A|$ constant within one orbit? Compute $d\mathbf A/dt$ and predict $|\mathbf A|/\mu k$ at perihelion and aphelion.

## C. Challenge (do at least one)

**C1.** (General-relativistic precession.) For a planet, GR modifies Binet's equation to $u'' + u = GM/h^2 + (3GM/c^2)u^2$, where $h = \ell/\mu$.
(a) Put $u = (GM/h^2)w$ and show $w'' + w = 1 + \delta w^2$ with $\delta = 3(GM)^2/(h^2c^2)$.
(b) Write $w = w_0 + \delta w_1$ with $w_0 = 1 + e\cos\varphi$. Find the equation for $w_1$ and identify the resonant term.
(c) Solve it, show $w \approx 1 + e\cos[(1 - \delta)\varphi]$, and derive (3.23).
(d) The project integrates $w'' + w = 1 + \delta w^2$ for $\delta = 10^{-3}$, $e = 0.2$. Explain the size of the discrepancy it finds from $2\pi\delta$.

**C2.** (Eccentric anomaly and Kepler's equation.) Define $\psi$ by $x = a(\cos\psi - e)$, $y = b\sin\psi$, with the force center at the origin.
(a) Show this traces the ellipse (3.13) and that $r = a(1 - e\cos\psi)$.
(b) Use the area law to derive Kepler's equation $2\pi t/T = \psi - e\sin\psi$, with $t$ measured from perihelion.
(c) Show the time averages $\langle 1/r\rangle = 1/a$ and $\langle r\rangle = a(1 + e^2/2)$.
(d) What fraction of each period does the Molniya satellite of B3 spend with $r > a$? Why is that useful?

**C3.** (Scattering by a repulsive $1/r^2$ potential.) $V = \alpha/r^2$, $\alpha > 0$.
(a) Use Binet's equation to show the orbit is $u = A\sin\gamma\varphi$ with $\gamma^2 = 1 + 2\mu\alpha/\ell^2$, and that $\theta = \pi(1 - 1/\gamma)$.
(b) Find $\theta(b)$ and show $\dfrac{d\sigma}{d\Omega} = \dfrac{\alpha\pi^2}{E}\,\dfrac{\pi - \theta}{\theta^2(2\pi - \theta)^2\sin\theta}$.
(c) Show $d\sigma/d\Omega \propto \theta^{-3}$ at small angles. Then use the impulse approximation to show that $V \propto r^{-s}$ gives $d\sigma/d\Omega \propto \theta^{-2-2/s}$ at small angles, and check $s = 1$ against Rutherford.

**C4.** (L4 and L5.) Masses $m_1, m_2$ circle their center of mass (the origin) at separation $d$ with $\Omega^2 = G(m_1 + m_2)/d^3$. In the co-rotating frame a massless particle in the orbital plane has $\Phi_{\rm eff}(\mathbf r) = -Gm_1/|\mathbf r - \mathbf r_1| - Gm_2/|\mathbf r - \mathbf r_2| - \tfrac12\Omega^2|\mathbf r|^2$.
(a) Show $\nabla\Phi_{\rm eff} = \mathbf 0$ at the two points that form equilateral triangles with $m_1$ and $m_2$.
(b) Show these are the only equilibria off the line through $m_1$ and $m_2$.
(c) Routh's criterion (stated, not derived) says L4 is linearly stable if $m_2/(m_1 + m_2) < \tfrac12(1 - \sqrt{23/27})$. Evaluate the bound and test Sun–Jupiter and Earth–Moon.
