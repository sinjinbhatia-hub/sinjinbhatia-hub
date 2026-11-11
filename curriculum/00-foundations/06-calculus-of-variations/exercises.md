# F6 Exercises

Conventions as in the notes: $J[y] = \int_a^b L(x, y, y')\,dx$; $L_y$, $L_{y'}$ are partial derivatives of $L$ with its three slots treated as independent; $\frac{d}{dx}$ is the total derivative along the curve; $\delta J = \frac{d}{d\varepsilon}J[y + \varepsilon\eta]|_{\varepsilon = 0}$ with $\eta(a) = \eta(b) = 0$. In the brachistochrone, $y$ is measured downward and $g = 9.81\ \mathrm{m/s^2}$.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Write the Euler–Lagrange equation (3.4) for each Lagrangian and find the general solution where asked.
(a) $L = y'^2 + y^2$ (solve).
(b) $L = x\,y'^2$ (solve; note which special case of §3.6 applies).
(c) $L = y'^2 - 2y\sin x$ (solve).
(d) $L = x^2 y'^2 + y^2$ (write the equation only; it is an Euler-type equation, and $y = x^r$ is a good guess).

**A2.** Let $J[y] = \int_0^1 (y'^2 + 2xy)\,dx$, $y(x) = x^2$, $\eta(x) = x(1 - x)$.
(a) Compute $\Phi(\varepsilon) = J[y + \varepsilon\eta]$ explicitly as a polynomial in $\varepsilon$, and read off $\delta J = \Phi'(0)$.
(b) Compute $\delta J$ again from formula (3.3) and confirm the two agree.
(c) Write the Euler–Lagrange equation for this $L$. Is $y = x^2$ a solution? Should it be, given (a)?

**A3.** For each Lagrangian, say which special case of §3.6 applies (no $y$, no $y'$, no explicit $x$, or none), and write the corresponding first integral or the full E–L equation.
(a) $L = \sqrt{1 + y'^2}\,/\,x$. (b) $L = y^2(1 + y'^2)$. (c) $L = y'^2 + xy$. (d) $L = e^{y}\sqrt{1 + y'^2}$.

**A4.** (a) Show that the Euler–Lagrange equation of $L = y\sqrt{1 + y'^2}$ (the catenary with $\lambda = 0$, and the catenoid) simplifies to $y\,y'' = 1 + y'^2$. (b) Verify by direct differentiation that $y = c\cosh((x - x_0)/c)$ satisfies it.

**A5.** Minimize $J[y] = \int_0^1\big(\tfrac12 y'^2 + y\big)\,dx$ with $y(0) = 0$ and the right endpoint **free**.
(a) Write the E–L equation and its general solution. (b) Impose the natural boundary condition (3.15) at $x = 1$ and find $y$. (c) Compute $J$ at your solution. Then solve the same problem with $y(1) = 0$ imposed instead, compute $J$ again, and explain in one sentence why the free-end value must be the smaller of the two.

## B. Core (do all)

**B1.** (Fermat's principle in a graded medium.) Light travels in the upper half-plane $y > 0$ with refractive index $n(y)$, so the ray $y(x)$ makes $\int n(y)\sqrt{1 + y'^2}\,dx$ stationary.
(a) For $n(y) = 1/y$, use the Beltrami identity to show that every ray is a semicircle centered on the $x$-axis, $(x - x_0)^2 + y^2 = R^2$. (These are the geodesics of the hyperbolic plane.) What are the rays in the limit $R \to \infty$?
(b) For $n(y) = n_0(1 + \kappa y)$ (air above a hot road, with $\kappa > 0$ meaning the index increases with height), show the rays are catenaries $1 + \kappa y = (C/n_0)\cosh\big(n_0\kappa(x - x_0)/C\big)$. Which way do they curve, and why does that produce a mirage?

**B2.** (Geodesics.)
(a) On a cylinder of radius $R$ with coordinates $(\varphi, z)$, $ds^2 = R^2 d\varphi^2 + dz^2$. Describe a path by $z(\varphi)$, write the arc-length Lagrangian, and show that the geodesics are helices $z = m\varphi + z_0$. Explain the result in one sentence by unrolling the cylinder.
(b) On the unit sphere, parametrize a path by $(\theta(t), \varphi(t))$ with $L = \sqrt{\dot\theta^2 + \sin^2\theta\,\dot\varphi^2}$. Write the two E–L equations (3.13). Show that $L$ contains no $\varphi$ and write the resulting first integral. Show that the equator $\theta = \pi/2$, $\varphi = \omega t$ satisfies both E–L equations.

**B3.** (The du Bois-Reymond lemma, and E–L without assuming $y''$ exists.)
(a) Let $f$ be continuous on $[a, b]$ and suppose $\int_a^b f(x)\eta'(x)\,dx = 0$ for every $\eta \in C^1[a, b]$ with $\eta(a) = \eta(b) = 0$. Prove that $f$ is constant. (Hint: let $\bar f$ be the mean of $f$ and try $\eta(x) = \int_a^x (f(t) - \bar f)\,dt$.)
(b) Use (a) to show that if $y \in C^1$ (not necessarily $C^2$) makes $J$ stationary, then $L_{y'}(x, y, y') = \int_a^x L_y(t, y, y')\,dt + \text{const}$. (Hint: in (3.3), integrate the $L_y\eta$ term by parts instead of the $L_{y'}\eta'$ term.) Explain why this is (3.4) in integrated form and what it buys you.

**B4.** (The discrete derivation with the trapezoid rule.) Replace (3.7) by the trapezoid-rule discretization
$$
J_h = \frac h2\sum_{i=0}^{N}\Big[L(x_i, y_i, p_i) + L(x_{i+1}, y_{i+1}, p_i)\Big], \qquad p_i = \frac{y_{i+1} - y_i}{h}.
$$
(a) Compute $\partial J_h/\partial y_i$ for an interior node (there are now four terms containing $y_i$) and write it in the form $h\,[\,\cdots\,]$ with a difference quotient inside the bracket.
(b) Explain why the bracket is now *centered* at $x_i$ and argue that $\partial J_h/\partial y_i = h\,\mathrm{EL}(x_i) + O(h^3)$, one order better than (3.9).
(c) (Code.) Implement this discretization in `project.py` Part 8 and verify the rate: the error should drop by a factor of about $8$ when $h$ is halved.

**B5.** (Harmonic oscillator by Beltrami.) Let $L = \tfrac12\dot y^2 - \tfrac12\omega^2 y^2$ with independent variable $t$.
(a) Write the E–L equation. (b) Write the Beltrami quantity and identify it physically. (c) Starting from the first integral of (b), separate variables and derive $y(t) = A\sin(\omega t + \phi)$ with $A$ expressed through the Beltrami constant. (d) Compute the second variation $\delta^2 S$ for the test function $\eta = \sin(\pi t/T)$ on $[0, T]$ and show that the extremal is *not* a minimum of the action when $T > \pi/\omega$. (So "least action" is literally false for an oscillator followed for more than half a period; it is stationary action.)

**B6.** (An isoperimetric problem with a computable multiplier.) Minimize $J[y] = \int_0^1 y'^2\,dx$ subject to $\int_0^1 y\,dx = 1$ and $y(0) = y(1) = 0$.
(a) Check the hypothesis of Theorem 3.7 (that $y$ is not an extremal of the constraint functional). (b) Solve the E–L equation of $L + \lambda G$ with the boundary conditions, then fix $\lambda$ from the constraint. (c) Compute the minimum value of $J$. (d) Explain why this stationary point is the global minimum.

**B7.** (Brachistochrone to a wall.) A bead starts at rest at $(0, 0)$ and must reach the vertical wall $x = 1$ (any depth) in the least time; $y$ is measured downward.
(a) Write the natural boundary condition (3.15) at $x = 1$ for $L = \sqrt{1 + y'^2}/\sqrt{2gy}$ and show it says the wire must arrive at the wall horizontally. (b) Using the cycloid (4.1), find $\theta_1$, $a$, the arrival depth, and the time $T$, with numbers. (c) Compare $T$ with the time to the fixed point $(1, 1)$ from Example 4.2 and explain the direction of the inequality.

**B8.** (Code: the catenoid.) In `project.py`, after Part 4, implement the minimal surface of revolution between the rings $y(\pm d) = 1$: minimize the polygon area $2\pi\sum_i \Delta s_i\,(y_i + y_{i+1})/2$ with no constraint. (a) For $d = 0.5$, solve $c\cosh(d/c) = 1$ numerically (two roots), compute the area $2\pi c\,[\,d + \tfrac c2\sinh(2d/c)\,]$ of each catenoid, and check which one the minimizer finds. Compare both with the area $2\pi$ of two flat disks. (b) For $d = 0.8 > 0.6627$, describe what the minimizer does and why no catenoid can exist.

## C. Challenge (do at least one)

**C1.** (Global minimality by convexity.) Let $f(p) = \sqrt{1 + p^2}$.
(a) Show $f$ is convex and hence $f(p + q) \ge f(p) + f'(p)\,q$ for all real $p, q$.
(b) Let $y$ be the straight line through $(a, A)$ and $(b, B)$ and let $\tilde y$ be any $C^1$ function with the same endpoints. Write $\tilde y = y + \eta$ and use (a) pointwise under the integral to prove $J[\tilde y] \ge J[y]$, with equality only if $\tilde y = y$. This proves that the line is the *global* minimizer, not only a stationary point.
(c) Generalize: if $L = L(x, y')$ does not contain $y$ and is convex in $y'$ for each $x$, every extremal with the given endpoints is a global minimizer.
(d) Explain why this argument does not apply directly to the brachistochrone Lagrangian $\sqrt{1 + y'^2}/\sqrt{2gy}$.

**C2.** (The brachistochrone that comes back up.)
(a) The bead starts at rest at $(0, 0)$ and must reach $(X, 0)$, a point at the *same height* a distance $X$ away. Show that the fastest wire is a full cycloid arch with $\theta_1 = 2\pi$, find $a$, the maximum depth, and the time $T$, and evaluate for $X = 1\ \mathrm{m}$. Why does a horizontal wire never get there?
(b) Prove that equation (4.2) has exactly one solution $\theta_1 \in (0, 2\pi)$ for every ratio $x_1/y_1 > 0$, by showing $f(\theta) = (\theta - \sin\theta)/(1 - \cos\theta)$ is strictly increasing on $(0, 2\pi)$ with $f \to 0$ at $0$ and $f \to \infty$ at $2\pi$.

**C3.** (Fields: Laplace's equation as an Euler–Lagrange equation.) Let $u(x, y)$ be defined on a bounded region $\Omega$ with $u$ prescribed on the boundary $\partial\Omega$, and let
$$
J[u] = \iint_\Omega L(x, y, u, u_x, u_y)\,dx\,dy .
$$
(a) Derive the Euler–Lagrange equation $\dfrac{\partial L}{\partial u} - \dfrac{\partial}{\partial x}\dfrac{\partial L}{\partial u_x} - \dfrac{\partial}{\partial y}\dfrac{\partial L}{\partial u_y} = 0$, using the divergence theorem in place of integration by parts and stating the two-dimensional version of the fundamental lemma you need.
(b) Apply it to the Dirichlet energy $L = \tfrac12(u_x^2 + u_y^2)$ and to $L = \tfrac12(u_x^2 + u_y^2) - \rho(x, y)\,u$. Name the two equations you get and say where they appear in EM-2.

**C4.** (Maximum entropy gives the Gaussian: a variational problem over densities.) Among probability densities $q$ on $\mathbb{R}$ with mean $\mu$ and variance $\sigma^2$, find the one maximizing the entropy $H[q] = -\int q\ln q\,dx$.
(a) Set up $L^* = -q\ln q + \lambda_0 q + \lambda_1 xq + \lambda_2 x^2 q$ with three multipliers for the constraints $\int q = 1$, $\int xq = \mu$, $\int x^2 q = \sigma^2 + \mu^2$, and explain why the E–L equation is algebraic (§3.6, Case B).
(b) Solve it for $q$, determine the sign $\lambda_2$ must have, and find all three multipliers in terms of $\mu$ and $\sigma$.
(c) Check Legendre's condition in the form appropriate here ($\partial^2 L^*/\partial q^2$) and say whether you have a maximum or a minimum.
(d) In one paragraph, connect this to the maximum-entropy derivation of the Boltzmann distribution in SM-4 and to the optimal variational distribution in ML-10.
