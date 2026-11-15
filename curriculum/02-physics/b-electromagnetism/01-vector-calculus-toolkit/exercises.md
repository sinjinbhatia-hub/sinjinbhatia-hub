# EM-1 Exercises

Conventions as in the notes: index notation with summed repeated indices, $\varepsilon_{123} = +1$; spherical $(r, \theta, \phi)$ with $\theta$ from $+z$; cylindrical $(s, \phi, z)$; $d\mathbf a$ outward on closed surfaces and right-hand-rule oriented on open ones; $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Let $f = x^2y + yz^3$ and $\mathbf F = (xy,\ yz,\ zx)$.
(a) Compute $\nabla f$, $\nabla^2f$, $\nabla\cdot\mathbf F$, $\nabla\times\mathbf F$, and evaluate each at $(1, 1, 1)$.
(b) Verify $\nabla\times(\nabla f) = \mathbf 0$ and $\nabla\cdot(\nabla\times\mathbf F) = 0$ by direct computation.

**A2.** Let $\mathbf A = (1, 2, 0)$, $\mathbf B = (0, 1, 3)$, $\mathbf C = (2, 0, 1)$.
(a) Compute $\mathbf A\cdot(\mathbf B\times\mathbf C)$. What does it measure? Check that $\mathbf B\cdot(\mathbf C\times\mathbf A)$ gives the same number.
(b) Compute $\mathbf A\times(\mathbf B\times\mathbf C)$ directly, then again with BAC-CAB (3.2).
(c) Compute $(\mathbf A\times\mathbf B)\times\mathbf C$ and compare with (b).

**A3.** The temperature in a room is $T = x^2 + 2y^2 + 3z^2 + xy$ (°C, positions in metres). At $(1, -1, 2)$ find: $\nabla T$; the unit vector along which $T$ rises fastest and that maximum rate; the rate of change along $(1, 1, 1)/\sqrt3$.

**A4.** (a) Use (3.28) to compute $\nabla\cdot(r^n\hat{\mathbf r})$ for $r \ne 0$.
(b) Use the cylindrical divergence to compute $\nabla\cdot(s\,\hat{\mathbf s})$, and check against the Cartesian form of the same field, $(x, y, 0)$.
(c) Use (3.29) to compute $\nabla^2(r^n)$. For which $n$ is $r^n$ harmonic ($\nabla^2 = 0$) away from the origin?

**A5.** Evaluate:
(a) $\int_{\text{all space}}(x^2 + y^2 + z^2)\,\delta^3(\mathbf r - \mathbf a)\,d\tau$ with $\mathbf a = (1, 2, 1)$.
(b) $\int_V(x + y + z)\,\delta^3(\mathbf r - 2\hat{\mathbf x})\,d\tau$, where $V$ is the ball of radius 1 centred at the origin.
(c) $\int_{\text{all space}}(r^2 + 1)\,\nabla^2\!\left(\dfrac{1}{|\mathbf r - \hat{\mathbf z}|}\right)d\tau$.

**A6.** Let $\mathbf F = (y, -x, z)$.
(a) Compute $\int\mathbf F\cdot d\mathbf l$ along the helix $\mathbf r(t) = (\cos t, \sin t, t/2\pi)$, $0 \le t \le 2\pi$.
(b) Compute it along the straight segment from $(1, 0, 0)$ to $(1, 0, 1)$ (same endpoints).
(c) Compute $\nabla\times\mathbf F$ and use Stokes' theorem to explain the difference between (a) and (b).

## B. Core (do all)

**B1.** Prove $\nabla\times(f\nabla g) = \nabla f\times\nabla g$ twice: once from (P5) and Theorem 3.5, once directly in index notation. Deduce that $\oint_P f\,\nabla g\cdot d\mathbf l = \int_S(\nabla f\times\nabla g)\cdot d\mathbf a$, and that $\oint_P f\,\nabla f\cdot d\mathbf l = 0$ for every loop.

**B2.** (The vortex.) Let $\mathbf F = \hat{\boldsymbol\phi}/s$ on $\mathbb R^3$ minus the $z$-axis; in Cartesian form this is $(-y, x, 0)/(x^2 + y^2)$.
(a) Using the cylindrical curl from the §3.10 table, compute $\nabla\times(s\,\hat{\boldsymbol\phi})$ (the whirlpool) and $\nabla\times(\hat{\boldsymbol\phi}/s)$.
(b) Compute $\oint\mathbf F\cdot d\mathbf l$ counterclockwise around the circle of radius $a$ about the $z$-axis in the plane $z = 0$, and around the unit circle centred at $(2, 0, 0)$ in the same plane.
(c) Explain why (b) does not contradict Stokes' theorem, and argue that in the sense of distributions $\nabla\times\mathbf F = 2\pi\,\delta(x)\delta(y)\,\hat{\mathbf z}$.
(d) Which hypothesis of Theorem 3.16 fails for this field?

**B3.** (Divergence theorem on a half-ball.) Let $\mathbf F = (xz,\ yz,\ z^2 + 1)$ and $V$ the half-ball $r \le R$, $z \ge 0$. Compute $\int_V\nabla\cdot\mathbf F\,d\tau$, the flux through the curved hemisphere, and the flux through the flat base, all with outward normals, and verify (3.18). Give the three numbers for $R = 2$.

**B4.** (Stokes on a triangle.) Let $P$ be the triangle $(1,0,0)\to(0,1,0)\to(0,0,1)\to(1,0,0)$ and $\mathbf F = (y, z, x)$.
(a) Compute $\oint_P\mathbf F\cdot d\mathbf l$ edge by edge.
(b) Compute $\nabla\times\mathbf F$ and its flux through the flat triangle, with the normal the right-hand rule assigns to $P$. Which way does that normal point?
(c) What changes if $P$ is traversed in the opposite order?

**B5.** Prove (P5) and (P6) in index notation. Check (P6) on $\mathbf A = \boldsymbol\omega$ (constant), $\mathbf B = \mathbf r$ against Example 4.1.

**B6.** (Where the divergence theorem breaks.) Let $\mathbf F_n = r^n\hat{\mathbf r}$ and $V$ the ball of radius $R$ about the origin.
(a) For $n = -1$, compute both sides of (3.18).
(b) For general $n > -2$, compute both sides and show they agree.
(c) For $n = -2$, compute both sides. Which theorem of the notes repairs the mismatch, and what is $\nabla\cdot\mathbf F_{-2}$ as a distribution?
(d) For $n = -3$, compute the flux through the sphere $r = R$ and through a small sphere $r = \epsilon$. Why can no multiple of $\delta^3(\mathbf r)$ repair (3.18) here?

**B7.** (Building a vector potential.) Let $\mathbf F = (y^2, z^2, x^2)$.
(a) Show $\nabla\cdot\mathbf F = 0$.
(b) Construct $\mathbf W$ with $\nabla\times\mathbf W = \mathbf F$ from (3.36), and check the curl.
(c) Compute $\nabla\cdot\mathbf W$. Find a $\lambda$ such that $\mathbf W' = \mathbf W + \nabla\lambda$ is in Coulomb gauge ($\nabla\cdot\mathbf W' = 0$), and verify $\nabla\times\mathbf W' = \mathbf F$.
(d) Is $\mathbf W'$ the only Coulomb-gauge potential for $\mathbf F$? If not, describe all of them.

**B8.** (Code: extend `project.py`.)
(a) Using the grid of Part 3 (`fd_div_curl`), take $f = e^{-(x^2+y^2+z^2)}$ and $\mathbf A = (xy,\ yz,\ zx)$ on $[-1, 1]^3$. Compute $\nabla\cdot(f\mathbf A)$ by finite-differencing the product, and $f\,\nabla\cdot\mathbf A + \mathbf A\cdot\nabla f$ from separately differenced pieces. Report the largest difference for $n = 21$ and $n = 41$ grid points per side and explain the ratio. Why would $\mathbf A = (y,\ z^2,\ \sin x)$ have been a useless test?
(b) Add a fourth path of your own to Example 4.6's `paths` and confirm the work is still 17. Then integrate the vortex of B2 from $(1, 0, 0)$ to $(-1, 0, 0)$ along the upper and the lower unit half-circles. Explain the two answers using B2 and Theorem 3.16.

## C. Challenge (do at least one)

**C1.** (Solid angle.) The solid angle subtended at the origin by a surface $S$ (crossed at most once by each ray from the origin) is $\Omega = \int_S\hat{\mathbf r}\cdot d\mathbf a/r^2$.
(a) Explain why this is the area of the radial projection of $S$ onto the unit sphere, so that a closed surface around the origin has $\Omega = 4\pi$.
(b) For the disk of radius $a$ in the plane $z = z_0 > 0$, centred on the $z$-axis, show $\Omega = 2\pi\left(1 - z_0/\sqrt{z_0^2 + a^2}\right)$. Evaluate for $a = z_0 = 1$, and take $a\to\infty$.
(c) (Code.) Write `square_solid_angle(z0)`, the flux of $\hat{\mathbf r}/r^2$ through the square $[-\tfrac12, \tfrac12]^2$ in the plane $z = z_0$, by 2D Gauss–Legendre quadrature. Predict the value at $z_0 = \tfrac12$ from symmetry before running it.

**C2.** (Uniqueness in Helmholtz.) Let $\mathbf G$ be $C^1$ on $\mathbb R^3$ with $\nabla\cdot\mathbf G = 0$, $\nabla\times\mathbf G = \mathbf 0$, and $|\mathbf G| \le K/r^2$ for $r \ge 1$.
(a) Show $\mathbf G = -\nabla U$ for a harmonic $U$ with $|U| \le K/r$ for $r \ge 1$.
(b) Apply Green's first identity (3.22) on a ball of radius $R$ and let $R\to\infty$ to show $\mathbf G = \mathbf 0$.
(c) Give a nonzero field with zero divergence and curl everywhere. Which hypothesis does it violate, and where does the same field appear in the FFT decomposition of §3.12?

**C3.** (Screened Coulomb potential.) Let $u = e^{-\mu r}/r$ with $\mu > 0$.
(a) Show $\nabla^2u = \mu^2u$ for $r \ne 0$.
(b) Show that $(\nabla^2 - \mu^2)u = -4\pi\delta^3(\mathbf r)$ in the sense of distributions: for every smooth $\psi$ vanishing outside a bounded region, $\int u\,(\nabla^2 - \mu^2)\psi\,d\tau = -4\pi\psi(\mathbf 0)$. (Excise a ball and use Green's second identity.)
(c) The project evaluates the left side of (b) for $\psi = e^{-r^2}$, $\mu = 1.5$. Set up that integral as a single integral in $r$ and say what it should equal.

**C4.** (Why the curl is a pseudovector.) Let $R$ be orthogonal, $\mathbf x' = R\mathbf x$, and let a vector field transform as $\mathbf F'(\mathbf x') = R\,\mathbf F(\mathbf x)$.
(a) Prove $\varepsilon_{lmn}R_{li}R_{mj}R_{nk} = (\det R)\,\varepsilon_{ijk}$.
(b) Show $\nabla'\cdot\mathbf F' = \nabla\cdot\mathbf F$ and $(\nabla'\times\mathbf F')_i = (\det R)\,R_{il}(\nabla\times\mathbf F)_l$.
(c) Check (b) by hand for the whirlpool $(-y, x, 0)$ under the reflection $z\to-z$, and explain the sign physically.
