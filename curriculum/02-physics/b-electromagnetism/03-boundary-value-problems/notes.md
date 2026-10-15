# EM-3: Boundary-Value Problems: Separation of Variables, Images, Multipoles, Legendre Polynomials

**Time:** ~14 hours · **Prerequisites:** EM-2, F2 (Fourier series) · **Unlocks:** EM-4; QM-5 and QM-6 reuse the Legendre equation and the spherical harmonics

> **Conventions in this module.** SI units; $k \equiv 1/4\pi\varepsilon_0 = 8.988\times10^9$ N m²/C². Field point $\mathbf r$, source point $\mathbf r'$, separation $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$ (EM-2). Spherical coordinates $(r, \theta, \phi)$ with $\theta$ the polar angle **measured from the symmetry axis of the problem**, which is called $z$ even when it is not vertical; $x \equiv \cos\theta \in [-1, 1]$. Cylindrical coordinates $(s, \phi, z)$, $s$ the distance from the axis. Legendre polynomials are normalized by $P_\ell(1) = 1$ (not unit norm). Multipole moments follow Griffiths: monopole $Q = \int\rho\,d\tau'$, dipole $\mathbf p = \int\mathbf r'\rho\,d\tau'$, quadrupole $Q_{ij} = \int(3x'_ix'_j - r'^2\delta_{ij})\rho\,d\tau'$ (symmetric, traceless). Spherical harmonics $Y_\ell^m$ use Griffiths' table (Condon–Shortley phase, same as SciPy). At a conductor's surface $\hat{\mathbf n}$ points out of the metal into the field region and $\sigma = -\varepsilon_0\,\partial V/\partial n$ (EM-2 Theorem 3.4(v)). A "harmonic" function is one with $\nabla^2V = 0$. "Grounded" means $V = 0$.

## 0. Where this sits

You need from EM-2: Poisson's equation $\nabla^2V = -\rho/\varepsilon_0$ (3.9), the boundary conditions at a surface charge (3.10)–(3.12), the five properties of a conductor (Theorem 3.4), the field energy and surface pressure (3.16)–(3.17), and the dipole field (4.3). From F2: separation of variables for the heat equation (§3.12), the sine-series coefficient formula, and the fact that orthogonality of eigenfunctions is the spectral theorem for a self-adjoint operator (Theorem 5.2 and the Sturm–Liouville paragraph of §5). From EM-1: Green's first and second identities (Theorem 3.11).

This module answers the question EM-2 left open: given conductors at known potentials and charges you cannot see, find $V$. The three methods (images, separation of variables, multipoles) are the working tools of EM-4 and EM-5; the Legendre equation and the spherical harmonics are the angular half of hydrogen (QM-5, QM-6); the relaxation solvers are LA-13's iterative methods on a real problem.

## 1. The problem that forced this idea

Laplace's equation $\nabla^2V = 0$ was written down for gravity. Laplace, working in the 1780s on the shape of the planets and the rings of Saturn, showed that the gravitational potential of any body obeys it outside the body. Adrien-Marie Legendre, attacking the same problem in 1782 (the attraction of a spheroid), expanded the inverse distance $1/|\mathbf r - \mathbf r'|$ in powers of $r'/r$ and found that the coefficients were polynomials in the cosine of the angle between $\mathbf r$ and $\mathbf r'$. Those are the polynomials of this module, and that expansion is the multipole expansion of §3.9. The same equation governs steady heat flow, an incompressible irrotational fluid, a stretched membrane, and the electrostatic potential in empty space. Whatever is learned about one is learned about all.

Electrostatics added a new kind of difficulty. Coulomb's law tells you $V$ when you know every charge. Put a charge near a metal plate and you do not know the charges: they move until the plate is an equipotential, and where they end up is the answer, not the input. Poisson's 1812 computation of the charge on two touching conducting spheres was a heroic special case. The general statement of the problem is: find $V$ satisfying Poisson's equation in a region, taking prescribed values on the boundary (Dirichlet) or having prescribed normal derivative there (Neumann). Riemann later named the first form Dirichlet's problem, after lectures Dirichlet gave in the 1850s.

Three people made it solvable. George Green's 1828 essay introduced the "potential function," the identities that bear his name, and the idea of a function (the Green's function) that solves the problem for a point source with the boundary conditions built in. William Thomson (later Kelvin), who rescued Green's essay from obscurity in 1845, in the same year invented the method of electrical images: a point charge near a grounded sphere produces the same exterior field as the charge plus a fictitious charge inside the sphere. And Fourier's separation of variables, built for heat in 1807 and completed in F2, applies unchanged to Laplace's equation whenever the boundaries are coordinate surfaces.

One more historical fact matters here. Samuel Earnshaw proved in 1842 that no arrangement of static charges can hold a test charge in stable equilibrium. That theorem is why a classical atom cannot be a static cluster of electrons and protons, and so why the ion traps that hold single ions today use oscillating fields (Paul) or magnetic fields (Penning). It is a two-line consequence of the first property derived below.

## 2. Intuition

**Harmonic functions are as boring as the boundary allows.** Picture $V$ as the height of a rubber sheet stretched over a wire frame (two dimensions), or as the steady temperature in a block whose surface temperatures are fixed. There are no bumps in the interior: a bump would mean a point hotter than all its neighbours, and heat would flow away from it. The height at any point is the *average* of the heights on any sphere drawn around it. That is the mean-value property, and everything in §3.1–3.2 follows from it: no interior maxima or minima, no stable equilibrium for a charge, and uniqueness (two solutions with the same boundary values differ by a harmonic function that is zero on the boundary, hence zero everywhere).

**Images: put a fake charge where it makes the boundary condition true.** A charge $q$ above a grounded plane. The plane is at $V = 0$. What else is at $V = 0$? The midplane between $q$ and $-q$. So in the region above the plane, "$q$ plus a grounded plane" and "$q$ plus a $-q$ below the plane" have the same charge *in the region* and the same boundary values. Once you know the solution is unique, they are the same problem, and the fake charge tells you the force, the surface charge, everything, for the half of space you care about. Guess: for a sphere the image is also a single charge, smaller than $q$ and inside.

**Separation of variables: guess a product, get ODEs, superpose.** You did this for the heat equation in F2. $V(x, y) = X(x)Y(y)$ turns one PDE into two ODEs linked by a separation constant. Laplace's equation has a twist: $X''/X = -Y''/Y$. If $V$ oscillates in $y$ ($Y'' < 0$) it must grow or decay exponentially in $x$. A harmonic function cannot wiggle in two directions at once (that would make bumps). Boundary conditions pick out the allowed wavenumbers; superposing the products with coefficients chosen by orthogonality matches the remaining boundary. In spherical coordinates the sines become Legendre polynomials and the exponentials become $r^\ell$ and $r^{-(\ell+1)}$.

**Multipoles: a Taylor expansion of the far field, organized by angular shape.** Far from a lump of charge, the potential looks like that of a point charge $Q/r$. If $Q = 0$, look closer: the next term falls like $1/r^2$ and has the angular shape $\cos\theta$ (a dipole). If that vanishes too, $1/r^3$ with the shape $(3\cos^2\theta - 1)/2$ (a quadrupole). The shapes are the Legendre polynomials because they are the only angular functions that, multiplied by $r^{-(\ell+1)}$, satisfy Laplace's equation. Guess the general term: $r^{-(\ell+1)}P_\ell(\cos\theta)$ times a number that measures the "$\ell$-th moment" of the charge.

```
   V = 0 plane             sphere in E0                multipole far field
   ------------------      ~~~~~~~~~~                  Q/r      (monopole)
      q  (real)             (  )  --> E0            + p.r/r^2   (dipole)
   ------------------      ~~~~~~~~~~               + Q_ij../r^3 (quadrupole)
     -q  (image)          V=0 on sphere; V -> -E0 z   + ...     valid for r > size
```

By the end of §3 you should be able to say which of the three tools fits a geometry in one glance: a conductor with a simple shape and a point charge means images; conductors that are coordinate surfaces mean separation of variables; a compact charge distribution seen from far away means multipoles.

## 3. Formal development

### 3.1 Harmonic functions: the mean-value property, no extrema, Earnshaw

▶ **Try it first.** Let $A(\rho)$ be the average of $V$ over the sphere of radius $\rho$ centred at $\mathbf r_0$. Differentiate with respect to $\rho$ and recognize a flux.

**Theorem 3.1 (Mean-value property).** If $\nabla^2V = 0$ throughout a ball of radius $R$ centred at $\mathbf r_0$, then
$$V(\mathbf r_0) = \frac{1}{4\pi R^2}\oint_{|\mathbf r - \mathbf r_0| = R}V\,da. \tag{3.1}$$

**Derivation.**

1. Write the average on the sphere of radius $\rho$ in terms of the solid angle: $A(\rho) = \dfrac{1}{4\pi}\displaystyle\oint V(\mathbf r_0 + \rho\hat{\mathbf n})\,d\Omega$.
2. Differentiate under the integral: $\dfrac{dA}{d\rho} = \dfrac{1}{4\pi}\displaystyle\oint\hat{\mathbf n}\cdot\nabla V\,d\Omega = \dfrac{1}{4\pi\rho^2}\oint_{S_\rho}\nabla V\cdot d\mathbf a$, since $da = \rho^2\,d\Omega$.
3. Divergence theorem (EM-1 Theorem 3.8): $\displaystyle\oint_{S_\rho}\nabla V\cdot d\mathbf a = \int_{\text{ball}}\nabla^2V\,d\tau = 0$. So $dA/d\rho = 0$.
4. $A$ is constant in $\rho$, and $A(\rho)\to V(\mathbf r_0)$ as $\rho\to0$ by continuity. Hence $A(R) = V(\mathbf r_0)$. ∎

■ **What just happened.** Step 3 is Gauss's law with no enclosed charge: the flux of $\mathbf E = -\nabla V$ through any sphere in a charge-free region is zero, so the average of $V$ over the sphere does not change as the sphere grows. Exercise A1 checks (3.1) by hand on $x^2 + y^2 - 2z^2$, and shows that for the non-harmonic $x^2 + y^2 + z^2$ the average exceeds the centre value by $R^2$.

**Corollary (no interior extrema).** A harmonic function on a connected region attains its maximum and minimum only on the boundary, unless it is constant. *Proof.* If $V(\mathbf r_0)$ were a maximum at an interior point, then $V\le V(\mathbf r_0)$ on every small sphere around $\mathbf r_0$ with average equal to $V(\mathbf r_0)$ by (3.1), which forces $V = V(\mathbf r_0)$ on every such sphere, so $V$ is constant on a ball; repeating the argument from every point of that ball spreads the constant over the connected region. ∎

**Theorem 3.2 (Earnshaw).** A charge cannot be held in stable equilibrium by electrostatic forces alone.

*Proof.* Let $V$ be the potential of all the *other* charges; near the test charge's position $\mathbf r_0$ there is no other charge, so $V$ is harmonic there. Stable equilibrium of $q > 0$ requires $qV$ to have a strict local minimum at $\mathbf r_0$, i.e. $V$ has an interior minimum. The corollary forbids it. (For $q<0$, replace minimum by maximum.) ∎

Two other ways to see the same thing. Flux: a restoring force in every direction means $\mathbf E$ points inward all over a small sphere around $\mathbf r_0$, so the flux is negative and Gauss's law demands negative charge inside, but there is none. Second-derivative test (LA-8 Theorem 3.19): $\nabla^2V = \operatorname{tr}(\text{Hessian}) = 0$, so the Hessian's eigenvalues cannot all be positive; the best you can do is a saddle. Consequences: a classical atom cannot be a static arrangement of charges, so classical electrons must orbit, so (EM-8) they radiate and the atom collapses in $10^{-11}$ s, which is QM-1's starting problem. And ion traps must cheat: the Paul trap uses a radio-frequency saddle that flips faster than the ion can fall off it; the Penning trap adds a magnetic field.

### 3.2 Uniqueness theorems

▶ **Try it first.** Let $W = V_1 - V_2$ be the difference of two solutions with the same boundary values. Apply Green's first identity (EM-1 (3.22)) with $\varphi = \psi = W$.

**Theorem 3.3 (Dirichlet uniqueness).** Let $\mathcal V$ be a bounded region with boundary $S$ (possibly several pieces), and let $\rho$ be given in $\mathcal V$. If $V_1$ and $V_2$ both satisfy $\nabla^2V = -\rho/\varepsilon_0$ in $\mathcal V$ and $V_1 = V_2$ on $S$, then $V_1 = V_2$ throughout $\mathcal V$. The same holds for an unbounded region if both solutions tend to $0$ at infinity.

**Derivation.**

1. $W = V_1 - V_2$ satisfies $\nabla^2W = 0$ in $\mathcal V$ (the $\rho$ cancels) and $W = 0$ on $S$.
2. Green's first identity with $\varphi = \psi = W$: $\displaystyle\int_{\mathcal V}\left(W\nabla^2W + |\nabla W|^2\right)d\tau = \oint_SW\,\nabla W\cdot d\mathbf a$.
3. The first term of the integrand is zero by step 1; the right side is zero because $W = 0$ on $S$. So $\int_{\mathcal V}|\nabla W|^2\,d\tau = 0$.
4. A non-negative continuous integrand with zero integral is zero: $\nabla W = 0$, so $W$ is constant, and the constant is its boundary value $0$. For an unbounded region, take $S$ a large sphere: $W\sim 1/r$, $\nabla W\sim1/r^2$, area $\sim r^2$, so the surface term $\to0$. ∎

**Theorem 3.4 (Uniqueness with conductors).** Let $\mathcal V$ be the region outside a set of conductors (and inside an outer boundary on which $V$ is given, or extending to infinity where $V\to0$). Given $\rho$ in $\mathcal V$ and the **total charge** $Q_i$ on each conductor, the field $\mathbf E$ in $\mathcal V$ is unique.

**Derivation.**

1. $W = V_1 - V_2$ is harmonic in $\mathcal V$, zero on the outer boundary, and on each conductor $i$ it is some constant $W_i$ (each solution is an equipotential there, EM-2 Theorem 3.4(iv), but we do not know the values).
2. Green's first identity as before: $\displaystyle\int_{\mathcal V}|\nabla W|^2\,d\tau = \sum_i\oint_{S_i}W\,\nabla W\cdot d\mathbf a = \sum_iW_i\oint_{S_i}\nabla W\cdot d\mathbf a$, pulling the constant $W_i$ out of each surface integral.
3. $\displaystyle\oint_{S_i}\nabla W\cdot d\mathbf a = -\oint_{S_i}(\mathbf E_1 - \mathbf E_2)\cdot d\mathbf a = \mp\frac{Q_i - Q_i}{\varepsilon_0} = 0$ by Gauss's law, whichever way $d\mathbf a$ points: both solutions enclose the same charge $Q_i$.
4. So $\int|\nabla W|^2 = 0$ and $\mathbf E_1 = \mathbf E_2$. ∎

■ **What just happened.** You never needed to know how the charge arranges itself on each conductor. The theorem says: if you can find *any* field that has the right charge in the region, the right total charge on each conductor, and the right behaviour at the outer boundary, it is *the* field. That licence is the method of images (§3.3), the "guess the form of the answer" moves of §3.5–3.7, EM-2's cavity theorem (Theorem 3.5) and off-centre cavity (Example 4.6), and EM-1's uniqueness step in the Helmholtz theorem. Exercise B1 proves Green's reciprocity theorem from the second identity and uses it to get induced charges with no integration.

### 3.3 The method of images

**Point charge above a grounded plane.** Charge $q$ at $(0, 0, d)$; the plane $z = 0$ is a grounded conductor. Find $V$ for $z\ge0$.

▶ **Try it first.** Write the potential of $q$ and of an imaginary $-q$ at $(0,0,-d)$. Check that it is zero on $z = 0$ and that in $z>0$ its only charge is $q$. Say which theorem makes it the answer.

1. **Guess.** For $z\ge0$,
   $$V(x, y, z) = kq\left[\frac{1}{\sqrt{x^2 + y^2 + (z-d)^2}} - \frac{1}{\sqrt{x^2 + y^2 + (z+d)^2}}\right]. \tag{3.2}$$
2. **Check.** In the region $z>0$ the only singularity of (3.2) is $q$ at $(0,0,d)$, so it solves Poisson's equation with the right source. On $z = 0$ the two square roots are equal, so $V = 0$. As $r\to\infty$, $V\to0$. By Theorem 3.3 this is the unique solution for $z\ge0$. For $z<0$, inside or beyond the conductor, the real potential is $0$; the image is a device for the upper region only.
3. **Induced surface charge.** $\sigma = -\varepsilon_0\partial V/\partial z$ at $z = 0$ ($\hat{\mathbf n} = +\hat{\mathbf z}$ points out of the metal). Differentiating (3.2), with $s^2 = x^2 + y^2$,
   $$\frac{\partial V}{\partial z}\Big|_{z=0} = kq\left[\frac{d}{(s^2 + d^2)^{3/2}} + \frac{d}{(s^2+d^2)^{3/2}}\right] = \frac{2kqd}{(s^2+d^2)^{3/2}},\qquad
   \sigma(s) = -\frac{qd}{2\pi(s^2 + d^2)^{3/2}}. \tag{3.3}$$
   Negative everywhere, largest directly below the charge ($\sigma(0) = -q/2\pi d^2$).
4. **Total induced charge.** $\displaystyle Q_{\text{ind}} = \int_0^\infty\sigma(s)\,2\pi s\,ds = -qd\int_0^\infty\frac{s\,ds}{(s^2+d^2)^{3/2}} = -qd\left[-\frac{1}{\sqrt{s^2+d^2}}\right]_0^\infty = -q.$ Every field line from $q$ ends on the plane. Stopping the integral at $s = S$ gives $-q\left(1 - d/\sqrt{S^2 + d^2}\right)$: half the induced charge lies inside $s = \sqrt3\,d$.
5. **Force on $q$.** In $z>0$, $V = V_q + V_{\text{induced}}$ and also $V = V_q + V_{\text{image}}$, so $V_{\text{induced}} = V_{\text{image}}$ everywhere in $z>0$, including at the position of $q$. The force on $q$ is therefore the Coulomb force of the image, at distance $2d$:
   $$\mathbf F = -\frac{kq^2}{(2d)^2}\hat{\mathbf z} = -\frac{kq^2}{4d^2}\hat{\mathbf z}. \tag{3.4}$$
   Cross-check: integrate the surface pressure $\sigma^2/2\varepsilon_0$ (EM-2 (3.17)) over the plane; the project does it and gets the same number.
6. **Energy, with the $\tfrac12$.** Two real charges $\pm q$ a distance $2d$ apart have interaction energy $-kq^2/2d$ (EM-2 (3.14)). The image configuration has **half** that. Honest derivation: the work to bring $q$ in from infinity to height $d$ against (3.4) is
   $$W = \int_\infty^d\frac{kq^2}{4z^2}\,dz = -\frac{kq^2}{4d}. \tag{3.5}$$
   Why half: with two real charges, both move and the field fills all space; here only $q$ moves (the induced charge arrives from ground at $V = 0$, costing nothing) and the field exists only in $z>0$. In field-energy language, $\tfrac{\varepsilon_0}{2}\int E^2$ over the half-space is half the integral for the symmetric pair, so the interaction part is halved. Self-energies (EM-2 §3.8) are excluded, as always for point charges. ∎

■ **What just happened.** One fake charge gave the potential, the induced charge density, the total induced charge, the force, and the energy, and the only justification needed was Theorem 3.3. The recipe is general: find a charge configuration *outside the region* that makes the boundary an equipotential at the right value. Example 4.1 puts numbers on all of this.

**Point charge outside a grounded sphere.** Charge $q$ at distance $a$ from the centre of a grounded sphere of radius $R < a$. Try a single image $q'$ on the line from the centre to $q$, at distance $b<R$ from the centre.

▶ **Try it first.** Impose $V = 0$ at the two points of the sphere on the axis (nearest and farthest from $q$). Two equations, two unknowns $q', b$.

1. With $\theta$ measured from the axis through $q$, a point of the sphere is at distance $\sqrt{a^2 + R^2 - 2aR\cos\theta}$ from $q$ and $\sqrt{b^2 + R^2 - 2bR\cos\theta}$ from $q'$ (law of cosines). On the axis, $\theta = 0$ and $\theta = \pi$:
   $$\frac{q}{a - R} + \frac{q'}{R - b} = 0,\qquad \frac{q}{a + R} + \frac{q'}{R + b} = 0.$$
2. Divide: $\dfrac{a+R}{a-R} = \dfrac{R+b}{R-b}$. Cross-multiply: $aR - ab + R^2 - Rb = aR + ab - R^2 - Rb$, so $2R^2 = 2ab$, i.e. $b = R^2/a$.
3. Then $q' = -q\dfrac{R - b}{a - R} = -q\dfrac{R(1 - R/a)}{a - R} = -\dfrac{qR}{a}$:
   $$q' = -\frac{R}{a}\,q\quad\text{at}\quad b = \frac{R^2}{a}. \tag{3.6}$$
4. **Verify for every $\theta$.** With $b = R^2/a$: $b^2 + R^2 - 2bR\cos\theta = \dfrac{R^2}{a^2}\left(R^2 + a^2 - 2aR\cos\theta\right)$. So the image's distance is $(R/a)$ times $q$'s distance, and
   $$\frac{kq'}{(R/a)\sqrt{a^2 + R^2 - 2aR\cos\theta}} = -\frac{kq}{\sqrt{a^2 + R^2 - 2aR\cos\theta}},$$
   which cancels $q$'s potential at every point of the sphere. (The project checks $V$ at 1000 random points of the sphere: $|V| < 10^{-13}$ V for $q = 1$ nC.) By Theorem 3.3, (3.6) gives the exterior potential.
5. **Force.** Coulomb's law between $q$ and $q'$ at separation $a - b = (a^2 - R^2)/a$:
   $$F = \frac{kqq'}{(a-b)^2} = -\frac{kq^2Ra}{(a^2 - R^2)^2}. \tag{3.7}$$
   Near contact ($a\to R$) this is $-kq^2/4(a-R)^2$, the plane result (3.4) with $d = a - R$; far away it is $-kq^2R/a^3$, attractive and falling faster than Coulomb because the induced charge $-qR/a$ shrinks as $q$ recedes. ∎

■ **What just happened.** The image of a point in a sphere is the geometric inversion of the point ($ab = R^2$), with charge scaled by $-R/a$. Exercise B2 derives $\sigma(\theta)$, integrates it to $q'$, and handles the neutral (ungrounded) sphere by adding a second image $+qR/a$ at the centre, which Theorem 3.4 licenses because it restores the total charge to zero while keeping the sphere an equipotential. Exercise C1 chains images between two planes; EM-2 Example 4.6 (charge in a spherical cavity) is the same construction turned inside out.

### 3.4 Separation of variables in Cartesian coordinates

**The slot.** Two grounded plates at $y = 0$ and $y = a$ extend to $x\to\infty$ (and infinitely in $z$); the strip at $x = 0$ is held at $V_0(y)$. Find $V(x,y)$ in $0<y<a$, $x>0$.

▶ **Try it first.** Substitute $V = X(x)Y(y)$ into $V_{xx} + V_{yy} = 0$, divide by $XY$, and argue that each piece is constant. Decide the sign of the constant from the conditions at $y = 0$ and $y = a$.

1. $X''Y + XY'' = 0 \Rightarrow X''/X = -Y''/Y$. The left side depends only on $x$, the right only on $y$; a function of $x$ equal to a function of $y$ is a constant, call it $k^2$ (the sign is decided next).
2. $Y'' = -k^2Y$ must vanish at $y = 0$ and $y = a$. If the constant were negative, $Y$ would be a combination of real exponentials, which vanishes at two points only if it is identically zero. So $k^2>0$, $Y = A\sin ky + B\cos ky$, and $X'' = +k^2X$ gives $X = Ce^{kx} + De^{-kx}$. **This is why hyperbolic and exponential functions appear:** a harmonic function that oscillates in one direction must grow or decay in the other (no interior extrema).
3. $Y(0) = 0 \Rightarrow B = 0$; $Y(a) = 0\Rightarrow ka = n\pi$, $n = 1, 2, \dots$; $V\to0$ as $x\to\infty\Rightarrow C = 0$.
4. Superpose (Laplace's equation is linear and the boundary conditions so far are homogeneous):
   $$V(x,y) = \sum_{n=1}^\infty C_n\,e^{-n\pi x/a}\sin\frac{n\pi y}{a}. \tag{3.8}$$
5. The last condition: $\sum_nC_n\sin(n\pi y/a) = V_0(y)$ on $0<y<a$. This is the half-range sine series of F2 §3.5. Multiply by $\sin(n'\pi y/a)$, integrate over $(0, a)$, and use $\int_0^a\sin\frac{n\pi y}{a}\sin\frac{n'\pi y}{a}dy = \frac a2\delta_{nn'}$:
   $$C_n = \frac2a\int_0^aV_0(y)\sin\frac{n\pi y}{a}\,dy. \tag{3.9}$$
∎

■ **What just happened.** Two properties of the family $\{\sin(n\pi y/a)\}$ did the work. **Orthogonality** let us solve for each $C_n$ separately (F2 Theorem 3.2, and behind it F2 Theorem 5.2: they are eigenfunctions of the self-adjoint operator $d^2/dy^2$ with Dirichlet conditions). **Completeness** (F2 took it on faith; Stein & Shakarchi ch. 3 proves it) guarantees that *any* reasonable $V_0$ can be matched, so the method never gets stuck. Then Theorem 3.3 says what we built is the solution. Example 4.3 does $V_0 = $ const; Exercise B8 does a parabola.

**Finite box.** If the slot is closed by a wall at $x = b$ held at $V_0(y)$, with $x = 0$ grounded, the $x$ solutions must vanish at $x = 0$: $X = \sinh(n\pi x/a)$, and
$$V = \sum_nC_n\sinh\frac{n\pi x}{a}\sin\frac{n\pi y}{a},\qquad C_n = \frac{2}{a\sinh(n\pi b/a)}\int_0^aV_0(y)\sin\frac{n\pi y}a\,dy. \tag{3.10}$$
Exercise B3 works this out and finds the centre of a square pipe with one hot wall sits at exactly $V_0/4$. A three-dimensional box separates into $X(x)Y(y)Z(z)$ with two oscillating directions and one hyperbolic one, $k_z^2 = k_x^2 + k_y^2$; the double Fourier series of F2 does the rest.

### 3.5 Spherical coordinates with azimuthal symmetry: the Legendre equation and why $\ell$ is an integer

For problems symmetric about the $z$ axis, $V = V(r, \theta)$ and Laplace's equation (EM-1 §3.10, dropping the $\phi$ term) is
$$\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial V}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial V}{\partial\theta}\right) = 0. \tag{3.11}$$

▶ **Try it first.** Put $V = R(r)\Theta(\theta)$, multiply through by $r^2/R\Theta$, and separate. Solve the radial ODE with a power of $r$.

1. Separation gives $\dfrac1R\dfrac{d}{dr}\left(r^2\dfrac{dR}{dr}\right) = -\dfrac{1}{\Theta\sin\theta}\dfrac{d}{d\theta}\left(\sin\theta\dfrac{d\Theta}{d\theta}\right) = $ const. Call the constant $\ell(\ell+1)$; for now $\ell$ is any number that makes this true (every real constant $\ge -\tfrac14$ is of this form with real $\ell$), and step 5 of the next derivation will force $\ell$ to be a non-negative integer.
2. **Radial equation.** $r^2R'' + 2rR' - \ell(\ell+1)R = 0$ is an Euler equation: try $R = r^p$. Then $p(p-1) + 2p - \ell(\ell+1) = p^2 + p - \ell(\ell+1) = (p - \ell)(p + \ell + 1) = 0$, so
   $$R(r) = A\,r^\ell + \frac{B}{r^{\ell+1}}. \tag{3.12}$$
   (For $\ell = 0$: a constant and $1/r$, the point charge.)
3. **Angular equation.** Change variable to $x = \cos\theta$, so $\dfrac{d}{d\theta} = -\sin\theta\dfrac{d}{dx}$ and $\dfrac{1}{\sin\theta}\dfrac{d}{d\theta}\left(\sin\theta\dfrac{d\Theta}{d\theta}\right) = \dfrac{d}{dx}\left[(1 - x^2)\dfrac{dP}{dx}\right]$ with $P(x) = \Theta(\theta)$. The angular equation becomes **Legendre's equation**
   $$\frac{d}{dx}\left[(1 - x^2)\frac{dP}{dx}\right] + \ell(\ell+1)P = 0,\qquad\text{i.e.}\qquad (1-x^2)P'' - 2xP' + \ell(\ell+1)P = 0. \tag{3.13}$$

**Solving Legendre's equation by power series (Frobenius), and the quantization of $\ell$.**

▶ **Try it first.** Substitute $P = \sum_{k\ge0}a_kx^k$ into (3.13) and collect the coefficient of $x^k$. You should get a two-term recursion.

1. Substituting, the coefficient of $x^k$ is $(k+2)(k+1)a_{k+2} - \left[k(k-1) + 2k - \ell(\ell+1)\right]a_k$. Since $k(k+1) - \ell(\ell+1) = (k - \ell)(k + \ell + 1)$,
   $$a_{k+2} = \frac{(k - \ell)(k + \ell + 1)}{(k+1)(k+2)}\,a_k. \tag{3.14}$$
2. The recursion links $a_0\to a_2\to a_4\to\cdots$ and $a_1\to a_3\to\cdots$ separately. So there are two independent solutions: an **even** series ($a_0 = 1$, $a_1 = 0$) and an **odd** series ($a_0 = 0$, $a_1 = 1$). The point $x = 0$ is $\theta = \pi/2$ and nothing special happens there; the danger is at the poles $x = \pm1$.
3. For large $k$, $\dfrac{a_{k+2}}{a_k} = \dfrac{(k-\ell)(k+\ell+1)}{(k+1)(k+2)} = 1 - \dfrac{2}{k} + O(k^{-2})$. Hence $\dfrac{(k+2)a_{k+2}}{k\,a_k} = \left(1 + \dfrac2k\right)\left(1 - \dfrac2k + O(k^{-2})\right) = 1 + O(k^{-2})$, whose product converges: $k\,a_k\to C\ne0$. The coefficients decay only like $C/k$.
4. At $x = \pm1$ every term is $\pm a_k$, of size $C/k$. The partial sum through $x^N$ of either series behaves like $\frac C2\ln N$ (a sum of $1/k$ over every other integer): **the series diverges at the poles**. So a non-terminating series solution is infinite along the whole $z$ axis, where there is no charge to produce an infinity. It must be rejected. (Project Part 1 watches this: for $\ell = 2.5$ the partial sums at $x = 1$ grow like $0.4919\ln N$, and $C/2 = 0.4918$.)
5. The only escape is for the series to **terminate**: (3.14) gives $a_{k+2} = 0$ for all later $k$ exactly when $k = \ell$ for some $k$ in the series, i.e. when $\ell$ is a non-negative integer with the parity of the series. The other-parity series does not terminate and must be switched off. Conclusion: $\ell\in\{0, 1, 2, \dots\}$, and the acceptable solution $P_\ell$ is a polynomial of degree $\ell$ and parity $(-1)^\ell$, normalized by $P_\ell(1) = 1$. ∎

■ **What just happened.** Nothing about electrostatics forced $\ell$ to be an integer. The requirement that a solution be finite at the poles did. **This is exactly the argument QM-5 uses** on the same equation (and its $m\ne0$ generalization) to show that the orbital angular momentum operator $L^2$ has eigenvalues $\hbar^2\ell(\ell+1)$ with integer $\ell$: there the Frobenius argument on functions of $\theta$ is what excludes the half-integer values that the ladder-operator algebra alone would allow. Do this derivation until you can reproduce steps 3–5 from memory; it is the one piece of this module that students skip and then need twice later. Exercise B4 runs (3.14) by hand for $\ell = 3$ and watches the rejected $\ell = 1$ even series sum to $1 - x\operatorname{artanh}x$, which blows up at $x = \pm1$.

**The first few.** From (3.14): $\ell = 0$: $P_0 = 1$. $\ell = 1$: $P_1 = x$. $\ell = 2$ (even, $a_0 = 1$): $a_2 = \frac{(0-2)(0+3)}{1\cdot2}a_0 = -3$, $a_4 = 0$; so $P\propto 1 - 3x^2$, and normalizing to $P(1) = 1$ gives
$$P_0 = 1,\quad P_1 = x,\quad P_2 = \tfrac12(3x^2 - 1),\quad P_3 = \tfrac12(5x^3 - 3x),\quad P_4 = \tfrac18(35x^4 - 30x^2 + 3).$$
Each $P_\ell$ has $\ell$ zeros in $(-1, 1)$, so $P_\ell(\cos\theta)$ changes sign $\ell$ times between the poles: these are the angular "shapes" of §2. LA-8 Example 4.1 produced the same polynomials by Gram–Schmidt on $1, x, x^2, \dots$: the same objects, reached from orthogonality instead of from an ODE.

**General solution** with azimuthal symmetry, by superposition:
$$V(r, \theta) = \sum_{\ell=0}^\infty\left(A_\ell\,r^\ell + \frac{B_\ell}{r^{\ell+1}}\right)P_\ell(\cos\theta). \tag{3.15}$$

### 3.6 The Legendre toolkit: Rodrigues, Bonnet, the generating function, orthogonality

Three formulas, two stated and one derived in part, and then the theorem that makes (3.15) usable.

**Rodrigues' formula** (stated; proof in Arfken & Weber §12.1 or Exercise C2 route):
$$P_\ell(x) = \frac{1}{2^\ell\,\ell!}\frac{d^\ell}{dx^\ell}\left(x^2 - 1\right)^\ell. \tag{3.17}$$
*Verify for $\ell = 2$.* $(x^2 - 1)^2 = x^4 - 2x^2 + 1$; second derivative $12x^2 - 4$; divide by $2^2\cdot2! = 8$: $\tfrac12(3x^2 - 1) = P_2$. ✓ (Project Part 0 checks (3.17) against the recursion for $\ell\le10$.)

**Bonnet's recursion** (stated; derived in Exercise C2(a)):
$$(\ell+1)P_{\ell+1}(x) = (2\ell+1)\,x\,P_\ell(x) - \ell\,P_{\ell-1}(x). \tag{3.18}$$
Check at $\ell = 1$: $2P_2 = 3x\cdot x - 1$. ✓ This is how the project generates the $P_\ell$ in three lines.

**The generating function.** For $|t| < 1$ and $|x|\le1$,
$$\frac{1}{\sqrt{1 - 2xt + t^2}} = \sum_{\ell=0}^\infty P_\ell(x)\,t^\ell. \tag{3.19}$$
*Where it comes from.* Expand with the binomial series $(1+u)^{-1/2} = 1 - \tfrac12u + \tfrac38u^2 - \cdots$ and $u = t^2 - 2xt$:
$$1 - \tfrac12(t^2 - 2xt) + \tfrac38(t^2 - 2xt)^2 + O(t^3) = 1 + x\,t + \left(\tfrac32x^2 - \tfrac12\right)t^2 + O(t^3) = P_0 + P_1t + P_2t^2 + \cdots$$
The full proof that every coefficient is a Legendre polynomial is Exercise C2: the function on the left satisfies a PDE whose $t^\ell$ coefficient is Legendre's equation (3.13), and its coefficients obey (3.18). Legendre found his polynomials as these coefficients in 1782; §3.9 shows why, since with $x = \cos\alpha$ and $t = r'/r$ the left side is $r/\mathscr r$.

▶ **Try it first.** Write (3.13) as $\mathcal LP = -\ell(\ell+1)P$ with $\mathcal L = \frac{d}{dx}(1-x^2)\frac{d}{dx}$. Integrate $\int_{-1}^1f\,\mathcal Lg\,dx$ by parts once and look at the boundary term.

**Theorem 3.5 (Orthogonality and normalization).**
$$\int_{-1}^1P_\ell(x)P_{\ell'}(x)\,dx = \frac{2}{2\ell+1}\,\delta_{\ell\ell'},\qquad\text{equivalently}\qquad\int_0^\pi P_\ell(\cos\theta)P_{\ell'}(\cos\theta)\sin\theta\,d\theta = \frac{2}{2\ell+1}\delta_{\ell\ell'}. \tag{3.16}$$

**Derivation.**

1. Legendre's equation is $\mathcal LP_\ell = -\ell(\ell+1)P_\ell$ with $\mathcal L = \dfrac{d}{dx}\left[(1-x^2)\dfrac{d}{dx}\right]$: a Sturm–Liouville operator with $p(x) = 1 - x^2$, $q = 0$, weight $w = 1$ (F2 §5, "Sturm–Liouville").
2. $\mathcal L$ is self-adjoint on functions finite at $\pm1$. Integrate by parts:
   $$\int_{-1}^1f\,\mathcal Lg\,dx = \Big[f\,(1-x^2)\,g'\Big]_{-1}^{1} - \int_{-1}^1(1-x^2)f'g'\,dx = -\int_{-1}^1(1-x^2)f'g'\,dx.$$
   The boundary term dies because $p(\pm1) = 0$: no boundary conditions had to be imposed, which is the role the factor $(1-x^2)$ plays. The right side is symmetric in $f, g$, so $\int f\,\mathcal Lg = \int(\mathcal Lf)\,g$.
3. F2 Theorem 5.2 (LA-8 Theorem 3.10): eigenfunctions of a self-adjoint operator with distinct eigenvalues are orthogonal. Explicitly, $-\ell(\ell+1)\int P_\ell P_{\ell'} = \int(\mathcal LP_\ell)P_{\ell'} = \int P_\ell(\mathcal LP_{\ell'}) = -\ell'(\ell'+1)\int P_\ell P_{\ell'}$, and $\ell(\ell+1)\ne\ell'(\ell'+1)$ for distinct non-negative integers, so the integral vanishes.
4. **Normalization from the generating function.** Square (3.19) and integrate over $x$; cross terms vanish by step 3:
   $$\int_{-1}^1\frac{dx}{1 - 2xt + t^2} = \sum_\ell t^{2\ell}\int_{-1}^1P_\ell^2\,dx.$$
   The left side is $\left[-\dfrac{\ln(1 - 2xt + t^2)}{2t}\right]_{-1}^{1} = \dfrac{1}{2t}\left[\ln(1+t)^2 - \ln(1-t)^2\right] = \dfrac1t\ln\dfrac{1+t}{1-t} = \dfrac2t\left(t + \dfrac{t^3}3 + \dfrac{t^5}5 + \cdots\right) = \sum_\ell\dfrac{2}{2\ell+1}t^{2\ell}$. Match coefficients of $t^{2\ell}$: $\int P_\ell^2 = 2/(2\ell+1)$. ∎

■ **What just happened.** Orthogonality was the spectral theorem, not a computation (step 3 never evaluated an integral). The normalization used the generating function as a bookkeeping device. The project asserts (3.16) symbolically for $\ell,\ell'\le6$ and checks the series in step 4. With (3.16), expanding boundary data in Legendre polynomials is a projection, exactly as (3.9) was: $f(\theta) = \sum A_\ell P_\ell(\cos\theta)$ with $A_\ell = \frac{2\ell+1}{2}\int_0^\pi f(\theta)P_\ell(\cos\theta)\sin\theta\,d\theta$. Exercise A4 does it for $\cos^2\theta$.

**Three recipes from (3.15) and (3.16).**

*Sphere of radius $R$ with $V(R, \theta) = V_0(\theta)$ specified; find $V$ inside.* Finite at the origin kills $B_\ell$. At $r = R$: $\sum A_\ell R^\ell P_\ell = V_0(\theta)$, so
$$A_\ell = \frac{2\ell+1}{2R^\ell}\int_0^\pi V_0(\theta)P_\ell(\cos\theta)\sin\theta\,d\theta. \tag{3.20}$$
*Outside the same sphere.* $V\to0$ kills $A_\ell$, and
$$B_\ell = \frac{(2\ell+1)R^{\ell+1}}{2}\int_0^\pi V_0(\theta)P_\ell(\cos\theta)\sin\theta\,d\theta. \tag{3.21}$$
*Spherical shell of radius $R$ carrying $\sigma(\theta)$, nothing else.* Inside $V_{\text{in}} = \sum A_\ell r^\ell P_\ell$, outside $V_{\text{out}} = \sum B_\ell r^{-(\ell+1)}P_\ell$. Two matching conditions (EM-2 §3.7): $V$ continuous at $R$, so $B_\ell = A_\ell R^{2\ell+1}$; and the normal derivative jumps by $-\sigma/\varepsilon_0$ (EM-2 (3.12)): $\sum\left[-(\ell+1)B_\ell R^{-(\ell+2)} - \ell A_\ell R^{\ell-1}\right]P_\ell = -\sigma/\varepsilon_0$, i.e. $\sum(2\ell+1)A_\ell R^{\ell-1}P_\ell = \sigma/\varepsilon_0$. Project with (3.16):
$$A_\ell = \frac{1}{2\varepsilon_0R^{\ell-1}}\int_0^\pi\sigma(\theta)P_\ell(\cos\theta)\sin\theta\,d\theta,\qquad B_\ell = A_\ell R^{2\ell+1}. \tag{3.22}$$
Examples 4.4 and 4.5 use these; in both the boundary data is a single $P_1$, so a single term survives.

### 3.7 Cylindrical coordinates with no $z$ dependence

For long objects with translational symmetry along $z$, $V = V(s, \phi)$ and
$$\frac1s\frac{\partial}{\partial s}\left(s\frac{\partial V}{\partial s}\right) + \frac{1}{s^2}\frac{\partial^2V}{\partial\phi^2} = 0. \tag{3.23}$$
Separate $V = S(s)\Phi(\phi)$ and multiply by $s^2/S\Phi$: $\dfrac sS(sS')' = -\dfrac{\Phi''}{\Phi} = n^2$. Periodicity $\Phi(\phi + 2\pi) = \Phi(\phi)$ forces $n$ to be an integer (the analogue of regularity at the poles). For $n\ge1$, $\Phi = \cos n\phi, \sin n\phi$ and the radial Euler equation $s^2S'' + sS' - n^2S = 0$ has $S = s^n, s^{-n}$. For $n = 0$, $\Phi = $ const and $S = a_0 + b_0\ln s$:
$$V(s,\phi) = a_0 + b_0\ln s + \sum_{n=1}^\infty\left(a_ns^n + \frac{b_n}{s^n}\right)\left(c_n\cos n\phi + d_n\sin n\phi\right). \tag{3.24}$$
The $\ln s$ is EM-2 Example 4.2's line charge; $s^{-n}\cos n\phi$ are the two-dimensional multipoles; $s^n\cos n\phi$ are the fields regular on the axis ($s\cos\phi = x$ is a uniform field).

**Conducting cylinder in a uniform field.** Grounded cylinder of radius $a$, axis along $z$, in $\mathbf E_0 = E_0\hat{\mathbf x}$. Far away $V\to -E_0s\cos\phi$, so only $n = 1$, $\cos\phi$ terms can appear: $V = (-E_0s + b_1/s)\cos\phi$. $V = 0$ at $s = a$ gives $b_1 = E_0a^2$:
$$V = -E_0\left(s - \frac{a^2}{s}\right)\cos\phi,\qquad \sigma = -\varepsilon_0\frac{\partial V}{\partial s}\Big|_{s=a} = 2\varepsilon_0E_0\cos\phi. \tag{3.25}$$
Compare the sphere's $3\varepsilon_0E_0\cos\theta$ (Example 4.4): the cylinder's induced "two-dimensional dipole" $b_1 = E_0a^2$ gives a surface field $2E_0$ at the sides. EM-2's project Part 6 found $0.989\times2\varepsilon_0E_0\cos\phi$ by relaxation on a square grid; now you know the exact answer it was approaching. (Project Part 0 checks (3.25) symbolically.) With $z$ dependence the radial functions become Bessel functions, which this curriculum meets in QM-6; Arfken & Weber ch. 11 is the reference.

### 3.8 Without azimuthal symmetry: spherical harmonics (stated)

Drop the symmetry: $V = R(r)Y(\theta,\phi)$. The radial equation is unchanged, (3.12), and the angular part becomes
$$\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial Y}{\partial\theta}\right) + \frac{1}{\sin^2\theta}\frac{\partial^2Y}{\partial\phi^2} = -\ell(\ell+1)\,Y. \tag{3.26}$$
Separating once more, $Y = \Theta(\theta)e^{im\phi}$ with $m$ an integer (periodicity in $\phi$), and $\Theta$ obeys the *associated* Legendre equation, (3.13) with an extra $-m^2P/(1-x^2)$ term. The same Frobenius argument gives: regular solutions exist only for integer $\ell\ge0$ and $|m|\le\ell$, and they are $P_\ell^m(x) = (1-x^2)^{|m|/2}\,d^{|m|}P_\ell/dx^{|m|}$ up to normalization. The normalized products are the **spherical harmonics** $Y_\ell^m(\theta,\phi)$.

**Theorem 3.6 (Properties of $Y_\ell^m$; stated, proved in QM-5 and Arfken & Weber ch. 12).** (i) $Y_\ell^m$ is an eigenfunction of the angular Laplacian (the left side of (3.26)) with eigenvalue $-\ell(\ell+1)$. (ii) Orthonormality: $\oint Y_\ell^{m*}Y_{\ell'}^{m'}\,d\Omega = \delta_{\ell\ell'}\delta_{mm'}$. (iii) Completeness: any reasonable function on the sphere is $\sum_{\ell m}c_{\ell m}Y_\ell^m$ with $c_{\ell m} = \oint Y_\ell^{m*}f\,d\Omega$. (iv) For each $\ell$ there are $2\ell+1$ of them, $m = -\ell, \dots, \ell$. (v) $Y_\ell^0 = \sqrt{(2\ell+1)/4\pi}\,P_\ell(\cos\theta)$, and $Y_\ell^{-m} = (-1)^m\,(Y_\ell^m)^*$.

For $\ell\le2$ (Griffiths' table, Condon–Shortley phase):
$$Y_0^0 = \frac{1}{\sqrt{4\pi}},\quad Y_1^0 = \sqrt{\frac{3}{4\pi}}\cos\theta,\quad Y_1^{\pm1} = \mp\sqrt{\frac{3}{8\pi}}\sin\theta\,e^{\pm i\phi},$$
$$Y_2^0 = \sqrt{\frac{5}{16\pi}}(3\cos^2\theta - 1),\quad Y_2^{\pm1} = \mp\sqrt{\frac{15}{8\pi}}\sin\theta\cos\theta\,e^{\pm i\phi},\quad Y_2^{\pm2} = \sqrt{\frac{15}{32\pi}}\sin^2\theta\,e^{\pm2i\phi}.$$
The general solution of Laplace's equation is then $V = \sum_{\ell m}\left(A_{\ell m}r^\ell + B_{\ell m}r^{-(\ell+1)}\right)Y_\ell^m$, and (3.15) is its $m = 0$ slice. Project Part 6 checks (i), (ii), and (v) numerically against SciPy and plots $|Y_\ell^m|^2$. **The promise:** QM-5 derives everything in Theorem 3.6 from the commutation relations of angular momentum, $[L_x, L_y] = i\hbar L_z$, without solving a single differential equation, and then identifies the eigenvalue $-\ell(\ell+1)$ here with $-L^2/\hbar^2$. Conventions differ across books by signs and factors (the Condon–Shortley $(-1)^m$ is in Griffiths, SciPy, Jackson; some chemistry texts omit it, and some use real combinations). We follow Griffiths; QM-5 §8 flags the difference again.

### 3.9 The multipole expansion

▶ **Try it first.** With $\alpha$ the angle between $\mathbf r$ and $\mathbf r'$, write $\mathscr r^2 = r^2 + r'^2 - 2rr'\cos\alpha$, factor out $r^2$, and compare with the left side of (3.19).

1. $\dfrac1{\mathscr r} = \dfrac1r\left[1 - 2\left(\dfrac{r'}r\right)\cos\alpha + \left(\dfrac{r'}r\right)^2\right]^{-1/2}$. This is the generating function (3.19) with $t = r'/r$ and $x = \cos\alpha$, valid for $r' < r$:
   $$\frac1{\mathscr r} = \frac1r\sum_{\ell=0}^\infty\left(\frac{r'}{r}\right)^\ell P_\ell(\cos\alpha)\qquad(r > r'). \tag{3.27}$$
   For $r < r'$ swap the roles: $\frac1{\mathscr r} = \frac1{r'}\sum(r/r')^\ell P_\ell(\cos\alpha)$ (used in "Your turn" for the inside of a ring).
2. Insert into EM-2 (3.8), $V(\mathbf r) = k\int\rho(\mathbf r')/\mathscr r\,d\tau'$, and interchange sum and integral (legitimate when $r$ exceeds the largest $r'$ in the source):
   $$V(\mathbf r) = k\sum_{\ell=0}^\infty\frac{1}{r^{\ell+1}}\int(r')^\ell P_\ell(\cos\alpha)\,\rho(\mathbf r')\,d\tau'. \tag{3.28}$$
   **This converges only for $r$ larger than the source.** Inside or near the source it is not an expansion of anything.
3. **Monopole, $\ell = 0$:** $V_0 = kQ/r$, $Q = \int\rho\,d\tau'$.
4. **Dipole, $\ell = 1$:** $P_1(\cos\alpha) = \cos\alpha = \hat{\mathbf r}\cdot\hat{\mathbf r}'$, so $\int r'\cos\alpha\,\rho\,d\tau' = \hat{\mathbf r}\cdot\int\mathbf r'\rho\,d\tau'$:
   $$V_1 = \frac{k\,\mathbf p\cdot\hat{\mathbf r}}{r^2},\qquad \mathbf p = \int\mathbf r'\rho(\mathbf r')\,d\tau'. \tag{3.29}$$
   This is EM-2 (4.1) for a general source; its field is EM-2 (4.3): $\mathbf E_{\text{dip}} = k\left[3(\mathbf p\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf p\right]/r^3$. For point charges, $\mathbf p = \sum q_i\mathbf r_i$.
5. **Quadrupole, $\ell = 2$:** $r'^2P_2(\cos\alpha) = \tfrac12\left(3r'^2\cos^2\alpha - r'^2\right) = \tfrac12\left(3(\hat{\mathbf r}\cdot\mathbf r')^2 - r'^2\right) = \tfrac12\sum_{ij}\hat r_i\hat r_j\left(3x'_ix'_j - r'^2\delta_{ij}\right)$, using $\hat{\mathbf r}\cdot\hat{\mathbf r} = \sum_{ij}\hat r_i\hat r_j\delta_{ij} = 1$. So
   $$V_2 = \frac{k}{2r^3}\sum_{i,j}\hat r_i\,Q_{ij}\,\hat r_j,\qquad Q_{ij} = \int\left(3x'_ix'_j - r'^2\delta_{ij}\right)\rho\,d\tau'. \tag{3.30}$$
∎

**The quadrupole tensor.** $Q_{ij}$ is symmetric, and **traceless**: $\sum_iQ_{ii} = \int(3r'^2 - 3r'^2)\rho\,d\tau' = 0$. A symmetric $3\times3$ matrix has $6$ independent entries; tracelessness removes one, leaving **5**. That is $2\ell+1$ for $\ell = 2$: the five entries of $Q_{ij}$ are five linear combinations of the five $\int r'^2\,Y_2^m(\theta',\phi')^*\rho\,d\tau'$, and $\sum\hat r_iQ_{ij}\hat r_j$ is a combination of the five $Y_2^m(\theta,\phi)$. In LA-12's language (§3.18), $Q_{ij}$ is a symmetric traceless tensor of type $(0,2)$; under rotations its five components mix only among themselves (an irreducible piece), while the trace of a general symmetric tensor would be a separate scalar ($\ell = 0$) piece. The pattern continues: $\ell = 1$ has $3$ components ($\mathbf p$, a vector), $\ell = 0$ has $1$ ($Q$, a scalar), $\ell = 3$ has $7$. Project Part 5 computes $Q_{ij}$ for a random cluster and asserts its trace is zero to machine precision; Exercise B6 does a square quadrupole by hand.

**Origin dependence.** Move the origin by $\mathbf a$, so $\mathbf r'' = \mathbf r' - \mathbf a$. Then
$$\mathbf p'' = \int(\mathbf r' - \mathbf a)\rho\,d\tau' = \mathbf p - Q\,\mathbf a. \tag{3.31}$$
If $Q = 0$ the dipole moment is intrinsic to the distribution. If $Q\ne0$ it is not: choosing $\mathbf a = \mathbf p/Q$ (the "centre of charge") makes $\mathbf p'' = 0$, just as the centre of mass makes the first mass moment vanish. Likewise $Q_{ij}$ is origin-independent only when $Q = 0$ and $\mathbf p = 0$; in general the first non-vanishing moment is the intrinsic one. Exercise A5 and Part 5 check (3.31).

**Physical dipole versus pure dipole.** Two charges $\pm q$ a distance $d$ apart form a *physical* dipole with $\mathbf p = q\mathbf d$, but they also carry an octupole ($\ell = 3$) and all higher odd moments: Example 4.6(a) finds that on the axis at $r = 5d$ the dipole term is $4\%$ low and the octupole term repairs most of that. The *pure* (point) dipole is the limit $d\to0$, $q\to\infty$ with $qd$ fixed; its potential is exactly (3.29) and its field is EM-2 (4.3) everywhere except at the origin, where the field carries an extra $-\tfrac{1}{3\varepsilon_0}\mathbf p\,\delta^3(\mathbf r)$ (Griffiths §3.4.4 and Problem 3.48; stated here, needed in EM-4 for the field inside a polarized medium). The rule: the dipole moment of a distribution is a number you compute from $\rho$; the dipole *field* is the first term of its far-field expansion, exact only for the idealized limit.

### 3.10 A dipole in an external field; dipole–dipole interaction

Let a physical dipole, $+q$ at $\mathbf r + \mathbf d/2$ and $-q$ at $\mathbf r - \mathbf d/2$, sit in an external field $\mathbf E$ (not its own) with potential $V$.

▶ **Try it first.** Write the energy, torque, and net force as differences between the two charges, then Taylor-expand to first order in $\mathbf d$.

1. **Energy.** $U = qV(\mathbf r + \mathbf d/2) - qV(\mathbf r - \mathbf d/2) \approx q\,\mathbf d\cdot\nabla V$:
   $$U = -\mathbf p\cdot\mathbf E. \tag{3.32}$$
   Lowest when $\mathbf p\parallel\mathbf E$, highest when antiparallel; the difference is $2pE$.
2. **Torque** about the dipole's centre. In a uniform field, $\mathbf N = \frac{\mathbf d}2\times q\mathbf E + \left(-\frac{\mathbf d}2\right)\times(-q\mathbf E) = q\mathbf d\times\mathbf E$:
   $$\mathbf N = \mathbf p\times\mathbf E. \tag{3.33}$$
   It turns $\mathbf p$ toward $\mathbf E$, consistent with (3.32): $N = -\partial U/\partial\theta = -pE\sin\theta$ about the axis perpendicular to both. In a non-uniform field there is an additional torque of order $d^2$ which vanishes for a pure dipole.
3. **Force.** $\mathbf F = q\mathbf E(\mathbf r + \mathbf d/2) - q\mathbf E(\mathbf r - \mathbf d/2) \approx q(\mathbf d\cdot\nabla)\mathbf E$:
   $$\mathbf F = (\mathbf p\cdot\nabla)\mathbf E. \tag{3.34}$$
   Zero in a uniform field: a dipole feels a net force only from a *gradient*. For constant $\mathbf p$ and curl-free $\mathbf E$, the identity $\nabla(\mathbf p\cdot\mathbf E) = (\mathbf p\cdot\nabla)\mathbf E + \mathbf p\times(\nabla\times\mathbf E)$ gives $\mathbf F = \nabla(\mathbf p\cdot\mathbf E) = -\nabla U$, as it must. (Exercise B7 does the point-charge field; Part 7 checks all three formulas against a physical dipole with errors falling as $d^2$.)
4. **Dipole–dipole.** Put $\mathbf p_2$ at displacement $\mathbf r$ from $\mathbf p_1$ and use (3.32) with $\mathbf E_1$ from EM-2 (4.3):
   $$U = -\mathbf p_2\cdot\mathbf E_1 = \frac{k}{r^3}\left[\mathbf p_1\cdot\mathbf p_2 - 3(\mathbf p_1\cdot\hat{\mathbf r})(\mathbf p_2\cdot\hat{\mathbf r})\right]. \tag{3.35}$$
   Head-to-tail along $\hat{\mathbf r}$: $U = -2kp_1p_2/r^3$. Side by side and parallel: $+kp_1p_2/r^3$; side by side antiparallel: $-kp_1p_2/r^3$. Symmetric in $1\leftrightarrow2$, as an interaction energy must be. ∎

■ **What just happened.** Every formula is "difference of two charges, expanded to first order," and the first-order term is where $\mathbf p = q\mathbf d$ appears. The neglected terms are $O(d^2)$ because the symmetric placement $\pm\mathbf d/2$ kills the odd corrections.

**Where this shows up.** A molecule has a permanent dipole moment when its bonds are polar and its shape does not cancel them: water ($p = 1.85$ D $= 6.17\times10^{-30}$ C m, bent), HCl ($1.08$ D); CO$_2$ and CH$_4$ have $\mathbf p = 0$ by symmetry (CO$_2$ has a quadrupole). Two water molecules head-to-tail at $0.3$ nm have, by (3.35), $U = -2kp^2/r^3 = -2.5\times10^{-20}$ J $= -0.158$ eV $\approx -6\,k_BT$ at $300$ K: the right order of magnitude for a hydrogen bond from a two-line model. Freely rotating dipoles average (3.35) to zero, but Boltzmann weighting favours attractive orientations and gives a net $\langle U\rangle\approx-\tfrac23(kp^2/r^3)^2/k_BT\propto -1/(Tr^6)$ (the Keesom term, Exercise C3). A nonpolar molecule in a field acquires an *induced* moment $\mathbf p = \alpha\mathbf E$ (EM-4), which gives the Debye ($1/r^6$) and London ($1/r^6$, quantum) terms. Together these are the van der Waals attraction, the $-1/r^6$ in the Lennard-Jones potential that SM-6 puts into a gas.

## 4. Worked examples

Every number here is printed or asserted by `project.py`; the part is named in each example.

### Example 4.1: A charge above a grounded plane, with every number (Part 2)

$q = 1$ nC at $d = 1$ cm above a grounded plane; $kq = 8.988$ V m.

- *Field at the plane directly below $q$.* Both $q$ and the image point the field downward with $kq/d^2 = 8.988/10^{-4} = 8.988\times10^4$ V/m each: $E_z = -2kq/d^2 = -1.7975\times10^5$ V/m. Check against (3.3): $\sigma = \varepsilon_0E_z = 8.854\times10^{-12}\times(-1.7975\times10^5) = -1.592\ \mu$C/m², and (3.3) at $s = 0$ gives $-q/2\pi d^2 = -10^{-9}/(2\pi\times10^{-4}) = -1.592\ \mu$C/m². ✓ At $s = d$ the density is $2^{-3/2}$ of that: $-0.563\ \mu$C/m². The project also gets $\sigma(s)$ by numerically differentiating (3.2) and matches (3.3) to $10^{-6}$.
- *Induced charge.* Inside $s<d$: $-q(1 - 1/\sqrt2) = -0.293$ nC. Half of the total is inside $s = \sqrt3d = 1.73$ cm. Total: $-1.000000$ nC by numerical integration of (3.3).
- *Force.* (3.4): $F = kq^2/4d^2 = 8.988\times10^{-9}/(4\times10^{-4}) = 2.2469\times10^{-5}$ N toward the plane. The surface-pressure integral $\int(\sigma^2/2\varepsilon_0)\,da$ gives $2.2469\times10^{-5}$ N. ✓
- *Energy.* (3.5): $W = -kq^2/4d = -2.2469\times10^{-7}$ J. The "two real charges" value $-kq^2/2d = -4.4938\times10^{-7}$ J is twice as large and wrong for this problem (§3.3 step 6).

### Example 4.2: A charge outside a grounded sphere, and the force (Part 2)

$R = 10$ cm, $q = 1$ nC at $a = 20$ cm from the centre.

- *Image.* (3.6): $q' = -qR/a = -0.50$ nC at $b = R^2/a = 5.0$ cm. Numerically, $\max|V|$ over 1000 random points of the sphere is $5.9\times10^{-14}$ V, against $kq/R\approx 90$ V scale.
- *Force.* (3.7): $F = -kq^2Ra/(a^2 - R^2)^2 = -8.988\times10^{-9}\times0.02/(0.03)^2 = -1.9972\times10^{-7}$ N (attractive). The plane formula with the same gap $a - R = 10$ cm would give $-2.2469\times10^{-7}$ N, $12\%$ larger: the sphere curves away from the charge, so less induced charge sits close to it.
- *Surface charge* (Exercise B2 derives it): $\sigma(\theta) = -\dfrac{q(a^2 - R^2)}{4\pi R\,(R^2 + a^2 - 2aR\cos\theta)^{3/2}}$. At the near pole $\sigma(0) = -q(a^2-R^2)/4\pi R(a-R)^3 = -23.873$ nC/m²; at the far pole $-0.8842$ nC/m²; ratio $\left(\frac{a+R}{a-R}\right)^3 = 27.00$. Integrated over the sphere: $-0.500000$ nC $= q'$. ✓ The project also recovers $\sigma$ from $-\varepsilon_0\partial V/\partial r$ numerically.

### Example 4.3: The slot with $V_0$ constant: the first three coefficients (Part 3)

$V_0(y) = V_0$. By (3.9), $C_n = \dfrac{2V_0}{a}\displaystyle\int_0^a\sin\frac{n\pi y}a\,dy = \frac{2V_0}{n\pi}(1 - \cos n\pi) = \begin{cases}4V_0/n\pi & n\text{ odd}\\ 0 & n\text{ even.}\end{cases}$
So $C_1 = 4V_0/\pi = 1.2732\,V_0$, $C_2 = 0$, $C_3 = 4V_0/3\pi = 0.4244\,V_0$ (the project computes them by numerical quadrature and gets these to $10^{-12}$).

*At the centre line $y = a/2$, $x = a/2$:* $\sin(n\pi/2) = \pm1$ for odd $n$, so $V/V_0 = \frac4\pi\left[e^{-\pi/2} - \tfrac13e^{-3\pi/2} + \tfrac15e^{-5\pi/2} - \cdots\right]$: one term $0.26468$, two $0.26087$, three $0.26097$. The series can be summed exactly. With $z = e^{-\pi(x - iy)/a}$, $\sum_{n\text{ odd}}z^n/n = \tfrac12\ln\frac{1+z}{1-z}$, and taking the imaginary part gives
$$V(x,y) = \frac{2V_0}{\pi}\arctan\frac{\sin(\pi y/a)}{\sinh(\pi x/a)}, \tag{4.1}$$
which at the centre is $\frac2\pi\arctan(1/\sinh(\pi/2)) = 0.26096\,V_0$. ✓ Three terms suffice here, but the series slows down near $x = 0$, where $V_0(y)$ jumps at the corners: the project needs $N = 3$, $19$, $125$ terms for $10^{-4}$ accuracy at $x = 0.5a$, $0.1a$, $0.01a$. The jump costs one power of $n$ in the coefficients (F2 Theorem 3.9) and the exponential factor is what rescues convergence away from the wall. A finite-difference relaxation solution (red-black SOR, LA-13 §3.9, $h = a/40$, 213 sweeps) gives $V(a/2,a/2) = 0.26098\,V_0$ and agrees with (4.1) to $4.4\times10^{-4}V_0$ for $x\ge a/4$; `slot.png` overlays the equipotentials.

### Example 4.4: A conducting sphere in a uniform field (Part 4)

A grounded (equivalently, neutral) conducting sphere of radius $R$ in $\mathbf E_0 = E_0\hat{\mathbf z}$. Far away $V\to -E_0z = -E_0r\cos\theta$, and $V = 0$ on $r = R$.

Use (3.15) outside the sphere. The condition at infinity allows only one growing term, $A_1 = -E_0$ with $P_1 = \cos\theta$, and $A_\ell = 0$ otherwise. At $r = R$: $(-E_0R + B_1/R^2)\cos\theta + \sum_{\ell\ne1}B_\ell R^{-(\ell+1)}P_\ell = 0$. By (3.16) each coefficient vanishes separately: $B_\ell = 0$ for $\ell\ne1$ and $B_1 = E_0R^3$:
$$V(r,\theta) = -E_0\left(r - \frac{R^3}{r^2}\right)\cos\theta\qquad(r\ge R). \tag{4.2}$$
The second term is a pure dipole potential $kp\cos\theta/r^2$ with the **induced dipole moment**
$$p = 4\pi\varepsilon_0R^3E_0. \tag{4.3}$$
Surface charge from $\sigma = -\varepsilon_0\partial V/\partial r$ at $R$: $\partial V/\partial r = -E_0(1 + 2R^3/r^3)\cos\theta\to -3E_0\cos\theta$, so
$$\sigma(\theta) = 3\varepsilon_0E_0\cos\theta. \tag{4.4}$$
Checks: total charge $\propto\int_0^\pi\cos\theta\sin\theta\,d\theta = 0$ ✓; dipole moment of the surface charge $\int(R\cos\theta)\sigma\,da = 3\varepsilon_0E_0R\cdot2\pi R^2\cdot\tfrac23 = 4\pi\varepsilon_0R^3E_0 = p$ ✓; the field at the poles is $\sigma/\varepsilon_0 = 3E_0$, three times the applied field. (Project Part 0 verifies (4.2) solves (3.11), vanishes at $R$, and gives (4.4) and (4.3), all symbolically.)

*Numbers.* $R = 1$ cm, $E_0 = 10^5$ V/m: $\sigma_{\max} = 2.656\ \mu$C/m²; pole field $3.0\times10^5$ V/m; charge on the upper hemisphere $3\pi\varepsilon_0E_0R^2 = 0.834$ nC (and $-0.834$ nC on the lower); $p = 1.1127\times10^{-11}$ C m.

*Relaxation check.* Part 4 solves (3.11) by SOR on a $72\times161$ grid in $(\theta, \ln r)$ out to $r = 50R$, with $V = 0$ on the sphere and $-E_0r\cos\theta$ on the outer boundary. After 524 sweeps: $|V - (4.2)| \le 2.7\times10^{-3}E_0R$ for $r\le3R$ (relative error $1.1\times10^{-3}$ everywhere), $\sigma$ within $5.6\times10^{-4}$ of (4.4) in units of $3\varepsilon_0E_0$, induced $p = 1.0008\times4\pi\varepsilon_0R^3E_0$, total charge $10^{-13}$. The grid never knew about Legendre polynomials; `sphere_uniform_field.png` overlays the two.

### Example 4.5: A shell with $\sigma(\theta) = \sigma_0\cos\theta$ (Part 0)

A spherical shell of radius $R$ carries $\sigma = \sigma_0\cos\theta$ and nothing else is present. In (3.22) only $\ell = 1$ survives, since $\cos\theta = P_1$ and (3.16) kills every other integral: $A_1 = \dfrac{1}{2\varepsilon_0}\displaystyle\int_0^\pi\sigma_0\cos^2\theta\sin\theta\,d\theta = \frac{\sigma_0}{2\varepsilon_0}\cdot\frac23 = \frac{\sigma_0}{3\varepsilon_0}$, $B_1 = A_1R^3$:
$$V_{\text{in}} = \frac{\sigma_0}{3\varepsilon_0}\,r\cos\theta = \frac{\sigma_0}{3\varepsilon_0}\,z,\qquad V_{\text{out}} = \frac{\sigma_0R^3}{3\varepsilon_0}\frac{\cos\theta}{r^2}. \tag{4.5}$$
*Matching.* At $r = R$ both equal $\frac{\sigma_0R}{3\varepsilon_0}\cos\theta$ ✓. Normal derivatives: $\partial V_{\text{out}}/\partial r - \partial V_{\text{in}}/\partial r = -\frac{2\sigma_0}{3\varepsilon_0}\cos\theta - \frac{\sigma_0}{3\varepsilon_0}\cos\theta = -\frac{\sigma_0\cos\theta}{\varepsilon_0} = -\sigma/\varepsilon_0$ ✓ (EM-2 (3.12)). Project Part 0 checks both and that each piece is harmonic.

*Reading.* Inside, the field is **uniform**: $\mathbf E_{\text{in}} = -\frac{\sigma_0}{3\varepsilon_0}\hat{\mathbf z}$. Outside, it is a pure dipole with $p = 4\pi\varepsilon_0B_1 = \tfrac{4\pi}3R^3\sigma_0$, which is also $\int z\,\sigma\,da$ directly. Now compare with Example 4.4: there $\sigma = 3\varepsilon_0E_0\cos\theta$, so $\sigma_0 = 3\varepsilon_0E_0$ and the induced charge alone produces $\mathbf E_{\text{in}} = -E_0\hat{\mathbf z}$ inside the sphere, exactly cancelling the applied field, as EM-2 Theorem 3.4(i) demands. EM-4 reuses this shell as the surface charge of a uniformly polarized sphere.

*Numbers.* $R = 5$ cm, $\sigma_0 = 1\ \mu$C/m²: $E_{\text{in}} = \sigma_0/3\varepsilon_0 = 3.7647\times10^4$ V/m along $-\hat{\mathbf z}$; $p = 5.2360\times10^{-10}$ C m.

### Example 4.6: Two charge pairs to quadrupole order (Part 5)

*(a) $+q$ at $(0,0,d)$ and $-q$ at $(0,0,-d)$.* Moments about the origin: $Q = 0$; $\mathbf p = q(d\hat{\mathbf z}) + (-q)(-d\hat{\mathbf z}) = 2qd\,\hat{\mathbf z}$; quadrupole $Q_{zz} = q(3d^2 - d^2) + (-q)(3d^2 - d^2) = 0$, and likewise every $Q_{ij} = 0$. For charges on the axis, $\cos\alpha = \pm1$ and (3.28) becomes $V = k\sum_\ell\left[q\,d^\ell + (-q)(-d)^\ell\right]P_\ell(\cos\theta)/r^{\ell+1}$: only **odd** $\ell$ survive. On the positive axis ($\theta = 0$) at $r = 5d$, in units of $kq/d$: exact $\frac1{5-1} - \frac1{5+1} = \frac1{12} = 0.083333$; dipole term $2/25 = 0.080000$; adding the octupole $2/5^4 = 0.0032$ gives $0.083200$; adding $\ell = 5$ gives $0.083328$. The dipole term alone is $4\%$ low at five dipole lengths, and the error is second order in $d/r$ (EM-2 Example 4.3 said the same).

*(b) Linear quadrupole: $+q$ at $(0,0,\pm d)$ and $-2q$ at the origin.* $Q = 0$; $\mathbf p = qd\hat{\mathbf z} - qd\hat{\mathbf z} + 0 = 0$: **the dipole term vanishes.** Quadrupole: $Q_{zz} = 2q(3d^2 - d^2) = 4qd^2$, $Q_{xx} = Q_{yy} = 2q(0 - d^2) = -2qd^2$, off-diagonals zero; trace $4 - 2 - 2 = 0$ ✓. By (3.30), $V_2 = \frac{k}{2r^3}\hat r_iQ_{ij}\hat r_j = \frac{kqd^2}{r^3}\left(2\cos^2\theta - \sin^2\theta\right) = \frac{2kqd^2}{r^3}P_2(\cos\theta)$, which is also what (3.28) gives directly: $k\left[2q\,d^2\right]P_2/r^3$. On the axis at $r = 5d$: exact $\frac1{4} + \frac16 - \frac25 = \frac1{60} = 0.016667$; quadrupole term $2/125 = 0.016000$; with $\ell = 4$, $0.016640$; with $\ell = 6$, $0.016666$. Only even $\ell$ appear, by the symmetry $z\to-z$.

For a random cluster of six charges inside a ball of radius $d$, Part 5 compares the exact potential with (3.28) truncated after $\ell = 0, 1, 2$ at $r = 3d$ to $100d$; the RMS errors fall with fitted slopes $-2.002$, $-3.003$, $-4.001$ on a log-log plot (`multipole_errors.png`), since the first omitted term is $\ell+1$.

## 5. Why it matters and where it connects

- **EM-4.** A uniformly polarized sphere has bound surface charge $\sigma_b = P\cos\theta$: Example 4.5 *is* its field, and (4.2)–(4.4) with the sphere's response replaced by a dielectric's are EM-4's central example. The delta-function term in the pure dipole field (§3.9) is what makes the average field in a polarized medium come out right.
- **EM-5.** The vector potential of a current loop is expanded exactly as (3.28); the magnetic dipole $\mathbf m$ plays $\mathbf p$'s role, and $\mathbf N = \mathbf m\times\mathbf B$, $U = -\mathbf m\cdot\mathbf B$, $\mathbf F = \nabla(\mathbf m\cdot\mathbf B)$ copy (3.32)–(3.34). Images work for magnetic problems with superconductors.
- **EM-8.** The radiation from a time-dependent source is organized by the same multipole expansion; electric-dipole radiation ($\propto\omega^4$) is the $\ell = 1$ term.
- **QM-5.** The angular part of the Laplacian (3.26) is $-L^2/\hbar^2$; its eigenfunctions are the $Y_\ell^m$ of Theorem 3.6 and its eigenvalue $\ell(\ell+1)$ is quantized by the Frobenius argument of §3.5. QM-5 proves Theorem 3.6 by operator algebra.
- **QM-6.** Hydrogen separates as $R(r)Y_\ell^m$; the radial equation carries the same $\ell(\ell+1)/r^2$ term that produced $r^\ell$ and $r^{-(\ell+1)}$ in (3.12), and the $2\ell+1$ degeneracy is Theorem 3.6(iv).
- **QM-8.** The energy shift of an atom in a field uses $U = -\mathbf p\cdot\mathbf E$ with $\mathbf p$ an operator (Stark effect); the quadrupole tensor's five components become the $\ell = 2$ perturbation.
- **LA-8, F2.** Theorem 3.5 is LA-8 Theorem 3.10 applied to $\mathcal L = \frac d{dx}(1-x^2)\frac d{dx}$; (3.9) and (3.20) are projections onto an orthogonal basis (LA-5, F2 Theorem 3.5). LA-8 Example 4.1 built the same $P_\ell$ by Gram–Schmidt.
- **LA-12.** $Q_{ij}$ is the first symmetric traceless tensor you meet with a physical job; CM-7's inertia tensor is its mass analogue (not traceless), and the stress tensor in EM-6 is the next.
- **LA-13.** Project Parts 3 and 4 are red-black SOR for the five-point Laplacian; LA-13 Theorem 3.19 is why it converges and (3.24) there is why Jacobi would take fifty times longer.
- **SM-6, ST.** The Keesom average of Exercise C3 is the first appearance of a thermally averaged interaction; its $1/r^6$ goes into the Lennard-Jones and van der Waals equations of state.
- **CM-5.** Legendre's original problem: the gravitational potential of an oblate planet is $-\frac{GM}r\left[1 - J_2\left(\frac Rr\right)^2P_2(\cos\theta) - \cdots\right]$, and $J_2$ (Earth: $1.08\times10^{-3}$) is a quadrupole moment that precesses satellite orbits.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 drill the mean-value property, image numbers, choosing a method, a Legendre expansion, origin-dependence of moments, and water in a field. B1–B8 are the module: B1 (Green's reciprocity), B2 (image in a sphere, complete), B4 (Frobenius by hand, including the rejected solution), B5 (charged sphere in a field), B7 (dipole in a point-charge field), and B8 (`project.py` extension) are essential. C1 (images between two planes), C2 (the generating function proves Bonnet, Legendre, and the multipole expansion), C3 (Keesom's thermally averaged dipole–dipole force).

## 7. Coding project

`project.py` (about 8 s) builds the module's tools and checks every number above:

0. **Legendre, symbolically** (SymPy). `legendre_recurrence` (Bonnet (3.18)), `legendre_rodrigues` (3.17), `frobenius_series` (3.14): all three agree with each other and with `sympy.legendre` for $\ell\le10$; each $P_\ell$ satisfies (3.13), has $P_\ell(1) = 1$ and parity $(-1)^\ell$; orthogonality (3.16) is asserted exactly for $\ell,\ell'\le6$; the generating function (3.19) is expanded to $t^{10}$ and its coefficients compared; the series of Theorem 3.5 step 4 is checked; the two identities of Exercise C2 are verified. Then (4.2)–(4.5) and (3.25) are shown to be harmonic with the right boundary values and surface charges, and the sphere averages of Exercise A1 and the coefficients of A4 are computed.
1. **Why $\ell$ is an integer.** `frob_coeffs` runs (3.14) in floating point to $k = 2\times10^5$ for $\ell = 2, 2.5, 3.5$. The $\ell = 2$ series terminates; for the others $k\,a_k\to C$ and the partial sums at $x = 1$ grow like $(C/2)\ln N$ (fitted slopes $0.4919$, $0.5665$ against $C/2 = 0.4918$, $0.5663$). `frobenius_divergence.png` shows the $\ell = 2.5$ "polynomial" blowing up at the pole.
2. **Images.** `V_pts`, `E_pts` evaluate point-charge sums by broadcasting. For the plane: $\sigma$ from $-\varepsilon_0\partial V/\partial z$ matches (3.3); total $-q$; force by (3.4) and by the pressure integral; energy (3.5). For the sphere: $V = 0$ on the surface, (3.7), $\sigma(\theta)$ and its integral, the neutral-sphere force of B2(c) and the limits of B2(d). `images.png`: field lines from the image construction and $\sigma(s)$ against (3.3).
3. **The slot.** `sine_coeffs` computes (3.9) by oscillatory quadrature; `slot_series` sums (3.8); `sor` is red-black successive over-relaxation for any five-point stencil $\sum_{\text{nb}}c_{\text{nb}}(V_{\text{nb}} - V) = 0$; `solve_slot` applies it. Example 4.3's numbers, the terms-needed table, B3's box. `slot.png`.
4. **Sphere in a uniform field** by relaxation on an $(\ln r, \theta)$ grid: multiplying $r^2\times$(3.11) by $e^{\rho}\sin\theta$ ($\rho = \ln r$) gives a symmetric five-point stencil whose coefficients are the `cE, cW, cN, cS` arrays; $\theta$ sits at cell centres so the poles need no special case. Example 4.4's relaxation numbers; B5's charged sphere. `sphere_uniform_field.png`.
5. **Multipoles.** `moments` returns $Q$, $\mathbf p$ (3.29), $Q_{ij}$ (3.30) for point charges; `V_multipole` evaluates the truncated expansion in tensor form; the $\ell = 2$ term is also computed in Legendre form (3.28) and asserted equal. Error slopes $-2, -3, -4$ asserted within $0.05$; $\operatorname{tr}Q_{ij} = 0$ to $10^{-12}$ relative; the origin shift (3.31); Example 4.6, A5, B6. `multipole_errors.png`.
6. **Spherical harmonics.** Griffiths' table for $\ell\le2$ against `scipy.special.sph_harm_y` (they agree, including the Condon–Shortley sign), eigenvalue $-\ell(\ell+1)$ symbolically, the $9\times9$ Gram matrix on a Gauss–Legendre $\times$ uniform-$\phi$ grid equal to $I$, and $|Y_\ell^m|^2$ surfaces (`spherical_harmonics.png`, the QM-5 preview).
7. **Dipoles in fields.** A physical dipole in a point charge's field: relative errors of force, torque, energy against (3.32)–(3.34) fall from $4.4\times10^{-4}$ to $2.7\times10^{-5}$ as $d$ halves twice (factor $4$ each time, asserted). Four point charges reproduce (3.35) to $10^{-6}$. Water numbers; A6's Langevin average; C3's exact Boltzmann average against the leading-order formula.

**Build it yourself** before reading the code: `sor(V, fixed, cE, cW, cN, cS, omega)` (update the two colours alternately; stop when the largest change is below `tol`), `moments(qs, pos)` (one `einsum` for $Q_{ij}$), `V_multipole`, and `sine_coeffs`. Then compare. **Your turn** asks for: the image series between two grounded planes (C1) summed numerically; the potential of a charged ring expanded in Legendre polynomials inside and outside, with the on-axis coefficients $P_{2n}(0) = 1, -\tfrac12, \tfrac38, -\tfrac5{16}, \dots$; a water molecule as three point charges against its dipole approximation; and B8's parabolic slot. Reference solutions follow the prompts; write yours first.

## 8. Where you'll struggle

- **This module is pattern recognition plus patience.** There is no algorithm that tells you "images" versus "separation" versus "multipoles"; you learn it by doing fifteen problems and noticing that conductors with one simple shape plus a point charge want images, conductors that are coordinate surfaces want separation, and "far away from a compact source" wants multipoles. Griffiths' chapter 3 problems are the right fifteen. Do not skip them because the notes felt clear.
- **The Frobenius argument is the one derivation here that pays off later, and the one students skip.** Steps 3–5 of §3.5 (coefficients $\sim C/k$, logarithmic divergence at the poles, termination forces integer $\ell$) are reused verbatim in QM-5 for the quantization of angular momentum and in QM-6 for the hydrogen radial equation (where termination of a series gives the energy levels). If you can only reproduce one derivation from this module closed-book, make it this one.
- **The $\tfrac12$ in the image energy.** $-kq^2/4d$, not $-kq^2/2d$. Derive it as the work integral (3.5) every time, or remember: the field exists in half of space and only one charge is real. The same $\tfrac12$ caught you in EM-2's $W = \tfrac12\sum q_iV_i$.
- **$\theta$ is measured from the symmetry axis, which is often not the $z$ axis of the problem as drawn.** A sphere in a field along $x$ has $\sigma = 3\varepsilon_0E_0\cos\theta'$ with $\theta'$ measured from $\hat{\mathbf x}$. Before writing any $P_\ell(\cos\theta)$, say out loud which axis $\theta$ is measured from.
- **The multipole expansion converges only for $r$ larger than the source.** Inside the source, or at a point closer than the farthest charge, (3.28) is not an expansion of anything; use the $r < r'$ form (3.27) for the interior of a shell or ring. The sign that you have forgotten: a "potential" that grows as you approach the origin.
- **The dipole moment of a distribution is not the field of a physical dipole.** $\mathbf p = \int\mathbf r'\rho\,d\tau'$ is a number you compute; EM-2 (4.3) is exact only for the point-dipole limit. Example 4.6(a) shows the physical dipole's field is $4\%$ off at $r = 5d$ and that the correction is an octupole, not a dipole error. Conversely, $\mathbf p$ depends on the origin unless $Q = 0$ (3.31); report the origin with the moment.
- **Signs of induced charge.** $\sigma = -\varepsilon_0\partial V/\partial n$ with $\hat{\mathbf n}$ out of the metal. For the plane that is $+\hat{\mathbf z}$; for a sphere seen from outside it is $+\hat{\mathbf r}$; for a cavity wall it is $-\hat{\mathbf r}$. Check every result against "the induced charge has the opposite sign to the inducing charge and is largest nearest to it."
- **Dropping terms at infinity or at the origin.** Outside a sphere you may keep an $r^\ell$ term only if the field at infinity demands it (Example 4.4 keeps exactly one); inside, no $r^{-(\ell+1)}$ at all. Students write (3.15) with all terms and then cannot determine the coefficients.
- **Spherical harmonic conventions.** Sign ($(-1)^m$ or not), normalization, and whether $Y_\ell^{-m}$ is $(Y_\ell^m)^*$ or $(-1)^m(Y_\ell^m)^*$ differ between books and software. Griffiths, Jackson, SciPy and these notes agree; Arfken & Weber and many chemistry texts do not. When two sources disagree, check $Y_1^1$ at $\theta = \phi = \pi/2$.
- **Slow Fourier convergence near a jump.** Example 4.3 needs $125$ terms for $10^{-4}$ at $x = 0.01a$. That is F2's Gibbs phenomenon, not a bug; smooth boundary data (B8) converge like $n^{-3}$.

## 9. Check yourself

1. State the mean-value property and derive it in four lines from Gauss's law.
2. Why can a charge not be in stable electrostatic equilibrium? Give the flux argument in two sentences.
3. Prove that the potential in a region is determined by its values on the boundary. Which identity does the proof use, and where does the boundary condition enter?
4. What does the conductor version of the uniqueness theorem take as given on each conductor, and why does that make the method of images legitimate?
5. Image charge for $q$ at distance $a$ from the centre of a grounded sphere of radius $R$: value and position. What changes if the sphere is neutral and isolated?
6. Why is the image-charge energy half the two-charge energy?
7. In the slot problem, why is the $x$ dependence exponential and the $y$ dependence sinusoidal, and not the other way round?
8. Separate Laplace's equation in spherical coordinates with azimuthal symmetry and write the two radial solutions for a given $\ell$.
9. Write the Frobenius recursion for Legendre's equation. Why must $\ell$ be a non-negative integer?
10. Write $\int_{-1}^1P_\ell P_{\ell'}\,dx$. Which theorem of LA-8 makes it zero for $\ell\ne\ell'$?
11. Write the potential outside a conducting sphere in a uniform field and the induced surface charge.
12. What is the first term of the multipole expansion of a neutral molecule's potential, and how does it depend on the choice of origin?
13. How many independent components does the quadrupole tensor have, and why does that number equal $2\ell+1$?
14. Write the energy, torque, and force on a dipole in an external field. In which of the three does the field gradient appear?

## 10. Pair with

- **Griffiths, *Introduction to Electrodynamics*, ch. 3.** Same order and notation as these notes. §3.3.2 (spherical separation) has four more worked boundary-value problems; §3.4 (multipoles) is slower and gentler on the origin-dependence point. Do Problems 3.9–3.12 (images), 3.18–3.25 (separation), 3.27–3.33 (multipoles): those are the fifteen.
- **Jackson, *Classical Electrodynamics*, ch. 2–3** (reference only; **this is a graduate text** and assumes you already know Griffiths cold). Ch. 2 has the Green's function formulation of images and the full theory of the sphere; ch. 3 has the associated Legendre functions, the addition theorem for spherical harmonics (which is why the $\ell = 2$ term of (3.28) is five $Y_2^m$'s), and Bessel functions. Read §2.1–2.5 and §3.1–3.5 only when you want to know *why* a formula is true at a level these notes omit.
- **MIT 8.07 lecture notes (Guth), lectures 5–9** (OCW). The uniqueness theorems and the Frobenius argument with more care about convergence than Griffiths; the multipole expansion in Cartesian tensor form, which is how (3.30) is usually written in research papers.
- **Arfken & Weber, *Mathematical Methods for Physicists*, ch. 11–12** (Bessel and Legendre functions). Proofs of Rodrigues' formula, Bonnet's recursion, and the generating function in full, the associated Legendre functions and spherical harmonics with every normalization written out, and the Sturm–Liouville theory behind Theorem 3.5 (ch. 10). Its $Y_\ell^m$ convention omits the Condon–Shortley phase; be warned.
