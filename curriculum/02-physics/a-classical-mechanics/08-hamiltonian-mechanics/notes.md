# CM-8: Hamiltonian Mechanics: Legendre Transform, Phase Space, Poisson Brackets, Liouville, Canonical Transformations, Hamilton–Jacobi

**Time:** ~12 hours · **Prerequisites:** CM-4 (generalized momenta, cyclic coordinates, the energy function) · **Unlocks:** SM-3 (Liouville and ensembles), QM-4 (commutators from Poisson brackets), QM-9 (path integrals, the classical limit), QC-9 (Hamiltonian simulation)

> **Conventions in this module.** Generalized coordinates are $q = (q_1, \dots, q_n)$ and their conjugate momenta are $p_i = \partial L/\partial \dot q_i$, exactly as in CM-4. The Hamiltonian is $H(q, p, t) = \sum_i p_i \dot q_i - L$. Hamilton's equations are $\dot q_i = \partial H/\partial p_i$, $\dot p_i = -\partial H/\partial q_i$. The Poisson bracket is $\{F, G\} = \sum_i \left(\frac{\partial F}{\partial q_i}\frac{\partial G}{\partial p_i} - \frac{\partial F}{\partial p_i}\frac{\partial G}{\partial q_i}\right)$, so that $\{q, p\} = +1$ (some books, including older editions of Landau, use the opposite sign; check before copying a formula). Phase space is the set of points $z = (q, p) \in \mathbb{R}^{2n}$. The symplectic matrix is written in sans-serif, $\mathsf{J} = \begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}$, because the italic $J = \oint p\,dq$ is reserved for the action variable; both names are standard and I will not rename either. The angle variable $w$ conjugate to $J$ advances by exactly $1$ per period, so the frequency $\nu = \partial H/\partial J$ is in cycles per unit time and $\omega = 2\pi\nu$. $L$ with no index is the Lagrangian; $L_x, L_y, L_z$ are components of angular momentum and $L^2 = L_x^2 + L_y^2 + L_z^2$. $\epsilon_{ijk}$ is the Levi-Civita symbol with $\epsilon_{xyz} = 1$. Sums are written out; no summation convention. The pendulum's angle $\theta$ is measured from the downward vertical, so $V = -mgl\cos\theta$ and the hanging rest point is $\theta = 0$. $\hbar = h/2\pi$.

## 0. Where this sits

You can write a Lagrangian, get Euler–Lagrange equations out of it, read off conserved momenta from cyclic coordinates, and you know the energy function $h = \sum_i \dot q_i\,\partial L/\partial\dot q_i - L$ from CM-4. That is everything this module needs from mechanics. From mathematics it needs multivariable calculus (the chain rule for partial derivatives, equality of mixed partials, the divergence theorem) and the ability to compute a $2\times 2$ determinant.

What it makes possible: SM-3 builds the microcanonical ensemble on Liouville's theorem from §3.5. QM-4 builds the commutator $[\hat x, \hat p] = i\hbar$ and Heisenberg's equation of motion on the Poisson bracket of §3.6. QM-8's WKB approximation is the Hamilton–Jacobi equation of §3.10 with a correction term. QM-9's path integral starts from the phase-space action of §3.3. QC-9's Trotterization is the splitting idea behind the symplectic integrators of §3.11. And ST-6's Hamiltonian Monte Carlo is a Verlet integrator run on a probability distribution.

## 1. The problem that forced this idea

Lagrange's equations are $n$ second-order ordinary differential equations in the $n$ configuration variables $q_i$. They are coordinate-free, they absorb constraints, and by 1800 they were the best tool anyone had. But they have a defect that bothered the astronomers of the early nineteenth century: when you perturb a solved problem (a planet, whose orbit you know exactly, now tugged by Jupiter), the perturbation equations for the orbital elements come out in a strange symmetric form that Lagrange's framework does not explain. Lagrange (1808) and Poisson (1809) both found this form while working on planetary perturbations, and Poisson, in the course of it, defined a bracket of two functions of the coordinates and momenta that reappears throughout this module. Neither of them saw it as the foundation of mechanics. They saw it as a trick for perturbation theory.

William Rowan Hamilton came to it from optics. In 1827, at age twenty-one, he had shown that all of geometrical optics could be packed into a single function (his *characteristic function*), from which every ray could be recovered by differentiation. Light rays are the paths that make an integral stationary (Fermat's principle, CM-2). Mechanical trajectories are also paths that make an integral stationary (Hamilton's principle, CM-2). So Hamilton asked: is there a characteristic function for mechanics? In two papers, *On a General Method in Dynamics* (1834) and its *Second Essay* (1835), he answered yes, and along the way rewrote Lagrange's $n$ second-order equations as $2n$ first-order equations,
$$
\dot q_i = \frac{\partial H}{\partial p_i}, \qquad \dot p_i = -\frac{\partial H}{\partial q_i},
$$
in which the coordinates and the momenta appear on exactly the same footing. He called the function $H$ the "principal function" of the system; Jacobi would call it the Hamiltonian.

Three things happened to this reformulation within a decade, and a fourth took ninety years.

**(a) It made the geometry of motion visible.** A state of the system is now a point $(q, p)$ in a $2n$-dimensional space, and Hamilton's equations are a vector field on that space. The whole history of a system is a curve; the whole set of possible histories is a flow, like a fluid. Questions that are invisible in configuration space (where trajectories cross each other constantly) become obvious in phase space (where they never do).

**(b) The flow preserves a structure.** Joseph Liouville, in an 1838 note on a different subject (the variation of arbitrary constants in a system of ordinary differential equations), proved that for equations of Hamilton's form the flow carries any region of phase space to a region of the same volume. Jacobi noticed what it meant for mechanics. Boltzmann (1868–1871), Maxwell (1879), and finally Gibbs (1902, who named it the "conservation of density-in-phase") turned it into the foundation of statistical mechanics: it is the reason you are allowed to say "all accessible microstates are equally likely" and have the statement stay true as time passes. Liouville's theorem is a corollary of something stronger, which Hamilton's equations also preserve and which was only named in the twentieth century: the *symplectic form*.

**(c) It gave Jacobi a method that solved problems Lagrange could not.** Carl Gustav Jacob Jacobi (1837) simplified Hamilton's optics-inspired construction into a single first-order partial differential equation for a generating function $S(q, t)$, now the Hamilton–Jacobi equation. When that equation can be separated, the problem is solved by quadratures. Jacobi used it to find the geodesics on an ellipsoid (1838) and to handle the motion of a particle attracted by two fixed centers, problems where the Lagrangian equations give you nothing you can integrate. He also proved, in 1842, the identity that now carries his name, and with it the theorem Poisson had stated in 1809: the bracket of two conserved quantities is conserved.

**(d) It turned out to be the classical shadow of quantum mechanics.** In September 1925 Paul Dirac read Heisenberg's first paper on the new mechanics, in which the "quantum" variables $x$ and $p$ fail to commute: $xp - px \neq 0$. By his own account, during a Sunday walk he remembered a structure from his classical mechanics course that had the same shape, could not recall its definition, and had to wait for the libraries to open on Monday to check it. It was the Poisson bracket. His paper *The Fundamental Equations of Quantum Mechanics* (November 1925) proposed that every Poisson bracket of classical mechanics becomes a commutator divided by $i\hbar$. Hamilton's equations become Heisenberg's. Every conservation law of this module is a conservation law of QM-4 with hats on. Jacobi's identity, Poisson's theorem, and the angular-momentum algebra of §3.6 are inherited unchanged.

This is the module where mechanics stops being about forces and starts being about structure. The equations of motion do not change. What changes is what you are able to see in them.

## 2. Intuition

### Phase space is the space of states

Newton's second law needs two numbers per degree of freedom to start: a position and a velocity. Hamilton's version uses a position $q$ and a momentum $p$. The pair $(q, p)$ is the complete state: give it to me at one instant and the entire future and past are determined. So the natural place to draw the motion is not the $q$-axis (configuration space), where the particle slides back and forth over the same ground and its trajectory crosses itself, but the $(q, p)$-plane, where each state is one point and each history is one curve.

Here is the pendulum, $H = \frac{p^2}{2ml^2} - mgl\cos\theta$, drawn this way:

```
   p
   ^      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      rotation (E > mgl): the pendulum
   |     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      goes over the top, p never changes sign
   |   __----                                 ----__
   |  /       .-'''''''''-.                         \   separatrix (E = mgl)
   | X       /   .-----.   \       X                 X  = hyperbolic points (theta = +-pi, p = 0)
   |  \      \   '-----'   /       /                 /  libration (E < mgl): closed loops
   |   --____                                 ____--    around the elliptic point (0, 0)
   |     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   |      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   +---------------------------------------------------> theta
          -pi          0           pi
```

(The project draws the real thing as `phase_portraits.png`.) Three features you will meet again and again. A **fixed point** where the vector field vanishes: the bottom of the swing $(0, 0)$ is *elliptic* (nearby orbits circle it), the top $(\pm\pi, 0)$ is *hyperbolic* (nearby orbits approach along one direction and leave along another). A **separatrix**, the curve through the hyperbolic point that divides qualitatively different motions (swinging back and forth versus going over the top). And the fact that **trajectories never cross**: each point has one velocity vector, so one curve through it. That last fact is the uniqueness theorem for ODEs, and it is what makes the picture a *flow*.

### The flow is a fluid, and the fluid is incompressible

Think of the vector field $(\dot q, \dot p)$ as the velocity field of a fluid filling phase space. Put a drop of ink in (a small cloud of nearby initial conditions) and watch it be carried along. It stretches, it folds, it wraps around the elliptic point into a long thin spiral. Here is what you should guess it does not do: change its area.

Why guess that? Try the two systems you can solve in your head. The free particle, $H = p^2/2m$: $\dot q = p/m$, $\dot p = 0$. Points with large $p$ move right faster; the drop is sheared into a parallelogram. A shear preserves area. The harmonic oscillator, $H = \frac{p^2}{2m} + \frac12 m\omega^2 q^2$: in units where the ellipses are circles, the flow is a rigid rotation at angular speed $\omega$. A rotation preserves area. Now look at the general equations: $\dot q = \partial H/\partial p$ and $\dot p = -\partial H/\partial q$. The rate at which the fluid expands is its divergence, $\partial\dot q/\partial q + \partial\dot p/\partial p$. Compute it with Hamilton's equations and watch what happens to the mixed partial derivatives. You should be able to guess the answer before §3.5 states it, and it is the whole of Liouville's theorem.

### Conserved quantities are functions constant along the flow

A conserved quantity $F(q, p)$ is a function whose level sets contain whole trajectories: the flow is tangent to the surface $F = \text{const}$. For one degree of freedom, $H$ itself is conserved and the level sets of $H$ *are* the trajectories (that is why the phase portrait above is a contour plot). How fast does an arbitrary $F$ change along the flow? By the chain rule, $\dot F = \frac{\partial F}{\partial q}\dot q + \frac{\partial F}{\partial p}\dot p = \frac{\partial F}{\partial q}\frac{\partial H}{\partial p} - \frac{\partial F}{\partial p}\frac{\partial H}{\partial q}$. That combination is going to get a name: the Poisson bracket $\{F, H\}$.

### The Poisson bracket is "how fast does $F$ change along the flow generated by $G$"

Nothing in the expression above cared that $H$ was *the* Hamiltonian. Any function $G(q, p)$ defines its own vector field $(\partial G/\partial p, -\partial G/\partial q)$ and its own flow, and $\{F, G\}$ is the rate of change of $F$ along $G$'s flow. Take $G = p$: the flow is $\dot q = 1$, $\dot p = 0$, a uniform translation in $q$, and $\{F, p\} = \partial F/\partial q$, the rate of change of $F$ under translation. Take $G = L_z = xp_y - yp_x$: the flow (you will check this in §3.6) is a rotation about the $z$-axis. Take $G = H$: the flow is time evolution. Momentum *generates* translations, angular momentum *generates* rotations, the Hamiltonian *generates* time evolution.

Now use the one obvious property of the bracket, antisymmetry: $\{F, G\} = -\{G, F\}$. It says: $F$ is constant along $G$'s flow if and only if $G$ is constant along $F$'s flow. With $G = H$: a quantity $F$ is conserved under time evolution if and only if $H$ is unchanged by the flow $F$ generates, that is, if and only if $F$ generates a symmetry of $H$. That is Noether's theorem (CM-4), and in this language it is a one-line consequence of a minus sign. You should guess that every structural fact about symmetry and conservation is going to look like this: short, symmetric, and about brackets.

### A canonical transformation is a change of coordinates that keeps the fluid incompressible

In Lagrangian mechanics you can change coordinates $q \to Q(q)$ freely. In Hamiltonian mechanics there are twice as many coordinates, and a change $(q, p) \to (Q, P)$ that mixes them is allowed too, as long as Hamilton's equations keep their form in the new variables. Which changes are those? The guess from the fluid picture: for one degree of freedom, exactly the area-preserving maps of the plane. Checking whether a map is area-preserving is a $2\times 2$ determinant. For $n$ degrees of freedom the right condition is a little stronger than volume-preserving, and it is the condition that the bracket $\{Q, P\}$, computed in the old coordinates, still equals $1$.

### The Legendre transform describes a curve by its slopes instead of its points

A convex curve $y = f(v)$ can be described by listing its points, or by listing its tangent lines. A tangent line is determined by its slope $p$ and its $y$-intercept, so "the intercept as a function of the slope" is a complete description of the curve. The Legendre transform is that function (with a sign). Passing from $L(q, \dot q)$ to $H(q, p)$ is exactly this: you stop describing the system by a point on the curve $L$-versus-$\dot q$ and start describing it by the slope $p = \partial L/\partial\dot q$. One identity of the Legendre transform is that the derivative of the new function gives back the old variable, and that identity, applied to $H$, is the first of Hamilton's equations: $\partial H/\partial p = \dot q$. The second one is Euler–Lagrange in disguise.

If you can guess, by now, that the divergence of the Hamiltonian vector field is zero, that $\dot F = \{F, H\}$, that $\{q, p\} = 1$ is what canonical transformations preserve, and that $\partial H/\partial p = \dot q$ is a property of the Legendre transform rather than of physics, you have the module. Everything below is making those guesses precise and then using them.
