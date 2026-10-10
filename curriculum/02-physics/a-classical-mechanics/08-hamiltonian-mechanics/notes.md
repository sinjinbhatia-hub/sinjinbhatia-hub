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
in which the coordinates and the momenta appear on exactly the same footing. The letter $H$ is his; the name "Hamiltonian" for the function came later. What Hamilton himself called the *principal function* is a different object, the $S$ of §3.10, and it is the one he cared about most.

Three things happened to this reformulation within a decade, and a fourth took ninety years.

**(a) It made the geometry of motion visible.** A state of the system is now a point $(q, p)$ in a $2n$-dimensional space, and Hamilton's equations are a vector field on that space. The whole history of a system is a curve; the whole set of possible histories is a flow, like a fluid. Questions that are invisible in configuration space (where trajectories cross each other constantly) become obvious in phase space (where they never do).

**(b) The flow preserves a structure.** Joseph Liouville, in an 1838 note on a different subject (the variation of arbitrary constants in a system of ordinary differential equations), proved that for equations of Hamilton's form the flow carries any region of phase space to a region of the same volume. Jacobi noticed what it meant for mechanics. Boltzmann (1868–1871), Maxwell (1879), and finally Gibbs (1902, who named it the "conservation of density-in-phase") turned it into the foundation of statistical mechanics: it is the reason you are allowed to say "all accessible microstates are equally likely" and have the statement stay true as time passes. Liouville's theorem is a corollary of something stronger, which Hamilton's equations also preserve and which was only named in the twentieth century: the *symplectic form*.

**(c) It gave Jacobi a method that solved problems Lagrange could not.** Carl Gustav Jacob Jacobi (1837) simplified Hamilton's optics-inspired construction into a single first-order partial differential equation for a generating function $S(q, t)$, now the Hamilton–Jacobi equation. When that equation can be separated, the problem is solved by quadratures. Jacobi used it to find the geodesics on an ellipsoid (1838) and to handle the motion of a particle attracted by two fixed centers, problems where the Lagrangian equations give you nothing you can integrate. In his Königsberg lectures on dynamics (1842–43) he also proved the identity that now carries his name, and with it the theorem Poisson had stated in 1809: the bracket of two conserved quantities is conserved.

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

## 3. Formal development

Throughout, $n$ is the number of degrees of freedom, $q = (q_1, \dots, q_n)$, $p = (p_1, \dots, p_n)$, and $z = (q, p) \in \mathbb{R}^{2n}$. Partial derivatives of $L$ are taken holding $(q, \dot q, t)$ fixed; partial derivatives of $H$ hold $(q, p, t)$ fixed. Keeping track of *which variables are held fixed* is most of the work in §3.1–3.2, and most of the errors.

### 3.1 The Legendre transform

Let $f(v)$ be a smooth convex function of one real variable ($f'' > 0$).

**Definition 3.1 (Legendre transform).**
$$
g(p) = \max_v\,\big[\,pv - f(v)\,\big]. \tag{3.1}
$$

▶ **Try it first.** Differentiate the bracket with respect to $v$ to find the maximizer, then differentiate $g$ with respect to $p$ and see what comes out.

**Derivation.**

1. The bracket $pv - f(v)$ is concave in $v$ (second derivative $-f'' < 0$), so it has a unique maximum, where its $v$-derivative vanishes:
   $$p = f'(v^*). \tag{3.2}$$
   Since $f'$ is strictly increasing, (3.2) can be solved for $v^*(p)$, and $g(p) = p\,v^*(p) - f(v^*(p))$.
2. Differentiate $g$ by the chain rule:
   $$g'(p) = v^* + p\frac{dv^*}{dp} - f'(v^*)\frac{dv^*}{dp} = v^* + \big(p - f'(v^*)\big)\frac{dv^*}{dp} = v^*. \tag{3.3}$$
   The bracket vanishes by (3.2). The derivative of the new function is the old variable.
3. Differentiate (3.3) once more: $g''(p) = dv^*/dp = 1/f''(v^*) > 0$. So $g$ is convex too, and the transform can be applied to it.

**Theorem 3.1 (Involution).** The Legendre transform of $g$ is $f$.

*Proof.* The transform of $g$ is $\max_p\,[vp - g(p)]$. The maximizer satisfies $v = g'(p^*)$, which by (3.3) says $v = v^*(p^*)$, i.e. $p^* = f'(v)$ by (3.2). The value is $v\,p^* - g(p^*) = v f'(v) - \big[p^* v - f(v)\big] = v f'(v) - v f'(v) + f(v) = f(v)$. ∎

■ **What just happened.** (3.2) and (3.3) are one relation read in two directions: $p = f'(v)$ and $v = g'(p)$. Nothing is lost passing from $f$ to $g$, which is what makes it legitimate to describe a mechanical system by $H(q, p)$ instead of $L(q, \dot q)$.

**The tangent-line picture.** The tangent line to $y = f(v)$ at $v^*$ has slope $p = f'(v^*)$ and $y$-intercept $f(v^*) - pv^* = -g(p)$. So $-g(p)$ is the intercept of the tangent line of slope $p$. A convex curve is the envelope of its tangent lines, so "intercept as a function of slope" determines the curve: that is Theorem 3.1 as a picture. The project draws it (`legendre_tangent_lines.png`): for $f(v) = \tfrac12 v^2$, the tangent of slope $p$ hits the $y$-axis at $-p^2/2$. Until you can see this picture, the Legendre transform is a formula you will keep getting the sign of wrong.

For several variables, $f(v_1, \dots, v_n)$ convex, the same construction gives $g(p) = \max_v\,[\sum_i p_iv_i - f(v)]$ with $p_i = \partial f/\partial v_i$ and $v_i = \partial g/\partial p_i$. You may also transform only some of the variables and carry the others along as parameters; their partial derivatives then flip sign, as (3.6) will show.

**Thermodynamics uses the same transform.** In SM-2 the internal energy $U(S, V)$ has $\partial U/\partial S = T$. The Helmholtz free energy $F(T, V) = U - TS$ is $-g$ in the notation above (thermodynamics writes $U - TS$ rather than $TS - U$, so every sign is reversed), and (3.3) becomes $\partial F/\partial T = -S$. Describing a system by its temperature instead of its entropy is describing the curve $U(S)$ by its slopes. Each thermodynamic potential in SM-2 (enthalpy, Gibbs) is a Legendre transform in one more variable, and each Maxwell relation is equality of mixed partials of one of them.

### 3.2 From $L$ to $H$: Hamilton's equations

Fix $q$ and $t$ and regard $L(q, \dot q, t)$ as a function of the velocities. Assume it is convex in $\dot q$, which holds whenever $T$ is a positive-definite quadratic form in $\dot q$ (CM-3's mass matrix $M_{ij}$ positive definite). Its Legendre transform in the velocities, with $p_i = \partial L/\partial\dot q_i$ as the new variables, is the Hamiltonian:
$$
H(q, p, t) = \sum_{i=1}^n p_i\,\dot q_i - L(q, \dot q, t), \qquad \dot q = \dot q(q, p, t)\ \text{ solved from }\ p_i = \frac{\partial L}{\partial\dot q_i}. \tag{3.4}
$$
As a function of $(q, \dot q, t)$ this is the energy function $h$ of CM-4 (3.8). The only change is the variables it is written in, and that change is everything.

▶ **Try it first.** Write $dH$ two ways: from the left of (3.4) as $\sum_i (\partial H/\partial q_i)\,dq_i + (\partial H/\partial p_i)\,dp_i + (\partial H/\partial t)\,dt$, and from the right by the product rule. Compare coefficients.

**Theorem 3.2 (Hamilton's equations).** Along every solution of the Euler–Lagrange equations,
$$
\dot q_i = \frac{\partial H}{\partial p_i}, \qquad \dot p_i = -\frac{\partial H}{\partial q_i}, \qquad i = 1, \dots, n. \tag{3.5}
$$

**Derivation.**

1. Differential of the right side of (3.4), with $\dot q$ regarded as a function of $(q, p, t)$:
   $$dH = \sum_i \dot q_i\,dp_i + \sum_i p_i\,d\dot q_i - \sum_i\frac{\partial L}{\partial q_i}dq_i - \sum_i\frac{\partial L}{\partial\dot q_i}d\dot q_i - \frac{\partial L}{\partial t}dt.$$
2. The second and fourth sums cancel because $p_i = \partial L/\partial\dot q_i$. This is the step that makes the Legendre transform work: all dependence on $d\dot q$ drops out exactly, so
   $$dH = \sum_i \dot q_i\,dp_i - \sum_i\frac{\partial L}{\partial q_i}dq_i - \frac{\partial L}{\partial t}dt.$$
3. $H$ is a function of $(q, p, t)$, so also $dH = \sum_i\frac{\partial H}{\partial q_i}dq_i + \sum_i\frac{\partial H}{\partial p_i}dp_i + \frac{\partial H}{\partial t}dt$. The differentials $dq_i, dp_i, dt$ are independent, so coefficients match:
   $$\frac{\partial H}{\partial p_i} = \dot q_i, \qquad \frac{\partial H}{\partial q_i} = -\frac{\partial L}{\partial q_i}, \qquad \frac{\partial H}{\partial t} = -\frac{\partial L}{\partial t}. \tag{3.6}$$
4. The first of (3.6) is the first Hamilton equation, and it is an identity: it holds whether or not the motion obeys any law, because it is (3.3) for this transform. The second Hamilton equation needs the equations of motion: Euler–Lagrange says $\dot p_i = \partial L/\partial q_i$, and with the second of (3.6) that is $\dot p_i = -\partial H/\partial q_i$. ∎

■ **What just happened.** The Legendre transform turned $n$ second-order equations into $2n$ first-order ones and gave $q$ and $p$ symmetric roles up to a sign. Half of Hamilton's equations are a property of the transform, not of physics; the other half are Euler–Lagrange. The third identity in (3.6), $\partial H/\partial t = -\partial L/\partial t$, is the sign flip for the untransformed variable promised in §3.1.

**Energy.** Differentiate $H$ along a solution:
$$
\frac{dH}{dt} = \sum_i\left(\frac{\partial H}{\partial q_i}\dot q_i + \frac{\partial H}{\partial p_i}\dot p_i\right) + \frac{\partial H}{\partial t} = \sum_i\left(\frac{\partial H}{\partial q_i}\frac{\partial H}{\partial p_i} - \frac{\partial H}{\partial p_i}\frac{\partial H}{\partial q_i}\right) + \frac{\partial H}{\partial t} = \frac{\partial H}{\partial t}.
$$
So $H$ is conserved iff it has no explicit time dependence, which by (3.6) is iff $L$ has none: CM-4 Theorem 3.6 in new clothes.

**Theorem 3.3 ($H = T + V$).** If the constraints are fixed ($\mathbf{r}_a = \mathbf{r}_a(q)$, no $t$) and $V = V(q, t)$ is velocity-independent, then $H = T + V$.

*Proof.* Under these hypotheses $T = \tfrac12\sum_{ij}M_{ij}(q)\dot q_i\dot q_j$ is homogeneous of degree 2 in $\dot q$ (CM-3 (3.13) with $N_i = T_0 = 0$). Euler's theorem on homogeneous functions (CM-4 (3.10)) gives $\sum_i\dot q_i\,\partial T/\partial\dot q_i = 2T$, so $\sum_i p_i\dot q_i - L = 2T - (T - V) = T + V$. To write it in $p$: $p = M\dot q$, so $\dot q = M^{-1}p$ and $T = \tfrac12 p^\top M^{-1}p$. ∎

When the hypotheses fail (rotating constraints, CM-4 Example 4.3; magnetic forces, §3.4 below), $H$ is still conserved whenever $\partial L/\partial t = 0$, but it is not $T + V$. Do not write $H = T + V$ by reflex. Write $\sum p\dot q - L$ and see what you get.

### 3.3 The phase-space variational principle

Hamilton's principle (CM-2) says the true path makes $\int L\,dt$ stationary among paths $q(t)$ with fixed endpoints. Substitute $L = \sum_i p_i\dot q_i - H$:
$$
S[q, p] = \int_{t_1}^{t_2}\Big(\sum_i p_i\dot q_i - H(q, p, t)\Big)dt. \tag{3.7}
$$

▶ **Try it first.** Vary $q_i(t)$ and $p_i(t)$ *independently*, with $\delta q_i = 0$ at the endpoints and no condition at all on $\delta p_i$. What is needed for $\delta S = 0$?

**Derivation.**

1. $\displaystyle\delta S = \int\sum_i\Big(\delta p_i\,\dot q_i + p_i\,\delta\dot q_i - \frac{\partial H}{\partial q_i}\delta q_i - \frac{\partial H}{\partial p_i}\delta p_i\Big)dt.$
2. Integrate $p_i\,\delta\dot q_i$ by parts: $\int p_i\,\delta\dot q_i\,dt = \big[p_i\,\delta q_i\big]_{t_1}^{t_2} - \int\dot p_i\,\delta q_i\,dt$; the boundary term vanishes because $\delta q_i = 0$ there.
3. Collect: $\displaystyle\delta S = \int\sum_i\Big[\Big(\dot q_i - \frac{\partial H}{\partial p_i}\Big)\delta p_i - \Big(\dot p_i + \frac{\partial H}{\partial q_i}\Big)\delta q_i\Big]dt.$
4. $\delta q_i$ and $\delta p_i$ are independent and arbitrary, so both brackets vanish: Hamilton's equations (3.5). ∎

■ **What just happened.** In (3.7) the momenta are not derived from the velocities; they are independent variables, and $\dot q = \partial H/\partial p$ comes *out* of the variation instead of going in. The integrand contains no $\dot p$, which is why no endpoint condition on $p$ was needed. This is the form of the action that QM-9's path integral is built on: inserting complete sets of position and momentum eigenstates between short time steps produces $\int\mathcal{D}q\,\mathcal{D}p\;\exp\!\big(\tfrac{i}{\hbar}\int(p\dot q - H)\,dt\big)$, and the configuration-space form $\int\mathcal{D}q\,e^{iS[q]/\hbar}$ is what remains after the Gaussian $p$-integrals are done.

### 3.4 Phase space, the flow, and five examples

Write $z = (q_1, \dots, q_n, p_1, \dots, p_n) \in \mathbb{R}^{2n}$ and $\nabla H = (\partial H/\partial q, \partial H/\partial p)$. Hamilton's equations become one vector equation,
$$
\dot z = \mathsf{J}\,\nabla H, \qquad \mathsf{J} = \begin{pmatrix} 0 & I \\ -I & 0\end{pmatrix}, \qquad \mathsf{J}^\top = -\mathsf{J}, \qquad \mathsf{J}^2 = -I. \tag{3.8}
$$
The right side is a vector field on phase space, the **Hamiltonian vector field** of $H$. Its integral curves are the trajectories; the map $\Phi_t$ carrying each $z(0)$ to $z(t)$ is the **flow**. Since $\nabla H\cdot\mathsf{J}\nabla H = 0$ ($\mathsf{J}$ is antisymmetric), the field is tangent to the level sets of $H$. For $n = 1$ the level sets are curves, and the phase portrait is a contour plot of $H$.

**Free particle.** $H = p^2/2m$: $\dot q = p/m$, $\dot p = 0$. Trajectories are horizontal lines traversed at speed $p/m$, and the flow $\Phi_t(q, p) = (q + pt/m,\ p)$ is a shear.

**Harmonic oscillator.** $H = \frac{p^2}{2m} + \frac12 m\omega^2q^2$: $\dot q = p/m$, $\dot p = -m\omega^2q$, so $\ddot q = -\omega^2q$. The level sets are ellipses with semi-axes $\sqrt{2E/m\omega^2}$ along $q$ and $\sqrt{2mE}$ along $p$, of area $\pi\sqrt{2E/m\omega^2}\,\sqrt{2mE} = 2\pi E/\omega$, a number you will meet again in (3.25). In the scaled coordinates $(\sqrt{m\omega}\,q,\ p/\sqrt{m\omega})$ the ellipses are circles and the flow is a rigid clockwise rotation at angular speed $\omega$.

**Central force, reduced.** A particle in a plane with potential $V(r)$: $L = \tfrac12 m(\dot r^2 + r^2\dot\phi^2) - V(r)$, $p_r = m\dot r$, $p_\phi = mr^2\dot\phi$, and by Theorem 3.3
$$
H = \frac{p_r^2}{2m} + \frac{p_\phi^2}{2mr^2} + V(r).
$$
$\phi$ is cyclic, so $\dot p_\phi = -\partial H/\partial\phi = 0$ and $p_\phi = \ell$ is constant. Freeze it: $H$ is a one-degree-of-freedom Hamiltonian in $(r, p_r)$ with effective potential $V_{\text{eff}} = V + \ell^2/2mr^2$, the CM-5 reduction. The angle is recovered afterwards from $\dot\phi = \partial H/\partial p_\phi = \ell/mr^2$. In the Hamiltonian picture every cyclic coordinate drops out this cleanly; in the Lagrangian picture $\dot\phi$ still sits inside $L$ and substituting $\ell$ for it before varying gives the wrong equations (Goldstein's Routhian, §8.3, is the fix, and it is a partial Legendre transform). Example 4.2 does this for the spherical pendulum.

**Pendulum.** $H = \frac{p^2}{2ml^2} - mgl\cos\theta$, $\theta$ from the downward vertical. Fixed points, where $\nabla H = 0$: $(0, 0)$, a minimum of $H$ (elliptic), and $(\pm\pi, 0)$, saddles (hyperbolic). The separatrix is the level set through the saddle, $E_s = H(\pi, 0) = mgl$. Below it, libration (closed curves around the origin); above it, rotation ($p$ never changes sign). Example 4.3 computes the separatrix curve, the area inside it, and the motion along it; Exercise C1 computes the period, which diverges logarithmically as $E \to E_s$.

**Charged particle: canonical versus kinetic momentum.** For charge $q_e$ in potentials $\phi(\mathbf{r}, t)$ and $\mathbf{A}(\mathbf{r}, t)$ the Lagrangian is $L = \tfrac12 m|\dot{\mathbf{r}}|^2 + q_e\,\dot{\mathbf{r}}\cdot\mathbf{A} - q_e\phi$ (CM-3 §3.5). The canonical momentum is
$$
\mathbf{p} = \frac{\partial L}{\partial\dot{\mathbf{r}}} = m\dot{\mathbf{r}} + q_e\mathbf{A},
$$
and it is not $m\mathbf{v}$. Then $\dot{\mathbf{r}} = (\mathbf{p} - q_e\mathbf{A})/m$ and
$$
H = \mathbf{p}\cdot\dot{\mathbf{r}} - L = \dot{\mathbf{r}}\cdot(\mathbf{p} - q_e\mathbf{A}) - \tfrac12 m|\dot{\mathbf{r}}|^2 + q_e\phi = \tfrac12 m|\dot{\mathbf{r}}|^2 + q_e\phi = \frac{|\mathbf{p} - q_e\mathbf{A}|^2}{2m} + q_e\phi. \tag{3.9}
$$
The value is kinetic plus electrostatic energy, but as a function on phase space it is not $p^2/2m + V$: the magnetic field enters only through the shift $\mathbf{p} \to \mathbf{p} - q_e\mathbf{A}$. This is the form that survives into QM-4, where $\hat{\mathbf{p}} = -i\hbar\nabla$ is the *canonical* momentum and a magnetic field enters the Schrödinger equation as $(-i\hbar\nabla - q_e\mathbf{A})^2/2m$. Exercise B1 checks that Hamilton's equations for (3.9) give the Lorentz force.

### 3.5 Liouville's theorem

**Theorem 3.4 (Liouville).** The Hamiltonian vector field is divergence-free, and the flow preserves phase-space volume: for any region $D \subset \mathbb{R}^{2n}$, $\operatorname{vol}(\Phi_t(D)) = \operatorname{vol}(D)$.

▶ **Try it first.** Compute $\sum_i\big(\partial\dot q_i/\partial q_i + \partial\dot p_i/\partial p_i\big)$ from (3.5).

**Derivation.**

1. Divergence of the field $(\dot q, \dot p) = (\partial H/\partial p, -\partial H/\partial q)$:
   $$\nabla\cdot\dot z = \sum_{i=1}^n\left(\frac{\partial}{\partial q_i}\frac{\partial H}{\partial p_i} - \frac{\partial}{\partial p_i}\frac{\partial H}{\partial q_i}\right) = \sum_i\left(\frac{\partial^2H}{\partial q_i\,\partial p_i} - \frac{\partial^2H}{\partial p_i\,\partial q_i}\right) = 0. \tag{3.10}$$
   Equality of mixed partials, and nothing else.
2. Let $D_t = \Phi_t(D)$. Its volume changes at the rate at which the boundary moves outward: $\frac{d}{dt}\operatorname{vol}(D_t) = \oint_{\partial D_t}\dot z\cdot\hat{\mathbf{n}}\,dA$ (each boundary patch sweeps out volume at the rate normal velocity times area). By the divergence theorem this equals $\int_{D_t}\nabla\cdot\dot z\,dV$, which is zero by step 1. ∎

Step 2 is the transport argument that gives the continuity equation in fluid mechanics; here the fluid has velocity field $\mathsf{J}\nabla H$ and is incompressible. The project transports a disk of 5000 points under the pendulum flow and watches it stretch into a spiral arm while its area stays put (`liouville_blob.png`). Note what is *not* claimed: the shape is not preserved, and the diameter of the blob grows without bound. Only the volume is fixed.

**The Liouville equation.** Put a density on phase space: $\rho(z, t)\,d^{2n}z$ is the probability of finding the system in the volume element $d^{2n}z$. Conservation of probability is the continuity equation $\partial\rho/\partial t + \nabla\cdot(\rho\dot z) = 0$. Expand: $\nabla\cdot(\rho\dot z) = \dot z\cdot\nabla\rho + \rho\,\nabla\cdot\dot z = \dot z\cdot\nabla\rho$ by (3.10). So
$$
\frac{\partial\rho}{\partial t} + \sum_i\left(\frac{\partial\rho}{\partial q_i}\frac{\partial H}{\partial p_i} - \frac{\partial\rho}{\partial p_i}\frac{\partial H}{\partial q_i}\right) = 0, \qquad\text{i.e.}\qquad \frac{\partial\rho}{\partial t} + \{\rho, H\} = 0 \tag{3.11}
$$
in the bracket notation of (3.12) below. Equivalently $d\rho/dt = 0$ along trajectories: the density is carried along the flow unchanged, like dye in an incompressible fluid.

**Why statistical mechanics is possible.** A stationary ensemble has $\partial\rho/\partial t = 0$, i.e. $\{\rho, H\} = 0$. Any $\rho = f(H)$ qualifies, because $\{f(H), H\} = f'(H)\{H, H\} = 0$. The microcanonical ensemble of SM-3, $\rho \propto \delta(H - E)$, uniform on the energy surface, is such an $f$; so is SM-4's canonical $\rho \propto e^{-\beta H}$. Liouville's theorem is what makes "all accessible microstates are equally likely" a *consistent* assumption: if it is true at $t = 0$ it stays true. Without volume preservation a uniform density would bunch up as it flowed and the assumption would contradict itself. Nothing here says the uniform density is *reached* from a non-uniform start; that is the ergodic problem, which is hard, and SM-3 is honest about it.

**Poincaré recurrence (stated).** If the accessible region of phase space has finite volume (bounded motion at fixed energy), then for almost every initial point and every neighbourhood $U$ of it, the trajectory returns to $U$ infinitely often. Poincaré's proof (1890) uses only Liouville: the images $\Phi_{k\tau}(U)$, $k = 1, 2, \dots$, all have the same volume, so they cannot all be disjoint inside a region of finite volume; two overlap, and running the flow backward shows that $U$ meets a later image of itself. Zermelo used it in 1896 against Boltzmann's $H$-theorem. The resolution (recurrence times for $10^{23}$ particles that dwarf the age of the universe) is SM-3's business.
