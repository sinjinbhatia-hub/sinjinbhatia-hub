# CM-8: Hamiltonian Mechanics: Legendre Transform, Phase Space, Poisson Brackets, Liouville, Canonical Transformations, Hamilton–Jacobi

**Time:** ~12 hours · **Prerequisites:** CM-4 (generalized momenta, cyclic coordinates, the energy function) · **Unlocks:** SM-3 (Liouville and ensembles), QM-4 (commutators from Poisson brackets), QM-9 (path integrals, the classical limit), QC-9 (Hamiltonian simulation)

> **Conventions in this module.** Generalized coordinates are $q = (q_1, \dots, q_n)$ and their conjugate momenta are $p_i = \partial L/\partial \dot q_i$, exactly as in CM-4. The Hamiltonian is $H(q, p, t) = \sum_i p_i \dot q_i - L$. Hamilton's equations are $\dot q_i = \partial H/\partial p_i$, $\dot p_i = -\partial H/\partial q_i$. The Poisson bracket is $\{F, G\} = \sum_i \left(\frac{\partial F}{\partial q_i}\frac{\partial G}{\partial p_i} - \frac{\partial F}{\partial p_i}\frac{\partial G}{\partial q_i}\right)$, so that $\{q, p\} = +1$ (some books, including older editions of Landau, use the opposite sign; check before copying a formula). Phase space is the set of points $z = (q, p) \in \mathbb{R}^{2n}$. The symplectic matrix is written in sans-serif, $\mathsf{J} = \begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}$, because the italic $J = \oint p\,dq$ is reserved for the action variable; both names are standard and I will not rename either. The angle variable $w$ conjugate to $J$ advances by exactly $1$ per period, so the frequency $\nu = \partial H/\partial J$ is in cycles per unit time and $\omega = 2\pi\nu$. $L$ with no index is the Lagrangian; $L_x, L_y, L_z$ are components of angular momentum and $L^2 = L_x^2 + L_y^2 + L_z^2$. $\epsilon_{ijk}$ is the Levi-Civita symbol with $\epsilon_{xyz} = 1$. Sums are written out; no summation convention. The pendulum's angle $\theta$ is measured from the downward vertical, so $V = -mgl\cos\theta$ and the hanging rest point is $\theta = 0$. $\hbar = h/2\pi$.

## 0. Where this sits

You can write a Lagrangian, get Euler–Lagrange equations out of it, read off conserved momenta from cyclic coordinates, and you know the energy function $h = \sum_i \dot q_i\,\partial L/\partial\dot q_i - L$ from CM-4. That is everything this module needs from mechanics. From mathematics it needs multivariable calculus (the chain rule for partial derivatives, equality of mixed partials, the divergence theorem) and the ability to compute a $2\times 2$ determinant.

What it makes possible: SM-3's ensembles rest on Liouville's theorem (§3.5); QM-4's $[\hat x, \hat p] = i\hbar$ and Heisenberg equation are the Poisson brackets of §3.6; QM-8's WKB is the Hamilton–Jacobi equation of §3.10; QM-9's path integral starts from the action of §3.3; QC-9's Trotterization and ST-6's Hamiltonian Monte Carlo are the integrators of §3.11.

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

**Harmonic oscillator.** $H = \frac{p^2}{2m} + \frac12 m\omega^2q^2$: $\dot q = p/m$, $\dot p = -m\omega^2q$, so $\ddot q = -\omega^2q$. The level sets are ellipses with semi-axes $\sqrt{2E/m\omega^2}$ along $q$ and $\sqrt{2mE}$ along $p$, of area $\pi\sqrt{2E/m\omega^2}\,\sqrt{2mE} = 2\pi E/\omega$, a number you will meet again in (3.26). In the scaled coordinates $(\sqrt{m\omega}\,q,\ p/\sqrt{m\omega})$ the ellipses are circles and the flow is a rigid clockwise rotation at angular speed $\omega$.

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

### 3.6 Poisson brackets

**Definition 3.2 (Poisson bracket).** For functions $F(q, p, t)$ and $G(q, p, t)$,
$$
\{F, G\} = \sum_{i=1}^n\left(\frac{\partial F}{\partial q_i}\frac{\partial G}{\partial p_i} - \frac{\partial F}{\partial p_i}\frac{\partial G}{\partial q_i}\right), \qquad \{q, p\} = +1. \tag{3.12}
$$
Properties, for constants $a, b$: (i) antisymmetry $\{F, G\} = -\{G, F\}$; (ii) bilinearity $\{aF + bG, K\} = a\{F, K\} + b\{G, K\}$; (iii) Leibniz $\{FG, K\} = F\{G, K\} + \{F, K\}G$; (iv) Jacobi
$$
\{F, \{G, K\}\} + \{G, \{K, F\}\} + \{K, \{F, G\}\} = 0. \tag{3.13}
$$
(i) and (ii) are read off (3.12); (iii) is the product rule.

▶ **Try it first.** In the expansion of (3.13), each term is a product of three factors. How many of them are second derivatives?

**Proof of (3.13)** (Landau §42, any $n$).

1. Each term of $\{A, \{B, C\}\}$ is a first derivative of $A$ times a derivative of (first derivative of $B$ × first derivative of $C$). So every term of (3.13) contains exactly one second derivative, of $F$, $G$ or $K$.
2. Second derivatives of $F$ come only from the last two brackets. Write $D_G = \{\,\cdot\,, G\} = \sum_i\big(\frac{\partial G}{\partial p_i}\partial_{q_i} - \frac{\partial G}{\partial q_i}\partial_{p_i}\big)$, a first-order differential operator. By (i), $\{G, \{K, F\}\} = D_GD_KF$ and $\{K, \{F, G\}\} = -D_KD_GF$, so these terms are $(D_GD_K - D_KD_G)F$.
3. For first-order operators $X = \sum_\alpha a_\alpha\partial_\alpha$, $Y = \sum_\beta b_\beta\partial_\beta$, the second-derivative part of $XYF$ is $\sum_{\alpha\beta}a_\alpha b_\beta\,\partial_\alpha\partial_\beta F$, symmetric under $X \leftrightarrow Y$ by equality of mixed partials. It cancels in $XY - YX$. So (3.13) has no second derivatives of $F$.
4. (3.13) is invariant under $F \to G \to K \to F$, so it has none of $G$ or $K$ either. By step 1, every term cancels. ∎

**Theorem 3.5 (Time evolution).** Along a solution of (3.5), for any $F(q, p, t)$,
$$
\frac{dF}{dt} = \{F, H\} + \frac{\partial F}{\partial t}. \tag{3.14}
$$
*Proof.* Chain rule, $\dot F = \sum_i\big(\frac{\partial F}{\partial q_i}\dot q_i + \frac{\partial F}{\partial p_i}\dot p_i\big) + \frac{\partial F}{\partial t}$, then substitute (3.5). ∎

So $\dot q_i = \{q_i, H\}$, $\dot p_i = \{p_i, H\}$, and $F$ with no explicit $t$ is conserved iff $\{F, H\} = 0$. The **fundamental brackets** are
$$
\{q_i, q_j\} = 0, \qquad \{p_i, p_j\} = 0, \qquad \{q_i, p_j\} = \delta_{ij}, \tag{3.15}
$$
and with (i)–(iii) they compute the bracket of any two polynomials without differentiating.

**Theorem 3.6 (Poisson).** If $F$ and $G$ have no explicit $t$ and are conserved, so is $\{F, G\}$.

*Proof.* Jacobi with $K = H$: $\{H, \{F, G\}\} = -\{F, \{G, H\}\} - \{G, \{H, F\}\} = 0$. ∎

**Angular momentum.** For $\mathbf{L} = \mathbf{r}\times\mathbf{p}$,
$$
L_x = yp_z - zp_y, \qquad L_y = zp_x - xp_z, \qquad L_z = xp_y - yp_x, \tag{3.16}
$$
$$
\{L_i, L_j\} = \sum_k\epsilon_{ijk}L_k, \qquad \{L^2, L_i\} = 0. \tag{3.17}
$$
Example 4.4 derives $\{L_x, L_y\} = L_z$; cyclic relabelling $x \to y \to z$ maps (3.16) to itself and gives the rest. By Leibniz, $\{L^2, L_z\} = 2L_x\{L_x, L_z\} + 2L_y\{L_y, L_z\} = -2L_xL_y + 2L_yL_x = 0$. By Poisson's theorem, conserving $L_x$ and $L_y$ conserves $L_z$.

**Generators.** Any $G(q, p)$ defines a flow $dz/d\lambda = \mathsf{J}\nabla G$ along which $dF/d\lambda = \{F, G\}$. So
$$
\{F, p_j\} = \frac{\partial F}{\partial q_j}, \qquad \{x, L_z\} = -y, \quad \{y, L_z\} = x, \qquad \{F, H\} = \dot F: \tag{3.18}
$$
$p_j$ generates translation in $q_j$, $L_z$ generates counterclockwise rotation about $z$, $H$ generates time evolution. By antisymmetry, $\{F, H\} = 0$ says both "$F$ is conserved" and "$H$ is invariant under $F$'s flow": Noether's theorem (CM-4) without a Lagrangian.

**The quantum dictionary.** Dirac's 1925 rule replaces phase-space functions by operators and
$$
\{A, B\} \;\longrightarrow\; \frac{1}{i\hbar}[\hat A, \hat B], \qquad [\hat A, \hat B] = \hat A\hat B - \hat B\hat A. \tag{3.19}
$$

| Classical | Quantum (QM-4) |
|---|---|
| $\{q, p\} = 1$ | $[\hat x, \hat p] = i\hbar$ |
| $\dot F = \{F, H\}$ | $d\hat A/dt = \frac{i}{\hbar}[\hat H, \hat A]$ (Heisenberg) |
| $\{L_x, L_y\} = L_z$ | $[\hat L_x, \hat L_y] = i\hbar\hat L_z$ |
| conserved iff $\{F, H\} = 0$ | conserved iff $[\hat A, \hat H] = 0$ |
| $p$ generates translations | $e^{-ia\hat p/\hbar}$ translates by $a$ |

It is consistent because the commutator also obeys (i)–(iv), with Leibniz in the order-sensitive form $[\hat A\hat B, \hat C] = \hat A[\hat B, \hat C] + [\hat A, \hat C]\hat B$ (which is why (iii) was written that way). Every theorem above that used only (i)–(iv) holds in quantum mechanics with no new proof. What fails is applying (3.19) to *every* function: $qp$ has no unique operator, and Groenewold (1946) proved no rule satisfies (3.19) for all polynomials.

### 3.7 Canonical transformations

**Definition 3.3.** An invertible, time-independent map $z = (q, p) \mapsto Z = (Q, P)$ is **canonical** if, for every $H$, the motion obeys $\dot Q_i = \partial K/\partial P_i$, $\dot P_i = -\partial K/\partial Q_i$ with $K(Z) = H(z(Z))$. Let $M = \partial Z/\partial z$ be its $2n\times 2n$ Jacobian.

▶ **Try it first.** Use $\dot Z = M\dot z$, $\nabla_zH = M^\top\nabla_ZK$ and $\dot z = \mathsf{J}\nabla_zH$. What must $M$ satisfy?

**Theorem 3.7.** Equivalent: (a) the map is canonical; (b) $(Q, P)$ satisfy (3.15) when the brackets are computed in $(q, p)$; (c)
$$
M^\top\mathsf{J}M = \mathsf{J}. \tag{3.20}
$$

**Derivation.**

1. $\dot Z = M\mathsf{J}\nabla_zH = M\mathsf{J}M^\top\nabla_ZK$. This equals $\mathsf{J}\nabla_ZK$ for every $K$ iff $M\mathsf{J}M^\top = \mathsf{J}$. So (a) ⇔ $M\mathsf{J}M^\top = \mathsf{J}$.
2. By (3.12), $\{Z_a, Z_b\} = (\nabla_zZ_a)^\top\mathsf{J}\nabla_zZ_b = (M\mathsf{J}M^\top)_{ab}$, and (3.15) says this matrix is $\mathsf{J}$. So (b) ⇔ $M\mathsf{J}M^\top = \mathsf{J}$.
3. Invert $M\mathsf{J}M^\top = \mathsf{J}$ with $\mathsf{J}^{-1} = -\mathsf{J}$: $M^{-\top}\mathsf{J}M^{-1} = \mathsf{J}$; multiply by $M^\top$ on the left and $M$ on the right to get (3.20). Reversible. ∎

**One degree of freedom.** For $M = \begin{pmatrix} a & b\\ c & d\end{pmatrix}$, multiplying out gives $M^\top\mathsf{J}M = (ad - bc)\mathsf{J}$, and $ad - bc = \frac{\partial Q}{\partial q}\frac{\partial P}{\partial p} - \frac{\partial Q}{\partial p}\frac{\partial P}{\partial q} = \{Q, P\}$. So canonical ⇔ $\{Q, P\} = 1$ ⇔ $\det M = 1$ ⇔ area-preserving, the guess of §2. For $n > 1$, (3.20) implies $|\det M| = 1$ but is stronger: it preserves the sum of oriented areas projected onto the $(q_i, p_i)$ planes (the symplectic form).

**The flow is canonical.** Let $M(t) = \partial\Phi_t/\partial z$ and $S$ the Hessian of $H$. Differentiating $\dot z = \mathsf{J}\nabla H$ gives $\dot M = \mathsf{J}SM$, $M(0) = I$, so
$$
\tfrac{d}{dt}\big(M^\top\mathsf{J}M\big) = M^\top S\mathsf{J}^\top\mathsf{J}M + M^\top\mathsf{J}^2SM = M^\top SM - M^\top SM = 0.
$$
This is the structure stronger than volume promised in §1(b); Liouville is its determinant.

**Examples.** *Exchange* $Q = p$, $P = -q$: $\det M = 1$ ($Q = p$, $P = +q$ has $\det = -1$ and fails). *Point transformations* $Q = f(q)$, $P = p/f'(q)$: $\det M = 1$. *Oscillator polar coordinates*:
$$
q = \sqrt{\frac{2P}{m\omega}}\,\sin Q, \qquad p = \sqrt{2Pm\omega}\,\cos Q, \qquad H = \frac{2Pm\omega\cos^2Q}{2m} + \frac{m\omega^2}{2}\frac{2P}{m\omega}\sin^2Q = \omega P. \tag{3.21}
$$
Once Example 4.5 shows (3.21) is canonical, $K = \omega P$ has $Q$ cyclic: $P$ constant, $\dot Q = \omega$. The oscillator is solved by a change of variables.

### 3.8 Generating functions

To *build* canonical maps, use (3.7): two action integrands that differ by a total derivative give the same variational equations. So a map is canonical if
$$
\sum_ip_i\,dq_i - H\,dt = \sum_iP_i\,dQ_i - K\,dt + dF, \tag{3.22}
$$
where $F$ depends on $t$ and on one old and one new variable from each pair. For $F_1(q, Q, t)$, matching $dq$, $dQ$, $dt$ gives $p = \partial F_1/\partial q$, $P = -\partial F_1/\partial Q$, $K = H + \partial F_1/\partial t$.

▶ **Try it first.** Put $F = F_2(q, P, t) - \sum_iQ_iP_i$ into (3.22) and match coefficients.

**Derivation (type 2).**

1. $dF = \sum_i\frac{\partial F_2}{\partial q_i}dq_i + \sum_i\frac{\partial F_2}{\partial P_i}dP_i + \frac{\partial F_2}{\partial t}dt - \sum_iP_i\,dQ_i - \sum_iQ_i\,dP_i$.
2. The right side of (3.22) becomes $\sum_i\frac{\partial F_2}{\partial q_i}dq_i + \sum_i\big(\frac{\partial F_2}{\partial P_i} - Q_i\big)dP_i + \big(\frac{\partial F_2}{\partial t} - K\big)dt$; the $P\,dQ$ terms cancel.
3. The left side, in the independent variables $(q, P, t)$, is $\sum_ip_i\,dq_i + 0\cdot dP - H\,dt$. Matching:
$$
p_i = \frac{\partial F_2}{\partial q_i}, \qquad Q_i = \frac{\partial F_2}{\partial P_i}, \qquad K = H + \frac{\partial F_2}{\partial t}. \tag{3.23}
$$
∎

■ **What just happened.** Any function of one old and one new variable per pair generates a canonical map through its partial derivatives. Subtracting $QP$ is a Legendre transform from $F_1$ to $F_2$.

| Type | Variables | Old | New |
|---|---|---|---|
| $F_1$ | $q, Q$ | $p = \partial F_1/\partial q$ | $P = -\partial F_1/\partial Q$ |
| $F_2$ | $q, P$ | $p = \partial F_2/\partial q$ | $Q = \partial F_2/\partial P$ |
| $F_3$ | $p, Q$ | $q = -\partial F_3/\partial p$ | $P = -\partial F_3/\partial Q$ |
| $F_4$ | $p, P$ | $q = -\partial F_4/\partial p$ | $Q = \partial F_4/\partial P$ |

Always $K = H + \partial F/\partial t$. *Identity:* $F_2 = \sum_iq_iP_i$. *Point transformations:* $F_2 = \sum_if_i(q)P_i$ gives $Q_i = f_i(q)$ and $p_j = \sum_iP_i\,\partial f_i/\partial q_j$. *Exchange:* $F_1 = qQ$. Exercise B5 shows $F_1 = \frac12m\omega q^2\cot Q$ generates (3.21).

### 3.9 Action-angle variables and adiabatic invariants

One degree of freedom, $H(q, p)$ with no explicit $t$, periodic orbit at energy $E$. The **action** is
$$
J(E) = \oint p\,dq, \tag{3.24}
$$
over one cycle; for a libration it is the enclosed phase-plane area. We want canonical $(w, J)$ with $K = K(J)$.

▶ **Try it first.** If $K = K(J)$, what are $\dot J$ and $\dot w$? How far must $w$ advance per period?

**Derivation.**

1. Type-2 function $W(q, J)$: $p = \partial W/\partial q$, $w = \partial W/\partial J$, $K = H$. On the orbit, solve $H(q, p) = E(J)$ for $p(q, E)$; then $W = \int^q p(q', E(J))\,dq'$.
2. $\dot J = -\partial K/\partial w = 0$ and $\dot w = dK/dJ \equiv \nu$, constant: $w = \nu t + w_0$.
3. Over one cycle, $\Delta w = \oint\frac{\partial^2W}{\partial J\,\partial q}dq = \frac{d}{dJ}\oint p\,dq = 1$.
4. One cycle takes the period $T$, so $\nu T = 1$:
$$
\nu = \frac{dE}{dJ} = \frac{1}{T}, \qquad T = \frac{dJ}{dE}. \tag{3.25}
$$
∎

Check without canonical machinery: differentiating $H(q, p(q, E)) = E$ gives $\partial p/\partial E = 1/\dot q$, so $dJ/dE = \oint dq/\dot q = \oint dt = T$.

■ **What just happened.** The orbit's area as a function of energy knows the period. In action-angle variables every bounded 1-DOF motion is uniform rotation of $w$.

**Oscillator.** The §3.4 ellipse gives
$$
J = \frac{2\pi E}{\omega}, \qquad \nu = \frac{dE}{dJ} = \frac{\omega}{2\pi}, \tag{3.26}
$$
independent of amplitude. In (3.21), $P = J/2\pi$ and $Q = 2\pi w$.

**Pendulum (stated; derived in Exercise C1).** For libration with $k = \sin(\theta_0/2)$, $k^2 = (1 + E/mgl)/2$: $J = 16\,ml\sqrt{gl}\,[\mathrm{E}(k) - (1 - k^2)\mathrm{K}(k)]$ and $T = 4\sqrt{l/g}\,\mathrm{K}(k)$, with $\mathrm{K}(k) = \int_0^{\pi/2}(1 - k^2\sin^2\varphi)^{-1/2}d\varphi$ and $\mathrm{E}(k) = \int_0^{\pi/2}(1 - k^2\sin^2\varphi)^{1/2}d\varphi$. At $E = 0$ ($m = l = g = 1$) the project gets $J = 6.77770$, $T = 7.41630$, and $dJ/dE = T$ to $10^{-6}$. As $E \to mgl$, $J \to 16\,ml\sqrt{gl}$ (Example 4.3) and $T \to \infty$.

**Theorem 3.8 (Adiabatic invariance, stated).** If $H(q, p; \lambda)$ depends on a parameter that changes slowly ($T|\dot\lambda| \ll |\lambda|$), the total change in $J$ tends to zero as the rate tends to zero, while $E$ changes by a finite amount. (Landau §49–51 proves it by averaging over the angle.)

We derive the oscillator case, $H = \frac{p^2}{2m} + \frac12m\omega(t)^2q^2$.

1. By §3.2, $\dot E = \partial H/\partial t = m\omega\dot\omega q^2$.
2. Average over one period with $\omega$ frozen (the slowness assumption). The virial theorem gives $\langle\frac12m\omega^2q^2\rangle = E/2$, so $\langle q^2\rangle = E/m\omega^2$.
3. $\langle\dot E\rangle = (\dot\omega/\omega)E$, i.e. $d\ln E = d\ln\omega$:
$$
\frac{E}{\omega} = \text{const}, \qquad J = \frac{2\pi E}{\omega} = \text{const}. \tag{3.27}
$$
∎

At the 1911 Solvay conference Lorentz asked whether a pendulum whose string is slowly shortened stays in a quantum state; Einstein answered that $E/\nu$ is unchanged. Ehrenfest (1916) made it a principle: only adiabatic invariants may be quantized. Hence the old quantum theory's **Bohr–Sommerfeld rule** $\oint p\,dq = nh$ (QM-1). For the oscillator it gives $E = n\hbar\omega$; the correct $(n + \frac12)\hbar\omega$ needs QM-8's turning-point correction.

### 3.10 The Hamilton–Jacobi equation

The best canonical map makes $K = 0$: every new variable is then constant.

▶ **Try it first.** Take a type-2 $S(q, P, t)$ in (3.23) and demand $K = 0$. What equation must $S$ satisfy?

**Derivation.**

1. By (3.23): $p_i = \partial S/\partial q_i$, $Q_i = \partial S/\partial P_i$, $K = H + \partial S/\partial t$.
2. $K = 0$ gives $\dot Q = \dot P = 0$: $P_i = \alpha_i$, $Q_i = \beta_i$, constants.
3. Substituting $p = \partial S/\partial q$ into $K = 0$:
$$
H\Big(q, \frac{\partial S}{\partial q}, t\Big) + \frac{\partial S}{\partial t} = 0. \tag{3.28}
$$
4. Given a **complete integral** $S(q, \alpha, t)$ with $n$ independent constants, solve $\beta_i = \partial S/\partial\alpha_i$ for $q(t; \alpha, \beta)$; then $p = \partial S/\partial q$. ∎

**Separation.** If $H$ has no explicit $t$, $S = W(q, \alpha) - Et$ turns (3.28) into
$$
H\Big(q, \frac{\partial W}{\partial q}\Big) = E, \tag{3.29}
$$
and $W$ is the function of §3.9. If $H$ is a sum of pieces each involving one $q_i$, $W = \sum_iW_i(q_i)$ reduces it to $n$ single integrals.

**$S$ is the action.** Along a solution, $\frac{dS}{dt} = \sum_i\frac{\partial S}{\partial q_i}\dot q_i + \frac{\partial S}{\partial t} = \sum_ip_i\dot q_i - H = L$. So $S$ is CM-2's action $\int L\,dt$, viewed as a function of the endpoint, and $\partial S/\partial q = p$: momentum is the gradient of the action. For $H = p^2/2m + V$ the velocity $\nabla S/m$ is normal to the surfaces $S = \text{const}$: trajectories are rays of a wavefront, Hamilton's optical picture. Example 4.6 solves the oscillator this way.

**The capstone: Schrödinger's equation contains (3.28).**

▶ **Try it first.** Put $\psi = A\,e^{iS/\hbar}$ ($A$, $S$ real) into $i\hbar\psi_t = -\frac{\hbar^2}{2m}\psi_{qq} + V\psi$ and separate real and imaginary parts.

**Derivation** (subscripts are partial derivatives).

1. $\psi_t = \big(A_t + \frac{i}{\hbar}AS_t\big)e^{iS/\hbar}$, $\psi_{qq} = \big(A_{qq} + \frac{2i}{\hbar}A_qS_q + \frac{i}{\hbar}AS_{qq} - \frac{1}{\hbar^2}AS_q^2\big)e^{iS/\hbar}$.
2. Substituting and dividing by $e^{iS/\hbar}$: $i\hbar A_t - AS_t = -\frac{\hbar^2}{2m}A_{qq} - \frac{i\hbar}{2m}(2A_qS_q + AS_{qq}) + \frac{1}{2m}AS_q^2 + VA$.
3. Real part, divided by $-A$: $S_t + \frac{S_q^2}{2m} + V = \frac{\hbar^2}{2m}\frac{A_{qq}}{A}$.
4. Imaginary part, times $2A/\hbar$: $\partial_t(A^2) + \partial_q(A^2S_q/m) = 0$.
5. When $\hbar^2|A_{qq}/A| \ll S_q^2 = p^2$ (wavelength short compared with the scale on which $A$ varies), step 3 is
$$
\frac{\partial S}{\partial t} + \frac{1}{2m}\Big(\frac{\partial S}{\partial q}\Big)^2 + V = 0, \tag{3.30}
$$
which is (3.28) for $H = p^2/2m + V$. ∎

■ **What just happened.** To leading order the phase of the wavefunction, in units of $\hbar$, is Hamilton's principal function. Step 4 is a continuity equation for $\rho = A^2$ moving at $p/m$; for a stationary state it gives $A^2p = \text{const}$, so $|\psi|^2 \propto 1/p(q)$, the classical fraction of time spent near $q$. The project compares this with the exact $n = 20$ oscillator state: the probability of $|q| < A/2$ ($A$ the classical amplitude) is $\frac{2}{\pi}\arcsin\frac12 = \frac13$ classically and $0.3435$ exactly. The approximation fails at turning points, where $p \to 0$; QM-8 repairs it.

### 3.11 Symplectic integrators, integrability, chaos

RK4 is accurate per step but its step map is not canonical, so the energy error accumulates over many periods. A **symplectic integrator** is a canonical map; backward error analysis (Hairer, Lubich & Wanner, *Geometric Numerical Integration*, ch. IX) shows it is, up to exponentially small terms, the exact flow of a modified Hamiltonian $\tilde H = H + O(h^r)$, $r$ the order. It conserves $\tilde H$, so $H$ oscillates within $O(h^r)$ forever.

For $H = T(p) + V(q)$, **symplectic Euler** is kick then drift:
$$
p_{n+1} = p_n - h\,V'(q_n), \qquad q_{n+1} = q_n + h\,T'(p_{n+1}). \tag{3.31}
$$
The kick is the exact time-$h$ flow of $V$ alone, Jacobian $\begin{pmatrix}1 & 0\\ -hV'' & 1\end{pmatrix}$; the drift is the exact flow of $T$ alone, $\begin{pmatrix}1 & hT''\\ 0 & 1\end{pmatrix}$. Each is a shear of determinant 1, each is canonical because flows are, and products of matrices obeying (3.20) obey it. For the oscillator the composite exactly conserves $\tilde H = H - \frac{h}{2}\omega^2qp$ (checked symbolically in the project). Half kick, drift, half kick is **Störmer–Verlet**, second order and time-reversible:
$$
p_{n+\frac12} = p_n - \tfrac{h}{2}V'(q_n), \qquad q_{n+1} = q_n + h\,T'(p_{n+\frac12}), \qquad p_{n+1} = p_{n+\frac12} - \tfrac{h}{2}V'(q_{n+1}). \tag{3.32}
$$
Pendulum from $\theta_0 = 2$ rad, 50 steps per period, $10^4$ periods: RK4's relative energy error grows roughly linearly to $0.257$; Verlet's never exceeds $0.018$ (and is $0.0019$ at the end). Verlet's 1967 argon simulations used (3.32); Wisdom and Holman (1991) built solar-system integrators on the same splitting; it is also inside ST-6's Hamiltonian Monte Carlo and QC-9's Trotter formula.

**Integrability and chaos.** $n$ degrees of freedom with $n$ independent conserved quantities whose mutual brackets vanish make a system *integrable*. **Liouville–Arnold theorem (stated; Arnold ch. 10):** each compact connected level set is then an $n$-torus carrying action-angle variables with $H = H(J_1, \dots, J_n)$. Every 1-DOF system, Kepler and the symmetric top are integrable. The double pendulum has two degrees of freedom but only one conserved quantity, and at high energy it is chaotic. The KAM theorem (Kolmogorov 1954, Arnold 1963, Moser 1962): under small perturbations most tori survive, deformed, while those with nearly rational frequency ratios break into island chains and chaotic layers (Exercise C3).

## 4. Worked examples

### Example 4.1: Legendre transforms of $\frac12mv^2$ and $e^v$

(a) $f = \frac12mv^2$: (3.2) gives $p = mv^*$, $v^* = p/m$, and $g = p\cdot\frac{p}{m} - \frac{m}{2}\frac{p^2}{m^2} = \frac{p^2}{2m}$, with $g' = p/m = v^*$. This is $T \to p^2/2m$ inside $H$. With $m = 2$, $p = 3$: $v^* = 1.5$, $g = 2.25$.

(b) $f = e^v$: $p = e^{v^*}$, so $v^* = \ln p$, defined only for $p > 0$ (the slopes $e^v$ actually has). $g(p) = p\ln p - p$, $g' = \ln p = v^*$, $g'' = 1/p > 0$. Involution: $vp - (p\ln p - p)$ is maximized where $v = \ln p$, i.e. $p = e^v$, with value $ve^v - ve^v + e^v = e^v$. At $p = 2$: $v^* = 0.6931$, $g = 2\ln 2 - 2 = -0.6137$.

### Example 4.2: The spherical pendulum, reduced by $p_\phi$

Mass $m$, rod $l$, polar angle $\theta$ from the downward vertical, azimuth $\phi$: $T = \frac12ml^2(\dot\theta^2 + \sin^2\theta\,\dot\phi^2)$, $V = -mgl\cos\theta$, $p_\theta = ml^2\dot\theta$, $p_\phi = ml^2\sin^2\theta\,\dot\phi$, and by Theorem 3.3
$$
H = \frac{p_\theta^2}{2ml^2} + \frac{p_\phi^2}{2ml^2\sin^2\theta} - mgl\cos\theta.
$$
$\phi$ is cyclic, $p_\phi$ is constant, and $\theta$ moves in $V_{\text{eff}} = \frac{p_\phi^2}{2ml^2\sin^2\theta} - mgl\cos\theta$.

*Conical orbit:* $V_{\text{eff}}'(\theta_0) = -\frac{p_\phi^2\cos\theta_0}{ml^2\sin^3\theta_0} + mgl\sin\theta_0 = 0$ gives $p_\phi^2 = \frac{m^2gl^3\sin^4\theta_0}{\cos\theta_0}$, and $\dot\phi = \frac{p_\phi}{ml^2\sin^2\theta_0}$ gives $\dot\phi^2 = \frac{g}{l\cos\theta_0}$.

*Small oscillations about it:* $V_{\text{eff}}'' = \frac{p_\phi^2(\sin^2\theta + 3\cos^2\theta)}{ml^2\sin^4\theta} + mgl\cos\theta$. The equilibrium condition gives $\frac{p_\phi^2}{ml^2\sin^4\theta_0} = \frac{mgl}{\cos\theta_0}$, so $V_{\text{eff}}''(\theta_0) = \frac{mgl(1 + 3\cos^2\theta_0)}{\cos\theta_0}$ and $\Omega^2 = \frac{V_{\text{eff}}''}{ml^2} = \frac{g}{l}\frac{1 + 3\cos^2\theta_0}{\cos\theta_0}$.

*Numbers* ($m = 1$ kg, $l = 1$ m, $g = 9.81$ m/s², $\theta_0 = \pi/6$): $p_\phi = \sqrt{9.81\times 0.0625/0.8660} = 0.841413$ kg m²/s; $\dot\phi = \sqrt{9.81/0.8660} = 3.36565$ rad/s (one revolution in $1.86686$ s); $\Omega = 6.0675$ rad/s; $\Omega/\dot\phi = \sqrt{3.25} = 1.8028$. The ratio is irrational, so a perturbed conical orbit traces a rosette. The project integrates Hamilton's equations: $\theta$ stays at $\theta_0$ to $10^{-12}$, and the measured $\Omega$ is $6.0673$.

### Example 4.3: The pendulum's separatrix

For $H = \frac{p^2}{2ml^2} - mgl\cos\theta$, $\nabla H = 0$ at $(0, 0)$ and $(\pm\pi, 0)$. The Hessian $\operatorname{diag}(mgl\cos\theta, 1/ml^2)$ is positive definite at $(0, 0)$ (elliptic) and indefinite at $(\pi, 0)$ (hyperbolic). The separatrix energy is $E_s = H(\pi, 0) = mgl$. On it, $\frac{p^2}{2ml^2} = mgl(1 + \cos\theta) = 2mgl\cos^2\frac{\theta}{2}$, so
$$
p = \pm 2ml\sqrt{gl}\,\cos\tfrac{\theta}{2},
$$
two arcs meeting at $\theta = \pm\pi$: the eye of the §2 sketch, drawn thick in `phase_portraits.png`. Enclosed area: $\int_{-\pi}^{\pi}4ml\sqrt{gl}\cos\frac{\theta}{2}\,d\theta = 16\,ml\sqrt{gl}$. With $m = l = 1$, $g = 9.81$: $E_s = 9.81$ J, area $50.113$ J·s.

*Motion on it:* with $\omega_0 = \sqrt{g/l}$, $\dot\theta = 2\omega_0\cos\frac{\theta}{2}$. With $u = \theta/2$, $du/\cos u = \omega_0\,dt$, and $\frac{d}{du}\operatorname{artanh}(\sin u) = 1/\cos u$, so $\sin\frac{\theta}{2} = \tanh(\omega_0t)$. The top is reached only as $t \to \infty$, which is why the period diverges as $E \to E_s$ (Exercise C1).

### Example 4.4: $\{L_x, L_y\} = L_z$

Bilinearity: $\{yp_z - zp_y, zp_x - xp_z\} = \{yp_z, zp_x\} - \{yp_z, xp_z\} - \{zp_y, zp_x\} + \{zp_y, xp_z\}$. By (3.15) a bracket of monomials is nonzero only if a conjugate pair is split between the two sides.

- $\{yp_z, zp_x\}$: $(z, p_z)$ split; Leibniz gives $yp_x\{p_z, z\} = -yp_x$.
- $\{yp_z, xp_z\}$ and $\{zp_y, zp_x\}$: no split pair; zero.
- $\{zp_y, xp_z\}$: $(z, p_z)$ split; $xp_y\{z, p_z\} = xp_y$.

Total $xp_y - yp_x = L_z$. The project verifies all nine $\{L_i, L_j\}$, $\{L^2, L_i\} = 0$, and $\{L_z, H\} = \{L^2, H\} = 0$ for $H = |\mathbf{p}|^2/2m + V(r)$ in SymPy.

### Example 4.5: The oscillator's polar map is canonical

Invert (3.21): $\tan Q = m\omega q/p$, $P = (p^2 + m^2\omega^2q^2)/2m\omega$. With $D = p^2 + m^2\omega^2q^2$: $\partial_qQ = m\omega p/D$, $\partial_pQ = -m\omega q/D$, $\partial_qP = m\omega q$, $\partial_pP = p/m\omega$, so
$$
\{Q, P\} = \frac{m\omega p}{D}\cdot\frac{p}{m\omega} + \frac{m\omega q}{D}\cdot m\omega q = \frac{D}{D} = 1.
$$
In the other direction, the Jacobian of $(Q, P) \mapsto (q, p)$ has determinant $\cos^2Q + \sin^2Q = 1$. By Theorem 3.7 the map is canonical, $K = \omega P$, and $P = E/\omega$, $Q = \omega t + Q_0$. For contrast, $(Q, P) = (q, p^2)$ has $\det M = 2p$: not canonical. The project's checker prints $M^\top\mathsf{J}M - \mathsf{J}$ for both.

### Example 4.6: The oscillator by Hamilton–Jacobi

(3.28) reads $\frac{1}{2m}S_q^2 + \frac12m\omega^2q^2 + S_t = 0$.

1. Separate, $S = W(q) - Et$, with $E$ the new constant momentum: $W' = \sqrt{2mE - m^2\omega^2q^2}$.
2. $\beta = \partial S/\partial E = \int\frac{m\,dq}{\sqrt{2mE - m^2\omega^2q^2}} - t = \frac{1}{\omega}\arcsin\frac{q}{A} - t$, with $A = \sqrt{2E/m\omega^2}$.
3. Invert: $q = A\sin\omega(t + \beta)$, and $p = \partial S/\partial q = W' = \sqrt{2mE}\cos\omega(t + \beta)$.
4. In full, $W = \frac{q}{2}\sqrt{2mE - m^2\omega^2q^2} + \frac{E}{\omega}\arcsin\frac{q}{A}$.

No ODE was solved: one separation, one integral, one inversion. The project checks symbolically that $S = W - Et$ satisfies the PDE, then with $m = 1.3$, $\omega = 0.7$, $E = 2.0$, $\beta = 0.3$ checks along the trajectory that $\partial S/\partial q = p(t)$ to $10^{-9}$, that $\partial S/\partial E$ stays at $0.300000$, and that $S - S(0) = \int L\,dt$ to $10^{-5}$.

### Example 4.7: A pendulum whose length halves

$m = g = 1$, $l_0 = 1$, amplitude $\theta_a = 0.10$ rad; the string is pulled up to $l_0/2$ over 100 periods. Final amplitude?

1. Small angles: $\omega = \sqrt{g/l}$ and $E = mgl(1 - \cos\theta_a) \approx \frac12mgl\theta_a^2$.
2. By (3.27), $E/\omega = \frac12mg^{1/2}l^{3/2}\theta_a^2$ is constant, so $\theta_a \propto l^{-3/4}$.
3. Halving $l$: $\theta_a \to 0.10\times 2^{3/4} = 0.1682$ rad. $E$ rises by $\sqrt2$: you do work pulling the string against the tension.

The project integrates the exact equations with $l$ decreasing linearly. Slow: amplitude ratio $1.6823$ against $1.6818$ predicted, $E/\omega$ constant to $0.14\%$. Fast (half a period): ratio $1.5669$, $E/\omega$ off by $13.3\%$ (`adiabatic_invariant.png`).

## 5. Why it matters and where it connects

- **SM-3, SM-4:** $\{f(H), H\} = 0$ in (3.11) makes $\delta(H - E)$ and $e^{-\beta H}$ stationary; $d^{2n}z$ is the measure microstates are counted in. **SM-2:** $F = U - TS$ is §3.1's transform.
- **QM-4:** the dictionary (3.19), Heisenberg's equation from (3.14), the algebra (3.17), and $(\hat{\mathbf{p}} - q_e\mathbf{A})^2/2m$ from (3.9). Unitary transformations are the quantum canonical transformations.
- **QM-1, QM-8, QM-9:** Bohr–Sommerfeld (§3.9); WKB is (3.30) plus $\hbar^2A_{qq}/A$; the path integral starts from (3.7).
- **QC-9, ST-6:** Trotterization and Hamiltonian Monte Carlo are the splitting (3.31)–(3.32). HMC accepts with $\min(1, e^{-\Delta H})$ and no Jacobian factor because Verlet preserves volume.
- **CM-5, CM-7:** Kepler and the symmetric top are integrable; C2 shows why Kepler orbits close.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills. B1–B8 are the module; B1 (Lorentz force from (3.9)), B3 (hidden symmetry of the 2D oscillator), B5 (generating function for (3.21)) and B8 (project extension) are essential. C1 derives the pendulum's action and period, C2 Kepler's actions, C3 explores the standard map.

## 7. Coding project

`project.py` (about 15 s) checks every number in §3–§4 and the solutions, and saves eight figures: phase portraits (oscillator, pendulum with separatrix, double well $V = -x^2/2 + x^4/4$); the Legendre tangent lines; the Liouville blob (5000 interior points and a 4000-point boundary under the pendulum flow; at $t = 24$ the shoelace area ratio is $0.999998$ for symplectic Euler, $0.960$ for RK45, $1.202$ for explicit Euler); the signature figure `energy_drift.png` (RK4 vs Verlet, $10^4$ periods); the pendulum's $J(E)$ and $T = dJ/dE$; adiabatic invariance; $S(q, t)$; and WKB. In between, a SymPy `poisson_bracket` checks (3.15), (3.17), Jacobi and conservation laws, and a canonical checker evaluates $M^\top\mathsf{J}M - \mathsf{J}$.

**Build it yourself:** before reading Parts 4–5, write `symplectic_euler_step` and `verlet_pendulum` from (3.31)–(3.32) and reproduce `energy_drift.png`. If your Verlet curve drifts, you have a full kick at both ends or the kick and drift out of order. "Your turn" asks for a 4th-order symplectic method, Kepler action-angle variables, and the standard map.

## 8. Where you'll struggle

- **The Legendre transform is a formula until you see the tangent lines.** Draw Example 4.1(a) by hand: three tangents, their intercepts at $-g(p)$. Then §3.2 stops being symbol-pushing.
- **Canonical momentum is not $mv$.** In a magnetic field $\mathbf{p} = m\mathbf{v} + q_e\mathbf{A}$ (3.9). On a rotating hoop (B4) $H \ne T + V$. Compute $p = \partial L/\partial\dot q$ and $H = \sum p\dot q - L$ every time.
- **Signs flip between books.** Landau's bracket is minus (3.12); generating-function tables differ. Anchor everything to $\{q, p\} = +1$ and $\dot p = -\partial H/\partial q$, and re-derive a table entry from (3.22) rather than trusting memory.
- **Canonical transformations and generating functions stay abstract for most students.** The realistic target: brackets and Liouville fluent, a canonical map recognizable by $\{Q, P\} = 1$ or (3.20), Hamilton–Jacobi done once (Example 4.6, B7). Learn type 2 by deriving (3.23); look up the rest.
- **Liouville is about volume, not shape.** The blob becomes a long filament and looks as if it grew. Measure it.
- **The WKB page is a preview.** Follow its five steps and keep the slogan ($S$ is the phase); the turning-point matching is QM-8's.

## 9. Check yourself

1. Define $g(p)$ for convex $f(v)$; which two relations make it an involution?
2. Which Hamilton equation is an identity of the Legendre transform, and which is Euler–Lagrange?
3. Give a system where $H$ is conserved but is not $T + V$.
4. Compute $\nabla\cdot(\mathsf{J}\nabla H)$ in one line.
5. Write the Liouville equation; show $\rho = e^{-\beta H}$ is stationary.
6. Prove that conserved $L_x$ and $L_y$ force conserved $L_z$.
7. What is $\{F, p_x\}$, and what does $p_x$ generate?
8. Write the quantum versions of $\{x, p\} = 1$ and $\dot F = \{F, H\}$.
9. For one degree of freedom, what single number decides whether a map is canonical?
10. Derive the type-2 relations from (3.22).
11. What is $J$ for an oscillator, and what is $dE/dJ$?
12. Write the Hamilton–Jacobi equation; why is $\partial S/\partial q = p$?

## 10. Pair with

- **Taylor, *Classical Mechanics*, ch. 13.** Gentle; Liouville done carefully; good phase-portrait problems.
- **Goldstein, Poole & Safko, ch. 8–10.** The reference for canonical transformations, generating functions, Hamilton–Jacobi and action-angle variables (including Kepler). Use it as the lookup table.
- **Landau & Lifshitz, *Mechanics*, §40–52.** The shortest complete treatment: Jacobi (§42), Liouville (§46), Hamilton–Jacobi (§47–48), adiabatic invariants (§49–51). Opposite bracket sign.
- **MIT 8.09 (OCW), Hamiltonian lectures.** Canonical transformations and Hamilton–Jacobi at course pace.
- **Arnold, *Mathematical Methods of Classical Mechanics*, ch. 8 (optional).** The coordinate-free symplectic picture; why (3.20) is the right condition.
- **Susskind & Hrabovsky, *The Theoretical Minimum*, lectures 8–10.** Phase-space fluid, Liouville and brackets with almost no algebra; good for §2.
