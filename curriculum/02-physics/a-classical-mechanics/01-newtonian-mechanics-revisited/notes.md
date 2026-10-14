# CM-1: Newtonian Mechanics Revisited: Conservation Laws, and the Problem Constraints Create

**Time:** ~6 hours · **Prerequisites:** ODEs (and a first mechanics course, which this module assumes you passed) · **Unlocks:** CM-2, and through it CM-3 to CM-8

> **Conventions in this module.** SI units throughout; $g = 9.81\ \mathrm{m/s^2}$. Vectors are bold: $\mathbf{r} = (x, y, z)$, with $y$ pointing **up** wherever gravity appears, so gravitational potential energy is $V = mgy$. A dot is a time derivative: $\dot{x} = dx/dt$, $\ddot{x} = d^2x/dt^2$. Momentum is $\mathbf{p} = m\mathbf{v} = m\dot{\mathbf{r}}$ and Newton's second law is written $\mathbf{F} = \dot{\mathbf{p}}$. Potential energy is $V$ (not $U$), with $\mathbf{F} = -\nabla V$; kinetic energy is $T = \tfrac12 m v^2$; total energy $E = T + V$. Angular momentum about the origin is $\mathbf{L} = \mathbf{r}\times\mathbf{p}$ and torque is $\boldsymbol{\tau} = \mathbf{r}\times\mathbf{F}$. The pendulum angle $\theta$ is measured from the **downward** vertical. For the damped oscillator we write $\ddot{x} + \gamma\dot{x} + \omega_0^2 x = (F_0/m)\cos\omega t$ with $\gamma = b/m$ and $\omega_0^2 = k/m$; Taylor's $\beta$ is our $\gamma/2$. The quality factor is $Q = \omega_0/\gamma$.

## 0. Where this sits

You have done a mechanics course. You can write $F = ma$ for a block on a ramp, you know energy is conserved when there is no friction, and you have solved the harmonic oscillator. This module is not that course again. It is a six-hour pass over the same material with three purposes: to state the laws carefully enough that you notice what you were never told (what the first law actually says, when the third law fails, what "conservative" means), to derive the conservation laws instead of quoting them, and to make you feel, concretely and by hand, how badly Newton's formulation handles constraints.

That last item is the point. CM-2 (least action) and CM-3 (Lagrangian mechanics) exist to solve a problem, and you cannot appreciate the solution until you have suffered the problem. §3.10 is three worked constraint problems in raw Newtonian form, done in full. They are supposed to hurt. Everything you need from ODEs is here: linear constant-coefficient equations, separation of variables, and the idea that a second-order equation is a first-order system in $(x, \dot{x})$.

## 1. The problem that forced this idea

Newton's *Principia* (1687) gave three laws and the law of gravitation, and with them solved the problem that had defeated everyone: why the planets move on Kepler's ellipses. In principle the laws are complete. Given every force on every particle, $\mathbf{F} = \dot{\mathbf{p}}$ determines the future. Nothing in classical mechanics after 1687 adds a new law of nature. Everything after 1687 is about *formulation*: how to write the same physics so that problems you can state, you can also solve.

The trouble starts the moment something is not free. A bead slides on a curved wire. A pendulum bob hangs from a rod. A wheel rolls without slipping. In each case the body is forced to stay on a surface or a curve, and whatever holds it there (the wire pushing on the bead, the rod pulling the bob) is a force you did not choose and do not know in advance. It is whatever it has to be to keep the body on the surface. Newton's method requires you to introduce this unknown **constraint force**, write $\mathbf{F} = \dot{\mathbf{p}}$ in Cartesian components with the constraint force included, and then solve simultaneously for the motion *and* for the force, in coordinates ($x$, $y$, $z$) that have nothing to do with the geometry of the problem. For a pendulum you want the angle. Newton hands you $x(t)$, $y(t)$, and a tension $T(t)$ you were never interested in, coupled by the algebraic condition $x^2 + y^2 = \ell^2$.

Newton himself handled such problems (Book I, Section X of the *Principia* treats bodies moving on given surfaces and the pendulum), but by geometric arguments tailored to each case. Huygens had already found in 1673 that a pendulum swinging on a cycloid is exactly isochronous, again by geometry. There was no general machine. Making one took a century. Jacob Bernoulli's analysis of the compound pendulum (1703) and Johann Bernoulli's principle of virtual velocities (1717) contained the key idea: the forces of constraint do no work in any displacement the constraint allows, so if you project Newton's equations onto the allowed directions, the constraint forces vanish from the problem. D'Alembert made this systematic in his *Traité de dynamique* (1743): write $\mathbf{F} - \dot{\mathbf{p}} = 0$ and treat $-\dot{\mathbf{p}}$ as a fictitious extra force, so that every dynamics problem becomes a statics problem, and apply virtual work. Euler, in parallel, was the one who first wrote Newton's second law as the differential equations $m\ddot{x} = F_x$ and so on (1750) and who worked out rigid-body motion.

Lagrange finished it. His *Mécanique analytique* (1788) opens by boasting that the book contains no diagrams: the whole of mechanics has been reduced to analysis. You pick *any* coordinates that describe the configuration (an angle, an arc length, a relative separation), as many as there are genuine freedoms, and the equations of motion follow from a single scalar function by a fixed recipe. Constraint forces never appear. Hamilton (1834–35) then showed that the recipe is equivalent to a single statement about the whole path: the motion makes a certain integral stationary. That is where CM-2 begins.

So the arc of this module is: Newton's laws, done properly; the conservation laws they imply, derived; and then the pendulum, the bead, and the Atwood machine in Newtonian form, where you will watch an unknown force eat half the page.

## 2. Intuition

### Force is what changes momentum

Forget "force causes acceleration." Think: **momentum is a conserved substance, and force is its flow rate.** $\mathbf{F} = \dot{\mathbf{p}}$ says a force on a body is a rate at which momentum is being delivered to it. The third law then says momentum delivered to one body is taken from another: momentum is never created, only moved. The center-of-mass theorem in §3.2 is nothing more than this bookkeeping: internal transfers cancel in the total, so the total momentum of a system changes only by what flows in from outside.

This picture survives into field theory, where the naive third law does not (§3.1): two moving charges can push on each other with forces that are *not* equal and opposite, and the "missing" momentum is carried by the electromagnetic field. Force as momentum flow is the version of Newton that survives; "action equals reaction between pairs of bodies" is not.

### Conservation laws are symmetries in disguise

Three quantities are conserved in isolated systems: momentum, angular momentum, and energy. In this module each is derived from Newton's laws with a specific assumption: momentum needs the third law, angular momentum needs the internal forces to be central, energy needs the forces to come from a potential. It looks like three separate facts. It is not. CM-4 (Noether's theorem) shows that each conservation law is the shadow of a symmetry: if the physics does not care *where* you are, momentum is conserved; if it does not care which *direction* you face, angular momentum; if it does not care *when* you start, energy. Hold onto that as you go through §3; the assumptions in each derivation are exactly the places where a symmetry is being used without being named.

### The state of a system lives in phase space

Newton's second law is second order in time. To predict the future of a particle you need its position *and* its velocity now. So the "state" of a one-dimensional particle is not a point on a line; it is a point $(x, \dot{x})$ in a plane, and the motion is a curve in that plane. This plane (or, with $p = m\dot{x}$ in place of $\dot{x}$, the $(x, p)$ plane) is **phase space**. Through every point of it passes exactly one trajectory, because the ODE has a unique solution from each initial condition. Trajectories never cross.

For a conservative one-dimensional system the trajectories are the level curves of $E = \tfrac12 m\dot{x}^2 + V(x)$. That means you can read the entire motion off a sketch of $V(x)$ without solving anything. Here is the pendulum:

```
 V(θ)                             θ-dot
   |      _.-"""-._                 |   ~~~~~~~~  (E > 2mgl: over the top)
2mgl .-"           "-.        ------+------  separatrix (E = 2mgl)
   |"                 "        (  ( | )  )   closed loops (E < 2mgl)
   |                           (   (o)   )
   +-------+--------+-- θ     ------+------
  -π       0        π         -π    0    π   θ
```

A minimum of $V$ is a point the system circles around (closed loops: oscillation). A maximum is a point trajectories approach and veer away from (a saddle). The curve that passes exactly through the maximum, with $E$ equal to the maximum of $V$, divides the oscillating motions from the rotating ones; it is the **separatrix**. Guess now what §3.6 proves: the frequency of the small loops around a minimum is set by the *curvature* of $V$ there, and nothing else.

### A constraint is a surface the motion is glued to

Picture a bead on a wire bent into a curve in the plane. The bead has two Cartesian coordinates but only one freedom: where along the wire it is. The wire enforces this by pushing on the bead with a force $\mathbf{N}$ perpendicular to the wire (frictionless wire: no tangential force). Two things follow immediately. First, $\mathbf{N}$ does no work, because it is perpendicular to the velocity, which is tangent to the wire. So energy conservation holds for the bead *without ever knowing $\mathbf{N}$*. Second, if you project Newton's equation onto the tangent direction, $\mathbf{N}$ drops out, and what is left is one equation in one unknown, the position along the wire.

That is the whole idea of CM-3 in a sentence: **use coordinates that live on the constraint surface, and the constraint forces disappear.** Newton's formulation does not let you do this directly, because $\mathbf{F} = \dot{\mathbf{p}}$ is a vector equation in Cartesian components and the constraint is an algebraic side condition. §3.10 shows what it costs to do the projection by hand. The guess you should make: there ought to be a formulation that does the projection automatically, for any coordinates, from a single scalar function. There is. It is the Lagrangian.

## 3. Formal development

### 3.1 Newton's laws, stated carefully

**First law.** There exist reference frames, called **inertial frames**, in which a body subject to no force moves in a straight line at constant velocity.

Read that again. It is not a special case of the second law ($\mathbf{F} = 0 \Rightarrow \mathbf{a} = 0$). It is the statement that frames exist in which the second law *holds*. In a rotating frame, a free body does not move in a straight line, and the second law is false unless you add fictitious forces. The first law tells you which frames the other two laws are about.

**Second law.** In an inertial frame, the net force on a particle equals the rate of change of its momentum:
$$\mathbf{F} = \dot{\mathbf{p}} = \frac{d}{dt}(m\mathbf{v}). \tag{3.1}$$
For constant mass this is $\mathbf{F} = m\ddot{\mathbf{r}}$. We write it as (3.1) for two reasons: it is the form that survives into relativity and quantum mechanics, and it is the form you need when mass changes (§3.9), where the naive $\mathbf{F} = m\mathbf{a}$ produces wrong answers.

**Third law.** If body 1 exerts a force $\mathbf{F}_{21}$ on body 2, then body 2 exerts $\mathbf{F}_{12} = -\mathbf{F}_{21}$ on body 1. In its **strong form** the two forces also lie along the line joining the bodies (they are *central*). Gravity and electrostatics obey the strong form. Magnetic forces between moving charges obey neither: two charges moving at right angles to each other can have $\mathbf{F}_{12}$ and $\mathbf{F}_{21}$ that are not even antiparallel. Momentum is still conserved, because the field carries momentum (EM-6, Poynting's theorem). For this module all internal forces obey the strong third law.

**Galilean relativity.** Let frame $S'$ move at constant velocity $\mathbf{u}$ relative to inertial frame $S$: $\mathbf{r}' = \mathbf{r} - \mathbf{u}t$, $t' = t$.

▶ **Try it first.** Differentiate twice. What happens to $\ddot{\mathbf{r}}$?

1. $\dot{\mathbf{r}}' = \dot{\mathbf{r}} - \mathbf{u}$. Velocities differ by the constant $\mathbf{u}$.
2. $\ddot{\mathbf{r}}' = \ddot{\mathbf{r}}$. Accelerations are identical.
3. If the force depends only on relative positions $\mathbf{r}_i - \mathbf{r}_j$ and relative velocities $\dot{\mathbf{r}}_i - \dot{\mathbf{r}}_j$ (true for every fundamental force), it is also unchanged. So $\mathbf{F}' = m\ddot{\mathbf{r}}'$ holds in $S'$ exactly when it holds in $S$. ∎

■ **What just happened.** Every frame moving uniformly relative to an inertial frame is inertial. There is no experiment within mechanics that detects absolute velocity. Energy and momentum are *not* invariant under this change (a body at rest in $S'$ has kinetic energy in $S$), but the *laws* are, and in particular the work–energy theorem and the conservation laws hold in every inertial frame.

### 3.2 Momentum of a system and the center-of-mass theorem

Consider $N$ particles with masses $m_i$ and positions $\mathbf{r}_i$. Particle $i$ feels an external force $\mathbf{F}_i^{\text{ext}}$ and internal forces $\mathbf{F}_{ij}$ from each other particle $j \neq i$. Define the total mass $M = \sum_i m_i$, the total momentum $\mathbf{P} = \sum_i m_i\dot{\mathbf{r}}_i$, and the **center of mass**
$$\mathbf{R} = \frac{1}{M}\sum_i m_i\mathbf{r}_i.$$

**Theorem 3.1 (Center-of-mass theorem).** If the internal forces obey the third law, then
$$\dot{\mathbf{P}} = \mathbf{F}^{\text{ext}} \equiv \sum_i \mathbf{F}_i^{\text{ext}}, \qquad\text{equivalently}\qquad M\ddot{\mathbf{R}} = \mathbf{F}^{\text{ext}}. \tag{3.2}$$

▶ **Try it first.** Write (3.1) for each particle, add them all up, and look at the double sum of internal forces. Pair the terms.

**Derivation.**

1. For particle $i$: $m_i\ddot{\mathbf{r}}_i = \mathbf{F}_i^{\text{ext}} + \sum_{j\neq i}\mathbf{F}_{ij}$.
2. Sum over $i$: $\sum_i m_i\ddot{\mathbf{r}}_i = \sum_i\mathbf{F}_i^{\text{ext}} + \sum_i\sum_{j\neq i}\mathbf{F}_{ij}$.
3. The double sum runs over every ordered pair $(i, j)$ with $i \neq j$. Each unordered pair $\{i, j\}$ appears twice, once as $\mathbf{F}_{ij}$ and once as $\mathbf{F}_{ji}$. Group them: $\sum_{i}\sum_{j\neq i}\mathbf{F}_{ij} = \sum_{i<j}(\mathbf{F}_{ij} + \mathbf{F}_{ji})$.
4. By the third law each bracket is $\mathbf{0}$. The double sum vanishes.
5. The left side of step 2 is $\frac{d}{dt}\sum_i m_i\dot{\mathbf{r}}_i = \dot{\mathbf{P}}$, and also $M\ddot{\mathbf{R}}$ by the definition of $\mathbf{R}$ (the masses are constant). ∎

■ **What just happened.** The center of mass of any system moves exactly as a single particle of mass $M$ would under the net external force alone. Internal forces, however violent (an explosion, Example 4.1), cannot move the center of mass. If $\mathbf{F}^{\text{ext}} = \mathbf{0}$, $\mathbf{P}$ is conserved. Step 4 is where "internal forces cancel" lives; it is one line, and it uses the third law and nothing else.

### 3.3 Angular momentum and torque

For a single particle, $\mathbf{L} = \mathbf{r}\times\mathbf{p}$ and $\boldsymbol{\tau} = \mathbf{r}\times\mathbf{F}$, both about the origin.

▶ **Try it first.** Differentiate $\mathbf{L}$ with the product rule. One of the two terms is a cross product of a vector with itself.

1. $\dot{\mathbf{L}} = \dot{\mathbf{r}}\times\mathbf{p} + \mathbf{r}\times\dot{\mathbf{p}}$.
2. $\dot{\mathbf{r}}\times\mathbf{p} = \dot{\mathbf{r}}\times m\dot{\mathbf{r}} = \mathbf{0}$.
3. $\mathbf{r}\times\dot{\mathbf{p}} = \mathbf{r}\times\mathbf{F} = \boldsymbol{\tau}$ by (3.1). So $\dot{\mathbf{L}} = \boldsymbol{\tau}$. ∎

**Theorem 3.2 (Angular momentum of a system).** For $N$ particles with $\mathbf{L} = \sum_i\mathbf{r}_i\times\mathbf{p}_i$, if the internal forces obey the **strong** third law (equal, opposite, and central), then
$$\dot{\mathbf{L}} = \boldsymbol{\tau}^{\text{ext}} \equiv \sum_i\mathbf{r}_i\times\mathbf{F}_i^{\text{ext}}. \tag{3.3}$$

**Derivation.**

1. Apply the one-particle result to each particle and sum: $\dot{\mathbf{L}} = \sum_i\mathbf{r}_i\times\mathbf{F}_i^{\text{ext}} + \sum_i\sum_{j\neq i}\mathbf{r}_i\times\mathbf{F}_{ij}$.
2. Pair the internal terms as in §3.2: $\sum_{i<j}\left(\mathbf{r}_i\times\mathbf{F}_{ij} + \mathbf{r}_j\times\mathbf{F}_{ji}\right)$.
3. Use $\mathbf{F}_{ji} = -\mathbf{F}_{ij}$: each bracket is $(\mathbf{r}_i - \mathbf{r}_j)\times\mathbf{F}_{ij}$.
4. If $\mathbf{F}_{ij}$ is central, it is parallel to $\mathbf{r}_i - \mathbf{r}_j$, and the cross product is $\mathbf{0}$. ∎

■ **What just happened.** The weak third law was enough for momentum. For angular momentum you need more: the internal forces must point along the line between the particles. Step 3 shows exactly why. If an internal force pair were equal and opposite but *offset* (a couple), it would exert a net internal torque and spin the system up from nothing. Nature's fundamental forces do not do this, but the assumption is real and you should know you made it.

### 3.4 Work, kinetic energy, and the work–energy theorem

Define $T = \tfrac12 m\,\dot{\mathbf{r}}\cdot\dot{\mathbf{r}}$.

▶ **Try it first.** Differentiate $T$ and use (3.1).

1. $\dot{T} = \tfrac12 m(\ddot{\mathbf{r}}\cdot\dot{\mathbf{r}} + \dot{\mathbf{r}}\cdot\ddot{\mathbf{r}}) = m\ddot{\mathbf{r}}\cdot\dot{\mathbf{r}} = \mathbf{F}\cdot\dot{\mathbf{r}}$.
2. Integrate from $t_1$ to $t_2$ along the actual path: $T_2 - T_1 = \int_{t_1}^{t_2}\mathbf{F}\cdot\dot{\mathbf{r}}\,dt = \int_{\mathbf{r}_1}^{\mathbf{r}_2}\mathbf{F}\cdot d\mathbf{r} \equiv W$. ∎

$$\Delta T = W = \int\mathbf{F}\cdot d\mathbf{r}\quad\text{(along the path actually taken, net force).} \tag{3.4}$$

■ **What just happened.** This is Newton's second law dotted with velocity and integrated. It holds for **every** force: friction, drag, magnetic, whatever. The common misstatement is "the work–energy theorem requires conservative forces." It does not. What conservative forces buy you is the ability to compute $W$ *without knowing the path*, which is the next section.

### 3.5 Conservative forces: three equivalent conditions

Let $\mathbf{F}(\mathbf{r})$ be a force field that depends on position only (not on time, not on velocity) defined on a region $D$ of space. Consider three statements:

- **(a) Path independence.** $\int_{\mathbf{r}_1}^{\mathbf{r}_2}\mathbf{F}\cdot d\mathbf{r}$ depends only on the endpoints, not the path. Equivalently, $\oint\mathbf{F}\cdot d\mathbf{r} = 0$ around every closed loop in $D$.
- **(b) Curl-free.** $\nabla\times\mathbf{F} = \mathbf{0}$ everywhere in $D$.
- **(c) Potential exists.** There is a scalar function $V(\mathbf{r})$ with $\mathbf{F} = -\nabla V$.

**Theorem 3.3.** (c) $\Rightarrow$ (a), (a) $\Rightarrow$ (c), and (c) $\Rightarrow$ (b) hold on any region. (b) $\Rightarrow$ (a) holds if $D$ is **simply connected** (every loop can be shrunk to a point without leaving $D$). A force satisfying them is **conservative**.

▶ **Try it first.** For (c) $\Rightarrow$ (a), parametrize a path as $\mathbf{r}(s)$ and use the chain rule on $V(\mathbf{r}(s))$. For (a) $\Rightarrow$ (c), *define* $V$ as minus the work from a fixed base point, then compute $\partial V/\partial x$ by moving a little in $x$.

**Derivation of (c) $\Rightarrow$ (a).**

1. Along any path $\mathbf{r}(s)$, $s \in [0, 1]$, from $\mathbf{r}_1$ to $\mathbf{r}_2$: $\int\mathbf{F}\cdot d\mathbf{r} = -\int_0^1\nabla V\cdot\frac{d\mathbf{r}}{ds}\,ds$.
2. By the chain rule, $\nabla V\cdot\frac{d\mathbf{r}}{ds} = \frac{d}{ds}V(\mathbf{r}(s))$.
3. So the integral is $-[V(\mathbf{r}(1)) - V(\mathbf{r}(0))] = V(\mathbf{r}_1) - V(\mathbf{r}_2)$, which depends only on the endpoints. ∎

**Derivation of (a) $\Rightarrow$ (c).**

1. Fix a base point $\mathbf{r}_0 \in D$ and define $V(\mathbf{r}) = -\int_{\mathbf{r}_0}^{\mathbf{r}}\mathbf{F}\cdot d\mathbf{r}'$ along any path. By (a) this is well defined.
2. Compute $\partial V/\partial x$ at $\mathbf{r}$. By path independence, $V(\mathbf{r} + h\hat{\mathbf{x}}) - V(\mathbf{r}) = -\int_{\mathbf{r}}^{\mathbf{r}+h\hat{\mathbf{x}}}\mathbf{F}\cdot d\mathbf{r}'$ along the straight segment, which is $-\int_0^h F_x(\mathbf{r} + s\hat{\mathbf{x}})\,ds$.
3. Divide by $h$ and let $h \to 0$: by continuity of $F_x$, the limit is $-F_x(\mathbf{r})$. Likewise for $y$ and $z$. So $\nabla V = -\mathbf{F}$. ∎

**Derivation of (c) $\Rightarrow$ (b).**

1. The $z$-component of $\nabla\times\mathbf{F}$ is $\partial_x F_y - \partial_y F_x = -\partial_x\partial_y V + \partial_y\partial_x V$.
2. For $V$ with continuous second derivatives, mixed partials commute (Clairaut), so this is $0$. Same for the other two components. ∎

**About (b) $\Rightarrow$ (a).** Stokes' theorem says $\oint_C\mathbf{F}\cdot d\mathbf{r} = \iint_S(\nabla\times\mathbf{F})\cdot d\mathbf{A}$ for any surface $S$ in $D$ bounded by $C$. If $D$ is simply connected such a surface always exists, and (b) makes the right side zero. Stokes' theorem itself is proved in EM-1 (the vector calculus toolkit), and the full statement that a curl-free field on a simply connected region is a gradient is the first half of EM-1's Helmholtz theorem. Without simple connectivity it is false: the field $\mathbf{F} = (-y, x, 0)/(x^2 + y^2)$ on $\mathbb{R}^3$ minus the $z$-axis has $\nabla\times\mathbf{F} = \mathbf{0}$ everywhere it is defined, but $\oint\mathbf{F}\cdot d\mathbf{r} = 2\pi$ around the unit circle (exercise B2). The hole matters.

■ **What just happened.** Three different-looking statements are one property. The useful direction in practice is (b): compute a curl to decide whether a potential exists, then find $V$ by integrating. Two things that are *not* conservative in this sense, and are routinely misdescribed: drag (depends on velocity) and the magnetic force (depends on velocity; it does no work, but it has no $V(\mathbf{r})$). Time-dependent forces are also excluded.

**Theorem 3.4 (Energy conservation).** If a particle moves under a conservative force $\mathbf{F} = -\nabla V$ plus possibly other forces $\mathbf{F}_{\text{nc}}$, then $E = T + V$ obeys
$$\dot{E} = \mathbf{F}_{\text{nc}}\cdot\dot{\mathbf{r}}. \tag{3.5}$$
In particular $E$ is constant when $\mathbf{F}_{\text{nc}} = \mathbf{0}$, and also when $\mathbf{F}_{\text{nc}}$ is always perpendicular to the velocity (a frictionless constraint force, a magnetic force).

**Derivation.**
1. $\dot{T} = \mathbf{F}\cdot\dot{\mathbf{r}} = (-\nabla V + \mathbf{F}_{\text{nc}})\cdot\dot{\mathbf{r}}$ by §3.4.
2. $\dot{V} = \nabla V\cdot\dot{\mathbf{r}}$ by the chain rule.
3. Add. The $\nabla V$ terms cancel. ∎

■ **What just happened.** Energy conservation is the work–energy theorem plus the existence of $V$. The second sentence of the theorem is the one that makes constrained problems tractable at all: a frictionless wire or an inextensible rod exerts a force you do not know, but that force is perpendicular to the motion, so it cannot change $E$.

### 3.6 One-dimensional motion from energy alone

For a particle on a line with $E = \tfrac12 m\dot{x}^2 + V(x)$ conserved:
$$\dot{x} = \pm\sqrt{\frac{2}{m}\left(E - V(x)\right)}. \tag{3.6}$$
The motion is confined to where $V(x) \le E$. Points where $V(x) = E$ are **turning points**: $\dot{x} = 0$ there and the sign flips.

**Period of oscillation.** If the particle is trapped between turning points $x_1 < x_2$, the time from $x_1$ to $x_2$ is, from (3.6) with the $+$ sign, $\int dx/\dot{x}$, and the full period is twice that:
$$T_{\text{period}} = 2\int_{x_1}^{x_2}\frac{dx}{\sqrt{\tfrac{2}{m}(E - V(x))}} = \sqrt{2m}\int_{x_1}^{x_2}\frac{dx}{\sqrt{E - V(x)}}. \tag{3.7}$$
The integrand blows up at both ends like $1/\sqrt{x - x_1}$, which is integrable. For $V = \tfrac12 kx^2$ you can check that (3.7) gives $2\pi\sqrt{m/k}$ independent of $E$ (exercise A4 route); for anything else, the period depends on amplitude. Example 4.2 does the pendulum.

**Reading the phase portrait from $V$.** Phase-space trajectories are the curves $\tfrac12 m\dot{x}^2 + V(x) = E$. At each $x$, $\dot{x} = \pm\sqrt{2(E - V)/m}$: the curve is widest where $V$ is lowest. A minimum $x_0$ of $V$ is an **equilibrium**; nearby trajectories are closed loops around $(x_0, 0)$. A maximum $x_m$ is also an equilibrium, but trajectories near it are hyperbola-like: the curve with $E = V(x_m)$ passes through $(x_m, 0)$ and is the separatrix, dividing trajectories that turn back from those that pass over. Everything else about the motion (direction of travel: clockwise, since $\dot{x} > 0$ in the upper half-plane means $x$ increases) follows. The project draws this for the pendulum, where $V(\theta) = mg\ell(1 - \cos\theta)$, the maxima are at $\theta = \pm\pi$, and the separatrix is $E = 2mg\ell$.

**Small oscillations.** Let $x_0$ be a minimum of $V$, so $V'(x_0) = 0$ and $V''(x_0) > 0$.

▶ **Try it first.** Taylor-expand $V$ about $x_0$ to second order and write Newton's equation for $\xi = x - x_0$.

1. $V(x_0 + \xi) = V(x_0) + V'(x_0)\xi + \tfrac12 V''(x_0)\xi^2 + O(\xi^3) = V(x_0) + \tfrac12 V''(x_0)\xi^2 + O(\xi^3)$.
2. $m\ddot{\xi} = -\dfrac{dV}{dx} = -V''(x_0)\xi + O(\xi^2)$.
3. Dropping the higher-order terms: $\ddot{\xi} + \omega^2\xi = 0$ with
$$\omega^2 = \frac{V''(x_0)}{m}. \tag{3.8}$$ ∎

■ **What just happened.** Near any stable equilibrium, every system is a harmonic oscillator, with frequency set by the curvature of the potential. This single idea becomes CM-6 (normal modes: the same thing for many coordinates, where $V''$ is a matrix and the frequencies are generalized eigenvalues) and the harmonic approximation in QM-3 and SM-4.

### 3.7 Drag: linear and quadratic

A body of mass $m$ falling through a fluid, $y$ up, with velocity $v = \dot{y}$ (so $v < 0$ when falling).

**Linear drag** $\mathbf{F}_{\text{drag}} = -b\mathbf{v}$. Newton: $m\dot{v} = -mg - bv$. Define $\tau = m/b$ and $v_\infty = mg/b = g\tau$.

1. Rewrite: $\dot{v} = -(v + v_\infty)/\tau$. This is linear with constant coefficients.
2. The solution with $v(0) = v_0$: $v(t) = -v_\infty + (v_0 + v_\infty)e^{-t/\tau}$. (Check: $\dot{v} = -(v_0 + v_\infty)e^{-t/\tau}/\tau = -(v + v_\infty)/\tau$. ✓)
3. As $t \to \infty$, $v \to -v_\infty$: the **terminal velocity**, where drag balances weight. The approach is exponential with time constant $\tau$.

**Projectile with linear drag.** Horizontal: $m\dot{v}_x = -bv_x$, so $v_x = v_{x0}e^{-t/\tau}$ and $x(t) = v_{x0}\tau(1 - e^{-t/\tau})$. Vertical: integrate step 2 with $v_0 = v_{y0}$:
$$x(t) = v_{x0}\tau\left(1 - e^{-t/\tau}\right), \qquad y(t) = (v_{y0} + v_\infty)\tau\left(1 - e^{-t/\tau}\right) - v_\infty t. \tag{3.9}$$
The horizontal range can never exceed $x_\infty = v_{x0}\tau$: the trajectory has a vertical asymptote. Example 4.3 puts numbers in.

**Quadratic drag** $\mathbf{F}_{\text{drag}} = -c|\mathbf{v}|\mathbf{v}$. For a body falling from rest (take speed $s = -v > 0$ downward to avoid sign clutter): $m\dot{s} = mg - cs^2$.

1. Terminal speed where $\dot{s} = 0$: $s_\infty = \sqrt{mg/c}$. Write $m\dot{s} = mg(1 - s^2/s_\infty^2)$.
2. Separate: $\displaystyle\int\frac{ds}{1 - s^2/s_\infty^2} = g\,t$. The left side is $s_\infty\operatorname{artanh}(s/s_\infty)$.
3. So $s(t) = s_\infty\tanh(t/\tau)$ with $\tau = s_\infty/g$. ∎

Unlike linear drag, quadratic drag does not separate into independent $x$ and $y$ equations for a projectile ($|\mathbf{v}| = \sqrt{v_x^2 + v_y^2}$ couples them), and there is no closed form. The project integrates it numerically.

### 3.8 The damped driven oscillator: a complete worked review

$$m\ddot{x} + b\dot{x} + kx = F_0\cos\omega t \quad\Longleftrightarrow\quad \ddot{x} + \gamma\dot{x} + \omega_0^2 x = \frac{F_0}{m}\cos\omega t, \qquad \gamma = \frac{b}{m},\ \omega_0^2 = \frac{k}{m}.$$

This is a linear ODE with constant coefficients. Its general solution is (homogeneous solution) + (one particular solution). You know this from ODEs; here is the physics reading.

**Transient.** The homogeneous equation $\ddot{x} + \gamma\dot{x} + \omega_0^2 x = 0$ has solutions $e^{rt}$ with $r^2 + \gamma r + \omega_0^2 = 0$, so $r = -\gamma/2 \pm \sqrt{\gamma^2/4 - \omega_0^2}$. For weak damping ($\gamma < 2\omega_0$, **underdamped**), $r = -\gamma/2 \pm i\omega_1$ with $\omega_1 = \sqrt{\omega_0^2 - \gamma^2/4}$, and
$$x_h(t) = e^{-\gamma t/2}\left(C_1\cos\omega_1 t + C_2\sin\omega_1 t\right).$$
Whatever the initial conditions, this dies out on a time scale $2/\gamma$. It is called the transient.

**Steady state.** For the particular solution use complex exponentials.

▶ **Try it first.** Replace $\cos\omega t$ by $e^{i\omega t}$, look for $\tilde{x} = \tilde{A}e^{i\omega t}$, solve for the complex $\tilde{A}$, and take the real part.

1. Substitute $\tilde{x} = \tilde{A}e^{i\omega t}$: $(-\omega^2 + i\gamma\omega + \omega_0^2)\tilde{A} = F_0/m$.
2. So $\tilde{A} = \dfrac{F_0/m}{\omega_0^2 - \omega^2 + i\gamma\omega}$.
3. Write $\tilde{A} = Ae^{-i\delta}$ with $A = |\tilde{A}|$. Then $\operatorname{Re}[\tilde{A}e^{i\omega t}] = A\cos(\omega t - \delta)$ with
$$A(\omega) = \frac{F_0/m}{\sqrt{(\omega_0^2 - \omega^2)^2 + \gamma^2\omega^2}}, \qquad \tan\delta = \frac{\gamma\omega}{\omega_0^2 - \omega^2},\quad \delta \in [0, \pi]. \tag{3.10}$$
4. Since the real part of a solution of a real linear ODE with real-part forcing is a solution, $x_p = A\cos(\omega t - \delta)$ solves the original equation. ∎

■ **What just happened.** After the transient dies, the oscillator moves at the *driving* frequency $\omega$, not its own, with amplitude $A$ and lagging the drive by the phase $\delta$. The complex number $\tilde{A}/F_0 = 1/[m(\omega_0^2 - \omega^2 + i\gamma\omega)]$ is the system's **transfer function** $H(\omega)$: in F3 you will see that for *any* input $F(t)$, the steady-state output's Fourier transform is $H(\omega)\hat{F}(\omega)$. What we just computed by hand is one Fourier component.

**Resonance.** $A(\omega)$ is maximized where the denominator $(\omega_0^2 - \omega^2)^2 + \gamma^2\omega^2$ is minimized. Set $u = \omega^2$; the derivative $-2(\omega_0^2 - u) + \gamma^2 = 0$ gives $\omega_{\max}^2 = \omega_0^2 - \gamma^2/2$, slightly below $\omega_0$. At $\omega = \omega_0$ exactly:
$$A(\omega_0) = \frac{F_0}{m\gamma\omega_0} = \frac{F_0}{k}\cdot\frac{\omega_0}{\gamma} = Q\,\frac{F_0}{k}, \qquad \delta(\omega_0) = \frac{\pi}{2}.$$
The static displacement under a constant force $F_0$ is $F_0/k$; at resonance the response is $Q$ times larger. The phase goes from $0$ (far below resonance: the mass follows the force) through $\pi/2$ (at resonance: velocity in phase with force, maximum power absorbed) to $\pi$ (far above: the mass moves opposite to the force). The average power absorbed, $\langle F\dot{x}\rangle = \tfrac12 F_0\omega A\sin\delta$, has a peak of full width at half maximum exactly $\gamma$ (exercise B5), so $Q = \omega_0/\gamma = \omega_0/\Delta\omega_{\text{FWHM}}$. Example 4.4 gives the numbers for $m = 1$ kg, $k = 100$ N/m, $b = 2$ kg/s.

### 3.9 Variable mass: the rocket equation

A rocket of instantaneous mass $m(t)$ moves at velocity $v$ (one dimension, no gravity for now) and ejects exhaust backward at speed $u$ *relative to the rocket*.

The sign error everyone makes: write $F = \dot{p}$ for the rocket alone as $\dot{m}v + m\dot{v} = 0$ and conclude $m\dot{v} = -\dot{m}v$. This is wrong, and it is wrong in a way that makes the thrust depend on $v$, which would mean a rocket in deep space would not accelerate from rest. The mistake is that $\mathbf{F} = \dot{\mathbf{p}}$ applies to a *fixed collection of particles*. "The rocket" at time $t + dt$ is not the same collection as at time $t$; some of it has left. You must track the system that includes the exhaust.

▶ **Try it first.** Take the system to be the rocket plus the fuel it is about to eject in the next $dt$. Write its momentum before and after. Keep track of the sign of $dm$.

**Derivation.**

1. At time $t$ the system is the rocket, mass $m$, velocity $v$. Momentum $p(t) = mv$.
2. In time $dt$ the rocket's mass changes by $dm$ (which is **negative**), so it ejects a mass $-dm > 0$. The exhaust leaves at velocity $v - u$ in the inertial frame (the rocket moves at $v$; the exhaust goes backward at $u$ relative to it). The rocket is now at $v + dv$.
3. Momentum after: $p(t + dt) = (m + dm)(v + dv) + (-dm)(v - u)$.
4. Expand and drop the second-order term $dm\,dv$: $p(t + dt) = mv + m\,dv + v\,dm - v\,dm + u\,dm = mv + m\,dv + u\,dm$.
5. No external force, so $p(t + dt) = p(t)$: $m\,dv = -u\,dm$. Divide by $dt$:
$$m\dot{v} = -u\dot{m}. \tag{3.11}$$
Since $\dot{m} < 0$, the right side is positive: the thrust is $u|\dot{m}|$, independent of $v$.
6. Separate and integrate: $dv = -u\,dm/m$, so $v(t) - v_0 = u\ln\dfrac{m_0}{m(t)}$. ∎

■ **What just happened.** The exhaust velocity in step 2 is $v - u$, not $-u$; the terms $\pm v\,dm$ then cancel in step 4, which is what kills the spurious $\dot{m}v$ term. With gravity, add $-mg$ to the right side of (3.11): $m\dot{v} = -u\dot{m} - mg$, giving $v = v_0 + u\ln(m_0/m) - gt$ for a vertical burn. Lift-off requires $u|\dot{m}| > m_0 g$ at the start. Example 4.6 has numbers.

### 3.10 The problem constraints create

Three problems, each done in Newtonian form with the constraint force as an unknown. Read them for the shape of the work, not just the answer.

#### (a) The simple pendulum in Cartesian coordinates

Pivot at the origin, rod of length $\ell$, bob of mass $m$ at $(x, y)$ with $y$ up. Forces: gravity $(0, -mg)$ and the rod's tension $T$, directed from the bob toward the pivot, so $\mathbf{T} = -T(x, y)/\ell$. Unknowns: $x(t)$, $y(t)$, $T(t)$. Equations:
$$m\ddot{x} = -\frac{T}{\ell}x, \qquad m\ddot{y} = -\frac{T}{\ell}y - mg, \qquad x^2 + y^2 = \ell^2.$$
Two differential equations and one algebraic equation for three unknowns. The algebraic one must be differentiated to be useful.

1. Differentiate the constraint: $x\dot{x} + y\dot{y} = 0$. Again: $\dot{x}^2 + \dot{y}^2 + x\ddot{x} + y\ddot{y} = 0$, i.e. $x\ddot{x} + y\ddot{y} = -v^2$ where $v^2 = \dot{x}^2 + \dot{y}^2$.
2. Multiply the $x$-equation by $x$, the $y$-equation by $y$, and add: $m(x\ddot{x} + y\ddot{y}) = -\dfrac{T}{\ell}(x^2 + y^2) - mgy = -T\ell - mgy$.
3. Insert step 1: $-mv^2 = -T\ell - mgy$, so
$$T = \frac{mv^2}{\ell} - \frac{mgy}{\ell}. \tag{3.12}$$
We now know the tension in terms of the state, but we still have two coupled second-order equations for $x$ and $y$ with $T$ substituted, and they are ugly.
4. To get something usable, introduce the angle by hand: $x = \ell\sin\theta$, $y = -\ell\cos\theta$. Then $\dot{x} = \ell\cos\theta\,\dot{\theta}$, $\ddot{x} = \ell\cos\theta\,\ddot{\theta} - \ell\sin\theta\,\dot{\theta}^2$, $\dot{y} = \ell\sin\theta\,\dot{\theta}$, $\ddot{y} = \ell\sin\theta\,\ddot{\theta} + \ell\cos\theta\,\dot{\theta}^2$, and $v^2 = \ell^2\dot{\theta}^2$.
5. Substitute into the two Newton equations: $m\ell(\cos\theta\,\ddot{\theta} - \sin\theta\,\dot{\theta}^2) = -T\sin\theta$ and $m\ell(\sin\theta\,\ddot{\theta} + \cos\theta\,\dot{\theta}^2) = T\cos\theta - mg$.
6. Multiply the first by $\cos\theta$, the second by $\sin\theta$, and add: the $T$ terms cancel, the $\dot{\theta}^2$ terms cancel, and $m\ell\ddot{\theta} = -mg\sin\theta$:
$$\ddot{\theta} = -\frac{g}{\ell}\sin\theta. \tag{3.13}$$
7. Multiply the first by $-\sin\theta$, the second by $\cos\theta$, and add: $m\ell\dot{\theta}^2 = T - mg\cos\theta$, so $T = m\ell\dot{\theta}^2 + mg\cos\theta$, which agrees with (3.12) since $y = -\ell\cos\theta$ and $v^2 = \ell^2\dot{\theta}^2$. ∎

Seven steps, two of them pure algebra to make an unknown we never wanted cancel out, and the decisive move (step 4) was to abandon Cartesian coordinates for the angle, which Newton's formulation gave us no hint to do. **CM-3 does this in three lines.**

#### (b) A bead on a frictionless parabolic wire

Wire $y = kx^2$, $k > 0$, gravity down, bead of mass $m$. The wire exerts a normal force $\mathbf{N} = N\hat{\mathbf{n}}$ with $\hat{\mathbf{n}} = (-2kx, 1)/\sqrt{1 + 4k^2x^2}$ (perpendicular to the tangent $(1, 2kx)$). Unknowns: $x(t)$, $y(t)$, $N(t)$.

1. Newton: $m\ddot{x} = -\dfrac{2kx\,N}{\sqrt{1 + 4k^2x^2}}$, $\quad m\ddot{y} = \dfrac{N}{\sqrt{1 + 4k^2x^2}} - mg$. Constraint: $y = kx^2$.
2. Abbreviate $\lambda = N/\sqrt{1 + 4k^2x^2}$ so that $m\ddot{x} = -2kx\lambda$ and $m\ddot{y} = \lambda - mg$.
3. Differentiate the constraint twice: $\dot{y} = 2kx\dot{x}$, $\ddot{y} = 2k\dot{x}^2 + 2kx\ddot{x}$.
4. Put step 3 into the $y$-equation and solve for $\lambda$: $\lambda = mg + 2mk(\dot{x}^2 + x\ddot{x})$.
5. Put that into the $x$-equation: $m\ddot{x} = -2kx\left[mg + 2mk(\dot{x}^2 + x\ddot{x})\right]$. Collect the $\ddot{x}$ terms: $\ddot{x}(1 + 4k^2x^2) = -2kx(g + 2k\dot{x}^2)$:
$$\ddot{x} = -\frac{2kx\left(g + 2k\dot{x}^2\right)}{1 + 4k^2x^2}. \tag{3.14}$$
6. The normal force, from steps 2 and 4: $N = \lambda\sqrt{1 + 4k^2x^2} = m\left[g + 2k(\dot{x}^2 + x\ddot{x})\right]\sqrt{1 + 4k^2x^2}$, with $\ddot{x}$ from (3.14). At the bottom ($x = 0$, so $\ddot{x} = 0$):
$$N_{\text{bottom}} = m\left(g + 2k\dot{x}^2\right). \tag{3.15}$$ ∎

Check (3.15): the radius of curvature of $y = kx^2$ at the origin is $1/(2k)$, so $N - mg = mv^2/\rho = 2kmv^2$. ✓ Equation (3.14) is a one-dimensional equation in $x$, which is what we wanted, but the route there went through an unknown $N$, a unit normal, and an auxiliary $\lambda$. **CM-3 does this in three lines**, and the $\dot{x}^2$ term, which looks mysterious here, comes out automatically from the kinetic energy $\tfrac12 m(1 + 4k^2x^2)\dot{x}^2$.

#### (c) The Atwood machine

Masses $m_1$ and $m_2$ hang from a massless string over a massless frictionless pulley. Positions $y_1$, $y_2$ measured upward. Unknowns: $\ddot{y}_1$, $\ddot{y}_2$, and the tension $T$ (the same on both sides because the pulley is massless and frictionless).

1. Newton for each mass: $m_1\ddot{y}_1 = T - m_1 g$, $\quad m_2\ddot{y}_2 = T - m_2 g$.
2. The string is inextensible: $y_1 + y_2 = \text{const}$, so $\ddot{y}_2 = -\ddot{y}_1$.
3. Subtract the equations after substituting: $m_1\ddot{y}_1 + m_2\ddot{y}_1 = (T - m_1 g) - (T - m_2 g)$, i.e.
$$\ddot{y}_1 = \frac{m_2 - m_1}{m_1 + m_2}\,g, \qquad T = m_1(g + \ddot{y}_1) = \frac{2m_1m_2}{m_1 + m_2}\,g. \tag{3.16}$$ ∎

This one is easy, but notice the pattern: one more unknown ($T$) than there are freedoms, and one algebraic step to eliminate it. With two pulleys and three masses you have two tensions, three constraint relations, and a $5\times 5$ linear system. **CM-3 does this in three lines**, for any number of pulleys.

#### Generalized coordinates, informally

In each problem the body had two or more Cartesian coordinates but fewer genuine freedoms: one for the pendulum ($\theta$), one for the bead ($x$, or arc length), one for the Atwood machine ($y_1$). A set of numbers that fixes the configuration, with no constraints left over, is a set of **generalized coordinates**. Their count is the number of **degrees of freedom**: $3N$ ($2N$ in a plane) minus the number of independent constraint equations. The pendulum has $2 - 1 = 1$. Newton's laws do not know about generalized coordinates; you have to introduce them by hand (step 4 of the pendulum) and grind. CM-3's recipe takes any generalized coordinates as input.

### 3.11 Numerical integration of Newton's equations

Write $\ddot{x} = a(x)$ as the first-order system $\dot{x} = v$, $\dot{v} = a(x)$, and step forward with time step $h$.

**Euler's method:** $x_{n+1} = x_n + hv_n$, $v_{n+1} = v_n + ha(x_n)$. Simple and, for oscillators, wrong in a specific way.

▶ **Try it first.** Apply one Euler step to the harmonic oscillator $a(x) = -\omega^2 x$ and compute the new energy $E = \tfrac12(v^2 + \omega^2x^2)$ (per unit mass) in terms of the old.

1. $v_{n+1}^2 = v_n^2 - 2h\omega^2x_nv_n + h^2\omega^4x_n^2$.
2. $\omega^2x_{n+1}^2 = \omega^2x_n^2 + 2h\omega^2x_nv_n + h^2\omega^2v_n^2$.
3. Add: the cross terms cancel, and $v_{n+1}^2 + \omega^2x_{n+1}^2 = (1 + h^2\omega^2)(v_n^2 + \omega^2x_n^2)$. So
$$E_{n+1} = (1 + h^2\omega^2)\,E_n, \qquad E_N = (1 + h^2\omega^2)^N E_0. \tag{3.17}$$ ∎

■ **What just happened.** Euler's method multiplies the energy by a fixed factor greater than one *every step*, exactly, no matter how small $h$ is. The project uses $h = T_0/1000$ with $T_0 = 2\pi/\omega$, so $h\omega = 2\pi/1000 = 0.006283$ and the factor per step is $1 + 3.948\times 10^{-5}$: after 10 periods ($10^4$ steps) $E/E_0 = 1.484$, after 100 periods $E/E_0 = 51.8$. For a short run it looks fine. It is not fine. Geometrically, the Euler map is the linear map with matrix $\begin{pmatrix}1 & h\\ -h\omega^2 & 1\end{pmatrix}$, whose determinant is $1 + h^2\omega^2$: it inflates phase-space area by that factor each step. The true flow preserves phase-space area exactly (CM-8, Liouville's theorem), so Euler is wrong about the one thing that matters most over long times.

**Euler on a nonlinear force.** For general $a(x) = -V'(x)$ (per unit mass), expand: $\tfrac12(v_n + ha_n)^2 = \tfrac12 v_n^2 + hv_na_n + \tfrac12h^2a_n^2$ and $V(x_n + hv_n) = V(x_n) - hv_na_n + \tfrac12h^2V''(x_n)v_n^2 + O(h^3)$. The $O(h)$ terms cancel, leaving
$$E_{n+1} - E_n = \tfrac12h^2\left[a(x_n)^2 + V''(x_n)\,v_n^2\right] + O(h^3), \tag{3.18}$$
which is exactly $h^2\omega^2E_n$ when $V = \tfrac12\omega^2x^2$. For the pendulum $V'' = (g/\ell)\cos\theta$ shrinks as the swing grows, so the gain is slower than geometric. It is still fatal. In the project ($\theta_0 = 1$ rad, $h = T_0/1000$) Euler's energy is up 41% after 10 periods, a plot that looks like a slightly growing sine wave. It crosses the separatrix at period 61, and by period 100 ($E/E_0 = 5.9$) the "pendulum" is spinning over the top, dozens of full turns. Euler did not make a small error; it changed the kind of motion.

**RK4** (classical fourth-order Runge–Kutta) evaluates $a$ four times per step and has local error $O(h^5)$. The default of `scipy.integrate.solve_ivp`, RK45 (Dormand–Prince), is an adaptive cousin. It is excellent for a few hundred periods, but it does not preserve phase-space area either; the energy drifts, slowly and steadily, and over long enough runs the drift wins.

**Symplectic Euler:** $v_{n+1} = v_n + ha(x_n)$, then $x_{n+1} = x_n + hv_{n+1}$. One tiny change: use the *updated* velocity for the position step. The map is now $(x, v) \mapsto (x + hv + h^2a(x),\ v + ha(x))$, whose Jacobian has determinant exactly 1 for *any* $a(x)$ (exercise C2). It preserves phase-space area, like the true flow. Its energy error is $O(h)$ but does *not* accumulate: the energy oscillates around the true value and stays bounded forever. This is the simplest **symplectic integrator**; CM-8 explains why "preserves the symplectic structure" is the right property to ask for, and the Störmer–Verlet (leapfrog) method you will meet there is the second-order version, used in every molecular dynamics code and in the Hamiltonian Monte Carlo of ST-6. The project runs all three on the pendulum for 100 periods.

## 4. Worked examples

### Example 4.1: An exploding shell and the center-of-mass theorem

A shell is fired at $v_0 = 50$ m/s at $60°$ above horizontal. At the top of its flight it explodes into two equal fragments. One fragment's velocity immediately after the explosion is zero. Where does the other land?

*Before the explosion.* $v_x = 50\cos 60° = 25.000$ m/s, $v_y = 50\sin 60° = 43.301$ m/s. Time to apex $v_y/g = 4.414$ s, apex height $v_y^2/2g = 95.57$ m, horizontal position at apex $25\times 4.414 = 110.35$ m. The unexploded shell would land at $R = v_0^2\sin 120°/g = 2500\times 0.8660/9.81 = 220.70$ m.

*The explosion.* Internal forces only, so by Theorem 3.1 the center of mass continues on the original parabola and lands at $R = 220.70$ m. Momentum conservation across the explosion (it is brief, so gravity's impulse is negligible): $(2m)(25, 0) = m(0, 0) + m\mathbf{v}_2$, so $\mathbf{v}_2 = (50, 0)$ m/s.

*After.* Fragment 1 falls straight down from the apex and lands at $x_1 = 110.35$ m. Fragment 2 falls for the same time (same height, same zero vertical velocity), $4.414$ s, while moving at 50 m/s horizontally: $x_2 = 110.35 + 50\times 4.414 = 331.05$ m. Check the CM: $(110.35 + 331.05)/2 = 220.70$ m $= R$. ✓

You could have skipped the velocities entirely: the CM lands at $R$, fragment 1 at $R/2$, so fragment 2 lands at $2R - R/2 = 3R/2 = 331.05$ m. The project simulates both fragments and asserts this.

### Example 4.2: The pendulum period as an elliptic integral

From (3.13), energy conservation reads $\tfrac12\ell^2\dot{\theta}^2 + g\ell(1 - \cos\theta) = g\ell(1 - \cos\theta_0)$ for release from rest at $\theta_0$. Apply (3.7) with $m \to m\ell^2$ and $V \to mg\ell(1 - \cos\theta)$:
$$T_{\text{period}} = 4\sqrt{\frac{\ell}{2g}}\int_0^{\theta_0}\frac{d\theta}{\sqrt{\cos\theta - \cos\theta_0}}.$$
(The factor 4: a quarter period from $0$ to $\theta_0$.)

1. Use $\cos\theta = 1 - 2\sin^2(\theta/2)$: $\cos\theta - \cos\theta_0 = 2\left(\sin^2\tfrac{\theta_0}{2} - \sin^2\tfrac{\theta}{2}\right)$. Set $k = \sin(\theta_0/2)$.
2. Substitute $\sin(\theta/2) = k\sin\phi$, so $\phi$ runs from $0$ to $\pi/2$ as $\theta$ runs from $0$ to $\theta_0$. Then $\tfrac12\cos(\theta/2)\,d\theta = k\cos\phi\,d\phi$, i.e. $d\theta = \dfrac{2k\cos\phi\,d\phi}{\sqrt{1 - k^2\sin^2\phi}}$.
3. The square root in the denominator becomes $\sqrt{2(k^2 - k^2\sin^2\phi)} = \sqrt{2}\,k\cos\phi$.
4. Everything cancels except
$$T_{\text{period}} = 4\sqrt{\frac{\ell}{g}}\int_0^{\pi/2}\frac{d\phi}{\sqrt{1 - k^2\sin^2\phi}} = 4\sqrt{\frac{\ell}{g}}\,K(k), \qquad k = \sin\frac{\theta_0}{2}, \tag{4.1}$$
where $K$ is the complete elliptic integral of the first kind.
5. Small-angle limit and first correction. Expand $(1 - k^2\sin^2\phi)^{-1/2} = 1 + \tfrac12 k^2\sin^2\phi + O(k^4)$ and use $\int_0^{\pi/2}\sin^2\phi\,d\phi = \pi/4$: $K(k) = \tfrac{\pi}{2}\left(1 + \tfrac14 k^2 + O(k^4)\right)$. With $T_0 = 2\pi\sqrt{\ell/g}$ and $k^2 = \sin^2(\theta_0/2) = \theta_0^2/4 + O(\theta_0^4)$:
$$T_{\text{period}} = T_0\left(1 + \frac{\theta_0^2}{16} + O(\theta_0^4)\right). \tag{4.2}$$ ∎

*Numbers.* For $\theta_0 = 90°$: $k^2 = 0.5$, $K = 1.85407$, $T/T_0 = 2K/\pi = 1.18034$. The correction (4.2) gives $1 + (\pi/2)^2/16 = 1.15421$, low by 2.2%. For $\theta_0 = 30°$: $k^2 = 0.06699$, $T/T_0 = 1.01741$ versus $1.01713$ from (4.2). The next term, $11\theta_0^4/3072$ (exercise C1), brings $90°$ to $1.17601$, low by 0.37%. Note for the project: `scipy.special.ellipk(m)` takes the *parameter* $m = k^2$, not $k$. Passing $k$ is a classic bug.

### Example 4.3: A projectile with linear drag

$m = 0.1$ kg, $b = 0.05$ kg/s, so $\tau = 2$ s and $v_\infty = g\tau = 19.62$ m/s. Launch at 30 m/s, $45°$: $v_{x0} = v_{y0} = 21.213$ m/s.

From (3.9): $x(t) = 42.426(1 - e^{-t/2})$ m and $y(t) = 81.666(1 - e^{-t/2}) - 19.62t$ m. The horizontal asymptote is $x_\infty = 42.43$ m.

*Maximum height.* $v_y = 0$ when $(v_{y0} + v_\infty)e^{-t/\tau} = v_\infty$, i.e. $t = \tau\ln\frac{v_{y0} + v_\infty}{v_\infty} = 2\ln(2.0812) = 1.4659$ s, and $y_{\max} = 13.67$ m (no drag: 22.94 m).

*Time of flight and range.* $y(t) = 0$ is transcendental; solve numerically (the project uses `brentq`): $t_f = 3.4032$ s, and the range is $x(t_f) = 34.69$ m. Without drag: $t_f = 2v_{y0}/g = 4.325$ s, range $v_0^2/g = 91.74$ m. Linear drag with $\tau = 2$ s cut the range by 62%. The project overlays the closed form on a `solve_ivp` integration and asserts they agree to $10^{-6}$.

### Example 4.4: The driven oscillator at resonance, with numbers

$m = 1$ kg, $k = 100$ N/m, $b = 2$ kg/s, $F_0 = 1$ N. Then $\omega_0 = 10$ rad/s, $\gamma = 2$ s$^{-1}$, $Q = \omega_0/\gamma = 5$, $\omega_1 = \sqrt{100 - 1} = 9.950$ rad/s.

*Static response.* $F_0/k = 0.010$ m.

*At $\omega = \omega_0 = 10$.* From (3.10), $A = (F_0/m)/(\gamma\omega_0) = 1/20 = 0.050$ m $= Q\times 0.010$. Phase lag $\delta = \pi/2$: the displacement lags the force by a quarter cycle; velocity is in phase with force.

*The true peak.* $\omega_{\max} = \sqrt{\omega_0^2 - \gamma^2/2} = \sqrt{98} = 9.8995$ rad/s, where $A = 1/(\gamma\sqrt{\omega_0^2 - \gamma^2/4}) = 1/(2\sqrt{99}) = 0.050252$ m. Half a percent above $A(\omega_0)$; for $Q = 5$ the distinction barely matters.

*Off resonance.* At $\omega = 8$: $A = 1/\sqrt{36^2 + 16^2} = 1/\sqrt{1552} = 0.02538$ m, $\delta = \arctan(16/36) = 0.4182$ rad $= 23.96°$. At $\omega = 12$: $A = 1/\sqrt{44^2 + 24^2} = 0.01995$ m, $\delta = \pi - \arctan(24/44) = 2.6422$ rad $= 151.4°$ (second quadrant, because $\omega_0^2 - \omega^2 < 0$).

*Transient.* Decays as $e^{-\gamma t/2} = e^{-t}$; after 5 s it is $0.7\%$ of its initial size, after 15 s it is $3\times 10^{-7}$. The project integrates from rest for 20 s, measures the steady-state amplitude over the last 2 s, and overlays it on (3.10) across $\omega \in [2, 20]$. The power curve's full width at half maximum comes out as $2.000$ s$^{-1}$ $= \gamma$.

### Example 4.5: The bead on the parabolic wire, with the normal force at the bottom

$y = kx^2$ with $k = 0.5$ m$^{-1}$ (radius of curvature $1/2k = 1$ m at the bottom), $m = 0.2$ kg, released from rest at $x_0 = 1$ m (so $y_0 = 0.5$ m).

*Speed at the bottom.* The normal force does no work (Theorem 3.4), so $\tfrac12 mv^2 = mgy_0$: $v^2 = 2gkx_0^2 = 9.81$ m$^2$/s$^2$, $v = 3.132$ m/s.

*Normal force at the bottom.* By (3.15), $N = m(g + 2kv^2) = 0.2(9.81 + 2\times 0.5\times 9.81) = 0.2\times 19.62 = 3.924$ N, which is exactly $2mg$. In general $N_{\text{bottom}} = mg(1 + 4k^2x_0^2)$; for $2kx_0 = 1$ that is $2mg$.

*Normal force at release.* At $x = x_0$, $\dot{x} = 0$, so (3.14) gives $\ddot{x} = -2kx_0 g/(1 + 4k^2x_0^2) = -9.81/2 = -4.905$ m/s$^2$, and $N = m[g + 2kx_0\ddot{x}]\sqrt{1 + 4k^2x_0^2} = 0.2[9.81 - 4.905]\sqrt{2} = 1.3873$ N. Check: at rest the wire supports only the normal component of gravity; the wire's slope is $2kx_0 = 1$, i.e. $45°$, so $N = mg\cos 45° = 1.962/\sqrt2 = 1.3873$ N. ✓

The project integrates (3.14), computes $N(t)$ along the trajectory from step 6 of §3.10(b), and asserts both numbers.

### Example 4.6: Atwood numbers and a rocket burn

*Atwood.* $m_1 = 1.0$ kg, $m_2 = 1.5$ kg. From (3.16): $\ddot{y}_1 = 0.5\times 9.81/2.5 = 1.962$ m/s$^2$ (the lighter mass rises), $T = 2\times 1.5\times 9.81/2.5 = 11.772$ N. Sanity: $T$ lies between $m_1 g = 9.81$ N and $m_2 g = 14.715$ N, as it must (it has to lift $m_1$ and fail to hold $m_2$).

*Rocket.* $m_0 = 1000$ kg, final mass $300$ kg, exhaust speed $u = 2500$ m/s. From step 6 of §3.9: $\Delta v = 2500\ln(1000/300) = 2500\times 1.2040 = 3009.9$ m/s. If the burn is vertical and takes 60 s, subtract $gt = 588.6$ m/s: $2421.3$ m/s. Lift-off check: $|\dot{m}| = 700/60 = 11.67$ kg/s, thrust $u|\dot{m}| = 29{,}167$ N versus initial weight $9810$ N. ✓ If instead the burn took 5 minutes, the thrust would be $5833$ N and the rocket would sit on the pad burning fuel.

## 5. Why it matters and where it connects

- **CM-2** takes the hint from §3.10: there should be a formulation that projects out constraint forces automatically. Hamilton's principle is that formulation, and the Euler–Lagrange equation it produces is derived in F6 (calculus of variations) and applied in CM-2.
- **CM-3** redoes the pendulum (3.13), the bead (3.14), and the Atwood machine (3.16) from $L = T - V$ in generalized coordinates. Keep this module's derivations; you will compare line counts. Lagrange multipliers in CM-3 recover the tension (3.12) and the normal force (3.15) when you *do* want them.
- **CM-4** turns the three assumptions in §3.2–3.5 (third law, central forces, existence of $V$) into three symmetries (translation, rotation, time translation) via Noether's theorem.
- **CM-5** uses $\dot{\mathbf{L}} = \boldsymbol{\tau}$ with a central force, (3.3), to reduce the two-body problem to one dimension, then applies (3.6) and (3.7) to the effective potential to get Kepler orbits.
- **CM-6** generalizes (3.8): for $n$ coordinates, $V''(x_0)$ becomes the Hessian matrix $K$ and the kinetic energy gives a mass matrix $M$; normal-mode frequencies solve $\det(K - \omega^2 M) = 0$, a generalized eigenproblem (LA-8).
- **CM-8** is phase space taken seriously: Liouville's theorem is the exact version of "the true flow preserves area" from §3.11, and symplectic integrators are the numerical methods that respect it.
- **EM-1** proves Stokes' theorem and the Helmholtz theorem, closing the (b) $\Rightarrow$ (a) gap in Theorem 3.3. **EM-6** shows where the momentum goes when the third law fails for magnetic forces.
- **F3**: the complex amplitude $\tilde{A}/F_0$ of (3.10) is the transfer function; the response to an arbitrary force is its convolution with the oscillator's Green's function.
- **LA-7**: the Euler map in §3.11 is a $2\times 2$ matrix; its eigenvalues have modulus $\sqrt{1 + h^2\omega^2} > 1$, which is a second way to see (3.17).
- **ST-6**: Hamiltonian Monte Carlo integrates Hamilton's equations with leapfrog precisely because of the area-preservation property introduced in §3.11.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills: center of mass, curl tests for conservative forces, terminal velocity, small-oscillation frequency, Atwood numbers. B1–B8 are the module: B1 (splitting $\mathbf{L}$ and $T$ into CM and relative parts), B2 (the vortex field: curl-free but not conservative), B3 (phase portrait of a double well), B5 (power absorbed and the FWHM), B7 (Coulomb friction, code), and B8 (general bead on a wire) are essential. C1 (next term in the pendulum period), C2 (why symplectic Euler preserves area), and C3 (pendulum with a driven pivot, Newtonian style) each take an hour.

## 7. Coding project

`project.py` checks every number in the notes and solutions and makes §3.11 visible:

1. **Three integrators on the pendulum** for 100 periods: Euler, RK4 (hand-written, checked against `solve_ivp`), symplectic Euler; energy error versus time. Asserts (3.17) and the §3.11 numbers on the harmonic oscillator ($E/E_0 = 1.484$, $51.8$); on the pendulum, Euler up 41% at 10 periods and past the separatrix at period 61, RK4's error below $10^{-6}$ (it is $7\times10^{-11}$; at $h = T_0/100$ it drifts down to $7\times10^{-6}$), symplectic Euler bounded with no trend.
2. **Phase portrait** with the separatrix $E = 2mg\ell$.
3. **Period versus amplitude:** `solve_ivp` events versus (4.1), (4.2), and the C1 term.
4. **Drag:** (3.9) versus numerics (Example 4.3); quadratic drag versus none.
5. **Driven oscillator:** numerical steady state overlaid on (3.10); Example 4.4; the power FWHM.
6. **Bead on the wire:** integrates (3.14), computes $N(t)$, asserts Example 4.5 and energy conservation.
7. **Atwood, the Newtonian way:** solves the $3\times 3$ system for $(\ddot{y}_1, \ddot{y}_2, T)$ at every step; asserts (3.16).
8. **Shell and rocket:** Examples 4.1 and 4.6.
9. **SymPy and the rest:** re-derives (3.14), (3.17), the C1 series, the C2 identities, and the curl tests; a last cell asserts the remaining numbers in `solutions.md`.

**Build it yourself:** before reading the code, write `euler_step`, `rk4_step`, and `symplectic_euler_step` for a generic $\ddot{x} = a(x)$, and `pendulum_period_elliptic(theta0, l)`. Then compare. The "Your turn" section asks for Coulomb friction on the oscillator (linear amplitude decay; B7 derives the rate), the period of $V = \alpha x^4$ versus amplitude (B4 predicts $T \propto 1/A$), and a rocket simulation against the closed form of §3.9.

## 8. Where you'll struggle

- **Admitting you don't know this.** You passed an exam on it. The hardest thing in this module is reading "the work–energy theorem holds for all forces" or "a conservative force must depend on position only" and realizing your working definition was wrong. Test yourself: can you state, without looking, what hypothesis each of §3.2, §3.3, §3.5 needs? If not, that is the module's content, not a detail.
- **"Internal forces cancel."** It is one line (step 4 of §3.2), and most courses wave at it. The point is *which* law does the cancelling: the weak third law for momentum, the strong (central) third law for angular momentum. If you cannot say why Theorem 3.2 needs more than Theorem 3.1, redo both derivations.
- **Conservative forces.** Three statements, several implications, one of which (curl-free implies potential) needs a topological hypothesis you have probably never been told about. The vortex field in B2 is the counterexample; do it by hand. Separately: "does no work" is not the same as "conservative." The magnetic force does no work and has no potential.
- **The rocket sign.** You will write $\dot{m}v + m\dot{v} = F$ for the rocket and get a thrust proportional to $v$. The fix is to say out loud what system $\mathbf{F} = \dot{\mathbf{p}}$ is being applied to (rocket plus exhaust), and to write the exhaust velocity as $v - u$, not $-u$.
- **Euler's method looks fine.** Ten periods, the plot looks like a sine wave with a slightly growing amplitude, and you ship it. At a hundred periods the oscillator's energy is up fifty-fold and the project's pendulum is spinning over the top. (3.17) says the growth is *exact* and *geometric*; there is no $h$ small enough to fix it over long times. Use RK4 for accuracy over short runs, a symplectic method for long ones, and never Euler for an oscillator.
- **The constraint examples hurt. That is the point.** If §3.10(a) felt like busywork, you have understood it correctly: seven steps of algebra to eliminate a force you never wanted, and the only good move (switching to $\theta$) was one Newton's formulation does not suggest. Do not try to "get through" these. Sit with how bad it is. CM-3 will feel earned.
- **The phase lag $\delta$.** It lives in $[0, \pi]$, and $\arctan$ alone gives $(-\pi/2, \pi/2)$. Above resonance you must add $\pi$ (or use `arctan2`). Example 4.4 at $\omega = 12$ shows the correct value, $151°$, not $-29°$.
- **`ellipk` takes $m = k^2$.** Not $k$. Every student who codes (4.1) gets the wrong period once.

## 9. Check yourself

1. State Newton's first law so that it is not a special case of the second.
2. Why is $\mathbf{F} = \dot{\mathbf{p}}$ preferred over $\mathbf{F} = m\mathbf{a}$ in this module? Give one situation where they differ.
3. Which form of the third law does the center-of-mass theorem need? Which does $\dot{\mathbf{L}} = \boldsymbol{\tau}^{\text{ext}}$ need, and where in the derivation?
4. Does the work–energy theorem hold for a block sliding to rest under friction? Does energy conservation?
5. Give a force that is curl-free but not conservative, and say what hypothesis of Theorem 3.3 it violates.
6. Write the period of oscillation between turning points $x_1, x_2$ in $V(x)$ as an integral. Why is it finite despite the singular integrand?
7. From a sketch of $V(x)$ with two minima and one maximum, sketch the phase portrait. Where is the separatrix?
8. Derive $\omega^2 = V''(x_0)/m$ in three lines.
9. At resonance, what is the amplitude in terms of $Q$ and the static displacement $F_0/k$? What is the phase lag?
10. Derive the rocket equation $m\dot{v} = -u\dot{m}$, saying explicitly what system you apply $\mathbf{F} = \dot{\mathbf{p}}$ to.
11. For the pendulum in Cartesian coordinates, what are the unknowns and how many equations do you have? Which step eliminates the tension?
12. By what exact factor does Euler's method change the energy of a harmonic oscillator per step? Why does symplectic Euler not do this?

## 10. Pair with

- **Taylor, *Classical Mechanics*, ch. 1–4.** Ch. 1 states the laws with the care this module aims for; §1.5 on the third law and its failure for magnetism is the source of the §3.1 teaser. Ch. 2 on drag is more thorough than §3.7 (quadratic drag with both components, Taylor's Fig. 2.9). Ch. 3 has the rocket (§3.2) and the CM theorem; ch. 4 has the conservative-force conditions and the one-dimensional energy method. Taylor's ch. 5 (oscillations) is the full version of §3.8 if you want the overdamped and critically damped cases too.
- **Morin, *Introduction to Classical Mechanics*, ch. 1–5.** For problems. Morin's are harder and more inventive than Taylor's; do the starred ones in ch. 3 (oscillations) and ch. 5 (conservation of energy and momentum) if you want to be sure you are not fooling yourself. His §5.8 does the rocket with the same sign warning.
- **Feynman Lectures, vol. I, ch. 9–14.** Ch. 9 ("Newton's Laws of Dynamics") integrates a planetary orbit by hand with Euler-like steps and is the best thing ever written on what a differential equation *is*. Ch. 10 on momentum conservation, ch. 13–14 on work and potential energy, say in words what §3.4–3.5 say in symbols.
- **Hairer, Lubich & Wanner, *Geometric Numerical Integration*, ch. I** (only if §3.11 hooked you). The first chapter does exactly the Euler/symplectic-Euler/Störmer–Verlet comparison on the pendulum and the Kepler problem, with pictures.
