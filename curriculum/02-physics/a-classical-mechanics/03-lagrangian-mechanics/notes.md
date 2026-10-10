# CM-3: Lagrangian Mechanics: Generalized Coordinates, Constraints, Multipliers, Worked Systems

**Time:** ~12 hours · **Prerequisites:** CM-2 (Hamilton's principle, the Euler–Lagrange equation), F6 (constrained variational problems), CM-1 (constraints in Newtonian form) · **Unlocks:** CM-4 (Noether, cyclic coordinates, the energy function), CM-6 (normal modes: the mass matrix $M$ from §3.6), CM-7 (rigid bodies), ML-2 (Lagrange multipliers and KKT)

> **Conventions in this module.** The Lagrangian is $L = T - V$ (kinetic minus potential), in SI units. A system of $N$ particles has positions $\mathbf{r}_a$, $a = 1, \dots, N$, masses $m_a$. Generalized coordinates are $q_i$, $i = 1, \dots, n$, with velocities $\dot q_i = dq_i/dt$; $q$ without an index means the whole list. Generalized momentum is $p_i = \partial L/\partial \dot q_i$. Holonomic constraints are written $f_k(q, t) = 0$, $k = 1, \dots, m$, and their multipliers are $\lambda_k$. **Sign convention:** the constraint force associated with $f_k$ is $\lambda_k \nabla f_k$ (so $\lambda_k$ carries the sign of $f_k$: rewrite $f_k \to -f_k$ and $\lambda_k$ flips). Generalized forces are $Q_i = \sum_a \mathbf{F}_a \cdot \partial\mathbf{r}_a/\partial q_i$. A **virtual displacement** $\delta q_i$ is a displacement at frozen time ($\delta t = 0$); a real displacement is $dq_i$ and takes $dt$ with it. Angles of pendulums are measured from the downward vertical. $\boldsymbol\omega$ is an angular velocity vector; $\omega_0 = \sqrt{g/R}$ or $\sqrt{g/\ell}$ is a system's natural frequency. $g = 9.81\ \mathrm{m/s^2}$ in every numerical check. Repeated indices are **not** summed unless a $\sum$ is written.

## 0. Where this sits

CM-2 gave you the principle: the path a system takes makes $S = \int L\,dt$ stationary, and the Euler–Lagrange equation $\frac{d}{dt}\frac{\partial L}{\partial \dot q_i} = \frac{\partial L}{\partial q_i}$ is equivalent to Newton's second law. That module proved the equivalence for unconstrained particles in Cartesian coordinates. It did not tell you what to do with a bead on a wire, a block on a moving wedge, or two pendulums hanging from each other. This module does. It turns the principle into a **method**: choose coordinates that live on the constraint surface, write $T$ and $V$ in them, and the constraint forces vanish from the problem. When you want those forces back, you reintroduce them as Lagrange multipliers.

What you build here is used everywhere downstream. CM-4 takes the Lagrangians you write here and reads conservation laws off their symmetries. CM-6 takes the mass matrix $M_{ij}(q)$ derived in §3.6 and the small-angle double pendulum of Example 4.4 and turns them into the generalized eigenproblem $K\mathbf{a} = \Omega^2 M \mathbf{a}$. CM-7 writes rigid-body Lagrangians in Euler angles using exactly the procedure boxed in §3.10. ML-2 uses the multipliers of §3.8 under the name "dual variables": the KKT conditions of constrained optimization are the statics case of this module's equation (3.11).

## 1. The problem that forced this idea

Newton's laws are complete and, for constrained systems, almost unusable. Consider the oldest textbook example after the falling stone: a bead sliding on a bent wire. Newton says $m\ddot{\mathbf{r}} = m\mathbf{g} + \mathbf{N}$, where $\mathbf{N}$ is the force the wire exerts. But $\mathbf{N}$ is not given. It is whatever it needs to be to keep the bead on the wire. You have three equations, three unknown coordinates, and two more unknowns (the components of $\mathbf{N}$ normal to the wire), and the only extra information is geometric: the bead stays on the wire. CM-1 made you feel this pain on purpose. For a single bead you can push through. For a mechanism with ten links, or the solar system with its tidally locked moons, or a spinning top, the bookkeeping of unknown reaction forces is hopeless.

The way out came in three steps across seventy years.

**Virtual work (Johann Bernoulli, 1717).** In a letter to Varignon, Bernoulli stated the principle of "virtual velocities" for statics: a system of bodies is in equilibrium exactly when the total work done by the applied forces in any infinitesimal displacement *compatible with the constraints* is zero. The genius of the statement is what it leaves out: the constraint forces. A rigid rod, a smooth surface, a taut string: none of them do work when you move the system in a way they permit. So you never have to know them. Statics becomes a question about applied forces and allowed motions only.

**d'Alembert's principle (1743).** Jean le Rond d'Alembert, in the *Traité de dynamique*, saw how to make dynamics into statics. Write Newton's law as $\mathbf{F}_a - m_a\ddot{\mathbf{r}}_a = \mathbf{0}$ and read $-m_a\ddot{\mathbf{r}}_a$ as one more force (the "inertial force"). Then every instant of a motion is an equilibrium of applied, constraint, and inertial forces, and Bernoulli's principle applies: the total virtual work of the applied plus inertial forces vanishes, because the constraint forces contribute nothing. This is the physical content that makes the whole Lagrangian method go. Everything in §3 is an unpacking of the sentence *constraint forces do no virtual work*.

**Lagrange's equations (1788).** Joseph-Louis Lagrange's *Mécanique analytique* combined the two principles with a piece of pure calculus: if the configuration is described by any set of independent parameters $q_1, \dots, q_n$ (not necessarily Cartesian coordinates of anything), then d'Alembert's principle can be rewritten as $n$ ordinary differential equations in the $q_i$, with the kinetic energy $T$ and the potential $V$ as the only inputs. Lagrange was proud that the book contained no diagrams at all: "On ne trouvera point de Figures dans cet Ouvrage." Mechanics had become algebra. He also introduced, in the statics part of the same book, the *méthode des multiplicateurs*: instead of eliminating a constraint by choosing coordinates, keep the redundant coordinates and add to the equations a term $\lambda\,\partial f/\partial q_i$ for each constraint $f = 0$. Lagrange already read these terms as the forces the constraints exert. Over the nineteenth century, as the equations were applied to machinery and to celestial mechanics, this became the standard way to recover reaction forces from a Lagrangian calculation: the multiplier is the force you threw away, returned to you on demand. Today it is the backbone of every multibody dynamics code, and (ML-2) of constrained optimization.

One more piece of vocabulary came late. Heinrich Hertz, in his *Principles of Mechanics* (1894), coined the word **holonomic** for constraints that can be written as equations among the coordinates, $f(q, t) = 0$, to distinguish them from constraints like rolling without slipping, which restrict velocities and cannot be integrated to a condition on positions. Lagrange's method in its simple form needs holonomic constraints. §3.2 explains why, and §4.10 is honest about what to do otherwise.

## 2. Intuition

### A constrained system lives on a surface

Take a single particle in $\mathbb{R}^3$ confined to a surface $f(x, y, z) = 0$, say a sphere. Its **configuration space** before the constraint is all of $\mathbb{R}^3$. After the constraint it is the sphere: a two-dimensional surface inside the three-dimensional space. The particle has two **degrees of freedom**, and any two numbers that locate a point on the sphere (latitude and longitude, say) are **generalized coordinates** for it.

```
                 N = λ ∇f   (constraint force: normal to the surface)
                   ↑
          . - - - -●- - - - .        the surface f(x,y,z) = 0
       .      ←— δr —→         .     δr: a virtual displacement, tangent
     .   allowed motions lie     .
    .    IN the surface           .
     .                           .
       .                       .
          ' - - - - - - - - '
```

Now the key physical fact. A smooth surface can push on the particle only **perpendicular** to itself. It has no grip along the surface (that would be friction). So the constraint force is $\mathbf{N} = \lambda\nabla f$ for some scalar $\lambda$, because $\nabla f$ is the normal direction. And any displacement $\delta\mathbf{r}$ that keeps the particle on the surface is tangent to it, so $\mathbf{N}\cdot\delta\mathbf{r} = 0$. **That is d'Alembert's principle for one particle on one surface**: the constraint force does no work in displacements the constraint allows.

### Why choosing coordinates on the surface kills the constraint force

Here is the guess you should make before §3 proves it. Newton's law, dotted with any tangent direction $\partial\mathbf{r}/\partial q_i$ (the direction the particle moves when you nudge the coordinate $q_i$):

$$
m\ddot{\mathbf{r}}\cdot\frac{\partial\mathbf{r}}{\partial q_i} = \mathbf{F}^{\text{applied}}\cdot\frac{\partial\mathbf{r}}{\partial q_i} + \underbrace{\mathbf{N}\cdot\frac{\partial\mathbf{r}}{\partial q_i}}_{= 0}.
$$

The constraint force has disappeared, not because it is small but because we only ever look at the equations of motion *along the surface*. Two tangent directions for a sphere give two equations in two unknowns $q_1, q_2$. The system is closed with no reaction force in sight. The whole technical content of §3.4 is showing that the left side, $m\ddot{\mathbf{r}}\cdot\partial\mathbf{r}/\partial q_i$, can be computed from the kinetic energy alone, as $\frac{d}{dt}\frac{\partial T}{\partial\dot q_i} - \frac{\partial T}{\partial q_i}$. Once you believe that, Lagrange's equations are Newton's law projected onto the constraint surface.

### Why the multiplier is a price

Suppose you refuse to choose surface coordinates and insist on $(x, y, z)$. Then you must carry $\mathbf{N} = \lambda\nabla f$ in the equations, with $\lambda$ unknown, and you need one more equation to find it: $f(x,y,z) = 0$ itself. Four unknowns ($x, y, z, \lambda$), four equations (three Newton, one constraint). The multiplier $\lambda$ is the size of the force needed to enforce the constraint. In optimization (ML-2) the same $\lambda$ is called a *shadow price*: how much the objective would improve if you relaxed the constraint by one unit. Here it is literally the force the wall exerts. If you want the tension in a pendulum string or the load on a bearing, this is how you get it.

So there are two ways to use the method, and you will use both:

| | Eliminate the constraint | Keep the constraint |
|---|---|---|
| Coordinates | $n$ independent $q_i$ on the surface | $n + m$ redundant coordinates |
| Unknowns | $q_1, \dots, q_n$ | $q_1, \dots, q_{n+m}$ and $\lambda_1, \dots, \lambda_m$ |
| Equations | $n$ Lagrange equations | $n + m$ Lagrange equations with $\sum_k \lambda_k\partial f_k/\partial q_i$ on the right, plus $m$ constraint equations |
| Constraint forces | gone | $\lambda_k\nabla f_k$ |
| When | you want the motion | you want the forces too (engineering), or the surface coordinates are ugly |

### What kinetic energy looks like in generalized coordinates

If $\mathbf{r}_a = \mathbf{r}_a(q)$, then by the chain rule $\dot{\mathbf{r}}_a = \sum_i \frac{\partial\mathbf{r}_a}{\partial q_i}\dot q_i$, so $T = \frac12\sum_a m_a|\dot{\mathbf{r}}_a|^2$ is a **quadratic form in the velocities** $\dot q_i$ with coefficients that depend on position: $T = \frac12\sum_{ij} M_{ij}(q)\dot q_i\dot q_j$. You should guess that $M$ is symmetric and positive definite (kinetic energy is never negative, and is zero only at rest), that it reduces to $m$ times the identity for a free particle in Cartesian coordinates, and that for a pendulum in the angle coordinate it is the single number $m\ell^2$, the moment of inertia. $M$ is the **mass matrix**, and in CM-6 it is one of the two matrices in $K\mathbf{a} = \Omega^2 M\mathbf{a}$. It is also, for the geometrically minded, a metric on configuration space (LA-12): the "length" of a velocity $\dot q$ is $\sqrt{2T}$.

### A preview of the showpiece

A bead slides on a vertical circular hoop that spins about its vertical diameter at angular speed $\omega$. Slow spin: the bead rests at the bottom, as a pendulum would. Fast spin: the bottom becomes unstable (think of a centrifuge), and the bead rides up to a position where gravity pulling it down the hoop balances the centrifugal tendency pushing it out. Guess: there is a critical $\omega$ at which the bottom stops being stable, and above it two new equilibria appear symmetrically on either side. The critical speed should be the one where the hoop's rotation rate matches the pendulum frequency $\sqrt{g/R}$, because that is the only frequency in the problem. Example 4.3 derives this exactly: it is a **pitchfork bifurcation** at $\omega^2 = g/R$, and the small-oscillation frequency about the new equilibria is computable in two lines from an effective potential.

## 3. Formal development

### 3.1 Configuration space and degrees of freedom

**Definition 3.1 (Configuration space).** For $N$ particles with positions $\mathbf{r}_1, \dots, \mathbf{r}_N \in \mathbb{R}^3$, the **configuration space** is $\mathbb{R}^{3N}$: one point in it specifies where every particle is. A **constraint** is a condition that restricts the system to a subset of this space. The **number of degrees of freedom** $n$ is the number of independent coordinates needed to specify a configuration once all constraints are imposed.

For $m$ independent constraint equations on $N$ particles in three dimensions, $n = 3N - m$. Two particles joined by a rigid rod: $3 \cdot 2 - 1 = 5$. A rigid body made of $N$ particles has $3N$ coordinates and so many rod constraints that only $6$ remain (three for position, three for orientation): CM-7. A bead on a wire: $3 - 2 = 1$. A planar double pendulum: $2 \cdot 2 - 2 = 2$ (each bob is in a plane, and each rod fixes a distance).

### 3.2 Holonomic and non-holonomic constraints

**Definition 3.2 (Holonomic constraint).** A constraint is **holonomic** if it can be written as an equation among the coordinates and possibly time:
$$
f_k(\mathbf{r}_1, \dots, \mathbf{r}_N, t) = 0, \qquad k = 1, \dots, m. \tag{3.1}
$$
If $t$ does not appear, the constraint is **scleronomic** (a fixed wire); if it does, **rheonomic** (a hoop that someone is spinning). Anything else is **non-holonomic**: inequalities ($|\mathbf{r}| \le R$, a particle inside a box), and conditions on velocities that cannot be integrated to a condition on positions.

The standard non-holonomic example is **rolling without slipping**. A vertical disk of radius $R$ rolling on the $xy$-plane has coordinates $(x, y)$ for its contact point, $\theta$ for the direction it faces, and $\phi$ for how far it has rotated. No slipping means the contact point moves at speed $R\dot\phi$ in the direction $\theta$:
$$
\dot x - R\dot\phi\cos\theta = 0, \qquad \dot y - R\dot\phi\sin\theta = 0.
$$
These are two conditions on velocities. If they came from a holonomic $f(x, y, \theta, \phi) = 0$, then every closed path in the configuration variables would return $\phi$ to its starting value. It does not: roll the disk forward a distance $d$ ($\phi$ increases by $d/R$), turn it in place by $90°$ ($\phi$ unchanged, since the contact point does not move), roll forward $d$ again, turn $90°$, and repeat twice more. The disk is back at $(x, y, \theta)$ but $\phi$ has increased by $4d/R$. No function of $(x, y, \theta, \phi)$ can be constant along that loop while vanishing at both ends, so the constraint is not holonomic. (Exercise C1 makes this a proof.) By contrast, a disk rolling along a *line* has $\dot x = R\dot\phi$, which integrates to $x - R\phi = \text{const}$: holonomic. Whether rolling is holonomic depends on the geometry, not on the word "rolling."

Why this matters: the method of §3.3–3.4 works by **solving the constraints** to express every $\mathbf{r}_a$ as a function of $n$ free coordinates. You can do that for (3.1). You cannot do it for a velocity constraint that does not integrate, because there is no surface to put coordinates on: the disk can reach every point of its four-dimensional configuration space, it just cannot get there along arbitrary paths. Lagrange's equations in the simple form need holonomic constraints. §4.10 says what the honest alternatives are.

### 3.3 Generalized coordinates and virtual displacements

**Definition 3.3 (Generalized coordinates).** Suppose $m$ holonomic constraints (3.1) leave $n = 3N - m$ degrees of freedom. A set of **generalized coordinates** is any $n$ quantities $q_1, \dots, q_n$ such that the positions can be written
$$
\mathbf{r}_a = \mathbf{r}_a(q_1, \dots, q_n, t), \qquad a = 1, \dots, N, \tag{3.2}
$$
with the constraints satisfied identically for every value of the $q_i$, and such that the $n$ vectors $\partial\mathbf{r}_a/\partial q_i$ (stacked over $a$ into one vector in $\mathbb{R}^{3N}$) are linearly independent at every configuration you care about.

The second requirement is what "good coordinates" means. Latitude and longitude on a sphere fail at the poles (longitude is undefined there, and $\partial\mathbf{r}/\partial\phi = \mathbf{0}$). That failure shows up later as a mass matrix that is not invertible (§3.6). Choosing generalized coordinates is the skill of this module and §8 says bluntly how it is acquired.

The $q_i$ need not have units of length. Angles, areas, a "center of mass plus relative displacement" combination, the arc length along a wire: anything that parametrizes the constraint surface smoothly.

**Definition 3.4 (Virtual displacement).** A **virtual displacement** $\delta\mathbf{r}_a$ is an infinitesimal change of configuration **at a fixed instant** ($\delta t = 0$) that is consistent with the constraints as they are at that instant. In generalized coordinates,
$$
\delta\mathbf{r}_a = \sum_{i=1}^n \frac{\partial\mathbf{r}_a}{\partial q_i}\,\delta q_i, \tag{3.3}
$$
with the $\delta q_i$ arbitrary and independent. Compare the **real** displacement in time $dt$:
$$
d\mathbf{r}_a = \sum_{i=1}^n \frac{\partial\mathbf{r}_a}{\partial q_i}\,dq_i + \frac{\partial\mathbf{r}_a}{\partial t}\,dt.
$$
The two differ by the $\partial\mathbf{r}_a/\partial t\,dt$ term, which is present only for rheonomic constraints. For a bead on a spinning hoop, a real displacement includes the hoop carrying the bead around; a virtual one slides the bead along the hoop with the hoop frozen. This distinction is where d'Alembert's principle is slippery, and Example 4.3 shows it biting.

### 3.4 d'Alembert's principle and the second derivation of Lagrange's equations

Split the total force on particle $a$ into the **applied** force $\mathbf{F}_a$ (gravity, springs, fields, anything you can write down) and the **constraint** force $\mathbf{F}_a^{c}$ (whatever the wire, rod, or surface exerts). Newton: $m_a\ddot{\mathbf{r}}_a = \mathbf{F}_a + \mathbf{F}^c_a$.

**Principle 3.5 (d'Alembert).** For ideal constraints, the constraint forces do no work under any virtual displacement:
$$
\sum_{a=1}^N \mathbf{F}^c_a\cdot\delta\mathbf{r}_a = 0 \quad\text{for every virtual displacement } \delta\mathbf{r}_a. \tag{3.4}
$$
Consequently, along the actual motion,
$$
\sum_{a=1}^N \left(m_a\ddot{\mathbf{r}}_a - \mathbf{F}_a\right)\cdot\delta\mathbf{r}_a = 0 \quad\text{for every virtual displacement.} \tag{3.5}
$$

(3.5) follows from (3.4) by substituting $\mathbf{F}^c_a = m_a\ddot{\mathbf{r}}_a - \mathbf{F}_a$. The content is (3.4), and it is a physical assumption about the constraints, not a theorem. It holds for: smooth surfaces and wires (force normal to surface, displacement tangent), rigid rods and inextensible strings (the two end forces are equal and opposite along the rod and the ends move the same amount along it), rolling without slipping (the contact point does not move, so the contact force does no work), and frictionless hinges. It fails for sliding friction, which is why friction has to be put in by hand as an applied force (§3.5, §4.10). Note that (3.4) says nothing about *real* work: a moving constraint can and does do real work on the system. The virtual displacement freezes time precisely to exclude that.

Now the derivation. The goal is to rewrite (3.5) entirely in terms of the generalized coordinates and the kinetic energy $T = \frac12\sum_a m_a|\dot{\mathbf{r}}_a|^2$.

▶ **Try it first.** Substitute (3.3) into (3.5) and collect the coefficient of each $\delta q_i$. The force term becomes $\sum_a \mathbf{F}_a\cdot\partial\mathbf{r}_a/\partial q_i$. For the acceleration term you need $\sum_a m_a\ddot{\mathbf{r}}_a\cdot\partial\mathbf{r}_a/\partial q_i$. Try to write it as a total time derivative of something minus a correction, using the product rule backwards: $\ddot{\mathbf{r}}\cdot\mathbf{u} = \frac{d}{dt}(\dot{\mathbf{r}}\cdot\mathbf{u}) - \dot{\mathbf{r}}\cdot\dot{\mathbf{u}}$. Then try to recognize both pieces as derivatives of $T$.

**Derivation (Lagrange's equations from d'Alembert's principle).**

1. **Generalized forces.** Substitute (3.3) into the force term of (3.5):
   $$\sum_a \mathbf{F}_a\cdot\delta\mathbf{r}_a = \sum_a\sum_i \mathbf{F}_a\cdot\frac{\partial\mathbf{r}_a}{\partial q_i}\,\delta q_i = \sum_i Q_i\,\delta q_i, \qquad Q_i := \sum_{a=1}^N \mathbf{F}_a\cdot\frac{\partial\mathbf{r}_a}{\partial q_i}. \tag{3.6}$$
   $Q_i$ is the **generalized force** conjugate to $q_i$. It has whatever units make $Q_i\,\delta q_i$ an energy: a torque if $q_i$ is an angle, a force if $q_i$ is a length.

2. **Two identities about $\mathbf{r}_a(q, t)$.** Differentiate (3.2) along the motion by the chain rule:
   $$\dot{\mathbf{r}}_a = \sum_j \frac{\partial\mathbf{r}_a}{\partial q_j}\dot q_j + \frac{\partial\mathbf{r}_a}{\partial t}.$$
   The right side depends on $q$, $\dot q$, and $t$, and it is *linear* in the $\dot q_j$ with coefficients $\partial\mathbf{r}_a/\partial q_j$ that do not depend on $\dot q$. Treating $q$, $\dot q$, $t$ as independent variables (as the Lagrangian formalism always does) and differentiating with respect to $\dot q_i$ picks out the coefficient:
   $$\frac{\partial\dot{\mathbf{r}}_a}{\partial\dot q_i} = \frac{\partial\mathbf{r}_a}{\partial q_i}. \tag{3.7}$$
   This is the **cancellation of dots**. Next, differentiate $\partial\mathbf{r}_a/\partial q_i$, itself a function of $(q, t)$, along the motion:
   $$\frac{d}{dt}\frac{\partial\mathbf{r}_a}{\partial q_i} = \sum_j \frac{\partial^2\mathbf{r}_a}{\partial q_j\partial q_i}\dot q_j + \frac{\partial^2\mathbf{r}_a}{\partial t\,\partial q_i} = \frac{\partial}{\partial q_i}\left(\sum_j\frac{\partial\mathbf{r}_a}{\partial q_j}\dot q_j + \frac{\partial\mathbf{r}_a}{\partial t}\right) = \frac{\partial\dot{\mathbf{r}}_a}{\partial q_i}, \tag{3.8}$$
   using equality of mixed partials and the fact that $\dot q_j$ is independent of $q_i$. This is the **interchange** of $d/dt$ and $\partial/\partial q_i$. Both identities need $\mathbf{r}_a$ to be a twice-differentiable function of $(q, t)$, nothing more.

3. **The acceleration term.** Product rule, then (3.8):
   $$\sum_a m_a\ddot{\mathbf{r}}_a\cdot\frac{\partial\mathbf{r}_a}{\partial q_i} = \frac{d}{dt}\left(\sum_a m_a\dot{\mathbf{r}}_a\cdot\frac{\partial\mathbf{r}_a}{\partial q_i}\right) - \sum_a m_a\dot{\mathbf{r}}_a\cdot\frac{d}{dt}\frac{\partial\mathbf{r}_a}{\partial q_i} = \frac{d}{dt}\left(\sum_a m_a\dot{\mathbf{r}}_a\cdot\frac{\partial\mathbf{r}_a}{\partial q_i}\right) - \sum_a m_a\dot{\mathbf{r}}_a\cdot\frac{\partial\dot{\mathbf{r}}_a}{\partial q_i}.$$

4. **Recognize $T$.** With $T = \frac12\sum_a m_a\,\dot{\mathbf{r}}_a\cdot\dot{\mathbf{r}}_a$:
   $$\frac{\partial T}{\partial\dot q_i} = \sum_a m_a\dot{\mathbf{r}}_a\cdot\frac{\partial\dot{\mathbf{r}}_a}{\partial\dot q_i} \overset{(3.7)}{=} \sum_a m_a\dot{\mathbf{r}}_a\cdot\frac{\partial\mathbf{r}_a}{\partial q_i}, \qquad \frac{\partial T}{\partial q_i} = \sum_a m_a\dot{\mathbf{r}}_a\cdot\frac{\partial\dot{\mathbf{r}}_a}{\partial q_i}.$$
   These are exactly the two pieces in step 3. Therefore
   $$\boxed{\;\sum_{a=1}^N m_a\ddot{\mathbf{r}}_a\cdot\frac{\partial\mathbf{r}_a}{\partial q_i} = \frac{d}{dt}\frac{\partial T}{\partial\dot q_i} - \frac{\partial T}{\partial q_i}\;} \tag{3.9}$$
   This identity is pure calculus. It holds for any smooth $\mathbf{r}_a(q, t)$, whether or not the $q_i$ are independent, and whether or not any law of physics is true.

5. **Assemble.** Put (3.6) and (3.9) into (3.5):
   $$\sum_{i=1}^n\left[\frac{d}{dt}\frac{\partial T}{\partial\dot q_i} - \frac{\partial T}{\partial q_i} - Q_i\right]\delta q_i = 0.$$
   The $\delta q_i$ are independent and arbitrary (Definition 3.3: the $q_i$ are free coordinates on the constraint surface). Take $\delta q_1 \ne 0$ and all others zero, then $\delta q_2$, and so on. Each bracket must vanish:
   $$\frac{d}{dt}\frac{\partial T}{\partial\dot q_i} - \frac{\partial T}{\partial q_i} = Q_i, \qquad i = 1, \dots, n. \tag{3.10}$$

6. **Conservative forces.** If the applied forces come from a potential, $\mathbf{F}_a = -\nabla_a V(\mathbf{r}_1, \dots, \mathbf{r}_N, t)$, then by the chain rule
   $$Q_i = -\sum_a \nabla_a V\cdot\frac{\partial\mathbf{r}_a}{\partial q_i} = -\frac{\partial V}{\partial q_i},$$
   where $V$ is now regarded as a function of $(q, t)$ through (3.2). Since $V$ does not depend on $\dot q$, $\partial V/\partial\dot q_i = 0$, and (3.10) becomes, with $L = T - V$,
   $$\frac{d}{dt}\frac{\partial L}{\partial\dot q_i} - \frac{\partial L}{\partial q_i} = 0. \qquad\blacksquare$$

■ **What just happened.** Newton's law was dotted into each tangent direction $\partial\mathbf{r}_a/\partial q_i$ of the constraint surface. The constraint forces dropped out at the start because they are perpendicular to those directions (that is (3.4)). The mass-times-acceleration side turned into derivatives of $T$ by two chain-rule identities. The result is the same Euler–Lagrange equation CM-2 obtained from Hamilton's principle, now derived without any variational calculus and with constraints built in from the first line. You have two independent routes to the same equation. The variational route (CM-2) explains *why* $L = T - V$ is the right object and generalizes to fields and to quantum mechanics; this route explains *what happened to the constraint forces* and extends to forces that have no potential.

**Theorem 3.6 (Lagrange's equations).** For a system with holonomic constraints satisfying d'Alembert's principle, described by independent generalized coordinates $q_1, \dots, q_n$, the motion satisfies (3.10). If all applied forces are conservative, it satisfies $\frac{d}{dt}\frac{\partial L}{\partial\dot q_i} = \frac{\partial L}{\partial q_i}$ with $L = T - V$. ∎ (Proved above.)

### 3.5 Non-conservative forces and velocity-dependent potentials

Split the applied force into a part with a potential and a part without: $\mathbf{F}_a = -\nabla_a V + \mathbf{F}_a^{\text{nc}}$. Step 6 above goes through for the first part and the second part stays as a generalized force:
$$
\frac{d}{dt}\frac{\partial L}{\partial\dot q_i} - \frac{\partial L}{\partial q_i} = Q_i^{\text{nc}}, \qquad Q_i^{\text{nc}} = \sum_a \mathbf{F}_a^{\text{nc}}\cdot\frac{\partial\mathbf{r}_a}{\partial q_i}. \tag{3.11}
$$
This is how driving forces, drag, and friction enter. Compute $Q_i^{\text{nc}}$ directly from (3.6) or, for linear drag, from Rayleigh's function (§3.9).

One class of velocity-dependent forces still fits inside $L$. If there is a function $U(q, \dot q, t)$ with
$$
Q_i = -\frac{\partial U}{\partial q_i} + \frac{d}{dt}\frac{\partial U}{\partial\dot q_i},
$$
then (3.10) with $Q_i$ moved to the left is the Euler–Lagrange equation of $L = T - U$. The Lorentz force is the important case: $U = q_e\Phi - q_e\dot{\mathbf{r}}\cdot\mathbf{A}$ (Example 4.8). Friction is not of this form. Rayleigh's function (§3.9) is the next best thing.

The generalized momentum is defined as
$$
p_i := \frac{\partial L}{\partial\dot q_i}. \tag{3.12}
$$
For a free particle in Cartesian coordinates it is $m\dot x_i$. For an angle it is an angular momentum. For a charged particle it is $m\dot{\mathbf{r}} + q_e\mathbf{A}$, which is not $m\mathbf{v}$ (Example 4.8). Lagrange's equations read $\dot p_i = \partial L/\partial q_i + Q_i^{\text{nc}}$: the rate of change of the generalized momentum equals the generalized force. If $q_i$ does not appear in $L$ (a **cyclic** coordinate) and $Q_i^{\text{nc}} = 0$, then $p_i$ is conserved. CM-4 makes this into Noether's theorem; here it appears as a working tool in Examples 4.6 and 4.7.

### 3.6 The structure of $T$: the mass matrix

▶ **Try it first.** Substitute $\dot{\mathbf{r}}_a = \sum_j (\partial\mathbf{r}_a/\partial q_j)\dot q_j + \partial\mathbf{r}_a/\partial t$ into $T = \frac12\sum_a m_a|\dot{\mathbf{r}}_a|^2$ and sort the result by powers of $\dot q$.

**Derivation.**

1. Expand the square:
   $$T = \frac12\sum_a m_a\left|\sum_j\frac{\partial\mathbf{r}_a}{\partial q_j}\dot q_j + \frac{\partial\mathbf{r}_a}{\partial t}\right|^2 = \frac12\sum_{i,j}M_{ij}\dot q_i\dot q_j + \sum_i N_i\dot q_i + T_0, \tag{3.13}$$
   with
   $$M_{ij}(q,t) = \sum_a m_a\frac{\partial\mathbf{r}_a}{\partial q_i}\cdot\frac{\partial\mathbf{r}_a}{\partial q_j}, \qquad N_i = \sum_a m_a\frac{\partial\mathbf{r}_a}{\partial q_i}\cdot\frac{\partial\mathbf{r}_a}{\partial t}, \qquad T_0 = \frac12\sum_a m_a\left|\frac{\partial\mathbf{r}_a}{\partial t}\right|^2.$$
   Write $T = T_2 + T_1 + T_0$ for the three pieces (quadratic, linear, and constant in $\dot q$). For scleronomic constraints $\partial\mathbf{r}_a/\partial t = \mathbf{0}$ and $T = T_2$ exactly.

2. **$M$ is symmetric**, since the dot product is.

3. **$M$ is positive semidefinite.** For any $\mathbf{u} \in \mathbb{R}^n$,
   $$\mathbf{u}^\top M\mathbf{u} = \sum_{ij}u_i M_{ij}u_j = \sum_a m_a\left|\sum_i\frac{\partial\mathbf{r}_a}{\partial q_i}u_i\right|^2 \ge 0.$$

4. **$M$ is positive definite exactly when the coordinates are good.** Equality in step 3 requires $\sum_i(\partial\mathbf{r}_a/\partial q_i)u_i = \mathbf{0}$ for every $a$, i.e. a linear dependence among the tangent vectors $\partial\mathbf{r}_a/\partial q_i$ with coefficients $u_i$. Definition 3.3 forbids that unless $\mathbf{u} = \mathbf{0}$. So $M$ is positive definite, hence invertible, wherever the coordinates are non-degenerate. ∎

■ **What just happened.** Stack the tangent vectors into the $3N\times n$ Jacobian $J = \partial\mathbf{r}/\partial q$ and let $\mathcal{M} = \operatorname{diag}(m_1, m_1, m_1, m_2, \dots)$. Then $M = J^\top\mathcal{M}J$: a **Gram matrix** (LA-5, LA-8) of the tangent vectors in the mass-weighted inner product. Gram matrices are symmetric positive semidefinite, and positive definite exactly when the vectors are independent. That is the whole proof in one line of linear algebra. $M$ is the mass matrix of CM-6 and the Riemannian metric on configuration space in LA-12's language: $2T_2 = \dot q^\top M\dot q$ is the squared speed in that metric.

Invertibility of $M$ is what lets you solve Lagrange's equations for the accelerations. Expanding $\frac{d}{dt}\partial L/\partial\dot q_i$ for $L = T_2 + T_1 + T_0 - V$, the only place $\ddot q$ appears is $\sum_j M_{ij}\ddot q_j$, so Lagrange's equations have the form
$$
\sum_j M_{ij}(q, t)\,\ddot q_j = (\text{terms in } q, \dot q, t), \tag{3.14}
$$
and $M^{-1}$ turns this into $\ddot q = (\ldots)$, which is what `solve_ivp` wants. The project's `lagrangian_to_eom` does exactly this: differentiate to find $M$, solve the linear system, hand the result to the integrator.

### 3.7 The energy function (brief; CM-4 does it properly)

Define
$$
h(q, \dot q, t) := \sum_i \dot q_i\,\frac{\partial L}{\partial\dot q_i} - L. \tag{3.15}
$$

▶ **Try it first.** Differentiate $h$ along a solution of (3.11) and simplify. You should find that everything cancels except $\partial L/\partial t$ and the non-conservative power.

**Derivation.** With $p_i = \partial L/\partial\dot q_i$,
$$
\frac{dh}{dt} = \sum_i(\ddot q_i p_i + \dot q_i\dot p_i) - \sum_i\left(\frac{\partial L}{\partial q_i}\dot q_i + \frac{\partial L}{\partial\dot q_i}\ddot q_i\right) - \frac{\partial L}{\partial t} = \sum_i\dot q_i\left(\dot p_i - \frac{\partial L}{\partial q_i}\right) - \frac{\partial L}{\partial t} = \sum_i\dot q_i Q_i^{\text{nc}} - \frac{\partial L}{\partial t}. \tag{3.16}
$$
∎

■ **What just happened.** If $L$ has no explicit time dependence and there are no non-conservative forces, $h$ is conserved. When the constraints are scleronomic and $V$ does not depend on $\dot q$, $T = T_2$ is homogeneous quadratic in $\dot q$, so $\sum_i\dot q_i\partial T/\partial\dot q_i = 2T$ (Euler's theorem on homogeneous functions) and $h = 2T - (T - V) = T + V = E$, the mechanical energy. When the constraints are rheonomic but $L$ still has no explicit $t$ (the rotating hoop at constant $\omega$), $h$ is still conserved but is *not* $T + V$: it is $T_2 - T_0 + V$. Example 4.3 lives on this distinction and the project checks it numerically. CM-4 names $h$ the energy function, turns it into the Hamiltonian in CM-8, and explains the symmetry behind its conservation.

### 3.8 Constraint forces via Lagrange multipliers

Now the second use of the method. Suppose you describe the system with $n' = n + m$ **redundant** coordinates $q_1, \dots, q_{n'}$ (for instance the Cartesian coordinates of every particle, $n' = 3N$), connected by $m$ holonomic constraints $f_k(q, t) = 0$ whose gradients $\nabla f_k = (\partial f_k/\partial q_1, \dots, \partial f_k/\partial q_{n'})$ are linearly independent. The $\delta q_i$ are no longer independent: a virtual displacement must keep every $f_k$ zero at frozen time,
$$
\delta f_k = \sum_{i=1}^{n'}\frac{\partial f_k}{\partial q_i}\,\delta q_i = 0, \qquad k = 1, \dots, m, \tag{3.17}
$$
that is, $\delta q \perp \nabla f_k$ for every $k$. (No $\partial f_k/\partial t\,\delta t$ term: $\delta t = 0$.)

▶ **Try it first.** Steps 1–5 of §3.4 used the independence of the $\delta q_i$ only in step 5. Redo step 5 when the $\delta q$ are restricted to the subspace $\mathcal{T} = \{\delta q : \delta q\cdot\nabla f_k = 0\ \forall k\}$. What can you conclude about the vector of brackets $E_i = \frac{d}{dt}\partial T/\partial\dot q_i - \partial T/\partial q_i - Q_i$? Use LA-5: if $E\cdot\delta q = 0$ for all $\delta q\in\mathcal{T} = W^\perp$, then $E\in (W^\perp)^\perp = W$.

**Derivation (the multiplier equations).**

1. The identity (3.9) and the generalized-force formula (3.6) never used independence of the $q_i$: they hold for any smooth parametrization $\mathbf{r}_a(q, t)$, redundant or not. So d'Alembert's principle (3.5) still gives
   $$\sum_{i=1}^{n'}E_i\,\delta q_i = 0, \qquad E_i := \frac{d}{dt}\frac{\partial T}{\partial\dot q_i} - \frac{\partial T}{\partial q_i} - Q_i,$$
   but now only for $\delta q \in \mathcal{T}$, the subspace of $\mathbb{R}^{n'}$ orthogonal to $W := \operatorname{span}\{\nabla f_1, \dots, \nabla f_m\}$.

2. $E\cdot\delta q = 0$ for all $\delta q\in W^\perp$ says $E\in(W^\perp)^\perp$. In finite dimensions $(W^\perp)^\perp = W$ (LA-5). So $E$ is a linear combination of the constraint gradients:
   $$E_i = \sum_{k=1}^m\lambda_k\frac{\partial f_k}{\partial q_i}$$
   for some scalars $\lambda_1, \dots, \lambda_m$, unique because the $\nabla f_k$ are independent.

3. Split $Q_i = -\partial V/\partial q_i + Q_i^{\text{nc}}$ as in §3.5 and write $L = T - V$:
   $$\boxed{\;\frac{d}{dt}\frac{\partial L}{\partial\dot q_i} - \frac{\partial L}{\partial q_i} = Q_i^{\text{nc}} + \sum_{k=1}^m\lambda_k\frac{\partial f_k}{\partial q_i}, \qquad i = 1, \dots, n',\;} \tag{3.18}$$
   together with the $m$ constraint equations $f_k(q, t) = 0$. That is $n' + m$ equations for the $n' + m$ unknown functions $q_i(t)$ and $\lambda_k(t)$.

4. **The multipliers are the constraint forces.** Go back to Newton with the constraint force kept: $m_a\ddot{\mathbf{r}}_a = \mathbf{F}_a + \mathbf{F}^c_a$. Dot with $\partial\mathbf{r}_a/\partial q_i$, sum over $a$, and use (3.9) on the left:
   $$\frac{d}{dt}\frac{\partial T}{\partial\dot q_i} - \frac{\partial T}{\partial q_i} = Q_i + Q_i^c, \qquad Q_i^c := \sum_a\mathbf{F}^c_a\cdot\frac{\partial\mathbf{r}_a}{\partial q_i}.$$
   Compare with step 2: $E_i = Q_i^c$, so
   $$Q_i^c = \sum_{k=1}^m\lambda_k\frac{\partial f_k}{\partial q_i}. \tag{3.19}$$
   The generalized constraint force conjugate to $q_i$ is the multiplier term. In particular, if the redundant coordinates are the Cartesian coordinates of the particles, then $\partial\mathbf{r}_a/\partial q_i$ is a unit vector and (3.19) reads
   $$\mathbf{F}^c_a = \sum_{k=1}^m\lambda_k\nabla_a f_k. \tag{3.20}$$
   The force of constraint on particle $a$ is $\lambda_k$ times the gradient of the constraint function: perpendicular to the constraint surface, with magnitude set by $\lambda_k$. ∎

■ **What just happened.** D'Alembert's principle says the constraint force is orthogonal to every allowed displacement. The allowed displacements are the orthogonal complement of the constraint gradients. So the constraint force is in the span of the gradients: $\mathbf{F}^c = \sum\lambda_k\nabla f_k$. The multipliers are the coefficients. Everything else is bookkeeping. The variational version of the same result: F6's multiplier rule says that extremizing $\int L\,dt$ subject to $f_k(q,t) = 0$ is the same as extremizing $\int\big(L + \sum_k\lambda_k(t)f_k\big)dt$ freely in $q$ and $\lambda$; the Euler–Lagrange equations for $q_i$ are (3.18) with $Q^{\text{nc}} = 0$, and those for $\lambda_k$ are $f_k = 0$. Same equations, same sign.

**Sign convention, stated once more because it will bite you.** The force is $+\lambda_k\nabla f_k$. If you write the pendulum constraint as $f = x^2 + y^2 - \ell^2$, then $\nabla f$ points radially outward, and the tension (inward) comes out as a negative $\lambda$. If you write $f = \ell^2 - x^2 - y^2$, $\lambda$ is positive. Neither is wrong; what is wrong is reading $|\lambda|$ as the force without checking which way $\nabla f$ points. Example 4.2 does this for the Atwood machine, where $\lambda$ comes out as *minus* the tension because the coordinates are measured downward.

**Solving (3.18) in practice.** The system is a mix of differential equations (3.18) and algebraic ones ($f_k = 0$): a **differential-algebraic equation** (DAE). The standard way to integrate it is to differentiate each constraint twice along the motion. For scleronomic $f_k(q)$ and the case where $T = \frac12\dot q^\top M\dot q$ with constant $M$ (always true in Cartesian coordinates),
$$
\frac{d^2 f_k}{dt^2} = \sum_i\frac{\partial f_k}{\partial q_i}\ddot q_i + \sum_{ij}\frac{\partial^2 f_k}{\partial q_i\partial q_j}\dot q_i\dot q_j = 0,
$$
and (3.18) is $M\ddot q = F + J^\top\lambda$ with $J_{ki} = \partial f_k/\partial q_i$ and $F_i = Q_i^{\text{nc}} - \partial V/\partial q_i$. Together:
$$
\begin{pmatrix} M & -J^\top \\ J & 0\end{pmatrix}\begin{pmatrix}\ddot q \\ \lambda\end{pmatrix} = \begin{pmatrix} F \\ -\dot q^\top(\nabla^2 f)\dot q\end{pmatrix}. \tag{3.21}
$$
At each time step, solve this linear system for $\ddot q$ and $\lambda$, then advance $q$ and $\dot q$. The project does this for the Atwood machine and checks $\lambda$ against the hand-computed tension. The block matrix in (3.21) is the **KKT matrix** of ML-2, where the same structure (a positive definite block, a constraint Jacobian, a zero block) appears in equality-constrained quadratic programming. A warning: integrating (3.21) naively lets the constraint drift, since you enforce $\ddot f = 0$, not $f = 0$; the project measures the drift, and Exercise B2 asks you to fix it.

### 3.9 Dissipation: Rayleigh's function

Linear drag on each particle, $\mathbf{F}_a^{\text{nc}} = -c_a\dot{\mathbf{r}}_a$ with $c_a > 0$, is the one non-conservative force with a clean Lagrangian bookkeeping.

▶ **Try it first.** Compute $Q_i^{\text{nc}}$ from (3.6) for this force. Use the cancellation of dots (3.7) to turn $\partial\mathbf{r}_a/\partial q_i$ into $\partial\dot{\mathbf{r}}_a/\partial\dot q_i$, and recognize a derivative.

**Derivation.**

1. $Q_i^{\text{nc}} = -\sum_a c_a\dot{\mathbf{r}}_a\cdot\dfrac{\partial\mathbf{r}_a}{\partial q_i} \overset{(3.7)}{=} -\sum_a c_a\dot{\mathbf{r}}_a\cdot\dfrac{\partial\dot{\mathbf{r}}_a}{\partial\dot q_i}$.

2. The right side is $-\dfrac{\partial}{\partial\dot q_i}\left(\dfrac12\sum_a c_a\,\dot{\mathbf{r}}_a\cdot\dot{\mathbf{r}}_a\right)$. Define **Rayleigh's dissipation function**
   $$\mathcal{F}(q, \dot q, t) := \frac12\sum_a c_a|\dot{\mathbf{r}}_a|^2, \qquad\text{so that}\qquad Q_i^{\text{nc}} = -\frac{\partial\mathcal{F}}{\partial\dot q_i}. \tag{3.22}$$

3. Exactly as in §3.6, substituting $\dot{\mathbf{r}}_a = \sum_j(\partial\mathbf{r}_a/\partial q_j)\dot q_j$ (scleronomic case) gives
   $$\mathcal{F} = \frac12\sum_{ij}c_{ij}(q)\dot q_i\dot q_j, \qquad c_{ij} = \sum_a c_a\frac{\partial\mathbf{r}_a}{\partial q_i}\cdot\frac{\partial\mathbf{r}_a}{\partial q_j},$$
   a symmetric positive semidefinite quadratic form in the velocities, by the same Gram-matrix argument.

4. **Energy loss.** Put $Q_i^{\text{nc}} = -\partial\mathcal{F}/\partial\dot q_i$ into (3.16) with $\partial L/\partial t = 0$. Since $\mathcal{F}$ is homogeneous of degree 2 in $\dot q$, Euler's theorem gives $\sum_i\dot q_i\,\partial\mathcal{F}/\partial\dot q_i = 2\mathcal{F}$, so
   $$\frac{dE}{dt} = -2\mathcal{F}. \tag{3.23}$$
   ∎

■ **What just happened.** Drag enters Lagrange's equations as $\frac{d}{dt}\frac{\partial L}{\partial\dot q_i} - \frac{\partial L}{\partial q_i} + \frac{\partial\mathcal{F}}{\partial\dot q_i} = 0$: one extra function to differentiate, in the same coordinates as $L$, with no need to resolve drag forces along the constraint surface by hand. The energy drains at twice the rate $\mathcal{F}$; Rayleigh (*Theory of Sound*, 1877) defined $\mathcal{F}$ with the $\frac12$ so that it would be half the dissipated power, matching the $\frac12$ in $T$. Exercise B7 adds $\mathcal{F}$ to the double pendulum and checks (3.23) numerically.

### 3.10 The procedure

Here is the method, boxed. Every worked example in §4 follows it in this order, and you should too, in writing, until it is automatic.

> **Lagrangian procedure for a constrained system**
>
> 1. **Count degrees of freedom.** $n = 3N - m$ (or $2N - m$ in a plane). Check: can you specify the configuration with $n$ numbers and no fewer?
> 2. **Choose generalized coordinates** $q_1, \dots, q_n$ on the constraint surface. Prefer coordinates in which the constraint is *identically* satisfied (angles for pendulums, arc length for wires, distance along the incline for a block). Check they are non-degenerate where the motion goes.
> 3. **Write every position** $\mathbf{r}_a(q, t)$ explicitly. Include any imposed time dependence (a driven pivot, a spinning hoop) here.
> 4. **Compute $T$** by differentiating the $\mathbf{r}_a$ and summing $\frac12 m_a|\dot{\mathbf{r}}_a|^2$. Sort it as $T_2 + T_1 + T_0$ and read off the mass matrix $M_{ij}$. Check $M$ is symmetric and positive definite.
> 5. **Compute $V$** from the positions. Drop additive constants and pure functions of $t$ (they do not affect the equations).
> 6. **Form $L = T - V$.** If there are non-conservative forces, compute $Q_i^{\text{nc}}$ from (3.6) or $\mathcal{F}$ from (3.22). If you want a constraint force, keep the corresponding coordinate and add $\lambda\,\partial f/\partial q_i$.
> 7. **Write the Euler–Lagrange equations**, one per coordinate. Expand $\frac{d}{dt}\frac{\partial L}{\partial\dot q_i}$ fully (watch for $q$-dependence inside $M$: it produces $\dot q^2$ terms).
> 8. **Check limits.** Set masses to zero or infinity, angles to zero, $\omega$ to zero; each should reduce to something you already know. Check that cyclic coordinates give conserved momenta. Check units. Then check the whole thing with SymPy.

## 4. Worked examples

Eight systems, each run through the procedure of §3.10 with every step shown. The project (`project.py`) derives each Lagrangian's equations of motion symbolically and checks every number quoted here.

### Example 4.1: Simple pendulum (calibration)

Bob of mass $m$ on a massless rod of length $\ell$, pivot fixed. **DOF:** $2 - 1 = 1$. **Coordinate:** $\theta$ from the downward vertical. **Position:** $\mathbf{r} = (\ell\sin\theta, -\ell\cos\theta)$. **Velocity:** $\dot{\mathbf{r}} = \ell\dot\theta(\cos\theta, \sin\theta)$, so $T = \frac12 m\ell^2\dot\theta^2$ and $M = m\ell^2$ (the moment of inertia). **Potential:** $V = mgy = -mg\ell\cos\theta$.
$$
L = \tfrac12 m\ell^2\dot\theta^2 + mg\ell\cos\theta, \qquad \frac{d}{dt}(m\ell^2\dot\theta) = -mg\ell\sin\theta \quad\Longrightarrow\quad \ddot\theta = -\frac{g}{\ell}\sin\theta. \tag{4.1}
$$
The rod tension never appeared. **Limits:** small $\theta$ gives $\ddot\theta = -(g/\ell)\theta$, frequency $\omega_0 = \sqrt{g/\ell}$. One line, as promised. Everything else in this section is this line with more coordinates.

### Example 4.2: The Atwood machine, with and without the tension

Masses $m_1$ and $m_2$ hang from a massless, inextensible string of length $\ell$ over a massless frictionless pulley. Measure $x_1$ and $x_2$ **downward** from the pulley axle to each mass.

**Without multipliers.** **DOF:** 1 (the string length ties $x_1 + x_2 = \ell$, up to the constant half-circumference, which we absorb into $\ell$). **Coordinate:** $x := x_1$, so $x_2 = \ell - x$. **Velocities:** $\dot x_1 = \dot x$, $\dot x_2 = -\dot x$. $T = \frac12(m_1 + m_2)\dot x^2$. $V = -m_1 g x_1 - m_2 g x_2 = -m_1 g x - m_2 g(\ell - x) = -(m_1 - m_2)gx + \text{const}$.
$$
L = \tfrac12(m_1 + m_2)\dot x^2 + (m_1 - m_2)gx, \qquad (m_1 + m_2)\ddot x = (m_1 - m_2)g, \qquad \ddot x = \frac{m_1 - m_2}{m_1 + m_2}\,g. \tag{4.2}
$$
**Limits:** $m_1 = m_2$ gives no acceleration; $m_2 = 0$ gives free fall. The tension is invisible.

**With multipliers.** Keep both coordinates $x_1, x_2$ and the constraint $f(x_1, x_2) = x_1 + x_2 - \ell = 0$, so $\partial f/\partial x_1 = \partial f/\partial x_2 = 1$. The unconstrained Lagrangian is
$$
L = \tfrac12 m_1\dot x_1^2 + \tfrac12 m_2\dot x_2^2 + m_1 g x_1 + m_2 g x_2.
$$
Equations (3.18):
$$
m_1\ddot x_1 - m_1 g = \lambda, \qquad m_2\ddot x_2 - m_2 g = \lambda, \qquad x_1 + x_2 = \ell. \tag{4.3}
$$
Differentiate the constraint twice: $\ddot x_2 = -\ddot x_1$. Then from the two Lagrange equations, $\ddot x_1 = g + \lambda/m_1$ and $-\ddot x_1 = g + \lambda/m_2$. Add: $0 = 2g + \lambda(1/m_1 + 1/m_2)$, so
$$
\lambda = -\frac{2m_1 m_2}{m_1 + m_2}\,g, \qquad \ddot x_1 = g + \frac{\lambda}{m_1} = g - \frac{2m_2 g}{m_1 + m_2} = \frac{m_1 - m_2}{m_1 + m_2}\,g. \tag{4.4}
$$
Same acceleration as (4.2), and now a value for $\lambda$. **Interpretation via (3.20):** the constraint force on $m_1$ is $\lambda\,\partial f/\partial x_1 = \lambda$ in the $+x_1$ direction, which is downward. Since $\lambda < 0$, the string pulls *up* on $m_1$ with magnitude
$$
T_{\text{string}} = -\lambda = \frac{2m_1m_2}{m_1 + m_2}\,g.
$$
Check against Newton: $m_1\ddot x_1 = m_1 g - T$ gives $T = m_1(g - \ddot x_1) = m_1 g\cdot\frac{2m_2}{m_1 + m_2}$. ✓ **Numbers:** $m_1 = 3$ kg, $m_2 = 2$ kg: $\ddot x_1 = 9.81/5 = 1.962\ \mathrm{m/s^2}$, $T = 2\cdot3\cdot2\cdot9.81/5 = 23.544$ N, $\lambda = -23.544$ N. The project integrates (4.3) as the DAE (3.21) and asserts that the computed $\lambda$ equals $-23.544$ N at every step.

**The sign lesson.** $\lambda$ is negative because $f$ was written so that $\nabla f$ points along increasing $x_1, x_2$ (downward on both masses) and the string pulls the other way. Write $f = \ell - x_1 - x_2$ and $\lambda = +T$. Decide which way $\nabla f$ points *before* you interpret $\lambda$.

### Example 4.3: Bead on a rotating hoop (the showpiece)

A bead of mass $m$ slides without friction on a circular hoop of radius $R$ in a vertical plane. The hoop is driven to rotate about its vertical diameter at constant angular speed $\omega$. **DOF:** the bead is on a curve: $1$. **Coordinate:** $\theta$, the angle from the bottom of the hoop. **Position** (the hoop's plane makes azimuth $\phi = \omega t$ with the $x$-axis; this is the imposed time dependence):
$$
\mathbf{r} = \big(R\sin\theta\cos\omega t,\; R\sin\theta\sin\omega t,\; -R\cos\theta\big).
$$
**Velocity.** $\dot{\mathbf{r}} = R\dot\theta\,(\cos\theta\cos\omega t, \cos\theta\sin\omega t, \sin\theta) + R\omega\sin\theta\,(-\sin\omega t, \cos\omega t, 0)$. The two vectors are orthogonal (first is tangent to the hoop, second is azimuthal), with lengths $R\dot\theta$ and $R\omega\sin\theta$. So
$$
T = \tfrac12 mR^2\dot\theta^2 + \tfrac12 mR^2\omega^2\sin^2\theta, \qquad T_2 = \tfrac12 mR^2\dot\theta^2,\ T_1 = 0,\ T_0 = \tfrac12 mR^2\omega^2\sin^2\theta.
$$
**Potential:** $V = mgz = -mgR\cos\theta$.
$$
L = \tfrac12 mR^2\dot\theta^2 + \tfrac12 mR^2\omega^2\sin^2\theta + mgR\cos\theta. \tag{4.5}
$$
**Euler–Lagrange.** $\partial L/\partial\dot\theta = mR^2\dot\theta$; $\partial L/\partial\theta = mR^2\omega^2\sin\theta\cos\theta - mgR\sin\theta$. So
$$
\ddot\theta = \sin\theta\left(\omega^2\cos\theta - \frac{g}{R}\right). \tag{4.6}
$$
**Limits:** $\omega = 0$ gives the pendulum (4.1) with $\ell = R$. ✓

**Effective potential.** (4.6) has the form $mR^2\ddot\theta = -dV_{\text{eff}}/d\theta$ with
$$
V_{\text{eff}}(\theta) = -mgR\cos\theta - \tfrac12 mR^2\omega^2\sin^2\theta = V - T_0. \tag{4.7}
$$
Check: $-V_{\text{eff}}' = -mgR\sin\theta + mR^2\omega^2\sin\theta\cos\theta$. ✓ The $-T_0$ term is the centrifugal potential: the hoop's rotation acts on the bead's $\theta$-motion like an extra potential energy pulling it toward the rim ($\theta = \pi/2$). The conserved energy function (3.15) is $h = T_2 + V_{\text{eff}}$, **not** $T + V$. $L$ has no explicit $t$, so $h$ is conserved; $T + V = h + 2T_0$ is not, because the motor driving the hoop does real work on the bead. The project integrates (4.6) and checks: $h$ constant to $10^{-9}$, $T + V$ not constant.

**Equilibria.** $V_{\text{eff}}'(\theta) = \sin\theta\,(mgR - mR^2\omega^2\cos\theta) = 0$:
- $\sin\theta = 0$: $\theta = 0$ (bottom) and $\theta = \pi$ (top), for every $\omega$.
- $\cos\theta^* = \dfrac{g}{R\omega^2} = \dfrac{\omega_0^2}{\omega^2}$ with $\omega_0 := \sqrt{g/R}$. This requires $\omega_0^2/\omega^2 \le 1$, i.e. $\omega \ge \omega_0$. For $\omega > \omega_0$ there are two such points, $\pm\theta^*$, symmetric about the bottom.

**Stability.** $V_{\text{eff}}'' = mgR\cos\theta - mR^2\omega^2\cos 2\theta$. A stable equilibrium has $V_{\text{eff}}'' > 0$, and the frequency of small oscillations about it is $\Omega^2 = V_{\text{eff}}''/(mR^2)$ (since $mR^2\ddot{\epsilon} = -V_{\text{eff}}''\epsilon$ for $\theta = \theta_{\text{eq}} + \epsilon$).

1. **Bottom, $\theta = 0$:** $V_{\text{eff}}'' = mgR - mR^2\omega^2 = mR^2(\omega_0^2 - \omega^2)$. Stable iff $\omega < \omega_0$, with $\Omega^2 = \omega_0^2 - \omega^2$. The rotation *softens* the pendulum; at $\omega = \omega_0$ the restoring force vanishes and above it the bottom is a hilltop.
2. **Top, $\theta = \pi$:** $V_{\text{eff}}'' = -mgR - mR^2\omega^2 < 0$. Unstable for every $\omega$.
3. **Off-axis, $\cos\theta^* = \omega_0^2/\omega^2$:** use $\cos 2\theta = 2\cos^2\theta - 1$:
   $$V_{\text{eff}}''(\theta^*) = mgR\cos\theta^* - mR^2\omega^2(2\cos^2\theta^* - 1) = mR^2\omega^2\cos^2\theta^* - 2mR^2\omega^2\cos^2\theta^* + mR^2\omega^2 = mR^2\omega^2\sin^2\theta^* > 0,$$
   where $mgR\cos\theta^* = mR^2\omega^2\cos^2\theta^*$ by the equilibrium condition. Stable whenever it exists, with
   $$\Omega^2 = \omega^2\sin^2\theta^* = \omega^2\left(1 - \frac{\omega_0^4}{\omega^4}\right) = \omega^2 - \frac{\omega_0^4}{\omega^2}. \tag{4.8}$$

**The pitchfork.** Plot the equilibria against $\omega$: for $\omega < \omega_0$ one stable branch at $\theta = 0$; at $\omega = \omega_0$ it loses stability and two stable branches $\pm\arccos(\omega_0^2/\omega^2)$ grow out of it, approaching $\pm\pi/2$ as $\omega\to\infty$. The picture is a fork with a handle, hence the name **supercritical pitchfork bifurcation**. The symmetry $\theta\to-\theta$ of $L$ is what makes the two new branches appear as a pair; a system without that symmetry would bifurcate differently. The project draws the diagram from these formulas and confirms stability on each branch by integration.

**Numerical example.** Take $R = 0.5$ m, $g = 9.81\ \mathrm{m/s^2}$, so $\omega_0 = \sqrt{19.62} = 4.4294$ rad/s, and spin at $\omega = 2\omega_0 = 8.8589$ rad/s. Then $\cos\theta^* = 1/4$, $\theta^* = \arccos(0.25) = 1.3181$ rad $= 75.52°$, $\sin\theta^* = \sqrt{15}/4 = 0.96825$, and
$$
\Omega = \omega\sin\theta^* = \frac{\sqrt{15}}{2}\,\omega_0 = 1.9365\times4.4294 = 8.5776\ \mathrm{rad/s}, \qquad \text{period } \frac{2\pi}{\Omega} = 0.7325\ \text{s}.
$$
For comparison, the stationary hoop's pendulum period is $2\pi/\omega_0 = 1.4185$ s. The project perturbs the bead by $0.02$ rad from $\theta^*$, integrates (4.6), measures the period from zero crossings of $\theta - \theta^*$, and asserts agreement with $0.7325$ s to better than $0.5\%$ (the residual is the amplitude dependence of a nonlinear oscillator).

**Where d'Alembert is slippery here.** The hoop exerts a normal force on the bead with a component along $\hat{\boldsymbol\phi}$ (that is what drags the bead around). Under a *real* displacement $d\mathbf{r}$, which includes $R\omega\sin\theta\,dt\,\hat{\boldsymbol\phi}$, this force does work: the bead's $T + V$ changes. Under a *virtual* displacement, $\delta\mathbf{r} = (\partial\mathbf{r}/\partial\theta)\delta\theta$ at frozen $t$, the displacement is tangent to the hoop and the normal force does no work. Principle 3.5 holds, the derivation goes through, and (4.6) is correct. If you let "virtual" slide into "real," you will think d'Alembert fails for moving constraints. It does not; the frozen time is the whole point.

### Example 4.4: The double pendulum

Bob 1 (mass $m_1$) hangs from a fixed pivot on a rod of length $\ell_1$; bob 2 (mass $m_2$) hangs from bob 1 on a rod of length $\ell_2$. Planar motion. **DOF:** $2\cdot2 - 2 = 2$. **Coordinates:** $\theta_1, \theta_2$, each from the downward vertical. **Positions:**
$$
x_1 = \ell_1\sin\theta_1,\quad y_1 = -\ell_1\cos\theta_1,\qquad x_2 = x_1 + \ell_2\sin\theta_2,\quad y_2 = y_1 - \ell_2\cos\theta_2.
$$
**Velocities:** $\dot x_1 = \ell_1\dot\theta_1\cos\theta_1$, $\dot y_1 = \ell_1\dot\theta_1\sin\theta_1$; $\dot x_2 = \dot x_1 + \ell_2\dot\theta_2\cos\theta_2$, $\dot y_2 = \dot y_1 + \ell_2\dot\theta_2\sin\theta_2$. Then $|\dot{\mathbf{r}}_1|^2 = \ell_1^2\dot\theta_1^2$ and
$$
|\dot{\mathbf{r}}_2|^2 = \ell_1^2\dot\theta_1^2 + \ell_2^2\dot\theta_2^2 + 2\ell_1\ell_2\dot\theta_1\dot\theta_2(\cos\theta_1\cos\theta_2 + \sin\theta_1\sin\theta_2) = \ell_1^2\dot\theta_1^2 + \ell_2^2\dot\theta_2^2 + 2\ell_1\ell_2\dot\theta_1\dot\theta_2\cos(\theta_1 - \theta_2).
$$
**Kinetic and potential energy.** Writing $\Delta := \theta_1 - \theta_2$,
$$
T = \tfrac12(m_1 + m_2)\ell_1^2\dot\theta_1^2 + \tfrac12 m_2\ell_2^2\dot\theta_2^2 + m_2\ell_1\ell_2\dot\theta_1\dot\theta_2\cos\Delta, \qquad V = -(m_1 + m_2)g\ell_1\cos\theta_1 - m_2 g\ell_2\cos\theta_2. \tag{4.9}
$$
**Mass matrix:** $M = \begin{pmatrix}(m_1 + m_2)\ell_1^2 & m_2\ell_1\ell_2\cos\Delta\\ m_2\ell_1\ell_2\cos\Delta & m_2\ell_2^2\end{pmatrix}$, with $\det M = m_2\ell_1^2\ell_2^2\left[(m_1 + m_2) - m_2\cos^2\Delta\right] = m_2\ell_1^2\ell_2^2(m_1 + m_2\sin^2\Delta) > 0$. Positive definite everywhere, as §3.6 promised: these are good coordinates.

**Euler–Lagrange for $\theta_1$.** $\partial L/\partial\dot\theta_1 = (m_1 + m_2)\ell_1^2\dot\theta_1 + m_2\ell_1\ell_2\dot\theta_2\cos\Delta$. Its time derivative (note $\dot\Delta = \dot\theta_1 - \dot\theta_2$):
$$
(m_1 + m_2)\ell_1^2\ddot\theta_1 + m_2\ell_1\ell_2\ddot\theta_2\cos\Delta - m_2\ell_1\ell_2\dot\theta_2(\dot\theta_1 - \dot\theta_2)\sin\Delta.
$$
$\partial L/\partial\theta_1 = -m_2\ell_1\ell_2\dot\theta_1\dot\theta_2\sin\Delta - (m_1 + m_2)g\ell_1\sin\theta_1$. Subtracting, the $\dot\theta_1\dot\theta_2\sin\Delta$ terms cancel ($-m_2\ell_1\ell_2\dot\theta_1\dot\theta_2\sin\Delta$ from the first, $+m_2\ell_1\ell_2\dot\theta_1\dot\theta_2\sin\Delta$ from the second), leaving, after dividing by $\ell_1$,
$$
(m_1 + m_2)\ell_1\ddot\theta_1 + m_2\ell_2\ddot\theta_2\cos\Delta + m_2\ell_2\dot\theta_2^2\sin\Delta + (m_1 + m_2)g\sin\theta_1 = 0. \tag{4.10a}
$$
**Euler–Lagrange for $\theta_2$.** $\partial L/\partial\dot\theta_2 = m_2\ell_2^2\dot\theta_2 + m_2\ell_1\ell_2\dot\theta_1\cos\Delta$; time derivative $m_2\ell_2^2\ddot\theta_2 + m_2\ell_1\ell_2\ddot\theta_1\cos\Delta - m_2\ell_1\ell_2\dot\theta_1(\dot\theta_1 - \dot\theta_2)\sin\Delta$. $\partial L/\partial\theta_2 = +m_2\ell_1\ell_2\dot\theta_1\dot\theta_2\sin\Delta - m_2 g\ell_2\sin\theta_2$ (the sign flips because $\partial\Delta/\partial\theta_2 = -1$). Again the $\dot\theta_1\dot\theta_2$ terms cancel; divide by $m_2\ell_2$:
$$
\ell_2\ddot\theta_2 + \ell_1\ddot\theta_1\cos\Delta - \ell_1\dot\theta_1^2\sin\Delta + g\sin\theta_2 = 0. \tag{4.10b}
$$
(4.10) is a linear system for $(\ddot\theta_1, \ddot\theta_2)$ with matrix $M$ (divided through by $\ell_1$ and $m_2\ell_2$); solving it by hand produces half a page of trigonometry. The project lets SymPy do that and asserts the result agrees with (4.10) symbolically. **Limits:** $m_2 \to 0$ in (4.10a) gives the simple pendulum. $\theta_1 = \theta_2$ identically gives a rigid rod of length $\ell_1 + \ell_2$ with two masses, and (4.10) collapses correctly (check it).

**Small oscillations (CM-6 preview).** Take $m_1 = m_2 = m$, $\ell_1 = \ell_2 = \ell$, and linearize about $\theta_1 = \theta_2 = 0$: $\cos\Delta\to1$, $\sin\theta_i\to\theta_i$, and the $\dot\theta^2\sin\Delta$ terms are second order and drop. (4.10) becomes
$$
M_0\begin{pmatrix}\ddot\theta_1\\\ddot\theta_2\end{pmatrix} + K\begin{pmatrix}\theta_1\\\theta_2\end{pmatrix} = 0, \qquad M_0 = m\ell^2\begin{pmatrix}2 & 1\\1 & 1\end{pmatrix}, \quad K = mg\ell\begin{pmatrix}2 & 0\\0 & 1\end{pmatrix}.
$$
($M_0$ is the mass matrix at equilibrium and $K_{ij} = \partial^2V/\partial\theta_i\partial\theta_j$ at equilibrium.) Try $\boldsymbol\theta = \mathbf{a}\cos\Omega t$: $(K - \Omega^2M_0)\mathbf{a} = 0$, the **generalized eigenproblem** of CM-6. With $\mu := \Omega^2\ell/g$,
$$
\det\begin{pmatrix}2 - 2\mu & -\mu\\ -\mu & 1 - \mu\end{pmatrix} = (2 - 2\mu)(1 - \mu) - \mu^2 = \mu^2 - 4\mu + 2 = 0 \quad\Longrightarrow\quad \mu = 2 \mp \sqrt2.
$$
So $\Omega_\pm^2 = (2\mp\sqrt2)\,g/\ell$. The mode shapes come from the first row, $a_2/a_1 = (2 - 2\mu)/\mu$: for $\mu = 2 - \sqrt2$, $a_2/a_1 = (2\sqrt2 - 2)/(2 - \sqrt2) = \sqrt2$ (both bobs swing the same way, the lower one $\sqrt2$ times further); for $\mu = 2 + \sqrt2$, $a_2/a_1 = -\sqrt2$ (opposite directions). **Numbers** ($m = 1$, $\ell = 1$, $g = 9.81$): $\Omega_-^2 = 0.585786\times9.81 = 5.7466$, $\Omega_- = 2.3972$ rad/s; $\Omega_+^2 = 3.414214\times9.81 = 33.4934$, $\Omega_+ = 5.7874$ rad/s. The project computes $M_0$ and $K$ as SymPy Hessians and solves the eigenproblem with `scipy.linalg.eigh(K, M0)`.

**Chaos.** At large amplitude the double pendulum is chaotic: two trajectories that start $10^{-8}$ rad apart separate exponentially, $\delta(t)\sim\delta_0 e^{\lambda_L t}$, until the separation saturates at the size of the accessible region. The exponent $\lambda_L$ is the **Lyapunov exponent**. There is no closed-form solution to find; the statement is numerical. The project integrates two such trajectories from $\theta_1 = \theta_2 = 2$ rad at rest, plots $\log\delta$ against $t$, and fits the slope: the estimate is about $1.2\ \mathrm{s^{-1}}$ for $m = \ell = 1$, meaning the trajectories lose one decimal digit of agreement every $\ln 10/1.2 \approx 1.9$ s. After 30 s, $10^{-8}$ has become order 1. A snapshot grid in the project shows the two pendulums agreeing, then not.

### Example 4.5: Pendulum with a horizontally oscillating support

The pivot of a pendulum (mass $m$, length $\ell$) is driven along the $x$-axis: $x_p(t) = a\cos\gamma t$. **DOF:** 1 (the pivot's motion is imposed, not free). **Coordinate:** $\theta$. **Position:** $x = a\cos\gamma t + \ell\sin\theta$, $y = -\ell\cos\theta$. **Velocity:** $\dot x = -a\gamma\sin\gamma t + \ell\dot\theta\cos\theta$, $\dot y = \ell\dot\theta\sin\theta$.
$$
T = \tfrac12 m\left(\ell^2\dot\theta^2 - 2a\gamma\ell\,\dot\theta\sin\gamma t\cos\theta + a^2\gamma^2\sin^2\gamma t\right), \qquad V = -mg\ell\cos\theta.
$$
Here $T_2 = \frac12 m\ell^2\dot\theta^2$, $T_1 = -ma\gamma\ell\dot\theta\sin\gamma t\cos\theta$ (a rheonomic cross term), and $T_0 = \frac12 ma^2\gamma^2\sin^2\gamma t$, a pure function of $t$ that drops out of the equations (step 5 of the procedure).

**Euler–Lagrange.** $\partial L/\partial\dot\theta = m\ell^2\dot\theta - ma\gamma\ell\sin\gamma t\cos\theta$. Time derivative:
$$
m\ell^2\ddot\theta - ma\gamma^2\ell\cos\gamma t\cos\theta + ma\gamma\ell\,\dot\theta\sin\gamma t\sin\theta.
$$
$\partial L/\partial\theta = ma\gamma\ell\,\dot\theta\sin\gamma t\sin\theta - mg\ell\sin\theta$. Subtract; the $\dot\theta\sin\gamma t\sin\theta$ terms cancel:
$$
\ddot\theta + \frac{g}{\ell}\sin\theta = \frac{a\gamma^2}{\ell}\cos\gamma t\,\cos\theta. \tag{4.11}
$$
**Interpretation.** $-\ddot x_p = a\gamma^2\cos\gamma t$ is the pivot's acceleration reversed. In the pivot's frame there is an effective horizontal gravity $-\ddot x_p$, and (4.11) is the pendulum equation with total effective gravity $(−\ddot x_p, -g)$ projected onto the tangent direction. **Limit:** $a = 0$ gives (4.1). ✓ The drive multiplies $\cos\theta$, a state-dependent coefficient, which is what "parametric" means; at small angles $\cos\theta\approx1$ and (4.11) is an ordinary driven oscillator, resonant at $\gamma = \omega_0$. If instead the pivot is shaken **vertically**, the drive multiplies $\sin\theta$, the linearization is $\ddot\theta + (\omega_0^2 - (a\gamma^2/\ell)\cos\gamma t)\theta = 0$, the **Mathieu equation**, and its stability chart contains Kapitza's result that fast vertical shaking stabilizes the inverted pendulum. Exercise C2 derives that equation and tests Kapitza's criterion numerically. The project checks (4.11) against SymPy.

### Example 4.6: Block on a sliding wedge

A wedge of mass $M$ and angle $\alpha$ sits on a frictionless floor. A block of mass $m$ slides without friction down its incline. **DOF:** wedge in 1D ($1$), block on the incline relative to the wedge ($1$): $n = 2$. **Coordinates:** $X$, the wedge's horizontal position; $s$, the block's distance down the incline from the wedge's top. **Positions:** wedge $(X, 0)$; block $(X + s\cos\alpha,\; h - s\sin\alpha)$ with $h$ the wedge height. **Velocities:** $\dot X$; block $(\dot X + \dot s\cos\alpha, -\dot s\sin\alpha)$.
$$
T = \tfrac12 M\dot X^2 + \tfrac12 m\left[(\dot X + \dot s\cos\alpha)^2 + \dot s^2\sin^2\alpha\right] = \tfrac12(M + m)\dot X^2 + m\cos\alpha\,\dot X\dot s + \tfrac12 m\dot s^2,
$$
$$
V = mg(h - s\sin\alpha) = -mgs\sin\alpha + \text{const}.
$$
**Mass matrix:** $\begin{pmatrix}M + m & m\cos\alpha\\ m\cos\alpha & m\end{pmatrix}$, determinant $m(M + m) - m^2\cos^2\alpha = m(M + m\sin^2\alpha) > 0$. (Exercise A4 writes it as $J^\top\mathcal{M}J$.)

**Euler–Lagrange.** $X$ is cyclic ($L$ does not contain $X$), so
$$
p_X = (M + m)\dot X + m\cos\alpha\,\dot s = \text{const}, \tag{4.12}
$$
which is conservation of total horizontal momentum (the floor pushes only vertically). Differentiating, $(M + m)\ddot X + m\cos\alpha\,\ddot s = 0$. For $s$: $\frac{d}{dt}(m\dot s + m\cos\alpha\,\dot X) = mg\sin\alpha$, i.e. $\ddot s + \cos\alpha\,\ddot X = g\sin\alpha$. Substitute $\ddot X = -\frac{m\cos\alpha}{M + m}\ddot s$:
$$
\ddot s\left(1 - \frac{m\cos^2\alpha}{M + m}\right) = g\sin\alpha \quad\Longrightarrow\quad \ddot s = \frac{(M + m)\,g\sin\alpha}{M + m\sin^2\alpha}, \qquad \ddot X = -\frac{m\,g\sin\alpha\cos\alpha}{M + m\sin^2\alpha}. \tag{4.13}
$$
**Limits:** $M\to\infty$: $\ddot s\to g\sin\alpha$, $\ddot X\to0$, the fixed incline. ✓ $\alpha\to\pi/2$: $\ddot s\to g$ (free fall beside a wall that is not pushed). ✓ $\alpha\to0$: nothing moves. ✓ Note $\ddot s > g\sin\alpha$ when the wedge can move: the incline retreats from under the block. **Numbers** ($M = 4$ kg, $m = 1$ kg, $\alpha = 30°$): $\ddot s = 5\times9.81\times0.5/(4 + 0.25) = 5.7706\ \mathrm{m/s^2}$, $\ddot X = -9.81\times0.5\times0.86603/4.25 = -0.9995\ \mathrm{m/s^2}$. Both constant: the motion is uniformly accelerated, and the project checks that the integrated $s(t)$ is exactly quadratic, that (4.12) holds, and that energy is conserved.

### Example 4.7: The spherical pendulum (CM-5 preview)

A bob of mass $m$ on a rigid massless rod of length $\ell$, free to swing in any direction. **DOF:** $3 - 1 = 2$. **Coordinates:** polar angle $\theta$ from the downward vertical and azimuth $\phi$. **Position:** $\mathbf{r} = \ell(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ -\cos\theta)$. **Velocity:** $\dot{\mathbf{r}} = \ell\dot\theta\,\hat{\boldsymbol\theta} + \ell\sin\theta\,\dot\phi\,\hat{\boldsymbol\phi}$ with orthogonal unit vectors, so
$$
T = \tfrac12 m\ell^2\left(\dot\theta^2 + \sin^2\theta\,\dot\phi^2\right), \qquad V = -mg\ell\cos\theta, \qquad L = \tfrac12 m\ell^2(\dot\theta^2 + \sin^2\theta\,\dot\phi^2) + mg\ell\cos\theta. \tag{4.14}
$$
**Mass matrix:** $m\ell^2\operatorname{diag}(1, \sin^2\theta)$: degenerate at $\theta = 0, \pi$, where $\phi$ is undefined. These coordinates are good everywhere else.

**Euler–Lagrange.** $\phi$ is cyclic:
$$
p_\phi = m\ell^2\sin^2\theta\,\dot\phi = \text{const}, \tag{4.15}
$$
the vertical component of angular momentum. For $\theta$: $m\ell^2\ddot\theta = m\ell^2\sin\theta\cos\theta\,\dot\phi^2 - mg\ell\sin\theta$.

**Reduction to one degree of freedom.** Eliminate $\dot\phi$ using (4.15), $\dot\phi = p_\phi/(m\ell^2\sin^2\theta)$:
$$
m\ell^2\ddot\theta = \frac{p_\phi^2\cos\theta}{m\ell^2\sin^3\theta} - mg\ell\sin\theta = -\frac{d}{d\theta}V_{\text{eff}}(\theta), \qquad V_{\text{eff}}(\theta) = \frac{p_\phi^2}{2m\ell^2\sin^2\theta} - mg\ell\cos\theta. \tag{4.16}
$$
(Check: $\frac{d}{d\theta}\frac{p_\phi^2}{2m\ell^2\sin^2\theta} = -\frac{p_\phi^2\cos\theta}{m\ell^2\sin^3\theta}$. ✓) The $\theta$-motion is a one-dimensional problem in the effective potential (4.16), whose first term, the **centrifugal barrier**, blows up at $\theta\to0$: a bob with any angular momentum can never pass through the bottom. The energy $E = \frac12 m\ell^2\dot\theta^2 + V_{\text{eff}}$ is conserved and equals the true $T + V$ here (scleronomic constraint, so $h = E$). This is exactly the reduction CM-5 performs for a planet: conserve angular momentum, fold it into an effective potential, solve a one-dimensional problem.

**Conical pendulum.** $V_{\text{eff}}'(\theta_0) = 0$ gives $\frac{p_\phi^2\cos\theta_0}{m\ell^2\sin^3\theta_0} = mg\ell\sin\theta_0$, i.e. with (4.15), $\dot\phi_0^2 = \dfrac{g}{\ell\cos\theta_0}$: the bob circles at constant $\theta_0$ and constant rate. Exercise B6 shows the small oscillations of $\theta$ about $\theta_0$ have frequency $\Omega = \dot\phi_0\sqrt{1 + 3\cos^2\theta_0}$, which is between $\dot\phi_0$ and $2\dot\phi_0$, so the orbit does not close in general: it precesses. The project checks (4.15), energy conservation, the conical solution, and this frequency.

**A trap, flagged here and again in §8.** Do not substitute (4.15) into $L$ and vary the result. $\frac12 m\ell^2\sin^2\theta\,\dot\phi^2$ becomes $+\frac{p_\phi^2}{2m\ell^2\sin^2\theta}$ *in $L$*, where it acts like $-V$: the centrifugal term comes out with the wrong sign. Substitute into the **equation of motion** (or into $E$), as done above. The object you may legitimately substitute into is the Routhian $L - p_\phi\dot\phi$, which CM-8 introduces alongside the Hamiltonian.

### Example 4.8: Charged particle in electric and magnetic fields (CM-2 recap, with the momentum highlighted)

Charge $q_e$, mass $m$, in fields given by potentials $\Phi(\mathbf{r}, t)$ and $\mathbf{A}(\mathbf{r}, t)$ with $\mathbf{E} = -\nabla\Phi - \partial_t\mathbf{A}$, $\mathbf{B} = \nabla\times\mathbf{A}$. No constraints; coordinates $x_1, x_2, x_3$. CM-2 showed that
$$
L = \tfrac12 m|\dot{\mathbf{r}}|^2 - q_e\Phi + q_e\,\dot{\mathbf{r}}\cdot\mathbf{A} \tag{4.17}
$$
reproduces the Lorentz force. Here is the computation again, because the generalized momentum is the point.

**Generalized momentum:** $p_i = \partial L/\partial\dot x_i = m\dot x_i + q_eA_i$, i.e.
$$
\mathbf{p} = m\mathbf{v} + q_e\mathbf{A}. \tag{4.18}
$$
This is **not** $m\mathbf{v}$. It is what is conserved when a coordinate is cyclic, it is what becomes $-i\hbar\nabla$ in QM-2, and confusing it with $m\mathbf{v}$ is the single most common error in this topic.

**Euler–Lagrange.** $\frac{d}{dt}(m\dot x_i + q_eA_i) = \partial L/\partial x_i = -q_e\partial_i\Phi + q_e\sum_j\dot x_j\partial_iA_j$. Along the motion $\frac{dA_i}{dt} = \partial_tA_i + \sum_j\dot x_j\partial_jA_i$. So
$$
m\ddot x_i = -q_e\partial_i\Phi - q_e\partial_tA_i + q_e\sum_j\dot x_j\left(\partial_iA_j - \partial_jA_i\right) = q_eE_i + q_e\,(\dot{\mathbf{r}}\times\mathbf{B})_i,
$$
using the identity $(\mathbf{v}\times(\nabla\times\mathbf{A}))_i = \sum_j v_j\partial_iA_j - \sum_j v_j\partial_jA_i$ (EM-1). ∎ The velocity-dependent potential of §3.5 is $U = q_e\Phi - q_e\dot{\mathbf{r}}\cdot\mathbf{A}$.

**Worked case.** Uniform $\mathbf{B} = B\hat{\mathbf{z}}$, $\Phi = 0$, in the **Landau gauge** $\mathbf{A} = (0, Bx, 0)$ (check: $\nabla\times\mathbf{A} = \partial_x A_y\,\hat{\mathbf{z}} = B\hat{\mathbf{z}}$). Then $L = \frac12 m(\dot x^2 + \dot y^2 + \dot z^2) + q_eBx\dot y$. Both $y$ and $z$ are cyclic: $p_y = m\dot y + q_eBx$ and $p_z = m\dot z$ are conserved. The $x$ equation: $m\ddot x = q_eB\dot y$. Combine with $\dot y = (p_y - q_eBx)/m$: $m\ddot x = \frac{q_eB}{m}(p_y - q_eBx)$, a harmonic oscillator about $x_c = p_y/(q_eB)$ with frequency $\omega_c = q_eB/m$, the **cyclotron frequency**. The project integrates this with $q_e = m = 1$, $B = 2$, $\mathbf{v}_0 = (1, 0, 0)$ and asserts $p_y$ is conserved while $m\dot y$ is not, that the orbit radius is $mv/(q_eB) = 0.5$, and that the period is $2\pi m/(q_eB) = \pi$.

### Example 4.9: Non-inertial frames from the Lagrangian

A particle of mass $m$ in a potential $V$, described by an observer in a frame rotating with angular velocity $\boldsymbol\omega(t)$ about the origin. Taylor ch. 9 derives the fictitious forces by differentiating vectors in a rotating basis twice. The Lagrangian route is shorter and does the bookkeeping automatically.

**Setup.** Let $\mathbf{r}$ be the inertial position and $\mathbf{r}'$ the same vector's components in the rotating frame: $\mathbf{r} = \mathcal{R}(t)\mathbf{r}'$ with $\mathcal{R}(t)$ a rotation matrix ($\mathcal{R}^\top\mathcal{R} = I$). The rotating observer measures $\mathbf{r}'$, $\dot{\mathbf{r}}'$, $\ddot{\mathbf{r}}'$.

▶ **Try it first.** Differentiate $\mathcal{R}^\top\mathcal{R} = I$ to show $\mathcal{R}^\top\dot{\mathcal{R}}$ is antisymmetric. An antisymmetric $3\times3$ matrix acts as a cross product with some vector; call it $\boldsymbol\omega$. Use this to express the inertial velocity in rotating components, then write $L$ and turn the crank.

**Derivation.**

1. **Angular velocity.** From $\mathcal{R}^\top\mathcal{R} = I$: $\dot{\mathcal{R}}^\top\mathcal{R} + \mathcal{R}^\top\dot{\mathcal{R}} = 0$, so $W := \mathcal{R}^\top\dot{\mathcal{R}}$ satisfies $W^\top = -W$. Every antisymmetric $3\times3$ matrix has the form $W\mathbf{u} = \boldsymbol\omega\times\mathbf{u}$ for a unique vector $\boldsymbol\omega$ (its three independent entries), which **defines** the angular velocity, in rotating-frame components. (LA-12's appendix on $SO(3)$ is this fact in general.)

2. **Inertial velocity.** $\dot{\mathbf{r}} = \dot{\mathcal{R}}\mathbf{r}' + \mathcal{R}\dot{\mathbf{r}}' = \mathcal{R}\left(\mathcal{R}^\top\dot{\mathcal{R}}\,\mathbf{r}' + \dot{\mathbf{r}}'\right) = \mathcal{R}\left(\dot{\mathbf{r}}' + \boldsymbol\omega\times\mathbf{r}'\right)$. Rotations preserve lengths, so
   $$|\dot{\mathbf{r}}|^2 = |\dot{\mathbf{r}}' + \boldsymbol\omega\times\mathbf{r}'|^2.$$

3. **The Lagrangian in rotating coordinates.** $V$ is a scalar, so $V(\mathbf{r}) = V(\mathcal{R}\mathbf{r}') =: V'(\mathbf{r}', t)$. Then
   $$L = \tfrac12 m|\dot{\mathbf{r}}' + \boldsymbol\omega\times\mathbf{r}'|^2 - V' = \tfrac12 m|\dot{\mathbf{r}}'|^2 + m\,\dot{\mathbf{r}}'\cdot(\boldsymbol\omega\times\mathbf{r}') + \tfrac12 m|\boldsymbol\omega\times\mathbf{r}'|^2 - V'. \tag{4.19}$$
   This is $T_2 + T_1 + T_0 - V$ in the notation of §3.6: the rotating coordinates are a rheonomic change of variables, and $T_1$ (the Coriolis term) and $T_0$ (the centrifugal term) are the $N_i$ and $T_0$ of (3.13).

4. **Momentum.** $\mathbf{p}' = \partial L/\partial\dot{\mathbf{r}}' = m\dot{\mathbf{r}}' + m\boldsymbol\omega\times\mathbf{r}'$ (this is $\mathcal{R}^\top$ applied to the inertial momentum). Its time derivative, keeping $\dot{\boldsymbol\omega}$:
   $$\frac{d\mathbf{p}'}{dt} = m\ddot{\mathbf{r}}' + m\dot{\boldsymbol\omega}\times\mathbf{r}' + m\boldsymbol\omega\times\dot{\mathbf{r}}'.$$

5. **Position derivative.** Two vector-calculus facts. First, $\dot{\mathbf{r}}'\cdot(\boldsymbol\omega\times\mathbf{r}') = \mathbf{r}'\cdot(\dot{\mathbf{r}}'\times\boldsymbol\omega)$ (cyclic permutation of the triple product), whose gradient with respect to $\mathbf{r}'$ (holding $\dot{\mathbf{r}}'$ fixed) is $\dot{\mathbf{r}}'\times\boldsymbol\omega$. Second, $|\boldsymbol\omega\times\mathbf{r}'|^2 = \omega^2|\mathbf{r}'|^2 - (\boldsymbol\omega\cdot\mathbf{r}')^2$, whose gradient is $2\omega^2\mathbf{r}' - 2(\boldsymbol\omega\cdot\mathbf{r}')\boldsymbol\omega = -2\,\boldsymbol\omega\times(\boldsymbol\omega\times\mathbf{r}')$, by the BAC–CAB identity $\boldsymbol\omega\times(\boldsymbol\omega\times\mathbf{r}') = (\boldsymbol\omega\cdot\mathbf{r}')\boldsymbol\omega - \omega^2\mathbf{r}'$. So
   $$\frac{\partial L}{\partial\mathbf{r}'} = m\,\dot{\mathbf{r}}'\times\boldsymbol\omega - m\,\boldsymbol\omega\times(\boldsymbol\omega\times\mathbf{r}') - \nabla' V'.$$

6. **Assemble** $d\mathbf{p}'/dt = \partial L/\partial\mathbf{r}'$ and use $\dot{\mathbf{r}}'\times\boldsymbol\omega = -\boldsymbol\omega\times\dot{\mathbf{r}}'$:
   $$\boxed{\;m\ddot{\mathbf{r}}' = -\nabla'V' \;\underbrace{-\;2m\,\boldsymbol\omega\times\dot{\mathbf{r}}'}_{\text{Coriolis}}\;\underbrace{-\;m\,\boldsymbol\omega\times(\boldsymbol\omega\times\mathbf{r}')}_{\text{centrifugal}}\;\underbrace{-\;m\,\dot{\boldsymbol\omega}\times\mathbf{r}'}_{\text{Euler}}\;} \tag{4.20}$$
   ∎

■ **What just happened.** The three fictitious forces fell out of the Euler–Lagrange equation of (4.19): the Coriolis force from the cross term $T_1$ (half of it from $d\mathbf{p}'/dt$, half from $\partial L/\partial\mathbf{r}'$, which is why the factor is 2), the centrifugal force from $T_0$, and the Euler force from differentiating $\boldsymbol\omega$ inside $\mathbf{p}'$. No basis vectors were differentiated. If you set $\dot{\boldsymbol\omega} = 0$ at step 4 you lose the Euler term silently, which is the §8 warning. Example 4.3 is the special case of (4.19) with $\boldsymbol\omega = \omega\hat{\mathbf{z}}$ constant and $\mathbf{r}'$ confined to a hoop: its $T_0 = \frac12 mR^2\omega^2\sin^2\theta$ is exactly $\frac12 m|\boldsymbol\omega\times\mathbf{r}'|^2$, and the effective potential (4.7) is $V - T_0$, the centrifugal potential. The project integrates a projectile from (4.19) with $\boldsymbol\omega = \omega\hat{\mathbf{z}}$ and $V = mgz$, integrates the same projectile in the inertial frame (a parabola), rotates the inertial solution into the rotating frame with $\mathcal{R}(t)^\top$, and asserts the two agree to $10^{-8}$.

### 4.10 When the Lagrangian method is the wrong tool

Honest list. The method is built on three things: a configuration space with coordinates, forces that are either conservative or expressible as $Q_i$, and constraints that are ideal and holonomic. Break one and the method bends or breaks.

- **Impulsive forces and collisions.** A collision is an infinite force for zero time. Lagrange's equations hold in the integrated form $\Delta p_i = \int Q_i\,dt$, but you still need an independent rule (coefficient of restitution, conservation of momentum) for what happens *during* the impulse, and that rule is Newtonian. Use momentum balance across the impact and the Lagrangian between impacts.
- **Friction with an unknown normal force.** Sliding friction is $\mu|\mathbf{N}|$ along $-\mathbf{v}$, and $|\mathbf{N}|$ is a constraint force, which the simple method has thrown away. You can recover $\mathbf{N}$ as a multiplier (§3.8), then feed $\mu|\lambda\nabla f|$ back in as $Q_i^{\text{nc}}$, and solve the coupled system numerically; it works but it is Newton with extra steps. Also, the sign of friction depends on the direction of motion, so the equations change form at every turning point.
- **Non-holonomic constraints.** For a velocity constraint $\sum_i a_{ki}(q)\dot q_i = 0$ that does not integrate, there is no surface and no generalized coordinates on it. The correct procedure is d'Alembert with the virtual displacements restricted by $\sum_i a_{ki}\delta q_i = 0$, which gives (3.18) with $a_{ki}$ in place of $\partial f_k/\partial q_i$ (Exercise C1 does the rolling disk). What does **not** work is substituting the constraint into $L$ and varying: that answers a different question (the "vakonomic" problem) and gives equations of motion that experiments contradict for rolling bodies. If you take one thing from this paragraph: a velocity constraint restricts $\delta q$, not $L$.
- **Dissipation beyond linear drag.** Rayleigh's function handles forces linear in velocity. Quadratic drag, dry friction, hysteresis: compute $Q_i^{\text{nc}}$ from (3.6) directly, and expect no conserved quantity to help you.
- **When Newton is faster.** One particle, Cartesian coordinates, no constraints: $\mathbf{F} = m\mathbf{a}$ is one line and $L$ is three. The method pays for itself exactly when there are constraints or curvilinear coordinates. Knowing when *not* to reach for it is part of knowing it.

## 5. Why it matters and where it connects

- **CM-4** starts from (3.12) and the observation in §3.5 that a cyclic coordinate has a conserved momentum, and proves Noether's theorem: every continuous symmetry of $L$ gives a conserved quantity. The energy function $h$ of (3.15) is its central example, and the distinction $h \ne T + V$ for the rotating hoop (Example 4.3) is a CM-4 exercise in disguise.
- **CM-5** repeats the reduction of Example 4.7 for a planet: conserve $p_\phi$, fold it into $V_{\text{eff}}(r) = p_\phi^2/(2mr^2) + V(r)$, and the two-body problem becomes one-dimensional. The centrifugal barrier, the circular orbit as a minimum of $V_{\text{eff}}$, and the small-oscillation frequency about it (Exercise B6) all reappear with $\sin\theta \to r$.
- **CM-6** takes $M_{ij}$ from (3.13) and $K_{ij} = \partial^2V/\partial q_i\partial q_j$ at an equilibrium and solves $K\mathbf{a} = \Omega^2M\mathbf{a}$ in general, which is Example 4.4's small-oscillation calculation with the spectral theorem (LA-8) doing the work. Positive definiteness of $M$ (§3.6, step 4) is what makes the generalized eigenproblem well posed.
- **CM-7** writes the Lagrangian of a rigid body in Euler angles by the procedure of §3.10 and derives Euler's equations. The rotating-frame Lagrangian (4.19) is the bridge: a rigid body is a rotating frame with mass attached.
- **CM-8** applies a Legendre transform to $L$ in the velocities, turning $h(q, \dot q)$ into $H(q, p)$. The generalized momenta (3.12), including $\mathbf{p} = m\mathbf{v} + q_e\mathbf{A}$ from (4.18), become the coordinates of phase space. The Routhian mentioned in Example 4.7 is a Legendre transform in some variables only.
- **QM-2, QM-4**: canonical quantization replaces the generalized momentum (3.12), not $m\mathbf{v}$, by $-i\hbar\partial/\partial q$. For a charged particle that gives $(-i\hbar\nabla - q_e\mathbf{A})^2/2m$, which is why (4.18) matters. CM-8 → QM-4 is the formal route; the content starts here.
- **EM-5, EM-8**: the gauge change $\mathbf{A}\to\mathbf{A} + \nabla\chi$, $\Phi\to\Phi - \partial_t\chi$ changes (4.17) by the total derivative $q_e\,d\chi/dt$, which does not change the equations of motion. That is the mechanical statement of gauge invariance.
- **LA-5, LA-8, LA-12**: the multiplier derivation (§3.8, step 2) is $(W^\perp)^\perp = W$. The mass matrix is a Gram matrix $J^\top\mathcal{M}J$, positive definite exactly when the tangent vectors are independent. In LA-12's language $M$ is a Riemannian metric on configuration space, and the free motion ($V = 0$) follows its geodesics.
- **ML-2**: the stationarity conditions of $\min f(x)$ subject to $g_k(x) = 0$ are $\nabla f = \sum_k\lambda_k\nabla g_k$: equation (3.18) with the kinetic terms removed (statics). The KKT matrix of equality-constrained quadratic programming is (3.21) with $M$ the Hessian. "The multiplier is the constraint force" and "the dual variable is the shadow price" are the same sentence in two dialects.
- **For the engineering courses**: multibody dynamics software (the kind that simulates linkages and vehicles) integrates (3.21) with constraint stabilization, and vibration analysis is CM-6 applied to the $M$ and $K$ you learned to compute here.

## 6. Exercises

See [exercises.md](exercises.md). **A1–A5** are drills: counting degrees of freedom, classifying constraints, one generalized force, one mass matrix, one dissipation check. **B1–B8** are the module. B1 (the two chain-rule identities, verified on the rotating hoop), B2 (pendulum tension by multipliers, with the DAE code), B4 (the parabolic wire, which bifurcates differently from the hoop), and B7 (Rayleigh damping on the double pendulum) are essential; B3, B5, B6, B8 extend the worked examples. **C1–C3**: the rolling disk's non-holonomic constraint and the correct equations for it, Kapitza's inverted pendulum, and the adiabatic invariant of a pendulum whose length changes slowly.

## 7. Coding project

`project.py` builds a symbolic-to-numeric pipeline and uses it to check every claim in §4.

1. **`lagrangian_to_eom(L, qs, t, Q=None, params=None)`.** Takes a SymPy Lagrangian in terms of `sp.Function` coordinates, forms the Euler–Lagrange expressions, reads off the mass matrix $M$ as the Jacobian of those expressions with respect to the accelerations (equation (3.14)), solves the linear system for $\ddot q$, and returns both the symbolic accelerations and a `rhs(t, y)` function ready for `scipy.integrate.solve_ivp`. The optional `Q` adds generalized forces (3.11). A companion `energy_function` builds $h$ from (3.15).
2. **Simple pendulum**: the pipeline reproduces (4.1) symbolically.
3. **Atwood machine**: the reduced Lagrangian gives (4.2); then `dae_accelerations` solves the KKT system (3.21) at every time step and the recorded $\lambda$ is asserted to equal $-23.544$ N. The constraint drift is measured.
4. **Rotating hoop**: (4.6) checked symbolically; $h$ conserved and $T + V$ not; the period about $\theta^*$ at $\omega = 2\omega_0$ measured from zero crossings and compared with $0.7325$ s; the bifurcation diagram drawn from the analytic branches with integrated trajectories confirming stability on each one.
5. **Double pendulum**: SymPy's accelerations are asserted equal to the hand-expanded (4.10); the small-oscillation eigenproblem is solved with `scipy.linalg.eigh(K, M0)` and gives $(2\mp\sqrt2)g/\ell$ with mode ratios $\pm\sqrt2$; two trajectories $10^{-8}$ apart are integrated, their separation plotted on a log scale, and the Lyapunov exponent fitted; a $3\times4$ grid of snapshots; energy conservation.
6. **Oscillating support, sliding wedge, spherical pendulum, charged particle**: (4.11), (4.13), (4.15)–(4.16), (4.18) checked symbolically or by integration, with the conserved momenta and energies verified.
7. **Rotating frame**: the equations of motion of (4.19) are asserted equal to (4.20) symbolically; a projectile is integrated in the rotating frame, the inertial parabola is rotated into that frame, and the two are overlaid and asserted to agree.
8. **Rayleigh**: $dE/dt = -2\mathcal{F}$ verified on a damped pendulum.

**Build it yourself first:** write `lagrangian_to_eom` before reading it. The key step is extracting $M$: differentiate each Euler–Lagrange expression with respect to each $\ddot q_j$. Then write `dae_accelerations(M, F, J, gamma)` as a single `np.linalg.solve` on the block matrix of (3.21). Compare with the file.

**Your turn** (at the end of the script): a pendulum whose length changes slowly (adiabatic invariant, Exercise C3), a bead on a rotating parabolic wire (Exercise B4, numerically), and Rayleigh damping on the double pendulum (Exercise B7).

## 8. Where you'll struggle

- **Choosing generalized coordinates is a skill, and the only way to get it is about twenty problems.** No rule generates the right coordinates. Angles for things that pivot, distances along constraints for things that slide, center-of-mass plus relative coordinates for pairs. You will choose badly several times (a Cartesian coordinate for a pendulum bob, which produces a Lagrangian with $\sqrt{\ell^2 - x^2}$ in it) and the Lagrangian itself will tell you, by becoming ugly. Do all of §4 with the notes closed, then the B exercises, then Taylor's chapter 7 problems until it stops hurting.
- **D'Alembert's principle is slippery, and "virtual" is the slippery word.** A virtual displacement is at frozen time. For a scleronomic constraint this makes no difference and you can be sloppy. For a rheonomic one (rotating hoop, oscillating pivot) the constraint force *does* do real work, the energy $T + V$ is not conserved, and students conclude the principle has failed. It has not: it says the constraint force does no *virtual* work, and the frozen time is exactly what excludes the real work done by the motor. If you find yourself writing "the hoop does work on the bead so d'Alembert does not apply," reread Definition 3.4 and Example 4.3's last paragraph.
- **The multiplier's sign.** $\mathbf{F}^c = +\lambda\nabla f$. The sign of $\lambda$ depends on how you wrote $f$ ($x_1 + x_2 - \ell$ versus $\ell - x_1 - x_2$). Students compute $\lambda$ correctly and then report the constraint force pointing the wrong way because they never checked which way $\nabla f$ points. Before interpreting $\lambda$, compute $\nabla f$ at a concrete configuration and draw the arrow. Example 4.2 and Exercise B2 both come out with $\lambda < 0$ for this reason.
- **Substituting a conserved momentum into $L$.** For the spherical pendulum, $p_\phi = m\ell^2\sin^2\theta\,\dot\phi$ is constant, and it is tempting to replace $\dot\phi$ in $L$ and vary the result. The centrifugal term then has the wrong sign. The reason: $L$ is a function of $(q, \dot q)$, and the variation that gives Lagrange's equations holds $\dot\phi$ independent of $\theta$; substituting makes $\dot\phi$ depend on $\theta$ and the variation is a different variation. Substitute into the equation of motion or into $h$ (Example 4.7), or learn the Routhian in CM-8.
- **The $\dot q^2$ terms from a $q$-dependent mass matrix.** When $M_{ij}$ depends on $q$, $\frac{d}{dt}\frac{\partial L}{\partial\dot q_i}$ contains $\sum_{jk}\frac{\partial M_{ij}}{\partial q_k}\dot q_k\dot q_j$, and $\frac{\partial L}{\partial q_i}$ contains $\frac12\sum_{jk}\frac{\partial M_{jk}}{\partial q_i}\dot q_j\dot q_k$. These partially cancel (Example 4.4) but not completely, and the survivors ($\dot\theta_2^2\sin\Delta$ in (4.10a)) are the terms students drop. Expand $\frac{d}{dt}$ with the full chain rule every time; never write $\frac{d}{dt}(M\dot q) = M\ddot q$ unless $M$ is constant.
- **The double pendulum's equations are long and you will make a typo.** Everyone does. (4.10) is right because SymPy says so, not because the author is careful. Use `lagrangian_to_eom` to check your hand derivation of every system with more than one coordinate, every time.
- **$T_0$ is not "just a constant."** For the rotating hoop, $T_0 = \frac12mR^2\omega^2\sin^2\theta$ depends on $\theta$ and is the whole bifurcation. Students who learned "constants in $L$ do not matter" sometimes drop it. Only terms that depend on $t$ alone (or are additive constants) can be dropped.
- **The rotating-frame derivation hides $\dot{\boldsymbol\omega}$.** If you assume $\boldsymbol\omega$ constant at step 4 of Example 4.9, the Euler force $-m\dot{\boldsymbol\omega}\times\mathbf{r}'$ never appears and you will not notice its absence. For Earth it is negligible; for a spinning-up centrifuge or a satellite reorienting itself it is not. Carry $\dot{\boldsymbol\omega}$ until the end and set it to zero last.
- **Non-holonomic constraints look holonomic.** "$\dot x = R\dot\phi\cos\theta$" looks like something you can substitute. You cannot, unless it integrates (§3.2), and whether it integrates depends on whether $\theta$ is a variable or a constant. Exercise C1(d) shows the wrong equations you get by substituting. If you cannot write the constraint as $f(q, t) = 0$ with no dots, stop and use (3.18) with the velocity coefficients in place of $\nabla f$.

## 9. Check yourself

1. How many degrees of freedom does a planar double pendulum have? A rigid body? A bead on a rotating hoop? Why does "degrees of freedom" need qualification for a sphere rolling on a plane?
2. State d'Alembert's principle in one sentence, using the word "virtual" correctly. Then say what a virtual displacement is *not*.
3. Derive the identity (3.9) in four lines, naming the two chain-rule facts you use.
4. Why do constraint forces disappear from Lagrange's equations in good generalized coordinates? Answer with a dot product.
5. Write the mass matrix $M_{ij}$ as a Gram matrix. Why is it positive definite? When does it fail to be?
6. Write the multiplier form of Lagrange's equations and say in one sentence what $\lambda_k\nabla f_k$ is. What determines the sign of $\lambda_k$?
7. For the rotating hoop at $\omega = 2\omega_0$, where are the equilibria and which are stable? What is conserved, and why is it not $T + V$?
8. Write the double pendulum's kinetic energy from memory. What is $\det M$ and why is it never zero?
9. The spherical pendulum: what is $p_\phi$, what is $V_{\text{eff}}(\theta)$, and what mistake must you avoid in deriving it?
10. What is the generalized momentum of a charged particle in a magnetic field, and why is it not $m\mathbf{v}$?
11. From $L = \frac12m|\dot{\mathbf{r}}' + \boldsymbol\omega\times\mathbf{r}'|^2 - V$, which term produces the Coriolis force, and why does it come with a factor of 2?
12. Give one system for which the Lagrangian method is the wrong tool, and say what to use instead.

## 10. Pair with

- **Taylor, *Classical Mechanics*, ch. 7** (§7.1–7.7, 7.9–7.10). The clearest undergraduate treatment; his examples overlap ours (Atwood, block on wedge, bead on hoop) with different coordinate choices, which is instructive. Read §7.4 for his proof of the equivalence with Newton, and the end-of-chapter problems for the twenty you need. **Ch. 9** derives the fictitious forces by the vector method; compare with Example 4.9.
- **Goldstein, Poole & Safko, *Classical Mechanics*, ch. 1–2.** §1.3–1.4 is the canonical source for d'Alembert's principle and the derivation in our §3.4, with the same notation. §1.5 does the velocity-dependent potential of the Lorentz force. §2.4 does multipliers and the rolling-disk example of §3.2. Denser than Taylor and more honest about non-holonomic constraints.
- **MIT 8.223 (Classical Mechanics II, IAP), lectures 1–6** on OCW. A fast, problem-driven pass through exactly this material by someone working at a board; good for seeing the procedure of §3.10 executed live on unfamiliar systems.
- **Landau & Lifshitz, *Mechanics*, §1–§5 and the problems after §5.** Nine pages that derive Lagrange's equations from Galilean invariance alone, then five worked problems that are a stiffer version of our §4. Read after, not before, this module. §30 Problem 1 is Kapitza's pendulum (Exercise C2).
- **Strogatz, *Nonlinear Dynamics and Chaos*, §3.4 (pitchfork bifurcations) and §6.5 (the bead on a rotating hoop).** If the bifurcation in Example 4.3 interested you, this is the general theory, and Strogatz uses our exact system as his example.
