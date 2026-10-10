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
