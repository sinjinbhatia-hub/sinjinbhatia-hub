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
