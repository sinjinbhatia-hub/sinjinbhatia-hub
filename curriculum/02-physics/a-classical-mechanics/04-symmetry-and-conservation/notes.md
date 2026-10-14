# CM-4: Symmetry and Conservation: Noether's Theorem, Cyclic Coordinates, the Energy Function

**Time:** ~8 hours · **Prerequisites:** CM-3 · **Unlocks:** CM-5, CM-8; later QM-5 (symmetry and degeneracy) and QM-9 (Noether in quantum mechanics)

> **Conventions in this module.** Generalized coordinates are $q = (q_1, \dots, q_n)$, velocities $\dot q_i = dq_i/dt$, and the Lagrangian is $L(q, \dot q, t) = T - V$ as in CM-3. The generalized momentum conjugate to $q_i$ is $p_i = \partial L/\partial \dot q_i$, and the Euler–Lagrange equations are $\dot p_i = \partial L/\partial q_i$. Particles are labeled $a = 1, \dots, N$ with masses $m_a$ and Cartesian positions $\mathbf{r}_a$ (bold is a 3-vector); $M = \sum_a m_a$ and $\mathbf{R} = \sum_a m_a \mathbf{r}_a / M$ is the center of mass. **Italic $L$ is the Lagrangian; bold $\mathbf{L} = \sum_a \mathbf{r}_a \times \mathbf{p}_a$ is angular momentum.** The energy function is $h(q, \dot q, t) = \sum_i \dot q_i\, \partial L/\partial \dot q_i - L$ (Goldstein's $h$; it becomes the Hamiltonian $H$ in CM-8 once written in terms of $p$). "Conserved" means $dQ/dt = 0$ along every solution of the Euler–Lagrange equations, never identically. $\varepsilon$ is an infinitesimal parameter and every transformation is kept to first order in $\varepsilon$. The Kepler potential is $V = -\kappa/r$ with $\kappa = GMm$ (we need $k$ for the degree of a homogeneous function, so the force constant is $\kappa$, not $k$ as in Goldstein). The unit vector along $\mathbf{r}$ is $\hat{\mathbf{r}} = \mathbf{r}/r$.

## 0. Where this sits

CM-3 gave you the Euler–Lagrange equations in any coordinates and a workshop of worked systems. You have already seen, in passing, that a coordinate missing from $L$ produces a conserved momentum. This module turns that observation into the organizing principle of mechanics: every continuous symmetry of $L$ produces a conserved quantity, and (for the quantities you care about) every conserved quantity comes from a symmetry.

You need CM-3's Euler–Lagrange equations and its rotating-hoop example, CM-2's fact that adding $dF/dt$ to $L$ does not change the equations of motion, and F6's Beltrami identity. What this unlocks: CM-5 uses $p_\theta$ and $E$ to reduce the two-body problem to one radial equation and uses the Laplace–Runge–Lenz vector to get the orbit without integrating; CM-8 turns $h$ into the Hamiltonian and turns "hidden symmetry" into a precise statement with Poisson brackets; QM-5 and QM-6 use the same symmetry $\to$ conservation $\to$ degeneracy chain for angular momentum and hydrogen.

## 1. The problem that forced this idea

Conservation laws were found one at a time, over two and a half centuries, and nobody knew why there were exactly the ones there were.

**Momentum.** Descartes, in the *Principles of Philosophy* (1644), declared that God had put a fixed "quantity of motion" into the world, size times speed, and that it was conserved. It was the right idea with the wrong quantity: speed is not velocity, and Descartes's scalar rule gave wrong answers for collisions. Huygens, Wren and Wallis fixed collisions in 1668–69, and Newton's *Principia* (1687) made the law vectorial: the total $\sum m\mathbf{v}$ of an isolated system is constant. That is Newton's third law in disguise, and Newton knew it.

**Energy.** Leibniz, in 1686, attacked Descartes and proposed that the conserved thing was *vis viva*, $mv^2$. The resulting argument between Cartesians and Leibnizians ran for sixty years; d'Alembert finally called it a dispute about words in 1743, and he was half right: both $m\mathbf{v}$ and $\tfrac12 mv^2$ are conserved in elastic collisions, and the question of what happens to vis viva in *inelastic* collisions, where it visibly disappears, was not settled until the 1840s. Mayer (1842), Joule (1843–50) and Helmholtz (1847, *On the Conservation of Force*) established that what disappears as motion reappears as heat, and the word "energy" in its modern sense dates from the 1850s. So energy conservation was found two centuries after momentum conservation, and by a completely different route.

**Angular momentum.** Kepler's second law (1609), equal areas in equal times, is angular momentum conservation for a planet, but nobody could say so until Newton showed in Proposition 1 of the *Principia* that the area law follows from any central force. Euler, in papers between 1750 and 1775, stated the law of the moment of momentum as a principle in its own right, independent of $\mathbf{F} = m\mathbf{a}$, which is what you need for rigid bodies and fluids.

By 1900, then, mechanics had ten conserved quantities for an isolated system: energy, three components of momentum, three of angular momentum, and the three components of the statement that the center of mass moves uniformly. Ten is also the number of parameters of the Galilean group (one time shift, three translations, three rotations, three boosts). Nobody connected the two numbers.

The connection was made by Emmy Noether in 1918, and for a reason that had nothing to do with mechanics. Hilbert and Klein in Göttingen were trying to understand energy conservation in Einstein's general relativity, where the would-be energy conservation law seemed to be an identity, true for any metric whatsoever, and therefore to say nothing. Hilbert had invited Noether to Göttingen in 1915 for her expertise in invariant theory (the university then refused to let her teach under her own name; Hilbert's reply to the objectors, that the faculty senate was not a bathhouse, is one of the better-documented remarks in the history of mathematics). Her paper *Invariante Variationsprobleme* contains two theorems. Theorem I: if the action is invariant under a continuous group of transformations with finitely many parameters, each parameter gives a conserved quantity, and conversely. Theorem II: if the group has infinitely many parameters (like the coordinate changes of general relativity), you get identities among the equations of motion instead, which is exactly the situation Hilbert had been puzzled by. Theorem II answered the question she had been asked. Theorem I is what this module is about.

Physicists ignored both for about forty years. The version in every textbook today comes from a 1951 *Reviews of Modern Physics* article by E. L. Hill, which simplified Theorem I to the case of point transformations and dropped the converse (Kosmann-Schwarzbach's *The Noether Theorems* documents the reception in detail). When gauge theories arrived in the 1960s and 70s, Noether's theorems became the way physics is organized: you write down the symmetries, and the conservation laws, the degeneracies, and often the interactions follow. This module is the finite-dimensional case, which is all mechanics needs, and you will see that the whole proof is three lines once the right question is asked.

## 2. Intuition

### A symmetry is a flat direction; the conserved thing is the momentum along it

Take a free particle on a line, $L = \tfrac12 m\dot x^2$. Slide the whole experiment to the right by $\varepsilon$: $x \to x + \varepsilon$. Nothing in $L$ changes, because $L$ does not contain $x$ at all. Euler–Lagrange says $\dot p_x = \partial L/\partial x = 0$. So $p_x = m\dot x$ is constant.

That is the whole theorem in its simplest case, and it is worth staring at. Think of $L$ as a landscape over configuration space. The statement "sliding the experiment left does not change $L$" means the landscape is **flat along the $x$ direction**. Euler–Lagrange says the rate of change of momentum is the slope of the landscape. Flat landscape, no slope, constant momentum.

```
   L(x, xdot)                     L(theta, thetadot) for a central force
   flat along x                   flat along theta

        |  .  .  .  .  .                 |  .  .  .  .  .
   xdot |  .  .  .  .  .     theta-dot   |  .  .  .  .  .
        |  .  .  .  .  .                 |  .  .  .  .  .
        +---------------> x              +---------------> theta
        slide: L same                    rotate: L same
        => p_x conserved                 => p_theta conserved
```

Now a particle in a central potential, $L = \tfrac12 m(\dot r^2 + r^2\dot\theta^2) - V(r)$. Rotate the whole experiment by $\varepsilon$: $\theta \to \theta + \varepsilon$. Again $L$ does not change, because $\theta$ does not appear. Flat along $\theta$, so $p_\theta = mr^2\dot\theta$ is constant. That is angular momentum.

### Two coins

Two particles on a line with $L = \tfrac12 m_1\dot x_1^2 + \tfrac12 m_2\dot x_2^2 - V(x_1 - x_2)$. Neither $x_1$ nor $x_2$ is missing from $L$, so neither $p_1$ nor $p_2$ is conserved. But slide *both* particles by the same $\varepsilon$: $x_1 \to x_1 + \varepsilon$, $x_2 \to x_2 + \varepsilon$. The difference $x_1 - x_2$ does not change, so $V$ does not change, and neither does $T$. The landscape is flat along the **diagonal direction** $K = (1, 1)$ in the $(x_1, x_2)$ plane, even though it is not flat along either axis.

What should be conserved? If the picture is right, the momentum "in the direction $K$," which is the component $p_1 K_1 + p_2 K_2 = p_1 + p_2$. Total momentum. (Check with Euler–Lagrange: $\dot p_1 = -V'(x_1 - x_2)$, $\dot p_2 = +V'(x_1 - x_2)$, sum zero.)

### The guess

So here is the shape of the theorem you should now be able to state. A **symmetry** is a direction $K(q)$ in configuration space (one $K_i$ per coordinate, possibly depending on where you are) such that moving every point a little along $K$ leaves $L$ unchanged. The **conserved quantity** is the momentum along that direction, $Q = \sum_i p_i K_i$. Rotation is the case $K = \hat{\mathbf{n}} \times \mathbf{r}$ (the velocity field of a rigid rotation about $\hat{\mathbf{n}}$), and the charge is $\mathbf{p} \cdot (\hat{\mathbf{n}} \times \mathbf{r}) = \hat{\mathbf{n}} \cdot (\mathbf{r} \times \mathbf{p})$: the component of angular momentum along the axis.

Two things the picture does not yet cover. First, time: "waiting a while and repeating the experiment gives the same result" is a symmetry, and it should give energy, but $t$ is not a coordinate and there is no $p_t$. Section 3.4 finds the right object, the energy function $h$, and it is *not* always $T + V$. Second, some transformations change $L$ but only by something harmless. Boosting to a moving frame changes a free particle's $L = \tfrac12 m\dot x^2$ to $\tfrac12 m(\dot x + \varepsilon v)^2 = L + \varepsilon m v\dot x + O(\varepsilon^2)$, and $mv\dot x = \tfrac{d}{dt}(mvx)$ is a total time derivative, which CM-2 showed does not affect the equations of motion. Such "symmetries up to a total derivative" also give conservation laws, with a correction term. The Galilean boost is the case that needs it.

## 3. Formal development

Throughout, $L(q, \dot q, t)$ with $n$ coordinates, $p_i = \partial L/\partial \dot q_i$, and the Euler–Lagrange equations

$$
\frac{d}{dt}\frac{\partial L}{\partial \dot q_i} = \frac{\partial L}{\partial q_i}, \qquad \text{i.e.}\quad \dot p_i = \frac{\partial L}{\partial q_i}. \tag{3.1}
$$

### 3.1 Cyclic coordinates

**Theorem 3.1 (Cyclic coordinates).** If $\partial L/\partial q_j = 0$ (the coordinate $q_j$ is *cyclic* or *ignorable*), then $p_j = \partial L/\partial \dot q_j$ is conserved.

*Proof.* By (3.1), $\dot p_j = \partial L/\partial q_j = 0$. ∎

$$
\frac{\partial L}{\partial q_j} = 0 \;\Longrightarrow\; p_j = \text{const}. \tag{3.2}
$$

Three examples, each one line:

- **Central force in polar coordinates.** $L = \tfrac12 m(\dot r^2 + r^2\dot\theta^2) - V(r)$. $\theta$ is cyclic, so $p_\theta = mr^2\dot\theta$ is conserved. In CM-5 this is $\ell$, and it is why the two-body problem reduces to one dimension.
- **Particle on a cylinder.** Radius $R$, coordinates $(\phi, z)$, $L = \tfrac12 m(R^2\dot\phi^2 + \dot z^2) - V(z)$. $\phi$ is cyclic, so $p_\phi = mR^2\dot\phi$, the angular momentum about the axis, is conserved. The particle spirals with constant angular speed whatever $V(z)$ does to it vertically.
- **Free particle.** $L = \tfrac12 m|\dot{\mathbf{r}}|^2$: all three Cartesian coordinates are cyclic, so $\mathbf{p} = m\dot{\mathbf{r}}$ is conserved.

Now the catch. Write the central-force problem in Cartesian coordinates: $L = \tfrac12 m(\dot x^2 + \dot y^2) - V(\sqrt{x^2 + y^2})$. *No* coordinate is cyclic. The rotational symmetry has not gone anywhere, and $m(x\dot y - y\dot x)$ is still conserved, but Theorem 3.1 cannot see it. Whether a coordinate is cyclic depends on the coordinates; whether $L$ has a symmetry does not. We need a statement that does not care which coordinates you chose.

### 3.2 Noether's theorem for point symmetries

**Definition 3.2 (Infinitesimal point transformation; invariance).** An infinitesimal point transformation is a family

$$
q_i \;\to\; q_i + \varepsilon K_i(q, t), \qquad i = 1, \dots, n,
$$

with $K_i$ given functions (the *generator*). It induces on velocities $\dot q_i \to \dot q_i + \varepsilon \dot K_i$, where $\dot K_i = \sum_j \frac{\partial K_i}{\partial q_j}\dot q_j + \frac{\partial K_i}{\partial t}$ is the total time derivative of $K_i$ along the path. The transformation is a **symmetry** of $L$ if $L$ is unchanged to first order in $\varepsilon$, for *every* $(q, \dot q, t)$, not only along solutions:

$$
\delta L \;\equiv\; \frac{d}{d\varepsilon}\Big|_{\varepsilon = 0} L\big(q + \varepsilon K,\ \dot q + \varepsilon\dot K,\ t\big) \;=\; \sum_{i=1}^n \left( \frac{\partial L}{\partial q_i} K_i + \frac{\partial L}{\partial \dot q_i}\dot K_i \right) \;=\; 0. \tag{3.3}
$$

The middle expression is the chain rule; (3.3) is the condition you check in practice.

**Theorem 3.3 (Noether's theorem, point symmetries).** If (3.3) holds, then

$$
Q \;=\; \sum_{i=1}^n p_i K_i \;=\; \sum_{i=1}^n \frac{\partial L}{\partial \dot q_i} K_i \tag{3.4}
$$

is conserved along every solution of the Euler–Lagrange equations.

▶ **Try it first.** Differentiate (3.4) with respect to time, use the product rule, then replace $\dot p_i$ using (3.1). Compare with (3.3).

**Derivation.**

1. Product rule: $\dfrac{dQ}{dt} = \sum_i \big(\dot p_i K_i + p_i \dot K_i\big)$.
2. On a solution, (3.1) gives $\dot p_i = \partial L/\partial q_i$, and $p_i = \partial L/\partial\dot q_i$ by definition. So $\dfrac{dQ}{dt} = \sum_i \left(\dfrac{\partial L}{\partial q_i}K_i + \dfrac{\partial L}{\partial \dot q_i}\dot K_i\right)$.
3. That is $\delta L$, which is zero by (3.3). ∎

■ **What just happened.** Step 1 is algebra. Step 2 is the *only* place the equations of motion enter, which is why $Q$ is conserved on solutions and not identically. Step 3 is the hypothesis, which must hold for all $(q, \dot q, t)$: the symmetry is a property of the landscape, the conservation a property of the paths. Theorem 3.1 is the case $K = (0, \dots, 1, \dots, 0)$, a unit step in the $q_j$ direction, for which (3.3) reads $\partial L/\partial q_j = 0$ and (3.4) reads $Q = p_j$.

Now the harmless failure. Suppose $L$ is not invariant, but changes by a total time derivative:

**Theorem 3.4 (Noether's theorem, quasi-symmetries).** If there is a function $F(q, t)$ such that

$$
\delta L = \sum_i \left( \frac{\partial L}{\partial q_i} K_i + \frac{\partial L}{\partial \dot q_i}\dot K_i \right) = \frac{dF}{dt} = \sum_j \frac{\partial F}{\partial q_j}\dot q_j + \frac{\partial F}{\partial t},
$$

then $Q = \sum_i p_i K_i - F$ is conserved.

$$
Q \;=\; \sum_i p_i K_i \;-\; F. \tag{3.5}
$$

*Proof.* As before, $\dfrac{d}{dt}\sum_i p_i K_i = \delta L$ on solutions, so $\dfrac{dQ}{dt} = \delta L - \dfrac{dF}{dt} = 0$. ∎

Why this is the right generalization: under the transformation the action $S = \int L\,dt$ changes by $\varepsilon\int (dF/dt)\,dt = \varepsilon\,[F]_{t_1}^{t_2}$, a boundary term, and CM-2 showed that Lagrangians differing by $dF/dt$ give the same Euler–Lagrange equations. So a quasi-symmetry of $L$ is a genuine symmetry of the *dynamics*. The Galilean boost (§3.3c) is the case where you cannot avoid this. Note the hypothesis: $F$ depends on $(q, t)$ only, so $dF/dt$ is linear in $\dot q$. Not every $\delta L$ of that form is a total derivative; the project's `total_derivative_potential` constructs $F$ when it exists and refuses when it does not.

**Theorem 3.5 (Converse, for charges linear in momentum).** Suppose $L$ is regular (the matrix $W_{ij} = \partial^2 L/\partial\dot q_i\partial\dot q_j$ is invertible, so that (3.1) can be solved for the accelerations and every point $(q, \dot q, t)$ lies on exactly one solution). If $Q = \sum_i p_i K_i(q, t) - F(q, t)$ is conserved, then $q_i \to q_i + \varepsilon K_i$ is a quasi-symmetry with $\delta L = dF/dt$.

*Proof.* On a solution, $\dot Q = \delta L - dF/dt$ by the computation in Theorem 3.4, and this is zero. The quantity $\delta L - dF/dt$ depends only on $(q, \dot q, t)$ (no accelerations, since (3.1) was used to remove $\dot p_i$). By regularity every $(q, \dot q, t)$ is a point of some solution, so $\delta L - dF/dt = 0$ everywhere. ∎

■ **What just happened.** For conserved quantities that are *linear in the momenta*, symmetry and conservation are the same thing. This is the part of Noether's Theorem I that Hill's 1951 simplification dropped. Conserved quantities that are *not* linear in $p$ (energy is quadratic; the Laplace–Runge–Lenz vector of §3.8 is quadratic) still come from symmetries, but the transformation then has to move velocities as well as positions. You will meet the first such case in §3.4 and the precise statement in CM-8.

### 3.3 The three Galilean symmetries, derived

**(a) Translation $\Rightarrow$ total momentum.** $N$ particles, $L = \sum_a \tfrac12 m_a|\dot{\mathbf{r}}_a|^2 - V$, with $V$ depending only on the differences $\mathbf{r}_a - \mathbf{r}_b$ (no external field). Translate everything along a unit vector $\hat{\mathbf{n}}$: $K_a = \hat{\mathbf{n}}$ for every $a$, so $\dot K_a = 0$.

1. $\partial L/\partial \mathbf{r}_a = -\nabla_a V$, so $\delta L = -\hat{\mathbf{n}}\cdot\sum_a \nabla_a V$.
2. $V$ depends on differences only, so translating every argument by the same vector leaves it unchanged: $V(\mathbf{r}_1 + \mathbf{s}, \dots, \mathbf{r}_N + \mathbf{s}) = V(\mathbf{r}_1, \dots, \mathbf{r}_N)$ for all $\mathbf{s}$. Differentiating with respect to $\mathbf{s}$ at $\mathbf{s} = 0$: $\sum_a \nabla_a V = \mathbf{0}$. (For two bodies with $V(\mathbf{r}_1 - \mathbf{r}_2)$: $\nabla_1 V = -\nabla_2 V$. This is Newton's third law.)
3. So $\delta L = 0$ and $Q = \sum_a \mathbf{p}_a\cdot\hat{\mathbf{n}} = \hat{\mathbf{n}}\cdot\mathbf{P}$ is conserved. Three independent $\hat{\mathbf{n}}$: the total momentum $\mathbf{P} = \sum_a \mathbf{p}_a$ is conserved. ∎

**(b) Rotation $\Rightarrow$ angular momentum about the axis.** An infinitesimal rotation by angle $\varepsilon$ about the unit vector $\hat{\mathbf{n}}$ moves every point by

$$
\mathbf{r}_a \to \mathbf{r}_a + \varepsilon\,\hat{\mathbf{n}}\times\mathbf{r}_a, \qquad K_a = \hat{\mathbf{n}}\times\mathbf{r}_a, \qquad \dot K_a = \hat{\mathbf{n}}\times\dot{\mathbf{r}}_a. \tag{3.6}
$$

(Why this is a rotation: $\hat{\mathbf{n}}\times\mathbf{r}$ is perpendicular to both the axis and $\mathbf{r}$, has length $r\sin(\text{angle to axis})$, which is the distance from the axis, and $\varepsilon$ times that is the arc length for a rotation by $\varepsilon$. It is the velocity field $\boldsymbol{\omega}\times\mathbf{r}$ of CM-7 with $\boldsymbol{\omega} = \hat{\mathbf{n}}$.)

▶ **Try it first.** Show that every dot product $\mathbf{u}\cdot\mathbf{v}$ between vectors that both transform by (3.6) is unchanged to first order. Then conclude $\delta L = 0$ for any $L$ built from dot products of the $\mathbf{r}_a$ and $\dot{\mathbf{r}}_a$ (which includes every central potential and every $V$ depending on inter-particle distances).

**Derivation.**

1. $(\mathbf{u} + \varepsilon\hat{\mathbf{n}}\times\mathbf{u})\cdot(\mathbf{v} + \varepsilon\hat{\mathbf{n}}\times\mathbf{v}) = \mathbf{u}\cdot\mathbf{v} + \varepsilon\big[\mathbf{u}\cdot(\hat{\mathbf{n}}\times\mathbf{v}) + \mathbf{v}\cdot(\hat{\mathbf{n}}\times\mathbf{u})\big] + O(\varepsilon^2)$.
2. The scalar triple product is unchanged by cyclic permutation and flips sign under a swap: $\mathbf{u}\cdot(\hat{\mathbf{n}}\times\mathbf{v}) = \hat{\mathbf{n}}\cdot(\mathbf{v}\times\mathbf{u})$ and $\mathbf{v}\cdot(\hat{\mathbf{n}}\times\mathbf{u}) = \hat{\mathbf{n}}\cdot(\mathbf{u}\times\mathbf{v}) = -\hat{\mathbf{n}}\cdot(\mathbf{v}\times\mathbf{u})$. The bracket vanishes.
3. $T = \sum_a \tfrac12 m_a\,\dot{\mathbf{r}}_a\cdot\dot{\mathbf{r}}_a$ and any $V$ depending on $|\mathbf{r}_a|$ and $|\mathbf{r}_a - \mathbf{r}_b|$ are built from such dot products, so $\delta L = 0$.
4. The charge: $Q = \sum_a \mathbf{p}_a\cdot(\hat{\mathbf{n}}\times\mathbf{r}_a) = \sum_a \hat{\mathbf{n}}\cdot(\mathbf{r}_a\times\mathbf{p}_a) = \hat{\mathbf{n}}\cdot\mathbf{L}$, again by cyclic permutation. ∎

■ **What just happened.** The conserved quantity attached to rotation about $\hat{\mathbf{n}}$ is the component of angular momentum along $\hat{\mathbf{n}}$. If $L$ is invariant under rotations about every axis, all of $\mathbf{L}$ is conserved; if only about $\hat{\mathbf{z}}$ (a particle in a vertical uniform field plus a central force, say), only $L_z$. Example 4.1 does step 4 in components, in case the triple-product manipulations felt too slick.

**(c) Galilean boost $\Rightarrow$ the center of mass moves uniformly.** Change to a frame moving at velocity $\varepsilon\mathbf{v}$: $\mathbf{r}_a \to \mathbf{r}_a + \varepsilon\mathbf{v}t$. So $K_a = \mathbf{v}t$ and $\dot K_a = \mathbf{v}$. This is the case where $K$ depends on $t$ and where $L$ is *not* invariant.

1. $\delta L = \sum_a\big[(-\nabla_a V)\cdot\mathbf{v}t + m_a\dot{\mathbf{r}}_a\cdot\mathbf{v}\big]$. The first sum vanishes as in (a). So $\delta L = \mathbf{v}\cdot\sum_a m_a\dot{\mathbf{r}}_a \ne 0$.
2. But $\sum_a m_a\dot{\mathbf{r}}_a = \dfrac{d}{dt}\sum_a m_a\mathbf{r}_a = \dfrac{d}{dt}(M\mathbf{R})$. So $\delta L = \dfrac{dF}{dt}$ with $F = M\mathbf{v}\cdot\mathbf{R}$, a function of positions only. Theorem 3.4 applies.
3. $Q = \sum_a\mathbf{p}_a\cdot\mathbf{v}t - M\mathbf{v}\cdot\mathbf{R} = \mathbf{v}\cdot\big(\mathbf{P}t - M\mathbf{R}\big)$. This holds for every $\mathbf{v}$, so the vector

$$
\mathbf{G} = \mathbf{P}\,t - M\mathbf{R} \tag{3.7}
$$

is conserved. ∎

■ **What just happened.** $\mathbf{G}$ constant means $M\mathbf{R}(t) = M\mathbf{R}(0) + \mathbf{P}t$: the center of mass moves in a straight line at velocity $\mathbf{P}/M$. Be honest with yourself about what is new here. Given that $\mathbf{P}$ is conserved and that $\mathbf{R}$ is *defined* as the mass-weighted centroid, this follows in one line. What Noether adds is *which* point moves uniformly: the one weighted by the masses that appear in the kinetic term. In relativistic mechanics the same boost symmetry picks out the energy-weighted centroid instead, and that is not obvious at all. Example 4.2 works this out for two bodies with all the terms visible.

### 3.4 Time translation and the energy function

Waiting is a symmetry when $L$ does not depend explicitly on $t$. But $t$ is not a coordinate, and the recipe (3.4) has no $p_t$. The right object is

$$
h(q, \dot q, t) \;=\; \sum_{i=1}^n \dot q_i\,\frac{\partial L}{\partial \dot q_i} \;-\; L. \tag{3.8}
$$

**Theorem 3.6 (Energy function).** Along every solution of (3.1),

$$
\frac{dh}{dt} = -\frac{\partial L}{\partial t}. \tag{3.9}
$$

In particular, if $L$ has no explicit time dependence, $h$ is conserved.

▶ **Try it first.** Differentiate (3.8) by the product rule; expand $dL/dt$ by the chain rule (remember the explicit $t$); subtract.

**Derivation.**

1. $\dfrac{d}{dt}\sum_i \dot q_i p_i = \sum_i\big(\ddot q_i p_i + \dot q_i\dot p_i\big)$.
2. $\dfrac{dL}{dt} = \sum_i\left(\dfrac{\partial L}{\partial q_i}\dot q_i + \dfrac{\partial L}{\partial \dot q_i}\ddot q_i\right) + \dfrac{\partial L}{\partial t} = \sum_i\left(\dfrac{\partial L}{\partial q_i}\dot q_i + p_i\ddot q_i\right) + \dfrac{\partial L}{\partial t}$.
3. Subtract: the $p_i\ddot q_i$ terms cancel, leaving $\dfrac{dh}{dt} = \sum_i \dot q_i\left(\dot p_i - \dfrac{\partial L}{\partial q_i}\right) - \dfrac{\partial L}{\partial t}$. The bracket is zero by (3.1). ∎

■ **What just happened.** For one coordinate and no explicit $t$, (3.9) says $\dot q\,\partial L/\partial\dot q - L = \text{const}$, which is exactly the **Beltrami identity** of F6 (the first integral of the Euler–Lagrange equation when the integrand has no explicit dependence on the independent variable). Same identity, same proof; there it saved you an integration in the brachistochrone, here it is energy conservation.

Time translation *is* a Noether symmetry, with a twist. Shift the path in time: $q_i(t) \to q_i(t + \varepsilon) = q_i + \varepsilon\dot q_i$, so the generator is $K_i = \dot q_i$, which depends on velocity. The proof of Theorem 3.3 never used that $K$ was velocity-independent (check: steps 1–3 go through with $\dot K_i = \ddot q_i$), and (3.3) gives $\delta L = \sum_i(\partial L/\partial q_i\,\dot q_i + \partial L/\partial\dot q_i\,\ddot q_i) = dL/dt - \partial L/\partial t$. When $\partial L/\partial t = 0$ this is $dL/dt$, a total derivative with $F = L$, and (3.5) gives $Q = \sum_i p_i\dot q_i - L = h$. The project checks this symbolically. This is your first velocity-dependent generator, and the reason $h$ is quadratic in the velocities rather than linear.

**When is $h$ the energy?** Suppose the Cartesian positions are functions of the generalized coordinates only, $\mathbf{r}_a = \mathbf{r}_a(q)$, with no explicit $t$ (the constraints do not move). Then

$$
\dot{\mathbf{r}}_a = \sum_j \frac{\partial\mathbf{r}_a}{\partial q_j}\dot q_j, \qquad T = \sum_a \tfrac12 m_a|\dot{\mathbf{r}}_a|^2 = \tfrac12\sum_{j,k} M_{jk}(q)\,\dot q_j\dot q_k, \quad M_{jk} = \sum_a m_a\frac{\partial\mathbf{r}_a}{\partial q_j}\cdot\frac{\partial\mathbf{r}_a}{\partial q_k},
$$

a **homogeneous quadratic** in the velocities: $T(q, \lambda\dot q) = \lambda^2 T(q, \dot q)$.

**Euler's theorem on homogeneous functions.** If $f(\lambda x_1, \dots, \lambda x_n) = \lambda^k f(x_1, \dots, x_n)$ for all $\lambda > 0$, then

$$
\sum_i x_i\frac{\partial f}{\partial x_i} = k\,f. \tag{3.10}
$$

*Proof.* Differentiate the defining identity with respect to $\lambda$: the chain rule on the left gives $\sum_i x_i\,\partial_i f(\lambda x)$, the right gives $k\lambda^{k-1}f(x)$. Set $\lambda = 1$. ∎

**Theorem 3.7 ($h = E$ for fixed constraints).** If $\mathbf{r}_a = \mathbf{r}_a(q)$ and $V = V(q, t)$ is velocity-independent, then $h = T + V$.

*Proof.* $T$ is homogeneous of degree 2 in $\dot q$, so by (3.10) $\sum_i\dot q_i\,\partial T/\partial\dot q_i = 2T$. $V$ has no $\dot q$. So $h = 2T - (T - V) = T + V$. ∎

If instead $\mathbf{r}_a = \mathbf{r}_a(q, t)$ (moving constraints: a rotating hoop, a driven support), then $\dot{\mathbf{r}}_a = \sum_j(\partial\mathbf{r}_a/\partial q_j)\dot q_j + \partial\mathbf{r}_a/\partial t$ and $T = T_2 + T_1 + T_0$ with $T_2$ quadratic, $T_1$ linear, $T_0$ independent of $\dot q$. Euler's theorem on each piece gives $\sum\dot q_i\partial T/\partial\dot q_i = 2T_2 + T_1$, so

$$
h = 2T_2 + T_1 - (T_2 + T_1 + T_0 - V) = T_2 - T_0 + V, \qquad E - h = T_1 + 2T_0. \tag{3.11}
$$

**The counterexample you must know: the bead on the rotating hoop** (CM-3). A hoop of radius $R$ rotates about its vertical diameter at fixed $\omega$; $\theta$ is the angle from the bottom. Then $T = \tfrac12 mR^2\dot\theta^2 + \tfrac12 mR^2\omega^2\sin^2\theta$ (so $T_1 = 0$, $T_0 = \tfrac12 mR^2\omega^2\sin^2\theta$) and $V = -mgR\cos\theta$. $L$ has no explicit $t$, so $h$ is conserved, but

$$
h = \tfrac12 mR^2\dot\theta^2 - \tfrac12 mR^2\omega^2\sin^2\theta - mgR\cos\theta, \qquad E = T + V = \tfrac12 mR^2\dot\theta^2 + \tfrac12 mR^2\omega^2\sin^2\theta - mgR\cos\theta,
$$

and $E - h = mR^2\omega^2\sin^2\theta = 2T_0$ changes as the bead moves. Interpretation: $h$ is the energy *as seen in the rotating frame*, kinetic energy $\tfrac12 mR^2\dot\theta^2$ plus the effective potential $V_{\text{eff}} = -mgR\cos\theta - \tfrac12 mR^2\omega^2\sin^2\theta$, whose second term is the centrifugal potential. That is conserved because in the rotating frame nothing depends on time. $E$ is the lab-frame energy, and it is *not* conserved because the motor holding $\omega$ fixed does work: $dE/dt = 2mR^2\omega^2\sin\theta\cos\theta\,\dot\theta$, which is the motor's torque $dL_z/dt$ times $\omega$. Example 4.3 puts numbers on this and the project plots both.

### 3.5 Symmetry of $L$ versus symmetry of the equations of motion

If $L$ is invariant (or quasi-invariant) under a transformation, the transformed path of any solution is a solution: the action is the same (up to boundary terms), so stationarity is preserved. The converse fails. Take the free particle $L = \tfrac12 m\dot x^2$ and the **scaling** $x \to \lambda x$. The equation of motion $\ddot x = 0$ is invariant. But $L \to \lambda^2 L$: multiplied by a constant, not shifted by a total derivative. No Noether charge, and indeed nothing new is conserved. Scaling is a symmetry of the physics that is not a symmetry of $L$, and Noether's theorem says nothing about it.

It does say something else, though, and that is the next section.

### 3.6 Mechanical similarity

**Theorem 3.8 (Mechanical similarity; Landau §10).** Let $V$ be homogeneous of degree $k$ in the coordinates: $V(\alpha\mathbf{r}_1, \dots, \alpha\mathbf{r}_N) = \alpha^k V(\mathbf{r}_1, \dots, \mathbf{r}_N)$. If $\mathbf{r}_a(t)$ solves the equations of motion, so does $\alpha\,\mathbf{r}_a(t/\beta)$ with

$$
\beta = \alpha^{1 - k/2}. \tag{3.12}
$$

Consequently, for geometrically similar motions, times scale as $t'/t = (l'/l)^{1 - k/2}$, velocities as $(l'/l)^{k/2}$, energies as $(l'/l)^k$ and angular momenta as $(l'/l)^{1 + k/2}$.

▶ **Try it first.** Substitute $\mathbf{r}_a = \alpha\tilde{\mathbf{r}}_a$, $t = \beta\tilde t$ into $L$ and ask what $\beta$ makes $T$ and $V$ pick up the *same* factor.

**Derivation.**

1. Under $\mathbf{r}_a \to \alpha\mathbf{r}_a$, $t \to \beta t$: $d\mathbf{r}_a/dt \to (\alpha/\beta)\,d\mathbf{r}_a/dt$, so $T \to (\alpha^2/\beta^2)\,T$ and $V \to \alpha^k V$.
2. Choose $\beta$ so that $\alpha^2/\beta^2 = \alpha^k$, i.e. $\beta^2 = \alpha^{2-k}$, $\beta = \alpha^{1-k/2}$. Then $L \to \alpha^k L$.
3. Multiplying $L$ by a constant does not change the Euler–Lagrange equations. So if $\tilde{\mathbf{r}}_a(\tilde t)$ is a solution of the scaled problem, it is a solution of the original one; unwinding the substitution, $\alpha\tilde{\mathbf{r}}_a(t/\beta)$ is a solution whenever $\tilde{\mathbf{r}}_a(t)$ is.
4. The scaled motion traverses a path $\alpha$ times larger in a time $\beta$ times longer: $t'/t = \alpha^{1-k/2}$. Velocities scale as $\alpha/\beta = \alpha^{k/2}$, energies as $\alpha^k$ (both $T$ and $V$ do), angular momenta $mrv$ as $\alpha\cdot\alpha^{k/2}$. ∎

■ **What just happened.** You got three laws for free. **$k = -1$ (gravity, Coulomb):** $t \propto l^{3/2}$, so for similar orbits $T^2 \propto a^3$, Kepler's third law, without solving anything (Example 4.4). **$k = 2$ (harmonic):** $t \propto l^0$, the period does not depend on amplitude; this is the whole reason pendulum clocks work and the reason small oscillations (CM-6) are tractable. **$k = 1$ (uniform gravity):** $t \propto \sqrt{l}$, Galileo's law for falling. For $k = 4$ the period *shrinks* as $1/A$ with amplitude. For $k = -2$ the angular momentum scales as $\alpha^0$: the scaling leaves $\ell$ unchanged, which is why $V \propto -1/r^2$ competes exactly with the centrifugal barrier $\ell^2/2mr^2$ and has no stable orbits (exercise C3). The project measures all these exponents.

### 3.7 The virial theorem

**Theorem 3.9 (Virial theorem).** For a system of particles whose motion is bounded (all $\mathbf{r}_a$ and $\mathbf{p}_a$ stay finite), with time averages $\langle\cdot\rangle$ over long times,

$$
2\langle T\rangle = -\Big\langle\sum_a \mathbf{F}_a\cdot\mathbf{r}_a\Big\rangle, \qquad\text{and if } V \text{ is homogeneous of degree } k:\qquad 2\langle T\rangle = k\langle V\rangle. \tag{3.13}
$$

For periodic motion the average over one period gives (3.13) exactly.

▶ **Try it first.** Let $G = \sum_a\mathbf{p}_a\cdot\mathbf{r}_a$. Compute $dG/dt$ and time-average it.

**Derivation.**

1. $\dfrac{dG}{dt} = \sum_a\big(\dot{\mathbf{p}}_a\cdot\mathbf{r}_a + \mathbf{p}_a\cdot\dot{\mathbf{r}}_a\big) = \sum_a\mathbf{F}_a\cdot\mathbf{r}_a + \sum_a m_a|\dot{\mathbf{r}}_a|^2 = \sum_a\mathbf{F}_a\cdot\mathbf{r}_a + 2T$.
2. Average over $[0, \tau]$: $\left\langle\dfrac{dG}{dt}\right\rangle = \dfrac{G(\tau) - G(0)}{\tau}$. If the motion is bounded, $G$ is bounded, and this tends to $0$ as $\tau\to\infty$. If the motion is periodic with period $\tau$, $G(\tau) = G(0)$ and it is exactly $0$.
3. So $2\langle T\rangle = -\langle\sum_a\mathbf{F}_a\cdot\mathbf{r}_a\rangle$. With $\mathbf{F}_a = -\nabla_a V$ and $V$ homogeneous of degree $k$ in the full set of coordinates, Euler's theorem (3.10) gives $\sum_a\mathbf{r}_a\cdot\nabla_a V = kV$, hence $2\langle T\rangle = k\langle V\rangle$. ∎

■ **What just happened.** For gravity, $k = -1$: $2\langle T\rangle = -\langle V\rangle$ and $E = \langle T\rangle + \langle V\rangle = -\langle T\rangle$. A bound gravitational system has negative total energy, and the more energy it *loses*, the *faster* its parts move: stars heat up as they radiate. For the harmonic oscillator, $k = 2$: $\langle T\rangle = \langle V\rangle$. Example 4.5 does both, including an order-of-magnitude temperature for the Sun's interior from nothing but $M$ and $R$. The name: Clausius (1870) called $-\tfrac12\sum\mathbf{F}\cdot\mathbf{r}$ the *virial*, from *vis*, force.

### 3.8 Hidden symmetry: the Laplace–Runge–Lenz vector

For the Kepler problem, $L = \tfrac12 m|\dot{\mathbf{r}}|^2 + \kappa/r$, the symmetries so far give $E$ and $\mathbf{L}$. There is one more conserved vector:

$$
\mathbf{A} = \mathbf{p}\times\mathbf{L} - m\kappa\,\hat{\mathbf{r}}. \tag{3.14}
$$

**Theorem 3.10.** For $\dot{\mathbf{p}} = -\kappa\hat{\mathbf{r}}/r^2$, $d\mathbf{A}/dt = 0$. Moreover $\mathbf{A}\cdot\mathbf{L} = 0$, $\mathbf{A}$ points from the force center to the pericenter, and $|\mathbf{A}| = m\kappa e$ with $e$ the eccentricity.

Example 4.6 proves conservation by direct differentiation (it is a BAC–CAB computation) and extracts the orbit equation from $\mathbf{A}\cdot\mathbf{r}$ in two lines. The history is a tangle: Jakob Hermann and Johann Bernoulli had it in 1710, Laplace in 1799, Hamilton rediscovered it in 1845 via the hodograph, Runge put it in a 1919 vector-analysis textbook, Lenz used it in 1924 for the old quantum theory of hydrogen, and Pauli used it in 1926 to derive the hydrogen spectrum from matrix mechanics before Schrödinger's equation existed.

Why "hidden"? $\mathbf{A}$ is quadratic in $\mathbf{p}$, so by the remark after Theorem 3.5 no point transformation $\mathbf{r}\to\mathbf{r} + \varepsilon\mathbf{K}(\mathbf{r}, t)$ generates it. The symmetry that does generate it moves positions *and* momenta together: it is a symmetry of phase space (the $(\mathbf{r}, \mathbf{p})$ space of CM-8), not of configuration space. There is nothing you can do to the experiment, no rotation or translation of the apparatus, that corresponds to it; it exists because the $1/r$ potential is special. Count what it buys: $E$, $\mathbf{L}$, $\mathbf{A}$ are seven conserved quantities with two relations ($\mathbf{A}\cdot\mathbf{L} = 0$ and $A^2 = m^2\kappa^2 + 2mEL^2$), so five independent constants in a six-dimensional phase space. The trajectory is pinned to a curve: every bound Kepler orbit *closes*. Add any other central term to $V$ (the project adds $-\epsilon/r^3$) and $\mathbf{A}$ is no longer conserved, the orbit precesses, and you can watch $\mathbf{A}$ turn while $E$ and $\mathbf{L}$ stay put. The same hidden symmetry is why hydrogen's energy levels depend only on $n$ and not on $\ell$: the quantum Runge–Lenz operator commutes with $\hat H$ and connects states of different $\ell$ at the same $n$ (QM-6). Until CM-8's Poisson brackets, "hidden symmetry" is a label, not an explanation; take it as that for now.

### 3.9 Symmetry and degeneracy (preview)

A symmetry does more than conserve: it forces *coincidences*. In CM-6 you will solve the generalized eigenproblem $K\mathbf{a} = \omega^2 M\mathbf{a}$ for normal modes, and whenever the configuration of masses and springs has a symmetry (three equal masses on an equilateral triangle, a square membrane) two or more modes share the same frequency: a rotation of the configuration maps a mode to a different mode with the same $\omega$, and by LA-8 the eigenvectors of a degenerate eigenvalue span an invariant subspace on which the symmetry acts. In QM-5 the same logic, with $[\hat H, \hat L_i] = 0$, forces the $2\ell + 1$ states of given $\ell$ to share one energy; in QM-6 the hidden Runge–Lenz symmetry raises that to the $n^2$ degeneracy of hydrogen. Degeneracy is what symmetry looks like in a spectrum.

### 3.10 Noether in field theory (one paragraph)

Replace $q_i(t)$ by a field $\phi(\mathbf{x}, t)$ and $L$ by a density $\mathcal{L}(\phi, \partial_\mu\phi)$. The same three-line proof produces not a conserved number but a conserved *current* $j^\mu = (\rho, \mathbf{j})$ with $\partial_\mu j^\mu = 0$, i.e. $\partial\rho/\partial t + \nabla\cdot\mathbf{j} = 0$: the charge inside any region changes only by flow through its boundary. The continuity equation for electric charge in EM-6 is Noether's theorem for the phase symmetry of the matter fields; Poynting's theorem (energy flow) and the Maxwell stress tensor (momentum flow) in EM-6 are Noether's theorem for translations in time and space. Every conservation law in this curriculum from here on is an instance of the same argument.

### 3.11 Discrete symmetries (one sentence)

Reflections and time reversal are symmetries without a continuous parameter $\varepsilon$ to differentiate with respect to, so they give no conservation law in classical mechanics; in quantum mechanics they give multiplicative quantum numbers (parity) and selection rules instead (QM-5, QM-8).

### 3.12 The table

> | Symmetry of $L$ | Generator $K$ | Conserved quantity | Where derived |
> |---|---|---|---|
> | Translation along $\hat{\mathbf{n}}$ | $K_a = \hat{\mathbf{n}}$ | $\hat{\mathbf{n}}\cdot\mathbf{P}$, total momentum | §3.3a |
> | Rotation about $\hat{\mathbf{n}}$ | $K_a = \hat{\mathbf{n}}\times\mathbf{r}_a$ | $\hat{\mathbf{n}}\cdot\mathbf{L}$, angular momentum | §3.3b, Ex. 4.1 |
> | Galilean boost along $\mathbf{v}$ (quasi-symmetry, $F = M\mathbf{v}\cdot\mathbf{R}$) | $K_a = \mathbf{v}t$ | $\mathbf{v}\cdot(\mathbf{P}t - M\mathbf{R})$: uniform motion of the center of mass | §3.3c, Ex. 4.2 |
> | Time translation ($\partial L/\partial t = 0$; quasi-symmetry, $F = L$) | $K_i = \dot q_i$ | $h = \sum\dot q_i p_i - L$ (equals $T + V$ iff constraints are fixed) | §3.4, Ex. 4.3 |
> | Cyclic coordinate $q_j$ | $K = \mathbf{e}_j$ | $p_j$ | §3.1 |
> | Scaling (symmetry of equations, *not* of $L$) | — | nothing conserved; instead similarity $t\propto l^{1-k/2}$ | §3.5–3.6 |
> | Hidden (phase-space) symmetry of $1/r$ | moves $\mathbf{r}$ and $\mathbf{p}$ | $\mathbf{A} = \mathbf{p}\times\mathbf{L} - m\kappa\hat{\mathbf{r}}$ | §3.8, Ex. 4.6 |
> | Discrete (parity, time reversal) | no $\varepsilon$ | nothing classically; selection rules in QM | §3.11 |

## 4. Worked examples

### Example 4.1: Angular momentum from rotational invariance, in components

A single particle, $L = \tfrac12 m(\dot x^2 + \dot y^2 + \dot z^2) - V(r)$, $r = \sqrt{x^2 + y^2 + z^2}$. Rotation about $\hat{\mathbf{z}}$: $K = \hat{\mathbf{z}}\times\mathbf{r} = (0, 0, 1)\times(x, y, z) = (-y, x, 0)$. So $\dot K = (-\dot y, \dot x, 0)$.

*Check (3.3).* $\partial L/\partial x = -V'(r)\,\partial r/\partial x = -V'(r)\,x/r$, and likewise for $y$, $z$. Then
$$
\sum_i\frac{\partial L}{\partial q_i}K_i = -\frac{V'(r)}{r}\big[x\cdot(-y) + y\cdot x + z\cdot 0\big] = 0, \qquad \sum_i\frac{\partial L}{\partial\dot q_i}\dot K_i = m\big[\dot x(-\dot y) + \dot y\dot x + 0\big] = 0.
$$
$\delta L = 0$: a symmetry.

*Charge (3.4).* $Q = p_xK_x + p_yK_y + p_zK_z = m\dot x(-y) + m\dot y(x) = m(x\dot y - y\dot x) = L_z$.

*General axis.* $\hat{\mathbf{n}} = (n_1, n_2, n_3)$, $K = \hat{\mathbf{n}}\times\mathbf{r} = (n_2z - n_3y,\ n_3x - n_1z,\ n_1y - n_2x)$. The potential term: $-\tfrac{V'}{r}\big[x(n_2z - n_3y) + y(n_3x - n_1z) + z(n_1y - n_2x)\big]$; every product appears twice with opposite signs ($n_2xz - n_2xz$, etc.), total zero. The kinetic term is the same expression with dots, also zero. The charge:
$$
Q = m\big[\dot x(n_2z - n_3y) + \dot y(n_3x - n_1z) + \dot z(n_1y - n_2x)\big] = n_1\,m(y\dot z - z\dot y) + n_2\,m(z\dot x - x\dot z) + n_3\,m(x\dot y - y\dot x) = \hat{\mathbf{n}}\cdot\mathbf{L}.
$$
Collecting the coefficients of $n_1, n_2, n_3$ is the component form of $\mathbf{p}\cdot(\hat{\mathbf{n}}\times\mathbf{r}) = \hat{\mathbf{n}}\cdot(\mathbf{r}\times\mathbf{p})$. The project does this with symbolic $n_1, n_2, n_3$ and an arbitrary $V$.

### Example 4.2: The Galilean boost and the center-of-mass theorem for two bodies

$L = \tfrac12 m_1|\dot{\mathbf{r}}_1|^2 + \tfrac12 m_2|\dot{\mathbf{r}}_2|^2 - V(|\mathbf{r}_1 - \mathbf{r}_2|)$. Boost: $\mathbf{r}_a \to \mathbf{r}_a + \varepsilon\mathbf{v}t$, $K_1 = K_2 = \mathbf{v}t$, $\dot K_1 = \dot K_2 = \mathbf{v}$.

*The variation.* Let $\mathbf{s} = \mathbf{r}_1 - \mathbf{r}_2$, $s = |\mathbf{s}|$. $\partial L/\partial\mathbf{r}_1 = -V'(s)\,\mathbf{s}/s$ and $\partial L/\partial\mathbf{r}_2 = +V'(s)\,\mathbf{s}/s$. So
$$
\delta L = -V'(s)\frac{\mathbf{s}}{s}\cdot\mathbf{v}t + V'(s)\frac{\mathbf{s}}{s}\cdot\mathbf{v}t + m_1\dot{\mathbf{r}}_1\cdot\mathbf{v} + m_2\dot{\mathbf{r}}_2\cdot\mathbf{v} = \mathbf{v}\cdot(m_1\dot{\mathbf{r}}_1 + m_2\dot{\mathbf{r}}_2).
$$
Not zero. But it equals $\dfrac{d}{dt}\big[\mathbf{v}\cdot(m_1\mathbf{r}_1 + m_2\mathbf{r}_2)\big]$, so $F = \mathbf{v}\cdot(m_1\mathbf{r}_1 + m_2\mathbf{r}_2) = M\mathbf{v}\cdot\mathbf{R}$ with $M = m_1 + m_2$.

*The charge (3.5).* $Q = \mathbf{p}_1\cdot\mathbf{v}t + \mathbf{p}_2\cdot\mathbf{v}t - M\mathbf{v}\cdot\mathbf{R} = \mathbf{v}\cdot(\mathbf{P}t - M\mathbf{R})$.

*Check it directly.* $\dfrac{d}{dt}(\mathbf{P}t - M\mathbf{R}) = \dot{\mathbf{P}}t + \mathbf{P} - M\dot{\mathbf{R}}$. Now $M\dot{\mathbf{R}} = m_1\dot{\mathbf{r}}_1 + m_2\dot{\mathbf{r}}_2 = \mathbf{P}$ identically, and $\dot{\mathbf{P}} = -V'\mathbf{s}/s + V'\mathbf{s}/s = \mathbf{0}$. So $\mathbf{G} = \mathbf{P}t - M\mathbf{R}$ is conserved, and $\mathbf{R}(t) = \mathbf{R}(0) + (\mathbf{P}/M)\,t$.

*Numbers.* $m_1 = 1$, $m_2 = 3$ (so $M = 4$), $\mathbf{r}_1(0) = (0, 0, 0)$, $\mathbf{r}_2(0) = (4, 0, 0)$, $\dot{\mathbf{r}}_1(0) = (2, 1, 0)$, $\dot{\mathbf{r}}_2(0) = (0, -1, 0)$. Then $\mathbf{P} = (2, 1, 0) + (0, -3, 0) = (2, -2, 0)$, $\mathbf{R}(0) = (0 + 12, 0, 0)/4 = (3, 0, 0)$, and $\mathbf{G} = -M\mathbf{R}(0) = (-12, 0, 0)$ at $t = 0$. Whatever $V$ is, at $t = 2$ the center of mass is at $\mathbf{R}(2) = (3, 0, 0) + (2, -2, 0)\cdot 2/4 = (4, -1, 0)$. The project's symbolic version (one dimension, arbitrary $V$) constructs $F$ automatically and checks $\dot Q = 0$ on shell.

### Example 4.3: $h$ versus $E$ for the bead on the rotating hoop

Take $m = 1$ kg, $R = 1$ m, $g = 9.81$ m/s², $\omega = 2$ rad/s, released from rest at $\theta_0 = 1$ rad. (Since $\omega^2 = 4 < g/R = 9.81$, the bottom is stable and the bead oscillates about $\theta = 0$; CM-3 found the critical $\omega^2 = g/R$.)

$$
h = \tfrac12\dot\theta^2 - \tfrac12\cdot 4\sin^2\theta - 9.81\cos\theta, \qquad E = \tfrac12\dot\theta^2 + \tfrac12\cdot 4\sin^2\theta - 9.81\cos\theta.
$$

At $t = 0$: $\sin^2 1 = 0.70807$, $\cos 1 = 0.54030$, so $2\sin^2 1 = 1.41615$ and $9.81\cos 1 = 5.30037$. Hence
$$
h = -1.41615 - 5.30037 = -6.7165\ \text{J}, \qquad E(0) = +1.41615 - 5.30037 = -3.8842\ \text{J}.
$$
$h$ stays at $-6.7165$ J forever (the project integrates for 10 s and finds it flat to $10^{-10}$). $E$ does not: each time the bead passes through $\theta = 0$, where $\sin\theta = 0$, $E = h = -6.7165$ J, and each time it reaches $\theta = \pm 1$, $E = -3.8842$ J. So $E$ oscillates between $-6.72$ and $-3.88$ J with the bead's period. The difference is $E - h = mR^2\omega^2\sin^2\theta = 4\sin^2\theta$, which is $2T_0$, twice the kinetic energy of the azimuthal motion the motor imposes.

Where does the energy come from? The bead's angular momentum about the vertical axis is $L_z = mR^2\sin^2\theta\,\omega$, which changes as $\theta$ changes; the motor supplies the torque $\tau = dL_z/dt = 2mR^2\omega\sin\theta\cos\theta\,\dot\theta$, at power $\tau\omega = 2mR^2\omega^2\sin\theta\cos\theta\,\dot\theta$. That is exactly $dE/dt$ computed from the expression above (the project checks the identity symbolically and numerically). Meanwhile $h$ is the energy an observer rotating with the hoop would write down: $\tfrac12 mR^2\dot\theta^2$ plus the effective potential $-mgR\cos\theta - \tfrac12 mR^2\omega^2\sin^2\theta$. In that frame the hoop is at rest, nothing depends on time, and energy is conserved. Both statements are right; they are about different quantities.

### Example 4.4: Kepler's third law from scaling

For $V = -\kappa/r$, $k = -1$ and (3.12) gives $\beta = \alpha^{3/2}$: scale an orbit up by $\alpha$ and the time to traverse it goes up by $\alpha^{3/2}$. Two geometrically similar orbits (same eccentricity) with semi-major axes $a$ and $a'$ therefore have periods with $T'/T = (a'/a)^{3/2}$, i.e. $T^2/a^3$ is the same for both. That is Kepler's third law, and we did not solve a single equation of motion.

Two cautions. First, mechanical similarity only compares *similar* orbits; that $T^2/a^3$ is also independent of eccentricity needs the actual solution (CM-5), where $T = 2\pi\sqrt{ma^3/\kappa}$. Second, $\kappa = GMm$ includes the planet's mass $m$, and so does $T$ in the Lagrangian $\tfrac12 m|\dot{\mathbf{r}}|^2$; the ratio $\kappa/m = GM$ is what sets the scale, so the constant in $T^2 = (4\pi^2/GM)\,a^3$ is the same for every planet around the same star.

*Numbers.* Mars has $a = 1.5237$ AU. Then $T_{\text{Mars}} = 1.5237^{3/2}$ yr. $\sqrt{1.5237} = 1.23438$, times $1.5237$ gives $1.8808$ yr. The observed sidereal period is $686.98$ d $= 1.8808$ yr. The project also integrates four Kepler orbits with $a = 0.5, 1, 2, 4$ (in units $m = \kappa = 1$) and finds $T^2/a^3 = 39.478418 = 4\pi^2$ for all four.

### Example 4.5: The virial theorem, for an oscillator and for a star

**(a) Harmonic oscillator.** $x = A\cos\omega t$, $V = \tfrac12 m\omega^2x^2$. Over a period, $\langle\cos^2\omega t\rangle = \langle\sin^2\omega t\rangle = \tfrac12$. So $\langle T\rangle = \tfrac12 m\omega^2A^2\cdot\tfrac12 = \tfrac14 m\omega^2A^2$ and $\langle V\rangle = \tfrac12 m\omega^2A^2\cdot\tfrac12 = \tfrac14 m\omega^2A^2$. Equal, as (3.13) with $k = 2$ demands. (The project finds $\langle T\rangle = \langle V\rangle = 0.25$ for $m = \omega = A = 1$, and for the quartic oscillator $V = x^4/4$ finds $2\langle T\rangle = 4\langle V\rangle = 0.3333$.)

**(b) The Sun.** Model the Sun as a self-gravitating ball of ideal gas. Assumptions: (i) uniform density (the real Sun is centrally concentrated; its binding energy is about $1.8\,GM^2/R$, three times the uniform-sphere value, which would triple the estimate below); (ii) the gas is fully ionized hydrogen and helium, mean particle mass $\bar m = \mu m_p$ with $\mu \approx 0.6$; (iii) classical ideal gas, so the thermal kinetic energy is $\tfrac32 Nk_BT$; (iv) the star is static, so the virial theorem holds for instantaneous values, not only time averages (this is the hydrostatic-equilibrium form of the theorem).

*Gravitational energy of a uniform sphere.* Assemble it shell by shell. When the sphere has radius $r$ its mass is $m(r) = Mr^3/R^3$; bringing in the next shell $dm = (3Mr^2/R^3)\,dr$ from infinity releases $Gm(r)\,dm/r$. So
$$
V = -\int_0^R \frac{G}{r}\cdot\frac{Mr^3}{R^3}\cdot\frac{3Mr^2}{R^3}\,dr = -\frac{3GM^2}{R^6}\int_0^R r^4\,dr = -\frac{3}{5}\frac{GM^2}{R}.
$$

*Virial.* $2\langle T\rangle = -\langle V\rangle$ gives $T_{\text{kin}} = \tfrac{3}{10}GM^2/R$. Set this equal to $\tfrac32 Nk_BT$ with $N = M/\bar m$:
$$
k_BT = \frac{2}{3}\cdot\frac{3}{10}\frac{GM^2}{R}\cdot\frac{\bar m}{M} = \frac{GM\bar m}{5R}.
$$

*Numbers.* $G = 6.674\times10^{-11}$, $M = 1.989\times10^{30}$ kg, $R = 6.957\times10^8$ m, $m_p = 1.6726\times10^{-27}$ kg, $k_B = 1.3807\times10^{-23}$ J/K. $GM/R = 1.908\times10^{11}$ J/kg; times $\bar m = 1.004\times10^{-27}$ kg gives $1.915\times10^{-16}$ J; divide by $5k_B$: $T \approx 2.8\times10^6$ K. The Sun's central temperature is $1.57\times10^7$ K and its mass-averaged temperature is a few million kelvin, so a one-line estimate from $M$ and $R$ alone lands within a factor of a few of the truth. The negative total energy $E = -\langle T\rangle$ is why a star contracting and radiating gets *hotter*, which is how it reaches ignition in the first place.

### Example 4.6: The Laplace–Runge–Lenz vector is conserved

$\mathbf{A} = \mathbf{p}\times\mathbf{L} - m\kappa\hat{\mathbf{r}}$ with $\dot{\mathbf{p}} = -\kappa\hat{\mathbf{r}}/r^2 = -\kappa\mathbf{r}/r^3$ and $\mathbf{L} = \mathbf{r}\times\mathbf{p}$ constant.

1. $\dfrac{d\mathbf{A}}{dt} = \dot{\mathbf{p}}\times\mathbf{L} + \mathbf{p}\times\dot{\mathbf{L}} - m\kappa\dfrac{d\hat{\mathbf{r}}}{dt} = \dot{\mathbf{p}}\times\mathbf{L} - m\kappa\dfrac{d\hat{\mathbf{r}}}{dt}$.
2. First term, with $\mathbf{L} = m\,\mathbf{r}\times\dot{\mathbf{r}}$ and BAC–CAB, $\mathbf{a}\times(\mathbf{b}\times\mathbf{c}) = \mathbf{b}(\mathbf{a}\cdot\mathbf{c}) - \mathbf{c}(\mathbf{a}\cdot\mathbf{b})$:
$$
\dot{\mathbf{p}}\times\mathbf{L} = -\frac{\kappa m}{r^3}\,\mathbf{r}\times(\mathbf{r}\times\dot{\mathbf{r}}) = -\frac{\kappa m}{r^3}\big[\mathbf{r}(\mathbf{r}\cdot\dot{\mathbf{r}}) - \dot{\mathbf{r}}\,r^2\big] = \kappa m\left[\frac{\dot{\mathbf{r}}}{r} - \frac{\mathbf{r}(\mathbf{r}\cdot\dot{\mathbf{r}})}{r^3}\right].
$$
3. Second term: $\dfrac{d}{dt}\dfrac{\mathbf{r}}{r} = \dfrac{\dot{\mathbf{r}}}{r} - \dfrac{\mathbf{r}\dot r}{r^2}$, and $\dot r = \dfrac{d}{dt}\sqrt{\mathbf{r}\cdot\mathbf{r}} = \dfrac{\mathbf{r}\cdot\dot{\mathbf{r}}}{r}$, so $\dfrac{d\hat{\mathbf{r}}}{dt} = \dfrac{\dot{\mathbf{r}}}{r} - \dfrac{\mathbf{r}(\mathbf{r}\cdot\dot{\mathbf{r}})}{r^3}$.
4. Steps 2 and 3 are the same vector, so $d\mathbf{A}/dt = \kappa m[\cdots] - m\kappa[\cdots] = \mathbf{0}$. ∎

*What $\mathbf{A}$ tells you.* $\mathbf{A}\cdot\mathbf{L} = (\mathbf{p}\times\mathbf{L})\cdot\mathbf{L} - m\kappa\hat{\mathbf{r}}\cdot\mathbf{L} = 0 - 0$ (the second because $\mathbf{L}\perp\mathbf{r}$), so $\mathbf{A}$ lies in the orbital plane. Dot with $\mathbf{r}$: $\mathbf{A}\cdot\mathbf{r} = \mathbf{r}\cdot(\mathbf{p}\times\mathbf{L}) - m\kappa r = \mathbf{L}\cdot(\mathbf{r}\times\mathbf{p}) - m\kappa r = L^2 - m\kappa r$. Writing $\mathbf{A}\cdot\mathbf{r} = Ar\cos\theta$ with $\theta$ measured from $\mathbf{A}$:
$$
r = \frac{L^2/m\kappa}{1 + (A/m\kappa)\cos\theta},
$$
the conic of CM-5 with eccentricity $e = A/(m\kappa)$, pericenter at $\theta = 0$, i.e. along $\mathbf{A}$. Squaring $\mathbf{A}$ (exercise B7) gives $A^2 = m^2\kappa^2 + 2mEL^2$, the second relation.

*Numbers.* The project integrates the orbit with $m = \kappa = 1$, $a = 1$, $e = 0.5$, tilted so every component of $\mathbf{L}$ is nonzero, for ten periods: $E = -0.5 = -\kappa/2a$, $|\mathbf{L}| = \sqrt{0.75} = 0.8660 = \sqrt{m\kappa a(1 - e^2)}$, $|\mathbf{A}| = 0.5 = m\kappa e$, all flat to about $10^{-10}$. With $-\epsilon/r^3$ added ($\epsilon = 0.002$), $E$ and $\mathbf{L}$ stay flat while the pericenter advances $3.95°$ per orbit; first-order perturbation theory (Landau §15, problem 3) predicts $6\pi\kappa\epsilon m^2/L^4 = 3.84°$.

## 5. Why it matters and where it connects

- **CM-5 (central forces).** $p_\theta = \ell$ (Theorem 3.1) and $h = E$ (Theorem 3.7) reduce the two-body problem to $\tfrac12\mu\dot r^2 + V_{\text{eff}}(r) = E$ with $V_{\text{eff}} = V + \ell^2/2\mu r^2$. The Runge–Lenz vector (Example 4.6) hands you the orbit equation without integrating; mechanical similarity (Theorem 3.8) hands you Kepler's third law.
- **CM-6 (normal modes).** Symmetric configurations have degenerate frequencies in $K\mathbf{a} = \omega^2M\mathbf{a}$ (§3.9); amplitude-independence of the period is the $k = 2$ case of (3.12), which is why small oscillations are linear.
- **CM-8 (Hamiltonian mechanics).** $h$ written in terms of $p$ is $H$; (3.9) becomes $dH/dt = \partial H/\partial t$. "$Q$ is conserved" becomes the Poisson bracket $\{Q, H\} = 0$, and *every* conserved $Q$, including the quadratic ones, generates a symmetry of phase space via $\{\,\cdot\,, Q\}$. That is where "hidden symmetry" becomes a theorem: the components of $\mathbf{L}$ and $\mathbf{A}$ close under Poisson brackets into the Lie algebra of 4-dimensional rotations.
- **QM-5, QM-6, QM-9.** $[\hat H, \hat L_z] = 0$ is Theorem 3.3 with commutators in place of Poisson brackets; the $2\ell + 1$ degeneracy is §3.9; hydrogen's $n^2$ degeneracy is §3.8; Ehrenfest's theorem in QM-9 is how Noether charges keep being conserved as expectation values.
- **EM-6.** Charge conservation $\partial_t\rho + \nabla\cdot\mathbf{j} = 0$ and Poynting's theorem are §3.10.
- **SM (statistical mechanics).** The virial theorem (3.13) is the classical tool for stellar structure and for the equation of state of interacting gases (the virial expansion is named for the same $\sum\mathbf{F}\cdot\mathbf{r}$); ergodic time averages are what $\langle\cdot\rangle$ means there.
- **F6.** Theorem 3.6 *is* the Beltrami identity.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: cyclic coordinates, computing $h$ for three Lagrangians, Euler's theorem, scaling exponents, checking $\delta L$ by hand, time translation as a generator. B1–B8 are the module: B1 (uniform gravity: one translation becomes a quasi-symmetry, two rotations seem lost, and the free-fall picture shows what really happened to them), B3 (charged particle in a magnetic field: the conserved momentum is neither $m\dot{\mathbf{r}}$ nor $\mathbf{p}$), B5 (bead on a rotating rod, $h \ne E$ with explicit solution), B7 ($A^2 = m^2\kappa^2 + 2mEL^2$) and B8 (extend `project.py`) are essential. C1–C4: the $SU(2)$ charges of the isotropic oscillator, a careful virial argument with negative heat capacity, why $V \propto -1/r^2$ has no stable orbits, and the damped oscillator's hidden conserved energy.

## 7. Coding project

`project.py` has two halves.

1. **Noether in SymPy.** `delta_L(L, q, v, K)` evaluates the invariance condition (3.3); `noether_charge` returns $Q = \sum p_iK_i$ or raises if $\delta L \ne 0$; `noether_charge_quasi` constructs $F$ from $\delta L$ by the Poincaré-lemma integral $F = \int_0^1[\sum_i a_i(sq, st)q_i + b(sq, st)t]\,ds$ and returns $Q = \sum p_iK_i - F$; `energy_function` builds $h$; `el_accelerations` solves the Euler–Lagrange equations for $\ddot q$ via the Hessian $W$; `on_shell_dt` differentiates any expression along solutions. Applied to: polar central force ($p_\theta$), two-body translation ($P$) and boost ($Pt - MX$, with $F$ found automatically), rotation about $\hat{\mathbf{z}}$ and about a symbolic $\hat{\mathbf{n}}$ for a 3D particle in an arbitrary central potential, the rotating hoop ($dh/dt = 0$, $dE/dt = $ motor power, time translation with $K = \dot\theta$), and the driven oscillator ($dh/dt = -\partial L/\partial t$).
2. **Numerical conservation monitoring.** Kepler for ten periods: $E$, $\mathbf{L}$, $\mathbf{A}$ flat to $10^{-10}$ (asserted). Add $-\epsilon/r^3$: $E$, $\mathbf{L}$ flat, $\mathbf{A}$ rotates, orbit precesses; the advance per orbit is measured from pericenter passages and compared to first-order perturbation theory. Rotating hoop: $h$ flat, $E$ not, $E - h = mR^2\omega^2\sin^2\theta$. Mechanical similarity: $T^2/a^3$ for four orbit sizes; period versus amplitude for $k = 2$ and $k = 4$ by quadrature with log–log slopes $0$ and $-1$, cross-checked against an ODE integration. Virial: time averages over one period for Kepler, the quartic and the harmonic oscillator; the Sun's temperature.

**Build it yourself:** before reading Part 1, write `delta_L` and `noether_charge` from (3.3)–(3.4) (about ten lines each) and test them on the polar central force. Then write the pericenter detector for Part 3 (an event on $\mathbf{r}\cdot\mathbf{v}$ changing sign from negative to positive). "Your turn" at the end asks for the isotropic-oscillator charges, the double pendulum (which conserves only $h$), the $1/r^2$ scaling, and the dependence of the precession rate on $\epsilon$. The final section, "Reference checks for the exercises", re-derives every symbolic claim and prints every number that appears in `solutions.md`; do each exercise on paper before reading it.

## 8. Where you'll struggle

- **You will remember "symmetry implies conservation" and forget both the converse and the hypothesis.** The hypothesis is invariance of $L$, for all paths, up to a total derivative $dF/dt$ with $F = F(q, t)$. It is not "the physics looks the same." Scaling leaves the free particle's physics the same and gives no conservation law (§3.5). The converse (Theorem 3.5) is what makes the theorem useful in reverse: if you *find* something conserved, there is a symmetry behind it, and that is how the Runge–Lenz vector eventually got explained. Say the full sentence: *a continuous quasi-symmetry of $L$ gives a conserved charge, and a conserved charge linear in the momenta gives a quasi-symmetry of $L$.*
- **$h$ versus $E$ only matters when constraints move, and that is exactly when you need it.** For every system with fixed constraints $h = T + V$ and you can stop reading. The moment there is a rotating frame, a driven support, or a time-dependent constraint, $h \ne E$, and students write down $T + V$, find it is not conserved, and conclude the problem is dissipative. It is not. Compute $h$ from (3.8), every time, and interpret the difference $E - h = T_1 + 2T_0$ as the work done by whatever is moving the constraint (Example 4.3).
- **The Galilean boost needs the total-derivative version, and most books skip it.** $\delta L = \mathbf{v}\cdot\mathbf{P} \ne 0$; if you stop there you conclude boosts are not symmetries, which is false (Galilean relativity is the oldest symmetry in mechanics). You must recognize $\mathbf{v}\cdot\mathbf{P}$ as $\tfrac{d}{dt}(M\mathbf{v}\cdot\mathbf{R})$ and subtract $F$. Exercise B1 (uniform gravity) and B3 (magnetic field) are two more cases where the $F$ term is unavoidable; do them until subtracting $F$ is reflex.
- **"Hidden symmetry" explains nothing until CM-8.** The phrase is a promissory note. What you can say now: $\mathbf{A}$ is conserved (Example 4.6), it is quadratic in $\mathbf{p}$ so no transformation of the apparatus generates it (Theorem 3.5's remark), and its existence is why Kepler orbits close. What you cannot say yet is what the symmetry *is*. Poisson brackets make it precise; until then do not pretend the word explains anything.
- **Sign and order in the rotation generator.** $K = \hat{\mathbf{n}}\times\mathbf{r}$, not $\mathbf{r}\times\hat{\mathbf{n}}$; the charge is $\mathbf{p}\cdot(\hat{\mathbf{n}}\times\mathbf{r}) = \hat{\mathbf{n}}\cdot(\mathbf{r}\times\mathbf{p})$, which is $+\hat{\mathbf{n}}\cdot\mathbf{L}$. Get one cross product backwards and you conserve $-L_z$, which is harmless, but get the triple-product permutation wrong and you will "prove" $\delta L \ne 0$ for a rotation. Example 4.1 in components is the antidote.
- **On shell versus off shell.** The invariance condition (3.3) must hold for *every* $(q, \dot q, t)$; conservation holds only *on solutions*. Students check $\delta L = 0$ along one solution and call it a symmetry, or demand $dQ/dt = 0$ identically and conclude nothing is ever conserved. Step 2 of Theorem 3.3 is the only place the equations of motion enter; locate it and keep the two levels apart.
- **The virial theorem is a statement about time averages of bounded motion.** Applying $2T = -V$ at an instant of a Kepler orbit is wrong (at pericenter $2T > -V$, at apocenter $2T < -V$; the project's running averages show the approach). For a static star the instantaneous and averaged forms coincide, which is why Example 4.5(b) is allowed to use it that way; say that assumption out loud when you do.
- **Exercise B3 is harder than it looks.** The canonical momentum $\mathbf{p} = m\dot{\mathbf{r}} + q\mathbf{A}$ is not $m\dot{\mathbf{r}}$, the Noether charge for translation is $p_x$ minus the $F$ correction, and what ends up conserved ($m\dot x - qBy$) is neither the kinetic nor the canonical momentum. Keep the three apart.

## 9. Check yourself

1. State Theorem 3.3 precisely, including what "invariant" means and where the equations of motion enter the proof.
2. Why does $Q = \sum p_iK_i$ need the correction $-F$ for a Galilean boost? What is $F$?
3. For $L = \tfrac12 m(\dot r^2 + r^2\dot\theta^2) - V(r)$, which coordinate is cyclic and what is conserved? Rewrite in Cartesian coordinates: is anything cyclic? Is anything still conserved?
4. Write the generator $K$ of a rotation about $\hat{\mathbf{n}}$ and show in one line that the charge is $\hat{\mathbf{n}}\cdot\mathbf{L}$.
5. Derive $dh/dt = -\partial L/\partial t$ in three lines.
6. When is $h = T + V$? Give the hypothesis and the one-line proof using Euler's theorem.
7. For the bead on the rotating hoop, which of $h$, $E$ is conserved, and what physical system supplies or absorbs the difference?
8. Scaling $x \to \lambda x$ is a symmetry of the free particle's equation of motion. Why does Noether's theorem give nothing?
9. From $V$ homogeneous of degree $k$, how do period, velocity and energy scale with size? Get Kepler's third law and the amplitude-independence of the harmonic period from this.
10. State the virial theorem with its hypothesis. What does it say for a bound gravitational system's total energy?
11. What is $\mathbf{A}$? Give two of its properties and say in one sentence why its conservation does not come from a point transformation.
12. Why does a symmetry produce degenerate eigenvalues (one sentence, CM-6 or QM-5)?

## 10. Pair with

- **Taylor, *Classical Mechanics*, §7.8–7.9.** The cleanest textbook treatment of cyclic coordinates and of $h$ versus $E$, with the same rotating examples as CM-3. Taylor does not state Noether's theorem as such; read him for the energy-function discussion.
- **Goldstein, Poole & Safko, §2.6–2.7.** The source of the energy function $h$ and the notation used here, with the homogeneous-quadratic argument for $h = T + V$ done carefully. Denser than these notes; good for a second pass.
- **Neuenschwander, *Emmy Noether's Wonderful Theorem*, ch. 1–4.** A whole short book on this module. Chapter 1 has the history (including what Hilbert and Klein actually asked), chapters 2–4 build the theorem with many worked examples, and it is the only book at this level that states the converse. Read it alongside §3.2.
- **Landau & Lifshitz, *Mechanics*, §6–10.** Seven pages that contain §3.3, §3.4 and §3.6 of these notes with nothing wasted. §10 (mechanical similarity) and the virial theorem are done there in half a page each; read them after you have done the derivations here and see how much they compress.
- **John Baez, "Noether's theorem in a nutshell"** (online essay). The one-paragraph proof of Theorem 3.3, written to make you see that it is only the product rule plus Euler–Lagrange. Read it after §3.2 and check that you could reproduce it on a napkin.
