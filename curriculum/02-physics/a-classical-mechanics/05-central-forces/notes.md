# CM-5: Central Forces and the Two-Body Problem

**Time:** ~10 hours · **Prerequisites:** CM-4 · **Unlocks:** QM-6 (hydrogen), the virial theorem in statistical mechanics (SM), QM-1 (the Rutherford atom)

> **Conventions in this module.** Two bodies have masses $m_1, m_2$, total mass $M = m_1 + m_2$, positions $\mathbf r_1, \mathbf r_2$. The **relative coordinate** is $\mathbf r = \mathbf r_1 - \mathbf r_2$ (it points *from* body 2 *to* body 1), $r = |\mathbf r|$, $\hat{\mathbf r} = \mathbf r/r$, and the center of mass is $\mathbf R = (m_1\mathbf r_1 + m_2\mathbf r_2)/M$. The **reduced mass** is $\mu = m_1m_2/(m_1+m_2)$. Angular momentum of the relative motion is $\mathbf L = \mu\,\mathbf r\times\dot{\mathbf r}$ with magnitude $\ell = |\mathbf L|$. The force on the reduced particle is $F(r)\,\hat{\mathbf r}$ with $F = -dV/dr$, so **$F < 0$ is attractive**. Gravity: $V = -k/r$ with $k = Gm_1m_2 > 0$. $E$ is the energy of the *relative* motion, so **bound orbits have $E < 0$**. In the orbital plane, $\varphi$ is the polar angle, measured from perihelion unless stated; $u = 1/r$ and a prime means $d/d\varphi$. Orbit constants: eccentricity $e \ge 0$, semi-latus rectum $c = \ell^2/\mu k$, semi-major axis $a$, semi-minor axis $b$, perihelion and aphelion distances $r_p, r_a$ (pericenter and apocenter for a non-solar center). For scattering (§3.11) the repulsive Coulomb potential is written $V = +k/r$ with $k = Z_1Z_2e^2/4\pi\varepsilon_0 > 0$; $\theta$ is the deflection angle ($0 \le \theta \le \pi$, not the polar angle $\varphi$), $b$ is the impact parameter. "Kepler units" means $\mu = k = 1$. CM-4 wrote the force constant as $\kappa$ and the mass as $m$; here they are $k$ and $\mu$, as in Goldstein and Taylor.

## 0. Where this sits

You need CM-4: a cyclic coordinate gives a conserved momentum (CM-4 Theorem 3.1), a time-independent Lagrangian gives a conserved energy function $h$ (CM-4 Theorem 3.7), and the Laplace–Runge–Lenz vector is conserved for $-k/r$ (CM-4 Example 4.6). You also need CM-3's polar-coordinate Lagrangian and its spherical-pendulum reduction (CM-3 Example 4.7), which is this module's method in miniature.

This module uses those tools on the problem they were invented for. Afterwards, QM-6 solves hydrogen by the same reduction with $\ell^2 \to \hbar^2\ell(\ell+1)$; QM-1 starts from Rutherford's nucleus, which this module derives from scattering data; CM-8 recasts the Kepler problem in action-angle variables; CM-6 generalizes the radial oscillation frequency of §3.9 to many degrees of freedom.

## 1. The problem that forced this idea

Tycho Brahe spent twenty years measuring planetary positions to a minute or two of arc. Johannes Kepler inherited the data in 1601 and spent most of a decade on Mars. His best circle-based model missed Tycho's positions by eight minutes of arc, and Kepler would not blame the data: in *Astronomia Nova* (1609) he wrote that those eight minutes "led the way to the reform of all astronomy." He replaced the circle by an ellipse with the Sun at a focus (the first law) and found that the Sun–planet line sweeps equal areas in equal times (the second law, which he actually found first). In *Harmonices Mundi* (1619) he added the third: the square of the period is proportional to the cube of the orbit's size, with the same constant for every planet.

These were rules without a cause. In August 1684 Edmond Halley visited Isaac Newton in Cambridge and asked what curve a planet would follow if the Sun attracted it with an inverse-square force. Newton answered "an ellipse" and within months sent Halley the proof, *De motu corporum in gyrum*, which grew into the *Principia* (1687). Book I opens by showing that equal areas follow from *any* force aimed at a fixed center (Proposition 1: angular momentum conservation), then that an ellipse with the center of force at a focus requires an inverse-square force (Proposition 11). Newton also treated two bodies pulling on each other (Book I, Section XI) and found that Kepler's third law is not exact: its constant depends on the sum of the masses. Mechanics as a subject, with $\mathbf F = m\mathbf a$ as its law, exists because it explained Kepler's three sentences.

The two-body problem is also the only interacting gravitational problem solvable in closed form for all initial conditions. Three bodies defeated Euler, Lagrange and Laplace, and Poincaré showed around 1890 that no general solution of the kind derived here exists. So the two-body solution became the template: use the symmetries (CM-4) to strip away every coordinate you can, then solve the one-dimensional problem that is left. Schrödinger did exactly this for hydrogen in 1926 (QM-6).

The same orbits returned in 1911 with a different job. Hans Geiger and Ernest Marsden, in Rutherford's Manchester laboratory, had reported in 1909 that about one alpha particle in 8000 bounced back from a platinum target, which J. J. Thomson's atom (positive charge spread over the whole atom) could not do. Rutherford computed the hyperbolic orbit of an alpha particle near a point charge, derived the angular distribution $\propto\sin^{-4}(\theta/2)$, and concluded that the positive charge sits in a nucleus thousands of times smaller than the atom. Geiger and Marsden confirmed the $\sin^{-4}$ law in 1913 from $5°$ to $150°$. You will derive it in §3.11.

## 2. Intuition

### Two bodies are one body

Sun and Earth both circle their common center of mass, which drifts in a straight line because no outside force acts. The interaction depends only on the vector $\mathbf r$ from Sun to Earth. So there are two independent problems: a free particle (the center of mass) and one particle (the relative vector) in the potential $V(r)$. The relative particle's mass turns out to be $\mu = m_1m_2/(m_1+m_2)$, smaller than either: Earth's mass less three parts per million for Sun–Earth, half of either for equal masses.

### Angular momentum flattens the problem twice

A central force pulls along $\mathbf r$, so it exerts no torque about the center: $\mathbf L$ is a fixed vector, and since $\mathbf L \perp \mathbf r$ at every instant, $\mathbf r$ stays in the fixed plane perpendicular to $\mathbf L$. Three dimensions become two. In that plane the polar angle $\varphi$ is cyclic, so $\ell = \mu r^2\dot\varphi$ is constant and $\varphi$ just follows $r$. Two dimensions become one: a particle on a half-line with energy $\tfrac12\mu\dot r^2 + V_{\rm eff}(r)$, $V_{\rm eff} = V + \ell^2/2\mu r^2$. The extra term is a repulsive **centrifugal barrier**, the kinetic energy of the angular motion, growing as $1/r^2$ because $\ell$ is fixed. CM-3's spherical pendulum had the same barrier with $\sin\theta$ in place of $r$.

### Reading the orbit off $V_{\rm eff}$

For gravity, $V_{\rm eff} = -k/r + \ell^2/2\mu r^2$. The barrier wins at small $r$, the attraction wins at large $r$, and between them is a well:

```
 V_eff
   |  |
   |  \
   |   \                                     E > 0 : in from infinity, bounce, out
 + |----*---------------------------------   (one turning point)
   |     \
 0 +------*-------------------------------   E = 0 : just barely escapes
   |       \              _______-------       (V_eff -> 0 from below like -k/r)
   |        \       __----
 - |---------*-----*----------------------   E < 0 : trapped between r_p and r_a
   |          \__ /
   |             *                           E = V_min : sits at the bottom r0
   +------------------------------------------> r
              r_p  r0     r_a
```

At $E = V_{\min}$ the radius never changes: a **circle**. For $V_{\min} < E < 0$ the radius oscillates between turning points $r_p$ and $r_a$ while $\varphi$ advances: a bounded orbit. For $E \ge 0$ the particle comes in from infinity, bounces off the barrier once, and leaves. Two quantitative guesses need almost no algebra.

**Energy fixes the size.** The turning points solve $E = -k/r + \ell^2/2\mu r^2$, that is $Er^2 + kr - \ell^2/2\mu = 0$. The sum of the roots is $r_p + r_a = -k/E$: the orbit's long diameter depends on $E$ alone, not on $\ell$. Call half of it $a$; then $E = -k/2a$.

**The bounded orbit closes.** Near the bottom of the well, $r$ oscillates with angular frequency $\omega_r = \sqrt{V_{\rm eff}''(r_0)/\mu}$ while $\varphi$ advances at $\omega_\varphi = \ell/\mu r_0^2$. For $-k/r$, the minimum is at $r_0 = \ell^2/\mu k$, and $V_{\rm eff}''(r_0) = 3\ell^2/\mu r_0^4 - 2k/r_0^3 = k/r_0^3$. So $\omega_r^2 = k/\mu r_0^3$, while $\omega_\varphi^2 = \ell^2/\mu^2r_0^4 = k/\mu r_0^3$. They are *equal*: the radius goes in and out exactly once per revolution, so the orbit closes, with one closest point per turn, which is what an ellipse with the force center at a focus does. For a force $\propto r^n$ the same computation gives $\omega_r/\omega_\varphi = \sqrt{n+3}$, irrational for almost every $n$.

So guess, before §3 derives it: bound inverse-square orbits are ellipses with the force center at a focus; their energy and (by CM-4 Example 4.4's scaling) their period depend only on $a$; circular orbits are stable only for forces falling more slowly than $1/r^3$; and in scattering, closer passes mean bigger deflections, so the fraction that bounces back measures how compact the target's charge is.

## 3. Formal development

### 3.1 The two-body Lagrangian and the reduced mass

Two particles interacting through a potential that depends only on their separation have

$$
L = \tfrac12 m_1|\dot{\mathbf r}_1|^2 + \tfrac12 m_2|\dot{\mathbf r}_2|^2 - V(|\mathbf r_1 - \mathbf r_2|). \tag{3.1}
$$

**Theorem 3.1 (Two-body reduction).** In the coordinates $\mathbf R, \mathbf r$, (3.1) becomes $L = \tfrac12 M|\dot{\mathbf R}|^2 + \tfrac12\mu|\dot{\mathbf r}|^2 - V(r)$. The center of mass moves uniformly and the relative coordinate obeys the equation of a single particle of mass $\mu$ in the potential $V$.

▶ **Try it first.** Solve the definitions of $\mathbf R$ and $\mathbf r$ for $\mathbf r_1$ and $\mathbf r_2$, substitute into the kinetic energy, and watch the cross term.

**Derivation.**

1. From $M\mathbf R = m_1\mathbf r_1 + m_2\mathbf r_2$ and $\mathbf r = \mathbf r_1 - \mathbf r_2$: substitute $\mathbf r_2 = \mathbf r_1 - \mathbf r$ into the first to get $M\mathbf R = M\mathbf r_1 - m_2\mathbf r$, so
$$
\mathbf r_1 = \mathbf R + \frac{m_2}{M}\mathbf r, \qquad \mathbf r_2 = \mathbf R - \frac{m_1}{M}\mathbf r. \tag{3.2}
$$
2. Differentiate and square:
$$
\tfrac12 m_1\Big|\dot{\mathbf R} + \tfrac{m_2}{M}\dot{\mathbf r}\Big|^2 + \tfrac12 m_2\Big|\dot{\mathbf R} - \tfrac{m_1}{M}\dot{\mathbf r}\Big|^2
= \tfrac12(m_1 + m_2)|\dot{\mathbf R}|^2 + \Big(\tfrac{m_1m_2}{M} - \tfrac{m_2m_1}{M}\Big)\dot{\mathbf R}\cdot\dot{\mathbf r} + \tfrac12\,\frac{m_1m_2^2 + m_2m_1^2}{M^2}|\dot{\mathbf r}|^2.
$$
3. The cross term cancels exactly. That is what choosing the *mass-weighted* center buys. The last coefficient is $m_1m_2(m_1 + m_2)/M^2 = m_1m_2/M = \mu$. So
$$
L = \tfrac12 M|\dot{\mathbf R}|^2 + \tfrac12\mu|\dot{\mathbf r}|^2 - V(r). \tag{3.3}
$$
4. $\mathbf R$ is cyclic, so $M\dot{\mathbf R}$ is conserved (CM-4 Theorem 3.1). The Euler–Lagrange equations for $\mathbf r$ give $\mu\ddot{\mathbf r} = -\nabla V(r) = -V'(r)\hat{\mathbf r}$:
$$
M\ddot{\mathbf R} = \mathbf 0, \qquad \mu\ddot{\mathbf r} = F(r)\,\hat{\mathbf r}, \quad F = -\frac{dV}{dr}. \tag{3.4}
$$
∎

■ **What just happened.** Exactly, with no approximation, the two-body problem split into a free particle of mass $M$ and a central-force problem for mass $\mu$. The orbit of $\mathbf r$ is body 1 *as seen from body 2*; by (3.2) each body's orbit about the center of mass is the same curve scaled by $m_2/M$ or $-m_1/M$. The project integrates the full eight-dimensional problem (Part 1): $\mathbf R(t)$ is a straight line to $10^{-9}$ and $\mathbf r(t)$ is the predicted conic. (CM-4 Exercise B2 reached (3.3) by counting conserved charges.)

### 3.2 Motion in a plane, and Kepler's second law

**Theorem 3.2 (Planarity and the area law).** For (3.4), $\mathbf L = \mu\,\mathbf r\times\dot{\mathbf r}$ is constant. If $\mathbf L \ne \mathbf 0$, the motion stays in the plane through the origin perpendicular to $\mathbf L$; if $\mathbf L = \mathbf 0$, it stays on a line through the origin. The area swept by $\mathbf r$ grows at the constant rate $dA/dt = \ell/2\mu$.

▶ **Try it first.** Differentiate $\mathbf L$ and use (3.4). Then compute $\mathbf r\cdot\mathbf L$.

**Derivation.**

1. $\dot{\mathbf L} = \mu\,\dot{\mathbf r}\times\dot{\mathbf r} + \mu\,\mathbf r\times\ddot{\mathbf r} = \mathbf 0 + \mathbf r\times F(r)\hat{\mathbf r} = \mathbf 0$: both cross products have parallel factors. (This is CM-4's Noether charge for rotations.)
2. $\mathbf r\cdot\mathbf L = \mu\,\mathbf r\cdot(\mathbf r\times\dot{\mathbf r}) = 0$ (a triple product with a repeated vector). So $\mathbf r$ is always perpendicular to the *fixed* vector $\mathbf L$: it lies in the plane $\{\mathbf x : \mathbf x\cdot\mathbf L = 0\}$.
3. If $\mathbf L = \mathbf 0$ then $\dot{\mathbf r}\parallel\mathbf r$ always, so $\hat{\mathbf r}$ never changes: radial motion on a line.
4. Area. In time $dt$ the radius vector sweeps a thin triangle of area $dA = \tfrac12|\mathbf r\times d\mathbf r| = \tfrac12|\mathbf r\times\dot{\mathbf r}|\,dt = (\ell/2\mu)\,dt$.
5. In plane polar coordinates, $\dot{\mathbf r} = \dot r\hat{\mathbf r} + r\dot\varphi\hat{\boldsymbol\varphi}$ and $\mathbf r\times\dot{\mathbf r} = r^2\dot\varphi\,\hat{\mathbf z}$, so
$$
\mathbf L = \mu\,\mathbf r\times\dot{\mathbf r}, \qquad \ell = \mu r^2\dot\varphi, \tag{3.5}
$$
$$
\frac{dA}{dt} = \frac{\ell}{2\mu} = \tfrac12 r^2\dot\varphi. \tag{3.6}
$$
∎

■ **What just happened.** Kepler's second law is angular momentum conservation, true for every central force, which is why Newton could put it first. The direction of $\mathbf L$ fixes the plane, its magnitude the sweep rate.

### 3.3 The radial equation and the effective potential

In the orbital plane,

$$
L = \tfrac12\mu(\dot r^2 + r^2\dot\varphi^2) - V(r).
$$

**Theorem 3.3 (Radial equation).** $p_\varphi = \partial L/\partial\dot\varphi = \mu r^2\dot\varphi = \ell$ is conserved, and

$$
\mu\ddot r = -\frac{dV}{dr} + \frac{\ell^2}{\mu r^3}, \tag{3.7}
$$

$$
E = \tfrac12\mu\dot r^2 + V_{\rm eff}(r), \qquad V_{\rm eff}(r) = V(r) + \frac{\ell^2}{2\mu r^2}. \tag{3.8}
$$

**Derivation.**

1. $\varphi$ does not appear in $L$, so $\dot p_\varphi = \partial L/\partial\varphi = 0$: $p_\varphi = \mu r^2\dot\varphi$ is constant, and by (3.5) it is $\ell$.
2. The $r$ equation: $\partial L/\partial\dot r = \mu\dot r$ and $\partial L/\partial r = \mu r\dot\varphi^2 - V'(r)$, so $\mu\ddot r = \mu r\dot\varphi^2 - V'(r)$. Substitute $\dot\varphi = \ell/\mu r^2$ *in the equation of motion*: $\mu r\dot\varphi^2 = \ell^2/\mu r^3$. That is (3.7).
3. $L$ has no explicit $t$ and $T$ is quadratic in velocities, so $E = T + V$ is conserved (CM-4 Theorem 3.7): $E = \tfrac12\mu\dot r^2 + \tfrac12\mu r^2\dot\varphi^2 + V = \tfrac12\mu\dot r^2 + \ell^2/2\mu r^2 + V$, which is (3.8). As a check, multiply (3.7) by $\dot r$: $\mu\dot r\ddot r = -\dot r\,dV_{\rm eff}/dr$, i.e. $\tfrac{d}{dt}[\tfrac12\mu\dot r^2 + V_{\rm eff}] = 0$. ∎

■ **What just happened.** $r$ moves like a particle on a line in $V_{\rm eff}$. A warning: substitute $\dot\varphi = \ell/\mu r^2$ into the *equation of motion* or the *energy*, never into the *Lagrangian*. If you do, $L \to \tfrac12\mu\dot r^2 + \ell^2/2\mu r^2 - V$, whose Euler–Lagrange equation $\mu\ddot r = -V' - \ell^2/\mu r^3$ has the centrifugal force with the wrong sign: varying $L$ holds $\dot\varphi$ fixed, not $\ell$. The correct reduced Lagrangian is Routh's, $L - p_\varphi\dot\varphi$ (CM-8). CM-7's spinning top has the same trap.

### 3.4 Qualitative orbits for the inverse-square force

For $V = -k/r$,

$$
V_{\rm eff}' = \frac{k}{r^2} - \frac{\ell^2}{\mu r^3} = 0 \;\Rightarrow\; r_0 = \frac{\ell^2}{\mu k}, \qquad V_{\rm eff}(r_0) = -\frac{\mu k^2}{\ell^2} + \frac{\mu k^2}{2\ell^2} = -\frac{\mu k^2}{2\ell^2}. \tag{3.9}
$$

$V_{\rm eff}''(r_0) = k/r_0^3 > 0$: a minimum, the circular orbit, where $\mu v^2/r_0 = k/r_0^2$. The classification of §2 is now exact:

| Energy | Radial motion | Orbit (proved in §3.6) |
|---|---|---|
| $E = -\mu k^2/2\ell^2$ | $r = r_0$ | circle |
| $-\mu k^2/2\ell^2 < E < 0$ | $r_p \le r \le r_a$ | ellipse |
| $E = 0$ | $r \ge r_p$, reaches $\infty$ with zero speed | parabola |
| $E > 0$ | $r \ge r_p$, reaches $\infty$ with speed $\sqrt{2E/\mu}$ | hyperbola |

The turning points solve $Er^2 + kr - \ell^2/2\mu = 0$, so for $E < 0$, $r_p + r_a = -k/E$ and $r_pr_a = -\ell^2/2\mu E$. Project Part 3 plots $V_{\rm eff}$ with four energy levels and the four orbits.

### 3.5 The orbit equation

Solving (3.8) gives $t(r)$; you want the *shape* $r(\varphi)$. The trick is to use $u = 1/r$ as the dependent variable and $\varphi$ as the independent one. Memorize it; nobody rediscovers it under exam pressure.

**Theorem 3.4 (Binet's equation).** For any central force $F(r)$, the orbit $u(\varphi) = 1/r(\varphi)$ satisfies

$$
\frac{d^2u}{d\varphi^2} + u = -\frac{\mu}{\ell^2u^2}\,F\!\left(\frac1u\right). \tag{3.10}
$$

▶ **Try it first.** Use $\dot\varphi = \ell u^2/\mu$ to turn every $d/dt$ into $d/d\varphi$. Compute $\dot r$, then $\ddot r$, then substitute into (3.7).

**Derivation.**

1. By (3.5), $\dot\varphi = \ell/\mu r^2 = \ell u^2/\mu$. For any quantity $Q$ along the orbit, the chain rule gives $\dot Q = \dfrac{dQ}{d\varphi}\dot\varphi = \dfrac{\ell u^2}{\mu}\dfrac{dQ}{d\varphi}$.
2. First derivative: $r = 1/u$, so $\dfrac{dr}{d\varphi} = -\dfrac{1}{u^2}u'$ and
$$
\dot r = \frac{\ell u^2}{\mu}\cdot\Big(-\frac{u'}{u^2}\Big) = -\frac{\ell}{\mu}u'.
$$
The $u^2$ cancels. That cancellation is the reason for choosing $u = 1/r$.
3. Second derivative: apply step 1 to $Q = \dot r = -(\ell/\mu)u'$, with $\ell, \mu$ constant:
$$
\ddot r = \frac{\ell u^2}{\mu}\cdot\frac{d}{d\varphi}\Big(-\frac{\ell}{\mu}u'\Big) = -\frac{\ell^2u^2}{\mu^2}u''.
$$
4. Substitute into (3.7) written as $\mu\ddot r - \ell^2/\mu r^3 = F$, with $1/r^3 = u^3$:
$$
-\frac{\ell^2u^2}{\mu}u'' - \frac{\ell^2u^3}{\mu} = F(1/u) \;\Longrightarrow\; -\frac{\ell^2u^2}{\mu}\,(u'' + u) = F(1/u).
$$
5. Divide by $-\ell^2u^2/\mu$. That is (3.10). ∎

■ **What just happened.** Steps 2–3 are chain rule, step 4 is the radial equation. For the inverse-square force the right side of (3.10) is a constant, so the orbit equation is a *harmonic oscillator* in $\varphi$ with a constant drive. Project Part 10 repeats steps 2–3 in SymPy.

### 3.6 The Kepler orbit: conic sections

**Theorem 3.5 (Kepler's first law).** For $F = -k/r^2$, every orbit with $\ell \ne 0$ is $r(\varphi) = c/[1 + e\cos(\varphi - \varphi_0)]$ with $c = \ell^2/\mu k$, a conic with eccentricity $e \ge 0$ and a focus at the origin. Bounded orbits ($e < 1$) are ellipses.

**Derivation.**

1. $F(1/u) = -ku^2$, so the right side of (3.10) is $-\dfrac{\mu}{\ell^2u^2}(-ku^2) = \dfrac{\mu k}{\ell^2}$:
$$
u'' + u = \frac{\mu k}{\ell^2}. \tag{3.11}
$$
2. This is linear with constant coefficients. A particular solution is the constant $\mu k/\ell^2$; the homogeneous solutions are $A\cos(\varphi - \varphi_0)$. Choose $\varphi_0$ so that $A \ge 0$ and write $A = e\,\mu k/\ell^2$ with $e \ge 0$:
$$
u = \frac{\mu k}{\ell^2}\big[1 + e\cos(\varphi - \varphi_0)\big] \quad\Longleftrightarrow\quad r(\varphi) = \frac{c}{1 + e\cos(\varphi - \varphi_0)}, \qquad c = \frac{\ell^2}{\mu k}. \tag{3.12}
$$
The two constants of the second-order equation are $e$ and $\varphi_0$; the third constant of the orbit is $\ell$, hidden in $c$.
3. *Why it is a conic.* Rotate axes so $\varphi_0 = 0$ and put $x = r\cos\varphi$. Then (3.12) reads $r + ex = c$, i.e. $r = e\,(c/e - x)$: the distance to the origin equals $e$ times the distance to the vertical line $x = c/e$. That is the focus–directrix definition of a conic of eccentricity $e$: circle ($e = 0$), ellipse ($0 < e < 1$), parabola ($e = 1$), hyperbola ($e > 1$, one branch, the one with $1 + e\cos\varphi > 0$).
4. *The ellipse explicitly.* Square $r = c - ex$: $x^2 + y^2 = c^2 - 2cex + e^2x^2$, so $(1 - e^2)x^2 + 2cex + y^2 = c^2$. For $e < 1$ complete the square:
$$
(1 - e^2)\Big(x + \frac{ce}{1 - e^2}\Big)^2 + y^2 = \frac{c^2}{1 - e^2} \;\Longrightarrow\; \frac{(x + ae)^2}{a^2} + \frac{y^2}{b^2} = 1,
$$
$$
a = \frac{c}{1 - e^2}, \qquad b = a\sqrt{1 - e^2}, \qquad r_p = \frac{c}{1+e} = a(1 - e), \qquad r_a = \frac{c}{1-e} = a(1 + e). \tag{3.13}
$$
The center is at $x = -ae$, so the origin (the force center) sits a distance $ae$ from the center: at a focus. ∎

■ **What just happened.** Binet's equation became a driven oscillator with natural frequency exactly 1 in $\varphi$, so $r$ repeats after $\Delta\varphi = 2\pi$: §2's "$\omega_r = \omega_\varphi$," now exact at every amplitude. The ellipse geometry is (3.13); $c = a(1 - e^2)$ is the radius at $\varphi = 90°$.

### 3.7 Energy, eccentricity, and Kepler's third law

**Theorem 3.6 (Energy and eccentricity).** For the Kepler orbit (3.12),

$$
e = \sqrt{1 + \frac{2E\ell^2}{\mu k^2}}, \tag{3.14}
$$

and for $E \ne 0$

$$
a = \frac{c}{1 - e^2} = -\frac{k}{2E}. \tag{3.15}
$$

▶ **Try it first.** Evaluate $E = V_{\rm eff}$ at perihelion, where $\dot r = 0$, using $1/r_p = (1 + e)/c$.

**Derivation.**

1. At perihelion $\dot r = 0$, so by (3.8) $E = -k/r_p + \ell^2/2\mu r_p^2$.
2. Put $1/r_p = (1 + e)\mu k/\ell^2$:
$$
E = -\frac{\mu k^2}{\ell^2}(1 + e) + \frac{\ell^2}{2\mu}\cdot\frac{\mu^2k^2}{\ell^4}(1 + e)^2 = \frac{\mu k^2}{2\ell^2}\big[(1 + e)^2 - 2(1 + e)\big] = \frac{\mu k^2}{2\ell^2}(e^2 - 1).
$$
3. Solve for $e$: $e^2 = 1 + 2E\ell^2/\mu k^2$, which is (3.14).
4. Then $1 - e^2 = -2E\ell^2/\mu k^2$ and $a = c/(1 - e^2) = (\ell^2/\mu k)\cdot(-\mu k^2/2E\ell^2) = -k/2E$. ∎

Signs: $E < 0 \Leftrightarrow e < 1 \Leftrightarrow a > 0$ (ellipse); $E = 0 \Leftrightarrow e = 1$ (parabola, $a = \infty$); $E > 0 \Leftrightarrow e > 1 \Leftrightarrow a < 0$ (hyperbola, $|a| = k/2E$). The minimum energy (3.9) gives $e = 0$. Since $k$ enters (3.14) squared, it also holds for a repulsive force (§3.11).

■ **What just happened.** Why (3.15) is remarkable: in a general potential a bound orbit has two independent sizes, $r_p$ and $r_a$, and $E$ depends on both. For $-k/r$ it depends only on their sum, so every orbit with the same $a$, from the circle to a needle-thin ellipse, has the same energy, and (next theorem) the same period. This degeneracy is the conserved Laplace–Runge–Lenz vector of CM-4 §3.8 at work; in QM-6 it becomes hydrogen's energy depending on $n$ but not on $\ell$.

**Theorem 3.7 (Kepler's third law).** A bound Kepler orbit has period

$$
T^2 = \frac{4\pi^2\mu}{k}\,a^3 = \frac{4\pi^2}{G(m_1 + m_2)}\,a^3. \tag{3.16}
$$

**Derivation.**

1. By (3.6) the area grows at the constant rate $\ell/2\mu$, and in one period it covers the ellipse, area $\pi ab$. So $T = \pi ab\big/(\ell/2\mu) = 2\pi\mu ab/\ell$.
2. From $c = \ell^2/\mu k$ and $c = a(1 - e^2)$: $\ell = \sqrt{\mu k a(1 - e^2)}$. With $b = a\sqrt{1 - e^2}$:
$$
T = \frac{2\pi\mu\,a^2\sqrt{1 - e^2}}{\sqrt{\mu k a}\,\sqrt{1 - e^2}} = 2\pi\sqrt{\frac{\mu a^3}{k}}.
$$
3. For gravity $\mu/k = \dfrac{m_1m_2/M}{Gm_1m_2} = \dfrac{1}{G(m_1 + m_2)}$. ∎

■ **What just happened.** The eccentricity cancelled, which CM-4's scaling argument could not show. Kepler's constant $T^2/a^3 = 4\pi^2/G(M_\odot + m)$ differs between planets by their mass fractions, about $10^{-3}$ for Jupiter (Exercise A4), and $a$ is the semi-major axis of the *relative* orbit. Project Part 2 measures ten orbits with random $a, e$: log–log slope $1.50000000$.

### 3.8 Speeds, vis-viva, and orbit transfers

At perihelion and aphelion the velocity is perpendicular to $\mathbf r$, so $\ell = \mu r_pv_p = \mu r_av_a$. With $\ell^2 = \mu k a(1 - e^2)$ and (3.13):

$$
v_p = \sqrt{\frac{k}{\mu a}}\sqrt{\frac{1 + e}{1 - e}}, \qquad v_a = \sqrt{\frac{k}{\mu a}}\sqrt{\frac{1 - e}{1 + e}}, \qquad \frac{v_p}{v_a} = \frac{1 + e}{1 - e}. \tag{3.17}
$$

(For $v_p$: $v_p^2 = \ell^2/\mu^2r_p^2 = \mu ka(1-e^2)/\mu^2a^2(1-e)^2 = (k/\mu a)(1+e)/(1-e)$.)

**Theorem 3.8 (Vis-viva).** On any Kepler orbit, at distance $r$ from the force center,

$$
v^2 = \frac{k}{\mu}\left(\frac{2}{r} - \frac{1}{a}\right), \tag{3.18}
$$

with $a = -k/2E$ (negative for hyperbolas, infinite for parabolas).

*Proof.* $E = \tfrac12\mu v^2 - k/r$ and $E = -k/2a$ by (3.15). Solve for $v^2$. ∎

It is the most-used formula in astrodynamics: speed depends only on where you are and on the orbit's size. $r = a$ on a circle gives $v_c = \sqrt{k/\mu r}$; $a \to \infty$ gives escape speed $\sqrt{2k/\mu r} = \sqrt2\,v_c$.

**Hohmann transfer.** To go from a circular orbit $r_1$ to a coplanar circular orbit $r_2 > r_1$ with two short burns, fly half of the ellipse touching both: $r_p = r_1$, $r_a = r_2$, $a_t = (r_1 + r_2)/2$. For a spacecraft $k/\mu \to GM$. By (3.18):

$$
\Delta v_1 = \sqrt{\frac{GM}{r_1}}\left(\sqrt{\frac{2r_2}{r_1 + r_2}} - 1\right), \qquad
\Delta v_2 = \sqrt{\frac{GM}{r_2}}\left(1 - \sqrt{\frac{2r_1}{r_1 + r_2}}\right), \qquad
t_{\rm transfer} = \pi\sqrt{\frac{a_t^3}{GM}}. \tag{3.19}
$$

(Burn 1 is the transfer ellipse's $v_p = \sqrt{GM(2/r_1 - 1/a_t)} = \sqrt{GM/r_1}\sqrt{2r_2/(r_1+r_2)}$ minus the circular speed at $r_1$; burn 2 is the circular speed at $r_2$ minus the ellipse's $v_a$.) Walter Hohmann published this in 1925. It is the cheapest two-burn transfer between coplanar circles; for $r_2/r_1 > 11.94$ a three-burn bi-elliptic transfer beats it.

### 3.9 Stability of circular orbits and closed orbits

For a general attractive force write $f(r) = -F(r) > 0$. A circular orbit at $r_0$ needs $V_{\rm eff}'(r_0) = 0$, i.e. $\ell^2 = \mu r_0^3f(r_0)$. Its stability and its small-oscillation frequency come from $V_{\rm eff}''$:

$$
\omega_r^2 = \frac{V_{\rm eff}''(r_0)}{\mu}, \qquad \frac{\omega_r^2}{\omega_\varphi^2} = 3 + \frac{r_0f'(r_0)}{f(r_0)}, \tag{3.20}
$$

where $\omega_\varphi = \ell/\mu r_0^2$. (Proof: $V_{\rm eff}'' = f' + 3\ell^2/\mu r_0^4 = f' + 3f/r_0$ and $\omega_\varphi^2 = \ell^2/\mu^2r_0^4 = f/\mu r_0$; divide.)

**Theorem 3.9 (Power-law forces).** For $F = -kr^n$ with $k > 0$, circular orbits exist at every radius, they are stable if and only if $n > -3$, and small radial oscillations about them have

$$
\omega_r^2 = \frac{(n + 3)\,k\,r_0^{n-1}}{\mu}, \qquad \frac{\omega_r}{\omega_\varphi} = \sqrt{n + 3}. \tag{3.21}
$$

▶ **Try it first.** Write $r = r_0 + x$ in (3.7), expand to first order in $x$, and use the circular-orbit condition to kill the zeroth-order term.

**Derivation.**

1. Circular orbit: $\ddot r = 0$ in (3.7) gives $\ell^2/\mu r_0^3 = kr_0^n$, i.e. $\ell^2 = \mu k r_0^{n+3}$. Any $r_0$ works by choosing $\ell$.
2. Put $r = r_0 + x$ with $\ell$ fixed. Expand $(r_0 + x)^{-3} \approx r_0^{-3}(1 - 3x/r_0)$ and $(r_0 + x)^n \approx r_0^n(1 + nx/r_0)$:
$$
\mu\ddot x = \frac{\ell^2}{\mu r_0^3}\Big(1 - \frac{3x}{r_0}\Big) - kr_0^n\Big(1 + \frac{nx}{r_0}\Big).
$$
3. The constant terms cancel by step 1. Replace $\ell^2/\mu r_0^3$ by $kr_0^n$ in the $x$ term: $\mu\ddot x = -(3 + n)\,k r_0^{n-1}\,x$.
4. If $n + 3 > 0$ this is a harmonic oscillator with $\omega_r^2 = (n+3)kr_0^{n-1}/\mu$; if $n + 3 < 0$ it is exponential growth: unstable.
5. $\omega_\varphi^2 = \ell^2/\mu^2r_0^4 = kr_0^{n-1}/\mu$ by step 1. Divide. ∎

■ **What just happened.** The centrifugal barrier, as a force, falls as $r^{-3}$. An attraction that falls faster loses to it when you move out and beats it when you move in: no restoring force. At $n = -3$ both have the same $r$-dependence, so whichever is stronger wins everywhere and orbits spiral in or out (Cotes's spirals). Gravity is safely stable, with $\omega_r = \omega_\varphi$ (Example 4.5).

**Closed orbits.** A near-circular orbit closes if $\omega_r/\omega_\varphi$ is rational. Perihelion to perihelion is $2\pi/\sqrt{n+3}$; perihelion to aphelion (the **apsidal angle**) is $\pi/\sqrt{n+3}$. $n = -2$: ratio 1, apsidal angle $180°$, ellipse with the center of force at a focus. $n = 1$ (isotropic oscillator): ratio 2, apsidal angle $90°$, ellipse with the center of force at its *center* (Exercise B5). $n = -1$: ratio $\sqrt2$, a rosette that never closes. Project Part 4 measures all of these.

**Theorem 3.10 (Bertrand, 1873; stated).** The only central forces for which *every* bounded orbit is closed are $F = -k/r^2$ and $F = -kr$.

*What is skipped.* Closure of near-circular orbits at every $r_0$ forces $3 + rf'/f$ in (3.20) to be a fixed rational square, so $f$ is a power law; expanding Binet's equation to third order in amplitude then allows only $n + 3 = 1$ or $4$. Full proof: Goldstein, Poole & Safko, Appendix A. Physically, closed planetary orbits are evidence for the inverse-square law, and *small departures from closure measure small departures from it*.

### 3.10 Precession under small perturbations

**Theorem 3.11 (Inverse-cube perturbation).** For $F = -k/r^2 - \lambda/r^3$ with $\mu\lambda < \ell^2$, the bound orbits are precessing ellipses, $u = \dfrac{\mu k}{\ell^2\beta^2}\big[1 + e\cos\beta\varphi\big]$ with $\beta^2 = 1 - \mu\lambda/\ell^2$, and the perihelion advances by

$$
\Delta = 2\pi\Big(\frac1\beta - 1\Big) \approx \frac{\pi\mu\lambda}{\ell^2} \quad\text{per orbit.} \tag{3.22}
$$

**Derivation.**

1. Right side of (3.10): $-\dfrac{\mu}{\ell^2u^2}(-ku^2 - \lambda u^3) = \dfrac{\mu k}{\ell^2} + \dfrac{\mu\lambda}{\ell^2}u$.
2. Move the $u$ term left: $u'' + \beta^2u = \mu k/\ell^2$. Still linear, constant coefficients, now with natural frequency $\beta$.
3. Solution: $u = \mu k/\ell^2\beta^2 + A\cos\beta\varphi$, the stated form.
4. Perihelia (maxima of $u$) occur at $\beta\varphi = 2\pi m$, i.e. $\varphi = 2\pi m/\beta$. Successive perihelia are $2\pi/\beta$ apart, so the advance per orbit is $2\pi/\beta - 2\pi$.
5. For small $x = \mu\lambda/\ell^2$, $1/\beta = (1 - x)^{-1/2} \approx 1 + x/2$, so $\Delta \approx \pi x$. ∎

■ **What just happened.** A small extra attraction ($\lambda > 0$) makes the radial oscillation slightly *slower* than the orbital motion, so perihelion arrives late each time and creeps forward. The near-circular formula (3.20) agrees to first order: with $f = k/r^2 + \lambda/r^3$, $3 + rf'/f \approx 1 - \lambda/kr_0$, so $\omega_\varphi/\omega_r \approx 1 + \lambda/2kr_0$ and the advance is $\pi\lambda/kr_0 = \pi\mu\lambda/\ell^2$. Project Part 5 measures $0.076163$ rad per orbit for $\lambda = 0.02$, $\ell^2 = 0.84$ (exact formula to $10^{-6}$); first order gives $0.074800$, $1.8\%$ low.

**Mercury and general relativity.** In 1859 Urbain Le Verrier found that Mercury's perihelion advanced faster than the other planets' pull explains; Simon Newcomb later put the excess at about $43''$ per century, on top of roughly $532''$ from the planets. General relativity (Einstein, November 1915) adds $(3GM/c^2)u^2$ to the right side of Binet's equation (written per unit mass, with $\ell/\mu$ the specific angular momentum). As a small perturbation (Exercise C1) it gives

$$
\delta\varphi = \frac{6\pi GM}{c^2a(1 - e^2)} \quad\text{per orbit.} \tag{3.23}
$$

We state the modified equation, which comes from the Schwarzschild geodesics in any GR text (Hartle's *Gravity*, Carroll's *Spacetime and Geometry*). For Mercury ($a = 5.7909\times10^{10}$ m, $e = 0.20563$, period $87.969$ d) the project gets $5.019\times10^{-7}$ rad per orbit, $42.98''$ per century.

### 3.11 Scattering and the cross-section

You cannot follow one alpha particle through a gold foil; you can count how many come out at each angle. The quantity connecting that count to the force law is the cross-section, a new kind of quantity: an *area per target* that turns a beam's flux into a counting rate.

**Setup.** A uniform beam of particles of energy $E$ travels along $+\hat{\mathbf x}$ toward a fixed force center with **intensity** $I$ (particles per unit area per unit time). A particle's **impact parameter** $b$ is the distance from the force center to its incoming straight-line path. $\ell = \mu v_\infty b$ and $E = \tfrac12\mu v_\infty^2$ fix the orbit, so the deflection is a function $\theta(b)$. The **differential cross-section** is defined by

$$
(\text{particles per unit time scattered into solid angle } d\Omega \text{ around } \theta) = I\,\frac{d\sigma}{d\Omega}\,d\Omega.
$$

**Theorem 3.12 (Cross-section from the deflection function).** If $\theta(b)$ is monotonic,

$$
\frac{d\sigma}{d\Omega} = \frac{b}{\sin\theta}\left|\frac{db}{d\theta}\right|. \tag{3.24}
$$

▶ **Try it first.** Count particles in an annulus of the beam, $b$ to $b + db$. Where do they all go?

**Derivation.**

1. Particles with impact parameter between $b$ and $b + db$ cross the annulus of area $2\pi b\,|db|$ per unit time at rate $I\cdot 2\pi b\,|db|$.
2. By symmetry about the beam axis, all of them leave at angles between $\theta(b)$ and $\theta(b) + d\theta$, into the cone shell of solid angle $d\Omega = 2\pi\sin\theta\,|d\theta|$.
3. No particles are created or lost: $I\cdot 2\pi b\,|db| = I\,\dfrac{d\sigma}{d\Omega}\,2\pi\sin\theta\,|d\theta|$. Cancel $2\pi I$. ∎

■ **What just happened.** $d\sigma/d\Omega$ is the Jacobian of the map from beam area $2\pi b\,db$ to solid angle $2\pi\sin\theta\,d\theta$. The **total cross-section** $\sigma = \int(d\sigma/d\Omega)\,d\Omega$ is the area of beam scattered at all. For a thin foil with $n$ targets per unit volume and thickness $t$, the probability of scattering into some range of angles is $n\,t\,\sigma_{\rm range}$, valid while that product is small (single scattering).

**Hard sphere.** A particle bouncing elastically off a fixed sphere of radius $R$ hits at angle of incidence $\alpha$ to the normal, with $b = R\sin\alpha$, and reflects symmetrically, so $\theta = \pi - 2\alpha$. Hence $b = R\cos(\theta/2)$, $|db/d\theta| = \tfrac12R\sin(\theta/2)$, and with $\sin\theta = 2\sin(\theta/2)\cos(\theta/2)$:

$$
\frac{d\sigma}{d\Omega} = \frac{R\cos(\theta/2)}{2\sin(\theta/2)\cos(\theta/2)}\cdot\frac{R\sin(\theta/2)}{2} = \frac{R^2}{4}, \qquad \sigma = 4\pi\cdot\frac{R^2}{4} = \pi R^2. \tag{3.25}
$$

Isotropic, and the total is the geometric shadow: half the scattered particles go backward.

**The repulsive Coulomb orbit.** Take $V = +k/r$, so $F = +k/r^2$ and $F(1/u) = ku^2$.

**Theorem 3.13 (Coulomb deflection).** A particle of energy $E > 0$ and impact parameter $b$ in the potential $+k/r$ moves on a hyperbola of eccentricity $e = \sqrt{1 + (2Eb/k)^2}$ and is deflected by

$$
\sin\frac\theta2 = \frac1e, \qquad\text{equivalently}\qquad b = \frac{k}{2E}\cot\frac\theta2. \tag{3.26}
$$

▶ **Try it first.** Solve Binet's equation with the sign of $k$ flipped. Then find the two directions in which $r \to \infty$, and convince yourself that the *velocity* at infinity points along the *position* direction at infinity.

**Derivation.**

1. Binet (3.10): $u'' + u = -\dfrac{\mu}{\ell^2u^2}ku^2 = -\dfrac{\mu k}{\ell^2}$. General solution $u = -\mu k/\ell^2 + A\cos\varphi$, with $\varphi = 0$ chosen at closest approach (maximum $u$), $A > 0$. Write $A = e\mu k/\ell^2$:
$$
u = \frac{\mu k}{\ell^2}\,(e\cos\varphi - 1), \qquad r = \frac{c}{e\cos\varphi - 1}, \quad c = \frac{\ell^2}{\mu k}.
$$
2. Since $u = 1/r > 0$, we need $e\cos\varphi > 1$. So $e > 1$ (a hyperbola) and the polar angle is confined to $-\varphi_\infty < \varphi < \varphi_\infty$ with $\cos\varphi_\infty = 1/e$. As $\varphi \to \pm\varphi_\infty$, $r \to \infty$. (This is the branch of the hyperbola that curves *away* from the focus.)
3. Energy. At closest approach $r_{\min} = c/(e - 1)$ and $\dot r = 0$, so $E = \ell^2/2\mu r_{\min}^2 + k/r_{\min}$. With $1/r_{\min} = (e-1)\mu k/\ell^2$: $E = \dfrac{\mu k^2}{2\ell^2}\big[(e-1)^2 + 2(e - 1)\big] = \dfrac{\mu k^2}{2\ell^2}(e^2 - 1)$. Same as (3.14).
4. Asymptotic directions, carefully. Far from the center, the transverse velocity is $r\dot\varphi = \ell/\mu r \to 0$, so the velocity becomes purely radial: $\dot{\mathbf r} \to \dot r\,\hat{\mathbf r}$. Coming in ($t \to -\infty$), the particle is at polar angle $-\varphi_\infty$ (say) with $\dot r < 0$, so its velocity points along $-\hat{\mathbf r}(-\varphi_\infty)$, the direction with polar angle $\pi - \varphi_\infty$. Going out, it is at $+\varphi_\infty$ with $\dot r > 0$, velocity along $\hat{\mathbf r}(\varphi_\infty)$, polar angle $\varphi_\infty$. The asymptotes do not pass through the origin (they are offset by $b$), but that offset does not matter for directions: $b/r \to 0$.
5. Deflection. The angle between the incoming direction $\pi - \varphi_\infty$ and outgoing direction $\varphi_\infty$ is $\theta = (\pi - \varphi_\infty) - \varphi_\infty = \pi - 2\varphi_\infty$. So $\theta/2 = \pi/2 - \varphi_\infty$ and $\sin(\theta/2) = \cos\varphi_\infty = 1/e$.
6. Impact parameter. At infinity $\ell = \mu v_\infty b$ and $E = \tfrac12\mu v_\infty^2$, so $\dfrac{2E\ell^2}{\mu k^2} = \dfrac{\mu v_\infty^2\cdot\mu^2v_\infty^2b^2}{\mu k^2} = \Big(\dfrac{\mu v_\infty^2b}{k}\Big)^2 = \Big(\dfrac{2Eb}{k}\Big)^2$, giving $e^2 = 1 + (2Eb/k)^2$. Then $\cot^2(\theta/2) = 1/\sin^2(\theta/2) - 1 = e^2 - 1 = (2Eb/k)^2$. ∎

■ **What just happened.** Most books skip step 4 with "by symmetry." The deflection is an angle between *velocities*; the orbit equation gives *positions*; they agree at infinity only because fixed $\ell = \mu r\cdot(r\dot\varphi)$ forces the transverse speed to zero. For an *attractive* hyperbola the same steps give $r = c/(1 + e\cos\varphi)$, $\cos\varphi_\infty = -1/e$, $\theta = 2\varphi_\infty - \pi$, and again $\sin(\theta/2) = 1/e$ (Example 4.3).

**Theorem 3.14 (Rutherford's formula).** For the Coulomb potential $\pm k/r$,

$$
\frac{d\sigma}{d\Omega} = \left(\frac{k}{4E}\right)^2\frac{1}{\sin^4(\theta/2)}. \tag{3.27}
$$

**Derivation.**

1. From (3.26), $b = (k/2E)\cot(\theta/2)$, so $\dfrac{db}{d\theta} = -\dfrac{k}{2E}\cdot\dfrac{1}{2\sin^2(\theta/2)}$.
2. Into (3.24), with $\sin\theta = 2\sin(\theta/2)\cos(\theta/2)$:
$$
\frac{d\sigma}{d\Omega} = \frac{(k/2E)\cos(\theta/2)/\sin(\theta/2)}{2\sin(\theta/2)\cos(\theta/2)}\cdot\frac{k}{4E\sin^2(\theta/2)} = \frac{k}{2E}\cdot\frac{k}{4E}\cdot\frac{1}{2\sin^4(\theta/2)} = \Big(\frac{k}{4E}\Big)^2\frac{1}{\sin^4(\theta/2)}.
$$
∎

Two consequences. Everything with $b < b(\theta_0)$ scatters beyond $\theta_0$ (smaller $b$, bigger deflection), so

$$
\sigma(\theta > \theta_0) = \pi\,b(\theta_0)^2 = \pi\Big(\frac{k}{2E}\Big)^2\cot^2\frac{\theta_0}{2}; \qquad \sigma(\theta > 90°) = \pi\Big(\frac{k}{2E}\Big)^2. \tag{3.28}
$$

As $\theta_0 \to 0$ this diverges: the Coulomb force has infinite range. In a real atom the electrons screen the nucleus beyond about $10^{-10}$ m and cut off the small angles.

■ **What just happened, and why it meant a nucleus.** The hard sphere scatters isotropically; the point charge scatters overwhelmingly forward, with a large-angle tail falling only as $\sin^{-4}(\theta/2)$ ($150°$ is three thousand times weaker than $15°$, not zero). Large angles need small $b$, hence close approaches, down to $r_{\min} = k/E$ head-on. Charge spread over the whole atom (Thomson) would deflect an alpha particle by only about $k/(ER_{\rm atom})$ per atom, a fraction of a milliradian (Example 4.6). So the $\sin^{-4}$ tail out to $150°$ says the charge sits inside the alpha particles' closest approach, tens of femtometres: a nucleus (QM-1). Two cautions: $\theta$ and $E$ in (3.27) are center-of-mass quantities ($E_{\rm cm} = 0.980\,E_{\rm lab}$ for alpha on gold), and the classical formula agrees exactly with the quantum one for pure Coulomb, a coincidence QM explains.

### 3.12 The Laplace–Runge–Lenz vector points at perihelion

CM-4 Example 4.6 showed that for $V = -k/r$

$$
\mathbf A = \mathbf p\times\mathbf L - \mu k\,\hat{\mathbf r}, \qquad \mathbf p = \mu\dot{\mathbf r}, \tag{3.29}
$$

is constant. Here is what it is.

**Theorem 3.15.** $\mathbf A$ points from the force center to the perihelion, and $|\mathbf A| = \mu k e$.

**Derivation.**

1. $\mathbf A$ is constant, so evaluate it at perihelion, where $\mathbf p \perp \mathbf r$. With $\hat{\mathbf z}$ along $\mathbf L$, $\mathbf p = p\,\hat{\boldsymbol\varphi}$ and $\hat{\boldsymbol\varphi}\times\hat{\mathbf z} = \hat{\mathbf r}$, so $\mathbf p\times\mathbf L = p\ell\,\hat{\mathbf r}$.
2. At perihelion $\ell = r_pp$, so $p\ell = \ell^2/r_p = \ell^2(1 + e)\mu k/\ell^2 = \mu k(1 + e)$.
3. $\mathbf A = [\mu k(1 + e) - \mu k]\,\hat{\mathbf r}_p = \mu k e\,\hat{\mathbf r}_p$. ∎

■ **What just happened.** A circle has $\mathbf A = \mathbf 0$; every other orbit carries a fixed arrow to its perihelion, and by (3.14) $|\mathbf A|^2 = \mu^2k^2 + 2\mu E\ell^2$. "The perihelion does not move" and "$\mathbf A$ is conserved" are the same statement, so a perturbation that makes the orbit precess (§3.10) makes $\mathbf A$ rotate (Exercise B8). Project Part 7 checks $|\mathbf A| = \mu ke$ to $10^{-8}$ for $e = 0$ to $1.5$.

### 3.13 Three bodies

Add a third body and, after the center of mass, energy and angular momentum are used, too many degrees of freedom remain; Poincaré showed (around 1890) that there are no further conserved quantities and the general motion is chaotic, so no closed form of this module's kind exists. Special solutions survive, notably Lagrange's (1772) equilateral L4 and L5 points, where Jupiter's Trojan asteroids sit (Exercise C4, project "Your turn").

## 4. Worked examples

### Example 4.1: Earth and Sun

*Reduced mass.* Gravitational parameters are measured far more precisely than masses: $GM_\odot = 1.32712440018\times10^{20}$ m³/s², $Gm_\oplus = 3.986004418\times10^{14}$ m³/s², ratio $3.003\times10^{-6}$. So
$$
\frac{\mu}{m_\oplus} = \frac{M_\odot}{M_\odot + m_\oplus} = 0.999996997 = 1 - 3.003\times10^{-6}.
$$
The Sun's distance from the barycenter is $a\,m_\oplus/(M_\odot + m_\oplus) = 449$ km, far inside the Sun.

*Orbit.* $a = 1$ AU $= 1.495978707\times10^{11}$ m, $e = 0.0167$. By (3.13), $r_p = a(1 - e) = 147.10$ Gm and $r_a = a(1 + e) = 152.10$ Gm.

*Speeds.* $k/\mu = G(M_\odot + m_\oplus)$. At $r = a$ (ends of the minor axis) $v = \sqrt{G(M_\odot + m_\oplus)/a} = 29.785$ km/s. By (3.17), multiply or divide by $\sqrt{(1 + e)/(1 - e)} = 1.01684$: $v_p = 30.286$ km/s, $v_a = 29.291$ km/s, and $v_p/v_a = 1.0340 = r_a/r_p$.

*Period.* By (3.16) with $G(M_\odot + m_\oplus)$: $T = 365.25635$ d; with $GM_\odot$ alone, $365.25690$ d. Earth's mass shortens the year by $47.4$ s. The measured sidereal year is $365.25636$ d; agreement this close is partly luck (the Moon and planets shift Earth's mean motion at the $10^{-6}$ level too), but the correction is real. And if you compute $GM_\odot$ as $G\times M_\odot$ ($6.67430\times10^{-11}\times1.9885\times10^{30}$), you get $365.24862$ d, wrong by $0.0083$ d because $G$ is known only to about $2\times10^{-5}$: 15 times the effect you were computing.

### Example 4.2: An orbit with its constants fixed by $E$ and $\ell$

In Kepler units ($\mu = k = 1$), a particle at $r = 1$, $\varphi = 0$ has transverse velocity $r\dot\varphi = 1$ and is moving outward with $E = -0.3$. Find the orbit.

*Constants.* $\ell = \mu r\cdot r\dot\varphi = 1$. From $E = \tfrac12\dot r^2 + \tfrac12(r\dot\varphi)^2 - 1/r$: $\dot r^2 = 2(-0.3 + 1) - 1 = 0.4$, so $\dot r = 0.6325$. By (3.12) and (3.14): $c = \ell^2/\mu k = 1$ and $e = \sqrt{1 + 2(-0.3)(1)} = \sqrt{0.4} = 0.6325$.

*Orientation.* $r(0) = 1 = c/(1 + e\cos\varphi_0)$ forces $\cos\varphi_0 = 0$. Moving outward means perihelion is behind us, so $\varphi_0 = -90°$:
$$
r(\varphi) = \frac{1}{1 + 0.6325\cos(\varphi + 90°)} = \frac{1}{1 - 0.6325\sin\varphi}.
$$
(Check $\dot r$: $\dot r = -(\ell/\mu)u' = -(-0.6325\cos\varphi) = 0.6325$ at $\varphi = 0$. ✓)

*Geometry and motion.* $r_p = 1/1.6325 = 0.6126$, $r_a = 1/0.3675 = 2.7208$, $a = 1/(1 - 0.4) = 1.6667 = -k/2E$ ✓, $b = a\sqrt{0.6} = 1.2910$, $T = 2\pi(1.6667)^{3/2} = 13.519$, $v_p = \ell/r_p = 1.6325$, $v_a = 0.3675$. Project Part 3 integrates from these initial conditions, fits a conic, and recovers $c = 1.000000$, $e = 0.632456$, perihelion at $-90.0000°$.

### Example 4.3: A hyperbolic visitor

An object arrives from interstellar space at $v_\infty = 26.4$ km/s and passes the Sun at perihelion distance $q = 0.255$ AU (numbers close to those of 1I/ʻOumuamua, 2017). Find its orbit and how much the Sun bends its path.

*Energy and size.* Per unit reduced mass, $E/\mu = \tfrac12v_\infty^2 = 348.5$ MJ/kg $> 0$: unbound. By (3.15), $|a| = k/2E = GM_\odot/v_\infty^2 = 1.904\times10^{11}$ m $= 1.273$ AU.

*Eccentricity.* For the attractive hyperbola $r_p = c/(1 + e)$ with $c = a(1 - e^2) = |a|(e^2 - 1)$, so $q = |a|(e - 1)$ and $e = 1 + q/|a| = 1 + 0.255/1.273 = 1.200$.

*Deflection.* By (3.26), $\sin(\theta/2) = 1/e = 0.833$, so $\theta = 112.8°$: more than a right angle.

*Impact parameter and perihelion speed.* $\cot(\theta/2) = 2Eb/k = b/|a|$ gives $b = |a|\sqrt{e^2 - 1} = 0.845$ AU: aimed $0.845$ AU from the Sun, it passes at $0.255$ AU. Vis-viva (3.18) with $a = -|a|$: $v_p = \sqrt{GM_\odot(2/q + 1/|a|)} = 87.5$ km/s. Check: $v_\infty b = v_pq = 3.34\times10^{15}$ m²/s (Part 9).

### Example 4.4: Hohmann transfer from low Earth orbit to geostationary orbit

$GM_\oplus = 3.986004418\times10^{14}$ m³/s², $r_1 = 6771$ km (400 km altitude), $r_2 = 42{,}164$ km (geostationary).

*Circular speeds.* $v_{c1} = \sqrt{GM/r_1} = 7672.6$ m/s, $v_{c2} = \sqrt{GM/r_2} = 3074.7$ m/s.

*Transfer ellipse.* $a_t = (r_1 + r_2)/2 = 24{,}467.5$ km. Vis-viva: $v_p = \sqrt{GM(2/r_1 - 1/a_t)} = 10{,}072.1$ m/s at perigee, $v_a = \sqrt{GM(2/r_2 - 1/a_t)} = 1617.4$ m/s at apogee. Check: $v_pr_1 = v_ar_2$ ✓.

*Burns.* $\Delta v_1 = 10{,}072.1 - 7672.6 = 2399.5$ m/s (prograde, at perigee). Coast half an ellipse: $t = \pi\sqrt{a_t^3/GM} = 5.290$ h. $\Delta v_2 = 3074.7 - 1617.4 = 1457.2$ m/s (prograde, at apogee). Total $3856.7$ m/s.

Both burns speed the craft up, yet it ends *slower* than it started ($3.07$ vs $7.67$ km/s): raising an orbit costs energy but, by vis-viva, lowers the speed. Project Part 8 flies both burns and finds the final orbit circular at $r_2$ to $10^{-7}$.

### Example 4.5: Radial oscillations about a circular orbit

For $n = -2$, (3.21) gives $\omega_r = \omega_\varphi = \sqrt{k/\mu r_0^3}$. Put numbers on it: a satellite on a circular orbit at $r_0 = 6771$ km gets a radial kick $\Delta v_r = 10$ m/s.

*Frequencies.* $v_c = 7672.6$ m/s, $\omega_r = \omega_\varphi = v_c/r_0 = 1.1332\times10^{-3}$ rad/s: one radial oscillation per $92.4$ min orbit.

*Linear theory.* After the kick $x(t) = (\Delta v_r/\omega_r)\sin\omega_rt$, amplitude $8.82$ km. Since $\omega_r = \omega_\varphi$, low and high points stay $180°$ apart: a slightly eccentric ellipse, $e \approx 8.82/6771$.

*Exact check.* A radial kick leaves $\ell = \mu r_0v_c$ unchanged, so $c = \ell^2/\mu k = r_0$, and adds $\tfrac12\mu\Delta v_r^2$ to $E = -\tfrac12\mu v_c^2$. By (3.14): $e^2 = 1 + 2(E/\mu)(r_0v_c)^2/(GM)^2 = 1 + (\Delta v_r^2 - v_c^2)/v_c^2 = (\Delta v_r/v_c)^2$, so $e = \Delta v_r/v_c = 1.3033\times10^{-3}$ exactly. Then $r_a - r_0 = r_0e/(1 - e) = 8.836$ km and $r_0 - r_p = r_0e/(1 + e) = 8.813$ km, bracketing the linear value.

*Other powers* (Part 4, $3\%$ amplitude). $n = -1$: measured $\omega_r/\omega_\varphi = 1.414319$ against $\sqrt2 = 1.414214$ (the gap is second order in amplitude), apsidal angle $127.27°$, never closes. $n = 1$: $2.000000$, apsidal angle $90°$. $n = -2.1$: $0.948627$ against $\sqrt{0.9} = 0.948683$, apsidal angle $189.75°$, slow forward precession.

### Example 4.6: Alpha particles on gold

$5.0$ MeV alpha particles ($Z_1 = 2$) hit a gold foil ($Z_2 = 79$, $A = 196.967$, density $19.32$ g/cm³) $0.4$ μm thick. Use $e^2/4\pi\varepsilon_0 = 1.43996$ MeV·fm.

*Strength.* $k = Z_1Z_2e^2/4\pi\varepsilon_0 = 227.5$ MeV fm; head-on $r_{\min} = k/E = 45.5$ fm, outside the gold nucleus (radius $\approx 1.2A^{1/3} = 7.0$ fm), so (3.27) should hold at all angles. We ignore the 2% recoil correction ($E_{\rm cm}/E_{\rm lab} = 0.9801$).

*Cross-section.* $k/4E = 11.38$ fm, $(k/4E)^2 = 129.4$ fm²/sr $= 1.29$ b/sr (1 barn $= 100$ fm²). At $90°$, $\sin^4 45° = 1/4$, so $d\sigma/d\Omega = 518$ fm²/sr $= 5.18$ b/sr. By (3.28), $\sigma(\theta > 90°) = \pi(k/2E)^2 = \pi(22.75\text{ fm})^2 = 1626$ fm² $= 16.3$ b.

*Fraction beyond $90°$.* $n = 19.32\times6.022\times10^{23}/196.967 = 5.907\times10^{22}$ atoms/cm³. The probability is $n\,t\,\sigma = 5.907\times10^{22}\times(0.4\times10^{-4}\text{ cm})\times(1.626\times10^{-23}\text{ cm}^2) = 3.84\times10^{-5}$: about one in 26,000 bounces back, the same order as Geiger and Marsden's 1 in 8000 from platinum (different target and thickness, so only the order compares).

*The shape.* $d\sigma/d\Omega$ falls by $1/2999$ from $15°$ to $150°$ and by $240{,}468$ from $5°$ to $150°$; Geiger and Marsden's 1913 counts tracked $\sin^{-4}(\theta/2)$ across that range. A $7.0$ fm hard sphere would scatter isotropically with total cross-section $153$ fm², a tenth of the Coulomb cross-section beyond $90°$.

*Why Thomson's atom fails.* Spread the charge over $R_{\rm atom} \approx 10^5$ fm. Crossing it, an alpha particle feels at most $\sim k/R^2$ for a time $\sim 2R/v$: deflection $\Delta p/p \sim 2k/(Rmv^2) = k/(ER) = 4.55\times10^{-4}$ rad. The foil is $t\,n^{1/3} = 1558$ atomic layers; random kicks add as a random walk to $\sqrt{1558}\times4.55\times10^{-4} = 0.0180$ rad ($1.03°$), and reaching $90°$ has probability of order $\exp[-(\pi/2)^2/0.018^2] \approx 10^{-3322}$. Observed: about $10^{-4}$. Only one large kick from a point charge does it.

Project Part 6 integrates 12,000 trajectories with $b$ uniform over a disk and recovers (3.27) with median relative error $2.7\%$.

## 5. Why it matters and where it connects

- **QM-6 (hydrogen).** The reduction (3.3) is the first line, with $\mu$ the electron–proton reduced mass, so every level is the infinite-mass value times $\mu/m_e = 0.99946$ (Exercise A1). The radial equation contains $V + \hbar^2\ell(\ell+1)/2\mu r^2$, the quantum (3.8). Energy depending only on $a$ (3.15) becomes energy depending only on $n$, for the same reason: the LRL vector (3.29).
- **QM-1.** The Rutherford atom of §3.11 is classically unstable (an orbiting electron radiates, EM-8's Larmor formula). QM-1 starts from that contradiction.
- **CM-4.** §3.12 is the payoff of CM-4 Example 4.6. CM-4's solution to Exercise B6(d) quotes $\langle r\rangle = a(1 + e^2/2)$, derived here in Exercise C2.
- **Statistical mechanics.** By CM-4's virial theorem a bound gravitating system has $E = -\langle T\rangle = -k/2a$ per orbit, so it *heats up* as it loses energy (negative heat capacity, CM-4 Exercise C2); self-gravitating systems therefore have no canonical equilibrium, a caveat for SM-4's canonical ensemble.
- **CM-8.** In action-angle variables, $\omega_r = \omega_\varphi$ is the degeneracy of the Kepler frequencies (CM-8 Exercise C2).
- **CM-6.** $\omega_r^2 = V_{\rm eff}''(r_0)/\mu$ (3.20) is a one-degree-of-freedom normal mode; CM-6 generalizes it to $K\mathbf a = \omega^2M\mathbf a$.
- **Engineering.** Vis-viva (3.18) and Hohmann (3.19) open every astrodynamics course; (3.20) is the stability test for orbits in any axisymmetric potential, from galaxies to storage rings.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills: reduced masses (hydrogen, deuterium, Earth–Moon, positronium), reading a conic, Mercury's speeds, the mass correction to Kepler III for Jupiter, and total cross-sections. B1–B8 are the module: B2 (find the force from the orbit with Binet's equation), B4 (stability in a Yukawa potential, where the golden ratio appears), B6 (closest approach and nuclear size) and B8 (extend `project.py` to watch the LRL vector precess) are essential. C1–C4 go further: the GR perihelion shift from the modified orbit equation, Kepler's equation and time averages, the cross-section of a $1/r^2$ potential, and the Lagrange points L4/L5.

## 7. Coding project

`project.py` integrates the two-body problem *without* the reduction and checks §3 against the raw numbers: (1) eight coupled ODEs for $m_1 = 1$, $m_2 = 0.3$, center of mass straight to $10^{-9}$, relative orbit the conic with $e = 0.230769$, $a = 0.8125$ from (3.14)–(3.15); (2) Kepler III over ten orbits, slope $1.50000000$; (3) $V_{\rm eff}$, turning points and the four conics; (4) power laws $n = -2, -1, 1, -2.1, -3.5$ against (3.21); (5) precession against (3.22) and Mercury from (3.23); (6) the Rutherford Monte Carlo against (3.27); (7) the LRL vector; (8) a simulated Hohmann transfer; (9) every worked-example number; (10) SymPy checks of Binet's chain rule, (3.11), (3.14), (3.16), (3.21), (3.22), (3.25), (3.27), (3.29), and the numbers in the solutions.

**Build it yourself:** before reading Part 1, write `two_body_rhs` and the least-squares conic fit `fit_conic` (fit $1/r$ as $A + B\cos\varphi + C\sin\varphi$; then $c = 1/A$, $e = \sqrt{B^2 + C^2}/A$). Before reading Part 6, write the Rutherford Monte Carlo: sample $b$ with density $\propto b$ (that is, $b = \sqrt{\text{uniform}}$ over the disk), integrate, histogram $\theta$, divide by $2\pi\sin\theta\,\Delta\theta$. "Your turn" asks for the restricted three-body problem at L4 and L1, orbital decay under drag, and Halley's comet.

## 8. Where you'll struggle

- **Skipping the reduction and getting the mass wrong in Kepler III.** Treat the Sun as fixed and write $T^2 = 4\pi^2a^3/GM_\odot$, and you are off by the planet's mass fraction: $3\times10^{-6}$ for Earth (47 s per year), $10^{-3}$ for Jupiter, order 1 for binary stars, where (3.16) with $m_1 + m_2$ is how stellar masses are measured. Fix: always write the problem in $\mathbf r = \mathbf r_1 - \mathbf r_2$ first and say out loud $k/\mu = G(m_1 + m_2)$. Related trap: $G\times M$ instead of the measured $GM$ (Example 4.1).
- **The $u = 1/r$ trick does not suggest itself.** Nobody invents it in an exam. Memorize the three lines of §3.5: $\dot\varphi = \ell u^2/\mu$, $\dot r = -(\ell/\mu)u'$, $\ddot r = -(\ell^2u^2/\mu^2)u''$. Then Binet's equation (3.10) is a one-line substitution. Practise it on Exercise B2, which runs it backward.
- **Signs of $E$ and $k$.** A bound orbit needs $E < 0$ in (3.14), and only then is $a = -k/2E$ positive. Students drop the minus sign and get a negative semi-major axis for Earth, or write $V = +k/r$ for gravity. For scattering we deliberately write $V = +k/r$ with $k > 0$: (3.14) and (3.27) survive because $k$ appears squared, but the orbit becomes $c/(e\cos\varphi - 1)$. Check every result against $e = 0$, the circle (3.9).
- **Substituting $\ell$ into the Lagrangian** flips the centrifugal sign (§3.3) and makes every circular orbit unstable. Substitute into the energy or the equations of motion only.
- **The asymptote step in Rutherford.** Most books write "by symmetry, $\theta = \pi - 2\varphi_\infty$." Redo step 4 of Theorem 3.13 until you can say *why* the velocity direction at infinity equals the position direction there ($r\dot\varphi = \ell/\mu r \to 0$), and draw it for both signs of $k$. Keep $\theta$ (deflection) and $\varphi$ (polar angle) apart; they are measured from different lines.
- **Cross-sections are a new kind of quantity.** $d\sigma/d\Omega$ is neither a probability nor a force: it is area per solid angle per target, a rate only after multiplying by beam intensity, a probability only after multiplying by $n\,t$. Write units on every cross-section (fm²/sr, barns) and name the multiplication each time. Exercise A5 (Rutherford's total cross-section is infinite) tests whether you know what it counts.
- **Rational ratios versus closure.** Closure needs $\omega_r/\omega_\varphi$ rational *at all amplitudes*; $\sqrt{n+3}$ is only first order (for $n = -1$ the measured ratio drifts off $\sqrt2$, Example 4.5). Bertrand's theorem says only two force laws pass at every amplitude.

## 9. Check yourself

1. Write $\mathbf r_1, \mathbf r_2$ in terms of $\mathbf R, \mathbf r$. Why does the cross term $\dot{\mathbf R}\cdot\dot{\mathbf r}$ cancel?
2. Prove in two lines that a central-force orbit lies in a plane. What happens when $\mathbf L = \mathbf 0$?
3. Write $V_{\rm eff}$ and the radial equation. What goes wrong if you substitute $\dot\varphi = \ell/\mu r^2$ into the Lagrangian?
4. Derive $\dot r = -(\ell/\mu)\,du/d\varphi$.
5. Write Binet's equation and solve it for $F = -k/r^2$. What is $c$?
6. Derive $e^2 = 1 + 2E\ell^2/\mu k^2$ from the energy at perihelion.
7. Why is it remarkable that $E = -k/2a$ for every eccentricity?
8. Derive Kepler's third law from the area law in three lines. Where does $m_1 + m_2$ enter?
9. State vis-viva. Use it to find the escape speed from radius $r$.
10. For $F = -kr^n$, when are circular orbits stable, and what is $\omega_r/\omega_\varphi$? Which two values of $n$ give closed orbits?
11. Derive $d\sigma/d\Omega = (b/\sin\theta)|db/d\theta|$. What is $d\sigma/d\Omega$ for a hard sphere?
12. Where does the Laplace–Runge–Lenz vector point, and how long is it?

## 10. Pair with

- **Taylor, *Classical Mechanics*, ch. 8.** The same route as these notes (reduced mass, $V_{\rm eff}$, $u = 1/r$, Kepler orbits, unbounded orbits, changes of orbit) with more pictures and more worked numbers; read the changes-of-orbit section alongside §3.8. **Taylor ch. 14** builds cross-sections up slowly and does the conversion between lab and center-of-mass angles properly, which these notes skip.
- **Goldstein, Poole & Safko, ch. 3.** Denser and more complete: the orbit integral by quadrature, Bertrand's theorem with the proof in Appendix A, the LRL vector (§3.9), and scattering in a central field (§3.10–3.11). Read it for Bertrand and for the attractive-scattering subtleties (orbiting and rainbow scattering) that these notes skip.
- **MIT 8.223 (OCW), lectures 7–10.** Lectures covering central-force motion and scattering; watching the effective-potential arguments drawn live is the best way to make §2 automatic.
