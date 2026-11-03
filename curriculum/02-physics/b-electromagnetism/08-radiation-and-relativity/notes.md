# EM-8: Potentials, Radiation, and Relativity: Retarded Potentials, Dipole Radiation, Lorentz Covariance, the Field Tensor

**Time:** ~10 hours · **Prerequisites:** EM-7, LA-12 · **Unlocks:** QM-1 (the stability-of-atoms problem), QM-8 (spontaneous emission context), SM-5 (blackbody context)

> **Conventions in this module.** SI units; $c = 1/\sqrt{\mu_0\varepsilon_0}$. Source point $\mathbf r'$, field point $\mathbf r$, separation $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$, $\mathscr r = |\boldsymbol{\mathscr r}|$, $\hat{\boldsymbol{\mathscr r}} = \boldsymbol{\mathscr r}/\mathscr r$ (EM-6's notation). **Retarded time** $t_r = t - \mathscr r/c$. Potentials $\mathbf E = -\nabla V - \partial_t\mathbf A$, $\mathbf B = \nabla\times\mathbf A$ in the **Lorenz gauge** $\nabla\cdot\mathbf A + c^{-2}\partial_tV = 0$; $\Box = \nabla^2 - c^{-2}\partial_t^2$ (EM-6 §3.12). A point charge $q$ moves on $\mathbf w(t)$ with velocity $\mathbf v$ and acceleration $\mathbf a$. $\beta = v/c$, $\gamma = 1/\sqrt{1 - \beta^2}$. **Relativity:** frame $\bar S$ moves at velocity $v\hat{\mathbf x}$ relative to $S$, axes parallel, origins coinciding at $t = \bar t = 0$. Events $x^\mu = (ct, x, y, z)$, $\mu = 0, 1, 2, 3$; Greek indices run 0–3, Latin 1–3; repeated upper–lower indices are summed (Einstein convention, LA-12 §3.15). **Metric signature $(-,+,+,+)$** as in Griffiths: $\eta_{\mu\nu} = \operatorname{diag}(-1, 1, 1, 1)$, $x_\mu = \eta_{\mu\nu}x^\nu = (-ct, x, y, z)$, $\partial_\mu = \partial/\partial x^\mu = (c^{-1}\partial_t, \nabla)$, $\partial^\mu = \eta^{\mu\nu}\partial_\nu$. **Field tensor (Griffiths' convention):** $F^{0i} = E_i/c$, $F^{ij} = \varepsilon_{ijk}B_k$; dual $G^{0i} = B_i$, $G^{ij} = -\varepsilon_{ijk}E_k/c$. With these, Maxwell's equations read $\partial_\nu F^{\mu\nu} = \mu_0J^\mu$ (derivative on the **second** index). Jackson and most QFT books use $(+,-,-,-)$ and $F^{0i} = -E_i$ (Gaussian or SI), which flips signs in several places (§8). Four-current $J^\mu = (c\rho, \mathbf J)$, four-potential $A^\mu = (V/c, \mathbf A)$.

## 0. Where this sits

EM-6 wrote Maxwell's equations as four wave equations for the potentials, $\Box V = -\rho/\varepsilon_0$ and $\Box\mathbf A = -\mu_0\mathbf J$ (EM-6 (3.37)), and promised their solution. EM-7 studied waves already in flight and teased the oscillating dipole (Rayleigh scattering, EM-7 §3.13). Here we solve (3.37), find out what makes waves in the first place, and then show that the theory you have been using for seven modules was relativistic all along: $\mathbf E$ and $\mathbf B$ are the six components of one antisymmetric tensor (LA-12 §3.17's transformation law, applied to a $4\times4$ matrix). You need LA-12's index notation, raising and lowering, and the rule $T'^{\mu\nu} = \Lambda^\mu{}_\alpha\Lambda^\nu{}_\beta T^{\alpha\beta}$. QM-1 opens with this module's Example 4.4 (the classical atom dies in $10^{-11}$ s); QM-8 needs the dipole power formula for spontaneous emission; SM-5 needs the $\omega$-dependence of radiation for the blackbody.

## 1. The problem that forced this idea

By 1865 Maxwell had equations; by 1888 Hertz had waves. Three questions remained, and each one broke something.

**How does a changing source's influence spread?** Coulomb's and Biot–Savart's laws are instantaneous: move a charge here and the field changes everywhere at once. Ludvig Lorenz showed in 1867 that the potentials could instead be written as integrals over the source evaluated at an earlier time, delayed by the light travel time. Alfred-Marie Liénard (1898) and Emil Wiechert (1900) independently worked out what this means for a single point charge in arbitrary motion, and found a strange factor that depends on whether the charge is approaching or receding.

**Why do accelerated charges radiate, and how much?** Joseph Larmor in 1897 computed the power radiated by an accelerating charge: proportional to the square of the acceleration. In 1911 Rutherford's scattering experiments showed that an atom is a small positive nucleus with electrons around it. An orbiting electron accelerates, so by Larmor it radiates, loses energy and spirals in. The computation (Example 4.4) gives about $1.6\times10^{-11}$ s. Atoms last forever. This contradiction is where QM-1 starts.

**In which frame do Maxwell's equations hold?** They predict light at speed $c$, and Galilean velocity addition says a moving observer should measure a different speed. So either Maxwell's equations hold only in one frame (the ether's) or Galilean relativity is wrong. Michelson and Morley (1887) could not detect the Earth's motion through any ether. Lorentz (1895–1904) and Poincaré found the transformation that leaves Maxwell's equations unchanged, treating it at first as a property of electrons and the ether. Einstein (1905) made it a statement about space and time: the laws of physics, Maxwell's included, are the same in every inertial frame, and space and time transform so as to make that true. His paper was titled "On the electrodynamics of moving bodies"; its opening example is a magnet and a conductor in relative motion, where the induced current is "electric" in one frame and "magnetic" in the other.

This module does the first two problems completely and the third to the point of showing that electric and magnetic fields are one object seen from different frames.

## 2. Intuition

**Information travels at $c$.** The potential at your location now is built from where the charges were when the light that is just arriving left them. So replace $\rho(\mathbf r')$ in Coulomb's integral by $\rho(\mathbf r', t - \mathscr r/c)$ and guess that you are done. For the potentials, you are (§3.2). For the fields, the naive guess fails, because differentiating the delayed source brings in its rate of change.

**A moving source is counted with a weight.** Picture a train coming straight at you, and suppose you photograph it with light that all arrives at once. Light from the back of the train must leave earlier than light from the front, and the train was further away then, so the photograph shows a longer train than the real one: if the train moves at $v$ toward you, its apparent length is $L/(1 - v/c)$. A charge distribution moving toward you is "seen" over a larger retarded volume, so its retarded integral weighs it by $1/(1 - \hat{\boldsymbol{\mathscr r}}\cdot\mathbf v/c)$. Shrink the train to a point and the factor survives. That is the Liénard–Wiechert factor.

**Radiation is a kink that runs away.** Purcell's picture (Purcell–Morin appendix H, after J. J. Thomson). A charge sits at rest; its field lines are radial. At $t = 0$ it is given a brief kick and then coasts. At a later time $t$, nothing beyond radius $ct$ has heard about the kick: out there the field lines still point radially from where the charge *would* be if it had stayed. Inside radius $c(t - \Delta t)$ the lines point from where the charge now is. In the thin shell between them the lines must connect, so they run nearly *sideways*, tangentially. The shell's thickness is $c\Delta t$, fixed; its radius grows. The tangential field in the shell is set by geometry: the radial field there is $\propto 1/r^2$ and the lines are stretched sideways by a factor proportional to $r$, so the tangential field is $\propto 1/r$. A $1/r$ field has energy density $\propto 1/r^2$ in a shell of area $\propto r^2$: the energy in the shell is the same at every radius. That energy has left for good. Guess: the power goes like (acceleration)$^2$, because the sideways component is proportional to the velocity change $\Delta v$ while the shell thickness is $c\Delta t$, and energy goes like field squared. It is largest broadside to the acceleration ($\sin\theta$ in the field, $\sin^2\theta$ in the power) and zero along it. The project's first figure draws this.

```
           old Coulomb field (charge at origin)
       \     |     /
        \____|____/      <- kink shell, radius ct, lines run sideways
        /    |    \
       /     *-->  \     charge, now coasting
```

**A moving charge's field is squashed, and some of $\mathbf E$ becomes $\mathbf B$.** In the frame of a charge at rest there is only a Coulomb field. In a frame where it moves, the field is compressed along the motion (length contraction squeezes the pattern) and strengthened broadside to it, and the moving charge is a current, so there is a $\mathbf B$. Neither observer is wrong; they disagree about how to split one thing into "electric" and "magnetic." Guess: the components of $\mathbf E$ and $\mathbf B$ along the boost are untouched (nothing is contracted in that direction transverse to the plates of a capacitor) and the transverse components mix, the way $x$ and $t$ mix.

**The field tensor.** $x$ and $ct$ mix under a boost the way $x$ and $y$ mix under a rotation. A vector $\mathbf E$ is not enough to hold an object that mixes with $\mathbf B$: you need something with six components. An antisymmetric $4\times4$ matrix has exactly six independent entries. Put $\mathbf E/c$ in the time–space row and $\mathbf B$ in the space–space block, and the boost rules are matrix conjugation by $\Lambda$, the same law LA-12 used for the stress tensor under rotation (LA-12 (3.19)).

## 3. Formal development

### 3.1 The problem: solve the wave equations for the potentials

In the Lorenz gauge, Maxwell's equations are (EM-6 (3.37))
$$\Box V = -\frac{\rho}{\varepsilon_0},\qquad\Box\mathbf A = -\mu_0\mathbf J. \tag{3.1}$$
In statics these are Poisson's equations, solved by $V = \frac{1}{4\pi\varepsilon_0}\int\rho/\mathscr r\,d\tau'$ and $\mathbf A = \frac{\mu_0}{4\pi}\int\mathbf J/\mathscr r\,d\tau'$ (EM-2, EM-5).

### 3.2 The retarded potentials

$$V(\mathbf r, t) = \frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\mathbf r', t_r)}{\mathscr r}\,d\tau',\qquad\mathbf A(\mathbf r, t) = \frac{\mu_0}{4\pi}\int\frac{\mathbf J(\mathbf r', t_r)}{\mathscr r}\,d\tau',\qquad t_r = t - \frac{\mathscr r}{c}. \tag{3.2}$$

▶ **Try it first.** Compute $\nabla V$. The integrand depends on $\mathbf r$ in two places: through $1/\mathscr r$ and through $t_r$ inside $\rho$. Write $\nabla t_r$ first.

**Theorem 3.1 (retarded potentials).** For $\rho$, $\mathbf J$ smooth and localized, (3.2) satisfies (3.1).

*Derivation (for $V$; $\mathbf A$ is component-by-component identical).* Write $\dot\rho = \partial\rho/\partial t$ evaluated at $(\mathbf r', t_r)$.
1. Since $\mathscr r = |\mathbf r - \mathbf r'|$, $\nabla\mathscr r = \hat{\boldsymbol{\mathscr r}}$, so
$$\nabla t_r = -\frac1c\hat{\boldsymbol{\mathscr r}},\qquad\nabla\rho(\mathbf r', t_r) = \dot\rho\,\nabla t_r = -\frac{\dot\rho}{c}\hat{\boldsymbol{\mathscr r}}. \tag{3.3}$$
This is the step people drop: $\rho(\mathbf r', t_r)$ depends on the field point $\mathbf r$ even though $\mathbf r'$ is the integration variable.
2. Product rule, with $\nabla(1/\mathscr r) = -\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$:
$$\nabla V = \frac{1}{4\pi\varepsilon_0}\int\Big[-\frac{\dot\rho}{c}\frac{\hat{\boldsymbol{\mathscr r}}}{\mathscr r} - \rho\frac{\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\Big]d\tau'.$$
3. Take the divergence of each term. The vector identities needed: $\nabla\cdot(\hat{\boldsymbol{\mathscr r}}/\mathscr r) = 1/\mathscr r^2$ and $\nabla\cdot(\hat{\boldsymbol{\mathscr r}}/\mathscr r^2) = 4\pi\delta^3(\boldsymbol{\mathscr r})$ (EM-1), and by (3.3) $\nabla\dot\rho = -(\ddot\rho/c)\hat{\boldsymbol{\mathscr r}}$, $\nabla\rho = -(\dot\rho/c)\hat{\boldsymbol{\mathscr r}}$.
4. First term: $\nabla\cdot(\dot\rho\,\hat{\boldsymbol{\mathscr r}}/\mathscr r) = \frac{\hat{\boldsymbol{\mathscr r}}}{\mathscr r}\cdot\nabla\dot\rho + \dot\rho\frac{1}{\mathscr r^2} = -\frac{\ddot\rho}{c\mathscr r} + \frac{\dot\rho}{\mathscr r^2}$.
5. Second term: $\nabla\cdot(\rho\,\hat{\boldsymbol{\mathscr r}}/\mathscr r^2) = \frac{\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\cdot\nabla\rho + 4\pi\rho\,\delta^3(\boldsymbol{\mathscr r}) = -\frac{\dot\rho}{c\mathscr r^2} + 4\pi\rho\,\delta^3(\boldsymbol{\mathscr r})$.
6. Assemble:
$$\nabla^2V = \frac{1}{4\pi\varepsilon_0}\int\Big[\frac{\ddot\rho}{c^2\mathscr r} - \frac{\dot\rho}{c\mathscr r^2} + \frac{\dot\rho}{c\mathscr r^2} - 4\pi\rho\,\delta^3(\boldsymbol{\mathscr r})\Big]d\tau' = \frac{1}{c^2}\frac{1}{4\pi\varepsilon_0}\int\frac{\ddot\rho}{\mathscr r}d\tau' - \frac{\rho(\mathbf r, t)}{\varepsilon_0}.$$
The $\dot\rho/\mathscr r^2$ terms cancel; the delta function picks out $\mathbf r' = \mathbf r$, where $t_r = t$.
7. Since $t_r$ depends on $t$ only through $t - \mathscr r/c$, $\partial_t^2$ of (3.2) is the same integral with $\ddot\rho$. So $\nabla^2V = c^{-2}\partial_t^2V - \rho/\varepsilon_0$, which is (3.1). ∎

The retarded potentials also satisfy the Lorenz gauge condition (a consequence of charge conservation; Griffiths Problem 10.8). Project Part 1 checks the point-source case symbolically (Example 4.1) and Part 5 checks the dipole's potentials.

■ **What just happened.** The guess "Coulomb with a delay" is exactly right for the potentials. The derivation had two sources of $\mathbf r$-dependence; the chain rule through $t_r$ produced terms with $\dot\rho$ and $\ddot\rho$, and those are what make $\nabla^2V$ equal $c^{-2}\partial_t^2V$ instead of zero.

**Advanced potentials.** Replace $t_r$ by $t_a = t + \mathscr r/c$ and every step goes through unchanged: the wave equation is symmetric under $t\to-t$, so it has a second solution in which the potential now depends on what the source will do later. We discard it by a physical choice, causality: effects follow causes. That is an input, not a theorem of Maxwell's equations, and Wheeler and Feynman (1945) built a whole theory using half-retarded, half-advanced fields to see how far the choice could be pushed.

### 3.3 Jefimenko's equations

▶ **Try it first.** You already have $\nabla V$ (step 2 above). Add $-\partial_t\mathbf A$ and convert $\mu_0$ to $1/(\varepsilon_0c^2)$.

**Theorem 3.2 (Jefimenko, 1966).**
$$\mathbf E(\mathbf r, t) = \frac{1}{4\pi\varepsilon_0}\int\Big[\frac{\rho(\mathbf r', t_r)}{\mathscr r^2}\hat{\boldsymbol{\mathscr r}} + \frac{\dot\rho(\mathbf r', t_r)}{c\mathscr r}\hat{\boldsymbol{\mathscr r}} - \frac{\dot{\mathbf J}(\mathbf r', t_r)}{c^2\mathscr r}\Big]d\tau', \tag{3.4}$$
$$\mathbf B(\mathbf r, t) = \frac{\mu_0}{4\pi}\int\Big[\frac{\mathbf J(\mathbf r', t_r)}{\mathscr r^2} + \frac{\dot{\mathbf J}(\mathbf r', t_r)}{c\mathscr r}\Big]\times\hat{\boldsymbol{\mathscr r}}\,d\tau'. \tag{3.5}$$

*Derivation of (3.4).*
1. From step 2 of Theorem 3.1, $-\nabla V = \frac{1}{4\pi\varepsilon_0}\int\big[\rho\hat{\boldsymbol{\mathscr r}}/\mathscr r^2 + \dot\rho\hat{\boldsymbol{\mathscr r}}/(c\mathscr r)\big]d\tau'$.
2. $\partial_t\mathbf A = \frac{\mu_0}{4\pi}\int\dot{\mathbf J}/\mathscr r\,d\tau'$, because $\partial t_r/\partial t = 1$ at fixed $\mathbf r$, $\mathbf r'$.
3. $\mathbf E = -\nabla V - \partial_t\mathbf A$ and $\mu_0 = 1/(\varepsilon_0c^2)$ give (3.4). ∎

(3.5) comes the same way from $\nabla\times\mathbf A$, using $\nabla\times[\mathbf J(\mathbf r', t_r)] = \nabla t_r\times\dot{\mathbf J} = -\frac1c\hat{\boldsymbol{\mathscr r}}\times\dot{\mathbf J}$ (Griffiths §10.2.2 writes it out).

■ **What just happened.** The fields are not "Coulomb and Biot–Savart with a delay." They carry extra terms in the retarded $\dot\rho$, $\mathbf J$ and $\dot{\mathbf J}$, falling off only as $1/\mathscr r$. Those $1/\mathscr r$ terms are radiation. In statics they vanish and (3.4)–(3.5) collapse to Coulomb and Biot–Savart. Jefimenko's equations are the direct answer to "what causes what": the fields are produced by the retarded charges and currents and their rates of change. Nothing in them says a changing $\mathbf B$ "causes" $\mathbf E$; Faraday's law is a relation between fields that both come from the sources.

### 3.4 The Liénard–Wiechert potentials

A point charge $q$ moves on $\mathbf w(t)$ with $|\mathbf v| < c$. The retarded time at $(\mathbf r, t)$ is the root of
$$c(t - t_r) = |\mathbf r - \mathbf w(t_r)|. \tag{3.6}$$
**There is exactly one root.** If $t_1 < t_2$ were both roots, subtract: $c(t_2 - t_1) = \big||\mathbf r - \mathbf w(t_1)| - |\mathbf r - \mathbf w(t_2)|\big|\le|\mathbf w(t_2) - \mathbf w(t_1)|$, so the charge would have averaged speed $\ge c$ between them. (At least one root exists because the left side grows without bound as $t_r\to-\infty$ while the right grows at most at rate $v < c$.) The project's solver uses exactly this: $g(t_r) = c(t - t_r) - |\mathbf r - \mathbf w(t_r)|$ is strictly decreasing, so bisection finds the root.

▶ **Try it first.** Write the point charge as $\rho(\mathbf r', t') = q\,\delta^3(\mathbf r' - \mathbf w(t'))$ and insert a delta function to make the retarded time an integration variable: $\rho(\mathbf r', t_r) = \int\rho(\mathbf r', t')\,\delta(t' - t + |\mathbf r - \mathbf r'|/c)\,dt'$. Do the $\mathbf r'$ integral first.

**Theorem 3.3 (Liénard–Wiechert potentials).**
$$V(\mathbf r, t) = \frac{q}{4\pi\varepsilon_0}\frac{1}{\mathscr r\,(1 - \hat{\boldsymbol{\mathscr r}}\cdot\mathbf v/c)},\qquad\mathbf A(\mathbf r, t) = \frac{\mathbf v}{c^2}V(\mathbf r, t), \tag{3.7}$$
with $\boldsymbol{\mathscr r} = \mathbf r - \mathbf w(t_r)$ and $\mathbf v = \dot{\mathbf w}(t_r)$, all at the retarded time.

*Derivation.*
1. $V = \frac{1}{4\pi\varepsilon_0}\int\!\!\int\frac{q\,\delta^3(\mathbf r' - \mathbf w(t'))\,\delta(t' - t + |\mathbf r - \mathbf r'|/c)}{|\mathbf r - \mathbf r'|}\,dt'\,d\tau'$.
2. The $\mathbf r'$ integral sets $\mathbf r' = \mathbf w(t')$: $V = \frac{q}{4\pi\varepsilon_0}\int\frac{\delta(g(t'))}{|\mathbf r - \mathbf w(t')|}dt'$ with $g(t') = t' - t + |\mathbf r - \mathbf w(t')|/c$.
3. $g$ has the single root $t' = t_r$ by (3.6). Use $\delta(g(t')) = \delta(t' - t_r)/|g'(t_r)|$.
4. $g'(t') = 1 + \frac1c\frac{d}{dt'}|\mathbf r - \mathbf w(t')| = 1 - \frac{(\mathbf r - \mathbf w)\cdot\mathbf v}{c|\mathbf r - \mathbf w|} = 1 - \hat{\boldsymbol{\mathscr r}}\cdot\mathbf v/c$, which is positive for $v < c$.
5. So $V = \frac{q}{4\pi\varepsilon_0}\big[\mathscr r(1 - \hat{\boldsymbol{\mathscr r}}\cdot\mathbf v/c)\big]^{-1}$. For $\mathbf A$, $\mathbf J = \rho\mathbf v$ adds a factor $\mathbf v(t_r)$ to the same integral, and $\mu_0\varepsilon_0 = 1/c^2$. ∎

It is convenient to define $\mathbf u = c\hat{\boldsymbol{\mathscr r}} - \mathbf v$, so that $\mathscr r(c - \hat{\boldsymbol{\mathscr r}}\cdot\mathbf v) = \boldsymbol{\mathscr r}\cdot\mathbf u$ and
$$V = \frac{1}{4\pi\varepsilon_0}\frac{qc}{\boldsymbol{\mathscr r}\cdot\mathbf u},\qquad\mathbf A = \frac{\mathbf v}{c^2}V. \tag{3.8}$$

■ **What just happened.** For a point charge the retarded integral is not $q/\mathscr r$: the Jacobian of the delta function is the factor $1/(1 - \hat{\boldsymbol{\mathscr r}}\cdot\mathbf v/c)$, the train of §2. A charge approaching you ($\hat{\boldsymbol{\mathscr r}}\cdot\mathbf v > 0$) counts for more; a receding one for less. It is not a relativistic effect (no $\gamma$ appears); it is pure retardation. Example 4.2 takes the head-on case.

### 3.5 The fields of a point charge in arbitrary motion

Differentiating (3.8) needs $\nabla t_r$ and $\partial t_r/\partial t$ for a moving source, which is where (3.3) no longer applies: now $t_r$ enters through $\mathbf w(t_r)$ as well. Implicit differentiation of (3.6) gives (Exercise B2)
$$\frac{\partial t_r}{\partial t} = \frac{\mathscr rc}{\boldsymbol{\mathscr r}\cdot\mathbf u},\qquad\nabla t_r = -\frac{\boldsymbol{\mathscr r}}{\boldsymbol{\mathscr r}\cdot\mathbf u},$$
and with these Griffiths §10.3.2 carries the differentiation through (two pages; we state the result):
$$\mathbf E(\mathbf r, t) = \frac{q}{4\pi\varepsilon_0}\frac{\mathscr r}{(\boldsymbol{\mathscr r}\cdot\mathbf u)^3}\Big[(c^2 - v^2)\mathbf u + \boldsymbol{\mathscr r}\times(\mathbf u\times\mathbf a)\Big],\qquad\mathbf B = \frac1c\hat{\boldsymbol{\mathscr r}}\times\mathbf E. \tag{3.9}$$
The first term is the **velocity field** (generalized Coulomb), falling as $1/\mathscr r^2$. The second is the **acceleration field**, falling as $1/\mathscr r$ and perpendicular to $\hat{\boldsymbol{\mathscr r}}$: the kink of §2. Project Part 2 evaluates (3.9) numerically on whole grids.

### 3.6 A charge in uniform motion: the squashed Coulomb field

▶ **Try it first.** Set $\mathbf a = 0$. The charge's *present* position is $\mathbf w(t) = \mathbf w(t_r) + \mathbf v(t - t_r)$. Show that $\mathbf R = \mathbf r - \mathbf w(t)$ equals $\boldsymbol{\mathscr r} - \mathscr r\mathbf v/c$, so $\mathscr r\mathbf u = c\mathbf R$.

**Theorem 3.4 (field of a uniformly moving charge).** With $\mathbf R$ from the present position and $\theta$ the angle between $\mathbf R$ and $\mathbf v$,
$$\mathbf E = \frac{q}{4\pi\varepsilon_0}\frac{1 - \beta^2}{(1 - \beta^2\sin^2\theta)^{3/2}}\frac{\hat{\mathbf R}}{R^2},\qquad\mathbf B = \frac{1}{c^2}\mathbf v\times\mathbf E. \tag{3.10}$$

*Derivation.*
1. $\mathbf R = \mathbf r - \mathbf w(t_r) - \mathbf v(t - t_r) = \boldsymbol{\mathscr r} - \mathbf v\mathscr r/c$, so $\mathscr r\mathbf u = c\boldsymbol{\mathscr r} - \mathscr r\mathbf v = c\mathbf R$.
2. Then (3.9) with $\mathbf a = 0$ is $\mathbf E = \frac{q}{4\pi\varepsilon_0}\frac{(c^2 - v^2)c\mathbf R}{(\boldsymbol{\mathscr r}\cdot\mathbf u)^3}$. The field points from the *present* position, not the retarded one.
3. To evaluate $\boldsymbol{\mathscr r}\cdot\mathbf u = c\mathscr r - \boldsymbol{\mathscr r}\cdot\mathbf v$: from $\boldsymbol{\mathscr r} = \mathbf R + \mathscr r\mathbf v/c$, $\boldsymbol{\mathscr r}\cdot\mathbf v = \mathbf R\cdot\mathbf v + \mathscr rv^2/c$, so $\boldsymbol{\mathscr r}\cdot\mathbf u = c\mathscr r(1 - \beta^2) - \mathbf R\cdot\mathbf v$.
4. Squaring $\boldsymbol{\mathscr r} = \mathbf R + \mathscr r\mathbf v/c$: $\mathscr r^2(1 - \beta^2) = R^2 + 2\mathscr r\,\mathbf R\cdot\mathbf v/c$.
5. Square step 3: $(\boldsymbol{\mathscr r}\cdot\mathbf u)^2 = c^2\mathscr r^2(1 - \beta^2)^2 - 2c\mathscr r(1 - \beta^2)\mathbf R\cdot\mathbf v + (\mathbf R\cdot\mathbf v)^2$. Multiply step 4 by $c^2(1 - \beta^2)$ and substitute for the first term: $(\boldsymbol{\mathscr r}\cdot\mathbf u)^2 = c^2(1 - \beta^2)R^2 + (\mathbf R\cdot\mathbf v)^2 = c^2R^2 - v^2R^2\sin^2\theta$.
6. So $\boldsymbol{\mathscr r}\cdot\mathbf u = cR\sqrt{1 - \beta^2\sin^2\theta}$, and step 2 becomes (3.10). For $\mathbf B$: $\hat{\boldsymbol{\mathscr r}} = (\mathbf R + \mathscr r\mathbf v/c)/\mathscr r$ and $\mathbf R\times\mathbf E = 0$, so $\mathbf B = \frac1c\hat{\boldsymbol{\mathscr r}}\times\mathbf E = \frac{1}{c^2}\mathbf v\times\mathbf E$. ∎

■ **What just happened.** All the retardation conspires to make the field point straight from where the charge is *now*, as if the field "knew" to extrapolate. It is weakened by $1 - \beta^2$ along the motion and strengthened by $\gamma$ broadside ($\theta = 90°$): a pancake. In §3.17 we get (3.10) a second way, by boosting the Coulomb field, and Project Part 4 checks the two against each other and against the numerical Liénard–Wiechert solver.

### 3.7 What radiation is

A localized source radiates if energy leaves it and never comes back. Surround it with a sphere of radius $r$; the power through the sphere at time $t$ came from the source at $t_0 = t - r/c$. Define
$$P_{\rm rad}(t_0) = \lim_{r\to\infty}\oint\mathbf S(\mathbf r, t_0 + r/c)\cdot d\mathbf a,\qquad\mathbf S = \frac{1}{\mu_0}\mathbf E\times\mathbf B. \tag{3.11}$$
The area grows as $r^2$, so only the parts of $\mathbf E$ and $\mathbf B$ that fall as $1/r$ can contribute: $1/r^2$ fields give $\mathbf S\sim1/r^4$ and the integral vanishes. **Radiation means the $1/r$ fields.** A changing field is not necessarily radiation: the near field of a slowly varying dipole changes in time but its energy sloshes back and forth (its time-averaged flux through a large sphere is zero). In (3.9) only the acceleration term is $1/r$, so a uniformly moving charge does not radiate, consistent with the fact that in its own frame it is a static charge.

### 3.8 Electric dipole radiation

Model: two point charges $\pm q(t) = \pm q_0\cos\omega t$ at $z = \pm d/2$ joined by a wire carrying current $I(t) = \dot q = -q_0\omega\sin\omega t$, so $\mathbf p(t) = p_0\cos\omega t\,\hat{\mathbf z}$ with $p_0 = q_0d$. Spherical coordinates $(r, \theta, \phi)$ about the centre.

The three approximations, which must be used in this order:
$$\text{(i) } d\ll r\ \text{(point dipole)},\qquad\text{(ii) } d\ll\lambda = \frac{2\pi c}{\omega}\ \text{(the source is small compared with a wavelength)},\qquad\text{(iii) } r\gg\lambda\ \text{(radiation zone)}. \tag{3.12}$$

▶ **Try it first.** Write $V$ as the sum of two retarded Coulomb potentials with $\mathscr r_\pm = \sqrt{r^2\mp rd\cos\theta + d^2/4}$. Expand to first order in $d$.

*Derivation.*
1. $V = \frac{1}{4\pi\varepsilon_0}\Big[\frac{q_0\cos\omega(t - \mathscr r_+/c)}{\mathscr r_+} - \frac{q_0\cos\omega(t - \mathscr r_-/c)}{\mathscr r_-}\Big]$.
2. By (i), $\mathscr r_\pm\approx r(1\mp\frac{d}{2r}\cos\theta)$ and $1/\mathscr r_\pm\approx\frac1r(1\pm\frac{d}{2r}\cos\theta)$.
3. $\cos\omega(t - \mathscr r_\pm/c)\approx\cos\big(\omega(t - r/c)\pm\frac{\omega d}{2c}\cos\theta\big)$. By (ii), $\omega d/c\ll1$, so expand: $\approx\cos\omega(t - r/c)\mp\frac{\omega d}{2c}\cos\theta\,\sin\omega(t - r/c)$.
4. Multiply steps 2 and 3, keep first order in $d$, subtract:
$$V = \frac{p_0\cos\theta}{4\pi\varepsilon_0r}\Big[-\frac{\omega}{c}\sin\omega(t - r/c) + \frac1r\cos\omega(t - r/c)\Big].$$
In the static limit $\omega\to0$ this is the dipole potential of EM-2.
5. By (iii), $1/r\ll\omega/c$: drop the second term.
$$V = -\frac{p_0\omega}{4\pi\varepsilon_0c}\frac{\cos\theta}{r}\sin\omega(t - r/c),\qquad\mathbf A = -\frac{\mu_0p_0\omega}{4\pi r}\sin\omega(t - r/c)\,\hat{\mathbf z}. \tag{3.13}$$
For $\mathbf A$: $\mathbf A = \frac{\mu_0}{4\pi}\int_{-d/2}^{d/2}\frac{I(t - \mathscr r/c)}{\mathscr r}dz\,\hat{\mathbf z}$; to zeroth order in $d$ the integrand is constant, giving $d\,I(t - r/c)/r$, and $dI = -p_0\omega\sin$.
6. Fields, keeping only $1/r$ terms (derivatives that hit $1/r$ or $\cos\theta$ produce $1/r^2$ and are dropped by (iii)): $\nabla V\approx\hat{\mathbf r}\,\partial_rV = \frac{p_0\omega^2}{4\pi\varepsilon_0c^2}\frac{\cos\theta}{r}\cos\omega(t - r/c)\hat{\mathbf r}$ and $\partial_t\mathbf A = -\frac{\mu_0p_0\omega^2}{4\pi r}\cos\omega(t - r/c)(\cos\theta\,\hat{\mathbf r} - \sin\theta\,\hat{\boldsymbol\theta})$, using $\hat{\mathbf z} = \cos\theta\,\hat{\mathbf r} - \sin\theta\,\hat{\boldsymbol\theta}$. In $-\nabla V - \partial_t\mathbf A$ the $\hat{\mathbf r}$ parts cancel ($\varepsilon_0c^2 = 1/\mu_0$):
$$\mathbf E = -\frac{\mu_0p_0\omega^2}{4\pi}\frac{\sin\theta}{r}\cos\omega(t - r/c)\,\hat{\boldsymbol\theta},\qquad\mathbf B = -\frac{\mu_0p_0\omega^2}{4\pi c}\frac{\sin\theta}{r}\cos\omega(t - r/c)\,\hat{\boldsymbol\phi}. \tag{3.14}$$
$\mathbf B = \nabla\times\mathbf A$ keeps only the $\partial_r$ of the phase: $\nabla\times(A_z\hat{\mathbf z})\approx\hat{\mathbf r}\times\hat{\mathbf z}\,\partial_rA_z = -\sin\theta\,\hat{\boldsymbol\phi}\,\partial_rA_z$.
7. Poynting vector: $\mathbf S = \frac{1}{\mu_0}\mathbf E\times\mathbf B = \frac{\mu_0}{c}\Big[\frac{p_0\omega^2}{4\pi}\frac{\sin\theta}{r}\cos\omega(t - r/c)\Big]^2\hat{\mathbf r}$ ($\hat{\boldsymbol\theta}\times\hat{\boldsymbol\phi} = \hat{\mathbf r}$). Average over a period ($\langle\cos^2\rangle = \tfrac12$):
$$\langle\mathbf S\rangle = \frac{\mu_0p_0^2\omega^4}{32\pi^2c}\frac{\sin^2\theta}{r^2}\hat{\mathbf r}. \tag{3.15}$$
8. Total power: $\int_0^\pi\sin^2\theta\,\sin\theta\,d\theta = 4/3$, times $2\pi$ from $\phi$:
$$\langle P\rangle = \frac{\mu_0p_0^2\omega^4}{12\pi c}. \tag{3.16}$$
∎

■ **What just happened.** In the radiation zone $\mathbf E$, $\mathbf B$ and $\hat{\mathbf r}$ form EM-7's plane-wave triad with $E/B = c$, locally a plane wave, with amplitude $\propto\sin\theta/r$. The pattern is a donut: nothing along the axis, maximum broadside. Power is independent of $r$, as (3.11) requires. The $\omega^4$ is the blue sky: EM-7 §3.13 said the molecular dipoles are driven far below resonance so $p_0$ does not depend on $\omega$; now (3.16) supplies the $\omega^4$, and $(700/450)^4\approx5.9$ times more blue than red is scattered (Exercise A3). Project Part 3 computes $\langle S_r\rangle$ from the exact numerical fields of an oscillating charge and confirms the $\sin^2\theta$ pattern, the radius-independence of $P$, and the exponent 4.

**Magnetic dipole radiation** (stated). A loop of radius $b$ with current $I_0\cos\omega t$ has $\mathbf m = m_0\cos\omega t\,\hat{\mathbf z}$, $m_0 = \pi b^2I_0$. The same three approximations give $\mathbf E$ along $\hat{\boldsymbol\phi}$ and $\mathbf B$ along $\hat{\boldsymbol\theta}$, the same $\sin^2\theta$ donut, and
$$\langle P\rangle_{\rm mag} = \frac{\mu_0m_0^2\omega^4}{12\pi c^3}. \tag{3.17}$$
For a loop and a dipole of comparable size $b\sim d$ and the same current, $P_{\rm mag}/P_{\rm el}\sim(\omega b/c)^2$, suppressed by $(b/\lambda)^2$ (Exercise B4). That is why electric dipole radiation dominates whenever it is allowed, and why "forbidden" atomic transitions (QM-8) are slow.

**Arbitrary localized source.** Repeat the derivation for any $\rho$, $\mathbf J$ of size $\ll\lambda$, expanding $t_r$ about $t_0 = t - r/c$ to first order in $\mathbf r'\cdot\hat{\mathbf r}/c$ (Griffiths §11.2.1). The result, to dipole order, depends only on the second derivative of the dipole moment at the retarded time:
$$\mathbf E = \frac{\mu_0}{4\pi r}\big[\hat{\mathbf r}\times(\hat{\mathbf r}\times\ddot{\mathbf p})\big]_{t_0},\qquad P_{\rm rad}(t_0) = \frac{\mu_0\ddot p(t_0)^2}{6\pi c}. \tag{3.18}$$
For $\ddot p = -p_0\omega^2\cos\omega t$, $\langle\ddot p^2\rangle = p_0^2\omega^4/2$ and (3.18) reproduces (3.16).

### 3.9 The Larmor formula

▶ **Try it first.** At an instant when the charge is at rest ($\mathbf v = 0$, $\mathbf a\ne0$), simplify the acceleration term of (3.9). What is $\mathbf u$?

**Theorem 3.5 (Larmor, 1897).** A point charge with $v\ll c$ radiates
$$P = \frac{\mu_0q^2a^2}{6\pi c}. \tag{3.19}$$

*Derivation.*
1. With $\mathbf v = 0$, $\mathbf u = c\hat{\boldsymbol{\mathscr r}}$ and $\boldsymbol{\mathscr r}\cdot\mathbf u = c\mathscr r$. The acceleration field of (3.9) is $\mathbf E_{\rm rad} = \frac{q}{4\pi\varepsilon_0}\frac{\mathscr r}{c^3\mathscr r^3}\boldsymbol{\mathscr r}\times(c\hat{\boldsymbol{\mathscr r}}\times\mathbf a) = \frac{\mu_0q}{4\pi\mathscr r}\big[\hat{\boldsymbol{\mathscr r}}\times(\hat{\boldsymbol{\mathscr r}}\times\mathbf a)\big]$.
2. Let $\theta$ be the angle between $\hat{\boldsymbol{\mathscr r}}$ and $\mathbf a$. $|\hat{\boldsymbol{\mathscr r}}\times(\hat{\boldsymbol{\mathscr r}}\times\mathbf a)| = a\sin\theta$. With $\mathbf B = \hat{\boldsymbol{\mathscr r}}\times\mathbf E/c$ perpendicular to $\mathbf E$, $\mathbf S = \frac{E^2}{\mu_0c}\hat{\boldsymbol{\mathscr r}} = \frac{\mu_0q^2a^2}{16\pi^2c}\frac{\sin^2\theta}{\mathscr r^2}\hat{\boldsymbol{\mathscr r}}$.
3. Integrate over the sphere: $\frac{\mu_0q^2a^2}{16\pi^2c}\cdot2\pi\cdot\frac43 = \frac{\mu_0q^2a^2}{6\pi c}$. ∎

The velocity field contributes nothing at large $\mathscr r$ (it is $1/\mathscr r^2$). The same result follows from (3.18) with $\mathbf p = q\mathbf w$, $\ddot{\mathbf p} = q\mathbf a$.

■ **What just happened.** Power goes as $a^2$, exactly the §2 guess. The radiated power is what the kink of §2 carries; the $\sin^2\theta$ is the sideways projection of $\mathbf a$.

**Liénard's generalization** (stated; Griffiths §11.2.1 derives it). For arbitrary speed,
$$P = \frac{\mu_0q^2\gamma^6}{6\pi c}\Big(a^2 - \frac{|\mathbf v\times\mathbf a|^2}{c^2}\Big). \tag{3.20}$$
The $\gamma^6$ makes ultrarelativistic radiation enormous. For circular motion ($\mathbf a\perp\mathbf v$) the bracket is $a^2/\gamma^2$ and $P\propto\gamma^4a^2$; the project's Your-turn 1 checks this by integrating Liénard's angular distribution.

### 3.10 Consequences: the classical atom, synchrotrons, bremsstrahlung, antennas

**The classical atom.** An electron in a circular orbit of radius $r$ around a proton has $a = e^2/(4\pi\varepsilon_0m_er^2)$ and energy $E = -e^2/(8\pi\varepsilon_0r)$. Setting $dE/dt = -P$ from (3.19) gives an ODE for $r(t)$ whose solution reaches $r = 0$ from the Bohr radius in $t_f = a_0^3/(4r_e^2c)\approx1.6\times10^{-11}$ s, where $r_e = e^2/(4\pi\varepsilon_0m_ec^2)$ is the classical electron radius (Example 4.4 does the derivation and the numbers). The emitted frequency would sweep continuously upward. Real atoms are stable and emit sharp lines. Classical electrodynamics plus Rutherford's atom is inconsistent with the existence of matter. QM-1 starts here.

**Synchrotron radiation.** A relativistic charge in a magnetic field circles with $\mathbf a\perp\mathbf v$ and radiates $P\propto\gamma^4$ by (3.20). The emission is beamed into a forward cone of half-angle $\sim1/\gamma$ (Exercise C2, Your-turn 1), so a distant observer sees a flash once per orbit lasting a fraction $\sim1/\gamma^3$ of the period, whose spectrum extends to very high harmonics. For electron storage rings this is an X-ray source; for electron–positron colliders it is the energy-loss problem that capped LEP and is why the highest-energy colliders accelerate protons ($m_p/m_e = 1836$, and $P\propto\gamma^4$ at fixed energy goes as $m^{-4}$).

**Bremsstrahlung** ("braking radiation"). An electron decelerating in the Coulomb field of a nucleus radiates by (3.19)–(3.20). An X-ray tube is this: electrons accelerated through tens of kV slam into a metal target, and the continuous part of the spectrum is bremsstrahlung, cut off at the photon energy equal to the electron's kinetic energy (the cutoff itself is quantum; classically nothing stops the spectrum).

**Antennas.** A centre-fed wire of length $L$ carrying a standing-wave current is a sum of short dipoles with retarded phases. The half-wave dipole ($L = \lambda/2$, current $I_0\cos kz$) has intensity $\propto\cos^2(\tfrac\pi2\cos\theta)/\sin^2\theta$, a slightly sharper donut than $\sin^2\theta$, and radiation resistance $R = \frac{Z_0}{2\pi}\int_0^\pi\frac{\cos^2(\frac\pi2\cos\theta)}{\sin\theta}d\theta = 73.1\ \Omega$, $Z_0 = \mu_0c$ (project Part 6). Matching that 73 Ω to a 50 or 75 Ω cable is a standard EE design problem. Short antennas ($L\ll\lambda$) have tiny radiation resistance $\propto(L/\lambda)^2$ and are inefficient, which is why AM broadcast towers are hundreds of metres tall.

### 3.11 Radiation reaction: where classical E&M breaks

A radiating charge loses energy, so some force must act back on it. Require that the work done by a force $\mathbf F_{\rm rad}$ over an interval in which the motion is periodic (or starts and ends with $\mathbf a = 0$) equals minus the Larmor energy:
$$\int_{t_1}^{t_2}\mathbf F_{\rm rad}\cdot\mathbf v\,dt = -\frac{\mu_0q^2}{6\pi c}\int_{t_1}^{t_2}\mathbf a\cdot\mathbf a\,dt = -\frac{\mu_0q^2}{6\pi c}\Big[\mathbf v\cdot\mathbf a\Big]_{t_1}^{t_2} + \frac{\mu_0q^2}{6\pi c}\int_{t_1}^{t_2}\dot{\mathbf a}\cdot\mathbf v\,dt,$$
integrating by parts. The boundary term vanishes for such intervals, and the simplest force that makes the integrands agree is the **Abraham–Lorentz force**
$$\mathbf F_{\rm rad} = \frac{\mu_0q^2}{6\pi c}\dot{\mathbf a} = m\tau\dot{\mathbf a},\qquad\tau = \frac{\mu_0q^2}{6\pi mc} = 6.27\times10^{-24}\ \text{s (electron)}. \tag{3.21}$$
This "derivation" only matches averages, so it is a plausibility argument; Griffiths §11.2.3 gets the same force from the self-force of a small charged shell, which shows the problem is structural, not an artifact of the argument.

The pathologies are real. With no external force, $m\mathbf a = m\tau\dot{\mathbf a}$ has the solution $\mathbf a = \mathbf a_0e^{t/\tau}$: a free electron spontaneously accelerating forever (**runaway**). Forbid runaways by a boundary condition at $t\to+\infty$, and a force switched on at $t = 0$ produces an acceleration $a(t) = (F_0/m)e^{t/\tau}$ for $t < 0$: the charge responds *before* the force acts (**pre-acceleration**, acausal over a time $\tau$). Project Your-turn 3 integrates both. The time $\tau$ corresponds to a length $c\tau = \tfrac23r_e\approx1.9\times10^{-15}$ m, where classical electrodynamics already disagrees with quantum mechanics by orders of magnitude, so the honest conclusion is that classical point-charge electrodynamics is not self-consistent at that scale. The Landau–Lifshitz reduction (replace $\dot{\mathbf a}$ by its value computed from the external force) is the practical fix. QED is the real one.

### 3.12 Special relativity: the Lorentz transformation

**Postulates (Einstein, 1905).** (1) The laws of physics take the same form in every inertial frame. (2) The speed of light in vacuum is $c$ in every inertial frame, independent of the motion of the source.

▶ **Try it first.** Assume the transformation from $S$ to $\bar S$ is linear and that the origin of $\bar S$ ($\bar x = 0$) moves as $x = vt$. Write $\bar x = \gamma(x - vt)$ with unknown $\gamma$. What does postulate 1 say about the inverse? What does postulate 2 say about a light pulse $x = ct$?

**Theorem 3.6 (Lorentz transformation).**
$$\bar t = \gamma\Big(t - \frac{v}{c^2}x\Big),\qquad\bar x = \gamma(x - vt),\qquad\bar y = y,\qquad\bar z = z,\qquad\gamma = \frac{1}{\sqrt{1 - v^2/c^2}}. \tag{3.22}$$

*Derivation.*
1. *Linearity.* Homogeneity of space and time (no point or instant is special) forces the map to be affine; the shared origin makes it linear.
2. *Transverse directions.* Let two rods of equal rest length, one in each frame, pass each other perpendicular to $\hat{\mathbf x}$ with paint brushes at their ends. If $\bar y = ky$ with $k\neq1$, one rod would mark the other, and the reverse situation (by postulate 1 the same with $k$) would have it marked the other way; the only consistent answer is $k = 1$. So $\bar y = y$, $\bar z = z$.
3. *Along $x$.* The point $\bar x = 0$ moves as $x = vt$, so $\bar x = \gamma(x - vt)$ for some $\gamma(v)$. By postulate 1, $S$ moves at $-v$ relative to $\bar S$, so $x = \gamma(\bar x + v\bar t)$ with the same $\gamma$ (isotropy: $\gamma(-v) = \gamma(v)$).
4. *Light.* A pulse leaves the origin at $t = \bar t = 0$. By postulate 2, $x = ct$ and $\bar x = c\bar t$. Substitute into step 3: $c\bar t = \gamma(c - v)t$ and $ct = \gamma(c + v)\bar t$.
5. Multiply: $c^2t\bar t = \gamma^2(c^2 - v^2)t\bar t$, so $\gamma = 1/\sqrt{1 - v^2/c^2}$ (positive, so that $\gamma\to1$ as $v\to0$).
6. *Time.* From $x = \gamma(\bar x + v\bar t)$ and $\bar x = \gamma(x - vt)$: $v\bar t = x/\gamma - \gamma(x - vt)$, so $\bar t = \gamma t + x(1/\gamma - \gamma)/v = \gamma t - \gamma\frac{v}{c^2}x$, using $1/\gamma^2 - 1 = -v^2/c^2$. ∎

In four-vector form, $\bar x^\mu = \Lambda^\mu{}_\nu x^\nu$ with
$$\Lambda = \begin{pmatrix}\gamma & -\gamma\beta & 0 & 0\\ -\gamma\beta & \gamma & 0 & 0\\ 0 & 0 & 1 & 0\\ 0 & 0 & 0 & 1\end{pmatrix}. \tag{3.23}$$

■ **What just happened.** Linearity plus "the same law in both directions" fixed the form; one light pulse fixed $\gamma$. Time had to transform too: the term $-\gamma vx/c^2$ is what makes $c$ frame-independent.

**Consequences.**
- *Relativity of simultaneity.* Two events at $t = 0$ at $x_1\ne x_2$ have $\bar t_2 - \bar t_1 = -\gamma v(x_2 - x_1)/c^2\ne0$. Simultaneity depends on the frame.
- *Time dilation.* A clock at rest in $\bar S$ ($\bar x$ fixed) ticks $\Delta\bar t$; from the inverse transformation $t = \gamma(\bar t + v\bar x/c^2)$, $\Delta t = \gamma\Delta\bar t$. Moving clocks run slow. Cosmic-ray muons ($\tau = 2.2\ \mu$s at rest) reach the ground because of this (Exercise A1).
- *Length contraction.* A rod at rest in $\bar S$ with length $\bar L$: in $S$ mark both ends at the same $t$; $\Delta\bar x = \gamma\Delta x$, so $L = \bar L/\gamma$. Only the length along the motion contracts.
- *Velocity addition.* Differentiating (3.22), an object moving at $u$ in $S$ moves at $\bar u = (u - v)/(1 - uv/c^2)$ in $\bar S$. A useful identity used twice below (project Part 1 checks it): for $w = (u + v)/(1 + uv/c^2)$, $\gamma_w = \gamma_u\gamma_v(1 + uv/c^2)$. Proof: $1 - w^2/c^2 = \frac{(1 + uv/c^2)^2 - (u + v)^2/c^2}{(1 + uv/c^2)^2} = \frac{(1 - u^2/c^2)(1 - v^2/c^2)}{(1 + uv/c^2)^2}$.

### 3.13 The interval, four-vectors, and $E^2 = p^2c^2 + m^2c^4$

The **interval** between two events separated by $\Delta x^\mu$ is
$$\Delta s^2 = \eta_{\mu\nu}\Delta x^\mu\Delta x^\nu = -c^2\Delta t^2 + \Delta x^2 + \Delta y^2 + \Delta z^2. \tag{3.24}$$
It is invariant: substitute (3.22), $-c^2\gamma^2(\Delta t - v\Delta x/c^2)^2 + \gamma^2(\Delta x - v\Delta t)^2 = \gamma^2(1 - \beta^2)(\Delta x^2 - c^2\Delta t^2)$. In matrix form, $\Lambda^\top\eta\Lambda = \eta$: Lorentz transformations are the linear maps that preserve the indefinite form $\eta$, as rotations preserve $\delta_{ij}$. By Sylvester's law (LA-12 §3.8) the signature $(1, 3)$ is basis-independent, which is why "timelike" ($\Delta s^2 < 0$), "spacelike" ($> 0$) and "lightlike" ($= 0$) mean the same in every frame.

A **four-vector** is any $a^\mu$ that transforms like $x^\mu$: $\bar a^\mu = \Lambda^\mu{}_\nu a^\nu$. Lowering with $\eta$ (LA-12 (3.16)) gives $a_\mu = (-a^0, \mathbf a)$, and $a_\mu b^\mu = -a^0b^0 + \mathbf a\cdot\mathbf b$ is invariant. **Proper time** along a worldline, $d\tau = \sqrt{-ds^2}/c = dt/\gamma$, is the time on a clock carried along. **Four-velocity** $u^\mu = dx^\mu/d\tau = \gamma(c, \mathbf v)$ has $u_\mu u^\mu = \gamma^2(-c^2 + v^2) = -c^2$.

▶ **Try it first.** Define $p^\mu = mu^\mu$. Its spatial part is $\gamma m\mathbf v$. Show that $\frac{d}{dt}(\gamma mc^2) = \mathbf F\cdot\mathbf v$ when $\mathbf F = d(\gamma m\mathbf v)/dt$, so that $p^0c$ is the energy.

**Theorem 3.7 (energy and momentum).** With relativistic momentum $\mathbf p = \gamma m\mathbf v$ and $\mathbf F = d\mathbf p/dt$, the energy (work–energy theorem with $E = mc^2$ at rest) is $E = \gamma mc^2$; $p^\mu = (E/c, \mathbf p)$; and
$$E^2 = p^2c^2 + m^2c^4. \tag{3.25}$$

*Derivation.*
1. $d\gamma/dt = \gamma^3\mathbf v\cdot\mathbf a/c^2$ (differentiate $(1 - v^2/c^2)^{-1/2}$).
2. $\mathbf F\cdot\mathbf v = \frac{d}{dt}(\gamma m\mathbf v)\cdot\mathbf v = m\gamma\,\mathbf a\cdot\mathbf v + m\gamma^3\frac{(\mathbf v\cdot\mathbf a)v^2}{c^2} = m\gamma^3(\mathbf v\cdot\mathbf a)\Big(\frac{1}{\gamma^2} + \frac{v^2}{c^2}\Big) = m\gamma^3\mathbf v\cdot\mathbf a$.
3. By step 1 this is $\frac{d}{dt}(\gamma mc^2)$. So the work done equals the change of $\gamma mc^2$, which is the energy up to a constant; the constant is fixed by requiring $p^\mu$ to be a four-vector, $E = p^0c = \gamma mc^2$.
4. $p_\mu p^\mu = m^2u_\mu u^\mu = -m^2c^2$. Also $p_\mu p^\mu = -E^2/c^2 + p^2$. Equate and multiply by $-c^2$. ∎

(Why $\mathbf p = \gamma m\mathbf v$ and not $m\mathbf v$: $m\mathbf v$ is not conserved in all frames if it is conserved in one; $\gamma m\mathbf v$ is. Taylor §15.9 does the collision argument.) Project Part 1 checks (3.25) symbolically.

### 3.14 Minkowski diagrams in one page

Plot $ct$ up and $x$ across. A light ray is a $45°$ line. A body at rest is vertical; one moving at $v$ is a line at slope $c/v$ from the horizontal, steeper than $45°$. The $\bar t$ axis ($\bar x = 0$) is the line $x = vt$; the $\bar x$ axis ($\bar t = 0$) is $ct = \beta x$. The two primed axes tilt symmetrically toward the light line, like closing scissors.

```
   ct          ct-bar
    |         /
    |        /      . light (45 deg)
    |       /    .
    |      /  .
    |     / .     _ _ x-bar  (ct = beta x)
    |    /.  _ -
    |   .- -
    |_-_______________ x
```

Read off the effects. *Simultaneity:* lines of constant $\bar t$ are parallel to the tilted $\bar x$ axis, so two events on a horizontal line ($t$ equal) are at different $\bar t$. *Time dilation:* the tick of a clock at rest in $\bar S$ is on the $\bar t$ axis; drop a horizontal line to read it in $S$ and the coordinate time is longer. *Causality:* events inside each other's light cone ($\Delta s^2 < 0$) have the same time order in every frame, since no boost tips a timelike line past $45°$. Spacelike-separated events can have either order, which is why nothing can travel faster than light without some frame seeing the effect before the cause. The hyperbolae $-c^2t^2 + x^2 = $ const are invariant and calibrate the tick marks on the tilted axes. Spend twenty minutes drawing these with Taylor §15 or Morin's relativity book; the diagrams are the fastest way to stop being confused by twin and pole-in-barn paradoxes.

### 3.15 Magnetism as relativity: Purcell's wire

▶ **Try it first.** A neutral wire: positive ions at rest with line density $\lambda_0$, conduction electrons with density $-\lambda_0$ drifting at speed $v$ in $-\hat{\mathbf x}$, so the current is $I = \lambda_0v$ in $+\hat{\mathbf x}$. A test charge $q > 0$ at distance $s$ moves at $u$ in $+\hat{\mathbf x}$. In its rest frame it feels no magnetic force. Where does a force come from?

**Theorem 3.8.** In the test charge's rest frame $\bar S$ the wire has line charge $\bar\lambda = -\gamma_uIu/c^2$, and the resulting electric force, transformed back to $S$, is exactly the magnetic force $quB$, $B = \mu_0I/(2\pi s)$, attractive.

*Derivation.*
1. *Lab force.* In $S$ the wire is neutral, $B = \mu_0I/(2\pi s)$ circles it, and $\mathbf F = q\mathbf u\times\mathbf B$ has magnitude $F = quB = \frac{\mu_0qIu}{2\pi s}$, toward the wire (parallel currents attract).
2. *Charge is invariant and length contracts.* A line of charges with density $\lambda_{\rm rest}$ in its own rest frame has density $\gamma\lambda_{\rm rest}$ in a frame where it moves at speed with factor $\gamma$: the same charges occupy a shorter length.
3. *Ions.* At rest in $S$, so $\lambda_{\rm rest} = \lambda_0$. In $\bar S$ they move at $-u$: $\bar\lambda_+ = \gamma_u\lambda_0$.
4. *Electrons.* They have density $-\lambda_0$ in $S$ where they move at $v$, so their rest density is $-\lambda_0/\gamma_v$. In $\bar S$ their speed is $w$ with $|w| = (u + v)/(1 + uv/c^2)$ (velocity addition), so $\bar\lambda_- = -\gamma_w\lambda_0/\gamma_v = -\gamma_u(1 + uv/c^2)\lambda_0$ by the identity of §3.12.
5. *Net charge in $\bar S$:* $\bar\lambda = \gamma_u\lambda_0 - \gamma_u(1 + uv/c^2)\lambda_0 = -\gamma_u\lambda_0uv/c^2 = -\gamma_uIu/c^2$. The wire is negatively charged in the test charge's frame: the electrons, moving faster there than the ions, are more contracted.
6. *Force in $\bar S$.* The charge is at rest, so only $\bar E$ acts: $\bar E = \bar\lambda/(2\pi\varepsilon_0s)$ (EM-2's line charge; $s$ is transverse, so not contracted), and $\bar F = q|\bar\lambda|/(2\pi\varepsilon_0s) = \gamma_u\frac{qIu}{2\pi\varepsilon_0c^2s} = \gamma_u\frac{\mu_0qIu}{2\pi s}$, toward the wire.
7. *Transforming forces.* Transverse momentum is unchanged by a boost along $x$ ($\bar y = y$), so $dp_y = d\bar p_y$; the time intervals for the particle are $dt = \gamma_ud\bar t$ (time dilation, the particle is at rest in $\bar S$). So $F_y = \bar F_y/\gamma_u = \frac{\mu_0qIu}{2\pi s}$, equal to step 1. ∎

■ **What just happened.** The magnetic force in one frame is the electric force in another. The imbalance $uv/c^2$ is fantastically small (about $10^{-14}$ in Example 4.5), but it multiplies the enormous charge density of conduction electrons (about $10^4$ C/m in a 1 mm² copper wire), and the product is the familiar magnetic force. Magnetism is the relativistic correction to electrostatics that you can feel with your hands, because matter is so precisely neutral that the electric part cancels. Exercise B6 does the same argument from the electrons' frame, where the test charge moves and the forces must cancel.

### 3.16 How $\mathbf E$ and $\mathbf B$ transform

Assume that $\mathbf E$ and $\mathbf B$ at an event in $\bar S$ depend only on $\mathbf E$ and $\mathbf B$ at that event in $S$, and that charge is invariant. Then we can learn the rules from any configuration that produces given fields.

▶ **Try it first.** A parallel-plate capacitor with plates perpendicular to $\hat{\mathbf y}$ moves along $x$. Which surface density changes, and why does a $B_z$ appear?

**Theorem 3.9 (field transformation for a boost along $x$).**
$$\bar E_x = E_x,\quad\bar E_y = \gamma(E_y - vB_z),\quad\bar E_z = \gamma(E_z + vB_y),\qquad\bar B_x = B_x,\quad\bar B_y = \gamma\Big(B_y + \frac{v}{c^2}E_z\Big),\quad\bar B_z = \gamma\Big(B_z - \frac{v}{c^2}E_y\Big). \tag{3.26}$$

*Derivation (Griffiths' route).*
1. *The capacitor, plates $\perp\hat{\mathbf y}$.* In its rest frame $S_0$ the plates carry $\pm\sigma_0$ (bottom $+$ at $y = 0$, top $-$), $\mathbf E_0 = (\sigma_0/\varepsilon_0)\hat{\mathbf y}$, no $\mathbf B$. In $S$ the plates move at $v_0\hat{\mathbf x}$; length contraction along $x$ gives $\sigma = \gamma_0\sigma_0$, so $E_y = \sigma/\varepsilon_0$. The moving plates are surface currents $\pm\sigma v_0\hat{\mathbf x}$; each sheet contributes $\frac{\mu_0}2\mathbf K\times\hat{\mathbf n}$ (EM-5), and between the plates both give $+\hat{\mathbf z}$: $B_z = \mu_0\sigma v_0$. So $E_y = \gamma_0\sigma_0/\varepsilon_0$, $B_z = \mu_0\gamma_0\sigma_0v_0$.
2. In $\bar S$ the plates move at $\bar v_0 = (v_0 - v)/(1 - vv_0/c^2)$, so by the same reasoning $\bar E_y = \bar\gamma_0\sigma_0/\varepsilon_0$ and $\bar B_z = \mu_0\bar\gamma_0\sigma_0\bar v_0$. By §3.12's identity (with $-v$), $\bar\gamma_0 = \gamma\gamma_0(1 - vv_0/c^2)$.
3. Then $\bar E_y = \gamma\gamma_0\sigma_0/\varepsilon_0 - \gamma v\,\frac{\gamma_0\sigma_0v_0}{\varepsilon_0c^2} = \gamma(E_y - vB_z)$, using $1/(\varepsilon_0c^2) = \mu_0$.
4. And $\bar B_z = \mu_0\gamma\gamma_0\sigma_0(1 - vv_0/c^2)\frac{v_0 - v}{1 - vv_0/c^2} = \gamma(\mu_0\gamma_0\sigma_0v_0 - \mu_0\gamma_0\sigma_0v) = \gamma\big(B_z - \frac{v}{c^2}E_y\big)$.
5. *Plates $\perp\hat{\mathbf z}$.* Rotate the set-up by $90°$ about $x$ ($y\to z$, $z\to-y$): the current sheets now give $B_y = -\mu_0\sigma v_0$, and the same algebra gives $\bar E_z = \gamma(E_z + vB_y)$, $\bar B_y = \gamma(B_y + \frac{v}{c^2}E_z)$.
6. *Plates $\perp\hat{\mathbf x}$.* Motion along $x$ does not contract the plates' area, so $\sigma$ and $E_x = \sigma/\varepsilon_0$ are unchanged (the separation contracts, but the field does not depend on it). Charges moving along the plate normal make no current sheet in the plane, so there is no $B$.
7. *A solenoid along $x$.* $B_x = \mu_0nI$. In $\bar S$ the turns are closer, $\bar n = \gamma n$, and the current, charge per unit time through a wire whose clocks run slow by $\gamma$, is $\bar I = I/\gamma$. So $\bar B_x = B_x$. ∎

In vector form, with $\parallel$ and $\perp$ relative to $\mathbf v$: $\bar{\mathbf E}_\parallel = \mathbf E_\parallel$, $\bar{\mathbf B}_\parallel = \mathbf B_\parallel$, $\bar{\mathbf E}_\perp = \gamma(\mathbf E + \mathbf v\times\mathbf B)_\perp$, $\bar{\mathbf B}_\perp = \gamma(\mathbf B - \mathbf v\times\mathbf E/c^2)_\perp$. The first of these is the Lorentz force law in disguise: a charge moving with $\bar S$ feels $q(\mathbf E + \mathbf v\times\mathbf B)$ in $S$ and $q\bar{\mathbf E}$ in its own frame.

■ **What just happened.** Fields along the boost are untouched; fields across it mix $\mathbf E$ with $\mathbf v\times\mathbf B$, with a $\gamma$, the way $x$ mixes with $t$. Check on the moving charge: start from a pure Coulomb field in the rest frame and apply the inverse rules ($v\to-v$, no $\bar{\mathbf B}$): $E_\parallel = \bar E_\parallel$, $\mathbf E_\perp = \gamma\bar{\mathbf E}_\perp$, $\mathbf B = \mathbf v\times\mathbf E/c^2$. Since $\bar x = \gamma R_x$ at $t = 0$, $\mathbf E = \frac{q\gamma\mathbf R}{4\pi\varepsilon_0(\gamma^2R_x^2 + R_\perp^2)^{3/2}}$, and $\gamma^2R_x^2 + R_\perp^2 = \gamma^2R^2(1 - \beta^2\sin^2\theta)$ gives (3.10) again. Project Part 4 checks this numerically.

**Theorem 3.10 (invariants).** $E^2 - c^2B^2$ and $\mathbf E\cdot\mathbf B$ are the same in every frame.

*Proof from (3.26).*
1. The $x$ components are unchanged, so consider $\perp$ parts. $\bar E_y^2 + \bar E_z^2 = \gamma^2[(E_y - vB_z)^2 + (E_z + vB_y)^2]$ and $c^2(\bar B_y^2 + \bar B_z^2) = \gamma^2[(cB_y + \beta E_z)^2 + (cB_z - \beta E_y)^2]$.
2. Subtract. The cross terms are $\gamma^2[-2vE_yB_z + 2vE_zB_y - 2c\beta B_yE_z + 2c\beta B_zE_y] = 0$ since $c\beta = v$. The squares give $\gamma^2[(1 - \beta^2)(E_y^2 + E_z^2) - c^2(1 - \beta^2)(B_y^2 + B_z^2)] = E_\perp^2 - c^2B_\perp^2$.
3. $\bar E_y\bar B_y + \bar E_z\bar B_z = \gamma^2[(E_y - vB_z)(B_y + vE_z/c^2) + (E_z + vB_y)(B_z - vE_y/c^2)]$. Expanding, the $vE_yE_z/c^2$ terms cancel, the $vB_yB_z$ terms cancel, and what remains is $\gamma^2(1 - \beta^2)(E_yB_y + E_zB_z)$. ∎

Consequences: if $\mathbf E\perp\mathbf B$ in one frame they are perpendicular in all; if $E > cB$ in one frame there is no frame with a pure magnetic field; a plane wave ($\mathbf E\perp\mathbf B$, $E = cB$) is a plane wave for every observer. Exercise B7 finds the frame in which a crossed $\mathbf E\perp\mathbf B$ with $E < cB$ becomes purely magnetic.

### 3.17 The field tensor

▶ **Try it first.** Write the antisymmetric matrix $F^{\mu\nu}$ with $F^{0i} = E_i/c$, $F^{ij} = \varepsilon_{ijk}B_k$. Compute the $(0, 2)$ entry of $\Lambda F\Lambda^\top$ with (3.23).

$$F^{\mu\nu} = \begin{pmatrix}0 & E_x/c & E_y/c & E_z/c\\ -E_x/c & 0 & B_z & -B_y\\ -E_y/c & -B_z & 0 & B_x\\ -E_z/c & B_y & -B_x & 0\end{pmatrix},\qquad G^{\mu\nu} = \begin{pmatrix}0 & B_x & B_y & B_z\\ -B_x & 0 & -E_z/c & E_y/c\\ -B_y & E_z/c & 0 & -E_x/c\\ -B_z & -E_y/c & E_x/c & 0\end{pmatrix}. \tag{3.27}$$
$G$ (the **dual tensor**) is obtained from $F$ by $\mathbf E/c\to\mathbf B$, $\mathbf B\to-\mathbf E/c$. Warning: Jackson and most field-theory texts use $(+,-,-,-)$ with $F^{0i} = -E_i$; the physics is the same and every sign in this section flips somewhere (§8).

**Theorem 3.11.** The rule for a rank-2 contravariant tensor (LA-12 (3.18)),
$$\bar F^{\mu\nu} = \Lambda^\mu{}_\alpha\Lambda^\nu{}_\beta F^{\alpha\beta},\quad\text{i.e. }\bar F = \Lambda F\Lambda^\top, \tag{3.28}$$
reproduces (3.26).

*Verification of three components.* $\Lambda$ has $\Lambda^0{}_0 = \Lambda^1{}_1 = \gamma$, $\Lambda^0{}_1 = \Lambda^1{}_0 = -\gamma\beta$, $\Lambda^2{}_2 = \Lambda^3{}_3 = 1$, all others zero.
1. $\bar F^{01} = \Lambda^0{}_\alpha\Lambda^1{}_\beta F^{\alpha\beta}$; only $\alpha, \beta\in\{0, 1\}$ contribute and $F^{00} = F^{11} = 0$: $\bar F^{01} = \Lambda^0{}_0\Lambda^1{}_1F^{01} + \Lambda^0{}_1\Lambda^1{}_0F^{10} = \gamma^2F^{01} - \gamma^2\beta^2F^{01} = F^{01}$. So $\bar E_x = E_x$.
2. $\bar F^{02} = \Lambda^0{}_\alpha\Lambda^2{}_2F^{\alpha2} = \gamma F^{02} - \gamma\beta F^{12} = \gamma(E_y/c - \beta B_z)$. So $\bar E_y = \gamma(E_y - vB_z)$.
3. $\bar F^{12} = \Lambda^1{}_\alpha F^{\alpha2} = -\gamma\beta F^{02} + \gamma F^{12} = \gamma(B_z - \beta E_y/c)$. So $\bar B_z = \gamma(B_z - vE_y/c^2)$. ∎

The other three components follow the same way; Project Part 4 checks all six for 20 random fields and boosts in random directions.

**Invariants as contractions.** Lowering both indices flips the sign of the $0i$ entries: $F_{0i} = -F^{0i}$, $F_{ij} = F^{ij}$. Then
$$F_{\mu\nu}F^{\mu\nu} = 2\Big(B^2 - \frac{E^2}{c^2}\Big),\qquad F_{\mu\nu}G^{\mu\nu} = -\frac4c\,\mathbf E\cdot\mathbf B. \tag{3.29}$$
(First: $2\sum_iF_{0i}F^{0i} = -2E^2/c^2$, $\sum_{ij}F_{ij}F^{ij} = \varepsilon_{ijk}\varepsilon_{ijl}B_kB_l = 2B^2$.) A full contraction of tensors is a scalar, so Theorem 3.10 is now one line.

■ **What just happened.** Six numbers that two observers disagree about are the components of one object they agree about, exactly as $x$ and $t$ are components of one event. The matrix form makes the mixing obvious: a boost along $x$ mixes rows/columns 0 and 1, so it leaves $F^{01}$ ($E_x$) and $F^{23}$ ($B_x$) alone and mixes the pairs $(F^{02}, F^{12})$ and $(F^{03}, F^{13})$.

### 3.18 Maxwell's equations in covariant form

**Four-current and continuity.** $J^\mu = (c\rho, \mathbf J) = \rho_0u^\mu$ for a charge density $\rho_0$ in its rest frame, so it is a four-vector (the $\gamma$ from length contraction multiplies $\rho_0$). Continuity, $\partial_t\rho + \nabla\cdot\mathbf J = 0$, is $\partial_\mu J^\mu = 0$, an invariant statement.

▶ **Try it first.** Write out $\partial_\nu F^{0\nu}$ and $\partial_\nu F^{1\nu}$ with (3.27) and $\partial_0 = c^{-1}\partial_t$.

**Theorem 3.12 (covariant Maxwell equations).**
$$\partial_\nu F^{\mu\nu} = \mu_0J^\mu,\qquad\partial_\nu G^{\mu\nu} = 0. \tag{3.30}$$
These are the four Maxwell equations.

*Derivation (components).*
1. $\mu = 0$: $\partial_\nu F^{0\nu} = \partial_iF^{0i} = \frac1c\nabla\cdot\mathbf E$. Set equal to $\mu_0c\rho$: $\nabla\cdot\mathbf E = \mu_0c^2\rho = \rho/\varepsilon_0$. Gauss.
2. $\mu = 1$: $\partial_0F^{10} + \partial_2F^{12} + \partial_3F^{13} = -\frac{1}{c^2}\partial_tE_x + \partial_yB_z - \partial_zB_y = \mu_0J_x$. That is the $x$ component of $\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\partial_t\mathbf E$. Ampère–Maxwell. ($\mu = 2, 3$ are the other components by cyclic relabelling.)
3. $G$, $\mu = 0$: $\partial_iG^{0i} = \nabla\cdot\mathbf B = 0$.
4. $G$, $\mu = 1$: $\partial_0G^{10} + \partial_2G^{12} + \partial_3G^{13} = -\frac1c\partial_tB_x - \frac1c\partial_yE_z + \frac1c\partial_zE_y = -\frac1c\big(\partial_t\mathbf B + \nabla\times\mathbf E\big)_x = 0$. Faraday. ∎

**Continuity follows.** Apply $\partial_\mu$ to the first of (3.30): $\partial_\mu\partial_\nu F^{\mu\nu} = 0$ because $\partial_\mu\partial_\nu$ is symmetric and $F^{\mu\nu}$ antisymmetric (swap the dummy names and get minus itself). So $\partial_\mu J^\mu = 0$: Maxwell's equations are consistent only if charge is conserved, the same observation that forced Maxwell's displacement current in EM-6, now in one line.

Note the index placement. The brief form you may see, $\partial_\mu F^{\mu\nu} = \mu_0J^\nu$, is correct in the $(+,-,-,-)$ convention with $F^{0i} = -E_i/c$; in Griffiths' convention it has the opposite sign, since $\partial_\mu F^{\mu\nu} = -\partial_\mu F^{\nu\mu}$. Check any source's convention by working out $\nu = 0$ and seeing whether you get $+\rho/\varepsilon_0$.

**The four-potential.** $A^\mu = (V/c, \mathbf A)$ and
$$F^{\mu\nu} = \partial^\mu A^\nu - \partial^\nu A^\mu,\qquad\partial^\mu = \eta^{\mu\nu}\partial_\nu = (-c^{-1}\partial_t, \nabla). \tag{3.31}$$
*Two components.* $F^{01} = \partial^0A^1 - \partial^1A^0 = -\frac1c\partial_tA_x - \frac1c\partial_xV = \frac1c(-\nabla V - \partial_t\mathbf A)_x = E_x/c$. ✓ $F^{12} = \partial_xA_y - \partial_yA_x = (\nabla\times\mathbf A)_z = B_z$. ✓ Project Part 5 checks all six with SymPy for the exact potentials of an oscillating dipole. With (3.31), $\partial_\nu G^{\mu\nu} = 0$ holds identically (it is $\nabla\cdot\nabla\times = 0$ and $\nabla\times\nabla = 0$ in disguise). The Lorenz gauge is $\partial_\mu A^\mu = \frac{1}{c^2}\partial_tV + \nabla\cdot\mathbf A = 0$, manifestly invariant, and in it $\partial_\nu F^{\mu\nu} = \partial^\mu(\partial_\nu A^\nu) - \partial_\nu\partial^\nu A^\mu = -\Box A^\mu$, so (3.30) becomes $\Box A^\mu = -\mu_0J^\mu$: (3.1) in one line, with $\Box = \partial_\nu\partial^\nu$ an invariant operator. The retarded potentials of §3.2 are its solution in every frame, which is why the Liénard–Wiechert potentials (3.8) form a four-vector: $A^\mu = \frac{\mu_0qc}{4\pi}\frac{u^\mu}{-x_\nu u^\nu}$ in disguise, with $x^\nu$ the separation to the retarded point (Exercise C3 asks you to check it).

**The Lorentz force.** The **Minkowski force** $K^\mu = dp^\mu/d\tau$ on a charge is
$$K^\mu = q\,u_\nu F^{\mu\nu}. \tag{3.32}$$
Check $\mu = 1$: $u_\nu = \gamma(-c, \mathbf v)$, so $K^1 = q[u_0F^{10} + u_2F^{12} + u_3F^{13}] = q\gamma[(-c)(-E_x/c) + v_yB_z - v_zB_y] = q\gamma(\mathbf E + \mathbf v\times\mathbf B)_x$. Since $K^1 = dp_x/d\tau = \gamma\,dp_x/dt$, this is the Lorentz force law. ($\mu = 0$ gives $q\gamma\mathbf v\cdot\mathbf E/c$, the rate of work.)

### 3.19 The Lagrangian density (one paragraph)

CM-2 obtained equations of motion from an action $\int L\,dt$; F6 did the calculus of variations for functions of one variable. A field theory uses an action $S = \int\mathcal L\,d^4x$ and Euler–Lagrange equations $\partial_\nu\frac{\partial\mathcal L}{\partial(\partial_\nu A_\mu)} = \frac{\partial\mathcal L}{\partial A_\mu}$, one for each field component (the same derivation as F6 with one integration by parts per coordinate). For electromagnetism in our signature,
$$\mathcal L = -\frac{1}{4\mu_0}F_{\mu\nu}F^{\mu\nu} + A_\mu J^\mu = \frac{\varepsilon_0}{2}E^2 - \frac{1}{2\mu_0}B^2 - \rho V + \mathbf J\cdot\mathbf A, \tag{3.33}$$
using (3.29) and $A_\mu J^\mu = -V\rho + \mathbf A\cdot\mathbf J$. (In the $(+,-,-,-)$ convention the coupling is written $-A_\mu J^\mu$; it is the same function.) The last two terms are the charged-particle Lagrangian $-qV + q\mathbf v\cdot\mathbf A$ of CM-2 §3.9 per unit volume. *The $\mu = 0$ equation* (vary $V$): $\mathcal L$ depends on $\partial_iV$ only through $\mathbf E = -\nabla V - \partial_t\mathbf A$, so $\partial\mathcal L/\partial(\partial_iV) = \varepsilon_0E_i\cdot(-1)$, and $\partial\mathcal L/\partial V = -\rho$. Euler–Lagrange: $\partial_i(-\varepsilon_0E_i) = -\rho$, that is, $\nabla\cdot\mathbf E = \rho/\varepsilon_0$. Varying $\mathbf A$ gives Ampère–Maxwell (Exercise C3); covariantly, $\partial\mathcal L/\partial(\partial_\nu A_\mu) = -\frac{1}{\mu_0}F^{\nu\mu}$ and $\partial\mathcal L/\partial A_\mu = J^\mu$ give $\partial_\nu F^{\mu\nu} = \mu_0J^\mu$ at once. The homogeneous equations are not Euler–Lagrange equations; they are automatic once $F = \partial A - \partial A$.

### 3.20 What comes next

Classical electrodynamics is a complete, relativistic field theory, and it fails in two places this module has met: radiating atoms (§3.10) and point-charge self-interaction (§3.11). QM-1 resolves the first by abandoning classical orbits; the blackbody spectrum (SM-5) and the photoelectric effect show that the energy in each mode of the field (3.30) comes in quanta $\hbar\omega$. The photon is the quantum of $A^\mu$. Quantum electrodynamics is the quantum field theory of (3.33) with the Dirac electron, and it computes the electron's magnetic moment to about one part in $10^{12}$.

## 4. Worked examples

### Example 4.1: The retarded potential of a point source solves the wave equation (project Part 1)

Take $V(\mathbf r, t) = \frac{1}{4\pi\varepsilon_0}\frac{f(t - r/c)}{r}$, the retarded potential of a point source of strength $f(t)$ at the origin. (A pure point charge cannot change, by charge conservation; think of $f$ as the charge on a small sphere fed by a thin wire, and look only away from the wire.) For $r > 0$ we use the radial Laplacian $\nabla^2V = \frac1r\partial_r^2(rV)$, valid for spherically symmetric functions.
1. $rV = \frac{1}{4\pi\varepsilon_0}f(t - r/c)$. $\partial_r(rV) = -\frac{1}{4\pi\varepsilon_0c}f'$, $\partial_r^2(rV) = \frac{1}{4\pi\varepsilon_0c^2}f''$.
2. $\nabla^2V = \frac{1}{4\pi\varepsilon_0c^2}\frac{f''}{r}$. And $\partial_t^2V = \frac{1}{4\pi\varepsilon_0}\frac{f''}{r}$. So $\nabla^2V - c^{-2}\partial_t^2V = 0$ for $r > 0$.
3. Near $r = 0$, $V\approx f(t)/(4\pi\varepsilon_0r)$ and $\nabla^2(1/r) = -4\pi\delta^3(\mathbf r)$, so $\Box V = -f(t)\delta^3(\mathbf r)/\varepsilon_0$: the source is $\rho = f(t)\delta^3(\mathbf r)$.
The project checks step 2 with SymPy in Cartesian coordinates for an arbitrary function $f$. This is the Green's function of $\Box$: (3.2) is a superposition of these with $f = \rho\,d\tau'$.

### Example 4.2: The retarded-volume factor for a charge coming straight at you (project Part 2)

A charge $q = 1$ nC moves along $+\hat{\mathbf x}$ at $\beta = 0.6$, passing the origin at $t = 0$. The observer is at $x = 1$ m. Find the retarded time, the retarded distance, $V$ and $E$ at $t = 0$.
1. Retarded condition (3.6): $c(0 - t_r) = 1 - vt_r$ (the charge is to the left), so $t_r = -\frac{1}{c - v} = -\frac{1}{0.4c} = -8.339$ ns.
2. $\mathscr r = -ct_r = 1/(1 - \beta) = 2.500$ m. The light arriving now left when the charge was 2.5 m away, not 1 m.
3. $\hat{\boldsymbol{\mathscr r}}\cdot\mathbf v/c = \beta = 0.6$, so the Liénard–Wiechert factor is $1/(1 - \beta) = 2.5$ and (3.7) gives $V = \frac{q}{4\pi\varepsilon_0}\frac{1}{2.5\times0.4} = \frac{q}{4\pi\varepsilon_0\cdot1\text{ m}} = 8.988$ V. The naive "Coulomb at the retarded position" would be $8.988/2.5 = 3.595$ V. The factor exactly undoes the retardation: $\mathscr r(1 - \beta) = R$, the present distance.
4. The field from (3.10) with $\theta = 0$ (on the line of motion), $R = 1$ m: $E_x = (1 - \beta^2)\frac{q}{4\pi\varepsilon_0R^2} = 0.64\times8.988 = 5.752$ V/m, pointing away from the charge's *present* position.
The project's numerical retarded-time solver reproduces $t_r$, $\mathscr r$, $V$ and $E_x$ to $10^{-12}$.

### Example 4.3: A 1 m antenna at 100 MHz carrying 1 A (project Part 6)

Treat a 1 m wire with uniform current amplitude $I_0 = 1$ A at $f = 100$ MHz as a Hertzian dipole.
1. $\omega = 2\pi f = 6.283\times10^8$ rad/s, $\lambda = c/f = 2.998$ m. So $d/\lambda = 0.334$: approximation (ii) is marginal and the answer is an estimate.
2. Charge at the ends: $I = \dot q$ gives $q_0 = I_0/\omega$, so $p_0 = q_0d = I_0d/\omega = 1/(6.283\times10^8) = 1.592\times10^{-9}$ C·m.
3. (3.16): $\langle P\rangle = \frac{\mu_0p_0^2\omega^4}{12\pi c} = \frac{\mu_0I_0^2d^2\omega^2}{12\pi c} = 43.90$ W.
4. Radiation resistance, defined by $\langle P\rangle = \tfrac12I_0^2R$: $R = 87.79\ \Omega$. Equivalently $R = \frac{2\pi}{3}Z_0(d/\lambda)^2 = 789\ \Omega\times0.1113$.
5. Broadside at 1 km, (3.15): $\langle S\rangle = 5.240\times10^{-6}$ W/m², and the field amplitude from (3.14) is $\mu_0p_0\omega^2/(4\pi r) = 0.0628$ V/m. FM radios work with microvolts per metre, so this is a strong signal.
6. Honesty check: a real centre-fed wire's current falls to zero at the ends. With a triangular current profile, $p_0$ halves and the power drops by 4. A half-wave dipole (1.5 m here) with its cosine current has $R = 73.1\ \Omega$ (§3.10). The order of magnitude, tens of ohms and tens of watts per ampere squared, is right.

### Example 4.4: The classical hydrogen atom collapses (project Part 6)

1. On a circular orbit of radius $r$, Newton gives $m_ev^2/r = e^2/(4\pi\varepsilon_0r^2)$, so $a = e^2/(4\pi\varepsilon_0m_er^2)$ and $E = \tfrac12m_ev^2 - \frac{e^2}{4\pi\varepsilon_0r} = -\frac{e^2}{8\pi\varepsilon_0r}$.
2. At $r = a_0 = 5.292\times10^{-11}$ m: $a = 9.044\times10^{22}$ m/s², $E = -13.61$ eV, and Larmor (3.19): $P = 4.669\times10^{-8}$ W. (Fine: $v/c = 1/137$, so (3.19) applies.)
3. $\frac{dE}{dt} = \frac{e^2}{8\pi\varepsilon_0r^2}\frac{dr}{dt} = -\frac{\mu_0e^2}{6\pi c}\Big(\frac{e^2}{4\pi\varepsilon_0m_er^2}\Big)^2$, assuming the orbit stays nearly circular (each orbit loses a tiny fraction of the energy, checked below).
4. Solve for $\dot r$ and use $\mu_0 = 1/(\varepsilon_0c^2)$ and $r_e = e^2/(4\pi\varepsilon_0m_ec^2) = 2.818\times10^{-15}$ m: $\dot r = -\frac{e^4}{12\pi^2\varepsilon_0^2m_e^2c^3r^2} = -\frac43\frac{r_e^2c}{r^2}$.
5. Separate: $r^2dr = -\frac43r_e^2c\,dt$, so $r^3 = a_0^3 - 4r_e^2ct$, which reaches zero at
$$t_f = \frac{a_0^3}{4r_e^2c} = 1.556\times10^{-11}\ \text{s}.$$
6. The initial orbital period is $T_0 = 2\pi a_0/v = 1.520\times10^{-16}$ s. Since $r = a_0(1 - t/t_f)^{1/3}$ and $T\propto r^{3/2}$, the number of orbits is $\int_0^{t_f}dt/T = 2t_f/T_0 = 2.05\times10^5$. Each orbit loses on average a fraction $\sim10^{-5}$ of the energy, so the "nearly circular" assumption holds until the very end.
The project integrates the ODE in step 4 numerically and recovers $t_f$ to $10^{-6}$.

### Example 4.5: Purcell's wire in numbers (project Part 6)

A wire carries $I = 10$ A. A proton-charge test particle ($q = e$) at $s = 1$ cm moves parallel to the current at $u = 10^6$ m/s.
1. Lab: $B = \mu_0I/(2\pi s) = 2.000\times10^{-4}$ T; $F = quB = 3.2044\times10^{-17}$ N, toward the wire.
2. Test charge's frame, Theorem 3.8 step 5: $\bar\lambda = -\gamma_uIu/c^2 = -1.113\times10^{-10}$ C/m.
3. $\bar E = \bar\lambda/(2\pi\varepsilon_0s) = -200.0$ V/m (pointing into the wire); $\bar F = e|\bar E| = 3.20437\times10^{-17}$ N.
4. $\gamma_u - 1 = 5.56\times10^{-6}$, and $\bar F/\gamma_u = 3.20435\times10^{-17}$ N, equal to step 1 (the project asserts agreement to $10^{-10}$; the residual is the difference between CODATA's $\mu_0$ and $1/(\varepsilon_0c^2)$).
5. Where the charge comes from, for a 1 mm² copper wire ($n = 8.5\times10^{28}$ m⁻³): $\lambda_0 = nAe = 1.362\times10^4$ C/m, drift speed $v = I/\lambda_0 = 7.34\times10^{-4}$ m/s, and the fractional imbalance $uv/c^2 = 8.2\times10^{-15}$. Fifteen orders of magnitude of cancellation leave a force you can measure on a bench.

### Example 4.6: Boosting a pure electric field (project Part 4)

In $S$, $\mathbf E = E_0\hat{\mathbf x}$ with $E_0 = 1.0\times10^6$ V/m and $\mathbf B = 0$. Find the fields in a frame moving at $\mathbf v = 0.6c\,\hat{\mathbf y}$.
1. The boost is along $y$, so $E_x$ is transverse and the vector form of (3.26) applies: $\bar{\mathbf E}_\perp = \gamma(\mathbf E + \mathbf v\times\mathbf B)_\perp = \gamma E_0\hat{\mathbf x}$, $\bar{\mathbf B} = -\gamma\mathbf v\times\mathbf E/c^2 = -\gamma\frac{vE_0}{c^2}(\hat{\mathbf y}\times\hat{\mathbf x}) = \gamma\frac{vE_0}{c^2}\hat{\mathbf z}$.
2. $\gamma = 1.25$: $\bar{\mathbf E} = 1.25\times10^6\,\hat{\mathbf x}$ V/m and $\bar B_z = 1.25\times0.6\times10^6/c = 2.502\times10^{-3}$ T. A magnetic field appears where there was none.
3. Invariants: $\bar E^2 - c^2\bar B^2 = \gamma^2E_0^2(1 - \beta^2) = E_0^2 = 1.000\times10^{12}$ V²/m², and $\bar{\mathbf E}\cdot\bar{\mathbf B} = 0 = \mathbf E\cdot\mathbf B$.
4. Interpretation: the source of $\mathbf E$ (say, capacitor plates perpendicular to $x$) moves at $-0.6c\,\hat{\mathbf y}$ in $\bar S$, so its charges are currents along $-\hat{\mathbf y}$ in sheets perpendicular to $x$, which produce $B_z$.
The project computes the same numbers as $\Lambda F\Lambda^\top$ with a boost along $y$.

## 5. Why it matters and where it connects

- **QM-1.** Example 4.4 is the crisis that opens quantum mechanics; (3.16) with $p_0$ replaced by a transition dipole matrix element is the classical skeleton of the emission rate.
- **QM-8.** Spontaneous emission rates are (3.18) with $\langle\ddot p^2\rangle\to$ quantum matrix elements; magnetic dipole and electric quadrupole transitions are suppressed by $(a_0/\lambda)^2$, the same suppression as (3.17).
- **SM-5.** Rayleigh–Jeans counts field modes and gives each $k_BT$; the $\omega^2$ mode density times $k_BT$ diverges, the ultraviolet catastrophe. The field being quantized is the $A^\mu$ of §3.18.
- **CM-2, F6.** (3.33) is the first field Lagrangian in the curriculum; the Euler–Lagrange machinery is F6's with four independent variables.
- **LA-12.** $\Lambda^\top\eta\Lambda = \eta$ is an isometry of an indefinite symmetric form (LA-12 §3.8, Sylvester); (3.28) is the type-(2,0) transformation law; the invariants (3.29) are full contractions.
- **EM-7.** The radiation-zone fields (3.14) are locally EM-7's plane waves; the Rayleigh $\omega^4$ is now derived.
- **Engineering.** Antenna radiation resistance and patterns (Example 4.3, Your-turn 2's phased array is the principle of radar beam steering and 5G beamforming); synchrotron light sources; EMC: any fast current edge on a circuit board radiates by (3.18), which is why rise times and loop areas are regulated.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: muons, a 1 GeV proton, the blue sky, Larmor for a slow electron, a field transformation, index gymnastics with $F_{\mu\nu}$. B1–B9 are the module: B1 (switched-on current: retarded potentials by hand), B2 (the retarded-time derivatives, the error magnet), B3 (squashed field numbers), B5 (the collapse derivation and orbit count), B6 (Purcell from the electrons' frame) and B8 (covariant Maxwell and continuity) are essential; B9 extends `project.py` to a circulating charge. C1 (Abraham–Lorentz), C2 (relativistic beaming via aberration) and C3 (Ampère–Maxwell from the Lagrangian, and $A^\mu$ of a moving charge) go further.

## 7. Coding project

`project.py` runs in about 4 seconds; figures go to `figures/`.

1. **Part 1, SymPy:** Example 4.1 for arbitrary $f$; $\Lambda^\top\eta\Lambda = \eta$; $u_\mu u^\mu = -c^2$; (3.25); the $\gamma$ addition identity.
2. **Part 2, the retarded-time solver:** `retarded_time` brackets and bisects (3.6) on whole arrays of field points; `lw_fields` evaluates (3.8)–(3.9). Example 4.2. A charge kicked from rest to $0.6c$: field lines with the kink shell, and a check that outside the light sphere the field is still the old Coulomb field (`radiation_kink.png`, which also shows $r|\mathbf E|$ for an oscillating charge).
3. **Part 3, dipole radiation:** period-averaged $S_r$ from the numerical fields on spheres of radius $0.3\lambda$ to $50\lambda$: the $\sin^2\theta$ pattern, the same power through every sphere, agreement with (3.16) and (3.19), and the fitted exponent $4.0000$ (`dipole_pattern.png`).
4. **Part 4, the field tensor:** (3.27)–(3.29) for random fields and boosts; Example 4.6; the uniformly moving charge's numerical Liénard–Wiechert field equals both (3.10) and the boosted Coulomb field.
5. **Part 5, covariant Maxwell:** SymPy builds $F^{\mu\nu}$ from (3.31) for the exact dipole potentials and matches $-\nabla V - \partial_t\mathbf A$ and $\nabla\times\mathbf A$; Lorenz gauge; fourth-order finite differences show $\partial_\nu F^{\mu\nu} = 0$ and $\partial_\nu G^{\mu\nu} = 0$ off the source.
6. **Part 6:** Examples 4.3–4.5 and the half-wave dipole's 73 Ω.
7. **Your turn** (reference solutions included): synchrotron beaming from Liénard's angular distribution, two-element beam steering, the Abraham–Lorentz runaway and pre-acceleration. Part 7 prints every number in `solutions.md`, including Exercise B9.

**Build it yourself:** before reading `retarded_time`, write your own: prove that $g(t_r) = c(t - t_r) - |\mathbf r - \mathbf w(t_r)|$ is decreasing, bracket the root, and bisect, vectorized over a grid. Then write `lw_fields` from (3.9) and reproduce Example 4.2 before running anything else. Then build $F^{\mu\nu}$ and $\Lambda$ and check (3.26) by matrix multiplication.

## 8. Where you'll struggle

- **Retarded-time derivatives.** This is the error magnet. $\rho(\mathbf r', t_r)$ depends on $\mathbf r$ through $t_r$, and for a moving point charge $t_r$ depends on $\mathbf r$ and $t$ implicitly through $\mathbf w(t_r)$, so $\nabla t_r = -\boldsymbol{\mathscr r}/(\boldsymbol{\mathscr r}\cdot\mathbf u)$, not $-\hat{\boldsymbol{\mathscr r}}/c$. The fix: write $t_r$ explicitly as a function before differentiating, and do Exercise B2 slowly once. If a derivation of the Liénard–Wiechert fields comes out without $(\boldsymbol{\mathscr r}\cdot\mathbf u)^3$ in the denominator, you dropped a $\nabla t_r$.
- **The three dipole approximations.** Apply them in the order (i), (ii), (iii). If you drop the $1/r^2$ term before expanding in $d$, or set $\omega d/c = 0$ before expanding the cosine, you lose the radiation term entirely (it is first order in $\omega d/c$). Also, (ii) fails for real antennas of size comparable to $\lambda$ (Example 4.3 has $d/\lambda = 1/3$); quote such numbers as estimates.
- **"Radiation" means $1/r$ fields carrying energy to infinity.** It does not mean "any changing field." The near field of an oscillating dipole changes in time and carries no net energy away; a uniformly moving charge does not radiate at all. Test: does the time-averaged flux through a sphere survive as $r\to\infty$?
- **Index notation and the metric will cost you a day.** $x_\mu$ and $x^\mu$ differ in the sign of the time component; $\partial_\mu$ has a lower index although it acts like $1/x^\mu$; $\partial^0 = -\partial_0$. Fix: when in doubt, write the component out. LA-12 is the prerequisite for exactly this reason; re-read LA-12 §3.16 if raising and lowering is not automatic.
- **Sign conventions for $F^{\mu\nu}$ differ.** Griffiths: $(-,+,+,+)$, $F^{0i} = E_i/c$, $\partial_\nu F^{\mu\nu} = \mu_0J^\mu$. Jackson (Gaussian) and Peskin–Schroeder: $(+,-,-,-)$, $F^{0i} = -E_i$, $\partial_\mu F^{\mu\nu} = \frac{4\pi}{c}J^\nu$ or $J^\nu$. Some books define $F^{ij} = -\varepsilon_{ijk}B_k$. Never copy a formula between books without checking one component (the $\nu = 0$ equation must give $+\rho/\varepsilon_0$).
- **Relativity is covered at the minimum.** Enough to see $\mathbf E$ and $\mathbf B$ unify. You have not done paradoxes, the Doppler effect, relativistic collisions or acceleration. If the Minkowski diagram does not feel solid, do Taylor ch. 15 or Morin's *Special Relativity for the Enthusiastic Beginner* before QM.
- **The Abraham–Lorentz pathologies are real.** Runaways and pre-acceleration are not a pedagogical simplification. The classical theory of a point charge interacting with its own field is inconsistent; the practical fix (Landau–Lifshitz) works because the effect matters only on time scales $\sim10^{-23}$ s where quantum mechanics has taken over anyway.
- **Which frame's force?** In Purcell's argument you must transform the force back (step 7), and forget that and the two frames disagree by $\gamma_u$. Forces are not components of a four-vector; $K^\mu = \gamma(\mathbf F\cdot\mathbf v/c, \mathbf F)$ is.

## 9. Check yourself

1. Derive (3.3) and explain why the $\dot\rho/\mathscr r^2$ terms cancel in Theorem 3.1.
2. Why do we discard the advanced potentials? Is that forced by Maxwell's equations?
3. Derive the factor $1/(1 - \hat{\boldsymbol{\mathscr r}}\cdot\mathbf v/c)$ from the delta-function integral. Is it relativistic?
4. Which term of (3.9) is radiation, and why is a uniformly moving charge's field not radiation?
5. List the three dipole approximations in the order used and say what each one drops.
6. Derive $\langle P\rangle = \mu_0p_0^2\omega^4/(12\pi c)$ from (3.15), and Larmor from the acceleration field.
7. Estimate the classical atom's lifetime from Larmor in three lines.
8. Derive the Lorentz transformation from the two postulates.
9. Derive $E^2 = p^2c^2 + m^2c^4$ from $p^\mu p_\mu$.
10. In Purcell's argument, why is the wire charged in the test charge's frame, and how does the force transform back?
11. Write $F^{\mu\nu}$ and show that $\Lambda F\Lambda^\top$ gives $\bar E_y = \gamma(E_y - vB_z)$.
12. Expand $\partial_\nu F^{\mu\nu} = \mu_0J^\mu$ for $\mu = 0$ and $\mu = 1$, and show it implies $\partial_\mu J^\mu = 0$.

## 10. Pair with

- **Griffiths, *Introduction to Electrodynamics*, ch. 10–12.** The same route with more examples: the full differentiation of the Liénard–Wiechert potentials (§10.3.2), the general localized source (§11.1.4), radiation reaction from the self-force (§11.2.3), and the capacitor–solenoid field transformations (§12.3.2) that §3.16 follows.
- **Purcell & Morin, *Electricity and Magnetism*, ch. 5 and appendix H.** Electrostatics of moving charges done the relativity-first way, with the wire argument in full, and the Thomson kink construction of a charge's radiation field drawn properly.
- **Feynman, *Lectures on Physics* Vol. I ch. 15–17 and 34, Vol. II ch. 21, 25–26.** I-15–17 is the clearest short special relativity; I-34 and II-21 do the fields of moving charges and radiation from Feynman's own formula; II-25–26 does four-vectors and the field tensor with physical commentary on every step.
- **Taylor, *Classical Mechanics*, ch. 15.** Special relativity at the level you need, with the collision derivation of $\mathbf p = \gamma m\mathbf v$, Minkowski diagrams, and many worked problems.
- **MIT 8.07 lecture notes, lectures 25–30.** Radiation and relativistic electrodynamics at a slightly higher level, with the covariant formulation used throughout.
