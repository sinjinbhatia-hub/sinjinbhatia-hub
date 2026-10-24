# EM-5: Magnetostatics: Biot–Savart, Ampère, the Vector Potential, Magnetic Materials

**Time:** ~12 hours · **Prerequisites:** EM-4 (and through it EM-1's vector calculus and EM-3's multipole expansion) · **Unlocks:** EM-6 (Faraday, Maxwell), and later QM-4 (gauge invariance) and the flux-quantization hardware of QC

> **Conventions in this module.** SI units. $\mu_0 = 4\pi\times10^{-7}$ N/A² to better than one part in $10^9$ (§3.4). $\mathbf B$ in tesla, flux $\Phi = \int\mathbf B\cdot d\mathbf a$ in webers (T m²). Field point $\mathbf r$, source point $\mathbf r'$, separation $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$ (EM-1, EM-2). Currents: $I$ (A), surface density $\mathbf K$ (A/m), volume density $\mathbf J$ (A/m²); $\mathbf J_f$ is free current. **The right-hand rule is used identically everywhere:** $d\mathbf l'\times\hat{\boldsymbol{\mathscr r}}$ in Biot–Savart, $d\mathbf l$ against $d\mathbf a$ in Stokes, and $\mathbf m = I\mathbf a$. $\mathbf m$ is a dipole moment (A m²), $\mathbf M$ the magnetization (A/m, dipole moment per volume), $\mathbf H = \mathbf B/\mu_0 - \mathbf M$, and bound currents are $\mathbf J_b = \nabla\times\mathbf M$, $\mathbf K_b = \mathbf M\times\hat{\mathbf n}$ with $\hat{\mathbf n}$ outward. Lorentz force $\mathbf F = q(\mathbf E + \mathbf v\times\mathbf B)$. Vector potentials are in the **Coulomb gauge** $\nabla\cdot\mathbf A = 0$ unless stated. At a surface, $\hat{\mathbf n}$ points from "below" to "above" and jumps are "above minus below" (EM-2 §3.7). $\mathbf B$ is a pseudovector (EM-1 §3.2).

## 0. Where this sits

You need EM-1's product rules (P4)–(P6), $\nabla\times\nabla f = 0$, $\nabla\cdot(\nabla\times\mathbf A) = 0$ and the delta function in $\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$; EM-3's expansion of $1/\mathscr r$ in Legendre polynomials and the dipole formulas; and EM-4's bound-charge derivation, which this module copies line by line for bound currents. This module gives you the second half of Maxwell's equations in their static form, $\nabla\cdot\mathbf B = 0$ and $\nabla\times\mathbf B = \mu_0\mathbf J$, plus the vector potential that quantum mechanics cannot do without. EM-6 adds the two time-dependent terms (Faraday's law and the displacement current) and the energy in the field.

## 1. The problem that forced this idea

In April–July 1820 Hans Christian Ørsted noticed a compass needle turn near a current-carrying wire. The force was unlike any known one: Coulomb's and Newton's act along the line between bodies, but the needle was pushed *around* the wire. Within weeks André-Marie Ampère showed that parallel currents attract and opposed currents repel, wrote a force law between current elements, and proposed that a magnet is a collection of tiny circulating currents. That last guess is the bound current $\mathbf J_b = \nabla\times\mathbf M$ of §3.13, which makes it quantitative. Jean-Baptiste Biot and Félix Savart measured the force of a long straight wire on a magnet and found it fell as $1/s$; the law for each current element is what we call Biot–Savart.

Three differences from electrostatics followed. The force $q\mathbf v\times\mathbf B$ depends on velocity, is perpendicular to it, and so does no work (Theorem 3.1). The field lines never start or end, because no magnetic charge has ever been found. And the field circulates around its sources, so it needs its own integral theorem (Ampère's) and its own potential, a *vector* $\mathbf A$ with $\mathbf B = \nabla\times\mathbf A$. Faraday's induction (1831) and Maxwell's displacement current (1861) complete the picture in EM-6.

The vector potential has a strange history. Faraday's "electrotonic state" became Maxwell's vector potential, and for a century many physicists treated it as scaffolding (Heaviside and Hertz rewrote Maxwell's theory with $\mathbf E$ and $\mathbf B$ only). Aharonov and Bohm (1959; Ehrenberg and Siday in 1949) showed that an electron passing a solenoid where $\mathbf B = 0$ still has its interference pattern shifted, because quantum mechanics couples to $\mathbf A$. §3.8 computes the classical shadow of that effect, and it is why you cannot skip the vector potential.

## 2. Intuition

**A wire makes the field curl around it.** Point your right thumb along the current; your fingers curl the way compass needles point. Far from a long wire the field falls as $1/s$, so $B\cdot2\pi s$ is constant, and that constant is proportional to the current.

**Why the field never diverges.** A line that goes up through a current loop must return around the outside; there is no magnetic charge for it to end on. So the flux out of any closed surface is zero: $\nabla\cdot\mathbf B = 0$.

**Guess Ampère's law.** For a long wire $B\cdot2\pi s = \mu_0I$ at every radius, so the circulation around *any* loop depends only on the current through it. The magnetic twin of Gauss's law is "circulation $=\mu_0\times$ enclosed current", i.e. $\nabla\times\mathbf B = \mu_0\mathbf J$. In electrostatics $\nabla\cdot\mathbf E\ne0$ and $\nabla\times\mathbf E = 0$; here the two trade places, and the potential changes with them: a curl-free field is a gradient, but a divergence-free field is a curl, $\mathbf B = \nabla\times\mathbf A$, which needs a vector.

**Why a magnet has surface current.** Fill a block with identical tiny loops all pointing the same way. Inside, each loop's current is cancelled by its neighbour's on the shared edge; only the surface loops have an uncancelled edge, so the block acts like one sheet of current around its side. For a long rod this is a solenoid with $nI = M$, so $B\approx\mu_0M$ inside. Carry this guess into §3.13: $\mathbf K_b = \mathbf M\times\hat{\mathbf n}$, the twin of EM-4's $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$.

## 3. Formal development

### 3.1 The Lorentz force, cyclotron motion, and crossed fields

A charge $q$ with velocity $\mathbf v$ in fields $\mathbf E$ and $\mathbf B$ feels
$$\mathbf F = q(\mathbf E + \mathbf v\times\mathbf B). \tag{3.1}$$
This defines $\mathbf B$ operationally.

**Theorem 3.1 (Magnetic forces do no work).** In a pure magnetic field a charge's kinetic energy is constant.

*Proof.* The power is $\mathbf F_{\rm mag}\cdot\mathbf v = q(\mathbf v\times\mathbf B)\cdot\mathbf v = 0$. In general $dK/dt = q\mathbf E\cdot\mathbf v$. ∎

A magnet lifting a paperclip does not violate this: the energy comes from whatever moved the magnet or sustains the currents (EM-6).

**Cyclotron motion.** Take $\mathbf E = 0$ and uniform $\mathbf B$. The component of $\mathbf v$ along $\mathbf B$ is unforced. The perpendicular part feels a force of constant magnitude $qv_\perp B$, perpendicular to $\mathbf v_\perp$ (Theorem 3.1 keeps $v_\perp$ fixed): uniform circular motion with $mv_\perp^2/r = qv_\perp B$,
$$r = \frac{mv_\perp}{qB},\qquad \omega = \frac{v_\perp}{r} = \frac{qB}{m},\qquad T = \frac{2\pi m}{qB}. \tag{3.2}$$
The period does not depend on speed (non-relativistically), which is why Lawrence's cyclotron works: a fixed-frequency voltage kicks ions each half-turn as they spiral outward. A 1 MeV proton in 1 T has $v = 1.384\times10^7$ m/s, $r = 14.45$ cm, $f = 15.25$ MHz (project Part 1; the relativistic correction $\gamma - 1 = 10^{-3}$ is ignored).

**Mass spectrometer.** Accelerate an ion through $V$ ($qV = \tfrac12mv^2$) and send it perpendicular into $B$. By (3.2),
$$r = \frac{1}{B}\sqrt{\frac{2mV}{q}}. \tag{3.3}$$
For singly charged ions at 2 kV and 0.5 T, ${}^{12}$C$^+$ has $r = 4.461$ cm and ${}^{13}$C$^+$ has $4.644$ cm: the landing points ($2r$ from the slit) differ by $3.655$ mm.

▶ **Try it first.** Take $\mathbf E = E\hat{\mathbf y}$, $\mathbf B = B\hat{\mathbf z}$, charge released from rest. Write the two equations of motion and look for a combination that is a harmonic oscillator.

**Derivation (the cycloid).** Let $\omega = qB/m$.
1. $\dot v_x = \omega v_y$ and $\dot v_y = \omega(E/B - v_x)$, since $\mathbf v\times\mathbf B = (v_yB, -v_xB, 0)$.
2. With $u = v_x - E/B$: $\dot u = \omega v_y$ and $\dot v_y = -\omega u$, a rotation in the $(u, v_y)$ plane at rate $\omega$.
3. From rest, $u = -(E/B)\cos\omega t$ and $v_y = (E/B)\sin\omega t$, so $v_x = (E/B)(1 - \cos\omega t)$.
4. Integrating from the origin,
$$x = \frac{E}{\omega B}(\omega t - \sin\omega t),\qquad y = \frac{E}{\omega B}(1 - \cos\omega t). \tag{3.4}$$
∎

■ **What just happened.** Step 2 subtracted the **drift velocity** $\mathbf v_d = \mathbf E\times\mathbf B/B^2$, independent of $q$ and $m$. In a frame moving at $\mathbf v_d$ the electric field vanishes (non-relativistically) and the motion is a pure circle; a circle rolling along at $v_d$ traces a cycloid, with arch height $2mE/(qB^2)$. For a proton with $E = 1$ kV/m, $B = 10$ mT: $v_d = 10^5$ m/s and the arch is $208.8$ mm (Your turn 1 checks (3.4) against an integrator).

### 3.2 Currents and the continuity equation

Charge density $\rho$ moving with velocity $\mathbf v$ gives volume current density; the same idea on a surface and on a line gives
$$\mathbf J = \rho\mathbf v\ \ (\text{A/m}^2),\qquad \mathbf K = \sigma\mathbf v\ \ (\text{A/m}),\qquad I = \lambda v\ \ (\text{A}). \tag{3.5}$$
$\mathbf K$ is current per unit *width* perpendicular to the flow; $\mathbf J$ is current per unit *area* perpendicular to the flow. In each case the element $I\,d\mathbf l$, $\mathbf K\,da$ or $\mathbf J\,d\tau$ is a charge $dq$ times its velocity, so summing (3.1)'s magnetic part over charges gives
$$\mathbf F = \int I\,d\mathbf l\times\mathbf B = \int\mathbf K\times\mathbf B\,da = \int\mathbf J\times\mathbf B\,d\tau. \tag{3.6}$$

**Theorem 3.2 (Continuity).** $\displaystyle\frac{\partial\rho}{\partial t} = -\nabla\cdot\mathbf J.$ (3.7)

*Proof.* The charge in a fixed volume $\mathcal V$ with boundary $\mathcal S$ changes only by flow through $\mathcal S$: $\frac{d}{dt}\int_{\mathcal V}\rho\,d\tau = -\oint_{\mathcal S}\mathbf J\cdot d\mathbf a$. The divergence theorem (EM-1 (3.18)) turns the right side into $-\int_{\mathcal V}\nabla\cdot\mathbf J\,d\tau$, so $\int_{\mathcal V}(\partial_t\rho + \nabla\cdot\mathbf J)\,d\tau = 0$ for every $\mathcal V$, which forces the integrand to vanish. ∎

**Steady currents** have $\partial_t\rho = 0$, hence
$$\nabla\cdot\mathbf J = 0. \tag{3.8}$$
Current lines have no beginning or end: they close. That is Kirchhoff's current law at a node, and the reason wires must form loops. Magnetostatics means steady currents. A charging capacitor violates (3.8) in the gap, and EM-6 will show that Ampère's law fails there until Maxwell fixes it.

### 3.3 The Biot–Savart law

The experimental input of magnetostatics is the field of a steady current, summed over its elements:
$$\mathbf B(\mathbf r) = \frac{\mu_0}{4\pi}\int\frac{\mathbf J(\mathbf r')\times\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\,d\tau' \ \longrightarrow\ \frac{\mu_0I}{4\pi}\oint\frac{d\mathbf l'\times\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\ \text{(wire)},\quad \frac{\mu_0}{4\pi}\int\frac{\mathbf K\times\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\,da'\ \text{(sheet)}. \tag{3.9}$$
Compare Coulomb (EM-2 (3.3)): the same $1/\mathscr r^2$ and the same superposition, but the source $\rho$ is replaced by a *vector* $\mathbf J$ and the direction $\hat{\boldsymbol{\mathscr r}}$ by a cross product with it. $d\mathbf l'\times\hat{\boldsymbol{\mathscr r}}$ is the right-hand rule: fingers along the current, curl toward the field point, thumb along $\mathbf B$. A single current element has no meaning on its own (steady currents close, (3.8)); the law is a statement about complete circuits, and for one moving point charge it is valid only for $v\ll c$.

### 3.4 Four fields from Biot–Savart

**The straight wire.** A wire along the $z$-axis from $z_1$ to $z_2$ carries $I$. Find $\mathbf B$ at distance $s$ from the axis, at $z = 0$.

▶ **Try it first.** Write $d\mathbf l'$ and $\boldsymbol{\mathscr r}$ in cylindrical components and do the cross product before integrating.

1. $d\mathbf l' = dz'\,\hat{\mathbf z}$ and $\boldsymbol{\mathscr r} = s\hat{\mathbf s} - z'\hat{\mathbf z}$, so $\mathscr r = \sqrt{s^2 + z'^2}$ and $d\mathbf l'\times\boldsymbol{\mathscr r} = s\,dz'\,\hat{\boldsymbol\phi}$.
2. $\mathbf B = \dfrac{\mu_0I}{4\pi}\hat{\boldsymbol\phi}\displaystyle\int_{z_1}^{z_2}\frac{s\,dz'}{(s^2 + z'^2)^{3/2}} = \frac{\mu_0I}{4\pi s}\hat{\boldsymbol\phi}\left[\frac{z'}{\sqrt{s^2 + z'^2}}\right]_{z_1}^{z_2}$.
3. With $\sin\theta = z'/\sqrt{s^2 + z'^2}$ (the angle from the perpendicular to a wire end),
$$\mathbf B = \frac{\mu_0I}{4\pi s}(\sin\theta_2 - \sin\theta_1)\,\hat{\boldsymbol\phi}. \tag{3.10}$$
4. Infinite wire ($\theta_{1,2}\to\mp\pi/2$):
$$\mathbf B = \frac{\mu_0I}{2\pi s}\hat{\boldsymbol\phi}. \tag{3.11}$$
∎

For a segment in any position, let $\hat{\boldsymbol\ell}$ be its direction, $L$ its length, $u = (\mathbf r - \mathbf a)\cdot\hat{\boldsymbol\ell}$ the offset of the foot of the perpendicular from the start $\mathbf a$, and $d$ the perpendicular distance. Step 2 is the same integral:
$$\mathbf B = \frac{\mu_0I}{4\pi}\frac{\hat{\boldsymbol\ell}\times(\mathbf r - \mathbf a)}{d^2}\left[\frac{L - u}{\sqrt{d^2 + (L-u)^2}} + \frac{u}{\sqrt{d^2 + u^2}}\right]. \tag{3.12}$$
This exact segment field is what the project's integrator sums over the sides of a polygon.

**The circular loop on axis.** Ring of radius $R$ in the plane $z = 0$, current counterclockwise from $+z$. At $(0,0,z)$ every element has $d\mathbf l'\perp\boldsymbol{\mathscr r}$ and $\mathscr r^2 = R^2 + z^2$, so each contributes $\mu_0I\,dl'/4\pi\mathscr r^2$ perpendicular to $\boldsymbol{\mathscr r}$ and $d\mathbf l'$. Horizontal parts cancel in opposite pairs; the vertical part carries a factor $R/\mathscr r$:
$$B_z = \frac{\mu_0I}{4\pi}\frac{R}{\mathscr r^3}\oint dl' = \frac{\mu_0IR^2}{2(R^2 + z^2)^{3/2}}. \tag{3.13}$$
At the centre this is $\mu_0I/2R$; far away it is $\mu_0IR^2/2z^3 = \mu_0m/2\pi z^3$ with $m = I\pi R^2$ (§3.10).

**The solenoid on axis.** A tube of radius $R$ and length $L$ with $n$ turns per length is a stack of rings; the slice $dz'$ carries $nI\,dz'$. With $w = z - z'$ and $\int R^2dw/(R^2 + w^2)^{3/2} = w/\sqrt{R^2 + w^2}$, (3.13) integrates to
$$B_z(z) = \frac{\mu_0nI}{2}\left[\frac{z + L/2}{\sqrt{R^2 + (z + L/2)^2}} - \frac{z - L/2}{\sqrt{R^2 + (z - L/2)^2}}\right]\xrightarrow{L\to\infty}\mu_0nI. \tag{3.14}$$
(The bracket is $\cos\theta_2 - \cos\theta_1$ for the angles the ends subtend.) For $n = 4000$ /m, $I = 3$ A, $R = 2$ cm, $L = 0.5$ m: $\mu_0nI = 15.080$ mT and the centre field is $15.032$ mT (Example 4.2).

**Helmholtz coils.** Two coaxial rings, radius $R$, current $I$, at $z = \pm d/2$. With $f(u) = (R^2 + u^2)^{-3/2}$, (3.13) gives $B_z = \tfrac12\mu_0IR^2[f(z - d/2) + f(z + d/2)]$. The pair is even in $z$, so odd derivatives vanish at the centre. Since $f'' = 3(4u^2 - R^2)(R^2 + u^2)^{-7/2}$ and both terms equal $3(d^2 - R^2)/(R^2 + d^2/4)^{7/2}$ at $u = \pm d/2$,
$$B_z''(0) = \frac{3\mu_0IR^2(d^2 - R^2)}{(R^2 + d^2/4)^{7/2}}, \tag{3.15}$$
which vanishes exactly at $d = R$. Then
$$B_z = B_0\left[1 - \tfrac{144}{125}(z/R)^4 + \cdots\right],\qquad B_0 = (4/5)^{3/2}\mu_0I/R, \tag{3.16}$$
with the coefficient derived symbolically in the project. With *opposite* currents and $d = \sqrt3R$ the third derivative cancels instead, giving a clean linear gradient (Exercise B3).

**Parallel wires.** Wire 1 (current $I_1$ along $+z$) makes $B = \mu_0I_1/2\pi d$ circling it, perpendicular to wire 2 at distance $d$. By (3.6) the force per length on wire 2 (current $I_2$ along $+z$) has magnitude $I_2B$ and points toward wire 1:
$$\frac{F}{L} = \frac{\mu_0I_1I_2}{2\pi d}\quad(\text{attractive for parallel currents}). \tag{3.17}$$
For $1$ A each at $1$ m this is $2\times10^{-7}$ N/m. Until 2019 that *defined* the ampere and fixed $\mu_0 = 4\pi\times10^{-7}$ N/A² exactly. Since 20 May 2019 the ampere is defined by fixing the elementary charge $e$; $\mu_0$ is measured and equals $4\pi\times10^{-7}$ to better than $10^{-9}$ (the project uses SciPy's value).

### 3.5 The divergence and curl of $\mathbf B$

This is the longest vector-calculus computation in the course. Do it once, on paper, with EM-1's product-rule table open.

**Theorem 3.3.** For a steady, localized current density ($\nabla\cdot\mathbf J = 0$, $\mathbf J = 0$ outside a bounded region), the Biot–Savart field (3.9) satisfies
$$\nabla\cdot\mathbf B = 0,\qquad\nabla\times\mathbf B = \mu_0\mathbf J. \tag{3.18}$$

▶ **Try it first.** Write $\mathbf F = \hat{\boldsymbol{\mathscr r}}/\mathscr r^2$. In $\nabla\cdot(\mathbf J\times\mathbf F)$ and $\nabla\times(\mathbf J\times\mathbf F)$, which derivatives hit $\mathbf J$? (None: $\nabla$ differentiates $\mathbf r$, and $\mathbf J$ depends on $\mathbf r'$.) Then pick the right line of EM-1 §3.5.

**Derivation.**
1. $\nabla$ acts on the field point $\mathbf r$ and passes inside the integral: $\mathbf B = \frac{\mu_0}{4\pi}\int\mathbf J(\mathbf r')\times\mathbf F\,d\tau'$ with $\mathbf F = \hat{\boldsymbol{\mathscr r}}/\mathscr r^2 = -\nabla(1/\mathscr r)$.
2. *Divergence.* By (P4), $\nabla\cdot(\mathbf J\times\mathbf F) = \mathbf F\cdot(\nabla\times\mathbf J) - \mathbf J\cdot(\nabla\times\mathbf F)$. Here $\nabla\times\mathbf J = 0$ (no $\mathbf r$ dependence) and $\nabla\times\mathbf F = -\nabla\times\nabla(1/\mathscr r) = 0$ by EM-1 (3.11). So $\nabla\cdot\mathbf B = 0$. ∎
3. *Curl.* By (P6) with $\mathbf J$ constant under $\nabla$: $\nabla\times(\mathbf J\times\mathbf F) = \mathbf J(\nabla\cdot\mathbf F) - (\mathbf J\cdot\nabla)\mathbf F$.
4. *The delta function.* By EM-1 (3.32), $\nabla\cdot\mathbf F = 4\pi\delta^3(\boldsymbol{\mathscr r})$, so the first term integrates to $\frac{\mu_0}{4\pi}\int\mathbf J(\mathbf r')\,4\pi\delta^3(\mathbf r - \mathbf r')\,d\tau' = \mu_0\mathbf J(\mathbf r)$.
5. *The surface term.* For a function of $\boldsymbol{\mathscr r}$, $\nabla = -\nabla'$. Component $i$ of the second term: $[(\mathbf J\cdot\nabla)\mathbf F]_i = -J_j\partial'_jF_i = -\partial'_j(J_jF_i) + F_i\,\partial'_jJ_j$. The last piece is $F_i\,\nabla'\cdot\mathbf J = 0$ by (3.8). The first integrates by the divergence theorem to $-\oint F_i\,\mathbf J\cdot d\mathbf a'$ over a surface enclosing all the current, where $\mathbf J = 0$. So $\int(\mathbf J\cdot\nabla)\mathbf F\,d\tau' = 0$.
6. Combining 3–5: $\nabla\times\mathbf B = \mu_0\mathbf J$. ∎

■ **What just happened.** The divergence vanishes because Biot–Savart's kernel is a gradient. The curl gives $\mu_0\mathbf J$ because the kernel's divergence is a delta function, the same $4\pi\delta^3$ that turned Coulomb into Gauss in EM-2. Step 5 needs $\nabla\cdot\mathbf J = 0$; without it a second term survives, which is Maxwell's problem with the charging capacitor (EM-6). One subtlety: $\partial_jF_i$ is a distribution, $\mathrm{PV}\,[(\delta_{ij} - 3\hat{\mathscr r}_i\hat{\mathscr r}_j)/\mathscr r^3] + \tfrac{4\pi}{3}\delta_{ij}\delta^3(\boldsymbol{\mathscr r})$, whose trace is step 4's $4\pi\delta^3$. Step 5 integrates by parts, which *defines* the distributional derivative (EM-1 §3.11), so no boundary term appears at $\boldsymbol{\mathscr r} = 0$. The same $\tfrac{4\pi}{3}$ gives the $\tfrac23\mu_0\mathbf m\,\delta^3$ in the dipole field (3.33). Project Part 3 checks both equations of (3.18) on grids.
### 3.6 Ampère's law

**Theorem 3.4 (Ampère's law).** For any closed path $\mathcal P$ and any surface $\mathcal S$ it bounds, oriented by the right-hand rule,
$$\oint_{\mathcal P}\mathbf B\cdot d\mathbf l = \mu_0I_{\rm enc},\qquad I_{\rm enc} = \int_{\mathcal S}\mathbf J\cdot d\mathbf a. \tag{3.19}$$

*Proof.* Stokes (EM-1 (3.19)) and (3.18). $I_{\rm enc}$ does not depend on the surface: two surfaces with the same boundary form a closed surface, through which the net current is $\int\nabla\cdot\mathbf J\,d\tau = 0$ by (3.8). ∎

Ampère's law is always true for steady currents. It *determines* $\mathbf B$ only when symmetry fixes the direction of $\mathbf B$ and makes $\mathbf B\cdot d\mathbf l$ constant along some path. That needs enough continuous translation or rotation symmetry that $\mathbf B$ depends on one coordinate: an infinite wire, an infinite sheet, an infinite solenoid, a toroid. A finite wire, a loop or a finite solenoid lacks it; the law still holds but does not fix $\mathbf B$, and you use Biot–Savart.

**Infinite wire.** By (3.9), $\mathbf B\propto d\mathbf l'\times\boldsymbol{\mathscr r}$ is perpendicular to the plane containing the wire and the field point, so $\mathbf B = B_\phi(s)\hat{\boldsymbol\phi}$. A circle of radius $s$: $2\pi sB_\phi = \mu_0I$, which is (3.11).

**Infinite sheet** $\mathbf K = K\hat{\mathbf x}$ at $z = 0$. Pair the strips at $\pm y'$: their $z$-components cancel (mirror images about $y = 0$), and the $x$-component vanishes because $\hat{\mathbf x}\times\boldsymbol{\mathscr r}$ has none. So $\mathbf B = B(z)\hat{\mathbf y}$, odd in $z$. A rectangle in the $yz$-plane with long sides (length $\ell$, along $\hat{\mathbf y}$) at $\pm z$ gives $2B\ell = \mu_0K\ell$:
$$\mathbf B = -\tfrac12\mu_0K\,\mathrm{sgn}(z)\,\hat{\mathbf y}. \tag{3.20}$$
(Above the sheet the field is along $\hat{\mathbf x}\times\hat{\mathbf z} = -\hat{\mathbf y}$.)

**Infinite solenoid** ($n$ turns per length, a cylindrical sheet with $K = nI$). $B_s = 0$, else a coaxial cylinder would have net flux $2\pi sLB_s\ne0$. $B_\phi = 0$: a coaxial circle encloses no axial current. So $\mathbf B = B(s)\hat{\mathbf z}$. A rectangle with both long sides outside encloses no current, so the outside field is constant, hence zero (it vanishes at infinity). A rectangle with one long side inside encloses $nI\ell$:
$$\mathbf B = \mu_0nI\,\hat{\mathbf z}\ \ (s < R),\qquad\mathbf B = 0\ \ (s > R). \tag{3.21}$$
**Toroid** ($N$ turns). Rotation symmetry and dense winding give $\mathbf B = B(s)\hat{\boldsymbol\phi}$. A circle inside the doughnut encloses $N$ wires of current $I$; one in the hole or outside encloses none (outside, the current goes in and comes out):
$$B = \frac{\mu_0NI}{2\pi s}\ \ (\text{inside the doughnut}),\qquad B = 0\ \ (\text{elsewhere}). \tag{3.22}$$
Unlike the solenoid, this varies across the section ("$\mu_0nI$ with $n = N/2\pi s_{\rm mean}$" is the thin-toroid approximation; Your turn 2).

| | Electrostatics | Magnetostatics |
|---|---|---|
| Source | charge $\rho$ (scalar) | current $\mathbf J$ (vector), $\nabla\cdot\mathbf J = 0$ |
| Force on a charge | $q\mathbf E$ | $q\mathbf v\times\mathbf B$, does no work |
| Field from sources | Coulomb: $\rho\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$ | Biot–Savart: $\mathbf J\times\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$ |
| Divergence, curl | $\nabla\cdot\mathbf E = \rho/\varepsilon_0$, $\nabla\times\mathbf E = 0$ | $\nabla\cdot\mathbf B = 0$, $\nabla\times\mathbf B = \mu_0\mathbf J$ |
| Integral law | Gauss (flux) | Ampère (circulation) |
| Symmetries that make it usable | spherical, cylindrical, planar | infinite line, plane, solenoid; toroid |
| Potential | $\mathbf E = -\nabla V$, $\nabla^2V = -\rho/\varepsilon_0$ | $\mathbf B = \nabla\times\mathbf A$, $\nabla^2\mathbf A = -\mu_0\mathbf J$ (Coulomb gauge) |
| Surface charge or current | $E^\perp$ jumps by $\sigma/\varepsilon_0$; $E^\parallel$ continuous | $B^\perp$ continuous; $B^\parallel$ jumps by $\mu_0K$ |
| Dipole | $\mathbf p$, $V = \mathbf p\cdot\hat{\mathbf r}/4\pi\varepsilon_0r^2$ | $\mathbf m$, $\mathbf A = \mu_0\mathbf m\times\hat{\mathbf r}/4\pi r^2$ |
| In matter | $\mathbf P$, $\rho_b = -\nabla\cdot\mathbf P$, $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$, $\mathbf D$ | $\mathbf M$, $\mathbf J_b = \nabla\times\mathbf M$, $\mathbf K_b = \mathbf M\times\hat{\mathbf n}$, $\mathbf H$ |
| Uniformly polarized or magnetized sphere | $\mathbf E_{\rm in} = -\mathbf P/3\varepsilon_0$ | $\mathbf B_{\rm in} = \tfrac23\mu_0\mathbf M$, $\mathbf H_{\rm in} = -\mathbf M/3$ |

### 3.7 The vector potential

**Theorem 3.5 (Vector potential).** For steady localized currents there is a field $\mathbf A$ with
$$\mathbf B = \nabla\times\mathbf A,\qquad \mathbf A\to\mathbf A + \nabla\lambda\ \text{changes nothing},\qquad \nabla\cdot\mathbf A = 0\ \text{is reachable}, \tag{3.23}$$
and in this **Coulomb gauge**
$$\nabla^2\mathbf A = -\mu_0\mathbf J,\qquad \mathbf A(\mathbf r) = \frac{\mu_0}{4\pi}\int\frac{\mathbf J(\mathbf r')}{\mathscr r}\,d\tau'. \tag{3.24}$$

▶ **Try it first.** Take the curl of $\mathbf A = \frac{\mu_0}{4\pi}\int\mathbf J/\mathscr r$. Which product rule gives you Biot–Savart back?

**Derivation.**
1. *Existence.* $\nabla\cdot\mathbf B = 0$ on $\mathbb R^3$ implies $\mathbf B = \nabla\times\mathbf A$ (EM-1 Theorem 3.17, with the explicit $\mathbf W$ of (3.36) there).
2. *Gauge freedom.* $\nabla\times\nabla\lambda = 0$ (EM-1 (3.11)), so $\mathbf A' = \mathbf A + \nabla\lambda$ has the same curl.
3. *Reaching the Coulomb gauge.* Given any $\mathbf A_0$, choose $\lambda$ with $\nabla^2\lambda = -\nabla\cdot\mathbf A_0$. This is Poisson's equation with source $\nabla\cdot\mathbf A_0$, solved by $\lambda = \frac{1}{4\pi}\int(\nabla'\cdot\mathbf A_0)/\mathscr r\,d\tau'$ (EM-1 (3.33)). Then $\nabla\cdot(\mathbf A_0 + \nabla\lambda) = 0$.
4. *The equation.* $\mu_0\mathbf J = \nabla\times\mathbf B = \nabla\times(\nabla\times\mathbf A) = \nabla(\nabla\cdot\mathbf A) - \nabla^2\mathbf A$ by EM-1 (3.13); in Coulomb gauge the first term vanishes. Each Cartesian component is Poisson's equation, whose solution vanishing at infinity is (EM-1 (3.33)) $\mathbf A = \frac{\mu_0}{4\pi}\int\mathbf J/\mathscr r\,d\tau'$.
5. *Consistency.* $\nabla\cdot\mathbf A = \frac{\mu_0}{4\pi}\int\mathbf J\cdot\nabla\frac1{\mathscr r}\,d\tau' = -\frac{\mu_0}{4\pi}\int\mathbf J\cdot\nabla'\frac{1}{\mathscr r}\,d\tau'$. By (P3) and the divergence theorem this is $\frac{\mu_0}{4\pi}\int\frac{\nabla'\cdot\mathbf J}{\mathscr r}\,d\tau' - \frac{\mu_0}{4\pi}\oint\frac{\mathbf J\cdot d\mathbf a'}{\mathscr r} = 0$: zero by (3.8) and because $\mathbf J = 0$ on the outer surface. And by (P5), $\nabla\times(\mathbf J/\mathscr r) = \nabla\frac1{\mathscr r}\times\mathbf J = \mathbf J\times\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$, which returns (3.9). ∎

■ **What just happened.** $\mathbf A$ is to $\mathbf B$ what $V$ is to $\mathbf E$ with div and curl swapped, and (3.24) is EM-2's $V = \frac{1}{4\pi\varepsilon_0}\int\rho/\mathscr r$ with $\rho/\varepsilon_0\to\mu_0\mathbf J$, component by component. $\mathbf A$ points *along* the current where $\mathbf B$ circles it. For infinite currents the integral diverges as in EM-2: an infinite wire has $A_z = -(\mu_0I/2\pi)\ln s + $ const, the log potential of EM-2 Example 4.2.

### 3.8 Worked potentials: the solenoid, the loop integral, the spinning shell

**Theorem 3.6 (The loop integral of $\mathbf A$ is the flux).** For a closed loop $\mathcal P$ bounding $\mathcal S$,
$$\oint_{\mathcal P}\mathbf A\cdot d\mathbf l = \int_{\mathcal S}(\nabla\times\mathbf A)\cdot d\mathbf a = \Phi_{\mathcal S}. \tag{3.25}$$
*Proof.* Stokes. ∎ It is gauge invariant, since $\oint\nabla\lambda\cdot d\mathbf l = 0$ for single-valued $\lambda$ (EM-1 (3.20)). So $\mathbf A$ itself changes with gauge and is not classically measurable, but its loop integral is: it is the flux. That is what survives in quantum hardware: flux quantization in a superconducting loop ($\Phi = n\,h/2e$) and the flux through a SQUID or flux-tunable qubit are statements about $\oint\mathbf A\cdot d\mathbf l$ (QM-4 for the phase).

**The infinite solenoid.** Radius $R$, $\mathbf B = \mu_0nI\hat{\mathbf z}$ inside and $0$ outside. By symmetry $\mathbf A = A_\phi(s)\hat{\boldsymbol\phi}$, and (3.25) on a circle of radius $s$ gives $2\pi sA_\phi = \Phi_{\rm enc}$, which is $\mu_0nI\pi s^2$ for $s<R$ and $\mu_0nI\pi R^2$ for $s>R$:
$$A_\phi = \frac{\mu_0nI}{2}s\ \ (s < R),\qquad A_\phi = \frac{\mu_0nIR^2}{2s}\ \ (s > R). \tag{3.26}$$
Checks: $\nabla\cdot\mathbf A = s^{-1}\partial_\phi A_\phi = 0$; $(\nabla\times\mathbf A)_z = s^{-1}\partial_s(sA_\phi)$ is $\mu_0nI$ inside and $0$ outside; $\mathbf A\to0$ at infinity; $A_\phi$ is continuous at $s = R$ and its slope jumps by $-\mu_0nI = -\mu_0K$, the boundary condition (3.29).

**Where $\mathbf B = 0$, $\mathbf A\ne0$.** For the example solenoid ($R = 2$ cm, $\mu_0nI = 15.08$ mT), $A_\phi(1.5R) = 100.5\ \mu$T m and $\oint\mathbf A\cdot d\mathbf l = 1.895\times10^{-5}$ Wb however large the circle. This is the Aharonov–Bohm setup: an electron's wavefunction acquires the phase $(q/\hbar)\int\mathbf A\cdot d\mathbf l$ along each path, so two paths on opposite sides of the solenoid differ by $q\Phi/\hbar$ with no field on either. The fringes shift by one each time $\Phi$ grows by $h/e = 4.136\times10^{-15}$ Wb. The project verifies for the *finite* solenoid that the numerical $\oint\mathbf A\cdot d\mathbf l$ equals the flux of the Biot–Savart $\mathbf B$ (to about $5\times10^{-4}$, at $1.5R$ and $3R$) and that a gauge change leaves it unchanged.

**The spinning charged shell.** A sphere of radius $R$ with uniform $\sigma$ rotates at $\boldsymbol\omega$: $\mathbf K = \sigma\boldsymbol\omega\times\mathbf r'$ by (3.5).

▶ **Try it first.** By symmetry $\mathbf A$ is $\boldsymbol\omega\times\hat{\mathbf r}$ times a function of $r$. Put the field point on the $z$-axis, $\boldsymbol\omega$ in the $xz$-plane at angle $\psi$ to it, and see which components survive the $\phi'$ integral.

**Derivation.**
1. $\mathbf A = \frac{\mu_0}{4\pi}\int\mathbf K/\mathscr r\,da'$ with $\mathbf r = r\hat{\mathbf z}$, so $\mathscr r^2 = r^2 + R^2 - 2rR\cos\theta'$ depends on $\theta'$ only.
2. With $\boldsymbol\omega = \omega(\sin\psi\,\hat{\mathbf x} + \cos\psi\,\hat{\mathbf z})$ and $\mathbf r' = R(\sin\theta'\cos\phi', \sin\theta'\sin\phi', \cos\theta')$, the $\hat{\mathbf y}$-component of $\boldsymbol\omega\times\mathbf r'$ is $\omega R(\cos\psi\sin\theta'\cos\phi' - \sin\psi\cos\theta')$; the $\hat{\mathbf x}$ and $\hat{\mathbf z}$ components are proportional to $\sin\phi'$.
3. $\int_0^{2\pi}\sin\phi'\,d\phi' = \int_0^{2\pi}\cos\phi'\,d\phi' = 0$ leaves only $-\sin\psi\cos\theta'$: $\mathbf A = -\tfrac12\mu_0\sigma\omega R^3\sin\psi\,\hat{\mathbf y}\int_{-1}^1u\,du/\mathscr r$, $u = \cos\theta'$.
4. EM-3 (3.27): $1/\mathscr r = \sum_\ell(r_<^\ell/r_>^{\ell+1})P_\ell(u)$. With $u = P_1$ and $\int_{-1}^1P_1P_\ell\,du = \tfrac23\delta_{\ell1}$, the integral is $\tfrac23r_</r_>^2$.
5. $\boldsymbol\omega\times\mathbf r = -\omega r\sin\psi\,\hat{\mathbf y}$, so the result is rotation covariant:
$$\mathbf A = \frac{\mu_0R\sigma}{3}\,\boldsymbol\omega\times\mathbf r\ \ (r < R),\qquad \mathbf A = \frac{\mu_0R^4\sigma}{3}\,\frac{\boldsymbol\omega\times\mathbf r}{r^3}\ \ (r > R). \tag{3.27}$$
6. Inside, $\nabla\times(\boldsymbol\omega\times\mathbf r) = 2\boldsymbol\omega$ (EM-1 (P6)), so $\mathbf B_{\rm in} = \tfrac23\mu_0R\sigma\boldsymbol\omega$, **uniform**. Outside, (3.27) is $\frac{\mu_0}{4\pi}\mathbf m\times\hat{\mathbf r}/r^2$ with $\mathbf m = \tfrac{4\pi}{3}R^4\sigma\boldsymbol\omega$, a pure dipole (§3.10). ∎

■ **What just happened.** The Legendre expansion reduced an ugly integral to one orthogonality relation. The field is uniform inside and dipolar outside, the structure of EM-3 Example 4.5, because $\mathbf K\propto\sin\theta\,\hat{\boldsymbol\phi}$ is a $\cos\theta$-type source. Example 4.3 gives numbers; project Part 5 builds the shell from 200 current rings and reproduces (3.27).

### 3.9 Boundary conditions at a surface current

Let a sheet carry $\mathbf K$, with $\hat{\mathbf n}$ from below to above.

**Theorem 3.7.**
$$B^\perp_{\rm above} = B^\perp_{\rm below},\qquad \mathbf B^\parallel_{\rm above} - \mathbf B^\parallel_{\rm below} = \mu_0\,\mathbf K\times\hat{\mathbf n}, \tag{3.28}$$
$$\mathbf A_{\rm above} = \mathbf A_{\rm below},\qquad \frac{\partial\mathbf A_{\rm above}}{\partial n} - \frac{\partial\mathbf A_{\rm below}}{\partial n} = -\mu_0\mathbf K. \tag{3.29}$$

*Proof.* (i) A flat pillbox straddling the sheet: $\nabla\cdot\mathbf B = 0$ gives $(B^\perp_{\rm above} - B^\perp_{\rm below})a = 0$ as the thickness $\to0$. (ii) A thin loop with long sides of length $\ell$ along a unit vector $\hat{\mathbf t}\perp\hat{\mathbf n}$, traversed along $+\hat{\mathbf t}$ above, has right-hand-rule normal $\hat{\mathbf n}\times\hat{\mathbf t}$. Its circulation is $(\mathbf B_{\rm above} - \mathbf B_{\rm below})\cdot\hat{\mathbf t}\,\ell$ and the enclosed current is $\ell\,\mathbf K\cdot(\hat{\mathbf n}\times\hat{\mathbf t}) = \ell\,\hat{\mathbf t}\cdot(\mathbf K\times\hat{\mathbf n})$; equate by (3.19) for every $\hat{\mathbf t}$. (iii) The flux through a thin loop straddling the sheet tends to zero, so (3.25) makes tangential $\mathbf A$ continuous; a pillbox with $\nabla\cdot\mathbf A = 0$ does the same for the normal part. (iv) With $\mathbf A$ continuous only normal derivatives can jump, and $\nabla\times\mathbf A$ jumps by $\hat{\mathbf n}\times[\partial_n\mathbf A]$. Setting this equal to $\mu_0\mathbf K\times\hat{\mathbf n} = -\mu_0\hat{\mathbf n}\times\mathbf K$ gives $[\partial_n\mathbf A] = -\mu_0\mathbf K$ (a normal component is excluded by $\nabla\cdot\mathbf A = 0$). ∎

Checks: for the plane (3.20), $\mathbf K\times\hat{\mathbf n} = -K\hat{\mathbf y}$ and (3.20) jumps by $-\mu_0K\hat{\mathbf y}$. At the solenoid wall $B_{\rm in} - B_{\rm out} = \mu_0nI$, confirmed on the computed finite solenoid ($15.078$ vs $15.080$ mT).

### 3.10 The multipole expansion and the magnetic dipole

Expand $1/\mathscr r$ for $r > r'$ with EM-3 (3.27) in $\mathbf A = \frac{\mu_0I}{4\pi}\oint d\mathbf l'/\mathscr r$:
$$\mathbf A = \frac{\mu_0I}{4\pi}\sum_{\ell=0}^\infty\frac{1}{r^{\ell+1}}\oint(r')^\ell P_\ell(\cos\alpha)\,d\mathbf l'. \tag{3.30}$$

**The monopole term vanishes.** $\ell = 0$ gives $\oint d\mathbf l' = 0$ for a closed loop. For any steady current, $\int\mathbf J\,d\tau' = 0$: component $i$ is $\int\mathbf J\cdot\nabla'x'_i\,d\tau' = \oint x'_i\mathbf J\cdot d\mathbf a' - \int x'_i\nabla'\cdot\mathbf J\,d\tau' = 0$. This is $\nabla\cdot\mathbf B = 0$ again: no magnetic charge supplies a $1/r$ term. **The dipole is the leading term.**

**The dipole term.** For $\ell = 1$, $P_1(\cos\alpha) = \hat{\mathbf r}\cdot\mathbf r'/r'$, so $\mathbf A_{\rm dip} = \frac{\mu_0I}{4\pi r^2}\oint(\hat{\mathbf r}\cdot\mathbf r')\,d\mathbf l'$.

▶ **Try it first.** Dot with a constant vector $\mathbf c$ and apply Stokes to the field $(\hat{\mathbf r}\cdot\mathbf r')\mathbf c$.

1. $\mathbf c\cdot\oint(\hat{\mathbf r}\cdot\mathbf r')d\mathbf l' = \oint\mathbf F\cdot d\mathbf l'$ with $\mathbf F = (\hat{\mathbf r}\cdot\mathbf r')\mathbf c$.
2. $\nabla'\times\mathbf F = \nabla'(\hat{\mathbf r}\cdot\mathbf r')\times\mathbf c = \hat{\mathbf r}\times\mathbf c$ by (P5).
3. Stokes: $\oint\mathbf F\cdot d\mathbf l' = (\hat{\mathbf r}\times\mathbf c)\cdot\mathbf a = \mathbf c\cdot(\mathbf a\times\hat{\mathbf r})$, where $\mathbf a = \int d\mathbf a'$ is the **vector area**. Since $\mathbf c$ is arbitrary,
$$\oint(\hat{\mathbf r}\cdot\mathbf r')\,d\mathbf l' = \mathbf a\times\hat{\mathbf r} = -\hat{\mathbf r}\times\mathbf a. \tag{3.31}$$
4. Define the **magnetic dipole moment** $\mathbf m = I\mathbf a$ ($I\times$ area along the right-hand normal for a planar loop):
$$\mathbf A_{\rm dip} = \frac{\mu_0}{4\pi}\frac{\mathbf m\times\hat{\mathbf r}}{r^2}. \tag{3.32}$$
∎

(The companion $\oint\mathbf r'(\hat{\mathbf r}\cdot d\mathbf l') = +\hat{\mathbf r}\times\mathbf a$ has the opposite sign, since the two integrands sum to the exact differential $d[(\hat{\mathbf r}\cdot\mathbf r')\mathbf r']$; this is a common sign error. For general currents $\mathbf m = \tfrac12\int\mathbf r'\times\mathbf J\,d\tau'$.)

**The dipole field.** For $\mathbf m = m\hat{\mathbf z}$, $A_\phi = \mu_0m\sin\theta/4\pi r^2$. The spherical curl (EM-1) gives $B_r = \frac{1}{r\sin\theta}\partial_\theta(\sin\theta A_\phi) = \frac{\mu_0m}{4\pi}\frac{2\cos\theta}{r^3}$ and $B_\theta = -\frac1r\partial_r(rA_\phi) = \frac{\mu_0m}{4\pi}\frac{\sin\theta}{r^3}$. In vector form, (P6) with $\mathbf m$ constant, $\nabla\cdot(\mathbf r/r^3) = 4\pi\delta^3$ and $(\mathbf m\cdot\nabla)(\mathbf r/r^3) = [\mathbf m - 3(\mathbf m\cdot\hat{\mathbf r})\hat{\mathbf r}]/r^3$ for $r\ne0$:
$$\mathbf B_{\rm dip} = \frac{\mu_0}{4\pi}\frac{3(\mathbf m\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf m}{r^3} + \frac{2\mu_0}{3}\,\mathbf m\,\delta^3(\mathbf r). \tag{3.33}$$
The delta term comes from the $\tfrac{4\pi}{3}\delta_{ij}\delta^3$ in $\partial_i(r_j/r^3)$ (§3.5); Exercise C1 checks it by integrating $\mathbf B$ over a sphere.

**Compare EM-3.** Away from the origin, (3.33) has exactly the form of $\mathbf E_{\rm dip} = \frac{1}{4\pi\varepsilon_0}[3(\mathbf p\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf p]/r^3$, with $\mathbf p/\varepsilon_0\leftrightarrow\mu_0\mathbf m$. The physics differs at the source. $\mathbf E$ lines run from $+$ to $-$ and *inside* the dipole $\mathbf E$ points against $\mathbf p$ (average $-\mathbf p/4\pi\varepsilon_0R^3$ over a sphere). $\mathbf B$ lines have no ends, thread the loop, and inside $\mathbf B$ points *along* $\mathbf m$ (the delta term has the opposite sign). Far from the loop you cannot tell a current loop from a pair of magnetic charges; at it you can.

### 3.11 A dipole in an external field

A loop carrying $I$ sits in a field $\mathbf B(\mathbf r)$ made elsewhere; by (3.6), $\mathbf F = I\oint d\mathbf l\times\mathbf B$.

▶ **Try it first.** Uniform $\mathbf B$: what is $\oint d\mathbf l$? For the torque, expand $\mathbf r\times(d\mathbf l\times\mathbf B)$.

**Derivation.**
1. *Uniform $\mathbf B$, force:* $\mathbf F = I(\oint d\mathbf l)\times\mathbf B = 0$.
2. *Torque:* $\mathbf r\times(d\mathbf l\times\mathbf B) = d\mathbf l\,(\mathbf r\cdot\mathbf B) - \mathbf B\,(\mathbf r\cdot d\mathbf l)$. The second integral is $\oint d(r^2/2) = 0$. The first is the identity in step 3 of §3.10, $\oint(\mathbf r\cdot\mathbf c)\,d\mathbf l = \mathbf a\times\mathbf c$ (with $\mathbf c = \mathbf B$), so
$$\mathbf N = I\,\mathbf a\times\mathbf B = \mathbf m\times\mathbf B. \tag{3.34}$$
3. *Force in a non-uniform field:* for a small loop expand $B_k\approx B_k(0) + x_l\partial_lB_k$; the constant part integrates to zero. Since $\oint d(x_lx_j) = 0$, $\oint x_l\,dl_j = \epsilon_{mlj}a_m$. Then $F_i = I\epsilon_{ijk}\partial_lB_k\,\epsilon_{mlj}a_m = Ia_m\partial_lB_k(\delta_{km}\delta_{il} - \delta_{kl}\delta_{im}) = \partial_i(\mathbf m\cdot\mathbf B) - m_i\nabla\cdot\mathbf B$, by EM-1 (3.1). With $\nabla\cdot\mathbf B = 0$,
$$\mathbf F = \nabla(\mathbf m\cdot\mathbf B). \tag{3.35}$$
∎

■ **What just happened.** (3.35) is $\mathbf F = -\nabla U$ with $U = -\mathbf m\cdot\mathbf B$, and (3.34) is $-\partial U/\partial\theta$ for the angle $\theta$ between $\mathbf m$ and $\mathbf B$: the dipole turns to align with $\mathbf B$ and is pulled toward stronger field, as a compass needle and an iron filing are, and these are EM-3 (3.32)–(3.34) with $\mathbf p\to\mathbf m$, $\mathbf E\to\mathbf B$. One caveat: magnetic forces do no work (Theorem 3.1), so $U$ is bookkeeping; the energy that appears as the loop moves comes from whatever holds its current constant, through Faraday's induced EMF (EM-6). Where $\nabla\times\mathbf B = 0$, (3.35) equals $(\mathbf m\cdot\nabla)\mathbf B$. A diamagnetic body has $\mathbf m$ induced *against* $\mathbf B$ and so is pushed toward weaker field.

### 3.12 Magnetism in matter: three behaviours and one honest warning

**Classical physics cannot produce any of them.** The **Bohr–van Leeuwen theorem** (Bohr 1911, van Leeuwen 1919): a classical system has $Z = \int e^{-\beta H(\mathbf p - q\mathbf A(\mathbf r),\,\mathbf r)}d^3p\,d^3r$; the shift $\mathbf p\to\mathbf p + q\mathbf A$ removes $\mathbf A$, so the free energy does not depend on $\mathbf B$ and $M = -\partial F/\partial B = 0$. Dia-, para- and ferromagnetism are all quantum (spin, orbital states, exclusion). The pictures below get signs and sizes right; they do not derive the effects.

**Diamagnetism** ($\chi_m<0$, present in all materials, weak): induced moments oppose $\mathbf B$. Orbit model: an electron (charge $-e$) circles at fixed radius $R$ and speed $v$, orbit normal along $\mathbf B$. The binding force is fixed; $B$ adds $\pm evB$ radially, so $m(v'^2 - v^2)/R\approx\pm evB$, $\Delta v = \pm eRB/2m$. The orbital moment has magnitude $evR/2$, so $|\Delta m| = e^2R^2B/4m$, and for either sense of circulation the change opposes $\mathbf B$:
$$\Delta\mathbf m = -\frac{e^2R^2}{4m}\mathbf B. \tag{3.36}$$
Averaging over orbit orientations ($\langle\rho^2\rangle = \tfrac23\langle r^2\rangle$) gives $\chi_m = -\mu_0nZe^2\langle r^2\rangle/6m_e$; for $n = 6\times10^{28}$ m$^{-3}$, $Z = 10$, $\langle r^2\rangle^{1/2} = 1$ Å this is $-3.5\times10^{-5}$, the observed order (table below). A superconductor is the extreme, $\chi_m = -1$.

**Paramagnetism** ($\chi_m>0$, weak): atoms with unpaired spins carry permanent moments $m_0$. A moment in $\mathbf B$ has energy $\mp m_0B$ (spin-½) and thermal agitation fights alignment; Boltzmann weights give SM-3's two-state paramagnet,
$$M = nm_0\tanh\frac{m_0B}{k_BT}\ \xrightarrow{\ m_0B\ll k_BT\ }\ \chi_m = \frac{\mu_0nm_0^2}{k_BT}\quad(\text{Curie's law}). \tag{3.37}$$
For one Bohr magneton per atom, $m_0B/k_BT = 2.2\times10^{-3}$ at 1 T and 300 K, and with $n = 6\times10^{28}$ m$^{-3}$, $\chi_m = 1.6\times10^{-3}$ (an upper end). Only a fraction of a percent of the moments align, which is why $\chi_m$ is small.

**Ferromagnetism** is cooperative alignment, the strong effect (§3.16).

### 3.13 Magnetization and bound currents

**Definition.** The **magnetization** $\mathbf M(\mathbf r)$ is the dipole moment per unit volume: a volume $d\tau'$ at $\mathbf r'$ carries moment $\mathbf M(\mathbf r')\,d\tau'$. Its vector potential is, from (3.32) with $\mathbf r\to\boldsymbol{\mathscr r}$,
$$\mathbf A(\mathbf r) = \frac{\mu_0}{4\pi}\int\frac{\mathbf M(\mathbf r')\times\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\,d\tau'. \tag{3.38}$$

**Theorem 3.8 (Bound currents).** $\mathbf A$ equals the vector potential of a volume current $\mathbf J_b$ and a surface current $\mathbf K_b$,
$$\mathbf A = \frac{\mu_0}{4\pi}\left[\int\frac{\mathbf J_b}{\mathscr r}\,d\tau' + \oint\frac{\mathbf K_b}{\mathscr r}\,da'\right],\qquad \mathbf J_b = \nabla\times\mathbf M,\quad\mathbf K_b = \mathbf M\times\hat{\mathbf n}. \tag{3.39}$$

▶ **Try it first.** This is EM-4's Theorem 3.3 with one integration by parts. Replace $\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$ by $\nabla'(1/\mathscr r)$ and look for the cross-product version of (P3).

**Derivation.**
1. By EM-4 (3.9), $\hat{\boldsymbol{\mathscr r}}/\mathscr r^2 = \nabla'(1/\mathscr r)$, so $\mathbf A = \frac{\mu_0}{4\pi}\int\mathbf M\times\nabla'\frac1{\mathscr r}\,d\tau'$.
2. (P5) with $f = 1/\mathscr r$: $\nabla'\times(\mathbf M/\mathscr r) = \frac1{\mathscr r}\nabla'\times\mathbf M - \mathbf M\times\nabla'\frac1{\mathscr r}$, so $\mathbf M\times\nabla'\frac1{\mathscr r} = \frac1{\mathscr r}\nabla'\times\mathbf M - \nabla'\times\frac{\mathbf M}{\mathscr r}$.
3. *Curl theorem:* $\int_{\mathcal V}\nabla\times\mathbf F\,d\tau = \oint_{\mathcal S}\hat{\mathbf n}\times\mathbf F\,da$. (Apply the divergence theorem to $\mathbf F\times\mathbf c$ for constant $\mathbf c$ and use (P4): $\oint(\mathbf F\times\mathbf c)\cdot d\mathbf a = \mathbf c\cdot\int\nabla\times\mathbf F\,d\tau$, and $(\mathbf F\times\mathbf c)\cdot\hat{\mathbf n} = \mathbf c\cdot(\hat{\mathbf n}\times\mathbf F)$.)
4. So $\mathbf A = \frac{\mu_0}{4\pi}\left[\int\frac{\nabla'\times\mathbf M}{\mathscr r}d\tau' - \oint\frac{\hat{\mathbf n}\times\mathbf M}{\mathscr r}da'\right]$, and $-\hat{\mathbf n}\times\mathbf M = \mathbf M\times\hat{\mathbf n}$. ∎

■ **What just happened.** A continuum of dipoles became a current distribution. Counting makes it concrete: take $\mathbf M = M_z(y)\hat{\mathbf z}$ and cubes of side $dx$, $dy$, $dz$ with loop current $I = M_zdz$ (so $m = M_zdx\,dy\,dz$). The edge shared by the cells at $y$ and $y + dy$ carries $+x$ current $[M_z(y + dy) - M_z(y)]dz$ net, i.e. $J_x = \partial_yM_z = (\nabla\times\mathbf M)_x$. This is Ampère's 1820 guess of molecular currents, made quantitative. Bound currents are real currents (they source $\mathbf B$) but not free ones (you cannot attach a battery). They close: $\nabla\cdot\mathbf J_b = \nabla\cdot(\nabla\times\mathbf M) = 0$.

**The uniformly magnetized sphere.** $\nabla\times\mathbf M = 0$ and $\mathbf K_b = M\hat{\mathbf z}\times\hat{\mathbf r} = M\sin\theta\,\hat{\boldsymbol\phi}$. That is the spinning shell with $\sigma R\omega = M$, so (3.27) gives
$$\mathbf B_{\rm in} = \tfrac23\mu_0\mathbf M,\qquad \mathbf B_{\rm out} = \text{dipole with }\mathbf m = \tfrac{4\pi}{3}R^3\mathbf M. \tag{3.40}$$
At the pole, $B_r$ is continuous: $\tfrac23\mu_0M$ inside and $2\mu_0m/4\pi R^3 = \tfrac23\mu_0M$ outside (consistent with (3.28)). The $\tfrac23$ is the twin of EM-4's $\tfrac13$: $\mathbf E_{\rm in} = -\mathbf P/3\varepsilon_0$ there, $\mathbf B_{\rm in} = +\tfrac23\mu_0\mathbf M$ here, the same relation as the delta terms of (3.33) and EM-3.

### 3.14 The auxiliary field $\mathbf H$

Total current is free plus bound: $\nabla\times\mathbf B = \mu_0(\mathbf J_f + \nabla\times\mathbf M)$. Collect the curls, $\nabla\times(\mathbf B/\mu_0 - \mathbf M) = \mathbf J_f$, and define
$$\mathbf H\equiv\frac{\mathbf B}{\mu_0} - \mathbf M,\qquad\nabla\times\mathbf H = \mathbf J_f,\qquad\oint\mathbf H\cdot d\mathbf l = I_{f,{\rm enc}}. \tag{3.41}$$
The right side contains only what you control. **What $\mathbf H$ is not:** $\mathbf B$ is the field: it fixes the force on charges (3.1), the flux, and $\mathbf A$. $\mathbf H$ is bookkeeping with a convenient curl, exactly as $\mathbf D$ was in EM-4. $\nabla\cdot\mathbf H = -\nabla\cdot\mathbf M\ne0$ in general, so $\mathbf H$ is *not* determined by $\mathbf J_f$ alone: a bar magnet has $\mathbf J_f = 0$ yet $\mathbf H\ne0$. Only when symmetry makes the $\mathbf M$-divergence vanish (the infinite solenoid, the toroid) does Ampère's law for $\mathbf H$ give it directly. Books disagree about which of $\mathbf B$ and $\mathbf H$ deserves "the magnetic field"; here $\mathbf B$ is the field.

**Boundary conditions.** From (3.28) and (3.41) with free surface current $\mathbf K_f$:
$$H^\perp_{\rm above} - H^\perp_{\rm below} = -(M^\perp_{\rm above} - M^\perp_{\rm below}),\qquad\mathbf H^\parallel_{\rm above} - \mathbf H^\parallel_{\rm below} = \mathbf K_f\times\hat{\mathbf n}. \tag{3.42}$$
The first is $B^\perp$ continuous rewritten; the second is the loop argument of §3.9 with only free current on the right.

**Linear media.** If $\mathbf M = \chi_m\mathbf H$ then $\mathbf B = \mu_0(1 + \chi_m)\mathbf H = \mu\mathbf H$, $\mu_r = 1 + \chi_m$ (note: $\chi_m$ multiplies $\mathbf H$, not $\mathbf B$). In a homogeneous linear medium
$$\mathbf J_b = \nabla\times(\chi_m\mathbf H) = \chi_m\mathbf J_f, \tag{3.43}$$
so bound current appears in the interior only where free current flows; the rest is on surfaces.

| Material (room temperature; Griffiths Table 6.1, approximate) | $\chi_m$ |
|---|---|
| Bismuth, gold, copper, water (diamagnetic) | $-1.6\times10^{-4}$, $-3.4\times10^{-5}$, $-9.7\times10^{-6}$, $-9.0\times10^{-6}$ |
| Oxygen (1 atm), aluminium, tungsten (paramagnetic) | $1.9\times10^{-6}$, $2.1\times10^{-5}$, $7.8\times10^{-5}$ |
| Soft iron and electrical steels (ferromagnetic, nonlinear) | $\mu_r\sim10^3$ to $10^5$, depends on $H$ and history |

Why $\chi_m$ is tiny except for ferromagnets: atomic moments are $\sim\mu_B = 9.3\times10^{-24}$ A m² and $\mu_BB/k_BT\sim2\times10^{-3}$ at 1 T (§3.12), and the diamagnetic response is $\sim10^{-5}$. Ferromagnets couple neighbouring moments with an exchange energy of order an eV $\gg k_BT$, so the alignment is cooperative.

### 3.15 The demagnetizing field and the bar magnet

For the uniformly magnetized sphere, (3.40) gives
$$\mathbf H_{\rm in} = \frac{\mathbf B_{\rm in}}{\mu_0} - \mathbf M = \tfrac23\mathbf M - \mathbf M = -\tfrac13\mathbf M. \tag{3.44}$$
Inside the magnet $\mathbf H$ points *against* $\mathbf M$ while $\mathbf B$ points with it: the **demagnetizing field**. The reason is EM-4's. Where $\mathbf J_f = 0$, $\nabla\times\mathbf H = 0$ so $\mathbf H = -\nabla\Phi_M$, and $\nabla\cdot\mathbf H = -\nabla\cdot\mathbf M$ is Poisson's equation with fictitious "magnetic charge" $\rho_M = -\nabla\cdot\mathbf M$, $\sigma_M = \mathbf M\cdot\hat{\mathbf n}$. So $\mathbf H$ is the *electric* field of EM-4 with $\mathbf P/\varepsilon_0\to\mathbf M$: (3.44) is EM-4 (3.12), and every EM-4 result for $\mathbf E$ transfers to $\mathbf H$. (This is bookkeeping, not a claim that magnetic charges exist.)

Shape matters. A thin plate magnetized *across* its faces has $\mathbf H = -\mathbf M$ inside (the slab of EM-2), so $\mathbf B = 0$ there, consistent with $B^\perp$ continuous and zero outside. A long rod magnetized *along* its axis has $\mathbf K_b = M$, a solenoid with $nI = M$, so $\mathbf B = \mu_0\mathbf M$ and $\mathbf H = 0$. The sphere is in between. A bar magnet about twice as long as it is wide sits closer to the rod: its $\mathbf B$ lines close through the magnet while its $\mathbf H$ lines leave the north face and re-enter the south face (project figure, bottom row; Example 4.7(c) for a short cylinder).

### 3.16 Ferromagnetism (qualitatively; not derivable classically)

Iron, cobalt, nickel and their alloys keep a magnetization with no applied field below the **Curie temperature** ($1043$ K for Fe, $1388$ K for Co, $627$ K for Ni). The cause is the **exchange interaction**, a quantum effect (exclusion plus Coulomb repulsion lowers the energy of parallel spins in these materials); by Bohr–van Leeuwen no classical model can supply it. SM-6's Ising model and mean-field theory reproduce the spontaneous magnetization and Curie point *given* the coupling.

A bulk sample is not uniformly magnetized: **domains** (Weiss, 1907) of aligned moments with different orientations lower the magnetostatic energy $\int B^2/2\mu_0$ by closing flux inside the material, at the cost of domain-wall energy. An applied field moves walls and rotates domains, so $\mathbf M(\mathbf H)$ is nonlinear, saturates at $M_s$ (iron: $\mu_0M_s = 2.15$ T), and is **hysteretic**: it depends on history. Cycling $H$ traces a loop with a **remanence** $B_r$ (what remains at $H = 0$) and a **coercivity** $H_c$ (the reverse field that brings $B$ to zero). The loop area is the energy dissipated per cycle per unit volume (EM-6 derives $dW = \mathbf H\cdot d\mathbf B$). **Soft** materials (silicon steel, permalloy; $H_c\lesssim100$ A/m, $\mu_r$ up to $10^5$) have thin loops and low loss: transformer and motor cores, inductors, disk-drive write heads. **Hard** materials (NdFeB, alnico; $H_c\sim10^5$ to $10^6$ A/m, $B_r\sim1$ T) keep their field: permanent magnets, recorded bits. The project simulates a loop with a Preisach model (bistable "hysterons" flipping up at $H > \alpha_i$ and down at $H < \beta_i$, averaged): remanence, coercivity and loop area appear qualitatively, and nothing explains *why* a material is ferromagnetic.

### 3.17 Magnetic circuits and the electromagnet with a gap

A core of high $\mu$ guides flux like a wire guides current. Take a core (cross-section $A_k$, length $l_k$, permeability $\mu_k$ for segment $k$) carrying $N$ turns of current $I$.

1. Ampère for $\mathbf H$ (3.41) around the core: $\sum_kH_kl_k = NI$.
2. $\nabla\cdot\mathbf B = 0$ makes the flux $\Phi = B_kA_k$ the same in every segment of a series circuit (neglecting leakage). With $H_k = B_k/\mu_k = \Phi/\mu_kA_k$:
$$NI = \Phi\sum_k\mathcal R_k,\qquad\mathcal R_k = \frac{l_k}{\mu_kA_k}. \tag{3.45}$$
This is Ohm's law with $NI$ (magnetomotive force) for the EMF, $\Phi$ for the current and the **reluctance** $\mathcal R$ for the resistance (Hopkinson's law); series reluctances add.
3. *Air gap.* A core of mean length $l$ and relative permeability $\mu_r$ with a gap $g$ in the same cross-section has $\mathcal R_{\rm core} = l/\mu_r\mu_0A$ and $\mathcal R_{\rm gap} = g/\mu_0A$, so
$$B_{\rm gap} = \frac{\mu_0NI}{g + l/\mu_r}. \tag{3.46}$$
A millimetre of air equals $\mu_r$ millimetres of iron, so a small gap dominates. That is a feature: the gap fixes $B$ by $NI$ and $g$ alone ($B\approx\mu_0NI/g$) and keeps the iron away from saturation, where $\mu$ collapses and the linear model fails. The model neglects leakage (flux that bypasses the gap, lowering $B$) and fringing (flux that bulges out of the gap, raising the effective area); Example 4.6 measures both with a 2D finite-difference solve of $\nabla\times(\nu\nabla\times\mathbf A) = \mathbf J_f$, $\nu = 1/\mu$, which for $\mathbf A = A_z(x,y)\hat{\mathbf z}$ is
$$-\nabla\cdot(\nu\nabla A_z) = J_z,\qquad\mathbf B = (\partial_yA_z,\ -\partial_xA_z). \tag{3.47}$$

## 4. Worked examples

### Example 4.1: Helmholtz coils (Part 2)
Two coils, $N = 100$ turns, $R = 10$ cm, $I = 2$ A. One coil's centre field is $\mu_0NI/2R = 1.2566$ mT (3.13). The pair at $d = R$ has $B_0 = (4/5)^{3/2}\mu_0NI/R = 1.7984$ mT, and the correction in (3.16) is $-1.8\times10^{-3}B_0$ at $z = 0.2R$. The integrator's largest deviation within $0.2R$ of the centre is $1.76\times10^{-3}$ at $d = R$ and $2.91\times10^{-2}$ at $d = 1.5R$.

### Example 4.2: A solenoid, finite and infinite, field and potential (Part 2–4)
$n = 4000$ /m, $I = 3$ A, $R = 2$ cm, $L = 0.5$ m. Ampère (3.21): $\mu_0nI = 15.080$ mT inside, $0$ outside. Biot–Savart on the finite solenoid: $15.032$ mT at the centre, matching (3.14) and $0.3\%$ below the infinite value ($L/R = 25$); $B_z = -46.1\ \mu$T at $s = 2R$ (return flux, $0.31\%$ of inside); wall jump $15.078$ vs $\mu_0K = 15.080$ mT. The flux is $\Phi_\infty = \mu_0nI\pi R^2 = 1.895\times10^{-5}$ Wb and (3.26) gives $A_\phi(1.5R) = 100.5\ \mu$T m. Numerically $\oint\mathbf A\cdot d\mathbf l$ is $0.9901\,\Phi_\infty$ at $s = 1.5R$ and $0.9697\,\Phi_\infty$ at $3R$, each equal to the numerical flux of $\mathbf B$ through the disk; the deficit is flux returning outside the finite solenoid.

### Example 4.3: The spinning shell (Part 5)
$R = 1$ cm, $\sigma = 10^{-6}$ C/m², $\omega = 100$ rad/s. (3.27): $B_{\rm in} = \tfrac23\mu_0\sigma R\omega = 0.838$ pT and $A_y(0.3R, 0, 0) = 1.257\times10^{-15}$ T m, both reproduced by 200 current rings. The surface field is $\sigma/\varepsilon_0 = 1.13\times10^5$ V/m and $\mu_0 = 1/\varepsilon_0c^2$, so $B_{\rm in} = \tfrac23(\omega R/c^2)(\sigma/\varepsilon_0)$: at a surface speed of $1$ m/s, magnetism is suppressed by $v/c^2$ relative to electricity (Purcell & Morin ch. 5–6 derive it as a relativistic effect).

### Example 4.4: Dipole field, close and far; the Earth (Part 2, 5)
A loop with $R = 10$ cm, $I = 2$ A has $m = I\pi R^2 = 0.0628$ A m². On axis, (3.33) gives $\mu_0m/2\pi z^3$ and (3.13) gives that times $(1 + R^2/z^2)^{-3/2}$: $0.9429$ at $z = 5R$ and $0.9852$ at $10R$. For the Earth take a dipole of about $8\times10^{22}$ A m² and $R_E = 6.371\times10^6$ m: $B_{\rm eq} = \mu_0m/4\pi R_E^3 = 30.9\ \mu$T and $B_{\rm pole} = 61.9\ \mu$T, inside the observed 25–65 $\mu$T.

### Example 4.5: The uniformly magnetized sphere (Part 5)
$R = 1$ cm, $M = 10^6$ A/m ($\mu_0M = 1.257$ T, NdFeB-like). (3.40): $B_{\rm in} = \tfrac23\mu_0M = 0.8378$ T; (3.44): $H_{\rm in} = -M/3 = -3.333\times10^5$ A/m; on the axis at $2R$, $B = \mu_0M/12 = 0.1047$ T. The ring model gives $0.8380$ T, $-3.331\times10^5$ A/m and $0.1046$ T. Just above the pole $B = 0.8379$ T equals the interior value ($B^\perp$ continuous), while $H^\perp$ jumps from $-3.33\times10^5$ to $+6.67\times10^5$ A/m: a jump of $M$, as (3.42) requires.

### Example 4.6: An electromagnet with an air gap (Part 7)
Square C-core, $20$ cm outside, $4$ cm legs, mean path $l = 0.636$ m (gap removed), gap $g = 4$ mm, $\mu_r = 1000$, $NI = 1000$ A (500 turns, 2 A). (3.46): $B = 1.2566\times10^{-6}\times1000/(0.004 + 0.000636) = 0.2711$ T, with $86.3\%$ of the magnetomotive force across the gap. Check with $H$: $H_{\rm core} = B/\mu = 215.7$ A/m, $H_{\rm gap} = B/\mu_0 = 2.157\times10^5$ A/m, and $137.2 + 862.8 = 1000.0$ A. The 2D solve of (3.47) on $102{,}400$ cells ($h = 1$ mm) gives $0.2658$ T at the gap centre and $0.2652$ T from gap flux$/39$ mm: $2\%$ below the model, because leakage across the window outweighs fringing for this shape. In the solved field $\oint\mathbf H\cdot d\mathbf l$ around the core centre-line is $1000.0$ A. With no gap the linear model would give $\mu_0\mu_rNI/l = 1.96$ T, at the edge of saturation; the gap caps the field at $\mu_0NI/g = 0.314$ T, and $\mu_r = 5000$ gives only $0.3045$ T.

### Example 4.7: Three problems with matter (Part 5)
(a) *Iron-filled solenoid.* $n = 1000$ /m, $I = 0.1$ A, $\mu_r = 2000$. $H = nI = 100$ A/m regardless of the filling; $B = \mu H = 0.2513$ T (empty: $0.126$ mT); $M = (\mu_r - 1)H = 1.999\times10^5$ A/m $= K_b$, which is $1999$ times the free $K$ and parallel to it. Aluminium ($\chi_m = 2.1\times10^{-5}$) gives $K_b = 2.1$ mA/m, same sign; a diamagnet's is antiparallel.
(b) *Toroid.* $N = 500$, $I = 0.5$ A, mean radius $10$ cm, $\mu_r = 2000$: $H = NI/2\pi s = 397.9$ A/m, $B = 1.000$ T, $M = 7.954\times10^5$ A/m. Across an $8$–$12$ cm section $B$ falls by $1.5$ from inner to outer edge, and the mean-radius flux estimate is $1.35\%$ below the exact $\ln1.5$ result (Your turn 2).
(c) *Magnetized cylinder.* $R = 5$ mm, $L = 20$ mm, $M = 10^6$ A/m. On the end faces $\mathbf M\parallel\hat{\mathbf n}$ so $\mathbf K_b = 0$; on the side $\mathbf K_b = M\hat{\boldsymbol\phi}$, a solenoid with $nI = M$. By (3.14), $B(\text{centre}) = \mu_0M(L/2)/\sqrt{R^2 + L^2/4} = 1.1240$ T (Biot–Savart: $1.1243$ T), below the long-rod limit $\mu_0M = 1.2566$ T; at the centre of an end face $B = \tfrac12\mu_0ML/\sqrt{R^2 + L^2} = 0.6096$ T, where $H$ jumps from $-5.149\times10^5$ (inside) to $+4.851\times10^5$ A/m (Exercise B6).

## 5. Why it matters and where it connects

- **EM-6.** Theorem 3.3 step 5 needs $\nabla\cdot\mathbf J = 0$; Maxwell's displacement current repairs it. Faraday's law couples $\mathbf E$ to $\partial_t\mathbf A$, and the energy density $\tfrac12\mathbf B\cdot\mathbf H$ with the induced EMF settles the "$U = -\mathbf m\cdot\mathbf B$" bookkeeping (§3.11) and the hysteresis energy (§3.16).
- **EM-7, EM-8.** $n = \sqrt{\varepsilon_r\mu_r}$, with $\mu_r\approx1$ at optical frequencies; retarded potentials generalize (3.24); the field tensor of LA-12 packages $\mathbf B$ with $\mathbf E$.
- **CM-3, QM-4, QC.** (3.1) follows from CM-3's Lagrangian $L = \tfrac12mv^2 - q\phi + q\mathbf v\cdot\mathbf A$ (its (4.17)); in QM-4, $\mathbf p\to\mathbf p - q\mathbf A$ and a gauge change is a local phase. Theorem 3.6 is why the flux through a superconducting loop (quantum $h/2e$) controls SQUIDs and flux-tunable qubits; a spin in $\mathbf B$ has levels $\mp m_0B$, QC-1's two-level system.
- **SM-3, SM-6.** (3.37) is SM-3's two-state paramagnet; SM-6's Ising model is the ferromagnet.
- **Mechanical engineering.** Magnetic circuits (§3.17) design motors, relays, solenoid actuators, transformers and bearings; cyclotrons and spectrometers are (3.2) at work.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5: cyclotron numbers, wires, a coil on axis, a toroid, units and the right-hand rule. B1–B8: B1 (surface independence of $I_{\rm enc}$), B3 (anti-Helmholtz), B5 (the sheet's $\mathbf A$), B6 (magnetized cylinder, $\mathbf B$ and $\mathbf H$) and B8 (`project.py`: force and torque on a small loop) are essential; B2 (non-uniform $\mathbf J$), B4 (two gauges) and B7 (a permeable wire's bound currents) are drills. C1–C3: the delta term in the dipole field, a saturating core, demagnetizing a Preisach ensemble.

## 7. Coding project

`project.py` (about ten seconds, figures in `figures/`) builds: (1) **Lorentz force**: a Boris pusher; asserts (3.2) and $|\mathbf v|$ constant to $10^{-12}$ (Theorem 3.1); spectrometer numbers. (2) **A Biot–Savart integrator** for any polygon using the exact segment formula (3.12) and its vector-potential twin; asserts (3.13), (3.14), the Helmholtz derivative and coefficient (SymPy), and the uniformity map. (3) **Ampère as a numerical loop integral** around and beside a wire and a ring; $\nabla\cdot\mathbf B$ and $\nabla\times\mathbf B$ near a ring and for Gaussian current filaments; the jump across the solenoid wall. (4) **The vector potential**: $\nabla\times\mathbf A = \mathbf B$, $\oint\mathbf A\cdot d\mathbf l$ against the flux, gauge invariance. (5) **Magnetized matter** from rings: sphere, spinning shell, cylinder; the numbers of Examples 4.3–4.7. (6) **Field-line plots** (`figures/field_lines.png`): loop, solenoid, a bar magnet's $\mathbf B$ and $\mathbf H$. (7) **Magnetic circuit** versus a 2D finite-difference solve of (3.47) (`figures/c_core.png`). (8) **Your turn**: cycloid, toroid, Preisach loop. (9) **The numbers quoted in `solutions.md`**, including the reference solution to Exercise B8.

**Build it yourself:** write `biot_savart(P, a, b, I)` for a polygon with the *midpoint rule* (one $d\mathbf l\times\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$ per segment) and watch its error against the exact formula fall as $1/N^2$ for a ring (`biot_savart_midpoint` in Part 2 is a reference). Then derive (3.12) and replace it.

## 8. Where you'll struggle

- **The right-hand rule appears in three places and you will flip a sign in each**: $d\mathbf l'\times\hat{\boldsymbol{\mathscr r}}$, the orientation of a loop against its surface (Stokes, Ampère), and $\mathbf m = I\mathbf a$. A wire with current along $+\hat{\mathbf x}$ has $\mathbf B$ along $-\hat{\mathbf y}$ above it ($\hat{\mathbf x}\times\hat{\mathbf z}$). Draw the picture and test one case before trusting any sign.
- **$\nabla\times\mathbf B = \mu_0\mathbf J$ from Biot–Savart (§3.5)** is the longest vector-calculus computation in the course. Do it once, on paper, and notice exactly where $\nabla\cdot\mathbf J = 0$ and the delta function enter.
- **Ampère's law is always true and only sometimes useful.** Students apply $B\cdot2\pi s = \mu_0I$ to a finite wire or a loop. The symmetry (§3.6) must make $\mathbf B\cdot d\mathbf l$ constant on the path; otherwise use Biot–Savart.
- **$\mathbf A$ feels unphysical and you will want to skip it. Don't.** QM-4 and QC need it. The classical content is (3.25): the loop integral is the flux, and $\mathbf A\ne0$ where $\mathbf B = 0$ outside a solenoid.
- **$\mathbf B$ versus $\mathbf H$.** Books disagree about which is "the magnetic field". Here $\mathbf B$ is the field and $\mathbf H$ is bookkeeping with $\nabla\cdot\mathbf H = -\nabla\cdot\mathbf M$. Inside a bar magnet $\mathbf H$ points against $\mathbf B$ (the $-\mathbf M/3$ of (3.44)); expect your intuition to resist, and use the project figure.
- **$\mu_0$ versus $1/\varepsilon_0$.** $\mu_0$ multiplies in (3.9) where $1/\varepsilon_0$ multiplies in Coulomb's law; $\mu_0\varepsilon_0 = 1/c^2$ is why Example 4.3 is tiny. And $I$ (A), $\mathbf K$ (A/m), $\mathbf J$ (A/m²) differ by a length each: a bound-current answer off by a length is a $\mathbf K$-versus-$\mathbf J$ slip.
- **Ferromagnetism is not derivable classically,** and this module says so (§3.12, §3.16). Learn the phenomenology and the magnetic-circuit use; the cause waits for SM-6 and quantum mechanics.

## 9. Check yourself

1. State the Lorentz force and prove in one line that a magnetic field does no work.
2. Derive $r = mv/qB$ and $\omega = qB/m$. Why is $\omega$ independent of speed?
3. What is the drift velocity in crossed $\mathbf E$ and $\mathbf B$, and why is it independent of the charge?
4. Derive $\partial_t\rho = -\nabla\cdot\mathbf J$. What does $\nabla\cdot\mathbf J = 0$ say about the path of current?
5. Write Biot–Savart for a volume current and obtain the infinite-wire field from it.
6. In deriving $\nabla\times\mathbf B = \mu_0\mathbf J$, where does the delta function appear, and where is $\nabla\cdot\mathbf J = 0$ used?
7. For which four geometries does Ampère's law determine $\mathbf B$, and what symmetry do they share?
8. Why does $\nabla\cdot\mathbf B = 0$ imply $\mathbf B = \nabla\times\mathbf A$, and how do you reach the Coulomb gauge?
9. What is $\oint\mathbf A\cdot d\mathbf l$ around a loop outside an infinite solenoid, and why is it gauge invariant?
10. State the boundary conditions for $\mathbf B$ and $\mathbf H$ at a surface current.
11. Why does the monopole term in the expansion of $\mathbf A$ vanish, and what is $\mathbf A_{\rm dip}$?
12. Give $\mathbf N$, $\mathbf F$ and $U$ for a dipole in a field. Give $\mathbf J_b$ and $\mathbf K_b$ for a magnetized object, and $\mathbf B$ and $\mathbf H$ inside a uniformly magnetized sphere.

## 10. Pair with

- **Griffiths, *Introduction to Electrodynamics*, ch. 5–6.** The same order and nearly the same derivations; ch. 5 is where §3.5's curl computation and the spinning shell come from, ch. 6 matches §3.13–3.15. Griffiths is better on worked problems and on the figures of the right-hand rule.
- **Purcell and Morin, *Electricity and Magnetism*, ch. 6 and ch. 11.** Magnetism derived as a relativistic consequence of the electric force on moving charges (read for intuition, e.g. why Example 4.3 is $c^{-2}$ small), and a long treatment of matter.
- **MIT 8.07 notes, lectures 10–14.** A more compact, more mathematical presentation of the vector potential, gauge choice and the multipole expansion; good as a second pass.
