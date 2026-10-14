# EM-2: Electrostatics: Coulomb → Gauss → Poisson; Potential, Energy, Conductors

**Time:** ~12 hours · **Prerequisites:** EM-1 · **Unlocks:** EM-3

> **Conventions in this module.** SI units. $\varepsilon_0 = 8.854\times10^{-12}$ F/m is the permittivity of free space and $k \equiv 1/4\pi\varepsilon_0 = 8.988\times10^{9}$ N m²/C² is Coulomb's constant. The **field point** $\mathbf r$ is where we want the field; the **source point** $\mathbf r'$ is where a charge sits; $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$ is the separation vector from source to field point, with length $\mathscr r$ and unit vector $\hat{\boldsymbol{\mathscr r}}$ (Griffiths' script r, as in EM-1). Charge densities: $\rho$ per volume, $\sigma$ per area, $\lambda$ per length. $d\tau'$ is a volume element at the source point; $d\mathbf a$ is an area element pointing **outward** from a closed surface. $\mathbf E$ is in V/m. The potential is $V$ (volts), with $\mathbf E = -\nabla V$ and $V(\infty) = 0$ for localized charge; for infinite distributions a finite reference point is named every time. Volumes are written $\mathcal V$ so they never clash with $V$. At a charged surface, $\hat{\mathbf n}$ points from side 1 ("below") to side 2 ("above"). Spherical coordinates $(r,\theta,\phi)$; cylindrical $(s,\phi,z)$ with $s$ the distance from the axis. "Conductor" means an ideal conductor in electrostatic equilibrium.

## 0. Where this sits

You need EM-1: the divergence and Stokes theorems, the delta-function identity $\nabla\cdot(\hat{\boldsymbol{\mathscr r}}/\mathscr r^2) = 4\pi\delta^3(\boldsymbol{\mathscr r})$ and its twin $\nabla^2(1/\mathscr r) = -4\pi\delta^3(\boldsymbol{\mathscr r})$, the theorem that a curl-free field on all of space is a gradient, and Green's first identity. CM-1 Theorem 3.3 (conservative forces) is used once.

This module turns Coulomb's law into the three tools every later E&M module uses: Gauss's law, the potential with Poisson's equation, and the field-energy formula. EM-3 solves Poisson's equation with conductors as boundaries (uniqueness, images, separation of variables, multipoles); EM-4 adds matter; EM-6 makes the field energy $\tfrac{\varepsilon_0}{2}E^2$ move.

## 1. The problem that forced this idea

By the 1760s experimenters knew that like charges repel and unlike attract, but not how the force depends on distance. Joseph Priestley, following up an observation of Benjamin Franklin's that a charged metal cup has no charge on its inside surface, argued in 1767 that this was exactly what Newton's shell theorem would predict for an inverse-square force. Henry Cavendish tested the idea directly around 1773 with two concentric metal spheres and found the exponent was 2 to within about 2 percent, but he never published; Maxwell edited his papers in 1879. The public result came from Charles-Augustin de Coulomb, who in 1785 hung a needle from a fine wire, measured the twist needed to hold two small charged balls apart at different separations, and reported (from a handful of readings) that the repulsion falls as the inverse square of the distance.

Coulomb's law describes action at a distance: a charge here pushes a charge there, instantly, through empty space. Newton had called the analogous idea for gravity "so great an absurdity" in a 1693 letter to Richard Bentley, but used it anyway because it worked. It also worked for electricity, as long as you knew where every charge was. That was the practical problem. In a conductor the charges move until they stop, and you do not know where they end up. Poisson tackled exactly this in 1812, computing how charge distributes itself over two touching conducting spheres. Summing Coulomb's law over unknown charges is not a calculation you can start.

The fix came from three directions. In 1813 Poisson showed that the gravitational potential obeys a local differential equation inside matter, extending Laplace's $\nabla^2 V = 0$ for empty space. In 1828 George Green, a miller's son in Nottingham with almost no formal schooling, self-published *An Essay on the Application of Mathematical Analysis to the Theories of Electricity and Magnetism*, which named the "potential function," proved the integral identities that bear his name, and treated conductors as boundary conditions. Almost no one read it until William Thomson (Kelvin) found a copy in 1845 and had it reprinted. Gauss stated the flux theorem for inverse-square forces in 1835 (published posthumously in 1867). Meanwhile Michael Faraday, with no mathematics to speak of, insisted from the 1830s on that the space between charges is filled with "lines of force" that are physically real, and built experiments (a foil-covered room in 1836, the ice-pail experiment in 1843) around the idea. Maxwell later turned Faraday's lines into the field equations.

So the field is not a convenience. It is what replaces "sum over every charge in the universe" with a local law, $\nabla\cdot\mathbf E = \rho/\varepsilon_0$, $\nabla\times\mathbf E = 0$, that you can solve with boundary conditions on conductors whose charge you never knew.

## 2. Intuition

**Field lines as fluid flow.** Picture each positive charge as a tap that emits an incompressible fluid at a rate proportional to its charge, and each negative charge as a drain. The flow lines are field lines; the speed of the fluid is $|\mathbf E|$. Around one tap in three dimensions, the fluid crosses spheres of area $4\pi r^2$, so its speed falls as $1/r^2$. That is Coulomb's law. The inverse square is not an arbitrary exponent: it is what "lines are neither created nor destroyed in empty space" means in three dimensions.

Now draw any closed surface and count the fluid leaving it. Fluid from taps outside the surface enters and leaves, net zero. Fluid from taps inside all leaves. So the net outflow equals the total tap strength inside, no matter what shape the surface is. You should be able to guess **Gauss's law**: flux of $\mathbf E$ out of a closed surface $= Q_{\text{enc}}/\varepsilon_0$.

**Potential as height.** Picture a landscape whose height at each point is $V$. The field points straight downhill, and its size is the steepness: $\mathbf E = -\nabla V$. A positive charge is a peak, a negative charge a pit. Walk around any closed loop on a landscape and you return to the same height, so the total "climb" $\oint\mathbf E\cdot d\mathbf l$ is zero. Guess: $\nabla\times\mathbf E = 0$, and the work to move a charge depends only on the endpoints.

**Laplace as averaging.** In a region with no charge, the landscape has no peaks or pits. The best way to say this: the height at a point equals the average height on a small circle around it. That is $\nabla^2 V = 0$, and it is literally the update rule the project uses to solve it on a grid. With charge present, the landscape is pushed up locally, like a membrane under pressure: $\nabla^2 V = -\rho/\varepsilon_0$.

**Conductors as flat ground.** In a metal the charges are free to roll. They roll until nothing pushes them, so inside a conductor the ground is flat: $\mathbf E = 0$ and $V$ is constant. The surface is a contour line of the landscape, and water runs downhill perpendicular to contour lines, so field lines leave a conductor at right angles. Any net charge must sit on the surface, because a bump of charge inside would make a peak.

**A jump in the slope.** A sheet of charge emits lines on both sides, so the field just above minus the field just below is the line density of the sheet: guess $\Delta E_\perp = \sigma/\varepsilon_0$. The height itself does not jump, only its slope.

**Energy.** Pushing two like charges together crowds field lines between them and costs work. Guess where that work goes: into the field, at a density that grows with $E^2$ (since $W \propto q^2$ and $E\propto q$).

## 3. Formal development

### 3.1 Coulomb's law and superposition

A point charge $q$ at $\mathbf r'$ exerts on a charge $Q$ at $\mathbf r$ the force
$$\mathbf F = \frac{1}{4\pi\varepsilon_0}\frac{qQ}{\mathscr r^2}\,\hat{\boldsymbol{\mathscr r}}. \tag{3.1}$$
Forces from several charges add as vectors (**superposition**, an experimental fact: the force between two charges is unchanged by the presence of a third). Dividing by the test charge $Q$ defines the **electric field** of charges $q_1,\dots,q_N$:
$$\mathbf E(\mathbf r) = \frac{1}{4\pi\varepsilon_0}\sum_{i=1}^N \frac{q_i}{\mathscr r_i^2}\hat{\boldsymbol{\mathscr r}}_i, \qquad \boldsymbol{\mathscr r}_i = \mathbf r - \mathbf r_i. \tag{3.2}$$
For a continuous distribution replace $q_i$ by $\rho(\mathbf r')\,d\tau'$ (or $\sigma\,da'$, $\lambda\,dl'$):
$$\mathbf E(\mathbf r) = \frac{1}{4\pi\varepsilon_0}\int \rho(\mathbf r')\frac{\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\,d\tau'. \tag{3.3}$$
The test charge is assumed small enough not to move the sources. In statics, $\mathbf E$ is a summary of what force a charge would feel; from EM-6 on it carries energy and momentum and becomes the primary object.

### 3.2 The field-line picture

A **field line** is a curve tangent to $\mathbf E$ at every point: a solution of $d\mathbf r/dt = \mathbf E(\mathbf r)$. Rules:

1. Lines start on positive charges and end on negative charges or at infinity. (In a charge-free region they cannot begin or end; this is $\nabla\cdot\mathbf E = 0$ there, derived in §3.3.)
2. Lines never cross, because $\mathbf E$ has one direction at each point (except where $\mathbf E = 0$).
3. If you draw a number of lines proportional to each charge, the number crossing a small area $da$ perpendicular to $\mathbf E$ is proportional to $|\mathbf E|\,da$.

Rule 3 is consistent only because of the inverse square: $N$ lines from a point charge cross a sphere of area $4\pi r^2$, giving density $\propto 1/r^2 \propto |\mathbf E|$. The project draws the lines of a dipole and a quadrupole; Exercise B4 shows dipole lines are the curves $r = r_0\sin^2\theta$.

### 3.3 The divergence of $\mathbf E$ and Gauss's law

▶ **Try it first.** Take the divergence of (3.3) with respect to $\mathbf r$. Only $\boldsymbol{\mathscr r}$ depends on $\mathbf r$. Then use EM-1's delta-function identity.

**Derivation.**

1. The derivatives act on the field point $\mathbf r$; $\rho(\mathbf r')$ is a constant for them. So
   $$\nabla\cdot\mathbf E = \frac{1}{4\pi\varepsilon_0}\int \rho(\mathbf r')\,\nabla\cdot\!\left(\frac{\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\right)d\tau'.$$
2. EM-1: $\nabla\cdot(\hat{\boldsymbol{\mathscr r}}/\mathscr r^2) = 4\pi\delta^3(\mathbf r - \mathbf r')$. Zero away from the source, with total flux $4\pi$ through any surface around it.
3. The delta function sifts out $\rho(\mathbf r)$:
   $$\nabla\cdot\mathbf E = \frac{\rho}{\varepsilon_0}. \tag{3.4}$$
4. Integrate (3.4) over a volume $\mathcal V$ with closed boundary $S$ and apply the divergence theorem:
   $$\oint_S \mathbf E\cdot d\mathbf a = \frac{1}{\varepsilon_0}\int_{\mathcal V}\rho\,d\tau = \frac{Q_{\text{enc}}}{\varepsilon_0}. \tag{3.5}$$
∎

**Theorem 3.1 (Gauss's law).** For every electrostatic field produced by charges obeying (3.3), and every closed surface $S$, (3.5) holds. Equivalently, (3.4) holds at every point.

■ **What just happened.** The only place the inverse-square law entered was step 2. $\nabla\cdot(\hat{\boldsymbol{\mathscr r}}/\mathscr r^n)$ vanishes away from the source only for $n = 2$; for a force $\propto 1/r^{2+\epsilon}$ the flux through a sphere would depend on its radius and (3.5) would be false. Gauss's law *is* the statement that the exponent is 2. That is why Cavendish's sphere-inside-a-sphere experiment measured the exponent: a nonzero $\epsilon$ would put charge on the inner sphere. (A 1971 version by Williams, Faller and Hill put $|\epsilon|$ below a few parts in $10^{16}$.)

**The converse: Gauss plus symmetry gives Coulomb.**

▶ **Try it first.** For a point charge $q$ at the origin, use rotational symmetry to decide which direction $\mathbf E$ can point and what it can depend on. Then apply (3.5) to a sphere.

1. Any rotation about the origin carries the charge to itself, so it must carry the field to itself. Take a point $\mathbf r$. Rotating by $180°$ about the axis through $\mathbf r$ fixes the point and the charge but reverses any component of $\mathbf E(\mathbf r)$ perpendicular to $\hat{\mathbf r}$. That component must therefore be zero: $\mathbf E = E\,\hat{\mathbf r}$.
2. Rotations carry any point at radius $r$ to any other, so $E = E(r)$.
3. On the sphere of radius $r$, $\mathbf E\cdot d\mathbf a = E(r)\,da$ with $E(r)$ constant, so the flux is $4\pi r^2E(r)$.
4. Gauss: $4\pi r^2E(r) = q/\varepsilon_0$, so $\mathbf E = q\hat{\mathbf r}/(4\pi\varepsilon_0r^2)$, and the force on $Q$ is (3.1). ∎

■ **What just happened.** Symmetry fixed the direction and the variable; Gauss supplied one number. That is the template for every Gauss's-law calculation. Gauss's law alone does not determine a field (you could add any divergence-free field); by the Helmholtz theorem (EM-1) you also need the curl, which §3.5 supplies, and a condition at infinity.

### 3.4 When Gauss's law is useful: three symmetries

Gauss's law is true for every closed surface and every charge distribution. It computes $\mathbf E$ only when you can find a surface on which $\mathbf E\cdot\hat{\mathbf n}$ is constant (on some pieces) and zero (on the rest), so that $|\mathbf E|$ comes out of the integral. That needs a symmetry that carries any point of the surface to any other point and leaves the charge unchanged. The surfaces that a continuous family of rigid motions can sweep out point-to-point in this way are **spheres** (rotations about a point), **infinite cylinders** (rotations about an axis plus translations along it), and **planes** (translations within the plane). Those are the three useful cases: spherical, cylindrical, planar. A finite rod, a cube, or a disk has no such surface; Gauss's law is still true for them but gives one equation for an unknown function (§8).

**Spherical shell** (radius $R$, charge $Q$ spread uniformly). By the argument of §3.3, $\mathbf E = E(r)\hat{\mathbf r}$. A Gaussian sphere of radius $r<R$ encloses nothing, so $E = 0$ inside. For $r>R$ it encloses $Q$, so $E = Q/(4\pi\varepsilon_0r^2)$: outside, the shell is indistinguishable from a point charge at its centre (Newton's shell theorem, CM-5). At $r = R$ the field jumps from $0$ to $Q/(4\pi\varepsilon_0R^2) = \sigma/\varepsilon_0$, previewing §3.8.

**Infinite plane** (charge $\sigma$ per area in $z = 0$). Translations in $x$ and $y$ leave everything unchanged, so $\mathbf E$ depends only on $z$. A $180°$ rotation about the vertical line through the field point reverses $E_x, E_y$, so $\mathbf E = E_z(z)\hat{\mathbf z}$. Reflection $z\to-z$ gives $E_z(-z) = -E_z(z)$. A pillbox with faces of area $A$ at heights $\pm z$ has flux $2E_z(z)A$ (none through the sides) and encloses $\sigma A$:
$$\mathbf E = \frac{\sigma}{2\varepsilon_0}\operatorname{sign}(z)\,\hat{\mathbf z}.$$
The field does not fall off with distance, because lines from an infinite sheet have nowhere to spread.

**Coaxial cylinders.** A thin cylindrical shell of radius $a$ carries $+\lambda$ per length; a coaxial shell of radius $b > a$ carries $-\lambda$. Translation along and rotation about the axis make $\mathbf E$ depend only on $s$; reflection through the plane $z = $ const kills $E_z$; reflection through the plane containing the axis and the field point kills $E_\phi$. So $\mathbf E = E(s)\hat{\mathbf s}$. A Gaussian cylinder of radius $s$ and length $\ell$ has flux $2\pi s\ell E(s)$ (the end caps contribute zero since $\mathbf E\perp\hat{\mathbf n}$ there) and encloses $0$ for $s<a$, $\lambda\ell$ for $a<s<b$, and $0$ for $s>b$:
$$\mathbf E = \frac{\lambda}{2\pi\varepsilon_0 s}\hat{\mathbf s}\quad(a<s<b),\qquad \mathbf E = 0\ \text{otherwise}.$$
The field is confined between the conductors. That is what a coaxial cable is for.

The infinite line and the solid sphere are Examples 4.2 and 4.1.

### 3.5 The curl of $\mathbf E$ and the potential

**Theorem 3.2.** Every electrostatic field satisfies
$$\nabla\times\mathbf E = 0. \tag{3.6}$$

▶ **Try it first.** Write $\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$ as the gradient (with respect to $\mathbf r$) of a function of $\mathscr r$.

**Derivation.**

1. $1/\mathscr r$ depends on $\mathbf r$ only through the distance $\mathscr r$, with $d(1/\mathscr r)/d\mathscr r = -1/\mathscr r^2$ in the direction $\hat{\boldsymbol{\mathscr r}}$. So $\nabla(1/\mathscr r) = -\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$.
2. Substitute into (3.3) and pull $\nabla$ (which acts on $\mathbf r$) outside the $\mathbf r'$ integral:
   $$\mathbf E(\mathbf r) = -\nabla\left[\frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\mathbf r')}{\mathscr r}\,d\tau'\right].$$
3. The curl of a gradient is zero (EM-1). ∎

Step 2 hands us the potential. Define
$$V(\mathbf r) = \frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\mathbf r')}{\mathscr r}\,d\tau', \qquad \mathbf E = -\nabla V. \tag{3.8}$$
For localized $\rho$ the integral tends to $0$ as $r\to\infty$, which is the convention $V(\infty) = 0$. By the fundamental theorem for gradients (EM-1),
$$V(\mathbf b) - V(\mathbf a) = -\int_{\mathbf a}^{\mathbf b}\mathbf E\cdot d\mathbf l, \tag{3.7}$$
independent of path, and $\oint\mathbf E\cdot d\mathbf l = 0$ for every loop. The electrostatic force $Q\mathbf E$ is conservative in the sense of CM-1 Theorem 3.3 with potential energy $QV$.

■ **What just happened.** One scalar function carries all the information in three field components, and it superposes (it is linear in $\rho$). Two routes to it: (3.7) needs only $\nabla\times\mathbf E = 0$ on all of space (EM-1's curl-free-implies-gradient theorem) and works for any distribution once you pick a reference point; (3.8) needs the charge to be localized, and fails for infinite lines and planes (Example 4.2). The project checks $\oint\mathbf E\cdot d\mathbf l = 0$ around a loop threading a quadrupole, and that central differences of $V$ reproduce $\mathbf E$.

### 3.6 Poisson's and Laplace's equations

Substitute $\mathbf E = -\nabla V$ into (3.4):
$$\nabla^2 V = -\frac{\rho}{\varepsilon_0}\quad\text{(Poisson)},\qquad \nabla^2V = 0 \ \text{ where } \rho = 0\quad\text{(Laplace)}. \tag{3.9}$$
Consistency check: apply $\nabla^2$ to (3.8) and use $\nabla^2(1/\mathscr r) = -4\pi\delta^3(\boldsymbol{\mathscr r})$ (EM-1): $\nabla^2V = \frac{1}{4\pi\varepsilon_0}\int\rho(\mathbf r')(-4\pi)\delta^3(\mathbf r-\mathbf r')\,d\tau' = -\rho/\varepsilon_0$. So (3.8) *is* the solution of (3.9) that vanishes at infinity, and $1/(4\pi\mathscr r)$ is the Green's function of $-\nabla^2$ (F3's impulse response, in space).

You now have six arrows between $\rho$, $\mathbf E$, $V$: $\mathbf E$ from $\rho$ by (3.3), $\rho$ from $\mathbf E$ by (3.4), $V$ from $\rho$ by (3.8), $\rho$ from $V$ by (3.9), $V$ from $\mathbf E$ by (3.7), $\mathbf E$ from $V$ by $-\nabla$. The integral formulas need every charge known. Poisson's equation needs only the charge you know plus **boundary conditions**, which is how conductors with unknown charge get handled (EM-3; project Part 5).

**Standard potentials.** Point charge: $V = q/(4\pi\varepsilon_0 r)$. Spherical shell, integrating $E$ in from infinity with (3.7): $V = Q/(4\pi\varepsilon_0r)$ outside, and inside, where $E = 0$, $V$ stays at its surface value $Q/(4\pi\varepsilon_0R)$: zero field does not mean zero potential. Solid sphere: Example 4.1. Line: Example 4.2. Dipole: Example 4.3. Infinite plane: $V = -\sigma|z|/(2\varepsilon_0)$ plus a constant, unbounded in both directions, so as with the line, $V(\infty) = 0$ is unavailable.

### 3.7 Boundary conditions at a surface charge

▶ **Try it first.** Apply (3.5) to a flat pillbox straddling the surface, and $\oint\mathbf E\cdot d\mathbf l = 0$ to a thin rectangle straddling it. Let the thickness go to zero.

**Derivation.** Let a surface carry density $\sigma$, with $\hat{\mathbf n}$ from side 1 (below) to side 2 (above), and let $\mathbf E_1,\mathbf E_2$ be the fields just below and above.

1. *Pillbox.* Faces of area $A$ on each side, thickness $\epsilon$. The flux is $(\mathbf E_2 - \mathbf E_1)\cdot\hat{\mathbf n}\,A$ plus a side contribution proportional to $\epsilon$, which vanishes as $\epsilon\to0$ because $\mathbf E$ is finite. The enclosed charge is $\sigma A$. So
   $$E_2^\perp - E_1^\perp = \frac{\sigma}{\varepsilon_0}. \tag{3.10}$$
2. *Loop.* A rectangle with sides of length $l$ parallel to the surface and height $\epsilon\to0$: $\oint\mathbf E\cdot d\mathbf l = (\mathbf E_2-\mathbf E_1)\cdot\hat{\mathbf t}\,l = 0$ by (3.6), for any tangent direction $\hat{\mathbf t}$. So
   $$\mathbf E_2^\parallel = \mathbf E_1^\parallel. \tag{3.11}$$
   Together: $\mathbf E_2 - \mathbf E_1 = (\sigma/\varepsilon_0)\hat{\mathbf n}$.
3. *Potential.* $V_2 - V_1 = -\int_1^2\mathbf E\cdot d\mathbf l$ over a path of length $\epsilon\to0$ with $\mathbf E$ finite, so $V$ is continuous.
4. *Normal derivative.* With $E^\perp = -\partial V/\partial n$, (3.10) becomes
   $$\frac{\partial V_2}{\partial n} - \frac{\partial V_1}{\partial n} = -\frac{\sigma}{\varepsilon_0}. \tag{3.12}$$
∎

■ **What just happened.** A surface charge puts a kink in $V$ (the slope jumps) but not a step. Check against §3.4: for the plane, $E_2^\perp - E_1^\perp = \sigma/2\varepsilon_0 - (-\sigma/2\varepsilon_0) = \sigma/\varepsilon_0$; for the shell, $\sigma/\varepsilon_0 - 0$. The sign convention is "above minus below, with $\hat{\mathbf n}$ pointing up." Get it backwards and every induced charge in EM-3 flips sign.

### 3.8 Work and energy

Moving a test charge $Q$ slowly from $\mathbf a$ to $\mathbf b$ takes force $-Q\mathbf E$ from you, so by (3.7) your work is
$$W = -Q\int_{\mathbf a}^{\mathbf b}\mathbf E\cdot d\mathbf l = Q\,[V(\mathbf b) - V(\mathbf a)]. \tag{3.13}$$
$V$ is potential energy per unit charge.

**Energy of point charges.**

▶ **Try it first.** Bring three charges in from infinity one at a time and add up the work. Then rewrite the result as a sum over all $i$ and all $j\ne i$.

1. The first charge costs nothing: there is no field yet.
2. The second costs $q_2V_1(\mathbf r_2) = kq_1q_2/\mathscr r_{12}$ by (3.13).
3. The third costs $kq_3(q_1/\mathscr r_{13} + q_2/\mathscr r_{23})$.
4. In general each pair is paid for once, when its second member arrives: $W = k\sum_{i<j}q_iq_j/\mathscr r_{ij}$. The order of assembly does not matter.
5. A sum over all ordered pairs $i\neq j$ counts each pair twice, so
   $$W = \frac12\sum_{i=1}^N q_i\,V(\mathbf r_i),\qquad V(\mathbf r_i) = k\sum_{j\ne i}\frac{q_j}{\mathscr r_{ij}}, \tag{3.14}$$
   where $V(\mathbf r_i)$ is the potential at $q_i$ due to all the *other* charges.

■ **What just happened.** The $\tfrac12$ is a double-counting correction: in $\sum_iq_iV(\mathbf r_i)$, the pair $(1,2)$ appears as $q_1\times$(potential of $q_2$) and again as $q_2\times$(potential of $q_1$). The restriction $j\ne i$ means (3.14) leaves out the energy needed to assemble each point charge itself. Hold that thought.

**Continuous charge.** Replace $q_i$ by $\rho\,d\tau$:
$$W = \frac12\int\rho\,V\,d\tau. \tag{3.15}$$
Now $V$ is the full potential. For a smooth $\rho$ the potential of an element $d\tau$ at its own location is of order $\rho\,(d\tau)^{2/3}$, which vanishes as $d\tau\to0$, so "other charges" and "all charges" agree.

**Theorem 3.3 (Energy in the field).** For a localized charge distribution,
$$W = \frac{\varepsilon_0}{2}\int_{\text{all space}}E^2\,d\tau. \tag{3.16}$$

▶ **Try it first.** Use (3.9) to replace $\rho$ by $-\varepsilon_0\nabla^2V$ in (3.15), then integrate by parts once (Green's first identity, EM-1).

**Derivation.**

1. Over any volume $\mathcal V$ containing all the charge: $W = -\dfrac{\varepsilon_0}{2}\displaystyle\int_{\mathcal V}V\,\nabla^2V\,d\tau$.
2. Product rule: $\nabla\cdot(V\nabla V) = |\nabla V|^2 + V\nabla^2V$.
3. Substitute and use the divergence theorem on the first term:
   $$W = \frac{\varepsilon_0}{2}\left[\int_{\mathcal V}|\nabla V|^2\,d\tau - \oint_{S}V\,\nabla V\cdot d\mathbf a\right].$$
4. Let $\mathcal V$ be a ball of radius $r\to\infty$. Far away $V\sim 1/r$, $|\nabla V|\sim1/r^2$, area $\sim r^2$, so the surface term is $\sim 1/r\to0$. The volume integral grows to the all-space integral.
5. $|\nabla V|^2 = E^2$. ∎

■ **What just happened.** The same number is an integral over where the charge is, (3.15), and an integral over where the field is, (3.16). The integrands differ point by point: in Example 4.4, five sixths of the field energy of a charged ball lies outside the ball, where $\rho = 0$ and the integrand of (3.15) is zero. Both bookkeepings are correct for the total (§8). Also note $W$ is quadratic, so energies do not superpose: $W[\mathbf E_1+\mathbf E_2] = W_1 + W_2 + \varepsilon_0\int\mathbf E_1\cdot\mathbf E_2\,d\tau$ (Exercise B5).

**The self-energy problem.** Apply (3.16) to a single point charge, $E = kq/r^2$:
$$W = \frac{\varepsilon_0}{2}\int_0^\infty\left(\frac{kq}{r^2}\right)^2 4\pi r^2\,dr = \frac{kq^2}{2}\int_0^\infty\frac{dr}{r^2} = \infty.$$
With the integral started at radius $a$ instead, $W = kq^2/(2a)$. This is the energy (3.14) left out. It is not an artifact of the formula. For an electron, the field energy outside $a = 1$ fm is already $1.41\,m_ec^2$, and outside $10^{-3}$ fm it is $1.41\times10^3\,m_ec^2$ (project Part 3), while scattering experiments show no electron structure at all down to far smaller distances. Classical electrodynamics has no fix: a point charge has infinite field energy. Quantum electrodynamics handles the divergence by absorbing it into the measured mass (renormalization), which is a working procedure, not a classical explanation. What saves you in practice: for a fixed set of point charges the self-energies are constants that cancel from every energy *difference*, so (3.14) is the right formula for forces and work.

### 3.9 Conductors

**Theorem 3.4 (Conductor in equilibrium).** For a conductor in electrostatic equilibrium: (i) $\mathbf E = 0$ inside the material; (ii) $\rho = 0$ inside; (iii) any net charge lies on the surface; (iv) the conductor is an equipotential; (v) just outside, $\mathbf E$ is perpendicular to the surface, $\mathbf E = (\sigma/\varepsilon_0)\hat{\mathbf n}$, so $\sigma = -\varepsilon_0\,\partial V/\partial n$.

▶ **Try it first.** Derive each from the one before, using (3.4), (3.7), (3.10), (3.11).

*Proof.* (i) Free charges in a field feel force $q\mathbf E$ and move; equilibrium means no motion, so $\mathbf E = 0$ wherever there are free charges, i.e., throughout the material. (That such an arrangement exists for any external field is part of EM-3's uniqueness theory; physically, the induced charges move until their field cancels the applied one, in a time of order $\varepsilon_0/\sigma_{\text{cond}} \approx 1.5\times10^{-19}$ s for copper, a naive estimate that the true electron collision time of about $10^{-14}$ s replaces; either way, effectively instantly.) (ii) $\rho = \varepsilon_0\nabla\cdot\mathbf E = 0$ by (3.4). (iii) By (ii) there is no charge in the interior, so net charge can only sit on the surface. (iv) For any two points in a connected conductor, (3.7) along a path inside gives $V(\mathbf b)-V(\mathbf a) = 0$. (v) Tangential $\mathbf E$ is continuous (3.11) and zero inside, so zero just outside. (3.10) with $E_1^\perp = 0$ gives $E_2^\perp = \sigma/\varepsilon_0$. ∎

**Induced charge.** Bring a positive charge near a neutral conductor. Negative charge gathers on the near side and positive on the far side, in exactly the pattern that cancels the external field inside. The near (negative) charge is closer, so the force is always attractive. The pattern itself is a boundary-value problem: EM-3 solves it by images, and project Part 6 solves it numerically for a conductor in a uniform field.

**Theorem 3.5 (Empty cavity).** If a cavity inside a conductor contains no charge, then $\mathbf E = 0$ throughout the cavity and the cavity wall carries no charge, whatever the fields outside the conductor.

*Proof.* Suppose $\mathbf E \ne 0$ at some point in the cavity. Follow the field line through it in both directions. With no charge in the cavity it cannot begin or end there, and it cannot close on itself (that would make $\oint\mathbf E\cdot d\mathbf l > 0$), so both ends are on the wall. Along the line $\mathbf E\cdot d\mathbf l > 0$, so by (3.7) the potential at its start exceeds that at its end. Both ends are on one conductor, an equipotential. Contradiction. So $\mathbf E = 0$ in the cavity, and by (3.10) $\sigma = 0$ on the wall. ∎ (EM-3's uniqueness theorem gives a second proof: $V = $ const solves Laplace's equation in the cavity with the right boundary value, so it is the solution.)

This is the **Faraday cage**: in 1836 Faraday built a foil-covered cube about 12 feet on a side, sat inside it with electroscopes while it was charged to the point of sparking outside, and detected nothing. A wire mesh works too, as long as the holes are small compared with the distances you care about.

**The force on a surface charge.** At a charged surface the field jumps (3.10). Which field pushes the charge?

▶ **Try it first.** Split the field near a small patch into the patch's own field and everything else.

1. Write $\mathbf E = \mathbf E_{\text{patch}} + \mathbf E_{\text{other}}$ near a small patch of the surface.
2. Seen from very close, the patch is an infinite plane: $\mathbf E_{\text{patch}} = \pm(\sigma/2\varepsilon_0)\hat{\mathbf n}$ just above and below (§3.4).
3. $\mathbf E_{\text{other}}$ has its sources elsewhere, so it is continuous across the patch.
4. Hence $\mathbf E_2 = \mathbf E_{\text{other}} + \frac{\sigma}{2\varepsilon_0}\hat{\mathbf n}$ and $\mathbf E_1 = \mathbf E_{\text{other}} - \frac{\sigma}{2\varepsilon_0}\hat{\mathbf n}$, so $\mathbf E_{\text{other}} = \tfrac12(\mathbf E_1+\mathbf E_2)$.
5. The patch cannot push itself, so the force per unit area is
   $$\mathbf f = \sigma\,\mathbf E_{\text{avg}} = \frac{\sigma}{2}(\mathbf E_1 + \mathbf E_2). \tag{3.17}$$
6. For a conductor $\mathbf E_1 = 0$, $\mathbf E_2 = (\sigma/\varepsilon_0)\hat{\mathbf n}$: $\ \mathbf f = \dfrac{\sigma^2}{2\varepsilon_0}\hat{\mathbf n} = \dfrac{\varepsilon_0}{2}E^2\,\hat{\mathbf n}$. ∎

■ **What just happened.** The surface charge feels the *average* of the fields on its two sides. For a conductor that is an outward pressure (for either sign of $\sigma$) equal to the energy density $\tfrac{\varepsilon_0}{2}E^2$ just outside. That equality is virtual work: push a patch of area $A$ outward by $dx$ at fixed charge; the field outside is unchanged, a slab of volume $A\,dx$ that held energy $\tfrac{\varepsilon_0}{2}E^2A\,dx$ now has none, and that energy drop is the work done by the force. Exercise B6 checks it on a sphere. Faraday's picture of field lines as stretched strings that also push sideways becomes the Maxwell stress tensor in EM-6.

### 3.10 Capacitors

Two conductors carry $+Q$ and $-Q$; let $V$ be the potential of the positive one minus the negative one. Doubling $Q$ doubles every field and every potential difference (everything is linear), so $Q/V$ depends only on geometry. The **capacitance** is $C = Q/V$, in farads.

▶ **Try it first.** For each geometry: find $\mathbf E$ between the conductors by Gauss, integrate it from one to the other with (3.7), divide.

**Parallel plates** (area $A$, gap $d$, edges ignored). Each plate is a sheet; between them the fields $\sigma/2\varepsilon_0$ of the two sheets add to $E = \sigma/\varepsilon_0 = Q/(\varepsilon_0A)$, and outside they cancel. $V = Ed = Qd/(\varepsilon_0A)$, so
$$C = \frac{\varepsilon_0A}{d}. \tag{3.18a}$$
**Spherical** (inner radius $a$ with $+Q$, outer radius $b$ with $-Q$). Between them $E = Q/(4\pi\varepsilon_0r^2)$; $V = \frac{Q}{4\pi\varepsilon_0}\left(\frac1a-\frac1b\right)$, so
$$C = 4\pi\varepsilon_0\frac{ab}{b-a}. \tag{3.18b}$$
As $b\to\infty$, $C\to4\pi\varepsilon_0a$: the capacitance of an isolated sphere.

**Cylindrical** (length $L$, radii $a<b$, $\lambda = Q/L$). Between them $E = \lambda/(2\pi\varepsilon_0s)$ (§3.4); $V = \frac{\lambda}{2\pi\varepsilon_0}\ln\frac ba$, so
$$C = \frac{2\pi\varepsilon_0L}{\ln(b/a)}. \tag{3.18c}$$

**Energy.** Charge the capacitor by moving charge $dq$ from the negative plate to the positive one when the charge is already $q$. By (3.13) this costs $dq\cdot(q/C)$. Integrating from $0$ to $Q$:
$$W = \int_0^Q\frac{q}{C}\,dq = \frac{Q^2}{2C} = \frac12CV^2 = \frac12QV. \tag{3.19}$$
Check with (3.16) for parallel plates: $\tfrac{\varepsilon_0}{2}E^2\cdot Ad = \tfrac{\varepsilon_0}{2}(V/d)^2Ad = \tfrac12(\varepsilon_0A/d)V^2$. ✓

**Combinations.** *Parallel:* all capacitors share $V$, charges add, $Q = \sum C_iV$, so $C = \sum C_i$. *Series:* the conductors between neighbouring capacitors are isolated and neutral, so each carries $\pm Q$ induced charge and every capacitor holds the same $Q$; voltages add, $V = \sum Q/C_i$, so $1/C = \sum 1/C_i$.

### 3.11 Charge in a cavity, and a preview of uniqueness

Put a charge $q$ inside a cavity in a neutral conductor. A Gaussian surface lying inside the metal and surrounding the cavity has $\mathbf E = 0$ on it, so it encloses zero charge: the cavity wall carries exactly $-q$. The conductor is neutral, so $+q$ sits on the outer surface. How the $-q$ is spread over the wall depends on where $q$ is and on the cavity's shape. How the $+q$ is spread over the outer surface does not depend on either. The reason is EM-3's **second uniqueness theorem**: in the region outside a set of conductors, the field is completely determined by the conductors' shapes and their total charges. The outside world sees "a conductor of this shape carrying total charge $q$," and that fixes the field outside. Faraday's ice-pail experiment (1843) is this statement: a charged ball lowered into a metal pail gives the same reading on an electroscope wired to the pail's outside wherever it hangs, and the same reading again after it touches the pail. Example 4.6 works it out for a sphere.

### 3.12 One equation, several subjects

Poisson's equation is the steady state of every "flux = $-$gradient, divergence of flux = source" law.

- **Gravity** (CM-5). $\mathbf g = -\nabla\Phi$, $\nabla\cdot\mathbf g = -4\pi G\rho_m$, so $\nabla^2\Phi = 4\pi G\rho_m$. Same equation with $1/(4\pi\varepsilon_0)\to -G$: like masses attract. The shell theorem is Gauss's law for gravity; Exercise B1 uses it.
- **Steady heat conduction.** Fourier's law $\mathbf q = -\kappa\nabla T$ with heat source density $\dot s$ in steady state gives $\nabla^2T = -\dot s/\kappa$. Temperature plays $V$, heat flux plays $\varepsilon_0\mathbf E$, and an isothermal boundary plays a conductor.
- **Diffusion as a solver.** The heat equation $\partial_tu = \alpha\nabla^2u$ relaxes to the solution of Laplace's equation with the same boundary values. Discretize it with the five-point Laplacian and explicit Euler at the largest stable step $\alpha\,\Delta t = h^2/4$ and the update becomes $u_{ij}\leftarrow$ average of its four neighbours: exactly the Jacobi iteration of LA-13. Relaxation methods are heat flow. They are slow for the reason heat flow is slow: the smoothest mode on a grid of $N$ points per side decays by a factor $\cos(\pi/(N+1))$ per sweep (LA-13 (3.24)), so the number of sweeps grows like $N^2$ (Exercise C2).
- **Variational form.** Laplace's equation is the Euler–Lagrange equation of $\int|\nabla V|^2\,d\tau$ (F6, Exercise C3, the Dirichlet principle). With conductor potentials fixed, the true potential minimizes the field energy (3.16). Exercise C3 here proves the companion statement for fixed charges (Thomson's theorem).

## 4. Worked examples

### Example 4.1: The uniformly charged solid sphere

Charge $Q$ is spread uniformly through a ball of radius $R$, so $\rho = 3Q/(4\pi R^3)$. Find $\mathbf E$ and $V$ everywhere.

*Field.* By §3.3's symmetry argument, $\mathbf E = E(r)\hat{\mathbf r}$. Gaussian sphere of radius $r$:

- $r\ge R$: $Q_{\text{enc}} = Q$, so $E = \dfrac{Q}{4\pi\varepsilon_0r^2}$.
- $r<R$: $Q_{\text{enc}} = \rho\cdot\tfrac43\pi r^3 = Q r^3/R^3$, so $4\pi r^2E = Qr^3/(\varepsilon_0R^3)$ and $E = \dfrac{Qr}{4\pi\varepsilon_0R^3}$.

*Potential*, by integrating (3.7) inward from $V(\infty) = 0$:

- $r\ge R$: $V(r) = -\displaystyle\int_\infty^r\frac{Q}{4\pi\varepsilon_0r'^2}dr' = \frac{Q}{4\pi\varepsilon_0r}$.
- $r<R$: continue from $V(R) = Q/(4\pi\varepsilon_0R)$:
  $$V(r) = \frac{Q}{4\pi\varepsilon_0R} - \int_R^r\frac{Qr'}{4\pi\varepsilon_0R^3}dr' = \frac{Q}{4\pi\varepsilon_0R} + \frac{Q(R^2-r^2)}{8\pi\varepsilon_0R^3} = \frac{Q}{8\pi\varepsilon_0R}\left(3 - \frac{r^2}{R^2}\right).$$

So $V(0) = \tfrac32V(R)$. Check against (3.9) in spherical coordinates: $\frac1{r^2}\frac{d}{dr}\left(r^2\frac{dV}{dr}\right) = \frac{1}{r^2}\frac{d}{dr}\left(-\frac{Qr^3}{4\pi\varepsilon_0R^3}\right) = -\frac{3Q}{4\pi\varepsilon_0R^3} = -\frac{\rho}{\varepsilon_0}$ ✓, and $V$ and $dV/dr$ are continuous at $R$ (no surface charge, (3.12)).

*The plots* (project, `charged_sphere.png`). $E$ rises linearly from zero at the centre, peaks at $r = R$, and falls as $1/r^2$, continuous but with a corner at $R$. $V$ is an inverted parabola inside joined smoothly to $1/r$ outside: the "hill" of §2 has a rounded top.

*Numbers.* $Q = 1$ nC, $R = 0.1$ m: $E(R) = 898.8$ V/m, $V(R) = 89.88$ V, $V(0) = 134.81$ V. The project integrates $E$ inward numerically and matches the closed form to $10^{-9}$.

### Example 4.2: The infinite line and the log potential

A straight line carries charge $\lambda$ per length. Find $\mathbf E$, then try to find $V$.

*Field.* Cylindrical symmetry (§3.4) gives $\mathbf E = E(s)\hat{\mathbf s}$. A Gaussian cylinder of radius $s$, length $\ell$: $2\pi s\ell E = \lambda\ell/\varepsilon_0$, so
$$\mathbf E = \frac{\lambda}{2\pi\varepsilon_0s}\hat{\mathbf s}.$$
*Potential with $V(\infty) = 0$.* $V(s) = -\int_\infty^s\frac{\lambda}{2\pi\varepsilon_0s'}ds' = \frac{\lambda}{2\pi\varepsilon_0}[\ln s']_s^{\infty} = \infty$. The integral diverges logarithmically. Formula (3.8) fails the same way: $\frac{1}{4\pi\varepsilon_0}\int_{-\infty}^{\infty}\frac{\lambda\,dz}{\sqrt{s^2+z^2}}$ diverges at large $|z|$.

*The fix.* Choose a finite reference radius $s_0$ with $V(s_0) = 0$:
$$V(s) = -\frac{\lambda}{2\pi\varepsilon_0}\ln\frac{s}{s_0}.$$
Changing $s_0$ adds a constant, which no field or force can detect.

*Why the failure is physical.* Take a finite line of half-length $L$. Then (3.8) converges: $V(s) = \frac{\lambda}{2\pi\varepsilon_0}\operatorname{asinh}(L/s) \approx \frac{\lambda}{2\pi\varepsilon_0}[\ln 2L - \ln s]$ for $L\gg s$. The potential grows without bound as $L\to\infty$, but the $s$-dependence, which is all that matters, converges. With $\lambda = 1$ nC/m ($\lambda/2\pi\varepsilon_0 = 17.975$ V), the project gets

| $L$ | $V(1\text{ cm})$ | $V(10\text{ cm})$ | difference |
|---|---|---|---|
| 1 m | 95.24 V | 53.89 V | 41.345 V |
| 1 km | 219.41 V | 178.02 V | 41.389 V |
| 1000 km | 343.57 V | 302.18 V | 41.389 V |

and the infinite-line prediction for the difference is $17.975\ln10 = 41.389$ V. A "potential" of an infinite object only exists as a difference.

### Example 4.3: The dipole to leading order

Charges $+q$ at $(0,0,d/2)$ and $-q$ at $(0,0,-d/2)$; dipole moment $\mathbf p = qd\,\hat{\mathbf z}$. Find $V$ and $\mathbf E$ for $r\gg d$.

*Potential.* The distances from the charges are $\mathscr r_\pm^2 = r^2 \mp rd\cos\theta + d^2/4$. For $r\gg d$, keeping first order in $d/r$:
$$\frac{1}{\mathscr r_\pm} = \frac1r\left(1 \mp \frac dr\cos\theta + \frac{d^2}{4r^2}\right)^{-1/2} \approx \frac1r\left(1 \pm \frac{d}{2r}\cos\theta\right).$$
So
$$V = \frac{q}{4\pi\varepsilon_0}\left(\frac1{\mathscr r_+} - \frac{1}{\mathscr r_-}\right) \approx \frac{qd\cos\theta}{4\pi\varepsilon_0r^2} = \frac{\mathbf p\cdot\hat{\mathbf r}}{4\pi\varepsilon_0r^2}. \tag{4.1}$$
The monopole term cancels, and the dipole potential falls as $1/r^2$, one power faster than a charge.

*Field in spherical components.* $\mathbf E = -\nabla V$ with $\nabla = \hat{\mathbf r}\partial_r + \hat{\boldsymbol\theta}\frac1r\partial_\theta$:
$$E_r = \frac{2p\cos\theta}{4\pi\varepsilon_0r^3},\qquad E_\theta = \frac{p\sin\theta}{4\pi\varepsilon_0r^3},\qquad E_\phi = 0. \tag{4.2}$$
*Coordinate-free form.* Write $\hat{\mathbf z} = \cos\theta\,\hat{\mathbf r} - \sin\theta\,\hat{\boldsymbol\theta}$, so $\mathbf p = p\cos\theta\,\hat{\mathbf r} - p\sin\theta\,\hat{\boldsymbol\theta}$ and $\mathbf p\cdot\hat{\mathbf r} = p\cos\theta$. Then $3(\mathbf p\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf p = 2p\cos\theta\,\hat{\mathbf r} + p\sin\theta\,\hat{\boldsymbol\theta}$, which is (4.2) times $4\pi\varepsilon_0r^3$:
$$\mathbf E_{\text{dip}} = \frac{1}{4\pi\varepsilon_0r^3}\left[3(\mathbf p\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf p\right]. \tag{4.3}$$
(Exercise B4 derives (4.3) directly in Cartesian coordinates; the project checks it symbolically.) On the axis the field is $2p/(4\pi\varepsilon_0r^3)$ along $\mathbf p$; in the equatorial plane it is half as strong and points against $\mathbf p$.

*How good is "leading order"?* On the axis the exact potential is $\frac{q}{4\pi\varepsilon_0}\left[\frac{1}{r-d/2} - \frac{1}{r+d/2}\right] = \frac{qd}{4\pi\varepsilon_0(r^2 - d^2/4)}$, so $V_{\text{dip}}/V_{\text{exact}} = 1 - d^2/4r^2$: exactly $0.99$ at $r = 5d$. The correction is second order in $d/r$, because the next term with the dipole's symmetry is the octupole (EM-3). The project compares (4.3) with the exact two-charge field in 50 random directions at $r = 20d$: worst relative error $2.2\times10^{-3}$.

### Example 4.4: Energy of a charged sphere, two ways (and a nucleus)

Find the energy of the uniformly charged ball of Example 4.1.

*Way 1, (3.15).* The charge is only inside, where $V = \frac{Q}{8\pi\varepsilon_0R}(3 - r^2/R^2)$:
$$W = \frac12\int_0^R\rho\,\frac{Q}{8\pi\varepsilon_0R}\left(3-\frac{r^2}{R^2}\right)4\pi r^2dr = \frac{\rho Q}{4\varepsilon_0R}\left[R^3 - \frac{R^3}{5}\right] = \frac{\rho QR^2}{5\varepsilon_0}.$$
With $\rho = 3Q/(4\pi R^3)$: $W = \dfrac{3Q^2}{20\pi\varepsilon_0R} = \dfrac35\,\dfrac{Q^2}{4\pi\varepsilon_0R}$.

*Way 2, (3.16).* Split at $R$:
$$W = \frac{\varepsilon_0}{2}\left(\frac{Q}{4\pi\varepsilon_0}\right)^2 4\pi\left[\int_0^R\frac{r^2}{R^6}r^2dr + \int_R^\infty\frac{1}{r^4}r^2dr\right] = \frac{Q^2}{8\pi\varepsilon_0}\left[\frac1{5R} + \frac1R\right] = \frac{3Q^2}{20\pi\varepsilon_0R}.\ \checkmark$$
The interior contributes $\frac{1/5}{6/5} = \frac16$ of the field energy; **five sixths sits outside the ball**, where there is no charge. Way 1 puts all of it inside. Both totals are right. (Exercise C1 gets the same $\tfrac35$ a third way, by building the ball shell by shell.)

Compare a shell of the same $Q$ and $R$: $E = 0$ inside, so only the outer integral survives and $W = \tfrac12\,Q^2/(4\pi\varepsilon_0R)$. Spreading charge through the volume costs $\tfrac35$ instead of $\tfrac12$ because some charge sits closer together. With $Q = 1$ nC, $R = 0.1$ m: $W = 5.3925\times10^{-8}$ J; the project gets this from both integrals to $10^{-10}$.

*From nuclear physics.* A nucleus with $Z$ protons and $A$ nucleons has radius $R \approx r_0A^{1/3}$ with $r_0\approx1.2$ fm. Treat its charge as a uniform ball. In nuclear units $e^2/4\pi\varepsilon_0 = 1.4400$ MeV fm, so the Coulomb energy is
$$E_C = \frac35\frac{e^2}{4\pi\varepsilon_0r_0}\frac{Z(Z-1)}{A^{1/3}} = a_C\frac{Z(Z-1)}{A^{1/3}},\qquad a_C = 0.720\ \text{MeV}.$$
Why $Z(Z-1)$ and not $Z^2$? Each proton does not repel itself: this is (3.14)'s $j\ne i$ restriction, the discrete-versus-continuum distinction, in a real formula. The "Coulomb term" of the semi-empirical mass formula has this form, with fitted coefficients close to $0.7$ MeV. For uranium-238 ($Z = 92$): $E_C = 973$ MeV. Split it into two palladium-119 nuclei ($Z = 46$) and the total is $2\times303 = 606$ MeV. The $367$ MeV difference drives fission; roughly half of it is paid back as the extra surface energy of two smaller drops, and the rest is the famous $\sim200$ MeV per fission.

The project also checks the discrete version: $3000$ equal point charges scattered uniformly through a ball have pair energy (3.14) equal to $0.9968$ of $\tfrac35kQ^2/R$.

### Example 4.5: A parallel-plate capacitor with numbers

Square plates $10$ cm on a side ($A = 0.01$ m²), gap $d = 1$ mm, connected to $V = 100$ V. Ignore edges.

- Capacitance (3.18a): $C = \varepsilon_0A/d = 8.854\times10^{-12}\times0.01/10^{-3} = 88.54$ pF.
- Charge: $Q = CV = 8.854$ nC. Surface density $\sigma = Q/A = 0.8854\ \mu$C/m².
- Field: $E = V/d = 1.0\times10^5$ V/m, and $\sigma/\varepsilon_0 = 0.8854\times10^{-6}/8.854\times10^{-12} = 1.0\times10^5$ V/m. ✓
- Energy (3.19): $W = \tfrac12CV^2 = 4.427\times10^{-7}$ J.
- Force between the plates. Each plate feels the *average* field at its surface (3.17): zero outside, $\sigma/\varepsilon_0$ inside, average $\sigma/2\varepsilon_0$. So $F = \sigma A\cdot\sigma/(2\varepsilon_0) = Q^2/(2\varepsilon_0A) = 4.427\times10^{-4}$ N, attractive, the weight of about 45 mg. The pressure is $\tfrac{\varepsilon_0}{2}E^2 = 0.0443$ Pa.
- Check by virtual work at fixed $Q$ (battery disconnected): $W = Q^2d/(2\varepsilon_0A)$ grows linearly with $d$, so pulling the plates apart takes $F = dW/dd = Q^2/(2\varepsilon_0A) = W/d = 4.427\times10^{-4}$ N. ✓ (At fixed $V$ instead, $W = \tfrac12CV^2$ *decreases* as $d$ grows, yet the force is the same and still attractive: the battery absorbs charge and does work $V\,dQ$, twice the change in $W$. Forgetting the battery gives the wrong sign.)

How good is "ignore edges"? The project solves Poisson's equation for two plates of width $a$ in two dimensions (per unit depth, inside a grounded box of half-width $1.5a$). The field at the centre equals $V/d$, and $\sigma$ on the gap side at the centre equals $\varepsilon_0V/d$, both to $10^{-3}$. But the total capacitance exceeds $\varepsilon_0a/d$:

| $d/a$ | 0.32 | 0.16 | 0.08 | 0.04 |
|---|---|---|---|---|
| $C/C_{\text{ideal}}$ | 1.487 | 1.271 | 1.153 | 1.088 |
| $\sigma_{\text{outer}}/\sigma_{\text{gap}}$ at centre | 0.197 | 0.102 | 0.052 | 0.027 |

The excess comes from charge piling up at the plate edges (where $\sigma$ diverges like (distance)$^{-1/2}$) and from charge on the *outer* faces, which the ideal model sets to zero. Both shrink as $d/a\to0$, but slowly: the excess $C/C_{\text{ideal}} - 1$ falls only somewhat less than proportionally to $d/a$, because of a slowly growing logarithm. For this example's $d/a = 0.01$ the edge correction is a few percent.

### Example 4.6: A point charge in an off-centre cavity

A neutral conducting sphere of radius $R$ has a spherical cavity of radius $a$, not centred. A charge $q$ sits in the cavity at distance $b$ from the cavity's centre. Find the charge on each surface and the field outside.

*Cavity wall.* Take a Gaussian surface inside the metal enclosing the cavity. $\mathbf E = 0$ on it, so $Q_{\text{enc}} = 0$, so the wall carries $-q$. It is not uniform: field lines from $q$ end on the wall, and more of them end on the near side. For a spherical cavity EM-3's image method gives
$$\sigma_{\text{wall}}(\theta) = -\frac{q\,(a^2-b^2)}{4\pi a\,(a^2+b^2-2ab\cos\theta)^{3/2}},$$
with $\theta$ measured from the direction of $q$'s offset. The density is largest at $\theta = 0$ and smallest at $\theta = \pi$, in the ratio $\left(\frac{a+b}{a-b}\right)^3$, which is $27$ for $b = a/2$. The project integrates the formula over the wall ($-q$ to $10^{-10}$) and confirms that the field of $q$ plus its image reproduces it.

*Outer surface.* The sphere is neutral, so the outer surface carries $+q$. How is it distributed? The field inside the metal is zero, so the charges on the outer surface must produce, together with everything else, zero field in the metal. The combination "$q$ plus the wall charge" already produces zero field everywhere outside the cavity: their total is zero, and by Theorem 3.4 the wall charge arranges itself to cancel $q$'s field in the metal, and a field that vanishes in the metal and has no charge beyond the wall vanishes beyond it too. That leaves the outer charge to produce zero field inside radius $R$ on its own and make the outer surface an equipotential. A uniform $\sigma = q/(4\pi R^2)$ does both (the shell of §3.4). By the uniqueness theorem (§3.11, proved in EM-3) it is the only distribution that does.

*Field outside.* A uniform shell of total charge $q$ produces
$$\mathbf E = \frac{q}{4\pi\varepsilon_0r^2}\hat{\mathbf r}\qquad (r > R),$$
with $r$ measured from the **centre of the conducting sphere**, not from the charge or the cavity. Outside, you cannot tell where the charge is or what shape the cavity has. Numbers: $q = 1$ nC, $R = 10$ cm gives $\sigma_{\text{outer}} = 7.96$ nC/m², uniform, and the conductor sits at $V = q/(4\pi\varepsilon_0R) = 89.88$ V. If the outer surface were not a sphere, the outer charge would not be uniform, but it would still be independent of where $q$ sits inside.

## 5. Why it matters and where it connects

- **EM-3** takes Poisson's equation (3.9) with conductors as boundary conditions. Its uniqueness theorems justify Theorem 3.5, §3.11, and Example 4.6; the method of images computes the induced charge of §3.9 exactly; the multipole expansion generalizes Example 4.3 to all orders; the conductor in a uniform field (project Part 6) is solved there in closed form.
- **EM-4** adds matter: polarization produces bound charge, and the boundary conditions (3.10)–(3.12) return with free charge in place of total charge.
- **EM-5** mirrors this module for magnetism: $\nabla\cdot\mathbf B = 0$ (no magnetic Gauss charge) and Ampère's law in place of $\nabla\times\mathbf E = 0$.
- **EM-6** makes $\tfrac{\varepsilon_0}{2}E^2$ part of the energy density $u = \tfrac{\varepsilon_0}{2}E^2 + \tfrac{1}{2\mu_0}B^2$ that flows with the Poynting vector, settling "where is the energy" in favour of the field; the surface pressure (3.17) becomes the Maxwell stress tensor.
- **LA-13** is the numerical engine of project Part 5: the five-point Laplacian, Jacobi and Gauss–Seidel (LA-13 Theorem 3.19 for when they converge; (3.24) for how slowly), and conjugate gradients.
- **CM-5** is the gravitational twin: the shell theorem is Gauss's law with $-G$ in place of $1/4\pi\varepsilon_0$, and the $1/r$ potential gives Kepler orbits.
- **QM-6** puts the Coulomb potential $-e^2/(4\pi\varepsilon_0r)$ into the Schrödinger equation for hydrogen. The finite size of the proton enters through exactly Example 4.1's interior potential $\propto(3-r^2/R^2)$.
- **F6** (Dirichlet principle) and **F3** (Green's functions as impulse responses) are the mathematical frames for §3.12 and (3.8).

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: superposition, flux through cube faces, a charged slab, a boundary-condition jump, a capacitor network, the energy of four charges. B1–B8 are the module: B1 (Gauss's law for gravity, the tunnel through the Earth), B3 (the shell potential by direct integration), B5 (why energies do not superpose, and the interaction term), B6 (surface force on a sphere), and B7 (conductor with two cavities) are essential; B2(d), B4(c), and B8 use `project.py`. C1–C3 go further: the self-energy and the classical electron radius, relaxation as diffusion (code), and Thomson's theorem.

## 7. Coding project

`project.py` (about 15 s) builds:

1. **Symbolic checks** with SymPy: Example 4.1's $V$ solves (3.9) with continuous $V$ and $dV/dr$; the energy $\tfrac35kQ^2/R$ three ways with interior share $\tfrac16$; (4.3) from $-\nabla$ of (4.1); the capacitor limits of Exercise B2; the shell potential of B3.
2. **Point charges.** `field_potential` evaluates (3.2) and the matching $V$ by broadcasting. Field lines and equipotentials of a dipole and a quadrupole (`dipole_quadrupole.png`); Example 4.3's numbers; $\oint\mathbf E\cdot d\mathbf l = 0$; a dipole field line traced numerically stays on $r = r_0\sin^2\theta$ to $4\times10^{-5}$.
3. **Gauss's-law checker.** Flux through spheres (Gauss–Legendre in $\cos\theta$) and cubes (Gauss–Legendre on each face) around twelve random charges, equal to $Q_{\text{enc}}/\varepsilon_0$ to $10^{-7}$; the corner-charge flux of A2.
4. **The charged sphere and the line** (Examples 4.1, 4.2, 4.4): numerical inward integration of $E$; energy both ways; the point-charge pair sum; self-energy; nuclear numbers.
5. **Poisson on a grid.** `solve_dirichlet` assembles the five-point Laplacian with conductors as fixed nodes and solves with `scipy.sparse`. The parallel-plate capacitor with fringing (`capacitor_fringing.png`), $C$ from energy equal to $C$ from charge (a discrete Green's identity: (3.15) and (3.16) agree to roundoff on the grid), and a Jacobi run whose sweep count quadruples ($3279\to13124$) when $h$ halves.
6. **Relaxation around a conductor.** Red-black SOR for a grounded cylinder in a uniform field (`conductor_uniform_field.png`): $\sigma$ from the discrete charge at the conductor's nodes is $0.989\times2\varepsilon_0E_0\cos\phi$; the ratio of tangential to normal field on circles just outside the conductor falls ($0.189, 0.105, 0.061$ at distances $0.2a, 0.1a, 0.05a$) toward zero, as Theorem 3.4(v) requires.
7. **The cavity** of Example 4.6, and every number in `solutions.md`.

**Build it yourself:** before reading the code, write `field_potential` (no Python loop over charges) and `solve_dirichlet` (the sparse matrix for the free nodes, with fixed neighbours moved to the right-hand side). Then compare. **Your turn** asks for the ring and disk on axis, a coaxial capacitor's $C$ against (3.18c), and a numerical demonstration that the field inside a charged shell vanishes; reference solutions follow the prompts.

## 8. Where you'll struggle

- **Using Gauss's law where it is true but useless.** Everyone tries it on a finite rod or a charged square. The law holds; the trouble is that $|\mathbf E|$ varies over every closed surface you can draw, so you cannot pull it out of the integral. Before writing a Gaussian surface, write the symmetry argument that fixes the direction of $\mathbf E$ and what it depends on (§3.3, the $180°$ rotation). If you cannot, use (3.3) or (3.8) instead.
- **The $\tfrac12$ in the energy.** Derive it every time from "each pair once" (§3.8) until it is automatic. The test: in $\tfrac12\sum q_iV(\mathbf r_i)$, can you say whether $V(\mathbf r_i)$ includes $q_i$'s own potential? (Discrete: no. Continuous: it doesn't matter.) The same $\tfrac12$ trips people in image-charge energies in EM-3.
- **The self-energy infinity is real.** Do not convince yourself that it is a mathematical trick that goes away. It is a genuine failure of classical point charges, unresolved within classical physics. The practical rule: use (3.14) for point charges, (3.15)/(3.16) for continuous distributions, and never mix a point charge's own field into (3.16).
- **"Where is the energy?"** (3.15) puts it where the charge is, (3.16) where the field is. People decide one must be wrong. In statics both are correct bookkeepings of the same total and no measurement distinguishes them. The field picture is forced only by dynamics (EM-6: radiation carries energy through empty space).
- **Infinite distributions.** $V(\infty) = 0$ fails for infinite lines and planes because they put charge at infinity. Pick a finite reference point and work only with differences (Example 4.2). The symptom is a $\ln$ or a linear term that will not go to zero.
- **Boundary-condition signs.** (3.10) is "above minus below, $\hat{\mathbf n}$ pointing from below to above," and (3.12) has a minus sign. Check every use against the infinite plane: $E = +\sigma/2\varepsilon_0$ above, $-\sigma/2\varepsilon_0$ below.
- **"$E = 0$ so $V = 0$."** Inside a shell or a conductor $E = 0$ and $V$ is *constant*, equal to its surface value. Zero slope, not zero height.
- **Finite-difference Poisson solvers are slow.** The project's Jacobi sweep count goes up fourfold when the grid spacing halves. That is not a bug: the smooth error modes decay at rate $\cos(\pi h)$ per sweep (LA-13 (3.24)). Use a sparse direct solver or conjugate gradients for real work, and SOR when you want relaxation to be tolerable.

## 9. Check yourself

1. State Gauss's law in integral and differential form. Which step of its derivation uses the inverse-square law?
2. Why does Gauss's law plus spherical symmetry give Coulomb's law, but Gauss's law alone does not determine $\mathbf E$?
3. Name the three symmetries for which Gauss's law determines $\mathbf E$, and the Gaussian surface for each.
4. Derive $\nabla\times\mathbf E = 0$ for a point charge in one line.
5. Why does $V(\infty) = 0$ fail for an infinite line charge, and what do you do instead?
6. Write the boundary conditions on $\mathbf E$ and on $\partial V/\partial n$ at a surface charge, with the direction of $\hat{\mathbf n}$.
7. Where does the $\tfrac12$ in $W = \tfrac12\sum q_iV(\mathbf r_i)$ come from? Does $V(\mathbf r_i)$ include $q_i$?
8. Derive $W = \tfrac{\varepsilon_0}{2}\int E^2\,d\tau$ from $W = \tfrac12\int\rho V\,d\tau$. Why does the surface term vanish?
9. What fraction of a uniformly charged ball's field energy lies inside the ball?
10. List the five properties of a conductor in equilibrium. Which one gives $\sigma = -\varepsilon_0\partial V/\partial n$?
11. What is the force per unit area on a conductor's surface, and why is it not $\sigma\cdot\sigma/\varepsilon_0$?
12. A charge sits off-centre in a cavity inside a neutral conducting sphere. What field does an observer outside measure, and centred where?

## 10. Pair with

- **Griffiths, *Introduction to Electrodynamics*, ch. 2.** The standard treatment, in the same order and notation (script $\mathscr r$). §2.4–2.5 on energy and conductors are especially clear; do its problems for drill, since these notes have fewer of the routine Gauss's-law geometries.
- **Purcell & Morin, *Electricity and Magnetism*, ch. 1–3.** Better physical pictures: field lines, the energy of a charge distribution built up piece by piece, and conductors. Its exercises on the surface force and on charges near conductors are excellent.
- **MIT 8.07 lecture notes (Guth), lectures 1–4** (OCW). A more mathematical pass at the delta function, Gauss, the potential, and energy, one level above Griffiths; read after the notes as a second derivation.
- **Feynman Lectures on Physics, vol. II, ch. 4–8.** Ch. 4–6 for intuition on fields of standard configurations; ch. 8 on electrostatic energy and where it is located. Ch. 28 (electromagnetic mass) is the most honest short account of the self-energy problem.
