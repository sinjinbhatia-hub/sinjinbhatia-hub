# SM-1: Thermodynamics I: State Variables, Temperature, the First Law, the Ideal Gas, Heat Capacities, Processes

**Time:** ~8 hours · **Prerequisites:** multivariable calculus (partial derivatives, chain rule, line integrals) · **Unlocks:** SM-2

> **Conventions in this module.** SI units throughout. $k_B = 1.380649\times10^{-23}$ J/K (Boltzmann's constant) is always written explicitly; $N$ is the number of particles, $n$ the number of moles, $N_A = 6.02214076\times10^{23}$ mol⁻¹, $R = N_Ak_B = 8.314$ J/(mol·K), so $Nk_B = nR$. State functions are capitals: $U$ (internal energy), $H$ (enthalpy), $P$, $V$, $T$. **First law: $dU = \delta Q + \delta W$, with $\delta Q$ the heat added TO the system and $\delta W$ the work done ON the system** (Schroeder's sign). Compression work is $\delta W = -P\,dV$. Engineering texts (Çengel & Boles) write $\Delta U = Q - W$ with $W$ the work done BY the system; their $W$ is minus ours. A $\delta$ (or đ) marks an inexact differential: a small amount of something that is not the change of any function. $(\partial x/\partial y)_z$ means "$z$ held fixed." $f$ is the number of quadratic degrees of freedom per molecule; $\gamma = C_P/C_V$. In 2D, "volume" is area $A$ and "pressure" is force per unit length.

## 0. Where this sits

You have taken an engineering thermodynamics course, so most of the formulas here will look familiar. This module redoes that material from a physicist's angle, for one purpose: to make SM-3's statistical derivation of thermodynamics land. That needs three things you may not have been given: a microscopic picture of $P$ and $T$ (kinetic theory), a precise statement of why $U$ is a state function and $Q$ is not (exact vs inexact differentials), and fluency with partial derivatives where the variable held fixed changes from line to line.

SM-2 adds entropy and the second law to the first law developed here. SM-3 derives $T$ from counting microstates. SM-4 proves the equipartition theorem stated here and derives the Maxwell distribution. SM-5 derives the Stefan–Boltzmann law and fixes the heat-capacity failures flagged in §3.3. SM-6 uses the van der Waals equation of §3.11 for the liquid–gas transition.

## 1. The problem that forced this idea

Steam engines came before the theory of steam engines. Newcomen's atmospheric engine pumped water out of mines from 1712. James Watt patented the separate condenser in 1769, and Boulton & Watt engines were in commercial service by 1776. Engineers improved them by trial and error. Nobody could say how much work a given amount of coal could possibly deliver.

Sadi Carnot asked exactly that in 1824 (*Réflexions sur la puissance motrice du feu*). He believed heat was a fluid, "caloric," that flowed from hot to cold without being used up, the way water drives a mill wheel without being consumed. On that wrong picture he still got the right answer: the maximum efficiency depends only on the two temperatures. That result belongs to SM-2. What matters here is that Carnot could not say what heat *was*.

The caloric picture had a known problem. Count Rumford (1798) reported that boring cannon barrels produced heat apparently without limit, as long as the horses kept turning the drill. Julius Robert Mayer (1842) argued that heat and work are interconvertible. James Joule spent the 1840s measuring it: a falling weight drove a paddle wheel in water, and the water warmed by an amount proportional to the work done. Heat is energy in transit. Rudolf Clausius (1850) and William Thomson, later Lord Kelvin (1851), reconciled Carnot with Joule and stated the two laws of thermodynamics in essentially their modern form.

Here is the problem that this module is about. Energy is conserved. Yet "heat" and "work" are not properties of a system. You cannot say how much heat a gas *contains*; you can only say how much heat flowed in during a particular process. Take a gas from state 1 to state 2 by two different routes and the heat absorbed differs, the work done differs, but the difference of the two is the same. The first law's content is precisely this: $Q$ and $W$ separately depend on the path; $Q + W$ does not, so it is the change of a function $U$ of the state. Making that sentence precise needs the mathematics of exact and inexact differentials (§3.6), and making $U$, $P$, and $T$ meaningful needs a picture of molecules (§3.2).

## 2. Intuition

### A box of molecules

A litre of air holds about $2.7\times10^{22}$ molecules. They fly in straight lines at a few hundred metres per second, bounce off each other every hundred nanometres or so, and bounce off the walls.

**Pressure** is what the walls feel. Each molecule that hits a wall reverses its perpendicular momentum, delivering an impulse $2mv_x$. Many impulses per unit time per unit area make a steady force per area. Double the number of molecules in the box: twice as many hits, twice the pressure. Double the molecular speed: each hit delivers twice the impulse *and* hits come twice as often, so four times the pressure. Pressure goes as (number density) × $m$ × (speed)².

**Temperature** is what is equal between two bodies once they have stopped exchanging energy. When fast molecules collide with slow ones, on average energy passes from the fast to the slow. Exchange stops when the average kinetic energy per molecule is the same on both sides. So temperature should measure the average kinetic energy per molecule.

Put those together: $PV \propto N\langle\tfrac12 mv^2\rangle$, and the empirical gas law $PV = Nk_BT$ then says $\langle\tfrac12 mv^2\rangle \propto k_BT$. You can guess the constant from the counting: there are three directions; the pressure on one wall involves only $v_x$; so expect $PV = N m\langle v_x^2\rangle = \tfrac13 Nm\langle v^2\rangle$, and therefore $\langle\tfrac12 mv^2\rangle = \tfrac32 k_BT$, or $\tfrac12k_BT$ per direction of motion.

### Heat versus work

Energy can enter the gas two ways.

- **Work**: a wall moves. A piston pushes in, every molecule that hits it bounces back faster (a ball hitting an advancing bat). The energy transfer is organized: one macroscopic coordinate (the piston position) changes, and the force on it is $PA$.
- **Heat**: no wall moves, but the wall is hotter than the gas. Wall atoms jiggle randomly and molecules leave the wall faster on average than they arrived. The energy transfer is disorganized: no macroscopic coordinate changes.

Both raise the molecules' kinetic energy. Once the energy is inside, it is just $U$, and there is no record of which route it took. That is why "the heat in the gas" is meaningless and "the energy in the gas" is not.

### States as points, processes as paths

For a fixed amount of gas in equilibrium, two numbers fix everything: say $(P, V)$. The state is a point on the $P$–$V$ plane. A slow process is a curve. The work done on the gas along a curve is minus the area under it:

```
   P
   |   2 *
   |    / \      path b: up, then right
   |   /   |     path a: right, then up
   |  /    |
   |1*-----+
   +-------------- V
```

Go from 1 to 2 by two routes and the area under each route differs, so $W$ differs. But $U$ is a function of where you are: $U(2) - U(1)$ is fixed. So $Q = \Delta U - W$ must differ too, by exactly the amount $W$ does.

**The guesses you should now be able to make.** (i) $PV = \tfrac23 N\langle\tfrac12mv^2\rangle$. (ii) For a monatomic gas, $U = \tfrac32 Nk_BT$, so heating at constant volume costs $\tfrac32 Nk_B$ per kelvin. (iii) Heating at constant pressure costs more, because the gas also pushes the atmosphere back: the extra should be $P\,\Delta V/\Delta T = Nk_B$. (iv) Squeeze a gas fast enough that no heat escapes and it heats up, so $P$ rises faster than the isothermal $1/V$.

## 3. Formal development

### 3.1 Systems, state variables, the zeroth law

A **system** is the part of the world we track; the rest is the **surroundings**. A system is in **equilibrium** when its macroscopic properties have stopped changing. In equilibrium a simple fluid of fixed $N$ is described by a few **state variables** (here $P$, $V$, $T$), related by an **equation of state** $F(P, V, T) = 0$. Two of them are independent.

State variables split into **extensive** (scale with system size: $V$, $N$, $U$) and **intensive** (do not: $P$, $T$).

**Zeroth law.** If $A$ is in thermal equilibrium with $C$, and $B$ is in thermal equilibrium with $C$, then $A$ is in thermal equilibrium with $B$.

This is what makes a thermometer possible: $C$ is the thermometer, and the zeroth law guarantees that "reads the same on $C$" is a transitive relation, so it sorts all bodies into classes. Temperature is a label on those classes. Any monotone property of $C$ (mercury length, resistance, gas pressure) gives a valid label; the problem is that different thermometers disagree between their calibration points.

**The ideal-gas scale.** For any gas at low enough density, $PV/N$ at fixed temperature is the same for all gases. So define $T$ by $T \propto \lim_{\rho\to0} PV/N$. This scale does not depend on which gas you use. Its zero is where $PV/N \to 0$, and there are no negative temperatures on it. Since 2019 the kelvin is defined by fixing $k_B = 1.380649\times10^{-23}$ J/K exactly, which makes the scale

$$PV = Nk_BT = nRT. \tag{3.1}$$

The thermodynamic (absolute) scale that Kelvin defined in 1848 from Carnot's engine coincides with this one; SM-2 proves it.

### 3.2 Kinetic theory: pressure from momentum transfer

Daniel Bernoulli derived Boyle's law from colliding particles in 1738, a century early. John Waterston's 1845 kinetic-theory paper was refused by the Royal Society and published only in 1892 by Rayleigh. Clausius (1857) and Maxwell (1860) made kinetic theory standard. Here is the derivation.

**Model.** $N$ molecules of mass $m$ in volume $V$, number density $\rho_N = N/V$. Velocities have a distribution: a fraction $g(\mathbf v)\,d^3v$ of molecules has velocity in $d^3v$ around $\mathbf v$, the same everywhere in the box, and isotropic ($g$ depends only on $|\mathbf v|$). Walls reflect molecules elastically on average (if they don't, the result is the same once the wall and gas are in equilibrium, because the outgoing flux then mirrors the incoming one).

▶ **Try it first.** Take a patch of wall of area $A$ perpendicular to $x$. In time $dt$, how many molecules with velocity component $v_x > 0$ hit it? (Draw the slanted cylinder they must start in.) Each delivers what impulse?

**Derivation.**

1. A molecule with velocity $\mathbf v$, $v_x > 0$, hits the patch in the next $dt$ if it starts inside a slanted cylinder with base $A$ and perpendicular height $v_x\,dt$. That cylinder has volume $Av_x\,dt$.
2. The number of such molecules is $\rho_N\,g(\mathbf v)\,d^3v \cdot A v_x\,dt$.
3. Each reverses $v_x$, so the wall receives impulse $2mv_x$ per molecule.
4. Total impulse in $dt$: $\displaystyle dp = \int_{v_x>0} 2mv_x \cdot \rho_N g(\mathbf v) A v_x\,dt\,d^3v = 2m\rho_N A\,dt\int_{v_x>0} v_x^2 g\,d^3v$.
5. By the symmetry $g(v_x,\dots) = g(-v_x,\dots)$, the half-space integral is half the full one: $\int_{v_x>0} v_x^2 g\,d^3v = \tfrac12\langle v_x^2\rangle$.
6. Pressure is impulse per time per area: $P = dp/(A\,dt) = \rho_N m\langle v_x^2\rangle$.
7. Isotropy: $\langle v_x^2\rangle = \langle v_y^2\rangle = \langle v_z^2\rangle = \tfrac13\langle v^2\rangle$. So

$$P = \frac13\,\frac{N}{V}\,m\langle v^2\rangle, \qquad PV = \frac23 N\left\langle\tfrac12 mv^2\right\rangle. \tag{3.2}$$

Collisions between molecules never entered. They redistribute velocities but, in equilibrium, leave $g$ unchanged, and step 2 needs only $g$. ∎

Now compare (3.2) with the gas law (3.1):

$$\left\langle\tfrac12 mv^2\right\rangle = \tfrac32 k_BT, \qquad \left\langle\tfrac12 mv_x^2\right\rangle = \tfrac12 k_BT. \tag{3.3}$$

■ **What just happened.** Pressure is momentum flux. Combining the mechanical result (3.2) with the empirical thermometer (3.1) tells you what temperature *is* for an ideal gas: $\tfrac23$ of the mean translational kinetic energy per molecule, divided by $k_B$. $k_B$ is a units conversion between kelvins and joules. $R = N_Ak_B$ is the same conversion per mole.

**In $d$ dimensions** step 7 becomes $\langle v_x^2\rangle = \langle v^2\rangle/d$, so $PV = \tfrac{2}{d}N\langle\tfrac12mv^2\rangle$ and $\langle\tfrac12mv^2\rangle = \tfrac d2 k_BT$. In 2D (the project's gas), $PA = Nk_BT$ with $\langle\tfrac12mv^2\rangle = k_BT$.

**The Maxwell speed distribution** (stated here; derived in SM-4 from the Boltzmann factor). The fraction of molecules with speed in $[v, v + dv]$ is $\mathcal D(v)\,dv$ with

$$\mathcal D_{3D}(v) = 4\pi\left(\frac{m}{2\pi k_BT}\right)^{3/2} v^2 e^{-mv^2/2k_BT}, \qquad \mathcal D_{2D}(v) = \frac{m}{k_BT}\,v\,e^{-mv^2/2k_BT}. \tag{3.4}$$

The prefactor $v^2$ (or $v$ in 2D) is the area of the sphere (circumference of the circle) of radius $v$ in velocity space; each Cartesian component is Gaussian with variance $k_BT/m$. Three speeds you will quote: the most probable $v_p = \sqrt{2k_BT/m}$, the mean $\bar v = \sqrt{8k_BT/\pi m}$, the rms $v_{\rm rms} = \sqrt{3k_BT/m}$ (3D). The 2D cumulative distribution is $F(v) = 1 - e^{-mv^2/2k_BT}$, which the project uses.

### 3.3 Equipartition

(3.3) says each translational direction carries $\tfrac12k_BT$. That is a special case of a general rule.

**Theorem 3.1 (Equipartition; stated, proved in SM-4).** For a classical system in thermal equilibrium at temperature $T$, every degree of freedom that enters the energy only as a quadratic term ($\tfrac12mv_x^2$, $\tfrac12I\omega^2$, $\tfrac12\kappa x^2$) has mean energy $\tfrac12k_BT$.

So if each molecule has $f$ quadratic degrees of freedom and there are no other energy contributions,

$$U = \frac f2\,Nk_BT. \tag{3.5}$$

Counting $f$:

| System | Quadratic terms | $f$ | $C_V/(Nk_B)$ predicted |
|---|---|---|---|
| Monatomic gas (He, Ar) | 3 translations | 3 | 3/2 |
| Diatomic gas (N₂, O₂), rigid | 3 translations + 2 rotations (about the two axes perpendicular to the bond) | 5 | 5/2 |
| Diatomic, vibrating | + 1 kinetic + 1 potential for the bond stretch | 7 | 7/2 |
| Solid, atoms on springs | 3 kinetic + 3 potential per atom | 6 | 3 (Dulong–Petit) |

Dulong and Petit observed in 1819 that many solid elements have molar heat capacity near $3R \approx 25$ J/(mol·K). Equipartition explains it.

**The failures, which are not footnotes.** Air at room temperature has $C_V \approx \tfrac52 R$, not $\tfrac72 R$: the vibration contributes nothing. Hydrogen gas has $C_V \approx \tfrac32R$ below about 100 K, $\tfrac52R$ near room temperature, and climbs toward $\tfrac72R$ only at thousands of kelvin. Diamond at room temperature has a molar heat capacity about a quarter of $3R$, and every solid's $C_V \to 0$ as $T\to0$. A classical degree of freedom cannot "switch off": Theorem 3.1 has no temperature dependence in it. Maxwell called the specific-heat problem the greatest difficulty facing molecular theory (1875). Kelvin's 1900 lecture on "two clouds" over nineteenth-century physics named equipartition's failure as one of them. The resolution is that energy levels are discrete: a mode whose level spacing exceeds $k_BT$ is frozen out. Einstein's 1907 model of a solid did this first; SM-3 counts its states, SM-5 does the Debye refinement. This is the first place in the curriculum where classical physics gives a wrong number and quantum mechanics is the fix.

### 3.4 Exact and inexact differentials

Everything about "state functions" reduces to one theorem from multivariable calculus.

A **differential form** in two variables is an expression $\omega = M(x,y)\,dx + N(x,y)\,dy$. Integrated along a path $C$ it gives $\int_C\omega = \int (M\dot x + N\dot y)\,dt$. It is **exact** if $\omega = dF$ for some function $F$, i.e., $M = \partial F/\partial x$ and $N = \partial F/\partial y$. Then $\int_C\omega = F(\text{end}) - F(\text{start})$ for every path $C$: path independent. If $\omega$ is not exact we write it as $\delta(\cdot)$ and its integral depends on the path.

**Theorem 3.2 (Exactness test).** Let $M, N$ be $C^1$ on a simply connected region $D$ of the $(x,y)$ plane. Then $\omega = M\,dx + N\,dy$ is exact on $D$ if and only if

$$\frac{\partial M}{\partial y} = \frac{\partial N}{\partial x}\quad\text{on } D. \tag{3.6}$$

▶ **Try it first.** One direction is Clairaut's theorem (equal mixed partials). For the other, think of $(M, N, 0)$ as a vector field in 3D. What is its curl?

**Derivation.**

1. (Exact ⇒ (3.6).) If $M = F_x$ and $N = F_y$ with $F$ twice continuously differentiable, then $M_y = F_{xy} = F_{yx} = N_x$ by Clairaut.
2. ((3.6) ⇒ exact.) Set $\mathbf F = (M, N, 0)$. Its curl is $(0, 0, N_x - M_y) = \mathbf 0$ by (3.6).
3. EM-1 Theorem 3.16 (curl-free fields on a simply connected region are gradients) gives $\Phi$ with $\nabla\Phi = \mathbf F$ (the theorem's $\mathbf F = -\nabla U$; set $\Phi = -U$). In particular $\Phi_x = M$, $\Phi_y = N$, and $\Phi$ does not depend on $z$ because $\Phi_z = 0$.
4. So $\omega = d\Phi$. Concretely, $\Phi(x,y) = \int_{(x_0,y_0)}^{(x,y)}\omega$ along any path in $D$, path independent by Stokes. ∎

■ **What just happened.** "Is there a function whose changes this is?" is the same question as "is this force conservative?" from CM-1 and EM-1. The test is mixed partials, and you can run it without ever finding the function.

In thermodynamics we will read Theorem 3.2 both ways. If a quantity is a function of the state, its differential passes test (3.6). If a differential fails (3.6), no state function has it as its change, and its integral depends on the path.

**The test on an ideal gas.** Use $(T, V)$ as coordinates for $N$ molecules with $U = C_VT$ ($C_V$ constant). The three differentials, with $\delta W = -P\,dV$ (justified in §3.5) and $\delta Q = dU - \delta W$:

- $dU = C_V\,dT + 0\,dV$. $\partial C_V/\partial V = 0 = \partial 0/\partial T$. Exact. ✓.
- $\delta W = 0\,dT - \dfrac{Nk_BT}{V}\,dV$. $\partial 0/\partial V = 0$ but $\partial(-Nk_BT/V)/\partial T = -Nk_B/V \ne 0$. **Inexact.**
- $\delta Q = C_V\,dT + \dfrac{Nk_BT}{V}\,dV$. $\partial C_V/\partial V = 0$ but $\partial(Nk_BT/V)/\partial T = Nk_B/V \ne 0$. **Inexact.**

Divide $\delta Q$ by $T$: $\delta Q/T = (C_V/T)\,dT + (Nk_B/V)\,dV$. Now both cross-derivatives are zero. **Exact.** The function is $C_V\ln T + Nk_B\ln V$ + const. That is the entropy of an ideal gas; $1/T$ is an integrating factor for heat. SM-2 shows this works for every system, and that is the second law. Keep this in your pocket.

### 3.5 The first law and compression work

**First law.** For any process between equilibrium states, there is a state function $U$, the internal energy, such that

$$\Delta U = Q + W, \qquad dU = \delta Q + \delta W, \tag{3.7}$$

with $Q$ the heat added to the system and $W$ the work done on it.

The physics content: $U$ exists as a state function, even though $Q$ and $W$ are not. Joule's experiments are the evidence: insulate a system (so $Q = 0$) and the work needed to take it between two given states is the same however you do the work (paddle wheel, electric heater, compression). That defines $\Delta U$ adiabatically; then (3.7) defines $Q$ for non-adiabatic processes.

> **Sign conventions (boxed, and repeated where used).**
> Physics (this curriculum, Schroeder, Kittel): $\Delta U = Q + W$, $W$ = work done **on** the system, $W = -\int P\,dV$.
> Engineering (Çengel & Boles, most ME courses): $\Delta U = Q - W$, $W$ = work done **by** the system, $W_b = +\int P\,dV$ ("boundary work").
> Map: $W_{\rm ours} = -W_{\rm Çengel}$. $Q$ is the same in both (heat **into** the system is positive). When reading an engineering text, translate every $W$ at the moment you read it; never mix in one calculation.

**Compression work for a quasistatic process.** A process is **quasistatic** if it is slow enough that the system passes through a continuous sequence of equilibrium states, so $P$ and $T$ are defined at every instant.

▶ **Try it first.** A gas in a cylinder of cross-section $A$ is pushed by a piston. The piston moves a distance $dx$ inward. What force does the piston exert on the gas, and how much work does it do? Express in terms of $dV$.

**Derivation.**

1. Let $x$ be the piston's position, with the gas occupying $V = Ax$. In a quasistatic process the gas pressure is uniform and equal to $P$, and the piston exerts force $PA$ on the gas, directed into the gas (toward decreasing $x$).
2. The piston moves from $x$ to $x + dx$. The work done on the gas is force times displacement along the force: $\delta W = (PA)(-dx)$, since the force points in $-x$.
3. $A\,dx = dV$. So $\delta W = -P\,dV$, and for a finite process

$$W = -\int_{V_i}^{V_f} P\,dV. \tag{3.8}$$

4. Compression ($dV < 0$) gives $W > 0$: work is done on the gas. Expansion gives $W < 0$. ∎

■ **What just happened.** Work is minus the area under the path on the $P$–$V$ diagram, taken in the direction of travel. Different paths have different areas, which is the inexactness of $\delta W$ seen geometrically. Around a closed cycle, $-\oint P\,dV$ is the enclosed area with a sign: clockwise cycles on a $P$–$V$ diagram do net work **on** the surroundings ($W < 0$), which is what an engine does.

**When (3.8) fails.** If the piston moves fast, the gas near it is compressed more than the gas far away; $P$ is not defined for the gas as a whole, and the force on the piston is not $P_{\rm eq}A$. In a free expansion into vacuum, nothing pushes back, so $W = 0$ even though $V$ changes. (3.8) is for quasistatic processes only. The project makes this concrete by moving its piston fast.

### 3.6 Heat capacities, enthalpy, and $C_P - C_V$

The heat capacity is the heat needed per unit temperature rise. Because $\delta Q$ is inexact, this depends on the path: you must say what is held fixed.

**At constant volume**, $\delta W = 0$, so $\delta Q = dU$ and

$$C_V = \left(\frac{\delta Q}{dT}\right)_V = \left(\frac{\partial U}{\partial T}\right)_V. \tag{3.9}$$

**At constant pressure**, $\delta Q = dU + P\,dV$. Since $P$ is fixed, $P\,dV = d(PV)$, so $\delta Q = d(U + PV)$. Define the **enthalpy**

$$H = U + PV, \qquad C_P = \left(\frac{\partial H}{\partial T}\right)_P. \tag{3.10}$$

At constant pressure, heat added equals the change in $H$. That is why chemists and engineers tabulate enthalpies: lab beakers and open flow systems sit at atmospheric pressure. $H$ is "energy plus the work you had to do to make room for the system against constant pressure." SM-2 shows $H$ is a Legendre transform of $U$.

For an ideal gas, $U$ depends on $T$ only. Joule's free-expansion experiment (1845) found no temperature change when a gas expanded into vacuum, so $(\partial U/\partial V)_T = 0$ to his accuracy; kinetic theory says the same thing, since (3.5) contains no $V$. Real gases show a small effect (SM-2's Joule–Thomson coefficient).

▶ **Try it first.** Write $U = U(T, V)$ and expand $dU$. Insert it in $\delta Q = dU + P\,dV$, then hold $P$ fixed and divide by $dT$.

**Derivation of $C_P - C_V$.**

1. $dU = (\partial U/\partial T)_V\,dT + (\partial U/\partial V)_T\,dV = C_V\,dT + (\partial U/\partial V)_T\,dV$.
2. $\delta Q = dU + P\,dV = C_V\,dT + \left[(\partial U/\partial V)_T + P\right]dV$.
3. Along a constant-$P$ path, $dV = (\partial V/\partial T)_P\,dT$. Divide by $dT$:
$$C_P = C_V + \left[\left(\frac{\partial U}{\partial V}\right)_T + P\right]\left(\frac{\partial V}{\partial T}\right)_P. \tag{3.11}$$
4. Ideal gas: $(\partial U/\partial V)_T = 0$ and $V = Nk_BT/P$ gives $(\partial V/\partial T)_P = Nk_B/P$. So

$$C_P - C_V = Nk_B = nR. \tag{3.12}$$

∎

■ **What just happened.** At constant pressure, part of the heat goes into pushing the surroundings back. For an ideal gas that extra is $P\,\Delta V = Nk_B\,\Delta T$, exactly. (3.11) holds for any substance; the bracket is "internal pressure plus external pressure." SM-2 rewrites it as $TV\beta^2/\kappa_T$ using a Maxwell relation.

With (3.5), $C_V = \tfrac f2Nk_B$, $C_P = \tfrac{f+2}{2}Nk_B$, and

$$\gamma \equiv \frac{C_P}{C_V} = \frac{f+2}{f}: \qquad \gamma = \tfrac53\ \text{(monatomic)}, \quad \tfrac75 = 1.4\ \text{(rigid diatomic)}, \quad 2\ \text{(2D monatomic, } f=2). \tag{3.13}$$

### 3.7 Latent heat

At a phase transition (boiling, melting) heat goes in and $T$ does not change. The **latent heat** $L$ is the heat per unit mass needed to convert the phase completely, normally quoted at constant pressure, so $L = \Delta h$ per unit mass (enthalpy, not energy). For water at 1 atm and 100 °C, $L_{\rm vap} = 2.257$ MJ/kg.

How much of that is energy and how much is pushing the atmosphere back? The vapour from 1 kg of water occupies (ideal gas) $V = (m/M)RT/P$; the liquid's volume is negligible. The work done by the vapour on the atmosphere is $P\,\Delta V \approx (m/M)RT = (1/0.018015)(8.314)(373.15) = 0.172$ MJ. So $\Delta U = 2.257 - 0.172 = 2.085$ MJ: about 92% of the latent heat breaks intermolecular bonds and 8% makes room. The Clausius–Clapeyron equation that relates $L$ to the slope of the boiling curve needs entropy (SM-2).

### 3.8 Ideal-gas processes

Take $N$ molecules of ideal gas with constant $C_V = \tfrac f2Nk_B$, from $(P_i, V_i, T_i)$ to $(P_f, V_f, T_f)$, quasistatically. In each case $\Delta U = C_V(T_f - T_i)$ regardless of process, $W$ from (3.8), and $Q = \Delta U - W$.

**Isochoric** ($V$ fixed). $W = 0$, $Q = \Delta U = C_V\Delta T$. $P \propto T$.

**Isobaric** ($P$ fixed). $W = -P(V_f - V_i) = -Nk_B\Delta T$. $Q = \Delta U - W = (C_V + Nk_B)\Delta T = C_P\Delta T = \Delta H$.

**Isothermal** ($T$ fixed, gas in contact with a reservoir). $\Delta U = 0$.

$$W = -\int_{V_i}^{V_f}\frac{Nk_BT}{V}\,dV = Nk_BT\ln\frac{V_i}{V_f}, \qquad Q = -W = Nk_BT\ln\frac{V_f}{V_i}. \tag{3.14}$$

Compressing isothermally ($V_f < V_i$) means you do positive work and the same amount of heat leaves into the reservoir. The isotherms are hyperbolas $P = Nk_BT/V$.

**Adiabatic** ($Q = 0$, insulated, quasistatic).

▶ **Try it first.** Set $\delta Q = 0$ in the first law. You have two expressions for $dU$: $C_V\,dT$ and $-P\,dV$. Eliminate $P$ with the gas law and separate variables.

**Derivation.**

1. $Q = 0$ in every step, so $dU = \delta W = -P\,dV$ (quasistatic, (3.8)).
2. For an ideal gas $dU = C_V\,dT = \tfrac f2Nk_B\,dT$ (holds along any path because $U = U(T)$).
3. Equate and use $P = Nk_BT/V$: $\tfrac f2Nk_B\,dT = -\dfrac{Nk_BT}{V}\,dV$.
4. Divide by $Nk_BT$: $\dfrac f2\dfrac{dT}{T} = -\dfrac{dV}{V}$.
5. Integrate: $\tfrac f2\ln T + \ln V = $ const, i.e.
$$VT^{f/2} = \text{const}, \qquad\text{equivalently}\qquad TV^{\gamma-1} = \text{const}, \tag{3.15}$$
since $2/f = \gamma - 1$ by (3.13).
6. Substitute $T = PV/(Nk_B)$: $PV\cdot V^{\gamma-1} = $ const:
$$PV^\gamma = \text{const}. \tag{3.16}$$
7. The work is $W = \Delta U = C_V(T_f - T_i) = \dfrac{P_fV_f - P_iV_i}{\gamma - 1}$. ∎

■ **What just happened.** With no heat flow, the work done on the gas goes entirely into raising $T$, so compressing it raises $P$ twice over: less volume *and* hotter. Adiabats are steeper than isotherms on the $P$–$V$ diagram by a factor $\gamma$ at every point: $(\partial P/\partial V)_{\rm adiabat} = -\gamma P/V$ versus $(\partial P/\partial V)_T = -P/V$. In 2D with $f = 2$, $\gamma = 2$ and $TA = $ const: halve the area, double the temperature. The project checks exactly that.

**Where the adiabatic formula does not apply.** Step 1 used $\delta W = -P\,dV$, which needs quasistatic compression. A free expansion is adiabatic ($Q = 0$) but not quasistatic: $W = 0$, so $\Delta U = 0$ and an ideal gas's $T$ does not change, while (3.15) would predict cooling. A fast piston compression delivers more work than $-\int P\,dV$ along the adiabat because the gas piles up against the piston; the final $T$ is higher than (3.15) predicts. The difference is entropy production (SM-2). The adiabatic formulas need two conditions: no heat flow, and slow.

**The $P$–$V$ picture of all four.**

```
   P
   |\
   | \  adiabat (steeper, P ~ V^-γ)
   |  \__
   |   \ \__   isotherm (P ~ 1/V)
   |    \   \___
   |     \      \_____
   |------\------------ isobar (horizontal)
   |       \
   |  isochore = vertical line
   +-------------------------- V
```

### 3.9 Partial derivatives in thermodynamics

Three variables $P, V, T$, one equation of state, two independent. Every partial derivative must say which variable is held fixed, because $(\partial U/\partial T)_V \ne (\partial U/\partial T)_P$ in general: the first is $C_V$, the second is $C_V + (\partial U/\partial V)_T(\partial V/\partial T)_P$. Most errors in thermodynamics are errors here.

Let $x, y, z$ be related by $F(x, y, z) = 0$, so any one is a function of the other two.

**Theorem 3.3 (Reciprocity and the cyclic rule).**

$$\left(\frac{\partial x}{\partial y}\right)_z = \frac{1}{(\partial y/\partial x)_z}, \tag{3.17}$$

$$\left(\frac{\partial x}{\partial y}\right)_z\left(\frac{\partial y}{\partial z}\right)_x\left(\frac{\partial z}{\partial x}\right)_y = -1. \tag{3.18}$$

▶ **Try it first.** Write $x = x(y, z)$ and $y = y(x, z)$, expand $dx$ and $dy$, substitute one into the other. You get an identity in $dx$ and $dz$ that must hold for arbitrary $dx$, $dz$.

**Derivation.**

1. $dx = \left(\dfrac{\partial x}{\partial y}\right)_z dy + \left(\dfrac{\partial x}{\partial z}\right)_y dz$ and $dy = \left(\dfrac{\partial y}{\partial x}\right)_z dx + \left(\dfrac{\partial y}{\partial z}\right)_x dz$.
2. Substitute the second into the first:
$$dx = \left(\frac{\partial x}{\partial y}\right)_z\left(\frac{\partial y}{\partial x}\right)_z dx + \left[\left(\frac{\partial x}{\partial y}\right)_z\left(\frac{\partial y}{\partial z}\right)_x + \left(\frac{\partial x}{\partial z}\right)_y\right]dz.$$
3. $x$ and $z$ are independent, so the coefficients must match on each side. The $dx$ coefficient: $(\partial x/\partial y)_z(\partial y/\partial x)_z = 1$. That is (3.17).
4. The $dz$ coefficient: $(\partial x/\partial y)_z(\partial y/\partial z)_x + (\partial x/\partial z)_y = 0$, so $(\partial x/\partial y)_z(\partial y/\partial z)_x = -(\partial x/\partial z)_y$.
5. Multiply both sides by $(\partial z/\partial x)_y$ and use (3.17) on the right: $(\partial x/\partial z)_y(\partial z/\partial x)_y = 1$. Result: (3.18). ∎

■ **What just happened.** The $-1$ is the surprise. Naively "cancelling" the $\partial$'s gives $+1$, and that is wrong because each factor holds a different variable fixed. Geometric reason: on a constant-$x$ contour in the $(y,z)$ plane, moving in $y$ forces a compensating move in $z$, and the minus sign records "compensating." Never cancel partials that hold different things fixed.

A third identity you will use constantly (the **change-of-constraint** rule). For a function $u$ of $(x, y)$, and any third variable $w$,

$$\left(\frac{\partial u}{\partial x}\right)_w = \left(\frac{\partial u}{\partial x}\right)_y + \left(\frac{\partial u}{\partial y}\right)_x\left(\frac{\partial y}{\partial x}\right)_w. \tag{3.19}$$

Proof: divide $du = u_x\,dx + u_y\,dy$ by $dx$ along a path of constant $w$. Equation (3.11) is (3.19) with $u = U$, $x = T$, $y = V$, $w = P$, plus the $P(\partial V/\partial T)_P$ term.

**Response functions.** Two partial derivatives are measured often enough to have names:

$$\beta = \frac1V\left(\frac{\partial V}{\partial T}\right)_P\ \ \text{(volume expansion coefficient)}, \qquad \kappa_T = -\frac1V\left(\frac{\partial V}{\partial P}\right)_T\ \ \text{(isothermal compressibility)}. \tag{3.20}$$

For an ideal gas $\beta = 1/T$ and $\kappa_T = 1/P$. The cyclic rule (3.18) with $(x,y,z) = (P, T, V)$ and (3.17) give the pressure rise when a substance is heated in a rigid container:

$$\left(\frac{\partial P}{\partial T}\right)_V = -\frac{(\partial V/\partial T)_P}{(\partial V/\partial P)_T} = \frac{\beta}{\kappa_T}. \tag{3.21}$$

Example 4.6 applies this to a sealed pipe full of water, where the answer is alarming.

### 3.10 Heat transfer in one page

The first law counts heat; it does not say how fast heat moves. Three mechanisms.

**Conduction.** Energy passes by collisions (and in metals, electrons) without bulk motion. **Fourier's law** (Fourier, *Théorie analytique de la chaleur*, 1822): the heat flux $\mathbf q$ (W/m²) is

$$\mathbf q = -k_{\rm th}\nabla T, \tag{3.22}$$

with $k_{\rm th}$ the thermal conductivity (W/(m·K)). It is an empirical linear-response law, like Ohm's.

▶ **Try it first.** Take a slab $[x, x + dx]$ of a rod with cross-section $A$, density $\rho$, specific heat $c$. Write "rate of change of energy inside = flux in − flux out" and take $dx \to 0$.

**Derivation of the heat equation.**

1. Energy in the slab: $\rho c A\,dx\,T(x,t)$ (relative to a reference temperature; heating at fixed volume of a solid, so $c \approx c_P \approx c_V$).
2. Rate of change: $\rho cA\,dx\,\partial_tT$.
3. Net inflow: $A[q(x) - q(x+dx)] = -A\,\partial_xq\,dx$.
4. Equate, cancel $A\,dx$, insert (3.22): $\rho c\,\partial_tT = k_{\rm th}\,\partial_x^2T$, i.e.

$$\frac{\partial T}{\partial t} = \alpha\frac{\partial^2T}{\partial x^2}, \qquad \alpha = \frac{k_{\rm th}}{\rho c}. \tag{3.23}$$

∎ In 3D the same argument with the divergence theorem (EM-1 Theorem 3.8) gives $\partial_tT = \alpha\nabla^2T$. This is the equation F2 §3.12 solves with sine series. In steady state with fixed end temperatures the profile is linear and the flux is constant: a 20 cm brick wall ($k_{\rm th} \approx 0.7$ W/(m·K)) with 20 K across it conducts $q = 0.7\times20/0.2 = 70$ W/m².

**Convection.** Bulk fluid motion carries energy. Engineering summarizes it with Newton's law of cooling $q = h(T_s - T_\infty)$, where $h$ is a heat transfer coefficient that packages the fluid mechanics; computing $h$ is a boundary-layer problem, not a thermodynamics one.

**Radiation.** Every surface emits electromagnetic radiation. A black body emits power per area

$$j = \sigma_{SB}T^4, \qquad \sigma_{SB} = 5.670\times10^{-8}\ \text{W/(m}^2\text{K}^4\text{)}. \tag{3.24}$$

Stefan found the $T^4$ empirically (1879); Boltzmann derived it from thermodynamics plus radiation pressure (1884); the value of $\sigma_{SB}$ in terms of $h$, $c$, and $k_B$ needs Planck, and SM-5 derives it. Real surfaces emit $\varepsilon\sigma_{SB}T^4$ with emissivity $0 \le \varepsilon \le 1$.

### 3.11 Real gases: van der Waals

The ideal gas ignores two facts: molecules take up room, and they attract each other at moderate distances. Johannes van der Waals (thesis, 1873) corrected both.

1. **Excluded volume.** Each mole of molecules occupies a volume $b$ the others cannot enter. Replace $V$ by $V - nb$.
2. **Attraction.** A molecule about to hit the wall is pulled back by the molecules behind it. The pull on each molecule is proportional to the density $n/V$, and the number of molecules hitting the wall is also proportional to $n/V$, so the pressure is reduced by $a(n/V)^2$.

$$\left(P + \frac{an^2}{V^2}\right)(V - nb) = nRT. \tag{3.25}$$

Below a critical temperature the isotherms of (3.25) develop a wiggle where $(\partial P/\partial V)_T > 0$, which is unphysical (squeeze it and the pressure drops) and signals the liquid–gas transition. The critical point is where the wiggle first appears, i.e. where $(\partial P/\partial V)_T = 0$ and $(\partial^2 P/\partial V^2)_T = 0$ together; exercise C1 solves these to get $T_c = 8a/(27Rb)$, $P_c = a/(27b^2)$, $V_c = 3nb$. For CO₂ ($a = 0.364$ Pa·m⁶/mol², $b = 4.267\times10^{-5}$ m³/mol) this gives $T_c = 304$ K, close to the measured 304.1 K; $P_c$ comes out 7.40 MPa versus the measured 7.38 MPa. SM-6 builds the transition (Maxwell construction, critical exponents) on this equation.

## 4. Worked examples

All numbers here are reproduced by `project.py` Part 1 ("Worked examples").

### Example 4.1: Kinetic theory for air at STP

Air at $T = 273.15$ K, $P = 101325$ Pa. Mean molar mass $M = 28.97$ g/mol, so $m = M/N_A = 4.811\times10^{-26}$ kg.

- Number density: $\rho_N = P/(k_BT) = 101325/(1.380649\times10^{-23}\times273.15) = 2.687\times10^{25}$ m⁻³. (Loschmidt's number. A litre holds $2.69\times10^{22}$ molecules, as promised in §2.)
- Mass density: $\rho = \rho_Nm = 1.292$ kg/m³.
- From (3.2): $\langle v^2\rangle = 3P/\rho = 3(101325)/1.2925 = 2.352\times10^5$ m²/s².
- $v_{\rm rms} = \sqrt{\langle v^2\rangle} = 485$ m/s. Equivalently $\sqrt{3k_BT/m}$, same number.
- Mean translational KE: $\tfrac32k_BT = 5.657\times10^{-21}$ J $= 0.0353$ eV.
- From (3.4): $v_p = \sqrt{2k_BT/m} = 396$ m/s, $\bar v = 447$ m/s. The ordering $v_p < \bar v < v_{\rm rms}$ is a property of the right-skewed distribution.
- Molar volume $RT/P = 22.41$ L.

Sanity check: the speed of sound in air at 0 °C is 331 m/s $= \sqrt{\gamma P/\rho} = v_{\rm rms}\sqrt{\gamma/3}$ with $\gamma = 1.4$. Sound is carried by molecules, so it is a bit slower than they are.

### Example 4.2: $C_P - C_V$ for an ideal gas, and the check against air

From (3.12), per mole, $C_P - C_V = R = 8.314$ J/(mol·K), independent of the gas. For air ($f = 5$): $C_V = \tfrac52R = 20.79$ and $C_P = \tfrac72R = 29.10$ J/(mol·K), so $\gamma = 1.400$. Tabulated air at 300 K has $c_p = 1.005$ kJ/(kg·K); times $M = 0.02897$ kg/mol gives 29.11 J/(mol·K). For helium ($f = 3$): $C_V = 12.47$, $C_P = 20.79$, $\gamma = 1.667$.

Directly from the definitions, to see (3.12) as bookkeeping. Heat 1 mol from 300 K to 301 K at constant pressure 1 atm. $\Delta U = C_V\Delta T = 20.79$ J whatever the path. The volume grows by $\Delta V = R\Delta T/P = 8.314/101325 = 8.21\times10^{-5}$ m³ (82 cm³), so the gas does $P\Delta V = 8.314$ J of work on the atmosphere. Heat needed: $20.79 + 8.31 = 29.10$ J $= C_P$.

### Example 4.3: Adiabatic compression of air by a factor of 10

Air at $T_i = 300$ K, $P_i = 1$ atm is compressed quasistatically and adiabatically to $V_f = V_i/10$. $\gamma = 1.4$.

- From (3.15): $T_f = T_i(V_i/V_f)^{\gamma - 1} = 300\times10^{0.4} = 300\times2.512 = 753.6$ K (480 °C).
- From (3.16): $P_f = P_i(V_i/V_f)^\gamma = 10^{1.4} = 25.12$ atm.
- Check with the gas law: $P_fV_f/(P_iV_i) = 25.12/10 = 2.512 = T_f/T_i$. ✓
- Work per mole: $W = C_V(T_f - T_i) = 20.79\times453.6 = 9428$ J, positive (done on the gas). $Q = 0$.
- An isothermal compression by the same ratio would end at 10 atm and 300 K and need $W = RT\ln10 = 5743$ J per mole, all of which would leave as heat.

The adiabatic result is why a Diesel engine (compression ratio 15–20) needs no spark plug: the air is hot enough after compression to ignite the injected fuel. Real compression is neither perfectly adiabatic nor quasistatic, but at engine speeds the heat loss in one stroke is small, so (3.15) is a good first estimate.

### Example 4.4: Path dependence of $Q$ and $W$, same $\Delta U$

One sample of monatomic ideal gas ($U = \tfrac32 PV$ by (3.1) and (3.5)) goes from state 1, $(P_1, V_1) = (100\text{ kPa}, 10\text{ L})$, to state 2, $(200\text{ kPa}, 20\text{ L})$. Then $T_2/T_1 = 4$ and

$$\Delta U = \tfrac32(P_2V_2 - P_1V_1) = \tfrac32(4000 - 1000)\text{ J} = 4500\text{ J}.$$

(kPa × L = J.) Three quasistatic paths:

| Path | Description | $W = -\int P\,dV$ | $Q = \Delta U - W$ |
|---|---|---|---|
| a | expand at 100 kPa to 20 L, then raise $P$ at 20 L | $-100\times10 = -1000$ J | 5500 J |
| b | raise $P$ at 10 L to 200 kPa, then expand at 200 kPa | $-200\times10 = -2000$ J | 6500 J |
| c | straight line on the $P$–$V$ diagram | $-\tfrac12(100+200)\times10 = -1500$ J | 6000 J |

$\Delta U = 4500$ J in all three. Going 1 → 2 by path b and back by path a reversed is a cycle with $\Delta U = 0$, $W_{\rm net} = -2000 + 1000 = -1000$ J (the gas does 1000 J of net work: the rectangle's area) and $Q_{\rm net} = +1000$ J. The cycle turns 1000 J of heat into 1000 J of work. That is an engine, and SM-2 asks how much heat it must also dump.

The project computes all three by integrating $P\,dV$ and $dU$ numerically along parametrized paths, and also checks Theorem 3.2: $\oint dU = 0$ and $\oint\delta Q = -\oint\delta W \neq 0$ around the cycle.

### Example 4.5: An Otto cycle, leg by leg

The air-standard Otto cycle models a spark-ignition engine. One cylinder: $V_1 = 500$ cm³ at bottom dead centre, compression ratio $r = V_1/V_2 = 8$, intake at $P_1 = 100$ kPa, $T_1 = 300$ K. Air, $\gamma = 1.4$, $C_V = \tfrac52R$ per mole. Combustion is modelled as constant-volume heat addition to $T_3 = 1800$ K.

Amount of gas: $n = P_1V_1/(RT_1) = (10^5)(5\times10^{-4})/(8.314\times300) = 0.02005$ mol. $nC_V = 0.4167$ J/K.

**1 → 2, adiabatic compression.** $T_2 = T_1r^{\gamma-1} = 300\times8^{0.4} = 689.2$ K. $P_2 = P_1r^\gamma = 100\times8^{1.4} = 1838$ kPa. $Q_{12} = 0$, $W_{12} = nC_V(T_2 - T_1) = +162.2$ J (on the gas).

**2 → 3, constant-volume heat addition.** $W_{23} = 0$. $Q_{23} = nC_V(T_3 - T_2) = 0.4167\times1110.8 = +462.8$ J. $P_3 = P_2T_3/T_2 = 4800$ kPa.

**3 → 4, adiabatic expansion.** $T_4 = T_3/r^{\gamma-1} = 1800/2.297 = 783.5$ K. $P_4 = P_3/r^\gamma = 261.2$ kPa. $Q_{34} = 0$, $W_{34} = nC_V(T_4 - T_3) = -423.5$ J (the gas does 423.5 J on the piston).

**4 → 1, constant-volume heat rejection** (exhaust modelled as cooling). $W_{41} = 0$. $Q_{41} = nC_V(T_1 - T_4) = -201.5$ J.

**Totals.** $\sum W = 162.2 - 423.5 = -261.4$ J, so the engine delivers 261.4 J per cycle. $\sum Q = 462.8 - 201.5 = 261.4$ J. $\Delta U_{\rm cycle} = 0$ ✓. The ratio of work out to heat in is $261.4/462.8 = 0.565$, which equals $1 - r^{1-\gamma} = 1 - 8^{-0.4} = 0.565$. Why this formula, and why no engine can beat a related bound, is SM-2's job. In Çengel's notation, $W_{\rm net,out} = +261.4$ J and $q_{\rm in} = Q_{23}/m$.

At 3000 rpm a four-stroke engine fires each cylinder 25 times a second, so this one cylinder gives about 6.5 kW in the ideal model. The project integrates $-\oint P\,dV$ numerically around the $P$–$V$ loop and matches $-261.4$ J.

### Example 4.6: The cyclic rule for the ideal gas, and water in a sealed pipe

**Verification for the ideal gas.** $P = Nk_BT/V$, $V = Nk_BT/P$, $T = PV/(Nk_B)$.

$$\left(\frac{\partial P}{\partial T}\right)_V = \frac{Nk_B}{V}, \quad \left(\frac{\partial T}{\partial V}\right)_P = \frac{P}{Nk_B}, \quad \left(\frac{\partial V}{\partial P}\right)_T = -\frac{Nk_BT}{P^2}.$$

Product: $\dfrac{Nk_B}{V}\cdot\dfrac{P}{Nk_B}\cdot\left(-\dfrac{Nk_BT}{P^2}\right) = -\dfrac{Nk_BT}{PV} = -1$ by the gas law. ✓ The project checks this symbolically with SymPy, and also for the van der Waals gas (3.25), where the algebra is not something you want to do by hand.

**A sealed pipe of water.** Liquid water at 20 °C has $\beta = 2.07\times10^{-4}$ K⁻¹ and $\kappa_T = 4.59\times10^{-10}$ Pa⁻¹. Heat a completely full, rigid, sealed pipe by 1 K. By (3.21), $(\partial P/\partial T)_V = \beta/\kappa_T = 4.51\times10^5$ Pa/K $= 4.45$ atm per kelvin. Ten degrees of warming in a closed hydronic loop with no expansion tank means 45 atm. This is why such systems have expansion tanks, and it came from a two-line partial-derivative identity.

## 5. Why it matters and where it connects

- **SM-2.** The first law plus "$\delta Q/T$ is exact" (seen for the ideal gas in §3.4) is the second law. The Otto cycle of Example 4.5 gets its efficiency proof, the Carnot bound, and the entropy accounting of the fast-piston heating. (3.11) becomes $C_P - C_V = TV\beta^2/\kappa_T$ via a Maxwell relation, which is Theorem 3.2's mixed-partial test applied to the potentials.
- **SM-3.** Temperature gets a microscopic definition, $1/T = \partial S/\partial U$, and the Einstein solid explains the heat-capacity freeze-out of §3.3.
- **SM-4.** Theorem 3.1 (equipartition) and the Maxwell distribution (3.4) are derived from the Boltzmann factor. The kinetic-theory identity $\langle\tfrac12mv_x^2\rangle = \tfrac12k_BT$ becomes a Gaussian integral.
- **SM-5.** Stefan–Boltzmann (3.24) and the Debye $T^3$ law.
- **SM-6.** The van der Waals critical point (3.25) and its failure near $T_c$.
- **F2 §3.12.** The heat equation (3.23) solved by Fourier series; the project's rod is that problem with Fourier's flux added.
- **CM-1, EM-1.** Theorem 3.2 is the same mathematics as "a force is conservative iff it is curl-free" (EM-1 Theorem 3.16). $U$ is to $\delta Q + \delta W$ what potential energy is to a conservative force.
- **ST and ML.** The Maxwell distribution is a Gaussian in each velocity component; SM-7 shows softmax is the Boltzmann distribution of SM-4, whose simplest case is (3.4).

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills (rms speeds, sign translation, isothermal work, exactness tests, counting degrees of freedom, response functions). B1–B8 are the module: kinetic theory in $d$ dimensions, adiabatic work, the cyclic rule for van der Waals, the Diesel cycle, a rectangular cycle, the free-expansion trap, the integrating factor $1/T$, and a project extension. B3, B5, B6, and B7 are essential. C1 (van der Waals critical point), C2 (the dry adiabatic lapse rate), C3 (Rüchardt's measurement of $\gamma$) go one step further.

## 7. Coding project

`project.py` has seven parts. (1) Every number in §4, recomputed and asserted. (2) A 2D hard-disk gas: $N$ elastic disks in a box, advanced in small time steps with pair collisions found by a k-d tree. Starting with every disk moving along $\pm x$ at the same speed, you watch $\langle v_x^2\rangle$ and $\langle v_y^2\rangle$ equalize (equipartition) and the speed histogram relax to the 2D Maxwell distribution (3.4). The pressure measured from wall impulses is compared with $Nk_BT/A$ and with its first finite-size correction. (3) A piston: a slow adiabatic compression gives $TA \approx$ const ($\gamma = 2$ in 2D, (3.15)); a slow compression with thermalizing walls gives $W \approx Nk_BT\ln(A_i/A_f)$ (3.14); a fast piston heats the gas beyond the adiabat (irreversibility, SM-2 preview). (4) The Otto cycle's $P$–$V$ diagram and $-\oint P\,dV$ by quadrature. (5) Path dependence: $Q$ and $W$ along three paths. (6) The heat equation on an aluminium rod by finite differences, with the Fourier-law flux at both ends checked against the energy balance. (7) Numbers used in the solutions.

**Build it yourself first:** the wall-impulse pressure estimator and the moving-piston reflection rule $v_x \to 2u_w - v_x$ (a ball bouncing off a moving wall, in the wall's frame). Write them before reading `step()`.

**Your turn** (at the end of the script): van der Waals isotherms for CO₂ with the critical point found numerically; the Diesel cycle; a 3D hard-sphere gas and $\gamma = 5/3$.

## 8. Where you'll struggle

- **The sign of $W$.** You learned $\Delta U = Q - W$ in your ME course and it is now wired in. You will get a sign wrong at least once this week. The fix: translate at the moment of reading, never mid-calculation, and do one sanity check every time: compressing a gas means work is done *on* it, so in this curriculum $W > 0$ for compression. If your formula gives $W < 0$ for a compression, you are using the other convention.
- **"Heat" is a verb.** Students write "the heat in the gas" or "the heat content." There is no such thing (§3.4: $\delta Q$ has no function behind it). Say "energy in the gas" and "heat flowed in." If a sentence would still make sense with "heat" replaced by "energy," and you meant a stock rather than a flow, you meant energy. This is not pedantry: the confusion is exactly the caloric theory, and it took fifty years to kill.
- **Applying $PV^\gamma$ to the wrong process.** (3.16) needs no heat flow *and* quasistatic change. A free expansion satisfies the first and not the second, and students routinely predict it cools an ideal gas. It doesn't (B6). Before using (3.16), say out loud: insulated? slow?
- **Equipartition's failure.** It is tempting to treat "vibrations frozen out" as a detail. It is the crack in classical physics that became quantum mechanics. If you are not bothered that a degree of freedom can contribute zero, reread §3.3 and look at a plot of $C_V(T)$ for H₂.
- **What is held fixed.** This is where most thermodynamic errors live: $(\partial U/\partial T)_P$ is not $C_V$; $\partial$'s holding different things fixed do not cancel; the cyclic rule has a $-1$. The fix is mechanical and non-negotiable: always write the subscript, and when a derivative holds the "wrong" variable fixed, convert it with (3.17)–(3.19). Do B3 by hand even though SymPy can do it. In SM-2 you will juggle eight such derivatives per line.
- **Kinetic theory is a model, thermodynamics isn't.** §3.2 assumes point particles and a specific velocity distribution; the first law assumes nothing about molecules. When the two seem to conflict (real gases, §3.11), the first law wins.

## 9. Check yourself

1. State the zeroth law. Why does it make a thermometer possible?
2. Derive $P = \tfrac13\rho_Nm\langle v^2\rangle$ in five lines. Where did you use isotropy?
3. What is $v_{\rm rms}$ for N₂ at 300 K, to one significant figure?
4. State equipartition. What does it predict for $C_V$ of a diatomic gas, and why is that prediction wrong at room temperature?
5. Write the first law with this curriculum's sign. What is $W$ in Çengel's convention in terms of ours?
6. Show $\delta W = -P\,dV$ is inexact in $(T, V)$ coordinates for an ideal gas.
7. Why is heat added at constant pressure equal to $\Delta H$?
8. Derive $C_P - C_V = Nk_B$ for an ideal gas.
9. Derive $PV^\gamma = $ const. Which two assumptions did you use?
10. Write the cyclic rule. Why isn't the right side $+1$?
11. A gas expands freely into vacuum, insulated. What are $Q$, $W$, $\Delta U$, $\Delta T$ for an ideal gas?
12. Write Fourier's law and derive the 1D heat equation from it.

## 10. Pair with

- **Schroeder, *An Introduction to Thermal Physics*, ch. 1.** The main pairing; same sign convention, same order of ideas. Read for the problems, many of which are better than mine, and §1.6 on heat capacities.
- **Fermi, *Thermodynamics*, ch. 1–2.** Short and classical. The cleanest statement of the first law as a postulate about adiabatic work, without molecules.
- **Feynman Lectures, vol. I, ch. 39–40.** Kinetic theory done the Feynman way, including why the piston speeds up molecules and a first look at the Boltzmann factor.
- **Kittel & Kroemer, *Thermal Physics*, ch. 1 (skim).** Starts from states and counting, which is where SM-3 goes.
- **Çengel & Boles, *Thermodynamics: An Engineering Approach*, ch. 1–4.** Your ME counterpart. Notation map: their $W$ is work by the system ($W_{\rm theirs} = -W_{\rm ours}$); lowercase $u, h, v$ are per unit mass; $R$ is the specific gas constant $R_u/M$; their specific-heat ratio is written $k$ (our $\gamma$; nothing to do with $k_B$); "boundary work" $W_b = \int P\,dV$. Read ch. 4 for real-gas tables and for the control-volume (open system) form of the first law, which these notes do not cover.
