# EM-6: Electrodynamics: Faraday, Maxwell's Equations, Energy and Momentum, Gauge Freedom

**Time:** ~12 hours · **Prerequisites:** EM-5 (and through it EM-1's product rules and EM-4's bound charges) · **Unlocks:** EM-7 (waves), EM-8 (radiation, relativity), QM-1 (photons carry energy and momentum), CM-4 (Noether's theorem for fields)

> **Conventions in this module.** SI units; $\mu_0\varepsilon_0 = 1/c^2$. EMF is $\mathcal E$ (volts), flux is $\Phi = \int_{\mathcal S}\mathbf B\cdot d\mathbf a$ (webers). **Orientation:** a loop and any surface spanning it are related by the right-hand rule (fingers along $d\mathbf l$, thumb along $d\mathbf a$), as in EM-1's Stokes theorem; a positive $\mathcal E$ drives current along $d\mathbf l$. Faraday's law is $\nabla\times\mathbf E = -\partial\mathbf B/\partial t$; Maxwell–Ampère is $\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\,\partial\mathbf E/\partial t$. Energy density $u = \tfrac12(\varepsilon_0E^2 + B^2/\mu_0)$, Poynting vector $\mathbf S = \mathbf E\times\mathbf B/\mu_0$, momentum density $\mathbf g = \varepsilon_0\mathbf E\times\mathbf B = \mathbf S/c^2$. **Stress tensor with Griffiths' sign:** $T_{ij} = \varepsilon_0(E_iE_j - \tfrac12\delta_{ij}E^2) + \mu_0^{-1}(B_iB_j - \tfrac12\delta_{ij}B^2)$, so $T_{ij}\,da_j$ is the force *on* the contents of a closed surface (tension positive) and $-T$ is the momentum flux. Repeated Latin indices are summed. Potentials: $\mathbf E = -\nabla V - \partial_t\mathbf A$, $\mathbf B = \nabla\times\mathbf A$; gauge change $V\to V - \partial_t\lambda$, $\mathbf A\to\mathbf A + \nabla\lambda$; **Lorenz** gauge (Ludvig Lorenz) $\nabla\cdot\mathbf A + c^{-2}\partial_tV = 0$; d'Alembertian $\Box\equiv\nabla^2 - c^{-2}\partial_t^2$ (Griffiths writes $\Box^2$). Source point $\mathbf r'$, separation $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$. Surface jumps are "above minus below" with $\hat{\mathbf n}$ pointing up (EM-2 §3.7).

## 0. Where this sits

EM-2 to EM-5 were statics: $\nabla\cdot\mathbf E = \rho/\varepsilon_0$, $\nabla\times\mathbf E = 0$, $\nabla\cdot\mathbf B = 0$, $\nabla\times\mathbf B = \mu_0\mathbf J$. You need EM-5's vector potential and $\oint\mathbf A\cdot d\mathbf l = \Phi$ (EM-5 (3.25)), EM-4's $\rho_b = -\nabla\cdot\mathbf P$, and EM-1's product rule (P4) and curl-of-curl (3.13). This module adds two time derivatives, which turn four static laws into one dynamical theory with waves, field energy and field momentum. EM-7 solves it for waves; EM-8 solves the potential equations of §3.12 with retarded potentials and packages everything into one tensor.

## 1. The problem that forced this idea

Ørsted (1820) showed that a current makes a magnetic field. Michael Faraday found the converse on 29 August 1831: with two coils on an iron ring, a galvanometer on the second kicked when the current in the first was *switched* on or off, not while it flowed. Within two months he got the same effect with a moving magnet and with a copper disk spinning between magnet poles. The law was about *change*. Joseph Henry found induction independently (published 1832, including self-induction); Heinrich Lenz stated the sign rule in 1834; Franz Neumann wrote the mutual inductance of two circuits as a double line integral in 1845.

Static theory could not survive motion. Ampère's law $\nabla\times\mathbf B = \mu_0\mathbf J$ implies, after taking a divergence, that $\nabla\cdot\mathbf J = 0$: charge can never pile up. But capacitors charge; EM-5 Exercise B1(c) found two surfaces on one loop that disagree about the enclosed current. James Clerk Maxwell, in *On Physical Lines of Force* (1861–62), added a "displacement current" to Ampère's law. He reached it through a mechanical model of the ether, but the requirement he later stated in the *Treatise* (1873) is the one we use: the total current always flows in closed loops. The corrected equations had waves travelling at $1/\sqrt{\mu_0\varepsilon_0}$. Weber and Kohlrausch had measured that combination of electrical constants in 1856 (about $3.1\times10^8$ m/s); Fizeau had measured the speed of light in 1849 (also about $3.1\times10^8$ m/s). Maxwell wrote in 1862: "we can scarcely avoid the inference that light consists in the transverse undulations of the same medium which is the cause of electric and magnetic phenomena."

Hertz generated and detected the waves in 1887–88. Poynting (1884) found that the energy a battery delivers to a resistor flows through the empty space around the wire. Heaviside compressed Maxwell's equations into the four we use. The forcing problem: the static equations contradict each other once anything moves, and the fix unified electricity, magnetism and optics.

## 2. Intuition

**A changing $\mathbf B$ makes $\mathbf E$ curl around it.** Push a magnet into a ring of wire and current flows; move the ring instead and you get the same current. In the first case nothing in the ring moves, so no magnetic force acts on its charges: the push is an electric field, and a circulating one (a conservative field cannot drive charge around a loop). In EM-5 a current made $\mathbf B$ curl around it with circulation $\mu_0I$. Guess the twin: **$\oint\mathbf E\cdot d\mathbf l = -d\Phi/dt$**, with $-\partial\mathbf B/\partial t$ playing the role of $\mu_0\mathbf J$. The minus sign is Lenz: the induced current makes flux that fights the change. With the other sign, a bar sliding on rails would push itself ever faster.

```
   B out of page, growing           induced E circulates clockwise
          . . .                            <----
        . . . . .                        |       ^
          . . .                          v       |
                                           ---->
```

**A changing $\mathbf E$ makes $\mathbf B$ curl around it.** Draw a loop around the wire feeding a charging capacitor. A flat disk on the loop is pierced by the current $I$; a bag-shaped surface on the same loop that passes between the plates is pierced by no current. Same loop, same $\oint\mathbf B\cdot d\mathbf l$. What does the bag catch? A growing electric field: the plate charge grows at rate $I$, so by Gauss the flux of $\varepsilon_0\mathbf E$ through the bag grows at rate $I$. Guess: $\varepsilon_0\,\partial\mathbf E/\partial t$ counts as current. That is Maxwell's term.

**Two curls make a wave.** A changing $\mathbf B$ makes a curling $\mathbf E$, whose change makes a curling $\mathbf B$, and so on, with no charges anywhere. The only speed built from $\mu_0$ (N/A²) and $\varepsilon_0$ (C²/N m²) is $1/\sqrt{\mu_0\varepsilon_0}$. Evaluate it before §3.8 does.

**Energy lives in the field and flows through it.** In a resistor $\mathbf E$ points along the wire and $\mathbf B$ circles it. The one direction built from both that is along neither is $\mathbf E\times\mathbf B$, which points *radially inward*. Guess: energy enters the resistor through its sides. And if field energy moves at $c$, relativity's $E = pc$ suggests a momentum density $u/c$ along the flow, i.e. $\mathbf g = \mathbf S/c^2$. Faraday pictured field lines as stretched strings that push sideways on their neighbours; the stress tensor makes that literal.

## 3. Formal development

### 3.1 Ohm's law, the Drude model, Joule heating

For most conductors the current density is proportional to the force per unit charge:
$$\mathbf J = \sigma\mathbf f\ \xrightarrow{\ \text{conductor at rest}\ }\ \mathbf J = \sigma\mathbf E. \tag{3.1}$$
$\sigma$ is the conductivity (S/m). This is a **constitutive relation**, a statement about a material like EM-4's $\mathbf P = \varepsilon_0\chi_e\mathbf E$, not a law of nature.

**Drude (1900).** Conduction electrons (density $n$, charge $q$, mass $m$) are accelerated by $\mathbf E$ and lose their directed momentum in collisions at mean interval $\tau$: $m\dot{\mathbf v} = q\mathbf E - m\mathbf v/\tau$. In steady state $\mathbf v = q\tau\mathbf E/m$, and $\mathbf J = nq\mathbf v$ gives
$$\sigma = \frac{nq^2\tau}{m}. \tag{3.2}$$
For copper ($n = 8.49\times10^{28}$ m$^{-3}$, $\sigma = 5.96\times10^7$ S/m), $\tau = 2.5\times10^{-14}$ s, and 1 A in a 1 mm² wire drifts at $7.4\times10^{-5}$ m/s (project Part 1). The order of magnitude is right; the temperature dependence needs the Fermi gas of SM-5.

**Joule heating from the force law.** A charge $dq$ with velocity $\mathbf v$ receives power $dq(\mathbf E + \mathbf v\times\mathbf B)\cdot\mathbf v = dq\,\mathbf E\cdot\mathbf v$; the magnetic part does no work (EM-5 Theorem 3.1). With $dq\,\mathbf v = \mathbf J\,d\tau$, for a uniform wire of length $L$, area $A$, $V = EL$, $I = JA$:
$$P = \int\mathbf E\cdot\mathbf J\,d\tau = VI = I^2R,\qquad R = \frac{L}{\sigma A}. \tag{3.3}$$

### 3.2 EMF and motional EMF: the flux rule

A circuit is driven by a source force $\mathbf f_s$ (chemistry, a moving conductor) plus the electrostatic $\mathbf E$ the source sets up. The **electromotive force** is
$$\mathcal E = \oint\mathbf f\cdot d\mathbf l = \oint\mathbf f_s\cdot d\mathbf l, \tag{3.4}$$
since static $\mathbf E$ has zero circulation. It is work per unit charge around the loop, not a force; for one resistive loop $\mathcal E = IR$.

**Motional EMF.** A charge in a moving wire has the wire's velocity $\mathbf v$ plus a drift $\mathbf u$ along the wire. Since $(\mathbf u\times\mathbf B)\cdot d\mathbf l = 0$, only $\mathbf v\times\mathbf B$ contributes to $\mathcal E$.

▶ **Try it first.** For a bar of length $\ell$ sliding at $v$ on rails in a perpendicular uniform $B$ (Example 4.1), compute $\oint(\mathbf v\times\mathbf B)\cdot d\mathbf l$ and $d\Phi/dt$. Then try a general moving loop.

**Theorem 3.1 (Flux rule).** Let a thin wire loop move and deform in a time-independent $\mathbf B$, each element $d\mathbf l$ moving with velocity $\mathbf v$. Then
$$\mathcal E = \oint(\mathbf v\times\mathbf B)\cdot d\mathbf l = -\frac{d\Phi}{dt}. \tag{3.5}$$

*Proof.*
1. In time $dt$ the element $d\mathbf l$ sweeps a parallelogram of vector area $(\mathbf v\,dt)\times d\mathbf l$. (Orientation check: a counterclockwise loop in the $xy$-plane expanding outward has, on its right edge, $\mathbf v\parallel\hat{\mathbf x}$, $d\mathbf l\parallel\hat{\mathbf y}$, and $\hat{\mathbf x}\times\hat{\mathbf y} = +\hat{\mathbf z}$: area added with the loop's own orientation.)
2. The old surface, the swept ribbon and the new surface bound a volume. $\nabla\cdot\mathbf B = 0$ makes the net flux out of it zero, so $\Phi(t + dt) = \Phi(t) + \int_{\rm ribbon}\mathbf B\cdot d\mathbf a$.
3. Hence $d\Phi/dt = \oint\mathbf B\cdot(\mathbf v\times d\mathbf l)$.
4. Cyclically, $\mathbf B\cdot(\mathbf v\times d\mathbf l) = d\mathbf l\cdot(\mathbf B\times\mathbf v) = -(\mathbf v\times\mathbf B)\cdot d\mathbf l$, so $d\Phi/dt = -\mathcal E$. ∎

■ **What just happened.** The flux rule for a moving loop is the magnetic Lorentz force plus $\nabla\cdot\mathbf B = 0$ (which makes the flux surface-independent). The magnetic force still does no work; whoever pulls the bar pays (Example 4.1).

**Where the flux rule needs care.** The proof assumed a thin wire whose material moves with the curve. **Faraday's disk** violates that: a conducting disk of radius $a$ spins at $\omega$ in uniform $B\hat{\mathbf z}$, with brushes at axle and rim. The material moves with $\mathbf v = \omega s\hat{\boldsymbol\phi}$, so $\mathbf v\times\mathbf B = \omega Bs\,\hat{\mathbf s} = \nabla(\tfrac12\omega Bs^2)$ and along any path from axle to rim
$$\mathcal E = \int_0^a\omega Bs\,ds = \tfrac12\omega Ba^2. \tag{3.6}$$
Draw "the circuit" as a fixed radius and its flux never changes, yet the EMF is steady. Feynman's rocking plates (*Lectures* II §17-2) do the reverse: a tiny rocking motion moves the contact point across a large field region, the flux through "the circuit" changes a lot, and the EMF is nearly zero. The rule that never fails is the force law: $\mathcal E = \oint(\mathbf E + \mathbf v\times\mathbf B)\cdot d\mathbf l$ with $\mathbf v$ the velocity of the material the current flows through.

### 3.3 Faraday's law

Faraday found the same $\mathcal E = -d\Phi/dt$ when the loop is at rest and the magnet moves. (Einstein opened his 1905 relativity paper with this coincidence: two different explanations, one measured current.)

▶ **Try it first.** For a stationary loop the charges have $\mathbf v = 0$. What force drives them? Turn $\mathcal E = -d\Phi/dt$ into a differential equation with Stokes.

**Theorem 3.2 (Faraday's law).** If for every fixed loop with spanning surface $\mathcal S$
$$\oint\mathbf E\cdot d\mathbf l = -\frac{d}{dt}\int_{\mathcal S}\mathbf B\cdot d\mathbf a = -\int_{\mathcal S}\frac{\partial\mathbf B}{\partial t}\cdot d\mathbf a, \tag{3.7}$$
then
$$\nabla\times\mathbf E = -\frac{\partial\mathbf B}{\partial t}. \tag{3.8}$$

*Proof.*
1. The loop is at rest, so no magnetic force acts; $\mathcal E = \oint\mathbf E\cdot d\mathbf l$. Experiment sets this equal to $-d\Phi/dt$, and for a fixed surface $d/dt$ passes inside: (3.7).
2. Stokes (EM-1 (3.19)): the left side is $\int_{\mathcal S}(\nabla\times\mathbf E)\cdot d\mathbf a$.
3. So $\int_{\mathcal S}(\nabla\times\mathbf E + \partial_t\mathbf B)\cdot d\mathbf a = 0$ for every surface. Shrink $\mathcal S$ to a small disk with any normal $\hat{\mathbf n}$ at any point: every component of the integrand vanishes. ∎

■ **What just happened.** Electrostatics' $\nabla\times\mathbf E = 0$ was a special case. Once $\mathbf B$ changes, $\int\mathbf E\cdot d\mathbf l$ between two points depends on the path, so a voltmeter's reading depends on how its leads run around the changing flux (Exercise A2(c)).

**The induced field is a magnetostatics problem.** Where $\rho = 0$ the induced field obeys $\nabla\cdot\mathbf E = 0$, $\nabla\times\mathbf E = -\partial_t\mathbf B$: EM-5's equations with $\mu_0\mathbf J\to-\partial_t\mathbf B$. So Ampère-type loops ($\oint\mathbf E\cdot d\mathbf l = -d\Phi/dt$ in place of $\oint\mathbf B\cdot d\mathbf l = \mu_0I_{\rm enc}$) find symmetric induced fields, and Biot–Savart translates to
$$\mathbf E_{\rm induced}(\mathbf r) = -\frac{1}{4\pi}\int\frac{\partial_t\mathbf B(\mathbf r')\times\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\,d\tau'. \tag{3.9}$$
(3.9) treats $\partial_t\mathbf B$ everywhere at the same instant: the **quasistatic** approximation, good when the system is small compared with $c$ times the time scale of change.

**The universal flux rule.** For a loop moving through a time-dependent field, $d\Phi/dt = \int\partial_t\mathbf B\cdot d\mathbf a + \oint\mathbf B\cdot(\mathbf v\times d\mathbf l)$ (field change plus Theorem 3.1's ribbon), so
$$\mathcal E = \oint(\mathbf E + \mathbf v\times\mathbf B)\cdot d\mathbf l = -\frac{d\Phi}{dt}, \tag{3.10}$$
with the exceptions of §3.2.

**A solenoid with a ramping current.** A long solenoid (radius $R$, $n$ turns per length) carries $I(t)$: $\mathbf B = \mu_0nI\hat{\mathbf z}$ inside, $0$ outside (EM-5 (3.21)). By symmetry $\mathbf E = E_\phi(s)\hat{\boldsymbol\phi}$. A circle of radius $s$, counterclockwise about $+\hat{\mathbf z}$, gives $2\pi sE_\phi = -d\Phi/dt$ with $\Phi = \mu_0nI\pi\min(s, R)^2$:
$$E_\phi = -\frac{\mu_0ns}{2}\frac{dI}{dt}\ \ (s < R),\qquad E_\phi = -\frac{\mu_0nR^2}{2s}\frac{dI}{dt}\ \ (s > R). \tag{3.11}$$
Outside, $\mathbf B$ is zero and never changes, yet $\mathbf E\ne0$. Compare EM-5 (3.26): $A_\phi = \mu_0nIs/2$ inside, $\mu_0nIR^2/2s$ outside. So (3.11) is $\mathbf E = -\partial_t\mathbf A$; §3.11 shows why.

### 3.4 Inductance

Current $I_1$ in loop 1 makes $\mathbf B_1\propto I_1$, so its flux through loop 2 is $\Phi_2 = M_{21}I_1$, and a changing $I_1$ induces $\mathcal E_2 = -M_{21}\,dI_1/dt$. $M_{21}$ (henries) is the **mutual inductance**.

▶ **Try it first.** Write $\Phi_2$ as $\oint_2\mathbf A_1\cdot d\mathbf l_2$, insert EM-5's wire formula for $\mathbf A_1$, and look for a symmetry.

**Theorem 3.3 (Neumann formula; reciprocity).**
$$M_{21} = \frac{\mu_0}{4\pi}\oint_1\oint_2\frac{d\mathbf l_1\cdot d\mathbf l_2}{|\mathbf r_1 - \mathbf r_2|} = M_{12}. \tag{3.12}$$

*Proof.*
1. Stokes: $\Phi_2 = \int_{\mathcal S_2}(\nabla\times\mathbf A_1)\cdot d\mathbf a_2 = \oint_2\mathbf A_1\cdot d\mathbf l_2$.
2. For a thin wire (Coulomb gauge), $\mathbf A_1(\mathbf r) = \frac{\mu_0I_1}{4\pi}\oint_1\frac{d\mathbf l_1}{|\mathbf r - \mathbf r_1|}$ (EM-5 (3.24)).
3. Insert and divide by $I_1$: the double integral in (3.12).
4. The integrand is symmetric under $1\leftrightarrow2$, so it is also $M_{12}$. ∎

■ **What just happened.** The flux through loop 2 per ampere in loop 1 equals the flux through loop 1 per ampere in loop 2, however different the loops; in the $\mathbf B$ picture that is far from obvious (Example 4.3(b)). Step 1 holds in any gauge ($\oint\nabla\lambda\cdot d\mathbf l = 0$). $M$ is pure geometry.

**Self-inductance.** A loop's own flux is $\Phi = LI$, and
$$\mathcal E = -L\frac{dI}{dt}. \tag{3.13}$$
(3.12) with loop 1 = loop 2 diverges logarithmically for a filament, so $L$ comes from the flux of a wire of finite size. A long solenoid of length $\ell$ ($N = n\ell$ turns, area $\pi R^2$) links $N\cdot\mu_0nI\pi R^2$. A toroid of $N$ turns with rectangular cross-section ($a < s < b$, height $h$) has $B = \mu_0NI/2\pi s$ (EM-5 (3.22)), so each turn links $\frac{\mu_0NIh}{2\pi}\ln\frac ba$:
$$L_{\rm sol} = \mu_0n^2\pi R^2\ell,\qquad L_{\rm tor} = \frac{\mu_0N^2h}{2\pi}\ln\frac ba. \tag{3.14}$$
For $n = 2000$/m, $R = 2$ cm, $\ell = 0.5$ m: $L_{\rm sol} = 3.158$ mH. For $N = 500$, $h = 1$ cm, $a = 4$ cm, $b = 6$ cm: $L_{\rm tor} = 0.2027$ mH.

### 3.5 Energy in the magnetic field

▶ **Try it first.** To raise the current from $0$ to $I$ the source pushes against the back-EMF (3.13). Compute the work. Then write $LI$ as $\oint\mathbf A\cdot d\mathbf l$.

**Theorem 3.4 (Magnetic energy).** The work to establish currents slowly from zero is
$$W = \tfrac12LI^2\ \ (\text{one loop}), \tag{3.15}$$
$$W = \tfrac12\int\mathbf A\cdot\mathbf J\,d\tau, \tag{3.16}$$
$$W = \frac{1}{2\mu_0}\int_{\text{all space}}B^2\,d\tau. \tag{3.17}$$

*Proof.* **Route 1, the back-EMF.**
1. Against $\mathcal E = -L\,dI/dt$ the source does work $L\,dI/dt$ per unit charge on charge flowing at rate $I$: $dW/dt = LI\,dI/dt$. Integrating gives (3.15).
2. $LI = \Phi = \oint\mathbf A\cdot d\mathbf l$, so $W = \tfrac12I\oint\mathbf A\cdot d\mathbf l$; for a distribution, $I\,d\mathbf l\to\mathbf J\,d\tau$ gives (3.16). (Directly: against the induced $\mathbf E = -\partial_t\mathbf A$ the source supplies $\int\partial_t\mathbf A\cdot\mathbf J\,d\tau$ per unit time; ramp every current by a common factor $\alpha$ from $0$ to $1$ and $\int_0^1\alpha\,d\alpha$ supplies the $\tfrac12$.)

**Route 2, the product rule.**
3. With $\mathbf J = \nabla\times\mathbf B/\mu_0$: $W = \frac{1}{2\mu_0}\int\mathbf A\cdot(\nabla\times\mathbf B)\,d\tau$.
4. (P4): $\mathbf A\cdot(\nabla\times\mathbf B) = \mathbf B\cdot(\nabla\times\mathbf A) - \nabla\cdot(\mathbf A\times\mathbf B) = B^2 - \nabla\cdot(\mathbf A\times\mathbf B)$.
5. Over a ball of radius $r$ containing all currents, $W = \frac{1}{2\mu_0}\left[\int B^2\,d\tau - \oint(\mathbf A\times\mathbf B)\cdot d\mathbf a\right]$. Far away $A\sim r^{-2}$, $B\sim r^{-3}$, so the surface term is $\sim r^{-5}r^2\to0$: (3.17). ∎

■ **What just happened.** The same energy is an integral over the currents and an integral over the field, the twin of EM-2 (3.15)–(3.16). Check on the solenoid: $(B^2/2\mu_0)\pi R^2\ell = \tfrac12\mu_0n^2I^2\pi R^2\ell = \tfrac12L_{\rm sol}I^2$; for the toroid the project integrates $B^2/2\mu_0$ and recovers $\tfrac12L_{\rm tor}I^2$.

**Electric versus magnetic storage.** At $1$ T, $u_B = 3.98\times10^5$ J/m³. Air breaks down near $3$ MV/m, where $u_E = \tfrac12\varepsilon_0E^2 = 39.8$ J/m³. Matching $u_B$ would need $E = cB = 3\times10^8$ V/m. That factor of $10^4$ is why bulk energy is stored in magnets (an MRI magnet holds megajoules) and why capacitors need thin gaps and dielectrics.

### 3.6 Maxwell's correction

▶ **Try it first.** Take the divergence of $\nabla\times\mathbf B = \mu_0\mathbf J$ and compare with continuity (EM-5 Theorem 3.2). Then use Gauss's law to find a vector whose divergence is $-\nabla\cdot\mathbf J$.

**Theorem 3.5 (Displacement current).** Ampère's law contradicts charge conservation whenever $\partial_t\rho\ne0$. The modification
$$\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\frac{\partial\mathbf E}{\partial t} \tag{3.18}$$
is consistent with both continuity and Gauss's law.

*Derivation.*
1. $\nabla\cdot(\nabla\times\mathbf B) = 0$ (EM-1 (3.12)), so Ampère's law forces $\nabla\cdot\mathbf J = 0$.
2. Continuity: $\nabla\cdot\mathbf J = -\partial_t\rho$. Steps 1–2 agree only if $\partial_t\rho = 0$. Ampère's law is magnetostatics and nothing more.
3. Differentiate Gauss's law: $\partial_t\rho = \nabla\cdot(\varepsilon_0\partial_t\mathbf E)$.
4. So $\nabla\cdot(\mathbf J + \varepsilon_0\partial_t\mathbf E) = 0$ always. Put this divergence-free combination where $\mathbf J$ was: (3.18). Its divergence now reproduces continuity. ∎

■ **What just happened.** $\mathbf J_d\equiv\varepsilon_0\partial_t\mathbf E$ is the **displacement current density**. Steps 1–4 show it is the minimal fix consistent with Gauss's law, not the only conceivable one; experiment decides, and the waves of §3.8 exist. In integral form, with $\Phi_E = \int\mathbf E\cdot d\mathbf a$,
$$\oint\mathbf B\cdot d\mathbf l = \mu_0I_{\rm enc} + \mu_0\varepsilon_0\frac{d\Phi_E}{dt}, \tag{3.19}$$
and now any two surfaces on one loop agree. Example 4.4 has the numbers. The term is negligible in low-frequency circuits because $\varepsilon_0$ is tiny (Exercise A4) and is the only source of $\mathbf B$ in empty space.

### 3.7 Maxwell's equations, in vacuum and in matter

$$\boxed{\ \nabla\cdot\mathbf E = \frac{\rho}{\varepsilon_0},\qquad\nabla\cdot\mathbf B = 0,\qquad\nabla\times\mathbf E = -\frac{\partial\mathbf B}{\partial t},\qquad\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\frac{\partial\mathbf E}{\partial t}\ } \tag{3.20}$$
With the Lorentz force this is all of classical electrodynamics; $\rho$ and $\mathbf J$ are *all* charges and currents. Helmholtz's theorem (EM-1 Theorem 3.15) says divergence and curl fix a field that vanishes at infinity, and (3.20) supplies both.

**In matter.** EM-4 wrote $\rho = \rho_f - \nabla\cdot\mathbf P$; EM-5 wrote $\mathbf J = \mathbf J_f + \nabla\times\mathbf M$ for steady currents. A *changing* polarization adds a third current.

▶ **Try it first.** If $\mathbf P$ in a slab grows, the surface charge $\mathbf P\cdot\hat{\mathbf n}$ grows. What current density delivers it? Check continuity with $\rho_b$.

*Derivation.*
1. Increasing $\mathbf P$ by $d\mathbf P$ moves bound charge $d\mathbf P\cdot d\mathbf a$ across an area element (EM-4 §3.5's count at a tilted surface). So the **polarization current** is
$$\mathbf J_p = \frac{\partial\mathbf P}{\partial t}. \tag{3.21}$$
2. Check: $\nabla\cdot\mathbf J_p = \partial_t(\nabla\cdot\mathbf P) = -\partial_t\rho_b$, continuity for bound charge. ($\nabla\times\mathbf M$ is divergence-free and needs no partner.)
3. Insert $\rho = \rho_f - \nabla\cdot\mathbf P$ and $\mathbf J = \mathbf J_f + \nabla\times\mathbf M + \partial_t\mathbf P$ into (3.20). Gauss: $\nabla\cdot(\varepsilon_0\mathbf E + \mathbf P) = \rho_f$. Maxwell–Ampère: $\nabla\times(\mathbf B/\mu_0 - \mathbf M) = \mathbf J_f + \partial_t(\varepsilon_0\mathbf E + \mathbf P)$. With $\mathbf D = \varepsilon_0\mathbf E + \mathbf P$ and $\mathbf H = \mathbf B/\mu_0 - \mathbf M$:
$$\nabla\cdot\mathbf D = \rho_f,\qquad\nabla\cdot\mathbf B = 0,\qquad\nabla\times\mathbf E = -\frac{\partial\mathbf B}{\partial t},\qquad\nabla\times\mathbf H = \mathbf J_f + \frac{\partial\mathbf D}{\partial t}. \tag{3.22}$$
∎ These need constitutive relations ($\mathbf D = \varepsilon\mathbf E$, $\mathbf H = \mathbf B/\mu$, $\mathbf J_f = \sigma\mathbf E$) before they can be solved.

**Boundary conditions.** At a surface with free charge $\sigma_f$ and free current $\mathbf K_f$:
$$D^\perp_{\rm above} - D^\perp_{\rm below} = \sigma_f,\quad B^\perp_{\rm above} = B^\perp_{\rm below},\quad\mathbf E^\parallel_{\rm above} = \mathbf E^\parallel_{\rm below},\quad\mathbf H^\parallel_{\rm above} - \mathbf H^\parallel_{\rm below} = \mathbf K_f\times\hat{\mathbf n}. \tag{3.23}$$
The divergence equations are unchanged, so the first two are EM-4 (3.16) and EM-5 (3.28). The curl equations gained time derivatives; we check that the last two survive.

1. Take a rectangle straddling the surface: sides of length $\ell$ along a unit vector $\hat{\boldsymbol\ell}$ in the surface, just above and just below; height $\epsilon$. Apply (3.7).
2. As $\epsilon\to0$ the short sides contribute nothing (finite $\mathbf E$) and the flux of $\mathbf B$ through the rectangle is at most $\ell\epsilon\max|\partial_t\mathbf B|\to0$.
3. What remains is $(\mathbf E_{\rm above} - \mathbf E_{\rm below})\cdot\hat{\boldsymbol\ell} = 0$ for every $\hat{\boldsymbol\ell}$: $\mathbf E^\parallel$ is continuous.
4. Repeat with $\oint\mathbf H\cdot d\mathbf l = I_{f,\rm enc} + \frac{d}{dt}\int\mathbf D\cdot d\mathbf a$. The $\mathbf D$ flux vanishes as in step 2; the free current through the rectangle is $\ell\,\mathbf K_f\cdot(\hat{\mathbf n}\times\hat{\boldsymbol\ell})$, exactly the static case EM-5 (3.42). ∎

■ **What just happened.** A flux of a finite field through a vanishing area vanishes, so only *surface* sources cause jumps. EM-7 uses (3.23) at every interface.

### 3.8 The wave equation and the speed of light

▶ **Try it first.** In vacuum take the curl of Faraday's law. Use EM-1 (3.13) on the left and Maxwell–Ampère on the right.

**Theorem 3.6 (Waves).** In vacuum ($\rho = 0$, $\mathbf J = 0$)
$$\nabla^2\mathbf E = \mu_0\varepsilon_0\frac{\partial^2\mathbf E}{\partial t^2},\qquad\nabla^2\mathbf B = \mu_0\varepsilon_0\frac{\partial^2\mathbf B}{\partial t^2}. \tag{3.24}$$

*Derivation.*
1. $\nabla\times(\nabla\times\mathbf E) = -\partial_t(\nabla\times\mathbf B)$.
2. Left: $\nabla(\nabla\cdot\mathbf E) - \nabla^2\mathbf E = -\nabla^2\mathbf E$, since $\nabla\cdot\mathbf E = 0$.
3. Right: $-\mu_0\varepsilon_0\partial_t^2\mathbf E$. Equate.
4. The curl of Maxwell–Ampère gives $-\nabla^2\mathbf B = \mu_0\varepsilon_0\partial_t(\nabla\times\mathbf E) = -\mu_0\varepsilon_0\partial_t^2\mathbf B$. ∎

■ **What just happened.** Every component obeys the wave equation with speed
$$c = \frac{1}{\sqrt{\mu_0\varepsilon_0}} = 2.998\times10^8\ \text{m/s}. \tag{3.25}$$
Constants measured with charged spheres and current balances give the speed of light. Without the displacement current step 3 has no $\partial_t^2$ and no wave. (Today $c$ is fixed by definition and $\varepsilon_0$ is defined as $1/\mu_0c^2$, so (3.25) holds by construction; in 1862 it was a discovery.)

The coupled equations say more. Try $\mathbf E = f(x - ct)\hat{\mathbf y}$: Faraday gives $\partial_tB_z = -\partial_xE_y = -f'$, so $B_z = f(x - ct)/c$. $\mathbf B$ is perpendicular to $\mathbf E$ and to the motion, $|\mathbf B| = |\mathbf E|/c$, and $\mathbf E\times\mathbf B$ points along the motion. The project launches its FDTD pulse from this pair.

**Magnetic charge.** Add $\rho_m$ and $\mathbf J_m$ via $\nabla\cdot\mathbf B = \mu_0\rho_m$ and $\nabla\times\mathbf E = -\mu_0\mathbf J_m - \partial_t\mathbf B$, and the equations become symmetric; in vacuum they are already invariant under the **duality rotation** $\mathbf E\to c\mathbf B$, $c\mathbf B\to-\mathbf E$. No monopole has been found. Dirac (1931) showed that if one exists, quantum mechanics requires every charge $q$ to satisfy $q\Phi_m = nh$ ($n$ an integer, $\Phi_m$ the monopole's flux): EM-5's Aharonov–Bohm phase $q\Phi_m/\hbar$ around the invisible tube bringing the flux in from infinity must be a multiple of $2\pi$. One monopole anywhere would explain why all charges are multiples of $e$.

### 3.9 Poynting's theorem

EM-5 Theorem 3.2 is the template for every conservation law in field theory: $\partial_t\rho + \nabla\cdot\mathbf J = 0$ says the amount in a fixed volume changes only by flow through its boundary. Energy and momentum of the fields follow the same pattern, with a source term for exchange with matter.

▶ **Try it first.** Start from the power delivered to charges, $\int\mathbf E\cdot\mathbf J\,d\tau$ (3.3). Eliminate $\mathbf J$ with Maxwell–Ampère; you will need (P4) and Faraday. Do it on paper; everyone loses a sign here once.

**Theorem 3.7 (Poynting).** With
$$u = \frac12\left(\varepsilon_0E^2 + \frac{B^2}{\mu_0}\right),\qquad\mathbf S = \frac{1}{\mu_0}\mathbf E\times\mathbf B, \tag{3.26}$$
the work done per unit time on the charges in a fixed volume $\mathcal V$ with boundary $\mathcal S$ is
$$\frac{dW_{\rm mech}}{dt} = \int_{\mathcal V}\mathbf E\cdot\mathbf J\,d\tau = -\frac{d}{dt}\int_{\mathcal V}u\,d\tau - \oint_{\mathcal S}\mathbf S\cdot d\mathbf a, \tag{3.27}$$
and pointwise
$$\frac{\partial u}{\partial t} + \nabla\cdot\mathbf S = -\mathbf E\cdot\mathbf J. \tag{3.28}$$

*Derivation.*
1. The magnetic force does no work, so $dW_{\rm mech}/dt = \int\mathbf E\cdot\mathbf J\,d\tau$.
2. Maxwell–Ampère: $\mathbf E\cdot\mathbf J = \mu_0^{-1}\mathbf E\cdot(\nabla\times\mathbf B) - \varepsilon_0\mathbf E\cdot\partial_t\mathbf E$.
3. (P4): $\mathbf E\cdot(\nabla\times\mathbf B) = \mathbf B\cdot(\nabla\times\mathbf E) - \nabla\cdot(\mathbf E\times\mathbf B)$.
4. Faraday: $\mathbf B\cdot(\nabla\times\mathbf E) = -\mathbf B\cdot\partial_t\mathbf B = -\tfrac12\partial_t(B^2)$; and $\mathbf E\cdot\partial_t\mathbf E = \tfrac12\partial_t(E^2)$.
5. Assemble: $\mathbf E\cdot\mathbf J = -\frac{1}{2\mu_0}\partial_tB^2 - \frac{\varepsilon_0}{2}\partial_tE^2 - \frac{1}{\mu_0}\nabla\cdot(\mathbf E\times\mathbf B) = -\partial_tu - \nabla\cdot\mathbf S$: (3.28).
6. Integrate over $\mathcal V$ and use the divergence theorem: (3.27). ∎

■ **What just happened.** Energy given to matter equals the drop in field energy inside plus the energy flowing *in* through the surface ($-\oint\mathbf S\cdot d\mathbf a$). $\mathbf S$ is the energy flux (W/m²), and (3.28) is the continuity template with source $-\mathbf E\cdot\mathbf J$. The theorem fixes only $\oint\mathbf S\cdot d\mathbf a$; reading $\mathbf S$ as the local flow is the standard interpretation, and it is the one radiation confirms. For a plane wave $|\mathbf S| = E^2/\mu_0c = cu$: the energy moves at $c$. The project checks (3.28) symbolically for fields from random polynomial potentials.

### 3.10 Momentum in the fields: the Maxwell stress tensor

Magnetic forces between moving charges violate Newton's third law (EM-5's current elements pull unequally). If momentum is conserved, the field must hold some. We find how much by eliminating the sources from $\mathbf f = \rho\mathbf E + \mathbf J\times\mathbf B$, as in §3.9. It is long; every step is a product rule.

▶ **Try it first.** Insert $\rho = \varepsilon_0\nabla\cdot\mathbf E$ and Maxwell–Ampère's $\mathbf J$ into $\mathbf f$ and try to write the $\mathbf E$ terms as $\partial_j(\text{something})_{ij}$.

**Theorem 3.8 (Momentum balance).** With
$$T_{ij} = \varepsilon_0\left(E_iE_j - \tfrac12\delta_{ij}E^2\right) + \frac{1}{\mu_0}\left(B_iB_j - \tfrac12\delta_{ij}B^2\right), \tag{3.29}$$
the force density on charges is $f_i = \partial_jT_{ij} - \varepsilon_0\mu_0\partial_tS_i$, and the force on the charges in $\mathcal V$ is
$$\mathbf F = \oint_{\mathcal S}\overleftrightarrow{T}\cdot d\mathbf a - \varepsilon_0\mu_0\frac{d}{dt}\int_{\mathcal V}\mathbf S\,d\tau. \tag{3.30}$$

*Derivation.*
1. $\mathbf f = \varepsilon_0(\nabla\cdot\mathbf E)\mathbf E + \mu_0^{-1}(\nabla\times\mathbf B)\times\mathbf B - \varepsilon_0(\partial_t\mathbf E)\times\mathbf B$.
2. $\partial_t(\mathbf E\times\mathbf B) = (\partial_t\mathbf E)\times\mathbf B + \mathbf E\times\partial_t\mathbf B = (\partial_t\mathbf E)\times\mathbf B - \mathbf E\times(\nabla\times\mathbf E)$ by Faraday. So the last term of step 1 is $-\varepsilon_0\partial_t(\mathbf E\times\mathbf B) - \varepsilon_0\mathbf E\times(\nabla\times\mathbf E)$.
3. Add $\mu_0^{-1}(\nabla\cdot\mathbf B)\mathbf B = 0$ for symmetry:
$$\mathbf f = \varepsilon_0\big[(\nabla\cdot\mathbf E)\mathbf E - \mathbf E\times(\nabla\times\mathbf E)\big] + \frac{1}{\mu_0}\big[(\nabla\cdot\mathbf B)\mathbf B - \mathbf B\times(\nabla\times\mathbf B)\big] - \varepsilon_0\partial_t(\mathbf E\times\mathbf B).$$
4. EM-1 (3.1): $[\mathbf E\times(\nabla\times\mathbf E)]_i = \varepsilon_{ijk}\varepsilon_{klm}E_j\partial_lE_m = E_j\partial_iE_j - E_j\partial_jE_i = \tfrac12\partial_i(E^2) - E_j\partial_jE_i$.
5. So the $\mathbf E$ bracket is $E_i\partial_jE_j + E_j\partial_jE_i - \tfrac12\partial_iE^2 = \partial_j(E_iE_j - \tfrac12\delta_{ij}E^2)$: the product rule read backwards. Same for $\mathbf B$.
6. Hence $f_i = \partial_jT_{ij} - \varepsilon_0\mu_0\partial_tS_i$; integrate and apply the divergence theorem row by row: (3.30). ∎

■ **What just happened.** $T_{ij}\,da_j$ is the $i$-th component of the force on a surface element whose normal is along $j$. That is Cauchy's traction $\mathbf t = \sigma\mathbf n$ of LA-12 §3.18: $T$ is a symmetric $(0,2)$ tensor, matrix times area vector gives force. For $\mathbf E = E\hat{\mathbf x}$ alone, $T = \tfrac12\varepsilon_0E^2\,\mathrm{diag}(1, -1, -1)$: a **tension** $\tfrac12\varepsilon_0E^2$ along the lines and a **pressure** of the same size across them. That is Faraday's picture of stretched, mutually repelling lines, made quantitative; at a conductor it is EM-2's surface pressure (3.17). At $1$ T the magnetic stress is $3.98\times10^5$ Pa, about 3.9 atmospheres on a magnet's windings.

**Field momentum.** With $\mathbf F = d\mathbf p_{\rm mech}/dt$, (3.30) reads
$$\frac{d\mathbf p_{\rm mech}}{dt} = -\frac{d}{dt}\int_{\mathcal V}\mathbf g\,d\tau + \oint_{\mathcal S}\overleftrightarrow{T}\cdot d\mathbf a,\qquad\mathbf g = \varepsilon_0\mathbf E\times\mathbf B = \frac{\mathbf S}{c^2},\qquad\frac{\partial\mathbf g}{\partial t} + \nabla\cdot(-\overleftrightarrow{T}) = -\mathbf f. \tag{3.31}$$
$\mathbf g$ is the field's momentum density and $-T$ its momentum flux: the continuity template once per component. For light, $g/u = 1/c$, the classical side of the photon's $p = E/c$ (QM-1).

**Angular momentum** density is $\mathbf r\times\mathbf g$, and static fields can hold it. Feynman's **disk paradox** (*Lectures* II §17-4, resolved in §27-6): a disk on a frictionless axle carries a coaxial solenoid with steady current and charged spheres near its rim. Switch the current off; the induced $\mathbf E$ (3.11) torques the charges and the disk turns, with nothing external acting. The angular momentum was in the static $\varepsilon_0\mathbf r\times(\mathbf E\times\mathbf B)$ all along.

### 3.11 The potentials in electrodynamics

$\nabla\cdot\mathbf B = 0$ still gives $\mathbf B = \nabla\times\mathbf A$, but $\nabla\times\mathbf E\ne0$ kills $\mathbf E = -\nabla V$.

▶ **Try it first.** Insert $\mathbf B = \nabla\times\mathbf A$ into Faraday's law and find a combination of $\mathbf E$ and $\mathbf A$ with zero curl.

*Derivation.*
1. $\nabla\times(\mathbf E + \partial_t\mathbf A) = 0$, so (EM-1 Theorem 3.16) $\mathbf E + \partial_t\mathbf A = -\nabla V$:
$$\mathbf E = -\nabla V - \frac{\partial\mathbf A}{\partial t},\qquad\mathbf B = \nabla\times\mathbf A. \tag{3.32}$$
2. The homogeneous equations now hold identically. Gauss gives $-\nabla^2V - \partial_t(\nabla\cdot\mathbf A) = \rho/\varepsilon_0$.
3. Maxwell–Ampère gives $\nabla\times(\nabla\times\mathbf A) = \mu_0\mathbf J - \mu_0\varepsilon_0\nabla\partial_tV - \mu_0\varepsilon_0\partial_t^2\mathbf A$; use EM-1 (3.13) and collect:
$$\nabla^2V + \frac{\partial}{\partial t}(\nabla\cdot\mathbf A) = -\frac{\rho}{\varepsilon_0},\qquad\Box\mathbf A - \nabla\left(\nabla\cdot\mathbf A + \frac{1}{c^2}\frac{\partial V}{\partial t}\right) = -\mu_0\mathbf J. \tag{3.33}$$
∎

■ **What just happened.** Six field components obeying eight equations became four potentials obeying four coupled equations; the freedom in the potentials will uncouple them. The solenoid of §3.3 is (3.32) with $V = 0$.

### 3.12 Gauge transformations; the Coulomb and Lorenz gauges

**Theorem 3.9 (Gauge invariance).** For any smooth $\lambda(\mathbf r, t)$,
$$\mathbf A' = \mathbf A + \nabla\lambda,\qquad V' = V - \frac{\partial\lambda}{\partial t} \tag{3.34}$$
give the same $\mathbf E$ and $\mathbf B$.

*Proof.* $\nabla\times\mathbf A' = \mathbf B + \nabla\times\nabla\lambda = \mathbf B$ (EM-1 (3.11)). $-\nabla V' - \partial_t\mathbf A' = \mathbf E + \nabla\partial_t\lambda - \partial_t\nabla\lambda = \mathbf E$, since mixed partials commute. ∎

Conversely, if two pairs give the same fields, $\nabla\times(\mathbf A' - \mathbf A) = 0$ makes $\mathbf A' - \mathbf A = \nabla\lambda$, and then $\nabla(V' - V + \partial_t\lambda) = 0$ leaves a function of $t$ alone, absorbable into $\lambda$: (3.34) is the whole freedom. The project runs the proof on a grid: central differences along different axes commute exactly, so a gauge change that moves $V$ and $\mathbf A$ by $O(1)$ moves the computed $\mathbf E$ and $\mathbf B$ by $2\times10^{-14}$.

**Coulomb gauge, $\nabla\cdot\mathbf A = 0$** (reachable as in EM-5 §3.7). The first of (3.33) becomes
$$\nabla^2V = -\frac{\rho}{\varepsilon_0},\qquad V(\mathbf r, t) = \frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\mathbf r', t)}{\mathscr r}\,d\tau'. \tag{3.35}$$
$V$ responds *instantly* to charge everywhere. That breaks nothing, because $V$ is not measurable: the instantaneous part of $\mathbf A$'s source ($\varepsilon_0\nabla\partial_tV$ in (3.33)) cancels the instantaneous part of $-\nabla V$ in $\mathbf E$. (Skipped: the proof that $\mathbf E$ comes out retarded; Jackson, *Classical Electrodynamics*, §6.3.) The Coulomb gauge suits statics and the quantized field.

▶ **Try it first.** Which condition kills the gradient term in the second of (3.33)? Insert it into the first.

**Lorenz gauge.**
$$\nabla\cdot\mathbf A + \frac{1}{c^2}\frac{\partial V}{\partial t} = 0. \tag{3.36}$$

*Derivation.*
1. *Reachable:* under (3.34), $\nabla\cdot\mathbf A + c^{-2}\partial_tV$ changes by $\nabla^2\lambda - c^{-2}\partial_t^2\lambda = \Box\lambda$. Choose $\lambda$ solving $\Box\lambda = -(\nabla\cdot\mathbf A_0 + c^{-2}\partial_tV_0)$, an inhomogeneous wave equation, which always has solutions (EM-8 writes them).
2. The second of (3.33) loses its gradient: $\Box\mathbf A = -\mu_0\mathbf J$.
3. In the first, $\partial_t(\nabla\cdot\mathbf A) = -c^{-2}\partial_t^2V$:
$$\Box V = -\frac{\rho}{\varepsilon_0},\qquad\Box\mathbf A = -\mu_0\mathbf J. \tag{3.37}$$
∎

■ **What just happened.** Four uncoupled copies of one wave equation, each potential driven by its own source. In statics they reduce to the Poisson equations of EM-2 and EM-5. EM-8 solves them with retarded potentials, $V(\mathbf r, t) = \frac{1}{4\pi\varepsilon_0}\int\rho(\mathbf r', t - \mathscr r/c)/\mathscr r\,d\tau'$. The project checks with SymPy that the exact potentials of an oscillating dipole satisfy (3.36), solve (3.37) off the origin, and give fields obeying all four vacuum equations. (3.36) leaves a residual freedom ($\Box\lambda = 0$) and is relativistically invariant ($\partial_\mu A^\mu = 0$ in EM-8). It is Ludvig Lorenz's (1867), not Hendrik Lorentz's.

### 3.13 Gauge invariance, the Lagrangian, and Noether

Gauge freedom is the deepest structural fact in the subject. In quantum mechanics the Schrödinger equation for a charge contains $\mathbf p - q\mathbf A$ and $qV$, and under (3.34) it keeps its form if $\psi\to e^{iq\lambda/\hbar}\psi$, a *local* phase change (QM-4). Read backwards: demanding that the phase of $\psi$ be freely choosable at each point forces a field that transforms like $(V, \mathbf A)$. In superconducting circuits the same phase quantizes the flux through a ring in units of $h/2e$, the basis of SQUIDs and flux-tunable qubits (the hardware behind QC-1's two-level systems).

Classically, the charged particle's Lagrangian (CM-3 (4.17)) changes under (3.34) by $q\,d\lambda/dt$, a total derivative that leaves the motion unchanged (CM-3 §5). The fields have a Lagrangian density too:
$$\mathcal L = \frac{\varepsilon_0}{2}E^2 - \frac{1}{2\mu_0}B^2 - \rho V + \mathbf J\cdot\mathbf A\quad\Big(= -\frac{1}{4\mu_0}F_{\mu\nu}F^{\mu\nu} - A_\mu J^\mu\ \text{in EM-8's notation}\Big), \tag{3.38}$$
with $\mathbf E$, $\mathbf B$ written through (3.32). Its Euler–Lagrange equations are Maxwell's. The Gauss half: $\mathcal L$ contains $V$ through $-\rho V$ and $\tfrac{\varepsilon_0}{2}|\nabla V + \partial_t\mathbf A|^2$, so $\partial\mathcal L/\partial V - \nabla\cdot\big(\partial\mathcal L/\partial(\nabla V)\big) = -\rho + \varepsilon_0\nabla\cdot\mathbf E = 0$. Exercise C3 varies $\mathbf A$ and gets Maxwell–Ampère.

Now apply CM-4's Noether theorem in its field form (CM-4 §3.10). Time-translation invariance of $\mathcal L$ gives a conserved energy current: Poynting's theorem (3.28). Space translations give momentum with flux $-T_{ij}$: (3.31). Rotations give field angular momentum: Feynman's disk. And the source term $-\rho V + \mathbf J\cdot\mathbf A$ changes under (3.34) by a total derivative exactly when $\partial_t\rho + \nabla\cdot\mathbf J = 0$: gauge invariance and charge conservation are tied, and charge conservation is what forced Maxwell's correction in §3.6. The four big results of this module are one theorem used four times.

## 4. Worked examples

### Example 4.1: The bar on rails, and the power balance (project Part 1)

*From mechanical engineering: a linear generator, and run backwards an eddy-current brake.* Rails $\ell = 0.5$ m apart are joined at one end by $R = 2\ \Omega$; a bar slides at $v = 3$ m/s in $B = 0.8$ T perpendicular to the plane. Rails and bar have negligible resistance.

The flux $\Phi = B\ell x$ grows at $B\ell v$, so $|\mathcal E| = B\ell v = 1.200$ V and $I = \mathcal E/R = 0.600$ A, circulating so its flux opposes the growth (Lenz). The bar feels $F = I\ell B = 0.240$ N opposing $v$. Keeping $v$ constant costs $Fv = 0.720$ W; the resistor dissipates $I^2R = 0.720$ W. In general $Fv = (B\ell v)^2/R = \mathcal E^2/R$.

*Coasting.* Release the bar (mass $m = 0.1$ kg): $m\dot v = -(B\ell)^2v/R$, so $v = v_0e^{-t/\tau}$ with $\tau = mR/(B\ell)^2 = 1.250$ s. Total heat $\int_0^\infty I^2R\,dt = \frac{(B\ell v_0)^2}{R}\frac\tau2 = \tfrac12mv_0^2 = 0.4500$ J: all the kinetic energy, as the project's ODE integration confirms. Eddy-current brakes on trains work this way, with no contact and no wear.

### Example 4.2: The induced $\mathbf E$ of a ramping solenoid (project Part 2)

$n = 2000$/m, $R = 2$ cm, $dI/dt = 100$ A/s, so $\mu_0n\,dI/dt = 0.2513$ T/s. By (3.11), $|E_\phi| = 1.257$ mV/m at $s = 1$ cm, $2.513$ mV/m at the wall, and $1.257$ mV/m at $4$ cm: linear inside, $1/s$ outside, the profile of $B$ around a uniform-current wire (EM-5 §3.6), as the magnetostatics analogy predicts. With $dI/dt > 0$ it circulates clockwise seen from $+z$.

*A harder check.* A circle of radius $3$ cm centred $2$ cm off axis encloses only part of the solenoid, and $\mathbf E$ is not tangent to it, so symmetry gives nothing. Faraday still says the circulation is $\mu_0n\,(dI/dt)$ times the overlap area of the two disks, $9.320\times10^{-4}$ m². The project integrates $\mathbf E\cdot d\mathbf l$ around the circle numerically: $234.2492\ \mu$V, equal to the flux rate $234.2492\ \mu$V.

### Example 4.3: Mutual inductance, two ways (project Part 3)

**(a) Coaxial loops.** Radii $R_1 = 5$ cm, $R_2 = 3$ cm, separation $z$. If loop 2 is small compared with the distance, it sits in loop 1's on-axis field $\mu_0I_1R_1^2/2(R_1^2 + z^2)^{3/2}$ (EM-5 (3.13)), so
$$M\approx\frac{\mu_0\pi R_1^2R_2^2}{2(R_1^2 + z^2)^{3/2}}\ \xrightarrow{z\gg R_1}\ \frac{\mu_0\pi R_1^2R_2^2}{2z^3},$$
the last form being symmetric in the two loops, as reciprocity requires. The exact value (Maxwell's elliptic-integral form of (3.12)) and the project's numerical Neumann integral agree to ten digits:

| $z$ | exact = Neumann | first approximation | $\mu_0\pi R_1^2R_2^2/2z^3$ |
|---|---|---|---|
| $10$ cm | $2.933\times10^{-9}$ H | $3.178\times10^{-9}$ (+8%) | $4.441\times10^{-9}$ (+51%) |
| $50$ cm | $3.482\times10^{-11}$ H | $3.500\times10^{-11}$ (+0.5%) | $3.553\times10^{-11}$ (+2.0%) |

The dipole formula needs $z$ several times the larger radius (`figures/mutual_coax.png`).

**(b) A small loop inside a solenoid.** Solenoid $n = 2000$/m, $R = 2$ cm, length $0.5$ m; small loop of radius $r = 5$ mm at its centre, normal tilted $\theta = 30°$ from the axis.

*Way 1 (current $I_1$ in the solenoid).* $\mathbf B_1 = \mu_0nI_1\hat{\mathbf z}$ has flux $\mu_0nI_1\pi r^2\cos\theta$ through the loop: $M_{21} = \mu_0n\pi r^2\cos\theta = 1.7095\times10^{-7}$ H.

*Way 2 (current $I_2$ in the small loop).* Count the small loop's flux through the solenoid's turns. Only its axial moment $m_z = I_2\pi r^2\cos\theta$ contributes (the transverse moment's flux through each coaxial turn cancels by symmetry). The dipole's $A_\phi = \mu_0m_z\sin\vartheta/4\pi\mathscr r^2$ (EM-5 (3.32)) gives, around a turn of radius $R$ at height $z$, $\oint\mathbf A\cdot d\mathbf l = 2\pi RA_\phi = \mu_0m_zR^2/2(R^2 + z^2)^{3/2}$. Summing over $n\,dz$ turns of an infinite solenoid,
$$\Phi_1 = n\int_{-\infty}^\infty\frac{\mu_0m_zR^2\,dz}{2(R^2 + z^2)^{3/2}} = \mu_0nm_z = \mu_0n\pi r^2\cos\theta\,I_2,$$
using $\int_{-\infty}^{\infty}R^2\,dz/(R^2 + z^2)^{3/2} = 2$. So $M_{12} = M_{21}$ from two calculations that share nothing. For the finite solenoid the project does both independently (way 1: exact elliptic fields of 1000 rings integrated over the small disk; way 2: the loop integral of the small loop's own $\mathbf A$ around each turn) and gets $1.70402\times10^{-7}$ H both times, $0.32\%$ below the infinite value: the finite solenoid's centre-field factor $(\ell/2)/\sqrt{(\ell/2)^2 + R^2} = 0.99682$.

### Example 4.4: The charging capacitor (project Part 4)

Circular plates of radius $a = 5$ cm, gap $d = 1$ mm, charged by $I = 2$ A through axial wires; ignore fringing.

*Displacement current.* $E = Q/\varepsilon_0\pi a^2$, so $\partial_tE = I/\varepsilon_0\pi a^2 = 2.876\times10^{13}$ V/(m s), and the total displacement current is $\varepsilon_0\partial_tE\cdot\pi a^2 = I$.

*$\mathbf B$ in the gap.* A coaxial circle of radius $s < a$ encloses displacement current $Is^2/a^2$; (3.19) gives $2\pi sB = \mu_0Is^2/a^2$, so $B = \mu_0Is/2\pi a^2$: $3.200\ \mu$T at $2$ cm, $8.000\ \mu$T at $s = a$. For $s > a$, $B = \mu_0I/2\pi s$, e.g. $4.000\ \mu$T at $10$ cm, the same as around the wire there. A bag surface on a circle around the wire that passes between the plates catches displacement current $I$ and gives the same answer; without Maxwell's term the two surfaces would disagree by $\mu_0I$.

*Energy flow.* At $V = 100$ V ($C = \varepsilon_0\pi a^2/d = 69.54$ pF, $E = 10^5$ V/m), at the rim $\mathbf E$ is axial and $\mathbf B$ azimuthal, so $\mathbf S$ points radially inward with $|\mathbf S| = EB(a)/\mu_0 = EI/2\pi a = 6.366\times10^5$ W/m². Through the band $2\pi ad$: $EId = VI = 200.0$ W $= d(\tfrac12CV^2)/dt$. The energy enters through the open sides of the gap, not along the wires.

### Example 4.5: Energy enters a resistor through its sides (project Part 5)

A nichrome wire (resistivity $1.10\times10^{-6}\ \Omega$ m), radius $a = 0.5$ mm, length $L = 1$ m, current $I = 2$ A: $R = 1.4006\ \Omega$, $V = 2.8011$ V, $P = I^2R = 5.6023$ W.

Inside, $\mathbf E = (V/L)\hat{\mathbf z} = 2.8011$ V/m; tangential $\mathbf E$ is continuous (3.23), so the same $E_z$ sits just outside. $B(a) = \mu_0I/2\pi a = 0.800$ mT, azimuthal. Since $\hat{\mathbf z}\times\hat{\boldsymbol\phi} = -\hat{\mathbf s}$, $\mathbf S$ points *into* the wire with
$$|\mathbf S(a)| = \frac{EB}{\mu_0} = \frac{EI}{2\pi a} = 1783.3\ \text{W/m}^2,$$
and over the side area $2\pi aL$ the inflow is $EIL = VI = 5.6023$ W $= I^2R$. The ends contribute nothing ($\mathbf S$ is radial). Inside, $B = \mu_0Is/2\pi a^2$, and the inflow through radius $s$ is $EIs^2L/a^2$, the Joule heat generated within $s$: energy soaks in and is consumed uniformly, $\nabla\cdot\mathbf S = -\mathbf E\cdot\mathbf J = -7.133\times10^6$ W/m³. The project builds $\mathbf S$ on a grid (`figures/poynting_resistor.png`), integrates its flux through a *square* around the wire (so symmetry does no work), gets $5.6023$ W, and checks $\nabla\cdot\mathbf S$ by finite differences. In a coaxial cable the energy travels in the insulator (Exercise B5).

### Example 4.6: Coulomb's law as a flux of stress (project Part 6)

Charges $q$ at $(0, 0, \pm d)$. Find the force on the upper one from fields on the midplane alone.

Enclose the upper charge by the midplane and a hemisphere at infinity (stresses fall as $r^{-4}$, area grows as $r^2$: no contribution). The outward normal on the midplane is $-\hat{\mathbf z}$, so for static fields $F_z = -\int T_{zz}\,da$.

*Like charges.* On the plane the $z$-components cancel: $\mathbf E = \frac{2qs}{4\pi\varepsilon_0(s^2 + d^2)^{3/2}}\hat{\mathbf s}$, and $T_{zz} = -\tfrac12\varepsilon_0E^2$, a pressure. With $w = s^2 + d^2$, $\int_0^\infty\frac{s^3\,ds}{(s^2 + d^2)^3} = \int_{d^2}^\infty\frac{w - d^2}{2w^3}dw = \frac1{4d^2}$, so
$$F_z = \frac{\varepsilon_0}{2}\int_0^\infty\left[\frac{2qs}{4\pi\varepsilon_0(s^2 + d^2)^{3/2}}\right]^22\pi s\,ds = \frac{q^2}{4\pi\varepsilon_0}\int_0^\infty\frac{s^3\,ds}{(s^2 + d^2)^3} = \frac{q^2}{4\pi\varepsilon_0(2d)^2},$$
upward: Coulomb repulsion at separation $2d$.

*Unlike charges.* Now the $\hat{\mathbf s}$ components cancel, $|E_z| = \frac{2qd}{4\pi\varepsilon_0(s^2 + d^2)^{3/2}}$, $T_{zz} = +\tfrac12\varepsilon_0E_z^2$ is a tension, and with $\int_0^\infty s\,ds/(s^2 + d^2)^3 = 1/4d^4$ the force is $-q^2/4\pi\varepsilon_0(2d)^2$: attraction. Like charges push through the pressure of lines lying in the plane; unlike charges pull through the tension of lines crossing it.

*Numbers.* $q = 1\ \mu$C, $d = 1$ cm: $22.4689$ N. The project builds the full $T_{ij}$ on a polar grid over the plane ($s = d\tan\vartheta$ maps it to a finite interval), contracts with $d\mathbf a$, and gets $\pm22.4689$ N for $F_z$ and zero for $F_x$, $F_y$.

## 5. Why it matters and where it connects

- **EM-7.** (3.24) becomes plane waves and polarization; (3.23) at interfaces gives Fresnel's equations (the project's $r = -1/3$, $t = 2/3$ at $n = 2$ is the normal-incidence case); Drude's model with a restoring force gives dispersion; (3.31) gives radiation pressure.
- **EM-8.** Retarded solutions of (3.37); radiated power is $\oint\mathbf S\cdot d\mathbf a$; $\mathbf E$ and $\mathbf B$ become the field tensor $F_{\mu\nu}$, and $u$, $\mathbf S$, $\mathbf g$, $T_{ij}$ one $4\times4$ stress–energy tensor with Lagrangian (3.38).
- **QM-1, QM-4, QC-1.** Light's momentum-to-energy ratio $1/c$ (§3.10); gauge invariance as local phase freedom (§3.13); flux quantization in superconducting qubit hardware.
- **CM-4, LA-12.** Poynting's theorem and (3.31) are Noether's theorem for translations of (3.38) (CM-4 §3.10). $T_{ij}$ is a symmetric $(0,2)$ tensor acting as Cauchy's traction (LA-12 §3.18).
- **Engineering.** Generators, transformers and wireless chargers (mutual inductance), induction motors and eddy-current brakes (Lenz), and FDTD itself, a workhorse of antenna and photonics design.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: Drude numbers, a loop in an ac field, solenoid energy, a tiny displacement current, gauge choices, stress values. B1–B8 are the module: B1 (Faraday's disk), B3 (induced $\mathbf E$ and $-\partial_t\mathbf A$), B4 (reciprocity and $M^2\le L_1L_2$), B5 (coaxial cable: $L$ two ways and Poynting), B6 (displacement current and the boundary conditions) and B8 (`project.py`: Courant and numerical dispersion) are essential; B2 (falling loop) and B7 (stress on capacitor plates) drill energy and momentum bookkeeping. C1–C3: the betatron condition, a lossy FDTD medium, Maxwell–Ampère from the Lagrangian.

## 7. Coding project

`project.py` runs in about 12 seconds, figures in `figures/`. Parts 1–6 check Examples 4.1–4.6 (with **Neumann's formula done numerically**, SymPy proofs of (3.28) and Theorem 3.8 step 6, and **the midplane stress integral**); Part 7 does **gauge invariance on a grid** and the Lorenz-gauge dipole; Part 8 is **1D FDTD** (`yee1d`); Part 9 a **2D TE FDTD** point source (`figures/fdtd_2d.png`); then Your turn (lossy medium, Faraday's disk, two coils, with reference solutions) and every number in `solutions.md`.

**The Yee scheme and the Courant condition.** For $\mathbf E = E_y(x, t)\hat{\mathbf y}$, $\mathbf B = B_z(x, t)\hat{\mathbf z}$ the curl equations are $\partial_tB_z = -\partial_xE_y$ and $\varepsilon\partial_tE_y = -\mu_0^{-1}\partial_xB_z$. Yee (1966) stores $E$ at grid points $x_i$ and integer times $t_n$, and $B$ at half points and half times, so every derivative is a centred difference (derived in the project's markdown cell). Why the time step is limited:

1. Insert $E^n_i = E_0e^{i(kx_i - \omega t_n)}$ and $B^{n+1/2}_{i+1/2} = B_0e^{i(kx_{i+1/2} - \omega t_{n+1/2})}$ into the $B$ update $B^{n+1/2}_{i+1/2} - B^{n-1/2}_{i+1/2} = -\frac{\Delta t}{\Delta x}(E^n_{i+1} - E^n_i)$. Each difference of exponentials is a sine: $B_0\sin\frac{\omega\Delta t}{2} = \frac{\Delta t}{\Delta x}E_0\sin\frac{k\Delta x}{2}$.
2. The $E$ update gives likewise $E_0\sin\frac{\omega\Delta t}{2} = \frac{c^2\Delta t}{\Delta x}B_0\sin\frac{k\Delta x}{2}$.
3. Multiply: $\sin\frac{\omega\Delta t}{2} = S\sin\frac{k\Delta x}{2}$, with **Courant number** $S = c\Delta t/\Delta x$.
4. If $S\le1$ every $k$ has real $\omega$. If $S > 1$ the shortest grid wave ($k\Delta x = \pi$) needs $\sin(\omega\Delta t/2) = S > 1$: $\omega$ is complex and that mode grows every step. ∎

At $S = 1$, $\omega = ck$ exactly in 1D; at $S = 0.5$ waves are slightly slow (the measured pulse speed is $2.3\times10^{-4}$ below $c$); at $S = 1.02$ round-off in the shortest mode grows to $6\times10^{87}$ in 600 steps. In 2D the limit is $S\le1/\sqrt2$. The scheme also conserves a discrete energy $\tfrac12\sum[\varepsilon(E^n)^2 + B^{n-1/2}B^{n+1/2}/\mu_0]\Delta x$ exactly (to $10^{-15}$ in the run), and that quantity is positive definite only for $S\le1$: the same condition seen through Poynting's theorem.

The runs: a Gaussian pulse moves at $c$; a perfect conductor ($E = 0$) reflects it inverted; at an interface to $n = 2$ the reflected and transmitted amplitudes are $-0.3331$ and $0.6667$ (theory $-1/3$, $2/3$), energy fractions $0.1110$ and $0.8890$; with absorbing (Mur) ends, $\int S_x\,dt$ through a plane near the end equals the pulse energy.

**Build it yourself:** before reading `yee1d`, write the two-line leapfrog for a vacuum box with $E = 0$ at both ends, launch $E = e^{-(x - x_0)^2/w^2}$ with $B = E/c$ shifted by half a time step, and check that the pulse moves right at $c$ and returns inverted. Then add the dielectric.

## 8. Where you'll struggle

- **The flux rule has real exceptions.** Faraday's disk has EMF without flux change; Feynman's rocking plates have flux change without EMF. Feynman *Lectures* II §17-2 is required reading. When in doubt compute $\oint(\mathbf E + \mathbf v\times\mathbf B)\cdot d\mathbf l$ with $\mathbf v$ the velocity of the conducting material; that never fails.
- **Signs in Faraday's law.** The minus sign means something only after the loop direction and the surface normal are tied by the right-hand rule. Fix the orientation, compute $\Phi$ with that normal, read the direction of $\mathcal E$ from the sign, and check against Lenz.
- **Poynting's theorem: everyone loses a sign.** The traps are (P4)'s minus and the outward $d\mathbf a$. Derive it twice, the second time closed-book, and test on Example 4.5: if $\mathbf S$ points out of the resistor, you lost one.
- **The stress tensor's indices.** This is the first time a tensor does real work in E&M, and §3.10 looks like a wall of indices. LA-12 is the cure: read $T_{ij}\,da_j$ as matrix times area vector, and step 5 as the product rule backwards. Then use tension along lines and pressure across them to get signs without computing.
- **The displacement current is not a current.** No charge crosses a vacuum gap; $\varepsilon_0\partial_t\mathbf E$ is a changing field that sources $\mathbf B$ as a current would. The polarization current $\partial_t\mathbf P$ in a dielectric *is* moving charge.
- **Lorenz, not Lorentz.** The gauge is Ludvig Lorenz's; many books printed it as "Lorentz" for decades.
- **FDTD stability.** Halve $\Delta x$ without halving $\Delta t$ and the code explodes: the Courant condition, not your update equations. Derive $\sin(\omega\Delta t/2) = S\sin(k\Delta x/2)$ once (§7) and always set `dt = S * dx / c`. In a dielectric the wave is slower, so the vacuum $S$ is the binding one; in 2D and 3D the limits are $1/\sqrt2$ and $1/\sqrt3$.
- **What does a voltmeter read?** With changing flux nearby, $\int\mathbf E\cdot d\mathbf l$ depends on the path, so the reading depends on where the leads run (Exercise A2(c)).

## 9. Check yourself

1. Derive $\sigma = nq^2\tau/m$. Why does $\mathbf B$ not contribute to Joule heating?
2. Prove the flux rule for a moving loop in static $\mathbf B$. Where is $\nabla\cdot\mathbf B = 0$ used?
3. Give the EMF of Faraday's disk and say why the flux rule seems not to apply.
4. Turn the integral form of Faraday's law into $\nabla\times\mathbf E = -\partial_t\mathbf B$. What does $-\partial_t\mathbf B$ play the role of?
5. Find $E(s)$ inside and outside a solenoid with ramping current. Why is $\mathbf E\ne0$ where $\mathbf B = 0$?
6. Derive Neumann's formula from $\Phi = \oint\mathbf A\cdot d\mathbf l$. Why does it prove $M_{12} = M_{21}$?
7. Get from $\tfrac12LI^2$ to $\int B^2/2\mu_0$. Which product rule, and which surface term vanishes?
8. Show that Ampère's law contradicts charge conservation and derive the displacement current.
9. Derive the polarization current and write Maxwell's equations in matter.
10. Derive the vacuum wave equation and the value of $c$.
11. Derive Poynting's theorem. Which way does energy flow into a resistor?
12. Write $T_{ij}$ and say what $T_{ij}\,da_j$ means; why do like charges repel in the tension/pressure picture? Prove that (3.34) leaves $\mathbf E$ and $\mathbf B$ unchanged and derive $\Box V = -\rho/\varepsilon_0$ in the Lorenz gauge.

## 10. Pair with

- **Griffiths, *Introduction to Electrodynamics*, ch. 7–8.** The same order and derivations (§7.1.3's flux rule, §7.3's displacement current, §8.1–8.2's Poynting theorem and stress tensor), with more worked problems and better pictures.
- **Purcell and Morin, *Electricity and Magnetism*, ch. 7 and 9.** More physical insight per page, and an electromagnetic wave built from a current sheet switched on at $t = 0$, which makes "the curls sustain each other" concrete.
- **Feynman, *Lectures on Physics* Vol. II, ch. 16–18 and 27.** Ch. 17 §2 (exceptions to the flux rule) is required; ch. 27 is the best discussion of where field energy and momentum go, including the disk paradox.
- **MIT 8.07 lecture notes, lectures 15–19.** A compact, more mathematical second pass through the conservation laws, the stress tensor and gauge choices.
