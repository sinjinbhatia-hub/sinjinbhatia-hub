# EM-6: Electrodynamics: Faraday, Maxwell's Equations, Energy and Momentum, Gauge Freedom

**Time:** ~12 hours · **Prerequisites:** EM-5 (and through it EM-1's product rules and EM-4's bound charges) · **Unlocks:** EM-7 (waves), EM-8 (radiation, relativity), QM-1 (photons carry energy and momentum), CM-4 (Noether's theorem for fields)

> **Conventions in this module.** SI units; $\mu_0\varepsilon_0 = 1/c^2$. EMF is $\mathcal E$ (volts), flux is $\Phi = \int_{\mathcal S}\mathbf B\cdot d\mathbf a$ (webers). **Orientation:** a loop $\mathcal P$ and any surface $\mathcal S$ spanning it are related by the right-hand rule (fingers along $d\mathbf l$, thumb along $d\mathbf a$), exactly as in EM-1's Stokes theorem; a positive $\mathcal E$ drives current along $d\mathbf l$. Faraday's law is $\nabla\times\mathbf E = -\partial\mathbf B/\partial t$; Maxwell–Ampère is $\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\,\partial\mathbf E/\partial t$. Energy density $u = \tfrac12(\varepsilon_0E^2 + B^2/\mu_0)$, Poynting vector $\mathbf S = \mathbf E\times\mathbf B/\mu_0$, momentum density $\mathbf g = \varepsilon_0\mathbf E\times\mathbf B = \mathbf S/c^2$. **Stress tensor with Griffiths' sign:** $T_{ij} = \varepsilon_0(E_iE_j - \tfrac12\delta_{ij}E^2) + \mu_0^{-1}(B_iB_j - \tfrac12\delta_{ij}B^2)$, so that $T_{ij}\,da_j$ is the force *on* the contents of a surface (tension positive), and $-T$ is the momentum flux. Repeated Latin indices are summed. Potentials: $\mathbf E = -\nabla V - \partial_t\mathbf A$, $\mathbf B = \nabla\times\mathbf A$; gauge change $V\to V - \partial_t\lambda$, $\mathbf A\to\mathbf A + \nabla\lambda$; **Lorenz** gauge (Ludvig Lorenz) $\nabla\cdot\mathbf A + c^{-2}\partial_tV = 0$; d'Alembertian $\Box\equiv\nabla^2 - c^{-2}\partial_t^2$ (Griffiths writes $\Box^2$). Field point $\mathbf r$, source point $\mathbf r'$, $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$. Surface jumps are "above minus below" with $\hat{\mathbf n}$ pointing up (EM-2 §3.7).

## 0. Where this sits

EM-2 to EM-5 were statics: $\nabla\cdot\mathbf E = \rho/\varepsilon_0$, $\nabla\times\mathbf E = 0$, $\nabla\cdot\mathbf B = 0$, $\nabla\times\mathbf B = \mu_0\mathbf J$. You need EM-5's vector potential and its loop integral $\oint\mathbf A\cdot d\mathbf l = \Phi$ (EM-5 (3.25)), EM-4's $\rho_b = -\nabla\cdot\mathbf P$, and EM-1's product rule (P4) and curl-of-curl (3.13). This module adds the two time-derivative terms, which turn four static laws into one dynamical theory with waves, field energy and field momentum. EM-7 solves it for waves; EM-8 solves the potential equations of §3.13 with retarded potentials and packages everything into one tensor.

## 1. The problem that forced this idea

Ørsted (1820) showed that a current makes a magnetic field. Michael Faraday spent years looking for the converse and found it on 29 August 1831: two coils wound on an iron ring, and a galvanometer on the second coil kicked when the current in the first was *switched* on or off, not while it flowed. Over the following weeks he produced the same effect by pushing a magnet into a coil and by spinning a copper disk between magnet poles. The law was about *change*. Joseph Henry in Albany found induction independently and published in 1832, including self-induction; Heinrich Lenz stated the sign rule in 1834; Franz Neumann wrote the mutual inductance of two circuits as a double line integral in 1845.

Static theory could not hold once things move. Ampère's law $\nabla\times\mathbf B = \mu_0\mathbf J$ says, after taking a divergence, that $\nabla\cdot\mathbf J = 0$: charge can never pile up anywhere. But a capacitor charges. EM-5 Exercise B1(c) found the symptom: two surfaces with the same rim give different "enclosed currents". James Clerk Maxwell, in *On Physical Lines of Force* (1861–62), added a "displacement current" $\varepsilon_0\partial\mathbf E/\partial t$ to Ampère's law. He reached it through a mechanical model of the ether, but in the *Treatise* (1873) the requirement he stated is the one we will use: the total current, conduction plus displacement, always flows in closed loops. With the new term the equations had wave solutions travelling at $1/\sqrt{\mu_0\varepsilon_0}$. Wilhelm Weber and Rudolf Kohlrausch had measured exactly that combination of electrical constants in 1856 (about $3.1\times10^8$ m/s), and Hippolyte Fizeau had measured the speed of light in 1849 (also about $3.1\times10^8$ m/s). Maxwell wrote in 1862: "we can scarcely avoid the inference that light consists in the transverse undulations of the same medium which is the cause of electric and magnetic phenomena."

Heinrich Hertz generated and detected the waves with spark-gap circuits in 1887–88. John Henry Poynting (1884) asked where the energy goes when a battery drives a resistor, and found the uncomfortable answer that it flows through the empty space around the wire. Oliver Heaviside rewrote Maxwell's twenty equations as the four vector equations we use. Ludvig Lorenz had already, in 1867, written the potentials in the gauge now named for him. The forcing problem in one sentence: the static equations contradict each other as soon as anything moves, and the fix unified electricity, magnetism and optics.

## 2. Intuition

**A changing $\mathbf B$ makes $\mathbf E$ curl around it.** Push a magnet into a ring of wire: current flows. Hold the magnet still and move the ring the other way: the same current. In the first case nothing in the ring moves, so no magnetic force acts on its charges; the push must be an electric field, and one that circulates (a conservative field cannot drive charge around a loop). Its circulation equals the rate of change of the flux through the ring. In EM-5 a current made $\mathbf B$ curl around it with circulation $\mu_0 I$. Guess the twin: **$\oint\mathbf E\cdot d\mathbf l = -d\Phi/dt$**, so $-\partial\mathbf B/\partial t$ plays the role of $\mu_0\mathbf J$. The minus sign is Lenz: the induced current makes flux that fights the change. With the opposite sign a sliding bar on rails would push itself faster and faster.

```
     changing B (out of page, growing)        induced E circulates clockwise
              . . .                                  <----
            . . . . .                              |       ^
            . . . . .                              v       |
              . . .                                  ---->
```

**A changing $\mathbf E$ makes $\mathbf B$ curl around it.** Draw a loop around the wire feeding a charging capacitor. A flat disk on that loop is pierced by the current $I$. A bag-shaped surface on the same loop that passes between the plates is pierced by nothing, yet the loop and its $\oint\mathbf B\cdot d\mathbf l$ are the same. What does the bag catch? A growing electric field: the plate charge $Q$ grows at rate $I$, so the flux of $\varepsilon_0\mathbf E$ through the bag grows at rate $I$ too (Gauss). Guess: $\varepsilon_0\,\partial\mathbf E/\partial t$ counts as current. That is Maxwell's term.

**Two curls make a wave.** Now there is a loop that needs no charges. A changing $\mathbf B$ makes a curling $\mathbf E$; that $\mathbf E$ changes, so it makes a curling $\mathbf B$; and so on. The disturbance carries itself through empty space. How fast? The only speed you can build from $\mu_0$ (N/A²) and $\varepsilon_0$ (C²/N m²) is $1/\sqrt{\mu_0\varepsilon_0}$. Guess that number before §3.8 computes it.

**Energy lives in the field and flows through it.** A battery drives current through a resistor. Inside the wire $\mathbf E$ points along the wire; $\mathbf B$ circles it. The only direction you can build from both that is not along either is $\mathbf E\times\mathbf B$, which points *radially inward*. If energy flows anywhere it flows along $\mathbf E\times\mathbf B$, so it should enter the resistor through its sides. And if a field carries energy at speed $c$, relativity's $E = pc$ says it carries momentum $u/c$ per volume: guess $\mathbf g = \mathbf S/c^2$. Faraday pictured field lines as stretched strings that also push sideways on their neighbours; the stress tensor of §3.11 makes that literal.

## 3. Formal development

### 3.1 Ohm's law, the Drude model, Joule heating

To keep a current flowing in a real conductor you must push. For most conductors the current density is proportional to the force per unit charge $\mathbf f$:
$$\mathbf J = \sigma\mathbf f\ \xrightarrow{\ \text{fields at rest}\ }\ \mathbf J = \sigma\mathbf E. \tag{3.1}$$
$\sigma$ is the conductivity (S/m), $1/\sigma$ the resistivity. This is a **constitutive relation**, a statement about a material, like $\mathbf P = \varepsilon_0\chi_e\mathbf E$ in EM-4; it is not a law of nature.

**Drude's justification (1900).** Treat conduction electrons (density $n$, charge $q$, mass $m$) as a gas that is accelerated by $\mathbf E$ and loses its directed momentum in collisions at mean interval $\tau$: $m\,d\mathbf v/dt = q\mathbf E - m\mathbf v/\tau$. In steady state $\mathbf v = q\tau\mathbf E/m$ and $\mathbf J = nq\mathbf v$, so
$$\sigma = \frac{nq^2\tau}{m}. \tag{3.2}$$
For copper ($n = 8.49\times10^{28}$ m$^{-3}$, $\sigma = 5.96\times10^7$ S/m) this gives $\tau = 2.5\times10^{-14}$ s, and a current of 1 A in a 1 mm² wire has drift speed $7.4\times10^{-5}$ m/s (project Part 1). The classical picture gets the order of magnitude right and the temperature dependence wrong; the quantum theory of metals (SM-7) keeps (3.2) with a different meaning for $\tau$.

**Resistance and Joule heating.** For a uniform wire of length $L$ and cross-section $A$ carrying $I$: $J = I/A$, $E = J/\sigma$, $V = EL$, so $V = IR$ with $R = L/(\sigma A)$. The power the field delivers to the charges follows from the force law alone: a charge $dq$ moving with $\mathbf v$ receives power $dq\,(\mathbf E + \mathbf v\times\mathbf B)\cdot\mathbf v = dq\,\mathbf E\cdot\mathbf v$ (the magnetic part does no work, EM-5 Theorem 3.1). Summing, with $dq\,\mathbf v = \mathbf J\,d\tau$,
$$P = \int\mathbf E\cdot\mathbf J\,d\tau = (EL)(JA) = VI = I^2R. \tag{3.3}$$
In a resistor that power ends as heat.

### 3.2 EMF and motional EMF: the flux rule

Two forces drive a circuit: the source force $\mathbf f_s$ (chemistry in a battery, a moving conductor in a generator) and the electrostatic $\mathbf E$ that the source sets up to distribute the charge. The **electromotive force** is the line integral of the total force per unit charge around the loop,
$$\mathcal E = \oint\mathbf f\cdot d\mathbf l. \tag{3.4}$$
The static $\mathbf E$ has zero circulation (EM-2), so $\mathcal E = \oint\mathbf f_s\cdot d\mathbf l$. For a single resistive loop, $\mathcal E = IR$. "EMF" is not a force; it is work per unit charge around the loop.

**Motional EMF.** Move a conducting loop through a static $\mathbf B$. A charge in the wire moves with the wire's velocity $\mathbf v$ plus a drift $\mathbf u$ along the wire, and feels $q(\mathbf v + \mathbf u)\times\mathbf B$. Since $\mathbf u\parallel d\mathbf l$, $(\mathbf u\times\mathbf B)\cdot d\mathbf l = 0$, and only $\mathbf v\times\mathbf B$ contributes to $\mathcal E$.

▶ **Try it first.** A rectangular loop with one side of length $\ell$ sliding at speed $v$ (the bar on rails, Example 4.1) in uniform $\mathbf B$ perpendicular to the loop. Compute $\oint(\mathbf v\times\mathbf B)\cdot d\mathbf l$ and $d\Phi/dt$ and compare. Then try the general loop.

**Theorem 3.1 (Flux rule).** Let a thin wire loop $\mathcal P(t)$ move and deform in a time-independent $\mathbf B$, each element $d\mathbf l$ moving with velocity $\mathbf v$. Then
$$\mathcal E = \oint(\mathbf v\times\mathbf B)\cdot d\mathbf l = -\frac{d\Phi}{dt}. \tag{3.5}$$

*Proof.*
1. In time $dt$ the element $d\mathbf l$ sweeps a parallelogram with vector area $(\mathbf v\,dt)\times d\mathbf l$. (Check the orientation: a loop in the $xy$-plane traversed counterclockwise and expanding has, on its right edge, $\mathbf v = \hat{\mathbf x}$, $d\mathbf l\propto\hat{\mathbf y}$, and $\hat{\mathbf x}\times\hat{\mathbf y} = +\hat{\mathbf z}$: added area with the loop's own orientation.)
2. The old surface $\mathcal S(t)$, the swept ribbon, and the new surface $\mathcal S(t + dt)$ bound a volume. Since $\nabla\cdot\mathbf B = 0$, the net flux out of it is zero, so with the orientations of step 1, $\Phi(t + dt) = \Phi(t) + \int_{\text{ribbon}}\mathbf B\cdot d\mathbf a$.
3. Hence $d\Phi/dt = \oint\mathbf B\cdot(\mathbf v\times d\mathbf l)$.
4. The triple product is cyclic: $\mathbf B\cdot(\mathbf v\times d\mathbf l) = d\mathbf l\cdot(\mathbf B\times\mathbf v) = -(\mathbf v\times\mathbf B)\cdot d\mathbf l$. So $d\Phi/dt = -\oint(\mathbf v\times\mathbf B)\cdot d\mathbf l = -\mathcal E$. ∎

■ **What just happened.** The flux rule for a moving loop is not a new law. It is the magnetic Lorentz force plus $\nabla\cdot\mathbf B = 0$, which makes "the flux through the loop" independent of which spanning surface you choose. Which force does the work? Not the magnetic one, which never does. In the sliding bar, $\mathbf v\times\mathbf B$ pushes charges along the bar, and as they then move along the bar the magnetic force on them gains a component opposing the bar's motion; whoever pulls the bar does the work (Example 4.1 balances it).

**Where the flux rule needs care.** The proof assumed a thin wire whose material moves with the curve. **Faraday's disk** breaks this. A conducting disk of radius $a$ spins at $\omega$ in a uniform $\mathbf B = B\hat{\mathbf z}$, with brushes on the axle and the rim. The material velocity is $\mathbf v = \omega s\hat{\boldsymbol\phi}$, so $\mathbf v\times\mathbf B = \omega Bs\,\hat{\mathbf s}$, and along any path from axle to rim
$$\mathcal E = \int_0^a\omega Bs\,ds = \tfrac12\omega Ba^2. \tag{3.6}$$
(The integrand $\omega Bs\,\hat{\mathbf s}$ is $\nabla(\tfrac12\omega Bs^2)$, so the path does not matter; the project integrates along a radius and along a spiral.) Draw the circuit as a fixed radius through the disk and its flux never changes, yet there is a steady EMF. Feynman's rocking plates (*Lectures* II §17-2) do the reverse: two curved plates rock slightly so that the contact point, and with it the circuit, sweeps across a large field region; the flux through "the circuit" changes a lot and the EMF is nearly zero. The rule that never fails is the Lorentz force: $\mathcal E = \oint(\mathbf E + \mathbf v\times\mathbf B)\cdot d\mathbf l$, with $\mathbf v$ the velocity of the material the current actually flows through.

### 3.3 Faraday's law

Faraday's experimental result: when the loop is at rest and the magnet moves (or a nearby current changes), the galvanometer reads the same $\mathcal E = -d\Phi/dt$ as when the loop moves. Einstein opened his 1905 relativity paper with this coincidence: the two explanations differ (magnetic force in one frame, electric in the other), the measured current does not.

▶ **Try it first.** For a stationary loop the charges have $\mathbf v = 0$. What force drives them? Turn $\mathcal E = -d\Phi/dt$ into a differential equation using Stokes.

**Theorem 3.2 (Faraday's law).** If for every fixed loop $\mathcal P$ with spanning surface $\mathcal S$
$$\oint_{\mathcal P}\mathbf E\cdot d\mathbf l = -\frac{d}{dt}\int_{\mathcal S}\mathbf B\cdot d\mathbf a = -\int_{\mathcal S}\frac{\partial\mathbf B}{\partial t}\cdot d\mathbf a, \tag{3.7}$$
then
$$\nabla\times\mathbf E = -\frac{\partial\mathbf B}{\partial t}. \tag{3.8}$$

*Proof.*
1. The loop is at rest, so no magnetic force acts; the force per unit charge that drives current is $\mathbf E$, and $\mathcal E = \oint\mathbf E\cdot d\mathbf l$. The experiment says this equals $-d\Phi/dt$; for a fixed surface the time derivative passes inside the integral. That is (3.7).
2. Stokes (EM-1 (3.19)) turns the left side into $\int_{\mathcal S}(\nabla\times\mathbf E)\cdot d\mathbf a$.
3. So $\int_{\mathcal S}(\nabla\times\mathbf E + \partial_t\mathbf B)\cdot d\mathbf a = 0$ for every surface. Shrink $\mathcal S$ to a small disk with normal $\hat{\mathbf n}$ around any point: the integrand's $\hat{\mathbf n}$ component vanishes there, for every $\hat{\mathbf n}$. ∎

■ **What just happened.** Electrostatics' $\nabla\times\mathbf E = 0$ was a special case. A changing magnetic field produces an electric field with curl, and its line integral around a loop is not zero, so in electrodynamics "the potential difference between two points" depends on the path. A voltmeter's reading depends on how its leads are routed around a changing flux (Exercise A2).

**The induced field is a magnetostatics problem.** Split $\mathbf E$ into the part made by charges and the induced part. Where $\rho = 0$, the induced field obeys
$$\nabla\cdot\mathbf E = 0,\quad\nabla\times\mathbf E = -\partial_t\mathbf B\qquad\longleftrightarrow\qquad\nabla\cdot\mathbf B = 0,\quad\nabla\times\mathbf B = \mu_0\mathbf J.$$

| Magnetostatics (EM-5) | Faraday-induced $\mathbf E$ |
|---|---|
| source $\mu_0\mathbf J$ | source $-\partial\mathbf B/\partial t$ |
| $\oint\mathbf B\cdot d\mathbf l = \mu_0I_{\rm enc}$ | $\oint\mathbf E\cdot d\mathbf l = -d\Phi/dt$ |
| Biot–Savart $\dfrac{\mu_0}{4\pi}\displaystyle\int\frac{\mathbf J\times\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}d\tau'$ | $-\dfrac{1}{4\pi}\displaystyle\int\frac{\partial_t\mathbf B\times\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}d\tau'$ |

$$\mathbf E_{\rm induced}(\mathbf r) = -\frac{1}{4\pi}\int\frac{\partial_t\mathbf B(\mathbf r')\times\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}\,d\tau'. \tag{3.9}$$
Every Ampère-law trick of EM-5 §3.6 now finds induced electric fields. One caveat: (3.9) treats $\partial_t\mathbf B$ as known everywhere at the same instant. That is the **quasistatic** approximation, valid when the system is small compared with $c$ times the time over which things change; §3.6's correction and EM-8's retardation fix it otherwise.

**The universal flux rule.** For a loop that moves through a time-dependent field, $d\Phi/dt$ has two parts: $\int\partial_t\mathbf B\cdot d\mathbf a$ from the field changing (Theorem 3.2) and $\oint\mathbf B\cdot(\mathbf v\times d\mathbf l)$ from the loop moving (Theorem 3.1, steps 1–3, applied at each instant). Adding,
$$\mathcal E = \oint(\mathbf E + \mathbf v\times\mathbf B)\cdot d\mathbf l = -\frac{d\Phi}{dt}, \tag{3.10}$$
with the same exceptions as in §3.2.

**The solenoid with a ramping current.** A long solenoid (radius $R$, $n$ turns per length) carries $I(t)$, so inside $\mathbf B = \mu_0nI\hat{\mathbf z}$ and outside $\mathbf B = 0$ (EM-5 (3.21)). By symmetry $\mathbf E = E_\phi(s)\hat{\boldsymbol\phi}$ (a radial component would need charge; an axial one would make $\oint\mathbf E\cdot d\mathbf l\ne0$ around a rectangle where $\partial_t\mathbf B = 0$). Take a circle of radius $s$ oriented counterclockwise about $+\hat{\mathbf z}$: $2\pi sE_\phi = -d\Phi/dt$, with $\Phi = \mu_0nI\pi\min(s, R)^2$:
$$E_\phi = -\frac{\mu_0ns}{2}\frac{dI}{dt}\ \ (s < R),\qquad E_\phi = -\frac{\mu_0nR^2}{2s}\frac{dI}{dt}\ \ (s > R). \tag{3.11}$$
Outside, $\mathbf B$ is zero and never changes, yet $\mathbf E\ne0$. Compare EM-5 (3.26): $A_\phi = \mu_0nIs/2$ inside and $\mu_0nIR^2/2s$ outside, so (3.11) is exactly $\mathbf E = -\partial\mathbf A/\partial t$. §3.12 shows why.

### 3.4 Inductance

Current $I_1$ in loop 1 makes $\mathbf B_1\propto I_1$ (Biot–Savart is linear), so the flux of $\mathbf B_1$ through loop 2 is $\Phi_2 = M_{21}I_1$. The constant $M_{21}$ (henries, Wb/A) is the **mutual inductance**. If $I_1$ changes, $\mathcal E_2 = -M_{21}\,dI_1/dt$.

▶ **Try it first.** Write $\Phi_2$ as a loop integral of $\mathbf A_1$ around loop 2 (EM-5 (3.25)), then insert the wire form of EM-5 (3.24). Look at the result for a symmetry.

**Theorem 3.3 (Neumann formula; reciprocity).**
$$M_{21} = \frac{\mu_0}{4\pi}\oint_1\oint_2\frac{d\mathbf l_1\cdot d\mathbf l_2}{|\mathbf r_1 - \mathbf r_2|} = M_{12}. \tag{3.12}$$

*Proof.*
1. Stokes: $\Phi_2 = \int_{\mathcal S_2}(\nabla\times\mathbf A_1)\cdot d\mathbf a_2 = \oint_2\mathbf A_1\cdot d\mathbf l_2$.
2. For a thin wire in the Coulomb gauge, $\mathbf A_1(\mathbf r) = \frac{\mu_0I_1}{4\pi}\oint_1\frac{d\mathbf l_1}{|\mathbf r - \mathbf r_1|}$ (EM-5 (3.24)).
3. Insert: $\Phi_2 = \frac{\mu_0I_1}{4\pi}\oint_2\oint_1\frac{d\mathbf l_1\cdot d\mathbf l_2}{|\mathbf r_2 - \mathbf r_1|}$; divide by $I_1$.
4. The integrand is unchanged by swapping the labels 1 and 2, so the same expression is $M_{12}$. ∎

■ **What just happened.** The flux through loop 2 per ampere in loop 1 equals the flux through loop 1 per ampere in loop 2, however different the loops are. From the $\mathbf B$ picture this is not obvious at all: in Example 4.3(b) one side is "a uniform field through a small loop" and the other is "a dipole's flux threading a thousand turns". Step 1 holds in any gauge, because a gauge term $\nabla\lambda$ has zero loop integral. $M$ depends only on geometry.

**Self-inductance.** A changing current also changes the flux through its own loop: $\Phi = LI$ and
$$\mathcal E = -L\frac{dI}{dt}. \tag{3.13}$$
$L > 0$ always (the back-EMF opposes the change). Setting loop 1 = loop 2 in (3.12) diverges logarithmically for an infinitely thin wire, so $L$ is computed from the flux of a wire of finite size. A long solenoid of length $\ell$, $N = n\ell$ turns, cross-section $\pi R^2$: each turn links $\mu_0nI\pi R^2$, so $\Phi = N\mu_0nI\pi R^2$. A toroid of $N$ turns with rectangular cross-section (inner radius $a$, outer $b$, height $h$) has $B = \mu_0NI/2\pi s$ (EM-5 (3.22)); one turn links $\int_a^b\frac{\mu_0NI}{2\pi s}h\,ds$. So
$$L_{\rm sol} = \mu_0n^2\pi R^2\ell,\qquad L_{\rm tor} = \frac{\mu_0N^2h}{2\pi}\ln\frac ba. \tag{3.14}$$
For $n = 2000$/m, $R = 2$ cm, $\ell = 0.5$ m: $L_{\rm sol} = 3.158$ mH. For $N = 500$, $h = 1$ cm, $a = 4$ cm, $b = 6$ cm: $L_{\rm tor} = 0.2027$ mH (project Part 3).

### 3.5 Energy in the magnetic field

▶ **Try it first.** To raise the current in a loop from $0$ to $I$, the source must push against the back-EMF (3.13). Compute the total work. Then rewrite $LI$ as $\oint\mathbf A\cdot d\mathbf l$.

**Theorem 3.4 (Magnetic energy).** The work needed to establish a current distribution $\mathbf J$ (starting from zero, slowly, with no losses) is
$$W = \tfrac12LI^2\ \ (\text{one loop}), \tag{3.15}$$
$$W = \tfrac12\int\mathbf A\cdot\mathbf J\,d\tau, \tag{3.16}$$
$$W = \frac{1}{2\mu_0}\int_{\text{all space}}B^2\,d\tau. \tag{3.17}$$

*Proof.* **Route 1 (the circuit and the back-EMF).**
1. The source does work per unit charge $-\mathcal E = L\,dI/dt$ on charge flowing at rate $I$: $dW/dt = LI\,dI/dt$.
2. Integrate from $0$ to $I$: $W = \tfrac12LI^2$. That is (3.15).
3. $LI = \Phi = \oint\mathbf A\cdot d\mathbf l$, so $W = \tfrac12I\oint\mathbf A\cdot d\mathbf l$. For a distribution, $I\,d\mathbf l\to\mathbf J\,d\tau$, giving (3.16). (The same result follows directly at field level: the source does work $-\int\mathbf E\cdot\mathbf J\,d\tau = \int\partial_t\mathbf A\cdot\mathbf J\,d\tau$ per unit time against the induced field $\mathbf E = -\partial_t\mathbf A$; raise every current by the same factor $\alpha$ from $0$ to $1$, so $\mathbf A\propto\alpha$ and the integral over $\alpha$ supplies the $\tfrac12$.)

**Route 2 (the product rule).**
4. Replace $\mathbf J$ by $\nabla\times\mathbf B/\mu_0$ (steady currents): $W = \frac{1}{2\mu_0}\int\mathbf A\cdot(\nabla\times\mathbf B)\,d\tau$.
5. By (P4) of EM-1, $\nabla\cdot(\mathbf A\times\mathbf B) = \mathbf B\cdot(\nabla\times\mathbf A) - \mathbf A\cdot(\nabla\times\mathbf B)$, so $\mathbf A\cdot(\nabla\times\mathbf B) = B^2 - \nabla\cdot(\mathbf A\times\mathbf B)$.
6. $W = \frac{1}{2\mu_0}\left[\int_{\mathcal V}B^2\,d\tau - \oint_{\mathcal S}(\mathbf A\times\mathbf B)\cdot d\mathbf a\right]$ for any volume $\mathcal V$ containing all the current.
7. Let $\mathcal V$ be a ball of radius $r\to\infty$. Far away $A\sim r^{-2}$ and $B\sim r^{-3}$ (dipole), so the surface term is $\sim r^{-5}\cdot r^2\to0$. That is (3.17). ∎

■ **What just happened.** The same energy is $\tfrac12\int\mathbf A\cdot\mathbf J$, an integral over where the current is, and $\int B^2/2\mu_0$, an integral over where the field is: the exact twin of EM-2's $\tfrac12\int\rho V\,d\tau = \tfrac{\varepsilon_0}{2}\int E^2\,d\tau$ (EM-2 (3.15)–(3.16)). Check on the solenoid: $B^2/2\mu_0$ times the volume $\pi R^2\ell$ is $\tfrac12\mu_0n^2I^2\pi R^2\ell = \tfrac12L_{\rm sol}I^2$. For the toroid the project integrates $B^2/2\mu_0$ over the cross-section and recovers $\tfrac12L_{\rm tor}I^2$.

**Electric versus magnetic storage.** At $B = 1$ T, $u_B = B^2/2\mu_0 = 3.98\times10^5$ J/m³. Air breaks down near $E = 3$ MV/m, where $u_E = \tfrac12\varepsilon_0E^2 = 39.8$ J/m³. Matching the magnetic density would need $E = cB = 3\times10^8$ V/m. That factor of $10^4$ is why large energies are stored in magnets (a clinical MRI magnet holds megajoules) and why capacitors need thin gaps and good dielectrics.

### 3.6 Maxwell's correction

▶ **Try it first.** Take the divergence of $\nabla\times\mathbf B = \mu_0\mathbf J$ and compare with the continuity equation (EM-5 Theorem 3.2). Then use Gauss's law to find a vector whose divergence is $-\nabla\cdot\mathbf J$.

**Theorem 3.5 (The displacement current).** Ampère's law $\nabla\times\mathbf B = \mu_0\mathbf J$ is inconsistent with charge conservation whenever $\partial_t\rho\ne0$. The modification
$$\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\frac{\partial\mathbf E}{\partial t} \tag{3.18}$$
is consistent with both continuity and Gauss's law.

*Derivation.*
1. The divergence of a curl vanishes (EM-1 (3.12)), so Ampère's law implies $\mu_0\nabla\cdot\mathbf J = 0$.
2. Continuity says $\nabla\cdot\mathbf J = -\partial_t\rho$. Steps 1 and 2 agree only if $\partial_t\rho = 0$: Ampère's law is a law of *magnetostatics*, and nothing more.
3. Gauss's law $\rho = \varepsilon_0\nabla\cdot\mathbf E$, differentiated: $\partial_t\rho = \nabla\cdot(\varepsilon_0\partial_t\mathbf E)$.
4. So $\nabla\cdot(\mathbf J + \varepsilon_0\partial_t\mathbf E) = 0$ always. The vector $\mathbf J + \varepsilon_0\partial_t\mathbf E$ is divergence-free and reduces to $\mathbf J$ in statics, so put it where $\mathbf J$ was. That is (3.18), and its divergence now reproduces continuity instead of contradicting it. ∎

■ **What just happened.** The term $\mathbf J_d\equiv\varepsilon_0\partial_t\mathbf E$ is the **displacement current density**. It is the minimal fix consistent with Gauss's law; nothing in steps 1–4 proves it is the *only* fix. Experiment decides: the waves of §3.8 exist and move at the predicted speed. In integral form, with $\Phi_E = \int\mathbf E\cdot d\mathbf a$,
$$\oint_{\mathcal P}\mathbf B\cdot d\mathbf l = \mu_0I_{\rm enc} + \mu_0\varepsilon_0\frac{d\Phi_E}{dt}, \tag{3.19}$$
and now any two surfaces on the same loop give the same right-hand side (by step 4 and the divergence theorem). Example 4.4 works the capacitor with numbers. The term is invisible in circuits at low frequency because $\varepsilon_0$ is tiny (Exercise A4); it dominates in empty space, where it is the only source of $\mathbf B$.

### 3.7 Maxwell's equations, in vacuum and in matter

$$\boxed{\ \nabla\cdot\mathbf E = \frac{\rho}{\varepsilon_0},\qquad\nabla\cdot\mathbf B = 0,\qquad\nabla\times\mathbf E = -\frac{\partial\mathbf B}{\partial t},\qquad\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\frac{\partial\mathbf E}{\partial t}\ } \tag{3.20}$$
With the Lorentz force $\mathbf F = q(\mathbf E + \mathbf v\times\mathbf B)$ and Newton's second law, this is all of classical electrodynamics. Here $\rho$ and $\mathbf J$ are *all* charges and currents. By Helmholtz's theorem (EM-1 Theorem 3.15) a field that vanishes at infinity is fixed by its divergence and curl, and (3.20) supplies both for each field.

**In matter.** EM-4 split charge into free and bound, $\rho = \rho_f + \rho_b$ with $\rho_b = -\nabla\cdot\mathbf P$; EM-5 split steady current into $\mathbf J_f + \mathbf J_b$ with $\mathbf J_b = \nabla\times\mathbf M$. A *changing* polarization adds a third current.

▶ **Try it first.** A slab with polarization $\mathbf P(t)$ has surface bound charge $\sigma_b = \mathbf P\cdot\hat{\mathbf n}$. If $P$ grows, charge must move to the surface. Find the current density that does it, and check continuity with $\rho_b$.

*Derivation of the polarization current.*
1. $\mathbf P$ is dipole moment per volume. Increasing it by $d\mathbf P$ shifts bound charge so that $d\mathbf P\cdot d\mathbf a$ crosses an area element $d\mathbf a$ (EM-4 §3.5's tilted-surface count). The current density is therefore
$$\mathbf J_p = \frac{\partial\mathbf P}{\partial t}. \tag{3.21}$$
2. Check: $\nabla\cdot\mathbf J_p = \partial_t(\nabla\cdot\mathbf P) = -\partial_t\rho_b$, which is continuity for bound charge.
3. $\mathbf J_b = \nabla\times\mathbf M$ has zero divergence, so it never violates continuity and needs no partner. ∎

Now insert $\rho = \rho_f - \nabla\cdot\mathbf P$ and $\mathbf J = \mathbf J_f + \nabla\times\mathbf M + \partial_t\mathbf P$ into (3.20). Gauss: $\nabla\cdot(\varepsilon_0\mathbf E + \mathbf P) = \rho_f$. Maxwell–Ampère: $\nabla\times(\mathbf B/\mu_0 - \mathbf M) = \mathbf J_f + \partial_t(\varepsilon_0\mathbf E + \mathbf P)$. With $\mathbf D = \varepsilon_0\mathbf E + \mathbf P$ (EM-4) and $\mathbf H = \mathbf B/\mu_0 - \mathbf M$ (EM-5),
$$\nabla\cdot\mathbf D = \rho_f,\qquad\nabla\cdot\mathbf B = 0,\qquad\nabla\times\mathbf E = -\frac{\partial\mathbf B}{\partial t},\qquad\nabla\times\mathbf H = \mathbf J_f + \frac{\partial\mathbf D}{\partial t}. \tag{3.22}$$
These need constitutive relations ($\mathbf D = \varepsilon\mathbf E$, $\mathbf H = \mathbf B/\mu$, $\mathbf J_f = \sigma\mathbf E$ for linear media) before they can be solved.

**Boundary conditions.** At a surface carrying free charge $\sigma_f$ and free current $\mathbf K_f$:
$$D^\perp_{\rm above} - D^\perp_{\rm below} = \sigma_f,\quad B^\perp_{\rm above} = B^\perp_{\rm below},\quad\mathbf E^\parallel_{\rm above} = \mathbf E^\parallel_{\rm below},\quad\mathbf H^\parallel_{\rm above} - \mathbf H^\parallel_{\rm below} = \mathbf K_f\times\hat{\mathbf n}. \tag{3.23}$$
The divergence equations did not change, so the first two are EM-4 (3.16) and EM-5 (3.28) word for word. The curl equations gained time derivatives; we must check the last two survive.

*Derivation of the two curl conditions.*
1. Take a rectangle straddling the surface: two sides of length $\ell$ parallel to the surface (one just above, one just below), two sides of height $\epsilon$. Apply the integral form of Faraday's law (3.7): $\oint\mathbf E\cdot d\mathbf l = -\frac{d}{dt}\int\mathbf B\cdot d\mathbf a$.
2. As $\epsilon\to0$ the short sides contribute nothing (finite $\mathbf E$), and the flux through the rectangle is at most $\ell\epsilon\max|\partial_t\mathbf B|\to0$ (finite $\partial_t\mathbf B$).
3. What remains is $(\mathbf E^\parallel_{\rm above} - \mathbf E^\parallel_{\rm below})\cdot\hat{\boldsymbol\ell}\,\ell = 0$ for every direction $\hat{\boldsymbol\ell}$ in the surface: $\mathbf E^\parallel$ is continuous.
4. Repeat with $\oint\mathbf H\cdot d\mathbf l = I_{f,\rm enc} + \frac{d}{dt}\int\mathbf D\cdot d\mathbf a$. The $\mathbf D$-flux term vanishes as in step 2. The free current through the shrinking rectangle does not vanish if there is a surface current: it is $\ell\,\mathbf K_f\cdot(\hat{\mathbf n}\times\hat{\boldsymbol\ell})$, which is exactly the static case of EM-5 (3.42), giving the last condition. ∎

■ **What just happened.** Displacement current and $\partial_t\mathbf B$ are fluxes of finite fields, and a flux through a vanishing area vanishes. Only *surface* sources ($\sigma_f$, $\mathbf K_f$) cause jumps. EM-7 uses (3.23) at every interface; the project's dielectric interface obeys it by construction.

### 3.8 The wave equation and the speed of light

▶ **Try it first.** In vacuum ($\rho = 0$, $\mathbf J = 0$) take the curl of Faraday's law. Use EM-1 (3.13) on the left and Maxwell–Ampère on the right.

**Theorem 3.6 (Waves).** In vacuum every component of $\mathbf E$ and $\mathbf B$ satisfies
$$\nabla^2\mathbf E = \mu_0\varepsilon_0\frac{\partial^2\mathbf E}{\partial t^2},\qquad\nabla^2\mathbf B = \mu_0\varepsilon_0\frac{\partial^2\mathbf B}{\partial t^2}. \tag{3.24}$$

*Derivation.*
1. Curl of Faraday: $\nabla\times(\nabla\times\mathbf E) = -\partial_t(\nabla\times\mathbf B)$ (space and time derivatives commute).
2. Left: $\nabla(\nabla\cdot\mathbf E) - \nabla^2\mathbf E = -\nabla^2\mathbf E$, since $\nabla\cdot\mathbf E = 0$ in vacuum.
3. Right: $-\partial_t(\mu_0\varepsilon_0\partial_t\mathbf E) = -\mu_0\varepsilon_0\partial_t^2\mathbf E$. Equate: the first of (3.24).
4. Curl of Maxwell–Ampère: $-\nabla^2\mathbf B = \mu_0\varepsilon_0\partial_t(\nabla\times\mathbf E) = -\mu_0\varepsilon_0\partial_t^2\mathbf B$. ∎

■ **What just happened.** Each equation in (3.24) is the wave equation $\nabla^2f = v^{-2}\partial_t^2f$ with speed
$$c = \frac{1}{\sqrt{\mu_0\varepsilon_0}} = 2.998\times10^8\ \text{m/s}. \tag{3.25}$$
Two constants measured with charged spheres and current balances give the speed of light. That was Maxwell's inference, and it is why step 3 needs the displacement current: without it there is no $\partial_t^2$ and no wave. (Since 2019 the SI fixes $c$ and $\mu_0$ is measured, so (3.25) now holds by construction; in 1862 it was a discovery.)

The coupled equations say more than (3.24) alone. Try $\mathbf E = f(x - ct)\hat{\mathbf y}$. Faraday: $\partial_tB_z = -\partial_xE_y = -f'$, so $B_z = f(x - ct)/c$. The magnetic field is perpendicular to $\mathbf E$ and to the direction of travel, with $|\mathbf B| = |\mathbf E|/c$, and $\mathbf E\times\mathbf B$ points along the motion. The project launches its FDTD pulse from exactly this pair; EM-7 develops plane waves fully.

**Magnetic charge (one paragraph).** Nothing in (3.20) forbids a magnetic charge density $\rho_m$ and current $\mathbf J_m$: $\nabla\cdot\mathbf B = \mu_0\rho_m$ and $\nabla\times\mathbf E = -\mu_0\mathbf J_m - \partial_t\mathbf B$ would make the equations symmetric. In vacuum the equations are already invariant under the **duality rotation** $\mathbf E\to c\mathbf B$, $c\mathbf B\to-\mathbf E$ (check it on (3.20) with no sources). No monopole has been found. Dirac (1931) showed that if even one exists, quantum mechanics requires every electric charge $q$ to satisfy $q\,\Phi_m = nh$ ($n$ an integer), where $\Phi_m$ is the monopole's total magnetic flux. The argument is EM-5's Aharonov–Bohm phase applied to the invisible flux tube that must carry the monopole's flux in from infinity: the phase $q\Phi_m/\hbar$ must be a multiple of $2\pi$. One monopole anywhere would explain why all charges are integer multiples of $e$.

### 3.9 Conservation laws: the continuity equation as template

EM-5 Theorem 3.2 is the pattern for every conservation law in field theory. A density $\rho$ and a current $\mathbf J$ satisfy $\partial_t\rho + \nabla\cdot\mathbf J = 0$; integrated over a fixed volume, the amount inside changes only by flow through the boundary. When the quantity can also be exchanged with something else, the right side gains a source term. Energy and momentum of the fields follow this pattern, with the source being the exchange with charged matter.

### 3.10 Poynting's theorem

▶ **Try it first.** Start with the power delivered to charges, $\int\mathbf E\cdot\mathbf J\,d\tau$ from (3.3). Eliminate $\mathbf J$ with Maxwell–Ampère. You will need (P4) and Faraday's law. Every student loses a sign here at least once; do it on paper.

**Theorem 3.7 (Poynting).** Let
$$u = \frac12\left(\varepsilon_0E^2 + \frac{1}{\mu_0}B^2\right),\qquad\mathbf S = \frac{1}{\mu_0}\mathbf E\times\mathbf B. \tag{3.26}$$
Then for any fixed volume $\mathcal V$ with boundary $\mathcal S$, the work done per unit time on the charges inside is
$$\frac{dW_{\rm mech}}{dt} = \int_{\mathcal V}\mathbf E\cdot\mathbf J\,d\tau = -\frac{d}{dt}\int_{\mathcal V}u\,d\tau - \oint_{\mathcal S}\mathbf S\cdot d\mathbf a, \tag{3.27}$$
and pointwise
$$\frac{\partial u}{\partial t} + \nabla\cdot\mathbf S = -\mathbf E\cdot\mathbf J. \tag{3.28}$$

*Derivation.*
1. Power on charges: the magnetic force does no work, so $dW_{\rm mech}/dt = \int\mathbf E\cdot\mathbf J\,d\tau$ (§3.1).
2. Maxwell–Ampère: $\mathbf J = \mu_0^{-1}\nabla\times\mathbf B - \varepsilon_0\partial_t\mathbf E$. So $\mathbf E\cdot\mathbf J = \mu_0^{-1}\mathbf E\cdot(\nabla\times\mathbf B) - \varepsilon_0\mathbf E\cdot\partial_t\mathbf E$.
3. (P4) with $\mathbf A = \mathbf E$: $\nabla\cdot(\mathbf E\times\mathbf B) = \mathbf B\cdot(\nabla\times\mathbf E) - \mathbf E\cdot(\nabla\times\mathbf B)$. So $\mathbf E\cdot(\nabla\times\mathbf B) = \mathbf B\cdot(\nabla\times\mathbf E) - \nabla\cdot(\mathbf E\times\mathbf B)$.
4. Faraday: $\mathbf B\cdot(\nabla\times\mathbf E) = -\mathbf B\cdot\partial_t\mathbf B = -\tfrac12\partial_t(B^2)$. Also $\mathbf E\cdot\partial_t\mathbf E = \tfrac12\partial_t(E^2)$.
5. Assemble: $\mathbf E\cdot\mathbf J = -\frac{1}{2\mu_0}\partial_t(B^2) - \frac{1}{\mu_0}\nabla\cdot(\mathbf E\times\mathbf B) - \frac{\varepsilon_0}{2}\partial_t(E^2) = -\partial_tu - \nabla\cdot\mathbf S$. That is (3.28).
6. Integrate over $\mathcal V$ and use the divergence theorem on the last term: (3.27). ∎

■ **What just happened.** Read (3.27) right to left: the energy delivered to matter equals the decrease of field energy inside plus the energy flowing *in* through the surface ($-\oint\mathbf S\cdot d\mathbf a$). So $\mathbf S$ is the energy flux density (W/m²), and (3.28) is the continuity equation of §3.9 with the source $-\mathbf E\cdot\mathbf J$: energy handed from field to charges. The project checks (3.28) symbolically for fields built from random polynomial potentials. Two cautions. First, only $\oint\mathbf S\cdot d\mathbf a$ is fixed by the theorem; adding any divergence-free field to $\mathbf S$ changes nothing in (3.27), so "the energy flows exactly along $\mathbf S$ at each point" is an interpretation (the standard one, and the one that matches radiation). Second, a static charge next to a magnet has $\mathbf S\ne0$ circulating forever with $\nabla\cdot\mathbf S = 0$: harmless for energy, but it carries momentum and angular momentum (§3.11). For a plane wave, $|\mathbf S| = E^2/\mu_0c = c\,u$: the energy moves at $c$.

### 3.11 Momentum in the fields: the Maxwell stress tensor

The force on the charges in a volume is $\mathbf F = \int(\rho\mathbf E + \mathbf J\times\mathbf B)\,d\tau$. Newton's third law fails for magnetic forces between moving charges (EM-5's two current elements pull unequally), so if momentum is conserved, the field must hold some. We find out how much by eliminating the sources, exactly as in §3.10. This is long; every step is a product rule.

▶ **Try it first.** Write the force density $\mathbf f = \rho\mathbf E + \mathbf J\times\mathbf B$ with $\rho = \varepsilon_0\nabla\cdot\mathbf E$ and $\mathbf J$ from Maxwell–Ampère. Try to write the $\mathbf E$ terms as $\partial_j(\text{something})_{ij}$.

**Theorem 3.8 (Momentum balance).** With $T_{ij}$ as in the conventions box,
$$T_{ij} = \varepsilon_0\left(E_iE_j - \tfrac12\delta_{ij}E^2\right) + \frac{1}{\mu_0}\left(B_iB_j - \tfrac12\delta_{ij}B^2\right), \tag{3.29}$$
the force density on charges is $f_i = \partial_jT_{ij} - \varepsilon_0\mu_0\,\partial_tS_i$, and the total force on the charges in $\mathcal V$ is
$$\mathbf F = \oint_{\mathcal S}\overleftrightarrow{T}\cdot d\mathbf a - \varepsilon_0\mu_0\frac{d}{dt}\int_{\mathcal V}\mathbf S\,d\tau. \tag{3.30}$$

*Derivation.*
1. $\mathbf f = \varepsilon_0(\nabla\cdot\mathbf E)\mathbf E + \mu_0^{-1}(\nabla\times\mathbf B)\times\mathbf B - \varepsilon_0(\partial_t\mathbf E)\times\mathbf B$.
2. Move the time derivative off $\mathbf E$: $\partial_t(\mathbf E\times\mathbf B) = (\partial_t\mathbf E)\times\mathbf B + \mathbf E\times\partial_t\mathbf B = (\partial_t\mathbf E)\times\mathbf B - \mathbf E\times(\nabla\times\mathbf E)$ (Faraday). So $-\varepsilon_0(\partial_t\mathbf E)\times\mathbf B = -\varepsilon_0\partial_t(\mathbf E\times\mathbf B) - \varepsilon_0\mathbf E\times(\nabla\times\mathbf E)$.
3. Add $\mu_0^{-1}(\nabla\cdot\mathbf B)\mathbf B$, which is zero, to make the expression symmetric:
$$\mathbf f = \varepsilon_0\big[(\nabla\cdot\mathbf E)\mathbf E - \mathbf E\times(\nabla\times\mathbf E)\big] + \frac{1}{\mu_0}\big[(\nabla\cdot\mathbf B)\mathbf B - \mathbf B\times(\nabla\times\mathbf B)\big] - \varepsilon_0\partial_t(\mathbf E\times\mathbf B).$$
4. In index notation, using EM-1 (3.1): $[\mathbf E\times(\nabla\times\mathbf E)]_i = \varepsilon_{ijk}\varepsilon_{klm}E_j\partial_lE_m = E_j\partial_iE_j - E_j\partial_jE_i = \tfrac12\partial_i(E^2) - E_j\partial_jE_i$.
5. So the $\mathbf E$ bracket is $E_i\partial_jE_j + E_j\partial_jE_i - \tfrac12\partial_i(E^2) = \partial_j(E_iE_j) - \partial_j(\tfrac12\delta_{ij}E^2)$, the product rule read backwards. The $\mathbf B$ bracket is the same with $\mathbf B$.
6. Hence $f_i = \partial_jT_{ij} - \varepsilon_0\partial_t(\mathbf E\times\mathbf B)_i = \partial_jT_{ij} - \varepsilon_0\mu_0\partial_tS_i$. Integrate over $\mathcal V$ and apply the divergence theorem to each row $i$: (3.30). ∎

■ **What just happened.** $T_{ij}$ is a symmetric $3\times3$ array that eats an area vector and returns a force: $F_i = \oint T_{ij}\,da_j$. That is exactly Cauchy's traction $\mathbf t = \sigma\mathbf n$ of LA-12 §3.18, so $T$ is a symmetric $(0,2)$ tensor (equivalently $(1,1)$ in an orthonormal frame) and $T_{ij}$ is the $i$-th component of force per area on a surface whose normal is along $j$. The project checks the identity in step 6 symbolically. For a field $\mathbf E = E\hat{\mathbf x}$ alone, $T = \tfrac12\varepsilon_0E^2\,\mathrm{diag}(1, -1, -1)$: a **tension** $\tfrac12\varepsilon_0E^2$ along the field lines and a **pressure** of the same size across them. That is Faraday's picture of stretched, mutually repelling lines, made quantitative; at a conductor it is EM-2's surface pressure (3.17). At $1$ MV/m the stress is $4.43$ Pa; at $1$ T the magnetic stress $B^2/2\mu_0$ is $3.98\times10^5$ Pa, about 3.9 atmospheres, which is why high-field magnets need steel bandages.

**Field momentum.** Newton's law for the matter in $\mathcal V$ is $\mathbf F = d\mathbf p_{\rm mech}/dt$. With (3.30),
$$\frac{d\mathbf p_{\rm mech}}{dt} = -\frac{d}{dt}\int_{\mathcal V}\mathbf g\,d\tau + \oint_{\mathcal S}\overleftrightarrow{T}\cdot d\mathbf a,\qquad\mathbf g = \varepsilon_0\mathbf E\times\mathbf B = \frac{\mathbf S}{c^2},\qquad\frac{\partial\mathbf g}{\partial t} + \nabla\cdot(-\overleftrightarrow{T}) = -\mathbf f. \tag{3.31}$$
$\mathbf g$ is the field's momentum density and $-T$ its momentum flux: the same template as (3.28), one equation per component. If no field crosses $\mathcal S$, $\mathbf p_{\rm mech} + \int\mathbf g\,d\tau$ is constant. For light, momentum per energy is $g/u = S/(c^2u) = 1/c$, the classical form of the photon's $p = E/c$ (QM-1).

**Angular momentum.** By the same logic the field carries angular momentum density $\boldsymbol\ell = \mathbf r\times\mathbf g$. Static configurations can hold it. Feynman's **disk paradox** (*Lectures* II §17-4, resolved in §27-6): a disk free to turn on an axle carries a coaxial solenoid with a steady current and a ring of charged spheres near its rim. Switch the current off. The falling flux induces a circulating $\mathbf E$ (3.11) that torques the charges, and the disk starts to turn. Nothing external acted, so where did the angular momentum come from? It was stored all along in the static $\varepsilon_0\mathbf r\times(\mathbf E\times\mathbf B)$ of the charges' $\mathbf E$ crossed with the solenoid's $\mathbf B$; the total is conserved.

### 3.12 The potentials in electrodynamics

$\nabla\cdot\mathbf B = 0$ still holds, so $\mathbf B = \nabla\times\mathbf A$ as in EM-5. But $\nabla\times\mathbf E\ne0$, so $\mathbf E = -\nabla V$ fails.

▶ **Try it first.** Insert $\mathbf B = \nabla\times\mathbf A$ into Faraday's law and find a combination of $\mathbf E$ and $\mathbf A$ whose curl vanishes.

*Derivation.*
1. $\nabla\times\mathbf E = -\partial_t(\nabla\times\mathbf A)$, so $\nabla\times(\mathbf E + \partial_t\mathbf A) = 0$.
2. A curl-free field is a gradient (EM-1 Theorem 3.16): $\mathbf E + \partial_t\mathbf A = -\nabla V$. Hence
$$\mathbf E = -\nabla V - \frac{\partial\mathbf A}{\partial t},\qquad\mathbf B = \nabla\times\mathbf A. \tag{3.32}$$
3. The two homogeneous equations ($\nabla\cdot\mathbf B = 0$, Faraday) now hold identically. Insert (3.32) into the two with sources. Gauss: $-\nabla^2V - \partial_t(\nabla\cdot\mathbf A) = \rho/\varepsilon_0$.
4. Maxwell–Ampère: $\nabla\times(\nabla\times\mathbf A) = \mu_0\mathbf J - \mu_0\varepsilon_0\nabla(\partial_tV) - \mu_0\varepsilon_0\partial_t^2\mathbf A$. Use EM-1 (3.13) on the left and collect:
$$\nabla^2V + \frac{\partial}{\partial t}(\nabla\cdot\mathbf A) = -\frac{\rho}{\varepsilon_0},\qquad\left(\nabla^2\mathbf A - \frac{1}{c^2}\frac{\partial^2\mathbf A}{\partial t^2}\right) - \nabla\left(\nabla\cdot\mathbf A + \frac{1}{c^2}\frac{\partial V}{\partial t}\right) = -\mu_0\mathbf J. \tag{3.33}$$
∎

■ **What just happened.** Six field components obeying eight equations became four potentials obeying four equations. The price is that (3.33) is coupled and ugly. The cure is the freedom in the potentials. The solenoid of §3.3 already showed (3.32) at work: with $\rho = 0$ and $V = 0$, $\mathbf E = -\partial_t\mathbf A$.

### 3.13 Gauge transformations; the Coulomb and Lorenz gauges

**Theorem 3.9 (Gauge invariance).** For any smooth $\lambda(\mathbf r, t)$, the potentials
$$\mathbf A' = \mathbf A + \nabla\lambda,\qquad V' = V - \frac{\partial\lambda}{\partial t} \tag{3.34}$$
give the same $\mathbf E$ and $\mathbf B$ as $(V, \mathbf A)$.

*Proof.* $\nabla\times\mathbf A' = \nabla\times\mathbf A + \nabla\times\nabla\lambda = \mathbf B$ by EM-1 (3.11). $-\nabla V' - \partial_t\mathbf A' = -\nabla V + \nabla\partial_t\lambda - \partial_t\mathbf A - \partial_t\nabla\lambda = \mathbf E$, because mixed partial derivatives commute. ∎

Conversely, if $(V, \mathbf A)$ and $(V', \mathbf A')$ give the same fields, then $\nabla\times(\mathbf A' - \mathbf A) = 0$, so $\mathbf A' - \mathbf A = \nabla\lambda$; then $\nabla(V' - V + \partial_t\lambda) = 0$, so $V' - V + \partial_t\lambda$ depends on $t$ only and can be absorbed into $\lambda$. Equation (3.34) is the *whole* freedom. The project does the proof on a grid: central differences along different axes commute exactly, so a gauge change that moves $V$ and $\mathbf A$ by $O(1)$ moves the computed $\mathbf E$ and $\mathbf B$ by $2\times10^{-14}$.

**Coulomb gauge, $\nabla\cdot\mathbf A = 0$.** It is reachable (EM-5 §3.7: solve $\nabla^2\lambda = -\nabla\cdot\mathbf A$). The first of (3.33) becomes
$$\nabla^2V = -\frac{\rho}{\varepsilon_0},\qquad V(\mathbf r, t) = \frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\mathbf r', t)}{\mathscr r}\,d\tau'. \tag{3.35}$$
$V$ responds *instantly* to charge anywhere in the universe. That is not a violation of relativity, because $V$ is not measurable; $\mathbf E$ also contains $-\partial_t\mathbf A$, and in this gauge the equation for $\mathbf A$ is sourced by $\mathbf J - \varepsilon_0\nabla\partial_tV$, whose instantaneous part cancels the instantaneous part of $-\nabla V$ in $\mathbf E$. (What is skipped: the full demonstration that $\mathbf E$ comes out retarded uses the transverse/longitudinal split of $\mathbf J$; see Jackson, *Classical Electrodynamics*, §6.3.) The Coulomb gauge is the natural one for statics and for quantizing the field (QM, quantum optics), where only the two transverse components of $\mathbf A$ are dynamical.

▶ **Try it first.** Look at the second of (3.33). Which condition on $(V, \mathbf A)$ kills its last bracket? Insert the same condition into the first equation.

**Lorenz gauge.** Impose
$$\nabla\cdot\mathbf A + \frac{1}{c^2}\frac{\partial V}{\partial t} = 0. \tag{3.36}$$

*Derivation.*
1. **It is reachable.** Given any $(V_0, \mathbf A_0)$, transform with $\lambda$. Then $\nabla\cdot\mathbf A + c^{-2}\partial_tV = (\nabla\cdot\mathbf A_0 + c^{-2}\partial_tV_0) + \nabla^2\lambda - c^{-2}\partial_t^2\lambda$. Choose $\lambda$ with $\Box\lambda = -(\nabla\cdot\mathbf A_0 + c^{-2}\partial_tV_0)$, an inhomogeneous wave equation, which always has solutions (EM-8 writes them down).
2. Second of (3.33): the gradient term vanishes, leaving $\Box\mathbf A = -\mu_0\mathbf J$.
3. First of (3.33): $\partial_t(\nabla\cdot\mathbf A) = -c^{-2}\partial_t^2V$, so $\nabla^2V - c^{-2}\partial_t^2V = -\rho/\varepsilon_0$.
$$\Box V = -\frac{\rho}{\varepsilon_0},\qquad\Box\mathbf A = -\mu_0\mathbf J,\qquad\Box\equiv\nabla^2 - \frac{1}{c^2}\frac{\partial^2}{\partial t^2}. \tag{3.37}$$
∎

■ **What just happened.** In the Lorenz gauge each potential obeys its own wave equation, sourced by its own source, with the same operator. Four uncoupled copies of one equation. In statics, (3.37) reduces to Poisson's equations of EM-2 and EM-5. EM-8 solves (3.37) with the retarded Green's function: $V(\mathbf r, t) = \frac{1}{4\pi\varepsilon_0}\int\rho(\mathbf r', t - \mathscr r/c)/\mathscr r\,d\tau'$. The project checks with SymPy that the exact potentials of an oscillating dipole satisfy (3.36), solve (3.37) away from the origin, and give fields satisfying all four vacuum equations. The Lorenz condition leaves a residual freedom ($\lambda$ with $\Box\lambda = 0$), and it is relativistically invariant ($\partial_\mu A^\mu = 0$ in EM-8), which is why relativists prefer it. It is named for Ludvig Lorenz (Danish, 1867), not Hendrik Lorentz.

### 3.14 Gauge invariance, the Lagrangian, and Noether

Gauge invariance looks like bookkeeping. It is the deepest structural fact in the subject. In quantum mechanics the Schrödinger equation for a charged particle contains $\mathbf A$ and $V$ through $\mathbf p\to\mathbf p - q\mathbf A$ and $V$; under (3.34) it stays the same equation if the wavefunction changes by a *local* phase, $\psi\to e^{iq\lambda/\hbar}\psi$ (QM-4). Read backwards: demanding that the phase of $\psi$ be freely choosable at each point *forces* the existence of a field that transforms like $(V, \mathbf A)$. In superconducting hardware the same phase forces the flux through a superconducting ring to be quantized in units of $h/2e$, which is how SQUIDs and flux qubits work (QC).

Classically, the Lagrangian of a charge (CM-3 (4.17)), $L = \tfrac12mv^2 - qV + q\mathbf v\cdot\mathbf A$, changes under (3.34) by $q\,d\lambda/dt$, a total time derivative, which does not change the equations of motion (CM-3 §5). The fields have a Lagrangian too. The **Lagrangian density**
$$\mathcal L = \frac{\varepsilon_0}{2}E^2 - \frac{1}{2\mu_0}B^2 - \rho V + \mathbf J\cdot\mathbf A\quad\Big(= -\frac{1}{4\mu_0}F_{\mu\nu}F^{\mu\nu} - A_\mu J^\mu\ \text{in EM-8's notation}\Big), \tag{3.38}$$
with $\mathbf E$ and $\mathbf B$ written through (3.32), has Maxwell's equations as its Euler–Lagrange equations. The Gauss half takes three lines: $\mathcal L$ depends on $V$ through $-\rho V$ and through $\tfrac{\varepsilon_0}{2}|\nabla V + \partial_t\mathbf A|^2$, so the Euler–Lagrange equation $\partial\mathcal L/\partial V - \nabla\cdot\big(\partial\mathcal L/\partial(\nabla V)\big) = 0$ reads $-\rho - \nabla\cdot(\varepsilon_0(\nabla V + \partial_t\mathbf A)) = -\rho + \varepsilon_0\nabla\cdot\mathbf E = 0$. Exercise C3 does $\mathbf A$ and gets Maxwell–Ampère.

Now CM-4's Noether theorem (CM-4 §3.10, the field version) applies. The field Lagrangian does not depend explicitly on time: the conserved current is **energy**, and its continuity equation is Poynting's theorem (3.28). It does not depend on position: three conserved currents, **momentum**, with flux $-T_{ij}$, which is (3.31). It is rotation invariant: field **angular momentum**, Feynman's disk. And gauge invariance of the $-\rho V + \mathbf J\cdot\mathbf A$ term (it changes by a total derivative exactly when $\partial_t\rho + \nabla\cdot\mathbf J = 0$) is tied to **charge conservation**, the very condition that forced Maxwell's correction in §3.6. The four big results of this module are one theorem applied four times.

## 4. Worked examples

### Example 4.1: The bar on rails, and the power balance (project Part 1)

*From mechanical engineering: this is a linear generator and, run backwards, an eddy-current brake.* Two parallel rails $\ell = 0.5$ m apart are joined at one end by a resistor $R = 2\ \Omega$. A bar slides along them at $v = 3$ m/s in a uniform $B = 0.8$ T perpendicular to the plane. Neglect the resistance of rails and bar.

*EMF.* The flux through the circuit is $\Phi = B\ell x$, with $x$ the distance from resistor to bar, growing at $dx/dt = v$. By (3.5), $|\mathcal E| = B\ell v = (0.8)(0.5)(3) = 1.200$ V. The current is $I = \mathcal E/R = 0.600$ A, directed so that its own flux opposes the growth (Lenz).

*Force.* The current $I$ flows along the bar in field $B$, so the bar feels $F = I\ell B = (0.6)(0.5)(0.8) = 0.240$ N, opposing $v$ (if it aided $v$, the bar would accelerate itself with no energy source).

*Power balance.* To keep $v$ constant, the agent pulling the bar supplies $Fv = (0.240)(3) = 0.720$ W. The resistor dissipates $I^2R = (0.36)(2) = 0.720$ W. Equal, as they must be: $Fv = (B\ell)^2v^2/R = \mathcal E^2/R$.

*Coasting.* Let go of the bar (mass $m = 0.1$ kg). Now $m\,dv/dt = -(B\ell)^2v/R$, so $v = v_0e^{-t/\tau}$ with $\tau = mR/(B\ell)^2 = (0.1)(2)/(0.4)^2 = 1.250$ s. The total heat is $\int_0^\infty I^2R\,dt = \frac{(B\ell)^2}{R}v_0^2\frac{\tau}{2} = \tfrac12mv_0^2 = 0.4500$ J: all the kinetic energy, and the project's ODE integration confirms it. That is how eddy-current brakes on trains and roller coasters stop things with no contact and no wear.

### Example 4.2: The induced $\mathbf E$ of a ramping solenoid (project Part 2)

A long solenoid with $n = 2000$ turns/m and radius $R = 2$ cm carries a current rising at $dI/dt = 100$ A/s. Then $\mu_0n\,dI/dt = 0.2513$ T/s. By (3.11):

| $s$ | $1$ cm | $2$ cm ($=R$) | $4$ cm |
|---|---|---|---|
| $\lvert E_\phi\rvert$ | $\tfrac12(0.2513)(0.01) = 1.257$ mV/m | $2.513$ mV/m | $\tfrac12(0.2513)(0.02)^2/0.04 = 1.257$ mV/m |

$|E|$ grows linearly inside, peaks at the wall, and falls as $1/s$ outside: the profile of $B$ around a wire of uniform current (EM-5 Example), as the analogy table predicts. With $dI/dt > 0$ the field circulates clockwise seen from $+z$.

*A harder check.* Take a circle of radius $3$ cm centred $2$ cm off the axis. It encloses only part of the solenoid's cross-section, and $\mathbf E$ is not tangent to it, so you cannot read off its circulation from symmetry. Faraday says the circulation is still $\mu_0n\,(dI/dt)\times(\text{overlap area})$, where the overlap of the two disks is $9.320\times10^{-4}$ m². The project integrates $\mathbf E\cdot d\mathbf l$ around the circle numerically, $234.2492\ \mu$V, and the flux rate is $234.2492\ \mu$V.

### Example 4.3: Mutual inductance, two ways (project Part 3)

**(a) Two coaxial loops far apart.** Radii $R_1 = 5$ cm, $R_2 = 3$ cm, centres a distance $z$ apart on a common axis. If $z\gg R_1, R_2$, loop 2 sits in loop 1's on-axis field $B_1 = \mu_0I_1R_1^2/2(R_1^2 + z^2)^{3/2}$ (EM-5 (3.13)), nearly uniform over loop 2's area, so
$$M\approx\frac{\mu_0\pi R_1^2R_2^2}{2(R_1^2 + z^2)^{3/2}}\approx\frac{\mu_0\pi R_1^2R_2^2}{2z^3}.$$
The second form is the dipole–dipole result: by reciprocity it is also the flux of loop 2's dipole field through loop 1. The exact answer (Maxwell's, from (3.12) in elliptic integrals) and the project's numerical Neumann double integral agree to ten digits:

| $z$ | exact = Neumann | $\mu_0\pi R_1^2R_2^2/2(R_1^2+z^2)^{3/2}$ | $\mu_0\pi R_1^2R_2^2/2z^3$ |
|---|---|---|---|
| $10$ cm | $2.933\times10^{-9}$ H | $3.178\times10^{-9}$ (+8%) | $4.441\times10^{-9}$ (+51%) |
| $50$ cm | $3.482\times10^{-11}$ H | $3.500\times10^{-11}$ (+0.5%) | $3.553\times10^{-11}$ (+2.0%) |

The dipole formula needs $z$ several times the larger radius; `figures/mutual_coax.png` shows the approach.

**(b) A small loop inside a solenoid, both ways.** Solenoid: $n = 2000$/m, $R = 2$ cm, length $0.5$ m (1000 turns). Small loop: radius $r = 5$ mm at the centre, its normal tilted $\theta = 30°$ from the axis.

*Way 1 (current $I_1$ in the solenoid).* Inside, $\mathbf B_1 = \mu_0nI_1\hat{\mathbf z}$; its flux through the small loop is $\mu_0nI_1\pi r^2\cos\theta$, so $M_{21} = \mu_0n\pi r^2\cos\theta = 1.7095\times10^{-7}$ H.

*Way 2 (current $I_2$ in the small loop).* Now count how much of the small loop's flux threads the solenoid's turns. Only the axial part of its moment, $m_z = I_2\pi r^2\cos\theta$, contributes (the transverse part's flux through each coaxial turn cancels by symmetry). A dipole $m_z$ at the origin has $A_\phi = \mu_0m_z\sin\vartheta/(4\pi\mathscr r^2)$ (EM-5 (3.32)); around a turn of radius $R$ at height $z$, $\oint\mathbf A\cdot d\mathbf l = 2\pi R A_\phi = \mu_0m_zR^2/2(R^2 + z^2)^{3/2}$. Sum over turns ($n\,dz$ per length) on an infinitely long solenoid:
$$\Phi_1 = n\int_{-\infty}^\infty\frac{\mu_0m_zR^2}{2(R^2 + z^2)^{3/2}}dz = n\frac{\mu_0m_z}{2}\cdot2 = \mu_0nm_z = \mu_0n\pi r^2\cos\theta\,I_2.$$
So $M_{12} = M_{21}$, as Theorem 3.3 promised, from two calculations that share nothing. For the finite solenoid the project computes both ways independently (way 1 from the exact elliptic fields of 1000 rings integrated over the small disk; way 2 as the loop integral of the small loop's own $\mathbf A$ around each of the 1000 turns): $1.70402\times10^{-7}$ H both times, $0.32\%$ below the infinite-solenoid value, which is the finite solenoid's centre-field factor $(\ell/2)/\sqrt{(\ell/2)^2 + R^2} = 0.99681$.

### Example 4.4: The charging capacitor (project Part 4)

Circular plates of radius $a = 5$ cm, gap $d = 1$ mm, charged by a steady current $I = 2$ A through thin wires on the axis. Ignore fringing.

*Displacement current.* $E = Q/(\varepsilon_0\pi a^2)$ between the plates, so $\partial_tE = I/(\varepsilon_0\pi a^2) = 2.876\times10^{13}$ V/(m s), and the total displacement current $\varepsilon_0\,\partial_tE\cdot\pi a^2 = I$ exactly: the conduction current stops at the plate and continues as displacement current.

*$\mathbf B$ between the plates.* Take a circle of radius $s$ coaxial with the plates, in the gap. A flat disk on it carries no conduction current and displacement current $I s^2/a^2$ (for $s < a$). By (3.19), $2\pi sB = \mu_0Is^2/a^2$:
$$B(s) = \frac{\mu_0Is}{2\pi a^2}:\quad B(2\ \text{cm}) = 3.200\ \mu\text{T},\quad B(a) = 8.000\ \mu\text{T}.$$
For $s > a$ the circle encloses all of it and $B = \mu_0I/2\pi s$, e.g. $4.000\ \mu$T at $10$ cm: the same as around the wire at that distance. Use a bag-shaped surface on the circle around the *wire* that passes between the plates: no conduction current, displacement current $I$, same answer. Without Maxwell's term the two surfaces would disagree by $\mu_0I$.

*Energy flow.* At the instant $V = 100$ V ($C = \varepsilon_0\pi a^2/d = 69.54$ pF, $E = V/d = 10^5$ V/m), at the rim $s = a$, $\mathbf E$ is axial and $\mathbf B$ azimuthal, so $\mathbf S$ points radially *inward* with $|\mathbf S| = EB(a)/\mu_0 = EI/2\pi a = 6.366\times10^5$ W/m². Through the band of area $2\pi ad$: $EId = VI = 200.0$ W $= d(\tfrac12CV^2)/dt$. The energy enters the capacitor through the open sides of the gap, not along the wires.

### Example 4.5: Energy flows into a resistor through its sides (project Part 5)

A nichrome wire (resistivity $1.10\times10^{-6}\ \Omega$ m) of radius $a = 0.5$ mm and length $L = 1$ m carries $I = 2$ A. Then $R = 1.4006\ \Omega$, $V = 2.8011$ V, $P = I^2R = 5.6023$ W.

*Fields at the surface.* Inside, $\mathbf E = (V/L)\hat{\mathbf z} = 2.8011$ V/m along the wire; tangential $\mathbf E$ is continuous (3.23), so the same $E_z$ sits just outside. Ampère: $B(a) = \mu_0I/2\pi a = 0.800$ mT, azimuthal. $\hat{\mathbf z}\times\hat{\boldsymbol\phi} = -\hat{\mathbf s}$: $\mathbf S$ points *into* the wire, with
$$|\mathbf S(a)| = \frac{EB}{\mu_0} = \frac{EI}{2\pi a} = 1783.3\ \text{W/m}^2.$$
Over the side area $2\pi aL$: $EIL = VI = 5.6023$ W $= I^2R$. The end faces contribute nothing ($\mathbf S$ is radial).

*Inside.* At radius $s < a$, $B = \mu_0Is/2\pi a^2$ and the inward flux through a cylinder of radius $s$ is $EIs^2L/a^2$: exactly the Joule heat $\sigma E^2\pi s^2L$ generated inside it. Energy soaks in from the surface and is consumed uniformly, $\nabla\cdot\mathbf S = -\mathbf E\cdot\mathbf J = -7.133\times10^6$ W/m³. The project computes $\mathbf S$ on a grid (`figures/poynting_resistor.png`), integrates the flux through a *square* around the wire (not a circle, so the symmetry does no work for you), gets $5.6023$ W, and checks $\nabla\cdot\mathbf S$ by finite differences. The battery's energy travels through the space around the circuit, guided by the surface charges that make $\mathbf E$ outside the wire; in a coaxial cable it travels in the insulator between the conductors (Exercise B5).

### Example 4.6: Coulomb's law as a flux of stress (project Part 6)

Two charges $q$ at $(0, 0, \pm d)$. Find the force on the upper one using only the fields on the midplane $z = 0$.

*Setup.* Enclose the upper charge by the midplane and a hemisphere at infinity (fields there fall as $r^{-2}$, stresses as $r^{-4}$, area as $r^2$: no contribution). The outward normal of the upper half-space on the midplane is $-\hat{\mathbf z}$, so by (3.30) with static fields $F_z = -\int T_{zz}\,da$.

*Like charges.* On the plane the $z$-components of the two fields cancel and the $\hat{\mathbf s}$ components add: $\mathbf E = \frac{2qs}{4\pi\varepsilon_0(s^2 + d^2)^{3/2}}\hat{\mathbf s}$. Then $T_{zz} = \varepsilon_0(0 - \tfrac12E^2) = -\tfrac12\varepsilon_0E^2$, a *pressure*. So
$$F_z = \frac{\varepsilon_0}{2}\int_0^\infty\left[\frac{2qs}{4\pi\varepsilon_0(s^2 + d^2)^{3/2}}\right]^22\pi s\,ds = \frac{q^2}{4\pi\varepsilon_0}\int_0^\infty\frac{s^3\,ds}{(s^2 + d^2)^3} = \frac{q^2}{4\pi\varepsilon_0}\cdot\frac{1}{4d^2} = \frac{q^2}{4\pi\varepsilon_0(2d)^2},$$
upward: repulsion, Coulomb's law at separation $2d$. (The integral: put $w = s^2 + d^2$; $\int_{d^2}^\infty\frac{(w - d^2)}{2w^3}dw = \frac12\left[\frac1{d^2} - \frac{d^2}{2d^4}\right] = \frac1{4d^2}$.)

*Unlike charges.* Now the $\hat{\mathbf s}$ components cancel and $\mathbf E = E_z\hat{\mathbf z}$ with $|E_z| = \frac{2qd}{4\pi\varepsilon_0(s^2 + d^2)^{3/2}}$; $T_{zz} = +\tfrac12\varepsilon_0E_z^2$, a *tension*, and $F_z = -\tfrac{\varepsilon_0}{2}\int E_z^2\,da = -q^2/(4\pi\varepsilon_0(2d)^2)$ using $\int_0^\infty s\,ds/(s^2 + d^2)^3 = 1/(4d^4)$: attraction. Like charges push on each other through the pressure of field lines lying in the plane; unlike charges pull through the tension of lines crossing it.

*Numbers.* $q = 1\ \mu$C, $d = 1$ cm: $F = 22.4689$ N. The project builds the full tensor $T_{ij}$ at every point of a polar grid on the plane (with $s = d\tan\vartheta$ mapping the infinite plane to a finite interval), contracts it with $d\mathbf a$, and gets $\pm22.4689$ N for the $z$-component and zero for $x$ and $y$.

## 5. Why it matters and where it connects

- **EM-7.** (3.24) becomes plane waves, polarization (the Jones vector, a qubit in QC-1), reflection with (3.23) at every interface (the project's $r = -1/3$, $t = 2/3$ at $n = 2$ previews the Fresnel equations), dispersion from Drude-like oscillators (§3.1's model with a restoring force), and waveguides. Radiation pressure is (3.31) applied to an absorbing wall.
- **EM-8.** (3.37) is solved with retarded potentials; dipole radiation's energy flux is $\mathbf S$; $\mathbf E$ and $\mathbf B$ become one antisymmetric tensor $F_{\mu\nu}$, $u$, $\mathbf S$, $\mathbf g$ and $T_{ij}$ become one symmetric $4\times4$ stress–energy tensor, and (3.38) is its Lagrangian.
- **QM-1, QM-4, QC.** Light carries energy and momentum in the ratio $1/c$ (§3.11), the classical side of photons. Gauge invariance is the local phase freedom of $\psi$ (§3.14); flux quantization in superconducting loops underlies SQUIDs and flux-tunable qubits.
- **CM-4.** Poynting's theorem and the stress tensor are Noether's theorem for time and space translations of (3.38) (CM-4 §3.10). **LA-12.** $T_{ij}$ is a symmetric $(0,2)$ tensor whose contraction with $d\mathbf a$ is a force, Cauchy's traction $\mathbf t = \sigma\mathbf n$ (LA-12 §3.18).
- **Mechanical and electrical engineering.** Generators (motional EMF), transformers and wireless chargers (mutual inductance), induction motors and eddy-current brakes (Lenz), inductive position sensors (LVDTs), and every numerical electromagnetics code you will meet (FDTD is the project; finite elements extend EM-3's methods).

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: Drude numbers, a loop in an ac field, solenoid energy, a tiny displacement current, gauge choices, stress values. B1–B8 are the module: B1 (Faraday's disk and the flux rule), B3 (induced $\mathbf E$ and $-\partial_t\mathbf A$), B4 (reciprocity and $M^2\le L_1L_2$), B5 (the coaxial cable: $L$ two ways and Poynting), B6 (Maxwell's correction and the boundary conditions) and B8 (`project.py`: Courant and numerical dispersion) are essential; B2 (the falling loop) and B7 (stress on capacitor plates) are drills in energy and momentum bookkeeping. C1–C3: the betatron condition, a lossy FDTD medium, and Maxwell–Ampère from the Lagrangian.

## 7. Coding project

`project.py` runs in about 12 seconds and writes figures to `figures/`. It builds: (1) Drude numbers and the energy balance of the sliding and falling bars (ODE integration); (2) Faraday's law on an off-centre loop, and a SymPy check that (3.11) has the right curl and zero divergence; (3) **Neumann's formula numerically** (`neumann`), checked against the exact coaxial-loop formula, and the two independent computations of Example 4.3(b); magnetic energy of the toroid; (4) the capacitor of Example 4.4; (5) the Poynting field of a resistor on a grid; (6) SymPy proofs that (3.28) and $f_i = \partial_jT_{ij} - \varepsilon_0\mu_0\partial_tS_i$ hold for fields from random polynomial potentials, and **the midplane stress integral** of Example 4.6; (7) **gauge invariance on a grid** and the Lorenz-gauge dipole in SymPy; (8) **a 1D FDTD solver** (`yee1d`); (9) **a 2D TE FDTD** point source (`figures/fdtd_2d.png`); Your turn (lossy medium, Faraday's disk, two coils, with reference solutions); (10) every number in `solutions.md`.

**The Yee scheme and the Courant condition.** For $\mathbf E = E_y(x, t)\hat{\mathbf y}$, $\mathbf B = B_z(x, t)\hat{\mathbf z}$, Maxwell's curl equations are $\partial_tB_z = -\partial_xE_y$ and $\varepsilon\partial_tE_y = -\mu_0^{-1}\partial_xB_z$. Yee (1966) stores $E$ at grid points $x_i$ and integer times $t_n$, $B$ at half points $x_{i+1/2}$ and half times $t_{n+1/2}$, so every derivative is a centred difference (the derivation is in the project's markdown cell). Why the time step is limited:

1. Insert $E^n_i = E_0e^{i(kx_i - \omega t_n)}$, $B^{n+1/2}_{i+1/2} = B_0e^{i(kx_{i+1/2} - \omega t_{n+1/2})}$ into the $B$ update $B^{n+1/2}_{i+1/2} - B^{n-1/2}_{i+1/2} = -\frac{\Delta t}{\Delta x}(E^n_{i+1} - E^n_i)$. Each difference of exponentials gives a sine: $B_0\sin\frac{\omega\Delta t}{2} = \frac{\Delta t}{\Delta x}E_0\sin\frac{k\Delta x}{2}$.
2. The $E$ update in vacuum gives, the same way, $E_0\sin\frac{\omega\Delta t}{2} = \frac{c^2\Delta t}{\Delta x}B_0\sin\frac{k\Delta x}{2}$.
3. Multiply: $\sin\frac{\omega\Delta t}{2} = S\sin\frac{k\Delta x}{2}$, with the **Courant number** $S = c\Delta t/\Delta x$.
4. If $S\le1$, every $k$ gives a real $\omega$: waves neither grow nor decay. If $S > 1$, the shortest wave on the grid ($k\Delta x = \pi$) needs $\sin(\omega\Delta t/2) = S > 1$, so $\omega$ is complex and that mode grows by a factor of order $e^{|\text{Im}\,\omega|\Delta t}$ every step. ∎

At $S = 1$ the relation is $\omega = ck$ exactly in 1D ("magic time step"); at $S = 0.5$ waves are slightly slow, and the project measures the pulse speed $2.3\times10^{-4}$ below $c$. At $S = 1.02$, round-off in the shortest mode grows to $6\times10^{87}$ in 600 steps. In 2D the limit is $S\le1/\sqrt2$. The scheme also conserves a discrete energy $\tfrac12\sum[\varepsilon(E^n)^2 + B^{n-1/2}B^{n+1/2}/\mu_0]\Delta x$ exactly (to $10^{-15}$ in the run), which is positive definite only when $S\le1$: the same condition seen through Poynting's theorem.

The 1D runs: a Gaussian pulse moves at $c$; a perfect conductor ($E = 0$ at the wall) reflects it inverted; at an interface to $n = 2$ the reflected and transmitted amplitudes are $-0.3331$ and $0.6667$ (theory $-1/3$, $2/3$) and the energy fractions $0.1110$ and $0.8890$ ($1/9$, $8/9$); with absorbing (Mur) ends, the time integral of $S_x$ through a plane near the end equals the pulse's energy to $10^{-5}$ and the domain empties. In 2D a point source radiates circular waves whose amplitude falls as $1/\sqrt r$.

**Build it yourself:** before reading `yee1d`, write the two-line leapfrog for a vacuum box with $E = 0$ at both ends, launch $E = e^{-(x - x_0)^2/w^2}$ with $B = E/c$ (shifted by half a time step), and check that the pulse moves right at $c$ and comes back inverted. Then add the dielectric.

## 8. Where you'll struggle

- **The flux rule has real exceptions.** Faraday's disk has an EMF with no flux change; Feynman's rocking plates have a flux change with no EMF. Read Feynman *Lectures* II §17-2; it is required. When in doubt, compute $\oint(\mathbf E + \mathbf v\times\mathbf B)\cdot d\mathbf l$ with $\mathbf v$ the velocity of the actual conducting material; that never fails.
- **Signs in Faraday's law.** The minus sign only means something once the loop's direction and the surface's normal are tied by the right-hand rule. Fix the orientation first, compute $\Phi$ with that normal, then read the direction of $\mathcal E$ from the sign. Check every answer against Lenz.
- **Poynting's theorem: every student loses a sign.** The traps are (P4)'s minus sign and the direction of $\oint\mathbf S\cdot d\mathbf a$ (outward). Do the derivation twice, the second time without looking, and test the result on Example 4.5: if your $\mathbf S$ points out of the resistor, you have lost one.
- **The stress tensor's index notation.** This is the first time a tensor does real work in E&M, and the derivation of §3.11 looks like a wall of indices. The cure is LA-12: once you read $T_{ij}\,da_j$ as "matrix times area vector gives force" and the bracket identity of step 5 as the product rule, it is a page of routine. Then use the tension/pressure picture to get signs without computing.
- **The displacement current is not a current.** No charge moves in a vacuum capacitor gap. The name is Maxwell's (his ether "displaced"); the term is a rate of change of field that sources $\mathbf B$ exactly as a current would. The polarization current $\partial_t\mathbf P$ inside a dielectric *is* moving charge; $\varepsilon_0\partial_t\mathbf E$ is not.
- **Lorenz, not Lorentz.** The gauge is Ludvig Lorenz's (1867); the force and transformation are Hendrik Lorentz's. Most books, Jackson's early editions included, have printed it wrong.
- **FDTD stability.** Change $\Delta x$ without changing $\Delta t$ and the code explodes; the bug is the Courant condition, not your update equations. Derive $\sin(\omega\Delta t/2) = S\sin(k\Delta x/2)$ once (§7) and always set $\Delta t$ from $S$: `dt = S * dx / c`. In a dielectric the local speed is lower, so the vacuum $S$ is the binding one; in 2D and 3D the limits are $1/\sqrt2$ and $1/\sqrt3$.
- **Which $\mathbf E$ does a voltmeter read?** With changing flux around, $\int\mathbf E\cdot d\mathbf l$ between two points depends on the path, so "the voltage across the inductor" depends on where the leads go. Kirchhoff's voltage law survives in circuit theory only because inductors' fields are confined to their own interior (Exercise A2(c)).

## 9. Check yourself

1. Derive $\sigma = nq^2\tau/m$ in two lines. Why does a magnetic field not contribute to Joule heating?
2. Prove the flux rule for a moving loop in a static $\mathbf B$. Where is $\nabla\cdot\mathbf B = 0$ used?
3. Give the EMF of Faraday's disk and say why the flux rule does not obviously apply.
4. Turn the integral form of Faraday's law into $\nabla\times\mathbf E = -\partial_t\mathbf B$. What source does $-\partial_t\mathbf B$ play the role of?
5. Find $E(s)$ inside and outside a solenoid with ramping current. Why is $\mathbf E\ne0$ where $\mathbf B = 0$?
6. Derive Neumann's formula from $\Phi = \oint\mathbf A\cdot d\mathbf l$. Why does it prove $M_{12} = M_{21}$?
7. Get from $\tfrac12LI^2$ to $\int B^2/2\mu_0$. Which product rule, and which surface term vanishes?
8. Show that Ampère's law contradicts charge conservation and derive the displacement current.
9. Derive the polarization current and write Maxwell's equations in matter.
10. Derive the vacuum wave equation and the value of $c$.
11. Derive Poynting's theorem. Which way does energy flow into a resistor?
12. Write $T_{ij}$ and say what $T_{ij}\,da_j$ means. Why do like charges repel in the tension/pressure picture? Write the gauge transformation, prove $\mathbf E$ and $\mathbf B$ are unchanged, and get $\Box V = -\rho/\varepsilon_0$ from the Lorenz condition.

## 10. Pair with

- **Griffiths, *Introduction to Electrodynamics*, ch. 7–8.** The same order and the same derivations (§7.1.3's flux-rule proof, §7.3's displacement current, §8.1–8.2's Poynting theorem and stress tensor). Griffiths has more worked problems and better pictures of the solenoid and capacitor.
- **Purcell and Morin, *Electricity and Magnetism*, ch. 7 and 9.** Induction and Maxwell's equations with more physical insight per page, and the electromagnetic wave built from a sheet of current switched on at $t = 0$, which makes "the curls sustain each other" concrete.
- **Feynman, *Lectures on Physics* Vol. II, ch. 16–18 and 27.** Ch. 17 §2 (exceptions to the flux rule) is required; ch. 18 builds the equations; ch. 27 is the best discussion of where field energy and momentum go, including the disk paradox.
- **MIT 8.07 lecture notes, lectures 15–19.** A more compact, more mathematical pass through the conservation laws, the stress tensor and the gauge choices; read after this module as a second view.
