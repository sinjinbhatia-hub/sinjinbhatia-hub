# EM-6 Exercises

Equation numbers refer to [notes.md](notes.md). Numerical answers are reproduced by `project.py` (Part 10 unless another part is named). Do A1–A6 and B1–B8 before moving on; do at least one C.

**A. Warm-up** (do all; 5–10 minutes each)

**A1. Drude numbers.** An aluminium wire of diameter $2$ mm and length $10$ m carries $5$ A. Aluminium has $\sigma = 3.77\times10^7$ S/m and $n = 1.81\times10^{29}$ conduction electrons per m³. Find (a) $R$ and the Joule power; (b) the field $E$ inside the wire; (c) the drift speed; (d) the Drude collision time $\tau$ from (3.2).

**A2. A loop in an ac field.** A $10\times20$ cm rectangular loop of resistance $0.5\ \Omega$ lies in a uniform field perpendicular to its plane, $B(t) = B_0\sin\omega t$ with $B_0 = 50$ mT and $f = 60$ Hz.
(a) Find the peak EMF and peak current.
(b) Take $\mathbf B$ pointing toward you. Which way does the current circulate at $t = 0$? Justify with Lenz's law and with the sign in (3.7).
(c) A uniform resistive ring encircles a long solenoid whose changing flux gives $\oint\mathbf E\cdot d\mathbf l = \mathcal E$ around the ring. A voltmeter is connected to diametrically opposite points $P$, $Q$, its leads running outside the ring either on the left or on the right. Show that the readings are $\mp\mathcal E/2$. Which assumption of circuit theory fails?

**A3. Solenoid energy.** A solenoid of $400$ turns, length $20$ cm, radius $1.5$ cm carries $3$ A. (a) Find $L$ from (3.14) and $W = \tfrac12LI^2$. (b) Find $B$ and $u_B = B^2/2\mu_0$, and check that $u_B$ times the volume equals $W$. (c) What electric field has the same energy density? Compare with the breakdown field of air ($\approx3$ MV/m).

**A4. A tiny displacement current.** A capacitor with circular plates of radius $3$ cm and gap $2$ mm is driven by $V(t) = V_0\sin\omega t$, $V_0 = 50$ V, $f = 1$ kHz. Find the peak displacement current, the peak displacement current density, and the peak $B$ at the rim of the plates. Compare with the Earth's field ($\sim50\ \mu$T) and explain why Ampère got away without Maxwell's term.

**A5. Two gauges for one field.** Take $V = -E_0x$, $\mathbf A = 0$. (a) Find $\mathbf E$ and $\mathbf B$. (b) Apply (3.34) with $\lambda = -E_0xt$ and find $V'$, $\mathbf A'$. Check that $\mathbf E$, $\mathbf B$ are unchanged. (c) Which of the two pairs satisfies the Coulomb condition? The Lorenz condition (3.36)?

**A6. Reading the stress tensor.** (a) Write $T$ for $\mathbf E = E_0\hat{\mathbf x}$ and for $\mathbf B = B_0\hat{\mathbf z}$. Give the eigenvalues and eigenvectors and say which are tensions and which pressures. Evaluate for $E_0 = 1$ MV/m and $B_0 = 1$ T. (b) Rotate: $\mathbf E = E_0(\hat{\mathbf x} + \hat{\mathbf y})/\sqrt2$. Write the $3\times3$ matrix $T_{ij}$ and find its eigenvectors. (c) Show that $T_{ii} = -u$ for any fields.

**B. Core** (do all; these are the module)

**B1. Faraday's disk.** A copper disk of radius $a = 10$ cm spins at $3000$ rpm in a uniform $B = 0.5$ T along its axis; brushes touch the axle and the rim.
(a) Derive the EMF (3.6) and evaluate it. What is the rim speed?
(b) Why does the flux rule (3.5) not apply directly? Rescue it: find the rate at which a radius moving with the material sweeps flux.
(c) Why does the EMF not depend on the path from axle to rim? Run the project's Your turn 2.
(d) A current $I = 10$ A flows radially outward from axle to rim. Find the magnetic torque on the disk and show that the power needed to keep it turning is $\mathcal EI$.

**B2. A falling loop.** A rectangular loop of width $\ell = 5$ cm, mass $m = 5$ g and resistance $R = 0.01\ \Omega$ hangs in a vertical plane with its top edge inside a region of uniform horizontal field $B = 1$ T (perpendicular to the loop) and its bottom edge below the region. It is released from rest.
(a) While the top edge is still in the field, show $m\dot v = mg - (B\ell)^2v/R$ and find the terminal speed.
(b) Solve for $v(t)$. Find the time constant and $v$ at $t = 50$ ms.
(c) At terminal speed find $I$ and show that $mgv_t = I^2R$. Where does the gravitational energy go?

**B3. Inside and outside a ramping solenoid.** A long solenoid with $n = 1000$ turns/m and radius $R = 3$ cm has $dI/dt = 500$ A/s.
(a) Find $|\mathbf E|$ at $s = 1$ cm and $s = 5$ cm, and its direction relative to the current.
(b) An electron is at rest at $s = 1$ cm. Find its initial acceleration.
(c) A coaxial ring of radius $5$ cm and resistance $0.2\ \Omega$ encircles the solenoid. Find the EMF and current. Do they change if the ring is moved off-centre or enlarged? Why?
(d) Show that (3.11) equals $-\partial\mathbf A/\partial t$ with EM-5 (3.26). Use (3.32) to explain why $\mathbf E\ne0$ outside although $\mathbf B = 0$ there.

**B4. Reciprocity and energy of two circuits.**
(a) A long straight wire and a coplanar rectangular loop: the loop's sides of length $a = 20$ cm are parallel to the wire; its other sides have length $b = 10$ cm, and its near side is $s = 5$ cm from the wire. Show $M = \frac{\mu_0a}{2\pi}\ln\frac{s + b}{s}$ and evaluate it. If the wire carries $10\sin(2\pi\cdot60\,t)$ A, find the peak EMF in the loop.
(b) By (3.12) the same $M$ gives the EMF induced "in the wire" when the loop's current changes. What does that mean once the wire is closed by a distant return path?
(c) Ramp $I_1$ from $0$ with $I_2 = 0$, then ramp $I_2$ with $I_1$ fixed. Compute the work against all back-EMFs and show $W = \tfrac12L_1I_1^2 + MI_1I_2 + \tfrac12L_2I_2^2$. Why must the other order give the same $W$?
(d) Use (3.17) to argue $W\ge0$ for every $(I_1, I_2)$ and deduce $M^2\le L_1L_2$. Which fact about $2\times2$ symmetric matrices did you use (LA-8)?

**B5. The coaxial cable.** Inner conductor radius $a = 0.5$ mm, outer conductor (thin) radius $b = 2$ mm, current $I$ in along the inner and back along the outer.
(a) Find $L$ per unit length from the flux between the conductors.
(b) Find it again from $\int B^2/2\mu_0\,d\tau$ per unit length (assume the current flows on the surface of the inner conductor). Evaluate.
(c) The cable joins a battery to a resistor: $V = 100$ V, $I = 2$ A. With $\mathbf E = \frac{V}{s\ln(b/a)}\hat{\mathbf s}$ between the conductors, find $\mathbf S$, its values at $s = a$ and $b$, and its flux through the annulus. Interpret.
(d) If the conductors have a small resistance, which way does $\mathbf S$ tilt near their surfaces, and why?

**B6. Displacement current, polarization current, and a radiating sheet.**
(a) Take the divergence of (3.18) and show it reproduces continuity only together with Gauss's law. Take the divergence of Faraday's law: what does it say about $\nabla\cdot\mathbf B$ at later times?
(b) A parallel-plate capacitor of area $A$ is filled with a linear dielectric $\varepsilon_r$ and charged by current $I$. In the dielectric, find $\partial_t\mathbf D$, the polarization current density $\partial_t\mathbf P$ and the vacuum term $\varepsilon_0\partial_t\mathbf E$. What fraction of the "total current" in the gap is moving bound charge?
(c) A current sheet $\mathbf K = K_0\cos(\omega t)\,\hat{\mathbf x}$ fills the plane $z = 0$ in vacuum. Show that $\mathbf E = E_0\cos(\omega(t - |z|/c))\hat{\mathbf x}$, with $\mathbf B$ chosen so each side's wave moves away from the sheet, solves Maxwell's equations for $z\ne0$; find $E_0$ from (3.23).
(d) For (c), compute the time-averaged power per area the source must supply against $\mathbf E$ (from $-\mathbf E\cdot\mathbf K$) and the power per area carried away by $\mathbf S$ on both sides. Evaluate for $K_0 = 1$ A/m.

**B7. Stress on capacitor plates and solenoid windings.**
(a) A parallel-plate capacitor (area $A$, charge $\pm Q$, fringing neglected). Enclose the top plate by a box whose bottom face is in the gap. Use (3.30) to find the force on the plate. Compare with EM-2's surface pressure (3.17). Evaluate for $Q = 1$ nC, $A = 1$ cm².
(b) Use $T$ to find the pressure on the windings of the solenoid of A3. Is it inward or outward? Check it with virtual work: $F = -\partial W/\partial R$ at fixed flux.

**B8. Extend `project.py`: numerical dispersion and the Courant limit.** Use `yee1d` with $\Delta x = 1$ mm.
(a) Launch $E = e^{-(x - x_0)^2/w^2}\cos(k_0(x - x_0))$ with $w = 40$ mm, $k_0\Delta x = 0.5$ and the matching right-moving $B$. For $S = 0.25, 0.5, 0.9, 1.0$ measure the speed of the centroid of $E^2$ over $600$ mm of travel.
(b) From $\sin(\omega\Delta t/2) = S\sin(k\Delta x/2)$ (§7) derive the numerical group velocity $v_g = c\cos(k\Delta x/2)/\sqrt{1 - S^2\sin^2(k\Delta x/2)}$ and compare with (a). Which $S$ is most accurate, and why is the answer special to 1D?
(c) Run $S = 1.02$ and measure the growth factor per step of $\max|E|$ between steps 400 and 600. Predict it from the dispersion relation at $k\Delta x = \pi$.

**C. Challenge** (do at least one)

**C1. The betatron.** Electrons circle at fixed radius $R$ in a field $B(s, t)$ along the axis (symmetric about it), accelerated by the induced $\mathbf E$ of the growing flux $\Phi(t)$ inside the orbit.
(a) Show that $dp/dt = (e/2\pi R)\,d\Phi/dt$ and that a circular orbit needs $p = eRB(R)$. Deduce the **betatron condition** $B(R) = \tfrac12\langle B\rangle$, with $\langle B\rangle = \Phi/\pi R^2$. Does the argument survive relativity?
(b) For $R = 0.5$ m and a final $B(R) = 0.4$ T, find the final momentum ($pc$ in MeV), the kinetic energy and the flux.
(c) If the flux ramps linearly in $4.0$ ms, find the energy gained per turn and the number of turns. Check against the distance light travels in $4$ ms.

**C2. A lossy medium in FDTD.** (a) Derive $c_a$, $c_b$ in `yee1d` by evaluating $\sigma E$ at $t_{n+1/2}$ as the average of $E^n$ and $E^{n+1}$. (b) From $\varepsilon_0\partial_t^2E + \sigma\partial_tE = \mu_0^{-1}\partial_x^2E$ show that for $\sigma\ll\omega\varepsilon_0$ a wave decays as $e^{-\alpha x}$, $\alpha = \sigma Z_0/2$ ($Z_0 = \mu_0c$), independent of $\omega$. (c) Run Your turn 1 and compare. Why only percent-level agreement? (d) What changes when $\sigma\gg\omega\varepsilon_0$?

**C3. Maxwell–Ampère from the Lagrangian.** (a) With $\mathcal L$ of (3.38) written in terms of $V$, $\mathbf A$ and their derivatives, write the Euler–Lagrange equation for $A_i$ and show that it is $\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\partial_t\mathbf E$. (b) Show that under (3.34), $\mathcal L$ changes by a total derivative plus $-\lambda(\partial_t\rho + \nabla\cdot\mathbf J)$. Conclude that gauge invariance of the action requires charge conservation. (c) Why does $V$'s Euler–Lagrange equation contain no $\partial_t^2V$? What does that say about $V$ as a dynamical variable?
