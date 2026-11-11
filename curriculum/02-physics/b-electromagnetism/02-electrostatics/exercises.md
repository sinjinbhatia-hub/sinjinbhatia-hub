# EM-2 Exercises

Conventions as in the notes: SI units, $k = 1/4\pi\varepsilon_0$, $V(\infty) = 0$ for localized charge, and at a charged surface $\hat{\mathbf n}$ points from side 1 (below) to side 2 (above).

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Superposition.) A charge $q_1 = +2$ nC sits at the origin and $q_2 = -1$ nC at $(0.3\text{ m}, 0, 0)$.
(a) Find $\mathbf E$ and $V$ at $P = (0, 0.4\text{ m}, 0)$ using (3.2).
(b) Find the point on the $x$-axis where $\mathbf E = 0$. First explain why it cannot lie between the charges or left of $q_1$.

**A2.** (Flux through cube faces.)
(a) A point charge $q$ sits at the centre of a cube. What is the flux of $\mathbf E$ through one face?
(b) Now $q$ sits at one corner of a cube. What is the flux through each of the six faces? (Hint: how many identical cubes share that corner?)
(c) In (a), why can't you use Gauss's law to find $\mathbf E$ at the centre of a face?

**A3.** (A charged slab.) Charge fills the slab $|z| < t$ with uniform density $\rho$; there is no charge elsewhere.
(a) Use symmetry and a pillbox centred on $z = 0$ to find $\mathbf E(z)$ everywhere.
(b) Find $V(z)$ with the reference $V(0) = 0$. Why can't you use $V(\infty) = 0$?
(c) Check that $E_z$ is continuous at $z = \pm t$ and say which equation of §3.7 requires it.

**A4.** (Boundary-condition jump.) The plane $z = 0$ carries $\sigma = 5.0\ \mu$C/m². Just below it the field is $\mathbf E_1 = (3.0,\ 0,\ 2.0)\times10^5$ V/m. Find the field $\mathbf E_2$ just above it, and the jump in $\partial V/\partial z$ across the plane.

**A5.** (Capacitor network.) A $2\ \mu$F and a $3\ \mu$F capacitor are in series; that pair is in parallel with a $4\ \mu$F capacitor; the combination is connected to a 12 V battery. Find the equivalent capacitance, the total charge, the charge on and voltage across each capacitor, and the stored energy. Check that the three capacitors' energies add to the total.

**A6.** (Energy of four charges.) Four charges of $+1$ nC sit at the corners of a square of side $a = 10$ cm.
(a) Compute the work needed to assemble them from infinity, bringing them in one at a time.
(b) Compute $\tfrac12\sum_iq_iV(\mathbf r_i)$, with $V(\mathbf r_i)$ the potential at corner $i$ due to the other three, and confirm it agrees with (a).
(c) Explain in one sentence why $\sum_iq_iV(\mathbf r_i)$ without the $\tfrac12$ is wrong.

## B. Core (do all)

**B1.** (Gauss's law for gravity; the tunnel through the Earth.) Newton's law gives the gravitational field of a point mass $M$ at the origin as $\mathbf g = -GM\hat{\mathbf r}/r^2$.
(a) Repeat the derivation of §3.3 for $\mathbf g$ to get $\oint\mathbf g\cdot d\mathbf a = -4\pi GM_{\text{enc}}$ and $\nabla\cdot\mathbf g = -4\pi G\rho_m$.
(b) Prove the shell theorem: a uniform spherical shell of mass attracts an outside particle as if all its mass sat at the centre, and exerts no force on a particle inside.
(c) Model the Earth as a uniform ball ($GM = 3.986\times10^{14}$ m³/s², $R = 6371$ km). Find $\mathbf g$ inside. Drill a straight frictionless tunnel through the centre and drop a stone in. Show the motion is simple harmonic, and find the period, the one-way travel time, and the speed at the centre.
(d) Now drill the tunnel along a chord that misses the centre. Show the motion along the chord has the same period. Compare the period with that of a satellite in a circular orbit skimming the surface.

**B2.** (Capacitor limits.)
(a) In the spherical capacitor (3.18b), set $b = a + d$ and let $d \ll a$. Show $C \to \varepsilon_0A/d$ with $A = 4\pi a^2$. Do the same for the cylindrical capacitor (3.18c) with $A = 2\pi aL$. Why must both reduce to the parallel-plate formula?
(b) Let $b\to\infty$ in (3.18b) to get the capacitance of an isolated sphere. Evaluate it for the Earth ($R = 6371$ km).
(c) In fair weather the field just above the ground is about 100 V/m, pointing down. What is the surface charge density on the ground, and what total charge would the Earth carry if this held everywhere?
(d) (Code.) In `project.py` ("Your turn" 2), solve Laplace's equation with `solve_dirichlet` for a coaxial cross-section, $a = 0.25$, $b = 1$, and compute the capacitance per length $C'$ from `edge_energy` at $h = 0.02, 0.01, 0.005$. Compare with $2\pi\varepsilon_0/\ln(b/a)$. Which way is the error, why, and how does it scale with $h$?

**B3.** (The shell potential by direct integration.) A spherical shell of radius $R$ carries uniform surface charge $\sigma$, total $Q$.
(a) Compute $V(r) = \dfrac{1}{4\pi\varepsilon_0}\displaystyle\int\frac{\sigma\,da'}{\mathscr r}$ at distance $r$ from the centre, for $r < R$ and $r > R$, without using Gauss's law. Put the field point on the $z$-axis and use the law of cosines for $\mathscr r$.
(b) Identify the step where "inside" and "outside" give different answers.
(c) Use (a) to compute $V(0)$ for the uniformly charged solid ball by adding up shells, and check against Example 4.1.

**B4.** (The dipole field and its field lines.)
(a) Derive (4.3) directly in Cartesian form: take $-\nabla$ of $V = \mathbf p\cdot\mathbf r/(4\pi\varepsilon_0r^3)$ using $\nabla(\mathbf p\cdot\mathbf r) = \mathbf p$ and $\nabla(r^{-3}) = -3r^{-5}\mathbf r$.
(b) Using (4.2), show that the field lines of an ideal dipole are the curves $r = r_0\sin^2\theta$.
(c) (Code.) Project Part 1 traces a field line of the physical dipole ($\pm1$ nC, separation $d$) from $(1\text{ m}, 0, 0)$ in to $r = 0.1$ m and reports $\max|r/\sin^2\theta - 1|$ for $d = 1$ and $2$ mm. Explain its size and its scaling with $d$.

**B5.** (Energies do not superpose.)
(a) Charge distributions 1 and 2 produce fields $\mathbf E_1, \mathbf E_2$ and have energies $W_1, W_2$ by (3.16). Show the energy of both together is $W = W_1 + W_2 + \varepsilon_0\int\mathbf E_1\cdot\mathbf E_2\,d\tau$.
(b) Show the interaction term equals $\int\rho_1V_2\,d\tau$. For two point charges, conclude that it is the pair energy $kq_1q_2/\mathscr r_{12}$ of (3.14), and say what $W_1$ and $W_2$ are.
(c) A thin shell of radius $a$ carries $+q$, and a concentric thin shell of radius $b > a$ carries $-q$. Compute $W_1$, $W_2$, and the interaction term separately, add them, and check the total against $(\varepsilon_0/2)\int E^2\,d\tau$ computed directly and against $Q^2/2C$ for the spherical capacitor.

**B6.** (The pressure on a charged conductor.) A conducting sphere of radius $R$ carries charge $Q$.
(a) Use (3.17) to find the outward force per unit area on its surface.
(b) Check (a) by virtual work: compute $dW/dR$ at fixed $Q$ and compare with the total outward force.
(c) Find the net force that the southern hemisphere exerts on the northern one. Evaluate for $Q = 1\ \mu$C, $R = 10$ cm.

**B7.** (A conductor with two cavities.) Two spherical cavities, of radii $a$ and $b$, are hollowed out of a neutral conducting sphere of radius $R$. A point charge $q_a$ sits at the centre of the first cavity and $q_b$ at the centre of the second.
(a) Find the surface charge densities on the two cavity walls and on the outer surface.
(b) Find $\mathbf E$ outside the conductor and inside each cavity.
(c) What is the force on $q_a$ and on $q_b$?
(d) A third charge $q_c$ is brought up from outside. Which answers to (a)–(c) change?
(e) Instead, move $q_a$ off-centre in its cavity, with $q_c$ absent. Which answers change?

**B8.** (Ring and disk on the axis; code.)
(a) A ring of radius $R$ carries charge $Q$ uniformly. Find $\mathbf E$ on its axis at height $z$. Where on the axis is $|\mathbf E|$ largest?
(b) A disk of radius $R$ carries uniform $\sigma$. Find $\mathbf E$ on its axis for $z > 0$ by adding rings. Check the limits $z \to 0^+$ and $z \gg R$.
(c) Why does Gauss's law not give you (b)?
(d) In `project.py` ("Your turn" 1), build the disk ($Q = 1$ nC, $R = 10$ cm) from point charges on a polar grid, sum their fields with `field_potential`, and compare with (b) at $z = 0.01, 0.05, 0.1, 1$ m.

## C. Challenge (do at least one)

**C1.** (Self-energy and the classical electron radius.)
(a) Get the energy of the uniformly charged ball a third way: build it shell by shell, bringing charge $dq$ from infinity to the surface of a ball of radius $r$ that already holds $q(r)$.
(b) Suppose the electron's rest energy $m_ec^2$ were entirely the electrostatic energy of its charge. Find its radius if the charge is a uniform ball and if it is a thin shell. Express both in terms of the classical electron radius $r_e = e^2/(4\pi\varepsilon_0m_ec^2)$ and evaluate.
(c) Why does no one believe either answer? Compare $r_e$ with the reduced Compton wavelength $\hbar/m_ec$ and say which physics takes over first as you look at smaller distances.

**C2.** (Relaxation as diffusion; code.)
(a) Discretize the heat equation $\partial_tu = \alpha\nabla^2u$ on a square grid of spacing $h$ with the five-point Laplacian and explicit Euler. Show that at $\alpha\Delta t = h^2/4$ the update is exactly the Jacobi iteration for Laplace's equation, and that this is the largest stable step.
(b) On the unit square with $N$ interior points per side ($h = 1/(N+1)$) and zero boundary values, show that $e_{mn} = \sin(p\pi mh)\sin(q\pi nh)$ is an eigenvector of one Jacobi sweep with eigenvalue $\tfrac12[\cos(p\pi h) + \cos(q\pi h)]$. Conclude that the error's slowest component shrinks by $\cos(\pi/(N+1))$ per sweep, and estimate the sweeps needed to reduce it by $10^6$ for $N = 30$ and $N = 100$.
(c) Project Part 8 runs Jacobi on the unit square with the top side at 1 V and the others grounded, $N = 30$, from $V = 0$, until $\max|V - V_{\text{exact}}| < 10^{-6}$ V. It measures fewer sweeps than (b) predicts. Explain the difference quantitatively.

**C3.** (Thomson's theorem and the Dirichlet principle.)
(a) Several conductors carry fixed total charges $Q_1, Q_2, \dots$. Prove that among all ways of placing those charges on the conductors (any distribution, surface or volume, as long as each conductor's total is right), the equilibrium distribution, the one that makes each conductor an equipotential, has the least energy (3.16).
(b) Now fix the potentials of the conductors instead, and prove that among all functions $V$ taking those boundary values, the solution of Laplace's equation minimizes $\int|\nabla V|^2\,d\tau$ (F6 Exercise C3 derived the Euler–Lagrange equation; here prove the minimum directly).
(c) Use (a) to explain in one sentence why charge placed inside a conductor moves to the surface.
