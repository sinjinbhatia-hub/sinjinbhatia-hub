# EM-5 Exercises

Equation numbers refer to [notes.md](notes.md). Numerical answers are reproduced by `project.py` (Part 9). Do A1–A5 and B1–B8 before moving on; do at least one C.

**A. Warm-up** (do all; 5–10 minutes each)

**A1. Cyclotron numbers.** A 5 MeV $\alpha$ particle ($q = 2e$, $m = 6.645\times10^{-27}$ kg) enters a uniform field $B = 0.8$ T perpendicular to its velocity.
(a) Find its speed, the radius of its circle and its orbital frequency. Check that treating it non-relativistically is fair.
(b) What kinetic energy must a proton have to follow the *same* circle?
(c) In one line, why does the field not change the particle's speed?

**A2. Wires.** (a) Find $B$ at $5$ cm from a long straight wire carrying $30$ A. (b) Find $B$ at the same distance from the midpoint of a $20$ cm straight segment carrying $30$ A, using (3.10). What fraction of the long-wire value is it? (c) Two long parallel wires $2$ cm apart carry $10$ A and $25$ A in *opposite* directions. Find the force per unit length and say whether it attracts or repels.

**A3. A coil on its axis.** A coil of $N = 50$ turns, radius $4$ cm, carries $1.5$ A. (a) Find $B$ at its centre and on the axis $3$ cm from the centre. (b) At what distance along the axis does $B$ fall to half its centre value?

**A4. A toroid.** A toroid of $N = 1200$ turns carries $0.8$ A; its cross-section runs from inner radius $6$ cm to outer radius $10$ cm. Find $B$ at $s = 6, 8, 10$ cm. What is $B$ in the hole and outside the winding? Justify both with an Ampère path.

**A5. Signs and units.** (a) A wire carries current along $+\hat{\mathbf x}$. Give the direction of $\mathbf B$ at $(0,0,1)$ m and at $(0,1,0)$ m. (b) A circular loop in the $xy$-plane carries current *clockwise* as seen from $+z$. Give the directions of $\mathbf m$ and of $\mathbf B$ on the axis. A second loop with counterclockwise current sits coaxially above it: do they attract or repel? (c) Show that $\mu_0I/2\pi s$ is in tesla when $\mu_0$ is in N/A$^2$, that $\oint\mathbf A\cdot d\mathbf l$ is in webers if $\mathbf A$ is in T m, and that $1/\sqrt{\mu_0\varepsilon_0}$ is a speed. Evaluate it.

**B. Core** (do all; these are the module)

**B1. Surface independence.** (a) Prove that for steady currents $\int_{\mathcal S_1}\mathbf J\cdot d\mathbf a = \int_{\mathcal S_2}\mathbf J\cdot d\mathbf a$ whenever $\mathcal S_1$ and $\mathcal S_2$ have the same boundary. (b) A long wire of radius $a$ carries uniform $\mathbf J = J_0\hat{\mathbf z}$. The loop $\mathcal P$ is a circle of radius $b > a$ in the plane $z = 0$, centred on the axis. Compute $I_{\rm enc}$ directly through (i) the flat disk and (ii) the hemisphere of radius $b$ bulging into $z > 0$ (use spherical coordinates and the cap where the wire pierces the sphere). (c) Which step of your proof in (a) fails for a charging capacitor with $\mathcal S_2$ passing between the plates, and what does that tell you about (3.19) there?

**B2. Ampère with a non-uniform $\mathbf J$.** A long cylindrical wire of radius $a$ has $\mathbf J = J_0(s/a)\hat{\mathbf z}$ for $s < a$ and $0$ outside. (a) Find the total current $I$. (b) Use Ampère's law to find $\mathbf B$ for $s < a$ and $s > a$. (c) Verify $(\nabla\times\mathbf B)_z = \mu_0J_z$ inside, with the cylindrical curl. (d) For $a = 1$ cm and $I = 10$ A find $B(a)$ and $B(a/2)$ and compare $B(a/2)$ with the value for uniform $\mathbf J$ of the same $I$. At what radius is half of $I$ enclosed?

**B3. Anti-Helmholtz pair.** Two coaxial rings of radius $R$ sit at $z = \pm d/2$. The ring at $z = +d/2$ carries $+I$ (counterclockwise from $+z$), the other $-I$. (a) Show $B_z$ is odd in $z$ and find $B_z'(0)$ as a function of $d$. (b) Show $B_z'''(0) = 0$ if and only if $d = \sqrt3R$. (c) At $d = \sqrt3R$ show $B_z'(0) = (48\sqrt3/7^{5/2})\,\mu_0I/R^2\approx0.641\,\mu_0I/R^2$. (d) Use $\nabla\cdot\mathbf B = 0$ near the axis to show $B_s\approx-\tfrac12s\,B_z'(0)$. (e) For $R = 5$ cm and $I = 1000$ A (total ampere-turns per ring) find the gradient in T/m, and use `biot_savart` to find how far $B_z$ departs from linear at $z = R/4$ and $R/2$; compare with the same pair at $d = R$.

**B4. Two gauges.** A uniform field $\mathbf B = B_0\hat{\mathbf z}$. (a) Show that $\mathbf A_1 = \tfrac12B_0(-y, x, 0)$ and $\mathbf A_2 = B_0(0, x, 0)$ both have curl $\mathbf B$ and zero divergence. (b) Find $\lambda$ with $\mathbf A_2 = \mathbf A_1 + \nabla\lambda$ and show $\nabla^2\lambda = 0$. What does this say about the Coulomb gauge? (c) For $B_0 = 0.5$ T evaluate $\oint\mathbf A\cdot d\mathbf l$ for both around a $10$ cm square in the $xy$-plane, and compare with the flux. (d) Which gauge would you choose for a problem with rotational symmetry about $z$?

**B5. The sheet's potential.** A sheet $\mathbf K = K\hat{\mathbf x}$ lies in $z = 0$. (a) Starting from (3.20), find $\mathbf A = A_x(z)\hat{\mathbf x}$ with $\nabla\times\mathbf A = \mathbf B$. (b) Check $\nabla\cdot\mathbf A = 0$, $\nabla^2\mathbf A = -\mu_0\mathbf J$, continuity of $\mathbf A$, and the jump $-\mu_0\mathbf K$ in $\partial\mathbf A/\partial z$ ((3.29)). (c) Why can $\mathbf A$ not vanish at infinity? Which result of EM-2 is the analogue? (d) Evaluate $|\mathbf B|$ for $K = 100$ A/m.

**B6. A magnetized cylinder.** A cylinder of radius $5$ mm and length $20$ mm is uniformly magnetized along its axis, $\mathbf M = 10^6\hat{\mathbf z}$ A/m. (a) Find $\mathbf J_b$ and $\mathbf K_b$. (b) Using (3.14) find $\mathbf B$ at the centre and at the centre of an end face. (c) Find $\mathbf H$ at the centre, and just inside and just outside the end face; show that $H^\perp$ jumps by $M$ there. (d) Compare $H$ at the centre with the sphere value $-M/3$ and the long-rod value $0$, and explain the order.

**B7. A permeable wire.** A long straight wire of radius $a = 5$ mm and linear permeability $\mu_r = 200$ carries a uniform free current $I = 20$ A. (a) Find $\mathbf H$, $\mathbf B$, $\mathbf M$ everywhere. (b) Find $\mathbf J_b$ and $\mathbf K_b$ and show that the total bound current through a cross-section vanishes. Why must it? (c) Compare $B$ at $s = 2a$ with the vacuum value.

**B8. Extend `project.py`: force and torque on a small loop.** Write, before reading Part 9, a function that takes a source (`a, b, I` arrays), a loop centre `c`, a normal `n`, a side length and a current, builds the four edges of a small square loop, evaluates `biot_savart` at the edge midpoints, and returns $\mathbf F = I\sum\mathbf{dl}\times\mathbf B$ and $\mathbf N = I\sum(\mathbf r - \mathbf c)\times(\mathbf{dl}\times\mathbf B)$ (from (3.6)). Take a $2$ mm loop with $1$ A, centred at $(0.03, 0, 0.10)$ m with its normal tilted $30^\circ$ from $z$ toward $x$, in the field of the $10$ cm, $2$ A ring. Compare $\mathbf F$ with $\nabla(\mathbf m\cdot\mathbf B)$ (central differences) and $\mathbf N$ with $\mathbf m\times\mathbf B$ ((3.34), (3.35)). Report the relative differences and how they change when you halve the side.

**C. Challenge** (do at least one)

**C1. The delta function in the dipole field.** (a) Using $\int\nabla\times\mathbf A\,d\tau = \oint\hat{\mathbf n}\times\mathbf A\,da$ (§3.13, step 3) and the dipole potential (3.32), show that $\int_{r<R}\mathbf B_{\rm dip}\,d\tau = \tfrac23\mu_0\mathbf m$ for any $R$. (b) Show that the non-delta part of (3.33) integrates to zero over the shell $\epsilon < r < R$. (c) Conclude that the delta term in (3.33) is $\tfrac{2\mu_0}{3}\mathbf m\,\delta^3(\mathbf r)$. (d) Check (a) against the uniformly magnetized sphere of Example 4.5, and say why the same integral of $\mathbf E_{\rm dip}$ has the *opposite* sign relative to $\mathbf p$.

**C2. A saturating core.** Replace the linear iron of Example 4.6 by $B = B_s\tanh(\mu_0\mu_rH/B_s)$ with $B_s = 1.6$ T and $\mu_r = 1000$, same geometry ($l = 0.636$ m, $g = 4$ mm). (a) Write the single equation for the common $B$ in core and gap given $NI$. (b) Solve it numerically (`scipy.optimize.brentq`) for $NI = 1000, 3000, 5000, 10^4$ A and compare with (3.46). (c) How many ampere-turns give $B = 1.0$ T, versus the linear model? (d) What is the largest gap field this model can reach, and why is that the right physical conclusion for a design that needs $2$ T?

**C3. Demagnetizing a Preisach ensemble.** Use the ensemble of Your turn 3. (a) Saturate at $H = +3$, then set $H = 0$: what remanence remains? (b) Instead apply $H_k = 3(-0.9)^k$ for $k = 1, \dots, 80$ and then $H = 0$. Report $M$ after $k = 10, 20, 40, 80$ and plot $M$ against $H$. (c) Repeat with factor $0.5$. (d) Explain, in terms of hysterons, why decaying alternation erases the remanence, and why coarser steps leave more.
