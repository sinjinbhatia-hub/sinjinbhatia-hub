# EM-8 Exercises

Conventions as in the notes: SI, Lorenz gauge, $t_r = t - \mathscr r/c$, metric $(-,+,+,+)$, $x^\mu = (ct, \mathbf r)$, $F^{0i} = E_i/c$, $F^{ij} = \varepsilon_{ijk}B_k$. Numbers: $c = 2.998\times10^8$ m/s, $\varepsilon_0 = 8.854\times10^{-12}$ F/m, $\mu_0 = 4\pi\times10^{-7}$ H/m, $e = 1.602\times10^{-19}$ C, $m_e = 9.109\times10^{-31}$ kg, $m_pc^2 = 938.27$ MeV, $a_0 = 5.292\times10^{-11}$ m.

## A. Warm-up (do all; 5–10 minutes each)

**A1. Muons.** Muons made 10 km up have rest-frame mean lifetime $\tau = 2.197\ \mu$s and move at $\beta = 0.995$.
(a) Find $\gamma$, the mean lifetime in the Earth frame, and the mean distance travelled. Compare with $\beta c\tau$.
(b) In the muon's frame, how thick is the 10 km of atmosphere? Explain why both frames agree that many muons reach the ground.

**A2. A 1 GeV proton.** A proton has kinetic energy 1.000 GeV. Find $\gamma$, $\beta$ and $pc$ (in MeV). Check (3.25) and $pc = \gamma\beta mc^2$.

**A3. Blue sky.** By what factor does the power scattered by an air molecule at 450 nm exceed that at 700 nm? State the two physical facts your answer uses.

**A4. Larmor for a slow electron.** An electron moves at $v = 1.0\times10^6$ m/s perpendicular to a uniform $B = 1.0$ T.
(a) Find its acceleration and the Larmor power.
(b) Show that its kinetic energy decays as $e^{-t/T}$ and find $T$. Is radiation loss important for this electron on laboratory time scales?

**A5. A moving observer in a magnet.** In the lab, $\mathbf E = 0$ and $\mathbf B = 1.0\,\hat{\mathbf z}$ T. Find $\bar{\mathbf E}$ and $\bar{\mathbf B}$ in a frame moving at $1.0\times10^3\,\hat{\mathbf x}$ m/s. Check both invariants. What has this to do with the motional EMF of EM-6?

**A6. Index gymnastics.** (a) Write $F_{\mu\nu} = \eta_{\mu\alpha}\eta_{\nu\beta}F^{\alpha\beta}$ as a matrix. (b) Compute $F^\mu{}_\nu$ and show it is not antisymmetric. (c) Verify $F_{\mu\nu}F^{\mu\nu} = 2(B^2 - E^2/c^2)$ component by component.

## B. Core (do all)

**B1. A switched-on current.** An infinite straight wire on the $z$ axis carries $I(t) = 0$ for $t < 0$ and $I(t) = I_0$ for $t\ge0$. The wire stays neutral.
(a) Using (3.2), show that at distance $s$ from the wire, $\mathbf A = 0$ for $t < s/c$ and, for $t > s/c$,
$$\mathbf A(s, t) = \frac{\mu_0I_0}{2\pi}\ln\Big(\frac{ct + \sqrt{c^2t^2 - s^2}}{s}\Big)\hat{\mathbf z}.$$
(b) Find $\mathbf E$ and $\mathbf B$. (c) Show that as $t\to\infty$ they become the static fields of EM-5. Why is there an electric field at all?

**B2. Retarded-time derivatives (the error magnet).** For a point charge on $\mathbf w(t)$, $t_r(\mathbf r, t)$ is defined by (3.6). With $\mathbf u = c\hat{\boldsymbol{\mathscr r}} - \mathbf v$:
(a) Differentiate $c(t - t_r) = |\mathbf r - \mathbf w(t_r)|$ with respect to $t$ and show $\partial t_r/\partial t = \mathscr rc/(\boldsymbol{\mathscr r}\cdot\mathbf u)$.
(b) Take the gradient and show $\nabla t_r = -\boldsymbol{\mathscr r}/(\boldsymbol{\mathscr r}\cdot\mathbf u)$.
(c) Explain why these reduce to $1$ and $-\hat{\boldsymbol{\mathscr r}}/c$ for a charge at rest, and why (3.3) was nevertheless right in Theorem 3.1.
(d) Use (a) to show that $\partial V/\partial t$ for the Liénard–Wiechert potential of a charge in uniform motion is not equal to the naive $\partial/\partial t$ of $q/(4\pi\varepsilon_0\mathscr r)$ with $\mathscr r$ treated as fixed. (One line of reasoning is enough.)

**B3. The pancake in numbers.** A charge moves at $\beta = 0.9$. At a fixed distance $R$ from its present position, find the field strength relative to Coulomb's (a) directly ahead, (b) broadside. (c) Find their ratio. (d) Show that $\mathbf B = \mathbf v\times\mathbf E/c^2$ vanishes ahead and is largest broadside.

**B4. Magnetic versus electric dipole.** (a) The radiation-zone fields of an oscillating magnetic dipole $\mathbf m = m_0\cos\omega t\,\hat{\mathbf z}$ are $\mathbf E = \frac{\mu_0m_0\omega^2}{4\pi c}\frac{\sin\theta}{r}\cos\omega(t - r/c)\hat{\boldsymbol\phi}$, $\mathbf B = -\frac{\mu_0m_0\omega^2}{4\pi c^2}\frac{\sin\theta}{r}\cos\omega(t - r/c)\hat{\boldsymbol\theta}$. Compute $\langle\mathbf S\rangle$ and derive (3.17).
(b) A loop of radius $b$ and a straight dipole of length $d = b$ carry the same current amplitude $I_0$ at the same $\omega$. Show $P_{\rm mag}/P_{\rm el} = (\pi b\omega/c)^2$ and evaluate for $b = 1$ cm at 100 MHz.

**B5. The classical atom, carefully.** (a) Derive $r(t)^3 = a_0^3 - 4r_e^2ct$ from Larmor, stating where you assume the orbit stays circular. (b) Evaluate $t_f$. (c) Show that the number of revolutions before collapse is $2t_f/T_0$, with $T_0$ the initial period, and evaluate it. (d) What does the answer to (c) say about the "nearly circular" assumption?

**B6. Purcell from the electrons' frame.** Use the wire of §3.15, but let the test charge $q$ be *at rest* in the lab. In the lab there is no force (neutral wire, $q$ at rest). Go to the electrons' rest frame $\tilde S$.
(a) Find the line charge density of the wire in $\tilde S$ and the electric force on $q$.
(b) In $\tilde S$, $q$ moves. Find the current, the magnetic field and the magnetic force on $q$.
(c) Show the two forces cancel exactly, to all orders in $v/c$.

**B7. Crossed fields.** In the lab, $\mathbf E = E\hat{\mathbf y}$ and $\mathbf B = B\hat{\mathbf z}$ with $E < cB$.
(a) Show that in a frame moving at $\mathbf v = \mathbf E\times\mathbf B/B^2$ the electric field vanishes, and find $\bar{\mathbf B}$.
(b) Use this to explain the $\mathbf E\times\mathbf B$ drift of a charged particle in crossed fields (EM-5) with no calculation of trajectories.
(c) For $E = 1.0\times10^3$ V/m and $B = 1.0$ T, find the drift velocity. Why is there no such frame when $E > cB$?

**B8. Covariant Maxwell and continuity.** (a) Expand $\partial_\nu F^{\mu\nu} = \mu_0J^\mu$ for $\mu = 2$ and identify the equation. (b) Expand $\partial_\nu G^{\mu\nu} = 0$ for $\mu = 3$. (c) Prove $\partial_\mu\partial_\nu F^{\mu\nu} = 0$ for any antisymmetric $F$, and state what it implies. (d) Show that $\partial_\lambda F_{\mu\nu} + \partial_\mu F_{\nu\lambda} + \partial_\nu F_{\lambda\mu} = 0$ holds identically when $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$, and that the case $(\lambda, \mu, \nu) = (1, 2, 3)$ is $\nabla\cdot\mathbf B = 0$.

**B9. (Code) A circulating charge.** In `project.py`, add a trajectory $\mathbf w = R(\cos\omega t, \sin\omega t, 0)$ with $R\omega/c = 0.005$ and use the Part 3 machinery to compute the period-averaged power and $\langle S_r\rangle(\theta)$ on a sphere of radius $5\lambda$.
(a) Predict the total power from Larmor, and its ratio to that of a single charge oscillating linearly with amplitude $R$.
(b) Predict the angular pattern by writing the motion as two perpendicular dipoles in quadrature. Check both numerically.

## C. Challenge (do at least one)

**C1. Abraham–Lorentz.** (a) For an electron, evaluate $\tau = \mu_0e^2/(6\pi m_ec)$ and show $c\tau = \tfrac23r_e$. (b) A force $F_0$ is switched on at $t = 0$ and stays on. Solve $ma = F + m\tau\dot a$ on $t < 0$ and $t > 0$, require $a$ continuous at $t = 0$ (why continuous?) and no runaway as $t\to\infty$, and find $a(t)$ everywhere. (c) Explain in two sentences why this is unphysical and why it does not matter for any experiment.

**C2. Beaming by aberration.** A particle moving at $\beta\hat{\mathbf x}$ emits light in its rest frame at angle $\bar\theta$ from $\hat{\mathbf x}$. (a) Use the Lorentz transformation of the photon's four-momentum to show the lab angle obeys $\cos\theta = \frac{\cos\bar\theta + \beta}{1 + \beta\cos\bar\theta}$, $\sin\theta = \frac{\sin\bar\theta}{\gamma(1 + \beta\cos\bar\theta)}$. (b) Show that light emitted at $\bar\theta = 90°$ appears at $\sin\theta = 1/\gamma$, so half of all rest-frame emission lands in a forward cone of half-angle $\approx1/\gamma$. Evaluate for $\gamma = 10$ and $1000$. (c) Why does this explain the narrowing seen in the project's synchrotron figure?

**C3. Ampère–Maxwell from the Lagrangian.** Write (3.33) in terms of $V$ and $\mathbf A$. Vary $A_x$: compute $\partial\mathcal L/\partial A_x$, $\partial\mathcal L/\partial(\partial_tA_x)$ and $\partial\mathcal L/\partial(\partial_jA_x)$, and show that the Euler–Lagrange equation is the $x$ component of $\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\partial_t\mathbf E$. (Hint: $\mathbf B\cdot\mathbf B$ depends on $\partial_yA_x$ through $B_z = \partial_xA_y - \partial_yA_x$ and on $\partial_zA_x$ through $B_y$.)
