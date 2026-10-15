# EM-3 Exercises

Conventions as in the notes: SI units, $k = 1/4\pi\varepsilon_0$, $V(\infty) = 0$ for localized charge, $P_\ell(1) = 1$, $x = \cos\theta$ with $\theta$ measured from the symmetry axis, and at a surface $\hat{\mathbf n}$ points from "below" (side 1) to "above" (side 2). Every number asked for below is printed by `project.py`; the solutions say where.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Mean-value property by hand.) Let $f = x^2 + y^2 - 2z^2$ and $g = x^2 + y^2 + z^2$.
(a) Show $\nabla^2 f = 0$. Compute the average of $f$ over the sphere of radius $R_s$ centred at $(x_0, y_0, z_0)$ by writing a point of the sphere as $(x_0 + R_s\sin u\cos v,\ y_0 + R_s\sin u\sin v,\ z_0 + R_s\cos u)$ and using that the three squared direction cosines each average to $\tfrac13$ and the cross terms to $0$. Confirm Theorem 3.1.
(b) Repeat for $g$. Show the average is $g(x_0, y_0, z_0) + R_s^2$ and explain the sign of the extra term from $\nabla^2 g = 6$.

**A2.** (Image numbers.) A charge $q = 3$ nC sits $d = 2$ cm above a grounded conducting plane. Using (3.3)–(3.5), find $\sigma$ directly below the charge and at $s = d$, the induced charge inside the circle $s < 2d$, the force on $q$, and the work needed to pull $q$ off to infinity.

**A3.** (Which method, and the first two Legendre solutions.)
(a) For each configuration name the method of this module you would reach for and write the first equation you would write down: (i) a point charge near a grounded conducting sphere; (ii) a long square pipe with three sides grounded and one at $V_0$; (iii) a long grounded conducting cylinder in a uniform transverse field; (iv) the field of a neutral molecule at distances much larger than the molecule; (v) a point charge inside a grounded spherical cavity.
(b) Write out the four functions $r^\ell P_\ell(\cos\theta)$ and $r^{-(\ell+1)}P_\ell(\cos\theta)$ for $\ell = 0, 1$ in Cartesian coordinates and say which familiar electrostatic potential each one is.

**A4.** (Expanding boundary data.) Use orthogonality (3.16) to expand $\cos^2\theta$ in Legendre polynomials: find $A_0, \dots, A_3$ in $\cos^2\theta = \sum_\ell A_\ell P_\ell(\cos\theta)$. Then write the potential outside a sphere of radius $R$ whose surface is held at $V(R, \theta) = V_0\cos^2\theta$, and identify the monopole and quadrupole moments it implies.

**A5.** (Moments depend on the origin.) Charges $3q$ at $(0, 0, a)$ and $-q$ at $(0, 0, -a)$. Find the monopole moment $Q$ and the dipole moment $\mathbf p$ about the origin, then $\mathbf p$ about the point $(0, 0, a)$, using (3.32). Find the origin for which $\mathbf p = 0$ (the "centre of charge"). Write the leading two terms of $V$ on the positive $z$ axis for $z \gg a$ with the origin at $(0,0,0)$.

**A6.** (Water in a field.) A water molecule has $p = 1.85$ D $= 6.17\times10^{-30}$ C m. Put it in $E = 10^6$ V/m (a strong laboratory field). Find the maximum torque (3.34), the energy difference between the aligned and anti-aligned orientations (3.33), and compare that difference with $k_BT$ at $300$ K. The thermal average of $\cos\theta$ for a dipole in a field is the Langevin function $\langle\cos\theta\rangle = \coth x - 1/x \approx x/3$ with $x = pE/k_BT$; evaluate it. How aligned is water in a strong field?

## B. Core (do all)

**B1.** (Green's reciprocity.) Let $\rho$ and $\rho'$ be two localized charge distributions in the same geometry (conductors included), with potentials $V$ and $V'$.
(a) Prove $\int\rho\,V'\,d\tau = \int\rho'\,V\,d\tau$ from Green's second identity (EM-1 Theorem 3.11, (3.23)) and Poisson's equation.
(b) Use it to show the total charge induced on a grounded infinite plane by a point charge $q$ is exactly $-q$, with no integral over $\sigma$: compare with the system "no point charge, plane held at $V_0$."
(c) Use it to show the total charge induced on a grounded sphere of radius $R$ by a charge $q$ at distance $a$ from its centre is $-qR/a$, confirming (3.6) by a different route. (Comparison system: the sphere alone at potential $V_0$.)

**B2.** (Image in a sphere, completely.) A charge $q$ sits at distance $a > R$ from the centre of a grounded conducting sphere of radius $R$; the image is (3.6).
(a) Verify $V = 0$ on the sphere by the law of cosines (not only at the two poles).
(b) Derive $\sigma(\theta)$ on the sphere, with $\theta$ measured from the direction of $q$, and integrate it over the sphere to show the total is $q' = -qR/a$. Do the integral.
(c) The sphere is now neutral and isolated instead of grounded. Find the extra image needed, the potential of the sphere, and the force on $q$.
(d) With $R = 10$ cm and $q = 1$ nC, compare the grounded and neutral forces for $a = 1.001R$, $1.01R$, $10R$ against the two limiting formulas (plane image near contact; $-2kq^2R^3/a^5$ far away). Which limit does each case reach, and why does the neutral force diverge near contact even though the sphere is neutral?

**B3.** (Box with a $\sinh$.) A rectangular pipe, infinite in $z$, has grounded walls at $y = 0$, $y = a$, $x = 0$ and the wall $x = b$ at potential $V_0$.
(a) Derive $V(x, y)$ by separation of variables (3.10): the $x$ dependence is $\sinh$, not $e^{-k x}$. Explain why.
(b) For $b = a$ evaluate $V$ at the centre with the first $1, 2, 3, 4$ odd terms.
(c) Give a symmetry argument (superposition plus uniqueness) that the exact centre value is $V_0/4$.

**B4.** (Frobenius by hand.) (a) Run the recursion (3.14) for $\ell = 3$ starting from $a_1 = 1$; show it terminates and normalize to $P_3(1) = 1$.
(b) For $\ell = 1$ take the *even* series ($a_0 = 1$, $a_1 = 0$). Compute the coefficients through $x^6$, guess the general coefficient, and show the series sums to $1 - x\operatorname{artanh}x$ for $|x| < 1$. What happens at $x = \pm1$? This is the second solution $Q_1$ of Legendre's equation; explain in one sentence why it is rejected for the sphere and what would change if the region excluded the $z$ axis.

**B5.** (Charged sphere in a field.) A conducting sphere of radius $R$ carries net charge $Q$ and sits in a uniform field $E_0\hat{\mathbf z}$.
(a) Write $V$ outside by superposing Example 4.4 with a point-charge potential, and say which theorem of §3.2 guarantees this is *the* solution.
(b) Find $\sigma(\theta)$ and show $\sigma = 0$ where $\cos\theta_0 = -Q/(12\pi\varepsilon_0R^2E_0)$. What is the threshold $|Q|$ above which $\sigma$ has one sign everywhere?
(c) Numbers: $R = 1$ cm, $E_0 = 10^5$ V/m, $Q = 0.1$ nC. Find the threshold charge, $\theta_0$, and the (radial) field just outside the north and south poles.

**B6.** (A square quadrupole.) Charges $+q$ at $(\pm a, 0, 0)$ and $-q$ at $(0, \pm a, 0)$.
(a) Show $Q = 0$ and $\mathbf p = 0$, and compute the quadrupole tensor (3.31). Check that it is traceless.
(b) Write $V$ to quadrupole order and evaluate it at $(10a, 0, 0)$ in units of $kq/a$. Compute the exact potential there from the four charges and the ratio of the two. Why is the quadrupole term off by about $0.4\%$ and not by $10\%$?

**B7.** (Dipole in the field of a point charge.) A point charge $Q$ sits at the origin; an ideal dipole $\mathbf p$ sits at $\mathbf r$.
(a) From (3.35) derive $\mathbf F = kQ\left[\mathbf p - 3(\mathbf p\cdot\hat{\mathbf r})\hat{\mathbf r}\right]/r^3$, and write the torque and energy from (3.34), (3.33).
(b) Project Part 7 builds a physical dipole $\pm q$ at $\mathbf r \pm \mathbf d/2$ with $qd = p$ fixed and halves $d$ three times. Explain why the relative errors of force, torque, and energy all fall by a factor $4$ per halving and not $2$: which term of the Taylor expansion of $\mathbf E$ about $\mathbf r$ is the first to survive?
(c) Why is the force on a dipole *not* $\nabla(\mathbf p\cdot\mathbf E)$ in general, and when are the two equal?

**B8.** (Code.) In project Part 3 change the slot's boundary data to $V_0(y) = 4V_0\,(y/a)(1 - y/a)$.
(a) Derive the coefficients (3.9): $C_n = 32V_0/(n\pi)^3$ for odd $n$, $0$ for even $n$.
(b) Compare the decay $n^{-3}$ with the step case $n^{-1}$ and explain it with F2 Theorems 3.9 and 3.10 (smoothness of the odd periodic extension).
(c) Compute $V(a/2, a/2)$ from one term and from the full series, then run the SOR solver on the new data and compare. The reference solution is the last "Your turn" cell.

## C. Challenge (do at least one)

**C1.** (Charge between two grounded planes.) A charge $q$ sits at height $z_0$ between grounded planes $z = 0$ and $z = L$.
(a) Construct the infinite set of images (positions and signs) and show $V = 0$ on both planes.
(b) Prove, using Green's reciprocity (B1) with the comparison system "no charge, bottom plane at $0$, top plane at $V_0$" (whose potential is $V' = V_0z/L$), that the induced charges are $Q_{\text{top}} = -qz_0/L$ and $Q_{\text{bottom}} = -q(L - z_0)/L$.
(c) Sum the image series numerically (project "Your turn 1" has a reference) and check (b) for $z_0 = 0.3L$.

**C2.** (The generating function proves everything.) Let $G(x, t) = (1 - 2xt + t^2)^{-1/2} = \sum_{\ell\ge0}c_\ell(x)t^\ell$ for $|t| < 1$.
(a) Show $(1 - 2xt + t^2)\,\partial G/\partial t = (x - t)\,G$. Match powers of $t$ to get $(\ell + 1)c_{\ell+1} = (2\ell+1)x\,c_\ell - \ell\,c_{\ell-1}$, Bonnet's recursion (3.18). With $c_0 = 1$, $c_1 = x$, conclude $c_\ell = P_\ell$.
(b) Show $\dfrac{\partial}{\partial x}\left[(1 - x^2)\dfrac{\partial G}{\partial x}\right] + t\,\dfrac{\partial^2}{\partial t^2}\big(tG\big) = 0$ and deduce that every $c_\ell$ satisfies Legendre's equation (3.13).
(c) Use (a) and (b) to justify the multipole expansion (3.28) of $1/\mathscr r$ for $r > r'$, and say what happens for $r < r'$.

**C3.** (Keesom: the thermally averaged dipole–dipole force.) Two water molecules ($p = 6.17\times10^{-30}$ C m) are a fixed distance $r$ apart and free to rotate at temperature $T$. Let $U_0 = kp^2/r^3$ and write (3.36) as $U = U_0\,g$ with $g = \hat{\mathbf p}_1\cdot\hat{\mathbf p}_2 - 3(\hat{\mathbf p}_1\cdot\hat{\mathbf r})(\hat{\mathbf p}_2\cdot\hat{\mathbf r})$.
(a) Show the unweighted average of $U$ over all orientations of both dipoles is zero.
(b) Boltzmann-weight the orientations and expand $e^{-U/k_BT}\approx 1 - U/k_BT$ to show $\langle U\rangle \approx -\langle U^2\rangle/k_BT = -\tfrac23 U_0^2/k_BT$. (You need $\langle g^2\rangle = \tfrac23$; compute it.)
(c) Evaluate $U_0$ and $\langle U\rangle$ at $r = 0.5, 1, 2$ nm and $T = 300$ K, and compare with the exact weighted average from project Part 7. Why does the approximation fail at $0.5$ nm?
(d) In one paragraph: why is this a $1/r^6$ attraction, and how do the Keesom (this problem), Debye (permanent–induced), and London (induced–induced) pieces of the van der Waals force differ? (SM-6 uses the result.)
