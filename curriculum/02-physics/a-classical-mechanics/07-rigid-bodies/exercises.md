# CM-7 Exercises

Conventions as in the notes: body axes $\hat{\mathbf e}_k$, space axes $\hat{\mathbf x}, \hat{\mathbf y}, \hat{\mathbf z}$; $\mathbf r_{\rm space} = R\,\mathbf r_{\rm body}$; $[\mathbf v]_\times\mathbf u = \mathbf v\times\mathbf u$ (3.3); $I_{ij}$ as in (3.8); Euler angles $(\phi, \theta, \psi)$ are $z$-$x'$-$z''$, $R = R_z(\phi)R_x(\theta)R_z(\psi)$ (3.21); for the heavy top $I_1 = I_2$ and $I_3$ are about the pivot. $g = 9.81$ m/s². Every number the solutions quote is reproduced in Part 9 of `project.py`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Moments by integration.) A uniform solid cylinder has mass $M$, radius $R$, height $h$. Slice it into disks and use the disk results of §3.7 ($\frac12mr^2$ about the axis, $\frac14mr^2$ about a diameter) with the parallel-axis theorem to find $I_3$ about its axis and $I_1 = I_2$ about a perpendicular axis through the CM. For what ratio $h/R$ is the cylinder a spherical top about its CM? Find $I_1$ about a diameter of one end face as well.

**A2.** (The hat map.) Let $\boldsymbol\omega = (1, 2, 3)$ and $\mathbf u = (4, 5, 6)$.
(a) Write $[\boldsymbol\omega]_\times$ and compute $[\boldsymbol\omega]_\times\mathbf u$. Compare with $\boldsymbol\omega\times\mathbf u$.
(b) Show that $\mathbf u^\top[\boldsymbol\omega]_\times\mathbf u = 0$ for every $\mathbf u$. What does this say about the velocity of a point of a rotating body relative to its position vector?
(c) Show $[\boldsymbol\omega]_\times^2 = \boldsymbol\omega\boldsymbol\omega^\top - |\boldsymbol\omega|^2\mathbb 1$.
(d) With $R = R_z(90°)$ and $\mathbf a = \hat{\mathbf x}$, check $R[\mathbf a]_\times R^\top = [R\mathbf a]_\times$ by multiplying out.

**A3.** (A point-mass tensor.) Redo the bent-shaft dumbbell of §2 with the tensor. Two masses $m$ sit at $\pm b\hat{\mathbf n}$ with $\hat{\mathbf n} = (\sin\alpha, 0, \cos\alpha)$.
(a) Show $I = 2mb^2(\mathbb 1 - \hat{\mathbf n}\hat{\mathbf n}^\top)$ and find its principal moments and axes.
(b) Compute $\mathbf L = I\boldsymbol\omega$ for $\boldsymbol\omega = \omega\hat{\mathbf z}$ and check it against §2.
(c) With $m = 0.5$ kg, $b = 0.2$ m, $\alpha = 30°$, $\omega = 10$ rad/s, find $\mathbf L$, the angle between $\mathbf L$ and $\boldsymbol\omega$, and the magnitude of the bearing torque.

**A4.** (One evaluation of Euler's equations.) A body has principal moments $(I_1, I_2, I_3) = (2, 3, 4)$ kg m². At one instant $\boldsymbol\omega = (1, -2, 0.5)$ rad/s and $\boldsymbol\tau = (0, 0, 1)$ N m, both in body components. Find $\dot{\boldsymbol\omega}$ from (3.17). Then check $\dot T = \boldsymbol\omega\cdot\boldsymbol\tau$ and $\frac{d}{dt}|\mathbf L|^2 = 2\mathbf L\cdot\boldsymbol\tau$ with your numbers.

**A5.** (Gimbal lock in words.) (a) Show that $(\phi, \theta, \psi) = (0.7, 0, 0.4)$ and $(0.2, 0, 0.9)$ give the same $R$. (b) A body starts at $R = \mathbb 1$ and rotates about the space $\hat{\mathbf x}$ axis, $R(t) = R_x(\Omega t)$; then, from $R = \mathbb 1$ again, about $\hat{\mathbf y}$, $R(t) = R_y(\Omega t)$. Find the Euler angles for $t > 0$ in each case. (c) Use (b) to explain why no smooth choice of $(\phi, \psi)$ exists near $\theta = 0$, and why an integrator that stores Euler angles fails there while one that stores $R$ does not.

## B. Core (do all; these are the module)

**B1.** (Euler's equations by matrices; the bounds on $|\mathbf L|^2$.)
(a) Write $\mathbf L_{\rm space} = R\,\mathbf L_{\rm body}$ and use $\dot R = R[\boldsymbol\omega_{\rm body}]_\times$ (§3.2) to derive $\dot{\mathbf L}_{\rm body} + \boldsymbol\omega_{\rm body}\times\mathbf L_{\rm body} = \boldsymbol\tau_{\rm body}$ without the basis-vector argument of §3.8.
(b) With $I_1 \le I_2 \le I_3$, prove $2TI_1 \le |\mathbf L|^2 \le 2TI_3$, and say when each equality holds.
(c) Define $D = |\mathbf L|^2/2T$. Explain in terms of $D$ which polhode (loop around $\hat{\mathbf e}_1$, loop around $\hat{\mathbf e}_3$, or separatrix) a torque-free body follows.

**B2.** (The cube about an edge.) A uniform cube of mass $M$, side $a$, occupies $[0, a]^3$.
(a) Find its moment about the edge along the $z$-axis twice: from Example 4.1's tensor, and from the scalar parallel-axis theorem with $I^{\rm cm} = \frac16Ma^2$.
(b) Find the inertia tensor about the midpoint of that edge, $(0, 0, a/2)$, and its principal moments and axes.
(c) The answer to (b) has a degenerate pair, as Example 4.1 did. Prove that the inertia tensor of *any* spherical top about *any* point $O$ is a symmetric top whose axis points from $O$ to the CM.
(d) The cube hangs from a horizontal edge (a frictionless hinge) and swings as a physical pendulum. Find the small-oscillation period for $a = 0.1$ m.

**B3.** (Prolate bodies and the outside cone.) A free symmetric top has $I_3 < I_1$ (prolate: a rugby ball, a pencil).
(a) Let $\theta_\omega$ and $\theta_L$ be the angles of $\boldsymbol\omega$ and $\mathbf L$ from $\hat{\mathbf e}_3$. Show $\tan\theta_L = (I_1/I_3)\tan\theta_\omega$, so $\boldsymbol\omega$ lies between $\hat{\mathbf e}_3$ and $\mathbf L$. Conclude that the body cone rolls on the outside of the space cone.
(b) Take $I_1 = 2$, $I_3 = 1$ and $\theta_\omega = 10°$. Find $\theta_L$, the half-angles of both cones, $\Omega_b/\omega_3$ and $\Omega_s/\omega_3$.
(c) What does the sign of $\Omega_b$ tell an observer riding the body?

**B4.** (Your phone's flip rate.) Model a phone as a uniform box, $15 \times 7 \times 0.8$ cm.
(a) Find the three principal moments about the CM in units of $M\,\text{cm}^2/12$, and name the intermediate axis.
(b) Toss it at 2 rev/s about each axis in turn. Compute the small-wobble frequencies about the two stable axes and the growth rate $\sigma$ about the intermediate one. How long does a 1% misalignment take to grow to order one?
(c) Does the answer depend on the phone's mass? On its overall size, at fixed proportions? Explain.
(d) Do the experiment (over a bed).

**B5.** (Euler's third equation from the Lagrangian.) Let $I_1, I_2, I_3$ be arbitrary.
(a) Write $T$ in Euler angles using (3.23).
(b) Show $\partial\omega_1/\partial\psi = \omega_2$ and $\partial\omega_2/\partial\psi = -\omega_1$, hence $\partial T/\partial\psi = (I_1 - I_2)\,\omega_1\omega_2$.
(c) Write the Euler–Lagrange equation for $\psi$ with a potential $V(\phi, \theta, \psi)$ and show it is the third of Euler's equations (3.17), with $\tau_3 = -\partial V/\partial\psi$. Why is $-\partial V/\partial\psi$ the torque about $\hat{\mathbf e}_3$?
(d) Why does this computation show that $p_\psi$ is conserved only for a symmetric top? Which heavy asymmetric tops are still integrable, and who found the hard case?

**B6.** (The gyroscope's precession quadratic.)
(a) Derive (3.31) from the $\theta$ equation (3.28) by setting $\theta = \theta_0$.
(b) Derive it again without a Lagrangian: in steady precession $\mathbf L$ turns rigidly about $\hat{\mathbf z}$ at $\dot\phi$, so $d\mathbf L/dt = \dot\phi\,\hat{\mathbf z}\times\mathbf L$. Write $\mathbf L$ in terms of $\hat{\mathbf z}$ and $\hat{\mathbf e}_3$ and equate to the gravitational torque.
(c) For the toy gyroscope of Example 4.5 at $\theta_0 = 60°$, find the minimum spin $\omega_3$ for which steady precession exists, and the single precession rate at that spin.
(d) A bicycle wheel (a hoop of mass 2 kg, radius 0.33 m) spins at 5 rev/s on a light axle held at one end by a string, with the wheel's center 0.15 m from the string, axle horizontal. Find the precession rate and period. Why is (3.31) exact here, with no fast-top approximation?
(e) Show that for $\theta_0 > 90°$ (axle pointing below the horizontal) steady precession exists at every spin, with one forward and one backward rate.

**B7.** (Separatrices, and the flip period.) Take $I = (1, 2, 3)$ and $|\mathbf L| = 1$.
(a) On the separatrix $2T = |\mathbf L|^2/I_2$. Subtract the sphere from the ellipsoid and show the separatrices lie on the two planes $L_3 = \pm\sqrt3\,L_1$. So they are great circles. Give the general slope for any $I_1 < I_2 < I_3$.
(b) Predict, from the linear growth rate alone, how the period $P(\delta)$ of $\omega_2(t)$ for a start at $\boldsymbol\omega = (\delta, 1, \delta)$ grows as $\delta\to0$.
(c) (Code.) Add a function to `project.py` that measures $P(\delta)$ from the zero crossings of $\omega_2(t)$ for $\delta = 10^{-2}, 10^{-3}, 10^{-4}$. Check the prediction of (b).

**B8.** (Monte Carlo hemisphere; code.) A uniform solid hemisphere has radius $R$, flat face in the $xy$ plane, dome toward $+z$.
(a) By hand: show the CM is at $\bar z = \frac38R$, that every principal moment about the center of the flat face is $\frac25MR^2$, and that $I_1^{\rm cm} = I_2^{\rm cm} = \frac{83}{320}MR^2$.
(b) In `project.py`, sample the hemisphere (reflect `sample_sphere`), build $I$ about the face center and about the sample's own CM, assert the parallel-axis theorem exactly on the sample and the closed forms to 1.5%.
(c) The hemisphere rests dome-down on a rough table and rocks without slipping. Find the small-oscillation period for $R = 0.1$ m. (Rocking in a plane is a holonomic rolling constraint, so CM-3's method applies.)

## C. Challenge (do at least one)

**C1.** (The major-axis rule.) (a) Prove that among all motions with a given $|\mathbf L|$, the kinetic energy is smallest for spin about the axis of largest moment, and that $T_{\max}/T_{\min} = I_3/I_1$. (b) A spacecraft has flexible parts that dissipate energy but exert only internal torques. Argue that its spin ends up about the axis of largest moment, whatever axis it started about. (c) Explorer 1 (1958) was spun about its long axis. What happened, and why? (d) Dual-spin satellites spin a rotor about the minimum axis while the platform is nearly despun. In two sentences, explain why (b) does not doom them.

**C2.** (A satellite deploys a boom; code.) A satellite with principal moments $(400, 500, 600)$ kg m² spins at 1 rad/s about $\hat{\mathbf e}_3$, with a 1% transverse perturbation. Between $t = 10$ s and $t = 20$ s it extends two booms symmetrically along $\pm\hat{\mathbf e}_3$, which raises $I_1$ and $I_2$ by 150 kg m² each and leaves $I_3$ alone.
(a) From $d(I\boldsymbol\omega)/dt + \boldsymbol\omega\times(I\boldsymbol\omega) = \mathbf 0$ in body components, write Euler's equations for time-dependent $I$. Why does the booms' own motion relative to the body add nothing to $\mathbf L$?
(b) Show $|\mathbf L|$ is still conserved but $\dot T = -\frac12\boldsymbol\omega^\top\dot I\boldsymbol\omega$.
(c) Which axis is intermediate after deployment? Predict the growth rate.
(d) Integrate it (modify `euler_rhs`): show steady spin before deployment and a flip after, and check (b) numerically.

**C3.** (The gyrocompass.) A rotor with spin moment $I_3$ spins at $\omega_s$ about an axis held horizontal by a gimbal that is free to turn about the local vertical; $I_v$ is the moment about the vertical. The Earth turns at $\Omega_E$, and you are at latitude $\lambda$. Let $\alpha$ be the angle of the spin axis from north, measured counterclockwise seen from above.
(a) Use $(d\mathbf L/dt)_{\rm inertial} = (d\mathbf L/dt)_{\rm Earth} + \boldsymbol\Omega_E\times\mathbf L$ and the fact that the gimbal exerts no torque about the vertical to derive $I_v\ddot\alpha = -I_3\omega_s\Omega_E\cos\lambda\,\sin\alpha$ (neglect terms of order $\Omega_E^2$).
(b) Which way does $\mathbf L$ point at the stable equilibrium? Where on Earth does the compass fail?
(c) Find the oscillation period for $I_3 = 0.02$ kg m², $I_v = 0.01$ kg m², 20,000 rpm, $\lambda = 45°$, $\Omega_E = 7.2921\times10^{-5}$ rad/s.
