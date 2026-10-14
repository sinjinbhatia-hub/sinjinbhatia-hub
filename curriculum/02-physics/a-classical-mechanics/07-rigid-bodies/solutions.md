# CM-7 Solutions

Conventions as in the notes and `exercises.md`. Every number below is asserted or printed in Part 9 of `project.py`.

## A. Warm-up

**A1.** Put the CM at the origin, axis along $z$, and slice into disks of thickness $dz$, mass $dm = (M/h)\,dz$, radius $R$, for $-h/2\le z\le h/2$.

*Axis.* Each disk contributes $\frac12R^2\,dm$, so $I_3 = \frac12MR^2$.

*Perpendicular axis through the CM ($x$-axis).* A disk at height $z$ has $\frac14R^2\,dm$ about its own diameter parallel to $x$ (§3.7, from (3.14)), plus $z^2\,dm$ from the parallel-axis theorem (3.13), since its center is $z$ from the $x$-axis. So
$$I_1 = \tfrac14MR^2 + \frac Mh\int_{-h/2}^{h/2}z^2\,dz = \tfrac14MR^2 + \tfrac1{12}Mh^2.$$
*Spherical top.* $I_1 = I_3$ requires $\frac14R^2 + \frac1{12}h^2 = \frac12R^2$, so $h^2 = 3R^2$ and $h/R = \sqrt3 = 1.732$. A cylinder whose height is $0.87$ of its diameter spins like a sphere about its CM. Part 9 samples one and finds all three eigenvalues equal to $0.5$ (units $MR^2$) within 1%.

*End face.* The diameter of an end face is parallel to the $x$-axis and $h/2$ from it: $I_1^{\rm end} = \frac14MR^2 + \frac1{12}Mh^2 + M\frac{h^2}4 = \frac14MR^2 + \frac13Mh^2$.

**A2.** (a) From (3.3),
$$[\boldsymbol\omega]_\times = \begin{pmatrix}0 & -3 & 2\\ 3 & 0 & -1\\ -2 & 1 & 0\end{pmatrix},\qquad [\boldsymbol\omega]_\times\mathbf u = \begin{pmatrix}0 - 15 + 12\\ 12 + 0 - 6\\ -8 + 5 + 0\end{pmatrix} = \begin{pmatrix}-3\\ 6\\ -3\end{pmatrix}.$$
Directly, $\boldsymbol\omega\times\mathbf u = (2\cdot6 - 3\cdot5,\ 3\cdot4 - 1\cdot6,\ 1\cdot5 - 2\cdot4) = (-3, 6, -3)$. Same.

(b) $\mathbf u^\top[\boldsymbol\omega]_\times\mathbf u$ is a number, so it equals its transpose $\mathbf u^\top[\boldsymbol\omega]_\times^\top\mathbf u = -\mathbf u^\top[\boldsymbol\omega]_\times\mathbf u$. A number equal to its negative is zero. With $\mathbf u = \mathbf r$, this says $\mathbf r\cdot\dot{\mathbf r} = 0$: every point's velocity is perpendicular to its position vector from the rotation center, so $|\mathbf r|$ is constant. That is rigidity.

(c) $[\boldsymbol\omega]_\times^2\mathbf u = \boldsymbol\omega\times(\boldsymbol\omega\times\mathbf u) = \boldsymbol\omega(\boldsymbol\omega\cdot\mathbf u) - \mathbf u|\boldsymbol\omega|^2$ by BAC-CAB, which is $(\boldsymbol\omega\boldsymbol\omega^\top - |\boldsymbol\omega|^2\mathbb 1)\mathbf u$. True for every $\mathbf u$, so the matrices are equal. (Compare (3.8): for one particle, $I = -m[\mathbf r]_\times^2$.)

(d) $R = \begin{pmatrix}0&-1&0\\1&0&0\\0&0&1\end{pmatrix}$, $[\hat{\mathbf x}]_\times = \begin{pmatrix}0&0&0\\0&0&-1\\0&1&0\end{pmatrix}$. Then $R[\hat{\mathbf x}]_\times = \begin{pmatrix}0&0&1\\0&0&0\\0&1&0\end{pmatrix}$ and
$$R[\hat{\mathbf x}]_\times R^\top = \begin{pmatrix}0&0&1\\0&0&0\\0&1&0\end{pmatrix}\begin{pmatrix}0&1&0\\-1&0&0\\0&0&1\end{pmatrix} = \begin{pmatrix}0&0&1\\0&0&0\\-1&0&0\end{pmatrix} = [\hat{\mathbf y}]_\times,$$
and $R\hat{\mathbf x} = \hat{\mathbf y}$. ✓

**A3.** (a) For each mass, $r^2 = b^2$ and $\mathbf r\mathbf r^\top = b^2\hat{\mathbf n}\hat{\mathbf n}^\top$ (the sign of $\pm b$ squares away). By (3.8), $I = 2m(b^2\mathbb 1 - b^2\hat{\mathbf n}\hat{\mathbf n}^\top) = 2mb^2(\mathbb 1 - \hat{\mathbf n}\hat{\mathbf n}^\top)$. $\mathbb 1 - \hat{\mathbf n}\hat{\mathbf n}^\top$ is the projector onto the plane perpendicular to $\hat{\mathbf n}$: it sends $\hat{\mathbf n}\mapsto\mathbf 0$ and fixes every $\mathbf v\perp\hat{\mathbf n}$. So the principal moments are $0$ along the rod and $2mb^2$ (doubly degenerate) for every direction perpendicular to it. $I$ is singular, as Theorem 3.2(b) says it must be for a line of masses.

(b) $I\omega\hat{\mathbf z} = 2mb^2\omega(\hat{\mathbf z} - \cos\alpha\,\hat{\mathbf n}) = 2mb^2\omega(-\sin\alpha\cos\alpha,\ 0,\ 1 - \cos^2\alpha) = 2mb^2\omega\sin\alpha(-\cos\alpha, 0, \sin\alpha)$, which is §2's result.

(c) $2mb^2\omega\sin\alpha = 2(0.5)(0.04)(10)(0.5) = 0.2$, so $\mathbf L = (-0.17321, 0, 0.1)$ kg m²/s and $|\mathbf L| = 0.2$. $\cos\angle(\mathbf L, \hat{\mathbf z}) = 0.1/0.2$, so the angle is $60° = 90° - \alpha$, as it must be since $\mathbf L\perp$ rod. Bearing torque $|\boldsymbol\omega\times\mathbf L| = \omega|\mathbf L|\sin60° = 10(0.2)(0.86603) = 1.7321$ N m.

**A4.** From (3.17), $\dot\omega_k = [\tau_k + (I_i - I_j)\omega_i\omega_j]/I_k$ with $(k, i, j)$ cyclic:
$\dot\omega_1 = [0 + (3 - 4)(-2)(0.5)]/2 = 1/2 = 0.5$,
$\dot\omega_2 = [0 + (4 - 2)(0.5)(1)]/3 = 1/3 = 0.33333$,
$\dot\omega_3 = [1 + (2 - 3)(1)(-2)]/4 = 3/4 = 0.75$ rad/s².

*Energy.* $\dot T = \sum I_k\omega_k\dot\omega_k = 2(1)(0.5) + 3(-2)(\frac13) + 4(0.5)(0.75) = 1 - 2 + 1.5 = 0.5$ W, and $\boldsymbol\omega\cdot\boldsymbol\tau = 0.5\cdot1 = 0.5$. ✓ (The gyroscopic terms do no work.)

*Angular momentum.* $\mathbf L = (2, -6, 2)$ and $I\dot{\boldsymbol\omega} = (1, 1, 3)$, so $\frac d{dt}|\mathbf L|^2 = 2\mathbf L\cdot I\dot{\boldsymbol\omega} = 2(2 - 6 + 6) = 4$, and $2\mathbf L\cdot\boldsymbol\tau = 2\cdot2\cdot1 = 4$. ✓ In body components $\dot{\mathbf L} = \boldsymbol\tau - \boldsymbol\omega\times\mathbf L$, and $\mathbf L\cdot(\boldsymbol\omega\times\mathbf L) = 0$.

**A5.** (a) At $\theta = 0$, $R_x(0) = \mathbb 1$ and $R = R_z(\phi)R_z(\psi) = R_z(\phi + \psi)$. Both triples have $\phi + \psi = 1.1$, so both give $R_z(1.1)$.

(b) *About $\hat{\mathbf x}$:* $R_x(\Omega t) = R_z(0)R_x(\Omega t)R_z(0)$, so $(\phi, \theta, \psi) = (0, \Omega t, 0)$. *About $\hat{\mathbf y}$:* the third column of (3.21) is $\hat{\mathbf e}_3 = (\sin\phi\sin\theta, -\cos\phi\sin\theta, \cos\theta)$, and $R_y(\Omega t)\hat{\mathbf z} = (\sin\Omega t, 0, \cos\Omega t)$. So $\theta = \Omega t$, $\cos\phi = 0$, $\sin\phi = 1$: $\phi = 90°$. Then $R_z(90°)R_x(\theta)R_z(\psi)$ must equal $R_y(\theta)$. Since $R_z(90°)R_x(\theta)R_z(-90°)$ is the rotation by $\theta$ about $R_z(90°)\hat{\mathbf x} = \hat{\mathbf y}$, $\psi = -90°$. So $(\phi, \theta, \psi) = (90°, \Omega t, -90°)$.

(c) Both motions pass smoothly through $R = \mathbb 1$, yet for $t > 0$ one has $\phi = 0$ and the other $\phi = 90°$, for every small $t$. So $\phi$ (and $\psi$) as functions of $R$ jump depending on the direction from which $R$ leaves $\theta = 0$: no continuous choice exists there, and a tumble that passes near $\theta = 0$ needs $\dot\phi, \dot\psi\sim1/\theta$, which is the $\det J = -\sin\theta$ singularity of §3.11. An integrator storing $R$ (or a quaternion) evolves $\dot R = R[\boldsymbol\omega_{\rm body}]_\times$, whose right side is smooth everywhere, so it never notices.

## B. Core

**B1.** (a) Write $\mathbf L_{\rm space} = R\mathbf L_{\rm body}$. From §3.2, $R^\top\dot R = [\boldsymbol\omega_{\rm body}]_\times$, so $\dot R = R[\boldsymbol\omega_{\rm body}]_\times$. Differentiate:

1. $\dot{\mathbf L}_{\rm space} = \dot R\mathbf L_{\rm body} + R\dot{\mathbf L}_{\rm body} = R\big([\boldsymbol\omega_{\rm body}]_\times\mathbf L_{\rm body} + \dot{\mathbf L}_{\rm body}\big)$.
2. Newton–Euler: $\dot{\mathbf L}_{\rm space} = \boldsymbol\tau_{\rm space} = R\boldsymbol\tau_{\rm body}$.
3. Multiply by $R^\top$: $\dot{\mathbf L}_{\rm body} + \boldsymbol\omega_{\rm body}\times\mathbf L_{\rm body} = \boldsymbol\tau_{\rm body}$.
4. $\mathbf L_{\rm body} = I_{\rm body}\boldsymbol\omega_{\rm body}$ with $I_{\rm body}$ constant gives (3.16), and principal axes give (3.17).

This is (3.15) with the basis vectors hidden in $R$: the columns of $R$ are the $\hat{\mathbf e}_k$, and $\dot R = R[\boldsymbol\omega_{\rm body}]_\times$ is $\dot{\hat{\mathbf e}}_k = \boldsymbol\omega\times\hat{\mathbf e}_k$ for all three at once.

(b) In principal axes $|\mathbf L|^2 = \sum_kI_k^2\omega_k^2$ and $2T = \sum_kI_k\omega_k^2$. Then
$$|\mathbf L|^2 - 2TI_1 = \sum_kI_k(I_k - I_1)\,\omega_k^2 \ge 0,\qquad 2TI_3 - |\mathbf L|^2 = \sum_kI_k(I_3 - I_k)\,\omega_k^2 \ge 0,$$
because every term is a product of non-negative factors. Equality on the left holds iff $\omega_k = 0$ for every $k$ with $I_k > I_1$: spin about the axis of smallest moment (or within the degenerate plane if $I_1 = I_2$). Equality on the right iff the spin is about the axis of largest moment.

(c) $D = |\mathbf L|^2/2T$ is constant (both are conserved) and lies in $[I_1, I_3]$ by (b). On the motion, $\sum_kL_k^2/I_k = 2T = |\mathbf L|^2/D = \sum_kL_k^2/D$, so
$$L_1^2\Big(\frac1{I_1} - \frac1D\Big) = L_2^2\Big(\frac1D - \frac1{I_2}\Big) + L_3^2\Big(\frac1D - \frac1{I_3}\Big).$$
Take $I_1 < D < I_2 < I_3$. Both brackets on the right are positive, so the right side is positive unless $L_2 = L_3 = 0$, which would make $D = I_1$. Hence $L_1^2 > 0$ throughout: $L_1$ never vanishes, the polhode stays on one side of the plane $L_1 = 0$, and it loops around $+\hat{\mathbf e}_1$ or $-\hat{\mathbf e}_1$. For $I_2 < D < I_3$ rearrange the same identity as
$$L_3^2\Big(\frac1D - \frac1{I_3}\Big) = L_1^2\Big(\frac1{I_1} - \frac1D\Big) + L_2^2\Big(\frac1{I_2} - \frac1D\Big).$$
Now $I_1 < I_2 < D < I_3$ makes all three brackets positive, so the right side is positive unless $L_1 = L_2 = 0$ (which would make $D = I_3$). Hence $L_3$ never vanishes: a loop around $\pm\hat{\mathbf e}_3$. $D = I_2$ is the separatrix (B7) or steady spin about $\hat{\mathbf e}_2$.

**B2.** (a) *From the tensor.* Example 4.1 gives $I_{zz} = \frac23Ma^2$ about the corner; $I_{zz}$ is the moment about the $z$-axis, which is the edge. *Scalar parallel axis.* The edge is parallel to a CM axis at perpendicular distance $d_\perp = |(a/2, a/2)| = a/\sqrt2$, so $I = \frac16Ma^2 + M\frac{a^2}2 = \frac23Ma^2$. ✓

(b) $O = (0, 0, a/2)$, so the CM is at $\mathbf d = (a/2, a/2, 0)$ relative to $O$, $d^2 = a^2/2$, $\mathbf d\mathbf d^\top = \frac{a^2}4\begin{pmatrix}1&1&0\\1&1&0\\0&0&0\end{pmatrix}$. By (3.13) with $I^{\rm cm} = \frac16Ma^2\mathbb 1$,
$$I^{(O)} = Ma^2\begin{pmatrix}\frac16 + \frac12 - \frac14 & -\frac14 & 0\\ -\frac14 & \frac16 + \frac12 - \frac14 & 0\\ 0 & 0 & \frac16 + \frac12\end{pmatrix} = Ma^2\begin{pmatrix}\frac5{12} & -\frac14 & 0\\ -\frac14 & \frac5{12} & 0\\ 0 & 0 & \frac23\end{pmatrix}.$$
The $z$ direction is an eigenvector with moment $\frac23Ma^2$. The $2\times2$ block $\begin{pmatrix}p & q\\ q & p\end{pmatrix}$ has eigenvectors $(1, \pm1)/\sqrt2$ with eigenvalues $p\pm q$: $\frac5{12} - \frac14 = \frac16$ along $(1, 1, 0)/\sqrt2$ and $\frac5{12} + \frac14 = \frac23$ along $(1, -1, 0)/\sqrt2$. Principal moments $\frac16, \frac23, \frac23$ (units $Ma^2$): the small one points from $O$ to the CM.

(c) For a spherical top $I^{\rm cm} = I_0\mathbb 1$. By (3.13), $I^{(O)} = (I_0 + Md^2)\mathbb 1 - M\mathbf d\mathbf d^\top$. Then $I^{(O)}\mathbf d = (I_0 + Md^2 - Md^2)\mathbf d = I_0\mathbf d$, and for every $\mathbf v\perp\mathbf d$, $I^{(O)}\mathbf v = (I_0 + Md^2)\mathbf v$. So the principal moments are $I_0$ (along $\mathbf d$) and $I_0 + Md^2$ (twice): a symmetric top with its axis along $O\to$CM. ∎ Example 4.1 ($d^2 = \frac34a^2$: $\frac16, \frac{11}{12}, \frac{11}{12}$) and (b) ($d^2 = \frac12a^2$: $\frac16, \frac23, \frac23$) are two cases.

(d) Physical pendulum about the edge: $I\ddot\vartheta = -Mgd_\perp\sin\vartheta$ with $I = \frac23Ma^2$, $d_\perp = a/\sqrt2$. Small oscillations: $\omega^2 = \dfrac{Mga/\sqrt2}{\frac23Ma^2} = \dfrac{3g}{2\sqrt2\,a} = \dfrac{29.43}{0.28284} = 104.05$ s⁻², $\omega = 10.201$ rad/s, period $2\pi/\omega = 0.6160$ s.

**B3.** (a) Write $\boldsymbol\omega = \omega_\perp\hat{\mathbf u} + \omega_3\hat{\mathbf e}_3$ with $\hat{\mathbf u}\perp\hat{\mathbf e}_3$. Then $\mathbf L = I_1\omega_\perp\hat{\mathbf u} + I_3\omega_3\hat{\mathbf e}_3$ (§3.9 step 4), so both vectors lie in the plane of $\hat{\mathbf e}_3$ and $\hat{\mathbf u}$, on the same side of $\hat{\mathbf e}_3$, with $\tan\theta_\omega = \omega_\perp/\omega_3$ and $\tan\theta_L = I_1\omega_\perp/(I_3\omega_3) = (I_1/I_3)\tan\theta_\omega$. If $I_3 < I_1$ then $\tan\theta_L > \tan\theta_\omega$, so $\theta_L > \theta_\omega$: going from $\hat{\mathbf e}_3$ outward in that plane you meet $\boldsymbol\omega$ first, then $\mathbf L$. The body cone has axis $\hat{\mathbf e}_3$ and half-angle $\theta_\omega$; the space cone has axis $\mathbf L$ and half-angle $\theta_L - \theta_\omega$. They touch along $\boldsymbol\omega$, and their axes lie on opposite sides of the contact line, so neither cone contains the other: the body cone rolls on the outside of the fixed space cone. (For $I_3 > I_1$, $\theta_L < \theta_\omega$, both axes lie on the same side of $\boldsymbol\omega$, and the small space cone sits inside the body cone, as in §3.9.)

(b) $\tan\theta_L = 2\tan10° = 2(0.17633) = 0.35265$, $\theta_L = 19.425°$. Body cone half-angle $10°$; space cone half-angle $19.425° - 10° = 9.425°$. By (3.18), $\Omega_b/\omega_3 = (I_3 - I_1)/I_1 = -0.5$. By (3.19), $|\mathbf L| = \sqrt{(I_1\omega_\perp)^2 + (I_3\omega_3)^2} = \omega_3\sqrt{4\tan^210° + 1} = 1.0604\,\omega_3$, and $\Omega_s = |\mathbf L|/I_1 = 0.5302\,\omega_3$.

(c) $\Omega_b < 0$: an observer on the body sees $\boldsymbol\omega$ (and $\mathbf L$) circle $\hat{\mathbf e}_3$ in the sense *opposite* to the spin. For an oblate body ($\Omega_b > 0$) they circle in the same sense as the spin. This is the sign test for which cone is outside.

**B4.** (a) For a box with edges $a, b, c$, the moment about the axis along $a$ is $\frac M{12}(b^2 + c^2)$ (§3.7's plate result plus the third dimension). With $a = 15$, $b = 7$, $c = 0.8$ cm, in units of $M\,\text{cm}^2/12$:
long axis $49 + 0.64 = 49.64$; in-plane short axis (along the 7 cm width) $225 + 0.64 = 225.64$; axis through the screen $225 + 49 = 274$. The intermediate axis is the in-plane short axis, as §3.10 predicted.

(b) $\omega_0 = 2\cdot2\pi = 12.566$ rad/s. Order $I_1 = 49.64 < I_2 = 225.64 < I_3 = 274$. From (3.20):
- About the long axis: $\sqrt{-\lambda}/\omega_0 = \sqrt{\dfrac{(I_2 - I_1)(I_3 - I_1)}{I_2I_3}} = \sqrt{\dfrac{176\times224.36}{225.64\times274}} = \sqrt{0.63869} = 0.79918$, i.e. $10.04$ rad/s.
- About the screen normal: $\sqrt{\dfrac{(I_3 - I_1)(I_3 - I_2)}{I_1I_2}} = \sqrt{\dfrac{224.36\times48.36}{49.64\times225.64}} = \sqrt{0.96869} = 0.98422$, i.e. $12.37$ rad/s.
- About the short in-plane axis (Theorem 3.7): $\sigma/\omega_0 = \sqrt{\dfrac{(I_2 - I_1)(I_3 - I_2)}{I_1I_3}} = \sqrt{\dfrac{176\times48.36}{49.64\times274}} = \sqrt{0.62578} = 0.79106$, so $\sigma = 9.941$ s⁻¹.

A factor of 100 takes $\ln100/\sigma = 4.6052/9.941 = 0.463$ s, about $0.93$ of a turn. The phone flips before you catch it.

(c) Mass cancels from every ratio of moments. Scaling all lengths by $s$ multiplies every moment by $s^2$, which also cancels. So $\sigma/\omega_0$ depends only on the shape (the ratios $15 : 7 : 0.8$), and $\sigma$ depends on the toss only through $\omega_0$.

(d) Spin about the long axis or about the screen normal: it comes back as it left, with a small wobble. Spin end over end about the width axis: it comes back having also made a half-turn about its long axis, so the screen faces the other way. Tossing faster does not help: the flip takes a fixed number of turns, $\ln(1/\delta)/(2\pi\times0.791)$. Only a cleaner toss (smaller $\delta$) delays it, and only logarithmically.

**B5.** (a) From (3.23),
$$T = \tfrac12I_1(\dot\phi\sin\theta\sin\psi + \dot\theta\cos\psi)^2 + \tfrac12I_2(\dot\phi\sin\theta\cos\psi - \dot\theta\sin\psi)^2 + \tfrac12I_3(\dot\phi\cos\theta + \dot\psi)^2.$$
(b) $\partial\omega_1/\partial\psi = \dot\phi\sin\theta\cos\psi - \dot\theta\sin\psi = \omega_2$; $\partial\omega_2/\partial\psi = -\dot\phi\sin\theta\sin\psi - \dot\theta\cos\psi = -\omega_1$; $\omega_3$ does not contain $\psi$. So $\partial T/\partial\psi = I_1\omega_1\omega_2 + I_2\omega_2(-\omega_1) = (I_1 - I_2)\omega_1\omega_2$. (Geometrically: a change of $\psi$ turns the body about $\hat{\mathbf e}_3$, which rotates $(\omega_1, \omega_2)$ in the body's 1–2 plane.)

(c) Only $\omega_3$ contains $\dot\psi$, with $\partial\omega_3/\partial\dot\psi = 1$, so $\partial T/\partial\dot\psi = I_3\omega_3$. The Euler–Lagrange equation $\frac d{dt}\partial\mathcal L/\partial\dot\psi = \partial\mathcal L/\partial\psi$ with $\mathcal L = T - V$ reads
$$I_3\dot\omega_3 - (I_1 - I_2)\,\omega_1\omega_2 = -\frac{\partial V}{\partial\psi},$$
the third line of (3.17) with $\tau_3 = -\partial V/\partial\psi$. *Why that is the torque about $\hat{\mathbf e}_3$:* $R(\phi, \theta, \psi + \delta\psi) = R\,R_z(\delta\psi)$, and $R\,R_z(\delta\psi) = (RR_z(\delta\psi)R^\top)R$ is $R$ followed by a rotation by $\delta\psi$ about $R\hat{\mathbf z} = \hat{\mathbf e}_3$. A torque does work $\boldsymbol\tau\cdot\hat{\mathbf e}_3\,\delta\psi$ in that rotation, and for a conservative torque that work is $-\delta V$. So the generalized force $Q_\psi = -\partial V/\partial\psi = \tau_3$.

(d) $p_\psi = \partial\mathcal L/\partial\dot\psi = I_3\omega_3$ is conserved for every motion iff $\partial\mathcal L/\partial\psi = (I_1 - I_2)\omega_1\omega_2 - \partial V/\partial\psi$ vanishes identically. The first term depends on the velocities and the second does not, so each must vanish separately: $I_1 = I_2$ and a potential independent of $\psi$ (CM on the symmetry axis). For an asymmetric body $\omega_3$ is driven by $\omega_1\omega_2$, so there is no third integral of this kind. The heavy top with a fixed point is known to be integrable for general initial data in three cases: Euler's (no torque, any moments), Lagrange's (§3.12), and Sofia Kovalevskaya's (1888): $I_1 = I_2 = 2I_3$ with the CM in the body's equatorial plane, which won her the Prix Bordin.

**B6.** (a) Steady precession: $\theta \equiv \theta_0$, so $\ddot\theta = 0$, and (3.27) makes $\dot\phi$ constant too. (3.28) becomes $0 = I_1\dot\phi^2\sin\theta_0\cos\theta_0 - p_\psi\dot\phi\sin\theta_0 + Mg\ell\sin\theta_0$. Divide by $\sin\theta_0 \ne 0$: $I_1\cos\theta_0\,\dot\phi^2 - p_\psi\dot\phi + Mg\ell = 0$, which is (3.31).

(b) With $\dot\theta = 0$, (3.22) gives $\boldsymbol\omega = \dot\phi\hat{\mathbf z} + \dot\psi\hat{\mathbf e}_3$. Its part along $\hat{\mathbf e}_3$ is $\omega_3 = \dot\psi + \dot\phi\cos\theta$, and its part perpendicular to $\hat{\mathbf e}_3$ is $\dot\phi(\hat{\mathbf z} - \cos\theta\,\hat{\mathbf e}_3)$. With $I_1 = I_2$,
$$\mathbf L = I_1\dot\phi(\hat{\mathbf z} - \cos\theta\,\hat{\mathbf e}_3) + p_\psi\hat{\mathbf e}_3 = I_1\dot\phi\,\hat{\mathbf z} + (p_\psi - I_1\dot\phi\cos\theta)\,\hat{\mathbf e}_3.$$
In steady precession the whole configuration turns about $\hat{\mathbf z}$ at $\dot\phi$, so $d\mathbf L/dt = \dot\phi\hat{\mathbf z}\times\mathbf L = \dot\phi(p_\psi - I_1\dot\phi\cos\theta)\,\hat{\mathbf z}\times\hat{\mathbf e}_3$. The gravitational torque about the pivot is $\ell\hat{\mathbf e}_3\times(-Mg\hat{\mathbf z}) = Mg\ell\,\hat{\mathbf z}\times\hat{\mathbf e}_3$. Equate coefficients: $\dot\phi\,p_\psi - I_1\cos\theta\,\dot\phi^2 = Mg\ell$, the same quadratic.

(c) Real roots need $p_\psi^2 \ge 4I_1Mg\ell\cos\theta_0$, i.e.
$$\omega_3 \ge \frac{\sqrt{4I_1Mg\ell\cos\theta_0}}{I_3} = \frac{\sqrt{4(1.825\times10^{-4})(0.03924)(0.5)}}{4.5\times10^{-5}} = \frac{3.7845\times10^{-3}}{4.5\times10^{-5}} = 84.10\ \text{rad/s} = 803.1\ \text{rpm}.$$
At that spin the double root is $\dot\phi = p_\psi/(2I_1\cos\theta_0) = 3.7845\times10^{-3}/1.825\times10^{-4} = 20.74$ rad/s (equivalently $2Mg\ell/p_\psi$). It is $1/\sqrt2$ of the sleeping threshold 118.9 rad/s because $\cos60° = \frac12$.

(d) $I_3 = MR^2 = 2(0.33)^2 = 0.2178$ kg m² (a hoop), $\omega_3 = 2\pi\cdot5 = 31.416$ rad/s, $p_\psi = 6.8424$ kg m²/s, $Mg\ell = 2(9.81)(0.15) = 2.943$ N m. At $\theta_0 = 90°$ the $\dot\phi^2$ term of (3.31) has coefficient $I_1\cos90° = 0$, so $\dot\phi = Mg\ell/p_\psi = 0.43011$ rad/s exactly, period $2\pi/0.43011 = 14.61$ s. It is exact because in (b) $\mathbf L = I_1\dot\phi\hat{\mathbf z} + p_\psi\hat{\mathbf e}_3$: the transverse part of $\mathbf L$ is vertical and does not rotate, so only the spin part needs torque. (At $90°$, $\omega_3 = \dot\psi$, so "5 rev/s relative to the axle" is $\omega_3$.)

(e) For $\cos\theta_0 < 0$ the discriminant $p_\psi^2 - 4I_1Mg\ell\cos\theta_0 > 0$ for every $p_\psi$, so two real roots always exist. Their product is $Mg\ell/(I_1\cos\theta_0) < 0$: one positive, one negative. At $p_\psi = 0$ they are $\pm\sqrt{Mg\ell/(I_1|\cos\theta_0|)}$, a conical pendulum going either way round. A hanging gyroscope needs no minimum spin.

**B7.** (a) With $|\mathbf L| = 1$ and $2T = 1/I_2 = \frac12$: ellipsoid $L_1^2 + \frac12L_2^2 + \frac13L_3^2 = \frac12$, sphere $L_1^2 + L_2^2 + L_3^2 = 1$. Subtract $\frac12\times$ sphere from the ellipsoid: $\frac12L_1^2 - \frac16L_3^2 = 0$, so $L_3 = \pm\sqrt3\,L_1$. Each is a plane through the origin (it contains $\hat{\mathbf e}_2$), and a plane through the center meets the sphere in a great circle. The planes make $\arctan\sqrt3 = 60°$ with $\hat{\mathbf e}_1$ in the 1–3 plane. In general, subtracting $1/I_2$ times the sphere: $L_1^2(\frac1{I_1} - \frac1{I_2}) = L_3^2(\frac1{I_2} - \frac1{I_3})$, so
$$\frac{L_3}{L_1} = \pm\sqrt{\frac{I_3(I_2 - I_1)}{I_1(I_3 - I_2)}}.$$
Part 9 runs a trajectory from near $\hat{\mathbf e}_2$ and checks $|L_3/L_1|\to\sqrt3$ along it once it leaves the saddle.

(b) Near $\pm\hat{\mathbf e}_2$ the linearized motion (3.20) is a saddle: one mode grows as $e^{\sigma t}$, one decays as $e^{-\sigma t}$. Call their amplitudes $y$ and $x$. The product $xy$ is constant in the linear flow (it is $e^{\sigma t}e^{-\sigma t}$ times a constant). A trajectory enters a fixed small box $|x|, |y|\le\varepsilon$ at $x = \varepsilon$, $y = c/\varepsilon$ and leaves at $y = \varepsilon$, $x = c/\varepsilon$, which takes $\frac1\sigma\ln(\varepsilon^2/c)$. The constant $c$ is set by how far the conserved quantities sit from the separatrix: $|\mathbf L|^2 - 2TI_2 = \sum I_k(I_k - I_2)\omega_k^2 = (-1 + 3)\delta^2 = 2\delta^2$ for the start $(\delta, 1, \delta)$, so $c\propto\delta^2$. Each passage near a saddle then takes $\frac2\sigma\ln(1/\delta) + $ const. Outside the boxes the trajectory follows the separatrix closely and the time is independent of $\delta$. One full period of $\omega_2$ passes near $+\hat{\mathbf e}_2$ and near $-\hat{\mathbf e}_2$ once each, so
$$P(\delta) \approx \frac4\sigma\ln\frac1\delta + C,\qquad P(\delta/10) - P(\delta) \approx \frac{4\ln10}\sigma = 4\sqrt3\ln10 = 15.953.$$
(c) The code (Part 9, `flip_period`) uses `solve_ivp` with an event on $\omega_2 = 0$ and takes the time between the first and third zero crossings:
```python
def flip_period(delta, I=I123):
    ev = lambda t, w, I: w[1]
    sol = solve_ivp(euler_rhs, (0, 400), [delta, 1.0, delta], args=(I,), method="DOP853",
                    rtol=1e-12, atol=1e-14, events=ev)
    tz = sol.t_events[0]
    return tz[2] - tz[0]
```
It gives $P = 39.106, 55.062, 71.015$ for $\delta = 10^{-2}, 10^{-3}, 10^{-4}$. Differences: $15.956$ and $15.953$, against the prediction $15.953$. The logarithm is why the Dzhanibekov flip looks periodic and regular: halving the misalignment adds only $\frac4\sigma\ln2$ to the period.

**B8.** (a) *CM.* The slice at height $z$ is a disk of area $\pi(R^2 - z^2)$; volume $V = \frac23\pi R^3$. So $\bar z = \frac1V\int_0^R z\,\pi(R^2 - z^2)\,dz = \frac{\pi(R^4/2 - R^4/4)}{2\pi R^3/3} = \frac38R$.
*Moments about the face center.* The full sphere of mass $2M$ has $\frac25(2M)R^2$ about every diameter. The integrands $y^2 + z^2$ etc. are even under $z\to -z$, so each hemisphere contributes half: $I_{xx} = I_{yy} = I_{zz} = \frac25MR^2$. The products vanish: $xy$ and $xz$ are odd under $x\to -x$, $yz$ under $y\to -y$.
*About the CM.* $\mathbf d = (0, 0, \frac38R)$. (3.13) run backwards: $I^{\rm cm} = I^{(O)} - M(d^2\mathbb 1 - \mathbf d\mathbf d^\top)$, so $I_3^{\rm cm} = \frac25MR^2$ (the $zz$ correction is $d^2 - d_z^2 = 0$) and $I_1^{\rm cm} = \frac25MR^2 - \frac9{64}MR^2 = \frac{128 - 45}{320}MR^2 = \frac{83}{320}MR^2 = 0.259375\,MR^2$.

(b) In Part 9:
```python
hemi = sample_sphere(1.0, N); hemi[:, 2] = np.abs(hemi[:, 2])   # reflect the lower half up
I_face = I_per_mass(hemi); c_h = hemi.mean(axis=0); I_hcm = I_per_mass(hemi - c_h)
assert np.allclose(I_face, parallel_axis(I_hcm, 1.0, -c_h), atol=1e-12)
assert np.allclose(np.diag(I_face), 0.4, rtol=0.015) and np.isclose(c_h[2], 3/8, rtol=0.01)
```
Reflection keeps the sample uniform because $z\mapsto|z|$ maps the uniform ball two-to-one onto the uniform hemisphere. The sample gives $I_1^{\rm cm}\approx0.259$ (units $MR^2$).

(c) The sphere's center $O$ stays at height $R$ while rolling (the dome is part of a sphere about $O$), and the CM sits $\frac38R$ below $O$ at rest. At tilt $\vartheta$ the CM height is $R - \frac38R\cos\vartheta$, so $V = -\frac38MgR\cos\vartheta + $ const. Rolling without slipping makes the contact point the instantaneous axis, which is $\frac58R$ from the CM at small $\vartheta$, so $T = \frac12I_P\dot\vartheta^2$ with $I_P = I_1^{\rm cm} + M(\frac58R)^2 = (\frac{83}{320} + \frac{125}{320})MR^2 = 0.65\,MR^2$ (to leading order in $\vartheta$). Small oscillations: $\omega^2 = \dfrac{\frac38MgR}{0.65MR^2} = 0.57692\,g/R = 56.596$ s⁻² for $R = 0.1$ m, $\omega = 7.523$ rad/s, period $0.8352$ s.

## C. Challenge

**C1.** (a) In principal axes $T = \sum_kL_k^2/(2I_k)$. With $I_1 \le I_2 \le I_3$, $\frac1{I_k}\ge\frac1{I_3}$, so $T \ge \sum L_k^2/(2I_3) = |\mathbf L|^2/(2I_3)$, with equality iff $L_k = 0$ whenever $I_k < I_3$: spin about the major axis. Likewise $T \le |\mathbf L|^2/(2I_1)$, with equality for spin about the minor axis. So $T_{\max}/T_{\min} = I_3/I_1$ at fixed $|\mathbf L|$.

(b) Internal forces exert no net torque, so $\mathbf L$ is conserved exactly. Dissipation lowers $T$ monotonically while any internal motion (flexing, sloshing) is excited, and a wobbling body excites internal motion at the wobble frequency. Energy keeps draining until the body sits where no wobble is left, and the only state with given $|\mathbf L|$ that cannot lose more energy is the minimum of $T$: spin about the major axis. Spin about the minor axis is stable for a rigid body (Theorem 3.7) but is the energy *maximum* at fixed $|\mathbf L|$, so dissipation drives it away.

(c) Explorer 1 (1958) was a long, thin cylinder spun about its long axis, which is its axis of *least* moment. Its flexible whip antennas dissipated energy, and it was soon tumbling about a transverse axis. Bracewell and Garriott's analysis (*Nature*, 1958) is the standard account and the origin of the "major-axis rule" in spacecraft design.

(d) The rule assumes one body with dissipation. A dual-spin craft is two bodies spinning at different rates; when the dissipation is on the slowly turning platform, the energy-sink argument has to be redone with both spin rates, and it shows a fast minor-axis rotor can be stable. The design rule becomes "put the dampers on the despun part" instead of "always spin about the major axis."

**C2.** (a) In body components, $\mathbf L_{\rm body} = I(t)\boldsymbol\omega$, and B1(a) gives $\dot{\mathbf L}_{\rm body} + \boldsymbol\omega\times\mathbf L_{\rm body} = \mathbf 0$. So
$$I\dot{\boldsymbol\omega} = -\dot I\boldsymbol\omega - \boldsymbol\omega\times(I\boldsymbol\omega).$$
The full angular momentum about the CM is $I\boldsymbol\omega + \sum m\,\mathbf r\times\mathbf v_{\rm rel}$, where $\mathbf v_{\rm rel}$ is boom velocity relative to the body. Each boom element moves along its own position vector (radially along $\pm\hat{\mathbf e}_3$ through the CM), so $\mathbf r\times\mathbf v_{\rm rel} = \mathbf 0$. Deploying symmetrically also keeps the CM fixed.

(b) $\frac d{dt}|\mathbf L|^2 = 2\mathbf L\cdot\dot{\mathbf L}_{\rm body} = -2\mathbf L\cdot(\boldsymbol\omega\times\mathbf L) = 0$. For $T = \frac12\boldsymbol\omega\cdot\mathbf L$: $\dot T = \frac12(\dot{\boldsymbol\omega}\cdot\mathbf L + \boldsymbol\omega\cdot\dot{\mathbf L})$. The second term is $-\frac12\boldsymbol\omega\cdot(\boldsymbol\omega\times\mathbf L) = 0$. In the first, $\dot{\boldsymbol\omega}\cdot\mathbf L = \dot{\boldsymbol\omega}^\top I\boldsymbol\omega = (I\dot{\boldsymbol\omega})\cdot\boldsymbol\omega = -\boldsymbol\omega^\top\dot I\boldsymbol\omega - \boldsymbol\omega\cdot(\boldsymbol\omega\times I\boldsymbol\omega) = -\boldsymbol\omega^\top\dot I\boldsymbol\omega$. So $\dot T = -\frac12\boldsymbol\omega^\top\dot I\boldsymbol\omega$. Here $\dot I = \dot s\cdot150\,\operatorname{diag}(1, 1, 0)$, so only the small transverse $\omega_1, \omega_2$ feed it, and $T$ barely changes. The deployment mechanism absorbs that energy.

(c) After deployment $I = (550, 650, 600)$: $\hat{\mathbf e}_3$ ($600$) is now intermediate. Spin about it at $\omega_3 = 1$ rad/s ($L_3 = I_3\omega_3$ and $I_3$ is unchanged, so $\omega_3$ stays 1 to first order) is unstable with
$$\sigma = \omega_3\sqrt{\frac{(600 - 550)(650 - 600)}{550\times650}} = \sqrt{0.0069930} = 0.08362\ \text{s}^{-1}.$$
Before deployment, $\hat{\mathbf e}_3$ is the major axis and the wobble frequency is $\sqrt{(600 - 400)(600 - 500)/(400\times500)} = \sqrt{0.1} = 0.3162$ rad/s.

(d) Part 9's `boom_rhs` uses a smooth ramp $s(t)$ from 0 to 1 on $[10, 20]$ s:
```python
def boom_I(t):
    u = np.clip((t - 10) / 10, 0, 1)
    s, ds = u * u * (3 - 2 * u), 6 * u * (1 - u) / 10                  # smoothstep and its derivative
    return np.array([400 + 150 * s, 500 + 150 * s, 600.0]), np.array([150 * ds, 150 * ds, 0.0])

def boom_rhs(t, w):
    I, Idot = boom_I(t)
    L = I * w
    return (-Idot * w - np.cross(w, L)) / I
```
Started at $\boldsymbol\omega = (0.01, 0.01, 1)$ and run to 250 s: during $0 \le t \le 10$ the transverse components stay below $0.017$ (a steady wobble), after deployment $\omega_3$ reverses sign (the flip). $|\mathbf L|$ stays constant to $10^{-9}$, and $T(t) - T(0)$ matches $-\int\frac12\boldsymbol\omega^\top\dot I\boldsymbol\omega\,dt$ to $10^{-9}$. A straight-line fit of $\ln|\omega_1|$ over the window $0.02 < |\omega_1| < 0.1$ after deployment gives $0.087$ s⁻¹ against the predicted $0.0836$; the 4% excess is the decaying mode, which has not fully died out at the start of the window.

**C3.** (a) Use local axes east $\hat{\mathbf x}$, north $\hat{\mathbf y}$, up $\hat{\mathbf u} = \hat{\mathbf z}$ (right-handed). The Earth's angular velocity is $\boldsymbol\Omega_E = \Omega_E(0, \cos\lambda, \sin\lambda)$. With $\alpha$ counterclockwise from north seen from above, the spin axis is $\hat{\mathbf n} = (-\sin\alpha, \cos\alpha, 0)$ and the gimbal turns about $\hat{\mathbf u}$ at $\dot\alpha$ relative to the Earth.

1. Angular momentum: $\mathbf L = I_3\omega_s\hat{\mathbf n} + I_v(\dot\alpha + \Omega_E\sin\lambda)\hat{\mathbf u} + O(I\Omega_E)$ horizontal terms, where $\omega_s$ is the total spin about $\hat{\mathbf n}$. The dropped terms come from $\boldsymbol\Omega_E$'s horizontal component and contribute only at order $\Omega_E^2$ below.
2. No torque about the vertical: $\hat{\mathbf u}\cdot(d\mathbf L/dt)_{\rm inertial} = 0$. Write $(d\mathbf L/dt)_{\rm inertial} = (d\mathbf L/dt)_{\rm Earth} + \boldsymbol\Omega_E\times\mathbf L$ (this is (3.15) with the Earth as the "body").
3. $\hat{\mathbf u}$ is fixed in the Earth frame, so $\hat{\mathbf u}\cdot(d\mathbf L/dt)_{\rm Earth} = \frac d{dt}(\hat{\mathbf u}\cdot\mathbf L) = I_v\ddot\alpha$.
4. $\hat{\mathbf u}\cdot(\boldsymbol\Omega_E\times\mathbf L) = \mathbf L\cdot(\hat{\mathbf u}\times\boldsymbol\Omega_E)$. Now $\hat{\mathbf u}\times\boldsymbol\Omega_E = \Omega_E\cos\lambda\,(\hat{\mathbf z}\times\hat{\mathbf y}) = -\Omega_E\cos\lambda\,\hat{\mathbf x}$, and $\mathbf L\cdot\hat{\mathbf x} = -I_3\omega_s\sin\alpha$ (the vertical part of $\mathbf L$ has no $\hat{\mathbf x}$ component). So this term is $I_3\omega_s\Omega_E\cos\lambda\sin\alpha$.
5. Adding: $I_v\ddot\alpha + I_3\omega_s\Omega_E\cos\lambda\,\sin\alpha = 0$. ∎

(b) For $I_3\omega_s > 0$ this is a pendulum equation with stable equilibrium $\alpha = 0$: the spin axis, and $\mathbf L$, point north, parallel to the horizontal component of $\boldsymbol\Omega_E$. ($\alpha = \pi$ is the inverted, unstable position.) The restoring coefficient vanishes at the poles, $\cos\lambda = 0$, where "north" is undefined, and weakens at high latitude.

(c) $\omega_s = 20000\cdot2\pi/60 = 2094.4$ rad/s. $\omega_{\rm gc}^2 = I_3\omega_s\Omega_E\cos\lambda/I_v = (0.02)(2094.4)(7.2921\times10^{-5})(0.70711)/0.01 = 0.21599$ s⁻², $\omega_{\rm gc} = 0.4647$ rad/s, period $13.52$ s. Undamped, it would swing about north forever; real compasses add damping, and marine gyrocompasses are additionally made pendulous and tuned to the 84-minute Schuler period so that the ship's accelerations do not deflect them.
