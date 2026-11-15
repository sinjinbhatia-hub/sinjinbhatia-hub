# CM-7: Rigid Bodies — the Inertia Tensor, Euler's Equations, Euler Angles, the Spinning Top

**Time:** ~12 hours · **Prerequisites:** CM-3, LA-8, LA-12 · **Unlocks:** QM-5 (angular momentum parallels), CM-8 (the top as a Hamiltonian system), ME dynamics (rotor balancing, gyroscopic loads, attitude control)

> **Conventions in this module.** The **space frame** has fixed orthonormal axes $\hat{\mathbf x}, \hat{\mathbf y}, \hat{\mathbf z}$; components in it carry subscripts $x, y, z$. The **body frame** has axes $\hat{\mathbf e}_1, \hat{\mathbf e}_2, \hat{\mathbf e}_3$ fixed in the body (principal axes whenever we say so); components in it carry subscripts $1, 2, 3$. The **attitude** is the rotation matrix $R(t) \in SO(3)$ with $\mathbf r_{\rm space} = R\,\mathbf r_{\rm body}$, so the columns of $R$ are the body axes written in space components. $\boldsymbol\omega$ is the angular velocity, $[\mathbf v]_\times$ is the antisymmetric matrix with $[\mathbf v]_\times\mathbf u = \mathbf v\times\mathbf u$. $I$ (entries $I_{ij}$) is the inertia tensor, $I_1, I_2, I_3$ the principal moments. $\mathbf L$ is angular momentum, $T$ kinetic energy, $\boldsymbol\tau$ torque, $\mathcal L$ the Lagrangian. $I$, $\mathbf L$, $\boldsymbol\tau$ are about the center of mass (CM) or a fixed pivot. **Euler angles** $(\phi, \theta, \psi)$ use the **$z$-$x'$-$z''$ convention** of Goldstein, $R = R_z(\phi)R_x(\theta)R_z(\psi)$; some texts rotate about $y'$ in the second step ($z$-$y'$-$z''$, which shifts $\phi$ by $+90°$ and $\psi$ by $-90°$ and leaves the top's Lagrangian unchanged), and aerospace yaw-pitch-roll ($z$-$y'$-$x''$) differs again. Goldstein's $A$ is our $R^\top$. SI units.

## 0. Where this sits

You need CM-3 (writing $\mathcal L = T - V$ in generalized coordinates and grinding out Euler–Lagrange equations, plus SymPy's `lagrangian_to_eom`), CM-4 Theorem 3.1 (cyclic coordinates give conserved momenta), LA-8 Theorem 3.11 (the spectral theorem for real symmetric matrices) and LA-12 (the $(0,2)$ tensor law (3.19) and Appendix A on $SO(3)$, especially (A.1)). LA-12 Example 4.6 already computed one inertia tensor with `einsum`; here it becomes dynamics.

This module makes possible: QM-5, where the same algebra of rotations, $[L_i, L_j] = \varepsilon_{ijk}L_k$, becomes quantum angular momentum and the symmetric top becomes the molecular rotor; CM-8, whose §3.11 lists the symmetric top among the integrable systems (as many conserved quantities as degrees of freedom); and every ME course that rotates something.

## 1. The problem that forced this idea

By the 1740s Newton's laws handled point masses and bodies on a fixed axle. Nobody could predict how a free, lopsided solid tumbles: throw a book spinning about a tilted axis and it wobbles in a way point-particle intuition does not predict. Euler worked on this from about 1750 to 1765 (*Theoria motus corporum solidorum seu rigidorum*, 1765). He found that angular momentum is generally **not parallel** to angular velocity, that every body nonetheless has three perpendicular axes about which it can spin without wobbling (J. A. Segner reached these "free axes" around 1755), and he wrote down the equations of the tumbling and the three angles we still use for orientation.

The payoffs came in a long line. Newton had explained the **precession of the equinoxes** (found by Hipparchus around 130 BCE) in the *Principia* (1687) as the Sun and Moon twisting Earth's equatorial bulge. Euler predicted that a rigid Earth spinning slightly off its figure axis would wobble with a period of about ten months; Seth Chandler found the wobble in 1891 at about fourteen months, and the mismatch measures the Earth's elasticity. Lagrange solved the heavy symmetric top (1788), Poinsot gave the geometric picture of free rotation (1834), and Foucault built and named the **gyroscope** (1852) to show the Earth's rotation. In 1985 Vladimir Dzhanibekov watched a wing nut on Salyut 7 flip end over end every few seconds: the **tennis-racket instability**, predicted by Euler's equations two centuries earlier.

For a mechanical engineer this is not history. Every rotating part loads its bearings through $\mathbf L \nparallel \boldsymbol\omega$ unless it is dynamically balanced; every satellite attitude, robot-arm joint torque and IMU is a rigid-body problem. Example 4.2 puts a number on the bearing load of a plate spun about the wrong axis; it is larger than you will guess.

## 2. Intuition

### A dumbbell on a bent shaft

Fix a massless rod with masses $m$ at distance $b$ either side of its midpoint to a vertical shaft at that midpoint, tilted at angle $\alpha$ from the shaft. Spin the shaft at $\boldsymbol\omega = \omega\hat{\mathbf z}$. At an instant, the upper mass is at $\mathbf r = b(\sin\alpha, 0, \cos\alpha)$ and moves with $\mathbf v = \boldsymbol\omega\times\mathbf r = \omega b\sin\alpha\,\hat{\mathbf y}$. Its angular momentum is $m\mathbf r\times\mathbf v = mb^2\omega\sin\alpha\,(-\cos\alpha, 0, \sin\alpha)$. The lower mass, at $-\mathbf r$ with velocity $-\mathbf v$, contributes the same. So
$$\mathbf L = 2mb^2\omega\sin\alpha\,(-\cos\alpha,\ 0,\ \sin\alpha).$$
This is not along $\hat{\mathbf z}$ unless $\alpha = 0$ or $90°$. It cannot be: each $m\mathbf r\times\mathbf v$ is perpendicular to $\mathbf r$, so $\mathbf L$ is perpendicular to the rod. As the rod turns, $\mathbf L$ turns with it, so $d\mathbf L/dt \ne \mathbf 0$ even though $\omega$ is constant. Something must supply that torque: the bearings. That is the whole problem of dynamic imbalance, in one sketch.

### $\mathbf L$ is a linear function of $\boldsymbol\omega$, so it is a matrix

Each velocity $\mathbf v_a = \boldsymbol\omega\times\mathbf r_a$ is linear in $\boldsymbol\omega$, and $\mathbf L = \sum_a m_a\mathbf r_a\times\mathbf v_a$ is linear in the velocities. A linear map $\mathbb R^3\to\mathbb R^3$ is a $3\times 3$ matrix (LA-4). So there is a matrix $I$, depending only on where the mass is, with $\mathbf L = I\boldsymbol\omega$. The dumbbell shows its off-diagonal entries are not zero.

Now the kinetic energy. Using the scalar triple product,
$$T = \tfrac12\sum_a m_a\mathbf v_a\cdot(\boldsymbol\omega\times\mathbf r_a) = \tfrac12\,\boldsymbol\omega\cdot\sum_a m_a\mathbf r_a\times\mathbf v_a = \tfrac12\,\boldsymbol\omega\cdot\mathbf L = \tfrac12\,\boldsymbol\omega^\top I\boldsymbol\omega.$$
$T > 0$ whenever the body actually moves. So $I$ defines a positive quadratic form. If $I$ is also symmetric (you will see in §3.3 that it is, by inspection), then you already know what LA-8 says about symmetric positive definite matrices.

**Make the guess now.** What does the spectral theorem promise about $I$? What happens to $\mathbf L$ when $\boldsymbol\omega$ points along one of the eigenvectors? What torque does the bearing need then? For the dumbbell, which directions are these?

(Answer, to check after you guess: three mutually perpendicular axes, the **principal axes**, along which $\mathbf L = I_k\boldsymbol\omega$ is parallel to $\boldsymbol\omega$, so a body spinning steadily about one needs no torque at all. For the dumbbell: along the rod, and any direction perpendicular to it.)

### Two more guesses

**The phone toss.** Toss your phone spinning about each of its three symmetry axes in turn (long axis, short in-plane axis, axis through the screen). In one toss it comes back flipped over. Guess which before §3.10. (Hint: order the three moments; the odd one out is neither largest nor smallest.)

**Why the top doesn't fall.** A fast top with $\mathbf L\approx I_3\omega_3\hat{\mathbf e}_3$ tilted at $\theta$ feels a gravitational torque $Mg\ell\sin\theta$, horizontal and perpendicular to $\mathbf L$. Since $d\mathbf L/dt = \boldsymbol\tau$, the tip of $\mathbf L$ moves sideways, not down. The horizontal part of $\mathbf L$ has length $I_3\omega_3\sin\theta$ and is pushed at rate $Mg\ell\sin\theta$, so it turns at $\Omega = Mg\ell/(I_3\omega_3)$, independent of tilt. §3.12 shows this is the slow root of a quadratic.

## 3. Formal development

### 3.1 Six degrees of freedom; Chasles

A **rigid body** is a set of particles whose mutual distances $|\mathbf r_a - \mathbf r_b|$ never change. Three non-collinear particles fix the body: 9 coordinates minus 3 fixed distances leaves **6 degrees of freedom**, 3 for the position of a reference point $O$ and 3 for the orientation. Write every particle as
$$\mathbf r_a(t) = \mathbf R_O(t) + R(t)\,\mathbf b_a, \tag{3.1}$$
where $\mathbf b_a$ is the particle's constant position in body coordinates and $R(t)\in SO(3)$ ($R^\top R = I$, $\det R = 1$; LA-12 A.1). Every rotation has an axis (LA-12 A.4), which gives **Chasles' theorem**: the most general rigid displacement is a translation plus a rotation about an axis, and in fact a *screw*, a rotation about some line combined with a translation along that same line. We will not need the screw form; we need $R(t)$ and its derivative.

### 3.2 The angular velocity vector exists

Nothing guarantees in advance that a tumbling body has *a* vector $\boldsymbol\omega$. Each particle has its own velocity. The claim is that one vector serves for all of them.

▶ **Try it first.** Differentiate $RR^\top = I$. What can you say about the matrix $\dot RR^\top$? Then write $\dot{\mathbf r}$ for a body-fixed point $\mathbf r = R\mathbf b$ using it.

**Theorem 3.1 (Angular velocity).** Let $R(t)\in SO(3)$ be differentiable. Then $\Omega := \dot RR^\top$ is antisymmetric, so $\Omega = [\boldsymbol\omega]_\times$ for a unique vector $\boldsymbol\omega(t)$, and every body-fixed vector $\mathbf r(t) = R(t)\mathbf b$ satisfies $\dot{\mathbf r} = \boldsymbol\omega\times\mathbf r$.

*Derivation.*

1. $RR^\top = I$ for all $t$. Differentiate with the product rule: $\dot RR^\top + R\dot R^\top = 0$.
2. $R\dot R^\top = (\dot RR^\top)^\top$, so step 1 reads
$$\Omega + \Omega^\top = 0. \tag{3.2}$$
3. An antisymmetric $3\times3$ matrix has three independent entries. Name them so that
$$\Omega = \begin{pmatrix}0 & -\omega_z & \omega_y\\ \omega_z & 0 & -\omega_x\\ -\omega_y & \omega_x & 0\end{pmatrix} =: [\boldsymbol\omega]_\times. \tag{3.3}$$
Multiplying out, $[\boldsymbol\omega]_\times\mathbf u = (\omega_yu_z - \omega_zu_y,\ \omega_zu_x - \omega_xu_z,\ \omega_xu_y - \omega_yu_x) = \boldsymbol\omega\times\mathbf u$.
4. For $\mathbf r = R\mathbf b$ with $\mathbf b$ constant: $\dot{\mathbf r} = \dot R\mathbf b = \dot RR^\top R\mathbf b = \Omega\mathbf r$, so
$$\dot{\mathbf r} = \boldsymbol\omega\times\mathbf r. \tag{3.4}$$ ∎

■ **What just happened.** Orthogonality alone forces the velocity field of a rotating body to be a cross product with one vector. This is LA-12 (A.1), "the tangent space of $SO(3)$ is the antisymmetric matrices," moved from $R = I$ to the current $R$ and given a physical name. In $n$ dimensions $\Omega$ has $n(n-1)/2$ entries, which equals $n$ only for $n = 3$; that is why angular velocity is a vector only in 3D. Note that $\boldsymbol\omega$ is *not* the time derivative of any angle vector.

**Body components.** $\boldsymbol\omega_{\rm body} = R^\top\boldsymbol\omega_{\rm space}$. Because a rotation preserves cross products, $R(\mathbf a\times\mathbf b) = R\mathbf a\times R\mathbf b$, which says $R[\mathbf a]_\times R^\top = [R\mathbf a]_\times$. Hence $R^\top\dot R = R^\top\Omega R = [\boldsymbol\omega_{\rm body}]_\times$. The project extracts $\boldsymbol\omega$ both ways from finite differences of $R(t)$.

**Finite rotations do not commute; infinitesimal ones do.** A $90°$ turn about $\hat{\mathbf x}$ then $90°$ about $\hat{\mathbf y}$ differs from the reverse order (the 3D version of LA-1 Example 4.3; the project checks it). But by Rodrigues' formula (LA-12 (A.4)), $R(\mathbf n, \varepsilon) = I + \varepsilon[\mathbf n]_\times + O(\varepsilon^2)$, so
$R(\mathbf n_1,\varepsilon)R(\mathbf n_2,\varepsilon) = I + \varepsilon([\mathbf n_1]_\times + [\mathbf n_2]_\times) + O(\varepsilon^2)$, the same in either order. Small rotations add like vectors; that is why rates of rotation do. Precisely, the **addition rule**: if $R = R_AR_B$ then
$$\dot RR^\top = \dot R_AR_A^\top + R_A(\dot R_BR_B^\top)R_A^\top = [\boldsymbol\omega_A]_\times + [R_A\boldsymbol\omega_B]_\times,$$
so the total angular velocity is $\boldsymbol\omega_A$ plus $\boldsymbol\omega_B$ carried into the space frame by $R_A$. §3.11 runs this rule three times.

### 3.3 Kinetic energy, angular momentum, the inertia tensor

Take the reference point $O$ to be the CM, or a pivot fixed in space. Write $\mathbf r_a$ for particle $a$'s position *relative to $O$*, and $\mathbf V$ for the velocity of $O$. By (3.4), $\mathbf v_a = \mathbf V + \boldsymbol\omega\times\mathbf r_a$.

▶ **Try it first.** Expand $T = \frac12\sum m_a|\mathbf V + \boldsymbol\omega\times\mathbf r_a|^2$ and kill the cross term. Then expand $|\boldsymbol\omega\times\mathbf r|^2$ as a quadratic form in $\boldsymbol\omega$.

*Derivation of $T$.*

1. $T = \frac12M|\mathbf V|^2 + \mathbf V\cdot\big(\boldsymbol\omega\times\sum_am_a\mathbf r_a\big) + \frac12\sum_am_a|\boldsymbol\omega\times\mathbf r_a|^2$.
2. The middle term vanishes: either $O$ is the CM, so $\sum m_a\mathbf r_a = \mathbf 0$, or $O$ is a fixed pivot, so $\mathbf V = \mathbf 0$. So $T = \frac12M|\mathbf V|^2 + T_{\rm rot}$ with
$$T_{\rm rot} = \tfrac12\sum_am_a|\boldsymbol\omega\times\mathbf r_a|^2. \tag{3.5}$$
3. Lagrange's identity $|\boldsymbol\omega\times\mathbf r|^2 = |\boldsymbol\omega|^2|\mathbf r|^2 - (\boldsymbol\omega\cdot\mathbf r)^2 = \sum_{i,j}\omega_i\big(r^2\delta_{ij} - x_ix_j\big)\omega_j$, with $x_i$ the components of $\mathbf r$. So
$$T_{\rm rot} = \tfrac12\,\boldsymbol\omega^\top I\boldsymbol\omega, \tag{3.6}$$
with $I$ given by (3.8) below.

*Derivation of $\mathbf L$.* About $O$,
$$\mathbf L = \sum_am_a\mathbf r_a\times(\boldsymbol\omega\times\mathbf r_a). \tag{3.7}$$
(For $O$ = CM this is the spin angular momentum; the orbital part $M\mathbf R_O\times\mathbf V$ separates off as in CM-1.) BAC-CAB: $\mathbf r\times(\boldsymbol\omega\times\mathbf r) = \boldsymbol\omega(\mathbf r\cdot\mathbf r) - \mathbf r(\mathbf r\cdot\boldsymbol\omega)$, whose $i$-th component is $\sum_j(r^2\delta_{ij} - x_ix_j)\omega_j$. So $\mathbf L = I\boldsymbol\omega$ with the **same** matrix
$$I_{ij} = \sum_am_a\big(r_a^2\,\delta_{ij} - x_{a,i}x_{a,j}\big), \qquad I_{ij} = \int\rho(\mathbf r)\,\big(r^2\delta_{ij} - x_ix_j\big)\,d^3r \ \text{ for a continuum of density }\rho. \tag{3.8}$$ ∎

■ **What just happened.** One matrix does both jobs: $\mathbf L = I\boldsymbol\omega$ and $T_{\rm rot} = \frac12\boldsymbol\omega^\top I\boldsymbol\omega = \frac12\boldsymbol\omega\cdot\mathbf L$. The diagonal entries $I_{xx} = \sum m(y^2 + z^2)$ are the **moments of inertia** about the coordinate axes (mass times squared distance from the axis). The off-diagonal entries $I_{xy} = -\sum m\,xy$ are the **products of inertia**, with the minus sign built in. (Some ME texts define products without the minus and put it in the matrix; check before copying a table.) In code, (3.8) is one `einsum` per term: `'a,ai,aj->ij'` is literally $\sum_am_ax_{a,i}x_{a,j}$ (LA-12 Example 4.6).

### 3.4 Symmetric, positive definite, and a tensor

**Theorem 3.2.** (a) $I = I^\top$. (b) $\boldsymbol\omega^\top I\boldsymbol\omega \ge 0$, with equality for $\boldsymbol\omega\ne\mathbf 0$ only if all the mass lies on the line through $O$ along $\boldsymbol\omega$; so $I$ is positive definite for every body that is not a line. (c) If the body is rotated by $R$ (every $\mathbf r_a\mapsto R\mathbf r_a$), its inertia tensor becomes
$$I' = RIR^\top. \tag{3.9}$$

*Proof.* (a) (3.8) is symmetric in $i, j$. (b) By (3.5)–(3.6), $\boldsymbol\omega^\top I\boldsymbol\omega = \sum m_a|\boldsymbol\omega\times\mathbf r_a|^2 \ge 0$, and it is zero iff every $\mathbf r_a\parallel\boldsymbol\omega$. (c) $|R\mathbf r|^2 = |\mathbf r|^2$ and $(R\mathbf r)(R\mathbf r)^\top = R\mathbf r\mathbf r^\top R^\top$, so $\sum m(r^2I - R\mathbf r\mathbf r^\top R^\top) = R\big[\sum m(r^2I - \mathbf r\mathbf r^\top)\big]R^\top$. ∎

(3.9) is the active form of LA-12's $(0,2)$ tensor law (3.19), $\sigma' = R^\top\sigma R$ for rotated *axes*. The consequence that drives everything below: as the body tumbles, its inertia tensor in space components is $I_{\rm space}(t) = R(t)\,I_{\rm body}\,R(t)^\top$, which changes in time; in body components it is **constant**. That is why we will write the equations of motion in the body frame.

### 3.5 Principal axes, the inertia ellipsoid, the triangle inequalities

**Theorem 3.3 (Principal axes).** For every rigid body and every reference point there is an orthonormal, right-handed body frame $\hat{\mathbf e}_1, \hat{\mathbf e}_2, \hat{\mathbf e}_3$ in which $I = \operatorname{diag}(I_1, I_2, I_3)$ with all $I_k \ge 0$. In that frame
$$\mathbf L = (I_1\omega_1,\ I_2\omega_2,\ I_3\omega_3), \qquad T_{\rm rot} = \tfrac12\big(I_1\omega_1^2 + I_2\omega_2^2 + I_3\omega_3^2\big), \tag{3.10}$$
and $\mathbf L\parallel\boldsymbol\omega$ exactly when $\boldsymbol\omega$ is an eigenvector of $I$.

*Proof.* $I$ is real symmetric, so LA-8 Theorem 3.11 gives $I = Q\Lambda Q^\top$ with $Q$ orthogonal and $\Lambda$ real diagonal; flip the sign of one column if needed to make $\det Q = +1$. The columns of $Q$ are the $\hat{\mathbf e}_k$, and $I_k \ge 0$ by Theorem 3.2(b). $I\boldsymbol\omega = \lambda\boldsymbol\omega$ is the definition of an eigenvector. ∎

So the guess of §2 is right: however irregular the body, three perpendicular axes exist about which it can spin steadily with no torque. Two useful facts for finding them without computation. A **mirror plane** $z\to -z$ makes $I_{xz} = I_{yz} = 0$ (the integrand is odd), so the normal to a symmetry plane is principal. An **$n$-fold symmetry axis** with $n\ge 3$ is principal with the two perpendicular moments equal, because the $2\times2$ block must commute with a rotation by $2\pi/n \ne \pi$, and the only symmetric $2\times2$ matrices that do are multiples of the identity. A body with $I_1 = I_2 \ne I_3$ is a **symmetric top**; with all three equal (sphere, cube about its center) a **spherical top**, for which every axis is principal.

**The inertia ellipsoid.** The surface $\boldsymbol\omega^\top I\boldsymbol\omega = 1$ is, in principal axes, $I_1\omega_1^2 + I_2\omega_2^2 + I_3\omega_3^2 = 1$: an ellipsoid with semi-axes $1/\sqrt{I_k}$ along the principal axes. Since
$$\nabla_{\boldsymbol\omega}\big(\tfrac12\boldsymbol\omega^\top I\boldsymbol\omega\big) = I\boldsymbol\omega = \mathbf L, \tag{3.11}$$
**$\mathbf L$ is normal to the ellipsoid at the point where $\boldsymbol\omega$ pierces it.** The normal of an ellipsoid is along the radius only at the ends of its axes. That is Theorem 3.3 as a picture, and it is the picture Poinsot used (§3.10).

**Theorem 3.4 (Triangle inequalities).** The principal moments about any point satisfy
$$I_1 + I_2 \ge I_3,\quad I_2 + I_3 \ge I_1,\quad I_3 + I_1 \ge I_2, \tag{3.12}$$
with equality in the first iff the body is a lamina in the $\hat{\mathbf e}_1\hat{\mathbf e}_2$ plane.

*Derivation.*

1. In principal axes the diagonal entries of (3.8) are the principal moments: $I_1 = \sum m(x_2^2 + x_3^2)$, $I_2 = \sum m(x_1^2 + x_3^2)$, $I_3 = \sum m(x_1^2 + x_2^2)$, with $x_k = \mathbf r\cdot\hat{\mathbf e}_k$.
2. $I_1 + I_2 - I_3 = 2\sum m\,x_3^2 \ge 0$. The other two are the same with indices permuted.
3. Equality iff every $x_3 = 0$. ∎

■ **What just happened.** Not every triple of positive numbers is a possible set of moments; $(1, 2, 4)$ is not. And the triangle inequality's equality case is a theorem you already use: §3.6's perpendicular-axis theorem. A bonus from LA-8 Theorem 3.17 (Rayleigh quotient): the moment about any axis $\hat{\mathbf n}$, $I_{\hat{\mathbf n}} = \hat{\mathbf n}^\top I\hat{\mathbf n}$, lies between $I_{\min}$ and $I_{\max}$.

### 3.6 The parallel-axis and perpendicular-axis theorems

▶ **Try it first.** Put $\mathbf r_a = \mathbf d + \mathbf s_a$, where $\mathbf d$ points from $O$ to the CM and $\sum m_a\mathbf s_a = \mathbf 0$. Expand $r^2\delta_{ij} - x_ix_j$ and sum.

**Theorem 3.5 (Parallel-axis).** With $M$ the total mass and $\mathbf d$ the position of the CM relative to $O$,
$$I^{(O)}_{ij} = I^{\rm (cm)}_{ij} + M\big(d^2\delta_{ij} - d_id_j\big). \tag{3.13}$$

*Derivation.*

1. $r_a^2 = d^2 + 2\mathbf d\cdot\mathbf s_a + s_a^2$ and $x_{a,i}x_{a,j} = d_id_j + d_is_{a,j} + s_{a,i}d_j + s_{a,i}s_{a,j}$.
2. Multiply by $m_a$ and sum. Every term linear in $\mathbf s_a$ is killed by $\sum m_a\mathbf s_a = \mathbf 0$.
3. What remains is $\sum m_a(s_a^2\delta_{ij} - s_{a,i}s_{a,j}) + M(d^2\delta_{ij} - d_id_j)$. The first sum is $I^{\rm (cm)}$. ∎

■ **What just happened.** The CM is the one point about which the cross terms die, so moving away from it adds a point-mass tensor $M(d^2\delta_{ij} - d_id_j)$ and nothing else. The scalar version you know from statics is the diagonal: about any axis, $I = I_{\rm cm} + Md_\perp^2$ with $d_\perp$ the distance between the two parallel axes. The tensor version also shifts the products of inertia, which is what makes Example 4.1 work. It only goes **from the CM**; to move between two non-CM points, go through the CM.

**Perpendicular-axis theorem.** For a lamina in the $xy$ plane ($z = 0$ for every particle), about any origin in its plane,
$$I_{zz} = I_{xx} + I_{yy}, \tag{3.14}$$
because $I_{xx} = \sum my^2$, $I_{yy} = \sum mx^2$ and $I_{zz} = \sum m(x^2 + y^2)$. This is the equality case of (3.12).

### 3.7 A short table, with the derivations

| Body (mass $M$) | Principal moments |
|---|---|
| Thin rod, length $\ell$ | $\frac1{12}M\ell^2$ about center, $\frac13M\ell^2$ about an end (perpendicular axes) |
| Rectangular plate $a\times b$ in the $xy$ plane, about its center | $I_{xx} = \frac1{12}Mb^2$, $I_{yy} = \frac1{12}Ma^2$, $I_{zz} = \frac1{12}M(a^2 + b^2)$ |
| Solid sphere, radius $R$ | $\frac25MR^2$ (spherical top) |
| Solid cube, side $a$, about its center | $\frac16Ma^2$ (spherical top) |
| Solid cone, height $h$, base radius $R$, about its apex | $I_3 = \frac3{10}MR^2$, $I_1 = I_2 = \frac3{20}MR^2 + \frac35Mh^2$ |
| Same cone about its CM ($\frac34h$ from apex) | $I_1 = I_2 = \frac3{20}MR^2 + \frac3{80}Mh^2$ |

*Rod:* $\int_{-\ell/2}^{\ell/2}(M/\ell)x^2\,dx = M\ell^2/12$; parallel axis with $d = \ell/2$ gives $M\ell^2/12 + M\ell^2/4 = M\ell^2/3$. *Plate:* slice into strips parallel to $y$; each is a rod of length $b$ whose midpoint lies on the $x$-axis, so each contributes (its mass)$\cdot b^2/12$ about that axis and $I_{xx} = Mb^2/12$; likewise $I_{yy} = Ma^2/12$; then (3.14). *Sphere:* by symmetry $I_{xx} = I_{yy} = I_{zz} =: I_s$, and their sum is $2\int\rho r^2\,dV$, so $3I_s = 2\cdot\frac{3M}{4\pi R^3}\int_0^R r^2\cdot4\pi r^2\,dr = \frac65MR^2$ and $I_s = \frac25MR^2$. *Cube:* $I_{xx} = \sum m(y^2 + z^2) = 2\cdot\frac{Ma^2}{12}$, the products vanish by mirror symmetry. *Cone:* the slice at height $z$ above the apex is a disk of radius $Rz/h$ and mass $dm = 3Mz^2\,dz/h^3$. A disk has $\frac12mr^2$ about its axis and, by (3.14), $\frac14mr^2$ about a diameter. So $I_3 = \int_0^h\frac12\frac{R^2z^2}{h^2}\frac{3Mz^2}{h^3}dz = \frac3{10}MR^2$ and $I_1 = \int_0^h\big(\frac14\frac{R^2z^2}{h^2} + z^2\big)\frac{3Mz^2}{h^3}dz = \frac35M\big(\frac{R^2}4 + h^2\big)$. The CM is at $\bar z = \int z\,dm/M = \frac34h$, and (3.13) gives $I_1^{\rm cm} = \frac3{20}MR^2 + \frac35Mh^2 - \frac9{16}Mh^2 = \frac3{20}MR^2 + \frac3{80}Mh^2$. The project Monte Carlo samples each body (rod, plate, sphere, cube, cone) and reproduces every entry to within 1%.

### 3.8 Euler's equations

The law is $d\mathbf L/dt = \boldsymbol\tau$ (about the CM or a fixed pivot; CM-1). The difficulty: in the space frame $I$ changes as the body turns; in the body frame $I$ is constant but the frame rotates. We take the second deal and pay for the rotating frame.

▶ **Try it first.** Write $\mathbf A = \sum_kA_k\hat{\mathbf e}_k$ in body axes and differentiate, using (3.4) on each $\hat{\mathbf e}_k$.

*Derivation.*

1. $\dfrac{d\mathbf A}{dt} = \sum_k\dot A_k\hat{\mathbf e}_k + \sum_kA_k\dot{\hat{\mathbf e}}_k$.
2. Each $\hat{\mathbf e}_k$ is body-fixed, so by (3.4) $\dot{\hat{\mathbf e}}_k = \boldsymbol\omega\times\hat{\mathbf e}_k$, and the second sum is $\boldsymbol\omega\times\mathbf A$. Writing $(d\mathbf A/dt)_{\rm body} := \sum_k\dot A_k\hat{\mathbf e}_k$ (the rate an observer riding the body sees),
$$\left(\frac{d\mathbf A}{dt}\right)_{\rm space} = \left(\frac{d\mathbf A}{dt}\right)_{\rm body} + \boldsymbol\omega\times\mathbf A. \tag{3.15}$$
3. Apply (3.15) to $\mathbf A = \mathbf L = I\boldsymbol\omega$ with $I$ constant in body components:
$$I\dot{\boldsymbol\omega} + \boldsymbol\omega\times(I\boldsymbol\omega) = \boldsymbol\tau. \tag{3.16}$$
4. In principal axes $\mathbf L = (I_1\omega_1, I_2\omega_2, I_3\omega_3)$ and $(\boldsymbol\omega\times\mathbf L)_1 = \omega_2L_3 - \omega_3L_2 = (I_3 - I_2)\omega_2\omega_3$, and cyclically. So:

**Theorem 3.6 (Euler's equations).** In principal body axes,
$$I_1\dot\omega_1 - (I_2 - I_3)\,\omega_2\omega_3 = \tau_1,\qquad I_2\dot\omega_2 - (I_3 - I_1)\,\omega_3\omega_1 = \tau_2,\qquad I_3\dot\omega_3 - (I_1 - I_2)\,\omega_1\omega_2 = \tau_3. \tag{3.17}$$ ∎

■ **What just happened.** Three first-order nonlinear ODEs for the body components of $\boldsymbol\omega$. Every symbol is a body component, including $\boldsymbol\tau$; the derivative on the left is still the space-frame derivative of $\mathbf L$, written in body components. The quadratic terms are the price of the rotating frame, and they are the source of everything surprising below. With $\boldsymbol\tau = \mathbf 0$ (free body), two quantities are conserved: $\dot T = \boldsymbol\omega\cdot I\dot{\boldsymbol\omega} = -\boldsymbol\omega\cdot(\boldsymbol\omega\times\mathbf L) = 0$, and since in body components $\dot{\mathbf L}_{\rm body} = -\boldsymbol\omega\times\mathbf L$, $\frac{d}{dt}|\mathbf L|^2 = -2\mathbf L\cdot(\boldsymbol\omega\times\mathbf L) = 0$. Knowing $\boldsymbol\omega_{\rm body}(t)$ is not knowing the orientation; for that you also integrate $\dot R = R[\boldsymbol\omega_{\rm body}]_\times$ (the project carries both).

### 3.9 The free symmetric top, solved exactly

Let $I_1 = I_2 \ne I_3$, $\boldsymbol\tau = \mathbf 0$.

▶ **Try it first.** Put $I_1 = I_2$ into (3.17). Show $\omega_3$ is constant, then write the first two equations as one equation for $\omega_1 + i\omega_2$ and solve it.

1. Third equation of (3.17): $I_3\dot\omega_3 = (I_1 - I_2)\omega_1\omega_2 = 0$, so $\omega_3$ is constant.
2. First two: $\dot\omega_1 = -\Omega_b\,\omega_2$, $\dot\omega_2 = \Omega_b\,\omega_1$ with
$$\Omega_b = \frac{I_3 - I_1}{I_1}\,\omega_3. \tag{3.18}$$
3. So $\omega_1 + i\omega_2 = \omega_\perp e^{i(\Omega_bt + \alpha)}$: in the body, the perpendicular part of $\boldsymbol\omega$ goes around $\hat{\mathbf e}_3$ at rate $\Omega_b$ with constant length. $\boldsymbol\omega$ sweeps the **body cone** about $\hat{\mathbf e}_3$.
4. In space, $\mathbf L$ is constant. Write $\mathbf L = I_1\boldsymbol\omega_\perp + I_3\omega_3\hat{\mathbf e}_3 = I_1\boldsymbol\omega + (I_3 - I_1)\omega_3\hat{\mathbf e}_3$. So $\boldsymbol\omega$, $\mathbf L$, $\hat{\mathbf e}_3$ are coplanar, and $\boldsymbol\omega = \mathbf L/I_1 - (\ldots)\hat{\mathbf e}_3$.
5. By (3.4), $\dot{\hat{\mathbf e}}_3 = \boldsymbol\omega\times\hat{\mathbf e}_3 = (\mathbf L/I_1)\times\hat{\mathbf e}_3$. The symmetry axis precesses about the fixed $\mathbf L$ at
$$\Omega_s = \frac{|\mathbf L|}{I_1}, \tag{3.19}$$
at fixed angle $\theta_L$ with $\tan\theta_L = I_1\omega_\perp/(I_3\omega_3)$. Since $\boldsymbol\omega$ is a fixed combination of $\mathbf L$ and $\hat{\mathbf e}_3$, it precesses with it, sweeping the **space cone** about $\mathbf L$.

■ **What just happened.** The motion is the body cone rolling without slipping on the space cone; their line of contact is $\boldsymbol\omega$, the instantaneous axis, along which body points are momentarily at rest. For an oblate body ($I_3 > I_1$, a disk, the Earth) $\theta_L$ is smaller than the angle of $\boldsymbol\omega$, so $\mathbf L$ lies between $\hat{\mathbf e}_3$ and $\boldsymbol\omega$ and the space cone sits inside the body cone. For a prolate body (a rugby ball, a rod) $\boldsymbol\omega$ lies between them and the body cone rolls on the outside of the space cone (Exercise B3). Two rates are in play and they are different: $\Omega_b$ is what an observer on the body sees; $\Omega_s$ is the wobble an outside observer sees. In Euler angles with $\hat{\mathbf z}\parallel\mathbf L$ (§3.11), $\dot\phi = \Omega_s$ and $\dot\psi = \omega_3 - \dot\phi\cos\theta_L = \omega_3 - I_3\omega_3/I_1 = -\Omega_b$ exactly, since $|\mathbf L|\cos\theta_L = I_3\omega_3$.

**The Chandler wobble, honestly.** The Earth is oblate with $(I_3 - I_1)/I_1 \approx 1/305$. If its rotation axis is slightly off its figure axis, (3.18) says the pole wanders around the figure axis with period $2\pi/\Omega_b = 305$ sidereal days. This is Euler's prediction. Chandler's observed wobble, a pole motion of a few meters, has a period of about **433 days**. The rigid model is wrong by 40%, and the error is physics: the Earth yields elastically and its oceans slosh, which lengthens the period (Newcomb's explanation, 1892). Use the rigid number as an order-of-magnitude statement, not a prediction.

### 3.10 The asymmetric body: stability, and Poinsot

Now $I_1 < I_2 < I_3$, torque-free. Steady spin $\boldsymbol\omega = \omega_0\hat{\mathbf e}_k$ about any principal axis is an exact solution of (3.17). Is it stable?

▶ **Try it first.** Spin about $\hat{\mathbf e}_1$: set $\boldsymbol\omega = (\omega_0 + \delta_1, \delta_2, \delta_3)$ with small $\delta$s, drop products of $\delta$s, and get a second-order equation for $\delta_2$.

*Derivation (spin about $\hat{\mathbf e}_k$; label the other two axes $i, j$).*

1. Linearize (3.17) about $\omega_0\hat{\mathbf e}_1$: $I_1\dot\delta_1 = (I_2 - I_3)\delta_2\delta_3 \approx 0$, $I_2\dot\delta_2 = (I_3 - I_1)\omega_0\delta_3$, $I_3\dot\delta_3 = (I_1 - I_2)\omega_0\delta_2$.
2. Differentiate the second and substitute the third: $\ddot\delta_2 = \dfrac{(I_3 - I_1)(I_1 - I_2)}{I_2I_3}\,\omega_0^2\,\delta_2$.
3. The same steps about any axis $k$ give, for each transverse component,
$$\ddot\delta = \lambda\,\delta, \qquad \lambda = \omega_0^2\,\frac{(I_k - I_i)(I_j - I_k)}{I_iI_j}. \tag{3.20}$$
4. If $I_k$ is the largest or the smallest moment, the two factors have opposite signs, $\lambda < 0$, and $\delta$ oscillates at $\sqrt{-\lambda}$: stable. If $I_k$ is the intermediate moment, $\lambda > 0$ and $\delta$ grows like $e^{\sigma t}$ with $\sigma = \sqrt\lambda$.

**Theorem 3.7 (Tennis-racket theorem).** Free rotation about the axes of largest and smallest moment is (linearly) stable; about the intermediate axis it is unstable, with growth rate
$$\sigma = \omega_0\sqrt{\frac{(I_2 - I_1)(I_3 - I_2)}{I_1I_3}}. \qquad ∎$$

■ **What just happened.** You can predict the phone toss: the flipping axis is the in-plane short axis, whose moment is in the middle. Linear theory only says the deviation grows. What happens next is not chaos; the conservation laws confine the motion, and the body flips over, flips back, and repeats periodically (Dzhanibekov's wing nut; the project shows $\omega_2$ reversing sign).

**Poinsot's picture: the global version.** In body components, $\mathbf L$ must satisfy two conservation laws at once:
$$L_1^2 + L_2^2 + L_3^2 = |\mathbf L|^2 \quad\text{(a sphere)}, \qquad \frac{L_1^2}{I_1} + \frac{L_2^2}{I_2} + \frac{L_3^2}{I_3} = 2T \quad\text{(an ellipsoid)}.$$
So $\mathbf L_{\rm body}(t)$ runs along an intersection curve (a **polhode**). Think of $2T$ as a function on the sphere $|\mathbf L| = $ const. It is largest at $\pm\hat{\mathbf e}_1$ (smallest $I$), smallest at $\pm\hat{\mathbf e}_3$, and has a saddle at $\pm\hat{\mathbf e}_2$. Level curves near a maximum or minimum are small closed loops: stable. Level curves through a saddle are two great-circle-like curves crossing there (the separatrices), and a curve that starts near $\hat{\mathbf e}_2$ runs all the way to near $-\hat{\mathbf e}_2$: the flip. The project draws these curves (`poinsot_sphere.png`). Poinsot's own construction is the space-frame version: by (3.11) the energy ellipsoid $\boldsymbol\omega^\top I\boldsymbol\omega = 2T$ has normal $\mathbf L$ at $\boldsymbol\omega$, and $\boldsymbol\omega\cdot\mathbf L = 2T$ is constant, so the ellipsoid rolls without slipping on a fixed plane perpendicular to $\mathbf L$ (Goldstein §5.6; Landau & Lifshitz §37).

### 3.11 Euler angles

To write a Lagrangian you need three coordinates for $R$. Euler's choice: start with body axes equal to space axes, then (1) rotate by $\phi$ about $\hat{\mathbf z}$; (2) rotate by $\theta$ about the *new* $x$-axis, called the **line of nodes** $\hat{\mathbf n}$; (3) rotate by $\psi$ about the *new* $z$-axis, which is the final $\hat{\mathbf e}_3$.

```
        z (space)    e3 (body)
         |  θ      /
         | <-->   /        θ: tilt of e3 from z
         |       /         φ: azimuth of the line of nodes n in the xy plane
         |      /          ψ: spin of the body about e3, measured from n
         +--------- y
        / \
       /   \ n (line of nodes, angle φ from x)
      x
```

Rotations about successively *new* axes compose by right multiplication: after the first rotation $R_1 = R_z(\phi)$, a rotation by $\theta$ about $\hat{\mathbf n} = R_1\hat{\mathbf x}$ is $R_1R_x(\theta)R_1^\top$, and applying it to $R_1$ gives $R_1R_x(\theta)$; once more for the third step. Hence
$$R(\phi,\theta,\psi) = R_z(\phi)\,R_x(\theta)\,R_z(\psi) = \begin{pmatrix} c_\phi c_\psi - s_\phi c_\theta s_\psi & -c_\phi s_\psi - s_\phi c_\theta c_\psi & s_\phi s_\theta\\ s_\phi c_\psi + c_\phi c_\theta s_\psi & -s_\phi s_\psi + c_\phi c_\theta c_\psi & -c_\phi s_\theta\\ s_\theta s_\psi & s_\theta c_\psi & c_\theta\end{pmatrix}, \tag{3.21}$$
with $c_\theta = \cos\theta$ and so on (LA-4: composition is multiplication). The third column is $\hat{\mathbf e}_3 = (\sin\phi\sin\theta, -\cos\phi\sin\theta, \cos\theta)$: $\theta$ is the tilt of the symmetry axis from vertical and $\phi$ is its azimuth. That is why this convention fits the top.

Now the computation every student dreads. Do it once, carefully, and never again.

▶ **Try it first.** Apply the addition rule of §3.2 to the product of three factors in (3.21). Each factor rotates about a single fixed axis, so each $\dot R_kR_k^\top$ is a rate times a hat matrix.

**Theorem 3.8 (Angular velocity in Euler angles).** For $R = R_z(\phi)R_x(\theta)R_z(\psi)$, the angular velocity is $\boldsymbol\omega = \dot\phi\,\hat{\mathbf z} + \dot\theta\,\hat{\mathbf n} + \dot\psi\,\hat{\mathbf e}_3$, with $\hat{\mathbf n} = R_z(\phi)\hat{\mathbf x}$ the line of nodes; its body components are (3.23) and its space components (3.24).

*Derivation.*

1. A rotation about a fixed axis $\hat{\mathbf a}$ by angle $\alpha(t)$ has $\dot R R^\top = \dot\alpha[\hat{\mathbf a}]_\times$ (differentiate Rodrigues' formula, or note that Theorem 3.1's $\boldsymbol\omega$ is $\dot\alpha\hat{\mathbf a}$). So $\dot R_z(\phi)R_z(\phi)^\top = \dot\phi[\hat{\mathbf z}]_\times$, $\dot R_x(\theta)R_x(\theta)^\top = \dot\theta[\hat{\mathbf x}]_\times$, $\dot R_z(\psi)R_z(\psi)^\top = \dot\psi[\hat{\mathbf z}]_\times$.
2. Apply the addition rule twice to $R = R_z(\phi)\cdot R_x(\theta)\cdot R_z(\psi)$:
$\dot RR^\top = \dot\phi[\hat{\mathbf z}]_\times + R_z(\phi)\,\dot\theta[\hat{\mathbf x}]_\times R_z(\phi)^\top + R_z(\phi)R_x(\theta)\,\dot\psi[\hat{\mathbf z}]_\times\,(R_z(\phi)R_x(\theta))^\top.$
3. Use $R[\mathbf a]_\times R^\top = [R\mathbf a]_\times$: the second term is $\dot\theta[\hat{\mathbf n}]_\times$ with $\hat{\mathbf n} = R_z(\phi)\hat{\mathbf x}$, the third is $\dot\psi[\hat{\mathbf e}_3]_\times$ with $\hat{\mathbf e}_3 = R_z(\phi)R_x(\theta)\hat{\mathbf z}$ (the last $R_z(\psi)$ does not move $\hat{\mathbf z}$). So
$$\boldsymbol\omega = \dot\phi\,\hat{\mathbf z} + \dot\theta\,\hat{\mathbf n} + \dot\psi\,\hat{\mathbf e}_3. \tag{3.22}$$
4. **Body components.** We need the three unit vectors in body components, $R^\top\mathbf a$. $\hat{\mathbf e}_3\to(0, 0, 1)$. $\hat{\mathbf z}\to$ the third *row* of (3.21), $(s_\theta s_\psi, s_\theta c_\psi, c_\theta)$. $\hat{\mathbf n}\to R_z(\psi)^\top R_x(\theta)^\top R_z(\phi)^\top R_z(\phi)\hat{\mathbf x} = R_z(-\psi)R_x(-\theta)\hat{\mathbf x} = R_z(-\psi)\hat{\mathbf x} = (\cos\psi, -\sin\psi, 0)$. Adding with the weights in (3.22):
$$\omega_1 = \dot\phi\sin\theta\sin\psi + \dot\theta\cos\psi,\qquad \omega_2 = \dot\phi\sin\theta\cos\psi - \dot\theta\sin\psi,\qquad \omega_3 = \dot\phi\cos\theta + \dot\psi. \tag{3.23}$$
5. **Space components.** $\hat{\mathbf z} = (0,0,1)$, $\hat{\mathbf n} = (\cos\phi, \sin\phi, 0)$, $\hat{\mathbf e}_3$ = third column of (3.21):
$$\omega_x = \dot\theta\cos\phi + \dot\psi\sin\theta\sin\phi,\qquad \omega_y = \dot\theta\sin\phi - \dot\psi\sin\theta\cos\phi,\qquad \omega_z = \dot\phi + \dot\psi\cos\theta. \tag{3.24}$$ ∎

■ **What just happened.** $\boldsymbol\omega$ is the sum of three angular velocities about three non-orthogonal axes; (3.23) and (3.24) are that sum resolved two ways. Checks you should run: $\hat{\mathbf z}\cdot\hat{\mathbf e}_3 = \cos\theta$, $\hat{\mathbf n}\perp$ both, so $|\boldsymbol\omega|^2 = \dot\phi^2 + \dot\theta^2 + \dot\psi^2 + 2\dot\phi\dot\psi\cos\theta$ from either form; and $R\,\boldsymbol\omega_{\rm body} = \boldsymbol\omega_{\rm space}$. The project checks (3.23) and (3.24) against $\boldsymbol\omega$ extracted from finite differences of $R(t)$ along random trajectories.

**Gimbal lock.** Write (3.23) as $\boldsymbol\omega_{\rm body} = J(\theta,\psi)\,(\dot\phi, \dot\theta, \dot\psi)^\top$ with
$J = \begin{pmatrix}s_\theta s_\psi & c_\psi & 0\\ s_\theta c_\psi & -s_\psi & 0\\ c_\theta & 0 & 1\end{pmatrix}$. Expanding along the third column, $\det J = s_\theta s_\psi(-s_\psi) - c_\psi s_\theta c_\psi = -\sin\theta$. At $\theta = 0$ or $\pi$ the map is singular. Geometrically, $\hat{\mathbf e}_3 = \pm\hat{\mathbf z}$, the first and third rotations are about the same axis, and (3.21) depends only on $\phi\pm\psi$: one degree of freedom has vanished from the coordinates, though not from the body. Near $\theta = 0$ a smooth tumble needs $\dot\phi$ and $\dot\psi$ of order $1/\theta$, and the condition number of $J$ grows like $1/\theta$ (project, `gimbal_lock_conditioning.png`). Physical three-gimbal platforms have the same disease; the Apollo guidance system had a warning for approaching gimbal lock. Every three-parameter description of $SO(3)$ is singular somewhere (a topological fact we do not prove). In software, integrate $R$ itself (and re-orthonormalize) or a unit quaternion, which is an element of $SU(2)$ and covers $SO(3)$ twice with no singularity (LA-12 A.5–A.6). Use Euler angles when a symmetry makes two of them cyclic, which is exactly the top.

**The symmetric top's kinetic energy.** With $I_1 = I_2$ (about the CM or the pivot),
$\omega_1^2 + \omega_2^2 = (\dot\phi s_\theta s_\psi + \dot\theta c_\psi)^2 + (\dot\phi s_\theta c_\psi - \dot\theta s_\psi)^2 = \dot\phi^2\sin^2\theta + \dot\theta^2$ (the cross terms cancel and $\psi$ drops out), so
$$T = \tfrac12I_1\big(\dot\theta^2 + \dot\phi^2\sin^2\theta\big) + \tfrac12I_3\big(\dot\psi + \dot\phi\cos\theta\big)^2. \tag{3.25}$$

### 3.12 The heavy symmetric top

A symmetric top spins on a fixed frictionless pivot on its symmetry axis, with its CM a distance $\ell$ up the axis. $I_1 = I_2$ and $I_3$ are taken **about the pivot** ($I_1 = I_1^{\rm cm} + M\ell^2$ by (3.13)). The CM is at height $\ell\cos\theta$. The Lagrangian is
$$\mathcal L = \tfrac12I_1\big(\dot\theta^2 + \dot\phi^2\sin^2\theta\big) + \tfrac12I_3\big(\dot\psi + \dot\phi\cos\theta\big)^2 - Mg\ell\cos\theta. \tag{3.26}$$

**Two cyclic coordinates.** $\phi$ and $\psi$ do not appear, so by CM-4 Theorem 3.1
$$p_\psi = \frac{\partial\mathcal L}{\partial\dot\psi} = I_3(\dot\psi + \dot\phi\cos\theta) = I_3\omega_3, \qquad p_\phi = \frac{\partial\mathcal L}{\partial\dot\phi} = I_1\dot\phi\sin^2\theta + p_\psi\cos\theta \tag{3.27}$$
are conserved. Physically, $p_\psi = L_3$ (spin along the axis) and $p_\phi = L_z$ (the vertical component); gravity's torque $\ell\hat{\mathbf e}_3\times(-Mg\hat{\mathbf z})$ is perpendicular to both $\hat{\mathbf z}$ and $\hat{\mathbf e}_3$. $\mathcal L$ has no explicit $t$, so the energy function $h$ is conserved (CM-4 Theorem 3.6), and $h = E = T + V$ because the pivot is fixed (CM-4 Theorem 3.7). Three constants for three coordinates: the problem is integrable, which is why Lagrange could solve it.

**The $\theta$ equation.** $\frac{d}{dt}(I_1\dot\theta) = \partial\mathcal L/\partial\theta$:
$$I_1\ddot\theta = I_1\dot\phi^2\sin\theta\cos\theta - p_\psi\dot\phi\sin\theta + Mg\ell\sin\theta. \tag{3.28}$$

▶ **Try it first.** Solve (3.27) for $\dot\phi$ and $\dot\psi$ in terms of $\theta$ and the two momenta. Substitute into $E$, not into $\mathcal L$, and collect what depends on $\theta$ alone.

*Derivation of the effective potential.*

1. From (3.27): $\dot\phi = \dfrac{p_\phi - p_\psi\cos\theta}{I_1\sin^2\theta}$ and $\dot\psi = \dfrac{p_\psi}{I_3} - \dot\phi\cos\theta$.
2. $E = \frac12I_1\dot\theta^2 + \frac12I_1\dot\phi^2\sin^2\theta + \dfrac{p_\psi^2}{2I_3} + Mg\ell\cos\theta$, using $I_3(\dot\psi + \dot\phi\cos\theta) = p_\psi$.
3. Substitute step 1 into the second term: $\frac12I_1\sin^2\theta\cdot\dfrac{(p_\phi - p_\psi\cos\theta)^2}{I_1^2\sin^4\theta} = \dfrac{(p_\phi - p_\psi\cos\theta)^2}{2I_1\sin^2\theta}$.
4. Move the constant $p_\psi^2/2I_3$ to the left. With $E' := E - p_\psi^2/(2I_3)$,
$$E' = \tfrac12I_1\dot\theta^2 + V_{\rm eff}(\theta), \qquad V_{\rm eff}(\theta) = \frac{(p_\phi - p_\psi\cos\theta)^2}{2I_1\sin^2\theta} + Mg\ell\cos\theta. \tag{3.29}$$ ∎

■ **What just happened.** The top is now a particle of "mass" $I_1$ in a one-dimensional potential. $V_{\rm eff}\to+\infty$ at $\theta = 0, \pi$ unless $p_\phi = \pm p_\psi$, so $\theta$ oscillates between two turning points: **nutation**. Why $E$ and not $\mathcal L$: substituting conserved momenta into $\mathcal L$ and then varying flips the sign of the centrifugal term (the same trap as CM-5's radial equation; the correct reduced Lagrangian is Routh's, which CM-5 and CM-8 §3.4 point to). The project confirms (3.29) symbolically.

**Turning points as roots of a cubic.** Put $u = \cos\theta$, so $\dot u^2 = \sin^2\theta\,\dot\theta^2$. Multiply (3.29) by $2(1 - u^2)/I_1$:
$$\dot u^2 = f(u) := (1 - u^2)(\alpha - \beta u) - (b - au)^2, \qquad \alpha = \frac{2E'}{I_1},\ \beta = \frac{2Mg\ell}{I_1},\ a = \frac{p_\psi}{I_1},\ b = \frac{p_\phi}{I_1}. \tag{3.30}$$
$f(\pm1) = -(b\mp a)^2 \le 0$, $f\to+\infty$ as $u\to+\infty$, and real motion needs $f \ge 0$ somewhere in $[-1,1]$. So $f$ has two roots $u_1\le u_2$ in $[-1, 1]$ bracketing the motion and a third root above 1. Meanwhile $\dot\phi = (b - au)/(1 - u^2)$. If $b - au$ keeps one sign on $[u_1, u_2]$ the axis traces a wavy curve; if it vanishes at a turning point, the curve has **cusps** (the top released from rest at its upper turning point); if it changes sign, **loops**. The project draws a wavy path (a kicked steady precession, $\dot\phi$ between 1.90 and 3.74 rad/s) and a cusped one on the unit sphere (`heavy_top_axis_path.png`).

**Steady precession.** $\theta = \theta_0$ constant requires $V_{\rm eff}'(\theta_0) = 0$. Differentiating (3.29) and writing the result in terms of $\dot\phi$ from step 1 gives $V_{\rm eff}' = \sin\theta\,\big[p_\psi\dot\phi - I_1\dot\phi^2\cos\theta - Mg\ell\big]$ (project, SymPy). So
$$I_1\cos\theta_0\,\dot\phi^2 - p_\psi\,\dot\phi + Mg\ell = 0, \qquad \dot\phi_\pm = \frac{p_\psi \pm\sqrt{p_\psi^2 - 4I_1Mg\ell\cos\theta_0}}{2I_1\cos\theta_0}. \tag{3.31}$$
For $\theta_0 < 90°$ steady precession needs $p_\psi^2 \ge 4I_1Mg\ell\cos\theta_0$: a minimum spin. For a **fast top**, $p_\psi^2 \gg 4I_1Mg\ell\cos\theta_0$, expand the root: $\dot\phi_-\approx Mg\ell/p_\psi = Mg\ell/(I_3\omega_3)$ (the slow precession you guessed in §2) and $\dot\phi_+\approx p_\psi/(I_1\cos\theta_0)$ (fast precession, essentially a free top's wobble). At $\theta_0 = 90°$ the quadratic degenerates and $\dot\phi = Mg\ell/p_\psi$ exactly.

**Small nutation.** Near a minimum $\theta_0$ of $V_{\rm eff}$, $\theta$ oscillates at $\omega_{\rm nut} = \sqrt{V_{\rm eff}''(\theta_0)/I_1}$. Differentiate $V_{\rm eff}' = \sin\theta\cdot[\ldots]$ once more; the term with the bracket undifferentiated vanishes at $\theta_0$. Using $d\dot\phi/d\theta = (p_\psi - 2I_1\dot\phi\cos\theta)/(I_1\sin\theta)$ from step 1,
$$\omega_{\rm nut}^2 = \frac{(p_\psi - 2I_1\dot\phi_0\cos\theta_0)^2 + I_1^2\dot\phi_0^2\sin^2\theta_0}{I_1^2} \;\approx\; \Big(\frac{I_3\omega_3}{I_1}\Big)^2 \ \text{(fast top, slow precession)}. \tag{3.32}$$
A fast top nutates fast and with small amplitude, so the eye sees smooth slow precession.

**The sleeping top.** A top spinning upright, $\theta = 0$, has $p_\phi = p_\psi$ by (3.27). Then $V_{\rm eff} = \dfrac{p_\psi^2}{2I_1}\dfrac{(1 - \cos\theta)^2}{\sin^2\theta} + Mg\ell\cos\theta = \dfrac{p_\psi^2}{2I_1}\tan^2\dfrac\theta2 + Mg\ell\cos\theta$. For small $\theta$, $\tan^2(\theta/2)\approx\theta^2/4$ and $\cos\theta\approx1 - \theta^2/2$:
$V_{\rm eff}\approx Mg\ell + \Big(\dfrac{p_\psi^2}{8I_1} - \dfrac{Mg\ell}{2}\Big)\theta^2$. Upright is stable iff the bracket is positive:
$$\omega_3^2 > \frac{4Mg\ell\,I_1}{I_3^2}. \tag{3.33}$$
As friction slows a sleeping top through this threshold, it "wakes up" and begins to wobble.

### 3.13 Gyroscopes, the equinoxes, and rolling

**Gyroscopes.** For a fast rotor, $\mathbf L\approx I_3\omega_3\hat{\mathbf e}_3$ and $d\mathbf L/dt = \boldsymbol\tau$ turns $\mathbf L$ toward $\boldsymbol\tau$: the axis moves perpendicular to the push. The engineering corollary runs the other way: if the casing forces a rotor's axis to turn at $\boldsymbol\Omega$, the bearings must supply $\boldsymbol\tau = \boldsymbol\Omega\times\mathbf L$ (a ship's turbine in a turning hull, a car wheel in a corner, a drone's rotors during a roll). A **gyrocompass** is a rotor whose axis is held horizontal but free to turn about the vertical; the Earth's rotation then produces a torque that swings the axis toward north (Anschütz-Kaempfe, 1908; Exercise C3 derives the restoring equation).

**Precession of the equinoxes.** The Earth is oblate, $H := (I_3 - I_1)/I_3\approx 1/305$, its axis is tilted $\varepsilon = 23.44°$ to the ecliptic, and the Sun and Moon pull harder on the near side of the equatorial bulge. Torque on a body from a point mass $GM$ at $\mathbf r$, origin at the body's CM:

1. Force on element $a$: $\mathbf F_a = GMm_a(\mathbf r - \mathbf r_a)/|\mathbf r - \mathbf r_a|^3$. To first order in $r_a/r$: $(\mathbf r - \mathbf r_a)/|\mathbf r - \mathbf r_a|^3 \approx \mathbf r/r^3 - \mathbf r_a/r^3 + 3\mathbf r(\hat{\mathbf r}\cdot\mathbf r_a)/r^4$.
2. $\boldsymbol\tau = \sum\mathbf r_a\times\mathbf F_a$. The $\mathbf r_a\times\mathbf r_a$ term is zero and $\sum m_a\mathbf r_a\times\mathbf r = 0$ (CM). So $\boldsymbol\tau = \frac{3GM}{r^3}\sum m_a(\hat{\mathbf r}\cdot\mathbf r_a)\,\mathbf r_a\times\hat{\mathbf r}$.
3. $\sum m_a\mathbf r_a(\mathbf r_a\cdot\hat{\mathbf r}) = \big(\sum m\mathbf r\mathbf r^\top\big)\hat{\mathbf r} = \big(\tfrac12\operatorname{tr}I\cdot\mathbb 1 - I\big)\hat{\mathbf r}$ by (3.8). Crossing with $\hat{\mathbf r}$ kills the identity part:
$$\boldsymbol\tau = \frac{3GM}{r^3}\,\hat{\mathbf r}\times(I\hat{\mathbf r}). \tag{3.34}$$
4. For a symmetric body $I\hat{\mathbf r} = I_1\hat{\mathbf r} + (I_3 - I_1)(\hat{\mathbf e}_3\cdot\hat{\mathbf r})\hat{\mathbf e}_3$. Average over a circular orbit $\hat{\mathbf r} = (\cos\lambda, \sin\lambda, 0)$ in the ecliptic with $\hat{\mathbf e}_3 = (0, \sin\varepsilon, \cos\varepsilon)$: $\langle\boldsymbol\tau\rangle = \frac{3GM}{2r^3}(I_3 - I_1)\sin\varepsilon\cos\varepsilon\,\hat{\mathbf x}$. Setting $d\mathbf L/dt = \boldsymbol\Omega\times\mathbf L$ with $\mathbf L = I_3\omega_3\hat{\mathbf e}_3$ and $\boldsymbol\Omega = \Omega\hat{\mathbf z}$ gives $\hat{\mathbf z}\times\hat{\mathbf e}_3 = -\sin\varepsilon\,\hat{\mathbf x}$ and
$$|\Omega| = \frac32\,\frac{GM}{r^3\,\omega_3}\,H\cos\varepsilon, \quad\text{retrograde.} \tag{3.35}$$

The project checks (3.34) against a direct sum over 4000 point masses, checks the orbit average of step 4 numerically, and evaluates (3.35) in Example 4.7: about 51 arcseconds per year, a period near 26,000 years. Treat the agreement with the observed 50.3″/yr as order of magnitude; the lunar orbit's tilt and the 18.6-year nutation are ignored.

**Rolling is non-holonomic.** A coin, a wheel or a ball rolling without slipping has $\mathbf v_{\rm cm} + \boldsymbol\omega\times(\mathbf r_{\rm contact} - \mathbf r_{\rm cm}) = \mathbf 0$. That is a condition on velocities, and $\boldsymbol\omega$ is not the derivative of any coordinate (§3.2), so in general it does not integrate to a condition on positions (CM-3 §3.2's rolling disk). CM-3's method of substituting constraints into $\mathcal L$ therefore fails; you must use Newton–Euler equations with the contact force, or velocity constraints with Lagrange multipliers (CM-3 §3.8 and §4.10; Goldstein §2.4). The project's "Your turn" item 3 sets up the spinning coin this way.

## 4. Worked examples

### Example 4.1: The cube about a corner — a diagonal and a degenerate pair

A uniform cube of mass $M$, side $a$, occupies $[0, a]^3$. Find its inertia tensor about the corner at the origin, and the principal axes.

*Entries.* Density $M/a^3$. $I_{xx} = \frac M{a^3}\int_0^a\!\int_0^a\!\int_0^a(y^2 + z^2)\,dx\,dy\,dz = \frac M{a^3}\big(a\cdot\frac{a^3}3\cdot a + a\cdot a\cdot\frac{a^3}3\big) = \frac23Ma^2$. $I_{xy} = -\frac M{a^3}\int xy\,dV = -\frac M{a^3}\cdot\frac{a^2}2\cdot\frac{a^2}2\cdot a = -\frac14Ma^2$. By symmetry all diagonal entries are $\frac23Ma^2$ and all products $-\frac14Ma^2$:
$$I = Ma^2\begin{pmatrix}\frac23 & -\frac14 & -\frac14\\ -\frac14 & \frac23 & -\frac14\\ -\frac14 & -\frac14 & \frac23\end{pmatrix} = Ma^2\Big(\tfrac{11}{12}\mathbb 1 - \tfrac14 J\Big),$$
where $J$ is the all-ones matrix.

*Check by parallel axis.* The CM is at $\mathbf d = \frac a2(1,1,1)$, $d^2 = \frac34a^2$, $\mathbf d\mathbf d^\top = \frac{a^2}4J$. With $I^{\rm cm} = \frac16Ma^2\mathbb 1$, (3.13) gives $Ma^2\big[(\frac16 + \frac34)\mathbb 1 - \frac14J\big] = Ma^2\big[\frac{11}{12}\mathbb 1 - \frac14J\big]$. ✓

*Principal axes.* $J(1,1,1)^\top = 3(1,1,1)^\top$, and $J\mathbf v = \mathbf 0$ for every $\mathbf v\perp(1,1,1)$. So: along the body diagonal $\frac1{\sqrt3}(1,1,1)$, $I_1 = Ma^2(\frac{11}{12} - \frac34) = \frac16Ma^2$; for **every** direction perpendicular to the diagonal, $I_2 = I_3 = \frac{11}{12}Ma^2$. The smallest moment equals $I^{\rm cm}$, as it must: the diagonal passes through the CM. The degenerate pair means that about its corner the cube is a symmetric top, with the diagonal as its symmetry axis (it has 3-fold symmetry about it; §3.5). Triangle check: $\frac16 + \frac{11}{12} = \frac{13}{12}\ge\frac{11}{12}$. ✓ The project's Monte Carlo cube gives eigenvalues $0.167, 0.917, 0.917$ (units $Ma^2$) and the exact parallel-axis identity on the same sample.

### Example 4.2: A plate spun about its diagonal — the bearing torque

A rectangular steel plate, $M = 2$ kg, $a = 0.3$ m by $b = 0.1$ m, is welded to a shaft along its diagonal through its center and spun at 50 Hz. What does the shaft ask of its bearings?

*Inertia.* Axes along the edges ($x$ along $a$, $y$ along $b$, $z$ normal) are principal by the mirror symmetries. From §3.7, $I = \frac M{12}\operatorname{diag}(b^2, a^2, a^2 + b^2) = \operatorname{diag}(0.0016667,\ 0.015,\ 0.016667)$ kg m².

*Angular velocity.* $|\boldsymbol\omega| = 2\pi\cdot50 = 314.16$ rad/s along $\hat{\mathbf n} = (a, b, 0)/\sqrt{a^2 + b^2} = (0.94868, 0.31623, 0)$, so $\boldsymbol\omega = (298.04,\ 99.346,\ 0)$ rad/s.

*Angular momentum.* $\mathbf L = I\boldsymbol\omega = (0.49673,\ 1.4902,\ 0)$ kg m²/s, $|\mathbf L| = \frac{M}{12}ab\,|\boldsymbol\omega| = 1.5708$ kg m²/s. The angle between $\mathbf L$ and $\boldsymbol\omega$: $\cos\gamma = \frac{\mathbf L\cdot\boldsymbol\omega}{|\mathbf L||\boldsymbol\omega|} = \frac{2ab}{a^2 + b^2} = 0.6$, so $\gamma = 53.13°$. $\mathbf L$ is closer to the long axis's perpendicular than $\boldsymbol\omega$ is, because the plate is "easier" to turn about $x$ ($I_{xx}$ small) and so stores less angular momentum in that direction.

*Torque.* $\mathbf L$ is fixed in the body, so it rotates at $\boldsymbol\omega$, and $d\mathbf L/dt = \boldsymbol\omega\times\mathbf L$ (this is (3.16) with $\dot{\boldsymbol\omega} = 0$). Its magnitude is $|\boldsymbol\omega||\mathbf L|\sin\gamma = 314.16\times1.5708\times0.8 = 394.8$ N m, or in closed form $\omega^2\frac{M}{12}ab\,\frac{a^2 - b^2}{a^2 + b^2}$. Kinetic energy $T = \frac12\boldsymbol\omega\cdot\mathbf L = 148.04$ J.

*What it means.* A 2 kg plate demands almost 400 N m of torque, rotating with the shaft at 50 Hz. If the bearings are 0.4 m apart, each carries a radial force of $394.8/0.4 = 987$ N that sweeps around once per revolution. The plate is **statically** balanced (CM on the axis) but **dynamically** unbalanced (the axis is not principal). Balancing machines measure exactly this torque and add mass until the shaft is a principal axis.

### Example 4.3: The frisbee — wobble versus spin

A thin uniform disk ($I_3 = 2I_1$ by (3.14)) is thrown spinning at 5 rev/s with $\boldsymbol\omega$ tilted $5°$ from its symmetry axis. Take $I_1 = 1$ (units cancel). $|\boldsymbol\omega| = 31.416$ rad/s, so $\omega_3 = 31.416\cos5° = 31.296$ and $\omega_\perp = 31.416\sin5° = 2.7381$ rad/s.

*Body frame.* By (3.18), $\Omega_b = \omega_3(2I_1 - I_1)/I_1 = \omega_3 = 31.296$ rad/s.

*Space frame.* $|\mathbf L| = \sqrt{(I_1\omega_\perp)^2 + (I_3\omega_3)^2} = I_1\omega_3\sqrt{\tan^25° + 4}$, so by (3.19) $\Omega_s/\omega_3 = \sqrt{\tan^25° + 4} = 2.0019$ and $\Omega_s = 62.65$ rad/s. The symmetry axis makes $\theta_L = \arctan(\tan5°/2) = 2.505°$ with $\mathbf L$, half the tilt of $\boldsymbol\omega$.

*What you see.* The outside observer sees the disk's axis wobble at $\dot\phi = \Omega_s$, while the disk's markings turn relative to the wobbling frame at $\dot\psi = -\Omega_b = -\omega_3$ (§3.9). Wobble rate : spin rate $= 2.0019 : 1$, i.e. 2:1 for small tilt. This is the cafeteria-plate calculation Feynman recounts in *Surely You're Joking, Mr. Feynman!*; he stated the ratio the wrong way round, which is a good warning about $\Omega_b$ versus $\Omega_s$. The project integrates Euler's equations with the attitude and measures $\Omega_b = 31.300$ and $\Omega_s/\omega_3 = 2.002$.

### Example 4.4: The tennis racket, $I_1 : I_2 : I_3 = 1 : 2 : 3$

*Stable axes.* From (3.20), spin about $\hat{\mathbf e}_1$: $\sqrt{-\lambda} = \omega_0\sqrt{(I_2 - I_1)(I_3 - I_1)/(I_2I_3)} = \omega_0\sqrt{1\cdot2/6} = 0.5774\,\omega_0$. About $\hat{\mathbf e}_3$: $\omega_0\sqrt{(I_3 - I_1)(I_3 - I_2)/(I_1I_2)} = \omega_0\sqrt{2\cdot1/2} = \omega_0$. Small wobbles, bounded.

*Unstable axis.* About $\hat{\mathbf e}_2$: $\sigma = \omega_0\sqrt{(2 - 1)(3 - 2)/(1\cdot3)} = \omega_0/\sqrt3 = 0.5774\,\omega_0$. A perturbation grows by a factor 100 in $\ln100/\sigma = 7.98/\omega_0$. The growing mode has $\delta_3 = I_1\sigma\delta_1/((I_2 - I_3)\omega_0) = -\sigma\delta_1/\omega_0$ (from the first linearized equation).

*With numbers.* Toss it at $\omega_0 = 2$ rev/s $= 12.566$ rad/s: $\sigma = 7.255$ s⁻¹, and a 1% misalignment reaches order one in $0.635$ s, about one and a quarter turns. That is why you cannot throw a racket cleanly about that axis. The project fits $\sigma = 0.5773\,\omega_0$ from the simulation and shows $\omega_2$ reversing sign.

### Example 4.5: A toy gyroscope — slow precession with numbers

A disk, $M = 0.1$ kg, radius $R = 0.03$ m, sits on a light axle whose tip rests on a pivot, with the disk's center $\ell = 0.04$ m from the pivot. It spins at 3000 rpm with the axle $60°$ from vertical.

*Moments about the pivot.* $I_3 = \frac12MR^2 = 4.5\times10^{-5}$ kg m². $I_1 = \frac14MR^2 + M\ell^2 = 2.25\times10^{-5} + 1.6\times10^{-4} = 1.825\times10^{-4}$ kg m² (disk about a diameter, then (3.13)). $Mg\ell = 0.1\times9.81\times0.04 = 0.03924$ N m. $\omega_3 = 314.16$ rad/s, $p_\psi = I_3\omega_3 = 0.014137$ kg m²/s.

*Steady precession* from (3.31) with $\cos\theta_0 = 0.5$: discriminant $p_\psi^2 - 4I_1Mg\ell\cos\theta_0 = 1.9986\times10^{-4} - 1.4323\times10^{-5}$, so $\dot\phi_- = 2.827$ rad/s (period 2.222 s) and $\dot\phi_+ = 152.1$ rad/s. The fast-top approximations give $Mg\ell/p_\psi = 2.776$ (2% low) and $p_\psi/(I_1\cos\theta_0) = 154.9$.

*Nutation.* From (3.32) at $\dot\phi_0 = 2.827$: $\omega_{\rm nut} = 74.68$ rad/s, against the fast-top estimate $p_\psi/I_1 = 77.46$. Kicking the steadily precessing top with $\dot\theta = 0.8$ rad/s gives a nutation the simulation measures at 74.69 rad/s, swinging between $59.39°$ and $60.62°$.

*Sleeping.* (3.33): $\omega_3 > 2\sqrt{Mg\ell I_1}/I_3 = 118.9$ rad/s $= 1136$ rpm. The project releases the top from rest at $0.02$ rad and bisects on $\omega_3$ for "the lower turning point of (3.30) is beyond $0.3$ rad"; it finds 117.9 rad/s, about 1% lower. The threshold is not a cliff. Released at rest exactly upright, $p_\phi = p_\psi$ and $E' = Mg\ell$, so $b = a$, $\alpha = \beta$ in (3.30) and $f(u) = (1 - u)^2\big[\beta(1 + u) - a^2\big]$: below threshold the top falls only to $\cos\theta_{\max} = I_3^2\omega_3^2/(2I_1Mg\ell) - 1$, which grows continuously from $\theta_{\max} = 0$ at (3.33). The criterion $\theta_{\max} > 0.3$ rad therefore trips at $\omega_3 = 118.9\cos(0.15) = 117.6$ rad/s, and starting at $0.02$ rad instead of exactly upright moves that to 117.9.

*Released from rest* at $30°$ with the same spin, the axis drops to $32.26°$ and comes back, tracing cusps; its mean precession rate is 2.843 rad/s, close to the slow rate.

### Example 4.6: Euler-angle $\boldsymbol\omega$ at a specific configuration

Take $\phi = 30°$, $\theta = 60°$, $\psi = 45°$ and rates $(\dot\phi, \dot\theta, \dot\psi) = (1, 0.5, 2)$ rad/s.

*Body components* (3.23): $\omega_1 = 1\cdot\frac{\sqrt3}2\cdot\frac{\sqrt2}2 + 0.5\cdot\frac{\sqrt2}2 = 0.61237 + 0.35355 = 0.96593$; $\omega_2 = 0.61237 - 0.35355 = 0.25882$; $\omega_3 = 1\cdot0.5 + 2 = 2.5$. (The first two are $\cos15°$ and $\sin15°$; not a coincidence, since $\frac{\sqrt3}2\sin45° + \frac12\cos45° = \sin75°$.)

*Space components* (3.24): $\omega_x = 0.5\cdot\frac{\sqrt3}2 + 2\cdot\frac{\sqrt3}2\cdot\frac12 = 0.43301 + 0.86603 = 1.29904$; $\omega_y = 0.5\cdot\frac12 - 2\cdot\frac{\sqrt3}2\cdot\frac{\sqrt3}2 = 0.25 - 1.5 = -1.25$; $\omega_z = 1 + 2\cdot\frac12 = 2$.

*Checks.* $|\boldsymbol\omega_{\rm body}|^2 = 0.93301 + 0.06699 + 6.25 = 7.25$; $|\boldsymbol\omega_{\rm space}|^2 = 1.6875 + 1.5625 + 4 = 7.25$; formula: $1 + 0.25 + 4 + 2\cdot1\cdot2\cdot\frac12 = 7.25$. ✓ So $|\boldsymbol\omega| = 2.6926$ rad/s. The project also confirms $R\,\boldsymbol\omega_{\rm body} = \boldsymbol\omega_{\rm space}$ with (3.21).

### Example 4.7: Precession of the equinoxes

Data: $GM_\odot = 1.32712\times10^{20}$ m³/s², $r_\odot = 1.496\times10^{11}$ m; $GM_{\rm moon} = 4.9028\times10^{12}$ m³/s², $r_{\rm moon} = 3.844\times10^8$ m; $\omega_3 = 7.2921\times10^{-5}$ rad/s; $H = 1/305.45$; $\varepsilon = 23.44°$.

*Sun.* $GM_\odot/r_\odot^3 = 3.964\times10^{-14}$ s⁻². (3.35): $|\Omega_\odot| = 1.5\times3.964\times10^{-14}/7.2921\times10^{-5}\times(1/305.45)\times0.9175 = 2.45\times10^{-12}$ rad/s $= 15.94$″/yr.

*Moon.* $GM/r^3$ is 2.18 times larger: $34.71$″/yr.

*Total.* $50.7$″/yr, so the axis goes around once in $360\times3600/50.7\approx25{,}600$ years. Observed: $50.3$″/yr, 25,772 years. Two inputs did the work: the oblateness $H$, and the fact that the Moon's tidal torque beats the Sun's because $GM/r^3$, not $GM/r^2$, sets it.

## 5. Why it matters and where it connects

- **LA-8 → here.** The principal axes are the spectral theorem (Theorem 3.11) applied to $I$; moments about arbitrary axes are Rayleigh quotients (Theorem 3.17). CM-6's mass matrix $M$ is the same kind of object: principal axes are "normal modes with no stiffness."
- **LA-12 → here.** (3.9) is the $(0,2)$ tensor law; Theorem 3.1 is Appendix A's (A.1) made physical; the addition rule and gimbal lock are why quaternions ($SU(2)$, A.5–A.6) run every flight computer and game engine.
- **QM-5.** The space components of $\mathbf L$ satisfy $\{L_i, L_j\} = \varepsilon_{ijk}L_k$ (CM-8 (3.17)), the classical shadow of $[J_i, J_j] = i\hbar\varepsilon_{ijk}J_k$. The body components satisfy the same algebra with the opposite sign, $\{L_1, L_2\} = -L_3$ and cyclically, and with that bracket $H = \sum_kL_k^2/2I_k$ generates exactly the body-frame Euler equations $\dot{\mathbf L} = \mathbf L\times\boldsymbol\omega$ (the project checks all three statements in Euler-angle canonical coordinates). The sign flip survives quantization: body-fixed angular momentum operators have "anomalous" commutators, which is a standard trap in molecular spectroscopy. The quantum symmetric top, with energies built from $1/(2I_1)$ and $1/(2I_3)$, is the model for molecular rotational spectra. Spin precession in a magnetic field, $\dot{\mathbf S} = \gamma\mathbf S\times\mathbf B$, is the gyroscope's slow precession with a different torque.
- **CM-8.** The heavy top has three degrees of freedom and three conserved quantities in involution, so it is integrable in the sense of CM-8 §3.11 (Liouville–Arnold). CM-8 §3.4 explains why a cyclic coordinate drops out cleanly in the Hamiltonian picture, which is the clean version of our "substitute into $E$, not $\mathcal L$" rule; the Lagrangian-side fix, the Routhian, is pointed to there (Goldstein §8.3).
- **ME dynamics.** Dynamic balancing (Example 4.2), gyroscopic bearing loads $\boldsymbol\Omega\times\mathbf L$, spin-stabilized and dual-spin satellites (Exercise C1), the recursive Newton–Euler algorithm for robot-arm torques (which is (3.16) applied link by link), and IMU attitude propagation by integrating $\dot R = R[\boldsymbol\omega]_\times$.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills: moments by integration, the hat map, a point-mass tensor, one evaluation of Euler's equations, gimbal lock in words. B1–B8 are the module; B2 (cube about an edge), B4 (your phone's flip rate), B6 (the gyroscope's precession quadratic) and B8 (Monte Carlo hemisphere in `project.py`) are essential. C1 proves the major-axis rule behind Explorer 1's tumble, C2 deploys a boom on a spinning satellite (with code), C3 derives the gyrocompass.

## 7. Coding project

`project.py` builds every object in this module and asserts the numbers above:

1. **Inertia tensor from a point cloud** with `einsum`; symmetry, positivity, $\mathbf L = I\boldsymbol\omega$ against $\sum m\mathbf r\times(\boldsymbol\omega\times\mathbf r)$, $T = \frac12\boldsymbol\omega^\top I\boldsymbol\omega$, the triangle inequalities.
2. **Monte Carlo bodies** (rod, cube, plate, sphere, cone) against §3.7; parallel- and perpendicular-axis theorems exactly on each sample; Example 4.1.
3. **$I' = RIR^\top$** by rotating the cloud; principal axes from `eigh`; Example 4.2's bearing torque.
4. **$\boldsymbol\omega$ from $\dot RR^\top$** by finite differences; non-commuting finite rotations.
5. **Euler's equations torque-free** from near each principal axis (`euler_torque_free.png`), Poinsot polhodes (`poinsot_sphere.png`), the growth rate and wobble frequencies of Example 4.4, the frisbee of Example 4.3 with attitude, the Chandler period.
6. **Euler angles:** $R(\phi,\theta,\psi)$, orthogonality, (3.23)–(3.24) by finite differencing, $\det J = -\sin\theta$, gimbal-lock conditioning.
7. **Heavy top from its Lagrangian** with SymPy (cyclic momenta, (3.28), the $V_{\rm eff}$ reduction checked symbolically), Example 4.5: steady precession, nutation frequency and turning points from the cubic (3.30), the sleeping-top threshold found by bisection, the axis path on the unit sphere.
8. **Equinoxes:** (3.34) against a direct sum, the orbit average of step 4, Example 4.7.
9. **Exercise checks** for the numbers in `solutions.md`.

**Build it yourself first:** write `inertia_tensor(masses, positions)` with two `einsum` calls, and `omega_body_from_rates` from (3.23), before reading them. Then write `euler_rhs` and integrate from near the intermediate axis; you should see the flip before you see the code's figure.

## 8. Where you'll struggle

- **$\mathbf L\nparallel\boldsymbol\omega$ is the hardest intuition in undergraduate mechanics.** Your reflex from fixed-axle problems is $L = I\omega$ with scalars. Fix: every time, compute $\mathbf L = I\boldsymbol\omega$ as a matrix-vector product and *draw both arrows*. Redo the dumbbell of §2 until you can say why $\mathbf L\perp$ rod without computing.
- **Body versus space components.** Euler's equations are the space-frame derivative of $\mathbf L$ written in body components. The torque in (3.17) must be in body components too, and for the top it is easier to use the Lagrangian than to resolve gravity into a rotating frame. Label every vector with its frame in your scratch work.
- **$\Omega_b$ versus $\Omega_s$.** Two different precession rates, seen by two different observers. Feynman mixed them up (Example 4.3). Write which observer you mean.
- **Conventions.** Products of inertia with or without the minus sign; $R$ versus Goldstein's $A = R^\top$; $z$-$x'$-$z''$ versus $z$-$y'$-$z''$ versus aerospace $z$-$y'$-$x''$; intrinsic (new axes, right-multiply) versus extrinsic (fixed axes, left-multiply). Before using any formula from another book, test it on one configuration against Example 4.6.
- **Euler angles are objectively unpleasant;** don't judge yourself by them. They are a bad coordinate system used here because the top's symmetry makes two of them cyclic. They are singular at $\theta = 0$, which for a top is the *sleeping* configuration, the one you most want to simulate. In code, integrate $R$ or a quaternion.
- **The top's effective-potential algebra is error-prone.** Three specific traps: using $I_1$ about the CM instead of the pivot; substituting the momenta into $\mathcal L$ instead of $E$; dropping the constant $p_\psi^2/2I_3$ and then wondering why the energy doesn't match. Let SymPy check every line, as the project does.
- **Chandler and the equinoxes are order-of-magnitude.** 305 days against 433; 50.7″/yr against 50.3″/yr. The first misses by physics (elasticity), the second agrees partly by luck. Say so when you quote them.
- **The stability result is linear.** "Unstable" means exponential growth at first; the flip is then periodic, not a fall into chaos. Don't over-read Theorem 3.7.

## 9. Check yourself

1. Why must $\dot RR^\top$ be antisymmetric, and what vector does it define?
2. Write $I_{ij}$ for point masses, and show $T = \frac12\boldsymbol\omega\cdot\mathbf L$.
3. When is $I$ singular?
4. Why do principal axes exist for every body? Which theorem, which hypothesis?
5. Can $(2, 3, 6)$ be principal moments? Why not?
6. State the parallel-axis theorem in tensor form. Which point must one of the two tensors be about?
7. Derive $(d\mathbf A/dt)_{\rm space} = (d\mathbf A/dt)_{\rm body} + \boldsymbol\omega\times\mathbf A$ in two lines.
8. For a free symmetric top, what are $\Omega_b$ and $\Omega_s$, and which does an outside observer see?
9. Which principal axis is unstable, and what is the growth rate?
10. Write $\boldsymbol\omega$ as a sum of three Euler-angle rates times three unit vectors. Name the vectors.
11. What is $\det J$ for the Euler-rate map, and what goes wrong where it vanishes?
12. For the heavy top, which momenta are conserved, what is $V_{\rm eff}$, and what is the sleeping-top condition?

## 10. Pair with

- **Taylor, *Classical Mechanics*, ch. 10.** The gentlest full treatment: inertia tensor, principal axes, Euler's equations, the free symmetric top and Euler angles, with many worked numbers. Read it alongside §3.3–3.11.
- **Goldstein, Poole & Safko, *Classical Mechanics*, ch. 4–5.** Ch. 4 does rotation matrices, Euler angles (our convention, with $A = R^\top$), the rotating-frame derivative and Cayley–Klein parameters; ch. 5 does Poinsot's construction and the heavy top in full, including the cubic $f(u)$ and the pictures of the three nutation types. The reference for §3.12.
- **Landau & Lifshitz, *Mechanics*, §§31–39.** Forty pages that cover everything here plus the asymmetric top in elliptic functions and rotating frames. Terse; read after these notes.
- **MIT 8.09 (Classical Mechanics III) rigid-body lectures, OCW.** Lectures on the inertia tensor, Euler equations, the symmetric top and its stability, at a pace between Taylor and Goldstein.
- **Kleppner & Kolenkow, *An Introduction to Mechanics*, ch. 7–8.** The best physical intuition for gyroscopes and $\mathbf L\nparallel\boldsymbol\omega$, with the dumbbell and the gyroscope done carefully from $d\mathbf L/dt = \boldsymbol\tau$ before any tensor appears.
