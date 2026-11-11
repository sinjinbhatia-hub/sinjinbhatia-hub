# EM-2 Solutions

Every number below is printed by `project.py` (Part 8 and the "Your turn" cells unless another part is named). $k = 8.988\times10^9$ N m²/C², so $k\times(1\text{ nC}) = 8.988$ V m.

## A. Warm-up

**A1.** (a) From $q_1$: $\boldsymbol{\mathscr r}_1 = (0, 0.4, 0)$ m, $\mathscr r_1 = 0.4$ m, so $\mathbf E_1 = kq_1/\mathscr r_1^2\,\hat{\mathbf y} = 17.975/0.16 = 112.34$ V/m along $\hat{\mathbf y}$.
From $q_2$: $\boldsymbol{\mathscr r}_2 = (0, 0.4, 0) - (0.3, 0, 0) = (-0.3, 0.4, 0)$ m, $\mathscr r_2 = 0.5$ m, $\hat{\boldsymbol{\mathscr r}}_2 = (-0.6, 0.8, 0)$. So $\mathbf E_2 = (kq_2/0.25)\hat{\boldsymbol{\mathscr r}}_2 = -35.95\times(-0.6, 0.8, 0) = (21.57, -28.76, 0)$ V/m.
Sum: $\mathbf E = (21.57,\ 83.58,\ 0)$ V/m, $|\mathbf E| = 86.32$ V/m.
$V = k(q_1/0.4 + q_2/0.5) = 8.988\times(5 - 2) = 26.96$ V.

(b) On the $x$-axis both fields lie along $\hat{\mathbf x}$. Between the charges, $q_1$ pushes a test charge toward $+x$ and $q_2$ pulls it toward $+x$: the fields add. Left of $q_1$, the larger charge $q_1$ is also the nearer one, so $2/x^2 > 1/(x-0.3)^2$ and $q_1$ wins. Right of $q_2$ the two fields are opposite and the smaller charge is nearer, so they can balance:
$$\frac{2}{x^2} = \frac{1}{(x-0.3)^2}\ \Rightarrow\ x - 0.3 = \frac{x}{\sqrt2}\ \Rightarrow\ x = \frac{0.3\sqrt2}{\sqrt2 - 1} = 1.0243\text{ m}.$$
The project evaluates (3.2) there and gets zero.

**A2.** (a) The total flux is $q/\varepsilon_0$ by (3.5). A rotation of the cube about its centre carries the charge to itself and any face to any other, so all six fluxes are equal: $q/(6\varepsilon_0)$ each.

(b) Eight copies of the cube fit around the corner, making a cube of side $2L$ centred on $q$. The flux out of the big cube is $q/\varepsilon_0$, and by symmetry each small cube's outer faces carry $q/(8\varepsilon_0)$. Our cube's outer faces are its three far faces. The three faces that meet at the corner contain the charge in their planes, so $\mathbf E$ (radial from the corner) lies in them and $\mathbf E\cdot\hat{\mathbf n} = 0$: zero flux. The three far faces are permuted by the $120°$ rotation about the cube's diagonal through the corner, so each carries
$$\Phi = \frac{q}{24\varepsilon_0}.$$
The project integrates $E_x$ over a far face and gets $0.0416666667\,q/\varepsilon_0$, which is $1/24$ to ten digits.

(c) Gauss's law gives the integral of $E_\perp$ over a closed surface. On a face of the cube $\mathbf E$ changes in both size and direction from point to point, and no symmetry makes $E_\perp$ constant on any closed surface through the face's centre. One integral cannot give the value at one point.

**A3.** (a) Translations in $x$ and $y$ leave the slab unchanged, so $\mathbf E$ depends only on $z$. A $180°$ rotation about the vertical line through the field point reverses $E_x, E_y$, so $\mathbf E = E_z(z)\hat{\mathbf z}$. The reflection $z\to-z$ gives $E_z(-z) = -E_z(z)$. A pillbox with faces of area $A$ at $\pm z$ has flux $2E_z(z)A$ and encloses $2\rho A|z|$ for $|z|<t$, $2\rho At$ for $|z|\ge t$:
$$E_z = \frac{\rho z}{\varepsilon_0}\ (|z|<t),\qquad E_z = \frac{\rho t}{\varepsilon_0}\operatorname{sign}(z)\ (|z|\ge t).$$
(b) $V(z) = -\int_0^zE_z\,dz'$. Inside: $V = -\rho z^2/(2\varepsilon_0)$. For $z>t$: $V = -\rho t^2/(2\varepsilon_0) - \rho t(z-t)/\varepsilon_0 = -\rho t(2z - t)/(2\varepsilon_0)$, and the same with $|z|$ for $z<-t$. For $\rho>0$, $V\to-\infty$ linearly as $|z|\to\infty$. The slab is infinite, so its charge reaches infinity and $V(\infty) = 0$ cannot be imposed (§3.6, Example 4.2).

(c) At $z = t$: inside value $\rho t/\varepsilon_0$, outside value $\rho t/\varepsilon_0$. Continuous. By (3.10), $E_\perp$ jumps only by $\sigma/\varepsilon_0$, and there is no surface charge at $z = \pm t$: the density there is a finite volume density, and a pillbox of zero thickness encloses none of it. As $t\to0$ with $2\rho t = \sigma$ fixed, the outside field tends to $\sigma/2\varepsilon_0$, the plane of §3.4.

**A4.** $\hat{\mathbf n} = \hat{\mathbf z}$. By (3.11) the tangential part is continuous: $E_{2x} = 3.0\times10^5$, $E_{2y} = 0$. By (3.10), $E_{2z} = E_{1z} + \sigma/\varepsilon_0 = 2.0\times10^5 + 5.0\times10^{-6}/8.854\times10^{-12} = 2.0\times10^5 + 5.647\times10^5$ V/m. So
$$\mathbf E_2 = (3.0,\ 0,\ 7.647)\times10^5\ \text{V/m}.$$
By (3.12), $\partial V_2/\partial z - \partial V_1/\partial z = -\sigma/\varepsilon_0 = -5.647\times10^5$ V/m. Check: $\partial V/\partial z = -E_z$, and $-7.647\times10^5 - (-2.0\times10^5) = -5.647\times10^5$. ✓

**A5.** Series pair: $1/C_s = 1/2 + 1/3 = 5/6$, so $C_s = 1.20\ \mu$F. In parallel with $4\ \mu$F: $C_{\text{eq}} = 5.20\ \mu$F. Total charge $Q = C_{\text{eq}}V = 62.4\ \mu$C.
The series pair has the full 12 V across it, and both its capacitors carry the same charge (the wire between them is isolated and neutral, §3.10): $Q_s = C_s\times12 = 14.4\ \mu$C. Voltages: $14.4/2 = 7.2$ V and $14.4/3 = 4.8$ V, which add to 12. ✓ The $4\ \mu$F capacitor has 12 V and $48\ \mu$C. Charges: $14.4 + 48 = 62.4\ \mu$C. ✓
Energy: $\tfrac12C_{\text{eq}}V^2 = \tfrac12(5.2\times10^{-6})(144) = 374.4\ \mu$J. Individually, $Q^2/2C$ gives $51.84$, $34.56$, $288.00\ \mu$J, which sum to $374.4\ \mu$J. ✓

**A6.** $kq^2/a = 8.988\times10^{-9}/0.1 = 8.988\times10^{-8}$ J.
(a) Corners 1–2–3–4 around the square. Charge 1: free. Charge 2 (adjacent to 1): $kq^2/a$. Charge 3 (adjacent to 2, diagonal from 1): $kq^2/a + kq^2/(a\sqrt2)$. Charge 4 (adjacent to 1 and 3, diagonal from 2): $2kq^2/a + kq^2/(a\sqrt2)$. Total: four sides and two diagonals,
$$W = \frac{kq^2}{a}\left(4 + \frac{2}{\sqrt2}\right) = (4+\sqrt2)\frac{kq^2}{a} = 4.866\times10^{-7}\text{ J}.$$
(b) At any corner, the other three give $V = kq(2/a + 1/(a\sqrt2)) = 8.988\times(20 + 7.071) = 243.3$ V. Then $\tfrac12\sum q_iV_i = \tfrac12\times4\times10^{-9}\times243.3 = 4.866\times10^{-7}$ J. ✓
(c) $\sum_iq_iV(\mathbf r_i)$ contains each pair twice, once as $q_i$ sitting in $q_j$'s potential and once the other way round, so it is twice the assembly work.

## B. Core

**B1.** (a) Newton's law has Coulomb's form with $kq\to-GM$. For a mass distribution, $\mathbf g(\mathbf r) = -G\int\rho_m(\mathbf r')\hat{\boldsymbol{\mathscr r}}/\mathscr r^2\,d\tau'$. Follow §3.3: take the divergence inside the integral, use $\nabla\cdot(\hat{\boldsymbol{\mathscr r}}/\mathscr r^2) = 4\pi\delta^3(\mathbf r - \mathbf r')$, and the delta function sifts out $\rho_m(\mathbf r)$:
$$\nabla\cdot\mathbf g = -4\pi G\rho_m.$$
Integrate over a volume and apply the divergence theorem: $\oint\mathbf g\cdot d\mathbf a = -4\pi GM_{\text{enc}}$. The dictionary is $1/\varepsilon_0\to-4\pi G$.

(b) The argument of §3.3's converse applies unchanged: rotations about the centre fix the shell, so $\mathbf g = g(r)\hat{\mathbf r}$. On a sphere of radius $r$, $4\pi r^2g(r) = -4\pi GM_{\text{enc}}$. For $r>R$, $M_{\text{enc}} = M$ and $g = -GM/r^2$: the field of a point mass $M$ at the centre. For $r<R$, $M_{\text{enc}} = 0$ and $g = 0$. ∎ Newton proved both halves geometrically (*Principia*, Book I, Propositions 70–71); Gauss's law makes each a two-line argument.

(c) Inside a uniform ball, $M_{\text{enc}} = Mr^3/R^3$, so $\mathbf g = -(GM/R^3)\,\mathbf r$, the same shape as Example 4.1's interior field. In a tunnel along a diameter, with $x$ the position along the tunnel, $\ddot x = -(GM/R^3)\,x$. That is simple harmonic motion with $\omega = \sqrt{GM/R^3}$:
$$T = 2\pi\sqrt{R^3/GM} = 5061\text{ s} = 84.3\text{ min},\qquad\text{one way } T/2 = 42.2\text{ min}.$$
Speed at the centre: amplitude times $\omega$, $v = R\omega = \sqrt{GM/R} = 7.91$ km/s. (Check: surface gravity $GM/R^2 = 9.82$ m/s².)

(d) Let the chord pass at distance $D$ from the centre. Write $\mathbf r = D\hat{\mathbf n} + x\hat{\mathbf t}$, with $\hat{\mathbf t}$ along the chord, $\hat{\mathbf n}$ perpendicular to it, and $x$ measured from the chord's midpoint. The component of $\mathbf g = -(GM/R^3)\mathbf r$ along the chord is $-(GM/R^3)x$. The perpendicular component is cancelled by the frictionless wall. So $\ddot x = -(GM/R^3)x$ again, with the same $\omega$ and $T = 84.3$ min for every chord. A satellite skimming the surface has $v^2/R = GM/R^2$, so $v = \sqrt{GM/R} = 7.91$ km/s and period $2\pi R/v = 2\pi\sqrt{R^3/GM}$, the same 84.3 min. The reason: uniform circular motion projected onto any line is simple harmonic motion with the same $\omega$. The falling stone is the satellite's shadow. (The real Earth is denser toward the centre, so these numbers belong to the model, not the planet.)

**B2.** (a) With $b = a + d$:
$$C = 4\pi\varepsilon_0\frac{a(a+d)}{d} = \frac{4\pi\varepsilon_0a^2}{d}\left(1 + \frac da\right)\ \to\ \frac{\varepsilon_0(4\pi a^2)}{d}\quad(d/a\to0).$$
Cylinder: $\ln(1 + d/a) = d/a - d^2/2a^2 + \cdots$, so $C = 2\pi\varepsilon_0L/\ln(1+d/a) = (2\pi\varepsilon_0La/d)(1 + d/2a + \cdots)\to\varepsilon_0(2\pi aL)/d$.
Why both must agree: when the gap is small compared with the radius of curvature, each patch of one conductor faces a patch of the other across a nearly flat gap. Gauss's law with a pillbox gives $E = \sigma/\varepsilon_0$ in the gap, and $V = Ed = Qd/(\varepsilon_0A)$, which is the parallel-plate calculation.

(b) $C = 4\pi\varepsilon_0a/(1 - a/b)\to4\pi\varepsilon_0a$. For the Earth: $4\pi\times8.854\times10^{-12}\times6.371\times10^6 = 709\ \mu$F. A planet has the capacitance of a modest electrolytic capacitor, because $C$ grows only linearly with size and $\varepsilon_0$ is small in SI units.

(c) The ground is a conductor, so Theorem 3.4(v) gives $\sigma = \varepsilon_0\,\mathbf E\cdot\hat{\mathbf n}$ with $\hat{\mathbf n}$ pointing up, out of the conductor. With $\mathbf E = -100\,\hat{\mathbf n}$ V/m, $\sigma = -8.854\times10^{-10}$ C/m². Over the whole surface, $Q = 4\pi R^2\sigma = -4.52\times10^5$ C. The ground is negative; a roughly matching positive charge is spread through the atmosphere above it, which is why the field dies away with altitude.

(d) The project gives $2\pi\varepsilon_0/\ln4 = 40.13$ pF/m against numerical values $38.84$, $39.62$, $39.86$ pF/m at $h = 0.02, 0.01, 0.005$ (ratios $0.9678$, $0.9872$, $0.9933$). The numerical value is always low. The inner conductor is the set of grid nodes with $s\le a$, a staircase lying inside the circle $s = a$. The outer conductor is the set of nodes with $s\ge b$, lying outside the circle $s = b$. So the effective $a$ is a little smaller and the effective $b$ a little larger, which makes $\ln(b/a)$ larger and $C$ smaller. Each boundary moves by a fraction of $h$, so the error is first order in $h$: halving $h$ roughly halves it, though unevenly, because the set of nodes that fall inside the circles changes irregularly with $h$. The five-point stencil's own error is $O(h^2)$, so the staircase dominates. To do better, place nodes on the boundary (a polar grid for this geometry) or correct the stencil next to curved boundaries.

**B3.** (a) Put the field point at $(0,0,r)$. A source point at polar angle $\theta'$ on the shell is at distance $\mathscr r = \sqrt{R^2 + r^2 - 2Rr\cos\theta'}$, and $da' = R^2\sin\theta'\,d\theta'\,d\phi'$. With $u = \cos\theta'$:
$$V(r) = \frac{\sigma}{4\pi\varepsilon_0}\,2\pi R^2\int_{-1}^{1}\frac{du}{\sqrt{R^2+r^2-2Rru}} = \frac{\sigma R^2}{2\varepsilon_0}\cdot\frac{1}{Rr}\Big[\sqrt{R^2+r^2-2Rru}\,\Big]_{u=1}^{u=-1},$$
using $\frac{d}{du}\sqrt{R^2+r^2-2Rru} = -Rr/\sqrt{R^2+r^2-2Rru}$. The bracket is $\sqrt{(R+r)^2} - \sqrt{(R-r)^2} = (R+r) - |R-r|$.
- $r>R$: the bracket is $2R$, so $V = \sigma R^2/(\varepsilon_0r) = Q/(4\pi\varepsilon_0r)$, using $Q = 4\pi R^2\sigma$.
- $r<R$: the bracket is $2r$, so $V = \sigma R/\varepsilon_0 = Q/(4\pi\varepsilon_0R)$, a constant.

(b) The step $\sqrt{(R-r)^2} = |R-r|$. That absolute value is the only place where inside and outside differ. Writing $R - r$ for every $r$ would give the inside answer everywhere. Differentiating gives $E = kQ/r^2$ outside and $0$ inside, matching §3.4 without any symmetry argument. Project Part 0 checks both cases with SymPy.

(c) A ball of density $\rho = 3Q/(4\pi R^3)$ is a stack of shells of radius $r'$ and charge $dq = 4\pi\rho r'^2dr'$. The centre is inside every shell, so each contributes $k\,dq/r'$:
$$V(0) = \int_0^Rk\frac{4\pi\rho r'^2}{r'}dr' = 2\pi k\rho R^2 = \frac{3kQ}{2R} = \tfrac32V(R),$$
which is Example 4.1's $V(0) = 134.81$ V for $Q = 1$ nC, $R = 0.1$ m (Part 3).

**B4.** (a) Product rule: $\nabla V = \frac{1}{4\pi\varepsilon_0}\left[r^{-3}\nabla(\mathbf p\cdot\mathbf r) + (\mathbf p\cdot\mathbf r)\nabla(r^{-3})\right]$. With $\mathbf p$ constant, $\nabla(p_xx + p_yy + p_zz) = \mathbf p$. And $\nabla(r^{-3}) = -3r^{-4}\nabla r = -3r^{-4}\,\mathbf r/r = -3r^{-5}\mathbf r$. So
$$\mathbf E = -\nabla V = \frac{1}{4\pi\varepsilon_0}\left[\frac{3(\mathbf p\cdot\mathbf r)\mathbf r}{r^5} - \frac{\mathbf p}{r^3}\right] = \frac{3(\mathbf p\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf p}{4\pi\varepsilon_0r^3},$$
which is (4.3). Project Part 0 confirms it with SymPy.

(b) A field line in a plane containing the axis ($E_\phi = 0$) has line element $(dr,\ r\,d\theta)$ parallel to $(E_r, E_\theta)$, so by (4.2)
$$\frac{dr}{r\,d\theta} = \frac{E_r}{E_\theta} = \frac{2\cos\theta}{\sin\theta}\ \Rightarrow\ \frac{dr}{r} = 2\frac{\cos\theta}{\sin\theta}d\theta\ \Rightarrow\ \ln r = 2\ln\sin\theta + \text{const},$$
so $r = r_0\sin^2\theta$, where $r_0$ is where the line crosses the equator. Each line leaves the dipole near $\theta = 0$, swings out to $r_0$, and returns near $\theta = \pi$.

(c) The project gives a maximum deviation of $4.4\times10^{-5}$ for $d = 1$ mm and $1.76\times10^{-4}$ for $d = 2$ mm, a ratio of $4.00$. The traced line belongs to two charges, not an ideal dipole. Their potential is odd under $z\to-z$, so it contains only odd multipoles: the next term after the dipole is the octupole, smaller by a relative factor of order $(d/r)^2$ (Example 4.3 found exactly $1 - d^2/4r^2$ on the axis). The deviation is therefore proportional to $d^2$ and is largest at the innermost point, $r = 0.1$ m, where $(d/r)^2 = 10^{-4}$. The ratio of exactly 4 shows that integration error (tolerance $10^{-10}$) plays no part.

**B5.** (a) $E^2 = |\mathbf E_1 + \mathbf E_2|^2 = E_1^2 + E_2^2 + 2\mathbf E_1\cdot\mathbf E_2$. Multiply by $\varepsilon_0/2$ and integrate: $W = W_1 + W_2 + \varepsilon_0\int\mathbf E_1\cdot\mathbf E_2\,d\tau$. ∎ The energy is quadratic in the fields, so it does not superpose. Doubling every charge quadruples $W$.

(b) Write $\mathbf E_1 = -\nabla V_1$ and use $\nabla\cdot(V_1\mathbf E_2) = \nabla V_1\cdot\mathbf E_2 + V_1\nabla\cdot\mathbf E_2$:
$$\varepsilon_0\int\mathbf E_1\cdot\mathbf E_2\,d\tau = -\varepsilon_0\oint V_1\mathbf E_2\cdot d\mathbf a + \varepsilon_0\int V_1\nabla\cdot\mathbf E_2\,d\tau.$$
On a sphere of radius $r\to\infty$ the surface term is of order $(1/r)(1/r^2)r^2\to0$, and $\nabla\cdot\mathbf E_2 = \rho_2/\varepsilon_0$, so the interaction term is $\int\rho_2V_1\,d\tau$. Swapping the labels gives $\int\rho_1V_2\,d\tau$. For point charges, $\rho_2 = q_2\delta^3(\mathbf r - \mathbf r_2)$ and the term is $q_2V_1(\mathbf r_2) = kq_1q_2/\mathscr r_{12}$, the pair energy of (3.14). $W_1$ and $W_2$ are the self-energies of each point charge, both infinite (§3.8). So (3.16) for point charges equals (3.14) plus the self-energies, and that is the whole difference between the two formulas.

(c) Shell 1 ($+q$ at radius $a$) has field $kq/r^2$ for $r>a$, so
$$W_1 = \frac{\varepsilon_0}{2}\int_a^\infty\frac{k^2q^2}{r^4}4\pi r^2dr = \frac{kq^2}{2a},$$
using $4\pi\varepsilon_0k = 1$. Likewise $W_2 = kq^2/(2b)$. Interaction: shell 2's potential is $-kq/b$ everywhere inside radius $b$, including on shell 1, so $\int\rho_1V_2\,d\tau = q(-kq/b) = -kq^2/b$. (Directly: $\mathbf E_2 \ne 0$ only for $r>b$, where $\mathbf E_1\cdot\mathbf E_2 = -k^2q^2/r^4$, and $\varepsilon_0\int_b^\infty(-k^2q^2/r^4)4\pi r^2dr = -kq^2/b$. ✓) Sum:
$$W = \frac{kq^2}{2a} + \frac{kq^2}{2b} - \frac{kq^2}{b} = \frac{kq^2}{2}\left(\frac1a - \frac1b\right).$$
Directly, the field is $kq/r^2$ only for $a<r<b$: $\frac{\varepsilon_0}{2}\int_a^b\frac{k^2q^2}{r^4}4\pi r^2dr = \frac{kq^2}{2}\left(\frac1a - \frac1b\right)$. ✓ And with $C = 4\pi\varepsilon_0ab/(b-a)$: $\dfrac{q^2}{2C} = \dfrac{q^2(b-a)}{8\pi\varepsilon_0ab} = \dfrac{kq^2}{2}\left(\dfrac1a - \dfrac1b\right)$. ✓ The interaction term is negative and larger in size than $W_2$: opposite charges lower the energy by sitting close.

**B6.** (a) $\sigma = Q/(4\pi R^2)$, $\mathbf E_1 = 0$ inside, $\mathbf E_2 = (\sigma/\varepsilon_0)\hat{\mathbf r}$ outside. By (3.17),
$$\mathbf f = \frac{\sigma}{2}\cdot\frac{\sigma}{\varepsilon_0}\hat{\mathbf r} = \frac{\sigma^2}{2\varepsilon_0}\hat{\mathbf r} = \frac{Q^2}{32\pi^2\varepsilon_0R^4}\hat{\mathbf r},$$
outward for either sign of $Q$.

(b) The charged sphere has $W = \tfrac12QV = Q^2/(8\pi\varepsilon_0R)$. Let $R$ grow by $dR$ with $Q$ fixed: $dW = -Q^2/(8\pi\varepsilon_0R^2)\,dR$. The electric forces do work $f\cdot4\pi R^2\cdot dR = Q^2/(8\pi\varepsilon_0R^2)\,dR$. Energy conservation needs that work to equal $-dW$, and it does. ✓

(c) By symmetry the net force on the northern hemisphere is along $\hat{\mathbf z}$:
$$F_z = \int_0^{2\pi}\!\!\int_0^{\pi/2}f\cos\theta\,R^2\sin\theta\,d\theta\,d\phi = 2\pi fR^2\cdot\tfrac12 = \pi R^2f = \frac{Q^2}{32\pi\varepsilon_0R^2}.$$
($\pi R^2f$ is the pressure times the projected area, the equatorial disk.) The forces among the northern hemisphere's own charges cancel in pairs, so this is the force exerted by the southern hemisphere. For $Q = 1\ \mu$C, $R = 0.1$ m: $F = 0.1123$ N, which the project also gets by direct quadrature.

**B7.** (a) A Gaussian surface inside the metal enclosing cavity $a$ has $\mathbf E = 0$ on it, so it encloses zero charge, and the wall of cavity $a$ carries $-q_a$. Likewise the wall of cavity $b$ carries $-q_b$. The conductor is neutral, so the outer surface carries $q_a + q_b$. To find the distributions, guess and verify. Guess $\sigma_a = -q_a/(4\pi a^2)$, $\sigma_b = -q_b/(4\pi b^2)$, $\sigma_R = (q_a+q_b)/(4\pi R^2)$, all uniform. Check that the guess gives $\mathbf E = 0$ in the metal:
- $q_a$ plus its uniform wall form a spherically symmetric system with zero net charge, so by §3.4 they produce zero field outside the cavity. The same holds for $q_b$ and its wall.
- The uniform outer shell produces zero field inside radius $R$ (§3.4).

So the field in the metal is zero, each wall is an equipotential, and every surface carries its required total. By the second uniqueness theorem (§3.11, proved in EM-3) this is the distribution.

(b) Outside: only the outer shell contributes, so $\mathbf E = k(q_a+q_b)\hat{\mathbf r}/r^2$, with $r$ measured from the centre of the big sphere. In cavity $a$: the wall's uniform charge gives zero field inside it, and everything outside the cavity gives zero field there too, so $\mathbf E = kq_a\hat{\mathbf r}_a/r_a^2$, measured from the cavity's centre. Cavity $b$ is the same with $q_b$.

(c) Zero for both. The field at $q_a$, leaving out its own, is the field of its wall (zero at the centre of a uniform shell) plus that of everything else (zero inside cavity $a$). The two charges do not feel each other: the metal shields both.

(d) The outer surface charge redistributes (its total stays $q_a + q_b$, but it is no longer uniform), and the field outside changes. $\sigma_a$, $\sigma_b$, the cavity fields, and the forces on $q_a$, $q_b$ are unchanged, because the field in each cavity is fixed by the charge in it and by the wall being an equipotential (Theorem 3.5's reasoning). $q_c$ itself feels a force.

(e) $\sigma_a$ becomes nonuniform, larger on the side nearer $q_a$ (Example 4.6's formula). The field in cavity $a$ changes, and $q_a$ now feels a force toward the nearest part of the wall, its induced charge. Nothing else changes: $\sigma_b$, cavity $b$'s field, the uniform $\sigma_R$, the outside field, and the zero force on $q_b$ (§3.11).

**B8.** (a) Put the ring in the $xy$-plane, centred at the origin. Every element $dq$ is at distance $\sqrt{z^2+R^2}$ from $(0,0,z)$. Diametrically opposite elements cancel each other's horizontal components. The vertical component of each is $k\,dq\,z/(z^2+R^2)^{3/2}$, so
$$E_z = \frac{kQz}{(z^2+R^2)^{3/2}}.$$
Maximum: $\dfrac{d}{dz}\left[z(z^2+R^2)^{-3/2}\right] = (z^2+R^2)^{-5/2}(R^2 - 2z^2) = 0$ at $z = R/\sqrt2$. For $Q = 1$ nC, $R = 0.1$ m: the maximum is $345.9$ V/m at $z = 0.0707$ m, and $E(0.1\text{ m}) = 317.8$ V/m. The project finds the maximum numerically.

(b) A ring of radius $s$ and width $ds$ carries $dQ = 2\pi\sigma s\,ds$. Using (a) and $2\pi k = 1/(2\varepsilon_0)$:
$$E_z = 2\pi k\sigma z\int_0^R\frac{s\,ds}{(z^2+s^2)^{3/2}} = 2\pi k\sigma z\left[\frac1z - \frac{1}{\sqrt{z^2+R^2}}\right] = \frac{\sigma}{2\varepsilon_0}\left(1 - \frac{z}{\sqrt{z^2+R^2}}\right).$$
As $z\to0^+$, $E_z\to\sigma/2\varepsilon_0$: up close the disk looks like an infinite plane (§3.4). For $z\gg R$, $z/\sqrt{z^2+R^2} = (1 + R^2/z^2)^{-1/2} \approx 1 - R^2/2z^2$, so $E_z\approx(\sigma/2\varepsilon_0)(R^2/2z^2) = \sigma\pi R^2/(4\pi\varepsilon_0z^2) = kQ/z^2$: a point charge.

(c) The disk's symmetries are rotations about one axis and one reflection. No family of rigid motions carries every point of a closed surface to every other point while fixing the disk, so there is no surface on which $E_\perp$ is constant (§3.4). Gauss's law still holds but gives one equation for an unknown function.

(d) With $\sigma/2\varepsilon_0 = 1797.5$ V/m, (b) gives $1618.65$, $993.64$, $526.48$, $8.92$ V/m at $z = 0.01, 0.05, 0.1, 1$ m. At 1 m the point-charge value is $kQ/z^2 = 8.988$ V/m. The disk's field is lower by the next term of the expansion, a factor $1 - \tfrac34(R/z)^2$. The superposition of $160{,}000$ point charges agrees with (b) to relative errors $2.9\times10^{-5}$, $2.1\times10^{-6}$, $1.0\times10^{-6}$, $2.3\times10^{-8}$. The error is largest at $z = 1$ cm because there the integrand varies on the scale $z$, which is closest to the grid spacing.

## C. Challenge

**C1.** (a) When the ball has grown to radius $r$ it holds $q(r) = Qr^3/R^3$, and its potential at its surface is $kq(r)/r$ (spherical symmetry, §3.4). Bringing $dq$ from infinity and spreading it as a thin layer costs $k\,q(r)\,dq/r$; the layer's own energy is of order $dq^2$ and drops out. With $dq = 3Qr^2dr/R^3$:
$$W = \int_0^Rk\frac{Qr^3}{R^3}\cdot\frac{3Qr^2}{R^3}\cdot\frac{dr}{r} = \frac{3kQ^2}{R^6}\int_0^Rr^4dr = \frac35\frac{kQ^2}{R}.$$
No $\tfrac12$ appears, because each increment pays only for the charge already present (each pair once, as in (3.14)). Project Part 0 does this integral in SymPy.

(b) Set the electrostatic energy equal to $m_ec^2$. Ball: $\tfrac35ke^2/R = m_ec^2$, so $R = \tfrac35r_e$. Shell: $\tfrac12ke^2/R = m_ec^2$, so $R = \tfrac12r_e$. With $ke^2 = 1.4400$ MeV fm and $m_ec^2 = 0.51100$ MeV, $r_e = 2.8179$ fm, so the ball has radius $1.691$ fm and the shell $1.409$ fm.

(c) Three reasons.
1. Scattering experiments show no electron structure down to distances far below 1 fm. The electron behaves as a point.
2. Quantum mechanics takes over long before $r_e$. The reduced Compton wavelength is $\hbar/m_ec = 386.16$ fm, which is $137.036\,r_e$ (that factor is $1/\alpha$, the inverse fine-structure constant). Confining an electron to a region smaller than this takes momentum of order $m_ec$, enough energy to create electron–positron pairs, so a classical charged ball is the wrong picture two orders of magnitude before you reach $r_e$. In quantum electrodynamics the self-energy diverges only logarithmically as the cutoff shrinks, and renormalization absorbs it into the measured mass (§3.8).
3. Even classically the model is incomplete. The ball and shell give different radii for the same energy, so the number depends on an unobservable internal structure. And by B6 a charged ball pushes itself apart, so it needs some non-electric binding (Poincaré proposed such stresses in 1906) about which electrostatics says nothing.

**C2.** (a) Five-point Laplacian plus explicit Euler, with $r \equiv \alpha\Delta t/h^2$:
$$u^{k+1}_{mn} = (1 - 4r)\,u^k_{mn} + r\left(u^k_{m+1,n} + u^k_{m-1,n} + u^k_{m,n+1} + u^k_{m,n-1}\right).$$
At $r = \tfrac14$ the first term vanishes and $u^{k+1}_{mn}$ is the average of its four neighbours: the Jacobi update for the five-point Laplace equation. Stability: by (b) below, the mode $(p,q)$ is multiplied each step by
$$g_{pq} = 1 - 4r + 2r[\cos p\pi h + \cos q\pi h] = 1 - 4r\left[\sin^2\tfrac{p\pi h}{2} + \sin^2\tfrac{q\pi h}{2}\right].$$
The bracket is largest for $p = q = N$, where it equals $2\sin^2(N\pi h/2)$, which tends to 2 as the grid is refined. So $|g|\le1$ for all modes requires $r\le1/[4\sin^2(N\pi h/2)]$, a bound that tends to $\tfrac14$; $r = \tfrac14$ is stable on every grid. Jacobi is heat flow at the largest step that is stable on all grids.

(b) Let $e^k = V^k - V_{\text{exact}}$. The exact discrete solution is a fixed point of the sweep, so subtracting shows $e$ obeys the same update, with $e = 0$ on the boundary. Put $e_{mn} = \sin(p\pi mh)\sin(q\pi nh)$. From $\sin(\beta(m+1)) + \sin(\beta(m-1)) = 2\sin(\beta m)\cos\beta$, the sum of the four neighbours is $[2\cos(p\pi h) + 2\cos(q\pi h)]\,e_{mn}$, so one sweep multiplies $e$ by
$$\lambda_{pq} = \tfrac12[\cos(p\pi h) + \cos(q\pi h)].$$
The boundary is automatic: $\sin0 = 0$ and $\sin(p\pi(N+1)h) = \sin p\pi = 0$. The $N^2$ vectors with $p, q = 1,\dots,N$ are orthogonal (the discrete sine basis, F2's sine series on a grid), so they are a basis for the interior values, and every error is a sum of them, each shrinking at its own rate. The largest $|\lambda|$ is $\cos(\pi h)$, at $p = q = 1$ (and at $p = q = N$ with the opposite sign). Reducing that component by $10^6$ takes $k = \ln10^{-6}/\ln\cos(\pi/(N+1))$ sweeps: $2686$ for $N = 30$ and $28554$ for $N = 100$. Since $\cos x\approx1 - x^2/2$, $k\approx2(N+1)^2\ln(10^6)/\pi^2\propto N^2$. This is LA-13 (3.24), and it is the factor 4 the project measures when $h$ halves in Part 5.

(c) The project measures $2510$ sweeps. Prediction (b) assumed the slow mode starts with amplitude 1 V, the size of the initial error. It doesn't. The initial error is $-V_{\text{exact}}$, and the slow mode's share of it is the projection $c_{11}$. In the continuum, separation of variables (EM-3) gives $V_{\text{exact}} = \sum_{p\text{ odd}}\frac{4}{p\pi}\sin(p\pi x)\frac{\sinh(p\pi y)}{\sinh p\pi}$, and only $p = 1$ survives the $x$-projection:
$$c_{11} = 4\int_0^1\!\!\int_0^1V\sin\pi x\sin\pi y\,dx\,dy = 4\cdot\frac4\pi\cdot\frac12\cdot\frac{1}{\sinh\pi}\int_0^1\sinh(\pi y)\sin(\pi y)\,dy = \frac{8}{\pi\sinh\pi}\cdot\frac{\sinh\pi}{2\pi} = \frac{4}{\pi^2} = 0.4053.$$
(The $y$-integral: integrate $e^{\pm\pi y}\sin\pi y$ by the standard formula. The two pieces give $(e^\pi+1)/2\pi$ and $(e^{-\pi}+1)/2\pi$, so their half-difference is $\sinh\pi/2\pi$.) On the $30\times30$ grid the project's discrete projection is $0.4046$ in size. The slow mode's peak is $|c_{11}|$ (at the centre), so it falls below $10^{-6}$ V after $\ln(10^{-6}/0.4046)/\ln\cos(\pi/31) = 2510$ sweeps, exactly what is measured. Every other mode is gone by then: the next slowest, $\lambda_{12}\approx1 - 5\pi^2h^2/4$, decays 2.5 times faster in the exponent than $\lambda_{11}\approx1 - \pi^2h^2/2$.

**C3.** (a) Let $\rho$ be the equilibrium distribution, with potential $V$ equal to the constant $V_i$ on conductor $i$, and field $\mathbf E$. Let a trial distribution be $\rho + \delta\rho$, where $\delta\rho$ lives on the conductors and $\int_{\text{conductor }i}\delta\rho\,d\tau = 0$ for each $i$. Let $\delta\mathbf E$, $\delta V$ be the field and potential of $\delta\rho$. The trial field is $\mathbf E + \delta\mathbf E$, and by (3.16)
$$W' = \frac{\varepsilon_0}{2}\int|\mathbf E + \delta\mathbf E|^2d\tau = W + \varepsilon_0\int\mathbf E\cdot\delta\mathbf E\,d\tau + \frac{\varepsilon_0}{2}\int|\delta\mathbf E|^2d\tau.$$
By B5(b), the middle term is $\int V\,\delta\rho\,d\tau$. $\delta\rho$ is nonzero only on the conductors, where $V = V_i$, so it equals $\sum_iV_i\int_i\delta\rho\,d\tau = 0$. Hence $W' = W + \tfrac{\varepsilon_0}{2}\int|\delta\mathbf E|^2d\tau\ge W$, with equality only if $\delta\mathbf E = 0$ everywhere, which means $\delta\rho = \varepsilon_0\nabla\cdot\delta\mathbf E = 0$. ∎

(b) Let $V$ solve $\nabla^2V = 0$ in the region $\Omega$ between the conductors, with the given boundary values, and let $U = V + \eta$ be any function with the same boundary values, so $\eta = 0$ on $\partial\Omega$. Then
$$\int_\Omega|\nabla U|^2 = \int_\Omega|\nabla V|^2 + 2\int_\Omega\nabla V\cdot\nabla\eta + \int_\Omega|\nabla\eta|^2.$$
From $\nabla\cdot(\eta\nabla V) = \nabla\eta\cdot\nabla V + \eta\nabla^2V$ and the divergence theorem, the middle integral is $\oint_{\partial\Omega}\eta\nabla V\cdot d\mathbf a - \int_\Omega\eta\nabla^2V = 0 - 0$. So $\int|\nabla U|^2\ge\int|\nabla V|^2$, with equality only if $\nabla\eta = 0$, which with $\eta = 0$ on the boundary means $\eta = 0$. (If $\Omega$ extends to infinity, require $\eta\to0$ there fast enough for the surface term to vanish, as in §3.8.) ∎ At fixed conductor potentials the true field minimizes the field energy. That is why relaxation (C2) and energy-based finite-element solvers converge to the right answer.

(c) By (a), the all-on-the-surface equilibrium (Theorem 3.4) is the unique minimum of the energy at fixed total charge, so charge placed inside a conductor is in a higher-energy arrangement and the free charges move to the surface, dissipating the difference as heat.
