# EM-1: The Vector Calculus Toolkit: Grad, Div, Curl, the Integral Theorems, the Delta Function, Helmholtz

**Time:** ~10 hours · **Prerequisites:** multivariable calculus (LA-12 §3.19 on $\varepsilon_{ijk}$ helps) · **Unlocks:** EM-2 and every later EM module; completes CM-1 Theorem 3.3

> **Conventions in this module.** SI units, though almost nothing here carries units. Vectors are bold; unit vectors wear hats. $\nabla = (\partial_x, \partial_y, \partial_z)$ is a vector operator. **Index notation:** $i, j, k, \ldots \in \{1, 2, 3\} \leftrightarrow \{x, y, z\}$, $\partial_i = \partial/\partial x_i$, and a repeated index is summed. $\delta_{ij}$ is the Kronecker delta; $\varepsilon_{ijk}$ is the Levi-Civita symbol, $\varepsilon_{123} = +1$, antisymmetric under any swap of two indices. Coordinates are right-handed. $\mathbf r$ is the **field point**, $\mathbf r'$ the **source point**, $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$ the **separation vector**, $\mathscr r = |\boldsymbol{\mathscr r}|$, $\hat{\boldsymbol{\mathscr r}} = \boldsymbol{\mathscr r}/\mathscr r$; with the source at the origin, $\boldsymbol{\mathscr r} = \mathbf r$. $d\tau$ is a volume element. $d\mathbf a = \hat{\mathbf n}\,da$ is an area element whose normal is **outward** on a closed surface and tied to the boundary's direction by the **right-hand rule** on an open one. $d\mathbf l$ is a line element. Spherical $(r, \theta, \phi)$: $\theta$ is the polar angle from $+z$, $\phi$ the azimuth. Cylindrical $(s, \phi, z)$: $s$ is the distance from the $z$-axis. $\nabla^2 = \nabla\cdot\nabla$, and the sign matters: $\nabla^2(1/r) = -4\pi\delta^3(\mathbf r)$. $\delta^3(\mathbf r) = \delta(x)\delta(y)\delta(z)$.

## 0. Where this sits

You can differentiate and integrate in several variables. What you lack is the toolkit *as a toolkit*: which identity to reach for, how to set up a surface integral without freezing, why $\nabla\cdot(\hat{\mathbf r}/r^2)$ is both zero and not zero. Every later EM module is written in this language: (3.32) is EM-2's Gauss's law, Green's identities are EM-3's uniqueness theorems, Theorem 3.17 is EM-5's $\mathbf B = \nabla\times\mathbf A$, and Helmholtz (Theorem 3.15) is why EM-6's Maxwell equations are complete. Theorem 3.16 closes the gap left in CM-1 Theorem 3.3.

## 1. The problem that forced this idea

In 1865 Maxwell published "A Dynamical Theory of the Electromagnetic Field": twenty equations in twenty unknowns, written component by component, so the reader had to notice that several were one statement written three times. Maxwell saw the structure (in the early 1870s he proposed the names "curl" and "convergence," minus our divergence), and his 1873 *Treatise* used some of Hamilton's quaternion notation, which mixed scalar and vector parts awkwardly.

Two people cut the knot independently. Josiah Willard Gibbs printed *Elements of Vector Analysis* privately for his Yale students (1881 and 1884). Oliver Heaviside, a self-taught telegraph engineer, developed the same notation through the 1880s, collected in *Electromagnetic Theory* (1893), and rewrote Maxwell's equations as the four we use now. Written as $\nabla\cdot\mathbf E$ and $\nabla\times\mathbf B$, they say where each field's sources are and how it circulates.

The integral theorems came earlier, from people computing things. Gauss used a special case of the divergence theorem in his 1813 work on the attraction of ellipsoids; Ostrogradsky proved a general version in the 1820s. George Green, a Nottingham miller with almost no formal schooling, self-published his *Essay on the Application of Mathematical Analysis to the Theories of Electricity and Magnetism* in 1828; it introduced the identities of §3.9 and was nearly lost until William Thomson (Kelvin) found a copy in 1845. Kelvin stated what we call Stokes' theorem in an 1850 letter to Stokes; Stokes set it as a Cambridge Smith's Prize examination question in 1854, and his name stuck. Helmholtz (1858), studying vortices in fluids, showed that a flow is fixed by its sources and its vorticity (Theorem 3.15). The point source hiding in $\hat{\mathbf r}/r^2$ was handled with "impulse functions" long before Dirac's *Principles of Quantum Mechanics* (1930) made $\delta$ standard and Laurent Schwartz's theory of distributions (late 1940s) made it rigorous.

You have met most of these pieces. This module exists because E&M is a vector-calculus endurance test, and what fails in Griffiths ch. 2 and after is fluency: picking the identity, setting up the integral, trusting the delta function.

## 2. Intuition

### Four fields

**(i) Spreading** $\mathbf r = (x, y, z)$. Arrows point out and lengthen. A tiny box anywhere: the far faces see a stronger field than the near ones, so more flows out than in. A tiny paddlewheel: same flow on both sides, no spin.

**(ii) Whirlpool** $(-y, x, 0)$. Arrows circle the $z$-axis with speed equal to the distance from it: a turntable at angular velocity 1. A box: in equals out. A paddlewheel anywhere is carried around rigidly and turns. Uniform spin.

**(iii) Vortex** $(-y, x, 0)/(x^2 + y^2)$, speed $1/s$. Still circling. But a paddlewheel at distance $s$ feels faster flow on its inner side; that shear turns it clockwise, while being carried around turns it counterclockwise, and with speed $\propto 1/s$ the two cancel exactly. No spin off the axis, yet the flow plainly goes around. The spin is all on the axis.

**(iv) Coulomb** $\hat{\mathbf r}/r^2$. Arrows point out and *weaken*. A tiny box away from the origin: its outer face is bigger, the field there weaker, and with $1/r^2$ these compensate exactly. Zero net outflow. Yet a sphere of radius $R$ catches flux $(1/R^2)(4\pi R^2) = 4\pi$ for every $R$. All the source is at one point.

§3.4 makes the boxes and paddlewheels precise. **Divergence** is net outflow per unit volume of a tiny box: source density. **Curl** is circulation per unit area of a tiny loop: spin density, twice the local angular velocity. Field (i) has divergence 3, curl 0; field (ii) divergence 0, curl $2\hat{\mathbf z}$; fields (iii) and (iv) have both zero away from their singular sets.

### Adding up boxes

Fill a region with tiny boxes. Each box's net outflow is (divergence) × (volume). Add them. Where two boxes share a face, the flux out of one is the flux into the other: same face, opposite outward normals. Every interior face cancels; only the outer skin survives.

```
 +-----+-----+-----+        shared face: "out of A" = "into B",
 |  A  →  B  |     |        so the two contributions cancel.
 +-----+-----+-----+
 |     |     |     |   sum of every box's outflow
 +-----+-----+-----+   = outflow through the outer skin only
```

So the integral of divergence over a volume is the flux through its boundary: the divergence theorem, guessed. Tile a surface with tiny loops instead: neighbours traverse a shared edge in opposite directions, interior edges cancel, and the total circulation is the circulation around the rim: Stokes' theorem. Both are $\int_a^b f'\,dx = f(b) - f(a)$ in a bigger room.

### The paradox, and the guess

Add up boxes for field (iv) over a ball centred at the origin. Every box away from the origin contributes zero; the boundary carries $4\pi$. So the box containing the origin contributes $4\pi$, however small. $\nabla\cdot(\hat{\mathbf r}/r^2)$ is "zero except at the origin, infinite there, with integral $4\pi$": guess $4\pi\delta^3(\mathbf r)$. The vortex is the 2D version, with curl $2\pi\,\delta(x)\delta(y)\,\hat{\mathbf z}$.

### Sources and swirls fix the field

Suppose you know the divergence (sources) and curl (swirls) everywhere, and the field dies off far away. Each bit of source sprays a Coulomb-like field, each bit of swirl drives a vortex-like circulation; add them up. There seems to be nothing else the field could be. That is Helmholtz's theorem, and why Maxwell's equations, which give exactly the divergence and curl of $\mathbf E$ and $\mathbf B$, determine the fields.

## 3. Formal development

### 3.1 Vector algebra in index notation

Dot product: $\mathbf A\cdot\mathbf B = A_iB_i = AB\cos\theta$. Cross product: $(\mathbf A\times\mathbf B)_i = \varepsilon_{ijk}A_jB_k$ (LA-12 (3.23)), of magnitude $AB\sin\theta$ (the area of the parallelogram on $\mathbf A$, $\mathbf B$), direction by the right-hand rule. **Triple product:** $\mathbf A\cdot(\mathbf B\times\mathbf C) = \varepsilon_{ijk}A_iB_jC_k = \det[\mathbf A\ \mathbf B\ \mathbf C]$, the signed volume of the parallelepiped on the three vectors (LA-6 Theorem 3.15(a)); unchanged by cyclic permutation, negated by a swap.

**Theorem 3.1 ($\varepsilon$–$\delta$ identity).**
$$\varepsilon_{ijk}\varepsilon_{imn} = \delta_{jm}\delta_{kn} - \delta_{jn}\delta_{km}. \tag{3.1}$$

*Proof.* Fix $j, k, m, n$. A term of the sum over $i$ is nonzero only if $i, j, k$ are distinct and $i, m, n$ are distinct, forcing $j \ne k$ and $\{m, n\} = \{j, k\}$. Then the sum is $\varepsilon_{ijk}^2 = 1$ if $(m, n) = (j, k)$ and $\varepsilon_{ijk}\varepsilon_{ikj} = -1$ if $(m, n) = (k, j)$. The right side takes the same values, and for $j = k$ its two terms cancel. ∎

Every identity in this module is (3.1) plus bookkeeping. First, **BAC-CAB**: cycling $\varepsilon_{ijk} = \varepsilon_{kij}$ to put the summed index first,
$(\mathbf A\times(\mathbf B\times\mathbf C))_i = \varepsilon_{kij}\varepsilon_{kmn}A_jB_mC_n = (\delta_{im}\delta_{jn} - \delta_{in}\delta_{jm})A_jB_mC_n = B_iA_jC_j - C_iA_jB_j$, so
$$\mathbf A\times(\mathbf B\times\mathbf C) = \mathbf B(\mathbf A\cdot\mathbf C) - \mathbf C(\mathbf A\cdot\mathbf B). \tag{3.2}$$
The cross product is not associative; Exercise A2 shows $(\mathbf A\times\mathbf B)\times\mathbf C \ne \mathbf A\times(\mathbf B\times\mathbf C)$ with numbers.

**Position and separation vectors.** $\mathbf r = (x, y, z) = r\hat{\mathbf r}$. In E&M the field at $\mathbf r$ from a source at $\mathbf r'$ depends on $\boldsymbol{\mathscr r} = \mathbf r - \mathbf r'$. Two facts you will use weekly, with $\nabla$ acting on $\mathbf r$:
$$\nabla\mathscr r = \hat{\boldsymbol{\mathscr r}}, \qquad \nabla\frac{1}{\mathscr r} = -\frac{\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2}.$$
(Proof: $\partial_x\sqrt{(x-x')^2 + (y-y')^2 + (z-z')^2} = (x-x')/\mathscr r$; then the chain rule.) Since $\boldsymbol{\mathscr r}$ depends only on $\mathbf r - \mathbf r'$, on any function of $\boldsymbol{\mathscr r}$ the source-point gradient is $\nabla' = -\nabla$.

### 3.2 What makes three numbers a vector

Rotate the axes. The new unit vectors are related to the old by an orthogonal matrix $R_{ij} = \hat{\mathbf e}'_i\cdot\hat{\mathbf e}_j$, $R^\top R = I$, $\det R = +1$ (LA-12 Appendix A.1). Displacements transform as $x'_i = R_{ij}x_j$. **A vector is a triple that transforms the same way:**
$$A'_i = R_{ij}A_j. \tag{3.3}$$
(Temperature, pressure, density) is a triple but not a vector: rotating the axes changes none of its entries. LA-12 §3.17 has the general law; we need three consequences.

1. **Dot products are scalars:** $A'_iB'_i = R_{ij}R_{ik}A_jB_k = (R^\top R)_{jk}A_jB_k = A_jB_j$.
2. **$\nabla$ transforms as a vector.** Inverting, $x_j = R_{ij}x'_i$, so $\partial x_j/\partial x'_i = R_{ij}$ and the chain rule gives $\partial'_i = R_{ij}\partial_j$. So $\nabla f$ is a vector and, by 1, $\nabla\cdot\mathbf F$ is a scalar.
3. **Cross products are vectors under rotations but not reflections.** The determinant identity $\varepsilon_{lmn}R_{il}R_{jm}R_{kn} = (\det R)\,\varepsilon_{ijk}$ (the Leibniz formula, LA-6 Theorem 3.9, applied to the matrix whose rows are rows $i, j, k$ of $R$; Exercise C4(a) proves it) gives $(\mathbf A'\times\mathbf B')_i = (\det R)\,R_{il}(\mathbf A\times\mathbf B)_l$. For rotations this is (3.3). A reflection ($\det R = -1$) adds a minus sign: $\mathbf A\times\mathbf B$, $\nabla\times\mathbf F$, and the magnetic field of EM-5 are **pseudovectors**. Exercise C4 does the curl.

Theorem 3.4 proves a fact about the $z$-component of the curl and extends it to every direction by rotating the axes; that step needs 3.

### 3.3 The gradient

**Theorem 3.2 (Gradient).** Let $f$ be differentiable and $\nabla f = (\partial_xf, \partial_yf, \partial_zf)$. For a displacement $d\mathbf l = (dx, dy, dz)$,
$$df = \nabla f\cdot d\mathbf l. \tag{3.4}$$
Hence the rate of change of $f$ along a unit vector $\hat{\mathbf u}$ is $\nabla f\cdot\hat{\mathbf u}$; it is largest, equal to $|\nabla f|$, when $\hat{\mathbf u}$ points along $\nabla f$; and $\nabla f$ is perpendicular to the level surfaces of $f$.

▶ **Try it first.** Write the multivariable chain rule for $df$ and recognize a dot product. For steepest ascent, use $\mathbf a\cdot\mathbf b = ab\cos\theta$.

**Derivation.**
1. Differentiability means $f(\mathbf r + d\mathbf l) - f(\mathbf r) = \partial_xf\,dx + \partial_yf\,dy + \partial_zf\,dz$ to first order. That sum is $\nabla f\cdot d\mathbf l$.
2. Put $d\mathbf l = \hat{\mathbf u}\,dt$ and divide by $dt$: $df/dt = \nabla f\cdot\hat{\mathbf u} = |\nabla f|\cos\theta$, $\theta$ the angle between them. Largest at $\theta = 0$.
3. If $d\mathbf l$ is tangent to a level surface, $df = 0$, so $\nabla f\cdot d\mathbf l = 0$ for every tangent direction. ∎

■ **What just happened.** Keep (3.4) as the definition: the gradient is *whatever you dot with $d\mathbf l$ to get $df$*. §3.10 uses it to get the spherical gradient without transforming a derivative.

### 3.4 Divergence and curl

For a differentiable vector field $\mathbf F$,
$$\nabla\cdot\mathbf F = \partial_iF_i = \partial_xF_x + \partial_yF_y + \partial_zF_z, \tag{3.5}$$
$$(\nabla\times\mathbf F)_i = \varepsilon_{ijk}\partial_jF_k: \qquad \nabla\times\mathbf F = (\partial_yF_z - \partial_zF_y,\ \partial_zF_x - \partial_xF_z,\ \partial_xF_y - \partial_yF_x). \tag{3.6}$$
These are formulas. The geometry is a theorem, and the geometry is what you will use.

**Theorem 3.3 (Divergence is flux density).** Let $\mathbf F$ be $C^1$ near $\mathbf r_0$, and let $B$ be the box centred at $\mathbf r_0$ with sides $\Delta x, \Delta y, \Delta z$ and volume $\Delta\tau = \Delta x\,\Delta y\,\Delta z$. Then
$$\nabla\cdot\mathbf F(\mathbf r_0) = \lim_{\Delta x, \Delta y, \Delta z\to0}\frac{1}{\Delta\tau}\oint_{\partial B}\mathbf F\cdot d\mathbf a. \tag{3.7}$$

▶ **Try it first.** Take the two faces perpendicular to $x$. Their outward normals are $\pm\hat{\mathbf x}$, so together they contribute $\iint[F_x(\text{right face}) - F_x(\text{left face})]\,dy\,dz$. Write the bracket as an integral of a derivative.

**Derivation.**
1. On the face $x = x_0 + \Delta x/2$ the outward normal is $+\hat{\mathbf x}$; on $x = x_0 - \Delta x/2$ it is $-\hat{\mathbf x}$. Combined flux: $\iint\left[F_x(x_0 + \tfrac{\Delta x}{2}, y, z) - F_x(x_0 - \tfrac{\Delta x}{2}, y, z)\right]dy\,dz$.
2. By the one-variable fundamental theorem the bracket is $\int\partial_xF_x\,dx$ over the box's $x$-range, so this pair contributes exactly $\iiint_B\partial_xF_x\,d\tau$.
3. Likewise for the $y$ and $z$ pairs. Adding: $\oint_{\partial B}\mathbf F\cdot d\mathbf a = \iiint_B\nabla\cdot\mathbf F\,d\tau$, exactly, for any box.
4. $\nabla\cdot\mathbf F$ is continuous, so $\iiint_B\nabla\cdot\mathbf F\,d\tau = \nabla\cdot\mathbf F(\mathbf r_0)\,\Delta\tau + o(\Delta\tau)$. Divide by $\Delta\tau$ and shrink the box. ∎

■ **What just happened.** Divergence is outflow per unit volume. Step 3 is already the divergence theorem for one box. Flux doesn't care about coordinates, so (3.7) is a coordinate-free definition; §3.10 uses it that way.

**Theorem 3.4 (Curl is circulation density).** Let $\mathbf F$ be $C^1$ near $\mathbf r_0$, $\hat{\mathbf n}$ a unit vector, and $P$ a small planar loop around $\mathbf r_0$ in the plane perpendicular to $\hat{\mathbf n}$, enclosing area $\Delta a$, traversed counterclockwise as seen from the tip of $\hat{\mathbf n}$. Then
$$\hat{\mathbf n}\cdot(\nabla\times\mathbf F)(\mathbf r_0) = \lim_{\Delta a\to0}\frac{1}{\Delta a}\oint_P\mathbf F\cdot d\mathbf l. \tag{3.8}$$

▶ **Try it first.** Take $\hat{\mathbf n} = \hat{\mathbf z}$ and $P$ a rectangle with sides $\Delta x$, $\Delta y$. Pair the bottom edge with the top, the left with the right.

**Derivation.**
1. *Rectangle in the $xy$-plane.* Counterclockwise from above, the bottom edge ($y = y_0 - \Delta y/2$) runs in $+x$ and the top edge in $-x$. Together: $\int\left[F_x(x, y_0 - \tfrac{\Delta y}{2}) - F_x(x, y_0 + \tfrac{\Delta y}{2})\right]dx = -\iint\partial_yF_x\,dy\,dx$.
2. The right edge runs in $+y$, the left in $-y$: $\int\left[F_y(x_0 + \tfrac{\Delta x}{2}, y) - F_y(x_0 - \tfrac{\Delta x}{2}, y)\right]dy = \iint\partial_xF_y\,dx\,dy$.
3. Sum: $\oint_P\mathbf F\cdot d\mathbf l = \iint(\partial_xF_y - \partial_yF_x)\,dx\,dy = (\nabla\times\mathbf F)_z\,\Delta a + o(\Delta a)$. (Before the limit, this is Green's theorem in the plane for a rectangle.)
4. *Any direction.* Rotate the axes so $\hat{\mathbf n}$ is the new $z'$-axis. The circulation is built from dot products, so it is axis-independent, and step 3 in the new axes gives $(\nabla'\times\mathbf F')_{z'}\Delta a = \hat{\mathbf n}\cdot(\nabla\times\mathbf F)\,\Delta a$, because the curl transforms as a vector (§3.2).
5. *Any shape.* Fill the loop with tiny rectangles; shared edges cancel (§2). ∎

■ **What just happened.** The curl's component along $\hat{\mathbf n}$ is circulation per unit area around $\hat{\mathbf n}$; a paddlewheel with axle $\hat{\mathbf n}$ turns at $\tfrac12\hat{\mathbf n}\cdot(\nabla\times\mathbf F)$ (the whirlpool: curl 2, rate 1). The project shrinks a tilted parallelogram by factors of 10 and the error in (3.8) falls by about 100 each time.

### 3.5 Product rules

| | |
|---|---|
| (P1) | $\nabla(fg) = f\nabla g + g\nabla f$ |
| (P2) | $\nabla(\mathbf A\cdot\mathbf B) = \mathbf A\times(\nabla\times\mathbf B) + \mathbf B\times(\nabla\times\mathbf A) + (\mathbf A\cdot\nabla)\mathbf B + (\mathbf B\cdot\nabla)\mathbf A$ |
| (P3) | $\nabla\cdot(f\mathbf A) = f\,\nabla\cdot\mathbf A + \mathbf A\cdot\nabla f$ |
| (P4) | $\nabla\cdot(\mathbf A\times\mathbf B) = \mathbf B\cdot(\nabla\times\mathbf A) - \mathbf A\cdot(\nabla\times\mathbf B)$ |
| (P5) | $\nabla\times(f\mathbf A) = f\,\nabla\times\mathbf A - \mathbf A\times\nabla f$ |
| (P6) | $\nabla\times(\mathbf A\times\mathbf B) = (\mathbf B\cdot\nabla)\mathbf A - (\mathbf A\cdot\nabla)\mathbf B + \mathbf A(\nabla\cdot\mathbf B) - \mathbf B(\nabla\cdot\mathbf A)$ |

Here $(\mathbf A\cdot\nabla)\mathbf B$ has components $A_j\partial_jB_i$. Don't memorize these; derive the one you need in index notation.

▶ **Try it first.** Prove (P3) and (P4). For (P4) you need one cyclic shift of $\varepsilon$ and one swap.

**Derivation of (P3).** $\partial_i(fA_i) = f\,\partial_iA_i + A_i\,\partial_if$, the one-variable product rule. ∎

**Derivation of (P4).**
1. $\partial_i(\varepsilon_{ijk}A_jB_k) = \varepsilon_{ijk}(\partial_iA_j)B_k + \varepsilon_{ijk}A_j(\partial_iB_k)$.
2. First term: $\varepsilon_{ijk} = \varepsilon_{kij}$, so it is $B_k\,\varepsilon_{kij}\partial_iA_j = \mathbf B\cdot(\nabla\times\mathbf A)$.
3. Second term: $\varepsilon_{ijk} = -\varepsilon_{jik}$, so it is $-A_j\,\varepsilon_{jik}\partial_iB_k = -\mathbf A\cdot(\nabla\times\mathbf B)$. ∎

■ **What just happened.** The minus sign in (P4) is the antisymmetry of $\varepsilon$ and nothing else. Exercise B5 does (P5) and (P6); (P6) needs (3.1). One more, used in CM-3's Lagrangian for a charged particle: for a constant vector $\mathbf v$, (3.1) gives $[\mathbf v\times(\nabla\times\mathbf A)]_i = \varepsilon_{kij}\varepsilon_{klm}v_j\partial_lA_m = v_j\partial_iA_j - v_j\partial_jA_i$.

### 3.6 Second derivatives

Applying $\nabla$ twice gives five combinations. The **Laplacian** of a scalar is
$$\nabla^2f = \nabla\cdot(\nabla f) = \partial_i\partial_if, \tag{3.9}$$
and of a vector field, *in Cartesian components only*,
$$\nabla^2\mathbf A = (\nabla^2A_x,\ \nabla^2A_y,\ \nabla^2A_z). \tag{3.10}$$
In curvilinear coordinates $\nabla^2\mathbf A$ is *not* the Laplacian of each component, because the unit vectors change from point to point; (3.13) below is the safe definition. $\nabla(\nabla\cdot\mathbf A)$ has no name. The remaining two always vanish.

**Theorem 3.5.** For $f$ and $\mathbf A$ with continuous second derivatives,
$$\nabla\times(\nabla f) = \mathbf 0, \tag{3.11}$$
$$\nabla\cdot(\nabla\times\mathbf A) = 0. \tag{3.12}$$

▶ **Try it first.** Write each in index notation. You get $\varepsilon$ contracted with a pair of derivatives. What symmetry does each factor have?

**Derivation.**
1. *Lemma.* If $S_{jk} = S_{kj}$ then $\varepsilon_{ijk}S_{jk} = 0$. Rename the dummy indices $j \leftrightarrow k$, then use both symmetries: $\varepsilon_{ijk}S_{jk} = \varepsilon_{ikj}S_{kj} = -\varepsilon_{ijk}S_{jk}$. A number equal to its own negative is zero.
2. $(\nabla\times\nabla f)_i = \varepsilon_{ijk}\partial_j\partial_kf$. Mixed partials commute (Clairaut, which needs continuous second derivatives), so $\partial_j\partial_kf$ is symmetric. Zero by the lemma.
3. $\nabla\cdot(\nabla\times\mathbf A) = \partial_i\varepsilon_{ijk}\partial_jA_k = \varepsilon_{kij}(\partial_i\partial_jA_k)$, and for each $k$ the bracket is symmetric in $i, j$. Zero. ∎

■ **What just happened.** Both are "symmetric contracted with antisymmetric is zero." §3.13 proves the converses on suitable regions.

**Theorem 3.6 (Curl of a curl).** For $\mathbf A$ with continuous second derivatives,
$$\nabla\times(\nabla\times\mathbf A) = \nabla(\nabla\cdot\mathbf A) - \nabla^2\mathbf A. \tag{3.13}$$

▶ **Try it first.** This is BAC-CAB (3.2) with $\mathbf A, \mathbf B \to \nabla$ and $\mathbf C\to\mathbf A$, except that the derivatives must stay to the left of $\mathbf A$. Redo the index computation of (3.2) keeping them there.

**Derivation.**
1. $[\nabla\times(\nabla\times\mathbf A)]_i = \varepsilon_{ijk}\partial_j(\varepsilon_{klm}\partial_lA_m) = \varepsilon_{kij}\varepsilon_{klm}\,\partial_j\partial_lA_m$ (the constants $\varepsilon$ pass through $\partial_j$; then cycle).
2. By (3.1): $(\delta_{il}\delta_{jm} - \delta_{im}\delta_{jl})\,\partial_j\partial_lA_m = \partial_j\partial_iA_j - \partial_j\partial_jA_i$.
3. The first term is $\partial_i(\nabla\cdot\mathbf A)$ by Clairaut; the second is $(\nabla^2\mathbf A)_i$ by (3.10). ∎

■ **What just happened.** Three lines, versus a page of components: that is what LA-12's $\varepsilon_{ijk}$ buys. The curl of Faraday's law plus (3.13) is EM-6's wave equation for light.

### 3.7 Line, surface, and volume integrals: the recipe

Most "I can't do this problem" moments in E&M are "I can't set up this integral." The recipe:

**Line integral.** Parametrize the path as $\mathbf r(t)$ from $t_a$ (start) to $t_b$ (end). Then $d\mathbf l = \mathbf r'(t)\,dt$ and
$$\int_P\mathbf F\cdot d\mathbf l = \int_{t_a}^{t_b}\mathbf F(\mathbf r(t))\cdot\mathbf r'(t)\,dt. \tag{3.14}$$
Reversing the direction flips the sign; reparametrizing changes nothing (chain rule).

**Surface integral.** Parametrize the surface as $\mathbf r(u, v)$, $(u, v)\in D$. The tangents $\mathbf r_u = \partial\mathbf r/\partial u$ and $\mathbf r_v$ span a small parallelogram of vector area $(\mathbf r_u\times\mathbf r_v)\,du\,dv$:
$$d\mathbf a = (\mathbf r_u\times\mathbf r_v)\,du\,dv, \qquad \int_S\mathbf F\cdot d\mathbf a = \iint_D\mathbf F(\mathbf r(u, v))\cdot(\mathbf r_u\times\mathbf r_v)\,du\,dv. \tag{3.15}$$
Then **check the direction** of $\mathbf r_u\times\mathbf r_v$ (outward on a closed surface, right-hand rule with the boundary on an open one); if wrong, swap $u$ and $v$. $|\mathbf r_u\times\mathbf r_v|$ is LA-6 §3.11's area-scaling Jacobian. Two you'll use constantly: the sphere $\mathbf r = R(\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$ has $\mathbf r_\theta\times\mathbf r_\phi = R^2\sin\theta\,\hat{\mathbf r}$ (outward); the disk $\mathbf r = (s\cos\phi, s\sin\phi, 0)$ has $\mathbf r_s\times\mathbf r_\phi = s\,\hat{\mathbf z}$.

**Volume integral.** With $J$ the Jacobian matrix of $(u, v, w)\mapsto\mathbf r$,
$$d\tau = |\det J|\,du\,dv\,dw. \tag{3.16}$$

**Before integrating, check** the sign of $d\mathbf a$, which terms a symmetry kills, and the units. Example 4.3 runs the recipe three times.

### 3.8 The fundamental theorems

**Theorem 3.7 (Fundamental theorem for gradients).** If $f$ is $C^1$ and $P$ is a piecewise-smooth path from $\mathbf a$ to $\mathbf b$,
$$\int_P\nabla f\cdot d\mathbf l = f(\mathbf b) - f(\mathbf a). \tag{3.17}$$
*Proof.* Along $\mathbf r(t)$, (3.4) gives $\frac{d}{dt}f(\mathbf r(t)) = \nabla f\cdot\mathbf r'(t)$. Integrate from $t_a$ to $t_b$ and compare with (3.14). ∎

**Theorem 3.8 (Divergence theorem, Gauss).** Let $V$ be a bounded region with piecewise-smooth closed boundary $S$, normal outward, and let $\mathbf F$ be $C^1$ on $V$ and $S$. Then
$$\int_V(\nabla\cdot\mathbf F)\,d\tau = \oint_S\mathbf F\cdot d\mathbf a. \tag{3.18}$$

▶ **Try it first.** Step 3 of Theorem 3.3 proved this for one box. Glue boxes together. What happens on a face two boxes share?

**Derivation.**
1. Cover $V$ with a grid of cubes of side $h$. For each cube $B_\alpha$ inside $V$, step 3 of Theorem 3.3 gives exactly $\int_{B_\alpha}\nabla\cdot\mathbf F\,d\tau = \oint_{\partial B_\alpha}\mathbf F\cdot d\mathbf a$.
2. Sum over $\alpha$. The left sides add to the integral over the union of the cubes.
3. On the right, each interior face belongs to two cubes, and its outward normal for one is the inward normal for the other: its two contributions cancel. The sum telescopes to the flux through the outer surface of the union.
4. As $h\to0$ the union fills $V$ and its surface approaches $S$; the two sides tend to the two sides of (3.18). ∎

**What is being skipped.** Step 4 near a curved $S$, where cubes are cut: a staircase's *area* does not converge to that of $S$, but its *flux* does, since $\mathbf F\cdot d\mathbf a$ sees only projected area. The rigorous proof (Marsden and Tromba, *Vector Calculus*, ch. 8) uses step 2 of Theorem 3.3 once per coordinate direction on suitably simple regions, then glues.

■ **What just happened.** The box calculation, summed, with interior faces cancelling in pairs: the §2 picture made exact.

**Theorem 3.9 (Stokes' theorem).** Let $S$ be a piecewise-smooth oriented surface with boundary curve $P$, oriented by the **right-hand rule** (fingers along $d\mathbf l$, thumb along $d\mathbf a$), and $\mathbf F$ be $C^1$ on $S$. Then
$$\int_S(\nabla\times\mathbf F)\cdot d\mathbf a = \oint_P\mathbf F\cdot d\mathbf l. \tag{3.19}$$

▶ **Try it first.** Same as Theorem 3.8, with loops for boxes and edges for faces.

**Derivation.**
1. Divide $S$ into small, nearly flat patches with vector areas $\Delta\mathbf a_\alpha = \hat{\mathbf n}_\alpha\Delta a_\alpha$; orient each patch boundary $P_\alpha$ by the right-hand rule.
2. Theorem 3.4: $\oint_{P_\alpha}\mathbf F\cdot d\mathbf l = (\nabla\times\mathbf F)\cdot\Delta\mathbf a_\alpha + o(\Delta a_\alpha)$.
3. Sum over $\alpha$. An edge shared by two neighbouring patches is traversed once in each direction (both loops run counterclockwise seen from the same side), so interior edges cancel and only $P$ survives, in its own direction.
4. The right sides form a Riemann sum for $\int_S(\nabla\times\mathbf F)\cdot d\mathbf a$; the errors total $o(1)$ times the area of $S$. Refine. ∎

■ **What just happened.** "Add up the paddlewheels." The orientation rule is part of the theorem. Any surface with boundary $P$ works (Example 4.3 checks three).

**Theorem 3.10 (Corollaries).** (a) $\int_{\mathbf a}^{\mathbf b}\nabla f\cdot d\mathbf l$ is independent of path, and around any closed loop
$$\oint\nabla f\cdot d\mathbf l = 0. \tag{3.20}$$
(b) $\int_S(\nabla\times\mathbf A)\cdot d\mathbf a$ depends only on the boundary of $S$. (c) Over any closed surface,
$$\oint_S(\nabla\times\mathbf A)\cdot d\mathbf a = 0. \tag{3.21}$$
*Proof.* (a) Theorem 3.7; for a loop $\mathbf a = \mathbf b$. (b) Every such surface gives $\oint_P\mathbf A\cdot d\mathbf l$ by Stokes. (c) Divergence theorem then (3.12): $\int_V\nabla\cdot(\nabla\times\mathbf A)\,d\tau = 0$. ∎

### 3.9 Green's identities

Integration by parts in three dimensions; the workhorses of EM-2 (field energy) and EM-3 (uniqueness theorems, Green's functions).

**Theorem 3.11 (Green's identities).** Let $\varphi, \psi$ have continuous second derivatives on $V$ and $S = \partial V$. Then
$$\int_V\left(\varphi\nabla^2\psi + \nabla\varphi\cdot\nabla\psi\right)d\tau = \oint_S\varphi\,\nabla\psi\cdot d\mathbf a, \tag{3.22}$$
$$\int_V\left(\varphi\nabla^2\psi - \psi\nabla^2\varphi\right)d\tau = \oint_S\left(\varphi\nabla\psi - \psi\nabla\varphi\right)\cdot d\mathbf a. \tag{3.23}$$

▶ **Try it first.** Apply the divergence theorem to $\mathbf F = \varphi\nabla\psi$ and expand $\nabla\cdot\mathbf F$ with (P3).

**Derivation.**
1. (P3) with $f = \varphi$, $\mathbf A = \nabla\psi$: $\nabla\cdot(\varphi\nabla\psi) = \varphi\nabla^2\psi + \nabla\varphi\cdot\nabla\psi$.
2. Integrate over $V$; apply (3.18) to the left side. That is **Green's first identity** (3.22).
3. Write (3.22) with $\varphi \leftrightarrow \psi$ and subtract. The symmetric term $\nabla\varphi\cdot\nabla\psi$ cancels: **Green's second identity** (3.23). ∎

■ **What just happened.** In 1D, (3.22) is $\int_a^b(\varphi\psi'' + \varphi'\psi')\,dx = [\varphi\psi']_a^b$. With $\varphi = \psi$ it trades $\int|\nabla\psi|^2$ for a boundary term minus $\int\psi\nabla^2\psi$: EM-2's field energy $\tfrac{\varepsilon_0}{2}\int E^2\,d\tau$, EM-3's uniqueness theorems, and Exercise C2 all run on this. (3.23) with $\psi = 1/\mathscr r$ is EM-3's route to Green's functions.

### 3.10 Curvilinear coordinates

**Spherical coordinates:**
$$x = r\sin\theta\cos\phi, \quad y = r\sin\theta\sin\phi, \quad z = r\cos\theta, \qquad r \ge 0,\ 0\le\theta\le\pi,\ 0\le\phi<2\pi. \tag{3.24}$$

▶ **Try it first.** Differentiate $\mathbf r(r, \theta, \phi)$ with respect to each coordinate. How long is each tangent vector? Are they perpendicular?

**Derivation (scale factors).**
1. $\partial\mathbf r/\partial r = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$, length 1. Call it $\hat{\mathbf r}$.
2. $\partial\mathbf r/\partial\theta = r(\cos\theta\cos\phi, \cos\theta\sin\phi, -\sin\theta)$, length $r$; its direction is $\hat{\boldsymbol\theta}$.
3. $\partial\mathbf r/\partial\phi = r\sin\theta(-\sin\phi, \cos\phi, 0)$, length $r\sin\theta$; its direction is $\hat{\boldsymbol\phi}$.
4. Direct dot products show the three are mutually orthogonal, and $\hat{\mathbf r}\times\hat{\boldsymbol\theta} = \hat{\boldsymbol\phi}$. The lengths are the **scale factors** $h_r = 1$, $h_\theta = r$, $h_\phi = r\sin\theta$, so
$$d\mathbf l = dr\,\hat{\mathbf r} + r\,d\theta\,\hat{\boldsymbol\theta} + r\sin\theta\,d\phi\,\hat{\boldsymbol\phi}. \tag{3.25}$$
5. The tangents are orthogonal, so the coordinate box is rectangular to first order with volume
$$d\tau = h_rh_\theta h_\phi\,dr\,d\theta\,d\phi = r^2\sin\theta\,dr\,d\theta\,d\phi, \tag{3.26}$$
which equals $|\det J|$ in (3.16); the project checks this symbolically. ∎

■ **What just happened.** One unit of $\theta$ moves you $r$ metres. Every formula below divides a coordinate derivative by a scale factor. The unit vectors vary ($\partial\hat{\mathbf r}/\partial\theta = \hat{\boldsymbol\theta}$), so never integrate $\hat{\mathbf r}$ as a constant.

**Theorem 3.12 (Spherical gradient).**
$$\nabla f = \frac{\partial f}{\partial r}\hat{\mathbf r} + \frac1r\frac{\partial f}{\partial\theta}\hat{\boldsymbol\theta} + \frac{1}{r\sin\theta}\frac{\partial f}{\partial\phi}\hat{\boldsymbol\phi}. \tag{3.27}$$

*Derivation.*
1. Chain rule: $df = \partial_rf\,dr + \partial_\theta f\,d\theta + \partial_\phi f\,d\phi$.
2. By (3.4) with $\nabla f = g_r\hat{\mathbf r} + g_\theta\hat{\boldsymbol\theta} + g_\phi\hat{\boldsymbol\phi}$ and (3.25): $df = g_r\,dr + g_\theta\,r\,d\theta + g_\phi\,r\sin\theta\,d\phi$.
3. Both hold for all $dr, d\theta, d\phi$; match coefficients. ∎

**Theorem 3.13 (Spherical divergence).**
$$\nabla\cdot\mathbf F = \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2F_r\right) + \frac{1}{r\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\,F_\theta\right) + \frac{1}{r\sin\theta}\frac{\partial F_\phi}{\partial\phi}. \tag{3.28}$$

▶ **Try it first.** Apply (3.7) to the box bounded by the coordinate surfaces at $r$, $r + dr$, $\theta$, $\theta + d\theta$, $\phi$, $\phi + d\phi$. The two $r$-faces have different areas; that is where the $r^2$ inside the derivative comes from.

*Derivation.*
1. Volume: $r^2\sin\theta\,dr\,d\theta\,d\phi$, by (3.26).
2. *$r$-faces.* The face at radius $r$ has area $(r\,d\theta)(r\sin\theta\,d\phi)$ and outward normal $-\hat{\mathbf r}$; the face at $r + dr$ has $+\hat{\mathbf r}$. Net flux: $\left[r^2F_r\right]_r^{r+dr}\sin\theta\,d\theta\,d\phi = \partial_r(r^2F_r)\,\sin\theta\,dr\,d\theta\,d\phi$.
3. *$\theta$-faces.* Area $dr\cdot r\sin\theta\,d\phi$, normals $\pm\hat{\boldsymbol\theta}$. Net: $\partial_\theta(\sin\theta\,F_\theta)\,r\,dr\,d\theta\,d\phi$.
4. *$\phi$-faces.* Area $dr\cdot r\,d\theta$, normals $\pm\hat{\boldsymbol\phi}$. Net: $\partial_\phi F_\phi\,r\,dr\,d\theta\,d\phi$.
5. Divide the total by the volume. ∎

■ **What just happened.** The coordinate-free definitions (3.4) and (3.7) did all the work. With general scale factors the same box gives $\nabla\cdot\mathbf F = \frac{1}{h_1h_2h_3}\sum_i\partial_i\big(\tfrac{h_1h_2h_3}{h_i}F_i\big)$.

Composing (3.27) and (3.28) gives the **spherical Laplacian**:
$$\nabla^2f = \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial f}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial f}{\partial\theta}\right) + \frac{1}{r^2\sin^2\theta}\frac{\partial^2f}{\partial\phi^2}. \tag{3.29}$$
For $f = 1/r$: $r^2\partial_r(1/r) = -1$, whose $r$-derivative is $0$. So **$\nabla^2(1/r) = 0$ for $r \ne 0$**: the flux of $\nabla(1/r) = -\hat{\mathbf r}/r^2$ through a sphere is $-4\pi$ whatever the radius, so none is created between two spheres. Hold that thought.

**Cylindrical coordinates** $x = s\cos\phi$, $y = s\sin\phi$, $z$: the same steps give scale factors $(1, s, 1)$,
$$d\mathbf l = ds\,\hat{\mathbf s} + s\,d\phi\,\hat{\boldsymbol\phi} + dz\,\hat{\mathbf z}, \qquad d\tau = s\,ds\,d\phi\,dz. \tag{3.30}$$

**The table.** The curls come from (3.8) on small loops in coordinate surfaces, as (3.28) came from (3.7): a page of bookkeeping with nothing new, so we state them (as on Griffiths' inside front cover; the project checks every entry against Cartesian):

| | Spherical $(r, \theta, \phi)$ | Cylindrical $(s, \phi, z)$ |
|---|---|---|
| $\nabla f$ | (3.27) | $\partial_sf\,\hat{\mathbf s} + \frac1s\partial_\phi f\,\hat{\boldsymbol\phi} + \partial_zf\,\hat{\mathbf z}$ |
| $\nabla\cdot\mathbf F$ | (3.28) | $\frac1s\partial_s(sF_s) + \frac1s\partial_\phi F_\phi + \partial_zF_z$ |
| $\nabla\times\mathbf F$ | $\frac{\partial_\theta(\sin\theta F_\phi) - \partial_\phi F_\theta}{r\sin\theta}\hat{\mathbf r} + \frac1r\Big[\frac{\partial_\phi F_r}{\sin\theta} - \partial_r(rF_\phi)\Big]\hat{\boldsymbol\theta} + \frac{\partial_r(rF_\theta) - \partial_\theta F_r}{r}\hat{\boldsymbol\phi}$ | $\Big[\frac{\partial_\phi F_z}{s} - \partial_zF_\phi\Big]\hat{\mathbf s} + \big[\partial_zF_s - \partial_sF_z\big]\hat{\boldsymbol\phi} + \frac{\partial_s(sF_\phi) - \partial_\phi F_s}{s}\hat{\mathbf z}$ |
| $\nabla^2f$ | (3.29) | $\frac1s\partial_s(s\,\partial_sf) + \frac{1}{s^2}\partial_\phi^2f + \partial_z^2f$ |

Do not memorize the spherical curl. Look it up every time.

### 3.11 The Dirac delta function in three dimensions

**The paradox, precisely.** Let $\mathbf v = \hat{\mathbf r}/r^2$. By (3.28), $\nabla\cdot\mathbf v = \frac{1}{r^2}\partial_r(1) = 0$ for $r \ne 0$. By (3.15), the flux through a sphere of radius $R$ is $\int\frac{1}{R^2}R^2\sin\theta\,d\theta\,d\phi = 4\pi$. The divergence theorem would give $\int_V0\,d\tau = 4\pi$. Theorem 3.8 has not failed (its hypothesis fails at the origin), but physics needs an answer: $\mathbf v$ is a point charge's field up to $q/4\pi\varepsilon_0$.

**The 3D delta.** From F3 (Definition 3.12, Theorem 3.14), $\delta(x)$ is not a function but a rule, $\int f(x)\,\delta(x - a)\,dx = f(a)$, and the weak limit of narrowing unit-area bumps. In 3D, $\delta^3(\mathbf r) = \delta(x)\delta(y)\delta(z)$, with the **sifting property**
$$\int_Vf(\mathbf r)\,\delta^3(\mathbf r - \mathbf a)\,d\tau = \begin{cases}f(\mathbf a), & \mathbf a \text{ inside } V,\\ 0, & \mathbf a \text{ outside } V.\end{cases} \tag{3.31}$$
So $\int\delta^3\,d\tau = 1$ over all space, $\delta^3$ has units of 1/volume, and substituting $\mathbf u = c\mathbf r$ in (3.31) gives $\delta^3(c\mathbf r) = \delta^3(\mathbf r)/|c|^3$.

**Theorem 3.14 (The delta function in $\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$).**
$$\nabla\cdot\frac{\hat{\boldsymbol{\mathscr r}}}{\mathscr r^2} = 4\pi\delta^3(\boldsymbol{\mathscr r}), \tag{3.32}$$
$$\nabla^2\frac{1}{\mathscr r} = -4\pi\delta^3(\boldsymbol{\mathscr r}). \tag{3.33}$$
These are equations between distributions. (3.32) means: for every smooth test function $\psi$ that vanishes outside a bounded region, $-\int(\hat{\boldsymbol{\mathscr r}}/\mathscr r^2)\cdot\nabla\psi\,d\tau = 4\pi\psi(\mathbf r')$. The left side is what $\int\psi\,\nabla\cdot\mathbf v\,d\tau$ becomes after one integration by parts, which is how F3 §3.10 defines the derivative of a distribution.

▶ **Try it first.** Put the source at the origin. Cut out a small ball of radius $\epsilon$. Outside it everything is smooth and $\nabla\cdot\mathbf v = 0$. Apply the divergence theorem to $\psi\mathbf v$ on the region outside the ball, then let $\epsilon\to0$.

**Derivation.**
1. Take $\mathbf r' = \mathbf 0$. Let $V_\epsilon = \{\epsilon < r < R\}$, with $R$ so large that $\psi = 0$ for $r \ge R$.
2. On $V_\epsilon$, (P3) gives $\nabla\cdot(\psi\mathbf v) = \psi\,\nabla\cdot\mathbf v + \mathbf v\cdot\nabla\psi = \mathbf v\cdot\nabla\psi$.
3. Divergence theorem on $V_\epsilon$. Its boundary is two spheres. On $r = R$, $\psi = 0$. On $r = \epsilon$ the outward normal *of $V_\epsilon$* points toward the origin: $d\mathbf a = -\hat{\mathbf r}\,\epsilon^2\sin\theta\,d\theta\,d\phi$. So
$$\int_{V_\epsilon}\mathbf v\cdot\nabla\psi\,d\tau = -\oint_{r=\epsilon}\psi\,\frac{1}{\epsilon^2}\,\epsilon^2\sin\theta\,d\theta\,d\phi = -\oint_{r = \epsilon}\psi\,d\Omega,$$
with $d\Omega = \sin\theta\,d\theta\,d\phi$ the solid-angle element.
4. Let $\epsilon\to0$. The right side tends to $-4\pi\psi(\mathbf 0)$: $\psi$ is continuous and the total solid angle is $4\pi$. The left side tends to the integral over all space, since $|\mathbf v|\,d\tau = \sin\theta\,dr\,d\theta\,d\phi$ is bounded near the origin. So $-\int\mathbf v\cdot\nabla\psi\,d\tau = 4\pi\psi(\mathbf 0)$, which is (3.32).
5. A source at $\mathbf r'$: shift the origin; nothing above depended on where it was.
6. (3.33): $\nabla(1/\mathscr r) = -\hat{\boldsymbol{\mathscr r}}/\mathscr r^2$ (§3.1). Take the divergence and use (3.32). ∎

■ **What just happened.** The flux of $\hat{\mathbf r}/r^2$ through a sphere is $4\pi$ however small the sphere, and the divergence theorem turns that into (3.32). Working shorthand: $\int_V\nabla\cdot(\hat{\boldsymbol{\mathscr r}}/\mathscr r^2)\,d\tau = 4\pi$ if the source is inside $V$, $0$ if outside. And (3.33) says $G = 1/(4\pi\mathscr r)$ satisfies $-\nabla^2G = \delta^3(\mathbf r - \mathbf r')$: it is the **Green's function** of $-\nabla^2$ (the spatial cousin of F3 Theorem 3.17's impulse response), so $u = \int G\rho\,d\tau'$ solves $-\nabla^2u = \rho$. Example 4.5 watches a shrinking charge cloud's Laplacian become $-4\pi\delta^3$.

### 3.12 The Helmholtz theorem

**Theorem 3.15 (Helmholtz).** Let $D$ be a scalar field and $\mathbf C$ a vector field with $\nabla\cdot\mathbf C = 0$, both smooth and falling off faster than $1/r^2$ as $r\to\infty$. Then there is exactly one vector field $\mathbf F$ with
$$\nabla\cdot\mathbf F = D, \qquad \nabla\times\mathbf F = \mathbf C, \qquad \mathbf F\to\mathbf 0 \text{ as } r\to\infty,$$
namely
$$\mathbf F = -\nabla U + \nabla\times\mathbf W, \qquad U(\mathbf r) = \frac{1}{4\pi}\int\frac{D(\mathbf r')}{\mathscr r}\,d\tau', \qquad \mathbf W(\mathbf r) = \frac{1}{4\pi}\int\frac{\mathbf C(\mathbf r')}{\mathscr r}\,d\tau'. \tag{3.34}$$
The condition $\nabla\cdot\mathbf C = 0$ is not optional: by (3.12) the curl of anything is divergence-free.

▶ **Try it first.** Take the divergence of (3.34). You need (3.12) and (3.33), and the fact that $\nabla$ (acting on $\mathbf r$) passes inside the $\mathbf r'$ integral.

**Proof sketch.**
1. *Divergence.* $\nabla\cdot(\nabla\times\mathbf W) = 0$ by (3.12), so $\nabla\cdot\mathbf F = -\nabla^2U = -\frac{1}{4\pi}\int D(\mathbf r')\,\nabla^2\frac{1}{\mathscr r}\,d\tau' = \int D(\mathbf r')\,\delta^3(\mathbf r - \mathbf r')\,d\tau' = D(\mathbf r)$, by (3.33) and (3.31).
2. *Curl.* $\nabla\times\nabla U = \mathbf 0$ by (3.11), so by (3.13) $\nabla\times\mathbf F = \nabla(\nabla\cdot\mathbf W) - \nabla^2\mathbf W$. Each Cartesian component of $\mathbf W$ has the form of $U$, so step 1 gives $-\nabla^2\mathbf W = \mathbf C$.
3. *Gauge term.* $\nabla\cdot\mathbf W = \frac{1}{4\pi}\int\mathbf C(\mathbf r')\cdot\nabla\frac{1}{\mathscr r}\,d\tau' = -\frac{1}{4\pi}\int\mathbf C\cdot\nabla'\frac{1}{\mathscr r}\,d\tau'$ (since $\nabla = -\nabla'$ on functions of $\boldsymbol{\mathscr r}$). By (P3) and the divergence theorem in $\mathbf r'$ this equals $\frac{1}{4\pi}\int\frac{\nabla'\cdot\mathbf C}{\mathscr r}\,d\tau' - \frac{1}{4\pi}\oint\frac{\mathbf C}{\mathscr r}\cdot d\mathbf a'$. The volume term vanishes because $\nabla'\cdot\mathbf C = 0$; the surface term, on a sphere of radius $R'\to\infty$, is $o(R'^{-2})\cdot R'^{-1}\cdot R'^2\to0$. So $\nabla\cdot\mathbf W = 0$ and $\nabla\times\mathbf F = \mathbf C$.
4. *Uniqueness.* If $\mathbf F_1$ and $\mathbf F_2$ both qualify, $\mathbf G = \mathbf F_1 - \mathbf F_2$ has zero divergence and curl, so by (3.13) $\nabla^2\mathbf G = \mathbf 0$: each Cartesian component is harmonic and tends to $0$ at infinity. A harmonic function has no interior maximum or minimum (EM-3 derives this from the mean-value property), so each component is squeezed to $0$. Exercise C2 gives a self-contained proof from Green's first identity when $\mathbf G$ falls off like $1/r^2$. ∎

(Skipped: justifying differentiation under the integral sign; see Griffiths Appendix B.)

■ **What just happened.** The divergence fixes a "Coulomb" potential $U$ of the sources, the curl a "Biot–Savart" potential $\mathbf W$ of the swirls, and vanishing at infinity rules out anything else. Maxwell's equations give a divergence and a curl for each of $\mathbf E$ and $\mathbf B$, which is exactly enough. The compatibility condition $\nabla\cdot\mathbf C = 0$ applied to $\nabla\times\mathbf B = \mu_0\mathbf J + \mu_0\varepsilon_0\partial_t\mathbf E$ forces $\nabla\cdot\mathbf J + \partial_t\rho = 0$: charge conservation is built in (EM-6).

**Helmholtz in Fourier space.** On a periodic box write $\mathbf F(\mathbf r) = \sum_{\mathbf k}\hat{\mathbf F}_{\mathbf k}e^{i\mathbf k\cdot\mathbf r}$. Then $\nabla$ acts on each mode as multiplication by $i\mathbf k$ (F3's derivative rule (3.8) with angular wavevector $\mathbf k = 2\pi\boldsymbol\xi$), so $\nabla\cdot \to i\mathbf k\cdot$ and $\nabla\times \to i\mathbf k\times$: calculus becomes algebra. For $\mathbf k \ne \mathbf 0$, with $\hat{\mathbf k} = \mathbf k/k$,
$$\hat{\mathbf F}_{\mathbf k} = \underbrace{\hat{\mathbf k}(\hat{\mathbf k}\cdot\hat{\mathbf F}_{\mathbf k})}_{\text{curl-free}} + \underbrace{\big[\hat{\mathbf F}_{\mathbf k} - \hat{\mathbf k}(\hat{\mathbf k}\cdot\hat{\mathbf F}_{\mathbf k})\big]}_{\text{divergence-free}}, \qquad \hat U_{\mathbf k} = \frac{i\,\mathbf k\cdot\hat{\mathbf F}_{\mathbf k}}{k^2}, \qquad \hat{\mathbf W}_{\mathbf k} = \frac{i\,\mathbf k\times\hat{\mathbf F}_{\mathbf k}}{k^2}. \tag{3.35}$$
Checks: $i\mathbf k\times\hat{\mathbf k}(\cdots) = \mathbf 0$; $i\mathbf k\cdot[\hat{\mathbf F} - \hat{\mathbf k}(\hat{\mathbf k}\cdot\hat{\mathbf F})] = 0$; $-i\mathbf k\,\hat U = \hat{\mathbf k}(\hat{\mathbf k}\cdot\hat{\mathbf F})$; by (3.2), $i\mathbf k\times\hat{\mathbf W} = -\mathbf k\times(\mathbf k\times\hat{\mathbf F})/k^2 = \hat{\mathbf F} - \hat{\mathbf k}(\hat{\mathbf k}\cdot\hat{\mathbf F})$; and $\mathbf k\cdot\hat{\mathbf W} = 0$, so $\nabla\cdot\mathbf W = 0$. The two parts are perpendicular mode by mode, hence orthogonal as fields (Parseval). The $\mathbf k = \mathbf 0$ mode (the box average) is curl- and divergence-free and belongs to neither potential: it is the extra harmonic piece that "$\mathbf F\to\mathbf 0$ at infinity" excludes on $\mathbb R^3$. The project runs (3.35) through the FFT; Example 4.7 does it by hand.

### 3.13 Potentials

**Theorem 3.16 (Curl-free fields are gradients).** Let $\mathbf F$ be $C^1$ on an open connected region $D$. Consider: (a) $\oint_P\mathbf F\cdot d\mathbf l = 0$ for every closed loop $P$ in $D$ (equivalently, line integrals are path-independent); (b) $\nabla\times\mathbf F = \mathbf 0$ in $D$; (c) $\mathbf F = -\nabla U$ for some scalar $U$. Then (a) $\Leftrightarrow$ (c) and (c) $\Rightarrow$ (b) on any region, and (b) $\Rightarrow$ (a) when $D$ is **simply connected**.

These are CM-1 Theorem 3.3's conditions ($V$ renamed $U$); CM-1 proved all but (b) $\Rightarrow$ (a).

▶ **Try it first.** Given a loop in $D$, find a surface in $D$ that it bounds, and apply Stokes.

**Derivation.**
1. (c) $\Rightarrow$ (a) is (3.20). (a) $\Rightarrow$ (c): define $U(\mathbf r) = -\int_{\mathbf r_0}^{\mathbf r}\mathbf F\cdot d\mathbf l$, well defined by (a), and differentiate along a short straight step in $x$ to get $\partial_xU = -F_x$, likewise $y, z$ (CM-1 §3.5). (c) $\Rightarrow$ (b) is (3.11).
2. (b) $\Rightarrow$ (a). Suppose first that $D$ is **star-shaped** about $\mathbf r_0$: every segment from $\mathbf r_0$ to a point of $D$ lies in $D$ (convex regions, and all of $\mathbb R^3$, qualify). Given a loop $P$ in $D$, let $S$ be the cone of segments from $\mathbf r_0$ to the points of $P$. Then $S\subset D$ and $\partial S = P$.
3. Stokes (3.19) on $S$: $\oint_P\mathbf F\cdot d\mathbf l = \int_S(\nabla\times\mathbf F)\cdot d\mathbf a = 0$ by (b). ∎

**What is being skipped.** For simply connected but not star-shaped $D$, contracting $P$ to a point sweeps out a spanning surface; making it smooth enough for Stokes is topology we don't need. Holes are where the theorem fails: the vortex $\hat{\boldsymbol\phi}/s$ on $\mathbb R^3$ minus the $z$-axis has zero curl but circulation $2\pi$ (Exercise B2; CM-1 Exercise B2), and every cone spanning that loop is pierced by the axis.

■ **What just happened.** The useful direction is (b) $\Rightarrow$ (c): compute a curl to learn whether a potential exists, then build $U$ by a line integral along the easiest path (Example 4.6). EM-2 gets the electrostatic potential this way.

**Theorem 3.17 (Divergence-free fields are curls).** If $\mathbf F$ is $C^1$ on $\mathbb R^3$ (or on a box containing the origin) and $\nabla\cdot\mathbf F = 0$, then $\mathbf F = \nabla\times\mathbf W$ with
$$\mathbf W = \left(\int_0^zF_y(x, y, z')\,dz',\ \ \int_0^xF_z(x', y, 0)\,dx' - \int_0^zF_x(x, y, z')\,dz',\ \ 0\right). \tag{3.36}$$

▶ **Try it first.** Look for $\mathbf W$ with $W_z = 0$. Then $\nabla\times\mathbf W = (-\partial_zW_y,\ \partial_zW_x,\ \partial_xW_y - \partial_yW_x)$. Get the first two components by integrating in $z$, then repair the third with a $z$-independent term.

**Derivation.**
1. $-\partial_zW_y = F_x(x, y, z)$: the $x'$-integral doesn't depend on $z$, and the fundamental theorem handles the $z'$-integral.
2. $\partial_zW_x = F_y(x, y, z)$, the same way.
3. $\partial_xW_y - \partial_yW_x = F_z(x, y, 0) - \int_0^z(\partial_xF_x + \partial_yF_y)(x, y, z')\,dz' = F_z(x, y, 0) + \int_0^z\partial_{z'}F_z\,dz' = F_z(x, y, z)$, using $\nabla\cdot\mathbf F = 0$ in the middle. ∎

■ **What just happened.** With (3.12): on $\mathbb R^3$, divergence-free $\Leftrightarrow$ a curl $\Leftrightarrow$ zero flux through every closed surface. That is why $\nabla\cdot\mathbf B = 0$ gives $\mathbf B = \nabla\times\mathbf A$ in EM-5.

**Gauge freedom.** Potentials are not unique even when $\mathbf F$ is: $U \to U + c$ leaves $-\nabla U$ alone, and $\mathbf W\to\mathbf W + \nabla\lambda$ leaves $\nabla\times\mathbf W$ alone by (3.11). Choosing $\lambda$ with $\nabla^2\lambda = -\nabla\cdot\mathbf W$ (Poisson, solved by (3.33)) makes $\nabla\cdot\mathbf W = 0$: the **Coulomb gauge**. The $\mathbf W$ of (3.34) and (3.35) is already in it; that of (3.36) usually is not (Exercise B7 repairs one). EM-5 and EM-6 choose gauges on purpose; in QM the freedom becomes local phase invariance.

## 4. Worked examples

### Example 4.1: Div and curl of the four fields of §2

*(i) $\mathbf r = (x, y, z)$.* $\nabla\cdot\mathbf r = 1 + 1 + 1 = 3$. $\nabla\times\mathbf r = (\partial_yz - \partial_zy, \ldots) = \mathbf 0$. Uniform source density, no spin.

*(ii) Whirlpool $(-y, x, 0)$.* Divergence $\partial_x(-y) + \partial_y(x) = 0$. The curl's $x$, $y$ components need $\partial_z$ of $z$-free things: $0$. Its $z$ component is $\partial_x(x) - \partial_y(-y) = 2$. This is rigid rotation $\boldsymbol\omega\times\mathbf r$ with $\boldsymbol\omega = \hat{\mathbf z}$; in general (P6) gives $\nabla\times(\boldsymbol\omega\times\mathbf r) = \boldsymbol\omega(\nabla\cdot\mathbf r) - (\boldsymbol\omega\cdot\nabla)\mathbf r = 2\boldsymbol\omega$.

*(iii) Vortex $(-y, x, 0)/s^2$, $s^2 = x^2 + y^2$.* $\partial_x(-y/s^2) = 2xy/s^4$ and $\partial_y(x/s^2) = -2xy/s^4$: divergence $0$. Curl $z$-component: $\partial_x(x/s^2) = (s^2 - 2x^2)/s^4$, $\partial_y(-y/s^2) = -(s^2 - 2y^2)/s^4$, difference $(2s^2 - 2x^2 - 2y^2)/s^4 = 0$. In cylindrical form the field is $\hat{\boldsymbol\phi}/s$, so on a circle about the axis $\mathbf F\cdot d\mathbf l = d\phi$: circulation $2\pi$ despite zero curl (Exercise B2).

*(iv) Coulomb $\mathbf r/r^3$.* With $\partial_xr = x/r$, $\partial_x(x/r^3) = 1/r^3 - 3x^2/r^5$. Summing over $x, y, z$: $3/r^3 - 3r^2/r^5 = 0$ for $r \ne 0$. Curl $x$-component: $\partial_y(z/r^3) - \partial_z(y/r^3) = -3yz/r^5 + 3yz/r^5 = 0$; the others likewise.

### Example 4.2: Divergence theorem on the unit cube

$\mathbf F = (x^2y + z,\ 3y^2z,\ 2xz^2)$ on $[0,1]^3$.

*Volume side.* $\nabla\cdot\mathbf F = 2xy + 6yz + 4xz$. Over the unit cube, $\int xy\,d\tau = \tfrac12\cdot\tfrac12 = \tfrac14$, and the same for $yz$, $xz$. Total: $(2 + 6 + 4)/4 = 3$.

*Surface side*, face by face (outward normals):

| Face | Normal | $\mathbf F\cdot\hat{\mathbf n}$ on face | Flux |
|---|---|---|---|
| $x = 1$ | $+\hat{\mathbf x}$ | $y + z$ | $\tfrac12 + \tfrac12 = 1$ |
| $x = 0$ | $-\hat{\mathbf x}$ | $-z$ | $-\tfrac12$ |
| $y = 1$ | $+\hat{\mathbf y}$ | $3z$ | $\tfrac32$ |
| $y = 0$ | $-\hat{\mathbf y}$ | $0$ | $0$ |
| $z = 1$ | $+\hat{\mathbf z}$ | $2x$ | $1$ |
| $z = 0$ | $-\hat{\mathbf z}$ | $0$ | $0$ |

Sum: $1 - \tfrac12 + \tfrac32 + 1 = 3$. ✓ The project reproduces every face by quadrature, then repeats the check for a non-polynomial field on a lopsided box (agreement to 12 digits).

### Example 4.3: Stokes three ways, hemisphere vs disk vs loop

$\mathbf F = (z - y,\ x,\ xy)$, $S$ the upper hemisphere $r = R = 2$, $z \ge 0$, normal outward (up). By the right-hand rule the boundary, the equator, runs counterclockwise seen from above.

*Curl.* $(\partial_y(xy) - \partial_zx,\ \partial_z(z - y) - \partial_x(xy),\ \partial_xx - \partial_y(z - y)) = (x,\ 1 - y,\ 2)$.

*Loop.* $\mathbf r(t) = R(\cos t, \sin t, 0)$, $\mathbf r' = R(-\sin t, \cos t, 0)$, $t: 0\to2\pi$. On the equator $\mathbf F = (-R\sin t,\ R\cos t,\ R^2\sin t\cos t)$, so $\mathbf F\cdot\mathbf r' = R^2\sin^2t + R^2\cos^2t = R^2$. Circulation $2\pi R^2 = 8\pi$.

*Hemisphere.* By (3.15), $d\mathbf a = R^2\sin\theta\,\hat{\mathbf r}\,d\theta\,d\phi$. With $x = R\sin\theta\cos\phi$, $y = R\sin\theta\sin\phi$:
$(\nabla\times\mathbf F)\cdot\hat{\mathbf r} = x\sin\theta\cos\phi + (1 - y)\sin\theta\sin\phi + 2\cos\theta = R\sin^2\theta\cos2\phi + \sin\theta\sin\phi + 2\cos\theta.$
Integrating over $\phi\in[0, 2\pi]$ kills $\cos2\phi$ and $\sin\phi$, leaving $4\pi\cos\theta$. Then $R^2\int_0^{\pi/2}4\pi\cos\theta\sin\theta\,d\theta = 4\pi R^2\cdot\tfrac12 = 2\pi R^2 = 8\pi$.

*Flat disk* $z = 0$, $s \le R$, normal $+\hat{\mathbf z}$ (same boundary orientation): $(\nabla\times\mathbf F)\cdot\hat{\mathbf z} = 2$, area $\pi R^2$, flux $2\pi R^2 = 8\pi$.

All three agree, as Theorem 3.10(b) requires; the project's quadratures agree to 10 digits.

### Example 4.4: Flux of $\hat{\mathbf r}/r^2$ through spheres, and Gauss's law in embryo

(a) *Concentric spheres.* §3.11 gave $4\pi$ for every radius. That is forced: the shell $R_1 < r < R_2$ excludes the origin and has $\nabla\cdot\mathbf v = 0$, so by (3.18) the flux out through $r = R_2$ minus the flux out through $r = R_1$ (where the shell's normal points inward) is $0$.

(b) *Off-centre sphere enclosing the origin* (centre $(0, 0, \tfrac12)$, radius 1). Direct integration is ugly. Apply (3.18) to the region between it and a tiny sphere $r = \epsilon$: divergence $0$ there, so the fluxes are equal, and the tiny sphere's is $4\pi$. Equivalently, $\int_V4\pi\delta^3(\mathbf r)\,d\tau = 4\pi$ by (3.32) and (3.31).

(c) *Spheres not enclosing the origin* (centre $(0,0,3)$, radius 1; centre $(2,2,0)$, radius 2, which stays $2\sqrt2 - 2$ away from the origin). $\mathbf v$ is smooth inside, so (3.18) applies directly: flux $0$.

The project's quadratures give $12.5663706 = 4\pi$ for radii $0.01$ to $1000$ and for (b), and $0$ to eight decimals for (c). Times $q/4\pi\varepsilon_0$, this is Gauss's law (EM-2).

### Example 4.5: $\nabla^2(1/r)$ directly, and a delta function you can watch form

(a) *Cartesian brute force, $r \ne 0$.* $\partial_x(1/r) = -x/r^3$, so $\partial_x^2(1/r) = -1/r^3 + 3x^2/r^5$. Summing over $x, y, z$: $-3/r^3 + 3r^2/r^5 = 0$. (The spherical route (3.29) took one line.)

(b) *A regularized point charge.* Replace the point by a unit Gaussian cloud $g_\sigma(r) = (2\pi\sigma^2)^{-3/2}e^{-r^2/2\sigma^2}$, with potential $u_\sigma(r) = \operatorname{erf}\big(r/(\sqrt2\,\sigma)\big)/r\to1/r$ as $\sigma\to0$. Check with (3.29): write $E(r) = \operatorname{erf}\big(r/(\sqrt2\,\sigma)\big)$, so $u = E/r$ and $r^2u' = rE' - E$. Then $\partial_r(r^2u') = rE''$, so $\nabla^2u_\sigma = E''/r$. Since $E' = \sqrt{2/\pi}\,\sigma^{-1}e^{-r^2/2\sigma^2}$, $E'' = -(r/\sigma^2)E'$, and
$$\nabla^2u_\sigma = -\sqrt{2/\pi}\,\sigma^{-3}e^{-r^2/2\sigma^2} = -4\pi g_\sigma(r),$$
because $4\pi(2\pi)^{-3/2} = \sqrt{2/\pi}$. Now test against $\psi = e^{-r^2}$ (so $\psi(0) = 1$):
$$\int\psi\,\nabla^2u_\sigma\,d\tau = -4\pi(2\pi\sigma^2)^{-3/2}\int e^{-r^2(1 + 1/2\sigma^2)}\,d\tau = -4\pi(2\pi\sigma^2)^{-3/2}\left(\frac{2\pi\sigma^2}{1 + 2\sigma^2}\right)^{3/2} = -4\pi(1 + 2\sigma^2)^{-3/2},$$
using $\int e^{-ar^2}d\tau = (\pi/a)^{3/2}$. For $\sigma = 0.5, 0.2, 0.1, 0.05$: $-6.8403, -11.1963, -12.1986, -12.4727$, approaching $-4\pi\psi(0) = -12.5664$. The project confirms each by quadrature, checks $\int(1/r)\nabla^2\psi\,d\tau = -4\pi$ (the distributional meaning of (3.33)), and sums a finite-difference Laplacian of $1/r$ over a grid: $-12.568$, 99.9% of it within four grid spacings of the origin.

### Example 4.6: A potential by line integration (a CM-1 force)

A particle feels $\mathbf F = (2xy - z,\ x^2 + z^2,\ 2yz - x)$ (N, positions in m). Is $\mathbf F$ conservative? Find $U$ with $\mathbf F = -\nabla U$ and the work from the origin to $(1, 2, 3)$.

*Curl:* $(2z - 2z,\ -1 + 1,\ 2x - 2x) = \mathbf 0$ (components in the order of (3.6)). $\mathbb R^3$ is star-shaped, so Theorem 3.16 guarantees $U$.

*Construct $U$* along the axis path $(0,0,0)\to(x,0,0)\to(x,y,0)\to(x,y,z)$, one coordinate per leg:
- Leg 1 ($y = z = 0$): $F_x = 0$. Contributes $0$.
- Leg 2 ($x$ fixed, $z = 0$): $F_y = x^2$. Contributes $x^2y$.
- Leg 3 ($x, y$ fixed): $F_z = 2yz' - x$. Contributes $yz^2 - xz$.

Work from the origin: $x^2y + yz^2 - xz$, so $U = -x^2y - yz^2 + xz$ (with $U(\mathbf 0) = 0$). Check: $-\nabla U = (2xy - z,\ x^2 + z^2,\ 2yz - x)$. ✓

*Work:* $U(1,2,3) = -2 - 18 + 3 = -17$, so $W = U(\mathbf 0) - U(1,2,3) = 17$ J. The project gets $17$ along the axis path, the straight line, and the curve $(t^2, 2\sin(\pi t/2), 3t^3)$.

### Example 4.7: A Helmholtz split by hand

$\mathbf F = (\sin(x + y), 0, 0)$ on the periodic box $[0, 2\pi)^3$. Since $\sin(x + y) = (e^{i(x+y)} - e^{-i(x+y)})/2i$, only $\mathbf k = \pm(1, 1, 0)$ appear, with $\hat{\mathbf F}_{\mathbf k}\parallel\hat{\mathbf x}$ and $\hat{\mathbf k}(\hat{\mathbf k}\cdot\hat{\mathbf x}) = (1, 1, 0)/2$ for either sign, so by (3.35)
$$\mathbf F = \underbrace{\tfrac12\sin(x+y)(1, 1, 0)}_{\text{curl-free}} + \underbrace{\tfrac12\sin(x+y)(1, -1, 0)}_{\text{divergence-free}}.$$
Check: the first is $-\nabla U$ with $U = \tfrac12\cos(x + y)$, since $\nabla U = -\tfrac12\sin(x+y)(1,1,0)$. The second is $\nabla\times\mathbf W$ with $\mathbf W = (0, 0, w)$, $w = -\tfrac12\cos(x+y)$: $\nabla\times\mathbf W = (\partial_yw, -\partial_xw, 0) = \tfrac12\sin(x+y)(1, -1, 0)$. ✓ From (3.35) directly, for $\mathbf k = (1,1,0)$, $\hat{\mathbf F}_{\mathbf k} = \hat{\mathbf x}/2i$ gives $\hat U_{\mathbf k} = i(1/2i)/2 = \tfrac14$, and the same for $-\mathbf k$, so $U = \tfrac14e^{i(x+y)} + \tfrac14e^{-i(x+y)} = \tfrac12\cos(x+y)$. ✓ The parts are perpendicular pointwise. The project reproduces this with the FFT, then splits a random smooth field and checks every property to better than $10^{-13}$.

## 5. Why it matters and where it connects

- **EM-2:** (3.32) plus superposition is $\nabla\cdot\mathbf E = \rho/\varepsilon_0$, and Example 4.4 is Gauss's law; Theorem 3.16 turns $\nabla\times\mathbf E = \mathbf 0$ into a potential; (3.33) is Poisson's equation with its Green's function; (3.22) turns the energy into $\tfrac{\varepsilon_0}{2}\int E^2\,d\tau$.
- **EM-3:** uniqueness from (3.22), separation of variables in (3.29), Green's functions from (3.23).
- **EM-5:** $\nabla\cdot\mathbf B = 0$ and $\nabla\times\mathbf B = \mu_0\mathbf J$ from Biot–Savart via (P4), (3.11), (3.33); Theorem 3.17 gives $\mathbf B = \nabla\times\mathbf A$, and (3.34)'s $\mathbf W$ with $\mathbf C = \mu_0\mathbf J$ is the vector potential. **EM-6–7:** Theorem 3.15 makes Maxwell's equations complete; (3.13) gives the wave equation.
- **CM-1 Theorem 3.3** is completed by Theorem 3.16; **CM-3** uses §3.5's $\mathbf v\times(\nabla\times\mathbf A)$ identity.
- **F3:** $\delta^3$ is a distribution (F3 §3.10); (3.35) is F3's derivative rule. **QM-5, QM-6:** the angular part of (3.29) is $-\hat L^2/\hbar^2r^2$. **ML-10:** Hyvärinen's score-matching identity $\int p\,\nabla\cdot\mathbf s\,d\tau = -\int\mathbf s\cdot\nabla p\,d\tau$ is (P3) plus (3.18) with a vanishing surface term.
- **Fluids** (your ME courses): continuity $\partial_t\rho + \nabla\cdot(\rho\mathbf u) = 0$ is (3.7) read as conservation; vorticity is $\nabla\times\mathbf u$.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 drill the operators, triple products, a gradient, curvilinear divergences, delta integrals, and a line integral. B1–B8 are the module; B3 and B4 (set up a surface integral and get the orientation right), B6 (where the divergence theorem breaks and how $\delta^3$ repairs it), B7 (build a vector potential and change gauge), and B8 (code) are essential. C1–C4: solid angle, uniqueness in Helmholtz via Green's identity, the screened Coulomb potential, and why the curl is a pseudovector.

## 7. Coding project

`project.py` checks every claim and number above, in ten parts: SymPy operators, identities and product rules (Example 4.1, Theorems 3.5–3.6); curvilinear Jacobians, scale factors, and the whole §3.10 table against Cartesian; finite-difference div and curl on a grid vs SymPy (second order, observed 1.98); a divergence-theorem checker (Example 4.2, a random box, B3); a Stokes checker (Example 4.3, a tilted parallelogram, the Theorem 3.4 limit, B4); fluxes of $\hat{\mathbf r}/r^2$ (Example 4.4, C1); the delta function three ways (Example 4.5); Helmholtz by FFT via (3.35) on Example 4.7 and a random field; potentials (Example 4.6, B7); and the numbers in the solutions, with reference solutions for B8, C1(c), C4 after "Your turn."

**Build it yourself:** write `fd_div_curl` (Part 3) with `np.gradient`, and `flux_box` (Part 4: six faces, outward signs) and reproduce Example 4.2's table, before reading them. Wrong face signs show up immediately.

## 8. Where you'll struggle

- **The delta function feels like cheating** the first time, and it is the single most-used fact in E&M: every Gauss's-law calculation, every Green's function. Stop thinking of $\delta^3$ as a function. Rerun the §3.11 derivation (excise a ball, apply the divergence theorem, shrink) until you can do it closed-book, and watch Example 4.5's numbers converge in the project.
- **The sign of the Laplacian.** $\nabla\cdot(\hat{\mathbf r}/r^2) = +4\pi\delta^3$ but $\nabla^2(1/r) = -4\pi\delta^3$. Rederive the sign from $\nabla(1/r) = -\hat{\mathbf r}/r^2$ ($1/r$ decreases outward, so its gradient points in) instead of memorizing it. Poisson's $\nabla^2V = -\rho/\varepsilon_0$ inherits the minus.
- **Stokes orientation is wrong in half of first attempts.** Draw the loop's arrow and the normal's arrow and check the right-hand rule before integrating. Closed surface: normal outward. Open surface: normal *chosen* to match the loop. Drill with Example 4.3 and B4.
- **Freezing on a surface integral.** Knowing the pieces is not setting one up. Use the §3.7 recipe mechanically, writing each step: parametrize, cross the tangents, check direction, substitute, integrate. B3 and C1 exist for this.
- **The spherical curl is unlearnable by memory.** Derive (3.27) and (3.28); look up the rest. Never take a vector Laplacian component-by-component in curvilinear coordinates.
- **Index notation is slower for one line and much faster for three.** If (3.13) takes you more than five minutes, redo LA-12 Exercise B6, then B5 here.
- **"Curl-free implies gradient" is false on regions with holes.** Each time you use Theorem 3.16, name the region.

## 9. Check yourself

1. State the $\varepsilon$–$\delta$ identity and use it to derive BAC-CAB in one line.
2. What is the divergence of a field, as a limit? What is the component of the curl along $\hat{\mathbf n}$, as a limit?
3. Why are $\nabla\times\nabla f$ and $\nabla\cdot(\nabla\times\mathbf A)$ zero? One sentence, no components.
4. Derive $\nabla\times(\nabla\times\mathbf A) = \nabla(\nabla\cdot\mathbf A) - \nabla^2\mathbf A$ in three lines of index notation.
5. Write $d\mathbf a$ for a sphere of radius $R$ and for a disk in the $xy$-plane. Which way do they point?
6. State the divergence theorem and Stokes' theorem, including the orientation rule. Why do interior faces/edges cancel?
7. Derive Green's first identity from the divergence theorem.
8. What are the scale factors in spherical coordinates? Derive $\nabla f$ in spherical coordinates from $df = \nabla f\cdot d\mathbf l$.
9. Compute $\nabla\cdot(\hat{\mathbf r}/r^2)$ for $r \ne 0$ and the flux through a sphere. Resolve the contradiction.
10. What is $\nabla^2(1/\mathscr r)$, and in what sense is that equation true?
11. State the Helmholtz theorem with its hypotheses. Why does it make Maxwell's equations complete?
12. Give a curl-free field that is not a gradient, and the hypothesis it violates.

## 10. Pair with

- **Griffiths, *Introduction to Electrodynamics* (4th ed.), ch. 1.** The standard treatment, in the same notation (script $\mathscr r$, $d\tau$, $s$ for the cylindrical radius). Its problems are better drilled than ours, and §1.5 (the delta function) and §1.6 (Helmholtz, potentials) are short and clear. The inside front cover is where you look up the spherical curl.
- **Schey, *Div, Grad, Curl, and All That*.** Builds divergence and curl from flux and circulation with fluid pictures, slowly and physically. Read it if §3.4 felt fast.
- **Purcell & Morin, *Electricity and Magnetism* (3rd ed.), ch. 2 and the appendix on curvilinear coordinates.** Derives div and curl from Gauss's law and circulation in the physical setting, and the appendix derives the curvilinear formulas, including the curl we skipped.
- **3Blue1Brown, "Divergence and curl: The language of Maxwell's equations, fluid flow, and more."** Twenty minutes of animated paddlewheels and sources. Watch it after §2.
