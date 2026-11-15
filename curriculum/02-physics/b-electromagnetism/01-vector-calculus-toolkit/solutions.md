# EM-1 Solutions

Every number below is asserted or printed by `project.py` (Part 10 and the reference-solution cell).

## A. Warm-up

**A1.** (a) $\nabla f = (\partial_x, \partial_y, \partial_z)(x^2y + yz^3) = (2xy,\ x^2 + z^3,\ 3yz^2)$, which at $(1,1,1)$ is $(2, 2, 3)$.
$\nabla^2f = \partial_x(2xy) + \partial_y(x^2 + z^3) + \partial_z(3yz^2) = 2y + 0 + 6yz$, which is $8$ at $(1,1,1)$.
$\nabla\cdot\mathbf F = \partial_x(xy) + \partial_y(yz) + \partial_z(zx) = y + z + x$, which is $3$.
$\nabla\times\mathbf F = (\partial_y(zx) - \partial_z(yz),\ \partial_z(xy) - \partial_x(zx),\ \partial_x(yz) - \partial_y(xy)) = (-y,\ -z,\ -x)$, which is $(-1, -1, -1)$.

(b) $\nabla\times\nabla f = (\partial_y(3yz^2) - \partial_z(x^2 + z^3),\ \partial_z(2xy) - \partial_x(3yz^2),\ \partial_x(x^2 + z^3) - \partial_y(2xy)) = (3z^2 - 3z^2,\ 0 - 0,\ 2x - 2x) = \mathbf 0$.
$\nabla\cdot(\nabla\times\mathbf F) = \partial_x(-y) + \partial_y(-z) + \partial_z(-x) = 0$.

**A2.** (a) $\mathbf B\times\mathbf C = (1\cdot1 - 3\cdot0,\ 3\cdot2 - 0\cdot1,\ 0\cdot0 - 1\cdot2) = (1, 6, -2)$, so $\mathbf A\cdot(\mathbf B\times\mathbf C) = 1 + 12 + 0 = 13$. This is the volume of the parallelepiped with edges $\mathbf A, \mathbf B, \mathbf C$; it is positive because the three form a right-handed triple. Cyclic check: $\mathbf C\times\mathbf A = (0\cdot0 - 1\cdot2,\ 1\cdot1 - 2\cdot0,\ 2\cdot2 - 0\cdot1) = (-2, 1, 4)$ and $\mathbf B\cdot(\mathbf C\times\mathbf A) = 0 + 1 + 12 = 13$.

(b) Directly: $\mathbf A\times(1, 6, -2) = (2\cdot(-2) - 0\cdot6,\ 0\cdot1 - 1\cdot(-2),\ 1\cdot6 - 2\cdot1) = (-4, 2, 4)$. By (3.2): $\mathbf A\cdot\mathbf C = 2$ and $\mathbf A\cdot\mathbf B = 2$, so $\mathbf B(2) - \mathbf C(2) = (0, 2, 6) - (4, 0, 2) = (-4, 2, 4)$. ✓

(c) $\mathbf A\times\mathbf B = (2\cdot3 - 0\cdot1,\ 0\cdot0 - 1\cdot3,\ 1\cdot1 - 2\cdot0) = (6, -3, 1)$, and $(6, -3, 1)\times(2, 0, 1) = (-3\cdot1 - 1\cdot0,\ 1\cdot2 - 6\cdot1,\ 6\cdot0 - (-3)\cdot2) = (-3, -4, 6)$. This is not $(-4, 2, 4)$: the cross product is not associative. Via (3.2): $(\mathbf A\times\mathbf B)\times\mathbf C = -\mathbf C\times(\mathbf A\times\mathbf B) = \mathbf B(\mathbf A\cdot\mathbf C) - \mathbf A(\mathbf B\cdot\mathbf C) = (0,2,6) - 3(1,2,0) = (-3, -4, 6)$. ✓

**A3.** $\nabla T = (2x + y,\ 4y + x,\ 6z) = (1, -3, 12)$ at $(1, -1, 2)$. Fastest rise is along $(1, -3, 12)/\sqrt{154}$, at the rate $|\nabla T| = \sqrt{1 + 9 + 144} = \sqrt{154} = 12.4097$ °C/m (Theorem 3.2). Along $\hat{\mathbf u} = (1,1,1)/\sqrt3$: $\nabla T\cdot\hat{\mathbf u} = (1 - 3 + 12)/\sqrt3 = 10/\sqrt3 = 5.7735$ °C/m.

**A4.** (a) With $F_r = r^n$ and the other components zero, (3.28) gives $\frac{1}{r^2}\partial_r(r^{n+2}) = (n + 2)\,r^{n-1}$. Checks: $n = 1$ gives $3$ (Example 4.1); $n = -2$ gives $0$ (§3.11).

(b) $\frac1s\partial_s(s\cdot s) = \frac1s(2s) = 2$. Cartesian: $\partial_xx + \partial_yy + \partial_z0 = 2$. ✓

(c) $\partial_r(r^n) = nr^{n-1}$, so $\frac{1}{r^2}\partial_r(r^2\cdot nr^{n-1}) = \frac{n(n+1)r^n}{r^2} = n(n+1)\,r^{n-2}$. This vanishes exactly for $n = 0$ and $n = -1$: the only harmonic powers are constants and $1/r$.

**A5.** (a) By the sifting property (3.31), the integral is $|\mathbf a|^2 = 1 + 4 + 1 = 6$.

(b) The point $2\hat{\mathbf x}$ is at distance $2$ from the origin, outside the unit ball, so by (3.31) the integral is $0$.

(c) By (3.33) with $\mathbf r' = \hat{\mathbf z}$, $\nabla^2(1/|\mathbf r - \hat{\mathbf z}|) = -4\pi\delta^3(\mathbf r - \hat{\mathbf z})$. Sifting: $-4\pi(|\hat{\mathbf z}|^2 + 1) = -8\pi = -25.1327$. The project also replaces $\delta^3$ by a Gaussian of width $\sigma$ centred at $\hat{\mathbf z}$; then the integral is $-4\pi(2 + 3\sigma^2)$ (the mean of $|\mathbf r|^2$ under that Gaussian is $1 + 3\sigma^2$), giving $-26.6407$ at $\sigma = 0.2$ and $-25.2270$ at $\sigma = 0.05$, tending to $-8\pi$.

**A6.** (a) $\mathbf r' = (-\sin t, \cos t, 1/2\pi)$ and $\mathbf F(\mathbf r(t)) = (\sin t, -\cos t, t/2\pi)$, so $\mathbf F\cdot\mathbf r' = -\sin^2t - \cos^2t + t/4\pi^2 = -1 + t/4\pi^2$. Integrating over $[0, 2\pi]$: $-2\pi + \frac{(2\pi)^2}{2\cdot4\pi^2} = \tfrac12 - 2\pi = -5.7832$.

(b) $\mathbf r(t) = (1, 0, t)$, $\mathbf r' = (0, 0, 1)$, $\mathbf F = (0, -1, t)$, so $\mathbf F\cdot\mathbf r' = t$ and the integral is $\tfrac12$.

(c) $\nabla\times\mathbf F = (\partial_yz - \partial_z(-x),\ \partial_zy - \partial_xz,\ \partial_x(-x) - \partial_yy) = (0, 0, -2)$. Going up the helix and back down the segment is a closed loop $C$ with $\oint_C\mathbf F\cdot d\mathbf l = (\tfrac12 - 2\pi) - \tfrac12 = -2\pi$. By Stokes this equals $\int_S(\nabla\times\mathbf F)\cdot d\mathbf a = -2\int_S\hat{\mathbf z}\cdot d\mathbf a$ for any $S$ bounded by $C$. To evaluate $\int_S\hat{\mathbf z}\cdot d\mathbf a$, note $\hat{\mathbf z} = \nabla\times\mathbf A$ with $\mathbf A = \tfrac12(-y, x, 0)$, so by Stokes again it is $\oint_C\tfrac12(x\,dy - y\,dx)$. On the helix $x\,dy - y\,dx = (\cos^2t + \sin^2t)\,dt$, contributing $\tfrac12\cdot2\pi = \pi$; on the vertical segment $x = 1$, $y = 0$ are constant, contributing $0$. So the flux is $-2\pi$, matching. The two paths disagree because $\mathbf F$ has curl, and the discrepancy is exactly the flux of the curl through the surface between them. The project prints helix $-5.7832$, segment $0.5000$, difference $-6.2832$.

## B. Core

**B1.** *From (P5).* With $\mathbf A = \nabla g$: $\nabla\times(f\nabla g) = f\,\nabla\times\nabla g - \nabla g\times\nabla f = \mathbf 0 + \nabla f\times\nabla g$, using (3.11) and $\mathbf a\times\mathbf b = -\mathbf b\times\mathbf a$.

*Index notation.* $[\nabla\times(f\nabla g)]_i = \varepsilon_{ijk}\partial_j(f\,\partial_kg) = \varepsilon_{ijk}(\partial_jf)(\partial_kg) + f\,\varepsilon_{ijk}\partial_j\partial_kg$. The first term is $(\nabla f\times\nabla g)_i$. The second is zero by the lemma in Theorem 3.5 ($\partial_j\partial_kg$ is symmetric).

Stokes (3.19) with $\mathbf F = f\nabla g$ gives $\oint_Pf\nabla g\cdot d\mathbf l = \int_S(\nabla f\times\nabla g)\cdot d\mathbf a$. With $g = f$, $\nabla f\times\nabla f = \mathbf 0$, so $\oint f\nabla f\cdot d\mathbf l = 0$ (consistent with $f\nabla f = \nabla(f^2/2)$ and (3.20)). The project checks the identity symbolically for arbitrary $f, g$.

**B2.** (a) For $\mathbf F = F_\phi(s)\hat{\boldsymbol\phi}$ the cylindrical curl's $\hat{\mathbf s}$ and $\hat{\boldsymbol\phi}$ components vanish ($F_z = F_s = 0$, $\partial_zF_\phi = 0$), leaving $\frac1s\partial_s(sF_\phi)\,\hat{\mathbf z}$.
For $s\hat{\boldsymbol\phi}$: $\frac1s\partial_s(s^2) = 2$, so the curl is $2\hat{\mathbf z}$, as in Example 4.1.
For $\hat{\boldsymbol\phi}/s$: $\frac1s\partial_s(1) = 0$ for $s \ne 0$.

(b) On the circle of radius $a$ about the axis, $d\mathbf l = a\,d\phi\,\hat{\boldsymbol\phi}$, so $\mathbf F\cdot d\mathbf l = \frac1a\cdot a\,d\phi = d\phi$ and the circulation is $2\pi$ for every $a$. The unit circle centred at $(2, 0, 0)$ and the disk it bounds lie in $s \ge 1$, where $\mathbf F$ is smooth with zero curl, so Stokes gives circulation $0$. (Equivalently: on the half-space $x > 0$, $\mathbf F = \nabla\phi$ with the single-valued $\phi = \arctan(y/x)$, and (3.20) applies.) The project prints $6.283185$ for $a = 0.5$ and $a = 3$, and $-1.2\times10^{-14}$ for the off-axis circle.

(c) Stokes needs $\mathbf F$ to be $C^1$ on a surface spanning the loop. Every surface spanning a loop around the axis is pierced by the axis, where $\mathbf F$ is undefined and unbounded, so Stokes does not apply and there is no contradiction. To find the curl as a distribution, repeat §3.11 in two dimensions. Let $\psi(x, y)$ be a smooth test function vanishing for $s \ge R$. For smooth fields, integration by parts gives $\int\psi\,(\partial_xF_y - \partial_yF_x)\,dx\,dy = -\int(F_y\partial_x\psi - F_x\partial_y\psi)\,dx\,dy$, and the right side is the definition of the weak curl.
1. On the annulus $\epsilon < s < R$, $\partial_xF_y - \partial_yF_x = 0$, so $\partial_x(\psi F_y) - \partial_y(\psi F_x) = F_y\partial_x\psi - F_x\partial_y\psi$.
2. Green's theorem in the plane (Theorem 3.4, step 3, summed) on the annulus: the integral of the left side equals the counterclockwise circulation of $\psi\mathbf F$ around $s = R$ (zero, since $\psi = 0$ there) minus the counterclockwise circulation around $s = \epsilon$, which is $\oint\psi\,d\phi$.
3. As $\epsilon\to0$, $\oint\psi\,d\phi\to2\pi\psi(0)$. So $-\int(F_y\partial_x\psi - F_x\partial_y\psi)\,dx\,dy = 2\pi\psi(0)$: the weak $z$-curl is $2\pi\delta(x)\delta(y)$, and since nothing depends on $z$, $\nabla\times\mathbf F = 2\pi\delta(x)\delta(y)\hat{\mathbf z}$.

(d) Simple connectivity. $\mathbb R^3$ minus the $z$-axis has loops (around the axis) that cannot be shrunk to a point inside the region. Locally $\mathbf F = \nabla\phi$, but $\phi$ increases by $2\pi$ each time around, so no single-valued potential exists.

**B3.** *Volume.* $\nabla\cdot\mathbf F = z + z + 2z = 4z$. In spherical coordinates,
$\int_V4z\,d\tau = 4\int_0^R\!\int_0^{\pi/2}\!\int_0^{2\pi}r\cos\theta\,r^2\sin\theta\,d\phi\,d\theta\,dr = 4\cdot\frac{R^4}{4}\cdot\frac12\cdot2\pi = \pi R^4.$

*Hemisphere.* By (3.15), $d\mathbf a = R^2\sin\theta\,\hat{\mathbf r}\,d\theta\,d\phi$ (outward). On $r = R$,
$\mathbf F\cdot\hat{\mathbf r} = \frac{x\cdot xz + y\cdot yz + z(z^2 + 1)}{R} = \frac{z(x^2 + y^2 + z^2 + 1)}{R} = (R^2 + 1)\cos\theta,$
using $z = R\cos\theta$. Flux: $(R^2 + 1)R^2\int_0^{2\pi}\!\int_0^{\pi/2}\cos\theta\sin\theta\,d\theta\,d\phi = (R^2 + 1)R^2\cdot2\pi\cdot\frac12 = \pi R^4 + \pi R^2$.

*Base.* On $z = 0$ the outward normal is $-\hat{\mathbf z}$ and $F_z = 0 + 1 = 1$, so the flux is $-\pi R^2$.

*Total:* $\pi R^4 + \pi R^2 - \pi R^2 = \pi R^4$, equal to the volume integral. ✓ For $R = 2$: hemisphere $20\pi = 62.8319$, base $-4\pi = -12.5664$, total and volume integral $16\pi = 50.2655$.

**B4.** (a) Parametrize each edge as $\mathbf r = \mathbf a + (\mathbf b - \mathbf a)t$, $t\in[0, 1]$.
- $(1,0,0)\to(0,1,0)$: $\mathbf r = (1 - t, t, 0)$, $\mathbf r' = (-1, 1, 0)$, $\mathbf F = (t, 0, 1 - t)$, $\mathbf F\cdot\mathbf r' = -t$, integral $-\tfrac12$.
- $(0,1,0)\to(0,0,1)$: $\mathbf r = (0, 1 - t, t)$, $\mathbf r' = (0, -1, 1)$, $\mathbf F = (1 - t, t, 0)$, $\mathbf F\cdot\mathbf r' = -t$, integral $-\tfrac12$.
- $(0,0,1)\to(1,0,0)$: $\mathbf r = (t, 0, 1 - t)$, $\mathbf r' = (1, 0, -1)$, $\mathbf F = (0, 1 - t, t)$, $\mathbf F\cdot\mathbf r' = -t$, integral $-\tfrac12$.

Total $-\tfrac32$.

(b) $\nabla\times\mathbf F = (\partial_yx - \partial_zz,\ \partial_zy - \partial_xx,\ \partial_xz - \partial_yy) = (-1, -1, -1)$. The triangle is flat, so its vector area is $\tfrac12(\mathbf v_2 - \mathbf v_1)\times(\mathbf v_3 - \mathbf v_1)$ with the vertices in the order traversed (this ordering *is* the right-hand rule): $\tfrac12(-1, 1, 0)\times(-1, 0, 1) = \tfrac12(1, 1, 1)$. The normal points away from the origin, toward $(1,1,1)$. Since the curl is constant, the flux is $(-1,-1,-1)\cdot\tfrac12(1,1,1) = -\tfrac32$. ✓

(c) Reversing the traversal flips the circulation to $+\tfrac32$ and, by the right-hand rule, the normal to $-(1,1,1)/\sqrt3$, so the flux is also $+\tfrac32$. Stokes holds as long as the loop and the normal are paired by the right-hand rule; reversing only one of them produces the classic sign error.

**B5.** *(P5).* $[\nabla\times(f\mathbf A)]_i = \varepsilon_{ijk}\partial_j(fA_k) = f\varepsilon_{ijk}\partial_jA_k + \varepsilon_{ijk}(\partial_jf)A_k = f(\nabla\times\mathbf A)_i + (\nabla f\times\mathbf A)_i = f(\nabla\times\mathbf A)_i - (\mathbf A\times\nabla f)_i$.

*(P6).*
1. $[\nabla\times(\mathbf A\times\mathbf B)]_i = \varepsilon_{ijk}\partial_j(\varepsilon_{klm}A_lB_m) = \varepsilon_{kij}\varepsilon_{klm}\partial_j(A_lB_m)$.
2. By (3.1): $(\delta_{il}\delta_{jm} - \delta_{im}\delta_{jl})\partial_j(A_lB_m) = \partial_j(A_iB_j) - \partial_j(A_jB_i)$.
3. Product rule: $B_j\partial_jA_i + A_i\partial_jB_j - B_i\partial_jA_j - A_j\partial_jB_i$.
4. These are, in order, $[(\mathbf B\cdot\nabla)\mathbf A]_i$, $A_i(\nabla\cdot\mathbf B)$, $-B_i(\nabla\cdot\mathbf A)$, $-[(\mathbf A\cdot\nabla)\mathbf B]_i$, which is (P6). ∎

*Check.* $\mathbf A = \boldsymbol\omega$ constant, $\mathbf B = \mathbf r$: $(\mathbf r\cdot\nabla)\boldsymbol\omega = \mathbf 0$; $(\boldsymbol\omega\cdot\nabla)\mathbf r = \boldsymbol\omega$ (since $\omega_j\partial_jx_i = \omega_i$); $\boldsymbol\omega(\nabla\cdot\mathbf r) = 3\boldsymbol\omega$; $\mathbf r(\nabla\cdot\boldsymbol\omega) = \mathbf 0$. Total $0 - \boldsymbol\omega + 3\boldsymbol\omega - 0 = 2\boldsymbol\omega$, as in Example 4.1. The project verifies (P5) and (P6) for arbitrary symbolic fields.

**B6.** (a) $n = -1$: by A4(a), $\nabla\cdot\mathbf F_{-1} = 1/r^2$. Volume: $\int_0^R\frac{1}{r^2}4\pi r^2\,dr = 4\pi R$. Flux: $\frac1R\cdot4\pi R^2 = 4\pi R$. ✓ (The field is singular at the origin, but $1/r^2$ times $r^2\,dr$ is integrable.)

(b) Volume: $\int_0^R(n+2)r^{n-1}4\pi r^2\,dr = 4\pi(n + 2)\frac{R^{n+2}}{n+2} = 4\pi R^{n+2}$; the integral converges at $0$ because $n + 1 > -1$. Flux: $R^n\cdot4\pi R^2 = 4\pi R^{n+2}$. ✓ Why the singularity doesn't matter: on the shell $\epsilon < r < R$ the theorem applies, and the inner sphere's flux $4\pi\epsilon^{n+2}\to0$.

(c) $n = -2$: $\nabla\cdot\mathbf F_{-2} = 0$ for $r \ne 0$, so the naive volume integral is $0$, but the flux is $4\pi$. Theorem 3.14 repairs it: $\nabla\cdot(\hat{\mathbf r}/r^2) = 4\pi\delta^3(\mathbf r)$, whose integral over the ball is $4\pi$. In the shell argument, the inner flux is $4\pi\epsilon^0 = 4\pi$ for every $\epsilon$, and that leftover is the delta function.

(d) $n = -3$: flux through $r = R$ is $R^{-3}\cdot4\pi R^2 = 4\pi/R$; through $r = \epsilon$ (outward from the origin) it is $4\pi/\epsilon$. On the shell, $\nabla\cdot\mathbf F_{-3} = -r^{-4}$ and $\int_\epsilon^R-r^{-4}\,4\pi r^2\,dr = 4\pi/R - 4\pi/\epsilon$, consistent with (3.18) on the shell. As $\epsilon\to0$ the inner flux $4\pi/\epsilon$ diverges, so the "point contribution" would have to be infinite: no finite multiple of $\delta^3$ works. Underneath: $|\mathbf F_{-3}| = r^{-3}$ is not integrable near the origin ($\int r^{-3}r^2\,dr$ diverges), so $\mathbf F_{-3}$ does not even define a distribution by integration. Only $n = -2$ sits exactly on the boundary where the inner flux is finite and nonzero.

**B7.** (a) $\partial_x(y^2) + \partial_y(z^2) + \partial_z(x^2) = 0$.

(b) From (3.36): $W_x = \int_0^zF_y\,dz' = \int_0^zz'^2\,dz' = z^3/3$; $W_y = \int_0^xF_z(x', y, 0)\,dx' - \int_0^zF_x\,dz' = \int_0^xx'^2\,dx' - \int_0^zy^2\,dz' = x^3/3 - y^2z$; $W_z = 0$. So $\mathbf W = (z^3/3,\ x^3/3 - y^2z,\ 0)$. Curl: $(\partial_yW_z - \partial_zW_y,\ \partial_zW_x - \partial_xW_z,\ \partial_xW_y - \partial_yW_x) = (y^2,\ z^2,\ x^2)$. ✓

(c) $\nabla\cdot\mathbf W = \partial_y(x^3/3 - y^2z) = -2yz$. We need $\nabla\cdot(\mathbf W + \nabla\lambda) = 0$, i.e. $\nabla^2\lambda = 2yz$. Try $\lambda = y^3z/3$: $\nabla^2\lambda = \partial_y^2(y^3z/3) = 2yz$. ✓ Then $\nabla\lambda = (0,\ y^2z,\ y^3/3)$ and $\mathbf W' = (z^3/3,\ x^3/3,\ y^3/3)$. Check: $\nabla\cdot\mathbf W' = 0$ (each component is independent of its own variable), and $\nabla\times\mathbf W' = (\partial_y\tfrac{y^3}{3} - \partial_z\tfrac{x^3}{3},\ \partial_z\tfrac{z^3}{3} - \partial_x\tfrac{y^3}{3},\ \partial_x\tfrac{x^3}{3} - \partial_y\tfrac{z^3}{3}) = (y^2, z^2, x^2)$. ✓

(d) No. If $\nabla^2\chi = 0$, then $\mathbf W' + \nabla\chi$ has the same curl (3.11) and divergence $\nabla^2\chi = 0$. For example $\chi = xy$ gives $\mathbf W' + (y, x, 0)$. The Coulomb-gauge potentials of $\mathbf F$ are exactly $\mathbf W' + \nabla\chi$ with $\chi$ harmonic; pinning down one requires a condition at infinity or on a boundary (compare C2).

**B8.** (a) The reference cell computes

```python
G, dv_fA, _ = fd_div_curl(numfield(f_b8 * A_b8), n)          # D·(fA)
_, dv_A, _ = fd_div_curl(numfield(A_b8), n)                   # D·A
rule = fg * dv_A + sum(Ac * np.gradient(fg, h, axis=ax, edge_order=2)
                       for ax, Ac in enumerate(numfield(A_b8)(*G)))   # f D·A + A·Df
gap = np.abs(dv_fA - rule).max()
```

The gaps are $5.61\times10^{-3}$ at $n = 21$ ($h = 0.1$) and $1.41\times10^{-3}$ at $n = 41$ ($h = 0.05$), ratio $3.99$. The discrete product rule fails at order $h^2$. For central differences $D$, algebra gives
$$D(fg)_i - (f_iDg_i + g_iDf_i) = \frac{(f_{i+1} - f_i)(g_{i+1} - g_i) - (f_{i-1} - f_i)(g_{i-1} - g_i)}{2h}.$$
Each product is $h^2f'g'$ plus $O(h^3)$, and the leading terms cancel, so the gap is $O(h^3)/h = O(h^2)$; `np.gradient`'s second-order one-sided edge formulas behave the same way. Halving $h$ divides the gap by $4$.
With $\mathbf A = (y, z^2, \sin x)$, each component is independent of the variable it is differentiated by in $\nabla\cdot$ (for instance $A_x = y$ is constant along $x$), so $D_x(fA_x) = A_xD_xf$ holds exactly on the grid. The gap would be roundoff ($\sim10^{-15}$) and would test nothing.

(b) The path $(\sin(\pi t/2),\ 2t,\ 3t^2)$, $t\in[0,1]$, runs from the origin to $(1,2,3)$ with $\mathbf r' = (\tfrac\pi2\cos(\pi t/2),\ 2,\ 6t)$; the project gets $17.0000000000$. For the vortex, along the upper half-circle the azimuth $\phi$ goes $0\to\pi$ and along the lower $0\to-\pi$; since $\mathbf F\cdot d\mathbf l = d\phi$ (B2), the integrals are $+\pi$ and $-\pi$ (project: $3.141593$, $-3.141593$). Same endpoints, different answers, although $\nabla\times\mathbf F = \mathbf 0$ along both paths. Their difference, $2\pi$, is the circulation around the closed loop "upper path, then lower path reversed", which encircles the axis: B2's $2\pi$. Theorem 3.16's hypothesis of simple connectivity fails.

## C. Challenge

**C1.** (a) At a point of $S$ at distance $r$, $\hat{\mathbf r}\cdot d\mathbf a$ is the area of the patch projected onto the plane perpendicular to the ray, i.e. onto the sphere of radius $r$ through that point. Scaling a patch of the radius-$r$ sphere down to the unit sphere divides its area by $r^2$. So $\hat{\mathbf r}\cdot d\mathbf a/r^2 = d\Omega$, the area of the patch's radial shadow on the unit sphere, and $\Omega$ is the area of the whole shadow. A closed surface around the origin casts a shadow covering the unit sphere once: $\Omega = 4\pi$ (the same answer as (3.32) with (3.18)).

(b) On the plane $z = z_0$, with $d\mathbf a = \hat{\mathbf z}\,s\,ds\,d\phi$ and $r^2 = s^2 + z_0^2$: $\hat{\mathbf r}\cdot\hat{\mathbf z}/r^2 = z_0/r^3$. So
$$\Omega = \int_0^{2\pi}\!\!\int_0^a\frac{z_0\,s}{(s^2 + z_0^2)^{3/2}}\,ds\,d\phi = 2\pi z_0\left[-\frac{1}{\sqrt{s^2 + z_0^2}}\right]_0^a = 2\pi\left(1 - \frac{z_0}{\sqrt{z_0^2 + a^2}}\right).$$
For $a = z_0 = 1$: $2\pi(1 - 1/\sqrt2) = 1.840302$. As $a\to\infty$, $\Omega\to2\pi$: an infinite plane hides half the sky. The project prints $5.657985$, $6.220357$, $6.276902$ for $a = 10, 100, 1000$ ($z_0 = 1$), approaching $2\pi = 6.283185$.

(c) The six faces of the unit cube centred at the origin are congruent and placed symmetrically about the centre, so each subtends the same solid angle; together they enclose the origin, so the six add to $4\pi$. Each face: $4\pi/6 = 2\pi/3 = 2.0944$. The top face is the square at $z_0 = \tfrac12$. Code:

```python
def square_solid_angle(z0, n=40):
    (u, wu), (v, wv) = gl(-0.5, 0.5, n), gl(-0.5, 0.5, n)
    U_, V_ = np.meshgrid(u, v, indexing="ij")
    return np.sum(np.outer(wu, wv) * z0 / (U_**2 + V_**2 + z0**2) ** 1.5)
```

The project gets $2.0943951024 = 2\pi/3$.

**C2.** (a) $\mathbb R^3$ is star-shaped, so by Theorem 3.16 $\mathbf G = -\nabla U_0$ for some $U_0$.
1. *The limit at infinity exists along each ray.* For $1 \le R_1 < R_2$, $|U_0(R_2\hat{\mathbf n}) - U_0(R_1\hat{\mathbf n})| = |\int_{R_1}^{R_2}\mathbf G\cdot\hat{\mathbf n}\,dt| \le \int_{R_1}^{R_2}K/t^2\,dt \le K/R_1$, which is Cauchy as $R_1\to\infty$.
2. *It doesn't depend on the direction.* Two points at radius $R$ are joined by a great-circle arc of length at most $\pi R$, along which $|\mathbf G| \le K/R^2$, so their $U_0$ values differ by at most $\pi K/R\to0$.
3. Call the common limit $c$ and set $U = U_0 - c$. Then $\mathbf G = -\nabla U$, $\nabla^2U = -\nabla\cdot\mathbf G = 0$, and for $r \ge 1$, $|U(r\hat{\mathbf n})| = \lim_{R\to\infty}|\int_r^R\mathbf G\cdot\hat{\mathbf n}\,dt| \le K/r$.

(b) Green's first identity (3.22) with $\varphi = \psi = U$ on the ball $B_R$, using $\nabla^2U = 0$:
$$\int_{B_R}|\mathbf G|^2\,d\tau = \oint_{r=R}U\,\nabla U\cdot d\mathbf a, \qquad \left|\oint_{r=R}U\,\nabla U\cdot d\mathbf a\right| \le \frac KR\cdot\frac{K}{R^2}\cdot4\pi R^2 = \frac{4\pi K^2}{R}.$$
The left side is nonnegative and increases with $R$; the right side tends to $0$. So $\int_{\mathbb R^3}|\mathbf G|^2\,d\tau = 0$, and since $\mathbf G$ is continuous, $\mathbf G = \mathbf 0$. Applied to $\mathbf G = \mathbf F_1 - \mathbf F_2$ for two solutions of the Helmholtz problem that fall off like $1/r^2$ (as (3.34) does for localized $D$, $\mathbf C$), this is the uniqueness half of Theorem 3.15.

(c) A constant field, e.g. $\mathbf G = (1, 0, 0)$: zero divergence and curl everywhere, nonzero. It violates the decay hypothesis. In the FFT decomposition (3.35) it is the $\mathbf k = \mathbf 0$ mode, the box average, which belongs to neither the gradient part nor the curl part (the project prints it for the random field).

**C3.** (a) $u' = -e^{-\mu r}(1/r^2 + \mu/r)$, so $r^2u' = -e^{-\mu r}(1 + \mu r)$. Differentiating: $\mu e^{-\mu r}(1 + \mu r) - \mu e^{-\mu r} = \mu^2re^{-\mu r}$. Dividing by $r^2$ as in (3.29): $\nabla^2u = \mu^2e^{-\mu r}/r = \mu^2u$ for $r \ne 0$. The project checks this symbolically.

(b) Let $\psi$ vanish for $r \ge R$ and $V_\epsilon = \{\epsilon < r < R\}$.
1. Green's second identity (3.23) with the pair $(u, \psi)$: $\int_{V_\epsilon}(u\nabla^2\psi - \psi\nabla^2u)\,d\tau = \oint_{\partial V_\epsilon}(u\nabla\psi - \psi\nabla u)\cdot d\mathbf a$.
2. On $V_\epsilon$, $\nabla^2u = \mu^2u$ by (a), so the left side is $\int_{V_\epsilon}u\,(\nabla^2 - \mu^2)\psi\,d\tau$.
3. On $r = R$, $\psi$ and $\nabla\psi$ vanish. On $r = \epsilon$, $d\mathbf a = -\hat{\mathbf r}\,\epsilon^2d\Omega$, so the right side is $-\oint(u\,\partial_r\psi - \psi\,\partial_ru)\,\epsilon^2\,d\Omega$.
4. As $\epsilon\to0$: $u\,\epsilon^2 = \epsilon e^{-\mu\epsilon}\to0$ and $\partial_r\psi$ is bounded, so the first term vanishes; $-\psi\,\partial_ru\,\epsilon^2 = \psi\,e^{-\mu\epsilon}(1 + \mu\epsilon)\to\psi(\mathbf 0)$. The right side tends to $-4\pi\psi(\mathbf 0)$.
5. The left side tends to the integral over all space, since $u\,r^2 = re^{-\mu r}$ is bounded. Hence $\int u\,(\nabla^2 - \mu^2)\psi\,d\tau = -4\pi\psi(\mathbf 0)$, i.e. $(\nabla^2 - \mu^2)u = -4\pi\delta^3(\mathbf r)$. ∎

Setting $\mu = 0$ recovers (3.33). This is the screened Coulomb potential: the Debye-screened potential of a charge in an electrolyte or plasma, and Yukawa's 1935 potential for the nuclear force (with $\mu$ the inverse range).

(c) For $\psi = e^{-r^2}$: $\psi' = -2re^{-r^2}$, $r^2\psi' = -2r^3e^{-r^2}$, and $\partial_r(r^2\psi') = (4r^4 - 6r^2)e^{-r^2}$, so $\nabla^2\psi = (4r^2 - 6)e^{-r^2}$. With $d\tau = 4\pi r^2dr$ for radial integrands:
$$\int u\,(\nabla^2 - \mu^2)\psi\,d\tau = 4\pi\int_0^\infty r\,e^{-\mu r}\,(4r^2 - 6 - \mu^2)\,e^{-r^2}\,dr,$$
which should equal $-4\pi\psi(0) = -4\pi$. The project gets $-12.5663706144 = -4\pi$ for $\mu = 1.5$.

**C4.** (a) Let $T_{ijk} = \varepsilon_{lmn}R_{li}R_{mj}R_{nk}$. Swapping $i$ and $j$ and then renaming the dummies $l \leftrightarrow m$: $T_{jik} = \varepsilon_{lmn}R_{lj}R_{mi}R_{nk} = \varepsilon_{mln}R_{mj}R_{li}R_{nk} = -\varepsilon_{lmn}R_{li}R_{mj}R_{nk} = -T_{ijk}$. The same argument works for any pair, so $T$ is totally antisymmetric and therefore $T_{ijk} = c\,\varepsilon_{ijk}$ (a totally antisymmetric array in three indices has one independent entry). Setting $(i,j,k) = (1,2,3)$: $c = \varepsilon_{lmn}R_{l1}R_{m2}R_{n3}$, the signed sum over permutations that defines the determinant of the matrix with columns $R_{\cdot1}, R_{\cdot2}, R_{\cdot3}$. So $c = \det R$.

(b) From §3.2, $\partial'_a = R_{aj}\partial_j$. *Divergence:* $\partial'_aF'_a = R_{aj}\partial_j(R_{ab}F_b) = (R^\top R)_{jb}\partial_jF_b = \partial_bF_b$.
*Curl:* $(\nabla'\times\mathbf F')_i = \varepsilon_{iab}\partial'_aF'_b = \varepsilon_{iab}R_{aj}R_{bc}\,\partial_jF_c$. To simplify $\varepsilon_{iab}R_{aj}R_{bc}$, write (a) as $\varepsilon_{lab}R_{lp}R_{aj}R_{bc} = (\det R)\varepsilon_{pjc}$, multiply by $R_{ip}$ and sum over $p$. Since $R_{lp}R_{ip} = (RR^\top)_{li} = \delta_{li}$, this gives $\varepsilon_{iab}R_{aj}R_{bc} = (\det R)R_{ip}\varepsilon_{pjc}$. Hence $(\nabla'\times\mathbf F')_i = (\det R)R_{ip}\varepsilon_{pjc}\partial_jF_c = (\det R)R_{ip}(\nabla\times\mathbf F)_p$. ∎

(c) $R = \operatorname{diag}(1, 1, -1)$, $\det R = -1$. The old coordinates are $\mathbf x = R^\top\mathbf x' = (x', y', -z')$, so $\mathbf F'(\mathbf x') = R\,(-y', x', 0) = (-y', x', 0)$ and $\nabla'\times\mathbf F' = (0, 0, 2)$. The formula: $(\det R)\,R\,(0, 0, 2) = -(0, 0, -2) = (0, 0, 2)$. ✓ A true vector would have become $R(0,0,2) = (0, 0, -2)$. Physically, the mirror $z\to-z$ leaves the turntable spinning counterclockwise as seen from above, because the motion lies in the $xy$-plane, which the mirror does not touch; but the right-hand-rule arrow we attach to that rotation "should" have been reflected to point down. Quantities defined through a cross product (angular velocity, angular momentum, $\mathbf B$) are pseudovectors. The project checks (b) symbolically for a rotation by $\arccos\tfrac35$ about $z$ and for this reflection, on two fields.
