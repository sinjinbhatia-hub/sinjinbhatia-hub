# F6 Solutions

Conventions as in the notes and exercises: $J[y] = \int_a^b L(x, y, y')\,dx$; $L_y$, $L_{y'}$ treat the three slots of $L$ as independent; $\frac{d}{dx}$ is the total derivative along the curve; $\delta J = \frac{d}{d\varepsilon}J[y + \varepsilon\eta]\big|_{\varepsilon=0}$ with $\eta(a) = \eta(b) = 0$. "E–L" is the Euler–Lagrange equation (3.4), "Case A/B/C" are the special cases of §3.6, "Beltrami" is (3.11). In the brachistochrone problems $y$ is measured downward and $g = 9.81\ \mathrm{m/s^2}$. Every number below is asserted in `project.py` (Parts 9 and 10 unless another part is named).

## A. Warm-up

**A1.** In each case compute $L_y$ and $L_{y'}$, then write $L_y - \frac{d}{dx}L_{y'} = 0$.

(a) $L = y'^2 + y^2$. $L_y = 2y$, $L_{y'} = 2y'$, $\frac{d}{dx}L_{y'} = 2y''$. E–L: $2y - 2y'' = 0$, i.e. $y'' = y$. General solution $y = C_1 e^{x} + C_2 e^{-x}$ (equivalently $A\cosh x + B\sinh x$). No special case is needed, though $L$ has no explicit $x$, so Beltrami gives the first integral $L - y'L_{y'} = y'^2 + y^2 - 2y'^2 = y^2 - y'^2 = \text{const}$, which is $y'' = y$ integrated once.

(b) $L = x\,y'^2$. No $y$ in the formula: **Case A**. $L_{y'} = 2x\,y' = C$, so $y' = C/(2x)$ and $y = \tfrac C2\ln|x| + D$; rename: $y = A\ln|x| + B$. Note that Beltrami does *not* apply here: $L$ contains $x$ explicitly, even though $L_y = 0$. (This is the warning in §8.)

(c) $L = y'^2 - 2y\sin x$. $L_y = -2\sin x$, $L_{y'} = 2y'$. E–L: $-2\sin x - 2y'' = 0$, i.e. $y'' = -\sin x$. Integrate twice: $y' = \cos x + C_1$, $y = \sin x + C_1 x + C_2$.

(d) $L = x^2y'^2 + y^2$. $L_y = 2y$, $L_{y'} = 2x^2y'$, $\frac{d}{dx}L_{y'} = 4xy' + 2x^2y''$. E–L: $2y - 4xy' - 2x^2y'' = 0$, i.e.
$$
x^2y'' + 2xy' - y = 0 .
$$
This is an Euler (equidimensional) equation. Trying $y = x^r$: $x^2\cdot r(r-1)x^{r-2} + 2x\cdot rx^{r-1} - x^r = x^r\,[r(r-1) + 2r - 1] = x^r\,(r^2 + r - 1)$, so $r^2 + r - 1 = 0$ and $r = \frac{-1 \pm\sqrt5}{2}$, i.e. $r_1 = 0.618$, $r_2 = -1.618$. General solution $y = C_1x^{r_1} + C_2x^{r_2}$ for $x > 0$.

Part 10 substitutes each general solution into its E–L equation with SymPy and gets $0$.

**A2.** $y = x^2$, $\eta = x(1 - x)$, so $y' = 2x$, $\eta' = 1 - 2x$, $y + \varepsilon\eta = x^2 + \varepsilon x(1 - x)$, $y' + \varepsilon\eta' = 2x + \varepsilon(1 - 2x)$.

(a) Expand the integrand of $\Phi(\varepsilon) = \int_0^1\big[(y' + \varepsilon\eta')^2 + 2x(y + \varepsilon\eta)\big]dx$:
$$
(2x + \varepsilon(1 - 2x))^2 = 4x^2 + 4\varepsilon x(1 - 2x) + \varepsilon^2(1 - 2x)^2, \qquad 2x\big(x^2 + \varepsilon x(1-x)\big) = 2x^3 + 2\varepsilon x^2(1 - x).
$$
Integrate term by term over $[0, 1]$:
- $\varepsilon^0$: $\int_0^1 4x^2\,dx + \int_0^1 2x^3\,dx = \tfrac43 + \tfrac12 = \tfrac{11}{6}$.
- $\varepsilon^1$: $\int_0^1 4x(1 - 2x)\,dx + \int_0^1 2x^2(1 - x)\,dx = 4\big(\tfrac12 - \tfrac23\big) + 2\big(\tfrac13 - \tfrac14\big) = -\tfrac23 + \tfrac16 = -\tfrac12$.
- $\varepsilon^2$: $\int_0^1(1 - 2x)^2\,dx = \int_0^1(1 - 4x + 4x^2)\,dx = 1 - 2 + \tfrac43 = \tfrac13$.

So $\Phi(\varepsilon) = \tfrac{11}{6} - \tfrac12\varepsilon + \tfrac13\varepsilon^2$ and $\delta J = \Phi'(0) = -\tfrac12$.

(b) $L = y'^2 + 2xy$ gives $L_y = 2x$ and $L_{y'} = 2y' = 4x$ on the curve. Formula (3.3):
$$
\delta J = \int_0^1\big[L_y\eta + L_{y'}\eta'\big]dx = \int_0^1\big[2x\cdot x(1 - x) + 4x(1 - 2x)\big]dx = \int_0^1\big(4x - 6x^2 - 2x^3\big)dx = 2 - 2 - \tfrac12 = -\tfrac12 . \quad\checkmark
$$

(c) E–L: $L_y - \frac{d}{dx}L_{y'} = 2x - 2y'' = 0$, i.e. $y'' = x$, with general solution $y = x^3/6 + C_1x + C_2$. For $y = x^2$, $y'' = 2 \ne x$: **not** a solution. It should not be: if $y$ satisfied E–L then by (3.5) $\delta J = \int(\text{E–L})\eta\,dx = 0$ for *every* $\eta$, and we found one $\eta$ with $\delta J = -\tfrac12$. (Also useful: $\Phi$ is a parabola with its minimum at $\varepsilon = \tfrac34$, so $y + \tfrac34\eta$ has a strictly smaller $J$ than $y$; $y = x^2$ is not even a local minimizer in this one direction.)

**A3.** Look at the formula for $L$, not at the curve.

(a) $L = \sqrt{1 + y'^2}\,/\,x$. Contains $x$ and $y'$, no $y$: **Case A** (Beltrami does not apply because $x$ is explicit). First integral (3.10): $L_{y'} = \dfrac{y'}{x\sqrt{1 + y'^2}} = C$. Solving it: $y'^2 = C^2x^2(1 + y'^2)$, $y' = \dfrac{Cx}{\sqrt{1 - C^2x^2}}$, $y = -\dfrac1C\sqrt{1 - C^2x^2} + D$, i.e. $x^2 + (y - D)^2 = 1/C^2$: circles centred on the $y$-axis (B1(a) with the roles of $x$ and $y$ swapped).

(b) $L = y^2(1 + y'^2)$. No explicit $x$: **Case C**. Beltrami: $L_{y'} = 2y^2y'$, so $L - y'L_{y'} = y^2(1 + y'^2) - 2y^2y'^2 = y^2(1 - y'^2) = C$.

(c) $L = y'^2 + xy$. Contains $x$, $y$ and $y'$: **none** of the special cases. Full E–L: $L_y = x$, $L_{y'} = 2y'$, so $x - 2y'' = 0$, $y'' = x/2$, $y = x^3/12 + C_1x + C_2$.

(d) $L = e^{y}\sqrt{1 + y'^2}$. No explicit $x$: **Case C**. $L_{y'} = e^y y'/\sqrt{1 + y'^2}$, so $L - y'L_{y'} = e^y\dfrac{(1 + y'^2) - y'^2}{\sqrt{1 + y'^2}} = \dfrac{e^y}{\sqrt{1 + y'^2}} = C$. (This is Fermat's principle (4.4) with $n(y) = e^y$.) Part 10 checks (b) and (d) symbolically.

**A4.** (a) $L = y\sqrt{1 + y'^2}$. $L_y = \sqrt{1 + y'^2}$ and $L_{y'} = \dfrac{y\,y'}{\sqrt{1 + y'^2}}$. Total derivative of $L_{y'}$ along the curve (product rule on $y\cdot y'$, chain rule on the root):
$$
\frac{d}{dx}L_{y'} = \frac{y'^2}{\sqrt{1 + y'^2}} + y\left[\frac{y''}{\sqrt{1 + y'^2}} - \frac{y'\cdot y'y''}{(1 + y'^2)^{3/2}}\right] = \frac{y'^2}{\sqrt{1 + y'^2}} + \frac{y\,y''}{(1 + y'^2)^{3/2}},
$$
because $y''(1 + y'^2) - y'^2y'' = y''$. E–L:
$$
\sqrt{1 + y'^2} - \frac{y'^2}{\sqrt{1 + y'^2}} - \frac{y\,y''}{(1 + y'^2)^{3/2}} = \frac{1}{\sqrt{1 + y'^2}} - \frac{y\,y''}{(1 + y'^2)^{3/2}} = 0 .
$$
Multiply by $(1 + y'^2)^{3/2}$: $(1 + y'^2) - y\,y'' = 0$, i.e. $y\,y'' = 1 + y'^2$. ∎ (Part 10: SymPy confirms $\mathrm{EL}\cdot(1 + y'^2)^{3/2} = 1 + y'^2 - yy''$.)

(b) $y = c\cosh u$ with $u = (x - x_0)/c$, $u' = 1/c$. Then $y' = c\sinh u\cdot\tfrac1c = \sinh u$ and $y'' = \cosh u\cdot\tfrac1c$. So $y\,y'' = c\cosh u\cdot\dfrac{\cosh u}{c} = \cosh^2u = 1 + \sinh^2u = 1 + y'^2$. ✓ The constant $c$ cancels, as it must: $yy''$ and $1 + y'^2$ are both unchanged when the curve is scaled.

**A5.** $L = \tfrac12y'^2 + y$, $L_y = 1$, $L_{y'} = y'$.

(a) E–L: $1 - y'' = 0$, so $y'' = 1$ and $y = \tfrac12x^2 + C_1x + C_2$.

(b) $y(0) = 0$ gives $C_2 = 0$. The natural boundary condition (3.15) at $x = 1$ is $L_{y'}(1) = y'(1) = 1 + C_1 = 0$, so $C_1 = -1$ and
$$
y = \tfrac12x^2 - x, \qquad y(1) = -\tfrac12 .
$$

(c) With $y' = x - 1$:
$$
J = \int_0^1\Big[\tfrac12(x - 1)^2 + \tfrac12x^2 - x\Big]dx = \int_0^1\big(x^2 - 2x + \tfrac12\big)dx = \tfrac13 - 1 + \tfrac12 = -\tfrac16 .
$$
With $y(1) = 0$ imposed instead: $\tfrac12 + C_1 = 0$, $C_1 = -\tfrac12$, $y = \tfrac12x(x - 1)$, $y' = x - \tfrac12$, and
$$
J = \int_0^1\Big[\tfrac12\big(x - \tfrac12\big)^2 + \tfrac12(x^2 - x)\Big]dx = \int_0^1\big(x^2 - x + \tfrac18\big)dx = \tfrac13 - \tfrac12 + \tfrac18 = -\tfrac1{24} .
$$
$-\tfrac16 < -\tfrac1{24}$. It must be: every function admissible in the fixed-end problem ($y(0) = 0$, $y(1) = 0$) is also admissible in the free-end problem ($y(0) = 0$ only), so the free-end minimum is taken over a larger set and cannot be larger. (Both are genuine minima: $L_{y'y'} = 1 > 0$ and $J$ is convex, being quadratic in $y'$ and linear in $y$.) Part 10 asserts both values.

## B. Core

**B1.** $L = n(y)\sqrt{1 + y'^2}$ has no explicit $x$, so Beltrami (4.4) gives $n(y)/\sqrt{1 + y'^2} = C$ along every ray.

(a) $n = 1/y$: $\dfrac{1}{y\sqrt{1 + y'^2}} = C$. Write $R := 1/C$; squaring, $y^2(1 + y'^2) = R^2$, so
$$
y'^2 = \frac{R^2 - y^2}{y^2}, \qquad y' = \pm\frac{\sqrt{R^2 - y^2}}{y} .
$$
Separate: $\dfrac{y\,dy}{\sqrt{R^2 - y^2}} = \pm dx$. The left side is $d\big(-\sqrt{R^2 - y^2}\big)$, so $-\sqrt{R^2 - y^2} = \pm(x - x_0)$. Square: $R^2 - y^2 = (x - x_0)^2$, i.e.
$$
(x - x_0)^2 + y^2 = R^2, \qquad y > 0 .
$$
A semicircle of radius $R = 1/C$ centred at $(x_0, 0)$ on the $x$-axis; the two signs are the rising and falling halves of the same circle. (Part 10 checks that the semicircle satisfies both $y^2(1 + y'^2) = R^2$ and the full E–L equation.) The limit $R \to\infty$ means $C \to 0$, i.e. $\sqrt{1 + y'^2} \to \infty$, i.e. $y' \to\infty$: the rays become **vertical lines** $x = x_0$. Semicircles centred on the boundary and vertical lines are exactly the geodesics of the Poincaré half-plane, whose metric is $ds = \sqrt{dx^2 + dy^2}/y$: the same functional.

(b) $n = n_0(1 + \kappa y)$: Beltrami gives $\dfrac{n_0(1 + \kappa y)}{\sqrt{1 + y'^2}} = C$, so $1 + y'^2 = \big(n_0(1 + \kappa y)/C\big)^2$. Put $w := 1 + \kappa y$, so $w' = \kappa y'$ and $y' = \sqrt{(n_0w/C)^2 - 1}$. Since $n_0w/C = \sqrt{1 + y'^2} \ge 1$ we may substitute $n_0w/C = \cosh u$; then $y' = \sinh u$, and differentiating the substitution, $n_0w'/C = \sinh u\cdot u'$, i.e. $n_0\kappa y'/C = \sinh u\cdot u'$, i.e. $n_0\kappa\sinh u/C = \sinh u\cdot u'$. Hence $u' = n_0\kappa/C$, $u = n_0\kappa(x - x_0)/C$, and
$$
1 + \kappa y = \frac{C}{n_0}\cosh\!\Big(\frac{n_0\kappa(x - x_0)}{C}\Big).
$$
(Part 10 substitutes this into the E–L equation of $L$ and gets $0$.) Curvature: $y'' = \dfrac{n_0\kappa}{C}\cosh(\cdots) > 0$ for $\kappa > 0$, so the ray is **concave up**, with its lowest point at $x = x_0$: it bends *toward* increasing $n$. Mirage: over a hot road the air at the surface is hottest, least dense and has the lowest index, so $n$ increases with height ($\kappa > 0$). A ray from the sky coming down toward the road flattens out and turns back up before touching it, then enters your eye travelling slightly *upward*, i.e. from a direction below the true horizon. Your eye extrapolates the final straight direction and places the image of the sky on the road surface, which therefore looks like a reflecting puddle. In the language of (4.5), $n\sin\theta = \text{const}$ with $\theta$ measured from the vertical: as the ray descends $n$ falls, so $\sin\theta$ must rise toward $1$, the ray turns horizontal, and then $\theta$ passes through $90°$ and the ray goes back up. Total internal reflection, done continuously.

**B2.** (a) With $z = z(\varphi)$ and $z' = dz/d\varphi$, $ds = \sqrt{R^2\,d\varphi^2 + dz^2} = \sqrt{R^2 + z'^2}\,d\varphi$, so $L = \sqrt{R^2 + z'^2}$. It contains neither $z$ nor $\varphi$: Case A. $L_{z'} = z'/\sqrt{R^2 + z'^2} = C$, so $z'^2 = C^2(R^2 + z'^2)$, $z'^2(1 - C^2) = C^2R^2$, and $z' = CR/\sqrt{1 - C^2} =: m$, a constant. So $z = m\varphi + z_0$: a **helix** (a circle when $m = 0$, with pitch $2\pi m$ per turn otherwise). Unrolling: the map $(\varphi, z) \mapsto (s, z) := (R\varphi, z)$ sends $ds^2 = R^2d\varphi^2 + dz^2$ to $ds^2 + dz^2$, the flat plane metric, so it preserves lengths; in the plane the shortest path is the straight line $z = (m/R)s + z_0$, which rolls back up into the helix. (Caution: on a cylinder two points are joined by helices with different winding numbers, all of which are extremals; only one is the shortest. Stationary is not minimum.)

(b) $L = \sqrt{\dot\theta^2 + \sin^2\theta\,\dot\varphi^2}$. Partial derivatives (treating $\theta, \varphi, \dot\theta, \dot\varphi$ as independent slots):
$$
L_\theta = \frac{\sin\theta\cos\theta\,\dot\varphi^2}{L}, \qquad L_{\dot\theta} = \frac{\dot\theta}{L}, \qquad L_\varphi = 0, \qquad L_{\dot\varphi} = \frac{\sin^2\theta\,\dot\varphi}{L}.
$$
The two E–L equations (3.13):
$$
\theta:\quad \frac{\sin\theta\cos\theta\,\dot\varphi^2}{L} - \frac{d}{dt}\Big(\frac{\dot\theta}{L}\Big) = 0, \qquad
\varphi:\quad -\frac{d}{dt}\Big(\frac{\sin^2\theta\,\dot\varphi}{L}\Big) = 0 .
$$
$L$ contains no $\varphi$ (Case A for the $\varphi$ equation), so the first integral is
$$
\frac{\sin^2\theta\,\dot\varphi}{L} = \text{const}.
$$
(With arc-length parametrization, $L = 1$, this is Clairaut's relation $\sin^2\theta\,\dot\varphi = \text{const}$; physically it is the $z$-component of angular momentum of a particle moving freely on the sphere, which is conserved because the sphere looks the same under rotation about the $z$-axis. CM-4 makes this Noether's theorem.) Equator: $\theta = \pi/2$, $\varphi = \omega t$, so $\dot\theta = 0$, $\sin\theta = 1$, $\cos\theta = 0$, $L = \sqrt{0 + 1\cdot\omega^2} = \omega$. The $\theta$ equation: $\dfrac{1\cdot 0\cdot\omega^2}{\omega} - \dfrac{d}{dt}\Big(\dfrac0\omega\Big) = 0$. ✓ The $\varphi$ equation: $\dfrac{d}{dt}\Big(\dfrac{1\cdot\omega}{\omega}\Big) = \dfrac{d}{dt}(1) = 0$. ✓ Part 10 does the same substitution in SymPy.

**B3.** (a) Let $\bar f := \dfrac{1}{b - a}\displaystyle\int_a^b f\,dx$ and define $\eta(x) := \displaystyle\int_a^x\big(f(t) - \bar f\big)dt$. Since $f$ is continuous, $\eta \in C^1[a, b]$ with $\eta' = f - \bar f$. Also $\eta(a) = 0$ and $\eta(b) = \int_a^b f\,dt - (b - a)\bar f = 0$. So $\eta$ is admissible, and the hypothesis gives
$$
0 = \int_a^b f\,\eta'\,dx = \int_a^b f\,(f - \bar f)\,dx .
$$
Also $\displaystyle\int_a^b\bar f\,(f - \bar f)\,dx = \bar f\int_a^b(f - \bar f)\,dx = \bar f\cdot 0 = 0$. Subtracting,
$$
\int_a^b (f - \bar f)^2\,dx = 0 .
$$
The integrand is continuous and $\ge 0$, so it is identically zero (if it were positive at one point it would be positive on an interval and the integral would be positive). Hence $f \equiv \bar f$, a constant. ∎

(b) Let $M(x) := \displaystyle\int_a^x L_y\big(t, y(t), y'(t)\big)dt$. Since $y \in C^1$ and $L \in C^2$, the integrand is continuous, so $M \in C^1$ with $M' = L_y$. Integrate the *first* term of (3.3) by parts:
$$
\int_a^b L_y\,\eta\,dx = \int_a^b M'\eta\,dx = \big[M\eta\big]_a^b - \int_a^b M\eta'\,dx = -\int_a^b M\eta'\,dx ,
$$
using $\eta(a) = \eta(b) = 0$. Then $\delta J = 0$ reads
$$
\int_a^b\big(L_{y'} - M\big)\eta'\,dx = 0 \quad\text{for every admissible } \eta .
$$
$L_{y'}(x, y, y') - M(x)$ is continuous (it involves only $y$ and $y'$, not $y''$), so by (a) it is a constant:
$$
L_{y'}\big(x, y(x), y'(x)\big) = \int_a^x L_y\big(t, y(t), y'(t)\big)dt + C . \qquad ∎
$$
This is (3.4) in integrated form: if $y$ happens to be $C^2$, differentiate both sides to get $\frac{d}{dx}L_{y'} = L_y$. What it buys: (i) **No $y''$ was assumed.** In §3.4 we needed $y \in C^2$ twice, to make $\frac{d}{dx}L_{y'}$ exist and to integrate the $L_{y'}\eta'$ term by parts; here we integrated the other term, which only needs $L_y$ continuous. (ii) **Regularity comes out instead of going in.** The right side is a $C^1$ function of $x$, so $L_{y'}(x, y(x), y'(x))$ is $C^1$ along the extremal even though $y'$ was only assumed continuous. Where $L_{y'y'} \ne 0$, the implicit function theorem then lets you solve $L_{y'}(x, y, y') = (\text{a } C^1 \text{ function})$ for $y'$ as a $C^1$ function, so $y \in C^2$ after all (Hilbert's differentiability theorem). The hypothesis of §3.4 is a conclusion. (iii) The same argument applied to piecewise-$C^1$ $y$ shows that $L_{y'}$ is continuous across a corner of $y$ even though $y'$ jumps there (the first Weierstrass–Erdmann corner condition). Gelfand & Fomin ch. 1 prove the lemma in (a) and derive (3.4) from it in exactly this way.

**B4.** (a) In $J_h = \frac h2\sum_{i=0}^N\big[L(x_i, y_i, p_i) + L(x_{i+1}, y_{i+1}, p_i)\big]$, the variable $y_i$ appears in the summands $i$ and $i - 1$, in four places:
- summand $i$, first term $L(x_i, y_i, p_i)$: in the $y$ slot and in $p_i$ (with $\partial p_i/\partial y_i = -1/h$);
- summand $i$, second term $L(x_{i+1}, y_{i+1}, p_i)$: in $p_i$ only;
- summand $i - 1$, first term $L(x_{i-1}, y_{i-1}, p_{i-1})$: in $p_{i-1}$ only (with $\partial p_{i-1}/\partial y_i = +1/h$);
- summand $i - 1$, second term $L(x_i, y_i, p_{i-1})$: in the $y$ slot and in $p_{i-1}$.

Chain rule on each:
$$
\frac{\partial J_h}{\partial y_i} = \frac h2\Big[L_y(x_i, y_i, p_i) - \tfrac1h L_{y'}(x_i, y_i, p_i) - \tfrac1h L_{y'}(x_{i+1}, y_{i+1}, p_i)\Big] + \frac h2\Big[\tfrac1h L_{y'}(x_{i-1}, y_{i-1}, p_{i-1}) + L_y(x_i, y_i, p_{i-1}) + \tfrac1h L_{y'}(x_i, y_i, p_{i-1})\Big].
$$
Define the **segment average** of $L_{y'}$ on segment $i$ (the one from $x_i$ to $x_{i+1}$, with slope $p_i$):
$$
\tilde Q_i := \tfrac12\Big[L_{y'}(x_i, y_i, p_i) + L_{y'}(x_{i+1}, y_{i+1}, p_i)\Big].
$$
Then the derivative regroups as
$$
\frac{\partial J_h}{\partial y_i} = h\left[\ \tfrac12\Big(L_y(x_i, y_i, p_i) + L_y(x_i, y_i, p_{i-1})\Big) - \frac{\tilde Q_i - \tilde Q_{i-1}}{h}\ \right].
$$

(b) Why centred, and why $O(h^3)$. Let $Q(x) := L_{y'}(x, y(x), y'(x))$ with $y$ smooth and $y_i = y(x_i)$. The slope $p_i = (y_{i+1} - y_i)/h$ is a *central* difference about the midpoint $x_{i+1/2}$, so $p_i = y'(x_{i+1/2}) + O(h^2)$; and the average of a smooth function at the two ends of an interval of length $h$ equals its value at the midpoint plus $O(h^2)$. Both facts together give $\tilde Q_i = Q(x_{i+1/2}) + h^2 r(x_{i+1/2}) + O(h^3)$, with $r$ a smooth function (it is built from derivatives of $L$ and $y$). Therefore
$$
\frac{\tilde Q_i - \tilde Q_{i-1}}{h} = \frac{Q(x_{i+1/2}) - Q(x_{i-1/2})}{h} + h^2\,\frac{r(x_{i+1/2}) - r(x_{i-1/2})}{h} + O(h^2) = Q'(x_i) + O(h^2),
$$
because the first quotient is a central difference of $Q$ about $x_i$ (error $O(h^2)$) and the second is $h^2$ times a bounded difference quotient. For the $L_y$ term: $p_i - y'(x_i) = +\tfrac h2y''(x_i) + O(h^2)$ and $p_{i-1} - y'(x_i) = -\tfrac h2y''(x_i) + O(h^2)$, so expanding $L_y(x_i, y_i, p)$ about $p = y'(x_i)$, the first-order terms cancel in the average and $\tfrac12\big(L_y(x_i, y_i, p_i) + L_y(x_i, y_i, p_{i-1})\big) = L_y(x_i, y(x_i), y'(x_i)) + O(h^2)$. Hence the bracket is $\mathrm{EL}(x_i) + O(h^2)$ and
$$
\frac{\partial J_h}{\partial y_i} = h\,\mathrm{EL}(x_i) + O(h^3).
$$
In (3.8), by contrast, $L_{y'}$ is evaluated at the *left* node of each segment with the *forward* slope, so the difference quotient is centred at $x_{i-1/2}$, not $x_i$, which costs one power of $h$: that is the $O(h^2)$ in (3.9).

(c) Reference implementation (Part 8b of `project.py`; same $L = y\sqrt{1 + y'^2} + xy^2$ and the same smooth non-solution $y = 1 + x^2/2 + \sin(2x)/3$ as Part 8):

```python
def J_trap(y_int, x):
    y = with_ends(y_int, y_f(x[0]), y_f(x[-1])); hh = x[1] - x[0]; p = np.diff(y) / hh
    return 0.5 * hh * np.sum(L_f(x[:-1], y[:-1], p) + L_f(x[1:], y[1:], p))

def grad_J_trap(y_int, x):
    y = with_ends(y_int, y_f(x[0]), y_f(x[-1])); hh = x[1] - x[0]; p = np.diff(y) / hh
    Q_seg = 0.5 * (Lp_f(x[:-1], y[:-1], p) + Lp_f(x[1:], y[1:], p))                 # Q~_i, one per segment
    Ly_avg = 0.5 * (Ly_f(x[1:-1], y[1:-1], p[1:]) + Ly_f(x[1:-1], y[1:-1], p[:-1]))  # one per interior node
    return hh * (Ly_avg - (Q_seg[1:] - Q_seg[:-1]) / hh)
```

The script first checks `grad_J_trap` against finite differences of `J_trap` (so the formula in (a) really is the gradient), then measures $\max_i|\partial J_h/\partial y_i - h\,\mathrm{EL}(x_i)|$ for $N = 50, 100, 200, 400$: it gets $3.82\times10^{-6}$, $4.95\times10^{-7}$, $6.29\times10^{-8}$, $7.94\times10^{-9}$, with successive ratios $7.72$, $7.86$, $7.93 \to 8$. Compare Part 8 with the left-Riemann discretization (3.7): $2.75\times10^{-4}$, $7.16\times10^{-5}$, $1.83\times10^{-5}$, $4.62\times10^{-6}$, ratios $3.83$, $3.92$, $3.96 \to 4$. At $N = 50$ the trapezoid gradient is already two orders of magnitude closer to the continuum E–L expression.

**B5.** $L = \tfrac12\dot y^2 - \tfrac12\omega^2y^2$, independent variable $t$.

(a) $L_y = -\omega^2y$, $L_{\dot y} = \dot y$. E–L: $-\omega^2y - \ddot y = 0$, i.e. $\ddot y + \omega^2y = 0$.

(b) Beltrami quantity: $H := L - \dot yL_{\dot y} = \tfrac12\dot y^2 - \tfrac12\omega^2y^2 - \dot y^2 = -\big(\tfrac12\dot y^2 + \tfrac12\omega^2y^2\big) = -E$, minus the total energy per unit mass (kinetic plus potential). $L$ has no explicit $t$, so $E$ is constant on solutions.

(c) From $\tfrac12\dot y^2 + \tfrac12\omega^2y^2 = E$ (with $E = -H \ge 0$): $\dot y = \pm\sqrt{2E - \omega^2y^2}$. Put $A := \sqrt{2E}/\omega$, so $2E = \omega^2A^2$ and $\dot y = \pm\omega\sqrt{A^2 - y^2}$. Separate:
$$
\frac{dy}{\sqrt{A^2 - y^2}} = \pm\omega\,dt \quad\Longrightarrow\quad \arcsin\frac yA = \pm\omega t + c \quad\Longrightarrow\quad y = A\sin(\pm\omega t + c).
$$
Since $\sin(-\omega t + c) = \sin(\omega t - c + \pi)$, both signs are of the form $y = A\sin(\omega t + \phi)$, with
$$
A = \frac{\sqrt{2E}}{\omega} = \frac{\sqrt{-2H}}{\omega}.
$$
Check: $\tfrac12\dot y^2 + \tfrac12\omega^2y^2 = \tfrac12A^2\omega^2(\cos^2 + \sin^2) = \tfrac12A^2\omega^2 = E$. ✓

(d) By (3.17) with $L_{yy} = -\omega^2$, $L_{y\dot y} = 0$, $L_{\dot y\dot y} = 1$:
$$
\delta^2S = \int_0^T\big(\dot\eta^2 - \omega^2\eta^2\big)dt .
$$
For $\eta = \sin(\pi t/T)$ (admissible: it vanishes at $0$ and $T$), $\dot\eta = (\pi/T)\cos(\pi t/T)$, and $\int_0^T\sin^2(\pi t/T)\,dt = \int_0^T\cos^2(\pi t/T)\,dt = T/2$ (one half-period of each). So
$$
\delta^2S = \Big(\frac\pi T\Big)^2\frac T2 - \omega^2\frac T2 = \frac T2\Big(\frac{\pi^2}{T^2} - \omega^2\Big) = \frac{\pi^2 - \omega^2T^2}{2T},
$$
which is **negative exactly when $T > \pi/\omega$**. Then $\Phi(\varepsilon) = S[y + \varepsilon\eta]$ has $\Phi'(0) = 0$ and $\Phi''(0) < 0$: the extremal is a local *maximum* of the action along this direction, hence not a minimum. (It is a saddle, not a maximum: for rapidly oscillating $\eta$ the $\dot\eta^2$ term dominates and $\delta^2S > 0$, as in the proof of Theorem 3.8.) $T = \pi/\omega$ is half a period, and it is precisely the conjugate point of Jacobi's condition in §3.11: every oscillator solution through $y(0) = 0$ returns to $y = 0$ at $t = \pi/\omega$, so neighbouring extremals through the start cross there. "Least action" is literally false for an oscillator watched for more than half a period; "stationary action" is what Hamilton's principle says. Part 10 evaluates $\delta^2S$ with SymPy.

**B6.** $J[y] = \int_0^1y'^2\,dx$, constraint $K[y] = \int_0^1y\,dx = 1$, $y(0) = y(1) = 0$. So $L = y'^2$ and $G = y$.

(a) $E_K = G_y - \frac{d}{dx}G_{y'} = 1 - 0 = 1 \ne 0$ for every $y$. No admissible function is an extremal of $K$, so the hypothesis of Theorem 3.7 holds automatically.

(b) $L^* = y'^2 + \lambda y$: $L^*_y = \lambda$, $L^*_{y'} = 2y'$. E–L: $\lambda - 2y'' = 0$, so $y'' = \lambda/2$ and $y = \tfrac\lambda4x^2 + C_1x + C_2$. $y(0) = 0$ gives $C_2 = 0$; $y(1) = 0$ gives $C_1 = -\lambda/4$. So $y = \tfrac\lambda4(x^2 - x) = \tfrac\lambda4\,x(x - 1)$. Constraint:
$$
\int_0^1\frac\lambda4(x^2 - x)\,dx = \frac\lambda4\Big(\frac13 - \frac12\Big) = -\frac{\lambda}{24} = 1 \quad\Longrightarrow\quad \lambda = -24,
$$
and $y = -6x(x - 1) = 6x(1 - x)$, a downward parabola with peak $3/2$ at $x = 1/2$.

(c) $y' = 6 - 12x$, so $J = \int_0^1(6 - 12x)^2\,dx = \int_0^1(36 - 144x + 144x^2)\,dx = 36 - 72 + 48 = 12$.

(d) Take any admissible competitor $\tilde y = y + \eta$ with $\eta(0) = \eta(1) = 0$ and $\int_0^1\eta\,dx = 0$ (so that $\tilde y$ also satisfies the constraint). Then
$$
J[\tilde y] = \int_0^1(y' + \eta')^2dx = J[y] + 2\int_0^1y'\eta'\,dx + \int_0^1\eta'^2\,dx .
$$
Integrate the cross term by parts: $2\int_0^1y'\eta'\,dx = 2[y'\eta]_0^1 - 2\int_0^1y''\eta\,dx = -2\cdot\frac\lambda2\int_0^1\eta\,dx = -\lambda\cdot 0 = 0$, using $y'' = \lambda/2$ and the constraint on $\eta$. So $J[\tilde y] = J[y] + \int_0^1\eta'^2\,dx \ge J[y]$, with equality only if $\eta' \equiv 0$, i.e. $\eta \equiv \eta(0) = 0$. The stationary point is the unique **global** minimizer. ∎ (The structure: $J$ is convex and the constraint set is affine, so the first-order condition is sufficient, the same reason the Lagrange multiplier conditions are sufficient for convex programs in ML-2.) Part 9 asserts $\int_0^1 6x(1 - x)\,dx = 1$ and $J = 12$; Part 10 checks that $6x(1 - x)$ solves the E–L equation of $y'^2 - 24y$.

**B7.** $L = \sqrt{1 + y'^2}/\sqrt{2gy}$, start $(0, 0)$, free arrival on the wall $x = 1$.

(a) $L_{y'} = \dfrac{y'}{\sqrt{1 + y'^2}\,\sqrt{2gy}}$. The natural boundary condition (3.15) at $x = 1$ is $L_{y'}(1) = 0$; since $y(1) > 0$ is finite, this forces $y'(1) = 0$: the wire must arrive at the wall **horizontally**.

(b) On the cycloid (4.1), $dy/dx = \sin\theta/(1 - \cos\theta) = \cot(\theta/2)$, which first vanishes at $\theta_1 = \pi$. Then $x_1 = a(\pi - \sin\pi) = \pi a = 1$ gives
$$
a = \frac1\pi = 0.3183\ \mathrm{m}, \qquad y_1 = a(1 - \cos\pi) = 2a = \frac2\pi = 0.6366\ \mathrm{m},
$$
so the wire arrives at the wall at the lowest point of its arch, $0.6366\ \mathrm{m}$ below the start. The time (4.3):
$$
T = \theta_1\sqrt{a/g} = \pi\sqrt{\frac{1}{\pi g}} = \sqrt{\frac\pi g} = \sqrt{\frac{\pi}{9.81}} = 0.5659\ \mathrm{s}.
$$

(c) Example 4.2 gives $T = 0.5829\ \mathrm{s}$ to the fixed point $(1, 1)$. The wall time $0.5659\ \mathrm{s}$ is **shorter**, and it must be: every wire from $(0, 0)$ to $(1, 1)$ is also a wire from $(0, 0)$ to the wall, so the wall problem minimizes over a larger set and its minimum cannot be larger. It is strictly smaller because the $(1, 1)$-cycloid arrives at $20.9°$ below horizontal, violating (3.15), so it is not even a stationary point of the wall problem. Equivalently: among all arrival depths $y_1$ on the wall, the descent time is minimized at $y_1 = 2/\pi$, and $y_1 = 1$ is past that optimum (the extra speed gained by going deeper no longer pays for the extra path). Part 9 asserts $\theta_1 = \pi$ and $T = \sqrt{\pi/g}$; Part 10 asserts $a$, $2a$, $T = 0.5659$ and $T < 0.5829$.

**B8.** Reference implementation (end of Part 10 in `project.py`). Each polygon segment, rotated about the $x$-axis, sweeps a conical frustum of exact lateral area $\pi(y_i + y_{i+1})\Delta s_i$, so the polygon functional $2\pi\sum_i\Delta s_i\,(y_i + y_{i+1})/2$ is the exact area of the rotated polygon, not an approximation.

```python
def polygon_area(y_int, x, ya, yb):
    y = with_ends(y_int, ya, yb)
    return float(2 * np.pi * np.sum(np.hypot(np.diff(x), np.diff(y)) * 0.5 * (y[1:] + y[:-1])))

def polygon_area_grad(y_int, x, ya, yb):
    y = with_ends(y_int, ya, yb); dy = np.diff(y); ds = np.hypot(np.diff(x), dy); ym = 0.5 * (y[1:] + y[:-1])
    return 2 * np.pi * ((dy[:-1] / ds[:-1]) * ym[:-1] + 0.5 * ds[:-1] - (dy[1:] / ds[1:]) * ym[1:] + 0.5 * ds[1:])

def catenoid_area(c, d):
    return 2 * np.pi * c * (d + 0.5 * c * np.sinh(2 * d / c))

Nk = 40; xk = np.linspace(-0.5, 0.5, Nk + 2)
res = minimize(polygon_area, np.ones(Nk), args=(xk, 1.0, 1.0), jac=polygon_area_grad,
               method="L-BFGS-B", bounds=[(0.0, None)] * Nk, options=dict(maxiter=5000, ftol=1e-15, gtol=1e-12))
```

The gradient: $y_i$ enters segments $i - 1$ and $i$; $\partial\Delta s_{i-1}/\partial y_i = \Delta y_{i-1}/\Delta s_{i-1}$, $\partial\Delta s_i/\partial y_i = -\Delta y_i/\Delta s_i$, and $\partial\bar y/\partial y_i = \tfrac12$ for both segment means. The script checks it against finite differences.

(a) *The two catenoids.* $y(\pm d) = 1$ with $y = c\cosh(x/c)$ requires $c\cosh(d/c) = 1$. The function $c \mapsto c\cosh(d/c)$ has a single minimum, at $u = d/c$ solving $u\tanh u = 1$ (Example 4.6), where it equals $d/0.6627$; for $d = 0.5$ that minimum is $0.7546 < 1$, so there are **two** roots: $c = 0.23509$ (thin neck) and $c = 0.84834$ (fat neck). The area of a catenoid: $y\sqrt{1 + y'^2} = c\cosh(x/c)\cdot\cosh(x/c)$, so
$$
A = 2\pi\int_{-d}^{d}c\cosh^2\frac xc\,dx = 2\pi c\int_{-d}^{d}\frac{1 + \cosh(2x/c)}{2}\,dx = 2\pi c\Big[\frac x2 + \frac c4\sinh\frac{2x}c\Big]_{-d}^{d} = 2\pi c\Big[d + \frac c2\sinh\frac{2d}c\Big].
$$
Numbers: thin neck $A = 6.8457$; fat neck $A = 5.9918$; two flat disks of radius $1$: $2\pi = 6.2832$. So $5.9918\ (\text{fat}) < 6.2832\ (\text{disks}) < 6.8457\ (\text{thin})$. Started from the cylinder $y \equiv 1$ (whose area is $2\pi\cdot 2d\cdot 1 = 2\pi$, the same as the two disks), the minimizer finds the **fat-neck** catenoid: area $5.9920$, maximum deviation $1.3\times10^{-5}$ from $0.84834\cosh(x/0.84834)$. Started exactly *on* the thin-neck catenoid it does not stay there: it slides to the same fat-neck solution with area $5.9920$. The thin catenoid satisfies E–L (it is an extremal) but is a saddle of the area functional, not a minimum. Stationary is not minimum, once more.

(b) *$d = 0.8$.* Now $\min_c c\cosh(0.8/c) = 0.8/0.6627 = 1.207 > 1$, so $c\cosh(d/c) = 1$ has **no** root and no catenoid spans the rings. The minimizer drives every interior node to the axis (the script asserts $\max_i y_i < 10^{-6}$, with the bound $y \ge 0$ active): the polygon becomes two vertical segments from the rings down to the axis plus a segment along the axis, with area $2\pi\sqrt{1 + h^2} = 6.2880 \to 2\pi$ as $h \to 0$. That is the Goldschmidt two-disk configuration of Example 4.6. It is not a smooth graph with $y > 0$, so it is not an admissible function: E–L has no solution in the admissible class, the infimum of the area is $2\pi$, and it is attained only outside the class. A soap film between rings pulled this far apart breaks into two flat disks.

## C. Challenge

**C1.** (a) $f(p) = \sqrt{1 + p^2}$: $f'(p) = p/\sqrt{1 + p^2}$, $f''(p) = (1 + p^2)^{-3/2} > 0$, so $f$ is strictly convex. Taylor's theorem with Lagrange remainder: $f(p + q) = f(p) + f'(p)q + \tfrac12f''(\xi)q^2$ for some $\xi$ between $p$ and $p + q$, and the remainder is $\ge 0$, with equality iff $q = 0$. So $f(p + q) \ge f(p) + f'(p)q$ for all real $p, q$. ∎

(b) Let $y$ be the line, $y' \equiv m = (B - A)/(b - a)$, and $\tilde y = y + \eta$ with $\eta(a) = \eta(b) = 0$, so $\tilde y' = m + \eta'$. Apply (a) pointwise with $p = m$, $q = \eta'(x)$ and integrate:
$$
J[\tilde y] = \int_a^b f(m + \eta')\,dx \ \ge\ \int_a^b\big[f(m) + f'(m)\eta'\big]dx = (b - a)f(m) + f'(m)\big[\eta(b) - \eta(a)\big] = J[y].
$$
(Check: $(b - a)\sqrt{1 + m^2} = \sqrt{(b - a)^2 + (B - A)^2}$, the distance between the endpoints.) Equality requires $\tfrac12f''(\xi)\eta'^2 = 0$ at every $x$ (the integrand of the difference is continuous and $\ge 0$), i.e. $\eta' \equiv 0$, i.e. $\eta \equiv \eta(a) = 0$, i.e. $\tilde y = y$. The line is the unique **global** minimizer. ∎

(c) Let $L = L(x, p)$ be convex in $p$ for each $x$, so $L(x, p + q) \ge L(x, p) + L_p(x, p)\,q$. Along an extremal, Case A gives $L_p(x, y'(x)) = C$, a constant. For any $\tilde y = y + \eta$ with the same endpoints,
$$
J[\tilde y] = \int_a^b L(x, y' + \eta')\,dx \ \ge\ \int_a^b\big[L(x, y') + L_p(x, y')\eta'\big]dx = J[y] + C\int_a^b\eta'\,dx = J[y] + C\big[\eta\big]_a^b = J[y]. \qquad ∎
$$
The step that makes it work: because $L_p$ is *constant* along the extremal, the linear term integrates to a boundary term, which the endpoint conditions kill. Convexity then makes the remainder nonnegative.

(d) Two things fail for $L = \sqrt{1 + y'^2}/\sqrt{2gy}$. First, $L$ depends on $y$, so (c) does not apply: $L_{y'}$ is not constant along the cycloid (the first integral is $y(1 + y'^2) = 2a$, which involves $y$), and the linear term $\int(L_y\eta + L_{y'}\eta')\,dx$ vanishes only after integrating by parts and using E–L, not by itself. One could try to rescue the argument with the *joint* convexity of $L$ in $(y, y')$, which would give $L(y + \eta, y' + \eta') \ge L + L_y\eta + L_{y'}\eta'$ pointwise and then $J[\tilde y] \ge J[y]$ after the integration by parts. But second, $L$ is **not** jointly convex: with $p = y'$,
$$
\det\begin{pmatrix}L_{yy} & L_{yp}\\ L_{yp} & L_{pp}\end{pmatrix} = \frac{3 - p^2}{8g\,y^3(1 + p^2)},
$$
(Part 10 computes this Hessian with SymPy), which is negative wherever $|y'| > \sqrt3$, i.e. wherever the wire is steeper than $60°$. Every brachistochrone starts vertically, so every one has a stretch where the Hessian is indefinite, and the pointwise inequality is false there. (Also $1/\sqrt y$ blows up at the start.) The cycloid *is* the global minimizer, but proving it needs the field theory of Weierstrass (Gelfand & Fomin ch. 6), not pointwise convexity.

**C2.** (a) The endpoint $(X, 0)$ has $y_1 = 0$, so (4.1) requires $a(1 - \cos\theta_1) = 0$ with $\theta_1 > 0$: $\theta_1 = 2\pi$, a full arch. Then $X = a(2\pi - \sin 2\pi) = 2\pi a$, so
$$
a = \frac{X}{2\pi}, \qquad \text{maximum depth } 2a = \frac X\pi \ \text{(reached at } \theta = \pi,\ x = X/2), \qquad T = \theta_1\sqrt{a/g} = 2\pi\sqrt{\frac{X}{2\pi g}} = \sqrt{\frac{2\pi X}{g}}.
$$
For $X = 1\ \mathrm{m}$: $a = 0.1592\ \mathrm{m}$, depth $0.3183\ \mathrm{m}$, $T = \sqrt{2\pi/9.81} = 0.8003\ \mathrm{s}$. A horizontal wire at depth $0$ never gets there: $v = \sqrt{2gy} = 0$ all along it, so a bead released from rest does not move at all ($T = \infty$; the integrand $1/\sqrt{2gy}$ is infinite everywhere). Every wire that arrives must first dip to gain speed; the cycloid is the optimal trade between dipping and distance. Note the scaling $T \propto\sqrt X$: the fastest time to a point at the same height grows like the square root of the distance, exactly as for free fall through a height $\propto X$.

(b) $f(\theta) = \dfrac{\theta - \sin\theta}{1 - \cos\theta}$ on $(0, 2\pi)$. Quotient rule:
$$
f'(\theta) = \frac{(1 - \cos\theta)(1 - \cos\theta) - (\theta - \sin\theta)\sin\theta}{(1 - \cos\theta)^2} = \frac{N(\theta)}{(1 - \cos\theta)^2}, \qquad N = 1 - 2\cos\theta + \cos^2\theta - \theta\sin\theta + \sin^2\theta = 2 - 2\cos\theta - \theta\sin\theta .
$$
Half angles, with $u := \theta/2 \in (0, \pi)$: $2 - 2\cos\theta = 4\sin^2u$ and $\theta\sin\theta = 2u\cdot 2\sin u\cos u$, so
$$
N = 4\sin u\,(\sin u - u\cos u).
$$
On $(0, \pi)$, $\sin u > 0$. Let $g(u) := \sin u - u\cos u$: $g(0) = 0$ and $g'(u) = \cos u - \cos u + u\sin u = u\sin u > 0$ on $(0, \pi)$, so $g > 0$ there. Hence $N > 0$ and $f' > 0$: $f$ is strictly increasing on $(0, 2\pi)$. Limits: as $\theta \to 0^+$, $\theta - \sin\theta = \theta^3/6 + O(\theta^5)$ and $1 - \cos\theta = \theta^2/2 + O(\theta^4)$, so $f \approx \theta/3 \to 0$; as $\theta \to 2\pi^-$, the numerator tends to $2\pi$ and the denominator to $0^+$, so $f \to +\infty$. A continuous strictly increasing function from $0$ to $\infty$ is a bijection $(0, 2\pi) \to (0, \infty)$, so (4.2) has exactly one solution $\theta_1$ for every $x_1/y_1 > 0$. ∎ (Part 10 checks monotonicity on a grid of $5000$ points and $f(\pi) = \pi/2$, the dip threshold of Example 4.2.)

**C3.** (a) Let $\eta(x, y)$ be $C^1$ on $\bar\Omega$ with $\eta = 0$ on $\partial\Omega$, and $\Phi(\varepsilon) = J[u + \varepsilon\eta]$. Differentiating under the integral with the chain rule on the three slots that depend on $\varepsilon$ (the slots $u, u_x, u_y$ receive $\eta, \eta_x, \eta_y$):
$$
\delta J = \Phi'(0) = \iint_\Omega\big[L_u\,\eta + L_{u_x}\eta_x + L_{u_y}\eta_y\big]dx\,dy .
$$
Define the vector field $\mathbf F := (L_{u_x}, L_{u_y})$ evaluated along $u$. Then $L_{u_x}\eta_x + L_{u_y}\eta_y = \mathbf F\cdot\nabla\eta = \nabla\cdot(\eta\mathbf F) - \eta\,\nabla\cdot\mathbf F$ (product rule for the divergence). By the divergence theorem,
$$
\iint_\Omega\nabla\cdot(\eta\mathbf F)\,dx\,dy = \oint_{\partial\Omega}\eta\,\mathbf F\cdot\mathbf n\,ds = 0,
$$
because $\eta = 0$ on $\partial\Omega$. This is the two-dimensional integration by parts, and the boundary term dies for the same reason as in §3.4. So
$$
\delta J = \iint_\Omega\Big[L_u - \frac{\partial}{\partial x}L_{u_x} - \frac{\partial}{\partial y}L_{u_y}\Big]\eta\,dx\,dy = 0 \quad\text{for every admissible } \eta .
$$
**Fundamental lemma in two dimensions.** *If $f$ is continuous on $\bar\Omega$ and $\iint_\Omega f\eta\,dx\,dy = 0$ for every $\eta \in C^1(\bar\Omega)$ vanishing on $\partial\Omega$, then $f \equiv 0$ on $\Omega$.* Proof: if $f(P) > 0$ at an interior point $P$, continuity gives a closed disk $D_r(P) \subset \Omega$ on which $f > \tfrac12f(P)$. Take the bump $\eta(X) = (r^2 - |X - P|^2)^2$ on $D_r(P)$ and $0$ outside; it is $C^1$ (the square kills the kink at the rim), vanishes on $\partial\Omega$, and is positive inside the disk, so $\iint f\eta \ge \tfrac12f(P)\iint_{D_r}\eta > 0$, a contradiction. ($f < 0$: replace $f$ by $-f$.) Applying the lemma to the bracket (continuous when $L \in C^2$, $u \in C^2$) gives
$$
\frac{\partial L}{\partial u} - \frac{\partial}{\partial x}\frac{\partial L}{\partial u_x} - \frac{\partial}{\partial y}\frac{\partial L}{\partial u_y} = 0 . \qquad ∎
$$

(b) $L = \tfrac12(u_x^2 + u_y^2)$: $L_u = 0$, $L_{u_x} = u_x$, $L_{u_y} = u_y$, so the equation is $-(u_{xx} + u_{yy}) = 0$: **Laplace's equation** $\nabla^2u = 0$. (Part 10 does this substitution symbolically.) $L = \tfrac12(u_x^2 + u_y^2) - \rho u$: $L_u = -\rho$, so $-\rho - \nabla^2u = 0$: **Poisson's equation** $\nabla^2u = -\rho$. In EM-2 these are the equations for the electrostatic potential $V$: $\nabla^2V = -\rho/\varepsilon_0$ where there is charge density $\rho$ and $\nabla^2V = 0$ where there is none; the functional is the field energy $\tfrac{\varepsilon_0}2\iint|\nabla V|^2 = \tfrac{\varepsilon_0}2\iint|\mathbf E|^2$ minus the interaction energy $\iint\rho V$ of the charges with the potential. So the electrostatic potential with prescribed values on the conductors is the one that makes the energy stationary (and, since the Dirichlet energy is convex, minimal): Dirichlet's principle, which is also why finite-element solvers for electrostatics and for the beam equation of §3.8 are minimizers. EM-3 is about solving these two equations with given boundary values.

**C4.** (a) Three isoperimetric constraints, $\int q = 1$, $\int xq = \mu$, $\int x^2q = \sigma^2 + \mu^2$, so three multipliers by the same argument as Theorem 3.7 with three fixed directions $\eta_2, \eta_3, \eta_4$ and the implicit function theorem in four variables. The augmented Lagrangian
$$
L^* = -q\ln q + \lambda_0q + \lambda_1xq + \lambda_2x^2q
$$
contains $q$ and $x$ but **not $q'$**, so $L^*_{q'} = 0$ and the E–L equation $L^*_q - \frac{d}{dx}L^*_{q'} = 0$ collapses to $L^*_q = 0$ (Case B): an algebraic equation at each $x$, with no ODE and no boundary conditions,
$$
-\ln q - 1 + \lambda_0 + \lambda_1x + \lambda_2x^2 = 0 .
$$

(b) Solve: $q(x) = \exp\big(\lambda_0 - 1 + \lambda_1x + \lambda_2x^2\big)$. For $q$ to be a density it must be integrable over $\mathbb R$, which requires **$\lambda_2 < 0$** (for $\lambda_2 > 0$ the exponent $\to +\infty$ at both ends; for $\lambda_2 = 0$ it is an exponential, infinite on one side, or a constant, not normalizable). Complete the square:
$$
\lambda_1x + \lambda_2x^2 = \lambda_2\Big(x + \frac{\lambda_1}{2\lambda_2}\Big)^2 - \frac{\lambda_1^2}{4\lambda_2},
$$
so $q$ is a Gaussian with mean $-\lambda_1/(2\lambda_2)$ and variance $-1/(2\lambda_2)$. Matching the constraints: $\sigma^2 = -1/(2\lambda_2)$ and $\mu = -\lambda_1/(2\lambda_2)$, i.e.
$$
\lambda_2 = -\frac{1}{2\sigma^2}, \qquad \lambda_1 = \frac{\mu}{\sigma^2}.
$$
Normalization fixes $\lambda_0$: $q = \dfrac{1}{\sigma\sqrt{2\pi}}\exp\Big(-\dfrac{(x - \mu)^2}{2\sigma^2}\Big)$ requires $\exp\big(\lambda_0 - 1 - \lambda_1^2/(4\lambda_2)\big) = 1/(\sigma\sqrt{2\pi})$. With $\lambda_1^2/(4\lambda_2) = \dfrac{\mu^2/\sigma^4}{-2/\sigma^2} = -\dfrac{\mu^2}{2\sigma^2}$, this is $\lambda_0 - 1 + \mu^2/(2\sigma^2) = -\tfrac12\ln(2\pi\sigma^2)$, so
$$
\lambda_0 = 1 - \frac{\mu^2}{2\sigma^2} - \frac12\ln\big(2\pi\sigma^2\big).
$$
The mean and variance constraints are then met automatically (they are the Gaussian's moments). The maximum entropy is $H = \tfrac12\ln(2\pi e\sigma^2)$. Part 10 checks with SymPy that this $q$ with these three multipliers satisfies the algebraic equation identically, has moments $1, \mu, \mu^2 + \sigma^2$, and has that entropy.

(c) $\partial^2L^*/\partial q^2 = -1/q < 0$ for every $x$ (since $q > 0$). The second variation here is $\delta^2 = \int L^*_{qq}\,\eta^2\,dx = -\int\eta^2/q\,dx < 0$ for every $\eta \not\equiv 0$: in Legendre's language (Theorem 3.8, "for a maximum, $\le 0$") the stationary density is a strict **maximum**, not a minimum. Because $-q\ln q$ is strictly concave in $q$ and all three constraints are linear in $q$, it is the global maximum over the constraint set: the Gaussian is *the* maximum-entropy density with given mean and variance (F5's entropy module states this; this is the proof).

(d) **SM-4** does the identical calculation with different constraints: maximize $S = -k_B\sum_ip_i\ln p_i$ subject to $\sum_ip_i = 1$ and $\sum_ip_iE_i = U$. Two multipliers, an algebraic stationarity condition $-\ln p_i - 1 + \lambda_0 - \beta E_i = 0$, and the solution $p_i \propto e^{-\beta E_i}$ is the **Boltzmann distribution**, with the energy constraint's multiplier $\beta$ identified as $1/(k_BT)$. That is Jaynes's 1957 "stat mech as inference" view, which SM-4 pairs with: the canonical ensemble is the least-committal distribution consistent with a known mean energy, exactly as the Gaussian is the least-committal density consistent with a known mean and variance. **ML-10**'s variational inference maximizes the ELBO $\mathcal L[q] = \mathbb E_q[\ln p(x, z)] - \mathbb E_q[\ln q(z)] = \mathbb E_q[\ln p(x, z)] + H[q]$ over densities $q(z)$: the same entropy functional, with the moment constraints replaced by the expected log-joint. Mean-field VI sets the functional derivative (3.18) of the ELBO plus a multiplier for $\int q_j = 1$ to zero; since no $q_j'$ appears, the condition is algebraic (Case B again) and gives $\ln q_j^*(z_j) = \mathbb E_{-j}[\ln p(x, z)] + \text{const}$, i.e. $q_j^* \propto \exp\big(\mathbb E_{-j}[\ln p(x, z)]\big)$. When $\ln p$ is quadratic in $z_j$, as in Gaussian models and in the VAE's Gaussian encoder, this is a Gaussian, which is why a Gaussian variational family is the default: by (b), the exponential of a quadratic is the only density a quadratic set of constraints can produce.
