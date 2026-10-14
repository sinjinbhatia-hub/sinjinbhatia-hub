# CM-1 Solutions

## A. Warm-up

**A1.** $M = 2 + 1 + 3 = 6$ kg.

$\mathbf{R} = \frac{1}{6}\left[2(1, 0, 0) + 1(0, 2, 0) + 3(0, 0, -1)\right] = \frac{1}{6}(2, 2, -3) = (0.3333, 0.3333, -0.5)$ m.

$\mathbf{P} = 2(0, 3, 0) + 1(-1, 0, 0) + 3(2, 1, 0) = (0, 6, 0) + (-1, 0, 0) + (6, 3, 0) = (5, 9, 0)$ kg m/s.

$\dot{\mathbf{R}} = \mathbf{P}/M = (5/6, 3/2, 0) = (0.8333, 1.5, 0)$ m/s, and $M\dot{\mathbf{R}} = (5, 9, 0) = \mathbf{P}$. ✓

Angular momentum, term by term with $\mathbf{p}_i = m_i\mathbf{v}_i$:
$\mathbf{r}_1\times\mathbf{p}_1 = (1, 0, 0)\times(0, 6, 0) = (0\cdot 0 - 0\cdot 6,\ 0\cdot 0 - 1\cdot 0,\ 1\cdot 6 - 0\cdot 0) = (0, 0, 6)$.
$\mathbf{r}_2\times\mathbf{p}_2 = (0, 2, 0)\times(-1, 0, 0) = (2\cdot 0 - 0\cdot 0,\ 0\cdot(-1) - 0\cdot 0,\ 0\cdot 0 - 2\cdot(-1)) = (0, 0, 2)$.
$\mathbf{r}_3\times\mathbf{p}_3 = (0, 0, -1)\times(6, 3, 0) = (0\cdot 0 - (-1)\cdot 3,\ (-1)\cdot 6 - 0\cdot 0,\ 0\cdot 3 - 0\cdot 6) = (3, -6, 0)$.
$\mathbf{L} = (3, -6, 8)$ kg m$^2$/s.

**A2.** Use $\nabla\times\mathbf{F} = (\partial_yF_z - \partial_zF_y,\ \partial_zF_x - \partial_xF_z,\ \partial_xF_y - \partial_yF_x)$.

(a) $\mathbf{F} = (y, x, 0)$: $z$-component $\partial_x(x) - \partial_y(y) = 1 - 1 = 0$; the other two involve only $\partial_z$ of $F_x, F_y$ (zero) and derivatives of $F_z = 0$. Curl-free on all of $\mathbb{R}^3$ (simply connected), so conservative. $V$: need $-\partial_xV = y$, so $V = -xy + f(y, z)$; $-\partial_yV = x - f_y = x$ gives $f_y = 0$; $-\partial_zV = -f_z = 0$. $V = -xy$.

(b) $\mathbf{F} = (-y, x, 0)$: $z$-component $\partial_x(x) - \partial_y(-y) = 1 + 1 = 2 \neq 0$. Not conservative. Directly: around a circle of radius $r$ counterclockwise, $\mathbf{F}$ is tangent with magnitude $r$, so $\oint\mathbf{F}\cdot d\mathbf{r} = r\cdot 2\pi r = 2\pi r^2 \neq 0$.

(c) $\mathbf{F} = (2xy, x^2 + z^2, 2yz)$: $x$-component $\partial_y(2yz) - \partial_z(x^2 + z^2) = 2z - 2z = 0$; $y$-component $\partial_z(2xy) - \partial_x(2yz) = 0 - 0 = 0$; $z$-component $\partial_x(x^2 + z^2) - \partial_y(2xy) = 2x - 2x = 0$. Conservative. $V$: $-\partial_xV = 2xy \Rightarrow V = -x^2y + f(y, z)$; $-\partial_yV = x^2 - f_y = x^2 + z^2 \Rightarrow f_y = -z^2 \Rightarrow f = -yz^2 + h(z)$; $-\partial_zV = 2yz - h' = 2yz \Rightarrow h' = 0$. $V = -(x^2y + yz^2)$.

(d) $\mathbf{F} = -b\mathbf{v}$ depends on velocity, not position, so it is not a force field $\mathbf{F}(\mathbf{r})$ at all and Theorem 3.3 does not apply. It is not conservative in any sense: along any closed path traversed at nonzero speed, $\oint\mathbf{F}\cdot d\mathbf{r} = -b\int|\mathbf{v}|^2dt < 0$.

**A3.** $s_\infty = \sqrt{mg/c} = \sqrt{0.145\times 9.81/1.225\times 10^{-3}} = \sqrt{1161.2} = 34.08$ m/s. $\tau = s_\infty/g = 3.474$ s. From $s = s_\infty\tanh(t/\tau)$, 90% of terminal when $\tanh(t/\tau) = 0.9$, $t = \tau\operatorname{artanh}(0.9) = 3.474\times 1.4722 = 5.114$ s. With no drag a ball dropped for 5 s would be at $gt = 49.05$ m/s, already above the terminal speed; drag matters for a baseball after about a second.

**A4.** $V'(x) = V_0\left(\dfrac{1}{a} - \dfrac{a}{x^2}\right) = 0 \Rightarrow x^2 = a^2 \Rightarrow x_0 = a$ (the only root with $x > 0$). $V''(x) = 2V_0a/x^3$, so $V''(a) = 2V_0/a^2 > 0$: stable. By (3.8), $\omega^2 = V''(a)/m = 2V_0/(ma^2)$, $\omega = \sqrt{2V_0/(ma^2)}$.

Harmonic check of (3.7): with $E = \tfrac12kA^2$, turning points $\pm A$,
$$T = \sqrt{2m}\int_{-A}^{A}\frac{dx}{\sqrt{\tfrac12k(A^2 - x^2)}} = \sqrt{2m}\sqrt{\frac{2}{k}}\int_{-A}^{A}\frac{dx}{\sqrt{A^2 - x^2}} = 2\sqrt{\frac{m}{k}}\left[\arcsin\frac{x}{A}\right]_{-A}^{A} = 2\sqrt{\frac{m}{k}}\cdot\pi = 2\pi\sqrt{\frac{m}{k}}.$$
Independent of $A$, as it must be.

**A5.** From (3.16) with $m_1 = 0.8$, $m_2 = 1.2$: $\ddot{y}_1 = \dfrac{0.4}{2.0}\times 9.81 = 1.962$ m/s$^2$ (the 0.8 kg mass rises), $T = \dfrac{2\times 0.8\times 1.2}{2.0}\times 9.81 = 0.96\times 9.81 = 9.4176$ N.

Limits. $m_1 = m_2 = m$: $\ddot{y}_1 = 0$, $T = mg$. The masses balance; the string holds each one's weight. $m_2 \gg m_1$: $\ddot{y}_1 \to g$ and $T \to 2m_1g$. The heavy mass falls freely, dragging the light one up at $g$, which needs $T - m_1g = m_1g$.

Axle force: the pulley is massless so the net force on it is zero; the string pulls down on it with $T$ on each side, so the axle supplies $2T = 18.835$ N upward. This is less than the total weight $(m_1 + m_2)g = 19.62$ N by $0.785$ N, which is $M\ddot{Y}$ for the center of mass: $\ddot{Y} = (m_1\ddot{y}_1 + m_2\ddot{y}_2)/M = (0.8 - 1.2)(1.962)/2 = -0.392$ m/s$^2$, and $2.0\times 0.392 = 0.785$ N. ✓ The CM accelerates downward, so the support is less than the weight (Theorem 3.1).

## B. Core

**B1.** $\sum_im_i\mathbf{r}_i' = \sum_im_i\mathbf{r}_i - \left(\sum_im_i\right)\mathbf{R} = M\mathbf{R} - M\mathbf{R} = \mathbf{0}$. Differentiating, $\sum_im_i\dot{\mathbf{r}}_i' = \mathbf{0}$ as well.

(a) $\mathbf{L} = \sum_i(\mathbf{R} + \mathbf{r}_i')\times m_i(\dot{\mathbf{R}} + \dot{\mathbf{r}}_i')$. Expand into four sums:
$\sum_im_i\,\mathbf{R}\times\dot{\mathbf{R}} = \mathbf{R}\times M\dot{\mathbf{R}} = \mathbf{R}\times\mathbf{P}$;
$\mathbf{R}\times\sum_im_i\dot{\mathbf{r}}_i' = \mathbf{R}\times\mathbf{0} = \mathbf{0}$;
$\left(\sum_im_i\mathbf{r}_i'\right)\times\dot{\mathbf{R}} = \mathbf{0}\times\dot{\mathbf{R}} = \mathbf{0}$;
$\sum_i\mathbf{r}_i'\times m_i\dot{\mathbf{r}}_i'$, the angular momentum about the CM.
So $\mathbf{L} = \mathbf{R}\times\mathbf{P} + \sum_i\mathbf{r}_i'\times m_i\dot{\mathbf{r}}_i'$. ∎

(b) $T = \tfrac12\sum_im_i(\dot{\mathbf{R}} + \dot{\mathbf{r}}_i')\cdot(\dot{\mathbf{R}} + \dot{\mathbf{r}}_i') = \tfrac12M\dot{\mathbf{R}}^2 + \dot{\mathbf{R}}\cdot\sum_im_i\dot{\mathbf{r}}_i' + \tfrac12\sum_im_i\dot{\mathbf{r}}_i'^2 = \tfrac12M\dot{\mathbf{R}}^2 + 0 + \tfrac12\sum_im_i\dot{\mathbf{r}}_i'^2$. ∎

The cross terms vanish in both because $\sum_im_i\mathbf{r}_i' = \mathbf{0}$ and its time derivative is zero: that is what "center of mass" means. Nothing about forces was used.

**B2.** (a) $F_x = -y/(x^2 + y^2)$, $F_y = x/(x^2 + y^2)$, $F_z = 0$, none depending on $z$. The $x$- and $y$-components of the curl involve $\partial_zF_x$, $\partial_zF_y$ (zero) and derivatives of $F_z$ (zero). The $z$-component:
$$\partial_xF_y = \frac{(x^2 + y^2) - x\cdot 2x}{(x^2 + y^2)^2} = \frac{y^2 - x^2}{(x^2 + y^2)^2}, \qquad \partial_yF_x = \frac{-(x^2 + y^2) + y\cdot 2y}{(x^2 + y^2)^2} = \frac{y^2 - x^2}{(x^2 + y^2)^2}.$$
Their difference is $0$. So $\nabla\times\mathbf{F} = \mathbf{0}$ wherever $x^2 + y^2 \neq 0$.

(b) Parametrize $\mathbf{r}(t) = (\cos t, \sin t, 0)$, $t \in [0, 2\pi]$, $d\mathbf{r} = (-\sin t, \cos t, 0)\,dt$. On the circle $x^2 + y^2 = 1$ so $\mathbf{F} = (-\sin t, \cos t, 0)$. Then $\mathbf{F}\cdot d\mathbf{r} = (\sin^2t + \cos^2t)\,dt = dt$ and $\oint\mathbf{F}\cdot d\mathbf{r} = 2\pi$.

(c) The region ($\mathbb{R}^3$ minus the $z$-axis) is not simply connected: the unit circle loops around the missing axis and cannot be shrunk to a point without crossing it. So (b) $\Rightarrow$ (a) does not apply, and indeed (a) fails by part (b).

On $x > 0$ take $V = -\arctan(y/x)$. Check: $\partial_xV = -\dfrac{1}{1 + y^2/x^2}\cdot\left(-\dfrac{y}{x^2}\right) = \dfrac{y}{x^2 + y^2}$, so $-\partial_xV = -y/(x^2 + y^2) = F_x$ ✓; $\partial_yV = -\dfrac{1}{1 + y^2/x^2}\cdot\dfrac{1}{x} = -\dfrac{x}{x^2 + y^2}$, so $-\partial_yV = F_y$ ✓. In polar coordinates $V = -\varphi$, minus the polar angle. Going once around the origin, $\varphi$ increases by $2\pi$, so any $V$ with $\mathbf{F} = -\nabla V$ on the whole punctured plane would have to satisfy $V(\text{end}) - V(\text{start}) = -\oint\mathbf{F}\cdot d\mathbf{r} = -2\pi$ for a loop that returns to its starting point, i.e. $V(\mathbf{r}) = V(\mathbf{r}) - 2\pi$. Impossible. $V$ exists locally everywhere, but not globally; it is multivalued, exactly like the polar angle.

**B3.** $V = -\tfrac12x^2 + \tfrac14x^4$, $m = 1$.

(a) $V'(x) = -x + x^3 = x(x^2 - 1) = 0$ at $x = 0, \pm 1$. $V''(x) = -1 + 3x^2$: $V''(0) = -1 < 0$ (unstable, a maximum, $V(0) = 0$); $V''(\pm 1) = 2 > 0$ (stable, minima, $V(\pm 1) = -\tfrac14$).

(b) By (3.8), $\omega^2 = V''(\pm 1)/m = 2$, $\omega = \sqrt{2}$ about each well.

(c) The separatrix passes through the unstable equilibrium $(0, 0)$ in phase space, so its energy is $E = V(0) = 0$. Phase portrait: two centers at $(\pm 1, 0)$, each surrounded by closed loops (oscillation in one well, $-\tfrac14 < E < 0$); the separatrix $\tfrac12\dot{x}^2 = \tfrac12x^2 - \tfrac14x^4$, i.e. $\dot{x} = \pm x\sqrt{1 - x^2/2}$, is a figure-eight through the origin with lobes reaching $x = \pm\sqrt{2}$; outside it ($E > 0$) are large closed loops that visit both wells, with $\dot{x} \neq 0$ at $x = 0$. All loops run clockwise ($\dot{x} > 0$ in the upper half-plane). Sketch:

```
   x-dot
     |      ____________
     |    /   ______     \        E > 0
     |   /  /  ___  \     \
     |  |  |  (   )  |     |      inner loops: E < 0
     |   \  \_/   \_/     /       figure-eight: E = 0
     |    \______________/
  ---+---------+----+----+------ x
   -sqrt2     -1    0    1    sqrt2
```

(d) $E = V(2) = -2 + 4 = 2$. At $x = 0$: $\tfrac12v^2 + 0 = 2$, $v = 2$ m/s. Turning point: $V(x) = 2 \Rightarrow x^4 - 2x^2 - 8 = 0 \Rightarrow x^2 = \dfrac{2 \pm\sqrt{4 + 32}}{2} = 4$ (the other root is negative), so $x = -2$. Symmetric, as expected from $V(-x) = V(x)$.

**B4.** (a) $E = \alpha A^4$ (release from rest at $A$), turning points $\pm A$. From (3.7), using symmetry,
$$T = \sqrt{2m}\int_{-A}^{A}\frac{dx}{\sqrt{\alpha(A^4 - x^4)}} = 2\sqrt{2m}\int_0^A\frac{dx}{\sqrt{\alpha(A^4 - x^4)}}.$$
Substitute $x = Au$, $dx = A\,du$, $\sqrt{A^4 - A^4u^4} = A^2\sqrt{1 - u^4}$:
$$T = \frac{2\sqrt{2m}\,A}{\sqrt{\alpha}\,A^2}\int_0^1\frac{du}{\sqrt{1 - u^4}} = \frac{2\sqrt{2}}{A}\sqrt{\frac{m}{\alpha}}\,I = \frac{4}{A}\sqrt{\frac{m}{2\alpha}}\,I,$$
since $2\sqrt2 = 4/\sqrt2$. ∎ The same substitution in $V = \alpha|x|^n$ gives $T \propto A\cdot A^{-n/2} = A^{1 - n/2}$: constant for $n = 2$, $1/A$ for $n = 4$.

(b) $I = \dfrac{\Gamma(1/4)\Gamma(1/2)}{4\Gamma(3/4)} = \dfrac{3.6256\times 1.7725}{4\times 1.2254} = 1.31103$ (the project confirms by quadrature). Why faster at larger amplitude: $x^4$ is flatter than a parabola near the origin and much steeper far out, so a large-amplitude particle spends its time where the restoring force per unit displacement ($V'/x = \alpha x^2$) is large, and it is whipped back faster than a harmonic oscillator would be. The pendulum's $V = mg\ell(1 - \cos\theta)$ is the opposite: *flatter* than a parabola at large angles (the restoring force $mg\sin\theta$ saturates), so large swings are slower.

(c) Reference code (goes under "Your turn" in `project.py`):

```python
alpha, m = 1.0, 1.0
I_quartic = quad(lambda u: 1/np.sqrt(1 - u**4), 0, 1)[0]      # 1.31103
for A in (0.5, 1.0, 2.0, 4.0):
    sol = solve_ivp(lambda t, y: [y[1], -4*alpha*y[0]**3/m], [0, 50/A], [A, 0],
                    events=lambda t, y: y[0], rtol=1e-10, atol=1e-12)
    te = sol.t_events[0]; T_num = 2*(te[1] - te[0])
    T_th = 4/A*np.sqrt(m/(2*alpha))*I_quartic
    assert abs(T_num/T_th - 1) < 1e-6, (A, T_num, T_th)
```
Doubling $A$ halves $T$: $T(1) = 3.708$ s, $T(2) = 1.854$ s, $T(4) = 0.927$ s (with $m = \alpha = 1$).

**B5.** (a) $\omega A = \dfrac{(F_0/m)\,\omega}{\sqrt{(\omega_0^2 - \omega^2)^2 + \gamma^2\omega^2}} = \dfrac{F_0/m}{\sqrt{(\omega_0^2 - \omega^2)^2/\omega^2 + \gamma^2}} = \dfrac{F_0/m}{\sqrt{(\omega_0^2/\omega - \omega)^2 + \gamma^2}}$. The denominator is smallest when $\omega_0^2/\omega = \omega$, i.e. $\omega = \omega_0$, whatever $\gamma$ is. ∎

(b) $F = F_0\cos\omega t$ and $\dot{x} = -\omega A\sin(\omega t - \delta) = -\omega A(\sin\omega t\cos\delta - \cos\omega t\sin\delta)$. So $F\dot{x} = -F_0\omega A\left[\cos\omega t\sin\omega t\cos\delta - \cos^2\omega t\sin\delta\right]$. Over a period $\langle\cos\omega t\sin\omega t\rangle = 0$ and $\langle\cos^2\omega t\rangle = \tfrac12$, giving $\langle F\dot{x}\rangle = \tfrac12F_0\omega A\sin\delta$. Dissipation: $\langle b\dot{x}^2\rangle = b\omega^2A^2\langle\sin^2(\omega t - \delta)\rangle = \tfrac12b\omega^2A^2$. They are equal iff $F_0\sin\delta = b\omega A$. Using $A = (F_0/m)/\sqrt{D}$ with $D = (\omega_0^2 - \omega^2)^2 + \gamma^2\omega^2$: $b\omega A = (b/m)\omega F_0/\sqrt{D} = F_0\gamma\omega/\sqrt{D} = F_0\sin\delta$. ✓ In steady state the drive supplies exactly what the damper eats, as it must, since $\langle dE/dt\rangle = 0$.

(c) $\langle P\rangle = \tfrac12b\omega^2A^2 = \dfrac{F_0^2\gamma}{2m}\cdot\dfrac{\omega^2}{(\omega_0^2 - \omega^2)^2 + \gamma^2\omega^2}$ (using $b/m^2 = \gamma/m$). At $\omega = \omega_0$ this is $P_{\max} = F_0^2/(2m\gamma)$. Half maximum: $\dfrac{\omega^2}{(\omega_0^2 - \omega^2)^2 + \gamma^2\omega^2} = \dfrac{1}{2\gamma^2}$, i.e. $(\omega_0^2 - \omega^2)^2 = \gamma^2\omega^2$, i.e. $\omega_0^2 - \omega^2 = \pm\gamma\omega$. These are two quadratics, $\omega^2 \pm\gamma\omega - \omega_0^2 = 0$, with positive roots $\omega_\pm = \dfrac{\mp\gamma + \sqrt{\gamma^2 + 4\omega_0^2}}{2}$. Their difference is $\omega_+ - \omega_- = \gamma$ **exactly**, for any damping; the weak-damping approximation is not even needed for the width (it is needed to say the half-power points are *symmetric* about $\omega_0$: $\omega_\pm \approx \omega_0 \pm\gamma/2$ when $\gamma \ll \omega_0$, from $\sqrt{\gamma^2 + 4\omega_0^2} \approx 2\omega_0$). Alternatively, with the approximation $\omega_0^2 - \omega^2 \approx 2\omega_0(\omega_0 - \omega)$ near the peak, the half-power condition becomes $2\omega_0|\omega_0 - \omega| = \gamma\omega_0$, so $|\omega_0 - \omega| = \gamma/2$. The project's `brentq` on Example 4.4 finds the half-power frequencies at $9.0499$ and $11.0499$ rad/s: width $2.0000 = \gamma$. ✓

**B6.** (a) $m(t) = m_0 - \beta t$ and $\dot{m} = -\beta$. From (3.11) with gravity: $m\dot{v} = u\beta - mg$, so $\dot{v} = \dfrac{u\beta}{m_0 - \beta t} - g$. Integrate from $0$ with $v(0) = 0$: $v = -u\ln(m_0 - \beta t) + u\ln m_0 - gt = u\ln\dfrac{m_0}{m_0 - \beta t} - gt$. ∎

(b) At $t = 0$, $\dot{v} = u\beta/m_0 - g$, positive only if $u\beta > m_0g$: thrust must exceed initial weight. If not, the rocket sits on the pad (the pad supplies the normal force that makes the net force zero) and burns fuel; $m(t)$ decreases until $u\beta = m(t)g$, and only then does it lift off, having wasted the fuel burned so far. The formula in (a) does not apply during the sitting phase.

(c) $\beta = 700/60 = 11.667$ kg/s, so $m_f = m_0 - 60\beta = 300$ kg. $v(60) = 2500\ln\dfrac{1000}{300} - 9.81\times 60 = 3009.9 - 588.6 = 2421.3$ m/s, matching Example 4.6.

Height. $h = \int_0^{t_b}v\,dt = u\int_0^{t_b}\ln\dfrac{m_0}{m_0 - \beta t}\,dt - \tfrac12gt_b^2$. Substitute $s = m_0 - \beta t$, $dt = -ds/\beta$: $\int_0^{t_b}\ln\dfrac{m_0}{m_0 - \beta t}\,dt = \dfrac{1}{\beta}\int_{m_f}^{m_0}\ln\dfrac{m_0}{s}\,ds$. Since $\dfrac{d}{ds}\left[s\ln\dfrac{m_0}{s} + s\right] = \ln\dfrac{m_0}{s} - 1 + 1 = \ln\dfrac{m_0}{s}$, the integral is $\dfrac{1}{\beta}\left[(m_0 - m_f) - m_f\ln\dfrac{m_0}{m_f}\right] = t_b - \dfrac{m_f}{\beta}\ln\dfrac{m_0}{m_f}$. So
$$h = ut_b - \frac{um_f}{\beta}\ln\frac{m_0}{m_f} - \tfrac12gt_b^2 = 150{,}000 - 64{,}285.7\times 1.20397 - 17{,}658 = 150{,}000 - 77{,}398.3 - 17{,}658 = 54{,}944\ \text{m},$$
about 54.9 km at burnout. The project integrates $v(t)$ numerically and gets the same value.

**B7.** (a) Over the half-cycle from rest at $x = A_1$ to rest at $x = -A_2$, $\Delta T = 0$. The work–energy theorem (3.4) with the net force: spring work $= -\Delta V = \tfrac12kA_1^2 - \tfrac12kA_2^2$; friction work $= -\mu mg\times(\text{distance}) = -\mu mg(A_1 + A_2)$ (friction opposes motion throughout, since the velocity does not change sign within a half-cycle). So $\tfrac12k(A_1^2 - A_2^2) = \mu mg(A_1 + A_2)$. Divide by $(A_1 + A_2) > 0$: $\tfrac12k(A_1 - A_2) = \mu mg$, i.e. $A_1 - A_2 = 2\mu mg/k$. Two half-cycles per period: the amplitude drops by $4\mu mg/k$ each period. The decrement is independent of the amplitude, so the decay is linear; the amplitude reaches zero in finite time. Linear damping gives $A_{n+1}/A_n = e^{-\gamma T/2}$, a constant *ratio*, hence exponential decay that never reaches zero. ∎

(b) Reference implementation (fixed-step RK4 with an explicit sticking rule; an adaptive integrator chatters on the sign flip when the mass stops):

```python
m, k, mu, x0 = 1.0, 100.0, 0.1, 0.5
def a_coulomb(x, v):
    return (-k*x - mu*m*g*math.copysign(1.0, v))/m
def rk4_xv(a, x, v, h):
    k1x, k1v = v, a(x, v)
    k2x, k2v = v + h/2*k1v, a(x + h/2*k1x, v + h/2*k1v)
    k3x, k3v = v + h/2*k2v, a(x + h/2*k2x, v + h/2*k2v)
    k4x, k4v = v + h*k3v, a(x + h*k3x, v + h*k3v)
    return x + h/6*(k1x+2*k2x+2*k3x+k4x), v + h/6*(k1v+2*k2v+2*k3v+k4v)
h, t, x, v = 1e-4, 0.0, x0, -1e-12          # nudge so sgn(v) is defined
amps, times = [x0], [0.0]
while t < 20:
    xn, vn = rk4_xv(a_coulomb, x, v, h); t += h
    if vn*v < 0:                                 # turning point
        amps.append(abs(xn)); times.append(t)
        if k*abs(xn) <= mu*m*g: x = xn; break    # spring can't beat static friction: stuck
    x, v = xn, vn
amps = np.array(amps)
assert np.allclose(np.diff(amps)[:-1], -2*mu*m*g/k, rtol=1e-2)   # 0.01962 per half-cycle
```
Output: turning-point amplitudes $0.5, 0.48038, 0.46076, 0.44114, \dots$, decreasing by $0.019620$ per half-cycle, i.e. $4\mu mg/k = 0.03924$ per period. The mass sticks when $k|x| \le \mu mg$, i.e. $|x| \le 0.00981$ m, at a turning point. After $n$ half-cycles $A_n = 0.5 - 0.01962n$; $A_{24} = 0.0291 > 0.00981$ (still moves), $A_{25} = 0.0095 \le 0.00981$: it stops after 25 half-cycles at $x = -0.0095$ m, at $t = 7.854$ s. Why that time: during each half-cycle the force is $-kx \mp\mu mg$, which is simple harmonic motion about the shifted center $x_c = \mp\mu mg/k$ at the undamped frequency $\omega_0 = 10$ rad/s, so every half-cycle lasts exactly $\pi/\omega_0 = 0.31416$ s, and $25\times 0.31416 = 7.854$ s. Coulomb friction changes the amplitude but not the period.

**B8.** (a) Tangent to the wire: $(1, f')$; unit normal $\hat{\mathbf{n}} = (-f', 1)/\sqrt{1 + f'^2}$. Let $\lambda = N/\sqrt{1 + f'^2}$. Newton: $m\ddot{x} = -f'\lambda$, $m\ddot{y} = \lambda - mg$. Constraint $y = f(x)$: $\dot{y} = f'\dot{x}$, $\ddot{y} = f''\dot{x}^2 + f'\ddot{x}$. Substitute into the $y$-equation: $\lambda = mg + m(f''\dot{x}^2 + f'\ddot{x})$. Substitute into the $x$-equation: $m\ddot{x} = -f'\left[mg + mf''\dot{x}^2 + mf'\ddot{x}\right]$, so $\ddot{x}(1 + f'^2) = -f'(g + f''\dot{x}^2)$:
$$\ddot{x} = -\frac{f'(g + f''\dot{x}^2)}{1 + f'^2}, \qquad N = \lambda\sqrt{1 + f'^2} = m\left(g + f''\dot{x}^2 + f'\ddot{x}\right)\sqrt{1 + f'^2}.$$ ∎

(b) $f = kx^2$: $f' = 2kx$, $f'' = 2k$. $\ddot{x} = -\dfrac{2kx(g + 2k\dot{x}^2)}{1 + 4k^2x^2}$, which is (3.14). $N = m(g + 2k\dot{x}^2 + 2kx\ddot{x})\sqrt{1 + 4k^2x^2}$, which is step 6 of §3.10(b); at $x = 0$, $\ddot{x} = 0$ and $N = m(g + 2k\dot{x}^2)$, which is (3.15). ✓

(c) $E = \tfrac12m(1 + f'^2)\dot{x}^2 + mgf$. Differentiate, using $\frac{d}{dt}f'^2 = 2f'f''\dot{x}$ and $\frac{d}{dt}f = f'\dot{x}$:
$$\dot{E} = m(1 + f'^2)\dot{x}\ddot{x} + mf'f''\dot{x}^3 + mgf'\dot{x} = m\dot{x}\left[(1 + f'^2)\ddot{x} + f'f''\dot{x}^2 + gf'\right].$$
Insert $(1 + f'^2)\ddot{x} = -f'(g + f''\dot{x}^2) = -gf' - f'f''\dot{x}^2$: the bracket is $-gf' - f'f''\dot{x}^2 + f'f''\dot{x}^2 + gf' = 0$. So $\dot{E} = 0$. ∎ (In CM-3, $\tfrac12m(1 + f'^2)\dot{x}^2$ is the kinetic energy in the coordinate $x$, and energy conservation falls out of the Lagrangian without this computation.)

## C. Challenge

**C1.** Binomial series: $(1 - z)^{-1/2} = 1 + \tfrac12z + \tfrac{3}{8}z^2 + O(z^3)$. With $z = k^2\sin^2\phi$:
$$K(k) = \int_0^{\pi/2}\left[1 + \tfrac12k^2\sin^2\phi + \tfrac38k^4\sin^4\phi + O(k^6)\right]d\phi = \frac{\pi}{2} + \frac{k^2}{2}\cdot\frac{\pi}{4} + \frac{3k^4}{8}\cdot\frac{3\pi}{16} + O(k^6) = \frac{\pi}{2}\left(1 + \frac{k^2}{4} + \frac{9k^4}{64} + O(k^6)\right).$$
Now $k = \sin(\theta_0/2) = \dfrac{\theta_0}{2} - \dfrac{\theta_0^3}{48} + O(\theta_0^5)$, so $k^2 = \dfrac{\theta_0^2}{4} - 2\cdot\dfrac{\theta_0}{2}\cdot\dfrac{\theta_0^3}{48} + O(\theta_0^6) = \dfrac{\theta_0^2}{4} - \dfrac{\theta_0^4}{48} + O(\theta_0^6)$ and $k^4 = \dfrac{\theta_0^4}{16} + O(\theta_0^6)$. Then
$$\frac{T}{T_0} = \frac{2K}{\pi} = 1 + \frac14\left(\frac{\theta_0^2}{4} - \frac{\theta_0^4}{48}\right) + \frac{9}{64}\cdot\frac{\theta_0^4}{16} + O(\theta_0^6) = 1 + \frac{\theta_0^2}{16} + \theta_0^4\left(-\frac{1}{192} + \frac{9}{1024}\right) + O(\theta_0^6).$$
Common denominator 3072: $-\dfrac{16}{3072} + \dfrac{27}{3072} = \dfrac{11}{3072}$. So $T = T_0\left(1 + \dfrac{\theta_0^2}{16} + \dfrac{11\theta_0^4}{3072} + O(\theta_0^6)\right)$. ∎

At $\theta_0 = \pi/2$: $\theta_0^2 = 2.4674$, $\theta_0^4 = 6.0881$. One term: $1 + 0.15421 = 1.15421$. Two terms: $1.15421 + 11\times 6.0881/3072 = 1.15421 + 0.02180 = 1.17601$. Exact: $1.18034$. The relative error drops from 2.2% to 0.37%. SymPy (project, Part 9) reproduces the series from `elliptic_k`.

**C2.** (a) $\Phi(x, v) = (x + hv,\ v + ha(x))$. Jacobian $\begin{pmatrix}\partial_x(x + hv) & \partial_v(x + hv)\\ \partial_x(v + ha) & \partial_v(v + ha)\end{pmatrix} = \begin{pmatrix}1 & h\\ ha'(x) & 1\end{pmatrix}$, determinant $1 - h^2a'(x)$. For $a = -\omega^2x$, $a' = -\omega^2$ and the determinant is $1 + h^2\omega^2$: the area inflation factor of (3.17). ∎

(b) $\Psi(x, v) = (x + hv + h^2a(x),\ v + ha(x))$. Jacobian $\begin{pmatrix}1 + h^2a'(x) & h\\ ha'(x) & 1\end{pmatrix}$, determinant $(1 + h^2a') - h\cdot ha' = 1$. ∎ Exactly area-preserving for any force law, any step.

(c) $v' = v - h\omega^2x$, $x' = x + hv' = (1 - h^2\omega^2)x + hv$. Compute the three pieces of $\tilde{E}(x', v')$:
$v'^2 = v^2 - 2h\omega^2xv + h^2\omega^4x^2$.
$x'^2 = (1 - h^2\omega^2)^2x^2 + 2h(1 - h^2\omega^2)xv + h^2v^2$.
$x'v' = (1 - h^2\omega^2)xv - h\omega^2(1 - h^2\omega^2)x^2 + hv^2 - h^2\omega^2xv = (1 - 2h^2\omega^2)xv - h\omega^2(1 - h^2\omega^2)x^2 + hv^2$.
Now $\tilde{E}' = \tfrac12v'^2 + \tfrac12\omega^2x'^2 - \tfrac12h\omega^2x'v'$. Collect coefficients:
$v^2$: $\tfrac12 + \tfrac12\omega^2h^2 - \tfrac12h\omega^2\cdot h = \tfrac12$.
$xv$: $-h\omega^2 + h\omega^2(1 - h^2\omega^2) - \tfrac12h\omega^2(1 - 2h^2\omega^2) = -h\omega^2 + h\omega^2 - h^3\omega^4 - \tfrac12h\omega^2 + h^3\omega^4 = -\tfrac12h\omega^2$.
$x^2$: $\tfrac12h^2\omega^4 + \tfrac12\omega^2(1 - 2h^2\omega^2 + h^4\omega^4) + \tfrac12h^2\omega^4(1 - h^2\omega^2) = \tfrac12\omega^2 + h^2\omega^4(\tfrac12 - 1 + \tfrac12) + h^4\omega^6(\tfrac12 - \tfrac12) = \tfrac12\omega^2$.
So $\tilde{E}(x', v') = \tfrac12v^2 + \tfrac12\omega^2x^2 - \tfrac12h\omega^2xv = \tilde{E}(x, v)$. ∎

Boundedness: $E - \tilde{E} = \tfrac12h\omega^2xv$, and $|xv| \le \tfrac12(\omega x^2 + v^2/\omega) = E/\omega$, so $|E - \tilde{E}| \le \tfrac12h\omega E$. Since $\tilde{E}$ is exactly constant, $E$ is trapped between $\tilde{E}/(1 + h\omega/2)$ and $\tilde{E}/(1 - h\omega/2)$ for all time (for $h\omega < 2$): a band of relative width $O(h\omega)$ that never grows. Euler, by contrast, multiplies $E$ by $1 + h^2\omega^2$ every step. The project checks the identity symbolically and watches $E$ oscillate numerically.

**C3.** (a) Bob position $(x, y) = (X + \ell\sin\theta,\ -\ell\cos\theta)$. Differentiating twice: $\ddot{x} = \ddot{X} + \ell\cos\theta\,\ddot{\theta} - \ell\sin\theta\,\dot{\theta}^2$, $\ddot{y} = \ell\sin\theta\,\ddot{\theta} + \ell\cos\theta\,\dot{\theta}^2$. Tension points from bob to pivot, along $-(\sin\theta, -\cos\theta)$: $\mathbf{T} = (-T\sin\theta, T\cos\theta)$. Newton:
$$m\left(\ddot{X} + \ell\cos\theta\,\ddot{\theta} - \ell\sin\theta\,\dot{\theta}^2\right) = -T\sin\theta, \qquad m\left(\ell\sin\theta\,\ddot{\theta} + \ell\cos\theta\,\dot{\theta}^2\right) = T\cos\theta - mg.$$

(b) Multiply the first by $\cos\theta$, the second by $\sin\theta$, add. The $T$ terms cancel, the $\dot{\theta}^2$ terms cancel ($-\sin\theta\cos\theta + \cos\theta\sin\theta$), and $\ddot{\theta}$ picks up $\cos^2\theta + \sin^2\theta = 1$:
$$m\left(\ddot{X}\cos\theta + \ell\ddot{\theta}\right) = -mg\sin\theta \quad\Longrightarrow\quad \ell\ddot{\theta} = -g\sin\theta - \ddot{X}\cos\theta = -g\sin\theta + A\Omega^2\cos(\Omega t)\cos\theta,$$
using $\ddot{X} = -A\Omega^2\cos\Omega t$. ∎ For $T$: multiply the first by $-\sin\theta$, the second by $\cos\theta$, add: $m\left(-\ddot{X}\sin\theta + \ell\dot{\theta}^2\right) = T - mg\cos\theta$, so
$$T = m\ell\dot{\theta}^2 + mg\cos\theta - m\ddot{X}\sin\theta = m\ell\dot{\theta}^2 + mg\cos\theta + mA\Omega^2\cos(\Omega t)\sin\theta.$$
At $A = 0$ this is the fixed-pivot result of §3.10(a) step 7.

(c) Small $\theta$ ($\sin\theta \approx\theta$, $\cos\theta \approx 1$): $\ddot{\theta} + \dfrac{g}{\ell}\theta = \dfrac{A\Omega^2}{\ell}\cos\Omega t$. This is §3.8 with $\omega_0^2 = g/\ell$, $\gamma = 0$, and $F_0/m \to A\Omega^2/\ell$. From (3.10), the steady-state amplitude is $\Theta = \dfrac{A\Omega^2/\ell}{|\omega_0^2 - \Omega^2|}$, with phase $0$ below resonance and $\pi$ above. As $\Omega\to\omega_0$ this diverges. Two things the model omits prevent that: damping (air drag, pivot friction) caps the amplitude at $Q$ times the static value, and the small-angle approximation fails first anyway: at large $\theta$ the pendulum's own frequency drops (Example 4.2: the period grows with amplitude), so the pendulum detunes itself from the drive. The real response is finite and, for strong driving, not even periodic.

(d) CM-3 writes $L = T - V = \tfrac12m\left[(\dot{X} + \ell\cos\theta\,\dot{\theta})^2 + \ell^2\sin^2\theta\,\dot{\theta}^2\right] + mg\ell\cos\theta$ and takes $\frac{d}{dt}\frac{\partial L}{\partial\dot{\theta}} = \frac{\partial L}{\partial\theta}$; the tension never appears, and (b) is the result.
