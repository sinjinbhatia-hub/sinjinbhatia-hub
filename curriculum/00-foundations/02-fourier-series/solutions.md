# F2 Solutions

Conventions as in the notes. Every number below is printed or asserted in `project.py`: Part 9 ("Exercise checks"), plus the B8 reference block just before it.

## A. Warm-up

**A1.** The product-to-sum formulas (§3.2 step 2): $2\cos A\cos B = \cos(A-B) + \cos(A+B)$, $2\sin A\sin B = \cos(A-B) - \cos(A+B)$, $2\sin A\cos B = \sin(A-B) + \sin(A+B)$.

(a) $\cos2x\cos5x = \tfrac12\big(\cos3x + \cos7x\big)$. Both have nonzero integer frequency, so $\int_{-\pi}^\pi\cos2x\cos5x\,dx = 0$.

(b) $\sin^23x = \tfrac12\big(\cos0 - \cos6x\big) = \tfrac12 - \tfrac12\cos6x$. The constant integrates to $\tfrac12\cdot2\pi = \pi$, the cosine to $0$. So $\int_{-\pi}^{\pi}\sin^23x\,dx = \pi$.

(c) $\cos4x\sin4x = \tfrac12\sin8x$, and $\int_{-\pi}^\pi\sin8x\,dx = \big[-\tfrac{\cos8x}{8}\big]_{-\pi}^{\pi} = 0$ (cosine is even). So the integral is $0$.

(d) $e^{-2ix}e^{3ix} = e^{ix}$, and $\int_{-\pi}^{\pi}e^{ix}dx = \big[\tfrac{e^{ix}}{i}\big]_{-\pi}^{\pi} = \tfrac{e^{i\pi} - e^{-i\pi}}{i} = \tfrac{(-1) - (-1)}{i} = 0$.

(e) $\langle e^{3ix}, e^{3ix}\rangle = \int_{-\pi}^{\pi}\overline{e^{3ix}}e^{3ix}\,dx = \int_{-\pi}^\pi1\,dx = 2\pi$.

Norms, from Theorems 3.2 and 3.3: $\|1\|^2 = 2\pi$, $\|\cos7x\|^2 = \pi$, $\|\sin7x\|^2 = \pi$, $\|e^{7ix}\|^2 = 2\pi$.

**A2.** (a) Match $f$ against $\frac{a_0}{2} + \sum(a_n\cos nx + b_n\sin nx)$: the constant term is $2 = a_0/2$, so $a_0 = 4$; $a_1 = 3$, $b_1 = -1$; $a_2 = 0$, $b_2 = 4$; $a_3 = -5$, $b_3 = 0$; everything else $0$. By (3.7), $c_0 = a_0/2 = 2$, $c_{\pm1} = \frac{a_1\mp ib_1}{2} = \frac{3\pm i}{2}$ (so $c_1 = \frac{3+i}{2}$, $c_{-1} = \frac{3-i}{2}$), $c_{\pm2} = \frac{0\mp4i}{2} = \mp2i$, $c_{\pm3} = -\frac52$.

Why no integral: $f\in\mathcal{T}_3$, and Theorem 3.5 says $S_3f$ is the *unique* closest point of $\mathcal{T}_3$ to $f$; since $f$ itself is at distance $0$, $S_3f = f$, and the coefficients of a vector in an orthogonal basis are unique. The integrals (3.4) would reproduce the numbers you can read off.

(b) Real form, step 1 of the Bessel derivation: $\|f\|^2 = \pi\Big(\frac{a_0^2}{2} + \sum_n(a_n^2 + b_n^2)\Big) = \pi\big(\tfrac{16}{2} + 9 + 1 + 16 + 25\big) = 59\pi\approx185.354$. Complex form: $\|f\|^2 = 2\pi\sum|c_n|^2 = 2\pi\big(4 + 2\cdot\tfrac{10}{4} + 2\cdot4 + 2\cdot\tfrac{25}{4}\big) = 2\pi\big(4 + 5 + 8 + 12.5\big) = 59\pi$. ✓

(c) Theorem 3.5: the closest element of $\mathcal{T}_1$ is $S_1f = 2 + 3\cos x - \sin x$, and by Pythagoras $\|f - S_1f\|^2 = \pi\big(b_2^2 + a_3^2\big) = \pi(16 + 25) = 41\pi\approx128.805$. For degree $\le2$: $S_2f = 2 + 3\cos x - \sin x + 4\sin2x$ and $\|f - S_2f\|^2 = \pi a_3^2 = 25\pi\approx78.540$.

**A3.** (a) Invert (3.7): $a_n = c_n + c_{-n}$, $b_n = i(c_n - c_{-n})$. $a_0 = 2c_0 = 6$; $a_1 = (1-2i) + (1+2i) = 2$, $b_1 = i\big((1-2i) - (1+2i)\big) = i(-4i) = 4$; $a_2 = \tfrac i2 - \tfrac i2 = 0$, $b_2 = i\big(\tfrac i2 + \tfrac i2\big) = i\cdot i = -1$. So $f(x) = 3 + 2\cos x + 4\sin x - \sin2x$. Check that it is real: $c_{-n} = \overline{c_n}$ holds for every $n$. ✓

(b) $\cos^2x = \tfrac12 + \tfrac12\cos2x$ (from $\cos2x = 2\cos^2x - 1$). Real: $a_0 = 1$ (constant term $\tfrac12$), $a_2 = \tfrac12$, all else $0$. Complex: $c_0 = \tfrac12$, $c_{\pm2} = \tfrac14$. Directly: $\cos^2x = \tfrac14(e^{ix} + e^{-ix})^2 = \tfrac14(e^{2ix} + 2 + e^{-2ix})$. ✓

$\sin^3x$: from $\sin3x = 3\sin x - 4\sin^3x$, $\sin^3x = \tfrac34\sin x - \tfrac14\sin3x$. Real: $b_1 = \tfrac34$, $b_3 = -\tfrac14$. Complex, by (3.7): $c_{\pm1} = \mp\tfrac{i}{2}b_1 = \mp\tfrac{3i}{8}$, $c_{\pm3} = \mp\tfrac i2b_3 = \pm\tfrac i8$. Directly: $\sin^3x = \Big(\frac{e^{ix} - e^{-ix}}{2i}\Big)^3 = \frac{e^{3ix} - 3e^{ix} + 3e^{-ix} - e^{-3ix}}{-8i} = \frac{i}{8}\big(e^{3ix} - 3e^{ix} + 3e^{-ix} - e^{-3ix}\big)$, so $c_3 = \tfrac i8$, $c_1 = -\tfrac{3i}{8}$, $c_{-1} = \tfrac{3i}{8}$, $c_{-3} = -\tfrac i8$. ✓

(c) $a_2 = \frac1\pi\int_{-\pi}^\pi\cos^2x\cos2x\,dx = \frac1\pi\int_{-\pi}^{\pi}\big(\tfrac12 + \tfrac12\cos2x\big)\cos2x\,dx = \frac1\pi\big(0 + \tfrac12\cdot\pi\big) = \tfrac12$, using $\int\cos2x = 0$ and $\int\cos^22x = \pi$ (Theorem 3.2). ✓

**A4.** By §3.5 with $\|\sin(n\pi x/L)\|^2 = L/2$ on $[0, L]$:
$$b_n = \frac2L\int_0^L\sin\frac{n\pi x}{L}\,dx = \frac2L\cdot\frac{L}{n\pi}\Big[-\cos\frac{n\pi x}{L}\Big]_0^L = \frac{2}{n\pi}\big(1 - \cos n\pi\big) = \frac{2\big(1 - (-1)^n\big)}{n\pi} = \begin{cases}\dfrac{4}{n\pi} & n\text{ odd}\\ 0 & n\text{ even.}\end{cases}$$
So $1\sim\frac4\pi\sum_{n\text{ odd}}\frac1n\sin\frac{n\pi x}{L}$ on $(0, L)$.

At $x = L/2$: $\sin(n\pi/2) = (-1)^{(n-1)/2}$ for odd $n$, so the series is $\frac4\pi\big(1 - \tfrac13 + \tfrac15 - \cdots\big) = \frac4\pi\cdot\frac\pi4 = 1$ (Leibniz). ✓ At $x = 0$ every term is $0$, so the series gives $0$, not $1$.

The sine series is the Fourier series of the **odd extension** of $f$ to $(-L, L)$: $g(x) = \operatorname{sgn}(x)$. That is the square wave of Example 4.1 with $x$ replaced by $\pi x/L$ (compare: $\frac4\pi\sum_{\text{odd}}\frac{\sin nx}{n}$). The value $0$ at $x = 0$ is Dirichlet's theorem giving the midpoint of the jump from $-1$ to $+1$.

$S_5f(L/2) = \frac4\pi\big(1 - \tfrac13 + \tfrac15\big) = \frac4\pi\cdot\frac{13}{15} = \frac{52}{15\pi}\approx1.1035$. (Above $1$; the alternating series overshoots and undershoots its limit in turn.)

**A5.** Direct, by parts with $u = t$, $dv = e^{-2\pi int}dt$, $v = \frac{e^{-2\pi int}}{-2\pi in}$, for $n\ne0$:
$$\hat c_n = \Big[t\,\frac{e^{-2\pi int}}{-2\pi in}\Big]_{-1/2}^{1/2} + \frac{1}{2\pi in}\int_{-1/2}^{1/2}e^{-2\pi int}dt.$$
The last integral is $0$ (a full period of a nonconstant exponential). The boundary term: $e^{\mp\pi in} = (-1)^n$, so it is $\frac{(-1)^n}{-2\pi in}\Big(\tfrac12 - \big(-\tfrac12\big)\Big) = -\frac{(-1)^n}{2\pi in} = \frac{i(-1)^n}{2\pi n}$, using $-\frac1i = i$. And $\hat c_0 = \int_{-1/2}^{1/2}t\,dt = 0$. So $\hat c_n = \dfrac{i(-1)^n}{2\pi n}$ for $n\ne0$; for instance $\hat c_1 = -\dfrac{i}{2\pi}\approx-0.1592\,i$.

By substitution: (3.9) with $T = 1$ is $\hat c_n = \int_{-1/2}^{1/2}f(t)e^{-2\pi int}dt$. Put $x = 2\pi t$, so $t = x/2\pi$, $dt = dx/2\pi$, and $f(t) = t = x/2\pi$:
$$\hat c_n = \int_{-\pi}^{\pi}\frac{x}{2\pi}\,e^{-inx}\,\frac{dx}{2\pi} = \frac{1}{2\pi}\cdot\underbrace{\frac{1}{2\pi}\int_{-\pi}^{\pi}x\,e^{-inx}dx}_{c_n\text{ of Example 4.2}} = \frac{1}{2\pi}\cdot\frac{i(-1)^n}{n}.$$
Same answer. The lesson: Osgood's sawtooth is our sawtooth scaled by $1/2\pi$ in height, so all coefficients scale by $1/2\pi$; the change of period itself costs nothing because (3.9) carries the $1/T$.

Real form: $\hat c_{-n} = -\hat c_n$, so $a_n = \hat c_n + \hat c_{-n} = 0$ and $b_n = i(\hat c_n - \hat c_{-n}) = 2i\hat c_n = 2i\cdot\frac{i(-1)^n}{2\pi n} = \frac{(-1)^{n+1}}{\pi n}$. Hence $t\sim\frac1\pi\sum_{n\ge1}\frac{(-1)^{n+1}}{n}\sin2\pi nt$. At $t = \tfrac14$: $\sin(n\pi/2)$ is $0$ for even $n$ and $(-1)^{(n-1)/2}$ for odd $n$, so the series is $\frac1\pi\big(1 - \tfrac13 + \tfrac15 - \cdots\big) = \frac1\pi\cdot\frac\pi4 = \tfrac14$. ✓

## B. Core

**B1.** (a) $b_n = \frac2\pi\int_0^\pi x(\pi - x)\sin nx\,dx = \frac2\pi\big(\pi I_1 - I_2\big)$ with $I_1 = \int_0^\pi x\sin nx\,dx$ and $I_2 = \int_0^\pi x^2\sin nx\,dx$.

From Example 4.2, $I_1 = -\dfrac{\pi(-1)^n}{n}$. For $I_2$, by parts with $u = x^2$, $dv = \sin nx\,dx$:
$$I_2 = \Big[-\frac{x^2\cos nx}{n}\Big]_0^\pi + \frac2n\int_0^\pi x\cos nx\,dx = -\frac{\pi^2(-1)^n}{n} + \frac2n\cdot\Big(-\frac{1 - (-1)^n}{n^2}\Big),$$
using $\int_0^\pi x\cos nx\,dx = -\frac{1 - (-1)^n}{n^2}$ from Example 4.3. So
$$\pi I_1 - I_2 = -\frac{\pi^2(-1)^n}{n} + \frac{\pi^2(-1)^n}{n} + \frac{2\big(1 - (-1)^n\big)}{n^3} = \frac{2\big(1 - (-1)^n\big)}{n^3},$$
$$b_n = \frac2\pi\cdot\frac{2\big(1 - (-1)^n\big)}{n^3} = \frac{4\big(1 - (-1)^n\big)}{\pi n^3} = \begin{cases}\dfrac{8}{\pi n^3} & n\text{ odd}\\ 0 & n\text{ even.}\end{cases}$$
$b_1 = \frac8\pi\approx2.5465$, $b_3 = \frac{8}{27\pi}\approx0.0943$, $b_5 = \frac{8}{125\pi}\approx0.0204$. So $x(\pi - x)\sim\frac8\pi\Big(\sin x + \frac{\sin3x}{27} + \frac{\sin5x}{125} + \cdots\Big)$ on $[0,\pi]$.

(b) At $x = \pi/2$: $f(\pi/2) = \frac\pi2\cdot\frac\pi2 = \frac{\pi^2}{4}$ and $\sin(n\pi/2) = (-1)^{(n-1)/2}$ for odd $n$. The odd extension is continuous and piecewise smooth, so Theorem 3.14 gives equality:
$$\frac{\pi^2}{4} = \frac8\pi\sum_{n\text{ odd}}\frac{(-1)^{(n-1)/2}}{n^3} \quad\Longrightarrow\quad 1 - \frac1{27} + \frac1{125} - \frac1{343} + \cdots = \frac{\pi^3}{32}\approx0.96895.$$

(c) The odd extension is $g(x) = x(\pi - |x|)$ on $[-\pi,\pi]$ (for $x < 0$: $-f(-x) = -(-x)(\pi + x) = x(\pi + x)$). Sketch: a smooth hump on $(0,\pi)$ of height $\pi^2/4$, its mirror image below the axis on $(-\pi, 0)$, repeated. $g(\pm\pi) = 0$, so $g$ is continuous across the seam. $g'(x) = \pi - 2|x|$, continuous everywhere, and $g'(\pm\pi) = -\pi$ on both sides, so $g'$ is continuous across the seam. $g''(x) = -2\operatorname{sgn}x$ jumps: at $x = 0$ from $+2$ to $-2$ ($J = -4$) and at the seam from $-2$ (at $\pi^-$) to $+2$ (at $-\pi^+$) ($J = +4$). So $g, g'$ are continuous and periodic and $g''$ is piecewise smooth: Theorem 3.10 with $k = 1$ gives $|c_n|\le C/n^3$. The $1/n^3$ is the price of the jump in the *second* derivative.

Recovering $b_n$ without integrals: $g'' = -2\operatorname{sgn}x$ is $-2$ times the square wave, so by Example 4.1 $c_n(g'') = -2\cdot\big(-\frac{2i}{\pi n}\big) = \frac{4i}{\pi n}$ for odd $n$, $0$ for even $n\ne0$. Theorem 3.18 twice ($g$ and $g'$ are continuous periodic): $c_n(g'') = (in)^2c_n(g) = -n^2c_n(g)$, so $c_n(g) = -\frac{4i}{\pi n^3}$ for odd $n$. Then $b_n = i(c_n - c_{-n}) = i\cdot2c_n = i\cdot2\cdot\Big(-\frac{4i}{\pi n^3}\Big) = \frac{8}{\pi n^3}$ for odd $n$. ✓ (Check also via Theorem 3.9 on $g''$: jumps $-4$ at $0$ and $+4$ at $\pi$, $g''' = 0$ on the pieces, so $c_n(g'') = \frac{1}{in}\cdot\frac{1}{2\pi}\big(-4 + 4(-1)^n\big) = \frac{-4(1 - (-1)^n)}{2\pi in} = \frac{2i(1 - (-1)^n)}{\pi n}$, which is $\frac{4i}{\pi n}$ for odd $n$. ✓)

(d) $f(\pi - x) = (\pi - x)x = f(x)$: $f$ is symmetric about $\pi/2$. And $\sin\big(n(\pi - x)\big) = \sin n\pi\cos nx - \cos n\pi\sin nx = (-1)^{n+1}\sin nx$: symmetric about $\pi/2$ for odd $n$, antisymmetric for even $n$. For even $n$ the integrand $f(x)\sin nx$ is antisymmetric about $\pi/2$ and its integral over $[0,\pi]$ vanishes. Same mechanism as Example 4.7.

**B2.** (a) By (3.8), $\langle e^{imx}, f\rangle = 2\pi c_m$. By linearity in the second slot and Theorem 3.3, $\langle e^{imx}, S_Nf\rangle = \sum_{|n|\le N}c_n\langle e^{imx}, e^{inx}\rangle = \sum_n c_n\cdot2\pi\delta_{mn} = 2\pi c_m$ for $|m|\le N$. Subtract: $\langle e^{imx}, f - S_Nf\rangle = 0$. Then, by conjugate-linearity in the first slot, $\langle S_Nf, f - S_Nf\rangle = \sum_{|m|\le N}\overline{c_m}\,\langle e^{imx}, f - S_Nf\rangle = 0$.

(b) $\|S_Nf\|^2 = \Big\langle\sum_m c_me^{imx}, \sum_nc_ne^{inx}\Big\rangle = \sum_{m,n}\overline{c_m}c_n\langle e^{imx}, e^{inx}\rangle = \sum_{m,n}\overline{c_m}c_n\cdot2\pi\delta_{mn} = 2\pi\sum_{|n|\le N}|c_n|^2$.

(c) Write $f = S_Nf + (f - S_Nf)$. The two pieces are orthogonal by (a), so Pythagoras gives $\|f\|^2 = \|S_Nf\|^2 + \|f - S_Nf\|^2$ (the cross term $2\operatorname{Re}\langle S_Nf, f - S_Nf\rangle$ is $0$). Hence
$$\|f - S_Nf\|^2 = \|f\|^2 - 2\pi\sum_{|n|\le N}|c_n|^2\ \ge 0, \qquad\text{i.e.}\qquad \sum_{|n|\le N}|c_n|^2\le\frac{1}{2\pi}\|f\|^2. \ \blacksquare$$

(d) Square wave: $\frac{1}{2\pi}\|f\|^2 = \frac{1}{2\pi}\int_{-\pi}^\pi1\,dx = 1$. Its coefficients are $c_n = -\frac{2i}{\pi n}$ for odd $n$ (both signs), so $\sum_{|n|\le N}|c_n|^2 = 2\sum_{\text{odd }k\le N}\frac{4}{\pi^2k^2} = \frac{8}{\pi^2}\sum_{\text{odd }k\le N}\frac1{k^2}$. The relative mean-square error is
$$\frac{1}{2\pi}\|f - S_Nf\|^2 = 1 - \frac{8}{\pi^2}\sum_{\text{odd }k\le N}\frac{1}{k^2}:$$
$N = 1$: $1 - \frac{8}{\pi^2} = 0.1894$; $N = 5$: $1 - \frac{8}{\pi^2}\big(1 + \tfrac19 + \tfrac1{25}\big) = 0.0669$; $N = 11$: $0.0337$; $N = 51$: $0.0078$. It goes to zero like $\frac{8}{\pi^2}\sum_{\text{odd }k > N}k^{-2}\approx\frac{4}{\pi^2N}$ (Parseval: $\frac{8}{\pi^2}\sum_{\text{odd}}k^{-2} = \frac{8}{\pi^2}\cdot\frac{\pi^2}{8} = 1$). Meanwhile the overshoot $\max_xS_Nf(x) - 1$ is $0.2732, 0.1884, 0.1813, 0.1791$ for the same $N$ (project, Part 2) and tends to $0.1790$, not to $0$. The worst-case error is worse still: $\sup_x|f(x) - S_Nf(x)| = 1$ for every $N$, because $S_Nf$ is continuous with $S_Nf(0) = 0$ while $f = 1$ on $(0,\pi)$, so $|f - S_Nf|\to1$ as $x\to0^+$. Theorem 3.5 controls the mean-square error only: it says nothing about either maximum, and §3.9 says neither goes to zero. The facts are consistent because the region where $S_Nf$ is far from $f$ has width $\sim1/N$ (§3.9), so it contributes $O(1/N)$ to $\frac{1}{2\pi}\|f - S_Nf\|^2$.

**B3.** (a) Triangle wave (Example 4.3): $a_0 = \pi$, $a_n = -\frac{4}{\pi n^2}$ for odd $n$, $0$ for even $n\ge2$, all $b_n = 0$. Left side of (3.11): $\frac1\pi\int_{-\pi}^\pi x^2\,dx = \frac1\pi\cdot\frac{2\pi^3}{3} = \frac{2\pi^2}{3}$. Right side: $\frac{a_0^2}{2} + \sum a_n^2 = \frac{\pi^2}{2} + \frac{16}{\pi^2}\sum_{n\text{ odd}}\frac1{n^4}$. So
$$\frac{16}{\pi^2}\sum_{n\text{ odd}}\frac1{n^4} = \frac{2\pi^2}{3} - \frac{\pi^2}{2} = \frac{\pi^2}{6} \quad\Longrightarrow\quad \sum_{n\text{ odd}}\frac{1}{n^4} = \frac{\pi^4}{96}\approx1.014678.$$
Even terms: $\sum_{n\text{ even}}\frac1{n^4} = \sum_{m\ge1}\frac{1}{(2m)^4} = \frac1{16}\sum_m\frac1{m^4} = \frac{\pi^4}{1440}\approx0.067645$. Check: $\frac{\pi^4}{90} - \frac{\pi^4}{96} = \pi^4\cdot\frac{96 - 90}{8640} = \frac{6\pi^4}{8640} = \frac{\pi^4}{1440}$. ✓

(b) Let $g$ be the odd extension of $f$ to $[-\pi,\pi]$. Then $g$ has only sine coefficients, and they are the half-range $b_n$ of §3.5 (both equal $\frac2\pi\int_0^\pi f\sin nx$). Also $\|g\|^2 = \int_{-\pi}^\pi g^2 = 2\int_0^\pi f^2$ because $g^2$ is even. Parseval (3.11) for $g$: $\frac1\pi\|g\|^2 = \sum b_n^2$, i.e. $\frac2\pi\int_0^\pi f^2\,dx = \sum_{n\ge1}b_n^2$.

Apply to $f = x(\pi - x)$ with $b_n = \frac{8}{\pi n^3}$ (odd $n$):
$$\int_0^\pi x^2(\pi - x)^2\,dx = \int_0^\pi\big(\pi^2x^2 - 2\pi x^3 + x^4\big)dx = \frac{\pi^5}{3} - \frac{\pi^5}{2} + \frac{\pi^5}{5} = \pi^5\cdot\frac{10 - 15 + 6}{30} = \frac{\pi^5}{30}.$$
Left side: $\frac2\pi\cdot\frac{\pi^5}{30} = \frac{\pi^4}{15}$. Right side: $\frac{64}{\pi^2}\sum_{n\text{ odd}}\frac1{n^6}$. So $\sum_{n\text{ odd}}\dfrac{1}{n^6} = \dfrac{\pi^4}{15}\cdot\dfrac{\pi^2}{64} = \dfrac{\pi^6}{960}\approx1.001447$.

Full sum: $\sum_{n\ge1}\frac1{n^6} = \sum_{\text{odd}} + \sum_{\text{even}} = \frac{\pi^6}{960} + \frac1{64}\sum_{n\ge1}\frac1{n^6}$, so $\big(1 - \frac1{64}\big)\sum\frac1{n^6} = \frac{\pi^6}{960}$ and $\sum_{n\ge1}\frac{1}{n^6} = \frac{64}{63}\cdot\frac{\pi^6}{960} = \frac{\pi^6}{945}\approx1.017343$ (Euler's $\zeta(6)$).

(c) For a function that is a polynomial on each piece, repeated integration by parts (Theorem 3.9) gives coefficients that are finite sums of terms $\frac{\text{rational}\times\pi^j}{n^k}$ with integer $k$; Parseval *squares* them, producing $\sum1/n^{2k}$ with $\pi^{2k}$ on the other side. Only even powers can arise this way. B1(b) came from *evaluating* a series at a point, not squaring it, and gave the odd power $n^{-3}$, but only in the *alternating* combination $\sum_{\text{odd}}(-1)^{(n-1)/2}/n^3 = \pi^3/32$. This route never reaches $\zeta(3) = \sum1/n^3$. For a polynomial on $[-\pi,\pi]$ the odd powers $n^{-3}, n^{-5},\dots$ arrive attached to $\sin nx$ (odd functions), and evaluating $\sum\sin(nx)/n^3$ at a point gives sums like B1(b)'s, never the plain $\sum 1/n^3$, which would need $\sum\cos(nx)/n^3$ at $x = 0$. Nobody knows a closed form for $\zeta(3)$. Apéry proved it irrational (announced 1978), but whether $\zeta(3)/\pi^3$ is rational is an open problem. Numerics only rule out fractions with small denominators.

**B4.** (a) $f(x) = \pi^4 - 2\pi^2x^2 + x^4$, an even polynomial, so continuity across the seam is the question $f^{(j)}(\pi^-) = f^{(j)}(-\pi^+)$, and for even $f^{(j)}$ that is automatic while for odd $f^{(j)}$ it requires $f^{(j)}(\pi) = 0$.
- $f(\pm\pi) = 0$: continuous. ✓
- $f'(x) = -4\pi^2x + 4x^3 = -4x(\pi^2 - x^2)$, $f'(\pm\pi) = 0$: continuous. ✓
- $f''(x) = -4\pi^2 + 12x^2$, even, $f''(\pm\pi) = 8\pi^2$: continuous. ✓
- $f'''(x) = 24x$, odd: $f'''(\pi^-) = 24\pi$, $f'''(-\pi^+) = -24\pi$. Jump at the seam $J = f'''(-\pi^+) - f'''(\pi^-) = -48\pi$.

(b) $f, f', f''$ are continuous, periodic and piecewise smooth, so Theorem 3.18 applies three times: $c_n(f') = in\,c_n(f)$, $c_n(f'') = in\,c_n(f')$, $c_n(f''') = in\,c_n(f'')$, giving $c_n(f''') = (in)^3c_n(f) = -in^3c_n(f)$. And $f''' = 24x$ on $(-\pi,\pi)$, so $c_n(f''') = 24\,c_n(x) = 24\cdot\frac{i(-1)^n}{n}$ (Example 4.2). Therefore
$$c_n(f) = \frac{c_n(f''')}{-in^3} = \frac{24i(-1)^n}{n\cdot(-in^3)} = -\frac{24(-1)^n}{n^4}\quad(n\ne0),$$
and since $f$ is even, $b_n = 0$ and $a_n = 2c_n = -\dfrac{48(-1)^n}{n^4}$.

Check with Theorem 3.9 on $g = f'''$: one jump $J = -48\pi$ at $x_j = \pi$, and $g' = 24$ has $c_n(g') = 0$ for $n\ne0$. (3.12): $c_n(g) = \frac{1}{in}\cdot\frac{1}{2\pi}(-48\pi)e^{-in\pi} = \frac{-24(-1)^n}{in} = \frac{24i(-1)^n}{n}$. ✓ Same as $24c_n(x)$.

(c) $\frac{a_0}{2} = \frac{1}{2\pi}\int_{-\pi}^\pi f\,dx = \frac1\pi\int_0^\pi\big(\pi^4 - 2\pi^2x^2 + x^4\big)dx = \frac1\pi\Big(\pi^5 - \frac{2\pi^5}{3} + \frac{\pi^5}{5}\Big) = \pi^4\cdot\frac{15 - 10 + 3}{15} = \frac{8\pi^4}{15}\approx51.9515$.

So $(\pi^2 - x^2)^2 = \dfrac{8\pi^4}{15} - 48\sum_{n\ge1}\dfrac{(-1)^n}{n^4}\cos nx$, with equality everywhere by Theorem 3.14 (continuous, piecewise smooth), indeed uniformly by Theorem 3.16.

(d) At $x = \pi$: $\cos n\pi = (-1)^n$, $(-1)^n(-1)^n = 1$, $f(\pi) = 0$: $0 = \frac{8\pi^4}{15} - 48\sum\frac1{n^4}$, so $\sum_{n\ge1}\frac1{n^4} = \frac{8\pi^4}{15\cdot48} = \frac{\pi^4}{90}$. ✓ At $x = 0$: $f(0) = \pi^4 = \frac{8\pi^4}{15} - 48\sum\frac{(-1)^n}{n^4}$, so $\sum_{n\ge1}\frac{(-1)^{n+1}}{n^4} = \frac{1}{48}\Big(\pi^4 - \frac{8\pi^4}{15}\Big) = \frac{7\pi^4}{720}\approx0.947033$.

(e) $f, f', f''$ continuous and periodic, $f'''$ piecewise smooth: Theorem 3.10 with $k = 2$, predicting $|c_n|\le C/n^{4}$ with $C = \frac{1}{2\pi}\Big(\sum|J_j| + \int_{-\pi}^\pi|f^{(4)}|dx\Big) = \frac{1}{2\pi}\big(48\pi + 24\cdot2\pi\big) = 48$. The actual value is $|c_n|n^4 = 24$, half the bound. The bound is not tight because it uses $|e^{-inx_j}| = 1$ and $|\int g'e^{-inx}|\le\int|g'|$ separately, whereas here the second term is exactly zero for $n\ne0$.

**B5.** (a) $D_N(x) = \sum_{n=-N}^Ne^{inx} = 1 + \sum_{n=1}^N\big(e^{inx} + e^{-inx}\big) = 1 + 2\sum_{n=1}^N\cos nx$. Cosine is even, so $D_N$ is even. Integrating term by term, $\int_{-\pi}^\pi e^{inx}dx = 2\pi\delta_{n0}$ (Theorem 3.3 with $m = 0$), so $\int_{-\pi}^\pi D_N = 2\pi$ and $\frac{1}{2\pi}\int D_N = 1$.

(b) From (3.13), $D_N(x) = 0$ where $\sin\big((N+\tfrac12)x\big) = 0$ but $\sin(x/2)\ne0$: $(N+\tfrac12)x = k\pi$, i.e. $x_k = \dfrac{2k\pi}{2N+1}$ for $k = 1, 2, \dots, 2N$ (at $k = 0$ and $k = 2N+1$ both sines vanish and $D_N = 2N+1$). So $D_N$ has $2N$ zeros in $(0, 2\pi)$. The central lobe runs from $-\frac{2\pi}{2N+1}$ to $\frac{2\pi}{2N+1}$: width $\frac{4\pi}{2N+1}\approx\frac{2\pi}{N}$, height $2N+1$, so its area is of order $2\pi$, consistent with (a). The side lobes alternate in sign with heights decaying like $1/\sin(x/2)$.

(c) $D_N$ is real, so $D_N^2 = \overline{D_N}D_N = \sum_{m,n}e^{-imx}e^{inx}$ and $\int_{-\pi}^\pi D_N^2\,dx = \sum_{|m|,|n|\le N}2\pi\delta_{mn} = 2\pi(2N+1)$. (Equivalently: Parseval for the trigonometric polynomial $D_N$, whose $2N+1$ coefficients are all $1$: $\frac{1}{2\pi}\|D_N\|^2 = 2N+1$.) For $N = 10$: $2\pi\cdot21 = 42\pi\approx131.95$.

(d) By evenness, $L_N = \frac{1}{2\pi}\int_{-\pi}^\pi|D_N| = \frac1\pi\int_0^\pi\frac{|\sin((N+\frac12)x)|}{\sin(x/2)}\,dx$. For $0 < x\le\pi$, $0 < \sin(x/2)\le x/2$, so $\frac{1}{\sin(x/2)}\ge\frac2x$ and
$$L_N\ge\frac1\pi\int_0^\pi\frac{2\,|\sin((N+\frac12)x)|}{x}\,dx = \frac2\pi\int_0^{(N+\frac12)\pi}\frac{|\sin u|}{u}\,du \qquad (u = (N+\tfrac12)x,\ dx/x = du/u).$$
Drop the last half-period and split the rest into the $N$ intervals $[(k-1)\pi, k\pi]$, on each of which $\frac1u\ge\frac{1}{k\pi}$:
$$\frac2\pi\int_0^{N\pi}\frac{|\sin u|}{u}du\ \ge\ \frac2\pi\sum_{k=1}^N\frac{1}{k\pi}\int_{(k-1)\pi}^{k\pi}|\sin u|\,du = \frac2\pi\sum_{k=1}^N\frac{2}{k\pi} = \frac{4}{\pi^2}\sum_{k=1}^N\frac1k,$$
using $\int_{(k-1)\pi}^{k\pi}|\sin u|du = 2$. Finally $\sum_{k=1}^N\frac1k\ge\int_1^{N+1}\frac{du}{u} = \ln(N+1)$ (each term $\frac1k$ is the height of a rectangle over $[k, k+1]$ lying above the curve $1/u$). So $L_N\ge\frac{4}{\pi^2}\ln(N+1)\to\infty$. ∎ Numerically (project): $L_1 = 1.436$, $L_{10} = 2.223$, $L_{100} = 3.139$, each above the bound $0.281, 0.972, 1.870$. A matching upper bound also holds, $L_N = \frac{4}{\pi^2}\ln N + O(1)$, so the lower bound has the right growth rate (the project's values exceed $\frac{4}{\pi^2}\ln N$ by $1.29$ at $N = 10$ and $1.27$ at $N = 100$, a nearly constant gap).

(e) From (3.14), $|S_Nf(x)| = \Big|\frac{1}{2\pi}\int f(x-t)D_N(t)\,dt\Big|\le\max|f|\cdot\frac{1}{2\pi}\int|D_N| = L_N\max|f|$. This is the best bound of its form: the function $f = \operatorname{sgn}D_N$ (or a continuous function close to it) gives $S_Nf(0)\approx L_N$. So the operator $S_N$ can magnify a function of size $1$ into a value of size $L_N\to\infty$ at a point, and the uniform boundedness principle (Banach–Steinhaus) turns "the norms $L_N$ are unbounded" into "there is a single continuous $f$ with $\sup_N|S_Nf(0)| = \infty$". That is the kind of function du Bois-Reymond constructed explicitly in 1873, long before the abstract principle existed, and it is the reason Dirichlet's theorem needs more than continuity. The culprit is the sign changes of $D_N$: $\int D_N = 2\pi$ is fixed but $\int|D_N|$ grows. A kernel $K_N\ge0$ with $\frac{1}{2\pi}\int K_N = 1$ has $\frac{1}{2\pi}\int|K_N| = 1$, so the corresponding averages satisfy $|\sigma_Nf|\le\max|f|$ for every $N$, and nothing can blow up. That is C3.

**B6.** (a) Theorem 3.17 (with $a_0 = 0$ for the sawtooth, so there is no $\frac{a_0x}{2}$ term):
$$\frac{x^2}{2} = \int_0^xt\,dt = 2\sum_{n\ge1}(-1)^{n+1}\int_0^x\frac{\sin nt}{n}\,dt = 2\sum_{n\ge1}(-1)^{n+1}\frac{1 - \cos nx}{n^2} = 2\sum_{n\ge1}\frac{(-1)^{n+1}}{n^2} - 2\sum_{n\ge1}\frac{(-1)^{n+1}\cos nx}{n^2}.$$
With $\sum(-1)^{n+1}/n^2 = \pi^2/12$: $\frac{x^2}{2} = \frac{\pi^2}{6} + 2\sum\frac{(-1)^n\cos nx}{n^2}$, so $x^2 = \frac{\pi^2}{3} + 4\sum_{n\ge1}\frac{(-1)^n}{n^2}\cos nx$. ✓ (Example 4.4.) The constant $2\sum(-1)^{n+1}/n^2 = \sum_nb_n/n$ is exactly what step 4 of Theorem 3.17 produces (it evaluates at $x = 0$). It could also have been found as the mean of $x^2/2$, namely $\frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{x^2}{2}dx = \frac{\pi^2}{6}$; doing it that way would *prove* $\sum(-1)^{n+1}/n^2 = \pi^2/12$ instead of using it.

(b) $x^2 - \frac{\pi^2}{3}$ has mean zero (its $a_0$ is $0$), so its integral from $0$ is periodic:
$$\int_0^x\Big(t^2 - \frac{\pi^2}{3}\Big)dt = \frac{x^3}{3} - \frac{\pi^2x}{3} = 4\sum_{n\ge1}\frac{(-1)^n}{n^2}\int_0^x\cos nt\,dt = 4\sum_{n\ge1}\frac{(-1)^n\sin nx}{n^3}.$$
Multiply by $-3$: $\pi^2x - x^3 = x(\pi^2 - x^2) = -12\sum\frac{(-1)^n\sin nx}{n^3} = 12\sum_{n\ge1}\frac{(-1)^{n+1}}{n^3}\sin nx$. ✓ This is $b_n = \frac{12(-1)^{n+1}}{n^3}$, the cubic row of the §3.7 table, obtained there by differentiating and here by integrating. (The project checks the partial sum against $x(\pi^2 - x^2)$ on a grid.)

(c) At $x = \pi/2$: left side $\frac\pi2\Big(\pi^2 - \frac{\pi^2}{4}\Big) = \frac{3\pi^3}{8}$. Right side: $\sin(n\pi/2) = 0$ for even $n$ and $(-1)^{(n-1)/2}$ for odd $n$, and $(-1)^{n+1} = 1$ for odd $n$, so it is $12\sum_{n\text{ odd}}\frac{(-1)^{(n-1)/2}}{n^3}$. Hence $\sum_{n\text{ odd}}\frac{(-1)^{(n-1)/2}}{n^3} = \frac{3\pi^3}{8\cdot12} = \frac{\pi^3}{32}$, the same value as B1(b). ✓ Two different functions, $x(\pi - x)$ on $[0,\pi]$ and $x(\pi^2 - x^2)$ on $[-\pi,\pi]$, evaluated at the same point, give the same sum because both series reduce there to the same alternating odd cubic sum with a different prefactor ($\frac8\pi\cdot\frac{\pi^3}{32} = \frac{\pi^2}{4}$ and $12\cdot\frac{\pi^3}{32} = \frac{3\pi^3}{8}$); consistency is a check that both coefficient computations are right.

(d) Termwise: $\frac{d}{dx}\Big[\frac{\pi^2}{3} + 4\sum\frac{(-1)^n\cos nx}{n^2}\Big] = -4\sum\frac{(-1)^n\sin nx}{n} = 2\cdot2\sum\frac{(-1)^{n+1}\sin nx}{n}$, which is $2\times$ the sawtooth series, i.e. the Fourier series of $2x$. Legitimate: $x^2$ is continuous as a periodic function ($\pi^2 = \pi^2$ at the seam) and piecewise smooth, and its derivative $2x$ is piecewise smooth, so Theorem 3.18 applies. By Theorem 3.14 the differentiated series converges to $2x$ on $(-\pi,\pi)$ and at $x = \pi$ to the midpoint of the jump of the periodic $2x$, $\frac{2\pi + (-2\pi)}{2} = 0$ (every term $\sin n\pi$ is $0$ there, consistent). The sawtooth fails the hypothesis: it is *not* continuous as a periodic function ($\pi\ne-\pi$ at the seam). Its termwise derivative $2\sum(-1)^{n+1}\cos nx$ has terms that do not tend to zero and diverges at every $x$; §3.11 explains that it is "trying to be" $1$ minus a train of delta functions at the jumps, which F3 makes precise.

**B7.** (a) $u = X(x)T(t)$: $XT' = kX''T$, so $\frac{T'}{kT} = \frac{X''}{X} = -\lambda$, a constant (a function of $t$ alone equal to a function of $x$ alone). The boundary conditions $u_x(0,t) = X'(0)T(t) = 0$ for all $t$ force $X'(0) = 0$ (else $T\equiv0$), and likewise $X'(L) = 0$. The $x$-problem is $X'' + \lambda X = 0$, $X'(0) = X'(L) = 0$.
- $\lambda < 0$, $\mu = \sqrt{-\lambda}$: $X = A\cosh\mu x + B\sinh\mu x$, $X' = \mu(A\sinh\mu x + B\cosh\mu x)$. $X'(0) = \mu B = 0$ gives $B = 0$; $X'(L) = \mu A\sinh\mu L = 0$ gives $A = 0$. Nothing.
- $\lambda = 0$: $X = Ax + B$, $X' = A$; $X'(0) = 0$ gives $A = 0$ and $X = B$ is a nonzero solution. **Eigenvalue $\lambda_0 = 0$, eigenfunction $X_0 = 1$.** (In the Dirichlet case of §3.12 this case died because $X(0) = 0$ killed $B$.)
- $\lambda > 0$: $X = A\cos\sqrt\lambda x + B\sin\sqrt\lambda x$, $X' = \sqrt\lambda(-A\sin\sqrt\lambda x + B\cos\sqrt\lambda x)$. $X'(0) = \sqrt\lambda B = 0$ gives $B = 0$; $X'(L) = -\sqrt\lambda A\sin(\sqrt\lambda L) = 0$ with $A\ne0$ gives $\sqrt\lambda L = n\pi$. So $\lambda_n = (n\pi/L)^2$, $X_n = \cos(n\pi x/L)$, $n\ge1$.
Then $T_n' = -k\lambda_nT_n$ gives $T_n = e^{-k(n\pi/L)^2t}$, and $T_0 = $ const.

(b) Superpose:
$$u(x,t) = \frac{a_0}{2} + \sum_{n\ge1}a_n\cos\frac{n\pi x}{L}\,e^{-k(n\pi/L)^2t}, \qquad a_n = \frac2L\int_0^L\phi(x)\cos\frac{n\pi x}{L}\,dx\ \ (n\ge0),$$
the cosine series of §3.5 (projection rule (3.5) with $\|\cos(n\pi x/L)\|^2 = L/2$ for $n\ge1$ and $\|1\|^2 = L$, which is why the constant term is $\frac{a_0}{2} = \frac1L\int_0^L\phi$). As $t\to\infty$ every term with $n\ge1$ decays and $u(x,t)\to\frac{a_0}{2} = \frac1L\int_0^L\phi(x)\,dx$, the **mean** of the initial temperature, uniformly in $x$.

Conservation: $\frac{d}{dt}\int_0^Lu\,dx = \int_0^Lu_t\,dx = k\int_0^Lu_{xx}\,dx = k\big[u_x\big]_0^L = k\big(u_x(L,t) - u_x(0,t)\big) = 0$. The total heat $\int_0^Lu\,dx$ is constant because insulated ends let none escape. The limit *must* then be the constant with the same integral as $\phi$, which is the mean. (Mode $n = 0$ is exactly the conserved quantity: $\frac{a_0}{2}\cdot L = \int\phi$, and it is the only mode with no decay.)

(c) $\frac{a_0}{2} = \frac1L\int_0^{L/2}u_0\,dx = \frac{u_0}{2}$. For $n\ge1$, split at the break $L/2$:
$$a_n = \frac2L\int_0^{L/2}u_0\cos\frac{n\pi x}{L}\,dx = \frac{2u_0}{L}\cdot\frac{L}{n\pi}\Big[\sin\frac{n\pi x}{L}\Big]_0^{L/2} = \frac{2u_0}{n\pi}\sin\frac{n\pi}{2}.$$
$\sin\frac{n\pi}{2} = 1, 0, -1, 0, 1, 0$ for $n = 1,\dots,6$, so
$$a_1 = \frac{2u_0}{\pi}\approx0.6366u_0,\ \ a_2 = 0,\ \ a_3 = -\frac{2u_0}{3\pi}\approx-0.2122u_0,\ \ a_4 = 0,\ \ a_5 = \frac{2u_0}{5\pi}\approx0.1273u_0,\ \ a_6 = 0.$$
$$u(x,t)\approx\frac{u_0}{2} + \frac{2u_0}{\pi}\Big[\cos\frac{\pi x}{L}e^{-\pi^2kt/L^2} - \frac13\cos\frac{3\pi x}{L}e^{-9\pi^2kt/L^2} + \frac15\cos\frac{5\pi x}{L}e^{-25\pi^2kt/L^2}\Big].$$
Even modes vanish because $\phi - \frac{u_0}{2}$ is **antisymmetric** about $x = L/2$ ($+\frac{u_0}{2}$ on the left, $-\frac{u_0}{2}$ on the right) while $\cos\frac{n\pi x}{L}$ is symmetric about $L/2$ for even $n$ and antisymmetric for odd $n$ ($\cos\frac{n\pi(L-x)}{L} = (-1)^n\cos\frac{n\pi x}{L}$). In Example 4.6 the basis was $\sin\frac{n\pi x}{L}$, which is symmetric for *odd* $n$ and antisymmetric for *even* $n$, and the same data had both parities relative to that basis (the mean was not removed by the basis), so even modes appeared there. Same initial data, different fate: with cold ends (Example 4.6) all the heat drains out and $u\to0$; with insulated ends it redistributes and $u\to\frac{u_0}{2}$. The slowest transient in both cases is mode $1$ with the same $\tau_1 = L^2/(\pi^2k)$ (for the copper rod of Example 4.6, $9.21$ s).

(d) Theorem 5.1's first integration by parts gives $\langle f, Lf\rangle = \big[\overline{f}f'\big]_0^L - \|f'\|^2$, and the boundary term vanishes under either Dirichlet ($f = 0$ at the ends) or Neumann ($f' = 0$ at the ends) conditions. So for an eigenfunction, $\lambda\|f\|^2 = -\|f'\|^2\le0$, with $\lambda = 0$ if and only if $f'\equiv0$, i.e. $f$ is constant. A nonzero constant satisfies the Neumann conditions, so $0$ is a Neumann eigenvalue; it violates the Dirichlet conditions, so $0$ is not a Dirichlet eigenvalue. Physically: a uniform temperature is a steady state with insulated ends, and is not a steady state with the ends held at $0$.

**B8.** Reference implementation (also at the end of `project.py`, after "Your turn", so you can compare after trying it yourself).

(a)
```python
def quartic(x):
    """(pi² - x²)²: f, f', f'' continuous and periodic, f''' jumps at ±pi -> |c_n| ~ 1/n^4 (B4)."""
    x = np.asarray(x, float)
    return (PI**2 - x**2) ** 2

n, c = fourier_coeffs(quartic, 40, M=2**16)
mask = n >= 1
slope4 = np.polyfit(np.log(n[mask]), np.log(np.abs(c[mask])), 1)[0]
assert abs(slope4 + 4) < 0.05
assert np.allclose(c[mask], -24 * (-1.0) ** n[mask] / n[mask] ** 4, atol=1e-8)   # B4(b)
assert np.isclose(c[n == 0][0].real, 8 * PI**4 / 15, atol=1e-8)                   # B4(c)
```
Output: fitted slope $-4.000$; $c_0 = 8\pi^4/15 = 51.9515$. The trapezoid rule with $M = 2^{16}$ points is accurate to $10^{-8}$ here because the aliasing error $\sum_{k\ne0}c_{n+kM}$ is of size $24/M^4$ (§5, "The DFT").

(b) The $O(N\times\text{grid})$ idea: both $x^2$ and $|x|$ are even, so $S_Nf = c_0 + \sum_{n\le N}a_n\cos nx$. Build the matrix whose row $n$ is $a_n\cos(nx)$ on the grid once, then `np.cumsum` along $n$ gives every partial sum at once.
```python
def terms_for_uniform_error(name, tol=0.01, Nmax=600):
    """Smallest N with max_x |f - S_N f| < tol, for the even waves 'parabola' and 'triangle'."""
    xf = np.linspace(-PI, PI, 8001)
    fx = WAVES[name](xf)
    nn_ = np.arange(1, Nmax + 1)
    a_n = 2 * exact_c(name, nn_).real                     # a_n = c_n + c_{-n} = 2 Re c_n
    terms = a_n[:, None] * np.cos(np.outer(nn_, xf))      # row n-1: a_n cos(n x)
    S = exact_c(name, [0])[0].real + np.cumsum(terms, axis=0)   # row N-1: S_N f on the grid
    err = np.abs(S - fx).max(axis=1)
    hits = np.nonzero(err < tol)[0]
    return int(nn_[hits[0]]) if hits.size else None
```
Output: $N = 400$ for $x^2$ and $N = 63$ for $|x|$. Theorem 3.16's tail bounds: for $x^2$, $|c_n| = 2/n^2$ for all $n\ne0$, so $\sum_{|n|>N}|c_n| = 4\sum_{n>N}n^{-2} < 4/N$, giving $N = 400$ for $0.01$: the bound is essentially sharp, because the maximum error of $x^2$ sits at the corner $x = \pm\pi$ where every term $\frac{4(-1)^n}{n^2}\cos n\pi = \frac{4}{n^2}$ has the same sign, so the tail adds up with no cancellation. For $|x|$, $|c_n| = \frac{2}{\pi n^2}$ for odd $n$ only, so $\sum_{|n|>N}|c_n| = \frac4\pi\sum_{\text{odd }n>N}n^{-2}\approx\frac{4}{\pi}\cdot\frac{1}{2N} = \frac{2}{\pi N}$, giving $N\approx64$; the code finds $63$, again at the corners $x = 0, \pm\pi$ where the odd cosines all have the same sign.

## C. Challenge

**C1.** (a) $\sin\big(n(\pi - t)\big) = \sin n\pi\cos nt - \cos n\pi\sin nt = -(-1)^n\sin nt = (-1)^{n+1}\sin nt$. So $S_Nf(\pi - t) = 2\sum_{n=1}^N(-1)^{n+1}\frac{(-1)^{n+1}\sin nt}{n} = 2\sum_{n=1}^N\frac{\sin nt}{n} = T_N(t)$.

(b) $T_N'(t) = 2\sum_{n=1}^N\cos nt = D_N(t) - 1$ by B5(a). Insert (3.13): $D_N(t) - 1 = \dfrac{\sin\big((N+\frac12)t\big) - \sin\frac t2}{\sin\frac t2}$, and $\sin A - \sin B = 2\cos\frac{A+B}{2}\sin\frac{A-B}{2}$ with $A = (N+\tfrac12)t$, $B = \tfrac t2$ gives $\frac{A+B}{2} = \frac{(N+1)t}{2}$, $\frac{A-B}{2} = \frac{Nt}{2}$:
$$T_N'(t) = \frac{2\cos\frac{(N+1)t}{2}\,\sin\frac{Nt}{2}}{\sin\frac t2}.$$
On $(0,\pi)$ the denominator is positive. $\sin\frac{Nt}{2} > 0$ for $0 < t < \frac{2\pi}{N}$, and $\cos\frac{(N+1)t}{2} > 0$ for $0 < t < \frac{\pi}{N+1}$, negative immediately after. Since $\frac{\pi}{N+1} < \frac{2\pi}{N}$, the first sign change of $T_N'$ is from $+$ to $-$ at $t_N = \dfrac{\pi}{N+1}$: the first maximum.

(c) $T_N(t_N) = 2\sum_{n=1}^N\frac{\sin(n\pi/(N+1))}{n}$. Put $\Delta u = \frac{\pi}{N+1}$ and $u_n = n\Delta u$; then $\frac{\sin(n\Delta u)}{n} = \frac{\sin u_n}{u_n}\Delta u$, so
$$T_N(t_N) = 2\sum_{n=1}^N\frac{\sin u_n}{u_n}\,\Delta u,$$
a Riemann sum for $2\int_0^\pi\frac{\sin u}{u}du$ with nodes $u_1, \dots, u_N$ spaced $\Delta u$ apart (it omits one cell of width $\Delta u$ at the right end, whose contribution is $O(1/N)$). The integrand $\frac{\sin u}{u}$ (equal to $1$ at $u = 0$) is continuous on $[0,\pi]$, so the Riemann sums converge: $\lim_NT_N(t_N) = 2\operatorname{Si}(\pi) = 2\times1.8519370 = 3.7038741$.

(d) $T_N(t_N) - (\pi - t_N)$: for $N = 11$, $T_{11}(\pi/12) = 3.4384$ and $\pi - \pi/12 = 2.8798$, difference $0.5586$; for $N = 51$, $T_{51}(\pi/52) = 3.6433$, $\pi - \pi/52 = 3.0812$, difference $0.5621$. Both tend to $2\operatorname{Si}(\pi) - \pi = 0.5623 = 0.0895\times2\pi$: nine percent of the jump $2\pi$, exactly as for the square wave's jump of $2$.

The project's grid maxima of $S_Nf - f$ are $0.5643$ ($N = 11$) and $0.5624$ ($N = 51$), slightly different. The reason: $t_N$ is where $T_N$ itself peaks ($T_N' = 0$), but the overshoot *above $f$* is $T_N(t) - (\pi - t)$, which peaks where $T_N'(t) = -1$, a little to the right of $t_N$ (where $T_N'$ has turned negative). Near $t_N$, $T_N'' = O(N^2)$ (the kernel $D_N$ varies on the scale $1/N$ with height $N$), so the shift is $\delta\sim1/N^2$ and the change in $T_N$ is $\sim N^2\delta^2\sim1/N^2$, while the change in $\pi - t$ is $\delta\sim1/N^2$: the two measurements differ by $O(1/N^2)$. Check: $0.0057$ at $N = 11$ versus $0.0003$ at $N = 51$, ratio $\approx19$, and $(51/11)^2\approx21$. ✓ Either way the limit is $0.5623$.

(e) The square wave is the wrong tool for removing jumps: $\operatorname{sgn}$ jumps twice per period (by $+2$ at $0$ and $-2$ at $\pi$), so subtracting a shifted copy at one jump creates a new jump half a period away. Use instead the sawtooth with a single jump per period,
$$\sigma(x) = \frac{\pi - x}{2\pi}\ \text{ on } (0, 2\pi),\ \text{extended } 2\pi\text{-periodically}: \qquad \sigma(0^+) = \tfrac12,\ \ \sigma(0^-) = -\tfrac12,$$
a jump of $+1$ at each multiple of $2\pi$, linear with slope $-\frac{1}{2\pi}$ everywhere else.

1. **Partial sums of $\sigma$.** For $t\in(0,2\pi)$, $\sigma(t) = \frac{f(\pi - t)}{2\pi}$ with $f$ the sawtooth $x$. Substituting $y = \pi - x$ in (3.8) gives $c_n\big(f(\pi - \cdot)\big) = (-1)^nc_{-n}(f)$, hence $S_N[f(\pi - \cdot)](t) = \sum_{|n|\le N}(-1)^nc_{-n}(f)e^{int} = \sum_{|m|\le N}c_m(f)e^{im(\pi - t)} = S_Nf(\pi - t)$ (put $m = -n$ and use $(-1)^m = e^{im\pi}$). So by (a), $S_N\sigma(t) = \frac{T_N(t)}{2\pi}$, and the error of the series for $2\pi\sigma$ is
   $$E_N(t) = T_N(t) - (\pi - t), \qquad 0 < t < 2\pi.$$
2. **A formula for $E_N$.** $T_N(0) = 0$ and $T_N' = D_N - 1$ by (b), so $T_N(t) = \int_0^tD_N(s)\,ds - t$ and $E_N(t) = \int_0^tD_N(s)\,ds - \pi$. Since $D_N$ is even with $\frac{1}{2\pi}\int_{-\pi}^{\pi}D_N = 1$, $\int_0^\pi D_N = \pi$, so $E_N(t) = -\int_t^\pi D_N(s)\,ds$ for $0 < t\le\pi$.
3. **Uniform convergence away from the jump.** Fix $\delta\in(0,\pi)$ and $t\in[\delta,\pi]$. Integrate by parts, integrating $\sin((N+\frac12)s)$ and differentiating $\frac{1}{\sin(s/2)}$:
   $$\int_t^\pi\frac{\sin((N+\frac12)s)}{\sin\frac s2}ds = \Big[-\frac{\cos((N+\frac12)s)}{(N+\frac12)\sin\frac s2}\Big]_t^\pi - \frac{1}{N+\frac12}\int_t^\pi\cos\big((N+\tfrac12)s\big)\frac{\cos\frac s2}{2\sin^2\frac s2}\,ds.$$
   The boundary term at $\pi$ is $0$ ($\cos((N+\frac12)\pi) = 0$) and at $t$ is at most $\frac{1}{(N+\frac12)\sin(\delta/2)}$; the integral is at most $\frac{1}{N+\frac12}\int_\delta^\pi\frac{\cos(s/2)}{2\sin^2(s/2)}ds = \frac{1}{N+\frac12}\Big(\frac{1}{\sin(\delta/2)} - 1\Big)$. So $|E_N(t)| < \frac{2}{(N+\frac12)\sin(\delta/2)}$ on $[\delta,\pi]$. $T_N$ is odd and $2\pi$-periodic, so $E_N(2\pi - t) = -E_N(t)$ and the same bound holds on $[\pi, 2\pi - \delta]$. Hence $S_N\sigma\to\sigma$ uniformly on every closed set at distance $\ge\delta$ from the jumps.
4. **The peak next to the jump.** On $(0,\pi]$, $E_N' = D_N$. The zeros of $D_N$ there are $t_k = \frac{2k\pi}{2N+1}$ (B5(b)); $D_N > 0$ on $(0, t_1)$ and alternates in sign on successive intervals $[t_k, t_{k+1}]$. With $u = (N+\frac12)s$ the area of lobe $k$ is $A_k = \int_{k\pi}^{(k+1)\pi}\frac{|\sin u|\,du}{(N+\frac12)\sin\frac{u}{2N+1}}$, and since $\sin\frac{u}{2N+1}$ increases for $\frac{u}{2N+1}\le\frac\pi2$, the areas decrease: $A_1 > A_2 > \cdots$. After $t_1$, $E_N$ first drops by $A_1$, then rises by at most $A_2 < A_1$, and so on (a partial lobe contributes less than the full one), so $E_N(t)\le E_N(t_1)$ on $(0,\pi]$: the maximum is at the first zero $t_1 = \frac{2\pi}{2N+1}$, the location found in (d). There, with $u = (N+\frac12)s$,
   $$E_N(t_1) = \int_0^\pi\frac{\sin u}{(N+\frac12)\sin\frac{u}{2N+1}}\,du - \pi\ \longrightarrow\ 2\int_0^\pi\frac{\sin u}{u}\,du - \pi = 2\operatorname{Si}(\pi) - \pi,$$
   because $(N+\frac12)\sin\frac{u}{2N+1}\to\frac u2$ uniformly on $[0,\pi]$ (step 5 of Theorem 3.15).
5. **Assemble.** Let $x_0 = x_1, x_2, \dots, x_m$ be the jumps of $f$ in one period, with sizes $J_1 = J > 0, J_2, \dots, J_m$, and set
   $$h(x) = f(x) - \sum_{j=1}^mJ_j\,\sigma(x - x_j).$$
   Each $J_j\sigma(\cdot - x_j)$ has exactly one jump per period, of size $J_j$ at $x_j$, and is linear elsewhere, so $h$ is continuous, $2\pi$-periodic and piecewise smooth, and $S_Nh\to h$ uniformly by Theorem 3.16. Choose $\delta_0 > 0$ less than the distance from $x_0$ to every other jump and let $I = [x_0 - \delta_0, x_0 + \delta_0]$. By step 3 (partial sums commute with shifts, since $c_n(\sigma(\cdot - a)) = e^{-ina}c_n(\sigma)$), $S_N[\sigma(\cdot - x_j)]\to\sigma(\cdot - x_j)$ uniformly on $I$ for $j\ge2$. So $G_N = S_Nh + \sum_{j\ge2}J_jS_N[\sigma(\cdot - x_j)]$ converges uniformly on $I$ to $g = h + \sum_{j\ge2}J_j\sigma(\cdot - x_j)$, which is continuous on $I$, and $f = g + J\sigma(\cdot - x_0)$, so $f(x_0^+) = g(x_0) + \frac J2$. For $x = x_0 + t$ with $0 < t < \delta\le\delta_0$, using $S_N[\sigma(\cdot - x_0)](x) = \frac{T_N(t)}{2\pi}$ and $\frac{J}{2\pi}T_N(t) - \frac J2 = \frac{J}{2\pi}E_N(t) - \frac{Jt}{2\pi}$:
   $$S_Nf(x) - f(x_0^+) = \underbrace{\big[G_N(x) - g(x)\big]}_{\to0\text{ uniformly}} + \underbrace{\Big[g(x) - g(x_0) - \frac{Jt}{2\pi}\Big]}_{|\cdot|\,\le\,\omega(\delta)} + \frac{J}{2\pi}E_N(t),$$
   where $\omega(\delta) = \max_{0\le t\le\delta}|g(x_0 + t) - g(x_0)| + \frac{J\delta}{2\pi}\to0$ as $\delta\to0$ by continuity of $g$. Since $J > 0$, step 4 gives $\max_{0<t<\delta}\frac{J}{2\pi}E_N(t) = \frac{J}{2\pi}E_N(t_1)\to\frac{J}{2\pi}\big(2\operatorname{Si}(\pi) - \pi\big)$ once $t_1 < \delta$. So both $\limsup_N$ and $\liminf_N$ of $\max_{x_0<x<x_0+\delta}\big(S_Nf(x) - f(x_0^+)\big)$ lie within $\omega(\delta)$ of
   $$\frac{J}{2\pi}\big(2\operatorname{Si}(\pi) - \pi\big) = \Big(\frac{\operatorname{Si}(\pi)}{\pi} - \frac12\Big)J\approx0.0895\,J,$$
   and letting $\delta\to0$ proves the claim. ∎

On the left of $x_0$ the same computation with $t\in(-\delta, 0)$ uses $\sigma(t) = \frac{-\pi - t}{2\pi}$ there and the oddness of $T_N$: the error term becomes $-\frac{J}{2\pi}E_N(|t|)$, an undershoot of $f(x_0^-)$ by the same $0.0895J$. A downward jump is the case $J > 0$ for $-f$. For the sawtooth itself, $J = -2\pi$ at $x_0 = \pi$: applying the result to $-f$, $S_Nf$ overshoots $f(\pi^-) = \pi$ from the left by $0.0895\times2\pi = 0.5623$, as in (d). Note what the order of limits buys: for fixed $\delta$ the bracket $\omega(\delta)$ need not vanish (the function itself varies on $(x_0, x_0 + \delta)$), which is why the statement sends $\delta\to0$ after $N\to\infty$.

**C2.** (a) Coefficients are linear: $c_n(f + g) = c_n(f) + c_n(g)$. Parseval (3.11) in the form $\|u\|^2 = 2\pi\sum|c_n(u)|^2$, applied to $u = f + g$:
$$\|f\|^2 + 2\operatorname{Re}\langle f, g\rangle + \|g\|^2 = 2\pi\sum_n\Big(|c_n(f)|^2 + 2\operatorname{Re}\big(\overline{c_n(f)}c_n(g)\big) + |c_n(g)|^2\Big).$$
Subtract Parseval for $f$ and for $g$: $\operatorname{Re}\langle f, g\rangle = 2\pi\sum\operatorname{Re}\big(\overline{c_n(f)}c_n(g)\big)$. Now replace $g$ by $ig$: $c_n(ig) = ic_n(g)$, $\langle f, ig\rangle = i\langle f, g\rangle$ (linearity in the second slot), and $\operatorname{Re}(iz) = -\operatorname{Im}z$. The same identity becomes $-\operatorname{Im}\langle f, g\rangle = -2\pi\sum\operatorname{Im}\big(\overline{c_n(f)}c_n(g)\big)$. Real and imaginary parts both agree, so $\langle f, g\rangle = 2\pi\sum_n\overline{c_n(f)}\,c_n(g)$. ∎ (This is the polarization identity of LA-8 Exercise C2: an inner product is determined by its norm.)

(b) $x$ and $y$ are continuous, $2\pi$-periodic and piecewise smooth, so by Theorem 3.18 (the continuous case of (3.12)) $c_n(x') = in\alpha_n$ and $c_n(y') = in\beta_n$. Parseval for $x'$ and $y'$, added:
$$\frac{1}{2\pi}\int_0^{2\pi}\big(x'^2 + y'^2\big)ds = \sum_n\big(|in\alpha_n|^2 + |in\beta_n|^2\big) = \sum_nn^2\big(|\alpha_n|^2 + |\beta_n|^2\big).$$
The left side is $\frac{1}{2\pi}\int_0^{2\pi}1\,ds = 1$ because the parametrization is by arc length. ∎

(c) Green's theorem for the region $\Omega$ bounded by $\gamma$ (positively oriented): $\oint_\gamma P\,dx + Q\,dy = \iint_\Omega(Q_x - P_y)\,dA$. With $P = 0$, $Q = x$: $\oint x\,dy = \iint1\,dA = A$. In the parameter, $\oint x\,dy = \int_0^{2\pi}x(s)y'(s)\,ds$. Since $x$ is real, $\int xy'\,ds = \int\overline{x}\,y'\,ds = \langle x, y'\rangle$, and by (a) with $c_n(y') = in\beta_n$:
$$A = \langle x, y'\rangle = 2\pi\sum_n\overline{\alpha_n}\,(in\beta_n) = 2\pi\sum_nin\,\overline{\alpha_n}\beta_n.$$

(d) Triangle inequality, then $|n|\le n^2$ (true for every integer, with equality only for $n\in\{0, \pm1\}$), then $2|\alpha||\beta|\le|\alpha|^2 + |\beta|^2$ (i.e. $(|\alpha| - |\beta|)^2\ge0$):
$$|A|\le2\pi\sum_n|n|\,|\alpha_n||\beta_n|\le2\pi\sum_nn^2|\alpha_n||\beta_n|\le2\pi\sum_nn^2\,\frac{|\alpha_n|^2 + |\beta_n|^2}{2} = \pi\cdot1 = \pi,$$
using (b) at the last step. ∎

Equality analysis. Equality in the middle inequality requires, for each $n$, either $|n| = n^2$ or $\alpha_n\beta_n = 0$; equality in the last requires $|\alpha_n| = |\beta_n|$ for every $n$ with $n\ne0$. For $|n|\ge2$ both together force $\alpha_n = \beta_n = 0$. So $x(s) = \alpha_0 + \alpha_1e^{is} + \overline{\alpha_1}e^{-is} = x_0 + 2\rho\cos(s + \varphi)$ and $y(s) = y_0 + 2\rho'\cos(s + \psi)$ (writing $\alpha_1 = \rho e^{i\varphi}$, $\beta_1 = \rho'e^{i\psi}$, with $\alpha_{-1} = \overline{\alpha_1}$ since $x$ is real), with $\rho = \rho'$. From (b): $1\cdot(|\alpha_1|^2 + |\beta_1|^2) + 1\cdot(|\alpha_{-1}|^2 + |\beta_{-1}|^2) = 4\rho^2 = 1$, so $2\rho = 1$. Arc length: $x'^2 + y'^2 = \sin^2(s + \varphi) + \sin^2(s + \psi) = 1$ for all $s$ forces $\psi = \varphi\pm\frac\pi2$, i.e. $y - y_0 = \mp\sin(s + \varphi)$. Then $(x - x_0)^2 + (y - y_0)^2 = \cos^2(s+\varphi) + \sin^2(s + \varphi) = 1$: a circle of radius $1$, area $\pi$, length $2\pi$. ∎

(e) If $\gamma$ has length $L$, the curve $\frac{2\pi}{L}\gamma$ has length $2\pi$ and area $\big(\frac{2\pi}{L}\big)^2A\le\pi$, so $A\le\dfrac{L^2}{4\pi}$, with equality only for circles.

Ellipse $x = 2\cos\theta$, $y = \sin\theta$: $A = \pi ab = 2\pi\approx6.2832$. $L = \int_0^{2\pi}\sqrt{4\sin^2\theta + \cos^2\theta}\,d\theta = 9.6884$ (no closed form; `scipy.integrate.quad`), so $L^2/4\pi = 7.4696 > 6.2832$. ✓ With $\theta$ as parameter: $\alpha_{\pm1} = 1$ ($2\cos\theta = e^{i\theta} + e^{-i\theta}$) and $\beta_{\pm1} = \mp\frac i2$ ($\sin\theta = \frac{e^{i\theta} - e^{-i\theta}}{2i}$). Then $2\pi\sum in\overline{\alpha_n}\beta_n = 2\pi\big[i\cdot1\cdot1\cdot(-\tfrac i2) + i\cdot(-1)\cdot1\cdot(\tfrac i2)\big] = 2\pi\big[\tfrac12 + \tfrac12\big] = 2\pi$. ✓ (The project does this with the FFT of the sampled curve.) Formula (c) holds for any parametrization because $\oint x\,dy$ is a line integral, independent of how the curve is traversed in time, and the polarized Parseval (a) is an identity for any two periodic functions; the arc-length condition entered only in (b), where it turned $\frac{1}{2\pi}\int(x'^2 + y'^2)$ into the number $1$. With the $\theta$-parametrization, $\sum n^2(|\alpha_n|^2 + |\beta_n|^2) = 2\cdot1 + 2\cdot\frac14 = \frac52\ne1$, and the inequality chain in (d) would give $|A|\le\frac52\pi$, true but useless; arc length is what makes the bound sharp.

**C3.** (a) $S_kf = \sum_{|n|\le k}c_ne^{inx}$. In $\sum_{k=0}^NS_kf$, the term $c_ne^{inx}$ appears once for each $k$ with $|n|\le k\le N$, which is $N + 1 - |n|$ values of $k$. So $\sum_{k=0}^NS_kf = \sum_{|n|\le N}(N+1-|n|)c_ne^{inx}$, and dividing by $N+1$ gives $\sigma_Nf = \sum_{|n|\le N}\Big(1 - \frac{|n|}{N+1}\Big)c_ne^{inx}$: the Dirichlet coefficients with a triangular taper.

(b) Average (3.14) over $k = 0, \dots, N$: $\sigma_Nf(x) = \frac{1}{2\pi}\int_{-\pi}^\pi f(x-t)\,\frac{1}{N+1}\sum_{k=0}^ND_k(t)\,dt$. With $D_k(t) = \frac{\sin((k+\frac12)t)}{\sin(t/2)}$:
$$\sum_{k=0}^ND_k(t) = \frac{1}{\sin\frac t2}\sum_{k=0}^N\sin\big((k+\tfrac12)t\big).$$
Multiply the sum by $2\sin\frac t2$ and use $2\sin A\sin B = \cos(A-B) - \cos(A+B)$ with $A = (k+\frac12)t$, $B = \frac t2$: $2\sin\frac t2\sin\big((k+\tfrac12)t\big) = \cos kt - \cos(k+1)t$. The sum telescopes: $\sum_{k=0}^N\big(\cos kt - \cos(k+1)t\big) = 1 - \cos(N+1)t = 2\sin^2\frac{(N+1)t}{2}$. So $\sum_{k=0}^N\sin\big((k+\tfrac12)t\big) = \dfrac{\sin^2\frac{(N+1)t}{2}}{\sin\frac t2}$ and
$$F_N(t) = \frac{1}{N+1}\sum_{k=0}^ND_k(t) = \frac{1}{N+1}\cdot\frac{\sin^2\frac{(N+1)t}{2}}{\sin^2\frac t2}, \qquad F_N(0) = \frac{1}{N+1}\sum_k(2k+1) = N+1.$$

(c) $F_N\ge0$ as a quotient of squares (and $N+1 > 0$ at $t = 0$). $\frac{1}{2\pi}\int F_N = \frac{1}{N+1}\sum_{k=0}^N\frac{1}{2\pi}\int D_k = \frac{1}{N+1}\cdot(N+1) = 1$ by Theorem 3.12.

(d) By (b), $\sigma_Nf(x) = \frac{1}{2\pi}\int f(x-t)F_N(t)\,dt$ is an average of values of $f$ with the nonnegative weights $\frac{1}{2\pi}F_N(t)\,dt$ that sum to $1$. An average of numbers in $[\min f, \max f]$ lies in $[\min f, \max f]$: $\min f\le\sigma_Nf(x)\le\max f$ for every $x$ and $N$. ∎ For the square wave, $|\sigma_Nf|\le1$ everywhere; the project finds $\max_x\sigma_{11}f = 0.9473$ and $\max_x\sigma_{51}f = 0.9878$, an *undershoot* of $0.0527$ and $0.0122$ that shrinks like $1/N$, compared with $S_N$'s permanent overshoot of $0.18$. Reading (a): Fejér summation multiplies $c_n$ by $1 - \frac{|n|}{N+1}$, which damps exactly the high harmonics whose abrupt truncation makes the Gibbs ears.

(e) *Fejér's theorem.* Let $f$ be continuous and $2\pi$-periodic, hence bounded by some $B$ and uniformly continuous. Given $\varepsilon > 0$ pick $\delta > 0$ with $|f(x - t) - f(x)| < \varepsilon$ whenever $|t| < \delta$. Using (c) to write $f(x) = \frac{1}{2\pi}\int f(x)F_N(t)\,dt$,
$$|\sigma_Nf(x) - f(x)| = \Big|\frac{1}{2\pi}\int_{-\pi}^\pi\big(f(x-t) - f(x)\big)F_N(t)\,dt\Big| \le \frac{1}{2\pi}\int_{|t|<\delta}\varepsilon F_N\,dt + \frac{1}{2\pi}\int_{\delta\le|t|\le\pi}2B\,F_N\,dt.$$
The first term is at most $\varepsilon$ (unit mass, positivity). For the second, on $\delta\le|t|\le\pi$ we have $\sin^2\frac t2\ge\sin^2\frac\delta2$ and $\sin^2\frac{(N+1)t}{2}\le1$, so $F_N(t)\le\frac{1}{(N+1)\sin^2(\delta/2)}$, and the second term is at most $\frac{2B}{(N+1)\sin^2(\delta/2)}$, which is $< \varepsilon$ for $N$ large. Neither bound depends on $x$. So $\max_x|\sigma_Nf(x) - f(x)| < 2\varepsilon$ for all large $N$: uniform convergence. ∎ (The three properties used, positivity, unit mass, and mass concentrating at $0$, define a *good kernel*; $D_N$ fails the first, which is B5(d)–(e).)

*Completeness.* Each $\sigma_Nf$ is a trigonometric polynomial of degree $\le N$. Uniform convergence implies mean-square convergence: $\|f - \sigma_Nf\|^2\le2\pi\max|f - \sigma_Nf|^2\to0$. By Theorem 3.5, $S_Nf$ is the *closest* element of $\mathcal{T}_N$ to $f$, so $\|f - S_Nf\|\le\|f - \sigma_Nf\|\to0$. That is mean-square convergence of the Fourier series for every continuous $f$, which by the proof of Theorem 3.8 is Parseval for continuous $f$. For $f\in PC$ with jumps, replace $f$ near each jump by a steep linear ramp over an interval of width $\eta$; the modified $g$ is continuous, $\|f - g\|^2\le(\text{number of jumps})\times(\max|f|)^2\times4\eta\to0$ as $\eta\to0$, and $\|f - S_Nf\|\le\|f - S_Ng\| \le\|f - g\| + \|g - S_Ng\|$ (using $S_Nf$ is the best approximation and $S_Ng\in\mathcal{T}_N$), both pieces small. So the trigonometric system is complete in $PC$, and Parseval (3.11) holds for every $f\in PC$: the one thing §3.6 took on faith, now proved in half a page given Fejér. Stein & Shakarchi ch. 2 §5 and ch. 3 §1 present the same argument for Riemann-integrable functions.

(f) Add to `project.py` (the reference is in Part 9):
```python
def fejer_sum(n, c, x):
    """sigma_N f(x) = sum_{|n|<=N} (1 - |n|/(N+1)) c_n e^{inx}, C3(a)."""
    N = n.max()
    return np.real(np.exp(1j * np.outer(x, n)) @ ((1 - np.abs(n) / (N + 1)) * c))

N = 51; n = np.arange(-N, N + 1)
sig = fejer_sum(n, exact_c("square", n), xg)
assert sig.max() <= 1 + 1e-12 and sig.min() >= -1 - 1e-12    # no overshoot, C3(d)
```
Plot it next to $S_{51}f$: the ears are gone, and the price is a slightly wider transition at the jump (the Fejér kernel's central lobe is about twice as wide as the Dirichlet kernel's for the same $N$).
