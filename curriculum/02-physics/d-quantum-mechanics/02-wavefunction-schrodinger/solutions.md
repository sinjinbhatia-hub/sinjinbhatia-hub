# QM-2 Solutions

Every number below is printed or asserted in Part 11 of `project.py` (or the part named).

## A. Warm-up

**A1.** (a) The time factor has modulus 1. $\int|A|^2e^{-2\lambda|x|}dx = 2|A|^2\int_0^\infty e^{-2\lambda x}dx = |A|^2/\lambda = 1$, so $A = \sqrt\lambda$.

(b) $\langle x\rangle = 0$: $x e^{-2\lambda|x|}$ is odd. $\langle x^2\rangle = 2\lambda\int_0^\infty x^2e^{-2\lambda x}dx = 2\lambda\cdot\frac{2}{(2\lambda)^3} = \frac{1}{2\lambda^2}$, using $\int_0^\infty x^ne^{-bx}dx = n!/b^{n+1}$. So $\sigma_x = \frac{1}{\sqrt2\,\lambda}$.

(c) $P(|x| > \sigma_x) = 2\int_{\sigma_x}^\infty\lambda e^{-2\lambda x}dx = e^{-2\lambda\sigma_x} = e^{-\sqrt2} = 0.2431$. (A Gaussian would give $0.3173$; the exponential has more mass near the center and in the far tails, less in between.)

**A2.** (a) For $\Psi = e^{i(kx-\omega t)}$: $i\hbar\partial_t\Psi = \hbar\omega\Psi$ and $-\frac{\hbar^2}{2m}\partial_x^2\Psi = \frac{\hbar^2k^2}{2m}\Psi$. The free equation holds for all $x, t$ iff $\hbar\omega\Psi = \frac{\hbar^2k^2}{2m}\Psi$, and since $\Psi \ne 0$, iff $\omega = \hbar k^2/2m$.

(b) $k = 2\pi/\lambda = 2\pi/10^{-10} = 6.283\times10^{10}$ m$^{-1}$. $v_g = \hbar k/m_e = (1.0546\times10^{-34})(6.283\times10^{10})/(9.109\times10^{-31}) = 7.274\times10^6$ m/s, and $v_p = v_g/2 = 3.637\times10^6$ m/s. $E = (\hbar k)^2/2m_e = 2.410\times10^{-17}$ J $= 150.4$ eV. (This is why electron diffraction from crystals, with spacings of a few Å, uses electrons of order 100 eV, as in Davisson–Germer, QM-1.)

**A3.** (a) $|e^{i\theta}\Psi|^2 = |e^{i\theta}|^2|\Psi|^2 = |\Psi|^2$. In $\langle Q\rangle = \int(e^{i\theta}\Psi)^*Q(x, -i\hbar\partial_x)(e^{i\theta}\Psi)\,dx$, the constant $e^{i\theta}$ passes through $\partial_x$, so the integrand is $e^{-i\theta}e^{i\theta}\Psi^*Q\Psi = \Psi^*Q\Psi$.

(b) $|e^{i\theta(x)}|^2 = 1$ still, so $|\Psi|^2$ and every $\langle f(x)\rangle$ are unchanged. But $\partial_x(e^{i\theta}\Psi) = e^{i\theta}(\Psi' + i\theta'\Psi)$, so
$\langle p\rangle_{\rm new} = -i\hbar\int\Psi^*(\Psi' + i\theta'\Psi)\,dx = \langle p\rangle + \hbar\int\theta'|\Psi|^2dx = \langle p\rangle + \hbar\langle\theta'\rangle.$
For $\theta = k_0x$ every momentum component is shifted by $\hbar k_0$ (in momentum space $\Phi(p)\to\Phi(p - \hbar k_0)$): the state is boosted by velocity $\hbar k_0/m$. That is how (3.23) gives the packet (3.16) a mean momentum. The phase of $\Psi$ is invisible point by point; its gradient is momentum.

**A4.** (a) $\hbar/m_ec = (1.0546\times10^{-34})/((9.109\times10^{-31})(2.998\times10^8)) = 3.862\times10^{-13}$ m. Confining to $\sigma_x < \hbar/m_ec$ forces $\sigma_p > \hbar/2\sigma_x > m_ec/2$ by (3.15), so the momentum spread is relativistic, kinetic energies reach $\sim m_ec^2$, and pair creation becomes possible: a one-particle, non-relativistic equation cannot describe that.

(b) $m = 4.0026\times10^{-3}/N_A = 6.646\times10^{-27}$ kg. $\lambda_{\rm th} = h/\sqrt{2\pi mk_BT} = 6.626\times10^{-34}/\sqrt{2\pi(6.646\times10^{-27})(1.381\times10^{-23})(4.2)} = 4.26\times10^{-10}$ m. $n = 125/m = 1.88\times10^{28}$ m$^{-3}$, spacing $n^{-1/3} = 3.76\times10^{-10}$ m. $n\lambda_{\rm th}^3 = 1.45$. Not classical: the thermal wavelength exceeds the spacing, so the atoms' wavefunctions overlap. Helium-4 becomes a superfluid below 2.17 K, and it never freezes at atmospheric pressure, because zero-point motion (Theorem 3.10(b)) outweighs the weak attraction.

**A5.** (a) For nice $f, g$, using Hermiticity of $\hat x$ and then of $\hat p$ (Theorem 3.4): $\langle f, \hat x\hat pg\rangle = \langle\hat xf, \hat pg\rangle = \langle\hat p\hat xf, g\rangle$. So $(\hat x\hat p)^\dagger = \hat p\hat x$. And $(\hat p\hat x)f = -i\hbar(xf)' = -i\hbar f - i\hbar xf' = (\hat x\hat p - i\hbar)f$. So $(\hat x\hat p)^\dagger = \hat x\hat p - i\hbar \ne \hat x\hat p$.

(b) $\langle\Psi, \hat x\hat p\Psi\rangle = -i\hbar\int x\Psi\Psi'dx = -\frac{i\hbar}{2}\int x\,(\Psi^2)'dx = -\frac{i\hbar}{2}\Big([x\Psi^2]_{-\infty}^{\infty} - \int\Psi^2dx\Big) = \frac{i\hbar}{2}.$
It is purely imaginary. A measured quantity cannot average to an imaginary number, so $\hat x\hat p$ is not an observable, confirming (a) by Theorem 3.4(4).

(c) $\big(\frac12(\hat x\hat p + \hat p\hat x)\big)^\dagger = \frac12\big((\hat x\hat p)^\dagger + (\hat p\hat x)^\dagger\big) = \frac12(\hat p\hat x + \hat x\hat p)$, using $(\hat p\hat x)^\dagger = \hat x^\dagger\hat p^\dagger = \hat x\hat p$ by the argument of (a) with the roles swapped. For a real $\Psi$ its expectation is $\frac12(\frac{i\hbar}{2} - \frac{i\hbar}{2}) = 0$, since $\langle\hat p\hat x\rangle = \langle\hat x\hat p\rangle - i\hbar$. This symmetrized product is the right quantum version of the classical $xp$ (rule (3.9) needs this choice whenever $x$ and $p$ are multiplied).

## B. Core

**B1.** Differentiate under the integral:
$$\frac{d\langle p\rangle}{dt} = -i\hbar\int\big(\partial_t\Psi^*\,\partial_x\Psi + \Psi^*\,\partial_x\partial_t\Psi\big)dx.$$
Integrate the second term by parts: $\int\Psi^*\partial_x(\partial_t\Psi)dx = [\Psi^*\partial_t\Psi]_{-\infty}^\infty - \int\partial_x\Psi^*\,\partial_t\Psi\,dx$. Dropping the boundary term (drop 1: $\Psi\to0$ at infinity while $\partial_t\Psi$ stays bounded),
$$\frac{d\langle p\rangle}{dt} = -i\hbar\int\big(\partial_t\Psi^*\,\partial_x\Psi - \partial_x\Psi^*\,\partial_t\Psi\big)dx.$$
Substitute $\partial_t\Psi = \frac{i\hbar}{2m}\Psi'' - \frac{i}{\hbar}V\Psi$ and $\partial_t\Psi^* = -\frac{i\hbar}{2m}\Psi^{*\prime\prime} + \frac{i}{\hbar}V\Psi^*$ ($V$ real). The integrand becomes
$$-\frac{i\hbar}{2m}\big(\Psi^{*\prime\prime}\Psi' + \Psi^{*\prime}\Psi''\big) + \frac{i}{\hbar}V\big(\Psi^*\Psi' + \Psi^{*\prime}\Psi\big) = -\frac{i\hbar}{2m}\partial_x\big(\Psi^{*\prime}\Psi'\big) + \frac{i}{\hbar}V\,\partial_x|\Psi|^2.$$
The first piece integrates to $-\frac{i\hbar}{2m}\big[|\Psi'|^2\big]_{-\infty}^\infty = 0$ (drop 2: $\Psi'\to0$). The second, by parts: $\frac{i}{\hbar}\int V\partial_x|\Psi|^2dx = \frac{i}{\hbar}\big[V|\Psi|^2\big]_{-\infty}^\infty - \frac{i}{\hbar}\int V'|\Psi|^2dx$, and the boundary term vanishes (drop 3: $V|\Psi|^2\to0$, true when $\Psi$ decays faster than $V$ grows, e.g. Gaussian tails against a polynomial $V$). So
$$\frac{d\langle p\rangle}{dt} = -i\hbar\cdot\Big(-\frac{i}{\hbar}\Big)\int V'|\Psi|^2dx.$$
Do the constants carefully: $-i\hbar\times\big(-\frac{i}{\hbar}\big) = i^2 = -1$. So $\frac{d\langle p\rangle}{dt} = -\int V'|\Psi|^2dx = \langle -\partial_xV\rangle$. ∎ The three dropped terms are the price of this route; the commutator route of §3.4 hides the same drops inside "$\hat H$ is Hermitian."

**B2.** (a) With $\int_{-\infty}^\infty e^{-\lambda|x|}e^{-ipx/\hbar}dx = \int_0^\infty e^{-(\lambda + ip/\hbar)x}dx + \int_0^\infty e^{-(\lambda - ip/\hbar)x}dx = \frac{1}{\lambda + ip/\hbar} + \frac{1}{\lambda - ip/\hbar} = \frac{2\lambda}{\lambda^2 + p^2/\hbar^2}$:
$$\Phi(p) = \frac{\sqrt\lambda}{\sqrt{2\pi\hbar}}\cdot\frac{2\lambda\hbar^2}{\hbar^2\lambda^2 + p^2} = \sqrt{\frac2\pi}\,\frac{(\hbar\lambda)^{3/2}}{\hbar^2\lambda^2 + p^2}.$$
A Lorentzian. Normalization: substitute $p = \hbar\lambda u$: $\int|\Phi|^2dp = \frac{2}{\pi}(\hbar\lambda)^3\cdot\frac{\hbar\lambda}{(\hbar\lambda)^4}\int\frac{du}{(1+u^2)^2} = \frac2\pi\cdot\frac\pi2 = 1$, using $\int_{-\infty}^\infty(1+u^2)^{-2}du = \pi/2$. The project integrates it numerically: $1.00000000$.

(b) $x$-space: $\Psi' = -\lambda\,\mathrm{sgn}(x)\Psi$ (away from $0$), so $|\Psi'|^2 = \lambda^2|\Psi|^2$ and $\langle p^2\rangle = \hbar^2\lambda^2$. $p$-space: $\int p^2|\Phi|^2dp = \frac2\pi(\hbar\lambda)^3\cdot\frac{(\hbar\lambda)^3}{(\hbar\lambda)^4}\int\frac{u^2du}{(1+u^2)^2} = \frac2\pi\hbar^2\lambda^2\cdot\frac\pi2 = \hbar^2\lambda^2$, using $\int u^2(1+u^2)^{-2}du = \int(1+u^2)^{-1}du - \int(1+u^2)^{-2}du = \pi - \frac\pi2$. They agree, as (3.18) says. $\langle p\rangle = 0$ ($\Psi$ real).

(c) $\sigma_x\sigma_p = \frac{1}{\sqrt2\lambda}\cdot\hbar\lambda = \frac{\hbar}{\sqrt2} = 0.7071\hbar$, $41\%$ above $\hbar/2$. Through (3.14), $\sigma_p = 2\pi\hbar\sigma_\xi$, so F3's $\sigma_x\sigma_\xi = \sqrt2/4\pi$ becomes $\sigma_x\sigma_p = 2\pi\hbar\cdot\sqrt2/4\pi = \hbar/\sqrt2$. Same number, as it must be.

(d) $|\Phi|^2\sim p^{-4}$ for large $|p|$. Then $p^2|\Phi|^2\sim p^{-2}$ is integrable, but $p^4|\Phi|^2\to$ const is not, so $\langle p^4\rangle = \infty$. In $x$-space: $\Psi'$ jumps from $+\lambda^{3/2}$ to $-\lambda^{3/2}$ at the kink, so $\Psi''$ contains $-2\lambda^{3/2}\delta(x)$, and $\langle p^4\rangle = \hbar^4\|\Psi''\|^2$ would need $\int\delta^2$. A kink costs the fourth moment; a jump in $\Psi$ itself would cost the second (F3's rect function).

**B3.** (a) $\Psi' = ik(Ae^{ikx} - Be^{-ikx})$, so
$\Psi^*\Psi' = ik\big(A^*e^{-ikx} + B^*e^{ikx}\big)\big(Ae^{ikx} - Be^{-ikx}\big) = ik\big(|A|^2 - |B|^2 + (z^* - z)\big)$ with $z = A^*Be^{-2ikx}$.
$z^* - z = -2i\operatorname{Im}z$ is purely imaginary, so $ik(z^* - z) = 2k\operatorname{Im}z$ is real and contributes nothing to $\operatorname{Im}(\Psi^*\Psi')$. Hence $\operatorname{Im}(\Psi^*\Psi') = k(|A|^2 - |B|^2)$ and $J = \frac{\hbar k}{m}(|A|^2 - |B|^2)$. The density $|\Psi|^2$ does contain the interference term $2\operatorname{Re}z$ (standing-wave ripples), but the current does not.

(b) For a stationary state $\partial_t\rho = 0$, so (3.5) gives $\partial_xJ = 0$: $J$ is the same constant on both sides. Left: $J = \frac{\hbar k}{m}(|A|^2 - |B|^2)$, the incident current minus the reflected one. Right: $J = \frac{\hbar q}{m}|C|^2$ by the plane-wave result. Define $R$ and $T$ as reflected and transmitted current divided by incident current: $R = \frac{(\hbar k/m)|B|^2}{(\hbar k/m)|A|^2} = \frac{|B|^2}{|A|^2}$, $T = \frac{(\hbar q/m)|C|^2}{(\hbar k/m)|A|^2} = \frac qk\frac{|C|^2}{|A|^2}$. Equating the two expressions for $J$ and dividing by the incident current gives $1 - R = T$. The factor $q/k$ is there because probability per unit length moves at different speeds on the two sides. QM-3 computes $B$ and $C$ from continuity of $\Psi$ and $\Psi'$; for $k = 2$, $V_0 = 1.5$ ($\hbar = m = 1$) it gives $R = 1/9$, the number the project's step uses.

**B4.** (a) Using $\sin\alpha\sin\beta = \frac12[\cos(\alpha-\beta) - \cos(\alpha+\beta)]$:
$\int_0^Lx\sin\frac{\pi x}{L}\sin\frac{2\pi x}{L}dx = \frac12\int_0^Lx\Big(\cos\frac{\pi x}{L} - \cos\frac{3\pi x}{L}\Big)dx.$
By parts, $\int_0^Lx\cos\frac{n\pi x}{L}dx = \frac{L^2}{n^2\pi^2}\big((-1)^n - 1\big)$, which is $-\frac{2L^2}{\pi^2}$ for $n = 1$ and $-\frac{2L^2}{9\pi^2}$ for $n = 3$. So the integral is $\frac12\big(-\frac{2L^2}{\pi^2} + \frac{2L^2}{9\pi^2}\big) = -\frac{8L^2}{9\pi^2}$, and $x_{12} = \frac2L\cdot\big(-\frac{8L^2}{9\pi^2}\big) = -\frac{16L}{9\pi^2}$.
With real $\psi_n$, $\langle x\rangle = \frac12x_{11} + \frac12x_{22} + x_{12}\cos\omega_{21}t$ (same algebra as (3.31)). By symmetry of $\psi_n^2$ about $L/2$, $x_{11} = x_{22} = L/2$. So $\langle x\rangle = \frac L2 - \frac{16L}{9\pi^2}\cos\omega_{21}t$.

(b) Route 1: $m\frac{d\langle x\rangle}{dt} = m\frac{16L}{9\pi^2}\omega_{21}\sin\omega_{21}t = m\frac{16L}{9\pi^2}\cdot\frac{3\pi^2\hbar}{2mL^2}\sin\omega_{21}t = \frac{8\hbar}{3L}\sin\omega_{21}t$.
Route 2: $p_{12} = -i\hbar\frac2L\cdot\frac{2\pi}{L}\int_0^L\sin\frac{\pi x}{L}\cos\frac{2\pi x}{L}dx$. With $\sin\alpha\cos\beta = \frac12[\sin(\alpha+\beta) + \sin(\alpha - \beta)]$, the integral is $\frac12\int_0^L\big(\sin\frac{3\pi x}{L} - \sin\frac{\pi x}{L}\big)dx = \frac12\big(\frac{2L}{3\pi} - \frac{2L}{\pi}\big) = -\frac{2L}{3\pi}$. So $p_{12} = -i\hbar\frac{4\pi}{L^2}\big(-\frac{2L}{3\pi}\big) = \frac{8i\hbar}{3L}$, and $p_{21} = \overline{p_{12}}$ (Hermiticity) $= -\frac{8i\hbar}{3L}$; $p_{11} = p_{22} = 0$ (real states). Then
$\langle p\rangle = \frac12\big(p_{12}e^{-i\omega_{21}t} + p_{21}e^{i\omega_{21}t}\big) = \frac12\cdot\frac{8i\hbar}{3L}\big(e^{-i\omega_{21}t} - e^{i\omega_{21}t}\big) = \frac12\cdot\frac{8i\hbar}{3L}(-2i\sin\omega_{21}t) = \frac{8\hbar}{3L}\sin\omega_{21}t.$
The two routes agree, which is Theorem 3.2 in action. (Here the $t$-dependence of the cross term is $e^{i(E_1 - E_2)t/\hbar} = e^{-i\omega_{21}t}$ multiplying $\overline{c_1}c_2p_{12}$.)

(c) $\omega_{21} = 3\pi^2/2 = 14.80$, period $2\pi/\omega_{21} = 4/3\pi = 0.4244$; $\langle x\rangle$ amplitude $16/9\pi^2 = 0.1801$; $\langle p\rangle$ amplitude $8/3 = 2.667$. Part 10 of the project checks $\langle x\rangle$ at $t = 0$, half a period and a full period by direct integration.

**B5.** (a) $|\Psi|^2$ does not see the chirp, so it is the Gaussian of (3.16): $\langle x\rangle = 0$, $\sigma_x = a$. $\Psi' = \big(-\frac{x}{2a^2} + 2icx\big)\Psi$. $\langle p\rangle = -i\hbar\int\big(-\frac{x}{2a^2} + 2icx\big)|\Psi|^2dx = 0$ (odd integrand). $|\Psi'|^2 = x^2\big(\frac{1}{4a^4} + 4c^2\big)|\Psi|^2$, so $\langle p^2\rangle = \hbar^2\big(\frac{1}{4a^4} + 4c^2\big)a^2 = \hbar^2\big(\frac{1}{4a^2} + 4c^2a^2\big)$. Then $\sigma_x\sigma_p = \hbar a\sqrt{\frac{1}{4a^2} + 4c^2a^2} = \frac\hbar2\sqrt{1 + 16c^2a^4}$. Equality with $\hbar/2$ iff $c = 0$: a Gaussian modulus is not enough.

(b) With $k_0 = 0$, the exponent of (3.25) is $-\frac{x^2}{4a^2(1+i\beta)} = -\frac{x^2(1 - i\beta)}{4a^2(1+\beta^2)} = -\frac{x^2}{4A^2} + iCx^2$ with $A^2 = a^2(1+\beta^2)$ and $C = \frac{\beta}{4a^2(1+\beta^2)}$. The prefactor $(2\pi a^2)^{-1/4}(1+i\beta)^{-1/2}$ has modulus $(2\pi A^2)^{-1/4}$ times a constant phase. So it is the chirped Gaussian of (a) with width $A$ and chirp $C$. Then $16C^2A^4 = 16\frac{\beta^2}{16a^4(1+\beta^2)^2}a^4(1+\beta^2)^2 = \beta^2$ and (a) gives $\frac\hbar2\sqrt{1+\beta^2}$, matching §3.8.

(c) From (b), the evolved state at "time" $\beta$ has $C = \beta/(4A^2)$. Given $A = 1$ and $C = -\frac14$: $\beta = 4CA^2 = -1$. Then $a^2 = A^2/(1+\beta^2) = \frac12$ and $t = 2ma^2\beta/\hbar = -1$. So $\Psi(x,0)$ is the minimum-uncertainty Gaussian with $a = 1/\sqrt2$ evolved to $t = -1$. Evolving forward, the width shrinks to its minimum $a = 0.7071$ at $t = +1$, where $\sigma_x\sigma_p = \hbar/2$, and grows afterwards. Check: $\sigma_p = \sqrt{\frac14 + 4\cdot\frac1{16}} = 0.7071 = \hbar/2a_{\min}$, constant in time. The project finds the minimum width $0.70711$ at $t = 1.000$. Physically, a negative chirp means the phase gradient (local momentum $\hbar S_x\propto 2\hbar Cx$) points inward on both sides: the packet is focused.

**B6.** (a) Write $x = x_0 + y$. For the Gaussian, $\langle y\rangle = \langle y^3\rangle = 0$ and $\langle y^2\rangle = a^2$. Then $\langle x^3\rangle = x_0^3 + 3x_0^2\langle y\rangle + 3x_0\langle y^2\rangle + \langle y^3\rangle = x_0^3 + 3x_0a^2$, so $\langle -V'\rangle = -\langle x^3\rangle = -(8 + 6a^2)$ for $x_0 = 2$, against $-V'(\langle x\rangle) = -8$.

(b) Let $\Delta(t) = \langle x\rangle - x_{\rm cl}(t)$. At $t = 0$, $\Delta = 0$ and $\dot\Delta = \langle p\rangle/m - 0 = 0$ ($\Psi$ real). By Ehrenfest, $m\ddot\Delta = \langle -V'\rangle + V'(x_{\rm cl}) = -6a^2$ at $t = 0$. So $\Delta\approx -3a^2t^2$, and $|\Delta|$ reaches $0.1$ at $t\approx\sqrt{0.1/3a^2} = 0.183/a$: $0.61$, $0.26$, $0.18$ for $a = 0.3$, $0.707$, $1$.

(c) The code in Part 11:
```python
for aw_ in (0.3, 1 / np.sqrt(2), 1.0):
    r_ = run_ehrenfest(Vq, dVq, gaussian(gh.x, 2.0, 0.0, aw_), dte, 2000)
    t_dev[aw_] = r_["t"][np.argmax(np.abs(r_["x"] - cl.y[0]) > 0.1)]
```
It measures $0.56$, $0.30$, $0.20$, close to (b). Narrower is better here, but not indefinitely: a packet of width $a$ has $\sigma_p = \hbar/2a$ and (roughly, as for a free packet) $\sigma_x(t)^2\approx a^2 + \hbar^2t^2/4m^2a^2$. The deviation grows like $\sigma_x(t)^2$, not $a^2$, and $a^2 + \hbar^2t^2/4m^2a^2$ is smallest at $a^2 = \hbar t/2m$. For any target time there is an optimal width; a very narrow packet spreads so fast that it leaves the classical orbit sooner. That is the uncertainty principle setting a limit on how classical a quantum particle can look in an anharmonic well.

**B7.** (a) (3.11) used only Hermiticity of $\hat H$, so it holds for any $\hat Q$, Hermitian or not. With $[\hat H, \hat x] = -\frac{i\hbar}{m}\hat p$ and $[\hat H, \hat p] = i\hbar V'$ (Theorem 3.5) and the product rule $[\hat H, \hat A\hat B] = [\hat H, \hat A]\hat B + \hat A[\hat H, \hat B]$:
$[\hat H, \hat x\hat p] = -\frac{i\hbar}{m}\hat p^2 + i\hbar\,xV'.$
So $\frac{d}{dt}\langle xp\rangle = \frac i\hbar\big(-\frac{i\hbar}{m}\langle p^2\rangle + i\hbar\langle xV'\rangle\big) = \frac{\langle p^2\rangle}{m} - \langle xV'\rangle = 2\langle T\rangle - \langle xV'\rangle$.

(b) In a stationary state every expectation value is constant (Theorem 3.9(a)), so the left side is $0$: $2\langle T\rangle = \langle xV'\rangle$. For $V = \alpha|x|^n$, $xV' = x\cdot n\alpha|x|^{n-1}\mathrm{sgn}(x) = n\alpha|x|^n = nV$, so $2\langle T\rangle = n\langle V\rangle$. The oscillator ($n = 2$): $\langle T\rangle = \langle V\rangle$, the quantum version of CM-8 Exercise A3's result.

(c) Example 4.6: $\langle T\rangle = \frac14$, $\langle xV'\rangle = \langle x^2\rangle = x_0^2 + \sigma_x^2 = 4.5$. So $\frac{d}{dt}\langle xp\rangle = \frac12 - 4.5 = -4$. Classically $xp = (2\cos t)(-2\sin t) = -2\sin2t$, whose derivative at $0$ is $-4$. They agree. ($\langle xp\rangle$ itself carries a constant imaginary part $i\hbar/2$ by A5(b), which has zero time derivative.)

**B8.** (a) Let $A = -iV\,dt/\hbar$, $B = -iT\,dt/\hbar$, both $O(dt)$. BCH: $e^Ae^B = \exp\big(A + B + \frac12[A,B] + O(dt^3)\big)$. Since the exact step is $e^{A+B}$, the unsymmetric split errs by $\frac12[A,B] = -\frac{dt^2}{2\hbar^2}[V, T]$ per step: $O(dt^2)$ locally, $O(dt)$ after $t/dt$ steps. For the symmetric split, first $e^{A/2}e^B = \exp\big(\frac A2 + B + \frac14[A,B] + O(dt^3)\big)$. Multiply by $e^{A/2}$ on the right and apply BCH again with $X = \frac A2 + B + \frac14[A,B]$, $Y = \frac A2$: $\exp\big(X + Y + \frac12[X, Y] + O(dt^3)\big)$, and $\frac12[X,Y] = \frac12[\frac A2 + B, \frac A2] + O(dt^3) = \frac14[B, A] + O(dt^3) = -\frac14[A,B]$. The $\frac14[A,B]$ terms cancel: $e^{A/2}e^Be^{A/2} = e^{A+B+O(dt^3)}$. Local error $O(dt^3)$, global $O(dt^2)$.

(b) The project's errors (norm of the state difference at $t = 2$) are $1.26\times10^{-3}$, $3.12\times10^{-4}$, $7.42\times10^{-5}$, with ratios $4.05$ and $4.20$. Halving $dt$ divides the error by $4 = 2^2$: second order, as (a) predicts.

(c) $p_{\max} = \pi\hbar/dx = 67.0$ (in units $\hbar = 1$). Each factor of (7.1) is a multiplication by unit-modulus phases in some basis, hence unitary for *any* $dt$: the norm is conserved to rounding even if $dt$ is absurdly large. Accuracy is separate. The splitting error is set by the commutators acting on the actual state, roughly by $p_{\rm typ}^2dt/m\hbar$ and $V_{\rm typ}dt/\hbar$ where the state has weight. A conserved norm is stability, not accuracy (LA-13's distinction); check energy conservation and halve $dt$.

## C. Challenge

**C1.** (a) Divide by $i\hbar$: $\partial_t\Phi + F\partial_p\Phi = -\frac{i}{\hbar}\frac{p^2}{2m}\Phi$. Along characteristics $p(t) = p_s + Ft$ (the classical momentum under force $F$), $\frac{d}{dt}\Phi(p(t),t) = \partial_t\Phi + F\partial_p\Phi = -\frac{i}{2m\hbar}(p_s + Ft)^2\Phi$. Integrating from $0$ to $t$:
$\Phi(p_s + Ft, t) = \Phi_0(p_s)\exp\Big(-\frac{i}{2m\hbar}\int_0^t(p_s + Ft')^2dt'\Big) = \Phi_0(p_s)\exp\Big(-\frac{i}{6m\hbar F}\big[(p_s+Ft)^3 - p_s^3\big]\Big).$
Set $p = p_s + Ft$: $\Phi(p,t) = \Phi_0(p - Ft)\exp\big(-\frac{i}{6m\hbar F}[p^3 - (p - Ft)^3]\big)$. Check: at $t = 0$ it is $\Phi_0(p)$; as $F\to0$, $p^3 - (p - Ft)^3\approx3p^2Ft$ and the phase becomes $-ip^2t/2m\hbar$, the free result.

(b) $|\Phi(p,t)|^2 = |\Phi_0(p - Ft)|^2$: the distribution slides rigidly by $Ft$. Hence $\langle p\rangle = p_0 + Ft$ and $\sigma_p$ is constant. For $\sigma_x$: let $C = \frac12\langle xp + px\rangle - \langle x\rangle\langle p\rangle$. By (3.11) with $V' = -F$: $\frac{d\langle x^2\rangle}{dt} = \frac1m\langle xp + px\rangle$, so $\frac{d\sigma_x^2}{dt} = \frac{2C}{m}$; and $[\hat H, \hat x\hat p + \hat p\hat x] = -\frac{2i\hbar}{m}\hat p^2 - 2i\hbar F\hat x$, so $\frac{d}{dt}\frac12\langle xp + px\rangle = \frac{\langle p^2\rangle}{m} + F\langle x\rangle$, and $\frac{dC}{dt} = \frac{\langle p^2\rangle}{m} + F\langle x\rangle - \frac{\langle p\rangle^2}{m} - \langle x\rangle F = \frac{\sigma_p^2}{m}$. These equations ($\dot{\sigma_x^2} = 2C/m$, $\dot C = \sigma_p^2/m$, $\sigma_p$ constant) are identical to the free ones ($F = 0$), with the same initial data. So $\sigma_x(t)$ is exactly the free width. "Your turn" 1 in the project asserts this.

(c) The project builds $\Phi(k, 6)$ from (a) for $a = 1$ ($\phi_0 = (2/\pi)^{1/4}e^{-k^2}$, $F = 0.1$), transforms to $x$, and finds $|\langle\Psi_{\rm analytic}, \Psi_{\rm numeric}\rangle| = 1.000000000000$. Exactness: the BCH series for $e^{A/2}e^Be^{A/2}$ contains only nested commutators of $T\propto p^2$ and $V\propto x$. $[p^2, x]\propto p$; then $[p^2, p] = 0$ and $[x, p] = i\hbar$, a number. Every higher commutator is zero or a number, so the splitting error is a multiple of the identity: a global phase per step. The overlap's modulus ignores it.

**C2.** (a) With $y = x - vt$ and $f = e^{-y^2/4\sigma^2}$: $f' = -\frac{y}{2\sigma^2}f$, $f'' = \big(\frac{y^2}{4\sigma^4} - \frac{1}{2\sigma^2}\big)f$. The normalization constant cancels in $f''/f$, so
$Q = -\frac{\hbar^2}{2m}\frac{f''}{f} = \frac{\hbar^2}{2m}\Big(\frac{1}{2\sigma^2} - \frac{y^2}{4\sigma^4}\Big), \qquad -\partial_xQ = \frac{\hbar^2y}{4m\sigma^4}.$
$Q$ is an inverted parabola centered on the packet; the force is proportional to $y$ and points away from the center: outward.

(b) $S/\hbar$ is the imaginary part of the exponent in (3.25) plus $x$-independent terms. From B5(b), that imaginary part is $k_0x + \frac{\beta y^2}{4a^2(1+\beta^2)}$, so $S_x = \hbar k_0 + \frac{\hbar\beta y}{2a^2(1+\beta^2)}$ and the velocity field is $\frac{S_x}{m} = v + \beta'\frac{\beta y}{1+\beta^2}$ with $\beta' = \frac{d\beta}{dt} = \frac{\hbar}{2ma^2}$. Along $x(t) = vt + y_0\sigma(t)/a$ we have $y = y_0\sqrt{1+\beta^2}$, and
$\dot x = v + y_0\frac{d}{dt}\sqrt{1+\beta^2} = v + y_0\frac{\beta\beta'}{\sqrt{1+\beta^2}} = v + \beta'\frac{\beta y}{1+\beta^2} = \frac{S_x}{m}.$
Second derivative: $\ddot x = y_0\frac{d}{dt}\frac{\beta\beta'}{\sqrt{1+\beta^2}} = y_0\beta'^2\frac{(1+\beta^2) - \beta^2}{(1+\beta^2)^{3/2}} = \frac{y_0\beta'^2}{(1+\beta^2)^{3/2}} = \frac{y_0\hbar^2}{4m^2a^4(1+\beta^2)^{3/2}}$. The force per mass from (a) is $\frac{\hbar^2y}{4m^2\sigma^4} = \frac{\hbar^2y_0\sqrt{1+\beta^2}}{4m^2a^4(1+\beta^2)^2}$, the same. So $m\ddot x = -\partial_xQ$. The project checks both identities in SymPy.

(c) In the Madelung picture probability is a fluid moving with velocity $S_x/m$ and pushed by $-\partial_x(V + Q)$. With $V = 0$, the quantum potential of a localized blob is a hill at its center, so it pushes the fluid outward and the blob expands; the narrower the blob, the larger $Q$ ($\propto1/\sigma^2$), which is the spreading rate $\tau\propto a^2$ in another language.

**C3.** (a) $\langle x\rangle = L/2$ by symmetry. With $\sin^2\theta = \frac12(1 - \cos2\theta)$: $\langle x^2\rangle = \frac1L\int_0^Lx^2dx - \frac1L\int_0^Lx^2\cos\frac{2\pi x}{L}dx = \frac{L^2}{3} - \frac{L^2}{2\pi^2}$, using $\int_0^Lx^2\cos kx\,dx = \frac{2L}{k^2}$ for $kL = 2\pi$ (two integrations by parts; the sine terms vanish at both ends). So $\sigma_x^2 = L^2\big(\frac13 - \frac1{2\pi^2} - \frac14\big) = L^2\big(\frac1{12} - \frac1{2\pi^2}\big)$. $\langle p\rangle = 0$ ($\psi_1$ real). $\langle p^2\rangle = \hbar^2\int|\psi_1'|^2 = \hbar^2\frac{\pi^2}{L^2}\cdot\frac2L\int_0^L\cos^2\frac{\pi x}{L}dx = \frac{\pi^2\hbar^2}{L^2}$. So $\sigma_x\sigma_p = \pi\hbar\sqrt{\frac1{12} - \frac1{2\pi^2}} = \hbar\sqrt{\frac{\pi^2}{12} - \frac12} = 0.5679\hbar$.

(b) No. By Theorem 3.7, $\sigma_x\sigma_p = \hbar/2$ only for $Ce^{-(x-x_0)^2/4a^2 + ip_0x/\hbar}$, which is nonzero for every $x$. A state confined to $0 < x < L$ vanishes outside, so it is not of that form, and the inequality is strict. The bound is still the infimum: a Gaussian of width $a\ll L$ centered at $L/2$, smoothly cut off far in its tails, has $\sigma_x\sigma_p$ as close to $\hbar/2$ as you like. (Its energy, $\approx\hbar^2/8ma^2$, is then far above $E_1$; the ground state trades uncertainty for energy.)
