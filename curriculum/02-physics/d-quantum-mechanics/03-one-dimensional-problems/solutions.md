# QM-3 Solutions

Every number below is printed or asserted in `project.py` (Part 8 and the parts cited).

## A. Warm-up

**A1.** $E_1 = \pi^2\hbar^2/2m_ea^2$. With $\hbar = 1.0546\times10^{-34}$ J s, $m_e = 9.109\times10^{-31}$ kg, $a = 10^{-9}$ m: $E_1 = 6.025\times10^{-20}$ J $= 0.3760$ eV. $E_2 = 4E_1$, so $E_2 - E_1 = 3E_1 = 1.1281$ eV. The photon wavelength is $\lambda = hc/\Delta E = 1239.84\text{ eV nm}/1.1281\text{ eV} = 1099$ nm: near infrared, just beyond the visible. (This is why quantum dots a few nm across emit in the visible: $E\propto 1/a^2$, and the effective mass in a semiconductor is smaller than $m_e$, which raises $E$ further.)

**A2.** $\langle\psi_1|\psi_3\rangle = \frac2a\int_0^a\sin\frac{\pi x}{a}\sin\frac{3\pi x}{a}dx$. Product to sum: $\sin A\sin B = \frac12[\cos(A - B) - \cos(A + B)]$, so the integrand is $\frac1a[\cos\frac{2\pi x}{a} - \cos\frac{4\pi x}{a}]$. Each cosine integrates over $[0, a]$ to $\frac{a}{2\pi}\sin2\pi = 0$ and $\frac{a}{4\pi}\sin4\pi = 0$. So $\langle\psi_1|\psi_3\rangle = 0$.

$\langle\psi_2|\psi_2\rangle = \frac2a\int_0^a\sin^2\frac{2\pi x}{a}dx = \frac2a\int_0^a\frac12\big(1 - \cos\frac{4\pi x}{a}\big)dx = \frac2a\cdot\frac a2 = 1$.

**A3.** As an operator, $a_+ = \frac{1}{\sqrt{2\hbar m\omega}}\big(-\hbar\frac{d}{dx} + m\omega x\big)$ (since $-i\hat p = -\hbar\,d/dx$). With $\psi_0 = Ne^{-m\omega x^2/2\hbar}$, $N = (m\omega/\pi\hbar)^{1/4}$: $\psi_0' = -\frac{m\omega x}{\hbar}\psi_0$, so
$$a_+\psi_0 = \frac{1}{\sqrt{2\hbar m\omega}}(m\omega x + m\omega x)\psi_0 = \frac{2m\omega}{\sqrt{2\hbar m\omega}}x\psi_0 = \sqrt{\frac{2m\omega}{\hbar}}\,x\,\psi_0.$$
Normalization: $\int x^2\psi_0^2dx = \sigma_x^2 = \hbar/2m\omega$ (QM-2 (3.16) with $a^2 = \hbar/2m\omega$), so $\int|a_+\psi_0|^2 = \frac{2m\omega}{\hbar}\cdot\frac{\hbar}{2m\omega} = 1$, consistent with $\sqrt{0 + 1} = 1$ in (3.13).

Now $a_- = \frac{1}{\sqrt{2\hbar m\omega}}\big(\hbar\frac{d}{dx} + m\omega x\big)$. Let $\psi_1 = cx\psi_0$, $c = \sqrt{2m\omega/\hbar}$. Then $\psi_1' = c\psi_0 + cx\psi_0' = c\psi_0 - c\frac{m\omega x^2}{\hbar}\psi_0$, and
$$a_-\psi_1 = \frac{c}{\sqrt{2\hbar m\omega}}\Big(\hbar\psi_0 - m\omega x^2\psi_0 + m\omega x^2\psi_0\Big) = \frac{c\hbar}{\sqrt{2\hbar m\omega}}\psi_0 = \sqrt{\frac{2m\omega}{\hbar}}\cdot\sqrt{\frac{\hbar}{2m\omega}}\,\psi_0 = \psi_0.$$
Coefficient $1 = \sqrt1$, as (3.13) says.

**A4.** $k' = \sqrt{2m(E - V_0)}/\hbar = \sqrt{2m\cdot\frac13V_0}/\hbar$ and $k = \sqrt{2m\cdot\frac43V_0}/\hbar$, so $k'/k = \sqrt{1/4} = \frac12$. With $r = \frac12$ in (3.23):
$R = \big(\frac{1 - 1/2}{1 + 1/2}\big)^2 = \big(\frac13\big)^2 = \frac19 = 0.1111$, and $T = \frac{4\cdot\frac12}{(3/2)^2} = \frac{2}{9/4} = \frac89 = 0.8889$. $R + T = 1$. Naive $|F/A|^2 = \big(\frac{2k}{k + k'}\big)^2 = \big(\frac{2}{3/2}\big)^2 = \frac{16}{9} = 1.7778$: impossible as a probability. The project asserts all three.

**A5.** $|\psi|^2 = \kappa e^{-2\kappa|x|}$. $P(|x| > 1/\kappa) = 2\int_{1/\kappa}^\infty\kappa e^{-2\kappa x}dx = 2\kappa\cdot\frac{e^{-2}}{2\kappa} = e^{-2} = 0.1353$. The sketch: a symmetric tent of two exponentials, peak $\sqrt\kappa$ at $x = 0$, with a corner there. $\psi'$ jumps from $+\kappa\sqrt\kappa$ to $-\kappa\sqrt\kappa$, a jump of $-2\kappa\sqrt\kappa = -\frac{2m\alpha}{\hbar^2}\psi(0)$, (3.19). $\psi'$ is continuous everywhere else.

## B. Core

**B1.** (a) By (3.14), $\hat x\psi_n = \sqrt{\hbar/2m\omega}\,(a_+\psi_n + a_-\psi_n) = \sqrt{\hbar/2m\omega}\,(\sqrt{n+1}\,\psi_{n+1} + \sqrt n\,\psi_{n-1})$ by (3.13). Take the inner product with $\psi_m$ and use orthonormality: $\langle\psi_m|\hat x|\psi_n\rangle = \sqrt{\hbar/2m\omega}\,(\sqrt{n+1}\,\delta_{m,n+1} + \sqrt n\,\delta_{m,n-1})$.

(b) An oscillating electric field couples to the dipole $\propto x$, and the transition rate is $\propto|\langle\psi_m|\hat x|\psi_n\rangle|^2$ (QM-8's golden rule). Only $m = n\pm1$ is nonzero, and every such transition has $|\Delta E| = \hbar\omega$, so absorption and emission occur at the single frequency $\omega$, whatever $n$ is.

(c) $\hat x^4 = (\hbar/2m\omega)^2(a_+ + a_-)^4$. Expanding the fourth power gives 16 words in $a_\pm$. Each $a_+$ raises $n$ by one and each $a_-$ lowers it, so to return to $\psi_0$ a word needs two of each (6 words). A word whose rightmost letter is $a_-$ gives $0$ (since $a_-\psi_0 = 0$); a word whose leftmost letter is $a_+$ ends in a state with $n \ge 1$, orthogonal to $\psi_0$. That leaves two words:
- $a_-a_-a_+a_+\psi_0$: $a_+a_+\psi_0 = \sqrt1\sqrt2\,\psi_2$; then $a_-a_-\psi_2 = \sqrt2\sqrt1\,\psi_0$. Contribution $2$.
- $a_-a_+a_-a_+\psi_0$: $a_+\psi_0 = \psi_1$, $a_-\psi_1 = \psi_0$, again $\psi_1$, again $\psi_0$. Contribution $1$.

Total $\langle\psi_0|\hat x^4|\psi_0\rangle = 3(\hbar/2m\omega)^2$. Check: $|\psi_0|^2$ is a Gaussian with variance $\sigma^2 = \hbar/2m\omega$, whose fourth moment is $3\sigma^4$. Agreement. (The project computes $(X^4)_{00} = 0.75$ with $\hbar = m = \omega = 1$.)

**B2.** (a) $\int_0^{a/2}\frac2a\,dx = 1$.

(b) $c_n = \int_0^{a/2}\sqrt{\frac2a}\sin\frac{n\pi x}{a}\sqrt{\frac2a}\,dx = \frac2a\cdot\frac{a}{n\pi}\Big[-\cos\frac{n\pi x}{a}\Big]_0^{a/2} = \frac{2}{n\pi}\Big(1 - \cos\frac{n\pi}{2}\Big)$.
So $c_1 = \frac2\pi = 0.6366$, $c_2 = \frac{2}{2\pi}(1 + 1) = \frac2\pi = 0.6366$, $c_3 = \frac{2}{3\pi} = 0.2122$, $c_4 = 0$.

(c) $P(E_1) = P(E_2) = 4/\pi^2 = 0.4053$.

(d) $|c_n|^2E_n = \frac{4(1 - \cos\frac{n\pi}{2})^2}{n^2\pi^2}\cdot\frac{n^2\pi^2\hbar^2}{2ma^2} = \frac{2\hbar^2}{ma^2}\Big(1 - \cos\frac{n\pi}{2}\Big)^2$. The factor $(1 - \cos\frac{n\pi}{2})^2$ cycles through $1, 4, 1, 0$, average $\frac32$. The terms do not decay, so the partial sums grow linearly, $\approx\frac{3\hbar^2}{ma^2}N$ (project with $\hbar = m = a = 1$: 34.0, 300.0, 3000.0 for $N = 10, 100, 1000$). $\langle H\rangle = \infty$.

Physically: kinetic energy is $\frac{\hbar^2}{2m}\int|\Psi'|^2dx$, and $\Psi$ jumps from $\sqrt{2/a}$ to $0$ at $x = a/2$, so $\Psi'$ contains a delta function, whose square is not integrable. Each truncation $\sum_{n\le N}c_n\psi_n$ is smooth (it rounds the jump over a width $\sim a/N$), with $\int|\Psi_N'|^2 \sim N$ finite, but the bound grows without limit. A real preparation never produces an exact jump; it produces a steep but smooth edge, and the energy is then large but finite.

**B3.** With $c_0 = c_1 = 1/\sqrt2$, $\Psi(t) = \frac{1}{\sqrt2}(\psi_0e^{-i\omega t/2} + \psi_1e^{-3i\omega t/2})$. By B1(a), the only nonzero matrix elements of $\hat x$ between $\psi_0, \psi_1$ are $\langle\psi_0|\hat x|\psi_1\rangle = \langle\psi_1|\hat x|\psi_0\rangle = \sqrt{\hbar/2m\omega}$. So
$$\langle x\rangle = \frac12\sqrt{\frac{\hbar}{2m\omega}}\big(e^{i\omega t/2}e^{-3i\omega t/2} + e^{3i\omega t/2}e^{-i\omega t/2}\big) = \sqrt{\frac{\hbar}{2m\omega}}\cos\omega t.$$
For $\hat p = i\sqrt{\hbar m\omega/2}(a_+ - a_-)$: $\langle\psi_1|\hat p|\psi_0\rangle = i\sqrt{\hbar m\omega/2}$, $\langle\psi_0|\hat p|\psi_1\rangle = -i\sqrt{\hbar m\omega/2}$. Then
The two cross terms of $\langle\Psi|\hat p|\Psi\rangle$ carry the phases $e^{3i\omega t/2}e^{-i\omega t/2} = e^{i\omega t}$ (for $\langle\psi_1|\hat p|\psi_0\rangle$) and $e^{-i\omega t}$ (for $\langle\psi_0|\hat p|\psi_1\rangle$), each with weight $\frac12$:
$$\langle p\rangle = \frac12\sqrt{\frac{\hbar m\omega}{2}}\big(ie^{i\omega t} - ie^{-i\omega t}\big) = \frac12\sqrt{\frac{\hbar m\omega}{2}}\,(i)(2i\sin\omega t) = -\sqrt{\frac{\hbar m\omega}{2}}\sin\omega t.$$ Check: $m\,d\langle x\rangle/dt = -m\omega\sqrt{\hbar/2m\omega}\sin\omega t = -\sqrt{\hbar m\omega/2}\sin\omega t = \langle p\rangle$, and $d\langle p\rangle/dt = -\omega\sqrt{\hbar m\omega/2}\cos\omega t = -m\omega^2\sqrt{\hbar/2m\omega}\cos\omega t = -m\omega^2\langle x\rangle$. Ehrenfest holds exactly (the force is linear). $\langle H\rangle = \frac12\cdot\frac12\hbar\omega + \frac12\cdot\frac32\hbar\omega = \hbar\omega$. The project asserts $\langle x\rangle(t) = \cos t/\sqrt2$ in units $\hbar = m = \omega = 1$.

**B4.** (a) Steps 1–2 hold on the ring: $W = \psi_1'\psi_2 - \psi_2'\psi_1$ is constant. Step 3 fails: there is no infinity where $\psi$ vanishes, so nothing forces $W = 0$. For $\psi_1 = e^{2\pi inx/L}$, $\psi_2 = e^{-2\pi inx/L}$, with $q = 2\pi n/L$: $W = iq\,e^{iqx}e^{-iqx} - (-iq)e^{-iqx}e^{iqx} = 2iq \ne 0$. Two independent solutions coexist at the same energy.

(b) Same failure: $e^{\pm ikx}$ do not vanish at infinity, $W = 2ik \ne 0$. Theorem 3.2 is about bound states only.

(c) Step 3 concludes $\chi = c\psi$ from Theorem 3.2. Without nondegeneracy, $\chi = \psi(-x)$ is a solution but need not be a multiple of $\psi$. Example: $V = 0$ is even, and $\psi = e^{ikx}$ is an eigenfunction, but $\psi(-x) = e^{-ikx}$ is not a multiple of it; $e^{ikx}$ is neither even nor odd. Step 2 still applies: $\cos kx$ and $\sin kx$ are definite-parity eigenfunctions spanning the same eigenspace. "Can be chosen with definite parity" survives; "is automatically of definite parity" does not.

**B5.** (a) $z_0 = a\sqrt{2m_eV_0}/\hbar$ with $a = 0.5$ nm, $V_0 = 10$ eV: $\sqrt{2m_eV_0} = \sqrt{2\cdot9.109\times10^{-31}\cdot1.602\times10^{-18}} = 1.708\times10^{-24}$ kg m/s; dividing by $\hbar$ gives $1.620\times10^{10}$ m$^{-1}$; times $a$, $z_0 = 8.100$. Number of bound states $\lfloor2z_0/\pi\rfloor + 1 = \lfloor5.157\rfloor + 1 = 6$.

(b) The infinite well of width $2a = 1$ nm, measured from the bottom, has $E_1 + V_0 \approx 0.3760$ eV (A1), so $E_1 \approx -9.6240$ eV. The true value (project, root of (3.22)) is $z_1 = 1.3974$, $E_1 = V_0(z_1^2/z_0^2 - 1) = -9.7024$ eV: **lower**. The true wavefunction leaks into the walls, so it lives in an effectively wider box and has less kinetic energy.

(c) The even condition (3.22) is $\kappa a = z\tan z$ (multiply both sides by $z$; $\kappa a = \sqrt{z_0^2 - z^2}$). For $z_0 \ll 1$ the root has $z < z_0 \ll 1$, so $\tan z \approx z$ and $\kappa a \approx z^2$. Since $z^2 = z_0^2 - (\kappa a)^2$ and $\kappa a$ is small (it is $\approx z^2 \le z_0^2$), $z^2 = z_0^2 - O(z_0^4)$, so $\kappa a \approx z_0^2$. Then $|E| = \hbar^2\kappa^2/2m \approx \frac{\hbar^2z_0^4}{2ma^2}$. Since $z_0^2 = 2mV_0a^2/\hbar^2$, this is $|E| \approx \frac{2mV_0^2a^2}{\hbar^2} = \frac{m}{2\hbar^2}(2aV_0)^2$, which is the delta-well binding $m\alpha^2/2\hbar^2$ with $\alpha = 2aV_0$, the well's "area." A shallow narrow well binds like a delta of the same strength.

**B6.** (a) In §3.8, $M_{22} = \cosh\kappa a + i\frac{\kappa^2 - k^2}{2k\kappa}\sinh\kappa a$ came from the propagator (3.24) with $q = i\kappa$. For $E > V_0$ the propagator has real $q = l$: $p_{11} = p_{22} = \cos la$, $p_{12} = \sin(la)/l$, $p_{21} = -l\sin la$. The same formula for $M_{22}$ (step 2 of §3.8) gives $M_{22} = \cos la - \frac{ik}{2l}\sin la - \frac{il}{2k}\sin la = \cos la - i\frac{k^2 + l^2}{2kl}\sin la$. Then $|M_{22}|^2 = \cos^2la + \frac{(k^2 + l^2)^2}{4k^2l^2}\sin^2la = 1 + \frac{(k^2 - l^2)^2}{4k^2l^2}\sin^2la$, using $(k^2 + l^2)^2 - 4k^2l^2 = (k^2 - l^2)^2$. Now $k^2 - l^2 = 2mV_0/\hbar^2$ and $k^2l^2 = 4m^2E(E - V_0)/\hbar^4$, so the coefficient is $V_0^2/4E(E - V_0)$ and
$$T = \Big[1 + \frac{V_0^2\sin^2(la)}{4E(E - V_0)}\Big]^{-1}.$$
(The substitution $\kappa\to il$ in (3.26) gives the same: $\sinh(i la) = i\sin la$, so $\sinh^2 \to -\sin^2$, and $V_0 - E \to -(E - V_0)$; the two minus signs cancel.) The project asserts this against `transfer()` for $1 < E \le 6$.

(b) $T = 1$ when $\sin la = 0$, $la = n\pi$: $E = V_0 + \frac{n^2\pi^2\hbar^2}{2ma^2}$, $n = 1, 2, \dots$

(c) As $E\to\infty$, $V_0^2/4E(E - V_0)\to0$ and $\sin^2\le1$, so $T\to1$. As $E\to V_0^+$, $l\to0$ and $\sin^2(la)\approx l^2a^2 = 2m(E - V_0)a^2/\hbar^2$; the $(E - V_0)$ cancels: $T\to\big[1 + \frac{V_0^2\cdot2ma^2}{4V_0\hbar^2}\big]^{-1} = \big[1 + \frac{mV_0a^2}{2\hbar^2}\big]^{-1}$. From below, $\sinh^2(\kappa a)\approx\kappa^2a^2 = 2m(V_0 - E)a^2/\hbar^2$ in (3.26), and the same cancellation gives the same limit. $T(E)$ is continuous through $E = V_0$, as it must be: nothing special happens to the wave there.

**B7.** (a) Far left $\psi = Ae^{ik_0x} + Be^{-ik_0x}$, far right $\psi = Fe^{ik_Nx}$ (with $k_0, k_N$ real; if $k_N$ is imaginary the right-hand solution is real up to a constant and carries $J = 0$, so $T = 0$). $J$ is independent of $x$ because $\partial_t\rho = 0$ in a stationary state and $\partial_xJ = -\partial_t\rho$ (QM-2 (3.5)), and that needs only $V$ real (used in deriving (3.5)). Evaluate on the left: $\psi^*\psi' = ik_0(|A|^2 - |B|^2) + ik_0(\bar BAe^{2ik_0x} - \bar AB e^{-2ik_0x})$. The second bracket is $z - \bar z$ for $z = \bar BAe^{2ik_0x}$, which is purely imaginary, so after multiplying by $ik_0$ it is real and drops out of $\operatorname{Im}$. So $J_L = \frac{\hbar k_0}{m}(|A|^2 - |B|^2)$. On the right, $J_R = \frac{\hbar k_N}{m}|F|^2$. Setting $J_L = J_R$ and dividing by $\hbar k_0|A|^2/m$: $1 - R = T$.

(b) From (3.21): $R + T = \frac{\beta^2}{1 + \beta^2} + \frac{1}{1 + \beta^2} = 1$.

(c) $M = W_{k_N}^{-1}\Phi W_{k_0}$ with $\det\Phi = 1$ (each $P_q$ has $\det = \cos^2 + \sin^2 = 1$; each $D$ is lower triangular with unit diagonal) and $\det W_k = -2ik$, so $\det M = (-2ik_0)/(-2ik_N) = k_0/k_N$. To connect with (a): because $V$ is real, if $\psi$ is a solution so is $\bar\psi$. Conjugation maps the amplitudes $(P, Q)$ of $Pe^{ikx} + Qe^{-ikx}$ to $(\bar Q, \bar P)$ (for real $k$). Applying $M$ to both $\psi$ and $\bar\psi$ forces $M_{11} = \bar M_{22}$ and $M_{12} = \bar M_{21}$. Then $\det M = |M_{22}|^2 - |M_{21}|^2 = k_0/k_N$. Divide by $|M_{22}|^2$: $1 - |M_{21}/M_{22}|^2 = (k_0/k_N)/|M_{22}|^2$, which by (3.31) is $1 - R = T$. The determinant identity is current conservation.

**B8.** Code (add below Part 5 of `project.py`):

```python
Er = np.linspace(1.01, 5, 400)                  # E / V0 with V0 = 1
Tst, Rst, Mst = transfer(Er, [0], [0, 1])
k, kp = np.sqrt(2 * Er), np.sqrt(2 * (Er - 1))
assert np.allclose(Rst, ((k - kp) / (k + kp))**2) and np.allclose(Tst, 4 * k * kp / (k + kp)**2)
naive = (2 * k / (k + kp))**2
fig, ax = plt.subplots()
ax.plot(Er, Rst, label="R"); ax.plot(Er, Tst, label="T"); ax.plot(Er, naive, "--", label="|F/A|^2")
ax.axhline(1, c="k", lw=0.5); ax.legend(); fig.savefig(FIG / "step.png"); plt.close(fig)
print(transfer([1.0], [0], [0, -3])[1])         # step down: R = 1/9
```

(a) The curves match (3.23): $R$ falls from 1 at $E = V_0$ (where $k' = 0$) toward 0, and $T$ rises toward 1. (b) $|F/A|^2 = 4/(1 + k'/k)^2 > 1$ for **every** $E > V_0$, because $k' < k$; it tends to 4 at threshold. It is not a probability. (c) For $V = 0\to-3E$: $k' = \sqrt{2m(E + 3E)}/\hbar = 2k$, so $R = \big(\frac{1 - 2}{1 + 2}\big)^2 = \frac19$ (the project asserts $0.1111$). Reflection is caused by the abrupt change of wavenumber, not by the sign of the change: matching $\psi$ and $\psi'$ cannot be done with a transmitted wave alone when $k' \ne k$ (step 2 of §3.7 forces $B \propto k - k'$). Classically nothing reflects off a drop.

**B9.** (a) Let $(\Pi f)(x) = f(-x)$. Then $(\Pi\hat x\Pi f)(x) = \Pi[x f(-x)] = (-x)f(x)$, so $\Pi\hat x\Pi = -\hat x$, and $\Pi\frac{d}{dx}\Pi = -\frac{d}{dx}$ (chain rule). So $\Pi a_+\Pi = -a_+$, i.e. $\Pi a_+ = -a_+\Pi$. $\psi_0$ is even: $\Pi\psi_0 = \psi_0$. Then $\Pi\psi_n = \frac{1}{\sqrt{n!}}\Pi a_+^n\psi_0 = \frac{(-1)^n}{\sqrt{n!}}a_+^n\Pi\psi_0 = (-1)^n\psi_n$.

(b) By Theorem 3.4, a bound eigenstate of an even $V$ is even or odd, so $|\psi(x)|^2$ is even. Then $x|\psi|^2$ is odd and integrates to $0$.

(c) No: B3 gives $\langle x\rangle = \sqrt{\hbar/2m\omega}\cos\omega t$. The superposition is not an eigenstate; it mixes an even and an odd state, and $|\Psi|^2$ contains the cross term $\psi_0\psi_1\cos\omega t$, which is odd in $x$. (b) applies only to eigenstates.

## C. Challenge

**C1.** (a) $\hbar = m = 1$; the jump condition (3.19) at each delta reads $\Delta\psi' = -2\alpha\psi$. Energy $E = -\kappa^2/2$. *Even:* $\psi = A\cosh\kappa x$ for $|x| < b$, $\psi = Ce^{-\kappa|x|}$ for $|x| > b$. At $x = b$: continuity $A\cosh\kappa b = Ce^{-\kappa b}$; jump $-\kappa Ce^{-\kappa b} - \kappa A\sinh\kappa b = -2\alpha A\cosh\kappa b$. Substitute $Ce^{-\kappa b} = A\cosh\kappa b$ and divide by $A\cosh\kappa b$: $\kappa + \kappa\tanh\kappa b = 2\alpha$. *Odd:* $\psi = A\sinh\kappa x$ inside, $\pm Ce^{-\kappa|x|}$ outside. At $x = b$: $A\sinh\kappa b = Ce^{-\kappa b}$ and $-\kappa Ce^{-\kappa b} - \kappa A\cosh\kappa b = -2\alpha A\sinh\kappa b$; dividing by $A\sinh\kappa b$: $\kappa + \kappa\coth\kappa b = 2\alpha$.

(b) Let $g(\kappa) = \kappa(1 + \coth\kappa b)$. As $\kappa\to0^+$, $\kappa\coth\kappa b\to1/b$, so $g\to1/b$. $g$ is increasing ($\kappa\coth\kappa b$ is increasing in $\kappa$: its derivative is $\coth u - u/\sinh^2u$ with $u = \kappa b$, and $\sinh u\cosh u > u$) and $g\to\infty$. So $g(\kappa) = 2\alpha$ has a root iff $2\alpha > 1/b$, i.e. $\alpha b > \frac12$. (The even equation's left side starts at $0$, so the even state always exists, the two-delta version of Theorem 3.1.)

(c) For $\kappa b\gg1$: $\tanh u = 1 - 2e^{-2u} + \dots$ and $\coth u = 1 + 2e^{-2u} + \dots$ Even: $\kappa(2 - 2e^{-2\kappa b}) = 2\alpha$, so $\kappa_+ = \alpha(1 + e^{-2\kappa b} + \dots) \approx \alpha + \alpha e^{-2\alpha b}$ (to leading order replace $\kappa$ by $\alpha$ in the exponent). Odd: $\kappa_- \approx \alpha - \alpha e^{-2\alpha b}$. The splitting $\Delta = E_- - E_+ = \frac12(\kappa_+^2 - \kappa_-^2) = \frac12(\kappa_+ + \kappa_-)(\kappa_+ - \kappa_-) \approx \frac12(2\alpha)(2\alpha e^{-2\alpha b}) = 2\alpha^2e^{-2\alpha b}$. With $|E_b| = \alpha^2/2$, $\kappa_0 = \alpha$, $d = 2b$: $\Delta = 4|E_b|e^{-\kappa_0d}$. The exponential in the separation is the double-well law (3.32).

(d) Exact (project, Brent's method on the two equations) vs $2e^{-2b}$:

| $b$ | $E_{\rm even}$ | $E_{\rm odd}$ | $\Delta$ exact | $2e^{-2b}$ |
|---|---|---|---|---|
| 1 | −0.614783 | −0.317455 | 0.2973 | 0.2707 |
| 2 | −0.517251 | −0.480369 | 0.03688 | 0.03663 |
| 3 | −0.502446 | −0.497487 | 0.004959 | 0.004958 |

The asymptotic formula is 10% off at $b = 1$ and exact to four figures by $b = 3$. Both levels move symmetrically away from the single-delta $-\frac12$.

**C2.** (a) For $x > 0$: $-\frac{\hbar^2}{2m}\psi'' + Fx\psi = E\psi$, i.e. $\psi'' = \frac{2mF}{\hbar^2}\big(x - \frac EF\big)\psi$. Put $x - E/F = \ell y$, so $\epsilon = E/F\ell$ in the exercise's notation, and $d^2/dx^2 = \ell^{-2}d^2/dy^2$: $\frac{d^2\psi}{dy^2} = \frac{2mF\ell^3}{\hbar^2}y\psi = y\psi$ when $\ell^3 = \hbar^2/2mF$. Airy's equation has solutions $\operatorname{Ai}(y)$, decaying as $y\to+\infty$, and $\operatorname{Bi}(y)$, growing; normalizability keeps $\psi = \operatorname{Ai}(x/\ell - E/F\ell)$. The wall requires $\psi(0) = \operatorname{Ai}(-E/F\ell) = 0$, so $-E/F\ell = a_n$, where $a_1 > a_2 > \dots$ are the (negative) zeros of $\operatorname{Ai}$. Hence $E_n = -a_nF\ell = -a_n\big(\frac{\hbar^2F^2}{2m}\big)^{1/3}$.

(b) Code (added in Part 8 of `project.py`):

```python
an = ai_zeros(3)[0]                        # -2.33811, -4.08795, -5.52056
Etri = -an * 0.5**(1 / 3)                  # hbar = m = F = 1: (F^2/2)^(1/3)
xt = np.linspace(0, 20, 8002)[1:-1]
assert np.allclose(fd_eigs(xt, xt, 3)[0], Etri, rtol=1e-5)
```

$E_1 = 2.33811\times0.79370 = 1.8558$, $E_2 = 3.2446$, $E_3 = 4.3817$, and `fd_eigs` on $[0, 20]$ agrees to $10^{-5}$. The spacing shrinks with $n$ ($E_n \sim n^{2/3}$ for large $n$), between the box ($n^2$) and the oscillator ($n$): the softer the walls, the more slowly the levels grow.

**C3.** (a) Write $V - E = E(r_2/r - 1)$, since $2Zk_ee^2 = Er_2$. So $\gamma = \frac{\sqrt{2mE}}{\hbar}\int_{r_1}^{r_2}\sqrt{\frac{r_2}{r} - 1}\,dr$. Substitute $r = r_2\sin^2u$, $dr = 2r_2\sin u\cos u\,du$, $\sqrt{r_2/r - 1} = \sqrt{1/\sin^2u - 1} = \cos u/\sin u$ (for $u\in(0, \pi/2]$). The limits are $u_1 = \arcsin\sqrt{r_1/r_2}$ and $\pi/2$:
$$\int_{r_1}^{r_2}\sqrt{\frac{r_2}{r} - 1}\,dr = 2r_2\int_{u_1}^{\pi/2}\cos^2u\,du = r_2\Big[u + \sin u\cos u\Big]_{u_1}^{\pi/2} = r_2\Big[\frac\pi2 - u_1 - \sin u_1\cos u_1\Big].$$
With $\sin u_1 = \sqrt{r_1/r_2}$, $\cos u_1 = \sqrt{1 - r_1/r_2}$ and $\frac\pi2 - u_1 = \arccos\sqrt{r_1/r_2}$, this is (3.30). The project checks it against numerical quadrature (agreement $2\times10^{-5}$, limited by the square-root endpoint).

(b) For $\epsilon = \sqrt{r_1/r_2}\ll1$: $u_1 = \arcsin\epsilon = \epsilon + O(\epsilon^3)$ and $\sin u_1\cos u_1 = \epsilon\sqrt{1 - \epsilon^2} = \epsilon + O(\epsilon^3)$. So the bracket is $\frac\pi2 - 2\epsilon + O(\epsilon^3)$ and $\gamma \approx \frac{\sqrt{2mE}}{\hbar}\big(\frac\pi2r_2 - 2\sqrt{r_1r_2}\big)$. The first term: $\sqrt{2mE}\,r_2 = \sqrt{2mE}\cdot\frac{2Zk_ee^2}{E} = \frac{2\sqrt{2m}\,Zk_ee^2}{\sqrt E}$, so $\gamma_{\rm lead} = \frac{\pi\sqrt{2m}\,Zk_ee^2}{\hbar\sqrt E}\propto\frac{Z}{\sqrt E}$.

(c) The half-life is $\propto e^{2\gamma}$ and $\gamma\propto Z/\sqrt E$ is large (tens), so a modest change in $\sqrt E$ shifts $2\gamma$ by tens of units. For U-238 versus Po-212, $\gamma$ drops from $44.4$ to $16.3$ and $e^{2(44.4 - 16.3)} = e^{56.1} = 10^{24.4}$, essentially the whole observed range ($10^{23.7}$) (the $2r_1/v$ prefactor changes by less than a factor of 2).
