# F3 Solutions

## A. Warm-up

**A1.** $\hat f(\xi) = \int_0^\infty e^{-3x}e^{-2\pi i x\xi}dx = \int_0^\infty e^{-(3+2\pi i\xi)x}dx = \left[\dfrac{e^{-(3+2\pi i\xi)x}}{-(3+2\pi i\xi)}\right]_0^\infty = \dfrac{1}{3+2\pi i\xi}$, the upper limit vanishing because $|e^{-(3+2\pi i\xi)x}| = e^{-3x}\to0$.

Check: $\hat f(0) = 1/3$ and $\int_0^\infty e^{-3x}dx = 1/3$. ✓

At $\xi = 1/(2\pi)$: $\hat f = 1/(3+i)$, $|\hat f| = 1/\sqrt{9+1} = 1/\sqrt{10} = 0.3162$.

**A2.** Start from $\hat\Pi = \operatorname{sinc}$ (3.3). Stretch (3.7) with $a = 1/4$: $\widehat{\Pi(x/4)}(\xi) = \frac{1}{|1/4|}\operatorname{sinc}(\xi/(1/4)) = 4\operatorname{sinc}(4\xi)$. Shift (3.6) by $2$: $\hat g(\xi) = e^{-2\pi i\cdot2\xi}\cdot4\operatorname{sinc}(4\xi) = 4e^{-4\pi i\xi}\operatorname{sinc}(4\xi)$.

Check: $\hat g(0) = 4$ and $\int g = $ width $= 4$. ✓ First zero of $|\hat g| = 4|\operatorname{sinc}(4\xi)|$: $\sin(4\pi\xi) = 0$ at $\xi = 1/4$. (A box of width $4$ has its first null at $1/4$ cycle per unit: first null $= 1/\text{width}$, always.)

**A3.** (a) The text's convention is the angular asymmetric one: $F(\omega) = \hat f(\omega/2\pi)$, i.e. $\hat f(\xi) = F(2\pi\xi) = \dfrac{2}{1+4\pi^2\xi^2}$. Equation (3.5) with $a=1$ gives $\dfrac{2\cdot1}{1+4\pi^2\xi^2}$. ✓

(b) Symmetric convention: $F_{\text{sym}}(k) = \frac{1}{\sqrt{2\pi}}\hat f(k/2\pi) = \frac{1}{\sqrt{2\pi}}\cdot\frac{2}{1+k^2}$.

(c) At zero frequency: ours $\hat f(0) = \int f = 2$; angular asymmetric $F(0) = \int f = 2$ (same, since the $1/2\pi$ is on the inverse); symmetric $F_{\text{sym}}(0) = \frac{1}{\sqrt{2\pi}}\int f = 2/\sqrt{2\pi} = 0.7979$. The zero-frequency check is the fastest way to detect which convention a formula is in.

**A4.** $N = 4$, $\omega = -i$, $F[k] = \sum_{n=0}^3f[n](-i)^{kn}$.

$f = (1,-1,1,-1)$: $F[k] = 1 - (-i)^k + (-i)^{2k} - (-i)^{3k} = 1-(-i)^k+(-1)^k-(-i)^{3k}$.
$k=0$: $1-1+1-1 = 0$. $k=1$: $1-(-i)+(-1)-(i) = 1+i-1-i = 0$. $k=2$: $1-(-1)+1-(-1) = 4$. $k=3$: $(-i)^3 = i$, $(-i)^9 = (-i)^1 = -i$: $1-i-1+i = 0$.
$F = (0,0,4,0)$. The input is $\cos(2\pi\cdot2n/4) = (-1)^n$: the tone with $k=2$ cycles per record (the Nyquist tone). Its coefficient is $N = 4$ because $e_2 = e_{-2}$ coincide, so the two halves $\tfrac12(e_2+e_{-2})$ add up to one full spike.

$g = (0,1,0,-1)$: $F[k] = (-i)^k - (-i)^{3k}$.
$k=0$: $0$. $k=1$: $-i - i = -2i$. $k=2$: $1-1 = 0$. $k=3$: $i-(-i) = 2i$.
$F = (0,-2i,0,2i)$. The input is $\sin(2\pi n/4)$, one cycle per record. Since $\sin\theta = (e^{i\theta}-e^{-i\theta})/(2i)$, $g = \frac{1}{2i}(e_1-e_3)$, giving $F[1] = N/(2i) = -2i$ and $F[3] = -N/(2i) = 2i$. Purely imaginary and odd: the real-odd rule of §3.4(g), discrete form.

**A5.** Folding formula $\nu_{\text{alias}} = |\nu - f_s\operatorname{round}(\nu/f_s)|$.
(a) $60$ at $50$: $\operatorname{round}(1.2) = 1$, $|60-50| = 10$ Hz.
(b) $60$ at $100$: $\operatorname{round}(0.6) = 1$, $|60-100| = 40$ Hz.
(c) $60$ at $130$: $\operatorname{round}(0.46) = 0$, $60$ Hz: no aliasing, since $130>120 = 2\cdot60$.
(d) $13$ rev/s at $12$ fps: $13-12 = +1$: appears to turn forward at $1$ rev/s (between frames it turns $1\tfrac1{12}$ turns, indistinguishable from $\tfrac1{12}$ of a turn forward). Compare $7$ at $10$: $7-10 = -3$, backward.
(e) $5$ Hz at $10$ Hz sits exactly at Nyquist. Samples $\sin(2\pi\cdot5\cdot n/10) = \sin(\pi n) = 0$ for every $n$. The signal is invisible. (With a phase shift, $\sin(\pi n+\phi) = (-1)^n\sin\phi$: amplitude depends on the unknown phase; see B8.)

**A6.** (a) $\langle x\delta',\varphi\rangle = \langle\delta', x\varphi\rangle = -\frac{d}{dx}(x\varphi)\big|_{x=0} = -[\varphi(0)+x\varphi'(x)]_{x=0} = -\varphi(0) = \langle-\delta,\varphi\rangle$. So $x\delta' = -\delta$. (Consistent with differentiating $x\delta = 0$ by the product rule: $\delta + x\delta' = 0$.)

(b) Take $a = 2\pi$ in (3.5): $\widehat{e^{-2\pi|x|}}(\xi) = \dfrac{4\pi}{4\pi^2+4\pi^2\xi^2} = \dfrac{1}{\pi(1+\xi^2)}$. So $g(x) := \dfrac{1}{\pi(1+x^2)}$ is the transform of $f(x) = e^{-2\pi|x|}$. Duality (3.10): $\hat g(\xi) = \hat{\hat f}(\xi) = f(-\xi) = e^{-2\pi|\xi|}$. By linearity, $\widehat{1/(1+x^2)} = \pi\hat g = \pi e^{-2\pi|\xi|}$.

Check at $0$: $\pi e^0 = \pi$, and $\int\frac{dx}{1+x^2} = \pi$. ✓ (The Cauchy density $\frac{1}{\pi(1+x^2)}$ has characteristic function $e^{-|t|}$ in the angular convention: the transform of a Lorentzian is a two-sided exponential, which has a kink at $0$, which is why the Cauchy distribution has no mean. F4 returns to this.)

## B. Core

**B1.** (a) $\widehat{e^{2\pi i\nu x}f(x)}(\xi) = \int f(x)e^{2\pi i\nu x}e^{-2\pi i x\xi}dx = \int f(x)e^{-2\pi i x(\xi-\nu)}dx = \hat f(\xi-\nu)$. ∎ With $\cos(2\pi\nu x) = \frac12(e^{2\pi i\nu x}+e^{-2\pi i\nu x})$ and linearity: $\widehat{\cos(2\pi\nu x)f} = \frac12[\hat f(\xi-\nu)+\hat f(\xi+\nu)]$.

(b) $\widehat{\Pi(x/\tau)} = \tau\operatorname{sinc}(\tau\xi)$ (Example 4.1). So
$$\hat f(\xi) = \frac\tau2\Big[\operatorname{sinc}\big(\tau(\xi-\nu_0)\big)+\operatorname{sinc}\big(\tau(\xi+\nu_0)\big)\Big].$$

(c) Two sinc lobes centered at $\pm\nu_0$, each of peak height $\tau/2$ (at $\nu_0\tau = 10$ the cross-talk from the other lobe at $\xi = \nu_0$ is $\frac\tau2\operatorname{sinc}(20) = 0$, so the peak is exactly $\tau/2$). Each main lobe extends from its center to the first zeros at $\xi = \nu_0\pm1/\tau$: width $2/\tau$. Sidelobes decay like $1/|\xi-\nu_0|$, the first one at about $22\%$ of the peak ($-13$ dB). Longer observation $\tau$: taller, narrower lobes. The total area under each lobe is $\frac12$ (= $f$'s "amount" at that frequency), independent of $\tau$.

(d) With the record of length $\tau$ sampled and DFT'd, the frequency grid is $\xi_k = k/\tau$ (3.28), and by Theorem 3.20 the DFT reads (up to the factor $\Delta$) $\hat f$ at those points: $F[k]\propto\frac\tau2[\operatorname{sinc}(k-\nu_0\tau)+\operatorname{sinc}(k+\nu_0\tau)]$. If $\nu_0\tau = m\in\mathbb Z$, then $\operatorname{sinc}(k-m) = \delta_{km}$ and $\operatorname{sinc}(k+m) = \delta_{k,-m}$: exactly two nonzero bins, $k = \pm m$ (i.e. $m$ and $N-m$). Conversely, if two bins are nonzero and the rest vanish, $\operatorname{sinc}(k-\nu_0\tau) = 0$ for all but one $k$, forcing $\nu_0\tau$ to be that integer. (Exactly, the DFT replaces each sinc by its periodic version, the Dirichlet kernel $\sin(\pi x)/(N\sin(\pi x/N))$ up to a phase factor, because the DFT samples the periodized spectrum (3.33). It has the same zeros, at the integers $x\not\equiv0\pmod N$, so the argument is unchanged.) If $\nu_0\tau\notin\mathbb Z$, the zeros of the sinc (at $\xi = \nu_0+j/\tau$) fall *between* grid points and no bin is exactly zero: energy leaks everywhere, largest in the two bins adjacent to $\nu_0\tau$, decaying like $1/|k-\nu_0\tau|$. That is the $(1, 1.192, 2.798, 2.398,\dots)$ of Example 4.5, where $\nu_0\tau = 2.5$.

**B2.** (a) $\|f\|^2 = a\int e^{-2a|x|}dx = 2a\int_0^\infty e^{-2ax}dx = 2a\cdot\frac{1}{2a} = 1$. ✓

(b) $\sigma_x^2 = a\int x^2e^{-2a|x|}dx = 2a\int_0^\infty x^2e^{-2ax}dx = 2a\cdot\frac{2!}{(2a)^3} = \frac{4a}{8a^3} = \frac{1}{2a^2}$. So $\sigma_x = \frac{1}{a\sqrt2}$.

(c) *Directly.* By (3.5), $\hat f(\xi) = \sqrt a\cdot\dfrac{2a}{a^2+4\pi^2\xi^2}$, so $|\hat f|^2 = \dfrac{4a^3}{(a^2+4\pi^2\xi^2)^2}$ and $\sigma_\xi^2 = 4a^3\displaystyle\int\frac{\xi^2\,d\xi}{(a^2+4\pi^2\xi^2)^2}$. Substitute $\xi = \frac{a}{2\pi}t$, $d\xi = \frac{a}{2\pi}dt$, $\xi^2 = \frac{a^2t^2}{4\pi^2}$, $a^2+4\pi^2\xi^2 = a^2(1+t^2)$:
$$\sigma_\xi^2 = 4a^3\cdot\frac{a^2}{4\pi^2}\cdot\frac{a}{2\pi}\cdot\frac{1}{a^4}\int\frac{t^2}{(1+t^2)^2}dt = \frac{a^2}{2\pi^3}\int\frac{t^2\,dt}{(1+t^2)^2}.$$
Partial fractions: $\dfrac{t^2}{(1+t^2)^2} = \dfrac{1}{1+t^2}-\dfrac{1}{(1+t^2)^2}$. $\int\frac{dt}{1+t^2} = \pi$; with $t = \tan\theta$, $\int\frac{dt}{(1+t^2)^2} = \int_{-\pi/2}^{\pi/2}\cos^2\theta\,d\theta = \frac\pi2$. So the integral is $\pi-\frac\pi2 = \frac\pi2$, and $\sigma_\xi^2 = \frac{a^2}{2\pi^3}\cdot\frac\pi2 = \frac{a^2}{4\pi^2}$, $\sigma_\xi = \frac{a}{2\pi}$.

*Via $\|f'\|$.* For $x\ne0$, $f'(x) = -a\operatorname{sgn}(x)f(x)$, so $|f'|^2 = a^2|f|^2$ and $\|f'\|^2 = a^2\|f\|^2 = a^2$ (the kink at $0$ is a single point and does not affect the integral; $f'$ has a jump there but no delta, so $f'\in L^2$). Step 1 of Theorem 3.11: $\sigma_\xi = \|f'\|/(2\pi) = a/(2\pi)$. Same answer, four lines instead of fourteen: that is what the derivative rule plus Plancherel buys you.

(d) $\sigma_x\sigma_\xi = \dfrac{1}{a\sqrt2}\cdot\dfrac{a}{2\pi} = \dfrac{1}{2\sqrt2\,\pi} = \dfrac{\sqrt2}{4\pi} = 0.11254$. The Gaussian minimum is $\frac{1}{4\pi} = 0.07958$; the exponential is $\sqrt2\approx1.414$ times worse. Independent of $a$, as the stretch theorem requires (scaling $x$ by $a$ scales $\sigma_x$ by $1/a$ and $\sigma_\xi$ by $a$).

(e) $\hat\Pi = \operatorname{sinc}$, so $\xi^2|\hat\Pi|^2 = \xi^2\frac{\sin^2\pi\xi}{\pi^2\xi^2} = \frac{\sin^2\pi\xi}{\pi^2}$, which does not decay: its integral over $|\xi|<L$ is $\approx L/\pi^2\to\infty$. Equivalently, $\Pi' = \delta_{-1/2}-\delta_{1/2}$ is not in $L^2$, so $\|f'\| = \infty$ in step 1 of the theorem. The hypothesis "$f'\in L^2$" is exactly what excludes jumps. The theorem's conclusion $\sigma_x\sigma_\xi\ge1/(4\pi)$ is trivially true for the rect ($\infty\ge$ anything), but the rect shows that a sharp edge in one domain costs *infinite* spread in the other, not merely a large one.

**B3.** (a) $C_{jk} = c[(j-k)\bmod4]$:
$$C = \begin{pmatrix}2&-1&0&-1\\-1&2&-1&0\\0&-1&2&-1\\-1&0&-1&2\end{pmatrix}.$$
Row $n$ has $2$ on the diagonal, $-1$ at $k = n-1$ and $k = n+1$ (mod $4$): $(Cx)[n] = 2x[n]-x[n-1]-x[n+1]$. This is $-(x[n+1]-2x[n]+x[n-1])$, the negative of the second difference, which approximates $-h^2u''$ on a grid of spacing $h$, with periodic (ring) boundary conditions.

(b) $\hat c[k] = \sum_mc[m]\omega^{km} = 2-\omega^k-\omega^{3k}$ with $\omega = e^{-2\pi i/4}$. Since $\omega^{3k} = \omega^{-k}$ (as $\omega^{4k}=1$), $\hat c[k] = 2-(\omega^k+\omega^{-k}) = 2-2\cos(2\pi k/4)$. Values: $k=0$: $0$; $k=1$: $2$; $k=2$: $4$; $k=3$: $2$. Eigenvalues $(0,2,4,2)$, all $\ge0$: $C$ is positive semidefinite, as a "minus Laplacian" should be.

(c) $e_1 = (1,i,-1,-i)$ (that is $e^{2\pi i\cdot1\cdot n/4} = i^n$). $Ce_1$: row 0: $2(1)-(i)+0-(-i) = 2$. Row 1: $-(1)+2(i)-(-1)+0 = 2i$. Row 2: $0-(i)+2(-1)-(-i) = -2$. Row 3: $-(1)+0-(-1)+2(-i) = -2i$. So $Ce_1 = (2,2i,-2,-2i) = 2e_1$. ✓ Eigenvalue $2 = \hat c[1]$.

(d) By (3.30), $C = \mathbf F^{-1}\Lambda\mathbf F$ with $\Lambda = \operatorname{diag}(0,2,4,2)$. Write $\mathbf u(0) = \sum_k\alpha_ke_k$ (with $\alpha_k = U[k]/N$, the inverse DFT). Each mode evolves independently: $\mathbf u(t) = \sum_k\alpha_k\,e^{-\kappa\lambda_kt/h^2}\,e_k$ with $\lambda_k = 2-2\cos(2\pi k/4)$. The $k=0$ mode (constant temperature, $e_0 = (1,1,1,1)$) has $\lambda_0 = 0$ and never decays: heat is conserved on the ring. The $k=2$ mode ($e_2 = (1,-1,1,-1)$, alternating hot/cold) has $\lambda_2 = 4$ and decays fastest, at rate $4\kappa/h^2$: the sharpest spatial variation smooths out first. Modes $k=1,3$ (one cycle around the ring) decay at $2\kappa/h^2$. This is the discrete version of F2's heat-equation result that the $n$-th Fourier mode decays like $e^{-n^2\kappa t}$.

(e) With fixed ends (Dirichlet), the matrix is tridiagonal Toeplitz but **not** circulant: the corner entries $C_{0,N-1} = C_{N-1,0} = -1$ are missing. It no longer commutes with the cyclic shift $S$, so Theorem 3.24 does not apply and the complex exponentials are not eigenvectors. Its eigenvectors are $\sin(\pi jk/(N+1))$, which vanish at the (virtual) endpoints, matching the boundary condition; the eigenvalues are $2-2\cos(\pi k/(N+1))$ (LA-8 Example 4.6 verifies this for $N = 3$). Same formula, different allowed wavenumbers: the ring quantizes $2\pi k/N$, the fixed bar quantizes $\pi k/(N+1)$. (This is the sine-series versus full-series distinction from F2.)

**B4.** (a) $M(N) = 2M(N/2)+N/2$. Unroll: $M(N) = 2M(N/2)+\frac N2 = 2[2M(N/4)+\frac N4]+\frac N2 = 4M(N/4)+\frac N2+\frac N2 = \dots = 2^jM(N/2^j)+j\frac N2$. At $j = m = \log_2N$: $M(N) = NM(1)+\frac N2\log_2N = \frac N2\log_2N$.
$A(N) = 2A(N/2)+N = 2^jA(N/2^j)+jN$; at $j = \log_2N$: $A(N) = N\log_2N$.

(b) $N = 2^{20}$: $M = 2^{19}\cdot20 = 10{,}485{,}760\approx1.05\times10^7$. $N^2 = 2^{40} = 1{,}099{,}511{,}627{,}776\approx1.10\times10^{12}$. Ratio $= 2^{40}/(2^{19}\cdot20) = 2^{21}/20 = 104{,}857.6$. About $10^5$: a one-second FFT would take a day directly.

(c) Level 1 (four length-2 DFTs): $4$ butterflies, each with twiddle $\omega_2^0 = 1$: $4$ trivial multiplications. Level 2 (two length-4 DFTs): $2\times2 = 4$ butterflies with twiddles $\omega_4^0 = 1$ and $\omega_4^1 = -i$: all trivial. Level 3: $4$ butterflies with $\omega_8^0 = 1$ (trivial), $\omega_8^1$, $\omega_8^2 = -i$ (trivial), $\omega_8^3$: $2$ non-trivial. Total $12 = \frac82\log_28$ as (3.32) says; non-trivial $2$.

(d) The DFT of the zero-padded sequence is $F[k] = \sum_{n=0}^3f[n]e^{-2\pi i kn/8}$: the same trigonometric polynomial $P(\xi) = \sum_{n=0}^3f[n]e^{-2\pi i n\xi}$ that the length-4 DFT evaluates at $\xi = 0,\frac14,\frac12,\frac34$, now evaluated at $\xi = 0,\frac18,\dots,\frac78$. The four new values are determined by the original four samples (they are values of a degree-3 polynomial in $e^{-2\pi i\xi}$ fixed by its four coefficients), so no information is added. What changes is the frequency grid (3.28): spacing $1/(N\Delta)$ halves. That interpolates the spectrum between the old bins, which helps locate a peak (the project's oscillator section) but does not sharpen it: the lobe width is set by the true record length $4\Delta$, not by $N$.

**B5.** (a) Let $\Phi(x) = \sum_nf(x+nT)$; it has period $T$ and (since $f$ decays rapidly) is smooth, so F2 applies with $\Phi(x) = \sum_kc_ke^{2\pi i kx/T}$, $c_k = \frac1T\int_0^T\Phi(x)e^{-2\pi i kx/T}dx$. Then
$$c_k = \frac1T\sum_n\int_0^Tf(x+nT)e^{-2\pi i kx/T}dx = \frac1T\sum_n\int_{nT}^{(n+1)T}f(y)e^{-2\pi i k(y-nT)/T}dy = \frac1T\sum_n\int_{nT}^{(n+1)T}f(y)e^{-2\pi i ky/T}dy,$$
using $e^{2\pi i kn} = 1$. The intervals $[nT,(n+1)T]$ tile $\mathbb R$, so $c_k = \frac1T\int_{-\infty}^\infty f(y)e^{-2\pi i ky/T}dy = \frac1T\hat f(k/T)$. Substituting back: $\sum_nf(x+nT) = \frac1T\sum_k\hat f(k/T)e^{2\pi i kx/T}$. ∎

(b) $x=0$, $T=1$: $\sum_nf(n) = \sum_k\hat f(k)$, which is (3.20). With $f(x) = e^{-\pi tx^2} = g(\sqrt t\,x)$, the stretch theorem and (3.11) give $\hat f(\xi) = t^{-1/2}e^{-\pi\xi^2/t}$. Hence $\sum_ne^{-\pi n^2t} = t^{-1/2}\sum_ke^{-\pi k^2/t}$. At $t = 1/2$: left $= 1+2e^{-\pi/2}+2e^{-2\pi}+2e^{-9\pi/2}+\cdots = 1+0.415759+0.003735+0.000001 = 1.419495$. Right $= \sqrt2(1+2e^{-2\pi}+2e^{-8\pi}+\cdots) = 1.414214\times(1+0.0037349+0.0000000\ldots) = 1.414214\times1.0037349 = 1.419495$. ✓ (The project computes both to $10^{-13}$ and checks these digits.)

(c) The Fourier series coefficients of the $T$-periodization of $f$ are the samples of $\hat f$ at the harmonics $k/T$, scaled by $1/T$. That is Theorem 3.27 with the roles of the two domains exchanged: periodizing in $x$ (convolving with $\mathrm{III}_T$) samples in $\xi$ (multiplying by $\frac1T\mathrm{III}_{1/T}$), just as sampling in $x$ periodizes in $\xi$. One theorem, read in both directions.

**B6.** (a) $y(t) = (h*x)(t) = \int h(\tau)e^{2\pi i\nu(t-\tau)}d\tau = e^{2\pi i\nu t}\int h(\tau)e^{-2\pi i\nu\tau}d\tau = H(\nu)e^{2\pi i\nu t}$. ∎ For real $h$, the response to the real input $\cos(\omega t) = \operatorname{Re}e^{i\omega t}$ is the real part of the response to $e^{i\omega t}$ (linearity with real coefficients): $\operatorname{Re}[H(\omega)e^{i\omega t}] = \operatorname{Re}[|H|e^{i\arg H}e^{i\omega t}] = |H(\omega)|\cos(\omega t+\arg H(\omega))$. (Here $H(\omega)$ means $\hat h$ at $\xi = \omega/2\pi$, as in (3.23).) "Steady state" because the convolution integral assumes the input has been on forever; the transient from switching on decays like $e^{-ct/2m}$ and is not part of this.

(b) $H(\omega) = 1/(4-\omega^2+0.2i\omega)$.
$\omega=1$: denominator $3+0.2i$, $|H| = 1/\sqrt{9.04} = 0.3326$, $\arg H = -\arctan(0.2/3) = -3.814°$.
$\omega=2$: denominator $0.4i$, $|H| = 2.500$, $\arg H = -90.00°$.
$\omega=3$: denominator $-5+0.6i$, $|H| = 1/\sqrt{25.36} = 0.1986$, $\arg H = -(180°-\arctan(0.6/5)) = -(180°-6.843°) = -173.16°$.
So the steady-state responses are $0.3326\cos(t-0.0666)$, $2.500\cos(2t-\pi/2) = 2.5\sin2t$, and $0.1986\cos(3t-3.0222)$.

(c) Below resonance the spring dominates ($k\gg m\omega^2$): the mass follows the force quasi-statically, in phase. At $\omega_0$ the spring and inertia forces cancel exactly ($k = m\omega_0^2$) and only the damper resists; the damper force is $c\dot x$, which is $90°$ ahead of $x$, so for it to balance the applied force the displacement must lag the force by $90°$: velocity in phase with force, maximum power transfer. Above resonance inertia dominates: almost all of the applied force goes into accelerating the mass, $f\approx m\ddot x$, so the force is in phase with the acceleration. For sinusoidal motion $\ddot x = -\omega^2x$, so the acceleration, and with it the force, is opposite to the displacement: the phase approaches $-180°$.

**B7.** Reference implementation (also in `project.py`, "Your turn" reference cells):

```python
def ifft_via_conjugation(X):
    """Inverse DFT from F^{-1} = conj(F)/N: f = conj(F conj(X))/N (eq. 3.27)."""
    X = np.asarray(X, complex)
    return np.conj(fft_recursive(np.conj(X))) / len(X)

X = np.fft.fft(rng.standard_normal(256) + 1j * rng.standard_normal(256))
assert np.max(np.abs(ifft_via_conjugation(X) - np.fft.ifft(X))) < 1e-10

# (b) Parseval for the DFT: sum |f|^2 = (1/N) sum |F|^2
f = rng.standard_normal(128); F = np.fft.fft(f)
assert np.isclose(np.sum(np.abs(f)**2), np.sum(np.abs(F)**2) / 128)

# (c) frequency-domain low-pass
fs, T = 100.0, 4.0
t = np.arange(0, T, 1 / fs)
x = np.sin(2*np.pi*2*t) + 0.5*np.sin(2*np.pi*5*t) + 0.3*rng.standard_normal(t.size)
X = np.fft.fft(x); xi = np.fft.fftfreq(t.size, d=1/fs)
X[np.abs(xi) > 3.0] = 0.0
y = np.fft.ifft(X).real
```

Why (a) works: $\mathbf F$ is symmetric, so $\mathbf F^* = \overline{\mathbf F}$ and $\mathbf F^{-1} = \overline{\mathbf F}/N$ (Theorem 3.19). Then $\mathbf F^{-1}X = \overline{\mathbf F}X/N = \overline{\mathbf F\bar X}/N$: conjugate the input, forward-transform, conjugate the output, divide by $N$. One FFT routine serves both directions.

(c) After filtering, the DFT bin at $5$ Hz ($k = 20$) is exactly zero and the $2$ Hz bin ($k=8$) is untouched, so $y$ contains the $2$ Hz tone plus the part of the noise below $3$ Hz; the project asserts $|Y[20]| = 0$ and $|Y[8]| \approx |X[8]|$. At the ends of the record: zeroing coefficients is multiplying $X$ by an ideal rect in frequency, which by Theorem 3.22 is **circular** convolution of $x$ with the inverse DFT of that rect (a periodic sinc, the Dirichlet kernel). Circular means the convolution wraps around: the filter at $t\approx0$ averages in samples from $t\approx T$. Since the record's two ends do not join smoothly (the signal is not exactly periodic with period $T$), the output shows ringing near both ends. Fixes: zero-pad the record before filtering, or taper the ends with a window. This is the discrete echo of the fact that an ideal low-pass filter's impulse response (a sinc) is infinitely long and non-causal.

**B8.** (a) $t_n = n/f_s$: $\cos(2\pi\frac{f_s}{2}\frac{n}{f_s}+\phi) = \cos(\pi n+\phi) = \cos(\pi n)\cos\phi - \sin(\pi n)\sin\phi = (-1)^n\cos\phi$. The samples are the alternating sequence scaled by $\cos\phi$. The phase $\phi$ and the amplitude are confounded: an amplitude-$1$ tone with $\phi = 60°$ gives the same samples as an amplitude-$0.5$ tone with $\phi = 0$, and $\phi = 90°$ gives all zeros. The sine component at Nyquist is invisible.

(b) At $x = m\Delta$: $\operatorname{sinc}\!\left(\frac{m\Delta-n\Delta}{\Delta}\right) = \operatorname{sinc}(m-n)$, which is $1$ if $m = n$ and $0$ otherwise (since $\sin(\pi(m-n)) = 0$ for integers $m\ne n$). The sum (3.34) collapses to the single term $f(m\Delta)\cdot1$. ✓ The sincs are the "Lagrange basis" of band-limited interpolation: each is $1$ at its own node and $0$ at all others.

(c) By (3.35), the samples of $e^{2\pi i\nu x}$ coincide with those of $e^{2\pi i(\nu-k/\Delta)x}$ for every integer $k$. Among all these candidate frequencies the sampler can only report one; by convention (and by the DFT's frequency grid (3.28), which covers $(-\frac{1}{2\Delta},\frac1{2\Delta}]$) it is the one with $|\nu-k/\Delta|$ smallest, i.e. $k = \operatorname{round}(\nu\Delta)$. For a real sinusoid, $\pm\nu'$ are indistinguishable (the cosine is even in frequency), so the apparent frequency is $|\nu-\frac1\Delta\operatorname{round}(\nu\Delta)|$. With $\nu = 9$, $1/\Delta = 8$: $\operatorname{round}(9/8) = 1$, $|9-8| = 1$ Hz. A $9$ Hz tone sampled at $8$ Hz looks like $1$ Hz.

(d) Theorem 3.28 as stated requires $\hat f(\xi) = 0$ for $|\xi|\ge B$ and $1/\Delta\ge2B$, i.e. the spectrum must vanish *at* $\pm B$ and beyond. If the spectrum is supported in the open interval $|\xi|<B$ and $1/\Delta = 2B$ exactly, the copies in (3.33) touch at $\pm B$ without overlapping, and the theorem applies. But if there is any spectral mass *at* $\xi = \pm B$ (a tone exactly at Nyquist), the $k=0$ and $k=\pm1$ copies overlap on that point and the proof's step 1 fails. Part (a) shows what goes wrong: $\cos(2\pi Bt+\phi)$ with $1/\Delta = 2B$ has samples $(-1)^n\cos\phi$, from which $\phi$ cannot be recovered. In practice one samples strictly above $2B$ and uses an anti-aliasing filter to guarantee $\hat f = 0$ near Nyquist.

## C. Challenge

**C1.** (a) Let $s_n(x) = \operatorname{sinc}\!\left(\frac{x-n\Delta}{\Delta}\right)$. From step 3 of Theorem 3.28, $\widehat{\operatorname{sinc}(x/\Delta)} = \Delta\Pi(\Delta\xi)$; by the shift rule, $\hat s_n(\xi) = \Delta\Pi(\Delta\xi)e^{-2\pi i n\Delta\xi}$. Plancherel (3.15):
$$\langle s_n, s_m\rangle = \langle\hat s_n,\hat s_m\rangle = \Delta^2\int\Pi(\Delta\xi)^2\,\overline{e^{-2\pi i n\Delta\xi}}\,e^{-2\pi i m\Delta\xi}d\xi = \Delta^2\int_{-1/(2\Delta)}^{1/(2\Delta)}e^{2\pi i(n-m)\Delta\xi}d\xi.$$
If $n = m$ the integrand is $1$ and the integral is $1/\Delta$, giving $\Delta$. If $n\ne m$, substitute $\eta = \Delta\xi$, $d\xi = d\eta/\Delta$: $\Delta^2\cdot\frac1\Delta\int_{-1/2}^{1/2}e^{2\pi i(n-m)\eta}d\eta = \Delta\cdot\frac{\sin(\pi(n-m))}{\pi(n-m)} = 0$. So $\langle s_n,s_m\rangle = \Delta\delta_{nm}$. ∎ (The project checks $\langle s_0,s_1\rangle\approx0$ and $\|s_0\|^2\approx\Delta$ numerically.)

(b) By (3.34), $f = \sum_nf(n\Delta)s_n$ with the $s_n$ orthogonal. Then $\|f\|^2 = \sum_{n,m}\overline{f(n\Delta)}f(m\Delta)\langle s_n,s_m\rangle = \Delta\sum_n|f(n\Delta)|^2$. ∎ The energy of a band-limited signal is $\Delta$ times the energy of its samples: a sampling-theorem Parseval, and the continuous justification for computing signal power from samples.

(c) $\sin(\pi x/\Delta) = \sin(2\pi\cdot\frac{1}{2\Delta}x)$, so by (3.19) its transform is $\frac{1}{2i}(\delta_{1/(2\Delta)}-\delta_{-1/(2\Delta)})$, supported exactly at $|\xi| = \frac{1}{2\Delta} = B$. Its samples are $\sin(\pi n) = 0$ for all $n$: identical to the samples of the zero function, which is also band-limited to $B$. So two different functions "band-limited to $B$" share the same samples, and (3.34) returns the wrong one (zero). This does not contradict Theorem 3.28 because the theorem requires $\hat f(\xi) = 0$ for $|\xi|\ge B$, and this $\hat f$ is nonzero *at* $|\xi| = B$. It shows the inequality $1/\Delta\ge2B$ cannot be relaxed: the Nyquist rate is sharp.

**C2.** (a) $(\mathbf F^2)_{kn} = \sum_m\mathbf F_{km}\mathbf F_{mn} = \sum_{m=0}^{N-1}\omega^{km}\omega^{mn} = \sum_m\omega^{m(k+n)}$. If $k+n\equiv0\pmod N$, every term is $1$: sum $N$. Otherwise $\omega^{k+n}$ is a nontrivial $N$-th root of unity and the geometric sum is $0$ (as in Theorem 3.19). So $\mathbf F^2 = NP$ with $P_{kn} = 1$ iff $k+n\equiv0$: $P$ reverses the vector, $(Pf)[n] = f[-n\bmod N]$, the discrete form of $\mathcal F^2f(x) = f(-x)$. $P^2 = I$, so $\mathbf F^4 = N^2P^2 = N^2I$, i.e. $(\mathbf F/\sqrt N)^4 = I$. If $(\mathbf F/\sqrt N)v = \lambda v$, $v\ne0$, then $v = \lambda^4v$ so $\lambda^4 = 1$: $\lambda\in\{1,-1,i,-i\}$. ∎ (With our sign convention the multiplicities of $1, -1, -i, i$ are $\lfloor N/4\rfloor+1,\ \lfloor(N+2)/4\rfloor,\ \lfloor(N+1)/4\rfloor,\ \lfloor(N-1)/4\rfloor$, a classical result we do not need; for $N = 8$ that is $3, 2, 2, 1$. The project checks the formula for $N = 1,\dots,16$. With the opposite sign, $e^{+2\pi i/N}$, the roles of $i$ and $-i$ swap.)

(b) $Rg = 2\pi xg - g' = 2\pi xg+2\pi xg = 4\pi xg$. Then $R^2g = R(4\pi xg) = 2\pi x\cdot4\pi xg - (4\pi xg)' = 8\pi^2x^2g - 4\pi g - 4\pi xg' = 8\pi^2x^2g-4\pi g+8\pi^2x^2g$ (using $g' = -2\pi xg$). So $h_2 = (16\pi^2x^2-4\pi)\,e^{-\pi x^2}$.

(c) By (3.9) applied twice, $\widehat{x^2g} = \left(\frac{i}{2\pi}\right)^2\hat g'' = -\frac{1}{4\pi^2}\frac{d^2}{d\xi^2}e^{-\pi\xi^2}$. Compute: $\frac{d}{d\xi}e^{-\pi\xi^2} = -2\pi\xi e^{-\pi\xi^2}$, $\frac{d^2}{d\xi^2}e^{-\pi\xi^2} = (-2\pi+4\pi^2\xi^2)e^{-\pi\xi^2}$. So $\widehat{x^2g} = -\frac{1}{4\pi^2}(4\pi^2\xi^2-2\pi)e^{-\pi\xi^2} = \left(\frac{1}{2\pi}-\xi^2\right)e^{-\pi\xi^2}$. Then
$$\hat h_2(\xi) = 16\pi^2\left(\frac1{2\pi}-\xi^2\right)e^{-\pi\xi^2} - 4\pi e^{-\pi\xi^2} = (8\pi-16\pi^2\xi^2-4\pi)e^{-\pi\xi^2} = -(16\pi^2\xi^2-4\pi)e^{-\pi\xi^2} = -h_2(\xi).$$
Eigenvalue $-1 = (-i)^2$. ✓ ($h_2$ is, up to normalization, the second Hermite function; in QM-3 it is the $n=2$ oscillator state, which is even in $x$ and has one node on each side.)

**C3.** Reference implementation (also in `project.py`):

```python
def hann(N):
    """Hann window w[n] = sin^2(pi n / N), the discrete form of cos^2(pi x / T) on |x| < T/2."""
    n = np.arange(N)
    return np.sin(np.pi * n / N) ** 2

def first_sidelobe_db(A):
    """Level (dB) of the first local maximum after the main lobe of a magnitude spectrum A, A[0] = peak."""
    i = 1
    while A[i] < A[i - 1]: i += 1          # descend the main lobe to the first zero
    while A[i + 1] >= A[i]: i += 1         # climb to the first sidelobe peak
    return 20 * np.log10(A[i] / A[0]), i

N, pad = 4096, 64
W = np.abs(np.fft.fft(hann(N), pad * N)); W /= W.max()
R = np.abs(np.fft.fft(np.ones(N), pad * N)); R /= R.max()
```

Measured: rect first sidelobe $-13.26$ dB, Hann $-31.47$ dB. Main lobe (peak to first zero): rect $1$ bin ($1/N$ cycles per sample, the first zero of $\operatorname{sinc}(N\xi)$), Hann $2$ bins. Up to a linear phase factor, the Hann transform is $\frac12\operatorname{sinc}(N\xi)+\frac14\operatorname{sinc}(N\xi-1)+\frac14\operatorname{sinc}(N\xi+1)$ (times $N$; this is §3.16's three sincs, and for $N = 4096$ the sinc form agrees with the exact zero-padded DFT to better than $10^{-11}$ of the peak over the first $200$ bins, which the project checks). Writing $s = \sin(\pi N\xi)/\pi$ and $u = N\xi$, it equals $s\left[\frac{1}{2u}-\frac{1}{4(u-1)}-\frac{1}{4(u+1)}\right] = \frac{-s}{2u(u^2-1)}$, because $\operatorname{sinc}(u\mp1) = -s/(u\mp1)$. At $u = 1$ the factor $s$ and the pole $1/(u-1)$ cancel to a nonzero value, so the first zero is at $u = 2$, i.e. $\xi = 2/N$; and the $1/u^3$ is the faster sidelobe decay. On the $2.5$-cycle cosine at $N = 64$, the rect-windowed DFT decays only like $1/|k-2.5|$ away from the tone, and the tails of the two lobes (at $+2.5$ and at $-2.5\equiv61.5$) add: the project measures about $10\%$ of the peak still leaking past bin $10$ and about $4.5\%$ at bin $30$. The Hann-windowed DFT is at $0.1\%$ of the peak by bin $10$ and below $10^{-5}$ by bin $30$. The trade-off in one sentence: a window with smoother edges buys sidelobes that decay like $1/\xi^3$ instead of $1/\xi$ (less leakage) at the price of a main lobe twice as wide (less resolution, since two tones closer than $2/T$ merge).

**C4.** (a) $\widehat{\operatorname{sgn}(x)e^{-a|x|}}(\xi) = \int_0^\infty e^{-(a+2\pi i\xi)x}dx - \int_{-\infty}^0e^{(a-2\pi i\xi)x}dx = \frac{1}{a+2\pi i\xi}-\frac{1}{a-2\pi i\xi} = \frac{(a-2\pi i\xi)-(a+2\pi i\xi)}{a^2+4\pi^2\xi^2} = \frac{-4\pi i\xi}{a^2+4\pi^2\xi^2}$. For fixed $\xi\ne0$, as $a\to0^+$ this tends to $\frac{-4\pi i\xi}{4\pi^2\xi^2} = \frac{-i}{\pi\xi} = \frac{1}{\pi i\xi}$.

(b) Split the integral at $|\xi| = 1$. On $|\xi|>1$, $\varphi(\xi)/\xi$ is integrable because $\varphi$ decays rapidly, so that part needs no limit. On $|\xi|<1$, write $\varphi(\xi) = \varphi(0)+\xi\psi(\xi)$ with $\psi(\xi) = (\varphi(\xi)-\varphi(0))/\xi$ for $\xi\ne0$ and $\psi(0) = \varphi'(0)$; by Taylor's theorem $\psi$ is smooth and bounded by $\sup_{|\xi|<1}|\varphi'|$. Then on $\epsilon<|\xi|<1$: $\int_{\epsilon<|\xi|<1}\frac{\varphi(\xi)}{\xi}d\xi = \varphi(0)\int_{\epsilon<|\xi|<1}\frac{d\xi}{\xi}+\int_{\epsilon<|\xi|<1}\psi(\xi)d\xi$. The first integral is exactly $0$ by oddness (the two sides cancel for every $\epsilon$); the second has a bounded integrand and converges as $\epsilon\to0$ to $\int_{-1}^1\psi$. So the limit exists. ∎ This is the **principal value** distribution $\mathrm{PV}\frac1\xi$; it is the odd symmetry that saves it, and the same argument fails for $1/|\xi|$ or $1/\xi^2$.

(c) The functions $\operatorname{sgn}(x)e^{-a|x|}$ converge to $\operatorname{sgn}(x)$ weakly as $a\to0^+$ (dominated convergence, paired with any test function). The transform is continuous on distributions, so $\widehat{\operatorname{sgn}} = \lim_{a\to0}\frac{-4\pi i\xi}{a^2+4\pi^2\xi^2}$ in the weak sense. Pair with $\varphi$ and split the integral at $|\xi| = 1$, as in (b). On $|\xi|>1$ the integrand $\frac{-4\pi i\xi\varphi(\xi)}{a^2+4\pi^2\xi^2}$ tends to $\frac{-i\varphi(\xi)}{\pi\xi}$ and is dominated by $|\varphi(\xi)|/(\pi|\xi|)$, which is integrable there, so that piece converges to $\frac{1}{\pi i}\int_{|\xi|>1}\frac{\varphi(\xi)}{\xi}d\xi$. On $|\xi|<1$ write $\varphi = \varphi(0)+\xi\psi$: the $\varphi(0)$ part integrates an odd function over a symmetric interval and vanishes for every $a$; the $\xi\psi$ part has integrand $\frac{-4\pi i\xi^2\psi(\xi)}{a^2+4\pi^2\xi^2}\to\frac{-i\psi(\xi)}{\pi}$, dominated by $|\psi|/\pi$, so it converges to $\frac{1}{\pi i}\int_{-1}^1\psi\,d\xi$, which by (b) is $\frac{1}{\pi i}\mathrm{PV}\!\int_{|\xi|<1}\frac{\varphi(\xi)}{\xi}d\xi$. Adding the two pieces, the limit is $\frac{1}{\pi i}\mathrm{PV}\!\int\frac{\varphi(\xi)}{\xi}d\xi$, i.e. $\widehat{\operatorname{sgn}} = \mathrm{PV}\frac{1}{\pi i\xi}$. Then $u = \frac12+\frac12\operatorname{sgn}$ and $\hat1 = \delta$ give
$$\hat u = \frac12\delta+\mathrm{PV}\frac{1}{2\pi i\xi}.$$

(d) $u' = \delta$ (§3.10), so by (3.8) extended to distributions, $2\pi i\xi\,\hat u = \hat\delta = 1$. Check with (c): $2\pi i\xi\cdot\frac12\delta = \pi i\,(\xi\delta) = 0$ since $\xi\delta = 0$; and $2\pi i\xi\cdot\mathrm{PV}\frac{1}{2\pi i\xi} = 1$ (multiplying $\mathrm{PV}\frac1\xi$ by the smooth function $\xi$ gives the constant $1$: $\langle\xi\,\mathrm{PV}\tfrac1\xi,\varphi\rangle = \mathrm{PV}\!\int\frac{\xi\varphi(\xi)}{\xi}d\xi = \int\varphi$). Sum: $0+1 = 1 = \hat\delta$. ✓ Note that the derivative rule alone could not have found the $\frac12\delta$: multiplying by $\xi$ kills any multiple of $\delta$, so (3.8) determines $\hat u$ only up to $C\delta$. The constant $C = \frac12$ is the mean value of $u$, fixed by the limiting argument in (c). This is the distributional version of "integration determines a function up to a constant."
