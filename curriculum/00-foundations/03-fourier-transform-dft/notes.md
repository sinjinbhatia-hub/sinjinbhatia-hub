# F3: The Fourier Transform, Delta Functions, Convolution, the DFT and FFT, and Sampling

**Time:** ~12 hours · **Prerequisites:** F2 · **Unlocks:** QM-2, EM-7, ST-8, ML-8, QC-5

> **Conventions in this module.** Fourier transform: $\hat f(\xi) = \int_{-\infty}^{\infty} f(x)\,e^{-2\pi i x\xi}\,dx$ and inverse $f(x) = \int_{-\infty}^{\infty} \hat f(\xi)\,e^{2\pi i x\xi}\,d\xi$ (the Osgood/Stein convention). No constants in front. $\xi$ is in **cycles per unit of $x$** (hertz if $x$ is seconds); the angular frequency is $\omega = 2\pi\xi$. We also write $\mathcal{F}f$ for $\hat f$. $\operatorname{sinc}(\xi) = \sin(\pi\xi)/(\pi\xi)$ with $\operatorname{sinc}(0) = 1$. $\Pi(x)$ is the unit rectangle: $1$ for $|x| < \tfrac12$, $0$ for $|x| > \tfrac12$. $\Lambda(x) = 1 - |x|$ for $|x| \le 1$ and $0$ otherwise (triangle). $u(x)$ is the Heaviside step ($0$ for $x<0$, $1$ for $x>0$). $\delta$ is the Dirac delta, $\delta_a(x) = \delta(x-a)$, and $\mathrm{III}(x) = \sum_{n\in\mathbb{Z}}\delta(x-n)$ is the comb. Convolution: $(f*g)(x) = \int f(y)\,g(x-y)\,dy$. The $L^2$ inner product is $\langle f, g\rangle = \int \overline{f(x)}\,g(x)\,dx$ (conjugate on the **first** slot, as in LA-8 and QM), $\|f\|^2 = \langle f,f\rangle$. The pairing of a distribution $T$ with a test function $\varphi$ is written $\langle T, \varphi\rangle$ and is **bilinear** (no conjugate); for an ordinary function $\langle f, \varphi\rangle = \int f\varphi\,dx$. DFT of $f\in\mathbb{C}^N$: $F[k] = \sum_{n=0}^{N-1} f[n]\,e^{-2\pi i kn/N}$ (NumPy's `np.fft.fft`, no $1/N$); $\omega = e^{-2\pi i/N}$; the Fourier matrix has entries $F_{kn} = \omega^{kn}$ for $0 \le k, n \le N-1$; the unitary version is $F/\sqrt N$. Circular convolution: $(f\circledast g)[n] = \sum_{m=0}^{N-1} f[m]\,g[(n-m)\bmod N]$. $\bar z$ is the complex conjugate. **Other books put the $2\pi$ elsewhere; the conversion table is in §3.1 and the sanity check is always $\hat f(0) = \int f$.**

## 0. Where this sits

F2 gave you Fourier series: a function on an interval of length $T$ is an orthogonal projection onto the basis $e^{2\pi i n x/T}$, with coefficients $c_n = \frac1T\int_{-T/2}^{T/2} f(x)e^{-2\pi i n x/T}dx$. That needs a period. This module removes the period ($T\to\infty$), which turns the sum into an integral and the coefficient list into a function $\hat f(\xi)$. It then adds four things you cannot do without: the delta function (made rigorous), convolution and why it becomes multiplication, the discrete transform that computers actually run (and why it is fast), and the sampling theorem that says when a list of numbers is as good as a continuous signal.

You need from before: complex exponentials and roots of unity (F1), inner products and orthogonal projection (LA-5, LA-8), Fourier series and the coefficient formula (F2), ODEs with constant coefficients, and integration by parts. Afterwards: QM-2 builds the momentum-space wavefunction and the Heisenberg principle directly on §3.8–3.9; EM-7 uses the transform for dispersion and the far-field diffraction pattern; ST-8's periodogram is $|F[k]|^2/N$ from §3.12; ML-8's convolution layer is the circulant matrix of §3.13; QC-5's quantum Fourier transform is the unitary matrix $F/\sqrt N$ of Theorem 3.18 with $N = 2^n$.

## 1. The problem that forced this idea

**The infinite rod.** Fourier's 1807 memoir and his 1822 *Théorie analytique de la chaleur* solved the heat equation on a finite bar by expanding the initial temperature in sines and cosines. The bar's length set the period. But heat also flows in a rod so long that its ends do not matter, and a single hot spot in the middle of it does not repeat. There is no period to expand in. Fourier's answer, which also appears in the 1822 book (Cauchy and Poisson derived versions of it independently around 1815–1816), was to let the period grow without bound. The discrete frequencies $n/T$ crowd together, the sum over $n$ becomes an integral over a continuous frequency $\xi$, and the coefficients become a function $\hat f(\xi)$. That limit is the opening derivation of §3.

**A function that is zero except at one point.** Engineers needed an idealized hammer blow: a force that acts for no time at all but delivers a finite impulse. Oliver Heaviside used such objects freely in his operational calculus in the 1890s. Paul Dirac introduced the symbol $\delta(x)$ in 1927 and made it a tool in *The Principles of Quantum Mechanics* (1930), with the defining property $\int \delta(x)\varphi(x)\,dx = \varphi(0)$. Mathematicians objected, correctly: no function has that property (we prove this in Theorem 3.13). John von Neumann called it a "fiction" in 1932. The resolution came from Sergei Sobolev (1936) and, in full generality, Laurent Schwartz (1945–1950): $\delta$ is not a function but a *linear functional*, a rule that eats a test function and returns a number. Once you accept that definition, every manipulation Dirac did becomes a theorem, and $\delta$ acquires a Fourier transform: the constant $1$.

**Computing it fast.** In 1805 Gauss wanted to interpolate the orbits of the asteroids Pallas and Juno from a handful of observations, which meant computing trigonometric coefficients from $12$ samples. He noticed that splitting the $12$ points into $3$ groups of $4$ (or $4$ of $3$) and recombining saved most of the arithmetic. He never published it; it appears in his collected works (1866), in Latin, and stayed unnoticed. In 1965 James Cooley and John Tukey published the same idea for $N = 2^m$ samples, motivated by the need to analyze seismometer arrays (part of detecting Soviet nuclear tests). Their algorithm computes the $N$-point transform in about $N\log_2 N$ operations instead of $N^2$. For $N = 2^{20}$ that is a factor of $50{,}000$. Spectral analysis, digital filtering, MRI reconstruction, JPEG, and every fast convolution in machine learning run on it. Heideman, Johnson, and Burrus traced the history back to Gauss in 1984.

Three problems, one object. The transform is the common language.

## 2. Intuition

### A signal is a sum of pure tones, with a continuous dial

In F2 a periodic signal was a sum of tones at frequencies $0, \pm1/T, \pm2/T, \dots$. Take a non-periodic signal, say one click. Which tones are in it? All of them. A single click has to be built from tones at *every* frequency, because tones at a discrete set of frequencies would add up to something that repeats. So the "list of coefficients" $c_n$ must become a function of a continuous variable: $\hat f(\xi)$ tells you how much of the tone $e^{2\pi i\xi x}$ is in $f$.

How do you measure "how much of $e^{2\pi i \xi x}$ is in $f$"? The same way F2 measured $c_n$: take the inner product with the tone. $\langle e^{2\pi i\xi x}, f\rangle = \int \overline{e^{2\pi i \xi x}} f(x)\,dx = \int f(x) e^{-2\pi i x\xi}\,dx$. That is the definition of $\hat f(\xi)$. The minus sign in the exponent is the conjugate on the first slot of the inner product. Nothing more.

### Sharp in one domain means spread in the other

Hit a drum: a short, sharp sound, and you hear a thud with no definite pitch. Strike a tuning fork: a long tone with one definite pitch. A short signal needs many frequencies to cancel each other everywhere except in a brief window; a long sinusoid needs only one. The transform makes this quantitative:

```
   f(x)                          |f^(xi)|
   ___                            __________
  |   |  narrow rect              /          \  wide sinc
__|   |__                     ___/            \___

   _______                         _
  |       |  wide rect            | |  narrow sinc
__|       |__                 ____| |____
```

Squeeze $f$ horizontally by a factor $a$ and $\hat f$ stretches by $a$. You should be able to guess that a *product* of widths is bounded below. It is: $\sigma_x\sigma_\xi \ge 1/(4\pi)$, with equality for Gaussians (§3.9). In QM-2 this exact theorem, with $\hbar$ pasted on, is the Heisenberg uncertainty principle.

### Convolution is smearing, and in frequency it is multiplication

Convolution $(f*g)(x) = \int f(y)g(x-y)\,dy$ is a moving weighted average: to get the output at $x$, center a flipped copy of $g$ at $x$, multiply by $f$, integrate. Convolve a rectangle with itself and you get a triangle:

```
  _____        _____          /\
 |     |   *  |     |   =    /  \
_|     |_    _|     |_    __/    \__
```

Every blur, every low-pass filter, every linear time-invariant system is a convolution. Now guess what convolution does to frequencies. Averaging smooths, and smoothing kills fast wiggles, so the high-frequency content of $f*g$ should be *reduced by a factor that depends on $g$*. The exact statement is the convolution theorem: $\widehat{f*g} = \hat f\,\hat g$. Smearing in $x$ is pointwise multiplication in $\xi$. This one fact is why the FFT matters: convolution costs $N^2$ directly and $N\log N$ through the transform.

### Sampling multiplies by a comb, which makes copies in frequency

Measuring $f$ every $\Delta$ seconds is multiplying $f$ by a comb of spikes $\sum_n \delta(x - n\Delta)$. The comb is periodic with period $\Delta$, so (F2) it is a sum of tones at frequencies $k/\Delta$. Multiplying $f$ by a tone at frequency $k/\Delta$ shifts $\hat f$ by $k/\Delta$ (that is the modulation rule, §3.4). So the spectrum of the sampled signal is the original spectrum **plus shifted copies at every multiple of $1/\Delta$**:

```
  |f^|       _/\_                            original
             ___________________________
  |f_s^|  _/\_   _/\_   _/\_   _/\_   _/\_    copies every 1/Delta
            -1/D    0    1/D   2/D
```

If the original spectrum lives in $|\xi| < B$ and the copies are spaced $1/\Delta > 2B$ apart, they do not overlap, and you can cut out the middle copy and recover $f$ exactly from its samples. If they overlap, high frequencies masquerade as low ones and nothing warns you. That is the Nyquist–Shannon theorem and aliasing, in one picture. The wagon wheel in a western that seems to turn backwards is this picture.

### The finite version

With $N$ samples there are only $N$ distinguishable frequencies, so the transform becomes an $N\times N$ matrix $F$ with entries $\omega^{kn}$, $\omega = e^{-2\pi i/N}$. Its columns are the $N$ discrete tones. They are orthogonal (F1's geometric sum of roots of unity), so $F/\sqrt N$ is unitary: the DFT is a rotation of $\mathbb{C}^N$, and it preserves lengths (Parseval). Every statement above has a finite twin: convolution becomes circular convolution, the convolution theorem survives, and "convolution operator" becomes "circulant matrix," which $F$ diagonalizes.

## 3. Formal development

### 3.1 From Fourier series to the transform: the $T\to\infty$ limit

▶ **Try it first.** Write the F2 series for a function on $[-T/2, T/2]$ with $\xi_n = n/T$ and $\Delta\xi = 1/T$. Rearrange it to look like a Riemann sum in $\xi$.

Take $f$ to be integrable, zero outside $[-T/2, T/2]$ for some large $T$ (any decaying $f$ is approximately like this). On that interval, F2 gives

$$
f(x) = \sum_{n=-\infty}^{\infty} c_n\,e^{2\pi i n x/T}, \qquad c_n = \frac1T\int_{-T/2}^{T/2} f(y)\,e^{-2\pi i n y/T}\,dy.
$$

**Derivation.**

1. Set $\xi_n = n/T$ and $\Delta\xi = \xi_{n+1}-\xi_n = 1/T$. Define $\hat f_T(\xi) = \int_{-T/2}^{T/2} f(y)e^{-2\pi i y\xi}\,dy$. Then $c_n = \Delta\xi\cdot\hat f_T(\xi_n)$.
2. Substitute into the series: $f(x) = \sum_n \hat f_T(\xi_n)\,e^{2\pi i\xi_n x}\,\Delta\xi$.
3. This is a Riemann sum for $\int \hat f_T(\xi)e^{2\pi i\xi x}\,d\xi$ on a grid of spacing $\Delta\xi = 1/T$.
4. Since $f$ vanishes outside $[-T/2, T/2]$, $\hat f_T(\xi) = \int_{-\infty}^{\infty} f(y)e^{-2\pi i y\xi}dy =: \hat f(\xi)$, independent of $T$. Let $T\to\infty$: the grid spacing $\to 0$ and the Riemann sum becomes the integral. ∎

**Definition 3.1 (Fourier transform and inverse).**
$$
\hat f(\xi) = \int_{-\infty}^{\infty} f(x)\,e^{-2\pi i x\xi}\,dx, \tag{3.1}
$$
$$
f(x) = \int_{-\infty}^{\infty} \hat f(\xi)\,e^{2\pi i x\xi}\,d\xi. \tag{3.2}
$$

(3.1) is **analysis** (take $f$ apart into tones); (3.2) is **synthesis** (put it back). Step 3 is a heuristic: a Riemann sum of an improper integral over a shrinking grid need not converge to the integral. The honest statement is the **Fourier inversion theorem**: if $f$ and $\hat f$ are both integrable, then (3.2) holds at every point where $f$ is continuous. Stein & Shakarchi prove it in ch. 5 §1 for rapidly decreasing functions by a Gaussian-approximation argument; Osgood ch. 2 gives the engineering version. We take it as proved and use it freely.

■ **What just happened.** The Fourier coefficient $c_n$ is $\hat f(n/T)$ times the grid spacing $1/T$. The series is a Riemann sum for the inverse transform. Letting the period go to infinity turns "which harmonics" into "how much at each frequency."

**Conventions, once.** Three conventions are in common use. Ours puts $2\pi$ in the exponent and nothing in front. With angular frequency $\omega = 2\pi\xi$:

| Convention | Forward | Inverse | Relation to ours |
|---|---|---|---|
| **Ours** (Osgood, Stein, signal processing, NumPy's `fftfreq`) | $\hat f(\xi) = \int f(x)e^{-2\pi i x\xi}dx$ | $f(x) = \int \hat f(\xi) e^{2\pi i x\xi}d\xi$ | — |
| **Angular, asymmetric** (most engineering ODE/PDE texts) | $F(\omega) = \int f(x) e^{-i\omega x}dx$ | $f(x) = \frac{1}{2\pi}\int F(\omega)e^{i\omega x}d\omega$ | $F(\omega) = \hat f(\omega/2\pi)$ |
| **Angular, symmetric** (QM-2, Griffiths, most physics) | $F(k) = \frac{1}{\sqrt{2\pi}}\int f(x)e^{-ikx}dx$ | $f(x) = \frac{1}{\sqrt{2\pi}}\int F(k)e^{ikx}dk$ | $F(k) = \frac{1}{\sqrt{2\pi}}\hat f(k/2\pi)$ |

The rules in §3.4 change constants when you change convention (the derivative rule becomes $i\omega$, the convolution theorem acquires a $\sqrt{2\pi}$ in the symmetric convention, and so on). **The one check that never lies:** in every convention, the forward transform at zero frequency is a known multiple of $\int f$. In ours, $\hat f(0) = \int f(x)\,dx$ exactly. Check it on every transform you compute.

### 3.2 Existence, boundedness, continuity, decay

**Theorem 3.2 (Integrable $f$ has a bounded continuous transform).** If $\int|f|\,dx < \infty$, then $\hat f(\xi)$ exists for every $\xi$, $|\hat f(\xi)| \le \int|f|\,dx$, and $\hat f$ is continuous.

**Derivation.**

1. Existence and bound: $|f(x)e^{-2\pi i x\xi}| = |f(x)|$, so the integral (3.1) converges absolutely and $|\hat f(\xi)| \le \int|f|$.
2. Continuity: $|\hat f(\xi+h) - \hat f(\xi)| = \left|\int f(x)e^{-2\pi i x\xi}\left(e^{-2\pi i x h}-1\right)dx\right| \le \int |f(x)|\,|e^{-2\pi i x h}-1|\,dx$.
3. For each fixed $x$, $|e^{-2\pi i xh}-1|\to 0$ as $h\to0$, and the integrand is bounded by $2|f(x)|$, which is integrable. Passing the limit inside the integral is the dominated convergence theorem (an analysis fact we take as given). So the right side $\to0$. ∎

**Theorem 3.3 (Riemann–Lebesgue).** If $\int|f| < \infty$ then $\hat f(\xi)\to0$ as $|\xi|\to\infty$.

▶ **Try it first.** Shifting $x$ by half a period of $e^{-2\pi i x\xi}$ flips the sign of the exponential. Use that to write $\hat f(\xi)$ as half the transform of a *difference* $f(x) - f(x-h)$.

**Derivation.**

1. Let $h = 1/(2\xi)$. Then $e^{-2\pi i (x+h)\xi} = e^{-2\pi i x\xi}e^{-\pi i} = -e^{-2\pi i x\xi}$.
2. So $\hat f(\xi) = -\int f(x)e^{-2\pi i(x+h)\xi}dx = -\int f(y-h)e^{-2\pi i y\xi}dy$ (substitute $y = x+h$).
3. Average this with (3.1): $\hat f(\xi) = \frac12\int\big[f(y) - f(y-h)\big]e^{-2\pi i y\xi}dy$, so $|\hat f(\xi)| \le \frac12\int|f(y)-f(y-h)|\,dy$.
4. As $|\xi|\to\infty$, $h\to0$, and $\int|f(y)-f(y-h)|dy\to0$ (translation is continuous in the $L^1$ norm: true for continuous compactly supported $f$ by uniform continuity, and for general integrable $f$ by approximating with such functions). ∎

■ **What just happened.** Integrability of $f$ buys you a bounded, continuous transform that decays. It does *not* buy you a rate of decay. The rate is controlled by the smoothness of $f$: a jump gives $1/|\xi|$ decay (the rect, next), a kink gives $1/\xi^2$ (the triangle and two-sided exponential), and infinite smoothness with decay gives faster-than-any-power decay (the Gaussian). The derivative rule in §3.4 is why.

### 3.3 First examples, each by direct integration

**Example 3.A (rect $\to$ sinc).** $\Pi(x) = 1$ for $|x|<\tfrac12$, else $0$.

1. $\hat\Pi(\xi) = \int_{-1/2}^{1/2} e^{-2\pi i x\xi}\,dx = \left[\dfrac{e^{-2\pi i x\xi}}{-2\pi i\xi}\right]_{-1/2}^{1/2} = \dfrac{e^{-\pi i\xi} - e^{\pi i\xi}}{-2\pi i\xi}$ for $\xi\neq0$.
2. $e^{-\pi i\xi}-e^{\pi i\xi} = -2i\sin(\pi\xi)$. So $\hat\Pi(\xi) = \dfrac{-2i\sin\pi\xi}{-2\pi i\xi} = \dfrac{\sin\pi\xi}{\pi\xi}$.
3. At $\xi=0$, directly: $\hat\Pi(0) = \int_{-1/2}^{1/2}1\,dx = 1 = \operatorname{sinc}(0)$. ✓ (This is the $\hat f(0) = \int f$ check.)

$$
\widehat{\Pi}(\xi) = \operatorname{sinc}(\xi) = \frac{\sin\pi\xi}{\pi\xi}. \tag{3.3}
$$

The zeros are at the nonzero integers. $\operatorname{sinc}$ decays like $1/|\xi|$ and is **not** absolutely integrable: the jump in $\Pi$ costs you slow decay. (The transform exists anyway, and $\int\operatorname{sinc} = 1$ as an improper integral, consistent with $\Pi(0)=1$ by (3.2).)

**Example 3.B (one-sided exponential).** $f(x) = e^{-ax}u(x)$, $a>0$.

1. $\hat f(\xi) = \int_0^\infty e^{-ax}e^{-2\pi i x\xi}dx = \int_0^\infty e^{-(a+2\pi i\xi)x}dx$.
2. The antiderivative is $e^{-(a+2\pi i\xi)x}/(-(a+2\pi i\xi))$; at $\infty$ it vanishes because $|e^{-(a+2\pi i\xi)x}| = e^{-ax}\to0$; at $0$ it is $-1/(a+2\pi i\xi)$.

$$
\widehat{e^{-ax}u(x)}(\xi) = \frac{1}{a+2\pi i\xi}. \tag{3.4}
$$

Check: $\hat f(0) = 1/a = \int_0^\infty e^{-ax}dx$. ✓ Note $|\hat f(\xi)|^2 = 1/(a^2+4\pi^2\xi^2)$ and $\hat f$ decays like $1/|\xi|$ (again a jump at $0$).

**Example 3.C (two-sided exponential $\to$ Lorentzian).** $f(x) = e^{-a|x|}$, $a>0$.

1. Split at $0$: $\hat f(\xi) = \int_0^\infty e^{-(a+2\pi i\xi)x}dx + \int_{-\infty}^0 e^{(a-2\pi i\xi)x}dx = \dfrac{1}{a+2\pi i\xi} + \dfrac{1}{a-2\pi i\xi}$.
2. Common denominator: $\dfrac{(a-2\pi i\xi)+(a+2\pi i\xi)}{a^2+4\pi^2\xi^2}$.

$$
\widehat{e^{-a|x|}}(\xi) = \frac{2a}{a^2+4\pi^2\xi^2}. \tag{3.5}
$$

Check: $\hat f(0) = 2/a = \int e^{-a|x|}dx$. ✓ This is the **Lorentzian** (the Cauchy density of F4, up to normalization). The kink at $0$ gives $1/\xi^2$ decay; real and even $f$ gave real and even $\hat f$ (§3.4(g)).

### 3.4 Properties, each derived

**Theorem 3.4 (Rules of the transform).** Let $f, g$ be integrable, $a, \nu\in\mathbb{R}$, $a\ne0$, $c_1, c_2\in\mathbb{C}$.

**(a) Linearity.** $\widehat{c_1 f + c_2 g} = c_1\hat f + c_2\hat g$. *Derivation:* the integral (3.1) is linear in $f$. ∎

**(b) Shift.**
$$
\widehat{f(x-a)}(\xi) = e^{-2\pi i a\xi}\,\hat f(\xi). \tag{3.6}
$$
*Derivation:* $\int f(x-a)e^{-2\pi i x\xi}dx$; put $y = x-a$, $dx = dy$: $\int f(y)e^{-2\pi i(y+a)\xi}dy = e^{-2\pi i a\xi}\int f(y)e^{-2\pi i y\xi}dy$. ∎ A delay changes the phase linearly in $\xi$ and leaves $|\hat f|$ alone.

**(c) Modulation.** $\widehat{e^{2\pi i\nu x}f(x)}(\xi) = \hat f(\xi-\nu)$. *Derivation:* $\int f(x)e^{2\pi i\nu x}e^{-2\pi i x\xi}dx = \int f(x)e^{-2\pi i x(\xi-\nu)}dx$. ∎ Consequence: $\widehat{\cos(2\pi\nu x)f(x)} = \tfrac12[\hat f(\xi-\nu)+\hat f(\xi+\nu)]$.

**(d) Scaling (the stretch theorem).**
$$
\widehat{f(ax)}(\xi) = \frac{1}{|a|}\,\hat f\!\left(\frac{\xi}{a}\right). \tag{3.7}
$$
*Derivation:* put $y = ax$, $dx = dy/a$. For $a>0$: $\int f(y)e^{-2\pi i (y/a)\xi}\,dy/a = \frac1a\hat f(\xi/a)$. For $a<0$ the limits of integration flip, which contributes a sign that turns $1/a$ into $1/|a|$. ∎ Compress $f$ (large $|a|$) and $\hat f$ spreads out and drops in height. The area $\int\hat f\,d\xi = f(0)$ is conserved.

**(e) Derivative.** If $f$ is differentiable, $f'$ is integrable, and $f(x)\to0$ as $|x|\to\infty$:
$$
\widehat{f'}(\xi) = 2\pi i\xi\,\hat f(\xi). \tag{3.8}
$$
*Derivation:* integrate by parts. $\int f'(x)e^{-2\pi i x\xi}dx = \big[f(x)e^{-2\pi i x\xi}\big]_{-\infty}^{\infty} - \int f(x)\,(-2\pi i\xi)e^{-2\pi i x\xi}dx$. The boundary term is $0$ because $f\to0$ and $|e^{-2\pi i x\xi}| = 1$. ∎ Differentiation multiplies by $2\pi i\xi$: it amplifies high frequencies. Repeating: $\widehat{f^{(n)}} = (2\pi i\xi)^n\hat f$. Read backwards, this is why smooth $f$ has a fast-decaying $\hat f$: if $f^{(n)}$ is integrable then $|\hat f(\xi)| \le \|f^{(n)}\|_1/(2\pi|\xi|)^n$ by Theorem 3.2 applied to $f^{(n)}$.

**(e′) Multiplication by $x$ (the reverse rule).** If $\int|x f(x)|dx<\infty$:
$$
\widehat{-2\pi i x\,f(x)}(\xi) = \frac{d}{d\xi}\hat f(\xi), \quad\text{equivalently}\quad \widehat{x f(x)}(\xi) = \frac{i}{2\pi}\,\hat f'(\xi). \tag{3.9}
$$
*Derivation:* differentiate (3.1) under the integral sign (allowed because $|x f(x)|$ is integrable and dominates the $\xi$-derivative of the integrand): $\hat f'(\xi) = \int f(x)\,(-2\pi i x)\,e^{-2\pi i x\xi}dx$. Divide by $-2\pi i$, using $1/(-2\pi i) = i/(2\pi)$. ∎

**(f) Duality.** Applying the transform twice reflects:
$$
\hat{\hat f}(x) = f(-x). \tag{3.10}
$$
*Derivation:* $\hat{\hat f}(x) = \int\hat f(\xi)e^{-2\pi i\xi x}d\xi = \int\hat f(\xi)e^{2\pi i\xi(-x)}d\xi = f(-x)$ by (3.2). ∎ Consequence: every transform pair reads both ways. From (3.3): $\widehat{\operatorname{sinc}}(\xi) = \Pi(-\xi) = \Pi(\xi)$. From (3.5) with $a = 2\pi$: $\widehat{1/(1+x^2)}(\xi) = \pi e^{-2\pi|\xi|}$ (exercise A6).

**(g) Conjugation and symmetry.** $\widehat{\bar f}(\xi) = \overline{\hat f(-\xi)}$. *Derivation:* $\int\overline{f(x)}e^{-2\pi i x\xi}dx = \overline{\int f(x)e^{2\pi i x\xi}dx} = \overline{\hat f(-\xi)}$. ∎ So if $f$ is **real**, $\hat f(-\xi) = \overline{\hat f(\xi)}$ (**Hermitian** symmetry): $|\hat f|$ is even, $\arg\hat f$ is odd, and the negative frequencies carry no new information. Real and even $f$ gives real and even $\hat f$; real and odd $f$ gives purely imaginary odd $\hat f$.

**(h) Reflection.** $\widehat{f(-x)}(\xi) = \hat f(-\xi)$: (3.7) with $a=-1$.

■ **What just happened.** Every rule is a one-line substitution in (3.1). The two you will use most: shift ↔ phase (3.6), and derivative ↔ multiply by $2\pi i\xi$ (3.8). Together they turn constant-coefficient ODEs into algebra (§3.11). The stretch theorem (3.7) is the uncertainty principle in embryo.

### 3.5 The Gaussian, by the ODE method

▶ **Try it first.** $g(x) = e^{-\pi x^2}$ satisfies $g' = -2\pi x\,g$. Transform both sides of that ODE using (3.8) and (3.9). What ODE does $\hat g$ satisfy?

**Theorem 3.5 (The Gaussian is its own transform).** $\widehat{e^{-\pi x^2}} = e^{-\pi\xi^2}$.

**Derivation.**

1. $g(x) = e^{-\pi x^2}$ satisfies $g'(x) = -2\pi x\,g(x)$, with $g\to0$ at $\pm\infty$ and $xg$ integrable, so (3.8) and (3.9) apply.
2. Transform the left side with (3.8): $\widehat{g'} = 2\pi i\xi\,\hat g$.
3. Transform the right side with (3.9): $\widehat{-2\pi x g} = -2\pi\cdot\frac{i}{2\pi}\hat g' = -i\,\hat g'$.
4. Equate: $2\pi i\xi\,\hat g = -i\,\hat g'$. Divide by $-i$: $\hat g'(\xi) = -2\pi\xi\,\hat g(\xi)$. This is the **same ODE** that $g$ satisfies, so $\hat g(\xi) = C e^{-\pi\xi^2}$.
5. Fix $C$ with the check $\hat g(0) = \int g$. Let $I = \int e^{-\pi x^2}dx$. Then $I^2 = \iint e^{-\pi(x^2+y^2)}dx\,dy = \int_0^{2\pi}\!\!\int_0^\infty e^{-\pi r^2}r\,dr\,d\theta = 2\pi\cdot\frac{1}{2\pi} = 1$. So $I = 1$, $C = 1$. ∎

$$
\widehat{e^{-\pi x^2}}(\xi) = e^{-\pi\xi^2}. \tag{3.11}
$$

By the stretch theorem (3.7), for $a>0$, writing $e^{-ax^2} = g\big(\sqrt{a/\pi}\,x\big)$:
$$
\widehat{e^{-ax^2}}(\xi) = \sqrt{\frac{\pi}{a}}\,e^{-\pi^2\xi^2/a}. \tag{3.12}
$$
A Gaussian of standard deviation $\sigma$ in $x$ (that is $a = 1/(2\sigma^2)$) transforms to a Gaussian of standard deviation $1/(2\pi\sigma)$ in $\xi$. Widths multiply to the constant $1/(2\pi)$. Hold onto that.

■ **What just happened.** We never integrated $e^{-\pi x^2}e^{-2\pi i x\xi}$. We noticed that the Gaussian is characterized by a first-order ODE, that the transform converts that ODE into itself, and that the normalization is fixed by $\hat g(0) = \int g$. This is LA-7 thinking: the Gaussian is an eigenvector of $\mathcal F$ with eigenvalue $1$. §3.8 finds the rest of the eigenvectors.

### 3.6 Convolution and the convolution theorem

**Definition 3.6.** $(f*g)(x) = \int_{-\infty}^\infty f(y)\,g(x-y)\,dy$.

Substituting $y\mapsto x-y$ shows $f*g = g*f$. Convolution is linear in each argument and associative (Fubini).

*Interpretation 1 (moving weighted average).* Take $g(x) = \frac1w\Pi(x/w)$, a box of width $w$ and area $1$. Then $(f*g)(x) = \frac1w\int_{x-w/2}^{x+w/2}f(y)\,dy$: the running mean of $f$ over a window of width $w$ centered at $x$. A general $g$ is a general weighting. $g$ is called the **kernel**.

*Interpretation 2 (output of a linear time-invariant system).* §3.11 proves that any linear system whose behavior does not change with time has the form $y = h*x$, where $h$ is the response to a unit impulse.

**Theorem 3.7 (Convolution theorem).** For integrable $f, g$:
$$
\widehat{f*g}(\xi) = \hat f(\xi)\,\hat g(\xi). \tag{3.13}
$$

▶ **Try it first.** Write out $\widehat{f*g}$ as a double integral, swap the order, and recognize the shift rule (3.6) in the inner integral.

**Derivation.**

1. $\widehat{f*g}(\xi) = \int\left(\int f(y)g(x-y)\,dy\right)e^{-2\pi i x\xi}\,dx$.
2. Swap the order of integration (Fubini: the double integral of $|f(y)||g(x-y)|$ is $\|f\|_1\|g\|_1<\infty$): $= \int f(y)\left(\int g(x-y)e^{-2\pi i x\xi}\,dx\right)dy$.
3. The inner integral is the transform of $g$ shifted by $y$; by (3.6) it equals $e^{-2\pi i y\xi}\hat g(\xi)$. (Or substitute $u = x-y$ directly.)
4. $= \hat g(\xi)\int f(y)e^{-2\pi i y\xi}dy = \hat g(\xi)\hat f(\xi)$. ∎

**The mirror statement.** Multiplication in $x$ is convolution in $\xi$:
$$
\widehat{fg}(\xi) = (\hat f*\hat g)(\xi). \tag{3.14}
$$
*Derivation.* Write $f$ by (3.2): $f(x)g(x) = \int\hat f(\eta)e^{2\pi i\eta x}d\eta\;g(x)$. Transform: $\widehat{fg}(\xi) = \int\hat f(\eta)\left(\int g(x)e^{-2\pi i x(\xi-\eta)}dx\right)d\eta = \int\hat f(\eta)\hat g(\xi-\eta)\,d\eta = (\hat f*\hat g)(\xi)$. ∎

*Example: $\Pi*\Pi = \Lambda$.* For $|x|\le1$ the overlap of two unit boxes offset by $x$ has length $1-|x|$; for $|x|>1$ they do not overlap. So $\Pi*\Pi = \Lambda$ and by (3.13),
$$
\hat\Lambda(\xi) = \operatorname{sinc}^2(\xi).
$$
Worked example 4.1 confirms this by direct integration. A kink in $x$, $1/\xi^2$ in frequency, consistent with the smoothness rule.

■ **What just happened.** Convolution tangles two functions together; the transform untangles them into a product, because the exponential $e^{-2\pi i x\xi}$ factors across a shift: $e^{-2\pi i(u+y)\xi} = e^{-2\pi i u\xi}e^{-2\pi i y\xi}$. That one factorization is the whole theorem.

### 3.7 Parseval and Plancherel: the transform preserves inner products

**Theorem 3.8 (Parseval/Plancherel).** For $f, g$ with $f, g, \hat f, \hat g$ integrable,
$$
\langle f, g\rangle = \langle\hat f, \hat g\rangle, \qquad\text{in particular}\qquad \int|f(x)|^2dx = \int|\hat f(\xi)|^2d\xi. \tag{3.15}
$$

▶ **Try it first.** Let $h(x) = \overline{f(-x)}$. Compute $\hat h$ with rules (g) and (h). Then evaluate $(g*h)(0)$ two ways: by the definition, and by writing $g*h$ as the inverse transform of $\widehat{g*h}$.

**Derivation.**

1. Let $h(x) = \overline{f(-x)}$. By (g), $\widehat{\bar f}(\xi) = \overline{\hat f(-\xi)}$; then by (h), reflecting the argument flips $\xi$ again: $\hat h(\xi) = \overline{\hat f(\xi)}$.
2. By the definition of convolution, $(g*h)(0) = \int g(y)h(0-y)\,dy = \int g(y)\overline{f(y)}\,dy = \langle f, g\rangle$.
3. By (3.2) at $x=0$ and the convolution theorem (3.13), $(g*h)(0) = \int\widehat{g*h}(\xi)\,d\xi = \int\hat g(\xi)\hat h(\xi)\,d\xi = \int\hat g(\xi)\overline{\hat f(\xi)}\,d\xi = \langle\hat f, \hat g\rangle$.
4. Equate 2 and 3. Set $g = f$ for the norm statement. ∎

The integrability hypotheses are needed for the derivation as written. **Plancherel's theorem** extends (3.15) to all of $L^2$ (all $f$ with $\int|f|^2<\infty$) by a density argument, and shows $\mathcal F$ is a bijection of $L^2$ onto itself. We use the extension without proof; Stein & Shakarchi ch. 5 §1.7 has it.

■ **What just happened.** In LA-8 language: $\mathcal F$ is a **unitary operator** on $L^2(\mathbb{R})$. It preserves lengths and angles. Energy computed in time equals energy computed in frequency. The DFT version, Theorem 3.18, is a literal unitary matrix.

### 3.8 The transform as an operator: $\mathcal F^4 = I$ and the Hermite eigenfunctions

**Theorem 3.9.** On $L^2$, $\mathcal F^2 f(x) = f(-x)$ and $\mathcal F^4 = I$. Hence every eigenvalue of $\mathcal F$ is a fourth root of unity: $\lambda\in\{1, -i, -1, i\}$.

*Proof.* $\mathcal F^2 f = f(-\,\cdot)$ is (3.10). Reflecting twice is the identity, so $\mathcal F^4 = I$. If $\mathcal F h = \lambda h$ with $h\ne0$ then $h = \mathcal F^4h = \lambda^4 h$, so $\lambda^4 = 1$. ∎

Which functions are eigenfunctions? Theorem 3.5 says $g = e^{-\pi x^2}$ has eigenvalue $1$. Define the operators $X f(x) = x f(x)$ and $D f = f'$. Rules (3.8) and (3.9) say
$$
\mathcal F D = 2\pi i\,X\,\mathcal F, \qquad \mathcal F X = \frac{i}{2\pi}\,D\,\mathcal F.
$$

**Theorem 3.10 (Hermite eigenfunctions).** Let $R = 2\pi X - D$, so $(Rf)(x) = 2\pi x f(x) - f'(x)$. Then $\mathcal F R = -i\,R\,\mathcal F$, and consequently $h_n := R^n g$ satisfies $\mathcal F h_n = (-i)^n h_n$. Each $h_n$ is a polynomial of degree $n$ times $e^{-\pi x^2}$ (these are the **Hermite functions**, the harmonic-oscillator eigenstates of QM-3).

**Derivation.**

1. $\mathcal F R f = \mathcal F(2\pi X f) - \mathcal F(Df) = 2\pi\cdot\frac{i}{2\pi}D\mathcal Ff - 2\pi i X\mathcal F f = i\,D\mathcal F f - 2\pi i\,X\mathcal Ff = -i\,(2\pi X - D)\mathcal F f = -i\,R\,\mathcal F f$.
2. Induction on $n$: $\mathcal F h_0 = \mathcal F g = g = (-i)^0 h_0$. If $\mathcal F h_n = (-i)^n h_n$, then $\mathcal F h_{n+1} = \mathcal F R h_n = -i R\mathcal F h_n = -i(-i)^n R h_n = (-i)^{n+1}h_{n+1}$.
3. Form: $R(p(x)g(x)) = 2\pi x p g - p' g - p g' = (2\pi x p - p' + 2\pi x p)g = (4\pi x p - p')g$. If $p$ has degree $n$, $4\pi xp - p'$ has degree $n+1$ with nonzero leading coefficient $4\pi\cdot(\text{lead of }p)$. So $h_n$ is a degree-$n$ polynomial times $g$. ∎

Concretely $h_1 = 4\pi x e^{-\pi x^2}$, and you can check it directly: by (3.9), $\widehat{x e^{-\pi x^2}} = \frac{i}{2\pi}\frac{d}{d\xi}e^{-\pi\xi^2} = \frac{i}{2\pi}(-2\pi\xi)e^{-\pi\xi^2} = -i\,\xi e^{-\pi\xi^2}$. Eigenvalue $-i$. ✓ The $h_n$ are orthogonal (different eigenvalues of a unitary operator) and in fact form an orthogonal basis of $L^2$; completeness is proved in Stein & Shakarchi ch. 5, exercise 7, and we do not need it here.

■ **What just happened.** The transform is a unitary operator with four eigenvalues, and its eigenvectors are generated from the Gaussian by a "raising operator" $R$ that the transform turns into $-iR$. This is the same ladder-operator argument QM-3 uses for the harmonic oscillator, and the reason the oscillator's energy eigenstates look identical in position and momentum space.

### 3.9 The uncertainty principle

**Definition.** For $f\in L^2$ with $\|f\| = 1$, so that $|f(x)|^2$ and $|\hat f(\xi)|^2$ are probability densities (they have the same total mass by (3.15)), define the spreads
$$
\sigma_x^2 = \int x^2|f(x)|^2dx, \qquad \sigma_\xi^2 = \int\xi^2|\hat f(\xi)|^2d\xi.
$$
(These are variances about $0$. For a general $f$, shift it by (3.6) and modulate by (c) so both means are $0$; neither operation changes $|f|$ or $|\hat f|$.)

**Theorem 3.11 (Uncertainty principle).** If $\|f\|=1$, $xf$ and $f'$ are in $L^2$, and $x|f(x)|^2\to0$ as $|x|\to\infty$, then
$$
\sigma_x\,\sigma_\xi \;\ge\; \frac{1}{4\pi}, \tag{3.16}
$$
with equality if and only if $f(x) = Ce^{-\pi a x^2}$ for some $a>0$ and $|C|^2 = \sqrt{2a}$.

▶ **Try it first.** Two ingredients: (i) write $\sigma_\xi^2$ in terms of $\|f'\|$ using (3.8) and (3.15); (ii) integrate $\int|f|^2 = 1$ by parts against $x$ to produce $\int x\,\overline{f}\,f'$. Then Cauchy–Schwarz (LA-1 Theorem 3.11, in its LA-8 form for functions).

**Derivation.**

1. *Frequency spread as a derivative.* By (3.8), $\widehat{f'}(\xi) = 2\pi i\xi\hat f(\xi)$, so $|\widehat{f'}(\xi)|^2 = 4\pi^2\xi^2|\hat f(\xi)|^2$. By Plancherel (3.15) applied to $f'$: $\|f'\|^2 = \int|\widehat{f'}|^2 = 4\pi^2\sigma_\xi^2$. So $\sigma_\xi = \|f'\|/(2\pi)$.
2. *Integration by parts.* $1 = \int|f|^2dx = \big[x|f(x)|^2\big]_{-\infty}^{\infty} - \int x\,\frac{d}{dx}|f|^2\,dx$. The boundary term vanishes **by hypothesis** $x|f|^2\to0$. (For a non-decaying $f$ this step fails; see §8.) And $\frac{d}{dx}|f|^2 = \frac{d}{dx}(\bar f f) = \bar f' f + \bar f f' = 2\operatorname{Re}(\bar f f')$. So
$$
1 = -2\operatorname{Re}\int x\,\overline{f(x)}\,f'(x)\,dx = -2\operatorname{Re}\langle xf, f'\rangle.
$$
3. *Cauchy–Schwarz.* $1 = |{-2}\operatorname{Re}\langle xf, f'\rangle| \le 2|\langle xf, f'\rangle| \le 2\|xf\|\,\|f'\| = 2\sigma_x\cdot2\pi\sigma_\xi = 4\pi\sigma_x\sigma_\xi$. ∎
4. *Equality.* Both inequalities in step 3 must be equalities. Cauchy–Schwarz is an equality iff $f' = \lambda\,xf$ for a constant $\lambda\in\mathbb C$. Then $\langle xf, f'\rangle = \lambda\|xf\|^2$, and $|\operatorname{Re}(\cdot)| = |\cdot|$ with the sign from step 2 forces $\lambda$ real and negative: $\lambda = -2\pi a$, $a>0$. The ODE $f' = -2\pi a x f$ gives $f = Ce^{-\pi a x^2}$; normalizing, $|C|^2\int e^{-2\pi a x^2}dx = |C|^2/\sqrt{2a} = 1$. ∎

*Check on the Gaussian.* Take $a=1$: $f = 2^{1/4}e^{-\pi x^2}$. Then $\sigma_x^2 = \sqrt2\int x^2e^{-2\pi x^2}dx = \sqrt2\cdot\frac{\sqrt\pi}{2}(2\pi)^{-3/2} = \frac{1}{4\pi}$. By (3.11) $\hat f = 2^{1/4}e^{-\pi\xi^2}$ has the same shape, so $\sigma_\xi^2 = \frac{1}{4\pi}$ and $\sigma_x\sigma_\xi = \frac{1}{4\pi}$ exactly. The project reproduces $0.0795775$.

*Two more cases.* For $f = \sqrt a\,e^{-a|x|}$ (normalized two-sided exponential), exercise B3 computes $\sigma_x = 1/(a\sqrt2)$, $\sigma_\xi = a/(2\pi)$, product $\sqrt2/(4\pi)\approx0.1125$: $41\%$ above the minimum. For the rect $f = \Pi$, $\sigma_x^2 = 1/12$ but $\sigma_\xi^2 = \int\xi^2\operatorname{sinc}^2\xi\,d\xi = \int\frac{\sin^2\pi\xi}{\pi^2}d\xi = \infty$: a jump has infinite frequency spread (its derivative is a pair of deltas, $\|f'\| = \infty$ in step 1). On a finite frequency grid $|\xi|<L$ you get $\sigma_\xi^2\approx L/\pi^2$, growing without bound; the project shows it.

**QM-2 in one line.** A wavefunction $\psi(x)$ with $\|\psi\|=1$ has momentum-space wavefunction proportional to $\hat\psi$ with $p = \hbar k = 2\pi\hbar\xi$. So $\sigma_p = 2\pi\hbar\sigma_\xi$, and (3.16) reads $\sigma_x\sigma_p\ge2\pi\hbar/(4\pi) = \hbar/2$. **The Heisenberg uncertainty principle is Theorem 3.11.** There is no additional physics in it; the physics is in the claim that momentum is the Fourier variable of position.

■ **What just happened.** The derivative rule converts frequency spread to $\|f'\|$; integration by parts relates $\|f\|$ to the overlap of $xf$ and $f'$; Cauchy–Schwarz bounds that overlap. The same three moves, with $X$ and $P$ as abstract operators, give the general uncertainty relation in QM-4. Equality is the Gaussian because the Gaussian is the function whose derivative is proportional to $xf$.

### 3.10 Distributions: what $\delta$ is, what it is not, and its transform

**Definition 3.12 (Test functions and distributions).** A **test function** is a smooth $\varphi:\mathbb{R}\to\mathbb{C}$ such that $\varphi$ and all its derivatives decay faster than any power of $1/|x|$ (the **Schwartz space** $\mathcal S$; Gaussians times polynomials are the model examples). A **(tempered) distribution** is a linear functional $T:\mathcal S\to\mathbb C$, written $\varphi\mapsto\langle T,\varphi\rangle$, that is continuous (if $\varphi_n\to\varphi$ in $\mathcal S$ then $\langle T,\varphi_n\rangle\to\langle T,\varphi\rangle$). Every ordinary function $f$ of moderate growth defines a distribution by $\langle f,\varphi\rangle = \int f\varphi\,dx$; we identify $f$ with that functional.

**The Dirac delta** is the distribution $\langle\delta, \varphi\rangle = \varphi(0)$. It is linear ($\varphi\mapsto\varphi(0)$ is) and continuous. **The sifting property** $\int\delta(x)\varphi(x)\,dx = \varphi(0)$ is not a theorem about an integral; it is the *definition* written in the notation of an integral. Likewise $\delta_a$: $\langle\delta_a,\varphi\rangle = \varphi(a)$.

**Theorem 3.13 ($\delta$ is not a function).** There is no integrable function $d$ with $\int d(x)\varphi(x)\,dx = \varphi(0)$ for all test functions $\varphi$.

*Proof.* Suppose there were. Take any test function $\varphi$ that vanishes on an interval $(-\varepsilon,\varepsilon)$ around $0$; then $\int d\varphi = 0$. Since this holds for every such $\varphi$, $d = 0$ almost everywhere outside every $(-\varepsilon,\varepsilon)$, hence $d=0$ almost everywhere. But then $\int d\varphi = 0$ for *every* $\varphi$, including one with $\varphi(0) = 1$. Contradiction. ∎

So "$\delta(0) = \infty$, $\delta(x) = 0$ elsewhere" is a picture, not a definition, and in QM-2 it will mislead you (a plane wave is not normalizable; $\langle k|k'\rangle = \delta(k-k')$ is a statement about functionals, not a number). The honest picture is a **limit**:

**Theorem 3.14 ($\delta$ as a weak limit).** Let $g$ be integrable with $\int g = 1$, and $g_\varepsilon(x) = \frac1\varepsilon g(x/\varepsilon)$. Then $\langle g_\varepsilon,\varphi\rangle\to\varphi(0)$ for every test function $\varphi$. We write $g_\varepsilon\to\delta$ **weakly** (or "in the sense of distributions").

**Derivation.**
1. $\int\frac1\varepsilon g(x/\varepsilon)\varphi(x)\,dx = \int g(y)\varphi(\varepsilon y)\,dy$ (substitute $x = \varepsilon y$).
2. As $\varepsilon\to0$, $\varphi(\varepsilon y)\to\varphi(0)$ for each $y$; the integrand is bounded by $|g(y)|\sup|\varphi|$, integrable; dominated convergence gives $\int g(y)\varphi(0)\,dy = \varphi(0)$. ∎

Narrow Gaussians $\frac1\varepsilon e^{-\pi x^2/\varepsilon^2}$ and narrow rects $\frac1\varepsilon\Pi(x/\varepsilon)$ are both approximate deltas: same area, shrinking width. Nothing converges pointwise at $0$ (the heights blow up); only the *pairings* converge.

**Operations on distributions** are defined by what they would be if $T$ were a function, moved onto $\varphi$ by the corresponding identity:

- **Derivative.** For a function, $\int T'\varphi = -\int T\varphi'$ (integration by parts; $\varphi$ kills the boundary terms). So *define*
$$
\langle T',\varphi\rangle = -\langle T,\varphi'\rangle. \tag{3.17}
$$
Then $\langle\delta',\varphi\rangle = -\varphi'(0)$: $\delta'$ samples the negative slope. And the Heaviside step: $\langle u',\varphi\rangle = -\int_0^\infty\varphi'(x)\,dx = -[\varphi]_0^\infty = \varphi(0)$, so $u' = \delta$. Every distribution is infinitely differentiable.
- **Multiplication by a smooth function $m$.** $\langle mT,\varphi\rangle = \langle T, m\varphi\rangle$. So $\langle x\delta,\varphi\rangle = \langle\delta, x\varphi\rangle = 0\cdot\varphi(0) = 0$: $x\,\delta = 0$. And $m\,\delta_a = m(a)\,\delta_a$.
- **Scaling.** For a function, $\int T(ax)\varphi(x)dx = \frac1{|a|}\int T(y)\varphi(y/a)dy$. So $\langle\delta(a\,\cdot),\varphi\rangle = \frac1{|a|}\varphi(0)$: $\delta(ax) = \delta(x)/|a|$.
- **Convolution with a function.** $(f*\delta)(x) = \int f(y)\delta(x-y)dy = f(x)$ and $(f*\delta_a)(x) = f(x-a)$, by sifting. $\delta$ is the identity for convolution; $\delta_a$ is the delay operator.
- **Fourier transform.** For functions, $\int\hat f(\xi)\varphi(\xi)\,d\xi = \iint f(x)e^{-2\pi i x\xi}\varphi(\xi)\,dx\,d\xi = \int f(x)\hat\varphi(x)\,dx$. So *define*
$$
\langle\hat T,\varphi\rangle = \langle T,\hat\varphi\rangle. \tag{3.18}
$$
(The Schwartz space is exactly the class of $\varphi$ for which $\hat\varphi$ is again a test function; that is why it is the right class.) Every rule of Theorem 3.4 carries over, with the same derivations moved onto $\varphi$.

**Theorem 3.15.** $\hat\delta = 1$ and $\hat 1 = \delta$.

**Derivation (distributional).** $\langle\hat\delta,\varphi\rangle = \langle\delta,\hat\varphi\rangle = \hat\varphi(0) = \int\varphi(x)\,dx = \langle1,\varphi\rangle$. So $\hat\delta = 1$. And $\langle\hat1,\varphi\rangle = \langle1,\hat\varphi\rangle = \int\hat\varphi(\xi)\,d\xi = \varphi(0)$ by the inversion formula (3.2) at $x=0$; that is $\langle\delta,\varphi\rangle$. So $\hat1=\delta$. ∎

**Derivation (limiting).** By Theorem 3.14, $\delta$ is the weak limit of $g_\varepsilon(x) = \frac1\varepsilon e^{-\pi x^2/\varepsilon^2}$. By (3.11) and the stretch theorem, $\hat g_\varepsilon(\xi) = e^{-\pi\varepsilon^2\xi^2}\to1$ pointwise and (dominated) weakly. The transform is continuous on distributions, so $\hat\delta = \lim\hat g_\varepsilon = 1$. Conversely, $1$ is the weak limit of the wide Gaussians $e^{-\pi\varepsilon^2x^2}$, whose transforms $\frac1\varepsilon e^{-\pi\xi^2/\varepsilon^2}$ are narrow Gaussians of area $1$, which tend to $\delta$. ∎

$$
\hat\delta = 1,\qquad \hat1 = \delta,\qquad \widehat{\delta_a}(\xi) = e^{-2\pi i a\xi},\qquad \widehat{e^{2\pi i\nu x}} = \delta_\nu,\qquad \widehat{\cos2\pi\nu x} = \tfrac12(\delta_\nu+\delta_{-\nu}). \tag{3.19}
$$
(The last three follow from the first two by shift and modulation.) A pure tone has a spectrum that is a single spike. $\hat f(0)=\int f$ in the form $\hat 1 = \delta$: the constant has infinite total mass, all at frequency $0$.

**The comb and Poisson summation.** $\mathrm{III}(x) = \sum_{n\in\mathbb Z}\delta(x-n)$, so $\langle\mathrm{III},\varphi\rangle = \sum_n\varphi(n)$ (converges because $\varphi$ decays fast).

**Theorem 3.16 (Poisson summation).** For a test function $\varphi$,
$$
\sum_{n\in\mathbb Z}\varphi(n) = \sum_{k\in\mathbb Z}\hat\varphi(k). \tag{3.20}
$$

▶ **Try it first.** Periodize $\varphi$: $\Phi(x) = \sum_n\varphi(x+n)$ has period $1$. Compute its F2 Fourier coefficients. Then evaluate the Fourier series at $x=0$.

**Derivation.**
1. $\Phi(x) = \sum_n\varphi(x+n)$ converges absolutely and is smooth with period $1$, so F2 applies: $\Phi(x) = \sum_k c_k e^{2\pi i kx}$, $c_k = \int_0^1\Phi(x)e^{-2\pi i kx}dx$.
2. $c_k = \sum_n\int_0^1\varphi(x+n)e^{-2\pi i kx}dx = \sum_n\int_n^{n+1}\varphi(y)e^{-2\pi i k(y-n)}dy$ (substitute $y = x+n$).
3. $e^{-2\pi i k(y-n)} = e^{-2\pi i ky}e^{2\pi i kn} = e^{-2\pi i ky}$ because $e^{2\pi i kn} = 1$. The intervals $[n,n+1]$ tile the line, so $c_k = \int_{-\infty}^\infty\varphi(y)e^{-2\pi i ky}dy = \hat\varphi(k)$.
4. The Fourier series of a smooth periodic function converges at every point (F2), so $\Phi(0) = \sum_k\hat\varphi(k)$. And $\Phi(0) = \sum_n\varphi(n)$. ∎

**Corollary (the comb is its own transform).** $\widehat{\mathrm{III}} = \mathrm{III}$. *Proof:* $\langle\widehat{\mathrm{III}},\varphi\rangle = \langle\mathrm{III},\hat\varphi\rangle = \sum_k\hat\varphi(k) = \sum_n\varphi(n) = \langle\mathrm{III},\varphi\rangle$ by (3.20). ∎

Scaled: define $\mathrm{III}_p(x) = \sum_n\delta(x-np)$, spikes every $p$. Since $\delta(x/p - n) = \delta((x-np)/p) = p\,\delta(x-np)$, we have $\mathrm{III}_p(x) = \frac1p\mathrm{III}(x/p)$, and by the stretch theorem (3.7),
$$
\widehat{\mathrm{III}_p}(\xi) = \mathrm{III}(p\xi) = \frac1p\,\mathrm{III}_{1/p}(\xi) = \frac1p\sum_k\delta\!\left(\xi-\frac kp\right). \tag{3.21}
$$
Spikes every $p$ in $x$ become spikes every $1/p$ in $\xi$, with height $1/p$.

**The F2 link, made explicit.** Periodizing $f$ with period $p$ is convolving with the comb: $(f*\mathrm{III}_p)(x) = \sum_n f(x-np)$. By (3.13) and (3.21), its transform is $\hat f(\xi)\cdot\frac1p\mathrm{III}_{1/p}(\xi) = \frac1p\sum_k\hat f(k/p)\,\delta(\xi-k/p)$: a spike at each harmonic $k/p$ with weight $\frac1p\hat f(k/p)$. Those weights are exactly the F2 coefficients $c_k$ of the periodized function. Fourier series are the Fourier transform of a periodic function; the transform of a periodic function is a comb of deltas weighted by the series coefficients.

**The Heaviside step.** Write $u = \frac12 + \frac12\operatorname{sgn}(x)$. By (3.4) and (3.5)'s computation with the sign flipped, $\widehat{\operatorname{sgn}(x)e^{-a|x|}} = \frac{1}{a+2\pi i\xi} - \frac{1}{a-2\pi i\xi} = \frac{-4\pi i\xi}{a^2+4\pi^2\xi^2}$. As $a\to0^+$ this tends to $\frac{-4\pi i\xi}{4\pi^2\xi^2} = \frac{1}{\pi i\xi}$ for $\xi\ne0$; the limit holds in the sense of distributions with the singularity at $0$ interpreted as a **principal value** $\mathrm{PV}\frac{1}{\pi i\xi}$ (the odd singularity's two sides cancel when paired with a test function; Osgood ch. 4 does the bookkeeping). Hence
$$
\hat u(\xi) = \frac12\delta(\xi) + \mathrm{PV}\,\frac{1}{2\pi i\xi}. \tag{3.22}
$$
The $\delta/2$ is the mean value $\tfrac12$ of $u$; the $1/(2\pi i\xi)$ is the jump, decaying like $1/\xi$ as every jump does.

■ **What just happened.** A distribution is a rule for producing numbers from test functions; $\delta$ is "evaluate at $0$." Every operation (derivative, scaling, transform) is defined by moving it onto the test function via the identity that holds for ordinary functions. With that, $\hat\delta=1$ is a one-line computation, and Poisson summation (which is F2's coefficient formula applied to a periodized function) says the comb is its own transform. §3.15 turns (3.21) into the sampling theorem.

### 3.11 Linear time-invariant systems and the damped driven oscillator

A **system** $L$ maps an input signal $x(t)$ to an output $y(t) = L[x](t)$. It is **linear** if $L[c_1x_1+c_2x_2] = c_1L[x_1]+c_2L[x_2]$ and **time-invariant** if delaying the input delays the output: $L[x(\cdot-\tau)](t) = y(t-\tau)$. The **impulse response** is $h = L[\delta]$.

**Theorem 3.17 (LTI systems are convolutions).** If $L$ is linear, time-invariant, and continuous (so that it commutes with the integral below), then $L[x] = h*x$ for every input $x$.

**Derivation.**
1. Write the input as a superposition of shifted impulses: $x(t) = (x*\delta)(t) = \int x(\tau)\,\delta(t-\tau)\,d\tau$.
2. By time invariance, $L[\delta(\cdot-\tau)] = h(\cdot-\tau)$.
3. By linearity (extended to the integral by continuity), $L[x](t) = \int x(\tau)\,L[\delta(\cdot-\tau)](t)\,d\tau = \int x(\tau)h(t-\tau)\,d\tau = (h*x)(t)$. ∎

Taking transforms with (3.13): $\hat y(\xi) = H(\xi)\hat x(\xi)$, with $H = \hat h$ the **transfer function** (or frequency response). Feed in a pure tone $x = e^{2\pi i\nu t}$: $y = h*x = \int h(\tau)e^{2\pi i\nu(t-\tau)}d\tau = H(\nu)e^{2\pi i\nu t}$. Pure tones are **eigenfunctions** of every LTI system, with eigenvalue $H(\nu)$: a sinusoid in gives the same sinusoid out, scaled by $|H(\nu)|$ and delayed in phase by $\arg H(\nu)$. This is the frequency-domain picture of LA-7: the exponentials diagonalize all shift-invariant operators at once. Theorem 3.22 is the finite version.

**The damped driven oscillator (ME vibration analysis).** Mass $m$, damping $c$, stiffness $k$, forcing $f(t)$:
$$
m\,x'' + c\,x' + k\,x = f(t).
$$
Transform both sides with (3.8) twice, writing $\omega = 2\pi\xi$ so that $\widehat{x'} = i\omega\hat x$ and $\widehat{x''} = -\omega^2\hat x$: $(-m\omega^2 + ic\omega + k)\hat x = \hat f$. So $\hat x = H\hat f$ with
$$
H(\omega) = \frac{1}{k - m\omega^2 + ic\omega}, \qquad |H(\omega)|^2 = \frac{1}{(k-m\omega^2)^2 + c^2\omega^2}. \tag{3.23}
$$
(Keeping $\omega$ here matches every vibration text; $H$ as a function of $\xi$ is $H(2\pi\xi)$.) Three regimes, which together make the **Bode plot** shape (log $|H|$ against log $\omega$):

- $\omega\to0$: $|H|\to1/k$, the static deflection per unit force. Flat.
- $\omega\to\infty$: $|H|\approx1/(m\omega^2)$. A straight line of slope $-2$ on log-log axes ($-40$ dB per decade): the mass cannot follow fast forcing.
- Near $\omega_0 = \sqrt{k/m}$: **resonance.** Maximize $|H|^2$ by minimizing the denominator $D(\omega^2) = (k-m\omega^2)^2 + c^2\omega^2$ in the variable $s = \omega^2$: $D'(s) = -2m(k-ms)+c^2 = 0$ gives
$$
\omega_{\text{peak}}^2 = \frac km - \frac{c^2}{2m^2}, \qquad |H|_{\max} = \frac{1}{c\sqrt{k/m - c^2/(4m^2)}} \tag{3.24}
$$
(exists when $c^2<2mk$; substitute $k - m\omega_{\text{peak}}^2 = c^2/(2m)$ into $D$ to get the height). Lighter damping: higher, narrower peak, closer to $\omega_0$. At $\omega=\omega_0$ exactly, $H = 1/(ic\omega_0)$: amplitude $1/(c\omega_0)$ and phase $-90°$, the response lagging the force by a quarter cycle. The phase of $H$ runs from $0$ (below resonance, in phase) through $-90°$ to $-180°$ (above resonance, opposite).

**The impulse response.** Integrate the ODE across $t=0$ with $f = \delta$: the $kx$ and $cx'$ terms contribute nothing over an interval of length $\to0$ (they are bounded), so $m[x'(0^+)-x'(0^-)] = 1$. Starting from rest, the impulse sets $x(0) = 0$, $x'(0^+) = 1/m$, and for $t>0$ the motion is free: $x = e^{-ct/(2m)}(A\cos\omega_d t + B\sin\omega_d t)$ with $\omega_d = \sqrt{k/m - c^2/(4m^2)}$ (the characteristic roots of $m r^2+cr+k=0$ are $-c/(2m)\pm i\omega_d$). The initial conditions give $A = 0$, $B = 1/(m\omega_d)$:
$$
h(t) = \frac{1}{m\omega_d}\,e^{-ct/(2m)}\sin(\omega_d t)\,u(t). \tag{3.25}
$$
Its transform must be (3.23); the project checks this numerically with `solve_ivp` and an FFT.

■ **What just happened.** An ODE with constant coefficients is an LTI system. Transforming turns $d/dt$ into $i\omega$ and the ODE into division by a polynomial. The transfer function (3.23) is the oscillator's whole personality: static stiffness, resonance, and mass-dominated roll-off, read off a single rational function.

### 3.12 The discrete Fourier transform

**Definition 3.18 (DFT).** For $f = (f[0],\dots,f[N-1])\in\mathbb C^N$, with $\omega = e^{-2\pi i/N}$,
$$
F[k] = \sum_{n=0}^{N-1} f[n]\,\omega^{kn} = \sum_{n=0}^{N-1}f[n]\,e^{-2\pi i kn/N}, \qquad k = 0,\dots,N-1. \tag{3.26}
$$
In matrix form $F = \mathbf F f$ where $\mathbf F$ is the $N\times N$ **Fourier matrix** with entries $\mathbf F_{kn} = \omega^{kn}$. (We write $\mathbf F$ bold for the matrix and $F[k]$ for the output vector.) This is NumPy's `np.fft.fft`: no $1/N$, no $1/\sqrt N$. Because $\omega^N = 1$, both $f$ and $F$ are naturally **periodic with period $N$**: $F[k+N] = F[k]$, and we read $f[n]$ with $n$ taken mod $N$.

Compare with (3.1): the integral became a sum, $x = n$, $\xi = k/N$. And compare with F2: $F[k] = \langle e_k, f\rangle$ where $e_k[n] = e^{2\pi i kn/N}$ is the $k$-th **discrete tone** and $\langle\cdot,\cdot\rangle$ is the $\mathbb C^N$ inner product with the conjugate on the first slot. The DFT is the inner product with the tones, exactly as the Fourier coefficient was.

**Theorem 3.19 ($\mathbf F/\sqrt N$ is unitary; inverse DFT).** The columns of $\mathbf F$ are mutually orthogonal, each of squared length $N$: $\mathbf F^*\mathbf F = N I$. Hence $\mathbf F/\sqrt N$ is unitary, $\mathbf F^{-1} = \overline{\mathbf F}/N$, and
$$
f[n] = \frac1N\sum_{k=0}^{N-1}F[k]\,e^{2\pi i kn/N}. \tag{3.27}
$$

▶ **Try it first.** The inner product of columns $m$ and $n$ is $\sum_k\overline{\omega^{km}}\omega^{kn} = \sum_k\omega^{k(n-m)}$. That is a geometric series in the root of unity $\omega^{n-m}$. Sum it (F1).

**Derivation.**
1. $(\mathbf F^*\mathbf F)_{mn} = \sum_{k=0}^{N-1}\overline{\omega^{km}}\,\omega^{kn} = \sum_{k=0}^{N-1}\omega^{k(n-m)}$, since $\bar\omega = \omega^{-1}$.
2. If $n = m$, every term is $1$ and the sum is $N$.
3. If $n\ne m$, put $r = \omega^{n-m}$. Then $r\ne1$ (because $0<|n-m|<N$ and $\omega$ is a *primitive* $N$-th root of unity) and $r^N = \omega^{N(n-m)} = 1$. Geometric sum: $\sum_{k=0}^{N-1}r^k = \frac{1-r^N}{1-r} = 0$.
4. So $\mathbf F^*\mathbf F = NI$, i.e. $(\mathbf F/\sqrt N)^*(\mathbf F/\sqrt N) = I$. A square matrix with orthonormal columns is unitary. $\mathbf F$ is symmetric ($\omega^{kn} = \omega^{nk}$), so $\mathbf F^* = \overline{\mathbf F}$ and $\mathbf F^{-1} = \overline{\mathbf F}/N$. Reading off the entries of $f = \overline{\mathbf F}F/N$ gives (3.27). ∎

Consequences, each immediate: **Parseval** $\sum_n|f[n]|^2 = \frac1N\sum_k|F[k]|^2$ (unitary matrices preserve length; the $1/N$ is the missing normalization). **Duality** $\mathbf F^2 = N\,P$ where $P$ is the reversal $n\mapsto-n\bmod N$ (exercise C2), so $\mathbf F^4 = N^2I$ and the eigenvalues of $\mathbf F/\sqrt N$ are fourth roots of unity, as in Theorem 3.9. **Hermitian symmetry** for real $f$: $F[N-k] = \overline{F[k]}$.

**Which frequency is index $k$?** The tone $e_k[n] = e^{2\pi i kn/N}$ completes $k$ cycles in $N$ samples: $k$ is a frequency in **cycles per record**. Since $e_{k}=e_{k-N}$, indices $k > N/2$ are the *negative* frequencies $k-N$. NumPy stores them in the second half of the array; `np.fft.fftshift` moves them to the left where a plot wants them, and `np.fft.fftfreq(N, d=Δ)` gives the physical frequencies. If the samples are $\Delta$ apart in physical units, index $k$ means
$$
\xi_k = \frac{k}{N\Delta}\ \text{cycles per unit},\qquad k = -\tfrac N2,\dots,\tfrac N2-1, \tag{3.28}
$$
a grid of spacing $1/(N\Delta) = 1/(\text{record length})$ extending to $\pm1/(2\Delta)$, the **Nyquist frequency**.

**Theorem 3.20 (The DFT approximates samples of $\hat f$).** If $f$ is negligible outside $[0, N\Delta)$ and $f[n] = f(n\Delta)$, then
$$
\hat f\!\left(\frac{k}{N\Delta}\right)\approx\Delta\,F[k]. \tag{3.29}
$$
*Derivation.* Riemann sum of (3.1) with spacing $\Delta$: $\hat f(\xi)\approx\sum_n f(n\Delta)e^{-2\pi i n\Delta\xi}\Delta$. At $\xi = k/(N\Delta)$ the exponent is $-2\pi i nk/N$, and the sum is $F[k]$. ∎ Two errors hide in "$\approx$": the grid $\Delta$ (which produces aliasing, §3.15: the DFT actually computes samples of the *periodized* spectrum) and the finite record (which produces leakage, §3.16). If the record starts at $x_0$ instead of $0$, multiply by $e^{-2\pi i x_0\xi_k}$ by the shift rule; the project does this.

**Definition 3.21 (Circular convolution).** $(f\circledast g)[n] = \sum_{m=0}^{N-1}f[m]\,g[(n-m)\bmod N]$.

**Theorem 3.22 (Circular convolution theorem).** $\mathrm{DFT}(f\circledast g)[k] = F[k]\,G[k]$.

**Derivation.**
1. $\sum_n(f\circledast g)[n]\omega^{kn} = \sum_n\sum_m f[m]g[n-m]\omega^{kn}$ (indices mod $N$).
2. Write $\omega^{kn} = \omega^{km}\omega^{k(n-m)}$; this is legitimate mod $N$ because $\omega^{kN}=1$, so $\omega^{k(n-m)}$ depends only on $(n-m)\bmod N$.
3. Swap sums: $= \sum_m f[m]\omega^{km}\sum_n g[n-m]\omega^{k(n-m)}$. As $n$ runs over $0,\dots,N-1$, so does $j = (n-m)\bmod N$; the inner sum is $\sum_j g[j]\omega^{kj} = G[k]$.
4. $= G[k]\sum_m f[m]\omega^{km} = F[k]G[k]$. ∎

**Linear convolution by zero padding.** If $f$ is nonzero only on $0\le n<p$ and $g$ only on $0\le n<q$, their ordinary (linear) convolution has length $p+q-1$. Pad both with zeros to length $N\ge p+q-1$; then the circular convolution cannot wrap around and equals the linear one. So linear convolution $=$ IFFT(FFT$(f)\cdot$FFT$(g))$ in $O(N\log N)$.

### 3.13 Circulant matrices: convolution is a matrix, and $\mathbf F$ diagonalizes it

**Definition 3.23.** The **cyclic shift** $S$ on $\mathbb C^N$ is $(Sx)[n] = x[(n-1)\bmod N]$ (every entry moves down one slot, the last wraps to the top). A **circulant matrix** with first column $c$ is $C = \sum_{m=0}^{N-1}c[m]\,S^m$; equivalently $C_{jk} = c[(j-k)\bmod N]$: each column is the previous one shifted down by one.

Since $(S^mx)[n] = x[n-m]$, we get $(Cx)[n] = \sum_m c[m]\,x[n-m] = (c\circledast x)[n]$: **circular convolution with $c$ is multiplication by the circulant matrix $C$.** For $N=4$ and $c = (c_0,c_1,c_2,c_3)$:
$$
C = \begin{pmatrix}c_0&c_3&c_2&c_1\\c_1&c_0&c_3&c_2\\c_2&c_1&c_0&c_3\\c_3&c_2&c_1&c_0\end{pmatrix}.
$$

**Theorem 3.24 (Every circulant is diagonalized by the Fourier matrix).** The discrete tones $e_k$ (the columns of $\overline{\mathbf F}$, i.e. of $N\mathbf F^{-1}$) are eigenvectors of $S$ with eigenvalues $\omega^k$, and of every circulant $C$ with eigenvalues $\hat c[k] = \mathrm{DFT}(c)[k]$. Hence
$$
C = \mathbf F^{-1}\,\Lambda\,\mathbf F,\qquad \Lambda = \operatorname{diag}(\mathrm{DFT}(c)). \tag{3.30}
$$

▶ **Try it first.** Apply $S$ to $e_k[n] = e^{2\pi i kn/N}$ and factor out $e_k[n]$.

**Derivation.**
1. $(Se_k)[n] = e_k[n-1] = e^{2\pi i k(n-1)/N} = e^{-2\pi i k/N}e_k[n] = \omega^k e_k[n]$. So $Se_k = \omega^ke_k$.
2. $Ce_k = \sum_mc[m]S^me_k = \sum_mc[m]\omega^{km}e_k = \hat c[k]\,e_k$: a polynomial in $S$ has the same eigenvectors as $S$, with the polynomial applied to the eigenvalues.
3. The $e_k$ are the columns of $\overline{\mathbf F}$, which is $N\mathbf F^{-1}$; collecting the eigenvector equations: $C\,\overline{\mathbf F} = \overline{\mathbf F}\Lambda$, so $C = \overline{\mathbf F}\Lambda\overline{\mathbf F}^{-1} = \mathbf F^{-1}\Lambda\mathbf F$ (the factors $N$ cancel). ∎

**Theorem 3.25 (Shift-invariant $\Leftrightarrow$ circulant).** A matrix $C$ commutes with $S$ if and only if it is circulant.

*Proof.* $(SC)_{jk} = C_{j-1,k}$ (shifting rows) and $(CS)_{jk} = \sum_mC_{jm}S_{mk} = C_{j,k+1}$ (since $S_{mk} = 1$ iff $m = k+1$). So $SC = CS$ iff $C_{j-1,k} = C_{j,k+1}$ for all $j,k$, i.e. $C_{j,k} = C_{j+1,k+1}$: constant along wrapped diagonals, which is the definition of circulant. Conversely a polynomial in $S$ commutes with $S$. ∎

This is Theorem 3.17 in finite dimensions: a linear operator that commutes with shifts is a convolution, and (3.30) says the tones diagonalize it. Numerically, (3.30) is how convolution is done fast: $Cx = \mathbf F^{-1}(\Lambda(\mathbf Fx))$ costs two FFTs and $N$ multiplications instead of $N^2$. **ML-8 preview:** a 1D convolution layer with kernel $c$ and circular padding *is* the matrix $C$ (with zero padding it is the closely related banded Toeplitz matrix); weight sharing is the statement that the matrix is circulant, and its spectrum is the DFT of the kernel. Strang's 18.065 lectures 31–32 are built around this.

### 3.14 The fast Fourier transform

Computing (3.26) directly costs $N^2$ complex multiplications. The FFT exploits the structure of $\omega$.

▶ **Try it first.** Split the sum (3.26) into even $n = 2m$ and odd $n = 2m+1$. Each half should look like a DFT of length $N/2$.

**Theorem 3.26 (Cooley–Tukey radix-2 step).** Let $N = 2M$. Let $E = \mathrm{DFT}_M(f[0], f[2], \dots, f[N-2])$ and $O = \mathrm{DFT}_M(f[1], f[3],\dots,f[N-1])$. Then for $k = 0,\dots,M-1$:
$$
F[k] = E[k] + \omega_N^k\,O[k],\qquad F[k+M] = E[k] - \omega_N^k\,O[k], \tag{3.31}
$$
where $\omega_N = e^{-2\pi i/N}$ and the factors $\omega_N^k$ are the **twiddle factors**.

**Derivation.**
1. $F[k] = \sum_{m=0}^{M-1}f[2m]\,\omega_N^{2mk} + \sum_{m=0}^{M-1}f[2m+1]\,\omega_N^{(2m+1)k}$.
2. $\omega_N^2 = e^{-2\pi i\cdot2/N} = e^{-2\pi i/M} = \omega_M$. So the first sum is $\sum_mf[2m]\omega_M^{mk} = E[k]$ and the second is $\omega_N^k\sum_mf[2m+1]\omega_M^{mk} = \omega_N^kO[k]$. This holds for all $k$, with $E$ and $O$ understood as $M$-periodic.
3. For the upper half, $E[k+M] = E[k]$ and $O[k+M] = O[k]$ by periodicity, while $\omega_N^{k+M} = \omega_N^k\,\omega_N^{N/2} = \omega_N^k e^{-\pi i} = -\omega_N^k$. ∎

Each pair $(E[k], O[k])\mapsto(F[k], F[k+M])$ is a **butterfly**: one complex multiplication ($\omega_N^kO[k]$) and two additions. If $N = 2^m$, recurse on $E$ and $O$ down to length $1$ (where the DFT is the identity). Worked example 4.6 does $N=8$ by hand.

**Cost.** Let $T(N)$ be the number of complex multiplications. The combine step needs $N/2$ (one per butterfly), so $T(N) = 2T(N/2) + N/2$ with $T(1) = 0$. Unroll: $T(N) = 2^jT(N/2^j) + j\cdot N/2$. At $j = \log_2N$, $T(N/2^j) = T(1) = 0$:
$$
T(N) = \frac N2\log_2N\quad\text{complex multiplications (and $N\log_2N$ additions), versus $N^2$ directly.} \tag{3.32}
$$
For $N = 2^{10} = 1024$: $5120$ vs $1{,}048{,}576$, a factor of $200$. For $N = 2^{20}$: a factor of $10^5$. (Many twiddles are trivial, $\pm1,\pm i$, so real implementations do even better.) The project times both.

**Why powers of two.** The radix-2 split needs $N$ even at every level. For $N = 3\cdot4$ Gauss split into $3$ groups of $4$; general mixed-radix FFTs exist for any $N$ with small prime factors (NumPy's `fft` handles all $N$, fastest for highly composite $N$). A prime $N$ has no split and needs other tricks (Bluestein, Rader). In practice: **zero-pad** to the next power of two. Zero padding does not add information; it evaluates the same trigonometric polynomial on a finer frequency grid (3.28), which is useful for locating peaks (the project's oscillator section needs it) but does not improve resolution, which is set by the record length.

### 3.15 Sampling: the Nyquist–Shannon theorem and aliasing

Sample $f$ every $\Delta$. As a distribution the sampled signal is
$$
f_s(x) = f(x)\,\mathrm{III}_\Delta(x) = \sum_nf(n\Delta)\,\delta(x-n\Delta).
$$

**Theorem 3.27 (Sampling periodizes the spectrum).**
$$
\hat f_s(\xi) = \frac1\Delta\sum_{k\in\mathbb Z}\hat f\!\left(\xi-\frac k\Delta\right). \tag{3.33}
$$
**Derivation.** By the mirror convolution theorem (3.14), $\hat f_s = \hat f*\widehat{\mathrm{III}_\Delta}$. By (3.21), $\widehat{\mathrm{III}_\Delta} = \frac1\Delta\sum_k\delta_{k/\Delta}$. Convolution with $\delta_a$ is a shift: $\hat f*\delta_{k/\Delta} = \hat f(\xi-k/\Delta)$. Sum. ∎

The spectrum of the samples is the original spectrum **repeated every $1/\Delta$** (the sampling rate), scaled by $1/\Delta$. This is the §2 picture, and (3.29) is the statement that the DFT computes samples of (3.33), not of $\hat f$.

**Theorem 3.28 (Nyquist–Shannon sampling theorem).** Suppose $f$ is **band-limited**: $\hat f(\xi) = 0$ for $|\xi|\ge B$. If the sampling rate exceeds twice the bandwidth, $1/\Delta\ge2B$, then $f$ is determined by its samples, and
$$
f(x) = \sum_{n\in\mathbb Z}f(n\Delta)\,\operatorname{sinc}\!\left(\frac{x-n\Delta}{\Delta}\right). \tag{3.34}
$$

▶ **Try it first.** Draw the copies in (3.33). Under the hypothesis they do not overlap. Multiply by a rect that keeps only the $k=0$ copy, and inverse-transform.

**Derivation.**
1. The $k$-th copy in (3.33) occupies $(k/\Delta-B,\,k/\Delta+B)$. The $k=0$ copy is inside $(-1/(2\Delta), 1/(2\Delta))$ and the $k=\pm1$ copies start at $\pm(1/\Delta-B)$, which is $\ge\pm1/(2\Delta)$ exactly when $B\le1/(2\Delta)$. So the copies do not overlap, and on $|\xi|<1/(2\Delta)$ only the $k=0$ term survives.
2. Therefore $\hat f(\xi) = \Delta\,\Pi(\Delta\xi)\,\hat f_s(\xi)$: the rect of width $1/\Delta$ (which equals $1$ for $|\xi|<1/(2\Delta)$) cuts out the central copy, and $\Delta$ undoes the $1/\Delta$.
3. Inverse-transform with the convolution theorem (3.13) read backwards: $f = \mathcal F^{-1}[\Delta\Pi(\Delta\xi)]*f_s$. By (3.3), duality, and the stretch theorem, $\mathcal F^{-1}[\Pi(\Delta\xi)](x) = \frac1\Delta\operatorname{sinc}(x/\Delta)$, so $\mathcal F^{-1}[\Delta\Pi(\Delta\xi)] = \operatorname{sinc}(x/\Delta)$.
4. $f(x) = \big(\operatorname{sinc}(\cdot/\Delta)*f_s\big)(x) = \sum_nf(n\Delta)\operatorname{sinc}\!\left(\frac{x-n\Delta}{\Delta}\right)$, using $g*\delta_a = g(\cdot-a)$. ∎

$2B$ is the **Nyquist rate**; $1/(2\Delta)$ is the **Nyquist frequency** of a given sampler. The interpolating function $\operatorname{sinc}((x-n\Delta)/\Delta)$ equals $1$ at $x = n\Delta$ and $0$ at every other sample point, so (3.34) passes through the samples; what the theorem adds is that it is the *only* band-limited function that does.

**Aliasing.** If $1/\Delta<2B$ the copies overlap and (3.33) mixes them irreversibly. For a single tone this is explicit: the samples of $e^{2\pi i\nu x}$ are $e^{2\pi i\nu n\Delta}$, and
$$
e^{2\pi i\nu n\Delta} = e^{2\pi i(\nu - k/\Delta)n\Delta}\quad\text{for every integer }k, \tag{3.35}
$$
because $e^{-2\pi i kn} = 1$. The sampler cannot distinguish $\nu$ from $\nu-k/\Delta$. The frequency it will report is the **alias**: the unique $\nu' \equiv \nu\pmod{1/\Delta}$ in $(-1/(2\Delta), 1/(2\Delta)]$; for a real sinusoid the apparent frequency is $|\nu'|$. **Folding formula:** $\nu_{\text{alias}} = \left|\nu - \tfrac1\Delta\operatorname{round}(\nu\Delta)\right|$. Frequencies above Nyquist fold back like a sheet of paper creased at $1/(2\Delta)$: $\nu = 1/(2\Delta)+\epsilon$ appears at $1/(2\Delta)-\epsilon$.

*The wagon wheel.* Film at $10$ frames per second a wheel turning at $7$ revolutions per second. Between frames it turns $0.7$ of a turn, which looks identical to turning $-0.3$ of a turn. On screen it turns backwards at $3$ rev/s. A stroboscope is the same effect used deliberately: flash at the rotation rate and the shaft appears to stand still (alias $0$). Worked example 4.7 does the arithmetic. **Nothing in the sampled data warns you**; the samples of a $7$ Hz sine at $10$ Hz are, to the last digit, the samples of a $3$ Hz sine with a sign flip. The only defenses are knowing the bandwidth or filtering before sampling (an **anti-aliasing filter**, which every ADC has in front of it for exactly this reason).

### 3.16 Windowing and spectral leakage

A finite record of length $T$ is $f(x)\Pi(x/T)$, so by (3.14) its spectrum is $\hat f*T\operatorname{sinc}(T\xi)$: every spectral line is smeared into a sinc of main-lobe width $2/T$ with sidelobes decaying like $1/\xi$ (the first sidelobe is $13$ dB down). For the DFT this means: a tone at frequency $\nu$ gives a single nonzero bin only if $\nu T$ is an integer (an integer number of cycles in the record), because then the zeros of the sinc land exactly on the other grid points $k/T$. Otherwise its energy **leaks** into every bin (worked example 4.5). A **window** $w$ with gentler edges, such as the **Hann window** $w(x) = \cos^2(\pi x/T) = \frac12(1+\cos(2\pi x/T))$ on $|x|<T/2$, has a transform that is three sincs ($\frac12\operatorname{sinc}$ at $0$ and $\frac14\operatorname{sinc}$ at $\pm1/T$, by the modulation rule) whose sidelobes nearly cancel: first sidelobe $31$ dB down, at the cost of a main lobe twice as wide. The trade is always sidelobe level against resolution; exercise C3 measures it.

### 3.17 Two dimensions

For $f(x,y)$, $\hat f(\xi,\eta) = \iint f(x,y)e^{-2\pi i(x\xi+y\eta)}dx\,dy$: the kernel factors as $e^{-2\pi i x\xi}e^{-2\pi i y\eta}$, so the 2D transform is a 1D transform of each row followed by a 1D transform of each column (**separable**), and a separable function $g(x)h(y)$ transforms to $\hat g(\xi)\hat h(\eta)$. Everything in §3.4–3.6 holds with vectors in place of scalars; in particular convolution is still multiplication in frequency, so **blurring an image** (convolving with a Gaussian or box kernel) is multiplying its spectrum by a low-pass filter that suppresses high spatial frequencies, which is where edges live. The physics link: in **Fraunhofer (far-field) diffraction**, EM-7 shows that the field at angle $\theta$ from an aperture is the Fourier transform of the aperture function evaluated at $\xi = \sin\theta/\lambda$. A slit of width $a$ (a rect) gives the intensity $\operatorname{sinc}^2(a\sin\theta/\lambda)$ and a circular aperture gives the Airy pattern; the diffraction pattern *is* the transform, and narrow apertures give wide patterns by the stretch theorem.

### 3.18 Time–frequency: the short-time Fourier transform

To see how the frequency content of a signal changes with time, slide a window $w$ of length $T_w$ along it and transform each piece: $S(\tau,\xi) = \int f(x)\,w(x-\tau)\,e^{-2\pi i x\xi}dx$. The plot of $|S(\tau,\xi)|^2$ is the **spectrogram**. Each column is the spectrum of $f\cdot w(\cdot-\tau)$, which by (3.14) is $\hat f$ smeared by $\hat w$ (width $\sim1/T_w$). So a short window locates events in time to within $T_w$ but blurs frequencies by $1/T_w$; a long window does the reverse. You cannot win both, and the reason is Theorem 3.11 applied to the window: $\sigma_\tau\sigma_\xi\ge1/(4\pi)$, with the Gaussian window as the best compromise. The project builds a spectrogram of a chirp with two window lengths and you will see the trade-off directly.
