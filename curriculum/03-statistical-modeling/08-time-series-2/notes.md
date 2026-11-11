# ST-8: Time Series II: Spectral Analysis, State-Space Models and the Kalman Filter, Volatility Models

**Time:** ~12 hours · **Prerequisites:** ST-7 (ARMA, ACF, Yule–Walker, Ljung–Box, exponential smoothing); F3 (Fourier transform, DFT, convolution theorem, windowing); F5 (MVN conditionals); LA-10 (Cholesky) · **Unlocks:** ML-4 (Gaussian processes as state-space models); control-theory pointers for an ME

> **Conventions in this module.** Frequency $\omega$ is in **cycles per sample**, $\omega\in[-\tfrac12,\tfrac12]$ (F3's convention); $\lambda=2\pi\omega$ is the angular frequency in radians per sample when it shortens a formula. The **spectral density** is $f(\omega)=\sum_h\gamma(h)e^{-2\pi i\omega h}$, so $\gamma(h)=\int_{-1/2}^{1/2}f(\omega)e^{2\pi i\omega h}\,d\omega$ and $\int f=\gamma(0)$: no factors of $2\pi$ anywhere. A filter $y_t=\sum_j\psi_jx_{t-j}$ has transfer polynomial $\psi(z)=\sum_j\psi_jz^j$ and **frequency response** $\psi(e^{-2\pi i\omega})$. The **periodogram** is $I(\omega_j)=|\sum_{t=1}^nx_te^{-2\pi i\omega_jt}|^2/n$ at the Fourier frequencies $\omega_j=j/n$. ARMA follows ST-7: $\phi(B)x_t=\theta(B)w_t$, $\phi(z)=1-\phi_1z-\dots-\phi_pz^p$, $\theta(z)=1+\theta_1z+\dots+\theta_qz^q$, $w_t$ white noise with variance $\sigma^2$. **State space:** state $\alpha_t\in\mathbb R^k$, observation $y_t\in\mathbb R^p$, $\alpha_{t+1}=T\alpha_t+R\eta_t$, $\eta_t\sim\mathcal N(0,Q)$; $y_t=Z\alpha_t+\varepsilon_t$, $\varepsilon_t\sim\mathcal N(0,H)$. $Y_t=(y_1,\dots,y_t)$. Filtered moments $a_{t|t},P_{t|t}$ (given $Y_t$), predicted $a_{t|t-1},P_{t|t-1}$ (given $Y_{t-1}$), smoothed $a_{t|n},P_{t|n}$. Innovation $v_t=y_t-Za_{t|t-1}$ with variance $F_t$; gain $K_t=P_{t|t-1}Z^\top F_t^{-1}$, so $a_{t|t}=a_{t|t-1}+K_tv_t$. **Volatility:** returns $\varepsilon_t=\sigma_tz_t$, $z_t$ iid mean 0 variance 1; GARCH(1,1) is $\sigma_t^2=\omega+\alpha\varepsilon_{t-1}^2+\beta\sigma_{t-1}^2$. In §3.15–3.17 only, $\omega,\alpha,\beta$ are GARCH parameters, not a frequency and a state. $\mathcal N(m,V)$ takes the covariance.

## 0. Where this sits

ST-7 described a stationary series by its autocovariance $\gamma(h)$ and fitted ARMA models to it. This module gives you two more views of the same object. The **frequency domain** says how the variance of the series is spread over frequencies, which is what an engineer reading a vibration record wants. **State space** says the series is a noisy readout of a hidden dynamical state, which handles missing data, unobserved components and nonlinear physics, and it is how ST-7's exact ARMA likelihood is computed in practice. A short third part models the variance itself (GARCH). You need ST-7's ARMA toolkit, F3's DFT and convolution theorem, and F5's Gaussian conditioning formula; you will derive that formula's consequences for the third time here. ML-4 uses §3.14's result that a Gaussian process with a Matérn kernel is a linear state-space model.

## 1. The problem that forced this idea

Arthur Schuster had a question about the weather and the Sun: is there a 26-day period (the solar rotation) in magnetic and meteorological records? In 1898 he proposed computing $|\sum_tx_te^{-2\pi i\omega t}|^2$ over a grid of trial frequencies and looking for a peak. He called it the periodogram, and he used it on earthquake records and, later, sunspots. It found peaks everywhere. Some were real; many were noise, and nobody could yet say which. G. Udny Yule's 1927 sunspot paper, the origin of ST-7's AR(2), was partly a reaction: a damped oscillator kicked by random shocks produces "cycles" without any strict period, and a periodogram of such a series is ragged. Eugen Slutsky showed the same year that moving-averaging pure noise creates convincing business cycles.

The theory came next. Norbert Wiener (1930) and Aleksandr Khinchin (1934) proved that for a stationary process, the autocovariance and a spectral measure are a Fourier pair. ST-7's time-domain theory and Schuster's frequency domain are the same theory written in two bases. Making it practical took John Tukey: with Ralph Blackman (1958) he explained why the raw periodogram must be smoothed, and with James Cooley (1965) he published the FFT (F3) that made spectral estimation cheap.

The second problem was navigation. In 1960 Rudolf Kalman published a recursive solution to "estimate the current state of a linear system from all noisy measurements so far, and update the estimate when a new one arrives, without reprocessing the past." Stanley Schmidt's group at NASA Ames adapted it, in its linearized (extended) form, for the Apollo midcourse navigation computer. The same recursion now runs in every GPS receiver, phone IMU and autopilot. For a statistician, the punchline is that ST-7's ARMA models are a special case: an ARMA series is a linear readout of a finite-dimensional state, and the Kalman filter gives its exact likelihood with missing values, time-varying parameters and unobserved components for free.

The third problem was finance. ST-7 modelled the conditional mean and assumed constant variance. Returns on assets have a nearly unpredictable mean and a strongly predictable variance: quiet weeks follow quiet weeks and turbulent ones follow turbulent ones. Robert Engle (1982) introduced ARCH (he first applied it to UK inflation uncertainty, not stock returns); Tim Bollerslev (1986) generalized it to GARCH. Engle shared the 2003 Nobel prize in economics for it.

Three forcing facts, then. (i) Oscillatory data (vibrations, tides, EEG, business cycles) is described naturally by frequency content, and the ACF hides it. (ii) "Estimate a hidden state from noisy observations" is the central problem of navigation, control and sensor fusion, and ARMA cannot even state it. (iii) Variance can be a time series too.

## 2. Intuition

**A spectrum is a variance budget.** Pass a stationary series through a bank of narrow band-pass filters, one per frequency, and measure the variance of each output. Those variances, per unit bandwidth, are $f(\omega)$; they add up to the total variance $\gamma(0)$. White noise has every lag uncorrelated, so no frequency is special: the spectrum is flat. An AR(1) with $\phi=0.7$ moves slowly (neighbouring values are similar), so its variance lives at low frequencies: a low-pass spectrum. Flip the sign to $\phi=-0.7$ and the series zigzags: the same spectrum mirrored to high frequencies. A damped oscillator kicked by noise (ST-7's AR(2) with complex roots) puts a bump near its ringing frequency. A pure sinusoid puts all its variance at one frequency: a spike.

```
 f(w)                       f(w)                     f(w)
  |\                         |     _                  |
  | \_  AR(1), phi>0         |    / \  AR(2),         |           white noise
  |   \___                   |   /   \__  complex     |-----------------------
  |       ''''-----          |__/       ''---         |
  +-------------------> w    +-------------------> w  +-------------------> w
  0                  1/2     0                  1/2   0                  1/2
```

**Filtering multiplies spectra.** If $y$ is $x$ passed through a linear filter, then at each frequency the filter multiplies the amplitude by $|\psi|$ and the variance by $|\psi|^2$. Guess what an ARMA spectrum is: white noise has flat spectrum $\sigma^2$, and ARMA is white noise through a rational filter, so $f$ must be $\sigma^2|\theta|^2/|\phi|^2$ evaluated on the unit circle. You will derive this in §3.3.

**The periodogram does not settle down.** The periodogram at frequency $\omega_j$ is the squared length of one DFT coefficient. With more data you get *more* Fourier frequencies, not better estimates at each one: each ordinate is still the squared modulus of one roughly Gaussian complex number, an exponential random variable with standard deviation equal to its mean. To estimate $f$ you must average neighbouring ordinates, trading resolution for stability. Guess the variance of one ordinate of white noise with variance 1: it is 1, for every $n$.

**The Kalman filter is predict-then-correct.** You track a cart on a rail. Your belief about its position is a Gaussian. Between measurements you push the belief forward with the dynamics (position += velocity × Δt) and it spreads because the dynamics are noisy. A measurement arrives, itself a Gaussian around the true position. The new belief is the product of the two Gaussians: a precision-weighted average whose weight on the measurement is

$$\text{gain}=\frac{\text{prior variance}}{\text{prior variance}+\text{measurement variance}}.$$

A sharp sensor and a vague prior give gain near 1 (trust the sensor); a noisy sensor and a sharp prior give gain near 0. That one ratio is the whole filter. You have seen it twice: as F5's conditional mean of one Gaussian given another, and as ST-5's Normal–Normal posterior mean. The Kalman filter applies it once per time step. Guess the steady-state gain when the process noise and measurement noise have equal variance: it is the golden ratio's reciprocal, 0.618 (§3.10).

**GARCH is ARMA for the squares.** Big returns follow big returns. So regress today's squared return on yesterday's: $\varepsilon_t^2\approx\omega+\alpha\varepsilon_{t-1}^2+\text{noise}$. That is an AR(1) in the squares, and adding $\beta\sigma^2_{t-1}$ makes it an ARMA(1,1). Everything ST-7 knows about ARMA (stationarity, mean reversion, forecasts) transfers to the variance. Guess the long-run variance of $\sigma_t^2=\omega+\alpha\varepsilon_{t-1}^2+\beta\sigma_{t-1}^2$: take expectations, $\bar\sigma^2=\omega+(\alpha+\beta)\bar\sigma^2$.

## 3. Formal development

### 3.1 Periodic components with random phase

Let $x_t=A\cos(2\pi\omega_0t+\Phi)$ with $A>0$ and $\omega_0\in(0,\tfrac12)$ fixed and $\Phi\sim\mathrm{Uniform}(0,2\pi)$.

▶ **Try it first.** Compute $\mathbb E x_t$ and $\mathbb E x_{t+h}x_t$ using $\cos a\cos b=\tfrac12[\cos(a-b)+\cos(a+b)]$.

1. $\mathbb E x_t=\frac{A}{2\pi}\int_0^{2\pi}\cos(2\pi\omega_0t+\phi)\,d\phi=0$: a full period of cosine.
2. $x_{t+h}x_t=\tfrac{A^2}{2}[\cos(2\pi\omega_0h)+\cos(2\pi\omega_0(2t+h)+2\Phi)]$. The second term integrates to zero over $\Phi$ (two full periods).
3. So $\gamma(h)=\tfrac{A^2}{2}\cos(2\pi\omega_0h)$, which depends on $h$ only: the process is weakly stationary. Writing the cosine as exponentials,
$$\gamma(h)=\tfrac{A^2}{4}e^{2\pi i\omega_0h}+\tfrac{A^2}{4}e^{-2\pi i\omega_0h}.\tag{3.1}$$
4. A sum $x_t=\sum_kA_k\cos(2\pi\omega_kt+\Phi_k)$ with independent phases has $\gamma(h)=\sum_k\tfrac{A_k^2}{2}\cos(2\pi\omega_kh)$ because cross terms have mean zero.

■ **What just happened.** A random-phase sinusoid is stationary, its ACF is a cosine that never decays, and its variance $A^2/2$ sits as two point masses of $A^2/4$ at $\pm\omega_0$. Every stationary process is a limit of such sums: **Herglotz's theorem** says $\gamma$ is nonnegative definite (ST-7) if and only if $\gamma(h)=\int_{(-1/2,1/2]}e^{2\pi i\omega h}\,dF(\omega)$ for a unique nondecreasing **spectral distribution** $F$ with total mass $\gamma(0)$. In symbols,
$$\gamma(h)=\int_{(-1/2,1/2]}e^{2\pi i\omega h}\,dF(\omega),\tag{3.2}$$
and (3.1) is the case where $F$ has two jumps. (Proof: Brockwell & Davis ch. 4; we prove the case with a density below.)

### 3.2 The spectral density and the Wiener–Khinchin theorem

**Definition 3.1.** If $\sum_h|\gamma(h)|<\infty$, the spectral density is
$$f(\omega)=\sum_{h=-\infty}^{\infty}\gamma(h)e^{-2\pi i\omega h},\qquad\omega\in[-\tfrac12,\tfrac12].\tag{3.3}$$

Note what kind of object this is. $h$ is an integer, so $f$ is a **Fourier series**, periodic in $\omega$ with period 1, and $\gamma(h)$ are its Fourier coefficients. Frequency lives on a circle. This is F2's setting (Fourier series of a periodic function), not F3's Fourier transform on $\mathbb R$, and the inversion formula is F2's coefficient formula.

**Theorem 3.2 (Wiener–Khinchin, summable case).** If $\sum|\gamma(h)|<\infty$ then $f$ in (3.3) is continuous, real, even, nonnegative, and
$$\gamma(h)=\int_{-1/2}^{1/2}f(\omega)e^{2\pi i\omega h}\,d\omega.\tag{3.4}$$

▶ **Try it first.** Multiply (3.3) by $e^{2\pi i\omega k}$ and integrate over $[-\tfrac12,\tfrac12]$. What is $\int_{-1/2}^{1/2}e^{2\pi i\omega m}d\omega$ for an integer $m$?

*Proof.*
1. $|\gamma(h)e^{-2\pi i\omega h}|=|\gamma(h)|$, and $\sum|\gamma(h)|<\infty$, so by the Weierstrass M-test the series converges uniformly in $\omega$ and $f$ is continuous.
2. Orthogonality on the circle: $\int_{-1/2}^{1/2}e^{2\pi i\omega m}\,d\omega=1$ if $m=0$ and $\frac{\sin(\pi m)}{\pi m}=0$ otherwise.
3. Uniform convergence lets us integrate term by term: $\int f(\omega)e^{2\pi i\omega k}d\omega=\sum_h\gamma(h)\int e^{2\pi i\omega(k-h)}d\omega=\gamma(k)$. That is (3.4).
4. Real and even: $\gamma(-h)=\gamma(h)$ is real, so pairing $h$ and $-h$ gives $f(\omega)=\gamma(0)+2\sum_{h\ge1}\gamma(h)\cos(2\pi\omega h)$, which is real and even in $\omega$.
5. Nonnegative: for any $n$, $0\le\frac1n\mathbb E\big|\sum_{t=1}^nx_te^{-2\pi i\omega t}\big|^2=\frac1n\sum_{s,t}\gamma(t-s)e^{-2\pi i\omega(t-s)}$. The number of pairs $(s,t)$ in $\{1..n\}^2$ with $t-s=h$ is $n-|h|$, so this equals
$$\sum_{|h|<n}\Big(1-\frac{|h|}{n}\Big)\gamma(h)e^{-2\pi i\omega h}.\tag{3.5}$$
As $n\to\infty$ each term tends to $\gamma(h)e^{-2\pi i\omega h}$ and the terms are dominated by $|\gamma(h)|$, so (3.5) $\to f(\omega)$ by dominated convergence. A limit of nonnegative numbers is nonnegative. ∎

Setting $h=0$ in (3.4):
$$\gamma(0)=\int_{-1/2}^{1/2}f(\omega)\,d\omega=2\int_0^{1/2}f(\omega)\,d\omega.\tag{3.6}$$

■ **What just happened.** The ACF and the spectral density carry the same information; (3.3) and (3.4) translate between them. (3.6) is the variance budget: $2f(\omega)\,d\omega$ is the variance contributed by frequencies in $[\omega,\omega+d\omega]$, counting $\pm\omega$ together. Step 5 also showed that (3.5) is the expected periodogram (§3.5), which is why the periodogram estimates $f$ at all.

### 3.3 The spectra of ARMA models

**White noise.** $\gamma(0)=\sigma^2$ and $\gamma(h)=0$ otherwise, so (3.3) gives
$$f_w(\omega)=\sigma^2\quad\text{for all }\omega.\tag{3.7}$$
Flat: every frequency carries the same variance, which is where the name "white" comes from.

**Theorem 3.3 (Linear filter theorem).** If $x_t$ is stationary with spectral density $f_x$ and $y_t=\sum_j\psi_jx_{t-j}$ with $\sum|\psi_j|<\infty$, then
$$f_y(\omega)=|\psi(e^{-2\pi i\omega})|^2f_x(\omega).\tag{3.8}$$

▶ **Try it first.** Write $\gamma_y(h)=\mathrm{Cov}(y_{t+h},y_t)$ as a double sum over $j,k$, then substitute into (3.3) and change the summation variable to $m=h-j+k$.

*Proof.*
1. $\gamma_y(h)=\sum_{j,k}\psi_j\psi_k\,\mathrm{Cov}(x_{t+h-j},x_{t-k})=\sum_{j,k}\psi_j\psi_k\gamma_x(h-j+k)$.
2. $f_y(\omega)=\sum_h\sum_{j,k}\psi_j\psi_k\gamma_x(h-j+k)e^{-2\pi i\omega h}$. Write $e^{-2\pi i\omega h}=e^{-2\pi i\omega j}\,e^{+2\pi i\omega k}\,e^{-2\pi i\omega(h-j+k)}$.
3. For fixed $j,k$, sum over $h$ with $m=h-j+k$: $\sum_m\gamma_x(m)e^{-2\pi i\omega m}=f_x(\omega)$. (Absolute summability of $\psi$ and $\gamma_x$ justifies reordering.)
4. What remains is $\big(\sum_j\psi_je^{-2\pi i\omega j}\big)\big(\sum_k\psi_ke^{2\pi i\omega k}\big)f_x=\psi(e^{-2\pi i\omega})\overline{\psi(e^{-2\pi i\omega})}f_x$, using real $\psi_k$. ∎

■ **What just happened.** This is F3's convolution theorem in second-moment form: convolution in time ($y=\psi*x$) is multiplication in frequency, and variance picks up the squared modulus. Steps 1–3 are the whole proof of "filtering multiplies spectra."

**ARMA.** $\phi(B)x_t=\theta(B)w_t$ says that the filter $\phi$ applied to $x$ equals the filter $\theta$ applied to $w$. Apply (3.8) to both sides and use (3.7):
$$f_x(\omega)=\sigma^2\frac{|\theta(e^{-2\pi i\omega})|^2}{|\phi(e^{-2\pi i\omega})|^2}.\tag{3.9}$$
A ratio of trigonometric polynomials: a **rational spectrum**. Roots of $\phi$ near the unit circle are poles that make peaks; roots of $\theta$ near the circle make troughs.

**AR(1).** With $\lambda=2\pi\omega$, $|1-\phi e^{-i\lambda}|^2=(1-\phi\cos\lambda)^2+\phi^2\sin^2\lambda=1-2\phi\cos\lambda+\phi^2$, so
$$f(\omega)=\frac{\sigma^2}{1-2\phi\cos(2\pi\omega)+\phi^2}.\tag{3.10}$$
The denominator is increasing in $\omega\in[0,\tfrac12]$ when $\phi>0$ (because $-2\phi\cos\lambda$ increases), so the spectrum falls from $\sigma^2/(1-\phi)^2$ at $\omega=0$ to $\sigma^2/(1+\phi)^2$ at $\omega=\tfrac12$: low-pass. For $\phi<0$ it rises: high-pass. Replacing $\phi$ by $-\phi$ maps $\cos\lambda\to-\cos\lambda$, i.e. $\omega\to\tfrac12-\omega$: the two spectra are mirror images.

**MA(1).** $f(\omega)=\sigma^2|1+\theta e^{-i\lambda}|^2=\sigma^2(1+2\theta\cos\lambda+\theta^2)$. Low-pass for $\theta>0$, high-pass for $\theta<0$, and never zero unless $|\theta|=1$ (the non-invertible boundary of ST-7).

**AR(2) and the location of the peak.** Let $c=\cos\lambda$.

▶ **Try it first.** Expand $D(\lambda)=|1-\phi_1e^{-i\lambda}-\phi_2e^{-2i\lambda}|^2$ as a polynomial in $c$ and minimize it.

1. Multiply by the conjugate: $D=1+\phi_1^2+\phi_2^2-2\phi_1\cos\lambda-2\phi_2\cos2\lambda+2\phi_1\phi_2\cos\lambda$. (The cross term $\phi_1\phi_2(e^{-i\lambda}e^{2i\lambda}+e^{i\lambda}e^{-2i\lambda})=2\phi_1\phi_2\cos\lambda$.)
2. Use $\cos2\lambda=2c^2-1$: $D(c)=1+\phi_1^2+\phi_2^2+2\phi_2+2\phi_1(\phi_2-1)c-4\phi_2c^2$.
3. $dD/dc=2\phi_1(\phi_2-1)-8\phi_2c=0$ gives
$$\cos(2\pi\omega^*)=\frac{\phi_1(\phi_2-1)}{4\phi_2}.\tag{3.11}$$
4. $d^2D/dc^2=-8\phi_2$. This is positive, so $\omega^*$ minimizes $D$ and maximizes $f$, exactly when $\phi_2<0$; and it lies in $(0,\tfrac12)$ only if the right side of (3.11) is in $(-1,1)$.
5. Compare with the roots. ST-7 showed that complex roots of $\phi$ have reciprocals $re^{\pm i\theta}$ with $\phi_1=2r\cos\theta$, $\phi_2=-r^2$. Substitute into (3.11): $\cos(2\pi\omega^*)=\cos\theta\cdot\frac{1+r^2}{2r}$. Since $\frac{1+r^2}{2r}\ge1$ (AM–GM), $\cos2\pi\omega^*\ge\cos\theta$, so $\omega^*\le\theta/2\pi$, with equality only as $r\to1$.

■ **What just happened.** Yule's quasi-period $2\pi/\theta$ from the roots and the spectral peak period $1/\omega^*$ are not the same number. Damping ($r<1$) broadens the resonance and pushes its peak toward lower frequency, the discrete-time version of CM-1's result that a damped oscillator's amplitude resonance sits below its natural frequency. As the roots approach the unit circle the two agree and the peak sharpens.

### 3.4 Linear filters

Every moving average and difference is a filter, and (3.8) says what it does to each frequency.

**First difference** $\nabla x_t=x_t-x_{t-1}$: $\psi(z)=1-z$ and
$$|1-e^{-2\pi i\omega}|^2=2-2\cos2\pi\omega=4\sin^2(\pi\omega).$$
Zero at $\omega=0$, 4 at $\omega=\tfrac12$: a high-pass filter. Differencing a trending series (ST-7) kills its lowest frequencies, which is why over-differencing leaves an MA spectrum with a hole at zero.

**The $k$-point moving average** $\psi_j=1/k$ for $j=0,\dots,k-1$: a geometric sum,
$$\psi(e^{-i\lambda})=\frac1k\frac{1-e^{-ik\lambda}}{1-e^{-i\lambda}}=\frac1ke^{-i(k-1)\lambda/2}\frac{\sin(k\lambda/2)}{\sin(\lambda/2)},\qquad|\psi|^2=\Big(\frac{\sin(\pi k\omega)}{k\sin(\pi\omega)}\Big)^2.\tag{3.12}$$
This is F3's Dirichlet kernel (the discrete sinc): 1 at $\omega=0$, zeros at $\omega=1/k,2/k,\dots$, and side lobes in between. A moving average is a low-pass filter with a leaky stopband: frequencies between its zeros pass through, attenuated but not removed. The phase factor $e^{-i(k-1)\lambda/2}$ is a delay of $(k-1)/2$ samples; centring the average removes it.

**Filter design intuition.** Sharp cutoffs need long filters (a short filter has a smooth frequency response: F3's time–frequency uncertainty). Cascading filters multiplies responses. Every practical smoother is a compromise between cutoff sharpness, side-lobe height and length.

**The Slutsky effect.** Apply $(1+B)^s$ and then $(1-B)^d$ to white noise. By (3.8) the output spectrum is $\sigma^2(4\cos^2\pi\omega)^s(4\sin^2\pi\omega)^d$. Take logs and differentiate in $\omega$: $-2s\pi\tan\pi\omega+2d\pi\cot\pi\omega=0$, so $\tan^2\pi\omega^*=d/s$ and $\cos2\pi\omega^*=\frac{1-\tan^2}{1+\tan^2}=\frac{s-d}{s+d}$. The product of a low-pass and a high-pass is a band-pass, and its output looks like a cycle of period $1/\omega^*$ that no mechanism produced. Slutsky's 1927 point, made quantitative: smoothing and differencing economic data manufactures "cycles." For $s=6$, $d=2$ the cycle has period exactly 6.

### 3.5 The periodogram

**Definition 3.4.** For data $x_1,\dots,x_n$ and $\omega_j=j/n$, the DFT is $d(\omega_j)=n^{-1/2}\sum_{t=1}^nx_te^{-2\pi i\omega_jt}$ and the periodogram is
$$I(\omega_j)=|d(\omega_j)|^2=\frac1n\Big|\sum_{t=1}^nx_te^{-2\pi i\omega_jt}\Big|^2.\tag{3.13}$$

For $j\neq0$, $\sum_te^{-2\pi i\omega_jt}=0$, so subtracting the mean changes nothing except $I(0)$; we always demean. By (3.5), $\mathbb EI(\omega)\to f(\omega)$: the periodogram is asymptotically unbiased.

**Theorem 3.5 (Distribution of the periodogram, Gaussian white noise).** If $x_t$ are iid $\mathcal N(0,\sigma^2)$ and $0<\omega_j<\tfrac12$, then $I(\omega_j)\sim\tfrac{\sigma^2}{2}\chi^2_2$, so $\mathbb EI=\sigma^2$ and $\mathrm{Var}\,I=\sigma^4$; ordinates at distinct such $\omega_j$ are independent.

▶ **Try it first.** Split $d(\omega_j)=C_j-iS_j$ into cosine and sine sums. What are their variances and covariance? You need F3's discrete orthogonality: $\sum_t\cos^2(2\pi\omega_jt)=\sum_t\sin^2(2\pi\omega_jt)=n/2$ and $\sum_t\cos(2\pi\omega_jt)\sin(2\pi\omega_kt)=0$ for $0<\omega_j,\omega_k<\tfrac12$.

*Proof.*
1. $C_j=n^{-1/2}\sum_tx_t\cos(2\pi\omega_jt)$ and $S_j=n^{-1/2}\sum_tx_t\sin(2\pi\omega_jt)$ are linear in Gaussian $x$, so all $C_j,S_j$ are jointly Gaussian with mean 0.
2. $\mathrm{Var}\,C_j=\frac{\sigma^2}{n}\sum_t\cos^2(2\pi\omega_jt)=\frac{\sigma^2}{2}$; likewise $\mathrm{Var}\,S_j=\frac{\sigma^2}2$.
3. All covariances among $\{C_j,S_j\}$ for distinct cosines/sines are $\frac{\sigma^2}{n}\times$ (an orthogonality sum) $=0$. Jointly Gaussian and uncorrelated means independent.
4. $I(\omega_j)=C_j^2+S_j^2=\frac{\sigma^2}{2}(U^2+V^2)$ with $U,V$ iid $\mathcal N(0,1)$, i.e. $\frac{\sigma^2}2\chi^2_2$. A $\chi^2_2$ has mean 2 and variance 4, so $\mathbb EI=\sigma^2$, $\mathrm{Var}\,I=\frac{\sigma^4}{4}\cdot4=\sigma^4$. ∎

For a general linear process with $\sum|\psi_j|\,|j|^{1/2}<\infty$ the same holds asymptotically with $\sigma^2$ replaced by $f(\omega_j)$:
$$\frac{2I(\omega_j)}{f(\omega_j)}\ \xrightarrow{d}\ \chi^2_2,\qquad\text{independently across fixed }\omega_j\in(0,\tfrac12).\tag{3.14}$$
(The proof writes $d_x(\omega)\approx\psi(e^{-2\pi i\omega})d_w(\omega)$ and applies Theorem 3.5 to $d_w$; Brockwell & Davis, *Time Series: Theory and Methods*, ch. 10, has the details.)

■ **What just happened, and the key fact.** $\mathrm{Var}\,I(\omega_j)\approx f(\omega_j)^2$ **does not shrink as $n$ grows.** The periodogram is not a consistent estimator of $f$. Its standard deviation equals its mean at every sample size; a raw periodogram's 95% interval runs from $0.27I$ to $39.5I$ (two-sided $\chi^2_2$ quantiles). Longer data buys more ordinates, each equally noisy. Schuster's spurious peaks were this.

**The fix: smooth across frequencies.** The ordinates at neighbouring Fourier frequencies are approximately independent and, if $f$ is smooth, have nearly the same mean. Average $L=2m+1$ of them (the **Daniell** estimator):
$$\hat f(\omega_j)=\frac1L\sum_{k=-m}^{m}I(\omega_{j+k}).\tag{3.15}$$

▶ **Try it first.** If $f$ is constant across the band, what is the distribution of $\hat f/f$? What do the mean and variance of $\hat f$ become?

1. By (3.14), $2I(\omega_{j+k})/f\approx\chi^2_2$ independently, so $2L\hat f/f\approx\chi^2_{2L}$: a sum of $L$ independent $\chi^2_2$.
2. $\mathbb E\hat f=f$ and $\mathrm{Var}\,\hat f=f^2\cdot\frac{2\cdot2L}{(2L)^2}=f^2/L$. The variance now falls as the span grows.
3. Define the **equivalent degrees of freedom** $\nu=2L$. For general weights $h_k\ge0$ summing to 1, the same moment matching ($\mathbb E\hat f=f$, $\mathrm{Var}\,\hat f=f^2\sum h_k^2$, matched to a scaled $\chi^2_\nu$ with variance $2f^2/\nu$) gives $\nu=2/\sum_kh_k^2$. The **modified Daniell** kernel gives the two end weights half weight, which smooths the jump at the band edge.
4. Since $\nu\hat f/f\approx\chi^2_\nu$, $P(\chi^2_{\nu,0.025}\le\nu\hat f/f\le\chi^2_{\nu,0.975})=0.95$. Solve for $f$:
$$\Big[\frac{\nu\hat f(\omega)}{\chi^2_{\nu,0.975}},\ \frac{\nu\hat f(\omega)}{\chi^2_{\nu,0.025}}\Big].\tag{3.16}$$
The interval is multiplicative, so on a log-scale plot it is a bar of constant length: plot spectra on log axes.
5. **Bias.** If $f$ curves, the average in (3.15) is biased. Taylor-expand $f(\omega_j+k/n)\approx f+f'k/n+\tfrac12f''k^2/n^2$; the linear terms cancel symmetrically and $\sum_{k=-m}^mk^2=m(m+1)(2m+1)/3$, so $\mathbb E\hat f-f\approx f''(\omega_j)\frac{m(m+1)}{6n^2}$. With bandwidth $B=L/n$, bias $\propto B^2$ and variance $\propto1/L=1/(nB)$.
6. $\mathrm{MSE}\approx c_1B^4+c_2/(nB)$, minimized where $4c_1B^3=c_2/(nB^2)$, i.e. $B^*\propto n^{-1/5}$.

■ **What just happened.** Smoothing turns an inconsistent estimator into a consistent one if the span grows with $n$ but more slowly than $n$ ($L\to\infty$, $L/n\to0$). Wide spans flatten peaks (bias), narrow spans leave noise (variance); there is no free choice, only a trade-off, and sharp peaks need narrow spans. The project's bias–variance plot shows it at a sharp AR(2) peak.

### 3.6 Tapering, parametric estimates, and hidden periodicities

**Leakage and tapers.** Observing $x_t$ only for $t=1..n$ multiplies the infinite series by a rectangular window. By F3's convolution theorem the DFT of the product is the DFT of the series convolved with the window's DFT, the Dirichlet kernel (3.12) with its slowly decaying side lobes (power falling like $1/(n\Delta)^2$ at distance $\Delta$ from a peak). A strong sinusoid's power **leaks** across the whole spectrum and can bury a weak component far away. A **taper** $h_t$ that goes smoothly to zero at both ends, such as the Hann taper $h_t=\sin^2(\pi t/n)$, has a kernel with side lobes falling like $1/\Delta^6$ in power. Compute $|\sum h_tx_te^{-2\pi i\omega t}|^2/\sum h_t^2$ (the normalization keeps $\mathbb EI\approx f$). The cost is a wider main lobe and fewer effective degrees of freedom, because the ends of the record are down-weighted. The standard recipe is **taper, then smooth**: tapering controls leakage (bias far from peaks), smoothing controls variance. In the project, a sinusoid of amplitude 0.01 sitting 0.2 cycles away from one of amplitude 10 is invisible without a taper and stands three orders of magnitude above the background with one.

**Parametric spectral estimation.** Fit an AR($p$) by Yule–Walker (ST-7) with $p$ chosen by AIC, and plug the estimates into (3.9). This is smooth by construction and resolves sharp peaks well, because a fitted pole near the unit circle produces a narrow peak with very few parameters. It is also the reason for caution: an over-fitted AR can create sharp peaks that are not there, and its intervals are model-based.

**Fisher's test for a hidden periodicity.** Under Gaussian white noise, Theorem 3.5 says the ordinates $I_1,\dots,I_m$ ($m=\lfloor(n-1)/2\rfloor$) are iid exponential. Fisher (1929) used $g=\max_jI_j/\sum_jI_j$ and derived its exact null distribution, $P(g>x)=\sum_{k=1}^{\lfloor1/x\rfloor}(-1)^{k-1}\binom mk(1-kx)^{m-1}$. Dividing by the sum makes the test free of the unknown $\sigma^2$. The test assumes the background is white; for coloured backgrounds compare each ordinate with a smoothed $\hat f$ instead. The project finds a sinusoid of amplitude 0.4 in unit-variance noise with $p\approx10^{-7}$ and checks that the test rejects 5% of the time under the null.

### 3.7 Cross-spectrum, coherence and transfer functions

For two jointly stationary series, the cross-covariance $\gamma_{yx}(h)=\mathrm{Cov}(y_{t+h},x_t)$ has Fourier series $f_{yx}(\omega)=\sum_h\gamma_{yx}(h)e^{-2\pi i\omega h}$, the **cross-spectrum**. It is complex because $\gamma_{yx}$ need not be even. Its argument is the **phase** (a lead or lag at each frequency) and the **squared coherence** is
$$\rho^2_{yx}(\omega)=\frac{|f_{yx}(\omega)|^2}{f_{xx}(\omega)f_{yy}(\omega)}\in[0,1].\tag{3.17}$$

▶ **Try it first.** Suppose $y_t=\sum_j\psi_jx_{t-j}+v_t$ with $v$ stationary and independent of $x$. Compute $\gamma_{yx}(h)$ and then $f_{yx}$.

1. $\gamma_{yx}(h)=\sum_j\psi_j\mathrm{Cov}(x_{t+h-j},x_t)=\sum_j\psi_j\gamma_{xx}(h-j)$.
2. The same reindexing as Theorem 3.3 gives $f_{yx}(\omega)=\psi(e^{-2\pi i\omega})f_{xx}(\omega)$, so the frequency response is
$$\psi(e^{-2\pi i\omega})=\frac{f_{yx}(\omega)}{f_{xx}(\omega)}.\tag{3.18}$$
3. By (3.8) and independence, $f_{yy}=|\psi|^2f_{xx}+f_{vv}$. Put both into (3.17): $\rho^2_{yx}=\frac{|\psi|^2f_{xx}}{|\psi|^2f_{xx}+f_{vv}}$, the fraction of $y$'s variance at $\omega$ that is linearly explained by $x$.
4. Coherence is a squared correlation at each frequency: $\mathbb E[d_y\bar d_x]\to f_{yx}$, $\mathbb E|d_x|^2\to f_{xx}$, $\mathbb E|d_y|^2\to f_{yy}$ (the cross version of (3.5)), so (3.17) is $|\mathrm{corr}(d_y(\omega),d_x(\omega))|^2$ for the DFT coefficients. By Cauchy–Schwarz it lies in $[0,1]$.

■ **What just happened.** Regression, done one frequency at a time. Estimate $f_{yx}$, $f_{xx}$, $f_{yy}$ by smoothing $d_y\bar d_x$, $|d_x|^2$, $|d_y|^2$ as in (3.15) (smoothing is essential: unsmoothed, $|d_y\bar d_x|^2/(|d_x|^2|d_y|^2)=1$ identically). Then (3.18) identifies an unknown linear system from its input and output, which is how mechanical engineers measure a structure's frequency response: shake it with broadband noise, record the response, divide spectra. Coherence tells you at which frequencies the estimate is trustworthy. A pure delay $y_t=x_{t-d}$ has $\psi=e^{-2\pi i\omega d}$: unit gain and phase $-2\pi\omega d$, linear in $\omega$. The project identifies a damped oscillator ($\zeta=0.05$, resonance at 0.078 cycles/sample) from $2^{15}$ samples of white-noise forcing with median gain error 1.2%.

### 3.8 State-space models

**Definition 3.6.** A linear Gaussian state-space model is
$$\alpha_{t+1}=T\alpha_t+R\eta_t,\quad\eta_t\sim\mathcal N(0,Q);\qquad y_t=Z\alpha_t+\varepsilon_t,\quad\varepsilon_t\sim\mathcal N(0,H),\tag{3.19}$$
with $\alpha_1\sim\mathcal N(a_{1|0},P_{1|0})$ and all of $\alpha_1,\{\eta_t\},\{\varepsilon_t\}$ independent. $T$ ($k\times k$) is the transition matrix, $Z$ ($p\times k$) the observation matrix, $R$ ($k\times r$) selects which states receive noise. Any of them may depend on $t$.

The state is Markov (F5): given $\alpha_t$, the future is independent of the past. Each observation depends only on the current state. The examples, which recur through the module:

- **AR($p$)** via the companion matrix (LA-7, ST-7): $\alpha_t=(x_t,x_{t-1},\dots,x_{t-p+1})^\top$, first row of $T$ is $(\phi_1,\dots,\phi_p)$, the subdiagonal is the identity, $R=e_1$, $Z=e_1^\top$, $Q=\sigma^2$, $H=0$. For AR(2): $T=\begin{pmatrix}\phi_1&\phi_2\\1&0\end{pmatrix}$ (Example 4.5).
- **ARMA($p,q$)** with $r=\max(p,q+1)$: $T$ has $(\phi_1,\dots,\phi_r)$ (zero-padded) in its first **column** and ones on the superdiagonal, $R=(1,\theta_1,\dots,\theta_{r-1})^\top$, $Z=e_1^\top$. For ARMA(1,1), $T=\begin{pmatrix}\phi&1\\0&0\end{pmatrix}$, $R=\binom1\theta$: then $\alpha_{t+1}^{(2)}=\theta w_{t+1}$ and $x_{t+1}=\phi x_t+\theta w_t+w_{t+1}$. Other books use other, equivalent forms (§8).
- **Local level:** $\mu_{t+1}=\mu_t+\eta_t$, $y_t=\mu_t+\varepsilon_t$. A random walk seen through noise. $T=Z=R=1$.
- **Local linear trend:** level $\mu$ and slope $\nu$, $\mu_{t+1}=\mu_t+\nu_t+\eta_t$, $\nu_{t+1}=\nu_t+\zeta_t$. $T=\begin{pmatrix}1&1\\0&1\end{pmatrix}$, $R=I$, $Q=\mathrm{diag}(\sigma_\eta^2,\sigma_\zeta^2)$.
- **Seasonal** of period $s$: $\gamma_{t+1}=-\sum_{j=0}^{s-2}\gamma_{t-j}+\omega_t$, so the $s$ seasonal effects sum to zero up to noise. State $(\gamma_t,\dots,\gamma_{t-s+2})$.
- **Constant velocity (the ME example).** Position and velocity, $\alpha=(p,u)^\top$, sampled every $\Delta t$, acceleration a white-noise shock $a_t$ held constant over each interval: $p_{t+1}=p_t+u_t\Delta t+\tfrac12a_t\Delta t^2$, $u_{t+1}=u_t+a_t\Delta t$. So $T=\begin{pmatrix}1&\Delta t\\0&1\end{pmatrix}$, $R=\binom{\Delta t^2/2}{\Delta t}$, $Q=\sigma_a^2$, and a position sensor has $Z=(1\ 0)$, $H=\sigma_y^2$.
- **Time-varying regression:** $y_t=x_t^\top\beta_t+\varepsilon_t$, $\beta_{t+1}=\beta_t+\eta_t$: $Z_t=x_t^\top$, $T=I$.

Why the form is general: any linear recursion of finite order can be stacked into a first-order vector recursion (the same trick that turns an $n$th-order ODE into a first-order system), and independent components combine by stacking states block-diagonally and adding their $Z$ rows. Section 3.13 builds trend + seasonal this way.

### 3.9 The Kalman filter derived

We need the conditional distribution of $\alpha_t$ given $Y_t$. Because everything in (3.19) is jointly Gaussian, every conditional distribution is Gaussian and is fixed by its mean and covariance. We carry $(a_{t|t},P_{t|t})$ forward in two steps.

You are about to derive F5's MVN conditioning formula for the third time: once in F5 itself, once in ST-5 as the Normal–Normal posterior, and now as a recursion. If it still feels like algebra, do the derivation by hand; this is usually the time it stops feeling like algebra.

**Step A: prediction.**

▶ **Try it first.** Given $\alpha_{t-1}\mid Y_{t-1}\sim\mathcal N(a_{t-1|t-1},P_{t-1|t-1})$ and $\alpha_t=T\alpha_{t-1}+R\eta_{t-1}$ with $\eta_{t-1}$ independent of $(\alpha_{t-1},Y_{t-1})$, find the distribution of $\alpha_t\mid Y_{t-1}$.

1. A linear map of a Gaussian plus independent Gaussian noise is Gaussian (F5).
2. Mean: $\mathbb E[T\alpha_{t-1}+R\eta_{t-1}\mid Y_{t-1}]=Ta_{t-1|t-1}$.
3. Covariance: the two terms are independent given $Y_{t-1}$, so variances add: $TP_{t-1|t-1}T^\top+RQR^\top$.
$$a_{t|t-1}=Ta_{t-1|t-1},\qquad P_{t|t-1}=TP_{t-1|t-1}T^\top+RQR^\top.\tag{3.20}$$

**Step B: update.**

▶ **Try it first.** Write the joint distribution of $(\alpha_t,y_t)$ given $Y_{t-1}$, then condition on $y_t$ using F5's formula: for $\binom{u}{v}\sim\mathcal N\big(\binom{m_u}{m_v},\begin{pmatrix}S_{uu}&S_{uv}\\S_{vu}&S_{vv}\end{pmatrix}\big)$, $u\mid v\sim\mathcal N(m_u+S_{uv}S_{vv}^{-1}(v-m_v),\ S_{uu}-S_{uv}S_{vv}^{-1}S_{vu})$.

4. $y_t=Z\alpha_t+\varepsilon_t$ with $\varepsilon_t$ independent of $(\alpha_t,Y_{t-1})$. So given $Y_{t-1}$: $\mathbb Ey_t=Za_{t|t-1}$; $\mathrm{Var}\,y_t=ZP_{t|t-1}Z^\top+H$; $\mathrm{Cov}(\alpha_t,y_t)=\mathrm{Cov}(\alpha_t,Z\alpha_t)=P_{t|t-1}Z^\top$.
5. Name the pieces:
$$v_t=y_t-Za_{t|t-1},\qquad F_t=ZP_{t|t-1}Z^\top+H.\tag{3.21}$$
$v_t$ is the **innovation**, the part of $y_t$ the past could not predict, and $F_t$ is its variance.
6. Apply F5's formula with $u=\alpha_t$, $v=y_t$: the conditional mean is $a_{t|t-1}+P_{t|t-1}Z^\top F_t^{-1}v_t$. Define the gain
$$K_t=P_{t|t-1}Z^\top F_t^{-1}.\tag{3.22}$$
7. The conditional covariance is $P_{t|t-1}-P_{t|t-1}Z^\top F_t^{-1}ZP_{t|t-1}$: the **Schur complement** of $F_t$ in the joint covariance (LA-10). Since $K_tF_tK_t^\top=P_{t|t-1}Z^\top F_t^{-1}F_tF_t^{-1}ZP_{t|t-1}$, this is $P_{t|t-1}-K_tF_tK_t^\top$.
$$a_{t|t}=a_{t|t-1}+K_tv_t,\qquad P_{t|t}=P_{t|t-1}-K_tF_tK_t^\top.\tag{3.23}$$

**The Kalman filter** is: start from $(a_{1|0},P_{1|0})$; for $t=1,\dots,n$ compute (3.21), (3.22), (3.23), then predict with (3.20) for $t+1$. The project's `kalman_filter(y, T, Z, R, Q, H, a0, P0)` is these lines. Each step costs $O(k^3+p^3)$, so the whole filter is linear in $n$.

Note that $P_{t|t}$, $P_{t|t-1}$ and $K_t$ never depend on the data. They can be computed in advance; only the means use $y$.

■ **What just happened.** Prediction is "push a Gaussian through linear dynamics"; update is "condition a Gaussian on a linear observation." Neither step needs anything except F5. The gain $K_t$ is a regression coefficient: the coefficient of the regression of $\alpha_t$ on $y_t$ given the past.

### 3.10 What the filter means: innovations, precision, stability, steady state

**Theorem 3.7 (Innovations).** The innovations $v_1,\dots,v_n$ are independent, $v_t\sim\mathcal N(0,F_t)$.

*Proof.*
1. $v_t=y_t-\mathbb E[y_t\mid Y_{t-1}]$, so $\mathbb E[v_t\mid Y_{t-1}]=0$, and given $Y_{t-1}$, $v_t\sim\mathcal N(0,F_t)$ by step 4 above. Since $F_t$ does not depend on $Y_{t-1}$, this is also the unconditional distribution.
2. For $s<t$, $v_s$ is a function of $Y_s\subseteq Y_{t-1}$. Tower property: $\mathbb E[v_tv_s^\top]=\mathbb E\big[\mathbb E[v_t\mid Y_{t-1}]v_s^\top\big]=0$.
3. The $v$'s are linear in the jointly Gaussian $y$'s, hence jointly Gaussian; uncorrelated jointly Gaussian vectors are independent. ∎

So a correctly specified filter turns the correlated $y$'s into white noise. That is the diagnostic: standardized innovations $v_t/\sqrt{F_t}$ should pass ST-7's Ljung–Box test and look $\mathcal N(0,1)$.

**The gain as relative precision.** For the local level model, $Z=1$ and (3.22) gives
$$K_t=\frac{P_{t|t-1}}{P_{t|t-1}+H}=\frac{1/H}{1/P_{t|t-1}+1/H},\qquad a_{t|t}=(1-K_t)a_{t|t-1}+K_ty_t.$$
The filtered mean is a precision-weighted average of prior mean and observation. This is ST-5's Normal–Normal posterior mean, and the weight is the shrinkage factor that ST-5 also found in the Beta–Binomial model. In matrix form the matrix inversion lemma rewrites (3.23) as $P_{t|t}^{-1}=P_{t|t-1}^{-1}+Z^\top H^{-1}Z$: precisions add. The filtered distribution **is** the Bayesian posterior of $\alpha_t$ given $Y_t$, and the predicted distribution is the next step's prior. The Kalman filter is recursive Bayesian updating, and with $T=I$, $Q=0$ it is recursive least squares for a fixed regression coefficient.

**Numerical stability.** In floating point, $P_{t|t-1}-K_tF_tK_t^\top$ subtracts two nearly equal matrices when the observation is precise, and rounding can leave $P_{t|t}$ asymmetric or with a negative eigenvalue (LA-13: the subtraction is ill-conditioned). The **Joseph form** $P_{t|t}=(I-K_tZ)P_{t|t-1}(I-K_tZ)^\top+K_tHK_t^\top$ is algebraically equal (expand it, using (3.22)) and is a sum of PSD terms, so it stays PSD under rounding. Production filters go further and propagate a Cholesky factor of $P$ (LA-10), the **square-root filter**, which also doubles the effective precision. The project symmetrizes $P$ after each update, the cheapest defence.

**The steady state and the Riccati equation.** Combining (3.20) and (3.23), the predicted covariance obeys the **discrete Riccati recursion** $P_{t+1|t}=T\big(P_{t|t-1}-P_{t|t-1}Z^\top F_t^{-1}ZP_{t|t-1}\big)T^\top+RQR^\top$. For a time-invariant, observable, stable-enough model it converges to a fixed point, and the gain converges to a constant: the steady-state filter.

▶ **Try it first.** For the local level model with $Q,H>0$, set $P_{t+1|t}=P_{t|t-1}=P$ and solve.

1. The recursion is $P_{t+1|t}=P_{t|t-1}-\frac{P_{t|t-1}^2}{P_{t|t-1}+H}+Q=\frac{P_{t|t-1}H}{P_{t|t-1}+H}+Q$.
2. A fixed point satisfies $\frac{P^2}{P+H}=Q$, i.e.
$$P^2-QP-QH=0.\tag{3.24}$$
3. The positive root is $P=\tfrac12\big(Q+\sqrt{Q^2+4QH}\big)$. With $q=Q/H$ (the signal-to-noise ratio) and $p=P/H$: $p=\tfrac12(q+\sqrt{q^2+4q})$ and $K=p/(p+1)$.
4. Convergence: $g(P)=PH/(P+H)+Q$ is increasing and concave with $g'(P)=H^2/(P+H)^2<1$ for $P>0$, so iterating $g$ from any $P\ge0$ converges monotonically to the unique positive fixed point.

■ **What just happened.** The steady-state gain depends only on $q$: $K\to0$ as $q\to0$ (a constant level; average everything), $K\to1$ as $q\to\infty$ (the level jumps; believe the latest point). With $K$ fixed, $a_{t|t}=a_{t-1|t-1}+K(y_t-a_{t-1|t-1})$ is ST-7's **simple exponential smoothing** with smoothing constant $K$. The ad hoc forecasting rule is the optimal filter for a local level model, and the model tells you how to choose the constant. For the constant-velocity model the steady-state filter is the tracking engineer's "$\alpha$–$\beta$ filter."

### 3.11 Smoothing, forecasting and missing data

The filter uses $Y_t$. After all the data are in, you want $\alpha_t\mid Y_n$, which also uses the future.

**Theorem 3.8 (Rauch–Tung–Striebel smoother).** With $J_t=P_{t|t}T^\top P_{t+1|t}^{-1}$, running backward from $a_{n|n},P_{n|n}$,
$$a_{t|n}=a_{t|t}+J_t(a_{t+1|n}-a_{t+1|t}),\qquad P_{t|n}=P_{t|t}+J_t(P_{t+1|n}-P_{t+1|t})J_t^\top.\tag{3.26}$$

▶ **Try it first.** Given $\alpha_{t+1}$, is $\alpha_t$ informed by $y_{t+1},\dots,y_n$? Use the Markov structure.

*Proof.*
1. Given $\alpha_{t+1}$, the future observations $y_{t+1},\dots,y_n$ depend on $\alpha_t$ only through $\alpha_{t+1}$ (they are functions of $\alpha_{t+1}$ and of noise after time $t$). So $p(\alpha_t\mid\alpha_{t+1},Y_n)=p(\alpha_t\mid\alpha_{t+1},Y_t)$.
2. Given $Y_t$, $(\alpha_t,\alpha_{t+1})$ is jointly Gaussian with means $(a_{t|t},Ta_{t|t})$, variances $P_{t|t}$ and $P_{t+1|t}$, and covariance $\mathrm{Cov}(\alpha_t,T\alpha_t+R\eta_t)=P_{t|t}T^\top$.
3. F5's formula: $\alpha_t\mid\alpha_{t+1},Y_t\sim\mathcal N\big(a_{t|t}+J_t(\alpha_{t+1}-a_{t+1|t}),\ P_{t|t}-J_tP_{t+1|t}J_t^\top\big)$.
4. Average over $\alpha_{t+1}\mid Y_n\sim\mathcal N(a_{t+1|n},P_{t+1|n})$. Law of total expectation: $a_{t|n}=a_{t|t}+J_t(a_{t+1|n}-a_{t+1|t})$. Law of total variance: $P_{t|n}=P_{t|t}-J_tP_{t+1|t}J_t^\top+J_tP_{t+1|n}J_t^\top$, which is (3.26). ∎

■ **What just happened.** One backward pass corrects each filtered estimate by how much the future revised its successor. Filtering is what a controller can use in real time; smoothing is what an analyst uses afterwards, and it is always at least as accurate ($P_{t|n}\preceq P_{t|t}$ because $P_{t+1|n}\preceq P_{t+1|t}$). In the project's tracker, filtering cuts position RMSE from 2.85 (raw sensor) to 1.76 and smoothing to 0.91.

**Forecasting.** With no new data, iterate (3.20): $a_{n+h|n}=T^ha_{n|n}$, $P_{n+h|n}=TP_{n+h-1|n}T^\top+RQR^\top$, and $y_{n+h}$ has forecast $Za_{n+h|n}$ and variance $ZP_{n+h|n}Z^\top+H$. For an AR($p$) in companion form this reproduces ST-7's forecasts and their growing variance.

**Missing data.** If $y_t$ is missing, there is nothing to condition on: skip (3.21)–(3.23), set $a_{t|t}=a_{t|t-1}$, $P_{t|t}=P_{t|t-1}$, and contribute nothing to the likelihood. The smoother then interpolates across the gap. Irregular sampling is handled the same way, with $T$ and $Q$ depending on the time step. ARMA in ST-7's form cannot do this without exactly this machinery.

### 3.12 The likelihood, and estimating the hyperparameters

**Theorem 3.9 (Prediction error decomposition).** For model (3.19) with parameters $\theta$ (entries of $T,Z,R,Q,H$),
$$\log p(y_1,\dots,y_n;\theta)=-\frac12\sum_{t=1}^n\Big[p\log2\pi+\log|F_t|+v_t^\top F_t^{-1}v_t\Big].\tag{3.27}$$

*Proof.*
1. Chain rule of densities: $p(y_1,\dots,y_n)=\prod_tp(y_t\mid Y_{t-1})$.
2. From step 4 of §3.9, $y_t\mid Y_{t-1}\sim\mathcal N(Za_{t|t-1},F_t)$, whose log density at $y_t$ is $-\tfrac12[p\log2\pi+\log|F_t|+v_t^\top F_t^{-1}v_t]$.
3. Sum the logs. ∎

(A second view: the map $(y_1,\dots,y_n)\mapsto(v_1,\dots,v_n)$ is lower block-triangular with identity diagonal blocks, so its Jacobian is 1, and Theorem 3.7 says the $v$'s are independent $\mathcal N(0,F_t)$.)

■ **What just happened.** The Kalman filter computes the exact Gaussian likelihood in $O(n)$, where the direct formula needs the $n\times n$ covariance and its Cholesky factor (LA-10) in $O(n^3)$. For ST-7's ARMA, put the model in the form of §3.8, start the filter at the stationary distribution $P_{1|0}=\Sigma$ solving the discrete Lyapunov equation $\Sigma=T\Sigma T^\top+RQR^\top$, and (3.27) is the exact likelihood. The project checks this against the direct multivariate Normal log density for an ARMA(1,1) with $n=200$: $-276.614139$ both ways. This is how statistical software computes ARMA likelihoods.

**Maximum likelihood.** Maximize (3.27) numerically over $\theta$, parametrizing variances by their logs so the optimizer is unconstrained, and get standard errors from the inverse of the numerical Hessian (ST-1's observed information; the delta method converts to the variance scale). For nonstationary states (random walks, trends), no stationary distribution exists; the usual approximation is a "diffuse" start with huge $P_{1|0}$, whose first few innovations contribute an almost parameter-free constant (Durbin & Koopman ch. 5 treat the exact diffuse likelihood). The project fits the local level model to $n=800$ simulated points and recovers $\hat\sigma^2_\eta=0.58\pm0.07$ (true 0.5) and $\hat\sigma^2_\varepsilon=0.90\pm0.08$ (true 1.0), and the local linear trend's three variances within their standard errors. The alternative is the **EM algorithm** (ML-6): the E-step is the smoother, the M-step has closed forms for $Q$ and $H$.

### 3.13 Structural time series

A **structural model** writes $y_t=\mu_t+\gamma_t+\varepsilon_t$: a trend (local linear trend), a seasonal (§3.8) and an irregular. Stack the states, $\alpha_t=(\mu_t,\nu_t,\gamma_t,\gamma_{t-1},\dots,\gamma_{t-s+2})$, make $T$ block-diagonal from the components' transition matrices, and set $Z=(1,0,1,0,\dots,0)$. The smoother (3.26) then delivers $a_{t|n}$, whose components are estimates of the trend and the seasonal with standard errors from $P_{t|n}$. This is what "seasonal adjustment" means when it is done with a model. Each component's variance is estimated by ML, and setting one to zero makes that component deterministic.

The connection to ST-7's forecasting folklore is exact. The steady-state filter for the local level model is simple exponential smoothing (§3.10); for the local linear trend it is Holt's linear method; adding the seasonal gives additive Holt–Winters. The Kalman filter supplies what the folklore lacks: a likelihood to choose the smoothing constants, prediction intervals, and the freedom to handle missing values and a non-steady start.

### 3.14 Nonlinear and non-Gaussian extensions; the GP connection

Real systems are nonlinear: $\alpha_{t+1}=f(\alpha_t)+\eta_t$, $y_t=h(\alpha_t)+\varepsilon_t$. The Gaussian conditioning argument no longer closes, because a nonlinear function of a Gaussian is not Gaussian.

**The extended Kalman filter (EKF)** linearizes. With Jacobians (LA-14) $T_t=\partial f/\partial\alpha$ at $a_{t|t}$ and $Z_t=\partial h/\partial\alpha$ at $a_{t|t-1}$,
$$a_{t+1|t}=f(a_{t|t}),\quad P_{t+1|t}=T_tP_{t|t}T_t^\top+Q;\qquad v_t=y_t-h(a_{t|t-1}),\quad F_t=Z_tP_{t|t-1}Z_t^\top+H,\tag{3.28}$$
with gain and update as in (3.22)–(3.23) using $Z_t$. The means go through the true nonlinearity; only the covariances use the linearization. It works when the nonlinearity is mild on the scale of the state uncertainty and fails, sometimes catastrophically, when it is not. For the project's pendulum, $\ddot\theta=-(g/L)\sin\theta$ discretized by semi-implicit Euler with step $\Delta t$,
$$f(\theta,u)=\big(\theta+\Delta t(u-\Delta t\tfrac gL\sin\theta),\ u-\Delta t\tfrac gL\sin\theta\big),\quad T_t=\begin{pmatrix}1-\Delta t^2\frac gL\cos\theta&\Delta t\\-\Delta t\frac gL\cos\theta&1\end{pmatrix}.\tag{3.29}$$
Released from $\theta_0=2$ rad (115°), the small-angle linear filter ($\sin\theta\approx\theta$) has the wrong period and its angle error is 1.48 rad RMS, worse than the raw sensor; the EKF's is 0.054 rad.

The **unscented Kalman filter** replaces the Jacobian with a deterministic set of "sigma points" pushed through $f$ exactly, matching the mean and covariance to second order without derivatives. **Particle filters** drop the Gaussian assumption: represent $p(\alpha_t\mid Y_t)$ by weighted samples, propagate each through the dynamics, reweight by $p(y_t\mid\alpha_t)$, resample. This is ST-6's importance sampling applied sequentially, and the average weight at each step estimates $p(y_t\mid Y_{t-1})$, so the product gives a likelihood estimate. GPS/INS fusion in aircraft and phones is an EKF or UKF; robot SLAM (simultaneous localization and mapping) began as an EKF over robot pose plus landmark positions and moved to particle and graph-optimization methods as maps grew.

**The Gaussian-process connection (ML-4).** An Ornstein–Uhlenbeck process with covariance $k(\tau)=s^2e^{-|\tau|/\ell}$ (the Matérn-½ kernel), sampled at unit spacing, is an AR(1) with $\phi=e^{-1/\ell}$ and innovation variance $s^2(1-\phi^2)$. GP regression with this kernel and Gaussian noise is therefore the RTS smoother on a scalar state, and costs $O(n)$ instead of the $O(n^3)$ of the kernel-matrix formula. Every Matérn kernel with half-integer smoothness $\nu$ is likewise the output of a linear SDE with a $(\nu+\tfrac12)$-dimensional state (Hartikainen and Särkkä, 2010). Exercise C2 verifies the scalar case to machine precision.

### 3.15 Volatility: the stylized facts and ARCH(1)

Daily returns on stocks, indices and exchange rates share four stylized facts. (i) The ACF of returns is near zero at all lags (the mean is close to unpredictable). (ii) The ACF of squared or absolute returns is positive and decays slowly. (iii) The marginal distribution has heavy tails (kurtosis well above 3). (iv) Volatility clusters: large moves, of either sign, follow large moves. Facts (i) and (ii) together say the returns are uncorrelated but not independent: a white noise that is not iid, which ST-7's ARMA cannot represent.

**ARCH(1)** (Engle 1982): $\varepsilon_t=\sigma_tz_t$, $z_t$ iid with mean 0 and variance 1, independent of the past, and
$$\sigma_t^2=\omega+\alpha\varepsilon_{t-1}^2,\qquad\omega>0,\ \alpha\ge0.\tag{3.30}$$

▶ **Try it first.** Show $\varepsilon_t$ is uncorrelated, find its unconditional variance, and write $\varepsilon_t^2$ as an AR(1).

1. $\mathbb E[\varepsilon_t\mid\text{past}]=\sigma_t\mathbb Ez_t=0$, so $\mathbb E\varepsilon_t\varepsilon_{t-h}=\mathbb E[\varepsilon_{t-h}\mathbb E(\varepsilon_t\mid\text{past})]=0$ for $h\ge1$. White noise.
2. Variance: $\mathbb E\varepsilon_t^2=\mathbb E\sigma_t^2=\omega+\alpha\mathbb E\varepsilon_{t-1}^2$. In the stationary state $\mu_2=\omega+\alpha\mu_2$, so $\mu_2=\omega/(1-\alpha)$, finite iff $\alpha<1$.
3. Let $u_t=\varepsilon_t^2-\sigma_t^2=\sigma_t^2(z_t^2-1)$. It has conditional mean zero (a martingale difference), so it is uncorrelated white noise if it has finite variance. Then $\varepsilon_t^2=\omega+\alpha\varepsilon_{t-1}^2+u_t$: **the squared ARCH(1) is an AR(1)** with coefficient $\alpha$, and its ACF is $\alpha^h$.
4. Fourth moment, Gaussian $z$ ($\mathbb Ez^4=3$): $\mu_4=\mathbb E\varepsilon_t^4=3\mathbb E\sigma_t^4=3(\omega^2+2\omega\alpha\mu_2+\alpha^2\mu_4)$. Solve: $\mu_4(1-3\alpha^2)=3\omega^2\big(1+\frac{2\alpha}{1-\alpha}\big)=3\omega^2\frac{1+\alpha}{1-\alpha}$.
5. Divide by $\mu_2^2=\omega^2/(1-\alpha)^2$:
$$\kappa=\frac{\mu_4}{\mu_2^2}=\frac{3(1-\alpha^2)}{1-3\alpha^2}>3,\qquad\text{finite iff }3\alpha^2<1.\tag{3.31}$$

■ **What just happened.** Gaussian shocks with a variance that depends on the past produce a heavy-tailed, uncorrelated, volatility-clustering series: all four stylized facts from one line. The excess kurtosis comes from mixing Normals of different variances (Jensen: $\mathbb E\sigma^4>(\mathbb E\sigma^2)^2$).

### 3.16 GARCH(1,1)

ARCH(1)'s squared returns have ACF $\alpha^h$; data need slow decay with small lag-1 correlation, which ARCH can only get with many lags. Bollerslev (1986) added the lagged variance:
$$\sigma_t^2=\omega+\alpha\varepsilon_{t-1}^2+\beta\sigma_{t-1}^2,\qquad\omega>0,\ \alpha,\beta\ge0.\tag{3.32}$$

▶ **Try it first.** Take expectations of (3.32) in the stationary state; then substitute $\sigma_{t-1}^2=\varepsilon_{t-1}^2-u_{t-1}$.

1. $\bar\sigma^2=\omega+\alpha\bar\sigma^2+\beta\bar\sigma^2$, so
$$\bar\sigma^2=\mathbb E\varepsilon_t^2=\frac{\omega}{1-\alpha-\beta},\qquad\text{finite iff }\alpha+\beta<1.\tag{3.33}$$
$\alpha+\beta$ is the **persistence**.
2. With $u_t=\varepsilon_t^2-\sigma_t^2$ as before: $\varepsilon_t^2-u_t=\omega+\alpha\varepsilon_{t-1}^2+\beta(\varepsilon_{t-1}^2-u_{t-1})$, i.e.
$$\varepsilon_t^2=\omega+(\alpha+\beta)\varepsilon_{t-1}^2+u_t-\beta u_{t-1}.\tag{3.34}$$
**The squared GARCH(1,1) is an ARMA(1,1)** in ST-7's sense, with AR coefficient $\alpha+\beta$ and MA coefficient $-\beta$. Its ACF is ST-7's ARMA(1,1) ACF with those values: $\rho(1)=\frac{\alpha(1-\alpha\beta-\beta^2)}{1-2\alpha\beta-\beta^2}$ and $\rho(h)=(\alpha+\beta)\rho(h-1)$ (exercise B9). With $\alpha$ small and $\alpha+\beta$ near 1, $\rho(1)$ is small and the decay is slow: the empirical pattern.
3. Kurtosis for Gaussian $z$: $\sigma^2_{t+1}=\omega+\sigma_t^2(\alpha z_t^2+\beta)$, and $\mathbb E(\alpha z^2+\beta)^2=3\alpha^2+2\alpha\beta+\beta^2$. Taking expectations of the square, $\mathbb E\sigma^4=\frac{\omega^2+2\omega(\alpha+\beta)\bar\sigma^2}{1-3\alpha^2-2\alpha\beta-\beta^2}$, and $\kappa=3\mathbb E\sigma^4/\bar\sigma^4$ simplifies (as in step 5 of §3.15) to $\kappa=\frac{3[1-(\alpha+\beta)^2]}{1-(\alpha+\beta)^2-2\alpha^2}$, finite iff $(\alpha+\beta)^2+2\alpha^2<1$. For $\alpha=0.1$, $\beta=0.85$: $\kappa=3.774$.

■ **What just happened.** ST-7's ARMA reappears, applied to the squares. Stationarity, the ACF and mean reversion of the variance all come from ST-7 with $\varepsilon_t^2$ in place of $x_t$. The one subtlety: $u_t$ is uncorrelated but heteroskedastic and non-Gaussian, so ARMA estimation on $\varepsilon_t^2$ is inefficient. Use the likelihood.

### 3.17 Estimation, forecasting, and uses

**Quasi-maximum likelihood.** Given $\varepsilon_1,\dots,\varepsilon_n$ and a start $\sigma_1^2$ (the sample variance), (3.32) computes every $\sigma_t^2(\theta)$ recursively, and the Gaussian log-likelihood is
$$\ell(\theta)=-\frac12\sum_{t=1}^n\Big[\log2\pi+\log\sigma_t^2(\theta)+\frac{\varepsilon_t^2}{\sigma_t^2(\theta)}\Big].\tag{3.35}$$
This is a prediction error decomposition again, with the conditional variance playing the role of $F_t$. Maximize it numerically under $\omega>0$, $\alpha,\beta\ge0$, $\alpha+\beta<1$. If $z_t$ is not Normal (it usually has heavier tails), maximizing (3.35) is **quasi**-ML: the estimator is still consistent and asymptotically Normal (Bollerslev–Wooldridge 1992), but the inverse Hessian is the wrong variance and you need ST-3's sandwich $A^{-1}BA^{-1}$, with $A$ the Hessian and $B$ the outer product of scores. The project fits $n=6000$ simulated returns and recovers $\hat\alpha=0.097$, $\hat\beta=0.856$ (true 0.1, 0.85).

**Forecasting volatility.**

▶ **Try it first.** Given information at time $t$, $\sigma_{t+1}^2$ is known. Compute $\mathbb E_t\sigma_{t+2}^2$ using $\mathbb E_t\varepsilon_{t+1}^2=\sigma_{t+1}^2$.

1. For $h\ge2$: $\mathbb E_t\sigma_{t+h}^2=\omega+\alpha\mathbb E_t\varepsilon_{t+h-1}^2+\beta\mathbb E_t\sigma_{t+h-1}^2=\omega+(\alpha+\beta)\mathbb E_t\sigma_{t+h-1}^2$, since $\mathbb E_t\varepsilon_{t+h-1}^2=\mathbb E_t[\sigma_{t+h-1}^2\mathbb E z^2]$.
2. Subtract $\bar\sigma^2=\omega+(\alpha+\beta)\bar\sigma^2$: the deviation $d_h=\mathbb E_t\sigma^2_{t+h}-\bar\sigma^2$ satisfies $d_h=(\alpha+\beta)d_{h-1}$.
3. Iterate from $d_1=\sigma_{t+1}^2-\bar\sigma^2$:
$$\sigma^2_{t+h|t}=\bar\sigma^2+(\alpha+\beta)^{h-1}\big(\sigma^2_{t+1}-\bar\sigma^2\big).\tag{3.36}$$

■ **What just happened.** Volatility forecasts revert geometrically to the long-run level at rate $\alpha+\beta$, with half-life $\log\tfrac12/\log(\alpha+\beta)$: the AR(1) forecast of ST-7, applied to variance. Persistence near 1 (typical daily estimates are 0.97–0.99) means shocks to volatility take months to fade. The forecast of the variance of an $h$-day return is the sum of (3.36) over the horizon.

**Value-at-risk.** A bank's one-day 99% VaR is the loss it exceeds with probability 1%: with conditionally Normal returns, $\mathrm{VaR}_{t+1}=-\big(\mu+\sigma_{t+1|t}\,\Phi^{-1}(0.01)\big)\approx2.33\,\sigma_{t+1|t}$ (in return units, $\mu\approx0$). GARCH makes VaR rise after turbulent days, which is the point and also the criticism (it reacts after the fact). Heavier-tailed $z_t$ (Student $t$) and backtesting the frequency of exceedances are standard.

**Stochastic volatility as state space.** An alternative treats log-variance as a latent AR(1): $h_{t+1}=\mu+\phi(h_t-\mu)+\eta_t$, $y_t=e^{h_t/2}z_t$. This is a nonlinear, non-Gaussian state-space model (3.19 with $Z\alpha$ replaced by a variance), so the Kalman filter does not apply directly. Taking $\log y_t^2=h_t+\log z_t^2$ makes it linear with non-Gaussian noise (a Kalman-filter QML), and a particle filter (§3.14) gives the exact likelihood. GARCH is the special case in which next period's variance is a deterministic function of the past, which is why its likelihood is so cheap.

**The three views.** *Time domain* (ST-7): ACF, ARMA, forecasting; right when the goal is prediction and the dynamics are short-memory. *Frequency domain* (§3.1–3.7): right when the physics is oscillatory or the question is "which frequencies," when you identify a linear system, or when filtering is the operation. *State space* (§3.8–3.14): right when there is a physical hidden state, missing or irregular data, several sensors, time-varying structure, or a need for exact likelihoods and smoothing. They are three coordinate systems on one object: (3.9) translates ARMA to the frequency domain, and §3.8 translates it to state space.

## 4. Worked examples

### Example 4.1: AR(1) spectra with $\phi=\pm0.7$

$\sigma^2=1$. By (3.10), $f(\omega)=1/(1-2\phi\cos2\pi\omega+\phi^2)$ with $1+\phi^2=1.49$.

$\phi=0.7$: at $\omega=0$, $\cos0=1$, $f=1/(1.49-1.4)=1/0.09=11.1111$. At $\omega=\tfrac14$, $\cos\tfrac\pi2=0$, $f=1/1.49=0.6711$. At $\omega=\tfrac12$, $\cos\pi=-1$, $f=1/2.89=0.3460$.

$\phi=-0.7$: $-2\phi\cos\lambda$ flips sign, so $f(0)=0.3460$, $f(\tfrac14)=0.6711$, $f(\tfrac12)=11.1111$.

Check the budget (3.6): $\gamma(0)=1/(1-0.49)=1.961$, the area under each curve. The positive-$\phi$ process has 32 times more spectral density at $\omega=0$ than at $\omega=\tfrac12$ ($11.11/0.346=(1.7/0.3)^2$): low-pass. At $\omega=\tfrac14$ the two agree, as the mirror symmetry requires.

### Example 4.2: The AR(2) peak for $\phi_1=1$, $\phi_2=-0.5$

This is ST-7's complex-root example. By (3.11), $\cos2\pi\omega^*=\frac{1\cdot(-0.5-1)}{4\cdot(-0.5)}=\frac{-1.5}{-2}=0.75$, inside $(-1,1)$, and $\phi_2<0$, so it is a peak. $\omega^*=\arccos(0.75)/2\pi=0.7227/6.2832=0.1150$ cycles per sample, a spectral period $1/\omega^*=8.69$.

The peak height: $D(c)=1+1+0.25-1+2(1)(-1.5)c+2c^2=1.25-3c+2c^2$; $D(0.75)=1.25-2.25+1.125=0.125$, so $f(\omega^*)=8.000$. Compare $f(0)=1/D(1)=1/0.25=4.000$.

From the roots: reciprocal roots have $r=\sqrt{-\phi_2}=0.7071$ and $\cos\theta=\phi_1/(2r)=0.7071$, $\theta=\pi/4$, quasi-period $2\pi/\theta=8.00$. Check step 5 of §3.3: $\cos\theta\cdot\frac{1+r^2}{2r}=0.7071\times\frac{1.5}{1.4142}=0.75$. ✓ The spectral peak period (8.69) exceeds the ACF's quasi-period (8) because $r=0.71$ is far from 1; the peak is broad, only twice $f(0)$.

### Example 4.3: Frequency responses of the 3-point moving average and the first difference

Centred 3-point MA: $\psi(e^{-i\lambda})=\tfrac13(e^{i\lambda}+1+e^{-i\lambda})=\tfrac13(1+2\cos\lambda)$, real (no phase shift, because it is centred).

| $\omega$ | 0 | 1/6 | 1/4 | 1/3 | 1/2 |
|---|---|---|---|---|---|
| $\psi$ (3-pt MA) | 1 | 2/3 | 1/3 | 0 | −1/3 |
| $\lvert\psi\rvert^2$ (3-pt MA) | 1 | 0.4444 | 0.1111 | 0 | 0.1111 |
| $\lvert\psi\rvert^2=4\sin^2\pi\omega$ (difference) | 0 | 1 | 2 | 3 | 4 |

The MA passes the low frequencies, kills $\omega=\tfrac13$ exactly (a period-3 cycle averages to zero over 3 points), and lets 11% of the variance at $\omega=\tfrac12$ through with its sign flipped: the side lobe of (3.12). The difference is the opposite: it removes $\omega=0$ and quadruples the variance at $\omega=\tfrac12$. White noise through the difference has variance $\int_{-1/2}^{1/2}4\sin^2(\pi\omega)\,d\omega=2$, which matches $\mathrm{Var}(w_t-w_{t-1})=2\sigma^2$. The project also shows the 5-point MA, whose $|\psi|^2$ equals $0.04=1/25$ at all four of $\omega=\tfrac16,\tfrac14,\tfrac13,\tfrac12$: the side lobes of (3.12) are not small.

### Example 4.4: The periodogram's variance and a smoothed-spectrum interval

Raw periodogram of Gaussian white noise, $\sigma^2=1$: by Theorem 3.5, $I(\omega_j)\sim\tfrac12\chi^2_2$, i.e. Exponential with mean 1: mean 1, variance 1, for every $n$. The project's 2000 replications at $\omega=\tfrac18$ give variances 1.05, 1.11, 1.03 for $n=128,1024,8192$: no decrease.

Daniell smoothing with $m=2$ averages $L=5$ ordinates, so $\nu=2L=10$. From $\chi^2_{10}$ tables, $\chi^2_{10,0.025}=3.247$ and $\chi^2_{10,0.975}=20.483$. By (3.16), the 95% interval is
$$\Big[\frac{10\hat f}{20.483},\ \frac{10\hat f}{3.247}\Big]=[0.488\hat f,\ 3.080\hat f].$$
If $\hat f(0.2)=2.5$, the interval is $[1.22,7.70]$. It is still wide (a factor of 6.3 from end to end) compared with the raw periodogram's factor of 146 ($0.271I$ to $39.5I$). The project checks that this interval covers the true AR(1) spectrum in 94.8% of 3000 simulated series.

### Example 4.5: AR(2) in state-space form

$x_{t+1}=\phi_1x_t+\phi_2x_{t-1}+w_{t+1}$. State $\alpha_t=(x_t,x_{t-1})^\top$. Then
$$\alpha_{t+1}=\underbrace{\begin{pmatrix}\phi_1&\phi_2\\1&0\end{pmatrix}}_{T}\alpha_t+\underbrace{\begin{pmatrix}1\\0\end{pmatrix}}_{R}\eta_t,\quad\eta_t=w_{t+1},\ Q=\sigma^2;\qquad y_t=\underbrace{(1\ \ 0)}_{Z}\alpha_t,\quad H=0.$$
The first row of $T$ is the AR recursion and the second row copies $x_t$ down. The eigenvalues of $T$ solve $\mu^2-\phi_1\mu-\phi_2=0$; they are the reciprocal roots of $\phi(z)$ (LA-7 and ST-7), so the model is stationary iff they lie inside the unit circle. For $\phi_1=1$, $\phi_2=-0.5$ they are $0.5\pm0.5i$, modulus $0.7071$. $H=0$ means the observation is exact; the filter then reproduces the ST-7 one-step forecasts, and its innovation variance converges to $\sigma^2$ after $p=2$ steps.

### Example 4.6: One predict–update cycle of the constant-velocity tracker

$\Delta t=1$: $T=\begin{pmatrix}1&1\\0&1\end{pmatrix}$, $R=\binom{0.5}{1}$, $Q=1$, $Z=(1\ 0)$, $H=1$. Start from $a_{0|0}=(0,1)^\top$ (at the origin, moving at 1 unit/step), $P_{0|0}=I$. The sensor reads $y_1=2.0$.

*Predict* (3.20): $a_{1|0}=Ta_{0|0}=(0+1,\ 1)^\top=(1,1)^\top$.
$TP_{0|0}T^\top=TT^\top=\begin{pmatrix}2&1\\1&1\end{pmatrix}$ and $RQR^\top=\begin{pmatrix}0.25&0.5\\0.5&1\end{pmatrix}$, so $P_{1|0}=\begin{pmatrix}2.25&1.5\\1.5&2\end{pmatrix}$. Position uncertainty grew from 1 to 2.25, and position and velocity are now correlated (+1.5): if the cart is further along than we think, it is probably faster.

*Update* (3.21)–(3.23): $v_1=2.0-1=1.0$. $F_1=P_{11}+H=3.25$. $K_1=P_{1|0}Z^\top/F_1=(2.25,1.5)^\top/3.25=(0.6923,0.4615)^\top$.
$a_{1|1}=(1,1)^\top+1.0\cdot(0.6923,0.4615)^\top=(1.6923,1.4615)^\top$.
$K_1F_1K_1^\top=\frac{1}{3.25}\begin{pmatrix}2.25^2&2.25\cdot1.5\\2.25\cdot1.5&1.5^2\end{pmatrix}=\begin{pmatrix}1.5577&1.0385\\1.0385&0.6923\end{pmatrix}$, so $P_{1|1}=\begin{pmatrix}0.6923&0.4615\\0.4615&1.3077\end{pmatrix}$.

Read it. The position estimate moved 69% of the way to the measurement because prior position variance (2.25) was 69% of the total (3.25). The **velocity** also moved, by 0.46, though velocity is never measured: the prior correlation carried the information. Velocity variance dropped from 2 to 1.31. The project reproduces every entry, and exercise B5 continues for a second step.

### Example 4.7: Steady-state gain of the local level model, $q=1$

$Q=H=1$. The Riccati fixed point (3.24): $P^2-P-1=0$, $P=\frac{1+\sqrt5}2=1.6180$, the golden ratio. $K=P/(P+1)=1.6180/2.6180=0.6180$, and the filtered variance is $P_{t|t}=PH/(P+H)=0.6180$ as well (so $P=P_{t|t}+Q$, as (3.20) says). Started from a diffuse prior, the gains are $K_1\approx1$, $K_2=2/3$, $K_3=5/8$, $K_4=13/21$, ratios of Fibonacci numbers converging to $0.6180$; the project prints $1,0.6667,0.625$ and $0.618$ at $t=10$. With signal and noise variances equal, each new observation gets 62% of the weight, and simple exponential smoothing with $\alpha=0.618$ is the optimal filter, which the project confirms to $10^{-8}$.

### Example 4.8: GARCH(1,1) with $\omega=0.01$, $\alpha=0.1$, $\beta=0.85$

Unconditional variance (3.33): $\bar\sigma^2=0.01/(1-0.95)=0.2$. Persistence $\alpha+\beta=0.95$, half-life $\log0.5/\log0.95=13.5$ periods. Suppose today's one-step variance is low, $\sigma^2_{t+1}=0.05$. By (3.36) the 10-step forecast is
$$\sigma^2_{t+10|t}=0.2+0.95^{9}(0.05-0.2)=0.2-0.15\times0.6302=0.1055.$$
Nine steps of reversion closed 37% of the gap to 0.2. The project checks the formula against 20,000 simulated paths, and the theoretical kurtosis is 3.774 (sample: 3.58 in $n=6000$; fourth moments converge slowly).

## 5. Why it matters and where it connects

- **ML-4 (Gaussian processes).** §3.14 and C2: a Matérn GP is a linear state-space model, and the RTS smoother (3.26) computes its posterior in $O(n)$. Temporal GPs on long series are done this way.
- **ML-6 (EM).** EM for state-space models uses the smoother as its E-step; the hyperparameters of §3.12 are its M-step.
- **ST-6 (MCMC).** Bayesian state-space models are sampled by Gibbs, drawing the whole state path given parameters (forward filtering, backward sampling, which is (3.26) with draws in place of means) and parameters given the path. Particle filters are ST-6's importance sampling made sequential.
- **CM-1 and F3.** An AR(2) with complex roots is CM-1's damped oscillator sampled in time; (3.11) is its resonance; §3.7 is F3's LTI-system frequency response measured from data.
- **Control (ME).** The Kalman filter is the estimator half of LQG control: the steady-state gain of §3.10 solves a Riccati equation that is the dual of the LQR controller's Riccati equation. Observability (can the state be recovered from $y$?) decides whether the Riccati recursion converges.
- **LA-10, LA-13.** The direct ARMA likelihood is a Toeplitz Cholesky; the Kalman filter is a fast, recursive version of the same factorization (the innovations are the Cholesky-whitened data), and the square-root filter is LA-10 applied to $P$.
- **ST-3.** GARCH's QML standard errors need the sandwich; the Ljung–Box test on innovations is ST-7's.

## 6. Exercises

See [exercises.md](exercises.md): 18 problems. A1–A6 are drills (spectra at given frequencies, a Daniell interval, state-space forms, a steady-state gain, a GARCH forecast). B1–B9 are the module: the random-phase ACF, the $k$-point MA response, the periodogram's distribution, the Slutsky peak, a second Kalman step, innovation independence and the likelihood, the local level as ARIMA(0,1,1), ARCH kurtosis, the GARCH ACF of squares. C1–C3: sensor fusion (modifies `project.py`), GP regression as Kalman smoothing (code), and a particle filter for stochastic volatility (code). Essential: A1, A3, A5, B3, B5, B6, B7, B9.

## 7. Coding project

`project.py` (about 10 seconds) builds every tool from numpy and scipy and checks each against theory. (1) ARMA spectra three ways: closed form (3.9), the Fourier series of the theoretical ACF (equal to $10^{-8}$), and the averaged periodogram of 200 simulations (`spectra.png`); Examples 4.1–4.2. (2) Filters on white noise against $|\psi|^2$ and the Slutsky peak. (3) The periodogram's variance at three sample sizes, the $\chi^2_{10}$ interval's coverage, and the bias–variance trade-off at a sharp AR(2) peak (`bias_variance.png`). (4) Leakage with and without a Hann taper; Fisher's test, with its null calibration. (5) A damped oscillator's frequency response and coherence from the cross-spectrum (`transfer.png`). (6) Generic `kalman_filter` and `rts_smoother`. (7) The tracker: Example 4.6, filtered vs smoothed vs true, innovation whiteness (Ljung–Box), a 50-step gap, and the gain's convergence to Example 4.7 (`tracker.png`, `gain.png`). (8) ML for the local level and local linear trend. (9) The ARMA(1,1) likelihood two ways. (10) A trend + seasonal decomposition. (11) EKF vs linear KF on a pendulum (`ekf.png`). (12) GARCH: stylized facts, QML, the forecast curve against Monte Carlo (`garch.png`). (13) Exercise numbers.

**Build it yourself first:** write `kalman_filter` (about 25 lines) from (3.20)–(3.23) and (3.27), with the missing-data branch, and `rts_smoother` from (3.26), before reading Part 6. Test on Example 4.6 by hand. Then predict, before running Part 7, how the smoothed standard deviation behaves across the 50-step gap (it peaks in the middle; why?). "Your turn": coherence and phase for a delayed common component, two-sensor fusion, and a particle filter for stochastic volatility.

## 8. Where you'll struggle

- **Factors of $2\pi$.** The spectral density is a Fourier *series* in $\omega$ with the ACF as coefficients, not a Fourier transform. Books use $\omega$ in cycles (this module, Shumway & Stoffer), in radians on $[-\pi,\pi]$ with a $1/2\pi$ in front (Brockwell & Davis), or one-sided spectra that double everything. Every factor-of-$2\pi$ error lives here. Keep F3's convention table open; check any formula against two facts: white noise with variance $\sigma^2$ has $f=\sigma^2$ in our convention, and $\int f=\gamma(0)$.
- **Reading structure into the periodogram's noise.** Almost everyone has computed a raw periodogram and found a "significant" peak. Ordinates are exponential: in 500 ordinates of white noise the largest is on average $\sum_{k\le500}1/k\approx6.8$ times the mean. Never interpret an unsmoothed periodogram; smooth it, plot on a log axis with the (3.16) bar, and use Fisher's test (or a smoothed background) before claiming a period.
- **Smoothing span is a choice with consequences.** There is no automatic bandwidth that is right for both a sharp peak and a broad background. Try several spans and say which features survive all of them.
- **The Kalman filter is F5's conditioning formula, and the gain is relative precision.** If (3.22) looks like an incantation, redo §3.9 with scalars ($k=p=1$) first: $K=P/(P+H)$. Then redo Example 4.6 by hand. Students who derive it the third time usually see it; students who memorize the five equations get the transposes wrong and cannot debug their code.
- **Timing conventions.** Some books (and the project's `a0`) start from $a_{1|0}$, others from $a_{0|0}$; some write the transition as $\alpha_t=T\alpha_{t-1}+\eta_t$ with the noise indexed by $t$. Off-by-one errors in the filter usually come from mixing these. Write down which moment you hold at each line of the loop.
- **The state-space form of ARMA is not unique.** Harvey's form (§3.8), Akaike's, and the "innovations form" all give the same likelihood with different states. When your book's $T$ and $R$ differ from these notes, check that the implied ACF of $y$ agrees; that is the only thing that must match.
- **Covariance updates in floating point.** $P-KFK^\top$ can lose symmetry and positive definiteness when measurements are precise (LA-13). Symmetrize, use the Joseph form, or a square-root filter. A filter whose $P$ goes negative reports negative variances and then nonsense.
- **Diffuse starts and the likelihood.** A huge $P_{1|0}$ works but adds a near-constant to (3.27); comparing likelihoods across models with different numbers of diffuse states is wrong without the exact diffuse treatment.
- **The EKF can diverge silently.** If the linearization is poor on the scale of $P$, the filter becomes overconfident, $P$ shrinks, the gain drops, and it stops listening to data. Watch the standardized innovations: if they are not $\mathcal N(0,1)$, the filter is lying.
- **GARCH standard errors.** Returns are heavy-tailed even after dividing by $\sigma_t$, so the inverse-Hessian SEs from (3.35) are too small. Use the sandwich (ST-3). And $\hat\alpha+\hat\beta$ near 1 is common and makes the likelihood flat in one direction.
- **Volatility is never observed.** You cannot plot "true $\sigma_t$" for real data, so every model check is indirect: standardized residuals $\varepsilon_t/\hat\sigma_t$ should have no ARCH left (Ljung–Box on their squares) and forecasts should be evaluated against squared returns, which are a very noisy proxy.

## 9. Check yourself

1. Write the spectral density of white noise, AR(1) and MA(1). Which is low-pass for $\phi>0$, and why is $\phi\to-\phi$ a reflection $\omega\to\tfrac12-\omega$?
2. Why is the inversion formula (3.4) a Fourier-series coefficient formula and not an inverse Fourier transform?
3. State the linear filter theorem and use it to derive the ARMA spectrum in one line.
4. For an AR(2), when does the spectrum have an interior peak, and why is the peak frequency below $\theta/2\pi$?
5. What is the distribution of a periodogram ordinate of Gaussian white noise, and why does that make the periodogram inconsistent?
6. Give the equivalent degrees of freedom of a Daniell smoother of span $L$ and the 95% interval for $f$.
7. What does a taper fix, and what does it cost?
8. Define squared coherence and say what it measures for $y=\psi*x+v$.
9. Derive the Kalman update from the joint distribution of $(\alpha_t,y_t)\mid Y_{t-1}$. Where is the Schur complement?
10. Why are the innovations independent, and how does that give the likelihood (3.27)?
11. Find the steady-state gain of the local level model as a function of $q$ and name the forecasting method it is equivalent to.
12. Derive $\bar\sigma^2$ for GARCH(1,1), show that $\varepsilon_t^2$ is ARMA(1,1), and write the $h$-step variance forecast.

## 10. Pair with

- **Shumway & Stoffer, *Time Series Analysis and Its Applications*, 4th ed., ch. 4 and 6.** The closest match to these notes' conventions (cycles per sample). Ch. 4 has more on spectral estimation, multitapers and coherence with real data; ch. 6 covers state space with R code, EM and the diffuse start.
- **Brockwell & Davis, *Introduction to Time Series and Forecasting*, ch. 4 and 8.** Herglotz's theorem and the periodogram's asymptotics proved carefully (radian convention: convert). Ch. 8 develops state space with a clean treatment of ARMA forms.
- **Durbin & Koopman, *Time Series Analysis by State Space Methods*, 2nd ed., ch. 2–4.** The local level model in full in ch. 2, then the general filter, smoother and likelihood with every derivation; the reference for structural models and diffuse initialization.
- **Särkkä, *Bayesian Filtering and Smoothing* (free PDF), ch. 4–8.** The Bayesian derivation of the filter and smoother, and the clearest account of the EKF, UKF and particle filters. Read it for §3.14.
- **Welch & Bishop, "An Introduction to the Kalman Filter" (free).** A short engineering-flavoured walk-through; useful to see the same equations in control notation.
- **Tsay, *Analysis of Financial Time Series*, 3rd ed., ch. 3.** ARCH, GARCH and their variants on real return data, with diagnostics and the stochastic-volatility comparison.
