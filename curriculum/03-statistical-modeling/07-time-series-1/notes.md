# ST-7: Time Series I: Stationarity, Autocovariance, AR/MA/ARMA, Yule–Walker, Forecasting, ARIMA, Unit Roots

**Time:** ~12 hours · **Prerequisites:** ST-3 (regression, FWL Theorem 3.4, AIC/BIC (3.36)); F5 (covariance matrices Theorem 3.19, the MVN density Theorem 3.21, conditional expectation as projection Theorem 3.29, Markov chains §3.12); LA-7 (eigenvalues, the linear recurrence (3.8)) · **Unlocks:** ST-8

> **Conventions in this module.** A time series is $\{x_t\}$ with $t\in\mathbb Z$ (observed for $t=1,\dots,n$). Mean $\mu=\mathbb E x_t$; **autocovariance** $\gamma(h)=\mathrm{Cov}(x_{t+h},x_t)$; **autocorrelation** $\rho(h)=\gamma(h)/\gamma(0)$. **White noise** $w_t\sim\mathrm{WN}(0,\sigma^2)$ means uncorrelated, mean 0, variance $\sigma^2$; "iid noise" and "Gaussian white noise" are said when meant. The **backshift** is $Bx_t=x_{t-1}$ and $\nabla=1-B$, $\nabla^d=(1-B)^d$, $\nabla_s=1-B^s$. **AR($p$):** $\phi(B)x_t=w_t$ with $\phi(z)=1-\phi_1z-\dots-\phi_pz^p$ (note the minus signs). **MA($q$):** $x_t=\theta(B)w_t$ with $\theta(z)=1+\theta_1z+\dots+\theta_qz^q$ (plus signs; some books, including older Box–Jenkins, use minus). ARMA($p,q$): $\phi(B)x_t=\theta(B)w_t$; ARIMA($p,d,q$): $\phi(B)\nabla^dx_t=\theta(B)w_t$. Hats are sample quantities: $\hat\gamma(h)$ uses divisor $n$, $\hat\rho(h)=\hat\gamma(h)/\hat\gamma(0)$. $\phi_{hh}$ is the PACF at lag $h$. $\Gamma_n=[\gamma(i-j)]_{i,j=1}^n$ is the $n\times n$ autocovariance matrix and $\boldsymbol\gamma_n=(\gamma(1),\dots,\gamma(n))^\top$. The **companion matrix** of an AR($p$) is $A$ in (3.11); its eigenvalues are $\lambda$, the roots of $\phi$ are $z$, and $\lambda=1/z$. $\psi_j$ are the MA($\infty$) weights, $\psi_0=1$. $\hat x_{n+h}$ is the $h$-step forecast from data up to $n$; $P_{n+h}$ its mean squared error. "Stationary" without a qualifier means weakly stationary.

## 0. Where this sits

You need ST-3's regression (OLS, FWL, standard errors and the assumption of independent errors that this module breaks), F5's covariance matrices, the MVN density and conditional expectation as projection, and LA-7's linear recurrences: a system $\mathbf u_k=A\mathbf u_{k-1}$ decays iff every eigenvalue has $|\lambda|<1$. ST-1's method of moments and asymptotic normality of the MLE are used without comment.

This module models dependence in one series indexed by time. ST-8 moves to the frequency domain (the spectral density is the Fourier transform of $\gamma$), to state-space models and the Kalman filter (which computes this module's exact likelihood for any ARMA), and to volatility models (GARCH). ML-8's convolutions are this module's linear filters with learned weights.

## 1. The problem that forced this idea

In 1927 George Udny Yule looked at Wolfer's sunspot numbers, which rise and fall with a period of roughly eleven years but never repeat exactly. The standard approach was harmonic analysis: fit sines with fixed periods and treat the rest as measurement error. Yule pointed out that the irregularity did not look like errors added to a clock. It looked like a pendulum being hit by boys with peashooters: the oscillation itself was disturbed, and each shock changed the phase and amplitude of everything that followed. He modeled the sunspot number as a linear combination of its own two previous values plus a random shock, the first autoregression, and estimated the coefficients by least squares. In the same year Eugen Slutsky showed the converse surprise: take pure random numbers (he used lottery draws), form moving sums, and you get series with apparent cycles that resemble business cycles. Cycles can be manufactured by filtering noise.

Herman Wold (1938) proved that the two ideas cover essentially everything: any weakly stationary process is the sum of a deterministic part and an infinite moving average of white noise. Box and Jenkins (1970) turned the theory into a procedure that engineers and economists could follow: identify a model from the autocorrelations, estimate it, check the residuals, forecast. The ARIMA model and that procedure are still the baseline every forecasting method is compared against.

Then came the failure that made the subject indispensable. Granger and Newbold (1974) simulated pairs of *independent* random walks and regressed one on the other. The textbook $t$ test reported a significant relationship most of the time. Much of the applied macroeconomics of the time regressed trending series on trending series, so this was not a curiosity. Dickey and Fuller (1979) gave a test for whether a series has a unit root; Nelson and Plosser (1982) could not reject one for most US macroeconomic series; Engle and Granger (1987) showed when regressions between random walks are meaningful (cointegration). Engle and Granger shared the 2003 Nobel prize in economics for this and for ARCH.

The forcing problem is ST-3's assumption. Regression standard errors assume the errors are independent. Data indexed by time almost never are: today's temperature, price, or vibration amplitude is close to yesterday's. Ignoring the dependence gives standard errors that are too small, forecasts with intervals that are too narrow, and, when the series wander, relationships that do not exist. The dependence has to be modeled, and the natural models are linear filters driven by white noise.

## 2. Intuition

**What stays the same.** Draw two hundred points of three series. White noise: a fuzzy band with no memory. Daily temperature anomalies: still a band around a fixed level, but smooth, because a warm day tends to be followed by a warm day. A stock price: no band at all, it wanders and its spread grows. The first two are **stationary**: their statistical character does not drift. Concretely, they have a level they return to, a spread that does not grow, and a correlation between today and $h$ days later that depends on $h$ but not on which day you started from. That last property is the whole definition, and you should be able to say it in one breath: *mean constant, autocovariance a function of the lag only*. The third series is not stationary; this module spends its second half on what to do about that.

**Memory as a function.** For a stationary series, the autocorrelation function $\rho(h)$ is its memory: $\rho(1)$ is how much today tells you about tomorrow, $\rho(10)$ how much about ten days out. Everything we do is reading and modeling this function.

**Two ways to build memory from noise.** Let $w_t$ be independent shocks.

- **Moving average:** today is a weighted sum of *recent shocks*, $x_t=w_t+0.5w_{t-1}$. Today and yesterday share the shock $w_{t-1}$, so they are correlated; today and two days ago share nothing. Guess: the ACF of an MA($q$) is exactly zero beyond lag $q$.
- **Autoregression:** today is a weighted sum of *recent days* plus a new shock, $x_t=0.7x_{t-1}+w_t$. Substitute repeatedly: $x_t=w_t+0.7w_{t-1}+0.49w_{t-2}+\dots$, an infinitely long moving average with geometric weights. Every shock is remembered forever but fades. Guess: $\rho(h)=0.7^h$.

**Stability is an eigenvalue condition.** The AR recursion $x_t=\phi_1x_{t-1}+\phi_2x_{t-2}+w_t$ is LA-7's linear recurrence with a forcing term. Stack $(x_t,x_{t-1})$ into a vector and it reads $\mathbf s_t=A\mathbf s_{t-1}+\text{shock}$. If every eigenvalue of $A$ is inside the unit circle, old shocks decay and the process settles into a stationary state; if one is on or outside the circle, a shock never dies or is amplified. Complex eigenvalues mean a damped rotation: the series oscillates with a drifting phase. That is Yule's pendulum and peashooters, and it is also CM-1 §3.8's underdamped oscillator sampled at discrete times and kicked by noise. Guess: an AR(2) is stationary iff the eigenvalues of a $2\times2$ matrix lie in the unit disk, and it looks periodic iff they are complex.

**The PACF is the AR detector.** For an AR(1), $x_t$ and $x_{t-2}$ are correlated ($0.49$), but only *through* $x_{t-1}$. If you already know $x_{t-1}$, $x_{t-2}$ adds nothing. The correlation that remains after removing the intermediate values (ST-3's "controlling for", FWL) is the partial autocorrelation. Guess: the PACF of an AR($p$) is zero beyond lag $p$, the mirror of the MA's ACF cutoff.

**Forecasting is conditional expectation.** Tomorrow's best guess for the AR(1) is $0.7x_n$; the day after, $0.49x_n$; far out, the mean. The uncertainty starts at $\sigma^2$ (one unknown shock) and grows to the series' full variance $\gamma(0)$, because far enough ahead you know nothing beyond the marginal distribution. For a random walk ($\phi=1$) the forecast is $x_n$ forever and the uncertainty grows without bound. That is the practical meaning of non-stationarity.

**Differencing.** A random walk $x_t=x_{t-1}+w_t$ has $x_t-x_{t-1}=w_t$: differencing turns it into noise. The "I" in ARIMA (integrated) is "model the differences, then add them back up". The danger: if you test whether $\phi=1$ with the usual $t$ statistic, you are testing on the boundary where the regressor $x_{t-1}$ itself wanders, and the usual distribution theory fails. And if you regress one wandering series on another, the residuals wander too, the standard errors assume they do not, and the $t$ statistic lies.

## 3. Formal development

### 3.1 Stationarity and the autocovariance function

**Definition 3.1.** $\{x_t\}$ is **strictly stationary** if $(x_{t_1+h},\dots,x_{t_k+h})$ has the same joint distribution as $(x_{t_1},\dots,x_{t_k})$ for all $k$, times and shifts $h$. It is **weakly (second-order) stationary** if $\mathbb Ex_t^2<\infty$, $\mathbb Ex_t=\mu$ for all $t$, and
$$\mathrm{Cov}(x_{t+h},x_t)=\gamma(h)\quad\text{depends on }h\text{ only.}\tag{3.1}$$

Strict implies weak when second moments exist. For Gaussian processes the two coincide, because the mean and covariances determine all joint distributions (F5 Theorem 3.21). Weak stationarity is what every result below uses.

**Theorem 3.2 (Properties of $\gamma$).** For a stationary process: (a) $\gamma(0)\ge0$; (b) $|\gamma(h)|\le\gamma(0)$; (c) $\gamma(-h)=\gamma(h)$; (d) $\gamma$ is **nonnegative definite**: for every $n$ and $\mathbf a\in\mathbb R^n$, $\sum_{i,j}a_ia_j\gamma(i-j)\ge0$, that is, the Toeplitz matrix
$$\Gamma_n=\begin{pmatrix}\gamma(0)&\gamma(1)&\cdots&\gamma(n-1)\\\gamma(1)&\gamma(0)&\cdots&\gamma(n-2)\\\vdots&&\ddots&\vdots\\\gamma(n-1)&\cdots&\gamma(1)&\gamma(0)\end{pmatrix}\succeq0.\tag{3.2}$$

▶ **Try it first.** Each property is one line from something you know: a variance, Cauchy–Schwarz (F5 Theorem 3.17), symmetry of covariance, and the variance of a linear combination.

*Proof.*
1. $\gamma(0)=\mathrm{Var}(x_t)\ge0$.
2. $|\mathrm{Cov}(x_{t+h},x_t)|\le\sqrt{\mathrm{Var}(x_{t+h})\mathrm{Var}(x_t)}=\gamma(0)$ by Cauchy–Schwarz.
3. $\gamma(-h)=\mathrm{Cov}(x_{t-h},x_t)=\mathrm{Cov}(x_t,x_{t-h})$, and with $s=t-h$ this is $\mathrm{Cov}(x_{s+h},x_s)=\gamma(h)$.
4. $0\le\mathrm{Var}(\sum_ia_ix_i)=\sum_{i,j}a_ia_j\mathrm{Cov}(x_i,x_j)=\sum_{i,j}a_ia_j\gamma(i-j)=\mathbf a^\top\Gamma_n\mathbf a$. ∎

■ **What just happened.** $\Gamma_n$ is the covariance matrix of $(x_1,\dots,x_n)$, so it is PSD by F5 Theorem 3.19 (a Gram matrix in LA-8 Theorem 3.16's sense, in the inner product $\mathbb E[UV]$ of centered variables). Stationarity adds the structure: constant along diagonals, a **Toeplitz** matrix. Not every even function with $|\gamma(h)|\le\gamma(0)$ is an autocovariance; (d) is a real constraint (exercise B1 shows $\rho(1)=0.6$, $\rho(h)=0$ for $h\ge2$ is impossible). Toeplitz systems will appear three times: Yule–Walker, the PACF, and the likelihood.

**The basic examples.**
- *White noise* $w_t$: $\gamma(0)=\sigma^2$, $\gamma(h)=0$ for $h\ne0$. Stationary.
- *Moving average* $x_t=w_t+\theta w_{t-1}$: stationary (§3.4).
- *Random walk* $x_t=x_{t-1}+w_t$, $x_0=0$: $x_t=\sum_{s=1}^tw_s$, so
$$\mathbb Ex_t=0,\qquad\mathrm{Var}(x_t)=t\sigma^2,\qquad\mathrm{Cov}(x_t,x_{t+h})=\mathrm{Var}(x_t)=t\sigma^2\ (h\ge0).\tag{3.3}$$
The variance grows with $t$, so (3.1) fails: not stationary. The correlation $\mathrm{corr}(x_t,x_{t+h})=t/\sqrt{t(t+h)}=\sqrt{t/(t+h)}$ is near 1 for $h\ll t$.
- *Random walk with drift* $x_t=\delta+x_{t-1}+w_t$: mean $\delta t$, not stationary either.

### 3.2 Estimating $\gamma$ from one realization

You see one path, not many. Estimating $\mathbb Ex_t$ by a time average requires that the time average converges to the ensemble average: **ergodicity**. For the linear processes of §3.3 with $\sum|\psi_j|<\infty$, it holds in mean square ($\mathrm{Var}(\bar x)=\frac1n\sum_{|h|<n}(1-|h|/n)\gamma(h)\to0$ when $\sum|\gamma(h)|<\infty$); we assume it throughout. The estimators are
$$\bar x=\frac1n\sum_{t=1}^nx_t,\qquad\hat\gamma(h)=\frac1n\sum_{t=1}^{n-h}(x_{t+h}-\bar x)(x_t-\bar x),\qquad\hat\rho(h)=\frac{\hat\gamma(h)}{\hat\gamma(0)}.\tag{3.4}$$
The divisor is $n$, not $n-h$. This biases $\hat\gamma(h)$ toward 0 by a factor about $1-h/n$ (and estimating $\mu$ adds a further $O(1/n)$ bias), but it guarantees that $\hat\Gamma_n$ is PSD: writing $\hat\Gamma_n=\frac1nDD^\top$ with $D$ the $n\times(2n-1)$ matrix whose $i$-th row is the centered data shifted right by $i-1$ places (zeros elsewhere) shows it is a Gram matrix. With $n-h$ the estimated matrix can be indefinite and Yule–Walker can return a nonstationary model.

**How big is $\hat\rho(h)$ when there is no correlation?** This gives the $\pm2/\sqrt n$ bands every software package draws.

▶ **Try it first.** Let $x_t$ be iid with mean 0 (known) and variance $\sigma^2$. For $h\ge1$, compute $\mathbb E[x_{t+h}x_t]$ and $\mathrm{Cov}(x_{t+h}x_t,\ x_{s+h}x_s)$ for $s\ne t$.

1. With the mean known, $\hat\gamma(h)=\frac1n\sum_{t=1}^{n-h}x_{t+h}x_t$. Each term has mean $\mathbb E x_{t+h}\mathbb Ex_t=0$.
2. For $s\ne t$, $\mathbb E[x_{t+h}x_tx_{s+h}x_s]=0$: among the four indices, the smallest one appears once (if $s<t$, $s$ is less than $t$, $s+h$ and $t+h$), and the factor with it has mean zero and is independent of the rest. So the terms are uncorrelated.
3. $\mathrm{Var}(x_{t+h}x_t)=\mathbb E x_{t+h}^2\mathbb Ex_t^2=\sigma^4$. Hence $\mathrm{Var}\hat\gamma(h)=(n-h)\sigma^4/n^2\approx\sigma^4/n$.
4. $\hat\gamma(0)\to\sigma^2$ (law of large numbers), so $\hat\rho(h)\approx\hat\gamma(h)/\sigma^2$ has variance $\approx1/n$. The terms $x_{t+h}x_t$ are $h$-dependent (terms more than $h$ apart are independent), and a CLT holds for such sequences.
5. The same computation with $\hat\gamma(h)$ and $\hat\gamma(k)$, $h\ne k$, gives covariance 0: the products $x_{t+h}x_t$ and $x_{s+k}x_s$ never have all indices paired.

$$\text{Under iid noise:}\quad\sqrt n\,\hat\rho(h)\xrightarrow{d}\mathcal N(0,1),\ \text{independently across }h=1,\dots,K.\tag{3.5}$$

■ **What just happened.** A sample autocorrelation of white noise is an average of $n$ uncorrelated products, so it is $O(1/\sqrt n)$; the $\pm1.96/\sqrt n\approx\pm2/\sqrt n$ band is a 95% band for **one** lag under the null of white noise. With 20 lags you expect one outside by chance. For a non-white process, Bartlett's formula gives $\mathrm{Var}\,\hat\rho(h)\approx\frac1n\sum_{k=1}^\infty[\rho(k+h)+\rho(k-h)-2\rho(h)\rho(k)]^2$ (stated; Brockwell & Davis §2.4); the project uses it to check sample against theoretical ACFs. Mean estimation adds only an $O(1/n)$ bias.

**The random-walk diagnostic.** For a random walk, (3.3) gives $\mathrm{corr}(x_t,x_{t+h})=\sqrt{t/(t+h)}$, close to 1 for every small lag. The sample ACF of one path of length 2000 (project Part 1) is $0.996,0.983,0.965,0.924$ at lags $1,5,10,20$: far above $2/\sqrt n=0.045$ at every lag and decaying almost linearly. This picture means "nonstationary, difference it", not "rich structure at every lag".

### 3.3 Linear processes

**Definition 3.3.** A **linear process** is $x_t=\mu+\sum_{j=-\infty}^\infty\psi_jw_{t-j}$ with $\sum|\psi_j|<\infty$. It is **causal** if $\psi_j=0$ for $j<0$ (depends only on present and past shocks).

Absolute summability makes the sum converge in mean square and almost surely, and lets us exchange sums and expectations.

**Theorem 3.4.** A linear process is stationary with mean $\mu$ and
$$\gamma(h)=\sigma^2\sum_{j=-\infty}^\infty\psi_j\psi_{j+h}.\tag{3.6}$$

▶ **Try it first.** Write $\mathrm{Cov}(x_{t+h},x_t)$ as a double sum and use $\mathbb E[w_aw_b]=\sigma^2\mathbb 1[a=b]$.

1. $\mathbb Ex_t=\mu+\sum\psi_j\mathbb Ew_{t-j}=\mu$.
2. $\mathrm{Cov}(x_{t+h},x_t)=\sum_j\sum_k\psi_j\psi_k\mathbb E[w_{t+h-j}w_{t-k}]$.
3. The expectation is $\sigma^2$ when $t+h-j=t-k$, i.e. $j=k+h$, and 0 otherwise. So $\gamma(h)=\sigma^2\sum_k\psi_{k+h}\psi_k$, independent of $t$. ∎

■ **What just happened.** Covariance at lag $h$ is the overlap between the weight sequence and itself shifted by $h$: an autocorrelation of the filter. Wold's theorem (stated; Brockwell & Davis ch. 2) says every stationary, purely nondeterministic process has this form with causal $\psi$ and uncorrelated $w_t$ (it guarantees $\sum\psi_j^2<\infty$, slightly weaker than the absolute summability of Definition 3.3), so (3.6) covers the whole subject. ST-8 Fourier-transforms it: at frequency $\nu$ in cycles per sample, the spectral density is $f(\nu)=\sum_h\gamma(h)e^{-2\pi i\nu h}=\sigma^2|\psi(e^{-2\pi i\nu})|^2$ (ST-8 (3.3) and (3.8); ST-8 writes $\omega$ for this $\nu$ and has no factors of $2\pi$ in front).

### 3.4 AR(1) and MA(1)

**AR(1).** $x_t=\phi x_{t-1}+w_t$.

▶ **Try it first.** Substitute the equation into itself $k$ times. What happens to the leftover term as $k\to\infty$ when $|\phi|<1$, $|\phi|=1$, $|\phi|>1$?

1. Iterating backward, $x_t=\phi^kx_{t-k}+\sum_{j=0}^{k-1}\phi^jw_{t-j}$.
2. If $|\phi|<1$ and $x$ is stationary (bounded variance), $\mathbb E\big(x_t-\sum_{j<k}\phi^jw_{t-j}\big)^2=\phi^{2k}\gamma(0)\to0$. So
$$x_t=\sum_{j=0}^\infty\phi^jw_{t-j}\qquad(|\phi|<1),\tag{3.7}$$
a causal linear process with $\psi_j=\phi^j$, $\sum|\psi_j|=1/(1-|\phi|)<\infty$. Conversely (3.7) satisfies the recursion, so it is the unique stationary solution.
3. Autocovariance from (3.6): $\gamma(h)=\sigma^2\sum_{j\ge0}\phi^j\phi^{j+h}=\sigma^2\phi^h/(1-\phi^2)$ for $h\ge0$. Or directly: multiply $x_t=\phi x_{t-1}+w_t$ by $x_{t-h}$ ($h\ge1$) and take expectations; $w_t$ is uncorrelated with the past, so $\gamma(h)=\phi\gamma(h-1)$; for $h=0$, $\gamma(0)=\phi^2\gamma(0)+\sigma^2$. Both give
$$\gamma(h)=\frac{\sigma^2\phi^{|h|}}{1-\phi^2},\qquad\rho(h)=\phi^{|h|}.\tag{3.8}$$
4. If $|\phi|=1$, step 1 gives $\mathrm{Var}(x_t-\phi^kx_{t-k})=k\sigma^2$; for a stationary process the left side is at most $4\gamma(0)$, a contradiction for large $k$. No stationary solution: the random walk (3.3).
5. If $|\phi|>1$, $\phi^kx_{t-k}$ explodes backward, but iterating *forward*, $x_t=\phi^{-1}x_{t+1}-\phi^{-1}w_{t+1}=\dots=-\sum_{j\ge1}\phi^{-j}w_{t+j}$, is a stationary solution that depends on the *future*. It is not causal and useless for modeling; started from a fixed value, the forward simulation explodes.

■ **What just happened.** Stationarity of AR(1) is $|\phi|<1$, the stability of the scalar recurrence $u_k=\phi u_{k-1}$. The process is a Gaussian Markov chain on $\mathbb R$ when $w_t$ is Gaussian (F5 §3.12: the transition kernel is $\mathcal N(\phi x,\sigma^2)$, the stationary distribution $\mathcal N(0,\sigma^2/(1-\phi^2))$). It is also the discrete-time sampling of the Ornstein–Uhlenbeck process $dX=-\kappa X\,dt+s\,dW$, the overdamped particle in a harmonic well kicked by noise: sampled every $\Delta$, it is exactly AR(1) with $\phi=e^{-\kappa\Delta}$ (exercise C3 derives the noise variance). The ACF decays geometrically and never reaches zero.

**MA(1).** $x_t=w_t+\theta w_{t-1}$. By (3.6) with $\psi_0=1,\psi_1=\theta$:
$$\gamma(0)=\sigma^2(1+\theta^2),\quad\gamma(\pm1)=\sigma^2\theta,\quad\gamma(h)=0\ (|h|\ge2);\qquad\rho(1)=\frac{\theta}{1+\theta^2}.\tag{3.9}$$
$|\rho(1)|\le\frac12$, with equality at $\theta=\pm1$ (since $1+\theta^2-2|\theta|=(1-|\theta|)^2\ge0$). An MA(1) is stationary for every $\theta$.

**Non-identifiability and invertibility.**

▶ **Try it first.** Replace $\theta$ by $1/\theta$ and $\sigma^2$ by $\theta^2\sigma^2$ in (3.9).

1. $\gamma(0)$ becomes $\theta^2\sigma^2(1+\theta^{-2})=\sigma^2(\theta^2+1)$: unchanged.
2. $\gamma(1)$ becomes $\theta^2\sigma^2\cdot\theta^{-1}=\sigma^2\theta$: unchanged.
3. All higher lags are 0 in both. The two models have identical autocovariances, hence identical Gaussian distributions (a Gaussian process is determined by its mean and $\gamma$).

So the data cannot distinguish $\theta=0.5$ from $\theta=2$. We pick the one in which the shocks can be recovered from the observations. From $w_t=x_t-\theta w_{t-1}$, iterate:
$$w_t=x_t-\theta x_{t-1}+\theta^2x_{t-2}-\dots=\sum_{j\ge0}(-\theta)^jx_{t-j},$$
which converges iff $|\theta|<1$. Such an MA(1) is **invertible**: it has an AR($\infty$) representation, $w_t$ is a function of the observed past, and the forecasting formulas of §3.9 work. Convention: always report the invertible $\theta$.

■ **What just happened.** Identifiability fails because (3.9) depends on $\theta$ only through $\theta/(1+\theta^2)$, which is invariant under $\theta\to1/\theta$; invertibility is the tie-breaker. In code this is the trap: a likelihood optimizer started at the wrong place happily converges to $|\theta|>1$, with the same likelihood value and a different $\sigma^2$. Restrict or reflect.

**MA($q$).** $x_t=\theta(B)w_t$, so $\psi_j=\theta_j$ for $0\le j\le q$ (with $\theta_0=1$) and $\psi_j=0$ otherwise. By (3.6),
$$\gamma(h)=\sigma^2\sum_{j=0}^{q-|h|}\theta_j\theta_{j+|h|}\ (|h|\le q),\qquad\gamma(h)=0\ (|h|>q).\tag{3.10}$$
The ACF **cuts off** after lag $q$. Invertibility: all roots of $\theta(z)$ outside the unit circle (same argument as §3.5 applied to $\theta$).

### 3.5 AR($p$): stationarity is an eigenvalue condition

Write the AR($p$) $x_t=\phi_1x_{t-1}+\dots+\phi_px_{t-p}+w_t$ in state form with $\mathbf s_t=(x_t,x_{t-1},\dots,x_{t-p+1})^\top$:
$$\mathbf s_t=A\mathbf s_{t-1}+w_t\mathbf e_1,\qquad A=\begin{pmatrix}\phi_1&\phi_2&\cdots&\phi_{p-1}&\phi_p\\1&0&\cdots&0&0\\0&1&\cdots&0&0\\\vdots&&\ddots&&\vdots\\0&0&\cdots&1&0\end{pmatrix}.\tag{3.11}$$

**Theorem 3.5 (Causal stationarity of AR($p$)).** The AR($p$) has a unique causal stationary solution iff every eigenvalue of $A$ satisfies $|\lambda|<1$, iff every root of $\phi(z)$ satisfies $|z|>1$. The eigenvalues are exactly the reciprocals of the roots:
$$\det(\lambda I-A)=\lambda^p-\phi_1\lambda^{p-1}-\dots-\phi_p=\lambda^p\,\phi(1/\lambda).\tag{3.12}$$

▶ **Try it first.** Iterate (3.11) backward as you did for AR(1), and use LA-7 (3.8) for the leftover $A^k\mathbf s_{t-k}$. For (3.12), guess an eigenvector: if $A\mathbf v=\lambda\mathbf v$, rows 2 to $p$ say $v_{i-1}=\lambda v_i$.

1. **Eigenvectors.** Rows $2,\dots,p$ of $A\mathbf v=\lambda\mathbf v$ read $v_{i-1}=\lambda v_i$. Setting $v_p=1$ forces $\mathbf v=(\lambda^{p-1},\lambda^{p-2},\dots,\lambda,1)^\top$ (and $v_p=0$ would force $\mathbf v=0$ when $\lambda\ne0$; $\lambda=0$ is not an eigenvalue when $\phi_p\ne0$, since $\det A=\pm\phi_p$).
2. **Row 1** then reads $\phi_1\lambda^{p-1}+\phi_2\lambda^{p-2}+\dots+\phi_p=\lambda\cdot\lambda^{p-1}$, i.e. $\lambda^p-\phi_1\lambda^{p-1}-\dots-\phi_p=0$. Dividing by $\lambda^p$: $1-\phi_1\lambda^{-1}-\dots-\phi_p\lambda^{-p}=\phi(1/\lambda)=0$. So $\lambda$ is an eigenvalue iff $z=1/\lambda$ is a root of $\phi$. To get (3.12) as an identity of polynomials, multiplicities included, expand $\det(\lambda I-A)$ along its first row $(\lambda-\phi_1,-\phi_2,\dots,-\phi_p)$ (cofactor expansion, LA-6): deleting row 1 and column $k$ leaves a block-triangular matrix with $k-1$ diagonal entries $-1$ and $p-k$ diagonal entries $\lambda$, so the cofactor is $(-1)^{1+k}(-1)^{k-1}\lambda^{p-k}=\lambda^{p-k}$, and $\det(\lambda I-A)=(\lambda-\phi_1)\lambda^{p-1}-\sum_{k\ge2}\phi_k\lambda^{p-k}$.
3. **Iterate.** $\mathbf s_t=A^k\mathbf s_{t-k}+\sum_{j=0}^{k-1}A^j\mathbf e_1w_{t-j}$. If the spectral radius $r(A)<1$, then $\|A^k\|\le Cr^k$ for some $r<1$ (LA-7 (3.8) for diagonalizable $A$; LA-11's Jordan form in general), so for a stationary $\mathbf s$ the first term goes to 0 in mean square and $x_t=\sum_{j\ge0}\psi_jw_{t-j}$ with $\psi_j=(A^j)_{11}$, absolutely summable because $|\psi_j|\le Cr^j$. Causal and stationary.
4. **Conversely**, if some eigenvalue has $|\lambda|\ge1$, let $\mathbf u^\top$ be a left eigenvector, $\mathbf u^\top A=\lambda\mathbf u^\top$, and $y_t=\mathbf u^\top\mathbf s_t$. Then $y_t=\lambda y_{t-1}+(\mathbf u^\top\mathbf e_1)w_t$, and $u_1\ne0$ (else $\mathbf u^\top A=\lambda\mathbf u^\top$ forces $\mathbf u=0$ by the shift structure). This is a scalar AR(1) with $|\lambda|\ge1$, which by §3.4 steps 4–5 has no causal stationary solution (for complex $\lambda$ the same variance argument applies to $|y_t|^2$). ∎

■ **What just happened.** The AR($p$) is a $p$-dimensional linear system driven by noise, and stationarity is the stability of that system: LA-7's statement that $A^k\to0$ iff all $|\lambda|<1$. The textbook phrase "roots of $\phi$ outside the unit circle" is the same condition read through $\lambda=1/z$. The project asserts the reciprocity for five AR(2)s. The MA($\infty$) weights are the first entry of $A^j\mathbf e_1$; equivalently they solve $\phi(B)\psi_j=0$ for $j\ge1$ with $\psi_0=1$, the recursion in (3.20) below.

**Yule–Walker equations.**

▶ **Try it first.** Multiply the AR($p$) equation by $x_{t-h}$ for $h\ge1$, take expectations, and use causality: $w_t$ is uncorrelated with $x_{t-h}$.

1. $\mathbb E[x_tx_{t-h}]=\sum_j\phi_j\mathbb E[x_{t-j}x_{t-h}]+\mathbb E[w_tx_{t-h}]$ (mean zero for simplicity).
2. Causality: $x_{t-h}=\sum_k\psi_kw_{t-h-k}$ involves only shocks before $t$, so $\mathbb E[w_tx_{t-h}]=0$ for $h\ge1$, and $=\sigma^2$ for $h=0$ (only $\psi_0w_t$ contributes).
3. Hence
$$\gamma(h)=\sum_{j=1}^p\phi_j\gamma(h-j)\ (h\ge1),\qquad\sigma^2=\gamma(0)-\sum_{j=1}^p\phi_j\gamma(j).\tag{3.13}$$
In matrix form for $h=1,\dots,p$: $\Gamma_p\boldsymbol\phi=\boldsymbol\gamma_p$.

■ **What just happened.** The autocovariance of an AR($p$) satisfies the same recurrence as the process itself, so the ACF is a combination of $\lambda_i^h$: geometric decay, with oscillation when eigenvalues are complex. Read forward, (3.13) computes $\gamma$ from $\boldsymbol\phi$; read backward, it is a $p\times p$ Toeplitz linear system (LA-2) giving $\boldsymbol\phi$ from $\gamma(0),\dots,\gamma(p)$. $\Gamma_p$ is positive definite whenever $\sigma^2>0$, so the solution is unique.

**AR(2) in full.** $\phi(z)=1-\phi_1z-\phi_2z^2$; eigenvalues solve $f(\lambda)=\lambda^2-\phi_1\lambda-\phi_2=0$, so $\lambda_1+\lambda_2=\phi_1$, $\lambda_1\lambda_2=-\phi_2$.

**Theorem 3.6 (The stationarity triangle).** An AR(2) is causal and stationary iff
$$\phi_1+\phi_2<1,\qquad\phi_2-\phi_1<1,\qquad|\phi_2|<1.\tag{3.14}$$
The roots are complex iff $\phi_1^2+4\phi_2<0$, the region under the parabola $\phi_2=-\phi_1^2/4$.

▶ **Try it first.** Evaluate $f(1)$ and $f(-1)$ in terms of the eigenvalues. What does $f(1)>0$ say about where the real roots can be?

1. $f(1)=1-\phi_1-\phi_2=(1-\lambda_1)(1-\lambda_2)$ and $f(-1)=1+\phi_1-\phi_2=(1+\lambda_1)(1+\lambda_2)$.
2. *Complex case.* $\lambda_2=\bar\lambda_1$, so $|\lambda_1|^2=\lambda_1\lambda_2=-\phi_2$; both are inside the circle iff $-\phi_2<1$. Here $f(\pm1)=|1\mp\lambda_1|^2>0$ automatically, and $\phi_2<0$ (since $\phi_1^2+4\phi_2<0$), so (3.14) reduces to $\phi_2>-1$. Consistent.
3. *Real case, necessity.* If $\lambda_1,\lambda_2\in(-1,1)$, then all factors in step 1 are positive and $|\phi_2|=|\lambda_1\lambda_2|<1$.
4. *Real case, sufficiency.* $f(1)>0$ means 1 is not between the roots: both are $<1$ or both $>1$. $f(-1)>0$: both $>-1$ or both $<-1$. "Both $>1$" or "both $<-1$" would give $\lambda_1\lambda_2>1$, contradicting $|\phi_2|<1$. So both lie in $(-1,1)$. ∎

■ **What just happened.** Two real eigenvalues are trapped in $(-1,1)$ by the sign of $f$ at $\pm1$ together with a bound on their product; a complex pair is trapped by its modulus $\sqrt{-\phi_2}$ alone. Each inequality in (3.14) is one edge of a triangle, and each edge is one way to lose stationarity: on $\phi_1+\phi_2=1$ an eigenvalue sits at $+1$ (a unit root, the random-walk boundary of §3.11), on $\phi_2-\phi_1=1$ one sits at $-1$, and on $\phi_2=-1$ a complex pair sits on the unit circle (an oscillation that never damps).

The vertices are $(-2,-1)$, $(2,-1)$, $(0,1)$; the project shades it (figure `ar2_triangle.png`).

**Complex roots and quasi-periodicity.** Write $\lambda=re^{\pm i\omega}$ with $r=\sqrt{-\phi_2}$ and $\cos\omega=\phi_1/(2r)$ (from $\lambda_1+\lambda_2=2r\cos\omega=\phi_1$). The ACF solves the recurrence (3.13), so $\rho(h)=c\lambda^h+\bar c\bar\lambda^h$, i.e.
$$\rho(h)=a\,r^h\cos(\omega h+b),\qquad\text{quasi-period }2\pi/\omega,\quad\omega=\arccos\frac{\phi_1}{2\sqrt{-\phi_2}},\tag{3.15}$$
with $a,b$ fixed by $\rho(0)=1$ and $\rho(1)=\phi_1/(1-\phi_2)$ (the $h=1$ Yule–Walker equation). A damped cosine: Yule's sunspots. The variance (exercise B2) is
$$\gamma(0)=\frac{(1-\phi_2)\,\sigma^2}{(1+\phi_2)\big[(1-\phi_2)^2-\phi_1^2\big]}.\tag{3.16}$$

### 3.6 The partial autocorrelation function

**Definition 3.7.** For $h\ge1$, let $\hat x_t^{(h)}=\phi_{h1}x_{t-1}+\dots+\phi_{hh}x_{t-h}$ be the best linear predictor of $x_t$ from $x_{t-1},\dots,x_{t-h}$ (mean zero). The **partial autocorrelation** at lag $h$ is the last coefficient,
$$\phi_{hh},\qquad\text{where }\Gamma_h(\phi_{h1},\dots,\phi_{hh})^\top=\boldsymbol\gamma_h.\tag{3.17}$$
The linear system is the normal equations (LA-5) of the projection of $x_t$ onto $\mathrm{span}(x_{t-1},\dots,x_{t-h})$ in the inner product $\mathbb E[UV]$. By FWL (ST-3 Theorem 3.4), $\phi_{hh}$ equals the coefficient from regressing the residual of $x_t$ (after removing $x_{t-1},\dots,x_{t-h+1}$) on the residual of $x_{t-h}$ (same removal). By time-reversibility of $\gamma$ the two residuals have equal variance, so $\phi_{hh}$ is their **correlation**: the correlation between $x_t$ and $x_{t-h}$ with the intermediate values partialled out. $\phi_{11}=\rho(1)$.

**Theorem 3.8 (PACF cutoff).** For a causal AR($p$), $\phi_{pp}=\phi_p$ and $\phi_{hh}=0$ for $h>p$.

▶ **Try it first.** For $h>p$, guess the best linear predictor of $x_t$ from $h$ lags, and check that the error is orthogonal to every predictor.

1. Candidate: $\hat x_t=\phi_1x_{t-1}+\dots+\phi_px_{t-p}+0\cdot x_{t-p-1}+\dots+0\cdot x_{t-h}$.
2. Error: $x_t-\hat x_t=w_t$. By causality $w_t$ is uncorrelated with $x_{t-1},\dots,x_{t-h}$ (step 2 of the Yule–Walker derivation).
3. The best linear predictor is characterized by orthogonality of the error to the span (F5 Theorem 3.29; LA-5), and the solution is unique since $\Gamma_h\succ0$. So the candidate is it; its last coefficient is 0 for $h>p$ and $\phi_p$ for $h=p$. ∎
$$\text{AR}(p):\ \phi_{hh}=0\ \text{ for }h>p.\tag{3.18}$$

■ **What just happened.** Once you condition on the $p$ most recent values, an AR($p$) has no further memory: the Markov property of the state $\mathbf s_t$. The PACF cuts off for AR exactly as the ACF does for MA. For an MA($q$) the PACF tails off, because an invertible MA is an AR($\infty$).

**Durbin–Levinson recursion (stated; Brockwell & Davis §2.5).** Solving (3.17) for each $h$ costs $O(h^3)$; exploiting the Toeplitz structure,
$$\phi_{hh}=\frac{\rho(h)-\sum_{j=1}^{h-1}\phi_{h-1,j}\,\rho(h-j)}{1-\sum_{j=1}^{h-1}\phi_{h-1,j}\,\rho(j)},\qquad\phi_{hj}=\phi_{h-1,j}-\phi_{hh}\,\phi_{h-1,h-j}\ (j<h),\tag{3.19}$$
gives all of them in $O(h^2)$, along with the prediction error variance $v_h=v_{h-1}(1-\phi_{hh}^2)$, $v_0=\gamma(0)$. So the PACF at lag $h$ **is** the last Yule–Walker coefficient of the AR($h$) fit. The sample PACF plugs $\hat\rho$ into (3.19); under an AR($p$) model, $\sqrt n\,\hat\phi_{hh}\approx\mathcal N(0,1)$ for $h>p$, so the same $\pm2/\sqrt n$ bands apply. The project asserts (3.19) against direct Toeplitz solves of (3.17) for four models.

### 3.7 ARMA($p,q$)

**Definition 3.9.** $\phi(B)x_t=\theta(B)w_t$, with $\phi$ and $\theta$ having no common roots. It is **causal** if all roots of $\phi$ lie outside the unit circle (Theorem 3.5's argument, now with the MA part as forcing), and **invertible** if all roots of $\theta$ do.

**$\psi$-weights.** For a causal ARMA, $x_t=\psi(B)w_t$ with $\psi(z)=\theta(z)/\phi(z)$. Equate coefficients in $\phi(z)\psi(z)=\theta(z)$:
$$\psi_j=\theta_j+\sum_{k=1}^{p}\phi_k\psi_{j-k},\qquad\psi_0=1,\ \theta_j=0\ (j>q),\ \psi_{j}=0\ (j<0).\tag{3.20}$$
For $j>q$ this is the AR recurrence again, so the $\psi_j$ (and $\gamma(h)$ for $h>q$) decay like $|\lambda_{\max}|^j$.

**ARMA(1,1)** $x_t=\phi x_{t-1}+w_t+\theta w_{t-1}$.

▶ **Try it first.** Apply (3.20) with $p=q=1$. Then sum (3.6) using geometric series.

1. $\psi_0=1$; $\psi_1=\theta+\phi\psi_0=\theta+\phi$; for $j\ge2$ the $\theta_j$ term is gone and $\psi_j=\phi\psi_{j-1}$. So for all $j\ge1$:
$$\psi_j=(\phi+\theta)\phi^{j-1}\ (j\ge1).\tag{3.21}$$
2. $\gamma(0)=\sigma^2\big[1+(\phi+\theta)^2\sum_{j\ge1}\phi^{2(j-1)}\big]=\sigma^2\Big[1+\frac{(\phi+\theta)^2}{1-\phi^2}\Big]=\sigma^2\frac{1+2\phi\theta+\theta^2}{1-\phi^2}$.
3. $\gamma(1)=\sigma^2\big[\psi_0\psi_1+\sum_{j\ge1}\psi_j\psi_{j+1}\big]=\sigma^2\Big[(\phi+\theta)+\frac{(\phi+\theta)^2\phi}{1-\phi^2}\Big]=\sigma^2\frac{(\phi+\theta)(1+\phi\theta)}{1-\phi^2}$.
4. For $h\ge2$, multiply the model by $x_{t-h}$: both $w_t$ and $w_{t-1}$ are uncorrelated with $x_{t-h}$, so $\gamma(h)=\phi\gamma(h-1)$.
$$\rho(1)=\frac{(\phi+\theta)(1+\phi\theta)}{1+2\phi\theta+\theta^2},\qquad\rho(h)=\phi^{h-1}\rho(1)\ (h\ge1).\tag{3.22}$$

■ **What just happened.** The MA part sets the starting value $\rho(1)$; the AR part sets the decay rate thereafter. Neither the ACF nor the PACF cuts off. If $\theta=-\phi$, then $\psi_j=0$ for $j\ge1$ and $x_t=w_t$: the factor $(1-\phi B)$ on both sides cancels (**parameter redundancy**). This is why the definition demands no common roots, and why overfitting an ARMA($p+1,q+1$) to an ARMA($p,q$) gives an unstable, nearly cancelling pair of parameters with huge standard errors.

**The identification table.**

| | ACF | PACF |
|---|---|---|
| AR($p$) | tails off (damped exponential/cosine) | **cuts off after $p$** |
| MA($q$) | **cuts off after $q$** | tails off |
| ARMA($p,q$) | tails off after lag $q-p$ | tails off after lag $p-q$ |
| Random walk / unit root | decays slowly, nearly linearly, from near 1 | $\phi_{11}\approx1$ |

The project's figure `identification_table.png` draws this for nine simulated processes.

### 3.8 Estimation and diagnostics

Assume mean zero (subtract $\bar x$; the effect on everything below is $O(1/n)$).

**Yule–Walker** (ST-1's method of moments): replace $\gamma$ by $\hat\gamma$ in (3.13),
$$\hat{\boldsymbol\phi}=\hat\Gamma_p^{-1}\hat{\boldsymbol\gamma}_p,\qquad\hat\sigma^2=\hat\gamma(0)-\hat{\boldsymbol\phi}^\top\hat{\boldsymbol\gamma}_p.\tag{3.23}$$
Because $\hat\Gamma_{p+1}$ is positive definite (divisor $n$), the fitted model is always causal (Brockwell & Davis §5.1). For MA or ARMA the moment equations are nonlinear and inefficient (exercise B5); use likelihood.

**Conditional least squares (CLS):** regress $x_t$ on $x_{t-1},\dots,x_{t-p}$ for $t=p+1,\dots,n$ by OLS. This maximizes the likelihood *conditional* on the first $p$ values. It can return a nonstationary fit (exercise B6).

**Exact Gaussian likelihood.** Any joint density factors into one-step predictive densities, $p(x_1,\dots,x_n)=\prod_{t=1}^np(x_t\mid x_1,\dots,x_{t-1})$. For a Gaussian process each factor is Normal with mean the best linear predictor $\hat x_t$ and variance its MSE $P_t$, so
$$-2\log L=\sum_{t=1}^n\Big[\log(2\pi P_t)+\frac{(x_t-\hat x_t)^2}{P_t}\Big]\tag{3.24}$$
(the **innovations form**). Equivalently, $(x_1,\dots,x_n)\sim\mathcal N(\mathbf 0,\Gamma_n)$ and $L$ is F5 Theorem 3.21's density with the Toeplitz covariance (3.2). The two are the same number: (3.24) is what you get from a Cholesky factorization $\Gamma_n=LL^\top$ (LA-10 Theorem 3.5), whose rows produce the innovations. For a general ARMA, $\hat x_t$ and $P_t$ come from the Kalman filter (ST-8). For AR($p$) they are explicit: for $t>p$, $\hat x_t=\sum\phi_jx_{t-j}$, $P_t=\sigma^2$; the first $p$ values are $\mathcal N(\mathbf 0,\Gamma_p)$.

▶ **Try it first.** Write the AR(1) exact likelihood using $x_1\sim\mathcal N(0,\sigma^2/(1-\phi^2))$ and $x_t\mid x_{t-1}\sim\mathcal N(\phi x_{t-1},\sigma^2)$.

1. $p(x_1)=\big(2\pi\sigma^2/(1-\phi^2)\big)^{-1/2}\exp\{-(1-\phi^2)x_1^2/(2\sigma^2)\}$.
2. $p(x_t\mid x_{1:t-1})=p(x_t\mid x_{t-1})=(2\pi\sigma^2)^{-1/2}\exp\{-(x_t-\phi x_{t-1})^2/(2\sigma^2)\}$ by the Markov property.
3. Multiply and take logs:
$$\log L(\phi,\sigma^2)=-\frac n2\log(2\pi\sigma^2)+\frac12\log(1-\phi^2)-\frac{(1-\phi^2)x_1^2+\sum_{t=2}^n(x_t-\phi x_{t-1})^2}{2\sigma^2}.\tag{3.25}$$

■ **What just happened.** The exact likelihood is the conditional (CLS) likelihood plus the density of the first observation under the stationary distribution. That term contributes $\frac12\log(1-\phi^2)$, which goes to $-\infty$ as $|\phi|\to1$: the exact MLE is always stationary. Comparing (3.25) with the MVN density shows $\det\Gamma_n=\sigma^{2n}/(1-\phi^2)$ and that $\Gamma_n^{-1}$ is tridiagonal (an AR(1) is Markov, so its precision matrix has zeros beyond the first off-diagonal; this sparsity is what makes ST-8's Kalman filter cost $O(n)$). The project checks (3.25) against the 200-dimensional Toeplitz MVN evaluated by Cholesky.

**Asymptotic normality.** For a causal AR($p$), Yule–Walker, CLS and exact ML all satisfy
$$\sqrt n(\hat{\boldsymbol\phi}-\boldsymbol\phi)\xrightarrow{d}\mathcal N(\mathbf 0,\ \sigma^2\Gamma_p^{-1}).\tag{3.26}$$
*Heuristic for AR(1) CLS.* $\hat\phi-\phi=\sum x_{t-1}w_t/\sum x_{t-1}^2$. The numerator is a sum of uncorrelated terms (a martingale difference sequence) with variance $n\sigma^2\gamma(0)$; the denominator is $\approx n\gamma(0)$. So $\mathrm{Var}(\hat\phi)\approx\sigma^2/(n\gamma(0))=(1-\phi^2)/n$, and a martingale CLT gives normality. Note that the precision improves as $|\phi|\to1$, where the regressor varies more. In the project (AR(2) $\phi=(1,-0.5)$, $n=200$, 500 replications), the asymptotic sd is 0.0612 for each coefficient; CLS and ML give sds 0.0622/0.0608 and means within 0.012 of the truth; Yule–Walker has a slightly larger finite-sample bias (mean $\hat\phi_1=0.978$) from the divisor-$n$ shrinkage. ML and CLS differ by 0.004 on average.

**Model selection.** AIC $=-2\log L+2k$ and BIC $=-2\log L+k\log n$ with $k=p+q+1$ parameters (ST-3 (3.36), same reasoning). Compare models fitted on the same data span.

**Residual diagnostics and Ljung–Box.** If the model is right, the residuals $\hat w_t$ are approximately white noise and (3.5) applies to their ACF. Summing squares of approximately independent $\mathcal N(0,1)$ variables gives a portmanteau statistic: Box–Pierce $n\sum_{h=1}^K\hat\rho^2(h)\approx\chi^2_K$. The exact variance of $\hat\rho(h)$ under iid noise is closer to $(n-h)/(n(n+2))$, which motivates the **Ljung–Box** correction
$$Q=n(n+2)\sum_{h=1}^K\frac{\hat\rho_{\hat w}^2(h)}{n-h}\ \approx\ \chi^2_{K-p-q}\ \text{ under a correct ARMA}(p,q).\tag{3.27}$$
The loss of $p+q$ degrees of freedom is because fitting forces the first few residual autocorrelations to be small (stated; Box & Pierce 1970). The project simulates the size: with a correctly fitted AR(2), $K=10$, nominal 5%, the rejection rate is 0.052; fitting AR(1) to the same data rejects every time.

### 3.9 Forecasting

**Best linear predictor.** Given $x_1,\dots,x_n$ (mean zero), the best linear predictor of $x_{n+h}$ is $\hat x_{n+h}=\sum_{j=1}^na_jx_{n+1-j}$ with the error orthogonal to every $x_i$ (F5 Theorem 3.29; LA-5 projection):
$$\Gamma_n\mathbf a=\big(\gamma(h),\gamma(h+1),\dots,\gamma(h+n-1)\big)^\top.\tag{3.28}$$
For Gaussian processes, linear prediction is the conditional expectation $\mathbb E[x_{n+h}\mid x_{1:n}]$ (F5 Theorem 3.24), so nothing nonlinear can do better.

**AR(1).**

▶ **Try it first.** Iterate $x_{n+h}=\phi x_{n+h-1}+w_{n+h}$ down to $x_n$. Which terms are known at time $n$?

1. $x_{n+h}=\phi^hx_n+\sum_{j=0}^{h-1}\phi^jw_{n+h-j}$.
2. The future shocks $w_{n+1},\dots,w_{n+h}$ are uncorrelated with $x_1,\dots,x_n$, so the prediction is the first term and the error is the sum: with mean $\mu$,
$$\hat x_{n+h}=\mu+\phi^h(x_n-\mu),\qquad P_{n+h}=\sigma^2\sum_{j=0}^{h-1}\phi^{2j}=\sigma^2\frac{1-\phi^{2h}}{1-\phi^2}.\tag{3.29}$$
3. As $h\to\infty$, $\hat x_{n+h}\to\mu$ exponentially fast and $P_{n+h}\to\sigma^2/(1-\phi^2)=\gamma(0)$.

**ARMA via $\psi$-weights.** For a causal invertible ARMA with the infinite past known, the shocks $w_s$, $s\le n$, are known (invertibility), and $x_{n+h}=\sum_{j=0}^{h-1}\psi_jw_{n+h-j}+\sum_{j\ge h}\psi_jw_{n+h-j}$. The first sum is future noise, the second is known. So
$$\hat x_{n+h}=\sum_{j\ge h}\psi_jw_{n+h-j},\qquad P_{n+h}=\sigma^2\sum_{j=0}^{h-1}\psi_j^2.\tag{3.30}$$
In practice you iterate the ARMA recursion with future $w$'s set to 0 and past $w$'s replaced by residuals. With only $x_1,\dots,x_n$ available, residuals before $t=1$ are set to 0 (the **truncated** predictor); the error this introduces decays like $|z_{\min}|^{-n}$ for an invertible model, where $z_{\min}$ is the root of $\theta(z)$ closest to the unit circle; it is negligible except for short series or near-unit MA roots (the exact finite-sample version is the Kalman filter, ST-8).

■ **What just happened.** Forecasting is deleting the unknown future shocks. The forecast reverts to the mean at the rate of the largest eigenvalue; the error variance accumulates one $\psi_j^2$ per step and saturates at $\gamma(0)$, since $\sum_{j\ge0}\psi_j^2\sigma^2=\gamma(0)$ by (3.6). A 95% interval is $\hat x_{n+h}\pm1.96\sqrt{P_{n+h}}$ under Gaussian shocks. It ignores parameter uncertainty, an $O(1/n)$ effect: the project's empirical coverage for a fitted AR(2) with $n=200$ is between 0.92 and 0.96 for $h=1,\dots,10$.

### 3.10 Nonstationarity: trends, differencing, ARIMA

Two models for a series that trends upward:
- **Trend-stationary:** $x_t=a+bt+y_t$ with $y_t$ stationary. Shocks are temporary; the series returns to the line. Remove the trend by regression (ST-3).
- **Difference-stationary:** $x_t=\delta+x_{t-1}+y_t$. Shocks are permanent; the "trend" is $\delta t$ plus accumulated shocks. Remove it by differencing: $\nabla x_t=\delta+y_t$.

They forecast differently. For a random walk with drift ($y_t=w_t$), $x_{n+h}=x_n+h\delta+\sum_{j=1}^hw_{n+j}$, so
$$\hat x_{n+h}=x_n+h\hat\delta,\qquad P_{n+h}=h\sigma^2.\tag{3.31}$$
The trend model gives a constant-width interval around the fitted line. Fit the trend model to a random walk with drift and its intervals are badly too narrow at long horizons: in the project ($n=200$, $\delta=0.2$, $h=20$, 1000 replications) the trend-regression interval covers 0.638 of the time, the ARIMA interval 0.930.

**Definition 3.10.** $x_t$ is **ARIMA($p,d,q$)** if $\nabla^dx_t$ is a causal invertible ARMA($p,q$):
$$\phi(B)(1-B)^dx_t=\theta(B)w_t.\tag{3.32}$$
The AR polynomial of $x_t$ itself, $\phi(z)(1-z)^d$, has a root at $z=1$ of multiplicity $d$: a **unit root**. Forecast $\nabla^dx$ as an ARMA, then sum back. The $\psi$-weights of $x$ (from (3.20) with $\phi(z)(1-z)^d$ in place of $\phi$) do not decay, so $P_{n+h}$ grows without bound (exercise B8).

**Over-differencing.** Differencing a stationary series introduces an MA unit root. If $x_t=a+bt+w_t$, then $\nabla x_t=b+w_t-w_{t-1}$, an MA(1) with $\theta=-1$: non-invertible, variance doubled ($2\sigma^2$ instead of $\sigma^2$), and lag-1 autocorrelation $-0.5$. A sample ACF of the differenced series with $\hat\rho(1)$ near $-0.5$ is the signature. Difference only as often as needed, usually $d\le1$ (occasionally 2), and prefer a test (§3.11) to eyeballing.

**Seasonality.** Monthly data with a yearly pattern use the seasonal difference $\nabla_{12}=1-B^{12}$ and seasonal ARMA polynomials in $B^{12}$: SARIMA$(p,d,q)\times(P,D,Q)_s$, $\Phi(B^s)\phi(B)\nabla_s^D\nabla^dx_t=\Theta(B^s)\theta(B)w_t$. The standard example is the "airline model" $(0,1,1)\times(0,1,1)_{12}$ of Box and Jenkins for log monthly airline passengers: $\nabla\nabla_{12}\log x_t=(1+\theta B)(1+\Theta B^{12})w_t$, whose ACF is nonzero only at lags 1, 11, 12, 13. Project "Your turn" 1 builds one.

> **The Box–Jenkins procedure.**
> 1. **Identify.** Plot the series. Transform (log) if the variance grows with the level. Difference until stationary (ACF decays quickly; unit-root test). Read candidate $(p,q)$ from the ACF/PACF table.
> 2. **Estimate.** Fit by exact ML (or CLS for AR). Check causality, invertibility, no near-cancelling roots.
> 3. **Diagnose.** Residual ACF inside $\pm2/\sqrt n$ bands, Ljung–Box (3.27), residual normality. Compare candidates by AIC/BIC. If the diagnostics fail, return to 1.
> 4. **Forecast.** Point forecasts and intervals by (3.30), undoing the differencing and the transformation.

### 3.11 Unit roots and the Dickey–Fuller test

Is $\phi=1$? Write the AR(1) as a regression of the difference on the lagged level,
$$\nabla x_t=\alpha+\varrho\,x_{t-1}+e_t,\qquad\varrho=\phi-1,\qquad H_0:\varrho=0\ \text{(unit root)},\ H_1:\varrho<0.\tag{3.33}$$
The natural statistic is the OLS $t$ ratio $\tau=\hat\varrho/\mathrm{se}(\hat\varrho)$. Under $H_0$ it does not have a $t$ distribution, even asymptotically.

▶ **Try it first.** Under $H_0$ (no constant, $x_0=0$), $\hat\phi-1=\sum x_{t-1}w_t/\sum x_{t-1}^2$. Find the order of magnitude in $n$ of the numerator and the denominator, using (3.3).

1. **Denominator.** $\mathbb E\sum_{t=1}^nx_{t-1}^2=\sigma^2\sum_{t=1}^n(t-1)\approx\sigma^2n^2/2$. It is $O_p(n^2)$, not $O_p(n)$ as for stationary $x$. ($Y_n=O_p(a_n)$ means $Y_n/a_n$ stays bounded in probability: "of size $a_n$", for random quantities.) More precisely, with $x_{\lfloor sn\rfloor}\approx\sigma\sqrt nW(s)$ for a standard Brownian motion $W$ (Donsker's theorem, stated), $n^{-2}\sum x_{t-1}^2\to\sigma^2\int_0^1W(s)^2ds$.
2. **Numerator.** Square $x_t=x_{t-1}+w_t$: $x_t^2=x_{t-1}^2+2x_{t-1}w_t+w_t^2$. Sum over $t$ and telescope:
$$\sum_{t=1}^nx_{t-1}w_t=\tfrac12\Big(x_n^2-\sum_{t=1}^nw_t^2\Big).$$
So $n^{-1}\sum x_{t-1}w_t\to\frac{\sigma^2}2(W(1)^2-1)$: $O_p(n)$.
3. **Ratio.**
$$n(\hat\phi-1)\xrightarrow{d}\frac{\frac12(W(1)^2-1)}{\int_0^1W(s)^2ds},\qquad\tau\xrightarrow{d}\frac{\frac12(W(1)^2-1)}{\sqrt{\int_0^1W(s)^2ds}}.\tag{3.34}$$
(For $\tau$, the standard error is $\hat\sigma/\sqrt{\sum x_{t-1}^2}=O_p(1/n)$.)

■ **What just happened.** Three things differ from the stationary case. (i) $\hat\phi$ converges at rate $n$, not $\sqrt n$ (**superconsistency**): a random-walk regressor has variance growing like $t$, so it is very informative about its own coefficient. (ii) The limit is not Normal: it is a ratio of functionals of one Brownian path. (iii) It is skewed left: $W(1)^2\sim\chi^2_1$ is below its mean 1 with probability 0.68, so the numerator is usually negative, and $\hat\phi$ is biased below 1. The derivation of (3.34) needs the functional CLT (Hamilton ch. 17 has it in full); what you need here is the conclusion that the $t$ table is wrong.

**The Dickey–Fuller distribution.** Critical values (asymptotic; Fuller 1976, refined by MacKinnon):

| regression (3.33) includes | 1% | 5% | 10% |
|---|---|---|---|
| no constant | $-2.58$ | $-1.95$ | $-1.62$ |
| constant $\alpha$ | $-3.43$ | $-2.86$ | $-2.57$ |
| constant and trend $\beta t$ | $-3.96$ | $-3.41$ | $-3.12$ |

Compare the $t$ (or Normal) 5% one-sided value $-1.645$. The project simulates $10^4$ random walks of length 100: the 5% quantile of $\tau$ is $-1.94$ without a constant and $-2.89$ with one, against $-1.66$ for $t_{98}$. Including a constant (needed whenever the alternative has a nonzero mean) shifts the distribution further left, and a trend further still. Using the $t$ table rejects a true unit root far too often.

**Augmented Dickey–Fuller (ADF).** If $\nabla x_t$ is itself autocorrelated (an ARIMA($p,1,0$)), the errors in (3.33) are not white. Add lagged differences: $\nabla x_t=\alpha+\varrho x_{t-1}+\sum_{j=1}^k\beta_j\nabla x_{t-j}+e_t$. The $t$ statistic on $\varrho$ has the same limiting distribution (3.34), so the same table applies. Choose $k$ by AIC or a rule like $k\approx(n-1)^{1/3}$.

**Low power.** Near the null, a stationary AR(1) with $\phi=0.95$ and a random walk look alike over 100 observations. The project's power (5% test with constant, $n=100$): $0.059$ at $\phi=0.99$, $0.082$ at $0.97$, $0.129$ at $0.95$, $0.34$ at $0.9$, $0.873$ at $0.8$. "Failed to reject a unit root" is weak evidence for one. Tests with stationarity as the null (KPSS) exist for the reverse framing.

**The economic stakes.** Nelson and Plosser (1982) applied Dickey–Fuller tests to fourteen long annual US macroeconomic series (real GNP, employment, prices, wages and others) and could not reject a unit root for any except the unemployment rate. If real output is difference-stationary, recessions are permanent losses that are never made up, the "business cycle" is not a cycle around a fixed trend, and forecast uncertainty grows without bound; if trend-stationary, output returns to its trend. Perron (1989) showed that allowing a one-time break in the trend (the 1929 crash, the 1973 oil shock) reverses many of the conclusions. Given the power numbers above, the data alone settle little; the debate is a lesson in what unit-root tests can and cannot tell you.

### 3.12 Spurious regression and cointegration

**Theorem 3.11 (Spurious regression, heuristic).** Let $x_t,y_t$ be independent random walks and regress $y_t=a+bx_t+e_t$ by OLS. Then $\hat b$ does not converge to 0, and the OLS $t$ statistic for $b=0$ grows like $\sqrt n$:
$$P(|t_{\hat b}|>2)\to1\qquad(n\to\infty).\tag{3.35}$$

▶ **Try it first.** By (3.3), $x_t$ and $y_t$ are of size $\sqrt t$. What is the size of $\sum(x_t-\bar x)^2$, of the residuals, and of $s^2=\sum e_t^2/(n-2)$?

1. $\sum(x_t-\bar x)^2=O_p(n^2)$ and $\sum(x_t-\bar x)(y_t-\bar y)=O_p(n^2)$ (step 1 of §3.11, with Brownian limits). So $\hat b=O_p(1)$: it converges to a ratio of functionals of two independent Brownian motions, a nondegenerate random variable, not to 0.
2. The residual $e_t=y_t-\hat a-\hat bx_t$ is a combination of two random walks: it is itself $O_p(\sqrt n)$ and wanders. So $s^2=O_p(n)$, not $O_p(1)$.
3. OLS reports $\mathrm{se}(\hat b)=s/\sqrt{\sum(x_t-\bar x)^2}=O_p(\sqrt n/n)=O_p(n^{-1/2})$.
4. $t=\hat b/\mathrm{se}(\hat b)=O_p(1)/O_p(n^{-1/2})=O_p(\sqrt n)$, and $\sqrt n\cdot(\text{a nonzero random variable})$ exceeds 2 in absolute value with probability tending to 1. ∎

■ **What just happened.** The OLS standard error formula assumes the errors are uncorrelated with constant variance. Here the "errors" are a random walk, extremely autocorrelated, so the formula understates the real sampling variability of $\hat b$ by a factor of order $\sqrt n$. The project's simulation, $10^4$ pairs of independent random walks of length 50, finds $|t|>2$ in 67% of regressions; at $n=100,200,400$ (4000 pairs each, Part 11) the fraction is 0.75, 0.83, 0.88. Two independent stationary AR(1)s with $\phi=0.5$ give 12% (still inflated, because the errors are autocorrelated; ST-3's sandwich idea extended to time series is the Newey–West HAC standard error); two iid series give 5.2%. Granger and Newbold's rule of thumb: an $R^2$ higher than the Durbin–Watson statistic signals a spurious regression. Regression on time series without checking for unit roots is malpractice. The fix is to difference both (regress $\nabla y$ on $\nabla x$) unless they are cointegrated.

**Cointegration.** Two I(1) series (each a unit-root process) are **cointegrated** if some linear combination $y_t-\beta x_t=u_t$ is stationary. Example: $x_t$ a random walk, $y_t=1+2x_t+u_t$ with $u_t$ AR(1). Both wander, but they wander together; $u_t$ is the stationary deviation from a long-run equilibrium (prices of the same commodity in two markets, consumption and income, a stock's spot and futures prices).

*Engle–Granger two-step.* (1) Regress $y_t$ on $x_t$ by OLS; under cointegration $\hat\beta$ is superconsistent (in the project, mean $\hat\beta=1.9993$ for $\beta=2$). (2) Apply a DF test to the residuals $\hat u_t$. The critical values are **not** the DF table: OLS chooses $\hat\beta$ to make $\hat u$ as small as possible, so the residuals look more stationary than a fixed combination would. Simulated under no cointegration ($n=100$), the 5% critical value is $-3.39$ (MacKinnon's asymptotic value for two variables is $-3.34$); using the DF value $-2.89$ gives a size of 0.149 instead of 0.05. On the cointegrated pair the test rejects in essentially every replication.

*Error-correction form (Granger representation theorem, stated; Engle & Granger 1987).* Cointegrated series satisfy
$$\nabla y_t=\alpha\,(y_{t-1}-\beta x_{t-1})+\gamma\,\nabla x_t+\varepsilon_t,\qquad\alpha<0,\tag{3.36}$$
plus possible lagged differences (in this equation only, $\alpha$ and $\gamma$ are regression coefficients, not the constant of (3.33) or the autocovariance): changes in $y$ respond to last period's disequilibrium, pulling the pair back together. Exercise C2 derives it for the AR(1)-deviation example. With more than two series, Johansen's method estimates several cointegrating vectors (Hamilton ch. 20).

### 3.13 Exponential smoothing: the practitioner's ARIMA

**Simple exponential smoothing (SES):** a level $\ell_t=\alpha x_t+(1-\alpha)\ell_{t-1}$, $0<\alpha<1$, forecast $\hat x_{t+h}=\ell_t$ for all $h$. Unrolling, $\ell_t=\alpha\sum_{j\ge0}(1-\alpha)^jx_{t-j}$: a weighted average with geometrically decaying weights. **Holt's** method adds a slope $b_t$ updated the same way (forecast $\ell_t+hb_t$); **Holt–Winters** adds seasonal indices. They are popular because they are cheap, robust, and need no identification step. They are ARIMA models in disguise.

**Theorem 3.12.** SES with parameter $\alpha$ gives the same one-step forecasts as the ARIMA(0,1,1) $x_t=x_{t-1}+w_t+\theta w_{t-1}$ with $\theta=\alpha-1$ (and the same initialization).

▶ **Try it first.** For ARIMA(0,1,1), write $\hat x_{t+1}$ given data to $t$ in terms of $x_t$ and $w_t$, then replace $w_t$ by the one-step forecast error.

1. $x_{t+1}=x_t+w_{t+1}+\theta w_t$. The future shock $w_{t+1}$ is unpredictable, so $\hat x_{t+1}=x_t+\theta w_t$.
2. The shock is the one-step forecast error (invertibility, $|\theta|<1$): $w_t=x_t-\hat x_t$.
3. Substitute: $\hat x_{t+1}=x_t+\theta(x_t-\hat x_t)=(1+\theta)x_t-\theta\hat x_t$.
4. With $\alpha=1+\theta$: $\hat x_{t+1}=\alpha x_t+(1-\alpha)\hat x_t$, which is the SES recursion for $\ell_t=\hat x_{t+1}$. ∎
$$\text{SES}(\alpha)\equiv\text{ARIMA}(0,1,1)\text{ with }\theta=\alpha-1;\quad P_{n+h}=\sigma^2\big(1+(h-1)\alpha^2\big).\tag{3.37}$$

■ **What just happened.** SES is optimal exactly when the series is a random walk observed with noise (the ARIMA(0,1,1) is equivalent to that local-level model; ST-8 shows it via the Kalman filter). Invertibility $|\theta|<1$ allows $0<\alpha<2$; the usual $\alpha\in(0,1)$ is $\theta\in(-1,0)$. The ARIMA view adds what SES alone does not give you: a likelihood for $\alpha$ and the prediction variance in (3.37), from $\psi_0=1$, $\psi_j=1+\theta=\alpha$ for $j\ge1$ and (3.30). Stated without proof: Holt's linear method is an ARIMA(0,2,2), and additive Holt–Winters is a seasonal ARIMA (Hyndman & Athanasopoulos ch. 8–9).

### 3.14 What comes next

ST-8 Fourier-transforms $\gamma$: the spectral density $f(\nu)=\sum_h\gamma(h)e^{-2\pi i\nu h}$ ($\nu$ in cycles per sample; ST-8 (3.3)) turns convolution by $\psi$ into multiplication, a complex-root AR(2) into a peak near the frequency $\omega/2\pi$ of (3.15), and the sample ACF into the periodogram. The peak is not exactly at $\omega/2\pi$: damping pulls it to a slightly lower frequency (ST-8 (3.11); for Example 4.2 the spectral peak has period 8.69 against the ACF's quasi-period 8). It writes every ARMA in state-space form (the companion form (3.11) is one), so that the Kalman filter computes (3.24) exactly for any ARMA with missing data. And it drops the assumption that the shocks have constant variance: GARCH models the clustering of large shocks seen in financial returns (project "Your turn" 3).

## 4. Worked examples

All numbers are asserted or printed in Part 0 of `project.py`.

### Example 4.1: ACF of an AR(1) and an MA(1)

AR(1) with $\phi=0.7$, $\sigma^2=1$. By (3.8), $\gamma(0)=1/(1-0.49)=1/0.51=1.96078$, and $\rho(h)=0.7^h$:

| $h$ | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| $\rho(h)$ | 1 | 0.7 | 0.49 | 0.343 | 0.2401 | 0.16807 |

MA(1) with $\theta=0.5$, $\sigma^2=1$. By (3.9), $\gamma(0)=1.25$, $\gamma(1)=0.5$, $\gamma(h)=0$ for $h\ge2$, so $\rho(1)=0.5/1.25=0.4$ and $\rho(h)=0$ beyond. The MA(1) with $\theta=2$, $\sigma^2=0.25$ has the same autocovariances: $\gamma(0)=0.25\cdot5=1.25$, $\gamma(1)=0.25\cdot2=0.5$.

### Example 4.2: AR(2) with complex roots, $\phi_1=1.0$, $\phi_2=-0.5$

*Triangle.* $\phi_1+\phi_2=0.5<1$, $\phi_2-\phi_1=-1.5<1$, $|\phi_2|=0.5<1$: stationary. $\phi_1^2+4\phi_2=1-2=-1<0$: complex roots.

*Roots.* $\phi(z)=1-z+0.5z^2=0\Rightarrow z=\frac{1\pm\sqrt{1-2}}{1}=1\pm i$, $|z|=\sqrt2>1$.

*Companion matrix.* $A=\begin{pmatrix}1&-0.5\\1&0\end{pmatrix}$, $\det(\lambda I-A)=\lambda^2-\lambda+0.5$, $\lambda=0.5\pm0.5i$, $|\lambda|=1/\sqrt2=0.7071$. Check (3.12): $1/(1+i)=(1-i)/2=0.5-0.5i$. ✓

*Quasi-period.* $r=\sqrt{0.5}$, $\cos\omega=1/(2\sqrt{0.5})=0.7071$, $\omega=\pi/4$, period $2\pi/\omega=8$.

*ACF.* $\rho(1)=\phi_1/(1-\phi_2)=1/1.5=0.6667$, then $\rho(h)=\rho(h-1)-0.5\rho(h-2)$: $\rho(2)=0.6667-0.5=0.1667$, $\rho(3)=0.1667-0.3333=-0.1667$, $\rho(4)=-0.1667-0.0833=-0.25$, $\rho(5)=-0.1667$, $\rho(6)=-0.0417$, $\rho(7)=0.0417$, $\rho(8)=0.0625$. Solving for $a,b$ in (3.15): $\rho(h)=2^{-h/2}\big[\cos(\pi h/4)+\tfrac13\sin(\pi h/4)\big]$; at $h=8$, $2^{-4}=0.0625$. ✓ The variance by (3.16): $\gamma(0)=1.5/(0.5\cdot(2.25-1))=2.4$.

*Physics reading.* An underdamped oscillator $\ddot x+\gamma\dot x+\omega_0^2x=$ noise (CM-1 §3.8) sampled every $\Delta$ has homogeneous solutions $e^{(-\gamma/2\pm i\omega_1)t}$, i.e. eigenvalues $e^{-\gamma\Delta/2}e^{\pm i\omega_1\Delta}$ per step. Matching: $e^{-\gamma\Delta/2}=0.7071$ and $\omega_1\Delta=\pi/4$: an oscillator that loses half its energy ($r^2=0.5$) and advances one eighth of a cycle per sample.

### Example 4.3: Yule–Walker for AR(2)

Given $\gamma(0)=2$, $\gamma(1)=1.2$, $\gamma(2)=0.4$, solve $\Gamma_2\boldsymbol\phi=\boldsymbol\gamma_2$:
$$\begin{pmatrix}2&1.2\\1.2&2\end{pmatrix}\begin{pmatrix}\phi_1\\\phi_2\end{pmatrix}=\begin{pmatrix}1.2\\0.4\end{pmatrix}.$$
$\det=4-1.44=2.56$. $\phi_1=(2\cdot1.2-1.2\cdot0.4)/2.56=1.92/2.56=0.75$; $\phi_2=(2\cdot0.4-1.2\cdot1.2)/2.56=-0.64/2.56=-0.25$. Then $\sigma^2=\gamma(0)-\phi_1\gamma(1)-\phi_2\gamma(2)=2-0.9+0.1=1.2$. Check: $\phi_1+\phi_2=0.5$, $\phi_2-\phi_1=-1$, $|\phi_2|=0.25$, inside the triangle; $\phi_1^2+4\phi_2=-0.4375<0$, complex eigenvalues with modulus $\sqrt{0.25}=0.5$. The project confirms that this AR(2) reproduces $\gamma(0),\gamma(1),\gamma(2)$ exactly.

### Example 4.4: PACF of an AR(1) by the regression definition

$\phi=0.7$, so $\gamma\propto(1,0.7,0.49,\dots)$ (the constant cancels). Lag 1: $\phi_{11}=\rho(1)=0.7$. Lag 2: solve (3.17),
$$\begin{pmatrix}1&0.7\\0.7&1\end{pmatrix}\begin{pmatrix}\phi_{21}\\\phi_{22}\end{pmatrix}=\begin{pmatrix}0.7\\0.49\end{pmatrix}.$$
Try $(0.7,0)$: first row $0.7$, second row $0.49$. ✓ So $\phi_{22}=0$. In the partial-correlation form: the residuals after removing $x_{t-1}$ are $x_t-0.7x_{t-1}$ and $x_{t-2}-0.7x_{t-1}$ (the backward predictor has the same coefficient by symmetry). Their covariance, in units of $\gamma(0)$, is $\rho(2)-0.7\rho(1)-0.7\rho(1)+0.49\rho(0)=0.49-0.49-0.49+0.49=0$. Lag 3 is also 0 (project).

### Example 4.5: ARMA(1,1) with $\phi=0.5$, $\theta=0.4$, $\sigma^2=1$

$\psi$-weights by (3.21): $\psi_0=1$, $\psi_1=0.9$, $\psi_2=0.45$, $\psi_3=0.225$, $\psi_4=0.1125$.

$\gamma(0)=(1+2(0.5)(0.4)+0.16)/(1-0.25)=1.56/0.75=2.08$. $\gamma(1)=(0.9)(1.2)/0.75=1.44$. $\gamma(2)=0.72$, $\gamma(3)=0.36$.

$\rho(1)=1.44/2.08=1.08/1.56=0.69231$, $\rho(2)=0.34615$, $\rho(3)=0.17308$. Ratio $\rho(h+1)/\rho(h)=0.5=\phi$ from lag 1 on. An AR(1) with the same $\rho(1)=0.692$ would have $\rho(2)=0.479$; the ARMA's ACF starts higher (the MA term) and decays faster relative to $\rho(1)$.

### Example 4.6: Forecasting an AR(1), $\phi=0.8$, $\sigma^2=1$, $\mu=0$, $x_n=2$

By (3.29): one step, $\hat x_{n+1}=1.6$, $P=1$, 95% interval $1.6\pm1.96=[-0.36,\ 3.56]$.

Five steps: $\hat x_{n+5}=0.8^5\cdot2=0.32768\cdot2=0.65536$. $P_{n+5}=(1-0.8^{10})/(1-0.64)=(1-0.10737)/0.36=2.47952$, $\sqrt P=1.57465$, interval $0.65536\pm3.08631=[-2.4309,\ 3.7417]$. The limit as $h\to\infty$: forecast 0, variance $\gamma(0)=1/0.36=2.778$; after five steps the interval is already 94% as wide as it will ever get ($\sqrt{2.4795/2.7778}=0.945$).

### Example 4.7: The exact AR(1) likelihood for three observations

$\phi=0.6$, $\sigma^2=1$, $x=(0.5,\ 1.0,\ -0.2)$. By (3.25), term by term ($\log\varphi$ is the Normal log density):

- $x_1\sim\mathcal N(0,1/0.64=1.5625)$: $-\tfrac12\log(2\pi\cdot1.5625)-0.25/(2\cdot1.5625)=-1.14208-0.08=-1.22208$.
- $x_2\mid x_1\sim\mathcal N(0.3,1)$, residual $0.7$: $-0.91894-0.245=-1.16394$.
- $x_3\mid x_2\sim\mathcal N(0.6,1)$, residual $-0.8$: $-0.91894-0.32=-1.23894$.

Total $\log L=-3.62496$. The same number comes from the $3\times3$ MVN with Toeplitz covariance $\Gamma_3=\frac1{0.64}\begin{pmatrix}1&0.6&0.36\\0.6&1&0.6\\0.36&0.6&1\end{pmatrix}$ (project, by Cholesky). The CLS likelihood would drop the first term and condition on $x_1$.

### Example 4.8: SES is ARIMA(0,1,1)

Data $x=(10,12,11,13,12)$, $\alpha=0.3$, $\ell_1=x_1=10$. SES: $\ell_2=0.3(12)+0.7(10)=10.6$; $\ell_3=0.3(11)+0.7(10.6)=10.72$; $\ell_4=0.3(13)+0.7(10.72)=11.404$; $\ell_5=0.3(12)+0.7(11.404)=11.5828$.

ARIMA(0,1,1) with $\theta=-0.7$, starting with $\hat x_2=x_1=10$: $\hat x_3=x_2+\theta(x_2-\hat x_2)=12-0.7(2)=10.6$; $\hat x_4=11-0.7(11-10.6)=10.72$; $\hat x_5=13-0.7(13-10.72)=11.404$; $\hat x_6=12-0.7(12-11.404)=11.5828$. Identical, as Theorem 3.12 says. The forecast for every horizon is $11.5828$; the ARIMA view adds the variance $\sigma^2(1+0.09(h-1))$.

## 5. Why it matters and where it connects

- **ST-8.** The spectral density is the Fourier transform (F3) of (3.6); the companion form (3.11) is a state-space model and the Kalman filter is the innovations algorithm (3.24) for it; GARCH extends the white-noise assumption.
- **ST-3.** Every regression on time-indexed data needs this module: autocorrelated errors invalidate the standard errors (§3.12), and the cure is either modeling the errors as ARMA (regression with ARMA errors, fit by (3.24)), HAC standard errors, or differencing.
- **ST-6.** An MCMC chain is a stationary time series; the effective sample size (ST-6 (3.18)) is $n/\big(1+2\sum_{h\ge1}\rho(h)\big)$, and for an AR(1)-like chain with $\rho(h)=\phi^h$ it is $n(1-\phi)/(1+\phi)$ (ST-6 Example 4.5). Its $\hat R$ is a nonstationarity diagnostic.
- **LA-7, LA-10.** Stationarity is spectral radius < 1 of the companion matrix; the likelihood is a Toeplitz Cholesky; Durbin–Levinson is the $O(n^2)$ Toeplitz solver (Levinson 1947, for Wiener's filtering problem).
- **CM-1, SM.** AR(1) is the sampled Ornstein–Uhlenbeck process (Brownian motion in a harmonic trap, the Langevin dynamics of SM-7); AR(2) with complex roots is the sampled noisy underdamped oscillator (Example 4.2). Fluctuation–dissipation is a statement about $\gamma(h)$.
- **ML-8.** A causal 1-D convolution is a finite MA filter; a recurrent network is a nonlinear AR. Evaluating any learned forecaster requires the time-ordered splits and the ARIMA baseline from this module.

## 6. Exercises

See `exercises.md`: 6 warm-ups (ACFs, roots, random-walk correlation, redundancy, bands, Ljung–Box), 8 core problems (valid ACFs, the AR(2) ACF, Durbin–Levinson by hand, ARMA forecasting, MA estimation, YW vs CLS, spurious regression in code, ARIMA forecast variance), 3 challenges (the Dickey–Fuller identity and limit, the error-correction form, OU sampling). Essential: A2, B1, B3, B4, B6, B7.

## 7. Coding project

`project.py` builds the subject from scratch with numpy/scipy (statsmodels is not used): Part 0 reproduces every worked example; Part 1 draws the identification table for nine processes and checks sample ACFs against theory with Bartlett's formula; Part 2 plots the AR(2) triangle and asserts the root/eigenvalue reciprocity; Part 3 compares Yule–Walker, CLS and exact ML over 500 replications; Part 4 asserts PACF = Durbin–Levinson; Part 5 checks forecast coverage; Part 6 compares ARIMA and trend-regression intervals on a drifting random walk; Part 7 calibrates Ljung–Box; Part 8 simulates the Dickey–Fuller distribution and power; Part 9 does spurious regression and Engle–Granger; Part 10 checks SES = ARIMA(0,1,1) on a long series.

**Build it yourself first:** `sacf`, `durbin_levinson`, `ar_loglik_exact` and `df_t`. Write each from the equation, then compare with the script. The "Your turn" section asks for a seasonal ARIMA, Yule's sunspot model with its implied period, and a GARCH(1,1) whose returns are white noise while their squares are not.

## 8. Where you'll struggle

**Stationarity itself.** After a month of time series most students still cannot state the definition. It is two conditions: constant mean, and $\mathrm{Cov}(x_{t+h},x_t)$ depends only on $h$. Drill it until it is automatic, then use it as a test on every model you meet: write down $\mathbb Ex_t$ and $\mathrm{Cov}(x_{t+h},x_t)$ and see whether $t$ survives. A random walk fails because $t$ survives in $\mathrm{Var}(x_t)=t\sigma^2$; a trend fails because $t$ survives in the mean.

**Reading the ACF of a random walk as structure.** A slowly decaying sample ACF, significant at every lag, is the signature of a unit root, not of a rich AR model. If you fit an AR(5) to a random walk you get $\hat\phi_1\approx1$ and noise elsewhere. Rule: if the ACF has not dropped below about $0.5$ by lag 10, difference before identifying anything.

**AR stationarity as an eigenvalue condition.** "Roots outside the unit circle" is memorized and inverted constantly (inside or outside? of which polynomial?). Do not memorize it. Rederive it from the companion matrix in two lines: the recurrence is stable iff eigenvalues are inside the circle (LA-7), and eigenvalues are reciprocal roots (3.12). The sign convention $\phi(z)=1-\phi_1z-\dots$ also bites: with $x_t=1.5x_{t-1}-0.75x_{t-2}+w_t$, $\phi(z)=1-1.5z+0.75z^2$, not $1+1.5z-0.75z^2$.

**MA invertibility and the silent wrong $\theta$.** The Gaussian likelihood of an MA(1) is the same at $\theta$ and $1/\theta$ (with $\sigma^2$ rescaled). An optimizer can return $\hat\theta=2.1$, a CLS recursion $w_t=x_t-\theta w_{t-1}$ explodes for $|\theta|>1$, and two packages report different "answers" for the same data. Constrain $|\theta|<1$ or reflect the roots of $\hat\theta(z)$ outside the unit circle after fitting. Near $\theta=-1$ (an over-differenced series) everything is unstable.

**Not recognizing F5's formula.** The exact likelihood is $\mathcal N(\mathbf 0,\Gamma_n)$ with a Toeplitz covariance; that is all. Students see (3.25) or the innovations form and think it is new. Write the $3\times3$ case of Example 4.7 both ways once, and the connection stays.

**The Dickey–Fuller distribution.** This module only sketches why the limit is non-standard; the honest summary is: the regressor is a random walk, so sums of its squares scale like $n^2$, the estimator scales like $1/n$, and the limit is a ratio of Brownian functionals. Know the table (with constant: 5% value $-2.86$), know which regression you ran (constant? trend?), and know that the test has little power within 0.1 of a unit root.

**Spurious regression is the most important practical lesson.** It is not exotic: any two trending series (GDP and the number of mobile phones, temperature and stock prices) give a "significant" regression. Before regressing levels on levels, test each series for a unit root; if both have one, either difference or test for cointegration with the Engle–Granger critical values, not the DF ones.

**Box–Jenkins identification is an art.** Real sample ACFs and PACFs are noisy; with $n=100$ the bands are $\pm0.2$ wide and a true ARMA(1,1) can look like an AR(2) or an MA(3). Fit two or three candidates, compare AIC and residual diagnostics, and prefer the simpler model when forecasts barely differ. The module can only start this skill; Shumway & Stoffer ch. 3 has many real examples to practice on.

## 9. Check yourself

1. State weak stationarity. Is $x_t=w_t+w_{t-1}+t$ stationary? Is $\nabla x_t$?
2. Why must $\Gamma_n$ be PSD, and why does $\hat\gamma$ use divisor $n$?
3. What is the approximate distribution of $\hat\rho(5)$ for white noise with $n=400$?
4. Write $\gamma(h)$ for a linear process.
5. Derive $\gamma(h)$ for AR(1) in two lines.
6. Why are $\theta=0.5$ and $\theta=2$ indistinguishable for an MA(1), and which do you report?
7. Write the companion matrix of an AR(3) and the relation between its eigenvalues and the roots of $\phi$.
8. Is $x_t=0.9x_{t-1}+0.2x_{t-2}+w_t$ stationary?
9. Which of ACF/PACF cuts off for AR($p$), and why?
10. Give the $h$-step forecast and MSE for AR(1), and their limits.
11. Why does the DF $t$ statistic not have a $t$ distribution, and what is the 5% critical value with a constant?
12. Why does regressing one independent random walk on another give $|t|>2$ most of the time?

## 10. Pair with

- **Shumway & Stoffer, *Time Series Analysis and Its Applications*, ch. 1–3** (free "EZ" edition online). The applied core: real data sets (sunspots, global temperature, SOI, GNP) for every method here, and R code. Better than these notes for practice at identification.
- **Brockwell & Davis, *Introduction to Time Series and Forecasting*, ch. 1–6.** Cleaner theory: the Bartlett formula, Durbin–Levinson and the innovations algorithm derived in full (ch. 2–3), Wold's decomposition. Read it for the proofs this module states.
- **Hyndman & Athanasopoulos, *Forecasting: Principles and Practice*, ch. 8–9** (free online). Exponential smoothing in its state-space form and ARIMA from a forecaster's view, with time-series cross-validation. The practical complement.
- **Hamilton, *Time Series Analysis*, ch. 1–5 (and 17–19 for unit roots and cointegration).** The reference: difference equations via companion matrices (ch. 1–2, the same derivation as §3.5), exact likelihoods for every ARMA (ch. 5), and the Brownian-motion derivation of (3.34) (ch. 17).
