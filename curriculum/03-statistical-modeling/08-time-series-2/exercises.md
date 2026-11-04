# ST-8 Exercises

Conventions as in the notes: $\omega$ in cycles per sample, $f(\omega)=\sum_h\gamma(h)e^{-2\pi i\omega h}$, $I(\omega_j)=|\sum_tx_te^{-2\pi i\omega_jt}|^2/n$; state space $\alpha_{t+1}=T\alpha_t+R\eta_t$, $y_t=Z\alpha_t+\varepsilon_t$; GARCH(1,1) $\sigma_t^2=\omega+\alpha\varepsilon_{t-1}^2+\beta\sigma_{t-1}^2$. White-noise variance $\sigma^2=1$ unless stated.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Compute the spectral density of the MA(1) $x_t=w_t-0.5w_{t-1}$ at $\omega=0,\tfrac14,\tfrac12$. Is it low-pass or high-pass? Check that $\int_{-1/2}^{1/2}f=\gamma(0)$ by integrating the formula.

**A2.** For the AR(2) with $\phi_1=0.5$, $\phi_2=0.3$: (a) check stationarity and say whether the roots are real or complex; (b) evaluate (3.11). Does the spectrum have an interior peak? (c) Compute $f(0)$, $f(\tfrac12)$ and $f$ at the stationary point you found, and sketch the spectrum.

**A3.** A Daniell smoother with $m=3$. What is $L$, what are the equivalent degrees of freedom, and what are the 95% interval multipliers in (3.16)? Compare their ratio with the $m=2$ case of Example 4.4.

**A4.** Write in state-space form, identifying $T,Z,R,Q,H$ and the state vector: (a) ARMA(1,1) with $\phi=0.6$, $\theta=0.4$, $\sigma^2=1$, and verify that your form reproduces $x_{t+1}=0.6x_t+w_{t+1}+0.4w_t$; (b) a local linear trend observed with noise of variance 2, level-noise variance 0.1, slope-noise variance 0.01.

**A5.** Local level model with $q=Q/H=0.25$ and $H=1$. Compute the steady-state predicted variance $P$ and gain $K$. Then compute $K$ for $q=0.01$ and $q=100$ and explain each in one sentence.

**A6.** GARCH(1,1) with $\omega=0.02$, $\alpha=0.08$, $\beta=0.9$. Compute $\bar\sigma^2$, the persistence, the half-life of a volatility shock, and the 5-step variance forecast $\sigma^2_{t+5|t}$ from $\sigma^2_{t+1}=2$.

## B. Core (do all)

**B1.** Let $x_t=A\cos(2\pi\omega_0t)+B\sin(2\pi\omega_0t)$ with $A,B$ independent $\mathcal N(0,s^2)$.
(a) Show $x_t$ is stationary and find $\gamma(h)$.
(b) Show that $x_t$ can be written as $C\cos(2\pi\omega_0t+\Phi)$ with $\Phi$ uniform on $(0,2\pi)$ independent of $C$, and identify the distribution of $C^2$.
(c) For the sum of independent such components at $\omega=0.1$ (variance $s_1^2$) and $\omega=0.3$ (variance $s_2^2$), write $\gamma(h)$ and describe the spectral distribution $F$. Why does $f$ not exist?

**B2.** (a) Derive (3.12) for the $k$-point moving average and locate its zeros in $(0,\tfrac12]$.
(b) For $k=4$, compute $|\psi|^2$ at $\omega=\tfrac18,\tfrac14,\tfrac12$.
(c) Explain, using (b)'s pattern, why a 12-term moving average of monthly data removes a fixed seasonal pattern of period 12 and all its harmonics, but also attenuates a 3-year cycle.

**B3.** (a) Let $x_t$ be white noise with variance $\sigma^2$, not necessarily Gaussian. Prove $\mathbb EI(\omega_j)=\sigma^2$ for $0<\omega_j<\tfrac12$.
(b) For $n$ even and Gaussian white noise, find the distribution of $I(\tfrac12)$, and its variance. Why does it differ from Theorem 3.5?
(c) In two sentences, explain why taking $n$ ten times larger does not shrink $\mathrm{Var}\,I(\omega)$ at a fixed frequency, and what does improve.

**B4.** (a) Derive the peak frequency of the filter $(1+B)^s(1-B)^d$ applied to white noise.
(b) Where is the peak when $s=d$? Interpret.
(c) Slutsky's point is that the output looks cyclical. Compute the output's lag-1 autocorrelation for $s=1$, $d=1$ (the filter $1-B^2$) directly from its coefficients, and reconcile it with (b).

**B5.** Continue Example 4.6 for one more step. The sensor reads $y_2=2.5$. Compute $a_{2|1}$, $P_{2|1}$, $F_2$, $K_2$, $v_2$ and $a_{2|2}$. Why did the velocity estimate go **down**?

**B6.** Local level model with $a_{1|0}=0$, $P_{1|0}=P_0$, and two observations $y_1,y_2$.
(a) Run the filter symbolically to get $v_1,F_1,v_2,F_2$.
(b) Write the joint covariance $\Sigma$ of $(y_1,y_2)$ directly from the model.
(c) Show that $F_1F_2=\det\Sigma$ and that $v_1^2/F_1+v_2^2/F_2=y^\top\Sigma^{-1}y$, so that (3.27) equals the bivariate Normal log density.

**B7.** Local level model with $q=Q/H$.
(a) Show that $\nabla y_t=y_t-y_{t-1}$ is an MA(1) and that its lag-1 autocorrelation is $-1/(q+2)$.
(b) Let $K$ be the steady-state gain. Prove that $\theta=K-1$ is the invertible MA(1) coefficient with that autocorrelation.
(c) Conclude that simple exponential smoothing with constant $K$ is the optimal forecast for ARIMA(0,1,1) with $\theta=K-1$ (ST-7's equivalence), and check it numerically for $q=1$.

**B8.** ARCH(1) with $\mathbb Ez_t^4=\kappa_z$ (not necessarily 3).
(a) Derive the kurtosis of $\varepsilon_t$ and the condition for it to be finite.
(b) Evaluate for Gaussian $z$ and $\alpha=0.5$. What happens for $\alpha\ge1/\sqrt3$?
(c) Explain in two sentences why ARCH returns have heavier tails than their Gaussian shocks.

**B9.** (a) Using (3.34) and ST-7's ARMA(1,1) autocorrelation $\rho(1)=\frac{(1+\phi\theta)(\phi+\theta)}{1+2\phi\theta+\theta^2}$, derive $\rho(1)$ for $\varepsilon_t^2$ in GARCH(1,1) and the recursion for $\rho(h)$.
(b) Evaluate $\rho(1)$, $\rho(2)$, $\rho(10)$ for $\alpha=0.1$, $\beta=0.85$.
(c) The project's single simulated path of 6000 gives a lag-1 sample autocorrelation of 0.145. Is this evidence against the formula? What finite-moment condition does a sample ACF of squares need to behave well?

## C. Challenge (do at least one)

**C1. (Sensor fusion; modifies `project.py`.)** (a) A constant scalar $x$ with prior variance $\infty$ is measured once by each of two sensors with error variances $\sigma_1^2=1$ and $\sigma_2^2=4$. Stack them as $y=Zx+\varepsilon$, $Z=(1,1)^\top$, $H=\mathrm{diag}(\sigma_1^2,\sigma_2^2)$, and derive the posterior mean and variance from the Kalman update. (b) Prove the identity $K_t=P_{t|t}Z^\top H^{-1}$ and use it to show that, at every step of any filter with diagonal $H$ and identical rows of $Z$, the gain columns are in the ratio of the sensor precisions. (c) In `project.py`, track the constant-velocity object with two position sensors, $H=\mathrm{diag}(1,16)$, and print the steady-state gain.

**C2. (GP regression is Kalman smoothing; code.)** Let $f$ be a zero-mean Gaussian process on the integers with $k(\tau)=s^2e^{-|\tau|/\ell}$, observed as $y_t=f(t)+\varepsilon_t$, $\varepsilon_t\sim\mathcal N(0,\sigma_n^2)$, $t=1,\dots,n$.
(a) Show that the AR(1) $\alpha_{t+1}=\phi\alpha_t+\eta_t$ with $\phi=e^{-1/\ell}$, $\eta_t\sim\mathcal N(0,s^2(1-\phi^2))$ and $\alpha_1\sim\mathcal N(0,s^2)$ has exactly the covariance $k$.
(b) Conclude that the GP posterior mean $K(K+\sigma_n^2I)^{-1}y$ equals the RTS smoother's $a_{t|n}$, and compare the costs.
(c) Verify numerically for $n=50$, $\ell=3$, $s^2=1$, $\sigma_n^2=0.25$.

**C3. (A particle filter; code.)** (a) Write a bootstrap particle filter and test it on a local level model with $Q=H=1$, where the exact log-likelihood is known from the Kalman filter. (b) Explain why $\sum_t\log\big(\frac1N\sum_iw_t^{(i)}\big)$ estimates the log-likelihood and whether it is unbiased. (c) Apply it to the stochastic volatility model $h_{t+1}=\mu+\phi(h_t-\mu)+\eta_t$, $y_t=e^{h_t/2}z_t$ (simulated with $\mu=-1$, $\phi=0.95$, $\mathrm{Var}\,\eta=0.05$) and profile the likelihood over $\phi$.
