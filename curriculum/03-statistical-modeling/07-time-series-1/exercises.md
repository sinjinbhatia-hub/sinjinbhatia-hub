# ST-7 Exercises

Conventions as in the notes: $\phi(z)=1-\phi_1z-\dots-\phi_pz^p$, $\theta(z)=1+\theta_1z+\dots+\theta_qz^q$, $w_t\sim\mathrm{WN}(0,\sigma^2)$, $\sigma^2=1$ unless stated, mean zero unless stated. Equation numbers refer to `notes.md`. Numerical answers are printed by Part 11 of `project.py`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (MA(2) ACF.) For $x_t=w_t+0.6w_{t-1}-0.3w_{t-2}$, compute $\gamma(0),\gamma(1),\gamma(2),\gamma(3)$ and $\rho(1),\rho(2)$. Is the process invertible? (Find the roots of $\theta$.)

**A2.** (Stationary or not.) For each AR(2), find the roots of $\phi(z)$ and the companion eigenvalues, decide whether it is causal and stationary, and say whether its ACF oscillates. If it oscillates, give the quasi-period.
(a) $x_t=1.2x_{t-1}-0.2x_{t-2}+w_t$. (b) $x_t=0.4x_{t-1}+0.45x_{t-2}+w_t$. (c) $x_t=1.5x_{t-1}-0.75x_{t-2}+w_t$.
For (a), factor $\phi(z)$ and say what model $\nabla x_t$ follows.

**A3.** (Random walk.) For $x_t=\sum_{s=1}^tw_s$, derive $\mathrm{Cov}(x_s,x_t)$ for $s\le t$ and $\mathrm{corr}(x_{100},x_{110})$. Explain in two sentences why the sample ACF of a random walk stays near 1 for many lags even though "$\rho(h)$" is not defined for it.

**A4.** (Parameter redundancy.) (a) Show that $x_t=0.5x_{t-1}+w_t-0.5w_{t-1}$ is white noise. (b) Simplify $x_t=0.3x_{t-1}+0.4x_{t-2}+w_t+0.5w_{t-1}$ by factoring both polynomials. What model is it really?

**A5.** (Bands.) A residual series of length $n=400$ has $\hat\rho(1)=0.08$ and $\hat\rho(2)=0.12$. (a) Which lies outside the $\pm2/\sqrt n$ band? (b) You inspect 20 lags of a true white-noise series. What is the expected number of lags outside the band, and the probability that at least one is outside (treat the lags as independent, by (3.5))?

**A6.** (Ljung–Box by hand.) Residuals of length $n=100$ have $\hat\rho(1)=0.15$, $\hat\rho(2)=-0.10$, $\hat\rho(3)=0.05$. Compute $Q$ with $K=3$ and its $p$-value (a) if the residuals are claimed to be white noise, (b) if they come from a fitted AR(1).

## B. Core (do all)

**B1.** (Valid autocorrelations.) Let $\rho(0)=1$, $\rho(\pm1)=r$, $\rho(h)=0$ for $|h|\ge2$.
(a) For $r=0.6$, show that the $5\times5$ Toeplitz matrix $R_5=[\rho(i-j)]$ is not PSD, so no stationary process has this ACF. (Hint: the eigenvalues of the $n\times n$ tridiagonal Toeplitz matrix with diagonal 1 and off-diagonal $r$ are $1+2r\cos(k\pi/(n+1))$, $k=1,\dots,n$.) Show that $R_4$ is PD.
(b) Prove that this $\rho$ is a valid ACF iff $|r|\le\frac12$.

**B2.** (AR(2) ACF and variance.) (a) Derive the variance formula (3.16) from the Yule–Walker equations (3.13) for $h=1,2$ and the $\sigma^2$ equation. (b) For $x_t=1.5x_{t-1}-0.75x_{t-2}+w_t$, compute $\gamma(0)$ and $\rho(1),\dots,\rho(4)$, and the modulus $r$ and frequency $\omega$ in (3.15).

**B3.** (Durbin–Levinson by hand.) A series has $\hat\rho(1)=0.8$, $\hat\rho(2)=0.5$, $\hat\rho(3)=0.3$. (a) Use (3.19) to compute $\hat\phi_{11},\hat\phi_{22},\hat\phi_{33}$ and the AR(3) Yule–Walker coefficients. (b) Compute the prediction error variances $v_1,v_2,v_3$ (in units of $\hat\gamma(0)$). (c) With $n=100$, which AR order does the PACF suggest?

**B4.** (Forecasting an AR(2).) For $x_t=x_{t-1}-0.5x_{t-2}+w_t$ (Example 4.2) with $x_{n-1}=2$, $x_n=1$: compute $\hat x_{n+1},\hat x_{n+2},\hat x_{n+3}$, the $\psi$-weights $\psi_0,\psi_1,\psi_2$, the MSEs (3.30), and 95% intervals. What do the forecasts and the MSE converge to as $h\to\infty$?

**B5.** (Estimating an MA(1) by moments.) (a) Show that the moment equation $\hat\rho(1)=\theta/(1+\theta^2)$ has a real solution iff $|\hat\rho(1)|\le\frac12$, and that its two solutions are reciprocals. (b) Find the invertible solution for $\hat\rho(1)=0.4$. (c) What should an estimation routine do if $\hat\rho(1)=0.55$? (d) Explain why the Gaussian likelihood itself cannot break the tie between $\theta$ and $1/\theta$.

**B6.** (Yule–Walker vs CLS for AR(1).) With known mean 0, $\hat\phi_{\rm YW}=\sum_{t=2}^nx_tx_{t-1}/\sum_{t=1}^nx_t^2$ and $\hat\phi_{\rm CLS}=\sum_{t=2}^nx_tx_{t-1}/\sum_{t=2}^nx_{t-1}^2$. (a) Prove $|\hat\phi_{\rm YW}|<1$ for any data not identically zero. (b) Compute both for $x=(1,2,4)$. (c) Show $\hat\phi_{\rm CLS}/\hat\phi_{\rm YW}=1+x_n^2/\sum_{t=1}^{n-1}x_t^2$ and explain why the two agree for a long stationary series but not for a short or explosive one.

**B7.** (Spurious regression, in code.) Extend Part 9 of `project.py`. (a) Report the fraction of $|t|>2$ for independent random walks at $n=50,100,200,400$. (b) Regress $\nabla y$ on $\nabla x$ instead, at $n=400$, and report the fraction. (c) Use the heuristic of Theorem 3.11 to predict how the fraction behaves as $n$ grows, and check it against (a).

**B8.** (ARIMA forecasts.) $y_t$ is ARIMA(1,1,0): $\nabla y_t=0.5\nabla y_{t-1}+w_t$. With $y_{n-1}=9$, $y_n=10$: (a) compute $\hat y_{n+1},\hat y_{n+2},\hat y_{n+3}$ and their limit; (b) find the $\psi$-weights of $y_t$ (from $(1-0.5z)(1-z)\psi(z)=1$) and the MSEs for $h=1,2,3$; (c) show $P_{n+h}$ grows linearly in $h$ for large $h$ and find the slope.

## C. Challenge (do at least one)

**C1.** (Superconsistency.) Let $x_t=x_{t-1}+w_t$, $x_0=0$, and $\hat\phi=\sum_{t=1}^nx_{t-1}x_t/\sum_{t=1}^nx_{t-1}^2$. (a) Prove $\sum_{t=1}^nx_{t-1}w_t=\frac12(x_n^2-\sum_{t=1}^nw_t^2)$. (b) Using (a) and $\mathbb E\sum x_{t-1}^2$, explain why $n(\hat\phi-1)$ has a nondegenerate limit and $\sqrt n(\hat\phi-1)\to0$. (c) Why is the limit distribution of $n(\hat\phi-1)$ skewed to the left? (d) Using Part 8 of the project, report the mean and 5% quantile of $n(\hat\phi-1)$ at $n=100$.

**C2.** (Error correction.) Let $x_t=x_{t-1}+v_t$ and $y_t=\beta x_t+u_t$ with $u_t=\phi u_{t-1}+\varepsilon_t$, $|\phi|<1$, $v,\varepsilon$ white noise. (a) Show both $x$ and $y$ are I(1) and $y-\beta x$ is stationary. (b) Derive the error-correction form (3.36) and identify $\alpha$ and $\gamma$. (c) In Part 9 the Engle–Granger statistic was compared with a simulated critical value. What is the size of the test if you use the ordinary DF 5% critical value instead, and why is it too large?

**C3.** (Ornstein–Uhlenbeck sampled.) The OU process $dX=-\kappa X\,dt+s\,dW$ has solution $X_{t+\Delta}=e^{-\kappa\Delta}X_t+s\int_0^\Delta e^{-\kappa(\Delta-u)}dW_{t+u}$. (a) Show that $X$ sampled every $\Delta$ is an AR(1) and find $\phi$ and $\sigma^2$ (use the Itô isometry $\mathbb E(\int f\,dW)^2=\int f^2du$). (b) Check that $\sigma^2/(1-\phi^2)$ equals the OU stationary variance $s^2/(2\kappa)$. (c) Evaluate for $\kappa=2$, $\Delta=0.1$, $s=1$. (d) What happens to $\phi$ as $\Delta\to0$, and what does that say about fitting an AR(1) to very finely sampled data and testing for a unit root?
