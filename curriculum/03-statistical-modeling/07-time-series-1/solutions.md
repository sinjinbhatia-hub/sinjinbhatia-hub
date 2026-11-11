# ST-7 Solutions

Numbers are printed by Part 11 of `project.py` unless another part is named. $\sigma^2=1$ unless stated.

## A. Warm-up

**A1.** $\psi_0=1$, $\psi_1=0.6$, $\psi_2=-0.3$. By (3.10):
- $\gamma(0)=1+0.36+0.09=1.45$;
- $\gamma(1)=\psi_0\psi_1+\psi_1\psi_2=0.6-0.18=0.42$;
- $\gamma(2)=\psi_0\psi_2=-0.3$;
- $\gamma(3)=0$ (no overlap at lag 3 for an MA(2)).

$\rho(1)=0.42/1.45=0.2897$, $\rho(2)=-0.3/1.45=-0.2069$. Invertibility: $\theta(z)=1+0.6z-0.3z^2=0\iff0.3z^2-0.6z-1=0\iff z=(0.6\pm\sqrt{0.36+1.2})/0.6=(0.6\pm1.2490)/0.6$, so $z=3.0817$ or $z=-1.0817$. Both have $|z|>1$: invertible.

**A2.** Roots of $\phi(z)=1-\phi_1z-\phi_2z^2$; eigenvalues are their reciprocals (3.12), the roots of $\lambda^2-\phi_1\lambda-\phi_2$.

(a) $\lambda^2-1.2\lambda+0.2=(\lambda-1)(\lambda-0.2)$: $\lambda=1,0.2$, roots $z=1,5$. An eigenvalue on the unit circle: not stationary (the triangle's edge $\phi_1+\phi_2=1$). Factoring, $\phi(z)=(1-z)(1-0.2z)$, so $(1-0.2B)(1-B)x_t=w_t$: $\nabla x_t$ is a stationary AR(1) with coefficient $0.2$, and $x_t$ is ARIMA(1,1,0). Its sample ACF would decay slowly.

(b) $\lambda^2-0.4\lambda-0.45=(\lambda-0.9)(\lambda+0.5)$: $\lambda=0.9,-0.5$, roots $z=1.1111$ and $z=-2$. Both $|\lambda|<1$: stationary. Real eigenvalues, so no oscillation in the usual sense; the ACF is $c_1(0.9)^h+c_2(-0.5)^h$, dominated by the slowly decaying positive term, with a small alternating component at short lags. (Triangle check: $0.85<1$, $0.05<1$, $0.45<1$.)

(c) $\lambda^2-1.5\lambda+0.75=0$: $\lambda=0.75\pm i\sqrt{0.75-0.5625}=0.75\pm0.4330i$, $|\lambda|=\sqrt{0.75}=0.8660$; roots $|z|=1/0.8660=1.1547$. Stationary, complex: the ACF oscillates. By (3.15), $\cos\omega=1.5/(2\cdot0.8660)=0.8660$, $\omega=\pi/6$, quasi-period $12$ (asserted in the project). A monthly series with this AR(2) would look roughly annual.

**A3.** For $s\le t$, $x_t=x_s+\sum_{u=s+1}^tw_u$ and the second sum is independent of $x_s$, so $\mathrm{Cov}(x_s,x_t)=\mathrm{Var}(x_s)=s\sigma^2$. Then $\mathrm{corr}(x_{100},x_{110})=100/\sqrt{100\cdot110}=\sqrt{100/110}=0.9535$.

The sample ACF compares the series with a shifted copy of itself. For a path that wanders slowly, $x_{t+h}$ and $x_t$ are close for all $h\ll n$ (correlation $\sqrt{t/(t+h)}$), so the sample products are large and positive and $\hat\rho(h)$ stays near 1, decaying roughly like $1-h/n$ times a constant. No population $\rho(h)$ exists because $\mathrm{Cov}(x_t,x_{t+h})$ depends on $t$; $\hat\rho$ is still computable, and what it measures is the slow wandering.

**A4.** (a) In backshift form, $(1-0.5B)x_t=(1-0.5B)w_t$. The causal solution is $x_t=(1-0.5B)^{-1}(1-0.5B)w_t=w_t$. Check with (3.20): $\psi_1=\theta_1+\phi_1\psi_0=-0.5+0.5=0$, and then $\psi_j=0.5\psi_{j-1}=0$ for all $j\ge1$. White noise.
(b) $\phi(z)=1-0.3z-0.4z^2=(1-0.8z)(1+0.5z)$ (expand: $1+0.5z-0.8z-0.4z^2$). $\theta(z)=1+0.5z$. The common factor $(1+0.5z)$ cancels: $(1-0.8B)x_t=w_t$, an AR(1) with $\phi=0.8$. Fitting the ARMA(2,1) would give large, strongly correlated standard errors for the three coefficients.

**A5.** (a) $2/\sqrt{400}=0.1$. $\hat\rho(1)=0.08$ is inside; $\hat\rho(2)=0.12$ is outside.
(b) Each lag is outside with probability 0.05, so the expected count is $20\times0.05=1$, and $P(\text{at least one})=1-0.95^{20}=0.6415$. One lag at 0.12 among twenty is not evidence of misfit; that is what the Ljung–Box test is for.

**A6.** $Q=100\cdot102\Big(\frac{0.0225}{99}+\frac{0.01}{98}+\frac{0.0025}{97}\Big)=10200(2.2727+1.0204+0.2577)\times10^{-4}=3.6219$.
(a) White-noise null: compare with $\chi^2_3$; the 95% point is 7.815, $p=0.3053$. Do not reject.
(b) After an AR(1) fit, $K-p-q=2$ degrees of freedom: 95% point 5.991, $p=0.1635$. Still no evidence of misfit. The smaller degrees of freedom make the test stricter, because fitting has already made the residual autocorrelations smaller than they would be for raw noise.

## B. Core

**B1.** (a) $R_n$ is tridiagonal Toeplitz with diagonal 1 and off-diagonals $r=0.6$. Its eigenvalues are $1+1.2\cos(k\pi/(n+1))$, and the smallest is at $k=n$: $1+1.2\cos(n\pi/(n+1))=1-1.2\cos(\pi/(n+1))$.
- $n=4$: $1-1.2\cos36^\circ=1-1.2(0.80902)=0.0292>0$, so $R_4$ is PD.
- $n=5$: $1-1.2\cos30^\circ=1-1.2(0.86603)=-0.0392<0$.

So $R_5$ has a negative eigenvalue, and for its eigenvector $\mathbf a$ we would get $\mathrm{Var}(\sum a_ix_i)=\gamma(0)\mathbf a^\top R_5\mathbf a<0$. By Theorem 3.2(d) no stationary process has this ACF. (The project prints $-0.0812$ at $n=6$ and $-0.1514$ at $n=10$; the minimum tends to $1-2r=-0.2$.)

(b) *Necessity.* The smallest eigenvalue of $R_n$ is $1-2|r|\cos(\pi/(n+1))$ (for $r<0$ the smallest is at $k=1$, same value). As $n\to\infty$ this tends to $1-2|r|$, and it is decreasing in $n$. If $|r|>\frac12$, it is negative for large enough $n$, so $\rho$ is invalid. *Sufficiency.* If $|r|\le\frac12$, solve $\theta/(1+\theta^2)=r$ (possible by B5(a)); the MA(1) with that $\theta$ has exactly this ACF by (3.9). So the ACF is valid iff $|r|\le\frac12$.

**B2.** (a) The Yule–Walker equations (3.13) at $h=1,2$, divided by $\gamma(0)$:
1. $\rho(1)=\phi_1+\phi_2\rho(1)\Rightarrow\rho(1)=\phi_1/(1-\phi_2)$.
2. $\rho(2)=\phi_1\rho(1)+\phi_2=\phi_1^2/(1-\phi_2)+\phi_2$.
3. The $\sigma^2$ equation: $\sigma^2=\gamma(0)\big[1-\phi_1\rho(1)-\phi_2\rho(2)\big]$.
4. Substitute: $1-\frac{\phi_1^2}{1-\phi_2}-\phi_2\Big(\frac{\phi_1^2}{1-\phi_2}+\phi_2\Big)=1-\phi_2^2-\frac{\phi_1^2(1+\phi_2)}{1-\phi_2}=\frac{(1+\phi_2)\big[(1-\phi_2)^2-\phi_1^2\big]}{1-\phi_2}$, using $1-\phi_2^2=(1+\phi_2)(1-\phi_2)$.
5. Hence $\gamma(0)=\dfrac{(1-\phi_2)\sigma^2}{(1+\phi_2)[(1-\phi_2)^2-\phi_1^2]}$, which is (3.16). Inside the triangle each factor is positive: $1+\phi_2>0$, and $(1-\phi_2)^2-\phi_1^2=(1-\phi_2-\phi_1)(1-\phi_2+\phi_1)>0$. On the edges $\gamma(0)\to\infty$.

(b) $\phi_1=1.5$, $\phi_2=-0.75$: $\gamma(0)=1.75/(0.25\cdot(3.0625-2.25))=1.75/0.203125=8.6154$. $\rho(1)=1.5/1.75=0.8571$; $\rho(2)=1.5(0.8571)-0.75=0.5357$; $\rho(3)=1.5(0.5357)-0.75(0.8571)=0.1607$; $\rho(4)=1.5(0.1607)-0.75(0.5357)=-0.1607$. Modulus $r=\sqrt{0.75}=0.8660$, $\omega=\arccos(0.8660)=\pi/6=0.5236$, period 12. The variance is large ($8.6$ times $\sigma^2$) because the eigenvalues are close to the circle.

**B3.** (a) By (3.19):
- $\hat\phi_{11}=\hat\rho(1)=0.8$.
- $\hat\phi_{22}=\dfrac{\hat\rho(2)-\hat\phi_{11}\hat\rho(1)}{1-\hat\phi_{11}\hat\rho(1)}=\dfrac{0.5-0.64}{1-0.64}=\dfrac{-0.14}{0.36}=-0.3889$. Update $\hat\phi_{21}=\hat\phi_{11}-\hat\phi_{22}\hat\phi_{11}=0.8(1+0.3889)=1.1111$.
- $\hat\phi_{33}=\dfrac{\hat\rho(3)-\hat\phi_{21}\hat\rho(2)-\hat\phi_{22}\hat\rho(1)}{1-\hat\phi_{21}\hat\rho(1)-\hat\phi_{22}\hat\rho(2)}=\dfrac{0.3-0.5556+0.3111}{1-0.8889+0.1944}=\dfrac{0.0556}{0.3056}=0.1818$.
- AR(3) coefficients: $\hat\phi_{31}=\hat\phi_{21}-\hat\phi_{33}\hat\phi_{22}=1.1111+0.1818(0.3889)=1.1818$, $\hat\phi_{32}=\hat\phi_{22}-\hat\phi_{33}\hat\phi_{21}=-0.3889-0.1818(1.1111)=-0.5909$, $\hat\phi_{33}=0.1818$. The project checks these against a direct $3\times3$ Toeplitz solve.

(b) $v_h=v_{h-1}(1-\hat\phi_{hh}^2)$ with $v_0=1$: $v_1=1-0.64=0.36$, $v_2=0.36(1-0.1512)=0.3056$, $v_3=0.3056(1-0.0331)=0.2955$. Adding lag 2 reduces the error variance by 15%; lag 3 by about 3%.
(c) Bands $\pm2/\sqrt{100}=\pm0.2$. $\hat\phi_{22}=-0.39$ is outside, $\hat\phi_{33}=0.18$ is inside (barely). The PACF suggests AR(2); fit AR(2) and AR(3) and compare AIC and residual diagnostics before deciding.

**B4.** Forecasts iterate the recursion with future shocks set to 0:
- $\hat x_{n+1}=x_n-0.5x_{n-1}=1-1=0$;
- $\hat x_{n+2}=\hat x_{n+1}-0.5x_n=0-0.5=-0.5$;
- $\hat x_{n+3}=\hat x_{n+2}-0.5\hat x_{n+1}=-0.5-0=-0.5$.

$\psi$-weights by (3.20): $\psi_0=1$, $\psi_1=\phi_1=1$, $\psi_2=\phi_1\psi_1+\phi_2\psi_0=1-0.5=0.5$. MSEs by (3.30): $P_{n+1}=1$, $P_{n+2}=1+1=2$, $P_{n+3}=2+0.25=2.25$. Intervals $\hat x\pm1.96\sqrt P$: $[-1.96,1.96]$, $[-3.2719,2.2719]$, $[-3.44,2.44]$. As $h\to\infty$ the forecasts go to the mean 0 as a damped cosine with period 8 (the forecast function satisfies the same recurrence as the ACF), and $P_{n+h}\to\gamma(0)=2.4$ (Example 4.2).

**B5.** (a) $\theta/(1+\theta^2)=r$ is $r\theta^2-\theta+r=0$. For $r\ne0$ the discriminant is $1-4r^2$, nonnegative iff $|r|\le\frac12$. The product of the roots is $r/r=1$, so they are reciprocals: $\theta$ and $1/\theta$.
(b) $\theta=\dfrac{1\pm\sqrt{1-4(0.16)}}{2(0.4)}=\dfrac{1\pm0.6}{0.8}$, i.e. $2$ or $0.5$. The invertible solution is $\hat\theta=0.5$.
(c) No MA(1) has $|\rho(1)|>\frac12$ (B1). Either the sample ACF is noisy (with $n=50$, the sd of $\hat\rho(1)$ is about $0.1$, so $0.55$ is compatible with $\rho(1)=0.45$), or the model is wrong (an AR or a higher-order MA). A routine should report the boundary value $\hat\theta=1$ with a warning, or better, switch to maximum likelihood. Moment estimators for MA models are also inefficient; their asymptotic variance exceeds the MLE's.
(d) The Gaussian likelihood depends on the parameters only through $\Gamma_n$, i.e. through $\gamma(0)$ and $\gamma(1)$ (F5 Theorem 3.21 with the Toeplitz covariance). By §3.4, $(\theta,\sigma^2)$ and $(1/\theta,\theta^2\sigma^2)$ give the same $\gamma$, hence the same $\Gamma_n$ and the same likelihood at every data set. The data carry no information to break the tie; only the invertibility convention does.

**B6.** (a) By Cauchy–Schwarz, $\big|\sum_{t=2}^nx_tx_{t-1}\big|\le\Big(\sum_{t=2}^nx_t^2\Big)^{1/2}\Big(\sum_{t=1}^{n-1}x_t^2\Big)^{1/2}$. Each factor is at most $S=\sum_{t=1}^nx_t^2$, so the bound is $\le S$. For equality we need both $\sum_{t=2}^nx_t^2=S$ and $\sum_{t=1}^{n-1}x_t^2=S$, i.e. $x_1=0$ and $x_n=0$; and equality in Cauchy–Schwarz needs $(x_2,\dots,x_n)\propto(x_1,\dots,x_{n-1})$. With $x_1=0$ the proportionality forces $x_2=c\,x_1=0$, then $x_3=0$, and so on: the data are identically zero. Otherwise the inequality is strict and $|\hat\phi_{\rm YW}|<1$. This is the $p=1$ case of the general fact that Yule–Walker fits are causal.
(b) $x=(1,2,4)$: $\sum x_tx_{t-1}=2+8=10$. $\hat\phi_{\rm CLS}=10/(1+4)=2$; $\hat\phi_{\rm YW}=10/(1+4+16)=10/21=0.4762$.
(c) $\hat\phi_{\rm CLS}/\hat\phi_{\rm YW}=\sum_{t=1}^nx_t^2/\sum_{t=1}^{n-1}x_t^2=1+x_n^2/\sum_{t=1}^{n-1}x_t^2$. For a stationary series, $x_n^2=O_p(1)$ while the sum is $\approx(n-1)\gamma(0)$, so the ratio is $1+O_p(1/n)$ and the two estimators agree to first order (both satisfy (3.26)). For a short series, or an explosive one where $x_n^2$ is comparable to all earlier squares combined (as in (b), $16$ vs $5$), the last observation dominates and they differ a lot. YW then shrinks toward zero, which is the bias seen in project Part 3 (mean $\hat\phi_1=0.978$ against CLS's $0.988$).

**B7.** (a) The project (4000 replications per $n$): $0.657$, $0.752$, $0.830$, $0.876$ at $n=50,100,200,400$. (Part 9's $10^4$ replications at $n=50$ give $0.669$; the difference is Monte Carlo error, about $\pm0.007$ per run.)
(b) Differencing makes both series iid noise, independent of each other, so the regression is the textbook case: the fraction is $0.049$, close to the nominal $0.05$ for $|t|>2$ (which is $P(|t_{398}|>2)=0.046$).
(c) By Theorem 3.11, $t=\sqrt n\,Z_n$ where $Z_n$ converges to a nondegenerate random variable $Z$ (a ratio of Brownian functionals). So $P(|t|>2)=P(|Z_n|>2/\sqrt n)\to P(Z\ne0)=1$, with the gap $P(|Z|\le2/\sqrt n)$ shrinking roughly like $1/\sqrt n$ if $Z$ has a density at 0. Check: the fraction not significant is $0.343,0.248,0.170,0.124$; the ratios for each doubling of $n$ are $0.72,0.69,0.73$, close to $1/\sqrt2=0.71$. The fraction never stops growing, so a larger sample makes the spurious regression worse, not better.

**B8.** Let $d_t=\nabla y_t$, an AR(1) with $\phi=0.5$; $d_n=y_n-y_{n-1}=1$.
(a) $\hat d_{n+h}=0.5^h$, so $\hat y_{n+1}=10+0.5=10.5$, $\hat y_{n+2}=10.5+0.25=10.75$, $\hat y_{n+3}=10.75+0.125=10.875$. The limit is $10+\sum_{h\ge1}0.5^h=11$: the forecast levels off because the expected future differences decay, but it does not revert to any fixed mean.
(b) $\psi(z)=1/[(1-0.5z)(1-z)]=\sum_k0.5^kz^k\cdot\sum_jz^j$, so $\psi_j=\sum_{k=0}^j0.5^k=2(1-0.5^{j+1})$: $\psi_0=1$, $\psi_1=1.5$, $\psi_2=1.75$. By (3.30): $P_{n+1}=1$, $P_{n+2}=1+2.25=3.25$, $P_{n+3}=3.25+3.0625=6.3125$.
(c) $\psi_j\to2$, so each extra step adds $\psi_{h-1}^2\to4$: $P_{n+h}\approx4h+\text{const}$ for large $h$. The slope $4\sigma^2=\sigma^2/(1-\phi)^2$ is the long-run variance of the differences, $\sum_h\gamma_d(h)$: a unit-root process accumulates the long-run variance of its increments once per step.

## C. Challenge

**C1.** (a) Square $x_t=x_{t-1}+w_t$: $x_t^2=x_{t-1}^2+2x_{t-1}w_t+w_t^2$. Sum over $t=1,\dots,n$: the $x^2$ terms telescope to $x_n^2-x_0^2=x_n^2$, so $x_n^2=2\sum x_{t-1}w_t+\sum w_t^2$, which rearranges to the claim. (The project checks the identity on $10^4$ paths.)
(b) $\hat\phi-1=\sum x_{t-1}w_t/\sum x_{t-1}^2$. By (a), the numerator is $\frac12(x_n^2-\sum w_t^2)$, and $x_n^2\sim n\sigma^2\chi^2_1$ (Gaussian case) while $\sum w_t^2=n\sigma^2+O_p(\sqrt n)$: numerator $=\frac{n\sigma^2}2(\chi^2_1-1)+O_p(\sqrt n)$, of order $n$. The denominator has mean $\sigma^2\sum_{t=1}^n(t-1)=\sigma^2n(n-1)/2$ and is of order $n^2$ (it is a sum of $n$ terms of size $t$; it is not concentrated around its mean, since it converges after dividing by $n^2$ to $\sigma^2\int_0^1W^2$, a random variable). So $\hat\phi-1=O_p(n)/O_p(n^2)=O_p(1/n)$: $n(\hat\phi-1)$ is $O_p(1)$ with a nondegenerate limit (3.34), and $\sqrt n(\hat\phi-1)=O_p(n^{-1/2})\to0$. The usual $\sqrt n$ asymptotics of (3.26) would give the degenerate limit $\mathcal N(0,1-\phi^2)=\mathcal N(0,0)$ at $\phi=1$, consistent with this but useless for testing.
(c) The numerator's limit $\frac12(W(1)^2-1)$ is negative whenever $W(1)^2<1$, which has probability $P(\chi^2_1<1)=0.6827$ (printed in Part 8), and it is bounded below by $-\frac12$ but unbounded above. The denominator is positive and is small exactly on paths that stay near 0, which are paths where $W(1)^2$ is small. Small negative numerators divided by small denominators produce the long left tail. So $\hat\phi$ is usually below 1, and the distribution of $n(\hat\phi-1)$ has its mass on the negative side.
(d) Part 8 prints mean $-1.772$ and 5% quantile $-7.833$ for $n(\hat\phi-1)$ at $n=100$ (no constant). Fuller's tables give about $-8$ for the 5% point: a test that rejects when $n(\hat\phi-1)<-8$ is an alternative to the $\tau$ test.

**C2.** (a) $x_t$ is a random walk: I(1). $y_t-\beta x_t=u_t$ is a stationary AR(1) by construction. $\nabla y_t=\beta v_t+u_t-u_{t-1}$ is stationary (a sum of stationary processes), and $y_t$ itself is not, since $\mathrm{Var}(y_t)=\beta^2\mathrm{Var}(x_t)+\gamma_u(0)+2\beta\,\mathrm{Cov}(x_t,u_t)$ grows like $t\beta^2\sigma_v^2$ (for $\beta\ne0$, taking $v$ and $\varepsilon$ uncorrelated; the covariance term is bounded by Cauchy–Schwarz by $O(\sqrt t)$). So $y$ is I(1) and $(1,-\beta)$ is a cointegrating vector.
(b) Difference $y_t=\beta x_t+u_t$: $\nabla y_t=\beta\nabla x_t+\nabla u_t$. Now $\nabla u_t=u_t-u_{t-1}=(\phi-1)u_{t-1}+\varepsilon_t$, and $u_{t-1}=y_{t-1}-\beta x_{t-1}$. So
$$\nabla y_t=(\phi-1)(y_{t-1}-\beta x_{t-1})+\beta\nabla x_t+\varepsilon_t,$$
which is (3.36) with $\alpha=\phi-1<0$ and $\gamma=\beta$. With the project's $\phi=0.5$, half of last period's deviation is corrected each period. $|\alpha|$ is the speed of adjustment; $\alpha=0$ would mean $u$ is a random walk and no cointegration.
(c) Part 9 prints a size of $0.149$ when the Engle–Granger statistic is compared with the DF 5% value $-2.89$ instead of the simulated $-3.39$. The first-stage OLS chooses $\hat\beta$ to minimize $\sum\hat u_t^2$, so even when $y$ and $x$ are independent random walks it finds the combination that looks most stationary. The residuals therefore look more mean-reverting than a fixed combination would, the DF statistic is pushed left, and the DF critical value rejects too often. The correct critical values (MacKinnon) depend on the number of variables in the cointegrating regression.

**C3.** (a) Sampling at $t=k\Delta$: $X_{(k+1)\Delta}=e^{-\kappa\Delta}X_{k\Delta}+\eta_k$ with $\eta_k=s\int_0^\Delta e^{-\kappa(\Delta-u)}dW_{k\Delta+u}$. The $\eta_k$ are Gaussian, mean zero, and independent across $k$ (they use Brownian increments on disjoint intervals), and $\eta_k$ is independent of $X_{k\Delta}$. That is an AR(1) with $\phi=e^{-\kappa\Delta}$ and, by the Itô isometry,
$$\sigma^2=s^2\int_0^\Delta e^{-2\kappa(\Delta-u)}du=s^2\frac{1-e^{-2\kappa\Delta}}{2\kappa}.$$
(b) $\sigma^2/(1-\phi^2)=s^2\frac{1-e^{-2\kappa\Delta}}{2\kappa}\Big/(1-e^{-2\kappa\Delta})=s^2/(2\kappa)$, the OU stationary variance, independent of $\Delta$ as it must be.
(c) $\kappa=2$, $\Delta=0.1$, $s=1$: $\phi=e^{-0.2}=0.81873$, $\sigma^2=(1-e^{-0.4})/4=0.08242$. Check: $0.08242/(1-0.67032)=0.25=1/(2\cdot2)$.
(d) $\phi=e^{-\kappa\Delta}=1-\kappa\Delta+O(\Delta^2)\to1$. Sampling a mean-reverting process very finely makes it look like a random walk at the sampling scale: over $n$ samples you see a time span $n\Delta$, and if $n\Delta\ll1/\kappa$ the series has not had time to revert. A Dickey–Fuller test then has power close to its size, since the local alternative $\phi=1-\kappa\Delta$ is within $O(1/n)$ of 1 when $n\Delta$ is fixed. Power for unit-root tests comes from the length of the time span, not the number of observations; daily data over two years are less informative about mean reversion than annual data over a century. This is the same low-power problem as in §3.11, seen from the continuous-time side.
