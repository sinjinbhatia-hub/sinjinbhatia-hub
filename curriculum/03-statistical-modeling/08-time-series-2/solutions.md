# ST-8 Solutions

Every number below is printed or asserted by `project.py` (Part 13 and the parts cited).

## A. Warm-up

**A1.** $\theta=-0.5$, $\lambda=2\pi\omega$: $f(\omega)=|1-0.5e^{-i\lambda}|^2=1-\cos\lambda+0.25=1.25-\cos(2\pi\omega)$.
- $\omega=0$: $1.25-1=0.25$.
- $\omega=\tfrac14$: $1.25-0=1.25$.
- $\omega=\tfrac12$: $1.25+1=2.25$.

The density rises with $\omega$, so the filter is high-pass. That fits: $x_t=w_t-0.5w_{t-1}$ has $\rho(1)=-0.5/1.25=-0.4$, so neighbouring values tend to have opposite signs. Integral: $\int_{-1/2}^{1/2}(1.25-\cos2\pi\omega)\,d\omega=1.25-0=1.25$, and $\gamma(0)=1+\theta^2=1.25$. ✓

**A2.** (a) The roots of $\phi(z)=1-0.5z-0.3z^2$ come from $0.3z^2+0.5z-1=0$: $z=\frac{-0.5\pm\sqrt{0.25+1.2}}{0.6}$, i.e. $z=1.1736$ and $z=-2.8403$. Both are real and lie outside the unit circle, so the model is stationary. (Triangle check from ST-7: $\phi_1+\phi_2=0.8<1$, $\phi_2-\phi_1=-0.2<1$, $|\phi_2|<1$.)

(b) $\cos2\pi\omega^*=\frac{0.5(0.3-1)}{4(0.3)}=\frac{-0.35}{1.2}=-0.2917$. This lies inside $(-1,1)$, but $\phi_2=0.3>0$, so $d^2D/dc^2=-8\phi_2<0$. The stationary point maximizes $D$, so it **minimizes** $f$. There is no interior peak.

(c) $f(0)=1/(1-0.5-0.3)^2=1/0.04=25.0$. $f(\tfrac12)=1/(1+0.5-0.3)^2=1/1.44=0.6944$. The stationary point is at $\omega=\arccos(-0.2917)/2\pi=0.2971$, where $f=0.4897$. Sketch: a tall spike at $\omega=0$ (the root at 1.17 is close to the unit circle, so the process is persistent), falling to a trough of 0.49 near $\omega=0.30$ and rising slightly to 0.69 at $\omega=\tfrac12$. The negative root $-2.84$ is the small high-frequency bump.

**A3.** $m=3$ gives $L=2m+1=7$ and $\nu=2L=14$. With $\chi^2_{14,0.025}=5.629$ and $\chi^2_{14,0.975}=26.119$, the multipliers are $14/26.119=0.536$ and $14/5.629=2.487$, so the interval is $[0.536\hat f,2.487\hat f]$. Its end-to-end ratio is $4.64$, against $6.31$ for $m=2$. The interval is narrower, and in exchange the estimate averages over a bandwidth $7/n$ instead of $5/n$, which costs bias at peaks.

**A4.** (a) Use Harvey's form with $r=\max(1,2)=2$. State $\alpha_t=(x_t,\ \theta w_t)^\top$.
$$T=\begin{pmatrix}0.6&1\\0&0\end{pmatrix},\quad R=\begin{pmatrix}1\\0.4\end{pmatrix},\quad Q=1,\quad Z=(1\ \ 0),\quad H=0,$$
with $\eta_t=w_{t+1}$. Check: $\alpha_{t+1}^{(2)}=0\cdot\alpha_t+0.4w_{t+1}=0.4w_{t+1}$, and $\alpha_{t+1}^{(1)}=0.6\alpha_t^{(1)}+\alpha_t^{(2)}+w_{t+1}=0.6x_t+0.4w_t+w_{t+1}=x_{t+1}$. ✓ (Project Part 9 uses exactly this form.)

(b) State $\alpha_t=(\mu_t,\nu_t)^\top$:
$$T=\begin{pmatrix}1&1\\0&1\end{pmatrix},\quad R=I_2,\quad Q=\begin{pmatrix}0.1&0\\0&0.01\end{pmatrix},\quad Z=(1\ \ 0),\quad H=2.$$

**A5.** With $H=1$, $P$ solves $P^2-0.25P-0.25=0$, so $P=\tfrac12(0.25+\sqrt{0.0625+1})=\tfrac12(0.25+1.0308)=0.6404$ and $K=P/(P+1)=0.6404/1.6404=0.3904$. Each observation gets 39% of the weight.
- $q=0.01$: $P=\tfrac12(0.01+\sqrt{0.0001+0.04})=0.1051$, so $K=0.0951$. The level barely moves, so the filter averages about the last $1/K\approx10$ observations.
- $q=100$: $P=\tfrac12(100+\sqrt{10000+400})=100.99$, so $K=0.9902$. The level jumps far more than the noise, so the filter believes almost only the latest point.

**A6.** $\bar\sigma^2=0.02/(1-0.98)=1.0$. Persistence $\alpha+\beta=0.98$. Half-life $\log0.5/\log0.98=34.31$ periods. By (3.36), $\sigma^2_{t+5|t}=1+0.98^4(2-1)=1+0.9224=1.9224$. After four steps of reversion only 8% of the excess has decayed, because persistence is so close to 1.

## B. Core

**B1.** (a) $\mathbb Ex_t=0$. Since $\mathbb EAB=0$,
$$\mathbb Ex_{t+h}x_t=s^2[\cos(2\pi\omega_0(t+h))\cos(2\pi\omega_0t)+\sin(2\pi\omega_0(t+h))\sin(2\pi\omega_0t)]=s^2\cos(2\pi\omega_0h),$$
by the cosine subtraction formula. It depends only on $h$, so the process is stationary, with $\gamma(h)=s^2\cos(2\pi\omega_0h)$.

(b) Write $(A,B)=(C\cos\psi,\ C\sin\psi)$ in polar coordinates. Then $A\cos a+B\sin a=C\cos(a-\psi)$, so $x_t=C\cos(2\pi\omega_0t+\Phi)$ with $\Phi=-\psi$. The density of $(A,B)$ is $\frac{1}{2\pi s^2}e^{-(A^2+B^2)/2s^2}$. It depends only on the radius, so in polar coordinates it factors as $\frac{1}{2\pi}\cdot\frac{c}{s^2}e^{-c^2/2s^2}$. Therefore $\psi$, and hence $\Phi$ mod $2\pi$, is uniform and independent of $C$, and $C^2/s^2\sim\chi^2_2$. Conditional on $C$, §3.1 gives $\gamma(h)=\tfrac{C^2}{2}\cos(2\pi\omega_0h)$. Averaging over $C$ with $\mathbb EC^2=2s^2$ recovers (a).

(c) Independence makes the autocovariances add: $\gamma(h)=s_1^2\cos(0.2\pi h)+s_2^2\cos(0.6\pi h)$. The spectral distribution $F$ is a step function with jumps of $s_1^2/2$ at $\omega=\pm0.1$ and $s_2^2/2$ at $\pm0.3$, total mass $s_1^2+s_2^2=\gamma(0)$. There is no density, because $\gamma$ does not decay ($\sum|\gamma(h)|=\infty$), so Theorem 3.2 does not apply, and a function cannot put positive mass at a single point. The periodogram of such a series has spikes that grow like $n$.

**B2.** (a) $\psi(z)=\frac1k\sum_{j=0}^{k-1}z^j=\frac1k\frac{1-z^k}{1-z}$ for $z\ne1$. At $z=e^{-i\lambda}$, factor out half-angles:
$$1-e^{-ik\lambda}=e^{-ik\lambda/2}(e^{ik\lambda/2}-e^{-ik\lambda/2})=2ie^{-ik\lambda/2}\sin(k\lambda/2),$$
and likewise $1-e^{-i\lambda}=2ie^{-i\lambda/2}\sin(\lambda/2)$. The ratio is $e^{-i(k-1)\lambda/2}\frac{\sin(k\lambda/2)}{\sin(\lambda/2)}$. With $\lambda=2\pi\omega$, $|\psi|^2=\big(\frac{\sin\pi k\omega}{k\sin\pi\omega}\big)^2$, which is (3.12). Zeros: $\sin(\pi k\omega)=0$ with $\sin\pi\omega\ne0$, i.e. $\omega=j/k$ for $j=1,\dots,\lfloor k/2\rfloor$.

(b) $k=4$. At $\omega=\tfrac18$: $\sin(\pi/2)/(4\sin(\pi/8))=1/(4\times0.38268)=0.6533$, so $|\psi|^2=0.4268$. At $\omega=\tfrac14$: $\sin\pi=0$, so $|\psi|^2=0$. At $\omega=\tfrac12$: $\sin2\pi=0$, so $|\psi|^2=0$. Both zeros of the 4-term MA lie at $\tfrac14$ and $\tfrac24$.

(c) A fixed seasonal pattern of period 12 that sums to zero over a year is a combination of the frequencies $j/12$ for $j=1,\dots,6$ (its Fourier series, F2). By (a), these are exactly the zeros of the 12-term MA, so the MA annihilates the pattern and all its harmonics. (Equivalently, the MA averages over one full year, and a zero-sum pattern averages to zero.) A 3-year cycle has $\omega=1/36$, where $|\psi|^2=\big(\frac{\sin(\pi/3)}{12\sin(\pi/36)}\big)^2=0.6857$. Its variance is cut by about 31% and its amplitude by 17%. Seasonal adjustment by moving averages distorts slow cycles too.

**B3.** (a) $\mathbb EI(\omega_j)=\frac1n\sum_{s,t}\mathbb E[x_sx_t]e^{-2\pi i\omega_j(t-s)}$. White noise has $\mathbb Ex_sx_t=\sigma^2$ if $s=t$ and 0 otherwise, so only the $n$ diagonal terms survive and each equals $\sigma^2e^0$. Hence $\mathbb EI=\frac1n\cdot n\sigma^2=\sigma^2$. Only uncorrelatedness was used, not Gaussianity or independence.

(b) At $\omega=\tfrac12$, $e^{-2\pi i t/2}=(-1)^t$ is real, so the sine sum is identically zero and $d(\tfrac12)=n^{-1/2}\sum_t(-1)^tx_t\sim\mathcal N(0,\sigma^2)$. Hence $I(\tfrac12)=d^2=\sigma^2\chi^2_1$, with mean $\sigma^2$ and variance $2\sigma^4$. The cosine/sine pair that gives two degrees of freedom in Theorem 3.5 collapses to one real component at $\omega=0$ and $\omega=\tfrac12$. These two ordinates are twice as variable.

(c) $I(\omega_j)$ is the squared modulus of a single DFT coefficient. More data makes that coefficient a sum of more terms, but it is normalized by $n^{-1/2}$, so its distribution stays $\mathcal N_{\mathbb C}$ with the same variance: the ordinate stays $\tfrac{f}{2}\chi^2_2$. What improves is resolution: the Fourier frequencies become $1/n$ apart, so there are more nearly independent ordinates per unit bandwidth to average, which is what smoothing exploits.

**B4.** (a) The filter's squared response is $|1+e^{-i\lambda}|^{2s}|1-e^{-i\lambda}|^{2d}$. Here $|1+e^{-i\lambda}|^2=2+2\cos\lambda=4\cos^2(\pi\omega)$, and from §3.4, $|1-e^{-i\lambda}|^2=4\sin^2(\pi\omega)$. So $\log f_y=\text{const}+2s\log\cos\pi\omega+2d\log\sin\pi\omega$. Differentiate: $\pi(-2s\tan\pi\omega+2d\cot\pi\omega)=0$, which gives $\tan^2\pi\omega^*=d/s$. Then $\cos2\pi\omega^*=\frac{1-\tan^2\pi\omega^*}{1+\tan^2\pi\omega^*}=\frac{s-d}{s+d}$. It is a maximum because $f_y$ vanishes at both ends ($\omega=0$ when $d\ge1$, $\omega=\tfrac12$ when $s\ge1$) and is positive in between. For $s=6$, $d=2$: $\cos2\pi\omega^*=0.5$, so $\omega^*=1/6$, period 6 (project Part 2: 6.00).

(b) $s=d$ gives $\cos2\pi\omega^*=0$, so $\omega^*=\tfrac14$, period 4. The low-pass and high-pass factors are mirror images ($\omega\leftrightarrow\tfrac12-\omega$), so their product peaks in the middle.

(c) $s=d=1$: $(1+B)(1-B)=1-B^2$, so $y_t=w_t-w_{t-2}$. Then $\gamma(0)=2$, $\gamma(1)=\mathrm{Cov}(w_{t+1}-w_{t-1},w_t-w_{t-2})=0$, $\gamma(2)=\mathrm{Cov}(w_{t+2}-w_t,w_t-w_{t-2})=-1$, and $\gamma(h)=0$ for $h\ge3$. So $\rho(1)=0$ and $\rho(2)=-\tfrac12$. A value tends to reverse sign two steps later and return four steps later: a period-4 oscillation. That is (b)'s peak at $\omega=\tfrac14$, seen in the time domain. Directly, $|1-e^{-2i\lambda}|^2=4\sin^2(2\pi\omega)$, which is maximal at $\omega=\tfrac14$.

**B5.** From Example 4.6, $a_{1|1}=(1.6923,1.4615)^\top$ and $P_{1|1}=\begin{pmatrix}0.6923&0.4615\\0.4615&1.3077\end{pmatrix}$.

*Predict.* $a_{2|1}=Ta_{1|1}=(1.6923+1.4615,\ 1.4615)^\top=(3.1538,1.4615)^\top$.
$TP_{1|1}=\begin{pmatrix}1.1538&1.7692\\0.4615&1.3077\end{pmatrix}$ (row 1 is row 1 plus row 2). Then $(TP_{1|1})T^\top$ adds column 2 to column 1: $\begin{pmatrix}2.9231&1.7692\\1.7692&1.3077\end{pmatrix}$. Adding $RQR^\top=\begin{pmatrix}0.25&0.5\\0.5&1\end{pmatrix}$ gives $P_{2|1}=\begin{pmatrix}3.1731&2.2692\\2.2692&2.3077\end{pmatrix}$.

*Update.* $v_2=2.5-3.1538=-0.6538$. $F_2=3.1731+1=4.1731$. $K_2=(3.1731,2.2692)^\top/4.1731=(0.7604,0.5438)^\top$.
$a_{2|2}=(3.1538-0.7604\times0.6538,\ 1.4615-0.5438\times0.6538)^\top=(2.6567,1.1060)^\top$.

Why velocity went down: the filter predicted the cart at 3.15, but the sensor saw 2.5, so the cart is behind where the velocity estimate put it. Position and velocity errors are positively correlated in $P_{2|1}$ (2.27), so a negative position surprise implies the cart is slower than believed, and velocity drops by $0.5438\times0.6538=0.36$. Velocity is never measured; it is inferred from how positions change. The project checks that running both observations through `kalman_filter` in one call gives the same $a_{2|2}$.

**B6.** (a) $t=1$: $v_1=y_1$, $F_1=P_0+H$, $K_1=P_0/(P_0+H)$, $a_{1|1}=K_1y_1$, $P_{1|1}=P_0H/(P_0+H)$. Predict: $a_{2|1}=K_1y_1$, $P_{2|1}=P_0H/(P_0+H)+Q$. $t=2$: $v_2=y_2-K_1y_1$, $F_2=P_0H/(P_0+H)+Q+H$.

(b) $y_1=\mu_1+\varepsilon_1$ and $y_2=\mu_1+\eta_1+\varepsilon_2$, so $\Sigma=\begin{pmatrix}P_0+H&P_0\\P_0&P_0+Q+H\end{pmatrix}$.

(c) Determinant: $F_1F_2=(P_0+H)\big(\frac{P_0H}{P_0+H}+Q+H\big)=P_0H+(Q+H)(P_0+H)=P_0Q+2P_0H+QH+H^2$. Also $\det\Sigma=(P_0+H)(P_0+Q+H)-P_0^2=P_0Q+2P_0H+QH+H^2$. Equal. ✓

Quadratic form: $v=Ly$ with $L=\begin{pmatrix}1&0\\-K_1&1\end{pmatrix}$. Compute $\mathrm{Cov}(v)=L\Sigma L^\top$:
- $\mathrm{Var}\,v_1=P_0+H=F_1$.
- $\mathrm{Cov}(v_1,v_2)=\mathrm{Cov}(y_1,y_2)-K_1\mathrm{Var}\,y_1=P_0-\frac{P_0}{P_0+H}(P_0+H)=0$.
- $\mathrm{Var}\,v_2=(P_0+Q+H)-2K_1P_0+K_1^2(P_0+H)=P_0+Q+H-\frac{P_0^2}{P_0+H}=\frac{P_0H}{P_0+H}+Q+H=F_2$.

So $L\Sigma L^\top=D=\mathrm{diag}(F_1,F_2)$, hence $\Sigma^{-1}=L^\top D^{-1}L$ and $y^\top\Sigma^{-1}y=(Ly)^\top D^{-1}(Ly)=v_1^2/F_1+v_2^2/F_2$. Combined with the determinant, $-\tfrac12[2\log2\pi+\log F_1+\log F_2+v_1^2/F_1+v_2^2/F_2]=\log\mathcal N(y;0,\Sigma)$. This is (3.27) for $n=2$, and it shows what the filter is doing: $L$ is the inverse of the unit lower-triangular Cholesky factor of $\Sigma$ (LA-10), computed one row at a time. Project Part 9 checks the same identity for an ARMA(1,1) with $n=200$.

**B7.** (a) $y_t-y_{t-1}=\mu_t-\mu_{t-1}+\varepsilon_t-\varepsilon_{t-1}=\eta_{t-1}+\varepsilon_t-\varepsilon_{t-1}$. This is a sum of white noises at lags 0 and 1 only, so $\gamma(h)=0$ for $|h|\ge2$: an MA(1). $\gamma(0)=Q+2H$ and $\gamma(1)=\mathrm{Cov}(\eta_t+\varepsilon_{t+1}-\varepsilon_t,\ \eta_{t-1}+\varepsilon_t-\varepsilon_{t-1})=-H$. So $\rho(1)=-H/(Q+2H)=-1/(q+2)$.

(b) An MA(1) with coefficient $\theta$ has $\rho(1)=\theta/(1+\theta^2)$ (ST-7). Let $p=P/H$ be the scaled steady-state predicted variance, so that $p^2=qp+q$ by (3.24), and $K=p/(p+1)$. Set $\theta=K-1=-1/(p+1)$. Then
$$\frac{\theta}{1+\theta^2}=\frac{-1/(p+1)}{1+1/(p+1)^2}=\frac{-(p+1)}{(p+1)^2+1}.$$
We need this to equal $-1/(q+2)$, i.e. $(p+1)^2+1=(q+2)(p+1)$. Expanding, $p^2+2p+2=qp+q+2p+2$, i.e. $p^2=qp+q$, which is the Riccati equation. ✓ Since $p>0$, $|\theta|=1/(p+1)<1$: this is the invertible root.

(c) The ARIMA(0,1,1) $y_t=y_{t-1}+e_t+\theta e_{t-1}$ has one-step forecast $\hat y_{t+1}=y_t+\theta e_t$ with $e_t=y_t-\hat y_t$, so $\hat y_{t+1}=\hat y_t+(1+\theta)(y_t-\hat y_t)$. That is simple exponential smoothing with constant $1+\theta=K$. The steady-state local level filter's forecast is $a_{t+1|t}=a_{t|t}=a_{t|t-1}+K(y_t-a_{t|t-1})$, the same recursion. For $q=1$: $K=0.6180$, $\theta=-0.3820$, and $\theta/(1+\theta^2)=-0.3333=-1/3$. ✓ Project Part 7 also checks numerically that the filter and SES coincide.

**B8.** (a) As in §3.15, $\mu_2=\omega/(1-\alpha)$. Now $\mu_4=\mathbb E\sigma_t^4z_t^4=\kappa_z\mathbb E\sigma_t^4$, because $z_t$ is independent of $\sigma_t$. Then $\mathbb E\sigma_t^4=\mathbb E(\omega+\alpha\varepsilon_{t-1}^2)^2=\omega^2+2\omega\alpha\mu_2+\alpha^2\mu_4$. So $\mu_4=\kappa_z(\omega^2+2\omega\alpha\mu_2)+\kappa_z\alpha^2\mu_4$, giving
$$\mu_4(1-\kappa_z\alpha^2)=\kappa_z\omega^2\Big(1+\frac{2\alpha}{1-\alpha}\Big)=\kappa_z\omega^2\frac{1+\alpha}{1-\alpha}.$$
Divide by $\mu_2^2=\omega^2/(1-\alpha)^2$: $\kappa=\kappa_z\frac{1-\alpha^2}{1-\kappa_z\alpha^2}$, finite iff $\kappa_z\alpha^2<1$. Since $1-\alpha^2>1-\kappa_z\alpha^2$ when $\kappa_z>1$, we get $\kappa>\kappa_z$.

(b) $\kappa_z=3$, $\alpha=0.5$: $\kappa=3\times0.75/0.25=9.0$. For $\alpha\ge1/\sqrt3=0.5774$ the fourth moment is infinite, although the variance $\omega/(1-\alpha)$ is still finite for $\alpha<1$. Sample kurtosis then never settles and grows with $n$.

(c) Conditionally on the past, each return is Normal, but with a random variance $\sigma_t^2$. A mixture of Normals with different variances has more mass near zero (low-variance periods) and in the tails (high-variance periods) than a single Normal with the same overall variance. Algebraically, this is $\mathbb E\sigma^4>(\mathbb E\sigma^2)^2$ (Jensen).

**B9.** (a) From (3.34), $\varepsilon_t^2-\bar\sigma^2$ is an ARMA(1,1) with $\phi=\alpha+\beta$ and $\theta=-\beta$, driven by the white noise $u_t$ (assuming $\mathbb E\varepsilon^4<\infty$). Substitute:
- $1+\phi\theta=1-(\alpha+\beta)\beta=1-\alpha\beta-\beta^2$;
- $\phi+\theta=\alpha$;
- $1+2\phi\theta+\theta^2=1-2\alpha\beta-2\beta^2+\beta^2=1-2\alpha\beta-\beta^2$.

So $\rho(1)=\frac{\alpha(1-\alpha\beta-\beta^2)}{1-2\alpha\beta-\beta^2}$. For $h\ge2$ the MA part has dropped out, so $\rho(h)=(\alpha+\beta)\rho(h-1)$ (ST-7's ARMA(1,1) ACF).

(b) $\alpha=0.1$, $\beta=0.85$: numerator $0.1(1-0.085-0.7225)=0.1\times0.1925=0.01925$, denominator $1-0.17-0.7225=0.1075$, so $\rho(1)=0.1791$. Then $\rho(2)=0.95\times0.1791=0.1701$ and $\rho(10)=0.95^9\times0.1791=0.1129$.

(c) No. Sample autocorrelations of $\varepsilon_t^2$ are very noisy because their variance involves $\mathbb E\varepsilon^8$. Averaged over 40 independent paths of 4000, the project gets 0.191, close to 0.179, while single paths scatter widely. The ACF of $\varepsilon^2$ exists iff $\mathbb E\varepsilon^4<\infty$, i.e. $(\alpha+\beta)^2+2\alpha^2<1$; here $0.9225$. For the sample ACF to be $\sqrt n$-consistent and asymptotically Normal you need $\mathbb E\varepsilon^8<\infty$, which for GARCH(1,1) means $\mathbb E(\alpha z^2+\beta)^4<1$. With $\mathbb Ez^4=3$, $\mathbb Ez^6=15$, $\mathbb Ez^8=105$, this is $105\alpha^4+60\alpha^3\beta+18\alpha^2\beta^2+4\alpha\beta^3+\beta^4=0.9592$, which holds but only barely. The estimator therefore converges slowly. For real returns with $\alpha+\beta\approx0.99$ the condition usually fails, and sample ACFs of squared returns should be read as qualitative.

## C. Challenge

**C1.** (a) A diffuse prior means $P_{1|0}^{-1}=0$. In information form (§3.10), $P_{1|1}^{-1}=0+Z^\top H^{-1}Z=1/1+1/4=1.25$, so the posterior variance is $0.80$. The mean is $a_{1|1}=P_{1|1}Z^\top H^{-1}y=0.8(y_1/1+y_2/4)=0.80y_1+0.20y_2$. The weights are in the ratio $4:1$ of the precisions, and the fused variance $0.8$ is below the better sensor's 1.

(b) Start from $P_{t|t}=P^--P^-Z^\top F^{-1}ZP^-$ with $P^-=P_{t|t-1}$ and $F=ZP^-Z^\top+H$:
$$P_{t|t}Z^\top H^{-1}=P^-Z^\top H^{-1}-P^-Z^\top F^{-1}ZP^-Z^\top H^{-1}=P^-Z^\top F^{-1}\big[F-ZP^-Z^\top\big]H^{-1}=P^-Z^\top F^{-1}HH^{-1}=K_t.$$
(We wrote $P^-Z^\top H^{-1}=P^-Z^\top F^{-1}FH^{-1}$ to factor.) If every row of $Z$ equals the same row vector $z$, then $Z^\top=z^\top(1,\dots,1)$ and $K_t=P_{t|t}z^\top(1,\dots,1)H^{-1}$. Column $j$ is $P_{t|t}z^\top/\sigma_j^2$. All columns are parallel, scaled by sensor precision, at every $t$, not only in steady state.

(c) Add to Part 13: `Z2 = [[1,0],[1,0]]`, `H = diag(1, 16)`, run `kalman_filter` with the tracker's $T$, $R$, $Q=1$ on any two-column data. The project prints the steady-state gain $\begin{pmatrix}0.7106&0.0444\\0.4802&0.0300\end{pmatrix}$ and the column ratio $16.0000$ (asserted). The noisy sensor still improves the estimate, but only by 1/16 as much per reading.

**C2.** (a) Stationary variance: if $\mathrm{Var}\,\alpha_t=s^2$ then $\mathrm{Var}\,\alpha_{t+1}=\phi^2s^2+s^2(1-\phi^2)=s^2$. With $\alpha_1\sim\mathcal N(0,s^2)$, every $\alpha_t$ has variance $s^2$. Covariance: iterating, $\alpha_{t+h}=\phi^h\alpha_t+(\text{noise after }t)$, so $\mathrm{Cov}(\alpha_{t+h},\alpha_t)=\phi^hs^2=s^2e^{-h/\ell}=k(h)$. The AR(1) is zero-mean Gaussian with the GP's covariance on the integers, so the two have the same joint law.

(b) Both $K(K+\sigma_n^2I)^{-1}y$ (F5's conditional mean for the joint Gaussian of $(f,y)$) and the smoother's $a_{t|n}=\mathbb E[\alpha_t\mid Y_n]$ (Theorem 3.8) are the conditional expectation of the latent process given the same data under the same joint distribution, so they are equal. Cost: the kernel formula solves an $n\times n$ dense system, $O(n^3)$ time and $O(n^2)$ memory. The filter and smoother do $O(1)$ work per time step, $O(n)$ total. The predictive variance $P_{t|n}$ equals the GP posterior variance by the same argument.

(c) The project builds both and prints a maximum difference of $1.4\times10^{-15}$ (asserted below $10^{-10}$).

**C3.** (a) The project's `bootstrap_pf` works as follows. Draw $N$ particles from the initial distribution. At each $t$: compute log-weights $\log p(y_t\mid x^{(i)})$; add $\log\frac1N\sum_iw^{(i)}$ to the running log-likelihood (subtract the maximum log-weight first for numerical stability, then add it back); resample indices with probabilities proportional to $w$; propagate each particle through the transition. On a local level model with $Q=H=1$, $n=100$ and $N=2000$, five runs give $-190.51$ (sd 0.44) against the Kalman value $-190.35$.

(b) Before resampling, the particles are (approximately) draws from $p(\alpha_t\mid Y_{t-1})$. Hence $\frac1N\sum_ip(y_t\mid\alpha_t^{(i)})\approx\int p(y_t\mid\alpha_t)p(\alpha_t\mid Y_{t-1})d\alpha_t=p(y_t\mid Y_{t-1})$: importance sampling with the prior as proposal (ST-6). Multiplying over $t$ gives $p(Y_n)$ by the chain rule. The product of these averages is an **unbiased** estimate of the likelihood (a theorem of Del Moral; the resampling steps do not break it). Its log is biased **downward** by Jensen's inequality, by roughly half the variance of the log-estimate. That matches the project's PF mean sitting slightly below the Kalman value.

(c) Code sketch, to add under "Your turn":
```python
mu, phi_true, s2eta, n = -1.0, 0.95, 0.05, 1000
h = np.empty(n); h[0] = mu + np.sqrt(s2eta / (1 - phi_true**2)) * rng.standard_normal()
for t in range(1, n):
    h[t] = mu + phi_true * (h[t-1] - mu) + np.sqrt(s2eta) * rng.standard_normal()
y = np.exp(h / 2) * rng.standard_normal(n)
def sv_ll(phi, N=2000):
    r = np.random.default_rng(0)                      # common random numbers across phi
    x0 = mu + np.sqrt(s2eta / (1 - phi**2)) * r.standard_normal(N)
    prop = lambda x, rr: mu + phi * (x - mu) + np.sqrt(s2eta) * rr.standard_normal(x.size)
    obs = lambda yt, x: -0.5 * (np.log(2 * np.pi) + x + yt**2 * np.exp(-x))
    return bootstrap_pf(y, prop, obs, x0, N, r)
grid = np.linspace(0.80, 0.99, 20)
prof = [sv_ll(p) for p in grid]
```
Here `obs` is $\log\mathcal N(y_t;0,e^{x})$. Fixing the seed inside `sv_ll` (common random numbers) makes the profile a smooth-looking curve instead of Monte Carlo noise. Expect the maximum near the true 0.95 and a profile that is flat over a few hundredths: with $n=1000$, $\phi$ in a stochastic volatility model is weakly identified, because volatility is never observed (§8). For proper inference you would use particle MCMC (ST-6's Metropolis–Hastings with the PF likelihood plugged in) rather than a grid.
