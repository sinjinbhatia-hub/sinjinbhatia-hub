# ST-2 Solutions

Quantiles are upper-tail throughout, as in the notes. Every number below is printed or asserted in Part 10 (or the reference cells) of `project.py`.

## A. Warm-up

**A1.** (a) The standard error is $s/\sqrt n = 4.8/4 = 1.2$ and $t_{15,0.025} = 2.1314$, so the half-width is $2.1314\times1.2 = 2.558$ and the interval is $[49.742,\ 54.858]$.

(b) $z_{0.025}\sigma/\sqrt n = 1.960\times5/4 = 2.450$, giving $[49.850,\ 54.750]$. It is narrower because the multiplier matters more than the small difference between 5 and 4.8: $2.131\times4.8 = 10.23$ against $1.960\times5 = 9.80$. The $t$ interval pays for estimating σ with heavier-tailed quantiles.

(c) $(n-1)s^2 = 15\times23.04 = 345.6$. With $\chi^2_{15,0.025} = 27.488$ and $\chi^2_{15,0.975} = 6.262$, (3.4) gives $\sigma^2\in[345.6/27.488,\ 345.6/6.262] = [12.573,\ 55.19]$, so $\sigma\in[3.546,\ 7.429]$. The large χ² quantile gives the lower endpoint.

**A2.** *Wald:* $\hat p = 0$, so the estimated standard error is $\sqrt{0\cdot1/10} = 0$ and the interval is the single point $[0,0]$. It claims 95% confidence that $p$ is exactly 0. That is absurd: if $p = 0.1$, zero successes in 10 trials has probability $0.9^{10} = 0.349$.

*Wilson:* with $z^2 = 3.8415$, $z^2/n = 0.38415$. The centre is $(0+z^2/(2n))/(1+z^2/n) = 0.19208/1.38415 = 0.13877$. The half-width is $z\sqrt{0+z^2/(4n^2)}/(1+z^2/n) = (z^2/(2n))/(1+z^2/n) = 0.13877$. So the interval is $[0,\ 0.2775]$; its upper end is $z^2/(n+z^2) = 3.8415/13.8415 = 0.2775$. The rule of three, $3/n = 0.30$, comes from solving $(1-p)^n = 0.05$ approximately, a one-sided 95% bound; Wilson's two-sided upper end is a little tighter, and both say the same thing: $p$ up to roughly 0.3 is compatible with no successes in 10.

**A3.** (a) Size $= P_{1/2}(X\ge8) = \big[\binom{10}{8}+\binom{10}{9}+\binom{10}{10}\big]/2^{10} = (45+10+1)/1024 = 56/1024 = 0.0547$.

(b) At $p = 0.7$: $P(X = 8) = 45(0.7)^8(0.3)^2 = 0.2335$, $P(X = 9) = 10(0.7)^9(0.3) = 0.1211$, $P(X = 10) = 0.7^{10} = 0.0282$. Power $= 0.3828$; type II error rate $= 1-0.3828 = 0.6172$.

(c) For $X\ge9$: size $= 11/1024 = 0.0107$; power at 0.7 $= 0.1211+0.0282 = 0.1493$; type II error rate 0.8507. Cutting the false-alarm rate from 5.5% to 1.1% cuts the power from 38% to 15%: with $n$ fixed, the two errors trade against each other.

**A4.** Sorted, $m = 8$: 0.001, 0.012, 0.020, 0.024, 0.040, 0.046, 0.30, 0.70.

*Bonferroni:* threshold $0.05/8 = 0.00625$. Only 0.001 passes: **1 rejection**.

*Holm:* $0.001\le0.05/8 = 0.00625$ ✓; $0.012\le0.05/7 = 0.00714$ ✗, stop. **1 rejection.**

*BH:* thresholds $jq/m = 0.00625j$: 0.00625, 0.0125, 0.01875, 0.025, 0.03125, 0.0375, 0.04375, 0.05. Compare: $0.001\le0.00625$ ✓; $0.012\le0.0125$ ✓; $0.020>0.01875$ ✗; $0.024\le0.025$ ✓; $0.040>0.03125$ ✗; $0.046>0.0375$ ✗; $0.30$ ✗; $0.70$ ✗. The largest $j$ that passes is $j = 4$, so BH rejects the four smallest: 0.001, 0.012, 0.020, 0.024. **4 rejections.**

The p-value 0.020 exceeds its own threshold 0.01875 and is still rejected. BH is step-up: it finds the *largest* $j$ with $p_{(j)}\le jq/m$ and rejects everything at or below $p_{(j)}$. Because $p_{(4)} = 0.024$ passes, $p_{(3)}$ comes along.

**A5.** (a) False. The true difference is fixed; [0.4, 3.1] is fixed; there is nothing random left to have probability 0.95. The 95% describes the procedure.
(b) True. That is the definition of coverage.
(c) False. 0.02 is $P(\text{data at least this extreme}\mid H_0)$, not $P(H_0\mid\text{data})$; the latter needs a prior.
(d) True. 0 lies outside the 95% interval, so by duality (Theorem 3.5) the level-0.05 test rejects; consistently, $p = 0.02<0.05$.
(e) False. 3 lies inside [0.4, 3.1], so the level-0.05 test of "difference = 3" does not reject (duality again, assuming the interval inverts the same test).
(f) True. That is the definition of the p-value for this two-sided test.

**A6.** By (3.21) with $z_{0.025} = 1.960$ (two-sided) and $z_{0.10} = 1.2816$: $n = \big((1.960+1.2816)/0.3\big)^2 = (3.2415/0.3)^2 = 10.805^2 = 116.75$, so $n = 117$. For two groups, $\operatorname{Var}(\bar X-\bar Y) = 2\sigma^2/n$, which is the one-sample variance with $n$ replaced by $n/2$, so $n/2 = 116.75$: $n = 233.5$, i.e. 234 per group.

## B. Core

**B1.** (a) $\bar X = (X_1+X_2)/2$, so $X_1-\bar X = (X_1-X_2)/2$ and $X_2-\bar X = -(X_1-X_2)/2$. Then $S^2 = \frac{1}{2-1}\big[(X_1-X_2)^2/4+(X_1-X_2)^2/4\big] = (X_1-X_2)^2/2$.

(b) $(U,V)^\top = Q^\top\mathbf X$ with $Q = \frac{1}{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}$, and $Q^\top Q = I$. Since $\mathbf X\sim\mathcal N(\mu\mathbf 1,\sigma^2I)$, F5 Theorem 3.22 gives $(U,V)\sim\mathcal N\big(Q^\top\mu\mathbf 1,\ \sigma^2Q^\top Q\big) = \mathcal N\big((\sqrt2\mu,0)^\top,\ \sigma^2I\big)$. The covariance is diagonal, so $U$ and $V$ are independent (F5 Theorem 3.23). Now $\bar X = U/\sqrt2$ is a function of $U$ alone, and $S^2 = (X_1-X_2)^2/2 = V^2$ is a function of $V$ alone, so $\bar X$ and $S^2$ are independent. Finally $V/\sigma\sim\mathcal N(0,1)$, so $S^2/\sigma^2 = (V/\sigma)^2\sim\chi^2_1$, which is (3.1) with $n-1 = 1$.

(c) $\mathbf q_1 = (1,1)^\top/\sqrt2 = \mathbf 1/\sqrt2$ and $\mathbf q_2 = (1,-1)^\top/\sqrt2$, the only unit vector orthogonal to $\mathbf 1$ in $\mathbb R^2$ up to sign.

**B2.** (a) Each $X_i\sim\mathrm{Expo}(\lambda) = \mathrm{Gamma}(1,\lambda)$, so $Y = \sum X_i\sim\mathrm{Gamma}(n,\lambda)$ (F5 Theorem 3.8), with density $f_Y(y) = \lambda^ny^{n-1}e^{-\lambda y}/\Gamma(n)$. Let $W = 2\lambda Y$, so $y = w/(2\lambda)$ and $dy = dw/(2\lambda)$:
$$f_W(w) = \frac{\lambda^n}{\Gamma(n)}\Big(\frac{w}{2\lambda}\Big)^{n-1}e^{-w/2}\cdot\frac{1}{2\lambda} = \frac{(1/2)^n}{\Gamma(n)}w^{n-1}e^{-w/2},$$
which is the $\mathrm{Gamma}(n,\frac12) = \chi^2_{2n}$ density (F5 §3.5: $\chi^2_k = \mathrm{Gamma}(k/2,\frac12)$). It involves λ yet has a distribution free of λ: a pivot.

(b) $\sum x_i = 123.0$, $n = 8$, so $\hat\lambda = n/\sum x_i = 8/123 = 0.06504$. With $\chi^2_{16,0.975} = 6.908$ and $\chi^2_{16,0.025} = 28.845$: $P(6.908\le2\lambda\cdot123\le28.845) = 0.95$, so
$$\lambda\in\Big[\frac{6.908}{246},\ \frac{28.845}{246}\Big] = [0.02808,\ 0.11726].$$
Since $\lambda\mapsto1/\lambda$ is decreasing and one-to-one, $\{1/\lambda:\lambda\in C\}$ has the same coverage: mean lifetime $\in[246/28.845,\ 246/6.908] = [8.53,\ 35.61]$ hours.

(c) The interval is equal-tailed in the pivot, and $\chi^2_{16}$ is right-skewed, so equal tail areas sit at unequal distances from the centre: $\hat\lambda-0.02808 = 0.037$ but $0.11726-\hat\lambda = 0.052$. Nothing requires a confidence interval to be symmetric about the estimate; coverage is the only requirement.

**B3.** (a) The likelihood is $f_\lambda(\mathbf x) = \lambda^ne^{-\lambda s}$ with $s = \sum x_i$, so
$$\frac{f_{\lambda_1}(\mathbf x)}{f_{\lambda_0}(\mathbf x)} = \Big(\frac{\lambda_1}{\lambda_0}\Big)^ne^{(\lambda_0-\lambda_1)s}.$$
Because $\lambda_0-\lambda_1>0$ this is increasing in $s$, so $\{f_{\lambda_1}/f_{\lambda_0}>k\} = \{s>c\}$ for some $c$. By Theorem 3.7 the most powerful level-α test rejects when $\sum X_i>c$, with $c$ fixed by $P_{\lambda_0}(\sum X_i>c) = \alpha$. By B2(a), $2\lambda_0\sum X_i\sim\chi^2_{2n}$ under $H_0$, so $P(\chi^2_{2n}>2\lambda_0c) = \alpha$ and $c = \chi^2_{2n,\alpha}/(2\lambda_0)$.

(b) The threshold $c$ does not depend on $\lambda_1$. So one test is most powerful against every $\lambda_1<\lambda_0$, i.e. uniformly most powerful for $H_1:\lambda<\lambda_0$ (the argument of Theorem 3.8). If the null is the composite $\lambda\ge\lambda_0$, the power function $\beta(\lambda) = P(\chi^2_{2n}>2\lambda c)$ decreases in λ, so its supremum over the null is at $\lambda_0$ and the size is still α.

(c) $\chi^2_{16,0.05} = 26.296$, so $c = 26.296/(2\times0.1) = 131.48$. Observed $s = 123.0<131.48$: do not reject. The p-value is $P(\chi^2_{16}\ge2\times0.1\times123) = P(\chi^2_{16}\ge24.6) = 0.0772$. Power at $\lambda_1 = 0.05$ (true mean 20 hours): $P_{0.05}(\sum X_i>131.48) = P(\chi^2_{16}>2\times0.05\times131.48) = P(\chi^2_{16}>13.148) = 0.6619$. Even if the components lasted twice as long as claimed, eight of them would miss it a third of the time. "Not rejected" is not "claim confirmed."

**B4.** Here $\ell(\lambda) = n\log\lambda-\lambda s$, $\ell'(\lambda) = n/\lambda-s$, $\ell''(\lambda) = -n/\lambda^2$, so the total information is $nI_1(\lambda) = n/\lambda^2$ and $\hat\lambda = n/s = 0.06504$. Note $\lambda_0s = 12.3$.

(a) *LRT:* $2[\ell(\hat\lambda)-\ell(\lambda_0)] = 2[n\log(\hat\lambda/\lambda_0)-(\hat\lambda-\lambda_0)s]$. Here $\log(\hat\lambda/\lambda_0) = \log0.65041 = -0.43016$ and $(\hat\lambda-\lambda_0)s = n-\lambda_0s = 8-12.3 = -4.3$, so LRT $= 2[8(-0.43016)+4.3] = 2(0.8587) = 1.7175$.
*Wald:* $(\hat\lambda-\lambda_0)^2n/\hat\lambda^2 = n(1-\lambda_0/\hat\lambda)^2 = 8(1-12.3/8)^2 = 8(0.5375)^2 = 2.3112$.
*Score:* $\ell'(\lambda_0)^2/(n/\lambda_0^2) = (80-123)^2/800 = 1849/800 = 2.3112$. (Here Wald and score coincide, because both equal $(n-\lambda_0s)^2/n$ in this parametrization; that is special to this model.)

(b) With $\theta = 1/\lambda$: $\ell(\theta) = -n\log\theta-s/\theta$, $\hat\theta = s/n = 15.375$, $\ell'(\theta) = -n/\theta+s/\theta^2$, and $nI_1(\theta) = \mathbb E[-\ell''(\theta)] = \mathbb E[-n/\theta^2+2S/\theta^3] = -n/\theta^2+2n\theta/\theta^3 = n/\theta^2$.
*LRT:* $\ell(\theta)$ is $\ell(\lambda)$ evaluated at $\lambda = 1/\theta$, so the maximum value and the value at $\theta_0 = 10$ are the same numbers as before: LRT $= 1.7175$, unchanged.
*Score:* $\ell'(10) = -0.8+123/100 = 0.43$ and $nI_1(10) = 8/100 = 0.08$, so $S = 0.43^2/0.08 = 2.3112$, unchanged. (In general $\ell'_\theta = \ell'_\lambda\,d\lambda/d\theta$ and $I_\theta = I_\lambda(d\lambda/d\theta)^2$, so the ratio $\ell'^2/I$ is invariant.)
*Wald:* $(\hat\theta-\theta_0)^2n/\hat\theta^2 = 8(5.375/15.375)^2 = 8(0.34959)^2 = 0.9777$. Changed, from 2.311 to 0.978. Wald measures distance in the coordinate you chose, scaled by curvature at $\hat\theta$, and a nonlinear change of coordinates distorts that distance differently on the two sides.

(c) $\chi^2_1$ p-values: LRT 0.190; Wald in λ (and score) 0.128; Wald in θ 0.323. Exact: $v = 2\lambda_0s = 24.6$, $P(\chi^2_{16}\ge24.6) = 0.0772$, so the exact two-sided p-value is $2\times0.0772 = 0.154$. The LRT and score p-values bracket the exact value; the Wald p-value is 0.128 or 0.323 depending on an arbitrary choice of parametrization. Report the exact p-value when a pivot exists; among the asymptotic ones, prefer the LRT (or score), never a Wald statistic whose answer depends on how you wrote the model.

**B5.** (a) Under resampling, $X^*_1,\dots,X^*_n$ are iid, each uniform on $\{x_1,\dots,x_n\}$ (with multiplicity). So $\mathbb E^*[X^*_j] = \frac1n\sum_ix_i = \bar x$ and $\operatorname{Var}^*(X^*_j) = \frac1n\sum_i(x_i-\bar x)^2 = \hat\sigma^2$. For an average of $n$ iid draws, $\mathbb E^*[\bar X^*] = \bar x$ and $\operatorname{Var}^*(\bar X^*) = \hat\sigma^2/n$ (F4 Theorem 3.14). ∎

(b) $\mathbf x = (1,2,6)$: $\bar x = 3$, $\hat\sigma^2 = (4+1+9)/3 = 14/3$, so (a) predicts $\operatorname{Var}^*(\bar X^*) = 14/9 = 1.5556$. Count over the 27 ordered triples $(a,b,c)$: each value appears in each position in 9 triples, so $\sum(a+b+c) = 3\times9\times(1+2+6) = 243$ and the mean of $\bar X^* = (a+b+c)/3$ is $243/(27\times3) = 3$ ✓. For the second moment, $\sum(a+b+c)^2 = \sum(a^2+b^2+c^2)+2\sum(ab+bc+ca)$. The first sum is $3\times9\times(1+4+36) = 1107$. For one pair of positions, $\sum ab = 3\times(1+2+6)^2 = 243$ (the third position is free, 3 ways), and there are three pairs, so the second sum is $2\times729 = 1458$. Total $2565$. Then $\mathbb E^*[(\bar X^*)^2] = 2565/(27\times9) = 10.5556$ and $\operatorname{Var}^* = 10.5556-9 = 1.5556 = 14/9$ ✓. The project enumerates all 27.

(c) $\widehat{\rm se}_{\rm boot} = \sqrt{14/9} = 1.2472$, while $s^2 = 14/2 = 7$ gives $s/\sqrt3 = 1.5275$. The ratio is $\hat\sigma/s = \sqrt{(n-1)/n} = \sqrt{2/3} = 0.8165$: the bootstrap uses the plug-in variance, divisor $n$. For large $n$ the factor tends to 1; at $n = 3$ it is an 18% underestimate, one reason the bootstrap is unreliable for tiny samples.

**B6.** Let $\mathcal I_0$ be the true nulls, $m_0 = |\mathcal I_0|$. If $m_0 = 0$, FWER = 0. Otherwise let $i^*$ be the true null with the smallest p-value.

1. Holm examines hypotheses in increasing order of p-value and stops at the first failure. If it rejects any true null, it rejects $H_{i^*}$, because no true null comes before $i^*$ in the order.
2. Suppose $H_{i^*}$ is examined at step $j$. The $j-1$ hypotheses examined earlier have smaller p-values than $p_{i^*}$, so none is a true null; there are only $m-m_0$ false nulls, so $j-1\le m-m_0$, i.e. $m-j+1\ge m_0$.
3. Rejection at step $j$ requires $p_{i^*}\le\alpha/(m-j+1)\le\alpha/m_0$.
4. Hence $\{V\ge1\}\subseteq\{p_{i^*}\le\alpha/m_0\}\subseteq\bigcup_{i\in\mathcal I_0}\{p_i\le\alpha/m_0\}$, and by the union bound and Theorem 3.11, ${\rm FWER}\le m_0\cdot\alpha/m_0 = \alpha$. ∎ No independence was used.

*Holm dominates Bonferroni.* If Bonferroni rejects $H_i$, then $p_i\le\alpha/m$. Every hypothesis Holm examines before or at $H_i$ has $p_{(j)}\le p_i\le\alpha/m\le\alpha/(m-j+1)$, so Holm passes every one of those steps and rejects $H_i$. ∎

**B7.** (a) $p\le u\iff1-\Phi(Z)\le u\iff\Phi(Z)\ge1-u\iff Z\ge z_u$. With $Z\sim\mathcal N(\Delta,1)$, $P(Z\ge z_u) = 1-\Phi(z_u-\Delta) = \Phi(\Delta-z_u)$. At $u = \alpha$ this is the power $\beta = \Phi(\Delta-z_\alpha)$, the same as (3.13) with $\Delta = \sqrt n(\mu-\mu_0)/\sigma$.

(b) $G(u) = \Phi(\Delta-z_u)$ with $z_u = \Phi^{-1}(1-u)$, so $dz_u/du = -1/\varphi(z_u)$. Then
$$g(u) = G'(u) = \varphi(\Delta-z_u)\cdot\frac{1}{\varphi(z_u)} = \exp\Big(-\frac{(z_u-\Delta)^2}{2}+\frac{z_u^2}{2}\Big) = \exp\Big(\Delta z_u-\frac{\Delta^2}{2}\Big).$$
At $\Delta = 0$, $g\equiv1$: the Uniform density, Theorem 3.11 again.

(c) $z_{0.05} = 1.645$, $z_{0.04} = 1.751$, $z_{0.001} = 3.090$. Power $= \Phi(2.5-1.645) = \Phi(0.855) = 0.804$. $g(0.04) = \exp(2.5\times1.751-3.125) = \exp(1.252) = 3.50$; $g(0.001) = \exp(2.5\times3.090-3.125) = \exp(4.601) = 99.5$. Since the null density of $p$ is 1, $g(u)$ is the likelihood ratio of $H_1$ to $H_0$ given "$p$ landed near $u$." A well-powered study reporting $p = 0.04$ has data only 3.5 times more likely under $H_1$: with even prior odds, $P(H_0\mid p\approx0.04) = 1/(1+3.50) = 0.22$. That is weak evidence, much weaker than "4%" suggests. A p-value near 0.001 is a likelihood ratio near 100, which is strong.

**B8.** Reference implementation (bottom of `project.py`):

```python
def slope(a):
    x = a[..., 0] - a[..., 0].mean(1, keepdims=True)
    y = a[..., 1] - a[..., 1].mean(1, keepdims=True)
    return (x * y).sum(1) / (x * x).sum(1)

for _ in range(1000):
    xr = rng.uniform(0, 10, 25); yr = 1 + 0.5 * xr + rng.normal(0, 2, 25)
    ci = bootstrap_ci(np.c_[xr, yr], slope, 1000, rng, kinds=("percentile", "basic"))
    lr = stats.linregress(xr, yr); hw = stats.t.ppf(.975, 23) * lr.stderr
    # count whether 0.5 lies in each interval
```

Coverage over 1,000 datasets (Monte Carlo standard error about 0.007): percentile 0.931, basic 0.927, classical $t$ 0.952. The classical interval is exact here because the errors are Normal (ST-3 derives it with the projection of §3.1 in two dimensions). The bootstrap intervals undercover by about two points at $n = 25$: they are built from $\hat F_n$, whose spread is too small (the B5 effect) and whose quantiles are noisy. That is the usual small-sample price of the bootstrap; it buys validity when the errors are not Normal or the variance is not constant, where the classical interval can fail.

*Why pairs:* each observation is a pair $(x_i,y_i)$ drawn from a joint distribution, and the slope is a property of that joint distribution. Resampling $x$'s and $y$'s separately would pair each $y$ with a random $x$, destroying the relation; the bootstrap slopes would scatter around 0, simulating the null of no relationship rather than the sampling variability of $\hat b$. (The alternative *residual bootstrap* keeps the $x_i$ fixed and resamples the fitted residuals; Efron and Tibshirani ch. 9 compares the two.)

## C. Challenge

**C1.** The joint density of the independent pair $(Z,V)$ is $f(z,v) = \varphi(z)\,\dfrac{v^{k/2-1}e^{-v/2}}{2^{k/2}\Gamma(k/2)}$ for $v>0$.

1. Change variables $(z,v)\mapsto(t,v)$ with $t = z/\sqrt{v/k}$. The inverse is $z = t\sqrt{v/k}$, $v = v$, with Jacobian determinant $|\partial(z,v)/\partial(t,v)| = \sqrt{v/k}$ (the matrix is triangular). By F5 Theorem 3.6, $f_{T,V}(t,v) = \varphi\big(t\sqrt{v/k}\big)\sqrt{v/k}\,f_V(v)$.
2. Integrate out $v$:
$$f_T(t) = \frac{1}{\sqrt{2\pi k}\,2^{k/2}\Gamma(k/2)}\int_0^\infty v^{\frac{k+1}{2}-1}\exp\Big(-\frac v2\Big(1+\frac{t^2}{k}\Big)\Big)dv.$$
3. Use $\int_0^\infty v^{a-1}e^{-bv}dv = \Gamma(a)/b^a$ with $a = (k+1)/2$ and $b = (1+t^2/k)/2$:
$$f_T(t) = \frac{\Gamma\big(\frac{k+1}{2}\big)\,2^{(k+1)/2}}{\sqrt{2\pi k}\,2^{k/2}\Gamma(k/2)}\Big(1+\frac{t^2}{k}\Big)^{-(k+1)/2}.$$
4. The powers of 2 cancel: $2^{(k+1)/2}/(\sqrt2\cdot2^{k/2}) = 1$. So $f_T(t) = \dfrac{\Gamma(\frac{k+1}{2})}{\sqrt{k\pi}\,\Gamma(\frac k2)}\Big(1+\dfrac{t^2}{k}\Big)^{-(k+1)/2}$. ∎

*$k = 1$:* $\Gamma(1)/(\sqrt\pi\,\Gamma(1/2)) = 1/\pi$ (F4 Theorem 3.19(c)), giving $\frac{1}{\pi(1+t^2)}$, the Cauchy density. *$k\to\infty$:* $(1+t^2/k)^{-(k+1)/2}\to e^{-t^2/2}$, and $\Gamma(x+\frac12)/\Gamma(x)\sim\sqrt x$ as $x\to\infty$ (Stirling), so with $x = k/2$ the constant tends to $\sqrt{k/2}/\sqrt{k\pi} = 1/\sqrt{2\pi}$. The limit is $\varphi(t)$. The project checks step 2's integral against `scipy.stats.t.pdf` at three points for $k = 4$.

**C2.** Write $\ell = \sum_{j=1}^kN_j\log p_j$ with $p_k = 1-\sum_{j<k}p_j$, so $\partial p_k/\partial p_j = -1$ for $j<k$.

1. *Score:* $u_j = \partial\ell/\partial p_j = N_j/p_j-N_k/p_k$ for $j<k$.
2. *Information:* $\partial^2\ell/\partial p_j\partial p_l = -\delta_{jl}N_j/p_j^2-N_k/p_k^2$. Taking $-\mathbb E[\cdot]$ with $\mathbb E N_j = np_j$: $I_{jl} = \delta_{jl}\,n/p_j+n/p_k$, i.e. $I = n\big(\operatorname{diag}(1/p_j)+\frac1{p_k}\mathbf 1\mathbf 1^\top\big)$.
3. *Inverse:* let $D = \operatorname{diag}(p_1,\dots,p_{k-1})$ and $\mathbf p = (p_1,\dots,p_{k-1})^\top$. Then $\big(D^{-1}+\frac1{p_k}\mathbf 1\mathbf 1^\top\big)(D-\mathbf p\mathbf p^\top) = I-\mathbf 1\mathbf p^\top+\frac1{p_k}\mathbf 1\mathbf p^\top-\frac{1-p_k}{p_k}\mathbf 1\mathbf p^\top = I$, using $D^{-1}\mathbf p = \mathbf 1$, $\mathbf 1^\top D = \mathbf p^\top$ and $\mathbf 1^\top\mathbf p = 1-p_k$. So $I^{-1} = \frac1n(D-\mathbf p\mathbf p^\top)$.
4. *The quadratic form:* $\mathbf u^\top I^{-1}\mathbf u = \frac1n\Big[\sum_{j<k}p_ju_j^2-\big(\sum_{j<k}p_ju_j\big)^2\Big]$. Put $a_j = N_j/p_j$ for all $j = 1,\dots,k$, so $u_j = a_j-a_k$, and the $j = k$ term $p_k(a_k-a_k)$ is zero. Both sums may therefore run over all $j\le k$. With $A$ a random variable equal to $a_j$ with probability $p_j$, the bracket is $\mathbb E[(A-a_k)^2]-(\mathbb E[A-a_k])^2 = \operatorname{Var}(A) = \sum_jp_ja_j^2-\big(\sum_jp_ja_j\big)^2$.
5. $\sum_jp_ja_j = \sum_jN_j = n$ and $\sum_jp_ja_j^2 = \sum_jN_j^2/p_j$. So the score statistic is $\frac1n\big(\sum_jN_j^2/p_j-n^2\big) = \sum_j\frac{N_j^2}{np_j}-n$.
6. Pearson: $\sum_j\frac{(N_j-np_j)^2}{np_j} = \sum_j\frac{N_j^2}{np_j}-2\sum_jN_j+n\sum_jp_j = \sum_j\frac{N_j^2}{np_j}-n$. The same. ∎ (Everything is evaluated at $\mathbf p = \mathbf p^0$; the project checks it on a 4-cell example.)

**C3.** Reference implementation (bottom of `project.py`): draw a $2\times10^4\times100$ array of $\mathcal N(0,1)$, form $z_j = \sum_{i\le j}X_i/\sqrt j$ at $j = 10,20,\dots,100$, and record whether $|z_j|>1.960$ at any look.

Results: type I error **0.198** with the per-look threshold 0.05, and **0.027** with 0.05/10 = 0.005 per look.

Each look is a valid 5% test, but the analyst rejects if *any* of 10 looks rejects, and the probability of a union of ten 5% events is far above 5%. It is below the $1-0.95^{10} = 0.401$ of independent looks because successive $z_j$ are strongly correlated ($\operatorname{Corr}(z_j,z_l) = \sqrt{j/l}$ for $j<l$, since they share the first $j$ observations). Bonferroni restores control by the union bound (Theorem 3.12), but conservatively (0.027 instead of 0.05) because it ignores that correlation; group-sequential designs (Pocock, O'Brien–Fleming boundaries) spend the 0.05 exactly. With unlimited looks the type I error tends to 1: by the law of the iterated logarithm, $|z_j|$ exceeds any fixed threshold infinitely often under $H_0$.

Why it is p-hacking: the reported p-value is computed as though $n$ had been fixed in advance. The stopping rule makes $n$ depend on the data, so "$p<0.05$ at the final look" has probability 0.198 under $H_0$, not 0.05. Every single computation was correct; the procedure was not.

**C4.** *Derivation.* Under $H_0$, $\bar X\sim\mathcal N(0,1/n)$. Under $H_1$, $\bar X = \mu+e$ with $\mu\sim\mathcal N(0,\tau^2)$ independent of $e\sim\mathcal N(0,1/n)$, so $\bar X\sim\mathcal N(0,\tau^2+1/n)$ (independent Normals add, F5 §3.4). The Bayes factor is the ratio of the two densities at $\bar x$:
$$B_{01} = \frac{\sqrt n\,e^{-n\bar x^2/2}}{(\tau^2+1/n)^{-1/2}e^{-\bar x^2/(2(\tau^2+1/n))}} = \sqrt{1+n\tau^2}\,\exp\Big(-\frac{\bar x^2}{2}\Big(n-\frac{n}{1+n\tau^2}\Big)\Big).$$
Since $n-\frac{n}{1+n\tau^2} = \frac{n^2\tau^2}{1+n\tau^2}$ and $\bar x^2 = z^2/n$, the exponent is $-\frac{z^2}{2}\cdot\frac{n\tau^2}{1+n\tau^2}$, the formula of §3.11.

*Numbers* ($\tau = 1$, $z = 2.5$): $B_{01} = 0.1936$ at $n = 10$; $0.4554$ at $10^2$; $1.394$ at $10^3$; $4.395$ at $10^4$; $43.94$ at $10^6$. Setting $\frac12\log(1+n) = 3.125\,n/(1+n)$ and solving numerically gives $B_{01} = 1$ at $n = 510.7$. The two-sided p-value is $2(1-\Phi(2.5)) = 0.0124$ for every $n$.

*Why they diverge.* The p-value asks only how surprising $\bar x$ is under $H_0$, and a fixed $z$ is equally surprising at every $n$. The Bayes factor also asks how well $H_1$ predicted $\bar x$; at large $n$ the observed $\bar x = z/\sqrt n$ is tiny, which $H_1$ (spreading its prediction over $\mu\in\pm\tau$) predicted poorly, so the point null wins. ST-5 takes this up, including how much the answer depends on τ.
