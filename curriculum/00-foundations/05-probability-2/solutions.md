# F5 Solutions

Numbers are rounded to four places. Every one of them is asserted or printed in `project.py`, under the exercise's label.

## A. Warm-up

**A1.** (a) Row sums give $p_X = (0.4, 0.6)$. Column sums give $p_Y = (0.3, 0.4, 0.3)$. Conditionals by (3.3):
$p_{Y|X}(\cdot\mid 1) = (0.10, 0.20, 0.10)/0.4 = (0.25, 0.5, 0.25)$ and $p_{Y|X}(\cdot\mid 2) = (0.2, 0.2, 0.2)/0.6 = (\tfrac13, \tfrac13, \tfrac13)$.
The two slices have different shapes, so $X$ and $Y$ are dependent (Theorem 3.4(b)). One cell confirms it: $p(1,1) = 0.10 \ne p_X(1)p_Y(1) = 0.12$.

(b) $\mathbb{E}[X] = 0.4 + 2(0.6) = 1.6$. $\mathbb{E}[Y] = 0.3 + 2(0.4) + 3(0.3) = 2.0$.
$\mathbb{E}[XY] = 1\cdot(0.1 + 0.4 + 0.3) + 2\cdot(0.2 + 0.4 + 0.6) = 0.8 + 2.4 = 3.2$.
So $\mathrm{Cov}(X, Y) = 3.2 - 1.6\cdot2.0 = 0$.

(c) $\mathbb{E}[Y\mid X = 1] = 0.25 + 1 + 0.75 = 2$ and $\mathbb{E}[Y\mid X = 2] = (1 + 2 + 3)/3 = 2$. So $\mathbb{E}[Y\mid X] = 2$, a constant random variable.
$\mathbb{E}[Y^2\mid X = 1] = 0.25 + 2 + 2.25 = 4.5$, so $\mathrm{Var}(Y\mid X = 1) = 4.5 - 4 = 0.5$.
$\mathbb{E}[Y^2\mid X = 2] = 14/3$, so $\mathrm{Var}(Y\mid X = 2) = 14/3 - 4 = 2/3$.
Eve's law (3.21): $\mathbb{E}[\mathrm{Var}(Y\mid X)] = 0.4(0.5) + 0.6(\tfrac23) = 0.6$ and $\mathrm{Var}(\mathbb{E}[Y\mid X]) = \mathrm{Var}(2) = 0$. The sum is $0.6$. Directly, $\mathrm{Var}(Y) = (0.3 + 1.6 + 2.7) - 4 = 0.6$. ✓

The "between" term is zero. Knowing $X$ never changes your best guess of $Y$; it changes only how spread $Y$ is around that guess. This explains the zero covariance: $\mathrm{Cov}(X, Y) = \mathbb{E}[XY] - \mathbb{E}X\,\mathbb{E}Y$, and by Adam's law and taking out what is known, $\mathbb{E}[XY] = \mathbb{E}[X\,\mathbb{E}[Y\mid X]] = 2\,\mathbb{E}[X] = \mathbb{E}[X]\mathbb{E}[Y]$. The hierarchy is: independent ⇒ $\mathbb{E}[Y\mid X]$ constant ⇒ uncorrelated. This table shows that neither arrow reverses at the first step.

**A2.** Use $\mathrm{Var}(\mathbf{a}^\top\mathbf{X}) = \mathbf{a}^\top\Sigma\mathbf{a}$ and $\mathrm{Cov}(\mathbf{a}^\top\mathbf{X}, \mathbf{b}^\top\mathbf{X}) = \mathbf{a}^\top\Sigma\mathbf{b}$ (bilinearity).
- $\mathrm{Var}(X+Y) = 4 + 9 + 2(-3) = 7$.
- $\mathrm{Var}(X - Y) = 4 + 9 - 2(-3) = 19$.
- $\mathrm{Var}(2X - 3Y) = 4(4) + 9(9) - 2\cdot2\cdot3\cdot(-3) = 16 + 81 + 36 = 133$.
- The constant shift drops out: $\mathrm{Cov}(2X + 1, 3Y - X) = 6\,\mathrm{Cov}(X,Y) - 2\,\mathrm{Var}(X) = -18 - 8 = -26$.
- $\rho = -3/\sqrt{4\cdot9} = -0.5$.

Eigenvalues: $\det(\Sigma - \lambda I) = \lambda^2 - 13\lambda + 27 = 0$, so $\lambda = (13 \pm\sqrt{61})/2 = 10.405,\ 2.595$. Both are positive, so $\Sigma$ is positive definite, as Theorem 3.19 requires of a covariance matrix (no linear combination of $X, Y$ is constant).

**A3.** (a) $P(\min_i X_i > m) = \prod_i P(X_i > m) = (e^{-\lambda m})^n = e^{-n\lambda m}$ for $m > 0$. This is the survival function of $\mathrm{Exponential}(n\lambda)$. With $n = 5$, $\lambda = 2$, the rate is 10 and the mean is $0.1$.

(b) $F_{\max}(m) = (1 - e^{-m})^2$, so $f_{\max}(m) = 2(1 - e^{-m})e^{-m} = 2e^{-m} - 2e^{-2m}$ for $m > 0$. Its mean is $2\int_0^\infty me^{-m}dm - 2\int_0^\infty me^{-2m}dm = 2(1) - 2(\tfrac14) = \tfrac32$.
Story check: by (a) the minimum is $\mathrm{Exponential}(2)$ with mean $\tfrac12$. After it, the survivor's remaining life is $\mathrm{Exponential}(1)$ by memorylessness (F4), so $\mathbb{E}[\max] = \tfrac12 + 1 = \tfrac32$.

(c) $X + Y \sim \mathrm{Pois}(5)$ (§3.4), so $P(X + Y = 4) = e^{-5}5^4/4! = 625e^{-5}/24 = 0.1755$.

**A4.** (a) $\rho = -2/\sqrt{4\cdot4} = -0.5$. $\det(\Sigma - \lambda I) = (4 - \lambda)^2 - 4 = 0$ gives $\lambda = 6$ with eigenvector $(1, -1)/\sqrt2$ and $\lambda = 2$ with eigenvector $(1, 1)/\sqrt2$. The level curves are ellipses centered at $(1, 2)$. The long axis runs along the anti-diagonal $(1,-1)$ with semi-axis $\propto\sqrt6$, and the short axis along $(1,1)$ with semi-axis $\propto\sqrt2$. Negative correlation tilts the ellipse the other way from Example 4.4.

(b) By (3.19), the mean is $1 + \frac{-2}{4}(4 - 2) = 0$ and the variance is $4 - \frac{(-2)^2}{4} = 3$. So $X_1\mid X_2 = 4 \sim \mathcal{N}(0, 3)$. (Check: $\sigma_1^2(1 - \rho^2) = 4(0.75) = 3$.)

(c) Let $B = \begin{pmatrix}1 & 1\\1 & -1\end{pmatrix}$. By Theorem 3.22, $B\mathbf{X}$ is MVN with mean $B\boldsymbol\mu = (3, -1)$ and covariance $B\Sigma B^\top$:
- $\mathrm{Var}(X_1 + X_2) = 4 + 4 - 4 = 4$;
- $\mathrm{Var}(X_1 - X_2) = 4 + 4 + 4 = 12$;
- $\mathrm{Cov}(X_1 + X_2, X_1 - X_2) = \mathrm{Var}X_1 - \mathrm{Var}X_2 = 0$.

So $(X_1 + X_2, X_1 - X_2) \sim \mathcal{N}((3, -1), \mathrm{diag}(4, 12))$. They are jointly Gaussian and uncorrelated, hence independent by Theorem 3.23. The rows of $B$ are (unnormalized) eigenvectors of $\Sigma$, which is why the covariance came out diagonal: rotating to the principal axes decorrelates. The variances $4 = 2\cdot2$ and $12 = 2\cdot6$ are the eigenvalues times $\|\text{row}\|^2 = 2$.

**A5.** (a) A fair die has $H = \log_2 6 = 2.585$ bits. For $p = (\tfrac12, \tfrac14, \tfrac14)$, $H = \tfrac12\cdot1 + 2\cdot\tfrac14\cdot2 = 1.5$ bits.

(b) $D(p\|u) = \sum_ip_i\log_2\frac{p_i}{1/3} = \sum_ip_i\log_2p_i + \log_23 = \log_23 - H(p) = 1.585 - 1.5 = 0.0850$ bits.
$D(u\|p) = \tfrac13\big[\log_2\tfrac{1/3}{1/2} + 2\log_2\tfrac{1/3}{1/4}\big] = \tfrac13[-0.5850 + 0.8301] = 0.0817$ bits. Not equal to $D(p\|u)$.
$H(p, u) = -\sum_ip_i\log_2\tfrac13 = \log_23 = 1.585$ bits. This equals $H(p) + D(p\|u) = 1.5 + 0.085$, as (3.30) says.

The same computation for any $p$ on $m$ points gives $D(p\|u) = \log m - H(p)$. Gibbs' inequality (Theorem 3.45) says the left side is $\ge 0$, with equality iff $p = u$. So $H(p) \le \log m$, with equality iff $p$ is uniform: Theorem 3.43(b) again, with Gibbs in place of the direct Jensen argument.

## B. Core

**B1.** (a) The new variables are $(r, \theta)$ and the old ones $(z_1, z_2) = (r\cos\theta, r\sin\theta)$. The Jacobian of old with respect to new is
$$\frac{\partial(z_1, z_2)}{\partial(r, \theta)} = \begin{pmatrix}\cos\theta & -r\sin\theta\\ \sin\theta & r\cos\theta\end{pmatrix}, \qquad \det = r\cos^2\theta + r\sin^2\theta = r > 0.$$
The map is a bijection from $(0,\infty)\times(-\pi,\pi)$ onto the plane minus a ray, which has probability 0. By (3.5),
$$f_{R,\Theta}(r, \theta) = \varphi(r\cos\theta)\varphi(r\sin\theta)\cdot r = \frac{1}{2\pi}e^{-r^2/2}\,r = \underbrace{\frac{1}{2\pi}}_{\theta}\cdot\underbrace{re^{-r^2/2}}_{r}.$$
The density factors and the domain is a product, so $R \perp \Theta$ by Theorem 3.4(c). $\Theta$ has the constant density $\tfrac1{2\pi}$ on $(-\pi, \pi)$, so it is uniform. $R$ has density $re^{-r^2/2}$, which integrates to $[-e^{-r^2/2}]_0^\infty = 1$.
Mean: $\mathbb{E}[R] = \int_0^\infty r^2e^{-r^2/2}dr = \tfrac12\int_{-\infty}^\infty r^2e^{-r^2/2}dr = \tfrac12\sqrt{2\pi}\cdot\mathbb{E}[Z^2] = \sqrt{\pi/2} = 1.2533$. The project's sample mean is $1.2542$.
Also $P(R^2 > s) = P(R > \sqrt s) = e^{-s/2}$, so $R^2 \sim \mathrm{Exponential}(\tfrac12)$: Box–Muller read backwards.

(b) The new variables are $(t, w) = (z_1/z_2, z_2)$, defined where $z_2 \ne 0$ (probability 1). The old variables in terms of the new are $z_1 = tw$, $z_2 = w$, so
$$\frac{\partial(z_1, z_2)}{\partial(t, w)} = \begin{pmatrix} w & t\\ 0 & 1\end{pmatrix}, \qquad |\det| = |w|.$$
By (3.5), $f_{T,W}(t, w) = \varphi(tw)\varphi(w)|w| = \frac{1}{2\pi}|w|e^{-w^2(1+t^2)/2}$. Integrate out $w$. The integrand is even in $w$:
$$f_T(t) = 2\int_0^\infty\frac{1}{2\pi}we^{-w^2(1+t^2)/2}dw = \frac1\pi\left[\frac{-e^{-w^2(1+t^2)/2}}{1+t^2}\right]_0^\infty = \frac{1}{\pi(1+t^2)}.$$
The CDF is $F_T(t) = \tfrac12 + \tfrac1\pi\arctan t$. Setting it to $\tfrac14, \tfrac12, \tfrac34$ gives $\arctan t = -\tfrac\pi4, 0, \tfrac\pi4$, so the quartiles are $-1, 0, 1$. The project's sample quartiles are $-0.995, -0.006, 0.996$.

(c) $\int|t|\frac{dt}{\pi(1+t^2)}$ diverges logarithmically, so $T$ has no mean and no variance, and the CLT's hypotheses fail. Worse, by §3.11 the average of $n$ iid standard Cauchys is again standard Cauchy: averaging ratios of noisy quantities whose denominator can be near zero does not concentrate at all.

**B2.** Cross-covariance of random vectors is $\mathrm{Cov}(\mathbf{U}, \mathbf{V}) = \mathbb{E}[(\mathbf{U} - \mathbb{E}\mathbf{U})(\mathbf{V} - \mathbb{E}\mathbf{V})^\top]$. It is bilinear: $\mathrm{Cov}(A\mathbf{U} + \mathbf{c}, \mathbf{V}) = A\,\mathrm{Cov}(\mathbf{U}, \mathbf{V})$, by the same proof as Theorem 3.19.

(a) *Orthogonal implies optimal.* Suppose $\mathbf{R}^* = \mathbf{X}_1 - \mathbf{a}^* - B^*\mathbf{X}_2$ has $\mathbb{E}\mathbf{R}^* = \mathbf{0}$ and $\mathrm{Cov}(\mathbf{R}^*, \mathbf{X}_2) = 0$. Then every entry satisfies $\mathbb{E}[R^*_iX_{2j}] = \mathrm{Cov}(R^*_i, X_{2j}) + \mathbb{E}[R^*_i]\mathbb{E}[X_{2j}] = 0$. For any other $(\mathbf{a}, B)$,
$$\mathbf{X}_1 - \mathbf{a} - B\mathbf{X}_2 = \mathbf{R}^* + \mathbf{D}, \qquad \mathbf{D} = (\mathbf{a}^* - \mathbf{a}) + (B^* - B)\mathbf{X}_2 =: \mathbf{c} + E\mathbf{X}_2.$$
The cross term is $\mathbb{E}[\mathbf{R}^{*\top}\mathbf{D}] = \mathbf{c}^\top\mathbb{E}[\mathbf{R}^*] + \sum_{i,j}E_{ij}\mathbb{E}[R^*_iX_{2j}] = 0$. Hence
$$\mathbb{E}\|\mathbf{X}_1 - \mathbf{a} - B\mathbf{X}_2\|^2 = \mathbb{E}\|\mathbf{R}^*\|^2 + \mathbb{E}\|\mathbf{D}\|^2 \ge \mathbb{E}\|\mathbf{R}^*\|^2.$$
Equality holds iff $\mathbf{c} + E\mathbf{X}_2 = \mathbf{0}$ with probability 1. Each row of that says a linear combination of $\mathbf{X}_2$ is constant. That would make $\Sigma_{22}$ singular (Theorem 3.19), so $E = 0$ and then $\mathbf{c} = \mathbf{0}$. The orthogonal pair is therefore the *unique* minimizer, which gives the "only if" direction too.

*Solving the conditions.* By bilinearity, $\mathrm{Cov}(\mathbf{R}^*, \mathbf{X}_2) = \Sigma_{12} - B^*\Sigma_{22} = 0$, so $B^* = \Sigma_{12}\Sigma_{22}^{-1}$. And $\mathbb{E}\mathbf{R}^* = \boldsymbol\mu_1 - \mathbf{a}^* - B^*\boldsymbol\mu_2 = \mathbf{0}$, so $\mathbf{a}^* = \boldsymbol\mu_1 - B^*\boldsymbol\mu_2$. Such a pair exists, so the minimizer exists.

(b) Expand by bilinearity:
$$\mathrm{Cov}(\mathbf{X}_1 - B^*\mathbf{X}_2) = \Sigma_{11} - B^*\Sigma_{21} - \Sigma_{12}B^{*\top} + B^*\Sigma_{22}B^{*\top}.$$
With $B^* = \Sigma_{12}\Sigma_{22}^{-1}$ and $B^{*\top} = \Sigma_{22}^{-1}\Sigma_{21}$ ($\Sigma$ symmetric), each of the last three terms equals $\Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21}$, with signs $-, -, +$. The result is $\Sigma_{11} - \Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21}$, the Schur complement.

(c) These are exactly the mean and covariance in (3.19), and nothing above used Gaussianity. So for *any* distribution, (3.19)'s formulas give the best *linear* predictor and its error covariance. The Gaussian case adds one thing. $(\mathbf{R}^*, \mathbf{X}_2)$ is a linear transform of $\mathbf{X}$, hence jointly Gaussian (Theorem 3.22), and it is uncorrelated, hence independent (Theorem 3.23). So $\mathbf{X}_1 = \mathbf{a}^* + B^*\mathbf{X}_2 + \mathbf{R}^*$ with noise $\mathbf{R}^* \sim \mathcal{N}(\mathbf{0}, S)$ independent of $\mathbf{X}_2$. Given $\mathbf{X}_2 = \mathbf{x}_2$, then, $\mathbf{X}_1 \sim \mathcal{N}(\mathbf{a}^* + B^*\mathbf{x}_2, S)$. This is a three-line second proof of Theorem 3.24, without completing the square. For non-Gaussian data, $\mathbb{E}[\mathbf{X}_1\mid\mathbf{X}_2]$ can be a curve and the conditional covariance can depend on $\mathbf{x}_2$. The project's B2 block uses Uniform inputs and a quadratic term: least squares reproduces $B^*$, and the residual covariance matches the Schur complement anyway.

**B3.** (a) Adam: $\mathbb{E}[Y] = \mathbb{E}[\mathbb{E}[Y\mid Z]] = 0.7(0) + 0.3(4) = 1.2$.

(b) Eve has two terms.
- Within: $\mathbb{E}[\mathrm{Var}(Y\mid Z)] = 0.7(1) + 0.3(4) = 1.9$.
- Between: $\mathbb{E}[Y\mid Z] = 4Z$, so $\mathrm{Var}(\mathbb{E}[Y\mid Z]) = 16\,\mathrm{Var}(Z) = 16(0.3)(0.7) = 3.36$.

$\mathrm{Var}(Y) = 1.9 + 3.36 = 5.26$. Most of the variance comes from not knowing which component you are in. The project's simulation gives $1.901 + 3.366$.

(c) Condition on $Z$ (law of total probability): $P(Y > 3) = 0.7\,P(\mathcal{N}(0,1) > 3) + 0.3\,P(\mathcal{N}(4, 4) > 3) = 0.7(1 - \Phi(3)) + 0.3(1 - \Phi(-0.5))$. That is $0.7(0.00135) + 0.3(0.69146) = 0.00094 + 0.20744 = 0.2084$.

(d) Bayes' rule with a density for the evidence:
$$P(Z = 1\mid Y = 3) = \frac{0.3\cdot\tfrac12\varphi(-0.5)}{0.7\,\varphi(3) + 0.3\cdot\tfrac12\varphi(-0.5)} = \frac{0.3(0.17603)}{0.7(0.004432) + 0.3(0.17603)} = \frac{0.05281}{0.00310 + 0.05281} = 0.9445.$$
(The $\mathcal{N}(4, 2^2)$ density at 3 is $\tfrac12\varphi(\tfrac{3-4}{2})$.) This is the "responsibility" of component 1 for the point $y = 3$, the E-step of EM in ML-6.

(e) No. Every normal has third central moment 0. For $Y\mid Z = z \sim \mathcal{N}(m, s^2)$ and $d = m - 1.2$, $\mathbb{E}[(Y - 1.2)^3\mid Z = z] = d^3 + 3ds^2$. So $\mathbb{E}[(Y - 1.2)^3] = 0.7(-1.728 - 3.6) + 0.3(21.952 + 33.6) = 12.936 \ne 0$. A mixture of normals is not normal; it is right-skewed here. (The project's Monte Carlo value is $12.957$.)

**B4.** (a) By (3.12) the posterior is $\mathrm{Beta}(1 + 7, 1 + 3) = \mathrm{Beta}(8, 4)$, with mean $8/12 = 2/3$. For the variance, use $\mathbb{E}[p^2] = \frac{B(a+2, b)}{B(a,b)} = \frac{a(a+1)}{(a+b)(a+b+1)}$ (Theorem 3.11 and $\Gamma(x+1) = x\Gamma(x)$). Then
$$\mathrm{Var}(p) = \frac{a(a+1)}{(a+b)(a+b+1)} - \frac{a^2}{(a+b)^2} = \frac{ab}{(a+b)^2(a+b+1)} = \frac{32}{144\cdot13} = 0.01709,$$
and the standard deviation is $0.1307$.

(b) Let $w = \frac{a+b}{a+b+n}$. Then
$$w\cdot\frac{a}{a+b} + (1 - w)\cdot\frac kn = \frac{a}{a+b+n} + \frac{n}{a+b+n}\cdot\frac kn = \frac{a + k}{a+b+n}.$$
So the posterior mean is a weighted average of the prior mean and the sample proportion. The prior's weight is its pseudo-sample size $a + b$ relative to the total. Here $w = 2/12 = 1/6$ on the prior mean $\tfrac12$ and $5/6$ on $0.7$: $\tfrac16(0.5) + \tfrac56(0.7) = 0.6667$. As $n\to\infty$, $w \to 0$ and the data win.

(c) By Theorem 3.9, $\mathrm{Beta}(k, n - k + 1)$ is the $k$-th smallest of $n$ iid Uniforms. Here $k = 8$ and $n - k + 1 = 4$, so $n = 11$: $p$ has the law of $U_{(8)}$ from 11 Uniforms. The event $U_{(8)} > \tfrac12$ means fewer than 8 of the 11 Uniforms fell below $\tfrac12$. The count below $\tfrac12$ is $\mathrm{Bin}(11, \tfrac12)$, so
$$P(p > \tfrac12\mid\text{data}) = P(\mathrm{Bin}(11, \tfrac12) \le 7) = 1 - \frac{\binom{11}{8} + \binom{11}{9} + \binom{11}{10} + \binom{11}{11}}{2^{11}} = 1 - \frac{165 + 55 + 11 + 1}{2048} = \frac{227}{256} = 0.8867.$$

**B5.** (a) Condition on $S$: $P(Y \le y) = \tfrac12P(X \le y) + \tfrac12P(-X \le y) = \tfrac12\Phi(y) + \tfrac12\Phi(y) = \Phi(y)$, because $-X \sim \mathcal{N}(0,1)$ by symmetry.

(b) $\mathbb{E}[XY] = \mathbb{E}[SX^2] = \mathbb{E}[S]\,\mathbb{E}[X^2] = 0$ by independence of $S$ and $X$ (Corollary 3.5). With $\mathbb{E}X = \mathbb{E}Y = 0$, $\mathrm{Cov}(X, Y) = 0$.

(c) $|Y| = |X|$ always. So $P(|X| < 1, |Y| > 1) = 0$, while $P(|X| < 1)\,P(|Y| > 1) > 0$. That violates Definition 3.3, so $X$ and $Y$ are not independent.

(d) $X + Y = (1 + S)X$ is $0$ when $S = -1$ and $2X$ when $S = 1$. So $P(X + Y = 0) = \tfrac12$. A normal variable with positive variance has no atoms, and a zero-variance normal is a constant; $X + Y$ is neither. So the linear combination $X + Y$ is not normal, and $(X, Y)$ is not MVN, because every linear combination of an MVN vector is normal (corollary to Theorem 3.22). Theorem 3.23 needs *joint* Gaussianity. Gaussian marginals plus zero correlation are not enough. The project finds $P(X + Y = 0) = 0.5009$.

**B6.** (a) $\mathrm{Var}(\bar X) = \tfrac14/100 = 0.0025$. Chebyshev gives $P(|\bar X - \tfrac12| \ge 0.1) \le 0.0025/0.01 = 0.25$.

(b) The standard deviation of $\bar X$ is $0.05$, so $0.1$ is two standard deviations. The CLT gives $2(1 - \Phi(2)) = 0.0455$.
With the continuity correction, let $S = 100\bar X \sim \mathrm{Bin}(100, \tfrac12)$, with mean 50 and SD 5. The event is $S \le 40$ or $S \ge 60$, and
$$P \approx \Phi\!\left(\tfrac{40.5 - 50}{5}\right) + 1 - \Phi\!\left(\tfrac{59.5 - 50}{5}\right) = 2\Phi(-1.9) = 0.0574.$$
The exact Binomial value is $0.0569$. Chebyshev is off by more than a factor of 4. The corrected normal approximation is off by $0.0005$.

(c) Chebyshev needs $\frac{1/4}{0.01n} \le 0.05$, i.e. $n \ge 500$. The CLT needs $2(1 - \Phi(0.2\sqrt n)) \le 0.05$, i.e. $0.2\sqrt n \ge 1.96$, so $n \ge 96.04$ and $n = 97$ suffices. The difference is the price of a guarantee that assumes nothing about the distribution's shape against an approximation that uses it.

(d) For Bernoulli(0.1), $\sigma = 0.3$ and $\sigma^3 = 0.027$. Also $\rho_3 = \mathbb{E}|X - 0.1|^3 = 0.9(0.1)^3 + 0.1(0.9)^3 = 0.0009 + 0.0729 = 0.0738$. So $\rho_3/\sigma^3 = 2.7333$, and the bound (3.25) at $n = 100$ is $0.4748\cdot2.7333/\sqrt{100} = 0.1298$. The actual $\sup_x|F_{100} - \Phi|$, computed exactly from the Binomial CDF, is $0.0832$. The bound holds and is loose by a factor of about 1.6. It has the right *rate*: the project's log-log slope for this case is $-0.48$.

**B7.** (a) By linearity, $\mathbb{E}\bar X_n = \mu$. By (3.14) with all cross-covariances zero, $\mathrm{Var}(\bar X_n) = \frac1{n^2}\sum_i\mathrm{Var}(X_i) \le \frac{\sigma^2}{n}$. Chebyshev then gives $P(|\bar X_n - \mu| \ge \epsilon) \le \frac{\sigma^2}{n\epsilon^2} \to 0$. ∎ Independence was never used, only zero covariance.

(b) From §3.12, $\boldsymbol\pi = (\frac{b}{a+b}, \frac{a}{a+b}) = (0.6, 0.4)$ and $\lambda_2 = 1 - a - b = 0.5$. Let $q_t = P(X_t = 1) = 1 - p_t$. Step 3 of the two-state proof, $p_t - p^* = \lambda_2^t(p_0 - p^*)$, gives $q_t - \pi_1 = \lambda_2^t(q_0 - \pi_1)$. Start from $X_0 = 1$ ($q_0 = 1$):
$$P(X_k = 1\mid X_0 = 1) = \pi_1 + \lambda_2^k(1 - \pi_1) = \pi_1 + \lambda_2^k\pi_0.$$
Since $X_0X_k = 1$ exactly when both are 1, $\mathbb{E}[X_0X_k] = \pi_1(\pi_1 + \lambda_2^k\pi_0)$. Started in stationarity, $\mathbb{E}X_0 = \mathbb{E}X_k = \pi_1$. So $\mathrm{Cov}(X_0, X_k) = \lambda_2^k\pi_0\pi_1$. ∎ The chain started in $\boldsymbol\pi$ stays in $\boldsymbol\pi$, and $k$-step probabilities depend only on $k$, so $\mathrm{Cov}(X_s, X_t) = \lambda_2^{|t-s|}\sigma^2$ with $\sigma^2 = \pi_0\pi_1 = 0.24$. The project measures the ratios $0.503, 0.255, 0.129$ for $k = 1, 2, 3$, against $0.5^k$.

(c) By (3.14) applied to $\bar X_n = \frac1n\sum_{s=0}^{n-1}X_s$,
$$\mathrm{Var}(\bar X_n) = \frac{\sigma^2}{n^2}\sum_{s=0}^{n-1}\sum_{t=0}^{n-1}\lambda_2^{|t-s|} \le \frac{\sigma^2}{n^2}\sum_{s=0}^{n-1}\Big(1 + 2\sum_{k=1}^\infty|\lambda_2|^k\Big) = \frac{\sigma^2}{n}\left(1 + \frac{2|\lambda_2|}{1 - |\lambda_2|}\right) = \frac{\sigma^2}{n}\cdot\frac{1 + |\lambda_2|}{1 - |\lambda_2|}.$$
For each fixed $s$, the terms $t = s$ (one of them) and $t = s \pm k$ (at most two for each $k \ge 1$) were bounded by extending the sum to infinity. Here the constant is $0.24\cdot3 = 0.72$, so $\mathrm{Var}(\bar X_n) \le 0.72/n$. Chebyshev then gives $\bar X_n \xrightarrow{p} \pi_1 = 0.4$, and the project's running mean after $2\times10^5$ steps is $0.3995$. The variables are dependent, yet the law of large numbers holds. This is Markov's 1906 argument against Nekrasov's claim that the law requires independence.

**B8.** (a) The rows are $(0, 1, 0, 0)$, $(\tfrac12, 0, \tfrac12, 0)$, $(0, \tfrac12, 0, \tfrac12)$, $(0, 0, 1, 0)$. The degrees are $(1, 2, 2, 1)$ and $2|E| = 6$, so $\boldsymbol\pi = (\tfrac16, \tfrac13, \tfrac13, \tfrac16)$ by the degree formula of §3.12.

(b) The graph is symmetric under $1\leftrightarrow4$, $2\leftrightarrow3$, so look for eigenvectors that are symmetric or antisymmetric under it.
- Symmetric $\mathbf{v} = (x, y, y, x)$. Row 1 gives $y = \lambda x$ and row 2 gives $\tfrac12x + \tfrac12y = \lambda y$. Hence $\lambda^2 - \tfrac12\lambda - \tfrac12 = 0$, so $\lambda = 1$ or $-\tfrac12$.
- Antisymmetric $\mathbf{v} = (x, y, -y, -x)$. Row 1 gives $y = \lambda x$ and row 2 gives $\tfrac12x - \tfrac12y = \lambda y$. Hence $\lambda^2 + \tfrac12\lambda - \tfrac12 = 0$, so $\lambda = \tfrac12$ or $-1$.

The eigenvalues are $\{1, \tfrac12, -\tfrac12, -1\}$. The eigenvalue $-1$, with eigenvector $(1, -1, 1, -1)$, is the period-2 signature: its component in $\boldsymbol\pi_0P^t$ is multiplied by $(-1)^t$ and never decays.

(c) Split the vertices into $A = \{1, 3\}$ and $A^c = \{2, 4\}$. Every step crosses between them. Starting at 1, $\boldsymbol\pi_t$ is supported on $A$ at even $t$ and on $A^c$ at odd $t$, while $\boldsymbol\pi(A) = \tfrac16 + \tfrac13 = \tfrac12$. If $C$ is the class holding the support at time $t$, the triangle inequality gives
$$\|\boldsymbol\pi_t - \boldsymbol\pi\|_1 \ge \Big|\sum_{i\in C}(\pi_{t,i} - \pi_i)\Big| + \Big|\sum_{i\notin C}(\pi_{t,i} - \pi_i)\Big| = |1 - \tfrac12| + |0 - \tfrac12| = 1.$$
For example, $\boldsymbol\pi_2 = (\tfrac12, 0, \tfrac12, 0)$ is at distance exactly 1, and the project prints $1.000$ at $t = 20$. There is no contradiction with Theorem 3.40, because the chain is not aperiodic: returns to vertex 1 happen only at even times. $\boldsymbol\pi$ is still stationary. The long-run *fraction of time* in each vertex still converges to $\boldsymbol\pi$; that part needs only irreducibility (Levin & Peres, ch. 4). What fails is convergence of $\boldsymbol\pi_t$ itself.

(d) $\boldsymbol\pi\tilde P = \tfrac12(\boldsymbol\pi + \boldsymbol\pi P) = \boldsymbol\pi$, so $\tilde P$ has the same $\boldsymbol\pi$. Any eigenvector of $P$ with eigenvalue $\lambda$ is an eigenvector of $\tilde P$ with eigenvalue $\tfrac12(1 + \lambda)$, which gives $\{1, \tfrac34, \tfrac14, 0\}$. The $-1$ has become $0$. All other eigenvalues lie in $[0, 1)$, the positive diagonal makes the chain aperiodic, and $\|\boldsymbol\pi_t - \boldsymbol\pi\|_1$ decays like $0.75^t$. The project gets $0.0042$ after 20 lazy steps, against $1.000$ for the plain walk.

**B9.** (a) The marginals are $p_Y = (0.95, 0.05)$ and $p_X = (0.9, 0.1)$.
- $H(Y) = -0.95\log_20.95 - 0.05\log_20.05 = 0.0703 + 0.2161 = 0.2864$ bits.
- Given $X = 0$, $Y = 0$ surely, so $H(Y\mid X = 0) = 0$.
- Given $X = 1$, $Y$ is a fair coin, so $H(Y\mid X = 1) = 1$ bit.
- $H(Y\mid X) = 0.9(0) + 0.1(1) = 0.1$ bits.

(b) Learning the rare event $X = 1$ makes $Y$ *more* uncertain (1 bit against 0.2864). But averaged over $X$, conditioning cannot hurt:
$$I(X;Y) = \sum_{x,y}p(x,y)\log\frac{p(x,y)}{p(x)p(y)} = \sum_{x,y}p(x,y)\log\frac{p(y\mid x)}{p(y)} = -H(Y\mid X) + H(Y).$$
The middle step uses $p(x,y)/p(x) = p(y\mid x)$. The last step uses the definition of $H(Y\mid X)$ and, for the second term, $\sum_xp(x,y) = p(y)$ by (3.1). Since $I = D(p_{X,Y}\|p_Xp_Y) \ge 0$ by Gibbs (Theorem 3.45), $H(Y\mid X) \le H(Y)$, with equality iff $X \perp Y$. ∎

(c) $I(X;Y) = 0.2864 - 0.1 = 0.1864$ bits. $H(X, Y) = H(0.9, 0.05, 0.05) = 0.9(0.1520) + 2(0.05)(4.3219) = 0.5690$ bits. Chain rule: $H(X) + H(Y\mid X) = h(0.1) + 0.1 = 0.4690 + 0.1 = 0.5690$. ✓

**B10.** (a) For $k \ne j$, $P_{kj} = Q_{kj}\min(1, \pi_j/\pi_k)$. Multiply by $\pi_k > 0$:
$$\pi_kP_{kj} = Q_{kj}\min(\pi_k, \pi_j).$$
The right side is symmetric in $(k, j)$, because $Q_{kj} = Q_{jk}$ and $\min$ is symmetric. So $\pi_kP_{kj} = \pi_jP_{jk}$ for all $k \ne j$; for $k = j$ it is trivial. The diagonal entries absorb all rejections and play no role in (3.27). By Theorem 3.41, $\boldsymbol\pi$ is stationary. ∎ Only the ratio $\pi_j/\pi_k$ entered, so the normalizing constant of $\boldsymbol\pi$ is never needed.

(b) Here $\boldsymbol\pi = (1, 2, 3, 4, 5)/15$. Moving up is always accepted ($\pi_{k+1}/\pi_k > 1$), so $P_{k,k+1} = \tfrac12$ for $k \le 4$. Moving down from $k$ is accepted with probability $(k-1)/k$, so $P_{k,k-1} = \frac{k-1}{2k}$. From 1, the downward proposal leaves the space and is rejected; from 5, the upward one is. Hence
$$P = \begin{pmatrix} \tfrac12 & \tfrac12 & 0 & 0 & 0\\ \tfrac14 & \tfrac14 & \tfrac12 & 0 & 0\\ 0 & \tfrac13 & \tfrac16 & \tfrac12 & 0\\ 0 & 0 & \tfrac38 & \tfrac18 & \tfrac12\\ 0 & 0 & 0 & \tfrac25 & \tfrac35\end{pmatrix}.$$
Column 3 of $\boldsymbol\pi P$: $\tfrac{2}{15}\cdot\tfrac12 + \tfrac3{15}\cdot\tfrac16 + \tfrac4{15}\cdot\tfrac38 = \tfrac{1 + 0.5 + 1.5}{15} = \tfrac3{15} = \pi_3$. ✓ The project checks all five columns and the full detailed-balance matrix.

(c) The states are $0, \dots, 9$ in code, standing for $k = 1, \dots, 10$. The log-target may be off by any additive constant:
```python
k10 = np.arange(1, 11, dtype=float)
chain_10 = metropolis_discrete(2 * np.log(k10), 10, 400_000, rng)
freq_10 = np.bincount(chain_10, minlength=10) / len(chain_10)
assert np.abs(freq_10 - k10**2 / 385).max() < 0.01
```
The normalizing constant $\sum_{k=1}^{10}k^2 = 385$ appears only in the *check*, never in the sampler. The sampler sees only differences of log-targets, which is the point of (a). The project's maximum deviation is $0.0007$.

## C. Challenge

**C1.** *Lemma.* If $\mathbf{V}$ has mean $\mathbf{m}$ and covariance $C$, then $\mathbb{E}[\mathbf{V}^\top A\mathbf{V}] = \mathrm{tr}(AC) + \mathbf{m}^\top A\mathbf{m}$.
*Proof.* A scalar is its own trace, and $\mathrm{tr}(XY) = \mathrm{tr}(YX)$, so $\mathbf{V}^\top A\mathbf{V} = \mathrm{tr}(A\mathbf{V}\mathbf{V}^\top)$. The trace is linear, so it commutes with $\mathbb{E}$. Expanding (3.15) gives $\mathbb{E}[\mathbf{V}\mathbf{V}^\top] = C + \mathbf{m}\mathbf{m}^\top$. Hence $\mathbb{E}[\mathbf{V}^\top A\mathbf{V}] = \mathrm{tr}(AC) + \mathrm{tr}(A\mathbf{m}\mathbf{m}^\top) = \mathrm{tr}(AC) + \mathbf{m}^\top A\mathbf{m}$. ∎

*Derivation.* By (3.17), $\ln f_i(\mathbf{x}) = -\tfrac k2\ln2\pi - \tfrac12\ln\det\Sigma_i - \tfrac12(\mathbf{x} - \boldsymbol\mu_i)^\top\Sigma_i^{-1}(\mathbf{x} - \boldsymbol\mu_i)$. Take $\mathbf{X} \sim f_1$:
$$D = \mathbb{E}_1[\ln f_1 - \ln f_2] = \tfrac12\ln\frac{\det\Sigma_2}{\det\Sigma_1} - \tfrac12\mathbb{E}_1\big[(\mathbf{X} - \boldsymbol\mu_1)^\top\Sigma_1^{-1}(\mathbf{X} - \boldsymbol\mu_1)\big] + \tfrac12\mathbb{E}_1\big[(\mathbf{X} - \boldsymbol\mu_2)^\top\Sigma_2^{-1}(\mathbf{X} - \boldsymbol\mu_2)\big].$$
- First expectation: the lemma with $\mathbf{V} = \mathbf{X} - \boldsymbol\mu_1$ (mean $\mathbf{0}$, covariance $\Sigma_1$) gives $\mathrm{tr}(\Sigma_1^{-1}\Sigma_1) = \mathrm{tr}(I_k) = k$.
- Second: $\mathbf{V} = \mathbf{X} - \boldsymbol\mu_2$ has mean $\boldsymbol\mu_1 - \boldsymbol\mu_2$ and covariance $\Sigma_1$, so the lemma gives $\mathrm{tr}(\Sigma_2^{-1}\Sigma_1) + (\boldsymbol\mu_1 - \boldsymbol\mu_2)^\top\Sigma_2^{-1}(\boldsymbol\mu_1 - \boldsymbol\mu_2)$.

Substituting yields the stated formula. (The quadratic form is unchanged by flipping the sign of $\boldsymbol\mu_1 - \boldsymbol\mu_2$.)

*$k = 1$.* $\tfrac12\big[\tfrac{\sigma_1^2}{\sigma_2^2} + \tfrac{(\mu_1 - \mu_2)^2}{\sigma_2^2} - 1 + \ln\tfrac{\sigma_2^2}{\sigma_1^2}\big] = \ln\tfrac{\sigma_2}{\sigma_1} + \tfrac{\sigma_1^2 + (\mu_1 - \mu_2)^2}{2\sigma_2^2} - \tfrac12$, which is (3.32). ✓

*Numbers.* With $\Sigma_2 = I$: $\mathrm{tr}(\Sigma_1) = 4$, $\|\boldsymbol\mu_2 - \boldsymbol\mu_1\|^2 = 1$, $k = 2$ and $\ln(1/\det\Sigma_1) = \ln\tfrac13 = -1.0986$. So $D = \tfrac12(4 + 1 - 2 - 1.0986) = 0.9507$ nats. The project's Monte Carlo estimate is $0.9482$.

**C2.** *The identity.* Integrating out $y$ by (3.2) gives $\iint f(x,y)\ln f_X(x)\,dx\,dy = \int f_X(x)\ln f_X(x)\,dx = -h(X)$, and likewise for $Y$. So
$$I(X;Y) = \iint f\ln\frac{f(x,y)}{f_X(x)f_Y(y)} = \iint f\ln f - \iint f\ln f_X - \iint f\ln f_Y = -h(X,Y) + h(X) + h(Y).$$

*Evaluation.* By (3.31), $h(X) = h(Y) = \tfrac12\ln(2\pi e)$. From §3.13, $h(X, Y) = \tfrac12\ln\big((2\pi e)^2\det\Sigma\big)$ with $\det\Sigma = 1 - \rho^2$. So
$$I = \ln(2\pi e) - \ln(2\pi e) - \tfrac12\ln(1 - \rho^2) = -\tfrac12\ln(1 - \rho^2).$$
At $\rho = 0.5$: $-\tfrac12\ln0.75 = 0.1438$ nats $= 0.2075$ bits, against the project's Monte Carlo value of $0.1437$.

*Limits.* As $\rho \to \pm1$, $I \to \infty$. In the limit $Y = \pm X$ exactly, and a continuous variable pins down another continuous variable to infinitely many digits. At $\rho = 0$, $I = 0$, which matches Theorem 3.23 (uncorrelated jointly Gaussian variables are independent). For a Gaussian pair, $\rho$ captures *all* the dependence. For non-Gaussian pairs it does not ($X$ and $X^2$; "Your turn" item 3).

**C3.** (a) Maximize $H(p) = -\sum_ip_i\ln p_i$ subject to $\sum_ip_i = 1$. With $\mathcal{L} = H(p) - \alpha(\sum_ip_i - 1)$, setting $\partial\mathcal{L}/\partial p_i = -\ln p_i - 1 - \alpha = 0$ forces $p_i = e^{-1-\alpha}$, the same for every $i$. So $p_i = 1/m$.
*Global maximum.* Each $-p\ln p$ is concave (second derivative $-1/p < 0$), so $H$ is concave and lies below its tangent plane. For any PMF $q$,
$$H(q) \le H(p) + \sum_i(-\ln p_i - 1)(q_i - p_i) = H(p) + (\ln m - 1)\sum_i(q_i - p_i) = H(p),$$
because $\sum_iq_i = \sum_ip_i = 1$. So the maximum is $H = \ln m$.

(b) Add the constraint $\sum_iip_i = 4.5$ with multiplier $\beta$. Stationarity gives $-\ln p_i - 1 - \alpha - \beta i = 0$, so $p_i = e^{-\beta i}/Z(\beta)$ with $Z(\beta) = \sum_{i=1}^6e^{-\beta i}$. The mean is
$$m(\beta) = \sum_i ip_i = -\frac{d}{d\beta}\ln Z, \qquad m'(\beta) = -\frac{d^2}{d\beta^2}\ln Z = -\mathrm{Var}_\beta(i) < 0,$$
so $m$ is strictly decreasing. Since $m(0) = 3.5 < 4.5$, there is exactly one root, and it has $\beta < 0$. A root-finder (the project uses `brentq`) gives $\beta = -0.3710$. Then $p = (0.0544, 0.0788, 0.1142, 0.1654, 0.2398, 0.3475)$, each entry $e^{0.371} = 1.449$ times the previous one, and $H(p) = 2.328$ bits. That is below $\log_26 = 2.585$: the mean constraint costs a quarter of a bit.

(c) Let $q$ be any PMF on $\{1, \dots, 6\}$ with mean $4.5$. Since $\ln p_i = -\beta i - \ln Z$ is *linear* in $i$,
$$0 \le D(q\|p) = -H(q) - \sum_iq_i\ln p_i = -H(q) + 4.5\beta + \ln Z.$$
The same computation with $q = p$ gives $H(p) = 4.5\beta + \ln Z$. Hence $H(q) \le H(p)$, with equality iff $q = p$ (Gibbs). The project also tests 500 random mean-4.5 alternatives, and none beats $p$.
*SM-4.* Replace face $i$ by the energy $E_i$ of microstate $i$ and the mean $4.5$ by a fixed mean energy. The same two steps give the Boltzmann distribution $p_i = e^{-\beta E_i}/Z$. The multiplier is $\beta = 1/k_BT$, $Z$ is the partition function, and $\langle E\rangle = -\partial_\beta\ln Z$ is the formula for $m(\beta)$ above. Here $\beta < 0$ because the required mean is *above* the uniform mean: the distribution favors high "energies," the analogue of a negative temperature.
