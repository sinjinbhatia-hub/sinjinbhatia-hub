# ST-1 Solutions

Every number below is printed or asserted by `project.py` (the cell "Numbers used in the exercise solutions", Part 8, and the reference cells at the bottom).

## A. Warm-up

**A1.** Use linearity of expectation and, for independent variables, $\mathrm{Var}(aX_1 + bX_2) = (a^2 + b^2)\sigma^2$ (F5 (3.14)).

- $\hat\theta_a$: $\mathbb{E} = \mu$, bias 0; $\mathrm{Var} = (\frac14 + \frac14)\sigma^2 = \sigma^2/2$; $\mathrm{MSE} = \sigma^2/2$.
- $\hat\theta_b$: $\mathbb{E} = (\mu + 2\mu)/3 = \mu$, bias 0; $\mathrm{Var} = (\frac19 + \frac49)\sigma^2 = 5\sigma^2/9$; $\mathrm{MSE} = 5\sigma^2/9$.
- $\hat\theta_c$: $\mathbb{E} = 2\mu/3$, bias $-\mu/3$; $\mathrm{Var} = (\frac19 + \frac19)\sigma^2 = 2\sigma^2/9$; $\mathrm{MSE} = \mu^2/9 + 2\sigma^2/9$ by (3.2).

Of the unbiased pair, $\hat\theta_a$ is better ($\frac12 < \frac59$): equal weights minimize the variance of a weighted average of iid variables. $\hat\theta_c$ beats $\hat\theta_a$ iff $\frac{\mu^2}{9} + \frac{2\sigma^2}{9} < \frac{\sigma^2}{2}$, i.e., $\mu^2 + 2\sigma^2 < \frac92\sigma^2$, i.e., $\mu^2 < \frac52\sigma^2$. (The project checks this identity with SymPy.) The lesson: MSE comparisons usually depend on the unknown parameter, so "minimize MSE" alone does not pick an estimator. You need an extra principle: unbiasedness, a worst case over $\theta$, or an average over $\theta$ (a prior, §3.11).

**A2.** (a) $\ell(p) = \sum_i\big[(x_i - 1)\log(1 - p) + \log p\big] = (\sum x_i - n)\log(1-p) + n\log p$. With $n = 5$ and $\sum x_i = 15$: $\ell'(p) = -\frac{10}{1-p} + \frac5p = 0$ gives $5(1 - p) = 10p$, so $\hat p = 1/3$. In general $\hat p = n/\sum x_i = 1/\bar x$. Check: $\ell''(p) = -\frac{10}{(1-p)^2} - \frac{5}{p^2} < 0$.

(b) By invariance (Theorem 3.11): the MLE of $1/p$ is $1/\hat p = 3 = \bar x$, and the MLE of $P(X = 1) = p$ is $\hat p = 1/3$.

(c) $\log p(k;p) = (k - 1)\log(1-p) + \log p$, so $-\partial_p^2\log p = \frac{k-1}{(1-p)^2} + \frac{1}{p^2}$. With $\mathbb{E}[X - 1] = \frac1p - 1 = \frac{1-p}{p}$:
$$I_1(p) = \frac{1-p}{p(1-p)^2} + \frac1{p^2} = \frac{1}{p(1-p)} + \frac{1}{p^2} = \frac{p + (1 - p)}{p^2(1-p)} = \frac{1}{p^2(1-p)}.$$
Wald standard error: $\sqrt{\hat p^2(1 - \hat p)/n} = \sqrt{\frac19\cdot\frac23\cdot\frac15} = \sqrt{2/135} = 0.1217$.

**A3.** (a) $\mathbb{E}X = 0$ for every $\theta$, so $m_1(\theta) = \hat m_1$ either has no solution or holds for every $\theta$: the first moment carries no information about $\theta$.

(b) $\mathbb{E}X^2 = \int_{-\theta}^\theta\frac{x^2}{2\theta}dx = \frac{\theta^2}{3}$. Setting $\theta^2/3 = \hat m_2$ gives $\hat\theta_{\text{MoM}} = \sqrt{3\hat m_2}$.

(c) $\hat m_2 = (1.44 + 0.16 + 4.41 + 0.49 + 2.25)/5 = 8.75/5 = 1.75$, so $\hat\theta_{\text{MoM}} = \sqrt{5.25} = 2.2913$. The likelihood is $L(\theta) = (2\theta)^{-n}$ if $\theta \ge \max|x_i|$ and 0 otherwise (some point would be impossible). It is decreasing on $[\max|x_i|, \infty)$, so $\hat\theta_{\text{MLE}} = \max|x_i| = 2.1$. This is §3.5(e) with $|X_i| \sim \mathrm{Uniform}(0, \theta)$.

**A4.** (a) $\log p = \log x - \log\theta - \frac{x^2}{2\theta}$. Then
$$s = \partial_\theta\log p = -\frac1\theta + \frac{x^2}{2\theta^2}, \qquad \partial_\theta^2\log p = \frac{1}{\theta^2} - \frac{x^2}{\theta^3}.$$
From $X^2/(2\theta) \sim \mathrm{Exponential}(1)$: $\mathbb{E}X^2 = 2\theta$ and $\mathrm{Var}(X^2) = (2\theta)^2 = 4\theta^2$.
- Curvature form: $-\mathbb{E}[\partial_\theta^2\log p] = -\frac1{\theta^2} + \frac{2\theta}{\theta^3} = \frac{1}{\theta^2}$.
- Score form: $\mathbb{E}s = -\frac1\theta + \frac{2\theta}{2\theta^2} = 0$, so $\mathbb{E}s^2 = \mathrm{Var}(s) = \mathrm{Var}\big(\frac{X^2}{2\theta^2}\big) = \frac{4\theta^2}{4\theta^4} = \frac1{\theta^2}$.

Both give $I_1 = 1/\theta^2$ (the project confirms both integrals with SymPy).

(b) $\ell(\theta) = \sum\log x_i - n\log\theta - \frac{1}{2\theta}\sum x_i^2$, $\ell'(\theta) = -\frac n\theta + \frac{1}{2\theta^2}\sum x_i^2 = 0$, so $\hat\theta = \frac{1}{2n}\sum x_i^2$. It is unbiased: $\mathbb{E}\hat\theta = \frac{1}{2n}\cdot n\cdot 2\theta = \theta$. Its variance is $\frac{1}{4n^2}\cdot n\cdot 4\theta^2 = \frac{\theta^2}{n} = \frac{1}{nI_1}$: it attains the bound. (Theorem 3.17 predicts this: $s_n = -\frac n\theta + \frac{1}{2\theta^2}\sum x_i^2 = \frac{n}{\theta^2}(\hat\theta - \theta) = nI_1(\hat\theta - \theta)$.)

**A5.** $(1-p)^{k-1}p = \exp\big(k\log(1-p) + \log p - \log(1-p)\big)$. So $h(k) = 1$, $T(k) = k$, $\eta = \log(1-p)$ (which is negative), and $A = -\log\frac{p}{1-p}$. In terms of $\eta$: $1 - p = e^\eta$, $p = 1 - e^\eta$, so
$$A(\eta) = \log\frac{e^\eta}{1 - e^\eta} = \eta - \log(1 - e^\eta).$$
Then $A'(\eta) = 1 + \frac{e^\eta}{1 - e^\eta} = \frac{1}{1 - e^\eta} = \frac1p$, the mean. And $A''(\eta) = \frac{e^\eta}{(1 - e^\eta)^2} = \frac{1-p}{p^2}$, the variance. Both match the known values (the project checks both with SymPy).

**A6.** (a) $p_n = \prod(1-p)^{x_i - 1}p = (1-p)^{\sum x_i - n}p^n$: $g(t;p) = (1-p)^{t-n}p^n$, $h = 1$, so $T = \sum x_i$.

(b) $p_n = \theta^n\big(\prod x_i\big)^{\theta - 1}$ on $(0,1)^n$: $g(t;\theta) = \theta^nt^{\theta-1}$ with $t = \prod x_i$, $h = \mathbb{1}\{\text{all } x_i \in (0,1)\}$. So $T = \prod x_i$, equivalently $\sum\log x_i$.

(c) $p_n = \prod\mathbb{1}\{\theta \le x_i \le \theta + 1\} = \mathbb{1}\{\theta \le \min x_i\}\cdot\mathbb{1}\{\max x_i \le \theta + 1\}$, a function of $\theta$ and $(\min x_i, \max x_i)$ only, with $h = 1$. So $T = (X_{(1)}, X_{(n)})$: two dimensions for a one-dimensional parameter. (It is minimal: the likelihood ratio criterion of §3.9 needs both endpoints. One-parameter exponential families always have a one-dimensional sufficient statistic; this family is not one.)

## B. Core

**B1.** Let $\mathrm{SS} = \sum(X_i - \bar X)^2$. By Theorem 3.7, $\mathrm{SS}/\sigma^2 \sim \chi^2_{n-1}$, so $\mathbb{E}\,\mathrm{SS} = (n-1)\sigma^2$ and $\mathrm{Var}(\mathrm{SS}) = 2(n-1)\sigma^4$.

(a) $\mathbb{E}\hat\sigma^2_c = c(n-1)\sigma^2$, so the bias is $(c(n-1) - 1)\sigma^2$. $\mathrm{Var}(\hat\sigma^2_c) = 2c^2(n-1)\sigma^4$. By (3.2):
$$\frac{\mathrm{MSE}(c)}{\sigma^4} = \big(c(n-1) - 1\big)^2 + 2c^2(n-1).$$

(b) Differentiate in $c$: $2(n-1)\big(c(n-1) - 1\big) + 4c(n-1) = 0$. Divide by $2(n-1)$: $c(n-1) - 1 + 2c = 0$, so $c^* = 1/(n+1)$. The second derivative $2(n-1)^2 + 4(n-1) > 0$ confirms a minimum. At $c^*$: $c^*(n-1) - 1 = -\frac{2}{n+1}$, so
$$\frac{\mathrm{MSE}(c^*)}{\sigma^4} = \frac{4}{(n+1)^2} + \frac{2(n-1)}{(n+1)^2} = \frac{2(n+1)}{(n+1)^2} = \frac{2}{n+1}.$$

(c) At $n = 10$: $c = \frac19$ gives $(1 - 1)^2 + \frac{2\cdot 9}{81} = \frac29 = 0.2222$; $c = \frac1{10}$ gives $0.01 + 0.18 = 0.19$; $c = \frac1{11}$ gives $\frac{2}{11} = 0.1818$.

(d) Yes within this class: $\mathrm{MSE}(c)/\sigma^4$ does not involve $\sigma^2$, so $c^*$ is best at every $\sigma^2$. Nobody uses it because (1) the gain is small, $\frac{2}{n+1}$ against $\frac{2}{n-1}$, which matters only for tiny $n$; and (2) $c^*$ was derived from $\mathrm{Var}(\chi^2_{n-1}) = 2(n-1)$, which is a normal-data fact: for heavier-tailed data the variance of SS is larger (it involves the fourth moment), and the optimal divisor changes. $S^2$ is unbiased for any distribution with a variance (Theorem 3.6), which is why it is the convention.

**B2.** From Example 4.4: $\mathbb{E}X_{(n)} = \frac{n\theta}{n+1}$, $\mathbb{E}X_{(n)}^2 = \frac{n\theta^2}{n+2}$.

(a) $\mathrm{MSE}(c) = \mathbb{E}(cX_{(n)} - \theta)^2 = c^2\frac{n\theta^2}{n+2} - 2c\frac{n\theta^2}{n+1} + \theta^2$. Setting the $c$-derivative to zero: $2c\frac{n}{n+2} = 2\frac{n}{n+1}$, so $c^* = \frac{n+2}{n+1}$. Then
$$\frac{\mathrm{MSE}(c^*)}{\theta^2} = \frac{(n+2)^2}{(n+1)^2}\cdot\frac{n}{n+2} - 2\frac{n+2}{n+1}\cdot\frac{n}{n+1} + 1 = 1 - \frac{n(n+2)}{(n+1)^2} = \frac{1}{(n+1)^2}.$$

(b) At $n = 10$, in units of $\theta^2$:
- $2\bar X$: unbiased, variance $\frac{1}{3n} = 0.03333$.
- $X_{(n)}$: bias $-\frac{1}{n+1}$, variance $\frac{n}{(n+1)^2(n+2)}$, so $\mathrm{MSE} = \frac{(n+2) + n}{(n+1)^2(n+2)} = \frac{2}{(n+1)(n+2)} = \frac{2}{132} = 0.01515$.
- $\frac{n+1}{n}X_{(n)}$: unbiased, variance $\frac{1}{n(n+2)} = \frac{1}{120} = 0.00833$.
- $c^*X_{(n)}$ with $c^* = 12/11 = 1.0909$: $\frac{1}{121} = 0.00826$.

(c) $2\bar X$ is an average, and averages fluctuate on the CLT scale $1/\sqrt n$, so its MSE is order $1/n$. The maximum fluctuates on a different scale: the gap $\theta - X_{(n)}$ is the distance from the largest of $n$ uniform points to the edge, which is typically about $\theta/(n+1)$ (one of $n+1$ roughly equal gaps). An error of order $1/n$ squares to an MSE of order $1/n^2$. Regular models cannot do this because their information grows only linearly in $n$; the edge of the support is information that the smooth theory has no way to count.

**B3.** (a) With $S \sim \mathrm{Gamma}(n, \lambda)$ and $k < n$,
$$\mathbb{E}[S^{-k}] = \int_0^\infty s^{-k}\frac{\lambda^ns^{n-1}e^{-\lambda s}}{\Gamma(n)}ds = \frac{\lambda^n}{\Gamma(n)}\cdot\frac{\Gamma(n-k)}{\lambda^{n-k}} = \lambda^k\frac{\Gamma(n-k)}{\Gamma(n)},$$
using $\int_0^\infty s^{a-1}e^{-\lambda s}ds = \Gamma(a)/\lambda^a$ with $a = n - k > 0$. For $k = 1$: $\lambda\frac{(n-2)!}{(n-1)!} = \frac{\lambda}{n-1}$. For $k = 2$: $\lambda^2\frac{(n-3)!}{(n-1)!} = \frac{\lambda^2}{(n-1)(n-2)}$.

(b) $\hat\lambda = n/S$: $\mathbb{E}\hat\lambda = \frac{n\lambda}{n-1}$, bias $\frac{\lambda}{n-1}$. $\mathrm{Var}\hat\lambda = n^2\Big(\frac{\lambda^2}{(n-1)(n-2)} - \frac{\lambda^2}{(n-1)^2}\Big) = \frac{n^2\lambda^2}{(n-1)^2}\cdot\frac{(n-1) - (n-2)}{n-2} = \frac{n^2\lambda^2}{(n-1)^2(n-2)}$.
$\tilde\lambda = (n-1)/S$: $\mathbb{E}\tilde\lambda = \lambda$, unbiased. $\mathrm{Var}\tilde\lambda = (n-1)^2\Big(\frac{\lambda^2}{(n-1)(n-2)} - \frac{\lambda^2}{(n-1)^2}\Big) = \frac{\lambda^2}{n-2}$.

(c) $I_1 = 1/\lambda^2$, so the bound is $\lambda^2/n$. Efficiency of $\tilde\lambda$: $\frac{\lambda^2/n}{\lambda^2/(n-2)} = \frac{n-2}{n} < 1$, tending to 1.

(d) The score is $s_n = \frac n\lambda - S = -n\big(\bar X - \frac1\lambda\big)$. By Theorem 3.17 an unbiased estimator attains the bound only if the score is $nI_1(\theta)(\hat\theta - \theta)$ for the parameter being estimated. The score is affine in $\bar X$, so the bound is attainable exactly for estimating $\mathbb{E}\bar X = 1/\lambda$ (and affine functions of it), by $\bar X$: $\mathrm{Var}\bar X = \frac{1}{n\lambda^2}$, and the bound for $g(\lambda) = 1/\lambda$ is $\frac{g'(\lambda)^2}{nI_1} = \frac{\lambda^{-4}}{n\lambda^{-2}} = \frac{1}{n\lambda^2}$. Equal. For $\lambda$ itself, attaining would need $s_n = \frac{n}{\lambda^2}(\hat\lambda - \lambda)$, i.e., $\hat\lambda = \lambda + \frac{\lambda^2}{n}(\frac n\lambda - S) = 2\lambda - \frac{\lambda^2S}{n}$, which depends on the unknown $\lambda$ and so is not an estimator.

(e) $n = 10$: bias of $\hat\lambda$ $= \lambda/9 = 0.1111\lambda$; $\mathrm{Var}\hat\lambda = \frac{100}{81\cdot 8}\lambda^2 = 0.1543\lambda^2$; $\mathrm{Var}\tilde\lambda = \lambda^2/8 = 0.1250\lambda^2$; bound $0.1000\lambda^2$; efficiency $0.8$.

**B4.** (a) $p = \sigma(\eta) = 1/(1 + e^{-\eta})$, $dp/d\eta = p(1-p)$. By (3.13), $I_\eta = \big(p(1-p)\big)^2\cdot\frac{1}{p(1-p)} = p(1-p)$.

(b) From the table in §3.10, $A(\eta) = \log(1 + e^\eta)$, $A'(\eta) = \frac{e^\eta}{1 + e^\eta} = p$, and $A''(\eta) = \frac{e^\eta}{(1 + e^\eta)^2} = p(1-p)$. By Theorem 3.26, $I_1(\eta) = A''(\eta) = p(1-p)$. Same answer.

(c) By invariance $\hat\eta = \log\frac{0.3}{0.7} = -0.8473$. Wald: $\mathrm{se} = 1/\sqrt{nI_1(\hat\eta)} = 1/\sqrt{50\cdot 0.21} = 1/\sqrt{10.5} = 0.3086$.

(d) $g(p) = \log\frac{p}{1-p}$ has $g'(p) = \frac1p + \frac{1}{1-p} = \frac{1}{p(1-p)}$. Delta method: $\mathrm{se}(\hat\eta) \approx g'(\hat p)\,\mathrm{se}(\hat p) = \frac{1}{\hat p(1-\hat p)}\sqrt{\frac{\hat p(1-\hat p)}{n}} = \frac{1}{\sqrt{n\hat p(1-\hat p)}}$, the same $0.3086$. This is not a coincidence: (3.13) and (3.19) are the same chain rule.

(e) The interval $\hat\eta \pm 1.96\,\mathrm{se}$ is an interval of real numbers, and mapping its endpoints through $\sigma$ always lands inside $(0, 1)$. The interval $\hat p \pm 1.96\,\mathrm{se}$ can cross 0 or 1 when $T$ is small or close to $n$, and collapses to a point when $T = 0$ (§3.12). A parametrization with an unbounded range usually has a sampling distribution closer to normal, so the Wald approximation is better there.

**B5.** (a) $\mathbb{E}W = P(X_1 = 0) = e^{-\lambda}$. $W$ is Bernoulli, so $\mathrm{Var}W = e^{-\lambda}(1 - e^{-\lambda})$.

(b) Example 4.7 gives $X_1 \mid T = t \sim \mathrm{Binomial}(t, 1/n)$, so $\mathbb{E}[W \mid T = t] = P(X_1 = 0 \mid T = t) = (1 - 1/n)^t$. Hence $\tilde g = (1 - 1/n)^T$. By Theorem 3.23 it is unbiased with variance no larger than $\mathrm{Var}W$.

(c) $T \sim \mathrm{Poisson}(n\lambda)$, so $\mathbb{E}[z^T] = \sum_t z^te^{-n\lambda}\frac{(n\lambda)^t}{t!} = e^{-n\lambda}e^{n\lambda z} = e^{n\lambda(z-1)}$.
- With $z = 1 - \frac1n$: $\mathbb{E}\tilde g = e^{n\lambda(-1/n)} = e^{-\lambda}$. Unbiased.
- With $z = (1 - \frac1n)^2 = 1 - \frac2n + \frac1{n^2}$: $\mathbb{E}\tilde g^2 = e^{n\lambda(-2/n + 1/n^2)} = e^{-2\lambda + \lambda/n}$.
- $\mathrm{Var}\,\tilde g = e^{-2\lambda + \lambda/n} - e^{-2\lambda} = e^{-2\lambda}(e^{\lambda/n} - 1)$.

(d) With $I_1 = 1/\lambda$ and $g'(\lambda) = -e^{-\lambda}$, the bound for unbiased estimators of $g$ is $\frac{g'(\lambda)^2}{nI_1} = \frac{\lambda e^{-2\lambda}}{n}$. Since $e^x - 1 > x$ for every $x > 0$, $\mathrm{Var}\,\tilde g = e^{-2\lambda}(e^{\lambda/n} - 1) > e^{-2\lambda}\frac{\lambda}{n}$: never attained. Why: the score is $s_n = \frac T\lambda - n$, affine in $T$. By Theorem 3.17, only unbiased estimators affine in $T$ can attain the bound, so only $g$ of the form $a + b\lambda$ can be estimated efficiently. $(1 - 1/n)^T$ is not affine in $T$. (It is nonetheless the best unbiased estimator of $e^{-\lambda}$, by Lehmann–Scheffé, because $T$ is complete; the bound is not achievable.)

(e) At $\lambda = 3.2$, $n = 20$: $\mathrm{Var}\,W = 0.03910$; $\mathrm{Var}\,\tilde g = 2.883\times 10^{-4}$; bound $= 2.658\times 10^{-4}$; ratio $1.0844$. Rao–Blackwell cut the variance by a factor of more than 100. The MLE $e^{-\bar X} = e^{-T/n}$ is also a function of $T$, so the same pgf gives its exact moments. With $z = e^{-1/n}$: $\mathbb{E}e^{-\bar X} = e^{n\lambda(e^{-1/n} - 1)} = 0.044099$, so its bias is $3.337\times 10^{-3}$. With $z = e^{-2/n}$: $\mathbb{E}e^{-2\bar X} = e^{n\lambda(e^{-2/n} - 1)} = 0.0022645$, so $\mathrm{MSE} = \mathbb{E}e^{-2\bar X} - 2g\,\mathbb{E}e^{-\bar X} + g^2 = 3.309\times 10^{-4}$. (Part 8's simulation gives $3.57\times 10^{-3}$ and $3.358\times 10^{-4}$, within Monte Carlo error, and asserts both against the exact values.) The MLE is worse than the unbiased $\tilde g$ here: at $n = 20$ its bias (Jensen: $e^{-x}$ is convex, so $\mathbb{E}e^{-\bar X} > e^{-\lambda}$) is not negligible.

**B6.** (a) With $w = n/(n+2)$, $\tilde p = w\hat p + (1 - w)\frac12$. So $\mathbb{E}\tilde p = wp + \frac{1-w}{2}$, bias $(1 - w)(\frac12 - p)$, $\mathrm{Var}\,\tilde p = w^2\frac{p(1-p)}{n}$, and
$$\mathrm{MSE}(\tilde p) = w^2\frac{p(1-p)}{n} + (1 - w)^2\Big(\frac12 - p\Big)^2.$$

(b) Set equal to $\mathrm{MSE}(\hat p) = p(1-p)/n$: $(1 - w^2)\frac{p(1-p)}{n} = (1 - w)^2(\frac12 - p)^2$. For $n = 50$, $w = 50/52$, $1 - w^2 = \frac{52^2 - 50^2}{52^2} = \frac{204}{2704}$ and $(1 - w)^2 = \frac{4}{2704}$. Multiply by 2704:
$$\frac{204}{50}p(1-p) = 4\Big(\frac14 - p + p^2\Big) \iff 4.08p - 4.08p^2 = 1 - 4p + 4p^2 \iff 8.08p^2 - 8.08p + 1 = 0.$$
So $p = \frac{8.08 \pm\sqrt{8.08^2 - 32.32}}{16.16} = \frac{8.08 \pm 5.7416}{16.16}$, giving $p = 0.1447$ and $p = 0.8553$. Between them $\tilde p$ has smaller MSE; outside, $\hat p$ does. (At $p = 0.3$: $0.003942$ against $0.0042$.)

(c) If you have reason to believe $p$ is not extreme, or you care about average performance over a range of $p$ near $\frac12$, report $\tilde p$; it also never gives the absurd estimates 0 or 1. If $p$ could be near 0 or 1 (rare events), $\tilde p$'s pull toward $\frac12$ hurts, and $\hat p$ (with an honest interval) is safer. The answer depends on prior knowledge about $p$, which is exactly what ST-5 makes explicit.

**B7.** (a) By (3.22), $\mathbb{E}T = A'(\eta) = \mu$ and $\mathrm{Var}\,T = A''(\eta)$. Then $\mathbb{E}\bar T = \mu$ and, for iid terms, $\mathrm{Var}\,\bar T = A''(\eta)/n$.

(b) $I_\eta = A''(\eta)$ (Theorem 3.26). The map $\eta \mapsto \mu = A'(\eta)$ has derivative $A''(\eta) > 0$, so it is invertible with $d\eta/d\mu = 1/A''(\eta)$. By (3.13), $I_\mu = (d\eta/d\mu)^2I_\eta = \frac{A''}{(A'')^2} = \frac{1}{A''(\eta)}$.

(c) The bound for unbiased estimators of $\mu$ is $\frac{1}{nI_\mu} = \frac{A''(\eta)}{n} = \mathrm{Var}\,\bar T$. Attained.

(d) Poisson: $A = e^\eta$, $\mu = \lambda$, $A'' = \lambda$; $\mathrm{Var}\bar X = \lambda/n$ and $I_\lambda = 1/\lambda$, consistent with the table in §3.6. Exponential: $A = -\log(-\eta)$, $A' = -1/\eta = 1/\lambda = \mu$, $A'' = 1/\eta^2 = 1/\lambda^2$. So $\bar X$ is efficient for the mean $1/\lambda$, with variance $\frac{1}{n\lambda^2}$, and this is the parameter B3(d) identified. B3 is about $\lambda = 1/\mu$, a nonlinear function of the mean parameter, which is why no unbiased estimator of it is efficient.

**B8.** (a) Write $f_j(x) = \phi(x - m_j)$ with $\phi$ the standard normal density. $\ell(w) = \sum_i\log\big(wf_1(x_i) + (1 - w)f_2(x_i)\big)$. Differentiate twice:
$$\ell'(w) = \sum_i\frac{f_1(x_i) - f_2(x_i)}{wf_1(x_i) + (1-w)f_2(x_i)}, \qquad \ell''(w) = -\sum_i\frac{(f_1(x_i) - f_2(x_i))^2}{(wf_1(x_i) + (1-w)f_2(x_i))^2} \le 0.$$
So $\ell$ is concave in $w$. (Each term is the log of an affine function of $w$.)

(b) $wf_1 + (1-w)f_2$ with $(w, m_1, m_2)$ is the same function of $x$ as $(1 - w)f_2 + wf_1$, which is the density with parameters $(1 - w, m_2, m_1)$. So the likelihoods agree at every data set. Suppose the MLE $\boldsymbol\theta_A = (\hat w, \hat m_1, \hat m_2)$ has $\hat m_1 \ne \hat m_2$. Then $\boldsymbol\theta_B = (1 - \hat w, \hat m_2, \hat m_1)$ is a different point with the same, maximal, likelihood. If $\ell$ were concave, $\ell\big(\frac{\boldsymbol\theta_A + \boldsymbol\theta_B}{2}\big) \ge \frac12\ell(\boldsymbol\theta_A) + \frac12\ell(\boldsymbol\theta_B) = \ell_{\max}$, so the midpoint would also be a maximizer. The midpoint is $(\frac12, \bar m, \bar m)$ with $\bar m = \frac{\hat m_1 + \hat m_2}{2}$: both components equal, i.e., the single normal $\mathcal{N}(\bar m, 1)$. Whenever that single normal fits strictly worse than the mixture MLE, $\ell$ is not concave.

(c) The reference cell at the bottom of `project.py` does it. With the seeded data the two starts converge to $(0.717, 2.052, -1.006)$ and $(0.283, -1.006, 2.052)$, both with negative log-likelihood $748.40$; at their midpoint the negative log-likelihood is $1034.43$, far larger, so $\ell$ is not concave. The cell also checks that second differences of $\ell(w)$ on a grid of 197 values with the means fixed are all negative, as (a) says. The code:

```python
mix_nll = lambda q: -np.log(q[0] * stats.norm.pdf(xm, q[1]) + (1 - q[0]) * stats.norm.pdf(xm, q[2])).sum()
fA = optimize.minimize(mix_nll, [0.5, -2, 3], bounds=b); fB = optimize.minimize(mix_nll, [0.5, 3, -2], bounds=b)
mid = (fA.x + fB.x) / 2          # the single normal N(m_bar, 1)
assert mix_nll(mid) > fA.fun + 50
```

This is why EM (ML-6) finds *a* local maximum and why mixture fits are reported up to relabeling.

## C. Challenge

**C1.** (a) $\ell = -n\log(2\pi\sigma^2) - \frac{1}{2\sigma^2}\sum_i\sum_{j=1}^2(x_{ij} - \mu_i)^2$. For each $i$, $\partial_{\mu_i}\ell = \frac{1}{\sigma^2}\sum_j(x_{ij} - \mu_i) = 0$ gives $\hat\mu_i = \frac{x_{i1} + x_{i2}}{2}$. Then $\partial_{\sigma^2}\ell = -\frac{n}{\sigma^2} + \frac{1}{2\sigma^4}\sum_{i,j}(x_{ij} - \hat\mu_i)^2 = 0$ gives $\hat\sigma^2 = \frac{1}{2n}\sum_{i,j}(x_{ij} - \hat\mu_i)^2$. Since $x_{i1} - \hat\mu_i = \frac{x_{i1} - x_{i2}}{2} = -(x_{i2} - \hat\mu_i)$, each pair contributes $2\cdot\frac{(x_{i1} - x_{i2})^2}{4} = \frac{(x_{i1} - x_{i2})^2}{2}$. So $\hat\sigma^2 = \frac{1}{4n}\sum_i(x_{i1} - x_{i2})^2$.

(b) $D_i = X_{i1} - X_{i2} \sim \mathcal{N}(0, 2\sigma^2)$ are iid, whatever the $\mu_i$. By the weak law, $\frac1n\sum D_i^2 \xrightarrow{p} \mathbb{E}D^2 = 2\sigma^2$. So $\hat\sigma^2 = \frac14\cdot\frac1n\sum D_i^2 \xrightarrow{p} \frac{2\sigma^2}{4} = \frac{\sigma^2}{2}$. The project's simulation with 20,000 pairs and $\sigma^2 = 1$ gives $0.4986$.

(c) Theorem 3.18 and the consistency sketch treat a fixed, finite-dimensional $\theta$ and let $n \to \infty$ with $\frac1n\ell_n(\theta) \to M(\theta)$. Here the parameter $(\mu_1, \dots, \mu_n, \sigma^2)$ grows with $n$. Each $\mu_i$ is informed by two observations no matter how large $n$ is, so $\hat\mu_i$ never converges, and its error eats one of the two degrees of freedom in each pair (the $n - 1$ of Theorem 3.6 with $n = 2$). The divisor should be one per pair, not two.

(d) $\tilde\sigma^2 = \frac{1}{2n}\sum_i(X_{i1} - X_{i2})^2 = 2\hat\sigma^2$, which converges to $\sigma^2$ by the same weak-law step. Equivalently, it is the MLE based only on the differences $D_i \sim \mathcal{N}(0, 2\sigma^2)$, a likelihood with no nuisance parameters. It is also $\frac1n\sum_iS_i^2$ with $S_i^2$ the unbiased within-pair variance.

**C2.** (a) With $u = x - \theta$, $\log p = -\log\pi - \log(1 + u^2)$, so $s = \frac{2u}{1 + u^2}$, which is bounded by 1 in absolute value; the score has every moment even though $X$ has none. Then
$$I_1 = \mathbb{E}s^2 = \int_{-\infty}^\infty\frac{4u^2}{(1 + u^2)^2}\cdot\frac{1}{\pi(1 + u^2)}du = \frac4\pi\int_{-\infty}^\infty\frac{u^2}{(1+u^2)^3}du.$$
Substitute $t = u^2$ on $(0, \infty)$ ($du = dt/(2\sqrt t)$) and double: $\int_{-\infty}^\infty\frac{u^2}{(1+u^2)^3}du = \int_0^\infty\frac{t^{1/2}}{(1+t)^3}dt = B(\tfrac32, \tfrac32) = \frac{\Gamma(3/2)^2}{\Gamma(3)} = \frac{\pi/4}{2} = \frac\pi8$, using $\int_0^\infty\frac{t^{a-1}}{(1+t)^{a+b}}dt = B(a, b)$. So $I_1 = \frac4\pi\cdot\frac\pi8 = \frac12$. The project confirms $0.5000$ by quadrature.

(b) F5 §3.11 shows, via characteristic functions, that the average of $n$ iid standard Cauchy variables is again standard Cauchy. Shifting by $\theta$ changes nothing, so $\bar X \sim$ Cauchy$(\theta, 1)$ for every $n$: averaging 100 observations is no better than looking at one. Its interquartile range is 2 for every $n$ (the quartiles of a standard Cauchy are $\pm 1$).

(c) The reference cell computes the MLE by a grid search over median $\pm 1$ followed by five vectorized Newton steps (the median is consistent, and starting there picks out the right root among possibly several, §3.12). Results with $n = 100$, 2000 replications: $n\mathrm{Var}(\text{MLE}) = 2.029$, close to $1/I_1 = 2$; $n\mathrm{Var}(\text{median}) = 2.514$, close to its asymptotic value $\pi^2/4 = 2.467$ (the median's asymptotic variance is $1/(4f(\theta)^2n)$ with $f(\theta) = 1/\pi$); and the interquartile range of $\bar X$ is $2.02$, the same as for a single observation. The MLE beats the median by the factor $2/2.467$, as Theorem 3.19 predicts; the mean is useless.

**C3.** (a) Unbiasedness: $\int\hat\theta_i(\mathbf{x})\,p_n(\mathbf{x};\theta)\,d\mathbf{x} = \theta_i$ for each $i$. Differentiate in $\theta_j$ under the integral: $\int\hat\theta_i\,\partial_{\theta_j}p_n\,d\mathbf{x} = \delta_{ij}$, i.e., $\mathbb{E}[\hat\theta_i\,(s_n)_j] = \delta_{ij}$. In matrix form $\mathbb{E}[\hat\theta s_n^\top] = I_k$. Since $\mathbb{E}s_n = \mathbf{0}$, also $\mathrm{Cov}(\hat\theta, s_n) = I_k$.

(b) Let $U = \mathbf{a}^\top\hat\theta$ and $V = \mathbf{b}^\top s_n$. Then $\mathrm{Cov}(U, V) = \mathbf{a}^\top\mathrm{Cov}(\hat\theta, s_n)\mathbf{b} = \mathbf{a}^\top\mathbf{b}$, $\mathrm{Var}(U) = \mathbf{a}^\top\mathrm{Cov}(\hat\theta)\mathbf{a}$, and $\mathrm{Var}(V) = \mathbf{b}^\top(nI_1)\mathbf{b}$ (Theorems 3.12, 3.13). Cauchy–Schwarz (F5 Theorem 3.17):
$$(\mathbf{a}^\top\mathbf{b})^2 \le \mathbf{a}^\top\mathrm{Cov}(\hat\theta)\mathbf{a}\cdot\mathbf{b}^\top(nI_1)\mathbf{b}.$$
Choose $\mathbf{b} = (nI_1)^{-1}\mathbf{a}$ (the choice that makes the inequality tightest). Then $\mathbf{a}^\top\mathbf{b} = \mathbf{a}^\top(nI_1)^{-1}\mathbf{a}$ and $\mathbf{b}^\top(nI_1)\mathbf{b} = \mathbf{a}^\top(nI_1)^{-1}\mathbf{a}$. The inequality becomes $\big(\mathbf{a}^\top(nI_1)^{-1}\mathbf{a}\big)^2 \le \mathbf{a}^\top\mathrm{Cov}(\hat\theta)\mathbf{a}\cdot\mathbf{a}^\top(nI_1)^{-1}\mathbf{a}$. For $\mathbf{a} \ne \mathbf{0}$ the factor $\mathbf{a}^\top(nI_1)^{-1}\mathbf{a}$ is positive ($(nI_1)^{-1}$ is positive definite by (R4)); divide by it to get $\mathbf{a}^\top\mathrm{Cov}(\hat\theta)\mathbf{a} \ge \mathbf{a}^\top(nI_1)^{-1}\mathbf{a}$. Since this holds for every $\mathbf{a}$, $\mathrm{Cov}(\hat\theta) - (nI_1)^{-1} \succeq 0$. ∎

(c) $\langle\mathbf{u}, \mathbf{v}\rangle_M = \mathbf{u}^\top M\mathbf{v}$ is an inner product because $M$ is symmetric positive definite (LA-8). Apply Cauchy–Schwarz (LA-8 Theorem 3.4) to $\mathbf{u} = M^{-1}\mathbf{e}_j$ and $\mathbf{v} = \mathbf{e}_j$: $\langle\mathbf{u}, \mathbf{v}\rangle_M = \mathbf{e}_j^\top M^{-1}M\mathbf{e}_j = 1$, $\langle\mathbf{u}, \mathbf{u}\rangle_M = \mathbf{e}_j^\top M^{-1}\mathbf{e}_j = [M^{-1}]_{jj}$, $\langle\mathbf{v}, \mathbf{v}\rangle_M = M_{jj}$. So $1 \le [M^{-1}]_{jj}M_{jj}$, i.e., $[M^{-1}]_{jj} \ge 1/M_{jj}$. ∎ Interpretation: $1/[nI_1]_{jj}$ is the bound for $\theta_j$ when all other parameters are known; $[(nI_1)^{-1}]_{jj}$ is the bound when they must be estimated too. For the Gamma model at $\alpha = 2$, $I_{\alpha\alpha} = \psi'(2) = \pi^2/6 - 1 = 0.6449$. With $\beta$ known the asymptotic variance of $\hat\alpha$ would be $1/\psi'(2) = 1.5505$ per observation; with $\beta$ unknown it is $[I_1^{-1}]_{\alpha\alpha} = 6.900$ (Example 4.5). Not knowing the rate makes the shape more than four times harder to estimate, because the two parameters trade off along the ridge of the likelihood ($\alpha/\beta$ fixed by $\bar x$).
