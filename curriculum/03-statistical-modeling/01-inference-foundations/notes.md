# ST-1: Inference Foundations: Models, Estimators, Bias/Variance/MSE, Consistency, Sufficiency, MLE, Method of Moments, Fisher Information, Cramér–Rao

**Time:** ~12 hours · **Prerequisites:** F5 (and LA-14 for the multiparameter calculus) · **Unlocks:** ST-2, ML-1; through §3.10, also ST-4 and ST-5

> **Conventions in this module.** Data $X_1, \dots, X_n$ are iid from a density or PMF $p(x;\theta)$ with $\theta \in \Theta \subseteq \mathbb{R}^k$; lowercase $x_i$ are observed values; $\theta_0$ is the true value. An **estimator** $\hat\theta_n = T(X_1, \dots, X_n)$ is a function of the data alone, hence a random variable. $\mathbb{E}_\theta, \mathrm{Var}_\theta$ are computed with parameter $\theta$. Logs are natural. Log-likelihood $\ell_n(\theta) = \sum_i \log p(X_i;\theta)$; per-observation score $s(\theta;x) = \nabla_\theta \log p(x;\theta)$ (a column); total score $s_n = \nabla\ell_n$. Fisher information per observation $I_1(\theta) = \mathbb{E}_\theta[ss^\top] = -\mathbb{E}_\theta[\nabla^2\log p]$; for $n$ observations $nI_1$. Bias $b(\theta) = \mathbb{E}_\theta\hat\theta - \theta$; $\mathrm{MSE} = \mathbb{E}_\theta(\hat\theta - \theta)^2$. $\xrightarrow{p}$, $\xrightarrow{d}$ are F5 Definition 3.33. Parametrizations: $\mathrm{Exponential}(\lambda)$ has **rate** $\lambda$, density $\lambda e^{-\lambda x}$; $\mathrm{Gamma}(\alpha, \beta)$ has shape $\alpha$ and **rate** $\beta$, density $\beta^\alpha x^{\alpha-1}e^{-\beta x}/\Gamma(\alpha)$ (F5's $\mathrm{Gamma}(a, \lambda)$); $\mathcal{N}(\mu, \sigma^2)$ is parametrized by the **variance**. $S^2 = \frac{1}{n-1}\sum(X_i - \bar X)^2$ and $\hat\sigma^2 = \frac1n\sum(X_i - \bar X)^2$. $\psi = \Gamma'/\Gamma$ is the digamma function. $A \succeq B$ means $A - B$ is positive semidefinite (LA-8).

## 0. Where this sits

You need F5 fluently: independence as factorization (Theorem 3.4), the weak law (Theorem 3.34), the CLT (3.24), Adam's and Eve's laws (Theorems 3.26, 3.28), Cauchy–Schwarz for covariances (Theorem 3.17), and KL divergence with Gibbs' inequality (Theorem 3.45). LA-14 supplies gradients, Hessians, and §3.9's Gaussian log-likelihood gradients (3.38)–(3.39); only §3.5(f) and the matrix statements need it.

Probability went from a known model to data. Inference goes backwards. This module builds the vocabulary every later statistics and ML module assumes: estimator, standard error, likelihood, score, Fisher information, sufficiency. ST-2 turns sampling distributions into intervals and tests. ML-1 reuses the bias–variance decomposition and recognizes log-loss minimization as maximum likelihood. ST-4 (GLMs) and ST-5 (conjugate priors) are built on the exponential families of §3.10.

## 1. The problem that forced this idea

In *Theoria motus* (1809) Gauss justified least squares, the method behind his recovery of the asteroid Ceres, by asking which error distribution makes the arithmetic mean of repeated measurements the most probable value. The answer was the normal law. A year later Laplace connected this to his central limit theorem: errors built from many small independent causes are approximately normal. So by 1810 there was a method and a justification, but for one model only, and no account of how accurate the answer was.

In 1894 Karl Pearson was handed W. F. R. Weldon's measurements of 1,000 Naples crabs. The histogram of one body ratio was lopsided, and Pearson guessed it was a mixture of two normal populations, with five unknown parameters. He estimated them by setting the first five sample moments equal to the model's moments, which led to a ninth-degree polynomial he solved by hand. That is the **method of moments**. It gives *an* answer. It says nothing about whether the answer is good.

R. A. Fisher attacked that gap. As a Cambridge undergraduate he proposed (1912) choosing the parameters that make the observed sample most probable. In "On the mathematical foundations of theoretical statistics" (1922) he defined **consistency**, **efficiency**, and **sufficiency**, called his criterion **maximum likelihood**, and showed that Pearson's moments can waste much of the information in the data. In 1925 he defined the **amount of information** and claimed that maximum likelihood asymptotically extracts all of it.

The optimality theory followed. C. R. Rao (1945) and Harald Cramér (*Mathematical Methods of Statistics*, 1946) proved that no unbiased estimator can have variance below one over the information; Fréchet (1943) and Darmois (1945) found the same bound independently. Rao's 1945 paper and Blackwell (1947) showed how to improve any estimator by conditioning on a sufficient statistic. Wald (1949) gave conditions under which maximum likelihood is consistent, and Neyman and Scott (1948) showed one way it can fail.

The forcing question has three parts. Given data from a distribution with unknown parameters: how do you estimate them, how good is the estimate, and is there a best possible one? Answering the second part needs a conceptual move that probability alone does not supply. **The estimate is a number; the estimator is a random variable**, and its distribution over hypothetical repetitions of the experiment is what we study.

## 2. Intuition

### An estimator is a recipe, and recipes have distributions

Take ten draws from $\mathcal{N}(\mu, 1)$ with $\mu$ unknown. The recipe "average them" gives a number, say 4.83. Another ten draws would give 5.12. Repeat the experiment forever and the outputs pile up into a histogram: the **sampling distribution** of the estimator. Every statement about how good an estimator is, is a statement about that histogram, *not* about the one number you computed.

```
               truth
                 |
   A:        .:::|:::.            centered on the truth, wide
   B:            |   .:::.        off-center, narrow
                 |
```

Two numbers summarize the histogram relative to the truth. **Bias** is how far its center sits from the truth. **Variance** is how spread out it is. The overall quality is the average squared distance from the truth, the **mean squared error**. Think of the center as a point and the spread as fluctuation around that point: squared distance to the truth should then split like Pythagoras into (center offset)$^2$ + (spread)$^2$. That guess is right, and it means estimator B can beat estimator A: a little bias can buy a lot of variance.

### The likelihood is surprise, read backwards

A coin lands heads 7 times in 10. For each candidate bias $p$, ask how probable this exact outcome would be: $L(p) = p^7(1-p)^3$. A coin with $p = 0.1$ would be astonished by this data; one with $p = 0.7$ would not. The **likelihood** is the probability of the data with the data held fixed and the parameter varied, and the **maximum likelihood estimate** is the least-surprised parameter: here $p = 0.7$ (set the derivative of $7\log p + 3\log(1-p)$ to zero).

### Curvature is information

Now take $n$ flips with $T$ heads. The log-likelihood $\ell(p) = T\log p + (n-T)\log(1-p)$ peaks at $\hat p = T/n$, and its second derivative there is $-n/(\hat p(1-\hat p))$. Draw it for $n = 10$ and $n = 1000$ with the same proportion of heads:

```
 log L                                   log L
   |       .--.        n = 10              |         .        n = 1000
   |    .-'    '-.                         |        / \
   |  .'          '.                       |       /   \
   |_/______________\___ p                 |______/_____\_____ p
            ^ 0.7                                   ^ 0.7
```

A flat peak means many parameters explain the data about equally well; a sharp peak means the data pin the parameter down. Now compare with something you already know: $\mathrm{Var}(\hat p) = p(1-p)/n$, which is exactly **one over the curvature**. Guess the general rule: an estimator's variance is at least $1/(\text{expected curvature of the log-likelihood})$, and maximum likelihood achieves it for large $n$. The expected curvature per observation is the **Fisher information**. The first half of the guess is the Cramér–Rao bound (§3.7), the second half is the asymptotic efficiency of the MLE (§3.8).

So before §3 you should be able to guess three results: MSE = bias² + variance; the MLE's variance is about $1/(nI_1)$; and no unbiased estimator does better than $1/(nI_1)$.

## 3. Formal development

### 3.1 Models, estimators, sampling distributions

**Definition 3.1.** A **statistical model** is a family of distributions $\{p(\cdot;\theta) : \theta \in \Theta\}$. It is **parametric** if $\Theta \subseteq \mathbb{R}^k$ for some finite $k$, and **nonparametric** otherwise (for example "all distributions on $\mathbb{R}$ with finite variance"). We assume $X_1, \dots, X_n$ iid from $p(\cdot;\theta_0)$ for some unknown $\theta_0 \in \Theta$: the model is **correctly specified** (until §3.12). By F5 Theorem 3.4 the joint density is $p_n(\mathbf{x};\theta) = \prod_i p(x_i;\theta)$. A **statistic** is any function of the data not involving $\theta$; an **estimator** of $\theta$ is a statistic used to guess $\theta$. Its distribution under $p(\cdot;\theta)$ is its **sampling distribution**, and the **standard error** $\mathrm{se}(\hat\theta)$ is the standard deviation of that distribution.

### 3.2 Bias, variance, MSE, consistency

**Definition 3.2.** For an estimator $\hat\theta$ of a scalar $\theta$:
$$b(\theta) = \mathbb{E}_\theta[\hat\theta] - \theta, \qquad \mathrm{Var}_\theta(\hat\theta), \qquad \mathrm{MSE}(\theta) = \mathbb{E}_\theta\big[(\hat\theta - \theta)^2\big]. \tag{3.1}$$
All three are functions of $\theta$: an estimator can be good where the truth is one value and bad where it is another.

**Theorem 3.3 (Bias–variance decomposition).** $\mathrm{MSE}(\theta) = b(\theta)^2 + \mathrm{Var}_\theta(\hat\theta)$.

▶ **Try it first.** Add and subtract $m = \mathbb{E}_\theta\hat\theta$ inside the square.

**Derivation.**
1. $\hat\theta - \theta = (\hat\theta - m) + (m - \theta)$, and $m - \theta = b$ is a constant.
2. Square: $(\hat\theta - m)^2 + 2b(\hat\theta - m) + b^2$.
3. Take expectations. $\mathbb{E}[\hat\theta - m] = 0$ kills the middle term:
$$\mathrm{MSE} = \mathrm{Var}(\hat\theta) + b^2. \tag{3.2}$$ ∎

■ **What just happened.** Pythagoras, in F5's inner product $\langle U, V\rangle = \mathbb{E}[UV]$: the fluctuation $\hat\theta - m$ has mean zero, so it is orthogonal to every constant, including the offset $b$.

**Definition 3.4.** $\hat\theta_n$ is **consistent** if $\hat\theta_n \xrightarrow{p} \theta$ for every $\theta \in \Theta$.

**Theorem 3.5.** If $\mathrm{MSE}_n(\theta) \to 0$, then $\hat\theta_n$ is consistent. In particular, bias $\to 0$ and variance $\to 0$ suffice.

▶ **Try it first.** Apply Markov's inequality (F5 Theorem 3.30) to the nonnegative variable $(\hat\theta_n - \theta)^2$.

**Derivation.**
1. Fix $\epsilon > 0$. The events $\{|\hat\theta_n - \theta| \ge \epsilon\}$ and $\{(\hat\theta_n - \theta)^2 \ge \epsilon^2\}$ are the same.
2. Markov: $P\big((\hat\theta_n - \theta)^2 \ge \epsilon^2\big) \le \mathbb{E}(\hat\theta_n - \theta)^2/\epsilon^2 = \mathrm{MSE}_n/\epsilon^2 \to 0$.
3. By (3.2), $\mathrm{MSE}_n = b_n^2 + \mathrm{Var}_n$, so $b_n \to 0$ and $\mathrm{Var}_n \to 0$ give $\mathrm{MSE}_n \to 0$. ∎

■ **What just happened.** This is Chebyshev's inequality centered at the truth instead of the mean. It is how the weak law (F5 Theorem 3.34) was proved, and it is the first thing to try for any new estimator.

### 3.3 The running example: sample mean and sample variance

$\mathbb{E}\bar X = \mu$ and $\mathrm{Var}\,\bar X = \sigma^2/n$ (F5 §3.11), so $\bar X$ is unbiased with $\mathrm{MSE} = \sigma^2/n \to 0$: consistent by Theorem 3.5.

**Theorem 3.6.** If the $X_i$ are iid with variance $\sigma^2 < \infty$, then $\mathbb{E}[S^2] = \sigma^2$.

▶ **Try it first.** Write $X_i - \bar X = (X_i - \mu) - (\bar X - \mu)$ and expand $\sum(X_i - \bar X)^2$.

**Derivation.**
1. $\sum(X_i - \bar X)^2 = \sum(X_i - \mu)^2 - 2(\bar X - \mu)\sum(X_i - \mu) + n(\bar X - \mu)^2$.
2. $\sum(X_i - \mu) = n(\bar X - \mu)$, so the middle term is $-2n(\bar X - \mu)^2$ and, for any constant $\mu$,
$$\sum_i(X_i - \bar X)^2 = \sum_i(X_i - \mu)^2 - n(\bar X - \mu)^2. \tag{3.3}$$
3. Take expectations with $\mu$ the true mean: $n\sigma^2 - n\cdot\sigma^2/n = (n-1)\sigma^2$. Divide by $n - 1$. ∎

■ **What just happened.** $\bar X$ is the point closest to the data in squared distance, so deviations from it are systematically smaller than deviations from $\mu$. The shortfall is exactly one $\sigma^2$: one degree of freedom spent estimating $\mu$.

**Theorem 3.7 (Normal samples).** If $X_i \overset{\text{iid}}{\sim} \mathcal{N}(\mu, \sigma^2)$, then $\bar X$ and $S^2$ are independent, $(n-1)S^2/\sigma^2 \sim \chi^2_{n-1}$, and
$$\mathrm{Var}(S^2) = \frac{2\sigma^4}{n-1}. \tag{3.4}$$

▶ **Try it first.** Rotate the data vector $\mathbf{X} = (X_1, \dots, X_n)^\top$ by an orthogonal matrix whose first row is $\mathbf{1}^\top/\sqrt n$. What does F5 Theorem 3.22 say about the rotated vector?

**Derivation.**
1. Complete $\mathbf{1}/\sqrt n$ to an orthonormal basis (Gram–Schmidt, LA-5) and let $Q$ have these as rows. Put $\mathbf{Y} = Q\mathbf{X}$. Since $\mathbf{X} \sim \mathcal{N}(\mu\mathbf{1}, \sigma^2 I)$, F5 Theorem 3.22 gives $\mathbf{Y} \sim \mathcal{N}(\mu Q\mathbf{1}, \sigma^2 QQ^\top) = \mathcal{N}(\mu Q\mathbf{1}, \sigma^2 I)$: the $Y_i$ are independent normals with variance $\sigma^2$.
2. $(Q\mathbf{1})_1 = \mathbf{1}^\top\mathbf{1}/\sqrt n = \sqrt n$, and for $i \ge 2$ row $i$ is orthogonal to $\mathbf{1}$, so $(Q\mathbf{1})_i = 0$. Hence $Y_1 = \sqrt n\,\bar X$, and $Y_2, \dots, Y_n \sim \mathcal{N}(0, \sigma^2)$.
3. Orthogonal maps preserve length: $\sum Y_i^2 = \sum X_i^2$. So $\sum_{i \ge 2}Y_i^2 = \sum X_i^2 - n\bar X^2 = \sum(X_i - \bar X)^2$ by (3.3) with the constant $0$.
4. So $(n-1)S^2/\sigma^2 = \sum_{i=2}^n (Y_i/\sigma)^2$ is a sum of $n-1$ independent squared standard normals, $\chi^2_{n-1}$ by definition. It depends only on $Y_2, \dots, Y_n$, while $\bar X$ depends only on $Y_1$: independent.
5. For $Z \sim \mathcal{N}(0,1)$, $\mathbb{E}Z^4 = 3$, so $\mathrm{Var}(Z^2) = 2$ and $\mathrm{Var}(\chi^2_{n-1}) = 2(n-1)$. Then $\mathrm{Var}(S^2) = \frac{\sigma^4}{(n-1)^2}\cdot 2(n-1)$. ∎

■ **What just happened.** One rotation splits the data into the mean direction and $n-1$ directions orthogonal to it. The residuals live in that $(n-1)$-dimensional subspace, which is what "$n-1$ degrees of freedom" means. ST-3 does the same with $p$ regression directions instead of one; ST-2 uses the independence of $\bar X$ and $S^2$ to build the $t$ statistic.

**The MLE of $\sigma^2$ has lower MSE.** §3.5 shows the MLE is $\hat\sigma^2 = \frac{n-1}{n}S^2$. From Theorem 3.6 and (3.4):
$$b = -\frac{\sigma^2}{n}, \qquad \mathrm{Var}(\hat\sigma^2) = \frac{2(n-1)\sigma^4}{n^2}, \qquad \mathrm{MSE}(\hat\sigma^2) = \frac{(2n-1)\sigma^4}{n^2} < \frac{2\sigma^4}{n-1} = \mathrm{MSE}(S^2), \tag{3.5}$$
because $(2n-1)(n-1) = 2n^2 - 3n + 1 < 2n^2$. Shrinking $S^2$ by $(n-1)/n$ costs bias² $\sigma^4/n^2$ and saves more than that in variance. This is the bias–variance trade-off's first appearance; Exercise B1 finds the best divisor ($n+1$), and Example 4.1 has numbers. Unbiasedness is a convention, not a law of nature.

### 3.4 The method of moments

**Definition 3.8.** Let $m_j(\theta) = \mathbb{E}_\theta[X^j]$ and $\hat m_j = \frac1n\sum_i X_i^j$. For $\theta \in \mathbb{R}^k$, the **method-of-moments (MoM) estimator** solves
$$m_j(\hat\theta) = \hat m_j, \qquad j = 1, \dots, k. \tag{3.6}$$
(Matching central moments instead is equivalent.)

- **Normal.** $m_1 = \mu$, $m_2 = \sigma^2 + \mu^2$. So $\hat\mu = \bar X$ and $\hat\sigma^2 = \hat m_2 - \bar X^2 = \frac1n\sum(X_i - \bar X)^2$ (divisor $n$; use (3.3) with the constant 0).
- **Gamma$(\alpha, \beta)$.** Mean $\alpha/\beta$, variance $\alpha/\beta^2$ (F5 (3.9)). Setting $\alpha/\beta = \bar X$ and $\alpha/\beta^2 = \hat\sigma^2$ and dividing,
$$\hat\beta_{\text{MoM}} = \bar X/\hat\sigma^2, \qquad \hat\alpha_{\text{MoM}} = \bar X^2/\hat\sigma^2. \tag{3.7}$$
- **Uniform$(0, \theta)$.** $m_1 = \theta/2$, so $\hat\theta = 2\bar X$: unbiased, with variance $4\cdot\frac{\theta^2}{12}\cdot\frac1n = \frac{\theta^2}{3n}$. Its defect: $2\bar X$ can be smaller than the largest observation, an estimate the data have already ruled out.

**Theorem 3.9 (Consistency of MoM).** Suppose $\mathbb{E}[X^{2k}] < \infty$ and the moment map $\theta \mapsto (m_1(\theta), \dots, m_k(\theta))$ has an inverse $h$ that is continuous at the true moments. Then $\hat\theta = h(\hat m_1, \dots, \hat m_k) \xrightarrow{p} \theta_0$.

*Proof.* By the weak law applied to $X^j$ (which has finite variance because $\mathbb{E}X^{2j} < \infty$), $\hat m_j \xrightarrow{p} m_j(\theta_0)$ for each $j$. **Continuous mapping:** if $\mathbf{Y}_n \xrightarrow{p} \mathbf{c}$ and $h$ is continuous at $\mathbf{c}$, then $h(\mathbf{Y}_n) \xrightarrow{p} h(\mathbf{c})$, because for any $\epsilon$ there is a $\delta$ with $\|\mathbf{y} - \mathbf{c}\| < \delta \Rightarrow \|h(\mathbf{y}) - h(\mathbf{c})\| < \epsilon$, so $P(\|h(\mathbf{Y}_n) - h(\mathbf{c})\| \ge \epsilon) \le P(\|\mathbf{Y}_n - \mathbf{c}\| \ge \delta) \le \sum_j P(|Y_{n,j} - c_j| \ge \delta/\sqrt k) \to 0$. ∎

MoM is consistent for the cheapest possible reason: averages converge. It is often inefficient. For the Gamma shape at $\alpha = 2$, its asymptotic variance is $2\alpha(\alpha+1)/n = 12/n$ (delta method on $(\hat m_1, \hat m_2)$; the project derives it) against the MLE's $6.90/n$ (Example 4.5): an asymptotic efficiency of $0.575$. That was Fisher's 1922 complaint.

### 3.5 Maximum likelihood

**Definition 3.10.** Given observed data, the **likelihood** is $L(\theta) = \prod_i p(x_i;\theta)$, the log-likelihood is $\ell(\theta) = \sum_i\log p(x_i;\theta)$, and the **maximum likelihood estimator** is $\hat\theta = \arg\max_\theta L(\theta) = \arg\max_\theta\ell(\theta)$. $L$ is a function of $\theta$, not a density in $\theta$; it need not integrate to 1.

**Theorem 3.11 (Invariance).** If $\hat\theta$ is the MLE of $\theta$, then $g(\hat\theta)$ is the MLE of $\eta = g(\theta)$ for any function $g$, where for non-injective $g$ the likelihood of $\eta$ is the **profile** $L^*(\eta) = \sup\{L(\theta) : g(\theta) = \eta\}$.

▶ **Try it first.** Do the case of a bijection $g$ first.

**Derivation.**
1. *Bijection.* $L^*(\eta) = L(g^{-1}(\eta))$. As $\eta$ ranges over $g(\Theta)$, $g^{-1}(\eta)$ ranges over all of $\Theta$, so $L^*$ is maximized where $g^{-1}(\eta) = \hat\theta$, i.e., at $\eta = g(\hat\theta)$.
2. *General.* For every $\eta$, $L^*(\eta) \le \sup_\Theta L = L(\hat\theta)$.
3. $\hat\theta$ belongs to $\{\theta : g(\theta) = g(\hat\theta)\}$, so $L^*(g(\hat\theta)) \ge L(\hat\theta)$. With step 2, $L^*(g(\hat\theta)) = L(\hat\theta) \ge L^*(\eta)$ for all $\eta$. ∎

■ **What just happened.** Maximum likelihood does not care how you label the parameter: the MLE of $\sigma$ is $\sqrt{\hat\sigma^2}$, of $e^{-\lambda}$ is $e^{-\hat\lambda}$. Unbiasedness has no such property: $S^2$ is unbiased for $\sigma^2$, but $S$ is biased for $\sigma$ (Jensen, F5 Theorem 3.32).

**The standard MLEs.** Write $T = \sum x_i$.

(a) **Bernoulli$(p)$.** $\ell(p) = T\log p + (n-T)\log(1-p)$. $\ell'(p) = T/p - (n-T)/(1-p) = 0$ gives $\hat p = T/n$, and $\ell''(p) = -T/p^2 - (n-T)/(1-p)^2 < 0$, so it is the maximum when $0 < T < n$ (the other cases are in §3.12).

(b) **Poisson$(\lambda)$.** $\ell(\lambda) = T\log\lambda - n\lambda - \sum\log x_i!$. $\ell' = T/\lambda - n = 0$ gives $\hat\lambda = \bar x$; $\ell'' = -T/\lambda^2 < 0$.

(c) **Exponential$(\lambda)$.** $\ell(\lambda) = n\log\lambda - \lambda T$. $\ell' = n/\lambda - T = 0$ gives $\hat\lambda = n/T = 1/\bar x$; $\ell'' = -n/\lambda^2 < 0$.

(d) **Normal$(\mu, v)$, $v = \sigma^2$, both unknown.**
$$\ell(\mu, v) = -\frac n2\log(2\pi v) - \frac{1}{2v}\sum_i(x_i - \mu)^2. \tag{3.8}$$
The gradient (LA-14 (3.38)–(3.39) with $p = 1$) is $\partial_\mu\ell = \sum(x_i - \mu)/v$, $\partial_v\ell = -\frac{n}{2v} + \frac{1}{2v^2}\sum(x_i - \mu)^2$; its zero is $\hat\mu = \bar x$, $\hat v = \frac1n\sum(x_i - \bar x)^2$. The Hessian entries are $\partial^2_\mu\ell = -n/v$, $\partial_\mu\partial_v\ell = -\sum(x_i - \mu)/v^2$, $\partial^2_v\ell = \frac{n}{2v^2} - \frac{1}{v^3}\sum(x_i - \mu)^2$. At the MLE the cross term is zero and $\partial_v^2\ell = \frac{n}{2\hat v^2} - \frac{n\hat v}{\hat v^3} = -\frac{n}{2\hat v^2}$. The Hessian $\mathrm{diag}(-n/\hat v, -n/(2\hat v^2))$ is negative definite, so this is a strict local maximum (LA-8 Theorem 3.19); it is the only critical point and $\ell \to -\infty$ at the edges, so it is global.

(e) **Uniform$(0, \theta)$: where calculus fails.** $p(x;\theta) = 1/\theta$ on $[0, \theta]$. If $\theta < \max_i x_i$, some observation is impossible and $L(\theta) = 0$. For $\theta \ge \max_i x_i$, $L(\theta) = \theta^{-n}$, strictly decreasing. So $\hat\theta = X_{(n)} = \max_i X_i$. The derivative $\ell'(\theta) = -n/\theta$ never vanishes; the maximum sits at a discontinuity, where setting a derivative to zero finds nothing. Example 4.4 shows the regularity behind §3.6–§3.8 fails here too.

(f) **Multivariate normal $\mathcal{N}(\boldsymbol\mu, \Sigma)$ on $\mathbb{R}^d$.** LA-14 §3.9 derived, using $\nabla\log\det\Sigma = \Sigma^{-1}$ (LA-14 (3.17)) and $d(\Sigma^{-1}) = -\Sigma^{-1}(d\Sigma)\Sigma^{-1}$ (LA-14 (3.6)), that $\nabla_{\boldsymbol\mu}\ell = \Sigma^{-1}\sum_i(\mathbf{x}_i - \boldsymbol\mu)$ and $\nabla_\Sigma\ell = \frac12(\Sigma^{-1}S_{\boldsymbol\mu}\Sigma^{-1} - n\Sigma^{-1})$ with $S_{\boldsymbol\mu} = \sum_i(\mathbf{x}_i - \boldsymbol\mu)(\mathbf{x}_i - \boldsymbol\mu)^\top$. Setting both to zero and multiplying by $\Sigma$:
$$\hat{\boldsymbol\mu} = \bar{\mathbf{x}}, \qquad \hat\Sigma = \frac1n\sum_i(\mathbf{x}_i - \bar{\mathbf{x}})(\mathbf{x}_i - \bar{\mathbf{x}})^\top. \tag{3.9}$$
The divisor is $n$, biased by the factor $(n-1)/n$ entrywise; (d) is the case $d = 1$.

(g) **Gamma$(\alpha, \beta)$: no closed form.** $\ell = n\alpha\log\beta - n\log\Gamma(\alpha) + (\alpha - 1)\sum\log x_i - \beta\sum x_i$.
1. $\partial_\beta\ell = n\alpha/\beta - \sum x_i = 0$ gives $\beta = \alpha/\bar x$.
2. $\partial_\alpha\ell = n\log\beta - n\psi(\alpha) + \sum\log x_i = 0$. Substitute step 1:
$$\log\alpha - \psi(\alpha) = \log\bar x - \frac1n\sum_i\log x_i =: s. \tag{3.10}$$
3. The right side is $\ge 0$ by Jensen ($\log$ is concave), and $> 0$ unless all $x_i$ are equal. The left side is strictly decreasing (the standard inequality $\psi'(\alpha) > 1/\alpha$) from $+\infty$ as $\alpha \to 0$ to $0$ as $\alpha \to \infty$ (because $\psi(\alpha) = \log\alpha - \frac{1}{2\alpha} + O(\alpha^{-2})$). So (3.10) has exactly one root, found numerically. Example 4.5 and project Part 4 do it.

### 3.6 The score and Fisher information

We need **regularity conditions (R)**: (R1) the support $\{x : p(x;\theta) > 0\}$ does not depend on $\theta$; (R2) $p$ is twice differentiable in $\theta$ and $\int p\,dx$ may be differentiated twice under the integral sign; (R3) $\theta$ is an interior point of $\Theta$; (R4) $I_1(\theta)$ is finite and positive definite. Every family in this module satisfies (R) except Uniform$(0, \theta)$, which fails (R1).

**Theorem 3.12 (Score identities).** Under (R):
$$\mathbb{E}_\theta[s(\theta;X)] = \mathbf{0}, \qquad I_1(\theta) = \mathrm{Cov}_\theta(s) = \mathbb{E}_\theta[ss^\top] = -\mathbb{E}_\theta\big[\nabla^2_\theta\log p(X;\theta)\big]. \tag{3.11}$$

▶ **Try it first.** Differentiate $\int p(x;\theta)\,dx = 1$ in $\theta$ and use $\partial_\theta p = p\,\partial_\theta\log p$. Then differentiate once more.

**Derivation** (scalar $\theta$; the matrix version is the same entry by entry).
1. $\int p(x;\theta)\,dx = 1$ for every $\theta$. Differentiate, using (R1)–(R2): $\int\partial_\theta p\,dx = 0$.
2. $\partial_\theta p = p\,\partial_\theta\log p = p\,s$. So $\int s\,p\,dx = \mathbb{E}_\theta[s] = 0$.
3. Differentiate $\int s(\theta;x)\,p(x;\theta)\,dx = 0$ again: $\int(\partial_\theta s)\,p\,dx + \int s\,\partial_\theta p\,dx = 0$.
4. The second integral is $\int s^2 p\,dx = \mathbb{E}[s^2]$; the first is $\mathbb{E}[\partial^2_\theta\log p]$. So $\mathbb{E}[s^2] = -\mathbb{E}[\partial^2_\theta\log p]$.
5. Since $\mathbb{E}[s] = 0$, $\mathbb{E}[s^2] = \mathrm{Var}(s)$. ∎

■ **What just happened.** Two faces of one number: how much the data jostle the slope of the log-likelihood equals the average curvature of the log-likelihood. The curvature form is usually easier to compute. Step 1 is where (R1) matters: if the limits of integration move with $\theta$, differentiating the integral produces boundary terms (Leibniz's rule), and Example 4.4 shows the identity then fails.

**The standard families**, per observation, from $-\mathbb{E}[\partial^2_\theta\log p]$:

| Family | $\log p(x;\theta)$ | $-\partial^2_\theta\log p$ | $I_1(\theta)$ |
|---|---|---|---|
| Bernoulli$(p)$ | $x\log p + (1-x)\log(1-p)$ | $x/p^2 + (1-x)/(1-p)^2$ | $1/p + 1/(1-p) = 1/(p(1-p))$ |
| $\mathcal{N}(\mu, \sigma^2)$, $\sigma$ known | $-(x - \mu)^2/(2\sigma^2) + \text{const}$ | $1/\sigma^2$ | $1/\sigma^2$ |
| Poisson$(\lambda)$ | $x\log\lambda - \lambda - \log x!$ | $x/\lambda^2$ | $1/\lambda$ |
| Exponential$(\lambda)$ | $\log\lambda - \lambda x$ | $1/\lambda^2$ | $1/\lambda^2$ |

For $\mathcal{N}(\mu, v)$ with both unknown, take the Hessian from §3.5(d) with $n = 1$: $-\partial^2_\mu = 1/v$; $-\partial_\mu\partial_v = (x - \mu)/v^2$, mean $0$; $-\partial^2_v = -\frac{1}{2v^2} + \frac{(x-\mu)^2}{v^3}$, mean $-\frac{1}{2v^2} + \frac{1}{v^2} = \frac{1}{2v^2}$. So
$$I_1(\mu, \sigma^2) = \begin{pmatrix}1/\sigma^2 & 0\\ 0 & 1/(2\sigma^4)\end{pmatrix}, \qquad I_1^{-1} = \begin{pmatrix}\sigma^2 & 0\\ 0 & 2\sigma^4\end{pmatrix}. \tag{3.12}$$

**Theorem 3.13 (Information adds).** For $n$ independent observations from the model, the information is $nI_1(\theta)$.

*Proof.* $\ell_n = \sum_i\log p(X_i;\theta)$, so $\nabla^2\ell_n = \sum_i\nabla^2\log p(X_i;\theta)$; take $-\mathbb{E}$ of each of the $n$ identically distributed terms. ∎ (Equivalently: $s_n$ is a sum of $n$ independent mean-zero scores, and variances of independent sums add, F5 (3.14).)

**Theorem 3.14 (Reparametrization).** If $\theta = h(\boldsymbol\phi)$ with Jacobian $J = \partial\theta/\partial\boldsymbol\phi$ (LA-14 Theorem 3.2), the information about $\boldsymbol\phi$ is
$$I^{(\phi)}(\boldsymbol\phi) = J^\top I^{(\theta)}\big(h(\boldsymbol\phi)\big)\,J; \qquad \text{one parameter: } I_\phi = h'(\phi)^2 I_\theta. \tag{3.13}$$

*Proof.* By the chain rule (LA-14 Theorem 3.6), $\nabla_{\boldsymbol\phi}\log p(x;h(\boldsymbol\phi)) = J^\top\nabla_\theta\log p$, i.e., $s_\phi = J^\top s_\theta$. Then $\mathbb{E}[s_\phi s_\phi^\top] = J^\top\mathbb{E}[s_\theta s_\theta^\top]J$. ∎

Information transforms like a metric (a covariant 2-tensor, LA-12), not like a function value. Example: for $\mathcal{N}(0, \sigma^2)$ parametrized by $\sigma$, $v = \sigma^2$ and $dv/d\sigma = 2\sigma$, so $I_\sigma = 4\sigma^2\cdot\frac{1}{2\sigma^4} = 2/\sigma^2$. Exercise B4 does the Bernoulli log-odds.

### 3.7 The Cramér–Rao bound and efficiency

**Theorem 3.15 (Cramér–Rao).** Assume (R), and let $\hat\theta$ be an unbiased estimator of a scalar $\theta$ with finite variance, such that $\mathbb{E}_\theta[\hat\theta] = \int\hat\theta(\mathbf{x})\,p_n(\mathbf{x};\theta)\,d\mathbf{x}$ may be differentiated under the integral sign. Then
$$\mathrm{Var}_\theta(\hat\theta) \ge \frac{1}{nI_1(\theta)}. \tag{3.14}$$

▶ **Try it first.** Write unbiasedness as an integral and differentiate in $\theta$. You should find $\mathbb{E}[\hat\theta\,s_n] = 1$. Which inequality relates one covariance to two variances?

**Derivation.**
1. Unbiasedness: $\int\hat\theta(\mathbf{x})\,p_n(\mathbf{x};\theta)\,d\mathbf{x} = \theta$ for all $\theta$. The function $\hat\theta(\mathbf{x})$ does not involve $\theta$.
2. Differentiate: $\int\hat\theta(\mathbf{x})\,\partial_\theta p_n\,d\mathbf{x} = 1$. With $\partial_\theta p_n = p_n s_n$, this says $\mathbb{E}_\theta[\hat\theta\,s_n] = 1$.
3. $\mathbb{E}[s_n] = 0$ (Theorem 3.12, summed), so $\mathrm{Cov}(\hat\theta, s_n) = \mathbb{E}[\hat\theta s_n] - \mathbb{E}[\hat\theta]\mathbb{E}[s_n] = 1$. **This is the covariance that equals 1: estimator with total score.**
4. Cauchy–Schwarz for covariances (F5 Theorem 3.17; it is LA-8 Theorem 3.4 in the inner product $\mathbb{E}[UV]$ on centered variables): $\mathrm{Cov}(\hat\theta, s_n)^2 \le \mathrm{Var}(\hat\theta)\,\mathrm{Var}(s_n)$.
5. $\mathrm{Var}(s_n) = nI_1$ (Theorems 3.12 and 3.13). So $1 \le \mathrm{Var}(\hat\theta)\cdot nI_1$. ∎

■ **What just happened.** In the space of random variables with inner product $\mathbb{E}[UV]$, the score is the direction along which the data move when $\theta$ moves, and it has length $\sqrt{nI_1}$. Unbiasedness forces $\hat\theta - \theta$ to have inner product exactly 1 with that direction. A vector with a fixed inner product with a given direction is shortest when it points along that direction, and then its length is $1/\sqrt{nI_1}$. That is the whole proof.

**Variants** (same proof, different step 3):
- **Biased estimators.** If $\mathbb{E}\hat\theta = \theta + b(\theta)$, step 2 gives $\mathrm{Cov}(\hat\theta, s_n) = 1 + b'(\theta)$, so
$$\mathrm{Var}(\hat\theta) \ge \frac{(1 + b'(\theta))^2}{nI_1(\theta)}. \tag{3.15}$$
For $\hat\sigma^2$, $b = -\sigma^2/n$ and the bound is $(1 - 1/n)^2\cdot 2\sigma^4/n$; the actual variance $2(n-1)\sigma^4/n^2$ is above it by the factor $n/(n-1)$. So $\hat\sigma^2$ legitimately has variance *below* the unbiased bound $2\sigma^4/n$, by the factor $(n-1)/n$ (project Part 5).
- **Unbiased estimators of $g(\theta)$:** $\mathrm{Var} \ge g'(\theta)^2/(nI_1)$, since now $\mathrm{Cov}(\hat g, s_n) = g'(\theta)$.
- **Multiparameter.** If $\hat\theta \in \mathbb{R}^k$ is unbiased, $\mathrm{Cov}(\hat\theta) \succeq (nI_1)^{-1}$: $\mathbf{a}^\top\mathrm{Cov}(\hat\theta)\mathbf{a} \ge \mathbf{a}^\top(nI_1)^{-1}\mathbf{a}$ for all $\mathbf{a}$ (Exercise C3 proves it). In particular $\mathrm{Var}(\hat\theta_j) \ge [(nI_1)^{-1}]_{jj}$, which is at least $1/[nI_1]_{jj}$: having to estimate the other parameters can only cost you.

**Definition 3.16.** An unbiased estimator is **efficient** if it attains (3.14); its **efficiency** is $\frac{1/(nI_1)}{\mathrm{Var}(\hat\theta)} \le 1$.

- $\bar X$ for a normal mean: $\mathrm{Var} = \sigma^2/n = 1/(n\cdot\sigma^{-2})$. Efficient.
- $\hat p$ (Bernoulli), $\hat\lambda = \bar X$ (Poisson): variances $p(1-p)/n$ and $\lambda/n$ equal $1/(nI_1)$. Efficient, at every $n$.
- $S^2$ for $\sigma^2$ with $\mu$ unknown: the bound from (3.12) is $[I_1^{-1}]_{22}/n = 2\sigma^4/n$, while $\mathrm{Var}(S^2) = 2\sigma^4/(n-1)$. Efficiency $(n-1)/n$: not efficient at any finite $n$. (No unbiased estimator does better; $S^2$ is the best unbiased one, by the Lehmann–Scheffé theorem, Casella & Berger §7.3.)

**Theorem 3.17 (When the bound is attained).** Under (R), an unbiased $\hat\theta$ attains (3.14) at every $\theta$ if and only if
$$s_n(\theta;\mathbf{x}) = nI_1(\theta)\,\big(\hat\theta(\mathbf{x}) - \theta\big), \tag{3.16}$$
and then $p_n(\mathbf{x};\theta) = e^{C(\mathbf{x})}\exp\big(\eta(\theta)\hat\theta(\mathbf{x}) - B(\theta)\big)$ for some functions $C, \eta, B$: an exponential family (§3.10) with $\hat\theta$ as its statistic.

▶ **Try it first.** When is Cauchy–Schwarz an equality?

**Derivation.**
1. Equality in step 4 of Theorem 3.15 holds iff the two centered variables are proportional (F5 Theorem 3.17): $s_n = a(\theta)(\hat\theta - \theta)$ with probability 1.
2. Find $a$. $\mathrm{Cov}(s_n, \hat\theta) = a\,\mathrm{Var}\,\hat\theta = 1$ and $\mathrm{Var}(s_n) = a^2\mathrm{Var}\,\hat\theta = nI_1$. Dividing, $a = nI_1$: that is (3.16).
3. Integrate (3.16) in $\theta$: $\partial_\theta\log p_n = a(\theta)\hat\theta(\mathbf{x}) - a(\theta)\theta$. With $\eta(\theta) = \int a\,d\theta$ and $B(\theta) = \int a(\theta)\theta\,d\theta$, $\log p_n = \eta(\theta)\hat\theta(\mathbf{x}) - B(\theta) + C(\mathbf{x})$. ∎

■ **What just happened.** The bound is attainable only when the score is affine in the estimator. That forces exponential-family form and also decides *which* parameter can be estimated efficiently: the mean of the natural statistic. For Exponential$(\lambda)$ data, $s_n = n/\lambda - \sum X_i = -n(\bar X - 1/\lambda)$ is affine in $\bar X$, so $1/\lambda$ is efficiently estimable (by $\bar X$), and $\lambda$ itself is not, by any unbiased estimator (Exercise B3).

### 3.8 Asymptotics of the MLE

**Theorem 3.18 (KL identification).** Let $M(\theta) = \mathbb{E}_{\theta_0}[\log p(X;\theta)]$. Then
$$M(\theta_0) - M(\theta) = D\big(p_{\theta_0}\,\|\,p_\theta\big) \ge 0, \tag{3.17}$$
with equality iff $p_\theta = p_{\theta_0}$. If the model is **identifiable** ($p_\theta = p_{\theta_0} \Rightarrow \theta = \theta_0$), $\theta_0$ is the unique maximizer of $M$.

▶ **Try it first.** Write $M(\theta_0) - M(\theta)$ as one expectation of a log-ratio, and name it.

**Derivation.**
1. $M(\theta_0) - M(\theta) = \mathbb{E}_{\theta_0}\big[\log p(X;\theta_0) - \log p(X;\theta)\big] = \mathbb{E}_{\theta_0}\Big[\log\frac{p(X;\theta_0)}{p(X;\theta)}\Big]$.
2. That is $D(p_{\theta_0}\|p_\theta)$, F5 (3.29). Gibbs' inequality (F5 Theorem 3.45): $\ge 0$, with equality iff the distributions agree.
3. Identifiability turns "same distribution" into "same parameter". ∎

■ **What just happened.** This is the conceptual core of likelihood. By the weak law, $\frac1n\ell_n(\theta) \xrightarrow{p} M(\theta)$ for each $\theta$: the average log-likelihood is a noisy copy of a function whose peak is exactly the truth. Maximizing likelihood is minimizing an empirical estimate of the KL divergence from the truth (F5 (3.30)).

**Consistency (sketch, with the gap named).** We have $\frac1n\ell_n \to M$ pointwise and $\theta_0 = \arg\max M$. Pointwise convergence of functions does not make their maximizers converge; you need the convergence to be uniform in $\theta$ and the maximum of $M$ to be well separated. Both hold, for example, when $\Theta$ is compact and $\log p(x;\theta)$ is continuous in $\theta$ and dominated by an integrable function (Wald 1949; van der Vaart, *Asymptotic Statistics*, Theorem 5.7). Then $\hat\theta_n \xrightarrow{p} \theta_0$. We assume consistency from here on.

We need one tool from outside F5. **Slutsky's theorem** (cited; van der Vaart Lemma 2.8): if $A_n \xrightarrow{d} A$ and $B_n \xrightarrow{p} c$ for a constant $c$, then $A_n + B_n \xrightarrow{d} A + c$, $B_nA_n \xrightarrow{d} cA$, and $A_n/B_n \xrightarrow{d} A/c$ if $c \ne 0$.

**Theorem 3.19 (Asymptotic normality of the MLE).** Assume (R), a consistent MLE, and a third-derivative bound $|\partial^3_\theta\log p(x;\theta)| \le K(x)$ for $\theta$ near $\theta_0$ with $\mathbb{E}K(X) < \infty$. Then
$$\sqrt n\,(\hat\theta_n - \theta_0) \xrightarrow{d} \mathcal{N}\big(0,\; 1/I_1(\theta_0)\big) \qquad \big(\text{vector case: } \mathcal{N}(\mathbf{0}, I_1(\theta_0)^{-1})\big). \tag{3.18}$$

▶ **Try it first.** Taylor-expand $\ell_n'(\hat\theta)$ around $\theta_0$ to first order and solve for $\hat\theta - \theta_0$. Which limit theorem handles the numerator, and which the denominator?

**Derivation** (scalar). Each bracket names the assumption being used.
1. $\hat\theta$ maximizes $\ell_n$ at an interior point, so $\ell_n'(\hat\theta) = 0$. [(R3), consistency: $\hat\theta$ is near the interior point $\theta_0$ with probability $\to 1$.]
2. Taylor with the Lagrange remainder: $\ell_n'(\hat\theta) = \ell_n'(\theta_0) + \ell_n''(\tilde\theta)(\hat\theta - \theta_0)$ for some $\tilde\theta$ between $\hat\theta$ and $\theta_0$. [(R2).]
3. Set it to zero, rearrange, and multiply by $\sqrt n$:
$$\sqrt n(\hat\theta - \theta_0) = \frac{n^{-1/2}\,\ell_n'(\theta_0)}{-n^{-1}\,\ell_n''(\tilde\theta)}.$$
4. Numerator: $n^{-1/2}\sum_i s(\theta_0;X_i)$, a scaled sum of iid variables with mean 0 and variance $I_1(\theta_0)$ (Theorem 3.12). By the CLT (F5 (3.24)) it converges in distribution to $\mathcal{N}(0, I_1(\theta_0))$. [(R4): finite variance.]
5. Denominator: $-\frac1n\ell_n''(\theta_0) = \frac1n\sum_i\big(-\partial^2_\theta\log p(X_i;\theta_0)\big) \xrightarrow{p} I_1(\theta_0)$ by the weak law and (3.11). Moving from $\theta_0$ to $\tilde\theta$ costs at most $|\tilde\theta - \theta_0|\cdot\frac1n\sum_iK(X_i)$, which is $\xrightarrow{p} 0\cdot\mathbb{E}K = 0$ because $\tilde\theta \xrightarrow{p} \theta_0$. [Third-derivative bound; consistency.]
6. Slutsky: the ratio converges in distribution to $\mathcal{N}(0, I_1)/I_1 = \mathcal{N}(0, 1/I_1)$. [(R4): $I_1 > 0$.] ∎

■ **What just happened.** The MLE is the truth plus (score noise)/(curvature). The score noise has variance $nI_1$, the curvature is $nI_1$, so the error variance is $nI_1/(nI_1)^2 = 1/(nI_1)$: **the MLE attains the Cramér–Rao bound asymptotically**. It is asymptotically efficient, and unlike Theorem 3.15 this needs no unbiasedness. Every bracket above is a place where the conclusion can fail: Uniform$(0, \theta)$ breaks step 1 (no zero derivative at the maximum), a true value on the boundary of $\Theta$ breaks (R3), and a likelihood with several roots (§3.12) means step 1 must be applied to the right root.

**Theorem 3.20 (Delta method).** If $\sqrt n(\hat\theta_n - \theta_0) \xrightarrow{d} \mathcal{N}(0, v)$ and $g$ is continuously differentiable near $\theta_0$, then
$$\sqrt n\big(g(\hat\theta_n) - g(\theta_0)\big) \xrightarrow{d} \mathcal{N}\big(0,\; g'(\theta_0)^2 v\big) \qquad \big(\text{vector: } \mathcal{N}(0, \nabla g^\top V\nabla g)\big). \tag{3.19}$$

*Proof.* 1. Mean value theorem: $g(\hat\theta) = g(\theta_0) + g'(\tilde\theta)(\hat\theta - \theta_0)$ with $\tilde\theta$ between them. 2. $\hat\theta - \theta_0 = n^{-1/2}\cdot\sqrt n(\hat\theta - \theta_0) \xrightarrow{p} 0$ by Slutsky, so $\tilde\theta \xrightarrow{p} \theta_0$ and $g'(\tilde\theta) \xrightarrow{p} g'(\theta_0)$ by continuous mapping. 3. $\sqrt n(g(\hat\theta) - g(\theta_0)) = g'(\tilde\theta)\cdot\sqrt n(\hat\theta - \theta_0) \xrightarrow{d} g'(\theta_0)\,\mathcal{N}(0, v)$ by Slutsky. ∎

For the MLE, $v = 1/I_1$, so $g(\hat\theta)$ (the MLE of $g(\theta)$, by invariance) has asymptotic variance $g'(\theta_0)^2/(nI_1)$, the Cramér–Rao bound for $g$: asymptotically efficient too.

**Standard errors in practice.** $\theta_0$ is unknown, so plug in the estimate. The **Wald standard error** and the **observed-information** standard error are
$$\widehat{\mathrm{se}} = \frac{1}{\sqrt{nI_1(\hat\theta)}}, \qquad \widehat{\mathrm{se}}_{\text{obs}} = \frac{1}{\sqrt{J_n}}, \quad J_n = -\ell_n''(\hat\theta). \tag{3.20}$$
They coincide in many exponential families (Poisson: $J_n = T/\hat\lambda^2 = n/\hat\lambda = nI_1(\hat\lambda)$) and differ in general; Efron and Hinkley (1978) argued for the observed version. The interval $\hat\theta \pm 1.96\,\widehat{\mathrm{se}}$ is the Wald 95% interval of ST-2. Example 4.2 checks the normal approximation against the exact Binomial.

### 3.9 Sufficiency and Rao–Blackwell

**Definition 3.21.** A statistic $T = T(\mathbf{X})$ is **sufficient** for $\theta$ if the conditional distribution of $\mathbf{X}$ given $T$ does not depend on $\theta$.

Once you know $T$, the rest of the data is noise with a known distribution: you could generate a fake dataset with the same distribution as the real one without knowing $\theta$.

**Theorem 3.22 (Factorization, Fisher–Neyman).** $T$ is sufficient iff $p_n(\mathbf{x};\theta) = g(T(\mathbf{x});\theta)\,h(\mathbf{x})$ for some nonnegative functions $g, h$. (Proof for discrete $\mathbf{X}$.)

▶ **Try it first.** For "$\Leftarrow$", compute $P(\mathbf{X} = \mathbf{x} \mid T = t)$ and watch $g$ cancel.

*Proof.* ($\Leftarrow$) 1. Fix $t$ and $\mathbf{x}$ with $T(\mathbf{x}) = t$ (otherwise the conditional probability is 0, free of $\theta$). Since $\{\mathbf{X} = \mathbf{x}\} \subseteq \{T = t\}$, $P(\mathbf{X} = \mathbf{x} \mid T = t) = P(\mathbf{X} = \mathbf{x})/P(T = t)$. 2. $P(T = t) = \sum_{\mathbf{y}:T(\mathbf{y}) = t}g(t;\theta)h(\mathbf{y}) = g(t;\theta)\sum_{T(\mathbf{y}) = t}h(\mathbf{y})$. 3. The ratio is $h(\mathbf{x})/\sum_{T(\mathbf{y}) = t}h(\mathbf{y})$: no $\theta$.
($\Rightarrow$) 4. If $T$ is sufficient, $p_n(\mathbf{x};\theta) = P_\theta(T = T(\mathbf{x}))\cdot P(\mathbf{X} = \mathbf{x} \mid T = T(\mathbf{x}))$. The second factor is free of $\theta$; call it $h(\mathbf{x})$, and call the first $g(T(\mathbf{x});\theta)$. ∎

(The continuous case has the same statement; its proof needs conditional densities, Casella & Berger Theorem 6.2.6.) Factorizing the joint densities of §3.5: Bernoulli $p^{\sum x}(1-p)^{n-\sum x}$, Poisson $\lambda^{\sum x}e^{-n\lambda}\cdot\prod(x_i!)^{-1}$ and Exponential $\lambda^ne^{-\lambda\sum x}$ all have $T = \sum x_i$; the Normal density $(2\pi\sigma^2)^{-n/2}\exp\big(-\frac{\sum x^2 - 2\mu\sum x + n\mu^2}{2\sigma^2}\big)$ has $T = (\sum x_i, \sum x_i^2)$; Uniform$(0,\theta)$, $\theta^{-n}\mathbb{1}\{\max x \le \theta\}\cdot\mathbb{1}\{\min x \ge 0\}$, has $T = \max x_i$; Gamma, $(\beta^\alpha/\Gamma(\alpha))^n(\prod x)^{\alpha-1}e^{-\beta\sum x}$, has $T = (\sum\log x_i, \sum x_i)$. Every MLE in §3.5 is a function of these, as it must be: the likelihood sees the data only through $g(T(\mathbf{x});\theta)$.

**Minimal sufficiency.** The whole dataset is trivially sufficient. A sufficient $T$ is **minimal** if it is a function of every other sufficient statistic: maximal compression with no information about $\theta$ lost. Criterion (Lehmann and Scheffé, 1950): $T$ is minimal if $p_n(\mathbf{x};\theta)/p_n(\mathbf{y};\theta)$ is free of $\theta$ exactly when $T(\mathbf{x}) = T(\mathbf{y})$. For Poisson the ratio is $\lambda^{\sum x_i - \sum y_i}\prod y_i!/\prod x_i!$, so $\sum x_i$ is minimal; the same check works for each statistic above.

**Theorem 3.23 (Rao–Blackwell).** Let $T$ be sufficient and $\hat\theta$ an estimator with finite variance. Then $\tilde\theta = \mathbb{E}[\hat\theta \mid T]$ is an estimator, has the same bias as $\hat\theta$, and $\mathrm{Var}(\tilde\theta) \le \mathrm{Var}(\hat\theta)$; hence $\mathrm{MSE}(\tilde\theta) \le \mathrm{MSE}(\hat\theta)$, with equality iff $\hat\theta$ is (with probability 1) a function of $T$.

▶ **Try it first.** Adam's law for the mean, Eve's law for the variance. Where does sufficiency enter?

**Derivation.**
1. $\mathbb{E}[\hat\theta \mid T = t]$ is computed from the conditional distribution of $\mathbf{X}$ given $T = t$, which by sufficiency does not involve $\theta$. So $\tilde\theta$ is a statistic. (Without sufficiency it would depend on the unknown $\theta$ and be useless.)
2. Adam's law (F5 Theorem 3.26): $\mathbb{E}\tilde\theta = \mathbb{E}[\mathbb{E}[\hat\theta \mid T]] = \mathbb{E}\hat\theta$. Same bias.
3. Eve's law (F5 (3.21)): $\mathrm{Var}\,\hat\theta = \mathbb{E}[\mathrm{Var}(\hat\theta \mid T)] + \mathrm{Var}(\tilde\theta) \ge \mathrm{Var}(\tilde\theta)$.
4. Same bias and no larger variance, so by (3.2) no larger MSE. Equality iff $\mathbb{E}[\mathrm{Var}(\hat\theta \mid T)] = 0$, i.e., $\hat\theta$ is determined by $T$. ∎

■ **What just happened.** Averaging an estimator over the part of the data that $T$ certifies as pure noise can only help. In F5's projection picture (Theorem 3.29), $\tilde\theta$ is the orthogonal projection of $\hat\theta$ onto the functions of $T$. If $T$ is also *complete* (Casella & Berger §6.2), the Rao–Blackwellized unbiased estimator is the unique best unbiased estimator (Lehmann–Scheffé); we cite this without proof. Example 4.7 and Exercise B5 use the theorem.

### 3.10 Exponential families

**Definition 3.24.** A model is an **exponential family** if
$$p(x;\theta) = h(x)\exp\big(\boldsymbol\eta(\theta)^\top\mathbf{T}(x) - A(\boldsymbol\eta)\big), \qquad A(\boldsymbol\eta) = \log\int h(x)\,e^{\boldsymbol\eta^\top\mathbf{T}(x)}\,dx, \tag{3.21}$$
with **natural parameter** $\boldsymbol\eta \in \mathbb{R}^k$, **sufficient statistic** $\mathbf{T}(x) \in \mathbb{R}^k$, **base measure** $h \ge 0$, and **log-partition function** $A$ (a sum for discrete $x$). For $n$ iid observations, $p_n = \prod h(x_i)\cdot\exp\big(\boldsymbol\eta^\top\sum\mathbf{T}(x_i) - nA(\boldsymbol\eta)\big)$, so $\sum_i\mathbf{T}(x_i)$ is sufficient by Theorem 3.22.

| Family | $h(x)$ | $\mathbf{T}(x)$ | $\boldsymbol\eta$ | $A(\boldsymbol\eta)$ |
|---|---|---|---|---|
| Bernoulli$(p)$ | $1$ | $x$ | $\log\frac{p}{1-p}$ | $\log(1 + e^\eta)$ |
| Poisson$(\lambda)$ | $1/x!$ | $x$ | $\log\lambda$ | $e^\eta$ |
| Exponential$(\lambda)$ | $\mathbb{1}\{x > 0\}$ | $x$ | $-\lambda$ | $-\log(-\eta)$ |
| $\mathcal{N}(\mu, \sigma^2)$ | $1/\sqrt{2\pi}$ | $(x, x^2)$ | $(\mu/\sigma^2,\, -1/(2\sigma^2))$ | $-\frac{\eta_1^2}{4\eta_2} - \frac12\log(-2\eta_2)$ |
| Gamma$(\alpha, \beta)$ | $\mathbb{1}\{x > 0\}/x$ | $(\log x, x)$ | $(\alpha, -\beta)$ | $\log\Gamma(\eta_1) - \eta_1\log(-\eta_2)$ |

Check the Normal row: $-\frac{(x-\mu)^2}{2\sigma^2} - \frac12\log(2\pi\sigma^2) = \frac{\mu}{\sigma^2}x - \frac{x^2}{2\sigma^2} - \frac{\mu^2}{2\sigma^2} - \frac12\log\sigma^2 - \frac12\log 2\pi$, and with $\eta_1 = \mu/\sigma^2$, $\eta_2 = -1/(2\sigma^2)$, $\frac{\mu^2}{2\sigma^2} = -\frac{\eta_1^2}{4\eta_2}$ and $\frac12\log\sigma^2 = -\frac12\log(-2\eta_2)$. Uniform$(0,\theta)$ (moving support) and the Cauchy are not exponential families.

**Theorem 3.25 (Derivatives of the log-partition function).** For $\boldsymbol\eta$ in the interior of $\{\boldsymbol\eta : A(\boldsymbol\eta) < \infty\}$,
$$\nabla A(\boldsymbol\eta) = \mathbb{E}_{\boldsymbol\eta}[\mathbf{T}], \qquad \nabla^2A(\boldsymbol\eta) = \mathrm{Cov}_{\boldsymbol\eta}(\mathbf{T}) \succeq 0, \tag{3.22}$$
so $A$ is convex.

▶ **Try it first.** Exponentiate: $e^{A(\eta)} = \int h\,e^{\eta T}dx$. Differentiate both sides in $\eta$.

**Derivation** (one parameter, then the vector version).
1. $e^{A(\eta)} = \int h(x)e^{\eta T(x)}dx$. Differentiate under the integral (legitimate in the interior; Casella & Berger Theorem 2.4.3): $A'(\eta)e^{A(\eta)} = \int T\,h\,e^{\eta T}dx$.
2. Divide by $e^{A}$: $A'(\eta) = \int T\,h\,e^{\eta T - A}dx = \mathbb{E}_\eta[T]$.
3. Differentiate step 1 again: $\big(A'' + (A')^2\big)e^A = \int T^2h\,e^{\eta T}dx$, so $A'' + (A')^2 = \mathbb{E}[T^2]$ and $A'' = \mathrm{Var}(T)$.
4. With $\partial/\partial\eta_j$ and $\partial^2/\partial\eta_j\partial\eta_k$ the same steps give $\nabla^2A = \mathbb{E}[\mathbf{T}\mathbf{T}^\top] - \mathbb{E}\mathbf{T}\,\mathbb{E}\mathbf{T}^\top = \mathrm{Cov}(\mathbf{T})$, which is PSD (F5 Theorem 3.19). A function with PSD Hessian is convex (LA-14 Theorem 3.13). ∎

■ **What just happened.** $A$ is a cumulant generating function in disguise: its derivatives are the cumulants of $\mathbf{T}$. Check: Poisson $A = e^\eta$ gives $A' = A'' = e^\eta = \lambda$ (mean = variance); Bernoulli $A' = e^\eta/(1 + e^\eta) = p$, $A'' = p(1-p)$. **SM-4's canonical ensemble is this object.** $p_i = e^{-\beta E_i}/Z$ is an exponential family with $T = E$ (energy), $\eta = -\beta$, $A = \ln Z$. Then $\langle E\rangle = -\partial\ln Z/\partial\beta$ and $\mathrm{Var}(E) = \partial^2\ln Z/\partial\beta^2 = k_BT^2C_V$ are (3.22): same function, same derivatives, different names.

**Theorem 3.26 (MLE in an exponential family).** The log-likelihood $\ell(\boldsymbol\eta) = \sum\log h(x_i) + \boldsymbol\eta^\top\sum\mathbf{T}(x_i) - nA(\boldsymbol\eta)$ is concave in $\boldsymbol\eta$, and the MLE (when it exists) solves the moment-matching equation
$$\mathbb{E}_{\hat{\boldsymbol\eta}}[\mathbf{T}] = \bar{\mathbf{T}} := \frac1n\sum_i\mathbf{T}(x_i). \tag{3.23}$$
The Fisher information in $\boldsymbol\eta$ is $I_1(\boldsymbol\eta) = \nabla^2A(\boldsymbol\eta)$, and the observed information equals the expected information.

▶ **Try it first.** Compute $\nabla\ell$ and $\nabla^2\ell$. Notice that the Hessian contains no data.

**Derivation.**
1. $\nabla\ell = \sum\mathbf{T}(x_i) - n\nabla A(\boldsymbol\eta) = n\big(\bar{\mathbf{T}} - \mathbb{E}_{\boldsymbol\eta}[\mathbf{T}]\big)$ by (3.22). It vanishes iff (3.23) holds.
2. $\nabla^2\ell = -n\nabla^2A(\boldsymbol\eta) = -n\,\mathrm{Cov}_{\boldsymbol\eta}(\mathbf{T}) \preceq 0$: $\ell$ is concave. If $\mathrm{Cov}(\mathbf{T})$ is positive definite (no nontrivial linear combination of $\mathbf{T}$ is constant, a **minimal** family), $\ell$ is strictly concave and a critical point is the unique global maximum (LA-8 Theorem 3.19 plus concavity).
3. $\nabla^2\ell$ contains no data, so $-\nabla^2\ell$ equals its own expectation: observed information $= n\nabla^2A(\boldsymbol\eta) =$ expected information, and by (3.11) $I_1(\boldsymbol\eta) = \nabla^2A(\boldsymbol\eta)$. ∎

■ **What just happened.** In an exponential family, maximum likelihood *is* the method of moments, applied to $\mathbf{T}$. For the Gamma, $\mathbb{E}[\log X] = \partial A/\partial\eta_1 = \psi(\alpha) - \log\beta$, and matching it to $\frac1n\sum\log x_i$ (with $\mathbb{E}X = \bar x$) is exactly (3.10). Concavity is why GLM fitting (ST-4) has a unique answer and Newton's method on it (IRLS) is safe. It holds in $\boldsymbol\eta$; other parametrizations need not be concave, though by invariance they have the same maximizer.

**Theorem 3.27 (Maximum entropy).** Fix $\mathbf{T}$ and take $h = 1$ (counting or Lebesgue measure). Among all distributions $q$ with $\mathbb{E}_q[\mathbf{T}] = \boldsymbol\tau$, the entropy $H(q) = -\mathbb{E}_q[\log q]$ is maximized by $p_{\boldsymbol\eta}(x) = \exp(\boldsymbol\eta^\top\mathbf{T}(x) - A(\boldsymbol\eta))$ with $\boldsymbol\eta$ chosen so that $\mathbb{E}_{\boldsymbol\eta}\mathbf{T} = \boldsymbol\tau$ (when such $\boldsymbol\eta$ exists).

*Where the form comes from (outline).* Maximize $-\sum_xq(x)\log q(x)$ subject to $\sum q = 1$ and $\sum q\mathbf{T} = \boldsymbol\tau$ with multipliers $\lambda_0, \boldsymbol\eta$. Stationarity in each $q(x)$: $-\log q(x) - 1 + \lambda_0 + \boldsymbol\eta^\top\mathbf{T}(x) = 0$, so $q(x) \propto e^{\boldsymbol\eta^\top\mathbf{T}(x)}$. This is F5 Theorem 3.46's calculation with general $\mathbf{T}$.

*Proof that it is the maximum.* For any $q$ with $\mathbb{E}_q\mathbf{T} = \boldsymbol\tau$: $-\mathbb{E}_q[\log p_{\boldsymbol\eta}] = -\boldsymbol\eta^\top\boldsymbol\tau + A(\boldsymbol\eta)$, which is the same number as $-\mathbb{E}_{p_{\boldsymbol\eta}}[\log p_{\boldsymbol\eta}] = H(p_{\boldsymbol\eta})$, because both expectations of $\mathbf{T}$ equal $\boldsymbol\tau$. So $H(p_{\boldsymbol\eta}) - H(q) = -\mathbb{E}_q[\log p_{\boldsymbol\eta}] + \mathbb{E}_q[\log q] = D(q\,\|\,p_{\boldsymbol\eta}) \ge 0$ by Gibbs. ∎

Fix the mean on $(0, \infty)$ and you get the Exponential; mean and variance on $\mathbb{R}$, the Normal; mean energy, Boltzmann (SM-4). Setting $\boldsymbol\tau = \bar{\mathbf{T}}$ reproduces (3.23): maximum entropy and maximum likelihood pick the same member of the family.

### 3.11 Bayesian estimation, previewed

Treat $\theta$ as random with a **prior** $\pi(\theta)$. Bayes' rule gives the **posterior**
$$\pi(\theta \mid \mathbf{x}) = \frac{p_n(\mathbf{x};\theta)\,\pi(\theta)}{\int p_n(\mathbf{x};\theta')\,\pi(\theta')\,d\theta'} \propto L(\theta)\,\pi(\theta). \tag{3.24}$$
The **posterior mean** $\hat\theta_B = \mathbb{E}[\theta \mid \mathbf{x}]$ is a function of the data, so it is an estimator with a sampling distribution, a bias and a variance like any other. For Bernoulli data with a $\mathrm{Beta}(a, b)$ prior, F5 Theorem 3.12 gives $p \mid \mathbf{x} \sim \mathrm{Beta}(a + T, b + n - T)$ and
$$\hat p_B = \frac{a + T}{a + b + n} = w\,\hat p + (1 - w)\frac{a}{a + b}, \qquad w = \frac{n}{a + b + n}. \tag{3.25}$$
A weighted average of the MLE and the prior mean: **shrinkage**. At a fixed true $p$, its bias is $(1 - w)\big(\frac{a}{a+b} - p\big)$ and its variance is $w^2p(1-p)/n$, smaller than the MLE's by $w^2$. Laplace's rule of succession is $a = b = 1$, $\hat p_B = (T + 1)/(n + 2)$. With $n = 50$ and $T = 15$ it gives $0.3077$ instead of $0.3$, and if the truth is $p = 0.3$ its MSE is $0.003942$ against $0.0042$ for $\hat p$. It wins near $1/2$ and loses near 0 and 1 (Exercise B6 finds the crossovers $0.1447$ and $0.8553$). The posterior mean minimizes the prior-averaged MSE (it is a conditional expectation, F5 Theorem 3.29). ST-5 is the full subject.

### 3.12 What can go wrong

- **Non-identifiability.** If $p_\theta = p_{\theta'}$ for $\theta \ne \theta'$, (3.17) only locates the set $\{\theta : p_\theta = p_{\theta_0}\}$. In $X \sim \mathcal{N}(\theta_1 + \theta_2, 1)$ only the sum is estimable; mixture labels can be swapped (Exercise B8).
- **No MLE, or one on the boundary.** Bernoulli data with $T = 0$ give $\hat p = 0$ with Wald standard error 0: nonsense. Logistic regression on separable data (ST-4) has likelihood increasing forever as coefficients grow. A Gaussian mixture with free variances has unbounded likelihood: center a component on a data point and shrink its variance to 0.
- **Multiple maxima.** The Cauchy location likelihood can have several local maxima (Exercise C2); mixtures always do. Theorem 3.19 concerns the root near $\theta_0$: start from a consistent estimator (the median, for Cauchy) and take Newton steps.
- **Neyman–Scott (1948).** Pairs $X_{i1}, X_{i2} \sim \mathcal{N}(\mu_i, \sigma^2)$, $i = 1, \dots, n$, each with its own mean. The MLE of $\sigma^2$ converges to $\sigma^2/2$ (Exercise C1). The number of parameters grows with $n$, each $\mu_i$ is estimated from two points forever, and the argument behind Theorem 3.18 never applies to them. The fix is to eliminate the nuisance means, e.g. via $X_{i1} - X_{i2} \sim \mathcal{N}(0, 2\sigma^2)$ (REML in ST-3 generalizes this).
- **Misspecification.** The data come from some $q$ that is not in the model.

**Theorem 3.28 (The MLE targets the KL projection).** If $X_i$ are iid from $q$, then $\frac1n\ell_n(\theta) \xrightarrow{p} -H(q) - D(q\,\|\,p_\theta)$, so the MLE targets
$$\theta^* = \arg\min_\theta D\big(q\,\|\,p_\theta\big). \tag{3.26}$$

*Derivation.* 1. Weak law under $q$: $\frac1n\sum\log p(X_i;\theta) \xrightarrow{p} \mathbb{E}_q[\log p_\theta]$. 2. $\mathbb{E}_q[\log p_\theta] = \mathbb{E}_q[\log q] - \mathbb{E}_q\big[\log\frac{q}{p_\theta}\big] = -H(q) - D(q\|p_\theta)$. 3. $H(q)$ does not depend on $\theta$, so maximizing the limit means minimizing $D(q\|p_\theta)$; with the same uniform-convergence caveat as before, $\hat\theta \xrightarrow{p} \theta^*$. ∎ Theorem 3.18 is the special case $q = p_{\theta_0}$, where the minimum is 0.

*Example.* Fit $\mathcal{N}(\mu, \sigma^2)$ to Exponential$(\lambda)$ data. $-\mathbb{E}_q[\log p_{\mu,\sigma^2}] = \frac12\log(2\pi\sigma^2) + \mathbb{E}_q[(X - \mu)^2]/(2\sigma^2)$ is minimized at $\mu^* = \mathbb{E}_qX = 1/\lambda$ and $\sigma^{*2} = \mathrm{Var}_qX = 1/\lambda^2$ (moment matching again: (3.23) with the expectation under $q$). The gap is $D = \frac12\log(2\pi e\sigma^{*2}) - h(q) = \frac12\log(2\pi e) - 1 \approx 0.4189$ nats for every $\lambda$, using $h(\mathrm{Exp}(\lambda)) = 1 - \log\lambda$. Project Part 10 confirms both. Under misspecification (3.11) fails, since it needed the data to come from the model, and the asymptotic variance becomes the **sandwich** $J^{-1}KJ^{-1}$ with $J = -\mathbb{E}_q[\nabla^2\log p_{\theta^*}]$, $K = \mathrm{Cov}_q(\nabla\log p_{\theta^*})$ (Huber 1967; White 1982), used in ST-2 and ST-3.

## 4. Worked examples

### Example 4.1: Two variance estimators, $n = 10$

Normal data, $n = 10$, in units of $\sigma^4$ (set $\sigma^2 = 1$). From (3.4)–(3.5):
- $S^2$: bias 0, variance $2/9$, so $\mathrm{MSE} = 0.2222$.
- $\hat\sigma^2$: bias $-1/10$, variance $2\cdot 9/100 = 0.18$, so $\mathrm{MSE} = 0.18 + 0.01 = 0.19$.

The biased estimator wins. Exercise B1's divisor $n + 1$ does better still, $\mathrm{MSE} = 2/11 = 0.1818$. Project Part 2 simulates $10^4$ datasets and confirms each figure to within Monte Carlo error; Part 1 does the same at $n = 5, 20, 100$ with $\sigma^2 = 4$, where for example at $n = 5$ the theory gives $\mathrm{MSE}(S^2) = 8.0$ and $\mathrm{MSE}(\hat\sigma^2) = 5.76$. At $n = 100$ the three MSEs differ in the third digit: the trade-off matters only for small samples.

### Example 4.2: Bernoulli: MLE, information, Cramér–Rao attained, exact vs asymptotic

$n = 50$ trials, $T = 15$ successes. $\hat p = 0.3$. The curvature at the peak is $\ell''(\hat p) = -15/0.09 - 35/0.49 = -166.67 - 71.43 = -238.1$, which equals $-n/(\hat p(1 - \hat p)) = -50/0.21$. Per observation, $I_1(\hat p) = 1/0.21 = 4.7619$. The Wald standard error is $1/\sqrt{238.1} = 0.0648$ and the 95% Wald interval is $(0.173, 0.427)$.

Cramér–Rao is attained *exactly*, at every $n$: $\mathrm{Var}(\hat p) = p(1-p)/n = 1/(nI_1(p))$. Theorem 3.17 predicted this: $s_n = \frac{T}{p} - \frac{n - T}{1 - p} = \frac{T - np}{p(1-p)} = nI_1(p)(\hat p - p)$, affine in $\hat p$.

**Exact versus asymptotic.** Suppose the truth is $p = 0.3$. How likely is $\hat p \le 0.2$, i.e., $T \le 10$? Exactly, $P(T \le 10) = 0.0789$ (Binomial CDF). The asymptotic $\hat p \approx \mathcal{N}(0.3, 0.0042)$ gives $z = (0.2 - 0.3)/0.0648 = -1.543$ and probability $0.0614$, noticeably low. With the continuity correction (F5 §3.11), $\Phi\big((10.5 - 15)/\sqrt{10.5}\big) = 0.0825$, much closer. At $n = 50$ Theorem 3.19 is a good description of the center of the sampling distribution and a rough one of the tails.

### Example 4.3: Normal $(\mu, \sigma^2)$ jointly: MLE, Fisher matrix, inverse

Eight measurements of a nominal 10 Ω resistor: $9.8, 10.4, 10.1, 9.6, 10.3, 10.0, 9.7, 10.1$. The sum is $80.0$, so $\hat\mu = 10.0$. Deviations: $-0.2, 0.4, 0.1, -0.4, 0.3, 0.0, -0.3, 0.1$; their squares sum to $0.56$. So $\hat\sigma^2 = 0.56/8 = 0.07$ (and $S^2 = 0.08$).

By §3.5(d) the Hessian at the MLE is $\mathrm{diag}(-n/\hat\sigma^2, -n/(2\hat\sigma^4))$, and by (3.12) the Fisher matrix for the sample is
$$nI_1(\hat\theta) = \begin{pmatrix}8/0.07 & 0\\ 0 & 8/(2\cdot 0.0049)\end{pmatrix} = \begin{pmatrix}114.29 & 0\\ 0 & 816.33\end{pmatrix}, \qquad (nI_1)^{-1} = \begin{pmatrix}0.00875 & 0\\ 0 & 0.001225\end{pmatrix}.$$
Observed and expected information agree here. The standard errors are $\mathrm{se}(\hat\mu) = \sqrt{0.00875} = 0.0935$ Ω and $\mathrm{se}(\hat\sigma^2) = \sqrt{0.001225} = 0.0350$ Ω². The zero off-diagonal says the two estimates are asymptotically uncorrelated; Theorem 3.7 says they are exactly independent.

Project Part 3 draws the log-likelihood surface in $(\mu, \sigma)$ with the MLE and the ellipse $\frac12(\theta - \hat\theta)^\top J(\theta - \hat\theta) = c$ from the observed information $J = \mathrm{diag}(n/\hat\sigma^2, 2n/\hat\sigma^2)$ (Theorem 3.14 with $v = \sigma^2$). The exact contours are lopsided in $\sigma$, steep toward $\sigma \to 0$ and gentle toward large $\sigma$; the ellipse is symmetric. With $n = 8$ the quadratic approximation is only local. Here $\hat\sigma^2 \pm 1.96\,\mathrm{se}$ reaches down almost to zero, a sign to use the exact $\chi^2$ distribution of Theorem 3.7 instead (ST-2).

### Example 4.4: Uniform$(0, \theta)$: the maximum, its bias, the fix, and why Cramér–Rao does not apply

By F5 Theorem 3.9 (with $k = n$), $X_{(n)}$ has density $nx^{n-1}/\theta^n$ on $[0, \theta]$. So
$$\mathbb{E}X_{(n)} = \int_0^\theta x\,\frac{nx^{n-1}}{\theta^n}dx = \frac{n\theta}{n+1}, \qquad \mathbb{E}X_{(n)}^2 = \frac{n\theta^2}{n+2}, \qquad \mathrm{Var}X_{(n)} = \frac{n\theta^2}{(n+1)^2(n+2)}.$$
The MLE is biased low by $\theta/(n+1)$: the sample maximum can never exceed $\theta$. The unbiased correction $\tilde\theta = \frac{n+1}{n}X_{(n)}$ has variance $\frac{(n+1)^2}{n^2}\cdot\frac{n\theta^2}{(n+1)^2(n+2)} = \frac{\theta^2}{n(n+2)}$. At $n = 10$ that is $0.00833\,\theta^2$, against $0.0333\,\theta^2$ for the MoM estimator $2\bar X$: four times better, and the ratio $(n+2)/3$ grows with $n$.

Now compute the "Cramér–Rao bound" naively. $\log p = -\log\theta$ on $[0, \theta]$, so $s = -1/\theta$ and $\mathbb{E}s^2 = 1/\theta^2$, giving $\theta^2/n = 0.1\,\theta^2$ at $n = 10$. The unbiased $\tilde\theta$ beats this by the factor $n + 2 = 12$, and its variance decays like $1/n^2$, faster than any Cramér–Rao rate. There is no contradiction: (R1) fails. Step 1 of Theorem 3.12 breaks: $\mathbb{E}_\theta[s] = -1/\theta \ne 0$, because $\frac{d}{d\theta}\int_0^\theta\frac1\theta dx = 0$ while $\int_0^\theta\partial_\theta\frac1\theta\,dx = -\frac1\theta$; the difference is Leibniz's boundary term $p(\theta;\theta) = 1/\theta$. A sharp edge in the density pins the parameter down much harder than smooth curvature can. Project Part 1 checks all three estimators at $n = 5, 20, 100$; Exercise B2 finds the MSE-optimal multiple of $X_{(n)}$.

### Example 4.5: Gamma by moments and by likelihood

Ten waiting times (minutes): $0.8, 1.9, 2.4, 0.5, 3.1, 1.2, 2.0, 1.6, 0.9, 2.6$. Then $\bar x = 1.7$ and $\hat\sigma^2 = 0.654$ (divisor $n$; squared deviations sum to $6.54$).

*MoM*, (3.7): $\hat\alpha = 1.7^2/0.654 = 4.4190$ and $\hat\beta = 1.7/0.654 = 2.5994$.

*MLE*, (3.10): $\frac1n\sum\log x_i = 0.3928$, so $s = \log 1.7 - 0.3928 = 0.1378$. Solving $\log\alpha - \psi(\alpha) = 0.1378$ numerically gives $\hat\alpha = 3.7865$, and $\hat\beta = \hat\alpha/\bar x = 2.2273$. Project Part 4 gets the same root three ways: a bracketing root-finder on (3.10), direct numerical maximization of $\ell$, and `scipy.stats.gamma.fit`.

The two methods disagree by 15% in $\hat\alpha$. Which to trust? The Fisher matrix in $(\alpha, \beta)$, from $-\partial^2\log p$, is
$$I_1 = \begin{pmatrix}\psi'(\alpha) & -1/\beta\\ -1/\beta & \alpha/\beta^2\end{pmatrix}, \qquad [I_1^{-1}]_{\alpha\alpha} = \frac{\alpha/\beta^2}{(\alpha\psi'(\alpha) - 1)/\beta^2} = \frac{\alpha}{\alpha\psi'(\alpha) - 1}.$$
At $\alpha = 2$ this is $6.90$, against $2\alpha(\alpha + 1) = 12$ for MoM: the MLE needs about 58% as many observations for the same precision. In the project's simulation at $\alpha = 2$, $n = 50$, $n\mathrm{Var}(\hat\alpha_{\text{MoM}}) = 12.53$ and $n\mathrm{Var}(\hat\alpha_{\text{MLE}}) = 8.26$. The MLE is clearly better, though at $n = 50$ still about 20% above its asymptotic value: (3.18) is a large-sample statement.

### Example 4.6: Counting decays: asymptotic normality, a Wald interval, and the delta method

A Geiger counter records counts in twenty one-second windows: $3, 5, 2, 4, 3, 1, 4, 3, 6, 2, 3, 4, 2, 5, 3, 2, 4, 3, 2, 3$ (total 64). Decays in disjoint windows are independent Poisson counts (F4), so the model is Poisson$(\lambda)$ with $\lambda$ the mean rate per second.

$\hat\lambda = 64/20 = 3.2$. By (3.18) with $I_1 = 1/\lambda$, $\sqrt n(\hat\lambda - \lambda) \xrightarrow{d} \mathcal{N}(0, \lambda)$. The Wald standard error is $\sqrt{\hat\lambda/n} = \sqrt{0.16} = 0.4$; the observed information $J_n = T/\hat\lambda^2 = 64/10.24 = 6.25 = n/\hat\lambda$ gives the same. The 95% Wald interval is $3.2 \pm 1.96(0.4) = (2.416, 3.984)$ decays per second. In the project's simulation at $\lambda = 3.2$, $n = 20$, this recipe covers the truth $94.5\%$ of the time.

Now estimate the probability that a one-second window is empty, $g(\lambda) = e^{-\lambda}$. By invariance the MLE is $e^{-3.2} = 0.04076$. By the delta method (3.19), $\mathrm{se} \approx |g'(\hat\lambda)|\cdot 0.4 = e^{-3.2}\cdot 0.4 = 0.01630$. The simulated standard deviation of $e^{-\bar X}$ at $n = 20$ is $0.01779$ (ratio 1.091); at $n = 200$ the ratio is 1.014. The delta method linearizes $e^{-\lambda}$, which curves noticeably over $\pm 0.4$.

### Example 4.7: Rao–Blackwellizing $X_1$ for the Poisson mean

Start from a deliberately bad unbiased estimator of $\lambda$: the first observation, $\hat\theta = X_1$, with variance $\lambda$. Condition on the sufficient statistic $T = \sum X_i$. Using that $T \sim \mathrm{Poisson}(n\lambda)$ and $\sum_{i \ge 2}X_i \sim \mathrm{Poisson}((n-1)\lambda)$ independent of $X_1$ (sums of independent Poissons, F5 §3.4):
$$P(X_1 = k \mid T = t) = \frac{\frac{e^{-\lambda}\lambda^k}{k!}\cdot\frac{e^{-(n-1)\lambda}((n-1)\lambda)^{t-k}}{(t-k)!}}{\frac{e^{-n\lambda}(n\lambda)^t}{t!}} = \binom tk\Big(\frac1n\Big)^k\Big(1 - \frac1n\Big)^{t-k}.$$
$\lambda$ cancels, which is sufficiency at work: $X_1 \mid T = t \sim \mathrm{Binomial}(t, 1/n)$. So $\tilde\theta = \mathbb{E}[X_1 \mid T] = T/n = \bar X$. Rao–Blackwell turned the worst reasonable estimator into the efficient one, and the variance dropped from $\lambda$ to $\lambda/n$. With $\lambda = 3.2$ and $n = 20$, from $3.2$ to $0.16$; the project's $10^4$ simulated datasets give $3.181$ and $0.1593$, and it checks $\mathbb{E}[X_1 \mid T = t] = t/n$ for several $t$. Exercise B5 does the same for $e^{-\lambda}$, where the answer is not something you would guess.

## 5. Why it matters and where it connects

- **ST-2** turns sampling distributions into confidence intervals (pivots such as Theorem 3.7's $\chi^2$, Wald intervals from (3.20)) and tests; Wilks' theorem for likelihood-ratio tests is the Taylor expansion of Theorem 3.19 taken one order further; the bootstrap simulates the sampling distribution from the data.
- **ST-3**: least squares is the Normal MLE; Gauss–Markov is an efficiency theorem; ridge regression is Example 4.1's trade-off in $p$ dimensions; the residual variance with divisor $n - p$ is Theorem 3.7 with $p$ rotated-out directions.
- **ST-4**: a GLM is an exponential family with $\eta = \mathbf{x}^\top\boldsymbol\beta$; Theorem 3.26 gives concavity and $I = \nabla^2A$, so Fisher scoring = Newton = IRLS.
- **ST-5**: conjugate priors exist for every exponential family, generalizing (3.25); the posterior concentrates at the rate $1/\sqrt{nI_1}$ (Bernstein–von Mises).
- **ST-6, ML-10**: variational inference minimizes $D(q\|p)$, the KL divergence in the opposite direction from Theorem 3.28's projection.
- **ML-1**: the bias–variance decomposition of prediction error is (3.2) applied to a fitted function; minimizing cross-entropy is maximizing likelihood (F5 (3.30)).
- **ML-2**: the natural gradient is $I_1^{-1}\nabla\ell$, LA-14 §3.1's Riesz representative with the Fisher metric of Theorem 3.14.
- **ML-6**: EM maximizes mixture likelihoods, the non-concave case of §3.12 and Exercise B8.
- **SM-4**: $\ln Z$ is $A(\eta)$, and the heat capacity is the Fisher information about $\beta$ (Theorem 3.25).

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 drill bias/MSE, an MLE with invariance, method of moments, both forms of Fisher information, an exponential-family rewrite, and factorization. B1–B8 are the module: B1 (best divisor for $\sigma^2$), B3 (the Exponential rate cannot be estimated efficiently), B5 (Rao–Blackwell for $e^{-\lambda}$, which misses Cramér–Rao), and B8 (code: mixture likelihoods) are essential. C1 proves Neyman–Scott inconsistency, C2 (code) is the Cauchy location model, C3 proves the multiparameter Cramér–Rao bound.

## 7. Coding project

`project.py` makes "an estimator is a random variable" literal by drawing $10^4$ datasets for every setting:

1. **Sampling-distribution simulator.** Normal, Exponential, Uniform$(0, \theta)$ at $n = 5, 20, 100$; MoM, MLE and unbiased corrections; histograms of all sampling distributions; a table of bias, variance and MSE against the formulas of §3.3, Exercise B3 and Example 4.4, asserted within 4 Monte Carlo standard errors.
2. **Bias–variance plot** for the three variance estimators against $n$ (Example 4.1, Exercise B1).
3. **Likelihood surface** of Example 4.3 with the MLE and the Fisher ellipse; finite-difference Hessian checked against (3.12).
4. **Gamma MLE** from (3.10) three ways, plus a simulation comparing MoM and MLE with their asymptotic variances (Example 4.5).
5. **Cramér–Rao ratios** $\mathrm{Var}(\text{MLE})\cdot nI_1$ against $n$ for Bernoulli, Poisson (exactly 1), Exponential rate ($n^3/((n-1)^2(n-2))$), and Normal variance ($(n-1)/n$).
6. **QQ plots** of $\sqrt n(\hat\lambda - \lambda)/\lambda$ for the Exponential MLE at $n = 5, 20, 100, 500$, and Example 4.2's exact-vs-normal check.
7. **Wald interval and delta method** (Example 4.6), with simulated coverage.
8. **Rao–Blackwell** (Example 4.7, Exercise B5).
9. **Exponential-family identities**: $A' = \mathbb{E}T$, $A'' = \mathrm{Var}\,T$ by finite differences for Poisson and Bernoulli, and symbolically with SymPy.
10. **Misspecification**: a Normal fitted to Exponential(1) data converges to the KL projection $\mathcal{N}(1, 1)$, found independently by numerically minimizing the KL divergence.

**Build it yourself first:** write `mc_check` (Monte Carlo bias/variance/MSE with standard errors), the vectorized Newton solver `gamma_shape_mle` for (3.10), and the Part 8 Rao–Blackwell comparison before reading the file's versions. "Your turn" then asks for the mixture-weight likelihood (B8), the Cauchy location model (C2), and a bootstrap standard error for the Gamma MLE (a preview of ST-2); reference solutions for B8 and C2 sit at the bottom of the file.

## 8. Where you'll struggle

- **What is random.** The single hardest step in this module is that an estimator is a random variable, and "the standard error" is the standard deviation of its sampling distribution over imagined repetitions, not a property of your one dataset. People compute $\hat p \pm 1.96\sqrt{\hat p(1-\hat p)/n}$ without being able to say what varies. The fix: for every standard error you compute this month, write one sentence saying what is being repeated, then simulate it as in project Part 1.
- **MLE is not a guarantee.** It is not always unbiased ($\hat\sigma^2$, $1/\bar X$), not always consistent (Neyman–Scott), and not always computable (Gamma needs a root-finder; mixtures have many roots). The Uniform example exists to show that the regularity conditions are not decoration: drop (R1) and both Cramér–Rao and the $1/\sqrt n$ rate fail.
- **Which covariance equals 1.** The Cramér–Rao proof is three lines of Cauchy–Schwarz and everyone forgets the key fact. It is $\mathrm{Cov}(\hat\theta, s_n) = 1$, the estimator against the *total* score, from differentiating $\mathbb{E}_\theta\hat\theta = \theta$. If you can derive that line, you can rebuild the whole proof, including the biased and $g(\theta)$ variants.
- **The asymptotic-normality derivation hides an assumption at every step.** Read it with the brackets: interior maximum, Taylor expansion, CLT for the score, LLN for the Hessian, consistency to move from $\theta_0$ to $\tilde\theta$, Slutsky. When an exam or a paper says "under regularity conditions", this list is what they mean.
- **Sufficiency feels like a formality** until Rao–Blackwell (§3.9) and exponential families (§3.10) use it. Do Example 4.7 by hand; watching $\lambda$ cancel is what makes the definition stick.
- **The exponential-family section is dense** and pays off three modules later, in ST-4 and ST-5. Read §3.10 twice, and fill in the table yourself for one family before checking it.
- **Parametrization conventions.** Rate vs scale for Exponential and Gamma, variance vs standard deviation for the Normal: the Fisher information changes with the parametrization (3.13), and scipy uses scale. Write the density down before computing anything.

## 9. Check yourself

1. What is the sampling distribution of an estimator, and what is its standard deviation called?
2. Derive $\mathrm{MSE} = b^2 + \mathrm{Var}$ in two lines.
3. Why does $S^2$ divide by $n - 1$? Which estimator of $\sigma^2$ has the smaller MSE for normal data, $S^2$ or $\hat\sigma^2$?
4. Write the MoM and MLE estimators of $\theta$ for Uniform$(0, \theta)$. Which is better, and why does Cramér–Rao not settle it?
5. State the invariance property of the MLE. Is unbiasedness invariant?
6. State the two expressions for Fisher information and the condition that makes them equal.
7. In the Cramér–Rao proof, which covariance equals 1, and why?
8. When can an unbiased estimator attain the Cramér–Rao bound?
9. Why does the average log-likelihood peak at $\theta_0$? Name the inequality.
10. Write the asymptotic distribution of the MLE and the delta-method variance of $g(\hat\theta)$.
11. State the factorization theorem and use it to find a sufficient statistic for Poisson data.
12. In an exponential family, what are $A'(\eta)$ and $A''(\eta)$, and what equation does the MLE solve?

## 10. Pair with

- **Wasserman, *All of Statistics*, ch. 6–7 and 9.** The same material, terser and with more examples; ch. 9 on parametric inference is the closest match to §3.5–§3.8. (Ch. 8, the bootstrap, belongs to ST-2.)
- **MIT 18.650 (Rigollet), lecture slides and notes ch. 1–3** on OCW. A modern statistician's framing of models, identifiability, MLE and Fisher information, with the KL view of §3.8 front and center.
- **Casella & Berger, *Statistical Inference*, ch. 6–7.** The reference: full proofs of factorization in the continuous case, completeness, Lehmann–Scheffé, and many more exercises. Use it when these notes cite it.
- **Rice, *Mathematical Statistics and Data Analysis*, ch. 8.** Gentle and concrete, with real datasets for MoM and MLE and a clear treatment of the Gamma example.
- **Efron & Hastie, *Computer Age Statistical Inference*, ch. 1–4** (free PDF). The big picture: how frequentist, Bayesian and Fisherian inference relate, and why Fisher's information still organizes the field.
