# ST-1 Exercises

Conventions as in the notes: Exponential and Gamma use the **rate**; $\mathcal{N}(\mu, \sigma^2)$ is parametrized by the variance; $S^2$ divides by $n - 1$, $\hat\sigma^2$ by $n$; $I_1$ is the Fisher information per observation.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Let $X_1, X_2$ be iid with mean $\mu$ and variance $\sigma^2$. Three estimators of $\mu$: $\hat\theta_a = (X_1 + X_2)/2$, $\hat\theta_b = (X_1 + 2X_2)/3$, $\hat\theta_c = (X_1 + X_2)/3$. Compute the bias, variance and MSE of each. Which of the two unbiased ones is better? For which values of $\mu$ does $\hat\theta_c$ beat $\hat\theta_a$? What does the answer say about choosing estimators by MSE?

**A2.** The Geometric distribution counts trials up to and including the first success: $P(X = k) = (1 - p)^{k-1}p$, $k = 1, 2, \dots$. Data: $3, 1, 4, 2, 5$.
(a) Find the MLE of $p$.
(b) Find the MLEs of the mean $1/p$ and of $P(X = 1)$ without new calculus.
(c) Compute $I_1(p)$ and the Wald standard error of $\hat p$.

**A3.** Let $X_i \sim \mathrm{Uniform}(-\theta, \theta)$.
(a) Why is the first-moment equation useless here?
(b) Derive the method-of-moments estimator from the second moment.
(c) For the data $-1.2, 0.4, 2.1, -0.7, 1.5$, compute the MoM estimate and the MLE (sketch $L(\theta)$ to find the MLE).

**A4.** The Rayleigh density is $p(x;\theta) = \frac{x}{\theta}e^{-x^2/(2\theta)}$ for $x > 0$, $\theta > 0$. You may use that $X^2/(2\theta) \sim \mathrm{Exponential}(1)$.
(a) Compute $I_1(\theta)$ both ways, as $\mathbb{E}[s^2]$ and as $-\mathbb{E}[\partial_\theta^2\log p]$.
(b) Find the MLE, show it is unbiased, and show it attains the Cramér–Rao bound.

**A5.** Write the Geometric distribution of A2 in exponential-family form (3.21): identify $h$, $T$, $\eta$, and $A(\eta)$. Use Theorem 3.25 to find its mean and variance, and check them against $1/p$ and $(1 - p)/p^2$.

**A6.** Use the factorization theorem to find a sufficient statistic for an iid sample from:
(a) Geometric$(p)$;
(b) $p(x;\theta) = \theta x^{\theta - 1}$ on $(0, 1)$, $\theta > 0$;
(c) Uniform$(\theta, \theta + 1)$.
In (c), compare the dimension of the sufficient statistic with the dimension of the parameter.

## B. Core (do all; these are the module)

**B1.** (The best divisor.) For normal data with $\mu$ unknown, consider $\hat\sigma^2_c = c\sum_i(X_i - \bar X)^2$.
(a) Using Theorem 3.7, find the bias, variance and MSE of $\hat\sigma^2_c$ as functions of $c$.
(b) Show that the MSE is minimized at $c = 1/(n + 1)$, and find the minimum.
(c) Evaluate $\mathrm{MSE}/\sigma^4$ at $n = 10$ for $c = 1/(n-1)$, $1/n$, $1/(n+1)$.
(d) Is $c = 1/(n+1)$ best for every value of $\sigma^2$? Give two reasons nobody uses it.

**B2.** (Uniform estimators.) For Uniform$(0, \theta)$, consider estimators $cX_{(n)}$.
(a) Using Example 4.4's moments, find $\mathrm{MSE}(c)$ and the minimizing $c^*$.
(b) At $n = 10$, compute $\mathrm{MSE}/\theta^2$ for $2\bar X$, $X_{(n)}$, $\frac{n+1}{n}X_{(n)}$, and $c^*X_{(n)}$.
(c) Explain in words why the order-statistic estimators have MSE of order $1/n^2$ while $2\bar X$ has order $1/n$.

**B3.** (A rate cannot be estimated efficiently.) Let $X_i \sim \mathrm{Exponential}(\lambda)$ and $S = \sum X_i \sim \mathrm{Gamma}(n, \lambda)$ (F5 Theorem 3.8).
(a) Show $\mathbb{E}[1/S] = \lambda/(n - 1)$ for $n \ge 2$ and $\mathbb{E}[1/S^2] = \lambda^2/((n-1)(n-2))$ for $n \ge 3$.
(b) Find the bias and variance of the MLE $\hat\lambda = n/S$ and of $\tilde\lambda = (n-1)/S$.
(c) Compare $\mathrm{Var}(\tilde\lambda)$ with the Cramér–Rao bound and compute its efficiency.
(d) Use Theorem 3.17 to explain why no unbiased estimator of $\lambda$ attains the bound, while $\bar X$ attains it for $1/\lambda$.
(e) Give the numbers at $n = 10$ in units of $\lambda$ or $\lambda^2$.

**B4.** (Reparametrization and the delta method.) Bernoulli$(p)$ with log-odds $\eta = \log(p/(1-p))$.
(a) Compute $I_1(\eta)$ by the Jacobian rule (3.13).
(b) Compute it again as $A''(\eta)$ (§3.10).
(c) For $n = 50$, $T = 15$, find $\hat\eta$ and its Wald standard error.
(d) Get the same standard error from $\mathrm{se}(\hat p)$ by the delta method.
(e) Why can it be better to build a Wald interval for $\eta$ and transform its endpoints back to $p$?

**B5.** (Rao–Blackwell for $e^{-\lambda}$.) Poisson$(\lambda)$ data, $T = \sum X_i$, target $g(\lambda) = e^{-\lambda} = P(X = 0)$.
(a) Show $W = \mathbb{1}\{X_1 = 0\}$ is unbiased for $g$, and find its variance.
(b) Using Example 4.7's conditional distribution, show $\mathbb{E}[W \mid T] = (1 - 1/n)^T$.
(c) Verify directly, with the Poisson probability generating function $\mathbb{E}[z^T] = e^{n\lambda(z-1)}$, that $(1 - 1/n)^T$ is unbiased, and show its variance is $e^{-2\lambda}(e^{\lambda/n} - 1)$.
(d) Compute the Cramér–Rao bound for unbiased estimators of $g$. Show the Rao–Blackwellized estimator never attains it, and explain why using Theorem 3.17.
(e) Evaluate everything at $\lambda = 3.2$, $n = 20$, and compare with the simulated MSE of the MLE $e^{-\bar X}$ in project Part 8.

**B6.** (Shrinkage.) Bernoulli data, Laplace's estimator $\tilde p = (T + 1)/(n + 2)$ (the posterior mean under a uniform prior).
(a) Find its bias, variance and MSE as functions of $p$.
(b) For $n = 50$, find the values of $p$ at which $\mathrm{MSE}(\tilde p) = \mathrm{MSE}(\hat p)$.
(c) Which estimator would you report, and what does your answer depend on?

**B7.** (Exponential families attain the bound.) Let $p(x;\eta) = h(x)\exp(\eta T(x) - A(\eta))$ be a one-parameter exponential family and $\mu = A'(\eta)$ its mean parameter.
(a) Show $\bar T = \frac1n\sum T(X_i)$ is unbiased for $\mu$ with variance $A''(\eta)/n$.
(b) Use (3.13) to show the Fisher information about $\mu$ is $1/A''(\eta)$.
(c) Conclude that $\bar T$ attains the Cramér–Rao bound for $\mu$.
(d) Check (a)–(c) for Poisson and for Exponential, and reconcile the Exponential case with B3.

**B8.** (Code: mixture likelihoods.) Consider the model $w\,\mathcal{N}(m_1, 1) + (1 - w)\,\mathcal{N}(m_2, 1)$.
(a) Prove that for fixed $(m_1, m_2)$ the log-likelihood is concave in $w$.
(b) Show that $(w, m_1, m_2)$ and $(1 - w, m_2, m_1)$ give the same likelihood. Use this to show that if $\ell$ were concave in all three parameters, the midpoint of the two maximizers would also be a maximizer, and explain what that midpoint model is.
(c) In `project.py`, simulate $n = 400$ points from $0.3\,\mathcal{N}(-1, 1) + 0.7\,\mathcal{N}(2, 1)$, maximize the likelihood from the starts $(0.5, -2, 3)$ and $(0.5, 3, -2)$, and check (b) numerically, including the likelihood at the midpoint. Check (a) on a grid of $w$.

## C. Challenge (do at least one)

**C1.** (Neyman–Scott.) Let $X_{i1}, X_{i2} \overset{\text{iid}}{\sim} \mathcal{N}(\mu_i, \sigma^2)$ for $i = 1, \dots, n$, all independent, with $\mu_1, \dots, \mu_n, \sigma^2$ unknown.
(a) Find the joint MLE. Show $\hat\sigma^2 = \frac{1}{4n}\sum_i(X_{i1} - X_{i2})^2$.
(b) Show $\hat\sigma^2 \xrightarrow{p} \sigma^2/2$.
(c) Which assumption behind Theorem 3.18 and the consistency sketch fails?
(d) Give a consistent estimator of $\sigma^2$, and check (b) with a simulation.

**C2.** (Code: the Cauchy location model.) $p(x;\theta) = \frac{1}{\pi(1 + (x - \theta)^2)}$.
(a) Show $I_1(\theta) = 1/2$, although $\mathrm{Var}(X)$ does not exist. (Hint: $\int_{-\infty}^\infty \frac{u^2}{(1+u^2)^3}du$ is a Beta integral.)
(b) Explain why $\bar X$ is useless here (F5 §3.11, "What the CLT does not say").
(c) In `project.py`, for $n = 100$ and 2000 replications, compute the MLE (start at the median, then maximize), the median, and the mean. Compare $n\mathrm{Var}(\text{MLE})$ with $1/I_1 = 2$ and $n\mathrm{Var}(\text{median})$ with its asymptotic value $\pi^2/4$, and measure the spread of $\bar X$ by its interquartile range.

**C3.** (Multiparameter Cramér–Rao.) Let $\hat\theta \in \mathbb{R}^k$ be unbiased for $\theta$ under (R).
(a) Show $\mathbb{E}[\hat\theta\,s_n^\top] = I_k$ (the identity matrix) by differentiating unbiasedness.
(b) For fixed $\mathbf{a}, \mathbf{b} \in \mathbb{R}^k$, apply Cauchy–Schwarz to $\mathbf{a}^\top\hat\theta$ and $\mathbf{b}^\top s_n$, then choose $\mathbf{b}$ to prove $\mathbf{a}^\top\mathrm{Cov}(\hat\theta)\mathbf{a} \ge \mathbf{a}^\top(nI_1)^{-1}\mathbf{a}$.
(c) Prove that for a positive definite $M$, $[M^{-1}]_{jj} \ge 1/M_{jj}$. Interpret for the Gamma model of Example 4.5: compare the asymptotic variance of $\hat\alpha$ when $\beta$ is unknown with the value it would have if $\beta$ were known, at $\alpha = 2$.
