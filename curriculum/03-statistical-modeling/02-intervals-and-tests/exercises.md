# ST-2 Exercises

Conventions as in the notes: quantiles are upper-tail ($P(Z>z_a) = a$, so $z_{0.025} = 1.960$; $\chi^2_{k,0.975}$ is the small one); $S^2$ has divisor $n-1$ and $\hat\sigma^2$ divisor $n$; β(θ) is the power function and $\beta_{\rm II}$ a type II error rate. Use SciPy or tables for quantiles.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** A sample of $n = 16$ Normal observations has $\bar x = 52.3$ and $s = 4.8$.
(a) Give the 95% $t$ interval for μ.
(b) If σ were known to be 5, give the 95% $z$ interval. Why is it narrower than (a) even though $5>4.8$?
(c) Give the 95% interval for σ from (3.4).

**A2.** In $n = 10$ trials there are 0 successes. Compute the 95% Wald interval (3.5) and the 95% Wilson interval (3.6). Compare Wilson's upper end with the "rule of three" bound $3/n$. What does the Wald interval claim, and why is the claim absurd?

**A3.** You toss a coin 10 times to test $H_0:p = 1/2$ against $H_1:p>1/2$, rejecting when the number of heads $X\ge8$.
(a) Find the size of the test.
(b) Find the power at $p = 0.7$ and the type II error rate there.
(c) Repeat (a) and (b) for the rule $X\ge9$. Describe the trade-off in one sentence.

**A4.** Eight p-values: 0.001, 0.012, 0.020, 0.024, 0.040, 0.046, 0.30, 0.70. With α = $q$ = 0.05, which hypotheses do Bonferroni, Holm and Benjamini–Hochberg reject? One p-value is rejected by BH even though it exceeds its own threshold $jq/m$. Which, and why?

**A5.** A paper reports: "mean difference 1.75, 95% CI [0.4, 3.1], $p = 0.02$" for a two-sided test of no difference. True or false, with one sentence of reason each:
(a) There is a 95% probability that the true difference lies in [0.4, 3.1].
(b) If the study were repeated many times, about 95% of the intervals computed this way would contain the true difference.
(c) The probability that there is no true difference is 0.02.
(d) A level-0.05 test of $H_0$: difference = 0 rejects.
(e) A level-0.05 test of $H_0$: difference = 3 rejects.
(f) If there were no true difference, data at least this extreme would occur with probability 0.02.

**A6.** How many observations does a one-sample two-sided $z$ test at α = 0.05 need for power 0.90 against a shift of $0.3\sigma$? How many per group for a two-sample comparison of the same shift?

## B. Core (do all)

**B1.** (Projection for $n = 2$.) Let $X_1,X_2$ be iid $\mathcal N(\mu,\sigma^2)$.
(a) Show that $S^2 = (X_1-X_2)^2/2$.
(b) Let $U = (X_1+X_2)/\sqrt2$ and $V = (X_1-X_2)/\sqrt2$. Find the joint distribution of $(U,V)$ and conclude that $\bar X$ and $S^2$ are independent and that $S^2/\sigma^2\sim\chi^2_1$.
(c) Which orthonormal vectors $\mathbf q_1,\mathbf q_2$ of Theorem 3.1 did you use?

**B2.** (An exact pivot.) Let $X_1,\dots,X_n$ be iid Exponential with rate λ (density $\lambda e^{-\lambda x}$, mean $1/\lambda$).
(a) Show that $2\lambda\sum X_i\sim\chi^2_{2n}$, so it is a pivot.
(b) Eight component lifetimes (hours) are 12.1, 3.4, 27.9, 8.6, 15.2, 5.3, 41.0, 9.5. Give the MLE $\hat\lambda$, the exact 95% interval for λ, and the corresponding interval for the mean lifetime $1/\lambda$.
(c) Why is the interval for λ not centred at $\hat\lambda$?

**B3.** (Neyman–Pearson and UMP.) Same model as B2.
(a) For $H_0:\lambda = \lambda_0$ versus $H_1:\lambda = \lambda_1$ with $\lambda_1<\lambda_0$, use Theorem 3.7 to find the most powerful level-α test. Show that it rejects for large $\sum X_i$ and give the threshold using B2(a).
(b) Explain why the same test is uniformly most powerful for $H_1:\lambda<\lambda_0$.
(c) A supplier claims a mean lifetime of 10 hours ($\lambda_0 = 0.1$); you suspect it is longer. With α = 0.05 and the data of B2, find the threshold, the decision, the p-value, and the power at $\lambda_1 = 0.05$.

**B4.** (The trinity, and who is invariant.) Same data, $H_0:\lambda = 0.1$ against a two-sided alternative.
(a) Compute the LRT, Wald and score statistics (3.16) in the rate parametrization.
(b) Reparametrize by the mean $\theta = 1/\lambda$ ($H_0:\theta = 10$). Recompute all three. Which change?
(c) Give the $\chi^2_1$ p-values and compare them with the exact two-sided p-value from B2's pivot, $2\min\{P(\chi^2_{16}\ge v),P(\chi^2_{16}\le v)\}$ with $v = 2\lambda_0\sum x_i$. Which asymptotic statistic would you report, and why?

**B5.** (The bootstrap of the mean, exactly.) Fix data $x_1,\dots,x_n$ with mean $\bar x$ and $\hat\sigma^2 = \frac1n\sum(x_i-\bar x)^2$.
(a) Prove that under resampling, $\mathbb E^*[\bar X^*] = \bar x$ and $\operatorname{Var}^*(\bar X^*) = \hat\sigma^2/n$.
(b) For $\mathbf x = (1,2,6)$, list or count the 27 equally likely ordered resamples and verify (a) numerically.
(c) Compare $\widehat{\rm se}_{\rm boot}$ with $s/\sqrt n$. What factor relates them, and why does it not matter for large $n$?

**B6.** (Holm.) Prove that Holm's step-down procedure controls FWER at α under arbitrary dependence. Hint: let $i^*$ be the true null with the smallest p-value, and ask how many hypotheses can be rejected before Holm reaches it. Then prove that Holm rejects every hypothesis Bonferroni rejects.

**B7.** (What a p-value looks like when $H_0$ is false.) In a one-sided $z$ test, $Z = \sqrt n(\bar X-\mu_0)/\sigma\sim\mathcal N(\Delta,1)$ with $\Delta = \sqrt n\delta/\sigma$, and $p = 1-\Phi(Z)$.
(a) Show that $P(p\le u) = \Phi(\Delta-z_u)$. What is this at $u = \alpha$?
(b) Show that the density of $p$ is $g(u) = \exp(\Delta z_u-\Delta^2/2)$ on $(0,1)$. What is $g$ when $\Delta = 0$?
(c) With $\Delta = 2.5$, compute the power at α = 0.05, $g(0.04)$ and $g(0.001)$. A study with 80% power reports $p = 0.04$. By what factor is that p-value more likely under $H_1$ than under $H_0$? Is that strong evidence?

**B8.** (Code: bootstrap a regression slope, an ST-3 preview.) In `project.py`, write `slope(a)` returning the least-squares slope $\sum(x_i-\bar x)(y_i-\bar y)/\sum(x_i-\bar x)^2$ along axis 1 of an array of $(x,y)$ pairs, so that it plugs into `bootstrap_ci`. Generate 1,000 datasets with $n = 25$, $x_i\sim\mathrm{Unif}(0,10)$, $y_i = 1+0.5x_i+\varepsilon_i$, $\varepsilon_i\sim\mathcal N(0,2^2)$. For each, compute the percentile and basic bootstrap intervals ($B = 1000$) and the classical interval $\hat b\pm t_{23,0.025}\,{\rm se}(\hat b)$ from `scipy.stats.linregress`. Report the three coverages. Why do you resample $(x,y)$ pairs rather than $x$ and $y$ separately?

## C. Challenge (do at least one)

**C1.** (The $t$ density.) Let $Z\sim\mathcal N(0,1)$ and $V\sim\chi^2_k$ be independent, $T = Z/\sqrt{V/k}$. Use the change of variables $(Z,V)\mapsto(T,V)$ (F5 Theorem 3.6) and integrate out $V$ to show
$$f_T(t) = \frac{\Gamma\big(\frac{k+1}{2}\big)}{\sqrt{k\pi}\,\Gamma\big(\frac k2\big)}\Big(1+\frac{t^2}{k}\Big)^{-(k+1)/2}.$$
Check that $k = 1$ gives the Cauchy density and that $f_T(t)\to\varphi(t)$ as $k\to\infty$.

**C2.** (Pearson's $X^2$ is the score statistic.) For the multinomial with $k$ cells, use $\boldsymbol\theta = (p_1,\dots,p_{k-1})$ as free parameters, $p_k = 1-\sum_{j<k}p_j$. Compute the score vector $\mathbf u$ at $\mathbf p^0$, the Fisher information $I(\mathbf p^0) = n\big(\operatorname{diag}(1/p_j^0)+\frac{1}{p_k^0}\mathbf 1\mathbf 1^\top\big)$ (derive it), its inverse $\frac1n(\operatorname{diag}(p_j^0)-\mathbf p^0\mathbf p^{0\top})$ (verify it), and show that the score statistic $\mathbf u^\top I^{-1}\mathbf u$ equals $\sum_{j=1}^k(N_j-np_j^0)^2/(np_j^0)$ exactly.

**C3.** (Code: optional stopping.) Under $H_0:\mu = 0$ with $X_i\sim\mathcal N(0,1)$ and σ known, an analyst computes the two-sided $z$ test after every 10 observations, up to 100, and stops the first time $p<0.05$. Simulate $2\times10^4$ such analysts and estimate the type I error. Repeat with the per-look threshold $0.05/10$. Explain both numbers, and explain why "I'll stop when it's significant" is p-hacking even when every individual test is correct.

**C4.** (Jeffreys–Lindley in numbers.) For $\bar X\sim\mathcal N(\mu,1/n)$, $H_0:\mu = 0$, $H_1:\mu\sim\mathcal N(0,\tau^2)$, derive the Bayes factor formula of §3.11. With $\tau = 1$ and the observed $z = \sqrt n\,\bar x = 2.5$ held fixed, compute $B_{01}$ for $n = 10, 10^2, 10^3, 10^4, 10^6$, and find the $n$ at which $B_{01} = 1$. What is the p-value in every case? Explain in two sentences why the two answers diverge.
