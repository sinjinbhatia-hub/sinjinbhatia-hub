# ST-5: Bayesian Inference I: Priors, Posteriors, Conjugacy, Credible Intervals, Predictive Distributions, Bayesian Regression, Model Comparison

**Time:** ~12 hours · **Prerequisites:** ST-4 (and through it ST-1 §3.10, exponential families; F5, Beta–Binomial and the Gaussian conditional) · **Unlocks:** ST-6, ML-4, ML-10

> **Conventions in this module.** The prior is $p(\theta)$ (ST-1 wrote $\pi(\theta)$), the likelihood $p(y\mid\theta)$, the posterior $p(\theta\mid y)\propto p(y\mid\theta)p(\theta)$; "$\propto$" drops factors that do not involve the quantity being updated. The **evidence** (marginal likelihood) is $p(y)=\int p(y\mid\theta)p(\theta)\,d\theta$ and the **posterior predictive** is $p(\tilde y\mid y)=\int p(\tilde y\mid\theta)p(\theta\mid y)\,d\theta$, where $\tilde y$ is a future observation. Hyperparameters are Greek with subscript 0 for the prior ($\alpha_0,\beta_0,\mu_0,\kappa_0$) and subscript $n$ for the posterior after $n$ observations ($\alpha_n,\mu_n,\dots$). $\mathcal N(\mu,\sigma^2)$ takes the **variance**; **precision** is its reciprocal. $\mathrm{Gamma}(\alpha,\beta)$ has shape $\alpha$ and **rate** $\beta$ (as in ST-1); $\mathrm{IG}(a,b)$ has density $b^a x^{-a-1}e^{-b/x}/\Gamma(a)$ and mean $b/(a-1)$; $t_\nu(m,s^2)$ is the Student $t$ with location $m$ and scale $s$, density $\propto[1+(x-m)^2/(\nu s^2)]^{-(\nu+1)/2}$. The **Bayes factor** is $\mathrm{BF}_{10}=p(y\mid M_1)/p(y\mid M_0)$; ST-2 (3.23) wrote $B_{01}$, the reciprocal, so watch the subscripts. A **credible** interval is a posterior probability statement; a **confidence** interval (ST-2) is a procedure property. $\hat\theta$ is the MLE, $I_1(\theta)$ the per-observation Fisher information, $\sigma(t)=1/(1+e^{-t})$. Logs are natural.

## 0. Where this sits

You need: from ST-1, the exponential family (3.21) with $\mathbb E\,T=\nabla A$ (3.22), Fisher information and its reparametrization rule (3.13), and asymptotic normality of the MLE (Theorem 3.19); from ST-2, the $t$ interval (Example 4.2) and the one-page Bayes factor (3.23); from ST-3, ridge (3.38)–(3.39) and the lasso (3.43); from ST-4, IRLS (3.13) and the Gamma–Poisson mixture (3.31); from F5, Beta–Binomial conjugacy (Theorem 3.12), the Gaussian conditional (Theorem 3.24), and Adam's and Eve's laws (Theorems 3.26, 3.28); from LA-14, the Hessian and second-order Taylor expansion (Theorem 3.13).

This module replaces "the parameter is a fixed unknown number" by "the parameter has a distribution that data reshapes". Closed-form updating works for a small list of models (§§3.2–3.6) and Bayesian linear regression (§3.11). Everything else needs computation, which is ST-6. ML-4's Gaussian processes are §3.11 with infinitely many features; ML-10's ELBO is the evidence (3.36) bounded from below.

## 1. The problem that forced this idea

Thomas Bayes (died 1761) asked: if an event has happened $y$ times in $n$ trials, what can you say about its unknown chance? Richard Price communicated Bayes's essay to the Royal Society in 1763. Pierre-Simon Laplace reached the same inversion independently in 1774 and stated it generally: the probability of a cause given an effect is proportional to the probability the cause gives the effect (times the cause's prior). He used it for the proportion of male births, asking how sure one could be that a boy is more probable than a girl, and for astronomical quantities such as the mass of Saturn. For about a century "inverse probability", with uniform priors as the default, was ordinary practice.

Then it fell out of favor. A uniform prior is not uniform after reparametrization (§3.7 derives this), and a prior looked like an opinion smuggled into science. Fisher (maximum likelihood, 1920s) and Neyman (confidence intervals, 1930s) built methods that needed no prior. Jeffreys (*Theory of Probability*, 1939) kept the Bayesian program alive with priors built from invariance; de Finetti, Savage, Lindley and Jaynes defended it from coherence, decision theory and logic. What held it back in practice was the integral (3.2): closed forms existed only for conjugate models. Metropolis et al. (1953) had the sampling idea, and Gelfand and Smith (1990) showed that Gibbs sampling made posterior integrals routine. Since the 1990s it has been a working method in astrophysics, genomics, forecasting and machine learning.

The forcing problem is a question frequentist inference does not answer. ST-1 and ST-2 tell you how a procedure behaves over hypothetical repetitions. You often need instead: *given this data, what should I believe about $\theta$, and what should I forecast?* You want prior knowledge in (last year's measurements, physical bounds), and you want uncertainty about $\theta$ carried through to predictions rather than dropped by plugging in $\hat\theta$. If you accept the arguments of Cox, de Finetti and Savage that consistent degrees of belief obey the probability calculus, the only coherent object for that job is a probability distribution over $\theta$. The cost is that you must supply a prior. §8 is honest about what that costs.

## 2. Intuition

**Two curves multiply.** Take a coin with unknown chance $p$ of heads. Before data, your belief is a curve over $[0,1]$: say broad, centered at $\tfrac12$ (Beta(2,2)). After 7 heads in 10 flips, the likelihood is another curve, peaked at $0.7$. The posterior is their pointwise product, renormalized:

```
 belief
   ^        prior          likelihood           posterior (product)
   |         .-.               .-.                    .-.
   |       .'   '.            .' '.                  .' '.
   |     .'       '.        .'     '.              .'     '.
   +-----0---.5-----1    0----.7----1          0----.64---1
```

The product peaks between the two peaks, closer to the likelihood because 10 flips outweigh a prior worth about 4 flips, and it is narrower than either: data sharpen belief.

**Conjugacy is counting.** Multiply $p^{\alpha_0-1}(1-p)^{\beta_0-1}$ by $p^y(1-p)^{n-y}$ and the exponents add. Updating a Beta is adding successes to one hyperparameter and failures to the other: the prior acts as $\alpha_0+\beta_0$ imaginary flips. **Guess the result:** the posterior is a Beta with $\alpha_0+y$ and $\beta_0+n-y$; its mean is a weighted average of the prior mean and the observed fraction, with weights proportional to imaginary and real flips; its spread shrinks like $1/n$. Now guess the general rule. In ST-1's exponential family the data enter only through the sum $\sum T(y_i)$ and the count $n$. A prior that has the same shape as the likelihood in $\eta$ must carry a "prior sum" and a "prior count". The update should be: add the real sum to the prior sum, add $n$ to the prior count. §3.6 proves it, and it contains every closed form in this module.

**The predictive averages over ignorance.** To forecast the next flip, you do not use one value of $p$; you average the forecast over your belief. For one flip the answer is the posterior mean. For ten future flips the answer is *wider* than Binomial(10, any single $p$), because the flips share the unknown $p$: if it is high, many come up heads together.

**Evidence is a bet placed before the data.** A model must spread its probability over every dataset it might produce. A flexible model spreads thin; a rigid model bets big on a narrow set.

```
 p(data | model)
   ^   rigid:  /\
   |          /  \          flexible: _____
   |         /    \                  /       \
   +--------+------+--------------------------> data
```

If the data land where the rigid model bet, it wins by the ratio of heights. If they land elsewhere, the flexible model wins. Complexity is penalized automatically, with no penalty term added by hand.

## 3. Formal development

### 3.1 Bayes' theorem for parameters

F4 Theorem 3.6 says $P(A\mid B)=P(B\mid A)P(A)/P(B)$. With $A$ a statement about $\theta$ and $B$ the data, written with densities:
$$p(\theta\mid y)=\frac{p(y\mid\theta)\,p(\theta)}{p(y)},\qquad p(y)=\int p(y\mid\theta)\,p(\theta)\,d\theta, \tag{3.1, 3.2}$$
where the denominator is the law of total probability (F4 Theorem 3.7). Three distributions (prior, likelihood, posterior) and two integrals: the **evidence** (3.2), which is the density the model assigned to the data before seeing them, and the **posterior predictive**, for $\tilde y$ conditionally independent of $y$ given $\theta$,
$$p(\tilde y\mid y)=\int p(\tilde y\mid\theta)\,p(\theta\mid y)\,d\theta. \tag{3.3}$$
Three consequences, each one line. (i) Only the *shape in $\theta$* of $p(y\mid\theta)p(\theta)$ matters for the posterior; constants drop out. (ii) The evidence is the integral of that unnormalized shape, so if the posterior kernel has a known normalizing constant $Z_n$ and the prior has $Z_0$, then $p(y)=(\text{data factor})\cdot Z_n/Z_0$. We use this three times. (iii) **Sequential equals batch.** If $y_1,y_2$ are conditionally independent given $\theta$, then $p(\theta\mid y_1,y_2)\propto p(y_2\mid\theta)\,[p(y_1\mid\theta)p(\theta)]\propto p(y_2\mid\theta)\,p(\theta\mid y_1)$: yesterday's posterior is today's prior. The project checks this at every checkpoint of Part 1. By ST-1's factorization theorem the posterior depends on the data only through a sufficient statistic.

### 3.2 The Beta–Binomial model

**Theorem 3.1.** Let $p\sim\mathrm{Beta}(\alpha_0,\beta_0)$ and $y\mid p\sim\mathrm{Binomial}(n,p)$. Write $\alpha_n=\alpha_0+y$, $\beta_n=\beta_0+n-y$, $\hat p=y/n$, $m_n=\alpha_n/(\alpha_n+\beta_n)$. Then:

(a) $p\mid y\sim\mathrm{Beta}(\alpha_n,\beta_n)$. $\qquad(3.4)$

(b) $p(y)=\binom ny\,B(\alpha_n,\beta_n)/B(\alpha_0,\beta_0)$. $\qquad(3.5)$

(c) $\mathbb E[p\mid y]=w\,\hat p+(1-w)\dfrac{\alpha_0}{\alpha_0+\beta_0}$ with $w=\dfrac{n}{\alpha_0+\beta_0+n}$. $\qquad(3.6)$

(d) $\mathrm{Var}(p\mid y)=\dfrac{m_n(1-m_n)}{\alpha_0+\beta_0+n+1}$. $\qquad(3.7)$

(e) For $m$ future trials, $P(\tilde y=k\mid y)=\binom mk B(\alpha_n+k,\beta_n+m-k)/B(\alpha_n,\beta_n)$ (the **Beta–Binomial**), with mean $m\,m_n$ and variance $m\,m_n(1-m_n)\,\dfrac{\alpha_n+\beta_n+m}{\alpha_n+\beta_n+1}$. $\qquad(3.8)$

(f) With $\alpha_0=\beta_0=1$ and $m=1$: $P(\tilde y=1\mid y)=(y+1)/(n+2)$, **Laplace's rule of succession**. $\qquad(3.9)$

▶ **Try it first.** Multiply likelihood and prior, keep only powers of $p$, and recognize the kernel. For (c), write $\alpha_n/(\alpha_n+\beta_n)$ as a mixture of $y/n$ and $\alpha_0/(\alpha_0+\beta_0)$. For (e), the predictive integral is the evidence integral with a different prior.

**Derivation.**
1. The likelihood is $\binom ny p^y(1-p)^{n-y}$ and the prior is $p^{\alpha_0-1}(1-p)^{\beta_0-1}/B(\alpha_0,\beta_0)$. The product, as a function of $p$, is $p^{\alpha_0+y-1}(1-p)^{\beta_0+n-y-1}$ times constants: the kernel of $\mathrm{Beta}(\alpha_n,\beta_n)$. A density with this kernel is unique, giving (a) (F5 Theorem 3.12).
2. The evidence is the integral of that product: $\binom ny B(\alpha_n,\beta_n)/B(\alpha_0,\beta_0)$, the ratio of normalizing constants, which is (b). For $\alpha_0=\beta_0=1$ it equals $\binom ny\frac{y!(n-y)!}{(n+1)!}=\frac1{n+1}$ for every $y$: the uniform prior treats all outcome counts as equally likely.
3. $\dfrac{\alpha_0+y}{\alpha_0+\beta_0+n}=\dfrac{n}{\alpha_0+\beta_0+n}\cdot\dfrac yn+\dfrac{\alpha_0+\beta_0}{\alpha_0+\beta_0+n}\cdot\dfrac{\alpha_0}{\alpha_0+\beta_0}$. That is (c).
4. F5's moments give $\mathbb E[p^2]=\frac{\alpha(\alpha+1)}{(\alpha+\beta)(\alpha+\beta+1)}$ for $\mathrm{Beta}(\alpha,\beta)$, so $\mathrm{Var}=\frac{\alpha(\alpha+1)(\alpha+\beta)-\alpha^2(\alpha+\beta+1)}{(\alpha+\beta)^2(\alpha+\beta+1)}=\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}=\frac{m(1-m)}{\alpha+\beta+1}$. With $\alpha+\beta=\alpha_0+\beta_0+n$ this is (d).
5. By (3.3), $P(\tilde y=k\mid y)=\int\binom mk p^k(1-p)^{m-k}\,\mathrm{Beta}(p;\alpha_n,\beta_n)\,dp$, which is step 2 with $(\alpha_0,\beta_0,n,y)$ replaced by $(\alpha_n,\beta_n,m,k)$. For the variance use Eve's law (F5 Theorem 3.28) with $\tilde y\mid p\sim\mathrm{Bin}(m,p)$: $\mathrm{Var}(\tilde y)=m\,\mathbb E[p(1-p)]+m^2\mathrm{Var}(p)$. Since $\mathbb E[p(1-p)]=m_n(1-m_n)-\mathrm{Var}(p)$, this is $m\,m_n(1-m_n)+m(m-1)\mathrm{Var}(p\mid y)$, and substituting (d) gives the factor $1+\frac{m-1}{\alpha_n+\beta_n+1}$.
6. For $m=1$, $\alpha_0=\beta_0=1$: $P(\tilde y=1\mid y)=m_n=(y+1)/(n+2)$. ∎

■ **What just happened.** The posterior mean is a shrinkage estimator: it pulls $\hat p$ toward the prior mean by an amount fixed by the ratio of imaginary to real trials, $\alpha_0+\beta_0$ against $n$. Once $n\gg\alpha_0+\beta_0$ the data dominate ($w\to1$) and the variance (3.7) is $\approx\hat p(1-\hat p)/n=1/(nI_1)$, the Cramér–Rao scale of ST-1. The predictive variance has the extra factor $(\alpha_n+\beta_n+m)/(\alpha_n+\beta_n+1)>1$: future trials share the unknown $p$, so they are positively correlated. Plugging in $\hat p$ or $m_n$ throws that away. Laplace applied (3.9) to the sun rising tomorrow, and noted that other knowledge would make the probability far larger; the rule is only as good as the prior.

### 3.3 The Normal model with known variance

**Theorem 3.2.** Let $y_i\mid\mu\sim\mathcal N(\mu,\sigma^2)$ be independent with $\sigma^2$ known, and $\mu\sim\mathcal N(\mu_0,\tau_0^2)$. Then $\mu\mid y\sim\mathcal N(\mu_n,\tau_n^2)$ with
$$\frac1{\tau_n^2}=\frac1{\tau_0^2}+\frac n{\sigma^2},\qquad \mu_n=\tau_n^2\Big(\frac{\mu_0}{\tau_0^2}+\frac{n\bar y}{\sigma^2}\Big), \tag{3.10}$$
and $\tilde y\mid y\sim\mathcal N(\mu_n,\ \sigma^2+\tau_n^2)$. $\qquad(3.11)$

▶ **Try it first.** Replace the likelihood by $\exp[-n(\bar y-\mu)^2/(2\sigma^2)]$ (why is that legitimate?), multiply by the prior, and complete the square in $\mu$.

**Derivation.**
1. $\sum(y_i-\mu)^2=\sum(y_i-\bar y)^2+n(\bar y-\mu)^2$, and the first term does not involve $\mu$. So the likelihood in $\mu$ is $\propto\exp[-n(\bar y-\mu)^2/(2\sigma^2)]$: $\bar y$ is sufficient.
2. The exponent of likelihood times prior is $-\tfrac12\big[\tfrac n{\sigma^2}(\mu-\bar y)^2+\tfrac1{\tau_0^2}(\mu-\mu_0)^2\big]$. Expand: the coefficient of $\mu^2$ is $-\tfrac12(\tfrac n{\sigma^2}+\tfrac1{\tau_0^2})=-\tfrac1{2\tau_n^2}$, and the coefficient of $\mu$ is $\tfrac{n\bar y}{\sigma^2}+\tfrac{\mu_0}{\tau_0^2}=\mu_n/\tau_n^2$.
3. A quadratic $-\frac1{2\tau_n^2}\mu^2+\frac{\mu_n}{\tau_n^2}\mu$ equals $-\frac1{2\tau_n^2}(\mu-\mu_n)^2+\text{const}$. That is (3.10).
4. Write $\tilde y=\mu+\tilde\varepsilon$ with $\tilde\varepsilon\sim\mathcal N(0,\sigma^2)$ independent of $\mu\mid y$. A sum of independent Normals is Normal (F5 Theorem 3.22), with variances adding: (3.11). ∎

■ **What just happened.** **Precisions add** and the mean is precision-weighted: the data have precision $n/\sigma^2$, the prior $1/\tau_0^2$, and the posterior mean is the average of $\bar y$ and $\mu_0$ with those weights. Writing $n_0=\sigma^2/\tau_0^2$ gives $\mu_n=(n_0\mu_0+n\bar y)/(n_0+n)$, the same shape as (3.6) with $n_0$ imaginary observations. This is F5's Gaussian conditional in disguise: $(\mu,\bar y)$ is jointly Gaussian with $\mathrm{Cov}(\mu,\bar y)=\tau_0^2$ and $\mathrm{Var}(\bar y)=\tau_0^2+\sigma^2/n$, and F5 (3.19) gives $\mathbb E[\mu\mid\bar y]=\mu_0+\frac{\tau_0^2}{\tau_0^2+\sigma^2/n}(\bar y-\mu_0)$, which is $\mu_n$. The predictive variance is irreducible noise $\sigma^2$ plus parameter uncertainty $\tau_n^2\sim\sigma^2/n$. With $\tau_0\to\infty$ (flat prior) the posterior is $\mathcal N(\bar y,\sigma^2/n)$, whose 95% interval is ST-2's $z$ interval.

### 3.4 The Normal model with unknown variance

**Mean known, variance unknown.** The likelihood in $v=\sigma^2$ is $v^{-n/2}\exp(-\mathrm{SS}_\mu/(2v))$ with $\mathrm{SS}_\mu=\sum(y_i-\mu)^2$. Against a prior $v\sim\mathrm{IG}(a_0,b_0)\propto v^{-a_0-1}e^{-b_0/v}$ the exponents add:
$$v\mid y\sim\mathrm{IG}\big(a_0+\tfrac n2,\ b_0+\tfrac12\mathrm{SS}_\mu\big). \tag{3.12}$$
Read $2a_0$ as imaginary observations and $2b_0$ as their imaginary sum of squares. (In terms of the precision $1/v$ this is a Gamma prior, and BDA3 writes the same family as scaled inverse-$\chi^2$.)

**Both unknown.** Use the **Normal–Inverse-Gamma** (NIG) prior: $\sigma^2\sim\mathrm{IG}(\alpha_0,\beta_0)$ and $\mu\mid\sigma^2\sim\mathcal N(\mu_0,\sigma^2/\kappa_0)$. The prior on $\mu$ is scaled by $\sigma^2$: the prior is worth $\kappa_0$ observations *of the same noise level as the data*.

**Theorem 3.3 (NIG update).** With $S=\sum(y_i-\bar y)^2$ the posterior is $\mathrm{NIG}(\mu_n,\kappa_n,\alpha_n,\beta_n)$:
$$\kappa_n=\kappa_0+n,\quad \mu_n=\frac{\kappa_0\mu_0+n\bar y}{\kappa_n},\quad \alpha_n=\alpha_0+\frac n2,\quad \beta_n=\beta_0+\frac12\Big[S+\frac{\kappa_0 n}{\kappa_n}(\bar y-\mu_0)^2\Big]. \tag{3.13}$$

▶ **Try it first.** Write the joint density (likelihood × prior) as $(\sigma^2)^{-\text{something}}\exp\{-[\dots]/(2\sigma^2)\}$. Complete the square in $\mu$ inside the bracket. What is left over?

**Derivation.**
1. Likelihood: $(\sigma^2)^{-n/2}\exp\{-[S+n(\bar y-\mu)^2]/(2\sigma^2)\}$. Prior: $(\sigma^2)^{-\alpha_0-1}e^{-\beta_0/\sigma^2}\cdot(\sigma^2)^{-1/2}\exp\{-\kappa_0(\mu-\mu_0)^2/(2\sigma^2)\}$.
2. The product has exponent $-\frac1{2\sigma^2}\big[2\beta_0+S+n(\mu-\bar y)^2+\kappa_0(\mu-\mu_0)^2\big]$ and power of $\sigma^2$ equal to $-(\alpha_0+\frac n2)-1-\frac12$.
3. The two quadratics combine exactly as in Theorem 3.2 (precisions $n$ and $\kappa_0$, both scaled by $1/\sigma^2$): $n(\mu-\bar y)^2+\kappa_0(\mu-\mu_0)^2=\kappa_n(\mu-\mu_n)^2+\frac{n\kappa_0}{\kappa_n}(\bar y-\mu_0)^2$. (Check the constant: $\kappa_n\mu_n^2+\frac{n\kappa_0}{\kappa_n}(\bar y-\mu_0)^2=n\bar y^2+\kappa_0\mu_0^2$.)
4. The bracket is now $\kappa_n(\mu-\mu_n)^2+2\beta_n$ with $\beta_n$ as in (3.13), and the power of $\sigma^2$ is $-\alpha_n-1-\frac12$. That is the NIG kernel with the stated parameters. ∎

■ **What just happened.** The mean's imaginary count $\kappa_0$ and the variance's imaginary count $2\alpha_0$ each grow by the real count ($n$ and $n$). The extra term in $\beta_n$ is the **prior–data conflict**: if $\bar y$ is far from $\mu_0$, the posterior concludes the noise is larger.

**Lemma 3.4 (Normal scale mixtures are $t$).** If $x\mid v\sim\mathcal N(m,sv)$ and $v\sim\mathrm{IG}(a,b)$, then $x\sim t_{2a}\big(m,\ sb/a\big)$.

*Proof.* $\int(2\pi sv)^{-1/2}e^{-(x-m)^2/(2sv)}\frac{b^a}{\Gamma(a)}v^{-a-1}e^{-b/v}dv=\frac{b^a}{\Gamma(a)\sqrt{2\pi s}}\int v^{-(a+\frac12)-1}\exp\{-[b+\frac{(x-m)^2}{2s}]/v\}dv=\frac{b^a\,\Gamma(a+\frac12)}{\Gamma(a)\sqrt{2\pi s}}\big[b+\tfrac{(x-m)^2}{2s}\big]^{-(a+\frac12)}$, using $\int v^{-c-1}e^{-d/v}dv=\Gamma(c)d^{-c}$. As a function of $x$ this is $\propto[1+\frac{(x-m)^2}{2sb}]^{-(a+1/2)}$; match with the $t_\nu(m,c^2)$ density: $\nu=2a$ and $\nu c^2=2sb$. ∎

**Theorem 3.5 (Marginal and predictive are $t$).** Under the NIG posterior,
$$\mu\mid y\sim t_{2\alpha_n}\Big(\mu_n,\ \frac{\beta_n}{\alpha_n\kappa_n}\Big),\qquad \tilde y\mid y\sim t_{2\alpha_n}\Big(\mu_n,\ \frac{\beta_n(\kappa_n+1)}{\alpha_n\kappa_n}\Big). \tag{3.14, 3.15}$$

*Proof.* Given $\sigma^2$, $\mu\sim\mathcal N(\mu_n,\sigma^2/\kappa_n)$, and $\tilde y=\mu+\tilde\varepsilon\sim\mathcal N(\mu_n,\sigma^2(1+1/\kappa_n))$. Apply Lemma 3.4 with $a=\alpha_n$, $b=\beta_n$ and $s=1/\kappa_n$ or $s=1+1/\kappa_n$. ∎

**The uninformative case reproduces ST-2.** The density $p(\mu,\sigma^2)\propto1/\sigma^2$ ("flat in $\mu$ and in $\log\sigma$") is the formal limit $\kappa_0\to0$, $\alpha_0\to-\frac12$, $\beta_0\to0$ of the NIG density (power of $\sigma^2$: $-\alpha_0-\frac32\to-1$). It is improper, but the posterior is proper: $\kappa_n=n$, $\mu_n=\bar y$, $\alpha_n=\frac{n-1}2$, $\beta_n=\frac S2$, so (3.14) reads $\mu\mid y\sim t_{n-1}\big(\bar y,\ s^2/n\big)$ with $s^2=S/(n-1)$. The 95% credible interval is $\bar y\pm t_{n-1,0.975}\,s/\sqrt n$: **numerically identical to ST-2's $t$ interval**. Likewise $\sigma^2\mid y\sim\mathrm{IG}(\frac{n-1}2,\frac S2)$ gives ST-2's (3.4) interval for $\sigma^2$. The coincidence happens because $(\mu-\bar y)/(s/\sqrt n)$ has the same $t_{n-1}$ law whether the data are random with $\mu$ fixed (ST-2) or $\mu$ is random with the data fixed (here); that symmetry is special to this prior, which is unchanged by shifting and rescaling the data. The *meaning* differs: ST-2 says the procedure covers 95% of the time over repeated data; this says $P(\mu\in\text{interval}\mid y)=0.95$. The predictive (3.15) in this case has scale $s\sqrt{1+1/n}$, the same as the frequentist prediction interval of ST-3 (3.17).

### 3.5 The Poisson–Gamma model

**Theorem 3.6.** Let $y_i\mid\lambda\sim\mathrm{Poisson}(\lambda t_i)$ with known exposures $t_i$ (take $t_i=1$ for plain counts), and $\lambda\sim\mathrm{Gamma}(\alpha_0,\beta_0)$. Then
$$\lambda\mid y\sim\mathrm{Gamma}\Big(\alpha_0+\sum y_i,\ \beta_0+\sum t_i\Big), \tag{3.16}$$
and for a future count over exposure $\tilde t$, with $\alpha_n,\beta_n$ the posterior parameters,
$$P(\tilde y=k\mid y)=\frac{\Gamma(k+\alpha_n)}{k!\,\Gamma(\alpha_n)}\Big(\frac{\beta_n}{\beta_n+\tilde t}\Big)^{\alpha_n}\Big(\frac{\tilde t}{\beta_n+\tilde t}\Big)^k, \tag{3.17}$$
a Negative Binomial with mean $\alpha_n\tilde t/\beta_n$ and variance $\frac{\alpha_n\tilde t}{\beta_n}\big(1+\frac{\tilde t}{\beta_n}\big)$.

*Derivation.* (i) The likelihood is $\prod(\lambda t_i)^{y_i}e^{-\lambda t_i}/y_i!\propto\lambda^{\sum y_i}e^{-\lambda\sum t_i}$; times $\lambda^{\alpha_0-1}e^{-\beta_0\lambda}$ the exponents add: (3.16). (ii) $\int\frac{(\lambda\tilde t)^ke^{-\lambda\tilde t}}{k!}\frac{\beta_n^{\alpha_n}\lambda^{\alpha_n-1}e^{-\beta_n\lambda}}{\Gamma(\alpha_n)}d\lambda=\frac{\tilde t^k\beta_n^{\alpha_n}}{k!\,\Gamma(\alpha_n)}\cdot\frac{\Gamma(k+\alpha_n)}{(\beta_n+\tilde t)^{k+\alpha_n}}$, which is (3.17). (iii) Mean and variance by Adam's and Eve's laws: $\mathbb E[\lambda\tilde t]$ and $\mathbb E[\lambda\tilde t]+\mathrm{Var}(\lambda\tilde t)$. ∎

■ **What just happened.** This is ST-4 Theorem 3.8 (3.31) with $r=\alpha_n$ and mean $\alpha_n\tilde t/\beta_n$: the same algebra, a different meaning. There the Gamma modeled heterogeneity *across units*; here it is your belief about one rate. "Overdispersion" (variance/mean $=1+\tilde t/\beta_n>1$) is parameter uncertainty. The evidence is the ratio of Gamma normalizers: $p(y)=\prod\frac{t_i^{y_i}}{y_i!}\cdot\frac{\Gamma(\alpha_n)}{\Gamma(\alpha_0)}\frac{\beta_0^{\alpha_0}}{\beta_n^{\alpha_n}}$.

### 3.6 Conjugacy for every exponential family

Take ST-1's family (3.21), $p(y\mid\eta)=h(y)\exp\{\eta^\top T(y)-A(\eta)\}$, so $n$ independent observations have likelihood $\prod h(y_i)\cdot\exp\{\eta^\top\sum T(y_i)-nA(\eta)\}$. Define the **conjugate prior** with hyperparameters $(\chi,\nu)$:
$$p(\eta\mid\chi,\nu)=\frac{1}{Z(\chi,\nu)}\exp\{\eta^\top\chi-\nu A(\eta)\},\qquad Z(\chi,\nu)=\int\exp\{\eta^\top\chi-\nu A(\eta)\}d\eta. \tag{3.18}$$

**Theorem 3.7 (Conjugate update).** (a) The posterior is $p(\eta\mid\chi',\nu')$ with
$$\chi'=\chi+\sum_i T(y_i),\qquad \nu'=\nu+n. \tag{3.19}$$
(b) The evidence is $p(y)=\prod h(y_i)\cdot Z(\chi',\nu')/Z(\chi,\nu)$. (c) If the prior density vanishes at the boundary of the natural parameter space, then $\mathbb E[\nabla A(\eta)]=\chi/\nu$; since $\nabla A(\eta)=\mathbb E_\eta[T]$ is the mean parameter (ST-1 (3.22)),
$$\mathbb E\big[\mathbb E_\eta T\mid y\big]=\frac{\chi+\sum T(y_i)}{\nu+n}=\frac{\nu}{\nu+n}\cdot\frac\chi\nu+\frac n{\nu+n}\cdot\bar T. \tag{3.20}$$

▶ **Try it first.** (a) and (b) are one multiplication. For (c), differentiate the *unnormalized* prior with respect to $\eta_j$ and integrate over $\eta$.

**Derivation.**
1. Prior × likelihood $\propto\exp\{\eta^\top(\chi+\sum T)-(\nu+n)A(\eta)\}$: same form, new hyperparameters. The normalizer is $Z(\chi',\nu')$, which gives (a). The evidence is the ratio of normalizers by §3.1 (ii), which is (b).
2. $\partial_j\exp\{\eta^\top\chi-\nu A\}=(\chi_j-\nu\,\partial_jA)\exp\{\eta^\top\chi-\nu A\}$. Integrate over $\eta$: the left side is a boundary term, zero by assumption. So $0=Z\,(\chi_j-\nu\,\mathbb E[\partial_jA])$, hence $\mathbb E[\partial_jA]=\chi_j/\nu$.
3. Apply step 2 to the posterior, whose hyperparameters are $(\chi',\nu')$; the right side of (3.20) rearranges the fraction. ∎

■ **What just happened.** $\chi$ is a **prior sufficient statistic** and $\nu$ a **prior sample size**: the prior is the likelihood of $\nu$ imaginary observations whose $T$ sums to $\chi$. Updating adds real to imaginary. The posterior mean of the mean parameter is *always* a convex combination of prior mean $\chi/\nu$ and sample mean $\bar T$ with weights $\nu$ and $n$. Every shrinkage formula so far is this one:

| Family | $T$ | $\eta$, $A(\eta)$ | Prior in the usual coordinates | $(\chi,\nu)$ |
|---|---|---|---|---|
| Bernoulli | $y$ | $\mathrm{logit}\,p$, $\log(1+e^\eta)$ | $p^{\chi}(1-p)^{\nu-\chi}\cdot\frac{1}{p(1-p)}$: $\mathrm{Beta}(\chi,\nu-\chi)$ | $(\alpha_0,\ \alpha_0+\beta_0)$ |
| Poisson | $y$ | $\log\lambda$, $e^\eta$ | $\lambda^{\chi}e^{-\nu\lambda}\cdot\frac1\lambda$: $\mathrm{Gamma}(\chi,\nu)$ | $(\alpha_0,\ \beta_0)$ |
| Normal, $\sigma^2$ known | $y$ | $\mu/\sigma^2$, $\sigma^2\eta^2/2$ | $\exp\{\mu\chi/\sigma^2-\nu\mu^2/(2\sigma^2)\}$: $\mathcal N(\chi/\nu,\sigma^2/\nu)$ | $(\nu\mu_0,\ \sigma^2/\tau_0^2)$ |

The factors $1/(p(1-p))$ and $1/\lambda$ are the Jacobians $d\eta/dp$ and $d\eta/d\lambda$. Check the first row against (3.4): $\chi+\sum y=\alpha_0+y$ and $\nu+n-(\chi+\sum y)=\beta_0+n-y$. The third row's $\nu=\sigma^2/\tau_0^2$ is the $n_0$ of §3.3. The NIG family of §3.4 is the two-parameter Normal's conjugate family with one extra degree of freedom: the exponential-family form ties $\kappa_0=\nu$ and $\alpha_0=(\nu+3)/2$, and (3.13) preserves that tie, so letting $\kappa_0$ and $\alpha_0$ vary independently gives a larger conjugate family. Conjugacy exists for any exponential family, but the posterior is a named, easy distribution only for a few; for the rest (logistic regression, mixtures, hierarchical models) you need §3.11's Laplace approximation or ST-6.

### 3.7 Priors: informative, weakly informative, "noninformative"

**Flat is not a property of a distribution.** Let $p\sim\mathrm{Uniform}(0,1)$ and $\eta=\log\frac p{1-p}$, so $p=\sigma(\eta)$ and $dp/d\eta=p(1-p)$. By the change-of-variables rule (F5 Theorem 3.6),
$$f_\eta(\eta)=f_p(\sigma(\eta))\,\Big|\frac{dp}{d\eta}\Big|=\sigma(\eta)(1-\sigma(\eta))=\frac{e^\eta}{(1+e^\eta)^2}, \tag{3.21}$$
the logistic density: almost all its mass lies in $|\eta|<4$. "No opinion about $p$" is an opinion about the log-odds. Reverse the question: flat in $\eta$ means density $\propto d\eta/dp=1/(p(1-p))$ in $p$, the improper $\mathrm{Beta}(0,0)$ (Haldane). Three candidate "ignorance" priors for a probability are Beta(0,0), Beta($\frac12,\frac12$) and Beta(1,1). After ten failures in ten trials the last two give posterior means $0.0455$ and $0.0833$ (Example 4.5); the first gives an improper posterior. There is no neutral choice, only choices that behave well under particular transformations.

**Jeffreys' prior.** Define
$$p_J(\theta)\propto\sqrt{I_1(\theta)}\qquad(\text{vector }\theta:\ \sqrt{\det I_1(\theta)}). \tag{3.22}$$

**Theorem 3.8 (Invariance).** If $\phi=h(\theta)$ is a smooth one-to-one reparametrization, the Jeffreys prior computed in $\phi$ equals the Jeffreys prior in $\theta$ transformed as a density.

▶ **Try it first.** ST-1 (3.13) says how $I_1$ changes under reparametrization. Take square roots and compare with the change-of-variables rule for densities.

*Proof.* ST-1 (3.13), scalar case: $I_\phi(\phi)=I_\theta(\theta(\phi))\,(d\theta/d\phi)^2$. Take square roots: $\sqrt{I_\phi}=\sqrt{I_\theta}\,|d\theta/d\phi|$. The density transformation rule says $p_\phi(\phi)=p_\theta(\theta(\phi))\,|d\theta/d\phi|$. The two right sides have the same form, so $p_J$ in $\phi$ is $p_J$ in $\theta$ pushed through $h$. (Vector case: $\det I_\phi=\det(J^\top I_\theta J)=\det I_\theta\,(\det J)^2$.) ∎ $\qquad(3.23)$

■ **What just happened.** Fisher information transforms like the square of a Jacobian, which is exactly what makes $\sqrt I$ transform like a density. The rule is "uniform in arc length of the Fisher metric". Examples, from ST-1's table:

| Model | $I_1$ | Jeffreys prior | Posterior |
|---|---|---|---|
| Bernoulli($p$) | $\frac1{p(1-p)}$ | $p^{-1/2}(1-p)^{-1/2}$ = Beta($\frac12,\frac12$) | Beta($y+\frac12,\ n-y+\frac12$) |
| Poisson($\lambda$) | $\frac1\lambda$ | $\lambda^{-1/2}$ = Gamma($\frac12,0$), improper | Gamma($\sum y+\frac12,\ n$) |
| $\mathcal N(\mu,\sigma^2)$, $\sigma$ known | $\frac1{\sigma^2}$ | constant (flat on $\mathbb R$), improper | $\mathcal N(\bar y,\sigma^2/n)$ |
| $\mathcal N(\mu,\sigma^2)$, $\mu$ known, scale $\sigma$ | $\frac2{\sigma^2}$ | $\frac1\sigma$ (flat in $\log\sigma$); in $v=\sigma^2$: $\frac1v$ | IG($\frac n2,\ \frac12\mathrm{SS}_\mu$) |

For the last row, $I_\sigma=2/\sigma^2$ is ST-1's (3.12) transformed by (3.13) ($I_v=\frac1{2v^2}$, $dv/d\sigma=2\sigma$). In the Bernoulli row, the $\eta$-scale prior is $\frac1{2\pi\cosh(\eta/2)}$ (the project asserts this equals $\sqrt{I_\eta}/\pi$ with $I_\eta=p(1-p)$). The two Normal parameters together: the joint Jeffreys prior $\sqrt{\det I_1}$ is $\propto\sigma^{-2}$ in $(\mu,\sigma)$, but the standard reference prior treats location and scale separately and uses $p(\mu,\sigma^2)\propto\sigma^{-2}$, i.e. $\sigma^{-1}$ in $(\mu,\sigma)$, which is what gave ST-2's $t_{n-1}$ in §3.4. The joint Jeffreys prior gives $t_n$ (narrower: half-width $0.0279$ against $0.0298$ for the Example 4.3 data). "The" noninformative prior is a choice.

**Improper priors** do not integrate to 1. They are acceptable when the posterior is proper, which you must check: Haldane's Beta(0,0) with $y=0$ gives Beta($0,n$), not a distribution. They are *not* acceptable in Bayes factors (§3.12): the unknown constant multiplying an improper prior multiplies the evidence and does not cancel from a ratio against a different model.

**Weakly informative priors** are the modern default (Gelman et al.): proper distributions that rule out absurd values and leave plausible ones open. Check them with a **prior predictive check**: draw $\theta\sim p(\theta)$, then data $\sim p(y\mid\theta)$, and look at the datasets. Example: a "vague" $\mathcal N(0,10^2)$ prior on the log-odds says $P(p<0.01\text{ or }p>0.99)=2\Phi(-4.595/10)=0.646$: two thirds of your prior belief is that the coin is nearly deterministic. A $\mathcal N(0,1.5^2)$ prior puts $0.0022$ there. Maximum-entropy priors (ST-1 Theorem 3.27, SM-4) and Bernardo's reference priors are other principled defaults (BDA3 ch. 2).

### 3.8 Credible intervals, calibration, and Bernstein–von Mises

A set $C(y)$ is a $100(1-\alpha)\%$ **credible set** if $P(\theta\in C(y)\mid y)=1-\alpha$. Two standard choices: **equal-tailed** (the $\alpha/2$ and $1-\alpha/2$ posterior quantiles) and **highest posterior density** (HPD), $C=\{\theta:p(\theta\mid y)\ge c\}$ with $c$ set so the mass is $1-\alpha$.

**Theorem 3.9 (HPD is shortest).** Among all sets with posterior mass $1-\alpha$, the HPD set has the smallest volume.

*Proof.* Let $H=\{p\ge c\}$ and $C$ be any set of equal mass. Removing the overlap, $\int_{H\setminus C}p=\int_{C\setminus H}p$. On $H\setminus C$ the density is $\ge c$ and on $C\setminus H$ it is $<c$, so $c\,|H\setminus C|\le\int_{H\setminus C}p=\int_{C\setminus H}p\le c\,|C\setminus H|$. Hence $|H|\le|C|$. ∎

The HPD interval has equal density at its ends (a unimodal posterior), is not invariant under reparametrization (density changes with the Jacobian), and for a monotone posterior (a pile-up at 0) is one-sided. Equal-tailed intervals map to equal-tailed intervals under monotone maps.

**Interpretation.** $P(\theta\in C(y)\mid y)=1-\alpha$ is a statement about $\theta$ given the data you have: what ST-2 insisted a confidence interval is *not*. The price is the prior. What guarantee survives the dependence on the prior?

**Theorem 3.10 (Bayesian calibration).** If $\theta\sim p(\theta)$, $y\mid\theta\sim p(y\mid\theta)$ and $C(y)$ is a $(1-\alpha)$ credible set computed from that same prior and model, then $P(\theta\in C(y))=1-\alpha$, the probability taken over both $\theta$ and $y$.

*Proof.* $P(\theta\in C(y))=\mathbb E_y\big[P(\theta\in C(y)\mid y)\big]=\mathbb E_y[1-\alpha]=1-\alpha$ by Adam's law (F5 Theorem 3.26). ∎ $\qquad(3.24)$

This is coverage *averaged over the prior*. At a fixed $\theta_0$ nothing is promised: the project computes exact coverage of the Beta–Binomial interval at $\theta_0=0.1$. With a uniform prior it is $0.930$ at $n=10$, $0.936$ at $n=100$, $0.949$ at $n=1000$. With a Beta(10,10) prior centered in the wrong place it is $0.000$, $0.451$, $0.883$: a strong prior's damage fades slowly (its bias is $O(1/n)$ but the posterior sd is $O(n^{-1/2})$, so their ratio shrinks only like $n^{-1/2}$). At $\theta_0=0.5$ the same prior *over*-covers (1.000, 0.965, 0.954). And averaging the uniform-prior coverage over $\theta_0\sim\mathrm{Uniform}(0,1)$ at $n=20$ gives $0.951$, Theorem 3.10 in action.

**Theorem 3.11 (Bernstein–von Mises; stated, heuristic derivation).** Under ST-1's regularity conditions, a prior with positive continuous density at the truth $\theta_0$, and a correctly specified model, the posterior is asymptotically $\mathcal N\big(\hat\theta,\ 1/(nI_1(\theta_0))\big)$ in total variation, in probability. $\qquad(3.25)$

*Heuristic derivation* (full proof: van der Vaart, *Asymptotic Statistics*, ch. 10). Write $\log p(\theta\mid y)=\ell_n(\theta)+\log p(\theta)+\text{const}$. (1) Expand $\ell_n$ about the MLE: $\ell_n(\theta)\approx\ell_n(\hat\theta)-\frac12 n\hat J(\theta-\hat\theta)^2$, where $\hat J=-\ell_n''(\hat\theta)/n\to I_1(\theta_0)$ by the law of large numbers (ST-1 Theorem 3.19). (2) The posterior lives where $|\theta-\hat\theta|\sim n^{-1/2}$; there $\log p(\theta)$ varies by $O(n^{-1/2})$, negligible against the $O(1)$ quadratic. (3) So $p(\theta\mid y)\propto\exp\{-\frac12nI_1(\theta-\hat\theta)^2\}$. ∎

■ **What just happened.** ST-1's asymptotic normality ($\sqrt n(\hat\theta-\theta_0)\to\mathcal N(0,1/I_1)$) read from the parameter's side: the same Normal, but *centered at the data and describing belief*. Consequences: credible intervals are asymptotically Wald intervals $\hat\theta\pm z_{\alpha/2}/\sqrt{nI_1}$ (ST-2 §3.6), their frequentist coverage tends to $1-\alpha$, the posterior mean tends to the MLE, and the prior is forgotten at rate $1/n$ (compare $w\to1$ in (3.6)). It fails, or applies late, for small $n$, strong priors, parameters on a boundary, non-identifiable parameters, a dimension that grows with $n$, and a misspecified model (the posterior then concentrates at the KL projection, ST-1 Theorem 3.28, with the wrong width).

### 3.9 Point estimates and decisions

Choose an action $a$ (an estimate) to minimize the **posterior expected loss** $\rho(a\mid y)=\mathbb E[L(a,\theta)\mid y]$. The frequentist **risk** of a rule $\delta(y)$ is $R(\theta,\delta)=\mathbb E_{y\mid\theta}L(\delta(y),\theta)$, and its **Bayes risk** is
$$r(\delta)=\mathbb E_\theta R(\theta,\delta)=\mathbb E_y\,\rho(\delta(y)\mid y). \tag{3.26}$$
(Swap the order of integration: the Bayes-risk-minimizing rule is the rule that minimizes $\rho$ at every $y$.)

**Theorem 3.12.** (a) Squared loss: the minimizer is the posterior **mean** $m=\mathbb E[\theta\mid y]$, and the minimum is $\mathrm{Var}(\theta\mid y)$, so the Bayes risk of the posterior mean is $\mathbb E_y[\mathrm{Var}(\theta\mid y)]$. (b) Absolute loss: the posterior **median**. (c) 0–1 loss with tolerance $\varepsilon\to0$: the posterior **mode**.

*Proof.* (a) $\rho(a)=\mathbb E[(\theta-m+m-a)^2]=\mathrm{Var}(\theta\mid y)+(m-a)^2$, since the cross term has factor $\mathbb E[\theta-m]=0$. (b) For a continuous posterior with CDF $F$, $\rho(a)=\int_{-\infty}^a(a-\theta)p+\int_a^\infty(\theta-a)p$ and $\rho'(a)=F(a)-(1-F(a))$, which vanishes at $F(a)=\frac12$; $\rho''=2p(a)\ge0$. (c) $L=\mathbb 1\{|\theta-a|>\varepsilon\}$ gives $\rho(a)=1-P(|\theta-a|\le\varepsilon\mid y)\approx1-2\varepsilon\,p(a\mid y)$, minimized at the maximum of the density. ∎

For the Example 4.1 posterior Beta(9,5): mean $0.6429$, median $0.6498$, mode $0.6667$. Three questions, three answers (the project minimizes each loss numerically).

**MAP is penalized maximum likelihood.** The mode maximizes $\log p(y\mid\theta)+\log p(\theta)$. For $y\mid\beta\sim\mathcal N(X\beta,\sigma^2I)$:

*Ridge.* With $\beta\sim\mathcal N(0,\tau^2I)$, $-\log p(\beta\mid y)=\frac1{2\sigma^2}\|y-X\beta\|^2+\frac1{2\tau^2}\|\beta\|^2+c=\frac1{2\sigma^2}\big[\|y-X\beta\|^2+\lambda\|\beta\|^2\big]+c$ with
$$\lambda=\sigma^2/\tau^2. \tag{3.27}$$
This is ST-3's ridge (3.38). A tight prior (small $\tau^2$) is a large penalty.

*Lasso.* With independent Laplace priors $p(\beta_j)=e^{-|\beta_j|/b}/(2b)$, $-\log p(\beta\mid y)=\frac1{\sigma^2}\big[\frac12\|y-X\beta\|^2+\frac{\sigma^2}b\|\beta\|_1\big]+c$, which is ST-3's (3.43) with
$$\lambda=\sigma^2/b. \tag{3.28}$$

■ **What just happened, and the catch.** The MAP *is* the penalized fit, and ST-3's ridge-versus-lasso contrast is a contrast of priors. But the mode ignores the shape of the posterior and is not invariant: Beta(9,5) has mode $2/3$ in $p$, while its density on the logit scale peaks at $p=9/14$. The lasso's exact zeros are a property of the mode: for an orthonormal design with $\sigma^2=1$, $z=0.3$ and prior scale $b=2$ ($\lambda=\frac12$) the MAP is exactly $0$ but the posterior mean is $0.2047$; the posterior puts no point mass at zero. (Park and Casella's "Bayesian lasso" makes this point.) MAP estimates are computed with a prior but are not Bayesian in spirit.

### 3.10 Posterior predictive checks

The posterior predictive (3.3) is also a simulator: draw $\theta^{(s)}\sim p(\theta\mid y)$, then replicated data $y^{\rm rep(s)}\sim p(y\mid\theta^{(s)})$. If the model is adequate, $y^{\rm rep}$ should look like $y$. Pick a **test quantity** $T(y)$ for a feature you care about, and compute the **posterior predictive $p$-value**
$$p_B=P\big(T(y^{\rm rep})\ge T(y)\mid y\big). \tag{3.29}$$
Values near 0 or 1 flag a feature the model cannot reproduce. It is a diagnostic, not a calibrated test (the data are used twice, so $p_B$ is conservative), and the statistic should not be one the fit forces to match: a Normal fit by $\bar y,s^2$ cannot fail a check on those. Example 4.11 fits a Normal to skewed data. The *prior* predictive check of §3.7 is the same procedure with $p(\theta)$ in place of $p(\theta\mid y)$.

### 3.11 Bayesian linear regression

**Theorem 3.13 (Known noise variance).** Let $y=X\beta+\varepsilon$, $\varepsilon\sim\mathcal N(0,\sigma^2I_n)$ with $\sigma^2$ known, $X$ an $n\times p$ design, and $\beta\sim\mathcal N(\mu_0,\Sigma_0)$. Then $\beta\mid y\sim\mathcal N(\mu_n,\Sigma_n)$ with
$$\Sigma_n^{-1}=\Sigma_0^{-1}+\frac{X^\top X}{\sigma^2},\qquad \mu_n=\Sigma_n\Big(\Sigma_0^{-1}\mu_0+\frac{X^\top y}{\sigma^2}\Big). \tag{3.30}$$

▶ **Try it first.** Write the log posterior as a quadratic in $\beta$. You are repeating §3.3 with matrices.

**Derivation.**
1. $\log p(\beta\mid y)=-\frac1{2\sigma^2}\|y-X\beta\|^2-\frac12(\beta-\mu_0)^\top\Sigma_0^{-1}(\beta-\mu_0)+\text{const}$.
2. Expand in $\beta$. The quadratic part is $-\frac12\beta^\top\big(\Sigma_0^{-1}+X^\top X/\sigma^2\big)\beta=-\frac12\beta^\top\Sigma_n^{-1}\beta$; the linear part is $\beta^\top\big(\Sigma_0^{-1}\mu_0+X^\top y/\sigma^2\big)=\beta^\top\Sigma_n^{-1}\mu_n$.
3. $-\frac12\beta^\top\Sigma_n^{-1}\beta+\beta^\top\Sigma_n^{-1}\mu_n=-\frac12(\beta-\mu_n)^\top\Sigma_n^{-1}(\beta-\mu_n)+\text{const}$. A Gaussian kernel with this mean and precision. ∎

■ **What just happened.** Precision matrices add and the mean is precision-weighted, as in (3.10). This is the **fourth** time you have completed this square (LA-5's normal equations, F5 Theorem 3.24, §3.3, and now). It is also F5 (3.19) applied to the jointly Gaussian pair $(\beta,y)$ with $\mathrm{Cov}(\beta,y)=\Sigma_0X^\top$ and $\mathrm{Cov}(y)=X\Sigma_0X^\top+\sigma^2I$: $\mathbb E[\beta\mid y]=\mu_0+\Sigma_0X^\top(X\Sigma_0X^\top+\sigma^2I)^{-1}(y-X\mu_0)$, which equals $\mu_n$ by the matrix inversion lemma (the project asserts it). Two special cases:

*Ridge.* $\mu_0=0$, $\Sigma_0=\tau^2I$ gives
$$\mu_n=(X^\top X+\lambda I)^{-1}X^\top y,\quad \Sigma_n=\sigma^2(X^\top X+\lambda I)^{-1},\quad\lambda=\sigma^2/\tau^2, \tag{3.31}$$
ST-3's (3.39) for the mean, and now a covariance too: ST-3 could only say ridge shrinks; here the posterior also says how uncertain the shrunk estimate is. Mean and mode coincide because the posterior is Gaussian.

*Flat prior.* $\Sigma_0^{-1}\to0$ gives $\mu_n\to\hat\beta_{\rm OLS}$ and $\Sigma_n\to\sigma^2(X^\top X)^{-1}$, ST-3 Theorem 3.5(b).

**Predictive distribution.** For a new input $x_*$, $\tilde y=x_*^\top\beta+\tilde\varepsilon$. Given $y$, $x_*^\top\beta\sim\mathcal N(x_*^\top\mu_n,\ x_*^\top\Sigma_nx_*)$ (F5 Theorem 3.22) and $\tilde\varepsilon\sim\mathcal N(0,\sigma^2)$ is independent, so
$$\tilde y\mid y\sim\mathcal N\big(x_*^\top\mu_n,\ \sigma^2+x_*^\top\Sigma_nx_*\big). \tag{3.32}$$
With a flat prior $x_*^\top\Sigma_nx_*=\sigma^2h_*$ with $h_*=x_*^\top(X^\top X)^{-1}x_*$, the leverage of ST-3 (3.16): the predictive variance is $\sigma^2(1+h_*)$, ST-3's prediction interval (3.17) with $\sigma^2$ known. For simple regression $h_*=\frac1n+\frac{(x_*-\bar x)^2}{S_{xx}}$: **the variance grows quadratically with distance from the data's center.** The band is narrowest where you have data and flares outward. This is the picture of ML-4. The function value $f(x)=\phi(x)^\top\beta$ has prior covariance $k(x,x')=\phi(x)^\top\Sigma_0\phi(x')$ and posterior covariance $\phi(x)^\top\Sigma_n\phi(x')$: Bayesian linear regression is Gaussian-process regression with a finite feature map, and ML-4 lets the number of features grow without bound.

**Unknown $\sigma^2$.** Put $\beta\mid\sigma^2\sim\mathcal N(\mu_0,\sigma^2V_0)$ and $\sigma^2\sim\mathrm{IG}(a_0,b_0)$. Then $\beta,\sigma^2\mid y$ is NIG with
$$V_n=(V_0^{-1}+X^\top X)^{-1},\ \ \mu_n=V_n(V_0^{-1}\mu_0+X^\top y),\ \ a_n=a_0+\tfrac n2,\ \ b_n=b_0+\tfrac12\big(y^\top y+\mu_0^\top V_0^{-1}\mu_0-\mu_n^\top V_n^{-1}\mu_n\big),$$
the matrix form of Theorem 3.3. The derivation is the same two steps: complete the square in $\beta$ (the identity $\|y-X\beta\|^2+(\beta-\mu_0)^\top V_0^{-1}(\beta-\mu_0)=(\beta-\mu_n)^\top V_n^{-1}(\beta-\mu_n)+[y^\top y+\mu_0^\top V_0^{-1}\mu_0-\mu_n^\top V_n^{-1}\mu_n]$ is a direct expansion), then collect powers of $\sigma^2$: $(\sigma^2)^{-n/2-p/2-a_0-1}\exp\{-[2b_0+\text{bracket}]/(2\sigma^2)\}$. With $X=\mathbf 1$ and $V_0=1/\kappa_0$ it reduces to (3.13). By Lemma 3.4 the predictive is a $t$:
$$\tilde y\mid y\sim t_{2a_n}\Big(x_*^\top\mu_n,\ \frac{b_n}{a_n}\big(1+x_*^\top V_nx_*\big)\Big). \tag{3.33}$$
**Zellner's $g$-prior** takes $\mu_0=0$, $V_0=g\,(X^\top X)^{-1}$: the prior covariance has the shape of the data's, so $\mu_n=\frac g{1+g}\hat\beta_{\rm OLS}$, uniform shrinkage with one scalar to choose ($g=n$ is the "unit information" prior).

**Logistic regression has no conjugate prior.** With $\beta\sim\mathcal N(\mu_0,\Sigma_0)$ the log posterior is
$$f(\beta)=\sum_i\big[y_ix_i^\top\beta-\log(1+e^{x_i^\top\beta})\big]-\tfrac12(\beta-\mu_0)^\top\Sigma_0^{-1}(\beta-\mu_0)+\text{const}.$$
Both terms are concave (ST-4 §3.3; the prior term is a negative-definite quadratic), so $f$ has a unique maximizer $\hat\beta$, the MAP, found by Newton's method, which is ST-4's IRLS with a ridge term (3.13). It exists even under complete separation (ST-4 Theorem 3.4), where the MLE does not. The **Laplace approximation** expands $f$ to second order at $\hat\beta$ (LA-14 Theorem 3.13): the gradient vanishes there, and $\nabla^2f(\hat\beta)=-A$ with $A=X^\top\hat WX+\Sigma_0^{-1}$, $\hat W=\mathrm{diag}\,\hat p_i(1-\hat p_i)$ (ST-4 (3.6)). So
$$p(\beta\mid y)\approx\mathcal N\big(\hat\beta,\ A^{-1}\big),\qquad p(y)\approx e^{f(\hat\beta)}(2\pi)^{p/2}|A|^{-1/2}, \tag{3.34, 3.35}$$
where $f$ here includes the normalizing constants, so that $e^{f}$ is the joint density. This is a local fit at the peak: it cannot see skewness. Example 4.8 shows both its success (the evidence) and its failure (the mean). The same expansion underlies ML-4's GP classification (Rasmussen & Williams ch. 3). **Hierarchical priors** (an unknown $\tau^2$ in $\beta\sim\mathcal N(0,\tau^2I)$, learned across groups) are developed in ST-6.

### 3.12 Model comparison

**The evidence.** For a model $M$ with parameters $\theta$,
$$p(y\mid M)=\int p(y\mid\theta,M)\,p(\theta\mid M)\,d\theta, \tag{3.36}$$
the likelihood *averaged over the prior*. By the chain rule of probability (F4 Theorem 3.5) $p(y)=\prod_ip(y_i\mid y_{<i})$, the product of one-step-ahead predictive densities (3.3): the evidence is the score a model earns for predicting each observation from the earlier ones, with no held-out set. The project verifies this for the NIG regression model of §3.11. The **Bayes factor** is $\mathrm{BF}_{10}=p(y\mid M_1)/p(y\mid M_0)$ and, by Bayes' theorem in odds form (F4 Theorem 3.6), posterior odds = $\mathrm{BF}_{10}\times$ prior odds. Verbal scales for $\mathrm{BF}$ are conventions, not theorems.

**The Occam factor.** Suppose the posterior is near-Normal with width $\sigma_{\rm post}$ at its peak $\hat\theta$ (Laplace, one parameter). Then
$$p(y\mid M)\approx\underbrace{p(y\mid\hat\theta)}_{\text{best fit}}\times\underbrace{p(\hat\theta)\sqrt{2\pi}\,\sigma_{\rm post}}_{\text{Occam factor}}\ =\ p(y\mid\hat\theta)\cdot\frac{\sqrt{2\pi}\,\sigma_{\rm post}}{\Delta}\ \ \text{(prior flat on width }\Delta). \tag{3.37}$$
The Occam factor is the fraction of the prior's range that survives the data: posterior width over prior width. It is at most about 1, it shrinks with each parameter the data pin down, and a diffuse prior makes it smaller. MacKay's point: the penalty falls out of the integral; nobody adds it. For $d$ parameters the factors multiply. In the linear-Gaussian model this is exact, not approximate.

**Theorem 3.14 (Evidence of Bayesian linear regression).** For $y=\Phi\beta+\varepsilon$, $\varepsilon\sim\mathcal N(0,\sigma^2I_n)$, $\beta\sim\mathcal N(0,\alpha^{-1}I_d)$, let $A=\alpha I+\Phi^\top\Phi/\sigma^2$ and $m=A^{-1}\Phi^\top y/\sigma^2$. Then
$$\log p(y)=\underbrace{-\tfrac n2\log(2\pi\sigma^2)-\tfrac1{2\sigma^2}\|y-\Phi m\|^2}_{\text{fit}}\ \underbrace{-\tfrac\alpha2\|m\|^2+\tfrac d2\log\alpha-\tfrac12\log|A|}_{\text{log Occam factor}}, \tag{3.38}$$
and equivalently $y\sim\mathcal N(0,\ \sigma^2I+\alpha^{-1}\Phi\Phi^\top)$.

*Proof.* $g(\beta)=\log p(y\mid\beta)+\log p(\beta)$ is an exact quadratic in $\beta$ with Hessian $-A$ and maximizer $m$ (Theorem 3.13), so $g(\beta)=g(m)-\frac12(\beta-m)^\top A(\beta-m)$ with no remainder. Integrating $e^{g}$ gives $e^{g(m)}(2\pi)^{d/2}|A|^{-1/2}$. Write out $g(m)=-\frac n2\log(2\pi\sigma^2)-\frac1{2\sigma^2}\|y-\Phi m\|^2+\frac d2\log\frac\alpha{2\pi}-\frac\alpha2\|m\|^2$ and combine the $(2\pi)^{d/2}$: (3.38). The second form: $y=\Phi\beta+\varepsilon$ is a linear map of Gaussians (F5 Theorem 3.22). ∎

For unknown $\sigma^2$ with the NIG prior, the evidence is $\log p(y)=\log\frac{\Gamma(a_n)}{\Gamma(a_0)}+\frac12\log\frac{|V_n|}{|V_0|}+a_0\log b_0-a_n\log b_n-\frac n2\log2\pi$, again a ratio of normalizers (§3.1 (ii)). $\qquad(3.39)$

**The Jeffreys–Lindley paradox.** Let $\bar x\mid\mu\sim\mathcal N(\mu,\sigma^2/n)$, $H_0:\mu=0$, $H_1:\mu\sim\mathcal N(0,\tau^2)$. Under $H_0$, $\bar x\sim\mathcal N(0,\sigma^2/n)$; under $H_1$ it is $\mathcal N(0,\sigma^2/n+\tau^2)$ (sum of independent Normals). With $z=\bar x\sqrt n/\sigma$ and $u=n\tau^2/\sigma^2$, the ratio of the two densities at $\bar x$ is
$$\mathrm{BF}_{10}=\sqrt{\frac{\sigma^2/n}{\sigma^2/n+\tau^2}}\ \exp\Big\{\frac{\bar x^2}2\Big[\frac1{\sigma^2/n}-\frac1{\sigma^2/n+\tau^2}\Big]\Big\}=(1+u)^{-1/2}\exp\Big\{\frac{z^2}2\cdot\frac u{1+u}\Big\}. \tag{3.40}$$
Hold $z$, and with it the two-sided $p$-value, fixed and let $n\to\infty$. Then $u\to\infty$ and $\mathrm{BF}_{10}\approx u^{-1/2}e^{z^2/2}\to0$: evidence for $H_0$ grows like $\sqrt n$ while the $p$-value stays put. With $n=10^4$, $\bar x=0.02$, $\sigma=\tau=1$: $z=2$, $p=0.0455$, $\mathrm{BF}_{10}=0.0739$, so $H_0$ is favored $13.5:1$. Why: $H_1$ spreads its bets over $|\mu|\lesssim1$, so it predicted $\bar x$ to within about $\pm1$, while $H_0$ predicted $\pm0.01$; an observation at $0.02$ is 2 standard errors from $H_0$'s prediction but sits in a region where $H_1$'s density is a hundred times thinner than $H_0$'s peak. Estimation under $H_1$ says something else: $\mu\mid y\sim\mathcal N(0.0200,0.0100^2)$, interval $(0.0004,0.0396)$ excludes 0. All three statements are true. The $p$-value says the data are unusual under $H_0$; the interval says $\mu$ is probably positive and tiny; the Bayes factor says a diffuse alternative predicted the data worse than the point null did.

**Sensitivity is the honest problem.** $\mathrm{BF}_{10}$ depends on $\tau$ in a way posteriors do not. Maximize (3.40) over $u$: $\frac d{du}\log\mathrm{BF}_{10}=-\frac1{2(1+u)}+\frac{z^2}{2(1+u)^2}=0$ at $1+u=z^2$, so
$$\max_\tau\mathrm{BF}_{10}=\frac1{|z|}\,e^{(z^2-1)/2}\quad(|z|>1). \tag{3.41}$$
For $z=2$ that is $2.24$, reached at $\tau\approx0.017$. No Normal prior does better, and no prior of any shape does better than $e^{z^2/2}=7.39$ (put a point mass at $\bar x$: $\mathrm{BF}_{10}\le p(\bar x\mid\mu=\bar x)/p(\bar x\mid0)$). A $p$-value of $0.0455$ is, at best, between 2:1 and 7:1 evidence against the null, which is far weaker than the number suggests (Edwards, Lindman and Savage, 1963, and Berger and Sellke, 1987, made this point). With *improper* priors the constant is arbitrary and $\mathrm{BF}$ is meaningless; use proper priors and report $\mathrm{BF}$ across a range of scales.

**BIC is the large-$n$ evidence.** Take $d$ parameters and apply the Laplace expansion (3.35) to $\int e^{\ell_n(\theta)}p(\theta)d\theta$ with Hessian $-n\hat J$, $\hat J\to I_1$:
$$\log p(y\mid M)=\ell_n(\hat\theta)+\log p(\hat\theta)+\tfrac d2\log2\pi-\tfrac d2\log n-\tfrac12\log|\hat J|+o(1).$$
Drop the terms that do not grow with $n$ (the prior density at $\hat\theta$, $2\pi$, $|\hat J|$):
$$-2\log p(y\mid M)\ \approx\ \mathrm{BIC}=-2\ell_n(\hat\theta)+d\log n. \tag{3.42}$$
The error in $\log\mathrm{BF}$ is $O(1)$, not vanishing: those dropped terms are where the prior lives, so BIC is "prior-free" only by ignoring a constant that can matter. Kass and Wasserman showed it is the evidence under a unit-information prior up to a smaller error. **Information criteria** that target *predictive* accuracy instead (WAIC; leave-one-out cross-validation via Pareto-smoothed importance sampling) are computed from posterior draws in ST-6; they answer "which model will predict new data better after fitting", while the evidence answers "which model predicted these data better *before* fitting". **Bayesian model averaging** keeps all models: $p(\tilde y\mid y)=\sum_kP(M_k\mid y)\,p(\tilde y\mid y,M_k)$, with weights proportional to evidence times prior probability; it inherits the evidence's prior sensitivity.

### 3.13 Two answers to two questions

Frequentist and Bayesian inference are valid answers to different questions. The frequentist asks how a *procedure* behaves over repeated data at a fixed $\theta$; the Bayesian asks what to believe about $\theta$ given the *data in hand*, from a stated prior. They **agree** for large $n$ (Theorem 3.11) and for the flat-prior location–scale model of §3.4. They **disagree** for small $n$, for strong priors (the Beta(10,10) coverage above), for testing a point null (the paradox above), and on **stopping rules**. The likelihood principle says all the evidence about $\theta$ in $y$ is in the likelihood function $p(y\mid\theta)$ up to a constant. The posterior obeys it; $p$-values do not. Example 4.12: 9 heads and 3 tails give the same likelihood $\propto p^9(1-p)^3$ whether you planned 12 flips or planned to stop at the third tail, so the posterior is identical, but the one-sided $p$-values against $p=\frac12$ are $0.0730$ and $0.0327$: the frequentist answer depends on a plan that was not otherwise recorded in the data.

The position this course takes: use the method that answers the question you are asking. Where they coincide the choice is free; where they differ, compute both and ask which assumption (the prior, or the sampling plan) is doing the work. Defend every prior with a prior predictive check, check calibration by simulation (the project does), prefer predictive criteria (ST-6) over Bayes factors for model selection, and report Bayes factors only with proper priors across a range of scales.

## 4. Worked examples

Every number below is asserted or printed by `project.py` (the part is named in parentheses).

### Example 4.1: Beta–Binomial, Beta(2,2) prior, 7 successes in 10 (Part 1, 2)

By (3.4) the posterior is Beta$(2+7,\,2+3)=$ Beta(9,5). Mean $9/14=0.6429$; check against (3.6): $w=10/14=0.7143$ and $0.7143(0.7)+0.2857(0.5)=0.6429$. By (3.7), $\mathrm{Var}=\frac{9\cdot5}{14^2\cdot15}=0.01531$, sd $0.1237$. The 95% equal-tailed interval is $(0.3857,\,0.8614)$; the HPD interval is $(0.4013,\,0.8737)$ (equal density, $p(0.4013\mid y)=p(0.8737\mid y)$, slightly shorter and shifted toward the mode). The probability the next trial succeeds is the mean, $0.6429$; the uniform prior would give Laplace's $8/12=0.6667$. The evidence (3.5) is $p(y{=}7)=0.1119$, against $1/11=0.0909$ for a uniform prior: the Beta(2,2) prior put more mass near $\frac12$ and predicted this outcome slightly better. For ten future trials the predictive (3.8) has mean $6.43$ and variance $3.673$ (Monte Carlo from 400,000 draws of $p$ then $\tilde y$: $3.674$), against $2.10$ for Binomial(10, 0.7) and $2.30$ for Binomial(10, 0.6429): the plug-in forecasts are too confident by the factor $(14+10)/(14+1)=1.6$ relative to the second.

### Example 4.2: Normal mean, precisions added (Part 1)

Eight measurements $0.5,3.1,1.2,2.4,-0.3,2.9,0.8,1.4$ ($\bar y=1.5$), known $\sigma^2=4$, prior $\mathcal N(0,1)$. Precisions: $1/1+8/4=3$, so $\tau_n^2=\frac13$. By (3.10), $\mu_n=\frac13(0+8\cdot1.5/4)=1.0$. The data carry weight $2/3$ (their precision 2 out of 3), so the mean sits at $\frac23(1.5)+\frac13(0)=1.0$. Posterior $\mathcal N(1,\frac13)$, 95% interval $(-0.1316,\ 2.1316)$ against the frequentist $1.5\pm1.386$. A new observation: $\mathcal N(1,\ 4+\frac13)$ by (3.11), variance $4.333$ against the plug-in $4$; the extra $\frac13$ is the price of not knowing $\mu$.

### Example 4.3: Ten measurements of $g$, Normal–Inverse-Gamma, against ST-2 (Part 5)

ST-2's data: $\bar y=9.806$, $S=0.01564$, $n=10$. *Reference prior:* $\mu\mid y\sim t_9(9.806,\ s^2/10)$, interval $(9.7762,\ 9.8358)$, **identical to ST-2's $t$ interval to four decimals**; $\sigma^2\in(0.00082,\ 0.00579)$, $\sigma\in(0.0287,\ 0.0761)$, ST-2's $\chi^2$ interval. *Proper prior* $\mathrm{NIG}(\mu_0{=}9.5,\kappa_0{=}1,\alpha_0{=}2,\beta_0{=}0.02)$ (prior mean $\sigma^2=0.02$): (3.13) gives $\kappa_n=11$, $\mu_n=(9.5+98.06)/11=9.7782$, $\alpha_n=7$, $\beta_n=0.02+\frac12[0.01564+\frac{10}{11}(0.306)^2]=0.07038$. The marginal (3.14) is $t_{14}(9.7782,\ 0.07038/77)$, interval $(9.7133,\ 9.8430)$: the prior, 0.3 below the data, drags the mean down by $0.028$ and, through the conflict term in $\beta_n$ (which supplies $0.0426$ of the $0.0704$), inflates the variance, so the interval is *wider* than ST-2's ($\pm0.065$ against $\pm0.030$). The predictive (3.15) for the next measurement: $(9.5536,\ 10.0028)$. The closed-form marginal agrees with brute-force integration over a $(\mu,\sigma^2)$ grid to $0.3\%$ of the peak, and the joint-Jeffreys $t_{10}$ half-width is $0.0279$ against the reference $0.0298$.

### Example 4.4: Exponential-family update reproduces Poisson–Gamma (Part 6)

Counts $3,5,2,4,6$, prior Gamma(2,1). In (3.18) coordinates the Poisson prior is $(\chi,\nu)=(\alpha_0,\beta_0)=(2,1)$. By (3.19), $\chi'=2+20=22$ and $\nu'=1+5=6$: Gamma(22,6), exactly what (3.16) gives. The posterior mean $22/6=3.667$ is (3.20): $\frac{1}{6}(2)+\frac56(4)=3.667$. A grid posterior matches the Gamma(22,6) density to $10^{-4}$, and quadrature confirms $\mathbb E[e^\eta]=\chi/\nu$ under the prior $e^{\eta\chi-\nu e^\eta}$. The predictive (3.17) for the next count has mean $3.667$ and variance $\frac{22\cdot7}{36}=4.278$, against the Poisson plug-in's $3.667$; $P(\tilde y=0)=0.0337$ (Monte Carlo agrees to $0.002$ in every cell).

### Example 4.5: Jeffreys' prior for Bernoulli and Poisson (Part 4)

*Bernoulli.* $\log p=y\log p+(1-y)\log(1-p)$; $-\partial_p^2\log p=\frac y{p^2}+\frac{1-y}{(1-p)^2}$, expectation $\frac1p+\frac1{1-p}=\frac1{p(1-p)}$. So $p_J\propto p^{-1/2}(1-p)^{-1/2}$. *Poisson.* $\log p=y\log\lambda-\lambda-\log y!$; $-\partial_\lambda^2\log p=y/\lambda^2$, expectation $1/\lambda$, so $p_J\propto\lambda^{-1/2}$. (SymPy confirms both and $I_\mu=\sigma^{-2}$, $I_\sigma=2/\sigma^2$.) *Posterior after 0 successes in 10 trials:* Jeffreys gives Beta($\frac12,\frac{21}2$), mean $0.0455$, 95% equal-tailed $(0.00005,\ 0.2172)$; the uniform prior gives Beta(1,11), mean $0.0833$, HPD $(0,\ 0.2384)$ (the pile-up at 0 makes the HPD one-sided). *Poisson, zero counts in 10 intervals:* Gamma($\frac12,10$), mean $0.05$. *Invariance:* pushing the flat prior to the logit scale gives the logistic density (3.21), which the project verifies equals `stats.logistic.pdf`; pushing Beta($\frac12,\frac12$) gives $1/(2\pi\cosh(\eta/2))$, matching $\sqrt{I_\eta}/\pi$ with $I_\eta=p(1-p)$ as Theorem 3.8 requires.

### Example 4.6: MAP is ridge, with $\lambda=\sigma^2/\tau^2$ (Part 7)

ST-3's five points centered: $x_c=(-2,-1,0,1,2)$, $y_c=(-2,0,1,0,1)$, so $S_{xx}=10$, $S_{xy}=6$, OLS slope $0.6$. Take $\sigma^2=0.8$ and prior slope $\mathcal N(0,\tau^2=0.5)$: $\lambda=0.8/0.5=1.6$. The MAP slope is $6/(10+1.6)=0.5172$, the OLS slope times the shrinkage factor $10/11.6=0.8621$; a numerical minimization of $-\log$ posterior agrees. The posterior is Gaussian, so mean $=$ mode $=0.5172$, with variance $0.8/11.6=0.06897$, smaller than the OLS variance $0.8/10=0.08$: the prior adds information. (The lasso version, soft-thresholding, is in §3.9: with $z=1.2$ and Laplace scale 2, MAP $=0.7$ and posterior mean $0.8633$.)

### Example 4.7: Bayesian regression on five points (Part 8)

ST-3 Example 4.1: $x=1,\dots,5$, $y=(2,4,5,4,5)$, design $[\mathbf 1,x]$, treat $\sigma^2=0.8$ as known, prior $\beta\sim\mathcal N(0,10I)$ ($\lambda=0.08$). From (3.30): posterior mean $(2.0363,\ 0.6437)$ (OLS: $(2.2,\ 0.6)$), sd $(0.897,\ 0.272)$. The posterior mean equals the ridge solution (3.31), the minimizer of $-\log$ posterior, and the F5 conditional formula, all asserted. Predictive (3.32) at the interpolation point $x_*=3$: mean $3.967$, variance $0.958$ (sd $0.98$). At the extrapolation point $x_*=8$: mean $7.186$, variance $2.846$ (sd $1.69$), three times larger. In the flat-prior limit the variances are $0.8(1+h_*)$ with $h_*=\frac15+\frac{(x_*-3)^2}{10}=0.2,\ 2.7$, i.e. $0.96$ and $2.96$, ST-3's prediction variance. The 50 posterior lines fan out away from $x=3$ and the $\pm2$ sd band flares. With unknown $\sigma^2$ the predictive is a $t$ (3.33), and the project checks that the closed-form evidence (3.39) equals the product of the five one-step $t$ predictives, which tests both formulas at once.

### Example 4.8: Laplace approximation, logistic regression (Part 9)

*One parameter.* Data $x=(-2,-1,0.5,1,2,3)$, $y=(0,0,1,0,1,1)$, no intercept, prior $\beta\sim\mathcal N(0,2^2)$. Newton gives $\hat\beta=1.0508$ and $A=1/0.7407^2$, so the Laplace posterior is $\mathcal N(1.051,\ 0.741^2)$. A fine-grid posterior has mean $1.379$ and sd $0.843$: the true posterior is skewed right, and the Laplace approximation, a local fit at the peak, misses that (mode $\ne$ mean, sd $12\%$ low). The evidence (3.35) is much better: $\log p(y)=-3.519$ against $-3.4725$ exactly, a $5\%$ error in $p(y)$. *Two parameters,* $n=60$ (Part 9): MAP $(0.037,\ 1.670)$, Laplace sd $(0.307,\ 0.464)$, against grid mean $(0.046,\ 1.780)$, sd $(0.314,\ 0.481)$; $\log$ evidence $-35.847$ against $-35.822$. At $n=60$ it is good enough for most purposes; the six-point example shows it visibly off in the mean and spread.

### Example 4.9: Evidence picks the degree (Part 10)

Thirty points from $y=0.5+x-2x^2+\mathcal N(0,0.3^2)$ on $[-1,1]$; polynomial models of degree 0 to 7 with prior $\beta\sim\mathcal N(0,2I)$ and known $\sigma=0.3$. Closed form (3.38), cross-checked against the multivariate Normal density:

| degree | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| max log-lik (MLE) | $-126.80$ | $-75.81$ | $-7.26$ | $-6.82$ | $-4.97$ | $-4.87$ | $-4.63$ | $-4.48$ |
| fit at posterior mean | $-126.80$ | $-75.81$ | $-7.28$ | $-6.85$ | $-5.47$ | $-5.37$ | $-5.34$ | $-5.30$ |
| log Occam factor | $-3.26$ | $-6.21$ | $-9.36$ | $-10.95$ | $-12.64$ | $-13.27$ | $-13.66$ | $-14.14$ |
| log evidence (fit + Occam) | $-130.06$ | $-82.02$ | $\mathbf{-16.64}$ | $-17.80$ | $-18.12$ | $-18.64$ | $-19.00$ | $-19.44$ |
| $-\mathrm{BIC}/2$ | $-128.50$ | $-79.21$ | $\mathbf{-12.36}$ | $-13.62$ | $-13.48$ | $-15.07$ | $-16.54$ | $-18.09$ |

The maximized likelihood never decreases with degree (the models are nested) and the Occam factor only worsens, so the evidence peaks where extra flexibility stops paying: the true degree 2. From degree 2 to 3 the fit gains $0.43$ and the Occam factor loses $1.59$; the net is $\mathrm{BF}_{23}=e^{1.16}=3.2$. BIC (3.42) peaks at the same degree but penalizes more heavily here than this fairly diffuse prior does.

### Example 4.10: The Jeffreys–Lindley paradox (Part 11)

$n=10^4$, $\bar x=0.02$, $\sigma=1$: $z=0.02\sqrt{10^4}=2$, two-sided $p=0.0455$. With $\tau=1$, $u=10^4$, (3.40) gives $\mathrm{BF}_{10}=(10001)^{-1/2}\exp\{2\cdot\frac{10^4}{10001}\}=0.00999\times7.388=0.0739$ (quadrature of $\int\mathcal N(\bar x;\mu,\sigma^2/n)\,\mathcal N(\mu;0,1)\,d\mu$ agrees), so $\mathrm{BF}_{01}=13.5$. The most favorable $\tau$ for $H_1$ is $0.0173$, giving $\mathrm{BF}_{10}=2.241=\frac12e^{3/2}$ as (3.41) predicts, and the bound for any prior is $e^{2}=7.389$. Under $H_1$ the posterior of $\mu$ is $\mathcal N(0.0200,\ 0.0100^2)$ with 95% interval $(0.0004,\ 0.0396)$. The project plots $p$ and $\mathrm{BF}_{01}$ as $n$ grows with $z=2$ fixed: $\mathrm{BF}_{01}<1$ at $n=10$ and above $10^3$ by $n=10^8$.

### Example 4.11: Posterior predictive check (Part 12)

Sixty draws from a lognormal$(0,0.8^2)$ have sample skewness $2.712$. Fit a Normal (reference prior), draw 4,000 $(\mu,\sigma^2)$ from the posterior, simulate 4,000 replicated datasets of size 60, compute each skewness: none reaches $2.712$ ($p_B<1/4000$). The Normal model cannot produce data this skewed. Sixty actual Normal draws give $p_B=0.665$. The cure here is a different likelihood (lognormal, or Normal on $\log y$), not a different prior.

### Example 4.12: Same data, two stopping rules (Part 13)

Nine heads and three tails. *Plan A,* 12 flips: $P(Y\ge9\mid p=\frac12)=0.0730$. *Plan B,* flip until the third tail: $P(\ge9\text{ heads before the 3rd tail})=0.0327$. Both likelihoods are $\propto p^9(1-p)^3$, so the uniform-prior posterior is Beta(10,4) in both (the project checks the grid posteriors coincide). A 5% test rejects under B and not under A.

## 5. Why it matters and where it connects

- **ST-6 (Bayesian computation).** Everything here stops at the model whose integral (3.2) has a closed form. For the rest, MCMC draws $\theta^{(s)}\sim p(\theta\mid y)$ and the predictive (3.3) and checks (3.29) are computed from the draws. WAIC and leave-one-out (§3.12) are computed there; hierarchical models (the previewed $\tau^2$ with a prior) live there.
- **ML-4 (Gaussian processes).** §3.11's predictive (3.32) with $\phi(x)$ infinite-dimensional is GP regression; the kernel is $k(x,x')=\phi(x)^\top\Sigma_0\phi(x')$. ST-3's kernel ridge regression is the posterior *mean*; the variance flaring away from data is what the Bayesian view adds. GP classification uses the Laplace approximation (3.34)–(3.35).
- **ML-1, ML-2, ML-3.** An $\ell^2$ penalty is a Gaussian prior (3.27) and weight decay is its gradient; $\ell^1$ is a Laplace prior (3.28); cross-entropy loss is ST-4's negative log-likelihood.
- **ML-10 (VAEs, the ELBO).** The evidence (3.36) is the quantity a latent-variable model maximizes; when the integral is intractable one maximizes a lower bound on $\log p(y)$ and the gap is a KL divergence (F5 Definition 3.44) between an approximate and the true posterior.
- **SM-4 (canonical ensemble).** The conjugate normalizer $Z(\chi,\nu)$ in (3.18) is a partition function; the evidence is a ratio of partition functions, $Z(\chi',\nu')/Z(\chi,\nu)$, and $\ln Z$ differences are free energies. The Occam factor (3.37) is a ratio of phase-space volumes: the volume of parameter space consistent with the data over the volume the prior allowed.
- **ST-8 (state-space models).** The Kalman filter is (3.10) and (3.30) applied one observation at a time; sequential equals batch (§3.1) is why it works.
- **ST-9 (causal inference).** A prior does not remove confounding. Bayesian regression of $y$ on $x$ is exactly as causal as OLS.

## 6. Exercises

See [exercises.md](exercises.md): 16 problems. A1–A5 are drills (each family's update, a credible interval in closed form). B1–B8 are the module: predictive distributions, the exponential-family update for a new family, reference priors, decision theory, an A/B test (code), and the Bayes-factor bound. C1–C3 are the challenges: a minimax estimator from a constant-risk Bayes rule, the Jeffreys–Lindley crossover, and Laplace-versus-BIC on real code. Essential: A1, A2, B1, B2, B3, B5, B8. B7 and C3 modify `project.py`.

## 7. Coding project

`project.py` (about 15 seconds) has 13 parts, an exercise-check cell, and "Your turn". (1) Beta–Binomial and Normal–Normal updating as data arrive (figure grid), sequential asserted equal to batch. (2) Posterior predictive against plug-in, with Monte Carlo ($\theta$ then $\tilde y$) checked against the analytic predictives. (3) Calibration (Theorem 3.10) and exact frequentist coverage against $n$ and prior strength. (4) SymPy Fisher information, Jeffreys priors, induced logit-scale priors, the non-invariant mode, a prior predictive check. (5) NIG: the $t$ marginal against ST-2's interval and a brute-force grid. (6) Exponential-family conjugacy reproducing Beta, Gamma and Normal updates, and the Negative Binomial predictive. (7) Mean, median, mode by minimizing each loss; MAP equals ridge; lasso as soft-thresholding. (8) Bayesian regression: contours, 50 sampled lines, the predictive band, ridge and F5-conditional equivalences, and the chain-rule check of the NIG evidence. (9) Laplace against a 2-D grid posterior. (10) Polynomial evidence, Occam factors, BIC. (11) Jeffreys–Lindley. (12) Posterior predictive checks. (13) Stopping rules.

**Build it yourself first:** write `beta_update`, `normal_update` and `blr_posterior` (about ten lines) from (3.4), (3.10) and (3.30) before reading Parts 1 and 8. Predict which of the 50 posterior lines in Part 8 will fan out, and where. Exercise B7 adds a third arm to the A/B test (reference solution at the bottom of the script); C3 extends Part 9 with BIC.

## 8. Where you'll struggle

- **"The prior is subjective" is true and the wrong objection.** The likelihood is a modeling choice too (Normal errors, independence, a linear mean), and nobody calls that unscientific. What is actually at stake is narrower: with small $n$ the prior moves the answer (the Beta(10,10) coverage of 0.000 at $n=10$), and Bayes factors move with the prior even when posteriors do not (Example 4.10). The remedy is the same as for any assumption: state it, check it with a prior predictive check, and show a sensitivity analysis.
- **Credible and confidence intervals are different answers to different questions,** and the student who understood ST-2 must now hold both. §3.4 makes this harder because the numbers coincide: $(9.7762,\ 9.8358)$ is both. Test yourself: for each reading, name what is random and what is fixed. In ST-2 the interval is random and $g$ is fixed; here $g$ is random and the interval is fixed. The coincidence is a property of one prior (§3.4), not a law (Beta(10,10) breaks it).
- **You will complete the square for the fourth time.** §3.3, §3.4 and §3.11 are the same computation as F5 Theorem 3.24. If (3.30) still feels like a trick, stop and redo F5's block-LDU derivation; ST-8 and ML-4 will need it fluently. Work in *information form* (precision and precision×mean), where updates are additions.
- **Flat priors are not noninformative,** and Jeffreys' prior is a choice, not a cure: it depends on the sampling model, so it violates the likelihood principle (for the Binomial it is $p^{-1/2}(1-p)^{-1/2}$; for the Negative Binomial, whose Fisher information is $r/(p^2(1-p))$, it is $p^{-1}(1-p)^{-1/2}$, though the likelihood is the same kernel).
- **Bayes factors are the least settled part of the field.** They are sensitive to prior scale in a way no posterior is (Example 4.10), undefined for improper priors, and the prior scale can be tuned to produce almost any value within the bound (3.41). Model comparison by predictive accuracy (ST-6) is more robust; use the evidence only with defended proper priors.
- **MAP is not Bayesian in spirit.** It uses the prior but ignores the posterior's shape, is not invariant to reparametrization (Example 4.1's mode, $2/3$, against $9/14$ on the logit scale), and yields exact lasso zeros that the posterior does not have.
- **Parametrization bugs cost hours.** Gamma and Inverse-Gamma come in shape–rate and shape–scale forms; `scipy.stats.gamma(a, scale=1/rate)` and `scipy.stats.invgamma(a, scale=b)`. The predictive variance includes $\sigma^2$ as well as the parameter uncertainty (a common slip is plotting the band for $f(x)$ and calling it a prediction band).

## 9. Check yourself

1. State the posterior, the evidence, and the posterior predictive as integrals of the prior and likelihood.
2. Beta(α₀,β₀) prior, $y$ successes in $n$ trials: give the posterior, its mean, and the weight on the data. What is $\alpha_0+\beta_0$ in words? Why is the predictive for $m$ future trials more spread out than Binomial$(m,\hat p)$?
3. State the Normal–Normal update. What adds? What is $n_0$?
4. What is the posterior of $\mu$ under $p(\mu,\sigma^2)\propto\sigma^{-2}$, how does its 95% interval compare with ST-2's, and what differs in interpretation?
5. State the exponential-family conjugate update. What do $\chi$ and $\nu$ mean? Verify it for Poisson–Gamma.
6. Derive the logit-scale density of a uniform prior on $p$. What does it say about "flat"?
7. Define Jeffreys' prior, prove it is invariant under reparametrization, and give it for Bernoulli.
8. What is the guarantee in Theorem 3.10, and why is it not a statement about a fixed $\theta$?
9. Prove the posterior mean minimizes expected squared error. Which prior makes MAP equal ridge, and what is $\lambda$?
10. Give the Bayesian linear regression posterior and predictive variance. Why does the band widen away from the data?
11. What is the Occam factor, and why does the evidence penalize complexity without a penalty term?
12. State the Jeffreys–Lindley paradox and explain why a significant $p$-value and a Bayes factor favoring $H_0$ can both be right.

## 10. Pair with

- **Gelman, Carlin, Stern, Dunson, Vehtari and Rubin, *Bayesian Data Analysis*, 3rd ed. (free PDF), ch. 1–5 and 14.** The reference for conjugate derivations, the general NIG algebra, and regression; ch. 6–7 for posterior predictive checks and model comparison. Denser than these notes; better for worked applied examples.
- **McElreath, *Statistical Rethinking*, 2nd ed., ch. 1–4 and 7 (and his lecture videos).** The best way to build the "grid approximation first" intuition and to learn prior predictive simulation; ch. 7 on overfitting and information criteria prepares ST-6.
- **Hoff, *A First Course in Bayesian Statistical Methods*, ch. 3, 5, 9.** The clearest conjugate derivations in print; read it for the Normal–Gamma algebra if §3.4 felt rushed.
- **Murphy, *Probabilistic Machine Learning: An Introduction*, ch. 4.** The machine-learning view of the same models, including Bayesian regression and the Laplace approximation, in the notation ML-4 and ML-10 use.
- **Jaynes, *Probability Theory: The Logic of Science*, ch. 1–4.** The philosophical case for probability as extended logic; polemical and worth reading skeptically.
- **MacKay, *Information Theory, Inference, and Learning Algorithms* (free), ch. 28.** The Occam factor argument of §3.12 in its original form, with pictures.
