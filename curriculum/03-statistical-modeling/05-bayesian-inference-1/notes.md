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
