# ST-2: Intervals and Tests: Sampling Distributions, Pivots, the Bootstrap, Neyman–Pearson, Likelihood-Ratio Tests, p-Values and Their Abuse, Multiple Testing

**Time:** ~12 hours · **Prerequisites:** ST-1 · **Unlocks:** ST-3

> **Conventions in this module.** Data $X_1,\dots,X_n$ iid, stacked as $\mathbf X\in\mathbb R^n$; $\mathbf 1 = (1,\dots,1)^\top$. $\bar X = \frac1n\sum X_i$; $S^2 = \frac{1}{n-1}\sum(X_i-\bar X)^2$ (divisor $n-1$); $\hat\sigma^2 = \frac1n\sum(X_i-\bar X)^2$ (divisor $n$, the MLE). Confidence level $1-\alpha$: a $(1-\alpha)$ confidence set $C(\mathbf X)$ has $P_\theta(\theta\in C(\mathbf X))\ge 1-\alpha$ for every θ. Null $H_0$, alternative $H_1$, test statistic $T$, rejection region $R$, size $\alpha = \sup_{\theta\in H_0}P_\theta(\text{reject})$, **power function** $\beta(\theta) = P_\theta(\text{reject})$. Because β is taken, the type II error rate at a design alternative $\theta_1$ is written $\beta_{\rm II} = 1-\beta(\theta_1)$ (applied texts call it β). **Quantiles are upper-tail:** $P(Z>z_a) = a$ for $Z\sim\mathcal N(0,1)$, so $z_{0.025} = 1.960$ and $z_{0.05} = 1.645$; likewise $P(T>t_{k,a}) = a$ for $T\sim t_k$ and $P(V>\chi^2_{k,a}) = a$ for $V\sim\chi^2_k$, so $\chi^2_{k,0.975}$ is the *small* quantile. $\Phi,\varphi$ are the standard Normal CDF and density; $\log$ is natural. $\ell(\theta)$ is the log-likelihood and $I_1(\theta)$ the per-observation Fisher information (ST-1). $\hat F_n$ is the empirical CDF; bootstrap replicates are $\hat\theta^*_1,\dots,\hat\theta^*_B$.

## 0. Where this sits

ST-1 gave you estimators, their sampling distributions, Fisher information, and the asymptotic normality of the MLE. This module turns a sampling distribution into the two things people report: an interval and a test. You need ST-1's MLE and $I_1$, LA-5's projections, LA-8's spectral theorem, and F4/F5's Normal, $\chi^2$ and CLT. ST-3 (regression) uses all of it: $t$ tests on coefficients, $F$ tests for nested models, and §3.1's projection proof in $p$ dimensions. ST-5 returns to testing from the Bayesian side.

## 1. The problem that forced this idea

A point estimate alone is nearly useless. "The new alloy's mean fatigue life is 1.3 million cycles" means nothing until you know whether the number is good to ±1% or ±50%. And "does the treatment work?" needs a decision rule whose error rate you know before you see the data.

William Sealy Gosset hit the first problem at the Guinness brewery in Dublin. Brewing experiments gave him four or five measurements, and the method of the day (estimate σ by $s$, then use Normal tables) gave intervals that were too narrow. In "The probable error of a mean" (*Biometrika*, 1908), published as "Student" because Guinness restricted staff publication, he found the distribution of $(\bar X-\mu)/s$ for Normal samples. He could not prove every step, so he checked the result empirically, with 750 samples of size 4 drawn from cards recording measurements of 3,000 criminals. A complete proof came from the young R. A. Fisher, who in 1912 sent Gosset an argument using the geometry of $n$-dimensional space. It is the projection argument of §3.1.

Fisher's *Statistical Methods for Research Workers* (1925) made significance testing a laboratory tool. He tabulated the distributions, called $P = 0.05$ a convenient point for judging significance, and treated the p-value as graded evidence against a null hypothesis. He did not mean 0.05 as a rule: in 1956 he wrote that no scientific worker has a fixed significance level at which he rejects hypotheses in all circumstances.

Jerzy Neyman and Egon Pearson asked a different question: among procedures with a fixed false-alarm rate, which detects real effects most often? Their 1933 paper (*Phil. Trans. R. Soc. A*) introduced the alternative hypothesis, the two error types, power, and the lemma of §3.5. Neyman's 1937 paper defined the confidence interval by its long-run coverage. Fisher attacked the framework for the rest of his life. Textbooks then fused the two: Fisher's p-value compared with Neyman's fixed α, with power computed only when a grant application demands it. Neither man endorsed the hybrid.

Two later events matter. Bradley Efron's "Bootstrap methods: another look at the jackknife" (*Annals of Statistics*, 1979) replaced analytic sampling distributions by resampling on a computer, so an interval can be attached to a median, a correlation, or a 500-line pipeline. The replication crisis, from John Ioannidis's "Why most published research findings are false" (*PLoS Medicine*, 2005) to the Open Science Collaboration's 2015 attempt to replicate 100 psychology studies (97% of the originals were significant, 36% of the replications), showed what the hybrid does when combined with low power, flexible analysis and a publication filter. In 2016 the American Statistical Association published a statement on p-values (Wasserstein and Lazar), the first time it had taken a position on a specific matter of statistical practice.

The forcing problem: every estimate needs an honest uncertainty, and every decision between hypotheses needs a calibrated procedure. "Calibrated" means "over hypothetical repetitions of the experiment," and almost everyone misreads that.

## 2. Intuition

### An interval is a ring toss, and the stake does not move

Picture a stake in the ground (the true μ) and a machine that throws a ring each time you run the experiment. The ring lands somewhere different each time, because the data differ. A 95% procedure is a machine whose ring lands around the stake on 95% of throws. After a throw, the ring is where it is: around the stake or not. The randomness is in the ring, never in the stake.

### Pivots: how to build the machine

If $X_1,\dots,X_n\sim\mathcal N(\mu,1)$, then $\bar X-\mu\sim\mathcal N(0,1/n)$ whatever μ is. A quantity that involves the unknown parameter but whose distribution does not depend on it is a **pivot**. Once you have one, $P(|\bar X-\mu|\le 1.96/\sqrt n) = 0.95$ for every μ, and you rearrange the inequality to trap μ. If σ is unknown, divide by an estimate of it; the estimate's own noise fattens the tails. Guess: at $n = 5$ the 1.96 must grow. (It grows to 2.776.)

### A test asks how surprising the data would be if nothing were going on

Fix a null hypothesis. Compute a statistic. Ask: if $H_0$ were true, how often would a statistic at least this extreme appear? That frequency is the **p-value**. It is computed *assuming* $H_0$, so it cannot be the probability that $H_0$ is true. Guess its null distribution: "the fraction of null datasets more extreme than this one" is a rank, and a random rank is uniform. Under $H_0$ the p-value is equally likely to be anywhere in $[0,1]$.

### Choosing the rejection region is a budget problem

You have a false-alarm budget α. Putting a possible data point $x$ in the rejection region costs $f_0(x)$ (its probability under $H_0$) and earns $f_1(x)$ (its probability under $H_1$) in power. To earn the most on a fixed budget, buy the best value-per-cost items first: put $x$ in the region when $f_1(x)/f_0(x)$ is large. That is the Neyman–Pearson lemma, and the ratio is the likelihood ratio.

### The bootstrap: the sample stands in for the population

The sampling distribution of a median depends on the unknown population. The sample is your best estimate of that population. So draw new samples of size $n$ *from the sample*, with replacement, recompute the median each time, and watch it wobble. That wobble estimates the real one. Guess where it breaks: a resample can never exceed the observed maximum, so for a statistic like the maximum, resampling cannot see what lies beyond the data.

### Many tests, many false alarms

Run 20 independent tests at α = 0.05 with every null true. The chance of no false alarm is $0.95^{20} = 0.36$, so the chance of at least one is 0.64. Guess the fix: test each at $0.05/20$, and the union bound caps the chance of any false alarm at 0.05.

## 3. Formal development

### 3.1 Sampling distributions for Normal data

In this subsection $X_1,\dots,X_n$ are iid $\mathcal N(\mu,\sigma^2)$. Write $\mathbf X = \mu\mathbf 1+\sigma\mathbf Z$ with $\mathbf Z\sim\mathcal N(\mathbf 0,I_n)$.

**Theorem 3.1 (Sampling distributions of $\bar X$ and $S^2$).** (a) and (b):
$$\bar X\sim\mathcal N(\mu,\sigma^2/n),\qquad (n-1)S^2/\sigma^2\sim\chi^2_{n-1}. \tag{3.1}$$
(c) $\bar X$ and $S^2$ are independent.

▶ **Try it first.** $\bar X\mathbf 1$ is the projection of $\mathbf X$ onto the line through $\mathbf 1$ (LA-1 Example 4.5's formula with $\mathbf v = \mathbf 1$). The residual $\mathbf X-\bar X\mathbf 1$ lives in the $(n-1)$-dimensional orthogonal complement. Rotate to a basis adapted to that split and ask what a rotation does to $\mathcal N(\mathbf 0,I)$.

**Derivation.**
1. Choose an orthonormal basis $\mathbf q_1 = \mathbf 1/\sqrt n,\mathbf q_2,\dots,\mathbf q_n$ of $\mathbb R^n$ (Gram–Schmidt, LA-5 Theorem 3.19, applied to $\mathbf 1,\mathbf e_2,\dots,\mathbf e_n$). Let $Q = (\mathbf q_1\cdots\mathbf q_n)$, so $Q^\top Q = I$.
2. Let $\mathbf W = Q^\top\mathbf Z$. By F5 Theorem 3.22, $\mathbf W\sim\mathcal N(\mathbf 0,Q^\top Q) = \mathcal N(\mathbf 0,I)$, and by F5 Theorem 3.23 its coordinates are independent. A rotation of white noise is white noise.
3. The mean: $\mathbf q_1^\top\mathbf X = \sqrt n\,\bar X$, and also $\mathbf q_1^\top\mathbf X = \mu\sqrt n+\sigma W_1$. So $\bar X = \mu+\sigma W_1/\sqrt n$, which is (a).
4. The residual: $\mathbf X-\bar X\mathbf 1 = (I-\mathbf q_1\mathbf q_1^\top)\mathbf X = \sigma(I-\mathbf q_1\mathbf q_1^\top)\mathbf Z$, because $(I-\mathbf q_1\mathbf q_1^\top)\mathbf 1 = \mathbf 0$. Expanding $\mathbf Z = \sum_iW_i\mathbf q_i$, this is $\sigma\sum_{i\ge2}W_i\mathbf q_i$.
5. Its squared length, by Pythagoras for orthonormal vectors (LA-5 Theorem 3.2): $(n-1)S^2 = \|\mathbf X-\bar X\mathbf 1\|^2 = \sigma^2\sum_{i=2}^nW_i^2$. A sum of $n-1$ squared iid standard Normals is $\chi^2_{n-1}$ (F5 §3.5). That is (b).
6. $\bar X$ is a function of $W_1$ alone and $S^2$ of $W_2,\dots,W_n$ alone, so they are independent. That is (c). ∎

■ **What just happened.** The data vector splits into a component along $\mathbf 1$ (one dimension: the mean) and a component perpendicular to it ($n-1$ dimensions: the spread). Spherical Gaussian noise does not care which axes you use, so the components are independent. **Degrees of freedom are the dimension of the subspace the residual lives in.** You lose one because you projected out one direction; ST-3 projects out $p$ and gets $n-p$. Part (c) is special to the Normal: for Exponential data the project finds $\operatorname{corr}(\bar X,S^2)\approx0.69$.

**Definition 3.2 ($t$ and $F$).** If $Z\sim\mathcal N(0,1)$ and $V\sim\chi^2_k$ are independent,
$$T = \frac{Z}{\sqrt{V/k}}\sim t_k. \tag{3.2}$$
If $V_1\sim\chi^2_{k_1}$ and $V_2\sim\chi^2_{k_2}$ are independent, $F = (V_1/k_1)/(V_2/k_2)\sim F_{k_1,k_2}$. The $t_k$ density is proportional to $(1+t^2/k)^{-(k+1)/2}$ (exercise C1 derives it with F5 Theorem 3.6's change of variables). Its tails are polynomial, and $t_k\to\mathcal N(0,1)$ as $k\to\infty$ because $V/k\to1$ by the law of large numbers. From the definitions, $T\sim t_k$ implies $T^2\sim F_{1,k}$.

**Corollary 3.3 (The $t$ statistic).** For Normal data,
$$T = \frac{\bar X-\mu}{S/\sqrt n}\sim t_{n-1}. \tag{3.3}$$
*Proof.* $T = \dfrac{(\bar X-\mu)/(\sigma/\sqrt n)}{\sqrt{[(n-1)S^2/\sigma^2]/(n-1)}}$. The numerator is $\mathcal N(0,1)$ and the bracket is $\chi^2_{n-1}$ by (3.1); they are independent by Theorem 3.1(c). This is (3.2) with $k = n-1$. ∎ The unknown σ cancels; that is the point. For two independent Normal samples the same argument gives $(S_1^2/\sigma_1^2)/(S_2^2/\sigma_2^2)\sim F_{n_1-1,n_2-1}$, the pivot for a variance ratio.

### 3.2 Confidence intervals from pivots

**Definition 3.4.** A **pivot** is a function $g(\mathbf X,\theta)$ whose distribution under $P_\theta$ is the same for every θ. The **coverage** of a set-valued procedure $C$ is the function $\theta\mapsto P_\theta(\theta\in C(\mathbf X))$; $C$ is a $(1-\alpha)$ confidence set if its coverage is at least $1-\alpha$ everywhere. Coverage is the criterion an interval is judged by, and you can check it by simulation.

The recipe: find $a,b$ with $P(a\le g(\mathbf X,\theta)\le b) = 1-\alpha$ and solve the inequality for θ. Applied to (3.1) and (3.3):
$$\sigma\text{ known: }\bar X\pm z_{\alpha/2}\frac{\sigma}{\sqrt n};\qquad \sigma\text{ unknown: }\bar X\pm t_{n-1,\alpha/2}\frac{S}{\sqrt n};\qquad \sigma^2\in\left[\frac{(n-1)S^2}{\chi^2_{n-1,\alpha/2}},\frac{(n-1)S^2}{\chi^2_{n-1,1-\alpha/2}}\right]. \tag{3.4}$$
With the upper-tail convention, $\chi^2_{n-1,\alpha/2}$ is the large quantile and gives the lower endpoint. The variance interval is lopsided because $\chi^2$ is skewed.

**What a 95% interval means.** Said correctly: *before the data are collected, the procedure has probability 0.95 of producing an interval that contains the true μ.* Once you have computed $[9.776, 9.836]$, the sentence "$P(9.776\le\mu\le9.836) = 0.95$" is false in this framework: nothing in it is random. Also wrong: "95% of the data lie in the interval" (that is a tolerance interval, much wider); "a replication's mean lands in this interval with probability 0.95" (about 83% under ideal conditions, because both means are random); "values inside are equally plausible." Right and useful: "the interval is the set of $\mu_0$ a level-0.05 test would not reject" (Theorem 3.5).

**A proportion: Wald and Wilson.** Let $X\sim\mathrm{Bin}(n,p)$ and $\hat p = X/n$. The CLT (F5 Theorem 3.35) gives the approximate pivot $(\hat p-p)/\sqrt{p(1-p)/n}\approx\mathcal N(0,1)$. **Wald** replaces $p$ in the denominator by $\hat p$, then solves a linear inequality:
$$\hat p\pm z_{\alpha/2}\sqrt{\hat p(1-\hat p)/n}. \tag{3.5}$$
**Wilson** (1927) keeps $p$ in the denominator and solves $(\hat p-p)^2\le z^2p(1-p)/n$ for $p$, with $z = z_{\alpha/2}$.

▶ **Try it first.** Write the inequality as $Ap^2-Bp+C\le0$ and use the quadratic formula.

1. Expand: $p^2(1+z^2/n)-p(2\hat p+z^2/n)+\hat p^2\le0$. An upward parabola is $\le0$ between its roots.
2. Discriminant: $(2\hat p+z^2/n)^2-4(1+z^2/n)\hat p^2 = 4\hat pz^2/n+z^4/n^2-4\hat p^2z^2/n = 4z^2\big[\hat p(1-\hat p)/n+z^2/(4n^2)\big]$.
3. Hence
$$p\in\frac{\hat p+\dfrac{z^2}{2n}\pm z\sqrt{\dfrac{\hat p(1-\hat p)}{n}+\dfrac{z^2}{4n^2}}}{1+z^2/n}. \tag{3.6}$$
The project checks (3.6) against SymPy's roots.

■ **What just happened.** Wilson inverted the test exactly; Wald linearized first. Wilson's centre $(X+z^2/2)/(n+z^2)$ is pulled toward 1/2, as if about two successes and two failures were added, and the interval never leaves $[0,1]$. Wald collapses to the point $\{0\}$ when $X = 0$ (zero estimated variance) and spills below 0 when $X$ is small. Example 4.3 measures the damage: at $n = 20$, $p = 0.05$, the "95%" Wald interval covers 64% of the time. Brown, Cai and DasGupta (*Statistical Science*, 2001) showed that Wald's coverage oscillates erratically even when $np$ is large; the project plots it.

**Two means.** For independent samples $X_1..X_{n_1}\sim\mathcal N(\mu_1,\sigma^2)$ and $Y_1..Y_{n_2}\sim\mathcal N(\mu_2,\sigma^2)$ with a *common* σ², $\bar X-\bar Y\sim\mathcal N(\mu_1-\mu_2,\sigma^2(1/n_1+1/n_2))$, and $[(n_1-1)S_1^2+(n_2-1)S_2^2]/\sigma^2\sim\chi^2_{n_1+n_2-2}$ (Theorem 3.1 for each sample; independent $\chi^2$'s add). Dividing as in Corollary 3.3, with pooled variance $S_p^2 = [(n_1-1)S_1^2+(n_2-1)S_2^2]/(n_1+n_2-2)$:
$$\bar X-\bar Y\pm t_{n_1+n_2-2,\alpha/2}\,S_p\sqrt{1/n_1+1/n_2}. \tag{3.7}$$
With unequal variances there is no exact pivot (the Behrens–Fisher problem). **Welch's** interval (1947) is $\bar X-\bar Y\pm t_{\nu,\alpha/2}\sqrt{S_1^2/n_1+S_2^2/n_2}$, with $\nu = (S_1^2/n_1+S_2^2/n_2)^2\big/\big[\frac{(S_1^2/n_1)^2}{n_1-1}+\frac{(S_2^2/n_2)^2}{n_2-1}\big]$ matching two moments of the variance estimate to a scaled $\chi^2_\nu$ (Satterthwaite). Use Welch by default: it loses little when the variances are equal, while the pooled interval miscovers badly when they differ and $n_1\ne n_2$.

**Asymptotic (Wald) intervals.** ST-1's asymptotic normality $\sqrt n(\hat\theta-\theta)\xrightarrow{d}\mathcal N(0,1/I_1(\theta))$ makes $(\hat\theta-\theta)/\widehat{\rm se}$ an approximate pivot, with $\widehat{\rm se} = 1/\sqrt{nI_1(\hat\theta)}$, hence the interval $\hat\theta\pm z_{\alpha/2}\widehat{\rm se}$. Formula (3.5) is this interval for the Bernoulli. It fails near the boundary of the parameter space (poor Normal approximation; the interval escapes the space), it ties the width to the centre's error by plugging $\hat\theta$ into the standard error, and it depends on the parametrization (an interval for $\log\lambda$, exponentiated, differs from one for λ). Inverting the score test (as Wilson does) or the likelihood-ratio test (§3.6) usually does better.

### 3.3 Intervals and tests are the same object

A level-α test of $H_0:\theta = \theta_0$ is described by its **acceptance region** $A(\theta_0)$, the data for which it does not reject, with $P_{\theta_0}(\mathbf X\notin A(\theta_0))\le\alpha$.

**Theorem 3.5 (Duality).** (a) If for each $\theta_0$ there is a level-α test with acceptance region $A(\theta_0)$, then $C(\mathbf x) = \{\theta_0:\mathbf x\in A(\theta_0)\}$ is a $(1-\alpha)$ confidence set. (b) Conversely, if $C$ is a $(1-\alpha)$ confidence set, "reject $H_0:\theta = \theta_0$ iff $\theta_0\notin C(\mathbf X)$" is a level-α test.

*Proof.* (a) By construction $\theta\in C(\mathbf X)$ iff $\mathbf X\in A(\theta)$, so $P_\theta(\theta\in C(\mathbf X)) = P_\theta(\mathbf X\in A(\theta))\ge1-\alpha$ for every θ. (b) $P_{\theta_0}(\text{reject}) = P_{\theta_0}(\theta_0\notin C(\mathbf X))\le\alpha$ by the coverage of $C$ at $\theta_0$. ∎

The $t$ interval is the set of $\mu_0$ the two-sided $t$ test does not reject, and Wilson is the inverted score test: "is $\mu_0$ in the 95% interval?" and "is $p\ge0.05$?" are one question.

### 3.4 The bootstrap and the permutation test

**The plug-in principle.** Many parameters are functionals of the distribution, $\theta = T(F)$: the mean $\int x\,dF$, the median $F^{-1}(1/2)$, a correlation. The plug-in estimate replaces $F$ by $\hat F_n$, which puts mass $1/n$ on each observation: $\hat\theta = T(\hat F_n)$. The standard error of $\hat\theta$ is also a functional of $F$: ${\rm se}_F(\hat\theta) = \sqrt{\operatorname{Var}_F\hat\theta(X_1,\dots,X_n)}$ with the $X_i$ iid $F$. Plug in $\hat F_n$:
$$\widehat{\rm se}_{\rm boot} = \sqrt{\operatorname{Var}_{\hat F_n}\hat\theta(X_1^*,\dots,X_n^*)},\qquad X_i^*\text{ iid }\hat F_n. \tag{3.8}$$
Drawing iid from $\hat F_n$ is sampling the data with replacement. The variance in (3.8) is a finite sum over $n^n$ equally likely ordered resamples, far too many to enumerate, so estimate it by Monte Carlo.

**Algorithm (bootstrap standard error).** For $b = 1,\dots,B$, draw $X^*_{b1},\dots,X^*_{bn}$ with replacement from the data and compute $\hat\theta^*_b$. Return
$$\widehat{\rm se}_B = \Big(\tfrac{1}{B-1}\textstyle\sum_b(\hat\theta^*_b-\bar\theta^*)^2\Big)^{1/2}. \tag{3.9}$$

Keep two errors apart. The *Monte Carlo error* $\widehat{\rm se}_B-\widehat{\rm se}_{\rm boot}$ shrinks like $1/\sqrt B$ and is yours to control ($B = 2000$ is ample for a standard error). The *statistical error* $\widehat{\rm se}_{\rm boot}-{\rm se}_F$ comes from $\hat F_n\ne F$. Why it is small: the Glivenko–Cantelli theorem says $\sup_x|\hat F_n(x)-F(x)|\to0$ almost surely, so if ${\rm se}_F(\hat\theta)$ depends continuously on $F$ (true for smooth functionals such as means, variances, correlations, and quantiles of continuous distributions), plugging in a nearby $\hat F_n$ gives a nearby answer. The precise condition is a differentiability of $T$ in $F$; van der Vaart's *Asymptotic Statistics*, ch. 23, has the theory. For the mean you can check everything by hand: the bootstrap variance of $\bar X^*$ is exactly $\hat\sigma^2/n$ (exercise B5).

**Bootstrap intervals.** With $\hat\theta^*_{(a)}$ the $a$-quantile of the replicates:
$$\text{normal: }\hat\theta\pm z_{\alpha/2}\widehat{\rm se}_B;\qquad\text{percentile: }\big[\hat\theta^*_{(\alpha/2)},\,\hat\theta^*_{(1-\alpha/2)}\big];\qquad\text{basic: }\big[2\hat\theta-\hat\theta^*_{(1-\alpha/2)},\,2\hat\theta-\hat\theta^*_{(\alpha/2)}\big]. \tag{3.10}$$
The basic (pivotal) interval needs a derivation. Treat $\Delta = \hat\theta-\theta$ as an approximate pivot and approximate its distribution by that of $\Delta^* = \hat\theta^*-\hat\theta$, the same quantity in the bootstrap world where $\hat\theta$ plays the truth. (1) In the bootstrap world, $P^*(\hat\theta^*_{(\alpha/2)}-\hat\theta\le\Delta^*\le\hat\theta^*_{(1-\alpha/2)}-\hat\theta) = 1-\alpha$. (2) Transfer: $P(\hat\theta^*_{(\alpha/2)}-\hat\theta\le\hat\theta-\theta\le\hat\theta^*_{(1-\alpha/2)}-\hat\theta)\approx1-\alpha$. (3) Solve for θ to get (3.10)'s basic interval. The percentile interval is right when some monotone transform of $\hat\theta$ is a symmetric pivot, and it respects the parameter's range. Neither dominates: in the project (skewed data, $n = 30$) both cover the mean 90–93% of the time, and for the median the basic interval drops to 0.83. **BCa** (Efron 1987; Efron and Tibshirani ch. 14) corrects the percentile endpoints for bias and skewness. The **parametric bootstrap** resamples from a fitted model $p(x;\hat\theta)$ instead of $\hat F_n$.

**When it fails: the maximum.** Let $X_i\sim\mathrm{Unif}(0,\theta)$ and $\hat\theta = X_{(n)}$, the sample maximum (ST-1's MLE here).

▶ **Try it first.** What is the probability that a resample contains the observed maximum? What is the probability that $X_{(n)} = \theta$ in the real world?

1. Real world: $P\big(n(\theta-X_{(n)})/\theta>t\big) = P\big(\text{all }X_i<\theta(1-t/n)\big) = (1-t/n)^n\to e^{-t}$. The scaled error is asymptotically Exponential(1), a continuous law with no atom at 0.
2. Bootstrap world: $\hat\theta^* = \hat\theta$ exactly when at least one of the $n$ draws is the observation $X_{(n)}$. Each draw misses it with probability $1-1/n$, so
$$P^*(\hat\theta^* = \hat\theta) = 1-(1-1/n)^n\to1-e^{-1}\approx0.632. \tag{3.11}$$
3. The bootstrap law of $n(\hat\theta-\hat\theta^*)$ has an atom of mass ≈ 0.632 at 0; the true law of $n(\theta-\hat\theta)$ has none. They do not approach each other.
4. Worse, $\hat\theta^*\le\hat\theta<\theta$ always, so the percentile interval lies entirely below θ: its coverage is exactly **0**. The project gets 0.000 for percentile and 0.875 for basic. ∎

■ **What just happened.** The maximum depends on $F$ only through the edge of its support, and $\hat F_n$ has nothing beyond $X_{(n)}$. Smoothness in $F$ is what the bootstrap needs, and the maximum does not have it. Here an exact pivot exists: $X_{(n)}/\theta\sim\mathrm{Beta}(n,1)$, so $[X_{(n)},X_{(n)}/\alpha^{1/n}]$ has coverage exactly $1-\alpha$ (0.955 in the project's 400 replications). Other failures: **heavy tails** (with infinite variance the bootstrap of the mean is inconsistent), **dependent data** (resampling points destroys the dependence; block bootstraps resample stretches), and **small $n$** (Example 4.4's bootstrap distribution of a median sits on three or four values).

**The jackknife** (Quenouille 1949, Tukey 1958), the bootstrap's ancestor, recomputes $\hat\theta_{(-i)}$ without observation $i$ and estimates ${\rm se}^2$ by $\frac{n-1}{n}\sum_i(\hat\theta_{(-i)}-\bar\theta_{(\cdot)})^2$. It agrees with the bootstrap for smooth statistics and is inconsistent for the median, which takes only two values when one point is deleted.

**Permutation tests.** Two samples $X_1..X_m\sim F$, $Y_1..Y_n\sim G$; test $H_0:F = G$ with a statistic such as $T = \bar X-\bar Y$. Let $N = m+n$.

**Theorem 3.6 (The permutation test is exact).** Let $p$ be the fraction of the $\binom Nm$ relabelings of the pooled data whose statistic is at least as extreme as the observed one. Under $H_0$, $P(p\le\alpha)\le\alpha$ for every $N$, with no distributional assumption.

*Proof.* Under $H_0$ the $N$ pooled variables are iid, so their joint density $\prod_kf(z_k)$ is unchanged by permuting them: they are **exchangeable**. Condition on the multiset of observed values. Every assignment of those values to the $N$ positions is then equally likely, so the observed labeling is a uniform draw from the $\binom Nm$ labelings, and the observed $T$ is a uniform draw from the list of $\binom Nm$ values of $T$. The fraction of the list at least as extreme as a uniformly chosen entry is at most α with probability at most α (at most $\alpha\binom Nm$ entries can have that fraction $\le\alpha$). This holds conditionally on every multiset, hence unconditionally. ∎

The bootstrap resamples *with* replacement to mimic sampling from $F$; the permutation test reshuffles *without* replacement to exploit the null's symmetry (Fisher, *The Design of Experiments*, 1935: if the treatment does nothing, the labels are arbitrary). The project enumerates all 924 splits of two samples of 6 (p = 0.0173) and compares the pooled $t$ test (p = 0.0113).

### 3.5 The Neyman–Pearson framework

A test is a function $\phi(\mathbf x)\in[0,1]$, the probability of rejecting $H_0$ (usually 0 or 1). A type I error rejects a true $H_0$ (probability at most α); a type II error keeps a false one (probability $1-\beta(\theta)$). Shrinking the rejection region lowers α and power together, so at fixed $n$ the Neyman–Pearson move is: fix α, then maximize power. The asymmetry is deliberate: $H_0$ is the claim you need strong evidence to abandon.

**Theorem 3.7 (Neyman–Pearson lemma).** For simple $H_0:f = f_0$ versus $H_1:f = f_1$, let
$$\phi^*(\mathbf x) = \begin{cases}1, & f_1(\mathbf x)>kf_0(\mathbf x),\\ 0, & f_1(\mathbf x)<kf_0(\mathbf x),\end{cases} \tag{3.12}$$
with any value on $\{f_1 = kf_0\}$, and $k\ge0$ chosen so that $\mathbb E_0\phi^* = \alpha$. Then every test φ with $\mathbb E_0\phi\le\alpha$ has $\mathbb E_1\phi\le\mathbb E_1\phi^*$: the likelihood-ratio test is **most powerful** at level α.

▶ **Try it first.** Show that $(\phi^*(\mathbf x)-\phi(\mathbf x))(f_1(\mathbf x)-kf_0(\mathbf x))\ge0$ for every $\mathbf x$, then integrate.

**Derivation.**
1. Where $f_1>kf_0$: $\phi^* = 1\ge\phi$, so both factors are $\ge0$. Where $f_1<kf_0$: $\phi^* = 0\le\phi$, so both factors are $\le0$. Where $f_1 = kf_0$ the second factor is 0. The product is $\ge0$ everywhere.
2. Integrate (or sum) over $\mathbf x$: $\int(\phi^*-\phi)f_1\ge k\int(\phi^*-\phi)f_0$.
3. The left side is $\mathbb E_1\phi^*-\mathbb E_1\phi$. The right side is $k(\alpha-\mathbb E_0\phi)\ge0$ since $k\ge0$ and $\mathbb E_0\phi\le\alpha$.
4. So $\mathbb E_1\phi^*\ge\mathbb E_1\phi$. ∎

■ **What just happened.** Step 1 is §2's budget argument as an inequality: φ* spends its false-alarm budget on exactly the points with the highest $f_1/f_0$. Any other test moves budget from high-ratio points to low-ratio points and loses power. Randomizing on $\{f_1 = kf_0\}$ matters only for discrete data, to hit α exactly.

**Uniformly most powerful one-sided tests.** Let $X_i\sim\mathcal N(\mu,\sigma^2)$ with σ known, $H_0:\mu\le\mu_0$ versus $H_1:\mu>\mu_0$.

**Theorem 3.8 (UMP test for a Normal mean).** The test rejecting when $\bar X>\mu_0+z_\alpha\sigma/\sqrt n$ has size α and is most powerful against every $\mu_1>\mu_0$ at once.

*Proof.* (1) Fix $\mu_1>\mu_0$. Then $\log\dfrac{f_{\mu_1}(\mathbf x)}{f_{\mu_0}(\mathbf x)} = \dfrac{1}{2\sigma^2}\sum_i\big[(x_i-\mu_0)^2-(x_i-\mu_1)^2\big] = \dfrac{n(\mu_1-\mu_0)}{\sigma^2}\bar x-\dfrac{n(\mu_1^2-\mu_0^2)}{2\sigma^2}$, increasing in $\bar x$. So "$f_{\mu_1}/f_{\mu_0}>k$" is "$\bar x>c$". (2) Size α at $\mu_0$ forces $c = \mu_0+z_\alpha\sigma/\sqrt n$, which does not involve $\mu_1$. By Theorem 3.7 this single test is most powerful for $\mu_0$ against every $\mu_1>\mu_0$. (3) Its power function
$$\beta(\mu) = P_\mu(\bar X>c) = 1-\Phi\Big(z_\alpha-\frac{\sqrt n(\mu-\mu_0)}{\sigma}\Big) \tag{3.13}$$
increases in μ, so $\sup_{\mu\le\mu_0}\beta(\mu) = \beta(\mu_0) = \alpha$: size α for the composite null. (4) Any level-α test of the composite null has level α at $\mu_0$, so by (2) it is no more powerful at any $\mu_1>\mu_0$. ∎

The proof only used that $f_{\theta_1}/f_{\theta_0}$ increases in one statistic $T(\mathbf x)$ for every $\theta_1>\theta_0$, a **monotone likelihood ratio**. The **Karlin–Rubin theorem** (Casella and Berger §8.3.2) is the general statement; every one-parameter exponential family with increasing natural parameter qualifies. Against two-sided alternatives no UMP test exists: the best test against $\mu_1>\mu_0$ has almost no power against $\mu_1<\mu_0$.

### 3.6 Likelihood-ratio tests, Wilks, and the trinity

For composite hypotheses $H_0:\theta\in\Theta_0$ inside a model $\theta\in\Theta$, maximize each likelihood:
$$\lambda(\mathbf x) = \frac{\sup_{\theta\in\Theta_0}L(\theta)}{\sup_{\theta\in\Theta}L(\theta)} = \frac{L(\hat\theta_0)}{L(\hat\theta)}\in(0,1], \tag{3.14}$$
and reject for small λ. The threshold needs the null distribution of λ. Sometimes it is exact (Example 4.6: the one-sample $t$ test is an LRT). In general Wilks (1938) supplies it for large $n$.

**Theorem 3.9 (Wilks).** Under ST-1's regularity conditions (true θ interior to Θ, smooth log-likelihood, consistent MLE), for one-dimensional θ and $H_0:\theta = \theta_0$,
$$-2\log\lambda = 2[\ell(\hat\theta)-\ell(\theta_0)]\xrightarrow{d}\chi^2_1. \tag{3.15}$$
For $\theta\in\mathbb R^k$ and $H_0$ fixing all $k$ coordinates, the limit is $\chi^2_k$; if $H_0$ fixes $r$ of them, it is $\chi^2_r$, with $r = \dim\Theta-\dim\Theta_0$.

▶ **Try it first.** Taylor-expand $\ell(\theta_0)$ about $\hat\theta$ to second order. What is $\ell'(\hat\theta)$? Then use ST-1's asymptotic normality of $\sqrt n(\hat\theta-\theta_0)$.

**Derivation (one parameter).**
1. Taylor with Lagrange remainder: $\ell(\theta_0) = \ell(\hat\theta)+\ell'(\hat\theta)(\theta_0-\hat\theta)+\frac12\ell''(\tilde\theta)(\theta_0-\hat\theta)^2$ with $\tilde\theta$ between $\theta_0$ and $\hat\theta$.
2. $\ell'(\hat\theta) = 0$ at an interior maximum, so $-2\log\lambda = -\ell''(\tilde\theta)(\hat\theta-\theta_0)^2$.
3. Regroup: $-2\log\lambda = \big[-\tfrac1n\ell''(\tilde\theta)\big]\cdot\big[\sqrt n(\hat\theta-\theta_0)\big]^2$.
4. First factor: $-\frac1n\ell''(\theta_0) = -\frac1n\sum_i\partial_\theta^2\log p(X_i;\theta_0)\to I_1(\theta_0)$ in probability by the law of large numbers and ST-1's second form of the information. Since $\tilde\theta$ is squeezed between $\theta_0$ and the consistent $\hat\theta$, the same limit holds at $\tilde\theta$ (a regularity condition: $\ell''/n$ is continuous in θ uniformly near $\theta_0$).
5. Second factor: $\sqrt{I_1(\theta_0)}\sqrt n(\hat\theta-\theta_0)\xrightarrow{d}Z\sim\mathcal N(0,1)$ by ST-1.
6. Slutsky's theorem (if $A_n\to a$ in probability and $B_n\xrightarrow{d}B$, then $A_nB_n\xrightarrow{d}aB$; Wasserman §5.2) gives $-2\log\lambda = \frac{-\ell''(\tilde\theta)/n}{I_1(\theta_0)}\big[\sqrt{I_1(\theta_0)}\sqrt n(\hat\theta-\theta_0)\big]^2\xrightarrow{d}Z^2\sim\chi^2_1$. ∎

**$k$ parameters.** Steps 1–3 give $-2\log\lambda\approx\mathbf u^\top\hat J\mathbf u$ with $\mathbf u = \sqrt n(\hat\theta-\theta_0)\approx\mathcal N(\mathbf 0,I_1^{-1})$ and $\hat J\to I_1$, the $k\times k$ information matrix. $I_1$ is symmetric positive definite, so LA-8 Theorem 3.11 gives $I_1 = Q\Lambda Q^\top$ and a symmetric square root $I_1^{1/2} = Q\Lambda^{1/2}Q^\top$. Then $\mathbf v = I_1^{1/2}\mathbf u\approx\mathcal N(\mathbf 0,I_1^{1/2}I_1^{-1}I_1^{1/2}) = \mathcal N(\mathbf 0,I_k)$, and $\mathbf u^\top I_1\mathbf u = \|\mathbf v\|^2\approx\chi^2_k$. In the right coordinates the quadratic form is a squared length, the move of Theorem 3.1. With $r$ restrictions, the statistic is the squared length of $\mathbf v$'s projection onto an $r$-dimensional subspace (Casella and Berger §10.3).

■ **What just happened.** Near its peak every smooth log-likelihood is a parabola with curvature $nI_1$, and the MLE is Normal with variance $1/(nI_1)$. The drop from the peak to $\theta_0$ is ½ × curvature × distance², twice which is a squared standard Normal. The project simulates a Poisson rate test at $n = 30$ and gets null rejection rates within 0.01 of 0.05.

**The trinity.** The same parabola can be measured three ways:
$$\text{LRT: }2[\ell(\hat\theta)-\ell(\theta_0)],\qquad\text{Wald: }W = n I_1(\hat\theta)(\hat\theta-\theta_0)^2,\qquad\text{score (Rao): }S = \frac{\ell'(\theta_0)^2}{nI_1(\theta_0)}. \tag{3.16}$$
The LRT measures the vertical drop, Wald the horizontal distance scaled by curvature at $\hat\theta$, and the score test the slope at $\theta_0$ scaled by curvature there (no fit under $H_1$ needed). They agree to first order. $W$ is step 6's bracket squared with $I_1(\hat\theta)\to I_1(\theta_0)$. For $S$: expand the score, $\ell'(\theta_0) = \ell'(\hat\theta)+\ell''(\tilde\theta)(\theta_0-\hat\theta) = \ell''(\tilde\theta)(\theta_0-\hat\theta)$, so $S = [\ell''(\tilde\theta)/n]^2\,n(\hat\theta-\theta_0)^2/I_1(\theta_0)\approx nI_1(\theta_0)(\hat\theta-\theta_0)^2$. So all three are asymptotically $\chi^2_1$ under $H_0$. At finite $n$ the LRT and score test are invariant to reparametrization and Wald is not (exercise B4: Wald = 2.31 in the rate, 0.98 in the mean, same data). Inverting Wald gives the Wald interval; inverting the score test gives Wilson. For a Poisson mean, $\ell(\lambda) = n\bar x\log\lambda-n\lambda+\text{const}$, $\hat\lambda = \bar x$, $I_1 = 1/\lambda$, so $-2\log\lambda = 2n[\bar x\log(\bar x/\lambda_0)-(\bar x-\lambda_0)]$, $W = n(\bar x-\lambda_0)^2/\bar x$, $S = n(\bar x-\lambda_0)^2/\lambda_0$. The project's sample with $\bar x = 2.167$, $\lambda_0 = 2$, $n = 30$ gives 0.406, 0.385, 0.417.

### 3.7 The χ² tests and ANOVA

**Goodness of fit.** Counts $(N_1,\dots,N_k)\sim\mathrm{Multinomial}(n;p_1,\dots,p_k)$ and $H_0:p_j = p_j^0$ for all $j$. The unrestricted MLE is $\hat p_j = N_j/n$. With expected counts $E_j = np_j^0$,
$$G = -2\log\lambda = 2\sum_jN_j\log\frac{N_j}{E_j}. \tag{3.17}$$
The simplex $\{\sum p_j = 1\}$ has dimension $k-1$ and $H_0$ pins it down, so Wilks gives $G\approx\chi^2_{k-1}$.

▶ **Try it first.** Write $N_j = E_j+\delta_j$ (so $\sum\delta_j = 0$) and expand $G$ to second order in $\delta_j/E_j$.

1. $\log(N_j/E_j) = \delta_j/E_j-\delta_j^2/(2E_j^2)+O(\delta_j^3/E_j^3)$.
2. $N_j\log(N_j/E_j) = (E_j+\delta_j)\big(\delta_j/E_j-\delta_j^2/(2E_j^2)\big)+\cdots = \delta_j+\delta_j^2/(2E_j)+O(\delta_j^3/E_j^2)$.
3. Sum, multiply by 2, and use $\sum\delta_j = 0$:
$$G\approx\sum_j\frac{\delta_j^2}{E_j} = \sum_j\frac{(N_j-E_j)^2}{E_j} =: X^2. \tag{3.18}$$

■ **What just happened.** Pearson's $X^2$ (1900), the formula everyone memorizes, is the quadratic approximation to the likelihood-ratio statistic: the "every log-likelihood is a parabola" fact again. It is also *exactly* the score statistic (exercise C2). Since $\delta_j = O(\sqrt n)$, the neglected terms are $O(1/\sqrt n)$; the rule "every $E_j\ge5$" keeps them small.

**Independence in an $r\times c$ table.** Counts $N_{ij}$ with cell probabilities $p_{ij}$; $H_0:p_{ij} = p_{i\cdot}p_{\cdot j}$. The null MLEs are the margins, so the expected counts are $E_{ij} = N_{i\cdot}N_{\cdot j}/n$ (row total × column total / grand total), and $X^2 = \sum_{ij}(N_{ij}-E_{ij})^2/E_{ij}\approx\chi^2_\nu$ with $\nu = (rc-1)-[(r-1)+(c-1)] = (r-1)(c-1)$: free parameters under $H_1$ minus free parameters under $H_0$.

**One-way ANOVA as a projection.** $k$ groups, $y_{gj} = \mu_g+\varepsilon_{gj}$, $\varepsilon_{gj}$ iid $\mathcal N(0,\sigma^2)$, $n_g$ observations in group $g$, $N = \sum n_g$. Test $H_0:\mu_1 = \cdots = \mu_k$. Stack the data into $\mathbf y\in\mathbb R^N$ and take two nested subspaces: $V_0 = \operatorname{span}(\mathbf 1)$ (dimension 1, the mean vectors $H_0$ allows) and $V_1$ = span of the $k$ group-indicator vectors (dimension $k$, the mean vectors $H_1$ allows). $V_0\subset V_1$ because $\mathbf 1$ is the sum of the indicators. Let $P_0,P_1$ be the orthogonal projections (LA-5 Theorem 3.12): $P_1\mathbf y$ replaces each observation by its group mean and $P_0\mathbf y = \bar y\mathbf 1$. The pieces $(I-P_1)\mathbf y$, $(P_1-P_0)\mathbf y$, $P_0\mathbf y$ lie in the mutually orthogonal subspaces $V_1^\perp$, $V_1\cap V_0^\perp$, $V_0$ of dimensions $N-k$, $k-1$, $1$. Pythagoras gives the **sum-of-squares decomposition**
$$\underbrace{\|\mathbf y-\bar y\mathbf 1\|^2}_{\rm SST} = \underbrace{\textstyle\sum_gn_g(\bar y_g-\bar y)^2}_{{\rm SSB} = \|(P_1-P_0)\mathbf y\|^2}+\underbrace{\textstyle\sum_{g,j}(y_{gj}-\bar y_g)^2}_{{\rm SSW} = \|(I-P_1)\mathbf y\|^2}. \tag{3.19}$$
Under $H_0$, $\mathbf y = \mu\mathbf 1+\sigma\mathbf z$, and $P_1-P_0$ and $I-P_1$ both kill $\mu\mathbf 1$. Repeat Theorem 3.1's steps 2–6 in an orthonormal basis adapted to $V_0\subset V_1\subset\mathbb R^N$: ${\rm SSB}/\sigma^2\sim\chi^2_{k-1}$ and ${\rm SSW}/\sigma^2\sim\chi^2_{N-k}$, independent. Hence
$$F = \frac{{\rm SSB}/(k-1)}{{\rm SSW}/(N-k)}\sim F_{k-1,N-k}\quad\text{under }H_0. \tag{3.20}$$
It is the LRT: maximizing the Normal likelihood over σ² for a fixed mean vector gives $L_{\max}\propto({\rm RSS}/N)^{-N/2}$ (RSS = residual sum of squares), so $\lambda = ({\rm SSW}/{\rm SST})^{N/2}$ and $-2\log\lambda = N\log\big(1+\frac{k-1}{N-k}F\big)$, increasing in $F$. Under $H_1$, SSB picks up $\|(P_1-P_0)\boldsymbol\mu\|^2 = \sum_gn_g(\mu_g-\bar\mu)^2$ ($\bar\mu$ the $n_g$-weighted mean), which inflates $F$. With $k = 2$, $F = t^2$ for the pooled two-sample statistic behind (3.7): the two-sample $t$ test is the two-group LRT. ST-3 replaces the group indicators by any design matrix and gets every regression $F$ test from (3.19)–(3.20) unchanged.

### 3.8 p-values

**Definition 3.10.** For a family of level-α tests that reject for large $T$, the **p-value** is the smallest α at which the observed data would be rejected. With null CDF $F_0$ for $T$ (at the boundary point of a composite null), it is the tail probability $p = P_{H_0}(T\ge t_{\rm obs}) = 1-F_0(t_{\rm obs})$ for continuous $T$.

**Theorem 3.11 (p-values are Uniform under $H_0$).** If $T$ has a continuous, strictly increasing CDF $F_0$ under a simple $H_0$, then $p = 1-F_0(T)\sim\mathrm{Unif}(0,1)$. For discrete $T$, $P_{H_0}(p\le u)\le u$ for every $u$.

*Proof.* By F4 Theorem 3.15(b), $F_0(T)\sim\mathrm{Unif}(0,1)$, and $1-U\sim\mathrm{Unif}(0,1)$ when $U$ is. For discrete $T$, $p\le u$ forces $T$ into an upper tail whose null probability is at most $u$. ∎

So "reject when $p\le\alpha$" has size exactly α, and under $H_0$, $p = 0.8$ and $p = 0.08$ are equally ordinary (the project checks $10^4$ null p-values with a Kolmogorov–Smirnov test).

**What a p-value is not.** Greenland et al. (2016) list 25 misreadings of p-values, intervals and power. The ones that do the most damage:
1. *Not $P(H_0\mid\text{data})$.* It is $P(\text{data at least this extreme}\mid H_0)$. Turning one into the other needs a prior and Bayes' theorem (F4 Theorem 3.6), and the two can differ wildly (§3.11).
2. *Not "the probability the result is due to chance."* Same error: "due to chance" means "$H_0$ is true."
3. *Not the probability that a replication succeeds.* If the true effect equals the observed one and $p = 0.05$, a same-size replication is significant about half the time.
4. *A small p is not a large effect.* With $n = 10^6$, a mean shift of 0.003σ gives $p\approx0.003$. Significance measures detectability, which grows with $n$; importance is the effect size, which does not.
5. *A large p is not evidence for $H_0$.* An underpowered study returns large p-values for real effects most of the time. The confidence interval shows which effects remain compatible with the data.
6. *0.049 and 0.051 are not different findings.* The threshold is a decision device, not a property of nature.

Report an **effect size** in meaningful units with its **confidence interval**. By Theorem 3.5 the interval contains the test's verdict, plus magnitude and precision.

### 3.9 Power, sample size, and the winner's curse

The power function (3.13) is the design tool.

▶ **Try it first.** Set $\beta(\mu_0+\delta) = 1-\beta_{\rm II}$ in (3.13) and solve for $n$.

1. $1-\Phi(z_\alpha-\sqrt n\delta/\sigma) = 1-\beta_{\rm II}$, i.e. $\Phi(z_\alpha-\sqrt n\delta/\sigma) = \beta_{\rm II}$.
2. With upper-tail quantiles, $\Phi(-z_{\beta_{\rm II}}) = \beta_{\rm II}$, so $z_\alpha-\sqrt n\delta/\sigma = -z_{\beta_{\rm II}}$.
3. Solve:
$$n = \left(\frac{(z_\alpha+z_{\beta_{\rm II}})\sigma}{\delta}\right)^2. \tag{3.21}$$
For a two-sided test use $z_{\alpha/2}$ (neglecting the tiny chance of rejecting in the wrong tail). For two groups of $n$ each, ${\rm Var}(\bar X-\bar Y)$ doubles and so does $n$.

■ **What just happened.** Sample size grows like $1/\delta^2$: halving the detectable effect quadruples the experiment. That is why small effects need enormous studies, and why "we found no effect with $n = 20$" usually means "we could not have found one."

**The winner's curse.** Suppose an estimate $\hat\theta\sim\mathcal N(\theta,s^2)$ is published only if significant in the positive direction, $\hat\theta>c = z_{\alpha/2}s$. Published estimates are draws from a truncated Normal. With $a = (c-\theta)/s$:
1. $\mathbb E[\hat\theta\mid\hat\theta>c] = \theta+s\,\mathbb E[Z\mid Z>a]$, $Z\sim\mathcal N(0,1)$.
2. $\mathbb E[Z\mid Z>a] = \dfrac{\int_a^\infty z\varphi(z)\,dz}{1-\Phi(a)} = \dfrac{\varphi(a)}{1-\Phi(a)}$, because $\varphi'(z) = -z\varphi(z)$.
3. Therefore
$$\mathbb E[\hat\theta\mid\text{published}] = \theta+s\,\frac{\varphi(a)}{1-\Phi(a)}. \tag{3.22}$$
Take $\theta = 0.2\sigma$, $n = 20$: $s = 0.224\sigma$, $c = 0.438\sigma$, $a = 1.066$, power about 0.14, and (3.22) gives $0.553\sigma$, **2.8 times the truth**. The project simulates $10^5$ studies with a $t$ test and finds published estimates averaging 2.7× the truth. When power is low, the significance filter passes only the estimates that got lucky, so the literature overstates effects and same-size replications "fail." Gelman and Carlin (2014) call this a type M (magnitude) error; type S (sign) errors, significant results with the wrong sign, are 1.8% of the significant results here.

### 3.10 Multiple testing

Test $H_1,\dots,H_m$ with p-values $p_1,\dots,p_m$; the true nulls form the set $\mathcal I_0$, of size $m_0$. Let $V$ count true nulls rejected (false discoveries) and $R$ all rejections. The **family-wise error rate** is ${\rm FWER} = P(V\ge1)$; the **false discovery rate** is ${\rm FDR} = \mathbb E[V/\max(R,1)]$.

**Theorem 3.12 (Bonferroni).** Rejecting $H_i$ iff $p_i\le\alpha/m$ gives ${\rm FWER}\le\alpha$, whatever the dependence between tests.

*Proof.* ${\rm FWER} = P\big(\bigcup_{i\in\mathcal I_0}\{p_i\le\alpha/m\}\big)\le\sum_{i\in\mathcal I_0}P(p_i\le\alpha/m)\le m_0\alpha/m\le\alpha$, by the union bound (F4 Theorem 3.4) and Theorem 3.11. ∎

**Holm (1979)**, step-down: sort $p_{(1)}\le\cdots\le p_{(m)}$; reject $H_{(1)}$ if $p_{(1)}\le\alpha/m$, then $H_{(2)}$ if $p_{(2)}\le\alpha/(m-1)$, and so on, stopping at the first failure. It controls FWER at α under any dependence (exercise B6 proves it), and since its first threshold is Bonferroni's and the rest are larger, it rejects everything Bonferroni does. There is no reason to use plain Bonferroni except arithmetic.

FWER control fits when a single false claim is costly (a drug approval). In screening (10,000 genes), you can live with some false discoveries if their *proportion* is controlled.

**Theorem 3.13 (Benjamini–Hochberg, 1995).** Let $k = \max\{j:p_{(j)}\le jq/m\}$ and reject $H_{(1)},\dots,H_{(k)}$ (nothing if no $j$ qualifies). If the p-values are independent and each true-null p-value is $\mathrm{Unif}(0,1)$, then ${\rm FDR} = (m_0/m)\,q\le q$.

▶ **Try it first.** Write $V/\max(R,1) = \sum_{i\in\mathcal I_0}\mathbf 1\{H_i\text{ rejected}\}/\max(R,1)$ and handle one true null at a time. What happens to $R$ if you replace $p_i$ by 0?

**Derivation.**
1. Let $N(t) = \#\{j:p_j\le t\}$. Then $p_{(j)}\le jq/m$ iff $N(jq/m)\ge j$, so $R = \max\{j:N(jq/m)\ge j\}$, and BH rejects exactly the $H_i$ with $p_i\le Rq/m$.
2. Fix a true null $i$ and let $R^{(i)}$ be the number of rejections BH makes with $p_i$ replaced by 0. $R^{(i)}$ depends only on the other p-values, so it is independent of $p_i$, and $R^{(i)}\ge1$.
3. *Claim:* if $p_i\le kq/m$, then $R = k$ iff $R^{(i)} = k$. Replacing $p_i$ by 0 changes $N(t)$ only for $t<p_i$. For $j\ge k$, $jq/m\ge kq/m\ge p_i$, so $N(jq/m)$ is unchanged. "$R = k$" means $N(kq/m)\ge k$ and $N(jq/m)<j$ for all $j>k$, conditions involving only $j\ge k$.
4. By steps 1 and 3, $\{H_i\text{ rejected},R = k\} = \{p_i\le kq/m,\,R^{(i)} = k\}$. So
$$\mathbb E\Big[\frac{\mathbf 1\{H_i\text{ rej.}\}}{\max(R,1)}\Big] = \sum_{k=1}^m\frac1kP\Big(p_i\le\frac{kq}{m}\Big)P(R^{(i)} = k) = \sum_{k=1}^m\frac1k\cdot\frac{kq}{m}P(R^{(i)} = k) = \frac qm,$$
using independence, $P(p_i\le u) = u$, and $\sum_kP(R^{(i)} = k) = 1$.
5. Sum over the $m_0$ true nulls: ${\rm FDR} = m_0q/m\le q$. ∎

■ **What just happened.** Each true null contributes exactly $q/m$, because the step-up threshold $kq/m$ grows in proportion to the number $k$ of rejections it gets divided by. Independence was used once, in step 4. BH still controls FDR under positive dependence (Benjamini and Yekutieli 2001); under arbitrary dependence, replace $q$ by $q/\sum_{j\le m}j^{-1}$. In the project (80 true nulls, 20 real effects), BH has FDR 0.041 (theory 0.040) and power 0.72 against Holm's 0.40, while its FWER is 0.46: it buys power by accepting some false discoveries.

**The look-elsewhere effect.** A physicist scanning a mass spectrum in 100 bins runs 100 tests. Under background only, a local $3\sigma$ excess ($p = 1.35\times10^{-3}$, one-sided) somewhere among 100 independent bins has probability $1-(1-1.35\times10^{-3})^{100} = 0.126$. Particle physics adopted the $5\sigma$ discovery threshold ($p = 2.87\times10^{-7}$) partly because $3\sigma$ effects kept evaporating (Lyons, "Discovering the significance of 5σ," 2013), and the 2012 Higgs announcements quoted global significances, corrected for this effect, alongside the local ones.

**p-hacking and forking paths.** An analyst facing a true null who measures 20 independent outcomes and reports whichever is significant has false-positive rate $1-0.95^{20} = 0.64$; one who tries 5 subgroups, 0.23 (the project simulates both). Simmons, Nelson and Simonsohn (2011) showed that four common "researcher degrees of freedom" (optional stopping, an extra outcome, a covariate, dropping a condition) together push the false-positive rate above 60%. Gelman and Loken's "garden of forking paths" (2013) adds that no conscious fishing is needed: the analysis you *would* have run depends on the data you saw, while the p-value assumed a fixed analysis. **Pre-registration** (publicly fixing hypotheses, outcomes, sample size and analysis before collecting data) restores that assumption; registered reports, accepted before results exist, also remove §3.9's publication filter. Exploration stays legitimate when labelled as such.

### 3.11 Bayesian testing in one page (ST-5 develops it)

A Bayesian compares hypotheses by the **Bayes factor**
$$B_{01} = \frac{p(\mathbf x\mid H_0)}{p(\mathbf x\mid H_1)} = \frac{p(\mathbf x\mid\theta_0)}{\int p(\mathbf x\mid\theta)\pi(\theta)\,d\theta}, \tag{3.23}$$
the ratio of how well each hypothesis predicted the data, with $H_1$'s prediction averaged over its prior π; posterior odds = $B_{01}$ × prior odds (F4 Theorem 3.6). For $\bar X\sim\mathcal N(\mu,1/n)$, $H_0:\mu = 0$, $H_1:\mu\sim\mathcal N(0,\tau^2)$, the marginal of $\bar X$ under $H_1$ is $\mathcal N(0,\tau^2+1/n)$, and the ratio of the two Normal densities at $\bar x = z/\sqrt n$ is
$$B_{01} = \sqrt{1+n\tau^2}\,\exp\Big(-\frac{z^2}{2}\cdot\frac{n\tau^2}{1+n\tau^2}\Big).$$
**The Jeffreys–Lindley paradox** (Jeffreys 1939, Lindley 1957): hold $z$, and so the p-value, fixed and let $n\to\infty$. The exponential tends to $e^{-z^2/2}$ but the prefactor grows like $\sqrt n$, so $B_{01}\to\infty$: data that reject $H_0$ at any fixed α end up favouring $H_0$. At $z = 2.5$ ($p = 0.012$) and $\tau = 1$, $B_{01}$ rises from 0.19 at $n = 10$ to 44 at $n = 10^6$ (exercise C4). A fixed $z$ at large $n$ means a tiny $\bar x\approx z/\sqrt n$, which the point null predicts far better than a diffuse alternative. The p-value ignores what $H_1$ predicted; the Bayes factor depends on it, and on τ. They answer different questions.

### 3.12 How to report an analysis

1. State the question, the estimand, and (if testing) $H_0$, $H_1$ and α before seeing outcome data; say whether the analysis was pre-registered.
2. Lead with the estimate and its confidence interval, in meaningful units.
3. Give the exact p-value, the test, its assumptions, and how you checked them.
4. Report the sample size and the power calculation behind it, including the design effect δ.
5. Report every test run, not only the significant ones, and the multiplicity correction (FWER or FDR, and why).
6. Never write "no effect" because $p>0.05$; write the interval.
7. Label confirmatory and exploratory results separately.

## 4. Worked examples

### Example 4.1: The projection proof for $n = 3$, with numbers

Take $\mathbf x = (2,5,8)^\top$: $\bar x = 5$ and $(n-1)s^2 = 9+0+9 = 18$. An orthonormal basis adapted to $\mathbf 1$: $\mathbf q_1 = \frac{1}{\sqrt3}(1,1,1)^\top$, $\mathbf q_2 = \frac{1}{\sqrt2}(1,-1,0)^\top$, $\mathbf q_3 = \frac{1}{\sqrt6}(1,1,-2)^\top$ (each has length 1; each pair has dot product 0). The coordinates $w_i = \mathbf q_i^\top\mathbf x$ are $w_1 = 15/\sqrt3 = 8.660 = \sqrt3\,\bar x$, $w_2 = -3/\sqrt2 = -2.121$, $w_3 = (2+5-16)/\sqrt6 = -3.674$. Then $w_2^2+w_3^2 = 4.5+13.5 = 18 = (n-1)s^2$. ✓ The mean lives in coordinate 1, the spread in coordinates 2 and 3. For $\mathbf X\sim\mathcal N(\mu\mathbf 1,\sigma^2I_3)$, $W_2/\sigma$, $W_3/\sigma$ are iid standard Normals independent of $W_1$, so $2S^2/\sigma^2 = (W_2^2+W_3^2)/\sigma^2\sim\chi^2_2$, and
$$T = \frac{\bar X-\mu}{S/\sqrt3} = \frac{(W_1-\sqrt3\mu)/\sigma}{\sqrt{(W_2^2+W_3^2)/(2\sigma^2)}}$$
is a standard Normal over the root-mean-square of two independent standard Normals: $t_2$ by (3.2), with σ cancelled.

### Example 4.2: A $t$ interval for ten measurements of $g$

A pendulum lab gives ten measurements of $g$ (m/s²): 9.79, 9.83, 9.76, 9.85, 9.81, 9.78, 9.88, 9.80, 9.74, 9.82.
- $\bar x = 98.06/10 = 9.806$; $\sum(x_i-\bar x)^2 = 0.01564$; $s^2 = 0.01564/9 = 0.0017378$; $s = 0.04169$; standard error $s/\sqrt{10} = 0.01318$.
- $t_{9,0.025} = 2.2622$, so the half-width is $2.2622\times0.01318 = 0.0298$.
- 95% interval for μ: $[9.7762,\ 9.8358]$ m/s².

It came from a procedure that covers the true $g$ in 95% of repetitions; that is all the 95% says. Using 1.96 in place of 2.262 gives the "z with $S$" interval, whose coverage at $n = 10$ the project measures as 0.916.

For σ, (3.4) with $\chi^2_{9,0.025} = 19.023$ and $\chi^2_{9,0.975} = 2.700$ gives $\sigma^2\in[0.01564/19.023,\ 0.01564/2.700]$, so $\sigma\in[0.0287,\ 0.0761]$ m/s². Ten points pin down a mean far better than a spread.

### Example 4.3: Wald versus Wilson for 2 successes in 20 trials

$\hat p = 0.1$, $z = 1.960$.
- *Wald* (3.5): ${\rm se} = \sqrt{0.1\times0.9/20} = 0.0671$; interval $0.1\pm0.1315 = [-0.0315,\ 0.2315]$. The lower end is impossible.
- *Wilson* (3.6): $z^2/n = 0.1921$; centre $(0.1+0.0960)/1.1921 = 0.1645$; half-width $1.96\sqrt{0.0045+0.0024}/1.1921 = 0.1365$; interval $[0.0279,\ 0.3010]$.

*Coverage by enumeration.* At a true $p$, coverage is the finite sum $C(p) = \sum_{x=0}^{20}\binom{20}{x}p^x(1-p)^{20-x}\mathbf 1\{p\in{\rm CI}(x)\}$; no simulation needed. At $p = 0.10$: Wald 0.876, Wilson 0.957. At $p = 0.05$, Wald contains 0.05 only for $x\in\{1,2,3,4\}$ (at $x = 0$ the interval is $\{0\}$; at $x = 5$ it is $[0.060,0.440]$), so $C = 1-P(X = 0)-P(X\ge5) = 0.639$; Wilson contains it for $x\in\{0,1,2\}$, giving 0.925. Averaged over $p\in[0.01,0.99]$ at $n = 20$, Wilson's coverage is 0.954 and Wald's 0.858. At $p = 0.05$ the "95%" Wald interval is a 64% interval.

### Example 4.4: Bootstrap standard error of a median of 15

Repair times (hours): 3.1, 4.7, 2.2, 8.9, 5.0, 3.8, 12.4, 4.1, 6.3, 2.9, 5.6, 4.4, 7.2, 3.5, 19.8. Sorted, the 8th value is the median, 4.7. The asymptotic formula for a median's standard error needs the unknown density at the median, so bootstrap it.

*Algorithm output.* $B = 10^4$ resamples of size 15 (seed 0 in the project) give $\widehat{\rm se}_B = 0.802$.

*The ideal bootstrap, exactly.* A resample's median is $\le x_{(k)}$ iff at least 8 of its 15 draws are $\le x_{(k)}$, and each draw is with probability $k/15$. So $P^*(\text{median}^*\le x_{(k)}) = P(\mathrm{Bin}(15,k/15)\ge8)$; differencing in $k$ gives the exact bootstrap distribution, whose standard deviation, the $B = \infty$ answer, is 0.798. The Monte Carlo error at $B = 10^4$ is under 1%.

*What to notice.* The bootstrap distribution sits on a few values: median* = 4.7 with probability 0.207, and 4.4 or 5.0 with probability 0.183 each (`figures/bootstrap_median.png` is spiky). The standard error is usable; a percentile interval from it is crude at $n = 15$. For contrast, the mean's standard error $s/\sqrt n = 1.181$ is inflated by the 19.8, which barely touches the median.

### Example 4.5: Neyman–Pearson for a Normal mean

$X_1,\dots,X_{16}\sim\mathcal N(\mu,2^2)$; $H_0:\mu = 0$ versus $H_1:\mu = 1$; α = 0.05.
- Likelihood ratio (Theorem 3.8, step 1, with $n = 16$, $\sigma^2 = 4$): $\log\Lambda = 4\bar x-2$, increasing in $\bar x$. Reject for $\bar x>c$.
- Size: under $H_0$, $\bar X\sim\mathcal N(0,0.5^2)$, so $c = z_{0.05}\times0.5 = 1.6449\times0.5 = 0.8224$ (equivalently $k = e^{4(0.8224)-2} = 3.73$).
- Power: under $H_1$, $\bar X\sim\mathcal N(1,0.5^2)$, so $\beta(1) = 1-\Phi\big((0.8224-1)/0.5\big) = \Phi(0.355) = 0.639$.

No level-0.05 test of these hypotheses has power above 0.639. The type II error rate is 0.361: this design misses a half-σ effect more than a third of the time.

### Example 4.6: The one-sample $t$ test is a likelihood-ratio test

$X_i\sim\mathcal N(\mu,\sigma^2)$, both unknown; $H_0:\mu = \mu_0$ (σ² free) versus μ free.
1. For fixed μ, maximizing $\ell(\mu,\sigma^2) = -\frac n2\log(2\pi\sigma^2)-\frac{1}{2\sigma^2}\sum(x_i-\mu)^2$ over σ² gives $\hat\sigma^2(\mu) = \frac1n\sum(x_i-\mu)^2$ and $L_{\max}(\mu) = (2\pi e\,\hat\sigma^2(\mu))^{-n/2}$.
2. Under $H_1$, $\hat\mu = \bar x$; under $H_0$, $\mu = \mu_0$. So $\lambda = \big(\hat\sigma^2(\bar x)/\hat\sigma^2(\mu_0)\big)^{n/2} = \Big(\dfrac{\sum(x_i-\bar x)^2}{\sum(x_i-\mu_0)^2}\Big)^{n/2}$.
3. Pythagoras along $\mathbf 1$ again: $\sum(x_i-\mu_0)^2 = \sum(x_i-\bar x)^2+n(\bar x-\mu_0)^2$.
4. Divide by $\sum(x_i-\bar x)^2 = (n-1)s^2$: $\lambda^{-2/n} = 1+\dfrac{n(\bar x-\mu_0)^2}{(n-1)s^2} = 1+\dfrac{t^2}{n-1}$, with $t = (\bar x-\mu_0)/(s/\sqrt n)$. So
$$\lambda = \Big(1+\frac{t^2}{n-1}\Big)^{-n/2}. \tag{4.1}$$
λ decreases in $t^2$: "reject for small λ" is "reject for large $|t|$." The LRT *is* the two-sided $t$ test, and its exact null distribution comes from Corollary 3.3, not from Wilks.

*Numbers.* The $g$ data of Example 4.2 with $\mu_0 = 9.78$: $t = 0.026/0.01318 = 1.972$, $\lambda = (1+3.890/9)^{-5} = 0.166$, $-2\log\lambda = 10\log1.4322 = 3.592$. The exact $t_9$ p-value is 0.080; Wilks' $\chi^2_1$ approximation gives 0.058. At $n = 10$ the asymptotic answer is noticeably anti-conservative; use the exact one when you have it.

### Example 4.7: A $2\times2$ table

A trial puts 50 patients on treatment and 50 on control; 30 and 18 improve.

| | improved | not | total |
|---|---|---|---|
| treatment | 30 | 20 | 50 |
| control | 18 | 32 | 50 |
| total | 48 | 52 | 100 |

Expected counts under independence, $E_{ij} = (\text{row total})(\text{column total})/100$: 24 and 26 in each row. Every cell is off by 6, so
$X^2 = 36\big(\tfrac1{24}+\tfrac1{26}+\tfrac1{24}+\tfrac1{26}\big) = 36(0.08333+0.07692) = 5.769$ on $(2-1)(2-1) = 1$ degree of freedom: $p = 0.0163$. The LRT form (3.17) gives $G = 5.826$, $p = 0.0158$. $X^2$ also equals the square of the pooled two-proportion statistic $z = (0.60-0.36)/\sqrt{0.48\times0.52\times(2/50)} = 2.402$.

The report that matters: the improvement rate is higher by 0.24 under treatment, with 95% Wald interval $0.24\pm1.96\sqrt{0.6\times0.4/50+0.36\times0.64/50} = [0.050,\ 0.430]$ (odds ratio 2.67). The data are compatible with anything from a marginal benefit to a very large one, which "$p = 0.016$" alone hides.

### Example 4.8: Sample size to detect half a standard deviation

Two-sided α = 0.05, power 0.80 at $\delta = 0.5\sigma$. By (3.21) with $z_{0.025} = 1.960$ and $z_{0.20} = 0.8416$: $n = \big((1.960+0.8416)/0.5\big)^2 = 5.603^2 = 31.4$, so $n = 32$. Check: the $z$ test's power at $n = 32$ is 0.807. Since σ will be estimated, the real test is a $t$ test, with power from the noncentral $t$ distribution: 0.782 at $n = 32$; the smallest $n$ reaching 0.80 is 34 (power 0.808; the project's $2\times10^4$ simulated studies give 0.804). The formula is a good first answer; add a few observations for the $t$ correction and more for dropouts.

### Example 4.9: Bonferroni, Holm and BH on ten p-values

Sorted: 0.001, 0.004, 0.006, 0.012, 0.019, 0.033, 0.047, 0.21, 0.38, 0.73; $m = 10$, α = $q$ = 0.05.
- *Bonferroni*: threshold 0.005. Rejects 0.001, 0.004: **2**.
- *Holm*: $0.001\le0.05/10$ ✓; $0.004\le0.05/9 = 0.00556$ ✓; $0.006\le0.05/8 = 0.00625$ ✓; $0.012\le0.05/7 = 0.00714$ ✗, stop: **3**.
- *BH*: compare $p_{(j)}$ with $0.005j$: $0.001\le0.005$, $0.004\le0.010$, $0.006\le0.015$, $0.012\le0.020$, $0.019\le0.025$ all pass; $0.033>0.030$, $0.047>0.035$, and the rest fail. Largest passing $j = 5$: **5**.

Seven p-values fall below 0.05. They yield 2, 3 or 5 discoveries depending on whether you control the chance of *any* false discovery or their expected *proportion*. Choose before you look.

## 5. Why it matters and where it connects

- **ST-3**: coefficient $t$ statistics with $n-p$ degrees of freedom (Theorem 3.1 with $p$ directions projected out), nested-model $F$ tests from (3.19)–(3.20) with a general design matrix, and bootstrap intervals for coefficients (exercise B8).
- **ST-4** tests GLM coefficients with the trinity of §3.6: Wald $z$ statistics in every software printout, likelihood-ratio tests as differences of deviances (Wilks), score tests for adding terms.
- **ST-5** reopens §3.11: credible versus confidence intervals, Bayes factors, the Jeffreys–Lindley paradox in full, and why priors matter most for testing.
- **ST-7**'s Dickey–Fuller unit-root statistic is a $t$ statistic whose null distribution is *not* $t$: under the unit-root null the regressor is nonstationary and the CLT behind Wilks and the trinity fails.
- **ML-1**'s optimism of the best cross-validation score among many models is the winner's curse and the look-elsewhere effect of §3.9–3.10; nested cross-validation is its multiplicity correction.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 drill the pivot intervals, Wilson at $x = 0$, size and power of a discrete test, the three multiplicity procedures, the interpretation of an interval and a p-value, and a sample-size calculation. B1–B8 are the module: B1 (the $n = 2$ projection proof), B3 (Neyman–Pearson and UMP for an Exponential rate), B4 (the trinity and Wald's non-invariance), B6 (Holm controls FWER) and B8 (bootstrap a regression slope in `project.py`) are essential. C1–C4 go further: the $t$ density, Pearson's $X^2$ as an exact score statistic, optional stopping, and the Jeffreys–Lindley paradox in numbers.

## 7. Coding project

`project.py` checks every promise in this module by brute force, in nine parts: (1) Example 4.1's rotation and simulated $\chi^2_4$, $t_4$ against SciPy, with $\operatorname{corr}(\bar X,S^2)\approx0$ for Normal data and 0.69 for Exponential; (2) coverage of $z$, $t$ and "z with $S$" intervals ($10^4$ replications) and of Wald and Wilson by exact enumeration, with coverage curves over $p$; (3) a bootstrap library (`bootstrap_se`, `bootstrap_ci` with normal, percentile, basic) on the mean, median and correlation with coverage by simulation, Example 4.4 against the exact ideal bootstrap, and the maximum's 0.632 atom and zero coverage; (4) a permutation test over all 924 relabelings; (5) Example 4.5, power curves from the noncentral $t$, Example 4.8; (6) Example 4.6, Wilks for a Poisson rate (QQ plot, SymPy expansion) and the trinity; (7) Example 4.7 and ANOVA from explicit projection matrices, with traces equal to degrees of freedom and $F$ equal to `scipy.stats.f_oneway`; (8) uniform null p-values, the winner's curse, p-hacking; (9) Example 4.9 and FWER/FDR/power of Bonferroni, Holm and BH on 2,000 simulated families.

**Build it yourself:** before reading the code, write `wilson_prop` from (3.6), `bootstrap_ci` from (3.10), and `bh` from Theorem 3.13, about ten lines each. "Your turn" asks for a bootstrapped regression slope (B8), optional stopping (C3), a look-elsewhere bump hunt, and heavy-tailed data; reference solutions for B8 and C3 are at the bottom.

## 8. Where you'll struggle

- **The sentence about the interval.** "There is a 95% probability that μ is in [9.776, 9.836]" is the single most common misstatement in science, and you will catch yourself writing it. The fix is mechanical: for every interval in this module, say aloud what is random (the endpoints) and what is fixed (μ). The Bayesian sentence you want to say needs a prior and lives in ST-5.
- **Fisher versus Neyman–Pearson.** Fisher: the p-value is graded evidence, no alternative needed, no fixed α. Neyman–Pearson: fix α and power in advance, then the only output is reject or not, and the size of $p$ below α carries no extra meaning. The textbook hybrid (report $p$ to three decimals *and* star it against 0.05) is incoherent. Know which one you are doing; for a confirmatory study it should be Neyman–Pearson with a power calculation.
- **The $t$-test-as-LRT algebra** (Example 4.6) is skipped by most courses. Do it once by hand. After it, Wilks, the trinity and ANOVA's $F$ all look like the same calculation, because they are.
- **Upper-tail quantiles.** $\chi^2_{9,0.025} = 19.02$ is the large one in this module; software uses the CDF (`chi2.ppf(0.975, 9)` = 19.02). Mixing the conventions swaps the endpoints of (3.4). Write the tail probability next to every quantile.
- **The bootstrap is not magic.** It estimates the sampling distribution under $\hat F_n$, and it inherits every flaw of $\hat F_n$: no tail beyond the data (the maximum), lumpy distributions at small $n$ (Example 4.4), destroyed dependence. Before trusting it on a new statistic, run a coverage simulation like Part 3 of the project.
- **p-values are not probabilities of hypotheses** and never will be. If a sentence contains "probability" and "hypothesis" and no prior, it is wrong.
- **Multiple testing is where real research goes wrong**, and almost never because someone forgot Bonferroni. It goes wrong through unreported analyses: subgroups tried, outcomes switched, stopping after peeking. Run the simulations in Parts 8–9 and exercise C3.
- **Power analysis gets skipped** in practice, and with the winner's curse that is why literatures fill with exaggerated and false-positive effects. If you cannot state the effect size your study can detect with 80% power, you are not ready to collect data.

## 9. Check yourself

1. State, correctly, what "95%" means for a confidence interval. What is random?
2. Why does $(n-1)S^2/\sigma^2$ have $n-1$ degrees of freedom and not $n$? Answer in terms of a subspace.
3. Why are $\bar X$ and $S^2$ independent for Normal data? Name the step that fails for Exponential data.
4. What is a pivot? Give the pivot behind the $t$ interval and the one behind the variance interval.
5. Why does the Wald interval for a proportion fail at $X = 0$, and how does Wilson avoid it?
6. State the duality between tests and confidence sets in one line.
7. Write the percentile and basic bootstrap intervals. Why does the percentile interval for $\max_iX_i$ have coverage 0?
8. Prove the Neyman–Pearson lemma in four lines.
9. Write $-2\log\lambda$, Wald and score statistics for a one-parameter model. Which is not invariant to reparametrization?
10. Show that $p\sim\mathrm{Unif}(0,1)$ under a simple null with continuous $T$.
11. What $n$ detects a $0.5\sigma$ shift with power 0.8 at two-sided α = 0.05? Why does power 0.14 lead to exaggerated published effects?
12. Which procedure controls FWER under arbitrary dependence, and which controls FDR under independence? Where in the BH proof is independence used?

## 10. Pair with

- **Wasserman, *All of Statistics*, ch. 8 (bootstrap) and ch. 10 (testing, p-values).** Terse and complete: Wald tests, LRT, permutation tests and BH in a few pages each. A good second pass after §3.6.
- **MIT 18.650 (Rigollet), lecture notes ch. 4–5 and the OCW lectures.** Parametric tests and goodness of fit with a clean asymptotic emphasis.
- **Rice, *Mathematical Statistics and Data Analysis*, ch. 9 and 11.** More examples, including the rank-based two-sample tests skipped here.
- **Efron and Tibshirani, *An Introduction to the Bootstrap*, ch. 1–6 and 12–14.** The bootstrap from its inventor, through bootstrap-$t$ and BCa intervals.
- **The ASA statement on p-values (Wasserstein and Lazar, 2016) and Greenland et al. (2016), "Statistical tests, P values, confidence intervals, and power: a guide to misinterpretations."** Required reading; short, and better than any textbook at cataloguing the errors of §3.8.
- **Casella and Berger, *Statistical Inference*, ch. 8–9.** The proofs with full regularity conditions: Neyman–Pearson and Karlin–Rubin (§8.3), pivots and test inversion (§9.2).
