# ST-4: Generalized Linear Models: Exponential Families, Link Functions, IRLS, Deviance

**Time:** ~12 hours · **Prerequisites:** ST-3 (and ST-1 §3.10, exponential families) · **Unlocks:** ST-5, ST-9, ML-3

> **Conventions in this module.** Responses $y_1,\dots,y_n$ are independent with means $\mu_i = \mathbb{E}\,y_i$. $X$ is $n\times p$ with rows $x_i^\top$ (first column all ones), $\beta\in\mathbb{R}^p$. The **linear predictor** is $\eta_i = x_i^\top\beta + o_i$, where $o_i$ is a known **offset** (zero unless stated), and the **link** $g$ satisfies $g(\mu_i) = \eta_i$. The **exponential dispersion family** is $f(y;\theta,\phi) = \exp\{[y\theta - b(\theta)]/a_i(\phi) + c(y,\phi)\}$ with $a_i(\phi) = \phi/\omega_i$, where $\omega_i$ is a known **prior weight** ($\omega_i = m_i$ for a binomial proportion $y_i = Y_i/m_i$, otherwise 1). ST-1 called the natural parameter $\eta$; here it is $\theta$, because $\eta$ now means the linear predictor. $V(\mu) = b''(\theta(\mu))$ is the **variance function**. $w_i = \omega_i/[V(\mu_i)\,g'(\mu_i)^2]$ are the **working weights**, $W = \mathrm{diag}(w_i)$, and $z$ is the **working response** (3.10). The **unit deviance** is $d(y,\mu)$, the **deviance** $D = \sum_i\omega_i d(y_i,\hat\mu_i)$, so the *scaled* deviance is $D/\phi = 2[\ell(\text{saturated}) - \ell(\text{model})]$. $X^2 = \sum_i \omega_i(y_i - \hat\mu_i)^2/V(\hat\mu_i)$ is Pearson's statistic. $\sigma(t) = 1/(1+e^{-t})$ is the logistic function; $\Phi,\varphi$ are the standard Normal CDF and density. The Gamma is parametrized by its **mean** $\mu$ and **shape** $\nu = 1/\phi$ (so it is not ST-1's rate parametrization). Logs are natural. IRLS = iteratively reweighted least squares.

## 0. Where this sits

You need three things from earlier. From ST-1 §3.10: an exponential family has $A' = \mathbb{E}\,T$, $A'' = \mathrm{Cov}\,T$ (3.22), and its log-likelihood is concave with observed information equal to expected information (Theorem 3.26). From ST-2 §3.6: Wilks' theorem and the Wald/score/LRT trinity. From ST-3: least squares as projection, the hat matrix $H = X(X^\top X)^{-1}X^\top$, and weighted least squares. LA-14 supplies the chain rule and Newton's method (3.23).

This module takes linear regression and removes three restrictions: the response need not be Normal, its variance need not be constant, and its mean need not be linear in the predictors. What replaces them is one family of distributions, one link function, and one algorithm. ST-5 puts priors on the same likelihoods; ST-9 asks what the coefficients mean causally; ML-3 re-derives logistic and softmax regression as a learning problem; ML-7's sigmoid output unit is a GLM; SM-7's softmax is SM-4's Boltzmann distribution.

## 1. The problem that forced this idea

Linear regression assumes $y_i = x_i^\top\beta + \varepsilon_i$ with Normal errors of constant variance. Three common kinds of response violate all of it at once. The project's Part 1 shows each failure on simulated data.

**A binary outcome** ($y\in\{0,1\}$, with true $p(x) = \sigma(1.2x)$ for $x\in(-3,3)$, $n = 400$). The least-squares line is $0.493 + 0.197x$, and 16.0% of its fitted values fall outside $[0,1]$: it predicts probabilities of $-0.1$ and $1.1$. The residual variance is 0.031 where $|x|>2$ and 0.235 where $|x|<0.6$, because $\mathrm{Var}(y\mid x) = p(1-p)$ is a function of the mean; constant variance is false by construction. For $x>1.5$ the residual can only take two values ($-\hat p$ or $1-\hat p$); its skewness is $-2.26$ and a Shapiro test rejects Normality with $p\approx 10^{-14}$.

**A count** (Poisson with $\mu = e^{0.3+0.6x}$). 14.5% of the fitted means are negative, and the residual variance is 0.58 for $x<-1$ against 6.79 for $x>1$: the spread tracks the mean.

**A positive, skewed quantity** (waiting times, repair costs, rainfall). A Normal model puts mass below zero, and the standard deviation usually scales with the mean.

The pre-1970 remedies were transformations: arcsine-square-root for proportions, $\log(y+1)$ for counts. Each flattens one symptom and quietly changes the question. $\mathbb{E}[\log y] \ne \log\mathbb{E}[y]$ (Jensen, F5 Theorem 3.32), so you no longer model the mean; and $\log 0$ has no answer.

The way out came from several directions. Fisher's appendix to Bliss (1935) fitted a probit dose-response model by maximum likelihood with iterative weighting. Berkson (1944) proposed the logistic function as a simpler alternative (the logit). Nelder and Wedderburn (1972) saw that logistic, Poisson, Gamma and Normal regression are one model, with one fitting algorithm, and that the algorithm is weighted least squares repeated. The GLM is a unification, not an invention. The shared structure is ST-1's exponential family.

## 2. Intuition

**The link is a change of coordinates.** A line $\eta = x^\top\beta$ lives on all of $\mathbb{R}$. A probability lives in $(0,1)$ and a rate in $(0,\infty)$. Choose a monotone map $g$ from the mean's natural range onto $\mathbb{R}$ and put the line there: $g(\mu) = x^\top\beta$. For a probability the logit stretches $(0,1)$ onto $\mathbb{R}$, and the straight line in $\eta$-space becomes an S-curve in $\mu$-space that cannot leave $(0,1)$.

**The noise scale follows the mean.** Poisson counts near 100 scatter by about $\pm10$; counts near 1 scatter by about $\pm1$. Ordinary least squares weights all points equally, but the likelihood weights each by the reciprocal of its variance, which is $1/\mu_i$ on the count scale. The weights depend on the unknown mean, so you have a chicken-and-egg problem. The standard escape is to iterate: fit, recompute the weights from the fit, fit again.

**Linearize the link.** Around the current fit $\mu$, $g(y) \approx g(\mu) + (y-\mu)g'(\mu) = \eta + (y-\mu)g'(\mu)$. This **working response** is the observation, re-expressed on the $\eta$ scale. For a Poisson count with the log link it is $\eta + (y-\mu)/\mu$, which is $\log y$ corrected to first order about $\log\mu$ (and finite when $y = 0$). Its variance is $g'^2\cdot\mathrm{Var}(y)$, so its natural weight is the reciprocal (for Poisson with the log link, variance $(1/\mu)^2\mu = 1/\mu$ and weight $\mu$).

**A one-parameter experiment.** Poisson data $y = (2,4,6)$, intercept only, so every fitted mean is $e^\beta$ and the weights are all equal. Start at $\beta = 0$ and repeat "weighted-average the working responses":

$$\beta:\ 0 \to 3 \to 2.1991 \to 1.6427 \to 1.4165 \to 1.3867 \to 1.3863 = \log 4 = \log\bar y.$$

Far from the answer each pass removes about one unit; near it the error is squared ($4.5\times10^{-4}$, then $10^{-7}$). That is the signature of Newton's method (LA-14 (3.23)). At the end $\hat\mu = \bar y = 4$, so the fitted counts sum to the observed counts. **Guess the theorem:** each IRLS pass is a Newton step in disguise; the fixed point matches the data's sums; and with the right link the algorithm is safe. Section 3 turns the guesses into statements.

**The deviance generalizes the residual sum of squares.** In least squares, $(y-\mu)^2$ measures a gap in the units of a constant variance. When the variance changes with the mean, the right gap is $\approx (y-\mu)^2/V(\mu)$, and the exact version is a log-likelihood ratio against a model that fits each point perfectly (Section 3.6).

## 3. Formal development

### 3.1 The exponential dispersion family

**Definition 3.1.** $y$ has an **exponential dispersion family** distribution if its density (or PMF) is
$$f(y;\theta,\phi) = \exp\Big\{\frac{y\theta - b(\theta)}{a(\phi)} + c(y,\phi)\Big\},\qquad a(\phi) = \phi/\omega. \tag{3.1}$$
$\theta$ is the **natural parameter**, $b$ the **cumulant function**, $\phi>0$ the **dispersion**, $\omega$ a known weight. For fixed $\phi$ this is ST-1's family (3.21) with $T = y$, natural parameter $\theta/a$ and log-partition $b/a$.

**Theorem 3.1 (Mean and variance).** If $b$ is twice differentiable and the usual regularity (differentiation under the integral) holds,
$$\mathbb{E}\,y = \mu = b'(\theta),\qquad \mathrm{Var}(y) = a(\phi)\,b''(\theta) = a(\phi)\,V(\mu), \tag{3.2}$$
where, with $\theta(\mu) = (b')^{-1}(\mu)$,
$$V(\mu) := b''(\theta(\mu)).\tag{3.3}$$

▶ **Try it first.** Take the log of (3.1), differentiate once in $\theta$, and apply ST-1's two identities: the score has mean zero, and its variance equals minus the expected second derivative.

**Derivation.**
1. $\ell(\theta) = [y\theta - b(\theta)]/a + c$. The score is $\partial\ell/\partial\theta = [y - b'(\theta)]/a$, and $\partial^2\ell/\partial\theta^2 = -b''(\theta)/a$, which contains no $y$.
2. ST-1: $\mathbb{E}[\text{score}] = 0$. So $\mathbb{E}\,y = b'(\theta)$.
3. ST-1: $\mathbb{E}[\text{score}^2] = -\mathbb{E}[\partial^2\ell/\partial\theta^2]$. The left side is $\mathrm{Var}(y)/a^2$ (by step 2), the right side is $b''/a$. So $\mathrm{Var}(y) = a\,b''$.
4. $b'' = \mathrm{Var}(y)/a > 0$, so $b'$ is strictly increasing and invertible: $\theta$ is a function of $\mu$, and $V(\mu) = b''(\theta(\mu))$ is well defined. ∎

■ **What just happened.** One function $b$ holds the mean and the variance. The *mean* is its slope, and the *variance function* is its curvature rewritten as a function of the mean. This is ST-1's $A' = \mathbb{E}T$, $A'' = \mathrm{Var}\,T$ with a dispersion factor. Specifying $V(\mu)$ essentially specifies the family.

**The table.** Write each log-density in the form (3.1) and read off $\theta$, $b$ and $a$.

- **Normal** $\mathcal{N}(\mu,\sigma^2)$: $\log f = [y\mu - \mu^2/2]/\sigma^2 - y^2/(2\sigma^2) - \tfrac12\log(2\pi\sigma^2)$. So $\theta = \mu$, $b = \theta^2/2$, $\phi = \sigma^2$. Then $b' = \theta = \mu$ and $b'' = 1$.
- **Bernoulli**$(p)$: $\log f = y\log\frac{p}{1-p} + \log(1-p)$. So $\theta = \log\frac{p}{1-p}$ and, since $1-p = 1/(1+e^\theta)$, $b = -\log(1-p) = \log(1+e^\theta)$, $\phi = 1$. Then $b' = e^\theta/(1+e^\theta) = p$ and $b'' = p(1-p)$. For a **binomial proportion** $y = Y/m$, $\log f = m[y\theta - b(\theta)] + \log\binom{m}{my}$: same $\theta, b$ with $a = 1/m$.
- **Poisson**$(\mu)$: $\log f = y\log\mu - \mu - \log y!$. So $\theta = \log\mu$, $b = e^\theta$, $\phi = 1$, $b' = b'' = e^\theta = \mu$.
- **Gamma** (mean $\mu$, shape $\nu$): $\log f = \nu(-y/\mu - \log\mu) + \nu\log\nu + (\nu-1)\log y - \log\Gamma(\nu)$. So $\theta = -1/\mu$ (negative), $b = -\log(-\theta) = \log\mu$, $\phi = 1/\nu$. Then $b' = -1/\theta = \mu$ and $b'' = 1/\theta^2 = \mu^2$.
- **Inverse Gaussian** $(\mu,\lambda)$: $\log f = -\lambda(y-\mu)^2/(2\mu^2y) + \tfrac12\log\frac{\lambda}{2\pi y^3}$. Expanding the square, $\theta = -1/(2\mu^2)$, $b = -\sqrt{-2\theta} = -1/\mu$, $\phi = 1/\lambda$. Then $b' = (-2\theta)^{-1/2} = \mu$ and $b'' = (-2\theta)^{-3/2} = \mu^3$.

| Family | $\theta$ | $b(\theta)$ | $\phi$ | $\mu = b'(\theta)$ | $V(\mu)$ | Canonical link |
|---|---|---|---|---|---|---|
| Normal | $\mu$ | $\theta^2/2$ | $\sigma^2$ | $\theta$ | $1$ | identity |
| Bernoulli / Binomial | $\log\frac{\mu}{1-\mu}$ | $\log(1+e^\theta)$ | $1$ ($a=1/m$) | $\frac{e^\theta}{1+e^\theta}$ | $\mu(1-\mu)$ | logit |
| Poisson | $\log\mu$ | $e^\theta$ | $1$ | $e^\theta$ | $\mu$ | log |
| Gamma | $-1/\mu$ | $-\log(-\theta)$ | $1/\nu$ | $-1/\theta$ | $\mu^2$ | inverse ($-1/\mu$) |
| Inverse Gaussian | $-1/(2\mu^2)$ | $-\sqrt{-2\theta}$ | $1/\lambda$ | $(-2\theta)^{-1/2}$ | $\mu^3$ | $-1/(2\mu^2)$ |

The variance functions $1,\mu,\mu^2,\mu^3$ are the Tweedie family's integer powers (Jørgensen 1987; stated, not derived). Project Part 2 has SymPy re-derive every row of the table, check (3.2) against `scipy.stats` for each distribution, and recover the unit deviances of Section 3.6 by symbolic integration. Dispersion matters differently by family: Bernoulli and Poisson have $\phi = 1$ built in (the variance is *determined* by the mean, which is the root of overdispersion in Section 3.8); Normal and Gamma have a free $\phi$.

### 3.2 The three components, links, and the canonical link

**Definition 3.2 (GLM).** A generalized linear model has
1. a **random component**: $y_i$ independent, each from (3.1) with mean $\mu_i$ and a common $\phi$;
2. a **systematic component**: the linear predictor $\eta_i = x_i^\top\beta + o_i$;
3. a **link**: a strictly monotone, differentiable $g$ with
$$g(\mu_i) = \eta_i. \tag{3.4}$$

The **canonical link** is $g = (b')^{-1}$, i.e. $\eta_i = \theta_i$. Equivalently, since $\theta(\mu) = g(\mu)$, $g'(\mu) = d\theta/d\mu = 1/V(\mu)$ (by step 4 of Theorem 3.1). The last column of the table lists them: logit for Bernoulli, log for Poisson, identity for Normal, $-1/\mu$ for Gamma. (Software writes $1/\mu$; the sign only flips $\beta$.)

The link's job is to map the mean's range onto $\mathbb{R}$. Beyond the canonical choices:

- **Probit** $g = \Phi^{-1}$ and **complementary log-log** $g(\mu) = \log(-\log(1-\mu))$ for probabilities. Both arise from a latent-variable model with a different error distribution (Section 3.7).
- **Identity** for a Poisson mean (additive rates). It is legitimate only while $x^\top\beta > 0$, which the fit does not enforce.
- **Log** for a Gamma mean (multiplicative effects). The canonical $-1/\mu$ can produce negative means, so the log link is the practical default.

"Canonical" means algebraically convenient, not true. The next two sections show what the convenience buys.

### 3.3 The likelihood equations

**Theorem 3.2 (Score and information).** For a GLM with log-likelihood $\ell(\beta) = \sum_i\ell_i$, $\ell_i = \omega_i[y_i\theta_i - b(\theta_i)]/\phi + c(y_i,\phi)$,
$$\frac{\partial\ell}{\partial\beta_j} = \frac1\phi\sum_{i=1}^n\frac{\omega_i\,(y_i-\mu_i)\,x_{ij}}{V(\mu_i)\,g'(\mu_i)}, \tag{3.5}$$
$$I(\beta) = \mathbb{E}\Big[\frac{\partial\ell}{\partial\beta}\frac{\partial\ell}{\partial\beta^\top}\Big] = \frac1\phi X^\top WX,\qquad w_i = \frac{\omega_i}{V(\mu_i)\,g'(\mu_i)^2}. \tag{3.6}$$
For the canonical link the score collapses to
$$\frac{\partial\ell}{\partial\beta} = \frac1\phi X^\top\Omega(y-\mu),\qquad \Omega = \mathrm{diag}(\omega_i), \tag{3.7}$$
and the Hessian $-\partial^2\ell/\partial\beta\partial\beta^\top = X^\top WX/\phi$ contains no $y$: observed information equals expected information.

▶ **Try it first.** $\beta$ reaches $\ell_i$ through a chain $\beta\to\eta_i\to\mu_i\to\theta_i\to\ell_i$. Write $\partial\ell_i/\partial\beta_j$ as a product of four derivatives, one per arrow, and evaluate each. (LA-14's chain rule; this is the place everyone drops a factor.)

**Derivation.**
1. $\partial\ell_i/\partial\theta_i = \omega_i[y_i - b'(\theta_i)]/\phi = \omega_i(y_i-\mu_i)/\phi$.
2. $\mu_i = b'(\theta_i)$, so $d\mu_i/d\theta_i = b''(\theta_i) = V(\mu_i)$, and $d\theta_i/d\mu_i = 1/V(\mu_i)$.
3. $g(\mu_i) = \eta_i$ differentiates (inverse function theorem) to $d\mu_i/d\eta_i = 1/g'(\mu_i)$. **This is the factor that gets dropped.**
4. $\eta_i = x_i^\top\beta + o_i$, so $\partial\eta_i/\partial\beta_j = x_{ij}$.
5. Multiply: $\partial\ell_i/\partial\beta_j = \dfrac{\omega_i(y_i-\mu_i)}{\phi}\cdot\dfrac{1}{V(\mu_i)}\cdot\dfrac{1}{g'(\mu_i)}\cdot x_{ij}$. Summing over $i$ gives (3.5).
6. *Information.* The terms are independent with mean zero, so cross terms vanish in $\mathbb{E}[uu^\top]$ ($u$ = score) and $I_{jk} = \sum_i x_{ij}x_{ik}\,\omega_i^2\,\mathrm{Var}(y_i)/[\phi V g']^2$. With $\mathrm{Var}(y_i) = \phi V(\mu_i)/\omega_i$ from (3.2), each term is $\omega_i x_{ij}x_{ik}/[\phi\,V g'^2]$, i.e. (3.6).
7. *Canonical link.* $g' = 1/V$ gives $Vg' = 1$ in (3.5), hence (3.7). Differentiate (3.7): $-\partial^2\ell/\partial\beta_j\partial\beta_k = \frac1\phi\sum_i\omega_ix_{ij}x_{ik}\,d\mu_i/d\eta_i$. For the canonical link $\eta = \theta$, so $d\mu/d\eta = d\mu/d\theta = V$, and $w_i = \omega_i/(V g'^2) = \omega_i V$ too. The two agree, and no $y_i$ appears. ∎

■ **What just happened.** The score is "residual × predictor", each residual rescaled by the *variance* ($V$) and by the *link's stretch* ($g'$). Three consequences, all for the canonical link. (i) The likelihood equations are $X^\top\Omega y = X^\top\Omega\hat\mu$: **the fitted model reproduces the sufficient statistic $X^\top y$**, which is ST-1's moment matching (3.23). With an intercept, $\sum\omega_iy_i = \sum\omega_i\hat\mu_i$; for each indicator column, the group's fitted total equals its observed total. (ii) $\ell$ is concave: $X^\top WX = (W^{1/2}X)^\top(W^{1/2}X)\succeq 0$ (LA-8 Theorem 3.16), and positive definite when $X$ has full column rank and every $w_i>0$. A critical point is then the unique global maximum. This is ST-1 Theorem 3.26 with $\theta_i = x_i^\top\beta$ plugged in. (iii) Newton's method (negative Hessian) and Fisher scoring (expected negative Hessian) are the same algorithm. For non-canonical links none of this holds exactly: the score vanishes at $\hat\beta$, but it is not moment matching, and the Hessian has an extra term proportional to $(y-\mu)$ (3.11).

### 3.4 IRLS and Fisher scoring

**Theorem 3.3 (Fisher scoring is iteratively reweighted least squares).** The Fisher scoring step $\beta\leftarrow\beta + I(\beta)^{-1}u(\beta)$ equals
$$\beta^{\text{new}} = \beta + (X^\top WX)^{-1}X^\top WG(y-\mu) \tag{3.8}$$
$$\phantom{\beta^{\text{new}}}= (X^\top WX)^{-1}X^\top W z, \tag{3.9}$$
$$z_i = \eta_i - o_i + (y_i-\mu_i)\,g'(\mu_i),\qquad G = \mathrm{diag}(g'(\mu_i)), \tag{3.10}$$
with $W$, $\mu$, $\eta$ evaluated at the current $\beta$. (3.9) is the weighted least squares solution of regressing $z$ on $X$ with weights $w_i$.

▶ **Try it first.** Newton for maximizing $\ell$ is $\beta + (-\nabla^2\ell)^{-1}\nabla\ell$ (LA-14 (3.23)); Fisher scoring swaps $-\nabla^2\ell$ for $I$. Express $u$ through $W$ and $G$, then add and subtract $X^\top WX\beta$ to make the right side look like the normal equations.

**Derivation.**
1. By (3.5), $u = \frac1\phi X^\top WG(y-\mu)$, because $w_ig'_i = \omega_i/(V_ig'_i)$.
2. $I^{-1}u = \phi(X^\top WX)^{-1}\cdot\frac1\phi X^\top WG(y-\mu)$: **$\phi$ cancels**, which gives (3.8).
3. $\eta - o = X\beta$, so $X^\top Wz = X^\top WX\beta + X^\top WG(y-\mu)$.
4. Factor out $(X^\top WX)^{-1}$: $(X^\top WX)^{-1}X^\top Wz = \beta + (X^\top WX)^{-1}X^\top WG(y-\mu)$, which is (3.8). ∎

■ **What just happened.** Each Newton-type step is a weighted regression of a *constructed response* $z$ on $X$ (ST-3), so everything you know about weighted least squares applies to each pass. $z_i + o_i$ is the first-order expansion of $g(y_i)$ about $\mu_i$. Its variance is $g'^2\cdot\phi V/\omega_i = \phi/w_i$: the weights are exactly the inverse variances, which is when WLS is best (ST-3, Gauss–Markov with weights). $\phi$ never enters the iteration, so Normal and Gamma fits need no dispersion until you want standard errors.

| Model | $w_i$ | $z_i$ (with $o_i = 0$) |
|---|---|---|
| Normal, identity | $1$ | $y_i$ (one pass is OLS) |
| Bernoulli, logit | $p_i(1-p_i)$ | $\eta_i + (y_i-p_i)/[p_i(1-p_i)]$ |
| Poisson, log | $\mu_i$ | $\eta_i + (y_i-\mu_i)/\mu_i$ |
| Gamma, log | $1$ (Exercise B1) | $\eta_i + (y_i-\mu_i)/\mu_i$ |
| Bernoulli, probit | $\varphi(\eta_i)^2/[p_i(1-p_i)]$ | $\eta_i + (y_i-p_i)/\varphi(\eta_i)$ |

**The algorithm.** Initialize $\mu_i^{(0)}$ from the data, nudged off the boundary ($y_i + 0.5$ for counts, $(\omega_iy_i + 0.5)/(\omega_i+1)$ for proportions), set $\eta^{(0)} = g(\mu^{(0)})$. Repeat (3.9)–(3.10) until the relative change in deviance is below a tolerance (typically 4 to 8 passes). **Convergence.** For a canonical link this is Newton's method (LA-14 (3.23)) on a concave function: local convergence is quadratic (the rate is proved in ML-2; the trace below shows it), and there are no other local maxima to be trapped by. The project's logistic fit ($n = 300$) starts from $\beta = 0$ and its error $\|\beta_t - \hat\beta\|$ runs $1.1\times10^{-1},\ 4.7\times10^{-3},\ 8.5\times10^{-6},\ 2.8\times10^{-11},\ 5\times10^{-16}$: each error is about the square of the last. Undamped Newton can still overshoot from a poor start (Section 2's experiment); production code halves the step when the deviance rises.

**Non-canonical links.** Differentiate the $i$-th term of (3.5) with respect to $\eta_i$, using $d\mu/d\eta = 1/g'$ (write $V' = dV/d\mu$):
$$w_i^{\text{obs}} = w_i + \frac{\omega_i\,(y_i-\mu_i)\,[V'(\mu_i)g'(\mu_i) + V(\mu_i)g''(\mu_i)]}{V(\mu_i)^2\,g'(\mu_i)^3}. \tag{3.11}$$
The extra term has expectation zero, so $\mathbb{E}\,w^{\text{obs}} = w$. Newton uses $X^\top W^{\text{obs}}X$ on the left of (3.8), Fisher scoring uses $X^\top WX$. They reach the same $\hat\beta$ (the root of the same score) but at different rates. Your turn (c) fits an identity-link Poisson model ($n=100$): from $\beta=(1,1)$, Fisher's errors fall by a constant factor $\approx 0.04$ per pass ($1.8\times10^{-2},\ 6.6\times10^{-4},\ 2.5\times10^{-5}$), while Newton's go $1.1,\ 0.39,\ 4.8\times10^{-2},\ 6.7\times10^{-4},\ 1.1\times10^{-7}$: slow, then quadratic. For a good fit the extra term is small, and Fisher scoring is nearly as good.

**When the MLE does not exist.**

**Theorem 3.4 (Complete separation).** In logistic regression, suppose some $b\in\mathbb{R}^p$ satisfies $s_ix_i^\top b>0$ for all $i$, where $s_i = 2y_i - 1$ (a hyperplane classifies every point correctly). Then $\ell(cb)$ is strictly increasing in $c>0$, bounded above by $0$, and the likelihood has **no maximizer**.

*Proof.* With $a_i = s_ix_i^\top b>0$, $\ell(cb) = -\sum_i\log(1+e^{-ca_i})$ (since $\sigma(t)^y(1-\sigma(t))^{1-y} = \sigma(s\,t)$). Then
$$\frac{d}{dc}\ell(cb) = \sum_i\frac{a_ie^{-ca_i}}{1+e^{-ca_i}}>0\tag{3.12}$$
and $\ell(cb)\to0^-$ as $c\to\infty$. Since $\ell(\beta)<0$ for every finite $\beta$ (each $\sigma<1$) while $\ell(cb)\to0$, every $\beta$ is beaten by some multiple of $b$. ∎

(Quasi-complete separation, where $s_ix_i^\top b\ge0$ with some equalities, also destroys existence; Albert and Anderson (1984) give the full characterization.) **What IRLS does.** The working residual is $z_i - \eta_i = (y_i-\mu_i)/[\mu_i(1-\mu_i)]\to\pm1$ as $\mu_i\to y_i$, and the weights $\mu_i(1-\mu_i)$ concentrate on the points nearest the boundary. In the project's data ($x = 1,\dots,6$, $y = 000111$) the two points at $x=3$ and $x=4$ dominate, each with working residual $\pm1$, so every pass adds $2/(4-3) = 2$ to the slope. The measured slope is $1.03,\ 1.82,\ 2.88,\ 6.54$ (iteration 5), $16.6$ (10), $46.6$ (25), and the deviance falls to $10^{-3}$ at iteration 10 and $3\times10^{-10}$ at 25. The algorithm "converges" only because the deviance stops changing. Software that does not warn you hands you a slope of 46.6 with a standard error of $1.6\times10^5$, so the Wald statistic is $z = 0.0003$ for the strongest possible effect (the Hauck–Donner phenomenon, 1977).

**The ridge fix.** Maximize $\ell(\beta) - \frac\lambda2\|\beta_{\ne0}\|^2$ (intercept unpenalized). With $I_0 = \mathrm{diag}(0,1,\dots,1)$, the gradient is $u - \lambda I_0\beta$ and the negative Hessian is $X^\top WX + \lambda I_0$. The Newton step is $\beta + (X^\top WX+\lambda I_0)^{-1}(u-\lambda I_0\beta)$; add and subtract $\lambda I_0\beta$ exactly as in step 3 above:
$$\beta^{\text{new}} = (X^\top WX + \lambda I_0)^{-1}X^\top Wz. \tag{3.13}$$
For the separated data with $\lambda = 1$ the iteration converges in 5 passes to intercept $-3.922$, slope $1.121$ (boundary at $x = -\beta_0/\beta_1 = 3.922/1.121 = 3.5$, midway between $x=3$ and $x=4$). The penalty is a Normal prior on $\beta$ (ST-5: this is the MAP estimate); Firth's (1993) Jeffreys-prior penalty is the other standard repair.

### 3.5 Inference

ST-1's asymptotic normality of the MLE extends to independent, non-identically distributed GLM data under mild conditions on the design (Fahrmeir and Kaufmann 1985; stated, not proved here). With $I(\beta) = X^\top WX/\phi$ from (3.6),
$$\hat\beta\ \dot\sim\ \mathcal{N}\big(\beta,\ \phi\,(X^\top WX)^{-1}\big),\qquad \widehat{\mathrm{se}}(\hat\beta_j) = \sqrt{\phi\,[(X^\top WX)^{-1}]_{jj}}, \tag{3.14}$$
with $W$ evaluated at $\hat\beta$. The **Wald** statistic $z_j = \hat\beta_j/\widehat{\mathrm{se}}_j$ and interval $\hat\beta_j\pm1.96\,\widehat{\mathrm{se}}_j$ follow (ST-2 (3.16)). For Bernoulli and Poisson, $\phi = 1$ is an *assumption*. For Normal and Gamma, $\phi$ is estimated by the moment (Pearson) estimator, from $\mathbb{E}[\omega(y-\mu)^2/V(\mu)] = \phi$:
$$\hat\phi = \frac{X^2}{n-p}, \qquad X^2 = \sum_i\frac{\omega_i(y_i-\hat\mu_i)^2}{V(\hat\mu_i)}, \tag{3.15}$$
and then $z$ is replaced by $t_{n-p}$ (exact for the Normal: it is ST-3's $\hat\sigma^2 = \mathrm{RSS}/(n-p)$).

**Likelihood-ratio tests.** For nested models $M_0\subset M_1$ differing by $r$ parameters, $2[\ell(M_1)-\ell(M_0)] = (D_0-D_1)/\phi$, because both log-likelihoods contain the same saturated term (Section 3.6). Wilks (ST-2 Theorem 3.9) gives
$$\frac{D_0-D_1}{\phi}\ \dot\sim\ \chi^2_r,\qquad\text{or, with }\phi\text{ estimated,}\qquad F = \frac{(D_0-D_1)/r}{\hat\phi}\ \dot\sim\ F_{r,\,n-p_1}. \tag{3.16}$$
For the Normal this is exactly ST-3's extra-sum-of-squares $F$ test. An **analysis of deviance table** lists the sequential differences as terms enter (project Part 11: Poisson, $n = 150$; residual deviance $234.16\to182.44\to169.27\to168.25$ as $x_1,x_2,x_3$ enter; the changes $51.7,\ 13.2,\ 1.02$ have $\chi^2_1$ p-values $6.4\times10^{-13}$, $2.8\times10^{-4}$ and $0.31$; $x_3$ is pure noise). Under a true null the LRT statistic is $\chi^2_1$ to simulation accuracy: 2000 replications give a 5%-level rejection rate of 0.046 and a Kolmogorov–Smirnov $p = 0.83$. The Wald test is not invariant to reparametrization and fails near boundaries (the separation case); prefer the LRT in small samples (ST-2). **AIC** $= -2\ell(\hat\beta) + 2k$, where $k$ counts $\phi$ when it is estimated (ST-3).

### 3.6 Deviance and residuals

The **saturated model** has one free parameter per observation, so $\tilde\mu_i = y_i$ and $\tilde\theta_i = \theta(y_i)$. Writing $\hat\theta_i = \theta(\hat\mu_i)$,
$$\frac{D}{\phi} = 2[\ell(\text{saturated})-\ell(\hat\mu)] = \frac2\phi\sum_i\omega_i\big[y_i(\tilde\theta_i-\hat\theta_i) - b(\tilde\theta_i) + b(\hat\theta_i)\big]. \tag{3.17}$$
$c(y,\phi)$ cancels, so the deviance needs neither $\phi$ nor the normalizing constant. Define the **unit deviance**
$$d(y,\mu) = 2\big[y(\tilde\theta-\hat\theta) - b(\tilde\theta)+b(\hat\theta)\big],\qquad D = \sum_i\omega_i\,d(y_i,\hat\mu_i). \tag{3.18}$$

**Theorem 3.5 (Integral form).** $d(y,\mu) = 2\displaystyle\int_\mu^y\frac{y-t}{V(t)}\,dt$. Hence $d\ge0$ with equality iff $\mu=y$, and $d(y,\mu)\approx(y-\mu)^2/V(\mu)$ for $y$ near $\mu$. $\tag{3.19}$

▶ **Try it first.** Differentiate $d(y,\mu)$ with respect to $\mu$ only. What does $\partial\hat\theta/\partial\mu$ equal?

**Derivation.**
1. $\hat\theta = \theta(\mu)$ and $\hat\theta' = 1/V(\mu)$ (Theorem 3.1, step 4). $\partial d/\partial\mu = 2[-y\hat\theta' + b'(\hat\theta)\hat\theta'] = -2(y-\mu)\hat\theta' = -2(y-\mu)/V(\mu)$.
2. $d(y,y) = 0$. Integrate $\partial d/\partial t$ from $t=y$ down to $t=\mu$: $d(y,\mu) = \int_y^\mu\!-2(y-t)/V(t)\,dt = 2\int_\mu^y(y-t)/V(t)\,dt$.
3. The integrand $(y-t)/V(t)$ has the sign of $y-t$, so the integral is $\ge0$ for either order of $y,\mu$. Expanding $1/V(t)\approx1/V(\mu)$ gives $2\cdot(y-\mu)^2/(2V(\mu))$. ∎

■ **What just happened.** The deviance is the squared gap $(y-\mu)^2$ *integrated against the local precision* $1/V$. Where the variance is constant it is the squared gap itself; elsewhere it adapts, and to second order it equals Pearson's $X^2$ (3.15). Closed forms (each confirmed by quadrature in the project):
$$\begin{aligned}
\text{Normal: } & (y-\mu)^2, & \text{Poisson: } & 2\big[y\log\tfrac{y}{\mu} - (y-\mu)\big], \\
\text{Bernoulli: } & 2\big[y\log\tfrac y\mu + (1-y)\log\tfrac{1-y}{1-\mu}\big], & \text{Gamma: } & 2\big[-\log\tfrac y\mu + \tfrac{y-\mu}{\mu}\big],
\end{aligned}\tag{3.20}$$
and Inverse Gaussian $(y-\mu)^2/(y\mu^2)$ (Exercise B8). The Poisson one comes straight from (3.18): $\tilde\theta = \log y$, $\hat\theta = \log\mu$, $b = e^\theta$, so $d = 2[y(\log y - \log\mu) - y + \mu]$. With an intercept and the canonical link, (3.7) gives $\sum(y_i-\hat\mu_i)=0$, so $D = 2\sum y_i\log(y_i/\hat\mu_i)$. For the Normal, $D = \mathrm{RSS}$: the deviance is the generalized RSS, and the analysis of deviance generalizes ANOVA.

**Residuals.** Two definitions, both squaring to a goodness-of-fit statistic:
$$r_i^P = \frac{\sqrt{\omega_i}\,(y_i-\hat\mu_i)}{\sqrt{V(\hat\mu_i)}},\qquad r_i^D = \mathrm{sign}(y_i-\hat\mu_i)\sqrt{\omega_i\,d(y_i,\hat\mu_i)},\qquad \sum(r^P_i)^2 = X^2,\ \ \sum(r^D_i)^2 = D. \tag{3.21}$$
(The deviance residual is the signed square root of the unit deviance, by construction.) In practice, for small counts the deviance residuals look closer to Normal than the Pearson ones (stated, not derived). To account for leverage, divide by $\sqrt{\phi(1-h_{ii})}$.

**Hat matrix and influence.** At convergence $\hat\beta$ is the WLS fit of $z$ on $X$ with weights $W$. Put $\tilde X = W^{1/2}X$; then
$$H = W^{1/2}X(X^\top WX)^{-1}X^\top W^{1/2} = \tilde X(\tilde X^\top\tilde X)^{-1}\tilde X^\top \tag{3.22}$$
is the projection onto $C(\tilde X)$: symmetric, idempotent, $\mathrm{tr}\,H = p$, by ST-3's proofs applied to $\tilde X$. The leverage $h_{ii} = w_ix_i^\top(X^\top WX)^{-1}x_i$ now involves the weights, so a point can have extreme $x$ and *small* leverage (in a logistic fit, $\hat p_i\approx0$ or $1$ makes $w_i$ tiny). ST-3's leave-one-out argument on the final WLS fit gives a one-step Cook's distance (Pregibon 1981):
$$C_i\approx\frac{(r^P_i)^2\,h_{ii}}{p\,\phi\,(1-h_{ii})^2}. \tag{3.23}$$
Project Part 12 checks $\mathrm{tr}\,H = 3$, symmetry and idempotence, $\sum(r^D)^2 = D$ and $\sum(r^P)^2 = X^2$.

**Goodness of fit: where the deviance is $\chi^2$ and where it is not.** $D/\phi$ is a Wilks statistic against the saturated model, so $D/\phi\ \dot\sim\ \chi^2_{n-p}$ holds when the saturated model has a *fixed* dimension while the information per observation grows: grouped binomial data with $m_i\to\infty$, or Poisson counts with $\mu_i\to\infty$. In the project, 8 groups of $m=100$ give $\bar D = 5.99$ against $n-p = 6$, and 4.6% rejections at the 5% point. It fails for Bernoulli data:

**Theorem 3.6 (Bernoulli deviance is a function of the fitted values).** For logistic regression on 0/1 data, with canonical link and $\hat\eta_i = x_i^\top\hat\beta$,
$$D = -2\Big[\hat\beta^\top X^\top\hat p - \sum_i\log(1+e^{\hat\eta_i})\Big]. \tag{3.24}$$

*Proof.* (1) For $y\in\{0,1\}$ the saturated log-likelihood is 0, and $\ell(\hat\mu) = \sum[y_i\hat\eta_i - \log(1+e^{\hat\eta_i})]$. So $D = -2\sum[y_i\hat\eta_i - \log(1+e^{\hat\eta_i})]$. (2) $\sum y_i\hat\eta_i = \hat\beta^\top X^\top y$. (3) By (3.7), $X^\top y = X^\top\hat p$. Substitute. ∎

■ **What just happened.** Given the fitted probabilities, $D$ is a fixed number: the data touch it only through $X^\top y$, which the fit has already matched. It cannot measure misfit, and it is not $\chi^2_{n-p}$. For $n=100$ with $\beta_1=3$ and the model *exactly right*, the simulated mean of $D$ is 70.5 (sd 11.3), against $n-p = 98$ for a $\chi^2_{98}$ (sd 14.0). An intercept-only Bernoulli fit to 15 successes in 50 trials has $D = 61.09$ whatever the order of the data (Exercise B6). *Differences* of deviances between nested models remain valid (3.16): the saturated terms cancel. For Bernoulli data use grouping, the Hosmer–Lemeshow idea (3.28), or residual diagnostics.

### 3.7 Logistic regression

**Odds, odds ratios, marginal effects.** The model is $\log[p_i/(1-p_i)] = x_i^\top\beta$, i.e. odds $p_i/(1-p_i) = e^{x_i^\top\beta}$. Raising $x_{ij}$ by one unit multiplies the odds by $e^{\beta_j}$: $e^{\beta_j}$ is an **odds ratio**. The effect on the *probability* is not constant: since $\sigma' = \sigma(1-\sigma)$, $\partial p/\partial x_j = \beta_jp(1-p)$, so $|\partial p/\partial x_j|\le|\beta_j|/4$. The odds ratio is not the **relative risk** $p_1/p_0$. With baseline odds $o_0 = p_0/(1-p_0)$ and $p_1 = \mathrm{OR}\,o_0/(1+\mathrm{OR}\,o_0)$, using $1+o_0 = 1/(1-p_0)$,
$$\mathrm{RR} = \frac{p_1}{p_0} = \frac{\mathrm{OR}}{1-p_0+p_0\,\mathrm{OR}}. \tag{3.25}$$
OR $\approx$ RR only when $p_0$ is small (Example 4.3).

**Latent-variable derivation.** Let $Y^* = x^\top\beta + \varepsilon$ and $Y = \mathbb{1}\{Y^*>0\}$, with $\varepsilon$ having CDF $F$ symmetric about 0. Then
$$P(Y=1\mid x) = P(\varepsilon>-x^\top\beta) = 1 - F(-x^\top\beta) = F(x^\top\beta). \tag{3.26}$$
A logistic $\varepsilon$ ($F=\sigma$) gives the logit model; a standard Normal gives the probit. A Gumbel $F(t) = \exp(-e^{-t})$ is not symmetric: $1-F(-\eta) = 1-\exp(-e^\eta)$, which is the complementary log-log. The threshold fixes the location but the *scale* of $\varepsilon$ is not identified, only $\beta/\mathrm{sd}(\varepsilon)$. The logistic has sd $\pi/\sqrt3 = 1.814$, the Normal has 1. Which constant converts probit coefficients to logit ones? Matching CDF slopes at 0 gives $\sigma'(0)/\varphi(0)\ \Rightarrow\ c = 4\varphi(0) = 1.596$; matching sds gives 1.814; the least-squares fit of $\sigma(ct)$ to $\Phi(t)$ on the whole line gives $c = 1.701$. Simulating $n = 40{,}000$ from a latent-Normal model with slope 1, the probit fit returns 0.998, the logit fit 1.701, ratio 1.705. Use "logit $\approx 1.6$ to $1.7\times$ probit".

**Loss, physics, boundary.** $-\ell/n = -\frac1n\sum[y_i\log\hat p_i + (1-y_i)\log(1-\hat p_i)]$ is the **binary cross-entropy** of ML-3 (minimizing it minimizes the KL divergence from the empirical distribution to the model, ST-1 Theorem 3.28). A two-state system with energies $E_0,E_1$ has $P(1) = \sigma(-(E_1-E_0)/k_BT)$ (SM-4): the logit is minus an energy gap in units of $k_BT$, which is SM-7's sigmoid-as-Boltzmann. Classifying $\hat p\ge t$ is $x^\top\beta\ge\log\frac{t}{1-t}$: a **hyperplane** with normal $\beta$, whatever the data (the signed distance from the boundary is $x^\top\beta/\|\beta\|$, a projection onto $\beta$ as in LA-5). Curved boundaries need feature maps (ML-3).

**Inference.** Standard errors are $\sqrt{[(X^\top WX)^{-1}]_{jj}}$ with $w_i = \hat p_i(1-\hat p_i)$, then Wald and LRT as in Section 3.5. An interval for the odds ratio is $\exp$ of the interval for $\beta_j$ (this transformation is why one reports intervals, not standard errors, for ORs).

**Discrimination: ROC and AUC.** At threshold $t$ the **true positive rate** is $\mathrm{TPR}(t) = P(S>t\mid Y=1)$ and the **false positive rate** $\mathrm{FPR}(t) = P(S>t\mid Y=0)$ for a score $S$ (here $\hat p$). The **ROC curve** plots $\mathrm{TPR}$ against $\mathrm{FPR}$ as $t$ sweeps from $+\infty$ to $-\infty$; the **AUC** is the area beneath it.

**Theorem 3.7 (AUC is a probability).** Let $S_1,S_0$ be independent scores of a randomly chosen positive and negative case, with continuous distributions. Then
$$\mathrm{AUC} = P(S_1>S_0). \tag{3.27}$$
For a finite sample with $n_1$ positives and $n_0$ negatives, $\mathrm{AUC} = \frac1{n_1n_0}\sum_{i:y_i=1}\sum_{j:y_j=0}[\mathbb{1}\{s_i>s_j\}+\frac12\mathbb{1}\{s_i=s_j\}]$, the Mann–Whitney $U/(n_1n_0)$.

▶ **Try it first.** Write $\mathrm{TPR}$ and $\mathrm{FPR}$ with the CDFs $F_1,F_0$ of $S_1,S_0$, and integrate $\mathrm{TPR}\,d\mathrm{FPR}$.

**Derivation.**
1. $\mathrm{TPR}(t) = 1-F_1(t)$ and $\mathrm{FPR}(t) = 1-F_0(t)$, so $d\,\mathrm{FPR} = -f_0(t)\,dt$. As $t$ goes from $+\infty$ to $-\infty$ the curve runs from $(0,0)$ to $(1,1)$.
2. $\mathrm{AUC} = \int_{-\infty}^\infty[1-F_1(t)]\,f_0(t)\,dt = \int P(S_1>t)\,f_0(t)\,dt$.
3. This is $\mathbb{E}[P(S_1>S_0\mid S_0)] = P(S_1>S_0)$ by independence (Adam's law, F5 Theorem 3.26).
4. Empirical version: the ROC is a staircase; each (positive, negative) pair where the positive outranks the negative contributes a rectangle of area $1/(n_1n_0)$, and a tie contributes a diagonal trapezoid of half that. ∎

■ **What just happened.** AUC is the probability that a random positive outranks a random negative. It depends only on the *ranking* (it is unchanged by any increasing transformation of the score), so it says **nothing about calibration**. A model is **calibrated** if $P(Y=1\mid\hat p=q) = q$. Project Part 6: replacing $\hat p$ by $\hat p^3$ leaves the AUC at 0.8229 but wrecks the calibration plot. The **Hosmer–Lemeshow** idea (1980) sorts cases by $\hat p$, splits them into $g=10$ groups with observed positives $O_g$, expected positives $E_g=\sum_g\hat p_i$ and size $n_g$, and computes
$$C = \sum_{g}\frac{(O_g-E_g)^2}{E_g(1-E_g/n_g)}\ \approx\ \chi^2_{g-2}\ \text{(heuristic)}. \tag{3.28}$$
The fitted model gets $C = 9.60$ ($p = 0.29$); the same ranking through $\hat p^3$ gets $C = 2639$.

**Multinomial logistic regression.** Let $y\in\{0,1\}^K$ be a one-hot vector with probabilities $\pi_k$. Using $y_K = 1-\sum_{k<K}y_k$,
$$\prod_k\pi_k^{y_k} = \exp\Big\{\sum_{k<K}y_k\log\frac{\pi_k}{\pi_K} + \log\pi_K\Big\},$$
an exponential family with $\theta_k = \log(\pi_k/\pi_K)$ and $b(\theta) = -\log\pi_K = \log(1+\sum_{k<K}e^{\theta_k})$. Then $\partial b/\partial\theta_k = e^{\theta_k}/(1+\sum_je^{\theta_j}) = \pi_k$ and $\partial^2b/\partial\theta_k\partial\theta_l = \pi_k\delta_{kl}-\pi_k\pi_l$ (the multinomial covariance, Theorem 3.1 in vector form). Setting $\theta_K:=0$ gives the **softmax**
$$\pi_k = \frac{e^{\theta_k}}{\sum_{j=1}^Ke^{\theta_j}},\qquad \theta_k = x^\top\beta_k. \tag{3.29}$$
Only $K-1$ logit vectors are free: adding a common vector to every $\beta_k$ changes nothing (the project asserts this). The canonical-link score is again residual times predictor: $\partial(-\ell)/\partial\beta_k = \sum_i(s_{ik}-y_{ik})x_i$, which is LA-14 (3.21), and the Hessian $\sum_i(\mathrm{diag}\,s_i - s_is_i^\top)\otimes x_ix_i^\top$ is PSD, so the fit is a convex problem (ML-3). The decision regions $\{x: (\beta_k-\beta_j)^\top x\ge0\ \forall j\}$ are convex cells. In the project (3 Gaussian clusters, $n=600$), gradient descent with step 1.0 for 6000 steps reaches loss 0.36764, equal to L-BFGS's to five digits, gradient norm $10^{-16}$, accuracy 0.862. **Ordinal** responses use the same latent-variable story with cut points $\alpha_1<\dots<\alpha_{K-1}$: the proportional-odds model $\mathrm{logit}\,P(Y\le k\mid x) = \alpha_k - x^\top\beta$.

### 3.8 Poisson regression, overdispersion, and the Gamma GLM

**Poisson with a log link.** $\mu_i = \exp(x_i^\top\beta)$, so $e^{\beta_j}$ is a **rate ratio**: one unit of $x_j$ multiplies the expected count by $e^{\beta_j}$. **Offsets.** If $y_i$ counts events over an exposure $t_i$ (time, area, population) at rate $\lambda_i = e^{x_i^\top\beta}$ per unit exposure, then $y_i\sim\mathrm{Poisson}(t_i\lambda_i)$, because Poisson counts over disjoint pieces of exposure add (F5: a sum of independent Poissons is Poisson), so doubling $t_i$ doubles the mean. Hence
$$\log\mu_i = \log t_i + x_i^\top\beta, \tag{3.30}$$
an **offset** $o_i = \log t_i$ whose coefficient is fixed at 1; IRLS handles it by regressing $z-o$ (3.10). Entering $\log t_i$ as a free covariate *tests* the coefficient-1 assumption: on simulated sites with $t\in(1,10)$ the free coefficient is $1.127\pm0.079$ (interval covers 1); ignoring exposure altogether leaves deviance 651.7 against 399.9 with the offset. The Poisson deviance is (3.20).

**Overdispersion.** A Poisson variable has $\mathrm{Var}=\mu$, so $\phi=1$ is a theorem *inside* the model. Real counts vary more because rates differ between units for reasons not in $X$.

**Theorem 3.8 (Gamma–Poisson mixture is Negative Binomial).** Let $\lambda\sim\mathrm{Gamma}(\text{shape }r,\ \text{mean }\mu)$ and $y\mid\lambda\sim\mathrm{Poisson}(\lambda)$. Then
$$P(y=k) = \frac{\Gamma(k+r)}{k!\,\Gamma(r)}\Big(\frac{r}{r+\mu}\Big)^r\Big(\frac{\mu}{r+\mu}\Big)^k,\qquad \mathbb{E}\,y = \mu,\quad\mathrm{Var}(y) = \mu+\frac{\mu^2}{r}. \tag{3.31}$$

▶ **Try it first.** Integrate $\lambda$ out: multiply the two densities, collect the powers of $\lambda$ into a Gamma integral. For the variance use Adam's and Eve's laws (F5 Theorems 3.26, 3.28).

**Derivation.**
1. $P(y=k) = \int_0^\infty\frac{e^{-\lambda}\lambda^k}{k!}\cdot\frac{(r/\mu)^r\lambda^{r-1}e^{-r\lambda/\mu}}{\Gamma(r)}d\lambda = \frac{(r/\mu)^r}{k!\,\Gamma(r)}\int_0^\infty\lambda^{k+r-1}e^{-\lambda(\mu+r)/\mu}\,d\lambda$.
2. $\int_0^\infty\lambda^{a-1}e^{-c\lambda}d\lambda = \Gamma(a)/c^a$ (substitute $u=c\lambda$ in F4 (3.29); the same integral gives F5 §3.5's Gamma–Poisson posterior), with $a=k+r$ and $c=(\mu+r)/\mu$, gives $\Gamma(k+r)\big(\frac{\mu}{\mu+r}\big)^{k+r}$.
3. Multiply: $\big(\frac r\mu\big)^r\big(\frac\mu{\mu+r}\big)^r = \big(\frac r{r+\mu}\big)^r$, leaving the first claim of (3.31).
4. Adam: $\mathbb{E}\,y = \mathbb{E}[\mathbb{E}(y\mid\lambda)] = \mathbb{E}\lambda = \mu$. Eve: $\mathrm{Var}(y) = \mathbb{E}[\mathrm{Var}(y\mid\lambda)]+\mathrm{Var}(\mathbb{E}[y\mid\lambda]) = \mathbb{E}\lambda+\mathrm{Var}\lambda = \mu+\mu^2/r$, since a Gamma with shape $r$ and mean $\mu$ has variance $\mu^2/r$. ∎

■ **What just happened.** Heterogeneity in the rate adds a $\mu^2/r$ term to the variance: the variance grows *faster* than the mean. For fixed $r$ the NB is an exponential family: with $e^\theta = \mu/(\mu+r)$, so $1-e^\theta = r/(\mu+r)$, the PMF is $\exp\{k\theta - b(\theta)+c\}$ with $b = -r\log(1-e^\theta)$. Then $b' = re^\theta/(1-e^\theta) = \mu$ and $b'' = b'/(1-e^\theta) = \mu(1+\mu/r)$. So it is a GLM with $V(\mu) = \mu+\mu^2/r$; with the log link $w = \mu/(1+\mu/r)$. Fit $r$ by maximizing the profile log-likelihood (one IRLS per $r$).

**Detecting and handling it.** Compute $\hat\phi = X^2/(n-p)$ (3.15); far above 1 means overdispersion. Poisson standard errors are then too small by about $\sqrt{\hat\phi}$. **Quasi-likelihood** (Wedderburn 1974): the estimating equation (3.5) involves only $\mu$ and $V$. Its root $\hat\beta$ is consistent whenever the *mean* model is right, whatever the true variance, because $\mathbb{E}[\text{(3.5)}]=0$ term by term. So **quasi-Poisson** keeps the Poisson $\hat\beta$ and multiplies the standard errors by $\sqrt{\hat\phi}$ (no likelihood, hence no AIC; use the $F$ test in (3.16)). The **NB GLM** instead changes the weights: $\mu/(1+\mu/r)$ against Poisson's $\mu$ downweights the large counts (Exercise B5), which is more efficient when the variance really is quadratic. Project Part 10 (300 replications, $n=200$, true $r=2$, true slope 0.5):

| Method | 95% interval coverage for the slope |
|---|---|
| Poisson | 0.743 |
| quasi-Poisson | 0.927 |
| Negative Binomial ($r$ by profile likelihood) | 0.957 |

The mean $\hat\phi$ is 2.52 and the median $\hat r$ is 2.04. On one dataset the slope's standard error is 0.0460 (Poisson), 0.0633 (quasi: $\times\sqrt{\hat\phi} = 1.375$), 0.0666 (NB, $\hat r=3.00$), and AIC is 871.2 for the Poisson fit against 823.1 for the NB (counting $r$). **Zero inflation** is the other failure: structural zeros (a machine that was off) give $P(0) = \pi+(1-\pi)e^{-\mu}$ and $P(k) = (1-\pi)\,\mathrm{Pois}(k;\mu)$ for $k>0$. This is a mixture, not an exponential family, so you maximize it directly or by EM (Your turn (b): $\hat\pi = 0.307$ for truth 0.3; the data have 39.9% zeros where the Poisson GLM expects 22.3%).

**The Gamma GLM.** For positive skewed $y$, $V=\mu^2$ gives $\mathrm{sd} = \sqrt\phi\,\mu$: a **constant coefficient of variation** $1/\sqrt\nu$. With the log link, effects are multiplicative and $\mathbb{E}\,y = e^{x^\top\beta}$. The alternative of OLS on $\log y$ models $\mathbb{E}\log y$, not $\log\mathbb{E}\,y$: if $y = \mu G$ with $G\sim\mathrm{Gamma}(\nu,\text{mean }1)$, then $\mathbb{E}\log y = \log\mu+\psi(\nu)-\log\nu$ ($\psi$ the digamma function). Slopes are unaffected; the intercept shifts by $\psi(\nu)-\log\nu$, which is $-0.176$ at $\nu=3$, so back-transformed OLS predicts means 16% too low (Jensen). Your turn (a): GLM $\hat\beta = (0.978, 0.521)$ against OLS-on-log $(0.793, 0.517)$ for truth $(1, 0.5)$. With the log link $W=I$ (Exercise B1): IRLS is plain least squares on $z$.

### 3.9 GLMs as maximum-entropy models

ST-1 Theorem 3.27: among distributions with $\mathbb{E}[T]=\tau$, the entropy maximizer is the exponential-family member with that mean. The regression version (take $\omega_i=\phi=1$): choose each $p_i(y)$, relative to the base measure $h$, to maximize $-\sum_i\int p_i\log(p_i/h)$ subject to **one constraint per column of $X$**: $\sum_ix_i\mathbb{E}_{p_i}[y] = \sum_ix_iy_i^{\rm obs}$, i.e. the model must reproduce the observed $X^\top y$.

1. Lagrangian with multipliers $\beta\in\mathbb{R}^p$ and $\lambda_i$ (normalization). Stationarity in $p_i(y)$: $-\log\frac{p_i(y)}{h(y)}-1+\lambda_i+y\,x_i^\top\beta=0$.
2. So $p_i(y)\propto h(y)\exp(y\,x_i^\top\beta)$: an exponential family with $\theta_i = x_i^\top\beta$, the canonical-link GLM, with the multipliers as its coefficients.
3. The constraints then read $\sum_ix_i\mu_i(\beta) = X^\top y$, which is (3.7). That this stationary point is the maximum is ST-1's Gibbs argument (Theorem 3.27).

**A canonical-link GLM is the least-committal model reproducing the data's sufficient statistic $X^\top y$.** The score equation is the constraint, and maximum likelihood is the dual problem: minimizing $\sum_ib(x_i^\top\beta)-\beta^\top X^\top y = -\ell$ over $\beta$. This is the bridge to SM-4 (the canonical ensemble maximizes entropy at fixed $\langle E\rangle$) and SM-7. Non-canonical links have no such reading, which is the best argument for treating canonical as special.

## 4. Worked examples

All numbers below are asserted or printed by `project.py`.

### Example 4.1: The table, for Bernoulli and Poisson

*Bernoulli, $p=0.3$.* $\log f = y\log\frac p{1-p}+\log(1-p)$, so $\theta = \log(0.3/0.7) = -0.8473$. Then $1-p = 1/(1+e^\theta)$ gives $b(\theta) = \log(1+e^\theta) = \log(10/7) = 0.3567 = -\log0.7$. Check: $b'(\theta) = e^\theta/(1+e^\theta) = (3/7)/(10/7) = 0.3$ and $b''(\theta) = 0.3\times0.7 = 0.21 = V(\mu)$. As a binomial proportion with $m=10$, $a = 1/10$ and $\mathrm{Var}(\bar y) = 0.021$.
*Poisson, $\mu=3.7$.* $\theta = \log3.7 = 1.3083$, $b = e^\theta = 3.7$, $b' = b'' = 3.7$: mean and variance are both 3.7, and $V(\mu)=\mu$.

### Example 4.2: Six data points, one IRLS pass by hand

Data $x = (-2,-1,0,1,2,3)$, $y = (0,0,1,0,1,1)$, $X = [\mathbf 1\ x]$, logit link.
*Start $\beta^{(0)}=0$.* Every $p_i = \tfrac12$, so $w_i = \tfrac14$ and $z_i = 0+(y_i-\tfrac12)/\tfrac14 = 4y_i-2 = (-2,-2,2,-2,2,2)$. The score (3.7) is $X^\top(y-\tfrac12) = (0,\ 3.5)$ (since $\sum x_iy_i = 5$ and $\sum x_i = 3$). $X^\top X = \begin{pmatrix}6&3\\3&19\end{pmatrix}$, $\det = 105$, $X^\top z = (0,14)$.
*Pass 1 via (3.9).* The $\tfrac14$ in $W$ cancels: $\beta^{(1)} = (X^\top X)^{-1}X^\top z = \frac1{105}\begin{pmatrix}19&-3\\-3&6\end{pmatrix}\binom{0}{14} = (-0.4,\ 0.8)$. *Via (3.8):* $0 + 4(X^\top X)^{-1}(0,3.5)^\top$, the same vector.
*Pass 2.* $p_i = \sigma(-0.4+0.8x_i)$, new weights and working responses: $\beta^{(2)} = (-0.5549,\ 1.1098)$.
*Convergence.* Six passes reach $\hat\beta = (-0.6070,\ 1.2140)$ with $\widehat{\mathrm{se}} = (1.2181,\ 0.9126)$ from (3.14). Residual deviance 4.9560; null deviance $12\log2 = 8.3178$ (the intercept-only fit has $\hat p=\tfrac12$). Check: $X^\top(y-\hat p) = 0$ to $10^{-8}$. *Tests for the slope.* LRT: $8.3178-4.9560 = 3.3618$, $p = 0.067$. Wald: $z = 1.2140/0.9126 = 1.330$, $p = 0.18$. At $n=6$ the two disagree (ST-2's trinity), and the LRT is the one to trust.

### Example 4.3: An odds ratio, and why it is not a risk ratio

In Part 5's synthetic crack-detection data ($n=200$; $x_1$ a standardized stress amplitude), $\hat\beta_1 = 1.2057$ with Wald interval $[0.7871, 1.6243]$. One standard deviation more stress multiplies the odds by $e^{1.2057} = 3.34$, interval $[e^{0.7871}, e^{1.6243}] = [2.20,\ 5.08]$. What that does to a *probability* depends on where you start. Take $\beta=0.7$, $\mathrm{OR}=e^{0.7}=2.0138$, and (3.25):

| baseline $p_0$ | new $p_1$ | RR |
|---|---|---|
| 0.01 | 0.0199 | 1.9935 |
| 0.10 | 0.1828 | 1.8284 |
| 0.50 | 0.6682 | 1.3364 |

Doubling the odds nearly doubles a rare risk and does nothing like that to a common one. Reading "OR 2.0" as "twice as likely" is the standard misreading.

### Example 4.4: Complete separation

$x=1,\dots,6$, $y=000111$ is separated by any threshold in $(3,4)$. IRLS from $\beta=0$ gives (intercept, slope, deviance): iteration 1: $(-3.600,\ 1.029,\ 2.94)$; 2: $(-6.360,\ 1.817,\ 1.65)$; 5: $(-22.89,\ 6.539,\ 0.149)$; 10: $(-58.22,\ 16.64,\ 9.8\times10^{-4})$; 25: $(-163.2,\ 46.64,\ 3.0\times10^{-10})$. From iteration 9 on each pass adds $2.00$ to the slope (to within 0.05), as predicted. The ridge fit ($\lambda=1$) is $(-3.922,\ 1.121)$ and its first pass from 0 is $(-2.930,\ 0.837)$, from (3.13) with $\mathrm{diag}(0,1)$.

### Example 4.5: Poisson deviance with and without an offset

A Geiger counter records counts $y = (3,5,9,21)$ over dwell times $t = (1,2,4,8)$ s.
*With the offset* $\log t$ and an intercept: (3.7) gives $\sum\hat\mu_i = \sum y_i = 38$, and $\hat\mu_i = t_i\hat\lambda$, so $\hat\lambda = 38/15 = 2.5333$ s$^{-1}$ and $\hat\mu = (2.533,\ 5.067,\ 10.133,\ 20.267)$. Since $\sum(y_i-\hat\mu_i) = 0$, $D = 2\sum y_i\log(y_i/\hat\mu_i) = 2(0.5072-0.0662-1.0675+0.7464) = 0.2400$ on 3 df ($p = 0.97$): a perfect fit.
*Without the offset*, the intercept-only model has $\hat\mu_i = \bar y = 9.5$ and $D = 2[3\log\tfrac3{9.5}+5\log\tfrac5{9.5}+9\log\tfrac9{9.5}+21\log\tfrac{21}{9.5}] = 19.008$ ($X^2 = 20.53$; $p = 2.7\times10^{-4}$ against $\chi^2_3$). It treats a longer count as unrelated to a longer exposure.

### Example 4.6: Negative Binomial variance from the Gamma mixture

Take $r = 2$ and $\mu=4$. By Theorem 3.8, $\mathrm{Var}(y) = 4+16/2 = 12$: a dispersion index of 3. The $k=3$ probability from (3.31) is $\frac{\Gamma(5)}{3!\,\Gamma(2)}(\tfrac26)^2(\tfrac46)^3 = 4\cdot\tfrac19\cdot\tfrac8{27} = 0.13169$. Numerical integration of Poisson against the Gamma density gives 0.13169 too, and $10^6$ simulated mixtures have mean 4.006 and variance 12.014.

### Example 4.7: One deviance residual

In Example 4.5's no-offset fit, the $y=21$ observation has $\hat\mu = 9.5$. By (3.20), $d = 2[21\log(21/9.5)-(21-9.5)] = 2(16.6578-11.5) = 10.3157$, so $r^D = +\sqrt{10.3157} = 3.2118$. Pearson: $r^P = 11.5/\sqrt{9.5} = 3.7311$. The deviance residual is smaller because the Poisson's right skew makes a large positive gap less surprising than the symmetric $\sqrt V$ scaling suggests. All four deviance residuals are $(-2.467,\ -1.607,\ -0.164,\ 3.212)$ and their squares sum to $D = 19.008$.

### Example 4.8: AUC on six cases

Positives have scores $(0.8, 0.5, 0.4)$, negatives $(0.6,0.3,0.1)$. Of the nine pairs, $0.8$ beats all three negatives, $0.5$ beats two ($0.3$, $0.1$) and loses to $0.6$, and $0.4$ beats two. That is $3+2+2 = 7$ wins, so $\mathrm{AUC} = 7/9 = 0.7778$, by (3.27). The trapezoid area under the ROC staircase gives the same number.

## 5. Why it matters and where it connects

- **ST-5.** A Bayesian GLM keeps this likelihood. With a Normal prior the log-posterior is (3.13)'s penalized $\ell$; the Laplace approximation to the posterior uses $I(\hat\beta) = X^\top WX$; and Gamma–Poisson (F5 §3.5) is the conjugate special case.
- **ST-6.** Hierarchical (mixed) models add random effects to $\eta_i$; the Negative Binomial is the simplest such model, with the random effect integrated out.
- **ST-9.** An odds ratio from logistic regression is an association. Whether it estimates an effect is ST-9's question; and the odds ratio is not "collapsible" over covariates the way a risk difference is, so its value can change when you adjust for a variable that is not a confounder.
- **ML-3.** Logistic and softmax regression with regularization are this module's likelihoods with a different story: cross-entropy loss (Section 3.7), gradient $X^\top(s-y)$, convexity, ridge as (3.13). Newton and IRLS give a second-order alternative to gradient descent; ML-3 also adds Gaussian discriminant analysis, which produces the same sigmoid.
- **ML-7.** A single sigmoid unit with cross-entropy loss *is* logistic regression, and a softmax output layer is (3.29). Backpropagation begins with the factor $(s-y)$ of (3.7).
- **SM-4 / SM-7.** The log-partition function $b$ is SM-4's $\ln Z$; $b' = \mathbb{E}[y]$, $b'' = \mathrm{Var}(y)$ are its derivatives; the logit is an energy gap and the softmax is the Boltzmann distribution (Section 3.7); Section 3.9 is the maximum-entropy principle.
- **Beyond.** Splines or smooth functions in the linear predictor give generalized additive models (ML-4's kernels are a relative); random effects give mixed models (ST-6); both are still IRLS with a modified penalty or weights.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 drill the exponential-dispersion form, odds and rate ratios, deviance residuals, and a hand IRLS iteration. B1–B8 are the module: B1 (Gamma-log weights and the chain rule), B2 (AUC as a probability), B3 (separation and the ridge step), B6 (why the Bernoulli deviance cannot test fit) and B8 (`project.py`: add the Inverse Gaussian family) are essential. C1–C3 go further: Newton versus Fisher scoring for a non-canonical link, the Gamma GLM versus OLS on $\log y$, and the zero-inflated Poisson.

## 7. Coding project

`project.py` builds `glm_fit` and tests it from every side. Parts: (0) the engine, with links and families as data; (1) the three failures of linear regression, with numbers; (2) SymPy derivation of the family table, canonical links and unit deviances; (3) IRLS against `scipy.optimize` direct maximum likelihood for Normal/identity, Bernoulli/logit, Bernoulli/probit, Poisson/log, Gamma/log, against ST-3's OLS (one pass), and the printed quadratic convergence; (4) Example 4.2 by hand and by code; (5) the logistic inference table, Wald coverage over 1000 replications, odds ratios, the latent-variable constants 1.596, 1.814, 1.701; (6) ROC, AUC against the pairwise probability and `scipy.stats.mannwhitneyu`, ties, calibration, Hosmer–Lemeshow, the decision boundary; (7) separation and the ridge fix; (8) softmax regression by gradient descent against L-BFGS, with finite-difference gradient check; (9) Poisson with offsets, Examples 4.5 and 4.7; (10) overdispersion: Poisson, quasi-Poisson and Negative Binomial coverage; (11) the analysis of deviance, $\chi^2$ calibration under the null, and the Bernoulli failure of Theorem 3.6; (12) residuals, hat matrix, Cook's distance; and "Your turn".

**Build it yourself.** Before reading Part 0, write `glm_fit(X, y, family, link)` for the Poisson/log case from (3.9)–(3.10): about 15 lines. Test it on a Normal/identity problem first (one pass must reproduce `np.linalg.lstsq`), then on Poisson against `scipy.optimize.minimize` of the `scipy.stats.poisson` log-likelihood. Then add the logit link, then the Gamma family. Your turn (a)–(c) are Exercises C2, C3, C1; Exercise B8's reference solution sits at the bottom of the script.

## 8. Where you'll struggle

- **The table is the module.** If any row of Section 3.1's table feels like a memorized fact rather than three lines of algebra, redo it with the Gamma or Inverse Gaussian before moving on. Everything later (weights, deviances, residuals) is read off the table. If ST-1 §3.10 was skipped, go back.
- **The chain rule drops a factor.** Almost everyone forgets $d\mu/d\eta = 1/g'(\mu)$ in step 3 of Theorem 3.2, and the result is a score that is right only for the identity link, which is the one you test on first (Exercise B1(c)). Check yourself: for Poisson with the log link the score is $\sum(y_i-\mu_i)x_i$, with no $1/\mu_i$ weights.
- **IRLS looks like a trick.** It stops looking like one when you see it is Newton's method with the expected Hessian, applied to a log-likelihood whose derivatives you computed. If (3.8) to (3.9) feels like sleight of hand, redo it with Example 4.2's numbers.
- **"Odds ratio" is read as "relative risk" constantly,** including in published papers. Use (3.25) whenever the outcome is common.
- **Complete separation is common in small datasets, and software is quiet about it.** A coefficient of 25 with a standard error of 10⁵ is not a strong effect; it is Theorem 3.4. Check for it whenever a predictor "perfectly predicts" the outcome.
- **Deviance is not $\chi^2$ for Bernoulli data,** and people test with it anyway (Theorem 3.6). Differences of deviances are fine; the single residual deviance on 0/1 data is not a fit test.
- **Overdispersion is the default for real counts.** A Poisson model is a hypothesis; compute $\hat\phi$ every time. Quasi-Poisson fixes standard errors, not the weights.
- **AUC is a ranking metric.** A model with AUC 0.95 can be useless for estimating risks; check calibration separately.
- **Parametrizations.** The Gamma canonical link is $-1/\mu$ (software flips the sign); Negative Binomial software parametrizes with $r$, $1/r$ or $\alpha = 1/r$; "dispersion" sometimes means $\phi$ and sometimes $a(\phi)$. Write down $V(\mu)$ for whatever you are fitting and you will not be confused.

## 9. Check yourself

1. State $\mathbb{E}\,y$ and $\mathrm{Var}(y)$ for the density (3.1) in terms of $b$, and derive them from the score identities.
2. Give $\theta$, $b(\theta)$, $\mu(\theta)$ and $V(\mu)$ for the Poisson and for the Bernoulli.
3. What makes a link canonical, and what three things does it buy you?
4. Write $\partial\ell/\partial\beta_j$ for a general link as a product of four derivatives. Which one gets dropped?
5. What are the working weight and working response for Poisson/log? Why is the IRLS step a weighted least squares problem, and why does $\phi$ never enter?
6. Why is Newton's method equal to Fisher scoring for the canonical link, and what does that imply about concavity and uniqueness?
7. Prove that complete separation destroys the logistic MLE. What does IRLS do?
8. A logistic coefficient is 0.7. What does it do to the odds? To a probability of 0.5? To a probability of 0.01?
9. Derive $\mathrm{AUC} = P(S_1>S_0)$. Why can a model with AUC 0.9 be badly calibrated?
10. Why does a Poisson offset $\log t_i$ have coefficient exactly 1?
11. Show that a Gamma-mixed Poisson has $\mathrm{Var} = \mu+\mu^2/r$. What do quasi-Poisson and the Negative Binomial each do about it?
12. Write the unit deviance as an integral of $(y-t)/V(t)$. Why is the residual deviance on Bernoulli data not $\chi^2_{n-p}$?

## 10. Pair with

- **Agresti, *Foundations of Linear and Generalized Linear Models*, ch. 4–7.** The cleanest end-to-end treatment of exactly this module's logistic, Poisson, deviance and overdispersion material; stronger than these notes on categorical data and software output.
- **McCullagh and Nelder, *Generalized Linear Models*, ch. 2.** The original; dense, but the exponential-family and quasi-likelihood sections are the primary source.
- **Stanford STATS 305B notes.** The same development with a computational and asymptotic emphasis.
- **Dobson and Barnett, *An Introduction to Generalized Linear Models*.** A gentler route with many small worked datasets, if Section 3.3 felt too fast.
- **Gelman, Hill and Vehtari, *Regression and Other Stories*, ch. 13–15.** Applied logistic regression: how to interpret, check and communicate coefficients, which these notes treat only briefly.
