# ST-4 Exercises

Conventions as in the notes: $\theta,b,a(\phi),V(\mu)$ from (3.1)–(3.3); $w_i$ from (3.6) and $z_i$ from (3.10); $d(y,\mu)$ is the unit deviance (3.18); $r^P,r^D$ are (3.21). $\sigma(t)=1/(1+e^{-t})$. Use a calculator or Python for numbers; the project's helper functions are fair game.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Let $Y\sim\mathrm{Binomial}(m,p)$ and $y = Y/m$.
(a) Write the PMF of $y$ in the form (3.1) and identify $\theta$, $b(\theta)$, $a(\phi)$ and $c$.
(b) Use (3.2) to get $\mathbb{E}\,y$ and $\mathrm{Var}(y)$, and check $\mathrm{Var}(y)$ against $p(1-p)/m$ for $m=10$, $p=0.3$.

**A2.** A logistic fit has $\mathrm{logit}\,p = -3+0.8x$.
(a) Compute $p$ at $x=2$ and at $x=3$.
(b) From the two odds, verify that the odds ratio for one more unit of $x$ is $e^{0.8}$.
(c) Compute the relative risk $p(3)/p(2)$ and the marginal effect $\partial p/\partial x$ at $x=3$.
(d) At which $x$ is $p=\tfrac12$, and what is the largest marginal effect anywhere?

**A3.** Crash counts on road segments are Poisson with $\log\mu = \log t + \beta_0+\beta_1\,\mathrm{wet}$, where $t$ is traffic in million vehicle-km, $\beta_0=-2.0$, $\beta_1=0.5$. For $t=3.2$ find the expected count on a dry and on a wet segment, the rate ratio, and the probability of no crash in each case.

**A4.** Three counts $y=(0,2,7)$ are fitted by $\hat\mu=(1.2,3.1,4.4)$. Compute the unit deviances (use $0\log0=0$), $D$, the Pearson and deviance residuals, and $X^2$.

**A5.** Fit an intercept-only Poisson/log model to $y=(2,4,6)$ by IRLS.
(a) Show that one pass maps $\beta\mapsto\beta+\bar y\,e^{-\beta}-1$.
(b) Start at $\beta=0$ and compute three passes.
(c) What does a pass do when started at $\log\bar y$?

## B. Core (do all; these are the module)

**B1.** (The chain rule and the weights.)
(a) For the Gamma family with log link, find $V$, $g'$, the score (3.5), and $w_i$ (take $\omega_i=1$). What does IRLS reduce to?
(b) For $y=(1,2,6)$ and an intercept-only Gamma/log model, compute the first pass from $\beta=0$ and the fixed point.
(c) For a Poisson model with log link, (3.5) gives $\sum(y_i-\mu_i)x_i$. Now forget the factor $d\mu/d\eta = 1/g'$ of step 3 in Theorem 3.2 and write the resulting "score". Test it on one observation $y=3$ with $\mu=e^\beta$ at $\beta=1$: compare the true $\partial\ell/\partial\beta$ with the wrong one. For which link is the omission invisible?

**B2.** (AUC.) Positives have scores $(0.9,0.6,0.6,0.3)$, negatives $(0.7,0.6,0.2,0.1)$.
(a) Compute the AUC from the pairwise definition, ties counting $\tfrac12$.
(b) Sweep the threshold, list the ROC vertices, and compute the trapezoid area.
(c) Prove in general that the trapezoid area of an empirical ROC equals the Mann–Whitney $U/(n_1n_0)$.
(d) Show that the AUC is unchanged by $s\mapsto s^3$ for positive scores. Why does that make it blind to calibration?

**B3.** (Separation.) The data are $x=1,\dots,6$, $y=000111$, with $X=[\mathbf 1\ x]$.
(a) Prove the logistic MLE does not exist, using $b=(-3.5,\,1)^\top$.
(b) Compute the first IRLS pass from $\beta=0$ by hand.
(c) Compute the first pass of the ridge-penalized IRLS (3.13) with $\lambda=1$ on the slope.

**B4.** (Offsets.)
(a) Derive the coefficient 1 on $\log t$ from the additivity of Poisson counts over disjoint exposure.
(b) A Poisson fit with $\log t$ as a free covariate returns $\hat\beta = 0.93$ with standard error $0.06$. Test $H_0:\beta=1$ and state what it means for the offset model.
(c) If exposures vary but you omit them, name two things that go wrong.

**B5.** (Overdispersion.)
(a) Derive $\mathrm{Var}(y)=\mu+\mu^2/r$ for a Gamma-mixed Poisson using Adam's and Eve's laws.
(b) Counts have mean 5 and variance 12.5. Which $r$ matches?
(c) With that $r$, compare the working weights of the Poisson (and quasi-Poisson) and the Negative Binomial log-link GLM at $\mu=5$ and $\mu=50$. Which observations does the NB downweight relative to Poisson?
(d) By what factor does quasi-Poisson inflate the Poisson standard errors if $\hat\phi=2.5$?

**B6.** (The Bernoulli deviance.)
(a) An intercept-only Bernoulli fit to 50 trials with 15 successes: compute $D$ from (3.18).
(b) Compute $P(\chi^2_{49}>D)$. What would the "goodness-of-fit test" conclude?
(c) Shuffle the data. What is $D$ now? Explain via Theorem 3.6 and name two valid ways to check fit.

**B7.** (Latent variables.)
(a) If $\varepsilon$ is logistic (CDF $\sigma$), show $P(Y=1\mid x)=\sigma(x^\top\beta)$.
(b) If $F(t)=\exp(-e^{-t})$, show $P(Y=1\mid x) = 1-\exp(-e^{x^\top\beta})$ and invert it to find the link.
(c) Derive the slope-matching constant $4\varphi(0)$ between logit and probit coefficients.
(d) A probit coefficient is 0.9. What logit coefficient does the constant 1.7 predict? Compare $\Phi(1.5)$ with $\sigma(1.7\times1.5)$.

**B8.** (Code: extend `project.py`.) Add the Inverse Gaussian family to `FAMILIES`: $V=\mu^3$ and unit deviance $(y-\mu)^2/(y\mu^2)$ (first derive that deviance from (3.19)). Simulate $n=400$ points with $\mu_i=e^{0.5+0.4x_i}$, $\lambda=5$, using `scipy.stats.invgauss.rvs(mu/lam, scale=lam, random_state=rng)` with `rng = np.random.default_rng(23)`. Fit with the log link, and compare $\hat\beta$ with direct maximization of the `scipy.stats.invgauss` log-likelihood. Report $\hat\phi$ against $1/\lambda$.

## C. Challenge (do at least one)

**C1.** (Newton versus Fisher.) For a Poisson model with identity link:
(a) Use (3.11) to show $w_i^{\rm obs}=y_i/\mu_i^2$ while $w_i=1/\mu_i$.
(b) With `glm_fit(..., observed=True)` and `observed=False`, fit $n=100$ points ($x\sim U(0,2)$, $y\sim\mathrm{Poisson}(2+3x)$, `default_rng(22)`) from $\beta=(1,1)$. Print the error $\|\beta_t-\hat\beta\|$ per pass. Describe the two convergence rates.
(c) Compute standard errors from the expected and the observed information at $\hat\beta$. Why are they close here, and when would they differ?

**C2.** (Gamma GLM versus OLS on $\log y$.) Let $y=\mu G$ with $G\sim\mathrm{Gamma}(\nu,\text{mean }1)$.
(a) Show $\mathbb{E}\log y = \log\mu+\psi(\nu)-\log\nu$.
(b) What does OLS of $\log y$ on $x$ estimate when $\log\mu = \beta_0+\beta_1x$? What does $\exp(\text{fitted value})$ estimate?
(c) Simulate $n=500$, $\nu=3$, $\mu=e^{1+0.5x}$, $x\sim\mathcal{N}(0,1)$ (`default_rng(20)`), and compare the Gamma/log GLM with OLS on $\log y$, including $\overline{\exp(\text{OLS fit})}/\bar y$ against $e^{\psi(\nu)-\log\nu}$.
(d) Show that the ML shape solves $\log\nu-\psi(\nu) = D/(2n)$ and compute it.

**C3.** (Zero-inflated Poisson.)
(a) Write $\ell(\beta,\pi)$ for $P(0)=\pi+(1-\pi)e^{-\mu}$, $P(k)=(1-\pi)\mathrm{Pois}(k;\mu)$.
(b) Show $\mathbb{E}\,y=(1-\pi)\mu$ and $\mathrm{Var}(y)=(1-\pi)\mu(1+\pi\mu)$.
(c) Simulate $n=1500$ with $x\sim\mathcal{N}(0,1)$, $\mu=e^{0.8+0.5x}$ and $\pi=0.3$ (`default_rng(21)`), and fit by `scipy.optimize.minimize`. Compare the observed fraction of zeros with the Poisson GLM's expectation $\overline{e^{-\hat\mu}}$, and the two log-likelihoods.
