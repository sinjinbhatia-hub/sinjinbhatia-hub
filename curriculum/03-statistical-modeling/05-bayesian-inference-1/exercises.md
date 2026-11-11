# ST-5 Exercises

Conventions as in the notes: $\mathrm{Gamma}(\alpha,\beta)$ is shape–rate; $\mathrm{IG}(a,b)$ has mean $b/(a-1)$; "95% interval" means equal-tailed unless stated. Equation numbers refer to `notes.md`. Quantiles may be computed by software; set up the formula by hand first.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** A prior $p\sim\mathrm{Beta}(1,3)$, then 4 successes in 6 trials. (a) Give the posterior. (b) Give its mean as a weighted average of the prior mean and the sample fraction (3.6), with the weight on the data. (c) Find the probability that the next two trials both succeed, and compare it with (posterior mean)$^2$. Why do they differ?

**A2.** Normal data with known $\sigma^2=9$, prior $\mu\sim\mathcal N(20,16)$, and $n=6$ observations with $\bar y=23$. (a) Find the posterior (3.10) and its 95% credible interval. (b) Find the 95% predictive interval for a new observation (3.11). (c) How many "imaginary observations" is the prior worth ($n_0$)?

**A3.** Poisson counts $4,0,2$ in three unit-length intervals, prior $\lambda\sim\mathrm{Gamma}(3,2)$. (a) Give the posterior. (b) Give the mean and variance of the count in the next interval. (c) Find $P(\tilde y=0\mid y)$.

**A4.** For $y\sim\mathrm{Exponential}(\lambda)$ (rate $\lambda$): (a) Compute $I_1(\lambda)$ and Jeffreys' prior (3.22). (b) Show the posterior for $n$ observations with sum $s$ is $\mathrm{Gamma}(n,s)$. (c) For the data $2.1,\ 0.7,\ 3.4$, compute the posterior mean, compare it with the MLE, and give a 95% credible interval (software for the quantiles).

**A5.** Ten trials, zero successes, uniform prior. (a) Show the posterior is $\mathrm{Beta}(1,11)$ with CDF $1-(1-p)^{11}$. (b) Compute the 95% equal-tailed interval in closed form. (c) Compute the 95% HPD interval. Explain why the HPD interval has $0$ as an endpoint.

## B. Core (do all)

**B1.** (Predictive.) For the posterior Beta$(\alpha_n,\beta_n)$, derive the Beta–Binomial predictive pmf for $m$ future trials (3.8) and its variance, using Eve's law. Then for Beta(9,5) and $m=10$ compute the variance, the ratio to the variance of Binomial$(10,\ 9/14)$, and $P(\tilde y=10)$ (all ten succeed).

**B2.** (NIG and ST-2.) Under the reference prior $p(\mu,\sigma^2)\propto\sigma^{-2}$ and ST-2's ten measurements of $g$ (Example 4.2 there: $\bar y=9.806$, $s=0.04169$): (a) show the posterior predictive for a new measurement is $t_{n-1}\big(\bar y,\ s^2(1+1/n)\big)$. (b) Compute the 95% predictive interval and compare it with the plug-in $\bar y\pm1.96\,s$. (c) A colleague says "the 95% confidence interval $(9.7762,9.8358)$ has a 95% chance of containing $g$." Under which framework is that sentence correct, what is random in each framework, and what is the one fact that makes the two intervals numerically equal here?

**B3.** (Exponential-family conjugacy, a new family.) Write $\mathrm{Exponential}(\lambda)$ in ST-1's form (3.21): identify $T$, $\eta$, $A(\eta)$. (a) Write the conjugate prior (3.18) and show that in the $\lambda$ coordinate it is $\mathrm{Gamma}(\nu+1,\chi)$. (b) Apply (3.19) and check it against A4(b). (c) With prior $\mathrm{Gamma}(2,4)$ and the data of A4(c), find the posterior. Then use (3.20) to write the posterior mean of the *mean lifetime* $1/\lambda=\nabla A$ as a weighted average of a prior value and the sample mean lifetime, and name the weights. What is the posterior mean of $\lambda$ itself, and why is it not of that form? (d) Derive the predictive survival function $P(\tilde X>t\mid y)=\big(\beta_n/(\beta_n+t)\big)^{\alpha_n}$ and evaluate it at $t=1$.

**B4.** (Normal scale.) $y_1,\dots,y_{10}\sim\mathcal N(0,v)$ with $\sum y_i^2=16$. (a) Show Jeffreys' prior for the precision $\phi=1/v$ is $\propto1/\phi$, consistent with $1/v$ for the variance (invariance, Theorem 3.8). (b) With $p(v)\propto1/v$ find the posterior (3.12). (c) Compute its mean, mode, and 95% interval, and compare with $\sum y_i^2/n$.

**B5.** (Gamma–Poisson with exposure.) Events occur at an unknown rate $\lambda$ per hour: 3 events in 2.5 hours, then 9 events in 4 hours. (a) Derive the posterior with a $\mathrm{Gamma}(2,1)$ prior (3.16) and state it. (b) Give the predictive mean, variance and $P(\tilde y=0)$ for the next 1 hour. (c) Compare with a Poisson forecast using the MLE rate.

**B6.** (Decisions.) (a) Prove that the Bayes risk of the posterior mean under squared loss is $\mathbb E_y[\mathrm{Var}(\theta\mid y)]$. (b) With a uniform prior and $y\sim\mathrm{Binomial}(10,p)$, compute the Bayes risk of Laplace's estimator $(y+1)/12$ and of the MLE $y/10$ (average the squared error over the prior and the data). Which is smaller, and why does the prior-averaged MLE risk equal $\mathbb E[p(1-p)]/n$?

**B7.** (Code.) An A/B test: arm A converted 120 of 1000, arm B 145 of 1000, arm C 138 of 1000. Use uniform priors. (a) Give each posterior. (b) Using `prob_B_beats_A` from the bottom of `project.py`, compute $P(\theta_B>\theta_A\mid\text{data})$ and confirm by Monte Carlo. (c) Extend to three arms: estimate $P(\text{each arm is best})$ from 400,000 posterior draws. (d) Repeat (b) with a Beta(20,180) prior on each arm and explain the change.

**B8.** (Bayes factors.) (a) Derive (3.40) from the two marginal densities of $\bar x$. (b) Show that $\mathrm{BF}_{10}$ as a function of $u=n\tau^2/\sigma^2$ is maximized at $u=z^2-1$ with value $e^{(z^2-1)/2}/|z|$, and that no prior can exceed $e^{z^2/2}$. (c) Evaluate both bounds for $z=1.96$ ($p=0.05$) and $z=2.576$ ($p=0.01$). What does the first number say about "$p=0.05$ is strong evidence against $H_0$"?

## C. Challenge (do at least one)

**C1.** (Minimax from a constant-risk Bayes rule.) Let $y\sim\mathrm{Binomial}(n,p)$, squared-error loss. (a) For the prior $\mathrm{Beta}(a,a)$ with $a=\sqrt n/2$, show the posterior mean $\delta(y)=(y+a)/(n+2a)$ has risk $R(p,\delta)=\dfrac{np(1-p)+a^2(1-2p)^2}{(n+2a)^2}$ and that it does not depend on $p$; evaluate it for $n=100$. (b) Prove: if a Bayes rule has constant risk $r$, then it is minimax, i.e. no rule has smaller maximum risk. (c) Compare with the maximum risk of the MLE, $1/(4n)$, and say when each is preferable.

**C2.** (Jeffreys–Lindley, quantified.) With $z=2$, $\sigma=\tau=1$: (a) Find the $n$ at which $\mathrm{BF}_{01}$ first exceeds 1. (b) Show that for large $n$, $\mathrm{BF}_{01}\approx\sqrt u\,e^{-z^2/2}$ and check it against Example 4.10. (c) Repeat (a) for $\tau=0.1$. Extend Part 11 of `project.py` to plot $\mathrm{BF}_{01}$ against $\tau$ on a log grid for $n=10^4$, and mark the maximizer of $\mathrm{BF}_{10}$.

**C3.** (Laplace versus BIC.) For the 2-D logistic data of Part 9: (a) compute $\hat\ell$ at the MLE (call `laplace_logistic` with a huge prior variance) and $-\mathrm{BIC}/2$ (3.42). (b) Compare with the grid log evidence and the Laplace log evidence (3.35). (c) Write out the terms dropped in going from the Laplace expansion to BIC, evaluate each for this dataset, and check that their sum accounts for the gap between $-\mathrm{BIC}/2$ and the Laplace value.
