# ST-6: Bayesian Computation: Metropolis–Hastings, Gibbs Sampling, Hamiltonian Monte Carlo, Diagnostics, Hierarchical Models, Variational Inference

**Time:** ~14 hours · **Prerequisites:** ST-5 (posteriors, NIG, predictive checks); F5 (Markov chains §3.12, detailed balance Theorem 3.41, KL Definition 3.44); CM-8 (Hamilton's equations, Liouville Theorem 3.4, Störmer–Verlet §3.11) for §3.4 · **Unlocks:** SM-7, ML-10, ST-8

> **Conventions in this module.** The **target** is the posterior $\pi(\theta)\propto p(y\mid\theta)p(\theta)$, known only up to its normalizing constant; "unnormalized" is fine everywhere because only ratios $\pi(\theta')/\pi(\theta)$ appear. A **proposal** density is $q(\theta'\mid\theta)$ (move from $\theta$ to $\theta'$). The Metropolis–Hastings acceptance probability is $\alpha(\theta,\theta')=\min\!\big(1,\ \pi(\theta')q(\theta\mid\theta')/[\pi(\theta)q(\theta'\mid\theta)]\big)$. A chain is $\theta^{(1)},\dots,\theta^{(S)}$; **warm-up** (burn-in) is the discarded start; **thinning** keeps every $k$-th draw. $\rho_k$ (always with a subscript) is the lag-$k$ autocorrelation of a chain; $n_{\rm eff}$ is the effective sample size; $\hat R$ is the potential scale reduction. In HMC, $\theta$ is position, $\rho$ (no subscript) is momentum, $M$ the mass matrix (here $M=I$ unless stated), $U(\theta)=-\log\pi(\theta)$ the potential energy, $H(\theta,\rho)=U(\theta)+\tfrac12\rho^\top M^{-1}\rho$ the Hamiltonian, $\varepsilon$ the leapfrog step and $L$ the number of steps. In VI, $q_\phi$ is a member of the variational family with parameters $\phi$ and the ELBO is $\mathbb E_q[\log p(\theta,y)-\log q(\theta)]$. As in ST-5: $\mathcal N(m,v)$ takes the variance, $\mathrm{IG}(a,b)\propto x^{-a-1}e^{-b/x}$ has mean $b/(a-1)$, NIG$(\mu_0,\kappa_0,\alpha_0,\beta_0)$ means $\sigma^2\sim\mathrm{IG}(\alpha_0,\beta_0)$, $\mu\mid\sigma^2\sim\mathcal N(\mu_0,\sigma^2/\kappa_0)$. In the hierarchical sections $\tau$ is a between-group **standard deviation**. $\Phi,\varphi$ are the standard Normal CDF and density.

## 0. Where this sits

You need ST-5's machinery (Bayes' theorem for parameters, the NIG update Theorem 3.3, Bayesian regression (3.30), posterior predictive checks §3.10), F5's Markov chains (stationarity (3.26), Theorem 3.40 on convergence, Theorem 3.41 that detailed balance implies stationarity) and F5's KL divergence. For HMC you need CM-8: Hamilton's equations (3.5), Liouville's theorem (Theorem 3.4), and the kick–drift–kick integrator of §3.11. LA-14 supplies gradients; LA-13's preconditioning reappears as the mass matrix.

This module turns "the posterior is known up to a constant" into numbers. SM-7 reads it backwards (the Boltzmann distribution as a posterior, the ELBO as a free energy), ML-10 uses the ELBO and the reparametrization trick of §3.8 to train VAEs, and ST-8 samples latent states with the Gibbs and HMC samplers built here.

## 1. The problem that forced this idea

The first Markov chain Monte Carlo run was physics. Metropolis, Rosenbluth, Rosenbluth, Teller and Teller (1953) wanted equilibrium properties of a gas of hard disks, averages over the Boltzmann distribution $e^{-E/kT}/Z$ in a space of hundreds of coordinates, where $Z$ is hopeless to compute. On the MANIAC computer at Los Alamos they ran a random walk that moved one particle at a time and accepted moves with probability $\min(1,e^{-\Delta E/kT})$; $Z$ never appeared (SM-6 builds this for the Ising model). Hastings (1970) generalized the acceptance rule to asymmetric proposals and wrote it for statisticians, who mostly did not notice. Geman and Geman (1984) restored noisy images with a sampler that updated one pixel at a time from its exact conditional distribution and named it the Gibbs sampler. Duane, Kennedy, Pendleton and Roweth (1987) needed to sample gauge fields in lattice QCD and combined Hamiltonian dynamics with a Metropolis correction: "hybrid Monte Carlo". Neal brought it to statistics and neural networks in the 1990s.

Gelfand and Smith (1990) showed statisticians that the Gibbs sampler computed posteriors nobody could integrate, and applied Bayesian statistics changed within a decade: BUGS automated Gibbs sampling, and Stan (2012) made HMC with automatic differentiation the default. Variational inference came from the same physics: the mean-field approximation (SM-6, SM-7) entered machine learning in the 1990s through Jordan, Hinton and coauthors, and Blei and coauthors made it standard for large models.

The forcing problem is the one ST-5 stopped at. Conjugate models are the only ones whose posterior has a closed form. Logistic regression, any hierarchical model, any model with more than a handful of parameters, has a posterior you can evaluate pointwise up to the constant $p(y)$ and nothing more. Every quantity you report (a posterior mean, an interval, a predictive probability) is an integral of that function over a space of tens to thousands of dimensions, where grids are impossible ($100$ points per axis in 20 dimensions is $10^{40}$ evaluations). Two general escapes exist: draw samples from the posterior and average (MCMC), or replace the posterior by the closest member of a tractable family (VI). This module builds both, from scratch, and the tools that tell you whether they worked.

## 2. Intuition

**A walker that spends time in proportion to height.** Picture the posterior as a landscape. You want a walker whose long-run fraction of time near $\theta$ equals $\pi(\theta)\,d\theta$. Then the walker's path is a sample, and averaging any $f$ along the path gives $\mathbb E_\pi f$. Design rule: propose a step; if it goes uphill, take it; if it goes downhill by a factor $r=\pi(\theta')/\pi(\theta)<1$, take it with probability $r$. Why should that work? Look at two points $a,b$ with $\pi(b)=\frac12\pi(a)$ and a symmetric proposal. Flow $a\to b$: (time at $a$) × (propose $b$) × $\frac12$. Flow $b\to a$: (time at $b$) × (propose $a$) × $1$. If time at $a$ is twice time at $b$, the two flows match, and nothing changes. That balance, pair by pair, is F5's detailed balance (3.27). **Guess the result:** for an asymmetric proposal the acceptance must also correct for how much easier it is to propose $b$ from $a$ than $a$ from $b$; the ratio to use is $\pi(b)q(a\mid b)/[\pi(a)q(b\mid a)]$.

**Step size is a trade-off.** Tiny steps are almost always accepted but the walker diffuses: after $n$ steps of size $s$ it has moved about $s\sqrt n$. Huge steps land in the tails and get rejected; the walker stands still. The best is in between, at a moderate acceptance rate. In high dimension the right step shrinks, because a random direction in $d$ dimensions is almost always downhill.

**Gibbs: walk along the axes.** If you can sample each coordinate exactly given the others, do that: no rejection. In a round posterior this mixes fast. In a long diagonal ridge, each axis-aligned move is pinned between the ridge walls and the walker staircases along it in tiny steps:

```
  theta2                      ridge (correlation 0.99)
    ^                    _/
    |                 _/|_
    |              _/   |_|      each move is along an axis,
    |           _/|_|            limited by the ridge width
    |        _/                 
    +-------------------------> theta1
```

**HMC: give the walker momentum.** Treat $U=-\log\pi$ as a potential energy surface and the walker as a frictionless puck. Flick it with a random momentum and let Newton's laws run. It coasts along the valley, converting potential to kinetic energy and back, and ends far from where it started, at a point of the same total energy. Exact Hamiltonian flow conserves energy and phase-space volume (CM-8), so the endpoint would always be accepted; the only rejection comes from the integrator's energy error. **Guess the result:** the acceptance probability is $\min(1,e^{-\Delta H})$, and since the integrator is second order, $\Delta H=O(\varepsilon^2)$.

**Did it work?** Start several walkers far apart. If they have forgotten their starts, the spread between their averages matches the spread within each. Strongly correlated draws also carry less information than independent ones: 10,000 draws with lag-one correlation 0.9 are worth about 500. **Guess:** $n_{\rm eff}=S/(1+2\sum_k\rho_k)$.

**Hierarchical models are where this bites.** Eight schools each estimate a coaching effect with noise. Each school's true effect $\theta_j$ is drawn from a population with spread $\tau$. If $\tau$ is small, all $\theta_j$ must be near the common mean; if $\tau$ is large, they can be anywhere. In the $(\theta_j,\log\tau)$ plane the posterior is a funnel: wide at the top, a thin neck at the bottom. No single step size works in both places. The fix is a change of coordinates: write $\theta_j=\mu+\tau\eta_j$ and sample $\eta_j$, which is $\mathcal N(0,1)$ a priori whatever $\tau$ is. The funnel becomes a cylinder.

**VI: optimize instead of sample.** Pick a family of easy distributions (independent Gaussians, say) and find the member closest to $\pi$. "Closest" must be measured in a way you can compute without knowing $p(y)$; the direction $\mathrm{KL}(q\Vert\pi)$ works because its expectation is over $q$. That direction penalizes $q$ for putting mass where $\pi$ has none, so $q$ hides inside $\pi$: it underestimates spread.

## 3. Formal development

### 3.1 Why draws suffice; rejection and importance sampling

Given draws $\theta^{(s)}\sim\pi$, every posterior summary is an average:
$$\mathbb E_\pi[f(\theta)]\approx\frac1S\sum_{s=1}^Sf(\theta^{(s)}),\qquad P(\theta\in A\mid y)\approx\frac1S\#\{s:\theta^{(s)}\in A\},\qquad \tilde y^{(s)}\sim p(\tilde y\mid\theta^{(s)}). \tag{3.1}$$
Quantiles of the draws give credible intervals; the last item samples ST-5's predictive (3.3) by composition. For independent draws the error is $\mathrm{sd}_\pi(f)/\sqrt S$ in any dimension: no curse of dimension *once you have the draws*. The difficulty is producing them.

**Rejection sampling.** If $\pi(\theta)\le Mg(\theta)$ for a density $g$ you can sample, draw $\theta\sim g$ and accept with probability $\pi(\theta)/(Mg(\theta))$. Accepted draws are exactly $\pi$-distributed, and the acceptance rate is $1/M$ (for normalized $\pi$). **Importance sampling** keeps every draw from $g$ and weights it: $\mathbb E_\pi f\approx\sum w_sf(\theta^{(s)})/\sum w_s$ with $w_s=\pi(\theta^{(s)})/g(\theta^{(s)})$. The normalizing constant cancels in the ratio. Its effective sample size is $(\sum w_s)^2/\sum w_s^2\approx S\,(\mathbb E_gw)^2/\mathbb E_gw^2$. $\qquad(3.2)$

**The curse, quantified.** Take $\pi=\mathcal N(0,I_d)$ and $g=\mathcal N(0,s^2I_d)$ with $s>1$, a proposal only slightly too wide. The ratio $\pi/g$ peaks at $0$ with value $s^d$, so rejection accepts a fraction $s^{-d}$. For importance sampling, per coordinate $\mathbb E_gw^2=\int\pi^2/g=s/\sqrt{2-1/s^2}$ (a Gaussian integral), and coordinates multiply:
$$\text{acceptance}=s^{-d},\qquad \frac{\mathrm{ESS}}{S}\approx\Big(\frac{\sqrt{2-1/s^2}}{s}\Big)^d. \tag{3.3}$$
With $s=1.5$: at $d=10$, $0.0173$ and $0.158$; at $d=100$, $2.46\times10^{-18}$ and $9.66\times10^{-9}$ (project Part 1, which checks (3.3) by simulation at $d=1,5$). Both decay exponentially because the overlap of two slightly different high-dimensional distributions is exponentially small. Any method that proposes globally and blindly fails; MCMC proposes *locally* and lets the chain find the mass. (Importance sampling survives in one place in this module: PSIS-LOO, §3.7, where the proposal is the full posterior and the target is a leave-one-out posterior that differs from it only slightly.)

### 3.2 Markov chain Monte Carlo and Metropolis–Hastings

**The idea.** Build a Markov chain with transition kernel $K(\theta\to\theta')$ whose stationary distribution is $\pi$: $\int\pi(\theta)K(\theta\to\theta')\,d\theta=\pi(\theta')$. F5 Theorem 3.40 (for finite chains; the general-state-space version is in Robert and Casella ch. 6, cited not proved) says that if the chain is irreducible (it can reach any set of positive $\pi$-probability) and aperiodic, then the distribution of $\theta^{(t)}$ converges to $\pi$ from any start, and time averages converge to $\pi$-expectations (the ergodic theorem). Stationarity is the property we design for. The easiest sufficient condition is **detailed balance**,
$$\pi(\theta)K(\theta\to\theta')=\pi(\theta')K(\theta'\to\theta)\quad\text{for all }\theta\neq\theta'. \tag{3.4}$$
Integrating (3.4) over $\theta$ gives stationarity in one line (F5 Theorem 3.41, with sums replaced by integrals).

**Theorem 3.1 (Metropolis–Hastings).** Propose $\theta'\sim q(\cdot\mid\theta)$ and accept with probability
$$\alpha(\theta,\theta')=\min\left(1,\ \frac{\pi(\theta')\,q(\theta\mid\theta')}{\pi(\theta)\,q(\theta'\mid\theta)}\right); \tag{3.5}$$
on rejection stay at $\theta$. The resulting chain satisfies (3.4), so $\pi$ is stationary.

▶ **Try it first.** For $\theta'\neq\theta$ the kernel is $K(\theta\to\theta')=q(\theta'\mid\theta)\alpha(\theta,\theta')$ (the rejection mass sits at $\theta'=\theta$, where (3.4) is trivial). Write (3.4) with an unknown $\alpha$, and look for the largest acceptance probabilities that satisfy it.

**Derivation.**
1. Detailed balance for $\theta'\ne\theta$ requires $\pi(\theta)q(\theta'\mid\theta)\alpha(\theta,\theta')=\pi(\theta')q(\theta\mid\theta')\alpha(\theta',\theta)$.
2. So the ratio is forced: $\dfrac{\alpha(\theta,\theta')}{\alpha(\theta',\theta)}=\dfrac{\pi(\theta')q(\theta\mid\theta')}{\pi(\theta)q(\theta'\mid\theta)}=:r(\theta,\theta')$, and $r(\theta',\theta)=1/r(\theta,\theta')$.
3. Any pair of probabilities with ratio $r$ works. To move as often as possible, make the larger of the two equal to 1. If $r\ge1$: $\alpha(\theta,\theta')=1$ and $\alpha(\theta',\theta)=1/r$. If $r<1$: $\alpha(\theta,\theta')=r$ and $\alpha(\theta',\theta)=1$. Both cases are $\alpha(\theta,\theta')=\min(1,r)$.
4. Verify: if $r<1$, LHS of step 1 is $\pi(\theta)q(\theta'\mid\theta)\,r=\pi(\theta')q(\theta\mid\theta')$, and RHS is $\pi(\theta')q(\theta\mid\theta')\cdot1$. Equal. The case $r\ge1$ is the same with roles swapped. ∎

■ **What just happened.** Detailed balance fixes only the *ratio* of forward and backward acceptance. Choosing the larger one to be 1 is Peskun-optimal: any other valid choice rejects more. $\pi$ enters only through $\pi(\theta')/\pi(\theta)$, so $p(y)$ cancels. In code, compare $\log u<\log r$ to avoid underflow.

**Special cases.**
- **Metropolis** (symmetric $q$, $q(\theta'\mid\theta)=q(\theta\mid\theta')$): $\alpha=\min(1,\pi(\theta')/\pi(\theta))$. This is the 1953 rule with $\pi\propto e^{-E/kT}$.
- **Random-walk Metropolis** (RWM): $\theta'=\theta+s\,C z$, $z\sim\mathcal N(0,I)$, with scale $s$ and a fixed matrix $C$ (often a Cholesky factor of a covariance guess). Symmetric, so the Metropolis rule applies.
- **Independence sampler:** $q(\theta'\mid\theta)=g(\theta')$, so $r=w(\theta')/w(\theta)$ with §3.1's weights $w=\pi/g$; it inherits importance sampling's curse.
- **Proposals on a transformed scale need the Hastings factor.** Example 4.1 shows a multiplicative proposal for a positive parameter: forget the factor and you sample the wrong distribution (the project gets mean 1.98 instead of 3).

**Tuning the scale.** Small $s$: acceptance near 1, but the chain diffuses with step $s$, and the time to cross a posterior of width 1 grows like $1/s^2$. Large $s$: nearly every proposal lands in the tails and is rejected. Roberts, Gelman and Gilks (1997) proved, for targets that are products of $d$ identical one-dimensional densities, that as $d\to\infty$ the optimal scale is $s\approx2.38/\sqrt d$ (in units of the target sd) and the corresponding acceptance rate is $0.234$; in one dimension the optimum is near $0.44$. The proof is beyond this module (it is a diffusion limit; see the paper). Two consequences: the step must shrink like $d^{-1/2}$, so the chain needs $O(d)$ iterations per independent draw; and acceptance between about 0.15 and 0.4 is near-optimal, since the efficiency curve is flat at the top. The project's $d=20$ sweep peaks at acceptance $0.21$. In $d=2$ with an isotropic proposal on a $\rho=0.9$ Gaussian the curve is flat from 0.09 to 0.3: the narrow direction limits the step, the long one the crossing time. The fix is a proposal covariance shaped like the posterior.

**Adaptive scaling during warm-up.** Adjust $s$ (e.g. $\log s\leftarrow\log s+\gamma_t(\alpha_t-0.234)$, $\gamma_t$ decreasing) and estimate a proposal covariance from the draws so far; then *freeze* the tuning and sample. Adapting forever can break stationarity, since the kernel then depends on the chain's history.

### 3.3 Gibbs sampling

**Algorithm.** Split $\theta=(\theta_1,\dots,\theta_k)$ into blocks. One sweep: for $j=1,\dots,k$, draw
$$\theta_j\sim\pi(\theta_j\mid\theta_{-j}),\tag{3.6}$$
the **full conditional** given the current values of all other blocks.

**Theorem 3.2.** Each update (3.6) is a Metropolis–Hastings step with proposal $q(\theta'\mid\theta)=\pi(\theta'_j\mid\theta_{-j})$ (leaving $\theta_{-j}$ unchanged), and its acceptance probability is 1.

▶ **Try it first.** Put this $q$ into (3.5). Write $\pi(\theta)=\pi(\theta_j\mid\theta_{-j})\,\pi(\theta_{-j})$ and note that $\theta'_{-j}=\theta_{-j}$.

**Derivation.**
1. The proposal changes only block $j$: $\theta'=(\theta'_j,\theta_{-j})$. So $q(\theta'\mid\theta)=\pi(\theta'_j\mid\theta_{-j})$ and $q(\theta\mid\theta')=\pi(\theta_j\mid\theta_{-j})$.
2. Factor the target: $\pi(\theta')=\pi(\theta'_j\mid\theta_{-j})\pi(\theta_{-j})$ and $\pi(\theta)=\pi(\theta_j\mid\theta_{-j})\pi(\theta_{-j})$.
3. The ratio in (3.5) is $\dfrac{\pi(\theta'_j\mid\theta_{-j})\pi(\theta_{-j})\cdot\pi(\theta_j\mid\theta_{-j})}{\pi(\theta_j\mid\theta_{-j})\pi(\theta_{-j})\cdot\pi(\theta'_j\mid\theta_{-j})}=1$. ∎

■ **What just happened.** Gibbs is MH with a perfect proposal for one block at a time. Each block update leaves $\pi$ invariant, so the sweep (a composition) does too. A fixed-order sweep is not itself reversible, but stationarity is what matters; random-scan Gibbs is reversible.

**When conditionals are available.** The full conditional is the joint density viewed as a function of one block, so you read it off by keeping only the factors that contain that block. Non-conjugate models often have conjugate pieces.

**The Normal model, $\mu$ and $\sigma^2$ unknown.** Prior NIG$(\mu_0,\kappa_0,\alpha_0,\beta_0)$. The joint density in $(\mu,\sigma^2)$ is
$$p(\mu,\sigma^2,y)\propto(\sigma^2)^{-\alpha_0-1-\frac{n+1}2}\exp\Big\{-\frac{1}{2\sigma^2}\Big[2\beta_0+\sum_i(y_i-\mu)^2+\kappa_0(\mu-\mu_0)^2\Big]\Big\}. \tag{3.8}$$
1. *Fix $\sigma^2$.* As a function of $\mu$, the bracket is $\kappa_n(\mu-\mu_n)^2+$ const by ST-5's completing-the-square step (Theorem 3.3, step 3), with $\kappa_n=\kappa_0+n$, $\mu_n=(\kappa_0\mu_0+n\bar y)/\kappa_n$. So
$$\mu\mid\sigma^2,y\sim\mathcal N(\mu_n,\ \sigma^2/\kappa_n). \tag{3.9}$$
2. *Fix $\mu$.* As a function of $\sigma^2$, (3.8) is $(\sigma^2)^{-a-1}e^{-b/\sigma^2}$, an inverse-Gamma kernel:
$$\sigma^2\mid\mu,y\sim\mathrm{IG}\Big(\alpha_0+\tfrac{n+1}2,\ \beta_0+\tfrac12\Big[\textstyle\sum_i(y_i-\mu)^2+\kappa_0(\mu-\mu_0)^2\Big]\Big). \tag{3.10}$$
The $\frac{n+1}2$ counts $n$ data terms plus the prior term for $\mu$, which also carries $\sigma^2$. This model has a closed-form posterior (ST-5), which is why it is the test case: the project's Gibbs draws pass Kolmogorov–Smirnov tests against the exact $t$ marginal of $\mu$ and IG marginal of $\sigma^2$. Replace the NIG prior by independent priors $\mu\sim\mathcal N(m_0,v_0)$, $\sigma^2\sim\mathrm{IG}$ (the "semi-conjugate" prior) and the posterior has no closed form, but the conditionals are still Normal and IG.

**Bayesian linear regression with unknown $\sigma^2$.** With $\beta\mid\sigma^2\sim\mathcal N(0,\sigma^2V_0)$ and $\sigma^2\sim\mathrm{IG}(a_0,b_0)$, the same two steps give $\beta\mid\sigma^2,y\sim\mathcal N(V_nX^\top y,\ \sigma^2V_n)$ with $V_n=(X^\top X+V_0^{-1})^{-1}$ (ST-5 (3.30)), and $\sigma^2\mid\beta,y\sim\mathrm{IG}\big(a_0+\frac{n+p}2,\ b_0+\frac12[\|y-X\beta\|^2+\beta^\top V_0^{-1}\beta]\big)$, with $p$ the number of coefficients. Project Part 3 matches the exact posterior means to three decimals.

**Data augmentation: probit regression (Albert and Chib, 1993).** The probit model $P(y_i=1)=\Phi(x_i^\top\beta)$ has no conjugate prior. ST-4 (§3.7) wrote it as a latent-variable model: $z_i=x_i^\top\beta+\varepsilon_i$, $\varepsilon_i\sim\mathcal N(0,1)$, $y_i=\mathbb 1[z_i>0]$. Add the $z_i$ to the parameters and sample $(\beta,z)$ jointly; the marginal of $\beta$ is the posterior you want.

▶ **Try it first.** Write the joint $p(\beta,z,y)=p(\beta)\prod_i\varphi(z_i-x_i^\top\beta)\,\mathbb 1[\text{sign of }z_i\text{ agrees with }y_i]$. Read off each conditional.

**Derivation.**
1. *Given $\beta$ and $y$*, the factors containing $z_i$ are $\varphi(z_i-x_i^\top\beta)\mathbb 1[\cdot]$: a Normal truncated to the side $y_i$ dictates,
$$z_i\mid\beta,y\sim\mathcal N(x_i^\top\beta,1)\ \text{truncated to }(0,\infty)\text{ if }y_i=1,\ (-\infty,0]\text{ if }y_i=0, \tag{3.11}$$
independently across $i$.
2. *Given $z$*, $y$ is a function of $z$, so it carries no further information about $\beta$, and $z=X\beta+\varepsilon$ is ordinary linear regression with known variance 1. With prior $\beta\sim\mathcal N(0,B_0)$, ST-5 (3.30) gives
$$\beta\mid z,y\sim\mathcal N\big(BX^\top z,\ B\big),\qquad B=(X^\top X+B_0^{-1})^{-1}. \tag{3.12}$$
3. Sampling (3.11) by inversion: with $m=x_i^\top\beta$, if $y_i=1$ draw $u\sim\mathrm U(\Phi(-m),1)$, if $y_i=0$ draw $u\sim\mathrm U(0,\Phi(-m))$; set $z_i=m+\Phi^{-1}(u)$. ∎

■ **What just happened.** Adding latent variables turned a non-conjugate model into two conjugate ones. The price is a larger state space and correlation between $z$ and $\beta$, which slows mixing when the data are very informative. The same trick drives Gibbs samplers for mixtures (latent labels), missing data, and ST-8's state-space models (latent states).

**The failure mode: correlation.** Take a bivariate Normal target with unit variances and correlation $\rho$. The conditional is $\theta_1\mid\theta_2\sim\mathcal N(\rho\theta_2,1-\rho^2)$. One sweep gives $\theta_1^{(t+1)}=\rho\theta_2^{(t)}+\dots=\rho^2\theta_1^{(t)}+\text{noise}$: the $\theta_1$ chain is an AR(1) process with coefficient $\rho^2$. At $\rho=0.99$, $\rho_1=0.9801$ and (Example 4.5's formula) $n_{\rm eff}/S=(1-\rho^2)/(1+\rho^2)=1/99.5$. The project measures lag-one autocorrelation $0.9797$ and $n_{\rm eff}/S=0.0105$. Hierarchical models create exactly this correlation between group effects and their variance (§3.6). Remedies: block correlated parameters together, reparametrize, or switch to HMC.

**Metropolis-within-Gibbs.** If one conditional is not a known distribution, replace that block's exact draw with any MH step that leaves the conditional invariant (a random walk on that block alone). The composite still leaves $\pi$ invariant, by the argument after Theorem 3.2.

### 3.4 Hamiltonian Monte Carlo

**Augment with momentum.** Introduce $\rho\in\mathbb R^d$ independent of $\theta$ with $\rho\sim\mathcal N(0,M)$, and define
$$H(\theta,\rho)=U(\theta)+K(\rho),\qquad U=-\log\pi(\theta),\qquad K=\tfrac12\rho^\top M^{-1}\rho. \tag{3.13}$$
This is CM-8's Hamiltonian $T+V$ with potential energy $-\log\pi$: SM-7's dictionary "energy = negative log-probability" with temperature 1. The joint density $p(\theta,\rho)\propto e^{-H}=\pi(\theta)\,e^{-K(\rho)}$ factors, so integrating out $\rho$ leaves $\pi(\theta)$: **if we sample the joint, the $\theta$-coordinates are posterior draws**. Hamilton's equations (CM-8 (3.5)) are
$$\dot\theta=\partial H/\partial\rho=M^{-1}\rho,\qquad\dot\rho=-\partial H/\partial\theta=\nabla\log\pi(\theta).$$
The force is the gradient of the log posterior: the puck is pushed toward high density.

**Theorem 3.3 (Exact dynamics are always accepted).** Let $\Phi_T$ be the time-$T$ Hamiltonian flow and $F(\theta,\rho)=(\theta,-\rho)$ the momentum flip. The deterministic proposal $z\mapsto G(z)=F(\Phi_T(z))$, accepted with probability $\min\!\big(1,\ e^{-H(G(z))}/e^{-H(z)}\big)$, leaves $e^{-H}$ invariant, and the acceptance probability is 1.

▶ **Try it first.** You need three facts from CM-8: $H$ is constant along the flow; the flow preserves volume (Liouville); running the flow, flipping momentum, running again, and flipping returns you to the start (time-reversibility, true because $H$ is even in $\rho$). Use them to show $G\circ G=$ identity, $|\det DG|=1$, and $H\circ G=H$.

**Derivation.**
1. *Involution.* Reversibility says $F\circ\Phi_T\circ F=\Phi_T^{-1}$ (check: flip, run forward, flip is running backward, because $\dot\theta=M^{-1}\rho$ changes sign under the flip and $\dot\rho$ does not). So $G\circ G=F\Phi_TF\Phi_T=\Phi_T^{-1}\Phi_T=I$.
2. *Volume.* $|\det D\Phi_T|=1$ (Liouville, CM-8 Theorem 3.4), and $|\det DF|=1$. So $|\det DG|=1$.
3. *Detailed balance for a deterministic involution.* Let $p=e^{-H}$ and $a(z)=\min(1,p(G(z))/p(z))$. For sets $A,B$, the probability of being in $A$ and moving to $B$ is $\int_Ap(z)a(z)\mathbb 1[G(z)\in B]\,dz$. Substitute $z=G(w)$ ($dz=dw$ by step 2; $G(z)=w$ by step 1): it becomes $\int_Bp(G(w))a(G(w))\mathbb 1[G(w)\in A]\,dw$. Now $p(z)a(z)=\min(p(z),p(G(z)))$ is symmetric under $z\leftrightarrow G(z)$, so $p(G(w))a(G(w))=p(w)a(w)$, and the expression is the probability of being in $B$ and moving to $A$. Detailed balance holds.
4. *Acceptance.* Energy conservation gives $H(\Phi_T(z))=H(z)$, and $H$ is even in $\rho$, so $H(G(z))=H(z)$ and $a=1$.
5. The momentum refresh $\rho\sim\mathcal N(0,M)$ is a Gibbs step on $\rho$ given $\theta$, so it also preserves $e^{-H}$. Without it the chain would stay on one energy level. ∎

■ **What just happened.** The three properties of Hamiltonian flow each do a job: reversibility makes the proposal its own inverse (so a Metropolis test is possible), volume preservation removes the Jacobian from that test, and energy conservation makes the test always pass. The flip $F$ is bookkeeping: in practice $\rho$ is resampled next anyway, so you never apply it. If you replace the exact flow by an integrator that keeps the first two properties, only the third breaks, and the Metropolis test corrects for it exactly.

**The leapfrog integrator.** CM-8 §3.11 built Störmer–Verlet as a composition of exact flows. Split $H=U+K$. The flow of $U$ alone (a "kick") moves only $\rho$: $\rho\leftarrow\rho+t\nabla\log\pi(\theta)$. The flow of $K$ alone (a "drift") moves only $\theta$: $\theta\leftarrow\theta+tM^{-1}\rho$. Leapfrog is half kick, drift, half kick:
$$\rho_{1/2}=\rho+\tfrac\varepsilon2\nabla\log\pi(\theta),\qquad\theta'=\theta+\varepsilon M^{-1}\rho_{1/2},\qquad\rho'=\rho_{1/2}+\tfrac\varepsilon2\nabla\log\pi(\theta'). \tag{3.14}$$
1. *Volume.* Each sub-step is a shear: the kick has Jacobian $\begin{pmatrix}I&0\\ \frac\varepsilon2\nabla^2\log\pi&I\end{pmatrix}$, the drift $\begin{pmatrix}I&\varepsilon M^{-1}\\0&I\end{pmatrix}$. Each has determinant 1, so the composition does. (It is in fact symplectic, CM-8 (3.20), since each piece is an exact Hamiltonian flow.)
2. *Reversibility.* The composition is symmetric (kick–drift–kick reads the same backwards), so flipping $\rho'$ and applying (3.14) returns to $(\theta,-\rho)$. The project checks both: Jacobian determinant $1.00000000$ for 40 leapfrog steps on the banana density, against $2.62$ for explicit Euler; reversal error $1.4\times10^{-16}$.
3. *Energy is only approximately conserved.* Leapfrog is second order: the local error per step is $O(\varepsilon^3)$ and over a trajectory of fixed length $T=L\varepsilon$ the error in $(\theta,\rho)$ is $O(\varepsilon^2)$.

For a symplectic integrator the energy error does *not* accumulate over long trajectories. The cleanest way to see it is the Gaussian case, $U=\theta^2/2$ ($M=1$), where the leapfrog map exactly conserves the **shadow Hamiltonian**
$$\tilde H(\theta,\rho)=\tfrac12\rho^2+\tfrac12\big(1-\tfrac{\varepsilon^2}4\big)\theta^2. \tag{3.15}$$
(Expand one step of (3.14) and substitute; the project verifies $\tilde H(\theta',\rho')-\tilde H(\theta,\rho)=0$ in SymPy.) Then $H-\tilde H=\frac{\varepsilon^2}8\theta^2$, so $|H(\text{end})-H(\text{start})|=\frac{\varepsilon^2}8|\theta_{\rm end}^2-\theta_{\rm start}^2|$, bounded by $O(\varepsilon^2)$ for all $L$. For general smooth $U$, backward error analysis (Hairer, Lubich and Wanner, *Geometric Numerical Integration*, ch. IX; cited) shows a shadow Hamiltonian $\tilde H=H+O(\varepsilon^2)$ exists for symplectic methods, with the same consequence. The project measures the maximum $|\Delta H|$ over a trajectory of length $T=2$ on the $\rho=0.9$ Gaussian for $\varepsilon=0.2,0.1,0.05,0.025$: $0.556, 0.138, 0.0351, 0.00877$, log–log slope $1.994$.

The leapfrog step is stable only if $\varepsilon$ is below $2\times$ (the smallest standard deviation of the target, in the metric $M$): on the $\rho=0.9$ Gaussian the narrow direction has sd $\sqrt{1-\rho}=0.316$, so $\varepsilon<0.632$. Above that the energy error grows exponentially and every proposal is rejected.

**The algorithm.**
$$\boxed{\begin{array}{l}\text{for }s=1,\dots,S:\\ \quad\rho\sim\mathcal N(0,M);\ (\theta^*,\rho^*)\leftarrow L\text{ leapfrog steps (3.14) from }(\theta^{(s-1)},\rho)\\ \quad\text{with probability }\min\big(1,\ \exp[H(\theta^{(s-1)},\rho)-H(\theta^*,\rho^*)]\big)\text{ set }\theta^{(s)}=\theta^*,\ \text{else }\theta^{(s)}=\theta^{(s-1)}\end{array}} \tag{3.16}$$
Cost: $L$ gradient evaluations per iteration (the gradient at the end of one step is reused at the start of the next). The gradient comes from LA-14's rules by hand (the project does this, and checks each by finite differences) or, in practice, from automatic differentiation.

**Tuning.**
- *Step size $\varepsilon$.* Controls acceptance through $\Delta H=O(\varepsilon^2)$. Stan adapts it during warm-up to a target acceptance of about 0.8 (the high-dimensional optimum is about 0.65; Beskos et al. 2013). Jitter $\varepsilon$ by $\pm20\%$ to avoid resonance with periodic orbits.
- *Trajectory length $T=L\varepsilon$.* Too short and HMC is a random walk; too long and the trajectory turns around (on a Gaussian, the orbit is periodic) and comes back. The **No-U-Turn Sampler** (Hoffman and Gelman, 2014) chooses $L$ adaptively: it extends the trajectory forward and backward in time by doubling until the two ends start moving toward each other, $(\theta^+-\theta^-)\cdot\rho^\pm<0$, then samples a point from the trajectory in a way that preserves detailed balance. "Your turn" 2 asks you to implement the U-turn criterion in its simplest form.
- *Mass matrix $M$.* With $M\approx\mathrm{Cov}_\pi(\theta)^{-1}$, the dynamics in the coordinates $M^{1/2}\theta$ see an isotropic target. This is LA-13's preconditioning: the leapfrog stability limit is set by the smallest sd and the number of steps to cross by the largest, so the cost scales with the condition number of the posterior covariance, and $M$ reduces it to 1. Stan estimates a diagonal $M$ during warm-up.

**Why HMC beats random walks: scaling with dimension.** For a product target in $d$ dimensions, RWM needs step $s\propto d^{-1/2}$ (§3.2) and its diffusive motion needs $O(1/s^2)=O(d)$ steps to move an $O(1)$ distance, each costing one density evaluation. HMC's energy error over a trajectory sums $d$ independent per-coordinate errors of size $\varepsilon^2$; for $\Delta H$ to stay $O(1)$ in the exponent, with errors partly cancelling as in a random sum, $\varepsilon\propto d^{-1/4}$. A trajectory of fixed length $T$ then needs $L=T/\varepsilon\propto d^{1/4}$ steps, and because the motion is ballistic (distance $\propto T$, not $\sqrt T$) one trajectory gives a nearly independent draw. Neal (2011, §5.4.4) states the result: **the number of gradient evaluations per effectively independent draw grows like $d^{1/4}$ for HMC against $d$ for RWM** (each evaluation costs $O(d)$, so total work is $d^{5/4}$ against $d^2$). The project's scaling plot measures $n_{\rm eff}$ per evaluation on $\mathcal N(0,I_d)$: RWM $0.143, 0.030, 0.0065$ at $d=2,10,50$ (close to $0.3/d$), HMC $0.55, 0.34, 0.18$; the HMC advantage grows from $3.8\times$ to $28\times$.

**Divergences.** When the trajectory enters a region where the curvature of $U$ is much larger than $1/\varepsilon^2$, the leapfrog goes unstable and $\Delta H$ explodes. Stan flags a **divergence** when $\Delta H$ exceeds 1000 (the project uses the same threshold). Divergences mean there is a region the sampler cannot enter, and estimates are biased away from it; lowering $\varepsilon$ usually moves the problem deeper rather than removing it. The usual culprit is the neck of a hierarchical funnel (§3.6).

### 3.5 Diagnostics

**Effective sample size.** Let $\theta^{(1)},\dots,\theta^{(S)}$ be a stationary chain with variance $\sigma^2$ and autocorrelations $\rho_k=\mathrm{Corr}(\theta^{(t)},\theta^{(t+k)})$.

▶ **Try it first.** Expand $\mathrm{Var}(\bar\theta)=S^{-2}\sum_s\sum_t\mathrm{Cov}(\theta^{(s)},\theta^{(t)})$ and count how many pairs are at lag $k$.

**Derivation.**
1. $\mathrm{Var}(\bar\theta)=\frac{\sigma^2}{S^2}\sum_{s,t}\rho_{|s-t|}$.
2. There are $S$ pairs at lag $0$ and $2(S-k)$ pairs at lag $k\ge1$, so
$$\mathrm{Var}(\bar\theta)=\frac{\sigma^2}S\Big[1+2\sum_{k=1}^{S-1}\Big(1-\frac kS\Big)\rho_k\Big]. \tag{3.17}$$
3. If $\sum_k|\rho_k|<\infty$, the factor $(1-k/S)\to1$ for each fixed $k$ as $S\to\infty$ and the tail is negligible, so $\mathrm{Var}(\bar\theta)\approx\sigma^2/n_{\rm eff}$ with
$$n_{\rm eff}=\frac{S}{1+2\sum_{k\ge1}\rho_k}. \tag{3.18}$$ ∎

■ **What just happened.** $n_{\rm eff}$ is the number of independent draws that would give the same variance of the mean. It is per quantity: the chain for $\theta_1$ and the chain for $\theta_1^2$ or a tail probability have different $n_{\rm eff}$. The **Monte Carlo standard error** is $\mathrm{MCSE}=\hat\sigma/\sqrt{n_{\rm eff}}$; report posterior means to the digits it supports. $\qquad(3.19)$

**Estimating the sum.** Sample autocorrelations at large lags are pure noise, so the sum must be truncated. Geyer's *initial positive sequence* estimator sums the pairs $\rho_{2m}+\rho_{2m+1}$ (which are positive for reversible chains) until the first negative pair. The project's `ess` does this with FFT autocorrelations. Antithetic chains (negative $\rho_1$, which HMC can produce) can have $n_{\rm eff}>S$.

**Potential scale reduction $\hat R$ (Gelman and Rubin, 1992).** Run $m$ chains of length $n$ from starting points dispersed relative to the posterior. Let $\bar\theta_{\cdot j}$ be chain $j$'s mean, $s_j^2$ its sample variance, and
$$W=\frac1m\sum_js_j^2,\qquad B=\frac n{m-1}\sum_j(\bar\theta_{\cdot j}-\bar\theta_{\cdot\cdot})^2,\qquad\widehat{\mathrm{var}}^+=\frac{n-1}nW+\frac Bn,\qquad\hat R=\sqrt{\widehat{\mathrm{var}}^+/W}. \tag{3.20}$$

▶ **Try it first.** Suppose all chains are already independent draws from $\pi$ with variance $\sigma^2$. Compute $\mathbb E W$ and $\mathbb E[B/n]$.

**Derivation.**
1. Under that supposition $\mathbb E s_j^2=\sigma^2$, so $\mathbb EW=\sigma^2$.
2. The chain means are independent with variance $\sigma^2/n$, and $B/n$ is their sample variance, so $\mathbb E[B/n]=\sigma^2/n$.
3. Then $\mathbb E\,\widehat{\mathrm{var}}^+=\frac{n-1}n\sigma^2+\frac{\sigma^2}n=\sigma^2$. Both $W$ and $\widehat{\mathrm{var}}^+$ estimate $\sigma^2$, and $\hat R\to1$.
4. Now suppose the chains have not mixed. Each chain has explored only part of the posterior, so $W$ underestimates $\sigma^2$; the chain means differ by more than sampling noise, so $B/n$ is inflated, and with overdispersed starts $\widehat{\mathrm{var}}^+$ overestimates $\sigma^2$. Both push $\hat R$ above 1. ∎

■ **What just happened.** $\hat R$ compares between-chain and within-chain variability, which is analysis of variance (ST-3's $F$ statistic) applied to chains. $\hat R$ near 1 says the chains are indistinguishable. Modern practice: **split** each chain in half (so a chain that drifts slowly shows up as two disagreeing halves), require $\hat R<1.01$, and apply it to rank-normalized draws as well (Vehtari et al. 2021) so heavy tails do not hide disagreement. The project's split-$\hat R$: four dispersed RWM chains on the Gaussian give $1.004$; two chains each stuck in one mode of $\frac12\mathcal N(-4,1)+\frac12\mathcal N(4,1)$ give $3.65$.

**Trace plots.** Plot every chain against iteration on one axis. Mixed chains are overlapping stationary fuzzy bands; a trend, a level shift, or one chain in a different band means trouble. One chain cannot judge its own convergence: a chain stuck in one mode looks perfectly stationary.

**Warm-up, burn-in, thinning.**
- *Burn-in* (discarding the start) matters for finite $S$ when you start far from the mass, though not asymptotically. *Warm-up* in Stan also adapts $\varepsilon$ and $M$, and those draws must be discarded because the kernel was changing.
- *Thinning* (keeping every $k$-th draw) never reduces Monte Carlo variance for a fixed run: the thinned average uses less information (Geyer 1992). It is justified only by memory or the cost of downstream computation per draw.

**The honest statement.** No diagnostic proves convergence. $\hat R\approx1$ and large $n_{\rm eff}$ mean only that the chains you ran agree with each other. If every chain misses the same mode, or none enters the neck of a funnel, every diagnostic is clean and every answer is wrong. Diagnostics can disprove convergence; they cannot prove it. Simulation-based calibration (draw $\theta$ from the prior, simulate $y$, fit, check that the true $\theta$'s rank among the posterior draws is uniform) is the strongest check available for a sampler-plus-model combination.

### 3.6 Hierarchical models

**The eight schools.** Rubin (1981) analysed an experiment on SAT coaching in eight high schools; BDA3 §5.5 uses it throughout. School $j$ reports an estimated effect $y_j$ with known standard error $\sigma_j$:

| school | A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|---|
| $y_j$ | 28 | 8 | −3 | 7 | −1 | 1 | 18 | 12 |
| $\sigma_j$ | 15 | 10 | 16 | 11 | 9 | 11 | 10 | 18 |

The model:
$$y_j\mid\theta_j\sim\mathcal N(\theta_j,\sigma_j^2),\qquad\theta_j\mid\mu,\tau\sim\mathcal N(\mu,\tau^2),\qquad\mu\sim\mathcal N(0,5^2),\qquad\tau\sim\text{half-Cauchy}(0,5). \tag{3.21}$$
(The hyperpriors are those of the Stan case study; BDA3 uses a flat prior on $\mu$.) **Exchangeability** is the assumption behind the second line: before the data you have no reason to rank the schools, so your joint prior for $(\theta_1,\dots,\theta_8)$ is permutation-invariant, and de Finetti's theorem makes such a prior (for an extendable sequence) a mixture of iid models, which is the hierarchical form. If you knew something distinguishing a school, you would add a covariate.

**Three fits.**
- *No pooling* ($\tau=\infty$): each $\theta_j$ estimated from its own $y_j$; school A's effect is $28\pm15$.
- *Complete pooling* ($\tau=0$): one common effect, the precision-weighted mean $\bar y_w=\sum y_j\sigma_j^{-2}/\sum\sigma_j^{-2}=7.69$ with sd $4.07$.
- *Partial pooling*: $\tau$ is estimated from the data, and each $\theta_j$ lands in between.

**Shrinkage.**

▶ **Try it first.** Given $\mu$ and $\tau$, school $j$ has a Normal prior $\mathcal N(\mu,\tau^2)$ and one Normal observation $y_j$ with variance $\sigma_j^2$. You did this update in ST-5 (3.10).

**Derivation.**
1. Conditional on $(\mu,\tau)$ the schools are independent, and $\theta_j$'s posterior combines precisions: $\mathrm{prec}=1/\sigma_j^2+1/\tau^2$, mean $=(y_j/\sigma_j^2+\mu/\tau^2)/\mathrm{prec}$. So
$$\theta_j\mid\mu,\tau,y\sim\mathcal N\big((1-B_j)y_j+B_j\mu,\ (1-B_j)\sigma_j^2\big),\qquad B_j=\frac{\sigma_j^2}{\sigma_j^2+\tau^2}. \tag{3.22}$$
2. Integrating $\theta_j$ out of the model, $y_j\mid\mu,\tau\sim\mathcal N(\mu,\sigma_j^2+\tau^2)$ independently. With $\mu\sim\mathcal N(0,25)$, $\mu$'s posterior given $\tau$ is again precision-weighted:
$$\mu\mid\tau,y\sim\mathcal N\Big(\hat\mu(\tau),\ V_\mu(\tau)\Big),\qquad\hat\mu(\tau)=\frac{\sum_jw_jy_j}{\sum_jw_j+1/25},\quad V_\mu=\frac1{\sum_jw_j+1/25},\quad w_j=\frac1{\sigma_j^2+\tau^2}. \tag{3.23}$$
3. Averaging (3.22) over (3.23),
$$\mathbb E[\theta_j\mid\tau,y]=(1-B_j)\,y_j+B_j\,\hat\mu(\tau). \tag{3.24}$$
4. Finally $p(\tau\mid y)\propto p(\tau)\,\mathcal N\big(y;\,0,\ \mathrm{diag}(\sigma_j^2+\tau^2)+25\,\mathbf 1\mathbf 1^\top\big)$, since $y$ is jointly Normal given $\tau$ after integrating out $\mu$ and $\theta$ (F5 Theorem 3.22). This is a one-dimensional density you can put on a grid, and every posterior expectation is a grid integral of (3.24). The project does this to get exact answers to check the samplers against: $\mathbb E[\tau\mid y]=3.60$. ∎

■ **What just happened.** Each school's estimate is pulled from its own data toward the population mean by the **shrinkage factor** $B_j$, the fraction of school $j$'s total variance that is measurement noise. Noisy schools (large $\sigma_j$) are shrunk more. The amount of pooling is not chosen by you: $\tau$ is estimated from how much the $y_j$ vary beyond what their $\sigma_j$ explain. Here they barely do (the $y_j$ spread is consistent with noise alone), so $p(\tau\mid y)$ concentrates on small values and the posterior means $\mathbb E[\theta_j\mid y]$ are all between $3.6$ and $6.3$.

**Why partial pooling wins.** Each $y_j$ is unbiased for $\theta_j$ but noisy; $\bar y_w$ is biased for every $\theta_j$ but stable. Shrinkage trades a little bias for a large reduction in variance, the bias–variance trade-off of ST-3's ridge regression; the hierarchical model chooses the trade-off from the data ("borrowing strength"). Stein (1956) and James and Stein (1961) proved the frequentist version: for $k\ge3$ Normal means with known variance, the estimator that shrinks every $y_j$ toward a common point by a data-chosen factor has smaller total squared-error risk than $y$ itself, *for every* $\theta$. The hierarchical posterior mean is a smoothed version of that estimator (Efron and Morris, 1975).

**The funnel.** Look at the joint posterior of $(\theta_1,\log\tau)$. Given $\tau$, (3.22) says $\theta_j$'s spread around $\mu$ is about $\min(\tau,\sigma_j)$. As $\tau\to0$ the $\theta_j$ are forced together: the posterior is wide in $\theta$ at large $\tau$ and needle-thin at small $\tau$, and the prior term $\prod_j\tau^{-1}e^{-(\theta_j-\mu)^2/2\tau^2}$ pulls mass into the neck. The curvature of $-\log\pi$ in the $\theta$ directions is about $1/\tau^2$, so a step size that works at $\tau=5$ is unstable at $\tau=0.2$, and a step size small enough for the neck crawls at the top. RWM rejects everything in the neck. Gibbs fails the way §3.3 predicts: given $\tau$ small, $\theta_j$ cannot move; given $\theta_j$ close together, $\tau$ cannot grow; the chain freezes. HMC with fixed $\varepsilon$ diverges in the neck. With eight schools and noisy data the funnel is mild, and it is still enough. The project runs four HMC chains on the **centered** parametrization (3.21) with $\varepsilon=0.15$, $L=25$: one chain enters the neck at $\tau\approx0.16$ and never leaves (2500 divergences in 2500 iterations), the others record 13, 4 and 22 divergences, $n_{\rm eff}(\log\tau)=8$ and $\hat R=1.91$. Funnels are the normal case for hierarchical models with few groups or weak data, not an exotic pathology.

**The non-centered parametrization.** Write
$$\theta_j=\mu+\tau\eta_j,\qquad\eta_j\sim\mathcal N(0,1),\qquad y_j\sim\mathcal N(\mu+\tau\eta_j,\sigma_j^2), \tag{3.26}$$
and sample $(\mu,\log\tau,\eta)$.

▶ **Try it first.** Write the prior density of $(\theta,\tau)$ in the centered form and of $(\eta,\tau)$ in the non-centered form. Which one factors?

**Derivation.**
1. Centered prior: $p(\tau)\prod_j\tau^{-1}\varphi((\theta_j-\mu)/\tau)$. The $\tau^{-J}$ factor and the $\tau$ inside $\varphi$ couple $\tau$ to every $\theta_j$: the conditional sd of $\theta_j$ is $\tau$. That coupling is the funnel.
2. The map $\eta_j=(\theta_j-\mu)/\tau$ has Jacobian $\tau^{-J}$, which cancels the $\tau^{-J}$: the non-centered prior is $p(\tau)\prod_j\varphi(\eta_j)$. **In the prior, $\eta$ and $\tau$ are independent**, and the geometry of each $\eta_j$ is a unit Gaussian at every $\tau$.
3. The likelihood $\prod_j\mathcal N(y_j;\mu+\tau\eta_j,\sigma_j^2)$ couples them again, but weakly when the data are weak ($\sigma_j$ large relative to $\tau$): the posterior is then close to the prior, which is easy. ∎

■ **What just happened.** The two parametrizations describe the same model; they differ only in the coordinates the sampler sees. Non-centered is better when the data are weak relative to the prior (few observations per group, large $\sigma_j$); centered is better when each group has a lot of data, because then the likelihood pins $\theta_j$ and the non-centered coordinates become correlated instead. On the eight schools the non-centered HMC run has zero divergences, $n_{\rm eff}(\log\tau)=5969$, $\hat R=1.000$, and $\mathbb E[\tau\mid y]=3.62$ against the grid's $3.60$.

**Hierarchical logistic regression.** Groups $j=1,\dots,J$, binary outcomes with random intercepts:
$$P(y_i=1)=\sigma(\alpha_{j[i]}+\beta x_i),\qquad\alpha_j=\mu+\tau\eta_j,\quad\eta_j\sim\mathcal N(0,1). \tag{3.27}$$
No conditional is conjugate, so HMC on $(\mu,\log\tau,\beta,\eta)$ is the tool, and intercepts are shrunk toward $\mu$, more for small groups. This is the **generalized linear mixed model**: the random-effects models ST-4 pointed to are these hierarchical models, which frequentist software fits by maximizing the likelihood with the random effects integrated out (Laplace approximation) rather than by sampling. The project fits (3.27) to 8 synthetic groups of 15 to 59 observations: posterior intercepts sit between the raw group logits and $\mu$ (for the 15-observation group the raw logit $-1.01$ becomes $-0.63$, near $\hat\mu=-0.55$).

### 3.7 Model evaluation with draws

**Posterior predictive checks** (ST-5 §3.10) need no formula when you have draws: for each draw $\theta^{(s)}$ simulate $y^{\rm rep(s)}\sim p(y\mid\theta^{(s)})$, compute a test statistic $T$, and report $p_B=\frac1S\#\{s:T(y^{\rm rep(s)})\ge T(y)\}$. $\qquad(3.28)$ For the eight schools with $T=\max_jy_j-\min_jy_j=31$, the non-centered fit gives $p_B=0.72$: the model reproduces the observed spread.

**Predictive accuracy.** The quantity we want is the expected log predictive density for new data, $\mathrm{elpd}=\sum_i\mathbb E_{\tilde y_i}\log p(\tilde y_i\mid y)$. Evaluating on the training data overestimates it. Two estimates:

- **WAIC** (Watanabe, 2010). The log pointwise predictive density is computed from draws, and the overfit is estimated by the posterior variance of the log-likelihood:
$$\mathrm{lppd}=\sum_i\log\Big(\frac1S\sum_sp(y_i\mid\theta^{(s)})\Big),\qquad p_{\rm WAIC}=\sum_i\mathrm{Var}_s\big[\log p(y_i\mid\theta^{(s)})\big],\qquad\widehat{\mathrm{elpd}}_{\rm WAIC}=\mathrm{lppd}-p_{\rm WAIC}. \tag{3.29}$$
$p_{\rm WAIC}$ is an effective number of parameters: in a regular model with weak priors it is close to the parameter count (AIC's penalty), and it is smaller when a prior or pooling constrains parameters.

- **Leave-one-out cross-validation by importance sampling.** The LOO predictive is $p(y_i\mid y_{-i})=\int p(y_i\mid\theta)p(\theta\mid y_{-i})d\theta$. Since $p(\theta\mid y_{-i})\propto p(\theta\mid y)/p(y_i\mid\theta)$ (for conditionally independent data), full-posterior draws reweighted by $r_s=1/p(y_i\mid\theta^{(s)})$ are importance draws from $p(\theta\mid y_{-i})$:
$$p(y_i\mid y_{-i})\approx\frac{\sum_sr_s\,p(y_i\mid\theta^{(s)})}{\sum_sr_s},\qquad\widehat{\mathrm{elpd}}_{\rm LOO}=\sum_i\log p(y_i\mid y_{-i}). \tag{3.30}$$
The raw ratios can have infinite variance when observation $i$ is influential (§3.1's warning). **PSIS** (Vehtari, Gelman and Gabry, 2017) fits a generalized Pareto distribution to the largest $M\approx\min(S/5,3\sqrt S)$ ratios and replaces them by the fitted distribution's quantiles, which stabilizes the estimate. The fitted shape $\hat k$ is a diagnostic: $\hat k<0.5$ fine, $0.5$ to $0.7$ usable, above $0.7$ the estimate is unreliable and that observation needs an exact refit. $\qquad(3.31)$

The project compares the three eight-schools fits: $\widehat{\mathrm{elpd}}_{\rm LOO}=-37.1$ (no pooling, max $\hat k=1.23$), $-30.6$ (complete pooling), $-30.7$ (hierarchical, max $\hat k=0.71$); WAIC gives $-34.1$, $-30.6$, $-30.7$ with $p_{\rm WAIC}=4.02$, $0.67$, $0.83$. No pooling loses, and its $\hat k>1$ is honest: under no pooling, leaving out school $j$ leaves $\theta_j$ with only its flat prior, so the full posterior is a terrible proposal for that LOO posterior. Complete and partial pooling tie, within any reasonable standard error, because the data are consistent with $\tau$ near 0; the hierarchical model is preferred on grounds of not committing to $\tau=0$.

**Why these instead of Bayes factors.** ST-5 Example 4.10 showed the evidence moves by large factors when a prior's scale changes even when the posterior does not, and it is undefined for improper priors. WAIC and LOO depend on the posterior, so they are stable under prior changes that do not change the fit, they answer the usual question (which model predicts new data better after fitting), and they come with per-observation diagnostics. Use them to compare and posterior predictive checks to criticize.

### 3.8 Variational inference

**The idea.** Choose a family $\mathcal Q=\{q_\phi\}$ of tractable densities and solve
$$q^*=\arg\min_{q\in\mathcal Q}\mathrm{KL}\big(q\,\Vert\,p(\cdot\mid y)\big),\qquad\mathrm{KL}(q\Vert p)=\mathbb E_q\Big[\log\frac{q(\theta)}{p(\theta\mid y)}\Big]$$
(F5 Definition 3.44, with densities). This direction puts the expectation under $q$, which we can sample or integrate; the reverse, $\mathrm{KL}(p\Vert q)$, needs expectations under the posterior, which is the problem we cannot solve.

**Theorem 3.4 (ELBO).** For any density $q$ with $q\ll p(\cdot\mid y)$,
$$\log p(y)=\underbrace{\mathbb E_q[\log p(\theta,y)]-\mathbb E_q[\log q(\theta)]}_{\mathrm{ELBO}(q)}+\mathrm{KL}\big(q\Vert p(\cdot\mid y)\big). \tag{3.32}$$
Hence $\mathrm{ELBO}(q)\le\log p(y)$, and maximizing the ELBO over $q$ is the same as minimizing the KL.

▶ **Try it first.** Write $p(\theta\mid y)=p(\theta,y)/p(y)$ inside the KL and pull out the constant.

**Derivation.**
1. $\mathrm{KL}(q\Vert p(\cdot\mid y))=\mathbb E_q[\log q(\theta)]-\mathbb E_q[\log p(\theta,y)]+\log p(y)$, since $\log p(y)$ does not depend on $\theta$ and $\mathbb E_q$ of a constant is the constant.
2. Rearrange: $\log p(y)=\mathbb E_q[\log p(\theta,y)]-\mathbb E_q[\log q]+\mathrm{KL}$.
3. KL $\ge0$ (Gibbs' inequality, F5 §3.13), so the ELBO is a lower bound, and since $\log p(y)$ does not depend on $q$, ELBO and $-$KL differ by a constant. ∎

■ **What just happened.** The intractable normalizing constant was moved to the left side, where it is a constant. The ELBO needs only the *unnormalized* joint $p(\theta,y)$. In SM-7's language, $-\mathrm{ELBO}$ is a variational free energy (energy minus entropy) and $-\log p(y)$ the true free energy; the KL is the gap. ML-10 maximizes the same bound over model parameters as well.

**Mean-field and CAVI.** The mean-field family is $q(\theta)=\prod_jq_j(\theta_j)$, independent blocks with no other restriction.

**Theorem 3.5 (Coordinate-ascent update).** Holding $q_{-j}$ fixed, the ELBO is maximized over $q_j$ by
$$q_j^*(\theta_j)\propto\exp\big\{\mathbb E_{-j}[\log p(\theta,y)]\big\}, \tag{3.33}$$
where $\mathbb E_{-j}$ averages over $\prod_{i\ne j}q_i$.

▶ **Try it first.** Write the ELBO as a function of $q_j$ alone, collecting everything else into constants, and recognize a negative KL.

**Derivation.**
1. With $q=q_j\prod_{i\ne j}q_i$, $\mathbb E_q[\log p(\theta,y)]=\mathbb E_{q_j}\big[\mathbb E_{-j}\log p(\theta,y)\big]$, and $\mathbb E_q[\log q]=\mathbb E_{q_j}[\log q_j]+\sum_{i\ne j}\mathbb E_{q_i}[\log q_i]$.
2. Define $\tilde q_j(\theta_j)=\exp\{\mathbb E_{-j}\log p(\theta,y)\}/Z_j$. Then $\mathbb E_{-j}\log p=\log\tilde q_j+\log Z_j$ and
$$\mathrm{ELBO}=\mathbb E_{q_j}[\log\tilde q_j]-\mathbb E_{q_j}[\log q_j]+\text{const}=-\mathrm{KL}(q_j\Vert\tilde q_j)+\text{const}.$$
3. The KL is zero exactly when $q_j=\tilde q_j$. ∎

■ **What just happened.** Each factor is set to "the conditional, with the other blocks' uncertainty averaged in the exponent". Compare Gibbs (3.6): Gibbs *samples* $\theta_j$ from $\pi(\theta_j\mid\theta_{-j})\propto\exp\log p(\theta,y)$ at the current $\theta_{-j}$; CAVI averages the log over $q_{-j}$. This is the mean-field equation of the Ising model (SM-6: each spin sees the *average* field of its neighbours). The ELBO increases at every update, so CAVI converges to a local optimum. When the full conditionals are exponential-family (the conjugate pieces of §3.3), $q_j^*$ is in the same family and the update is a formula.

**CAVI for the Normal model.** Prior NIG$(\mu_0,\kappa_0,\alpha_0,\beta_0)$, mean-field $q(\mu,\sigma^2)=q(\mu)q(\sigma^2)$, and the log joint from (3.8):
$$\log p=-\big(\alpha_0+\tfrac{n+1}2+1\big)\log\sigma^2-\frac1{\sigma^2}\Big(\beta_0+\tfrac12\Big[\sum_i(y_i-\mu)^2+\kappa_0(\mu-\mu_0)^2\Big]\Big)+\text{const}.$$
1. *Update $q(\mu)$.* Take $\mathbb E_{\sigma^2}$ (only $1/\sigma^2$ appears with $\mu$): $\mathbb E_{\sigma^2}\log p=-\frac{\mathbb E[1/\sigma^2]}2\big[\sum(y_i-\mu)^2+\kappa_0(\mu-\mu_0)^2\big]+$const $=-\frac{\mathbb E[1/\sigma^2]\kappa_n}{2}(\mu-\mu_n)^2+$const. So
$$q(\mu)=\mathcal N\Big(\mu_n,\ v\Big),\qquad v=\frac1{\kappa_n\,\mathbb E_q[1/\sigma^2]}. \tag{3.34}$$
2. *Update $q(\sigma^2)$.* Take $\mathbb E_\mu$: the bracket's expectation under $\mathcal N(\mu_n,v)$ is $\sum(y_i-\mu_n)^2+nv+\kappa_0[(\mu_n-\mu_0)^2+v]$ (using $\mathbb E(c-\mu)^2=(c-\mu_n)^2+v$). The result is an IG kernel:
$$q(\sigma^2)=\mathrm{IG}(A,B),\quad A=\alpha_0+\tfrac{n+1}2,\quad B=\beta_0+\tfrac12\Big[\sum_i(y_i-\mu_n)^2+\kappa_0(\mu_n-\mu_0)^2+\kappa_nv\Big],\quad\mathbb E_q[1/\sigma^2]=\frac AB. \tag{3.35}$$
3. *Fixed point.* The first two terms in the bracket equal $2(\beta_n-\beta_0)$ (the identity of Example 4.2), so $B=\beta_n+\frac12\kappa_nv$ and $v=B/(\kappa_nA)$. Substituting: $B=\beta_n+B/(2A)$, so $B=\beta_n\frac{2A}{2A-1}$. With $A=\alpha_n+\frac12$: $B=\beta_n\frac{2\alpha_n+1}{2\alpha_n}$ and $\mathbb E_q[1/\sigma^2]=A/B=\alpha_n/\beta_n$, the exact posterior mean of the precision.
4. *Compare with the truth.* $q(\mu)$ has the exact mean $\mu_n$ and variance $v=\beta_n/(\kappa_n\alpha_n)$. The exact marginal is $t_{2\alpha_n}(\mu_n,\beta_n/(\alpha_n\kappa_n))$ (ST-5 (3.14)) with variance $\frac{\beta_n}{\alpha_n\kappa_n}\cdot\frac{2\alpha_n}{2\alpha_n-2}=\frac{\beta_n}{\kappa_n(\alpha_n-1)}$. The ratio is
$$\frac{\mathrm{Var}_q(\mu)}{\mathrm{Var}(\mu\mid y)}=\frac{\alpha_n-1}{\alpha_n}<1. \tag{3.36}$$ ∎

■ **What just happened.** Mean-field VI got the mean exactly and the variance too small. The missing variance is the part that comes from the *dependence* between $\mu$ and $\sigma^2$ in the true posterior: when $\sigma^2$ is large, $\mu$ is uncertain, and averaging over that gives the $t$'s heavier tails. A factorized $q$ cannot represent the dependence, so it replaces "Normal with random variance" by "Normal with average precision". With $\alpha_n=5$ (Example 4.7) VI reports 80% of the true variance.

**Gaussian VI and the reparametrization trick.** For non-conjugate models take $q_\phi=\mathcal N(m,LL^\top)$ with $L$ lower triangular and maximize the ELBO by gradient ascent. The entropy is $\log|\det L|+$const. The expectation term's gradient with respect to $\phi$ is not an expectation under a fixed distribution, because $q$ depends on $\phi$. Write $\theta=m+Lz$ with $z\sim\mathcal N(0,I)$; now the randomness is in $z$, which does not depend on $\phi$, and
$$\nabla_m\mathrm{ELBO}=\mathbb E_z[\nabla\log p(m+Lz,y)],\qquad\nabla_L\mathrm{ELBO}=\mathbb E_z[\nabla\log p(m+Lz,y)\,z^\top]+L^{-\top}\ (\text{lower triangle}), \tag{3.37}$$
estimated by averaging over a few draws of $z$ (chain rule, LA-14). This is the **reparametrization trick** that ML-10 uses to train VAEs. It needs the same gradient as HMC.

**VI's bias.** Minimizing $\mathrm{KL}(q\Vert\pi)=\mathbb E_q[\log q-\log\pi]$ makes $q$ pay infinitely wherever $q>0$ and $\pi\approx0$, and nothing for missing parts of $\pi$. So $q$ fits inside the high-density region: it is **mode-seeking** and underestimates variance; with a multimodal posterior it typically locks onto one mode. On the banana density (project Part 8), Gaussian VI centres at $\theta_2=0.22$ with variances $0.23,0.25$, against HMC's mean $1.05$ and variances $1.04,2.54$ (exact: $1$, and $1, 2.25$): it sits in the bend of the banana and ignores both arms. MCMC is slow and exact as $S\to\infty$; VI is fast and biased at any budget. Use VI for huge models, inside training loops (ML-10), or for initialization; use MCMC when you will report uncertainty.

### 3.9 The workflow, and software

Gelman et al., "Bayesian Workflow" (2020), summarized: (1) prior predictive check; (2) fit, first on simulated data with known parameters; (3) check the computation ($\hat R$, $n_{\rm eff}$, divergences, traces), fixing geometry by reparametrizing before tuning; (4) posterior predictive checks; (5) compare with LOO; (6) expand or simplify and repeat. Their "folk theorem": a model that is hard to compute is often wrong statistically.

Software: **Stan** (own language, NUTS, autodiff), **PyMC** (Python, NUTS and VI), **NumPyro** (JAX, NUTS on GPUs). All implement what this module builds, with careful adaptation and diagnostics.

## 4. Worked examples

Every number below is asserted or printed by `project.py` (the part is in parentheses).

### Example 4.1: The MH ratio, symmetric and not (Part 2)

*Gaussian random walk on $\pi=\mathcal N(0,1)$.* The proposal $\theta'=\theta+sz$ is symmetric, so (3.5) is $\min(1,e^{-(\theta'^2-\theta^2)/2})$. From $\theta=1$ to $\theta'=1.8$: $\theta'^2-\theta^2=2.24$, $\alpha=e^{-1.12}=0.3263$. The reverse move, $1.8\to1$, has $\alpha=1$. Check (3.4): $\pi(1)\,q\,(0.3263)=\varphi(1)(0.3263)q=0.2420\times0.3263\,q=0.0790\,q$ and $\pi(1.8)\,q\cdot1=0.0790\,q$.

*Multiplicative proposal on a positive parameter.* Target Gamma(3,1), $\pi(\theta)\propto\theta^2e^{-\theta}$. Propose $\theta'=\theta e^{z}$, $z\sim\mathcal N(0,0.8^2)$, i.e. $\log\theta'\sim\mathcal N(\log\theta,0.64)$. The density of $\theta'$ is $q(\theta'\mid\theta)=\varphi_{0.8}(\log\theta'-\log\theta)/\theta'$; the $\varphi$ part is symmetric, so $q(\theta\mid\theta')/q(\theta'\mid\theta)=\theta'/\theta$. From $\theta=2$ to $\theta'=3$: $\pi(3)/\pi(2)=9e^{-3}/(4e^{-2})=2.25e^{-1}=0.8277$; Hastings factor $3/2$; $r=1.2416$, so $\alpha=1$. Without the factor you would accept with probability $0.8277$ and, as the project shows, sample $\theta^2e^{-\theta}\cdot\theta^{-1}$, i.e. Gamma(2,1): sample mean $1.98$ instead of $3.01$.

### Example 4.2: Gibbs for the Normal model (Part 3)

Data $y=(4.2,5.1,3.8,6.0,5.4,4.7)$, $n=6$, $\bar y=4.8667$, $S=\sum(y_i-\bar y)^2=3.2333$. Prior NIG$(0,1,2,2)$. Then $\kappa_n=7$, $\mu_n=6(4.8667)/7=4.1714$, $\alpha_n=5$, $\beta_n=2+\frac12[3.2333+\frac67(4.8667)^2]=2+\frac12[3.2333+20.3010]=13.7671$.

The sampler alternates (3.9) and (3.10). At $\sigma^2=1$: $\mu\sim\mathcal N(4.1714,\ 1/7=0.1429)$. At $\mu=\mu_n$: $\sigma^2\sim\mathrm{IG}(2+3.5,\ 2+\frac12[\sum(y_i-4.1714)^2+(4.1714)^2])=\mathrm{IG}(5.5,\ 13.7671)$. The second parameter equals $\beta_n$ exactly, because $\sum(y_i-\mu_n)^2+\kappa_0(\mu_n-\mu_0)^2=S+n(\bar y-\mu_n)^2+\kappa_0(\mu_n-\mu_0)^2$ and, with $\bar y-\mu_n=\kappa_0(\bar y-\mu_0)/\kappa_n$ and $\mu_n-\mu_0=n(\bar y-\mu_0)/\kappa_n$, the last two terms sum to $\frac{n\kappa_0}{\kappa_n}(\bar y-\mu_0)^2$. Notice the conditional shape $\alpha_0+\frac{n+1}2=5.5$ differs from the marginal shape $\alpha_n=5$: conditioning on $\mu$ adds information. After 40,000 sweeps (thinned by 10 after 1000 warm-up), KS tests against the exact $t_{10}(4.1714,\ 13.7671/35)$ marginal of $\mu$ and $\mathrm{IG}(5,13.7671)$ marginal of $\sigma^2$ give $p=0.060$ and $p=0.937$.

### Example 4.3: Albert–Chib conditionals (Part 3)

Suppose at the current $\beta$ observation $i$ has $x_i^\top\beta=0.5$. If $y_i=1$, (3.11) draws $z_i$ from $\mathcal N(0.5,1)$ truncated to $(0,\infty)$, whose mean is $0.5+\varphi(0.5)/\Phi(0.5)=0.5+0.3521/0.6915=1.0092$. If $y_i=0$, truncated to $(-\infty,0]$: $0.5-\varphi(0.5)/(1-\Phi(0.5))=0.5-0.3521/0.3085=-0.6411$. (Truncated-Normal mean: $\mathbb E[Z\mid Z>-m]=\varphi(m)/\Phi(m)$ for $Z\sim\mathcal N(0,1)$, by integrating $z\varphi(z)$, whose antiderivative is $-\varphi(z)$.) A 1 pushes the latent value up past the threshold; a 0 pushes it below. Then (3.12) regresses $z$ on $X$. On 400 simulated observations with true $\beta=(-0.3,0.8)$ and prior $\mathcal N(0,100I)$, 3500 kept draws give posterior mean $(-0.262,\ 0.645)$, sd $(0.070,\ 0.075)$, against ST-4's probit MLE $(-0.262,\ 0.640)$: with a weak prior and $n=400$ the posterior mean and the MLE agree to within a tenth of a posterior sd, as Bernstein–von Mises (ST-5 §3.8) predicts.

### Example 4.4: One leapfrog step by hand (Part 4)

Target $\mathcal N(0,1)$: $U=\theta^2/2$, $\nabla\log\pi=-\theta$, $M=1$. Start $(\theta,\rho)=(1,0.5)$, $\varepsilon=0.1$. $H_0=\frac12+\frac18=0.625$.
1. Half kick: $\rho_{1/2}=0.5-0.05(1)=0.45$.
2. Drift: $\theta'=1+0.1(0.45)=1.045$.
3. Half kick: $\rho'=0.45-0.05(1.045)=0.39775$.

$H_1=\frac12(1.045)^2+\frac12(0.39775)^2=0.5460125+0.0791025=0.6251150$, so $\Delta H=1.1503\times10^{-4}$ and the acceptance probability is $e^{-\Delta H}=0.99988$. By (3.15), $\Delta H=\frac{\varepsilon^2}8(\theta'^2-\theta^2)=0.00125\times0.092025=1.1503\times10^{-4}$: the energy error is exactly the change in the shadow term. The exact flow would follow the circle $\theta^2+\rho^2=1.25$; leapfrog follows the ellipse $\rho^2+0.9975\,\theta^2=$ const.

### Example 4.5: $n_{\rm eff}$ for $\rho_k=0.9^k$ (Part 5)

By (3.18), $1+2\sum_{k\ge1}0.9^k=1+2\cdot\frac{0.9}{0.1}=19$, so $n_{\rm eff}=S/19$. In general for $\rho_k=\phi^k$, $1+\frac{2\phi}{1-\phi}=\frac{1+\phi}{1-\phi}$. An AR(1) chain $x_t=0.9x_{t-1}+\sqrt{0.19}\,e_t$ has exactly these autocorrelations. With $S=200{,}000$: theory $10{,}526$; Geyer's estimator gives $10{,}123$ ($4\%$ low, within its sampling error). The MCSE of the mean is $\sqrt{19}=4.4$ times the naive $\sigma/\sqrt S$.

### Example 4.6: Eight-schools shrinkage at $\tau=10$ (Part 6)

Take $\tau=10$, a value in the upper part of the posterior. Weights $w_j=1/(\sigma_j^2+100)$: $1/325,1/200,1/356,1/221,1/181,1/221,1/200,1/424$, with $\sum w_j=0.03282$ and $\sum w_jy_j=0.26670$. By (3.23), $\hat\mu(10)=0.26670/(0.03282+0.04)=3.663$. (With a flat prior on $\mu$ it would be $0.26670/0.03282=8.13$: the $\mathcal N(0,5^2)$ prior on $\mu$ is worth more than all eight schools at this $\tau$, which is worth knowing about the case-study prior.) Shrinkage factors $B_j=\sigma_j^2/(\sigma_j^2+100)$: $0.692, 0.500, 0.719, 0.548, 0.448, 0.548, 0.500, 0.764$. By (3.24), school A: $(1-0.692)(28)+0.692(3.663)=11.15$; school E (the most precise, $B=0.448$): $0.552(-1)+0.448(3.663)=1.09$; school C: $1.79$; school G: $10.83$. At the posterior mean $\tau=3.60$ the factors are $0.86$ to $0.96$, and the full-posterior means $\mathbb E[\theta_j\mid y]$, integrated over $p(\tau\mid y)$ on a grid, are $6.21, 4.94, 3.93, 4.76, 3.62, 4.04, 6.30, 4.85$: almost complete pooling, and school A's 28 is pulled to about 6.

### Example 4.7: CAVI on Example 4.2's data (Part 8)

Same data and prior: $\mu_n=4.1714$, $\kappa_n=7$, $\alpha_n=5$, $\beta_n=13.7671$, and $A=\alpha_0+\frac{n+1}2=5.5$ is fixed. Start with $\mathbb E[1/\sigma^2]=1$.
1. (3.34): $v=1/7=0.14286$. (3.35): $B=\beta_n+\frac12\kappa_nv=13.7671+0.5=14.2671$.
2. $\mathbb E[1/\sigma^2]=5.5/14.2671$, $v=B/(\kappa_nA)=14.2671/38.5=0.37058$, $B=13.7671+3.5(0.37058)=15.0642$.
3. $v=0.39128$, $B=15.1366$; and so on to the fixed point $B=13.7671\times\frac{11}{10}=15.1439$, $v=0.39335$.

So $q(\mu)=\mathcal N(4.1714,\ 0.39335)$ and $q(\sigma^2)=\mathrm{IG}(5.5,\ 15.1439)$, with $\mathbb E_q[1/\sigma^2]=0.3632=\alpha_n/\beta_n$. The exact $\mathrm{Var}(\mu\mid y)=\beta_n/(\kappa_n(\alpha_n-1))=13.7671/28=0.49168$. Ratio $0.8000=(\alpha_n-1)/\alpha_n$, as (3.36) predicts. The iteration contracts: the map $B\mapsto\beta_n+B/(2A)$ has slope $1/11$.

## 5. Why it matters and where it connects

- **SM-6, SM-7.** Metropolis is SM-6's Ising simulation; SM-7 reads a posterior as a Boltzmann distribution with energy $-\log\pi$, the evidence as a partition function, $-$ELBO (3.32) as variational free energy, (3.33) as the mean-field equation. Langevin dynamics, HMC with $L=1$, is the sampling step of diffusion models.
- **CM-8.** HMC is the module where Liouville's theorem (Theorem 3.4 there) and the symplecticity of Störmer–Verlet (§3.11 there) become load-bearing in a non-physics computation: they are exactly steps 1–2 of Theorem 3.3's derivation here.
- **ML-10.** The ELBO (3.32) and the reparametrization gradient (3.37), with a neural network producing $(m,L)$ from each data point ("amortized" VI), are the variational autoencoder.
- **ST-8.** State-space models have one latent state per time step; Gibbs samplers draw the whole state path given parameters (forward filtering, backward sampling) and the parameters given the path, the data-augmentation pattern of §3.3.
- **LA-13, LA-14.** The mass matrix is a preconditioner; HMC's cost is governed by the posterior covariance's condition number.

## 6. Exercises

See [exercises.md](exercises.md): 16 problems. A1–A5 are drills (the MH ratio, a leapfrog step, $n_{\rm eff}$ for AR(1) chains, Gibbs conditionals). B1–B9 are the module: MH for a non-symmetric proposal, Gibbs on a correlated Gaussian, HMC acceptance with a mass matrix, the non-centered Jacobian, shrinkage, $\hat R$ by hand, a funnel experiment in code, and the CAVI variance ratio. C1–C2 are challenges: Metropolis–Hastings with a deterministic involution (a general version of Theorem 3.3) and parallel tempering in code. Essential: A1, A2, A3, B1, B2, B5, B6, B9. B8 and C2 modify `project.py`.

## 7. Coding project

`project.py` (about 20 seconds) builds every sampler from numpy and scipy, each tested first on a target with a known answer. (1) The curse of dimension (3.3). (2) MH: Example 4.1, the missing-Hastings-factor bug, RWM traces, acceptance and $n_{\rm eff}$ against scale (`rwm_scale.png`). (3) Gibbs for the Normal model (KS tests against NIG), a $\rho=0.99$ Gaussian, regression, Albert–Chib against the probit MLE. (4) HMC: Example 4.4, the shadow Hamiltonian in SymPy, energy error slope, acceptance against $\varepsilon$, leapfrog's Jacobian and reversibility against Euler, banana trajectories, and the scaling plot (`scaling.png`, the signature figure). (5) $n_{\rm eff}$ and $\hat R$. (6) Eight schools: exact grid posterior, centered against non-centered HMC, funnel and shrinkage plots (`eight_schools.png`), a predictive check, WAIC and PSIS-LOO. (7) Hierarchical logistic regression. (8) CAVI and Gaussian VI against HMC (`vi_banana.png`). (9) Exercise numbers.

**Build it yourself first:** write `rwm`, `leapfrog`, `hmc` and `ess` (about 40 lines) from (3.5), (3.14), (3.16) and (3.18) before reading the shared-tools cell; then write the non-centered eight-schools log density and gradient and check it with `check_grad`. Predict why the non-centered run can reach $\tau\approx0$ while the centered one cannot. "Your turn": Ising denoising, the U-turn criterion, parallel tempering.

## 8. Where you'll struggle

- **Detailed balance is one line, and the whole subject.** Theorem 3.1, Theorem 3.2 and Theorem 3.3 are three applications of (3.4). If you cannot reproduce F5's proof that (3.27) implies stationarity, or do not see why the Metropolis ratio is the answer to "make the two flows equal", go back to F5 §3.12 and redo its exercise on the Metropolis chain. Everything else in this module is engineering on top of that line.
- **A sampler that runs is not a sampler that works.** A chain always produces numbers, and the numbers always have a mean. Every bug in this module (a missing Hastings factor, a wrong sign in a gradient, a Jacobian forgotten after a log transform) produces a chain that looks fine and samples the wrong distribution. The defences: test every sampler on a target with a known answer first (that is why the project starts each method on a conjugate model), check gradients by finite differences, run several chains from dispersed starts, and read $\hat R$, $n_{\rm eff}$ and the divergence count every time.
- **HMC genuinely requires CM-8.** If "leapfrog preserves volume" is a phrase rather than a computation (two shear matrices of determinant 1) and "reversible" is not "flip, run, flip is the inverse", Theorem 3.3 will be magic and divergences will be mysterious. Redo CM-8 §3.5 and §3.11 and the project's leapfrog-Jacobian check (Part 4), then Example 4.4 by hand.
- **The log-transform Jacobian.** Sampling $\log\tau$ requires adding $\log\tau$ to the log density (the `+ lt` in the project). Forget it and you silently sample a different prior.
- **The funnel is the normal case.** Every hierarchical model with weak data per group has one. The non-centered parametrization (3.26) is mandatory knowledge, and so is its limitation: with strong data per group, centered is better. Diagnose with divergences and the scatter of $\log\tau$ against one group parameter.
- **$\hat R\approx1$ does not prove convergence; it fails to disprove it.** Four chains that all miss a mode agree perfectly. The project's centered run looks bad because one chain got stuck; with a different seed all four might avoid the neck and report clean diagnostics while being biased against small $\tau$. Divergences are the extra signal HMC gives you; take every one seriously.
- **"Burn-in" and "thinning" are folklore with partial justification** (§3.5): discard adaptive warm-up; thin only for storage.
- **VI is fast and biased; MCMC is slow and exact in the limit.** Knowing which error you can tolerate is the craft. VI's means are often good, its variances too small (3.36), and its answer to a multimodal posterior is one mode.
- **Implement these samplers once from scratch.** Using Stan afterwards is fine; using it before means you cannot tell a modelling problem from a sampling problem when it reports divergences.

## 9. Check yourself

1. Why does $p(y)$ never appear in an MCMC algorithm? Where does it reappear in this module?
2. Derive the MH acceptance probability from detailed balance. What changes for a symmetric proposal?
3. Prove that a Gibbs update is an MH step with acceptance 1.
4. Write the two full conditionals for the Normal model with an NIG prior, and the two for Albert–Chib probit.
5. Why does Gibbs sampling slow down on a correlated posterior? Give the lag-one autocorrelation for a bivariate Normal with correlation $\rho$.
6. State the three properties of Hamiltonian flow and what each contributes to HMC's correctness.
7. Write the leapfrog step and show its Jacobian determinant is 1.
8. Why is HMC's energy error $O(\varepsilon^2)$ and not growing with $L$? What happens above the stability limit?
9. Derive $n_{\rm eff}=S/(1+2\sum\rho_k)$ and compute it for $\rho_k=\phi^k$.
10. Define $W$, $B$ and $\hat R$. Why can four chains with $\hat R=1.00$ still be wrong?
11. Give the shrinkage factor for school $j$ and the non-centered parametrization. Why does the latter remove the funnel?
12. Derive $\log p(y)=\mathrm{ELBO}+\mathrm{KL}$ and the CAVI update. Why does mean-field VI underestimate the variance of $\mu$ in the Normal model, and by what factor?

## 10. Pair with

- **Betancourt, "A Conceptual Introduction to Hamiltonian Monte Carlo" (arXiv 1701.02434). Required.** Typical sets and why random walks fail in high dimension, drawn far better than here; read §§2–5 after §3.4.
- **Gelman et al., *Bayesian Data Analysis*, 3rd ed., ch. 10–13 and ch. 5.** The reference for MH, Gibbs, HMC, diagnostics ($\hat R$, $n_{\rm eff}$ as Stan computes them) and the eight schools (§5.5). Denser, with more applied examples.
- **McElreath, *Statistical Rethinking*, 2nd ed., ch. 9 and 13.** The best intuition-first treatment of MCMC and multilevel models, with divergent transitions in practice.
- **Neal, "MCMC using Hamiltonian dynamics" (Handbook of MCMC, 2011).** The original careful account of HMC, its scaling (§5.4.4) and its tuning; the source of the $d^{1/4}$ result.
- **Blei, Kucukelbir and McAuliffe, "Variational Inference: A Review for Statisticians" (JASA 2017).** CAVI, the exponential-family case, and stochastic VI, with the clearest derivations in print.
- **Gelman and Hill, *Data Analysis Using Regression and Multilevel/Hierarchical Models*, ch. 11–13.** Hierarchical models as applied statisticians use them, and the link to mixed models.
- **Stan documentation (User's Guide, "Reparameterization" and "Efficiency tuning").** For practice once you have built the samplers here.
