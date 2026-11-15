# ST-6 Exercises

Conventions as in the notes: $\mathrm{IG}(a,b)$ has mean $b/(a-1)$; $H=U+\frac12\rho^\top M^{-1}\rho$ with $U=-\log\pi$; $n_{\rm eff}$ is (3.18) and $\hat R$ is (3.20) (unsplit unless stated). Equation numbers refer to `notes.md`. Numerical answers are printed by Part 9 of `project.py`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Curse of dimension.) Target $\mathcal N(0,I_d)$, rejection sampling from $\mathcal N(0,1.2^2I_d)$. (a) Show the smallest valid envelope constant is $M=1.2^d$ and the acceptance rate is $1.2^{-d}$. (b) Evaluate at $d=10$ and $d=50$. (c) In one sentence, why does random-walk Metropolis not suffer this exponential decay?

**A2.** (The MH ratio.) Target $\mathcal N(0,1)$, Gaussian random-walk proposal. (a) Compute the acceptance probability of the move $0.5\to-1.2$ and of the reverse move. (b) Verify detailed balance (3.4) for this pair of points, up to the common proposal density.

**A3.** (Leapfrog.) Target $\mathcal N(0,1)$, $M=1$, $\varepsilon=0.2$, start $(\theta,\rho)=(0,1)$. (a) Take two leapfrog steps (3.14) by hand. (b) Compute $\Delta H$ and the acceptance probability. (c) Check your $\Delta H$ against the shadow Hamiltonian (3.15).

**A4.** ($n_{\rm eff}$.) For a chain with $\rho_k=\phi^k$, compute $n_{\rm eff}/S$ for $\phi=0.5$, $\phi=0.99$ and $\phi=-0.5$. What does the last answer mean, and which sampler in this module can produce negative lag-one autocorrelation?

**A5.** (Reading off full conditionals.) The model $y_i\mid\lambda\sim\mathrm{Poisson}(\lambda)$ for $i=1,\dots,n$, $\lambda\mid b\sim\mathrm{Gamma}(a,b)$ (shape–rate), $b\sim\mathrm{Gamma}(c,d)$, with $a,c,d$ known. Write the joint density and derive both full conditionals for a Gibbs sampler.

## B. Core (do all)

**B1.** (The Hastings factor.) For a positive parameter, propose $\theta'=\theta e^z$ with $z\sim\mathcal N(0,s^2)$. (a) Derive $q(\theta'\mid\theta)$ and show $q(\theta\mid\theta')/q(\theta'\mid\theta)=\theta'/\theta$. (b) Suppose you forget the factor and accept with $\min(1,\pi(\theta')/\pi(\theta))$. Show the chain is a correct MH sampler for a different target $\tilde\pi$, and find $\tilde\pi$. (c) For $\pi=$ Gamma(3,1), what does the wrong chain sample, and what mean does it report? (Compare Example 4.1.) (d) Show that the same factor arises if you instead run a symmetric random walk on $\phi=\log\theta$ with target density for $\phi$ obtained by change of variables.

**B2.** (Gibbs on a ridge.) Target $(\theta_1,\theta_2)\sim\mathcal N\big(0,\begin{pmatrix}1&\rho\\\rho&1\end{pmatrix}\big)$. (a) Derive the full conditionals (F5 Theorem 3.24). (b) Show that the $\theta_1$-subsequence of a systematic-scan Gibbs sampler is an AR(1) process with coefficient $\rho^2$ and stationary variance 1. (c) Compute $n_{\rm eff}/S$ for $\rho=0.9$ and $\rho=0.99$. (d) Give a linear reparametrization in which Gibbs produces independent draws, and say why it is not available in general.

**B3.** (Mass matrix and stability.) Target $\mathcal N(0,\Sigma)$ with $\Sigma$ positive definite. (a) With $M=I$, show that leapfrog decouples into one-dimensional oscillators along the eigenvectors of $\Sigma$, with frequencies $\omega_i=\lambda_i^{-1/2}$ ($\lambda_i$ the eigenvalues). (b) For one oscillator, write one leapfrog step as a $2\times2$ matrix and show it is stable (bounded orbits) iff $\varepsilon\omega<2$. Conclude that $\varepsilon<2\sqrt{\lambda_{\min}}$, and evaluate it for $\rho=0.9$ in B2's matrix. (c) Show that with $M=\Sigma^{-1}$ every frequency is 1. Why does this make the number of steps per trajectory independent of the condition number of $\Sigma$?

**B4.** (Centered vs non-centered.) For one school, $y\mid\theta\sim\mathcal N(\theta,\sigma^2)$, $\theta=\mu+\tau\eta$, $\eta\sim\mathcal N(0,1)$, with $\mu,\tau$ fixed. (a) Show the Jacobian of $\eta\mapsto\theta$ cancels the $\tau^{-1}$ in the centered prior density. (b) Find the posterior of $\eta$ given $(\mu,\tau,y)$. (c) Show its sd is about 1 when $\sigma\gg\tau$ and about $\sigma/\tau$ when $\sigma\ll\tau$. (d) Use (c) to explain when the non-centered parametrization has its own funnel, and therefore when the centered one is better.

**B5.** (Shrinkage.) A group has $\sigma_j=4$. (a) Compute $B_j$ for $\tau=1,4,16$ and interpret each. (b) With $\hat\mu=10$ and $y_j=30$, compute $\mathbb E[\theta_j\mid\tau,y]$ at each $\tau$. (c) Show from (3.22) that the conditional posterior variance $(1-B_j)\sigma_j^2$ is less than both $\sigma_j^2$ and $\tau^2$.

**B6.** ($\hat R$ by hand.) Two chains of length 4: $(1.0,1.2,0.8,1.1)$ and $(3.0,3.2,2.8,3.1)$. (a) Compute $W$, $B$, $\widehat{\mathrm{var}}^+$ and $\hat R$ from (3.20). (b) Now suppose a single chain consists of the first sequence followed by the second. Explain why split-$\hat R$ catches this and unsplit $\hat R$ on one chain cannot be computed at all.

**B7.** (Thinning.) For an AR(1) chain with $\phi=0.9$, compare $n_{\rm eff}/S$ using all $S$ draws with $n_{\rm eff}/S$ (still divided by the original $S$) after keeping every 10th draw. Prove in general that for $\rho_k=\phi^k$ with $0<\phi<1$, thinning by $k\ge2$ reduces $n_{\rm eff}$. When is thinning justified?

**B8.** (Code: the funnel and the step size.) Modify Part 6 of `project.py` to run the centered model with $\varepsilon=0.05$ and $\varepsilon=0.3$, keeping $\varepsilon L\approx3.75$, four chains of 1000 iterations. Report divergences per chain, $n_{\rm eff}(\log\tau)$, $\hat R$, the smallest $\tau$ visited and $\mathbb E[\tau\mid y]$, and compare with the exact $3.60$ and with the non-centered run. Does a smaller step fix the funnel? What does it cost?

**B9.** (VI underestimates variance.) (a) Re-derive the CAVI fixed point (3.34)–(3.35) and the ratio (3.36). (b) Evaluate the ratio for $\alpha_n=2,5,50$. (c) Explain in words, using Eve's law (F5 Theorem 3.28), where the missing variance went. (d) Does $\mathrm{Var}_q(\sigma^2)$ also underestimate $\mathrm{Var}(\sigma^2\mid y)$? Compute both for Example 4.7.

## C. Challenge (do at least one)

**C1.** (Involutive MH.) Let $p$ be a density on $\mathbb R^m$ and $G:\mathbb R^m\to\mathbb R^m$ a smooth involution ($G\circ G=I$). (a) Prove that the deterministic proposal $z\to G(z)$ accepted with probability $\min\big(1,\ p(G(z))\,|\det DG(z)|/p(z)\big)$ satisfies detailed balance. (b) Show $|\det DG(G(z))|=1/|\det DG(z)|$, which is what makes (a) work when $G$ does not preserve volume. (c) Recover Theorem 3.3 and the HMC acceptance rule for leapfrog as special cases. (d) Why must HMC be combined with momentum resampling to be useful?

**C2.** (Code: parallel tempering.) Implement "Your turn" 3: chains at temperatures $T\in\{1,2,4,8\}$ targeting $\pi^{1/T}$ for $\pi=\frac12\mathcal N(-4,1)+\frac12\mathcal N(4,1)$, a random-walk step in each chain with scale $1.5\sqrt T$, and one proposed swap of a random adjacent pair per iteration. (a) Derive the swap acceptance probability from (3.5) applied to the product target $\prod_i\pi(x_i)^{1/T_i}$. (b) Run two replicas started in opposite modes; report $\hat R$ of the $T=1$ chain against the single-temperature value in Part 5, $P(x>0)$, and the swap rate. (c) Why does the $T=1$ chain still have exactly the right stationary distribution?
