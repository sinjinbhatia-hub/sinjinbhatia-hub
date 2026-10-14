# F5 Exercises

Conventions as in the notes: Gamma and Exponential use the **rate** $\lambda$; Markov chains use **row** vectors, $\boldsymbol\pi_{t+1} = \boldsymbol\pi_tP$; $\log$ is natural unless written $\log_2$. Every numerical answer is checked in `project.py`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** $X \in \{1, 2\}$ and $Y \in \{1, 2, 3\}$ have joint PMF

| $p(x,y)$ | $y=1$ | $y=2$ | $y=3$ |
|---|---|---|---|
| $x=1$ | 0.10 | 0.20 | 0.10 |
| $x=2$ | 0.20 | 0.20 | 0.20 |

(a) Find both marginals and both conditional PMFs $p_{Y|X}(\cdot\mid x)$. Are $X$ and $Y$ independent?
(b) Compute $\mathbb{E}[X]$, $\mathbb{E}[Y]$, $\mathbb{E}[XY]$ and $\mathrm{Cov}(X, Y)$.
(c) Compute $\mathbb{E}[Y\mid X]$ and $\mathrm{Var}(Y\mid X = x)$ for both $x$. Verify Eve's law. Which of its two terms is zero, and what does that say?

**A2.** $(X, Y)$ has covariance matrix $\Sigma = \begin{pmatrix} 4 & -3 \\ -3 & 9\end{pmatrix}$. Using only bilinearity (Theorem 3.14) or (3.14), compute $\mathrm{Var}(X + Y)$, $\mathrm{Var}(X - Y)$, $\mathrm{Var}(2X - 3Y)$, $\mathrm{Cov}(2X + 1,\, 3Y - X)$ and $\rho$. Then find the eigenvalues of $\Sigma$ and confirm Theorem 3.19.

**A3.** (a) $X_1, \dots, X_n$ iid $\mathrm{Exponential}(\lambda)$. Show $\min_i X_i \sim \mathrm{Exponential}(n\lambda)$. What is the mean for $n = 5$, $\lambda = 2$?
(b) For $X_1, X_2$ iid $\mathrm{Exponential}(1)$, find the density of $\max(X_1, X_2)$ and its mean.
(c) $X \sim \mathrm{Pois}(2)$ and $Y \sim \mathrm{Pois}(3)$ are independent. Find $P(X + Y = 4)$.

**A4.** $\mathbf{X} = (X_1, X_2) \sim \mathcal{N}\!\left(\begin{pmatrix}1\\2\end{pmatrix}, \begin{pmatrix}4 & -2\\-2 & 4\end{pmatrix}\right)$.
(a) Find $\rho$, the eigenvalues and eigenvectors of $\Sigma$, and describe the level curves.
(b) Find the distribution of $X_1 \mid X_2 = 4$.
(c) Find the joint distribution of $(X_1 + X_2,\; X_1 - X_2)$. Are these two independent? Which theorem?

**A5.** (a) Compute the entropy in bits of a fair six-sided die and of $p = (\tfrac12, \tfrac14, \tfrac14)$.
(b) With $u$ uniform on three outcomes, compute $D(p\|u)$ and $D(u\|p)$ in bits, and the cross-entropy $H(p, u)$. Show that $D(p\|u) = \log_2 3 - H(p)$, and explain why this gives a second proof of Theorem 3.43(b).

## B. Core (do all)

**B1.** (Box–Muller backwards.) Let $Z_1, Z_2$ be iid $\mathcal{N}(0,1)$ and write $(Z_1, Z_2) = (R\cos\Theta, R\sin\Theta)$ with $R > 0$, $\Theta \in (-\pi, \pi]$.
(a) Use Theorem 3.6 to find the joint density of $(R, \Theta)$. Show $R \perp \Theta$, $\Theta$ is uniform, and $R$ has density $re^{-r^2/2}$ (the Rayleigh distribution). Compute $\mathbb{E}[R]$.
(b) Let $T = Z_1/Z_2$. By changing variables from $(Z_1, Z_2)$ to $(T, W)$ with $W = Z_2$, show $f_T(t) = \dfrac{1}{\pi(1 + t^2)}$ (Cauchy). Find the quartiles of $T$.
(c) Why does part (b) mean the CLT cannot apply to averages of ratios of normals?

**B2.** (Best linear predictor, vector form.) Let $\mathbf{X}_1 \in \mathbb{R}^p$ and $\mathbf{X}_2 \in \mathbb{R}^q$ have *any* joint distribution with finite second moments, means $\boldsymbol\mu_1, \boldsymbol\mu_2$, and covariance blocks $\Sigma_{11}, \Sigma_{12}, \Sigma_{22}$ with $\Sigma_{22}$ invertible.
(a) Show that $\mathbf{a}^*$, $B^*$ minimize $\mathbb{E}\|\mathbf{X}_1 - \mathbf{a} - B\mathbf{X}_2\|^2$ if and only if the residual $\mathbf{R} = \mathbf{X}_1 - \mathbf{a}^* - B^*\mathbf{X}_2$ has mean zero and $\mathrm{Cov}(\mathbf{R}, \mathbf{X}_2) = 0$. Deduce $B^* = \Sigma_{12}\Sigma_{22}^{-1}$ and $\mathbf{a}^* = \boldsymbol\mu_1 - B^*\boldsymbol\mu_2$.
(b) Show $\mathrm{Cov}(\mathbf{R}) = \Sigma_{11} - \Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21}$.
(c) Compare with (3.19). What is special about the Gaussian case?

**B3.** (Adam and Eve on a mixture.) $Z \sim \mathrm{Bernoulli}(0.3)$. Given $Z = 0$, $Y \sim \mathcal{N}(0, 1)$; given $Z = 1$, $Y \sim \mathcal{N}(4, 2^2)$.
(a) Find $\mathbb{E}[Y]$ by Adam's law.
(b) Find $\mathrm{Var}(Y)$ by Eve's law. Name the "within" and "between" terms.
(c) Find $P(Y > 3)$.
(d) Find $P(Z = 1\mid Y = 3)$ by Bayes' rule for densities.
(e) Is $Y$ normal? Justify in one sentence.

**B4.** (Beta–Binomial.) A coin has unknown heads-probability $p \sim \mathrm{Beta}(1, 1)$. You observe 7 heads in 10 tosses.
(a) Give the posterior distribution, its mean and its standard deviation.
(b) Prove, for general $a, b, n, k$, that the posterior mean $\frac{a+k}{a+b+n}$ is a weighted average of the prior mean and $k/n$. Give the weights, and their values here.
(c) Compute $P(p > \tfrac12 \mid \text{data})$ *exactly*, without integrating, by recognizing $\mathrm{Beta}(8, 4)$ as an order statistic of Uniforms (Theorem 3.9).

**B5.** (Gaussian marginals are not enough.) Let $X \sim \mathcal{N}(0, 1)$ and, independently, $S = \pm1$ with probability $\tfrac12$ each. Put $Y = SX$.
(a) Show $Y \sim \mathcal{N}(0,1)$.
(b) Show $\mathrm{Cov}(X, Y) = 0$.
(c) Show $X$ and $Y$ are not independent.
(d) Show $(X, Y)$ is not bivariate normal. (Find a linear combination that is not normal.) Which hypothesis of Theorem 3.23 fails?

**B6.** (Chebyshev versus the CLT.) Toss a fair coin $n = 100$ times; let $\bar X$ be the fraction of heads.
(a) Bound $P(|\bar X - \tfrac12| \ge 0.1)$ by Chebyshev.
(b) Approximate it by the CLT, without and with the continuity correction. Compare with the exact Binomial value.
(c) How large must $n$ be for Chebyshev to *guarantee* $P(|\bar X - \tfrac12| \ge 0.1) \le 0.05$? How large does the CLT say is enough?
(d) For $\mathrm{Bernoulli}(0.1)$ summands, compute $\rho_3/\sigma^3$ and the Berry–Esseen bound (3.25) at $n = 100$. Compare with the actual sup error, which `project.py` prints.

**B7.** (A law of large numbers without independence.)
(a) Let $X_1, X_2, \dots$ be **uncorrelated** with common mean $\mu$ and $\mathrm{Var}(X_i) \le \sigma^2$. Prove $\bar X_n \xrightarrow{p} \mu$.
(b) Take the two-state chain on $\{0, 1\}$ with $P = \begin{pmatrix}1-a & a\\ b & 1-b\end{pmatrix}$, $a = 0.2$, $b = 0.3$, started from its stationary distribution. Show $\mathrm{Cov}(X_0, X_k) = \lambda_2^k\,\pi_0\pi_1$ with $\lambda_2 = 1 - a - b$.
(c) Show $\mathrm{Var}(\bar X_n) \le \dfrac{\pi_0\pi_1}{n}\cdot\dfrac{1 + |\lambda_2|}{1 - |\lambda_2|}$ for the average of $X_0, \dots, X_{n-1}$, and conclude the weak law holds for this dependent sequence. Evaluate the constant. Whose argument is this, against whom?

**B8.** (Periodicity.) Consider the simple random walk on the path graph $1 - 2 - 3 - 4$.
(a) Write $P$ and find $\boldsymbol\pi$ from the degrees.
(b) Find the eigenvalues of $P$. Which one is the signature of period 2?
(c) Starting at vertex 1, show $\|\boldsymbol\pi_t - \boldsymbol\pi\|_1 \ge 1$ for every $t$. Why doesn't this contradict Theorem 3.40?
(d) The lazy walk is $\tilde P = \tfrac12(I + P)$. Show it has the same $\boldsymbol\pi$, find its eigenvalues, and give its convergence rate.

**B9.** (Conditioning and entropy.) $X, Y \in \{0, 1\}$ with $p(0,0) = 0.9$, $p(0,1) = 0$, $p(1,0) = 0.05$, $p(1,1) = 0.05$.
(a) Compute $H(Y)$, $H(Y\mid X = 0)$, $H(Y\mid X = 1)$, $H(Y\mid X)$ in bits.
(b) Observe that $H(Y\mid X = 1) > H(Y)$. Prove nonetheless that $H(Y\mid X) \le H(Y)$ always.
(c) Compute $I(X;Y)$ and $H(X, Y)$, and check the chain rule.

**B10.** (Metropolis; needs code.) Target $\pi_k \propto k$ on $\{1, \dots, 5\}$. From state $k$, propose $k \pm 1$ with probability $\tfrac12$ each. Reject proposals that leave $\{1, \dots, 5\}$. Accept a proposal $j$ with probability $\min(1, \pi_j/\pi_k)$; otherwise stay.
(a) Prove that, for *any* target $\boldsymbol\pi > 0$ and any symmetric proposal $Q$ ($Q_{kj} = Q_{jk}$), the Metropolis chain satisfies detailed balance (3.27) with respect to $\boldsymbol\pi$.
(b) Write out the $5\times5$ transition matrix and verify $\boldsymbol\pi P = \boldsymbol\pi$ directly for one column.
(c) In `project.py`, `metropolis_discrete` implements this chain. Change the target to $\pi_k \propto k^2$ on $\{1, \dots, 10\}$ and assert that the visit frequencies after $4\times10^5$ steps are within $0.01$ of the normalized target. Where is the normalizing constant used?

## C. Challenge (do at least one)

**C1.** (KL between multivariate Gaussians.) Derive
$$D\big(\mathcal{N}(\boldsymbol\mu_1, \Sigma_1)\,\|\,\mathcal{N}(\boldsymbol\mu_2, \Sigma_2)\big) = \tfrac12\Big[\mathrm{tr}(\Sigma_2^{-1}\Sigma_1) + (\boldsymbol\mu_2 - \boldsymbol\mu_1)^\top\Sigma_2^{-1}(\boldsymbol\mu_2 - \boldsymbol\mu_1) - k + \ln\frac{\det\Sigma_2}{\det\Sigma_1}\Big].$$
(You will need $\mathbb{E}[\mathbf{V}^\top A\mathbf{V}] = \mathrm{tr}(A\,\mathrm{Cov}(\mathbf{V})) + \mathbb{E}[\mathbf{V}]^\top A\,\mathbb{E}[\mathbf{V}]$; prove it.) Check that it reduces to (3.32) for $k = 1$. Evaluate it for $\boldsymbol\mu_1 = \mathbf{0}$, $\Sigma_1 = \begin{pmatrix}2&1\\1&2\end{pmatrix}$, $\boldsymbol\mu_2 = (1, 0)^\top$, $\Sigma_2 = I$.

**C2.** (Mutual information of a Gaussian pair.) $(X, Y)$ is bivariate normal with unit variances and correlation $\rho$, $|\rho| < 1$. Show $I(X; Y) = -\tfrac12\ln(1 - \rho^2)$, using $I = h(X) + h(Y) - h(X,Y)$ (prove that identity for densities first). Evaluate at $\rho = 0.5$ in nats and bits. What happens as $\rho \to \pm1$, and why is that reasonable? Then do "Your turn" item 3 in `project.py`.

**C3.** (Maximum entropy with constraints.)
(a) Use a Lagrange multiplier to show that the uniform distribution maximizes $H(p)$ over PMFs on $m$ points. Use concavity to show the critical point is the maximum.
(b) (Jaynes' die.) Find the maximum-entropy PMF on $\{1, \dots, 6\}$ with mean $4.5$. Show it has the form $p_i \propto e^{-\beta i}$, find $\beta$ numerically, and give $p$ and $H(p)$ in bits.
(c) Prove your answer to (b) is the global maximum with the KL argument of Theorem 3.46. Relate $\beta$ to SM-4.
