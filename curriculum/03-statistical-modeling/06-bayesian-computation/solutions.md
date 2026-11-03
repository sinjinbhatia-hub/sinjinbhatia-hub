# ST-6 Solutions

Numbers are printed by Part 9 of `project.py` unless another part is named.

## A. Warm-up

**A1.** (a) The ratio of densities is $\pi(\theta)/g(\theta)=\dfrac{(2\pi)^{-d/2}e^{-|\theta|^2/2}}{(2\pi\cdot1.44)^{-d/2}e^{-|\theta|^2/2.88}}=1.2^d\exp\{-|\theta|^2(\tfrac12-\tfrac1{2.88})\}$. The exponent is $\le0$ with equality at $\theta=0$, so $\sup\pi/g=1.2^d$ and the smallest valid $M$ is $1.2^d$. The acceptance probability is $\mathbb E_g[\pi/(Mg)]=\frac1M\int\pi=1.2^{-d}$.
(b) $1.2^{-10}=0.1615$; $1.2^{-50}=1.1\times10^{-4}$.
(c) RWM never compares the whole target with a fixed global proposal; each proposal is a small local perturbation of a point that is already in the high-probability region, so the acceptance ratio compares two nearby points and stays $O(1)$ once the step is scaled like $d^{-1/2}$.

**A2.** (a) The proposal is symmetric, so $\alpha=\min(1,\pi(\theta')/\pi(\theta))=\min\big(1,e^{-(\theta'^2-\theta^2)/2}\big)$. For $0.5\to-1.2$: $\theta'^2-\theta^2=1.44-0.25=1.19$, $\alpha=e^{-0.595}=0.5516$. For $-1.2\to0.5$: the ratio is $e^{0.595}>1$, so $\alpha=1$.
(b) With $q$ the common proposal density $q(-1.2\mid0.5)=q(0.5\mid-1.2)$: $\pi(0.5)\,q\,\alpha(0.5,-1.2)=0.35207\times0.5516\,q=0.19419\,q$, and $\pi(-1.2)\,q\cdot1=\varphi(1.2)\,q=0.19419\,q$. Equal, as (3.4) requires.

**A3.** (a) $\nabla\log\pi(\theta)=-\theta$.
Step 1: $\rho_{1/2}=1-0.1(0)=1$; $\theta=0+0.2(1)=0.2$; $\rho=1-0.1(0.2)=0.98$.
Step 2: $\rho_{1/2}=0.98-0.1(0.2)=0.96$; $\theta=0.2+0.2(0.96)=0.392$; $\rho=0.96-0.1(0.392)=0.9208$.
(b) $H_0=\frac12(0^2+1^2)=0.5$. $H_2=\frac12(0.392^2+0.9208^2)=\frac12(0.153664+0.847873)=0.500768$. $\Delta H=7.6832\times10^{-4}$, acceptance $e^{-\Delta H}=0.99923$.
(c) By (3.15), $\tilde H=\frac12\rho^2+\frac12(1-\frac{\varepsilon^2}4)\theta^2$ is exactly conserved, so $H-\tilde H=\frac{\varepsilon^2}8\theta^2$ and $\Delta H=\frac{\varepsilon^2}8(\theta_2^2-\theta_0^2)=\frac{0.04}8(0.153664-0)=7.6832\times10^{-4}$. It matches to every digit; the script asserts it.

**A4.** $n_{\rm eff}/S=(1-\phi)/(1+\phi)$ (Example 4.5). $\phi=0.5$: $0.5/1.5=0.3333$. $\phi=0.99$: $0.01/1.99=0.0050$, so a chain of 100,000 draws is worth about 500. $\phi=-0.5$: $1.5/0.5=3$: the chain is worth three times as many independent draws. With alternating-sign correlation, consecutive errors in the running mean cancel (antithetic sampling), so the mean has *smaller* variance than an iid average. HMC with long trajectories on near-Gaussian targets can overshoot to the opposite side of the mode and produce negative $\rho_1$; this is why the Geyer estimator must allow $n_{\rm eff}>S$.

**A5.** The joint density is
$$p(y,\lambda,b)\propto\Big[\prod_i\lambda^{y_i}e^{-\lambda}\Big]\cdot\frac{b^a}{\Gamma(a)}\lambda^{a-1}e^{-b\lambda}\cdot b^{c-1}e^{-db}.$$
*Conditional of $\lambda$.* Keep factors with $\lambda$: $\lambda^{\sum y_i+a-1}e^{-(n+b)\lambda}$, a Gamma kernel: $\lambda\mid b,y\sim\mathrm{Gamma}(a+\sum_iy_i,\ n+b)$.
*Conditional of $b$.* Keep factors with $b$: $b^ae^{-b\lambda}\cdot b^{c-1}e^{-db}=b^{a+c-1}e^{-(\lambda+d)b}$: $b\mid\lambda,y\sim\mathrm{Gamma}(a+c,\ \lambda+d)$. It does not depend on $y$ given $\lambda$, because $y$ enters the model only through $\lambda$.

## B. Core

**B1.** (a) $\theta'=\theta e^z$ means $\log\theta'=\log\theta+z\sim\mathcal N(\log\theta,s^2)$. By the change of variables $u=\log\theta'$, $du/d\theta'=1/\theta'$:
$$q(\theta'\mid\theta)=\frac1{\theta'}\,\varphi_s(\log\theta'-\log\theta),\qquad\varphi_s(x)=\frac{e^{-x^2/2s^2}}{\sqrt{2\pi}s}.$$
Since $\varphi_s$ is even, $\varphi_s(\log\theta-\log\theta')=\varphi_s(\log\theta'-\log\theta)$, so $q(\theta\mid\theta')/q(\theta'\mid\theta)=(1/\theta)/(1/\theta')=\theta'/\theta$.
(b) The wrong chain accepts with $\min(1,\pi(\theta')/\pi(\theta))$. A correct MH sampler for $\tilde\pi$ with this $q$ accepts with $\min\big(1,\tilde\pi(\theta')\theta'/(\tilde\pi(\theta)\theta)\big)$. These coincide for all pairs iff $\tilde\pi(\theta')\theta'/(\tilde\pi(\theta)\theta)=\pi(\theta')/\pi(\theta)$, i.e. $\tilde\pi(\theta)\theta\propto\pi(\theta)$, so $\tilde\pi(\theta)\propto\pi(\theta)/\theta$. By Theorem 3.1 the wrong chain is a correct sampler for $\tilde\pi$.
(c) $\pi\propto\theta^2e^{-\theta}$ gives $\tilde\pi\propto\theta e^{-\theta}$, which is Gamma(2,1), mean 2. Part 2 of the project reports $1.983$ (and $3.006$ with the factor).
(d) Let $\phi=\log\theta$. The density of $\phi$ is $\pi_\phi(\phi)=\pi(e^\phi)\,|d\theta/d\phi|=\pi(e^\phi)e^\phi$. A symmetric random walk on $\phi$ accepts with $\pi_\phi(\phi')/\pi_\phi(\phi)=\frac{\pi(\theta')}{\pi(\theta)}\cdot\frac{\theta'}{\theta}$, the same ratio as in (a). The Hastings factor is the Jacobian of the transformation, seen from the other side. Forgetting it is the same bug as forgetting the `+ lt` term when sampling $\log\tau$.

**B2.** (a) By F5 Theorem 3.24 with unit variances, $\theta_1\mid\theta_2\sim\mathcal N(\rho\theta_2,1-\rho^2)$ and $\theta_2\mid\theta_1\sim\mathcal N(\rho\theta_1,1-\rho^2)$.
(b) Write a sweep as $\theta_2^{(t)}=\rho\theta_1^{(t)}+\sqrt{1-\rho^2}\,e_t$, then $\theta_1^{(t+1)}=\rho\theta_2^{(t)}+\sqrt{1-\rho^2}\,e'_t$ with independent standard Normal $e_t,e'_t$. Substituting,
$$\theta_1^{(t+1)}=\rho^2\theta_1^{(t)}+\rho\sqrt{1-\rho^2}\,e_t+\sqrt{1-\rho^2}\,e'_t,$$
an AR(1) with coefficient $\rho^2$ and innovation variance $\rho^2(1-\rho^2)+(1-\rho^2)=(1-\rho^2)(1+\rho^2)=1-\rho^4$. The stationary variance of an AR(1) with coefficient $\phi$ and innovation variance $v$ is $v/(1-\phi^2)=(1-\rho^4)/(1-\rho^4)=1$, the correct marginal variance.
(c) With $\phi=\rho^2$: $n_{\rm eff}/S=(1-\rho^2)/(1+\rho^2)$. $\rho=0.9$: $0.19/1.81=0.1050$. $\rho=0.99$: $0.0199/1.9801=0.01005$ (the project measures $0.0105$).
(d) Put $u=(\theta_1+\theta_2)/\sqrt2$, $v=(\theta_1-\theta_2)/\sqrt2$. This is a rotation onto the eigenvectors of the covariance (LA-8), so $\mathrm{Cov}(u,v)=0$, $\mathrm{Var}\,u=1+\rho$, $\mathrm{Var}\,v=1-\rho$, and for a Gaussian uncorrelated means independent (F5 Theorem 3.23). Then the full conditional of $u$ given $v$ is its marginal, and Gibbs on $(u,v)$ draws independent samples. In general you do not know the posterior covariance before sampling, and for non-Gaussian posteriors the dependence is not removable by a linear map; this is the job adaptive proposals and HMC's mass matrix do approximately.

**B3.** (a) With $M=I$, $\nabla\log\pi=-\Sigma^{-1}\theta$. Write $\Sigma=Q\Lambda Q^\top$ and $\theta=Qx$, $\rho=Qp$. Since $Q$ is orthogonal, $\frac12|\rho|^2=\frac12|p|^2$ and $\theta^\top\Sigma^{-1}\theta=x^\top\Lambda^{-1}x$; each leapfrog sub-step (3.14) commutes with the rotation, so in $(x,p)$ coordinates the steps are $p_i\leftarrow p_i-\frac\varepsilon2\lambda_i^{-1}x_i$, $x_i\leftarrow x_i+\varepsilon p_i$, independently for each $i$. Each pair is a harmonic oscillator $H_i=\frac12p_i^2+\frac12\omega_i^2x_i^2$ with $\omega_i^2=1/\lambda_i$.
(b) For one oscillator, with $w=\varepsilon\omega$, the three sub-steps compose to
$$\begin{pmatrix}x'\\p'\end{pmatrix}=\begin{pmatrix}1&0\\-\frac{\varepsilon\omega^2}2&1\end{pmatrix}\begin{pmatrix}1&\varepsilon\\0&1\end{pmatrix}\begin{pmatrix}1&0\\-\frac{\varepsilon\omega^2}2&1\end{pmatrix}\begin{pmatrix}x\\p\end{pmatrix}=\begin{pmatrix}1-\frac{w^2}2&\varepsilon\\-\varepsilon\omega^2(1-\frac{w^2}4)&1-\frac{w^2}2\end{pmatrix}\begin{pmatrix}x\\p\end{pmatrix}.$$
The determinant is $(1-\frac{w^2}2)^2+w^2(1-\frac{w^2}4)=1$ and the trace is $2-w^2$. The eigenvalues of a $2\times2$ matrix with determinant 1 are $\lambda,1/\lambda$ solving $\lambda^2-(2-w^2)\lambda+1=0$. If $|2-w^2|<2$, i.e. $0<w<2$, they are complex conjugates of modulus 1 and the orbits are bounded ellipses (this is the shadow ellipse of (3.15) for $\omega=1$). If $w>2$, the trace is below $-2$, the eigenvalues are real with one of modulus greater than 1, and the orbit grows exponentially: the energy error explodes, every proposal is rejected. So stability requires $\varepsilon\omega_i<2$ for all $i$, i.e. $\varepsilon<2/\omega_{\max}=2\sqrt{\lambda_{\min}}$. For $\rho=0.9$, $\lambda_{\min}=1-\rho=0.1$, so $\varepsilon<2\sqrt{0.1}=0.632$ (printed in Part 4, where acceptance falls from $0.999$ at $\varepsilon=0.05$ to $0.805$ at $0.45$).
(c) With mass matrix $M$, $\dot\theta=M^{-1}\rho$, $\dot\rho=-\Sigma^{-1}\theta$, so $\ddot\theta=-M^{-1}\Sigma^{-1}\theta$. With $M=\Sigma^{-1}$, $M^{-1}\Sigma^{-1}=\Sigma\Sigma^{-1}=I$: every mode has $\omega=1$. A trajectory must last long enough for the slowest mode to move (time $\sim1/\omega_{\min}$) with a step limited by the fastest ($\varepsilon<2/\omega_{\max}$), so $L\gtrsim\omega_{\max}/\omega_{\min}=\sqrt{\lambda_{\max}/\lambda_{\min}}$, the square root of the condition number. With all frequencies equal to 1 that ratio is 1, so $L$ no longer depends on $\Sigma$'s conditioning. This is LA-13's preconditioning.

**B4.** (a) Centered prior density of $\theta$: $\tau^{-1}\varphi((\theta-\mu)/\tau)$. With $\theta=\mu+\tau\eta$, $d\theta/d\eta=\tau$, so the density of $\eta$ is $\tau^{-1}\varphi(\eta)\cdot\tau=\varphi(\eta)$: free of $\tau$.
(b) The likelihood is $\mathcal N(y;\mu+\tau\eta,\sigma^2)$, i.e. in $\eta$ an observation $(y-\mu)/\tau$ with variance $\sigma^2/\tau^2$. Combining with the $\mathcal N(0,1)$ prior (precisions add, ST-5 (3.10)):
$$\eta\mid\mu,\tau,y\sim\mathcal N\Big(\frac{\tau(y-\mu)}{\sigma^2+\tau^2},\ \frac{\sigma^2}{\sigma^2+\tau^2}\Big).$$
(c) If $\sigma\gg\tau$ the variance is $\approx1$. If $\sigma\ll\tau$ it is $\approx\sigma^2/\tau^2$, sd $\approx\sigma/\tau$.
(d) In case (c)'s second regime the posterior scale of $\eta$ depends on $\tau$ as $1/\tau$: as $\tau$ varies over its posterior, $\eta$'s width changes, which is a funnel in $(\eta,\log\tau)$. Meanwhile in the centered coordinates $\theta\mid\mu,\tau,y$ has variance $(1-B)\sigma^2\approx\sigma^2$, independent of $\tau$: no funnel. So non-centered is right when the data per group are weak ($\sigma_j$ large relative to $\tau$, the eight schools) and centered is right when they are strong (many observations per group). Mixed cases can be parametrized group by group.

**B5.** (a) $B_j=16/(16+\tau^2)$: $\tau=1$: $16/17=0.9412$ (the population is tight, so the group's own noisy data count for 6%); $\tau=4$: $0.5$ (noise and population spread equal, so half-and-half); $\tau=16$: $16/272=0.0588$ (groups genuinely differ, so the group mostly keeps its own estimate).
(b) $\mathbb E[\theta_j\mid\tau,y]=(1-B_j)30+B_j10=30-20B_j$: $11.18$, $20.00$, $28.82$.
(c) $(1-B_j)\sigma_j^2=\frac{\tau^2}{\sigma_j^2+\tau^2}\sigma_j^2=\frac{\sigma_j^2\tau^2}{\sigma_j^2+\tau^2}=\big(\sigma_j^{-2}+\tau^{-2}\big)^{-1}$. Since $\sigma_j^{-2}+\tau^{-2}>\sigma_j^{-2}$, the result is $<\sigma_j^2$, and since it is $>\tau^{-2}$, the result is $<\tau^2$. Two sources of information beat either alone.

**B6.** (a) Chain 1: mean $1.025$, deviations $-0.025,0.175,-0.225,0.075$, squares sum $0.0875$, $s_1^2=0.0875/3=0.02917$. Chain 2 is chain 1 plus 2, so $s_2^2=0.02917$. $W=0.02917$. Grand mean $2.025$; $B=\frac{4}{2-1}[(1.025-2.025)^2+(3.025-2.025)^2]=4\times2=8$. $\widehat{\mathrm{var}}^+=\frac34(0.02917)+\frac84=2.0219$. $\hat R=\sqrt{2.0219/0.02917}=\sqrt{69.3}=8.33$. The chains clearly sit in different places.
(b) A single chain gives $m=1$ and $B$ has denominator $m-1=0$: unsplit $\hat R$ is undefined for one chain, and nothing within the chain is compared. Split-$\hat R$ halves the chain into $(1.0,1.2,0.8,1.1)$ and $(3.0,3.2,2.8,3.1)$, exactly the two chains of (a), and reports $8.33$. Splitting catches chains that drift or jump between regions over time.

**B7.** All draws: $n_{\rm eff}/S=0.1/1.9=0.05263$. Thinned by 10: the kept draws form an AR(1) with $\phi^{10}=0.3487$, so $n_{\rm eff}=\frac S{10}\cdot\frac{1-0.3487}{1+0.3487}$, and $n_{\rm eff}/S=0.04829$. Thinning lost 8% of the information.
*Proof for general $k\ge2$.* We need $\frac1k\frac{1-\phi^k}{1+\phi^k}<\frac{1-\phi}{1+\phi}$. Divide by $1-\phi>0$ and use $\frac{1-\phi^k}{1-\phi}=\sum_{i=0}^{k-1}\phi^i$: the claim is $(1+\phi)\sum_{i=0}^{k-1}\phi^i<k(1+\phi^k)$. The left side is $\sum_{i=0}^{k-1}(\phi^i+\phi^{i+1})=1+2\sum_{i=1}^{k-1}\phi^i+\phi^k$. So the claim is $2\sum_{i=1}^{k-1}\phi^i<(k-1)(1+\phi^k)$. For each $1\le i\le k-1$, $(1-\phi^i)(1-\phi^{k-i})>0$ gives $\phi^i+\phi^{k-i}<1+\phi^k$. Summing over $i=1,\dots,k-1$, the left side is $\sum_i\phi^i+\sum_i\phi^{k-i}=2\sum_{i=1}^{k-1}\phi^i$ and the right side is $(k-1)(1+\phi^k)$. ∎
Thinning is justified only when storing or post-processing each draw is expensive (memory for a high-dimensional model, an expensive function evaluated per draw); then keeping every $k$-th draw of a highly autocorrelated chain loses little.

**B8.** Reference implementation in Part 9 of the script (same $\varepsilon L\approx3.75$, four chains of 1000, 200 warm-up discarded):

| step | divergences per chain | $n_{\rm eff}(\log\tau)$ | $\hat R$ | min $\tau$ | $\mathbb E\tau$ |
|---|---|---|---|---|---|
| $\varepsilon=0.05$, $L=75$ | 0, 2, 1, 1 | 120 | 1.03 | 0.088 | 3.50 |
| $\varepsilon=0.15$, $L=25$ (Part 6, 2500 its) | 13, 2500, 4, 22 | 8 | 1.91 | 0.155 | 2.86 |
| $\varepsilon=0.3$, $L=12$ | 8, 17, 746, 30 | 9 | 2.51 | 0.255 | 2.02 |
| non-centered, $\varepsilon=0.15$ (Part 6) | 0, 0, 0, 0 | 5969 | 1.000 | $\approx0$ | 3.62 |

Exact $\mathbb E[\tau\mid y]=3.60$. A larger step cannot enter the neck: the smallest $\tau$ reached rises to $0.255$, one chain sticks there with hundreds of divergences, and $\mathbb E\tau$ is biased low because the stuck chain sits at small $\tau$ while the others avoid the neck (both errors at once). A smaller step goes deeper (min $\tau=0.088$) with few divergences and a nearly correct mean, so it does help here. The cost: three times as many gradients per iteration as $\varepsilon=0.15$, and still $n_{\rm eff}=120$ from 3200 draws, against 5969 from 8000 for the non-centered run at a third of the cost per iteration. Per gradient evaluation the non-centered sampler is roughly 100 times more efficient. And the neck extends to $\tau\to0$ without limit, so any fixed step eventually fails there; the remaining divergences say so. Reparametrization removes the problem; step size only negotiates with it.

**B9.** (a) Notes §3.8, steps 1–4: $q(\mu)=\mathcal N(\mu_n,v)$ with $v=1/(\kappa_n\mathbb E[1/\sigma^2])$ from collecting the $\mu$-quadratic; $q(\sigma^2)=\mathrm{IG}(A,B)$ with $A=\alpha_n+\frac12$, $B=\beta_n+\frac12\kappa_nv$ from taking $\mathbb E_\mu$ of the bracket; substituting $v=B/(\kappa_nA)$ gives $B=\beta_n+B/(2A)$, so $B=\beta_n\cdot2A/(2A-1)=\beta_n(2\alpha_n+1)/(2\alpha_n)$ and $v=B/(\kappa_nA)=\beta_n/(\kappa_n\alpha_n)$. The exact $t_{2\alpha_n}(\mu_n,\beta_n/(\alpha_n\kappa_n))$ has variance $\frac{\beta_n}{\alpha_n\kappa_n}\cdot\frac{2\alpha_n}{2\alpha_n-2}=\frac{\beta_n}{\kappa_n(\alpha_n-1)}$. Ratio $(\alpha_n-1)/\alpha_n$.
(b) $\alpha_n=2$: $0.5$. $\alpha_n=5$: $0.8$. $\alpha_n=50$: $0.98$. Since $\alpha_n=\alpha_0+n/2$, the ratio is $1-2/(2\alpha_0+n)$: VI's error is $O(1/n)$ and vanishes as data accumulate, but at small $n$ it halves the variance.
(c) Eve's law: $\mathrm{Var}(\mu\mid y)=\mathbb E[\mathrm{Var}(\mu\mid\sigma^2,y)]+\mathrm{Var}(\mathbb E[\mu\mid\sigma^2,y])$. Here $\mathbb E[\mu\mid\sigma^2,y]=\mu_n$ does not depend on $\sigma^2$, so the second term is 0, and the first is $\mathbb E[\sigma^2/\kappa_n]=\beta_n/(\kappa_n(\alpha_n-1))$. VI instead uses $1/(\kappa_n\mathbb E[1/\sigma^2])=\beta_n/(\kappa_n\alpha_n)$. By Jensen (F5 Theorem 3.32, $1/x$ convex), $\mathbb E[\sigma^2]=\mathbb E[1/(1/\sigma^2)]>1/\mathbb E[1/\sigma^2]$. The missing variance is the Jensen gap: mean-field VI averages the *precision* in the exponent instead of averaging the Normal over its random variance.
(d) $\mathrm{Var}_q(\sigma^2)=B^2/((A-1)^2(A-2))=15.1439^2/(4.5^2\times3.5)=3.2358$. Exact: $\beta_n^2/((\alpha_n-1)^2(\alpha_n-2))=13.7671^2/(16\times3)=3.9486$. Yes, about 18% too small; and $\mathbb E_q\sigma^2=3.3653$ against $3.4418$. VI's shape parameter $A=\alpha_n+\frac12$ is too large because it treats $\mu$'s contribution as half an extra observation with averaged uncertainty.

## C. Challenge

**C1.** (a) Let $a(z)=\min\big(1,\ p(G(z))J(z)/p(z)\big)$ with $J(z)=|\det DG(z)|$. The probability of starting in $A$ and moving to $B$ is $P(A\to B)=\int_Ap(z)a(z)\mathbb 1[G(z)\in B]\,dz$. Substitute $z=G(w)$; since $G$ is an involution $G(z)=w$, and $dz=J(w)\,dw$. Then $P(A\to B)=\int_Bp(G(w))a(G(w))J(w)\mathbb 1[G(w)\in A]\,dw$. Detailed balance requires this to equal $P(B\to A)=\int_Bp(w)a(w)\mathbb 1[G(w)\in A]\,dw$, which holds if $p(G(w))a(G(w))J(w)=p(w)a(w)$ for all $w$. Compute, using (b) $J(G(w))=1/J(w)$:
$$p(G(w))J(w)\,a(G(w))=p(G(w))J(w)\min\Big(1,\frac{p(w)J(G(w))}{p(G(w))}\Big)=\min\big(p(G(w))J(w),\ p(w)\big),$$
$$p(w)a(w)=\min\big(p(w),\ p(G(w))J(w)\big).$$
Equal. ∎
(b) Differentiate $G(G(z))=z$: $DG(G(z))\,DG(z)=I$. Take determinants: $\det DG(G(z))\det DG(z)=1$.
(c) HMC with exact flow: $G=F\circ\Phi_T$ is an involution (Theorem 3.3, step 1) with $J=1$ (Liouville), so the rule is $\min(1,p(G(z))/p(z))=\min(1,e^{-\Delta H})$, and $\Delta H=0$ gives acceptance 1. With $L$ leapfrog steps, $G=F\circ\Psi^L$ where $\Psi$ is (3.14): it is an involution because leapfrog is symmetric (flip, $L$ steps, flip is $\Psi^{-L}$), and $J=1$ because each sub-step is a shear. The rule is $\min(1,e^{-\Delta H})$ with the leapfrog's nonzero $\Delta H$: box (3.16).
(d) A deterministic involution alone visits only $z$ and $G(z)$: the chain alternates between two points and is not irreducible. Resampling $\rho\sim\mathcal N(0,M)$ is a Gibbs step that preserves $e^{-H}$ and moves the chain to a new energy level and a new trajectory each iteration; the composition of the two invariant kernels is invariant and, with the refresh, irreducible.

**C2.** (a) The state is $(x_1,\dots,x_K)$ with target $\prod_i\pi(x_i)^{1/T_i}$. Swapping $x_i$ and $x_{i+1}$ is a symmetric deterministic proposal (an involution with Jacobian 1), so by (3.5) the acceptance is $\min(1,r)$ with
$$r=\frac{\pi(x_{i+1})^{1/T_i}\,\pi(x_i)^{1/T_{i+1}}}{\pi(x_i)^{1/T_i}\,\pi(x_{i+1})^{1/T_{i+1}}},\qquad\log r=\big(\log\pi(x_{i+1})-\log\pi(x_i)\big)\Big(\frac1{T_i}-\frac1{T_{i+1}}\Big).$$
The within-chain moves are random-walk Metropolis on $\pi^{1/T_i}$, accepting with $\exp\{(\log\pi(y)-\log\pi(x_i))/T_i\}$.
(b) The reference implementation is `ptemper` in Part 9. Two replicas of 8000 iterations started at $-4$ and $+4$ (1000 discarded): $\hat R$ of the $T=1$ chains $1.002$, against $3.65$ for the single-temperature chains of Part 5; $P(x>0)=0.482$ (truth $0.5$); swap acceptance $0.81$. At $T=8$ the barrier between modes ($\log\pi$ drops by about 8 between them) is flattened to about 1 unit, so the hot chain crosses freely and swaps carry its mode changes down to $T=1$.
(c) Every move (within-chain Metropolis steps and swaps) is an MH step for the product target, so the product is stationary; its marginal for the first coordinate is $\pi^{1/1}=\pi$. The hot chains are auxiliary, like HMC's momentum.
