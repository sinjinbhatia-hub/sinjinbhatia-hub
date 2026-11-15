# ST-5 Solutions

Equation numbers refer to `notes.md`. Every number below is asserted or printed by Part 14 of `project.py` (or by the part named in the notes).

## A. Warm-up

**A1.** (a) By (3.4), Beta$(1+4,\ 3+2)=$ Beta(5,5).

(b) Prior mean $1/4$, sample fraction $4/6$, data weight $w=n/(\alpha_0+\beta_0+n)=6/10=0.6$. Then $0.6\cdot\frac23+0.4\cdot\frac14=0.4+0.1=0.5=\frac5{10}$, the posterior mean. ✓

(c) Given $p$, the two trials are independent, so $P(\text{both})=\mathbb E[p^2\mid y]$. For Beta$(\alpha,\beta)$, $\mathbb E[p^2]=\frac{\alpha(\alpha+1)}{(\alpha+\beta)(\alpha+\beta+1)}=\frac{5\cdot6}{10\cdot11}=0.2727$. The squared mean is $0.25$. They differ by $\mathrm{Var}(p\mid y)=\frac{25}{100\cdot11}=0.0227$, since $\mathbb E[p^2]=(\mathbb E p)^2+\mathrm{Var}(p)$. The two future trials share the unknown $p$, so they are positively correlated; the plug-in answer $0.25$ treats them as independent.

**A2.** (a) Precision: $1/16+6/9=0.0625+0.6667=0.72917$, so $\tau_n^2=1.3714$. Mean: $\mu_n=\tau_n^2\,(20/16+6\cdot23/9)=1.3714\,(1.25+15.3333)=22.7429$. The posterior is $\mathcal N(22.7429,\ 1.3714)$, sd $1.1711$, and the 95% interval is $22.7429\pm1.96(1.1711)=(20.4475,\ 25.0382)$.

(b) By (3.11), predictive variance $9+1.3714=10.3714$, sd $3.2205$; interval $22.7429\pm6.3122=(16.4307,\ 29.0550)$.

(c) $n_0=\sigma^2/\tau_0^2=9/16=0.5625$ imaginary observations. Check: $\mu_n=(0.5625\cdot20+6\cdot23)/6.5625=149.25/6.5625=22.7429$. ✓

**A3.** (a) $\sum y=6$, $n=3$: Gamma$(3+6,\ 2+3)=$ Gamma(9,5) by (3.16).

(b) For $\tilde t=1$, (3.17) gives mean $\alpha_n/\beta_n=9/5=1.8$ and variance $1.8\,(1+1/5)=2.16$. (The Poisson plug-in at $1.8$ would say variance $1.8$.)

(c) $P(\tilde y=0)=(\beta_n/(\beta_n+1))^{\alpha_n}=(5/6)^9=0.1938$.

**A4.** (a) $\log p=\log\lambda-\lambda x$, $\partial^2_\lambda\log p=-1/\lambda^2$, which does not involve $x$, so $I_1(\lambda)=1/\lambda^2$ and $p_J\propto\sqrt{I_1}=1/\lambda$.

(b) Likelihood $\prod\lambda e^{-\lambda x_i}=\lambda^ne^{-\lambda s}$; times $\lambda^{-1}$: $\lambda^{n-1}e^{-\lambda s}$, the kernel of Gamma$(n,s)$.

(c) $n=3$, $s=6.2$: Gamma(3, 6.2), mean $3/6.2=0.4839$, which equals the MLE $n/s$. This is not a coincidence: Jeffreys' prior is the Gamma$(0,0)$ limit, and the Gamma$(\alpha,\beta)$ posterior mean $(\alpha+n)/(\beta+s)$ reduces to $n/s$. The 95% interval is $(0.0998,\ 1.1653)$: wide and asymmetric, as three observations deserve.

**A5.** (a) Likelihood $(1-p)^{10}$ times the flat prior gives Beta$(1,11)$ with density $11(1-p)^{10}$. Integrating from $0$ to $q$: $1-(1-q)^{11}$.

(b) Lower endpoint: $1-(1-q)^{11}=0.025\Rightarrow q=1-0.975^{1/11}=0.0023$. Upper: $1-(1-q)^{11}=0.975\Rightarrow q=1-0.025^{1/11}=0.2849$. Interval $(0.0023,\ 0.2849)$.

(c) The density $11(1-p)^{10}$ is strictly decreasing, so $\{p:\ \text{density}\ge c\}=[0,q]$. Mass $0.95$: $1-(1-q)^{11}=0.95\Rightarrow q=1-0.05^{1/11}=0.2384$. HPD $=(0,\ 0.2384)$, shorter than the equal-tailed interval (length $0.2384$ against $0.2826$). Its left endpoint is $0$ because the mode of the posterior is at the boundary; excluding the values near $0$ (where density is highest) to balance tails would throw away the most probable region.

## B. Core

**B1.** *Pmf.* By (3.3), $P(\tilde y=k\mid y)=\int\binom mk p^k(1-p)^{m-k}\frac{p^{\alpha_n-1}(1-p)^{\beta_n-1}}{B(\alpha_n,\beta_n)}dp=\binom mk\frac{B(\alpha_n+k,\ \beta_n+m-k)}{B(\alpha_n,\beta_n)}$, using the Beta integral for the kernel $p^{\alpha_n+k-1}(1-p)^{\beta_n+m-k-1}$.

*Variance.* Eve's law: $\mathrm{Var}(\tilde y)=\mathbb E[m\,p(1-p)]+\mathrm{Var}(mp)$. Write $m_n=\mathbb E[p\mid y]$ and $V=\mathrm{Var}(p\mid y)=\frac{m_n(1-m_n)}{\alpha_n+\beta_n+1}$. Since $\mathbb E[p(1-p)]=m_n-(V+m_n^2)=m_n(1-m_n)-V$, we get $\mathrm{Var}(\tilde y)=m\,m_n(1-m_n)-mV+m^2V=m\,m_n(1-m_n)+m(m-1)V=m\,m_n(1-m_n)\big[1+\frac{m-1}{\alpha_n+\beta_n+1}\big]=m\,m_n(1-m_n)\frac{\alpha_n+\beta_n+m}{\alpha_n+\beta_n+1}$.

*Numbers.* $\alpha_n+\beta_n=14$, $m_n=9/14$, $m\,m_n(1-m_n)=10\cdot\frac9{14}\cdot\frac5{14}=2.2959$. Variance $=2.2959\cdot\frac{24}{15}=3.6735$; ratio to the Binomial$(10,9/14)$ variance $2.2959$ is $1.6$. For all ten successes, $\frac{B(19,5)}{B(9,5)}=\prod_{j=0}^{9}\frac{9+j}{14+j}=\frac{9\cdot10\cdots18}{14\cdot15\cdots23}=0.03825$, against $(9/14)^{10}=0.01205$ from the plug-in: the plug-in forecast is three times too pessimistic about a clean sweep, because it ignores that if $p$ is high then everything comes up together.

**B2.** (a) Reference-prior posterior: $\kappa_n=n$, $\mu_n=\bar y$, $\alpha_n=\frac{n-1}2$, $\beta_n=\frac S2$ with $S=(n-1)s^2$. By (3.15) the predictive is $t_{2\alpha_n}=t_{n-1}$ with scale$^2=\frac{\beta_n(\kappa_n+1)}{\alpha_n\kappa_n}=\frac{S/2\cdot(n+1)}{\frac{n-1}2\cdot n}=s^2\frac{n+1}n=s^2(1+\frac1n)$. ∎

(b) $s\sqrt{1.1}=0.04169\cdot1.04881=0.04373$; $t_{9,0.975}=2.2622$; half-width $0.09890$; interval $(9.7071,\ 9.9049)$. The plug-in $9.806\pm1.96(0.04169)=(9.7243,\ 9.8877)$ is narrower by 17%: it ignores the uncertainty in $\mu$ and uses $z$ where $t_9$ is needed. The first interval is also the frequentist 95% prediction interval (ST-3 (3.17) for the intercept-only model).

(c) The sentence is correct in the Bayesian framework with the reference prior: $P(g\in(9.7762,9.8358)\mid\text{data})=0.95$, with $g$ random (belief) and the interval fixed. In ST-2's framework $g$ is a fixed constant and the *interval* is the random object; "95%" describes the long-run frequency with which the procedure captures $g$, and once the data are in, the realized interval either contains $g$ or does not. What makes the numbers equal is that the quantity $(\mu-\bar y)/(s/\sqrt n)$ has a $t_{n-1}$ law whether the data are random with $\mu$ fixed or $\mu$ is random (reference prior) with the data fixed. With the proper NIG prior of Example 4.3 the interval is $(9.7133,\ 9.8430)$ and the equality is gone.

**B3.** $p(y\mid\lambda)=\lambda e^{-\lambda y}=\exp\{\eta y-A(\eta)\}$ with $T=y$, $\eta=-\lambda<0$, $A(\eta)=-\log(-\eta)=-\log\lambda$, $h=1$.

(a) Conjugate prior (3.18): $\exp\{\eta\chi-\nu A\}=\exp\{\eta\chi+\nu\log(-\eta)\}=(-\eta)^\nu e^{\eta\chi}$. Put $\lambda=-\eta$ (Jacobian $1$): $\lambda^\nu e^{-\lambda\chi}$, the kernel of Gamma$(\nu+1,\ \chi)$.

(b) (3.19): $\chi'=\chi+s$, $\nu'=\nu+n$, so the posterior is Gamma$(\nu+n+1,\ \chi+s)$. In Gamma terms, $\alpha_n=\alpha_0+n$, $\beta_n=\beta_0+s$ (with $\alpha_0=\nu+1$, $\beta_0=\chi$), the same as A4(b) in the limit $\alpha_0,\beta_0\to0$.

(c) Prior Gamma(2,4) means $\nu=1$, $\chi=4$. With $n=3$, $s=6.2$: Gamma$(5,\ 10.2)$. In (3.20) the quantity averaged is $\nabla A(\eta)=-1/\eta=1/\lambda$, the **mean lifetime**: $\mathbb E[1/\lambda\mid y]=\frac{\chi+s}{\nu+n}=\frac{4+6.2}{4}=2.55=\frac14(4)+\frac34(2.0667)$, weights $\nu=1$ and $n=3$ on the prior value $\chi/\nu=4$ and the sample mean lifetime $\bar y=2.0667$. (Check: for Gamma$(5,10.2)$, $\mathbb E[1/\lambda]=\beta/(\alpha-1)=10.2/4$. ✓) The posterior mean of $\lambda$ itself is $5/10.2=0.4902$; it is not a weighted average of the prior mean rate $0.5$ and the sample rate $0.4839$ with these weights, because $\mathbb E[1/\lambda]\ne1/\mathbb E[\lambda]$ (Jensen, F5 Theorem 3.32): conjugate shrinkage is linear in the *mean parameter* $\nabla A$, not in other transforms.

(d) $P(\tilde X>t\mid y)=\int e^{-\lambda t}\frac{\beta_n^{\alpha_n}\lambda^{\alpha_n-1}e^{-\beta_n\lambda}}{\Gamma(\alpha_n)}d\lambda=\frac{\beta_n^{\alpha_n}}{(\beta_n+t)^{\alpha_n}}$, since the integral is a Gamma kernel with rate $\beta_n+t$. At $t=1$: $(10.2/11.2)^5=0.6265$. (The MLE plug-in $e^{-3/6.2}=0.6165$ is lower: the posterior's tail is heavier.)

**B4.** (a) For $\mathcal N(0,v)$, ST-1's variance information is $I_v=\frac1{2v^2}$. With $\phi=1/v$, $dv/d\phi=-1/\phi^2$, so by (3.13) $I_\phi=I_v\,\phi^{-4}=\frac{\phi^2}2\cdot\frac1{\phi^4}=\frac1{2\phi^2}$ and $p_J(\phi)\propto1/\phi$. Pushing $p_J(v)\propto1/v$ to $\phi$: $p(\phi)=\frac1{1/\phi}\cdot\frac1{\phi^2}=\frac1\phi$. The two agree, as Theorem 3.8 requires.

(b) Likelihood $v^{-5}e^{-16/(2v)}$ times $v^{-1}$ is $v^{-6}e^{-8/v}$, the kernel of IG$(5,8)$ (power $-a-1=-6$). By (3.12) with $a_0=b_0=0$ this is IG$(\frac n2,\frac12\sum y^2)$. ✓

(c) Mean $b/(a-1)=8/4=2.0$; mode $b/(a+1)=8/6=1.333$; 95% interval $(0.7811,\ 4.9277)$. The MLE $\sum y^2/n=1.6$ lies between the mode and the mean; the posterior is right-skewed, so mean $>$ mode, and the Jeffreys MAP $\sum y^2/(n+2)$ is below the MLE.

**B5.** (a) $\sum y=3+9=12$ events, $\sum t=2.5+4=6.5$ hours. By (3.16), Gamma$(2+12,\ 1+6.5)=$ Gamma(14, 7.5), mean $1.867$ events/hour.

(b) For $\tilde t=1$ (3.17): mean $14/7.5=1.8667$; variance $1.8667(1+1/7.5)=2.1156$; $P(\tilde y=0)=(7.5/8.5)^{14}=0.1734$.

(c) MLE rate $12/6.5=1.846$. Poisson with that rate: variance $1.846$ (the index of dispersion is 1 by construction), $P(0)=e^{-1.846}=0.1578$. The plug-in understates the variance by 13% and the chance of a silent hour by 1.6 percentage points; both are the cost of ignoring uncertainty in $\lambda$.

**B6.** (a) By (3.26), $r(\delta)=\mathbb E_y\,\rho(\delta(y)\mid y)$. For $\delta=$ posterior mean, Theorem 3.12(a) gives $\rho(\delta(y)\mid y)=\mathrm{Var}(\theta\mid y)$. Hence $r=\mathbb E_y[\mathrm{Var}(\theta\mid y)]$. ∎

(b) By (3.5) with a uniform prior, $y$ is uniform on $\{0,\dots,10\}$ (each value has probability $1/11$). Posterior Beta$(y+1,\ 11-y)$, $\alpha_n+\beta_n+1=13$, so $\mathrm{Var}=m(1-m)/13$ with $m=(y+1)/12$. Averaging over $y$: $\frac1{13}\cdot\frac1{11}\sum_{k=1}^{11}\frac k{12}\big(1-\frac k{12}\big)=\frac1{13\cdot11\cdot144}\sum_{k=1}^{11}k(12-k)$. Since $\sum k=66$ and $\sum k^2=506$, $\sum k(12-k)=12(66)-506=286$. So $r=\frac{286}{13\cdot11\cdot144}=\frac{286}{20592}=\frac1{72}=0.01389$. For the MLE, $R(p,\hat p)=p(1-p)/n$ does not involve $y$ beyond the Binomial variance, so $r=\mathbb E[p(1-p)]/n$ with $\mathbb E[p(1-p)]=\frac12-\frac13=\frac16$ under the uniform prior: $r=\frac1{60}=0.01667$. Laplace's estimator is smaller (1/72 against 1/60), as it must be: it is the posterior mean, hence the Bayes rule for squared loss, and no rule has lower Bayes risk.

**B7.** (a) Beta(121,881), Beta(146,856), Beta(139,863).

(b) `prob_B_beats_A(120,1000,145,1000)` $=\int f_A(t)\,[1-F_B(t)]\,dt=0.9502$. Monte Carlo with 400,000 paired draws gives $0.9501$ (standard error $\approx0.0003$).

(c) Draw 400,000 samples from each posterior and take the arg-max per draw: $P(\text{best})\approx(0.022,\ 0.659,\ 0.319)$ for A, B, C (each $\pm0.001$). B is most likely best, but C has a one-in-three chance; the data do not separate B and C well.

(d) With Beta(20,180) priors (mean $0.10$, worth 200 trials) the posteriors are Beta(140,1060) and Beta(165,1035), with means $0.1167$ and $0.1375$ instead of $0.1208$ and $0.1457$: both are pulled toward $0.10$, the gap shrinks from $0.025$ to $0.021$, and $P(B>A)$ falls from $0.9502$ to $0.9376$. Neither is a magic number: a decision rule should compare expected loss, not threshold the probability at $0.95$.

**B8.** (a) Under $H_0$, $\bar x\sim\mathcal N(0,v_0)$ with $v_0=\sigma^2/n$. Under $H_1$, $\bar x=\mu+e$ with $\mu\sim\mathcal N(0,\tau^2)$ and $e\sim\mathcal N(0,v_0)$ independent, so $\bar x\sim\mathcal N(0,v_0+\tau^2)$ (F5 Theorem 3.22). The ratio of the densities at $\bar x$ is $\sqrt{\frac{v_0}{v_0+\tau^2}}\exp\{-\frac{\bar x^2}2[\frac1{v_0+\tau^2}-\frac1{v_0}]\}$. With $u=\tau^2/v_0=n\tau^2/\sigma^2$: $\sqrt{v_0/(v_0+\tau^2)}=(1+u)^{-1/2}$, and $\frac1{v_0}-\frac1{v_0+\tau^2}=\frac1{v_0}\cdot\frac u{1+u}$; with $z^2=\bar x^2/v_0$ the exponent is $\frac{z^2}2\frac u{1+u}$. That is (3.40).

(b) $\frac d{du}\log\mathrm{BF}_{10}=-\frac1{2(1+u)}+\frac{z^2}{2(1+u)^2}=\frac{z^2-(1+u)}{2(1+u)^2}$, positive for $1+u<z^2$ and negative after, so the maximum is at $1+u=z^2$ (requires $|z|>1$), where $\mathrm{BF}_{10}=(z^2)^{-1/2}\exp\{\frac{z^2}2\cdot\frac{z^2-1}{z^2}\}=\frac1{|z|}e^{(z^2-1)/2}$. For any prior $\pi$ on $\mu$, $p(\bar x\mid H_1)=\int\mathcal N(\bar x;\mu,v_0)\pi(\mu)d\mu\le\mathcal N(\bar x;\bar x,v_0)$, so $\mathrm{BF}_{10}\le\mathcal N(\bar x;\bar x,v_0)/\mathcal N(\bar x;0,v_0)=e^{\bar x^2/(2v_0)}=e^{z^2/2}$.

(c) $z=1.96$: $\frac1{1.96}e^{1.4208}=2.112$ and $e^{1.9208}=6.826$. $z=2.576$: $\frac1{2.576}e^{2.8179}=6.499$ and $e^{3.3179}=27.60$. At $p=0.05$, the best a Normal prior can do for the alternative is about 2.1 to 1, and the best any prior can do is about 6.8 to 1: "$p=0.05$ is strong evidence" is not supported. Even $p=0.01$ gives at most about 6.5 to 1 (Normal) or 27.6 to 1 (any prior).

## C. Challenge

**C1.** (a) $\delta=(y+a)/(n+2a)$. Variance: $\mathrm{Var}(y)/(n+2a)^2=np(1-p)/(n+2a)^2$. Bias: $\frac{np+a}{n+2a}-p=\frac{a(1-2p)}{n+2a}$. So $R=\frac{np(1-p)+a^2(1-2p)^2}{(n+2a)^2}$. With $a^2=n/4$ the numerator is $np-np^2+\frac n4(1-4p+4p^2)=np-np^2+\frac n4-np+np^2=\frac n4$, free of $p$. Hence $R=\frac{n/4}{(n+\sqrt n)^2}$; for $n=100$, $\frac{25}{12100}=0.002066$.

(b) Let $\delta$ be the Bayes rule for a prior $\pi$ with $R(p,\delta)=r$ for all $p$. Its Bayes risk is $\int R\,\pi=r$. Suppose another rule $\delta'$ had $\sup_pR(p,\delta')<r$. Its Bayes risk $\int R(p,\delta')\pi(p)dp\le\sup_pR(p,\delta')<r$, contradicting that $\delta$ minimizes Bayes risk among all rules. So $\sup_pR(p,\delta')\ge r$ for every rule: $\delta$ is minimax. ∎

(c) The MLE's maximum risk is $1/(4n)=0.0025$ (at $p=\frac12$), larger than $0.002066$. But $p(1-p)/n<0.002066$ whenever $p(1-p)<0.2066$, i.e. $p<0.2917$ or $p>0.7083$: outside that middle band the MLE is better. The minimax rule is preferable only if $p$ could be anywhere and the worst case is what you fear; if $p$ is known to be small (rare events), the MLE, or a Bayes rule with a matching prior, is far better.

**C2.** (a) $\mathrm{BF}_{01}=(1+u)^{1/2}\exp\{-\frac{z^2}2\frac u{1+u}\}$ (the reciprocal of (3.40)). Setting it to 1 with $z=2$: $\frac12\log(1+u)=\frac{2u}{1+u}$. Numerically $u=49.44$. With $\sigma=\tau=1$, $u=n$, so $H_0$ is favored for $n\ge50$ (at $\bar x=z/\sqrt n=0.28$).

(b) For large $u$, $(1+u)^{1/2}\approx\sqrt u$ and $\frac u{1+u}\to1$, so $\mathrm{BF}_{01}\approx\sqrt u\,e^{-z^2/2}$. At $u=10^4$, $z=2$: $100\,e^{-2}=13.53$, against the exact $13.54$ of Example 4.10.

(c) $u=n\tau^2/\sigma^2=n/100$ for $\tau=0.1$, so $n=100(49.44)\approx4944$: a narrower alternative takes more data before it loses. Extension for Part 11: `taus = np.logspace(-3, 1, 200); plt.loglog(taus, 1/bf10(0.02, 10**4, 1, taus)); plt.axvline(np.sqrt(3)/100)`. The curve starts at $1$ for tiny $\tau$ ($H_1\to H_0$), dips to $1/2.241=0.446$ at $\tau=\sqrt{z^2-1}\,\sigma/\sqrt n=0.0173$, crosses $1$ at $\tau=\sqrt{49.44}/100=0.0703$ (the same $u$ as in (a)), and then rises without bound. The data favor $H_1$ for $\tau<0.0703$ and $H_0$ above: the verdict is a function of the prior scale.

**C3.** (a) With the prior variance set to $10^{12}$, Newton returns the MLE $\hat\beta=(0.041,\ 1.711)$, with $\hat\ell=-31.533$. For $n=60$, $d=2$: $-\mathrm{BIC}/2=\hat\ell-\frac d2\log n=-31.533-4.094=-35.628$.

(b) The grid log evidence is $-35.822$ and the Laplace value is $-35.847$ (they differ by $0.025$). BIC is off by $0.19$ from the grid: a factor $e^{0.19}\approx1.2$ in the evidence.

(c) Expanding about the MLE with Hessian $-n\hat J$: $\log p(y)\approx\hat\ell+\log p(\hat\beta)+\frac d2\log2\pi-\frac d2\log n-\frac12\log|\hat J|$ (the derivation before (3.42)). BIC keeps the first and fourth terms and drops three: the log prior density at $\hat\beta$, $\log\mathcal N(\hat\beta;0,9I)=-4.198$; $\frac d2\log2\pi=1.838$; and $-\frac12\log|\hat J|=2.176$, where $\hat J=X^\top\hat WX/n$. Their sum is $-0.184$. Then $-35.628-0.184=-35.812$, within $0.010$ of the grid value and $0.035$ of the Laplace value. So the gap is exactly the dropped constants; none grows with $n$, but each depends on the prior or on the information per observation, and in a comparison between models the difference of these terms can exceed the BIC difference itself.
