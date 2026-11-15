# ST-9 Solutions

Every number below is printed or asserted by `project.py` (Part 13 unless another part is cited).

## A. Warm-up

**A1.** Individual effects $Y^1-Y^0=(2,1,3,0,2,1)$.
- ATE $=(2+1+3+0+2+1)/6=9/6=1.5$.
- ATT (units 1–3) $=(2+1+3)/3=2.0$.
- $\mathbb E[Y\mid A=1]=$ mean of $Y^1$ over units 1–3 $=(6+7+8)/3=7$. $\mathbb E[Y\mid A=0]=$ mean of $Y^0$ over units 4–6 $=(8+7+9)/3=8$. Association $=-1.0$.
- Selection bias $=\mathbb E[Y^0\mid A=1]-\mathbb E[Y^0\mid A=0]=(4+6+5)/3-8=5-8=-3.0$.
- Check: ATT $+$ selection bias $=2.0-3.0=-1.0$. ✓

The treated had worse untreated prospects by 3 points, which more than hides an effect of 2 on them. The ATT (2.0) also differs from the ATE (1.5) because the treated happen to have larger effects.

**A2.** Paths from $A$ to $Y$ (ignoring direction):
1. $A\to Y$ (causal).
2. $A\to M\to Y$ (causal).
3. $A\leftarrow L\to Y$ (backdoor: starts with an arrow into $A$).
4. $A\to C\leftarrow Y$ (not backdoor; $C$ is a collider, blocked unless conditioned on).

The only backdoor path is 3, and it is blocked by $L$. $M$ and $C$ are descendants of $A$, so any set containing them violates condition (i). Checking the eight subsets:
- $\varnothing$: path 3 open. Fails.
- $\{L\}$: blocks path 3, no descendants of $A$. **Valid.**
- $\{M\}$, $\{L,M\}$: contain a descendant of $A$ (a mediator); adjusting removes the indirect effect along path 2.
- $\{C\}$, $\{L,C\}$: contain a descendant of $A$ (a collider); adjusting opens path 4.
- $\{M,C\}$, $\{L,M,C\}$: both problems.

The project's `backdoor_ok` returns exactly $[\{L\}]$.

**A3.** Pooled outcomes $\{5,7,2,3\}$, total 17; observed difference $6-2.5=3.5$. There are $\binom42=6$ assignments. A treated pair with sum $s$ gives difference $(s-(17-s))/2=s-8.5$. The pairs and sums: $\{7,5\}=12$, $\{7,3\}=10$, $\{7,2\}=9$, $\{5,3\}=8$, $\{5,2\}=7$, $\{3,2\}=5$. Only the observed pair reaches $s\ge12$, so the one-sided $p=1/6=0.1667$.

That is also the smallest p-value the design can produce, whatever the data. With four units the experiment cannot reach significance at 0.05 even if the treatment is enormous. Compute the minimum attainable p-value, $1/\binom N{n_1}$ one-sided, before you run a small randomized experiment; with $N=6$ it is $0.05$, with $N=8$ (4 treated) it is $1/70=0.014$.

**A4.** ITT $=12.6-12.0=0.6$. First stage $=0.45-0.15=0.30$. Wald $=0.6/0.30=2.00$.

Under monotonicity, $P(A=1\mid Z=0)=0.15$ is the always-taker share, $1-0.45=0.55$ is the never-taker share, and the compliers are $1-0.15-0.55=0.30$, the first stage. The estimate 2.0 is the average effect of $A$ among those who take the treatment if and only if encouraged. It says nothing about the always-takers (15%) or never-takers (55%).

**A5.** Standard order $(--,+-,-+,++)=(20,30,24,26)$.
- $A$: $\frac12(-20+30-24+26)=6$.
- $B$: $\frac12(-20-30+24+26)=0$.
- $AB$: $\frac12(20-30-24+26)=-4$.

The effect of $A$ is $30-20=10$ at low $B$ and $26-24=2$ at high $B$; the interaction is half the difference, $(2-10)/2=-4$. $B$'s main effect is zero, but its effect is $24-20=+4$ at low $A$ and $26-30=-4$ at high $A$. "B has no effect" is false: $B$ has effects of opposite sign that cancel on average. With a large interaction, report the simple effects (or the interaction plot), not the main effects.

**A6.** DiD $=(15-10)-(11-8)=5-3=2$.

Assumption (parallel trends): without the treatment, the treated group's average would have risen by the same 3 points as the control group's. Failure example: the treated group is a set of schools that adopted a new curriculum *because* their scores were already improving (or because a new principal arrived the same year), so their untreated trend was steeper. Then part of the 5-point rise was trend, and DiD overstates the effect. The project's Part 9 shows a trend difference of 0.7 producing an estimate of 2.699 for a true effect of 2.

## B. Core

**B1.** Marginals: $n_{L=0}=200$, $n_{L=1}=300$, so $P(L=0)=0.4$, $P(L=1)=0.6$. Propensities $P(A=1\mid L=0)=40/200=0.2$, $P(A=1\mid L=1)=240/300=0.8$. Stratum risks: $(L,A)=(0,0):0.30$, $(0,1):0.50$, $(1,0):0.50$, $(1,1):0.70$.

(a) $\mathbb E[Y\mid A=1]=(20+168)/280=0.6714$, $\mathbb E[Y\mid A=0]=(48+30)/220=0.3545$; crude difference $0.3169$.

(b) $\mathbb E[Y^1]=0.50(0.4)+0.70(0.6)=0.20+0.42=0.62$. $\mathbb E[Y^0]=0.30(0.4)+0.50(0.6)=0.12+0.30=0.42$. ATE $=0.20$.

(c) Weights $1/P(A\mid L)$: $(0,1)$: 5; $(0,0)$: 1.25; $(1,1)$: 1.25; $(1,0)$: 5.
$\mathbb E[Y^1]=(20\cdot5+168\cdot1.25)/500=(100+210)/500=0.62$. $\mathbb E[Y^0]=(48\cdot1.25+30\cdot5)/500=(60+150)/500=0.42$. Same as (b), as the saturated-case argument in §3.4 guarantees.

(d) Weighted counts:

| $L$ | $A$ | weighted $n$ | weighted $Y=1$ |
|---|---|---|---|
| 0 | 0 | $160\times1.25=200$ | 60 |
| 0 | 1 | $40\times5=200$ | 100 |
| 1 | 0 | $60\times5=300$ | 150 |
| 1 | 1 | $240\times1.25=300$ | 210 |

In each arm, $P(L=1\mid A)=300/500=0.6=P(L=1)$: in the pseudo-population treatment is independent of $L$. Each arm is a full copy of the 500-patient population, and the crude difference in the pseudo-population, $310/500-210/500=0.20$, is the ATE.

**B2.** (a) $Y^1=(6,7,8,8,9,10)$: mean 8, squared deviations $4,1,0,0,1,4$, sum 10, $S_1^2=10/5=2.0$. $Y^0=(4,6,5,8,7,9)$: mean 6.5, squared deviations $6.25,0.25,2.25,2.25,0.25,6.25$, sum 17.5, $S_0^2=3.5$. $\tau=(2,1,3,0,2,1)$: mean 1.5, squared deviations $0.25,0.25,2.25,2.25,0.25,0.25$, sum 5.5, $S_\tau^2=1.1$.
$$\operatorname{Var}\hat\tau=\frac{2.0}3+\frac{3.5}3-\frac{1.1}6=0.6667+1.1667-0.1833=1.6500.$$

(b) The project's `neyman_exact(y0, y1, 3)` enumerates the 20 assignments; the variance of the 20 values of $\hat\tau$ is 1.6500, and their mean is 1.5.

(c) The sample variances $s_1^2,s_0^2$ within the treated and control groups are unbiased over the randomization for $S_1^2,S_0^2$ (a simple random sample without replacement from a finite population has $\mathbb Es^2=S^2$). So the conservative estimator has expectation $S_1^2/3+S_0^2/3=1.8333$. It overstates the variance by $S_\tau^2/N=1.1/6=0.1833$, about 11%.

(d) $\mathbb E[s_1^2/n_1+s_0^2/n_0]-\operatorname{Var}\hat\tau=S_\tau^2/N$ by (c) and (3.6). This is zero if and only if $S_\tau^2=0$, i.e. $\frac1{N-1}\sum_i(\tau_i-\bar\tau)^2=0$, i.e. $\tau_i=\bar\tau$ for every $i$: a constant effect. ∎

**B3.** (a) Treks with unit $\omega$'s. $\operatorname{Var}A$: top $L$ (path $L\to A$ twice, $1\cdot1$) and top $A$ (trivial, 1): $\operatorname{Var}A=2$. $\operatorname{Cov}(A,Y)$: top $L$: left path $L\to A$ (1); right paths to $Y$: $L\to Y$ (1) and $L\to A\to Y$ (0.5), total 1.5, product $1.5$; top $A$: trivial and $A\to Y$, $0.5$. Total $\operatorname{Cov}(A,Y)=2.0$.

(b) $Y$ on $A$ alone: $2.0/2=1.000$. $Y$ on $(A,L)$: $L$ satisfies the backdoor criterion, so the coefficient is the effect, $0.500$. Directly: $\operatorname{Cov}(A,L)=1$, $\operatorname{Var}L=1$, $\operatorname{Cov}(L,Y)=1+0.5=1.5$; solving $\begin{pmatrix}2&1\\1&1\end{pmatrix}\beta=\begin{pmatrix}2\\1.5\end{pmatrix}$ gives $\beta=(0.5,1)$. OVB: the short coefficient is $0.5+1\cdot\gamma$ with $\gamma=\operatorname{Cov}(A,L)/\operatorname{Var}A=1/2$, so $1.0$. ✓ The bias term $b\gamma=0.5$ is the confounding through the fork.

(c) Add $C=A+Y+\varepsilon_C$. Covariances: $\operatorname{Var}Y=0.25\cdot2+1\cdot1+2(0.5)(1)\operatorname{Cov}(A,L)+1=0.5+1+1+1=3.5$. $\operatorname{Cov}(A,C)=\operatorname{Var}A+\operatorname{Cov}(A,Y)=4$. $\operatorname{Cov}(L,C)=\operatorname{Cov}(L,A)+\operatorname{Cov}(L,Y)=1+1.5=2.5$. $\operatorname{Cov}(Y,C)=\operatorname{Cov}(A,Y)+\operatorname{Var}Y=5.5$. $\operatorname{Var}C=\operatorname{Var}A+\operatorname{Var}Y+2\operatorname{Cov}(A,Y)+1=2+3.5+4+1=10.5$. Regress $Y$ on $(A,L,C)$: by FWL first partial out $L$. Given $L$: $\operatorname{Var}(A\mid L)=1$, $\operatorname{Cov}(A,Y\mid L)=2-1.5=0.5$, $\operatorname{Var}(Y\mid L)=3.5-2.25=1.25$, $\operatorname{Cov}(A,C\mid L)=4-2.5=1.5$, $\operatorname{Cov}(Y,C\mid L)=5.5-3.75=1.75$, $\operatorname{Var}(C\mid L)=10.5-6.25=4.25$. These are exactly Example 4.5's numbers: within levels of $L$ the problem *is* Example 4.5, so the coefficient on $A$ is $-0.2500$. `SCM.pop_coef("Y", ["A","L","C"])` returns $-0.2500$. Adjusting for the confounder fixed the fork; adding the collider broke it again and flipped the sign.

**B4.** (a) $A\to S\leftarrow Y$, with no arrow between $A$ and $Y$. Analysing only $S=1$ conditions on the collider.

(b) Bayes within $A=a$ (since $A\perp Y$, $P(Y=1\mid A=a)=0.2$):
$$P(Y=1\mid A=a,S=1)=\frac{0.2\,P(S=1\mid a,1)}{0.2\,P(S=1\mid a,1)+0.8\,P(S=1\mid a,0)}.$$
- $a=1$: $P(S\mid1,1)=0.9$, $P(S\mid1,0)=0.5$: $\frac{0.18}{0.18+0.40}=0.3103$.
- $a=0$: $P(S\mid0,1)=0.5$, $P(S\mid0,0)=0.1$: $\frac{0.10}{0.10+0.08}=0.5556$.

Ratio $0.5586$: in the study sample, exposure appears to *halve* the risk of the outcome, although they are independent in the population.

(c) "Being exposed was itself a reason to enter the study, so the exposed patients in the sample did not need the disease to get in, and the unexposed ones mostly did: the sample is a subset chosen by the very variables you are comparing."

**B5.** (a) Direct $0.5$ plus indirect $1\times0.8$: total effect $1.3$.

(b) Backdoor paths start $A\leftarrow L_1$: the only one is $A\leftarrow L_1\to L_2\to Y$. It is blocked by $L_1$, by $L_2$, or by both. $\varnothing$ leaves it open. Valid: $\{L_1\}$, $\{L_2\}$, $\{L_1,L_2\}$, three distinct sets. `backdoor_ok` agrees.

(c) Output of Part 13:

| set | backdoor | coefficient | $n\operatorname{Var}(\hat\beta_A)$ |
|---|---|---|---|
| $\varnothing$ | no | 1.8000 | 1.5700 |
| $\{L_1\}$ | yes | 1.3000 | 2.6400 |
| $\{L_2\}$ | yes | 1.3000 | 1.0933 |
| $\{L_1,L_2\}$ | yes | 1.3000 | 1.6400 |
| $\{M\}$ | no | 1.0000 | 3.7500 |
| $\{L_1,C\}$ | no | $-0.3681$ | 1.7793 |

Two checks by hand. For $\{L_2\}$: substituting $M$, $Y=1.3A+L_2+0.8\varepsilon_M+\varepsilon_Y$, and $\varepsilon_M,\varepsilon_Y$ are independent of $(A,L_2)$, so the coefficient is 1.3 with residual variance $0.64+1=1.64$; $\operatorname{Var}A=\operatorname{Var}L_2=2$, $\operatorname{Cov}(A,L_2)=1$, so $\operatorname{Var}(A\mid L_2)=2-\frac12=1.5$, and $1.64/1.5=1.0933$. For $\{L_1\}$: $Y=1.3A+L_1+\varepsilon_{L_2}+0.8\varepsilon_M+\varepsilon_Y$, residual variance $1+1.64=2.64$, $\operatorname{Var}(A\mid L_1)=1$, so 2.64. The $\{M\}$ row is 1.0 by coincidence of two errors: adjusting for the mediator removes the indirect $0.8$ (leaving 0.5), and the unblocked backdoor path adds $\operatorname{Cov}(A,L_1\mid M)/\operatorname{Var}(A\mid M)=(1/3)/(2/3)=0.5$.

(d) $\{L_2\}$ is the most precise valid set. Adjusting for a variable that predicts $Y$ (and not $A$) shrinks the residual variance in the numerator; adjusting for a variable that predicts only $A$ shrinks $\operatorname{Var}(A\mid S)$ in the denominator, which *increases* the variance. Rule: among valid adjustment sets, prefer those made of causes of the outcome over causes of the treatment. (Adding $L_1$ to $\{L_2\}$ costs precision and buys nothing.) This is the start of the "optimal adjustment set" literature (Henckel, Perković & Maathuis 2022).

**B6.** (a) $\operatorname{Var}A=0.25+1+1=2.25$. $\operatorname{Cov}(A,Y)=\operatorname{Var}A+2\operatorname{Cov}(A,U)=2.25+2=4.25$. OLS plim $=4.25/2.25=1+2/2.25=1.8889$. The project's simulated OLS slope is 1.890 (Part 8).

(b) $\operatorname{Cov}(Z,A)=0.5$; $\operatorname{Cov}(Z,Y)=\operatorname{Cov}(Z,A)+2\operatorname{Cov}(Z,U)=0.5+0=0.5$. Ratio $1$. ✓

(c) Stage 1: $\hat A_i=\hat g_0+\hat g_1Z_i$ with $\hat g_1=\widehat{\operatorname{Cov}}(Z,A)/\widehat{\operatorname{Var}}Z$. Stage 2 slope: $\widehat{\operatorname{Cov}}(\hat A,Y)/\widehat{\operatorname{Var}}(\hat A)$. Since $\hat A$ is an affine function of $Z$, $\widehat{\operatorname{Cov}}(\hat A,Y)=\hat g_1\widehat{\operatorname{Cov}}(Z,Y)$ and $\widehat{\operatorname{Var}}\hat A=\hat g_1^2\widehat{\operatorname{Var}}Z$. The slope is
$$\frac{\hat g_1\widehat{\operatorname{Cov}}(Z,Y)}{\hat g_1^2\widehat{\operatorname{Var}}Z}=\frac{\widehat{\operatorname{Cov}}(Z,Y)}{\hat g_1\widehat{\operatorname{Var}}Z}=\frac{\widehat{\operatorname{Cov}}(Z,Y)}{\widehat{\operatorname{Cov}}(Z,A)}.\ ∎$$
For binary $Z$, $\widehat{\operatorname{Cov}}(Z,Y)/\widehat{\operatorname{Cov}}(Z,A)$ equals the ratio of differences in means (both covariances equal $\bar Z(1-\bar Z)$ times the difference in group means), which is (3.20).

(d) The first-stage residual variance is $\operatorname{Var}(A\mid Z)=1+1=2$ and the concentration parameter is $n\pi^2\operatorname{Var}Z/\operatorname{Var}(A\mid Z)=400\times0.25/2=50$. The expected $F$ with one instrument is about $1+50=51$; the project's 1000 simulations give a mean of 51.2. Far above 10: not weak.

(e) Stage-2 OLS computes residuals $y-\hat\beta\hat A$. The structural error is $y-\hat\beta A$, and $A-\hat A$ (the first-stage residual, which carries all of $U$'s influence) is in the stage-2 residual multiplied by $\hat\beta$. So the stage-2 residual variance estimates $\operatorname{Var}(\varepsilon+\beta v)$ instead of $\operatorname{Var}\varepsilon$, where $v$ is the first-stage error: wrong in general, in either direction. The correct estimator uses $y-X\hat\beta$ with the original $X$ in $\hat\sigma^2$, and $(X^\top P_ZX)^{-1}$ for the matrix (§3.8).

**B7.** (a) $\gamma\delta=\hat\beta$ gives $\delta=0.30/0.5=0.60$: the confounder would have to differ by 0.6 units (per unit of treatment) between the treated and untreated, after adjusting for $L$. The explain-away curve is the hyperbola $\delta=0.30/\gamma$; points above it reverse the sign. The project's `sensitivity.png` plots $\hat\beta-\gamma\delta$ against $\gamma$ for $\delta\in\{0.1,0.2,0.4,0.6\}$; the $\delta=0.6$ line crosses zero at $\gamma=0.5$. To also remove statistical significance you would need only $\gamma\delta=0.30-1.96(0.05)=0.202$.

(b) $E=1.5+\sqrt{1.5\times0.5}=1.5+0.866=2.366$. An unmeasured confounder associated with both treatment and outcome by a risk ratio of 2.37 each (above and beyond the measured covariates) could explain the observed risk ratio away; weaker confounding could not. Whether that is plausible is a subject-matter judgement: compare with the strongest measured confounder.

(c) The number of covariates says nothing about the ones you did not measure. (3.24) shows the bias depends only on the omitted variable's $\gamma$ and $\delta$. Forty weak covariates can leave one strong confounder untouched, and if any of the forty are colliders or mediators (§3.6), adjusting for them adds bias.

## C. Challenge

**C1.** (a) Under independence of $Z$ from $(A(0),A(1),Y^0,Y^1)$ and exclusion, $Y=Y^{A(Z)}$.
1. $\mathbb E[Y\mid Z=1]-\mathbb E[Y\mid Z=0]=\mathbb E[Y^{A(1)}\mid Z=1]-\mathbb E[Y^{A(0)}\mid Z=0]=\mathbb E[Y^{A(1)}-Y^{A(0)}]$.
2. $Y^{A(1)}-Y^{A(0)}=(A(1)-A(0))\tau$ with $\tau=Y^1-Y^0$: this is $\tau$ for compliers ($A(1)-A(0)=1$), $-\tau$ for defiers ($-1$), 0 for always- and never-takers.
3. Numerator $=\pi_c\mathbb E[\tau\mid c]-\pi_d\mathbb E[\tau\mid d]$.
4. First stage $=\mathbb E[A(1)-A(0)]=\pi_c-\pi_d$.

Dividing gives the stated ratio. ∎

(b) $\dfrac{0.4(2)-0.1(5)}{0.4-0.1}=\dfrac{0.8-0.5}{0.3}=1.00$. The compliers' effect is 2 and the defiers' is 5, and the Wald ratio is 1, below both. It is a weighted difference with a negative weight on defiers, not an average of anyone's effect; with larger defier effects it can even have the wrong sign.

(c) An encouragement letter for a programme that some recipients resent ("being told to do it makes me refuse") produces defiers. In judge-leniency designs, a judge may be harsher than average on one kind of case and more lenient on another, which violates monotonicity across defendant types.

**C2.** Markov factorization for $U\to A$, $U\to Y$, $A\to M$, $M\to Y$:
$$p(u,a,m,y)=p(u)\,p(a\mid u)\,p(m\mid a)\,p(y\mid m,u).$$
1. Truncated factorization (3.13): delete $p(a\mid u)$ and set $A=a$: $p(u,m,y\mid do(a))=p(u)\,p(m\mid a)\,p(y\mid m,u)$.
2. Sum over $u,m$: $p(y\mid do(a))=\sum_mp(m\mid a)\,\big[\sum_up(u)\,p(y\mid m,u)\big]$. The first factor is observable. The bracket involves the unobserved $U$.
3. Write $p(u)=\sum_{a'}p(u\mid a')p(a')$: the bracket is $\sum_{a'}p(a')\sum_up(u\mid a')p(y\mid m,u)$.
4. $M\perp U\mid A$ (the only path from $M$ to $U$ is $M\leftarrow A\leftarrow U$, a chain through $A$, plus $M\to Y\leftarrow U$, a collider not conditioned on), so $p(u\mid a')=p(u\mid a',m)$.
5. $Y\perp A\mid\{M,U\}$ (every path from $A$ to $Y$ goes through $M$ as a chain or through $U$ as a fork), so $p(y\mid m,u)=p(y\mid m,u,a')$.
6. Then $\sum_up(u\mid a',m)\,p(y\mid m,u,a')=p(y\mid m,a')$ by the law of total probability.
7. Combine: $p(y\mid do(a))=\sum_mp(m\mid a)\sum_{a'}p(y\mid m,a')\,p(a')$. ∎

Reading: the effect of $A$ on $M$ is unconfounded (backdoor $A\leftarrow U\to Y\leftarrow M$ is blocked at the collider $Y$), and the effect of $M$ on $Y$ is deconfounded by $A$ (backdoor $M\leftarrow A\leftarrow U\to Y$ is blocked by $A$). The formula chains the two.

**C3.** The project's Part 13 has the extension (1000 replications, $n=200$, $A=\pi\sum_{j=1}^kZ_j+U+\varepsilon_A$, $Y=A+2U+\varepsilon_Y$):

| $k$ | $\pi$ | mean $F$ | median bias | mean bias | OLS bias |
|---|---|---|---|---|---|
| 1 | 0.15 | 3.22 | 0.173 | 0.391 | 0.989 |
| 5 | 0.15 | 3.26 | 0.277 | 0.240 | 0.947 |

The OLS bias is $2/(k\pi^2+2)$, from $\operatorname{Var}A=k\pi^2+2$ and $\operatorname{Cov}(A,2U)=2$.

With $k=5$ the mean bias of 2SLS is $0.240$, a fraction $0.240/0.947=0.25$ of the OLS bias; the rule of thumb $1/F=1/3.26=0.31$ is in the right range. Both instrument sets have the same $F$, but the five-instrument estimator is pulled further toward OLS in median (0.277 vs 0.173): each added weak instrument adds a dimension to the projection $P_Z$ that captures noise correlated with $U$ (overfitting the first stage).

For $k=1$, the just-identified IV estimator is a ratio $\widehat{\operatorname{Cov}}(Z,Y)/\widehat{\operatorname{Cov}}(Z,A)$ whose denominator has positive density at zero. Its distribution has tails like a Cauchy's and no finite mean, so the sample mean of 1000 draws (0.391) is dominated by a few huge values and would change a lot with another seed. The median is the meaningful summary. With $k\ge2$ instruments and one regressor, 2SLS has a finite mean (moments exist up to order $k-1$; Kinal 1980), which is why the $k=5$ mean is stable.
