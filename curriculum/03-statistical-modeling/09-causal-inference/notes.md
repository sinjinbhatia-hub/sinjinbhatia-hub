# ST-9: Causal Inference and Experimental Design: Randomization, Confounding, DAGs, Adjustment, Instruments, the Limits of Regression

**Time:** ~8 hours · **Prerequisites:** ST-4 (logistic regression and IRLS, used for propensity scores); F4 (conditional probability, Simpson's paradox); ST-2 (permutation tests, power) and ST-3 (FWL, omitted-variable bias) are used throughout; ST-5 helpful · **Unlocks:** the honest use of every other module in this course and in 04-machine-learning

> **Conventions in this module.** $A$ is the treatment (binary, $A\in\{0,1\}$, unless stated), $Y$ the outcome, $L$ measured covariates, $U$ unmeasured variables, $Z$ an instrument. $Y^a$ (also $Y^{a=1}, Y^{a=0}$, or $Y^1,Y^0$) is the **potential outcome**: the value $Y$ would take if $A$ were set to $a$. The **individual causal effect** is $Y^1-Y^0$; the **average treatment effect** is $\mathrm{ATE}=\mathbb E[Y^1]-\mathbb E[Y^0]$; the **effect in the treated** is $\mathrm{ATT}=\mathbb E[Y^1-Y^0\mid A=1]$. The **association** is $\mathbb E[Y\mid A=1]-\mathbb E[Y\mid A=0]$. $P(Y\mid do(A=a))$ is the distribution of $Y$ when $A$ is *set* to $a$ by intervention, so $\mathbb E[Y\mid do(A=a)]=\mathbb E[Y^a]$. **Exchangeability** is $Y^a\perp A$; **conditional exchangeability** is $Y^a\perp A\mid L$. In a DAG, $X\to W$ means $X$ is a direct cause of $W$; a **fork** is $A\leftarrow L\to Y$ (confounder $L$), a **chain** $A\to M\to Y$ (mediator $M$), a **collider** $A\to C\leftarrow Y$. In linear DAGs, $B_{jk}$ is the coefficient on the edge $k\to j$ and $\omega_j$ the variance of node $j$'s own noise. $e(L)=P(A=1\mid L)$ is the **propensity score**. For finite populations of $N$ units, $S^2$ denotes a variance with divisor $N-1$.

## 0. Where this sits

Every model in ST-1 to ST-8 estimates features of a joint distribution: means, conditional means, regression coefficients, spectra. The questions people bring to data are usually about what would happen if something were *changed*: does this drug work, does this policy help, does this design change improve the product. Those are questions about a different distribution, the one after intervention. This module makes the gap precise and says exactly what buys it closed: randomization, or assumptions you state as a DAG and defend. You need F4's conditional probability and Simpson's paradox (F4 Trap 4), ST-2's permutation test (ST-2 Theorem 3.6) and power formula (ST-2 §3.9), ST-3's Frisch–Waugh–Lovell theorem (ST-3 Theorem 3.4) and omitted-variable bias (ST-3 Theorem 3.12), and ST-4's logistic regression. Afterwards, re-read every coefficient table in ST-3 and every feature importance in 04-machine-learning (ML-3, ML-5) with one question: what DAG would make this causal?

## 1. The problem that forced this idea

Agricultural field trials in the early twentieth century compared fertilizers on plots that differed in soil, drainage and shade. Experimenters arranged plots systematically (alternating treatments, "chessboard" layouts) to balance the soil, and argued endlessly about whether the balance worked. In *Statistical Methods for Research Workers* (1925) and "The arrangement of field experiments" (1926), R. A. Fisher at Rothamsted gave the answer that ended the argument: assign treatments **at random**. Randomization does not balance any particular plot; it makes the assignment independent of everything about the plots, measured or not, and it supplies the probability distribution that justifies the test. His 1935 *Design of Experiments* opens with the lady tasting tea, the first randomization test.

Two years before Fisher's book, Jerzy Neyman (writing as Splawa-Neyman, 1923, in Polish) had written down, for the same agricultural problem, one "potential yield" per plot per variety, of which only one is ever observed. The notation sat unused for half a century until Donald Rubin (1974) revived it for observational studies in education and medicine, where nobody randomized. Paul Holland (1986) named the obstacle the **fundamental problem of causal inference**: for each unit you see one potential outcome and never the other.

A separate line started with Sewall Wright, a geneticist studying coat colour in guinea pigs. In "Correlation and causation" (1921) he drew arrows between variables and showed that a correlation equals a sum over paths of products of path coefficients. Statisticians mostly ignored path diagrams for decades, while economists built simultaneous-equation models with instruments (Philip Wright, Sewall's father, used an instrumental variable in 1928 to estimate supply and demand for flaxseed). Judea Pearl, working in AI on probabilistic reasoning (1988), turned the drawings into a calculus: directed acyclic graphs, d-separation, the do-operator, the backdoor criterion (1993; the 1995 *Biometrika* paper), collected in *Causality* (2000).

The hardest test case was smoking. Richard Doll and Austin Bradford Hill's 1950 case-control study found smokers heavily over-represented among lung cancer patients. Fisher himself argued that a genotype could cause both smoking and cancer. Nobody could randomize people to smoke for thirty years. In 1965 Hill set out nine considerations for judging whether an association is causal; §3.12 lists them. The smoking case was won by the convergence of many imperfect studies, each with different weaknesses.

The forcing fact is this. An association is a property of one distribution, $P(Y\mid A)$. A causal effect compares two distributions, $P(Y\mid do(A=1))$ and $P(Y\mid do(A=0))$, and observational data come from neither. The difference between them is not a slogan; it is a term you can write down and compute (§3.1). Closing it takes randomization or untestable assumptions, stated out loud.

## 2. Intuition

**Two worlds, one observed.** A patient takes the drug and recovers. Did the drug help? The question compares the world where she took it with the world where she didn't. You only see one. So "the effect for this patient" is a missing-data problem with the other half always missing. You cannot fix this with more data on the same patient; you can only fix it by comparing *different* patients and arguing that the comparison is fair.

**When is comparing groups fair?** Suppose doctors give the drug to the sickest patients. The treated recover less often, not because the drug hurts but because they were sicker to begin with. The untreated group is a bad stand-in for "the treated, had they not been treated." The gap between what the untreated did and what the treated *would have done* untreated is the **selection bias**. Guess the decomposition: the observed difference in means is the effect on the treated plus how differently the two groups would have fared with no treatment at all.

**Randomization makes the groups interchangeable.** Flip a coin for each patient. Now the treated and untreated are two random halves of the same population. Each half's untreated outcome distribution is the same, so the selection bias is zero in expectation and the observed difference estimates the average effect. Coin flips achieve what no list of measured covariates can: balance on everything, including what you didn't think to measure. Guess what randomization buys formally: treatment is independent of the potential outcomes.

**A confounder is a common cause.** Severity causes treatment (doctors' choice) and causes recovery. That is a fork, $A\leftarrow L\to Y$. If you compare treated and untreated *within* each severity level, the fork is cut, and averaging those within-level comparisons over the population's severity distribution gives the effect. Guess the formula: $\mathbb E[Y^a]=\sum_l\mathbb E[Y\mid A=a,L=l]P(L=l)$.

**A DAG tells you what to adjust for, and what not to.** Three shapes cover every local situation:

```
  fork (confounder)       chain (mediator)        collider
       L                                              C
      / \                 A --> M --> Y              / \
     v   v                                          ^   ^
    A     Y                                        A     Y
  A, Y dependent;        A, Y dependent;        A, Y independent;
  given L: independent   given M: independent   given C: DEPENDENT
```

The first two behave the way intuition says: conditioning on the middle variable blocks the flow. The third is the one that hurts. Talent and looks are independent in the population. Film stars are selected for having enough of one or the other. Among film stars, the less talented are better looking: conditioning on the collider "became famous" created a negative correlation from nothing. Guess the consequence for regression: "control for everything you measured" can *create* bias if one of those variables is a collider, and it removes part of the effect if one is a mediator.

**Regression "controls for" in the FWL sense.** ST-3 Theorem 3.4 says the coefficient on $A$ in a regression on $A$ and $L$ is the slope of $Y$-residualized-on-$L$ against $A$-residualized-on-$L$. That is an operation on a joint distribution. It is the right operation exactly when $L$ blocks the confounding paths and opens none, which the DAG decides and the data cannot.

## 3. Formal development

### 3.1 Potential outcomes and the selection-bias decomposition

Each unit $i$ has two potential outcomes $Y_i^1,Y_i^0$. We observe $A_i$ and
$$Y_i=A_iY_i^1+(1-A_i)Y_i^0\qquad\text{(consistency)}.\tag{3.1}$$
Consistency says the treatment is well defined: "taking the drug" means one thing, and one unit's treatment does not change another's outcome (no interference; together with consistency, Rubin's SUTVA). Vaccines (herd immunity) and network effects violate the second part.

▶ **Try it first.** Write $\mathbb E[Y\mid A=1]-\mathbb E[Y\mid A=0]$ using (3.1), then add and subtract $\mathbb E[Y^0\mid A=1]$.

1. By (3.1), on $\{A=1\}$ we have $Y=Y^1$ and on $\{A=0\}$ we have $Y=Y^0$. So the association is $\mathbb E[Y^1\mid A=1]-\mathbb E[Y^0\mid A=0]$.
2. Add and subtract $\mathbb E[Y^0\mid A=1]$, the mean untreated outcome of the treated, which is never observed:
$$\underbrace{\mathbb E[Y\mid A=1]-\mathbb E[Y\mid A=0]}_{\text{association}}=\underbrace{\mathbb E[Y^1-Y^0\mid A=1]}_{\text{ATT}}+\underbrace{\mathbb E[Y^0\mid A=1]-\mathbb E[Y^0\mid A=0]}_{\text{selection bias}}.\tag{3.2}$$
3. The ATT itself differs from the ATE when effects are heterogeneous and the treated are not typical: $\mathrm{ATT}=\mathrm{ATE}+(\mathrm{ATT}-\mathrm{ATE})$.

■ **What just happened.** The observed comparison is the causal quantity plus a term built from an unobservable mean. Everything in this module is a way of making that term zero (randomization), estimating around it (adjustment, IV, RD, DiD), or bounding it (sensitivity analysis). Example 4.2 makes it a number.

### 3.2 Randomized experiments

**Theorem 3.1 (Randomization identifies the ATE).** If $A\perp(Y^0,Y^1)$ (exchangeability) and $0<P(A=1)<1$, then the association equals the ATE.

*Proof.*
1. $\mathbb E[Y\mid A=a]=\mathbb E[Y^a\mid A=a]$ by consistency (3.1).
2. $\mathbb E[Y^a\mid A=a]=\mathbb E[Y^a]$ because $Y^a$ is independent of $A$; positivity makes the conditioning event non-null.
3. Subtract the $a=0$ case from the $a=1$ case. ∎

In (3.2), exchangeability makes the selection bias $\mathbb E[Y^0\mid A=1]-\mathbb E[Y^0\mid A=0]=0$ and the ATT equal to the ATE. A coin, a random-number generator, or a lottery produces $A\perp(Y^0,Y^1)$ by construction; nothing else does.
$$Y^a\perp A\quad\Longrightarrow\quad\mathbb E[Y\mid A=1]-\mathbb E[Y\mid A=0]=\mathrm{ATE}.\tag{3.3}$$

**Fisher's randomization test.** Take $N$ units, $n_1$ assigned to treatment by drawing one of the $\binom N{n_1}$ subsets uniformly. Fisher's **sharp null** is
$$H_0:\ Y_i^1=Y_i^0\ \text{for every } i.\tag{3.4}$$
Under (3.4), each unit's outcome is the same whatever its assignment, so the observed vector $y$ is the outcome vector under *every* assignment. For any statistic $T(w,y)$ of an assignment $w$ (say, difference in means) we can compute its value under all $\binom N{n_1}$ assignments. The exact p-value is
$$p=\frac{\#\{w:\ T(w,y)\ge T(w_{\rm obs},y)\}}{\binom N{n_1}}.\tag{3.5}$$
It is exact because under $H_0$ the observed assignment is a uniform draw from that set, the same argument as ST-2 Theorem 3.6: the permutation test *is* Fisher's randomization test, with the randomization the experimenter did standing in for the exchangeability ST-2 assumed. The test needs no model for $Y$, no Normality, no large $N$. Its limitation is the null: (3.4) says no unit is affected, which is stronger than "the average effect is zero."

**Neyman's variance.** Neyman's target is the finite-population ATE $\tau=\frac1N\sum_i(Y_i^1-Y_i^0)$ and the estimator $\hat\tau=\bar Y_1-\bar Y_0$ (mean among treated minus mean among controls). The randomness is only in the assignment.

**Theorem 3.2 (Neyman, 1923).** Under complete randomization with $n_1$ treated and $n_0=N-n_1$ controls, $\mathbb E\hat\tau=\tau$ and
$$\operatorname{Var}(\hat\tau)=\frac{S_1^2}{n_1}+\frac{S_0^2}{n_0}-\frac{S_\tau^2}{N},\tag{3.6}$$
where $S_1^2,S_0^2,S_\tau^2$ are the finite-population variances (divisor $N-1$) of $Y_i^1$, $Y_i^0$ and $Y_i^1-Y_i^0$.

▶ **Try it first.** Write $\hat\tau$ as $\sum_iW_ic_i+\text{const}$ with $W_i$ the treatment indicator. What are $\operatorname{Var}W_i$ and $\operatorname{Cov}(W_i,W_j)$ when exactly $n_1$ of $N$ are treated?

1. With $p=n_1/N$: $\hat\tau=\frac1{n_1}\sum W_iY_i^1-\frac1{n_0}\sum(1-W_i)Y_i^0=\sum_iW_ic_i-\frac1{n_0}\sum_iY_i^0$, where $c_i=Y_i^1/n_1+Y_i^0/n_0$.
2. Unbiasedness: $\mathbb EW_i=p$, so $\mathbb E\hat\tau=\frac{p}{n_1}\sum Y_i^1+\frac p{n_0}\sum Y_i^0-\frac1{n_0}\sum Y_i^0=\bar Y^1+(\frac{n_1}{Nn_0}-\frac1{n_0})\sum Y_i^0=\bar Y^1-\bar Y^0=\tau$.
3. $\operatorname{Var}W_i=p(1-p)$. Since $\sum_iW_i=n_1$ is fixed, $0=\operatorname{Var}\sum W_i=Np(1-p)+N(N-1)\operatorname{Cov}(W_i,W_j)$, so $\operatorname{Cov}(W_i,W_j)=-p(1-p)/(N-1)$.
4. $\operatorname{Var}\sum W_ic_i=p(1-p)\big[\sum c_i^2-\frac1{N-1}\sum_{i\ne j}c_ic_j\big]$. Using $\sum_{i\neq j}c_ic_j=(\sum c_i)^2-\sum c_i^2$ this is $p(1-p)\frac N{N-1}\sum(c_i-\bar c)^2=\frac{n_1n_0}{N}S_c^2$.
5. $S_c^2=S_1^2/n_1^2+S_0^2/n_0^2+2S_{10}/(n_1n_0)$, with $S_{10}$ the finite-population covariance of $Y^1,Y^0$. Multiply by $n_1n_0/N$: $\operatorname{Var}\hat\tau=\frac{n_0}{Nn_1}S_1^2+\frac{n_1}{Nn_0}S_0^2+\frac2NS_{10}$.
6. Write $\frac{n_0}{Nn_1}=\frac1{n_1}-\frac1N$ and $\frac{n_1}{Nn_0}=\frac1{n_0}-\frac1N$, and use $S_\tau^2=S_1^2+S_0^2-2S_{10}$. This gives (3.6). ∎

■ **What just happened.** The first two terms are the textbook two-sample variance. The third is negative and involves $S_\tau^2$, the variance of individual effects, which needs both potential outcomes for each unit and so can never be estimated. Dropping it gives the **conservative** estimator $s_1^2/n_1+s_0^2/n_0$, unbiased exactly when the effect is constant. The ordinary Welch standard error is therefore valid for randomized experiments without any sampling-from-a-population story, and slightly too large.

**Blocking.** If you know before the experiment that units fall into groups with different outcome levels (sick and healthy, two machines, two plots), randomize separately within each block $b$ of size $N_b$ and combine $\hat\tau_B=\sum_b(N_b/N)\hat\tau_b$. Blocks are randomized independently, so
$$\operatorname{Var}\hat\tau_B=\sum_b\Big(\frac{N_b}N\Big)^2V_b,\tag{3.7}$$
with $V_b$ given by (3.6) inside block $b$. Why this helps, in the simplest case: constant effect, two blocks of equal size, half treated in each. Then (3.6) gives $\operatorname{Var}\hat\tau=4S^2/N$ under complete randomization and (3.7) gives $\operatorname{Var}\hat\tau_B=\frac14(8S_1^2/N+8S_2^2/N)=4\bar S_w^2/N$, where $\bar S_w^2=\frac12(S_{(1)}^2+S_{(2)}^2)$ is the mean within-block variance of $Y^0$. The ANOVA identity (ST-2 Example 4.7's decomposition) for two equal blocks with means $\mu_1,\mu_2$ is
$$(N-1)S^2=(N-2)\bar S_w^2+\tfrac N4(\mu_1-\mu_2)^2.\tag{3.8}$$
So blocking removes the between-block term $\frac N4(\mu_1-\mu_2)^2$ from the variance. The more the blocks differ, the bigger the gain. Fisher's slogan: block what you can, randomize what you cannot. Example 4.2's population shows a threefold reduction.

### 3.3 Factorial designs, design of experiments, A/B tests

**The $2^k$ factorial.** $k$ factors, each at a low ($-1$) and high ($+1$) level, every combination run ($2^k$ cells, $r$ replicates each, $N=r2^k$ runs). Let $x_{ji}=\pm1$ be factor $j$'s level in cell $i$ and $\bar y_i$ the cell mean. The **main effect** of factor $j$ is the mean response at its high level minus the mean at its low level, averaging over everything else; the **interaction** of $j$ and $k$ is half the difference between $j$'s effect at high $k$ and at low $k$. Both are contrasts:
$$\widehat{\mathrm{eff}}_S=\frac1{2^{k-1}}\sum_{i=1}^{2^k}\Big(\prod_{j\in S}x_{ji}\Big)\bar y_i,\qquad S\subseteq\{1,\dots,k\}.\tag{3.9}$$
The columns $\prod_{j\in S}x_{j\cdot}$ are mutually orthogonal $\pm1$ vectors (LA-5), which is why every effect is estimated with all $N$ runs. Each contrast has variance $2^k\sigma^2/(r\,4^{k-1})=4\sigma^2/N$. Compare one-factor-at-a-time (OFAT): to estimate three main effects, each as a difference of two averages of $m$ runs with the same variance $4\sigma^2/8=\sigma^2/2$ as a $2^3$ factorial with $N=8$, you need $2\sigma^2/m=\sigma^2/2$, so $m=4$, and $4+3\times4=16$ runs. OFAT costs twice as much and cannot see interactions at all. A regression of $y$ on the $\pm1$ columns gives coefficients equal to half the effects, because the regression slope is per unit of $x$ and the effect spans two units.

This is engineering **design of experiments (DOE)**. Fractional factorials run a $2^{-p}$ fraction of the cells and deliberately confound high-order interactions with each other (Box, Hunter & Hunter ch. 6). Taguchi methods add a noise-factor array to design products robust to uncontrollable variation; response-surface methods (Box & Wilson, 1951) fit quadratic models around an optimum. For an ME, these are how a process is tuned with a fixed budget of runs.

**A/B tests** are randomized experiments on users. Sample size comes from ST-2 §3.9: to detect a conversion rate change from $p_0=0.10$ to $p_1=0.11$ with a two-sided 5% test and 80% power, $n=(z_{0.975}+z_{0.8})^2[p_0(1-p_0)+p_1(1-p_1)]/\Delta^2=14{,}749$ per arm (the project computes it). The classic failure is **peeking**: checking the p-value as data arrive and stopping at the first $p<0.05$. With 10 equally spaced looks at an A/A test (no difference), the project finds a false-positive rate of 0.192, not 0.05. This is ST-2 exercise C3's optional stopping. Fix it with a fixed horizon, a group-sequential boundary, or an always-valid sequential test.

**Noncompliance.** People assigned to treatment may not take it. Comparing by treatment *received* breaks randomization (who complies is not random). Comparing by treatment *assigned* is the **intention-to-treat (ITT)** effect: unbiased for the effect of the assignment, by Theorem 3.1 applied to the assignment. To get the effect of the treatment itself among those who comply, treat assignment as an instrument (§3.8).

### 3.4 Observational studies: conditional exchangeability, standardization, IPW

Without randomization, the best hope is that within strata of measured covariates $L$, treatment is as good as random:
$$Y^a\perp A\mid L\quad\text{(conditional exchangeability, no unmeasured confounding)},\qquad 0<P(A=a\mid L=l)\ \text{for all } l,a\quad\text{(positivity).}$$

**Theorem 3.3 (Standardization, the g-formula).** Under consistency, conditional exchangeability and positivity,
$$\mathbb E[Y^a]=\sum_l\mathbb E[Y\mid A=a,L=l]\,P(L=l).\tag{3.10}$$

▶ **Try it first.** Condition $\mathbb E[Y^a]$ on $L$ and use the three assumptions one at a time.

1. Law of total expectation: $\mathbb E[Y^a]=\sum_l\mathbb E[Y^a\mid L=l]P(L=l)$.
2. Conditional exchangeability: $\mathbb E[Y^a\mid L=l]=\mathbb E[Y^a\mid A=a,L=l]$. Positivity makes the right side defined.
3. Consistency: on $\{A=a\}$, $Y^a=Y$. So $\mathbb E[Y^a\mid A=a,L=l]=\mathbb E[Y\mid A=a,L=l]$. ∎

**Theorem 3.4 (Inverse probability weighting).** Under the same assumptions,
$$\mathbb E[Y^a]=\mathbb E\!\left[\frac{\mathbb 1(A=a)\,Y}{P(A=a\mid L)}\right].\tag{3.11}$$

*Proof.*
1. Consistency: $\mathbb 1(A=a)Y=\mathbb 1(A=a)Y^a$.
2. Condition on $L$: $\mathbb E\big[\frac{\mathbb 1(A=a)Y^a}{P(A=a\mid L)}\big]=\mathbb E\big[\frac{\mathbb E[\mathbb 1(A=a)Y^a\mid L]}{P(A=a\mid L)}\big]$.
3. Conditional exchangeability factorizes the inner expectation: $\mathbb E[\mathbb 1(A=a)Y^a\mid L]=P(A=a\mid L)\,\mathbb E[Y^a\mid L]$.
4. The propensity cancels (positivity: it is not zero), leaving $\mathbb E[\mathbb E[Y^a\mid L]]=\mathbb E[Y^a]$. ∎

■ **What just happened.** Standardization models the outcome; IPW models the treatment. IPW builds a **pseudo-population**: each treated unit with propensity $0.25$ stands for 4 people, so in the reweighted population treatment is independent of $L$, as if randomized. If $P(A=a\mid L=l)$ is near 0 for some $l$, the weights explode: practical positivity violations are why IPW estimates can be wildly variable.

**Equivalence in the saturated case.** With discrete $L$ and the empirical probabilities plugged in ($\hat P(l)=n_l/n$, $\hat P(A=a\mid l)=n_{la}/n_l$), the IPW estimate is
$$\frac1n\sum_l\sum_{i:\,L_i=l,A_i=a}\frac{Y_i}{n_{la}/n_l}=\sum_l\frac{n_l}{n}\cdot\frac1{n_{la}}\sum_{i:\,L_i=l,A_i=a}Y_i=\sum_l\hat P(l)\,\hat{\mathbb E}[Y\mid A=a,L=l],$$
which is the plug-in g-formula. The two are the same number (Example 4.4). They differ only once you fit parametric models for the outcome or the propensity, and then each is right only if its model is. A **doubly robust** (augmented IPW) estimator combines both and is consistent if either model is correct (Robins, Rotnitzky & Zhao 1994; Hernán & Robins ch. 13).

**Theorem 3.5 (Propensity score, Rosenbaum & Rubin 1983).** Let $e(L)=P(A=1\mid L)$. (i) $A\perp L\mid e(L)$ (balance). (ii) If $Y^a\perp A\mid L$ then $Y^a\perp A\mid e(L)$.

*Proof.* (i) $P(A=1\mid L,e(L))=P(A=1\mid L)=e(L)$, since $e(L)$ is a function of $L$. Also $P(A=1\mid e(L))=\mathbb E[P(A=1\mid L)\mid e(L)]=\mathbb E[e(L)\mid e(L)]=e(L)$. The two agree, so given $e(L)$, $L$ carries no further information about $A$. (ii) $P(A=1\mid Y^a,e(L))=\mathbb E[P(A=1\mid Y^a,L)\mid Y^a,e(L)]=\mathbb E[e(L)\mid Y^a,e(L)]=e(L)$, using conditional exchangeability in the middle step. This does not depend on $Y^a$. ∎

So a high-dimensional $L$ can be replaced by a scalar: adjust for, stratify on, or match on $e(L)$. **Matching** pairs each treated unit with one or more controls with similar $L$ or $e(L)$ and averages the differences; it estimates the ATT, its bias comes from imperfect matches, and its advantage is transparency (you can inspect balance before seeing outcomes). The propensity score theorem does not create exchangeability; it only compresses the $L$ you already had. "We used propensity scores" says nothing about unmeasured confounding.

### 3.5 DAGs, the do-operator, and d-separation

A **DAG** $G$ on variables $V_1,\dots,V_p$ is a set of arrows with no directed cycles. A **path** is any sequence of distinct nodes joined by edges, ignoring direction; a **directed path** follows the arrows. $\mathrm{pa}(j)$ are the parents of $V_j$; descendants are nodes reachable by directed paths. A distribution is **Markov** to $G$ if
$$p(v_1,\dots,v_p)=\prod_jp(v_j\mid \mathrm{pa}_j).\tag{3.12}$$
A causal DAG says more: each factor is an autonomous mechanism, and an intervention $do(A=a)$ replaces $A$'s mechanism by the constant $a$ and leaves the others alone (the **truncated factorization**):
$$p(v\mid do(A=a))=\prod_{j:\,V_j\ne A}p(v_j\mid \mathrm{pa}_j)\Big|_{A=a}.\tag{3.13}$$
In a linear SCM, $V_j=\sum_kB_{jk}V_k+\varepsilon_j$, intervention deletes row $A$ of $B$ and sets $V_A=a$; the project's `SCM.sample(do=...)` does exactly that. Potential outcomes and do are two notations for one object: $Y^a$ is $Y$ in the intervened model.

**The three structures, derived.** Use linear Gaussian variables with unit noise variances and F5's conditional covariance: for jointly Gaussian $X,W,C$,
$$\operatorname{Cov}(X,W\mid C)=\Sigma_{XW}-\Sigma_{XC}\Sigma_{WC}/\Sigma_{CC}.\tag{3.14}$$

▶ **Try it first.** For each structure, compute $\Sigma_{AY}$ and then (3.14) with the middle node as $C$.

1. **Fork** $A\leftarrow L\to Y$, $A=\lambda L+\varepsilon_A$, $Y=\mu L+\varepsilon_Y$. $\Sigma_{AY}=\lambda\mu$, $\Sigma_{AL}=\lambda$, $\Sigma_{YL}=\mu$, $\Sigma_{LL}=1$. Conditional covariance: $\lambda\mu-\lambda\mu=0$. Dependent; independent given $L$.
2. **Chain** $A\to M\to Y$, $M=aA+\varepsilon_M$, $Y=bM+\varepsilon_Y$. $\Sigma_{AY}=ab$, $\Sigma_{AM}=a$, $\Sigma_{MM}=a^2+1$, $\Sigma_{MY}=b(a^2+1)$. Conditional: $ab-a\cdot b(a^2+1)/(a^2+1)=0$. Dependent; independent given $M$.
3. **Collider** $A\to C\leftarrow Y$, $C=\gamma A+\delta Y+\varepsilon_C$, $A,Y$ independent. $\Sigma_{AY}=0$, $\Sigma_{AC}=\gamma$, $\Sigma_{YC}=\delta$, $\Sigma_{CC}=\gamma^2+\delta^2+1$.
$$\operatorname{Cov}(A,Y\mid C)=-\frac{\gamma\delta}{\gamma^2+\delta^2+1}\neq0.\tag{3.15}$$
Independent; *dependent* given $C$, with sign opposite to $\gamma\delta$.
4. Conditioning on a **descendant** of a collider (a noisy copy $C'=C+\varepsilon'$) gives (3.15) with $\Sigma_{C'C'}$ in the denominator: smaller but still nonzero. A descendant of a collider partially opens it.

■ **What just happened.** Information flows through forks and chains unless you condition on the middle; it is blocked at a collider unless you condition on it or a descendant. (3.15) is Berkson's paradox (Berkson 1946, hospital patients): selecting on a common effect induces association between its causes. The project's "famous" sample shows correlation 0.000 in the population turning into $-0.683$ among the selected.

**Definition 3.6 (d-separation).** A path between $X$ and $W$ is **blocked** by a set $S$ if it contains a node $m$ such that either (a) $m$ is a non-collider on the path (fork or chain) and $m\in S$, or (b) $m$ is a collider on the path and neither $m$ nor any descendant of $m$ is in $S$. $X$ and $W$ are **d-separated** by $S$ if every path between them is blocked.

The **global Markov property** says: if $X$ and $W$ are d-separated by $S$ in $G$, then $X\perp W\mid S$ in every distribution Markov to $G$. Conversely, if they are not d-separated, they are dependent given $S$ in almost every such distribution (all except a measure-zero set of parameter values where effects cancel; "faithfulness"). The general proof (Lauritzen 1996, §3.2; Pearl 2009, §1.2.3) is a careful induction on the factorization; the three cases above are its local building blocks.

### 3.6 The backdoor criterion and the adjustment traps

A **backdoor path** from $A$ to $Y$ is a path that begins with an arrow *into* $A$. Such paths carry association that is not effect.

**Theorem 3.7 (Backdoor criterion, Pearl 1993).** If a set $L$ (i) contains no descendant of $A$ and (ii) blocks every backdoor path from $A$ to $Y$, then $Y^a\perp A\mid L$, and
$$P(y\mid do(A=a))=\sum_lP(y\mid A=a,L=l)P(L=l).$$

*Proof in the simplest case* (DAG $L\to A$, $L\to Y$, $A\to Y$, so $L$ satisfies the criterion).
1. Markov factorization (3.12): $p(l,a,y)=p(l)\,p(a\mid l)\,p(y\mid a,l)$.
2. Truncated factorization (3.13): $p(l,y\mid do(A=a))=p(l)\,p(y\mid a,l)$.
3. Sum over $l$: $p(y\mid do(a))=\sum_lp(y\mid a,l)p(l)$, the g-formula (3.10). ∎

The general proof (Pearl 2009, §3.3.1) adds an arrow $F_A\to A$ for the intervention and shows that conditions (i)–(ii) make $Y$ and $F_A$ d-separated given $\{A\}\cup L$. Condition (i) matters: a descendant of $A$ is either a mediator (adjusting removes effect) or downstream of a collider path (adjusting opens it).

**Applying it to five DAGs.**

| DAG | Correct set | Wrong choice and why |
|---|---|---|
| Fork: $L\to A$, $L\to Y$, $A\to Y$ | $\{L\}$ | $\varnothing$: path $A\leftarrow L\to Y$ open |
| Chain: $A\to M\to Y$, $A\to Y$ | $\varnothing$ (total effect) | $\{M\}$: descendant of $A$, removes the indirect effect |
| Collider: $A\to C\leftarrow Y$, $A\to Y$ | $\varnothing$ | $\{C\}$: descendant of $A$, opens $A\to C\leftarrow Y$ |
| M-bias: $A\leftarrow U_1\to C\leftarrow U_2\to Y$, $A\to Y$ | $\varnothing$ | $\{C\}$: pre-treatment, but a collider; opens the path |
| Example 4.6's five-node DAG | $\{Z_1,Z_3\}$, $\{Z_2,Z_3\}$, $\{Z_1,Z_2,Z_3\}$ | $\{Z_3\}$ alone blocks two paths and opens one |

**M-bias and the "adjust for everything" fallacy.** In the M-structure, $C$ is measured before treatment, correlated with both $A$ and $Y$, and looks like a textbook confounder. Adjusting for it opens $A\leftarrow U_1\to C\leftarrow U_2\to Y$. With all coefficients 1 and effect 0.5, the project finds the regression on $A$ alone gives 0.5000 and on $(A,C)$ gives 0.3000. The rule "a pre-treatment variable associated with treatment and outcome is a confounder" is false. In practice M-bias is usually smaller than the confounding a variable removes (Ding & Miratrix 2015), but the logic is the point: whether to adjust is a property of the DAG, not of correlations.

**Front door** (stated). If $A\to M\to Y$ with $M$ the only route, $A\leftarrow U\to Y$ confounded, and $M$ unconfounded with $A$ and with $Y$ given $A$, then $P(y\mid do(a))=\sum_mP(m\mid a)\sum_{a'}P(y\mid m,a')P(a')$ (Pearl 2009, §3.3.2). It identifies an effect with an unmeasured confounder, at the cost of strong assumptions about the mediator.

**Mediation and the post-treatment fallacy.** Linear chain with a direct arrow: $M=aA+\varepsilon_M$, $Y=cA+bM+\varepsilon_Y$. The total effect of $A$ is $c+ab$ (substitute $M$). Regressing $Y$ on $A$ and $M$ gives $c$: by FWL the coefficient on $A$ uses only the part of $A$ unrelated to $M$, so it is the **controlled direct effect**,
$$\beta_{A\mid M}=c=\text{total}-ab.\tag{3.16}$$
If the question is "does $A$ work," adjusting for $M$ answers a different question. Worse, mediators are usually confounded with the outcome. Add $U\to M$ (coefficient $h$) and $U\to Y$ (coefficient $g$), $U$ independent of $A$, unit variances. Then $M$ is a collider on $A\to M\leftarrow U\to Y$.

▶ **Try it first.** Compute $\Sigma_{AM},\Sigma_{MM},\Sigma_{AY},\Sigma_{MY}$, then the coefficient on $A$ in the regression of $Y$ on $(A,M)$.

1. $\Sigma_{AA}=1$, $\Sigma_{AM}=a$, $\Sigma_{MM}=a^2+h^2+1$, $\Sigma_{AY}=c+ab$, $\Sigma_{MY}=ac+b\Sigma_{MM}+gh$.
2. Two-regressor formula: $\beta_A=\dfrac{\Sigma_{MM}\Sigma_{AY}-\Sigma_{AM}\Sigma_{MY}}{\Sigma_{AA}\Sigma_{MM}-\Sigma_{AM}^2}$. The denominator is $h^2+1$.
3. Numerator: $\Sigma_{MM}(c+ab)-a(ac+b\Sigma_{MM}+gh)=c(\Sigma_{MM}-a^2)-agh=c(h^2+1)-agh$.
$$\beta_{A\mid M}=c-\frac{agh}{1+h^2}.\tag{3.17}$$

■ **What just happened.** Even the direct effect is biased once the mediator has its own confounder: conditioning on $M$ opens the path through $U$. With $a=0.8$, $b=0.5$, $c=0.3$, $g=h=1$ the project gets $0.3-0.4=-0.10$: the direct effect has the wrong sign. Mediation analysis needs no-confounding assumptions for $M$ as well as for $A$.

**Selection bias is collider bias.** If inclusion in the data $S$ depends on treatment and outcome (hospitalized patients, survivors, published studies, customers who did not churn), analysing only $S=1$ conditions on a collider $A\to S\leftarrow Y$. Survivorship bias is the same graph. Exercise B4 computes a case.

**Measurement error.** Adjusting for a noisy version $L^*$ of a confounder $L$ blocks the fork only partially (in a linear model, by the reliability ratio $\operatorname{Var}L/\operatorname{Var}L^*$), leaving residual confounding that does not shrink with $n$. Error in $A$ attenuates the effect toward zero (ST-3's errors-in-variables).

### 3.7 Regression and causation

A coefficient on $A$ in an OLS regression of $Y$ on $(A,L)$ estimates a causal effect when all of these hold:
1. **$L$ satisfies the backdoor criterion.** Otherwise there is open confounding or collider bias (§3.6).
2. **$L$ contains no descendants of $A$.** No mediators, no post-treatment colliders ("bad controls," Angrist & Pischke §3.2.3).
3. **The functional form is right** for $\mathbb E[Y\mid A,L]$, or else you accept the coefficient as an approximation. Linearity in $L$ that is wrong leaves confounding through the misfit.
4. **Effect homogeneity, or a reinterpretation.** If effects vary with $L$, the OLS coefficient is a weighted average of stratum effects with weights proportional to $\operatorname{Var}(A\mid L)$, not population shares (Angrist 1998). That is neither the ATE nor the ATT. Standardization (3.10) targets the ATE directly.

**The omitted-variable bias formula is the confounding formula.** ST-3 Theorem 3.12: omitting $L$ from a regression with true coefficient $b$ on $L$ shifts the coefficient on $A$ by $b\cdot\gamma$, $\gamma$ the slope of $L$ on $A$. In the fork $L\to A$ (coefficient $a$), $L\to Y$ ($b$), $A\to Y$ ($c$) with unit variances, $\gamma=\operatorname{Cov}(A,L)/\operatorname{Var}A=a/(a^2+1)$, so the short coefficient is $c+ab/(a^2+1)$. "Controlling for" $L$ via FWL removes exactly this.

**Wright's path rules.** Write the linear SCM as $V=BV+\varepsilon$ with $\operatorname{Cov}\varepsilon=\Omega=\operatorname{diag}(\omega_j)$.

▶ **Try it first.** Solve for $V$ and write $\Sigma=\operatorname{Cov}V$. Then expand $(I-B)^{-1}$ as a geometric series. Why does it terminate?

1. $V=(I-B)^{-1}\varepsilon$, so
$$\Sigma=(I-B)^{-1}\,\Omega\,(I-B)^{-\top}.\tag{3.18}$$
2. Order the nodes topologically. Then $B$ is strictly lower triangular, hence nilpotent ($B^p=0$), and $(I-B)^{-1}=I+B+B^2+\dots+B^{p-1}$ exactly.
3. $(B^k)_{it}=\sum B_{i,j_{k-1}}\cdots B_{j_1,t}$ over all node sequences, and a term is nonzero only if $t\to j_1\to\dots\to i$ is a directed path. So $[(I-B)^{-1}]_{it}$ is the sum over directed paths from $t$ to $i$ of the products of their coefficients (the empty path when $t=i$ contributes 1).
4. Entry $(i,j)$ of (3.18) is $\sum_t[(I-B)^{-1}]_{it}\,\omega_t\,[(I-B)^{-1}]_{jt}$. A pair of directed paths $t\rightsquigarrow i$, $t\rightsquigarrow j$ from a common source $t$ is a **trek** with top $t$:
$$\operatorname{Cov}(V_i,V_j)=\sum_{\text{treks } i\leftsquigarrow t\rightsquigarrow j}\omega_t\prod_{\text{edges on both sides}}B.\tag{3.19}$$

■ **What just happened.** Covariance is a sum over treks, and a trek is exactly a path with no collider: it goes up from $i$ to a top and down to $j$. Unconditional covariance never flows through a collider. Two cases, unit $\omega$'s: **3-node fork** $A\leftarrow L\to Y$ with coefficients $a,b$: one trek, top $L$, $\operatorname{Cov}(A,Y)=ab$. **3-node collider** $A\to C\leftarrow Y$: $A$ and $Y$ share no ancestor, no trek, $\operatorname{Cov}(A,Y)=0$, while $\operatorname{Cov}(A,C)=\gamma$ and $\operatorname{Cov}(Y,C)=\delta$ each have one trek. Wright's 1921 rule was the standardized version. The project's `trek_cov` enumerates treks and matches (3.18) and the sample covariance.

**The sign of collider bias.** For the pure collider, (3.15) and the regression coefficient on $A$ in $Y\sim A+C$ both have the sign of $-\gamma\delta$. Two positive causes of a selected outcome look like substitutes. Example 4.5 adds a real effect and shows the bad control can flip its sign.

**The Table 2 fallacy** (Westreich & Greenland 2013). Papers report a regression of $Y$ on the exposure $A$ and confounders $L$, then interpret every coefficient causally. But in the DAG $L\to A\to Y$, $L\to Y$, the coefficient on $L$ given $A$ is $L$'s *direct* effect with its mediator $A$ held fixed, not its total effect, and the set that deconfounds $A$ need not deconfound $L$. One regression, one causal question.

### 3.8 Instrumental variables

When the confounder $U$ is unmeasured, adjustment is impossible. An **instrument** $Z$ satisfies:
- **Relevance**: $Z$ affects $A$ ($\operatorname{Cov}(Z,A)\neq0$).
- **Exclusion**: $Z$ affects $Y$ only through $A$.
- **Independence**: $Z$ shares no causes with $Y$ ($Z\perp U$).

```
   Z ---> A ---> Y
          ^      ^
           \    /
             U        (no arrow Z -> Y, no arrow U -> Z)
```

**The Wald estimator.** Take the linear model $Y=\beta A+g(U)+\varepsilon$ with constant effect $\beta$ and binary $Z$.

▶ **Try it first.** Compute $\mathbb E[Y\mid Z=1]-\mathbb E[Y\mid Z=0]$ and see which terms the three assumptions kill.

1. $\mathbb E[Y\mid Z=z]=\beta\,\mathbb E[A\mid Z=z]+\mathbb E[g(U)\mid Z=z]+\mathbb E[\varepsilon\mid Z=z]$. Exclusion means $Z$ does not appear in the equation for $Y$.
2. Independence: $\mathbb E[g(U)\mid Z=z]$ and $\mathbb E[\varepsilon\mid Z=z]$ do not depend on $z$.
3. Subtract and divide (relevance makes the denominator nonzero):
$$\beta=\frac{\mathbb E[Y\mid Z=1]-\mathbb E[Y\mid Z=0]}{\mathbb E[A\mid Z=1]-\mathbb E[A\mid Z=0]}=\frac{\text{reduced form}}{\text{first stage}}.\tag{3.20}$$
For continuous $Z$ the same steps give $\beta=\operatorname{Cov}(Z,Y)/\operatorname{Cov}(Z,A)$, which for binary $Z$ is (3.20).

■ **What just happened.** $Z$ moves $A$ for reasons unrelated to $U$. The induced movement in $Y$, divided by the induced movement in $A$, is the effect. All of $A$'s confounded variation is thrown away, which is why IV estimates are noisy.

**LATE.** With heterogeneous effects, (3.20) estimates something specific. Each unit has potential treatments $A(1),A(0)$ under $Z=1,0$: **compliers** ($A(1)=1,A(0)=0$), **always-takers**, **never-takers**, **defiers** ($A(1)=0,A(0)=1$). Assume $Z\perp(A(0),A(1),Y^0,Y^1)$ and exclusion, so $Y=Y^{A(Z)}$.
1. $\mathbb E[Y\mid Z=1]-\mathbb E[Y\mid Z=0]=\mathbb E[Y^{A(1)}-Y^{A(0)}]$ by independence.
2. $Y^{A(1)}-Y^{A(0)}=(A(1)-A(0))(Y^1-Y^0)$, which is $+(Y^1-Y^0)$ for compliers, $-(Y^1-Y^0)$ for defiers, 0 otherwise.
3. So the numerator is $\pi_c\,\mathbb E[\tau\mid c]-\pi_d\,\mathbb E[\tau\mid d]$ and the first stage is $\pi_c-\pi_d$, with $\tau=Y^1-Y^0$ and $\pi$ the type shares.
4. **Monotonicity** (no defiers, $\pi_d=0$) gives
$$\text{Wald}=\mathbb E[Y^1-Y^0\mid \text{complier}]\quad(\text{LATE; Imbens \& Angrist 1994}).\tag{3.21}$$
IV estimates the effect for the people the instrument moves, who are not identifiable individually and may not be the people a policy would target.

**Two-stage least squares.** General linear model $y=X\beta+\varepsilon$ ($X$ is $n\times k$ and includes $A$ and exogenous controls), instruments $Z$ ($n\times m$, $m\ge k$, controls included as their own instruments) with $\mathbb E[Z^\top\varepsilon]=0$. In the just-identified case $m=k$, the sample analogue of $\mathbb E[Z^\top(y-X\beta)]=0$ is
$$\hat\beta_{\rm IV}=(Z^\top X)^{-1}Z^\top y.\tag{3.22}$$
**2SLS** works for $m\ge k$. Stage 1: regress each column of $X$ on $Z$, $\hat X=P_ZX$ with $P_Z=Z(Z^\top Z)^{-1}Z^\top$ the projection onto the column space of $Z$ (LA-5). Stage 2: regress $y$ on $\hat X$:
$$\hat\beta_{\rm 2SLS}=(\hat X^\top\hat X)^{-1}\hat X^\top y=(X^\top P_ZX)^{-1}X^\top P_Zy,\tag{3.23}$$
using $P_Z^\top P_Z=P_Z$. When $m=k$, $Z^\top X$ is square and invertible, and $(X^\top Z(Z^\top Z)^{-1}Z^\top X)^{-1}=(Z^\top X)^{-1}(Z^\top Z)(X^\top Z)^{-1}$. Multiplying by $X^\top Z(Z^\top Z)^{-1}Z^\top y$ collapses (3.23) to (3.22). **Consistency**: $\hat\beta_{\rm 2SLS}-\beta=(X^\top P_ZX/n)^{-1}(X^\top Z/n)(Z^\top Z/n)^{-1}(Z^\top\varepsilon/n)$ and $Z^\top\varepsilon/n\to0$. **Standard errors**: $\widehat{\operatorname{Var}}=\hat\sigma^2(X^\top P_ZX)^{-1}$ with $\hat\sigma^2$ from residuals $y-X\hat\beta$, **not** from the second-stage residuals $y-\hat X\hat\beta$. Running the two stages by hand and reading OLS output from stage 2 gives the wrong SE.

**Weak instruments.** When $Z$ barely moves $A$, the first stage $P_ZX$ is mostly the projection of noise. In the extreme $\operatorname{Cov}(Z,A)=0$, $\hat X$ is $X$ projected onto a random subspace; the part of $X$ that lands there is correlated with $\varepsilon$ the same way $X$ is, and 2SLS behaves like OLS. Staiger & Stock (1997) and Stock, Wright & Yogo (2002) show that the bias of 2SLS relative to OLS is approximately $1/F$, with $F$ the first-stage F statistic for the excluded instruments. Hence the rule of thumb **$F>10$** keeps the relative bias near 10% or less. Weak instruments also make the usual confidence intervals badly undercover; use Anderson–Rubin intervals. The project plots 2SLS median bias against first-stage $F$.

**Classic instruments, critiqued.**
- *Draft lottery* (Angrist 1990): random Vietnam draft numbers → military service → later earnings. Independence holds by design. Exclusion is arguable: men with low numbers may have stayed in school to get deferments, which affects earnings directly.
- *Distance to college* (Card 1995): growing up near a college → more schooling → earnings. Distance correlates with urban residence and family background, both of which affect earnings: independence is doubtful, and the paper conditions on region and background to defend it.
- *Rainfall* (Miguel, Satyanath & Sergenti 2004): rainfall → economic growth → civil conflict in African countries. Rainfall can affect conflict through other channels (roads, disease, the mobility of fighters), so exclusion is contested.

The pattern: relevance is testable, independence is sometimes guaranteed by design, and exclusion is almost never testable and almost always arguable.

### 3.9 Other designs: regression discontinuity, difference-in-differences, synthetic control

**Regression discontinuity** (Thistlethwaite & Campbell 1960). Treatment is assigned by a threshold on a running variable $R$: $A=\mathbb 1(R\ge c)$ (a scholarship for test scores above a cutoff). If units cannot precisely manipulate $R$, those just above and just below $c$ are as good as randomized, and the jump $\lim_{r\downarrow c}\mathbb E[Y\mid R=r]-\lim_{r\uparrow c}\mathbb E[Y\mid R=r]$ is the effect at $c$. Estimate each limit by a local linear regression within a bandwidth $h$ of $c$. Small $h$: little bias, few points, high variance. Large $h$: stable but the curvature of $\mathbb E[Y\mid R]$ leaks into the jump. Data-driven bandwidths (Imbens–Kalyanaraman; Calonico–Cattaneo–Titiunik) balance the two. Check for bunching of $R$ just above $c$ (McCrary's density test); bunching means manipulation. The project's simulated jump of 2 is recovered as 2.003 at $h=0.3$, while a naive comparison of all treated and untreated gives 2.802.

**Difference-in-differences.** Two groups ($G=1$ treated, $G=0$ control), two periods ($T=0$ before, $T=1$ after), treatment only for $G=1,T=1$. Estimator: $\hat\tau=(\bar Y_{11}-\bar Y_{10})-(\bar Y_{01}-\bar Y_{00})$. **Parallel trends** says $\mathbb E[Y^0_{T=1}-Y^0_{T=0}\mid G=1]=\mathbb E[Y^0_{T=1}-Y^0_{T=0}\mid G=0]$: absent treatment, both groups would have changed by the same amount. Then the treated group's change is ATT plus the common trend, the control's change is the common trend, and the difference is the ATT. The regression form is $Y=\beta_0+\beta_1G+\beta_2T+\beta_3G\cdot T+\varepsilon$, with $\hat\beta_3=\hat\tau$ exactly (four cells, four parameters). With several periods, an **event-study plot** shows the coefficients on leads and lags of treatment; leads (pre-treatment periods) near zero support parallel trends, though they cannot prove it for the post-treatment period. Card & Krueger's 1994 New Jersey/Pennsylvania minimum-wage study is the textbook case. In the project, diverging trends of 0.7 produce a DiD of 2.699 against a true 2. Staggered adoption with heterogeneous effects breaks the two-way fixed-effects regression (Goodman-Bacon 2021); use the newer estimators.

**Synthetic control** (Abadie & Gardeazabal 2003; Abadie, Diamond & Hainmueller 2010) builds the counterfactual for one treated unit (California after its 1988 tobacco tax) as a weighted average of untreated units chosen to match its pre-treatment path.

### 3.10 Sensitivity analysis

No observational analysis can rule out an unmeasured confounder. It can say how strong one would have to be.

**Linear-bias version.** Truth: $Y=\beta A+\theta^\top L+\gamma U+\varepsilon$. You regress $Y$ on $(A,L)$. By ST-3 Theorem 3.12, with $\delta$ the coefficient on $A$ in the population regression of $U$ on $(A,L)$,
1. $\hat\beta\to\beta+\gamma\delta$, so
$$\beta=\hat\beta-\gamma\delta.\tag{3.24}$$
2. The effect is explained away when $\gamma\delta=\hat\beta$: a hyperbola in the $(\gamma,\delta)$ plane. For $\hat\beta=0.30$, a confounder with $\gamma=0.5$ needs $\delta=0.6$.
3. Benchmark against measured covariates: if no observed $L_j$ comes close to that $(\gamma,\delta)$ pair, an unobserved one would have to be stronger than anything you measured (Cinelli & Hazlett 2020 formalize this with partial $R^2$).

For risk ratios, VanderWeele & Ding (2017) define the **E-value**: the minimum strength of association, on the risk-ratio scale, that $U$ must have with both $A$ and $Y$ to explain away an observed risk ratio $RR$:
$$E=RR+\sqrt{RR(RR-1)}.\tag{3.25}$$
For $RR=2$, $E=3.414$ (derivation from the Ding–VanderWeele bounding factor in their paper). **Rosenbaum bounds** ask the same question for matched designs: how large may the odds ratio $\Gamma$ of treatment between two matched units differ before the conclusion fails?

### 3.11 Causal discovery, and causal questions in the other modules

**Causal discovery.** From observational data alone you can learn at most the **Markov equivalence class** of the DAG: graphs with the same skeleton and the same **v-structures** (colliders $X\to C\leftarrow W$ with $X,W$ non-adjacent) imply the same conditional independences (Verma & Pearl 1990). A chain $A\to M\to Y$, the reverse chain and the fork $A\leftarrow M\to Y$ all say "$A\perp Y\mid M$" and nothing else, so no amount of data distinguishes them. The collider says "$A\perp Y$, but not given $C$," a different pattern, which is why arrow directions can be learned only at colliders (and propagated from there). The PC algorithm (Spirtes, Glymour & Scheines) does this by testing conditional independences; it needs faithfulness and no hidden confounders, and its output is only as good as hundreds of tests.

**ST-3's regressions** are associations: the coefficient on $x_j$ is the FWL partial slope, causal only if a stated DAG makes the other regressors a backdoor set for $x_j$ (§3.7).

**Granger causality** (Granger 1969). $X$ "Granger-causes" $Y$ if past $X$ improves the prediction of $Y_t$ beyond past $Y$, tested by an F test on added lags in ST-7's autoregressions. It is a statement about prediction. A common cause $Z$ that drives $X$ with lag 1 and $Y$ with lag 2 makes past $X$ a leading indicator of $Y$: the project finds $F=3214.6$ although there is no arrow $X\to Y$ and setting $X$ to anything leaves $Y$ unchanged. Granger causality is causal only under "no unmeasured common causes," the same assumption as everywhere else.

**Machine learning.** A predictive model estimates $\mathbb E[Y\mid X]$, which is the right object for prediction under the same distribution and the wrong object for deciding an action, because acting changes the distribution (the **shifted distribution** problem). The canonical real case: a pneumonia-risk model (Caruana et al. 2015) learned that asthma patients had *lower* mortality, because asthmatics were routinely sent to intensive care. Used to triage, it would have sent asthmatics home. The project simulates the same thing: treatment is given to the sicker, the model sees a proxy for severity, achieves $R^2=0.65$, and assigns treatment a coefficient of $+0.12$ (harmful) when the true effect is $-1.00$. The policy "treat nobody" it implies raises the mean bad outcome from $-1.01$ to $-0.01$. Prediction accuracy is not evidence of causal correctness. **Reinforcement learning** escapes this partly by experimenting: an agent that chooses actions and sees results is running a sequential randomized trial.

### 3.12 Hill's criteria and the honest summary

Hill (1965) listed nine considerations for judging whether an observed association is causal: **strength** (large associations are harder to explain by confounding; that is the E-value's logic), **consistency** across studies and populations, **specificity** (one exposure, one disease, now considered weak), **temporality** (cause before effect, the only one that is necessary), **biological gradient** (dose-response), **plausibility** (a mechanism), **coherence** with known facts, **experiment** (does removing the exposure reduce the outcome?), and **analogy**. Hill wrote that none of the nine can bring indisputable evidence and none can be required as a *sine qua non*; later writers single out temporality as the one that is logically necessary. They are a checklist for argument, not a test.

The honest summary is short. Observational causal inference is assumption-driven. Every method in §3.4–3.10 converts an untestable assumption (no unmeasured confounding, exclusion, continuity at a threshold, parallel trends) into an estimate. The assumptions must be stated, as a DAG or an equivalent list, so that others can attack them. "The data speak for themselves" is false: the same table supports opposite conclusions under different DAGs (Example 4.1), and only the DAG, which comes from outside the data, says which.

## 4. Worked examples

### Example 4.1: Simpson's paradox, both DAGs

Charig et al. (1986) compared open surgery and percutaneous nephrolithotomy (PCNL) for kidney stones (F4 Trap 4):

| | small stones | large stones | all |
|---|---|---|---|
| open | 81/87 = 0.9310 | 192/263 = 0.7300 | 273/350 = 0.7800 |
| PCNL | 234/270 = 0.8667 | 55/80 = 0.6875 | 289/350 = 0.8257 |

Open surgery wins in each stratum and loses overall, because it was given mostly to large stones (263 of 350), which are harder.

*DAG 1: stone size is a confounder,* size → treatment, size → recovery, treatment → recovery. The backdoor criterion says adjust for size. Standardize to the pooled size distribution (357 small, 343 large of 700):
- open: $0.9310\times\frac{357}{700}+0.7300\times\frac{343}{700}=0.8325$;
- PCNL: $0.8667\times\frac{357}{700}+0.6875\times\frac{343}{700}=0.7789$.

Open surgery is better by 5.4 percentage points.

*DAG 2: the stratifying variable is a mediator.* Keep the numbers and relabel the column variable as something measured *after* treatment that the treatment affects and that affects recovery: Pearl's example is blood pressure after a drug, treatment → BP → recovery. Now adjusting is conditioning on a mediator, and the crude comparison (0.7800 vs 0.8257) is the total effect: the second treatment is better by 4.6 points. Same numbers, opposite answers. Which DAG holds is a fact about the world (stone size exists before the doctor chooses), not about the table.

### Example 4.2: the selection-bias decomposition with eight units

Health scores (higher is better). The sicker four were treated.

| unit | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| $Y^0$ | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| $Y^1$ | 5 | 6 | 6 | 7 | 7 | 8 | 10 | 11 |
| $Y^1-Y^0$ | 3 | 3 | 2 | 2 | 1 | 1 | 2 | 2 |
| $A$ | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 |

- ATE $=16/8=2.0$. ATT $=(3+3+2+2)/4=2.5$.
- Observed: $\mathbb E[Y\mid A=1]=(5+6+6+7)/4=6.0$; $\mathbb E[Y\mid A=0]=(6+7+8+9)/4=7.5$. Association $=-1.5$: the treatment looks harmful.
- Selection bias $=\mathbb E[Y^0\mid A=1]-\mathbb E[Y^0\mid A=0]=3.5-7.5=-4.0$.
- Check (3.2): $2.5+(-4.0)=-1.5$. ✓

Now randomize instead: all $\binom84=70$ assignments of four treated. The project enumerates them: the mean of $\hat\tau$ is exactly 2.0, and its variance is 2.5000, matching (3.6) with $S_1^2=4.2857$, $S_0^2=6$, $S_\tau^2=0.5714$: $4.2857/4+6/4-0.5714/8=2.5000$. The conservative bound is 2.5714. Blocking on the sick/healthy split (units 1–4 and 5–8, two treated per block, 36 assignments) gives variance 0.8750 by enumeration and by (3.7), a reduction by a factor of 2.9, because the blocks differ in $Y^0$ by 4 points.

### Example 4.3: Fisher's exact test with six units

Three of six units are randomized to treatment. Outcomes: treated 7, 9, 5; control 4, 6, 3. Observed difference $7-4.333=2.667$.

Under the sharp null (3.4) the six outcomes are fixed; the difference for an assignment with treated sum $s$ is $(s-(34-s))/3=(2s-34)/3$. The observed $s=21$. Enumerate the $\binom63=20$ triples from $\{9,7,6,5,4,3\}$ with $s\ge21$: $\{9,7,6\}$ (22) and $\{9,7,5\}$ (21). So the one-sided $p=2/20=0.10$. For the two-sided test we also need $s\le13$: $\{3,4,5\}$ (12) and $\{3,4,6\}$ (13), so $p=4/20=0.20$. The smallest possible one-sided p-value with this design is $1/20=0.05$: a six-unit experiment can never be significant at a level below 0.05, which is worth knowing before running it.

### Example 4.4: standardization and IPW on a $2\times2\times2$ table

1000 patients, binary severity $L$, treatment $A$, recovery $Y$:

| $L$ | $A$ | $n$ | recovered | $\hat P(Y=1\mid A,L)$ |
|---|---|---|---|---|
| 0 | 0 | 450 | 90 | 0.20 |
| 0 | 1 | 150 | 45 | 0.30 |
| 1 | 0 | 100 | 50 | 0.50 |
| 1 | 1 | 300 | 210 | 0.70 |

$P(L=1)=0.4$, $P(A=1\mid L=0)=0.25$, $P(A=1\mid L=1)=0.75$.

*Crude:* $\mathbb E[Y\mid A=1]=255/450=0.5667$, $\mathbb E[Y\mid A=0]=140/550=0.2545$, difference 0.3121.

*Standardization (3.10):* $\mathbb E[Y^1]=0.30(0.6)+0.70(0.4)=0.46$; $\mathbb E[Y^0]=0.20(0.6)+0.50(0.4)=0.32$; ATE $=0.14$.

*IPW (3.11):* weights are $1/0.25=4$ for $(L{=}0,A{=}1)$, $1/0.75=4/3$ for $(0,0)$, $4/3$ for $(1,1)$, $4$ for $(1,0)$.
$\mathbb E[Y^1]=(45\cdot4+210\cdot\tfrac43)/1000=(180+280)/1000=0.46$; $\mathbb E[Y^0]=(90\cdot\tfrac43+50\cdot4)/1000=(120+200)/1000=0.32$.

Identical, as §3.4's saturated-case proof requires. In the pseudo-population, the 150 treated low-severity patients become 600, the full size of that stratum: everyone is represented under each treatment. More than half the crude difference was confounding.

### Example 4.5: a collider, by Wright's rules

Linear DAG $A\to Y$ (0.5), $A\to C$ (1), $Y\to C$ (1), unit noise variances. The true effect of $A$ on $Y$ is 0.5. An analyst "controls for" $C$.

*Covariances by treks (3.19).* $\operatorname{Var}A=1$. $\operatorname{Cov}(A,Y)$: one trek (top $A$, $A\to Y$), $=0.5$. $\operatorname{Var}Y$: top $A$ ($0.5^2$) plus top $Y$ (1) $=1.25$. $\operatorname{Cov}(A,C)$: top $A$ with directed paths $A$ (trivial) and $A\to C$ (1) or $A\to Y\to C$ (0.5): $1+0.5=1.5$. $\operatorname{Cov}(Y,C)$: top $A$: $(0.5)(1)+(0.5)(0.5)=0.75$; top $Y$: $1\cdot1=1$; total 1.75. $\operatorname{Var}C$: top $A$: $(1.5)^2=2.25$; top $Y$: 1; top $C$: 1; total 4.25. The project checks all six against (3.18) and a sample of 200,000.

*Regression of $Y$ on $(A,C)$.* $\begin{pmatrix}1&1.5\\1.5&4.25\end{pmatrix}^{-1}\begin{pmatrix}0.5\\1.75\end{pmatrix}=\frac12\begin{pmatrix}4.25&-1.5\\-1.5&1\end{pmatrix}\begin{pmatrix}0.5\\1.75\end{pmatrix}=\begin{pmatrix}-0.25\\0.5\end{pmatrix}$.

The coefficient on $A$ is $-0.25$. The true effect is $+0.5$. Adjusting for a common effect of treatment and outcome reversed the sign. Without the $A\to Y$ arrow (pure collider) the same calculation gives $-0.5$ from a true 0.

### Example 4.6: the backdoor criterion on a five-node DAG

DAG: $Z_1\to Z_3\leftarrow Z_2$, $Z_1\to A$, $Z_3\to A$, $Z_3\to Y$, $Z_2\to Y$, $A\to Y$ (Pearl's classic, without the mediator). Coefficients 0.8, 0.8, 0.7, 0.6, 0.5, 0.7 and effect 1.0.

Backdoor paths from $A$:
1. $A\leftarrow Z_3\to Y$ (open unless $Z_3\in L$);
2. $A\leftarrow Z_3\leftarrow Z_2\to Y$ (open unless $Z_3$ or $Z_2\in L$);
3. $A\leftarrow Z_1\to Z_3\to Y$ (open unless $Z_1$ or $Z_3\in L$);
4. $A\leftarrow Z_1\to Z_3\leftarrow Z_2\to Y$ ($Z_3$ is a collider: blocked unless $Z_3\in L$; then open unless $Z_1$ or $Z_2\in L$).

Path 1 forces $Z_3\in L$; then path 4 needs $Z_1$ or $Z_2$. Valid: $\{Z_1,Z_3\}$, $\{Z_2,Z_3\}$, $\{Z_1,Z_2,Z_3\}$. Invalid: $\varnothing,\{Z_1\},\{Z_2\},\{Z_1,Z_2\}$ (path 1 open) and $\{Z_3\}$ (path 4 opened). The project's d-separation code lists the same four paths and the same three valid sets, and the population regression coefficient on $A$ is exactly 1.0000 for each valid set, while the invalid ones give 1.4358, 1.5206, 1.2805, 1.2206 and, for $\{Z_3\}$, 0.8983. Adjusting for the one variable that "obviously" confounds ($Z_3$ causes both) makes the answer wrong in a new direction.

### Example 4.7: the Wald estimator from a $2\times2$ table of means

An encouragement design: $Z=1$ means a randomly chosen half of employees got a letter encouraging enrolment in a training programme ($A$); $Y$ is a later productivity score.

| | $\mathbb E[A\mid Z]$ | $\mathbb E[Y\mid Z]$ |
|---|---|---|
| $Z=1$ | 0.70 | 5.00 |
| $Z=0$ | 0.20 | 4.10 |

ITT effect on $Y$: $0.90$. First stage: $0.50$. Wald (3.20): $0.90/0.50=1.80$. Under monotonicity (nobody refuses *because* of the letter) the compliers are 50% of employees and the training raises their score by 1.8 (LATE, 3.21). The comparison of enrolled with non-enrolled would mix in whatever made people enrol.

### Example 4.8: a $2^2$ factorial

A process with temperature $A$ and pressure $B$ at two levels each, cell means of yield:

| | $B$ low | $B$ high |
|---|---|---|
| $A$ low | 60 | 54 |
| $A$ high | 72 | 80 |

By (3.9), in standard order $(--,+-,-+,++)=(60,72,54,80)$:
- $A$: $\frac12(-60+72-54+80)=19$;
- $B$: $\frac12(-60-72+54+80)=1$;
- $AB$: $\frac12(60-72-54+80)=7$.

Pressure has almost no main effect, but the effect of temperature is $12$ at low pressure and $26$ at high pressure; the interaction $7$ is half that difference. An OFAT experimenter who varied pressure at low temperature would see $54-60=-6$, conclude pressure hurts, and never find the best corner (80). The interaction plot (project, `interaction.png`) shows non-parallel lines.

## 5. Why it matters and where it connects

- **ST-3**: every coefficient table. The OVB formula (ST-3 Theorem 3.12) is (3.24); FWL (ST-3 Theorem 3.4) is what "adjusting" does, and §3.7's four conditions are when it means anything.
- **ST-2**: the permutation test (Theorem 3.6) is Fisher's randomization test; the power calculation (§3.9) sizes A/B tests; optional stopping (exercise C3) is the peeking problem.
- **ST-4**: propensity scores are logistic regressions fitted by IRLS; in the project the propensity is fitted by the same Newton step.
- **ST-5/ST-6**: Bayesian models are as causal as their DAG; a hierarchical model with the wrong adjustment set is wrong with a credible interval.
- **ST-7**: Granger tests are predictive (§3.11).
- **LA-5**: 2SLS is a projection followed by a regression (3.23); factorial contrasts are orthogonal columns.
- **F5**: the conditional-covariance formula (3.14) is the Gaussian conditioning formula; d-separation reads it off a graph.
- **04-machine-learning**: ML-1's generalization guarantees hold under the training distribution, and an intervention changes the distribution. Feature importances in ML-5 and coefficients in ML-3 are associations. Distribution shift, fairness and off-policy evaluation all need §3.4–3.6.
- **Engineering**: §3.3 is design of experiments for process tuning and reliability testing; blocking by machine, batch or day is (3.7).

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 drill the definitions with small tables: the decomposition, adjustment sets, Fisher's test, the Wald estimator, a factorial, a DiD. B1–B7 are the module: B1 (standardization = IPW), B3 (Wright's rules and a bad control), B5 (adjustment sets and their precision, in `project.py`) and B6 (2SLS) are essential. C1–C3: LATE with defiers, the propensity score theorem, and weak instruments with many instruments (code).

## 7. Coding project

`project.py` is a small causal simulation library in 13 parts: (1) a linear-Gaussian `SCM` class with `sample(do=...)`, the covariance (3.18), trek enumeration (3.19), d-separation and the backdoor check, verified on Example 4.5; (2) Example 4.1; (3) Example 4.2 with Neyman's variance and blocking by exact enumeration; (4) Fisher's test by full enumeration (Example 4.3) and on a 10-unit experiment against the $t$ test; (5) Example 4.4, then a simulation where the association is far from the ATE under confounding and equal to it under randomization, and the g-formula and IPW (logistic propensity by Newton) recover ATE $=2$; (6) Berkson's collider, the mediator traps (3.16)–(3.17), M-bias, and Example 4.6 over all eight adjustment sets; (7) factorial designs with interaction plots; (8) Example 4.7, OLS vs 2SLS, and 2SLS bias against first-stage $F$; (9) RD with bandwidths and DiD with and without parallel trends; (10) the sensitivity plot; (11) Granger from a common cause and the prediction-vs-intervention model; (12) A/B sizing and peeking; (13) exercise numbers.

**Build it yourself:** before reading the code, write `blocked` and `backdoor_ok` from Definition 3.6 and Theorem 3.7 (about 15 lines), and `neyman_exact` (enumerate assignments, about 8 lines). "Your turn" asks for propensity-score matching, a DAG of your own with three adjustment sets, Monty Hall as a collider, and a positivity stress test.

## 8. Where you'll struggle

- **Potential-outcomes notation looks like overkill.** $Y^1$, $Y^0$, $Y^{A(1)}$: it feels like decoration until (3.2), when the selection bias appears as a number you could not have named without it. Do Example 4.2 by hand with your own eight units. After that, write every causal claim you read as an equation in $Y^a$.
- **"Adjust for everything you have" is the most common and most harmful mistake.** It feels conservative. It is not: colliders (Example 4.5) and mediators (3.16) make adjusted estimates worse than crude ones, and you cannot tell from the data. Run the project's Part 6 and change the coefficients until the sign flip stops surprising you. Then draw the DAG *before* choosing covariates, every time.
- **DAGs encode assumptions you cannot test.** People find this unsatisfying and look for a data-driven way to pick adjustment sets. There isn't one (§3.11: equivalence classes). This is the point of the DAG, not its weakness: it puts the untestable part on paper where it can be argued about. A missing arrow is a strong claim; check each one.
- **Collider bias hides inside study design.** Selection into a sample (hospital patients, survivors, respondents, published results) is conditioning. Ask "how did a unit get into this dataset?" before anything else.
- **IV's exclusion restriction is almost never testable and almost always arguable.** The first stage you can check; independence sometimes comes from a lottery; exclusion is a story. Weak instruments add bias toward OLS and bad intervals: report the first-stage $F$ and use Anderson–Rubin intervals if it is below 10.
- **Reading regression output differently.** After this module a coefficient on $x_j$ is causal only when a stated DAG makes the other regressors a backdoor set for $x_j$ and none of them is downstream of it. That is rarely true for more than one coefficient in a table (Table 2 fallacy). This should change how you read ST-3's output permanently.
- **Scope.** This is 8 hours on a subject with semester-long courses and an active research frontier (heterogeneous effects, double machine learning, staggered DiD, synthetic control inference). It is here because without it every other module in the course is dangerous. Treat §10's Hernán & Robins as the continuation.

## 9. Check yourself

1. Write the selection-bias decomposition (3.2) and say which term randomization kills and why.
2. State Fisher's sharp null and explain why it makes the randomization distribution computable. How does it differ from Neyman's null?
3. In (3.6), which term cannot be estimated, and in which direction does dropping it err?
4. Derive the g-formula from consistency, conditional exchangeability and positivity.
5. Why do IPW and standardization give the same number for a saturated discrete $L$?
6. What does conditioning do in a fork, a chain and a collider? Give the linear-Gaussian formula for the collider.
7. State the backdoor criterion. Give a pre-treatment variable it forbids adjusting for.
8. Adjusting for a mediator gives which effect? What goes wrong if the mediator has its own confounder with $Y$?
9. State the trek rule and explain why covariance never flows through a collider.
10. Derive the Wald estimator and say what it estimates with heterogeneous effects under monotonicity.
11. Why does 2SLS drift toward OLS when the instrument is weak, and what is the rule of thumb?
12. What is identifiable about a DAG from observational data alone, and why are colliders special?

## 10. Pair with

- **Hernán & Robins, *Causal Inference: What If* (free PDF), ch. 1–7.** The potential-outcomes and DAG views merged, with exactly these conventions ($Y^a$, exchangeability, positivity, consistency); ch. 2–3 on randomization and observational studies, ch. 6–7 on graphs and confounding. The best continuation; Part II (ch. 11–18) does the modelling these notes only name.
- **McElreath, *Statistical Rethinking*, 2nd ed., ch. 5–6.** The clearest DAG introduction anywhere: the fork, pipe, collider and descendant ("the haunted DAG") with simulations; read ch. 6 for collider and post-treatment bias with more examples than §3.6.
- **Pearl, Glymour & Jewell, *Causal Inference in Statistics: A Primer*, ch. 1–3.** Short and formal: d-separation, the do-operator, backdoor and front-door with proofs these notes give only in the simplest case.
- **Angrist & Pischke, *Mostly Harmless Econometrics*, ch. 2–4.** The econometric view: selection bias, regression as a weighted average of effects, bad controls, IV, LATE and 2SLS with the standard errors done properly.
- **Cunningham, *Causal Inference: The Mixtape* (free online).** Worked examples with data and code for matching, RD, DiD and synthetic control, the designs §3.9 treats in a paragraph each.
- **Fisher, *The Design of Experiments*, ch. 2.** The lady tasting tea: the original randomization argument in Fisher's own words, still the best statement of why the physical act of randomizing is what licenses the test.
