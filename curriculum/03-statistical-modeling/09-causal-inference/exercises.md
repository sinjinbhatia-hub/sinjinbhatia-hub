# ST-9 Exercises

Conventions as in the notes: $A$ treatment, $Y$ outcome, $L$ measured covariates, $U$ unmeasured, $Z$ instrument; potential outcomes $Y^1,Y^0$; ATE $=\mathbb E[Y^1-Y^0]$, ATT $=\mathbb E[Y^1-Y^0\mid A=1]$. Linear DAGs have unit noise variances unless stated. Equation numbers refer to `notes.md`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Six units:

| unit | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| $Y^0$ | 4 | 6 | 5 | 8 | 7 | 9 |
| $Y^1$ | 6 | 7 | 8 | 8 | 9 | 10 |
| $A$ | 1 | 1 | 1 | 0 | 0 | 0 |

Compute the ATE, the ATT, the observed association and the selection bias, and check (3.2).

**A2.** DAG: $L\to A$, $L\to Y$, $A\to M\to Y$, $A\to Y$, $A\to C\leftarrow Y$. List every path from $A$ to $Y$, mark which are backdoor paths, and find every subset of $\{L,M,C\}$ that satisfies the backdoor criterion. For each subset that fails, say why.

**A3.** Four units, two randomized to treatment. Treated outcomes 5, 7; control outcomes 2, 3. Compute the exact one-sided Fisher p-value for the difference in means. What is the smallest p-value this design can ever produce, and what does that tell you about planning?

**A4.** An encouragement experiment gives $\mathbb E[A\mid Z=1]=0.45$, $\mathbb E[A\mid Z=0]=0.15$, $\mathbb E[Y\mid Z=1]=12.6$, $\mathbb E[Y\mid Z=0]=12.0$. Compute the ITT effect, the first stage, and the Wald estimate. Under monotonicity, what fraction of the population are compliers, and for whom is the Wald estimate an effect?

**A5.** A $2^2$ factorial has cell means $(A,B)=(-,-):20$, $(+,-):30$, $(-,+):24$, $(+,+):26$. Compute both main effects and the interaction with (3.9). What is the effect of $A$ at each level of $B$? Would you report "B has no effect"?

**A6.** Average outcomes: treated group 10 before, 15 after; control group 8 before, 11 after. Compute the difference-in-differences estimate. State the assumption in words and give one concrete way it could fail.

## B. Core (do all)

**B1.** Five hundred patients:

| $L$ | $A$ | $n$ | $Y=1$ |
|---|---|---|---|
| 0 | 0 | 160 | 48 |
| 0 | 1 | 40 | 20 |
| 1 | 0 | 60 | 30 |
| 1 | 1 | 240 | 168 |

(a) Compute the crude risk difference. (b) Compute $\mathbb E[Y^1]$, $\mathbb E[Y^0]$ and the ATE by standardization (3.10). (c) Compute them by IPW (3.11) and confirm they agree. (d) Build the IPW pseudo-population table (four cells, weighted counts) and show $A\perp L$ in it.

**B2.** Use A1's population with complete randomization of 3 of 6 units.
(a) Compute $S_1^2$, $S_0^2$, $S_\tau^2$ and Neyman's variance (3.6).
(b) Verify by enumerating all 20 assignments (by hand or with `neyman_exact` in `project.py`).
(c) Compute the expected value of the conservative estimator $S_1^2/3+S_0^2/3$ and say by how much it overstates the variance.
(d) Prove that the conservative estimator is unbiased for the true variance if and only if the treatment effect is the same for every unit.

**B3.** Linear DAG: $L\to A$ (coefficient 1), $L\to Y$ (1), $A\to Y$ (0.5).
(a) Use the trek rule (3.19) to compute $\operatorname{Var}A$ and $\operatorname{Cov}(A,Y)$.
(b) Compute the population coefficient of $Y$ on $A$ alone, and on $(A,L)$. Identify the omitted-variable bias term with ST-3 Theorem 3.12.
(c) Add a collider: $A\to C$ (1), $Y\to C$ (1). Compute the coefficient on $A$ in the regression of $Y$ on $(A,L,C)$. Check with `SCM.pop_coef`.

**B4.** (Selection as collider.) In a population, $A$ and $Y$ are independent binary variables with $P(Y=1)=0.2$ and $P(A=1)=0.3$. A unit enters the study with probability $P(S=1\mid A,Y)=0.1+0.4A+0.4Y$.
(a) Draw the DAG.
(b) Compute $P(Y=1\mid A=1,S=1)$ and $P(Y=1\mid A=0,S=1)$ and their ratio.
(c) Explain the result in one sentence to a clinician who thinks "the study sample is just a subset."

**B5.** (Code.) Build this DAG with `SCM` in `project.py`: $L_1\to A$ (1), $L_1\to L_2$ (1), $L_2\to Y$ (1), $A\to M$ (1), $M\to Y$ (0.8), $A\to Y$ (0.5), $A\to C$ (1), $Y\to C$ (1).
(a) What is the total effect of $A$ on $Y$?
(b) By hand, find all valid backdoor adjustment sets among subsets of $\{L_1,L_2\}$. Confirm with `backdoor_ok`.
(c) For $\varnothing$, $\{L_1\}$, $\{L_2\}$, $\{L_1,L_2\}$, $\{M\}$, $\{L_1,C\}$, compute the population coefficient on $A$ and $n\operatorname{Var}(\hat\beta_A)=\sigma^2_{Y\mid A,S}/\operatorname{Var}(A\mid S)$.
(d) Which valid set gives the most precise estimate, and why? Give a rule.

**B6.** $Z,U\sim\mathcal N(0,1)$ independent, $A=0.5Z+U+\varepsilon_A$, $Y=A+2U+\varepsilon_Y$, unit-variance noises.
(a) Compute the probability limit of the OLS slope of $Y$ on $A$.
(b) Show that $\operatorname{Cov}(Z,Y)/\operatorname{Cov}(Z,A)=1$.
(c) Prove that 2SLS with one instrument and one regressor (plus intercepts) equals $\widehat{\operatorname{Cov}}(Z,Y)/\widehat{\operatorname{Cov}}(Z,A)$.
(d) For $n=400$, approximate the expected first-stage $F$. Is the instrument weak?
(e) Explain why the residual variance from the stage-2 OLS output gives a wrong standard error.

**B7.** An observational study reports an adjusted coefficient $\hat\beta=0.30$ (SE 0.05) for the effect of an exercise programme on a fitness score.
(a) With (3.24), find the imbalance $\delta$ needed to explain the effect away if an unmeasured confounder has $\gamma=0.5$. Sketch the explain-away curve in the $(\gamma,\delta)$ plane.
(b) A different study reports a risk ratio of 1.5. Compute its E-value and interpret it.
(c) Why is "we adjusted for 40 covariates" not a substitute for (a)?

## C. Challenge (do at least one)

**C1.** (LATE with defiers.) (a) Repeat the derivation of (3.21) without assuming monotonicity and show that the Wald ratio is $\dfrac{\pi_c\,\mathbb E[\tau\mid c]-\pi_d\,\mathbb E[\tau\mid d]}{\pi_c-\pi_d}$. (b) Take $\pi_c=0.4$, $\mathbb E[\tau\mid c]=2$, $\pi_d=0.1$, $\mathbb E[\tau\mid d]=5$. Compute the Wald ratio and comment. (c) Give a real setting in which defiers are plausible.

**C2.** (Front door.) DAG: $U\to A$, $U\to Y$, $A\to M\to Y$, no other arrows, $U$ unobserved. Starting from the truncated factorization (3.13), derive
$$P(y\mid do(a))=\sum_mP(m\mid a)\sum_{a'}P(y\mid m,a')P(a').$$
Say at each step which conditional independence of the DAG you use.

**C3.** (Code: weak instruments.) Extend Part 8 of `project.py` so that $A$ depends on $k$ instruments, each with coefficient $\pi$. With $n=200$ and 1000 replications, compare $k=1,\pi=0.15$ with $k=5,\pi=0.15$: report the mean first-stage $F$, the median and mean bias of 2SLS, and the OLS bias. Compare the ratio of the 2SLS bias to the OLS bias with $1/F$. Why is the mean bias for $k=1$ not a meaningful number?
