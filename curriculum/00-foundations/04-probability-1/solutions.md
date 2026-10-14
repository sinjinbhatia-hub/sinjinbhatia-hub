# F4 Solutions

Conventions as in the notes and exercises. Every number below is asserted or printed in `project.py`, mostly in Part 10; the simulation figures quoted come from that seeded run.

## A. Warm-up

**A1.** (a) Complement: "at least one 6" fails only if all four rolls avoid 6. The rolls are independent, so $P(\text{no } 6) = (5/6)^4 = 625/1296$, and
$$P(\text{at least one } 6) = 1 - \tfrac{625}{1296} = \tfrac{671}{1296} = 0.5177.$$

(b) One roll of two dice gives double-six with probability $1/36$. For 24 independent rolls, $P(\text{at least one}) = 1 - (35/36)^{24} = 0.4914$.

(c) De Méré used the linear rule "chance $\approx$ number of tries $\times$ chance per try", which says both bets are worth $4 \cdot \tfrac16 = 24 \cdot \tfrac{1}{36} = \tfrac23$. The true value $1 - (1-p)^n$ depends on $n$ and $p$ separately, not only through $np$: with $np = 2/3$ held fixed it decreases as $p$ shrinks, toward the Poisson limit $1 - e^{-2/3} = 0.4866$ (§3.10). The linear rule is only the first-order term.

(d) We need $(35/36)^n < 1/2$, i.e. $n > \ln 2 / \ln(36/35) = 24.6$. So $n = 25$, giving $1 - (35/36)^{25} = 0.5055$.

**A2.** (a) Count the ways to choose a committee of $k$ people from $n$ and name one member as chair. Choosing the committee first ($\binom nk$ ways) and then the chair from its members ($k$ ways) gives $k\binom nk$. Choosing the chair first ($n$ ways) and then the other $k - 1$ members from the remaining $n - 1$ people gives $n\binom{n-1}{k-1}$. Both count the same set of (committee, chair) pairs. ∎

(b) A room holds $n$ men and $n$ women. The right side counts the $n$-person committees from all $2n$ people. Sort them by the number $k$ of men: there are $\binom nk$ ways to choose the men and $\binom{n}{n-k}$ ways to choose the $n - k$ women, and $\binom{n}{n-k} = \binom nk$ by Theorem 3.2(a). Summing over $k$ counts each committee once, so $\sum_k \binom nk^2 = \binom{2n}{n}$. ∎ (This is Vandermonde with $m = n$, $r = n$.) For $n = 5$: $1 + 25 + 100 + 100 + 25 + 1 = 252 = \binom{10}{5}$.

**A3.** From the table: $P(F) = 80/200 = 0.4$, $P(S) = 150/200 = 0.75$, $P(F \cap S) = 72/200 = 0.36$.
- $P(S \mid F) = P(F \cap S)/P(F) = 72/80 = 0.9$.
- $P(S \mid F^c) = 78/120 = 0.65$.
- $P(F \mid S) = 72/150 = 0.48$.
- $P(F \mid S^c) = 8/50 = 0.16$.

Bayes (3.10): $P(F \mid S) = P(S \mid F)P(F)/P(S) = 0.9 \times 0.4 / 0.75 = 0.48$. ✓

Independence test: $P(F \cap S) = 0.36$ but $P(F)P(S) = 0.4 \times 0.75 = 0.30$. Not equal, so $F$ and $S$ are dependent (equivalently $P(S \mid F) = 0.9 \neq 0.75 = P(S)$). These are observational counts, so the dependence says nothing yet about F4 *causing* passes (ST-9).

**A4.** (a) $\Omega$ is the 36 equally likely pairs $(i, j)$. The maximum is at most $k$ exactly when both dice are at most $k$, which is $k^2$ pairs. So $F_X(k) = k^2/36$ for $k = 1, \dots, 6$. The PMF is the jump: $P(X = k) = F_X(k) - F_X(k-1) = \frac{k^2 - (k-1)^2}{36} = \frac{2k - 1}{36}$, i.e. $\frac{1}{36}, \frac{3}{36}, \frac{5}{36}, \frac{7}{36}, \frac{9}{36}, \frac{11}{36}$. They sum to $36/36$. ✓

(b) $\mathbb{E}[X] = \frac{1}{36}\sum_k k(2k-1) = \frac{1 + 6 + 15 + 28 + 45 + 66}{36} = \frac{161}{36} = 4.4722$.
By LOTUS, $\mathbb{E}[X^2] = \frac{1}{36}\sum_k k^2(2k-1) = \frac{1 + 12 + 45 + 112 + 225 + 396}{36} = \frac{791}{36} = 21.9722$.
By (3.17), $\operatorname{Var}(X) = \frac{791}{36} - \left(\frac{161}{36}\right)^2 = \frac{791 \cdot 36 - 25921}{1296} = \frac{2555}{1296} = 1.9715$.

**A5.** (a) Standardize with $Z = (X - 100)/15$.
- $P(X > 130) = P(Z > 2) = 1 - \Phi(2) = 0.02275$.
- $P(85 < X < 115) = P(-1 < Z < 1) = 0.6827$ (the 68% of §3.12).
- 99th percentile: $x = 100 + 15\,\Phi^{-1}(0.99) = 100 + 15(2.3263) = 134.89$.

(b) Mean $1/\lambda = 10$ min, so $\lambda = 0.1$ per minute.
- $P(T > 20) = e^{-2} = 0.1353$.
- $P(T < 5) = 1 - e^{-0.5} = 0.3935$.
- $P(T > 30 \mid T > 10) = P(T > 20) = 0.1353$ by memorylessness (Theorem 3.16).
- Median: $1 - e^{-0.1 m} = 1/2$ gives $m = 10\ln 2 = 6.931$ min.
- $\sigma_T = 1/\lambda = 10$ min.

## B. Core

**B1.** (a) Fix $\omega$. If $\omega$ is in none of $A, B, C$, every indicator is 0, the product is $1 \cdot 1 \cdot 1 = 1$, and the right side is $0 = I_{A\cup B\cup C}(\omega)$. If $\omega$ is in at least one of them, say $A$, then the factor $1 - I_A(\omega) = 0$ kills the product, and the right side is $1 = I_{A\cup B\cup C}(\omega)$. The two functions agree at every $\omega$.

(b) Write $a, b, c$ for $I_A, I_B, I_C$. Expanding,
$$1 - (1-a)(1-b)(1-c) = a + b + c - ab - ac - bc + abc.$$
Products of indicators are indicators of intersections: $ab = I_{A \cap B}$, $abc = I_{A\cap B\cap C}$. Take expectations, use linearity (Theorem 3.8) and the bridge $\mathbb{E}[I_E] = P(E)$ (Theorem 3.9):
$$P(A \cup B \cup C) = P(A) + P(B) + P(C) - P(A\cap B) - P(A\cap C) - P(B\cap C) + P(A\cap B\cap C),$$
which is (3.6). ∎

(c) For $n$ events, $I_{\cup A_i} = 1 - \prod_{i=1}^n (1 - I_{A_i})$ by the same two-case check. Expanding the product gives one term for each nonempty subset $S \subseteq \{1, \dots, n\}$, equal to $(-1)^{|S|+1}\prod_{i \in S} I_{A_i} = (-1)^{|S|+1} I_{\cap_{i\in S} A_i}$; taking expectations term by term gives (3.7).

(d) $A_i$ = roll $i$ is a 6. Single terms: $3 \cdot \frac16$. Pairs: $3 \cdot \frac1{36}$. Triple: $\frac{1}{216}$. So $P = \frac{108 - 18 + 1}{216} = \frac{91}{216} = 0.4213$, and $1 - (5/6)^3 = 1 - \frac{125}{216} = \frac{91}{216}$. ✓

**B2.** $P(F) = 1/2$; given $F$ each toss is heads with probability $1/2$, given $F^c$ with probability $3/4$, and tosses of a given coin are independent.

(a) LOTP: $P(H_1) = \frac12\cdot\frac12 + \frac12\cdot\frac34 = \frac58$. Bayes:
$$P(F \mid H_1) = \frac{\frac12\cdot\frac12}{5/8} = \frac{1/4}{5/8} = \frac25, \qquad P(F \mid T_1) = \frac{\frac12\cdot\frac12}{3/8} = \frac23.$$

(b) LOTP again: $P(H_1 \cap H_2) = \frac12\cdot\frac14 + \frac12\cdot\frac{9}{16} = \frac{13}{32}$. So $P(H_2 \mid H_1) = \frac{13/32}{5/8} = \frac{13}{20} = 0.65$, while $P(H_2) = \frac58 = 0.625$ by the same LOTP as $H_1$.
- *Not independent:* $P(H_1 \cap H_2) = \frac{13}{32} \neq \frac{25}{64} = P(H_1)P(H_2)$.
- *Conditionally independent given $F$:* $P(H_1 \cap H_2 \mid F) = \frac14 = \frac12\cdot\frac12 = P(H_1\mid F)P(H_2\mid F)$, and given $F^c$, $\frac{9}{16} = \frac34\cdot\frac34$. This is the model's assumption that a fixed coin's tosses are independent.

The first head is evidence for the biased coin, which raises the chance of the second.

(c) *Directly:* $P(F \mid H_1 \cap H_2) = \frac{P(F)P(H_1 \cap H_2 \mid F)}{P(H_1 \cap H_2)} = \frac{\frac12\cdot\frac14}{13/32} = \frac{1/8}{13/32} = \frac{4}{13}$.

*Sequentially by (3.12) with $E = H_1$:* the prior is now $P(F \mid H_1) = \frac25$. By conditional independence $P(H_2 \mid F, H_1) = \frac12$ and $P(H_2 \mid F^c, H_1) = \frac34$. So
$$P(F \mid H_2, H_1) = \frac{\frac12\cdot\frac25}{\frac12\cdot\frac25 + \frac34\cdot\frac35} = \frac{1/5}{1/5 + 9/20} = \frac{4/20}{13/20} = \frac{4}{13}. \;✓$$

(d) $P(F \mid H_1 \cap \cdots \cap H_n) = \dfrac{\frac12(\frac12)^n}{\frac12(\frac12)^n + \frac12(\frac34)^n} = \dfrac{1}{1 + (3/2)^n}$. At $n = 5$: $1/(1 + 7.59375) = 0.1164$. As $n \to \infty$ it tends to 0. That is right: a fair coin produces $n$ heads in a row with probability $2^{-n}$, the biased one with probability $(3/4)^n$, and their ratio $(2/3)^n \to 0$; evidence accumulates against the hypothesis that explains the data worse.

**B3.** (a) $\Omega = \{HH, HT, TH, TT\}$, each $1/4$. $A = \{HH, HT\}$, $B = \{HH, TH\}$, $C = \{HH, TT\}$, each with probability $1/2$. $A\cap B = A \cap C = B \cap C = \{HH\}$, probability $\frac14 = \frac12\cdot\frac12$: all three pairs are independent. $A \cap B \cap C = \{HH\}$ has probability $\frac14 \neq \frac18$.

(b) $P(A) = P(B) = \frac12$, $P(C) = \frac46 = \frac23$. $A \cap B \cap C = \{3\}$, so $P(A\cap B\cap C) = \frac16 = \frac12\cdot\frac12\cdot\frac23$. But $A \cap B = \{3\}$ has probability $\frac16 \neq \frac14 = P(A)P(B)$. So the triple product condition holds without pairwise independence: mutual independence needs all four conditions.

(c) *Intersection.* Using the triple condition and then the $(B, C)$ pair condition:
$$P(A \cap (B \cap C)) = P(A)P(B)P(C) = P(A)\,P(B \cap C).$$
*Union.* Distribute, then use (3.5):
$$P(A \cap (B \cup C)) = P((A\cap B) \cup (A \cap C)) = P(A\cap B) + P(A\cap C) - P(A\cap B\cap C).$$
By mutual independence this is $P(A)P(B) + P(A)P(C) - P(A)P(B)P(C) = P(A)\big[P(B) + P(C) - P(B)P(C)\big]$. The bracket is $P(B) + P(C) - P(B\cap C) = P(B \cup C)$. So $P(A \cap (B\cup C)) = P(A)P(B \cup C)$. ∎

(d) In (a), $B \cap C = \{HH\}$ has probability $\frac14$, and $A \cap (B \cap C) = \{HH\}$ has probability $\frac14 \neq P(A)P(B\cap C) = \frac18$. Knowing $B \cap C$ (both tosses heads) determines $A$.

**B4.** (a) For $|q| < 1$, $\sum_{k\ge0} q^k = (1-q)^{-1}$. Differentiate twice in $q$: $\sum_{k \ge 2} k(k-1)q^{k-2} = 2(1-q)^{-3}$. Then
$$\mathbb{E}[X(X-1)] = \sum_{k\ge 2} k(k-1) q^k p = pq^2 \sum_{k\ge2}k(k-1)q^{k-2} = \frac{2pq^2}{p^3} = \frac{2q^2}{p^2}.$$
So $\mathbb{E}[X^2] = \frac{2q^2}{p^2} + \frac qp$ and
$$\operatorname{Var}(X) = \frac{2q^2}{p^2} + \frac qp - \frac{q^2}{p^2} = \frac{q^2 + pq}{p^2} = \frac{q(q+p)}{p^2} = \frac{q}{p^2}.$$

(b) $\mathbb{E}[Y] = \frac qp + 1 = \frac{q + p}{p} = \frac1p$. $\operatorname{Var}(Y) = \operatorname{Var}(X + 1) = \frac{q}{p^2}$ by Theorem 3.12: shifting every value by 1 moves the centroid, not the spread.

(c) $X \ge k$ exactly when the first $k$ trials all fail, so $P(X \ge k) = q^k$. (Or sum: $\sum_{j\ge k} q^j p = pq^k/(1-q) = q^k$.) Then for nonnegative integers $s, t$, since $\{X \ge s+t\} \subseteq \{X \ge s\}$,
$$P(X \ge s+t \mid X \ge s) = \frac{q^{s+t}}{q^s} = q^t = P(X \ge t).$$

(d) $p = 0.2$, $q = 0.8$: $\mathbb{E}[X] = 0.8/0.2 = 4$; $\operatorname{Var}(X) = 0.8/0.04 = 20$; $\sigma_X = \sqrt{20} = 4.4721$; $\mathbb{E}[Y] = 5$; $P(X \ge 10) = 0.8^{10} = 0.1074$.

**B5.** (a) White balls = aces: $w = 4$, $b = 48$, $N = 52$, $n = 5$. $X \sim \mathrm{HGeom}(4, 48, 5)$:
$$P(X = k) = \frac{\binom4k\binom{48}{5-k}}{\binom{52}{5}}, \qquad k = 0, 1, \dots, 4.$$

(b) Let $I_j$ indicate that the $j$-th card dealt is an ace. By symmetry each card dealt is equally likely to be any of the 52, so $P(I_j = 1) = \frac{4}{52} = \frac1{13}$. By linearity, $\mathbb{E}[X] = 5/13 = 0.3846$.

(c) $\mathbb{E}[I_1 I_2] = P(\text{first two are aces}) = \frac{4}{52}\cdot\frac{3}{51} = \frac{1}{221}$ by the chain rule. Then $\operatorname{Cov}(I_1, I_2) = \frac{1}{221} - \frac{1}{169} = -\frac{12}{8619} = -0.001392$. Formula: $p = \frac1{13}$, $q = \frac{12}{13}$, $N - 1 = 51$, so $-\frac{pq}{N-1} = -\frac{12/169}{51} = -\frac{12}{8619}$. ✓

(d) $\operatorname{Var}(X) = npq\,\frac{N-n}{N-1} = 5 \cdot \frac{12}{169}\cdot\frac{47}{51} = \frac{2820}{8619} = 0.3272$. (Check via (c): $5\cdot\frac{12}{169} + 20\cdot\left(-\frac{12}{8619}\right) = \frac{3060 - 240}{8619} = \frac{2820}{8619}$.)

(e) $P(X \ge 1) = 1 - \frac{\binom{48}{5}}{\binom{52}{5}} = 1 - \frac{1{,}712{,}304}{2{,}598{,}960} = 0.3412$. The $\mathrm{Bin}(5, 1/13)$ approximation gives $P(X \ge 1) = 1 - (12/13)^5 = 0.3298$ and $\operatorname{Var} = 60/169 = 0.3550$. The mean $5/13$ is exactly right in both, because linearity uses only the marginal $P(I_j = 1) = 1/13$, which is the same with or without replacement. The variance and $P(X \ge 1)$ depend on the dependence between draws, which the Binomial ignores: without replacement, an ace drawn makes the next ace less likely (negative covariance), so the count is less spread out and more likely to be at least 1.

**B6.** (a) Substitute $u = \lambda x$, $dx = du/\lambda$:
$$\mathbb{E}[X^k] = \frac{\lambda^a}{\Gamma(a)}\int_0^\infty x^{a+k-1}e^{-\lambda x}\,dx = \frac{\lambda^a}{\Gamma(a)}\int_0^\infty \frac{u^{a+k-1}}{\lambda^{a+k-1}}e^{-u}\,\frac{du}{\lambda} = \frac{\lambda^a}{\Gamma(a)\lambda^{a+k}}\Gamma(a + k) = \frac{\Gamma(a+k)}{\Gamma(a)\lambda^k}.$$

(b) By Theorem 3.19(a), $\Gamma(a+1) = a\Gamma(a)$ and $\Gamma(a+2) = (a+1)a\Gamma(a)$. So $\mathbb{E}[X] = a/\lambda$, $\mathbb{E}[X^2] = a(a+1)/\lambda^2$, and $\operatorname{Var}(X) = \frac{a^2 + a - a^2}{\lambda^2} = \frac{a}{\lambda^2}$. At $a = 1$: $1/\lambda$ and $1/\lambda^2$, the Exponential's values from §3.11. ✓

(c) For $t < \lambda$,
$$M_X(t) = \frac{\lambda^a}{\Gamma(a)}\int_0^\infty x^{a-1}e^{-(\lambda - t)x}\,dx.$$
This is the $k = 0$ integral of (a) with $\lambda - t > 0$ in place of $\lambda$, which equals $\Gamma(a)/(\lambda - t)^a$. So $M_X(t) = \left(\frac{\lambda}{\lambda - t}\right)^a$. For $t \ge \lambda$ the integral diverges.

(d) By Theorem 3.21, each $X_i$ has MGF $\frac{\lambda}{\lambda - t}$ for $t < \lambda$. By Theorem 3.22 applied $n - 1$ times (the partial sum $X_1 + \cdots + X_{j}$ is independent of $X_{j+1}$), the sum has MGF $\left(\frac{\lambda}{\lambda - t}\right)^n$ on $t < \lambda$, an open interval around 0. By (c) that is the MGF of $\mathrm{Gamma}(n, \lambda)$, and by Theorem 3.23 the sum has that distribution. ∎

(e) $\mathrm{Gamma}(3, 2)$: $\mathbb{E}[X] = 3/2 = 1.5$, $\operatorname{Var}(X) = 3/4 = 0.75$, $\mathbb{E}[X^3] = \frac{\Gamma(6)}{\Gamma(3)\,2^3} = \frac{120}{2 \cdot 8} = 7.5$.

**B7.** (a) $F(x) = \int_{-\infty}^x \frac{dt}{\pi(1+t^2)} = \frac1\pi\big[\arctan t\big]_{-\infty}^{x} = \frac1\pi\left(\arctan x + \frac\pi2\right) = \frac12 + \frac{\arctan x}{\pi}$. Solve $u = F(x)$: $\arctan x = \pi(u - \tfrac12)$, so $F^{-1}(u) = \tan\big(\pi(u - \tfrac12)\big)$. Reference implementation (it is in Part 10 of `project.py`):

```python
def cauchy_sample(size, rng):
    """F^{-1}(u) = tan(pi (u - 1/2)) applied to U ~ Unif(0, 1)."""
    return np.tan(np.pi * (rng.random(size) - 0.5))

xc = cauchy_sample(100_000, rng)
assert stats.kstest(xc, stats.cauchy.cdf).pvalue > 1e-3
```

The KS test passes; overlay the empirical CDF exactly as Part 4 does for the Exponential.

(b) $\int_{-\infty}^\infty |x| f(x)\,dx = \frac{2}{\pi}\int_0^\infty \frac{x}{1 + x^2}dx = \frac1\pi\lim_{R\to\infty}\ln(1 + R^2) = \infty$. So $\mathbb{E}[X]$ is undefined (the positive and negative halves are each infinite; $\infty - \infty$ has no value). Part 8's three seeds end at running means $-0.015$, $0.658$, $-1.853$ after $10^5$ draws. One of them happens to land near 0, and that is the trap: with no mean there is no target, and the next huge draw can move the average anywhere. (In fact the average of $n$ i.i.d. standard Cauchy variables is again standard Cauchy, so averaging $10^5$ draws is no more concentrated than one draw. The proof uses characteristic functions, which this curriculum does not develop; Feller vol. II has it.)

(c) Quantiles converge, because they are defined through the CDF, which every distribution has, and the fraction of samples below $x$ is a Bernoulli average converging to $F(x)$. $F^{-1}(1/2) = \tan 0 = 0$ and $F^{-1}(3/4) = \tan(\pi/4) = 1$. The run gives sample median $-0.0037$ and sample upper quartile $1.0017$.

(d) The ray leaves $(0, 1)$ at angle $\Theta$ from the downward vertical and meets the line $y = 0$ where $\tan\Theta = x/1$, so $X = \tan\Theta$. Since $\tan$ is increasing on $(-\pi/2, \pi/2)$,
$$P(X \le x) = P(\Theta \le \arctan x) = \frac{\arctan x + \pi/2}{\pi} = \frac12 + \frac{\arctan x}\pi = F(x).$$
With $\Theta = \pi(U - \frac12)$ for $U \sim \mathrm{Unif}(0, 1)$, $X = \tan(\pi(U - \frac12))$: this is the sampler from (a).

**B8.** (a) Complete the square: $tz - \frac{z^2}{2} = -\frac12(z - t)^2 + \frac{t^2}{2}$. So
$$M_Z(t) = \frac{1}{\sqrt{2\pi}}\int e^{tz - z^2/2}dz = e^{t^2/2}\cdot\frac{1}{\sqrt{2\pi}}\int e^{-(z-t)^2/2}dz = e^{t^2/2},$$
the last integral being $\sqrt{2\pi}$ by Theorem 3.18 after the shift $z \mapsto z - t$. Expand:
$$e^{t^2/2} = \sum_{k\ge0}\frac{(t^2/2)^k}{k!} = \sum_{k\ge 0}\frac{t^{2k}}{2^k k!}.$$
Theorem 3.20 says the coefficient of $t^n$ is $\mathbb{E}[Z^n]/n!$. Odd powers are absent, so odd moments are 0, and
$$\mathbb{E}[Z^{2k}] = \frac{(2k)!}{2^k k!} = 1\cdot3\cdot5\cdots(2k-1).$$
$\mathbb{E}[Z^4] = \frac{24}{4\cdot2} = 3$, $\mathbb{E}[Z^6] = \frac{720}{8 \cdot 6} = 15$.

(b) Binomial expansion and linearity, with $\mathbb{E}[Z] = \mathbb{E}[Z^3] = 0$, $\mathbb{E}[Z^2] = 1$, $\mathbb{E}[Z^4] = 3$:
- $\mathbb{E}[X^2] = \mu^2 + 2\mu\sigma\,\mathbb{E}Z + \sigma^2\mathbb{E}Z^2 = \mu^2 + \sigma^2$.
- $\mathbb{E}[X^3] = \mu^3 + 3\mu^2\sigma\,\mathbb{E}Z + 3\mu\sigma^2\,\mathbb{E}Z^2 + \sigma^3\mathbb{E}Z^3 = \mu^3 + 3\mu\sigma^2$.
- $\mathbb{E}[X^4] = \mu^4 + 6\mu^2\sigma^2\,\mathbb{E}Z^2 + \sigma^4\,\mathbb{E}Z^4 = \mu^4 + 6\mu^2\sigma^2 + 3\sigma^4$ (the odd-power terms vanish).

At $\mu = 1$, $\sigma = 2$: $\mathbb{E}[X^2] = 5$, $\mathbb{E}[X^3] = 1 + 12 = 13$, $\mathbb{E}[X^4] = 1 + 24 + 48 = 73$.

(c) $X \sim \mathrm{Bin}(10, \frac12)$, so $\operatorname{Var}(X) = 10\cdot\frac12\cdot\frac12 = 2.5$. $Y = 10 - X \sim \mathrm{Bin}(10, \frac12)$, $\operatorname{Var}(Y) = 2.5$. $X + Y = 10$ is constant, so $\operatorname{Var}(X + Y) = 0$. From (3.19), $0 = 2.5 + 2.5 + 2\operatorname{Cov}(X, Y)$, so $\operatorname{Cov}(X, Y) = -2.5$.

$\operatorname{Var}(X - Y)$ two ways:
- From (3.19) with $-Y$ in place of $Y$, using $\operatorname{Cov}(X, -Y) = -\operatorname{Cov}(X, Y)$: $2.5 + 2.5 - 2(-2.5) = 10$.
- Directly: $X - Y = 2X - 10$, so $\operatorname{Var}(X - Y) = 4\operatorname{Var}(X) = 10$. ✓

The hypothesis of Theorem 3.14 that fails is independence: $Y$ is a function of $X$, so the covariance term does not vanish.

## C. Challenge

**C1.** (a) Each child has one of $2 \times 7 = 14$ equally likely (sex, weekday) types, independently, so $\Omega$ is the $14^2 = 196$ ordered pairs of types, each with probability $1/196$. Let $T$ = at least one child is a Tuesday boy. The complement "neither child is a Tuesday boy" has $13 \times 13 = 169$ outcomes, so $|T| = 196 - 169 = 27$. Among the $7 \times 7 = 49$ boy-boy outcomes, those with no Tuesday are $6 \times 6 = 36$, so $|BB \cap T| = 13$. Hence
$$P(BB \mid T) = \frac{13}{27} = 0.4815.$$

(b) A child is a "boy with the property" with probability $r/2$. So $P(T) = 1 - (1 - r/2)^2 = r - r^2/4$. Given both are boys (probability $1/4$), at least one has the property with probability $1 - (1 - r)^2 = 2r - r^2$. So
$$P(BB \mid T) = \frac{\frac14(2r - r^2)}{r - \frac{r^2}{4}} = \frac{2r - r^2}{4r - r^2} = \frac{2 - r}{4 - r}.$$
At $r = 1/7$: $\frac{13/7}{27/7} = \frac{13}{27}$. ✓ As $r \to 1$ every boy has the property, the statement reduces to "at least one boy", and the answer tends to $1/3$, Trap 3(i). As $r \to 0$ the answer tends to $1/2$, Trap 3(ii)/(iii). Words: a rare property almost singles out one specific child, and once a specific child is identified as a boy, the other is a boy with probability $1/2$. A two-boy family has two chances to contain a Tuesday boy, so it is nearly twice as likely to produce the statement, which pushes $BB$ up from $1/3$ toward $1/2$.

(c) Part 10 enumerates the 196 outcomes and gets exactly $13/27$; a $2\times10^5$-family simulation gives $0.4858$ against $0.4815$, within 1.5 standard errors.

**C2.** (a) States: $S_0$ = no progress, $S_1$ = last toss was $H$. Let $a$ and $b$ be the expected numbers of further tosses to finish from $S_0$ and $S_1$; $\mathbb{E}[W] = a$.

*$HT$.* From $S_0$: one toss; $H$ moves to $S_1$, $T$ stays at $S_0$. From $S_1$: $T$ finishes, $H$ stays in $S_1$.
$$a = 1 + \tfrac12 b + \tfrac12 a, \qquad b = 1 + \tfrac12 b.$$
The second gives $b = 2$; the first gives $\frac12 a = 1 + 1$, so $a = 4$.

*$HH$.* From $S_1$: $H$ finishes, $T$ sends you back to $S_0$.
$$a = 1 + \tfrac12 b + \tfrac12 a, \qquad b = 1 + \tfrac12 a.$$
Substitute: $\frac12 a = 1 + \frac12(1 + \frac12 a)$, so $\frac14 a = \frac32$ and $a = 6$.

(b) The difference is what a failure costs once you are halfway. Waiting for $HT$, a wrong toss in $S_1$ is an $H$, which leaves you in $S_1$: still halfway. Waiting for $HH$, a wrong toss in $S_1$ is a $T$, which sends you back to the start. Equivalently, $HH$ overlaps itself ($HHH$ contains it twice), so its occurrences come in clumps; the long-run rate is $1/4$ per position for both patterns, but clumped occurrences mean longer gaps before the first one.

(c) Same states, $H$ with probability $p$, $T$ with probability $q$.

*$HT$:* $b = 1 + pb$, so $b = 1/q$. Then $a = 1 + pb + qa$ gives $pa = 1 + p/q$, so $a = \frac1p + \frac1q$.

*$HH$:* $b = 1 + qa$ and $a = 1 + pb + qa$. Substituting, $pa = 1 + p(1 + qa)$, so $pa(1 - q) = 1 + p$, i.e. $p^2 a = 1 + p$ and $a = \frac{1}{p^2} + \frac1p$.

At $p = 1/2$ these give 4 and 6. ✓ At $p = 0.6$: $\mathbb{E}[W_{HT}] = 1.6667 + 2.5 = 4.1667$; $\mathbb{E}[W_{HH}] = 1.6667 + 2.7778 = 4.4444$.

(d) Part 10 simulates $10^5$ runs of each. Results: $4.0047$ and $5.9916$ at $p = 1/2$, $4.1623$ and $4.4415$ at $p = 0.6$, all within one standard error of the formulas.

**C3.** (a) For nonnegative integers $s, t$: $\{X \ge s + t\} \subseteq \{X \ge s\}$, so
$$G(s + t) = P(X \ge s+t) = P(X \ge s + t \mid X \ge s)\,P(X \ge s) = G(t)\,G(s).$$
By induction on $k$, $G(k) = G(k-1)G(1) = G(1)^k$, with $G(0) = 1$. Let $p = P(X = 0)$, so $G(1) = 1 - p = q$. Then $P(X = k) = G(k) - G(k+1) = q^k - q^{k+1} = q^k p$, the $\mathrm{Geom}(p)$ PMF. ∎

Where $0 < p < 1$ was used:
- $p < 1$ gives $G(s) = q^s > 0$ for every $s$, so the conditional probabilities in the hypothesis are defined.
- $p > 0$ gives $q < 1$, so $G(k) \to 0$ and the PMF sums to 1. If $p = 0$ then $G(k) = 1$ for all $k$, which no integer-valued random variable can satisfy.

(b) For a nonnegative integer $k$, $\lfloor T \rfloor \ge k$ iff $T \ge k$. So $P(\lfloor T\rfloor \ge k) = P(T \ge k) = e^{-\lambda k} = (e^{-\lambda})^k$. That is $G(k) = q^k$ with $q = e^{-\lambda}$, so $\lfloor T \rfloor \sim \mathrm{Geom}(1 - e^{-\lambda})$ by (a). It had to be so: $T$ is memoryless at all times (Theorem 3.16), in particular at integer times, so $\lfloor T \rfloor$ is memoryless on the integers, and by (a) the only such distribution is Geometric.

(c) $\lambda = \ln 2$: $p = 1 - \frac12 = \frac12$; $\mathbb{E}[X] = q/p = 1$; $P(X \ge 3) = q^3 = \frac18$. Simulation in Part 10: $1.0053$ and $0.1261$.

**C4.** (a) For $i \neq j$, $I_iI_j = 1$ iff both $i$ and $j$ are fixed, which leaves $(n-2)!$ arrangements of the rest. So $\mathbb{E}[I_iI_j] = \frac{(n-2)!}{n!} = \frac{1}{n(n-1)}$.

(b) $X^2 = \left(\sum_j I_j\right)^2 = \sum_j I_j^2 + \sum_{i \ne j} I_iI_j = \sum_j I_j + \sum_{i\neq j}I_iI_j$, using $I_j^2 = I_j$. There are $n(n-1)$ ordered pairs $i \ne j$. By linearity,
$$\mathbb{E}[X^2] = n\cdot\frac1n + n(n-1)\cdot\frac{1}{n(n-1)} = 2, \qquad \operatorname{Var}(X) = 2 - 1^2 = 1.$$

(c) $n = 4$, $4! = 24$ permutations:
- $k = 0$: the derangements of 4, numbering 9 (Example 4.4).
- $k = 1$: choose the fixed point (4 ways), derange the other 3 (2 ways): 8.
- $k = 2$: choose the pair ($\binom42 = 6$), swap the other two (1 way): 6.
- $k = 3$: 0, since three fixed points force the fourth.
- $k = 4$: 1.

The total is $9 + 8 + 6 + 0 + 1 = 24$. ✓ $\mathbb{E}[X] = \frac{0\cdot9 + 1\cdot8 + 2\cdot6 + 4\cdot1}{24} = \frac{24}{24} = 1$. $\mathbb{E}[X^2] = \frac{8 + 4\cdot6 + 16\cdot1}{24} = \frac{48}{24} = 2$. ✓

(d) For $n = 10$:
- $P(X = 0) = \sum_{i=0}^{10}\frac{(-1)^i}{i!} = 0.3678795$.
- $P(X = 1) = \sum_{i=0}^{9}\frac{(-1)^i}{i!} = 0.3678792$.
- $P(X = 2) = \frac12\sum_{i=0}^{8}\frac{(-1)^i}{i!} = 0.1839410$.

The $\mathrm{Pois}(1)$ values are $e^{-1} = 0.3678794$, $e^{-1} = 0.3678794$, $e^{-1}/2 = 0.1839397$. They agree to about $10^{-6}$, because each sum is a truncated series for $e^{-1}$. Part 10 enumerates all $n!$ permutations for $n = 4, 5, 6$ and asserts $\mathbb{E}[X] = 1$ and $\mathbb{E}[X^2] = 2$ exactly in rational arithmetic, so $\operatorname{Var}(X) = 1$.
