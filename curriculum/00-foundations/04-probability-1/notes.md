# F4: Probability I: Foundations, Random Variables, Named Distributions, Expectation

**Time:** ~14 hours · **Prerequisites:** F1 · **Unlocks:** F5, SM-3, ST-1

> **Conventions in this module.** The sample space is $\Omega$, an outcome is $\omega \in \Omega$, events are subsets $A, B \subseteq \Omega$, and $P(A)$ is the probability of $A$. Random variables are capital letters $X, Y$; the values they take are lowercase $x, y$. $\mathbb{E}[X]$ is expectation, $\operatorname{Var}(X)$ variance, $\sigma_X = \sqrt{\operatorname{Var}(X)}$ standard deviation. The PMF of a discrete $X$ is $p_X(x) = P(X = x)$; the PDF of a continuous $X$ is $f_X(x)$; the CDF of any $X$ is $F_X(x) = P(X \le x)$. $I_A$ is the indicator of $A$ ($1$ on $A$, $0$ off it). $\binom{n}{k} = \frac{n!}{k!(n-k)!}$. $q$ always means $1 - p$. $\log$ is the natural logarithm. **Geometric$(p)$ counts failures before the first success**, so its support is $\{0, 1, 2, \dots\}$ (Blitzstein's convention; SciPy's `geom` counts trials instead, see the project). **NegBin$(r, p)$ counts failures before the $r$-th success.** Exponential is parametrized by rate $\lambda$ (mean $1/\lambda$). Normal is $\mathcal{N}(\mu, \sigma^2)$ with the *variance* second. $X \perp Y$ means independent.

## 0. Where this sits

You need F1: set notation, quantifiers, induction, and the habit of writing "let $A$ be the event that..." before computing anything. You need single-variable calculus fluently (series, integration by parts, a polar-coordinates double integral); all of that you have.

This module is the probability that everything downstream shares. F5 builds joint distributions, the multivariate Gaussian, the central limit theorem, and conditional expectation on top of it. SM-3 is this module's counting applied to microstates: the multiplicity of an Einstein solid is a binomial coefficient, and Boltzmann's entropy is the logarithm of one. ST-1 takes every PMF and PDF here and reads it backwards, as a likelihood in the parameter. If the named distributions, linearity of expectation, LOTUS, and Bayes' theorem are not reflexive by the end, F5 and ST-1 will both be twice as hard as they need to be.

## 1. The problem that forced this idea

Probability started as a question about gambling debts. In 1654 the Chevalier de Méré put two problems to Blaise Pascal. One was the **problem of points**: two players stake money on a sequence of fair rounds, the first to win (say) five rounds takes the pot, and the game is interrupted at 4–3. How should the pot be split? Pascal wrote to Pierre de Fermat, and over the summer of 1654 they solved it two different ways. Fermat enumerated the possible continuations of the game and counted which ones each player wins. Pascal reasoned backwards from the final round: the value of a position is the average of the values of the positions it leads to. That backward recursion is the first use of what this module calls *expectation*, and Fermat's enumeration is the first use of *equally likely outcomes*. Christiaan Huygens learned of the correspondence and in 1657 published *De ratiociniis in ludo aleae*, the first textbook on the subject, built around the idea of the fair price of a gamble.

The next pressure came from money that was not on a table. Governments sold annuities and needed to price them; Edmond Halley published a life table in 1693 for that purpose, and Abraham de Moivre's *Doctrine of Chances* (1718) was written partly for insurers. In 1713 Jacob Bernoulli's *Ars Conjectandi* appeared posthumously with the first limit theorem: if you repeat a trial with success probability $p$ many times, the observed fraction of successes will, with probability as close to 1 as you like, be as close to $p$ as you like. This is the **weak law of large numbers**, and it is the bridge between a probability (an abstract weight) and a frequency (something you can measure). De Moivre followed in 1733 with the observation that the binomial distribution, for large $n$, is shaped like $e^{-x^2}$: the first appearance of the Normal curve. Pierre-Simon Laplace's *Théorie analytique des probabilités* (1812) systematized all of it, including the rule for updating probabilities on evidence that we now call Bayes' theorem (Thomas Bayes's own version was published posthumously in 1763).

Then physics. In 1860 James Clerk Maxwell asked what the velocities of the molecules in a gas look like, and answered with a probability distribution: the first time a law of physics was *stated as* a distribution rather than as a trajectory. In 1877 Ludwig Boltzmann connected entropy to the number of microscopic arrangements compatible with a macroscopic state, so that thermodynamics became counting. SM-3 and SM-4 are this story; they need §3.2 (counting), §3.10 (the Binomial and its Normal limit), and §3.14 (generating functions, which become partition functions).

All of this happened without anyone agreeing on what probability *is*. The frequentist position (Richard von Mises, 1919) says a probability is the limiting frequency in a long run of repetitions. The degree-of-belief position (John Maynard Keynes 1921, Frank Ramsey 1926, Bruno de Finetti 1931) says a probability is a measure of how strongly a rational agent believes something, with coherence enforced by betting odds. The two camps still argue, and the argument matters for ST-5. But in 1933 Andrey Kolmogorov's *Grundbegriffe der Wahrscheinlichkeitsrechnung* ended the argument about the *mathematics* by laying down three axioms (§3.4) that both camps accept. Whatever you think a probability means, it is a non-negative number, the whole space has probability 1, and probabilities of disjoint events add. Everything in this module follows from those three sentences plus definitions.

The last problem, inference (given data, what produced it?) was there from Bayes and Laplace onward, but it got its modern form from Fisher in the 1920s and Neyman and Pearson in the 1930s. That is ST-1 and ST-2. This module gives them their vocabulary.

## 2. Intuition

### Probability is mass

Picture $\Omega$ as a table covered with pebbles. Each pebble is one possible outcome $\omega$. Each pebble has a mass, and the total mass is 1. An **event** is a set of pebbles, and its probability is the total mass of the pebbles in it. That is the whole of §3.1 and §3.4. Disjoint events have no pebbles in common, so their masses add. The complement of $A$ has mass $1 - P(A)$ because the total is 1. If $A \subseteq B$ then $P(A) \le P(B)$ because $B$ has all of $A$'s pebbles and maybe more. Inclusion–exclusion, $P(A \cup B) = P(A) + P(B) - P(A \cap B)$, says: if you add the masses of $A$ and $B$ you have counted the shared pebbles twice, so subtract them once.

```
   Omega (total mass 1)
   +-------------------------------------+
   |  . .   .    A           .   .       |
   |    . .---------.    .      .        |
   |  .  |  . . .  .|-------.            |
   |     |   . .  . |  . .  | B          |
   |     '-----------|-.----'   .   .    |
   |   .   .    .    '------'            |
   +-------------------------------------+
   P(A u B) = mass(A) + mass(B) - mass(overlap)
```

In the **naive** setting every pebble has the same mass $1/|\Omega|$, and probability is counting: $P(A) = |A|/|\Omega|$. Most of §3.2 and §3.3 is learning to count.

### Conditioning is restricting the universe and renormalizing

"Given that $B$ happened" means: sweep every pebble outside $B$ off the table. The pebbles that remain have total mass $P(B)$, which is less than 1, so rescale every remaining mass by $1/P(B)$ to make the total 1 again. The probability of $A$ in this new, smaller universe is the mass of $A$'s pebbles that survived, rescaled:

$$P(A \mid B) = \frac{P(A \cap B)}{P(B)}.$$

Everything about conditional probability in §3.5 is this picture. Bayes' theorem is the observation that $P(A \cap B)$ can be computed from either side: restrict to $B$ and look at $A$, or restrict to $A$ and look at $B$. Independence means restricting to $B$ does not change the *proportion* of mass on $A$.

### A random variable is a function, not a number

Now write a number on each pebble. That labeling is a **random variable** $X$: a function from $\Omega$ to $\mathbb{R}$. "The sum of two dice" is a function on the 36 pebbles $(i, j)$, namely $X(i, j) = i + j$. The **distribution** of $X$ is the list of how much total mass sits on each label: $P(X = 7) = 6/36$ because six pebbles carry the label 7. Two different functions can have the same distribution (the first die and the second die, for instance) and are still different random variables. Blitzstein's phrase: the distribution is the *blueprint*, the random variable is the *house*. You will be tempted to think of $X$ as a number that hasn't been revealed yet. Resist. It is a function, and $X + Y$, $X^2$, $e^{tX}$ are all functions too, defined pebble by pebble.

### Expectation is a center of mass; variance is a moment of inertia

Put the pebbles on a number line, each at its label $x$, with its mass $P(X = x)$. The **expectation** is

$$\mathbb{E}[X] = \sum_x x \, P(X = x),$$

which is the centroid of a system of point masses with total mass 1: $\bar{x} = \sum m_i x_i / \sum m_i$. You have computed centroids since first-year statics. The **variance** is

$$\operatorname{Var}(X) = \sum_x (x - \mathbb{E}X)^2 \, P(X = x),$$

which is the second moment of mass about the centroid: $I_{\mathrm{cm}} = \sum m_i r_i^2$ with $r_i$ the distance from the center of mass. Everything you know about moments of inertia transfers:

- Shifting all the masses by $b$ moves the centroid by $b$ and does not change $I_{\mathrm{cm}}$: $\mathbb{E}[X + b] = \mathbb{E}[X] + b$, $\operatorname{Var}(X + b) = \operatorname{Var}(X)$.
- Scaling all distances by $a$ scales the centroid by $a$ and the moment by $a^2$: $\mathbb{E}[aX] = a\mathbb{E}[X]$, $\operatorname{Var}(aX) = a^2 \operatorname{Var}(X)$.
- The **parallel axis theorem** $I_{\text{axis}} = I_{\mathrm{cm}} + M d^2$, with $M = 1$ and the axis at the origin ($d = \mathbb{E}X$), reads $\mathbb{E}[X^2] = \operatorname{Var}(X) + (\mathbb{E}X)^2$. That is equation (3.17) below, the formula you will use to compute every variance in this module.

Standard deviation is the radius of gyration: the distance at which you would concentrate all the mass to get the same second moment. It has the same units as $X$.

### What you should be able to guess now

1. **Averages add.** The centroid of the labels $X(\omega) + Y(\omega)$ should be the sum of the centroids, because summing over pebbles is linear. Nothing about how $X$ and $Y$ relate should matter. (True: Theorem 3.8. Variances do *not* add in general: Theorem 3.14.)
2. **Rare events, many chances, count is Poisson.** If $n$ independent things each happen with small probability $p$, and $np = \lambda$ is moderate, the number that happen should depend only on $\lambda$. (True: §3.10.)
3. **Many small independent contributions add up to a bell curve.** De Moivre saw it in 1733. (True: the CLT, previewed in §3.12, proved in F5.)
4. **Independence is not disjointness.** Disjoint events with positive probability are as dependent as events can be: if one happens, the other certainly didn't. (§3.6.)

## 3. Formal development

### 3.1 Sample spaces, events, and the logic of sets

**Definition 3.1.** A **sample space** $\Omega$ is the set of all possible outcomes of an experiment. An **event** is a subset $A \subseteq \Omega$. The event $A$ *occurs* if the realized outcome $\omega$ lies in $A$.

Set operations are logical operations on events:

| Set language | Event language |
|---|---|
| $A \cup B$ | $A$ or $B$ (or both) occurs |
| $A \cap B$ | $A$ and $B$ both occur |
| $A^c = \Omega \setminus A$ | $A$ does not occur |
| $A \subseteq B$ | $A$ implies $B$ |
| $A \cap B = \emptyset$ | $A$ and $B$ are **disjoint**: they cannot both occur |
| $\Omega$, $\emptyset$ | the certain event, the impossible event |

De Morgan's laws, $(A \cup B)^c = A^c \cap B^c$ and $(A \cap B)^c = A^c \cup B^c$, are the probability versions of "not (A or B) = (not A) and (not B)." You will use the first one constantly in the form "the probability that *at least one* happens is one minus the probability that *none* happen."

Example: two dice. $\Omega = \{(i, j) : 1 \le i, j \le 6\}$, $|\Omega| = 36$. The event "sum is 7" is $\{(1,6), (2,5), (3,4), (4,3), (5,2), (6,1)\}$. The event "first die is 3" is $\{(3, j)\}$. Their intersection is $\{(3, 4)\}$.

### 3.2 Counting

**Multiplication rule.** If an experiment has $r$ stages with $n_1$ outcomes at stage 1, $n_2$ at stage 2 (whatever happened at stage 1), and so on, there are $n_1 n_2 \cdots n_r$ outcomes in all. Each outcome is a path through a tree with $n_1$ branches, each splitting into $n_2$, and so on.

**Permutations.** The number of orderings of $n$ distinct objects is $n! = n(n-1)\cdots 1$: $n$ choices for the first slot, $n-1$ for the second, and so on. The number of ordered selections of $k$ objects from $n$ is $n(n-1)\cdots(n-k+1) = n!/(n-k)!$.

**Binomial coefficients.** The number of $k$-element subsets of an $n$-element set is written $\binom{n}{k}$.

▶ **Try it first.** Count the ordered selections of $k$ from $n$ two ways: directly, and as "choose a subset, then order it." Solve for $\binom{n}{k}$.

**Derivation.**

1. Ordered selections of $k$ from $n$: $n!/(n-k)!$ by the multiplication rule.
2. Alternatively: choose the subset ($\binom{n}{k}$ ways), then order its $k$ elements ($k!$ ways). By the multiplication rule, $\binom{n}{k} \cdot k!$.
3. These count the same thing, so $\binom{n}{k} k! = n!/(n-k)!$, i.e.
$$\binom{n}{k} = \frac{n!}{k!\,(n-k)!}, \qquad 0 \le k \le n, \tag{3.1}$$
and $\binom{n}{k} = 0$ for $k > n$ or $k < 0$. ∎

■ **What just happened.** Every ordered selection was counted exactly $k!$ times when we ignored order, so dividing by $k!$ fixes the overcount. "Count two ways, set equal" is the method of this whole section.

**Theorem 3.1 (Binomial theorem).** For any $x, y$ and integer $n \ge 0$,
$$(x + y)^n = \sum_{k=0}^n \binom{n}{k} x^k y^{n-k}. \tag{3.2}$$

▶ **Try it first.** Expand $(x+y)(x+y)(x+y)$ without collecting terms. How many of the eight terms are $x^2 y$? Why that number?

*Proof.* $(x+y)^n$ is a product of $n$ factors. Expanding without collecting, each term is formed by choosing $x$ or $y$ from each factor, so there are $2^n$ terms. A term equals $x^k y^{n-k}$ exactly when $x$ was chosen from $k$ of the $n$ factors. The number of ways to pick *which* $k$ factors is $\binom{n}{k}$. Collecting like terms gives (3.2). ∎

■ **What just happened.** The algebra of expanding a product is counting in disguise. Set $x = y = 1$: $\sum_k \binom{n}{k} = 2^n$, the number of subsets of an $n$-set. Set $x = p$, $y = 1-p$: the Binomial PMF in §3.10 sums to 1.

**Theorem 3.2 (Three story proofs).** A *story proof* shows two expressions are equal by exhibiting one set that both expressions count.

(a) $\binom{n}{k} = \binom{n}{n-k}$.

(b) **Vandermonde's identity.** $\displaystyle \binom{m + n}{r} = \sum_{k=0}^{r} \binom{m}{k}\binom{n}{r-k}$.

(c) **Hockey-stick identity.** $\displaystyle \sum_{i=k}^{n} \binom{i}{k} = \binom{n+1}{k+1}$.

▶ **Try it first.** For (a), choosing who is *on* a committee is the same as choosing who is *off* it. For (b), split the people into two groups. For (c), choose $k+1$ numbers from $\{1, \dots, n+1\}$ and ask what the largest one is.

*Proof.* (a) Choosing a $k$-subset of $n$ people to be on a committee determines, and is determined by, the $(n-k)$-subset left off. The two sides count the same committees. ∎

(b) Choose $r$ people from a room containing $m$ men and $n$ women. The left side counts directly. The right side sorts the committees by how many men they contain: a committee with $k$ men and $r - k$ women can be formed in $\binom{m}{k}\binom{n}{r-k}$ ways (multiplication rule), and every committee has exactly one value of $k$, so summing over $k$ counts each committee once. ∎

(c) Choose a $(k+1)$-subset of $\{1, 2, \dots, n+1\}$: $\binom{n+1}{k+1}$ ways. Sort by the largest element chosen. If the largest is $i + 1$ (where $i$ ranges from $k$ to $n$, since there must be $k$ smaller elements available), the remaining $k$ elements are a $k$-subset of $\{1, \dots, i\}$: $\binom{i}{k}$ ways. Sum over $i$. ∎

■ **What just happened.** No algebra. Each identity is one set counted two ways. The habit of asking "what does this expression *count*?" is worth more than any formula in this section.

**Theorem 3.3 (Sampling table).** Choose $k$ objects from $n$ distinguishable types.

| | Order matters | Order doesn't matter |
|---|---|---|
| With replacement | $n^k$ | $\binom{n + k - 1}{k}$ |
| Without replacement | $\dfrac{n!}{(n-k)!}$ | $\dbinom{n}{k}$ |

*Proof of the top-right entry (stars and bars).* An unordered selection with replacement is a record of how many of each type you took: $(k_1, \dots, k_n)$ with $k_i \ge 0$ and $\sum k_i = k$. Encode it as $k$ stars separated into $n$ groups by $n - 1$ bars: for $n = 3$, $k = 4$, the selection "two of type 1, none of type 2, two of type 3" is $\star\star \,|\, |\, \star\star$. Every arrangement of $k$ stars and $n-1$ bars in a row corresponds to exactly one selection, and there are $\binom{n+k-1}{k}$ such arrangements (choose which $k$ of the $n+k-1$ positions hold stars). ∎ The other three entries are the multiplication rule and (3.1).

The top-right entry is the one that surprises people. In SM-3 it counts the microstates of an Einstein solid: $k$ quanta of energy distributed among $n$ oscillators, $\binom{n+k-1}{k}$ ways.

### 3.3 Naive probability and the birthday problem

**Definition 3.2 (Naive probability).** If $\Omega$ is finite and every outcome is equally likely, $P(A) = |A| / |\Omega|$.

This is the only definition Pascal and Fermat had. It is correct when the symmetry is real (fair dice, well-shuffled cards, uniformly random permutations) and wrong when it isn't. "There are two possibilities, so each has probability 1/2" is the most common probability error in the world.

**The birthday problem.** There are $k$ people in a room, each with a birthday uniformly distributed over 365 days, independently. What is the probability that at least two share a birthday?

▶ **Try it first.** Compute the probability that *no two* share a birthday by filling in the birthdays one person at a time.

**Derivation.**

1. $\Omega$ = all lists of $k$ birthdays, $|\Omega| = 365^k$, equally likely by assumption.
2. The event "all different": the first person can have any of 365 days, the second any of the remaining 364, ..., the $k$-th any of $365 - k + 1$. So $|A| = 365 \cdot 364 \cdots (365 - k + 1)$.
3. $$P(\text{no match}) = \prod_{j=0}^{k-1} \frac{365 - j}{365} = \prod_{j=0}^{k-1}\left(1 - \frac{j}{365}\right), \qquad P(\text{match}) = 1 - P(\text{no match}).$$
4. **Approximation.** For every real $x$, $1 - x \le e^{-x}$ (the function $g(x) = e^{-x} - 1 + x$ has $g(0) = 0$ and $g'(x) = 1 - e^{-x}$, which is negative for $x < 0$ and positive for $x > 0$, so $g \ge 0$). For small $x$ the inequality is nearly equality. Hence
$$P(\text{no match}) \approx \prod_{j=0}^{k-1} e^{-j/365} = \exp\!\left(-\frac{1}{365}\sum_{j=0}^{k-1} j\right) = \exp\!\left(-\frac{k(k-1)}{730}\right). \tag{3.3}$$

■ **What just happened.** The exact answer is a product of $k$ factors; the approximation turns the product into a single exponential of a sum, and the sum is a triangular number. Note $k(k-1)/2 = \binom{k}{2}$ is the number of *pairs*: (3.3) says $P(\text{no match}) \approx e^{-(\text{number of pairs})/365}$. In §3.10 you will recognize $e^{-\lambda}$ as the Poisson probability of zero events with $\lambda = \binom{k}{2}/365$. Example 4.1 does $k = 23$: exact $0.5073$, approximate $0.5000$.

### 3.4 Kolmogorov's axioms and what follows from them

**Definition 3.3 (Probability).** A probability is a function $P$ assigning to each event $A \subseteq \Omega$ a number $P(A)$ such that

$$\text{(K1) } P(A) \ge 0; \qquad \text{(K2) } P(\Omega) = 1; \qquad \text{(K3) } P\!\left(\bigcup_{i=1}^\infty A_i\right) = \sum_{i=1}^\infty P(A_i) \text{ whenever the } A_i \text{ are pairwise disjoint.} \tag{3.4}$$

(For uncountable $\Omega$ not every subset can be an event; the collection of events is a $\sigma$-algebra. That restriction never bites in this curriculum. The pointer, if you want it, is Billingsley, *Probability and Measure*, ch. 1–2.)

**Theorem 3.4 (Consequences of the axioms).** For events $A, B, C$:

1. $P(\emptyset) = 0$.
2. $P(A^c) = 1 - P(A)$.
3. If $A \subseteq B$ then $P(A) \le P(B)$. In particular $P(A) \le 1$.
4. $P(A \cup B) = P(A) + P(B) - P(A \cap B). \tag{3.5}$
5. $P(A \cup B \cup C) = P(A) + P(B) + P(C) - P(A \cap B) - P(A \cap C) - P(B \cap C) + P(A \cap B \cap C). \tag{3.6}$

▶ **Try it first.** Each one is "write a set as a disjoint union, apply (K3)." Find the disjoint union for each.

**Derivation.**

1. $\Omega = \Omega \cup \emptyset \cup \emptyset \cup \cdots$, a disjoint union. By (K3), $1 = 1 + P(\emptyset) + P(\emptyset) + \cdots$, forcing $P(\emptyset) = 0$. (Finite additivity now follows from (K3) by padding with empty sets.)
2. $\Omega = A \cup A^c$ disjointly, so $1 = P(A) + P(A^c)$.
3. $B = A \cup (B \setminus A)$ disjointly, so $P(B) = P(A) + P(B \setminus A) \ge P(A)$ by (K1). Take $B = \Omega$ for the "in particular."
4. $A \cup B = A \cup (B \setminus A)$ disjointly, so $P(A \cup B) = P(A) + P(B \setminus A)$. Also $B = (A \cap B) \cup (B \setminus A)$ disjointly, so $P(B \setminus A) = P(B) - P(A \cap B)$. Substitute.
5. Apply (3.5) to the two events $A$ and $B \cup C$:
   $P(A \cup B \cup C) = P(A) + P(B \cup C) - P(A \cap (B \cup C))$.
   Expand $P(B \cup C)$ by (3.5). Write $A \cap (B \cup C) = (A \cap B) \cup (A \cap C)$ and apply (3.5) once more: $P((A\cap B) \cup (A \cap C)) = P(A \cap B) + P(A \cap C) - P(A \cap B \cap C)$. Collect terms. ∎

■ **What just happened.** Three axioms, and the pebble picture became theorems. (3.5) is "the overlap was counted twice." (3.6) is "pairwise overlaps were counted twice, so subtract them, but then the triple overlap was added three times and subtracted three times, so add it back once."

**The general inclusion–exclusion formula** (stated; the proof is an induction on $n$ using exactly the step in item 5, or a one-line indicator argument you will be able to do after §3.9):
$$P\!\left(\bigcup_{i=1}^n A_i\right) = \sum_i P(A_i) - \sum_{i<j} P(A_i \cap A_j) + \sum_{i<j<k} P(A_i \cap A_j \cap A_k) - \cdots + (-1)^{n+1} P(A_1 \cap \cdots \cap A_n). \tag{3.7}$$

**De Montmort's matching problem (1708).** A deck of $n$ cards labeled $1, \dots, n$ is shuffled uniformly. A *match* is a card whose label equals its position. What is the probability of at least one match?

▶ **Try it first.** Let $A_i$ be the event that card $i$ is in position $i$. Compute $P(A_i)$, $P(A_i \cap A_j)$, and in general $P(A_{i_1} \cap \cdots \cap A_{i_k})$, and count how many terms of each kind appear in (3.7).

**Derivation.**

1. All $n!$ orderings are equally likely. $A_i$ fixes card $i$ and leaves the other $n - 1$ free: $P(A_i) = (n-1)!/n! = 1/n$.
2. $A_i \cap A_j$ ($i \ne j$) fixes two cards: $P(A_i \cap A_j) = (n-2)!/n! = \frac{1}{n(n-1)}$. In general a $k$-fold intersection has probability $(n-k)!/n!$.
3. There are $\binom{n}{k}$ ways to choose which $k$ events to intersect, so the $k$-th sum in (3.7) is
$$\binom{n}{k}\frac{(n-k)!}{n!} = \frac{n!}{k!(n-k)!}\cdot\frac{(n-k)!}{n!} = \frac{1}{k!}.$$
4. Therefore
$$P(\text{at least one match}) = 1 - \frac{1}{2!} + \frac{1}{3!} - \cdots + (-1)^{n+1}\frac{1}{n!} = 1 - \sum_{k=0}^{n} \frac{(-1)^k}{k!} \;\xrightarrow[n \to \infty]{}\; 1 - e^{-1} \approx 0.6321. \tag{3.8}$$ ∎

■ **What just happened.** Every term in a sum of $2^n - 1$ terms collapsed to $1/k!$ by symmetry, and the alternating series is the Taylor series of $e^{-1}$. The convergence is fast: $n = 4$ gives $0.625$, $n = 7$ already matches $1 - 1/e$ to four decimals (Example 4.4). The probability that a random shuffle leaves *no* card in place is about $37\%$ whether the deck has 10 cards or a million.

### 3.5 Conditional probability

**Definition 3.4.** For $P(B) > 0$, the **conditional probability of $A$ given $B$** is
$$P(A \mid B) = \frac{P(A \cap B)}{P(B)}. \tag{3.9}$$

This is the "restrict and renormalize" picture of §2 made exact. Check that $P(\cdot \mid B)$ is itself a probability: $P(A \mid B) \ge 0$; $P(\Omega \mid B) = P(B)/P(B) = 1$; and if the $A_i$ are disjoint then so are the $A_i \cap B$, so (K3) for $P$ gives (K3) for $P(\cdot \mid B)$. **Consequently every theorem about probabilities is also a theorem about conditional probabilities**: you may append "$\mid E$" to every term of any identity in this module. That remark is used below as "extra conditioning."

**Theorem 3.5 (Chain rule).** $P(A_1 \cap A_2 \cap \cdots \cap A_n) = P(A_1)\,P(A_2 \mid A_1)\,P(A_3 \mid A_1 \cap A_2)\cdots P(A_n \mid A_1 \cap \cdots \cap A_{n-1})$.

*Proof.* Rearrange (3.9) as $P(A \cap B) = P(B)P(A \mid B)$; apply it repeatedly, peeling off one event at a time. For $n = 3$: $P(A_1 \cap A_2 \cap A_3) = P(A_1 \cap A_2)P(A_3 \mid A_1 \cap A_2) = P(A_1)P(A_2 \mid A_1)P(A_3 \mid A_1 \cap A_2)$. ∎

**Theorem 3.6 (Bayes' theorem).** If $P(A), P(B) > 0$,
$$P(A \mid B) = \frac{P(B \mid A)\,P(A)}{P(B)}. \tag{3.10}$$

▶ **Try it first.** Write (3.9) twice, once for $P(A \mid B)$ and once for $P(B \mid A)$. Both contain $P(A \cap B)$.

*Proof.* From (3.9), $P(A \cap B) = P(A \mid B)P(B)$ and also $P(A \cap B) = P(B \mid A)P(A)$. Set them equal and divide by $P(B)$. ∎

■ **What just happened.** Two lines. The theorem is famous not because it is deep but because it reverses the direction of conditioning, and the two directions feel completely different to a human brain. $P(\text{symptom} \mid \text{disease})$ is what a textbook knows; $P(\text{disease} \mid \text{symptom})$ is what a patient wants. §3.7 shows how badly intuition handles the reversal.

**Theorem 3.7 (Law of total probability, LOTP).** Let $B_1, \dots, B_n$ be a **partition** of $\Omega$ (disjoint, union is $\Omega$) with each $P(B_i) > 0$. Then for any event $A$,
$$P(A) = \sum_{i=1}^n P(A \mid B_i)\,P(B_i). \tag{3.11}$$

*Proof.* $A = A \cap \Omega = A \cap (B_1 \cup \cdots \cup B_n) = (A \cap B_1) \cup \cdots \cup (A \cap B_n)$, and these pieces are disjoint because the $B_i$ are. By (K3), $P(A) = \sum_i P(A \cap B_i) = \sum_i P(A \mid B_i)P(B_i)$ by (3.9). ∎

■ **What just happened.** You cut $A$ into slices along the partition and added the slices. This is the formal version of "condition on what you wish you knew": if knowing which $B_i$ happened would make $P(A)$ easy, pretend you know, then average over the $B_i$. Nine problems out of ten in this module are solved by choosing the right partition. The denominator of Bayes' theorem is usually computed this way: $P(B) = P(B \mid A)P(A) + P(B \mid A^c)P(A^c)$.

**Extra conditioning.** By the remark after Definition 3.4, Bayes and LOTP hold inside any conditioning event $E$:
$$P(A \mid B, E) = \frac{P(B \mid A, E)\,P(A \mid E)}{P(B \mid E)}, \tag{3.12}$$
$$P(A \mid E) = \sum_i P(A \mid B_i, E)\,P(B_i \mid E). \tag{3.13}$$
Here $P(A \mid B, E)$ means $P(A \mid B \cap E)$. You will use (3.12) whenever evidence arrives in two installments: the posterior after the first piece becomes the prior for the second.

### 3.6 Independence

**Definition 3.5.** Events $A$ and $B$ are **independent** if
$$P(A \cap B) = P(A)\,P(B). \tag{3.14}$$
If $P(B) > 0$ this is equivalent to $P(A \mid B) = P(A)$: learning $B$ does not change the probability of $A$. Events $A_1, \dots, A_n$ are **mutually independent** if for *every* subcollection $A_{i_1}, \dots, A_{i_k}$, $P(A_{i_1} \cap \cdots \cap A_{i_k}) = P(A_{i_1})\cdots P(A_{i_k})$. For three events that is four conditions: three pairwise and one triple.

**Pairwise is not mutual.** Toss two fair coins. $A$ = first is heads, $B$ = second is heads, $C$ = the two tosses agree. Each has probability $1/2$. $P(A \cap B) = P(HH) = 1/4$, $P(A \cap C) = P(HH) = 1/4$, $P(B \cap C) = P(HH) = 1/4$: all three pairs are independent. But $P(A \cap B \cap C) = P(HH) = 1/4 \ne 1/8$. Knowing any one of them tells you nothing about any other; knowing two tells you the third with certainty.

**Independent is not disjoint.** If $A \cap B = \emptyset$ and both have positive probability, then $P(A \cap B) = 0 < P(A)P(B)$, so they are *dependent*, and about as dependent as possible: $P(A \mid B) = 0$. Students confuse these constantly because both feel like "nothing to do with each other." Disjoint means *cannot both happen*. Independent means *one happening does not change the odds of the other*. A coin landing heads and the same coin landing tails are disjoint. The first toss landing heads and the second toss landing tails are independent.

**Definition 3.6 (Conditional independence).** $A$ and $B$ are **conditionally independent given $E$** if $P(A \cap B \mid E) = P(A \mid E)\,P(B \mid E)$.

**Conditional independence does not imply independence (the two-coins example).** A bag holds a fair coin and a coin with $P(H) = 3/4$. You draw one at random and toss it twice. Let $F$ = you drew the fair coin, $H_1$, $H_2$ = the tosses land heads. Given $F$ (or given $F^c$) the tosses are independent: $P(H_1 \cap H_2 \mid F) = 1/4 = P(H_1 \mid F)P(H_2 \mid F)$, and $P(H_1 \cap H_2 \mid F^c) = 9/16$. But marginally, by LOTP,
$$P(H_1) = \tfrac12\cdot\tfrac12 + \tfrac12\cdot\tfrac34 = \tfrac58, \qquad P(H_1 \cap H_2) = \tfrac12\cdot\tfrac14 + \tfrac12\cdot\tfrac{9}{16} = \tfrac{13}{32} = 0.40625 \ne \left(\tfrac58\right)^2 = 0.390625.$$
The first toss landing heads is evidence that you hold the biased coin, which raises the odds of a second head. The tosses are dependent *through the unknown coin*. This structure (observations independent given a hidden parameter, dependent without it) is the entire setup of Bayesian inference in ST-5.

**Independence does not imply conditional independence.** Return to the two fair coins: $A$ and $B$ are independent, but given $C$ (the tosses agree), $P(A \cap B \mid C) = P(HH)/P(C) = (1/4)/(1/2) = 1/2$, while $P(A \mid C)P(B \mid C) = (1/2)(1/2) = 1/4$. Conditioning on a common consequence creates dependence.

**Definition 3.7 (Independent random variables).** Discrete $X$ and $Y$ are independent if $P(X = x, Y = y) = P(X = x)P(Y = y)$ for all $x, y$; in general, if $P(X \le x, Y \le y) = P(X \le x)P(Y \le y)$ for all $x, y$. A collection is **i.i.d.** (independent and identically distributed) if mutually independent with the same distribution. If $X \perp Y$ then $g(X) \perp h(Y)$ for any functions $g, h$: the event $\{g(X) \in S\}$ is a union of events $\{X = x\}$, and products distribute over the sums.

### 3.7 The classic traps, with full solutions

Each of these is a correct application of §3.5. Each one feels wrong. Work every one with pencil; the roadmap's prediction is that you will get at least two wrong on first contact, and that is normal.

**Trap 1: Monty Hall.** Three doors, one car, two goats. You pick door 1. The host, who knows where the car is, opens one of the other two doors to reveal a goat; if both hide goats he picks between them with equal probability. He opens door 3. Should you switch to door 2?

*Solution by Bayes.* Let $C_i$ = car behind door $i$, $P(C_i) = 1/3$. Let $O_3$ = host opens door 3. The host's rule gives $P(O_3 \mid C_1) = 1/2$ (he had a choice), $P(O_3 \mid C_2) = 1$ (door 3 was his only option), $P(O_3 \mid C_3) = 0$ (he never reveals the car). By LOTP, $P(O_3) = \frac12\cdot\frac13 + 1\cdot\frac13 + 0\cdot\frac13 = \frac12$. By Bayes,
$$P(C_1 \mid O_3) = \frac{(1/2)(1/3)}{1/2} = \frac13, \qquad P(C_2 \mid O_3) = \frac{1\cdot(1/3)}{1/2} = \frac23.$$
Switch.

*Solution by conditioning on the car.* Consider the strategy "always switch." Condition on where the car is (LOTP): if $C_1$, switching loses; if $C_2$ or $C_3$, the host is forced to reveal the other goat and switching wins. $P(\text{switch wins}) = 0\cdot\frac13 + 1\cdot\frac13 + 1\cdot\frac13 = \frac23$. The switcher wins exactly when the first guess was wrong, and the first guess is wrong with probability $2/3$.

*Why the host's strategy matters.* Suppose instead the host opens one of the two unchosen doors uniformly at random, without looking, and it happens to show a goat. Now $P(O_3 \cap \text{goat} \mid C_1) = 1/2$, $P(O_3 \cap \text{goat} \mid C_2) = 1/2$, $P(O_3 \cap \text{goat} \mid C_3) = 0$, and Bayes gives $P(C_1 \mid \cdot) = P(C_2 \mid \cdot) = 1/2$. Same observation, different answer, because the observation carries information only if the host's choice depended on the car. **Always write down the mechanism that produced the evidence.** That sentence is the whole of ST-5 in miniature.

**Trap 2: The base-rate fallacy.** A disease affects 1 in 10,000 people. A test is 99% accurate in both directions: $P(+ \mid D) = 0.99$ and $P(- \mid D^c) = 0.99$. You test positive. What is the probability you have the disease?

*Solution.* $P(D) = 10^{-4}$, $P(+ \mid D^c) = 0.01$. By Bayes with LOTP in the denominator,
$$P(D \mid +) = \frac{0.99 \times 10^{-4}}{0.99 \times 10^{-4} + 0.01 \times 0.9999} = \frac{0.000099}{0.000099 + 0.009999} = \frac{0.000099}{0.010098} \approx 0.0098.$$
About 1%. In a population of a million: 100 are sick and 99 of them test positive; 999,900 are healthy and 1% of them, 9,999 people, also test positive. Of the $99 + 9999 = 10{,}098$ positives, 99 are sick. The test is good; the disease is rare; the false positives outnumber the true positives a hundred to one. The fix is always the same: multiply the likelihood by the prior before you believe anything.

**Trap 3: The two-children problem.** A family has two children, each independently a boy or girl with probability $1/2$. $\Omega = \{BB, BG, GB, GG\}$ (elder listed first), each with probability $1/4$.

(i) "At least one is a boy." $P(BB \mid \text{at least one } B) = P(BB)/P(\{BB, BG, GB\}) = (1/4)/(3/4) = 1/3$.

(ii) "The elder is a boy." $P(BB \mid \text{elder } B) = P(BB)/P(\{BB, BG\}) = (1/4)/(1/2) = 1/2$.

(iii) "You meet one of the two children, chosen at random, and it is a boy." Now the sample space includes which child you met. $P(\text{met a boy}) = 1/2$ by symmetry; $P(BB \cap \text{met a boy}) = P(BB) = 1/4$; so $P(BB \mid \text{met a boy}) = 1/2$.

Three phrasings, two answers. The information "at least one boy" in (i) is *less* than the information in (ii) or (iii), where a specific child has been identified; less information, more uncertainty, lower probability of $BB$. Exercise C1 adds "born on a Tuesday" and gets $13/27$.

**Trap 4: Simpson's paradox.** A treatment can be better in every subgroup and worse overall. The standard example is the kidney-stone data of Charig et al. (1986), as usually quoted:

| | Small stones | Large stones | All patients |
|---|---|---|---|
| Treatment A (open surgery) | $81/87 = 93\%$ | $192/263 = 73\%$ | $273/350 = 78\%$ |
| Treatment B (percutaneous) | $234/270 = 87\%$ | $55/80 = 69\%$ | $289/350 = 83\%$ |

A beats B on small stones (93 vs 87) and on large stones (73 vs 69), and loses overall (78 vs 83). No arithmetic error: it is LOTP. With $S$ = success and $L$ = large stone,
$$P(S \mid A) = P(S \mid A, L)P(L \mid A) + P(S \mid A, L^c)P(L^c \mid A),$$
and the *weights* $P(L \mid A) = 263/350 = 0.75$ versus $P(L \mid B) = 80/350 = 0.23$ differ enormously: surgeons sent the hard cases to A. The overall rate is a weighted average, and A's average is dragged down by being weighted toward the hard stratum. The fix: never compare marginal rates when the groups were assigned by something that also affects the outcome. ST-9 (causal inference) is about when the stratified comparison is the right one.

**Trap 5: The prosecutor's fallacy.** A forensic match has probability one in a million for a random innocent person. The prosecutor says: "the probability the defendant is innocent is one in a million." That confuses $P(\text{match} \mid \text{innocent}) = 10^{-6}$ with $P(\text{innocent} \mid \text{match})$. If the defendant was found by searching a database of a city of 10 million with no other evidence, the prior is $P(\text{guilty}) = 10^{-7}$ and Bayes gives
$$P(\text{guilty} \mid \text{match}) = \frac{1 \cdot 10^{-7}}{1\cdot 10^{-7} + 10^{-6}(1 - 10^{-7})} \approx \frac{10^{-7}}{10^{-7} + 10^{-6}} = \frac{1}{11}.$$
About ten innocent people in the city also match. The 1999 conviction of Sally Clark in England rested partly on a 1-in-73-million figure for two infant deaths in one family being presented as the probability of her innocence (the figure also wrongly assumed independence of the deaths); the Royal Statistical Society protested publicly in 2001 and the conviction was quashed in 2003. The fix is the same as Trap 2: the likelihood is not the posterior.
