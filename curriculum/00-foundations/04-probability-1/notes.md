# F4: Probability I: Foundations, Random Variables, Named Distributions, Expectation

**Time:** ~14 hours · **Prerequisites:** F1 · **Unlocks:** F5, SM-3, ST-1

> **Conventions in this module.** The sample space is $\Omega$, an outcome is $\omega \in \Omega$, events are subsets $A, B \subseteq \Omega$, and $P(A)$ is the probability of $A$. Random variables are capital letters $X, Y$; the values they take are lowercase $x, y$. $\mathbb{E}[X]$ is expectation, $\operatorname{Var}(X)$ variance, $\sigma_X = \sqrt{\operatorname{Var}(X)}$ standard deviation. The PMF of a discrete $X$ is $p_X(x) = P(X = x)$; the PDF of a continuous $X$ is $f_X(x)$; the CDF of any $X$ is $F_X(x) = P(X \le x)$. $I_A$ is the indicator of $A$ ($1$ on $A$, $0$ off it). $\binom{n}{k} = \frac{n!}{k!(n-k)!}$. $q$ always means $1 - p$. $\log$ is the natural logarithm. **Geometric$(p)$ counts failures before the first success**, so its support is $\{0, 1, 2, \dots\}$ (Blitzstein's convention; SciPy's `geom` counts trials instead, see the project). **NegBin$(r, p)$ counts failures before the $r$-th success.** Exponential is parametrized by rate $\lambda$ (mean $1/\lambda$). Normal is $\mathcal{N}(\mu, \sigma^2)$ with the *variance* second. $X \perp Y$ means independent.

## 0. Where this sits

You need F1: set notation, quantifiers, induction, and the habit of writing "let $A$ be the event that..." before computing anything. You need single-variable calculus fluently (series, integration by parts, a polar-coordinates double integral); all of that you have.

This module is the probability that everything downstream shares. F5 builds joint distributions, the multivariate Gaussian, the central limit theorem, and conditional expectation on top of it. SM-3 is this module's counting applied to microstates: the multiplicity of an Einstein solid is a binomial coefficient, and Boltzmann's entropy is the logarithm of one. ST-1 takes every PMF and PDF here and reads it backwards, as a likelihood in the parameter. If the named distributions, linearity of expectation, LOTUS, and Bayes' theorem are not reflexive by the end, F5 and ST-1 will both be twice as hard as they need to be.

## 1. The problem that forced this idea

Probability started as a question about gambling debts. In 1654 the Chevalier de Méré put two problems to Blaise Pascal. One was the **problem of points**: two players stake money on a sequence of fair rounds, the first to win (say) five rounds takes the pot, and the game is interrupted at 4–3. How should the pot be split? Pascal wrote to Pierre de Fermat, and over the summer of 1654 they solved it two different ways. Fermat enumerated the possible continuations of the game and counted which ones each player wins. Pascal reasoned backwards from the final round: the value of a position is the average of the values of the positions it leads to. That backward recursion is the first use of what this module calls *expectation*, and Fermat's enumeration is the first use of *equally likely outcomes*. Christiaan Huygens learned of the correspondence and in 1657 published *De ratiociniis in ludo aleae*, the first textbook on the subject, built around the idea of the fair price of a gamble.

The next pressure came from money that was not on a table. Governments sold annuities and needed to price them; Edmond Halley published a life table in 1693 for that purpose, and Abraham de Moivre, whose *Doctrine of Chances* (1718) was the standard text on games, turned the same methods on life annuities in *Annuities upon Lives* (1725). In 1713 Jacob Bernoulli's *Ars Conjectandi* appeared posthumously with the first limit theorem: if you repeat a trial with success probability $p$ many times, the observed fraction of successes will, with probability as close to 1 as you like, be as close to $p$ as you like. This is the **weak law of large numbers**, and it is the bridge between a probability (an abstract weight) and a frequency (something you can measure). De Moivre followed in 1733 with the observation that the binomial distribution, for large $n$, is shaped like $e^{-x^2}$: the first appearance of the Normal curve. Pierre-Simon Laplace's *Théorie analytique des probabilités* (1812) systematized all of it, including the rule for updating probabilities on evidence that we now call Bayes' theorem (Thomas Bayes's own version was published posthumously in 1763).

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

■ **What just happened.** Every term in a sum of $2^n - 1$ terms collapsed to $1/k!$ by symmetry, and the alternating series is the Taylor series of $e^{-1}$. The convergence is fast: $n = 4$ gives $0.625$, $n = 7$ already matches $1 - 1/e$ to four decimals (the project prints every $n$ up to 10). The probability that a random shuffle leaves *no* card in place is about $37\%$ whether the deck has 10 cards or a million.

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

### 3.8 Random variables, PMF, CDF

**Definition 3.8.** A **random variable** is a function $X : \Omega \to \mathbb{R}$. It is **discrete** if it takes values in a finite or countably infinite set. The event $\{X = x\}$ means $\{\omega \in \Omega : X(\omega) = x\}$.

**Definition 3.9.** The **probability mass function (PMF)** of a discrete $X$ is $p_X(x) = P(X = x)$. The **support** of $X$ is $\{x : p_X(x) > 0\}$. The **cumulative distribution function (CDF)** of any random variable is
$$F_X(x) = P(X \le x), \qquad x \in \mathbb{R}. \tag{3.15}$$

A function is a valid PMF iff $p(x) \ge 0$ for all $x$ and $\sum_x p(x) = 1$; the sum is over the support. Both conditions are (K1) and (K2) applied to the partition $\{X = x\}$. The PMF of a discrete $X$ determines its CDF ($F(x) = \sum_{t \le x} p(t)$, a staircase that jumps by $p(x)$ at each support point) and the CDF determines the PMF (the jump sizes).

**Properties of any CDF.** (i) $F$ is nondecreasing, because $x \le y$ gives $\{X \le x\} \subseteq \{X \le y\}$ and Theorem 3.4(3) applies. (ii) $F(x) \to 0$ as $x \to -\infty$ and $F(x) \to 1$ as $x \to +\infty$. (iii) $F$ is right-continuous: $F(x_n) \to F(x)$ as $x_n \downarrow x$. (The last two use the continuity of $P$ along monotone sequences of events, which follows from (K3); Blitzstein §3.6.) And $P(a < X \le b) = F(b) - F(a)$, since $\{X \le b\} = \{X \le a\} \cup \{a < X \le b\}$ disjointly.

**Blueprint versus house.** Let $X$ be the indicator of heads on one fair toss and $Y = 1 - X$ the indicator of tails. Both have PMF $p(0) = p(1) = 1/2$: the same distribution. But $X + Y = 1$ always and $X - Y = \pm 1$ always, so $X$ and $Y$ are different functions. A distribution (the blueprint) tells you the probabilities of the values. A random variable (the house) is the specific function on $\Omega$. Questions about one variable are answered by its distribution. Questions about two or more variables together need the joint distribution, which is F5's subject; in this module we only use the joint structure through independence.

### 3.9 Expectation

Throughout this subsection $\Omega$ is finite or countable and $X$, $Y$ are random variables on it; all sums are assumed to converge absolutely, so they may be rearranged. (For a random variable with countable range on an arbitrary $\Omega$, every argument below goes through verbatim with "$\sum_{\omega : X(\omega) = x} P(\{\omega\})$" replaced by "$P(X = x)$"; for continuous variables see §3.11.)

**Definition 3.10.** The **expectation** of $X$ is
$$\mathbb{E}[X] = \sum_{\omega \in \Omega} X(\omega)\,P(\{\omega\}) = \sum_{x} x\,P(X = x). \tag{3.16}$$

The two expressions are equal by **grouping**: collect the pebbles $\omega$ by their label $x$. Formally, $\sum_\omega X(\omega)P(\{\omega\}) = \sum_x \sum_{\omega : X(\omega) = x} x\,P(\{\omega\}) = \sum_x x \sum_{\omega : X(\omega) = x} P(\{\omega\}) = \sum_x x\,P(X = x)$, using (K3) for the inner sum. The first form (sum over pebbles) is the one that makes proofs easy. The second (sum over values) is the one you compute with. It is the center of mass of §2.

**Theorem 3.8 (Linearity of expectation).** For random variables $X, Y$ and constants $a, b$:
$$\mathbb{E}[X + Y] = \mathbb{E}[X] + \mathbb{E}[Y], \qquad \mathbb{E}[aX + b] = a\,\mathbb{E}[X] + b.$$
**No assumption of independence is made.**

▶ **Try it first.** Use the pebble form of (3.16). $X + Y$ is the function $\omega \mapsto X(\omega) + Y(\omega)$.

*Proof.* $\mathbb{E}[X + Y] = \sum_\omega \big(X(\omega) + Y(\omega)\big)P(\{\omega\}) = \sum_\omega X(\omega)P(\{\omega\}) + \sum_\omega Y(\omega)P(\{\omega\}) = \mathbb{E}[X] + \mathbb{E}[Y]$. Likewise $\sum_\omega (aX(\omega) + b)P(\{\omega\}) = a\sum_\omega X(\omega)P(\{\omega\}) + b\sum_\omega P(\{\omega\}) = a\mathbb{E}[X] + b$. ∎

■ **What just happened.** A sum of sums is a sum. The proof never mentioned how $X$ and $Y$ relate because it worked pebble by pebble, where $X$ and $Y$ are both fixed numbers. This is the single most-used theorem in the module: every "expected number of ..." in §3.9–3.10 is linearity plus indicators. Students who try to prove it from the value form (3.16) get tangled in joint PMFs and start believing independence is needed. It isn't.

**Theorem 3.9 (The fundamental bridge).** For any event $A$, $\mathbb{E}[I_A] = P(A)$.

*Proof.* $I_A$ takes the value 1 with probability $P(A)$ and 0 with probability $1 - P(A)$. By (3.16), $\mathbb{E}[I_A] = 1\cdot P(A) + 0\cdot(1 - P(A)) = P(A)$. ∎

Combined with linearity, this converts "the expected number of things that happen" into "the sum of the probabilities that each thing happens," whether or not the things are independent. Useful algebra: $I_A^2 = I_A$, $I_{A \cap B} = I_A I_B$, $I_{A^c} = 1 - I_A$, $I_{A \cup B} = 1 - (1 - I_A)(1 - I_B)$. (Expanding the last for $n$ events and taking expectations proves (3.7) in one line; Exercise B1 asks you to do it for $n = 3$.)

**Theorem 3.10 (LOTUS: Law of the Unconscious Statistician).** For any function $g : \mathbb{R} \to \mathbb{R}$,
$$\mathbb{E}[g(X)] = \sum_x g(x)\,P(X = x).$$

*Proof.* $g(X)$ is the random variable $\omega \mapsto g(X(\omega))$. By the pebble form, $\mathbb{E}[g(X)] = \sum_\omega g(X(\omega))P(\{\omega\})$. Group the pebbles by the value of $X$: $= \sum_x g(x)\sum_{\omega : X(\omega) = x} P(\{\omega\}) = \sum_x g(x) P(X = x)$. ∎ The same grouping argument gives the two-variable version $\mathbb{E}[g(X, Y)] = \sum_{x, y} g(x, y)\,P(X = x, Y = y)$.

■ **What just happened.** You do not need the distribution of $g(X)$ to compute its expectation. You weight $g(x)$ by the probabilities of $x$. The name is a joke about people who do this without noticing it needs proof. It needed proof, and the proof was grouping. The roadmap flags LOTUS as a skill to drill; the exercises do.

**Definition 3.11.** The **variance** of $X$ with mean $\mu = \mathbb{E}[X]$ is $\operatorname{Var}(X) = \mathbb{E}[(X - \mu)^2]$. The **standard deviation** is $\sigma_X = \sqrt{\operatorname{Var}(X)}$.

Variance has the units of $X^2$; standard deviation has the units of $X$. If $X$ is a length in meters, $\operatorname{Var}(X)$ is in $\mathrm{m}^2$ and reporting it as a spread is a units error. Always report $\sigma$.

**Theorem 3.11.**
$$\operatorname{Var}(X) = \mathbb{E}[X^2] - (\mathbb{E}X)^2. \tag{3.17}$$

▶ **Try it first.** Expand $(X - \mu)^2$ and apply linearity, remembering that $\mu$ is a constant.

*Proof.* $(X - \mu)^2 = X^2 - 2\mu X + \mu^2$. By linearity, $\mathbb{E}[(X-\mu)^2] = \mathbb{E}[X^2] - 2\mu\,\mathbb{E}[X] + \mu^2 = \mathbb{E}[X^2] - 2\mu^2 + \mu^2 = \mathbb{E}[X^2] - \mu^2$. ∎

■ **What just happened.** This is the parallel axis theorem of §2: the second moment about the origin, $\mathbb{E}[X^2]$, equals the second moment about the centroid, $\operatorname{Var}(X)$, plus (mass 1)$\times$(distance from origin to centroid)$^2$. In practice you compute $\mathbb{E}[X^2]$ by LOTUS and subtract $\mu^2$.

**Theorem 3.12.** $\operatorname{Var}(aX + b) = a^2\operatorname{Var}(X)$. In particular $\operatorname{Var}(X + b) = \operatorname{Var}(X)$ and $\operatorname{Var}(c) = 0$ for a constant $c$.

*Proof.* $\mathbb{E}[aX + b] = a\mu + b$, so $(aX + b) - (a\mu + b) = a(X - \mu)$ and $\operatorname{Var}(aX + b) = \mathbb{E}[a^2(X - \mu)^2] = a^2\operatorname{Var}(X)$. ∎

Also $\operatorname{Var}(X) \ge 0$ always, with equality iff $X$ equals $\mu$ with probability 1 (a sum of non-negative terms is zero iff each term is zero).

**Theorem 3.13 (Products of independent variables).** If $X \perp Y$ then $\mathbb{E}[XY] = \mathbb{E}[X]\,\mathbb{E}[Y]$, and more generally $\mathbb{E}[g(X)h(Y)] = \mathbb{E}[g(X)]\,\mathbb{E}[h(Y)]$.

*Proof.* By two-variable LOTUS and Definition 3.7, $\mathbb{E}[g(X)h(Y)] = \sum_{x,y} g(x)h(y)P(X = x, Y = y) = \sum_{x,y} g(x)h(y)P(X = x)P(Y = y) = \Big(\sum_x g(x)P(X = x)\Big)\Big(\sum_y h(y)P(Y = y)\Big)$. ∎

**Definition 3.12 (Covariance, minimally).** $\operatorname{Cov}(X, Y) = \mathbb{E}[(X - \mu_X)(Y - \mu_Y)]$. Expanding exactly as in Theorem 3.11,
$$\operatorname{Cov}(X, Y) = \mathbb{E}[XY] - \mathbb{E}[X]\,\mathbb{E}[Y]. \tag{3.18}$$
By Theorem 3.13, **independent variables have covariance zero.** (The converse is false: Part 9 of the project has $X$ uniform on $\{-1, 0, 1\}$ and $Y = X^2$, dependent with $\operatorname{Cov} = 0$. F5 treats covariance properly.)

**Theorem 3.14 (Variance of a sum).**
$$\operatorname{Var}(X + Y) = \operatorname{Var}(X) + \operatorname{Var}(Y) + 2\operatorname{Cov}(X, Y), \tag{3.19}$$
and if $X_1, \dots, X_n$ are pairwise independent, $\operatorname{Var}(X_1 + \cdots + X_n) = \operatorname{Var}(X_1) + \cdots + \operatorname{Var}(X_n)$.

▶ **Try it first.** Write $X + Y - (\mu_X + \mu_Y) = (X - \mu_X) + (Y - \mu_Y)$ and square.

*Proof.* $\operatorname{Var}(X+Y) = \mathbb{E}\big[((X - \mu_X) + (Y - \mu_Y))^2\big] = \mathbb{E}[(X-\mu_X)^2] + \mathbb{E}[(Y - \mu_Y)^2] + 2\mathbb{E}[(X-\mu_X)(Y-\mu_Y)]$ by linearity, which is (3.19). For $n$ variables the same expansion produces every $\operatorname{Var}(X_i)$ once and every $\operatorname{Cov}(X_i, X_j)$, $i \ne j$, twice; pairwise independence kills the covariances. ∎

■ **What just happened.** Expectations add always; variances add only when the cross terms vanish. In the moment-of-inertia picture: the second moment of a combined body is not the sum of the parts' second moments unless the parts are arranged symmetrically about the center. Exercise B8 has a pair with $\operatorname{Var}(X + Y) = 0$ although each variance is positive.

**The indicator trick: three applications.** Write the count as a sum of indicators, take expectations, use the bridge.

*(1) Fixed points of a random permutation.* Shuffle $n$ cards (de Montmort again). Let $X$ be the number of matches and $I_j$ the indicator that card $j$ is in position $j$. Then $X = \sum_{j=1}^n I_j$ and
$$\mathbb{E}[X] = \sum_{j=1}^n P(\text{card } j \text{ in position } j) = n \cdot \frac1n = 1 \quad \text{for every } n.$$
The $I_j$ are dependent (if $n - 1$ cards are in place, so is the last), and it does not matter.

*(2) Birthday matches.* With $k$ people, let $X$ be the number of pairs sharing a birthday and $I_{ij}$ the indicator that persons $i$ and $j$ share one. $P(I_{ij} = 1) = 1/365$ (whatever $i$'s birthday is, $j$ matches it with probability $1/365$). There are $\binom{k}{2}$ pairs, so $\mathbb{E}[X] = \binom{k}{2}/365$. For $k = 23$: $253/365 = 0.693$. This is the exponent in (3.3).

*(3) The Binomial mean.* $n$ independent trials each succeeding with probability $p$; $X = \sum_{j=1}^n I_j$ with $I_j$ the indicator of success on trial $j$. $\mathbb{E}[X] = \sum_j p = np$. No PMF needed. §3.10 derives it from the PMF too, so you can see how much work the bridge saves.

**The St. Petersburg paradox (Nicolas Bernoulli, 1713; Daniel Bernoulli, 1738).** Toss a fair coin until the first head. If it appears on toss $n$, you win $2^n$ dollars. The fair price of the game is its expectation:
$$\mathbb{E}[\text{payout}] = \sum_{n=1}^\infty 2^n \cdot \frac{1}{2^n} = \sum_{n=1}^\infty 1 = \infty.$$
Nobody would pay more than a few dollars to play. The resolution is not that the mathematics is wrong: the center of mass of this distribution really is at infinity, because the arm extends forever with masses that decay exactly as fast as the lever lengthens. The resolutions are (a) utility is not linear in money (Daniel Bernoulli proposed $\log$), and (b) any real casino has finite reserves: if the payout is capped at $2^{40}$ dollars (about a trillion), the expectation drops to $\sum_{n=1}^{40} 1 + \sum_{n > 40} 2^{40}/2^n = 40 + 1 = 41$ dollars. Expectation is a tool, and this is the example to remember when someone reports a mean of a heavy-tailed quantity as if it summarized anything. The Cauchy distribution in the project (Part 8) has no mean at all, and its running average never settles.

### 3.10 The named discrete distributions

Each distribution below comes with a **story**: an experiment that produces it. Learn the story first. The PMF, mean, and variance follow from the story, and so does the recognition of the distribution when it appears in SM-3 or ST-1 under a different name. A summary table closes the subsection.

#### Bernoulli

**Story.** One trial that succeeds with probability $p$. $X = 1$ on success, $0$ on failure. Written $X \sim \mathrm{Bern}(p)$. Every indicator $I_A$ is $\mathrm{Bern}(P(A))$.

**PMF.** $P(X = 1) = p$, $P(X = 0) = q = 1 - p$.

**Mean and variance.** $\mathbb{E}[X] = p$ (the bridge). $X^2 = X$, so $\mathbb{E}[X^2] = p$ and by (3.17) $\operatorname{Var}(X) = p - p^2 = pq$. Maximum variance $1/4$ at $p = 1/2$: a fair coin is the most unpredictable.

#### Binomial

**Story.** $n$ independent $\mathrm{Bern}(p)$ trials; $X$ = number of successes. $X \sim \mathrm{Bin}(n, p)$.

▶ **Try it first.** What is the probability of one *specific* sequence with $k$ successes and $n - k$ failures? How many such sequences are there?

**PMF derived.** A specific sequence of outcomes with $k$ successes and $n-k$ failures has probability $p^k q^{n-k}$ by independence (Definition 3.5, applied $n$-fold). The number of such sequences is the number of ways to choose which $k$ of the $n$ positions hold the successes, $\binom{n}{k}$. These sequences are disjoint events, so
$$P(X = k) = \binom{n}{k} p^k q^{n-k}, \qquad k = 0, 1, \dots, n. \tag{3.20}$$
It sums to 1 by the binomial theorem (3.2) with $x = p$, $y = q$: $(p + q)^n = 1$.

**Mean, two ways.** *Via indicators:* $\mathbb{E}[X] = np$ (indicator trick, application 3). *Via the PMF:* using the identity $k\binom{n}{k} = n\binom{n-1}{k-1}$ (story: choose a $k$-person committee and its chair; either pick the committee then the chair, $\binom{n}{k}k$ ways, or the chair then the other $k-1$ members, $n\binom{n-1}{k-1}$ ways),
$$\mathbb{E}[X] = \sum_{k=0}^n k\binom{n}{k}p^k q^{n-k} = \sum_{k=1}^n n\binom{n-1}{k-1}p^k q^{n-k} = np\sum_{j=0}^{n-1}\binom{n-1}{j}p^j q^{(n-1)-j} = np\,(p + q)^{n-1} = np,$$
where $j = k - 1$ and the last sum is (3.2) again. Same answer, five times the work. That is the lesson.

**Variance via independence.** $X = \sum_{j=1}^n I_j$ with independent $I_j \sim \mathrm{Bern}(p)$, so by Theorem 3.14, $\operatorname{Var}(X) = \sum_j pq = npq$.

**Facts you should be able to prove from the story.** If $X \sim \mathrm{Bin}(n, p)$ then $n - X \sim \mathrm{Bin}(n, q)$ (relabel success and failure). If $X \sim \mathrm{Bin}(n, p)$ and $Y \sim \mathrm{Bin}(m, p)$ are independent then $X + Y \sim \mathrm{Bin}(n + m, p)$ (concatenate the trials; §3.14 proves it again with MGFs).

#### Geometric

**Story.** Independent $\mathrm{Bern}(p)$ trials, $0 < p \le 1$, repeated until the first success; $X$ = number of *failures* before it. $X \sim \mathrm{Geom}(p)$, support $\{0, 1, 2, \dots\}$. (The number of *trials* including the success is $X + 1$, sometimes called the First Success distribution; it has mean $1/p$ and the same variance.)

**PMF derived.** $X = k$ means $k$ failures then a success, one specific sequence: $P(X = k) = q^k p$. It sums to 1 by the geometric series: $\sum_{k\ge0} q^k p = p/(1 - q) = 1$. $\tag{3.21}$

**Mean, two ways.** *Via the series.* Differentiate $\sum_{k \ge 0} q^k = \frac{1}{1-q}$ with respect to $q$ (allowed inside the radius of convergence): $\sum_{k\ge1} kq^{k-1} = \frac{1}{(1-q)^2}$. Then
$$\mathbb{E}[X] = \sum_{k \ge 0} k q^k p = pq\sum_{k\ge1} k q^{k-1} = \frac{pq}{(1-q)^2} = \frac{pq}{p^2} = \frac{q}{p}.$$
*Via first-step conditioning.* Condition on the first trial (LOTP for expectations, which is linearity applied to $X = X I_{S} + X I_{S^c}$ with $S$ = first trial succeeds). If it succeeds, $X = 0$. If it fails, we have used one failure and the remaining trials are a fresh sequence of independent $\mathrm{Bern}(p)$ trials, so the number of *further* failures is again $\mathrm{Geom}(p)$ and $X = 1 + X'$ with $X' \sim \mathrm{Geom}(p)$. Therefore
$$\mathbb{E}[X] = p\cdot 0 + q\,(1 + \mathbb{E}[X]) \;\Longrightarrow\; \mathbb{E}[X](1 - q) = q \;\Longrightarrow\; \mathbb{E}[X] = \frac{q}{p}.$$
The series proof is calculus. The first-step proof is one line of algebra once you see that the process restarts. First-step analysis is how Pascal solved the problem of points, how F5 computes hitting times of Markov chains, and how Exercise C2 finds the expected wait for $HH$.

**Variance.** Same first-step trick on $X^2$: $\mathbb{E}[X^2] = p\cdot 0 + q\,\mathbb{E}[(1 + X')^2] = q\,(1 + 2\mathbb{E}[X] + \mathbb{E}[X^2])$, so $\mathbb{E}[X^2](1 - q) = q + 2q^2/p$, giving $\mathbb{E}[X^2] = \frac{q}{p} + \frac{2q^2}{p^2} = \frac{pq + 2q^2}{p^2}$. By (3.17),
$$\operatorname{Var}(X) = \frac{pq + 2q^2}{p^2} - \frac{q^2}{p^2} = \frac{pq + q^2}{p^2} = \frac{q(p + q)}{p^2} = \frac{q}{p^2}.$$

#### Negative Binomial

**Story.** Independent $\mathrm{Bern}(p)$ trials until the $r$-th success; $X$ = number of failures before it. $X \sim \mathrm{NBin}(r, p)$. $\mathrm{NBin}(1, p) = \mathrm{Geom}(p)$.

**PMF derived.** $X = n$ means the $(n + r)$-th trial is a success and exactly $r - 1$ of the first $n + r - 1$ trials are successes. There are $\binom{n + r - 1}{r - 1}$ ways to place those $r-1$ successes, each full sequence having probability $p^r q^n$:
$$P(X = n) = \binom{n + r - 1}{r - 1} p^r q^n, \qquad n = 0, 1, 2, \dots \tag{3.22}$$

**Mean and variance.** Let $G_i$ be the number of failures between the $(i-1)$-th and $i$-th success. By the fresh-start argument each $G_i \sim \mathrm{Geom}(p)$, and they are independent (they depend on disjoint blocks of trials). $X = G_1 + \cdots + G_r$, so by linearity and Theorem 3.14,
$$\mathbb{E}[X] = \frac{rq}{p}, \qquad \operatorname{Var}(X) = \frac{rq}{p^2}.$$
Deriving these from (3.22) directly is miserable. The story makes them two lines.

#### Hypergeometric

**Story.** An urn holds $w$ white and $b$ black balls, $N = w + b$. Draw $n$ **without replacement**; $X$ = number of white balls drawn. $X \sim \mathrm{HGeom}(w, b, n)$. Quality control: $w$ defective parts in a batch of $N$, inspect $n$ (Example 4.9). Ecology: tag $w$ fish, recapture $n$.

**PMF derived.** All $\binom{N}{n}$ unordered samples are equally likely. The samples with exactly $k$ white balls are formed by choosing $k$ of the $w$ white and $n - k$ of the $b$ black:
$$P(X = k) = \frac{\binom{w}{k}\binom{b}{n-k}}{\binom{N}{n}}, \qquad \max(0, n - b) \le k \le \min(w, n). \tag{3.23}$$
It sums to 1 by Vandermonde's identity, Theorem 3.2(b), with $m = w$, $r = n$.

**Mean via indicators.** Let $I_j$ be the indicator that the $j$-th ball drawn is white. By symmetry the $j$-th ball drawn is equally likely to be any of the $N$ balls (imagine drawing all $N$ in a random order and reading off position $j$), so $P(I_j = 1) = w/N$. Then $\mathbb{E}[X] = \sum_{j=1}^n w/N = n\,\frac{w}{N}$. The indicators are dependent (drawing a white ball first makes a white second ball less likely), and linearity does not care.

**Variance and the finite-population correction.** Write $p = w/N$. By Theorem 3.14 with covariances kept,
$$\operatorname{Var}(X) = \sum_j \operatorname{Var}(I_j) + \sum_{i \ne j}\operatorname{Cov}(I_i, I_j) = n\,pq + n(n-1)\operatorname{Cov}(I_1, I_2),$$
since all pairs have the same covariance by symmetry. Now $\mathbb{E}[I_1 I_2] = P(\text{first two both white}) = \frac{w}{N}\cdot\frac{w - 1}{N - 1}$ (chain rule), so by (3.18)
$$\operatorname{Cov}(I_1, I_2) = \frac{w(w-1)}{N(N-1)} - \frac{w^2}{N^2} = \frac{w(w-1)N - w^2(N-1)}{N^2(N-1)} = \frac{w^2 - wN}{N^2(N-1)} = -\frac{w(N - w)}{N^2(N-1)} = -\frac{pq}{N - 1}.$$
Negative, as expected. Hence
$$\operatorname{Var}(X) = npq - \frac{n(n-1)pq}{N-1} = npq\,\frac{N - n}{N - 1}.$$
The factor $\frac{N-n}{N-1}$ is the **finite-population correction**. It is 1 when $n = 1$, it is 0 when $n = N$ (you drew everything; no randomness), and it tends to 1 as $N \to \infty$ with $n$ fixed. In that limit sampling without replacement is indistinguishable from sampling with replacement and $\mathrm{HGeom}(w, b, n) \to \mathrm{Bin}(n, p)$: each draw barely changes the urn. Example 4.9 compares the two for a batch of 50 parts, where the correction is far from 1.

#### Poisson

**Story.** Count the number of events when there are *many* opportunities, each with *small* probability, roughly independent. Raindrops on a tile in a second; typos on a page; radioactive decays in a minute; birthday-sharing pairs in a room. The parameter $\lambda$ is the expected count.

▶ **Try it first.** Take $X_n \sim \mathrm{Bin}(n, \lambda/n)$ and compute $\lim_{n\to\infty} P(X_n = k)$ for fixed $k$. You will need $\lim (1 - \lambda/n)^n = e^{-\lambda}$.

**PMF derived as a limit.** Fix $k$ and $\lambda > 0$. From (3.20) with $p = \lambda/n$,
$$P(X_n = k) = \binom{n}{k}\left(\frac{\lambda}{n}\right)^k\left(1 - \frac{\lambda}{n}\right)^{n-k} = \underbrace{\frac{n(n-1)\cdots(n-k+1)}{n^k}}_{\to\ 1}\cdot\frac{\lambda^k}{k!}\cdot\underbrace{\left(1 - \frac{\lambda}{n}\right)^n}_{\to\ e^{-\lambda}}\cdot\underbrace{\left(1 - \frac{\lambda}{n}\right)^{-k}}_{\to\ 1}.$$
The first factor is a product of $k$ terms $(n - j)/n \to 1$. The third is the definition of $e^{-\lambda}$. The fourth has a fixed number of factors each tending to 1. Hence
$$P(X = k) = \frac{e^{-\lambda}\lambda^k}{k!}, \qquad k = 0, 1, 2, \dots \tag{3.24}$$
It sums to 1 because $\sum_k \lambda^k/k! = e^\lambda$. We write $X \sim \mathrm{Pois}(\lambda)$. Example 4.6 checks the limit numerically for $n = 1000$, $p = 0.003$.

**Mean and variance.**
$$\mathbb{E}[X] = \sum_{k\ge0} k\,\frac{e^{-\lambda}\lambda^k}{k!} = e^{-\lambda}\lambda\sum_{k \ge 1}\frac{\lambda^{k-1}}{(k-1)!} = e^{-\lambda}\lambda e^{\lambda} = \lambda.$$
For the variance, compute the *factorial moment* $\mathbb{E}[X(X-1)]$ instead of $\mathbb{E}[X^2]$, because $k(k-1)/k! = 1/(k-2)!$ cancels cleanly:
$$\mathbb{E}[X(X-1)] = \sum_{k \ge 2} k(k-1)\frac{e^{-\lambda}\lambda^k}{k!} = e^{-\lambda}\lambda^2\sum_{k\ge2}\frac{\lambda^{k-2}}{(k-2)!} = \lambda^2.$$
Then $\mathbb{E}[X^2] = \mathbb{E}[X(X-1)] + \mathbb{E}[X] = \lambda^2 + \lambda$ and $\operatorname{Var}(X) = \lambda^2 + \lambda - \lambda^2 = \lambda$. Mean equals variance: a diagnostic. Count data with variance much larger than the mean is *overdispersed* and not Poisson (ST-4).

**The Poisson paradigm.** If $A_1, \dots, A_n$ are events with small probabilities $p_j$, and they are independent or only weakly dependent, then the number that occur is approximately $\mathrm{Pois}(\lambda)$ with $\lambda = \sum_j p_j$. The derivation above is the independent, equal-$p_j$ case; F5 sketches the general one. Applied to birthdays: $\binom{k}{2}$ pair-events each with probability $1/365$, weakly dependent, so the number of matching pairs is approximately $\mathrm{Pois}\big(\binom{k}{2}/365\big)$ and $P(\text{no match}) \approx e^{-\binom{k}{2}/365}$, which is (3.3) again, now with a reason.

#### Discrete Uniform

**Story.** $X$ is equally likely to be any of $1, 2, \dots, n$. One fair die is $n = 6$.

**PMF.** $P(X = k) = 1/n$.

**Mean and variance.** $\mathbb{E}[X] = \frac1n\sum_{k=1}^n k = \frac{n+1}{2}$. By LOTUS and $\sum_{k=1}^n k^2 = \frac{n(n+1)(2n+1)}{6}$, $\mathbb{E}[X^2] = \frac{(n+1)(2n+1)}{6}$, so
$$\operatorname{Var}(X) = \frac{(n+1)(2n+1)}{6} - \frac{(n+1)^2}{4} = \frac{(n+1)\big[2(2n+1) - 3(n+1)\big]}{12} = \frac{(n+1)(n-1)}{12} = \frac{n^2 - 1}{12}.$$
For a die: mean $3.5$, variance $35/12 \approx 2.917$, standard deviation $1.708$.

#### Summary

| Distribution | Story | PMF $P(X = k)$ | $\mathbb{E}[X]$ | $\operatorname{Var}(X)$ |
|---|---|---|---|---|
| $\mathrm{Bern}(p)$ | one trial | $p^k q^{1-k}$, $k \in \{0, 1\}$ | $p$ | $pq$ |
| $\mathrm{Bin}(n, p)$ | successes in $n$ trials | $\binom{n}{k}p^k q^{n-k}$ | $np$ | $npq$ |
| $\mathrm{Geom}(p)$ | failures before 1st success | $q^k p$ | $q/p$ | $q/p^2$ |
| $\mathrm{NBin}(r, p)$ | failures before $r$-th success | $\binom{k+r-1}{r-1}p^r q^k$ | $rq/p$ | $rq/p^2$ |
| $\mathrm{HGeom}(w, b, n)$ | white in $n$ draws, no replacement | $\binom{w}{k}\binom{b}{n-k}/\binom{w+b}{n}$ | $nw/N$ | $npq\frac{N-n}{N-1}$ |
| $\mathrm{Pois}(\lambda)$ | count of rare events | $e^{-\lambda}\lambda^k/k!$ | $\lambda$ | $\lambda$ |
| $\mathrm{DUnif}\{1..n\}$ | fair $n$-sided die | $1/n$ | $(n+1)/2$ | $(n^2-1)/12$ |

### 3.11 Continuous random variables

**Definition 3.13.** $X$ is **continuous** if there is a function $f_X \ge 0$, the **probability density function (PDF)**, with
$$F_X(x) = P(X \le x) = \int_{-\infty}^x f_X(t)\,dt, \qquad\text{so that}\qquad f_X = F_X' \text{ where } F_X \text{ is differentiable.} \tag{3.25}$$

Consequences. $\int_{-\infty}^\infty f_X = 1$ (that is $F(\infty) = 1$). $P(a < X \le b) = F(b) - F(a) = \int_a^b f_X$. And $P(X = x) = \int_x^x f_X = 0$ for every single $x$: a continuous variable hits any exact value with probability zero, so $P(X \le x) = P(X < x)$ and endpoints never matter.

**Density is not probability.** $f_X(x)$ is not $P(X = x)$; that is 0. The probability is the *area*: $P(x < X < x + \epsilon) \approx f_X(x)\,\epsilon$. A density can exceed 1 ($X \sim \mathrm{Unif}(0, 1/2)$ has $f = 2$ on its support). A density has units of $1/[X]$: if $X$ is a time in seconds, $f_X$ is "per second," and only $f_X(x)\,dx$ is dimensionless. Every time you write a density, ask what its units are. This catches most errors in change-of-variables (F5).

**Expectation and LOTUS.** $\mathbb{E}[X] = \int x f_X(x)\,dx$ and $\mathbb{E}[g(X)] = \int g(x) f_X(x)\,dx$, when the integrals converge absolutely. The heuristic: chop the line into bins of width $dx$, put mass $f_X(x)\,dx$ at $x$, and apply the discrete definitions; the Riemann sums converge to the integrals. A rigorous proof of continuous LOTUS needs the measure-theoretic definition of $\mathbb{E}$ and is the one place in this module where we point elsewhere (Billingsley §16, or F5's sketch). Theorems 3.8, 3.11–3.14 and Definition 3.12 then hold verbatim for continuous variables, with the same proofs read through LOTUS.

#### Uniform

**Story.** A point chosen "completely at random" in $[a, b]$: probability proportional to length. $X \sim \mathrm{Unif}(a, b)$.

**PDF and CDF.** $f(x) = \frac{1}{b - a}$ on $[a, b]$, else 0. $F(x) = \frac{x - a}{b - a}$ on $[a, b]$.

**Mean and variance.** $\mathbb{E}[X] = \int_a^b \frac{x}{b-a}dx = \frac{b^2 - a^2}{2(b-a)} = \frac{a + b}{2}$. $\mathbb{E}[X^2] = \frac{b^3 - a^3}{3(b-a)} = \frac{a^2 + ab + b^2}{3}$. So
$$\operatorname{Var}(X) = \frac{a^2 + ab + b^2}{3} - \frac{(a+b)^2}{4} = \frac{4a^2 + 4ab + 4b^2 - 3a^2 - 6ab - 3b^2}{12} = \frac{(b - a)^2}{12}.$$
For $\mathrm{Unif}(0, 1)$: mean $1/2$, variance $1/12$. (The same $1/12$ as a uniform rod's moment of inertia about its center per unit mass per unit length squared. It is the same integral.)

**Theorem 3.15 (Universality of the Uniform).** Let $F$ be a CDF that is continuous and strictly increasing on the support of its distribution, with inverse $F^{-1} : (0, 1) \to \mathbb{R}$.

(a) If $U \sim \mathrm{Unif}(0, 1)$ then $X = F^{-1}(U)$ has CDF $F$.

(b) If $X$ has CDF $F$ then $F(X) \sim \mathrm{Unif}(0, 1)$.

▶ **Try it first.** For (a), compute $P(F^{-1}(U) \le x)$ by applying $F$ to both sides of the inequality. You will need $P(U \le u) = u$.

*Proof.* (a) $F$ is increasing, so $F^{-1}(U) \le x$ iff $U \le F(x)$. Hence $P(X \le x) = P(U \le F(x)) = F(x)$, since $P(U \le u) = u$ for $u \in [0, 1]$ and $F(x) \in [0, 1]$. ∎

(b) For $u \in (0, 1)$: $F(X) \le u$ iff $X \le F^{-1}(u)$ (again because $F$ is increasing). So $P(F(X) \le u) = P(X \le F^{-1}(u)) = F(F^{-1}(u)) = u$, which is the $\mathrm{Unif}(0, 1)$ CDF. ∎

■ **What just happened.** (a) is **inverse-CDF sampling**: a uniform random number generator plus the inverse CDF of any distribution gives you samples from that distribution. Example 4.7 does the Exponential; the project does a custom piecewise CDF. (b) says any continuous random variable can be turned into a uniform by pushing it through its own CDF; ST-2 uses it (a $p$-value under the null is $\mathrm{Unif}(0,1)$), and it is why Q–Q plots work. Where $F$ has flat pieces or jumps, replace $F^{-1}$ by the quantile function $F^{-1}(u) = \min\{x : F(x) \ge u\}$ and (a) still holds; (b) needs continuity.

#### Exponential

**Story.** The waiting time until an event that happens "at a constant rate $\lambda$ with no memory": the time to the next radioactive decay, the next customer, the failure of a component with no wear-out. $X \sim \mathrm{Expo}(\lambda)$, $\lambda > 0$, in units of events per unit time.

**PDF, CDF, survival.** $f(x) = \lambda e^{-\lambda x}$ for $x > 0$. $F(x) = \int_0^x \lambda e^{-\lambda t}dt = 1 - e^{-\lambda x}$. The **survival function** is $P(X > x) = e^{-\lambda x}$. $\tag{3.26}$

**Mean and variance.** Let $Y \sim \mathrm{Expo}(1)$. Integrating by parts with $u = y$, $dv = e^{-y}dy$: $\mathbb{E}[Y] = \int_0^\infty y e^{-y}dy = [-ye^{-y}]_0^\infty + \int_0^\infty e^{-y}dy = 0 + 1 = 1$. Again with $u = y^2$: $\mathbb{E}[Y^2] = [-y^2 e^{-y}]_0^\infty + 2\int_0^\infty y e^{-y}dy = 2$. So $\operatorname{Var}(Y) = 2 - 1 = 1$. Now $X = Y/\lambda$ (check: $P(Y/\lambda > x) = P(Y > \lambda x) = e^{-\lambda x}$), so by Theorems 3.8 and 3.12,
$$\mathbb{E}[X] = \frac1\lambda, \qquad \operatorname{Var}(X) = \frac{1}{\lambda^2}, \qquad \sigma_X = \frac1\lambda.$$
Mean equals standard deviation: another diagnostic. Median: $1 - e^{-\lambda m} = 1/2$ gives $m = \ln 2/\lambda \approx 0.693/\lambda$, well below the mean because of the long right tail.

**Theorem 3.16 (Memorylessness).** If $X \sim \mathrm{Expo}(\lambda)$ then for all $s, t \ge 0$,
$$P(X > s + t \mid X > s) = P(X > t).$$

*Proof.* $P(X > s + t \mid X > s) = \dfrac{P(X > s + t,\ X > s)}{P(X > s)} = \dfrac{P(X > s+t)}{P(X > s)} = \dfrac{e^{-\lambda(s+t)}}{e^{-\lambda s}} = e^{-\lambda t}$. ∎

A component that has lasted 3000 hours is, as far as its remaining life goes, as good as new. This is a strong modeling assumption, correct for radioactive decay and wrong for anything that wears out; Example 4.10 shows what it implies numerically.

**Theorem 3.17 (Memorylessness characterizes the Exponential).** Let $X$ be a positive random variable with a continuous CDF, with $P(X > t) \to 0$ as $t \to \infty$, satisfying $P(X > s + t \mid X > s) = P(X > t)$ for all $s, t \ge 0$. Then $X \sim \mathrm{Expo}(\lambda)$ for some $\lambda > 0$.

▶ **Try it first.** Write $G(t) = P(X > t)$. The hypothesis says $G(s + t) = G(s)G(t)$. What functions satisfy that? Start with $G(2) = G(1)^2$, then $G(1/2)$, then $G(m/n)$.

**Derivation.**

1. Let $G(t) = P(X > t)$, so $G(0) = 1$, $G$ is continuous (since $F$ is), non-increasing, and $G(t) \to 0$. The hypothesis is $G(s + t) = G(s)\,G(t)$ for all $s, t \ge 0$.
2. *$G$ never vanishes.* If $G(t_0) = 0$ for some $t_0 > 0$, then $G(t_0/2)^2 = G(t_0) = 0$, so $G(t_0/2) = 0$, and by induction $G(t_0/2^n) = 0$ for all $n$. But $t_0/2^n \to 0$ and $G$ is continuous with $G(0) = 1$. Contradiction.
3. *Set $\lambda$.* $0 < G(1) < 1$: positive by step 2, and if $G(1) = 1$ then $G(n) = G(1)^n = 1$ for all integers $n$, contradicting $G \to 0$. Define $\lambda = -\ln G(1) > 0$, so $G(1) = e^{-\lambda}$.
4. *Rationals.* $G(1) = G(\tfrac1n + \cdots + \tfrac1n) = G(\tfrac1n)^n$, so $G(\tfrac1n) = e^{-\lambda/n}$ (positive root, by step 2). Then $G(\tfrac mn) = G(\tfrac1n)^m = e^{-\lambda m/n}$.
5. *Reals.* For irrational $t$, choose rationals $r_k \to t$. Continuity gives $G(t) = \lim G(r_k) = \lim e^{-\lambda r_k} = e^{-\lambda t}$.
6. So $P(X > t) = e^{-\lambda t}$ for all $t \ge 0$, which is (3.26). ∎

■ **What just happened.** The functional equation $G(s + t) = G(s)G(t)$ forces an exponential, exactly as $h(s + t) = h(s) + h(t)$ with continuity forces a linear function. Among continuous positive distributions, "no memory" and "exponential" are the same thing. Exercise C3 does the discrete version: no memory on $\{0, 1, 2, \dots\}$ forces the Geometric, which is why the Geometric is the discrete shadow of the Exponential.

**Exponential and the Poisson process (stated; F5 proves it).** Suppose events occur in time so that the number in any interval of length $t$ is $\mathrm{Pois}(\lambda t)$ and counts in disjoint intervals are independent. Then the waiting time $T_1$ to the first event satisfies $P(T_1 > t) = P(\text{no events in } [0, t]) = e^{-\lambda t}$, so $T_1 \sim \mathrm{Expo}(\lambda)$; and the gaps between successive events are i.i.d. $\mathrm{Expo}(\lambda)$. Conversely, i.i.d. exponential gaps produce Poisson counts. The rate $\lambda$ is the same number in both descriptions: events per unit time.

### 3.12 The Normal distribution

**Definition 3.14.** The **standard Normal** $Z \sim \mathcal{N}(0, 1)$ has PDF
$$\varphi(z) = \frac{1}{\sqrt{2\pi}}\,e^{-z^2/2}, \qquad z \in \mathbb{R}, \tag{3.27}$$
and CDF $\Phi(z) = \int_{-\infty}^z \varphi$.

Where does $1/\sqrt{2\pi}$ come from? It is whatever constant makes the area 1, and finding it requires the one genuine trick in this module.

**Theorem 3.18 (The Gaussian integral).** $\displaystyle\int_{-\infty}^\infty e^{-z^2/2}\,dz = \sqrt{2\pi}$.

▶ **Try it first.** You cannot do this integral directly ($e^{-z^2/2}$ has no elementary antiderivative). Instead call the integral $I$, write $I^2$ as a double integral over the plane, and change to polar coordinates. Try it before reading; then memorize the move, because you will not rediscover it under pressure.

**Derivation.**

1. Let $I = \int_{-\infty}^\infty e^{-x^2/2}dx$. Then $I^2 = \left(\int e^{-x^2/2}dx\right)\left(\int e^{-y^2/2}dy\right) = \iint_{\mathbb{R}^2} e^{-(x^2 + y^2)/2}\,dx\,dy$, because the integrand factors and the two integrals have independent variables.
2. Polar coordinates: $x = r\cos\theta$, $y = r\sin\theta$, $x^2 + y^2 = r^2$, $dx\,dy = r\,dr\,d\theta$. So
$$I^2 = \int_0^{2\pi}\!\int_0^\infty e^{-r^2/2}\,r\,dr\,d\theta.$$
3. The extra factor $r$ is what makes the inner integral elementary: with $u = r^2/2$, $du = r\,dr$, $\int_0^\infty e^{-r^2/2}r\,dr = \int_0^\infty e^{-u}du = 1$.
4. $I^2 = \int_0^{2\pi} 1\,d\theta = 2\pi$, and $I > 0$, so $I = \sqrt{2\pi}$. ∎

■ **What just happened.** One copy of the integral is impossible; the product of two copies is a rotationally symmetric integral over the plane, and rotational symmetry hands you the Jacobian factor $r$ that makes it trivial. Example 4.8 checks it numerically. The trick reappears when you compute $\Gamma(1/2)$ in §3.13 and when F5 normalizes the multivariate Gaussian.

**Mean and variance of $Z$.** $\mathbb{E}[Z] = \frac{1}{\sqrt{2\pi}}\int z e^{-z^2/2}dz = 0$: the integrand is odd and the integral converges absolutely ($\int |z|e^{-z^2/2}dz = 2$, by the substitution $u = z^2/2$ on each half). For the variance, $\operatorname{Var}(Z) = \mathbb{E}[Z^2] - 0$, and by LOTUS
$$\mathbb{E}[Z^2] = \frac{1}{\sqrt{2\pi}}\int_{-\infty}^\infty z\cdot\big(z e^{-z^2/2}\big)dz.$$
Integrate by parts with $u = z$ and $dv = ze^{-z^2/2}dz$, so $v = -e^{-z^2/2}$:
$$\mathbb{E}[Z^2] = \frac{1}{\sqrt{2\pi}}\Big(\big[-ze^{-z^2/2}\big]_{-\infty}^{\infty} + \int_{-\infty}^\infty e^{-z^2/2}dz\Big) = \frac{1}{\sqrt{2\pi}}\big(0 + \sqrt{2\pi}\big) = 1.$$
The boundary term vanishes because $e^{-z^2/2}$ beats any polynomial. So $\mathbb{E}[Z] = 0$, $\operatorname{Var}(Z) = 1$: the name is honest.

**The general Normal by location-scale.** Define $X = \mu + \sigma Z$ with $\sigma > 0$. By Theorems 3.8 and 3.12, $\mathbb{E}[X] = \mu$ and $\operatorname{Var}(X) = \sigma^2$. Its CDF is $F_X(x) = P(\mu + \sigma Z \le x) = P\!\left(Z \le \frac{x - \mu}{\sigma}\right) = \Phi\!\left(\frac{x-\mu}{\sigma}\right)$, and differentiating with the chain rule,
$$f_X(x) = \frac{1}{\sigma}\,\varphi\!\left(\frac{x - \mu}{\sigma}\right) = \frac{1}{\sigma\sqrt{2\pi}}\exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right). \tag{3.28}$$
We write $X \sim \mathcal{N}(\mu, \sigma^2)$. The $1/\sigma$ in front is the units: $f_X$ must have units of $1/[X]$, and $\varphi$ is dimensionless.

The same two moves (write the CDF, then differentiate) give the general one-variable **change of variables**. If $Y = g(X)$ with $g$ strictly increasing and differentiable, then $F_Y(y) = P(g(X) \le y) = F_X(g^{-1}(y))$, and the chain rule gives $f_Y(y) = f_X(g^{-1}(y))\,\frac{d}{dy}g^{-1}(y)$. For $g$ strictly decreasing the inequality flips, $F_Y(y) = 1 - F_X(g^{-1}(y))$, and the derivative picks up a minus sign that cancels the negative $\frac{d}{dy}g^{-1}$. Both cases together:
$$f_Y(y) = f_X\big(g^{-1}(y)\big)\left|\frac{d}{dy}g^{-1}(y)\right|. \tag{3.28a}$$
The absolute derivative is the factor by which $g^{-1}$ stretches a small interval, which is exactly what keeps $f_Y(y)\,dy = f_X(x)\,dx$. F5 replaces it by a Jacobian determinant in several variables.

**Standardization.** Reverse the construction: if $X \sim \mathcal{N}(\mu, \sigma^2)$ then $Z = (X - \mu)/\sigma \sim \mathcal{N}(0, 1)$. Every Normal probability reduces to $\Phi$: $P(a < X \le b) = \Phi\!\left(\frac{b-\mu}{\sigma}\right) - \Phi\!\left(\frac{a-\mu}{\sigma}\right)$.

**The 68–95–99.7 rule.** $P(|Z| \le 1) = \Phi(1) - \Phi(-1) = 2\Phi(1) - 1 = 0.6827$; $P(|Z| \le 2) = 0.9545$; $P(|Z| \le 3) = 0.9973$. For any Normal: 68% of the mass within one $\sigma$ of the mean, 95% within two, 99.7% within three. (Here $\Phi(-z) = 1 - \Phi(z)$ by the symmetry of $\varphi$.) The project checks these from $\Phi$ and from samples.

**$\Phi$ has no closed form.** Liouville proved (1835) that $e^{-z^2}$ has no elementary antiderivative, so $\Phi$ is a new special function. It is tabulated and computed numerically (`scipy.stats.norm.cdf`), and written in terms of the error function as $\Phi(z) = \tfrac12\big(1 + \operatorname{erf}(z/\sqrt2)\big)$. You will never "solve for $x$" in $\Phi(x) = 0.99$ by algebra; you invert numerically ($z_{0.99} = 2.326$).

**Why the Normal is everywhere (preview of F5).** The **central limit theorem**: if $X_1, X_2, \dots$ are i.i.d. with mean $\mu$ and finite variance $\sigma^2$, then $\dfrac{X_1 + \cdots + X_n - n\mu}{\sigma\sqrt n} \to \mathcal{N}(0, 1)$ in distribution as $n \to \infty$. Whatever the $X_i$ look like (dice, Bernoullis, exponentials), their standardized sum looks Normal. De Moivre saw the Binomial case in 1733. F5 proves the general case with the MGFs of §3.14: the MGF of the standardized sum tends to $e^{t^2/2}$, and uniqueness finishes it. Measurement errors, heights, and the thermal velocity of a gas molecule (Maxwell 1860) are all sums of many small independent contributions.

### 3.13 Gamma and Beta (densities and normalizing constants only; F5 does the rest)

**Definition 3.15 (The Gamma function).** For $a > 0$,
$$\Gamma(a) = \int_0^\infty x^{a-1}e^{-x}\,dx. \tag{3.29}$$
The integral converges: near 0 the integrand is $\sim x^{a-1}$, integrable for $a > 0$; at infinity $e^{-x}$ wins.

**Theorem 3.19.** (a) $\Gamma(a + 1) = a\,\Gamma(a)$. (b) $\Gamma(n) = (n-1)!$ for integers $n \ge 1$. (c) $\Gamma(\tfrac12) = \sqrt\pi$.

*Proof.* (a) Integrate by parts with $u = x^a$, $dv = e^{-x}dx$: $\Gamma(a+1) = \int_0^\infty x^a e^{-x}dx = [-x^a e^{-x}]_0^\infty + a\int_0^\infty x^{a-1}e^{-x}dx = 0 + a\Gamma(a)$. (b) $\Gamma(1) = \int_0^\infty e^{-x}dx = 1 = 0!$. If $\Gamma(n) = (n-1)!$ then by (a) $\Gamma(n+1) = n\Gamma(n) = n\cdot(n-1)! = n!$. Induction. (c) $\Gamma(\tfrac12) = \int_0^\infty x^{-1/2}e^{-x}dx$. Substitute $x = u^2/2$, $dx = u\,du$, $x^{-1/2} = \sqrt2/u$:
$$\Gamma(\tfrac12) = \int_0^\infty \frac{\sqrt2}{u}\,e^{-u^2/2}\,u\,du = \sqrt2\int_0^\infty e^{-u^2/2}du = \sqrt2\cdot\frac{\sqrt{2\pi}}{2} = \sqrt\pi,$$
using Theorem 3.18 and the evenness of $e^{-u^2/2}$. ∎

So $\Gamma$ interpolates the factorial, and $\Gamma(\tfrac12) = \sqrt\pi$ is the Gaussian integral in disguise. Also $\Gamma(\tfrac32) = \tfrac12\sqrt\pi$, $\Gamma(\tfrac52) = \tfrac34\sqrt\pi$, by (a).

**Gamma distribution.** $X \sim \mathrm{Gamma}(a, \lambda)$, $a, \lambda > 0$, has PDF
$$f(x) = \frac{\lambda^a}{\Gamma(a)}\,x^{a-1}e^{-\lambda x}, \qquad x > 0.$$
It integrates to 1: substitute $u = \lambda x$ to get $\frac{\lambda^a}{\Gamma(a)}\int_0^\infty (u/\lambda)^{a-1}e^{-u}\frac{du}{\lambda} = \frac{1}{\Gamma(a)}\Gamma(a) = 1$. $\mathrm{Gamma}(1, \lambda) = \mathrm{Expo}(\lambda)$. The story, proved in F5: $\mathrm{Gamma}(n, \lambda)$ is the sum of $n$ i.i.d. $\mathrm{Expo}(\lambda)$, i.e., the time of the $n$-th event in a Poisson process; and $\mathrm{Gamma}(n/2, 1/2)$ is the $\chi^2_n$ distribution of ST-2. The mean is $a/\lambda$ and the variance $a/\lambda^2$ (Exercise B6 derives them from (3.29) and Theorem 3.19(a)).

**Beta distribution.** $X \sim \mathrm{Beta}(a, b)$, $a, b > 0$, has PDF
$$f(x) = \frac{1}{B(a, b)}\,x^{a-1}(1-x)^{b-1}, \qquad 0 < x < 1, \qquad B(a, b) = \int_0^1 x^{a-1}(1-x)^{b-1}dx = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a + b)}.$$
The last equality is stated here and proved in F5 (it is a two-dimensional change of variables on the product of two Gamma integrals, the same move as Theorem 3.18). $\mathrm{Beta}(1, 1) = \mathrm{Unif}(0, 1)$. The Beta is the distribution of a probability: in ST-5 it is the conjugate prior for the $p$ of a Binomial, and in F5 it is the distribution of the $k$-th smallest of $n$ uniforms.

### 3.14 Moment generating functions

**Definition 3.16.** The **moment generating function (MGF)** of $X$ is
$$M_X(t) = \mathbb{E}\big[e^{tX}\big], \tag{3.30}$$
provided this is finite for all $t$ in some open interval $(-\epsilon, \epsilon)$ around 0. Not every distribution has one (the Cauchy of the project does not; $e^{tX}$ has infinite expectation for every $t \ne 0$).

**Theorem 3.20 (MGFs generate moments).** If $M_X$ exists then $\mathbb{E}[X^n] = M_X^{(n)}(0)$, the $n$-th derivative at 0, and
$$M_X(t) = \sum_{n=0}^\infty \mathbb{E}[X^n]\,\frac{t^n}{n!}.$$

▶ **Try it first.** Expand $e^{tX}$ as a power series in $t$ and take the expectation of each term.

**Derivation.**

1. $e^{tX} = \sum_{n\ge0} \frac{t^n X^n}{n!}$ for every outcome.
2. Take expectations and exchange $\mathbb{E}$ with the sum: $M_X(t) = \sum_n \frac{t^n}{n!}\mathbb{E}[X^n]$. (The exchange is legitimate when $M_X$ is finite near 0; the series then converges absolutely and the justification is the dominated convergence theorem, which is the second and last place this module points to measure theory.)
3. A power series $\sum_n c_n t^n$ convergent near 0 has $c_n = M^{(n)}(0)/n!$ (Taylor). Matching coefficients, $\mathbb{E}[X^n]/n! = M_X^{(n)}(0)/n!$. ∎

■ **What just happened.** All the moments of $X$ are packed into one function, and differentiation unpacks them. In practice the Taylor-coefficient route is often faster than differentiating: Exercise B8 reads $\mathbb{E}[Z^4] = 3$ straight off the series of $e^{t^2/2}$. The logarithm $K_X(t) = \ln M_X(t)$ is the *cumulant generating function*; $K'(0) = \mathbb{E}[X]$ and $K''(0) = \operatorname{Var}(X)$ (check: $K' = M'/M$, $K'' = (M''M - M'^2)/M^2$, evaluate at $t = 0$ where $M = 1$). In SM-4 the partition function is $\ln Z(\beta) = K_E(-\beta)$ for the energy, and $\langle E\rangle = -\partial_\beta \ln Z$, $\operatorname{Var}(E) = \partial^2_\beta \ln Z$ are exactly these two identities.

**Theorem 3.21 (Five MGFs).**

| $X$ | $M_X(t)$ | valid for |
|---|---|---|
| $\mathrm{Bern}(p)$ | $q + pe^t$ | all $t$ |
| $\mathrm{Bin}(n, p)$ | $(q + pe^t)^n$ | all $t$ |
| $\mathrm{Pois}(\lambda)$ | $\exp\big(\lambda(e^t - 1)\big)$ | all $t$ |
| $\mathcal{N}(\mu, \sigma^2)$ | $\exp\big(\mu t + \tfrac12\sigma^2 t^2\big)$ | all $t$ |
| $\mathrm{Expo}(\lambda)$ | $\dfrac{\lambda}{\lambda - t}$ | $t < \lambda$ |

*Derivations.* *Bernoulli:* LOTUS, $\mathbb{E}[e^{tX}] = e^{0}q + e^{t}p$. *Binomial:* by LOTUS and the binomial theorem, $\sum_k \binom nk p^k q^{n-k}e^{tk} = \sum_k \binom nk (pe^t)^k q^{n-k} = (pe^t + q)^n$. (Or: Theorem 3.22 below applied to $n$ independent Bernoullis.) *Poisson:* $\sum_k e^{tk}\frac{e^{-\lambda}\lambda^k}{k!} = e^{-\lambda}\sum_k \frac{(\lambda e^t)^k}{k!} = e^{-\lambda}e^{\lambda e^t}$. *Standard Normal:* complete the square in the exponent, $tz - z^2/2 = -\tfrac12(z - t)^2 + \tfrac12 t^2$, so
$$M_Z(t) = \frac{1}{\sqrt{2\pi}}\int e^{tz - z^2/2}dz = e^{t^2/2}\cdot\frac{1}{\sqrt{2\pi}}\int e^{-(z-t)^2/2}dz = e^{t^2/2},$$
the last integral being the total mass of a Normal centered at $t$. *General Normal:* $X = \mu + \sigma Z$, $M_X(t) = \mathbb{E}[e^{t\mu}e^{t\sigma Z}] = e^{\mu t}M_Z(\sigma t) = e^{\mu t + \sigma^2t^2/2}$. *Exponential:* $\int_0^\infty e^{tx}\lambda e^{-\lambda x}dx = \lambda\int_0^\infty e^{-(\lambda - t)x}dx = \frac{\lambda}{\lambda - t}$, which converges only for $t < \lambda$. ∎

Check against §3.10–3.12: for Poisson, $M'(t) = \lambda e^t M(t)$ so $M'(0) = \lambda$; $M''(t) = (\lambda e^t)^2 M(t) + \lambda e^t M(t)$ so $M''(0) = \lambda^2 + \lambda = \mathbb{E}[X^2]$. For the Exponential, $M(t) = (1 - t/\lambda)^{-1} = \sum_n (t/\lambda)^n$, so $\mathbb{E}[X^n] = n!/\lambda^n$: mean $1/\lambda$, second moment $2/\lambda^2$. For $\mathcal{N}(0,1)$, $e^{t^2/2} = \sum_k \frac{t^{2k}}{2^k k!}$: odd moments vanish and $\mathbb{E}[Z^{2k}] = \frac{(2k)!}{2^k k!}$, so $\mathbb{E}[Z^2] = 1$, $\mathbb{E}[Z^4] = 3$.

**Theorem 3.22 (MGF of a sum of independent variables).** If $X \perp Y$ then $M_{X+Y}(t) = M_X(t)\,M_Y(t)$.

*Proof.* $M_{X+Y}(t) = \mathbb{E}[e^{t(X+Y)}] = \mathbb{E}[e^{tX}e^{tY}] = \mathbb{E}[e^{tX}]\,\mathbb{E}[e^{tY}]$, the last step by Theorem 3.13 with $g(x) = e^{tx}$, $h(y) = e^{ty}$. ∎

**Theorem 3.23 (Uniqueness; stated).** If $M_X(t) = M_Y(t)$ for all $t$ in some open interval around 0, then $X$ and $Y$ have the same distribution.

The proof inverts a Laplace transform and is in Billingsley §30 or Feller vol. II; we use the theorem as a tool. Together, Theorems 3.22 and 3.23 turn "what is the distribution of $X + Y$?" into "multiply two functions and recognize the product."

**Application 1: Binomials add.** If $X \sim \mathrm{Bin}(n, p)$ and $Y \sim \mathrm{Bin}(m, p)$ are independent, $M_{X+Y}(t) = (q + pe^t)^n(q + pe^t)^m = (q + pe^t)^{n+m}$, which is the MGF of $\mathrm{Bin}(n + m, p)$. By uniqueness, $X + Y \sim \mathrm{Bin}(n+m, p)$. ∎ (Same for Poissons: $e^{\lambda(e^t-1)}e^{\mu(e^t-1)} = e^{(\lambda+\mu)(e^t-1)}$, so independent Poissons add to a Poisson. And for Normals, $e^{\mu_1 t + \sigma_1^2t^2/2}e^{\mu_2 t + \sigma_2^2 t^2/2}$ is the MGF of $\mathcal{N}(\mu_1 + \mu_2, \sigma_1^2 + \sigma_2^2)$: independent Normals add to a Normal, with variances adding, as Theorem 3.14 said they must.)

**Application 2: The Poisson limit, again.** With $p = \lambda/n$,
$$M_{\mathrm{Bin}(n, \lambda/n)}(t) = \left(1 - \frac\lambda n + \frac\lambda n e^t\right)^n = \left(1 + \frac{\lambda(e^t - 1)}{n}\right)^n \;\xrightarrow[n\to\infty]{}\; \exp\big(\lambda(e^t - 1)\big) = M_{\mathrm{Pois}(\lambda)}(t),$$
using $(1 + c/n)^n \to e^c$. The MGFs converge to the Poisson MGF, which (by a continuity theorem that F5 states alongside the CLT) means the distributions converge: the same conclusion as §3.10, in one line instead of four factors. This is the template for F5's proof of the CLT.

### 3.15 How to solve a probability problem

1. **Name things in symbols before computing anything.** "Let $D$ be the event the part is defective; let $X$ be the number of defective parts in the sample." Half of all errors in this subject are ambiguities in the English that symbols would have exposed (Trap 3).
2. **Write down exactly what is asked** as $P(\cdot)$, $P(\cdot \mid \cdot)$, or $\mathbb{E}[\cdot]$.
3. **Condition on what you wish you knew.** If knowing $B$ would make the problem easy, use LOTP (3.11) over a partition, or first-step analysis. If you know $B$ and want $A$ but the problem gives you $P(B \mid A)$, use Bayes (3.10).
4. **Look for symmetry and indicators.** "Expected number of ..." is almost always $\sum_j P(A_j)$. "Probability that the $j$-th draw is ..." is almost always the same as for the first draw.
5. **Recognize a named distribution from its story**, then look up the mean and variance instead of recomputing.
6. **Sanity-check.** Probabilities lie in $[0, 1]$. Take $p \to 0$, $p \to 1$, $n \to \infty$ and see if the answer does something sensible. Check units of any density. Check $\sum p(x) = 1$.
7. **Simulate.** Ten lines of NumPy with $10^5$ trials settle most disputes to two decimals. The project is a template. If the simulation and the formula disagree, the formula is usually wrong, but check that the simulation implements the *same* conditioning (Monty Hall with a random host is a different simulation from Monty Hall with a knowing host).

## 4. Worked examples

Every number below is reproduced by an `assert` or a `print` in `project.py` (the part number is given in brackets).

### Example 4.1: The birthday problem for 23 people [Part 1]

*Exact.* From §3.3 step 3 with $k = 23$,
$$P(\text{no match}) = \prod_{j=0}^{22}\left(1 - \frac{j}{365}\right) = \frac{365}{365}\cdot\frac{364}{365}\cdots\frac{343}{365} = 0.49270, \qquad P(\text{match}) = 0.50730.$$
(Compute it as a running product; the project does. By hand, take logs: $\sum_{j=0}^{22}\ln(1 - j/365) = -0.7078$, $e^{-0.7078} = 0.4927$.)

*Approximate.* From (3.3), $\binom{23}{2} = 253$ pairs, so $P(\text{no match}) \approx e^{-253/365} = e^{-0.69315} = 0.50000$, $P(\text{match}) \approx 0.5000$. The approximation is low by $0.007$ because $1 - x < e^{-x}$ strictly, so the product of the $(1 - j/365)$ is smaller than the product of the exponentials and the true no-match probability is *smaller*, the true match probability *larger*. For $k = 50$: exact $0.9704$, approximate $1 - e^{-1225/365} = 0.9651$. The number of pairs grows like $k^2$, which is why 23 people suffice: it is 253 chances at a 1-in-365 event, not 23.

### Example 4.2: The medical test [Part 1]

Prevalence $P(D) = 10^{-4}$; sensitivity $P(+\mid D) = 0.99$; specificity $P(-\mid D^c) = 0.99$.

1. Denominator by LOTP: $P(+) = P(+\mid D)P(D) + P(+\mid D^c)P(D^c) = 0.99(0.0001) + 0.01(0.9999) = 0.000099 + 0.009999 = 0.010098$.
2. Bayes: $P(D \mid +) = 0.000099 / 0.010098 = 0.009804$.

Just under 1%. Now suppose a second, independent test (same accuracy) is also positive. Use extra conditioning (3.12) with the first result as $E$: the prior is now $0.0098$,
$$P(D \mid ++) = \frac{0.99 \times 0.009804}{0.99\times0.009804 + 0.01\times0.990196} = \frac{0.009706}{0.009706 + 0.009902} = 0.4950.$$
Two positives: a coin flip. Three: $0.990$. Each positive multiplies the odds $P(D)/P(D^c)$ by the likelihood ratio $0.99/0.01 = 99$; the prior odds are $1/9999$, so after one test $99/9999 \approx 1/101$, after two $\approx 0.98$, after three $\approx 97$. Odds and likelihood ratios are the right currency for sequential evidence; ST-5 uses them throughout.

### Example 4.3: Monty Hall, two ways [Part 1]

Worked in full in §3.7 (Trap 1): $P(\text{car behind door 2} \mid \text{host opens 3}) = 2/3$ by Bayes, and $P(\text{switching wins}) = 2/3$ by conditioning on the car's location. The project simulates $10^5$ games with a knowing host and asserts that the switch-win frequency is within $3\sigma$ of $2/3$, where $\sigma = \sqrt{(2/3)(1/3)/10^5} = 0.00149$ is the standard deviation of a $\mathrm{Bin}(10^5, 2/3)$ proportion, and the stay-win frequency within $3\sigma$ of $1/3$. It also simulates the random-host variant, conditions on "a goat was revealed," and gets $1/2$. Same observation, different mechanism, different answer.

### Example 4.4: The matching problem for $n = 4$ and $n = 10$ [Part 1]

From (3.8): $n = 4$: $1 - \frac12 + \frac16 - \frac1{24} = \frac{24 - 12 + 4 - 1}{24} = \frac{15}{24} = 0.625$. $n = 10$: $0.632121$, versus $1 - 1/e = 0.632121$; they agree to six decimals because the omitted tail of the series is less than $1/11! = 2.5 \times 10^{-8}$. The project shuffles $10^5$ decks of 10 and finds the fraction with at least one match within $3\sigma$ of $0.6321$, and the mean number of matches within $3\sigma$ of 1.

Check the $n = 4$ value by brute force: of the $24$ permutations of $\{1, 2, 3, 4\}$, the derangements (no fixed point) number $9$, so $P(\text{at least one}) = 15/24$. ✓

### Example 4.5: Expected fixed points via indicators [Part 1]

$\mathbb{E}[X] = \sum_{j=1}^n P(\sigma(j) = j) = n\cdot\frac1n = 1$, for $n = 1, 2, 3, \dots$. Compare with the work of doing it from the distribution: $P(X = k) = \frac{1}{k!}\sum_{i=0}^{n-k}\frac{(-1)^i}{i!}$ (choose the $k$ fixed points, derange the rest), then $\sum_k k\,P(X = k)$. The indicator argument is one line and does not need the distribution at all. Exercise C4 pushes it to $\operatorname{Var}(X) = 1$ for $n \ge 2$, and the pair (mean 1, variance 1) is the fingerprint of $\mathrm{Pois}(1)$, which is in fact the limiting distribution of $X$.

### Example 4.6: The Poisson limit, numerically [Part 3]

$n = 1000$, $p = 0.003$, $\lambda = np = 3$.

| $k$ | $\mathrm{Bin}(1000, 0.003)$ | $\mathrm{Pois}(3)$ | difference |
|---|---|---|---|
| 0 | $0.049563$ | $0.049787$ | $-0.000224$ |
| 1 | $0.149137$ | $0.149361$ | $-0.000225$ |
| 2 | $0.224154$ | $0.224042$ | $+0.000112$ |
| 3 | $0.224379$ | $0.224042$ | $+0.000337$ |
| 4 | $0.168284$ | $0.168031$ | $+0.000253$ |
| 5 | $0.100869$ | $0.100819$ | $+0.000050$ |

Hand check of $k = 0$: $(0.997)^{1000} = e^{1000\ln 0.997} = e^{1000(-0.0030045)} = e^{-3.0045} = 0.04956$, versus $e^{-3} = 0.04979$. The exponent $1000\ln(1 - 0.003) = -3 - 0.0045 - \cdots$ is the whole error: $\ln(1 - p) = -p - p^2/2 - \cdots$, so $(1-p)^n = e^{-\lambda}e^{-np^2/2}\cdots = e^{-\lambda}(1 - \lambda p/2 + \cdots)$, and the relative error at $k = 0$ is about $\lambda p/2 = 0.0045$. The project plots the total variation distance $\tfrac12\sum_k |P_{\mathrm{Bin}}(k) - P_{\mathrm{Pois}}(k)|$ against $n$ for fixed $\lambda = 3$ and finds it decreasing like $1/n$.

### Example 4.7: Inverse-CDF sampling for the Exponential [Part 4]

$F(x) = 1 - e^{-\lambda x}$ on $x > 0$. Solve $u = F(x)$: $e^{-\lambda x} = 1 - u$, so $x = -\ln(1 - u)/\lambda$. By Theorem 3.15(a), if $U \sim \mathrm{Unif}(0,1)$ then
$$X = -\frac{\ln(1 - U)}{\lambda} \sim \mathrm{Expo}(\lambda).$$
Since $1 - U \sim \mathrm{Unif}(0, 1)$ too (its CDF is $P(1 - U \le u) = P(U \ge 1 - u) = u$), $X = -\ln(U)/\lambda$ works equally well and saves a subtraction. Direct check without the theorem: $P(X > x) = P(-\ln U > \lambda x) = P(U < e^{-\lambda x}) = e^{-\lambda x}$. ✓ The project draws $10^5$ samples with $\lambda = 2$, overlays the empirical CDF on $F$, checks the sample mean against $1/\lambda = 0.5$, and runs a Kolmogorov–Smirnov test (ST-2 explains the test; here it is a sanity check, asserted at the loose level $p > 0.001$). Then the same for a custom piecewise CDF, $F(x) = x/2$ on $[0, 1]$ and $\tfrac12 + (x-1)/4$ on $[1, 3]$, whose inverse is $F^{-1}(u) = 2u$ for $u < \tfrac12$ and $1 + 4(u - \tfrac12)$ for $u \ge \tfrac12$.

### Example 4.8: The Gaussian integral, numerically [Part 5]

Theorem 3.18 says $\int_{-\infty}^\infty e^{-z^2/2}dz = \sqrt{2\pi} = 2.506628$. The project computes the one-dimensional integral with `scipy.integrate.quad` (agreement to $10^{-10}$), computes the two-dimensional integral $\iint e^{-(x^2+y^2)/2}$ in Cartesian coordinates with `dblquad` and gets $2\pi = 6.283185$, and computes the polar form $\int_0^{2\pi}\int_0^\infty e^{-r^2/2}r\,dr\,d\theta$ to confirm step 2 of the proof. It then checks $\int \varphi = 1$, $\int z^2\varphi = 1$, $\int z^4\varphi = 3$, and $\Gamma(\tfrac12) = \sqrt\pi = 1.772454$ via `scipy.special.gamma`.

### Example 4.9: Defective parts in a batch (Hypergeometric vs Binomial) [Part 2]

A batch of $N = 50$ machined parts contains $w = 5$ out of tolerance. You inspect $n = 10$ chosen at random without replacement. Let $X$ be the number of defective parts in the sample.

*Exact (Hypergeometric).* $X \sim \mathrm{HGeom}(5, 45, 10)$.
$$P(X = 0) = \frac{\binom{5}{0}\binom{45}{10}}{\binom{50}{10}} = \prod_{i=0}^{9}\frac{45 - i}{50 - i} = \frac{40\cdot39\cdot38\cdot37\cdot36}{50\cdot49\cdot48\cdot47\cdot46} = \frac{78{,}960{,}960}{254{,}251{,}200} = 0.31056,$$
after cancelling the common factors $45, 44, 43, 42, 41$. So $P(\text{at least one defect found}) = 0.6894$. $\mathbb{E}[X] = nw/N = 10\cdot5/50 = 1$. $\operatorname{Var}(X) = npq\frac{N-n}{N-1} = 10(0.1)(0.9)\frac{40}{49} = 0.7347$.

*Binomial approximation* (pretend each part is independently defective with $p = 0.1$): $P(X = 0) = 0.9^{10} = 0.3487$, $\mathbb{E}[X] = 1$, $\operatorname{Var}(X) = 0.9$. The mean is exact (linearity never cared about replacement), but the approximation overstates $P(X = 0)$ by $0.038$ and the variance by $22.5\%$: with $n/N = 20\%$ of the batch inspected, the finite-population correction $40/49 = 0.816$ is far from 1. Each defect you draw depletes the pool of defects, which pulls the count toward its mean. The project checks both PMFs against `scipy.stats` and plots them side by side.

*Decision.* If the acceptance rule is "reject the batch if any sampled part is defective," this batch (10% defective) is accepted with probability $0.31$. To get the acceptance probability of a 10%-defective batch below $0.05$ you need $\binom{50 - 5}{n}/\binom{50}{n} < 0.05$, which first happens at $n = 22$ (the ratio equals $\binom{28}{5}/\binom{50}{5} = 98{,}280/2{,}118{,}760 = 0.0464$ there, against $0.0560$ at $n = 21$; the project checks).

### Example 4.10: Time to failure of a component (Exponential) [Part 2]

A bearing's time to failure $T$ (hours) is modeled as $\mathrm{Expo}(\lambda)$ with mean time to failure $1/\lambda = 2000$ h, so $\lambda = 5 \times 10^{-4}$ per hour.

1. $P(T > 3000) = e^{-3000/2000} = e^{-1.5} = 0.2231$.
2. $P(T \le 500) = 1 - e^{-0.25} = 0.2212$. More than a fifth of bearings fail in the first quarter of their mean life.
3. Median: $\ln 2/\lambda = 0.6931 \times 2000 = 1386$ h. Half of all bearings fail before $1386$ h even though the *mean* is $2000$ h; the mean is dragged up by the long tail.
4. Memorylessness: $P(T > 6000 \mid T > 3000) = P(T > 3000) = 0.2231$. A bearing that has run 3000 h has, under this model, exactly the survival prospects of a new one. For a real bearing that is false (wear), and the Weibull distribution of ST-4 and the Gamma of F5 exist to fix it. The Exponential is correct for failures caused by external random shocks at a constant rate, not for wear-out.
5. Standard deviation $= 1/\lambda = 2000$ h, equal to the mean. If your failure data has a standard deviation much smaller than the mean, the failure mechanism is not memoryless.

The project asserts items 1–4 and simulates $10^5$ lifetimes to confirm the conditional survival in item 4 to within sampling error.

## 5. Why it matters and where it connects

- **F5** starts from this module's Definition 3.7 and builds joint PMFs and PDFs, covariance and correlation (Definition 3.12 done properly, with $|\operatorname{Cov}(X,Y)| \le \sigma_X\sigma_Y$ as Cauchy–Schwarz from LA-1 Theorem 3.11 applied to the inner product $\langle X, Y\rangle = \mathbb{E}[XY]$), conditional expectation (first-step analysis from §3.10 made general), the change-of-variables formula that gives Theorem 3.19's Beta normalizing constant, the multivariate Normal, and the laws of large numbers and the CLT. The CLT proof is Theorem 3.22 plus Theorem 3.23 plus a Taylor expansion of the cumulant generating function.
- **SM-3** is §3.2 and §3.10. The multiplicity of an Einstein solid is the stars-and-bars count $\binom{n+k-1}{k}$ of Theorem 3.3; a two-state paramagnet has $\binom{N}{k}$ microstates with $k$ spins up, a Binomial; Boltzmann's $S = k_B\ln\Omega$ is the log of a binomial coefficient, and Stirling's formula turns it into the entropy of a Bernoulli. The fact that fluctuations scale as $\sqrt{N}$ while means scale as $N$ is $\operatorname{Var}(\mathrm{Bin}) = npq$ next to $\mathbb{E}[\mathrm{Bin}] = np$.
- **SM-4** reads §3.14 backwards: the partition function $Z(\beta) = \sum_i e^{-\beta E_i}$ is the MGF of the energy at $t = -\beta$ (up to the counting measure), $\ln Z$ is a cumulant generating function, and $\langle E\rangle = -\partial_\beta\ln Z$, $\operatorname{Var}(E) = \partial_\beta^2\ln Z = k_B T^2 C_V$ are the two identities after Theorem 3.20.
- **ST-1** takes every PMF and PDF in §3.10–3.13, fixes the data, and views it as a function of the parameter: the likelihood. Maximum likelihood for $\mathrm{Bin}$, $\mathrm{Pois}$, $\mathrm{Expo}$, $\mathcal{N}$ is one derivative of the formulas here. Fisher information is the variance of a score, computed by LOTUS. Unbiasedness is a statement about $\mathbb{E}$ of an estimator, and the bias–variance decomposition of MSE is equation (3.17) with a shifted origin.
- **ST-2** uses Theorem 3.15(b): under the null hypothesis a $p$-value is $\mathrm{Unif}(0, 1)$. The Kolmogorov–Smirnov test in the project is the first hypothesis test you will have run.
- **ST-5** is §3.5 and §3.7 at scale: Bayes' theorem with the prior a density ($\mathrm{Beta}(a, b)$ on a Binomial $p$, §3.13), the two-coins example (§3.6) as the model of "data independent given the parameter," and Monty Hall's lesson that the *mechanism* producing the evidence is part of the model.
- **ML-1** defines risk as $\mathbb{E}[\text{loss}]$, computed by LOTUS, and estimates it by a sample mean, justified by the law of large numbers previewed in Part 8 of the project. **ML-3**: logistic regression is a $\mathrm{Bern}(\sigma(\mathbf{w}^\top\mathbf{x}))$ model; the cross-entropy loss is its negative log PMF. **ML-10**: the ELBO is an expectation under an approximate posterior, estimated by sampling exactly as in §3.11.
- **QM-2**: $|\psi(x)|^2$ is a PDF and must have units of $1/[x]$ (the "density is not probability" warning of §3.11 verbatim); $\langle \hat{x}\rangle = \int x|\psi|^2dx$ is LOTUS; $\Delta x$ is a standard deviation; and the Gaussian wave packet's normalization is Theorem 3.18.
- **LA-1** Theorem 3.11 (Cauchy–Schwarz) is, in F5, the reason correlation lies in $[-1, 1]$. The "introduce a free parameter $t$, use positivity, read off the discriminant" move from that proof is how F5 proves it.

## 6. Exercises

See [exercises.md](exercises.md). **A1–A5** are drills: de Méré's dice, two story proofs, a conditional probability from a table, a PMF/CDF/expectation computation, Normal and Exponential numbers. **B1–B8** are the module: B1 (inclusion–exclusion by indicators), B2 (the two-coin Bayes problem and conditional independence), B3 (pairwise vs mutual), B4 (Geometric variance and the first-success shift), B5 (Hypergeometric by indicators), B6 (Gamma moments), B7 (universality and a Cauchy sampler; extends `project.py`), B8 (Normal MGF and a sum with zero variance). B2, B4, B7 are essential. **C1–C4**: the Tuesday-boy problem, expected waits for $HH$ versus $HT$, the discrete memorylessness characterization, and the variance of the number of fixed points.

## 7. Coding project

`project.py` is a `Simulator`-style script: a seeded `numpy.random.default_rng`, exact formulas implemented from the notes, and every claim checked against simulation or against `scipy.stats`.

1. **Monte Carlo of the traps.** Birthday problem: simulated match frequency against the exact product for $k = 1, \dots, 60$, with the $e^{-k(k-1)/730}$ approximation overlaid (figure). Monty Hall: $10^5$ games, switch vs stay, asserted within $3\sigma$ of $2/3$ and $1/3$; plus the random-host variant giving $1/2$. The matching problem: mean fixed points $\to 1$, $P(\text{none}) \to 1/e$. The base-rate example: simulate a million patients and read off the positive predictive value.
2. **Exact PMFs from the formulas** (3.20)–(3.24) for Binomial, Geometric, Negative Binomial, Hypergeometric, Poisson, asserted against `scipy.stats` (watch the convention differences, commented in the code); means and variances computed from the PMFs and asserted against the table in §3.10; empirical histograms from `rng` overlaid on the PMFs (figure). Examples 4.9 and 4.10 asserted.
3. **Poisson approximation error** (total variation distance) against $n$ for fixed $\lambda = 3$ (figure), and the table of Example 4.6.
4. **Inverse-CDF sampling** for $\mathrm{Expo}(2)$ and for the piecewise CDF of Example 4.7, with empirical-CDF overlays (figure) and `scipy.stats.kstest` assertions at $p > 0.001$; universality in the other direction, $F(X) \sim \mathrm{Unif}$.
5. **The Gaussian integral** by `quad`, by `dblquad` in Cartesian coordinates, and in polar form; $\Gamma(1/2)$; the Normal moments.
6. **68–95–99.7** from $\Phi$ and from $10^6$ samples.
7. **MGF numerically.** The sample MGF $\hat M(t) = \frac1N\sum_i e^{tX_i}$ for $\mathrm{Pois}(2)$, its first three derivatives at 0 by central finite differences, compared with $\mathbb{E}[X] = 2$, $\mathbb{E}[X^2] = 6$, $\mathbb{E}[X^3] = 22$ (the exact values are obtained by differentiating $e^{\lambda(e^t - 1)}$ with SymPy). The Binomial-sum and Poisson-limit MGF identities are checked symbolically.
8. **Law of large numbers** running means for $\mathrm{Bern}(0.3)$ (converges) and standard Cauchy (does not; figure). The Cauchy has no mean, so there is nothing to converge to.
9. **Linearity on dependent variables.** $X$ uniform on $\{-1, 0, 1\}$, $Y = X^2$: $\mathbb{E}[X + Y] = \mathbb{E}X + \mathbb{E}Y$ exactly and by simulation; also $\operatorname{Cov}(X, Y) = 0$ despite dependence; and a dependent pair (heads and tails counts) where $\operatorname{Var}(X + Y) = 0 \ne \operatorname{Var}X + \operatorname{Var}Y$.

**Build it yourself first:** before reading Part 2, write `binom_pmf(k, n, p)`, `hgeom_pmf(k, w, b, n)`, and `pois_pmf(k, lam)` from (3.20), (3.23), (3.24) using `math.comb` and `math.exp`, and `expo_sample(lam, size, rng)` from Example 4.7. Then compare.

**Your turn** (prompts at the end of the script): simulate the two-children problem under phrasings (i), (ii), (iii) of Trap 3 and recover $1/3$, $1/2$, $1/2$; estimate the probability that among 23 people at least one pair has birthdays within one day of each other (circularly, so Dec 31 and Jan 1 count); derive $\mathbb{E}[\text{rolls to see all six faces}] = 6\left(1 + \tfrac12 + \tfrac13 + \tfrac14 + \tfrac15 + \tfrac16\right) = 14.7$ as a sum of independent first-success waiting times with $p = 6/6, 5/6, \dots, 1/6$, and simulate it.

## 8. Where you'll struggle

- **Conditional probability.** Your intuition for it is broken. So is everyone's; the roadmap says fifty problems and means it. The five traps in §3.7 are not curiosities, they are the standard failure modes, and you will fall for at least two of them on first contact even having read the solutions. The fix is mechanical: (1) write the sample space, (2) write the conditioning event as a *set*, (3) compute $P(A \cap B)/P(B)$. Never reason in English past step 1. When a problem says "given," ask *how* the information arrived (Monty Hall with a knowing host vs a random host; "at least one boy" vs "the elder is a boy"). If you cannot write down the mechanism, you cannot write down $P(B)$.
- **Independent vs disjoint.** You will use one word when you mean the other. Disjoint: $P(A \cap B) = 0$. Independent: $P(A \cap B) = P(A)P(B)$. Disjoint events with positive probability are maximally dependent. Say "cannot both happen" for disjoint and "doesn't change the odds" for independent until the words stop sliding.
- **Forgetting that a random variable is a function.** The symptom is writing things like "$X = 3$" as if $X$ were a number, or asking "what is $X$?" The cure is to write $X(\omega)$ explicitly for a few problems: the sum of two dice is $X(i, j) = i + j$, the number of matches is $X(\sigma) = \#\{j : \sigma(j) = j\}$. The payoff is that linearity (Theorem 3.8) becomes obvious instead of mysterious, and "$X$ and $Y$ have the same distribution but are different" stops being a riddle.
- **Treating a density as a probability.** You will write $f(x)$ where you mean $f(x)\,dx$, be surprised that $f > 1$, or forget the $1/\sigma$ in (3.28). Check units. A density has units of $1/[X]$. If your answer has the wrong units it is wrong; this catches most change-of-variable errors in F5 and most normalization errors in QM-2.
- **The Gaussian integral is not derivable from scratch.** Nobody looks at $\int e^{-z^2/2}dz$ and thinks of squaring it and going to polar coordinates; that is a trick someone found, and the only way to have it is to have memorized it. Memorize: *square it, polar, the Jacobian $r$ makes it elementary.* The same for $\Gamma(1/2)$: *substitute $x = u^2/2$ and it becomes the Gaussian.* Do not spend an hour trying to find a direct route.
- **MGFs look like busywork.** You will compute five of them, prove two things you already believed (Binomials add; Poisson is the Binomial limit), and wonder why. The answer is that F5's proof of the central limit theorem is six lines *with* MGFs and does not exist at this level without them, and SM-4's partition function is one. Learn Theorem 3.20's expansion and Theorem 3.22's product now so that F5 can spend its time on the theorem instead of the tool.
- **Counting problems are hard and stay hard.** The story proofs in §3.2 look effortless because the right story is given. Finding the right story is a skill that improves slowly, and people who are excellent at probability are still routinely stuck on counting problems. Do not judge your grasp of probability by your speed at counting. If a counting problem resists for twenty minutes, compute small cases by brute force (the project shows how), guess the pattern, and move on; the probability content of this module does not depend on it.
- **Which Geometric.** Blitzstein and these notes count failures (support from 0, mean $q/p$); Ross and SciPy's `geom` count trials (support from 1, mean $1/p$). Both are in the literature. Before using any formula or library call, write down which one it assumes. The project's Part 2 comments mark every convention conversion.
- **The exercises that are harder than they look.** B4(a) (differentiating the geometric series twice: the powers of $q$ pulled out in front are easy to get wrong by one; check your answer at $p = 1/2$, where $\operatorname{Var}(X) = 2$). B7(b) (the Cauchy has no mean, so "check the sample mean" is a trap; the exercise wants you to see that). C2 (the $HH$ vs $HT$ asymmetry: after a failed attempt at $HH$ you are back at the start, but after a failed attempt at $HT$ you are halfway).

## 9. Check yourself

1. State Kolmogorov's three axioms and derive $P(A^c) = 1 - P(A)$ from them in one line.
2. Prove $\binom{n}{k} = \binom{n}{n-k}$ with a story, not algebra.
3. Write Bayes' theorem and derive it from the definition of conditional probability in two lines.
4. Two events are disjoint and each has probability $0.3$. Are they independent? Why not?
5. Give an example of events that are conditionally independent given $E$ but not independent.
6. What is a random variable? What is a distribution? Give two different random variables with the same distribution.
7. State linearity of expectation. What assumption about $X$ and $Y$ does it need? Prove it.
8. What is $\mathbb{E}[I_A]$? Use it to find the expected number of fixed points of a random permutation of $n$ objects.
9. Derive the Poisson PMF as a limit of Binomials. Which three limits did you use?
10. Write $\operatorname{Var}(X) = \mathbb{E}[X^2] - (\mathbb{E}X)^2$ and name the mechanics theorem it corresponds to.
11. If $U \sim \mathrm{Unif}(0, 1)$, how do you produce an $\mathrm{Expo}(\lambda)$ sample? Prove it works.
12. State the memorylessness property. Which continuous distributions have it?
13. Compute $\int_{-\infty}^\infty e^{-z^2/2}dz$. Outline the trick in one sentence.
14. What is $M_X(t)$? Why is $M_{X+Y} = M_X M_Y$ for independent $X, Y$? What theorem lets you conclude the distribution of $X + Y$ from that product?
15. A density takes the value $2.5$ at some point. Is that a problem? What are its units?

## 10. Pair with

- **Blitzstein & Hwang, *Introduction to Probability*, 2nd ed. (free at probabilitybook.net), ch. 1–6.** The source of the story-proof style used here, with far more worked problems than these notes. Ch. 1 counting and story proofs, ch. 2 conditioning and the traps (their treatment of Monty Hall and of the two-coins example is the model for §3.6–3.7), ch. 3–4 discrete distributions and expectation (indicator trick, LOTUS), ch. 5 continuous distributions (universality of the Uniform, their §5.3, is better explained there than anywhere), ch. 6 moments and MGFs. Do their end-of-chapter problems; the ones marked with a circle have solutions online.
- **Harvard Stat 110 (Blitzstein), lectures 1–15 on YouTube.** Lectures 1–3 counting, story proofs, axioms, birthdays; 4–6 conditioning, Bayes, Monty Hall, Simpson; 7–8 random variables, Binomial, Hypergeometric; 9–10 expectation and indicators; 11 Poisson and the Poisson paradigm; 12–13 continuous distributions, the Uniform and universality, the Normal; 14 location–scale and LOTUS. Lectures 16 (Exponential, memorylessness) and 17 (MGFs) cover §3.11 and §3.14. Watch 9 and 14 even if you skip the rest: the indicator trick and LOTUS are where his delivery adds the most over the book.
- **Ross, *A First Course in Probability*, ch. 1–5 and 7.** Use as a drill source: more routine exercises than Blitzstein, with answers to odd problems. Warning: Ross's Geometric counts trials, not failures. Translate.
- **Feller, *An Introduction to Probability Theory and Its Applications*, vol. I, ch. I–II and IV.** Only if you want the classical treatment of the matching problem and the birthday problem with their history; Feller's prose is the best in the subject and his problem sets are the hardest.
