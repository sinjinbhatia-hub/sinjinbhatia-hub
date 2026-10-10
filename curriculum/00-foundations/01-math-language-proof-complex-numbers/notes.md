# F1: Mathematical Language, Proof, and Complex Numbers

**Time:** ~6 hours · **Prerequisites:** none · **Unlocks:** LA-3, LA-7, F4, and every proof in every later module

> **Conventions in this module.** Sets are uppercase ($A$, $B$), elements lowercase. $\mathbb{N} = \{1, 2, 3, \dots\}$ (no zero; Hammack's convention), $\mathbb{Z}$ integers, $\mathbb{Q}$ rationals, $\mathbb{R}$ reals, $\mathbb{C}$ complex numbers. Logic: $\neg$ not, $\wedge$ and, $\vee$ or (inclusive), $\Rightarrow$ implies, $\Leftrightarrow$ if and only if, $\forall$ for all, $\exists$ there exists. A function is written $f: A \to B$. Sums are $\sum_{k=m}^{n} a_k$ with explicit bounds; $\delta_{ij}$ is the Kronecker delta. Complex numbers: $i^2 = -1$ (engineers write $j$; this curriculum uses $i$ everywhere). $z = x + iy$ with $x = \operatorname{Re} z$ and $y = \operatorname{Im} z$ real; conjugate $\bar z = x - iy$; modulus $|z| = \sqrt{x^2 + y^2}$; argument $\theta = \arg z$ in radians, counterclockwise from the positive real axis, principal value in $(-\pi, \pi]$. Time dependence: phasors in this module use $e^{+i\omega t}$ (the engineering convention). Physics modules (QM-2, EM-7) use $e^{-i\omega t}$. The sign flip is flagged in §3.8 and §8.

## 0. Where this sits

You can compute. You can differentiate, integrate, solve ODEs, and write Python that works. What nobody has taught you is what a *proof* is, what a *definition* does, or how to turn the words "for all" into a plan of action. That gap is invisible until LA-3, when the problems stop saying "compute" and start saying "show that," and it costs students weeks. This module closes it in advance.

The second half of the module is complex numbers, treated as geometry: multiplication is rotation and scaling, and $e^{i\theta}$ is the point on the unit circle at angle $\theta$. You need that fluency for complex eigenvalues (LA-7), Fourier series and transforms (F2, F3), every wavefunction in QM, and every phase in QC. It also happens to be how vibration analysis and AC circuits are done in your own department.

This module needs nothing. LA-3, LA-7, and F4 need it, and so does every "show that" after this point.

## 1. The problem that forced this idea

**Problem one: the student who computes perfectly and fails the first proof.** Here is what happens. A student who has aced three semesters of calculus meets the problem "Show that if $f$ and $g$ are injective, then $g \circ f$ is injective." They write $g(f(x)) = g(f(y))$, stare at it, and either write nothing or write a paragraph of true facts that never connects to the question. The student is not weak. They have never been told that a proof is a specific kind of object with specific moves, that "injective" is a definition that *unpacks into a sentence you can use*, and that "show that for all $x$" means "begin by writing *Let $x$ be arbitrary*." Those are three facts. They take an afternoon to learn and nobody teaches them.

The idea of proof is old. Euclid's *Elements* (around 300 BCE) is the template: a handful of definitions and postulates, then propositions, each justified only by definitions, postulates, and earlier propositions. For two thousand years that was what "proof" meant and nobody needed to formalize it further. What forced the modern precision was calculus. Cauchy's *Cours d'analyse* (1821) set out to make the subject rigorous and still contains a famous error: it asserts that a convergent sum of continuous functions is continuous. Abel exhibited a Fourier-series counterexample in 1826. The mistake is a *quantifier-order* mistake: whether a tolerance $\delta$ is allowed to depend on the point $x$ or must work for all $x$ at once. Fixing it took the field until Weierstrass's lectures in the 1860s, which gave the $\varepsilon$–$\delta$ definitions we use today. The lesson is that the order of "for every" and "there exists" is not pedantry; it is the entire content of some theorems, and getting it wrong fooled Cauchy.

**Problem two: a formula that gives the right answer through nonsense.** In 1545 Cardano's *Ars Magna* published the solution of the cubic $x^3 = px + q$ (del Ferro's and Tartaglia's method). Applied to $x^3 = 15x + 4$ the formula produces

$$
x = \sqrt[3]{2 + \sqrt{-121}} + \sqrt[3]{2 - \sqrt{-121}},
$$

which contains square roots of a negative number, even though $x = 4$ is a perfectly real solution (check: $64 = 60 + 4$). Cardano had met $\sqrt{-15}$ in another problem and called such quantities "as subtle as they are useless." Rafael Bombelli, in his *L'Algebra* (1572), decided to compute with them anyway. He wrote down the rules we would write as $(+i)(+i) = -1$ and showed that $\sqrt[3]{2 + 11i} = 2 + i$ (check: $(2+i)^3 = 8 + 12i - 6 - i = 2 + 11i$), so the formula gives $(2 + i) + (2 - i) = 4$. The nonsense cancelled and the right answer came out. Descartes (1637) called the quantities "imaginary" as a dismissal, and the name stuck.

For nearly two hundred years they were a trick that worked. Two things changed that. Euler, in *Introductio in analysin infinitorum* (1748), found $e^{i\theta} = \cos\theta + i\sin\theta$ (he wrote $\sqrt{-1}$; the symbol $i$ is his too, from 1777), which turned multiplication of complex numbers into addition of angles. And Wessel (1797), Argand (1806), and Gauss (1831) drew them as points in a plane, making the rotation visible. Gauss coined the name "complex number." Hamilton (1833–37) removed the last mystery by *defining* a complex number as an ordered pair of reals with a specified multiplication rule, so that $i$ is not a square root of anything mysterious; it is the pair $(0, 1)$. That is the definition we use.

## 2. Intuition

### A proof is a chain, and every link has a name

A computation is a chain of equalities: $(1+i)^2 = 1 + 2i + i^2 = 2i$. Each link is "$=$". A proof is also a chain, but the links come in more kinds:

```
  hypothesis ──(by definition)──► unpacked form ──(algebra)──► new fact
      ──(by Theorem 3.4)──► another fact ──(by definition)──► conclusion
```

Every link cites its justification: a definition, an earlier theorem, a rule of algebra, or a hypothesis you were given. The chain is written in *sentences*, because the links are not all "$=$" and the reader needs to know which kind each one is. "Since $n$ is odd, $n = 2k + 1$ for some integer $k$" is a link of type "by definition." "$n^2 = 4k^2 + 4k + 1$" is a link of type "algebra." A page of symbols with no words is not a proof, because the reader cannot see the links.

You should be able to guess, now, what the first sentence of a proof of "for all $x \in A$, $P(x)$" is. It is "Let $x \in A$." You cannot check every $x$, so you take one, refuse to assume anything about it beyond $x \in A$, and show $P(x)$. Because you assumed nothing special, the argument works for all of them. And the first move of a proof of "there exists $x$ with $P(x)$" is to produce one and check it.

### Quantifiers are moves in a two-player game

Read $\forall x\, \exists y\; R(x, y)$ as a game. Your opponent picks $x$. Then you, knowing their $x$, pick $y$. You win if $R(x, y)$ holds. The statement is true if you have a winning strategy.

Read $\exists y\, \forall x\; R(x, y)$ as the same game with the moves swapped: *you* pick $y$ first, then the opponent picks $x$ knowing your $y$. This is harder to win. "Everyone has a mother" ($\forall$ person $\exists$ mother) is true. "Someone is everyone's mother" ($\exists$ mother $\forall$ person) is not. The second always implies the first (if one $y$ works for every $x$, then certainly each $x$ has some $y$), and never the other way round.

You can guess the negation rule from the game. The negation of a statement is the statement that the *other* player wins. So to negate, swap who moves: $\forall \leftrightarrow \exists$, and negate the winning condition. $\neg(\forall x\, \exists y\; R)$ is $\exists x\, \forall y\; \neg R$. That is a mechanical rule, and §3.4 makes it official.

### Multiplication by $i$ is a quarter turn

Draw the complex number $x + iy$ as the point $(x, y)$. Multiplying by $i$ sends $1 \mapsto i$, i.e., $(1, 0) \mapsto (0, 1)$: a counterclockwise quarter turn. It sends $i \mapsto i^2 = -1$, i.e., $(0, 1) \mapsto (-1, 0)$: another quarter turn. So "$i^2 = -1$" is not mysterious. Two quarter turns make a half turn, and a half turn sends $1$ to $-1$.

```
            i
            |
     -1 ----+---- 1        multiply by i: rotate 90° counterclockwise
            |              i·i = half turn = −1
           -i
```

Now guess the general rule. Multiplying by $2i$ should double the length and turn by $90°$. Multiplying by $-1$ turns by $180°$. Multiplying by $1 + i$, which sits at $45°$ with length $\sqrt 2$, should turn by $45°$ and stretch by $\sqrt 2$. The guess: **multiplying by $z$ rotates by the angle of $z$ and scales by the length of $z$.** Angles add; lengths multiply.

If angles add under multiplication, then the function "point on the unit circle at angle $\theta$" turns multiplication into addition of $\theta$. The only functions that do that are exponentials. So that point should be $e^{\text{something}\cdot\theta}$, and the something, by the quarter-turn picture, is $i$. That is Euler's formula, and §3.8 derives it two ways.

## 3. Formal development

### 3.1 Sets

**Definition 3.1.** A **set** is a collection of objects; $x \in A$ means $x$ is an element of $A$. Sets are written by listing, $\{1, 2, 3\}$, or by **set-builder notation**: $\{x \in A : P(x)\}$ is the set of elements of $A$ for which the statement $P(x)$ is true. The **empty set** is $\varnothing$. For sets $A, B$ inside a universe $U$:

- $A \subseteq B$ (**subset**) means $\forall x\,(x \in A \Rightarrow x \in B)$.
- $A = B$ means $A \subseteq B$ and $B \subseteq A$. *This is why every proof of set equality has two halves.*
- $A \cup B = \{x : x \in A \vee x \in B\}$ (**union**), $A \cap B = \{x : x \in A \wedge x \in B\}$ (**intersection**), $A \setminus B = \{x : x \in A \wedge x \notin B\}$ (**difference**), $A^c = U \setminus A$ (**complement**).
- $A \times B = \{(a, b) : a \in A,\ b \in B\}$ (**Cartesian product**), the set of ordered pairs. $\mathbb{R}^2 = \mathbb{R} \times \mathbb{R}$. If $A$ has $m$ elements and $B$ has $n$, then $A \times B$ has $mn$.

Example: the even integers are $\{n \in \mathbb{Z} : n = 2k \text{ for some } k \in \mathbb{Z}\}$. The sentence after the colon is the *definition of even*, and when a proof says "let $n$ be even," the next line is "so $n = 2k$ for some integer $k$." Unpacking a definition means writing down the sentence after the colon.

### 3.2 Functions

**Definition 3.2.** A **function** $f: A \to B$ assigns to each $a \in A$ exactly one element $f(a) \in B$. $A$ is the **domain**, $B$ the **codomain**, and $f(A) = \{f(a) : a \in A\} \subseteq B$ the **image**. The codomain is part of the function: $x \mapsto x^2$ as $\mathbb{R} \to \mathbb{R}$ and as $\mathbb{R} \to [0, \infty)$ are different functions (the second is surjective; the first is not).

**Definition 3.3.** Let $f: A \to B$.
- $f$ is **injective** (one-to-one) if $\forall a_1, a_2 \in A:\ f(a_1) = f(a_2) \Rightarrow a_1 = a_2$. Distinct inputs give distinct outputs.
- $f$ is **surjective** (onto) if $\forall b \in B\ \exists a \in A:\ f(a) = b$. Every element of the codomain is hit.
- $f$ is **bijective** if both.

```
 injective, not surjective     surjective, not injective      bijective
   A        B                    A        B                   A        B
   1 ───► a                      1 ───► a                     1 ───► a
   2 ───► b                      2 ───► a                     2 ───► b
          c                      3 ───► b                     3 ───► c
```

Forget the "horizontal line test." It only works for functions $\mathbb{R} \to \mathbb{R}$ you can graph, and it hides the logic. Use the **definitional test**, which works for every function ever:

- To prove injective: write "Suppose $f(a_1) = f(a_2)$." Derive $a_1 = a_2$.
- To prove surjective: write "Let $b \in B$." Construct an $a \in A$ with $f(a) = b$ and check it.
- To disprove either: one counterexample. Not injective: exhibit $a_1 \ne a_2$ with $f(a_1) = f(a_2)$. Not surjective: exhibit a $b$ nobody hits.

The **composition** of $f: A \to B$ and $g: B \to C$ is $g \circ f: A \to C$, $(g \circ f)(a) = g(f(a))$. The **identity** $\mathrm{id}_A: A \to A$ is $\mathrm{id}_A(a) = a$. An **inverse** of $f: A \to B$ is a $g: B \to A$ with $g \circ f = \mathrm{id}_A$ and $f \circ g = \mathrm{id}_B$.

**Theorem 3.4 (Inverses exist exactly for bijections).** $f: A \to B$ has an inverse if and only if $f$ is bijective.

▶ **Try it first.** This is an "if and only if," so it is two proofs. For ($\Leftarrow$), you have to *define* $g(b)$: which $a$ should it be, and why is there exactly one candidate? For ($\Rightarrow$), apply $g$ to both sides of $f(a_1) = f(a_2)$.

**Derivation.**

1. ($\Leftarrow$) Assume $f$ is bijective. Let $b \in B$. By surjectivity there is some $a$ with $f(a) = b$; by injectivity there is only one (if $f(a) = f(a') = b$ then $a = a'$). Define $g(b)$ to be that unique $a$. This is a well-defined function $B \to A$.
2. Check $f \circ g = \mathrm{id}_B$: $f(g(b)) = b$ because $g(b)$ was chosen to satisfy exactly that.
3. Check $g \circ f = \mathrm{id}_A$: for $a \in A$, $g(f(a))$ is the unique element mapping to $f(a)$, and $a$ is such an element, so $g(f(a)) = a$.
4. ($\Rightarrow$) Assume $g$ is an inverse. Injective: if $f(a_1) = f(a_2)$, apply $g$: $g(f(a_1)) = g(f(a_2))$, so $a_1 = a_2$ since $g \circ f = \mathrm{id}_A$. Surjective: given $b$, set $a = g(b)$; then $f(a) = f(g(b)) = b$. ∎

■ **What just happened.** "Invertible" is a single word that unpacks into exactly two checks, existence (surjective) and uniqueness (injective) of a preimage. In LA-4 the same theorem reads: a linear map is invertible iff its kernel is $\{\mathbf 0\}$ (injective) and its image is the whole codomain (surjective).

### 3.3 Logic: connectives and implication

A **statement** is a sentence that is either true or false. From statements $P, Q$ we build $\neg P$, $P \wedge Q$, $P \vee Q$ (true if at least one holds), and the one that causes trouble, $P \Rightarrow Q$ ("if $P$ then $Q$"):

| $P$ | $Q$ | $P \Rightarrow Q$ |
|---|---|---|
| T | T | T |
| T | F | **F** |
| F | T | T |
| F | F | T |

The only way an implication fails is $P$ true and $Q$ false. When $P$ is false the implication is true by default: a promise whose condition never arose was not broken. This is the right convention because we want "if $n > 2$ then $n^2 > 4$" to be a true statement about *all* integers, including $n = 1$, where the hypothesis is false. Consequence: to prove $P \Rightarrow Q$ you may assume $P$; the case $\neg P$ needs no work.

**Definition 3.5 (The four relatives of an implication).** For $P \Rightarrow Q$:

| Name | Statement | Equivalent to $P \Rightarrow Q$? |
|---|---|---|
| (original) | $P \Rightarrow Q$ | yes |
| **contrapositive** | $\neg Q \Rightarrow \neg P$ | **yes** |
| **converse** | $Q \Rightarrow P$ | **no** |
| **inverse** | $\neg P \Rightarrow \neg Q$ | **no** (but it is equivalent to the converse) |

▶ **Try it first.** Build the truth table for $\neg Q \Rightarrow \neg P$ and compare with the table above. Then do $Q \Rightarrow P$.

**Derivation.** $\neg Q \Rightarrow \neg P$ is false only when $\neg Q$ is true and $\neg P$ is false, i.e., $Q$ false and $P$ true. That is the same single row where $P \Rightarrow Q$ is false. Same table, same statement. $Q \Rightarrow P$ is false only in the row $Q$ true, $P$ false, which is a different row, so it is a different statement. The inverse $\neg P \Rightarrow \neg Q$ is the contrapositive of the converse, hence equivalent to the converse. ∎

Example: "If $A$ is invertible then $\det A \neq 0$." Contrapositive: "If $\det A = 0$ then $A$ is not invertible" (same fact). Converse: "If $\det A \ne 0$ then $A$ is invertible" (a *different* fact, which happens to also be true, LA-6). For "if $n$ is divisible by 4 then $n$ is even," the converse is false ($n = 2$). Equivalent statements are always both true or both false; the converse can go either way, and you must check it separately.

$P \Leftrightarrow Q$ means $(P \Rightarrow Q) \wedge (Q \Rightarrow P)$. **A proof of an "if and only if" is two proofs**, one for each direction, and you label them ($\Rightarrow$) and ($\Leftarrow$) as in Theorem 3.4. Also useful: $P \Rightarrow Q$ is the same as $\neg P \vee Q$ (check the table), which is what makes the negation rule (3.3) below work. Phrasebook: "$P$ only if $Q$" and "$P$ is sufficient for $Q$" and "$Q$ is necessary for $P$" all mean $P \Rightarrow Q$.

### 3.4 Quantifiers

$\forall x \in A,\ P(x)$ means $P$ holds for every element of $A$. $\exists x \in A,\ P(x)$ means at least one element satisfies $P$. Read them as game moves (§2): $\forall$ is the opponent's move, $\exists$ is yours.

**Theorem 3.6 (Order matters, one way).** $\exists y\, \forall x\; R(x, y)$ implies $\forall x\, \exists y\; R(x, y)$. The converse is false in general.

*Proof.* Suppose $y_0$ satisfies $R(x, y_0)$ for all $x$. Let $x$ be arbitrary. Then $y = y_0$ works for this $x$. So $\forall x\, \exists y\; R(x, y)$. For the converse: on $\mathbb{Z}$ with $R(x, y) = (y > x)$, $\forall x\, \exists y\, (y > x)$ is true (take $y = x + 1$) but $\exists y\, \forall x\, (y > x)$ is false (no integer exceeds all integers). ∎

**The $\varepsilon$–$\delta$ example.** You have met these in calculus; here they are as quantifier strings.

$f$ is **continuous at $a$**:
$$\forall \varepsilon > 0\ \ \exists \delta > 0\ \ \forall x\ \ \big(|x - a| < \delta \Rightarrow |f(x) - f(a)| < \varepsilon\big).$$
Game: opponent names a tolerance $\varepsilon$; you answer with a $\delta$; the opponent then picks any $x$ within $\delta$ of $a$; you win if $f(x)$ is within $\varepsilon$ of $f(a)$.

$f$ is **continuous on $\mathbb{R}$**: $\forall a\ \forall \varepsilon\ \exists \delta\ \forall x\ (\dots)$. The $\delta$ comes *after* $a$, so it may depend on $a$.

$f$ is **uniformly continuous on $\mathbb{R}$**: $\forall \varepsilon\ \exists \delta\ \forall a\ \forall x\ (\dots)$. The $\delta$ comes *before* $a$: one $\delta$ must work everywhere.

By Theorem 3.6, uniform continuity implies continuity. The converse fails: $f(x) = x^2$ is continuous on $\mathbb{R}$ but not uniformly (exercise C1, where you win the game explicitly). Cauchy's 1821 error was exactly this: he used a $\delta$ chosen before $x$ when his hypothesis only supplied one chosen after.

**Negating quantified statements mechanically.** Push $\neg$ inward, one symbol at a time, with these rules:

$$
\neg(\forall x\, P(x)) \equiv \exists x\, \neg P(x), \qquad \neg(\exists x\, P(x)) \equiv \forall x\, \neg P(x), \tag{3.1}
$$
$$
\neg(P \wedge Q) \equiv \neg P \vee \neg Q, \qquad \neg(P \vee Q) \equiv \neg P \wedge \neg Q, \tag{3.2}
$$
$$
\neg(P \Rightarrow Q) \equiv P \wedge \neg Q. \tag{3.3}
$$

Why (3.1): "not every $x$ has $P$" says "some $x$ lacks $P$"; "no $x$ has $P$" says "every $x$ lacks $P$." Why (3.3): from $P \Rightarrow Q \equiv \neg P \vee Q$ and (3.2). The algorithm: the $\neg$ walks past each quantifier, flipping it, until it reaches the innermost predicate, where it negates an inequality or equation ($\neg(a < b)$ is $a \ge b$; $\neg(a = b)$ is $a \ne b$). Example 4.1 runs it on two real statements. The negation of continuity at $a$ is
$$\exists \varepsilon > 0\ \ \forall \delta > 0\ \ \exists x\ \ \big(|x - a| < \delta \ \wedge\ |f(x) - f(a)| \ge \varepsilon\big),$$
read: there is a tolerance that no $\delta$ achieves. To *disprove* a universal statement you prove its negation, which is existential: exhibit one counterexample. To disprove an existential statement you prove a universal one, which needs a general argument.

### 3.5 Proof techniques

Each technique below comes with a skeleton, one annotated example, and one bad proof that students actually write, corrected. The annotations in brackets are the kind of link (§2).

**(a) Direct proof.** To prove $P \Rightarrow Q$: assume $P$, unpack it, compute or reason, arrive at $Q$.

*Claim.* If $n$ is odd, then $n^2$ is odd.
*Proof.* Let $n$ be an odd integer. [The $\forall$ became "let."] Then $n = 2k + 1$ for some integer $k$. [Definition of odd.] So $n^2 = 4k^2 + 4k + 1 = 2(2k^2 + 2k) + 1$. [Algebra.] Since $2k^2 + 2k$ is an integer, $n^2$ is odd. [Definition of odd, read backwards.] ∎

*Bad proof.* "$n$ is odd, so $n^2 = n \cdot n$ is odd times odd, which is odd." — This uses "odd times odd is odd," which is the claim in disguise. Nothing was unpacked. *Fix:* write $n = 2k + 1$ and compute. The definition is the only tool you are allowed to start from.

**(b) Proof by contrapositive.** To prove $P \Rightarrow Q$, prove $\neg Q \Rightarrow \neg P$ instead (equivalent by Definition 3.5). Use it when $\neg Q$ is easier to unpack than $P$: "irrational," "not even," "nonzero" are negatives, and their negations are positive and concrete.

*Claim.* For real $x$, if $x^2$ is irrational then $x$ is irrational.
*Proof.* We prove the contrapositive: if $x$ is rational then $x^2$ is rational. [Announce the technique.] Let $x$ be rational, so $x = p/q$ with $p, q$ integers, $q \ne 0$. [Definition.] Then $x^2 = p^2/q^2$ with $p^2, q^2$ integers and $q^2 \ne 0$. [Algebra.] So $x^2$ is rational. [Definition.] ∎

*Bad proof.* "Suppose $x$ is irrational. Then $x$ is not $p/q$, so $x^2$ is not $p^2/q^2$…" — This is the *inverse* $\neg P \Rightarrow \neg Q$, not equivalent to the claim (and the step is false anyway: $\sqrt2$ is irrational, $(\sqrt2)^2 = 2$). *Fix:* negate the conclusion, not the hypothesis. Write the contrapositive out explicitly before you start.

**(c) Proof by contradiction.** To prove $P$: assume $\neg P$, derive something impossible. Use when $\neg P$ gives you something concrete to compute with.

**Theorem 3.7.** $\sqrt 2$ is irrational.

▶ **Try it first.** Assume $\sqrt 2 = p/q$ *in lowest terms*. Square. What does the equation say about the parity of $p$? Then of $q$?

*Proof.* Suppose, for contradiction, that $\sqrt 2$ is rational. Then $\sqrt 2 = p/q$ with $p, q \in \mathbb{Z}$, $q \ne 0$, and we may take the fraction in lowest terms, so $p$ and $q$ have no common factor. [Definition of rational, plus a free normalization.] Squaring, $2q^2 = p^2$. [Algebra.] So $p^2$ is even, hence $p$ is even (Example 4.4). [Earlier result.] Write $p = 2m$. Then $2q^2 = 4m^2$, so $q^2 = 2m^2$ is even, hence $q$ is even. [Same result again.] Now $p$ and $q$ are both even, so both are divisible by 2, contradicting lowest terms. [The contradiction.] Hence $\sqrt 2$ is irrational. ∎

*Bad proof.* The same argument without "in lowest terms." It ends with "$p$ and $q$ are both even," and stops, as if that were impossible. It is not: $\sqrt2 = 2\sqrt2/2$ has even numerator and denominator if you believe the premise. *Fix:* the contradiction has to contradict *something you assumed or know*. "Lowest terms" is what makes "both even" impossible. Every contradiction proof must say exactly what it contradicts.

**(d) Proof by cases.** Split the possibilities into cases that *together cover everything*, and prove the claim in each.

*Claim.* For every integer $n$, $n^2 + n$ is even.
*Proof.* Let $n \in \mathbb{Z}$. Every integer is even or odd. [The cases are exhaustive.] Case 1: $n = 2k$. Then $n^2 + n = 4k^2 + 2k = 2(2k^2 + k)$, even. Case 2: $n = 2k+1$. Then $n^2 + n = 4k^2 + 4k + 1 + 2k + 1 = 2(2k^2 + 3k + 1)$, even. In both cases $n^2 + n$ is even. ∎

*Bad proof.* "$n = 1$: $2$, even. $n = 2$: $6$, even. $n = 3$: $12$, even. So it is always even." — Examples are not cases. Cases must be a complete partition of the domain, each handled by a general argument. *Fix:* say what the cases are and why they exhaust the integers.

**(e) Proof by induction.** To prove $P(n)$ for all integers $n \ge n_0$: prove $P(n_0)$ (**base case**); then prove that for any $m \ge n_0$, $P(m) \Rightarrow P(m+1)$ (**inductive step**). The assumption $P(m)$ in the step is the **inductive hypothesis**, and the step must *use* it.

**Theorem 3.8.** For all $n \in \mathbb{N}$, $\displaystyle\sum_{k=1}^{n} k = \frac{n(n+1)}{2}$.

*Proof.* Base case $n = 1$: left side $1$, right side $1 \cdot 2/2 = 1$. ✓ Inductive step: fix $m \ge 1$ and assume $\sum_{k=1}^{m} k = m(m+1)/2$. [Inductive hypothesis, for *one* $m$.] Then
$$\sum_{k=1}^{m+1} k = \Big(\sum_{k=1}^{m} k\Big) + (m+1) = \frac{m(m+1)}{2} + (m+1) = \frac{(m+1)(m+2)}{2},$$
where the second equality *is* the inductive hypothesis. The right side is the formula with $n = m+1$. So $P(m) \Rightarrow P(m+1)$, and by induction $P(n)$ holds for all $n \ge 1$. ∎

*Bad proof.* "Assume the formula holds for all $n$. Then for $n+1$, $\sum_{k=1}^{n+1} k = (n+1)(n+2)/2$ by the formula. ∎" — Two errors. "For all $n$" assumes the conclusion. And the step never splits off the last term, so the hypothesis is never used; the formula for $n+1$ is just asserted. *Fix:* assume $P(m)$ for one fixed $m$, split $\sum^{m+1} = \sum^{m} + (\text{last term})$, and point to the line where the hypothesis enters. If you cannot point to it, you have not done induction.

**(f) Strong induction.** Same, but the inductive hypothesis is $P(k)$ for *all* $n_0 \le k \le m$, not just $P(m)$. Use it when the case $m+1$ depends on cases smaller than $m$.

**Theorem 3.9.** Every integer $n \ge 2$ is a product of primes (a prime counting as a product of one prime).

*Proof.* Base: $2$ is prime. Step: fix $m \ge 2$ and assume every $k$ with $2 \le k \le m$ is a product of primes. Consider $m+1$. If $m+1$ is prime, done. Otherwise $m+1 = ab$ with $2 \le a, b \le m$. [Definition of not prime.] By the strong hypothesis, $a$ and $b$ are each products of primes, so $ab$ is too. ∎

*Bad proof.* Ordinary induction: "Assume $m$ is a product of primes. Then $m + 1$…" — Nothing about the factorization of $m$ tells you anything about $m+1$ (consider $m = 8$, $m + 1 = 9$). The factors of $m+1$ are *smaller* numbers, not $m$. *Fix:* assume the claim for everything below $m+1$.

**(g) Existence proofs.** To prove $\exists x\, P(x)$: either exhibit $x$ and verify $P(x)$ (**constructive**), or show that no $x$ is impossible (**non-constructive**).

*Constructive.* There exist integers $x, y$ with $3x + 5y = 1$: take $x = 2$, $y = -1$; check $6 - 5 = 1$. ∎

*Non-constructive.* There exist irrational $a, b$ with $a^b$ rational. Consider $\sqrt2^{\sqrt2}$. Case 1: it is rational; take $a = b = \sqrt2$ (irrational by Theorem 3.7). Case 2: it is irrational; take $a = \sqrt2^{\sqrt2}$, $b = \sqrt2$; then $a^b = \sqrt2^{\sqrt2\cdot\sqrt2} = \sqrt2^2 = 2$. In either case such a pair exists, and we never learned which case holds. ∎

*Bad proof.* "Let $a$ be an irrational number such that $a^b$ is rational. Then…" — The first sentence assumes the thing to be proved exists. *Fix:* produce the object, or argue that its absence is impossible.

**(h) Uniqueness proofs.** To prove "there is at most one $x$ with $P(x)$": assume $x_1$ and $x_2$ both satisfy $P$, and show $x_1 = x_2$.

*Claim.* A function has at most one inverse.
*Proof.* Let $g$ and $h$ both be inverses of $f: A \to B$. Then
$g = g \circ \mathrm{id}_B = g \circ (f \circ h) = (g \circ f) \circ h = \mathrm{id}_A \circ h = h.$ ∎
(Each equality: definition of identity; $h$ is an inverse; composition is associative; $g$ is an inverse; definition of identity.)

*Bad proof.* "The inverse is unique because it is defined as *the* function that undoes $f$." — Using the word "the" is not a proof that there is only one. *Fix:* assume two, show equal.

**(i) Proving set equality.** Two inclusions, each an element-chase.

*Claim.* $A \cap (B \cup C) = (A \cap B) \cup (A \cap C)$.
*Proof.* ($\subseteq$) Let $x \in A \cap (B \cup C)$. Then $x \in A$ and $x \in B \cup C$, so $x \in B$ or $x \in C$. If $x \in B$ then $x \in A \cap B$; if $x \in C$ then $x \in A \cap C$. Either way $x \in (A \cap B) \cup (A \cap C)$. ($\supseteq$) Let $x \in (A \cap B) \cup (A \cap C)$. Then $x \in A \cap B$ or $x \in A \cap C$. In both cases $x \in A$, and $x \in B$ or $x \in C$, so $x \in B \cup C$. Hence $x \in A \cap (B \cup C)$. ∎

*Bad proof.* A Venn diagram with the same region shaded twice. — A diagram is how you *find* the proof, not the proof; it covers one configuration of three sets and cannot handle, say, $17$ sets or infinite ones. *Fix:* chase an arbitrary element through both inclusions.

**(j) Proving a function is injective or surjective.** Use the definitional test from §3.2.

*Claim.* $f: \mathbb{R} \to \mathbb{R}$, $f(x) = 2x + 3$, is bijective.
*Proof.* Injective: suppose $f(a_1) = f(a_2)$. Then $2a_1 + 3 = 2a_2 + 3$, so $a_1 = a_2$. Surjective: let $y \in \mathbb{R}$. Set $x = (y - 3)/2 \in \mathbb{R}$. Then $f(x) = 2 \cdot \tfrac{y-3}{2} + 3 = y$. ∎

*Bad proof (surjectivity).* "$f$ is surjective because for every $x$ there is a $y$, namely $y = 2x + 3$." — That is the definition of *function*, with the quantifiers in the wrong order. Surjectivity starts from an arbitrary $y$ and must find an $x$. *Fix:* start with "Let $y \in \mathbb{R}$," and solve for $x$.

### 3.6 How to read a definition, and how to start when stuck

**A definition is a checklist.** Unpack it into numbered conditions before trying to use it. Here is the one you will meet first in LA-3.

> *Definition (LA-3 preview).* A subset $W$ of a vector space $V$ is a **subspace** if (i) $\mathbf 0 \in W$; (ii) $\mathbf u, \mathbf w \in W \Rightarrow \mathbf u + \mathbf w \in W$; (iii) $\mathbf u \in W,\ c \in \mathbb{R} \Rightarrow c\mathbf u \in W$.

To prove $W = \{(x, y) \in \mathbb{R}^2 : x + 2y = 0\}$ is a subspace, you do three things, one per item: (i) $0 + 2\cdot 0 = 0$ ✓. (ii) Let $(x_1, y_1), (x_2, y_2) \in W$, so $x_1 + 2y_1 = 0$ and $x_2 + 2y_2 = 0$; then $(x_1 + x_2) + 2(y_1 + y_2) = 0 + 0 = 0$ ✓. (iii) $cx_1 + 2cy_1 = c(x_1 + 2y_1) = 0$ ✓. To show $W' = \{(x, y) : x + 2y = 1\}$ is *not* a subspace, you fail one item: $(0,0) \notin W'$ since $0 \ne 1$. One failure is enough; three checks are required for success.

**When stuck:**
1. Write *what you know*: every hypothesis, each unpacked by its definition into equations or sentences.
2. Write *what you want*: the conclusion, unpacked the same way. If it is "$\forall x$," write "Let $x$." If it is "$\exists x$," write "We need to find $x$ such that…". If it is "$A = B$ for sets," write the two inclusions.
3. Work from both ends until they meet. Most short proofs are two or three lines down from the top and two or three lines up from the bottom.
4. If the conclusion is a negative ("irrational," "not injective," "$\ne$"), try the contrapositive or contradiction.
5. Try $n = 1, 2, 3$ or a $2 \times 2$ case. Small cases show you the mechanism; the proof is the mechanism written for general $n$.
6. Name things. "Some integer" becomes $k$. "A vector in $W$" becomes $\mathbf w$ with its defining equation written out.

### 3.7 Summation notation and index manipulation

$\sum_{k=m}^{n} a_k = a_m + a_{m+1} + \dots + a_n$ (and is $0$ if $n < m$). The index $k$ is a **bound** (dummy) variable, like the $x$ in $\int f(x)\,dx$: it has no meaning outside the sum and can be renamed freely, $\sum_k a_k = \sum_j a_j$. Rules, each checkable by writing out the terms:

$$
\sum_{k=m}^{n} (c\,a_k + d\,b_k) = c\sum_{k=m}^{n} a_k + d\sum_{k=m}^{n} b_k, \qquad
\sum_{k=m}^{n} a_k = \sum_{k=m}^{p} a_k + \sum_{k=p+1}^{n} a_k, \qquad
\sum_{k=m}^{n} a_k = \sum_{j=m+s}^{n+s} a_{j-s}. \tag{3.4}
$$

The last is an **index shift**: substitute $k = j - s$ and shift both bounds. Constants that do not involve the index factor out; constants that *do* involve it do not.

**Double sums.** $\sum_{i=1}^{m}\sum_{j=1}^{n} a_{ij}$ adds up an $m \times n$ table. Adding by rows or by columns gives the same total:

$$
\sum_{i=1}^{m}\sum_{j=1}^{n} a_{ij} = \sum_{j=1}^{n}\sum_{i=1}^{m} a_{ij}. \tag{3.5}
$$

For a **triangular** region the bounds change when you swap:

$$
\sum_{i=1}^{n}\sum_{j=1}^{i} a_{ij} = \sum_{j=1}^{n}\sum_{i=j}^{n} a_{ij}. \tag{3.6}
$$

▶ **Try it first.** Draw the lattice points $(i, j)$ with $1 \le j \le i \le n$ for $n = 4$. Describe the same set as "for each $j$, which $i$?"

**Derivation of (3.6).** Both sides sum $a_{ij}$ over the same finite set of index pairs $T = \{(i, j) : 1 \le j \le i \le n\}$, so they are equal; the only question is the description of $T$. The left side reads $T$ row by row: for each $i$ from $1$ to $n$, $j$ runs $1, \dots, i$. The right side reads it column by column: for a fixed $j$, the condition $j \le i \le n$ says $i$ runs $j, \dots, n$. Same set, same sum. ∎

■ **What just happened.** Swapping a double sum is never algebra; it is redescribing a region. Draw the region, read off the new bounds. Example 4.5 does this with $a_{ij} = 1$.

**Product of two sums.**
$$
\Big(\sum_{i=1}^{m} a_i\Big)\Big(\sum_{j=1}^{n} b_j\Big) = \sum_{i=1}^{m}\sum_{j=1}^{n} a_i b_j. \tag{3.7}
$$
*Derivation.* $\sum_j b_j$ is a single number, call it $B$. Then $\big(\sum_i a_i\big) B = \sum_i (a_i B)$ by linearity (3.4), and $a_i B = a_i \sum_j b_j = \sum_j a_i b_j$ again by linearity (the index $j$ does not appear in $a_i$). ∎ The indices *must* be different letters: $\big(\sum_i a_i\big)^2 = \sum_i\sum_j a_i a_j$, which is **not** $\sum_i a_i^2$ (exercise B6(c)).

**Kronecker delta.** $\delta_{ij} = 1$ if $i = j$ and $0$ otherwise. Its **sifting property**:
$$
\sum_{j=1}^{n} \delta_{ij}\, a_j = a_i \qquad (1 \le i \le n). \tag{3.8}
$$
*Derivation.* In the sum, every term with $j \ne i$ is $0 \cdot a_j = 0$; the single term $j = i$ is $1 \cdot a_i$. ∎ The identity matrix has entries $(I)_{ij} = \delta_{ij}$, and the basis vector $\mathbf e_j$ has entries $(\mathbf e_j)_i = \delta_{ij}$. So LA-1 Theorem 3.5, $A\mathbf e_j = \mathbf a_j$, is one line of (3.8): $(A\mathbf e_j)_i = \sum_k a_{ik}\delta_{kj} = a_{ij}$.

**Telescoping.**
$$
\sum_{k=1}^{n} (a_{k+1} - a_k) = a_{n+1} - a_1. \tag{3.9}
$$
*Derivation.* By linearity and the index shift (3.4) with $s = 1$: $\sum_{k=1}^{n} a_{k+1} = \sum_{j=2}^{n+1} a_j$. Subtract $\sum_{k=1}^{n} a_k$: the terms $a_2, \dots, a_n$ appear in both and cancel, leaving $a_{n+1} - a_1$. ∎ Example: $\frac{1}{k(k+1)} = \frac1k - \frac1{k+1}$, so $\sum_{k=1}^{n}\frac{1}{k(k+1)} = 1 - \frac{1}{n+1}$.

■ **What just happened.** These six rules are the whole toolkit behind the proofs of LA-1's (3.1), (3.2), (3.3): "the $i$-th entry of $\sum_j x_j \mathbf a_j$ is $\sum_j x_j a_{ij}$" is linearity; "$(AB)_{ik} = \sum_j a_{ij} b_{jk}$ makes $(AB)\mathbf x = A(B\mathbf x)$" is a swap of two finite sums (exercise B6(d)). If you can do B6 you will never be stopped by an index again.

### 3.8 Complex numbers

**Definition 3.10 (Hamilton).** $\mathbb{C} = \mathbb{R} \times \mathbb{R}$ with
$$(a, b) + (c, d) = (a + c,\ b + d), \qquad (a, b)\cdot(c, d) = (ac - bd,\ ad + bc).$$
Write $(a, 0)$ as the real number $a$ (the reals sit inside $\mathbb{C}$ as the first axis), and set $i = (0, 1)$.

**Derivation of $i^2 = -1$ and of the notation $a + bi$.** $i^2 = (0,1)(0,1) = (0\cdot 0 - 1\cdot 1,\ 0\cdot 1 + 1\cdot 0) = (-1, 0) = -1$. And $a + bi = (a, 0) + (b, 0)(0, 1) = (a, 0) + (b\cdot 0 - 0\cdot 1,\ b\cdot 1 + 0\cdot 0) = (a, 0) + (0, b) = (a, b)$. ∎ So $(a, b) = a + bi$, and the multiplication rule is what you get by expanding $(a + bi)(c + di)$ with $i^2 = -1$:
$$(a + bi)(c + di) = ac + adi + bci + bd\,i^2 = (ac - bd) + (ad + bc)i. \tag{3.10}$$
Nothing imaginary was invented; a pair of reals was given a multiplication rule, and $i$ is the name of the pair $(0, 1)$. Commutativity, associativity, and distributivity follow from those of $\mathbb{R}$ (check (3.10) is symmetric in the two factors).

**Definition 3.11.** For $z = a + bi$: $\operatorname{Re} z = a$, $\operatorname{Im} z = b$, the **conjugate** $\bar z = a - bi$, the **modulus** $|z| = \sqrt{a^2 + b^2}$.

**Theorem 3.12 (Conjugate algebra).** For all $z, w \in \mathbb{C}$:
$$
|z|^2 = z\bar z, \qquad \overline{z + w} = \bar z + \bar w, \qquad \overline{zw} = \bar z\,\bar w, \qquad \bar{\bar z} = z, \qquad z \in \mathbb{R} \Leftrightarrow z = \bar z. \tag{3.11}
$$

▶ **Try it first.** Multiply out $(a + bi)(a - bi)$. Then compute $\overline{zw}$ and $\bar z \bar w$ from (3.10) and compare.

**Derivation.** $z\bar z = (a + bi)(a - bi) = a^2 - (bi)^2 = a^2 + b^2 = |z|^2$. Sum: immediate from the definition. Product: by (3.10), $\overline{zw} = (ac - bd) - (ad + bc)i$, while $\bar z\bar w = (a - bi)(c - di) = (ac - bd) + (-ad - bc)i$, the same. Double conjugate: flip the sign of $b$ twice. Real iff conjugate-fixed: $a - bi = a + bi \Leftrightarrow b = 0$. ∎

Two immediate consequences. **Division:** for $z \ne 0$,
$$\frac{1}{z} = \frac{\bar z}{|z|^2}, \qquad \text{so} \qquad \frac{w}{z} = \frac{w\bar z}{|z|^2}, \tag{3.12}$$
because $z \cdot \bar z/|z|^2 = |z|^2/|z|^2 = 1$. (Multiply top and bottom by the conjugate.) **Real and imaginary parts:** $\operatorname{Re} z = (z + \bar z)/2$, $\operatorname{Im} z = (z - \bar z)/(2i)$.

**The complex plane.** Plot $a + bi$ at $(a, b)$. Addition is vector addition (parallelogram rule), immediately from Definition 3.10. $|z|$ is the distance from the origin and $|z - w|$ the distance between $z$ and $w$. Conjugation is reflection in the real axis.

**Polar form.** Write $a = r\cos\theta$, $b = r\sin\theta$ with $r = |z|$ and $\theta = \arg z$ (determined up to multiples of $2\pi$; the principal value lies in $(-\pi, \pi]$). Then $z = r(\cos\theta + i\sin\theta)$. Computing $\theta$ by hand: $\theta = \operatorname{atan2}(b, a)$, *not* $\arctan(b/a)$, which loses the quadrant ($-1 - i$ and $1 + i$ have the same $b/a$).

**Theorem 3.13 (Euler's formula).** For every real $\theta$,
$$e^{i\theta} = \cos\theta + i\sin\theta. \tag{3.13}$$

Here $e^w$ for complex $w$ is *defined* by the power series $e^w = \sum_{n=0}^{\infty} w^n/n!$, the same series as for real $w$. We take two analysis facts for granted: this series converges absolutely for every complex $w$ (compare with $\sum |w|^n/n!$), and an absolutely convergent series may be regrouped freely. Both are proved in any real analysis text (Abbott, *Understanding Analysis*, ch. 2). Everything else we derive.

▶ **Try it first.** Compute $i^0, i^1, i^2, i^3, i^4, i^5$. Then write out the first six terms of $\sum (i\theta)^n/n!$ and sort them into those with and without an $i$.

**Derivation 1 (power series).**

1. The powers of $i$ cycle with period 4: $i^0 = 1$, $i^1 = i$, $i^2 = -1$, $i^3 = -i$, $i^4 = 1$, …. So $i^{2k} = (-1)^k$ and $i^{2k+1} = (-1)^k i$.
2. Split the series by parity of $n$ (allowed by absolute convergence):
$$e^{i\theta} = \sum_{k=0}^{\infty}\frac{(i\theta)^{2k}}{(2k)!} + \sum_{k=0}^{\infty}\frac{(i\theta)^{2k+1}}{(2k+1)!} = \sum_{k=0}^{\infty}\frac{(-1)^k\theta^{2k}}{(2k)!} + i\sum_{k=0}^{\infty}\frac{(-1)^k\theta^{2k+1}}{(2k+1)!}.$$
3. The first sum is the Taylor series of $\cos\theta$ and the second is that of $\sin\theta$. ∎

**Derivation 2 (the ODE argument).** Both sides satisfy the same first-order initial value problem, and that problem has a unique solution.

1. Let $f(\theta) = \cos\theta + i\sin\theta$. Then $f'(\theta) = -\sin\theta + i\cos\theta = i(\cos\theta + i\sin\theta) = i f(\theta)$, and $f(0) = 1$.
2. Let $g(\theta) = e^{i\theta}$. Differentiating the series term by term, $g'(\theta) = \sum_{n\ge1} i^n \theta^{n-1}/(n-1)! = i\sum_{m\ge0} (i\theta)^m/m! = i\,g(\theta)$, and $g(0) = 1$.
3. Uniqueness. Let $u = f - g$. Then $u' = iu$ and $u(0) = 0$. Compute $\frac{d}{d\theta}|u|^2 = \frac{d}{d\theta}(u\bar u) = u'\bar u + u\bar u'$. Now $\bar u' = \overline{u'} = \overline{iu} = -i\bar u$ by (3.11). So $\frac{d}{d\theta}|u|^2 = iu\bar u - iu\bar u = 0$. Hence $|u|^2$ is constant, equal to $|u(0)|^2 = 0$, so $u \equiv 0$ and $f = g$. ∎

■ **What just happened.** Derivation 1 is bookkeeping; Derivation 2 is the one to remember. $f' = if$ says: the velocity of the point $e^{i\theta}$ is $i$ times its position, i.e., the position rotated by a quarter turn (§2), with the same length. A point whose velocity is always perpendicular to its position vector and of equal magnitude moves on a circle of radius $|f(0)| = 1$ at unit angular speed. So $e^{i\theta}$ is the point on the unit circle at angle $\theta$. That picture is the whole theory of the complex exponential. Corollaries, each one line from (3.13): $e^{i\pi} = -1$; $|e^{i\theta}|^2 = \cos^2\theta + \sin^2\theta = 1$; $\overline{e^{i\theta}} = e^{-i\theta}$; $e^{i\theta} = 1 \Leftrightarrow \theta \in 2\pi\mathbb{Z}$.

So every complex number is $z = re^{i\theta}$ with $r = |z|$, $\theta = \arg z$.

**Theorem 3.14 (Multiplication is rotation and scaling).** If $z = re^{i\theta}$ and $w = se^{i\varphi}$, then $zw = rs\,e^{i(\theta + \varphi)}$. Hence $|zw| = |z||w|$ and $\arg(zw) = \arg z + \arg w$ (mod $2\pi$).

▶ **Try it first.** Multiply $(\cos\theta + i\sin\theta)(\cos\varphi + i\sin\varphi)$ with (3.10). Which trig identities do you see?

**Derivation.** It suffices to show $e^{i\theta}e^{i\varphi} = e^{i(\theta + \varphi)}$. By (3.13) and (3.10),
$$e^{i\theta}e^{i\varphi} = (\cos\theta\cos\varphi - \sin\theta\sin\varphi) + i(\sin\theta\cos\varphi + \cos\theta\sin\varphi) = \cos(\theta + \varphi) + i\sin(\theta + \varphi) = e^{i(\theta+\varphi)},$$
using the angle-addition formulas. Then $zw = rs\,e^{i\theta}e^{i\varphi} = rs\,e^{i(\theta+\varphi)}$, whose modulus is $rs$ (since $|e^{i\psi}| = 1$) and whose argument is $\theta + \varphi$. ∎

■ **What just happened.** The guess from §2 is now a theorem. Multiplying by $w = se^{i\varphi}$ rotates the plane by $\varphi$ and scales it by $s$. The argument can be run backwards: if you *trust* the exponent law $e^{i\theta}e^{i\varphi} = e^{i(\theta+\varphi)}$ (which the series satisfies), Euler's formula *proves* the angle-addition formulas by comparing real and imaginary parts. You never have to memorize them again.

**Theorem 3.15 (The matrix of multiplication by $z$).** Identify $w = x + iy$ with the column $(x, y)^\top \in \mathbb{R}^2$. The map $w \mapsto zw$, for fixed $z = a + bi$, is linear, and its matrix is
$$M(z) = \begin{pmatrix} a & -b \\ b & a \end{pmatrix} = r\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}, \qquad z = re^{i\theta}. \tag{3.14}$$
Moreover $M(z)M(w) = M(zw)$ and $\det M(z) = a^2 + b^2 = |z|^2$.

▶ **Try it first.** The columns of a matrix are the images of $\mathbf e_1 = 1$ and $\mathbf e_2 = i$ (LA-1 Theorem 3.5). Compute $z \cdot 1$ and $z \cdot i$.

**Derivation.** Linearity: $z(w_1 + w_2) = zw_1 + zw_2$ and $z(cw) = c(zw)$ for real $c$, by distributivity. First column: $z\cdot 1 = a + bi \leftrightarrow (a, b)^\top$. Second column: $z \cdot i = ai + bi^2 = -b + ai \leftrightarrow (-b, a)^\top$. That is (3.14); the polar form follows from $a = r\cos\theta$, $b = r\sin\theta$. $M(z)M(w) = M(zw)$ because both are the matrix of $u \mapsto z(wu) = (zw)u$ (LA-1 §3.3: composition is matrix multiplication). The determinant is $a\cdot a - (-b)\cdot b$. ∎

■ **What just happened.** Complex numbers *are* the $2 \times 2$ real matrices of the form (3.14), and multiplication of complex numbers is multiplication of these matrices. The matrix with $r = 1$ is the rotation matrix $R_\theta$ you met in LA-1 Example 4.3 (there with $\theta = 90°$, which is $M(i)$). In LA-7 you will prove that any real $2 \times 2$ matrix with complex eigenvalues $a \pm bi$ is, in a suitable basis, exactly $M(a + bi)$: complex eigenvalues of a real matrix mean rotation-and-scaling, which is why they produce oscillation in ODE systems.

**Theorem 3.16 (de Moivre).** For every integer $n$, $(\cos\theta + i\sin\theta)^n = \cos n\theta + i\sin n\theta$.

*Proof.* The claim is $(e^{i\theta})^n = e^{in\theta}$. For $n \ge 0$, induction: $n = 0$ gives $1 = 1$; if $(e^{i\theta})^m = e^{im\theta}$ then $(e^{i\theta})^{m+1} = e^{im\theta}e^{i\theta} = e^{i(m+1)\theta}$ by Theorem 3.14. For $n < 0$, $(e^{i\theta})^{n} = 1/(e^{i\theta})^{|n|} = 1/e^{i|n|\theta} = e^{-i|n|\theta}$, using $e^{i\psi}e^{-i\psi} = e^{0} = 1$. ∎

Expanding $n = 2$ with (3.10) gives $\cos 2\theta = \cos^2\theta - \sin^2\theta$ and $\sin 2\theta = 2\sin\theta\cos\theta$ for free. Example 4.6 uses de Moivre to compute $(1 + i)^{10}$ in one line.

**Theorem 3.17 ($n$-th roots of unity).** For $n \in \mathbb{N}$, the equation $z^n = 1$ has exactly $n$ solutions,
$$\omega_k = e^{2\pi i k/n}, \qquad k = 0, 1, \dots, n-1, \tag{3.15}$$
equally spaced on the unit circle. Writing $\omega = \omega_1 = e^{2\pi i/n}$, they are $1, \omega, \omega^2, \dots, \omega^{n-1}$.

▶ **Try it first.** Write $z = re^{i\theta}$ and use de Moivre on $z^n$. Compare moduli, then arguments (remembering arguments are only defined mod $2\pi$).

**Derivation.**

1. Let $z = re^{i\theta}$ with $r > 0$. Then $z^n = r^n e^{in\theta}$ by Theorem 3.16. We need $r^n e^{in\theta} = 1 = 1\cdot e^{i\cdot 0}$.
2. Moduli: $r^n = 1$, and $r$ is a positive real, so $r = 1$.
3. Arguments: $e^{in\theta} = 1$ iff $n\theta = 2\pi k$ for some $k \in \mathbb{Z}$ (last corollary of Euler), so $\theta = 2\pi k/n$.
4. Distinctness: $k$ and $k + n$ give the same $\theta$ mod $2\pi$, so there are at most $n$ distinct points, namely $k = 0, \dots, n-1$; these are distinct because their angles differ by less than $2\pi$. Conversely each $\omega_k$ satisfies $\omega_k^n = e^{2\pi i k} = 1$. Finally $\omega_k = (e^{2\pi i/n})^k = \omega^k$. ∎

**Theorem 3.18 (The roots of unity sum to zero).** For $n \ge 2$, $\displaystyle\sum_{k=0}^{n-1}\omega^k = 0$.

*Proof.* Let $S = \sum_{k=0}^{n-1}\omega^k$. Multiply by $\omega$ and shift the index (3.4): $\omega S = \sum_{k=0}^{n-1}\omega^{k+1} = \sum_{k=1}^{n}\omega^k = S - \omega^0 + \omega^n = S - 1 + 1 = S$. So $(\omega - 1)S = 0$. Since $n \ge 2$, $\omega = e^{2\pi i/n} \ne 1$, so $S = 0$. ∎

■ **What just happened.** Geometrically: the vertices of a regular polygon centered at the origin have centroid at the origin. Algebraically: $S$ is a geometric series, and the proof is the standard "multiply by the ratio and subtract." In F3 this becomes the orthogonality of the discrete Fourier basis, which is what makes the DFT invertible (exercise B7(d)).

**Roots of a general complex number.** To solve $z^n = w$ with $w = se^{i\varphi} \ne 0$, the same argument gives $r = s^{1/n}$ (the positive real root) and $n\theta = \varphi + 2\pi k$:
$$z_k = s^{1/n}\,e^{i(\varphi + 2\pi k)/n}, \qquad k = 0, \dots, n-1. \tag{3.16}$$
One root at angle $\varphi/n$, the rest obtained by multiplying by the roots of unity $\omega^k$. Example 4.7 finds the cube roots of $8i$.

**Theorem 3.19 (cos and sin from exponentials).**
$$\cos\theta = \frac{e^{i\theta} + e^{-i\theta}}{2}, \qquad \sin\theta = \frac{e^{i\theta} - e^{-i\theta}}{2i}. \tag{3.17}$$
*Derivation.* Write (3.13) for $\theta$ and for $-\theta$, using $\cos(-\theta) = \cos\theta$, $\sin(-\theta) = -\sin\theta$: $e^{i\theta} = \cos\theta + i\sin\theta$ and $e^{-i\theta} = \cos\theta - i\sin\theta$. Add and divide by 2; subtract and divide by $2i$. ∎ This is why F2 writes Fourier series in $e^{inx}$ rather than in sines and cosines: one family of functions instead of two, and multiplication of them is addition of frequencies.

**Phasors.** A sinusoid of angular frequency $\omega$, amplitude $A$, and phase $\varphi$ is the real part of a rotating complex number:
$$A\cos(\omega t + \varphi) = \operatorname{Re}\big[A e^{i\varphi}\, e^{i\omega t}\big], \tag{3.18}$$
because $Ae^{i(\omega t + \varphi)} = A\cos(\omega t + \varphi) + iA\sin(\omega t + \varphi)$ by Euler, and we take the real part. The complex number $\tilde A = Ae^{i\varphi}$ is the **phasor**: it stores amplitude and phase, and the time dependence $e^{i\omega t}$ is common to everything at that frequency.

**Theorem 3.20 (Adding sinusoids of equal frequency).** $A_1\cos(\omega t + \varphi_1) + A_2\cos(\omega t + \varphi_2) = A\cos(\omega t + \varphi)$ where $Ae^{i\varphi} = A_1e^{i\varphi_1} + A_2e^{i\varphi_2}$. In words: add the phasors as vectors.

*Derivation.* $\operatorname{Re}$ is additive, $\operatorname{Re}(z_1 + z_2) = \operatorname{Re} z_1 + \operatorname{Re} z_2$ (Definition 3.10). So the left side is $\operatorname{Re}[\tilde A_1 e^{i\omega t}] + \operatorname{Re}[\tilde A_2 e^{i\omega t}] = \operatorname{Re}[(\tilde A_1 + \tilde A_2)e^{i\omega t}]$. Write $\tilde A_1 + \tilde A_2 = Ae^{i\varphi}$ in polar form; by (3.18) this is $A\cos(\omega t + \varphi)$. ∎

■ **What just happened.** Trigonometric addition problems became vector addition in the plane. This is how a vibration engineer adds the responses of a structure and how an electrical engineer adds voltages around an AC circuit: the time factor $e^{i\omega t}$ is dropped, the phasors are added, and the time factor is restored at the end. Example 4.8 does $3\cos\omega t + 4\cos(\omega t + \pi/2)$ and gets $5\cos(\omega t + 0.927)$; exercise B8 does a forced mass–spring–damper.

> **Sign convention warning.** Engineering writes $e^{+j\omega t}$ (and $j$ for $i$, since $i$ is current). Physics writes time dependence as $e^{-i\omega t}$: a plane wave is $e^{i(kx - \omega t)}$ (EM-7), and a stationary state evolves as $e^{-iEt/\hbar}$ (QM-2). Under the physics convention a phasor $Ae^{i\varphi}$ represents $A\cos(\omega t - \varphi)$, so phase *lags* and *leads* swap sign when you move between textbooks. Neither is wrong. Check which exponent a source uses before you read a phase off it.

**The complex exponential and linear ODEs.** For complex $\lambda = \alpha + i\beta$, define $e^{\lambda t} = e^{\alpha t}e^{i\beta t} = e^{\alpha t}(\cos\beta t + i\sin\beta t)$. Then
$$\frac{d}{dt}e^{\lambda t} = \lambda e^{\lambda t}. \tag{3.19}$$
*Derivation.* By the product rule, $\frac{d}{dt}\big[e^{\alpha t}(\cos\beta t + i\sin\beta t)\big] = \alpha e^{\alpha t}(\cos\beta t + i\sin\beta t) + e^{\alpha t}(-\beta\sin\beta t + i\beta\cos\beta t) = e^{\alpha t}\big[\alpha(\cos\beta t + i\sin\beta t) + i\beta(\cos\beta t + i\sin\beta t)\big] = (\alpha + i\beta)e^{\lambda t}$, using $i\beta(\cos + i\sin) = -\beta\sin + i\beta\cos$. ∎

So $e^{\lambda t}$ is an **eigenfunction** of $d/dt$ with eigenvalue $\lambda$: differentiation does not change its shape, only multiplies it. That is the entire reason exponentials solve constant-coefficient linear ODEs.

**Theorem 3.21.** For constants $a, b, c$ with $a \ne 0$, the function $y = e^{\lambda t}$ solves $ay'' + by' + cy = 0$ if and only if $a\lambda^2 + b\lambda + c = 0$. If $a, b, c$ are real and $y$ is a complex solution, then $\operatorname{Re} y$ and $\operatorname{Im} y$ are real solutions.

*Derivation.* By (3.19), $y' = \lambda e^{\lambda t}$ and $y'' = \lambda^2 e^{\lambda t}$, so $ay'' + by' + cy = (a\lambda^2 + b\lambda + c)e^{\lambda t}$. Since $e^{\lambda t} \ne 0$ (its modulus is $e^{\alpha t} > 0$), this vanishes iff the quadratic does. For the second claim, write $L[y] = ay'' + by' + cy$. Because $a, b, c$ are real and differentiation acts on real and imaginary parts separately, $\operatorname{Re} L[y] = L[\operatorname{Re} y]$. If $L[y] = 0$ then $L[\operatorname{Re} y] = \operatorname{Re} 0 = 0$; likewise for $\operatorname{Im}$. ∎

For the damped oscillator $\ddot x + 2\gamma\dot x + \omega_0^2 x = 0$ with $0 < \gamma < \omega_0$, the roots are $\lambda = -\gamma \pm i\omega_d$ with $\omega_d = \sqrt{\omega_0^2 - \gamma^2}$, so $e^{\lambda t} = e^{-\gamma t}(\cos\omega_d t + i\sin\omega_d t)$ and the real solutions are $e^{-\gamma t}\cos\omega_d t$ and $e^{-\gamma t}\sin\omega_d t$: decaying oscillation, with the decay rate in $\operatorname{Re}\lambda$ and the frequency in $\operatorname{Im}\lambda$. Example 4.9 and project Part 7 do this with numbers. In LA-7 the same idea becomes $e^{At}$ for a system $\dot{\mathbf x} = A\mathbf x$, with the eigenvalues of $A$ playing the role of $\lambda$.

**Unit modulus and phases (QC-1 preview).** Since $|e^{i\theta}| = 1$, multiplying by $e^{i\theta}$ preserves modulus: $|e^{i\theta}z| = |z|$. In QC-1 a qubit state is a pair $(\alpha, \beta) \in \mathbb{C}^2$ with $|\alpha|^2 + |\beta|^2 = 1$, and measurement probabilities are $|\alpha|^2$, $|\beta|^2$. Multiplying the whole state by a **global phase** $e^{i\theta}$ changes neither, so it is physically invisible; multiplying only $\beta$ by $e^{i\theta}$ (a **relative phase**) changes the state. Operators that preserve all moduli are *unitary* (LA-8), and $e^{i\theta}$ is the $1 \times 1$ case.

**Theorem 3.22 (Triangle inequality).** $|z + w| \le |z| + |w|$ for all $z, w \in \mathbb{C}$.

▶ **Try it first.** Expand $|z + w|^2 = (z + w)\overline{(z + w)}$ with (3.11). You will meet $z\bar w + \bar z w$. What is it in terms of $\operatorname{Re}$?

**Derivation.**

1. $|z + w|^2 = (z + w)(\bar z + \bar w) = z\bar z + z\bar w + w\bar z + w\bar w = |z|^2 + (z\bar w + \overline{z\bar w}) + |w|^2$, since $w\bar z = \overline{z\bar w}$ by (3.11).
2. $u + \bar u = 2\operatorname{Re} u$, so $|z + w|^2 = |z|^2 + 2\operatorname{Re}(z\bar w) + |w|^2$.
3. For any $u$, $\operatorname{Re} u \le |u|$, because $|u|^2 = (\operatorname{Re} u)^2 + (\operatorname{Im} u)^2 \ge (\operatorname{Re} u)^2$. And $|z\bar w| = |z||\bar w| = |z||w|$ by Theorem 3.14.
4. So $|z + w|^2 \le |z|^2 + 2|z||w| + |w|^2 = (|z| + |w|)^2$. Take square roots. ∎

■ **What just happened.** Same proof as LA-1 Theorem 3.12, with $\operatorname{Re}(z\bar w) \le |z||w|$ playing the role of Cauchy–Schwarz. In LA-8 the expression $z\bar w$ summed over components becomes the complex inner product, and this derivation becomes the triangle inequality in $\mathbb{C}^n$.

**Theorem 3.23 (Fundamental Theorem of Algebra).** Every polynomial $p(z) = a_nz^n + \dots + a_1z + a_0$ with complex coefficients and $n \ge 1$, $a_n \ne 0$, has at least one root in $\mathbb{C}$. Consequently $p(z) = a_n(z - r_1)(z - r_2)\cdots(z - r_n)$ for some $r_1, \dots, r_n \in \mathbb{C}$ (not necessarily distinct): $p$ has exactly $n$ roots counted with multiplicity.

The existence of one root is **not proved here.** Every known proof uses analysis or topology (the series definition of $e^z$ is not enough): see Needham, *Visual Complex Analysis*, ch. 7, for the winding-number argument, or any complex analysis text for the proof via Liouville's theorem. What we *can* derive is the consequence.

**Derivation of the factorization from the existence of one root.**

1. (Factor theorem.) If $p(r) = 0$ then $p(z) = (z - r)q(z)$ for a polynomial $q$ of degree $n - 1$ with leading coefficient $a_n$. Proof: $p(z) = p(z) - p(r) = \sum_{k=1}^{n} a_k(z^k - r^k)$, and each $z^k - r^k = (z - r)(z^{k-1} + z^{k-2}r + \dots + r^{k-1})$ (multiply out: it telescopes, (3.9)). Factor $(z - r)$ out of the whole sum.
2. Induct on $n$. For $n = 1$, $p(z) = a_1(z + a_0/a_1)$. For $n > 1$, FTA gives a root $r_n$; by step 1, $p = (z - r_n)q$ with $\deg q = n - 1$; by the inductive hypothesis $q = a_n\prod_{k=1}^{n-1}(z - r_k)$. ∎

**Why this forces LA-7 to work over $\mathbb{C}$.** The eigenvalues of an $n \times n$ matrix $A$ are the roots of its characteristic polynomial $\det(A - \lambda I)$, which has degree $n$. Over $\mathbb{R}$ a polynomial may have no roots at all: the rotation $M(i) = \begin{pmatrix}0 & -1\\1 & 0\end{pmatrix}$ has characteristic polynomial $\lambda^2 + 1$, so it has no real eigenvalue, which makes sense (a rotation fixes no direction). Over $\mathbb{C}$, Theorem 3.23 guarantees $n$ eigenvalues with multiplicity, so *every* matrix has at least one eigenvector, and that single fact is what makes triangularization (Schur form, LA-10) and the spectral theorem (LA-8) possible.

**Theorem 3.24 (Conjugate roots).** If $p$ has *real* coefficients and $p(z) = 0$, then $p(\bar z) = 0$.

▶ **Try it first.** Take the conjugate of both sides of $a_nz^n + \dots + a_0 = 0$ and push the bar inside with (3.11).

*Proof.* Using (3.11) repeatedly, $\overline{p(z)} = \overline{\sum_k a_k z^k} = \sum_k \bar a_k\,\bar z^{\,k} = \sum_k a_k\bar z^{\,k} = p(\bar z)$, where $\bar a_k = a_k$ because the coefficients are real. So $p(\bar z) = \overline{p(z)} = \bar 0 = 0$. ∎

Consequences. Non-real roots of a real polynomial come in conjugate pairs $r, \bar r$, and $(z - r)(z - \bar r) = z^2 - 2\operatorname{Re}(r)\,z + |r|^2$ is a quadratic with *real* coefficients. So every real polynomial factors over $\mathbb{R}$ into linear factors and irreducible quadratics, and a real polynomial of odd degree has a real root (exercise C4). For LA-7: the eigenvalues of a *real* matrix come in conjugate pairs $a \pm bi$, and each pair acts on a two-dimensional real subspace as the rotation-scaling $M(a + bi)$ of (3.14).

## 4. Worked examples

### Example 4.1: Negating two quantified statements

**(a)** "Every bounded sequence has a convergent subsequence." Formalize: $\forall$ sequences $(x_n)$: $\big[(x_n) \text{ bounded} \Rightarrow \exists$ subsequence $(x_{n_k})$ that converges$\big]$. Apply (3.1) then (3.3) then (3.1):
$$\exists (x_n):\ (x_n) \text{ bounded} \ \wedge\ \forall \text{ subsequences } (x_{n_k}):\ (x_{n_k}) \text{ does not converge}.$$
In words: "There is a bounded sequence none of whose subsequences converges." (For real sequences this negation is false; the original is the Bolzano–Weierstrass theorem. To disprove the original you would have to exhibit such a sequence; to prove it you must handle an arbitrary bounded sequence.)

**(b)** "$\forall \varepsilon > 0\ \exists \delta > 0\ \forall x\ (|x - a| < \delta \Rightarrow |f(x) - f(a)| < \varepsilon)$." Push the negation past each quantifier, flipping it, then apply (3.3) at the core and negate the inequality:
$$\exists \varepsilon > 0\ \ \forall \delta > 0\ \ \exists x\ \ \big(|x - a| < \delta\ \wedge\ |f(x) - f(a)| \ge \varepsilon\big).$$
Game reading: *you* now name an $\varepsilon$ first; the opponent answers with any $\delta$; you must then find an $x$ within $\delta$ of $a$ whose image is at least $\varepsilon$ away. If you can always do that, $f$ is discontinuous at $a$. Note "$<\varepsilon$" negates to "$\ge \varepsilon$," not "$> \varepsilon$."

### Example 4.2: Induction for $\sum_{k=1}^{n} k^2 = \dfrac{n(n+1)(2n+1)}{6}$

*Base case* $n = 1$: left side $1^2 = 1$; right side $1\cdot 2\cdot 3/6 = 1$. ✓

*Inductive step.* Fix $m \ge 1$ and assume $\sum_{k=1}^{m} k^2 = m(m+1)(2m+1)/6$. Then
$$\sum_{k=1}^{m+1} k^2 = \frac{m(m+1)(2m+1)}{6} + (m+1)^2 \quad [\text{inductive hypothesis used here}]$$
$$= \frac{(m+1)\big[m(2m+1) + 6(m+1)\big]}{6} = \frac{(m+1)(2m^2 + 7m + 6)}{6} = \frac{(m+1)(m+2)(2m+3)}{6},$$
since $(m+2)(2m+3) = 2m^2 + 3m + 4m + 6 = 2m^2 + 7m + 6$. The target formula at $n = m+1$ is $\frac{(m+1)(m+2)(2(m+1)+1)}{6} = \frac{(m+1)(m+2)(2m+3)}{6}$. Same. So the formula holds for all $n \ge 1$. ∎

The project checks this identity symbolically (SymPy's `summation` against the closed form), checks the inductive step symbolically (closed form at $m+1$ minus closed form at $m$ minus $(m+1)^2$ simplifies to $0$), and brute-forces $n \le 50$.

### Example 4.3: The composition of two injective functions is injective

Let $f: A \to B$ and $g: B \to C$ be injective. Claim: $g \circ f: A \to C$ is injective.

*Proof.* Suppose $(g \circ f)(a_1) = (g \circ f)(a_2)$ for some $a_1, a_2 \in A$. [Definitional test: assume equal outputs.] By definition of composition, $g(f(a_1)) = g(f(a_2))$. Since $g$ is injective, $f(a_1) = f(a_2)$. [Injectivity of $g$ applied to the inputs $f(a_1), f(a_2) \in B$.] Since $f$ is injective, $a_1 = a_2$. [Injectivity of $f$.] ∎

Notice the proof uses the hypotheses *in the order the functions are undone*: $g$ first (outermost), then $f$. Exercise B2 shows the converse is half true: if $g \circ f$ is injective then $f$ must be, but $g$ need not be. Project Part 8 tests this on random finite functions.

### Example 4.4: If $n^2$ is even then $n$ is even (contrapositive)

The conclusion "$n$ is even" is easy to negate concretely ("$n$ is odd, so $n = 2k+1$"), while the hypothesis "$n^2$ is even" is awkward to use directly (knowing $n^2 = 2m$ tells you little about $n$). So we prove the contrapositive: **if $n$ is odd then $n^2$ is odd.** Let $n = 2k + 1$. Then $n^2 = 4k^2 + 4k + 1 = 2(2k^2 + 2k) + 1$, which is odd. Hence if $n^2$ is even, $n$ is even. ∎

This is the lemma used twice inside the proof of Theorem 3.7.

### Example 4.5: Swap a double sum and evaluate $\sum_{i=1}^{3}\sum_{j=1}^{i} 1$

*Directly.* For $i = 1$ the inner sum has one term; for $i = 2$, two; for $i = 3$, three. Total $1 + 2 + 3 = 6$.

*Swapped, by (3.6).* The region is $\{(i, j) : 1 \le j \le i \le 3\}$:

```
  j=3:            (3,3)
  j=2:       (2,2)(3,2)
  j=1:  (1,1)(2,1)(3,1)
         i=1  i=2  i=3
```

Reading by columns (fixed $j$, $i$ from $j$ to $3$): $\sum_{j=1}^{3}\sum_{i=j}^{3} 1 = \sum_{j=1}^{3}(3 - j + 1) = 3 + 2 + 1 = 6$. ✓

In general $\sum_{i=1}^{n}\sum_{j=1}^{i} 1 = \sum_{i=1}^{n} i = n(n+1)/2$ (Theorem 3.8): the triangular region has $n(n+1)/2$ lattice points, and the swapped version gives $\sum_{j=1}^{n}(n - j + 1) = n + (n-1) + \dots + 1$, the same sum backwards.

### Example 4.6: $(1 + i)^{10}$ via polar form

$|1 + i| = \sqrt{1 + 1} = \sqrt 2$ and $\arg(1 + i) = \operatorname{atan2}(1, 1) = \pi/4$. So $1 + i = \sqrt2\, e^{i\pi/4}$ and by de Moivre (Theorem 3.16)
$$(1 + i)^{10} = (\sqrt2)^{10} e^{i\cdot 10\pi/4} = 2^5 e^{i5\pi/2} = 32\, e^{i\pi/2} = 32i,$$
using $e^{i5\pi/2} = e^{i(2\pi + \pi/2)} = e^{i\pi/2}$. Check without polar form: $(1 + i)^2 = 1 + 2i + i^2 = 2i$, so $(1 + i)^{10} = (2i)^5 = 32 i^5 = 32i$. ✓ The project computes this a third way, as $M(1+i)^{10}$ applied to $(1, 0)^\top$, and gets $(0, 32)^\top$.

### Example 4.7: All cube roots of $8i$

$8i = 8e^{i\pi/2}$, so $s = 8$, $\varphi = \pi/2$, $n = 3$ in (3.16): $z_k = 8^{1/3}e^{i(\pi/2 + 2\pi k)/3} = 2e^{i(\pi/6 + 2\pi k/3)}$.

- $k = 0$: $2e^{i\pi/6} = 2(\cos\tfrac{\pi}{6} + i\sin\tfrac{\pi}{6}) = 2(\tfrac{\sqrt3}{2} + \tfrac{i}{2}) = \sqrt3 + i$.
- $k = 1$: $2e^{i5\pi/6} = 2(-\tfrac{\sqrt3}{2} + \tfrac{i}{2}) = -\sqrt3 + i$.
- $k = 2$: $2e^{i3\pi/2} = 2(0 - i) = -2i$.

Checks. $(-2i)^3 = -8i^3 = -8(-i) = 8i$. ✓ $(\sqrt3 + i)^2 = 3 + 2\sqrt3 i + i^2 = 2 + 2\sqrt3 i$; then $(2 + 2\sqrt3 i)(\sqrt3 + i) = 2\sqrt3 + 2i + 6i + 2\sqrt3 i^2 = (2\sqrt3 - 2\sqrt3) + 8i = 8i$. ✓ The three roots form an equilateral triangle of circumradius $2$ with one vertex at $-2i$; the other two are $-2i$ rotated by $\pm 120°$.

### Example 4.8: Adding $3\cos\omega t + 4\cos(\omega t + \pi/2)$ with phasors (ME / AC circuits)

Phasors (3.18): $\tilde A_1 = 3e^{i0} = 3$ and $\tilde A_2 = 4e^{i\pi/2} = 4i$. Add as vectors (Theorem 3.20): $\tilde A = 3 + 4i$. Polar form: $|3 + 4i| = \sqrt{9 + 16} = 5$ and $\varphi = \operatorname{atan2}(4, 3) = 0.9273$ rad ($= 53.13°$). Therefore
$$3\cos\omega t + 4\cos(\omega t + \tfrac{\pi}{2}) = 5\cos(\omega t + 0.9273).$$
Check at $t = 0$: left side $3 + 4\cos(\pi/2) = 3$; right side $5\cos(0.9273) = 5\cdot\tfrac35 = 3$. ✓ Check the identity another way: $4\cos(\omega t + \pi/2) = -4\sin\omega t$, and $5\cos(\omega t + \varphi) = 5(\cos\omega t\cos\varphi - \sin\omega t\sin\varphi) = 5(0.6\cos\omega t - 0.8\sin\omega t) = 3\cos\omega t - 4\sin\omega t$. ✓

Mechanically: two displacements at the same frequency, one in phase and one a quarter cycle ahead, combine into a single oscillation of amplitude $5$, not $7$: the 3–4–5 triangle. The AC version: the voltage across a resistor and across an inductor in series are $90°$ apart, so the total amplitude is $\sqrt{V_R^2 + V_L^2}$. Project Part 6 plots the two phasors, their sum, and the three sinusoids, and asserts the identity on a time grid.

### Example 4.9: A damped oscillator by complex exponential (mechanics)

$\ddot x + \dot x + 4x = 0$, so $\gamma = 1/2$, $\omega_0 = 2$, with $x(0) = 1$, $\dot x(0) = 0$.

Ansatz $x = e^{\lambda t}$; Theorem 3.21 gives $\lambda^2 + \lambda + 4 = 0$, so $\lambda = \tfrac{-1 \pm \sqrt{1 - 16}}{2} = -\tfrac12 \pm i\tfrac{\sqrt{15}}{2}$. Thus $\omega_d = \sqrt{15}/2 = 1.9365$, and the real solutions are $e^{-t/2}\cos\omega_d t$ and $e^{-t/2}\sin\omega_d t$. General real solution $x = e^{-t/2}(C_1\cos\omega_d t + C_2\sin\omega_d t)$. From $x(0) = 1$: $C_1 = 1$. Differentiate: $\dot x = e^{-t/2}\big[-\tfrac12(C_1\cos\omega_d t + C_2\sin\omega_d t) + \omega_d(-C_1\sin\omega_d t + C_2\cos\omega_d t)\big]$, so $\dot x(0) = -\tfrac12 C_1 + \omega_d C_2 = 0$ gives $C_2 = \tfrac{1}{2\omega_d} = \tfrac{1}{\sqrt{15}} = 0.2582$. So
$$x(t) = e^{-t/2}\Big(\cos\tfrac{\sqrt{15}}{2}t + \tfrac{1}{\sqrt{15}}\sin\tfrac{\sqrt{15}}{2}t\Big).$$
Project Part 7 integrates the ODE numerically with `scipy.integrate.solve_ivp` and asserts agreement with this formula to within $10^{-6}$ on $[0, 10]$.

## 5. Why it matters and where it connects

- **LA-1.** The proofs that (3.1), (3.2), and (3.3) agree are the summation rules of §3.7; LA-1 Theorem 3.5 is the sifting property (3.8). Exercise B6(d) is $(AB)\mathbf x = A(B\mathbf x)$ in index form.
- **LA-3.** "Subspace," "linearly independent," "basis," and "spans" are all definitions-as-checklists (§3.6). Proving a set is a subspace is three checks; proving independence is "suppose $\sum c_i\mathbf v_i = \mathbf 0$ and show all $c_i = 0$," a direct proof whose first line is dictated by the definition. Set equality by two inclusions is how you prove two subspaces are equal.
- **LA-4.** Injective, surjective, bijective (Definition 3.3) become kernel, image, isomorphism. Theorem 3.4 is the statement that a linear map is invertible iff it is bijective; Rank–Nullity is what replaces counting when the sets are infinite.
- **LA-7.** Complex eigenvalues of a real matrix come in conjugate pairs (Theorem 3.24) and act as rotation-scalings $M(a + bi)$ (Theorem 3.15); the FTA (Theorem 3.23) is why every matrix has an eigenvalue over $\mathbb{C}$; $e^{\lambda t}$ as eigenfunction of $d/dt$ (3.19) becomes $e^{At}$.
- **LA-8.** The triangle inequality proof (Theorem 3.22) with $\sum_k z_k\bar w_k$ in place of $z\bar w$ is the triangle inequality for the complex inner product. $|e^{i\theta}| = 1$ is the $1 \times 1$ unitary.
- **F2, F3.** Fourier series are expansions in $e^{inx}$ (Theorem 3.19 is why one family suffices); the DFT matrix has entries $\omega^{jk}$ and Theorem 3.18 (with exercise B7(d)) is its orthogonality.
- **F4.** Probability's limit theorems are $\forall\varepsilon\,\exists N\,\forall n \ge N$ statements; "almost surely" versus "in probability" is a quantifier-order distinction, exactly §3.4's continuous versus uniformly continuous.
- **CM-6.** Normal-mode analysis is phasors with matrices: the ansatz $\mathbf x = \operatorname{Re}[\mathbf v e^{i\omega t}]$ turns the ODE system into an eigenproblem.
- **EM-7, QM-2.** Plane waves $e^{i(kx - \omega t)}$ and stationary states $e^{-iEt/\hbar}$ use the physics sign convention (§3.8 warning box).
- **QC-1.** Global versus relative phase is the modulus-preservation remark after Theorem 3.21.
- **ST-7.** An AR process is stationary iff the roots of its characteristic polynomial lie outside the unit circle: Theorem 3.23 and the complex plane, used for data.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 drill sets, truth tables, negation, complex arithmetic, and sums (do them all; they take an hour). B1–B8 are the module: B1–B2 (functions, with full definitional proofs), B3 (induction plus extending `project.py`), B6 (index manipulation, the LA-1 bridge), B7 (complex algebra and roots), and B8 (forced oscillator by phasors) are essential. C1 (uniform continuity: win the quantifier game), C2 (strong induction plus code), C3 (geometry of roots of unity), and C4 (real polynomials and the odd-degree theorem) go one step further; do at least C1.

## 7. Coding project

`project.py` checks every numerical claim in these notes and builds the pictures:

1. **Induction identity checker.** For $\sum k$ and $\sum k^2$: SymPy confirms the closed form symbolically, confirms the inductive step symbolically, and brute-forces $n \le 50$ with exact rationals. B3 asks you to add $\sum k^3$.
2. **Quantifier game.** On a finite set with a relation $R$, evaluates $\forall x\,\exists y\,R$ and $\exists y\,\forall x\,R$ and shows they differ; checks Theorem 3.6 and the negation rule (3.1) on 200 random relations.
3. **Roots of unity** for $n = 3, \dots, 8$: plotted, each checked to satisfy $z^n = 1$, and their sum asserted $\approx 0$ (Theorem 3.18).
4. **Euler's formula**, checked numerically at random $\theta$, and via partial sums of the series with the error plotted against the number of terms; the ODE $f' = if$ is checked by finite differences.
5. **Complex multiplication as the $2 \times 2$ matrix** (3.14): $M(z)\,\mathbf w = zw$, $M(z)M(w) = M(zw)$, $\det M(z) = |z|^2$ on random $z, w$; Example 4.6 via $M(1+i)^{10}$.
6. **Phasor addition** (Example 4.8): rotating vectors and the resulting sinusoid plotted, identity asserted against `np.cos` on a time grid.
7. **Damped oscillator** (Example 4.9): complex-exponential solution compared with `solve_ivp`.
8. **Finite functions**: `is_injective`, `is_surjective`, `compose` on dictionaries; Example 4.3 and exercise B2 tested on random functions, including a counterexample where $g \circ f$ is injective but $g$ is not.
9. **Worked-example numbers**: cube roots of $8i$, the double sum of Example 4.5, and the forced-oscillator phasor of B8.

**Build it yourself first:** before reading Parts 3 and 8, write `roots_of_unity(n)` (one line with `np.exp`) and `is_injective(f)` / `is_surjective(f, codomain)` for a function stored as a `dict`. Then compare.

**Your turn** (prompts at the end of the script): prove and check $\sum_{k=0}^{n}\binom{n}{k} = 2^n$; show the $n$-th roots of unity form a regular polygon and compute its area $(n/2)\sin(2\pi/n)$; write the matrix for multiplication by $e^{i\pi/4}$ and apply it to the unit square.

## 8. Where you'll struggle

- **Writing a proof as a computation with no words.** This is the single biggest failure, and it is not cosmetic. A column of equations does not tell the reader which lines are assumptions, which are definitions unpacked, and which follow from which. Every proof needs sentences: "Let…", "Suppose…", "By definition…", "Since $g$ is injective…", "Therefore…". If your proof has no verbs, it is not finished. Rule of thumb: a reader should be able to cover the symbols and still follow the logic from the words.
- **Confusing $A \Rightarrow B$ with $B \Rightarrow A$.** Everyone does this, including in LA-3 ("the vectors are independent, so the matrix is invertible" when the problem gave you invertibility). The fix is mechanical: before you start, write the hypothesis and the conclusion on separate lines and label them. If you catch yourself using the conclusion, you have proved the converse, which is a different (possibly false) theorem.
- **Induction proofs where the step never uses the hypothesis.** If you cannot point to the exact line where "by the inductive hypothesis" applies, you have either proved the result directly (fine, but then say so) or asserted it. Splitting off the last term is almost always the move that brings the hypothesis in.
- **Negation errors.** "$<$" negates to "$\ge$", not "$>$". "$\forall x\,\exists y$" negates to "$\exists x\,\forall y$", with the quantifiers flipped *and* in the same order; students flip the order too. Run the algorithm of §3.4 one symbol at a time and do not skip steps.
- **"Let $x$ be arbitrary" versus picking an $x$.** When proving $\forall x$, you may not choose $x$; the opponent does. If your proof of a universal statement contains a specific number, you have proved one case.
- **Complex numbers feel fake until the phasor example.** If you are impatient with §3.8, go to Example 4.8 and project Part 6 first, see the 3–4–5 triangle turn into $5\cos(\omega t + 0.927)$, and then come back for Euler. The abstraction becomes tolerable once you have seen it save you a page of trigonometry.
- **The $e^{+i\omega t}$ versus $e^{-i\omega t}$ convention.** Engineering (and this module) uses $e^{+j\omega t}$; physics uses $e^{-i\omega t}$, so QM-2's time evolution is $e^{-iEt/\hbar}$ and EM-7's wave is $e^{i(kx - \omega t)}$. A phase that "leads" in one convention "lags" in the other. Check the exponent before reading off a phase, and expect to translate.
- **$\arg$ is multivalued, and $\arctan(b/a)$ is wrong half the time.** Use $\operatorname{atan2}(b, a)$ (which is what `np.angle` does), and when taking $n$-th roots remember there are $n$ of them because $\varphi$ and $\varphi + 2\pi k$ are the same angle.
- **Exercise C1(c) is harder than it looks.** Negating "uniformly continuous" is mechanical; *winning* the resulting game (choosing $a$ and $x$ as functions of $\delta$) takes a real idea: push $a$ far out so that a small step in $x$ makes a big change in $x^2$.

## 9. Check yourself

1. Write the first sentence of a proof of "for all $x \in A$, $P(x)$." Then the first sentence of a proof that $f: A \to B$ is surjective.
2. State the contrapositive, converse, and inverse of "if $z^n = 1$ then $|z| = 1$." Which are equivalent to the original?
3. Negate, mechanically: $\forall \varepsilon > 0\ \exists N\ \forall n \ge N\ |a_n - L| < \varepsilon$.
4. What is the difference between $\forall a\,\forall\varepsilon\,\exists\delta$ and $\forall\varepsilon\,\exists\delta\,\forall a$? Name the two properties of functions these define.
5. In an induction proof, what is the inductive hypothesis, and what must the inductive step do with it?
6. Why does a proof of $A = B$ for sets have two parts?
7. Prove $\sqrt2$ is irrational, and say exactly what the contradiction contradicts.
8. Swap the order of $\sum_{i=1}^{n}\sum_{j=i}^{n} a_{ij}$. What is $\sum_{j}\delta_{ij}a_j$?
9. Derive $e^{i\theta} = \cos\theta + i\sin\theta$ from $f' = if$, $f(0) = 1$. What does $f' = if$ say geometrically?
10. Write the matrix of multiplication by $3 - 4i$. What are its determinant and its rotation angle (as an `atan2`)?
11. Find the fourth roots of unity and show they sum to zero two ways (algebra and picture).
12. Express $A\cos(\omega t + \varphi)$ as the real part of a complex exponential. Why does equal-frequency addition reduce to vector addition?
13. Why must LA-7 work over $\mathbb{C}$? Give a $2 \times 2$ real matrix with no real eigenvalues.

## 10. Pair with

- **Hammack, *Book of Proof*, 3rd ed. (free online), ch. 1–2, 4–10, 12.** Sets, logic, direct/contrapositive/contradiction, iff and existence proofs, set proofs, disproof, induction, functions. Far more exercises than here, with odd-numbered solutions; this is where you drill. Skip ch. 3 (counting) and 11 (relations) for now.
- **Velleman, *How to Prove It*, 3rd ed., ch. 1–3, 5–6.** Alternative to Hammack with a more explicit "givens and goals" method, which is §3.6's "write what you know, write what you want" in full. Better than Hammack on quantifier structure.
- **Needham, *Visual Complex Analysis*, ch. 1.** The geometry of complex multiplication, Euler's formula as uniform circular motion (our Derivation 2), and the history of Bombelli's cubic, told better than anyone. Read §1.I–1.III after §3.8.
- **3Blue1Brown, "Euler's formula with introductory group theory."** The rotation-and-scaling picture animated; watch after Theorem 3.14.
- **Welch Labs, *Imaginary Numbers Are Real* (video series).** Episodes 1–3 for the Cardano–Bombelli story and the plane picture, if you want the history with pictures.
