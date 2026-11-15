# F1 Exercises

Conventions as in the notes: $\mathbb{N} = \{1, 2, 3, \dots\}$; $i^2 = -1$; $\bar z$ is the conjugate, $|z|$ the modulus, $\arg z \in (-\pi, \pi]$ the principal argument; phasors use $e^{+i\omega t}$. Every "prove" means a written proof in sentences, with the technique named in the first line.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Let $U = \{1, 2, 3, 4, 5, 6\}$, $A = \{1, 2, 3, 4\}$, $B = \{3, 4, 5\}$. Compute $A \cup B$, $A \cap B$, $A \setminus B$, $B^c$, and list every element of $A \times \{0, 1\}$. How many elements does $A \times B$ have? Finally, write the set $\{2, 4, 6, 8, \dots\}$ in set-builder notation in two different ways.

**A2.** For the statement $S$: "If a square matrix $A$ is invertible, then $\det A \ne 0$," write the converse, the inverse, and the contrapositive in words. Which of the four are logically equivalent to $S$? (You may take as known from LA-6 that $A$ is invertible if and only if $\det A \ne 0$; which of the four statements are *true*?) Repeat for "If $n$ is divisible by $4$, then $n$ is even," and this time identify which of the four are false, with a counterexample for each false one.

**A3.** Negate each statement mechanically using (3.1)–(3.3), in symbols and then in words.
(a) $\forall x \in \mathbb{R}\ \exists y \in \mathbb{R}:\ y^2 = x$. Then decide which is true, the statement or its negation, and give the witness.
(b) $\exists M\ \forall n \in \mathbb{N}:\ |a_n| \le M$ (the sequence is bounded).
(c) "Every nonzero vector $\mathbf v \in \mathbb{R}^2$ has a unique unit vector in its direction." (Write "unique" as $\exists \mathbf u\,[P(\mathbf u) \wedge \forall \mathbf u'\,(P(\mathbf u') \Rightarrow \mathbf u' = \mathbf u)]$ before negating.)

**A4.** Let $z = 2 + 3i$ and $w = 1 - i$. Compute $zw$, $z/w$ (in the form $a + bi$), $\bar z$, $|z|$, and verify $|z|^2 = z\bar z$ and $|zw| = |z||w|$ numerically. Then write $-1 + i$, $-2i$, and $1 - \sqrt3\, i$ in polar form $re^{i\theta}$ with $\theta$ the principal argument, and use de Moivre (Theorem 3.16) to compute $(1 - \sqrt3\, i)^6$. Explain in one sentence why the answer is real.

**A5.** (a) Evaluate $\sum_{i=1}^{4}\sum_{j=1}^{3}(i + j)$ twice: inner sum over $j$ first, then with the order swapped by (3.5). (b) Simplify $\sum_{k=1}^{n}\delta_{k3}\,a_k$ (for $n \ge 3$) and $\sum_{j=1}^{n}\delta_{ij}\delta_{jk}$. What matrix identity is the second one? (c) Evaluate $\sum_{k=1}^{10}\frac{1}{k(k+1)}$ by telescoping (3.9).

## B. Core (do all)

**B1.** (Functions.) (a) Prove that $f: \mathbb{R} \to \mathbb{R}$, $f(x) = 3x - 7$, is bijective, using the definitional test of §3.2, and find $f^{-1}$. (b) Show that $g: \mathbb{R} \to \mathbb{R}$, $g(x) = x^2$, is neither injective nor surjective, with explicit witnesses. (c) Choose a domain and codomain on which $x \mapsto x^2$ is bijective and *prove* it is (the injectivity half needs an argument, not a picture).

**B2.** (Composition.) Let $f: A \to B$ and $g: B \to C$. (a) Prove: if $f$ and $g$ are surjective, then $g \circ f$ is surjective. (b) Prove: if $g \circ f$ is injective, then $f$ is injective. (c) Give explicit finite sets and functions for which $g \circ f$ is injective but $g$ is not. (So the converse of Example 4.3 is only half true.) Project Part 8 tests (a)–(c) on random functions.

**B3.** (Induction, code.) Prove by induction that for all $n \in \mathbb{N}$,
$$\sum_{k=1}^{n} k^3 = \Big(\frac{n(n+1)}{2}\Big)^2.$$
Mark the line where the inductive hypothesis is used. Then add this identity to the dictionary `identities` in Part 1 of `project.py` so the checker verifies it symbolically and by brute force for $n \le 50$.

**B4.** (Cases and set equality.) (a) Prove that for every integer $n$, $6$ divides $n^3 - n$. (Hint: $n^3 - n = (n-1)n(n+1)$. Show $2 \mid n^3 - n$ by cases on the parity of $n$, show $3 \mid n^3 - n$ by cases on the remainder of $n$ mod $3$, then show that if $2 \mid x$ and $3 \mid x$ then $6 \mid x$.) (b) Prove De Morgan's law $(A \cup B)^c = A^c \cap B^c$ by two inclusions. (c) Prove $A \times (B \cap C) = (A \times B) \cap (A \times C)$ by two inclusions.

**B5.** (Existence and uniqueness.) (a) Prove that between any two distinct rational numbers $a < b$ there is an irrational number. (Hint: consider $a + (b - a)/\sqrt2$; you may use Theorem 3.7.) Is your proof constructive? (b) Prove that a function has at most one inverse (write the full argument, naming the justification for each equality). (c) Prove that the additive identity of $\mathbb{C}$ is unique: if $u \in \mathbb{C}$ satisfies $z + u = z$ for every $z \in \mathbb{C}$, then $u = 0$.

**B6.** (Index gymnastics; the LA-1 bridge.) Prove each identity, citing (3.4)–(3.9) at each step.
(a) $\displaystyle\sum_{i=1}^{n}\sum_{j=i}^{n} a_j = \sum_{j=1}^{n} j\,a_j$.
(b) $\displaystyle\sum_{i=1}^{n}\sum_{j=1}^{n}\delta_{ij}\,a_i b_j = \sum_{i=1}^{n} a_i b_i$.
(c) $\displaystyle\Big(\sum_{i=1}^{n} a_i\Big)^2 = \sum_{i=1}^{n} a_i^2 + 2\sum_{1 \le i < j \le n} a_i a_j$.
(d) With $(AB)_{ik} = \sum_j a_{ij}b_{jk}$ and $(A\mathbf x)_i = \sum_j a_{ij}x_j$ (LA-1 (3.2), (3.3)), prove $\big((AB)\mathbf x\big)_i = \big(A(B\mathbf x)\big)_i$ for every $i$, in index notation.

**B7.** (Complex algebra and roots.) (a) Prove $|zw| = |z||w|$ using only $|z|^2 = z\bar z$ and (3.11), without polar form. (b) Find all four fourth roots of $-16$ in the form $a + bi$ and check one of them by direct multiplication. (c) Prove that the set of $n$-th roots of unity is closed under multiplication and under taking reciprocals. (d) Let $\omega = e^{2\pi i/n}$ with $n \ge 2$. Prove that for every integer $m$,
$$\sum_{k=0}^{n-1}\omega^{km} = \begin{cases} n & \text{if } n \mid m, \\ 0 & \text{otherwise.}\end{cases}$$
(This is the orthogonality relation behind the discrete Fourier transform in F3.)

**B8.** (Forced oscillator by phasors; mechanical engineering.) A mass–spring–damper driven by a sinusoidal force obeys $m\ddot x + c\dot x + kx = F_0\cos\omega t$. Look for a steady-state solution $x_p(t) = \operatorname{Re}[X e^{i\omega t}]$ with complex amplitude $X$.
(a) Substitute $x = Xe^{i\omega t}$ into the *complex* equation $m\ddot x + c\dot x + kx = F_0e^{i\omega t}$, use (3.19), and derive
$$X = \frac{F_0}{k - m\omega^2 + i c\omega}.$$
Then explain, using the argument of Theorem 3.21, why taking real parts gives a solution of the original real equation.
(b) For $m = 1$, $c = 0.5$, $k = 4$, $F_0 = 1$, $\omega = 1.5$, compute $|X|$ and $\arg X$ to three decimals, and write $x_p(t)$ as $A\cos(\omega t + \varphi)$. Does the mass lead or lag the force? Check your numbers against project Part 9.
(c) In one or two sentences: why does this method depend on the ODE being linear *and* having constant coefficients?

## C. Challenge (do at least one; C1 recommended)

**C1.** ($\varepsilon$–$\delta$ and quantifier order.) (a) Prove from the definition that $f(x) = 2x + 1$ is continuous at every $a \in \mathbb{R}$. (b) Prove from the definition that $f(x) = x^2$ is continuous at every $a \in \mathbb{R}$. (Hint: $|x^2 - a^2| = |x - a||x + a|$; first restrict $|x - a| < 1$ to bound $|x + a|$.) Note where your $\delta$ depends on $a$. (c) Prove that $f(x) = x^2$ is *not* uniformly continuous on $\mathbb{R}$: write the negation of the definition mechanically, then win the game. (Take $\varepsilon = 1$; given $\delta$, choose $a$ and $x$, depending on $\delta$, with $|x - a| < \delta$ and $|x^2 - a^2| \ge 1$.) Say in one sentence how this relates to Cauchy's 1821 error.

**C2.** (Strong induction, code.) (a) Prove that every integer $n \ge 12$ can be written as $4a + 5b$ with $a, b$ nonnegative integers. (b) The Fibonacci numbers are $F_1 = F_2 = 1$, $F_n = F_{n-1} + F_{n-2}$. Prove $F_n < 2^n$ for all $n \ge 1$. Why does ordinary induction (one previous case) not suffice? (c) Add a brute-force check of (a) for all $12 \le n \le 200$ to `project.py`, and confirm that $n = 11$ fails.

**C3.** (Geometry of roots of unity.) (a) Prove the reverse triangle inequality $\big||z| - |w|\big| \le |z - w|$ from Theorem 3.22. (b) Let $\omega = e^{2\pi i/n}$. Prove that consecutive roots of unity are equidistant, $|\omega^{k+1} - \omega^k| = |\omega - 1| = 2\sin(\pi/n)$ for every $k$, so the roots are the vertices of a regular $n$-gon inscribed in the unit circle, and prove its area is $\frac{n}{2}\sin\frac{2\pi}{n}$. Check the area numerically for $n = 3, \dots, 8$ in the project's "Your turn" (the shoelace formula is one line). (c) What is $\lim_{n\to\infty}\frac{n}{2}\sin\frac{2\pi}{n}$, and why is that the right answer?

**C4.** (Real polynomials and the odd-degree theorem.) Using the Fundamental Theorem of Algebra (Theorem 3.23) and the conjugate-roots theorem (Theorem 3.24): (a) Prove that every polynomial with real coefficients factors over $\mathbb{R}$ into a constant times real linear factors $(z - r)$ and real quadratic factors $z^2 - 2\operatorname{Re}(r)\,z + |r|^2$ with $r \notin \mathbb{R}$. (Strong induction on the degree; you may use that dividing a real polynomial by a real polynomial gives a real quotient.) (b) Deduce that every real polynomial of odd degree has a real root. (c) Deduce that every real $n \times n$ matrix with $n$ odd has at least one real eigenvalue, given that its characteristic polynomial $\det(A - \lambda I)$ is a real polynomial of degree $n$ (LA-6). Give a geometric example in $\mathbb{R}^3$.
