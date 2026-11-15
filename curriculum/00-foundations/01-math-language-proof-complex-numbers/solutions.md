# F1 Solutions

## A. Warm-up

**A1.** $A \cup B = \{1, 2, 3, 4, 5\}$. $A \cap B = \{3, 4\}$. $A \setminus B = \{1, 2\}$. $B^c = U \setminus B = \{1, 2, 6\}$.

$A \times \{0, 1\} = \{(1,0), (1,1), (2,0), (2,1), (3,0), (3,1), (4,0), (4,1)\}$, eight elements. $A \times B$ has $4 \cdot 3 = 12$ elements.

Two set-builder forms for $\{2, 4, 6, \dots\}$: $\{n \in \mathbb{N} : n \text{ is even}\}$, or $\{2k : k \in \mathbb{N}\}$, or spelled out, $\{n \in \mathbb{N} : \exists k \in \mathbb{N},\ n = 2k\}$.

**A2.** For $S$ = "if $A$ is invertible then $\det A \ne 0$":
- Converse: "If $\det A \ne 0$ then $A$ is invertible."
- Inverse: "If $A$ is not invertible then $\det A = 0$."
- Contrapositive: "If $\det A = 0$ then $A$ is not invertible."

Only the contrapositive is logically equivalent to $S$; the converse and inverse are equivalent to each other but not to $S$ (Definition 3.5). Since invertibility is *equivalent* to $\det A \ne 0$ (LA-6), both directions hold, so all four statements happen to be true. That is a fact about matrices, not about logic: it does not make the converse equivalent to $S$.

For "if $4 \mid n$ then $n$ is even":
- Converse: "If $n$ is even then $4 \mid n$." **False:** $n = 2$ is even and $4 \nmid 2$.
- Inverse: "If $4 \nmid n$ then $n$ is odd." **False:** $n = 2$ again ($4 \nmid 2$ but $2$ is even). The same counterexample kills both, as it must, since they are equivalent.
- Contrapositive: "If $n$ is odd then $4 \nmid n$." **True** (equivalent to the original: a multiple of $4$ is $2 \cdot 2k$, even).

**A3.** (a) Negation: $\exists x \in \mathbb{R}\ \forall y \in \mathbb{R}:\ y^2 \ne x$. In words: some real number is not the square of any real number. The negation is true, with witness $x = -1$: for every real $y$, $y^2 \ge 0 > -1$, so $y^2 \ne -1$. Hence the original statement is false.

(b) Negation: $\forall M\ \exists n \in \mathbb{N}:\ |a_n| > M$. In words: for every bound $M$, some term of the sequence exceeds it, i.e., the sequence is unbounded. (Note $\neg(|a_n| \le M)$ is $|a_n| > M$, strict.)

(c) Let $P(\mathbf u)$ be "$|\mathbf u| = 1$ and $\mathbf u = c\mathbf v$ for some $c > 0$." The statement is
$$\forall \mathbf v \ne \mathbf 0\ \ \exists \mathbf u\ \big[P(\mathbf u) \wedge \forall \mathbf u'\,(P(\mathbf u') \Rightarrow \mathbf u' = \mathbf u)\big].$$
Negate: flip $\forall \mathbf v$ to $\exists \mathbf v$, flip $\exists \mathbf u$ to $\forall \mathbf u$, and negate the bracket with (3.2) and (3.3):
$$\exists \mathbf v \ne \mathbf 0\ \ \forall \mathbf u\ \big[\neg P(\mathbf u) \ \vee\ \exists \mathbf u'\,(P(\mathbf u') \wedge \mathbf u' \ne \mathbf u)\big].$$
In words: there is a nonzero vector $\mathbf v$ such that every candidate $\mathbf u$ either fails to be a unit vector in the direction of $\mathbf v$, or is not the only one. Equivalently: some nonzero $\mathbf v$ has either no unit vector in its direction or at least two. (The original is true: $\mathbf u = \mathbf v/|\mathbf v|$ exists, and if $c\mathbf v$ has length $1$ with $c > 0$ then $c = 1/|\mathbf v|$, so it is unique.)

**A4.** $zw = (2 + 3i)(1 - i) = 2 - 2i + 3i - 3i^2 = 2 + i + 3 = 5 + i$.

$z/w$: multiply top and bottom by $\bar w = 1 + i$, (3.12): $\dfrac{(2 + 3i)(1 + i)}{|1 - i|^2} = \dfrac{2 + 2i + 3i + 3i^2}{2} = \dfrac{-1 + 5i}{2} = -\tfrac12 + \tfrac52 i$.

$\bar z = 2 - 3i$. $|z| = \sqrt{4 + 9} = \sqrt{13}$. Check: $z\bar z = (2 + 3i)(2 - 3i) = 4 - 9i^2 = 13 = |z|^2$. ✓ And $|zw| = |5 + i| = \sqrt{26}$, while $|z||w| = \sqrt{13}\sqrt{2} = \sqrt{26}$. ✓

Polar forms: $-1 + i$: $r = \sqrt2$, $\theta = \operatorname{atan2}(1, -1) = 3\pi/4$, so $\sqrt2\,e^{3\pi i/4}$. $-2i$: $r = 2$, $\theta = -\pi/2$, so $2e^{-i\pi/2}$. $1 - \sqrt3\,i$: $r = \sqrt{1 + 3} = 2$, $\theta = \operatorname{atan2}(-\sqrt3, 1) = -\pi/3$, so $2e^{-i\pi/3}$.

De Moivre: $(1 - \sqrt3\,i)^6 = 2^6 e^{-6i\pi/3} = 64\,e^{-2\pi i} = 64$. It is real because six copies of the angle $-60°$ make $-360°$, a full turn, which lands back on the positive real axis.

**A5.** (a) Inner sum first: $\sum_{j=1}^{3}(i + j) = 3i + (1 + 2 + 3) = 3i + 6$. Then $\sum_{i=1}^{4}(3i + 6) = 3(1 + 2 + 3 + 4) + 24 = 30 + 24 = 54$. Swapped: $\sum_{i=1}^{4}(i + j) = (1 + 2 + 3 + 4) + 4j = 10 + 4j$, then $\sum_{j=1}^{3}(10 + 4j) = 30 + 4 \cdot 6 = 54$. ✓

(b) $\sum_{k=1}^{n}\delta_{k3}a_k = a_3$ by sifting (3.8): every term with $k \ne 3$ is zero. $\sum_{j=1}^{n}\delta_{ij}\delta_{jk} = \delta_{ik}$: sift with $a_j = \delta_{jk}$ to get $\delta_{ik}$. In matrix language $(I \cdot I)_{ik} = \sum_j I_{ij}I_{jk}$, so this says $I^2 = I$.

(c) $\dfrac{1}{k(k+1)} = \dfrac1k - \dfrac1{k+1}$ (check: $\frac{(k+1) - k}{k(k+1)}$). So $\sum_{k=1}^{10}\big(\frac1k - \frac1{k+1}\big) = -\sum_{k=1}^{10}\big(\frac1{k+1} - \frac1k\big) = -\big(\frac1{11} - \frac11\big) = 1 - \frac1{11} = \frac{10}{11}$, by (3.9) with $a_k = 1/k$.

## B. Core

**B1.** (a) *Injective.* Suppose $f(a_1) = f(a_2)$. Then $3a_1 - 7 = 3a_2 - 7$; add $7$ and divide by $3$: $a_1 = a_2$. *Surjective.* Let $y \in \mathbb{R}$. Set $x = (y + 7)/3$, a real number. Then $f(x) = 3\cdot\frac{y+7}{3} - 7 = y$. So $f$ is bijective, and by Theorem 3.4 it has an inverse, $f^{-1}(y) = (y + 7)/3$. Check: $f^{-1}(f(x)) = (3x - 7 + 7)/3 = x$ and $f(f^{-1}(y)) = 3\cdot\frac{y+7}{3} - 7 = y$. ∎

(b) *Not injective:* $g(1) = 1 = g(-1)$ with $1 \ne -1$. *Not surjective:* $-1 \in \mathbb{R}$ is the codomain, but for every real $x$, $x^2 \ge 0 \ne -1$, so nothing maps to $-1$.

(c) Take $h: [0, \infty) \to [0, \infty)$, $h(x) = x^2$. *Injective.* Suppose $a, b \ge 0$ and $a^2 = b^2$. Then $0 = a^2 - b^2 = (a - b)(a + b)$, so $a = b$ or $a + b = 0$. In the second case, since $a, b \ge 0$, both are $0$, so again $a = b$. *Surjective.* Let $y \ge 0$. Set $x = \sqrt y$, which exists and lies in $[0, \infty)$; then $h(x) = y$. ∎ (Other correct choices: $(-\infty, 0] \to [0, \infty)$.)

**B2.** (a) Let $c \in C$. Since $g$ is surjective, there is $b \in B$ with $g(b) = c$. Since $f$ is surjective, there is $a \in A$ with $f(a) = b$. Then $(g \circ f)(a) = g(f(a)) = g(b) = c$. So every $c$ is hit, and $g \circ f$ is surjective. ∎ (The hypotheses are used in the reverse order from Example 4.3: surjectivity is unwound from the outside in.)

(b) Suppose $f(a_1) = f(a_2)$. Apply $g$ to both sides: $g(f(a_1)) = g(f(a_2))$, i.e., $(g \circ f)(a_1) = (g \circ f)(a_2)$. Since $g \circ f$ is injective, $a_1 = a_2$. Hence $f$ is injective. ∎

(c) $A = \{1, 2\}$, $B = \{a, b, c\}$, $C = \{1, 2\}$. Let $f(1) = a$, $f(2) = b$ and $g(a) = 1$, $g(b) = 2$, $g(c) = 2$. Then $g$ is not injective ($g(b) = g(c)$ with $b \ne c$), but $(g \circ f)(1) = 1$, $(g \circ f)(2) = 2$, which is injective. The element $c$ where $g$ fails is never reached by $f$, so the composition never sees the failure.

**B3.** *Base case* $n = 1$: left side $1^3 = 1$; right side $(1 \cdot 2/2)^2 = 1$. ✓

*Inductive step.* Fix $m \ge 1$ and assume $\sum_{k=1}^{m}k^3 = \big(\frac{m(m+1)}{2}\big)^2 = \frac{m^2(m+1)^2}{4}$. Then
$$\sum_{k=1}^{m+1}k^3 = \frac{m^2(m+1)^2}{4} + (m+1)^3 \qquad [\text{inductive hypothesis used here}]$$
$$= \frac{(m+1)^2\big[m^2 + 4(m+1)\big]}{4} = \frac{(m+1)^2(m^2 + 4m + 4)}{4} = \frac{(m+1)^2(m+2)^2}{4} = \Big(\frac{(m+1)(m+2)}{2}\Big)^2,$$
which is the formula at $n = m + 1$. By induction it holds for all $n \ge 1$. ∎

Code: in Part 1 of `project.py`, add the line

```python
identities["sum k^3 = (n(n+1)/2)^2"] = (k**3, (n * (n + 1) / 2) ** 2)
```

before the `for name, (term, closed) in identities.items()` loop. The checker then reports base case, inductive step, symbolic sum, and brute force for $n \le 50$ all passing. (Observation: $\sum k^3 = (\sum k)^2$.)

**B4.** (a) Let $n \in \mathbb{Z}$ and $x = n^3 - n = (n - 1)n(n + 1)$, the product of three consecutive integers.

*$2 \mid x$, by cases on parity.* If $n$ is even, $n = 2q$, then $x = (n-1)(2q)(n+1)$ is even. If $n$ is odd, $n = 2q + 1$, then $n - 1 = 2q$ and $x = (2q)\,n\,(n+1)$ is even. The cases are exhaustive.

*$3 \mid x$, by cases on $n$ mod $3$.* Every integer is $3q$, $3q + 1$, or $3q + 2$ for some integer $q$. If $n = 3q$, then $3 \mid n$. If $n = 3q + 1$, then $n - 1 = 3q$. If $n = 3q + 2$, then $n + 1 = 3(q + 1)$. In each case one of the three factors is divisible by $3$, so $x$ is.

*$2 \mid x$ and $3 \mid x$ imply $6 \mid x$.* Write $x = 2a$ and $x = 3b$ with $a, b \in \mathbb{Z}$. Then $x = 3x - 2x = 3(2a) - 2(3b) = 6(a - b)$, and $a - b$ is an integer. ∎

(b) ($\subseteq$) Let $x \in (A \cup B)^c$. Then $x \notin A \cup B$, i.e., $\neg(x \in A \vee x \in B)$. By (3.2) this is $x \notin A \wedge x \notin B$, i.e., $x \in A^c$ and $x \in B^c$, so $x \in A^c \cap B^c$. ($\supseteq$) Let $x \in A^c \cap B^c$. Then $x \notin A$ and $x \notin B$, so by (3.2) $\neg(x \in A \vee x \in B)$, i.e., $x \notin A \cup B$, so $x \in (A \cup B)^c$. ∎

(c) ($\subseteq$) Let $(p, q) \in A \times (B \cap C)$. Then $p \in A$ and $q \in B \cap C$, so $q \in B$ and $q \in C$. Hence $(p, q) \in A \times B$ and $(p, q) \in A \times C$, so $(p, q) \in (A \times B) \cap (A \times C)$. ($\supseteq$) Let $(p, q) \in (A \times B) \cap (A \times C)$. Then $(p, q) \in A \times B$, so $p \in A$ and $q \in B$; and $(p, q) \in A \times C$, so $q \in C$. Hence $q \in B \cap C$ and $(p, q) \in A \times (B \cap C)$. ∎

**B5.** (a) Let $a < b$ be rational and set $c = a + \dfrac{b - a}{\sqrt2}$. *Between:* since $b - a > 0$ and $0 < 1/\sqrt2 < 1$, we have $0 < (b - a)/\sqrt2 < b - a$, so $a < c < b$. *Irrational:* suppose for contradiction $c$ is rational. Then $c - a$ is rational (difference of rationals: $\frac pq - \frac rs = \frac{ps - qr}{qs}$) and nonzero (as $c > a$), so $\sqrt2 = \dfrac{b - a}{c - a}$ is a quotient of rationals with nonzero denominator, hence rational, contradicting Theorem 3.7. So $c$ is irrational. ∎ The proof is constructive: it names the number.

(b) Let $g$ and $h$ both be inverses of $f: A \to B$, so $g \circ f = \mathrm{id}_A$, $f \circ g = \mathrm{id}_B$, $h \circ f = \mathrm{id}_A$, $f \circ h = \mathrm{id}_B$. Then
$$g = g \circ \mathrm{id}_B \quad(\text{identity does nothing}) = g \circ (f \circ h) \quad(h \text{ is an inverse}) = (g \circ f) \circ h \quad(\text{composition is associative: both sides send } b \mapsto g(f(h(b)))) = \mathrm{id}_A \circ h \quad(g \text{ is an inverse}) = h.$$
So any two inverses are equal. ∎

(c) Suppose $z + u = z$ for every $z \in \mathbb{C}$. Take $z = 0$: $0 + u = 0$. But $0 + u = u$ (Definition 3.10: $(0, 0) + (u_1, u_2) = (u_1, u_2)$). So $u = 0$. ∎ (Equivalently, with two candidate identities $u$ and $u'$: $u = u + u' = u' + u = u'$.)

**B6.** (a) The index set is $\{(i, j) : 1 \le i \le j \le n\}$. Reading it with $j$ outer: for fixed $j$, $i$ runs from $1$ to $j$. So by the triangular swap (the mirror image of (3.6)),
$$\sum_{i=1}^{n}\sum_{j=i}^{n}a_j = \sum_{j=1}^{n}\sum_{i=1}^{j}a_j = \sum_{j=1}^{n}a_j\sum_{i=1}^{j}1 = \sum_{j=1}^{n}j\,a_j,$$
since $a_j$ does not depend on $i$ and factors out of the inner sum by linearity (3.4), and $\sum_{i=1}^{j}1 = j$. ∎

(b) $\sum_i\sum_j\delta_{ij}a_ib_j = \sum_i a_i\Big(\sum_j\delta_{ij}b_j\Big)$ (linearity: $a_i$ is constant with respect to $j$) $= \sum_i a_i b_i$ by sifting (3.8). ∎

(c) By (3.7) with $b_j = a_j$: $\big(\sum_i a_i\big)^2 = \sum_{i=1}^{n}\sum_{j=1}^{n}a_ia_j$. Split the $n^2$ index pairs into three disjoint sets: $i = j$, $i < j$, $i > j$:
$$\sum_i\sum_j a_ia_j = \sum_{i} a_i^2 + \sum_{i<j}a_ia_j + \sum_{i>j}a_ia_j.$$
In the last sum rename the bound variables $i \leftrightarrow j$: $\sum_{i>j}a_ia_j = \sum_{j>i}a_ja_i = \sum_{i<j}a_ia_j$ (multiplication commutes). So the two off-diagonal sums are equal and the total is $\sum_i a_i^2 + 2\sum_{i<j}a_ia_j$. ∎

(d) $\big((AB)\mathbf x\big)_i = \sum_k (AB)_{ik}x_k = \sum_k\Big(\sum_j a_{ij}b_{jk}\Big)x_k = \sum_k\sum_j a_{ij}b_{jk}x_k$ (linearity: $x_k$ is constant with respect to $j$, so it distributes into the inner sum) $= \sum_j\sum_k a_{ij}b_{jk}x_k$ (rectangular swap (3.5)) $= \sum_j a_{ij}\Big(\sum_k b_{jk}x_k\Big)$ (linearity: $a_{ij}$ is constant with respect to $k$) $= \sum_j a_{ij}(B\mathbf x)_j = \big(A(B\mathbf x)\big)_i$. ∎

**B7.** (a) $|zw|^2 = (zw)\overline{(zw)} = zw\,\bar z\,\bar w$ by (3.11) $= (z\bar z)(w\bar w)$ (commutativity) $= |z|^2|w|^2$. Both $|zw|$ and $|z||w|$ are nonnegative reals with the same square, so they are equal. ∎

(b) $-16 = 16e^{i\pi}$, so $s = 16$, $\varphi = \pi$, $n = 4$ in (3.16): $z_k = 16^{1/4}e^{i(\pi + 2\pi k)/4} = 2e^{i\pi/4 + i\pi k/2}$ for $k = 0, 1, 2, 3$. Since $2e^{i\pi/4} = 2(\frac{\sqrt2}{2} + \frac{\sqrt2}{2}i) = \sqrt2 + \sqrt2\,i$, and each further root is the previous one rotated by $90°$ (multiplied by $i$):
$$\sqrt2 + \sqrt2\,i,\qquad -\sqrt2 + \sqrt2\,i,\qquad -\sqrt2 - \sqrt2\,i,\qquad \sqrt2 - \sqrt2\,i.$$
Check: $(\sqrt2 + \sqrt2\,i)^2 = 2(1 + i)^2 = 2 \cdot 2i = 4i$, and $(4i)^2 = -16$. ✓

(c) Let $\omega_j = e^{2\pi ij/n}$ and $\omega_k = e^{2\pi ik/n}$ be roots of unity. Then $(\omega_j\omega_k)^n = \omega_j^n\omega_k^n = 1 \cdot 1 = 1$, so the product is an $n$-th root of unity (explicitly, $\omega_j\omega_k = e^{2\pi i(j+k)/n} = \omega_{(j+k) \bmod n}$). And $(\omega_k^{-1})^n = (\omega_k^n)^{-1} = 1$, so the reciprocal is a root of unity too (explicitly $\omega_k^{-1} = e^{-2\pi ik/n} = \omega_{n-k} = \bar\omega_k$). ∎

(d) Let $S = \sum_{k=0}^{n-1}(\omega^m)^k$, a geometric series with ratio $\rho = \omega^m = e^{2\pi im/n}$. *Case $n \mid m$:* then $m = n q$ and $\rho = e^{2\pi iq} = 1$, so every term is $1$ and $S = n$. *Case $n \nmid m$:* then $2\pi m/n$ is not an integer multiple of $2\pi$, so $\rho \ne 1$ (last corollary of Theorem 3.13). Multiply $S$ by $\rho$ and shift the index as in Theorem 3.18: $\rho S = \sum_{k=1}^{n}\rho^k = S - 1 + \rho^n$, and $\rho^n = \omega^{mn} = (\omega^n)^m = 1$. So $\rho S = S$, i.e., $(\rho - 1)S = 0$, and since $\rho \ne 1$, $S = 0$. ∎ (Project Part 3 checks this for $n = 6$ and $-12 \le m \le 12$.)

**B8.** (a) Let $x(t) = Xe^{i\omega t}$ with $X \in \mathbb{C}$ constant. By (3.19) with $\lambda = i\omega$: $\dot x = i\omega Xe^{i\omega t}$ and $\ddot x = (i\omega)^2Xe^{i\omega t} = -\omega^2Xe^{i\omega t}$. Substitute into $m\ddot x + c\dot x + kx = F_0e^{i\omega t}$:
$$\big(-m\omega^2 + ic\omega + k\big)Xe^{i\omega t} = F_0e^{i\omega t}.$$
Since $e^{i\omega t} \ne 0$, divide it out: $(k - m\omega^2 + ic\omega)X = F_0$, so $X = \dfrac{F_0}{k - m\omega^2 + ic\omega}$. (The denominator is nonzero whenever $c > 0$ and $\omega > 0$, because its imaginary part $c\omega$ is then nonzero.)

Real parts: write $L[x] = m\ddot x + c\dot x + kx$. Because $m, c, k$ are real and differentiation acts separately on real and imaginary parts, $\operatorname{Re}L[x] = L[\operatorname{Re}x]$ (Theorem 3.21). Taking real parts of $L[Xe^{i\omega t}] = F_0e^{i\omega t}$ gives $L[\operatorname{Re}(Xe^{i\omega t})] = \operatorname{Re}(F_0e^{i\omega t}) = F_0\cos\omega t$. So $x_p = \operatorname{Re}[Xe^{i\omega t}]$ solves the real equation. ∎

(b) $k - m\omega^2 = 4 - 1 \cdot 2.25 = 1.75$ and $c\omega = 0.5 \cdot 1.5 = 0.75$. So
$$X = \frac{1}{1.75 + 0.75i} = \frac{1.75 - 0.75i}{1.75^2 + 0.75^2} = \frac{1.75 - 0.75i}{3.0625 + 0.5625} = \frac{1.75 - 0.75i}{3.625} = 0.4828 - 0.2069i.$$
$|X| = 1/\sqrt{3.625} = 1/1.9039 = 0.525$. $\arg X = -\operatorname{atan2}(0.75, 1.75) = -0.405$ rad ($-23.2°$). By (3.18), $x_p(t) = 0.525\cos(1.5t - 0.405)$. The phase is negative, so the displacement **lags** the force by $0.405$ rad: the mass reaches its maximum $0.405/1.5 = 0.27$ time units after the force does. (Project Part 9 prints $|X| = 0.5252$, $\arg X = -0.4049$ and checks that $x_p$ satisfies the ODE on a time grid.)

(c) Constant coefficients make $e^{i\omega t}$ an eigenfunction of the whole operator $L$ (each derivative only multiplies it by $i\omega$), so the ODE collapses to one algebraic equation for $X$; if the coefficients depended on $t$, $L[e^{i\omega t}]$ would no longer be a multiple of $e^{i\omega t}$. Linearity is what lets us take real parts at the end: $\operatorname{Re}L[x] = L[\operatorname{Re}x]$ fails for a nonlinear term such as $x^2$, since $\operatorname{Re}(x^2) \ne (\operatorname{Re}x)^2$.

## C. Challenge

**C1.** (a) Let $a \in \mathbb{R}$ and $\varepsilon > 0$ be given. Choose $\delta = \varepsilon/2$. If $|x - a| < \delta$, then $|f(x) - f(a)| = |(2x + 1) - (2a + 1)| = 2|x - a| < 2\delta = \varepsilon$. So $f$ is continuous at $a$. ∎ (Here $\delta$ does not depend on $a$: this $f$ is in fact uniformly continuous.)

(b) Let $a \in \mathbb{R}$ and $\varepsilon > 0$. We have $|x^2 - a^2| = |x - a|\,|x + a|$. If $|x - a| < 1$, then by the triangle inequality $|x + a| = |(x - a) + 2a| \le |x - a| + 2|a| < 1 + 2|a|$. Choose
$$\delta = \min\Big(1,\ \frac{\varepsilon}{1 + 2|a|}\Big) > 0.$$
If $|x - a| < \delta$, then $|x - a| < 1$ (so the bound on $|x + a|$ applies) and
$$|x^2 - a^2| = |x - a|\,|x + a| < \delta\,(1 + 2|a|) \le \frac{\varepsilon}{1 + 2|a|}(1 + 2|a|) = \varepsilon.$$
So $f$ is continuous at $a$. ∎ The $\delta$ depends on $a$: it shrinks like $1/(1 + 2|a|)$ as $|a|$ grows. That is the hint that no single $\delta$ will work for all $a$.

(c) Uniform continuity on $\mathbb{R}$ is $\forall\varepsilon > 0\ \exists\delta > 0\ \forall a\ \forall x\ \big(|x - a| < \delta \Rightarrow |x^2 - a^2| < \varepsilon\big)$. Negating with (3.1) and (3.3):
$$\exists\varepsilon > 0\ \ \forall\delta > 0\ \ \exists a\ \exists x\ \ \big(|x - a| < \delta\ \wedge\ |x^2 - a^2| \ge \varepsilon\big).$$
We win this game. Take $\varepsilon = 1$. Let $\delta > 0$ be arbitrary (the opponent's move). Choose $a = 1/\delta$ and $x = a + \delta/2$. Then $|x - a| = \delta/2 < \delta$, and
$$x^2 - a^2 = (x - a)(x + a) = \frac{\delta}{2}\Big(\frac{2}{\delta} + \frac{\delta}{2}\Big) = 1 + \frac{\delta^2}{4} \ge 1 = \varepsilon.$$
So for every $\delta$ there is a pair within $\delta$ whose squares differ by at least $1$: $x^2$ is not uniformly continuous on $\mathbb{R}$. ∎

Cauchy's error had the same shape, with the number of terms $N$ of his series in place of $\delta$: pointwise convergence gave him, for each point $x$ and tolerance $\varepsilon$, an $N$ ($\forall x\,\forall\varepsilon\,\exists N$), and his argument used one $N$ for all points at once ($\forall\varepsilon\,\exists N\,\forall x$). Part (b) versus part (c) is that gap, for continuity instead of convergence, on one function.

**C2.** (a) *Strong induction on $n$.* Base cases: $12 = 4 \cdot 3 + 5 \cdot 0$, $13 = 4 \cdot 2 + 5 \cdot 1$, $14 = 4 \cdot 1 + 5 \cdot 2$, $15 = 4 \cdot 0 + 5 \cdot 3$. Inductive step: let $n \ge 16$ and assume every integer $m$ with $12 \le m < n$ is of the form $4a + 5b$ with $a, b \ge 0$. Since $n - 4 \ge 12$ and $n - 4 < n$, the hypothesis gives $n - 4 = 4a + 5b$, so $n = 4(a + 1) + 5b$ with $a + 1, b \ge 0$. ∎ Four base cases are needed because the step reaches back four.

(b) *Strong induction.* Base cases: $F_1 = 1 < 2 = 2^1$ and $F_2 = 1 < 4 = 2^2$. Step: let $n \ge 3$ and assume $F_m < 2^m$ for all $1 \le m < n$. Then
$$F_n = F_{n-1} + F_{n-2} < 2^{n-1} + 2^{n-2} < 2^{n-1} + 2^{n-1} = 2^n,$$
using the hypothesis for $m = n - 1$ and $m = n - 2$. ∎ Ordinary induction would only supply $F_{n-1} < 2^{n-1}$; the recurrence also needs $F_{n-2}$, so you need the claim for two previous values, which is what strong induction (or an induction on the pair $(F_{n-1}, F_n)$) provides. Two base cases are needed because the step reaches back two.

(c) Add to `project.py` (for example at the end of Part 9):

```python
def is_4a_plus_5b(n):
    """True iff n = 4a + 5b for some integers a, b >= 0 (exercise C2a)."""
    return any((n - 5 * b) >= 0 and (n - 5 * b) % 4 == 0 for b in range(n // 5 + 1))

assert all(is_4a_plus_5b(n) for n in range(12, 201))
assert not is_4a_plus_5b(11)   # 11 = 4a + 5b has no solution: b=0 gives 11, b=1 gives 6, b=2 gives 1
```

**C3.** (a) By Theorem 3.22, $|z| = |(z - w) + w| \le |z - w| + |w|$, so $|z| - |w| \le |z - w|$. Swapping the roles of $z$ and $w$: $|w| - |z| \le |w - z| = |z - w|$. A real number and its negative are both $\le |z - w|$, so $\big||z| - |w|\big| \le |z - w|$. ∎

(b) *Equidistant.* $|\omega^{k+1} - \omega^k| = |\omega^k(\omega - 1)| = |\omega^k|\,|\omega - 1| = |\omega - 1|$ by Theorem 3.14, since $|\omega^k| = 1$. This does not depend on $k$, so all $n$ consecutive distances (including from $\omega^{n-1}$ back to $\omega^0 = \omega^n$) are equal. All vertices lie on the unit circle (Theorem 3.17). A closed polygon inscribed in a circle with all sides equal is regular. For the side length: $\omega - 1 = e^{i\pi/n}(e^{i\pi/n} - e^{-i\pi/n}) = e^{i\pi/n}\cdot 2i\sin(\pi/n)$ by (3.17), so $|\omega - 1| = 2\sin(\pi/n)$.

*Area.* Join the center to the vertices: the polygon is $n$ congruent isosceles triangles with two sides of length $1$ (radii) and apex angle $2\pi/n$. The area of a triangle with sides $1, 1$ and included angle $\alpha$ is $\frac12 \cdot 1 \cdot 1 \cdot \sin\alpha$. So the total area is $n \cdot \frac12\sin\frac{2\pi}{n} = \frac n2\sin\frac{2\pi}{n}$. ∎ Values: $n = 3$: $1.2990$; $n = 4$: $2$; $n = 6$: $2.5981$; $n = 8$: $2.8284$.

Code for the "Your turn" (shoelace formula on the vertices in order):

```python
for n_ in range(3, 9):
    w = roots_of_unity(n_)
    x, y = w.real, w.imag
    area = 0.5 * abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
    assert np.isclose(area, n_ / 2 * np.sin(2 * np.pi / n_))
```

(c) $\dfrac n2\sin\dfrac{2\pi}{n} = \pi\cdot\dfrac{\sin(2\pi/n)}{2\pi/n} \to \pi \cdot 1 = \pi$ as $n \to \infty$, since $\sin u/u \to 1$ as $u \to 0$. That is the area of the unit disk: the inscribed polygons fill the circle. (Archimedes' method.)

**C4.** (a) *Claim.* Every real polynomial $p$ of degree $n \ge 1$ with leading coefficient $a_n$ can be written as $a_n$ times a product of factors of the form $(z - r)$ with $r \in \mathbb{R}$ and $(z^2 - 2\operatorname{Re}(r)z + |r|^2)$ with $r \notin \mathbb{R}$. *Proof by strong induction on $n$.*

Base $n = 1$: $p = a_1 z + a_0 = a_1(z - r)$ with $r = -a_0/a_1 \in \mathbb{R}$.

Step: let $n \ge 2$ and assume the claim for all real polynomials of degree $< n$. By the FTA (Theorem 3.23), $p$ has a root $r \in \mathbb{C}$.

*Case 1: $r \in \mathbb{R}$.* By the factor theorem (step 1 in the derivation after Theorem 3.23), $p = (z - r)q$. Dividing the real polynomial $p$ by the real polynomial $z - r$ produces a real quotient $q$ of degree $n - 1$ with leading coefficient $a_n$. Apply the hypothesis to $q$.

*Case 2: $r \notin \mathbb{R}$.* By Theorem 3.24, $\bar r$ is also a root, and $\bar r \ne r$. Factor: $p = (z - r)q_1$. Then $0 = p(\bar r) = (\bar r - r)q_1(\bar r)$ and $\bar r - r \ne 0$, so $q_1(\bar r) = 0$ and $q_1 = (z - \bar r)q_2$. Hence $p = (z - r)(z - \bar r)q_2 = (z^2 - 2\operatorname{Re}(r)z + |r|^2)\,q_2$, where $(z - r)(z - \bar r) = z^2 - (r + \bar r)z + r\bar r$ and $r + \bar r = 2\operatorname{Re} r$, $r\bar r = |r|^2$ by (3.11). The quadratic has real coefficients, so dividing $p$ by it gives a real $q_2$ of degree $n - 2$ with leading coefficient $a_n$. If $n = 2$, $q_2 = a_2$ is a constant and we are done; otherwise apply the hypothesis to $q_2$. ∎

(b) In the factorization of (a), $\deg p = (\text{number of linear factors}) + 2\,(\text{number of quadratic factors})$. If $\deg p$ is odd, the number of linear factors is odd, hence at least $1$, and each linear factor $(z - r)$ with $r \in \mathbb{R}$ gives a real root $r$. ∎

(c) For a real $n \times n$ matrix $A$, the characteristic polynomial $\chi(\lambda) = \det(A - \lambda I)$ has real coefficients and degree $n$. If $n$ is odd, (b) gives a real root $\lambda_0$, so $\det(A - \lambda_0 I) = 0$ and $A$ has the real eigenvalue $\lambda_0$ (LA-7: a real eigenvector then exists by solving the real system $(A - \lambda_0 I)\mathbf v = \mathbf 0$). Example: a rotation of $\mathbb{R}^3$ about an axis. Its characteristic polynomial is cubic, so it has a real eigenvalue; geometrically that eigenvalue is $1$ and the eigenvector is the axis, the one direction the rotation leaves fixed. A rotation of $\mathbb{R}^2$ by $90°$ (Theorem 3.15 with $z = i$) has characteristic polynomial $\lambda^2 + 1$ and no real eigenvalue, which is allowed because $n = 2$ is even.
