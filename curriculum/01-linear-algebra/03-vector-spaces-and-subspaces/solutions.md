# LA-3 Solutions

## A. Warm-up

**A1.**
(a) Not a subspace. $(1, 0)$ is in the set but $(-1)(1, 0) = (-1, 0)$ has $x < 0$. (Closed under addition, not under scaling.)
(b) Not a subspace. $\mathbf{0} = (0, 0)$ is not on the line ($0 \ne 2\cdot 0 + 1$). Every subspace contains $\mathbf{0}$ (Theorem 3.4).
(c) Not a subspace. The set is the union of the two axes. $(1, 0)$ and $(0, 1)$ are in it but $(1, 0) + (0, 1) = (1, 1)$ has $xy = 1 \ne 0$. (Closed under scaling, not under addition; compare B1(b).)
(d) Subspace: it is $N(\mathbf{n}^\top)$ for the $1 \times 3$ matrix $\mathbf{n}^\top = (1, 2, -1)$, and null spaces are subspaces (Theorem 3.6). As a span: special solutions $(-2, 1, 0)^\top$ and $(1, 0, 1)^\top$.
(e) Subspace. Nonempty: the zero polynomial has $p(0) = p(1) = 0$. If $p(0) = p(1)$ and $q(0) = q(1)$, then $(p + q)(0) = p(0) + q(0) = p(1) + q(1) = (p + q)(1)$ and $(cp)(0) = c\,p(0) = c\,p(1) = (cp)(1)$. In coefficient coordinates it is the null space of the $1 \times 4$ matrix $(0, 1, 1, 1)$, since $p(1) - p(0) = a_1 + a_2 + a_3$.
(f) Not a subspace. $\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}$ and $\begin{pmatrix}0 & 0\\0 & 1\end{pmatrix}$ both have determinant $0$; their sum is $I$ with determinant $1$.

**A2.** (a) Pivot columns $1, 3$; free columns $2, 4$; $r = 2$.
(b) Row 1: $x_1 = -3x_2 - 2x_4$. Row 2: $x_3 = x_4$. With $x_2 = 1, x_4 = 0$: $\mathbf{s}_1 = (-3, 1, 0, 0)^\top$. With $x_2 = 0, x_4 = 1$: $\mathbf{s}_2 = (-2, 0, 1, 1)^\top$. Check $A\mathbf{s}_2$: row 1: $-2 + 0 + 0 + 2 = 0$; row 2: $0 + 0 + 1 - 1 = 0$. ✓ Basis of $N(A)$: $(\mathbf{s}_1, \mathbf{s}_2)$; $\dim N(A) = 2 = 4 - 2$.
(c) $A$ is already $R$ and $E = I$, so set $x_2 = x_4 = 0$: $x_1 = 4$, $x_3 = 1$. $\mathbf{x}_p = (4, 0, 1, 0)^\top$. Complete solution:
$$\mathbf{x} = (4, 0, 1, 0)^\top + c_1(-3, 1, 0, 0)^\top + c_2(-2, 0, 1, 1)^\top.$$
(d) $\dim C(A) = 2$ (so $C(A) = \mathbb{R}^2$: both pivot columns $(1, 0)^\top, (0, 1)^\top$), $\dim N(A) = 2$, $\dim C(A^\top) = 2$, $\dim N(A^\top) = m - r = 0$. Every $\mathbf{b} \in \mathbb{R}^2$ is reachable.

**A3.** (a) Columns of $A = \begin{pmatrix} 1 & 2 & 3 \\ 2 & 3 & 4 \\ 3 & 4 & 5 \end{pmatrix}$. $R_2 - 2R_1 = (0, -1, -2)$, $R_3 - 3R_1 = (0, -2, -4)$, $R_3 - 2R_2 = \mathbf{0}$, scale $R_2$ by $-1$: $(0, 1, 2)$, $R_1 - 2R_2 = (1, 0, -1)$. So $R = \begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \\ 0 & 0 & 0\end{pmatrix}$, $r = 2 < 3$: dependent. Special solution $x_3 = 1$, $x_1 = 1$, $x_2 = -2$: relation $\mathbf{v}_1 - 2\mathbf{v}_2 + \mathbf{v}_3 = \mathbf{0}$. Check: $(1 - 4 + 3, 2 - 6 + 4, 3 - 8 + 5) = \mathbf{0}$. ✓
(b) Two independent directions (columns 1, 2): a plane through $\mathbf{0}$. A normal $\mathbf{n}$ needs $\mathbf{n}\cdot\mathbf{v}_1 = \mathbf{n}\cdot\mathbf{v}_2 = 0$: $n_1 + 2n_2 + 3n_3 = 0$ and $2n_1 + 3n_2 + 4n_3 = 0$. Subtract twice the first from the second: $-n_2 - 2n_3 = 0$, so $n_2 = -2n_3$, then $n_1 = -2n_2 - 3n_3 = n_3$. Take $\mathbf{n} = (1, -2, 1)^\top$: the plane is $x - 2y + z = 0$. (Check $\mathbf{v}_3$: $3 - 8 + 5 = 0$. ✓ That $\mathbf{n}$ equals the relation vector of (a) is a coincidence of this matrix, whose rows and columns are both arithmetic progressions, so $N(A) = N(A^\top)$.)
(c) $0 - 2 + 2 = 0$, so yes. Solve $a\mathbf{v}_1 + b\mathbf{v}_2 = (0, 1, 2)^\top$: $a + 2b = 0$, $2a + 3b = 1$. From the first $a = -2b$; then $-4b + 3b = 1$, $b = -1$, $a = 2$. Check third entry: $2\cdot 3 - 4 = 2$. ✓ So $(0, 1, 2)^\top = 2\mathbf{v}_1 - \mathbf{v}_2$.

**A4.** (a) $a(1, 0, 1) + b(0, 1, 1) = (1, 1, 1)$ needs $a = 1$, $b = 1$, $a + b = 1$. The third fails ($2 \ne 1$). Not in the span.
(b) $\mathbf{u}_1, \mathbf{u}_2$ are independent (neither is a multiple of the other: $\mathbf{u}_1$ has a $0$ in slot 2 and a $1$ in slot 1, $\mathbf{u}_2$ the reverse, so $c\mathbf{u}_1 = \mathbf{u}_2$ forces $c = 0$ and $c = 1$). By (a) and Lemma 3.20, $\mathbf{u}_1, \mathbf{u}_2, \mathbf{u}_3$ is independent. Three independent vectors in the 3-dimensional $\mathbb{R}^3$ form a basis (Corollary 3.26(a)).
(c) $a\mathbf{u}_1 + b\mathbf{u}_2 + c\mathbf{u}_3 = (2, 3, 4)^\top$: $a + c = 2$, $b + c = 3$, $a + b + c = 4$. From the first two, $a = 2 - c$, $b = 3 - c$; the third gives $5 - 2c + c = 4$, so $c = 1$, $a = 1$, $b = 2$. Check: $(1, 0, 1) + 2(0, 1, 1) + (1, 1, 1) = (2, 3, 4)$. ✓ Coordinates $(1, 2, 1)^\top$. Exactly one answer because the list is a basis: independence is uniqueness of coefficients (Theorem 3.15 / 3.17).

**A5.** Reduce $[A \mid I]$:
$$
\left(\begin{array}{ccc|ccc} 1 & 1 & 2 & 1 & 0 & 0 \\ 2 & 1 & 3 & 0 & 1 & 0 \\ 3 & 0 & 3 & 0 & 0 & 1 \end{array}\right)
\xrightarrow{R_2 - 2R_1,\ R_3 - 3R_1}
\left(\begin{array}{ccc|ccc} 1 & 1 & 2 & 1 & 0 & 0 \\ 0 & -1 & -1 & -2 & 1 & 0 \\ 0 & -3 & -3 & -3 & 0 & 1 \end{array}\right)
\xrightarrow{R_3 - 3R_2}
\left(\begin{array}{ccc|ccc} 1 & 1 & 2 & 1 & 0 & 0 \\ 0 & -1 & -1 & -2 & 1 & 0 \\ 0 & 0 & 0 & 3 & -3 & 1 \end{array}\right)
$$
$$
\xrightarrow{R_2 \leftarrow -R_2}
\left(\begin{array}{ccc|ccc} 1 & 1 & 2 & 1 & 0 & 0 \\ 0 & 1 & 1 & 2 & -1 & 0 \\ 0 & 0 & 0 & 3 & -3 & 1 \end{array}\right)
\xrightarrow{R_1 - R_2}
\left(\begin{array}{ccc|ccc} 1 & 0 & 1 & -1 & 1 & 0 \\ 0 & 1 & 1 & 2 & -1 & 0 \\ 0 & 0 & 0 & 3 & -3 & 1 \end{array}\right).
$$
So $R = \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 0 \end{pmatrix}$, $E = \begin{pmatrix} -1 & 1 & 0 \\ 2 & -1 & 0 \\ 3 & -3 & 1 \end{pmatrix}$, $r = 2$, pivots $1, 2$, free $3$. Check row 3 of $EA$: $3(1, 1, 2) - 3(2, 1, 3) + (3, 0, 3) = (0, 0, 0)$. ✓

- $C(A)$: columns 1, 2 of $A$: $(1, 2, 3)^\top, (1, 1, 0)^\top$. $\dim 2$. (Column 3 $=$ column 1 $+$ column 2, as $R$'s third column $(1, 1, 0)^\top$ says.)
- $N(A)$: $x_3 = 1$, $x_1 = -1$, $x_2 = -1$: $(-1, -1, 1)^\top$. $\dim 1$. Check: $-(1, 2, 3) - (1, 1, 0) + (2, 3, 3) = \mathbf{0}$. ✓
- $C(A^\top)$: rows of $R$: $(1, 0, 1)^\top, (0, 1, 1)^\top$. $\dim 2$.
- $N(A^\top)$: row 3 of $E$: $(3, -3, 1)^\top$. $\dim 1$. Check $\mathbf{y}^\top A = (3 - 6 + 3,\ 3 - 3 + 0,\ 6 - 9 + 3) = (0, 0, 0)$. ✓

Orthogonality (3.6): $(1, 0, 1)\cdot(-1, -1, 1) = 0$ and $(0, 1, 1)\cdot(-1, -1, 1) = 0$; $(3, -3, 1)\cdot(1, 2, 3) = 0$ and $(3, -3, 1)\cdot(1, 1, 0) = 0$. ✓

```
          R^3                                        R^3
   +-----------------------+                 +-----------------------+
   | C(A^T)  dim 2         |        A        | C(A)    dim 2         |
   | (1,0,1), (0,1,1)      | --------------> | (1,2,3), (1,1,0)      |
   + - - - - -_|_- - - - - +                 + - - - - -_|_- - - - - +
   | N(A)    dim 1         |     A x = 0     | N(A^T)  dim 1         |
   | (-1,-1,1)             | --------------> | (3,-3,1)              |
   +-----------------------+                 +-----------------------+
          2 + 1 = 3                                 2 + 1 = 3
```

## B. Core

**B1.** (a) *Subspace.* Nonempty: $\mathbf{0} = \mathbf{0} + \mathbf{0} \in U + W$. Let $\mathbf{x} = \mathbf{u} + \mathbf{w}$ and $\mathbf{x}' = \mathbf{u}' + \mathbf{w}'$ with $\mathbf{u}, \mathbf{u}' \in U$, $\mathbf{w}, \mathbf{w}' \in W$, and $c \in \mathbb{R}$. Then $\mathbf{x} + \mathbf{x}' = (\mathbf{u} + \mathbf{u}') + (\mathbf{w} + \mathbf{w}')$ by (A1), (A2), with $\mathbf{u} + \mathbf{u}' \in U$ and $\mathbf{w} + \mathbf{w}' \in W$ by closure of each subspace; so $\mathbf{x} + \mathbf{x}' \in U + W$. And $c\mathbf{x} = c\mathbf{u} + c\mathbf{w}$ by (M1), with $c\mathbf{u} \in U$, $c\mathbf{w} \in W$; so $c\mathbf{x} \in U + W$. By Theorem 3.4, $U + W$ is a subspace.
*Smallest.* $U \subseteq U + W$ (take $\mathbf{w} = \mathbf{0}$) and $W \subseteq U + W$ (take $\mathbf{u} = \mathbf{0}$), so $U \cup W \subseteq U + W$. If $X$ is any subspace containing $U \cup W$, then for $\mathbf{u} \in U$, $\mathbf{w} \in W$ both are in $X$, so $\mathbf{u} + \mathbf{w} \in X$ by closure. Hence $U + W \subseteq X$. ∎

(b) $(\Leftarrow)$ If $U \subseteq W$ then $U \cup W = W$, a subspace; likewise if $W \subseteq U$.
$(\Rightarrow)$ Suppose $U \cup W$ is a subspace and, for contradiction, neither contains the other. Pick $\mathbf{u} \in U$ with $\mathbf{u} \notin W$, and $\mathbf{w} \in W$ with $\mathbf{w} \notin U$. Both lie in $U \cup W$, so $\mathbf{u} + \mathbf{w} \in U \cup W$. Case 1: $\mathbf{u} + \mathbf{w} \in U$. Then $\mathbf{w} = (\mathbf{u} + \mathbf{w}) - \mathbf{u} \in U$ (closure of $U$ under addition and scaling by $-1$), contradicting the choice of $\mathbf{w}$. Case 2: $\mathbf{u} + \mathbf{w} \in W$. Then $\mathbf{u} = (\mathbf{u} + \mathbf{w}) - \mathbf{w} \in W$, contradicting the choice of $\mathbf{u}$. Both cases are impossible, so one subspace contains the other. ∎ (A1(c) is the picture: two lines, and the sum of a vector from each leaves the union.)

**B2.** (a) $c\mathbf{0} = c(\mathbf{0} + \mathbf{0}) = c\mathbf{0} + c\mathbf{0}$ by (A3) then (M1). Add $-(c\mathbf{0})$ to both sides, using (A2), (A4), (A3) exactly as in Theorem 3.2: $\mathbf{0} = c\mathbf{0}$. ∎
(b) $\mathbf{v} + (-1)\mathbf{v} = 1\mathbf{v} + (-1)\mathbf{v} = (1 + (-1))\mathbf{v} = 0\mathbf{v} = \mathbf{0}$ by (M4), (M2), arithmetic, Theorem 3.2. So $(-1)\mathbf{v}$ is a negative of $\mathbf{v}$. Negatives are unique: if $\mathbf{v} + \mathbf{w} = \mathbf{0}$ and $\mathbf{v} + \mathbf{w}' = \mathbf{0}$, then $\mathbf{w} = \mathbf{w} + \mathbf{0} = \mathbf{w} + (\mathbf{v} + \mathbf{w}') = (\mathbf{w} + \mathbf{v}) + \mathbf{w}' = (\mathbf{v} + \mathbf{w}) + \mathbf{w}' = \mathbf{0} + \mathbf{w}' = \mathbf{w}' + \mathbf{0} = \mathbf{w}'$, using (A3), hypothesis, (A2), (A1), hypothesis, (A1), (A3). Hence $(-1)\mathbf{v} = -\mathbf{v}$. ∎
(c) Suppose $\mathbf{0}$ and $\mathbf{0}'$ both satisfy (A3). Then $\mathbf{0} = \mathbf{0} + \mathbf{0}'$ (using $\mathbf{0}'$ as the zero) $= \mathbf{0}' + \mathbf{0}$ (A1) $= \mathbf{0}'$ (using $\mathbf{0}$ as the zero). ∎
(d) Suppose $c\mathbf{v} = \mathbf{0}$ and $c \ne 0$. Then $\mathbf{v} = 1\mathbf{v} = (c^{-1}c)\mathbf{v} = c^{-1}(c\mathbf{v}) = c^{-1}\mathbf{0} = \mathbf{0}$ by (M4), arithmetic, (M3), hypothesis, part (a). So either $c = 0$ or $\mathbf{v} = \mathbf{0}$. ∎

**B3.** (a) $0^\top = 0$, so $0 \in \mathrm{Sym}_n$. If $A^\top = A$ and $B^\top = B$: $(A + B)^\top = A^\top + B^\top = A + B$ and $(cA)^\top = cA^\top = cA$ (transpose is entrywise: $(A + B)^\top_{ij} = (A + B)_{ji} = A_{ji} + B_{ji}$). Theorem 3.4. ∎
(b) Let $E_{ij}$ be the matrix with a $1$ in position $(i, j)$. Basis: $\mathcal{B} = (E_{11},\ E_{12} + E_{21},\ E_{22})$. Every symmetric $2 \times 2$ matrix is $\begin{pmatrix} a & b \\ b & c \end{pmatrix} = aE_{11} + b(E_{12} + E_{21}) + cE_{22}$, and the coefficients $(a, b, c)$ are forced: $a$ is the $(1,1)$ entry, $b$ the $(1,2)$ entry, $c$ the $(2,2)$ entry. Existence and uniqueness of coefficients for every element is Theorem 3.17's criterion for a basis. $\dim \mathrm{Sym}_2 = 3$.
(c) $\mathrm{Sym}_n$: the list of $E_{ii}$ ($i = 1..n$) and $E_{ij} + E_{ji}$ ($i < j$). A symmetric $A$ equals $\sum_i a_{ii}E_{ii} + \sum_{i<j} a_{ij}(E_{ij} + E_{ji})$ because the right side has $(i, i)$ entry $a_{ii}$, $(i, j)$ entry $a_{ij}$ and $(j, i)$ entry $a_{ij} = a_{ji}$; and the coefficients are forced by the entries, so the representation is unique. Count: $n + \binom{n}{2} = n + n(n-1)/2 = n(n+1)/2$.
$\mathrm{Skew}_n$: $A^\top = -A$ forces $a_{ii} = -a_{ii}$, so the diagonal is zero. The list $E_{ij} - E_{ji}$ ($i < j$) works the same way: $A = \sum_{i<j} a_{ij}(E_{ij} - E_{ji})$ uniquely. Count $\binom{n}{2} = n(n-1)/2$.
(d) *Existence.* $S = \tfrac12(A + A^\top)$, $K = \tfrac12(A - A^\top)$. Then $S + K = A$, $S^\top = \tfrac12(A^\top + A) = S$, $K^\top = \tfrac12(A^\top - A) = -K$. *Uniqueness.* If also $A = S' + K'$, then $S - S' = K' - K$. The left side is symmetric, the right side skew (differences stay in each subspace). So $X = S - S'$ satisfies $X^\top = X$ and $X^\top = -X$, hence $X = -X$, $2X = 0$, $X = 0$. So $S = S'$, $K = K'$. ∎ Dimensions: $n(n+1)/2 + n(n-1)/2 = n^2 = \dim\mathcal{M}_{n\times n}$. That is no accident. Existence and uniqueness of the split say exactly that $\mathcal{M}_{n\times n} = \mathrm{Sym}_n \oplus \mathrm{Skew}_n$ (Definition 3.33; the uniqueness step also showed $\mathrm{Sym}_n \cap \mathrm{Skew}_n = \{0\}$, the form of Theorem 3.34(b)), and Theorem 3.34(c) says the dimensions of the summands of a direct sum add to the whole dimension.

**B4.** (a) Nonempty: the zero polynomial vanishes at $1$. If $p(1) = q(1) = 0$ then $(p + q)(1) = 0 + 0 = 0$ and $(cp)(1) = c\cdot 0 = 0$. Theorem 3.4. In coordinates $p = a_0 + a_1x + a_2x^2 + a_3x^3$, $p(1) = a_0 + a_1 + a_2 + a_3$, so $W$ corresponds to $N\big((1, 1, 1, 1)\big)$. One pivot (column 1), free columns $2, 3, 4$, $a_0 = -a_1 - a_2 - a_3$. Special solutions $(-1, 1, 0, 0), (-1, 0, 1, 0), (-1, 0, 0, 1)$, i.e. the polynomials $x - 1,\ x^2 - 1,\ x^3 - 1$. Basis of $W$ by Theorem 3.18 (coordinates preserve independence and span, Theorem 3.17). $\dim W = 4 - 1 = 3$.
(b) The zero polynomial has $p(1) = 0 \ne 1$, so $\mathbf{0}$ is not in the set. (Also: $p = 1$ is in it but $2p$ is not.)
(c) Two equations: $a_0 + a_1 + a_2 + a_3 = 0$ and $a_0 + 2a_1 + 4a_2 + 8a_3 = 0$. Matrix $\begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & 2 & 4 & 8 \end{pmatrix}$; $R_2 - R_1 = (0, 1, 3, 7)$; $R_1 - R_2 = (1, 0, -2, -6)$. $R = \begin{pmatrix} 1 & 0 & -2 & -6 \\ 0 & 1 & 3 & 7 \end{pmatrix}$, $r = 2$, free columns $3, 4$. Special solutions: $a_2 = 1, a_3 = 0$: $a_0 = 2$, $a_1 = -3$, giving $2 - 3x + x^2 = (x - 1)(x - 2)$. $a_2 = 0, a_3 = 1$: $a_0 = 6$, $a_1 = -7$, giving $6 - 7x + x^3$; check $1 - 7 + 6 = 0$ and $8 - 14 + 6 = 0$ ✓; factor: $(x - 1)(x - 2)(x + 3)$ (expand: $(x^2 - 3x + 2)(x + 3) = x^3 - 7x + 6$ ✓). $\dim W_2 = 4 - 2 = 2$. Any polynomial in $W_2$ is $(x - 1)(x - 2)q(x)$ with $\deg q \le 1$, so $\big((x - 1)(x - 2),\ x(x - 1)(x - 2)\big)$ is another basis.

**B5.** Let $W$ be a subspace of $\mathbb{R}^3$. By Corollary 3.26(d), $W$ is finite-dimensional with $\dim W \le 3$, so $\dim W \in \{0, 1, 2, 3\}$. Let $\mathcal{B}$ be a basis of $W$ (Theorem 3.22).
- $\dim W = 0$: the basis is empty, $W = \operatorname{span}\emptyset = \{\mathbf{0}\}$.
- $\dim W = 1$: $\mathcal{B} = (\mathbf{d})$ with $\mathbf{d} \ne \mathbf{0}$ (Theorem 3.13(c)); $W = \{t\mathbf{d}\}$, a line through $\mathbf{0}$.
- $\dim W = 2$: $\mathcal{B} = (\mathbf{d}_1, \mathbf{d}_2)$ independent; $W = \{s\mathbf{d}_1 + t\mathbf{d}_2\}$, a plane through $\mathbf{0}$ (independent means not parallel, so the combinations fill a plane, not a line).
- $\dim W = 3$: $W = \mathbb{R}^3$ by Corollary 3.26(d). ∎

*Two descriptions of a plane.* The span description is the $\dim 2$ case above. For the null-space description, let $W = \operatorname{span}\{\mathbf{d}_1, \mathbf{d}_2\}$. Let $D$ be the $2 \times 3$ matrix with rows $\mathbf{d}_1^\top, \mathbf{d}_2^\top$; its rank is $2$ (independent rows, Theorem 3.29), so $\dim N(D) = 3 - 2 = 1$. Pick $\mathbf{n} \ne \mathbf{0}$ in $N(D)$: then $\mathbf{n}\cdot\mathbf{d}_1 = \mathbf{n}\cdot\mathbf{d}_2 = 0$, so every vector of $W$ is orthogonal to $\mathbf{n}$ (by (3.6) applied to $D$, or directly by linearity of the dot product): $W \subseteq N(\mathbf{n}^\top)$. Now $\mathbf{n}^\top$ is a $1 \times 3$ matrix of rank $1$, so $\dim N(\mathbf{n}^\top) = 3 - 1 = 2 = \dim W$. A subspace of a space of the same dimension equals it (Corollary 3.26(d)), so $W = N(\mathbf{n}^\top)$. ∎ Conversely every $N(\mathbf{n}^\top)$ with $\mathbf{n} \ne \mathbf{0}$ has dimension $2$ and is a span of its two special solutions.

**B6.** (a) Column $j$ of $C\tilde R$ is $\sum_{i=1}^r \tilde R_{ij}\,\mathbf{c}_i = \sum_{i=1}^r R_{ij}\,\mathbf{a}_{p_i}$ (LA-1 (3.1); the $i$-th column of $C$ is $\mathbf{a}_{p_i}$ and the $i$-th row of $\tilde R$ is row $i$ of $R$). If $j = p_k$ is a pivot column, $R_{i,p_k} = 1$ for $i = k$ and $0$ otherwise (RREF), so the sum is $\mathbf{a}_{p_k}$, which is column $j$ of $A$. If $j = f$ is a free column, Lemma 3.27 (via Theorem 3.28 step 3–4) gives $\mathbf{a}_f = \sum_i R_{i,f}\mathbf{a}_{p_i}$, again column $j$ of $A$. Every column agrees, so $A = C\tilde R$. ∎
(b) Example 4.2: $C = \begin{pmatrix} 1 & 0 \\ 2 & 1 \\ 1 & 1 \end{pmatrix}$, $\tilde R = \begin{pmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 2 \end{pmatrix}$. Row 1 of $C\tilde R$: $1\cdot(1, 2, 0, 1) + 0 = (1, 2, 0, 1)$. Row 2: $2(1, 2, 0, 1) + 1(0, 0, 1, 2) = (2, 4, 1, 4)$. Row 3: $(1, 2, 0, 1) + (0, 0, 1, 2) = (1, 2, 1, 3)$. That is $A$. ✓ (The project asserts `C42 @ Rt42 == A42` for this matrix in Part 3, and $A = C\tilde R$ on 300 random matrices in Part 2.)
(c) Row $i$ of $A = C\tilde R$ is $\sum_{k=1}^r C_{ik}(\text{row } k \text{ of } \tilde R)$: a combination of the $r$ rows of $\tilde R$. So $C(A^\top) \subseteq \operatorname{span}\{\text{rows of } \tilde R\}$, which has dimension at most $r$ (a spanning list of $r$ vectors contains a basis, Theorem 3.22). Hence $\dim C(A^\top) \le r = \dim C(A)$, where $r = \dim C(A)$ is Theorem 3.28. Apply the same inequality to $A^\top$: $\dim C(A) = \dim C((A^\top)^\top) \le \dim C(A^\top)$. Together: equality. ∎

**B7.** (a) $(\supseteq)$ If $A\mathbf{x} = \mathbf{0}$ then $A^\top A\mathbf{x} = A^\top\mathbf{0} = \mathbf{0}$. $(\subseteq)$ If $A^\top A\mathbf{x} = \mathbf{0}$ then $\mathbf{x}^\top A^\top A\mathbf{x} = 0$. But $\mathbf{x}^\top A^\top A\mathbf{x} = (A\mathbf{x})^\top(A\mathbf{x}) = \|A\mathbf{x}\|^2$ (LA-1 Theorem 3.8 and Definition 3.9). A vector of length $0$ is $\mathbf{0}$, so $A\mathbf{x} = \mathbf{0}$. ∎
(b) $A$ and $A^\top A$ both have $n$ columns and the same null space, so by Rank–Nullity (3.5), $\operatorname{rank}(A^\top A) = n - \dim N(A^\top A) = n - \dim N(A) = \operatorname{rank}(A)$. ∎
(c) A square matrix is invertible iff its null space is $\{\mathbf{0}\}$ (LA-2; equivalently iff it has $n$ pivots, Corollary 3.10). So $A^\top A$ invertible $\iff N(A^\top A) = \{\mathbf{0}\} \iff N(A) = \{\mathbf{0}\} \iff$ columns of $A$ independent (Theorem 3.12). ∎
(d) Each column of $A^\top A$ is $A^\top$ times a column of $A$, hence in $C(A^\top)$; so $C(A^\top A) \subseteq C(A^\top)$ (Theorem 3.5). By (b) and Theorem 3.29(c), $\dim C(A^\top A) = \operatorname{rank}(A) = \dim C(A^\top)$. A subspace of equal dimension is the whole space (Corollary 3.26(d)). ∎ (LA-5 uses (c) to invert $A^\top A$ in the normal equations and (d) to say the least-squares solution is unique in the row space.)

**B8.** Reference implementation (add after `extend_to_basis` in `project.py`):

```python
def basis_of_sum(U, W, tol=1e-10):
    """Basis of U + W: pivot columns of [U | W] (Theorem 3.28 applied to the spanning list)."""
    return column_space_basis(np.hstack([U, W]), tol)


def basis_of_intersection(U, W, tol=1e-10):
    """Basis of U ∩ W. x in U ∩ W  iff  x = U a = W b  iff  [U | -W] (a; b) = 0.
    For each special solution (a; b) of [U | -W], the vector U a is in both subspaces;
    these vectors span U ∩ W and are independent because U has independent columns."""
    p = U.shape[1]
    N = null_space_basis(np.hstack([U, -W]), tol)        # columns (a; b)
    return U @ N[:p, :]


TOL = 1e-8                                                # explicit rank threshold, see below
for _ in range(100):
    n = 6
    s, u_only, w_only = rng.integers(1, 4, size=3)        # numbers of shared, U-only, W-only columns
    shared = rng.standard_normal((n, s))
    U = np.column_stack([shared, rng.standard_normal((n, u_only))])
    W = np.column_stack([shared, rng.standard_normal((n, w_only))])
    if np.linalg.matrix_rank(np.hstack([U, W]), tol=TOL) < min(n, s + u_only + w_only):
        continue                                          # (measure-zero) degenerate draw
    dU, dW = U.shape[1], W.shape[1]
    dS = basis_of_sum(U, W).shape[1]
    I = basis_of_intersection(U, W)
    dI = np.linalg.matrix_rank(I, tol=TOL) if I.shape[1] else 0
    assert dS == dU + dW - dI, (dS, dU, dW, dI)                       # (3.8), Theorem 3.35
    assert np.linalg.matrix_rank(np.hstack([U, I]), tol=TOL) == dU    # I ⊆ U
    assert np.linalg.matrix_rank(np.hstack([W, I]), tol=TOL) == dW    # I ⊆ W
```

Why the explicit `tol`: `matrix_rank` with no `tol` calls a singular value zero only if it is below about $\sigma_{\max}\cdot\max(m, n)\cdot 2.2\times 10^{-16}$, a few rounding errors. $[W \mid I]$ has exact rank $\dim W$, but its "zero" singular values are computed from rounded data and can land slightly above that line. With `rng = np.random.default_rng(1)`, draw 91 has three shared columns, three U-only and two W-only, so $U = \mathbb{R}^6$, $U \cap W = W$, and the smallest computed singular value of $[W \mid I]$ is $3.7\times 10^{-14}$: the default threshold says rank $6$, not $5$, and the assert fails. A threshold of $10^{-8}$ sits far above rounding noise and far below the singular values of these well-conditioned Gaussian matrices.

Why the intersection method works: $\mathbf{x} \in U \cap W$ means $\mathbf{x} = U\mathbf{a}$ for some $\mathbf{a}$ (in $U$) and $\mathbf{x} = W\mathbf{b}$ for some $\mathbf{b}$ (in $W$); subtracting, $U\mathbf{a} - W\mathbf{b} = \mathbf{0}$, i.e. $[U \mid -W]\binom{\mathbf{a}}{\mathbf{b}} = \mathbf{0}$. Conversely any null vector $(\mathbf{a}; \mathbf{b})$ gives $\mathbf{x} = U\mathbf{a} = W\mathbf{b}$ in both. The map $(\mathbf{a}; \mathbf{b}) \mapsto U\mathbf{a}$ is linear and one-to-one on the null space: if $U\mathbf{a} = \mathbf{0}$ then $\mathbf{a} = \mathbf{0}$ (independent columns), and then $W\mathbf{b} = \mathbf{0}$ forces $\mathbf{b} = \mathbf{0}$. So it carries a basis of the null space to a basis of $U \cap W$, and, with $p = \dim U$ and $q = \dim W$ the column counts as in the exercise, $\dim(U \cap W) = (p + q) - \operatorname{rank}[U \mid -W] = \dim U + \dim W - \dim(U + W)$ (the columns of $[U \mid -W]$ span $U + W$). That is (3.8) again, this time proved by Rank–Nullity. Example: one shared column and two extra columns in each of $U$ and $W$ (so $\dim U = \dim W = 3$) gives $\dim(U + W) = 5$, $\dim(U \cap W) = 1$: $5 = 3 + 3 - 1$. ✓

## C. Challenge

**C1.** (a) $(\Rightarrow)$ If $V = U \oplus W$, Theorem 3.34(b) gives $U + W = V$ and $U \cap W = \{\mathbf{0}\}$, and then (3.8) gives $\dim V = \dim(U + W) = \dim U + \dim W - 0$. $(\Leftarrow)$ If $U \cap W = \{\mathbf{0}\}$ and $\dim U + \dim W = \dim V$, then (3.8) gives $\dim(U + W) = \dim U + \dim W - 0 = \dim V$. $U + W$ is a subspace of $V$ (B1(a)) with the full dimension, so $U + W = V$ (Corollary 3.26(d)). Together with $U \cap W = \{\mathbf{0}\}$, Theorem 3.34(b) gives $V = U \oplus W$. ∎

(b) $U$ has a basis $\mathbf{u}_1, \dots, \mathbf{u}_k$ (Corollary 3.26(d)). Extend it to a basis $\mathbf{u}_1, \dots, \mathbf{u}_k, \mathbf{w}_1, \dots, \mathbf{w}_m$ of $V$ (Theorem 3.25) and let $W = \operatorname{span}\{\mathbf{w}_1, \dots, \mathbf{w}_m\}$, a subspace by Theorem 3.5. *$U + W = V$:* every $\mathbf{v}$ is a combination of the full basis; the part using the $\mathbf{u}$'s is in $U$ and the part using the $\mathbf{w}$'s is in $W$. *$U \cap W = \{\mathbf{0}\}$:* if $\mathbf{x} = \sum_i a_i\mathbf{u}_i = \sum_j b_j\mathbf{w}_j$, then $\sum_i a_i\mathbf{u}_i - \sum_j b_j\mathbf{w}_j = \mathbf{0}$ is a relation among the basis of $V$, so every coefficient is $0$ and $\mathbf{x} = \mathbf{0}$. By Theorem 3.34(b), $V = U \oplus W$. ∎

*Not unique.* In $\mathbb{R}^2$ let $U = \operatorname{span}\{(1, 0)\}$. Both $W_1 = \operatorname{span}\{(0, 1)\}$ and $W_2 = \operatorname{span}\{(1, 1)\}$ are complements: each is a line different from $U$, and two different lines through $\mathbf{0}$ meet only in $\mathbf{0}$ (a nonzero $\mathbf{x}$ in both would make each line equal to $\operatorname{span}\{\mathbf{x}\}$); and $\dim U + \dim W_i = 2 = \dim\mathbb{R}^2$, so (a) applies. LA-5 singles out the orthogonal complement $U^\perp = W_1$: LA-5 Theorem 3.7(b) says $\mathbb{R}^n = S \oplus S^\perp$ for every subspace $S$, and $S^\perp$ is the only complement perpendicular to $S$ (a complement $W$ with $W \perp S$ lies in $S^\perp$ and has the same dimension $n - \dim S$, so equals it by Corollary 3.26(d)).

(c) Take the three lines $U_1 = \operatorname{span}\{(1, 0)\}$, $U_2 = \operatorname{span}\{(0, 1)\}$, $U_3 = \operatorname{span}\{(1, 1)\}$. Any two are different lines, so $U_i \cap U_j = \{\mathbf{0}\}$ for $i \ne j$ (as in (b)), and the triple intersection is $\{\mathbf{0}\}$. The sum contains $U_1 + U_2 = \mathbb{R}^2$, so it is $\mathbb{R}^2$, of dimension $2$. The proposed formula gives $1 + 1 + 1 - 0 + 0 = 3 \ne 2$. The sum is not direct: $(1, 0) + (0, 1) + (-1, -1) = \mathbf{0}$ with $(1, 0) \in U_1$, $(0, 1) \in U_2$, $(-1, -1) \in U_3$, all nonzero, which Theorem 3.34(a) forbids. (Equivalently, a direct sum would have dimension $1 + 1 + 1 = 3$ by Theorem 3.34(c), more than $\dim\mathbb{R}^2$.) The correct count applies (3.8) twice: $\dim(U_1 + U_2 + U_3) = \dim(U_1 + U_2) + \dim U_3 - \dim\big((U_1 + U_2) \cap U_3\big) = 2 + 1 - 1 = 2$. The intersection that matters is $(U_1 + U_2) \cap U_3 = U_3$, which the pairwise intersections do not see.

(d) $U + W$ is a subspace of $\mathbb{R}^7$, so $\dim(U + W) \le 7$ (Corollary 3.26(d)). By (3.8), $\dim(U \cap W) = \dim U + \dim W - \dim(U + W) \ge 4 + 5 - 7 = 2$. ∎ *Equality:* $U = \operatorname{span}\{\mathbf{e}_1, \dots, \mathbf{e}_4\}$, $W = \operatorname{span}\{\mathbf{e}_3, \dots, \mathbf{e}_7\}$. Then $U + W$ contains every $\mathbf{e}_i$, so $U + W = \mathbb{R}^7$. A vector in $U$ has entries $5, 6, 7$ equal to $0$, and a vector in $W$ has entries $1, 2$ equal to $0$, so a vector in both is a combination of $\mathbf{e}_3, \mathbf{e}_4$; both of these lie in $U$ and in $W$. So $U \cap W = \operatorname{span}\{\mathbf{e}_3, \mathbf{e}_4\}$, of dimension $2$. Check: $4 + 5 - 2 = 7$. ✓

**C2.** Induction on $k$. *Base $k = 1$:* $e^{\lambda_1 t}$ is never zero, so $c_1 e^{\lambda_1 t} = 0$ forces $c_1 = 0$ (Theorem 3.13(c)).
*Step.* Assume the statement for any $k - 1$ distinct exponents, and suppose $\sum_{i=1}^k c_i e^{\lambda_i t} = 0$ for all $t$. Multiply by $e^{-\lambda_k t}$ (nonzero, so the equation is equivalent): $\sum_{i=1}^{k-1} c_i e^{(\lambda_i - \lambda_k)t} + c_k = 0$ for all $t$. Differentiate: $\sum_{i=1}^{k-1} c_i(\lambda_i - \lambda_k)e^{(\lambda_i - \lambda_k)t} = 0$. The $k - 1$ exponents $\mu_i = \lambda_i - \lambda_k$ are distinct (the $\lambda_i$ are), so by the induction hypothesis every coefficient $c_i(\lambda_i - \lambda_k)$ is $0$; and $\lambda_i - \lambda_k \ne 0$ for $i < k$, so $c_i = 0$ for $i = 1, \dots, k - 1$. Substituting back, $c_k = 0$. ∎

*Complex case.* Nothing in the proof used that $\lambda_i$ or $c_i$ were real: $e^{\mu t} \ne 0$ for complex $\mu$, differentiation of $e^{\mu t}$ gives $\mu e^{\mu t}$, and "distinct" is all that was needed. So over the complex scalars, $e^{\lambda_1 t}, \dots, e^{\lambda_k t}$ with distinct complex $\lambda_i$ are independent in the complex vector space of smooth complex-valued functions.

*Trigonometric functions.* Suppose $a_0 + a_1\cos t + b_1\sin t + a_2\cos 2t + b_2\sin 2t = 0$ for all $t$, with real coefficients. Write $\cos kt = \tfrac12(e^{ikt} + e^{-ikt})$ and $\sin kt = \tfrac{1}{2i}(e^{ikt} - e^{-ikt})$. The relation becomes
$$a_0 e^{0t} + \tfrac{a_1 - ib_1}{2}e^{it} + \tfrac{a_1 + ib_1}{2}e^{-it} + \tfrac{a_2 - ib_2}{2}e^{2it} + \tfrac{a_2 + ib_2}{2}e^{-2it} = 0,$$
a relation among $e^{\lambda t}$ for the five distinct complex exponents $0, \pm i, \pm 2i$. By the complex case all five coefficients vanish: $a_0 = 0$, $a_1 \mp ib_1 = 0$ (so $a_1 = b_1 = 0$), $a_2 \mp ib_2 = 0$ (so $a_2 = b_2 = 0$). ∎ The same argument for all $k \le K$ shows the first $2K + 1$ Fourier basis functions are independent, which is why Fourier coefficients (F2) are well-defined coordinates. In LA-7, exponentials with distinct eigenvalues as exponents are the basis of solutions of $\dot{\mathbf{x}} = A\mathbf{x}$, and this exercise is why those solutions are independent.

**C3.** (a) Column $j$ of $A + B$ is $\mathbf{a}_j + \mathbf{b}_j \in C(A) + C(B)$, so $C(A + B) \subseteq C(A) + C(B)$ (Theorem 3.5 and B1(a)). A basis of $C(A)$ followed by a basis of $C(B)$ spans $C(A) + C(B)$ (every $\mathbf{u} + \mathbf{w}$ is a combination of the two lists), so $\dim(C(A) + C(B)) \le \operatorname{rank} A + \operatorname{rank} B$ by Theorem 3.22 (a spanning list contains a basis, which is no longer than the list). Then Corollary 3.26(d): $\operatorname{rank}(A + B) = \dim C(A + B) \le \dim(C(A) + C(B)) \le \operatorname{rank} A + \operatorname{rank} B$. ∎
(b) By (3.7) twice, $\operatorname{rank}(PAQ) = \operatorname{rank}\big((PA)Q\big) \le \operatorname{rank}(PA) \le \operatorname{rank} A$. Also $A = P^{-1}(PAQ)Q^{-1}$, so by the same two steps $\operatorname{rank} A \le \operatorname{rank}(PAQ)$. Equality. ∎ (Special case: row operations, $P = E$, $Q = I$, preserve rank, which is why "rank $=$ number of pivots" was legitimate in LA-2.)
(c) Code and the numbers it produces (run after the main script, or paste into a cell):

```python
for eps in (1e-6, 1e-12, 1e-16):
    A = np.array([[1, 1], [1, 1 + eps]])
    print(eps, len(rref(A, tol=1e-10)[1]), np.linalg.matrix_rank(A), f"{np.linalg.det(A):.6e}", 1 + eps == 1)
# Output (columns: eps, rref rank, matrix_rank, det, 1 + eps == 1):
# 1e-06 2 2 1.000000e-06 False
# 1e-12 1 2 1.000089e-12 False
# 1e-16 1 1 0.000000e+00 True
```

*$\varepsilon = 10^{-6}$.* Elimination leaves $\varepsilon$ in the $(2, 2)$ slot; $10^{-6} > 10^{-10}$, so `rref` counts two pivots. The singular values are about $2$ and $\varepsilon/2 = 5 \times 10^{-7}$, far above `matrix_rank`'s default threshold (about $2 \cdot 2.2 \times 10^{-16} \cdot 2 \approx 10^{-15}$), so it also says $2$. The determinant is $\varepsilon$. All agree: rank $2$, but the matrix is close (within $10^{-6}$) to a rank-$1$ matrix.
*$\varepsilon = 10^{-12}$.* Elimination leaves $10^{-12}$ in the $(2, 2)$ slot, below the tolerance $10^{-10}$, so `rref` declares it zero: rank $1$. The SVD sees $\sigma_2 \approx 5 \times 10^{-13}$, still above its threshold: rank $2$. The determinant is $1.000089 \times 10^{-12} \ne 0$. Notice the fifth digit. The stored entry is not $1 + 10^{-12}$ but the nearest double, which is off by about $8.9 \times 10^{-17}$ (doubles near $1$ are $2.2 \times 10^{-16}$ apart). That error is tiny next to $1$ and is $9 \times 10^{-5}$ of $\varepsilon$; the subtraction $(1 + \varepsilon) - 1$, which both elimination and the determinant perform, cancels the $1$'s and leaves the error exposed (catastrophic cancellation, LA-13 §3.2). So even the computed "$\varepsilon$" is only good to four digits. The tools disagree because "rank" of a floating-point matrix is a decision about which small numbers to call zero, and the two tools draw the line in different places ($10^{-10}$ versus about $10^{-15}$).
*$\varepsilon = 10^{-16}$.* In double precision $1 + 10^{-16}$ rounds to exactly $1$ (the spacing of doubles near $1$ is $2.2 \times 10^{-16}$), so the stored matrix is exactly $\begin{pmatrix} 1 & 1 \\ 1 & 1\end{pmatrix}$ and every method says rank $1$; the "true" matrix you typed never existed in the machine.

Which is right? If the entries are exact, the mathematical rank is $2$ for every $\varepsilon \ne 0$ and only the $10^{-16}$ case is wrong (and it is the *input*, not the algorithm, that was wrong). If the entries are measurements with noise around $10^{-9}$, the honest answer is "rank $1$ to within the noise," and `rref` with a tolerance matching the noise is the right tool, while `matrix_rank` with its default threshold is overconfident. The decision needs a number (the noise level) that the matrix does not carry. Underneath is the fact that rank is not a continuous function of the entries: $A_\varepsilon$ has rank $2$ for every $\varepsilon \ne 0$ and rank $1$ at $\varepsilon = 0$, so any floating-point rank needs a threshold. The SVD's singular values are the right way to see *how far* a matrix is from each lower rank: LA-9 Theorem 3.8 (Eckart–Young, spectral norm) proves that the distance from $A$ to the nearest matrix of rank $\le k$ is exactly $\sigma_{k+1}$. Here $\sigma_2 \approx \varepsilon/2$, so $A_\varepsilon$ is within $\varepsilon/2$ of a rank-$1$ matrix, and the honest report is "rank $2$, but only $\varepsilon/2$ away from rank $1$."
