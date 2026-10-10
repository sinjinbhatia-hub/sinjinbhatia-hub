# LA-6: Determinants: Signed Volume, Properties, Cofactors, Cramer; the Alternating-Multilinear View

**Time:** ~8 hours · **Prerequisites:** LA-4 (and LA-2 for elimination, LA-3 for independence) · **Unlocks:** LA-7

> **Conventions in this module.** $A$ is a real $n \times n$ matrix with entries $a_{ij}$ (row $i$, column $j$), rows $\mathbf{r}_1^\top, \dots, \mathbf{r}_n^\top$ and columns $\mathbf{a}_1, \dots, \mathbf{a}_n$. We write $\det A$ or $|A|$ for the determinant and define it by its action on **rows**; every row statement transfers to columns once we prove $\det A^\top = \det A$ (Theorem 3.11). A permutation $\sigma$ of $\{1, \dots, n\}$ is written in one-line notation $(\sigma(1), \dots, \sigma(n))$; its sign is $\operatorname{sgn}(\sigma) = (-1)^{N(\sigma)}$ where $N(\sigma)$ counts inversions (pairs $i < j$ with $\sigma(i) > \sigma(j)$). The minor $M_{ij}$ is the determinant of the $(n-1)\times(n-1)$ matrix obtained by deleting row $i$ and column $j$; the cofactor is $C_{ij} = (-1)^{i+j} M_{ij}$; the adjugate is $\operatorname{adj} A = C^\top$ (transpose of the cofactor matrix). Orientation: an ordered basis of $\mathbb{R}^n$ is **positively oriented** if the matrix with those vectors as columns has positive determinant; in $\mathbb{R}^2$ that means counterclockwise from the first to the second, in $\mathbb{R}^3$ right-handed. Volume means $n$-dimensional volume (area when $n = 2$).

## 0. Where this sits

You already know, from LA-2, that a square matrix is invertible exactly when elimination produces $n$ nonzero pivots, and from LA-4 that a matrix is the coordinate description of a linear map. This module attaches a single number to a square matrix that answers the invertibility question ("is it zero?"), measures how the map scales volume ("what is its absolute value?"), and records whether the map flips orientation ("what is its sign?").

You need: elimination with row swaps (LA-2), linear independence and dimension (LA-3), the idea that a linear map is determined by what it does to a basis (LA-4, Theorem 3.5 of LA-1). You will use: LA-7 defines eigenvalues as the roots of $\det(A - \lambda I)$ and cannot start without this module. LA-8 (positive definiteness via minors), LA-12 (alternating forms and tensors), F5 (the $\det\Sigma$ in the Gaussian density and the Jacobian in change of variables), and CM-8 (Liouville's theorem) all reach back here.

Honest sizing: the computations are easy. The two ideas that take real effort are the proof that the three defining properties pin down a unique function, and the alternating-multilinear characterization in §3.9. Budget your eight hours accordingly: two on §2–3.4, two on §3.5–3.8, two on §3.9–3.13, two on exercises and the project.

## 1. The problem that forced this idea

The determinant was discovered twice, independently, as a by-product of eliminating unknowns from linear systems. In Japan, Seki Takakazu's *Kai Fukudai no Hō* (1683) gave rules for eliminating variables from systems of equations, and those rules contain, for small systems, exactly the expressions we now call determinants. In Europe, Leibniz wrote to l'Hôpital in 1693 about three equations in two unknowns. Using a two-digit index notation ($10, 11, 12, 20, 21, 22, 30, 31, 32$ for the coefficients) he wrote down the condition for the three equations to be simultaneously solvable, and the condition is a $3 \times 3$ determinant set equal to zero. Neither Seki nor Leibniz had a name for the object. They had a number that comes out of elimination and tells you whether elimination will succeed.

The next step came from geometry. Gabriel Cramer, in an appendix to his 1750 *Introduction à l'analyse des lignes courbes algébriques*, asked how to find the algebraic curve of a given degree passing through given points. A conic $ax^2 + bxy + cy^2 + dx + ey + f = 0$ through five points gives five linear equations in the six unknown coefficients (one is a normalization). Cramer wrote the general solution of $n$ equations in $n$ unknowns as a ratio of two expressions, each a signed sum of products of coefficients: the rule that now carries his name. (Maclaurin had the $2\times 2$ and $3\times 3$ cases in his *Treatise of Algebra*, published posthumously in 1748.) Cramer's denominator is the determinant of the coefficient matrix, and the rule fails exactly when it vanishes. So the determinant entered mathematics as the **obstacle to solvability**: the one number whose vanishing means the system is degenerate.

Through the late eighteenth century the object acquired its own theory. Vandermonde (1771) treated these expressions as things to study in their own right, Laplace (1772) gave the expansion along a row in terms of smaller determinants, and Gauss used the word *determinant* in 1801 for the discriminant of a quadratic form. The decisive paper is Cauchy's 1812 memoir, read to the Institut de France, which uses "determinant" in the modern sense, arranges the coefficients in a square array, introduces the row-and-column notation, and proves the multiplication theorem $\det(AB) = \det A \cdot \det B$ (Binet announced a related theorem the same day). Jacobi's three 1841 papers in Crelle's journal systematized the theory and introduced the *functional determinant*, the object you call the Jacobian.

Cauchy's definition is the sum over permutations (our Theorem 3.9). It is correct, it is complete, and it is useless for understanding. The question "what *is* a determinant?" was answered differently, late in the nineteenth century, in the Berlin lectures of Weierstrass and Kronecker: the determinant is the **unique** function of $n$ vectors in $\mathbb{R}^n$ that is linear in each vector, changes sign when two vectors are swapped, and equals $1$ on the standard basis. Every property of determinants follows from those three sentences, and so does the formula. That is the view this module takes. Strang's chapter 5 and Axler's chapter 9 are the two modern books that take it seriously, from opposite ends: Strang starts from the three properties and computes, Axler starts from alternating multilinear maps and proves. We will do both.

The payoff for you, concretely: when LA-7 defines eigenvalues as the roots of $\det(A - \lambda I)$, you will know why that polynomial does not depend on the basis you wrote $A$ in. When F5 writes the Gaussian density with a $(\det\Sigma)^{-1/2}$ in front, you will know that it is a volume. When you change variables in a triple integral and write $r^2\sin\phi\,dr\,d\phi\,d\theta$, you will know you are computing a determinant.

## 2. Intuition

### Area as a function of two columns

Take two vectors $\mathbf{a} = (a_1, a_2)^\top$ and $\mathbf{b} = (b_1, b_2)^\top$ in $\mathbb{R}^2$. They span a parallelogram with vertices $\mathbf{0}, \mathbf{a}, \mathbf{a} + \mathbf{b}, \mathbf{b}$.

```
                 a+b
          b  ___________
            /          /
           /          /
          /          /
         /__________/
        0            a
```

Call its area $\operatorname{Area}(\mathbf{a}, \mathbf{b})$, a function of two vector arguments. We want to understand how it responds to the two vector operations from LA-1, scaling and adding, applied to one argument at a time.

**Scaling one side.** Replace $\mathbf{a}$ by $c\mathbf{a}$ with $c > 0$. The base stretches by $c$, the height (the perpendicular distance from the line through $\mathbf{a}$ to the tip of $\mathbf{b}$) is unchanged. So $\operatorname{Area}(c\mathbf{a}, \mathbf{b}) = c\operatorname{Area}(\mathbf{a}, \mathbf{b})$.

**Adding to one side.** Write the area as base times height, where the height is measured perpendicular to $\mathbf{b}$:
$$\operatorname{Area}(\mathbf{a}, \mathbf{b}) = \|\mathbf{b}\| \cdot h(\mathbf{a}), \qquad h(\mathbf{a}) = \text{component of } \mathbf{a} \text{ perpendicular to } \mathbf{b}.$$
Let $\mathbf{n}$ be the unit vector obtained by rotating $\mathbf{b}/\|\mathbf{b}\|$ by $90°$ counterclockwise. Then $h(\mathbf{a}) = \mathbf{a} \cdot \mathbf{n}$, and the dot product is linear in $\mathbf{a}$ (LA-1, §3.5). So
$$\operatorname{Area}(\mathbf{a} + \mathbf{a}', \mathbf{b}) = \|\mathbf{b}\|\,(\mathbf{a} + \mathbf{a}')\cdot\mathbf{n} = \operatorname{Area}(\mathbf{a}, \mathbf{b}) + \operatorname{Area}(\mathbf{a}', \mathbf{b}).$$
Area is **additive** in the first argument. You can see the same fact without the dot product: shear the parallelogram on $\mathbf{a} + \mathbf{a}'$ along the direction of $\mathbf{b}$ and it decomposes into the parallelogram on $\mathbf{a}$ and the one on $\mathbf{a}'$, each with the same base $\mathbf{b}$.

But there is a catch in the formula $h(\mathbf{a}) = \mathbf{a} \cdot \mathbf{n}$: it can be **negative**. If $\mathbf{a}$ is on the clockwise side of $\mathbf{b}$, the "height" comes out negative, and $\|\mathbf{b}\|\,h(\mathbf{a})$ is minus the geometric area. The formula is not computing area. It is computing **signed area**: positive when going from $\mathbf{a}$ to $\mathbf{b}$ is a counterclockwise turn, negative when it is clockwise.

```
       b                            a
       ^                            ^
       |   positive:                |   negative:
       |   a -> b is CCW            |   a -> b is CW
       +--------> a                 +--------> b
```

This is forced on us. If we want area to be linear in each argument, including the step $\operatorname{Area}(-\mathbf{a}, \mathbf{b}) = -\operatorname{Area}(\mathbf{a}, \mathbf{b})$ that linearity demands, area must be allowed to be negative. Signed area is the price of linearity, and it turns out to be a feature: the sign records **orientation**.

**Swapping the two sides.** The parallelogram on $(\mathbf{b}, \mathbf{a})$ is the same set of points as the one on $(\mathbf{a}, \mathbf{b})$, so the unsigned area is the same. But going from $\mathbf{b}$ to $\mathbf{a}$ is the opposite turn. So $\operatorname{Area}(\mathbf{b}, \mathbf{a}) = -\operatorname{Area}(\mathbf{a}, \mathbf{b})$.

**Normalization.** The unit square has area $1$ and $(\mathbf{e}_1, \mathbf{e}_2)$ is counterclockwise: $\operatorname{Area}(\mathbf{e}_1, \mathbf{e}_2) = 1$.

### The three properties already determine the formula

Now watch what these three facts (linear in each argument, sign flips under swap, $1$ on the standard basis) force. Expand each argument in the standard basis and use linearity twice:
$$\operatorname{Area}(\mathbf{a}, \mathbf{b}) = \operatorname{Area}(a_1\mathbf{e}_1 + a_2\mathbf{e}_2,\; b_1\mathbf{e}_1 + b_2\mathbf{e}_2) = \sum_{i=1}^2\sum_{j=1}^2 a_i b_j \operatorname{Area}(\mathbf{e}_i, \mathbf{e}_j).$$
Four terms. Two of them have $\mathbf{e}_i = \mathbf{e}_j$: a parallelogram on two copies of the same vector is flat, and the swap rule confirms it ($\operatorname{Area}(\mathbf{e}_1, \mathbf{e}_1) = -\operatorname{Area}(\mathbf{e}_1, \mathbf{e}_1)$, so it is $0$). The other two are $\operatorname{Area}(\mathbf{e}_1, \mathbf{e}_2) = 1$ and $\operatorname{Area}(\mathbf{e}_2, \mathbf{e}_1) = -1$. So
$$\operatorname{Area}(\mathbf{a}, \mathbf{b}) = a_1 b_2 - a_2 b_1.$$
That is the $2\times 2$ determinant, and we did not define it; we derived it from three properties any sensible signed area must have. Nothing in the argument used the number $2$ in an essential way. In $\mathbb{R}^n$, $n$ vectors span a parallelepiped, its signed volume should be linear in each edge, should flip sign when two edges are swapped, and should be $1$ on the unit cube, and the same expansion will produce a sum of $n!$ signed products. That is the main theorem of §3.

### A linear map scales every area by the same factor

One more thing you should be able to guess. Let $A$ be a $2\times 2$ matrix and think of it as the map $\mathbf{x} \mapsto A\mathbf{x}$. It sends the unit square (edges $\mathbf{e}_1, \mathbf{e}_2$) to the parallelogram with edges $A\mathbf{e}_1 = \mathbf{a}_1$ and $A\mathbf{e}_2 = \mathbf{a}_2$, the columns of $A$. So the image of the unit square has signed area $\det A = a_{11}a_{22} - a_{12}a_{21}$.

What about some other region? Chop it into tiny squares of side $h$. The map sends each tiny square to a tiny parallelogram, a translate of $h$ times the image of the unit square, with area $h^2|\det A|$. Every piece is scaled by the same factor. So the whole region is scaled by $|\det A|$, and if $\det A < 0$ the map has turned every counterclockwise loop into a clockwise one. The determinant is the **volume scaling factor** of the linear map, with a sign for orientation. §3.10 proves this.

If $\det A = 0$, the two columns are parallel (or one is zero), the parallelogram is flat, and the map squashes the plane onto a line. A squashed plane cannot be unsquashed by any linear map: $A$ is not invertible. So you should also expect $\det A = 0 \iff A$ singular. That is the obstacle Leibniz and Cramer found.

## 3. Formal development

### 3.1 The definition: three properties

**Definition 3.1 (Determinant).** The determinant is a function $\det : \mathbb{R}^{n\times n} \to \mathbb{R}$ satisfying three properties:

- **(D1) Normalization.** $\det I = 1$.
- **(D2) Alternating.** Swapping two rows of $A$ multiplies $\det A$ by $-1$.
- **(D3) Multilinear.** $\det$ is linear in each row separately, with the other rows held fixed:
  $$\det\begin{pmatrix} \vdots \\ c\,\mathbf{r}^\top + \mathbf{r}'^\top \\ \vdots \end{pmatrix} = c\det\begin{pmatrix} \vdots \\ \mathbf{r}^\top \\ \vdots \end{pmatrix} + \det\begin{pmatrix} \vdots \\ \mathbf{r}'^\top \\ \vdots \end{pmatrix}. \tag{3.1}$$

Two warnings about (D3). It is linearity in **one row at a time**: $\det(A + B) \neq \det A + \det B$ in general, because adding $B$ changes every row at once. And scaling the whole matrix scales every row, so $\det(cA) = c^n \det A$, not $c\det A$. Both mistakes are common; §8 returns to them.

The definition does not say that such a function exists or that there is only one. We prove both. The strategy: §3.2 derives consequences of (D1)–(D3), assuming a determinant exists. §3.3 shows that these consequences force a specific value on every matrix (so there is at most one determinant) and then exhibits a formula that satisfies (D1)–(D3) (so there is at least one). After that, every result about the determinant can be proved either from the three properties or from the formula, whichever is shorter.

### 3.2 Consequences of the three properties

Throughout this subsection, $\det$ is any function satisfying (D1)–(D3).

**Theorem 3.2 (Two equal rows).** If two rows of $A$ are equal, $\det A = 0$.

*Proof.* Swap the two equal rows. The matrix is unchanged, so its determinant is unchanged. But (D2) says the determinant was multiplied by $-1$. So $\det A = -\det A$, hence $\det A = 0$. ∎

(This one-line argument uses that $2 \neq 0$ in $\mathbb{R}$. Over fields where $1 + 1 = 0$ it fails, which is why Axler takes "two equal arguments gives zero" as the axiom and derives the sign flip from it; see §3.9 step 1. Over $\mathbb{R}$ the two formulations are equivalent.)

**Theorem 3.3 (Row operations).** Subtracting a multiple of one row from another row does not change the determinant.

*Proof.* Replace row $i$ by $\mathbf{r}_i - c\,\mathbf{r}_j$ with $j \neq i$. By (D3) in row $i$,
$$\det\begin{pmatrix}\vdots \\ \mathbf{r}_i^\top - c\,\mathbf{r}_j^\top \\ \vdots\end{pmatrix} = \det\begin{pmatrix}\vdots \\ \mathbf{r}_i^\top \\ \vdots\end{pmatrix} - c\det\begin{pmatrix}\vdots \\ \mathbf{r}_j^\top \\ \vdots\end{pmatrix},$$
where the last matrix has $\mathbf{r}_j^\top$ in row $i$ **and** in row $j$. By Theorem 3.2 its determinant is $0$. ∎

**Theorem 3.4 (Zero row).** If some row of $A$ is zero, $\det A = 0$.

*Proof.* The zero row equals $0 \cdot \mathbf{r}^\top$ for any $\mathbf{r}$. By (D3) with $c = 0$ and no $\mathbf{r}'$ term, $\det A = 0 \cdot \det(\cdots) = 0$. ∎

**Theorem 3.5 (Dependent rows).** If the rows of $A$ are linearly dependent, $\det A = 0$.

*Proof.* Some row is a combination of the others, say $\mathbf{r}_i = \sum_{j \neq i} c_j \mathbf{r}_j$ (LA-3). By (D3) in row $i$, $\det A = \sum_{j\neq i} c_j \det A^{(j)}$, where $A^{(j)}$ is $A$ with row $i$ replaced by $\mathbf{r}_j$. Each $A^{(j)}$ has two equal rows (rows $i$ and $j$), so each term is $0$ by Theorem 3.2. ∎

**Theorem 3.6 (Triangular matrices).** If $U$ is upper or lower triangular, $\det U = u_{11} u_{22} \cdots u_{nn}$, the product of the diagonal entries.

▶ **Try it first.** If every diagonal entry is nonzero, use Theorem 3.3 to clear everything off the diagonal, then (D3) to pull out each diagonal entry. What goes wrong if some $u_{kk} = 0$?

*Proof.* Take $U$ upper triangular (the lower case is identical with the roles of "above" and "below" exchanged).

1. *Case: all $u_{kk} \neq 0$.* Work from the bottom row up. Row $n$ is $(0, \dots, 0, u_{nn})$. Subtract $u_{k n}/u_{nn}$ times row $n$ from each row $k < n$: this clears column $n$ above the diagonal and, by Theorem 3.3, does not change the determinant. Now row $n-1$ is $(0, \dots, u_{n-1,n-1}, 0)$; use it to clear column $n-1$. Continue. After $n - 1$ rounds the matrix is diagonal, $D = \operatorname{diag}(u_{11}, \dots, u_{nn})$, with $\det D = \det U$.
2. Row $k$ of $D$ is $u_{kk}\mathbf{e}_k^\top$. Apply (D3) to each row in turn to pull out the scalar: $\det D = u_{11}u_{22}\cdots u_{nn}\det I = u_{11}\cdots u_{nn}$ by (D1).
3. *Case: some $u_{kk} = 0$.* Then rows $k, k+1, \dots, n$ (that is $n - k + 1$ rows) have zeros in columns $1, \dots, k$: rows below $k$ because $U$ is triangular, row $k$ because its first nonzero entry would have been $u_{kk}$. So these $n - k + 1$ rows all lie in $\operatorname{span}\{\mathbf{e}_{k+1}, \dots, \mathbf{e}_n\}$, a space of dimension $n - k$. More vectors than the dimension means dependent (LA-3), and Theorem 3.5 gives $\det U = 0$, which equals the product of the diagonal since one factor is $0$. ∎

■ **What just happened.** Theorems 3.3 and 3.6 together say: the determinant survives elimination untouched and is then read off the diagonal. That is the algorithm.

**Theorem 3.7 (Determinant by elimination; singular iff zero).** Run Gaussian elimination with row swaps on $A$ (LA-2), producing an upper triangular $U$ with $s$ row swaps along the way. Then
$$\det A = (-1)^s\, u_{11} u_{22} \cdots u_{nn} = (-1)^s \times (\text{product of the pivots}). \tag{3.2}$$
In particular, $\det A \neq 0$ if and only if $A$ is invertible.

*Proof.*

1. Elimination uses two kinds of operations: swapping two rows (each one multiplies the determinant by $-1$, by (D2)) and subtracting a multiple of a pivot row from a lower row (no change, Theorem 3.3). After $s$ swaps and any number of subtractions, $\det U = (-1)^s \det A$.
2. $\det U$ is the product of its diagonal (Theorem 3.6). The diagonal of $U$ holds the pivots, with $0$ in any position where elimination found no pivot. Multiply through by $(-1)^s$ to get (3.2).
3. From LA-2: $A$ is invertible iff elimination produces $n$ nonzero pivots iff $u_{11}\cdots u_{nn} \neq 0$ iff $\det A \neq 0$. ∎

■ **What just happened.** Two things. First, the practical algorithm: eliminate, count swaps, multiply the diagonal. For a $1000 \times 1000$ matrix this costs about $\tfrac23 n^3 \approx 7 \times 10^8$ flops, well under a second. Second, a uniqueness proof in disguise: step 1 used only (D2) and (D3), and step 2 used only (D1)–(D3). So **any** function satisfying the three properties must return $(-1)^s u_{11}\cdots u_{nn}$ on $A$. If two functions both satisfy (D1)–(D3), they agree on every matrix. There is at most one determinant.

What we have not yet shown is that there is at least one. Different elimination orders (different swap choices) could in principle give different values of $(-1)^s\prod u_{kk}$, and if they did, no function could satisfy (D1)–(D3) at all. The next subsection removes that worry by writing down a formula and checking the three properties directly.

### 3.3 The big formula: existence and uniqueness

**Permutations and their signs.** A permutation $\sigma$ of $\{1, \dots, n\}$ is a bijection from that set to itself; we write it as the list $(\sigma(1), \dots, \sigma(n))$. There are $n!$ of them. An **inversion** of $\sigma$ is a pair of positions $i < j$ with $\sigma(i) > \sigma(j)$; $N(\sigma)$ is the number of inversions, and $\operatorname{sgn}(\sigma) = (-1)^{N(\sigma)}$. For example $(2, 3, 1)$ has inversions $(2,1)$ and $(3,1)$, so $N = 2$ and $\operatorname{sgn} = +1$; $(3, 2, 1)$ has $N = 3$ and $\operatorname{sgn} = -1$.

**Lemma 3.8 (A swap flips the sign).** Swapping two entries of the list $(\sigma(1), \dots, \sigma(n))$ changes $N(\sigma)$ by an odd number, hence flips $\operatorname{sgn}(\sigma)$.

*Proof.*

1. *Adjacent swap.* Swap the entries in positions $k$ and $k+1$. For any pair of positions not equal to $\{k, k+1\}$, the two entries keep their relative order (a pair involving position $k$ and some $m \notin \{k, k+1\}$ compares the same two values, in the same left-to-right order, as before). Only the pair $(k, k+1)$ changes status: it becomes an inversion if it was not one, and stops being one if it was. So $N$ changes by exactly $\pm 1$.
2. *General swap.* Swap positions $i < j$, with $d = j - i$. Do it with adjacent swaps: move the entry at position $j$ leftward to position $i+1$ ($d - 1$ adjacent swaps), swap positions $i$ and $i + 1$ (one swap), then move the entry now at position $i+1$ (the one originally at $i$) rightward to position $j$ ($d - 1$ swaps). Total $2d - 1$ adjacent swaps, an odd number, each changing $N$ by $\pm 1$. The net change is odd. ∎

Let $P_\sigma$ be the **permutation matrix** with rows $\mathbf{e}_{\sigma(1)}^\top, \dots, \mathbf{e}_{\sigma(n)}^\top$. Sorting its rows into the order $\mathbf{e}_1^\top, \dots, \mathbf{e}_n^\top$ by $s$ swaps turns it into $I$ and turns the list $\sigma$ into $(1, \dots, n)$, which has $N = 0$. Each swap flips the parity of $N$ (Lemma 3.8), so $(-1)^s = (-1)^{N(\sigma)}$. By (D2) and (D1), if a determinant exists then
$$\det P_\sigma = (-1)^s = \operatorname{sgn}(\sigma).$$

**Theorem 3.9 (Leibniz formula; existence and uniqueness).** There is exactly one function satisfying (D1)–(D3), and it is
$$\det A = \sum_{\sigma} \operatorname{sgn}(\sigma)\, a_{1\sigma(1)}\, a_{2\sigma(2)} \cdots a_{n\sigma(n)}, \tag{3.3}$$
the sum running over all $n!$ permutations of $\{1, \dots, n\}$.

▶ **Try it first.** Write row $i$ as $\mathbf{r}_i^\top = \sum_j a_{ij}\mathbf{e}_j^\top$ and expand $\det A$ using (D3) in every row. You get $n^n$ terms. Which ones survive?

*Proof of uniqueness (the formula is forced).* Suppose $\det$ satisfies (D1)–(D3).

1. Expand row 1: $\mathbf{r}_1^\top = \sum_{j_1} a_{1j_1}\mathbf{e}_{j_1}^\top$. By (D3), $\det A = \sum_{j_1} a_{1j_1}\det(\mathbf{e}_{j_1}^\top; \mathbf{r}_2^\top; \dots; \mathbf{r}_n^\top)$, where the semicolons list the rows.
2. Expand row 2 inside each term, then row 3, and so on. After all $n$ rows:
   $$\det A = \sum_{j_1=1}^n\cdots\sum_{j_n=1}^n a_{1j_1}a_{2j_2}\cdots a_{nj_n}\,\det(\mathbf{e}_{j_1}^\top; \dots; \mathbf{e}_{j_n}^\top).$$
3. If two of the indices $j_1, \dots, j_n$ are equal, the matrix $(\mathbf{e}_{j_1}^\top; \dots; \mathbf{e}_{j_n}^\top)$ has two equal rows and its determinant is $0$ (Theorem 3.2). The surviving terms have $j_1, \dots, j_n$ all distinct, so $(j_1, \dots, j_n)$ is a permutation $\sigma$ and the matrix is $P_\sigma$.
4. $\det P_\sigma = \operatorname{sgn}(\sigma)$ from the computation above. Substituting gives (3.3). So any function with properties (D1)–(D3) equals the right side of (3.3). ∎

*Proof of existence (the formula has the three properties).* Define $f(A)$ to be the right side of (3.3).

5. *(D1).* For $A = I$, the product $\prod_i a_{i\sigma(i)}$ is $0$ unless $\sigma(i) = i$ for every $i$. Only the identity permutation survives, with sign $+1$ and product $1$. So $f(I) = 1$.
6. *(D3).* Each product $a_{1\sigma(1)}\cdots a_{n\sigma(n)}$ contains exactly one factor from row $i$, namely $a_{i\sigma(i)}$. So each term is a linear function of row $i$ (the other rows fixed), and a sum of linear functions is linear.
7. *(D2).* Let $A'$ be $A$ with rows $i$ and $j$ swapped, and let $\tau$ be the transposition that exchanges $i$ and $j$, so $a'_{k m} = a_{\tau(k) m}$. Then
   $$\prod_{k} a'_{k\sigma(k)} = \prod_{k} a_{\tau(k)\sigma(k)} = \prod_{m} a_{m\,\sigma(\tau^{-1}(m))} = \prod_m a_{m\,(\sigma\circ\tau)(m)},$$
   substituting $m = \tau(k)$ and using $\tau^{-1} = \tau$. The list for $\sigma\circ\tau$ is the list for $\sigma$ with the entries in positions $i$ and $j$ swapped, so $\operatorname{sgn}(\sigma\circ\tau) = -\operatorname{sgn}(\sigma)$ by Lemma 3.8. As $\sigma$ runs over all permutations so does $\sigma\circ\tau$. Hence
   $$f(A') = \sum_\sigma \operatorname{sgn}(\sigma)\prod_m a_{m(\sigma\tau)(m)} = \sum_{\rho}\operatorname{sgn}(\rho\tau)\prod_m a_{m\rho(m)} = -\sum_\rho\operatorname{sgn}(\rho)\prod_m a_{m\rho(m)} = -f(A),$$
   with $\rho = \sigma\tau$. ∎

■ **What just happened.** The three properties are consistent (there is a function with all three) and complete (there is only one). From now on "the determinant" is a well-defined object and we may prove things about it using whichever description is convenient: the three properties, the elimination formula (3.2), or the Leibniz formula (3.3). The formula also tells you what kind of object $\det A$ is: a polynomial in the $n^2$ entries, homogeneous of degree $n$, with each term using exactly one entry from every row and every column. For $n = 2$ it has $2$ terms, for $n = 3$ it has $6$, for $n = 10$ it has $3{,}628{,}800$. Nobody computes with (3.3) past $n = 4$; §3.13 says why, and §3.9 says what it is good for instead.

### 3.4 Products, inverses, transposes

**Theorem 3.10 (Product rule).** For $n \times n$ matrices $A$ and $B$, $\det(AB) = \det A \cdot \det B$.

▶ **Try it first.** Fix $B$ with $\det B \neq 0$ and consider $f(A) = \det(AB)/\det B$ as a function of $A$. Check (D1)–(D3) for $f$. Then what does uniqueness say?

*Proof.*

1. *Case $\det B \neq 0$.* Define $f(A) = \det(AB)/\det B$. Row $i$ of $AB$ is $\mathbf{r}_i^\top B$ (LA-1, row form of a product), so each row of $AB$ depends only on the corresponding row of $A$, and linearly.
2. *(D1):* $f(I) = \det(IB)/\det B = 1$.
3. *(D2):* swapping rows $i, j$ of $A$ swaps rows $i, j$ of $AB$, which flips $\det(AB)$ and hence $f$.
4. *(D3):* if row $i$ of $A$ is $c\,\mathbf{r}^\top + \mathbf{r}'^\top$, row $i$ of $AB$ is $c\,\mathbf{r}^\top B + \mathbf{r}'^\top B$, and (D3) for $\det$ applied to $AB$ gives (D3) for $f$.
5. By uniqueness (Theorem 3.9), $f = \det$. So $\det(AB)/\det B = \det A$.
6. *Case $\det B = 0$.* Then $B$ is singular (Theorem 3.7), so there is $\mathbf{x} \neq \mathbf{0}$ with $B\mathbf{x} = \mathbf{0}$. Then $(AB)\mathbf{x} = A\mathbf{0} = \mathbf{0}$, so $AB$ is singular too and $\det(AB) = 0 = \det A \cdot 0$. ∎

■ **What just happened.** The proof never touched a formula. It said: "$A \mapsto \det(AB)/\det B$ behaves like a determinant, so it is one." This is the characteristic move of the whole module, and it is why uniqueness matters more than the formula. (An alternative proof writes $A$ as a product of elementary matrices and checks the rule for each; it is longer and teaches less.)

**Corollary 3.10′.** (a) If $A$ is invertible, $\det(A^{-1}) = 1/\det A$. (b) $\det(A^k) = (\det A)^k$. (c) For any invertible $S$, $\det(S^{-1}AS) = \det A$. (d) $\det(cA) = c^n \det A$.

*Proof.* (a) $1 = \det I = \det(A A^{-1}) = \det A\det(A^{-1})$. (b) Induction on $k$ with the product rule. (c) $\det(S^{-1}AS) = \det(S^{-1})\det A\det S = \det A$ by (a). (d) $cA = (cI)A$ and $\det(cI) = c^n$ by applying (D3) to each of the $n$ rows of $cI$. ∎

Part (c) is the result LA-7 needs most. In LA-4, two matrices $A$ and $S^{-1}AS$ describe the **same linear map** in two bases. So the determinant is a property of the map, not of the coordinates. §3.9 gives a second, coordinate-free proof of the same fact.

**Theorem 3.11 (Transpose).** $\det A^\top = \det A$.

*Proof.*

1. By (3.3) applied to $A^\top$, whose $(i, j)$ entry is $a_{ji}$: $\det A^\top = \sum_\sigma \operatorname{sgn}(\sigma)\prod_i a_{\sigma(i)\,i}$.
2. Reindex the product by $m = \sigma(i)$, so $i = \sigma^{-1}(m)$: $\prod_i a_{\sigma(i)\,i} = \prod_m a_{m\,\sigma^{-1}(m)}$.
3. $\operatorname{sgn}(\sigma^{-1}) = \operatorname{sgn}(\sigma)$: an inversion of $\sigma$ is a pair of positions $i < j$ with values $\sigma(i) > \sigma(j)$; the same pair, read as positions $\sigma(j) < \sigma(i)$ of $\sigma^{-1}$ carrying values $j > i$, is an inversion of $\sigma^{-1}$. This is a bijection between the two inversion sets, so $N(\sigma^{-1}) = N(\sigma)$.
4. As $\sigma$ runs over all permutations so does $\sigma^{-1}$. So $\det A^\top = \sum_{\rho}\operatorname{sgn}(\rho)\prod_m a_{m\rho(m)} = \det A$ with $\rho = \sigma^{-1}$. ∎

■ **What just happened.** Every row statement in §3.2 is now a column statement too: the determinant is linear in each column, flips under a column swap, is zero when two columns are equal or the columns are dependent, and is unchanged by adding a multiple of one column to another. The signed-area function of §2 was a function of columns; Definition 3.1 was stated for rows; Theorem 3.11 says they are the same function.

### 3.5 Cofactor expansion

Group the terms of (3.3) according to which column the first row uses.

**Theorem 3.12 (Laplace expansion).** For any row $i$,
$$\det A = \sum_{j=1}^n a_{ij}\, C_{ij} = \sum_{j=1}^n (-1)^{i+j} a_{ij} M_{ij}, \tag{3.4}$$
where $M_{ij}$ is the minor (the determinant of $A$ with row $i$ and column $j$ deleted) and $C_{ij} = (-1)^{i+j}M_{ij}$ is the cofactor. The same holds down any column $j$: $\det A = \sum_{i=1}^n a_{ij}C_{ij}$.

▶ **Try it first.** In (3.3), collect the terms with $\sigma(1) = j$. Factor out $a_{1j}$. Argue that what remains is $\pm$ the determinant of the matrix with row 1 and column $j$ deleted, and find the sign by counting inversions involving position 1.

*Proof.*

1. *Row 1.* Split the sum (3.3) by the value $j = \sigma(1)$:
   $$\det A = \sum_{j=1}^n a_{1j}\sum_{\sigma:\,\sigma(1) = j}\operatorname{sgn}(\sigma)\,a_{2\sigma(2)}\cdots a_{n\sigma(n)}.$$
2. Fix $j$. A permutation with $\sigma(1) = j$ is determined by its restriction to positions $2, \dots, n$, which is a bijection $\{2, \dots, n\} \to \{1, \dots, n\}\setminus\{j\}$. Relabel the rows $2, \dots, n$ as $1, \dots, n-1$ and the columns $\{1, \dots, n\}\setminus\{j\}$ as $1, \dots, n-1$ in increasing order; this turns the restriction into a permutation $\tau$ of $\{1, \dots, n-1\}$, and $a_{2\sigma(2)}\cdots a_{n\sigma(n)}$ into the product $\prod_k b_{k\tau(k)}$ of entries of the $(n-1)\times(n-1)$ matrix $B$ obtained by deleting row 1 and column $j$. Every $\tau$ arises from exactly one $\sigma$.
3. *The sign.* Inversions of $\sigma$ split into those involving position 1 and those among positions $2, \dots, n$. Position 1 carries the value $j$, and it is inverted with every later position carrying a smaller value; the values $1, \dots, j-1$ all appear later, so there are exactly $j - 1$ such inversions. Inversions among positions $2, \dots, n$ are preserved by the order-preserving relabeling, so they are exactly the inversions of $\tau$. Hence $N(\sigma) = (j-1) + N(\tau)$ and $\operatorname{sgn}(\sigma) = (-1)^{j-1}\operatorname{sgn}(\tau) = (-1)^{1+j}\operatorname{sgn}(\tau)$.
4. Substituting into step 1: the inner sum is $(-1)^{1+j}\sum_\tau \operatorname{sgn}(\tau)\prod_k b_{k\tau(k)} = (-1)^{1+j}\det B = (-1)^{1+j}M_{1j}$. This is (3.4) for $i = 1$.
5. *Row $i$.* Move row $i$ to the top by $i - 1$ adjacent swaps (swap it with the row above, repeatedly). The determinant is multiplied by $(-1)^{i-1}$, the rows other than $i$ keep their relative order, and so the minors of the new first row are exactly the $M_{ij}$ of the original matrix. Step 4 applied to the new matrix gives $\det A = (-1)^{i-1}\sum_j (-1)^{1+j}a_{ij}M_{ij} = \sum_j(-1)^{i+j}a_{ij}M_{ij}$.
6. *Columns.* Apply the row result to $A^\top$ and use Theorem 3.11. ∎

The signs $(-1)^{i+j}$ form a checkerboard starting with $+$ in the top-left corner:
$$\begin{pmatrix} + & - & + & - \\ - & + & - & + \\ + & - & + & - \\ - & + & - & + \end{pmatrix}.$$

■ **What just happened.** The $n\times n$ determinant is a combination of $n$ determinants of size $n - 1$. Recursing all the way down costs about $e \cdot n!$ multiplications (Exercise A-level arithmetic: $T(n) = nT(n-1) + n$, so $T(n)/n! \to e$). Cofactor expansion is the right tool when a row has many zeros, when the matrix has symbolic entries, and for proofs (§3.6, §3.7, §3.12). It is the wrong tool for a numerical $10\times 10$ matrix; use (3.2).

### 3.6 Cramer's rule

**Theorem 3.13 (Cramer).** If $\det A \neq 0$, the unique solution of $A\mathbf{x} = \mathbf{b}$ has components
$$x_j = \frac{\det A_j}{\det A}, \qquad A_j = A \text{ with column } j \text{ replaced by } \mathbf{b}. \tag{3.5}$$

▶ **Try it first.** Column $j$ of $A_j$ is $\mathbf{b} = A\mathbf{x} = \sum_k x_k\mathbf{a}_k$. Expand $\det A_j$ using linearity in that column.

*Proof.*

1. $\det A_j = \det(\mathbf{a}_1, \dots, \mathbf{a}_{j-1}, \mathbf{b}, \mathbf{a}_{j+1}, \dots, \mathbf{a}_n)$, listing columns.
2. By the column picture (LA-1, (3.1)), $\mathbf{b} = A\mathbf{x} = \sum_{k=1}^n x_k\mathbf{a}_k$. Linearity in column $j$ (Theorem 3.11 plus (D3)) gives
   $$\det A_j = \sum_{k=1}^n x_k\det(\mathbf{a}_1, \dots, \mathbf{a}_{j-1}, \mathbf{a}_k, \mathbf{a}_{j+1}, \dots, \mathbf{a}_n).$$
3. For $k \neq j$ the matrix in the $k$-th term has column $\mathbf{a}_k$ twice, so its determinant is $0$. For $k = j$ it is $A$. So $\det A_j = x_j\det A$. Divide. ∎

■ **What just happened.** Cramer's rule is one line once you have column-linearity. It is a *formula* for the solution, which is exactly what you want when the entries are symbols (a $2\times 2$ with resistances as entries, Example 4.5) or when you want to prove something about how the solution depends on the data (every $x_j$ is a ratio of polynomials in the entries; it is continuous wherever $\det A \neq 0$). It is not an algorithm: solving an $n \times n$ system by Cramer costs $n + 1$ determinants, each $\sim \tfrac23 n^3$ by elimination, against $\sim \tfrac23 n^3$ total for eliminating once.

### 3.7 The inverse via the adjugate

**Theorem 3.14 (Adjugate formula).** Let $C$ be the matrix of cofactors, $C_{ij} = (-1)^{i+j}M_{ij}$, and $\operatorname{adj}A = C^\top$. Then
$$A \,\operatorname{adj}A = (\det A)\, I, \qquad\text{so}\qquad A^{-1} = \frac{\operatorname{adj}A}{\det A}\ \text{ when } \det A \neq 0. \tag{3.6}$$

▶ **Try it first.** The $(i, k)$ entry of $A\,C^\top$ is $\sum_j a_{ij}C_{kj}$. For $i = k$ this is (3.4). For $i \neq k$, which matrix's cofactor expansion is it?

*Proof.*

1. $(A\,C^\top)_{ik} = \sum_j a_{ij}(C^\top)_{jk} = \sum_j a_{ij}C_{kj}$.
2. *Diagonal, $i = k$:* $\sum_j a_{ij}C_{ij} = \det A$ by (3.4), expansion along row $i$.
3. *Off-diagonal, $i \neq k$:* Let $A'$ be $A$ with row $k$ replaced by a copy of row $i$. The cofactors $C_{kj}$ are computed from the matrix with row $k$ deleted, so they do not see what is in row $k$: $A'$ and $A$ have the same $C_{kj}$ for every $j$. Expanding $\det A'$ along its row $k$, whose entries are $a_{ij}$, gives $\det A' = \sum_j a_{ij}C_{kj}$, which is exactly our sum. But $A'$ has two equal rows, so $\det A' = 0$.
4. So $A\,C^\top$ has $\det A$ on the diagonal and $0$ elsewhere. If $\det A \neq 0$, divide by it. ∎

For $n = 2$: $\operatorname{adj}\begin{pmatrix} a & b \\ c & d\end{pmatrix} = \begin{pmatrix} d & -b \\ -c & a\end{pmatrix}$, giving the formula you memorized for the $2\times 2$ inverse. For $n = 3$ you need nine $2\times 2$ determinants; for $n = 4$, sixteen $3\times 3$ ones. Like Cramer, this is a formula for proofs and symbols (it shows, for instance, that the entries of $A^{-1}$ are polynomials in the entries of $A$ divided by $\det A$), not a way to compute inverses numerically. The project's "Your turn" asks you to measure how badly it does.

### 3.8 The cross product is a determinant

For $\mathbf{u}, \mathbf{v} \in \mathbb{R}^3$, define
$$\mathbf{u}\times\mathbf{v} = \det\begin{pmatrix} \mathbf{e}_1 & \mathbf{e}_2 & \mathbf{e}_3 \\ u_1 & u_2 & u_3 \\ v_1 & v_2 & v_3\end{pmatrix} = (u_2v_3 - u_3v_2)\,\mathbf{e}_1 - (u_1v_3 - u_3v_1)\,\mathbf{e}_2 + (u_1v_2 - u_2v_1)\,\mathbf{e}_3, \tag{3.7}$$
where the "determinant" with vectors in the first row is shorthand for the cofactor expansion along that row.

**Theorem 3.15 (Cross product facts).** For $\mathbf{u}, \mathbf{v}, \mathbf{w} \in \mathbb{R}^3$:
(a) $\mathbf{w}\cdot(\mathbf{u}\times\mathbf{v}) = \det\begin{pmatrix}\mathbf{w}^\top \\ \mathbf{u}^\top \\ \mathbf{v}^\top\end{pmatrix}$ (the scalar triple product is the determinant with those rows).
(b) $\mathbf{u}\times\mathbf{v}$ is orthogonal to both $\mathbf{u}$ and $\mathbf{v}$.
(c) $\|\mathbf{u}\times\mathbf{v}\|^2 = \|\mathbf{u}\|^2\|\mathbf{v}\|^2 - (\mathbf{u}\cdot\mathbf{v})^2$, hence $\|\mathbf{u}\times\mathbf{v}\| = \|\mathbf{u}\|\|\mathbf{v}\|\sin\theta$, the area of the parallelogram on $\mathbf{u}, \mathbf{v}$.

*Proof.* (a) Dotting (3.7) with $\mathbf{w}$ replaces $\mathbf{e}_k$ by $w_k$ in each cofactor term, which is the cofactor expansion of the stated determinant along its first row. (b) Put $\mathbf{w} = \mathbf{u}$ in (a): the matrix has two equal rows, so the triple product is $0$. Same for $\mathbf{w} = \mathbf{v}$. (c) Expand both sides:
$$\|\mathbf{u}\|^2\|\mathbf{v}\|^2 - (\mathbf{u}\cdot\mathbf{v})^2 = \sum_{i,j}u_i^2v_j^2 - \sum_{i,j}u_iv_iu_jv_j = \sum_{i < j}\left(u_i^2v_j^2 + u_j^2v_i^2 - 2u_iv_iu_jv_j\right) = \sum_{i<j}(u_iv_j - u_jv_i)^2,$$
where the $i = j$ terms cancelled between the two double sums. The three terms $i < j$ are the squares of the three components in (3.7). Then $\|\mathbf{u}\|^2\|\mathbf{v}\|^2(1 - \cos^2\theta) = \|\mathbf{u}\|^2\|\mathbf{v}\|^2\sin^2\theta$ by the angle formula (LA-1, Theorem 3.10). ∎

■ **What just happened.** Every cross-product identity you memorized in mechanics is a determinant property. Anticommutativity $\mathbf{v}\times\mathbf{u} = -\mathbf{u}\times\mathbf{v}$ is a row swap. $\mathbf{u}\times\mathbf{u} = \mathbf{0}$ is two equal rows. Linearity in each factor is (D3). And (a) with (c) says $|\mathbf{w}\cdot(\mathbf{u}\times\mathbf{v})|$ is base area times height, the volume of the parallelepiped on $\mathbf{u}, \mathbf{v}, \mathbf{w}$, which is the $n = 3$ case of §3.10. The cross product exists only in $\mathbb{R}^3$ because only there does deleting one row of an $n\times n$ determinant leave exactly one vector's worth of cofactors; the object that generalizes is the alternating form of §3.9.

### 3.9 The alternating-multilinear characterization (Axler's route)

Forget matrices for a moment. A function $f : \underbrace{\mathbb{R}^n \times \cdots \times \mathbb{R}^n}_{n \text{ copies}} \to \mathbb{R}$ is **multilinear** if it is linear in each argument when the others are held fixed, and **alternating** if $f(\mathbf{v}_1, \dots, \mathbf{v}_n) = 0$ whenever two of the arguments are equal. Signed area in §2 was such a function with $n = 2$. The determinant, viewed as a function of the $n$ columns, is such a function for every $n$. The theorem is that there are no others, up to a constant.

**Theorem 3.16 (Uniqueness of alternating multilinear forms).** If $f : (\mathbb{R}^n)^n \to \mathbb{R}$ is multilinear and alternating, then for all $\mathbf{v}_1, \dots, \mathbf{v}_n$,
$$f(\mathbf{v}_1, \dots, \mathbf{v}_n) = f(\mathbf{e}_1, \dots, \mathbf{e}_n)\cdot\det\begin{pmatrix}\mathbf{v}_1 & \cdots & \mathbf{v}_n\end{pmatrix}. \tag{3.8}$$
In words: the alternating multilinear functions of $n$ vectors in $\mathbb{R}^n$ form a one-dimensional space, spanned by $\det$.

▶ **Try it first.** First show that "alternating" forces "swapping two arguments flips the sign" by expanding $f(\dots, \mathbf{u} + \mathbf{v}, \dots, \mathbf{u} + \mathbf{v}, \dots) = 0$. Then repeat the expansion of §2 with $n$ arguments.

*Proof.*

1. *Alternating implies antisymmetric.* Put $\mathbf{u} + \mathbf{v}$ in two argument slots. By alternation the value is $0$; by multilinearity in those two slots it expands to
   $$0 = f(\dots\mathbf{u}\dots\mathbf{u}\dots) + f(\dots\mathbf{u}\dots\mathbf{v}\dots) + f(\dots\mathbf{v}\dots\mathbf{u}\dots) + f(\dots\mathbf{v}\dots\mathbf{v}\dots).$$
   The first and last terms are $0$ by alternation, so $f(\dots\mathbf{v}\dots\mathbf{u}\dots) = -f(\dots\mathbf{u}\dots\mathbf{v}\dots)$.
2. *Expand in the standard basis.* Write $\mathbf{v}_j = \sum_{i} v_{ij}\mathbf{e}_i$, so $v_{ij}$ is the $(i, j)$ entry of the matrix $V = (\mathbf{v}_1 \cdots \mathbf{v}_n)$. Multilinearity in all $n$ slots gives
   $$f(\mathbf{v}_1, \dots, \mathbf{v}_n) = \sum_{i_1, \dots, i_n} v_{i_1 1}v_{i_2 2}\cdots v_{i_n n}\, f(\mathbf{e}_{i_1}, \dots, \mathbf{e}_{i_n}).$$
3. Terms with a repeated index vanish by alternation. The rest are indexed by permutations $\sigma$ with $i_j = \sigma(j)$. Sorting the arguments of $f(\mathbf{e}_{\sigma(1)}, \dots, \mathbf{e}_{\sigma(n)})$ into standard order by swaps, each swap flipping the sign (step 1) and the parity of $N$ (Lemma 3.8), gives $f(\mathbf{e}_{\sigma(1)}, \dots, \mathbf{e}_{\sigma(n)}) = \operatorname{sgn}(\sigma)f(\mathbf{e}_1, \dots, \mathbf{e}_n)$.
4. Therefore $f(\mathbf{v}_1, \dots, \mathbf{v}_n) = f(\mathbf{e}_1, \dots, \mathbf{e}_n)\sum_\sigma \operatorname{sgn}(\sigma)\prod_j v_{\sigma(j)j}$. The sum is the Leibniz formula (3.3) for $\det V^\top$, which equals $\det V$ by Theorem 3.11. ∎

■ **What just happened.** The proof is the proof of Theorem 3.9 again, with "the matrix with rows $\mathbf{e}_{\sigma(i)}$" replaced by "the value of $f$ on the basis vectors in permuted order." The content is the same; the packaging is coordinate-free, and that buys two things.

**Corollary 3.16′ (The determinant of a linear map).** Let $T : \mathbb{R}^n \to \mathbb{R}^n$ be linear and let $f$ be any nonzero alternating multilinear form. Then $g(\mathbf{v}_1, \dots, \mathbf{v}_n) = f(T\mathbf{v}_1, \dots, T\mathbf{v}_n)$ is again alternating and multilinear, so by Theorem 3.16 there is a scalar $c$ with $g = c\,f$. The scalar $c$ does not depend on which $f$ you chose, and it equals $\det[T]_{\mathcal B}$ for the matrix of $T$ in **any** basis $\mathcal B$.

*Proof.* $g$ is multilinear because $T$ is linear and $f$ is multilinear; it is alternating because equal inputs to $T$ give equal outputs. Apply (3.8) to both $g$ and $f$: $g(\mathbf{v}) = g(\mathbf{e})\det V$ and $f(\mathbf{v}) = f(\mathbf{e})\det V$, so $c = g(\mathbf{e})/f(\mathbf{e}) = f(T\mathbf{e}_1, \dots, T\mathbf{e}_n)/f(\mathbf{e}_1, \dots, \mathbf{e}_n) = \det(T\mathbf{e}_1 \cdots T\mathbf{e}_n) = \det[T]$ in the standard basis, using (3.8) once more. For another basis, repeat the argument with $\mathbf{e}$ replaced by that basis: the ratio $g/f$ is the same scalar $c$ no matter which basis you evaluate it on, and evaluating on $\mathcal B$ gives $\det[T]_{\mathcal B}$. ∎

So $\det T$ is defined without choosing coordinates, and Corollary 3.10′(c) is the shadow of this fact. The second thing it buys is a one-line product rule: $f(AB\mathbf{v}_1, \dots, AB\mathbf{v}_n) = (\det A) f(B\mathbf{v}_1, \dots, B\mathbf{v}_n) = (\det A)(\det B) f(\mathbf{v}_1, \dots, \mathbf{v}_n)$, so $\det(AB) = \det A\det B$.

This is also the door to differential forms. In EM-1 and LA-12 you will meet $dx \wedge dy$, the alternating form that eats two vectors in the plane and returns the signed area of their parallelogram. It is $\det$. The Jacobian factor in §3.11 is the statement that a smooth map pulls back this form to $(\det J)\,du\wedge dv$. And in QM-7, the antisymmetric wavefunction of $n$ fermions is a determinant (the Slater determinant) precisely because "two particles in the same state gives zero" is the alternating property.

### 3.10 Volume scaling and orientation

To prove that $|\det|$ is a volume we need to say what volume is. We take the following facts about $n$-dimensional volume as known from multivariable calculus (they are theorems of measure theory, which is where their rigorous proofs live):

- (V1) The unit cube has volume $1$; a box with edges $|d_1|, \dots, |d_n|$ along the axes has volume $|d_1\cdots d_n|$.
- (V2) Volume is unchanged by translation and by permuting the coordinate axes.
- (V3) A set contained in a subspace of dimension less than $n$ has volume $0$.
- (V4) *Cavalieri / shear.* Replacing an edge $\mathbf{v}_i$ of a parallelepiped by $\mathbf{v}_i + c\,\mathbf{v}_j$ ($j \neq i$) does not change its volume: every cross-section parallel to the face spanned by the other edges is translated, not changed.
- (V5) Volume is additive on finite unions of sets that overlap only on volume-zero sets, and monotone under inclusion.

The **parallelepiped** spanned by $\mathbf{v}_1, \dots, \mathbf{v}_n$ is $P(\mathbf{v}_1, \dots, \mathbf{v}_n) = \{t_1\mathbf{v}_1 + \cdots + t_n\mathbf{v}_n : 0 \le t_i \le 1\}$.

**Theorem 3.17 (Volume and determinant).** Let $A$ be $n\times n$ and $T(\mathbf{x}) = A\mathbf{x}$.
(a) $\operatorname{vol}P(\mathbf{v}_1, \dots, \mathbf{v}_n) = |\det(\mathbf{v}_1 \cdots \mathbf{v}_n)|$.
(b) For every parallelepiped $P$, $\operatorname{vol}T(P) = |\det A|\operatorname{vol}P$.
(c) For every bounded region $\Omega$ that has a volume, $\operatorname{vol}T(\Omega) = |\det A|\operatorname{vol}\Omega$. \tag{3.9}

▶ **Try it first.** For (a), run elimination on the **columns** of $V = (\mathbf{v}_1\cdots\mathbf{v}_n)$ using only swaps and "add a multiple of one column to another." Both operations leave $|\det V|$ alone (Theorems 3.11, 3.3, (D2)) and leave $\operatorname{vol}P$ alone (V2, V4). Where do you end up?

*Proof.*

1. *(a), singular case.* If $\det V = 0$ the columns are dependent (Theorem 3.5 via Theorem 3.11), so $P$ lies in their span, a subspace of dimension $< n$, and has volume $0$ by (V3). Both sides are $0$.
2. *(a), invertible case.* Column operations of two kinds, swapping two columns and adding a multiple of one column to another, change neither $|\det V|$ nor $\operatorname{vol}P(\mathbf{v}_1, \dots, \mathbf{v}_n)$ (the swap does not change the set $P$ at all; the addition is (V4)). Gauss–Jordan elimination on the columns of an invertible matrix (LA-2, applied to $V^\top$) reduces $V$ to a diagonal matrix $D = \operatorname{diag}(d_1, \dots, d_n)$ with all $d_k \neq 0$ using only these two operations: swap to bring a nonzero entry to the pivot position, then clear its row to the left and right. (Clearing the row of column $k$ uses column $k$ to modify other columns, which is allowed.)
3. $P(d_1\mathbf{e}_1, \dots, d_n\mathbf{e}_n)$ is a box with edge lengths $|d_k|$, so its volume is $|d_1\cdots d_n|$ by (V1), and $|\det D| = |d_1\cdots d_n|$ by Theorem 3.6. Since neither quantity changed during the reduction, $\operatorname{vol}P(\mathbf{v}_1, \dots, \mathbf{v}_n) = |\det V|$.
4. *(b).* $T$ is linear, so $T(P(\mathbf{v}_1, \dots, \mathbf{v}_n)) = \{\sum t_i A\mathbf{v}_i\} = P(A\mathbf{v}_1, \dots, A\mathbf{v}_n)$. By (a) and the product rule, $\operatorname{vol}T(P) = |\det(AV)| = |\det A|\,|\det V| = |\det A|\operatorname{vol}P$. A translated parallelepiped maps to a translated parallelepiped, so (V2) extends this to parallelepipeds not based at the origin.
5. *(c), singular case.* If $\det A = 0$, $T(\Omega)$ lies in the column space of $A$, a proper subspace, and has volume $0$ by (V3).
6. *(c), invertible case.* "$\Omega$ has a volume" means: for every $\varepsilon > 0$ there is a grid of cubes of side $h$ such that the union $I_h$ of the grid cubes contained in $\Omega$ and the union $O_h$ of the grid cubes meeting $\Omega$ satisfy $I_h \subseteq \Omega \subseteq O_h$ and $\operatorname{vol}O_h - \operatorname{vol}I_h < \varepsilon$. (This is Jordan measurability, the hypothesis under which Riemann integrals of volumes make sense.) Since $T$ is a bijection, $T(I_h) \subseteq T(\Omega) \subseteq T(O_h)$, and $T$ carries the disjoint cubes of $I_h$ to parallelepipeds overlapping only on faces. By step 4 and (V5),
   $$\operatorname{vol}T(I_h) = |\det A|\operatorname{vol}I_h \le \operatorname{vol}T(\Omega) \le \operatorname{vol}T(O_h) = |\det A|\operatorname{vol}O_h.$$
   Both outer quantities are within $|\det A|\,\varepsilon$ of $|\det A|\operatorname{vol}\Omega$. Let $\varepsilon \to 0$. ∎

■ **What just happened.** The determinant was defined by three algebraic properties on rows. Theorem 3.17 says its absolute value is the volume-scaling factor of the linear map, for every region, not only for parallelepipeds. The approximation step is the same one that makes Riemann integration work; it is where "the linear map scales each tiny cube by $|\det A|$" becomes "it scales everything by $|\det A|$."

**Orientation.** The sign of $\det A$ is what is left. Say $A$ (invertible) **preserves orientation** if $\det A > 0$ and **reverses** it if $\det A < 0$. Rotations in $\mathbb{R}^2$ have determinant $\cos^2\theta + \sin^2\theta = 1$; the reflection $\operatorname{diag}(1, -1)$ has determinant $-1$. The two classes cannot be connected: if $A(t)$ is a continuous path of invertible matrices from $A(0)$ to $A(1)$, then $\det A(t)$ is a continuous function (it is a polynomial in the entries) that never vanishes, so by the intermediate value theorem it cannot change sign. A reflection cannot be reached from the identity by a continuous motion through invertible maps. That is the precise content of "you cannot turn a left hand into a right hand by rotating it."

### 3.11 The Jacobian determinant

You have changed variables in multiple integrals since calculus: $dx\,dy = r\,dr\,d\theta$ in polar coordinates, $dx\,dy\,dz = r^2\sin\phi\,dr\,d\phi\,d\theta$ in spherical. The general statement is:

**Theorem 3.18 (Change of variables).** Let $\mathbf{g} : U \to \mathbb{R}^n$ be continuously differentiable and injective on an open set $U \subseteq \mathbb{R}^n$, with Jacobian matrix $J_{\mathbf{g}}(\mathbf{u}) = \left(\partial g_i/\partial u_j\right)$ invertible on $U$. Then for integrable $f$,
$$\int_{\mathbf{g}(U)} f(\mathbf{x})\,d\mathbf{x} = \int_U f(\mathbf{g}(\mathbf{u}))\,\bigl|\det J_{\mathbf{g}}(\mathbf{u})\bigr|\,d\mathbf{u}. \tag{3.10}$$

*Why it is Theorem 3.17.* Near a point $\mathbf{u}_0$, the differentiable map is its linearization plus a small error: $\mathbf{g}(\mathbf{u}_0 + \mathbf{h}) = \mathbf{g}(\mathbf{u}_0) + J_{\mathbf{g}}(\mathbf{u}_0)\mathbf{h} + o(\|\mathbf{h}\|)$. So a small cube $Q$ of side $h$ at $\mathbf{u}_0$ is mapped to, approximately, a translate of the parallelepiped $J_{\mathbf{g}}(\mathbf{u}_0)\,Q$, whose volume is $|\det J_{\mathbf{g}}(\mathbf{u}_0)|\,h^n$ by (3.9). Summing $f(\mathbf{g}(\mathbf{u}_0))\times(\text{volume of the image cell})$ over a grid of cubes is a Riemann sum for the right side of (3.10); summing $f(\mathbf{x})\times(\text{volume of the cell})$ over the image cells is a Riemann sum for the left side. The two sums agree up to the $o(h^n)$ errors per cell, which total $o(1)$.

What is skipped here, and is the entire difficulty of the real proof, is controlling the error uniformly: showing the image of a small cube is contained in a slightly inflated parallelepiped and contains a slightly deflated one. That argument is in Spivak, *Calculus on Manifolds*, Theorem 3-13, and Rudin, *Principles of Mathematical Analysis*, Theorem 10.9. The linear-algebra content of the theorem is exactly (3.9), and the $r$ in $r\,dr\,d\theta$ is exactly a $2\times 2$ determinant (Example 4.4).

The absolute value is there because volume is positive. The sign of $\det J$ still carries information: it says whether $\mathbf{g}$ preserves or reverses orientation at that point, and when you integrate differential forms (EM-1, Stokes' theorem) rather than densities, the sign stays.

### 3.12 The Vandermonde determinant and interpolation

Given $n$ distinct numbers $x_1, \dots, x_n$, the **Vandermonde matrix** is
$$V = \begin{pmatrix} 1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\ 1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\ \vdots & & & & \vdots \\ 1 & x_n & x_n^2 & \cdots & x_n^{n-1}\end{pmatrix}.$$
Row $i$ is the monomials evaluated at $x_i$. So $V\mathbf{c}$, for $\mathbf{c} = (c_0, \dots, c_{n-1})^\top$, is the vector of values $p(x_i)$ of the polynomial $p(x) = c_0 + c_1x + \cdots + c_{n-1}x^{n-1}$. Finding the polynomial of degree $\le n - 1$ through the points $(x_i, y_i)$ is solving $V\mathbf{c} = \mathbf{y}$. That is Cramer's 1750 problem.

**Theorem 3.19 (Vandermonde).**
$$\det V = \prod_{1 \le i < j \le n}(x_j - x_i). \tag{3.11}$$
In particular $\det V \neq 0$ iff the $x_i$ are distinct, and then the interpolating polynomial exists and is unique.

▶ **Try it first.** Think of $\det V$ as a polynomial in $x_n$ alone (the other $x_i$ are constants). What is its degree? Where must it vanish? What does that force?

*Proof.* Induction on $n$. For $n = 1$, $\det V = 1$ and the empty product is $1$.

1. Expand $\det V$ along the last row (3.4). Each term is $\pm x_n^{j-1}$ times a minor that does not involve $x_n$. So $\det V = p(x_n)$ for a polynomial $p$ of degree $\le n - 1$ whose coefficient of $x_n^{n-1}$ is $(-1)^{n+n}M_{nn} = M_{nn}$, the determinant of the $(n-1)\times(n-1)$ Vandermonde matrix on $x_1, \dots, x_{n-1}$.
2. If $x_n$ is set equal to any of $x_1, \dots, x_{n-1}$, two rows of $V$ become equal, so $p(x_i) = 0$ for $i = 1, \dots, n-1$. A polynomial of degree $\le n - 1$ with the $n - 1$ distinct roots $x_1, \dots, x_{n-1}$ (distinct by hypothesis) is a constant times $\prod_{i < n}(x_n - x_i)$, and matching the coefficient of $x_n^{n-1}$ identifies the constant: $p(x_n) = M_{nn}\prod_{i<n}(x_n - x_i)$.
3. By induction $M_{nn} = \prod_{1\le i<j\le n-1}(x_j - x_i)$. Multiplying by the new factors $\prod_{i < n}(x_n - x_i)$ gives (3.11). ∎

■ **What just happened.** A determinant identity was proved by treating the determinant as a polynomial and counting roots. The same technique proves the characteristic polynomial facts in LA-7. And note the two ways distinctness enters: if two nodes coincide, two rows coincide, the determinant is $0$, and interpolation fails (two different $y$ values at the same $x$ is impossible; the same $y$ value gives a one-parameter family of solutions). Exercise B4 extracts a theorem about polynomials from this.

### 3.13 What the determinant is, numerically, and what it is not

**It is a polynomial.** From (3.3), $\det A$ is a homogeneous polynomial of degree $n$ in the $n^2$ entries. Consequences: it is continuous (used in §3.10 for orientation and in LA-7 for the continuity of eigenvalues); it is differentiable (Exercise C1 computes the derivative); and $\det(cA) = c^n\det A$.

**Cost.** Cofactor recursion: $T(n) = nT(n-1) + O(n)$, so $T(n) \approx e\cdot n!$. For $n = 20$ that is about $6.6 \times 10^{18}$ multiplications, centuries on a laptop. Elimination: about $\tfrac23 n^3$ flops, $5{,}300$ for $n = 20$. The project measures both; the cofactor curve goes vertical at $n \approx 9$. Every library computes determinants by LU factorization, which is (3.2) with partial pivoting.

**Scale.** Because $\det(cA) = c^n\det A$, the determinant says nothing about how close $A$ is to singular. The matrix $\tfrac{1}{10}I_{100}$ has determinant $10^{-100}$ and is as far from singular as a matrix can be: its inverse is $10 I$, perfectly behaved. The matrix $\operatorname{diag}(10^6, 10^{-6})$ has determinant $1$ and is numerically treacherous: perturbing its $(2,2)$ entry by $10^{-6}$ makes it singular. The right measure of "how singular" is the condition number $\kappa(A) = \sigma_{\max}/\sigma_{\min}$, the ratio of the largest to smallest singular values, defined in LA-9 and used throughout LA-13. The determinant is $\pm\prod_i\sigma_i$, a product that one tiny factor can dominate or a hundred moderate factors can shrink to nothing, and neither tells you about the ratio.

**Overflow.** For $n$ in the hundreds, $\det A$ routinely exceeds $10^{308}$ or falls below $10^{-308}$ and overflows double precision even when every entry is of order $1$. NumPy's `slogdet` returns the sign and $\log|\det A|$ separately for this reason; use it whenever the determinant is a quantity you actually need, as in the Gaussian log-likelihood of F5 and ML-6.

**What it is for.** Three things, none of them "deciding whether a matrix is invertible numerically." (1) Theory: it is the one scalar invariant of a linear map that detects invertibility and multiplies under composition, so it is how eigenvalues are defined (LA-7) and how volume enters integration (§3.11). (2) Formulas: Cramer, the adjugate, Vandermonde, the cross product, when entries are symbols. (3) Small or structured matrices: a $2\times 2$ or $3\times 3$ by hand, a triangular or block matrix by inspection.

## 4. Worked examples

### Example 4.1: One $3\times 3$, three ways

$$A = \begin{pmatrix} 0 & 2 & 1 \\ 1 & 1 & 3 \\ 2 & -1 & 1 \end{pmatrix}.$$

*By elimination, (3.2).* The $(1,1)$ entry is $0$, so swap rows 1 and 2 ($s = 1$):
$$\begin{pmatrix} 1 & 1 & 3 \\ 0 & 2 & 1 \\ 2 & -1 & 1\end{pmatrix} \xrightarrow{R_3 \leftarrow R_3 - 2R_1} \begin{pmatrix} 1 & 1 & 3 \\ 0 & 2 & 1 \\ 0 & -3 & -5\end{pmatrix} \xrightarrow{R_3 \leftarrow R_3 + \tfrac32 R_2} \begin{pmatrix} 1 & 1 & 3 \\ 0 & 2 & 1 \\ 0 & 0 & -\tfrac72\end{pmatrix}.$$
Pivots $1, 2, -\tfrac72$, product $-7$, one swap: $\det A = (-1)^1(-7) = 7$.

*By cofactors along row 1, (3.4).* $a_{11} = 0$ kills the first term.
$$\det A = 0\cdot C_{11} + 2\cdot(-1)^{1+2}\begin{vmatrix} 1 & 3 \\ 2 & 1\end{vmatrix} + 1\cdot(-1)^{1+3}\begin{vmatrix} 1 & 1 \\ 2 & -1\end{vmatrix} = -2(1 - 6) + (-1 - 2) = 10 - 3 = 7.$$

*By cofactors down column 1.* $\det A = 0\cdot C_{11} + 1\cdot C_{21} + 2\cdot C_{31}$ with $C_{21} = (-1)^{2+1}\begin{vmatrix}2 & 1\\-1 & 1\end{vmatrix} = -(2 + 1) = -3$ and $C_{31} = (-1)^{3+1}\begin{vmatrix}2 & 1 \\ 1 & 3\end{vmatrix} = 6 - 1 = 5$. Total $-3 + 10 = 7$.

*By the Leibniz formula, (3.3),* all six terms: $(1,2,3)$: $+0\cdot1\cdot1 = 0$; $(1,3,2)$: $-0\cdot3\cdot(-1) = 0$; $(2,1,3)$: $-2\cdot1\cdot1 = -2$; $(2,3,1)$: $+2\cdot3\cdot2 = 12$; $(3,1,2)$: $+1\cdot1\cdot(-1) = -1$; $(3,2,1)$: $-1\cdot1\cdot2 = -2$. Sum: $7$.

Same answer four times, as uniqueness guarantees. The project checks all four.

### Example 4.2: The unit square under shear-then-rotate

Shear by $2$ along $x$, then rotate by $30°$:
$$S = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}, \qquad R = \begin{pmatrix} \cos 30° & -\sin 30° \\ \sin 30° & \cos 30°\end{pmatrix} = \begin{pmatrix} \tfrac{\sqrt3}{2} & -\tfrac12 \\ \tfrac12 & \tfrac{\sqrt3}{2}\end{pmatrix}, \qquad M = RS = \begin{pmatrix} \tfrac{\sqrt3}{2} & \sqrt3 - \tfrac12 \\ \tfrac12 & 1 + \tfrac{\sqrt3}{2}\end{pmatrix} \approx \begin{pmatrix} 0.866 & 1.232 \\ 0.500 & 1.866\end{pmatrix}.$$
The image of the unit square is the parallelogram with edges $M\mathbf{e}_1 \approx (0.866, 0.500)^\top$ and $M\mathbf{e}_2 \approx (1.232, 1.866)^\top$. Its signed area is
$$\det M = \tfrac{\sqrt3}{2}\left(1 + \tfrac{\sqrt3}{2}\right) - \tfrac12\left(\sqrt3 - \tfrac12\right) = \tfrac{\sqrt3}{2} + \tfrac34 - \tfrac{\sqrt3}{2} + \tfrac14 = 1.$$
Numerically $0.866 \times 1.866 - 1.232 \times 0.500 = 1.616 - 0.616 = 1.000$. You could have skipped the arithmetic: $\det M = \det R\det S = 1\cdot 1$, because a shear has determinant $1\cdot 1 - 2\cdot 0 = 1$ and a rotation has $\cos^2 + \sin^2 = 1$. Shearing and rotating move area around without creating or destroying it, and the product rule says so.

Now reflect across the $x$-axis as a third step, $F = \operatorname{diag}(1, -1)$: $\det(FM) = -1$. Same area, opposite orientation: going from $FM\mathbf{e}_1 \approx (0.866, -0.5)^\top$ to $FM\mathbf{e}_2 \approx (1.232, -1.866)^\top$ is a clockwise turn. And if you stretch first by $D = \operatorname{diag}(3, \tfrac12)$, the image of the square under $RSD$ has area $|\det R\det S\det D| = \tfrac32$. The project draws all of these.

### Example 4.3: A $3\times 3$ Vandermonde, computed, factored, and used

Nodes $x = (1, 2, 4)$:
$$V = \begin{pmatrix} 1 & 1 & 1 \\ 1 & 2 & 4 \\ 1 & 4 & 16\end{pmatrix}.$$
*Computed by elimination:* $R_2 \leftarrow R_2 - R_1 = (0, 1, 3)$, $R_3 \leftarrow R_3 - R_1 = (0, 3, 15)$, $R_3 \leftarrow R_3 - 3R_2 = (0, 0, 6)$. Pivots $1, 1, 6$, no swaps: $\det V = 6$.

*Factored by (3.11):* $(x_2 - x_1)(x_3 - x_1)(x_3 - x_2) = (1)(3)(2) = 6$. ✓

*Factored by hand, in symbols,* to see where (3.11) comes from: subtract row 1 from rows 2 and 3 (no change in determinant), then expand down column 1:
$$\det V = \begin{vmatrix} 1 & x_1 & x_1^2 \\ 0 & x_2 - x_1 & x_2^2 - x_1^2 \\ 0 & x_3 - x_1 & x_3^2 - x_1^2\end{vmatrix} = \begin{vmatrix} x_2 - x_1 & (x_2 - x_1)(x_2 + x_1) \\ x_3 - x_1 & (x_3 - x_1)(x_3 + x_1)\end{vmatrix} = (x_2 - x_1)(x_3 - x_1)\begin{vmatrix} 1 & x_2 + x_1 \\ 1 & x_3 + x_1\end{vmatrix},$$
pulling the common factor out of each row by (D3). The last determinant is $(x_3 + x_1) - (x_2 + x_1) = x_3 - x_2$. ✓

*Used for interpolation.* Which parabola passes through $(1, 2), (2, 3), (4, 11)$? Solve $V\mathbf{c} = (2, 3, 11)^\top$. By Cramer (3.5), the $x^2$ coefficient is
$$c_2 = \frac{1}{6}\begin{vmatrix} 1 & 1 & 2 \\ 1 & 2 & 3 \\ 1 & 4 & 11\end{vmatrix} = \frac{1}{6}\cdot 6 = 1$$
(eliminate: rows become $(1,1,2), (0,1,1), (0,3,9) \to (0,0,6)$). Likewise $c_1 = \tfrac{1}{6}\begin{vmatrix} 1 & 2 & 1 \\ 1 & 3 & 4 \\ 1 & 11 & 16\end{vmatrix} = \tfrac{-12}{6} = -2$ and $c_0 = \tfrac{1}{6}\begin{vmatrix} 2 & 1 & 1 \\ 3 & 2 & 4 \\ 11 & 4 & 16\end{vmatrix} = \tfrac{18}{6} = 3$. The parabola is $p(x) = x^2 - 2x + 3$; check $p(1) = 2$, $p(2) = 3$, $p(4) = 11$. ✓ Had two nodes coincided, $\det V$ would be $0$ and Cramer would divide by zero: no parabola, or infinitely many.

### Example 4.4: Polar coordinates, or why it is $r\,dr\,d\theta$

The map $\mathbf{g}(r, \theta) = (r\cos\theta, r\sin\theta)$ has Jacobian matrix
$$J = \begin{pmatrix} \partial x/\partial r & \partial x/\partial\theta \\ \partial y/\partial r & \partial y/\partial\theta\end{pmatrix} = \begin{pmatrix} \cos\theta & -r\sin\theta \\ \sin\theta & r\cos\theta\end{pmatrix}, \qquad \det J = r\cos^2\theta + r\sin^2\theta = r.$$
Read the columns geometrically. The first column, $\partial\mathbf{g}/\partial r = (\cos\theta, \sin\theta)^\top$, is a unit vector pointing radially: moving $dr$ in the $r$ direction moves you $dr$ outward. The second, $\partial\mathbf{g}/\partial\theta = r(-\sin\theta, \cos\theta)^\top$, has length $r$ and points tangentially: moving $d\theta$ moves you $r\,d\theta$ along the circle. The two columns are orthogonal (dot product $-r\cos\theta\sin\theta + r\sin\theta\cos\theta = 0$), so the parallelogram they span is a rectangle with sides $1$ and $r$, area $r$. The small coordinate cell $[r, r + dr]\times[\theta, \theta + d\theta]$ in the $(r, \theta)$ plane, of area $dr\,d\theta$, is carried to a near-rectangle of area $r\,dr\,d\theta$. That is (3.10) with $|\det J| = r$. The sign is positive: $(\partial_r\mathbf{g}, \partial_\theta\mathbf{g})$ is counterclockwise, so polar coordinates preserve orientation (for $r > 0$).

A check with numbers: the quarter annulus $1 \le r \le 2$, $0 \le \theta \le \pi/2$ has area $\int_0^{\pi/2}\int_1^2 r\,dr\,d\theta = \tfrac{\pi}{2}\cdot\tfrac{4 - 1}{2} = \tfrac{3\pi}{4} \approx 2.356$. Forgetting the Jacobian gives $\tfrac{\pi}{2}\cdot 1 \approx 1.571$, wrong by a factor of $\tfrac32$, the average value of $r$ on the region. The project estimates the area by Monte Carlo and recovers $3\pi/4$ only with the factor $r$ included. Exercise B6 does the same for spherical coordinates and produces $r^2\sin\phi$.

### Example 4.5: Cramer's rule in a two-loop circuit, against elimination

Two loops share a resistor $R_2$; a battery $V$ drives loop 1. Kirchhoff's voltage law around each loop, with loop currents $I_1, I_2$:
$$\begin{pmatrix} R_1 + R_2 & -R_2 \\ -R_2 & R_2 + R_3\end{pmatrix}\begin{pmatrix} I_1 \\ I_2\end{pmatrix} = \begin{pmatrix} V \\ 0\end{pmatrix}.$$
*Cramer with symbols.* $\det = (R_1 + R_2)(R_2 + R_3) - R_2^2 = R_1R_2 + R_1R_3 + R_2R_3$. Then
$$I_1 = \frac{\begin{vmatrix} V & -R_2 \\ 0 & R_2 + R_3\end{vmatrix}}{\det} = \frac{V(R_2 + R_3)}{R_1R_2 + R_1R_3 + R_2R_3}, \qquad I_2 = \frac{\begin{vmatrix} R_1 + R_2 & V \\ -R_2 & 0\end{vmatrix}}{\det} = \frac{VR_2}{R_1R_2 + R_1R_3 + R_2R_3}.$$
This is what you want from Cramer: closed forms showing how each current depends on each resistance. You can read off that $I_2 \to 0$ as $R_3 \to \infty$ (open the second loop) and that $\det > 0$ whenever all resistances are positive, so the circuit equations are never singular. No amount of elimination with numbers would tell you that.

*Numbers.* $R_1 = R_2 = R_3 = 1\,\Omega$, $V = 3\,\text{V}$: $\det = 3$, $I_1 = 3\cdot 2/3 = 2\,\text{A}$, $I_2 = 3\cdot 1/3 = 1\,\text{A}$.

*Elimination with the same numbers.* $\begin{pmatrix} 2 & -1 & | & 3 \\ -1 & 2 & | & 0\end{pmatrix} \xrightarrow{R_2 \leftarrow R_2 + \frac12 R_1} \begin{pmatrix} 2 & -1 & | & 3 \\ 0 & \tfrac32 & | & \tfrac32\end{pmatrix}$, so $I_2 = 1$ and $2I_1 - 1 = 3$ gives $I_1 = 2$. Same answer, fewer operations: elimination on a $2\times 2$ costs one multiplier and a back-substitution; Cramer costs three $2\times 2$ determinants. For $n = 2$ nobody cares. For $n = 1000$ (a thousand-node circuit, a finite-element mesh), Cramer by elimination-computed determinants costs $1001$ eliminations instead of one, and Cramer by cofactors costs $1001!$ operations, which is not a number of operations anyone has performed.

## 5. Why it matters and where it connects

- **LA-7.** Eigenvalues are the $\lambda$ with $A - \lambda I$ singular, that is, the roots of the **characteristic polynomial** $p(\lambda) = \det(A - \lambda I)$ (Theorem 3.7). It is a polynomial of degree $n$ by (3.3) (Exercise C2), and by Corollary 3.10′(c) it does not depend on the basis, so eigenvalues are a property of the map. Also $\det A = \prod_i\lambda_i$ and the Wronskian test for independent solutions of a linear ODE system is a determinant.
- **LA-8.** Orthogonal matrices have $\det Q = \pm 1$ (from $Q^\top Q = I$ and Theorems 3.10, 3.11; Exercise B2). Sylvester's criterion tests positive definiteness through the signs of the leading principal minors.
- **LA-9 and LA-13.** $|\det A| = \prod_i\sigma_i$, the product of the singular values, while the condition number is their ratio $\sigma_{\max}/\sigma_{\min}$. §3.13 is the reason the latter and not the former measures nearness to singularity. Determinants are computed via the LU factorization, which is (3.2).
- **LA-12.** Theorem 3.16 says the space of alternating $n$-linear forms on an $n$-dimensional space is one-dimensional. LA-12 builds the whole exterior algebra ($k$-forms for $k < n$, wedge products) on this fact; $\det$ is the top-degree case, and $dx\wedge dy\wedge dz$ is the volume form.
- **EM-1.** Curvilinear coordinates and the Jacobian; the cross product in $\nabla\times$; the volume element in Gauss's law computations.
- **CM-7, CM-8.** The inertia tensor's principal axes come from LA-8, but the volume of the body is a triple integral with a Jacobian. Liouville's theorem (CM-8) says Hamiltonian flow preserves phase-space volume: the Jacobian determinant of the flow map is identically $1$, and Exercise C1's formula $\tfrac{d}{dt}\det = \det\cdot\operatorname{tr}(\ldots)$ is how one proves it.
- **F5, ST-1, ML-6.** The multivariate Gaussian density has $(\det\Sigma)^{-1/2}$ in front because the map $\mathbf{z} \mapsto \Sigma^{1/2}\mathbf{z}$ scales volume by $(\det\Sigma)^{1/2}$ (Theorem 3.17). Transforming a density under $\mathbf{y} = \mathbf{g}(\mathbf{x})$ multiplies by $|\det J|^{-1}$ (3.10). Gaussian mixture log-likelihoods carry $\log\det\Sigma_k$; compute it with `slogdet`.
- **QM-7.** The Slater determinant: the antisymmetric wavefunction of $n$ identical fermions is $\det[\psi_i(\mathbf{x}_j)]$. The Pauli exclusion principle is Theorem 3.2: two fermions in the same single-particle state means two equal rows, and the wavefunction vanishes.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 drill the three ways to compute, the row-operation rules, permutation signs, Cramer and the adjugate on $2\times 2$, and the cross product. B1–B8 are the module: B1 (block triangular determinants), B2 (similarity invariance and orthogonal matrices), B3 (recognizing an alternating multilinear form), B5 (tetrahedron volume), B6 (the spherical Jacobian), and B8 (the `project.py` extension) are essential. C1 (Jacobi's formula for the derivative of a determinant), C2 (the characteristic polynomial, a direct preview of LA-7), and C3 (Hadamard's inequality and the Gram determinant, with code) go one step beyond; do at least one.

## 7. Coding project

`project.py` implements the determinant three ways and tests every claim in the notes that has a number in it:

1. **Three algorithms.** `det_elimination` (partial pivoting, counts swaps, applies (3.2)), `det_leibniz` (all $n!$ permutations with signs from inversion counts, (3.3)), and `det_cofactor` (recursive (3.4)). All three are asserted against `np.linalg.det` on random matrices, and all four computations of Example 4.1 are reproduced, including the six Leibniz terms and the pivots $1, 2, -\tfrac72$.
2. **Cost.** Times cofactor against elimination for $n = 2, \dots, 9$ and plots both on a log axis. You will watch factorial growth happen.
3. **Volume.** Maps $100{,}000$ uniform points in the unit square and unit cube through random matrices and estimates the image's volume two ways (bounding-box rejection with the inverse map, and `scipy.spatial.ConvexHull`), checking against $|\det A|$ to within $1\%$ (Theorem 3.17).
4. **Pictures.** The unit square and its images under Example 4.2's shear, rotation, and reflection, with orientation shown by the order of the edge vectors.
5. **Identities.** $\det(AB) = \det A\det B$, $\det A^\top = \det A$, $\det A^{-1} = 1/\det A$, $\det(cA) = c^n\det A$, the adjugate identity (3.6), Cramer's rule, and the Vandermonde formula (3.11), numerically on random matrices and symbolically in SymPy. The polar and spherical Jacobians are computed symbolically and the quarter-annulus area $3\pi/4$ is recovered by Monte Carlo with the factor $r$ (Example 4.4).
6. **Scale.** `diag(1e6, 1e-6)` has determinant $1$ and condition number $10^{12}$; $0.1 I_{100}$ has determinant $10^{-100}$ and condition number $1$.

**Build it yourself first:** write `det_elimination(A)` (eliminate with partial pivoting, count swaps, multiply the diagonal) and `det_leibniz(A)` (`itertools.permutations`, sign from counting inversions) before reading the code. Then the "Your turn" cells: the adjugate inverse versus `np.linalg.inv` on an ill-conditioned $4\times 4$, the volume of a tetrahedron from its vertices, and a repeated-row matrix under random perturbation of its third row.

## 8. Where you'll struggle

- **Which row operation does what.** Three operations, three different effects: swapping two rows multiplies $\det$ by $-1$; multiplying a row by $c$ multiplies $\det$ by $c$; adding a multiple of one row to another does **nothing**. The swap sign is the one people lose. When you eliminate by hand, write "$s = $" in the margin and increment it every time you swap; do not try to remember afterwards. Example 4.1 has exactly one swap and the sign matters: without it you get $-7$.
- **$\det(A + B) \neq \det A + \det B$ and $\det(cA) \neq c\det A$.** Linearity is in **one row at a time**. $\det(2I_3) = 8$, not $2$. If you catch yourself writing $\det(A + B)$ as a sum, go back to (3.1) and note that it has one row changing and $n - 1$ rows fixed.
- **The checkerboard.** The cofactor $C_{ij}$ has sign $(-1)^{i+j}$, so $C_{12}$ is **negative** and $C_{21}$ is negative, while $C_{13}$, $C_{22}$, $C_{31}$ are positive. Expanding along row 2 starts with a minus. Write the sign pattern down every time until you stop needing to. The other classic slip is using the minor $M_{ij}$ where the cofactor $C_{ij}$ is needed, which drops all the signs; the adjugate formula (3.6) is wrong without them.
- **A small determinant does not mean nearly singular.** It means nothing without knowing the scale. $\det(0.1 I_{100}) = 10^{-100}$. The quantity that measures nearness to singularity is the condition number $\kappa$, defined in LA-9 and used in LA-13; until then, the habit to build is "never test singularity by comparing $\det$ to a tolerance." The project's Part 6 shows both failure modes with numbers.
- **The alternating-multilinear characterization feels like abstraction for its own sake.** It is not. It is what makes the determinant **coordinate-free**: Corollary 3.16′ defines $\det T$ for a linear map with no basis in sight, and that is what LA-7 needs when it says "the eigenvalues of $T$" rather than "the eigenvalues of this matrix." It also gives the product rule in one line and generalizes to $k$ vectors in $\mathbb{R}^n$ (LA-12, where $\det$ becomes one member of a family). If §3.9 feels empty, reread the §2 derivation of $a_1b_2 - a_2b_1$ and notice that it was the $n = 2$ case of Theorem 3.16, done before you had any formula.
- **Cofactors and Cramer are formulas, not algorithms.** Both are the right tool for symbols (Example 4.5) and proofs (§3.7, §3.12), and the wrong tool for numbers past $n = 3$. If you are reaching for cofactor expansion on a numerical $5\times 5$, stop and eliminate.
- **The volume theorem's approximation step** (Theorem 3.17, step 6) is the first place in this course where "has a volume" is a hypothesis rather than a given. If it bothers you, that is appropriate; the full story is measure theory, and the version here is exactly what Riemann integration assumes.
- **B1 is harder than it looks.** The natural attempt (expand along rows, track signs) drowns in indices. The clean proof uses uniqueness: show that the block determinant, as a function of the top-left block, satisfies (D1)–(D3). C1 needs you to trust the uniqueness theorem for linear ODEs; if that is rusty, do C2 instead.

## 9. Check yourself

1. State the three defining properties of the determinant. Which one is used to prove that a matrix with two equal rows has determinant zero, and how?
2. Elimination on a $4\times 4$ produced pivots $2, -1, 3, \tfrac12$ with two row swaps. What is the determinant?
3. $\det A = 5$ for a $3\times 3$ matrix. What are $\det(2A)$, $\det(A^\top)$, $\det(A^{-1})$, $\det(A^2)$, and $\det(-A)$?
4. Write the Leibniz formula for $n = 3$ with all six terms and their signs. Which term corresponds to the permutation $(3, 1, 2)$?
5. Why is $\det(AB) = \det A\det B$? Give the uniqueness-based proof in three sentences.
6. What is the sign of the cofactor $C_{23}$? Of $C_{32}$? Expand a general $3\times 3$ determinant along its second row.
7. Derive Cramer's rule from linearity in a column, in two lines.
8. Why is $\mathbf{u}\times\mathbf{v}$ perpendicular to $\mathbf{u}$? Answer with a property of determinants, not with a computation.
9. State Theorem 3.16. Why does it mean the determinant of a linear map is independent of the basis?
10. A linear map has determinant $-3$. What happens to the area of a disk of radius $1$ under it? To the orientation of a counterclockwise loop?
11. Where does the $r$ in $r\,dr\,d\theta$ come from? Compute it.
12. $\det A = 10^{-15}$. Is $A$ nearly singular? What would you need to know to answer?

## 10. Pair with

- **Strang, *Introduction to Linear Algebra*, ch. 5.** §5.1 is the three properties and their consequences in the same order as §3.2 here; §5.2 derives the big formula and cofactors; §5.3 does Cramer, the adjugate, and volume. Strang is better than these notes at the $3\times 3$ pictures of the volume argument.
- **MIT 18.06, lectures 18, 19, 20 (Strang) on OCW.** Lecture 18 is the properties (watch him derive the triangular case live); 19 is the formula and cofactors; 20 is Cramer, the inverse formula, and volume. Lecture 19's derivation of the big formula by expanding rows is the one §3.3 follows.
- **Axler, *Linear Algebra Done Right*, 4th ed., ch. 9** (open access). §9B defines alternating multilinear forms and proves the one-dimensionality (our Theorem 3.16) in full generality over any field; §9C defines the determinant of an operator the way Corollary 3.16′ does and proves the product rule from it. Axler is better than these notes at the coordinate-free side and at the connection to the characteristic polynomial. Read §9B after §3.9 here.
- **3Blue1Brown, *Essence of Linear Algebra*, episode 6 ("The determinant").** Ten minutes of animated area scaling and orientation flips. Watch it after §2, before §3.
- **Spivak, *Calculus on Manifolds*, Theorem 3-13**, if you want the full proof of the change-of-variables formula (3.10). Not required.
