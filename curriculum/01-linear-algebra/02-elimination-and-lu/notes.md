# LA-2: Solving Systems: Elimination, Echelon Forms, LU, Inverses, Rank

**Time:** ~8 hours · **Prerequisites:** LA-1 · **Unlocks:** LA-3 (and through it the whole LA core; LA-6 uses the pivots for determinants, LA-10 and LA-13 return to LU numerically)

> **Conventions in this module.** $A$ is $m\times n$ with entries $a_{ij}$ (row $i$, column $j$), rows $\mathbf{r}_1^\top,\dots,\mathbf{r}_m^\top$ and columns $\mathbf{a}_1,\dots,\mathbf{a}_n$, as in LA-1. Unless stated otherwise $A$ is square, $n\times n$. $\mathbf{e}_i$ is the $i$-th standard basis vector and $\mathbf{e}_i\mathbf{e}_j^\top$ is the matrix with a single $1$ in position $(i,j)$. The basic row operation is **"row $i$ ← row $i$ − $\ell$ · row $j$"**; $\ell$ is the **multiplier**, written $\ell_{ij}$: the first index is the row being changed, the second is the pivot row being subtracted. When every pivot sits on the diagonal (square invertible $A$), $\ell_{ij}$ is also the multiplier that clears position $(i,j)$; after a free column the pivot of row $j$ lies right of the diagonal, and $\ell_{ij}$ clears row $i$ in that pivot's column (Example 4.2). Sign convention: $L$ stores $+\ell_{ij}$; the elementary matrix $E_{ij}$ that performs the operation contains $-\ell_{ij}$. **Pivots** are the first nonzero entries of the rows of an echelon form; for square invertible $A$ they are the diagonal of $U$. A permutation matrix $P$ has rows $\mathbf{e}_{\pi(1)}^\top, \dots, \mathbf{e}_{\pi(n)}^\top$ for a bijection $\pi$ of $\{1, \dots, n\}$, so row $i$ of $PA$ is row $\pi(i)$ of $A$ and $P\mathbf{e}_j = \mathbf{e}_{\pi^{-1}(j)}$ (the project's `np.eye(n)[perm]`). Permutations are written $PA = LU$ (SciPy's `lu` returns $A = PLU$, the transposed convention; the project checks the translation). A **flop** is one floating-point add, subtract, multiply, or divide. A pivot is "zero" numerically if $|u_{kk}| \le \texttt{tol}$ for a stated tolerance (so $\texttt{tol} = 0$ flags exact zeros only).

## 0. Where this sits

LA-1 told you *when* $A\mathbf{x}=\mathbf{b}$ has a solution (when $\mathbf{b}$ is in the span of the columns) but gave you no way to find one. This module is the algorithm: Gaussian elimination, and the factorization $A = LU$ that elimination secretly computes. You need LA-1's column picture, the row form of $A\mathbf{x}$ (LA-1 eq. 3.2), the column rule for $AB$ (LA-1 eq. 3.3), and $(AB)^\top = B^\top A^\top$ (LA-1 Theorem 3.8).

What it makes possible: LA-3 defines the null space as the solution set of $A\mathbf{x}=\mathbf{0}$ and the rank as a dimension; both are read off the echelon forms you compute here. LA-6 shows the determinant is $\pm$ the product of the pivots. LA-10 and LA-13 ask how $LU$ behaves in floating point; this module gives the first, decisive example of why that question matters. Every later project that calls `np.linalg.solve` is calling this module.

## 1. The problem that forced this idea

Elimination is the oldest algorithm in this curriculum. Chapter 8 of the Chinese *Nine Chapters on the Mathematical Art* (compiled roughly two thousand years ago) is titled *Fangcheng*, "rectangular arrays." It lays coefficients out on a counting board and subtracts multiples of one column from another until the array is triangular, then substitutes back. That is Gaussian elimination, minus the name, minus the matrices.

The name came much later. In 1809 and 1810 Gauss was fitting the orbit of the asteroid Pallas to many observations by least squares. The normal equations of least squares (derived in LA-5) are a square symmetric system with as many unknowns as orbital elements, six or more, and Gauss needed a *procedure* a human computer could follow without thinking: eliminate the unknowns one at a time, in a fixed order, with tight bookkeeping. Geodesists and astronomers adopted it as "Gauss's method," and the surveyor Wilhelm Jordan published, in the 1888 edition of his *Handbuch der Vermessungskunde*, the variant that also eliminates above the pivots so the system ends as the identity: Gauss–Jordan. (Joseph Grcar's 2011 article "How ordinary elimination became Gaussian elimination," *Historia Mathematica*, traces the naming.)

Two twentieth-century events turned the recipe into mathematics. First, people noticed it is a **factorization**: record what elimination does to $A$ and the record is a lower-triangular $L$, the result is an upper-triangular $U$, and $A = LU$. Alan Turing's 1948 paper "Rounding-off errors in matrix processes" is where the name $LU$ and the clean matrix statement appear (Banachiewicz in 1938 and Dwyer in 1944 had matrix forms of the same idea). Once elimination is a factorization you can ask when it exists, when it is unique, what it costs, and what it costs to reuse.

Second, the first electronic computers made people afraid of it. Harold Hotelling argued in 1943 that rounding errors could grow like $4^n$, which would make a $100\times 100$ system hopeless. Von Neumann and Goldstine (1947) showed the fear was unfounded for symmetric positive definite matrices, and James Wilkinson's backward-error analysis in the early 1960s explained why elimination *with row swaps* works in practice even though a pathological worst case exists. That vocabulary (growth factor, backward stability, conditioning) is LA-13. Here you will meet the single example that makes the subject necessary: a $2\times 2$ matrix on which elimination without swaps returns an answer wrong in the first digit.

So the questions: given $A$ and $\mathbf{b}$, how do I find $\mathbf{x}$; how do I know in advance whether there is one, many, or none; how do I do it a thousand times for a thousand right-hand sides without paying a thousand times; and when does the arithmetic lie to me.

## 2. Intuition

### Subtracting equations doesn't change the answer

Take $2x + y = 5$ and $4x + 5y = 13$. Subtract $2$ times the first from the second: $3y = 3$. The new pair "$2x + y = 5$, $3y = 3$" has exactly the same solutions as the old pair. Anything satisfying both old equations satisfies the new one (it is built from them), and the operation is reversible: add $2$ times the first equation back and you recover the old second one.

In the row picture, you replaced the line $4x + 5y = 13$ by the line $3y = 3$. Both pass through the solution $(2, 1)$: you rotated one line about the solution until it was horizontal. Elimination is a sequence of such rotations; every hyperplane moves, but all keep passing through the common intersection.

The new system is **triangular**: the last equation has one unknown. Read it off ($y = 1$), substitute upward ($x = 2$). That is **back substitution**. The whole strategy: *make it triangular, then substitute.*

### The staircase

For $3\times 3$: clear below the first entry of column 1 (two operations), then below the second entry of column 2 (one operation):

```
[ ■  *  * ]
[ 0  ■  * ]      ■ = pivot (nonzero)
[ 0  0  ■ ]      * = anything
```

For a non-square or degenerate matrix the staircase can have wide steps:

```
[ ■  *  *  *  * ]
[ 0  0  ■  *  * ]      columns 1, 3, 4 have pivots
[ 0  0  0  ■  * ]      columns 2, 5 do not ("free")
[ 0  0  0  0  0 ]      a zero row
```

That is a **row echelon form**. The number of pivots will be the most important number attached to a matrix.

### Where it breaks, and what breaking means

Elimination wants to divide by the pivot. Two things can go wrong.

*Case 1: the pivot position holds a zero but something below it doesn't.* Swap two rows (swapping two equations changes nothing) and continue. A nuisance, not a defect.

*Case 2: the pivot position and everything below it are zero.* No swap helps. The column contributes no step to the staircase: it is free. This is information: that column is a combination of the columns before it. One of LA-1's ingredients was redundant.

Push this to the end and the three outcomes of LA-1 exercise B1 (two of which LA-1 Example 4.2 showed) get a mechanism. Eliminate on $[A\,|\,\mathbf{b}]$. A row reading $0\ 0\ 0\,|\,c$ with $c\ne 0$ says $0 = c$: no solution. Otherwise there is a solution, unique exactly when every column of $A$ got a pivot; a free column means its unknown can be anything, so infinitely many.

### The record of the operations is a factorization

Guess this before §3 proves it. Each step is "row $i$ ← row $i$ − $\ell_{ij}$ · row $j$." At the end you have $U$ and the multipliers. Run the steps backwards: row $i$ of $A$ equals row $i$ of $U$ plus $\ell_{ij}$ times row $j$ of $U$ for each $j$ you subtracted. In matrix language, $A = LU$ with $L$ unit lower triangular holding the $\ell_{ij}$. **The multipliers are the factorization.**

Why care? Solving $A\mathbf{x}=\mathbf{b}$ becomes two triangular solves, $L\mathbf{c}=\mathbf{b}$ then $U\mathbf{x}=\mathbf{c}$. The expensive part depends only on $A$. A thousand right-hand sides (a thousand load cases on one structure, a thousand steps of an implicit ODE solver): factor once, substitute a thousand times.

### Guess the cost

About $n$ steps, each updating roughly $n^2$ entries with a multiply and a subtract: about $n^3$ flops, and the constant will turn out to be $2/3$. A triangular solve touches each entry of $U$ once: about $n^2$. A factorization costs about $n$ substitutions. Remember that ratio when someone suggests computing $A^{-1}$.

### When an inverse exists

$A^{-1}$ undoes $A$: $A^{-1}A = I$. Since $AA^{-1} = I$ says $A$ sends column $j$ of $A^{-1}$ to $\mathbf{e}_j$, the columns of $A^{-1}$ are the solutions of $A\mathbf{x} = \mathbf{e}_j$: $n$ systems with the same $A$, which is why Gauss–Jordan on $[A\,|\,I]$ works. If elimination finds $n$ pivots, every $A\mathbf{x}=\mathbf{b}$ has exactly one solution and those columns exist. If it finds fewer, some $\mathbf{x}\ne\mathbf{0}$ has $A\mathbf{x}=\mathbf{0}$, and then $A^{-1}$ would have to send $\mathbf{0}$ to both $\mathbf{0}$ and $\mathbf{x}$. So: **$A$ is invertible iff elimination produces $n$ pivots.** That is the main theorem.

## 3. Formal development

### 3.1 Row operations are left multiplications

Three row operations:

- **(R1) Replacement.** Row $i$ ← row $i$ − $\ell$ · row $j$, $i \ne j$.
- **(R2) Swap.** Exchange rows $i$ and $j$.
- **(R3) Scaling.** Row $i$ ← $c$ · row $i$, $c \ne 0$.

Each is a multiplication *on the left* by a square matrix, because of a one-line consequence of the row-times-column rule.

**Theorem 3.1 (Left multiplication acts on rows; right multiplication acts on columns).** For $E$ ($m\times m$) and $A$ ($m\times n$), row $i$ of $EA$ is $\sum_k E_{ik}\,(\text{row } k \text{ of } A)$. For $F$ ($n\times n$), column $j$ of $AF$ is $\sum_k F_{kj}\,(\text{column } k \text{ of } A)$.

*Proof.* $(EA)_{ij} = \sum_k E_{ik}A_{kj}$ (LA-1 eq. 3.3); fix $i$ and let $j$ run: that is the row vector $\sum_k E_{ik}\mathbf{r}_k^\top$. The column statement is LA-1 eq. 3.1 applied to each column of $F$. ∎

**Definition 3.2 (Elementary matrices).** For $i\ne j$ and scalars $\ell$, $c\ne 0$:

$$
E_{ij}(\ell) = I - \ell\,\mathbf{e}_i\mathbf{e}_j^\top, \qquad
P_{ij} = I - \mathbf{e}_i\mathbf{e}_i^\top - \mathbf{e}_j\mathbf{e}_j^\top + \mathbf{e}_i\mathbf{e}_j^\top + \mathbf{e}_j\mathbf{e}_i^\top, \qquad
D_i(c) = I + (c-1)\,\mathbf{e}_i\mathbf{e}_i^\top. \tag{3.1}
$$

Read them with Theorem 3.1. Row $i$ of $E_{ij}(\ell)$ is $\mathbf{e}_i^\top - \ell\mathbf{e}_j^\top$, so row $i$ of $E_{ij}(\ell)A$ is (row $i$) − $\ell$(row $j$), and all other rows of $E_{ij}(\ell)$ are rows of $I$, so all other rows of $A$ are untouched: (R1). $P_{ij}$ is $I$ with rows $i, j$ exchanged, so $P_{ij}A$ exchanges rows $i, j$ of $A$: (R2). $D_i(c)$ is $I$ with $c$ in position $(i,i)$: (R3). For $n = 3$:

$$
E_{21}(2) = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \qquad
P_{23} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}, \qquad
D_2(5) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 5 & 0 \\ 0 & 0 & 1 \end{pmatrix}.
$$

Note the sign: $E_{21}(2)$ *subtracts* $2\times$ row 1 from row 2 and stores $-2$. This is the most common bookkeeping error in the module (§8).

▶ **Try it first.** Each row operation can be undone by one of the same type. Guess the inverse of each matrix in (3.1), then verify by multiplying. You will need the product rule for the unit matrices $\mathbf{e}_i\mathbf{e}_j^\top$.

**Derivation (inverses of elementary matrices).**

1. For any indices, $(\mathbf{e}_i\mathbf{e}_j^\top)(\mathbf{e}_k\mathbf{e}_l^\top) = \mathbf{e}_i(\mathbf{e}_j^\top\mathbf{e}_k)\mathbf{e}_l^\top = \delta_{jk}\,\mathbf{e}_i\mathbf{e}_l^\top$, where $\delta_{jk} = 1$ if $j = k$ and $0$ otherwise, because $\mathbf{e}_j^\top\mathbf{e}_k = \mathbf{e}_j\cdot\mathbf{e}_k$. $\tag{3.2}$
2. $E_{ij}(\ell)E_{ij}(-\ell) = (I - \ell\mathbf{e}_i\mathbf{e}_j^\top)(I + \ell\mathbf{e}_i\mathbf{e}_j^\top) = I - \ell^2(\mathbf{e}_i\mathbf{e}_j^\top)(\mathbf{e}_i\mathbf{e}_j^\top) = I - \ell^2\delta_{ji}\,\mathbf{e}_i\mathbf{e}_j^\top = I$, since $i\ne j$; the same in the other order. So $E_{ij}(\ell)^{-1} = E_{ij}(-\ell)$: *to undo "subtract $\ell$ times row $j$," add $\ell$ times row $j$.*
3. $P_{ij}P_{ij}$ swaps rows $i, j$ of $P_{ij}$, giving back $I$. So $P_{ij}^{-1} = P_{ij}$.
4. $D_i(c)D_i(1/c)$ scales row $i$ of $I$ by $c$ then $1/c$. So $D_i(c)^{-1} = D_i(1/c)$. ∎

■ **What just happened.** Every row operation is invertible, with an inverse of the same type. That is the only property of row operations the rest of the module uses.

**Theorem 3.3 (Row operations preserve the solution set).** If $E$ is invertible, then $A\mathbf{x} = \mathbf{b}$ if and only if $(EA)\mathbf{x} = E\mathbf{b}$.

*Proof.* Multiply on the left by $E$ for one direction and by $E^{-1}$ for the other. ∎

A product of invertible matrices is invertible, with the inverses in reverse order: $(E_2E_1)(E_1^{-1}E_2^{-1}) = E_2(E_1E_1^{-1})E_2^{-1} = I$, and likewise $(E_1^{-1}E_2^{-1})(E_2E_1) = I$. So any sequence of row operations applied to $[A\,|\,\mathbf{b}]$ produces a system with exactly the same solutions.

### 3.2 Echelon forms and substitution

**Definition 3.4 (Row echelon form, pivots).** $U$ is in **row echelon form** if (i) every all-zero row is below every nonzero row, and (ii) the first nonzero entry of each nonzero row, its **pivot**, lies strictly to the right of the pivot of the row above. A column containing a pivot is a **pivot column**; the others are **free columns**. $U$ is in **reduced row echelon form (RREF)** if in addition every pivot equals $1$ and is the only nonzero entry in its column.

**Algorithm 3.5 (Elimination).** Start with pivot row $k = 1$. For each column $j = 1, \dots, n$: look at the entries of column $j$ in rows $k, \dots, m$. If all are zero, column $j$ is free; move to the next column. Otherwise choose a nonzero entry (by hand: any; in floating point: the largest in absolute value, §3.6), swap its row into row $k$, and for each $i > k$ set $\ell_{ik} = a_{ij}/a_{kj}$ (the entry of row $i$ in column $j$ divided by the pivot $a_{kj}$, both from the current matrix) and apply (R1): row $i$ ← row $i$ − $\ell_{ik}$ · row $k$. Column $j$ now has its pivot in row $k$ and zeros below. Set $k \leftarrow k+1$. Stop when the columns or rows run out. Until the first free column, $j = k$ and the pivot is the diagonal entry $a_{kk}$; after one, $j > k$. In column 3 of Example 4.2 the pivot is $a_{23}$ and $\ell_{32} = a_{33}/a_{23}$.

The output is in echelon form: each pivot is in a later column than the previous one, and every row below the last pivot row is zero. By Theorem 3.3, $U = E_s\cdots E_1 A$ and $[U\,|\,\mathbf{c}]$ with $\mathbf{c} = E_s\cdots E_1\mathbf{b}$ has the same solutions as $[A\,|\,\mathbf{b}]$.

For square $U$ with all $n$ pivots on the diagonal, **back substitution** solves $U\mathbf{x} = \mathbf{c}$:

$$
x_n = \frac{c_n}{u_{nn}}, \qquad x_i = \frac{1}{u_{ii}}\Big(c_i - \sum_{j=i+1}^{n} u_{ij}x_j\Big), \quad i = n-1, \dots, 1. \tag{3.3}
$$

For unit lower triangular $L$ (ones on the diagonal), **forward substitution** solves $L\mathbf{c} = \mathbf{b}$:

$$
c_1 = b_1, \qquad c_i = b_i - \sum_{j=1}^{i-1}\ell_{ij}c_j, \quad i = 2, \dots, n. \tag{3.4}
$$

**Theorem 3.6 (Triangular systems).** If $U$ is upper triangular with nonzero diagonal, $U\mathbf{x} = \mathbf{c}$ has exactly one solution, given by (3.3). Likewise $L\mathbf{c} = \mathbf{b}$ for unit lower triangular $L$, given by (3.4).

*Proof.* Equation $i$ of $U\mathbf{x}=\mathbf{c}$ involves only $x_i, \dots, x_n$. Equation $n$ forces $x_n = c_n/u_{nn}$. Given $x_n$, equation $n-1$ forces $x_{n-1}$, and so on upward; at each stage the later unknowns are fixed and the equation is solved for $x_i$ by dividing by $u_{ii}\ne 0$. Any solution must take these values (uniqueness) and these values satisfy every equation (existence). The lower-triangular case is the same from the top down. ∎

### 3.3 Existence and uniqueness, read off the echelon form

**Theorem 3.7 (The three outcomes).** Let $A$ be $m\times n$ and $[U\,|\,\mathbf{c}]$ any echelon form of $[A\,|\,\mathbf{b}]$, with $r$ pivots in the $U$ part.

(a) $A\mathbf{x} = \mathbf{b}$ has **no solution** iff some row of $[U\,|\,\mathbf{c}]$ is zero in the $U$ part and nonzero in the $\mathbf{c}$ part (the augmented column is a pivot column).

(b) If (a) does not occur and $r = n$ (every column of $A$ is a pivot column): **exactly one** solution.

(c) If (a) does not occur and $r < n$: **infinitely many**. The $n - r$ free variables may be given any values; the pivot variables are then determined uniquely. Every solution has the form $\mathbf{x}_p + \sum_{f\ \mathrm{free}} t_f\,\mathbf{s}_f$, where $\mathbf{x}_p$ is the solution with all free variables $0$ and the **special solution** $\mathbf{s}_f$ solves $U\mathbf{x} = \mathbf{0}$ with $x_f = 1$ and the other free variables $0$.

▶ **Try it first.** Convince yourself of (a) in one sentence, then *construct* a solution when (a) doesn't happen: choose the free variables, then work up through the pivot rows.

**Derivation.**

1. By Theorem 3.3 we may work with $[U\,|\,\mathbf{c}]$.
2. A row $(0\ \cdots\ 0\,|\,c_i)$ with $c_i \ne 0$ is the equation $0 = c_i$; nothing satisfies it. That is the "if" of (a).
3. Suppose no such row exists. Assign arbitrary values to the free variables. Read the pivot rows from the bottom up: pivot row $k$ has pivot $u_{k,p_k}\ne 0$ in column $p_k$ and zeros to its left, so it reads $u_{k,p_k}x_{p_k} + (\text{terms in } x_j,\ j > p_k) = c_k$. Every $x_j$ with $j > p_k$ is free (assigned) or belongs to a lower pivot row (already determined). Solve for $x_{p_k}$. Zero rows read $0 = 0$. We have a solution: the "only if" of (a).
4. In (b) there are no free variables, so step 3 involves no choices and every solution must equal the one built. Unique.
5. In (c) pick a free $x_f$; different values give solutions differing in coordinate $f$, so infinitely many. For the structure: if $\mathbf{x}$ solves the system, $\mathbf{x} - \mathbf{x}_p$ solves $U\mathbf{x} = \mathbf{0}$ with free coordinates $t_f = x_f$. By step 3 applied to $U\mathbf{x}=\mathbf{0}$, a homogeneous solution is determined by its free coordinates, and $\sum_f t_f\mathbf{s}_f$ is a homogeneous solution (linearity, LA-1 Theorem 3.4) with those same free coordinates. So $\mathbf{x} - \mathbf{x}_p = \sum_f t_f\mathbf{s}_f$. ∎

■ **What just happened.** LA-1's three outcomes now have a mechanism. (a) is "$\mathbf{b}$ is not in the span of the columns": elimination exposed a combination of the equations reading $0 = c$. (b) is "no column is redundant." (c) is "some column is a combination of earlier ones," and $\mathbf{s}_f$ writes down which. In LA-3 the set $\{\sum t_f\mathbf{s}_f\}$ becomes the **null space** $N(A)$ and $n - r$ its dimension.

### 3.4 The theorem this module exists to prove

**Theorem 3.8 (Fewer than $n$ pivots gives a nonzero null vector).** If $A$ is $n\times n$ and elimination produces fewer than $n$ pivots, then $A\mathbf{x} = \mathbf{0}$ has a solution $\mathbf{x} \ne \mathbf{0}$.

▶ **Try it first.** If column $f$ of $U$ has no pivot, set $x_f = 1$ and the other free variables to $0$. Can you always fill in the pivot variables so that $U\mathbf{x} = \mathbf{0}$? Why does that give $A\mathbf{x} = \mathbf{0}$ and not just $U\mathbf{x} = \mathbf{0}$?

*Proof.*

1. Elimination gives $EA = U$ with $E$ a product of elementary matrices, hence invertible (§3.1).
2. $U$ has $r < n$ pivots in $n$ columns, so some column $f$ is free.
3. $U\mathbf{x} = \mathbf{0}$ is consistent (outcome (a) cannot occur with a zero right side). By Theorem 3.7(c), set $x_f = 1$, the other free variables $0$, and back-substitute for the pivot variables. The result satisfies $U\mathbf{x} = \mathbf{0}$ and $x_f = 1$, so $\mathbf{x} \ne \mathbf{0}$.
4. By Theorem 3.3, $A\mathbf{x} = E^{-1}U\mathbf{x} = \mathbf{0}$. ∎

**Corollary 3.9 (More unknowns than equations).** If $A$ is $m\times n$ with $n > m$, then $A\mathbf{x} = \mathbf{0}$ has a nonzero solution: there are at most $m$ pivots, so some column is free, and the same construction works. (Exercise B3 asks for every step.)

### 3.5 $A = LU$: the record of elimination

Assume for now that elimination on square $A$ never needs a row swap (§3.6 removes this). The operations are $E_{ij} := E_{ij}(\ell_{ij})$, $i > j$, applied column by column:

$$
E_{n,n-1}\cdots E_{32}\,E_{n1}\cdots E_{31}E_{21}\,A = U,
\qquad\text{so}\qquad
A = E_{21}^{-1}E_{31}^{-1}\cdots E_{n1}^{-1}\,E_{32}^{-1}\cdots E_{n,n-1}^{-1}\,U =: LU. \tag{3.5}
$$

**Theorem 3.10 ($L$ stores the multipliers).** The matrix $L$ in (3.5) is unit lower triangular with $L_{ij} = \ell_{ij}$ for $i > j$: the entry in position $(i,j)$ is exactly the multiplier used to clear position $(i,j)$. No arithmetic beyond the elimination itself is needed to form $L$.

▶ **Try it first.** For $3\times 3$: row 1 of $U$ is row 1 of $A$ (never changed). What happened to row 3 of $A$ on its way to becoming row 3 of $U$? Write row 3 of $A$ as a combination of rows of $U$.

**Derivation (first way: run elimination backwards).**

1. Row $j$ is never changed after it serves as the pivot row (we only modify rows below the pivot row), so row $j$ of $U$ is "the pivot row $j$," call it $\mathbf{u}_j^\top$.
2. Row $i$ of $A$ is changed exactly $i-1$ times: at step $j = 1, \dots, i-1$ we subtract $\ell_{ij}\mathbf{u}_j^\top$. After the last of these it is $\mathbf{u}_i^\top$:
   $$\mathbf{u}_i^\top = (\text{row } i \text{ of } A) - \ell_{i1}\mathbf{u}_1^\top - \cdots - \ell_{i,i-1}\mathbf{u}_{i-1}^\top.$$
3. Rearranging: $(\text{row } i \text{ of } A) = \ell_{i1}\mathbf{u}_1^\top + \cdots + \ell_{i,i-1}\mathbf{u}_{i-1}^\top + 1\cdot\mathbf{u}_i^\top$.
4. By Theorem 3.1, this says $A = LU$ with row $i$ of $L$ equal to $(\ell_{i1}, \dots, \ell_{i,i-1}, 1, 0, \dots, 0)$: unit lower triangular with the multipliers below the diagonal. ∎

**Derivation (second way: multiply the inverses).** This explains why the *inverse* product in (3.5) is clean while the forward product $E_{32}E_{31}E_{21}$ is not.

1. $E_{ij}^{-1} = I + \ell_{ij}\mathbf{e}_i\mathbf{e}_j^\top$ (§3.1, step 2).
2. For $n = 3$: $(I + \ell_{21}\mathbf{e}_2\mathbf{e}_1^\top)(I + \ell_{31}\mathbf{e}_3\mathbf{e}_1^\top)$ has cross term $\ell_{21}\ell_{31}\delta_{13}\mathbf{e}_2\mathbf{e}_1^\top = 0$ by (3.2). Multiplying the result by $(I + \ell_{32}\mathbf{e}_3\mathbf{e}_2^\top)$ gives cross terms with $\delta_{13}$ again, both zero. So $L = I + \ell_{21}\mathbf{e}_2\mathbf{e}_1^\top + \ell_{31}\mathbf{e}_3\mathbf{e}_1^\top + \ell_{32}\mathbf{e}_3\mathbf{e}_2^\top$: the multipliers, overlaid.
3. In general, build the product (3.5) from left to right, one factor at a time. The factors come column by column, $q = 1, 2, \dots$, so every factor $(I + \ell_{ab}\mathbf{e}_a\mathbf{e}_b^\top)$ already multiplied in has $b \le q$ when the factor $(I + \ell_{pq}\mathbf{e}_p\mathbf{e}_q^\top)$, $p > q$, arrives. Inductive hypothesis: the product so far is $I + \sum \ell_{ab}\mathbf{e}_a\mathbf{e}_b^\top$, summed over the pairs $(a,b)$ already processed (true for the first factor). Multiplying by the new factor gives, by (3.2),
   $$\Big(I + \sum \ell_{ab}\mathbf{e}_a\mathbf{e}_b^\top\Big)(I + \ell_{pq}\mathbf{e}_p\mathbf{e}_q^\top) = I + \sum \ell_{ab}\mathbf{e}_a\mathbf{e}_b^\top + \ell_{pq}\mathbf{e}_p\mathbf{e}_q^\top + \sum \ell_{ab}\ell_{pq}\,\delta_{bp}\,\mathbf{e}_a\mathbf{e}_q^\top,$$
   and every cross term vanishes because $b \le q < p$ gives $\delta_{bp} = 0$. The hypothesis holds with $(p,q)$ added. After the last factor, $L = I + \sum_{i>j}\ell_{ij}\mathbf{e}_i\mathbf{e}_j^\top$. ∎

*Contrast.* In $E_{32}E_{31}E_{21}$ the factor $(I - \ell_{32}\mathbf{e}_3\mathbf{e}_2^\top)$ stands left of $(I - \ell_{21}\mathbf{e}_2\mathbf{e}_1^\top)$ and the cross term $\ell_{32}\ell_{21}\delta_{22}\mathbf{e}_3\mathbf{e}_1^\top$ survives: the $(3,1)$ entry of $E_{32}E_{31}E_{21}$ is $-\ell_{31} + \ell_{32}\ell_{21}$, not $-\ell_{31}$. Example 4.1 computes both.

■ **What just happened.** $A = LU$ is elimination *remembered*. In code, overwrite each zeroed $a_{ij}$ with $\ell_{ij}$ and one array holds $L$ (strictly below the diagonal) and $U$ (on and above) at no extra cost. Solving $A\mathbf{x}=\mathbf{b}$ becomes $L\mathbf{c}=\mathbf{b}$ by (3.4) then $U\mathbf{x}=\mathbf{c}$ by (3.3). Forward substitution *is* elimination applied to $\mathbf{b}$: $c_i = b_i - \sum_{j<i}\ell_{ij}c_j$ is what happens to $b_i$ when you subtract $\ell_{ij}$ times the processed $j$-th equation.

Two facts about triangular matrices are needed for uniqueness and for §3.10.

**Lemma 3.11 (Triangular facts).** (i) The product of two lower (upper) triangular matrices is lower (upper) triangular, with diagonal entries the products of the diagonal entries. (ii) A triangular matrix is invertible iff its diagonal entries are all nonzero, and the inverse of a unit lower triangular matrix is unit lower triangular.

*Proof.* (i) For lower triangular $A, B$ and $i < j$: $(AB)_{ij} = \sum_k a_{ik}b_{kj}$, and each term has $k > i$ (so $a_{ik} = 0$) or $k \le i < j$ (so $b_{kj} = 0$). On the diagonal only $k = i$ survives. Transpose for the upper case. (ii) Take $U$ upper triangular (transpose for lower). If every $u_{ii} \ne 0$, Theorem 3.6 solves $U\mathbf{x} = \mathbf{e}_j$ for each $j$, giving $X$ with $UX = I$; applied to $U^\top$ (lower triangular, same nonzero diagonal, with a division by $u_{ii}$ at each step of (3.4)) it gives $Y$ with $U^\top Y = I$, i.e. $Y^\top U = I$. A matrix with a left inverse and a right inverse is invertible and they coincide (Theorem 3.15 below is this one-line fact). If some $u_{kk} = 0$: columns $1, \dots, k$ of $U$ are nonzero only in rows $1, \dots, k-1$ (column $j \le k$ has zeros below row $j$, and column $k$ also in row $k$), so by Corollary 3.9 applied to that $(k-1)\times k$ block there is $\mathbf{y} \ne \mathbf{0}$ in $\mathbb{R}^k$ killed by it; padding $\mathbf{y}$ with zeros gives $U\mathbf{x} = \mathbf{0}$ with $\mathbf{x} \ne \mathbf{0}$, and an invertible $U$ would force $\mathbf{x} = U^{-1}U\mathbf{x} = \mathbf{0}$. For unit lower $L$, solving $LX = I$ column by column with (3.4): column $j$ has $x_1 = \cdots = x_{j-1} = 0$ and $x_j = 1$, so $X$ is unit lower triangular. ∎

**Theorem 3.12 (Uniqueness of $LU$).** If $A$ is invertible and $A = L_1U_1 = L_2U_2$ with $L_i$ unit lower triangular and $U_i$ upper triangular, then $L_1 = L_2$ and $U_1 = U_2$.

*Proof.* $U_1 = L_1^{-1}A$ is invertible, likewise $U_2$. From $L_1U_1 = L_2U_2$: $L_2^{-1}L_1 = U_2U_1^{-1}$. The left side is unit lower triangular (Lemma 3.11), the right side upper triangular. A matrix that is both is $I$. ∎

### 3.6 Row swaps: $PA = LU$ and why pivoting is not optional

Two reasons to swap rows. The first is necessity: $\begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix}$ is invertible but has no $LU$ with unit lower $L$ (the $(1,1)$ entry of $LU$ is $u_{11}$ and the $(2,1)$ entry is $\ell_{21}u_{11}$, which cannot be $1$ when $u_{11} = 0$). The second is arithmetic, and it is the one that matters.

**The $10^{-20}$ example.** Let $\varepsilon = 10^{-20}$ and solve

$$
\begin{pmatrix} \varepsilon & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2\end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}.
$$

The exact solution is $x_1 = 1/(1-\varepsilon)$, $x_2 = 1 - \varepsilon x_1$: both equal $1$ to twenty digits. Eliminate without swapping in IEEE double precision (about 16 significant digits). The multiplier is $\ell_{21} = 1/\varepsilon = 10^{20}$. The new $(2,2)$ entry is $1 - 10^{20}$, which the computer rounds to $-10^{20}$: the $1$ is twenty orders of magnitude below the rounding threshold and vanishes. The right side becomes $2 - 10^{20} \to -10^{20}$. Back substitution: $x_2 = 1$ (correct, by luck), then $x_1 = (1 - x_2)/\varepsilon = 0$. **Wrong in the first digit.** The computed $L$ and $U$ multiply to $\begin{pmatrix} \varepsilon & 1 \\ 1 & 0 \end{pmatrix}$: the entry $a_{22} = 1$ has been replaced by $0$. The algorithm solved a different matrix.

Swap the rows first: multiplier $\ell_{21} = \varepsilon$, new $(2,2)$ entry $1 - \varepsilon \to 1$, right side $1 - 2\varepsilon \to 1$. Then $x_2 = 1$, $x_1 = 2 - 1 = 1$. Correct.

The lesson is not "small pivots are bad numbers"; the matrix is perfectly nice. The lesson is that **a large multiplier amplifies rounding error**: subtracting $\ell\times(\text{row } j)$ with $|\ell| \gg 1$ swamps the information in row $i$. **Partial pivoting** is the rule: at step $k$, choose as pivot the entry of largest absolute value in column $k$ on or below the diagonal. Then $|\ell_{ik}| \le 1$ and no entry of $U$ can grow by more than a factor of $2$ per step. Every production library (LAPACK's `getrf`, behind `np.linalg.solve` and `scipy.linalg.lu`) does this.

**Theorem 3.13 ($PA = LU$).** Elimination with row swaps on an invertible $A$ produces a permutation matrix $P$, a unit lower triangular $L$, and an upper triangular $U$ with

$$
PA = LU. \tag{3.6}
$$

Under partial pivoting, $|L_{ij}| \le 1$. $P$ records the net reordering of the rows of $A$.

▶ **Try it first.** For $n = 3$ the operations are $E_2P_2E_1P_1A = U$: $P_1$ swaps row 1 with some row, $E_1$ clears column 1, $P_2$ swaps rows 2 and 3 (or does nothing), $E_2$ clears $(3,2)$. The problem is the $P_2$ between $E_1$ and $E_2$. Show $P_2E_1 = E_1'P_2$ where $E_1'$ is $E_1$ with its two multipliers swapped.

**Derivation.**

1. $E_1 = I - \ell_{21}\mathbf{e}_2\mathbf{e}_1^\top - \ell_{31}\mathbf{e}_3\mathbf{e}_1^\top$ and $P_2 = P_{23}$. Since $P_2^2 = I$, $P_2E_1 = (P_2E_1P_2)P_2$. Compute $E_1' := P_2E_1P_2 = I - \sum_{m=2,3}\ell_{m1}(P_2\mathbf{e}_m)(\mathbf{e}_1^\top P_2)$.
2. $P_2\mathbf{e}_2 = \mathbf{e}_3$, $P_2\mathbf{e}_3 = \mathbf{e}_2$ ($P_{23}$ swaps $\mathbf{e}_2$ and $\mathbf{e}_3$: its columns 2 and 3 are $\mathbf{e}_3$ and $\mathbf{e}_2$), and $\mathbf{e}_1^\top P_2 = (P_2^\top\mathbf{e}_1)^\top = (P_2\mathbf{e}_1)^\top = \mathbf{e}_1^\top$ since $P_{23}$ is symmetric and fixes index $1$. So $E_1' = I - \ell_{21}\mathbf{e}_3\mathbf{e}_1^\top - \ell_{31}\mathbf{e}_2\mathbf{e}_1^\top$: the same multipliers in swapped rows, still unit lower triangular.
3. Therefore $U = E_2E_1'(P_2P_1)A$. Set $P = P_2P_1$ (a product of permutation matrices is one, §3.10) and $L = (E_2E_1')^{-1} = E_1'^{-1}E_2^{-1}$, unit lower triangular with the reordered multipliers by Theorem 3.10's second derivation. Then $PA = LU$.
4. For general $n$ the identity $P_kE_j = E_j'P_k$ holds whenever $P_k$ swaps two rows below row $j$, which is always true for a swap at a later step $k > j$: $P_k$ fixes indices $1, \dots, k-1 \ni j$, so $\mathbf{e}_j^\top P_k = \mathbf{e}_j^\top$, and $P_k$ permutes the column-$j$ multipliers among rows $k, \dots, n$. Moving every $P_k$ to the right gives $E_{n-1}'\cdots E_1'\,P_{n-1}\cdots P_1\,A = U$. ∎

■ **What just happened.** In code this is one line: *when you swap rows $k$ and $p$ of the working matrix, also swap rows $k$ and $p$ in the already-computed columns $1, \dots, k-1$ of $L$.* That is $E_j \to E_j'$. Forgetting it is the standard bug in a hand-written $LU$; the project's asserts against `scipy.linalg.lu` catch it.

Honesty about the guarantee: $|\ell| \le 1$ allows the entries of $U$ to grow by up to $2^{n-1}$ overall, and there are matrices (Wilkinson's, exercise C1) that achieve it. So partial pivoting is not provably stable in the worst case. It is used universally because, for matrices that arise in practice, the growth is almost always small. Why that is so is a research-level question; LA-13 states what is known (Trefethen & Bau lecture 22).

### 3.7 Operation counts

▶ **Try it first.** At step $k$ of elimination on an $n\times n$ matrix, how many rows change, and how many entries in each? Count one division per multiplier and one multiply plus one subtract per updated entry. Don't count the entries you already know become zero.

**Derivation (cost of $LU$).**

1. At step $k$ ($k = 1, \dots, n-1$) there are $m = n - k$ rows below the pivot row, each needing one division for its multiplier: $m$ flops.
2. Each of those rows has $m$ entries right of column $k$ that change (column $k$ itself we *know* becomes zero), each by one multiply and one subtract: $2m^2$ flops.
3. Total: $\displaystyle\sum_{m=1}^{n-1}(m + 2m^2) = \frac{n(n-1)}{2} + \frac{(n-1)n(2n-1)}{3} = \frac23 n^3 - \frac12 n^2 - \frac16 n$, using $\sum_{m=1}^{N}m^2 = N(N+1)(2N+1)/6$. The leading term is
   $$\text{flops}(LU) \approx \tfrac{2}{3}n^3. \tag{3.7}$$

**Derivation (cost of the substitutions).** Forward substitution (3.4): row $i$ has $i-1$ multiplies and $i-1$ subtracts, total $\sum_i 2(i-1) = n(n-1) \approx n^2$. Back substitution (3.3): row $i$ has $n-i$ multiplies, $n-i$ subtracts, one divide: $\sum_i(2(n-i)+1) = n^2$ exactly. One solve, given $L$ and $U$, costs about $2n^2$ flops.

■ **What just happened.** Three consequences you should be able to recite.

- **Factor once, solve many.** $k$ right-hand sides cost $\tfrac23 n^3 + 2kn^2$ with one $LU$, versus $k\cdot\tfrac23 n^3$ from scratch. For $n = k = 1000$: about $2.7\times 10^9$ flops versus $6.7\times 10^{11}$, a factor of $250$. The project times the same comparison at $k = 10$ and $n \le 400$ (`timing_lu_vs_solve.png`), where the predicted ratio $k/(1 + 3k/n)$ is about $9.3$ at $n = 400$; the gap grows with $k$.
- **The inverse costs three factorizations.** $A^{-1}$ means solving $A\mathbf{x} = \mathbf{e}_j$ for $n$ right-hand sides. Naively $\tfrac23 n^3 + 2n\cdot n^2 = \tfrac83 n^3$; exploiting the zeros in $\mathbf{e}_j$ (exercise B8), about $2n^3$. At least three times the cost of $LU$.
- **The inverse saves nothing per solve.** $A^{-1}\mathbf{b}$ is $n^2$ multiplies and $n^2$ adds: $2n^2$ flops, the same as two substitutions. So $A^{-1}$ costs $3\times$ up front and gains nothing afterward.

There is also an accuracy argument, which should be stated carefully. The $LU$ route is **backward stable**: the computed $\hat{\mathbf{x}}$ exactly solves a nearby problem $(A + \delta A)\hat{\mathbf{x}} = \mathbf{b}$ with $\|\delta A\|$ about machine precision times $\|A\|$ times the growth factor (Wilkinson; LA-13), so the residual $\mathbf{b} - A\hat{\mathbf{x}}$ is tiny relative to $\|A\|\|\hat{\mathbf{x}}\|$ however ill-conditioned $A$ is. The route $\hat{\mathbf{x}} = \widehat{A^{-1}}\mathbf{b}$ has no such guarantee; the best general bound on its residual carries an extra factor of the condition number $\kappa(A) = \|A\|\|A^{-1}\|$ (LA-9, LA-13). In practice, with a good library in double precision, the two residuals are often within a small factor of each other; the project prints both and you should look rather than take folklore on faith. The arguments that *always* hold are cost, and structure: $A^{-1}$ is dense even when $A$ is sparse. The $n\times n$ version of Example 4.3's tridiagonal matrix has $3n-2$ nonzeros and an $LU$ with the same count; its inverse has $n^2$ nonzeros, none of them zero. For the $10^6\times 10^6$ systems of EM-3 and CM-6 that is the difference between feasible and impossible. **Never form $A^{-1}$ to solve a system.** Call `solve`, or factor and substitute.

### 3.8 Inverses

**Definition 3.14.** A square matrix $A$ is **invertible** (nonsingular) if there is $B$ with $AB = BA = I$; then $B$ is $A^{-1}$. Otherwise $A$ is **singular**.

**Theorem 3.15 (Uniqueness; one side suffices).** If $BA = I$ and $AC = I$, then $B = C$. So a matrix has at most one inverse, and any left inverse of a square matrix equals any right inverse.

*Proof.* $B = B(AC) = (BA)C = C$. ∎

**Theorem 3.16 (Invertible $\iff$ $n$ pivots).** An $n\times n$ matrix $A$ is invertible if and only if elimination produces $n$ pivots.

▶ **Try it first.** One direction is Theorem 3.8 plus one line. For the other, *build* $B$ with $BA = I$: which row operations take $U$ (with $n$ nonzero pivots) all the way to $I$?

**Derivation.**

($\Leftarrow$) Suppose $n$ pivots: $EA = U$ with $E$ a product of elementary matrices and $U$ upper triangular with nonzero diagonal.

1. Clear *above* the pivots, last column first. For column $n$: for each $i < n$ apply $E_{in}(u_{in}/u_{nn})$. Row $n$ is $(0, \dots, 0, u_{nn})$, so this changes only column $n$ of the other rows. Repeat for column $n-1$ using row $n-1$, which is zero left of column $n-1$ and now also in column $n$. Continue to column $2$. The result is $\operatorname{diag}(u_{11}, \dots, u_{nn})$.
2. Apply $D_i(1/u_{ii})$ for each $i$. The result is $I$.
3. So $MA = I$ with $M$ the product of all elementary matrices used; $M$ is invertible.
4. Multiply $MA = I$ on the left by $M^{-1}$: $A = M^{-1}$. Then $AM = M^{-1}M = I$. So $A^{-1} = M$.

($\Rightarrow$) Suppose $A$ is invertible with fewer than $n$ pivots.

5. By Theorem 3.8 there is $\mathbf{x} \ne \mathbf{0}$ with $A\mathbf{x} = \mathbf{0}$.
6. Then $\mathbf{x} = A^{-1}A\mathbf{x} = A^{-1}\mathbf{0} = \mathbf{0}$. Contradiction. ∎

■ **What just happened.** With $n$ pivots, more row operations drive $U$ on to $I$, and the product of every operation used is $A^{-1}$; with fewer, there is a nonzero null vector, which no invertible matrix can have. So invertibility is decided by counting pivots, and the proof of $\Leftarrow$ is already an algorithm for $A^{-1}$ (Corollary 3.17).

**Corollary 3.17 (Gauss–Jordan computes the inverse).** Row-reduce $[A\,|\,I]$ until the left block is $I$. The right block is then $A^{-1}$:

$$
M[A\,|\,I] = [MA\,|\,M] = [I\,|\,M] = [I\,|\,A^{-1}]. \tag{3.8}
$$

If the left block cannot reach $I$ (fewer than $n$ pivots), $A$ is singular.

*Proof.* Row operations act on both blocks (Theorem 3.1), and the operations that turn $A$ into $I$ are the $M$ of step 3, which is $A^{-1}$ by step 4. ∎

**Corollary 3.18.** $A$ is singular iff $A\mathbf{x} = \mathbf{0}$ has a nonzero solution. If $A$ is invertible, $A\mathbf{x} = \mathbf{b}$ has exactly one solution for every $\mathbf{b}$, namely $A^{-1}\mathbf{b}$.

*Proof.* Singular $\iff$ fewer than $n$ pivots (Theorem 3.16) $\Rightarrow$ nonzero null vector (Theorem 3.8); conversely a nonzero null vector forces singularity by step 6. $A(A^{-1}\mathbf{b}) = \mathbf{b}$ gives existence; if $A\mathbf{x} = \mathbf{b}$ then $\mathbf{x} = A^{-1}A\mathbf{x} = A^{-1}\mathbf{b}$. ∎

**Theorem 3.19 (Algebra of inverses).** For invertible $n\times n$ $A$, $B$ and $c \ne 0$:
(a) $(AB)^{-1} = B^{-1}A^{-1}$; (b) $(A^\top)^{-1} = (A^{-1})^\top$; (c) $(A^{-1})^{-1} = A$; (d) $(cA)^{-1} = c^{-1}A^{-1}$.

*Proof.* By Theorem 3.15 a one-sided inverse suffices. (a) $(AB)(B^{-1}A^{-1}) = A(BB^{-1})A^{-1} = I$. Note the order reversal, as in $(AB)^\top = B^\top A^\top$: to undo "do $B$ then $A$," undo $A$ first. (b) $(A^{-1})^\top A^\top = (AA^{-1})^\top = I$ by LA-1 Theorem 3.8. (c) $AA^{-1} = I$ says $A$ is a right inverse of $A^{-1}$. (d) $(cA)(c^{-1}A^{-1}) = I$. ∎

The converse of (a) also holds: if $AB$ is invertible, so are $A$ and $B$ (exercise B5). For $2\times 2$ there is a closed form, to be verified once and never used for anything larger: $\begin{pmatrix} a & b \\ c & d\end{pmatrix}^{-1} = \dfrac{1}{ad - bc}\begin{pmatrix} d & -b \\ -c & a\end{pmatrix}$. The condition $ad - bc \ne 0$ is "two pivots": if $a \ne 0$ the second pivot is $d - (c/a)b = (ad - bc)/a$. LA-6 shows the determinant of any square matrix is $\pm$ the product of its pivots.

### 3.9 Rank, provisionally

**Definition 3.20 (Rank, provisional).** The **rank** of $A$, written $r$ or $\operatorname{rank}A$, is the number of pivots produced by elimination.

What we can say now: $r \le \min(m, n)$ (at most one pivot per row and per column). A square matrix is invertible iff $r = n$. $A\mathbf{x} = \mathbf{b}$ is solvable for *every* $\mathbf{b}$ iff $r = m$ (every row has a pivot, so no row of $U$ is zero and outcome (a) cannot occur). The solution set, when nonempty, has $n - r$ free parameters.

What we cannot yet say honestly: that this number is a property of $A$ rather than of the particular row operations you chose. Different swap choices give different echelon forms. The pivot count is in fact the same for all of them, and exercise C3 gives a direct proof: a column is a pivot column exactly when it is not a combination of the columns to its left, and that property survives invertible row operations (Theorem 3.3 applied to the homogeneous system). LA-3 does better: it proves rank $= \dim C(A) = \dim C(A^\top)$ (LA-3 Theorems 3.28 and 3.29), a statement with no algorithm in it. The RREF $R$ itself is unique, and that follows from C3 here. By C3 the pivot columns $p_1 < \cdots < p_r$ are fixed by $A$. In a free column $f$ of $R$, the entry in row $i$ is the coefficient of pivot column $p_i$ when column $f$ of $R$ is written as a combination of the pivot columns (which in $R$ are $\mathbf{e}_1, \dots, \mathbf{e}_r$). By Theorem 3.3 (as in C3(b)) the same coefficients write column $f$ of $A$ in terms of the pivot columns of $A$, and they are the only such coefficients, because a relation among the pivot columns of $A$ would be one among $\mathbf{e}_1, \dots, \mathbf{e}_r$. So every entry of $R$ is determined by $A$. Until LA-3, treat "rank" as "pivot count," and know that in floating point the count depends on a tolerance: the project shows a matrix of exact rank 2 whose computed third pivot is about $10^{-16}$, not $0$.

### 3.10 $LDU$, symmetric matrices, and the special matrices

**Theorem 3.21 ($A = LDU$).** If $A$ is invertible and $A = LU$ with no row swaps, let $D = \operatorname{diag}(u_{11}, \dots, u_{nn})$. Then

$$
A = L\,D\,U', \qquad U' = D^{-1}U \text{ unit upper triangular}, \tag{3.9}
$$

and this factorization into (unit lower)·(diagonal)·(unit upper) is unique.

*Proof.* $D$ is invertible (nonzero pivots); $D^{-1}U$ divides row $i$ of $U$ by $u_{ii}$, giving ones on the diagonal; $LDU' = LU = A$. Uniqueness: if $L_1D_1U_1' = L_2D_2U_2'$, these are two $LU$ factorizations with upper factors $D_iU_i'$, so by Theorem 3.12 $L_1 = L_2$ and $D_1U_1' = D_2U_2'$. The diagonal of $D_iU_i'$ is that of $D_i$ (Lemma 3.11(i)), so $D_1 = D_2$ and then $U_1' = U_2'$. ∎

**Theorem 3.22 (Symmetric $\Rightarrow$ $A = LDL^\top$).** If $A = A^\top$ is invertible and has an $LDU'$ factorization (3.9), then $U' = L^\top$:

$$
A = LDL^\top. \tag{3.10}
$$

▶ **Try it first.** Transpose both sides of $A = LDU'$ and use uniqueness.

*Proof.* $A = A^\top = (LDU')^\top = U'^\top D L^\top$ (LA-1 Theorem 3.8 twice; $D^\top = D$). Here $U'^\top$ is unit lower triangular, $D$ diagonal, $L^\top$ unit upper triangular: an $LDU'$-type factorization of $A$. By uniqueness, $U'^\top = L$ and $L^\top = U'$. ∎

■ **What just happened.** For a symmetric matrix you need only the lower half of the elimination: $L$ and the pivots determine everything, and the work halves to about $n^3/3$ flops. When all pivots are positive, $A = (LD^{1/2})(LD^{1/2})^\top$ is the **Cholesky factorization**, and "all pivots positive" is one of the equivalent definitions of *positive definite* in LA-8. The project's "Your turn" has you implement (3.10).

**Special matrices, collected.**

- **Diagonal.** Invertible iff every diagonal entry is nonzero; the inverse is the diagonal of reciprocals. On the left $D$ scales rows; on the right, columns (Theorem 3.1).
- **Triangular.** Lemma 3.11: products and inverses stay triangular; invertible iff no zero on the diagonal; the pivots are the diagonal.
- **Permutation.** $P$ has exactly one $1$ in each row and column; its rows are $\mathbf{e}_{\pi(1)}^\top, \dots, \mathbf{e}_{\pi(n)}^\top$ for a bijection $\pi$ of $\{1, \dots, n\}$, so by Theorem 3.1 row $i$ of $PA$ is row $\pi(i)$ of $A$. Read by columns, $P\mathbf{e}_j = \mathbf{e}_{\pi^{-1}(j)}$: column $j$ has its $1$ in the row $i$ with $\pi(i) = j$ (in A6, rows $\mathbf{e}_3^\top, \mathbf{e}_1^\top, \mathbf{e}_2^\top$ give $P\mathbf{e}_1 = \mathbf{e}_2$). There are $n!$ of them, and products of permutation matrices are permutation matrices: row $i$ of $P_1P_2$ is $\mathbf{e}_{\pi_1(i)}^\top P_2 = (\text{row } \pi_1(i) \text{ of } P_2) = \mathbf{e}_{\pi_2(\pi_1(i))}^\top$, so $P_1P_2$ is the permutation matrix of $\pi_2\circ\pi_1$.

**Theorem 3.23.** For every permutation matrix $P$, $PP^\top = P^\top P = I$; that is, $P^{-1} = P^\top$.

*Proof.* $(PP^\top)_{ij} = (\text{row } i)\cdot(\text{row } j) = \mathbf{e}_{\pi(i)}\cdot\mathbf{e}_{\pi(j)} = \delta_{ij}$, because $\pi$ is a bijection. So $PP^\top = I$, hence $P^\top = P^{-1}$ by Theorem 3.15, and $P^\top P = I$ too. ∎

- **Symmetric.** $A = A^\top$. For any $A$, $A^\top A$ and $AA^\top$ are symmetric: $(A^\top A)^\top = A^\top A$. A symmetric invertible $A$ has symmetric inverse: $(A^{-1})^\top = (A^\top)^{-1} = A^{-1}$ (Theorem 3.19(b)). Symmetric matrices are where physics lives: stiffness matrices (CM-6), conductance matrices (Example 4.4), inertia tensors (CM-7), Hessians (ML-2), covariances (F5), and (3.10) is the first of several factorizations (Cholesky in LA-10, spectral in LA-8) that exploit the symmetry.

## 4. Worked examples

### Example 4.1: A full $3\times 3$ elimination with $L$ and $U$

Solve $A\mathbf{x} = \mathbf{b}$ with

$$
A = \begin{pmatrix} 2 & 1 & 1 \\ 4 & 1 & 0 \\ -2 & 2 & 1 \end{pmatrix}, \qquad \mathbf{b} = \begin{pmatrix} 3 \\ 6 \\ 1 \end{pmatrix}.
$$

*Step 1, column 1.* Pivot $u_{11} = 2$. Multipliers $\ell_{21} = 4/2 = 2$, $\ell_{31} = -2/2 = -1$.
Row 2 ← row 2 − $2$·row 1 $= (4, 1, 0) - (4, 2, 2) = (0, -1, -2)$.
Row 3 ← row 3 − $(-1)$·row 1 $= (-2, 2, 1) + (2, 1, 1) = (0, 3, 2)$. (Subtracting a negative multiple is adding; the stored multiplier is $-1$.)

*Step 2, column 2.* Pivot $u_{22} = -1$. Multiplier $\ell_{32} = 3/(-1) = -3$.
Row 3 ← row 3 − $(-3)$·row 2 $= (0, 3, 2) + 3(0, -1, -2) = (0, 0, -4)$.

$$
U = \begin{pmatrix} 2 & 1 & 1 \\ 0 & -1 & -2 \\ 0 & 0 & -4 \end{pmatrix}, \qquad
L = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ -1 & -3 & 1 \end{pmatrix}.
$$

Pivots $2, -1, -4$: three of them, so $A$ is invertible (Theorem 3.16).

*Check $A = LU$* by Theorem 3.10's recipe, row 3: $\ell_{31}\mathbf{u}_1^\top + \ell_{32}\mathbf{u}_2^\top + \mathbf{u}_3^\top = -(2,1,1) - 3(0,-1,-2) + (0,0,-4) = (-2, 2, 1)$. ✓ Row 2: $2(2,1,1) + (0,-1,-2) = (4, 1, 0)$. ✓

*Forward substitution (3.4):* $c_1 = 3$; $c_2 = 6 - 2\cdot 3 = 0$; $c_3 = 1 - (-1)(3) - (-3)(0) = 4$. $\mathbf{c} = (3, 0, 4)^\top$, which is exactly what the two elimination steps do to $\mathbf{b}$.

*Back substitution (3.3):* $-4x_3 = 4 \Rightarrow x_3 = -1$. $-x_2 - 2(-1) = 0 \Rightarrow x_2 = 2$. $2x_1 + 2 - 1 = 3 \Rightarrow x_1 = 1$. $\mathbf{x} = (1, 2, -1)^\top$. *Check:* $A\mathbf{x} = (2 + 2 - 1,\ 4 + 2,\ -2 + 4 - 1) = (3, 6, 1)$. ✓

*The forward product is not "minus $L$."* With $E_{21} = E_{21}(2)$, $E_{31} = E_{31}(-1)$, $E_{32} = E_{32}(-3)$:

$$
E_{32}E_{31}E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ -5 & 3 & 1 \end{pmatrix}.
$$

The $(3,1)$ entry is $-\ell_{31} + \ell_{32}\ell_{21} = 1 + (-3)(2) = -5$, the surviving cross term of §3.5. Its inverse is $L$, whose $(3,1)$ entry is $-1$.

### Example 4.2: A singular $3\times 3$, and what the echelon form says

$$
A = \begin{pmatrix} 1 & 2 & 3 \\ 2 & 4 & 7 \\ 3 & 6 & 10 \end{pmatrix}.
$$

*Column 1.* Pivot $1$. $\ell_{21} = 2$: row 2 ← $(2,4,7) - 2(1,2,3) = (0, 0, 1)$. $\ell_{31} = 3$: row 3 ← $(3,6,10) - 3(1,2,3) = (0, 0, 1)$.

*Column 2.* Positions $(2,2)$ and $(3,2)$ are both $0$. No swap helps. **Column 2 is free.** This is the informative kind of breaking: column 2 is $2\times$ column 1.

*Column 3.* Pivot row is still row 2; pivot $1$ at $(2,3)$. $\ell_{32} = 1$: row 3 ← $(0,0,1) - (0,0,1) = \mathbf{0}$.

$$
U = \begin{pmatrix} 1 & 2 & 3 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}.
$$

Two pivots (columns 1, 3), one free column (2), one zero row. Rank $2$; singular by Theorem 3.16. By Theorem 3.8, $A\mathbf{x} = \mathbf{0}$ has a nonzero solution: $x_2 = 1$, then $x_3 = 0$ from row 2 and $x_1 = -2$ from row 1. $\mathbf{s} = (-2, 1, 0)^\top$, and $A\mathbf{s} = -2(\text{col } 1) + (\text{col } 2) = \mathbf{0}$.

*Right-hand side 1:* $\mathbf{b} = (1, 3, 4)^\top$. Same operations: $c_2 = 3 - 2 = 1$, $c_3 = 4 - 3 = 1$, then $c_3 \leftarrow 1 - 1 = 0$. Rows of the augmented echelon form: $(1, 2, 3\,|\,1)$, $(0, 0, 1\,|\,1)$, $(0,0,0\,|\,0)$. Consistent. Theorem 3.7(c): $x_3 = 1$; $x_1 = 1 - 2x_2 - 3$. With $x_2 = t$:

$$
\mathbf{x} = \begin{pmatrix} -2 \\ 0 \\ 1 \end{pmatrix} + t\begin{pmatrix} -2 \\ 1 \\ 0 \end{pmatrix}.
$$

Check: $A(-2, 0, 1)^\top = (-2 + 3, -4 + 7, -6 + 10) = (1, 3, 4)$. ✓

*Right-hand side 2:* $\mathbf{b} = (1, 3, 5)^\top$. Now $c_3 = 5 - 3 = 2$, then $2 - 1 = 1$. Last row $(0, 0, 0\,|\,1)$: $0 = 1$. **No solution.** The consistency condition, read off the last row, is $b_3 - 3b_1 - (b_2 - 2b_1) = 0$, i.e. $b_3 = b_1 + b_2$: $4 = 1 + 3$ ✓, $5 \ne 1 + 3$ ✗. In LA-1's language, the columns span only a plane and $(1,3,5)^\top$ is off it.

### Example 4.3: Gauss–Jordan inverse of a $3\times 3$, and its $LDL^\top$

$$
A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 2 & 1 \\ 0 & 1 & 2 \end{pmatrix}.
$$

(Up to signs, this is the second-difference matrix. With $S = \operatorname{diag}(1, -1, 1)$, flipping the sign of $x_2$ and of equation 2 turns $A$ into $SAS = \operatorname{tridiag}(-1, 2, -1)$, which discretizes $-d^2/dx^2$ in EM-3 and is the stiffness matrix of three beads on a string in CM-6. The flip changes no pivot, only the signs of some off-diagonal entries of $L$, $U$ and $A^{-1}$; the project checks this.) Form $[A\,|\,I]$.

*Down.* $\ell_{21} = \tfrac12$: row 2 ← row 2 − $\tfrac12$ row 1 $= (0, \tfrac32, 1\,|\,-\tfrac12, 1, 0)$. $\ell_{31} = 0$. $\ell_{32} = 1/\tfrac32 = \tfrac23$: row 3 ← row 3 − $\tfrac23$ row 2 $= (0, 0, \tfrac43\,|\,\tfrac13, -\tfrac23, 1)$.

$$
\left[\begin{array}{ccc|ccc} 2 & 1 & 0 & 1 & 0 & 0 \\ 0 & \tfrac32 & 1 & -\tfrac12 & 1 & 0 \\ 0 & 0 & \tfrac43 & \tfrac13 & -\tfrac23 & 1 \end{array}\right]
$$

Pivots $2, \tfrac32, \tfrac43$, product $4$ (the determinant, LA-6). Three pivots: invertible.

*Up.* Row 2 ← row 2 − $\tfrac{1}{4/3}$ row 3 $=$ row 2 − $\tfrac34$ row 3 $= (0, \tfrac32, 0\,|\,-\tfrac34, \tfrac32, -\tfrac34)$. Row 1 already has $0$ in column 3. Row 1 ← row 1 − $\tfrac{1}{3/2}$ row 2 $=$ row 1 − $\tfrac23$ row 2 $= (2, 0, 0\,|\,\tfrac32, -1, \tfrac12)$.

*Scale.* Row 1 $\div 2$, row 2 $\div\tfrac32$, row 3 $\div\tfrac43$:

$$
\left[\begin{array}{ccc|ccc} 1 & 0 & 0 & \tfrac34 & -\tfrac12 & \tfrac14 \\ 0 & 1 & 0 & -\tfrac12 & 1 & -\tfrac12 \\ 0 & 0 & 1 & \tfrac14 & -\tfrac12 & \tfrac34 \end{array}\right],
\qquad
A^{-1} = \frac14\begin{pmatrix} 3 & -2 & 1 \\ -2 & 4 & -2 \\ 1 & -2 & 3 \end{pmatrix}.
$$

*Check:* row 2 of $A$ times column 1 of $A^{-1}$: $\tfrac14(3 - 4 + 1) = 0$ ✓; times column 2: $\tfrac14(-2 + 8 - 2) = 1$ ✓. $A^{-1}$ is symmetric, as §3.10 predicts, and **dense**: $A$ has four zeros, $A^{-1}$ none. That is §3.7's structural argument in miniature.

*$LDL^\top$.* From the downward pass, $L = \begin{pmatrix} 1 & 0 & 0 \\ \tfrac12 & 1 & 0 \\ 0 & \tfrac23 & 1 \end{pmatrix}$ and $D = \operatorname{diag}(2, \tfrac32, \tfrac43)$. Theorem 3.22 predicts $U' = L^\top$, so $U = DL^\top$ should have $u_{23} = \tfrac32\cdot\tfrac23 = 1$; the downward pass gave $u_{23} = 1$. ✓ All pivots positive: positive definite (LA-8).

### Example 4.4: A Kirchhoff circuit, and what the pivots mean physically

Three nodes with unknown voltages $v_1, v_2, v_3$ relative to ground. Resistors: $R_1 = 1\,\Omega$ from node 1 to ground, $R_2 = 2\,\Omega$ between nodes 1 and 2, $R_3 = 2\,\Omega$ between nodes 2 and 3, $R_4 = 1\,\Omega$ from node 3 to ground. A current source injects $3\,\mathrm{A}$ into node 1.

```
        3 A
         |
   v1 ---+---[R2 = 2]--- v2 ---[R3 = 2]--- v3
   |                                       |
 [R1 = 1]                               [R4 = 1]
   |                                       |
  gnd                                     gnd
```

Kirchhoff's current law at each node, with current through a resistor $= g\cdot(\text{voltage drop})$, $g = 1/R$: conductances $g_1 = 1, g_2 = \tfrac12, g_3 = \tfrac12, g_4 = 1$. Node 1: $g_1v_1 + g_2(v_1 - v_2) = 3$. Node 2: $g_2(v_2 - v_1) + g_3(v_2 - v_3) = 0$. Node 3: $g_3(v_3 - v_2) + g_4v_3 = 0$. In matrix form $G\mathbf{v} = \mathbf{i}$:

$$
G = \begin{pmatrix} \tfrac32 & -\tfrac12 & 0 \\ -\tfrac12 & 1 & -\tfrac12 \\ 0 & -\tfrac12 & \tfrac32 \end{pmatrix}, \qquad \mathbf{i} = \begin{pmatrix} 3 \\ 0 \\ 0 \end{pmatrix}.
$$

$G$ is symmetric (the conductance between two nodes is the same from either end); $G_{kk}$ is the total conductance touching node $k$; $G_{jk}$ is minus the conductance between $j$ and $k$. Every conductance matrix has this form.

*Elimination.* Pivot $u_{11} = \tfrac32$. $\ell_{21} = -\tfrac12/\tfrac32 = -\tfrac13$, $\ell_{31} = 0$. Row 2 ← row 2 + $\tfrac13$ row 1 $= (0, 1 - \tfrac16, -\tfrac12) = (0, \tfrac56, -\tfrac12)$. Pivot $u_{22} = \tfrac56$. $\ell_{32} = -\tfrac12/\tfrac56 = -\tfrac35$. Row 3 ← row 3 + $\tfrac35$ row 2 $= (0, 0, \tfrac32 - \tfrac{3}{10}) = (0, 0, \tfrac65)$.

$$
L = \begin{pmatrix} 1 & 0 & 0 \\ -\tfrac13 & 1 & 0 \\ 0 & -\tfrac35 & 1 \end{pmatrix}, \qquad
U = \begin{pmatrix} \tfrac32 & -\tfrac12 & 0 \\ 0 & \tfrac56 & -\tfrac12 \\ 0 & 0 & \tfrac65 \end{pmatrix}.
$$

*Solve.* Forward: $c_1 = 3$, $c_2 = 0 + \tfrac13\cdot 3 = 1$, $c_3 = 0 + \tfrac35\cdot 1 = \tfrac35$. Back: $v_3 = \tfrac35/\tfrac65 = \tfrac12$; $\tfrac56v_2 = 1 + \tfrac12\cdot\tfrac12 = \tfrac54 \Rightarrow v_2 = \tfrac32$; $\tfrac32v_1 = 3 + \tfrac12\cdot\tfrac32 = \tfrac{15}{4} \Rightarrow v_1 = \tfrac52$.

$\mathbf{v} = (2.5, 1.5, 0.5)^\top$ volts. *Physical check:* $2.5\,\mathrm{A}$ flows through $R_1$ to ground; the remaining $0.5\,\mathrm{A}$ crosses $R_2$ (drop $1\,\mathrm{V}$, $v_2 = 1.5$ ✓), crosses $R_3$ (drop $1\,\mathrm{V}$, $v_3 = 0.5$ ✓), and exits through $R_4$ ($0.5\,\mathrm{V}/1\,\Omega = 0.5\,\mathrm{A}$ ✓). Total $3\,\mathrm{A}$. ✓

*What the pivots mean.* Pivot 1, $\tfrac32 = g_1 + g_2$, is the conductance seen at node 1 with nodes 2, 3 held at zero. Pivot 2, $\tfrac56 = \tfrac13 + \tfrac12$: the $\tfrac13$ is $R_1$ and $R_2$ *in series* ($3\,\Omega$), the $\tfrac12$ is $g_3$. So pivot 2 is the conductance seen at node 2 *after node 1 has been eliminated from the circuit* and its resistors merged into an equivalent. Pivot 3, $\tfrac65 = \tfrac15 + 1$: $R_1 + R_2 + R_3 = 5\,\Omega$ in series, plus $g_4$. Each elimination step removes one node and replaces its connections by their series–parallel equivalent. Electrical engineers call this **Kron reduction**; in structural mechanics the same operation on a stiffness matrix is **static condensation**; in linear algebra it is the **Schur complement** (exercise C2). And a **zero pivot means a node with no path to ground**: remove $R_1$ and $R_4$ and the conductance matrix $\begin{pmatrix} \tfrac12 & -\tfrac12 & 0 \\ -\tfrac12 & 1 & -\tfrac12 \\ 0 & -\tfrac12 & \tfrac12\end{pmatrix}$ has pivots $\tfrac12, \tfrac12, 0$ and null vector $(1, 1, 1)^\top$: with no reference to ground, all voltages can shift by a common constant. Singular matrix, physically indeterminate potential.

### Example 4.5: The $10^{-20}$ example, as the project reproduces it

Without pivoting: $\ell_{21} = 10^{20}$, $U = \begin{pmatrix} 10^{-20} & 1 \\ 0 & \mathrm{fl}(1 - 10^{20})\end{pmatrix} = \begin{pmatrix} 10^{-20} & 1 \\ 0 & -10^{20}\end{pmatrix}$, and $LU = \begin{pmatrix} 10^{-20} & 1 \\ 1 & 10^{20} - 10^{20} \end{pmatrix} = \begin{pmatrix} 10^{-20} & 1 \\ 1 & 0\end{pmatrix} \ne A$. The computed solution of $A\mathbf{x} = (1, 2)^\top$ is $(0, 1)^\top$; the true one is $(1, 1)^\top$ to 20 digits. With pivoting the computed solution is $(1, 1)^\top$. The project sweeps $\varepsilon$ from $10^{-1}$ to $10^{-20}$ and plots both errors (`pivoting_error.png`). The unpivoted error does not jump; it grows in proportion to $1/\varepsilon$, roughly one decade for every decade $\varepsilon$ shrinks (the computed errors lie between about $10^{-18}/\varepsilon$ and $10^{-16}/\varepsilon$, machine precision divided by $\varepsilon$): $4\times 10^{-16}$ at $\varepsilon = 10^{-1}$, $5\times 10^{-9}$ at $10^{-8}$, $2\times 10^{-5}$ at $10^{-12}$, $1.2$ at $10^{-16}$, and exactly $1$ for $\varepsilon \le 10^{-17}$, where $x_1$ is computed as $0$. On the semilog plot that is a straight line, not a step. The multiplier $1/\varepsilon$ amplifies the rounding error at every $\varepsilon$; machine precision only decides where the damage reaches $O(1)$. The pivoted error stays at rounding level ($\le 2.2\times 10^{-16}$) throughout.

## 5. Why it matters and where it connects

- **LA-3** takes the special solutions $\mathbf{s}_f$ of Theorem 3.7(c) and calls their span the null space $N(A)$; the pivot columns of $A$ become a basis of the column space $C(A)$; the nonzero rows of the RREF a basis of the row space. The provisional rank becomes the common dimension of $C(A)$ and the row space, and $n - r$ becomes $\dim N(A)$ (LA-3 Theorems 3.18, 3.28, 3.29). RREF uniqueness is proved here, in §3.9, from exercise C3.
- **LA-4**: Rank–Nullity, $r + \dim N(A) = n$, is Theorem 3.7(c) without the algorithm.
- **LA-5**: least squares solves $A^\top A\hat{\mathbf{x}} = A^\top\mathbf{b}$; $A^\top A$ is symmetric (§3.10) with positive pivots when the columns of $A$ are independent, so Cholesky ($LDL^\top$ with $D > 0$) solves it. $QR$ replaces $LU$ when the matrix is rectangular.
- **LA-6**: $\det A = \pm(\text{product of pivots})$, the sign from the number of swaps. Example 4.3's pivots multiply to $4$; the project checks `np.linalg.det`.
- **LA-8**: symmetric $A$ is positive definite iff all pivots are positive iff $A = LDL^\top$ with $D > 0$. The energy $\mathbf{x}^\top A\mathbf{x} = \sum d_iy_i^2$, $\mathbf{y} = L^\top\mathbf{x}$ (exercise B7), is the first step of that proof.
- **LA-10** (Cholesky; which factorization when) and **LA-13** (backward stability of $LU$, growth factor, conditioning; Example 4.5 done properly) both start from (3.6).
- **CM-6**: small oscillations give $M\ddot{\mathbf{q}} + K\mathbf{q} = \mathbf{0}$ with $K$ symmetric and patterned like Example 4.3 and Example 4.4's $G$. Static condensation of interior degrees of freedom is §3.5's elimination; the pivots are effective stiffnesses as in Example 4.4 they are effective conductances.
- **EM-3**: the finite-difference Laplacian is banded and symmetric; $LU$ of a banded matrix stays banded ($O(n)$ for tridiagonal), while its inverse is dense (Example 4.3).
- **ST-3**: the normal equations $(X^\top X)\boldsymbol\beta = X^\top\mathbf{y}$ are solved by Cholesky or $QR$; `inv(X.T @ X)` in a regression code is §3.7's mistake.
- **ML-2**: Newton's method solves $H\mathbf{d} = -\nabla f$ at every iteration; one $LDL^\top$ per step, never $H^{-1}$.
- **QC-2**: the reversible classical gates (SWAP, CNOT, Toffoli) are permutation matrices on the computational basis, and Theorem 3.23 is why they are unitary with no further argument.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: a $2\times 2$ elimination with $E$, $L$, $U$; elementary matrices and inverses; an RREF by hand; substitution with a flop count; a $2\times 2$ Gauss–Jordan; a permutation matrix. B1–B9 are the module: B2 (full $LU$ with the forward-product contrast), B3 (prove Corollary 3.9; $PA = LU$ by hand), B5 (algebra of inverses including the converse of Theorem 3.19(a)), B7 ($LDL^\top$ and $\mathbf{x}^\top A\mathbf{x}$), and B9 (implement $LDL^\top$ in `project.py`) are essential. C1 (Wilkinson's growth matrix), C2 (block $LU$ and the Schur complement, with the circuit), C3 (the pivot count is independent of the row operations), C4 (Sherman–Morrison) go one step beyond.

## 7. Coding project

`project.py` builds the algorithms of this module from scratch and checks every number in §4.

1. **Elementary matrices.** `E(n, i, j, l)`, `Pswap`, `D` from (3.1). Asserts: left multiplication acts on rows and right on columns (Theorem 3.1); the inverses of §3.1; for Example 4.1, $E_{32}E_{31}E_{21}$ has $(3,1)$ entry $-5$ while $L$ has $-1$.
2. **$LU$ with partial pivoting.** `lu_decompose(A, pivot=True, tol=0.0)` returning `P, L, U` with $PA = LU$, including the swap of the earlier columns of $L$ (Theorem 3.13); `forward_sub` (3.4), `back_sub` (3.3), `solve(A, b)`. Asserted against `scipy.linalg.lu` (translating its $A = PLU$ convention) and `np.linalg.solve` on random matrices. Examples 4.1 and 4.4 reproduced to the digit, including pivots-as-conductances.
3. **RREF with a tolerance, and classification.** `rref(A, tol)` returns the reduced form and pivot columns; `classify(A, b)` returns `'unique'`, `'none'`, or `'infinite'` with the free columns, $\mathbf{x}_p$, and the special solutions (Theorem 3.7). Example 4.2 with both right-hand sides; LA-1's three $2\times 2$ outcomes; random systems of known rank.
4. **Gauss–Jordan inverse.** `gauss_jordan_inverse(A)` from (3.8); Example 4.3 with its $LDL^\top$; Theorem 3.19 on random matrices; the case against `inv`: density of the inverse of a tridiagonal matrix, wall-clock of `inv` versus `solve` (the flop ratio is $3$; wall-clock ratios depend on the library's threading), and the residuals of both routes, printed for you to look at.
5. **The $10^{-20}$ example.** Elimination with and without pivoting for $\varepsilon = 10^{-1}, \dots, 10^{-20}$; figure `pivoting_error.png`, with the reference line $10^{-16}/\varepsilon$. Asserts: the unpivoted error tracks $10^{-16}/\varepsilon$ (slope about one decade per decade) until it saturates at exactly $1$; the pivoted error stays at rounding level.
6. **Costs.** A flop counter inside `lu_decompose`, asserted equal to the exact sum of §3.7 and plotted against $\tfrac23n^3$ (`flop_count.png`); wall-clock of "$LU$ once then $k$ substitutions" versus "$k$ solves from scratch" for growing $n$ (`timing_lu_vs_solve.png`).
7. **Exercise numbers**, so that `solutions.md` and the code agree.

**Build it yourself first.** Before reading Part 2, write `lu_decompose` without pivoting (a loop over columns, `np.outer` for the rank-one update), then `forward_sub` and `back_sub` from (3.3) and (3.4). Make them pass Example 4.1. Then add pivoting and discover which rows of $L$ need swapping. Expected result: `np.allclose(P @ A, L @ U)` on random $50\times 50$ matrices, and your `solve` agreeing with NumPy to about $10^{-13}$.

**Your turn** (end of the script): implement `ldlt(A)` (Theorem 3.22) and check $LDL^\top = A$ and that $D$ holds the pivots; estimate rank by counting pivots above a tolerance and watch a matrix of exact rank 2 report rank 6 when the tolerance is $0$; build a random singular matrix as $BC$ with $B$ being $n\times 2$ and see what `solve` does with it.

## 8. Where you'll struggle

- **Signs of the multipliers.** The operation is row $i$ ← row $i$ − $\ell_{ij}$ · row $j$ with $\ell_{ij}$ = (entry being cleared in row $i$)/(pivot of row $j$), both taken from the *current* matrix; that is $a_{ij}/a_{jj}$ when the pivots are on the diagonal. $L$ stores $+\ell_{ij}$; $E_{ij}$ stores $-\ell_{ij}$. When $\ell_{ij}$ is negative you *add* a multiple of row $j$ and $L$ still stores the negative number (Example 4.1: $\ell_{31} = -1$, you added row 1, $L_{31} = -1$). The three ways to get it wrong: storing the negative; computing $\ell$ from the original $A$ instead of the current matrix at step 2; reading $L$ off $E_{32}E_{31}E_{21}$, which is wrong by a cross term. The fix: after every hand elimination, check one row of $A = LU$ with Theorem 3.10's recipe. Thirty seconds, catches all three.
- **Row operations versus column operations.** Left multiplication acts on rows, right multiplication on columns (Theorem 3.1). Elimination is all left multiplication, which is why $[A\,|\,I] \to [I\,|\,A^{-1}]$ works: one $M$ hits both blocks from the left. Column operations on an augmented matrix change the unknowns, not the equations, and the answer is to a different problem. Say "rows" out loud.
- **"Zero pivot" panic.** A zero in the pivot position means *look below it*. If something below is nonzero, swap and continue. Only if the entire rest of the column is zero is the column free. People either give up at the first zero (declaring singular a matrix that isn't) or try to swap when there is nothing to swap in (Example 4.2).
- **The inverse reflex.** `x = np.linalg.inv(A) @ b` is in a lot of code you'll read. It costs three times `solve`, saves nothing per right-hand side, destroys sparsity, and loses the backward-stability guarantee (§3.7). The reflex to fight is algebraic: $\mathbf{x} = A^{-1}\mathbf{b}$ is a correct *statement* and a bad *instruction*. Write $A^{-1}$ in proofs; call `solve` in code.
- **"Rank" is provisional here.** Definition 3.20 counts an algorithm's output, and we have not proved that different legal choices give the same count (C3 does; LA-3 does it properly). In floating point the count depends on a tolerance: a matrix of exact rank 2 produces a third pivot near $10^{-16}$, and a `rank` without a `tol` says 3. Numerical rank genuinely needs the SVD (LA-9). Don't build arguments on "the rank" until LA-3 defines it without an algorithm.
- **The pivoting bug.** When you add row swaps to `lu_decompose`, the swap must also be applied to the filled columns of $L$ (Theorem 3.13, step 4). Everyone forgets it once; the symptom is $PA \ne LU$ while every step looked right.
- **The direction of "free."** Free variables are chosen *first*; pivot variables are solved *from the bottom up*. Draw the staircase, mark the pivot columns, read the pivot rows from the bottom.
- **Not every invertible matrix has an $LU$.** $\begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix}$ does not (§3.6). Every invertible matrix has a $PA = LU$; the permutation is part of the factorization.
- **Exercise B3(b)** is harder than it looks because the second swap reorders column 1 of $L$. Write the permutation vector after each swap.

## 9. Check yourself

1. Write the $3\times 3$ elementary matrix that subtracts $5$ times row 1 from row 3, and its inverse. Which one has $+5$ in it?
2. Why does left-multiplying by $E$ change rows and not columns? One sentence using the row-times-column rule.
3. State the three outcomes for $A\mathbf{x} = \mathbf{b}$ in terms of pivot columns and the augmented column.
4. Prove in four lines: an $n\times n$ matrix with fewer than $n$ pivots has a nonzero null vector.
5. Why is $L$ exactly the matrix of multipliers, with no further arithmetic? Give the "run elimination backwards" argument for row 3.
6. Why is $E_{32}E_{31}E_{21}$ *not* $I$ with $-\ell_{ij}$ below the diagonal? Which entry is off, and by how much?
7. Give the $2\times 2$ matrix on which elimination without pivoting fails in floating point, and say in one sentence what goes wrong.
8. About how many flops for $LU$ of a $2000\times 2000$ matrix? For one back substitution? For $A^{-1}$?
9. State and prove $(AB)^{-1} = B^{-1}A^{-1}$.
10. Why does Gauss–Jordan on $[A\,|\,I]$ leave $A^{-1}$ in the right block?
11. For symmetric $A$ with an $LDU$ factorization, why is $U = L^\top$? Which theorem does the proof lean on?
12. In the circuit of Example 4.4, what does a pivot measure? What does a zero pivot mean physically?

## 10. Pair with

- **Strang, *Introduction to Linear Algebra*, ch. 2** (§2.2–2.7: elimination, elimination using matrices, rules for matrix operations, inverses, $A = LU$, transposes and permutations). Strang's treatment of "$L$ is the record of the multipliers" is the clearest in print; his §2.6 is the source of the $E_{32}E_{31}E_{21}$-versus-$L$ comparison.
- **MIT 18.06, lectures 2–5 (Strang), on OCW.** Lecture 2 is elimination with matrices; lecture 3 is multiplication and inverses with Gauss–Jordan done live; lecture 4 is $A = LU$ and the operation count; lecture 5 is transposes, permutations, and the start of vector spaces. Watch lecture 4 after §3.5 and notice how he counts the $n^3/3$ multiplications.
- **Trefethen & Bau, *Numerical Linear Algebra*, lectures 20–21.** Lecture 20 states Gaussian elimination as $L^{-1}A = U$ with the flop count in their style; lecture 21 is pivoting, with the $PA = LU$ derivation of Theorem 3.13 and the worst-case growth factor. Lecture 22 (stability) is for LA-13, but its first page is readable now.
- **Grcar, "How ordinary elimination became Gaussian elimination," *Historia Mathematica* 38 (2011).** Optional. The real history of §1, with primary sources.
