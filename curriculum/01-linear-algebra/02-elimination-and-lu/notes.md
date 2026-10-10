# LA-2: Solving Systems: Elimination, Echelon Forms, LU, Inverses, Rank

**Time:** ~8 hours · **Prerequisites:** LA-1 · **Unlocks:** LA-3 (and through it the whole LA core; LA-6 uses the pivots for determinants, LA-10 and LA-13 return to LU numerically)

> **Conventions in this module.** $A$ is $m\times n$ with entries $a_{ij}$ (row $i$, column $j$), rows $\mathbf{r}_1^\top,\dots,\mathbf{r}_m^\top$ and columns $\mathbf{a}_1,\dots,\mathbf{a}_n$, as in LA-1. Unless stated otherwise $A$ is square, $n\times n$. $\mathbf{e}_i$ is the $i$-th standard basis vector and $\mathbf{e}_i\mathbf{e}_j^\top$ is the matrix with a single $1$ in position $(i,j)$. The basic row operation is **"row $i$ ← row $i$ − $\ell$ · row $j$"**; the number $\ell$ is the **multiplier**, and when it is used to clear position $(i,j)$ we call it $\ell_{ij}$. Sign convention: $L$ stores $+\ell_{ij}$, the elementary matrix $E_{ij}$ that performs the operation contains $-\ell_{ij}$. **Pivots** are the entries $u_{11},\dots$ on the diagonal of the triangular matrix $U$ that elimination produces (more precisely, the first nonzero entry of each nonzero row of an echelon form). Permutations are written as $PA = LU$ (SciPy's `lu` returns the transposed convention $A = PLU$; the project checks the translation). A **flop** is one floating-point add, subtract, multiply, or divide. A pivot is "zero" numerically if $|u_{kk}| < \texttt{tol}$ for a stated tolerance.

## 0. Where this sits

LA-1 told you *when* $A\mathbf{x}=\mathbf{b}$ has a solution (when $\mathbf{b}$ is in the span of the columns) but gave you no way to find one beyond staring at $2\times 2$ examples. This module is the algorithm: Gaussian elimination, and the matrix factorization $A = LU$ that elimination secretly computes. You need LA-1's column picture, the row form of $A\mathbf{x}$ (LA-1 eq. 3.2), the column rule for $AB$ (LA-1 eq. 3.3), and $(AB)^\top = B^\top A^\top$ (LA-1 Theorem 3.8).

What it makes possible: LA-3 defines the null space as the set of solutions of $A\mathbf{x}=\mathbf{0}$ and the rank as a dimension; both are read directly off the echelon forms you compute here. LA-6 will show the determinant is $\pm$ the product of the pivots. LA-10 and LA-13 take $LU$ and ask how it behaves in floating point; this module gives you the first, decisive example of why that question matters. Every project in every later course that calls `np.linalg.solve` is calling this module.

## 1. The problem that forced this idea

Elimination is the oldest algorithm in this curriculum. Chapter 8 of the Chinese *Nine Chapters on the Mathematical Art* (compiled roughly two thousand years ago) is titled *Fangcheng*, "rectangular arrays." It lays coefficients out in a grid on a counting board and subtracts multiples of one column from another until the array is triangular, then substitutes back. That is Gaussian elimination, minus the name, minus the matrices, done with rods.

The name came much later, and the reason it stuck is worth knowing. In 1809 and 1810 Gauss was fitting the orbit of the asteroid Pallas to many observations by least squares. The "normal equations" of least squares (you will derive them in LA-5) are a square symmetric system with as many unknowns as orbital elements, six or more. Gauss needed a *procedure* that a human computer could follow without thinking, and he wrote one down: eliminate the unknowns one at a time, in a fixed order, and keep the bookkeeping tight. Nineteenth-century geodesists and astronomers adopted it as "Gauss's method," and the surveyor Wilhelm Jordan published, in the 1888 edition of his *Handbuch der Vermessungskunde*, the variant that keeps eliminating above the pivots too, so the system ends as the identity. That variant is what we now call Gauss–Jordan. (Joseph Grcar's 2011 article "How ordinary elimination became Gaussian elimination" in *Historia Mathematica* traces the naming; it is a good read once you've finished this module.)

Then two things happened in the twentieth century that turned a hand recipe into a piece of mathematics.

First, people noticed that the recipe is a **factorization**. If you record what elimination does to $A$, the record is a lower-triangular matrix $L$, the result is an upper-triangular matrix $U$, and $A = LU$. Alan Turing's 1948 paper "Rounding-off errors in matrix processes" is where the name $LU$ and the clean matrix statement appear (Banachiewicz in 1938 and Dwyer in 1944 had matrix forms of the same idea). Once elimination is a factorization you can reason about it: when it exists, when it is unique, how much it costs, and what it costs to reuse.

Second, the first electronic computers made people afraid of it. Harold Hotelling argued in 1943 that rounding errors in elimination could grow like $4^n$, which would make a $100\times 100$ system hopeless. Von Neumann and Goldstine (1947) showed the fear was unfounded for symmetric positive definite matrices, and James Wilkinson in the early 1960s gave the general backward-error analysis that explains why elimination *with row swaps* works in practice even though a pathological worst case exists. That story, and its vocabulary (growth factor, backward stability, conditioning), is LA-13. Here you will see the single example that makes the whole subject necessary: a $2\times 2$ matrix on which elimination without swaps returns an answer that is wrong in the first digit.

So the concrete questions this module answers: given $A$ and $\mathbf{b}$, how do I actually find $\mathbf{x}$, how do I know in advance whether there is one or many or none, how do I do it a thousand times for a thousand right-hand sides without paying a thousand times, and when does the arithmetic lie to me.

## 2. Intuition

### Subtracting equations doesn't change the answer

Take two equations in two unknowns:

$$
\begin{aligned}
2x + y &= 5 \\
4x + 5y &= 13.
\end{aligned}
$$

Subtract $2$ times the first from the second: $4x + 5y - 2(2x + y) = 13 - 10$, i.e. $3y = 3$. The new pair "$2x + y = 5$, $3y = 3$" has exactly the same solutions as the old pair. Why? Anything that satisfies both old equations satisfies the new second equation (it's built from them). And the operation is reversible: add $2$ times the first equation back and you recover the old second one. Same solution set, both directions.

In the row picture, you replaced the line $4x + 5y = 13$ by the line $3y = 3$. Both pass through the solution $(2, 1)$. You rotated one line about the solution point until it was horizontal. Elimination is a sequence of such rotations: every row operation moves one hyperplane, but all of them keep passing through the common intersection, so the intersection never changes.

Now the new system is **triangular**: the last equation has one unknown. Read it off ($y = 1$), substitute upward ($2x + 1 = 5$, $x = 2$). That is **back substitution**. The entire strategy of this module is: *make it triangular, then substitute.*

### The staircase

For a $3\times 3$ system, eliminate below the first entry of column 1 (two operations), then below the second entry of column 2 (one operation). Three operations, and you have

```
[ ■  *  * ]
[ 0  ■  * ]      ■ = pivot (nonzero)
[ 0  0  ■ ]      * = anything
```

For a non-square or a degenerate matrix the picture is a staircase whose steps can be wider than one column:

```
[ ■  *  *  *  * ]
[ 0  0  ■  *  * ]      columns 1, 3, 4 have pivots
[ 0  0  0  ■  * ]      columns 2, 5 do not ("free")
[ 0  0  0  0  0 ]      a zero row
```

That is a **row echelon form**. The pivots are the corners of the stairs. The number of pivots is going to be the most important number attached to a matrix.

### Where it breaks, and what breaking means

Elimination wants to divide by the pivot. Two things can go wrong.

*Case 1: the pivot position holds a zero, but something below it doesn't.* Swap two rows (swapping two equations changes nothing) and continue. This is a nuisance, not a defect.

*Case 2: the pivot position and everything below it are zero.* No swap helps. That column contributes no step to the staircase; it is a free column. This is real information: it says that column is a combination of the columns before it. In the column picture of LA-1, one of your ingredients was redundant.

Push this to the end and you get the three outcomes you met in LA-1 Example 4.2, now with a mechanism. Suppose you eliminate on the augmented matrix $[A\,|\,\mathbf{b}]$. If you end with a row that reads $0\ 0\ 0\,|\,c$ with $c\ne 0$, the system says $0 = c$: no solution. Otherwise there is a solution, and it is unique exactly when every column of $A$ got a pivot; if a column is free, its unknown can be anything and there are infinitely many solutions.

### The record of the operations is a factorization

Here is the idea you should try to guess before §3 proves it. Each elimination step is "row $i$ ← row $i$ − $\ell_{ij}$ · row $j$." At the end you have $U$ and a list of multipliers $\ell_{ij}$. Run the steps backwards: row $i$ of $A$ equals row $i$ of $U$ plus $\ell_{ij}$ times row $j$ of $U$ for each $j$ you subtracted. In matrix language: $A = LU$ where $L$ has ones on the diagonal and the $\ell_{ij}$ below it. **The multipliers are the factorization.** Nothing more needs to be computed.

Why care? Because solving $A\mathbf{x}=\mathbf{b}$ then splits into two triangular solves: $L\mathbf{c}=\mathbf{b}$ (forward), $U\mathbf{x}=\mathbf{c}$ (backward). The expensive part, producing $L$ and $U$, depends only on $A$. If you have a thousand right-hand sides (a thousand load cases on one structure, a thousand time steps of one implicit ODE solver), you factor once and substitute a thousand times.

### Guess the cost

Elimination has about $n$ steps. Each step updates roughly an $n\times n$ block of entries, each update being a multiply and a subtract. So about $n\cdot n^2 = n^3$ flops, and the constant will turn out to be $2/3$. A triangular solve touches each entry of $U$ once: about $n^2$. So a factorization costs $n$ triangular solves. Keep that ratio in mind when someone suggests computing $A^{-1}$.

### What an inverse is, and when it exists

$A^{-1}$ is the matrix that undoes $A$: $A^{-1}A = I$. If you can solve $A\mathbf{x} = \mathbf{e}_j$ for every $j$, the solutions are the columns of $A^{-1}$ (because $AA^{-1} = I$ says column $j$ of $A^{-1}$ is sent to $\mathbf{e}_j$). Solving $n$ systems with the same $A$ is exactly the "many right-hand sides" situation, which is why Gauss–Jordan on $[A\,|\,I]$ works: it eliminates all $n$ systems at once.

And when does the inverse exist? If elimination finds $n$ pivots, every system $A\mathbf{x}=\mathbf{b}$ has exactly one solution, so the columns of $A^{-1}$ exist. If it finds fewer than $n$, there is a free column, so $A\mathbf{x}=\mathbf{0}$ has a nonzero solution $\mathbf{x}$; but then $A^{-1}$ would have to send $\mathbf{0}$ to both $\mathbf{0}$ and $\mathbf{x}$, which no matrix does. So: **$A$ is invertible if and only if elimination produces $n$ pivots.** That is the main theorem. §3 proves it.

## 3. Formal development

### 3.1 Row operations are left multiplications

There are three row operations:

- **(R1) Replacement.** Row $i$ ← row $i$ − $\ell$ · row $j$, for $i \ne j$.
- **(R2) Swap.** Exchange rows $i$ and $j$.
- **(R3) Scaling.** Row $i$ ← $c$ · row $i$, for $c \ne 0$.

Each is performed by multiplying $A$ *on the left* by a square matrix. The reason is a one-line consequence of the row-times-column rule.

**Theorem 3.1 (Left multiplication acts on rows; right multiplication acts on columns).** For $E$ ($m\times m$) and $A$ ($m\times n$), row $i$ of $EA$ is $\sum_k E_{ik}\,(\text{row } k \text{ of } A)$: a combination of the rows of $A$ with coefficients from row $i$ of $E$. For $F$ ($n\times n$), column $j$ of $AF$ is $\sum_k F_{kj}\,(\text{column } k \text{ of } A)$.

*Proof.* $(EA)_{ij} = \sum_k E_{ik}A_{kj}$ (LA-1 eq. 3.3). Fix $i$ and let $j$ run over the columns: that is the row vector $\sum_k E_{ik}\mathbf{r}_k^\top$. The column statement is LA-1 eq. 3.1 applied to each column of $F$. ∎

**Definition 3.2 (Elementary matrices).** For $i\ne j$ and scalars $\ell$, $c\ne 0$:

$$
E_{ij}(\ell) = I - \ell\,\mathbf{e}_i\mathbf{e}_j^\top, \qquad
P_{ij} = I - \mathbf{e}_i\mathbf{e}_i^\top - \mathbf{e}_j\mathbf{e}_j^\top + \mathbf{e}_i\mathbf{e}_j^\top + \mathbf{e}_j\mathbf{e}_i^\top, \qquad
D_i(c) = I + (c-1)\,\mathbf{e}_i\mathbf{e}_i^\top. \tag{3.1}
$$

Read them with Theorem 3.1. Row $i$ of $E_{ij}(\ell)$ is $\mathbf{e}_i^\top - \ell\mathbf{e}_j^\top$, so row $i$ of $E_{ij}(\ell)A$ is (row $i$ of $A$) − $\ell$(row $j$ of $A$), and every other row of $E_{ij}(\ell)$ is a row of $I$, so every other row of $A$ is untouched. That is (R1). $P_{ij}$ is $I$ with rows $i$ and $j$ exchanged, so $P_{ij}A$ exchanges rows $i$ and $j$ of $A$: (R2). $D_i(c)$ is $I$ with the $(i,i)$ entry replaced by $c$: (R3). For $3\times 3$, with $\ell = 2$:

$$
E_{21}(2) = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \qquad
P_{23} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}, \qquad
D_2(5) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 5 & 0 \\ 0 & 0 & 1 \end{pmatrix}.
$$

Note the sign: $E_{21}(2)$ *subtracts* $2\times$ row 1 from row 2 and stores $-2$. This is the single most common bookkeeping error in the module; see §8.

▶ **Try it first.** Each row operation can be undone by a row operation of the same type. Guess the inverse of each matrix in (3.1), then verify by multiplying. You will need the product rule for the "unit matrices" $\mathbf{e}_i\mathbf{e}_j^\top$.

**Derivation (inverses of elementary matrices).**

1. For any indices, $(\mathbf{e}_i\mathbf{e}_j^\top)(\mathbf{e}_k\mathbf{e}_l^\top) = \mathbf{e}_i(\mathbf{e}_j^\top\mathbf{e}_k)\mathbf{e}_l^\top = \delta_{jk}\,\mathbf{e}_i\mathbf{e}_l^\top$, where $\delta_{jk}$ is $1$ if $j = k$ and $0$ otherwise, because $\mathbf{e}_j^\top\mathbf{e}_k = \mathbf{e}_j\cdot\mathbf{e}_k$. $\tag{3.2}$
2. $E_{ij}(\ell)E_{ij}(-\ell) = (I - \ell\mathbf{e}_i\mathbf{e}_j^\top)(I + \ell\mathbf{e}_i\mathbf{e}_j^\top) = I + \ell\mathbf{e}_i\mathbf{e}_j^\top - \ell\mathbf{e}_i\mathbf{e}_j^\top - \ell^2(\mathbf{e}_i\mathbf{e}_j^\top)(\mathbf{e}_i\mathbf{e}_j^\top)$. By (3.2) the last product is $\delta_{ji}\mathbf{e}_i\mathbf{e}_j^\top = 0$ since $i\ne j$. So the product is $I$. The same computation in the other order gives $E_{ij}(-\ell)E_{ij}(\ell) = I$. Hence $E_{ij}(\ell)^{-1} = E_{ij}(-\ell)$: *to undo "subtract $\ell$ times row $j$," add $\ell$ times row $j$.*
3. $P_{ij}P_{ij}$ swaps rows $i, j$ of $P_{ij}$, which gives back $I$ (swap twice). So $P_{ij}^{-1} = P_{ij}$.
4. $D_i(c)D_i(1/c)$ scales row $i$ of $I$ by $c$ then by $1/c$: $I$. So $D_i(c)^{-1} = D_i(1/c)$. ∎

■ **What just happened.** Every row operation is invertible, and its inverse is a row operation of the same type. That is what makes the next theorem true, and it is the only property of row operations that the rest of the module uses.

**Theorem 3.3 (Row operations preserve the solution set).** If $E$ is invertible, then $A\mathbf{x} = \mathbf{b}$ if and only if $(EA)\mathbf{x} = E\mathbf{b}$.

*Proof.* If $A\mathbf{x}=\mathbf{b}$, multiply on the left by $E$: $EA\mathbf{x} = E\mathbf{b}$. If $EA\mathbf{x} = E\mathbf{b}$, multiply on the left by $E^{-1}$: $A\mathbf{x} = \mathbf{b}$. ∎

A product of invertible matrices is invertible: $(E_2E_1)(E_1^{-1}E_2^{-1}) = E_2(E_1E_1^{-1})E_2^{-1} = I$, and the same in the other order; the inverse of a product is the product of the inverses *in reverse order*. So any sequence of row operations, applied to the augmented matrix $[A\,|\,\mathbf{b}]$, produces a system with exactly the same solutions.

### 3.2 Echelon forms and substitution

**Definition 3.4 (Row echelon form, pivots).** A matrix $U$ is in **row echelon form** if (i) every all-zero row is below every nonzero row, and (ii) the first nonzero entry of each nonzero row, its **pivot**, lies strictly to the right of the pivot of the row above. A column containing a pivot is a **pivot column**; the others are **free columns**. $U$ is in **reduced row echelon form (RREF)** if in addition every pivot equals $1$ and is the only nonzero entry in its column.

**Algorithm 3.5 (Elimination to echelon form).** Start with pivot row $k = 1$. For each column $j = 1, \dots, n$ in turn: look at the entries in column $j$ in rows $k, k+1, \dots, m$. If all are zero, column $j$ is free; go to the next column. Otherwise choose one nonzero entry (by hand: any; in floating point: the largest in absolute value, §3.6), swap its row into row $k$ (R2), and for each $i > k$ compute the multiplier $\ell_{ik} = a_{ik}/a_{kk}$ and apply (R1): row $i$ ← row $i$ − $\ell_{ik}$ · row $k$. Column $j$ now has its pivot in row $k$ and zeros below. Set $k \leftarrow k+1$ and go to the next column. Stop when the columns or the rows run out.

The output is in echelon form: each pivot is placed in a later column than the previous one (we moved right), and every row below the last pivot row is zero (no nonzero entry was found for it). By Theorem 3.3, $U = E_s\cdots E_1 A$ for elementary $E$'s, and $[U\,|\,\mathbf{c}]$ with $\mathbf{c} = E_s\cdots E_1\mathbf{b}$ has the same solutions as $[A\,|\,\mathbf{b}]$.

For a square $U$ with all $n$ pivots on the diagonal, **back substitution** solves $U\mathbf{x} = \mathbf{c}$:

$$
x_n = \frac{c_n}{u_{nn}}, \qquad x_i = \frac{1}{u_{ii}}\Big(c_i - \sum_{j=i+1}^{n} u_{ij}x_j\Big), \quad i = n-1, \dots, 1. \tag{3.3}
$$

For a unit lower triangular $L$ (ones on the diagonal), **forward substitution** solves $L\mathbf{c} = \mathbf{b}$:

$$
c_1 = b_1, \qquad c_i = b_i - \sum_{j=1}^{i-1}\ell_{ij}c_j, \quad i = 2, \dots, n. \tag{3.4}
$$

**Theorem 3.6 (Triangular systems).** If $U$ is upper triangular with nonzero diagonal, $U\mathbf{x} = \mathbf{c}$ has exactly one solution, given by (3.3). Likewise $L\mathbf{c} = \mathbf{b}$ with $L$ unit lower triangular has exactly one solution, (3.4).

*Proof.* Equation $i$ of $U\mathbf{x}=\mathbf{c}$ reads $u_{ii}x_i + \sum_{j>i}u_{ij}x_j = c_i$; it involves only $x_i, \dots, x_n$. Equation $n$ forces $x_n = c_n/u_{nn}$ (division allowed because $u_{nn}\ne 0$). Given $x_n$, equation $n-1$ forces $x_{n-1}$, and so on up: at each stage the later unknowns are already fixed and the equation is solved for $x_i$ by dividing by $u_{ii}\ne 0$. Any solution must take these values (uniqueness), and these values satisfy every equation (existence). The lower-triangular case is the same argument from the top down, with no divisions. ∎

### 3.3 Existence and uniqueness, read off the echelon form

**Theorem 3.7 (The three outcomes).** Let $A$ be $m\times n$ and let $[U\,|\,\mathbf{c}]$ be any echelon form of $[A\,|\,\mathbf{b}]$, with $r$ pivots in the $U$ part.

(a) $A\mathbf{x} = \mathbf{b}$ has **no solution** if and only if some row of $[U\,|\,\mathbf{c}]$ is zero in the $U$ part and nonzero in the $\mathbf{c}$ part (equivalently: the augmented column is a pivot column).

(b) If (a) does not occur and $r = n$ (every column of $A$ is a pivot column), there is **exactly one** solution.

(c) If (a) does not occur and $r < n$, there are **infinitely many**. The $n - r$ free variables may be given any values; the pivot variables are then determined uniquely. Every solution has the form $\mathbf{x}_p + \sum_{f\ \mathrm{free}} t_f\,\mathbf{s}_f$, where $\mathbf{x}_p$ is the solution with all free variables zero and $\mathbf{s}_f$ (the **special solution** for free column $f$) solves $U\mathbf{x} = \mathbf{0}$ with $x_f = 1$ and the other free variables $0$.

▶ **Try it first.** Convince yourself of (a) in one sentence, then try to *construct* a solution when (a) doesn't happen: choose the free variables, then work up through the pivot rows.

**Derivation.**

1. By Theorem 3.3 the two systems have the same solutions, so we may work with $[U\,|\,\mathbf{c}]$.
2. A row $(0\ \cdots\ 0\,|\,c_i)$ with $c_i \ne 0$ is the equation $0 = c_i$. Nothing satisfies it. That proves the "if" of (a).
3. Suppose no such row exists. Assign arbitrary values to the free variables. Now read the pivot rows from the bottom up. Pivot row $k$ has its pivot $u_{k,p_k}\ne 0$ in column $p_k$ and zeros to the left, so it reads $u_{k,p_k}x_{p_k} + (\text{terms in } x_j,\ j > p_k) = c_k$. Every $x_j$ with $j > p_k$ is either a free variable (assigned) or a pivot variable belonging to a lower pivot row (already determined, because we work upward). Solve for $x_{p_k}$ by dividing by $u_{k,p_k}$. Zero rows read $0 = 0$. We have produced a solution, which proves the "only if" of (a).
4. In case (b) there are no free variables, so step 3 involves no choices: each pivot variable is forced. Any solution must satisfy the pivot equations, so it must equal the one we built. Unique.
5. In case (c) pick a free variable $x_f$. Two assignments differing in $x_f$ yield two solutions differing in coordinate $f$; there are infinitely many assignments. For the structure claim: if $\mathbf{x}$ is any solution, then $\mathbf{x} - \mathbf{x}_p$ solves $U\mathbf{x} = \mathbf{0}$ (subtract the equations) and has free coordinates $t_f = x_f$. By step 3 applied to $U\mathbf{x} = \mathbf{0}$, a solution of the homogeneous system is determined by its free coordinates; $\sum_f t_f\mathbf{s}_f$ is a solution (combinations of solutions of a homogeneous system are solutions, by linearity, LA-1 Theorem 3.4) with the same free coordinates. So $\mathbf{x} - \mathbf{x}_p = \sum_f t_f\mathbf{s}_f$. ∎

■ **What just happened.** The three outcomes from LA-1 Example 4.2 and exercise B1 now have a mechanism. (a) is "$\mathbf{b}$ is not in the span of the columns": elimination exposed a combination of the equations that reads $0 = c$. (b) is "no column is redundant." (c) is "some column is a combination of earlier ones," and the special solutions $\mathbf{s}_f$ write down exactly which combination. In LA-3 the set $\{\sum t_f\mathbf{s}_f\}$ gets its name: the **null space** $N(A)$, and $n - r$ becomes its dimension.

### 3.4 The theorem this module exists to prove

**Theorem 3.8 (Fewer than $n$ pivots gives a nonzero null vector).** Let $A$ be $n\times n$. If elimination produces fewer than $n$ pivots, then $A\mathbf{x} = \mathbf{0}$ has a solution $\mathbf{x} \ne \mathbf{0}$.

*Proof.*

1. Elimination gives $EA = U$ with $E = E_s\cdots E_1$ a product of elementary matrices, hence invertible (§3.1).
2. $U$ has $r < n$ pivots in $n$ columns, so at least one column $f$ is free.
3. The system $U\mathbf{x} = \mathbf{0}$ is consistent (its right side is $\mathbf{0}$, so outcome (a) cannot occur). Apply Theorem 3.7(c): set $x_f = 1$, set the other free variables to $0$, and determine the pivot variables by back substitution. The result $\mathbf{x}$ satisfies $U\mathbf{x} = \mathbf{0}$ and $x_f = 1$, so $\mathbf{x} \ne \mathbf{0}$.
4. By Theorem 3.3, $A\mathbf{x} = E^{-1}U\mathbf{x} = \mathbf{0}$. ∎

**Corollary 3.9 (More unknowns than equations).** If $A$ is $m\times n$ with $n > m$, then $A\mathbf{x} = \mathbf{0}$ has a nonzero solution. (There are at most $m$ pivots, one per row, so some column is free; the same construction works. Exercise B3 asks you to write it out.)

### 3.5 $A = LU$: the record of elimination

Assume for now that elimination on the square matrix $A$ never needs a row swap: at each step the pivot position is already nonzero. (§3.6 removes this assumption.) Then the operations are $E_{ij} := E_{ij}(\ell_{ij})$ for $i > j$, applied column by column:

$$
E_{n,n-1}\cdots E_{32}\,E_{n1}\cdots E_{31}E_{21}\,A = U,
\qquad\text{so}\qquad
A = E_{21}^{-1}E_{31}^{-1}\cdots E_{n1}^{-1}\,E_{32}^{-1}\cdots E_{n,n-1}^{-1}\,U =: LU. \tag{3.5}
$$

**Theorem 3.10 ($L$ stores the multipliers).** The matrix $L$ in (3.5) is unit lower triangular with $L_{ij} = \ell_{ij}$ for $i > j$: the entry below the diagonal in position $(i,j)$ is exactly the multiplier that was used to clear position $(i,j)$. No arithmetic beyond the elimination itself is needed to form $L$.

▶ **Try it first.** Take the $3\times 3$ case. Row 1 of $U$ is row 1 of $A$ (it was never changed). What happened to row 3 of $A$ on its way to becoming row 3 of $U$? Write row 3 of $A$ as a combination of rows of $U$.

**Derivation (first way: run elimination backwards).**

1. Row $j$ of $U$ is the state of row $j$ at the moment it was used as a pivot row; after that moment row $j$ is never changed again (we only modify rows *below* the pivot row). So "row $j$ of $U$" is well defined as "the pivot row $j$," call it $\mathbf{u}_j^\top$.
2. Row $i$ of $A$ is changed exactly $i-1$ times: at step $j = 1, 2, \dots, i-1$ we subtract $\ell_{ij}\mathbf{u}_j^\top$ from it. After the last of these it is $\mathbf{u}_i^\top$. So
   $$\mathbf{u}_i^\top = (\text{row } i \text{ of } A) - \ell_{i1}\mathbf{u}_1^\top - \ell_{i2}\mathbf{u}_2^\top - \cdots - \ell_{i,i-1}\mathbf{u}_{i-1}^\top.$$
3. Rearranging: $(\text{row } i \text{ of } A) = \ell_{i1}\mathbf{u}_1^\top + \cdots + \ell_{i,i-1}\mathbf{u}_{i-1}^\top + 1\cdot\mathbf{u}_i^\top$.
4. By Theorem 3.1, "row $i$ of $A$ is a combination of the rows of $U$ with coefficients $(\ell_{i1}, \dots, \ell_{i,i-1}, 1, 0, \dots, 0)$" says exactly that $A = LU$ where row $i$ of $L$ is $(\ell_{i1}, \dots, \ell_{i,i-1}, 1, 0, \dots, 0)$. That $L$ is unit lower triangular with the multipliers below the diagonal. ∎

**Derivation (second way: multiply the inverses).** This explains why the *inverse* product in (3.5) is clean while the forward product $E_{32}E_{31}E_{21}$ is not.

1. $E_{ij}^{-1} = I + \ell_{ij}\mathbf{e}_i\mathbf{e}_j^\top$ (§3.1 step 2).
2. For $n = 3$: $(I + \ell_{21}\mathbf{e}_2\mathbf{e}_1^\top)(I + \ell_{31}\mathbf{e}_3\mathbf{e}_1^\top) = I + \ell_{21}\mathbf{e}_2\mathbf{e}_1^\top + \ell_{31}\mathbf{e}_3\mathbf{e}_1^\top + \ell_{21}\ell_{31}\,\delta_{13}\,\mathbf{e}_2\mathbf{e}_1^\top$ by (3.2), and $\delta_{13} = 0$. Multiplying by $(I + \ell_{32}\mathbf{e}_3\mathbf{e}_2^\top)$ produces cross terms with $\delta_{13}$ again (from $\mathbf{e}_2\mathbf{e}_1^\top\cdot\mathbf{e}_3\mathbf{e}_2^\top$ and $\mathbf{e}_3\mathbf{e}_1^\top\cdot\mathbf{e}_3\mathbf{e}_2^\top$), both zero. So $L = I + \ell_{21}\mathbf{e}_2\mathbf{e}_1^\top + \ell_{31}\mathbf{e}_3\mathbf{e}_1^\top + \ell_{32}\mathbf{e}_3\mathbf{e}_2^\top$: the multipliers, overlaid.
3. In general, the factors in (3.5) are ordered so that a factor $(I + \ell_{ij}\mathbf{e}_i\mathbf{e}_j^\top)$ stands to the left only of factors $(I + \ell_{pq}\mathbf{e}_p\mathbf{e}_q^\top)$ with $q \ge j$ and $p > q$. The cross term is $\ell_{ij}\ell_{pq}\,\delta_{jp}\,\mathbf{e}_i\mathbf{e}_q^\top$, and $p > q \ge j$ forces $p \ne j$. Multiplying in from the right, by induction, no cross term ever survives. ∎
4. By contrast, in $E_{32}E_{31}E_{21}$ the factor $(I - \ell_{32}\mathbf{e}_3\mathbf{e}_2^\top)$ stands to the left of $(I - \ell_{21}\mathbf{e}_2\mathbf{e}_1^\top)$, and the cross term $\ell_{32}\ell_{21}\,\delta_{22}\,\mathbf{e}_3\mathbf{e}_1^\top$ survives: the $(3,1)$ entry of $E_{32}E_{31}E_{21}$ is $-\ell_{31} + \ell_{32}\ell_{21}$, not $-\ell_{31}$. Example 4.1 computes both and the project asserts the difference.

■ **What just happened.** $A = LU$ is not a new computation; it is elimination *remembered*. In code, you overwrite each zeroed entry $a_{ij}$ with $\ell_{ij}$, and the array holds $L$ (strictly below the diagonal) and $U$ (on and above) at no extra cost. Solving $A\mathbf{x}=\mathbf{b}$ becomes $L\mathbf{c}=\mathbf{b}$ by (3.4), then $U\mathbf{x}=\mathbf{c}$ by (3.3). Forward substitution *is* elimination applied to $\mathbf{b}$: $c_i = b_i - \sum_{j<i}\ell_{ij}c_j$ is what happens to $b_i$ when you subtract $\ell_{ij}$ times the (already processed) $j$-th equation.

Two facts about triangular matrices are needed for uniqueness and for §3.9.

**Lemma 3.11 (Triangular facts).** (i) The product of two lower (upper) triangular matrices is lower (upper) triangular, with diagonal entries the products of the diagonal entries. (ii) A triangular matrix is invertible iff its diagonal entries are all nonzero, and the inverse of a unit lower triangular matrix is unit lower triangular.

*Proof.* (i) For lower triangular $A, B$ and $i < j$: $(AB)_{ij} = \sum_k a_{ik}b_{kj}$; each term has $k > i$ (so $a_{ik} = 0$) or $k \le i < j$ (so $b_{kj} = 0$). On the diagonal only $k = i$ survives. Transpose for the upper case. (ii) A triangular matrix *is* its own echelon form (possibly after reordering nothing), so its pivots are its diagonal entries; it has $n$ pivots iff none is zero, and Theorem 3.16 below will make "$n$ pivots" equivalent to invertibility. For the unit lower case, solve $LX = I$ column by column: column $j$ of $X$ solves $L\mathbf{x} = \mathbf{e}_j$, and forward substitution (3.4) gives $x_1 = \cdots = x_{j-1} = 0$ (the right side is zero there and the earlier $x$'s are zero) and $x_j = 1$. So $X$ is unit lower triangular, and $LX = I$ makes $X = L^{-1}$ by Theorem 3.15 below. ∎

**Theorem 3.12 (Uniqueness of $LU$).** If $A$ is invertible and $A = L_1U_1 = L_2U_2$ with $L_1, L_2$ unit lower triangular and $U_1, U_2$ upper triangular, then $L_1 = L_2$ and $U_1 = U_2$.

*Proof.* $U_1 = L_1^{-1}A$ is invertible (product of invertibles), likewise $U_2$. From $L_1U_1 = L_2U_2$: $L_2^{-1}L_1 = U_2U_1^{-1}$. The left side is unit lower triangular (Lemma 3.11), the right side is upper triangular. A matrix that is both is diagonal with ones on the diagonal: $I$. So $L_1 = L_2$, and then $U_1 = U_2$. ∎

### 3.6 Row swaps: $PA = LU$ and why pivoting is not optional

Two reasons to swap rows. The first is necessity: $A = \begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix}$ is invertible but has a zero in the pivot position, so $A = LU$ with unit lower $L$ is impossible (the $(1,1)$ entry of $LU$ is $u_{11}$, and then $(2,1)$ is $\ell_{21}u_{11}$, which can't be $1$ when $u_{11} = 0$). The second is arithmetic, and it is the one that matters in practice.

**The $10^{-20}$ example.** Let $\varepsilon = 10^{-20}$ and solve

$$
\begin{pmatrix} \varepsilon & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2\end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}.
$$

The exact solution is $x_1 = 1/(1-\varepsilon) \approx 1$, $x_2 = 1 - \varepsilon x_1 \approx 1$: both entries are $1$ to twenty digits. Eliminate without swapping, in IEEE double precision (about 16 significant digits). The multiplier is $\ell_{21} = 1/\varepsilon = 10^{20}$. The new $(2,2)$ entry is $1 - 10^{20}\cdot 1$, and the computer rounds $1 - 10^{20}$ to $-10^{20}$: the $1$ is twenty orders of magnitude below the rounding threshold of $10^{20}$ and vanishes. The right side becomes $2 - 10^{20} \to -10^{20}$ likewise. Back substitution: $x_2 = (-10^{20})/(-10^{20}) = 1$ (correct, by luck), then $x_1 = (1 - x_2)/\varepsilon = 0/\varepsilon = 0$. **Wrong in the first digit.** Check the factorization: the computed $L$ and $U$ multiply to $\begin{pmatrix} \varepsilon & 1 \\ 1 & 0 \end{pmatrix}$. The $(2,2)$ entry of $A$, which was $1$, has been replaced by $0$. The algorithm solved a different matrix.

Swap the rows first: $\begin{pmatrix} 1 & 1 \\ \varepsilon & 1 \end{pmatrix}$, multiplier $\ell_{21} = \varepsilon$, new $(2,2)$ entry $1 - \varepsilon \to 1$, right side $1 - 2\varepsilon \to 1$. Then $x_2 = 1$, $x_1 = 2 - 1 = 1$. Correct to all digits.

The lesson is not "small pivots are bad numbers." The matrix is perfectly nice. The lesson is that **a large multiplier amplifies rounding error**: subtracting $\ell\times(\text{row } j)$ with $|\ell| \gg 1$ swamps the information in row $i$. **Partial pivoting** is the rule: at step $k$, choose as pivot the entry of largest absolute value in column $k$ on or below the diagonal. Then every multiplier satisfies $|\ell_{ik}| = |a_{ik}/a_{kk}| \le 1$ and no entry of $U$ can grow by more than a factor of $2$ per step. Every production library (LAPACK's `getrf`, which `np.linalg.solve` and `scipy.linalg.lu` call) does this.

**Theorem 3.13 ($PA = LU$).** Elimination with row swaps on an invertible $A$ produces a permutation matrix $P$, a unit lower triangular $L$, and an upper triangular $U$ with

$$
PA = LU. \tag{3.6}
$$

Under partial pivoting, $|L_{ij}| \le 1$ for all $i, j$. Here $P$ records the net reordering of the rows of $A$.

▶ **Try it first.** For $n = 3$ the operations are $E_2P_2E_1P_1A = U$, where $P_1$ swaps row 1 with some row, $E_1$ clears column 1, $P_2$ swaps rows 2 and 3 (or does nothing), $E_2$ clears position $(3,2)$. The problem is the $P_2$ sitting between $E_1$ and $E_2$. Show that $P_2E_1 = E_1'P_2$ where $E_1'$ is $E_1$ with its two multipliers swapped.

**Derivation.**

1. Write $E_1 = I - \ell_{21}\mathbf{e}_2\mathbf{e}_1^\top - \ell_{31}\mathbf{e}_3\mathbf{e}_1^\top$ and $P_2 = P_{23}$. Since $P_2^2 = I$, $P_2E_1 = (P_2E_1P_2)P_2$. Compute $E_1' := P_2E_1P_2 = P_2P_2 - \sum_{m=2,3}\ell_{m1}(P_2\mathbf{e}_m)(\mathbf{e}_1^\top P_2)$.
2. $P_2\mathbf{e}_2 = \mathbf{e}_3$ and $P_2\mathbf{e}_3 = \mathbf{e}_2$ (a permutation matrix sends $\mathbf{e}_m$ to $\mathbf{e}_{\pi(m)}$). And $\mathbf{e}_1^\top P_2 = (P_2^\top\mathbf{e}_1)^\top = (P_2\mathbf{e}_1)^\top = \mathbf{e}_1^\top$, since $P_{23}$ is symmetric and fixes index $1$. So $E_1' = I - \ell_{21}\mathbf{e}_3\mathbf{e}_1^\top - \ell_{31}\mathbf{e}_2\mathbf{e}_1^\top$: the same two multipliers, in swapped rows. Still unit lower triangular.
3. Therefore $U = E_2P_2E_1P_1A = E_2E_1'(P_2P_1)A$. Set $P = P_2P_1$ (a product of permutation matrices is a permutation matrix, §3.9) and $L = (E_2E_1')^{-1} = E_1'^{-1}E_2^{-1}$, which by Theorem 3.10's second derivation is unit lower triangular holding the (reordered) multipliers. Then $PA = LU$.
4. For general $n$ the same identity $P_kE_j = E_j'P_k$ holds whenever $P_k$ swaps two rows both below row $j$, which is always the case for a swap performed at a later step $k > j$: $P_k$ fixes indices $1, \dots, k-1 \supseteq \{j\}$, so $\mathbf{e}_j^\top P_k = \mathbf{e}_j^\top$, and $P_k$ permutes the column-$j$ multipliers among rows $k, \dots, n$. Moving every $P_k$ to the right one at a time gives $E_{n-1}'\cdots E_1'\,P_{n-1}\cdots P_1\,A = U$. ∎

■ **What just happened.** In code, this is one line: *when you swap rows $k$ and $p$ of the working matrix, also swap rows $k$ and $p$ in the already-computed columns $1, \dots, k-1$ of $L$.* That is exactly $E_j \to E_j'$. Forgetting it is the standard bug in a hand-written LU; the project's asserts against `scipy.linalg.lu` will catch it.

Honesty about what partial pivoting guarantees: the bound "$|\ell| \le 1$" allows entries of $U$ to grow by at most $2^{n-1}$ overall, and there are matrices (Wilkinson's, exercise C1) that achieve it. So partial pivoting is not provably stable in the worst case. It is used universally because, for matrices that arise in practice, the growth is almost always small (typically a modest multiple of $1$, not $2^{n-1}$). Why that is so is a research-level question; LA-13 states what is known (Trefethen & Bau lecture 22).

### 3.7 Operation counts

▶ **Try it first.** At step $k$ of elimination on an $n\times n$ matrix, how many rows are modified, and how many entries in each? Count one division per multiplier, and one multiply plus one subtract per updated entry. Don't count the entries you already know will be zero.

**Derivation (cost of $LU$).**

1. At step $k$ ($k = 1, \dots, n-1$) there are $m = n - k$ rows below the pivot row. Each needs one division to form $\ell_{ik}$: $m$ flops.
2. Each of those rows has $m$ entries to the right of column $k$ that change (columns $k+1, \dots, n$); column $k$ itself we *know* becomes zero and don't compute. Each update $a_{ij} \leftarrow a_{ij} - \ell_{ik}a_{kj}$ is one multiply and one subtract: $2m^2$ flops for the step.
3. Total: $\displaystyle\sum_{m=1}^{n-1}(m + 2m^2) = \frac{n(n-1)}{2} + \frac{(n-1)n(2n-1)}{3}$, using $\sum_{m=1}^{N}m^2 = N(N+1)(2N+1)/6$ with $N = n-1$. The leading term is
   $$\text{flops}(LU) \approx \tfrac{2}{3}n^3. \tag{3.7}$$

**Derivation (cost of the substitutions).** Forward substitution (3.4): row $i$ has $i - 1$ multiplies and $i-1$ subtracts, total $\sum_{i=1}^n 2(i-1) = n(n-1) \approx n^2$. Back substitution (3.3): row $i$ has $n - i$ multiplies, $n-i$ subtracts, one divide: $\sum_i (2(n-i) + 1) = n(n-1) + n = n^2$ exactly. So one solve, given $L$ and $U$, costs about $2n^2$ flops.

■ **What just happened.** Three consequences, each of which you should be able to recite.

- **Factor once, solve many.** $k$ right-hand sides cost $\tfrac23 n^3 + 2kn^2$ with one $LU$, versus $k\cdot\tfrac23 n^3$ from scratch. For $n = k = 1000$: about $2.7\times 10^9$ flops versus $6.7\times 10^{11}$, a factor of $250$. The project times this.
- **The inverse costs three factorizations.** Computing $A^{-1}$ means solving $n$ systems $A\mathbf{x} = \mathbf{e}_j$. Done naively that is $\tfrac23 n^3 + 2n\cdot n^2 = \tfrac83 n^3$; done carefully (exploiting the zeros in $\mathbf{e}_j$, which Gauss–Jordan does automatically) it is about $2n^3$. Either way: at least three times the cost of $LU$.
- **The inverse saves nothing per solve.** Multiplying $A^{-1}\mathbf{b}$ is $n^2$ multiplies and $n^2$ adds: $2n^2$ flops, the same as the two triangular solves. So $A^{-1}$ costs $3\times$ up front and gains nothing afterward.

There is also an accuracy argument, and it should be stated carefully. The $LU$ route is **backward stable**: the computed $\hat{\mathbf{x}}$ exactly solves a nearby problem, $(A + \delta A)\hat{\mathbf{x}} = \mathbf{b}$ with $\|\delta A\|$ about machine precision times $\|A\|$ times the growth factor (Wilkinson; LA-13). So the residual $\mathbf{b} - A\hat{\mathbf{x}}$ is tiny relative to $\|A\|\|\hat{\mathbf{x}}\|$ regardless of how ill-conditioned $A$ is. The route $\hat{\mathbf{x}} = \widehat{A^{-1}}\mathbf{b}$ has no such guarantee; the best general bound on its residual carries an extra factor of the condition number $\kappa(A) = \|A\|\|A^{-1}\|$ (LA-9, LA-13). In practice, with a good library in double precision, the two residuals are often within a small factor of each other; the project measures them and you should look at the numbers rather than take the folklore on faith. The arguments that *always* hold are the cost argument above and a structural one: $A^{-1}$ is dense even when $A$ is sparse. The tridiagonal matrix of Example 4.3 has $3n - 2$ nonzeros and an $LU$ with about the same; its inverse has $n^2$ nonzeros, none of them zero. For the $10^6\times 10^6$ systems in EM-3 and CM-6 that is the difference between feasible and impossible. **Never form $A^{-1}$ to solve a system.** Call `solve`, or factor and substitute.

### 3.8 Inverses

**Definition 3.14.** A square matrix $A$ is **invertible** (nonsingular) if there is a matrix $B$ with $AB = BA = I$. Then $B$ is called $A^{-1}$. A square matrix that is not invertible is **singular**.

**Theorem 3.15 (Uniqueness; one side suffices).** If $BA = I$ and $AC = I$, then $B = C$. Consequently a matrix has at most one inverse, and any left inverse of a square matrix equals any right inverse.

*Proof.* $B = BI = B(AC) = (BA)C = IC = C$. ∎

**Theorem 3.16 (Invertible $\iff$ $n$ pivots).** An $n\times n$ matrix $A$ is invertible if and only if elimination produces $n$ pivots.

▶ **Try it first.** One direction is Theorem 3.8 plus one line. For the other, you need to *build* $B$ with $BA = I$: think about which row operations take $U$ (with $n$ nonzero pivots) all the way to $I$.

**Derivation.**

($\Leftarrow$) Suppose there are $n$ pivots, so $EA = U$ with $E$ a product of elementary matrices and $U$ upper triangular with nonzero diagonal.

1. Clear *above* the pivots, last column first. For column $n$: for each $i < n$ apply $E_{in}(u_{in}/u_{nn})$, which subtracts a multiple of row $n$ from row $i$. Row $n$ is $(0, \dots, 0, u_{nn})$, so this changes only column $n$ of the other rows: column $n$ becomes $(0, \dots, 0, u_{nn})^\top$ and nothing to the left is disturbed. Repeat for column $n-1$ using row $n-1$, which is zero to the left of column $n-1$ and now also zero in column $n$. Continue to column $2$. The result is diagonal: $\operatorname{diag}(u_{11}, \dots, u_{nn})$.
2. Apply $D_i(1/u_{ii})$ for each $i$. The result is $I$.
3. So $MA = I$ where $M$ is the product of all elementary matrices used (the original $E$, then the upward clearing, then the scalings). $M$ is invertible as a product of invertible matrices.
4. From $MA = I$: multiply on the left by $M^{-1}$ to get $A = M^{-1}$. Then $AM = M^{-1}M = I$. So $AM = MA = I$ and $A^{-1} = M$.

($\Rightarrow$) Suppose $A$ is invertible but elimination gives fewer than $n$ pivots.

5. By Theorem 3.8 there is $\mathbf{x} \ne \mathbf{0}$ with $A\mathbf{x} = \mathbf{0}$.
6. Multiply by $A^{-1}$: $\mathbf{x} = A^{-1}A\mathbf{x} = A^{-1}\mathbf{0} = \mathbf{0}$. Contradiction. ∎

**Corollary 3.17 (Gauss–Jordan computes the inverse).** Apply row operations to the $n\times 2n$ matrix $[A\,|\,I]$ until the left block is $I$. The right block is then $A^{-1}$:

$$
M[A\,|\,I] = [MA\,|\,MI] = [I\,|\,M] = [I\,|\,A^{-1}]. \tag{3.8}
$$

If the left block cannot be brought to $I$ (fewer than $n$ pivots), $A$ is singular.

*Proof.* Row operations on a block matrix act on both blocks (Theorem 3.1), and the operations that turn $A$ into $I$ are the $M$ of step 3, which is $A^{-1}$ by step 4. ∎

**Corollary 3.18 (Singular $\iff$ nontrivial null vector; unique solutions).** $A$ is singular iff $A\mathbf{x} = \mathbf{0}$ has a nonzero solution. If $A$ is invertible, $A\mathbf{x} = \mathbf{b}$ has exactly one solution for every $\mathbf{b}$, namely $A^{-1}\mathbf{b}$.

*Proof.* Singular $\iff$ fewer than $n$ pivots (Theorem 3.16) $\Rightarrow$ nonzero null vector (Theorem 3.8); conversely a nonzero null vector forces singularity by step 6. For the second claim: $A(A^{-1}\mathbf{b}) = \mathbf{b}$ gives existence; if $A\mathbf{x} = \mathbf{b}$ then $\mathbf{x} = A^{-1}A\mathbf{x} = A^{-1}\mathbf{b}$ gives uniqueness. ∎

**Theorem 3.19 (Algebra of inverses).** If $A$ and $B$ are invertible $n\times n$ matrices, then
(a) $(AB)^{-1} = B^{-1}A^{-1}$; (b) $(A^\top)^{-1} = (A^{-1})^\top$; (c) $(A^{-1})^{-1} = A$; (d) $(cA)^{-1} = c^{-1}A^{-1}$ for $c \ne 0$.

*Proof.* By Theorem 3.15 it suffices to exhibit a one-sided inverse in each case. (a) $(AB)(B^{-1}A^{-1}) = A(BB^{-1})A^{-1} = AA^{-1} = I$. Note the order reversal, the same pattern as $(AB)^\top = B^\top A^\top$: to undo "do $B$ then $A$," undo $A$ first. (b) $(A^{-1})^\top A^\top = (AA^{-1})^\top = I^\top = I$ by LA-1 Theorem 3.8. (c) $AA^{-1} = I$ says $A$ is a right inverse of $A^{-1}$. (d) $(cA)(c^{-1}A^{-1}) = (cc^{-1})AA^{-1} = I$. ∎

The converse of (a) is also true and useful: if $AB$ is invertible then so are $A$ and $B$ (exercise B5).

For $2\times 2$ matrices there is a closed form you should verify once and then never use for anything larger: $\begin{pmatrix} a & b \\ c & d\end{pmatrix}^{-1} = \dfrac{1}{ad - bc}\begin{pmatrix} d & -b \\ -c & a\end{pmatrix}$. Multiply it out. The condition $ad - bc \ne 0$ is "two pivots": if $a \ne 0$ the second pivot is $d - (c/a)b = (ad - bc)/a$. The number $ad - bc$ is the determinant; LA-6 shows the determinant of any square matrix is $\pm$ the product of its pivots.

### 3.9 Rank, provisionally

**Definition 3.20 (Rank, provisional).** The **rank** of $A$, written $r$ or $\operatorname{rank}A$, is the number of pivots produced by elimination.

What we can say now: $r \le \min(m, n)$ (at most one pivot per row and per column). A square matrix is invertible iff $r = n$ (Theorem 3.16). $A\mathbf{x} = \mathbf{b}$ is solvable for *every* $\mathbf{b}$ iff $r = m$ (every row has a pivot, so no row of $U$ is zero, so outcome (a) of Theorem 3.7 cannot occur). The solution set, when nonempty, has $n - r$ free parameters.

What we cannot yet say honestly: that this number is a property of $A$ rather than of the particular sequence of row operations you chose. Different swap choices give different echelon forms. The number of pivots is in fact the same for all of them, and exercise C3 gives a direct proof: a column is a pivot column exactly when it is not a linear combination of the columns to its left, and that property is unchanged by invertible row operations (Theorem 3.3 applied to the homogeneous system). LA-3 does better: it proves rank $= \dim C(A) = \dim$ (row space), a statement with no algorithm in it, and that the RREF of $A$ is unique. Until then, treat "rank" as "pivot count," and be aware that in floating point the count depends on a tolerance (the project shows a matrix of exact rank 2 whose computed third pivot is $10^{-16}$, not $0$).

### 3.10 $LDU$, symmetric matrices, and the special matrices

**Theorem 3.21 ($A = LDU$).** If $A$ is invertible and $A = LU$ with no row swaps, let $D = \operatorname{diag}(u_{11}, \dots, u_{nn})$ be the diagonal matrix of pivots. Then

$$
A = L\,D\,U', \qquad U' = D^{-1}U \text{ unit upper triangular}, \tag{3.9}
$$

and this factorization into (unit lower)·(diagonal)·(unit upper) is unique.

*Proof.* $D$ is invertible (nonzero pivots), $D^{-1}U$ divides row $i$ of $U$ by $u_{ii}$, which puts ones on the diagonal and keeps it upper triangular; $LDU' = LDD^{-1}U = LU = A$. Uniqueness: if $L_1D_1U_1' = L_2D_2U_2'$ then $L_1(D_1U_1') = L_2(D_2U_2')$ are two $LU$ factorizations with unit lower $L$ and upper triangular $D_iU_i'$, so by Theorem 3.12 $L_1 = L_2$ and $D_1U_1' = D_2U_2'$. The diagonal of $D_iU_i'$ is the diagonal of $D_i$ (Lemma 3.11(i)), so $D_1 = D_2$, and then $U_1' = U_2'$. ∎

**Theorem 3.22 (Symmetric $\Rightarrow$ $A = LDL^\top$).** If $A = A^\top$ is invertible and has an $LDU'$ factorization (3.9), then $U' = L^\top$:

$$
A = LDL^\top. \tag{3.10}
$$

▶ **Try it first.** Transpose both sides of $A = LDU'$ and use uniqueness.

*Proof.* $A = A^\top = (LDU')^\top = U'^\top D^\top L^\top = U'^\top D L^\top$ (LA-1 Theorem 3.8 twice; $D^\top = D$). Now $U'^\top$ is unit lower triangular, $D$ is diagonal, $L^\top$ is unit upper triangular. So $U'^\top D L^\top$ is an $LDU'$-type factorization of $A$. By uniqueness, $U'^\top = L$ and $L^\top = U'$. ∎

■ **What just happened.** For a symmetric matrix you only need to compute the lower half of the elimination: $L$ and the pivots determine everything. The work halves to about $n^3/3$ flops. When all pivots are positive, $A = LDL^\top = (LD^{1/2})(LD^{1/2})^\top$ is the **Cholesky factorization**, and "all pivots positive" is one of the equivalent definitions of *positive definite* in LA-8. The project's "Your turn" has you implement (3.10).

**Special matrices, collected.**

- **Diagonal.** Invertible iff every diagonal entry is nonzero; inverse is the diagonal of reciprocals. Multiplying by $D$ on the left scales rows; on the right scales columns (Theorem 3.1).
- **Triangular.** Lemma 3.11: products and inverses stay triangular; invertible iff no zero on the diagonal; the pivots are the diagonal.
- **Permutation.** $P$ has exactly one $1$ in each row and each column, zeros elsewhere; its rows are $\mathbf{e}_{\pi(1)}^\top, \dots, \mathbf{e}_{\pi(n)}^\top$ for a bijection $\pi$ of $\{1, \dots, n\}$, and $PA$ reorders the rows of $A$ by $\pi$. There are $n!$ of them. Products of permutation matrices are permutation matrices ($P_1P_2\mathbf{e}_j = P_1\mathbf{e}_{\pi_2(j)} = \mathbf{e}_{\pi_1(\pi_2(j))}$, and $\pi_1\circ\pi_2$ is a bijection).

**Theorem 3.23.** For every permutation matrix $P$, $P^\top P = PP^\top = I$; that is, $P^{-1} = P^\top$.

*Proof.* $(PP^\top)_{ij} = (\text{row } i \text{ of } P)\cdot(\text{row } j \text{ of } P) = \mathbf{e}_{\pi(i)}\cdot\mathbf{e}_{\pi(j)} = \delta_{\pi(i)\pi(j)} = \delta_{ij}$, the last step because $\pi$ is a bijection. So $PP^\top = I$, and $P^\top = P^{-1}$ by Theorem 3.15; then also $P^\top P = I$. ∎

- **Symmetric.** $A = A^\top$. For any $A$ (even rectangular), $A^\top A$ and $AA^\top$ are symmetric: $(A^\top A)^\top = A^\top(A^\top)^\top = A^\top A$. If $A$ is symmetric and invertible, $A^{-1}$ is symmetric: $(A^{-1})^\top = (A^\top)^{-1} = A^{-1}$ by Theorem 3.19(b). Symmetric matrices are where physics lives: stiffness matrices (CM-6), conductance matrices (Example 4.4), inertia tensors (CM-7), Hessians (ML-2), covariance matrices (F5) are all symmetric, and (3.10) is the first of several factorizations (Cholesky in LA-10, spectral in LA-8) that exploit it.

## 4. Worked examples

### Example 4.1: A full $3\times 3$ elimination with $L$ and $U$

Solve $A\mathbf{x} = \mathbf{b}$ with

$$
A = \begin{pmatrix} 2 & 1 & 1 \\ 4 & 1 & 0 \\ -2 & 2 & 1 \end{pmatrix}, \qquad \mathbf{b} = \begin{pmatrix} 3 \\ 6 \\ 1 \end{pmatrix}.
$$

*Step 1, column 1.* Pivot $u_{11} = 2$. Multipliers: $\ell_{21} = 4/2 = 2$, $\ell_{31} = -2/2 = -1$.
Row 2 ← row 2 − $2$·row 1 $= (4, 1, 0) - (4, 2, 2) = (0, -1, -2)$.
Row 3 ← row 3 − $(-1)$·row 1 $= (-2, 2, 1) + (2, 1, 1) = (0, 3, 2)$. (Subtracting a negative multiple is adding. The multiplier stored is $-1$.)

$$
\begin{pmatrix} 2 & 1 & 1 \\ 0 & -1 & -2 \\ 0 & 3 & 2 \end{pmatrix}
$$

*Step 2, column 2.* Pivot $u_{22} = -1$. Multiplier $\ell_{32} = 3/(-1) = -3$.
Row 3 ← row 3 − $(-3)$·row 2 $= (0, 3, 2) + 3(0, -1, -2) = (0, 0, -4)$.

$$
U = \begin{pmatrix} 2 & 1 & 1 \\ 0 & -1 & -2 \\ 0 & 0 & -4 \end{pmatrix}, \qquad
L = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ -1 & -3 & 1 \end{pmatrix}.
$$

Pivots $2, -1, -4$. Three pivots: $A$ is invertible (Theorem 3.16).

*Check $A = LU$* by Theorem 3.10's recipe, row 3: $\ell_{31}\mathbf{u}_1^\top + \ell_{32}\mathbf{u}_2^\top + \mathbf{u}_3^\top = -(2,1,1) - 3(0,-1,-2) + (0,0,-4) = (-2, -1+3, -1+6-4) = (-2, 2, 1)$. ✓ Row 2: $2(2,1,1) + (0,-1,-2) = (4, 1, 0)$. ✓

*Forward substitution (3.4) on $\mathbf{b}$:* $c_1 = 3$; $c_2 = 6 - \ell_{21}c_1 = 6 - 6 = 0$; $c_3 = 1 - \ell_{31}c_1 - \ell_{32}c_2 = 1 + 3 + 0 = 4$. So $\mathbf{c} = (3, 0, 4)^\top$. (This is exactly what the two elimination steps do to $\mathbf{b}$.)

*Back substitution (3.3):* $-4x_3 = 4 \Rightarrow x_3 = -1$. $-x_2 - 2(-1) = 0 \Rightarrow x_2 = 2$. $2x_1 + 2 + (-1) = 3 \Rightarrow x_1 = 1$. $\mathbf{x} = (1, 2, -1)^\top$.

*Check:* $A\mathbf{x} = (2 + 2 - 1,\ 4 + 2 + 0,\ -2 + 4 - 1) = (3, 6, 1)$. ✓

*The forward product is not "minus $L$."* $E_{21} = E_{21}(2)$, $E_{31} = E_{31}(-1)$, $E_{32} = E_{32}(-3)$:

$$
E_{32}E_{31}E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ -5 & 3 & 1 \end{pmatrix}.
$$

The $(3,1)$ entry is $-\ell_{31} + \ell_{32}\ell_{21} = 1 + (-3)(2) = -5$, exactly the surviving cross term of §3.5. Its inverse is $L$, whose $(3,1)$ entry is $-1$. Never read $L$ off the forward product.

### Example 4.2: A singular $3\times 3$, and what the echelon form says

$$
A = \begin{pmatrix} 1 & 2 & 3 \\ 2 & 4 & 7 \\ 3 & 6 & 10 \end{pmatrix}.
$$

*Column 1.* Pivot $1$. $\ell_{21} = 2$: row 2 ← $(2,4,7) - 2(1,2,3) = (0, 0, 1)$. $\ell_{31} = 3$: row 3 ← $(3,6,10) - 3(1,2,3) = (0, 0, 1)$.

*Column 2.* Position $(2,2)$ is $0$, and so is $(3,2)$. No swap helps. **Column 2 is free.** Elimination breaks here, and this is the informative kind of breaking: column 2 is $2\times$ column 1.

*Column 3.* Move to column 3 with the pivot row still row 2. Pivot $1$ at $(2,3)$. $\ell_{32} = 1/1 = 1$: row 3 ← $(0,0,1) - (0,0,1) = (0,0,0)$.

$$
U = \begin{pmatrix} 1 & 2 & 3 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}.
$$

Two pivots (columns 1 and 3), one free column (2), one zero row. Rank $2$. By Theorem 3.16, $A$ is singular. By Theorem 3.8, $A\mathbf{x} = \mathbf{0}$ has a nonzero solution: set $x_2 = 1$, back-substitute: $x_3 = 0$ from row 2, $x_1 + 2 = 0$ from row 1. $\mathbf{s} = (-2, 1, 0)^\top$, and indeed $A\mathbf{s} = -2(\text{col } 1) + (\text{col } 2) = \mathbf{0}$.

*Right-hand side 1:* $\mathbf{b} = (1, 3, 4)^\top$. Apply the same operations: $c_2 = 3 - 2 = 1$, $c_3 = 4 - 3 = 1$, then $c_3 \leftarrow 1 - 1 = 0$. Augmented echelon form rows: $(1, 2, 3\,|\,1)$, $(0, 0, 1\,|\,1)$, $(0,0,0\,|\,0)$. Consistent. Theorem 3.7(c): $x_3 = 1$; $x_1 = 1 - 2x_2 - 3 = -2 - 2x_2$. With $x_2 = t$:

$$
\mathbf{x} = \begin{pmatrix} -2 \\ 0 \\ 1 \end{pmatrix} + t\begin{pmatrix} -2 \\ 1 \\ 0 \end{pmatrix}.
$$

Check the particular solution: $A(-2, 0, 1)^\top = (-2 + 3, -4 + 7, -6 + 10) = (1, 3, 4)$. ✓

*Right-hand side 2:* $\mathbf{b} = (1, 3, 5)^\top$. Now $c_3 = 5 - 3 = 2$, then $2 - 1 = 1$. Last row: $(0, 0, 0\,|\,1)$, i.e. $0 = 1$. **No solution.** In LA-1's language, $(1,3,5)^\top$ is not in the span of the columns, which is only a plane (two independent directions). The condition for consistency, read off the last row: $b_3 - 3b_1 - (b_2 - 2b_1) = 0$, i.e. $b_3 = b_1 + b_2$. Check: $4 = 1 + 3$ ✓, $5 \ne 1 + 3$ ✗.

### Example 4.3: Gauss–Jordan inverse of a $3\times 3$, and its $LDL^\top$

$$
A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 2 & 1 \\ 0 & 1 & 2 \end{pmatrix}.
$$

(This is a second-difference matrix; its relatives discretize $-d^2/dx^2$ in EM-3 and appear as the stiffness matrix of three masses on springs in CM-6.) Form $[A\,|\,I]$ and eliminate.

*Down.* $\ell_{21} = 1/2$: row 2 ← row 2 − $\tfrac12$ row 1 $= (0, \tfrac32, 1\,|\,-\tfrac12, 1, 0)$. $\ell_{31} = 0$: nothing. $\ell_{32} = 1/(3/2) = \tfrac23$: row 3 ← row 3 − $\tfrac23$ row 2 $= (0, 0, 2 - \tfrac23\,|\,0 + \tfrac13, -\tfrac23, 1) = (0, 0, \tfrac43\,|\,\tfrac13, -\tfrac23, 1)$.

$$
\left[\begin{array}{ccc|ccc} 2 & 1 & 0 & 1 & 0 & 0 \\ 0 & \tfrac32 & 1 & -\tfrac12 & 1 & 0 \\ 0 & 0 & \tfrac43 & \tfrac13 & -\tfrac23 & 1 \end{array}\right]
$$

Pivots $2, \tfrac32, \tfrac43$. Product $4$ (this will be the determinant in LA-6). Three pivots: invertible.

*Up.* Clear above pivot 3: row 2 ← row 2 − $\tfrac{1}{4/3}$ row 3 $=$ row 2 − $\tfrac34$ row 3 $= (0, \tfrac32, 0\,|\,-\tfrac12 - \tfrac14, 1 + \tfrac12, -\tfrac34) = (0, \tfrac32, 0\,|\,-\tfrac34, \tfrac32, -\tfrac34)$. Row 1 already has $0$ in column 3. Clear above pivot 2: row 1 ← row 1 − $\tfrac{1}{3/2}$ row 2 $=$ row 1 − $\tfrac23$ row 2 $= (2, 0, 0\,|\,1 + \tfrac12, -1, \tfrac12) = (2, 0, 0\,|\,\tfrac32, -1, \tfrac12)$.

*Scale.* Row 1 $\div 2$, row 2 $\div \tfrac32$, row 3 $\div \tfrac43$:

$$
\left[\begin{array}{ccc|ccc} 1 & 0 & 0 & \tfrac34 & -\tfrac12 & \tfrac14 \\ 0 & 1 & 0 & -\tfrac12 & 1 & -\tfrac12 \\ 0 & 0 & 1 & \tfrac14 & -\tfrac12 & \tfrac34 \end{array}\right],
\qquad
A^{-1} = \frac14\begin{pmatrix} 3 & -2 & 1 \\ -2 & 4 & -2 \\ 1 & -2 & 3 \end{pmatrix}.
$$

*Check:* row 2 of $A$ times column 1 of $A^{-1}$: $\tfrac14(1\cdot 3 + 2\cdot(-2) + 1\cdot 1) = 0$ ✓; row 2 times column 2: $\tfrac14(-2 + 8 - 2) = 1$ ✓. $A^{-1}$ is symmetric, as Theorem 3.19(b) predicts for a symmetric $A$. And it is **dense**: $A$ has four zeros, $A^{-1}$ has none. That is §3.7's structural argument in miniature.

*$LDL^\top$.* From the downward pass, $L = \begin{pmatrix} 1 & 0 & 0 \\ \tfrac12 & 1 & 0 \\ 0 & \tfrac23 & 1 \end{pmatrix}$ and $D = \operatorname{diag}(2, \tfrac32, \tfrac43)$. Theorem 3.22 says $U' = L^\top$; check one entry: $U = DU'$ should have $u_{23} = \tfrac32\cdot\tfrac23 = 1$, and the downward pass gave $u_{23} = 1$. ✓ All pivots positive: $A$ is positive definite (LA-8), and these three numbers are the squared lengths of the "new directions" in a sense LA-5 makes precise.

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

Kirchhoff's current law at each node: current in equals current out, with current through a resistor $= (\text{voltage drop})/R = g\cdot(\text{drop})$ where $g = 1/R$ is the conductance. Conductances $g_1 = 1, g_2 = \tfrac12, g_3 = \tfrac12, g_4 = 1$ (siemens).

Node 1: $g_1 v_1 + g_2(v_1 - v_2) = 3$. Node 2: $g_2(v_2 - v_1) + g_3(v_2 - v_3) = 0$. Node 3: $g_3(v_3 - v_2) + g_4 v_3 = 0$. In matrix form $G\mathbf{v} = \mathbf{i}$:

$$
G = \begin{pmatrix} \tfrac32 & -\tfrac12 & 0 \\ -\tfrac12 & 1 & -\tfrac12 \\ 0 & -\tfrac12 & \tfrac32 \end{pmatrix}, \qquad \mathbf{i} = \begin{pmatrix} 3 \\ 0 \\ 0 \end{pmatrix}.
$$

$G$ is symmetric (the conductance between nodes 1 and 2 is the same seen from either end), the diagonal entry $G_{kk}$ is the total conductance touching node $k$, and the off-diagonal $G_{jk}$ is minus the conductance between $j$ and $k$. Every conductance matrix has this form.

*Elimination.* Pivot 1: $u_{11} = \tfrac32$. $\ell_{21} = (-\tfrac12)/\tfrac32 = -\tfrac13$, $\ell_{31} = 0$. Row 2 ← row 2 + $\tfrac13$ row 1 $= (0, 1 - \tfrac16, -\tfrac12) = (0, \tfrac56, -\tfrac12)$. Pivot 2: $u_{22} = \tfrac56$. $\ell_{32} = (-\tfrac12)/\tfrac56 = -\tfrac35$. Row 3 ← row 3 + $\tfrac35$ row 2 $= (0, 0, \tfrac32 - \tfrac35\cdot\tfrac12) = (0, 0, \tfrac32 - \tfrac{3}{10}) = (0, 0, \tfrac65)$.

$$
L = \begin{pmatrix} 1 & 0 & 0 \\ -\tfrac13 & 1 & 0 \\ 0 & -\tfrac35 & 1 \end{pmatrix}, \qquad
U = \begin{pmatrix} \tfrac32 & -\tfrac12 & 0 \\ 0 & \tfrac56 & -\tfrac12 \\ 0 & 0 & \tfrac65 \end{pmatrix}.
$$

*Solve.* Forward: $c_1 = 3$, $c_2 = 0 + \tfrac13\cdot 3 = 1$, $c_3 = 0 + \tfrac35\cdot 1 = \tfrac35$. Back: $v_3 = \tfrac35/\tfrac65 = \tfrac12$; $\tfrac56 v_2 = 1 + \tfrac12\cdot\tfrac12 = \tfrac54 \Rightarrow v_2 = \tfrac32$; $\tfrac32 v_1 = 3 + \tfrac12\cdot\tfrac32 = \tfrac{15}{4} \Rightarrow v_1 = \tfrac52$.

$\mathbf{v} = (2.5, 1.5, 0.5)^\top$ volts. *Physical check:* $2.5\,\mathrm{A}$ flows through $R_1$ to ground; the remaining $0.5\,\mathrm{A}$ crosses $R_2$ (drop $1\,\mathrm{V}$, so $v_2 = 1.5$ ✓), crosses $R_3$ (drop $1\,\mathrm{V}$, $v_3 = 0.5$ ✓), and exits through $R_4$ ($0.5\,\mathrm{V}/1\,\Omega = 0.5\,\mathrm{A}$ ✓). Total $3\,\mathrm{A}$. ✓

*What the pivots mean.* Pivot 1, $\tfrac32 = g_1 + g_2$, is the conductance seen at node 1 if nodes 2 and 3 were held at zero. Pivot 2 is more interesting: $\tfrac56 = \tfrac13 + \tfrac12$. The $\tfrac13$ is the conductance of $R_1$ and $R_2$ *in series* ($1 + 2 = 3\,\Omega$), the $\tfrac12$ is $g_3$. So pivot 2 is the conductance seen at node 2 *after node 1 has been eliminated from the circuit* and its resistors merged into an equivalent. Pivot 3, $\tfrac65 = \tfrac15 + 1$: $R_3 + R_2 + R_1 = 5\,\Omega$ in series gives $\tfrac15$, plus $g_4 = 1$. Each elimination step removes one node from the circuit and replaces its connections by their series–parallel equivalent. Electrical engineers call this **Kron reduction**; in structural mechanics the same operation on a stiffness matrix is **static condensation**; in linear algebra it is the **Schur complement** (exercise C2). And a **zero pivot would mean a node with no path to ground**: disconnect $R_1$ and $R_4$ and the conductance matrix becomes $\begin{pmatrix} \tfrac12 & -\tfrac12 & 0 \\ -\tfrac12 & 1 & -\tfrac12 \\ 0 & -\tfrac12 & \tfrac12\end{pmatrix}$, elimination gives pivots $\tfrac12, \tfrac12, 0$, and the null vector is $(1, 1, 1)^\top$: with no reference to ground, all voltages can shift by a common constant. Singular matrix, physically indeterminate potential. Same fact, two languages.

### Example 4.5: The $10^{-20}$ example, by hand in floating point

Already done in §3.6; here is the arithmetic in the form the project reproduces. Without pivoting: $\ell_{21} = 10^{20}$, $U = \begin{pmatrix} 10^{-20} & 1 \\ 0 & \mathrm{fl}(1 - 10^{20})\end{pmatrix} = \begin{pmatrix} 10^{-20} & 1 \\ 0 & -10^{20}\end{pmatrix}$, and $LU = \begin{pmatrix} 10^{-20} & 1 \\ 1 & 10^{20}\cdot 1 + (-10^{20}) \end{pmatrix} = \begin{pmatrix} 10^{-20} & 1 \\ 1 & 0\end{pmatrix} \ne A$. The computed solution of $A\mathbf{x} = (1, 2)^\top$ is $(0, 1)^\top$; the true one is $(1, 1)^\top$ to 20 digits. With pivoting, the computed solution is $(1, 1)^\top$. The project sweeps $\varepsilon$ from $10^{-1}$ to $10^{-20}$ and plots the error of both methods; you will see the unpivoted error jump to $O(1)$ once $\varepsilon$ drops below about $10^{-16}$, the machine precision, and stay there.

## 5. Why it matters and where it connects

- **LA-3** takes the special solutions $\mathbf{s}_f$ of Theorem 3.7(c) and calls their span the null space $N(A)$; the pivot columns of $A$ become a basis of the column space $C(A)$; the nonzero rows of the RREF become a basis of the row space. The provisional rank of Definition 3.20 becomes the common dimension of $C(A)$ and the row space, and $n - r$ (the number of free columns) becomes $\dim N(A)$. The uniqueness of the RREF is proved there.
- **LA-4**: Rank–Nullity, $r + \dim N(A) = n$, is Theorem 3.7(c) restated without an algorithm.
- **LA-5**: least squares solves $A^\top A\hat{\mathbf{x}} = A^\top\mathbf{b}$; $A^\top A$ is symmetric (§3.10) with positive pivots when the columns of $A$ are independent, so Cholesky ($LDL^\top$ with $D > 0$) solves it. LA-5 also introduces $QR$, which replaces $LU$ when the matrix is rectangular.
- **LA-6**: $\det A = \pm(\text{product of pivots})$, sign from the number of row swaps. Example 4.3's pivots multiply to $4$; the project checks `np.linalg.det`.
- **LA-8**: a symmetric matrix is positive definite iff all its pivots are positive iff $A = LDL^\top$ with $D > 0$. Example 4.3's matrix qualifies. The energy $\mathbf{x}^\top A\mathbf{x} = \sum_i d_i y_i^2$ with $\mathbf{y} = L^\top\mathbf{x}$ (exercise B7) is the first step of that proof.
- **LA-10** (Cholesky, and when to use which factorization) and **LA-13** (backward stability of $LU$, the growth factor, conditioning; Example 4.5 done properly) both start from (3.6).
- **CM-6**: the equations of motion for small oscillations are $M\ddot{\mathbf{q}} + K\mathbf{q} = \mathbf{0}$ with $K$ a symmetric stiffness matrix of the same pattern as Example 4.3 and Example 4.4's $G$. Static condensation (eliminating interior degrees of freedom) is literally §3.5's elimination; the pivots are effective stiffnesses, as in Example 4.4 they are effective conductances.
- **EM-3**: the finite-difference Laplacian for a boundary-value problem is a banded symmetric matrix; $LU$ of a banded matrix stays banded, costing $O(n)$ for tridiagonal. Forming its inverse would cost $O(n^2)$ storage and destroy the band (Example 4.3).
- **ST-3**: the normal equations $(X^\top X)\boldsymbol\beta = X^\top\mathbf{y}$ are solved by Cholesky or $QR$; `inv(X.T @ X)` in a regression code is the mistake §3.7 warns about.
- **ML-2**: Newton's method solves $H\mathbf{d} = -\nabla f$ at every iteration with $H$ the Hessian; one $LDL^\top$ per step, never $H^{-1}$.
- **QC-2**: the classical reversible gates (SWAP, CNOT, Toffoli) are permutation matrices on the computational basis, and Theorem 3.23, $P^{-1} = P^\top$, is why they are unitary with no further argument.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are drills: a $2\times 2$ elimination with $E$, $L$, $U$; writing elementary matrices and their inverses; an RREF by hand; substitution with a flop count; a $2\times 2$ Gauss–Jordan; a permutation matrix. B1–B9 are the module: B2 (full $LU$ with the forward product contrast), B3 (prove Corollary 3.9 and $PA = LU$ by hand), B5 (algebra of inverses including the converse of Theorem 3.19(a)), B7 ($LDL^\top$ and the energy $\mathbf{x}^\top A\mathbf{x}$), and B9 (implement $LDL^\top$ in `project.py`) are essential. C1 (Wilkinson's growth matrix), C2 (block $LU$ and the Schur complement, with the circuit), C3 (pivot count is independent of the row operations), C4 (Sherman–Morrison rank-one update) go one step beyond.

## 7. Coding project

`project.py` builds the algorithms of this module from scratch and checks every number in §4.

1. **Elementary matrices.** `E(i, j, l)`, `P(i, j)`, `D(i, c)` from (3.1). Asserts: left multiplication acts on rows and right multiplication on columns (Theorem 3.1); each inverse is as derived in §3.1; for Example 4.1, $E_{32}E_{31}E_{21}$ has $(3,1)$ entry $-5$ while $L$ has $-1$.
2. **$LU$ with partial pivoting.** `lu_decompose(A, pivot=True)` returning `P, L, U` with $PA = LU$, including the swap-the-earlier-columns-of-$L$ step of Theorem 3.13; `forward_sub` (3.4), `back_sub` (3.3), `solve(A, b)`. Asserted against `scipy.linalg.lu` (translating its $A = PLU$ convention) and `np.linalg.solve` on random matrices of several sizes. Examples 4.1, 4.3, 4.4 reproduced to the digit, including the pivots-as-conductances check.
3. **RREF with a tolerance, and classification.** `rref(A, tol)` returns the reduced form and the pivot columns; `classify(A, b)` returns `'unique'`, `'none'`, or `'infinite'` plus the free columns, a particular solution, and the special solutions (Theorem 3.7). Example 4.2 with both right-hand sides; LA-1's three $2\times 2$ outcomes.
4. **Gauss–Jordan inverse.** `gauss_jordan_inverse(A)` from (3.8); Example 4.3; Theorem 3.19 checked on random matrices; the never-invert demonstration: timing of `inv` vs. `solve`, residuals of both routes, and the density of the inverse of a tridiagonal matrix.
5. **The $10^{-20}$ example.** Elimination with and without pivoting for $\varepsilon = 10^{-1}, \dots, 10^{-20}$; figure `pivoting_error.png`.
6. **Costs.** A flop counter in `lu_decompose` verified against the exact sum in §3.7 and plotted against $\tfrac23 n^3$ (`flop_count.png`); wall-clock timing of "$LU$ once then $k$ substitutions" versus "$k$ solves from scratch" for growing $n$ (`timing_lu_vs_solve.png`).

**Build it yourself first.** Before reading Part 2, write `lu_decompose` without pivoting (a loop over columns, `np.outer` for the rank-one update), then `forward_sub` and `back_sub` directly from (3.3) and (3.4). Make them pass on Example 4.1. Then add pivoting and discover for yourself which rows of $L$ need swapping. Expected result: `np.allclose(P @ A, L @ U)` on random $50\times 50$ matrices, and your `solve` agreeing with NumPy to about $10^{-13}$.

**Your turn** (at the end of the script): implement `ldlt(A)` for symmetric $A$ (Theorem 3.22) and check $LDL^\top = A$ and that $D$ holds the pivots; estimate rank by counting pivots above a tolerance, and watch what happens to a matrix of exact rank 2 when the tolerance is too small; construct a random singular matrix as $BC$ with $B$ being $n\times r$ and see what `solve` does with it.

## 8. Where you'll struggle

- **Signs of the multipliers.** The operation is row $i$ ← row $i$ − $\ell_{ij}$ · row $j$, with $\ell_{ij} = a_{ij}/a_{jj}$ computed from the *current* entries. $L$ stores $+\ell_{ij}$. $E_{ij}$ stores $-\ell_{ij}$. When $\ell_{ij}$ is negative you *add* a multiple of row $j$, and $L$ still stores the negative number (Example 4.1: $\ell_{31} = -1$, you added row 1, $L_{31} = -1$). People get this wrong by (i) storing the negative, (ii) computing $\ell$ from the original $A$ instead of the current matrix at step 2, (iii) reading $L$ off $E_{32}E_{31}E_{21}$, which is wrong by a cross term. The fix: after every hand elimination, check one row of $A = LU$ with Theorem 3.10's recipe (row $i$ of $A = \sum_j \ell_{ij}\mathbf{u}_j^\top$). It takes thirty seconds and catches all three errors.
- **Row operations versus column operations.** Left multiplication $EA$ acts on rows; right multiplication $AE$ acts on columns (Theorem 3.1). Elimination is all left multiplication, which is why $[A\,|\,I] \to [I\,|\,A^{-1}]$ works: the same $M$ hits both blocks from the left. If you ever find yourself doing column operations on an augmented matrix, you have changed the unknowns, not the equations, and the answer is for a different problem. Say "rows" out loud.
- **"Zero pivot" panic.** A zero in the pivot position means *look below it*. If something below is nonzero, swap and continue; nothing is wrong. Only if the entire rest of the column is zero is the column free. Students either give up at the first zero (declaring a singular matrix that isn't) or swap when they shouldn't (there's nothing to swap in). Example 4.2 is the second case.
- **The inverse reflex.** `x = np.linalg.inv(A) @ b` is in a lot of code you'll read, including code by people who should know better. It costs three times what `solve` costs, saves nothing per right-hand side, destroys sparsity, and loses the backward-stability guarantee (§3.7). It also hides singular and nearly singular matrices behind a warning you'll ignore. The reflex to fight is the algebraic one: $\mathbf{x} = A^{-1}\mathbf{b}$ is a correct *statement* and a bad *instruction*. Write $A^{-1}$ in proofs; call `solve` in code.
- **"Rank" is provisional here.** Definition 3.20 counts pivots, which is an algorithm's output, and we have not proved that different legal choices give the same count (C3 does; LA-3 does it properly). Worse, in floating point the count depends on a tolerance: a matrix of exact rank 2 produces a third pivot around $10^{-16}$, not $0$, and `rank` without a `tol` will say 3. Numerical rank is a genuinely subtle notion that needs the SVD (LA-9). Don't build arguments on "the rank" until LA-3 gives it a definition without an algorithm in it.
- **The pivoting bug.** When you add row swaps to your `lu_decompose`, the swap must be applied to the already-filled columns of $L$ too (Theorem 3.13, step 4). Everyone forgets it once. The symptom: $PA \ne LU$ while each individual step looked right. The project asserts against SciPy so you can't ship it.
- **Theorem 3.7(c) and the direction of "free."** Free variables are chosen *first*, then pivot variables are solved *from the bottom up*. Students try to solve for free variables, or work top-down. Draw the staircase, mark the pivot columns, and read the pivot rows from the bottom.
- **Thinking every invertible matrix has an $LU$.** $\begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix}$ does not (§3.6). Every invertible matrix has a $PA = LU$. The permutation is part of the factorization, not an afterthought.
- **Exercise B3(b)** ($PA = LU$ by hand with two swaps) is harder than it looks because the second swap reorders column 1 of $L$. Do it slowly, writing the permutation vector after each swap.

## 9. Check yourself

1. Write the $3\times 3$ elementary matrix that subtracts $5$ times row 1 from row 3. Write its inverse. Which of the two has $+5$ in it?
2. Why does left-multiplying by $E$ change rows and not columns? One sentence using the row-times-column rule.
3. State the three outcomes for $A\mathbf{x} = \mathbf{b}$ in terms of pivot columns and the augmented column.
4. Prove in four lines: an $n\times n$ matrix with fewer than $n$ pivots has a nonzero null vector.
5. Why is $L$ exactly the matrix of multipliers, with no further arithmetic? Give the "run elimination backwards" argument for row 3.
6. Why is $E_{32}E_{31}E_{21}$ *not* simply $I$ with $-\ell_{ij}$ below the diagonal? Which entry is wrong and by how much?
7. Give the $2\times 2$ matrix for which elimination without pivoting fails in floating point, and say in one sentence what goes wrong.
8. About how many flops for $LU$ of a $2000\times 2000$ matrix? For one back substitution? For $A^{-1}$?
9. If $BA = I$ and $AC = I$, prove $B = C$.
10. State and prove $(AB)^{-1} = B^{-1}A^{-1}$.
11. Why does Gauss–Jordan on $[A\,|\,I]$ produce $A^{-1}$ in the right block?
12. For a symmetric $A$ with an $LDU$ factorization, why is $U = L^\top$? What theorem does the proof lean on?
13. What is $P^{-1}$ for a permutation matrix, and why?
14. In the circuit of Example 4.4, what does a pivot measure? What does a zero pivot mean physically?

## 10. Pair with

- **Strang, *Introduction to Linear Algebra*, ch. 2** (§2.2–2.7: elimination, elimination using matrices, rules for matrix operations, inverses, $A = LU$, transposes and permutations). Strang's treatment of "$L$ is the record of the multipliers" is the clearest in print; his §2.6 is the source of the $E_{32}E_{31}E_{21}$-versus-$L$ comparison.
- **MIT 18.06, lectures 2–5 (Strang), on OCW.** Lecture 2 is elimination with matrices; lecture 3 is multiplication and inverses with Gauss–Jordan done live; lecture 4 is $A = LU$ and the operation count; lecture 5 is transposes, permutations, and the start of vector spaces. Watch lecture 4 after §3.5 and notice how he counts the $n^3/3$ multiplications.
- **Trefethen & Bau, *Numerical Linear Algebra*, lectures 20–21.** Lecture 20 states Gaussian elimination as $L^{-1}A = U$ and gives the flop count in their style; lecture 21 is pivoting, with the $PA = LU$ derivation of Theorem 3.13 and the worst-case growth factor. Lecture 22 (stability) is for LA-13, but the first page is readable now.
- **Grcar, "How ordinary elimination became Gaussian elimination," *Historia Mathematica* 38 (2011).** Optional. The real history of §1, with the primary sources.
