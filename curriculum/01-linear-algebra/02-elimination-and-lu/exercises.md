# LA-2 Exercises

Conventions as in the notes: row $i$ ← row $i$ − $\ell_{ij}$ · row $j$ with $\ell_{ij}$ the multiplier; $E_{ij}(\ell) = I - \ell\,\mathbf{e}_i\mathbf{e}_j^\top$; $L$ stores $+\ell_{ij}$; $PA = LU$; a flop is one floating-point operation.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Solve $\begin{pmatrix} 3 & 1 \\ 6 & 5 \end{pmatrix}\mathbf{x} = \begin{pmatrix} 4 \\ 11 \end{pmatrix}$ by elimination. Write down the multiplier $\ell_{21}$, the elementary matrix $E_{21}$, and $L$ and $U$. Verify $LU = A$ by multiplication. Then find $\mathbf{x}$ by forward substitution (3.4) followed by back substitution (3.3), and check it in the original equations.

**A2.** For $n = 3$, write the elementary matrix $E$ that subtracts $4$ times row 1 from row 3, the permutation $P$ that swaps rows 1 and 2, and the diagonal matrix $D$ that multiplies row 2 by $\tfrac15$. Write the inverse of each. For $A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix}$ compute $EA$ and $AE$. Describe in words what $AE$ did to $A$ (Theorem 3.1).

**A3.** Reduce $M = \begin{pmatrix} 1 & 2 & 1 & 3 \\ 2 & 4 & 0 & 2 \\ 1 & 2 & -1 & -1 \end{pmatrix}$ to row echelon form and then to RREF. List the pivot columns, the free columns, and the rank. Write every solution of $M\mathbf{x} = \mathbf{0}$ in the form $s\,\mathbf{s}_1 + t\,\mathbf{s}_2$ (Theorem 3.7(c)).

**A4.** Using $L$ and $U$ from Example 4.1, solve $A\mathbf{x} = \mathbf{b}$ for $\mathbf{b} = (2, 7, -5)^\top$ by forward then back substitution. Check the answer against $A$ directly. Then count the flops you actually performed in each substitution and compare with the counts $n(n-1)$ and $n^2$ from §3.7.

**A5.** Compute the inverse of $A = \begin{pmatrix} 2 & 3 \\ 1 & 4 \end{pmatrix}$ by Gauss–Jordan on $[A\,|\,I]$. Check it against the $2\times 2$ formula at the end of §3.8. Then compute $(A^\top)^{-1}$ by the formula and confirm it equals $(A^{-1})^\top$ (Theorem 3.19(b)).

**A6.** Let $P = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}$. For a matrix $A$ with rows $\mathbf{r}_1^\top, \mathbf{r}_2^\top, \mathbf{r}_3^\top$, what are the rows of $PA$? Compute $P^\top$, $P^2$, and $P^3$. Confirm $P^{-1} = P^\top$ and identify which power of $P$ it equals.

## B. Core (do all)

**B1.** (Elementary matrices, carefully.) Using Theorem 3.1, prove that $E_{ij}(\ell)A$, $P_{ij}A$, and $D_i(c)A$ perform operations (R1), (R2), (R3) on $A$. Then prove $P_{ij}^2 = I$ by *direct multiplication* of the formula in (3.1), using the product rule (3.2) for the matrices $\mathbf{e}_a\mathbf{e}_b^\top$ (not the "swap twice" argument). Hint: write $P_{ij} = I + Q$ and show $Q^2 = -2Q$.

**B2.** (Full $LU$.) Let $A = \begin{pmatrix} 1 & 2 & 0 \\ 3 & 4 & 1 \\ 2 & 5 & -1 \end{pmatrix}$.
(a) Eliminate, tracking every multiplier. Write $L$ and $U$ and check one row of $A = LU$ with the recipe of Theorem 3.10.
(b) Solve $A\mathbf{x} = (-1, 1, -5)^\top$ by the two substitutions.
(c) Compute the forward product $E_{32}E_{31}E_{21}$ and compare its $(3,1)$ entry with $L_{31}$. Express the difference in terms of the multipliers.
(d) Write the $LDU'$ factorization (3.9).

**B3.** (Free columns and row swaps.)
(a) Prove Corollary 3.9: if $A$ is $m\times n$ with $n > m$, then $A\mathbf{x} = \mathbf{0}$ has a nonzero solution. Write out every step; do not cite Theorem 3.8 as a black box.
(b) Compute $PA = LU$ for $A = \begin{pmatrix} 0 & 2 & 1 \\ 1 & 1 & 0 \\ 2 & 3 & 1 \end{pmatrix}$ using partial pivoting (largest entry in the column). Keep a permutation vector, and when you swap rows at step 2, swap the corresponding entries of column 1 of $L$ (Theorem 3.13). Write $P$, $L$, $U$ and verify $PA = LU$ by multiplication.

**B4.** (Reading a singular system.) Let $A = \begin{pmatrix} 1 & 1 & 1 \\ 1 & 2 & 3 \\ 2 & 3 & 4 \end{pmatrix}$.
(a) Find an echelon form, the rank, and a nonzero solution of $A\mathbf{x} = \mathbf{0}$.
(b) Carry a general $\mathbf{b} = (b_1, b_2, b_3)^\top$ through the same row operations and find the single linear condition on $\mathbf{b}$ for $A\mathbf{x} = \mathbf{b}$ to be solvable. Interpret it in LA-1's language.
(c) For $\mathbf{b} = (1, 2, 3)^\top$, write the complete solution as $\mathbf{x}_p + t\,\mathbf{s}$.

**B5.** (Algebra of inverses.) Let $A, B$ be $n\times n$.
(a) If $A$ is invertible, prove $A^k$ is invertible with $(A^k)^{-1} = (A^{-1})^k$ for every $k \ge 1$.
(b) Prove the converse of Theorem 3.19(a): if $AB$ is invertible then $A$ and $B$ are both invertible. (Hint: Corollary 3.18 for $B$; then write $A$ as a product of invertible matrices.)
(c) Prove that if $A$ is invertible and $AB = AC$ then $B = C$. Give $2\times 2$ matrices $A, B, C$ with $AB = AC$ but $B \ne C$.

**B6.** (Gauss–Jordan, $3\times 3$.) Let $A = \begin{pmatrix} 1 & 0 & 2 \\ 2 & 1 & 3 \\ 0 & 1 & 1 \end{pmatrix}$.
(a) Compute $A^{-1}$ by row-reducing $[A\,|\,I]$. Show the matrix after the downward pass, after the upward pass, and after scaling.
(b) Check $AA^{-1} = I$.
(c) From the downward pass, read off $L$ and $U$ and the pivots. What is their product? (LA-6 will call it $\det A$.)
(d) In two sentences, explain why the right block ends as $A^{-1}$ (Corollary 3.17).

**B7.** ($LDL^\top$ and energy.) Suppose $A = LDL^\top$ with $L$ unit lower triangular and $D = \operatorname{diag}(d_1, \dots, d_n)$.
(a) Prove: if every $d_i > 0$ then $\mathbf{x}^\top A\mathbf{x} > 0$ for every $\mathbf{x} \ne \mathbf{0}$. (Hint: substitute $\mathbf{y} = L^\top\mathbf{x}$ and use Lemma 3.11.)
(b) For the matrix of Example 4.3, write $\mathbf{x}^\top A\mathbf{x}$ explicitly as a weighted sum of three squares, and expand to confirm it equals $2x_1^2 + 2x_1x_2 + 2x_2^2 + 2x_2x_3 + 2x_3^2$.
(c) Prove the converse of (a): if some $d_k \le 0$, exhibit $\mathbf{x} \ne \mathbf{0}$ with $\mathbf{x}^\top A\mathbf{x} \le 0$.

**B8.** (Operation counts.)
(a) Derive the exact flop count for $LU$ without pivoting, $\sum_{m=1}^{n-1}(m + 2m^2) = \tfrac{n(n-1)}{2} + \tfrac{n(n-1)(2n-1)}{3}$, from the description of Algorithm 3.5. Evaluate it for $n = 3$ and confirm by counting the arithmetic in Example 4.1 by hand.
(b) With one $LU$ and $k$ right-hand sides the cost is $\tfrac23 n^3 + 2kn^2$; from scratch it is $k\cdot\tfrac23 n^3$. For $n = 500$, find the $k$ at which the from-scratch cost is exactly $100$ times the factor-once cost.
(c) Computing $A^{-1}$ means solving $A\mathbf{x} = \mathbf{e}_j$ for $j = 1, \dots, n$. Show that forward substitution with right-hand side $\mathbf{e}_j$ costs only about $(n-j)^2$ flops (why?), so the $n$ forward substitutions together cost about $n^3/3$, and the total for $A^{-1}$ via $LU$ is about $2n^3$. Compare with the $\tfrac23 n^3$ of a single factorization.

**B9.** (Code.) In `project.py`, the "Your turn" section asks for `ldlt(A)`. Implement it: for symmetric $A$ with no row swaps needed, return `L, d` with $A = L\operatorname{diag}(d)L^\top$ (Theorem 3.22), computing only the lower half. Assert (i) `L @ np.diag(d) @ L.T` reproduces random symmetric positive definite matrices $BB^\top + I$ to `1e-10`; (ii) for Example 4.3, `d == [2, 3/2, 4/3]` and `L[1,0] == 1/2`, `L[2,1] == 2/3`; (iii) `d` equals the diagonal of $U$ from `lu_decompose(A, pivot=False)`. Then make sure `lu_decompose` raises `np.linalg.LinAlgError` when a pivot has absolute value at most a tolerance `tol` (the script's version does; read how). Confirm that it raises on the matrix of Example 4.2 with `tol=1e-12`, and that on a random matrix of exact rank 2 (built as $BC$ with $B$ being $6\times 2$ and $C$ being $2\times 6$) it raises with `tol=1e-10` but returns a factorization with `tol=0`. Write a two-sentence comment on what that factorization's third pivot looks like and what `solve` then returns.

## C. Challenge (do at least one)

**C1.** (Wilkinson's growth matrix.) Let $W_n$ have $1$ on the diagonal, $-1$ everywhere below the diagonal, $1$ in the entire last column, and $0$ elsewhere; for example $W_4 = \begin{pmatrix} 1 & 0 & 0 & 1 \\ -1 & 1 & 0 & 1 \\ -1 & -1 & 1 & 1 \\ -1 & -1 & -1 & 1 \end{pmatrix}$.
(a) Eliminate $W_4$ with partial pivoting. Show that no row swap is ever needed (every candidate pivot has the same absolute value as the diagonal; take the diagonal) and find $U$.
(b) Prove by induction on $n$ that the last column of $U$ is $(1, 2, 4, \dots, 2^{n-1})^\top$. Conclude that the growth factor $\max_{ij}|u_{ij}| / \max_{ij}|a_{ij}|$ equals $2^{n-1}$, the worst case allowed by partial pivoting.
(c) Code: build $W_{60}$, set $\mathbf{b} = W_{60}\mathbf{1}$ so the exact solution is $\mathbf{1}$, and solve with `scipy.linalg.lu_factor`/`lu_solve` and with your `solve`. Report $\max_i|x_i - 1|$ and the relative residual $\|W\mathbf{x} - \mathbf{b}\|/\|\mathbf{b}\|$. Explain what you see using $2^{59} \approx 5.8\times 10^{17}$ and machine precision $\approx 1.1\times 10^{-16}$.

**C2.** (Block elimination and the Schur complement.) Partition a square matrix as $A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}$ with $A_{11}$ square and invertible, and define $S = A_{22} - A_{21}A_{11}^{-1}A_{12}$.
(a) Verify by block multiplication that $A = \begin{pmatrix} I & 0 \\ A_{21}A_{11}^{-1} & I \end{pmatrix}\begin{pmatrix} A_{11} & A_{12} \\ 0 & S \end{pmatrix}$. (This is $LU$ with blocks; when $A_{11}$ is $1\times 1$ it is one step of Algorithm 3.5.)
(b) Prove that $A$ is invertible if and only if $S$ is invertible.
(c) For the conductance matrix $G$ of Example 4.4, take $A_{11} = G_{11}$ ($1\times 1$) and compute $S$. Show its $(1,1)$ entry is $\tfrac56$ and interpret it as a series–parallel combination. Then take $A_{11}$ to be the top-left $2\times 2$ block and show $S = \tfrac65$.

**C3.** (The pivot count does not depend on the row operations.) Let $U$ be any echelon form of $A$ (any legal sequence of row operations, any choice of swaps).
(a) Prove: column $j$ of $U$ is a pivot column if and only if it is *not* a linear combination of columns $1, \dots, j-1$ of $U$. (For the "free" direction, show the coefficients can be found by back substitution on the pivot rows above.)
(b) Prove: for invertible $E$, column $j$ of $A$ is a combination of columns $1, \dots, j-1$ of $A$ if and only if the same holds for $EA$ with the same coefficients. (Hint: a column relation is a statement $A\mathbf{y} = \mathbf{0}$ for a particular $\mathbf{y}$; Theorem 3.3.)
(c) Conclude that the set of pivot columns, and hence the rank of Definition 3.20, is the same for every echelon form of $A$.

**C4.** (Sherman–Morrison: updating a solve without refactoring.) Let $A$ be invertible and $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$ with $\sigma := 1 + \mathbf{v}^\top A^{-1}\mathbf{u} \ne 0$.
(a) Verify by multiplication that $(A + \mathbf{u}\mathbf{v}^\top)^{-1} = A^{-1} - \dfrac{A^{-1}\mathbf{u}\mathbf{v}^\top A^{-1}}{\sigma}$.
(b) Given $L, U$ with $A = LU$, show that $(A + \mathbf{u}\mathbf{v}^\top)\mathbf{x} = \mathbf{b}$ can be solved in $O(n^2)$ flops using only two pairs of triangular substitutions and $O(n)$ extra work. Write the algorithm.
(c) Show that if $\sigma = 0$ then $A + \mathbf{u}\mathbf{v}^\top$ is singular, by exhibiting a nonzero null vector.
(d) Code: check (a) numerically on a random $6\times 6$ $A$ and random $\mathbf{u}, \mathbf{v}$, and time (b) against a fresh `solve` for $n = 400$.
