# LA-2 Solutions

Conventions as in the notes: row $i$ ← row $i$ − $\ell_{ij}$ · row $j$; $L$ stores $+\ell_{ij}$; $E_{ij}(\ell) = I - \ell\,\mathbf{e}_i\mathbf{e}_j^\top$; $PA = LU$. Every number below is reproduced in `project.py`, Part 7.

## A. Warm-up

**A1.** $\ell_{21} = 6/3 = 2$. Row 2 ← row 2 − $2$·row 1 $= (6, 5) - (6, 2) = (0, 3)$.

$$E_{21} = E_{21}(2) = \begin{pmatrix} 1 & 0 \\ -2 & 1 \end{pmatrix}, \qquad U = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}, \qquad L = E_{21}^{-1} = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}.$$

$LU$: row 1 $= 1\cdot(3, 1) = (3, 1)$; row 2 $= 2\cdot(3, 1) + (0, 3) = (6, 5)$. $= A$. ✓

Forward (3.4): $c_1 = 4$, $c_2 = 11 - 2\cdot 4 = 3$. Back (3.3): $x_2 = 3/3 = 1$, $x_1 = (4 - 1\cdot 1)/3 = 1$. $\mathbf{x} = (1, 1)^\top$. Check: $3 + 1 = 4$ ✓, $6 + 5 = 11$ ✓.

**A2.**
$$E = E_{31}(4) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -4 & 0 & 1 \end{pmatrix}, \quad E^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 4 & 0 & 1 \end{pmatrix}; \qquad
P = P_{12} = \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad P^{-1} = P;$$
$$D = D_2(\tfrac15) = \operatorname{diag}(1, \tfrac15, 1), \quad D^{-1} = \operatorname{diag}(1, 5, 1).$$

$EA$: rows 1, 2 unchanged; row 3 $= (7, 8, 9) - 4(1, 2, 3) = (3, 0, -3)$. So $EA = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 3 & 0 & -3 \end{pmatrix}$.

$AE$: by Theorem 3.1, column $j$ of $AE$ is $\sum_k E_{kj}\,(\text{column } k \text{ of } A)$. Column 1 of $E$ is $(1, 0, -4)^\top$, so column 1 of $AE$ is $(\text{col } 1) - 4(\text{col } 3) = (1 - 12, 4 - 24, 7 - 36)^\top = (-11, -20, -29)^\top$. Columns 2 and 3 of $E$ are $\mathbf{e}_2, \mathbf{e}_3$, so columns 2, 3 of $AE$ are columns 2, 3 of $A$. $AE = \begin{pmatrix} -11 & 2 & 3 \\ -20 & 5 & 6 \\ -29 & 8 & 9 \end{pmatrix}$. In words: $E$ on the left subtracted $4\times$row 1 from row 3; $E$ on the right subtracted $4\times$column 3 from column 1. (Check one entry: $(AE)_{11} = \sum_k a_{1k}E_{k1} = 1\cdot 1 + 2\cdot 0 + 3\cdot(-4) = -11$. ✓)

**A3.** Column 1: pivot $1$. $\ell_{21} = 2$: row 2 ← $(2, 4, 0, 2) - 2(1, 2, 1, 3) = (0, 0, -2, -4)$. $\ell_{31} = 1$: row 3 ← $(1, 2, -1, -1) - (1, 2, 1, 3) = (0, 0, -2, -4)$. Column 2: positions $(2,2), (3,2)$ are both $0$: free. Column 3: pivot $-2$ in row 2. $\ell_{32} = (-2)/(-2) = 1$: row 3 ← $(0,0,-2,-4) - (0,0,-2,-4) = \mathbf{0}$.

Echelon form: $\begin{pmatrix} 1 & 2 & 1 & 3 \\ 0 & 0 & -2 & -4 \\ 0 & 0 & 0 & 0 \end{pmatrix}$. To RREF: row 2 $\div(-2)$ gives $(0, 0, 1, 2)$; row 1 ← row 1 − row 2 gives $(1, 2, 0, 1)$.

$$\operatorname{RREF}(M) = \begin{pmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 \end{pmatrix}.$$

Pivot columns $1, 3$; free columns $2, 4$; rank $2$. For $M\mathbf{x} = \mathbf{0}$: set $x_2 = s$, $x_4 = t$. Row 2: $x_3 + 2t = 0 \Rightarrow x_3 = -2t$. Row 1: $x_1 + 2s + t = 0 \Rightarrow x_1 = -2s - t$.

$$\mathbf{x} = s\begin{pmatrix} -2 \\ 1 \\ 0 \\ 0 \end{pmatrix} + t\begin{pmatrix} -1 \\ 0 \\ -2 \\ 1 \end{pmatrix}.$$

Check $\mathbf{s}_2 = (-1, 0, -2, 1)^\top$ against the original rows: row 1: $-1 + 0 - 2 + 3 = 0$ ✓; row 2: $-2 + 0 + 0 + 2 = 0$ ✓; row 3: $-1 + 0 + 2 - 1 = 0$ ✓.

**A4.** $L = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ -1 & -3 & 1 \end{pmatrix}$, $U = \begin{pmatrix} 2 & 1 & 1 \\ 0 & -1 & -2 \\ 0 & 0 & -4 \end{pmatrix}$, $\mathbf{b} = (2, 7, -5)^\top$.

Forward: $c_1 = 2$; $c_2 = 7 - 2\cdot 2 = 3$; $c_3 = -5 - (-1)(2) - (-3)(3) = -5 + 2 + 9 = 6$. $\mathbf{c} = (2, 3, 6)^\top$.

Back: $x_3 = 6/(-4) = -\tfrac32$; $x_2 = (3 - (-2)(-\tfrac32))/(-1) = (3 - 3)/(-1) = 0$; $x_1 = (2 - 1\cdot 0 - 1\cdot(-\tfrac32))/2 = \tfrac{7/2}{2} = \tfrac74$. $\mathbf{x} = (\tfrac74, 0, -\tfrac32)^\top$.

Check against $A = \begin{pmatrix} 2 & 1 & 1 \\ 4 & 1 & 0 \\ -2 & 2 & 1 \end{pmatrix}$: row 1: $\tfrac72 + 0 - \tfrac32 = 2$ ✓; row 2: $7 + 0 + 0 = 7$ ✓; row 3: $-\tfrac72 + 0 - \tfrac32 = -5$ ✓.

Flops. Forward: row 1 none; row 2 one multiply, one subtract ($2$); row 3 two multiplies, two subtracts ($4$). Total $6 = n(n-1) = 3\cdot 2$. ✓ Back: row 3 one divide ($1$); row 2 one multiply, one subtract, one divide ($3$); row 1 two multiplies, two subtracts, one divide ($5$). Total $9 = n^2$. ✓

**A5.** $[A\,|\,I] = \left[\begin{array}{cc|cc} 2 & 3 & 1 & 0 \\ 1 & 4 & 0 & 1 \end{array}\right]$. $\ell_{21} = \tfrac12$: row 2 ← row 2 − $\tfrac12$ row 1 $= (0, \tfrac52\,|\,-\tfrac12, 1)$. Up: row 1 ← row 1 − $\tfrac{3}{5/2}$ row 2 $=$ row 1 − $\tfrac65$ row 2 $= (2, 0\,|\,1 + \tfrac35, -\tfrac65) = (2, 0\,|\,\tfrac85, -\tfrac65)$. Scale: row 1 $\div 2 = (1, 0\,|\,\tfrac45, -\tfrac35)$; row 2 $\div\tfrac52 = (0, 1\,|\,-\tfrac15, \tfrac25)$.

$$A^{-1} = \frac15\begin{pmatrix} 4 & -3 \\ -1 & 2 \end{pmatrix}.$$

Formula: $ad - bc = 8 - 3 = 5$, $\frac15\begin{pmatrix} 4 & -3 \\ -1 & 2 \end{pmatrix}$. ✓ Then $A^\top = \begin{pmatrix} 2 & 1 \\ 3 & 4 \end{pmatrix}$, determinant still $5$, $(A^\top)^{-1} = \frac15\begin{pmatrix} 4 & -1 \\ -3 & 2 \end{pmatrix}$, which is the transpose of $A^{-1}$. ✓

**A6.** The rows of $P$ are $\mathbf{e}_3^\top, \mathbf{e}_1^\top, \mathbf{e}_2^\top$, so by Theorem 3.1 the rows of $PA$ are $\mathbf{r}_3^\top, \mathbf{r}_1^\top, \mathbf{r}_2^\top$: the last row moves to the top and the others shift down (a cyclic shift).

$P^\top = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$. $P^2$: $P\mathbf{e}_1 = \mathbf{e}_2$, $P\mathbf{e}_2 = \mathbf{e}_3$, $P\mathbf{e}_3 = \mathbf{e}_1$, so $P^2\mathbf{e}_1 = \mathbf{e}_3$, $P^2\mathbf{e}_2 = \mathbf{e}_1$, $P^2\mathbf{e}_3 = \mathbf{e}_2$: $P^2 = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix} = P^\top$. $P^3\mathbf{e}_j = \mathbf{e}_j$ for each $j$, so $P^3 = I$. Hence $P^{-1} = P^2 = P^\top$ (Theorem 3.23 confirmed; shifting three times is the identity, so shifting twice undoes shifting once).

## B. Core

**B1.** *(R1)–(R3) via Theorem 3.1.* Row $i$ of $E_{ij}(\ell) = I - \ell\mathbf{e}_i\mathbf{e}_j^\top$ is $\mathbf{e}_i^\top - \ell\mathbf{e}_j^\top$ (the matrix $\mathbf{e}_i\mathbf{e}_j^\top$ has row $i$ equal to $\mathbf{e}_j^\top$ and all other rows zero). By Theorem 3.1, row $i$ of $E_{ij}(\ell)A$ is $1\cdot\mathbf{r}_i^\top - \ell\,\mathbf{r}_j^\top$; every other row $k$ of $E_{ij}(\ell)$ is $\mathbf{e}_k^\top$, giving row $k$ of $A$ unchanged. That is (R1). For $P_{ij}$: rows $i$ and $j$ of $I$ are replaced by $\mathbf{e}_j^\top$ and $\mathbf{e}_i^\top$ respectively ($-\mathbf{e}_i\mathbf{e}_i^\top + \mathbf{e}_i\mathbf{e}_j^\top$ changes row $i$ from $\mathbf{e}_i^\top$ to $\mathbf{e}_j^\top$, and likewise for row $j$), so row $i$ of $P_{ij}A$ is $\mathbf{r}_j^\top$ and row $j$ is $\mathbf{r}_i^\top$: (R2). For $D_i(c)$: row $i$ is $c\,\mathbf{e}_i^\top$, so row $i$ of $D_i(c)A$ is $c\,\mathbf{r}_i^\top$: (R3).

*$P_{ij}^2 = I$ by direct multiplication.* Write $A_{ab} := \mathbf{e}_a\mathbf{e}_b^\top$, so (3.2) reads $A_{ab}A_{cd} = \delta_{bc}A_{ad}$. Then $P_{ij} = I + Q$ with $Q = -A_{ii} - A_{jj} + A_{ij} + A_{ji}$, and $P_{ij}^2 = I + 2Q + Q^2$. Expand $Q^2$ term by term (sixteen products; with $i \ne j$, $\delta_{ij} = 0$):

- $(-A_{ii})(-A_{ii}) = A_{ii}$; $(-A_{ii})(-A_{jj}) = \delta_{ij}A_{ij} = 0$; $(-A_{ii})(A_{ij}) = -A_{ij}$; $(-A_{ii})(A_{ji}) = -\delta_{ij}A_{ii} = 0$.
- $(-A_{jj})(-A_{ii}) = 0$; $(-A_{jj})(-A_{jj}) = A_{jj}$; $(-A_{jj})(A_{ij}) = -\delta_{ji}A_{jj} = 0$; $(-A_{jj})(A_{ji}) = -A_{ji}$.
- $(A_{ij})(-A_{ii}) = -\delta_{ji}A_{ii} = 0$; $(A_{ij})(-A_{jj}) = -A_{ij}$; $(A_{ij})(A_{ij}) = \delta_{ji}A_{ij} = 0$; $(A_{ij})(A_{ji}) = A_{ii}$.
- $(A_{ji})(-A_{ii}) = -A_{ji}$; $(A_{ji})(-A_{jj}) = -\delta_{ij}A_{jj} = 0$; $(A_{ji})(A_{ij}) = A_{jj}$; $(A_{ji})(A_{ji}) = 0$.

Sum: $2A_{ii} + 2A_{jj} - 2A_{ij} - 2A_{ji} = -2Q$. Therefore $P_{ij}^2 = I + 2Q - 2Q = I$. ∎

**B2.** (a) Column 1, pivot $1$: $\ell_{21} = 3$, row 2 ← $(3, 4, 1) - 3(1, 2, 0) = (0, -2, 1)$; $\ell_{31} = 2$, row 3 ← $(2, 5, -1) - 2(1, 2, 0) = (0, 1, -1)$. Column 2, pivot $-2$: $\ell_{32} = 1/(-2) = -\tfrac12$, row 3 ← $(0, 1, -1) - (-\tfrac12)(0, -2, 1) = (0, 1, -1) + (0, -1, \tfrac12) = (0, 0, -\tfrac12)$.

$$L = \begin{pmatrix} 1 & 0 & 0 \\ 3 & 1 & 0 \\ 2 & -\tfrac12 & 1 \end{pmatrix}, \qquad U = \begin{pmatrix} 1 & 2 & 0 \\ 0 & -2 & 1 \\ 0 & 0 & -\tfrac12 \end{pmatrix}.$$

Row 3 of $LU$ by Theorem 3.10: $2(1, 2, 0) - \tfrac12(0, -2, 1) + (0, 0, -\tfrac12) = (2, 4 + 1, -\tfrac12 - \tfrac12) = (2, 5, -1)$. ✓

(b) Forward: $c_1 = -1$; $c_2 = 1 - 3(-1) = 4$; $c_3 = -5 - 2(-1) - (-\tfrac12)(4) = -5 + 2 + 2 = -1$. Back: $x_3 = (-1)/(-\tfrac12) = 2$; $x_2 = (4 - 1\cdot 2)/(-2) = -1$; $x_1 = -1 - 2(-1) - 0 = 1$. $\mathbf{x} = (1, -1, 2)^\top$. Check: $(1 - 2 + 0, 3 - 4 + 2, 2 - 5 - 2) = (-1, 1, -5)$. ✓

(c) $E_{21} = E_{21}(3)$, $E_{31} = E_{31}(2)$, $E_{32} = E_{32}(-\tfrac12)$. $E_{31}E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -3 & 1 & 0 \\ -2 & 0 & 1 \end{pmatrix}$ (column 1 holds $-\ell_{21}, -\ell_{31}$; no cross term because both act on column 1 of $I$). Then $E_{32}(E_{31}E_{21})$ adds $\tfrac12\times$row 2 to row 3: row 3 $= (-2, 0, 1) + \tfrac12(-3, 1, 0) = (-\tfrac72, \tfrac12, 1)$.

$$E_{32}E_{31}E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -3 & 1 & 0 \\ -\tfrac72 & \tfrac12 & 1 \end{pmatrix}.$$

Its $(3,1)$ entry is $-\tfrac72$, while $L_{31} = 2$. The difference from the naive guess $-\ell_{31} = -2$ is the cross term $\ell_{32}\ell_{21} = (-\tfrac12)(3) = -\tfrac32$: $-\ell_{31} + \ell_{32}\ell_{21} = -2 - \tfrac32 = -\tfrac72$. ✓

(d) $D = \operatorname{diag}(1, -2, -\tfrac12)$, $U' = D^{-1}U = \begin{pmatrix} 1 & 2 & 0 \\ 0 & 1 & -\tfrac12 \\ 0 & 0 & 1 \end{pmatrix}$. $A = LDU'$.

**B3.** (a) Let $A$ be $m\times n$ with $n > m$. Run Algorithm 3.5: $EA = U$ with $E$ a product of elementary matrices, hence invertible (each elementary matrix is invertible by §3.1 and a product of invertible matrices is invertible). Each pivot lies in a distinct row, so $U$ has at most $m$ pivots, hence at most $m$ pivot columns, hence at least $n - m \ge 1$ free columns. Choose a free column $f$. Set $x_f = 1$ and every other free variable to $0$. Determine the pivot variables from the bottom pivot row upward: pivot row $k$ reads $u_{k,p_k}x_{p_k} + \sum_{j > p_k}u_{kj}x_j = 0$, every $x_j$ with $j > p_k$ is already fixed (a free variable, or a pivot variable from a lower row), and $u_{k,p_k} \ne 0$, so $x_{p_k} = -\tfrac{1}{u_{k,p_k}}\sum_{j>p_k}u_{kj}x_j$ is determined. Zero rows read $0 = 0$. The vector $\mathbf{x}$ built this way satisfies $U\mathbf{x} = \mathbf{0}$ and has $x_f = 1$, so $\mathbf{x} \ne \mathbf{0}$. Finally $A\mathbf{x} = E^{-1}U\mathbf{x} = E^{-1}\mathbf{0} = \mathbf{0}$. ∎

(b) $A = \begin{pmatrix} 0 & 2 & 1 \\ 1 & 1 & 0 \\ 2 & 3 & 1 \end{pmatrix}$, permutation vector $\pi = (1, 2, 3)$.

*Step 1.* Column 1 entries $0, 1, 2$: largest is $2$ in row 3. Swap rows 1 and 3: working matrix $\begin{pmatrix} 2 & 3 & 1 \\ 1 & 1 & 0 \\ 0 & 2 & 1 \end{pmatrix}$, $\pi = (3, 2, 1)$. Multipliers $\ell_{21} = \tfrac12$, $\ell_{31} = 0$. Row 2 ← $(1, 1, 0) - \tfrac12(2, 3, 1) = (0, -\tfrac12, -\tfrac12)$. Row 3 unchanged. Column 1 of $L$ so far: $(1, \tfrac12, 0)^\top$.

*Step 2.* Column 2 at or below row 2: $-\tfrac12$ (row 2), $2$ (row 3). Largest is $2$: swap rows 2 and 3. Working matrix $\begin{pmatrix} 2 & 3 & 1 \\ 0 & 2 & 1 \\ 0 & -\tfrac12 & -\tfrac12 \end{pmatrix}$, $\pi = (3, 1, 2)$. **Swap the same rows in column 1 of $L$:** it becomes $(1, 0, \tfrac12)^\top$. Multiplier $\ell_{32} = (-\tfrac12)/2 = -\tfrac14$. Row 3 ← $(0, -\tfrac12, -\tfrac12) + \tfrac14(0, 2, 1) = (0, 0, -\tfrac14)$.

$$P = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\ (\text{rows } \mathbf{e}_{\pi(i)}^\top), \qquad L = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ \tfrac12 & -\tfrac14 & 1 \end{pmatrix}, \qquad U = \begin{pmatrix} 2 & 3 & 1 \\ 0 & 2 & 1 \\ 0 & 0 & -\tfrac14 \end{pmatrix}.$$

$PA$ has rows $3, 1, 2$ of $A$: $\begin{pmatrix} 2 & 3 & 1 \\ 0 & 2 & 1 \\ 1 & 1 & 0 \end{pmatrix}$. $LU$: row 1 $= (2, 3, 1)$; row 2 $= (0, 2, 1)$; row 3 $= \tfrac12(2, 3, 1) - \tfrac14(0, 2, 1) + (0, 0, -\tfrac14) = (1, \tfrac32 - \tfrac12, \tfrac12 - \tfrac14 - \tfrac14) = (1, 1, 0)$. $PA = LU$. ✓ (Had you forgotten to swap column 1 of $L$, row 3 of $LU$ would have been $(0, \dots)$ and row 2 would have been $(1, \dots)$: wrong rows.)

**B4.** (a) $\ell_{21} = 1$: row 2 ← $(1, 2, 3) - (1, 1, 1) = (0, 1, 2)$. $\ell_{31} = 2$: row 3 ← $(2, 3, 4) - 2(1, 1, 1) = (0, 1, 2)$. $\ell_{32} = 1$: row 3 ← $(0,1,2) - (0,1,2) = \mathbf{0}$.

$U = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \\ 0 & 0 & 0 \end{pmatrix}$, rank $2$, column 3 free. Null vector: $x_3 = 1$; row 2: $x_2 + 2 = 0 \Rightarrow x_2 = -2$; row 1: $x_1 - 2 + 1 = 0 \Rightarrow x_1 = 1$. $\mathbf{s} = (1, -2, 1)^\top$. Check: $A\mathbf{s} = (1 - 2 + 1, 1 - 4 + 3, 2 - 6 + 4) = \mathbf{0}$. ✓

(b) Same operations on $\mathbf{b}$: $c_1 = b_1$; $c_2 = b_2 - b_1$; $c_3 = (b_3 - 2b_1) - (b_2 - b_1) = b_3 - b_1 - b_2$. The last row of the augmented echelon form is $(0, 0, 0\,|\,b_3 - b_1 - b_2)$, so the system is solvable iff
$$b_3 = b_1 + b_2.$$
LA-1 language: the column space of $A$ is the plane $\{\mathbf{b} : b_3 = b_1 + b_2\}$ (each column satisfies it: $2 = 1 + 1$, $3 = 1 + 2$, $4 = 1 + 3$), and $A\mathbf{x} = \mathbf{b}$ is solvable iff $\mathbf{b}$ lies in that plane (LA-1 Theorem 3.6).

(c) $\mathbf{b} = (1, 2, 3)^\top$: $3 = 1 + 2$ ✓. $\mathbf{c} = (1, 1, 0)^\top$. Set $x_3 = t$. Row 2: $x_2 = 1 - 2t$. Row 1: $x_1 = 1 - (1 - 2t) - t = t$. So
$$\mathbf{x} = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix} + t\begin{pmatrix} 1 \\ -2 \\ 1 \end{pmatrix}.$$
Check $\mathbf{x}_p$: $A(0, 1, 0)^\top = $ column 2 $= (1, 2, 3)^\top$. ✓

**B5.** (a) Induction on $k$. $k = 1$ is the definition. If $A^{k-1}(A^{-1})^{k-1} = I$, then $A^k(A^{-1})^k = A^{k-1}(AA^{-1})(A^{-1})^{k-1} = A^{k-1}I(A^{-1})^{k-1} = I$. So $(A^{-1})^k$ is a right inverse of $A^k$, hence the inverse (Theorem 3.15). ∎

(b) Suppose $AB$ is invertible. If $B$ were singular, Corollary 3.18 gives $\mathbf{x} \ne \mathbf{0}$ with $B\mathbf{x} = \mathbf{0}$; then $(AB)\mathbf{x} = A\mathbf{0} = \mathbf{0}$ with $\mathbf{x} \ne \mathbf{0}$, so $AB$ is singular by Corollary 3.18, contradiction. So $B$ is invertible. Then $A = (AB)B^{-1}$ is a product of two invertible matrices, hence invertible (with inverse $B(AB)^{-1}$). ∎

(c) Multiply $AB = AC$ on the left by $A^{-1}$: $B = C$. Counterexample when $A$ is singular: $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$, $B = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$, $C = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$. Then $AB = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix} = AC$ but $B \ne C$. (Any $B - C$ in the null space of $A$, columnwise, works.)

**B6.** (a) $[A\,|\,I]$, downward pass: $\ell_{21} = 2$: row 2 ← $(2, 1, 3\,|\,0, 1, 0) - 2(1, 0, 2\,|\,1, 0, 0) = (0, 1, -1\,|\,-2, 1, 0)$. $\ell_{31} = 0$. $\ell_{32} = 1$: row 3 ← $(0, 1, 1\,|\,0, 0, 1) - (0, 1, -1\,|\,-2, 1, 0) = (0, 0, 2\,|\,2, -1, 1)$.

$$\left[\begin{array}{ccc|ccc} 1 & 0 & 2 & 1 & 0 & 0 \\ 0 & 1 & -1 & -2 & 1 & 0 \\ 0 & 0 & 2 & 2 & -1 & 1 \end{array}\right]$$

Upward pass: clear column 3 above the pivot $2$. Row 2 ← row 2 + $\tfrac12$ row 3 $= (0, 1, 0\,|\,-2 + 1, 1 - \tfrac12, 0 + \tfrac12) = (0, 1, 0\,|\,-1, \tfrac12, \tfrac12)$. Row 1 ← row 1 − $\tfrac22$ row 3 $=$ row 1 − row 3 $= (1, 0, 0\,|\,-1, 1, -1)$. Column 2 above its pivot is already $0$.

$$\left[\begin{array}{ccc|ccc} 1 & 0 & 0 & -1 & 1 & -1 \\ 0 & 1 & 0 & -1 & \tfrac12 & \tfrac12 \\ 0 & 0 & 2 & 2 & -1 & 1 \end{array}\right]$$

Scaling: row 3 $\div 2$. Result:

$$A^{-1} = \begin{pmatrix} -1 & 1 & -1 \\ -1 & \tfrac12 & \tfrac12 \\ 1 & -\tfrac12 & \tfrac12 \end{pmatrix}.$$

(b) $AA^{-1}$: row 1 of $A$ is $(1, 0, 2)$: against the columns of $A^{-1}$: $(-1 + 2, 1 - 1, -1 + 1) = (1, 0, 0)$ ✓. Row 2, $(2, 1, 3)$: $(-2 - 1 + 3, 2 + \tfrac12 - \tfrac32, -2 + \tfrac12 + \tfrac32) = (0, 1, 0)$ ✓. Row 3, $(0, 1, 1)$: $(-1 + 1, \tfrac12 - \tfrac12, \tfrac12 + \tfrac12) = (0, 0, 1)$ ✓.

(c) $L = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix}$, $U = \begin{pmatrix} 1 & 0 & 2 \\ 0 & 1 & -1 \\ 0 & 0 & 2 \end{pmatrix}$, pivots $1, 1, 2$, product $2$ ($= \det A$ in LA-6; check: `np.linalg.det` gives $2$).

(d) The row operations of the whole reduction are one invertible matrix $M$ acting on the left of the block matrix, so $M[A\,|\,I] = [MA\,|\,M]$. We stopped when $MA = I$, which forces $M = A^{-1}$ (Theorem 3.15: a left inverse of a square matrix is the inverse), and the right block is exactly $M$.

**B7.** (a) Let $\mathbf{y} = L^\top\mathbf{x}$. $L^\top$ is unit upper triangular, hence invertible (Lemma 3.11(ii), transposed), so $\mathbf{x} \ne \mathbf{0}$ implies $\mathbf{y} \ne \mathbf{0}$. Then
$$\mathbf{x}^\top A\mathbf{x} = \mathbf{x}^\top LDL^\top\mathbf{x} = (L^\top\mathbf{x})^\top D(L^\top\mathbf{x}) = \mathbf{y}^\top D\mathbf{y} = \sum_{i=1}^n d_i y_i^2,$$
using $(L^\top\mathbf{x})^\top = \mathbf{x}^\top L$ (LA-1 Theorem 3.8). Every term is $\ge 0$ and at least one $y_i \ne 0$ has $d_iy_i^2 > 0$. So $\mathbf{x}^\top A\mathbf{x} > 0$. ∎

(b) For Example 4.3, $L^\top = \begin{pmatrix} 1 & \tfrac12 & 0 \\ 0 & 1 & \tfrac23 \\ 0 & 0 & 1 \end{pmatrix}$, so $\mathbf{y} = (x_1 + \tfrac12x_2,\ x_2 + \tfrac23x_3,\ x_3)^\top$ and $D = \operatorname{diag}(2, \tfrac32, \tfrac43)$:
$$\mathbf{x}^\top A\mathbf{x} = 2\big(x_1 + \tfrac12x_2\big)^2 + \tfrac32\big(x_2 + \tfrac23x_3\big)^2 + \tfrac43x_3^2.$$
Expand: $2(x_1^2 + x_1x_2 + \tfrac14x_2^2) = 2x_1^2 + 2x_1x_2 + \tfrac12x_2^2$; $\tfrac32(x_2^2 + \tfrac43x_2x_3 + \tfrac49x_3^2) = \tfrac32x_2^2 + 2x_2x_3 + \tfrac23x_3^2$; plus $\tfrac43x_3^2$. Total: $2x_1^2 + 2x_1x_2 + (\tfrac12 + \tfrac32)x_2^2 + 2x_2x_3 + (\tfrac23 + \tfrac43)x_3^2 = 2x_1^2 + 2x_1x_2 + 2x_2^2 + 2x_2x_3 + 2x_3^2$. ✓ (Directly: $\mathbf{x}^\top A\mathbf{x} = \sum_{ij}a_{ij}x_ix_j = 2x_1^2 + 2x_2^2 + 2x_3^2 + 2x_1x_2 + 2x_2x_3$.)

(c) Suppose $d_k \le 0$. Take $\mathbf{y} = \mathbf{e}_k$ and $\mathbf{x} = (L^\top)^{-1}\mathbf{e}_k$, which is nonzero because $(L^\top)^{-1}$ is invertible. Then $L^\top\mathbf{x} = \mathbf{e}_k$ and $\mathbf{x}^\top A\mathbf{x} = \mathbf{e}_k^\top D\mathbf{e}_k = d_k \le 0$. ∎ (So "all pivots positive" is equivalent to "$\mathbf{x}^\top A\mathbf{x} > 0$ for all $\mathbf{x} \ne \mathbf{0}$," which LA-8 takes as the definition of positive definite.)

**B8.** (a) At step $k$ there are $m = n - k$ rows below the pivot row. Each needs one division for its multiplier: $m$ flops. Each has $m$ entries to the right of the pivot column that change, each by one multiply and one subtract: $2m^2$ flops. Summing over $k = 1, \dots, n-1$, i.e. $m = n-1, \dots, 1$:
$$\sum_{m=1}^{n-1}(m + 2m^2) = \frac{(n-1)n}{2} + 2\cdot\frac{(n-1)n(2n-1)}{6} = \frac{n(n-1)}{2} + \frac{n(n-1)(2n-1)}{3}.$$
For $n = 3$: $3 + \frac{3\cdot 2\cdot 5}{3} = 3 + 10 = 13$. By hand in Example 4.1: step 1 had $2$ divisions ($\ell_{21}, \ell_{31}$) and updated $2$ rows $\times$ $2$ entries $\times$ $2$ flops $= 8$, total $10$; step 2 had $1$ division and $1\times 1\times 2 = 2$, total $3$. $10 + 3 = 13$. ✓ (Expanding, the total is $\tfrac23n^3 - \tfrac12n^2 - \tfrac16n$, so the ratio to $n^3$ approaches $\tfrac23$ from below.)

(b) Set $k\cdot\tfrac23n^3 = 100\big(\tfrac23n^3 + 2kn^2\big)$. Then $k\big(\tfrac23n^3 - 200n^2\big) = \tfrac{200}{3}n^3$, so $k = \dfrac{200n^3/3}{2n^3/3 - 200n^2} = \dfrac{200n/3}{2n/3 - 200}$. For $n = 500$: numerator $\tfrac{200\cdot 500}{3} = 33333.3$, denominator $\tfrac{1000}{3} - 200 = 133.3$, so $k = 250$. Check: factor-once cost $\tfrac23(500)^3 + 2\cdot 250\cdot 500^2 = 8.33\times 10^7 + 1.25\times 10^8 = 2.08\times 10^8$; from scratch $250\times 8.33\times 10^7 = 2.08\times 10^{10}$; ratio $100$. ✓

(c) Forward substitution $L\mathbf{c} = \mathbf{e}_j$: for $i < j$, $c_i = 0 - \sum_{k<i}\ell_{ik}c_k = 0$ by induction (all earlier $c_k$ are zero), with no arithmetic needed once you know it; $c_j = 1$; for $i > j$, $c_i = -\sum_{k=j}^{i-1}\ell_{ik}c_k$ costs about $2(i - j)$ flops. Total for one $\mathbf{e}_j$: $\sum_{i=j+1}^n 2(i-j) \approx (n-j)^2$. Over all $j$: $\sum_{j=1}^n(n-j)^2 \approx n^3/3$. The $n$ back substitutions have no such zeros and cost $n\cdot n^2 = n^3$. Total for $A^{-1}$: $\tfrac23n^3 + \tfrac13n^3 + n^3 = 2n^3$, three times the $\tfrac23n^3$ of one factorization, with nothing gained afterwards since $A^{-1}\mathbf{b}$ costs $2n^2$, the same as two substitutions.

**B9.** Reference implementation (also in `project.py`, "Your turn"):

```python
def ldlt(A):
    """A = L diag(d) L^T for symmetric A without row swaps (Theorem 3.22)."""
    A = np.asarray(A, float); n = A.shape[0]
    assert np.allclose(A, A.T)
    L = np.eye(n); d = np.zeros(n)
    for j in range(n):
        d[j] = A[j, j] - np.sum(L[j, :j] ** 2 * d[:j])
        if d[j] == 0:
            raise np.linalg.LinAlgError("zero pivot")
        for i in range(j + 1, n):
            L[i, j] = (A[i, j] - np.sum(L[i, :j] * L[j, :j] * d[:j])) / d[j]
    return L, d
```

Why these formulas: entry $(i, j)$ of $LDL^\top$ with $i \ge j$ is $\sum_k L_{ik}d_kL_{jk} = \sum_{k<j}L_{ik}d_kL_{jk} + L_{ij}d_j$ (since $L_{jk} = 0$ for $k > j$ and $L_{jj} = 1$). Setting this equal to $a_{ij}$ and solving gives $d_j$ (case $i = j$) and then $L_{ij}$. Only $a_{ij}$ with $i \ge j$ are read. The checks: `L @ np.diag(d) @ L.T` reproduces $BB^\top + I$ to `1e-10`; for Example 4.3, `d = [2, 1.5, 1.333...]`, `L[1,0] = 0.5`, `L[2,1] = 0.6667`; and `d` equals `np.diag(U)` from `lu_decompose(A, pivot=False)` because Theorem 3.22 says the pivots of a symmetric matrix are its $D$.

Tolerance: `lu_decompose` raises when `abs(U[k,k]) <= tol`. On Example 4.2 with `tol=1e-12` it raises at $k = 1$ (the whole column below the staircase is zero to within rounding). On a random $BC$ of rank 2 it raises with `tol=1e-10` but not with `tol=0`: the third pivot is a rounding artifact of size about $10^{-16}$ (the script prints it), and `solve` divides by it, returning an $\mathbf{x}$ with entries of order $10^{15}$–$10^{16}$ whose residual $A\mathbf{x} - \mathbf{b}$ is *not* small. The code is "solving" $(A + \delta A)\mathbf{x} = \mathbf{b}$ for a $\delta A$ at rounding level that happens to make the matrix invertible; the answer means nothing.

## C. Challenge

**C1.** (a) Step 1: column 1 is $(1, -1, -1, -1)^\top$; all candidates have $|a_{i1}| = 1$, so partial pivoting keeps the diagonal (no swap). Multipliers $\ell_{i1} = -1/1 = -1$ for $i = 2, 3, 4$: row $i$ ← row $i$ + row 1. Row 1 is $(1, 0, 0, 1)$, so this zeroes column 1 below the diagonal, leaves columns 2, 3 alone, and doubles the last column: rows $2, 3, 4$ become $(0, 1, 0, 2)$, $(0, -1, 1, 2)$, $(0, -1, -1, 2)$. The trailing $3\times 3$ block is $W_3$ with last column $2\cdot(1,1,1)^\top$. Step 2: same pattern, multipliers $-1$, last column of rows 3, 4 becomes $4$: rows $(0, 0, 1, 4)$, $(0, 0, -1, 4)$. Step 3: row 4 ← row 4 + row 3 $= (0, 0, 0, 8)$.

$$U = \begin{pmatrix} 1 & 0 & 0 & 1 \\ 0 & 1 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 8 \end{pmatrix}, \qquad L = \begin{pmatrix} 1 & 0 & 0 & 0 \\ -1 & 1 & 0 & 0 \\ -1 & -1 & 1 & 0 \\ -1 & -1 & -1 & 1 \end{pmatrix}, \qquad P = I.$$

(b) Claim: after $k$ steps of elimination on $W_n$ (no swaps), rows $k+1, \dots, n$ form the matrix $W_{n-k}$ in columns $k+1, \dots, n$ except that its last column is $2^k(1, \dots, 1)^\top$, and the first $k$ rows of $U$ have last-column entries $1, 2, \dots, 2^{k-1}$. Base $k = 0$: $W_n$ itself. Inductive step: the pivot row $k+1$ is $(0, \dots, 0, 1, 0, \dots, 0, 2^k)$ (a $1$ on the diagonal, zeros until the last column). Every entry below the diagonal in column $k+1$ is $-1$, so all candidates tie at $|{-1}| = 1$ and no swap occurs; the multipliers are $-1$. Adding the pivot row to each lower row zeroes column $k+1$, changes nothing in columns $k+2, \dots, n-1$ (the pivot row is zero there), and adds $2^k$ to the last-column entries, turning $2^k$ into $2^{k+1}$. The trailing block is now $W_{n-k-1}$ with last column $2^{k+1}(1, \dots, 1)^\top$, and the pivot row contributed $u_{k+1,n} = 2^k$. After $n-1$ steps the last column of $U$ is $(1, 2, \dots, 2^{n-1})^\top$. ∎ All entries of $W_n$ have absolute value $\le 1$ and $\max|u_{ij}| = 2^{n-1}$, so the growth factor is $2^{n-1}$, which is the maximum that $|\ell_{ij}| \le 1$ permits (each step can at most double the largest entry).

(c) With `scipy.linalg.lu_factor`/`lu_solve` on $W_{60}$ and $\mathbf{b} = W_{60}\mathbf{1}$ the project prints: $\max|u_{ij}| = 2^{59} \approx 5.76\times 10^{17}$, $\max_i|x_i - 1| = 1$ (some components of the computed solution are off by order one), and relative residual $\|W\mathbf{x} - \mathbf{b}\|/\|\mathbf{b}\| \approx 3\times 10^{-2}$ rather than $10^{-16}$; your own `solve` gives numbers of the same size. The explanation: backward stability of $LU$ bounds the backward error by roughly (machine precision) $\times$ (growth factor) $\times \|A\|$, and here $1.1\times 10^{-16}\times 5.8\times 10^{17} \approx 64$. The rounding errors made while forming the last column of $U$ are of size $2^{59}\varepsilon_{\text{mach}} \approx 60$ in absolute terms, against entries of $A$ of size $1$. Elimination has effectively solved a system whose matrix differs from $W_{60}$ by $O(10)$: a different problem. For $n \le 50$ the growth $2^{49}\varepsilon_{\text{mach}} \approx 6\times 10^{-2}$ is still below $1$ and the solution comes out exact; the failure switches on between $n = 50$ and $n = 55$. This is the matrix that shows partial pivoting is not unconditionally stable (LA-13); it essentially never occurs in practice, which is why everyone uses partial pivoting anyway.

**C2.** (a) Block multiplication (rows of the first factor times columns of the second, with blocks in place of entries, which is legitimate because block sizes match):
$$\begin{pmatrix} I & 0 \\ A_{21}A_{11}^{-1} & I \end{pmatrix}\begin{pmatrix} A_{11} & A_{12} \\ 0 & S \end{pmatrix} = \begin{pmatrix} A_{11} & A_{12} \\ A_{21}A_{11}^{-1}A_{11} & A_{21}A_{11}^{-1}A_{12} + S \end{pmatrix} = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}. \ ∎$$

(b) Call the two factors $\mathcal{L}$ and $\mathcal{T}$. $\mathcal{L}$ is invertible with $\mathcal{L}^{-1} = \begin{pmatrix} I & 0 \\ -A_{21}A_{11}^{-1} & I \end{pmatrix}$ (multiply: the off-diagonal block is $A_{21}A_{11}^{-1} - A_{21}A_{11}^{-1} = 0$). So $A = \mathcal{L}\mathcal{T}$ is invertible iff $\mathcal{T}$ is (if $\mathcal{T}$ is invertible, $A^{-1} = \mathcal{T}^{-1}\mathcal{L}^{-1}$; if $A$ is invertible, $\mathcal{T} = \mathcal{L}^{-1}A$ is). If $S$ is invertible, then $\mathcal{T}^{-1} = \begin{pmatrix} A_{11}^{-1} & -A_{11}^{-1}A_{12}S^{-1} \\ 0 & S^{-1} \end{pmatrix}$: check $\mathcal{T}\mathcal{T}^{-1}$ has top-right block $-A_{12}S^{-1} + A_{12}S^{-1} = 0$ and diagonal blocks $I$. If $S$ is singular, there is $\mathbf{y} \ne \mathbf{0}$ with $S\mathbf{y} = \mathbf{0}$ (Corollary 3.18); set $\mathbf{x} = \begin{pmatrix} -A_{11}^{-1}A_{12}\mathbf{y} \\ \mathbf{y} \end{pmatrix} \ne \mathbf{0}$. Then $\mathcal{T}\mathbf{x} = \begin{pmatrix} -A_{12}\mathbf{y} + A_{12}\mathbf{y} \\ S\mathbf{y} \end{pmatrix} = \mathbf{0}$, so $\mathcal{T}$ is singular. Hence $A$ invertible $\iff$ $\mathcal{T}$ invertible $\iff$ $S$ invertible. ∎

(c) $G = \begin{pmatrix} \tfrac32 & -\tfrac12 & 0 \\ -\tfrac12 & 1 & -\tfrac12 \\ 0 & -\tfrac12 & \tfrac32 \end{pmatrix}$. With $A_{11} = \tfrac32$, $A_{21} = (-\tfrac12, 0)^\top$, $A_{12} = (-\tfrac12, 0)$, $A_{22} = \begin{pmatrix} 1 & -\tfrac12 \\ -\tfrac12 & \tfrac32 \end{pmatrix}$:
$$S = A_{22} - \frac{1}{3/2}\begin{pmatrix} -\tfrac12 \\ 0 \end{pmatrix}\begin{pmatrix} -\tfrac12 & 0 \end{pmatrix} = A_{22} - \frac23\begin{pmatrix} \tfrac14 & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 - \tfrac16 & -\tfrac12 \\ -\tfrac12 & \tfrac32 \end{pmatrix} = \begin{pmatrix} \tfrac56 & -\tfrac12 \\ -\tfrac12 & \tfrac32 \end{pmatrix}.$$
This $S$ is exactly the trailing $2\times 2$ block after step 1 of Example 4.4, and $S_{11} = \tfrac56 = \tfrac13 + \tfrac12$: eliminating node 1 merges $R_1$ (to ground) and $R_2$ (to node 2) into a single $1 + 2 = 3\,\Omega$ path from node 2 to ground (conductance $\tfrac13$), in parallel with $R_3$'s conductance $\tfrac12$ toward node 3. $S$ is the conductance matrix of the two-node circuit that remains.

With $A_{11} = \begin{pmatrix} \tfrac32 & -\tfrac12 \\ -\tfrac12 & 1 \end{pmatrix}$ (determinant $\tfrac32 - \tfrac14 = \tfrac54$, inverse $\tfrac45\begin{pmatrix} 1 & \tfrac12 \\ \tfrac12 & \tfrac32 \end{pmatrix}$), $A_{12} = (0, -\tfrac12)^\top$, $A_{21} = (0, -\tfrac12)$, $A_{22} = \tfrac32$:
$$A_{21}A_{11}^{-1}A_{12} = \tfrac45\cdot(-\tfrac12)\cdot\tfrac32\cdot(-\tfrac12) = \tfrac45\cdot\tfrac38 = \tfrac{3}{10}, \qquad S = \tfrac32 - \tfrac{3}{10} = \tfrac65.$$
That is the third pivot of Example 4.4: $\tfrac15 + 1$, the series path $R_3 + R_2 + R_1 = 5\,\Omega$ to ground in parallel with $R_4$. Each pivot is a $1\times 1$ Schur complement.

**C3.** (a) Let $U$ be in echelon form with pivots in columns $p_1 < p_2 < \cdots < p_r$, pivot $k$ in row $k$.

*Pivot column $\Rightarrow$ not a combination.* Suppose column $j = p_k$. Rows $k, k+1, \dots$ have their first nonzero entry in column $\ge p_k = j$ (row $k$'s pivot is at $j$; lower rows' pivots are further right or the rows are zero). So every column $1, \dots, j-1$ has a zero in row $k$, hence so does any linear combination of them, while column $j$ has $u_{kj} \ne 0$ in row $k$. Column $j$ is not a combination of the earlier columns.

*Free column $\Rightarrow$ a combination.* Suppose column $j$ is free, and let $q$ be the number of pivot columns to its left, $p_1 < \cdots < p_q < j$. Rows $q+1, q+2, \dots$ have their first nonzero in a pivot column $p_{q+1} > j$ (or are zero), so column $j$ and columns $1, \dots, j-1$ are all zero in rows $> q$. We seek $y_1, \dots, y_q$ with $\sum_i y_i\,(\text{column } p_i) = \text{column } j$. In rows $> q$ both sides are zero. In rows $1, \dots, q$ this is the $q\times q$ system $T\mathbf{y} = (u_{1j}, \dots, u_{qj})^\top$ with $T_{ai} = u_{a,p_i}$. $T$ is upper triangular with nonzero diagonal: $T_{aa} = u_{a,p_a}$ is a pivot, and for $a > i$, $u_{a,p_i} = 0$ because row $a$'s first nonzero is at $p_a > p_i$. Theorem 3.6 solves it by back substitution. Non-pivot columns to the left of $j$ get coefficient $0$. So column $j$ is a combination of columns $1, \dots, j-1$. ∎

(b) "Column $j$ of $A$ equals $\sum_{i<j}y_i\,(\text{column } i \text{ of } A)$" is the statement $A\mathbf{y}' = \mathbf{0}$ with $\mathbf{y}' = (y_1, \dots, y_{j-1}, -1, 0, \dots, 0)^\top$ (LA-1 eq. 3.1). By Theorem 3.3, $A\mathbf{y}' = \mathbf{0} \iff (EA)\mathbf{y}' = \mathbf{0}$, which is the same statement about the columns of $EA$ with the same coefficients. ∎

(c) Let $U = EA$ be any echelon form of $A$. By (a), the pivot columns of $U$ are exactly the columns $j$ that are not combinations of columns $1, \dots, j-1$ of $U$. By (b), that is exactly the set of columns $j$ of $A$ that are not combinations of columns $1, \dots, j-1$ of $A$, a set defined by $A$ alone with no reference to $E$. So every echelon form of $A$ has the same pivot columns, and in particular the same number of them. ∎ (LA-3 goes further: this count is $\dim C(A)$, and the columns of $A$ in the pivot positions form a basis of $C(A)$.)

**C4.** (a) Write $\mathbf{w} = A^{-1}\mathbf{u}$ and $\sigma = 1 + \mathbf{v}^\top\mathbf{w}$. Then
$$(A + \mathbf{u}\mathbf{v}^\top)\Big(A^{-1} - \frac{\mathbf{w}\mathbf{v}^\top A^{-1}}{\sigma}\Big) = I - \frac{\mathbf{u}\mathbf{v}^\top A^{-1}}{\sigma} + \mathbf{u}\mathbf{v}^\top A^{-1} - \frac{\mathbf{u}(\mathbf{v}^\top\mathbf{w})\mathbf{v}^\top A^{-1}}{\sigma},$$
using $A\mathbf{w} = \mathbf{u}$ in the second term and $\mathbf{v}^\top\mathbf{w}$ being a scalar in the fourth. Collect the three $\mathbf{u}\mathbf{v}^\top A^{-1}$ terms: coefficient $-\tfrac1\sigma + 1 - \tfrac{\mathbf{v}^\top\mathbf{w}}{\sigma} = \tfrac{-1 + \sigma - (\sigma - 1)}{\sigma} = 0$. The product is $I$, so by Theorem 3.15 the formula gives the inverse. ∎

(b) $\mathbf{x} = (A + \mathbf{u}\mathbf{v}^\top)^{-1}\mathbf{b} = A^{-1}\mathbf{b} - \dfrac{\mathbf{w}\,(\mathbf{v}^\top A^{-1}\mathbf{b})}{\sigma}$. Algorithm: (1) $\mathbf{z} = A^{-1}\mathbf{b}$ by forward then back substitution with $L, U$: $2n^2$ flops. (2) $\mathbf{w} = A^{-1}\mathbf{u}$ the same way: $2n^2$. (3) $\sigma = 1 + \mathbf{v}\cdot\mathbf{w}$, $\beta = \mathbf{v}\cdot\mathbf{z}$: $4n$. (4) $\mathbf{x} = \mathbf{z} - (\beta/\sigma)\mathbf{w}$: $2n$. Total about $4n^2$, versus $\tfrac23n^3$ for refactoring $A + \mathbf{u}\mathbf{v}^\top$. (If $\mathbf{w}$ is reused for many $\mathbf{b}$, each extra solve costs $2n^2 + O(n)$.)

(c) If $\sigma = 0$, let $\mathbf{w} = A^{-1}\mathbf{u}$, which is nonzero unless $\mathbf{u} = \mathbf{0}$ (and if $\mathbf{u} = \mathbf{0}$ then $\sigma = 1 \ne 0$, so $\mathbf{u} \ne \mathbf{0}$ here). Then $(A + \mathbf{u}\mathbf{v}^\top)\mathbf{w} = A\mathbf{w} + \mathbf{u}(\mathbf{v}^\top\mathbf{w}) = \mathbf{u} + \mathbf{u}(\sigma - 1) = \sigma\mathbf{u} = \mathbf{0}$. A nonzero null vector: singular by Corollary 3.18. ∎

(d)
```python
A = rng.standard_normal((6, 6)); u = rng.standard_normal(6); v = rng.standard_normal(6)
Ai = np.linalg.inv(A); sigma = 1 + v @ Ai @ u
assert np.allclose(np.linalg.inv(A + np.outer(u, v)), Ai - np.outer(Ai @ u, v @ Ai) / sigma)

n = 400; A = rng.standard_normal((n, n)); u, v, b = (rng.standard_normal(n) for _ in range(3))
P, L, U = lu_decompose(A)
t0 = time.perf_counter()
z = back_sub(U, forward_sub(L, P @ b)); w = back_sub(U, forward_sub(L, P @ u))
x_update = z - w * (v @ z) / (1 + v @ w)
t_update = time.perf_counter() - t0
t0 = time.perf_counter(); x_fresh = solve(A + np.outer(u, v), b); t_fresh = time.perf_counter() - t0
assert np.allclose(x_update, x_fresh)
print(t_update, t_fresh)   # a few milliseconds vs. tens of milliseconds with this script's LU
```
The update route is two pairs of substitutions ($O(n^2)$); the fresh route is a new factorization ($O(n^3)$). With the pure-NumPy `lu_decompose` of the project, the fresh solve of a $400\times 400$ system takes about $70\,\mathrm{ms}$ and the update about $3\,\mathrm{ms}$; the ratio grows linearly with $n$.
