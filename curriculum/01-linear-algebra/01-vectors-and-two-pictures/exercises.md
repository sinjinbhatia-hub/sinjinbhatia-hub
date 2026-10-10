# LA-1 Exercises

Conventions as in the notes: column vectors, $\mathbf{u}\cdot\mathbf{v} = \mathbf{u}^\top\mathbf{v}$, $\|\cdot\|$ is the Euclidean norm unless subscripted.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Let $A = \begin{pmatrix} 1 & 3 \\ -2 & 0 \\ 4 & 1 \end{pmatrix}$ and $\mathbf{x} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}$. Compute $A\mathbf{x}$ twice: once as a combination of columns (3.1), once row by row (3.2). Confirm they agree.

**A2.** Write the system $\;3x + y = 7,\; x - 2y = 0\;$ in the form $A\mathbf{x} = \mathbf{b}$, then write it in column form as $x\,\mathbf{a}_1 + y\,\mathbf{a}_2 = \mathbf{b}$. Solve it and verify the column equation numerically.

**A3.** Compute $\mathbf{u}\cdot\mathbf{v}$, $\|\mathbf{u}\|$, $\|\mathbf{v}\|$, and the angle between $\mathbf{u} = (1, 2, 2)^\top$ and $\mathbf{v} = (2, -1, 0)^\top$.

**A4.** Let $A = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$ and $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$. Compute $AB$ and $BA$ using the column rule (3.3). Are they equal? Describe geometrically what each of $A$, $B$, $AB$, $BA$ does to the plane.

**A5.** For $\mathbf{v} = (3, -4)^\top$, compute $\|\mathbf{v}\|_1$, $\|\mathbf{v}\|_2$, $\|\mathbf{v}\|_\infty$. Which is largest? Will that ordering always hold? (Just conjecture; C2 proves it.)

## B. Core (do all)

**B1.** Give a $2\times 2$ matrix $A$ and vector $\mathbf{b}$ such that $A\mathbf{x} = \mathbf{b}$ has (a) exactly one solution, (b) no solution, (c) infinitely many solutions. For each, draw (or describe) the row picture *and* the column picture.

**B2.** Prove Theorem 3.4 (linearity of $\mathbf{x} \mapsto A\mathbf{x}$) starting from the **row** form (3.2) instead of the column form. Then explain in one sentence why the two proofs had to give the same result.

**B3.** Let $A$ be $m \times n$. Prove that $A\mathbf{0} = \mathbf{0}$ and that if $A\mathbf{x} = \mathbf{b}$ and $A\mathbf{y} = \mathbf{b}$ then $A(\mathbf{x} - \mathbf{y}) = \mathbf{0}$. Conclude: if $A\mathbf{x} = \mathbf{b}$ has two different solutions, then $A\mathbf{z} = \mathbf{0}$ has a nonzero solution. (This is the seed of the null space in LA-3.)

**B4.** Let $\mathbf{v}_1 = (1, 0, 1)^\top$, $\mathbf{v}_2 = (0, 1, 1)^\top$, $\mathbf{v}_3 = (1, 1, 2)^\top$.
(a) Show that $\mathbf{v}_3 \in \operatorname{span}\{\mathbf{v}_1, \mathbf{v}_2\}$.
(b) Describe $\operatorname{span}\{\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3\}$ geometrically. Find a single equation $\mathbf{n}\cdot\mathbf{x} = 0$ that describes it.
(c) Is $(1, 1, 1)^\top$ in the span? Is $(2, 3, 5)^\top$?

**B5.** Prove that $\operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_k\}$ is closed under addition and scalar multiplication: if $\mathbf{p}, \mathbf{q}$ are in the span and $c \in \mathbb{R}$, then $\mathbf{p} + \mathbf{q}$ and $c\mathbf{p}$ are in the span.

**B6.** (Parallelogram law.) Prove that for all $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$,
$$\|\mathbf{u} + \mathbf{v}\|^2 + \|\mathbf{u} - \mathbf{v}\|^2 = 2\|\mathbf{u}\|^2 + 2\|\mathbf{v}\|^2.$$
Interpret geometrically (what are the four lengths?). Then show that the $\ell^1$ norm does **not** satisfy this identity by finding a counterexample in $\mathbb{R}^2$. (Why this matters: in LA-8 you'll learn that a norm comes from an inner product *if and only if* it satisfies the parallelogram law.)

**B7.** (Projection.) Let $\mathbf{v} \neq \mathbf{0}$. Define $f(t) = \|\mathbf{u} - t\mathbf{v}\|^2$.
(a) Expand $f(t)$ and find the $t^*$ that minimizes it by calculus.
(b) Show the residual $\mathbf{u} - t^*\mathbf{v}$ is orthogonal to $\mathbf{v}$.
(c) Show $\|\mathbf{u}\|^2 = \|t^*\mathbf{v}\|^2 + \|\mathbf{u} - t^*\mathbf{v}\|^2$ (Pythagoras for the projection).
(d) Compute the projection of $(2, 3, 6)^\top$ onto the line through $(1, 1, 1)^\top$ and verify (b) and (c) numerically.

**B8.** (Code.) In `project.py`, Part 3 implements `matvec_columns` and `matvec_rows`. Add a function `matmat_columns(A, B)` that computes $AB$ using (3.3), column by column, calling `matvec_columns`. Assert it agrees with `A @ B` on random $4\times 3$ and $3\times 5$ matrices. Then write a one-paragraph comment explaining why this definition makes $(AB)\mathbf{x} = A(B\mathbf{x})$ automatic.

## C. Challenge (do at least one)

**C1.** (Weighted dot product.) Let $w_1, \dots, w_n > 0$ and define $\langle \mathbf{u}, \mathbf{v}\rangle_w = \sum_i w_i u_i v_i$. Show this satisfies the three dot-product properties (symmetric, linear in each slot, positive definite), and prove the corresponding Cauchy–Schwarz inequality $|\langle \mathbf{u}, \mathbf{v}\rangle_w| \le \sqrt{\langle \mathbf{u}, \mathbf{u}\rangle_w}\sqrt{\langle \mathbf{v}, \mathbf{v}\rangle_w}$ using only those three properties (no geometry). Where do weighted inner products show up? (Hint: weighted least squares in ST-3; the Fisher information metric in ST-1.)

**C2.** (Norm inequalities.) Prove that for all $\mathbf{v} \in \mathbb{R}^n$:
$$\|\mathbf{v}\|_\infty \le \|\mathbf{v}\|_2 \le \|\mathbf{v}\|_1 \le \sqrt{n}\,\|\mathbf{v}\|_2 \le n\,\|\mathbf{v}\|_\infty.$$
For each inequality, give a vector that achieves equality. (The middle one, $\|\mathbf{v}\|_1 \le \sqrt{n}\|\mathbf{v}\|_2$, is Cauchy–Schwarz in disguise: find the disguise.)

**C3.** (A first matrix norm.) Define the Frobenius norm $\|A\|_F = \sqrt{\sum_{i,j} a_{ij}^2}$.
(a) Show $\|A\|_F^2 = \sum_j \|\mathbf{a}_j\|^2$ (sum of squared column lengths) $= \operatorname{tr}(A^\top A)$, where $\operatorname{tr}$ is the sum of the diagonal entries.
(b) Prove $\|A\mathbf{x}\| \le \|A\|_F\,\|\mathbf{x}\|$ using Cauchy–Schwarz on each row.
(c) Find a $2\times 2$ matrix and vector where (b) is strict and one where it's equality. Is $\|A\|_F$ the *smallest* constant $K$ with $\|A\mathbf{x}\| \le K\|\mathbf{x}\|$ for all $\mathbf{x}$? (The smallest such $K$ is the *operator norm*; LA-9 shows it's the largest singular value.)
