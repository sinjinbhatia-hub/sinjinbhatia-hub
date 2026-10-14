# LA-1: Vectors, the Two Pictures of $A\mathbf{x} = \mathbf{b}$, Dot Product, Norms, Span

**Time:** ~8 hours · **Prerequisites:** none (multivariable calculus helps) · **Unlocks:** LA-2, LA-3, and indirectly everything else

> **Conventions in this module.** Vectors are columns and written in bold lowercase: $\mathbf{v} \in \mathbb{R}^n$. Their entries are $v_1, \dots, v_n$. Matrices are uppercase: $A$. The dot product of $\mathbf{u}$ and $\mathbf{v}$ is $\mathbf{u} \cdot \mathbf{v} = \mathbf{u}^\top \mathbf{v} = \sum_i u_i v_i$. The length (Euclidean norm) is $\|\mathbf{v}\| = \sqrt{\mathbf{v}\cdot\mathbf{v}}$; $\|\mathbf{v}\|$ with no subscript is always this norm, and $\|\mathbf{v}\|_1$, $\|\mathbf{v}\|_\infty$ are defined in §3.6. The standard basis vectors are $\mathbf{e}_1, \dots, \mathbf{e}_n$: $\mathbf{e}_i$ has a 1 in position $i$ and zeros elsewhere. Products compose right to left: in $AB$, $B$ acts first.

## 0. Where this sits

You've used vectors in physics as "arrows with magnitude and direction" and in multivariable calculus as position vectors and gradients. This module deliberately starts over, because the picture you need for the rest of this curriculum is not the arrow picture. It is the picture of **a matrix acting on a vector**, and especially of **$A\mathbf{x}$ as a combination of the columns of $A$**.

If this module feels easy, good. Do the exercises anyway. The column picture is the single idea that LA-3 (subspaces), LA-4 (linear maps), LA-5 (projections), and LA-9 (SVD) are all built on, and most students who struggle in those modules struggle because they never actually internalized it here.

## 1. The problem that forced this idea

Linear equations are old. Old Babylonian scribes (roughly 1900–1600 BCE) solved problems that amount to two equations in two unknowns. Chapter 8 of the Chinese text *Nine Chapters on the Mathematical Art* (compiled by about the first century CE from older material) solves systems of linear equations, starting with three equations in three unknowns, by a method you would recognize as elimination, laid out on a counting board. But none of these cultures had *vectors* or *matrices*. They had recipes for numbers.

The thing that forced vectors into existence was geometry meeting algebra. Descartes (1637) showed that points in the plane can be coordinates $(x, y)$, so geometric questions become algebraic ones. But the reverse problem took two more centuries: given a system of linear equations, what does it *look like*? Where is the geometry?

Two answers emerged, and they are different pictures of the same equations.

**The row picture** (what you learned in high school). Each equation $a_1 x + a_2 y = b$ is a line in the plane. Two equations are two lines; the solution is where they cross. Three equations in three unknowns are three planes meeting at a point. This picture is fine in two or three dimensions and useless in ten: nobody can see ten hyperplanes in $\mathbb{R}^{10}$.

**The column picture** (what we will live in). Hermann Grassmann (*Die lineale Ausdehnungslehre*, 1844) made "linear combination of vectors" a basic operation in its own right. Arthur Cayley ("A Memoir on the Theory of Matrices," 1858) treated the array of coefficients as a single object with its own algebra. Put those two moves together and you can write the system

$$
\begin{aligned}
2x - y &= 1 \\
x + y &= 5
\end{aligned}
$$

not as two lines but as one vector equation:

$$
x \begin{pmatrix} 2 \\ 1 \end{pmatrix} + y \begin{pmatrix} -1 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ 5 \end{pmatrix}.
$$

Now the question is: **how much of the first column and how much of the second column do I need to build the right-hand side?** The unknowns $x$ and $y$ are *amounts*. The columns are *ingredients*. The right side is the *recipe target*. This question makes sense in any dimension, because it never asks you to picture anything but "combine these ingredients."

The matrix $A$ is nothing more than the bookkeeping device that holds the columns:

$$
A = \begin{pmatrix} 2 & -1 \\ 1 & 1 \end{pmatrix}, \qquad \mathbf{x} = \begin{pmatrix} x \\ y \end{pmatrix}, \qquad A\mathbf{x} = \mathbf{b}.
$$

That is the whole module in one sentence: $A\mathbf{x}$ is a **linear combination of the columns of $A$ with coefficients from $\mathbf{x}$.** Everything else here is making sure you can compute with it, measure with it (dot products, norms), and ask the right question about it (what can I build? That's *span*).

## 2. Intuition

### Vectors as lists and as arrows

A vector in $\mathbb{R}^2$ is a pair of numbers $\mathbf{v} = (v_1, v_2)^\top$. Draw it as an arrow from the origin to the point $(v_1, v_2)$. Both views are correct; the list view is how you compute, the arrow view is how you think.

```
        ^ y
        |
        |      * (3, 2)
        |     /
        |    /  v = (3, 2)
        |   /
        |  /
        | /
        +-----------------> x
```

Two operations, and only two, define everything:

- **Addition** $\mathbf{u} + \mathbf{v}$: add entrywise. Geometrically, put the tail of $\mathbf{v}$ at the head of $\mathbf{u}$; the sum goes from the origin to the new head. (Parallelogram rule.)
- **Scalar multiplication** $c\mathbf{v}$: multiply every entry by $c$. Geometrically, stretch by $|c|$, and flip if $c < 0$.

A **linear combination** is anything you can build with these two operations: $c_1 \mathbf{v}_1 + c_2 \mathbf{v}_2 + \cdots + c_k \mathbf{v}_k$.

### What can you build? Span

Take one vector $\mathbf{v} = (1, 2)^\top$. All its multiples $c\mathbf{v}$ form a line through the origin. That line is the **span** of $\mathbf{v}$.

Take two vectors $\mathbf{u} = (1, 0)^\top$ and $\mathbf{v} = (1, 2)^\top$. Every point in the plane is $c_1 \mathbf{u} + c_2 \mathbf{v}$ for some $c_1, c_2$. (Convince yourself: to reach $(a, b)$, use $c_2 = b/2$ to get the height right, then $c_1 = a - b/2$ to fix the horizontal.) Their span is all of $\mathbb{R}^2$.

Now take $\mathbf{u} = (1, 2)^\top$ and $\mathbf{w} = (2, 4)^\top$. Since $\mathbf{w} = 2\mathbf{u}$, every combination $c_1 \mathbf{u} + c_2 \mathbf{w} = (c_1 + 2c_2)\mathbf{u}$ is still on the line through $\mathbf{u}$. Two vectors, but their span is only a line. The second vector added nothing new.

This is the whole intuition for the rest of the course: **the span of a set of vectors is the set of everything you can build from them, and the interesting question is always whether a new vector adds a new direction or is already buildable.**

### The two pictures, side by side

For $A\mathbf{x} = \mathbf{b}$ with $A$ being $2 \times 2$:

| | Row picture | Column picture |
|---|---|---|
| What is an equation? | A line (one row of $A$ dotted with $\mathbf{x}$ equals one entry of $\mathbf{b}$) | Not the unit of thought |
| What is the system? | Two lines; find the intersection | One vector equation: combine the columns to hit $\mathbf{b}$ |
| When is there a unique solution? | Lines cross at one point | Columns point in different directions, so they span the plane and there's one recipe for $\mathbf{b}$ |
| When is there no solution? | Parallel lines | Columns are parallel and $\mathbf{b}$ is off their line |
| When are there infinitely many? | Same line twice | Columns are parallel and $\mathbf{b}$ is on their line, so many recipes work |
| Scales to $\mathbb{R}^{100}$? | No | Yes |

### Measuring: the dot product

The dot product $\mathbf{u} \cdot \mathbf{v} = u_1 v_1 + u_2 v_2 + \cdots$ is a single number that answers "how much do $\mathbf{u}$ and $\mathbf{v}$ point the same way?" Positive means roughly aligned, zero means perpendicular, negative means roughly opposite. By the end of §3 you'll have derived $\mathbf{u}\cdot\mathbf{v} = \|\mathbf{u}\|\|\mathbf{v}\|\cos\theta$, which makes this precise.

The guess you should be able to make now: the dot product of a vector with itself is its squared length (Pythagoras), and two vectors are perpendicular exactly when their dot product is zero.

## 3. Formal development

### 3.1 Vectors in $\mathbb{R}^n$

**Definition 3.1.** $\mathbb{R}^n$ is the set of all ordered lists of $n$ real numbers, written as columns. For $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$ and $c \in \mathbb{R}$:

$$
\mathbf{u} + \mathbf{v} = \begin{pmatrix} u_1 + v_1 \\ \vdots \\ u_n + v_n \end{pmatrix}, \qquad c\mathbf{v} = \begin{pmatrix} c v_1 \\ \vdots \\ c v_n \end{pmatrix}.
$$

These satisfy the rules you'd expect: addition is commutative and associative, there is a zero vector $\mathbf{0}$, every vector has a negative, scalar multiplication distributes over both kinds of addition: $c(\mathbf{u} + \mathbf{v}) = c\mathbf{u} + c\mathbf{v}$ and $(c + d)\mathbf{v} = c\mathbf{v} + d\mathbf{v}$, scalars multiply consistently, $(cd)\mathbf{v} = c(d\mathbf{v})$, and $1\mathbf{v} = \mathbf{v}$. In LA-3 these eight rules become the *definition* of a vector space, and you'll see that functions, polynomials, and matrices are vectors too. For now, they're facts about lists of numbers that you can check entry by entry.

**Definition 3.2 (Linear combination, span).** A linear combination of $\mathbf{v}_1, \dots, \mathbf{v}_k$ is any vector of the form $c_1 \mathbf{v}_1 + \cdots + c_k \mathbf{v}_k$ with $c_i \in \mathbb{R}$. The **span** of $\{\mathbf{v}_1, \dots, \mathbf{v}_k\}$ is the set of all such combinations:

$$
\operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_k\} = \{ c_1 \mathbf{v}_1 + \cdots + c_k \mathbf{v}_k : c_1, \dots, c_k \in \mathbb{R} \}.
$$

The span always contains $\mathbf{0}$ (take all $c_i = 0$), and it is closed under addition and scalar multiplication (a combination of combinations is a combination). Those two closure facts are what will make a span a *subspace* in LA-3.

### 3.2 Matrix times vector: the two definitions

Let $A$ be an $m \times n$ matrix ($m$ rows, $n$ columns) with entries $a_{ij}$ (row $i$, column $j$). Write its columns as $\mathbf{a}_1, \dots, \mathbf{a}_n \in \mathbb{R}^m$ and its rows as $\mathbf{r}_1^\top, \dots, \mathbf{r}_m^\top$.

**Definition 3.3 (Column form).** For $\mathbf{x} \in \mathbb{R}^n$,

$$
A\mathbf{x} = x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n. \tag{3.1}
$$

$A\mathbf{x}$ is a linear combination of the columns of $A$, with the entries of $\mathbf{x}$ as the coefficients.

**Definition 3.3′ (Row form).** The $i$-th entry of $A\mathbf{x}$ is the dot product of the $i$-th row of $A$ with $\mathbf{x}$:

$$
(A\mathbf{x})_i = \mathbf{r}_i \cdot \mathbf{x} = \sum_{j=1}^n a_{ij} x_j. \tag{3.2}
$$

▶ **Try it first.** Show that (3.1) and (3.2) give the same vector. Pick the $i$-th entry of the right side of (3.1) and expand.

**Derivation.**

1. The $i$-th entry of $x_j \mathbf{a}_j$ is $x_j a_{ij}$, because the $i$-th entry of column $j$ is $a_{ij}$.
2. The $i$-th entry of a sum of vectors is the sum of the $i$-th entries. So the $i$-th entry of $\sum_j x_j \mathbf{a}_j$ is $\sum_j x_j a_{ij}$.
3. That is exactly $\sum_j a_{ij} x_j = \mathbf{r}_i \cdot \mathbf{x}$, which is (3.2). ∎

■ **What just happened.** Two ways of organizing the same $mn$ multiplications. The row form is how you compute by hand. The column form is how you *think*: $A\mathbf{x}$ lives in the span of the columns of $A$.

Three immediate consequences:

**Theorem 3.4 (Linearity).** For all $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$ and $c \in \mathbb{R}$: $A(\mathbf{x} + \mathbf{y}) = A\mathbf{x} + A\mathbf{y}$ and $A(c\mathbf{x}) = c\,A\mathbf{x}$.

*Proof.* From (3.1): $A(\mathbf{x} + \mathbf{y}) = \sum_j (x_j + y_j)\mathbf{a}_j = \sum_j x_j \mathbf{a}_j + \sum_j y_j \mathbf{a}_j = A\mathbf{x} + A\mathbf{y}$, using the distributive rule for scalars times vectors. Similarly $A(c\mathbf{x}) = \sum_j (c x_j)\mathbf{a}_j = c\sum_j x_j \mathbf{a}_j$. ∎

**Theorem 3.5 (Columns are images of basis vectors).** $A\mathbf{e}_j = \mathbf{a}_j$: the matrix times the $j$-th standard basis vector is the $j$-th column.

*Proof.* In (3.1), $\mathbf{e}_j$ has $x_j = 1$ and all other $x_i = 0$, so the sum collapses to $\mathbf{a}_j$. ∎

This tiny theorem is the key to LA-4: to know what a matrix does to *every* vector, you only need to know what it does to the basis vectors, and those answers are its columns.

**Theorem 3.6 (Solvability of $A\mathbf{x} = \mathbf{b}$).** The system $A\mathbf{x} = \mathbf{b}$ has a solution if and only if $\mathbf{b} \in \operatorname{span}\{\mathbf{a}_1, \dots, \mathbf{a}_n\}$.

*Proof.* By (3.1), $A\mathbf{x} = \mathbf{b}$ says exactly that $\mathbf{b}$ equals some linear combination of the columns, and the solution $\mathbf{x}$ *is* the list of coefficients. ∎

The span of the columns is important enough to get a name in LA-3: the **column space** $C(A)$. Theorem 3.6 says: $A\mathbf{x} = \mathbf{b}$ is solvable exactly when $\mathbf{b}$ is in the column space.

### 3.3 Matrix times matrix

If $B$ is $n \times p$ with columns $\mathbf{b}_1, \dots, \mathbf{b}_p$, define $AB$ to be the $m \times p$ matrix whose $k$-th column is $A\mathbf{b}_k$:

$$
AB = \begin{pmatrix} A\mathbf{b}_1 & A\mathbf{b}_2 & \cdots & A\mathbf{b}_p \end{pmatrix}. \tag{3.3}
$$

So the $(i, k)$ entry of $AB$ is $(A\mathbf{b}_k)_i = \mathbf{r}_i \cdot \mathbf{b}_k = \sum_j a_{ij} b_{jk}$. This is the "row times column" rule, and (3.3) is why it's the right rule: it makes $(AB)\mathbf{x} = A(B\mathbf{x})$.

▶ **Try it first.** Prove $(AB)\mathbf{x} = A(B\mathbf{x})$ for any $\mathbf{x} \in \mathbb{R}^p$, using (3.1) and (3.3) and Theorem 3.4.

**Derivation.**

1. By (3.1) applied to $AB$, whose columns are $A\mathbf{b}_k$: $(AB)\mathbf{x} = \sum_k x_k (A\mathbf{b}_k)$.
2. By linearity (Theorem 3.4), $\sum_k x_k (A\mathbf{b}_k) = A\left(\sum_k x_k \mathbf{b}_k\right)$.
3. By (3.1) applied to $B$: $\sum_k x_k \mathbf{b}_k = B\mathbf{x}$. So $(AB)\mathbf{x} = A(B\mathbf{x})$. ∎

■ **What just happened.** Matrix multiplication is defined so that "do $B$, then do $A$" is a single matrix $AB$. In LA-4 this becomes: *composition of linear maps is matrix multiplication*, and it's the reason the row-times-column rule exists at all. Note the order: $AB$ means $B$ acts first.

Matrix multiplication is **not commutative** in general. Example 4.3 shows $AB \ne BA$.

### 3.4 Transpose

**Definition 3.7.** The transpose $A^\top$ of an $m \times n$ matrix is the $n \times m$ matrix with $(A^\top)_{ij} = a_{ji}$. Rows become columns. For a column vector $\mathbf{v}$, $\mathbf{v}^\top$ is a row.

**Theorem 3.8.** $(AB)^\top = B^\top A^\top$.

*Proof.* $((AB)^\top)_{ki} = (AB)_{ik} = \sum_j a_{ij} b_{jk}$. And $(B^\top A^\top)_{ki} = \sum_j (B^\top)_{kj} (A^\top)_{ji} = \sum_j b_{jk} a_{ij}$. Same sum. ∎

The order reverses. This is the first of many "reverse the order" rules you'll meet (inverses in LA-2, adjoints in LA-8, daggers in QM-4).

### 3.5 The dot product and length

**Definition 3.9.** For $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$, the **dot product** (also *inner product* or *scalar product*) is

$$
\mathbf{u} \cdot \mathbf{v} = \mathbf{u}^\top \mathbf{v} = \sum_{i=1}^n u_i v_i,
$$

and the **Euclidean norm** (length) is $\|\mathbf{v}\| = \sqrt{\mathbf{v} \cdot \mathbf{v}} = \sqrt{v_1^2 + \cdots + v_n^2}$.

Properties, each checkable entry by entry: $\mathbf{u}\cdot\mathbf{v} = \mathbf{v}\cdot\mathbf{u}$ (symmetric); $(c\mathbf{u} + d\mathbf{w})\cdot\mathbf{v} = c(\mathbf{u}\cdot\mathbf{v}) + d(\mathbf{w}\cdot\mathbf{v})$ (linear in each slot); $\mathbf{v}\cdot\mathbf{v} \ge 0$ with equality only for $\mathbf{v} = \mathbf{0}$ (positive definite). In LA-8 these three properties become the *definition* of an inner product on an abstract space.

Why is $\|\mathbf{v}\|$ the length? In $\mathbb{R}^2$, $\sqrt{v_1^2 + v_2^2}$ is the Pythagorean theorem. In $\mathbb{R}^3$, apply Pythagoras twice: first in the $xy$-plane to get $\sqrt{v_1^2 + v_2^2}$, then with the $z$-leg. In $\mathbb{R}^n$ we *define* length this way; it agrees with geometry wherever geometry exists.

**Theorem 3.10 (Angle formula).** For nonzero $\mathbf{u}, \mathbf{v}$ in $\mathbb{R}^2$ or $\mathbb{R}^3$, let $\theta \in [0, \pi]$ be the geometric angle between them. Then

$$
\mathbf{u} \cdot \mathbf{v} = \|\mathbf{u}\| \|\mathbf{v}\| \cos\theta. \tag{3.4}
$$

The hypothesis "$\mathbb{R}^2$ or $\mathbb{R}^3$" matters: those are the spaces where "angle" already has a meaning from geometry, and where we checked above that $\|\cdot\|$ is geometric length. In $\mathbb{R}^n$ for $n \ge 4$ there is no angle to start from, and (3.5) below handles that case.

▶ **Try it first.** The three vectors $\mathbf{u}$, $\mathbf{v}$, and $\mathbf{u} - \mathbf{v}$ form a triangle. Write the law of cosines for it, then expand $\|\mathbf{u} - \mathbf{v}\|^2$ using the dot product.

**Derivation.**

0. *Parallel case.* If $\mathbf{v} = c\mathbf{u}$, the "triangle" is flat, so check (3.4) directly. Here $\mathbf{u}\cdot\mathbf{v} = c\|\mathbf{u}\|^2$ and $\|\mathbf{v}\| = |c|\,\|\mathbf{u}\|$. If $c > 0$ then $\theta = 0$ and the right side is $c\|\mathbf{u}\|^2$. If $c < 0$ then $\theta = \pi$ and the right side is $-|c|\,\|\mathbf{u}\|^2 = c\|\mathbf{u}\|^2$. Either way (3.4) holds. From now on assume $\mathbf{u}, \mathbf{v}$ are not parallel, so the triangle is a genuine one.
1. Law of cosines, with sides $\|\mathbf{u}\|$, $\|\mathbf{v}\|$, and the opposite side $\|\mathbf{u} - \mathbf{v}\|$:
   $$\|\mathbf{u} - \mathbf{v}\|^2 = \|\mathbf{u}\|^2 + \|\mathbf{v}\|^2 - 2\|\mathbf{u}\|\|\mathbf{v}\|\cos\theta.$$
2. Expand the left side with the dot product, using linearity and symmetry:
   $$\|\mathbf{u} - \mathbf{v}\|^2 = (\mathbf{u} - \mathbf{v})\cdot(\mathbf{u} - \mathbf{v}) = \mathbf{u}\cdot\mathbf{u} - 2\,\mathbf{u}\cdot\mathbf{v} + \mathbf{v}\cdot\mathbf{v} = \|\mathbf{u}\|^2 - 2\,\mathbf{u}\cdot\mathbf{v} + \|\mathbf{v}\|^2.$$
3. Set the two expressions equal and cancel $\|\mathbf{u}\|^2 + \|\mathbf{v}\|^2$ from both sides: $-2\,\mathbf{u}\cdot\mathbf{v} = -2\|\mathbf{u}\|\|\mathbf{v}\|\cos\theta$. Divide by $-2$. ∎

■ **What just happened.** The dot product, defined purely algebraically, computes an angle. In particular $\mathbf{u}\cdot\mathbf{v} = 0$ exactly when $\cos\theta = 0$, i.e., the vectors are **perpendicular**. In any $\mathbb{R}^n$ we take $\mathbf{u}\cdot\mathbf{v} = 0$ as the definition of **orthogonal**. This is the start of LA-5.

**Theorem 3.11 (Cauchy–Schwarz inequality).** For all $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$,
$$|\mathbf{u}\cdot\mathbf{v}| \le \|\mathbf{u}\|\|\mathbf{v}\|,$$
with equality if and only if one vector is a scalar multiple of the other.

In $\mathbb{R}^2$ and $\mathbb{R}^3$ this follows from (3.4) because $|\cos\theta| \le 1$. But (3.4) leaned on geometry (the law of cosines), so that argument says nothing in $\mathbb{R}^4$, and in LA-8 we'll want Cauchy–Schwarz in spaces with no geometry at all, like spaces of functions. So here is a proof that uses only the three algebraic properties of the dot product, in every $\mathbb{R}^n$.

▶ **Try it first.** For any real $t$, $\|\mathbf{u} - t\mathbf{v}\|^2 \ge 0$. Expand it as a quadratic in $t$. A quadratic that's never negative has a non-positive discriminant.

**Derivation.**

1. If $\mathbf{v} = \mathbf{0}$ both sides are $0$ and we're done. Assume $\mathbf{v} \ne \mathbf{0}$.
2. For every $t \in \mathbb{R}$:
   $$0 \le \|\mathbf{u} - t\mathbf{v}\|^2 = \|\mathbf{u}\|^2 - 2t\,(\mathbf{u}\cdot\mathbf{v}) + t^2 \|\mathbf{v}\|^2.$$
3. The right side is a quadratic $q(t) = \|\mathbf{v}\|^2 t^2 - 2(\mathbf{u}\cdot\mathbf{v}) t + \|\mathbf{u}\|^2$ with positive leading coefficient that is never negative, so it has at most one real root, so its discriminant is $\le 0$:
   $$4(\mathbf{u}\cdot\mathbf{v})^2 - 4\|\mathbf{v}\|^2\|\mathbf{u}\|^2 \le 0.$$
4. Rearranging: $(\mathbf{u}\cdot\mathbf{v})^2 \le \|\mathbf{u}\|^2\|\mathbf{v}\|^2$. Take square roots.
5. Equality holds iff the discriminant is $0$ iff $q$ has a real root $t_0$ iff $\|\mathbf{u} - t_0\mathbf{v}\| = 0$ iff $\mathbf{u} = t_0\mathbf{v}$. ∎

■ **What just happened.** The quadratic $q$ in step 3 is smallest at its vertex: $q'(t) = 2\|\mathbf{v}\|^2 t - 2\,\mathbf{u}\cdot\mathbf{v} = 0$ gives $t^* = (\mathbf{u}\cdot\mathbf{v})/\|\mathbf{v}\|^2$, and $q(t^*) = \|\mathbf{u}\|^2 - (\mathbf{u}\cdot\mathbf{v})^2/\|\mathbf{v}\|^2 \ge 0$ is Cauchy–Schwarz again. The vector $t^*\mathbf{v}$ is the **projection** of $\mathbf{u}$ onto the line through $\mathbf{v}$: the closest point on that line to $\mathbf{u}$. We just proved Cauchy–Schwarz by projecting. Keep this: it's the first appearance of the formula that runs all of LA-5 and all of least squares.

**The angle in $\mathbb{R}^n$.** Now turn (3.4) around. For nonzero $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$ with any $n$, *define* the angle $\theta \in [0, \pi]$ between them by

$$
\cos\theta = \frac{\mathbf{u}\cdot\mathbf{v}}{\|\mathbf{u}\|\,\|\mathbf{v}\|}. \tag{3.5}
$$

Cauchy–Schwarz is what guarantees the right side lies in $[-1, 1]$, so that $\theta$ exists. In $\mathbb{R}^2$ and $\mathbb{R}^3$, Theorem 3.10 says this is the geometric angle; for $n \ge 4$, (3.5) is the definition. This is the order of logic: algebraic Cauchy–Schwarz first, angle second. Example 4.4 uses (3.5) in $\mathbb{R}^4$.

**Theorem 3.12 (Triangle inequality).** $\|\mathbf{u} + \mathbf{v}\| \le \|\mathbf{u}\| + \|\mathbf{v}\|$.

*Proof.* $\|\mathbf{u} + \mathbf{v}\|^2 = \|\mathbf{u}\|^2 + 2\,\mathbf{u}\cdot\mathbf{v} + \|\mathbf{v}\|^2 \le \|\mathbf{u}\|^2 + 2\|\mathbf{u}\|\|\mathbf{v}\| + \|\mathbf{v}\|^2 = (\|\mathbf{u}\| + \|\mathbf{v}\|)^2$, using Cauchy–Schwarz in the middle. Take square roots. ∎

### 3.6 Other norms (a preview you'll need)

The Euclidean norm isn't the only way to measure length. Two others you'll meet constantly:

- **$\ell^1$ norm:** $\|\mathbf{v}\|_1 = |v_1| + \cdots + |v_n|$ (Manhattan distance). Central to lasso regression in ST-3 and sparsity in ML.
- **$\ell^\infty$ norm:** $\|\mathbf{v}\|_\infty = \max_i |v_i|$. The natural norm for componentwise error bounds in numerical work; the matrix norm it induces (largest absolute row sum) gives the condition number $\kappa_\infty$ in LA-13 Exercise A3.

All three satisfy: $\|\mathbf{v}\| \ge 0$ with equality only at $\mathbf{0}$; $\|c\mathbf{v}\| = |c|\,\|\mathbf{v}\|$; and the triangle inequality. Any function with those three properties is called a **norm**. The first two are immediate from the formulas. The triangle inequality for $\|\cdot\|_2$ is Theorem 3.12. For the other two it comes from the scalar fact $|a + b| \le |a| + |b|$, applied entry by entry:

1. For each $i$, $|u_i + v_i| \le |u_i| + |v_i|$.
2. Sum over $i$: $\|\mathbf{u} + \mathbf{v}\|_1 \le \|\mathbf{u}\|_1 + \|\mathbf{v}\|_1$.
3. Instead bound each right side by its maximum: $|u_i + v_i| \le \|\mathbf{u}\|_\infty + \|\mathbf{v}\|_\infty$ for every $i$, so the largest left side is too: $\|\mathbf{u} + \mathbf{v}\|_\infty \le \|\mathbf{u}\|_\infty + \|\mathbf{v}\|_\infty$. ∎

The unit balls $\{\mathbf{v} : \|\mathbf{v}\| \le 1\}$ are a disk ($\ell^2$), a diamond ($\ell^1$), and a square ($\ell^\infty$) in $\mathbb{R}^2$. The project draws their boundaries, the unit spheres $\{\mathbf{v} : \|\mathbf{v}\| = 1\}$.

### 3.7 Lines and planes as spans plus a shift

The line through $\mathbf{p}$ in direction $\mathbf{d}$ is $\{\mathbf{p} + t\mathbf{d} : t \in \mathbb{R}\}$. The plane through $\mathbf{p}$ spanned by $\mathbf{d}_1, \mathbf{d}_2$ (not parallel) is $\{\mathbf{p} + s\mathbf{d}_1 + t\mathbf{d}_2\}$. Both are "a span, shifted by $\mathbf{p}$."

A plane in $\mathbb{R}^3$ through the origin can also be described by a single equation $\mathbf{n}\cdot\mathbf{x} = 0$, where $\mathbf{n}$ is the normal vector: all vectors perpendicular to $\mathbf{n}$. That's the row picture of a plane (one equation) versus the column picture (span of two directions). LA-3 makes the relationship between these two descriptions precise: it's the relationship between the **null space** and the **row space**.

## 4. Worked examples

### Example 4.1: Both pictures of a $2 \times 2$ system

Solve
$$
\begin{pmatrix} 2 & -1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 1 \\ 5 \end{pmatrix}.
$$

*Row picture.* Line 1: $2x - y = 1$. Line 2: $x + y = 5$. Add them: $3x = 6$, so $x = 2$, then $y = 3$. Intersection at $(2, 3)$.

*Column picture.* We need $x\begin{pmatrix}2\\1\end{pmatrix} + y\begin{pmatrix}-1\\1\end{pmatrix} = \begin{pmatrix}1\\5\end{pmatrix}$. Check $x = 2, y = 3$: $2\begin{pmatrix}2\\1\end{pmatrix} + 3\begin{pmatrix}-1\\1\end{pmatrix} = \begin{pmatrix}4\\2\end{pmatrix} + \begin{pmatrix}-3\\3\end{pmatrix} = \begin{pmatrix}1\\5\end{pmatrix}$. ✓

Two of the first column plus three of the second builds the target. The project plots both pictures.

### Example 4.2: A system with no solution, seen both ways

$$
A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}, \qquad \mathbf{b} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}.
$$

*Row picture.* Lines $x + 2y = 1$ and $2x + 4y = 3$. The second is $x + 2y = 1.5$. Parallel, distinct. No intersection.

*Column picture.* Columns $(1, 2)^\top$ and $(2, 4)^\top$ are parallel (the second is twice the first). Their span is the line $\{t(1, 2)^\top\}$. Is $\mathbf{b} = (1, 3)^\top$ on it? We'd need $t = 1$ from the first entry and $t = 3/2$ from the second. No. So $\mathbf{b} \notin \operatorname{span}$ and by Theorem 3.6 there's no solution.

Change $\mathbf{b}$ to $(1, 2)^\top$: now $\mathbf{b} = 1\cdot(1,2)^\top + 0\cdot(2,4)^\top$, but also $= 0\cdot(1,2)^\top + \tfrac12 (2,4)^\top$, and infinitely many other recipes. Same columns, different right side, infinitely many solutions.

**The lesson:** whether $A\mathbf{x} = \mathbf{b}$ is solvable depends on $\mathbf{b}$ *and* on the columns; whether the solution is *unique* depends only on the columns (are they "independent," LA-3's word for not-parallel-in-any-generalized-sense).

### Example 4.3: Matrix multiplication, order matters, and a mechanics preview

Let $R$ rotate the plane by $90°$ counterclockwise and $S$ stretch the $x$-axis by $2$:
$$
R = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \qquad S = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}.
$$
(Check with Theorem 3.5: $R\mathbf{e}_1 = (0, 1)^\top$, which is $\mathbf{e}_1$ rotated by $90°$. $R\mathbf{e}_2 = (-1, 0)^\top$. ✓)

Compute both products using (3.3), column by column:

$RS$: columns are $R(2, 0)^\top = (0, 2)^\top$ and $R(0, 1)^\top = (-1, 0)^\top$. So $RS = \begin{pmatrix} 0 & -1 \\ 2 & 0 \end{pmatrix}$.

$SR$: columns are $S(0, 1)^\top = (0, 1)^\top$ and $S(-1, 0)^\top = (-2, 0)^\top$. So $SR = \begin{pmatrix} 0 & -2 \\ 1 & 0 \end{pmatrix}$.

$RS \ne SR$. Geometrically: $RS$ means stretch first, then rotate; $SR$ means rotate first, then stretch. Stretching along $x$ and then rotating is not the same as rotating and then stretching along $x$. Apply both to $\mathbf{e}_1$: $RS\mathbf{e}_1 = (0, 2)^\top$ but $SR\mathbf{e}_1 = (0, 1)^\top$.

*Mechanics preview.* In CM-7, the orientation of a rigid body is a rotation matrix, and the fact that rotations about different axes don't commute is why Euler angles are such a pain and why angular velocity is a vector but finite rotation is not.

### Example 4.4: Dot products in $\mathbb{R}^4$, where there's nothing to draw

Let $\mathbf{u} = (1, 2, 0, -1)^\top$ and $\mathbf{v} = (3, -1, 2, 1)^\top$.

$\mathbf{u}\cdot\mathbf{v} = 3 - 2 + 0 - 1 = 0$. These vectors are orthogonal, in a space you cannot picture, and you know it with certainty.

$\|\mathbf{u}\| = \sqrt{1 + 4 + 0 + 1} = \sqrt{6}$, $\|\mathbf{v}\| = \sqrt{9 + 1 + 4 + 1} = \sqrt{15}$.

Cauchy–Schwarz check: $|\mathbf{u}\cdot\mathbf{v}| = 0 \le \sqrt{6}\sqrt{15} = \sqrt{90} \approx 9.49$. ✓ (Trivially, but the project checks thousands of random pairs.)

Now $\mathbf{w} = (1, 1, 1, 1)^\top$: $\mathbf{u}\cdot\mathbf{w} = 2$, $\|\mathbf{w}\| = 2$, so $\cos\theta = 2/(\sqrt{6}\cdot 2) = 1/\sqrt{6} \approx 0.408$, $\theta \approx 65.9°$. The "angle between two 4-vectors" is a real, computable thing.

### Example 4.5: Projection, the formula that keeps coming back

Project $\mathbf{u} = (3, 1)^\top$ onto the line through $\mathbf{v} = (1, 1)^\top$.

From the Cauchy–Schwarz proof, the closest point on the line is $t^*\mathbf{v}$ with $t^* = (\mathbf{u}\cdot\mathbf{v})/\|\mathbf{v}\|^2 = 4/2 = 2$. So the projection is $\mathbf{p} = 2(1, 1)^\top = (2, 2)^\top$.

Check: the error $\mathbf{u} - \mathbf{p} = (1, -1)^\top$, and $(1, -1)\cdot(1, 1) = 0$. The error is perpendicular to the line. That perpendicularity is the defining property of a projection and the whole content of LA-5 and of least squares in ST-3.

## 5. Why it matters and where it connects

- **LA-3** names $\operatorname{span}$ of the columns the *column space* and asks "which columns are redundant?" (independence, basis, dimension). Theorem 3.6 becomes the first of four fundamental subspaces.
- **LA-4** takes Theorem 3.5 ($A\mathbf{e}_j = \mathbf{a}_j$) and makes it the bridge between abstract linear maps and matrices: a matrix is what a linear map does to a basis.
- **LA-5** takes Example 4.5 and generalizes it: projection onto a subspace, least squares, orthogonal bases. The $t^* = \mathbf{u}\cdot\mathbf{v}/\|\mathbf{v}\|^2$ formula becomes $P = A(A^\top A)^{-1}A^\top$.
- **LA-8** takes the three properties of the dot product (symmetric, linear, positive definite) as the definition of an inner product, so that functions can have angles and Fourier series (F2) become projections.
- **CM-7**: rotation matrices and non-commutativity (Example 4.3).
- **QM-4**: the dot product becomes $\langle\phi|\psi\rangle$, probabilities are $|\langle\phi|\psi\rangle|^2$, and Cauchy–Schwarz becomes the uncertainty principle. Literally: the uncertainty principle is Cauchy–Schwarz applied to two particular vectors.
- **ML-1, ML-4**: cosine similarity (Theorem 3.10) is how embeddings are compared; the kernel trick is "replace the dot product."
- **ST-3**: regression is projection (Example 4.5 in $n$ dimensions).

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills on computing $A\mathbf{x}$ both ways and dot products. B1–B7 are the module: B2 (prove linearity from the row form), B4 (span of three vectors in $\mathbb{R}^3$), B6 (parallelogram law), and B7 (project.py extension) are essential. C1–C3 go further: Cauchy–Schwarz for weighted dot products, the $\ell^1$ triangle inequality, and a first look at matrix norms.

## 7. Coding project

`project.py` builds the geometric intuition this module is about:

1. **Two pictures.** Plots the row picture (lines) and column picture (vectors and their combination) for Example 4.1, then for a random $2 \times 2$ system. Verifies that the solution found by NumPy is the right recipe for the columns.
2. **Span explorer.** Samples thousands of random coefficients for two or three vectors in $\mathbb{R}^3$ and plots where the combinations land, for independent vectors (fills space), two vectors (a plane), and a dependent triple (a plane, not space). You will see what "span" looks like.
3. **$A\mathbf{x}$ two ways.** Implements `matvec_rows` and `matvec_columns` from (3.1) and (3.2) and asserts they agree with each other and with NumPy on random matrices.
4. **Cauchy–Schwarz and the triangle inequality** are checked on 10,000 random pairs. The angle formula is checked against `np.arccos`.
5. **Unit balls** for $\ell^1$, $\ell^2$, $\ell^\infty$ are drawn.

**Build it yourself:** before reading the code, write `matvec_columns(A, x)` (loop over columns, accumulate `x[j] * A[:, j]`) and `project_onto_line(u, v)`. Then compare.

## 8. Where you'll struggle

- **Believing the column picture is worth it.** It feels like a notational trick when the row picture gives you the answer just as fast for $2\times 2$. It isn't a trick. The row picture does not survive past $\mathbb{R}^3$, and almost every theorem from LA-3 onward is a statement about columns. If you find yourself defaulting to "rows dotted with $\mathbf{x}$" when thinking (not computing), stop and rewrite as a combination of columns.
- **$AB$ versus $BA$.** You will compute the wrong one at least once. Say it out loud: "$AB$ applies $B$ first." Check with Example 4.3 until it's automatic.
- **Shape errors.** An $m \times n$ matrix times an $n \times p$ matrix is $m \times p$. The inner dimensions must match and disappear. Write shapes above every product until you stop making mistakes. (This habit saves hours in ML-9.)
- **The Cauchy–Schwarz proof** is the first "trick" proof in the course: introduce a free variable $t$, use a positivity fact, extract a conclusion from the discriminant. It feels like magic the first time. It's a standard move you'll see a dozen more times (LA-8, F4's variance inequality, QM-4's uncertainty principle). Learn the move, not just the result.
- **Thinking $\mathbf{u}\cdot\mathbf{v} = 0$ means one of them is zero.** It doesn't. Example 4.4.

## 9. Check yourself

1. Write $A\mathbf{x}$ for a $3\times 2$ matrix $A$ as a linear combination of its columns. Which vector supplies the coefficients?
2. What is $A\mathbf{e}_2$? Why?
3. State in one sentence when $A\mathbf{x} = \mathbf{b}$ has a solution, using the word "span."
4. If the columns of a $2 \times 2$ matrix are parallel, how many solutions can $A\mathbf{x} = \mathbf{b}$ have? (Two possible answers; give both and say when each happens.)
5. Why is matrix multiplication defined by "row times column"? (Answer with (3.3) and the composition property.)
6. What is $(AB)^\top$?
7. Derive $\mathbf{u}\cdot\mathbf{v} = \|\mathbf{u}\|\|\mathbf{v}\|\cos\theta$ from the law of cosines in three lines.
8. Prove Cauchy–Schwarz without geometry. What is the geometric meaning of the minimizing $t$?
9. Give two vectors in $\mathbb{R}^5$, neither zero, with dot product zero.
10. What are the three properties that make a function a norm? Sketch the $\ell^1$ unit ball in $\mathbb{R}^2$.
11. What is the projection of $\mathbf{u}$ onto the line through $\mathbf{v}$? What is perpendicular to what?

## 10. Pair with

- **Strang, *Introduction to Linear Algebra*, ch. 1** (§1.1–1.3). Strang's exposition of the two pictures is the original and still the best; read it for the pictures of planes in $\mathbb{R}^3$.
- **MIT 18.06, lecture 1 (Strang)** on OCW. Watch him draw the row and column pictures live. Twenty minutes well spent.
- **3Blue1Brown, *Essence of Linear Algebra*, episodes 1–4.** Animated vectors, span, linear transformations as grid deformations, and matrix multiplication as composition. Watch episodes 3 and 4 after §3.3; they make (3.3) visual.
- **Axler, *Linear Algebra Done Right*, §1A.** Only if you want to see $\mathbb{R}^n$ and $\mathbb{C}^n$ defined with full rigor. Skip for now otherwise; we return to Axler in LA-3.
