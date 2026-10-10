# LA-4: Linear Maps: Matrices as Coordinates of a Map, Change of Basis, Kernel and Image, Rank–Nullity, Isomorphism

**Time:** ~12 hours · **Prerequisites:** LA-3 · **Unlocks:** LA-5, LA-6

> **Conventions in this module.** $V$ and $W$ are finite-dimensional real vector spaces (LA-3) with $\dim V = n$ and $\dim W = m$. A linear map is written $T: V \to W$; $T\mathbf{v}$ and $T(\mathbf{v})$ mean the same thing. A **basis** is an *ordered* list $\mathcal{B} = (\mathbf{b}_1, \dots, \mathbf{b}_n)$; reordering it changes every matrix below. The coordinate vector of $\mathbf{v}$ in $\mathcal{B}$ is $[\mathbf{v}]_{\mathcal{B}} \in \mathbb{R}^n$. The matrix of $T$ with respect to the input basis $\mathcal{B}$ and output basis $\mathcal{C}$ is $[T]_{\mathcal{B},\mathcal{C}}$, an $m \times n$ matrix whose $j$-th column is $[T\mathbf{b}_j]_{\mathcal{C}}$; when $V = W$ and $\mathcal{B} = \mathcal{C}$ we write $[T]_{\mathcal{B}}$. **Change of basis:** if $\mathcal{B}$ is the old basis and $\mathcal{B}'$ the new one, $P$ is the matrix whose columns are the *new* basis vectors written in *old* coordinates, so $P[\mathbf{v}]_{\mathcal{B}'} = [\mathbf{v}]_{\mathcal{B}}$ (new coordinates in, old coordinates out) and $[T]_{\mathcal{B}'} = P^{-1}[T]_{\mathcal{B}}P$. $\mathcal{P}_n$ is the space of polynomials of degree at most $n$ (dimension $n + 1$) and $D$ is differentiation. $\ker T$ is the kernel (null space), $\operatorname{im} T$ the image (range). Composition $S \circ T$ means "$T$ first, then $S$." Rotation angles are counterclockwise positive. $\operatorname{tr} A$ is the sum of the diagonal entries.

## 0. Where this sits

From LA-3 you need: vector space, subspace, linear independence, basis, dimension, the fact that every basis of $V$ has the same length, and the **basis extension theorem** (every linearly independent list in a finite-dimensional space extends to a basis). From LA-1 you need Theorem 3.5 ($A\mathbf{e}_j = \mathbf{a}_j$) and the column-by-column definition of $AB$, LA-1 eq. (3.3). From LA-2 you need pivots, free variables, and the inverse.

Until now a matrix has been the object and "the map $\mathbf{x} \mapsto A\mathbf{x}$" a thing you say about it. This module reverses that. The **map** is the object. A matrix is a *description* of the map in a chosen basis, the way $(3, 2)$ is a description of an arrow once you have chosen axes. Rotate the axes and the numbers change; the arrow doesn't. Rotate the basis and the matrix changes; the map doesn't.

Everything from here on uses this. LA-5 studies projections, which are linear maps with $P^2 = P$, and orthogonal changes of basis. LA-6 defines the determinant as the volume-scaling factor of a map and proves it is basis-independent. LA-7 is entirely about finding the basis in which an operator's matrix is diagonal, and "diagonalization" is the change of basis $P^{-1}AP$ of this module with a well-chosen $P$. CM-7's inertia tensor and principal axes are Example 4.4 in three dimensions.

## 1. The problem that forced this idea

Take the matrix $\operatorname{diag}(2, 1) = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$. What does it do? "Stretches by 2 along the first basis direction." If the first basis vector is $\mathbf{e}_1$, that is a horizontal stretch. If someone else's first basis vector is $(1, 1)/\sqrt{2}$, the same matrix describes a stretch along the diagonal. The matrix alone does not tell you the map. The matrix plus a basis does.

The reverse is equally unsettling. Rotation of the plane by $90°$ has matrix $\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$ in the standard basis. Describe the same rotation using the basis $(1, 1), (1, -1)$ and you get $\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ (Example 4.2). Same motion, different numbers, and in this case the numbers *look like* rotation the other way.

Physicists hit this wall first, and hit it hard, because the laws they wrote down are about objects that exist independently of anyone's axes. By 1765 Euler had shown that a rigid body has three perpendicular **principal axes**: a frame in which the relationship between angular velocity and angular momentum, $\mathbf{L} = I\boldsymbol{\omega}$, decouples into three separate equations, because the matrix $I$ is diagonal there. In any other frame the same physical relationship is described by a full $3 \times 3$ array with off-diagonal entries, and the equations of motion look far worse. The body hasn't changed. Cauchy, in the 1820s, found the same structure in continuum mechanics: at each point of a loaded solid there is a **stress tensor** $\sigma$, a linear map that sends the unit normal of any small cut surface to the force per area acting across it, and there is a frame (the principal stress directions) in which the shear components vanish. Otto Mohr's circle (1882), which every mechanical engineer learns, is a graphical way of reading off how the components of $\sigma$ change as you rotate the frame. Example 4.4 shows it is nothing but the change-of-basis formula of this module.

So the problem was: **separate the thing from its description.** Say what a linear map *is* without reference to a basis, say how its matrix transforms when the basis changes, and find the quantities that do not change.

The algebra caught up in the second half of the nineteenth century. Arthur Cayley's 1858 *Memoir on the Theory of Matrices* is usually cited as the birth of matrix algebra, and his starting point was not arrays of numbers but *linear substitutions*: rules like $x' = ax + by,\ y' = cx + dy$ for changing variables. Doing one substitution after another gives another substitution, and Cayley observed that the rule for the new coefficients is exactly row-times-column. The matrix was the bookkeeping; the substitution was the thing. Twenty years later Georg Frobenius, in a long 1878 paper on linear substitutions and bilinear forms, asked the question this module is built around: when you change variables, which features of the array survive? The notion of the **rank** of a matrix is usually credited to him, as is the systematic study of when two matrices are "the same substitution in different variables." The fully abstract point of view, in which a linear map between vector spaces is defined first and matrices come second, arrived with Peano's 1888 axioms for vector spaces and became standard only in the twentieth century. It is the point of view of Axler's book, and of this module.

## 2. Intuition

### A linear map is a grid deformation

Draw the usual square grid on the plane. Apply a linear map $T$ and watch where every grid point goes. Three things always happen and nothing else can:

```
   before (standard grid)              after a shear T(x, y) = (x + y, y)

   +---+---+---+---+                       /---/---/---/---/
   |   |   |   |   |                      /   /   /   /   /
   +---+---+---+---+                     /---/---/---/---/
   |   |   |   |   |                    /   /   /   /   /
   +---+---+---+---+                   /---/---/---/---/
   |   |   |   |   |                  /   /   /   /   /
   O---+---+---+---+                 O---/---/---/---/
```

1. Grid lines stay straight.
2. Parallel, evenly spaced lines stay parallel and evenly spaced.
3. The origin stays put.

All three are the single algebraic condition $T(x\mathbf{e}_1 + y\mathbf{e}_2) = x\,T\mathbf{e}_1 + y\,T\mathbf{e}_2$. The grid lines are the sets "$x$ = constant" and "$y$ = constant," and after the map they are "multiples of $T\mathbf{e}_1$ plus a fixed amount of $T\mathbf{e}_2$," which are again straight, parallel, evenly spaced lines. The project draws this for a shear, a rotation, a projection, and a random map (`figures/grid_deformations.png`).

The consequence that drives the whole module: **once you know where the two basis arrows $\mathbf{e}_1, \mathbf{e}_2$ land, you know where every point lands**, because every point is $x\mathbf{e}_1 + y\mathbf{e}_2$ and linearity does the rest. Two arrows in $\mathbb{R}^2$ are four numbers. Arrange them as columns and you have a $2 \times 2$ matrix. LA-1 Theorem 3.5 said "the columns of $A$ are $A\mathbf{e}_j$"; now read it backwards: *a matrix is a record of where a basis goes.* That's all a matrix ever is.

### What gets squashed, where you can land

Project the plane onto the line $y = 2x$. Every point slides perpendicularly onto the line. Two questions:

- **What gets squashed to zero?** Everything on the perpendicular line through the origin (direction $(2, -1)$) lands on the origin. That line is the **kernel**. It is a subspace: one dimension of the plane is destroyed.
- **Where can you land?** Only on the line $y = 2x$. That line is the **image**. One dimension survives.

Dimension bookkeeping: $2 = 1 + 1$. Started with two, destroyed one, kept one.

Try it on a rotation: nothing is squashed (the kernel is $\{\mathbf{0}\}$ alone), everything is reachable (image is the whole plane). $2 = 0 + 2$. A map from $\mathbb{R}^3$ to $\mathbb{R}^2$, like $(x, y, z) \mapsto (x, y)$: the $z$-axis is squashed, the image is all of $\mathbb{R}^2$. $3 = 1 + 2$. A map from $\mathbb{R}^2$ into $\mathbb{R}^3$, like $(x, y) \mapsto (x, y, 0)$: nothing squashed, the image is a plane inside $\mathbb{R}^3$. $2 = 0 + 2$.

The guess you should now be able to make: **(dimensions squashed) + (dimensions kept) = (dimensions you started with).** That is Rank–Nullity, Theorem 3.14. Note that the right-hand side is the dimension of the *input* space, not the output. The output space only constrains how big the image can be.

A second guess: if the input and output spaces have the same dimension, then "nothing squashed" and "everything reachable" must happen together, because $n = 0 + (\text{kept})$ forces kept $= n$. That is Corollary 3.15.

### Same map, different grid

Go back to rotation by $90°$. In the standard grid: $\mathbf{e}_1 \mapsto \mathbf{e}_2$ and $\mathbf{e}_2 \mapsto -\mathbf{e}_1$. Columns $(0, 1)$ and $(-1, 0)$.

Now describe the same rotation using the grid built from $\mathbf{b}_1 = (1, 1)$ and $\mathbf{b}_2 = (1, -1)$. Rotating $\mathbf{b}_1 = (1, 1)$ by $90°$ gives $(-1, 1)$, which is $-\mathbf{b}_2$. Rotating $\mathbf{b}_2 = (1, -1)$ gives $(1, 1) = \mathbf{b}_1$. So in this grid the record is "$\mathbf{b}_1 \mapsto 0\cdot\mathbf{b}_1 - 1\cdot\mathbf{b}_2$, $\mathbf{b}_2 \mapsto 1\cdot\mathbf{b}_1 + 0\cdot\mathbf{b}_2$": columns $(0, -1)$ and $(1, 0)$. Different matrix. Same rotation.

Euler's principal axes are the same idea: the inertia tensor is one map, and there is one grid in which its matrix is diagonal. What you should expect from §3: a mechanical recipe for converting the matrix in one grid to the matrix in another, and it should read as "convert coordinates into the old grid, apply the old matrix, convert the answer back." Three steps, three matrices, multiplied. That is $P^{-1}AP$.

## 3. Formal development

### 3.1 Linear maps and examples

**Definition 3.1 (Linear map).** Let $V$ and $W$ be vector spaces. A function $T: V \to W$ is **linear** if for all $\mathbf{u}, \mathbf{v} \in V$ and all $a, b \in \mathbb{R}$,

$$
T(a\mathbf{u} + b\mathbf{v}) = a\,T\mathbf{u} + b\,T\mathbf{v}. \tag{3.1}
$$

Equivalently, $T$ respects both operations: $T(\mathbf{u} + \mathbf{v}) = T\mathbf{u} + T\mathbf{v}$ and $T(c\mathbf{v}) = c\,T\mathbf{v}$. (Take $a = b = 1$, and $b = 0$.) When $V = W$ we call $T$ an **operator** on $V$. The set of all linear maps from $V$ to $W$ is written $\mathcal{L}(V, W)$.

Two immediate consequences. First, $T\mathbf{0} = T(0\cdot\mathbf{0}) = 0\cdot T\mathbf{0} = \mathbf{0}$: a linear map always sends zero to zero. Second, by induction on (3.1), $T\left(\sum_i c_i\mathbf{v}_i\right) = \sum_i c_i\,T\mathbf{v}_i$ for any finite combination.

Examples, each with the one-line reason it is linear:

- **(a) Matrix multiplication.** $T\mathbf{x} = A\mathbf{x}$ for a fixed $m \times n$ matrix $A$, from $\mathbb{R}^n$ to $\mathbb{R}^m$. Linear by LA-1 Theorem 3.4.
- **(b) Rotation of the plane** by angle $\theta$. Linear because rotating a parallelogram rotates its diagonal, and rotating a scaled arrow scales the rotated arrow. Its matrix is in §3.13.
- **(c) Projection onto a line** through a unit vector $\mathbf{u}$: $T\mathbf{v} = (\mathbf{u}\cdot\mathbf{v})\,\mathbf{u}$ (LA-1 Example 4.5). Linear because the dot product is linear in $\mathbf{v}$.
- **(d) Differentiation** $D: \mathcal{P}_n \to \mathcal{P}_n$, $Dp = p'$. Linear because $(af + bg)' = af' + bg'$.
- **(e) Integration** $J: \mathcal{P}_n \to \mathcal{P}_{n+1}$, $(Jp)(x) = \int_0^x p(t)\,dt$. Linear because integrals are.
- **(f) Evaluation at a point** $\operatorname{ev}_a: \mathcal{P}_n \to \mathbb{R}$, $p \mapsto p(a)$. Linear because $(af + bg)(a) = af(a) + bg(a)$.
- **(g) Transpose** $A \mapsto A^\top$ from $\mathbb{R}^{m\times n}$ to $\mathbb{R}^{n\times m}$. Linear because $(aA + bB)^\top = aA^\top + bB^\top$ entry by entry.
- **(h) Matrix multiplication as a map on matrices.** For fixed $M$, $L_M: X \mapsto MX$ on $\mathbb{R}^{n\times n}$. Linear because $M(aX + bY) = aMX + bMY$ (LA-1 Theorem 3.4 applied column by column). So is $X \mapsto MX - XM$, the commutator with $M$, which runs all of QM-4.
- **(i) Trace** $A \mapsto \operatorname{tr} A$ from $\mathbb{R}^{n\times n}$ to $\mathbb{R}$. Linear because the diagonal of $aA + bB$ is $a$ times the diagonal of $A$ plus $b$ times the diagonal of $B$.

Non-examples, each with the property that fails: translation $\mathbf{v} \mapsto \mathbf{v} + \mathbf{c}$ ($T\mathbf{0} \neq \mathbf{0}$ unless $\mathbf{c} = \mathbf{0}$); the norm $\mathbf{v} \mapsto \|\mathbf{v}\|$ ($\|-\mathbf{v}\| \neq -\|\mathbf{v}\|$); the determinant on $2\times 2$ matrices ($\det(2A) = 4\det A$); squaring a polynomial.

### 3.2 A linear map is determined by what it does to a basis

**Theorem 3.2 (Existence and uniqueness of the linear extension).** Let $(\mathbf{b}_1, \dots, \mathbf{b}_n)$ be a basis of $V$ and let $\mathbf{w}_1, \dots, \mathbf{w}_n$ be *any* vectors in $W$ (repeats and zeros allowed). Then there is exactly one linear map $T: V \to W$ with $T\mathbf{b}_j = \mathbf{w}_j$ for every $j$.

▶ **Try it first.** If such a $T$ exists, what must it do to $\mathbf{v} = \sum_j x_j\mathbf{b}_j$? Write that down as the definition, then ask why the definition is legal: could $\mathbf{v}$ have two different expressions in the basis?

**Derivation.**

1. *Coordinates are unique.* Every $\mathbf{v} \in V$ is $\sum_j x_j \mathbf{b}_j$ for some scalars because a basis spans. If also $\mathbf{v} = \sum_j y_j\mathbf{b}_j$, subtract: $\sum_j (x_j - y_j)\mathbf{b}_j = \mathbf{0}$, and independence forces every $x_j - y_j = 0$. So the scalars $x_j$ are determined by $\mathbf{v}$.
2. *Define* $T\mathbf{v} := \sum_j x_j \mathbf{w}_j$, where $x_j$ are the coordinates from step 1. Step 1 is what makes this a well-defined function: one input, one output.
3. *Linearity.* Let $\mathbf{v} = \sum_j x_j\mathbf{b}_j$ and $\mathbf{u} = \sum_j y_j\mathbf{b}_j$. Then $a\mathbf{v} + b\mathbf{u} = \sum_j (ax_j + by_j)\mathbf{b}_j$, and by step 1 these are *the* coordinates of $a\mathbf{v} + b\mathbf{u}$. So
   $$T(a\mathbf{v} + b\mathbf{u}) = \sum_j (ax_j + by_j)\mathbf{w}_j = a\sum_j x_j\mathbf{w}_j + b\sum_j y_j\mathbf{w}_j = a\,T\mathbf{v} + b\,T\mathbf{u}.$$
4. *It hits the targets.* $\mathbf{b}_j$ has coordinates $(0, \dots, 1, \dots, 0)$ with the $1$ in slot $j$, so $T\mathbf{b}_j = \mathbf{w}_j$.
5. *Uniqueness.* If $T'$ is linear with $T'\mathbf{b}_j = \mathbf{w}_j$ for all $j$, then for any $\mathbf{v} = \sum_j x_j\mathbf{b}_j$, linearity gives $T'\mathbf{v} = \sum_j x_j T'\mathbf{b}_j = \sum_j x_j \mathbf{w}_j = T\mathbf{v}$. So $T' = T$. ∎

■ **What just happened.** Linearity is a rigidity condition. A general function $V \to W$ needs infinitely many values specified; a linear one needs $n$. Choosing $n$ vectors in an $m$-dimensional space is $mn$ numbers. That is the size of an $m \times n$ matrix, and Theorem 3.19 will make the match exact.

### 3.3 Coordinate vectors

**Definition 3.3 (Coordinate vector).** For a basis $\mathcal{B} = (\mathbf{b}_1, \dots, \mathbf{b}_n)$ of $V$ and $\mathbf{v} = \sum_j x_j\mathbf{b}_j$, the **coordinate vector** of $\mathbf{v}$ is $[\mathbf{v}]_{\mathcal{B}} = (x_1, \dots, x_n)^\top \in \mathbb{R}^n$.

Examples. In $\mathbb{R}^n$ with the standard basis, $[\mathbf{v}]_{\text{std}} = \mathbf{v}$: the vector is its own coordinate vector, which is why you have never had to think about this. In $\mathcal{P}_3$ with the monomial basis $(1, x, x^2, x^3)$, $[\,2 - x + 5x^3\,] = (2, -1, 0, 5)^\top$. In the basis $(1, x, x^2/2, x^3/6)$ the same polynomial is $2\cdot 1 - 1\cdot x + 0\cdot\tfrac{x^2}{2} + 30\cdot\tfrac{x^3}{6}$, so its coordinate vector is $(2, -1, 0, 30)^\top$. In $\mathbb{R}^{2\times 2}$ with the basis $(E_{11}, E_{12}, E_{21}, E_{22})$ (the matrix units: a single $1$ in the named slot), $\left[\begin{pmatrix} a & b \\ c & d\end{pmatrix}\right] = (a, b, c, d)^\top$: a matrix "flattened row by row."

**Theorem 3.4 (The coordinate map is an isomorphism).** The map $\kappa_{\mathcal{B}}: V \to \mathbb{R}^n$, $\mathbf{v} \mapsto [\mathbf{v}]_{\mathcal{B}}$, is linear, injective, and surjective.

*Proof.* Linear: step 3 of the previous derivation showed that the coordinates of $a\mathbf{v} + b\mathbf{u}$ are $ax_j + by_j$, i.e. $[a\mathbf{v} + b\mathbf{u}]_{\mathcal{B}} = a[\mathbf{v}]_{\mathcal{B}} + b[\mathbf{u}]_{\mathcal{B}}$. Injective: if $[\mathbf{v}]_{\mathcal{B}} = [\mathbf{u}]_{\mathcal{B}}$ then $\mathbf{v}$ and $\mathbf{u}$ are the same combination of the $\mathbf{b}_j$, so $\mathbf{v} = \mathbf{u}$. Surjective: given $\mathbf{x} \in \mathbb{R}^n$, the vector $\mathbf{v} = \sum_j x_j\mathbf{b}_j$ has $[\mathbf{v}]_{\mathcal{B}} = \mathbf{x}$. ∎

The word *isomorphism* is defined properly in §3.10. For now: $\kappa_{\mathcal{B}}$ is a perfect dictionary between $V$ and $\mathbb{R}^n$ that respects addition and scaling. Anything you can do in $\mathbb{R}^n$ with lists of numbers, you can do in $V$ by translating back and forth. This is the license for the project to store polynomials and $2 \times 2$ matrices as flat NumPy arrays.

### 3.4 The matrix of a linear map

**Definition 3.5 (Matrix of a map).** Let $T: V \to W$ be linear, $\mathcal{B} = (\mathbf{b}_1, \dots, \mathbf{b}_n)$ a basis of $V$ and $\mathcal{C} = (\mathbf{c}_1, \dots, \mathbf{c}_m)$ a basis of $W$. The **matrix of $T$ with respect to $\mathcal{B}$ and $\mathcal{C}$** is the $m \times n$ matrix

$$
[T]_{\mathcal{B},\mathcal{C}} = \Big( [T\mathbf{b}_1]_{\mathcal{C}} \;\; [T\mathbf{b}_2]_{\mathcal{C}} \;\; \cdots \;\; [T\mathbf{b}_n]_{\mathcal{C}} \Big): \qquad \text{column } j = [T\mathbf{b}_j]_{\mathcal{C}}. \tag{3.2}
$$

In words: send each input basis vector through $T$, write the result in output coordinates, make that column $j$. The entry $M_{ij}$ is the amount of $\mathbf{c}_i$ in $T\mathbf{b}_j$: $T\mathbf{b}_j = \sum_i M_{ij}\mathbf{c}_i$. Note the shape: $m$ rows (one per output basis vector), $n$ columns (one per input basis vector).

This is LA-1 Theorem 3.5 turned into a definition. There, for $T\mathbf{x} = A\mathbf{x}$ on $\mathbb{R}^n$ with standard bases, the columns of $A$ were $A\mathbf{e}_j = T\mathbf{e}_j = [T\mathbf{e}_j]_{\text{std}}$. So $[T]_{\text{std},\text{std}} = A$: a matrix is its own matrix in the standard bases. Definition 3.5 extends the idea to any map and any bases.

**Theorem 3.6 (The matrix computes the map).** For every $\mathbf{v} \in V$,

$$
[T\mathbf{v}]_{\mathcal{C}} = [T]_{\mathcal{B},\mathcal{C}}\,[\mathbf{v}]_{\mathcal{B}}. \tag{3.3}
$$

▶ **Try it first.** Write $\mathbf{v} = \sum_j x_j\mathbf{b}_j$, apply $T$, take $\mathcal{C}$-coordinates of both sides (they are linear by Theorem 3.4), and recognize LA-1 eq. (3.1).

**Derivation.** Write $M = [T]_{\mathcal{B},\mathcal{C}}$ and $\mathbf{x} = [\mathbf{v}]_{\mathcal{B}}$.

1. $\mathbf{v} = \sum_j x_j\mathbf{b}_j$ by definition of coordinates.
2. $T\mathbf{v} = \sum_j x_j\,T\mathbf{b}_j$ by linearity of $T$.
3. $[T\mathbf{v}]_{\mathcal{C}} = \sum_j x_j\,[T\mathbf{b}_j]_{\mathcal{C}}$ by linearity of $\kappa_{\mathcal{C}}$ (Theorem 3.4).
4. $[T\mathbf{b}_j]_{\mathcal{C}}$ is column $j$ of $M$, so the right side is $\sum_j x_j\,(\text{column } j \text{ of } M)$, which is $M\mathbf{x}$ by the column form of matrix–vector multiplication, LA-1 eq. (3.1). ∎

■ **What just happened.** Here is the picture to hold on to:

```
            T
      V  ───────►  W
      │            │
  [·]_B│            │[·]_C
      ▼            ▼
     R^n ───────► R^m
          [T]_{B,C}
```

Go down then across, or across then down: same answer. An abstract map on an abstract space is computed by an ordinary matrix–vector product on coordinate vectors. Every linear map between finite-dimensional spaces is "secretly" multiplication by a matrix, once bases are fixed. And conversely, by Theorem 3.2, every $m \times n$ matrix $M$ defines exactly one linear map (send $\mathbf{b}_j$ to $\sum_i M_{ij}\mathbf{c}_i$). The dictionary is perfect in both directions; Theorem 3.19 makes it an isomorphism of vector spaces.

### 3.5 Composition is matrix multiplication

**Theorem 3.7 (Composition).** Let $T: V \to W$ and $S: W \to U$ be linear, with bases $\mathcal{B}$ of $V$, $\mathcal{C}$ of $W$, $\mathcal{D}$ of $U$. Then

$$
[S \circ T]_{\mathcal{B},\mathcal{D}} = [S]_{\mathcal{C},\mathcal{D}}\,[T]_{\mathcal{B},\mathcal{C}}. \tag{3.4}
$$

▶ **Try it first.** Compute column $j$ of the left side by applying (3.3) to $S$ with input $T\mathbf{b}_j$. Then use LA-1 eq. (3.3), the column-by-column definition of a matrix product.

**Derivation.**

1. By (3.2), column $j$ of $[S\circ T]_{\mathcal{B},\mathcal{D}}$ is $[(S\circ T)\mathbf{b}_j]_{\mathcal{D}} = [S(T\mathbf{b}_j)]_{\mathcal{D}}$.
2. By (3.3) applied to $S$ with the vector $T\mathbf{b}_j \in W$: $[S(T\mathbf{b}_j)]_{\mathcal{D}} = [S]_{\mathcal{C},\mathcal{D}}\,[T\mathbf{b}_j]_{\mathcal{C}}$.
3. By (3.2), $[T\mathbf{b}_j]_{\mathcal{C}}$ is column $j$ of $[T]_{\mathcal{B},\mathcal{C}}$.
4. So column $j$ of the left side is $[S]_{\mathcal{C},\mathcal{D}}$ times column $j$ of $[T]_{\mathcal{B},\mathcal{C}}$. By LA-1 eq. (3.3), that is column $j$ of the product $[S]_{\mathcal{C},\mathcal{D}}[T]_{\mathcal{B},\mathcal{C}}$. Columns agree, so the matrices agree. ∎

■ **What just happened.** This is *why* matrices multiply the way they do. LA-1 defined $AB$ column by column so that $(AB)\mathbf{x} = A(B\mathbf{x})$ and remarked that row-times-column is the entry-level consequence. Now the reason is visible: "do $T$, then do $S$" must be a single linear map, and (3.4) says its matrix is the product. Any other rule for multiplying matrices would fail to track composition. Three things follow at once:

- *Order.* $S \circ T$ means $T$ first; in $[S][T]$ the matrix of $T$ sits nearest the vector. Matrices don't commute because composition doesn't (LA-1 Example 4.3: stretch-then-rotate is not rotate-then-stretch).
- *Identity.* The identity operator $I_V$ has matrix $[I_V]_{\mathcal{B},\mathcal{B}} = I$, since $[\mathbf{b}_j]_{\mathcal{B}} = \mathbf{e}_j$. With *different* input and output bases, $[I_V]_{\mathcal{B}',\mathcal{B}}$ is not the identity matrix. That matrix is the change-of-basis matrix, next.
- *Inverses.* If $T: V \to W$ has an inverse map $T^{-1}$ (which is linear, Theorem 3.17), then $[T^{-1}]_{\mathcal{C},\mathcal{B}}[T]_{\mathcal{B},\mathcal{C}} = [I_V]_{\mathcal{B}} = I$, so the matrix of the inverse is the inverse of the matrix.

### 3.6 Change of basis

Now one space $V$ and two bases: the **old** basis $\mathcal{B} = (\mathbf{b}_1, \dots, \mathbf{b}_n)$ and the **new** basis $\mathcal{B}' = (\mathbf{b}'_1, \dots, \mathbf{b}'_n)$. Every vector has two coordinate vectors, $[\mathbf{v}]_{\mathcal{B}}$ and $[\mathbf{v}]_{\mathcal{B}'}$. How are they related?

**Definition (Change-of-basis matrix).** $P := [I_V]_{\mathcal{B}',\mathcal{B}}$, the matrix of the identity map with *input* basis $\mathcal{B}'$ and *output* basis $\mathcal{B}$. By (3.2) its columns are $[I_V\mathbf{b}'_j]_{\mathcal{B}} = [\mathbf{b}'_j]_{\mathcal{B}}$:

$$
P = \Big( [\mathbf{b}'_1]_{\mathcal{B}} \;\; \cdots \;\; [\mathbf{b}'_n]_{\mathcal{B}} \Big), \qquad P\,[\mathbf{v}]_{\mathcal{B}'} = [\mathbf{v}]_{\mathcal{B}}. \tag{3.5}
$$

The second statement is (3.3) applied to $T = I_V$. **Read it carefully: the columns of $P$ are the new basis vectors written in old coordinates, and $P$ converts new coordinates to old coordinates.** $P$ is invertible (its columns are the coordinates of a basis, hence independent), and $P^{-1} = [I_V]_{\mathcal{B},\mathcal{B}'}$ converts old coordinates to new, because $[I]_{\mathcal{B},\mathcal{B}'}[I]_{\mathcal{B}',\mathcal{B}} = [I]_{\mathcal{B}'} = I$ by (3.4).

▶ **Try it first.** Let $\mathcal{B}$ be the standard basis of $\mathbb{R}^2$ and $\mathcal{B}' = \big((1, 1), (1, -1)\big)$. Write $P$. Find the $\mathcal{B}'$-coordinates of $\mathbf{v} = (3, 1)$ twice: by solving $a(1,1) + b(1,-1) = (3,1)$ directly, and by computing $P^{-1}\mathbf{v}$. They had better agree.

**The concrete case.** The new basis vectors in standard coordinates are themselves, so

$$
P = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}, \qquad P^2 = \begin{pmatrix} 2 & 0 \\ 0 & 2\end{pmatrix} = 2I, \qquad P^{-1} = \tfrac12 P = \begin{pmatrix} \tfrac12 & \tfrac12 \\ \tfrac12 & -\tfrac12 \end{pmatrix}.
$$

Directly: $a + b = 3$ and $a - b = 1$ give $a = 2, b = 1$. Via the matrix: $P^{-1}(3, 1)^\top = \tfrac12(3 + 1,\; 3 - 1)^\top = (2, 1)^\top$. ✓ Check (3.5): $P(2, 1)^\top = (2 + 1,\; 2 - 1)^\top = (3, 1)^\top$. ✓

The error to avoid: computing $P(3, 1)^\top = (4, 2)^\top$ and calling it "the new coordinates." It is not. $(4, 2)$ is the vector whose *new* coordinates are $(3, 1)$, i.e. $3\mathbf{b}'_1 + 1\mathbf{b}'_2$, expressed in old coordinates. Same matrix, opposite question. The test that never fails: **$P\mathbf{e}_j$ must be the $j$-th new basis vector in old coordinates.** If it isn't, you have $P^{-1}$.

**Theorem 3.8 (Change of basis for an operator).** Let $T: V \to V$ be linear, $A = [T]_{\mathcal{B}}$ its matrix in the old basis and $A' = [T]_{\mathcal{B}'}$ its matrix in the new basis. With $P$ as in (3.5),

$$
A' = P^{-1} A P. \tag{3.6}
$$

▶ **Try it first.** Insert the identity twice: $T = I_V \circ T \circ I_V$. Now apply (3.4) with the bases chosen as $\mathcal{B}' \to \mathcal{B} \to \mathcal{B} \to \mathcal{B}'$.

**Derivation.**

1. $T = I_V \circ T \circ I_V$ as maps on $V$.
2. Apply (3.4) twice, choosing the intermediate bases to be $\mathcal{B}$:
   $$[T]_{\mathcal{B}',\mathcal{B}'} = [I_V]_{\mathcal{B},\mathcal{B}'}\;[T]_{\mathcal{B},\mathcal{B}}\;[I_V]_{\mathcal{B}',\mathcal{B}}.$$
   The right factor takes $\mathcal{B}'$-coordinates in and gives $\mathcal{B}$-coordinates out; the middle factor works in $\mathcal{B}$; the left factor converts back.
3. $[I_V]_{\mathcal{B}',\mathcal{B}} = P$ and $[I_V]_{\mathcal{B},\mathcal{B}'} = P^{-1}$ by the definition and the remark after (3.5). So $A' = P^{-1}AP$. ∎

A second derivation, worth seeing because it is the one you will reconstruct under exam pressure: for any $\mathbf{v}$, $[T\mathbf{v}]_{\mathcal{B}'} = P^{-1}[T\mathbf{v}]_{\mathcal{B}} = P^{-1}A[\mathbf{v}]_{\mathcal{B}} = P^{-1}AP[\mathbf{v}]_{\mathcal{B}'}$. So the matrix $P^{-1}AP$ does to $\mathcal{B}'$-coordinates what $T$ does to vectors, and by Theorem 3.6 the matrix with that property is unique: it is $[T]_{\mathcal{B}'}$.

■ **What just happened.** Read $P^{-1}AP$ from the right, the direction the vector enters: $P$ converts new coordinates to old; $A$ applies $T$ in old coordinates; $P^{-1}$ converts the result back to new coordinates. For a map $T: V \to W$ with a basis change on each side ($P$ in $V$, $Q$ in $W$), the same argument gives $[T]_{\mathcal{B}',\mathcal{C}'} = Q^{-1}[T]_{\mathcal{B},\mathcal{C}}\,P$: convert input coordinates with $P$, apply, convert output coordinates back with $Q^{-1}$.

Example 4.2 runs (3.6) with numbers for the $90°$ rotation. Example 4.4 runs it for the stress tensor.

### 3.7 Similarity and its invariants

**Definition.** Two $n \times n$ matrices $A$ and $A'$ are **similar** if $A' = P^{-1}AP$ for some invertible $P$.

**Theorem 3.9 (Similar means same operator, different basis).** $A$ and $A'$ are similar if and only if they are the matrices of one operator $T$ in two bases.

*Proof.* ($\Leftarrow$) is Theorem 3.8. ($\Rightarrow$) Given $A' = P^{-1}AP$, let $T$ be the operator $\mathbf{x} \mapsto A\mathbf{x}$ on $\mathbb{R}^n$, so $[T]_{\text{std}} = A$, and let $\mathcal{B}'$ be the list of columns of $P$. They are independent (because $P$ is invertible, Theorem 3.18), so they form a basis, and the change-of-basis matrix from the standard basis to $\mathcal{B}'$ is exactly $P$ by (3.5). By Theorem 3.8, $[T]_{\mathcal{B}'} = P^{-1}AP = A'$. ∎

So any property shared by all matrices of an operator, regardless of basis, is a property of the *operator*. Which properties are those?

**Theorem 3.10 (Trace is a similarity invariant).** For any $X \in \mathbb{R}^{m\times n}$ and $Y \in \mathbb{R}^{n\times m}$, $\operatorname{tr}(XY) = \operatorname{tr}(YX)$. Consequently $\operatorname{tr}(P^{-1}AP) = \operatorname{tr}(A)$.

▶ **Try it first.** Write $(XY)_{ii}$ as a sum over $k$, then sum over $i$. Do the same for $(YX)_{kk}$. Compare the double sums.

**Derivation.**

1. $(XY)_{ii} = \sum_{k=1}^n x_{ik}y_{ki}$, so $\operatorname{tr}(XY) = \sum_{i=1}^m\sum_{k=1}^n x_{ik}y_{ki}$.
2. $(YX)_{kk} = \sum_{i=1}^m y_{ki}x_{ik}$, so $\operatorname{tr}(YX) = \sum_{k=1}^n\sum_{i=1}^m y_{ki}x_{ik}$.
3. These are the same $mn$ products $x_{ik}y_{ki}$, added in a different order. Finite sums can be reordered, so the traces agree.
4. Apply step 3 with $X = P^{-1}A$ and $Y = P$: $\operatorname{tr}(P^{-1}AP) = \operatorname{tr}\big((P^{-1}A)P\big) = \operatorname{tr}\big(P(P^{-1}A)\big) = \operatorname{tr}(A)$. ∎

■ **What just happened.** "The trace of the operator $T$" is now a meaningful number, $\operatorname{tr}T := \operatorname{tr}[T]_{\mathcal{B}}$ for any basis $\mathcal{B}$, even though it was defined by summing entries that depend on the basis. Note that $\operatorname{tr}(XYZ) = \operatorname{tr}(ZXY)$ (cyclic) but $\operatorname{tr}(XYZ) \neq \operatorname{tr}(YXZ)$ in general: only cyclic rotations are free.

Two more invariants, previewed:

- **Determinant.** $\det(P^{-1}AP) = \det(P^{-1})\det(A)\det(P) = \det A$, because LA-6 proves $\det(XY) = \det X\det Y$ and $\det(P^{-1}) = 1/\det P$. §3.13 shows the $2\times 2$ case as an area; LA-6 does the rest.
- **Eigenvalues.** A number $\lambda$ is an eigenvalue of $A$ if $A\mathbf{v} = \lambda\mathbf{v}$ for some $\mathbf{v} \neq \mathbf{0}$ (LA-7). If so, then $(P^{-1}AP)(P^{-1}\mathbf{v}) = P^{-1}A\mathbf{v} = \lambda\,P^{-1}\mathbf{v}$, and $P^{-1}\mathbf{v} \neq \mathbf{0}$ since $P^{-1}$ is invertible. So similar matrices have the same eigenvalues, with the eigenvectors converted by $P^{-1}$ (old coordinates to new, as they should be). LA-7 builds on exactly this line.
- **Rank** is also invariant: it is $\dim\operatorname{im}T$, a property of the operator (§3.8).

Not invariant: the entries, being diagonal, being symmetric (unless $P$ is orthogonal, LA-5), being upper triangular. And the converse fails: equal trace and determinant do *not* imply similar. $I$ and $\begin{pmatrix}1 & 1\\ 0 & 1\end{pmatrix}$ have trace $2$ and determinant $1$, but $P^{-1}IP = I$ for every $P$, so $I$ is similar only to itself (Exercise B4). The project checks trace, determinant, and eigenvalue invariance on random $P$.

### 3.8 Kernel and image

**Definition 3.11.** For $T: V \to W$ linear,
$$\ker T = \{\mathbf{v} \in V : T\mathbf{v} = \mathbf{0}\} \subseteq V, \qquad \operatorname{im} T = \{T\mathbf{v} : \mathbf{v} \in V\} \subseteq W.$$
The **rank** of $T$ is $\dim\operatorname{im}T$ and the **nullity** is $\dim\ker T$. For $T\mathbf{x} = A\mathbf{x}$, $\ker T = N(A)$ (the null space) and $\operatorname{im} T = C(A)$ (the column space), two of LA-3's four fundamental subspaces, and the rank of $T$ is the rank of $A$.

**Theorem 3.12.** $\ker T$ is a subspace of $V$ and $\operatorname{im} T$ is a subspace of $W$.

*Proof.* $T\mathbf{0} = \mathbf{0}$, so $\mathbf{0} \in \ker T$. If $\mathbf{u}, \mathbf{v} \in \ker T$ and $a, b \in \mathbb{R}$ then $T(a\mathbf{u} + b\mathbf{v}) = a\,T\mathbf{u} + b\,T\mathbf{v} = \mathbf{0}$, so $a\mathbf{u} + b\mathbf{v} \in \ker T$. For the image: $\mathbf{0} = T\mathbf{0} \in \operatorname{im}T$; if $\mathbf{w}_1 = T\mathbf{v}_1$ and $\mathbf{w}_2 = T\mathbf{v}_2$ then $a\mathbf{w}_1 + b\mathbf{w}_2 = T(a\mathbf{v}_1 + b\mathbf{v}_2) \in \operatorname{im}T$. Closure under combinations is LA-3's subspace test. ∎

**Theorem 3.13 (Injectivity test).** $T$ is injective (one-to-one) if and only if $\ker T = \{\mathbf{0}\}$.

*Proof.* ($\Rightarrow$) $T\mathbf{0} = \mathbf{0}$ always. If $T\mathbf{v} = \mathbf{0} = T\mathbf{0}$, injectivity gives $\mathbf{v} = \mathbf{0}$. ($\Leftarrow$) Suppose $T\mathbf{u} = T\mathbf{v}$. Then $T(\mathbf{u} - \mathbf{v}) = \mathbf{0}$, so $\mathbf{u} - \mathbf{v} \in \ker T = \{\mathbf{0}\}$, so $\mathbf{u} = \mathbf{v}$. ∎

■ **What just happened.** To check that a linear map never sends two inputs to the same output, you check *one* output, $\mathbf{0}$. For matrices this is LA-1 Exercise B3: two solutions of $A\mathbf{x} = \mathbf{b}$ differ by a solution of $A\mathbf{z} = \mathbf{0}$.

### 3.9 Rank–Nullity

**Theorem 3.14 (Rank–Nullity; Axler calls it the Fundamental Theorem of Linear Maps).** Let $V$ be finite-dimensional and $T: V \to W$ linear. Then $\operatorname{im}T$ is finite-dimensional and

$$
\dim V = \dim\ker T + \dim\operatorname{im}T. \tag{3.7}
$$

▶ **Try it first.** Take a basis of $\ker T$ and extend it to a basis of $V$ using LA-3's basis extension theorem. Apply $T$ to the vectors you added. Guess what they form, then prove both halves of the guess. The independence half is the one that needs thought.

**Derivation.**

1. $\ker T$ is a subspace of the finite-dimensional space $V$, so it has a basis $(\mathbf{u}_1, \dots, \mathbf{u}_k)$ with $k = \dim\ker T$ (LA-3). If $\ker T = \{\mathbf{0}\}$ then $k = 0$ and the list is empty; the argument below still works.
2. By the basis extension theorem there are vectors $\mathbf{v}_1, \dots, \mathbf{v}_r$ such that $(\mathbf{u}_1, \dots, \mathbf{u}_k, \mathbf{v}_1, \dots, \mathbf{v}_r)$ is a basis of $V$. Hence $\dim V = k + r$. It remains to show $\dim\operatorname{im}T = r$, and we do it by showing that $(T\mathbf{v}_1, \dots, T\mathbf{v}_r)$ is a basis of $\operatorname{im}T$.
3. *Spanning.* Let $\mathbf{w} \in \operatorname{im}T$, so $\mathbf{w} = T\mathbf{v}$ for some $\mathbf{v} \in V$. Expand $\mathbf{v}$ in the basis of step 2: $\mathbf{v} = \sum_i a_i\mathbf{u}_i + \sum_j c_j\mathbf{v}_j$. Then
   $$\mathbf{w} = T\mathbf{v} = \sum_i a_i\,T\mathbf{u}_i + \sum_j c_j\,T\mathbf{v}_j = \mathbf{0} + \sum_j c_j\,T\mathbf{v}_j,$$
   because each $\mathbf{u}_i$ is in the kernel. So every vector of the image is a combination of the $T\mathbf{v}_j$.
4. *Independence.* Suppose $\sum_j c_j\,T\mathbf{v}_j = \mathbf{0}$. By linearity, $T\left(\sum_j c_j\mathbf{v}_j\right) = \mathbf{0}$, so $\sum_j c_j\mathbf{v}_j \in \ker T$. The kernel is spanned by the $\mathbf{u}_i$, so $\sum_j c_j\mathbf{v}_j = \sum_i a_i\mathbf{u}_i$ for some scalars $a_i$. Move everything to one side:
   $$\sum_j c_j\mathbf{v}_j - \sum_i a_i\mathbf{u}_i = \mathbf{0}.$$
   This is a linear relation among the vectors of the basis from step 2, so every coefficient is zero. In particular every $c_j = 0$.
5. By steps 3 and 4, $(T\mathbf{v}_1, \dots, T\mathbf{v}_r)$ is a basis of $\operatorname{im}T$, so $\dim\operatorname{im}T = r$ and $\dim V = k + r = \dim\ker T + \dim\operatorname{im}T$. ∎

■ **What just happened.** The basis of $V$ was split into a part $T$ destroys and a part $T$ keeps, and the kept part is carried to a basis of the image *without collisions*. Step 4 is the one that students get wrong: you cannot conclude independence of the $T\mathbf{v}_j$ from independence of the $\mathbf{v}_j$ directly (a linear map can destroy independence; that is what the kernel *is*). You have to pull the relation back into $V$, notice it lands in the kernel, and use that the $\mathbf{v}_j$ were chosen *independent of the kernel*. §8 discusses why this proof feels circular and why it isn't.

Three remarks.

- For an $m \times n$ matrix $A$, (3.7) reads $n = \dim N(A) + \operatorname{rank}A$. LA-3 got this by counting: each free variable gives one basis vector of $N(A)$, each pivot one basis vector of $C(A)$, and every column is one or the other. That argument needs elimination. This one needs no coordinates and works for $D$ on polynomials or $A \mapsto A - A^\top$ on matrices equally well.
- $\operatorname{rank}T \le \min(\dim V, \dim W)$: the image is a subspace of $W$, and by (3.7) its dimension is at most $\dim V$.
- The right side of (3.7) involves only the *input* dimension. $\dim W$ puts a ceiling on the rank, nothing more.

**Corollary 3.15 (Equal dimensions).** Suppose $\dim V = \dim W$ and $T: V \to W$ is linear. Then $T$ is injective $\iff$ $T$ is surjective $\iff$ $T$ is bijective.

*Proof.* $T$ injective $\iff \ker T = \{\mathbf{0}\}$ (Theorem 3.13) $\iff \dim\ker T = 0 \iff \dim\operatorname{im}T = \dim V$ (by (3.7)) $\iff \dim\operatorname{im}T = \dim W \iff \operatorname{im}T = W$ (a subspace with the full dimension is the whole space, LA-3) $\iff T$ surjective. Bijective is both. ∎

This fails the moment the dimensions differ. $D: \mathcal{P}_3 \to \mathcal{P}_3$ is neither injective (constants die) nor surjective ($x^3$ is not a derivative of a cubic). $J: \mathcal{P}_2 \to \mathcal{P}_3$ is injective but not surjective (nothing integrates to $1$ from $0$). $D: \mathcal{P}_3 \to \mathcal{P}_2$ is surjective but not injective. Exercise B2 computes all three.

### 3.10 Isomorphism

**Definition 3.16.** An **isomorphism** is a bijective linear map $T: V \to W$. If one exists, $V$ and $W$ are **isomorphic**, written $V \cong W$. A bijective operator $T: V \to V$ is called **invertible**.

**Theorem 3.17.** (a) The inverse of an isomorphism is linear. (b) Every $n$-dimensional real vector space is isomorphic to $\mathbb{R}^n$. (c) $V \cong W$ if and only if $\dim V = \dim W$.

*Proof.* (a) Let $S = T^{-1}: W \to V$, which exists as a function because $T$ is bijective. For $\mathbf{w}_1, \mathbf{w}_2 \in W$ and scalars $a, b$:
$$T\big(a\,S\mathbf{w}_1 + b\,S\mathbf{w}_2\big) = a\,TS\mathbf{w}_1 + b\,TS\mathbf{w}_2 = a\mathbf{w}_1 + b\mathbf{w}_2.$$
Apply $S$ to both sides: $a\,S\mathbf{w}_1 + b\,S\mathbf{w}_2 = S(a\mathbf{w}_1 + b\mathbf{w}_2)$. That is linearity of $S$.

(b) Pick a basis $\mathcal{B}$; $\kappa_{\mathcal{B}}: V \to \mathbb{R}^n$ is a linear bijection by Theorem 3.4.

(c) ($\Rightarrow$) If $T: V \to W$ is an isomorphism, $\ker T = \{\mathbf{0}\}$ and $\operatorname{im}T = W$, so (3.7) gives $\dim V = 0 + \dim W$. ($\Leftarrow$) If $\dim V = \dim W = n$, then $\kappa_{\mathcal{C}}^{-1}\circ\kappa_{\mathcal{B}}: V \to \mathbb{R}^n \to W$ is a composition of two linear bijections (using (a) for the inverse), hence a linear bijection. Concretely it sends $\mathbf{b}_j \mapsto \mathbf{c}_j$ and has matrix $I$. ∎

■ **What just happened.** Up to relabeling, there is exactly one real vector space of each finite dimension. $\mathcal{P}_3$, $\mathbb{R}^{2\times 2}$, and $\mathbb{R}^4$ are the same space wearing different clothes; any linear-algebra question about one can be moved to another and back. This is why the four fundamental subspaces of LA-3, defined for matrices, are enough to study every finite-dimensional linear map: pick bases, and the map *is* a matrix.

### 3.11 Invertibility of a square matrix: all the equivalences in one place

**Theorem 3.18 (Invertible matrix theorem).** For an $n \times n$ real matrix $A$, the following are equivalent.

- (a) $A$ is invertible: there is $B$ with $AB = BA = I$.
- (b) The map $\mathbf{x} \mapsto A\mathbf{x}$ is an isomorphism of $\mathbb{R}^n$.
- (c) $N(A) = \{\mathbf{0}\}$: the only solution of $A\mathbf{x} = \mathbf{0}$ is $\mathbf{x} = \mathbf{0}$.
- (d) The columns of $A$ are linearly independent.
- (e) The columns of $A$ span $\mathbb{R}^n$, i.e. $A\mathbf{x} = \mathbf{b}$ is solvable for every $\mathbf{b}$.
- (f) $\operatorname{rank}A = n$.
- (g) $A$ has $n$ pivots (its reduced row echelon form is $I$).
- (h) The rows of $A$ are linearly independent.
- (i) $A\mathbf{x} = \mathbf{b}$ has exactly one solution for every $\mathbf{b}$.

Previews, proved later: (j) $\det A \neq 0$ (LA-6). (k) $0$ is not an eigenvalue of $A$ (LA-7; (k)$\iff$(c) is immediate once "eigenvalue $0$" is unpacked as "$A\mathbf{v} = \mathbf{0}$ with $\mathbf{v} \neq \mathbf{0}$").

*Proof.* Let $T\mathbf{x} = A\mathbf{x}$.

(c)$\iff$(d): $A\mathbf{x} = \sum_j x_j\mathbf{a}_j$ (LA-1 eq. 3.1), so "$A\mathbf{x} = \mathbf{0}$ only for $\mathbf{x} = \mathbf{0}$" says exactly "the only combination of the columns that vanishes is the trivial one," which is LA-3's definition of independence.

(c)$\iff$(f): Rank–Nullity, $n = \dim N(A) + \operatorname{rank}A$.

(f)$\iff$(e): $C(A)$ is a subspace of $\mathbb{R}^n$ of dimension $\operatorname{rank}A$; it equals $\mathbb{R}^n$ iff its dimension is $n$ (LA-3). And $C(A) = \mathbb{R}^n$ means every $\mathbf{b}$ is a combination of the columns, i.e. $A\mathbf{x} = \mathbf{b}$ is solvable (LA-1 Theorem 3.6).

(c)$\iff$(e) directly: Corollary 3.15, since input and output both have dimension $n$.

(f)$\iff$(h): LA-3 proved row rank equals column rank, $\dim C(A^\top) = \dim C(A)$. The $n$ rows are independent iff they span an $n$-dimensional space iff row rank $= n$.

(c)$\iff$(g): LA-2. The solutions of $A\mathbf{x} = \mathbf{0}$ are parametrized by the free variables, one per non-pivot column. There are none iff every column has a pivot iff there are $n$ pivots.

(b)$\iff$(c) and (e): bijective means injective and surjective; injective is (c) by Theorem 3.13; surjective means $\operatorname{im}T = C(A) = \mathbb{R}^n$, which is (e).

(a)$\iff$(b): If $AB = BA = I$ then $\mathbf{x} \mapsto B\mathbf{x}$ undoes $\mathbf{x} \mapsto A\mathbf{x}$ in both orders, so $T$ is bijective. Conversely if $T$ is bijective, $T^{-1}$ is linear (Theorem 3.17a), hence by Theorem 3.6 it is $\mathbf{x} \mapsto B\mathbf{x}$ for $B = [T^{-1}]_{\text{std}}$, and (3.4) gives $AB = [T\circ T^{-1}] = I$ and $BA = [T^{-1}\circ T] = I$.

(i)$\iff$(c) and (e): existence for every $\mathbf{b}$ is (e). Given existence, uniqueness for every $\mathbf{b}$ is equivalent to (c) because two solutions of $A\mathbf{x} = \mathbf{b}$ differ by an element of $N(A)$ (LA-1 Exercise B3), and conversely $N(A) \neq \{\mathbf{0}\}$ gives $A\mathbf{x} = \mathbf{0}$ two solutions. ∎

A corollary with real content: **for square matrices, a one-sided inverse is two-sided.** If $AB = I$, then $B\mathbf{x} = \mathbf{0}$ implies $\mathbf{x} = AB\mathbf{x} = \mathbf{0}$, so $B$ satisfies (c), hence (a), and then $A = A(BB^{-1}) = (AB)B^{-1} = B^{-1}$, so $BA = I$ too. Exercise B5 asks you to write this out and to say exactly where Rank–Nullity was used.

### 3.12 The space of linear maps

**Theorem 3.19.** $\mathcal{L}(V, W)$ is a vector space under $(S + T)\mathbf{v} := S\mathbf{v} + T\mathbf{v}$ and $(cT)\mathbf{v} := c\,T\mathbf{v}$, and if $\dim V = n$, $\dim W = m$, then $\dim\mathcal{L}(V, W) = mn$.

▶ **Try it first.** Fix bases and consider $\Phi: \mathcal{L}(V, W) \to \mathbb{R}^{m\times n}$, $T \mapsto [T]_{\mathcal{B},\mathcal{C}}$. Show it is linear, injective, and surjective. What does Theorem 3.17(c) then say?

**Derivation.**

1. *Vector space.* $S + T$ is linear: $(S+T)(a\mathbf{u} + b\mathbf{v}) = S(a\mathbf{u} + b\mathbf{v}) + T(a\mathbf{u} + b\mathbf{v}) = a(S\mathbf{u} + T\mathbf{u}) + b(S\mathbf{v} + T\mathbf{v}) = a(S+T)\mathbf{u} + b(S+T)\mathbf{v}$. So is $cT$: $(cT)(a\mathbf{u} + b\mathbf{v}) = c(aT\mathbf{u} + bT\mathbf{v}) = a(cT)\mathbf{u} + b(cT)\mathbf{v}$. The eight axioms of LA-3 hold because they hold in $W$ at every $\mathbf{v}$ ("pointwise"); the zero map $\mathbf{v} \mapsto \mathbf{0}$ is the zero vector.
2. *$\Phi$ is linear.* Column $j$ of $[S + T]_{\mathcal{B},\mathcal{C}}$ is $[(S+T)\mathbf{b}_j]_{\mathcal{C}} = [S\mathbf{b}_j]_{\mathcal{C}} + [T\mathbf{b}_j]_{\mathcal{C}}$ by Theorem 3.4, which is column $j$ of $[S] + [T]$. And column $j$ of $[cT]$ is $[cT\mathbf{b}_j]_{\mathcal{C}} = c[T\mathbf{b}_j]_{\mathcal{C}}$, so $[cT] = c[T]$.
3. *$\Phi$ is injective.* If $[T]_{\mathcal{B},\mathcal{C}} = 0$ then every $[T\mathbf{b}_j]_{\mathcal{C}} = \mathbf{0}$, so every $T\mathbf{b}_j = \mathbf{0}$, so $T$ agrees with the zero map on a basis and equals it by the uniqueness half of Theorem 3.2. So $\ker\Phi = \{0\}$ and Theorem 3.13 applies.
4. *$\Phi$ is surjective.* Given any $M \in \mathbb{R}^{m\times n}$, set $\mathbf{w}_j = \sum_i M_{ij}\mathbf{c}_i$ and let $T$ be the linear extension with $T\mathbf{b}_j = \mathbf{w}_j$ (existence half of Theorem 3.2). Then $[T\mathbf{b}_j]_{\mathcal{C}}$ is column $j$ of $M$, so $\Phi(T) = M$.
5. $\Phi$ is an isomorphism, so by Theorem 3.17(c), $\dim\mathcal{L}(V, W) = \dim\mathbb{R}^{m\times n} = mn$. ∎

■ **What just happened.** The "matrix of a map" construction is itself a linear isomorphism between maps and matrices. Everything in this module is one sentence: *$\mathcal{L}(V, W) \cong \mathbb{R}^{m\times n}$, with composition going to multiplication, and the isomorphism depends on the choice of bases by the rule $P^{-1}AP$.* The case $W = \mathbb{R}$ gives $\mathcal{L}(V, \mathbb{R})$, the **dual space** $V^*$, of dimension $n$; LA-12 is about it.

### 3.13 Some operators worth knowing by sight

**Nilpotent operators.** $T$ is **nilpotent** if $T^k = 0$ for some $k$ (the smallest such $k$ is the *index*). Differentiation on $\mathcal{P}_3$ is the standard example: each application lowers the degree by one or kills the polynomial, so $D^4 = 0$ on $\mathcal{P}_3$, while $D^3x^3 = 6 \neq 0$, so the index is $4$. Its matrix is in Example 4.1. In the basis $(1, x, x^2/2, x^3/6)$ it is the **shift matrix** with $1$s immediately above the diagonal, which is the model for every nilpotent operator (LA-11, Jordan form). Exercise C2 proves the index of a nilpotent operator is at most $\dim V$.

**Projections (idempotents).** $P$ is a **projection** if $P^2 = P$. Geometric meaning: once a vector has landed in the image, projecting again does nothing. Precisely: if $\mathbf{w} \in \operatorname{im}P$, say $\mathbf{w} = P\mathbf{v}$, then $P\mathbf{w} = P^2\mathbf{v} = P\mathbf{v} = \mathbf{w}$. So $P$ acts as the identity on its image and as zero on its kernel, and every vector splits as $\mathbf{v} = (\mathbf{v} - P\mathbf{v}) + P\mathbf{v}$ with the first piece in the kernel ($P(\mathbf{v} - P\mathbf{v}) = P\mathbf{v} - P^2\mathbf{v} = \mathbf{0}$) and the second in the image. In a basis made of a basis of the image followed by a basis of the kernel, the matrix is $\operatorname{diag}(1, \dots, 1, 0, \dots, 0)$ with $\operatorname{rank}P$ ones. Example 4.3 shows this for the line $y = 2x$; Exercise C3 proves it in general. LA-5 adds the orthogonality condition $P = P^\top$, which pins down *which* projection onto a given image.

**Rotation, reflection, shear in $\mathbb{R}^2$.** By Definition 3.5 with standard bases, write down where $\mathbf{e}_1$ and $\mathbf{e}_2$ go and those are the columns.

- Rotation by $\theta$: $\mathbf{e}_1 \mapsto (\cos\theta, \sin\theta)$, $\mathbf{e}_2 \mapsto (-\sin\theta, \cos\theta)$, so
  $$R_\theta = \begin{pmatrix}\cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{pmatrix}.$$
  $R_\theta R_\varphi = R_{\theta+\varphi}$ by (3.4), and multiplying out the left side proves the angle-addition formulas (Exercise A5).
- Reflection across the line through the origin at angle $\varphi$: the mirror image of $\mathbf{e}_1$ (at angle $0$) is at angle $2\varphi$; the mirror image of $\mathbf{e}_2$ (at angle $90°$) is at angle $2\varphi - 90°$. So
  $$F_\varphi = \begin{pmatrix}\cos 2\varphi & \sin 2\varphi \\ \sin 2\varphi & -\cos 2\varphi\end{pmatrix}, \qquad F_\varphi^2 = I.$$
- Horizontal shear with factor $k$: $\mathbf{e}_1 \mapsto \mathbf{e}_1$, $\mathbf{e}_2 \mapsto (k, 1)$, so $S_k = \begin{pmatrix}1 & k \\ 0 & 1\end{pmatrix}$. Horizontal lines slide sideways by $k$ times their height.
- Projection onto the line through the unit vector $\mathbf{u}$: $P = \mathbf{u}\mathbf{u}^\top$, since $P\mathbf{v} = \mathbf{u}(\mathbf{u}^\top\mathbf{v}) = (\mathbf{u}\cdot\mathbf{v})\mathbf{u}$.

**Area: a determinant preview.** The unit square has corners $\mathbf{0}, \mathbf{e}_1, \mathbf{e}_1 + \mathbf{e}_2, \mathbf{e}_2$. A linear map with matrix $A = \begin{pmatrix} a & b \\ c & d\end{pmatrix}$ sends it to the parallelogram with sides $\mathbf{u} = A\mathbf{e}_1 = (a, c)$ and $\mathbf{v} = A\mathbf{e}_2 = (b, d)$ (the orange shape in `figures/grid_deformations.png`).

**Theorem 3.20 (Area of the image of the unit square).** The parallelogram spanned by $(a, c)$ and $(b, d)$ has area $|ad - bc|$.

▶ **Try it first.** Area is base times height, $\|\mathbf{u}\|\,\|\mathbf{v}\|\,|\sin\theta|$. Square it, replace $\sin^2$ by $1 - \cos^2$, and use LA-1 eq. (3.4) for $\cos\theta$.

**Derivation.**

1. Area $= \|\mathbf{u}\|\,\|\mathbf{v}\||\sin\theta|$ where $\theta$ is the angle between the sides (base $\|\mathbf{u}\|$, height $\|\mathbf{v}\||\sin\theta|$).
2. Area$^2 = \|\mathbf{u}\|^2\|\mathbf{v}\|^2(1 - \cos^2\theta) = \|\mathbf{u}\|^2\|\mathbf{v}\|^2 - (\mathbf{u}\cdot\mathbf{v})^2$, using $\mathbf{u}\cdot\mathbf{v} = \|\mathbf{u}\|\|\mathbf{v}\|\cos\theta$.
3. With $\mathbf{u} = (a, c)$, $\mathbf{v} = (b, d)$: $(a^2 + c^2)(b^2 + d^2) - (ab + cd)^2 = a^2b^2 + a^2d^2 + b^2c^2 + c^2d^2 - a^2b^2 - 2abcd - c^2d^2 = a^2d^2 - 2abcd + b^2c^2 = (ad - bc)^2$. Take the square root. ∎

$$
\operatorname{area}\big(A(\text{unit square})\big) = |ad - bc|. \tag{3.8}
$$

So: rotation, $\cos^2 + \sin^2 = 1$, area preserved. Reflection, $-\cos^2 2\varphi - \sin^2 2\varphi = -1$, area preserved but the sign is negative: the square has been flipped over. Shear, $1\cdot 1 - k\cdot 0 = 1$, area preserved even though the shape is badly distorted. Projection $\tfrac15\begin{pmatrix}1 & 2 \\ 2 & 4\end{pmatrix}$: $\tfrac{1}{25}(4 - 4) = 0$, the square is flattened to a segment. $\begin{pmatrix} 2 & 1 \\ 1 & 3\end{pmatrix}$: area $5$. The project computes each area from the polygon and compares.

The number $ad - bc$ is the **determinant**, LA-6's subject. LA-6 proves it scales the area of *every* region, not only the unit square; that its sign records whether orientation is preserved; that $\det(AB) = \det A\det B$ (compose the maps, multiply the scale factors); and hence that $\det$ is a similarity invariant. The sign is the explanation of Example 4.2: the basis $(1, 1), (1, -1)$ has $\det P = -2 < 0$, it is a left-handed frame, so a counterclockwise rotation looks clockwise in its coordinates.

## 4. Worked examples

### Example 4.1: The differentiation operator on $\mathcal{P}_3$ in two bases

$D: \mathcal{P}_3 \to \mathcal{P}_3$, $Dp = p'$. **Monomial basis** $\mathcal{B} = (1, x, x^2, x^3)$. Apply $D$ to each basis vector and write the result in the same basis (Definition 3.5):

$$
D(1) = 0 = (0, 0, 0, 0),\quad D(x) = 1 = (1, 0, 0, 0),\quad D(x^2) = 2x = (0, 2, 0, 0),\quad D(x^3) = 3x^2 = (0, 0, 3, 0).
$$

These are the columns:

$$
[D]_{\mathcal{B}} = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \\ 0 & 0 & 0 & 0 \end{pmatrix}.
$$

Check (3.3) on $p = 2 - x + 5x^3$, coordinates $(2, -1, 0, 5)^\top$: $[D]_{\mathcal{B}}(2, -1, 0, 5)^\top = (-1, 0, 15, 0)^\top$, and indeed $p' = -1 + 15x^2$. ✓

Powers: $[D]^2$ has entries $2$ and $6$ in positions $(1,3)$ and $(2,4)$ (second derivatives: $D^2x^2 = 2$, $D^2x^3 = 6x$); $[D]^3$ has the single entry $6$ in position $(1, 4)$; $[D]^4 = 0$. So $D$ is nilpotent of index $4$ on a $4$-dimensional space.

Kernel and image: $\ker D$ = constants, dimension $1$. $\operatorname{im}D = \mathcal{P}_2$ (every polynomial of degree $\le 2$ is the derivative of a cubic), dimension $3$. Rank–Nullity: $4 = 1 + 3$. ✓ $D$ is neither injective nor surjective as a map $\mathcal{P}_3 \to \mathcal{P}_3$, consistent with Corollary 3.15 only because it is neither.

**Factorial basis** $\mathcal{B}' = (1, x, x^2/2, x^3/6)$. Now
$$D(1) = 0,\qquad D(x) = 1 = 1\cdot\mathbf{b}'_1,\qquad D(x^2/2) = x = 1\cdot\mathbf{b}'_2,\qquad D(x^3/6) = x^2/2 = 1\cdot\mathbf{b}'_3,$$
so
$$
[D]_{\mathcal{B}'} = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 \end{pmatrix},
$$
the shift matrix: it moves each coordinate one slot down and drops the first. Same operator, cleaner matrix, because the basis was chosen so that $D$ carries each basis vector to the previous one. (This is the Taylor-series basis: $f = \sum_k f^{(k)}(0)\,x^k/k!$, so the $\mathcal{B}'$-coordinates of $f$ are its derivatives at $0$, and differentiation shifts them.)

Consistency with Theorem 3.8: the change-of-basis matrix has columns $[\mathbf{b}'_j]_{\mathcal{B}}$, i.e. $P = \operatorname{diag}(1, 1, \tfrac12, \tfrac16)$ and $P^{-1} = \operatorname{diag}(1, 1, 2, 6)$. Then $[D]_{\mathcal{B}}P$ has columns $\mathbf{0}, (1,0,0,0)^\top, \tfrac12(0,2,0,0)^\top = (0,1,0,0)^\top, \tfrac16(0,0,3,0)^\top = (0,0,\tfrac12,0)^\top$, and multiplying by $P^{-1}$ on the left scales row $3$ by $2$ (and row $4$ by $6$, which is all zeros), turning the $\tfrac12$ into a $1$. Result: the shift matrix. ✓ The project computes both matrices from the definition, checks them against SymPy, and checks $P^{-1}AP$.

### Example 4.2: Rotation by $90°$ in the basis $(1, 1), (1, -1)$

$T$ = rotation by $+90°$. Standard matrix $A = \begin{pmatrix} 0 & -1 \\ 1 & 0\end{pmatrix}$ (columns $T\mathbf{e}_1 = \mathbf{e}_2$, $T\mathbf{e}_2 = -\mathbf{e}_1$). New basis $\mathbf{b}'_1 = (1, 1)$, $\mathbf{b}'_2 = (1, -1)$. From §3.6, $P = \begin{pmatrix} 1 & 1 \\ 1 & -1\end{pmatrix}$ and $P^{-1} = \tfrac12\begin{pmatrix} 1 & 1 \\ 1 & -1\end{pmatrix}$.

*Step 1: $AP$.* Columns are $A\mathbf{b}'_1 = (0\cdot 1 - 1\cdot 1,\; 1\cdot 1 + 0\cdot 1) = (-1, 1)$ and $A\mathbf{b}'_2 = (0\cdot 1 - 1\cdot(-1),\; 1\cdot 1 + 0) = (1, 1)$:
$$AP = \begin{pmatrix} -1 & 1 \\ 1 & 1\end{pmatrix}.$$
(These are the rotated basis vectors in old coordinates: $(1,1)$ rotated by $90°$ is $(-1, 1)$. ✓)

*Step 2: $P^{-1}(AP)$.*
$$P^{-1}AP = \tfrac12\begin{pmatrix} 1 & 1 \\ 1 & -1\end{pmatrix}\begin{pmatrix} -1 & 1 \\ 1 & 1\end{pmatrix} = \tfrac12\begin{pmatrix} -1 + 1 & 1 + 1 \\ -1 - 1 & 1 - 1\end{pmatrix} = \tfrac12\begin{pmatrix} 0 & 2 \\ -2 & 0\end{pmatrix} = \begin{pmatrix} 0 & 1 \\ -1 & 0\end{pmatrix}.$$

*Check against the definition.* Column $1$ should be $[T\mathbf{b}'_1]_{\mathcal{B}'}$. $T\mathbf{b}'_1 = (-1, 1)$; solve $a(1,1) + b(1,-1) = (-1, 1)$: $a + b = -1$, $a - b = 1$, so $a = 0$, $b = -1$. Column $(0, -1)^\top$. ✓ Column $2$: $T\mathbf{b}'_2 = (1, 1) = 1\cdot\mathbf{b}'_1 + 0\cdot\mathbf{b}'_2$. Column $(1, 0)^\top$. ✓ (This is exactly the computation of §2.)

*What it means.* In the standard basis, "rotation by $90°$" reads $\mathbf{e}_1 \mapsto \mathbf{e}_2$, $\mathbf{e}_2 \mapsto -\mathbf{e}_1$. In the new basis it reads $\mathbf{b}'_1 \mapsto -\mathbf{b}'_2$, $\mathbf{b}'_2 \mapsto \mathbf{b}'_1$, which is the pattern of rotation by $-90°$. Nothing is wrong. The frame $(\mathbf{b}'_1, \mathbf{b}'_2)$ is left-handed ($\det P = -2 < 0$: going from $\mathbf{b}'_1$ to $\mathbf{b}'_2$ is a clockwise turn), so a counterclockwise motion of the plane has clockwise-looking coordinates. Invariants: $\operatorname{tr} = 0$ and $\det = 1$ for both matrices, as Theorem 3.10 and §3.13 require.

*The wrong direction, for the record.* For this particular $P$, $P^{-1} = \tfrac12 P$, so $PAP^{-1} = \tfrac12 PAP = P^{-1}AP$: the wrong-direction product gives the right answer by accident, and a direction error would go unnoticed. In general $PAP^{-1}$ is the matrix of $T$ in the basis whose vectors (in old coordinates) are the columns of $P^{-1}$, which is not what was asked. Exercise B3 has a basis where the two products differ visibly, and where even the trace and determinant cannot tell them apart.

### Example 4.3: Projection onto the line $y = 2x$

Unit direction $\mathbf{u} = (1, 2)/\sqrt5$. By §3.13, $P = \mathbf{u}\mathbf{u}^\top = \tfrac15\begin{pmatrix} 1 \\ 2\end{pmatrix}\begin{pmatrix} 1 & 2\end{pmatrix} = \tfrac15\begin{pmatrix} 1 & 2 \\ 2 & 4\end{pmatrix}$.

*$P^2 = P$:*
$$P^2 = \tfrac1{25}\begin{pmatrix} 1 & 2 \\ 2 & 4\end{pmatrix}\begin{pmatrix} 1 & 2 \\ 2 & 4\end{pmatrix} = \tfrac1{25}\begin{pmatrix} 1 + 4 & 2 + 8 \\ 2 + 8 & 4 + 16\end{pmatrix} = \tfrac1{25}\begin{pmatrix} 5 & 10 \\ 10 & 20\end{pmatrix} = \tfrac15\begin{pmatrix} 1 & 2 \\ 2 & 4\end{pmatrix} = P. \checkmark$$
Or without arithmetic: $P^2 = \mathbf{u}(\mathbf{u}^\top\mathbf{u})\mathbf{u}^\top = \mathbf{u}\cdot 1\cdot\mathbf{u}^\top = P$.

*Kernel.* $P\mathbf{v} = \mathbf{0}$ iff $\mathbf{u}\cdot\mathbf{v} = 0$ iff $v_1 + 2v_2 = 0$: the line through $(2, -1)$. Check: $P(2, -1)^\top = \tfrac15(2 - 2,\; 4 - 4)^\top = \mathbf{0}$. ✓ Nullity $1$.

*Image.* Every $P\mathbf{v}$ is a multiple of $\mathbf{u}$, and $P\mathbf{u} = \mathbf{u}$, so the image is the line $y = 2x$. Check: $P(1, 2)^\top = \tfrac15(1 + 4,\; 2 + 8)^\top = (1, 2)^\top$. ✓ Rank $1$. Rank–Nullity: $2 = 1 + 1$. ✓

*Adapted basis.* Take $\mathcal{B}' = \big((1, 2), (2, -1)\big)$, image direction first, kernel direction second. $T\mathbf{b}'_1 = \mathbf{b}'_1$, $T\mathbf{b}'_2 = \mathbf{0}$, so $[P]_{\mathcal{B}'} = \operatorname{diag}(1, 0)$ by Definition 3.5. Via (3.6) with $Q = \begin{pmatrix} 1 & 2 \\ 2 & -1\end{pmatrix}$: $PQ$ has columns $P(1,2)^\top = (1, 2)^\top$ and $P(2,-1)^\top = \mathbf{0}$, so $PQ = \begin{pmatrix} 1 & 0 \\ 2 & 0\end{pmatrix}$, and $Q^{-1} = \tfrac{1}{-5}\begin{pmatrix} -1 & -2 \\ -2 & 1\end{pmatrix} = \tfrac15\begin{pmatrix} 1 & 2 \\ 2 & -1\end{pmatrix}$, giving $Q^{-1}PQ = \tfrac15\begin{pmatrix} 1 + 4 & 0 \\ 2 - 2 & 0\end{pmatrix} = \operatorname{diag}(1, 0)$. ✓ Trace $1 = \tfrac15 + \tfrac45$ and determinant $0$ in both bases.

*Geometry.* $P$ is the identity on its image and zero on its kernel; the plane is the direct sum of the two lines. That is what $P^2 = P$ means. The grid picture (`figures/grid_deformations.png`, bottom left) shows the whole plane collapsing onto one line.

### Example 4.4: Plane stress in a rotated frame (Mohr's circle is a change of basis)

At a point in a loaded plate the stress state is the symmetric matrix
$$\sigma = \begin{pmatrix} \sigma_x & \tau_{xy} \\ \tau_{xy} & \sigma_y\end{pmatrix} = \begin{pmatrix} 80 & 40 \\ 40 & 20\end{pmatrix}\ \text{MPa}$$
in the $x, y$ frame. Cauchy's theorem says $\sigma$ is a linear map: for a cut with unit normal $\mathbf{n}$, the traction (force per area) across the cut is $\mathbf{t} = \sigma\mathbf{n}$. The entries are the components of that map in the $x, y$ frame. A different engineer uses axes rotated by $\theta$: $\mathbf{e}'_1 = (\cos\theta, \sin\theta)$, $\mathbf{e}'_2 = (-\sin\theta, \cos\theta)$. Same plate, same map, different numbers. What are they?

*The change-of-basis step.* The new basis vectors in old coordinates are the columns of $R_\theta$, so $P = R_\theta$ in (3.5). Since $R_\theta^\top R_\theta = I$, $P^{-1} = R_\theta^\top$. By (3.6),
$$\sigma' = R_\theta^\top\,\sigma\,R_\theta.$$
Write $c = \cos\theta$, $s = \sin\theta$. First $\sigma R_\theta$ has columns $\sigma(c, s)^\top = (\sigma_x c + \tau s,\; \tau c + \sigma_y s)^\top$ and $\sigma(-s, c)^\top = (-\sigma_x s + \tau c,\; -\tau s + \sigma_y c)^\top$. Then the $(1,1)$ entry of $R_\theta^\top(\sigma R_\theta)$ is row $(c, s)$ dotted with column one:
$$\sigma'_{x} = c(\sigma_x c + \tau s) + s(\tau c + \sigma_y s) = \sigma_x c^2 + \sigma_y s^2 + 2\tau cs,$$
and the $(1,2)$ entry is row $(c, s)$ dotted with column two:
$$\tau' = c(-\sigma_x s + \tau c) + s(-\tau s + \sigma_y c) = -(\sigma_x - \sigma_y)cs + \tau(c^2 - s^2).$$
Using $c^2 = \tfrac12(1 + \cos 2\theta)$, $s^2 = \tfrac12(1 - \cos 2\theta)$, $2cs = \sin 2\theta$:
$$\sigma'_x = \frac{\sigma_x + \sigma_y}{2} + \frac{\sigma_x - \sigma_y}{2}\cos 2\theta + \tau\sin 2\theta, \qquad \tau' = -\frac{\sigma_x - \sigma_y}{2}\sin 2\theta + \tau\cos 2\theta.$$
These are the stress-transformation equations of every mechanics-of-materials course. They are (3.6), nothing more.

*Mohr's circle.* As $\theta$ varies, $(\sigma'_x, \tau')$ is a point at distance $\sqrt{\big(\tfrac{\sigma_x - \sigma_y}{2}\big)^2 + \tau^2}$ from $\big(\tfrac{\sigma_x + \sigma_y}{2}, 0\big)$, because the two expressions are a $(\cos, \sin)$ pair of the angle $2\theta$ against a fixed vector. It traces a circle. With the numbers: centre $\tfrac{80 + 20}{2} = 50$, radius $\sqrt{30^2 + 40^2} = 50$.

*Principal frame.* $\tau' = 0$ when $\tan 2\theta_p = \dfrac{2\tau}{\sigma_x - \sigma_y} = \dfrac{80}{60} = \dfrac43$, so $\cos 2\theta_p = \tfrac35$, $\sin 2\theta_p = \tfrac45$, $\theta_p = 26.57°$. Then
$$\sigma'_x = 50 + 30\cdot\tfrac35 + 40\cdot\tfrac45 = 50 + 18 + 32 = 100, \qquad \tau' = -30\cdot\tfrac45 + 40\cdot\tfrac35 = -24 + 24 = 0,$$
and $\sigma'_y = \operatorname{tr}\sigma - \sigma'_x = 100 - 100 = 0$ (trace is invariant, Theorem 3.10). So in the principal frame $\sigma' = \operatorname{diag}(100, 0)$ MPa: pure tension of $100$ MPa along one axis, nothing along the other. The determinant checks: $\det\sigma = 80\cdot 20 - 40^2 = 0 = 100\cdot 0$. At $\theta = 45°$: $\sigma'_x = 50 + 0 + 40 = 90$, $\tau' = -30$, $\sigma'_y = 10$; the matrix is $\begin{pmatrix} 90 & -30 \\ -30 & 10\end{pmatrix}$, trace $100$ again.

In LA-8 the principal frame will be the eigenvector basis of the symmetric matrix $\sigma$, and the spectral theorem will say such a frame always exists and is orthogonal. For now the point is narrower and more important: the engineer's transformation formulas are $P^{-1}AP$ with $P$ a rotation. The project draws the circle (`figures/mohr_circle.png`) and checks every number above.

## 5. Why it matters and where it connects

- **LA-5.** Projections are the operators with $P^2 = P$ (§3.13); orthogonal projections add $P = P^\top$. Gram–Schmidt produces orthonormal bases, for which the change-of-basis matrix satisfies $P^{-1} = P^\top$, so (3.6) becomes $P^\top AP$, which is what Example 4.4 used.
- **LA-6.** The determinant is (3.8) generalized: the signed volume scaling of a map. $\det(AB) = \det A\det B$ is Theorem 3.7 read through volumes, and $\det(P^{-1}AP) = \det A$ makes "the determinant of an operator" meaningful, as Theorem 3.10 did for the trace.
- **LA-7.** Diagonalization is Theorem 3.8 with $P$ chosen so that $P^{-1}AP$ is diagonal; the columns of $P$ are the eigenvectors. The eigenvalue invariance in §3.7 is the first line of that module. The shift matrix of Example 4.1 is a Jordan block, the building block of LA-11.
- **LA-9.** The image of the unit circle under any $2\times 2$ map is an ellipse (`figures/unit_circle_images.png`); its semi-axes are the singular values and their directions come from the SVD. Rank is the number of nonzero singular values.
- **LA-12.** $\mathcal{L}(V, \mathbb{R}) = V^*$ is the dual space (Theorem 3.19 with $m = 1$). "A tensor is an array that transforms like $P^{-1}AP$" is the physicist's definition; LA-12 proves it equivalent to the multilinear-map definition.
- **CM-7.** The inertia tensor is a symmetric operator on $\mathbb{R}^3$; the principal axes are the basis in which its matrix is diagonal. Example 4.4 is the two-dimensional version with stress instead of inertia.
- **QM-4.** An observable is an operator; its matrix in the energy eigenbasis is diagonal; changing basis from position to momentum representation is a change of basis in the sense of (3.5), with $P$ unitary. The commutator $X \mapsto MX - XM$ of §3.1(h) is the map whose kernel is "observables that commute with $M$."
- **ML-8, ML-9.** A linear layer in a neural network is a linear map $\mathbb{R}^n \to \mathbb{R}^m$ plus a shift; composing layers is (3.4); the Jacobian of a composition (LA-14, backpropagation) is a product of matrices for the same reason.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 drill the definition of linearity, computing matrices from (3.2), coordinate vectors, kernels, and the rotation and differentiation matrices. B1–B8 are the module: B1 (isomorphisms carry bases to bases), B2 (integration and differentiation, when a one-sided inverse exists), B3 (change of basis in both directions, numbers), B5 ($AB = I$ implies $BA = I$, via Rank–Nullity), and B8 (extend `project.py`) are essential. C1–C3 go further: rank of a composition, nilpotent index, and the structure theorem for projections.

## 7. Coding project

`project.py` builds and checks everything in §3–4:

1. **`coords(v, basis)` and `matrix_of_map(T, basis_in, basis_out)`**, Definitions 3.3 and 3.5 implemented literally: solve for coordinates, build column $j$ from $T(\mathbf{b}_j)$. $T$ is any Python function on vectors. Polynomials are stored as coefficient vectors, matrices flattened, so a basis of $\mathcal{P}_3$ or of $\mathbb{R}^{2\times 2}$ is a list of length-4 arrays (Theorem 3.17b in code).
2. **The differentiation matrix** in the monomial and factorial bases (Example 4.1), checked against an independent SymPy construction, and $D^4 = 0 \neq D^3$.
3. **Theorems 3.6 and 3.7** checked on random polynomials, on $D\circ D$, and on $D\circ J$ versus $J\circ D$ (Exercise B2).
4. **`change_of_basis(P, A)`** returning $P^{-1}AP$; Example 4.2 computed and compared with the direct definition; trace, determinant, and eigenvalues shown invariant on random $P$ up to $5\times 5$.
5. **Projection onto $y = 2x$** (Example 4.3): $P^2 = P$, kernel, image, $\operatorname{diag}(1, 0)$ in the adapted basis.
6. **Grid deformations** for shear, rotation, projection, and a random map; the image of the unit circle (an ellipse; the code labels its semi-axes with `np.linalg.svd`, a preview of LA-9 and not something derived here); area of the image of the unit square versus $|ad - bc|$ by the shoelace formula.
7. **Rank–Nullity on 300 random maps** of chosen rank, with `np.linalg.matrix_rank` and `scipy.linalg.null_space`.
8. **The stress tensor in a rotated frame** (Example 4.4) and Mohr's circle.

**Build it yourself:** before reading the code, write `coords` (solve $B\mathbf{c} = \mathbf{v}$ with the basis vectors as the columns of $B$) and `matrix_of_map` (a loop over `basis_in`, one `coords` call per column). Then compare. Then do the "Your turn" cells: the $4\times 4$ matrix of the transpose map on $2\times 2$ matrices, the kernel and image of $A \mapsto A - A^\top$, rotation about the $z$-axis in a tilted basis, and a numerical linearity test. Reference checks for the first three sit at the bottom of the file.

## 8. Where you'll struggle

- **Which way does $P$ go.** This is the single most common error in the whole course, and you will make it. The convention here: columns of $P$ are the *new* basis vectors in *old* coordinates, so $P$ eats new coordinates and returns old ones, and $[T]_{\text{new}} = P^{-1}[T]_{\text{old}}P$. Some books define $P$ the other way around and write $PAP^{-1}$; both are correct *with their own definition of $P$*. The fix is not to memorize but to test: **apply $P$ to $\mathbf{e}_1$.** If you get the first new basis vector written in old coordinates, $P$ converts new to old. If you get nonsense, you are holding $P^{-1}$. Do this every single time, including in Example 4.4 where the convention $P = R_\theta$ versus $P = R_\theta^\top$ is the difference between $+\theta$ and $-\theta$ in the Mohr's circle formulas.
- **Active versus passive.** $R_\theta$ rotates vectors by $\theta$ (active). The same matrix, used as $P$, describes new axes rotated by $\theta$, and then the *components* of a fixed vector in the new frame are $P^{-1}\mathbf{v} = R_{-\theta}\mathbf{v}$: rotate the axes one way and the coordinates turn the other way. Physicists switch between these views constantly without saying so. When a sign is wrong, this is usually why.
- **The map is not the matrix.** $T$ is a function. $[T]_{\mathcal{B}}$ is a table of numbers that depends on $\mathcal{B}$. Writing "$T = A$" is only legal when $V = \mathbb{R}^n$ with the standard basis, where Definition 3.5 collapses to LA-1 Theorem 3.5. The statement "the trace of $T$" is meaningful only because of Theorem 3.10; the statement "the $(1,2)$ entry of $T$" is meaningless. If you find yourself asking whether a *map* is symmetric or upper triangular, stop: those are properties of a matrix in a basis, and (except for symmetry under orthogonal changes of basis, LA-5) they do not survive $P^{-1}AP$.
- **Thinking Rank–Nullity is about matrices.** It is about linear maps. The matrix version ($n = $ nullity $+$ rank, from counting pivots) is the shadow of (3.7) in one basis. The map version applies to $D$ on polynomials, to $A \mapsto A - A^\top$, to evaluation at a point, to the Laplacian on a finite-dimensional space of functions, with no coordinates in sight. When you reach for pivots to answer a question about $D$, you have missed the point of the module.
- **The proof of Rank–Nullity feels circular.** "Extend a basis of the kernel to a basis of $V$": how do you know you can? Because of LA-3's basis extension theorem, which was proved by a different argument (adding vectors not in the span until you cannot, and showing this stops after at most $\dim V$ steps). The proof here is *layered*, not circular: it uses a previously proven fact about bases to prove a new fact about maps. The second confusion is in step 4: "the $\mathbf{v}_j$ are independent, so the $T\mathbf{v}_j$ are independent" is *false* in general (the kernel is precisely the collection of independence-destroying directions), and the proof never says it. It says: a relation among the $T\mathbf{v}_j$ pulls back to a vector in the kernel, and the $\mathbf{v}_j$ were chosen independent *of the kernel's basis*. Reread step 4 until you can reproduce it from the sentence "pull the relation back into $V$."
- **Coordinate vectors in $\mathbb{R}^n$.** Because $[\mathbf{v}]_{\text{std}} = \mathbf{v}$, you have never had to distinguish a vector from its coordinates. The moment the basis is $(1, 1), (1, -1)$ the two come apart, and most errors in B3 are from adding a vector to a coordinate vector. Write the subscript $\mathcal{B}$ every time until the distinction is reflexive.
- **Injective iff surjective only when dimensions agree.** Corollary 3.15 is a statement about *square* situations. Students apply it to $D: \mathcal{P}_3 \to \mathcal{P}_2$ and conclude it must be injective. It is not.
- **Exercises that are harder than they look.** B5 ($AB = I \Rightarrow BA = I$) is three lines, but the three lines use Theorem 3.13, Corollary 3.15, and associativity in a specific order; most first attempts assume what they are proving. C2 needs a chain $\ker N \subsetneq \ker N^2 \subsetneq \cdots$ and the observation that strictly increasing subspaces have strictly increasing dimensions. "Your turn" item 3 in the project is mostly bookkeeping, but the non-orthogonal tilt will give you a non-orthogonal matrix for a rotation, which is correct and surprising: "orthogonal matrix" is a property of a matrix in an *orthonormal* basis, not of the map.

## 9. Check yourself

1. State the definition of a linear map. Give one linear map on $\mathcal{P}_2$ and one non-linear map on $\mathcal{P}_2$, each with a one-line justification.
2. Why is a linear map determined by its values on a basis? Which two properties of a basis are used, and where?
3. What is column $j$ of $[T]_{\mathcal{B},\mathcal{C}}$, in words? What is its shape?
4. Write the formula $[T\mathbf{v}]_{\mathcal{C}} = [T]_{\mathcal{B},\mathcal{C}}[\mathbf{v}]_{\mathcal{B}}$ and derive it in four lines.
5. Why is matrix multiplication defined by row-times-column? Answer with Theorem 3.7.
6. Let $\mathcal{B}' = \big((2, 1), (1, 1)\big)$ in $\mathbb{R}^2$. Write $P$. Which way does it convert? Find $[\,(3, 2)\,]_{\mathcal{B}'}$.
7. Derive $A' = P^{-1}AP$ by inserting the identity twice.
8. Prove $\operatorname{tr}(XY) = \operatorname{tr}(YX)$. Why does this make "the trace of an operator" well defined?
9. State and prove: $T$ is injective iff $\ker T = \{\mathbf{0}\}$.
10. State Rank–Nullity and sketch its proof. Where does independence of the $T\mathbf{v}_j$ come from?
11. $D: \mathcal{P}_4 \to \mathcal{P}_4$: what are the kernel, image, rank, nullity, and the smallest $k$ with $D^k = 0$?
12. List six conditions equivalent to "$A$ is invertible" for a square matrix, and for two of them say which theorem connects them.
13. What is $\dim\mathcal{L}(\mathcal{P}_3, \mathbb{R}^{2\times 2})$?

## 10. Pair with

- **Axler, *Linear Algebra Done Right*, 4th ed., ch. 3** (§3A–3D; §3E–3F are LA-12). Axler proves Theorems 3.2, 3.14, 3.17 and 3.19 in the same order as here, with cleaner notation and more examples of maps on function spaces. His proof of the Fundamental Theorem of Linear Maps is the one in §3.9. Read 3D for matrices of maps after you have done Example 4.1 by hand.
- **Strang, *Introduction to Linear Algebra*, §8.2–8.3** (the matrix of a linear transformation; the search for a good basis). Strang's §8.2 does the differentiation and integration matrices (Example 4.1 and Exercise B2) and is better than these notes on the "a different basis gives a different matrix" examples. §8.3 previews why diagonal (eigenvector) and orthogonal (SVD) bases are the ones worth searching for.
- **3Blue1Brown, *Essence of Linear Algebra*, episodes 9–13.** Episode 13 (change of basis) is the one to watch twice; it animates (3.5) and $P^{-1}AP$ and uses the same "whose grid are we speaking" language as §2. Episodes 9–12 (dot products and duality, cross products, Cramer's rule) are the bridge toward LA-6 and LA-12; watch them once now for the pictures. Rewatch episodes 3, 4, 7 and 8 (transformations, composition, null space and column space, nonsquare matrices) if the grid-deformation picture of §2 is not yet automatic.
