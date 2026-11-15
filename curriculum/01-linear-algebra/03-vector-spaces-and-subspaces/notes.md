# LA-3: Vector Spaces, Subspaces, Independence, Basis, Dimension, and the Four Fundamental Subspaces

**Time:** ~12 hours · **Prerequisites:** LA-2 (elimination, echelon forms, rank), F1 (proof writing) · **Unlocks:** LA-4 (linear maps), and through it LA-5, LA-6, LA-7

> **Conventions in this module.** Vectors are columns, bold lowercase: $\mathbf{v}$. Scalars are real unless said otherwise. Replacing $\mathbb{R}$ by $\mathbb{C}$ changes nothing in the definitions, the proofs, or the dimension counts (LA-7 needs that), with one exception: the statements about *orthogonality*. Over $\mathbb{C}$, reading (3.6) as "$\perp$", the LA-5 preview after it, and exercise B7 are false until $A^\top$ is replaced by the conjugate transpose $A^\dagger$ (LA-8). Example: $A = (1, i)^\top$ has $A^\top A = 1 + i^2 = 0$, so $N(A^\top A) = \mathbb{C}$ while $N(A) = \{\mathbf{0}\}$. $A$ is $m \times n$ ($m$ rows, $n$ columns) with columns $\mathbf{a}_1, \dots, \mathbf{a}_n \in \mathbb{R}^m$. $R$ is the reduced row echelon form (RREF) of $A$ from LA-2 (it is unique: LA-2 §3.9, and the corollary after Theorem 3.28 here), $r$ is its number of pivots, and $E$ is the invertible $m \times m$ matrix that records the row operations, so $EA = R$. Pivot columns are indexed $p_1 < \dots < p_r$, free columns $f_1 < \dots < f_{n-r}$. The four subspaces are $C(A)$ (column space), $N(A)$ (null space), $C(A^\top)$ (row space), $N(A^\top)$ (left null space). "Independent" always means linearly independent. A **basis is an ordered list** (so coordinates make sense); when a theorem says "set of vectors" it means a list without repeats. $\mathcal{P}_n$ is polynomials of degree $\le n$; $\mathcal{M}_{m\times n}$ is real $m \times n$ matrices; $C^\infty(\mathbb{R})$ is infinitely differentiable real functions. $\mathbf{0}$ is the zero vector of whatever space is in play; $0$ is the number. $\dim$ is dimension.

## 0. Where this sits

You can now solve $A\mathbf{x} = \mathbf{b}$ by elimination (LA-2) and you know that $A\mathbf{x}$ is a combination of the columns of $A$ (LA-1). This module asks the questions that elimination answers without saying so: *what is the set of all solutions, what shape does it have, and how big is it?* The answers (subspace, basis, dimension) are the vocabulary of every later linear algebra module. LA-4 reads this module's $C(A)$ and $N(A)$ as the image and kernel of a linear map; LA-5 adds right angles to the four subspaces; LA-9 gives them their best bases.

This is also the module where you stop computing and start proving. The roadmap's warning is accurate: LA-3 and LA-4 are the wall. Section 3.0 tells you how to climb it.

## 1. The problem that forced this idea

By 1850, three communities were using the same algebra without knowing it.

Geometers had arrows. You add arrows head to tail and stretch them by numbers. Grassmann's *Ausdehnungslehre* (1844) tried to build geometry in any number of dimensions on those two operations alone. The book was nearly unreadable and, for decades, nearly unread.

Analysts had functions. Eighteenth-century work on linear differential equations (Euler, d'Alembert, Lagrange) leaned on a fact everyone used and nobody named: if $y_1$ and $y_2$ solve $y'' + y = 0$, so does $c_1 y_1 + c_2 y_2$. Fourier (1807; in print 1822) went further and claimed a function on an interval could be written as $\sum_k a_k \cos kx + b_k \sin kx$, with coefficients computed by a formula that reads exactly like "dot the function against $\cos kx$ and divide by the squared length." He was treating functions as arrows with infinitely many components, without the vocabulary to say so.

Algebraists had $n$-tuples and matrices (Cayley, 1858; LA-1 §1). And the same facts kept being reproved in each setting. "A solution of the inhomogeneous equation plus any solution of the homogeneous equation is again a solution." "If a combination of these objects is zero then every coefficient is zero, so the coefficients of any expansion are unique." "Two of them suffice and a third is redundant." Each community proved these for its own objects, and the proofs never used what the objects *were*, only how they added and scaled.

Giuseppe Peano noticed. In 1888, in a short book whose title credits Grassmann (*Calcolo geometrico secondo l'Ausdehnungslehre di H. Grassmann*), he wrote down axioms for what he called a *sistema lineare*: a set of things with an addition and a scaling obeying a handful of rules, and nothing else assumed. He gave polynomials and functions as examples and remarked that such a system could have infinitely many dimensions. His list is, up to small changes, the eight axioms of §3.1. Almost nobody read that either.

Physics was one of the forces that made the abstraction unavoidable. Hilbert's work on integral equations (1904–1910) traded a function for its infinite sequence of coefficients against an orthonormal family of functions, and Erhard Schmidt (1908) worked out the geometry (lengths, angles, orthogonal projection) of the resulting space $\ell^2$ of square-summable sequences. Then quantum mechanics arrived twice within a year. Heisenberg's 1925 theory, which Born and Jordan recast within months as a calculus of infinite matrices, worked with arrays of numbers; Schrödinger's 1926 theory used a differential equation for a wave function. Schrödinger, and independently Carl Eckart, showed in 1926 that the two gave the same answers, but the arguments were formal, and it was murky what single object both theories were describing. Von Neumann's answer (1927, then the 1932 *Mathematische Grundlagen der Quantenmechanik*) was one abstract space written in two coordinate systems: Heisenberg's vectors were sequences of components, Schrödinger's were functions. One space, two bases. That sentence is only *sayable* if "vector space" does not care whether its elements are sequences or functions, and "basis" lets you change coordinates between them. Weyl's *Raum, Zeit, Materie* (1918) and Banach's thesis (1922) had already put the axiomatic definition in front of physicists and analysts. Its route into the standard curriculum ran at least as much through algebra (van der Waerden's *Moderne Algebra*, 1930–31) and functional analysis (Banach's *Théorie des opérations linéaires*, 1932) as through physics; G. H. Moore, "The axiomatization of linear algebra: 1875–1940" (*Historia Mathematica*, 1995), traces it.

So here is the deal this module offers. You already know $\mathbb{R}^n$. You will now learn to prove things about $\mathbb{R}^n$ *using only the eight rules of addition and scaling*. Every such proof is then automatically a proof about polynomials, matrices, solutions of your ODEs, and wave functions. The price is that you have to write proofs, not check examples. That is what the twelve hours are for.

## 2. Intuition

### 2.1 Subspaces: flat things through the origin

In $\mathbb{R}^3$, which sets of vectors are closed off under the two operations? Pick a set $W$ and ask two questions. If I add two vectors in $W$, am I still in $W$? If I scale a vector in $W$, am I still in $W$?

A line through the origin: yes and yes. A plane through the origin: yes and yes. All of $\mathbb{R}^3$: trivially. The single point $\{\mathbf{0}\}$: yes, since $\mathbf{0} + \mathbf{0} = \mathbf{0}$ and $c\mathbf{0} = \mathbf{0}$.

A line *not* through the origin: no. Take a point $\mathbf{p}$ on it and scale by $0$; you land at the origin, which is not on the line. A plane missing the origin fails the same way. The unit sphere fails (add two points on it). The first quadrant of $\mathbb{R}^2$ fails (multiply by $-1$). The union of two lines through the origin fails (add one vector from each; the sum is on neither line).

```
   y                                   y
   |        /  line y = x              |          /  line y = x + 1
   |      /    through 0               |        /    misses 0
   |    /                              |      /
---+--/--------- x                  ---+----/------- x
   |/                                  |  /
  /|   v, 2v, -v, v+v all on it       /|   0 is not on it; 0*v leaves it
```

The flat sets through the origin are the **subspaces**. Carry this guess forward: *the subspaces of $\mathbb{R}^3$ are exactly $\{\mathbf{0}\}$, the lines through $\mathbf{0}$, the planes through $\mathbf{0}$, and $\mathbb{R}^3$ itself.* Nothing curved, nothing shifted, nothing with a hole. Corollary 3.26 proves it.

There are two ways to manufacture a subspace, and they are the two pictures of LA-1. *Column picture:* span some vectors ($\operatorname{span}\{\mathbf{v}_1, \mathbf{v}_2\}$ is a plane). *Row picture:* impose homogeneous linear equations ($\{\mathbf{x} : \mathbf{n}\cdot\mathbf{x} = 0\}$ is a plane). Every subspace of $\mathbb{R}^n$ has both descriptions (exercise B5 proves it for planes in $\mathbb{R}^3$, by a dimension count that works in every $\mathbb{R}^n$). Converting between them is the row space / null space story of §3.8–§3.9, and LA-5 finishes it with orthogonal complements (LA-5 Theorem 3.7).

### 2.2 Independence: no redundant directions

Start in $\mathbb{R}^2$ with one nonzero vector $\mathbf{v}$. Its span is a line. Add a second vector $\mathbf{w}$. Either $\mathbf{w}$ lies on that line (it adds nothing; it is *redundant*) or it does not (the span is now the whole plane). Add a third vector $\mathbf{u}$. It is in the plane, so it is in the span of the first two already: redundant, no matter what it is. In $\mathbb{R}^2$ you can have at most two non-redundant directions. In $\mathbb{R}^3$, three.

"Redundant" is symmetric. If $\mathbf{u} = \mathbf{v} + \mathbf{w}$ you could equally throw out $\mathbf{v}$, since $\mathbf{v} = \mathbf{u} - \mathbf{w}$. What is intrinsic is not *which* vector is redundant but the existence of a *relation* $c_1\mathbf{v} + c_2\mathbf{w} + c_3\mathbf{u} = \mathbf{0}$ with coefficients not all zero. A list with no such relation is **independent**: the only way to combine its vectors to get $\mathbf{0}$ is with all coefficients zero.

The number of non-redundant directions needed to build a flat set is its **dimension**. The guess that matters: *it does not depend on which directions you choose.* Two people building the same plane with different pairs of vectors both use exactly two. Theorem 3.23 proves it, and the proof is the hardest in the module.

### 2.3 Same algebra, different costume

Polynomials of degree at most two: $a + bx + cx^2 \leftrightarrow (a, b, c)$. Adding polynomials adds the coefficient lists; scaling a polynomial scales the list. So $\mathcal{P}_2$ is $\mathbb{R}^3$ in disguise, and $1, x, x^2$ play the role of $\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3$.

Solutions of $y'' + y = 0$: every one is $a\cos t + b\sin t$ for exactly one pair $(a, b)$, namely $(y(0), y'(0))$. The solution set is $\mathbb{R}^2$ in disguise, with $\cos t$ and $\sin t$ as the axes.

Symmetric $2 \times 2$ matrices $\begin{pmatrix} a & b \\ b & c\end{pmatrix} \leftrightarrow (a, b, c)$. Again $\mathbb{R}^3$.

What survives the change of costume is the number of free parameters, the rule for adding, and the rule for scaling. That is all a vector space is. So the guesses: **dimension** is the number of free parameters, a **basis** is one choice of coordinate axes, and **coordinates** are the parameters themselves.

### 2.4 The two questions about $A\mathbf{x} = \mathbf{b}$, and the count

LA-1 Theorem 3.6: $A\mathbf{x} = \mathbf{b}$ is solvable iff $\mathbf{b}$ is in the span of the columns. LA-1 exercise B3: if $A\mathbf{z} = \mathbf{0}$ forces $\mathbf{z} = \mathbf{0}$, then a solution, when one exists, is unique (Corollary 3.10 below gives the converse). These are different questions. *Existence* depends on $\mathbf{b}$ and the columns. *Uniqueness* depends on the columns alone.

Picture a $3 \times 4$ matrix of rank $2$. Four columns in $\mathbb{R}^3$, but elimination found only two pivots, so only two non-redundant directions: the column space is a plane in $\mathbb{R}^3$, dimension $2 = r$. The two redundant columns each produce a relation $A\mathbf{z} = \mathbf{0}$, and those relations are themselves non-redundant: the null space is a plane in $\mathbb{R}^4$, dimension $2 = n - r$. Now the rows. Three rows in $\mathbb{R}^4$, but elimination turned one of them into zeros, so two non-redundant rows: the row space has dimension $2 = r$. The one way of combining rows to get zero gives a vector $\mathbf{y}$ with $\mathbf{y}^\top A = \mathbf{0}$: the left null space is a line in $\mathbb{R}^3$, dimension $1 = m - r$.

You should now be able to guess the main theorem of the module (Theorem 3.30): the column space and row space both have dimension $r$; the null space has dimension $n - r$; the left null space has dimension $m - r$. On the input side $r + (n - r) = n$. On the output side $r + (m - r) = m$. The same $r$ on both sides is the statement "row rank equals column rank," and it is not obvious.

## 3. Formal development

### 3.0 How to read and write a proof in this module

> Every theorem below has the form "if [hypotheses] then [conclusion]." A proof is a chain of sentences, each of which is one of: (a) a hypothesis restated, (b) a definition unpacked, (c) an axiom or earlier theorem applied, (d) a logical step from earlier sentences. Nothing else is admissible. Not "it's clear from the picture." Not "it worked for the example." When you *read* a proof, label every sentence with its letter. When you *write* one, put down the first line ("Let ... be arbitrary with [hypotheses]") and the last line ("Therefore [conclusion]") before anything else, then fill the middle.
>
> Three patterns cover most of this module:
>
> - **To show a set is a subspace:** take two arbitrary elements and an arbitrary scalar; show the sum and the multiple stay in the set; check the set is nonempty. (Theorem 3.4 says that is enough.)
> - **To show two sets are equal:** show each contains the other. Two paragraphs, labelled $(\subseteq)$ and $(\supseteq)$.
> - **To show something is unique:** suppose there are two, subtract, show the difference is zero.
>
> Two fully annotated proofs follow. The bracketed labels are not part of the proof; they are what should be going on in your head.

**Example proof A (pattern: subspace).** *Claim.* If $U$ and $W$ are subspaces of a vector space $V$, then $U \cap W$ is a subspace of $V$.

*Proof.* [Set-up: name arbitrary elements satisfying the hypothesis.] Let $\mathbf{x}, \mathbf{y} \in U \cap W$ and let $c \in \mathbb{R}$. [Unpack the definition of $\cap$.] Then $\mathbf{x}, \mathbf{y} \in U$ and $\mathbf{x}, \mathbf{y} \in W$. [Apply the hypothesis "$U$ is a subspace": closure.] Since $U$ is a subspace, $\mathbf{x} + \mathbf{y} \in U$ and $c\mathbf{x} \in U$. [Same hypothesis for $W$.] Since $W$ is a subspace, $\mathbf{x} + \mathbf{y} \in W$ and $c\mathbf{x} \in W$. [Re-pack the definition of $\cap$.] So $\mathbf{x} + \mathbf{y} \in U \cap W$ and $c\mathbf{x} \in U \cap W$. [Nonempty.] Every subspace contains $\mathbf{0}$ (Theorem 3.4), so $\mathbf{0} \in U$ and $\mathbf{0} \in W$, hence $\mathbf{0} \in U \cap W$. [Invoke the test.] By Theorem 3.4, $U \cap W$ is a subspace of $V$. ∎

Notice what the proof did *not* do. It did not take $V = \mathbb{R}^3$ and two planes. It works for subspaces of function spaces you will never draw. That is the entire value of writing it this way.

**Example proof B (pattern: uniqueness, inside an "if and only if").** *Claim.* Let $\mathbf{v}_1, \dots, \mathbf{v}_k$ be vectors in $V$. The list is linearly independent (Definition 3.11) if and only if every vector in $\operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_k\}$ can be written as $c_1\mathbf{v}_1 + \dots + c_k\mathbf{v}_k$ in exactly one way.

*Proof.* [An "iff" is two implications. Label them and do them separately.]

$(\Rightarrow)$ [Hypothesis: independent. Goal: unique coefficients.] Suppose the list is independent, and suppose some $\mathbf{w}$ in the span has two representations, $\mathbf{w} = \sum_i c_i\mathbf{v}_i$ and $\mathbf{w} = \sum_i d_i\mathbf{v}_i$. [The uniqueness move: subtract.] Subtracting, $\sum_i (c_i - d_i)\mathbf{v}_i = \mathbf{w} - \mathbf{w} = \mathbf{0}$. [Apply the hypothesis to this relation.] By independence, every coefficient $c_i - d_i$ is $0$, so $c_i = d_i$ for all $i$. [Conclude.] The two representations were the same one.

$(\Leftarrow)$ [Hypothesis: unique coefficients. Goal: the only relation is the trivial one.] Suppose every vector in the span has a unique representation, and suppose $\sum_i c_i\mathbf{v}_i = \mathbf{0}$. [Find two representations of one vector.] The vector $\mathbf{0}$ is in the span, and $\mathbf{0} = \sum_i 0\cdot\mathbf{v}_i$ is one representation of it; $\sum_i c_i\mathbf{v}_i$ is another. [Apply the hypothesis.] By uniqueness, $c_i = 0$ for every $i$. [Conclude.] That is the definition of independence. ∎

This claim is Theorem 3.15 below and is the most important sentence in the module: **independence is uniqueness of coefficients; span is existence of coefficients.** Keep them apart.

### 3.1 Vector spaces

**Definition 3.1 (Vector space).** A real vector space is a set $V$ with two operations, addition ($\mathbf{u} + \mathbf{v} \in V$ for $\mathbf{u}, \mathbf{v} \in V$) and scalar multiplication ($c\mathbf{v} \in V$ for $c \in \mathbb{R}$, $\mathbf{v} \in V$), such that for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and all $c, d \in \mathbb{R}$:

| | Axiom | Name |
|---|---|---|
| (A1) | $\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}$ | commutativity |
| (A2) | $(\mathbf{u} + \mathbf{v}) + \mathbf{w} = \mathbf{u} + (\mathbf{v} + \mathbf{w})$ | associativity |
| (A3) | there is a $\mathbf{0} \in V$ with $\mathbf{v} + \mathbf{0} = \mathbf{v}$ for all $\mathbf{v}$ | zero vector |
| (A4) | for each $\mathbf{v}$ there is a $-\mathbf{v} \in V$ with $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$ | negatives |
| (M1) | $c(\mathbf{u} + \mathbf{v}) = c\mathbf{u} + c\mathbf{v}$ | distributivity over vector addition |
| (M2) | $(c + d)\mathbf{v} = c\mathbf{v} + d\mathbf{v}$ | distributivity over scalar addition |
| (M3) | $(cd)\mathbf{v} = c(d\mathbf{v})$ | compatibility |
| (M4) | $1\mathbf{v} = \mathbf{v}$ | unit |

The elements of $V$ are called vectors whatever they actually are. The requirement that $\mathbf{u} + \mathbf{v}$ and $c\mathbf{v}$ land back in $V$ (closure) is part of the definition.

**Examples.** In each, all eight axioms are verified by one move: compute both sides entry by entry (or point by point) and observe that each entry reduces to the corresponding rule for real numbers.

1. $\mathbb{R}^n$ with entrywise operations (LA-1 Definition 3.1).
2. $\mathcal{M}_{m\times n}$, real $m \times n$ matrices, entrywise. This is $\mathbb{R}^{mn}$ with the entries arranged in a rectangle.
3. $\mathcal{P}_n$, polynomials of degree $\le n$: $(p + q)(x) = p(x) + q(x)$, $(cp)(x) = c\,p(x)$. Degree $\le n$ is preserved because the degree of a sum is at most the larger degree. (Polynomials of degree *exactly* $n$ do not form a vector space: $x^n + (-x^n) = 0$ has no degree $n$.)
4. $C^\infty(\mathbb{R})$, infinitely differentiable functions, pointwise operations. Sums and multiples of smooth functions are smooth. Axiom (A1), for instance: $(f + g)(t) = f(t) + g(t) = g(t) + f(t) = (g + f)(t)$ for every $t$, so $f + g = g + f$. The same works for all functions on an interval, or all continuous ones.
5. $S = \{y \in C^\infty(\mathbb{R}) : y'' + y = 0\}$, the solutions of the harmonic oscillator. If $y_1, y_2 \in S$ and $c \in \mathbb{R}$, then $(y_1 + y_2)'' + (y_1 + y_2) = (y_1'' + y_1) + (y_2'' + y_2) = 0 + 0 = 0$ and $(cy_1)'' + cy_1 = c(y_1'' + y_1) = 0$. So $S$ is closed under both operations, and the axioms hold in $S$ because they hold in $C^\infty(\mathbb{R})$ (this is Theorem 3.4 in action). Example 4.5 shows $S$ is two-dimensional with basis $\cos t, \sin t$.
6. $\mathbb{R}^\infty$, all real sequences $(a_1, a_2, a_3, \dots)$, termwise. Infinite-dimensional (§3.6).

A non-example, to see the axioms have teeth: $\mathbb{R}^2$ with ordinary addition but scaling defined by $c \odot (x, y) = (cx, y)$. Then (M2) fails: $(1 + 1) \odot (x, y) = (2x, y)$ but $1 \odot (x, y) + 1 \odot (x, y) = (2x, 2y)$.

Anything true in *every* vector space must follow from the eight axioms alone. Here is the first such fact.

▶ **Try it first.** Show $0\mathbf{v} = \mathbf{0}$ for every $\mathbf{v}$. Hint: $0 = 0 + 0$.

**Theorem 3.2.** In any vector space: $0\mathbf{v} = \mathbf{0}$, $c\mathbf{0} = \mathbf{0}$, and $(-1)\mathbf{v} = -\mathbf{v}$.

**Derivation** (first statement; the other two are exercise B2).

1. $0\mathbf{v} = (0 + 0)\mathbf{v} = 0\mathbf{v} + 0\mathbf{v}$ by (M2).
2. By (A4), $0\mathbf{v}$ has a negative; add it to both sides: $0\mathbf{v} + (-(0\mathbf{v})) = (0\mathbf{v} + 0\mathbf{v}) + (-(0\mathbf{v}))$.
3. The left side is $\mathbf{0}$ by (A4). The right side, by (A2), is $0\mathbf{v} + \big(0\mathbf{v} + (-(0\mathbf{v}))\big) = 0\mathbf{v} + \mathbf{0} = 0\mathbf{v}$, using (A4) then (A3).
4. So $\mathbf{0} = 0\mathbf{v}$. ∎

■ **What just happened.** Nothing about lists of numbers was used. The result holds for polynomials, matrices, functions, sequences, and spaces nobody has thought of yet. That is the point of axioms: prove once.

### 3.2 Subspaces

**Definition 3.3 (Subspace).** A subset $W$ of a vector space $V$ is a **subspace** if $W$ is itself a vector space under the addition and scalar multiplication of $V$.

That reads like eight axioms to check. It is two.

**Theorem 3.4 (Two-closure test).** A subset $W \subseteq V$ is a subspace if and only if

$$
\text{(i) } W \ne \emptyset, \qquad \text{(ii) } \mathbf{u}, \mathbf{w} \in W \Rightarrow \mathbf{u} + \mathbf{w} \in W, \qquad \text{(iii) } \mathbf{w} \in W,\ c \in \mathbb{R} \Rightarrow c\mathbf{w} \in W. \tag{3.1}
$$

Moreover, every subspace contains $\mathbf{0}$, so (i) can be replaced by "$\mathbf{0} \in W$."

▶ **Try it first.** Which of the eight axioms are *identities* (true for all vectors in $V$, hence automatically for all vectors in $W$)? Which two are *existence* statements that need an argument?

**Derivation.**

1. $(\Rightarrow)$ If $W$ is a vector space under the inherited operations, closure (ii) and (iii) is part of Definition 3.1, and $W$ contains its zero vector so is nonempty.
2. $(\Leftarrow)$ Assume (3.1). Axioms (A1), (A2), (M1)–(M4) are equations that hold for *all* vectors of $V$; vectors of $W$ are vectors of $V$, so they hold in $W$. Closure makes the two operations well-defined on $W$.
3. (A3): $W$ is nonempty, so pick some $\mathbf{w} \in W$. By (iii) with $c = 0$, $0\mathbf{w} \in W$, and by Theorem 3.2, $0\mathbf{w} = \mathbf{0}$. So $\mathbf{0} \in W$ and it serves as the zero of $W$.
4. (A4): for $\mathbf{w} \in W$, (iii) with $c = -1$ gives $(-1)\mathbf{w} \in W$, and by Theorem 3.2 this is $-\mathbf{w}$.
5. All eight axioms hold in $W$; it is a vector space. Step 3 also proves the "moreover." ∎

■ **What just happened.** The axioms that are identities are inherited for free. The two existence axioms were manufactured by scaling: $\mathbf{0} = 0\mathbf{w}$ and $-\mathbf{w} = (-1)\mathbf{w}$. This is why §2.1's shifted line fails: it is not closed under scaling by $0$.

**Theorem 3.5 (Span is a subspace).** For any $\mathbf{v}_1, \dots, \mathbf{v}_k \in V$, the set $\operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_k\}$ is a subspace of $V$. It is the smallest one containing all the $\mathbf{v}_i$: any subspace containing $\mathbf{v}_1, \dots, \mathbf{v}_k$ contains their span.

*Proof.* Nonempty: $\mathbf{0} = \sum_i 0\mathbf{v}_i$ is in the span. Closure: if $\mathbf{p} = \sum_i a_i\mathbf{v}_i$ and $\mathbf{q} = \sum_i b_i\mathbf{v}_i$ then $\mathbf{p} + \mathbf{q} = \sum_i (a_i + b_i)\mathbf{v}_i$ by (A1), (A2), (M2), and $c\mathbf{p} = \sum_i (ca_i)\mathbf{v}_i$ by (M1), (M3). Both are in the span. By Theorem 3.4 the span is a subspace. Smallest: if $W$ is a subspace containing every $\mathbf{v}_i$, then by (iii) each $a_i\mathbf{v}_i \in W$, and by (ii) applied $k - 1$ times, $\sum_i a_i\mathbf{v}_i \in W$. ∎

**Theorem 3.6 (Column space and null space).** For an $m \times n$ matrix $A$:
- $C(A) = \{A\mathbf{x} : \mathbf{x} \in \mathbb{R}^n\} = \operatorname{span}\{\mathbf{a}_1, \dots, \mathbf{a}_n\}$ is a subspace of $\mathbb{R}^m$.
- $N(A) = \{\mathbf{x} \in \mathbb{R}^n : A\mathbf{x} = \mathbf{0}\}$ is a subspace of $\mathbb{R}^n$.

*Proof.* The two descriptions of $C(A)$ agree by LA-1 (3.1): $A\mathbf{x} = \sum_j x_j\mathbf{a}_j$ runs over exactly the combinations of the columns as $\mathbf{x}$ runs over $\mathbb{R}^n$. So $C(A)$ is a span, hence a subspace by Theorem 3.5. For $N(A)$: $A\mathbf{0} = \mathbf{0}$, so $\mathbf{0} \in N(A)$. If $A\mathbf{x} = \mathbf{0}$ and $A\mathbf{y} = \mathbf{0}$ then by linearity (LA-1 Theorem 3.4) $A(\mathbf{x} + \mathbf{y}) = A\mathbf{x} + A\mathbf{y} = \mathbf{0}$ and $A(c\mathbf{x}) = cA\mathbf{x} = \mathbf{0}$. Theorem 3.4 finishes. ∎

Note where each lives: $C(A) \subseteq \mathbb{R}^m$ (the space of right-hand sides $\mathbf{b}$), $N(A) \subseteq \mathbb{R}^n$ (the space of unknowns $\mathbf{x}$). Mixing these up is the most common error in the module.

### 3.3 Computing $N(A)$ from the RREF, and the complete solution

From LA-2: elimination produces an invertible $E$ with $EA = R$, where $R$ is in reduced row echelon form. Row $i \le r$ of $R$ has a $1$ in pivot column $p_i$, zeros in every other pivot column, zeros to the left of $p_i$, and arbitrary entries $R_{i,f_j}$ in the free columns to its right. Rows $r + 1, \dots, m$ of $R$ are zero.

**Lemma 3.7 (Row operations preserve the null space).** If $E$ is invertible, $N(EA) = N(A)$. In particular $N(R) = N(A)$.

*Proof.* $(\supseteq)$ If $A\mathbf{x} = \mathbf{0}$ then $EA\mathbf{x} = E\mathbf{0} = \mathbf{0}$. $(\subseteq)$ If $EA\mathbf{x} = \mathbf{0}$, multiply by $E^{-1}$: $A\mathbf{x} = \mathbf{0}$. ∎

Read $A\mathbf{x} = \mathbf{0}$ as $\sum_j x_j\mathbf{a}_j = \mathbf{0}$: a vector in $N(A)$ is a *relation among the columns*. The lemma says $A$ and $R$ have identical column relations. We use this repeatedly.

**Definition (Special solutions).** For each free column $f_j$ ($j = 1, \dots, n - r$) the special solution $\mathbf{s}_j \in \mathbb{R}^n$ has $(\mathbf{s}_j)_{f_j} = 1$, $(\mathbf{s}_j)_{f_k} = 0$ for $k \ne j$, and pivot entries $(\mathbf{s}_j)_{p_i} = -R_{i,f_j}$ for $i = 1, \dots, r$.

▶ **Try it first.** Write out what row $i$ of $R\mathbf{x} = \mathbf{0}$ says, solving for the pivot variable $x_{p_i}$ in terms of the free variables. Then check that $\mathbf{s}_j$ satisfies every row.

**Theorem 3.8 (Special solutions describe the null space).** Every $\mathbf{s}_j \in N(A)$, and every $\mathbf{x} \in N(A)$ is

$$
\mathbf{x} = \sum_{j=1}^{n-r} x_{f_j}\,\mathbf{s}_j, \tag{3.2}
$$

with coefficients equal to its own free entries.

**Derivation.**

1. Row $i \le r$ of $R\mathbf{x} = \mathbf{0}$ reads $x_{p_i} + \sum_j R_{i,f_j}x_{f_j} = 0$ (the $1$ in column $p_i$, zeros in other pivot columns, zeros to the left). So $\mathbf{x} \in N(R)$ iff for every $i$: $x_{p_i} = -\sum_j R_{i,f_j}x_{f_j}$. The free entries can be anything; the pivot entries are then forced.
2. $\mathbf{s}_j$ has free entries $(0, \dots, 1, \dots, 0)$ (the $1$ in slot $j$), so step 1 demands $x_{p_i} = -R_{i,f_j}$, which is exactly how $\mathbf{s}_j$ was defined. So $R\mathbf{s}_j = \mathbf{0}$, and by Lemma 3.7, $\mathbf{s}_j \in N(A)$.
3. Let $\mathbf{x} \in N(A) = N(R)$ and set $\mathbf{z} = \mathbf{x} - \sum_j x_{f_j}\mathbf{s}_j$. Since $N(R)$ is a subspace, $\mathbf{z} \in N(R)$. Its free entry $f_k$ is $x_{f_k} - \sum_j x_{f_j}(\mathbf{s}_j)_{f_k} = x_{f_k} - x_{f_k} = 0$.
4. A vector in $N(R)$ whose free entries are all zero has, by step 1, pivot entries $-\sum_j R_{i,f_j}\cdot 0 = 0$. So $\mathbf{z} = \mathbf{0}$, which is (3.2). ∎

■ **What just happened.** The null space has exactly $n - r$ free parameters, one per free column. Pivot variables are not parameters; they are determined. Theorem 3.18 (§3.5) shows the $\mathbf{s}_j$ are independent, so they are a basis and $\dim N(A) = n - r$.

**Theorem 3.9 (Complete solution).** If $A\mathbf{x}_p = \mathbf{b}$, then

$$
\{\mathbf{x} : A\mathbf{x} = \mathbf{b}\} = \{\mathbf{x}_p + \mathbf{x}_n : \mathbf{x}_n \in N(A)\}. \tag{3.3}
$$

*Proof.* $(\supseteq)$ $A(\mathbf{x}_p + \mathbf{x}_n) = A\mathbf{x}_p + A\mathbf{x}_n = \mathbf{b} + \mathbf{0} = \mathbf{b}$. $(\subseteq)$ Let $A\mathbf{x} = \mathbf{b}$ and put $\mathbf{x}_n = \mathbf{x} - \mathbf{x}_p$. Then $A\mathbf{x}_n = \mathbf{b} - \mathbf{b} = \mathbf{0}$, so $\mathbf{x}_n \in N(A)$ and $\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n$. ∎

To find $\mathbf{x}_p$: reduce $[A \mid \mathbf{b}]$ to $[R \mid E\mathbf{b}]$. The system is consistent iff the last $m - r$ entries of $E\mathbf{b}$ are zero; then set all free variables to $0$ and read $x_{p_i} = (E\mathbf{b})_i$.

**Corollary 3.10.** $A\mathbf{x} = \mathbf{b}$ has $0$, $1$, or infinitely many solutions. It has at most one solution (for every $\mathbf{b}$) iff $N(A) = \{\mathbf{0}\}$ iff $r = n$ iff there are no free columns.

*Proof.* If no $\mathbf{x}_p$ exists: $0$ solutions. If $\mathbf{x}_p$ exists and $N(A) = \{\mathbf{0}\}$, (3.3) gives exactly $\mathbf{x}_p$. If $\mathbf{x}_p$ exists and some $\mathbf{x}_n \ne \mathbf{0}$ is in $N(A)$, then $\mathbf{x}_p + c\mathbf{x}_n$ are infinitely many distinct solutions. So $N(A) = \{\mathbf{0}\}$ gives at most one solution for every $\mathbf{b}$. Conversely, if $N(A) \ne \{\mathbf{0}\}$, take $\mathbf{b} = \mathbf{0}$: then $\mathbf{x}_p = \mathbf{0}$ exists and the third case gives infinitely many solutions. Finally $N(A) = \{\mathbf{0}\}$ iff there is no special solution iff no free column, by Theorem 3.8 (if $f_1$ exists then $\mathbf{s}_1 \ne \mathbf{0}$). ∎

The solution set (3.3) is a subspace *shifted by* $\mathbf{x}_p$: the shifted line of §2.1. It is a subspace only when $\mathbf{b} = \mathbf{0}$. Your ODE course called the same structure "general solution = particular + homogeneous" (Example 4.5).

### 3.4 Linear independence

**Definition 3.11 (Linear independence).** A list $\mathbf{v}_1, \dots, \mathbf{v}_k$ in $V$ is **linearly independent** if

$$
c_1\mathbf{v}_1 + \dots + c_k\mathbf{v}_k = \mathbf{0} \quad\Longrightarrow\quad c_1 = \dots = c_k = 0. \tag{3.4}
$$

Otherwise it is **dependent**: some relation $\sum_i c_i\mathbf{v}_i = \mathbf{0}$ holds with not all $c_i$ zero. The empty list is independent by convention.

**Theorem 3.12 (Independence test).** The columns $\mathbf{a}_1, \dots, \mathbf{a}_n$ of $A$ are independent iff $N(A) = \{\mathbf{0}\}$ iff $A$ has $n$ pivots ($r = n$) iff $A$ has no free columns.

*Proof.* $\sum_j c_j\mathbf{a}_j = A\mathbf{c}$ by LA-1 (3.1). So a nontrivial relation among the columns is exactly a nonzero vector in $N(A)$: the columns are independent iff $N(A) = \{\mathbf{0}\}$. The rest is Corollary 3.10. ∎

In practice: put the vectors as the columns of a matrix, row reduce, count pivots. Independent iff every column is a pivot column. The project's `is_independent` does exactly this.

**Theorem 3.13.** (a) Any list containing $\mathbf{0}$ is dependent. (b) Any list with a repeated vector is dependent. (c) A one-vector list $\mathbf{v}$ is independent iff $\mathbf{v} \ne \mathbf{0}$.

*Proof.* (a) Say $\mathbf{v}_1 = \mathbf{0}$. Then $1\cdot\mathbf{v}_1 + 0\cdot\mathbf{v}_2 + \dots + 0\cdot\mathbf{v}_k = \mathbf{0}$ is a relation with a nonzero coefficient. (b) If $\mathbf{v}_i = \mathbf{v}_j$ with $i \ne j$, then $1\cdot\mathbf{v}_i + (-1)\cdot\mathbf{v}_j = \mathbf{0}$, other coefficients zero. (c) If $\mathbf{v} = \mathbf{0}$, (a) applies. If $\mathbf{v} \ne \mathbf{0}$ and $c\mathbf{v} = \mathbf{0}$ with $c \ne 0$, then $\mathbf{v} = 1\mathbf{v} = (c^{-1}c)\mathbf{v} = c^{-1}(c\mathbf{v}) = c^{-1}\mathbf{0} = \mathbf{0}$ by (M4), (M3), Theorem 3.2, a contradiction; so $c = 0$. ∎

**Theorem 3.14.** Any $n + 1$ vectors in $\mathbb{R}^n$ are dependent. More generally, any list of more than $m$ vectors in $\mathbb{R}^m$ is dependent.

*Proof.* Put the $k > m$ vectors as columns of an $m \times k$ matrix $A$. Each pivot sits in its own row, so $r \le m < k$: at least one column is free. By Theorem 3.12 the columns are dependent. ∎

**Theorem 3.15 (Independence is uniqueness of coefficients).** $\mathbf{v}_1, \dots, \mathbf{v}_k$ is independent iff every vector in its span has exactly one representation $\sum_i c_i\mathbf{v}_i$.

*Proof.* Example proof B in §3.0. ∎

A useful reformulation: a list is dependent iff some $\mathbf{v}_j$ is a combination of the *other* vectors in the list. (If $\sum_i c_i\mathbf{v}_i = \mathbf{0}$ with $c_j \ne 0$, then $\mathbf{v}_j = -c_j^{-1}\sum_{i\ne j}c_i\mathbf{v}_i$. Conversely if $\mathbf{v}_j = \sum_{i \ne j}d_i\mathbf{v}_i$, then $\sum_{i\ne j}d_i\mathbf{v}_i + (-1)\mathbf{v}_j = \mathbf{0}$ is a nontrivial relation.) Lemma 3.19 sharpens this.

### 3.5 Spanning sets, basis, coordinates

**Definition 3.16.** A list $\mathbf{v}_1, \dots, \mathbf{v}_k$ **spans** $V$ if $\operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_k\} = V$. A **basis** of $V$ is a list that is independent and spans $V$. $V$ is **finite-dimensional** if some finite list spans it, and **infinite-dimensional** otherwise.

**Theorem 3.17 (Coordinates).** A list $\mathcal{B} = (\mathbf{b}_1, \dots, \mathbf{b}_n)$ is a basis of $V$ iff every $\mathbf{v} \in V$ can be written as $\mathbf{v} = \sum_i c_i\mathbf{b}_i$ in exactly one way. The vector $[\mathbf{v}]_\mathcal{B} = (c_1, \dots, c_n)^\top \in \mathbb{R}^n$ is the **coordinate vector** of $\mathbf{v}$ in $\mathcal{B}$, and
$$[\mathbf{u} + \mathbf{v}]_\mathcal{B} = [\mathbf{u}]_\mathcal{B} + [\mathbf{v}]_\mathcal{B}, \qquad [c\mathbf{v}]_\mathcal{B} = c[\mathbf{v}]_\mathcal{B}.$$

*Proof.* "Spans" is existence of a representation for every $\mathbf{v}$; "independent" is uniqueness (Theorem 3.15). For the coordinate rules: if $\mathbf{u} = \sum_i c_i\mathbf{b}_i$ and $\mathbf{v} = \sum_i d_i\mathbf{b}_i$, then $\mathbf{u} + \mathbf{v} = \sum_i (c_i + d_i)\mathbf{b}_i$ is *a* representation, hence by uniqueness *the* representation, so $[\mathbf{u} + \mathbf{v}]_\mathcal{B} = (c_i + d_i)_i$. Likewise $c\mathbf{v} = \sum_i (cc_i)\mathbf{b}_i$. ∎

■ **What just happened.** This is the precise form of §2.3's "$\mathcal{P}_2$ is $\mathbb{R}^3$ in disguise." The map $\mathbf{v} \mapsto [\mathbf{v}]_\mathcal{B}$ is a one-to-one correspondence $V \leftrightarrow \mathbb{R}^n$ that converts vector addition into $\mathbb{R}^n$ addition and scaling into scaling. LA-4 calls such a map an **isomorphism** and proves every $n$-dimensional space is isomorphic to $\mathbb{R}^n$. Standard examples: $\mathbf{e}_1, \dots, \mathbf{e}_n$ for $\mathbb{R}^n$; $1, x, \dots, x^n$ for $\mathcal{P}_n$; the matrices $E_{ij}$ (a single $1$ in position $(i, j)$) for $\mathcal{M}_{m\times n}$.

**Theorem 3.18 (Special solutions are a basis of $N(A)$).** The special solutions $\mathbf{s}_1, \dots, \mathbf{s}_{n-r}$ form a basis of $N(A)$.

*Proof.* They span $N(A)$ by Theorem 3.8. Independence: suppose $\sum_j c_j\mathbf{s}_j = \mathbf{0}$. Look at the free entry $f_k$ of both sides: $\sum_j c_j(\mathbf{s}_j)_{f_k} = c_k$ on the left, $0$ on the right. So every $c_k = 0$. ∎

### 3.6 Dimension

The goal: prove that all bases of a finite-dimensional space have the same length. Everything rests on one lemma and one inequality.

**Lemma 3.19 (Linear dependence lemma).** Suppose $\mathbf{v}_1, \dots, \mathbf{v}_m$ is dependent. Then there is an index $j$ such that (a) $\mathbf{v}_j \in \operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_{j-1}\}$ (for $j = 1$ this means $\mathbf{v}_1 = \mathbf{0}$), and (b) removing $\mathbf{v}_j$ from the list does not change the span.

*Proof.* Take a nontrivial relation $\sum_i a_i\mathbf{v}_i = \mathbf{0}$ and let $j$ be the *largest* index with $a_j \ne 0$. Then $a_j\mathbf{v}_j = -\sum_{i<j}a_i\mathbf{v}_i$, so $\mathbf{v}_j = \sum_{i<j}(-a_i/a_j)\mathbf{v}_i$, which is (a); if $j = 1$ the sum is empty and $\mathbf{v}_1 = \mathbf{0}$. For (b): any $\mathbf{u} = \sum_i c_i\mathbf{v}_i$ in the original span becomes, after substituting the expression for $\mathbf{v}_j$, a combination of the list without $\mathbf{v}_j$. So the original span is contained in the smaller list's span; the reverse containment is automatic because fewer vectors can only span less. ∎

**Lemma 3.20 (Adjoining a vector outside the span).** If $\mathbf{v}_1, \dots, \mathbf{v}_k$ is independent and $\mathbf{w} \notin \operatorname{span}\{\mathbf{v}_1, \dots, \mathbf{v}_k\}$, then $\mathbf{v}_1, \dots, \mathbf{v}_k, \mathbf{w}$ is independent.

*Proof.* Suppose $\sum_i c_i\mathbf{v}_i + d\mathbf{w} = \mathbf{0}$. If $d \ne 0$ then $\mathbf{w} = -d^{-1}\sum_i c_i\mathbf{v}_i$ is in the span, contradicting the hypothesis. So $d = 0$, leaving $\sum_i c_i\mathbf{v}_i = \mathbf{0}$, and independence of the $\mathbf{v}_i$ gives all $c_i = 0$. ∎

▶ **Try it first.** Suppose $\mathbf{u}_1, \mathbf{u}_2$ is independent and $\mathbf{w}_1, \mathbf{w}_2, \mathbf{w}_3$ spans $V$. Adjoin $\mathbf{u}_1$ to the front of the $\mathbf{w}$ list. The list of four is dependent (why?). Use Lemma 3.19 to throw out one vector; argue it must be a $\mathbf{w}$. Repeat with $\mathbf{u}_2$. What have you shown about $2$ versus $3$?

**Theorem 3.21 (Steinitz exchange lemma).** If $\mathbf{u}_1, \dots, \mathbf{u}_k$ is independent in $V$ and $\mathbf{w}_1, \dots, \mathbf{w}_n$ spans $V$, then $k \le n$.

**Derivation.** We show, for each $j = 0, 1, \dots, k$, that there is a spanning list of length exactly $n$ consisting of $\mathbf{u}_1, \dots, \mathbf{u}_j$ followed by $n - j$ of the $\mathbf{w}$'s; in particular $n - j \ge 0$, and taking $j = k$ gives $k \le n$.

1. $j = 0$: the list $\mathbf{w}_1, \dots, \mathbf{w}_n$ spans by hypothesis.
2. Suppose the claim holds for $j - 1$, where $1 \le j \le k$: the list $L = (\mathbf{u}_1, \dots, \mathbf{u}_{j-1}, \mathbf{w}_{i_1}, \dots, \mathbf{w}_{i_{n-j+1}})$ spans $V$ and has length $n$.
3. Since $L$ spans $V$, $\mathbf{u}_j \in \operatorname{span} L$. Insert $\mathbf{u}_j$ right after $\mathbf{u}_{j-1}$ to get $L' = (\mathbf{u}_1, \dots, \mathbf{u}_j, \mathbf{w}_{i_1}, \dots, \mathbf{w}_{i_{n-j+1}})$, of length $n + 1$. $L'$ is dependent, because $\mathbf{u}_j$ is a combination of the other vectors in it (reformulation after Theorem 3.15).
4. Apply Lemma 3.19 to $L'$: some vector in $L'$ lies in the span of the vectors *before* it, and removing it leaves the span, which is $V$, unchanged.
5. That vector is not one of $\mathbf{u}_1, \dots, \mathbf{u}_j$: if $\mathbf{u}_i \in \operatorname{span}\{\mathbf{u}_1, \dots, \mathbf{u}_{i-1}\}$ for some $i \le j$, that gives a nontrivial relation among $\mathbf{u}_1, \dots, \mathbf{u}_i$ (write $\mathbf{u}_i$ as the combination and move it to one side), contradicting independence of the $\mathbf{u}$'s. So the removed vector is one of the $\mathbf{w}$'s. In particular there *was* a $\mathbf{w}$ in $L'$, i.e. $n - j + 1 \ge 1$.
6. After removal the list has length $n$, spans $V$, and consists of $\mathbf{u}_1, \dots, \mathbf{u}_j$ followed by $n - j$ of the $\mathbf{w}$'s. The claim holds for $j$. By induction it holds for $j = k$, and step 5 at the last stage gives $n - k \ge 0$. ∎

■ **What just happened.** Independent vectors can be swapped into a spanning list one at a time, each one evicting a $\mathbf{w}$, never another $\mathbf{u}$. You never run out of $\mathbf{w}$'s to evict, so there were at least as many $\mathbf{w}$'s as $\mathbf{u}$'s. *The length of an independent list never exceeds the length of a spanning list.* Theorem 3.14 was the special case $V = \mathbb{R}^m$ with the standard basis as the spanning list.

**Theorem 3.22 (Every spanning list contains a basis).** If $\mathbf{w}_1, \dots, \mathbf{w}_n$ spans $V$, then some sublist of it is a basis of $V$. In particular every finite-dimensional vector space has a basis.

*Proof.* If the list is independent, it is a basis. If not, Lemma 3.19(b) removes one vector without changing the span. Repeat. Each step shortens the list, so after at most $n$ steps the list is independent (the empty list is independent and spans $\{\mathbf{0}\}$) and still spans $V$. ∎

**Theorem 3.23 (All bases have the same length).** If $\mathcal{B}_1$ (length $k$) and $\mathcal{B}_2$ (length $n$) are bases of a finite-dimensional $V$, then $k = n$.

*Proof.* $\mathcal{B}_1$ is independent and $\mathcal{B}_2$ spans, so $k \le n$ by Theorem 3.21. Swap the roles: $n \le k$. ∎

**Definition 3.24 (Dimension).** $\dim V$ is the length of any basis of $V$. So $\dim\{\mathbf{0}\} = 0$ (empty basis), $\dim\mathbb{R}^n = n$, $\dim\mathcal{P}_n = n + 1$, $\dim\mathcal{M}_{m\times n} = mn$, and by Theorem 3.18, $\dim N(A) = n - r$.

The claim $\dim\mathcal{P}_n = n + 1$ needs $1, x, \dots, x^n$ to be a basis. It spans $\mathcal{P}_n$ by definition. It is independent: suppose $\sum_{j=0}^{n} a_jx^j = 0$ for every $x$. Fix $k$, differentiate $k$ times, and set $x = 0$. Terms with $j < k$ were killed by the differentiation; terms with $j > k$ still carry a factor $x^{j-k}$ and vanish at $0$; the term $j = k$ leaves $k!\,a_k$. So $k!\,a_k = 0$, hence $a_k = 0$, for every $k$. (Example 4.4 is the case $n = 2$.)

**Theorem 3.25 (Every independent list extends to a basis).** Let $V$ be finite-dimensional and $\mathbf{u}_1, \dots, \mathbf{u}_k$ independent. Then there are vectors $\mathbf{u}_{k+1}, \dots, \mathbf{u}_d$ such that $\mathbf{u}_1, \dots, \mathbf{u}_d$ is a basis of $V$.

*Proof.* Let $\mathbf{w}_1, \dots, \mathbf{w}_n$ span $V$. The list $(\mathbf{u}_1, \dots, \mathbf{u}_k, \mathbf{w}_1, \dots, \mathbf{w}_n)$ spans $V$. Run the removal process of Theorem 3.22 on it. At each step Lemma 3.19 removes a vector lying in the span of the vectors *before* it; as in step 5 of Theorem 3.21, that vector is never a $\mathbf{u}_i$, so only $\mathbf{w}$'s are removed. What remains is $(\mathbf{u}_1, \dots, \mathbf{u}_k, \text{some } \mathbf{w}\text{'s})$, independent and spanning. ∎

The constructive version, which the project implements as `extend_to_basis`: go through $\mathbf{w}_1, \dots, \mathbf{w}_n$ in order and keep $\mathbf{w}_i$ iff it is not in the span of what you have so far (Lemma 3.20 guarantees the list stays independent).

**Corollary 3.26.** Let $\dim V = n$.
(a) Any independent list of $n$ vectors in $V$ is a basis.
(b) Any spanning list of $n$ vectors is a basis.
(c) Any list of more than $n$ vectors is dependent; any list of fewer than $n$ vectors does not span $V$.
(d) If $W$ is a subspace of $V$, then $W$ is finite-dimensional, $\dim W \le n$, and $\dim W = n$ implies $W = V$.

*Proof.* (a) Extend to a basis (Theorem 3.25). The basis has length $n$ (Theorem 3.23), so nothing was added. (b) Shrink to a basis (Theorem 3.22); it has length $n$, so nothing was removed. (c) More than $n$ independent vectors would violate Theorem 3.21 against a basis of $V$ (a spanning list of length $n$). A spanning list of fewer than $n$ vectors would contain a basis of length less than $n$ (Theorem 3.22), violating Theorem 3.23. (d) Build a basis of $W$ greedily: start with the empty list; while its span is not $W$, pick $\mathbf{w} \in W$ outside the span and adjoin it. By Lemma 3.20 the list stays independent; by Theorem 3.21 (against a basis of $V$) its length never exceeds $n$; so the process stops, with an independent list spanning $W$, of length $\le n$. If $\dim W = n$, a basis of $W$ is an independent list of $n$ vectors of $V$, hence a basis of $V$ by (a), so $W = \operatorname{span} = V$. ∎

Now §2.1's guess is a theorem: a subspace of $\mathbb{R}^3$ has dimension $0$, $1$, $2$, or $3$, so it is $\{\mathbf{0}\}$, the span of one nonzero vector (a line), the span of two independent vectors (a plane), or all of $\mathbb{R}^3$.

**Infinite dimension.** $\mathbb{R}^\infty$ has no finite spanning list: the sequences $\mathbf{e}_1, \mathbf{e}_2, \dots, \mathbf{e}_{N+1}$ (a single $1$ in position $i$) are independent for every $N$, so if a list of length $N$ spanned, Theorem 3.21 would give $N + 1 \le N$. The same argument with $1, x, x^2, \dots$ (independent for every $n$, by the paragraph after Definition 3.24) shows the space of *all* polynomials is infinite-dimensional, and so is $C^\infty(\mathbb{R})$, which contains them. This is where Hilbert's integral equations and F2's Fourier series live; the finite-dimensional theory is the model for it.

### 3.7 The column space: pivot columns form a basis

**Lemma 3.27 (Same column relations).** Write the columns of $R$ as $\mathbf{r}_1, \dots, \mathbf{r}_n$. For every $\mathbf{c} \in \mathbb{R}^n$: $\sum_j c_j\mathbf{a}_j = \mathbf{0}$ iff $\sum_j c_j\mathbf{r}_j = \mathbf{0}$.

*Proof.* The left relation says $\mathbf{c} \in N(A)$, the right says $\mathbf{c} \in N(R)$, and these are the same set by Lemma 3.7. ∎

In words: a column of $A$ is a given combination of other columns of $A$ exactly when the same combination holds among the corresponding columns of $R$. The *relations* survive elimination even though the *columns themselves* change.

▶ **Try it first.** What do the pivot columns of $R$ look like? What do the free columns of $R$ look like, in terms of the pivot columns of $R$?

**Theorem 3.28 (Pivot columns are a basis of $C(A)$).** The pivot columns $\mathbf{a}_{p_1}, \dots, \mathbf{a}_{p_r}$ of $A$ (of $A$ itself, not of $R$) form a basis of $C(A)$. Hence $\dim C(A) = r$.

**Derivation.**

1. In $R$, pivot column $p_i$ is $\mathbf{e}_i \in \mathbb{R}^m$: a $1$ in row $i$, zeros elsewhere (definition of RREF). So $\mathbf{r}_{p_1}, \dots, \mathbf{r}_{p_r} = \mathbf{e}_1, \dots, \mathbf{e}_r$, which are independent (a relation $\sum_i c_i\mathbf{e}_i = \mathbf{0}$ has $c_i$ as its $i$-th entry).
2. By Lemma 3.27, any relation among $\mathbf{a}_{p_1}, \dots, \mathbf{a}_{p_r}$ would be a relation among $\mathbf{e}_1, \dots, \mathbf{e}_r$. There is none, so the pivot columns of $A$ are independent.
3. Free column $f_j$ of $R$ has nonzero entries only in rows $1, \dots, r$ (rows below $r$ are zero), so $\mathbf{r}_{f_j} = \sum_{i=1}^r R_{i,f_j}\mathbf{e}_i = \sum_i R_{i,f_j}\mathbf{r}_{p_i}$: a combination of the pivot columns of $R$.
4. By Lemma 3.27, the same relation holds in $A$: $\mathbf{a}_{f_j} = \sum_i R_{i,f_j}\mathbf{a}_{p_i}$. So every column of $A$ is in $\operatorname{span}\{\mathbf{a}_{p_1}, \dots, \mathbf{a}_{p_r}\}$, hence (Theorem 3.5, smallest subspace) $C(A) \subseteq \operatorname{span}\{\mathbf{a}_{p_i}\} \subseteq C(A)$.
5. Independent and spanning: a basis. Its length is $r$. ∎

■ **What just happened.** $R$ tells you *which* columns of $A$ to keep, but you keep the columns of $A$. In general $C(R) \ne C(A)$: for $A = \begin{pmatrix}1 \\ 1\end{pmatrix}$, $R = \begin{pmatrix}1 \\ 0\end{pmatrix}$ and the column spaces are different lines. Row operations change the column space; they preserve the column *relations*. The rank of LA-2 (number of pivots) now has a meaning independent of the algorithm: it is $\dim C(A)$.

**Corollary (The RREF is unique).** If $EA = R$ and $E'A = R'$ with $E, E'$ invertible and $R, R'$ both in RREF, then $R = R'$. (LA-2 §3.9 proves this from its exercise C3; here is the same argument in this module's language.)

*Proof.* *Pivot positions.* In any RREF, column $j$ is a pivot column iff it is not a combination of the columns to its left. If $j = p_i$, column $j$ is $\mathbf{e}_i$, while every column left of $p_i$ is zero in rows $i, \dots, m$ (row $i$ is zero left of its pivot, and so are the rows below it, whose pivots are further right); so $\mathbf{e}_i$ is not a combination of them. If $j = f$ is free, step 3 of Theorem 3.28 writes it as a combination of pivot columns, and only those to its left appear (if $p_i > f$ then $R_{i,f} = 0$, since row $i$ is zero left of $p_i$). By Lemma 3.27, "$\mathbf{r}_j = \sum_{k<j}d_k\mathbf{r}_k$" holds iff "$\mathbf{a}_j = \sum_{k<j}d_k\mathbf{a}_k$" does (move everything to one side to get a relation). So column $j$ is a pivot column of $R$ iff $\mathbf{a}_j \notin \operatorname{span}\{\mathbf{a}_1, \dots, \mathbf{a}_{j-1}\}$. That condition mentions only $A$, so $R$ and $R'$ have the same pivot columns $p_1 < \dots < p_r$, hence the same pivot columns $\mathbf{e}_1, \dots, \mathbf{e}_r$ and the same zero rows $r + 1, \dots, m$. *Free entries.* For a free column $f$, step 4 of Theorem 3.28 gives $\mathbf{a}_f = \sum_i R_{i,f}\,\mathbf{a}_{p_i}$ and equally $\mathbf{a}_f = \sum_i R'_{i,f}\,\mathbf{a}_{p_i}$. The $\mathbf{a}_{p_i}$ are independent (step 2), so coefficients are unique (Theorem 3.15): $R_{i,f} = R'_{i,f}$. Every entry agrees. ∎

So "the" RREF in the conventions box is earned: the free entries of $R$ are the coordinates of the free columns of $A$ in the basis of pivot columns. The project checks it by reducing each random matrix with two different pivoting strategies and comparing the results.

### 3.8 The row space: row rank equals column rank

▶ **Try it first.** Show that each row of $EA$ is a combination of rows of $A$. Then: why are the nonzero rows of $R$ independent? (Look at the entries in the pivot columns.)

**Theorem 3.29 (Row space).** (a) If $E$ is invertible, $C((EA)^\top) = C(A^\top)$: row operations do not change the row space. (b) The nonzero rows of $R$ form a basis of $C(A^\top)$. (c) $\dim C(A^\top) = r = \dim C(A)$.

**Derivation.**

1. (a) Row $i$ of $EA$ is $\sum_k E_{ik}(\text{row } k \text{ of } A)$, by the row form of matrix multiplication (LA-2 Theorem 3.1; entrywise, $(EA)_{ij} = \sum_k E_{ik}A_{kj}$). So every row of $EA$ is in $C(A^\top)$, and $C((EA)^\top) \subseteq C(A^\top)$ by Theorem 3.5. Apply the same argument to $A = E^{-1}(EA)$: $C(A^\top) \subseteq C((EA)^\top)$.
2. (b) By (a) with $EA = R$, the rows of $R$ span $C(A^\top)$; the zero rows contribute nothing, so the $r$ nonzero rows span. Independence: suppose $\sum_{i=1}^r c_i(\text{row } i \text{ of } R) = \mathbf{0}^\top$. Look at entry $p_i$: row $i$ has a $1$ there and every other row has a $0$ there (pivot column of RREF). So $c_i = 0$ for each $i$.
3. (c) The basis in (b) has $r$ elements; Theorem 3.28 gave $r$ for the column space. ∎

■ **What just happened.** "Row rank $=$ column rank" is the statement that elimination, which works on rows, counts something about columns. The proof shows both counts are the number of pivots. Exercise B6 gives a second proof via the factorization $A = CR$.

### 3.9 The four fundamental subspaces

**Theorem 3.30 (Fundamental theorem of linear algebra, part 1).** Let $A$ be $m \times n$ with rank $r$. Then
$$\dim C(A) = r, \qquad \dim N(A) = n - r, \qquad \dim C(A^\top) = r, \qquad \dim N(A^\top) = m - r.$$
In particular (**Rank–Nullity**):
$$\dim C(A) + \dim N(A) = n. \tag{3.5}$$

*Proof.* The first three are Theorems 3.28, 3.18, 3.29. For the fourth, $A^\top$ is $n \times m$ with rank $r$ (Theorem 3.29(c)), so Theorem 3.18 applied to $A^\top$ gives $\dim N(A^\top) = m - r$. ∎

**Theorem 3.31 (A basis of $N(A^\top)$ from elimination).** Let $EA = R$ with $E$ invertible, and let $\mathbf{y}_1^\top, \dots, \mathbf{y}_m^\top$ be the rows of $E$. Then $\mathbf{y}_{r+1}, \dots, \mathbf{y}_m$ form a basis of $N(A^\top)$.

*Proof.* Row $i$ of $EA$ is $\mathbf{y}_i^\top A$, and it equals row $i$ of $R$, which is zero for $i > r$. So $\mathbf{y}_i^\top A = \mathbf{0}^\top$, i.e. $A^\top\mathbf{y}_i = \mathbf{0}$, for $i = r + 1, \dots, m$: these $m - r$ vectors lie in $N(A^\top)$. They are independent: the rows of $E$ are the columns of $E^\top$, which is invertible (LA-2), so $N(E^\top) = \{\mathbf{0}\}$ and Theorem 3.12 says its columns are independent; a sublist of an independent list is independent. An independent list of $m - r$ vectors in the $(m - r)$-dimensional space $N(A^\top)$ is a basis by Corollary 3.26(a). ∎

To get $E$ in practice, reduce the augmented matrix $[A \mid I_m]$; it becomes $[R \mid E]$. The last $m - r$ rows of the right block are your basis. (Alternatively: compute $N(A^\top)$ by special solutions of $A^\top$. Both are in the project, and it checks that they span the same space.)

**The big picture.** Strang's diagram, which you should be able to draw from memory with the dimensions filled in. Example 4.2 fills it in with actual vectors.

```
          R^n  (inputs x)                               R^m  (outputs b)
   +-----------------------------+                +-----------------------------+
   |                             |                |                             |
   |     row space  C(A^T)       |       A        |     column space  C(A)      |
   |     dimension  r            | -------------> |     dimension  r            |
   |                             |   x_row -> b   |                             |
   + - - - - - -  _|_  - - - - - +                + - - - - - -  _|_  - - - - - +
   |                             |                |                             |
   |     null space  N(A)        |   A x_null = 0 |   left null space  N(A^T)   |
   |     dimension  n - r        | -------------> |     dimension  m - r        |
   |                             |       0        |                             |
   +-----------------------------+                +-----------------------------+
          r + (n - r) = n                                r + (m - r) = m
```

The symbol $\perp$ between the boxes is the content of LA-5, but one line of it is available now.

**Preview: $N(A) \perp C(A^\top)$.** If $A\mathbf{x} = \mathbf{0}$, then row by row (LA-1 (3.2)) each row of $A$ has dot product $0$ with $\mathbf{x}$; so $\mathbf{x}$ is orthogonal to every row, hence (by linearity of the dot product) to every combination of rows. That is,
$$\mathbf{x} \in N(A),\ \mathbf{w} \in C(A^\top) \quad\Longrightarrow\quad \mathbf{x}\cdot\mathbf{w} = 0. \tag{3.6}$$
Applying this to $A^\top$: $N(A^\top) \perp C(A)$. LA-5 proves more: the dimensions $r + (n - r) = n$ force these pairs to be orthogonal *complements*, and every $\mathbf{x} \in \mathbb{R}^n$ splits uniquely as $\mathbf{x}_{\text{row}} + \mathbf{x}_{\text{null}}$. The map $A$ sends $\mathbf{x}_{\text{null}}$ to $\mathbf{0}$ and $\mathbf{x}_{\text{row}}$ to $A\mathbf{x}$: all the action is between the two top boxes.

All of this is about *real* matrices (conventions box). Over $\mathbb{C}$, (3.6) is still true as written, since $\mathbf{x}^\top A^\top\mathbf{y} = (A\mathbf{x})^\top\mathbf{y} = 0$, but $\mathbf{x}^\top\mathbf{w} = 0$ is not orthogonality there, and the complement statement fails: for $A = (1, i)$, $N(A) = \operatorname{span}\{(1, i)^\top\} = C(A^\top)$, so the two subspaces coincide instead of splitting $\mathbb{C}^2$. LA-8 repairs it with $A^\dagger$ in place of $A^\top$.

### 3.10 Rank of a product

▶ **Try it first.** Each column of $AB$ is $A$ times a column of $B$ (LA-1 (3.3)). Which subspace must it lie in? And if $B\mathbf{x} = \mathbf{0}$, what is $AB\mathbf{x}$?

**Theorem 3.32.** For $A$ ($m \times p$) and $B$ ($p \times n$),
$$\operatorname{rank}(AB) \le \min\big(\operatorname{rank} A, \operatorname{rank} B\big). \tag{3.7}$$

**Derivation.**

1. Column $k$ of $AB$ is $A\mathbf{b}_k \in C(A)$. So $C(AB) \subseteq C(A)$, and by Corollary 3.26(d), $\operatorname{rank}(AB) = \dim C(AB) \le \dim C(A) = \operatorname{rank} A$.
2. If $B\mathbf{x} = \mathbf{0}$ then $AB\mathbf{x} = A\mathbf{0} = \mathbf{0}$. So $N(B) \subseteq N(AB)$, and $\dim N(B) \le \dim N(AB)$.
3. $B$ and $AB$ both have $n$ columns. Rank–Nullity (3.5) for each: $n - \operatorname{rank} B \le n - \operatorname{rank}(AB)$, i.e. $\operatorname{rank}(AB) \le \operatorname{rank} B$. ∎

■ **What just happened.** Multiplying can only lose rank, never create it. Strict inequality happens: $A = \begin{pmatrix}1 & 0\\0 & 0\end{pmatrix}$, $B = \begin{pmatrix}0 & 0\\0 & 1\end{pmatrix}$ both have rank $1$ and $AB = 0$ has rank $0$. If $A$ is invertible, $\operatorname{rank}(AB) = \operatorname{rank} B$ (exercise C3).

### 3.11 Sums, intersections, direct sums

Two subspaces $U, W$ of $V$ give two new ones. $U \cap W$ is a subspace (Example proof A in §3.0). So is the **sum** $U + W = \{\mathbf{u} + \mathbf{w} : \mathbf{u} \in U,\ \mathbf{w} \in W\}$: it contains $\mathbf{0} = \mathbf{0} + \mathbf{0}$, and $(\mathbf{u} + \mathbf{w}) + (\mathbf{u}' + \mathbf{w}') = (\mathbf{u} + \mathbf{u}') + (\mathbf{w} + \mathbf{w}')$ and $c(\mathbf{u} + \mathbf{w}) = c\mathbf{u} + c\mathbf{w}$ stay in it (exercise B1(a) writes this out with the axioms cited and shows $U + W$ is the smallest subspace containing $U \cup W$). The same goes for $U_1 + \dots + U_k$.

**Definition 3.33 (Direct sum).** Let $U_1, \dots, U_k$ be subspaces of $V$. $V$ is the **direct sum** $V = U_1 \oplus \dots \oplus U_k$ if every $\mathbf{v} \in V$ can be written as $\mathbf{v} = \mathbf{u}_1 + \dots + \mathbf{u}_k$ with each $\mathbf{u}_i \in U_i$, in exactly one way.

Existence of the splitting says $U_1 + \dots + U_k = V$. Uniqueness is the new part, and, as in Theorem 3.15, it reduces to a statement about $\mathbf{0}$.

**Theorem 3.34 (Testing a direct sum).** (a) $V = U_1 \oplus \dots \oplus U_k$ iff $U_1 + \dots + U_k = V$ and the only way to write $\mathbf{0} = \mathbf{u}_1 + \dots + \mathbf{u}_k$ with $\mathbf{u}_i \in U_i$ is with every $\mathbf{u}_i = \mathbf{0}$. (b) For two subspaces: $V = U \oplus W$ iff $U + W = V$ and $U \cap W = \{\mathbf{0}\}$. (c) If $V = U_1 \oplus \dots \oplus U_k$, then a basis of $U_1$, followed by a basis of $U_2$, …, followed by a basis of $U_k$, is a basis of $V$. So $\dim V = \dim U_1 + \dots + \dim U_k$.

*Proof.* (a) $(\Rightarrow)$ $\mathbf{0} = \mathbf{0} + \dots + \mathbf{0}$ is one splitting of $\mathbf{0}$, so by uniqueness it is the only one. $(\Leftarrow)$ If $\mathbf{v} = \sum_i\mathbf{u}_i = \sum_i\mathbf{u}_i'$, subtract: $\mathbf{0} = \sum_i(\mathbf{u}_i - \mathbf{u}_i')$ with $\mathbf{u}_i - \mathbf{u}_i' \in U_i$, so every $\mathbf{u}_i - \mathbf{u}_i' = \mathbf{0}$. (b) By (a) it is enough to show that "$\mathbf{0} = \mathbf{u} + \mathbf{w}$ forces $\mathbf{u} = \mathbf{w} = \mathbf{0}$" holds iff $U \cap W = \{\mathbf{0}\}$. If the first holds and $\mathbf{x} \in U \cap W$, then $\mathbf{0} = \mathbf{x} + (-\mathbf{x})$ with $\mathbf{x} \in U$, $-\mathbf{x} \in W$, so $\mathbf{x} = \mathbf{0}$. If $U \cap W = \{\mathbf{0}\}$ and $\mathbf{0} = \mathbf{u} + \mathbf{w}$, then $\mathbf{u} = -\mathbf{w}$ lies in both, so $\mathbf{u} = \mathbf{0}$ and then $\mathbf{w} = \mathbf{0}$. (c) *Spans:* every $\mathbf{v}$ is $\sum_i\mathbf{u}_i$, and each $\mathbf{u}_i$ is a combination of the basis of $U_i$. *Independent:* given a relation among all the basis vectors, let $\mathbf{u}_i$ be the part of it that uses the basis of $U_i$; then $\mathbf{u}_i \in U_i$ and $\mathbf{0} = \mathbf{u}_1 + \dots + \mathbf{u}_k$. By (a) every $\mathbf{u}_i = \mathbf{0}$, and independence of the basis of $U_i$ makes each of its coefficients zero. ∎

For $k \ge 3$, pairwise intersections $\{\mathbf{0}\}$ are *not* enough. The lines $U_1 = \operatorname{span}\{(1, 0)\}$, $U_2 = \operatorname{span}\{(0, 1)\}$, $U_3 = \operatorname{span}\{(1, 1)\}$ in $\mathbb{R}^2$ meet pairwise only in $\mathbf{0}$, yet $(1, 0) + (0, 1) + (-1, -1) = \mathbf{0}$ is a nonzero splitting of $\mathbf{0}$.

The worked instance is exercise B3(d): $\mathcal{M}_{n\times n} = \mathrm{Sym}_n \oplus \mathrm{Skew}_n$. Every square matrix is $\tfrac12(A + A^\top) + \tfrac12(A - A^\top)$, and a matrix that is both symmetric and skew is $0$, so Theorem 3.34(b) applies, and (c) says the two dimensions add to $n^2$.

▶ **Try it first.** Let $U \ne W$ be two planes through $\mathbf{0}$ in $\mathbb{R}^3$. What are $U + W$ and $U \cap W$? Guess a formula linking $\dim(U + W)$, $\dim U$, $\dim W$, $\dim(U \cap W)$. Why does the intersection have to be subtracted?

**Theorem 3.35 (Dimension formula).** If $U$ and $W$ are subspaces of a finite-dimensional $V$, then
$$\dim(U + W) = \dim U + \dim W - \dim(U \cap W). \tag{3.8}$$

**Derivation.**

1. $U \cap W$ is a subspace of $V$, hence finite-dimensional (Corollary 3.26(d)). Let $\mathbf{z}_1, \dots, \mathbf{z}_k$ be a basis, $k = \dim(U \cap W)$.
2. This list is independent and lies in $U$. Extend it to a basis $\mathbf{z}_1, \dots, \mathbf{z}_k, \mathbf{u}_1, \dots, \mathbf{u}_p$ of $U$ (Theorem 3.25), so $\dim U = k + p$. Separately, extend it to a basis $\mathbf{z}_1, \dots, \mathbf{z}_k, \mathbf{w}_1, \dots, \mathbf{w}_q$ of $W$, so $\dim W = k + q$.
3. *Spans.* Let $L = (\mathbf{z}_1, \dots, \mathbf{z}_k, \mathbf{u}_1, \dots, \mathbf{u}_p, \mathbf{w}_1, \dots, \mathbf{w}_q)$. Every vector of $L$ lies in $U$ or in $W$, hence in $U + W$, so $\operatorname{span} L \subseteq U + W$ (Theorem 3.5). Conversely, $\mathbf{u} + \mathbf{w} \in U + W$ is a combination of the $\mathbf{z}$'s and $\mathbf{u}$'s plus a combination of the $\mathbf{z}$'s and $\mathbf{w}$'s, hence a combination of $L$.
4. *Independent, stage 1.* Suppose $\sum_i a_i\mathbf{z}_i + \sum_j b_j\mathbf{u}_j + \sum_l c_l\mathbf{w}_l = \mathbf{0}$ $(*)$. Then $\mathbf{v} = \sum_l c_l\mathbf{w}_l$ equals $-\sum_i a_i\mathbf{z}_i - \sum_j b_j\mathbf{u}_j$. The first expression is in $W$, the second in $U$, so $\mathbf{v} \in U \cap W$ and $\mathbf{v} = \sum_i d_i\mathbf{z}_i$ for some $d_i$. Then $\sum_i d_i\mathbf{z}_i - \sum_l c_l\mathbf{w}_l = \mathbf{0}$ is a relation among the basis of $W$ from step 2, so every $c_l = 0$.
5. *Stage 2.* With the $c_l$ gone, $(*)$ is a relation among the basis of $U$ from step 2, so every $a_i$ and every $b_j$ is $0$.
6. $L$ is a basis of $U + W$ of length $k + p + q = (k + p) + (k + q) - k$. ∎

■ **What just happened.** Adding $\dim U + \dim W$ counts the intersection twice; (3.8) subtracts it once. The real work is step 4: a relation among the combined list manufactures a vector lying in both $U$ and $W$, and the basis of $U \cap W$, built first, absorbs it. Two consequences get used later. (i) For two subspaces, $V = U \oplus W$ iff $U \cap W = \{\mathbf{0}\}$ and $\dim U + \dim W = \dim V$ (exercise C1(a)). (ii) In $\mathbb{R}^n$, if $\dim U + \dim W > n$ then $U \cap W$ contains a nonzero vector, because $\dim(U \cap W) = \dim U + \dim W - \dim(U + W) \ge \dim U + \dim W - n \ge 1$. For two distinct planes through $\mathbf{0}$ in $\mathbb{R}^3$: $U + W$ is strictly bigger than $U$ (otherwise $W \subseteq U$, and equal dimensions force $W = U$ by Corollary 3.26(d)), so $U + W = \mathbb{R}^3$ and $\dim(U \cap W) = 2 + 2 - 3 = 1$: they meet in a line. LA-8's Courant–Fischer theorem and the spectral-norm Eckart–Young theorem (LA-9 Theorem 3.8) both run on (ii). Exercise B8 checks (3.8) numerically.

## 4. Worked examples

### Example 4.1: Testing independence in $\mathbb{R}^4$

Are $\mathbf{v}_1 = (1, 2, 0, 1)^\top$, $\mathbf{v}_2 = (0, 1, 1, 1)^\top$, $\mathbf{v}_3 = (2, 5, 1, 3)^\top$ independent?

Theorem 3.12: make them the columns of $A$ and reduce.
$$
A = \begin{pmatrix} 1 & 0 & 2 \\ 2 & 1 & 5 \\ 0 & 1 & 1 \\ 1 & 1 & 3 \end{pmatrix}
\xrightarrow{\ R_2 - 2R_1,\ R_4 - R_1\ }
\begin{pmatrix} 1 & 0 & 2 \\ 0 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 1 & 1 \end{pmatrix}
\xrightarrow{\ R_3 - R_2,\ R_4 - R_2\ }
R = \begin{pmatrix} 1 & 0 & 2 \\ 0 & 1 & 1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.
$$
Two pivots (columns 1, 2), one free column (3). Dependent. The special solution: $x_3 = 1$, $x_1 = -2$, $x_2 = -1$, so $\mathbf{s} = (-2, -1, 1)^\top$ and the relation is $-2\mathbf{v}_1 - \mathbf{v}_2 + \mathbf{v}_3 = \mathbf{0}$, i.e. $\mathbf{v}_3 = 2\mathbf{v}_1 + \mathbf{v}_2$. Check: $2(1, 2, 0, 1) + (0, 1, 1, 1) = (2, 5, 1, 3)$. ✓ Notice that $R$'s free column $(2, 1, 0, 0)^\top$ is literally the recipe (Theorem 3.28 step 4).

Now replace $\mathbf{v}_3$ by $\mathbf{w} = (2, 5, 1, 4)^\top$. The same two steps give rows $(1, 0, 2), (0, 1, 1), (0, 1, 1), (0, 1, 2)$; then $R_3 - R_2 = \mathbf{0}$ and $R_4 - R_2 = (0, 0, 1)$. Three pivots: independent. Three independent vectors in $\mathbb{R}^4$ do not span it (Corollary 3.26(c)); the project's `extend_to_basis` finds a fourth.

### Example 4.2: All four subspaces of a $3 \times 4$ matrix

$$
A = \begin{pmatrix} 1 & 2 & 0 & 1 \\ 2 & 4 & 1 & 4 \\ 1 & 2 & 1 & 3 \end{pmatrix}, \qquad m = 3,\ n = 4.
$$

**Eliminate, recording $E$.** Reduce $[A \mid I_3]$:
$$
\left(\begin{array}{cccc|ccc} 1 & 2 & 0 & 1 & 1 & 0 & 0 \\ 2 & 4 & 1 & 4 & 0 & 1 & 0 \\ 1 & 2 & 1 & 3 & 0 & 0 & 1 \end{array}\right)
\xrightarrow{R_2 - 2R_1,\ R_3 - R_1}
\left(\begin{array}{cccc|ccc} 1 & 2 & 0 & 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & 2 & -2 & 1 & 0 \\ 0 & 0 & 1 & 2 & -1 & 0 & 1 \end{array}\right)
\xrightarrow{R_3 - R_2}
\left(\begin{array}{cccc|ccc} 1 & 2 & 0 & 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & 2 & -2 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 & -1 & 1 \end{array}\right).
$$
The left block is already reduced (the pivot in column 3 has a $0$ above it). So
$$
R = \begin{pmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 \end{pmatrix}, \qquad
E = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 1 & -1 & 1 \end{pmatrix}, \qquad r = 2.
$$
Pivot columns $p_1 = 1$, $p_2 = 3$. Free columns $f_1 = 2$, $f_2 = 4$. Check $EA = R$ on one entry: row 3 of $E$ times column 4 of $A$ is $1\cdot 1 - 1\cdot 4 + 1\cdot 3 = 0$. ✓

**Column space** (Theorem 3.28): columns 1 and 3 *of $A$*:
$$C(A) = \operatorname{span}\left\{\begin{pmatrix}1\\2\\1\end{pmatrix}, \begin{pmatrix}0\\1\\1\end{pmatrix}\right\} \subseteq \mathbb{R}^3, \qquad \dim = 2.$$
A plane in $\mathbb{R}^3$. The free columns are combinations with coefficients read off $R$: $\mathbf{a}_2 = 2\mathbf{a}_1$ and $\mathbf{a}_4 = 1\mathbf{a}_1 + 2\mathbf{a}_3 = (1, 2, 1) + (0, 2, 2) = (1, 4, 3)$. ✓

**Null space** (Theorem 3.18): two special solutions. Row 1 of $R\mathbf{x} = \mathbf{0}$: $x_1 = -2x_2 - x_4$. Row 2: $x_3 = -2x_4$.
$x_2 = 1, x_4 = 0$: $\mathbf{s}_1 = (-2, 1, 0, 0)^\top$. $x_2 = 0, x_4 = 1$: $\mathbf{s}_2 = (-1, 0, -2, 1)^\top$.
$$N(A) = \operatorname{span}\{(-2, 1, 0, 0)^\top, (-1, 0, -2, 1)^\top\} \subseteq \mathbb{R}^4, \qquad \dim = 2 = n - r.$$
Check $A\mathbf{s}_2$: row 1: $-1 + 0 + 0 + 1 = 0$; row 2: $-2 + 0 - 2 + 4 = 0$; row 3: $-1 + 0 - 2 + 3 = 0$. ✓

**Row space** (Theorem 3.29): nonzero rows of $R$:
$$C(A^\top) = \operatorname{span}\{(1, 2, 0, 1)^\top, (0, 0, 1, 2)^\top\} \subseteq \mathbb{R}^4, \qquad \dim = 2.$$
Check (3.6): $(1, 2, 0, 1)\cdot(-1, 0, -2, 1) = -1 + 0 + 0 + 1 = 0$; $(0, 0, 1, 2)\cdot(-1, 0, -2, 1) = -2 + 2 = 0$; both rows against $\mathbf{s}_1$: $-2 + 2 = 0$ and $0$. ✓ All four dot products vanish.

**Left null space** (Theorem 3.31): the last $m - r = 1$ row of $E$:
$$N(A^\top) = \operatorname{span}\{(1, -1, 1)^\top\} \subseteq \mathbb{R}^3, \qquad \dim = 1 = m - r.$$
Check $\mathbf{y}^\top A = (1 - 2 + 1,\ 2 - 4 + 2,\ 0 - 1 + 1,\ 1 - 4 + 3) = (0, 0, 0, 0)$. ✓ And $\mathbf{y} \perp C(A)$: $(1, -1, 1)\cdot(1, 2, 1) = 0$, $(1, -1, 1)\cdot(0, 1, 1) = 0$. ✓ Meaning: row 1 $-$ row 2 $+$ row 3 of $A$ is the zero row, which is exactly the dependency elimination found.

**The big picture, filled in:**

```
          R^4                                              R^3
   +---------------------------+                  +---------------------------+
   |  C(A^T)   dim 2           |        A         |  C(A)     dim 2           |
   |  (1,2,0,1), (0,0,1,2)     | ---------------> |  (1,2,1), (0,1,1)         |
   + - - - - - - _|_ - - - - - +                  + - - - - - - _|_ - - - - - +
   |  N(A)     dim 2           |     A x = 0      |  N(A^T)   dim 1           |
   |  (-2,1,0,0), (-1,0,-2,1)  | ---------------> |  (1,-1,1)                 |
   +---------------------------+                  +---------------------------+
             2 + 2 = 4                                      2 + 1 = 3
```

**Complete solution.** Take $\mathbf{b} = (1, 3, 2)^\top$. $E\mathbf{b} = (1,\ 3 - 2,\ 1 - 3 + 2)^\top = (1, 1, 0)^\top$; last entry $0$, consistent. Free variables $0$: $x_1 = 1$, $x_3 = 1$, so $\mathbf{x}_p = (1, 0, 1, 0)^\top$ (indeed $\mathbf{b} = \mathbf{a}_1 + \mathbf{a}_3$). By (3.3),
$$\mathbf{x} = \begin{pmatrix}1\\0\\1\\0\end{pmatrix} + c_1\begin{pmatrix}-2\\1\\0\\0\end{pmatrix} + c_2\begin{pmatrix}-1\\0\\-2\\1\end{pmatrix}, \qquad c_1, c_2 \in \mathbb{R}.$$
Now $\mathbf{b}' = (1, 3, 3)^\top$: $E\mathbf{b}' = (1, 1, 1)^\top$, last entry $1 \ne 0$. No solution. Equivalently $\mathbf{y}\cdot\mathbf{b}' = 1 - 3 + 3 = 1 \ne 0$, so $\mathbf{b}' \notin C(A)$: the left null space is the test for solvability.

### Example 4.3: A rank-one matrix

$$A = \mathbf{u}\mathbf{v}^\top = \begin{pmatrix}1\\2\\3\end{pmatrix}\begin{pmatrix}1 & -1 & 2\end{pmatrix} = \begin{pmatrix} 1 & -1 & 2 \\ 2 & -2 & 4 \\ 3 & -3 & 6 \end{pmatrix}.$$

Every column is a multiple of $\mathbf{u}$ (column $j$ is $v_j\mathbf{u}$), every row a multiple of $\mathbf{v}^\top$. Elimination: $R_2 - 2R_1$ and $R_3 - 3R_1$ kill rows 2 and 3:
$$R = \begin{pmatrix} 1 & -1 & 2 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}, \qquad E = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ -3 & 0 & 1 \end{pmatrix}, \qquad r = 1.$$

- $C(A) = \operatorname{span}\{(1, 2, 3)^\top\} = \operatorname{span}\{\mathbf{u}\}$. A line. $\dim 1$.
- $C(A^\top) = \operatorname{span}\{(1, -1, 2)^\top\} = \operatorname{span}\{\mathbf{v}\}$. A line. $\dim 1$.
- $N(A)$: free $x_2, x_3$; $x_1 = x_2 - 2x_3$. $\mathbf{s}_1 = (1, 1, 0)^\top$, $\mathbf{s}_2 = (-2, 0, 1)^\top$. $\dim 2 = 3 - 1$. This is the plane $\mathbf{v}\cdot\mathbf{x} = 0$ (check: $1 - 1 = 0$, $-2 + 2 = 0$), as (3.6) predicts: $A\mathbf{x} = \mathbf{u}(\mathbf{v}^\top\mathbf{x}) = \mathbf{0}$ iff $\mathbf{v}\cdot\mathbf{x} = 0$.
- $N(A^\top)$: rows 2, 3 of $E$: $(-2, 1, 0)^\top, (-3, 0, 1)^\top$. $\dim 2 = 3 - 1$. This is the plane $\mathbf{u}\cdot\mathbf{y} = 0$ (check: $-2 + 2 = 0$, $-3 + 3 = 0$).

Dimensions $1, 2, 1, 2$: $1 + 2 = 3$ on each side. A rank-one matrix is a line-to-line map that annihilates a plane. LA-9 proves every rank-$r$ matrix is a sum of $r$ such pieces $\sigma_i\mathbf{u}_i\mathbf{v}_i^\top$ with the $\mathbf{u}_i$ an orthonormal basis of $C(A)$ and the $\mathbf{v}_i$ an orthonormal basis of $C(A^\top)$.

### Example 4.4: Two bases of $\mathcal{P}_2$ and converting coordinates

**Claim.** $\mathcal{B} = (1, x, x^2)$ and $\mathcal{C} = (1,\ 1 + x,\ 1 + x + x^2)$ are both bases of $\mathcal{P}_2$.

*$\mathcal{B}$ spans* by definition of $\mathcal{P}_2$. *$\mathcal{B}$ is independent:* suppose $a + bx + cx^2 = 0$ as a function. At $x = 0$: $a = 0$. Differentiate: $b + 2cx = 0$; at $x = 0$: $b = 0$. Differentiate again: $2c = 0$. So $\dim\mathcal{P}_2 = 3$.

*$\mathcal{C}$:* write each vector in $\mathcal{B}$-coordinates and make them the columns of
$$M = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{pmatrix}.$$
A relation $c_1\cdot 1 + c_2(1 + x) + c_3(1 + x + x^2) = 0$ is, in coordinates (Theorem 3.17), $M\mathbf{c} = \mathbf{0}$. $M$ is upper triangular with nonzero diagonal, so it has three pivots and $N(M) = \{\mathbf{0}\}$ (Theorem 3.12): $\mathcal{C}$ is independent. Three independent vectors in a three-dimensional space form a basis (Corollary 3.26(a)). ✓

**Convert.** $p(x) = 2 + 3x - x^2$ has $[p]_\mathcal{B} = (2, 3, -1)^\top$. For $[p]_\mathcal{C} = (c_1, c_2, c_3)^\top$ we need $c_1 + c_2(1 + x) + c_3(1 + x + x^2) = p$, i.e. $M[p]_\mathcal{C} = [p]_\mathcal{B}$. Back-substitute: $c_3 = -1$; $c_2 + c_3 = 3 \Rightarrow c_2 = 4$; $c_1 + c_2 + c_3 = 2 \Rightarrow c_1 = -1$. So
$$p = -1\cdot 1 + 4(1 + x) - 1(1 + x + x^2), \qquad [p]_\mathcal{C} = (-1, 4, -1)^\top.$$
Check: $-1 + 4 + 4x - 1 - x - x^2 = 2 + 3x - x^2$. ✓ In general $[p]_\mathcal{B} = M[p]_\mathcal{C}$ and $[p]_\mathcal{C} = M^{-1}[p]_\mathcal{B}$ with
$$M^{-1} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}, \qquad M^{-1}\begin{pmatrix}2\\3\\-1\end{pmatrix} = \begin{pmatrix}2 - 3\\3 + 1\\-1\end{pmatrix} = \begin{pmatrix}-1\\4\\-1\end{pmatrix}. ✓$$
$M$ is a **change-of-basis matrix**: its columns are the $\mathcal{B}$-coordinates of the $\mathcal{C}$ vectors. LA-4 builds its whole theory of "a matrix is a linear map written in coordinates" on this.

### Example 4.5: The solution space of $y'' + y = 0$ (mechanics)

A unit mass on a unit spring: $y'' + y = 0$. Let $S$ be the set of all solutions on $\mathbb{R}$. Example 5 of §3.1 showed $S$ is a subspace of $C^\infty(\mathbb{R})$.

**Claim.** $(\cos t, \sin t)$ is a basis of $S$, so $\dim S = 2$, and $[y]_{(\cos, \sin)} = (y(0), y'(0))^\top$.

*Both are in $S$:* $(\cos t)'' = -\cos t$ and $(\sin t)'' = -\sin t$.

*Independent:* suppose $a\cos t + b\sin t = 0$ for all $t$. At $t = 0$: $a = 0$. At $t = \pi/2$: $b = 0$.

*Spanning (this is the uniqueness theorem for this ODE, proved from scratch):* let $y \in S$ and set $z(t) = y(t) - y(0)\cos t - y'(0)\sin t$. Then $z \in S$ (subspace), $z(0) = y(0) - y(0) = 0$, and $z'(0) = y'(0) - y'(0)\cos 0 = 0$ (using $(\cos)' = -\sin$, $(\sin)' = \cos$). Define the energy $\mathcal{E}(t) = z'(t)^2 + z(t)^2$. Then
$$\mathcal{E}'(t) = 2z'z'' + 2zz' = 2z'(z'' + z) = 0,$$
so $\mathcal{E}$ is constant, equal to $\mathcal{E}(0) = 0$. Hence $z(t)^2 \le \mathcal{E}(t) = 0$ for all $t$, so $z \equiv 0$ and $y = y(0)\cos t + y'(0)\sin t$. ∎

The coordinates of a solution are its initial conditions. Changing basis to $(\cos t + \sin t, \cos t - \sin t)$, say, is a $2 \times 2$ matrix exactly as in Example 4.4.

**The forced oscillator and (3.3).** $y'' + y = t$ has the particular solution $y_p = t$ (check: $0 + t = t$). Every solution is $t + a\cos t + b\sin t$: the proof is Theorem 3.9 word for word with "$A$" replaced by the operator $L[y] = y'' + y$, which is linear in the same sense ($L[y_1 + y_2] = L[y_1] + L[y_2]$, $L[cy] = cL[y]$). "General solution $=$ particular $+$ homogeneous" is the complete-solution theorem; $N(L) = S$ is a null space with basis $\cos t, \sin t$. LA-4 makes "linear operator" precise, and LA-7 handles $y'' + y = 0$ as a first-order system whose null-space basis comes from eigenvectors. The project checks this example symbolically with SymPy.

## 5. Why it matters and where it connects

- **LA-4** reads $N(A)$ as the kernel and $C(A)$ as the image of the linear map $\mathbf{x} \mapsto A\mathbf{x}$, and (3.5) becomes the rank–nullity theorem for maps between abstract spaces. The coordinate map of Theorem 3.17 is the isomorphism $V \cong \mathbb{R}^n$; the matrix $M$ of Example 4.4 is a change-of-basis matrix.
- **LA-5** completes the big picture: (3.6) upgrades to "orthogonal complements," and least squares projects $\mathbf{b}$ onto $C(A)$. Exercise B7's result, $N(A^\top A) = N(A)$, is why $A^\top A$ is invertible exactly when the columns are independent, which is the hypothesis of the normal equations.
- **LA-6**: $\det A \ne 0$ iff the columns are independent iff $N(A) = \{\mathbf{0}\}$ (Theorem 3.12).
- **LA-7**: an eigenspace is a null space, $N(A - \lambda I)$, and its dimension (geometric multiplicity) is $n - \operatorname{rank}(A - \lambda I)$. The solution space of $\dot{\mathbf{x}} = A\mathbf{x}$ is an $n$-dimensional vector space of functions, exactly as in Example 4.5.
- **LA-8, LA-9, LA-11** use §3.11. Consequence (ii) of the dimension formula (3.8) drives LA-8's Courant–Fischer theorem and LA-9 Theorem 3.8 (Eckart–Young in the spectral norm). LA-11 splits $\mathbb{C}^n$ as a direct sum of generalized eigenspaces (Definition 3.33) and builds the Jordan basis by concatenating bases of the pieces, which is Theorem 3.34(c).
- **LA-9**: the SVD produces orthonormal bases for all four subspaces at once, and rank-$r$ matrices are sums of $r$ rank-one pieces (Example 4.3). `scipy.linalg.null_space` and `orth`, used as black boxes in the project, are SVD-based, and LA-9 Theorem 3.8 is the reason to prefer them numerically: the distance from $A$ to the nearest matrix of rank $\le k$ is exactly $\sigma_{k+1}$, so singular values say how far $A$ is from each lower rank, which a single pivot compared against a tolerance cannot. Exercise C3(c) shows the two answers disagreeing, and LA-13 §3.2 (catastrophic cancellation) explains why the pivot $(1 + \varepsilon) - 1$ that elimination computes there is inaccurate for tiny $\varepsilon$.
- **F2**: $1, \cos kx, \sin kx$ are independent in the infinite-dimensional space of functions on $[-\pi, \pi]$, and Fourier coefficients are coordinates in that list. Peano's remark about infinite dimension comes due.
- **CM-6**: the solutions of $n$ coupled oscillators form a $2n$-dimensional space (two initial conditions per coordinate, as in Example 4.5). Each normal mode $\mathbf{q}_k$ with frequency $\omega_k$ contributes two basis solutions, $\mathbf{q}_k\cos\omega_k t$ and $\mathbf{q}_k\sin\omega_k t$ (or $\mathbf{q}_k$ and $t\,\mathbf{q}_k$ for a zero-frequency mode), giving a basis of $2n$ functions. **CM-7**: the inertia tensor is a symmetric $3 \times 3$ matrix, an element of the six-dimensional subspace of exercise B3.
- **QM-4**: a quantum state space is a vector space; superposition is closure under addition. Heisenberg's matrix mechanics and Schrödinger's wave mechanics are one space written in two bases (§1). Do not confuse this with the Heisenberg and Schrödinger *pictures* of time evolution in QM-4, which differ in whether the states or the operators carry the time dependence. **QC-1**: a qubit is a two-dimensional complex vector space, $n$ qubits a $2^n$-dimensional one.
- **ST-3**: the design matrix $X$ has $C(X)$ as the set of fittable mean vectors; multicollinearity is $N(X) \ne \{\mathbf{0}\}$, i.e. dependent predictor columns, and the resulting non-uniqueness of the fit is (3.3). **ML-6**: PCA is low-rank approximation of a data matrix.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 drill the subspace test, special solutions, the independence test, and the four subspaces of a small matrix. B1–B8 are the module: B1 (sums and unions of subspaces), B3 (symmetric matrices as a subspace), B5 (every subspace of $\mathbb{R}^3$ is $\{\mathbf{0}\}$, a line, a plane, or $\mathbb{R}^3$; the span and null-space descriptions of a plane), B6 (the $A = CR$ factorization and a second proof that row rank equals column rank), B7 ($N(A^\top A) = N(A)$), and B8 (extend `project.py`) are essential. C1 (direct sums, complements, and why the dimension formula (3.8) has no three-subspace version), C2 (independence of exponentials), and C3 (rank inequalities and numerical rank) go one step beyond.

## 7. Coding project

`project.py` builds everything in §3.3–§3.10 from scratch and checks it against SciPy:

1. **`rref_with_E(A, tol)`** (and the wrapper `rref`, which drops $E$): Gauss–Jordan with partial pivoting on $[A \mid I]$, returning $R$, the pivot columns, and $E$ with $EA = R$. Checked on random matrices (RREF structure, $EA = R$, $E$ invertible, and the same $R$ from two pivoting strategies, as the uniqueness corollary after Theorem 3.28 says).
2. **Four subspaces**: `column_space_basis` (Theorem 3.28), `null_space_basis` (Theorem 3.18), `row_space_basis` (Theorem 3.29), `left_null_space_basis` (Theorem 3.31). Random rank-$r$ matrices are built as products $BC$ with $B$ an $m \times r$ and $C$ an $r \times n$ random matrix, so the rank is known in advance; the dimensions $r, n - r, r, m - r$ are asserted, the orthogonality relations (3.6) and the factorization $A = C\tilde R$ of exercise B6 are checked, and each basis is compared to `scipy.linalg.null_space` / `orth` by checking the two bases span the same subspace.
3. **Examples 4.1–4.4 reproduced exactly**: the matrices $R$ and $E$, the special solutions, the left-null vector $(1, -1, 1)$, the complete solution, the coordinates $(-1, 4, -1)$. Example 4.5 is checked symbolically with SymPy (both basis functions solve the ODE, the Wronskian is $1$, the energy identity $\mathcal{E}' = 2z'(z'' + z)$ holds for an arbitrary function $z$, the coordinates of $a\cos t + b\sin t$ are its initial conditions, `dsolve` returns particular plus homogeneous).
4. **Pictures**: $N(A)$ and $C(A^\top)$ drawn in $\mathbb{R}^3$ for a $1 \times 3$ matrix (plane and line) and a $2 \times 3$ matrix (line and plane), saved to `figures/`.
5. **`is_independent`**, **`extend_to_basis`**, and the rank inequality (3.7) on random products, including a strict case.

**Build it yourself:** before reading the code, write `rref` (loop over columns; find the pivot row; swap; scale to $1$; eliminate above and below; move to the next column) and `null_space_basis` from the special-solution recipe of §3.3. Then compare. The "Your turn" section asks for the symmetric-matrix basis, the vanishing-at-$1$ polynomial basis, and extending $(1, 1, 0), (0, 1, 1)$ to a basis of $\mathbb{R}^3$.

## 8. Where you'll struggle

- **This is the wall.** You go from executing elimination to proving statements about every vector in a set you cannot list. Nobody's intuition carries them through this on its own. The fix is mechanical and boring: write one proof per day, in full sentences, following the §3.0 patterns, for the whole time you are in LA-3 and LA-4. Start with the ones in this module that are three lines long (Theorems 3.5, 3.6, 3.9, 3.13) and rewrite them closed-book until the first and last lines come automatically. Do not skip F1's exercises on quantifiers; "for all $\mathbf{x}$ in $N(A)$" and "there exists $\mathbf{x}_p$" are the whole content of (3.3).
- **Existence versus uniqueness.** "$\mathbf{b}$ is in the span" means *there exist* coefficients. "The set is independent" means coefficients, when they exist, are *unique*. Students conflate them in both directions: proving a set spans by showing the trivial relation is the only one, or proving independence by exhibiting one combination. Theorem 3.15 is the exact bridge, and it only connects uniqueness to independence, never existence. When you are stuck, ask out loud: is this an existence question or a uniqueness question?
- **A subspace must contain $\mathbf{0}$.** The line $y = x + 1$, the set of vectors with $\|\mathbf{x}\| = 1$, the solution set of $A\mathbf{x} = \mathbf{b}$ with $\mathbf{b} \ne \mathbf{0}$: none of them is a subspace, and checking $\mathbf{0}$ first is the fastest way to see it. If $\mathbf{0}$ passes, *then* check the two closures with arbitrary elements, not with two examples you happen to like.
- **Rank is not about rows.** You compute it with row operations, so it feels like a property of the rows. It is $\dim C(A)$, a property of the columns, and *equally* $\dim C(A^\top)$. Theorem 3.29 is the reason; if it feels like a coincidence, redo exercise B6, whose $A = CR$ factorization shows the same $r$ counting columns of $C$ and rows of $R$.
- **$\dim N(A) = n - r$, not $r$, and $n$ is the number of columns.** Three slips happen here. Using $m$ instead of $n$ (the null space lives in $\mathbb{R}^n$, the input space). Counting pivots instead of free columns. Thinking a wide matrix (more columns than rows) might have a trivial null space; it cannot, because $r \le m < n$ forces free columns (Theorem 3.14).
- **Taking the pivot columns of $R$ instead of $A$.** $C(R) \ne C(A)$. $R$ tells you *which* columns; you take them from $A$. The project's `column_space_basis` does this and asserts it against `orth(A)`; break it deliberately and watch the assertion fail.
- **Steinitz.** The exchange proof is the hardest in the module and the first one where the induction hypothesis has structure ("a spanning list of length $n$ containing $\mathbf{u}_1, \dots, \mathbf{u}_j$"). Do the ▶ Try-it-first case with $k = 2$, $n = 3$ by hand with explicit vectors before reading the general proof. If you only ever read it, it will not stick.
- **Theorem 3.35 (dimension of $U + W$)** is harder than it looks because you have to *construct* a basis of $U + W$ from bases of $U$, $W$, and $U \cap W$ and then prove independence by a two-stage argument. Close the notes and rebuild step 4 yourself; budget an hour. Then do exercise C1(c), which shows why the "inclusion–exclusion" pattern stops at two subspaces.

## 9. Check yourself

1. State the two-closure test. Why does closure under scaling by $0$ and $-1$ make axioms (A3) and (A4) automatic?
2. Give a set in $\mathbb{R}^2$ that is closed under addition but not under scalar multiplication, and one closed under scalar multiplication but not addition.
3. Why is $N(A)$ a subspace of $\mathbb{R}^n$ and not of $\mathbb{R}^m$? Where does $C(A)$ live?
4. Write down the special solutions for $R = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 5\end{pmatrix}$ without pencil.
5. Define linear independence in one line with quantifiers. Restate it as a uniqueness statement.
6. Why are any five vectors in $\mathbb{R}^4$ dependent? Two different proofs: one via pivots, one via Steinitz.
7. State the Steinitz exchange lemma. What does it say about the lengths of two bases?
8. The pivot columns of $R$ are $\mathbf{e}_1, \dots, \mathbf{e}_r$. Why does that make the pivot columns of $A$ a basis of $C(A)$, even though $C(A) \ne C(R)$?
9. Dimensions of the four subspaces of a $5 \times 7$ matrix of rank $3$. Which two live in $\mathbb{R}^7$?
10. Give the one-line reason $N(A) \perp C(A^\top)$.
11. Why is $\operatorname{rank}(AB) \le \operatorname{rank} B$? (Which null space contains which?)
12. Why does the energy argument in Example 4.5 prove that $\cos t, \sin t$ *span* the solution space, rather than that they are independent?
13. When is $V = U \oplus W$? Three lines through $\mathbf{0}$ in $\mathbb{R}^2$ meet pairwise only in $\mathbf{0}$; why is their sum not direct?

## 10. Pair with

- **Strang, *Introduction to Linear Algebra*, ch. 3** (§3.1–3.5). The four subspaces are Strang's organizing idea and his §3.5 big-picture figure is the one drawn in §3.9. His worked $3 \times 4$ examples and the "special solutions" recipe are better rehearsed than anything here.
- **MIT 18.06, lectures 5–10** (OCW). Lecture 6 (column space and null space), 7 (special solutions), 8 (complete solution $\mathbf{x}_p + \mathbf{x}_n$), 9 (independence, basis, dimension), 10 (the four subspaces). Watch 9 and 10 after §3.9; Strang draws the picture live and talks through the dimension count.
- **Axler, *Linear Algebra Done Right*, 4th ed., ch. 1–2.** The abstract side done properly: axioms, subspaces, span, independence, the linear dependence lemma (his 2.19), and the proof that independent lists are no longer than spanning lists (his 2.22, the model for Theorem 3.21). Read ch. 2 slowly and do his exercises; this is the proof style LA-4 continues. He uses "list" exactly as this module does, and in these chapters avoids matrices entirely, which is useful discomfort.
- **Hammack, *Book of Proof*, ch. 4–8** (F1's text) if the §3.0 patterns feel foreign: direct proof, contrapositive, contradiction, proving set equality by double inclusion.
