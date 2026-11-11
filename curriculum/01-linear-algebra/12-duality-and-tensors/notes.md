# LA-12: Duality, Bilinear Forms, Tensors and Tensor Products, Kronecker Products, Index Notation

**Time:** ~14 hours · **Prerequisites:** LA-8 (LA-9's SVD is used once, in Theorem 3.16) · **Unlocks:** CM-7, EM-8, QM-4, QM-7, QC-2

> **Conventions in this module.** $V, W$ are finite-dimensional vector spaces over $\mathbb{F} = \mathbb{R}$ or $\mathbb{C}$. **Basis vectors carry lower indices** $\mathbf{e}_1, \dots, \mathbf{e}_n$; **coordinates of a vector carry upper indices**, $\mathbf{v} = \sum_i v^i \mathbf{e}_i$. The dual basis is $\mathbf{e}^1, \dots, \mathbf{e}^n$ (upper), and components of a functional $\varphi$ are $\varphi_i = \varphi(\mathbf{e}_i)$ (lower). **Change of basis:** $\mathbf{e}'_j = \sum_i P_{ij}\mathbf{e}_i$, so the **columns of $P$ are the new basis vectors written in old coordinates.** Then coordinates of vectors change by $\mathbf{v}' = P^{-1}\mathbf{v}$, components of functionals by $\boldsymbol\varphi' = P^\top\boldsymbol\varphi$, matrices of operators by $A' = P^{-1}AP$, matrices of bilinear forms by $B' = P^\top B P$. **Einstein convention:** an index repeated once up and once down is summed; in an orthonormal basis of $\mathbb{R}^n$ we relax this and sum any repeated index, as engineers do. Complex inner products are conjugate-linear in the **first** slot: $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^\dagger\mathbf{v}$ (LA-8 convention). **Stars and daggers:** $V^*$ is the dual space and $T^*$ the dual map, neither of which involves an inner product; the adjoint, which LA-8 wrote as $A^*$, is written $T^\dagger$ here. **Kronecker product:** $(A\otimes B)_{(i,k),(j,l)} = A_{ij}B_{kl}$ with the pair $(i,k)$ ordered row-major, $(i,k)\mapsto i\cdot(\text{rows of }B) + k$; this is `np.kron`. **Qubit 0 is the leftmost tensor factor:** $|01\rangle = |0\rangle\otimes|1\rangle = (0,1,0,0)^\top$. **$\operatorname{vec}(X)$ stacks the columns of $X$** (column-major, NumPy `order="F"`). **Stress:** $\sigma_{ij}$ is the $i$-component of the traction on the face whose outward normal is $\mathbf{e}_j$, $\mathbf{t} = \sigma\mathbf{n}$, tension positive, units MPa. **Rotations** are counterclockwise; $R(\theta) = \exp(\theta K)$ with $K$ the antisymmetric matrix of the axis. $\varepsilon_{123} = +1$.

## 0. Where this sits

You know the spectral theorem for symmetric and Hermitian matrices (LA-8), quadratic forms and positive definiteness (LA-8), and the SVD (LA-9). You have used the inertia tensor and the stress tensor in ME courses as "matrices that rotate by $Q\sigma Q^\top$," and you have probably been told that a tensor is "something that transforms like a tensor," which is not a definition.

This module gives the definition. It builds the dual space $V^*$, then bilinear forms, then the tensor product $V\otimes W$, and shows that the three things called "tensor" in mechanics, differential geometry, and quantum mechanics are one object seen in three coordinate systems. It is the last piece of linear algebra the physics courses need: the inertia tensor and Euler's equations in CM-7, the electromagnetic field tensor in EM-8, Dirac notation and the adjoint in QM-4, composite systems and entanglement in QM-7, and multi-qubit gates in QC-2 all start from here. The appendix on $SO(3)$ and $SU(2)$ is read in full in QM-5 and QC-1.

## 1. The problem that forced this idea

Three communities invented the tensor, separately, because three different problems refused to be a vector.

**Stress, 1822.** Augustin-Louis Cauchy was trying to describe the internal forces in a deformed elastic solid. Pressure in a fluid is a number: the force per area on any small surface is $p$ times the area, pointing along the normal. In a solid under load this fails twice. The force per area (the *traction* $\mathbf{t}$) depends on which way the surface faces, and it is not along the normal: a shaft under torque has shear on its cross-section. Cauchy's result, presented to the Académie in 1822, is that the dependence on the surface is *linear*: there is a $3\times 3$ array $\sigma$ such that the traction on the face with unit normal $\mathbf{n}$ is $\mathbf{t} = \sigma\mathbf{n}$. The argument is the tetrahedron you saw in solid mechanics: balance forces on a shrinking tetrahedron and the only consistent possibility is linearity in $\mathbf{n}$. Nine numbers that take a direction and return a force. Not a vector. Not a scalar. The word for it did not exist yet; Woldemar Voigt gave it the name *tensor* in 1898, working on crystal elasticity.

**Coordinates, 1900.** Gregorio Ricci-Curbastro and Tullio Levi-Civita published *Méthodes de calcul différentiel absolu et leurs applications* in 1900. Their problem was Riemann's: do geometry on a curved space using whatever coordinates are convenient, and make sure the answers do not depend on the coordinates. Their tool was bookkeeping. Attach to each coordinate system an array of components, and specify exactly how the array must change when the coordinates change. An object that is specified this way, in every coordinate system at once, is coordinate-independent by construction. The transformation rule came in two flavors, which they called *covariant* and *contravariant*, and the distinction is the single most confusing thing in the subject, because in the orthonormal Cartesian coordinates engineers use, the two flavors are identical. Einstein learned this calculus from his friend Marcel Grossmann in 1912 and wrote general relativity in it.

**Composite quantum systems, 1932.** John von Neumann's *Mathematische Grundlagen der Quantenmechanik* set quantum mechanics in Hilbert space, and had to answer: if particle 1 lives in a space $V$ and particle 2 lives in $W$, where does the pair live? Not in $V\times W$ (pairs of states), and not in $V\oplus W$ (dimension adds). The right space is the **tensor product** $V\otimes W$, whose dimension *multiplies*, and in which most vectors are not of the form $v\otimes w$. Those non-product vectors are the entangled states. Schrödinger named entanglement in 1935; the mathematics that makes it unavoidable is a dimension count: the product vectors $\mathbf{v}\otimes\mathbf{w}$ fill only a surface of dimension $\dim V + \dim W - 1$ inside the $(\dim V\cdot\dim W)$-dimensional space $V\otimes W$ (three dimensions out of four for two qubits), so almost every vector is not a product. The decomposition that measures how entangled a state is was found by Erhard Schmidt in 1907 for integral equations, decades before anyone needed it for physics.

Three problems: a linear machine that turns normals into forces; a bookkeeping rule for arrays under coordinate change; a vector space whose elements are mostly not products. This module's job is to show you that these are one idea, and to make you fluent enough in each dialect that you can read a mechanics text, a relativity text, and a quantum computing text without translating.

## 2. Intuition

### A matrix is three different things

Take a $2\times 2$ array of numbers $M = \begin{pmatrix} 2 & 1 \\ 1 & 3\end{pmatrix}$. What is it?

**Reading 1: an operator.** $M$ eats a vector and returns a vector: $\mathbf{v}\mapsto M\mathbf{v}$. This is the LA-4 reading. In the language of this module, $M$ is a tensor of **type $(1,1)$**: one slot for a vector in, one vector out.

**Reading 2: a bilinear form.** $M$ eats *two* vectors and returns a *number*: $(\mathbf{u},\mathbf{v})\mapsto\mathbf{u}^\top M\mathbf{v} = 2u^1v^1 + u^1v^2 + u^2v^1 + 3u^2v^2$. A "machine with two input slots and a numerical readout." This is a tensor of **type $(0,2)$**. The kinetic energy $\tfrac12\boldsymbol\omega^\top I\boldsymbol\omega$ of a rigid body is a bilinear form eating $\boldsymbol\omega$ twice.

**Reading 3: a vector in a bigger space.** Reshape $M$ into the 4-vector $(2, 1, 1, 3)^\top$. That is an element of $\mathbb{R}^2\otimes\mathbb{R}^2 \cong \mathbb{R}^4$, the space of two-qubit states if you normalize it. A tensor of **type $(2,0)$**.

Same four numbers. Three readings. The readings differ in exactly one place: **what happens when you change basis.** Under $\mathbf{e}'_j = \sum_i P_{ij}\mathbf{e}_i$, the operator's matrix becomes $P^{-1}MP$, the form's matrix becomes $P^\top MP$, and the reshaped vector's coefficient matrix becomes $P^{-1}M P^{-\top}$. You should guess now that these three laws are the only possible ones, and that the whole covariant/contravariant business is the question of whether a given slot picks up a $P^{-1}$ or a $P^\top$.

### Rulers and arrows

Here is the picture that makes covariant versus contravariant physical. Draw a vector $\mathbf{v}$ in the plane with basis $\mathbf{e}_1, \mathbf{e}_2$. Now double both basis vectors: $\mathbf{e}'_i = 2\mathbf{e}_i$. The arrow $\mathbf{v}$ has not moved, but its coordinates *halve*, because you need half as many of the longer basis vectors to reach the same point. Coordinates of vectors go the *opposite* way to the basis. That is **contravariant.**

Now consider a *measuring device*: a linear functional $\varphi$ that reads off, say, "height above the $\mathbf{e}_1$ axis in units of $\mathbf{e}_2$." Its components are its readings on the basis vectors: $\varphi_i = \varphi(\mathbf{e}_i)$. Double the basis vectors and the readings *double*. Components of functionals go the *same* way as the basis. That is **covariant.**

```
   basis doubled          vector v: coordinates halve   (contravariant, P^{-1})
   e_i -> 2 e_i           functional phi: readings double  (covariant,  P^T)
```

In an orthonormal basis that you rotate (never stretch), $P^{-1} = P^\top$ and the two behaviors coincide. This is why you got through two years of engineering without noticing there were two kinds of index. It is also why the stress tensor "transforms by $Q\sigma Q^\top$" in every ME text with no mention of which index is which: for rotations it does not matter.

### A vector in $V\otimes W$ is a matrix, and its rank is the whole story

Two qubits. Qubit 0 has basis $|0\rangle, |1\rangle$; so does qubit 1. The pair has basis $|00\rangle, |01\rangle, |10\rangle, |11\rangle$, four states, and a general state is $\psi = \sum_{ij} C_{ij}\,|i\rangle\otimes|j\rangle$. Arrange the coefficients as a $2\times 2$ matrix $C$.

The state $|+\rangle\otimes|0\rangle = \tfrac{1}{\sqrt2}(|0\rangle + |1\rangle)\otimes|0\rangle = \tfrac1{\sqrt2}(|00\rangle + |10\rangle)$ has $C = \tfrac1{\sqrt2}\begin{pmatrix}1 & 0\\ 1 & 0\end{pmatrix}$, rank 1.

The Bell state $\tfrac1{\sqrt2}(|00\rangle + |11\rangle)$ has $C = \tfrac1{\sqrt2}\begin{pmatrix}1 & 0\\ 0 & 1\end{pmatrix}$, rank 2.

Guess: a state is a product $a\otimes b$ exactly when its coefficient matrix is an outer product $\mathbf{a}\mathbf{b}^\top$, which is exactly when it has rank 1. And the way to write *any* matrix as a sum of as few rank-1 pieces as possible is the SVD (LA-9). So the SVD of $C$ should tell you how entangled $\psi$ is. That is Theorems 3.15–3.16, and the number of nonzero singular values is the Schmidt rank.

### Three definitions, one object

| Definition | Who uses it | What a type-$(1,1)$ tensor is |
|---|---|---|
| Multilinear map $V^*\times V\to\mathbb{F}$ | Mathematicians (Axler, Halmos) | $(\varphi, \mathbf{v})\mapsto\varphi(T\mathbf{v})$ |
| Element of $V\otimes V^*$ | Quantum mechanics | $\sum_{ij} T^i{}_j\,\mathbf{e}_i\otimes\mathbf{e}^j$ |
| Array $T^i{}_j$ that transforms as $T'^i{}_j = (P^{-1})^i{}_k\,T^k{}_l\,P^l{}_j$ | Physics, engineering (Carroll, Schutz, every ME text) | the matrix of $T$ in each basis |

The first is basis-free. The second is basis-free too, but is built to make "dimension multiplies" obvious. The third is the first with a basis chosen and then the basis thrown away by demanding consistency across all bases. §3.17 proves they are the same, and after that you may use whichever dialect the book in front of you speaks.

## 3. Formal development

### 3.1 The dual space

**Definition 3.1 (Linear functional, dual space).** A **linear functional** on $V$ is a linear map $\varphi: V\to\mathbb{F}$. The set of all of them, with pointwise addition and scaling, is a vector space called the **dual space** $V^* = \mathcal{L}(V, \mathbb{F})$.

Examples you already know, in disguise:

- On $\mathbb{R}^n$: "take the second coordinate," $\varphi(\mathbf{v}) = v^2$.
- On $\mathbb{R}^n$ with a fixed $\mathbf{a}$: "dot with $\mathbf{a}$," $\varphi_{\mathbf{a}}(\mathbf{v}) = \mathbf{a}\cdot\mathbf{v}$. (The first example is this one with $\mathbf{a} = \mathbf{e}_2$.)
- On polynomials of degree $\le n$: "evaluate at $x = 3$," $\varphi(p) = p(3)$; and "integrate from 0 to 1," $\varphi(p) = \int_0^1 p$.
- On $n\times n$ matrices: trace.
- In quantum mechanics: the bra $\langle\phi|$, which eats $|\psi\rangle$ and returns the number $\langle\phi|\psi\rangle$ (§3.6).

A functional is a row vector, in the following precise sense. Fix a basis $\mathbf{e}_1, \dots, \mathbf{e}_n$ and write $\mathbf{v} = \sum_j v^j\mathbf{e}_j$. By linearity
$$\varphi(\mathbf{v}) = \sum_j v^j\,\varphi(\mathbf{e}_j) = \sum_j \varphi_j v^j, \qquad \varphi_j := \varphi(\mathbf{e}_j). \tag{3.1}$$
So $\varphi$ is determined by the $n$ numbers $\varphi_j$, and $\varphi(\mathbf{v})$ is the row $(\varphi_1, \dots, \varphi_n)$ times the column $(v^1, \dots, v^n)^\top$. Lower index on $\varphi_j$, upper on $v^j$, one of each, summed: that is the Einstein convention earning its keep already.

### 3.2 The dual basis

▶ **Try it first.** Given a basis $\mathbf{e}_1,\dots,\mathbf{e}_n$ of $V$, write down $n$ functionals $\mathbf{e}^1, \dots, \mathbf{e}^n$ such that $\mathbf{e}^i$ returns the $i$-th coordinate of its input. Then convince yourself these $n$ functionals are a basis of $V^*$.

**Theorem 3.2 (Dual basis).** Let $\mathbf{e}_1, \dots, \mathbf{e}_n$ be a basis of $V$. There exist unique functionals $\mathbf{e}^1, \dots, \mathbf{e}^n \in V^*$ with
$$\mathbf{e}^i(\mathbf{e}_j) = \delta^i_j = \begin{cases}1 & i = j\\ 0 & i\ne j,\end{cases} \tag{3.2}$$
and they form a basis of $V^*$. Consequently $\dim V^* = \dim V$, and every $\varphi\in V^*$ equals $\varphi = \sum_i \varphi(\mathbf{e}_i)\,\mathbf{e}^i$.

**Derivation.**

1. *Existence.* A linear map is determined by its values on a basis, and any choice of values extends uniquely to a linear map (LA-4). So define $\mathbf{e}^i$ by prescribing $\mathbf{e}^i(\mathbf{e}_j) = \delta^i_j$ on the basis and extending linearly. Explicitly, $\mathbf{e}^i(\sum_j v^j\mathbf{e}_j) = \sum_j v^j\delta^i_j = v^i$: the $i$-th coordinate functional.
2. *Uniqueness.* Two functionals agreeing on a basis agree everywhere, by (3.1).
3. *Independence.* Suppose $\sum_i c_i\mathbf{e}^i = 0$ in $V^*$. Apply both sides to $\mathbf{e}_j$: $\sum_i c_i\delta^i_j = c_j = 0$. All coefficients vanish.
4. *Spanning.* Let $\varphi\in V^*$ and set $\psi = \sum_i\varphi(\mathbf{e}_i)\mathbf{e}^i$. Then $\psi(\mathbf{e}_j) = \sum_i\varphi(\mathbf{e}_i)\delta^i_j = \varphi(\mathbf{e}_j)$ for every $j$, so $\psi = \varphi$ by step 2.
5. $n$ independent spanning vectors: a basis, and $\dim V^* = n$. ∎

■ **What just happened.** The dual basis is "the coordinate functions." Its existence is the LA-4 fact that you may assign a linear map's values on a basis freely. The formula $\varphi = \sum_i\varphi(\mathbf{e}_i)\mathbf{e}^i$ says: a functional's components in the dual basis are its readings on the original basis. That is (3.1) again.

**The dual basis is not the transpose.** Write the basis vectors as the columns of an $n\times n$ matrix $E = (\mathbf{e}_1\ \cdots\ \mathbf{e}_n)$ (coordinates in some reference basis, say the standard one). Write the functional $\mathbf{e}^i$ as a row vector $\mathbf{r}^i$, so that $\mathbf{e}^i(\mathbf{v}) = \mathbf{r}^i\mathbf{v}$. Stack the rows into a matrix $F$ with rows $\mathbf{r}^1, \dots, \mathbf{r}^n$. Condition (3.2) says $\mathbf{r}^i\mathbf{e}_j = \delta^i_j$, i.e.
$$FE = I, \qquad F = E^{-1}. \tag{3.3}$$
**The dual basis vectors are the rows of $E^{-1}$,** not the rows of $E^\top$. They coincide only when $E^{-1} = E^\top$, i.e. when the basis is orthonormal. Example 4.1 computes a case where the dual basis has the right directions but the wrong scale, and one where it has the wrong directions altogether. If you have ever been handed a non-orthogonal basis and tried to read off coordinates by dotting with the basis vectors, this is why it failed.

### 3.3 The double dual and the natural isomorphism

$V^*$ is a vector space, so it has its own dual $V^{**} = (V^*)^*$: functionals on functionals. By Theorem 3.2 applied twice, $\dim V^{**} = \dim V$, so they are isomorphic. But $V$ and $V^*$ are also isomorphic for the same reason, and that isomorphism is *not canonical*: to send $\mathbf{v}$ to a functional you must choose something (a basis, or an inner product). Between $V$ and $V^{**}$ there is an isomorphism that requires no choice at all.

**Theorem 3.3 (Natural isomorphism $V\cong V^{**}$).** For $\mathbf{v}\in V$ define $\Lambda_{\mathbf{v}}\in V^{**}$ by $\Lambda_{\mathbf{v}}(\varphi) = \varphi(\mathbf{v})$ ("evaluate at $\mathbf{v}$"). The map $\Lambda: V\to V^{**}$, $\mathbf{v}\mapsto\Lambda_{\mathbf{v}}$, is a linear isomorphism.

*Proof.* (i) $\Lambda_{\mathbf{v}}$ is linear in $\varphi$: $\Lambda_{\mathbf{v}}(a\varphi + b\psi) = (a\varphi + b\psi)(\mathbf{v}) = a\varphi(\mathbf{v}) + b\psi(\mathbf{v})$. (ii) $\Lambda$ is linear in $\mathbf{v}$: $\Lambda_{a\mathbf{u} + b\mathbf{v}}(\varphi) = \varphi(a\mathbf{u} + b\mathbf{v}) = a\Lambda_{\mathbf{u}}(\varphi) + b\Lambda_{\mathbf{v}}(\varphi)$. (iii) $\Lambda$ is injective: if $\Lambda_{\mathbf{v}} = 0$ then $\varphi(\mathbf{v}) = 0$ for every $\varphi$. If $\mathbf{v}\ne\mathbf{0}$, extend it to a basis $\mathbf{e}_1 = \mathbf{v}, \mathbf{e}_2, \dots$ and take $\varphi = \mathbf{e}^1$; then $\mathbf{e}^1(\mathbf{v}) = 1 \ne 0$. Contradiction, so $\mathbf{v} = \mathbf{0}$. (iv) Injective linear map between spaces of the same finite dimension is an isomorphism (rank–nullity, LA-4). ∎

This is why, in index notation, a vector is "a $(1,0)$ tensor": it is a linear machine that eats one functional and returns a number. The identification costs nothing. In infinite dimensions (QM-4's Hilbert spaces) the statement needs care and is called reflexivity; for $L^2$ it holds.

### 3.4 The dual map, and what the transpose is

**Definition 3.4 (Dual map).** Let $T: V\to W$ be linear. The **dual map** $T^*: W^*\to V^*$ is
$$T^*(\psi) = \psi\circ T, \qquad\text{i.e.}\qquad (T^*\psi)(\mathbf{v}) = \psi(T\mathbf{v}).$$
It goes *backwards*: a functional on the target is pulled back to a functional on the source by first applying $T$. Check it lands in $V^*$: $\psi\circ T$ is a composition of linear maps, hence linear. Check $T^*$ is linear in $\psi$: $(a\psi + b\chi)\circ T = a(\psi\circ T) + b(\chi\circ T)$.

▶ **Try it first.** Pick bases $\mathbf{v}_1..\mathbf{v}_n$ of $V$ and $\mathbf{w}_1..\mathbf{w}_m$ of $W$ with dual bases $\mathbf{v}^j$ and $\mathbf{w}^i$. Let $A$ be the $m\times n$ matrix of $T$: $T\mathbf{v}_j = \sum_i A_{ij}\mathbf{w}_i$. Compute the matrix of $T^*$ with respect to the dual bases by applying $T^*\mathbf{w}^i$ to $\mathbf{v}_j$.

**Theorem 3.5 (The matrix of the dual map is the transpose).** With bases as above, the matrix of $T^*: W^*\to V^*$ with respect to the dual bases $\{\mathbf{w}^i\}$, $\{\mathbf{v}^j\}$ is $A^\top$.

**Derivation.**

1. The matrix $C$ of $T^*$ is defined by $T^*\mathbf{w}^i = \sum_j C_{ji}\,\mathbf{v}^j$ (column $i$ holds the coordinates of the image of the $i$-th basis functional). $C$ is $n\times m$.
2. By Theorem 3.2, the coordinates of a functional in the dual basis are its values on the basis: $C_{ji} = (T^*\mathbf{w}^i)(\mathbf{v}_j)$.
3. By Definition 3.4, $(T^*\mathbf{w}^i)(\mathbf{v}_j) = \mathbf{w}^i(T\mathbf{v}_j) = \mathbf{w}^i\!\left(\sum_k A_{kj}\mathbf{w}_k\right) = \sum_k A_{kj}\,\delta^i_k = A_{ij}$.
4. So $C_{ji} = A_{ij}$, i.e. $C = A^\top$. ∎

■ **What just happened.** The transpose of a matrix has been a bookkeeping operation since LA-1 ("flip rows and columns"). It is the matrix of a specific linear map: *pull functionals back along $T$.* That is what $A^\top$ **is**, with no inner product in sight. Every property of the transpose follows: $(ST)^* = T^*S^*$ because $\psi\circ(S\circ T) = (\psi\circ S)\circ T$, which is $(BA)^\top = A^\top B^\top$ (LA-1 Theorem 3.8). In LA-8 you met the *adjoint*, written $A^*$ there and $T^\dagger$ here (this module reserves $^*$ for the dual map), defined through an inner product; on a real inner product space with orthonormal bases its matrix is also $A^\top$. The two coincide there because Riesz (§3.6) identifies $V^*$ with $V$ isometrically when the basis is orthonormal. On a non-orthonormal basis they differ, and the dual map is the one that is always the plain transpose.

### 3.5 Annihilators and the four subspaces in dual language (brief)

**Definition 3.6.** For a subspace $U\subseteq V$, the **annihilator** is $U^0 = \{\varphi\in V^*: \varphi(\mathbf{u}) = 0\ \forall\mathbf{u}\in U\}$.

**Theorem 3.7.** $\dim U + \dim U^0 = \dim V$.

*Proof.* Choose a basis $\mathbf{e}_1, \dots, \mathbf{e}_k$ of $U$ and extend to a basis $\mathbf{e}_1, \dots, \mathbf{e}_n$ of $V$; take the dual basis. Claim: $U^0 = \operatorname{span}\{\mathbf{e}^{k+1}, \dots, \mathbf{e}^n\}$. ($\supseteq$) Each $\mathbf{e}^i$ with $i > k$ kills $\mathbf{e}_1, \dots, \mathbf{e}_k$, hence kills $U$. ($\subseteq$) If $\varphi\in U^0$ then $\varphi = \sum_i\varphi(\mathbf{e}_i)\mathbf{e}^i$ and $\varphi(\mathbf{e}_i) = 0$ for $i\le k$, so $\varphi$ is a combination of $\mathbf{e}^{k+1}, \dots, \mathbf{e}^n$. These $n - k$ functionals are independent (part of a basis), so $\dim U^0 = n - k$. ∎

**Theorem 3.8.** For $T: V\to W$: $\operatorname{null}T^* = (\operatorname{range}T)^0$ and $\operatorname{range}T^* = (\operatorname{null}T)^0$.

*Proof of the first.* $\psi\in\operatorname{null}T^*$ iff $\psi\circ T = 0$ iff $\psi(T\mathbf{v}) = 0$ for all $\mathbf{v}$ iff $\psi$ kills $\operatorname{range}T$. The second: $\operatorname{range}T^*\subseteq(\operatorname{null}T)^0$ since $(T^*\psi)(\mathbf{v}) = \psi(T\mathbf{v}) = \psi(\mathbf{0}) = 0$ for $\mathbf{v}\in\operatorname{null}T$. The dimensions agree. With $n = \dim V$, $m = \dim W$: rank–nullity for $T^*$, the first statement, and Theorem 3.7 in $W$ give $\dim\operatorname{range}T^* = m - \dim\operatorname{null}T^* = m - \dim(\operatorname{range}T)^0 = \dim\operatorname{range}T$; Theorem 3.7 in $V$ and rank–nullity for $T$ give $\dim(\operatorname{null}T)^0 = n - \dim\operatorname{null}T = \dim\operatorname{range}T$ as well. ∎ (Read through Theorem 3.5, the first chain says $\operatorname{rank}A^\top = \operatorname{rank}A$: LA-3's Theorem 3.29(c), now with a reason.)

Once Riesz (next) identifies functionals with vectors, $(\cdot)^0$ becomes $(\cdot)^\perp$, $T^*$ becomes $A^\top$, and Theorem 3.8 reads $N(A^\top) = C(A)^\perp$ and $C(A^\top) = N(A)^\perp$: the four fundamental subspaces of LA-3, with the dual space as the reason they come in orthogonal pairs.

### 3.6 Riesz representation, and bras

On an inner product space there is a cheap way to make functionals: fix $\mathbf{v}$ and take $\mathbf{u}\mapsto\langle\mathbf{v},\mathbf{u}\rangle$. (Linear in $\mathbf{u}$ because our inner product is linear in the *second* slot.) Riesz says these are all of them.

▶ **Try it first.** Given $\varphi\in V^*$ and an orthonormal basis $\mathbf{e}_1..\mathbf{e}_n$, guess the vector $\mathbf{v}$ with $\varphi = \langle\mathbf{v},\cdot\rangle$ by asking what $\langle\mathbf{v},\mathbf{e}_i\rangle$ must equal. Mind the conjugate.

**Theorem 3.9 (Riesz representation, finite-dimensional).** Let $V$ be a finite-dimensional inner product space. For every $\varphi\in V^*$ there is a unique $\mathbf{v}\in V$ with $\varphi(\mathbf{u}) = \langle\mathbf{v},\mathbf{u}\rangle$ for all $\mathbf{u}$. In an orthonormal basis, $\mathbf{v} = \sum_i\overline{\varphi(\mathbf{e}_i)}\,\mathbf{e}_i$.

**Derivation.**

1. *Existence.* Set $\mathbf{v} = \sum_i\overline{\varphi(\mathbf{e}_i)}\,\mathbf{e}_i$. For $\mathbf{u} = \sum_j u^j\mathbf{e}_j$, conjugate-linearity in the first slot gives $\langle\mathbf{v},\mathbf{u}\rangle = \sum_{i,j}\overline{\overline{\varphi(\mathbf{e}_i)}}\,u^j\langle\mathbf{e}_i,\mathbf{e}_j\rangle = \sum_i\varphi(\mathbf{e}_i)u^i = \varphi(\mathbf{u})$ by (3.1).
2. *Uniqueness.* If $\langle\mathbf{v},\mathbf{u}\rangle = \langle\mathbf{v}',\mathbf{u}\rangle$ for all $\mathbf{u}$, then $\langle\mathbf{v} - \mathbf{v}',\mathbf{u}\rangle = 0$ for all $\mathbf{u}$; take $\mathbf{u} = \mathbf{v} - \mathbf{v}'$ to get $\|\mathbf{v} - \mathbf{v}'\|^2 = 0$. ∎

■ **What just happened.** An inner product is a *choice* of isomorphism $V\to V^*$, $\mathbf{v}\mapsto\langle\mathbf{v},\cdot\rangle$. It is linear over $\mathbb{R}$ and conjugate-linear over $\mathbb{C}$ (the conjugate in step 1 is the physics convention of LA-8 showing up exactly where the roadmap warned it would). Different inner products give different isomorphisms; no inner product, no canonical isomorphism. Compare Theorem 3.3, where $V\cong V^{**}$ needed nothing.

**Bras.** In Dirac notation, a vector is a ket $|\psi\rangle\in V$. The **bra** $\langle\phi|$ is the functional $|\psi\rangle\mapsto\langle\phi|\psi\rangle$, that is, the Riesz partner of $|\phi\rangle$. Bras live in $V^*$, kets in $V$; the bracket $\langle\phi|\psi\rangle$ is "functional applied to vector." In coordinates, $|\psi\rangle$ is a column, $\langle\phi|$ is the conjugate-transposed row $\boldsymbol\phi^\dagger$, and the bra–ket is the $1\times 1$ product. The outer product $|\psi\rangle\langle\phi|$ is a column times a row: an operator, $|\chi\rangle\mapsto|\psi\rangle\langle\phi|\chi\rangle$. In §3.12 you will recognize it as a simple tensor in $V\otimes V^*$. QM-4 builds everything on this.

### 3.7 Bilinear forms and their change of basis

**Definition 3.10.** A **bilinear form** on $V$ is a map $B: V\times V\to\mathbb{F}$ linear in each slot separately: $B(a\mathbf{u} + b\mathbf{u}', \mathbf{v}) = aB(\mathbf{u},\mathbf{v}) + bB(\mathbf{u}',\mathbf{v})$ and likewise in the second slot. Its **matrix** in a basis is $B_{ij} = B(\mathbf{e}_i,\mathbf{e}_j)$.

Expanding both inputs in the basis:
$$B(\mathbf{u},\mathbf{v}) = B\!\left(\sum_i u^i\mathbf{e}_i, \sum_j v^j\mathbf{e}_j\right) = \sum_{i,j}u^iv^j\,B(\mathbf{e}_i,\mathbf{e}_j) = \sum_{i,j}B_{ij}u^iv^j = \mathbf{u}^\top B\,\mathbf{v}, \tag{3.4}$$
where in the last expression $\mathbf{u},\mathbf{v}$ are coordinate columns. Two lower indices on $B_{ij}$, two upper on the coordinates, both summed.

▶ **Try it first.** Change basis by $\mathbf{e}'_j = \sum_i P_{ij}\mathbf{e}_i$. Compute $B'_{jl} = B(\mathbf{e}'_j, \mathbf{e}'_l)$ by bilinearity and read off the matrix equation. Then do the same for an operator $T$ with $T\mathbf{e}_j = \sum_i A_{ij}\mathbf{e}_i$ and see where the difference comes from.

**Theorem 3.11 (Change of basis: forms versus operators).** Under $\mathbf{e}'_j = \sum_i P_{ij}\mathbf{e}_i$,
$$B' = P^\top B P \quad\text{(bilinear form)}, \qquad A' = P^{-1} A P \quad\text{(operator)}. \tag{3.5}$$

**Derivation.**

1. *Form.* $B'_{jl} = B(\mathbf{e}'_j,\mathbf{e}'_l) = B\!\left(\sum_i P_{ij}\mathbf{e}_i, \sum_k P_{kl}\mathbf{e}_k\right) = \sum_{i,k}P_{ij}P_{kl}B_{ik} = \sum_{i,k}(P^\top)_{ji}B_{ik}P_{kl} = (P^\top BP)_{jl}$.
2. *Operator.* $T\mathbf{e}'_j = T\!\left(\sum_i P_{ij}\mathbf{e}_i\right) = \sum_i P_{ij}\sum_k A_{ki}\mathbf{e}_k$. This is expressed in the *old* basis; the matrix $A'$ needs it in the *new* basis. Invert the basis relation: $\mathbf{e}_k = \sum_l (P^{-1})_{lk}\mathbf{e}'_l$. Substitute: $T\mathbf{e}'_j = \sum_l\Big[\sum_{k,i}(P^{-1})_{lk}A_{ki}P_{ij}\Big]\mathbf{e}'_l$, so $A'_{lj} = (P^{-1}AP)_{lj}$. ∎

■ **What just happened.** Look at where the two derivations diverge. For the form, *both* inputs are vectors fed into slots, and each slot contributes a factor $P$ from expanding $\mathbf{e}'$ in terms of $\mathbf{e}$; nothing needs to be re-expressed, because the output is a number. For the operator, the input contributes a $P$ the same way, but the *output is a vector*, and re-expressing it in the new basis costs a $P^{-1}$. **A slot that eats a vector picks up $P$ (and is called covariant); a slot that produces a vector picks up $P^{-1}$ (contravariant).** That one sentence is the entire covariant/contravariant distinction. Everything in §3.15–3.17 is this sentence with more indices.

Consequences worth memorizing. A matrix $A'$ similar to $A$ has the same eigenvalues, trace, and determinant. A matrix $B' = P^\top BP$ **congruent** to $B$ has in general different eigenvalues and trace, and $\det B' = (\det P)^2\det B$. Example 4.2 does this with numbers: the same $2\times 2$ array, the same $P$, two different answers. When $P$ is orthogonal, $P^\top = P^{-1}$, the two laws coincide, and the eigenvalues of a form become meaningful (they are the principal stresses when $B = \sigma$).

### 3.8 Symmetric forms, quadratic forms, and Sylvester's law of inertia

A bilinear form is **symmetric** if $B(\mathbf{u},\mathbf{v}) = B(\mathbf{v},\mathbf{u})$, equivalently $B_{ij} = B_{ji}$. Its **quadratic form** is $q(\mathbf{v}) = B(\mathbf{v},\mathbf{v}) = \mathbf{v}^\top B\mathbf{v}$ (LA-8). The form is recoverable from $q$ by polarization, $B(\mathbf{u},\mathbf{v}) = \tfrac12[q(\mathbf{u}+\mathbf{v}) - q(\mathbf{u}) - q(\mathbf{v})]$, which you verify by expanding $q(\mathbf{u} + \mathbf{v})$ with (3.4) and symmetry.

The spectral theorem (LA-8) says a real symmetric $B$ is $Q D Q^\top$ with $Q$ orthogonal, $D$ diagonal. Since $Q^\top = Q^{-1}$, this is both a similarity *and* a congruence: in the orthonormal eigenbasis the form is $q = \sum_i\lambda_i (v'^i)^2$. But congruence allows more than orthogonal $P$. Scaling each eigenvector by $|\lambda_i|^{-1/2}$ (when $\lambda_i \ne 0$) gives $P = QS$ with $S = \operatorname{diag}(|\lambda_i|^{-1/2}\text{ or }1)$ and
$$P^\top BP = \operatorname{diag}(\underbrace{1,\dots,1}_{p}, \underbrace{-1,\dots,-1}_{q}, \underbrace{0,\dots,0}_{z}). \tag{3.6}$$
So the eigenvalues of a form are *not* invariant under change of basis, but something is:

**Theorem 3.12 (Sylvester's law of inertia).** If a real symmetric bilinear form has the diagonal matrix (3.6) in one basis and $\operatorname{diag}(1,\dots,1,-1,\dots,-1,0,\dots,0)$ with counts $(p',q',z')$ in another, then $p = p'$, $q = q'$, $z = z'$. The triple $(p, q, z)$ is the **signature** of the form.

▶ **Try it first.** Let $U_+$ be the span of the $p$ basis vectors on which the first diagonalization is $+1$, and $W_{\le 0}$ the span of the $q' + z'$ vectors on which the second is $-1$ or $0$. Show $U_+\cap W_{\le 0} = \{\mathbf{0}\}$ and count dimensions.

**Derivation.**

1. On $U_+$ the form is positive definite: for $\mathbf{u} = \sum_{i\le p}c^i\mathbf{e}_i \ne\mathbf{0}$, $q(\mathbf{u}) = \sum_{i\le p}(c^i)^2 > 0$.
2. On $W_{\le 0}$ it is negative semidefinite: for $\mathbf{w} = \sum_{j > p'}d^j\mathbf{e}'_j$, $q(\mathbf{w}) = -\sum_{p' < j\le p'+q'}(d^j)^2 \le 0$.
3. A vector in both would satisfy $q > 0$ and $q\le 0$ unless it is $\mathbf{0}$. So $U_+\cap W_{\le 0} = \{\mathbf{0}\}$, and $\dim U_+ + \dim W_{\le 0}\le n$: $p + (n - p')\le n$, i.e. $p\le p'$.
4. Swap the roles of the two bases: $p'\le p$. So $p = p'$.
5. $p + q = \operatorname{rank}B$ in either basis, and rank is invariant under congruence ($P$ invertible). So $q = q'$ and $z = z'$. ∎

■ **What just happened.** Under similarity, $n$ numbers (the eigenvalues) are invariant. Under congruence, only three integers survive. That is the right amount of information for a bilinear form: a form is determined up to change of basis by its signature. Positive definite (LA-8) is $(n, 0, 0)$. The Minkowski metric of special relativity (EM-8) has signature $(1,3,0)$ or $(3,1,0)$ depending on the author, and *no* change of coordinates can make it positive definite: Sylvester says so.

### 3.9 Sesquilinear forms (the complex case)

Over $\mathbb{C}$ a bilinear form cannot be positive definite ($B(i\mathbf{v}, i\mathbf{v}) = -B(\mathbf{v},\mathbf{v})$), so inner products are not bilinear. A **sesquilinear form** is conjugate-linear in the first slot and linear in the second (our convention), with matrix $B_{ij} = B(\mathbf{e}_i,\mathbf{e}_j)$ and $B(\mathbf{u},\mathbf{v}) = \mathbf{u}^\dagger B\mathbf{v}$. Re-run the derivation of Theorem 3.11: the first slot contributes $\overline{P_{ij}}$, so
$$B' = P^\dagger BP.$$
The form is **Hermitian** if $B(\mathbf{u},\mathbf{v}) = \overline{B(\mathbf{v},\mathbf{u})}$, i.e. $B = B^\dagger$. Sylvester's law holds verbatim with $P^\dagger$ in place of $P^\top$ (the same proof, since $q(\mathbf{v}) = \mathbf{v}^\dagger B\mathbf{v}$ is real for Hermitian $B$). An inner product is a Hermitian form of signature $(n,0,0)$.

### 3.10 Multilinear maps

**Definition 3.13.** A map $T: V_1\times\cdots\times V_k\to\mathbb{F}$ is **$k$-linear** (multilinear) if it is linear in each slot with the others held fixed. The set of them is a vector space $\mathcal{L}(V_1,\dots,V_k;\mathbb{F})$.

The same expansion as (3.4) shows $T$ is determined by its values on tuples of basis vectors, $T_{i_1\cdots i_k} = T(\mathbf{e}^{(1)}_{i_1}, \dots, \mathbf{e}^{(k)}_{i_k})$, and these $\prod_a\dim V_a$ numbers may be chosen freely (define $T$ on basis tuples and extend multilinearly). Hence
$$\dim\mathcal{L}(V_1,\dots,V_k;\mathbb{F}) = \prod_{a=1}^k\dim V_a. \tag{3.7}$$
A basis: the products of dual basis vectors $(\mathbf{v}_1,\dots,\mathbf{v}_k)\mapsto\mathbf{e}^{i_1}(\mathbf{v}_1)\cdots\mathbf{e}^{i_k}(\mathbf{v}_k)$. Dimension *multiplies*. This is the shadow of the tensor product, which is what the next section builds.

### 3.11 The tensor product $V\otimes W$

**The universal property, in one paragraph.** We want a vector space $V\otimes W$ and a bilinear map $\otimes: V\times W\to V\otimes W$, $(\mathbf{v},\mathbf{w})\mapsto\mathbf{v}\otimes\mathbf{w}$, such that *every* bilinear map $\beta: V\times W\to U$ (any target $U$) factors uniquely as $\beta = \tilde\beta\circ\otimes$ with $\tilde\beta: V\otimes W\to U$ **linear**. In words: $V\otimes W$ is the space that turns bilinear maps on $V\times W$ into linear maps on one space, and it does so with no extra junk (uniqueness) and no omissions (existence). Any two spaces with this property are isomorphic by a unique isomorphism compatible with $\otimes$, which is why one speaks of "the" tensor product. You will not feel the force of this definition until you need to define a map *out of* $V\otimes W$ (as in §3.14 for $A\otimes B$) and realize you only have to say what it does to pairs $(\mathbf{v},\mathbf{w})$, check bilinearity, and you are done.

**The concrete construction.** Let $\mathcal{F}$ be the (enormous) vector space of all finite formal linear combinations of symbols $[\mathbf{v},\mathbf{w}]$, one symbol for each pair $(\mathbf{v},\mathbf{w})\in V\times W$, with no relations between them. Let $\mathcal{R}\subseteq\mathcal{F}$ be the subspace spanned by all elements of the four forms
$$[\mathbf{v}_1 + \mathbf{v}_2, \mathbf{w}] - [\mathbf{v}_1,\mathbf{w}] - [\mathbf{v}_2,\mathbf{w}], \quad [c\mathbf{v},\mathbf{w}] - c[\mathbf{v},\mathbf{w}], \quad [\mathbf{v},\mathbf{w}_1 + \mathbf{w}_2] - [\mathbf{v},\mathbf{w}_1] - [\mathbf{v},\mathbf{w}_2], \quad [\mathbf{v}, c\mathbf{w}] - c[\mathbf{v},\mathbf{w}].$$
Define $V\otimes W = \mathcal{F}/\mathcal{R}$ and $\mathbf{v}\otimes\mathbf{w} = [\mathbf{v},\mathbf{w}] + \mathcal{R}$. Quotienting by $\mathcal{R}$ *forces* $\otimes$ to be bilinear: $(\mathbf{v}_1 + \mathbf{v}_2)\otimes\mathbf{w} = \mathbf{v}_1\otimes\mathbf{w} + \mathbf{v}_2\otimes\mathbf{w}$ because their difference lies in $\mathcal{R}$, and so on. Elements of the form $\mathbf{v}\otimes\mathbf{w}$ are **simple** (or product, or decomposable) tensors; a general element is a finite sum of them.

*Universal property holds for this construction.* Given bilinear $\beta: V\times W\to U$, define $\hat\beta: \mathcal{F}\to U$ on the basis of symbols by $\hat\beta([\mathbf{v},\mathbf{w}]) = \beta(\mathbf{v},\mathbf{w})$ and extend linearly. Bilinearity of $\beta$ says exactly that $\hat\beta$ vanishes on each of the four generator types of $\mathcal{R}$, hence on $\mathcal{R}$. So $\hat\beta$ descends to a linear $\tilde\beta: \mathcal{F}/\mathcal{R}\to U$ with $\tilde\beta(\mathbf{v}\otimes\mathbf{w}) = \beta(\mathbf{v},\mathbf{w})$. Uniqueness: the simple tensors span $V\otimes W$, so a linear map is determined by its values on them.

▶ **Try it first.** Let $\mathbf{e}_1..\mathbf{e}_n$ be a basis of $V$ and $\mathbf{f}_1..\mathbf{f}_m$ a basis of $W$. Show every simple tensor $\mathbf{v}\otimes\mathbf{w}$ is a combination of the $\mathbf{e}_i\otimes\mathbf{f}_j$. Then, to show the $nm$ elements $\mathbf{e}_i\otimes\mathbf{f}_j$ are independent, use the universal property with a cleverly chosen bilinear $\beta$ into $\mathbb{F}$.

**Theorem 3.14 (Basis of the tensor product).** $\{\mathbf{e}_i\otimes\mathbf{f}_j: 1\le i\le n, 1\le j\le m\}$ is a basis of $V\otimes W$. Hence $\dim(V\otimes W) = (\dim V)(\dim W)$.

**Derivation.**

1. *Spanning.* $\mathbf{v}\otimes\mathbf{w} = \left(\sum_i v^i\mathbf{e}_i\right)\otimes\left(\sum_j w^j\mathbf{f}_j\right) = \sum_{i,j}v^iw^j\,\mathbf{e}_i\otimes\mathbf{f}_j$ by bilinearity of $\otimes$. Simple tensors span, so the $\mathbf{e}_i\otimes\mathbf{f}_j$ span.
2. *Independence.* Fix $(k,l)$. The map $\beta_{kl}: V\times W\to\mathbb{F}$, $\beta_{kl}(\mathbf{v},\mathbf{w}) = \mathbf{e}^k(\mathbf{v})\,\mathbf{f}^l(\mathbf{w}) = v^kw^l$, is bilinear (product of a linear function of $\mathbf{v}$ and a linear function of $\mathbf{w}$). By the universal property there is a linear $\tilde\beta_{kl}: V\otimes W\to\mathbb{F}$ with $\tilde\beta_{kl}(\mathbf{e}_i\otimes\mathbf{f}_j) = \delta^k_i\delta^l_j$.
3. Suppose $\sum_{i,j}c^{ij}\,\mathbf{e}_i\otimes\mathbf{f}_j = 0$. Apply $\tilde\beta_{kl}$: $\sum_{i,j}c^{ij}\delta^k_i\delta^l_j = c^{kl} = 0$. Every coefficient vanishes. ∎

■ **What just happened.** Spanning was bilinearity. Independence needed the universal property, because the construction is a *quotient*, and in a quotient you must prove things are *not* zero; the way to do it is to exhibit a linear map that does not kill them. The maps $\tilde\beta_{kl}$ are the "coordinate functions" on $V\otimes W$: the dual basis of $V\otimes W$ is $\{\mathbf{e}^k\otimes\mathbf{f}^l\}$. Also note that $\dim(V\otimes W) = nm$, while the simple tensors $\mathbf{v}\otimes\mathbf{w}$ form only an $(n + m - 1)$-dimensional surface, not a subspace: the pair $(\mathbf{v},\mathbf{w})$ has $n + m$ coordinates, and one of them is wasted because $(c\mathbf{v})\otimes\mathbf{w} = \mathbf{v}\otimes(c\mathbf{w})$. For two qubits, $4 > 3$; for two five-qubit registers, $1024 \gg 63$. That gap is where entanglement lives. (The project checks the count by computing the rank of the Jacobian of $(\mathbf{v},\mathbf{w})\mapsto\mathbf{v}\otimes\mathbf{w}$, which is $n + m - 1$.)

### 3.12 The coefficient matrix; simple tensors are rank-1 matrices

By Theorem 3.14 every $t\in V\otimes W$ is uniquely
$$t = \sum_{i,j}C_{ij}\,\mathbf{e}_i\otimes\mathbf{f}_j, \tag{3.8}$$
and the $n\times m$ array $C$ is the **coefficient matrix** of $t$. The map $t\mapsto C$ is a linear isomorphism $V\otimes W\to\mathbb{F}^{n\times m}$. There is a basis-free version: $V\otimes W\cong\operatorname{Hom}(V^*, W)$ via $\mathbf{v}\otimes\mathbf{w}\mapsto(\varphi\mapsto\varphi(\mathbf{v})\,\mathbf{w})$. (Bilinear in $(\mathbf{v},\mathbf{w})$, so it extends by the universal property; it sends $\mathbf{e}_i\otimes\mathbf{f}_j$ to the map $\mathbf{e}^k\mapsto\delta^k_i\mathbf{f}_j$, whose matrix is the $(j,i)$ unit matrix, so the basis goes to a basis.) If $V$ has an inner product, Riesz turns $V^*$ into $V$ and $V\otimes W\cong\operatorname{Hom}(V,W)$: the tensor $\mathbf{v}\otimes\mathbf{w}$ becomes the operator $\mathbf{u}\mapsto\langle\mathbf{v},\mathbf{u}\rangle\mathbf{w}$, which is the outer product $\mathbf{w}\mathbf{v}^\dagger$. That is Dirac's $|w\rangle\langle v|$.

**Theorem 3.15 (Simple tensors have rank-1 coefficient matrices).** $t\in V\otimes W$ is a simple tensor $\mathbf{v}\otimes\mathbf{w}$ if and only if its coefficient matrix has rank $\le 1$.

*Proof.* ($\Rightarrow$) By step 1 of Theorem 3.14, $\mathbf{v}\otimes\mathbf{w}$ has $C_{ij} = v^iw^j$, so $C = \mathbf{v}\mathbf{w}^\top$ (coordinate columns), an outer product, rank $\le 1$. ($\Leftarrow$) A rank-$\le 1$ matrix has a one-dimensional (or zero) column space, so every column is a multiple of one column $\mathbf{a}$: $C = \mathbf{a}\mathbf{b}^\top$ for some $\mathbf{b}$ (LA-3). Then $t = \sum_{i,j}a^ib^j\,\mathbf{e}_i\otimes\mathbf{f}_j = \left(\sum_i a^i\mathbf{e}_i\right)\otimes\left(\sum_j b^j\mathbf{f}_j\right)$ by bilinearity. ∎

So the question "is this two-particle state a product state?" is the question "is this matrix rank 1?", and for $2\times 2$ matrices that is "is the determinant zero?"

### 3.13 The Schmidt decomposition

▶ **Try it first.** Take the SVD of the coefficient matrix, $C = U\Sigma V^\dagger = \sum_k\sigma_k\mathbf{u}_k\mathbf{v}_k^\dagger$ (LA-9), and substitute into (3.8). Collect terms.

**Theorem 3.16 (Schmidt decomposition).** Let $V, W$ be inner product spaces of dimensions $n, m$ and $t\in V\otimes W$. There exist orthonormal sets $\{\mathbf{a}_k\}\subset V$, $\{\mathbf{b}_k\}\subset W$ and numbers $\sigma_1\ge\sigma_2\ge\cdots\ge\sigma_r > 0$, $r\le\min(n,m)$, with
$$t = \sum_{k=1}^r\sigma_k\,\mathbf{a}_k\otimes\mathbf{b}_k. \tag{3.9}$$
The $\sigma_k$ are the nonzero singular values of the coefficient matrix of $t$ in orthonormal bases, $r = \operatorname{rank}C$ is the **Schmidt rank**, and $\sum_k\sigma_k^2 = \|t\|^2$.

**Derivation.**

1. Take orthonormal bases $\mathbf{e}_i$, $\mathbf{f}_j$ and the coefficient matrix $C$ from (3.8). SVD (LA-9): $C = U\Sigma V^\dagger$, i.e. $C_{ij} = \sum_{k=1}^r\sigma_k U_{ik}\overline{V_{jk}}$ with $U, V$ having orthonormal columns.
2. Substitute: $t = \sum_{i,j}\sum_k\sigma_kU_{ik}\overline{V_{jk}}\,\mathbf{e}_i\otimes\mathbf{f}_j = \sum_k\sigma_k\left(\sum_iU_{ik}\mathbf{e}_i\right)\otimes\left(\sum_j\overline{V_{jk}}\,\mathbf{f}_j\right)$ by bilinearity.
3. Set $\mathbf{a}_k = \sum_iU_{ik}\mathbf{e}_i$ (the $k$-th column of $U$ in the $\mathbf{e}$ basis) and $\mathbf{b}_k = \sum_j\overline{V_{jk}}\mathbf{f}_j$. Orthonormal columns of $U$ in an orthonormal basis give orthonormal vectors; the same for $\overline V$, since conjugating every entry preserves $\overline{V}^\dagger\overline V = \overline{V^\dagger V} = I$.
4. The inner product on $V\otimes W$ is defined on simple tensors by $\langle\mathbf{v}\otimes\mathbf{w},\mathbf{v}'\otimes\mathbf{w}'\rangle = \langle\mathbf{v},\mathbf{v}'\rangle\langle\mathbf{w},\mathbf{w}'\rangle$ and extended; with this, $\{\mathbf{e}_i\otimes\mathbf{f}_j\}$ is orthonormal and $\|t\|^2 = \sum_{ij}|C_{ij}|^2 = \|C\|_F^2 = \sum_k\sigma_k^2$ (LA-9). ∎

■ **What just happened.** The SVD, which in LA-9 was about low-rank approximation of data matrices, is the structure theorem for vectors in a tensor product. For a normalized two-particle state, $\sum_k\sigma_k^2 = 1$, and the $\sigma_k^2$ are the probabilities of finding the pair in the $k$-th matched product $\mathbf{a}_k\otimes\mathbf{b}_k$. **A two-qubit state is entangled (not a product) if and only if its Schmidt rank is 2,** by Theorem 3.15. For the Bell state, $C = I/\sqrt2$ has singular values $(1/\sqrt2, 1/\sqrt2)$, rank 2: entangled. For $|+\rangle|0\rangle$, $C = \tfrac1{\sqrt2}\begin{pmatrix}1&0\\1&0\end{pmatrix}$ has singular values $(1, 0)$, rank 1: a product, as it was built to be. Example 4.3 does a less obvious case. The quantity $S = -\sum_k\sigma_k^2\log\sigma_k^2$ is the **entanglement entropy**; it is $0$ for products and $\log 2$ for Bell states; QM-7 and QC-7 meet it again as the von Neumann entropy of the reduced density matrix $CC^\dagger$ (§5).

### 3.14 Kronecker products

If $A: V\to V'$ and $B: W\to W'$ are linear, the map $(\mathbf{v},\mathbf{w})\mapsto A\mathbf{v}\otimes B\mathbf{w}$ is bilinear, so by the universal property there is a unique linear $A\otimes B: V\otimes W\to V'\otimes W'$ with
$$(A\otimes B)(\mathbf{v}\otimes\mathbf{w}) = A\mathbf{v}\otimes B\mathbf{w}. \tag{3.10}$$
(This is the first place the universal property did real work: without it, you would have to check that defining $A\otimes B$ on simple tensors is consistent across the many ways of writing a tensor as a sum of simple ones.)

▶ **Try it first.** Apply (3.10) to a basis vector $\mathbf{e}_j\otimes\mathbf{f}_l$ and expand $A\mathbf{e}_j = \sum_i A_{ij}\mathbf{e}_i$, $B\mathbf{f}_l = \sum_k B_{kl}\mathbf{f}_k$. Read off the matrix entry in row $(i,k)$, column $(j,l)$.

**Derivation of the matrix.**

1. $(A\otimes B)(\mathbf{e}_j\otimes\mathbf{f}_l) = \left(\sum_i A_{ij}\mathbf{e}_i\right)\otimes\left(\sum_k B_{kl}\mathbf{f}_k\right) = \sum_{i,k}A_{ij}B_{kl}\,(\mathbf{e}_i\otimes\mathbf{f}_k)$.
2. The matrix of a linear map has, in column $(j,l)$, the coordinates of the image of basis vector $(j,l)$. So
$$(A\otimes B)_{(i,k),(j,l)} = A_{ij}B_{kl}. \tag{3.11}$$
3. Order the pairs row-major: $(i,k)\mapsto i\cdot m_B + k$ where $m_B$ is the number of rows of $B$ (zero-based), and likewise for columns. Then for fixed $(i,j)$, the entries with $k, l$ ranging form the block $A_{ij}B$, and
$$A\otimes B = \begin{pmatrix}A_{11}B & A_{12}B & \cdots\\ A_{21}B & A_{22}B & \cdots\\ \vdots & & \ddots\end{pmatrix}. \tag{3.12}$$

■ **What just happened.** The matrix of $A\otimes B$ is forced by (3.10) and the choice of how to order the pair index $(i,k)$. This is the **Kronecker product**, and it is exactly `np.kron(A, B)`. The ordering convention is the one that makes "qubit 0 is the leftmost factor" true: $|01\rangle = \mathbf{e}_0\otimes\mathbf{e}_1$ is the pair $(0,1)\mapsto 0\cdot 2 + 1 = 1$, the second standard basis vector $(0,1,0,0)^\top$. Reverse the ordering and every two-qubit gate in QC-2 comes out with its blocks swapped. The vectorized version of (3.8) is the same statement: `psi = C.reshape(-1)` (row-major) and `C = psi.reshape(n, m)`.

**Theorem 3.17 (Kronecker identities).** For matrices of compatible sizes:
(a) $(A\otimes B)(C\otimes D) = AC\otimes BD$; (b) $(A\otimes B)^\top = A^\top\otimes B^\top$; (c) $(A\otimes B)^{-1} = A^{-1}\otimes B^{-1}$ when both inverses exist; (d) $\operatorname{tr}(A\otimes B) = \operatorname{tr}A\cdot\operatorname{tr}B$; (e) if $A$ is $n\times n$ with eigenvalues $\lambda_1..\lambda_n$ and $B$ is $m\times m$ with eigenvalues $\mu_1..\mu_m$ (with multiplicity), the $nm$ eigenvalues of $A\otimes B$ are the products $\lambda_i\mu_k$; (f) $\det(A\otimes B) = (\det A)^m(\det B)^n$.

*Proofs.* (a) Both sides are linear maps on $V\otimes W$; it suffices to check them on simple tensors, which span: $(A\otimes B)(C\otimes D)(\mathbf{v}\otimes\mathbf{w}) = (A\otimes B)(C\mathbf{v}\otimes D\mathbf{w}) = AC\mathbf{v}\otimes BD\mathbf{w} = (AC\otimes BD)(\mathbf{v}\otimes\mathbf{w})$. (b) From (3.11): $((A\otimes B)^\top)_{(i,k),(j,l)} = (A\otimes B)_{(j,l),(i,k)} = A_{ji}B_{lk} = (A^\top)_{ij}(B^\top)_{kl}$. (c) By (a), $(A\otimes B)(A^{-1}\otimes B^{-1}) = AA^{-1}\otimes BB^{-1} = I\otimes I = I$. (d) The diagonal entries are $(i,k),(i,k)$: $\sum_{i,k}A_{ii}B_{kk} = \left(\sum_iA_{ii}\right)\left(\sum_kB_{kk}\right)$. (e) If $A\mathbf{x} = \lambda\mathbf{x}$ and $B\mathbf{y} = \mu\mathbf{y}$ then by (3.10), $(A\otimes B)(\mathbf{x}\otimes\mathbf{y}) = \lambda\mathbf{x}\otimes\mu\mathbf{y} = \lambda\mu\,(\mathbf{x}\otimes\mathbf{y})$. When $A$ and $B$ are diagonalizable with eigenbases $\{\mathbf{x}_i\}$, $\{\mathbf{y}_k\}$, Theorem 3.14 (with these as the bases) says the $\mathbf{x}_i\otimes\mathbf{y}_k$ form a basis of eigenvectors, so these $nm$ products are all the eigenvalues with the right multiplicities. For non-diagonalizable $A, B$, use Schur (LA-10): $A = Q T Q^\dagger$, $B = R S R^\dagger$ with $T, S$ upper triangular; by (a), $A\otimes B = (Q\otimes R)(T\otimes S)(Q\otimes R)^\dagger$, and $T\otimes S$ is upper triangular (check with (3.12): blocks below the diagonal are $T_{ij}S = 0$, and the diagonal blocks $T_{ii}S$ are upper triangular) with diagonal entries $T_{ii}S_{kk} = \lambda_i\mu_k$. (f) Product of the eigenvalues in (e): $\prod_{i,k}\lambda_i\mu_k = \prod_i\lambda_i^m\prod_k\mu_k^n$. ∎

**Theorem 3.18 (The vec trick).** With $\operatorname{vec}$ stacking columns,
$$\operatorname{vec}(AXB) = (B^\top\otimes A)\operatorname{vec}(X). \tag{3.13}$$

▶ **Try it first.** Column $j$ of $AXB$ is $A$ times column $j$ of $XB$, and column $j$ of $XB$ is $\sum_k B_{kj}\mathbf{x}_k$ where $\mathbf{x}_k$ is column $k$ of $X$. Write this as a block-row of a block matrix acting on the stacked columns.

*Proof.* Column $j$ of $AXB$ is $A(XB)\mathbf{e}_j = AX\mathbf{b}_j = A\sum_kB_{kj}\mathbf{x}_k = \sum_kB_{kj}A\mathbf{x}_k = \begin{pmatrix}B_{1j}A & B_{2j}A & \cdots\end{pmatrix}\operatorname{vec}(X)$. By (3.12), the $j$-th block-row of $B^\top\otimes A$ is $\begin{pmatrix}(B^\top)_{j1}A & (B^\top)_{j2}A & \cdots\end{pmatrix} = \begin{pmatrix}B_{1j}A & B_{2j}A & \cdots\end{pmatrix}$. Stacking the columns $j = 1, 2, \dots$ gives (3.13). ∎

■ **What just happened.** Left-multiplication by $A$ and right-multiplication by $B$ are linear maps on the space of matrices, and (3.13) is their matrix once you flatten $X$ into a column. Why you care: a matrix equation like Sylvester's $AX + XB = C$ is linear in $X$, and (3.13) turns it into an ordinary linear system $(I\otimes A + B^\top\otimes I)\operatorname{vec}X = \operatorname{vec}C$ that `np.linalg.solve` handles. LA-14 uses (3.13) to write the derivative of $X\mapsto AXB$ as a single matrix. **Trap:** NumPy's default `ravel()` is row-major. In row-major vectorization the identity reads $\operatorname{vec}_r(AXB) = (A\otimes B^\top)\operatorname{vec}_r(X)$; the project checks both.

### 3.15 Tensors of type $(p, q)$ and index notation

**Definition 3.19.** A **tensor of type $(p,q)$** on $V$ is a multilinear map
$$T: \underbrace{V^*\times\cdots\times V^*}_{p}\times\underbrace{V\times\cdots\times V}_{q}\to\mathbb{F}.$$
Its **components** in a basis are its values on basis vectors and dual basis vectors,
$$T^{i_1\cdots i_p}{}_{j_1\cdots j_q} = T(\mathbf{e}^{i_1},\dots,\mathbf{e}^{i_p},\mathbf{e}_{j_1},\dots,\mathbf{e}_{j_q}), \tag{3.14}$$
with the $p$ **contravariant** indices up and the $q$ **covariant** indices down.

Check this against what you know. A $(0,0)$ tensor is a scalar. A $(0,1)$ tensor is a functional, components $\varphi_j = \varphi(\mathbf{e}_j)$. A $(1,0)$ tensor eats a functional and returns a number: by Theorem 3.3 that is a vector $\mathbf{v}$ acting as $\varphi\mapsto\varphi(\mathbf{v})$, with components $v^i = \mathbf{e}^i(\mathbf{v})$, the coordinates. A $(0,2)$ tensor is a bilinear form, $B_{ij}$. A $(1,1)$ tensor $T(\varphi,\mathbf{v})$ is an operator in disguise: given $S: V\to V$, set $T(\varphi,\mathbf{v}) = \varphi(S\mathbf{v})$, components $T^i{}_j = \mathbf{e}^i(S\mathbf{e}_j) = A_{ij}$, the matrix of $S$ with the row index up and the column index down; conversely every $(1,1)$ tensor arises this way (fix $\mathbf{v}$; $T(\cdot,\mathbf{v})\in V^{**}\cong V$ is $S\mathbf{v}$). A $(2,0)$ tensor is an element of $V\otimes V$ (coefficient matrix $C^{ij}$). The space of $(p,q)$ tensors is $V^{\otimes p}\otimes(V^*)^{\otimes q}$ via (3.7) and Theorem 3.14: both have dimension $n^{p+q}$, and $\mathbf{e}_{i_1}\otimes\cdots\otimes\mathbf{e}^{j_q}$ corresponds to the multilinear map $(\varphi^1,\dots,\mathbf{v}_q)\mapsto\varphi^1(\mathbf{e}_{i_1})\cdots\mathbf{e}^{j_q}(\mathbf{v}_q)$. That is the first of the three reconciliations.

**Einstein summation.** An index appearing once up and once down in a product is summed over, and the summation sign is dropped: $\varphi(\mathbf{v}) = \varphi_iv^i$, $(A\mathbf{v})^i = A^i{}_jv^j$, $B(\mathbf{u},\mathbf{v}) = B_{ij}u^iv^j$, $(AB)^i{}_k = A^i{}_jB^j{}_k$. A repeated index is a **dummy** index: its letter is irrelevant, $A^i{}_jv^j = A^i{}_kv^k$. An unrepeated index is **free** and must match on both sides of an equation. Two rules that catch most errors: never repeat an index three times, and never have the same letter free on one side and dummy on the other.

**Contraction.** Setting an upper and a lower index equal and summing turns a $(p,q)$ tensor into a $(p-1,q-1)$ tensor. For a $(1,1)$ tensor, $T^i{}_i = \sum_iA_{ii} = \operatorname{tr}A$. For the $(2,1)$ tensor $A^i{}_jv^k$ (an operator tensored with a vector), contracting $j$ with $k$ gives $A^i{}_jv^j = (A\mathbf{v})^i$. The invariance of the trace under $A\mapsto P^{-1}AP$ is the statement that contraction produces a tensor (the $P$ and $P^{-1}$ cancel, §3.17). In `np.einsum`, contraction is the repeated letter: `einsum('ii', A)` is the trace, `einsum('ij,jk->ik', A, B)` is the matrix product, `einsum('ij,i,j', B, u, v)` is $B(\mathbf{u},\mathbf{v})$.

**Tensor product of tensors.** If $S$ is $(p,q)$ and $T$ is $(p',q')$, then $(S\otimes T)(\ldots) = S(\text{first }p+q\text{ arguments})\cdot T(\text{rest})$ is $(p+p', q+q')$ with components $S^{\cdots}{}_{\cdots}T^{\cdots}{}_{\cdots}$, products of components. Every tensor is a sum of such products of vectors and functionals; that is Theorem 3.14 again.

### 3.16 The metric tensor and raising/lowering indices

An inner product on a real $V$ is a symmetric positive-definite $(0,2)$ tensor, the **metric**, with components
$$g_{ij} = \langle\mathbf{e}_i,\mathbf{e}_j\rangle. \tag{3.15}$$
Riesz (Theorem 3.9) maps $\mathbf{v}\mapsto\langle\mathbf{v},\cdot\rangle$, a vector to a functional. In components: $\langle\mathbf{v},\mathbf{e}_j\rangle = \langle v^i\mathbf{e}_i,\mathbf{e}_j\rangle = v^ig_{ij}$, so the functional has components
$$v_j = g_{ij}v^i \qquad\text{(lowering an index)}. \tag{3.16}$$
The matrix $g$ is invertible (positive definite), and its inverse has components written $g^{ij}$, so $g^{ij}g_{jk} = \delta^i_k$ and $v^i = g^{ij}v_j$ (raising). Then $\langle\mathbf{u},\mathbf{v}\rangle = g_{ij}u^iv^j = u_jv^j = u^iv_i$: the inner product is a contraction once one index has been lowered.

In an **orthonormal** basis $g_{ij} = \delta_{ij}$, so $v_i = v^i$ numerically, and the whole up/down distinction is invisible. **This is why engineers never see it.** The stress tensor in an orthonormal frame has $\sigma_{ij} = \sigma^i{}_j = \sigma^{ij}$ as arrays, and ME texts write everything down and sum over anything repeated. The moment you use a non-orthonormal basis (crystallography, curvilinear coordinates in continuum mechanics, every chart in relativity), the metric comes back.

A non-orthonormal example where it matters, in $\mathbb{R}^2$ with the standard dot product. Basis $\mathbf{e}_1 = (1,0)^\top$, $\mathbf{e}_2 = (1,1)^\top$. Metric: $g_{11} = 1$, $g_{12} = g_{21} = 1$, $g_{22} = 2$, so $g = \begin{pmatrix}1&1\\1&2\end{pmatrix}$, $g^{-1} = \begin{pmatrix}2&-1\\-1&1\end{pmatrix}$. The vector $\mathbf{v} = (2,1)^\top$ (standard coordinates) is $\mathbf{v} = 1\cdot\mathbf{e}_1 + 1\cdot\mathbf{e}_2$, so $v^1 = v^2 = 1$. Lower: $v_1 = g_{11}v^1 + g_{21}v^2 = 2$, $v_2 = g_{12}v^1 + g_{22}v^2 = 3$. Check against the definition: $v_j = \langle\mathbf{v},\mathbf{e}_j\rangle$ gives $(2,1)\cdot(1,0) = 2$ and $(2,1)\cdot(1,1) = 3$. ✓ And $\|\mathbf{v}\|^2 = v^iv_i = 1\cdot2 + 1\cdot3 = 5 = 2^2 + 1^2$. ✓ The components $(1,1)$ and $(2,3)$ are different numbers describing the same arrow; neither is "the" components of $\mathbf{v}$ until you say up or down.

### 3.17 The transformation law, and the reconciliation theorem

Now the third definition. Change basis by $\mathbf{e}'_j = P^i{}_j\mathbf{e}_i$ (that is $P_{ij}$ with indices placed; columns of $P$ are the new basis vectors). The dual basis changes by the inverse: $\mathbf{e}'^i = (P^{-1})^i{}_k\mathbf{e}^k$. (Check: $\mathbf{e}'^i(\mathbf{e}'_j) = (P^{-1})^i{}_kP^l{}_j\,\mathbf{e}^k(\mathbf{e}_l) = (P^{-1})^i{}_kP^k{}_j = \delta^i_j$. ✓)

▶ **Try it first.** Put the primed basis vectors into (3.14) for a $(1,1)$ tensor and expand by multilinearity.

**Derivation of the transformation law.**

1. $(1,1)$ case: $T'^i{}_j = T(\mathbf{e}'^i,\mathbf{e}'_j) = T\big((P^{-1})^i{}_k\mathbf{e}^k, P^l{}_j\mathbf{e}_l\big) = (P^{-1})^i{}_k\,P^l{}_j\,T(\mathbf{e}^k,\mathbf{e}_l)$, so
$$T'^i{}_j = (P^{-1})^i{}_k\,T^k{}_l\,P^l{}_j, \qquad\text{i.e.}\qquad T' = P^{-1}TP. \tag{3.17}$$
2. $(0,2)$ case: $B'_{ij} = B(P^k{}_i\mathbf{e}_k, P^l{}_j\mathbf{e}_l) = P^k{}_iP^l{}_jB_{kl}$, i.e. $B' = P^\top BP$. Same as (3.5).
3. General case: each upper index contracts with a $P^{-1}$, each lower index with a $P$:
$$T'^{i_1\cdots i_p}{}_{j_1\cdots j_q} = (P^{-1})^{i_1}{}_{k_1}\cdots(P^{-1})^{i_p}{}_{k_p}\;P^{l_1}{}_{j_1}\cdots P^{l_q}{}_{j_q}\;T^{k_1\cdots k_p}{}_{l_1\cdots l_q}. \tag{3.18}$$
Each factor comes from expanding one argument by multilinearity. ∎

■ **What just happened.** (3.17) is $P^{-1}AP$ and step 2 is $P^\top BP$: the contrast from Theorem 3.11 is the special case $p + q = 2$ of (3.18). Upper index, $P^{-1}$; lower index, $P^\top$ (written with the indices in the order $P^l{}_j$, the transpose is automatic). Physics texts write $P^{-1}$ as $\partial x'^i/\partial x^k$ and $P$ as $\partial x^l/\partial x'^j$ because their bases are coordinate bases; the algebra is identical.

**Theorem 3.20 (Three definitions of a tensor agree).** The following are in natural bijection:
(i) multilinear maps $T: (V^*)^p\times V^q\to\mathbb{F}$;
(ii) elements of $V^{\otimes p}\otimes(V^*)^{\otimes q}$;
(iii) assignments, to every basis of $V$, of an array $T^{i_1\cdots i_p}{}_{j_1\cdots j_q}$ of $n^{p+q}$ numbers, such that the arrays for any two bases are related by (3.18).

*Proof.* (i)$\leftrightarrow$(ii) is §3.15: both spaces have dimension $n^{p+q}$ and the map sending $\mathbf{e}_{i_1}\otimes\cdots\otimes\mathbf{e}^{j_q}$ to the multilinear map $(\varphi^1,\dots,\mathbf{v}_q)\mapsto\varphi^1(\mathbf{e}_{i_1})\cdots\mathbf{e}^{j_q}(\mathbf{v}_q)$ sends a basis to a basis (the latter are the basis of (3.7)). (i)$\to$(iii): given $T$, define the array in each basis by (3.14); the derivation above shows the arrays satisfy (3.18). (iii)$\to$(i): given the assignment, pick one basis and define $T$ by multilinear extension of its array in that basis, $T(\varphi^1,\dots,\mathbf{v}_q) = \varphi^1_{i_1}\cdots v_q^{j_q}\,T^{i_1\cdots}{}_{\cdots j_q}$. This $T$ has, in *any other* basis, components given by (3.14), which by the derivation equal the (3.18)-transform of the original array, which by hypothesis is the array the assignment specifies for that basis. So $T$ reproduces the whole assignment, and does not depend on which basis was used to define it. The two constructions are inverse to each other. ∎

This is the theorem the roadmap promised. The third definition, "an array that transforms like a tensor," is not circular: it says *a tensor is a consistent family of arrays*, one per basis, and consistency is (3.18). You can now read Carroll's chapter 1 and Axler's chapter 9 and recognize the same object.

### 3.18 The physical tensors, fully worked

**Stress.** Cauchy's theorem: the traction on a face with unit outward normal $\mathbf{n}$ is $\mathbf{t} = \sigma\mathbf{n}$, $t^i = \sigma^i{}_jn^j$. Read this way, $\sigma$ is a $(1,1)$ tensor: eats a normal (vector), returns a force per area (vector). Balance of angular momentum on a small cube forces $\sigma_{ij} = \sigma_{ji}$. Alternatively, $\sigma(\mathbf{a},\mathbf{n}) = \mathbf{a}\cdot\sigma\mathbf{n}$ is the component along $\mathbf{a}$ of the traction on the $\mathbf{n}$-face: a $(0,2)$ tensor. In an orthonormal frame, these are the same array (§3.16). Under a rotation of axes $\mathbf{e}'_j = R_{ij}\mathbf{e}_i$ ($R$ orthogonal), both laws give
$$\sigma' = R^\top\sigma R, \tag{3.19}$$
because $R^{-1} = R^\top$. Your ME text writes $\sigma' = Q\sigma Q^\top$ with $Q = R^\top$ (rows of $Q$ are the new axes); it is the same equation. In 2D with $R = \begin{pmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{pmatrix}$, multiply out (3.19):
$$\sigma'_{x'x'} = \tfrac{\sigma_{xx}+\sigma_{yy}}2 + \tfrac{\sigma_{xx}-\sigma_{yy}}2\cos2\theta + \sigma_{xy}\sin2\theta, \qquad \sigma'_{x'y'} = -\tfrac{\sigma_{xx}-\sigma_{yy}}2\sin2\theta + \sigma_{xy}\cos2\theta. \tag{3.20}$$
(Derivation: $\sigma'_{x'x'} = \mathbf{e}'_1{}^\top\sigma\mathbf{e}'_1 = \sigma_{xx}c^2 + 2\sigma_{xy}sc + \sigma_{yy}s^2$ with $c = \cos\theta$, $s = \sin\theta$; use $c^2 = \tfrac{1+\cos2\theta}2$, $s^2 = \tfrac{1-\cos2\theta}2$, $2sc = \sin2\theta$. In the same way, $\sigma'_{x'y'} = \mathbf{e}'_1{}^\top\sigma\mathbf{e}'_2 = -\sigma_{xx}sc + \sigma_{xy}(c^2 - s^2) + \sigma_{yy}sc$.) Write $\tfrac{\sigma_{xx}-\sigma_{yy}}2 = r\cos2\theta_p$ and $\sigma_{xy} = r\sin2\theta_p$ with $r = \sqrt{\left(\tfrac{\sigma_{xx}-\sigma_{yy}}2\right)^2 + \sigma_{xy}^2}$. Then (3.20) becomes
$$\sigma'_{x'x'} = c_0 + r\cos(2\theta - 2\theta_p), \qquad \sigma'_{x'y'} = -r\sin(2\theta - 2\theta_p), \qquad c_0 = \tfrac{\sigma_{xx}+\sigma_{yy}}2. \tag{3.21}$$
As $\theta$ runs over $[0,\pi)$, the point $(\sigma'_{x'x'}, \sigma'_{x'y'})$ traces the circle of center $(c_0, 0)$ and radius $r$ once: **Mohr's circle.** The shear vanishes at $\theta = \theta_p$, where the normal stress is $c_0 + r$ (and at $\theta_p + 90°$, where it is $c_0 - r$). Those are the eigenvalues of $\sigma$: $\lambda = c_0\pm r$ is exactly the quadratic formula for $\det(\sigma - \lambda I) = 0$. **Mohr's circle is the 2D spectral theorem drawn as a picture**, with the orthogonal eigenvectors being the principal directions. The sign of the $\sin$ in (3.21) is why ME texts plot $\tau$ *downward*: with $\tau$ up, increasing $\theta$ moves the point clockwise. Example 4.5 has the numbers.

**Inertia.** For a body with density $\rho$,
$$I_{ij} = \int\rho\,(x_kx_k\,\delta_{ij} - x_ix_j)\,dV, \qquad I_{ij} = \sum_a m_a\big(|\mathbf{r}_a|^2\delta_{ij} - x_{a,i}x_{a,j}\big)\ \text{for point masses}. \tag{3.22}$$
$I$ is a symmetric $(0,2)$ tensor: $\mathbf{L} = I\boldsymbol\omega$ (or, as a form, $I(\boldsymbol\omega,\boldsymbol\omega) = 2\times$ kinetic energy). Derivation of (3.22) from $\mathbf{L} = \sum m_a\mathbf{r}_a\times(\boldsymbol\omega\times\mathbf{r}_a)$ and the identity $\mathbf{a}\times(\mathbf{b}\times\mathbf{c}) = \mathbf{b}(\mathbf{a}\cdot\mathbf{c}) - \mathbf{c}(\mathbf{a}\cdot\mathbf{b})$: $\mathbf{r}\times(\boldsymbol\omega\times\mathbf{r}) = \boldsymbol\omega|\mathbf{r}|^2 - \mathbf{r}(\mathbf{r}\cdot\boldsymbol\omega)$, whose $i$-component is $(|\mathbf{r}|^2\delta_{ij} - x_ix_j)\omega^j$. Because $I$ is real symmetric, the spectral theorem (LA-8) gives an orthonormal eigenbasis: the **principal axes**, in which $I$ is diagonal and $\mathbf{L}\parallel\boldsymbol\omega$. That is why principal axes exist for every rigid body, however irregular. The quadratic form $\boldsymbol\omega^\top I\boldsymbol\omega = \sum_am_a(|\mathbf{r}_a|^2|\boldsymbol\omega|^2 - (\mathbf{r}_a\cdot\boldsymbol\omega)^2) = \sum_am_a|\mathbf{r}_a\times\boldsymbol\omega|^2\ge 0$, so $I$ is positive semidefinite, with a zero eigenvalue only if all mass lies on one line. Example 4.6 computes one with `einsum`; CM-7 puts it into Euler's equations.

**Strain.** For a displacement field $\mathbf{u}(\mathbf{x})$, the small-strain tensor is $\varepsilon_{ij} = \tfrac12(\partial_iu_j + \partial_ju_i)$, the symmetric part of the displacement gradient $\partial_iu_j$ (the antisymmetric part is a local rotation, which produces no strain). It is a symmetric $(0,2)$ tensor, transforms by (3.19), has principal strains and a Mohr's circle, and Hooke's law for an isotropic solid is a linear map from one $(0,2)$ tensor to another, $\sigma_{ij} = \lambda\varepsilon_{kk}\delta_{ij} + 2\mu\varepsilon_{ij}$, which is itself a $(0,4)$ tensor $C_{ijkl}$ with the two Lamé constants as its only independent components in the isotropic case. **Trap:** the "engineering shear strain" $\gamma_{xy} = 2\varepsilon_{xy}$ is not a tensor component, and strain-transformation formulas written in terms of $\gamma$ carry factors of 2 that (3.20) does not.

### 3.19 Levi-Civita, cross products, and alternating forms (brief)

The **Levi-Civita symbol** $\varepsilon_{ijk}$ is $+1$ for $(ijk)$ an even permutation of $(123)$, $-1$ for odd, $0$ if any index repeats. The cross product is a contraction with it:
$$(\mathbf{a}\times\mathbf{b})_i = \varepsilon_{ijk}a_jb_k. \tag{3.23}$$
Check $i = 1$: $\varepsilon_{123}a_2b_3 + \varepsilon_{132}a_3b_2 = a_2b_3 - a_3b_2$. ✓ The identity $\varepsilon_{ijk}\varepsilon_{imn} = \delta_{jm}\delta_{kn} - \delta_{jn}\delta_{km}$ (exercise B6) turns every vector identity of EM-1 into two lines of index algebra. In $\mathbb{R}^3$ with orthonormal bases, $\varepsilon_{ijk}$ is the component array of the $(0,3)$ tensor $(\mathbf{a},\mathbf{b},\mathbf{c})\mapsto\det(\mathbf{a}\ \mathbf{b}\ \mathbf{c})$, the volume form; it is invariant under rotations ($\det R = 1$) and flips sign under reflections, which is why the cross product is a "pseudovector."

A $(0,k)$ tensor is **alternating** if it changes sign when any two arguments are swapped. Alternating $k$-forms on an $n$-dimensional space form a vector space of dimension $\binom nk$ (a basis is $\mathbf{e}^{i_1}\wedge\cdots\wedge\mathbf{e}^{i_k}$ for $i_1 < \cdots < i_k$, where the **wedge** $\alpha\wedge\beta = \alpha\otimes\beta - \beta\otimes\alpha$ for 1-forms), so for $k = n$ there is exactly one up to scale. LA-6 defined the determinant as the alternating $n$-linear function of the columns with $\det I = 1$; in this language, $\det$ is the unique top-degree alternating form normalized on the standard basis, and $\det(A) = $ the factor by which $A$ scales it. Axler ch. 9 does this in full; it is the cleanest proof that $\det(AB) = \det A\det B$. The electromagnetic field tensor $F_{\mu\nu}$ of EM-8 is an alternating $(0,2)$ tensor on 4-dimensional spacetime; its six independent components are $\mathbf{E}$ and $\mathbf{B}$. That is the only pointer to it in this module.

### 3.20 Outer products in data (brief)

The coefficient matrix of $\mathbf{x}\otimes\mathbf{y}$ is the outer product $\mathbf{x}\mathbf{y}^\top$, and the SVD $A = \sum_k\sigma_k\mathbf{u}_k\mathbf{v}_k^\top$ (LA-9) is the Schmidt decomposition of $A$ viewed as an element of $\mathbb{R}^m\otimes\mathbb{R}^n$. A **rank-1 update** $A + \mathbf{x}\mathbf{y}^\top$ is a tensor-product perturbation, and the Sherman–Morrison formula for its inverse (LA-2, exercise C4) is cheap for the same reason Theorem 3.15 is true. The **Khatri–Rao product** (column-wise Kronecker) appears in tensor decompositions of data arrays with three or more indices; it is a name to recognize, not something this module develops.

## 4. Worked examples

### Example 4.1: Dual basis of $\{(1,1), (1,-1)\}$, and why it is not the transpose

Basis of $\mathbb{R}^2$: $\mathbf{e}_1 = (1,1)^\top$, $\mathbf{e}_2 = (1,-1)^\top$. Basis matrix $E = \begin{pmatrix}1&1\\1&-1\end{pmatrix}$, $\det E = -2$.

By (3.3) the dual basis vectors are the rows of $E^{-1} = \tfrac{1}{-2}\begin{pmatrix}-1&-1\\-1&1\end{pmatrix} = \begin{pmatrix}\tfrac12&\tfrac12\\ \tfrac12&-\tfrac12\end{pmatrix}$:
$$\mathbf{e}^1 = (\tfrac12, \tfrac12), \qquad \mathbf{e}^2 = (\tfrac12, -\tfrac12).$$
Check (3.2): $\mathbf{e}^1(\mathbf{e}_1) = \tfrac12 + \tfrac12 = 1$, $\mathbf{e}^1(\mathbf{e}_2) = \tfrac12 - \tfrac12 = 0$, $\mathbf{e}^2(\mathbf{e}_1) = 0$, $\mathbf{e}^2(\mathbf{e}_2) = \tfrac12 + \tfrac12 = 1$. ✓

The rows of $E^\top$ are $(1,1)$ and $(1,-1)$: the right directions (this basis is orthogonal) but twice too long, because the basis vectors have length $\sqrt2$, not 1. Using them as coordinate functions on $\mathbf{v} = (3,1)^\top$ gives $(4, 2)$, but $\mathbf{v} = 2\mathbf{e}_1 + 1\mathbf{e}_2$, and the true dual basis gives $\mathbf{e}^1(\mathbf{v}) = 2$, $\mathbf{e}^2(\mathbf{v}) = 1$. ✓

Now a non-orthogonal basis: $\mathbf{e}_1 = (1,0)^\top$, $\mathbf{e}_2 = (1,1)^\top$, $E = \begin{pmatrix}1&1\\0&1\end{pmatrix}$, $E^{-1} = \begin{pmatrix}1&-1\\0&1\end{pmatrix}$. Dual basis: $\mathbf{e}^1 = (1,-1)$, $\mathbf{e}^2 = (0,1)$. The rows of $E^\top$ would be $(1,0)$ and $(1,1)$: wrong directions entirely. $\mathbf{e}^1$ is perpendicular to $\mathbf{e}_2$, not to $\mathbf{e}_1$: the $i$-th dual vector is orthogonal to all the *other* basis vectors, which is the geometric content of $\mathbf{e}^i(\mathbf{e}_j) = 0$ for $j\ne i$. In the orthonormal case "orthogonal to all the others" happens to mean "parallel to $\mathbf{e}_i$," and the distinction disappears.

### Example 4.2: Same numbers, same $P$, two different answers

Let $M = \begin{pmatrix}2&1\\1&3\end{pmatrix}$ and change basis by the shear $P = \begin{pmatrix}1&1\\0&1\end{pmatrix}$ (new basis $\mathbf{e}'_1 = \mathbf{e}_1$, $\mathbf{e}'_2 = \mathbf{e}_1 + \mathbf{e}_2$). $P^{-1} = \begin{pmatrix}1&-1\\0&1\end{pmatrix}$.

Common factor: $MP = \begin{pmatrix}2&1\\1&3\end{pmatrix}\begin{pmatrix}1&1\\0&1\end{pmatrix} = \begin{pmatrix}2&3\\1&4\end{pmatrix}$.

*As an operator* (Theorem 3.11, $A' = P^{-1}MP$):
$$P^{-1}MP = \begin{pmatrix}1&-1\\0&1\end{pmatrix}\begin{pmatrix}2&3\\1&4\end{pmatrix} = \begin{pmatrix}1&-1\\1&4\end{pmatrix}.$$
Trace $5 = \operatorname{tr}M$, determinant $4 + 1 = 5 = \det M$, eigenvalues $\tfrac{5\pm\sqrt5}2 \approx 1.382, 3.618$, the same as $M$'s. Not symmetric: similarity does not preserve symmetry unless $P$ is orthogonal.

*As a bilinear form* ($B' = P^\top MP$):
$$P^\top MP = \begin{pmatrix}1&0\\1&1\end{pmatrix}\begin{pmatrix}2&3\\1&4\end{pmatrix} = \begin{pmatrix}2&3\\3&7\end{pmatrix}.$$
Trace $9\ne 5$. Determinant $14 - 9 = 5 = (\det P)^2\det M$ with $\det P = 1$. Eigenvalues $\tfrac{9\pm\sqrt{61}}2\approx 0.595, 8.405$: different from $M$'s, but both still positive. Signature $(2,0,0)$ before and after: Sylvester (Theorem 3.12) in action. Still symmetric: congruence preserves symmetry.

Two matrices, $\begin{pmatrix}1&-1\\1&4\end{pmatrix}$ and $\begin{pmatrix}2&3\\3&7\end{pmatrix}$, both "the matrix of $M$ in the new basis." The question "which is right?" has no answer until you say what $M$ *is*. The project reproduces both.

### Example 4.3: Schmidt decomposition of $\tfrac12(|00\rangle + |01\rangle + |10\rangle - |11\rangle)$

Coefficient matrix (3.8), row index = qubit 0, column index = qubit 1:
$$C = \tfrac12\begin{pmatrix}1&1\\1&-1\end{pmatrix}.$$
$\det C = \tfrac14(-1-1) = -\tfrac12\ne 0$, so rank 2: entangled (Theorem 3.15). How entangled? $C^\dagger C = \tfrac14\begin{pmatrix}1&1\\1&-1\end{pmatrix}^2 = \tfrac14\begin{pmatrix}2&0\\0&2\end{pmatrix} = \tfrac12 I$, so both singular values are $\sigma_1 = \sigma_2 = 1/\sqrt2$. Schmidt coefficients $(1/\sqrt2, 1/\sqrt2)$, $\sigma_1^2 + \sigma_2^2 = 1$ ✓, entanglement entropy $-2\cdot\tfrac12\log\tfrac12 = \log 2$: **maximally entangled**, as entangled as the Bell state.

Explicit Schmidt vectors: since $C = \tfrac1{\sqrt2}H$ with $H = \tfrac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}$ orthogonal, one SVD is $U = I$, $\Sigma = \tfrac1{\sqrt2}I$, $V^\top = H$. By Theorem 3.16, $\mathbf{a}_k = $ columns of $U$: $|0\rangle, |1\rangle$; $\mathbf{b}_k = $ columns of $V$ = rows of $H$: $|+\rangle = \tfrac1{\sqrt2}(|0\rangle + |1\rangle)$ and $|-\rangle = \tfrac1{\sqrt2}(|0\rangle - |1\rangle)$. So
$$\psi = \tfrac1{\sqrt2}\big(|0\rangle\otimes|+\rangle + |1\rangle\otimes|-\rangle\big).$$
Check by expanding: $\tfrac1{\sqrt2}\cdot\tfrac1{\sqrt2}(|00\rangle + |01\rangle + |10\rangle - |11\rangle) = \tfrac12(\ldots)$. ✓ (Because the singular values are equal, the Schmidt basis is not unique: `np.linalg.svd` returns a different but equally valid pair, which the project prints.)

### Example 4.4: A Kronecker product with its eigenvalues checked

$A = \begin{pmatrix}2&1\\0&3\end{pmatrix}$ (eigenvalues $2, 3$), $B = \begin{pmatrix}1&2\\2&1\end{pmatrix}$ (eigenvalues $3, -1$). By (3.12),
$$A\otimes B = \begin{pmatrix}2B & 1B\\ 0B & 3B\end{pmatrix} = \begin{pmatrix}2&4&1&2\\4&2&2&1\\0&0&3&6\\0&0&6&3\end{pmatrix}.$$
Theorem 3.17(e): eigenvalues $\{2\cdot3, 2\cdot(-1), 3\cdot3, 3\cdot(-1)\} = \{6, -2, 9, -3\}$. Check (d): $\operatorname{tr} = 2 + 2 + 3 + 3 = 10 = 5\cdot 2$ ✓, and $6 - 2 + 9 - 3 = 10$ ✓. Check (f): $\det = (\det A)^2(\det B)^2 = 36\cdot 9 = 324$, and $6\cdot(-2)\cdot9\cdot(-3) = 324$ ✓. An eigenvector for $6$: $\mathbf{x} = (1,0)^\top$ for $A$'s eigenvalue 2, $\mathbf{y} = (1,1)^\top$ for $B$'s eigenvalue 3, so $\mathbf{x}\otimes\mathbf{y} = (1,1,0,0)^\top$; and $(A\otimes B)(1,1,0,0)^\top = (6,6,0,0)^\top$ ✓.

### Example 4.5: Stress at 30°, by $P^\top\sigma P$ and by Mohr's circle

$\sigma = \begin{pmatrix}80&30\\30&20\end{pmatrix}$ MPa, rotate axes by $\theta = 30°$: $c = \cos30° = 0.8660$, $s = \sin30° = 0.5$, $P = \begin{pmatrix}c&-s\\s&c\end{pmatrix}$.

*Matrix route (3.19).* $\sigma P = \begin{pmatrix}80c + 30s & -80s + 30c\\ 30c + 20s & -30s + 20c\end{pmatrix} = \begin{pmatrix}84.282 & -14.019\\ 35.981 & 2.321\end{pmatrix}$. Then $P^\top\sigma P$: entry $(1,1) = c(84.282) + s(35.981) = 72.990 + 17.990 = 90.981$; entry $(1,2) = c(-14.019) + s(2.321) = -12.141 + 1.160 = -10.981$; entry $(2,2) = -s(-14.019) + c(2.321) = 7.010 + 2.010 = 9.019$. So
$$\sigma' = \begin{pmatrix}90.98 & -10.98\\ -10.98 & 9.02\end{pmatrix}\text{ MPa}.$$
Trace check: $90.98 + 9.02 = 100 = 80 + 20$ ✓ (rotations are similarities too).

*Mohr route (3.21).* Center $c_0 = \tfrac{80+20}2 = 50$. Radius $r = \sqrt{30^2 + 30^2} = 30\sqrt2 = 42.43$. Principal angle: $\tan2\theta_p = \sigma_{xy}/\tfrac{\sigma_{xx}-\sigma_{yy}}2 = 30/30 = 1$, so $2\theta_p = 45°$, $\theta_p = 22.5°$. Principal stresses $50\pm42.43 = 92.43, 7.57$ MPa (the eigenvalues of $\sigma$: $\lambda^2 - 100\lambda + 700 = 0$ gives $50\pm\sqrt{1800}$ ✓). At $\theta = 30°$, $2\theta - 2\theta_p = 60° - 45° = 15°$:
$$\sigma'_{x'x'} = 50 + 42.43\cos15° = 50 + 40.98 = 90.98, \qquad \sigma'_{x'y'} = -42.43\sin15° = -10.98.$$
Agreement to the last digit. The project draws the circle with the $\theta = 0$ point $(80, 30)$, the $\theta = 30°$ point $(90.98, -10.98)$, and the principal points on the axis.

### Example 4.6: Inertia tensor of three point masses, in `einsum`, with principal axes

Masses $m = (1, 2, 1)$ at $\mathbf{r}_1 = (1,1,0)$, $\mathbf{r}_2 = (-1,0,1)$, $\mathbf{r}_3 = (0,1,1)$ (about the origin, not the center of mass, which is at $(-1,2,3)/4$). Each $|\mathbf{r}_a|^2 = 2$. From (3.22):

$a=1$ ($m=1$): $2I - \mathbf{r}\mathbf{r}^\top = \begin{pmatrix}2&0&0\\0&2&0\\0&0&2\end{pmatrix} - \begin{pmatrix}1&1&0\\1&1&0\\0&0&0\end{pmatrix} = \begin{pmatrix}1&-1&0\\-1&1&0\\0&0&2\end{pmatrix}$.

$a=2$ ($m=2$): $2\left(2I - \begin{pmatrix}1&0&-1\\0&0&0\\-1&0&1\end{pmatrix}\right) = \begin{pmatrix}2&0&2\\0&4&0\\2&0&2\end{pmatrix}$.

$a=3$ ($m=1$): $2I - \begin{pmatrix}0&0&0\\0&1&1\\0&1&1\end{pmatrix} = \begin{pmatrix}2&0&0\\0&1&-1\\0&-1&1\end{pmatrix}$.

Sum:
$$I = \begin{pmatrix}5&-1&2\\-1&6&-1\\2&-1&5\end{pmatrix}.$$
In `einsum`: `I = einsum('a,a,ij->ij', m, r2, eye(3)) - einsum('a,ai,aj->ij', m, X, X)` with `r2 = einsum('ai,ai->a', X, X)`. The string `'a,ai,aj->ij'` *is* $\sum_am_ax_{a,i}x_{a,j}$.

Principal moments: $\operatorname{tr}I = 16$, $\det I = 5(30 - 1) + 1(-5 + 2) + 2(1 - 12) = 145 - 3 - 22 = 120$. Try $\mathbf{v} = (1,0,-1)^\top$: $I\mathbf{v} = (5-2, -1+1, 2-5)^\top = 3\mathbf{v}$, eigenvalue $3$. The other two have sum $13$ and product $40$: $5$ and $8$. Eigenvectors: for $8$, $(I - 8I_3)\mathbf{v} = 0$ is solved by $(1,-1,1)^\top$ (check: $I(1,-1,1)^\top = (5+1+2, -1-6-1, 2+1+5)^\top = (8,-8,8)^\top$ ✓); for $5$, orthogonal to both: $(1,2,1)^\top$ (check: $I(1,2,1)^\top = (5-2+2, -1+12-1, 2-2+5)^\top = (5,10,5)^\top$ ✓). Principal axes $\tfrac1{\sqrt2}(1,0,-1)$, $\tfrac1{\sqrt6}(1,2,1)$, $\tfrac1{\sqrt3}(1,-1,1)$ with moments $3, 5, 8$. Spinning about the axes of smallest and largest moment, $(1,0,-1)$ and $(1,-1,1)$, is stable; spinning about the middle axis $(1,2,1)$ is not. That is the tennis-racket theorem, and CM-7 derives it from Euler's equations.

### Example 4.7: Rodrigues' formula for a 90° rotation about $(1,1,1)/\sqrt3$

From Appendix A, (A.4): $R = I + \sin\theta\,K + (1 - \cos\theta)K^2$ where $K = [\mathbf{n}]_\times$ is the antisymmetric matrix with $K\mathbf{v} = \mathbf{n}\times\mathbf{v}$. With $\mathbf{n} = \tfrac1{\sqrt3}(1,1,1)$ and $\theta = 90°$: $\sin\theta = 1$, $1 - \cos\theta = 1$, so $R = I + K + K^2$.
$$K = \tfrac1{\sqrt3}\begin{pmatrix}0&-1&1\\1&0&-1\\-1&1&0\end{pmatrix}, \qquad K^2 = \mathbf{n}\mathbf{n}^\top - I = \tfrac13\begin{pmatrix}1&1&1\\1&1&1\\1&1&1\end{pmatrix} - I = \begin{pmatrix}-\tfrac23&\tfrac13&\tfrac13\\ \tfrac13&-\tfrac23&\tfrac13\\ \tfrac13&\tfrac13&-\tfrac23\end{pmatrix}.$$
(The identity $K^2 = \mathbf{n}\mathbf{n}^\top - I$ for a unit $\mathbf{n}$ is $\mathbf{n}\times(\mathbf{n}\times\mathbf{v}) = \mathbf{n}(\mathbf{n}\cdot\mathbf{v}) - \mathbf{v}$.) Adding,
$$R = \begin{pmatrix}\tfrac13 & \tfrac13 - \tfrac1{\sqrt3} & \tfrac13 + \tfrac1{\sqrt3}\\ \tfrac13 + \tfrac1{\sqrt3} & \tfrac13 & \tfrac13 - \tfrac1{\sqrt3}\\ \tfrac13 - \tfrac1{\sqrt3} & \tfrac13 + \tfrac1{\sqrt3} & \tfrac13\end{pmatrix} \approx \begin{pmatrix}0.3333 & -0.2440 & 0.9107\\ 0.9107 & 0.3333 & -0.2440\\ -0.2440 & 0.9107 & 0.3333\end{pmatrix}.$$
Checks: $R\mathbf{n} = \mathbf{n}$ (each row sums to $1$, since $K\mathbf{n} = 0$ and $K^2\mathbf{n} = 0$) ✓; $\operatorname{tr}R = 1 = 1 + 2\cos90°$ ✓; $\det R = 1$ ✓. $R\mathbf{e}_1 = (0.333, 0.911, -0.244)^\top$ makes angle $\arccos(0.333) = 70.5°$ with $\mathbf{e}_1$, not $90°$, which is right: $\mathbf{e}_1$ is not perpendicular to the axis, so only its perpendicular component turns through the full $90°$. Since four quarter-turns are a full turn, $R^4 = I$; and $\theta = 120°$ about the same axis gives the cyclic permutation $\mathbf{e}_1\to\mathbf{e}_2\to\mathbf{e}_3$, which you can confirm from (A.4) with $\sin\theta = \tfrac{\sqrt3}2$, $1 - \cos\theta = \tfrac32$. The project compares this $R$ to `scipy.spatial.transform.Rotation.from_rotvec`.

## 5. Why it matters and where it connects

- **CM-7.** The inertia tensor (3.22) is the central object. Euler's equations are $\mathbf{L} = I\boldsymbol\omega$ differentiated in a rotating frame; principal axes exist by LA-8's spectral theorem applied to the symmetric $(0,2)$ tensor $I$; the rotation group of Appendix A is the configuration space of a rigid body, and Euler angles are coordinates on it.
- **EM-8.** The field tensor $F_{\mu\nu}$ is an alternating $(0,2)$ tensor on spacetime (§3.19); Lorentz transformations are the $P$ in (3.18) that preserve the Minkowski metric, a symmetric $(0,2)$ tensor of signature $(1,3,0)$ by Theorem 3.12; raising and lowering indices with $\eta_{\mu\nu}$ is (3.16) with a non-positive-definite $g$.
- **QM-4.** Bras are functionals (§3.6), $|\psi\rangle\langle\phi|$ is a simple tensor in $V\otimes V^*$ (§3.12), the adjoint is the dual map composed with Riesz (Theorem 3.5), and the completeness relation $\sum_i|i\rangle\langle i| = I$ is Theorem 3.2's formula $\mathbf{v} = \sum_i\mathbf{e}^i(\mathbf{v})\mathbf{e}_i$.
- **QM-7.** The state space of two particles is $V\otimes W$ (Theorem 3.14); entanglement is Schmidt rank $> 1$ (Theorems 3.15–3.16); the reduced density matrix of qubit 0 is $CC^\dagger$ and its eigenvalues are the $\sigma_k^2$.
- **QC-2.** Two-qubit gates are $4\times 4$ matrices; $U_1\otimes U_2$ is the Kronecker product (3.12) in the qubit-0-leftmost convention; CNOT is *not* a Kronecker product, which is why it creates entanglement. The density matrices of **QC-7** make the entanglement entropy $-\sum\sigma_k^2\log\sigma_k^2$ of §3.13 the von Neumann entropy of $CC^\dagger$.
- **QM-5, QC-1.** Appendix A: angular momentum operators are the generators of $SO(3)$; spin-$\tfrac12$ is the $SU(2)$ representation; single-qubit gates are $SU(2)$ elements and the Bloch sphere is $SO(3)$ acting on it through the 2-to-1 map.
- **LA-14.** The vec trick (3.13) is how the derivative of a matrix-valued function of a matrix is written as a matrix; the Hessian of a function of a matrix is a $(0,4)$ tensor flattened by $\operatorname{vec}$.
- **ML-9.** "Tensor" in PyTorch means "multi-index array" with no transformation law; the shapes bookkeeping of attention, `einsum('bhqd,bhkd->bhqk', Q, K)`, is §3.15's contraction with batch and head indices along for the ride.
- **LA-9, LA-10.** SVD is the Schmidt decomposition of a matrix; rank-1 updates are tensor-product perturbations.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 drill dual bases, the transformation laws, and small Kronecker products. B1–B7 are the module: B2 (prove $(ST)^* = T^*S^*$ and read off $(BA)^\top$), B3 (the Schmidt decomposition of a given state by hand), B4 (Sylvester on a $3\times 3$ form), B5 (the transformation law for a $(1,2)$ tensor), B6 (the $\varepsilon$–$\delta$ identity and two vector identities), and B7 (extend `project.py`: solve a Sylvester equation via the vec trick) are essential. C1–C3 go further: the universal property used to prove $V\otimes W\cong W\otimes V$ and $(V\otimes W)^*\cong V^*\otimes W^*$, the three-qubit state that has no Schmidt decomposition, and the Lie algebra isomorphism $\mathfrak{so}(3)\cong\mathfrak{su}(2)$ made explicit.

## 7. Coding project

`project.py` makes every theorem in §3 numerical:

1. **Dual bases.** `dual_basis(E)` returns the rows of $E^{-1}$ (3.3) and asserts $\mathbf{e}^i(\mathbf{e}_j) = \delta^i_j$. Example 4.1's two bases are run, and the dual basis is shown to differ from the transpose.
2. **Change of basis, two ways.** `change_basis_operator(A, P)` and `change_basis_form(B, P)` implement (3.5). Example 4.2's numbers are reproduced. On random symmetric forms and random invertible $P$, the signature is asserted invariant (Sylvester) while the eigenvalues are not; on random operators the eigenvalues are asserted invariant.
3. **Kronecker identities.** Theorem 3.17(a)–(f) checked on random matrices with `np.kron`; Example 4.4 reproduced; the vec trick (3.13) checked in both column-major and row-major form.
4. **Schmidt decomposition.** `schmidt_decomposition(psi, dims)` reshapes to the coefficient matrix and takes the SVD (Theorem 3.16). Run on the four Bell states, on product states, on random two-qubit states, on Example 4.3's state, and on a random $3\times 4$ bipartite state; the state is reconstructed from (3.9). `entanglement_entropy` asserts $0$ for products and $\log 2$ for Bell states.
5. **`einsum` and index notation.** Trace, matrix product, bilinear form, the inertia tensor of Example 4.6 with its principal axes, the cross product as $\varepsilon_{ijk}a_jb_k$ (3.23), the $\varepsilon$–$\delta$ identity, the metric example of §3.16, and the transformation law (3.18) applied to a $(1,1)$ and a $(0,2)$ tensor and compared with $P^{-1}AP$ and $P^\top BP$.
6. **Stress and Mohr's circle.** Example 4.5 by matrix and by (3.21), asserted equal; the circle plotted with the $0°$, $30°$, and principal points marked.
7. **Rotations.** `rodrigues(n, theta)` implements (A.4), compared with `scipy.spatial.transform.Rotation`; Example 4.7 reproduced. The tangent space of $SO(3)$ is shown to be antisymmetric by numerical differentiation; the commutators $[L_i, L_j] = \varepsilon_{ijk}L_k$ are asserted; $SU(2)\to SO(3)$ is checked via $U(\mathbf{n}'\cdot\boldsymbol\sigma)U^\dagger = (R\mathbf{n}')\cdot\boldsymbol\sigma$ for random axes and angles; and $\theta = 2\pi$ gives $U = -I$, $R = I$.

**Build it yourself first:** `dual_basis`, `change_basis_form`, `schmidt_decomposition` (three lines: reshape, SVD, return), and `rodrigues`. Write them from the equations before reading the code.

**Your turn** (in the script): compute the metric tensor of a skewed basis in $\mathbb{R}^3$ and lower an index; find the principal moments of a thin rectangular plate by sampling points and compare with $\tfrac{1}{12}m(a^2 + b^2)$; solve $AX + XB = C$ via (3.13) and `np.linalg.solve`.

## 8. Where you'll struggle

- **You will feel you are learning the same thing three times, and you are.** Multilinear map, element of a tensor product, array with a transformation law. The roadmap warned you. The fix is to pick one definition as home (I recommend the multilinear map, because components are then *literally* "plug in basis vectors," (3.14)) and to translate the other two into it every time, until the translation is reflexive. Theorem 3.20 is the dictionary. Do not try to hold all three as primary.

- **The direction of $P$.** Is it $P^\top BP$ or $PBP^\top$? $P^{-1}AP$ or $PAP^{-1}$? Both appear in books, depending on whether $P$'s columns are the new basis in old coordinates (ours) or the old in new. This is where index gymnastics gets its bad reputation, and it is not your fault. The fix: never memorize the formula; re-derive it from "$\mathbf{e}'_j = \sum_iP_{ij}\mathbf{e}_i$, plug into $B(\mathbf{e}'_j,\mathbf{e}'_l)$" in two lines, every time, for the first month. Theorem 3.11's derivation is the template. And when in doubt, test with $P = \operatorname{diag}(2, 1)$ on a $2\times 2$ example: vector components should halve.

- **Covariant vs contravariant feels arbitrary until you stretch a basis.** For rotations the two are identical, so every example from ME gives you no feedback. Do §3.16's skewed example and the metric "Your turn" by hand until the two sets of components $(1,1)$ and $(2,3)$ feel like two descriptions of one arrow. Then do the "rulers and arrows" picture of §2 once more.

- **Qubit ordering will burn you in QC.** `np.kron(A, B)` puts $A$ on qubit 0, the *left* factor, the *most significant* bit of the basis index: $|01\rangle$ is index $1$, $|10\rangle$ is index $2$. Qiskit uses the *opposite* (qubit 0 is the rightmost factor and least significant bit). Both are self-consistent; mixing them produces CNOTs that act on the wrong qubit and no error message. The fix: write a one-line `ket(bits)` helper that builds basis states by Kronecker product in your convention and *only ever* construct states through it. The project includes one.

- **The universal property feels like nothing until you need it.** On first reading, §3.11's paragraph seems to say "a tensor product is a thing with a property." You will get it the first time you try to define $A\otimes B$ on a sum like $\mathbf{v}_1\otimes\mathbf{w}_1 + \mathbf{v}_2\otimes\mathbf{w}_2$ and worry whether a different way of writing the same tensor gives a different answer. The universal property is the one-line guarantee that it does not. Exercise C1 makes you use it twice; after that it stops feeling empty.

- **The Schmidt decomposition proof is five lines and you will not believe it.** "Reshape and take the SVD" cannot possibly be the structure theorem for entanglement. It is. The content is in Theorem 3.14, which says the reshape is legitimate (the $\mathbf{e}_i\otimes\mathbf{f}_j$ are a basis), and in LA-9, which did the hard work. If it feels too easy, do C2: for three particles there is *no* Schmidt decomposition, and seeing why that fails tells you what the two-particle case was really using.

- **Mohr's circle sign conventions.** (3.21) has $\sigma'_{x'y'} = -r\sin(2\theta - 2\theta_p)$, so with shear plotted *upward* the point moves *clockwise* as $\theta$ increases, and the angle on the circle is $2\theta$, not $\theta$. Different texts flip the shear axis, or the sign of $\sigma_{xy}$ on one face, to make it counterclockwise. Do not memorize a direction; compute $P^\top\sigma P$ once and see which way your book's circle goes.

- **The appendix is dense.** It compresses a chapter of Stillwell into 1000 words. Read it once now for the shape (rotations form a group; their derivatives at the identity are antisymmetric; the exponential gets you back; $SU(2)$ does the same with $2\times 2$ complex matrices and covers $SO(3)$ twice). Skim the rest. Come back and read every line when QM-5 or QC-1 needs it; it will go much faster then.

- **Exercise B5 is harder than it looks.** Writing the transformation law for a $(1,2)$ tensor is one line of (3.18); *implementing* it in `einsum` and checking it against a composition of matrix products takes most people an hour of index-shuffling. Budget for it. The error is almost always one transposed $P$.

## 9. Check yourself

1. What is the dual basis to $\mathbf{e}_1 = (2,0)^\top$, $\mathbf{e}_2 = (0,1)^\top$? Is it the transpose of the basis matrix?
2. State the defining property of the dual map $T^*$ and say in one sentence why its matrix is $A^\top$.
3. Why is $V\cong V^{**}$ "natural" but $V\cong V^*$ not?
4. Riesz: given $\varphi(\mathbf{u}) = 3u^1 - u^2$ on $\mathbb{C}^2$ with the standard inner product, what vector $\mathbf{v}$ has $\varphi = \langle\mathbf{v},\cdot\rangle$? (Careful: our convention conjugates the first slot.)
5. Under $\mathbf{e}'_j = \sum_iP_{ij}\mathbf{e}_i$, how do the matrices of an operator and a bilinear form change? In one sentence, *why* do they differ?
6. What is preserved by congruence $B\mapsto P^\top BP$ that is not preserved by similarity, and what is preserved by similarity that is not preserved by congruence?
7. What is $\dim(V\otimes W)$, and what is the dimension of the subset of simple tensors (as a surface, not a subspace) when $\dim V = \dim W = 2$?
8. A two-qubit state has coefficient matrix $\begin{pmatrix}a&b\\c&d\end{pmatrix}$. Give the one-number test for entanglement.
9. Write $(A\otimes B)_{(i,k),(j,l)}$ and say which qubit $A$ acts on under `np.kron(A, B)`.
10. What does $\operatorname{vec}(AXB)$ equal? Column-major or row-major?
11. Write the transformation law for a $(2,1)$ tensor in index notation.
12. In an orthonormal basis of $\mathbb{R}^3$, why can you ignore the difference between $\sigma^i{}_j$ and $\sigma_{ij}$? Give a basis where you cannot.
13. What is the tangent space of $SO(3)$ at the identity, and what is $\exp$ of a small element of it?

## 10. Pair with

- **Axler, *Linear Algebra Done Right*, 4th ed., ch. 9 (Multilinear Algebra and Determinants), and §3F (Duality).** Axler's duality section is the cleanest proof that the dual map's matrix is the transpose, and his ch. 9 treats bilinear forms, quadratic forms, alternating forms, and the tensor product in exactly this module's order with full proofs. Read §3F alongside §3.1–3.5 and ch. 9 alongside §3.7–3.13. Open access.
- **Halmos, *Finite-Dimensional Vector Spaces*, the sections Dual Spaces, Brackets, Dual Bases, Reflexivity, Annihilators (§§13–17 in the 2nd ed.) and Bilinear Forms, Tensor Products, Product Bases (§§23–25).** Older notation, but Halmos is more careful than anyone about what is natural and what requires a choice. His Reflexivity and Tensor Products sections do what §3.3 and §3.11 do here, more slowly.
- **Carroll, *Spacetime and Geometry*, ch. 1 (§1.4–1.7, vectors through manipulating tensors, especially).** The physicist's definition of a tensor by its transformation law, with the metric raising and lowering indices from the start. Read it after §3.17 to see definition (iii) spoken natively, and to see the $(1,3)$ signature at work.
- **Schutz, *A First Course in General Relativity*, ch. 3 (Tensor analysis in special relativity).** Gentler than Carroll on one-forms as "machines that eat vectors," with the clearest pictures of covariant versus contravariant you will find. If §3.16 did not click, read Schutz §3.3–3.5.
- **Stillwell, *Naive Lie Theory*, ch. 1–2.** Rotations of the plane and space, quaternions, $SU(2)$ and the 2-to-1 map, at undergraduate pace with no manifolds. Appendix A is a compression of these two chapters; Stillwell has the pictures and the quaternion viewpoint that the appendix omits.
- **Nielsen & Chuang, §2.1.7 (tensor products) and §2.5 (Schmidt decomposition).** Two pages each, in qubit language, with the ordering convention this module uses. Read after §3.13.

---

## Appendix A: $SO(3)$ and $SU(2)$

*Read once now for shape; read closely when QM-5 or QC-1 sends you back.*

**A.1 The rotation group.** A rotation of $\mathbb{R}^3$ is a linear map preserving lengths and orientation: $R^\top R = I$ and $\det R = +1$. The set of such matrices is $SO(3)$ ("special orthogonal"). It is a **group**: products and inverses of rotations are rotations ($(R_1R_2)^\top R_1R_2 = R_2^\top R_1^\top R_1R_2 = I$, $\det(R_1R_2) = 1$; $R^{-1} = R^\top$ is also in it). It is not commutative (LA-1 Example 4.3). It is a 3-dimensional surface inside the 9-dimensional space of matrices: $R^\top R = I$ imposes 6 independent equations (symmetric matrix), $9 - 6 = 3$. Those 3 parameters are an axis (2) and an angle (1), or three Euler angles (CM-7).

**A.2 The tangent space at the identity.** Let $R(t)$ be a smooth curve in $SO(3)$ with $R(0) = I$. Differentiate $R(t)^\top R(t) = I$ at $t = 0$:
$$\dot R(0)^\top R(0) + R(0)^\top\dot R(0) = 0 \quad\Longrightarrow\quad \dot R(0)^\top + \dot R(0) = 0. \tag{A.1}$$
Every velocity vector of a curve through $I$ is an **antisymmetric** matrix. Conversely every antisymmetric $K$ is such a velocity: $R(t) = \exp(tK)$ satisfies $R^\top R = \exp(tK^\top)\exp(tK) = \exp(-tK)\exp(tK) = I$ (the two exponents commute, so the LA-7 rule $e^Xe^Y = e^{X+Y}$ applies) and $\det\exp(tK) = e^{t\operatorname{tr}K} = e^0 = 1$. So the tangent space is exactly the 3-dimensional space $\mathfrak{so}(3)$ of antisymmetric $3\times 3$ matrices, with basis
$$L_x = \begin{pmatrix}0&0&0\\0&0&-1\\0&1&0\end{pmatrix}, \quad L_y = \begin{pmatrix}0&0&1\\0&0&0\\-1&0&0\end{pmatrix}, \quad L_z = \begin{pmatrix}0&-1&0\\1&0&0\\0&0&0\end{pmatrix}, \qquad (L_i)_{jk} = -\varepsilon_{ijk}. \tag{A.2}$$
$L_z$ is the derivative at $\theta = 0$ of the rotation by $\theta$ about $z$; the others likewise. For any vector $\boldsymbol\omega$, $(\omega_xL_x + \omega_yL_y + \omega_zL_z)\mathbf{v} = \boldsymbol\omega\times\mathbf{v}$ (check with (3.23): $(\omega_iL_i)_{jk}v_k = -\varepsilon_{ijk}\omega_iv_k = \varepsilon_{jik}\omega_iv_k = (\boldsymbol\omega\times\mathbf{v})_j$). So the antisymmetric matrix $K = [\boldsymbol\omega]_\times$ *is* "cross product with $\boldsymbol\omega$," which is why angular velocity is a vector in 3D (and only in 3D: $\dim\mathfrak{so}(n) = n(n-1)/2$ equals $n$ only for $n = 3$).

**A.3 Commutators.** Direct multiplication gives $L_xL_y - L_yL_x = L_z$ and cyclically:
$$[L_i, L_j] = \varepsilon_{ijk}L_k. \tag{A.3}$$
(Check one: $L_xL_y = \begin{pmatrix}0&0&0\\1&0&0\\0&0&0\end{pmatrix}$, $L_yL_x = \begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}$, difference $= L_z$. ✓ The project checks all nine.) The commutator measures the failure of rotations to commute to second order: $e^{sL_x}e^{tL_y}e^{-sL_x}e^{-tL_y} = I + st[L_x, L_y] + O(3)$. In QM-5 the angular momentum operators satisfy $[J_i, J_j] = i\hbar\varepsilon_{ijk}J_k$, which is (A.3) with $J = i\hbar L$; the whole theory of angular momentum is the study of matrices satisfying (A.3).

**A.4 The exponential map and Rodrigues' formula.** For a unit axis $\mathbf{n}$ and $K = [\mathbf{n}]_\times$, the rotation by $\theta$ about $\mathbf{n}$ is $R = \exp(\theta K)$. To evaluate the series, note $K^2 = \mathbf{n}\mathbf{n}^\top - I$ (the identity $\mathbf{n}\times(\mathbf{n}\times\mathbf{v}) = \mathbf{n}(\mathbf{n}\cdot\mathbf{v}) - \mathbf{v}$), and then $K^3 = K(\mathbf{n}\mathbf{n}^\top - I) = (K\mathbf{n})\mathbf{n}^\top - K = -K$ since $K\mathbf{n} = \mathbf{n}\times\mathbf{n} = 0$. So the powers cycle: $K^3 = -K$, $K^4 = -K^2$, $K^5 = K$, and
$$\exp(\theta K) = I + \theta K + \tfrac{\theta^2}2K^2 - \tfrac{\theta^3}{6}K - \tfrac{\theta^4}{24}K^2 + \cdots = I + \Big(\theta - \tfrac{\theta^3}6 + \cdots\Big)K + \Big(\tfrac{\theta^2}2 - \tfrac{\theta^4}{24} + \cdots\Big)K^2,$$
$$R(\mathbf{n},\theta) = I + \sin\theta\,K + (1 - \cos\theta)\,K^2. \tag{A.4}$$
Applied to a vector: $R\mathbf{v} = \mathbf{v} + \sin\theta\,\mathbf{n}\times\mathbf{v} + (1-\cos\theta)(\mathbf{n}(\mathbf{n}\cdot\mathbf{v}) - \mathbf{v}) = \cos\theta\,\mathbf{v} + \sin\theta\,\mathbf{n}\times\mathbf{v} + (1 - \cos\theta)(\mathbf{n}\cdot\mathbf{v})\mathbf{n}$: the component along $\mathbf{n}$ is fixed and the perpendicular component turns by $\theta$. Example 4.7 evaluates it. Every element of $SO(3)$ is $\exp$ of some $\theta K$, because every rotation has an axis. The eigenvalues of $R$ have modulus $1$ ($R$ preserves length), the non-real ones come in conjugate pairs (LA-7 Theorem 3.16), a real cubic has a real root, and the product of all three eigenvalues is $\det R = 1$. So the real eigenvalues are $\pm1$ and $+1$ occurs: with one real eigenvalue it equals $\det R = 1$; with three, an even number are $-1$. A unit eigenvector $\mathbf{n}$ for $+1$ is the axis, $R$ maps the plane $\mathbf{n}^\perp$ to itself as a 2D rotation by some $\theta$, and that is exactly what (A.4) describes. So $\exp: \mathfrak{so}(3)\to SO(3)$ is onto, but not one-to-one: $\theta$ and $\theta + 2\pi$ give the same $R$.

**A.5 $SU(2)$.** The group of $2\times 2$ complex matrices with $U^\dagger U = I$ and $\det U = 1$. Writing $U = \begin{pmatrix}a&b\\-\bar b&\bar a\end{pmatrix}$ with $|a|^2 + |b|^2 = 1$ (exercise: this is the general form), $SU(2)$ is the unit sphere $S^3$ in $\mathbb{C}^2 = \mathbb{R}^4$: 3-dimensional, like $SO(3)$. Differentiating $U^\dagger U = I$ and $\det U = 1$ along a curve through $I$: the tangent space $\mathfrak{su}(2)$ consists of **anti-Hermitian, traceless** matrices ($\dot U^\dagger = -\dot U$, $\operatorname{tr}\dot U = 0$). A basis is $-\tfrac i2\sigma_j$ with the Pauli matrices
$$\sigma_x = \begin{pmatrix}0&1\\1&0\end{pmatrix}, \quad \sigma_y = \begin{pmatrix}0&-i\\i&0\end{pmatrix}, \quad \sigma_z = \begin{pmatrix}1&0\\0&-1\end{pmatrix}, \qquad \sigma_j\sigma_k = \delta_{jk}I + i\varepsilon_{jkl}\sigma_l. \tag{A.5}$$
From (A.5), $[\sigma_j,\sigma_k] = 2i\varepsilon_{jkl}\sigma_l$, so
$$\big[-\tfrac i2\sigma_j, -\tfrac i2\sigma_k\big] = -\tfrac14\cdot 2i\varepsilon_{jkl}\sigma_l = \varepsilon_{jkl}\big(-\tfrac i2\sigma_l\big). \tag{A.6}$$
**The same commutation relations as (A.3).** Two different groups, two different sizes of matrix, the same infinitesimal structure: $\mathfrak{so}(3)\cong\mathfrak{su}(2)$ as Lie algebras (exercise C3). The exponential is again explicit: for a unit $\mathbf{n}$, $(\mathbf{n}\cdot\boldsymbol\sigma)^2 = I$ by (A.5), so
$$U(\mathbf{n},\theta) = \exp\!\big(-\tfrac{i\theta}2\,\mathbf{n}\cdot\boldsymbol\sigma\big) = \cos\tfrac\theta2\,I - i\sin\tfrac\theta2\,\mathbf{n}\cdot\boldsymbol\sigma. \tag{A.7}$$
Note the $\theta/2$: this is the spin-$\tfrac12$ rotation operator of QM-5 and the general single-qubit gate of QC-1 (up to a global phase).

**A.6 The 2-to-1 map $SU(2)\to SO(3)$.** A real vector $\mathbf{v}\in\mathbb{R}^3$ becomes a traceless Hermitian matrix $\mathbf{v}\cdot\boldsymbol\sigma = \begin{pmatrix}v_z & v_x - iv_y\\ v_x + iv_y & -v_z\end{pmatrix}$, with $\det(\mathbf{v}\cdot\boldsymbol\sigma) = -|\mathbf{v}|^2$. For $U\in SU(2)$, the matrix $U(\mathbf{v}\cdot\boldsymbol\sigma)U^\dagger$ is again traceless Hermitian, so equals $\mathbf{v}'\cdot\boldsymbol\sigma$ for a unique real $\mathbf{v}'$, and $\det$ is preserved, so $|\mathbf{v}'| = |\mathbf{v}|$. The map $\mathbf{v}\mapsto\mathbf{v}'$ is linear and length-preserving: an orthogonal matrix $R(U)$, and by continuity from $U = I$ it has $\det = +1$. So $U\mapsto R(U)$ is a map $SU(2)\to SO(3)$, and it is a homomorphism: $R(U_1U_2) = R(U_1)R(U_2)$ since conjugating twice is conjugating by the product. **Claim:** $R(U(\mathbf{n},\theta))$ is the rotation by $\theta$ about $\mathbf{n}$. Demonstration for $\mathbf{n} = \mathbf{e}_z$: $U = \operatorname{diag}(e^{-i\theta/2}, e^{i\theta/2})$, and $U\sigma_xU^\dagger = \begin{pmatrix}0 & e^{-i\theta}\\ e^{i\theta} & 0\end{pmatrix} = \cos\theta\,\sigma_x + \sin\theta\,\sigma_y$, which is $\mathbf{e}_x\mapsto\cos\theta\,\mathbf{e}_x + \sin\theta\,\mathbf{e}_y$: rotation by $\theta$ (not $\theta/2$) about $z$. ✓ The general axis follows by conjugating with a $U$ that moves $\mathbf{e}_z$ to $\mathbf{n}$; the project verifies $U(\mathbf{n}'\cdot\boldsymbol\sigma)U^\dagger = (R\mathbf{n}')\cdot\boldsymbol\sigma$ numerically for random $\mathbf{n}, \theta, \mathbf{n}'$ with $R$ from (A.4). **2-to-1:** $U$ and $-U$ give the same conjugation, hence the same $R$. Concretely, $\theta = 2\pi$ in (A.7) gives $U = \cos\pi\,I = -I$, while (A.4) gives $R = I$. A spin-$\tfrac12$ particle rotated by $360°$ comes back to minus itself; it takes $720°$ to return. This is not a metaphor: it is (A.7) at $\theta = 2\pi$, it has been measured with neutron interferometers, and the relative phase it produces is observable.

**A.7 Why this matters.** QM-5: angular momentum is the representation theory of (A.3), and spin-$\tfrac12$ exists because $SU(2)$ has a 2-dimensional representation that $SO(3)$ lacks. QC-1: every single-qubit gate is (A.7) times a phase, and the Bloch sphere is $\mathbb{R}^3$ with $SU(2)$ acting through $R(U)$. CM-7: $\boldsymbol\omega = [\dot R R^\top]_\times^{-1}$ is the angular velocity, and (A.1) is why it exists. Stillwell ch. 1–2 does all of this with quaternions, which are $SU(2)$ in disguise.
