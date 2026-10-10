# LA-4 Exercises

Conventions as in the notes: bases are ordered lists; $[\mathbf{v}]_{\mathcal{B}}$ is the coordinate vector; column $j$ of $[T]_{\mathcal{B},\mathcal{C}}$ is $[T\mathbf{b}_j]_{\mathcal{C}}$; the change-of-basis matrix $P$ has the new basis vectors as columns in old coordinates, so $P[\mathbf{v}]_{\text{new}} = [\mathbf{v}]_{\text{old}}$ and $[T]_{\text{new}} = P^{-1}[T]_{\text{old}}P$. $\mathcal{P}_n$ uses the monomial basis $(1, x, \dots, x^n)$ unless stated. $\mathbb{R}^{2\times 2}$ uses the basis $(E_{11}, E_{12}, E_{21}, E_{22})$, i.e. a matrix is flattened row by row.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Linear or not? For each map, decide and justify in one line (either verify (3.1) or exhibit a specific failure).
(a) $T(x, y) = (x + 2y,\; 3x - y,\; y)$ from $\mathbb{R}^2$ to $\mathbb{R}^3$.
(b) $T(x, y) = (x + 1,\; y)$ on $\mathbb{R}^2$.
(c) $T(p) = p(2)$ from $\mathcal{P}_3$ to $\mathbb{R}$.
(d) $T(p) = x\,p'(x)$ on $\mathcal{P}_3$.
(e) $T(A) = A^\top A$ on $\mathbb{R}^{2\times 2}$.
(f) $T(A) = (\operatorname{tr}A)\,I$ on $\mathbb{R}^{2\times 2}$.

**A2.** For $T$ in A1(a), write $[T]$ with respect to the standard bases (what shape is it?). Compute $T(1, -1)$ directly and via (3.3). Find $\ker T$ and $\operatorname{im}T$, state rank and nullity, and say whether $T$ is injective, surjective, both, or neither.

**A3.** Let $\mathcal{B}' = \big((2, 1), (1, 1)\big)$ in $\mathbb{R}^2$.
(a) Write $P$ and $P^{-1}$.
(b) Find $[\,(3, 2)\,]_{\mathcal{B}'}$ twice: by solving $a(2,1) + b(1,1) = (3,2)$, and by one matrix–vector product. Say which matrix you used and why.
(c) Which vector of $\mathbb{R}^2$ has $\mathcal{B}'$-coordinates $(3, 2)$?

**A4.** Let $A = \begin{pmatrix} 1 & 2 & 3 \\ 2 & 4 & 6\end{pmatrix}$. Find a basis of $\ker A$ and a basis of $\operatorname{im}A$, check Rank–Nullity, and say whether $\mathbf{x} \mapsto A\mathbf{x}$ is injective or surjective.

**A5.** (a) Multiply out $R_\theta R_\varphi$ by row-times-column. Explain why Theorem 3.7 forces $R_\theta R_\varphi = R_{\theta+\varphi}$, and read off the angle-addition formulas for $\cos(\theta + \varphi)$ and $\sin(\theta + \varphi)$ from the comparison.
(b) Compute $F_\varphi^2$ for the reflection matrix of §3.13 and confirm it is $I$. What does $F_\varphi R_\theta$ do to the plane, and what is its determinant?

**A6.** (a) Write $[D]$ for $D: \mathcal{P}_2 \to \mathcal{P}_2$ in the monomial basis, then compute $[D]^2$ and $[D]^3$. Which coordinate vectors does $[D]^2$ send to zero, and what polynomials are those?
(b) Write the matrix of $T(p) = x\,p'(x)$ on $\mathcal{P}_2$ in the monomial basis. Why is it diagonal? (One sentence; the word you want arrives in LA-7.)

## B. Core (do all)

**B1.** (Isomorphisms preserve everything.) Let $T: V \to W$ be an isomorphism and $\mathbf{v}_1, \dots, \mathbf{v}_k \in V$.
(a) Prove: $(\mathbf{v}_i)$ is linearly independent in $V$ if and only if $(T\mathbf{v}_i)$ is linearly independent in $W$.
(b) Prove: $(\mathbf{v}_i)$ spans $V$ if and only if $(T\mathbf{v}_i)$ spans $W$.
(c) Apply (a) to $\kappa_{\mathcal{B}}$ to justify the rule "to test whether polynomials are independent, test their coefficient vectors." Then decide whether $1 + x,\; 1 - x,\; 1 + x^2$ are independent in $\mathcal{P}_2$.

**B2.** (Integration and differentiation.) Let $J: \mathcal{P}_2 \to \mathcal{P}_3$, $(Jp)(x) = \int_0^x p(t)\,dt$, and $D: \mathcal{P}_3 \to \mathcal{P}_2$, $Dp = p'$, both in monomial bases.
(a) Write $[J]$ ($4\times 3$) and $[D]$ ($3\times 4$).
(b) Compute $[D][J]$ and $[J][D]$. Interpret both products as statements about calculus. Which polynomials does $J\circ D$ fix, and which does it kill?
(c) Find $\ker J$, $\operatorname{im}J$, $\ker D$, $\operatorname{im}D$, and verify Rank–Nullity for each map.
(d) Explain using Rank–Nullity why no linear map $\mathcal{P}_3 \to \mathcal{P}_2$ can have a two-sided inverse, and why $J$ has a left inverse but cannot have a right inverse.

**B3.** (Change of basis, both directions.) Let $T$ have standard matrix $A = \begin{pmatrix} 2 & 1 \\ 1 & 2\end{pmatrix}$ and let $\mathcal{B}' = \big((2, 1), (1, 1)\big)$ as in A3.
(a) Compute $[T]_{\mathcal{B}'} = P^{-1}AP$, every intermediate matrix shown.
(b) Verify the answer from Definition 3.5: compute $T\mathbf{b}'_1$ and $T\mathbf{b}'_2$ and write each in $\mathcal{B}'$-coordinates.
(c) Compute the wrong-direction product $PAP^{-1}$. Compare its trace and determinant with those of (a). What does this tell you about using invariants to catch a direction error? What basis is $PAP^{-1}$ actually the matrix of $T$ in?
(d) Find a basis in which the matrix of $T$ is diagonal, and compute that diagonal matrix with (3.6). (Hint: try $(1, 1)$ and $(1, -1)$ and look at what $T$ does to each.)

**B4.** (Trace.)
(a) Prove $\operatorname{tr}(XYZ) = \operatorname{tr}(ZXY)$ for square $X, Y, Z$, and give $2\times 2$ matrices with $\operatorname{tr}(XYZ) \neq \operatorname{tr}(YXZ)$.
(b) Show that $I$ and $N = \begin{pmatrix} 1 & 1 \\ 0 & 1\end{pmatrix}$ have the same trace and determinant but are not similar.
(c) Prove there are no $n\times n$ real matrices with $XY - YX = I$. (In QM-4 the position and momentum operators satisfy $[\hat x, \hat p] = i\hbar\,I$; this exercise is why they cannot be finite matrices.)
(d) Define the trace of an operator. Compute $\operatorname{tr}D$ for $D$ on $\mathcal{P}_3$, the trace of the projection of Example 4.3, and the trace of the transpose map on $\mathbb{R}^{2\times 2}$.

**B5.** (One-sided inverses.)
(a) Let $A, B$ be $n\times n$ with $AB = I$. Prove $BA = I$. Mark exactly where Theorem 3.13, Corollary 3.15, and associativity are used.
(b) Give a $2\times 3$ matrix $A$ and a $3\times 2$ matrix $B$ with $AB = I_2$ but $BA \neq I_3$. Use Rank–Nullity (or the rank of a product, C1) to prove $BA = I_3$ is impossible for *any* such $A, B$.

**B6.** (Isomorphisms.)
(a) Prove that the composition of two isomorphisms is an isomorphism and that the inverse of an isomorphism is an isomorphism. (So "$\cong$" is an equivalence relation.)
(b) Write down an explicit isomorphism $\Phi: \mathcal{P}_3 \to \mathbb{R}^{2\times 2}$ and its inverse. Transport $D$ across it: compute $(\Phi\circ D\circ\Phi^{-1})\begin{pmatrix} a & b \\ c & d\end{pmatrix}$.
(c) Prove: a linear map $T: V \to W$ is an isomorphism if and only if it sends some basis of $V$ to a basis of $W$.

**B7.** (The transpose map, by hand.) On $\mathbb{R}^{2\times 2}$ with basis $(E_{11}, E_{12}, E_{21}, E_{22})$, let $T(A) = A^\top$.
(a) Write the $4\times 4$ matrix $[T]$. Show $[T]^2 = I$. What are $\ker T$ and $\operatorname{im}T$?
(b) Let $K(A) = A - A^\top$. Show $[K] = I - [T]$, find a basis of $\ker K$ and of $\operatorname{im}K$ as matrices (not only as coordinate vectors), and check $\dim\ker K + \dim\operatorname{im}K = 4$.
(c) Let $S = \tfrac12(I + T)$, i.e. $S(A) = \tfrac12(A + A^\top)$. Show $S^2 = S$, identify its image and kernel, and check that $\operatorname{tr}[S]$ equals its rank. (C3 explains the last coincidence.)

**B8.** (Code.) Extend `project.py`:
(a) Write `is_linear(T, basis, trials=20)` that draws random scalars $a, b$ and random combinations $\mathbf{u}, \mathbf{v}$ of the basis vectors and checks $T(a\mathbf{u} + b\mathbf{v}) = aT\mathbf{u} + bT\mathbf{v}$ numerically. Run it on `D`, on the transpose map, and on the coefficient-wise square $p \mapsto p^2$. Explain why `matrix_of_map` returns a matrix for the non-linear map anyway, and what that matrix predicts wrongly.
(b) Write `kernel_basis(T, basis_in, basis_out)` and `image_basis(T, basis_in, basis_out)` that return bases of $\ker T$ and $\operatorname{im}T$ as *ambient* vectors (polynomial coefficient vectors, flattened matrices), by computing `matrix_of_map`, taking `scipy.linalg.null_space` / `scipy.linalg.orth`, and converting coordinate vectors back to vectors. Apply both to $D$ on $\mathcal{P}_3$ and to $A \mapsto A - A^\top$. Assert that $T$ sends every returned kernel vector to zero and that the dimensions add to $\dim V$.
(c) In a comment, explain why converting back uses `basis_in` for the kernel and `basis_out` for the image.

## C. Challenge (do at least one)

**C1.** (Rank of a composition.) Let $T: V \to W$ and $S: W \to U$ be linear.
(a) Prove $\operatorname{rank}(S\circ T) \le \min(\operatorname{rank}S, \operatorname{rank}T)$.
(b) Prove $\dim\ker(S\circ T) \le \dim\ker T + \dim\ker S$. (Hint: apply Rank–Nullity to the restriction of $S$ to the subspace $\operatorname{im}T$.)
(c) Deduce Sylvester's inequality: $\operatorname{rank}(S\circ T) \ge \operatorname{rank}S + \operatorname{rank}T - \dim W$. Give matrices where it is an equality and where it is strict.

**C2.** (Nilpotent index.) Let $N$ be a nilpotent operator on $V$, $\dim V = n$, with index $k$ (so $N^k = 0$ and $N^{k-1} \neq 0$).
(a) Show $\ker N \subseteq \ker N^2 \subseteq \ker N^3 \subseteq \cdots$.
(b) Show that if $\ker N^j = \ker N^{j+1}$ for some $j$, then $\ker N^{j+1} = \ker N^{j+2}$, and hence the chain never grows again.
(c) Conclude that the chain is strictly increasing up to $\ker N^k = V$, and therefore $k \le n$.
(d) Show $D$ on $\mathcal{P}_n$ has index exactly $n + 1 = \dim\mathcal{P}_n$, and give a nilpotent operator on $\mathbb{R}^3$ of index $2$.
(e) Prove $I - N$ is invertible, with inverse $I + N + N^2 + \cdots + N^{k-1}$.

**C3.** (Structure of projections.) Let $P: V \to V$ satisfy $P^2 = P$.
(a) Prove every $\mathbf{v} \in V$ is *uniquely* $\mathbf{k} + \mathbf{w}$ with $\mathbf{k} \in \ker P$, $\mathbf{w} \in \operatorname{im}P$. (This is what $V = \ker P \oplus \operatorname{im}P$ means.)
(b) Prove that in a basis consisting of a basis of $\operatorname{im}P$ followed by a basis of $\ker P$, $[P] = \operatorname{diag}(1, \dots, 1, 0, \dots, 0)$ with $\operatorname{rank}P$ ones. Conclude $\operatorname{tr}P = \operatorname{rank}P$.
(c) Conversely, given subspaces $U, K$ with $V = U \oplus K$, prove there is exactly one projection with image $U$ and kernel $K$.
(d) Show $I - P$ is a projection and identify its image and kernel.
(e) Compute the standard matrix of the projection onto the line $y = 2x$ *along* the direction $(1, 0)$ (image $\operatorname{span}(1, 2)$, kernel $\operatorname{span}(1, 0)$). Verify $P^2 = P$, note $P \neq P^\top$, and compare with the orthogonal projection of Example 4.3. Which of the two does LA-5 single out, and by what extra condition?
