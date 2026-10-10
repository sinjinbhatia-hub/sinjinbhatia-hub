# LA-3 Exercises

Conventions as in the notes: $A$ is $m \times n$, $R$ its RREF with $r$ pivots, $EA = R$; the four subspaces are $C(A), N(A), C(A^\top), N(A^\top)$; a basis is an ordered list. "Prove" means a proof in the style of §3.0: full sentences, every step a hypothesis, a definition, an earlier result, or a logical step.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** For each set, decide whether it is a subspace of the stated space. If yes, say which description (span or null space) you'd use. If no, give an explicit witness: a vector in the set whose scalar multiple leaves it, or two vectors in the set whose sum leaves it, or the observation that $\mathbf{0}$ is missing.
(a) $\{(x, y) : x \ge 0\}$ in $\mathbb{R}^2$.
(b) $\{(x, y) : y = 2x + 1\}$ in $\mathbb{R}^2$.
(c) $\{(x, y) : xy = 0\}$ in $\mathbb{R}^2$.
(d) $\{(x, y, z) : x + 2y - z = 0\}$ in $\mathbb{R}^3$.
(e) $\{p \in \mathcal{P}_3 : p(0) = p(1)\}$ in $\mathcal{P}_3$.
(f) $\{A \in \mathcal{M}_{2\times 2} : \det A = 0\}$ in $\mathcal{M}_{2\times 2}$.

**A2.** Let $A = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & -1 \end{pmatrix}$ (already in RREF) and $\mathbf{b} = (4, 1)^\top$.
(a) Identify the pivot and free columns and $r$.
(b) Write the special solutions and a basis of $N(A)$.
(c) Find a particular solution $\mathbf{x}_p$ with free variables zero and write the complete solution (3.3).
(d) State the dimensions of all four fundamental subspaces of $A$.

**A3.** Let $\mathbf{v}_1 = (1, 2, 3)^\top$, $\mathbf{v}_2 = (2, 3, 4)^\top$, $\mathbf{v}_3 = (3, 4, 5)^\top$.
(a) Use the pivot test (Theorem 3.12) to decide whether they are independent; if not, exhibit the relation.
(b) Describe $\operatorname{span}\{\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3\}$ geometrically and give it as a null space $\{\mathbf{x} : \mathbf{n}\cdot\mathbf{x} = 0\}$.
(c) Is $(0, 1, 2)^\top$ in the span? If so, write it as a combination of $\mathbf{v}_1, \mathbf{v}_2$.

**A4.** Let $\mathbf{u}_1 = (1, 0, 1)^\top$, $\mathbf{u}_2 = (0, 1, 1)^\top$, $\mathbf{u}_3 = (1, 1, 1)^\top$.
(a) Is $\mathbf{u}_3 \in \operatorname{span}\{\mathbf{u}_1, \mathbf{u}_2\}$? (An existence question.)
(b) Is $(\mathbf{u}_1, \mathbf{u}_2, \mathbf{u}_3)$ a basis of $\mathbb{R}^3$? Answer using Lemma 3.20 and Corollary 3.26, without row reducing.
(c) Find the coordinates of $(2, 3, 4)^\top$ in this basis. (A uniqueness question: why is there exactly one answer?)

**A5.** Let $A = \begin{pmatrix} 1 & 1 & 2 \\ 2 & 1 & 3 \\ 3 & 0 & 3 \end{pmatrix}$. Reduce $[A \mid I]$ to $[R \mid E]$ by hand (no pivoting; first nonzero entry is the pivot). Give bases and dimensions for $C(A)$, $N(A)$, $C(A^\top)$, $N(A^\top)$, check the two orthogonality relations (3.6) on your bases, and draw the big picture with the numbers filled in.

## B. Core (do all; these are the module)

**B1.** Let $U$ and $W$ be subspaces of $V$.
(a) Prove that $U + W = \{\mathbf{u} + \mathbf{w} : \mathbf{u} \in U, \mathbf{w} \in W\}$ is a subspace of $V$, and that it is the smallest subspace containing $U \cup W$.
(b) Prove that $U \cup W$ is a subspace if and only if $U \subseteq W$ or $W \subseteq U$.

**B2.** From the eight axioms and Theorem 3.2 only, prove:
(a) $c\mathbf{0} = \mathbf{0}$ for every scalar $c$.
(b) $(-1)\mathbf{v} = -\mathbf{v}$. (First show $(-1)\mathbf{v}$ is *a* negative of $\mathbf{v}$; then show negatives are unique.)
(c) The zero vector is unique.
(d) If $c\mathbf{v} = \mathbf{0}$ then $c = 0$ or $\mathbf{v} = \mathbf{0}$.

**B3.** (Symmetric and skew-symmetric matrices; mechanics.) Let $\mathrm{Sym}_n = \{A \in \mathcal{M}_{n\times n} : A^\top = A\}$ and $\mathrm{Skew}_n = \{A : A^\top = -A\}$.
(a) Prove $\mathrm{Sym}_n$ is a subspace of $\mathcal{M}_{n\times n}$.
(b) Find a basis of $\mathrm{Sym}_2$ and prove it is one (use Theorem 3.17). State $\dim \mathrm{Sym}_2$.
(c) Prove $\dim \mathrm{Sym}_n = n(n+1)/2$ and $\dim \mathrm{Skew}_n = n(n-1)/2$.
(d) Prove every $A \in \mathcal{M}_{n\times n}$ is $S + K$ with $S \in \mathrm{Sym}_n$, $K \in \mathrm{Skew}_n$ in exactly one way. Check the dimensions add up. (In CM-7 the inertia tensor lives in $\mathrm{Sym}_3$; in continuum mechanics the velocity gradient splits as strain rate plus spin exactly this way.)

**B4.** (a) Prove $W = \{p \in \mathcal{P}_3 : p(1) = 0\}$ is a subspace of $\mathcal{P}_3$. Using coefficient coordinates $(a_0, a_1, a_2, a_3)$, write $W$ as the null space of a $1 \times 4$ matrix, find its special solutions, translate them back into polynomials, and state $\dim W$.
(b) Show $\{p \in \mathcal{P}_3 : p(1) = 1\}$ is not a subspace.
(c) Let $W_2 = \{p \in \mathcal{P}_3 : p(1) = 0 \text{ and } p(2) = 0\}$. Find a basis and $\dim W_2$ by the same method, and factor each basis polynomial.

**B5.** (The subspaces of $\mathbb{R}^3$.) Prove that every subspace $W$ of $\mathbb{R}^3$ is one of: $\{\mathbf{0}\}$, a line through $\mathbf{0}$, a plane through $\mathbf{0}$, or $\mathbb{R}^3$. Then prove that every plane through $\mathbf{0}$ has both descriptions of §2.1: it is $\operatorname{span}\{\mathbf{d}_1, \mathbf{d}_2\}$ for some independent $\mathbf{d}_1, \mathbf{d}_2$ and it is $N(\mathbf{n}^\top)$ for some nonzero $1 \times 3$ matrix $\mathbf{n}^\top$. (Hint for the second: what is $\dim N(\mathbf{n}^\top)$? Then use Corollary 3.26(d).)

**B6.** (Strang's $A = CR$.) Let $A$ be $m \times n$ of rank $r$, let $C$ be the $m \times r$ matrix of pivot columns of $A$, and let $\tilde R$ be the $r \times n$ matrix of nonzero rows of $R$.
(a) Prove $A = C\tilde R$. (Use Lemma 3.27 and the shape of $R$.)
(b) Verify it numerically for the matrix of Example 4.2.
(c) Use $A = C\tilde R$ to give a second proof that $\dim C(A^\top) = \dim C(A)$. (Show every row of $A$ is a combination of the $r$ rows of $\tilde R$, conclude $\dim C(A^\top) \le r$, and then apply the same inequality to $A^\top$.)

**B7.** (The matrix $A^\top A$; needed in LA-5 and ST-3.)
(a) Prove $N(A^\top A) = N(A)$. (For $\subseteq$: if $A^\top A\mathbf{x} = \mathbf{0}$, compute $\mathbf{x}^\top A^\top A\mathbf{x}$.)
(b) Deduce $\operatorname{rank}(A^\top A) = \operatorname{rank}(A)$.
(c) Deduce that the $n \times n$ matrix $A^\top A$ is invertible if and only if the columns of $A$ are independent.
(d) Prove $C(A^\top A) = C(A^\top)$.

**B8.** (Code.) In `project.py`, Part 5 implements `in_span` and `extend_to_basis`. Add two functions:
`basis_of_sum(U, W)` returning a basis of $U + W$ given bases $U$ ($n \times p$) and $W$ ($n \times q$) as columns; and `basis_of_intersection(U, W)` returning a basis of $U \cap W$. For the intersection, explain why $\mathbf{x} \in U \cap W$ iff $\mathbf{x} = U\mathbf{a} = W\mathbf{b}$ for some $\mathbf{a}, \mathbf{b}$, iff $(\mathbf{a}, \mathbf{b})$ is in the null space of the $n \times (p + q)$ matrix $[\,U \mid -W\,]$, and use `null_space_basis` on it. Test both on random subspaces of $\mathbb{R}^6$ (construct $U$ and $W$ sharing some columns so the intersection is nontrivial) and assert the dimension formula of C1, $\dim(U + W) = \dim U + \dim W - \dim(U \cap W)$.

## C. Challenge (do at least one)

**C1.** (Dimension formula.) Let $U, W$ be subspaces of a finite-dimensional $V$. Prove
$$\dim(U + W) = \dim U + \dim W - \dim(U \cap W).$$
Strategy: take a basis $\mathbf{z}_1, \dots, \mathbf{z}_k$ of $U \cap W$; extend it to a basis of $U$ and, separately, to a basis of $W$; show the union of the three lists is a basis of $U + W$. The independence step needs an argument in two stages. Then apply the formula to show two distinct planes through $\mathbf{0}$ in $\mathbb{R}^3$ meet in a line.

**C2.** (Exponentials are independent; LA-7 and F2.) Let $\lambda_1, \dots, \lambda_k$ be distinct real numbers. Prove that $e^{\lambda_1 t}, \dots, e^{\lambda_k t}$ are linearly independent in $C^\infty(\mathbb{R})$. (Induction on $k$: given a relation, divide by $e^{\lambda_k t}$, differentiate, and apply the induction hypothesis.) Then explain in one paragraph why the same proof works for distinct *complex* $\lambda_j$ if scalars are complex, and deduce that $1, \cos t, \sin t, \cos 2t, \sin 2t$ are independent over $\mathbb{R}$.

**C3.** (Rank inequalities and numerical rank.)
(a) Prove $\operatorname{rank}(A + B) \le \operatorname{rank} A + \operatorname{rank} B$.
(b) Prove that if $P$ ($m \times m$) and $Q$ ($n \times n$) are invertible, then $\operatorname{rank}(PAQ) = \operatorname{rank} A$. (Two applications of (3.7), one of them to $A = P^{-1}(PAQ)Q^{-1}$.)
(c) (Code.) For $A_\varepsilon = \begin{pmatrix} 1 & 1 \\ 1 & 1 + \varepsilon \end{pmatrix}$ with $\varepsilon = 10^{-6}, 10^{-12}, 10^{-16}$, compare the rank reported by `rref(A, tol=1e-10)` with `np.linalg.matrix_rank(A)` and with `np.linalg.det(A)`. Explain each disagreement, and say which answer is "right." (There is no single right answer; that is the point, and LA-13 is about it.)
