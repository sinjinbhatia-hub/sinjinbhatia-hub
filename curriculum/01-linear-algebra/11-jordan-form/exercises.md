# LA-11 Exercises

Conventions as in the notes: $J_s(\lambda) = \lambda I + N_s$ with ones on the superdiagonal, chains written bottom first, $d_k(\lambda) = \dim\ker(A-\lambda I)^k$, $p_A(x) = \det(xI - A)$ monic.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** Let $A = \begin{pmatrix}1&1\\-1&3\end{pmatrix}$. Find $p_A$, show there is only one eigenvector direction, find a Jordan chain $\mathbf{v}_2\to\mathbf{v}_1\to\mathbf{0}$, write $P$ and verify $P^{-1}AP = J_2(2)$. Then compute $e^{At}$ in closed form using (3.8) and (3.10), and check that its derivative at $t = 0$ is $A$.

**A2.** Let $N = \begin{pmatrix}0&a&b\\0&0&c\\0&0&0\end{pmatrix}$. Compute $N^2$ and $N^3$. For which $(a,b,c)$ is the index of nilpotency $3$, $2$, $1$? Give the Jordan form of $N$ in each case. Then prove that *every* strictly upper triangular $n\times n$ matrix satisfies $N^n = 0$, using the flag $V_k = \operatorname{span}(\mathbf{e}_1,\dots,\mathbf{e}_k)$.

**A3.** An eigenvalue $\lambda$ of some $6\times6$ matrix has kernel dimensions $d_0, d_1, d_2, d_3, d_4 = 0, 3, 5, 6, 6$. Find the Jordan block sizes for $\lambda$, the geometric and algebraic multiplicities, and the exponent of $(x-\lambda)$ in the minimal polynomial.

**A4.** For $J = J_3(0)\oplus J_1(0)\oplus J_2(5)$ (a $6\times6$ matrix), write down $p_J$, $m_J$, the sequences $d_k(0)$ and $d_k(5)$, and say whether $J$ is diagonalizable.

**A5.** Write out $X(t) = e^{tJ_3(\lambda)}$ from (3.9) and verify directly that $X'(t) = J_3(\lambda)X(t)$ and $X(0) = I$. Then evaluate $e^{J}$ for $J = \begin{pmatrix}2&1\\0&2\end{pmatrix}$ to four decimals.

**A6.** For $A = \begin{pmatrix}a&b\\c&d\end{pmatrix}$, write $p_A(x) = x^2 - (a+d)x + (ad-bc)$ and verify Cayley–Hamilton by multiplying out $A^2 - (a+d)A + (ad-bc)I$. Deduce $A^{-1} = \dfrac{(a+d)I - A}{ad-bc}$ and check it is the formula you already know.

## B. Core (do all)

**B1.** (The splitting lemma.) Let $T$ be any operator on $V$ with $\dim V = n$. Prove $V = \ker T^n\oplus\operatorname{im}T^n$. (Hint: for the intersection, if $\mathbf{v} = T^n\mathbf{u}$ and $T^n\mathbf{v} = \mathbf{0}$, what does Lemma 3.3(c) say about $\mathbf{u}$? For the sum, count.) Then apply it to $T = A - \lambda I$: conclude $V = G_\lambda\oplus U_\lambda$ with $U_\lambda = \operatorname{im}(A-\lambda I)^n$ invariant under $A$, and that $\lambda$ is *not* an eigenvalue of $A|_{U_\lambda}$. (This lemma is the engine of Axler's proof of Theorem 3.9, which needs no triangular form; see §10.)

**B2.** Let $A = \begin{pmatrix}2&1&-1\\0&1&0\\1&1&0\end{pmatrix}$. Find $p_A$, the sequence $d_k$ for each eigenvalue, the Jordan block sizes, and a Jordan basis $P$; verify $P^{-1}AP = J$ by checking the chain relations $A\mathbf{v}_j = \lambda\mathbf{v}_j + \mathbf{v}_{j-1}$ column by column. Write the general solution of $\mathbf{x}' = A\mathbf{x}$.

**B3.** (What the $d_k$ can be.) Let $d_k = \dim\ker N^k$ for a nilpotent $N$.
(a) Prove that the number of blocks of size exactly $k$ is $2d_k - d_{k-1} - d_{k+1}$.
(b) Prove that the jumps $d_k - d_{k-1}$ are non-increasing in $k$ (so the sequence $d_k$ is concave).
(c) Decide whether $(0, 2, 4, 5)$ and $(0, 1, 3, 4)$ can occur as $(d_0, d_1, d_2, d_3)$; where they can, give the block sizes.
(d) Give a $4\times4$ nilpotent matrix with $d = (0, 2, 4)$ that is *not* already in Jordan form, and find its Jordan basis.

**B4.** (Diagonalizable iff distinct roots, without the Jordan form.) Suppose $m_A(x) = \prod_{i=1}^m(x-\lambda_i)$ with the $\lambda_i$ distinct. Let $\ell_i(x) = \prod_{j\ne i}\dfrac{x-\lambda_j}{\lambda_i-\lambda_j}$ be the Lagrange polynomials.
(a) Show $\sum_i\ell_i(x) = 1$ identically. (Hint: degree and number of roots.)
(b) Show $(A-\lambda_iI)\,\ell_i(A) = 0$.
(c) Conclude that every $\mathbf{v}$ is a sum of eigenvectors, hence $A$ is diagonalizable.
(d) Prove the converse directly: if $A = PDP^{-1}$ with $D$ diagonal, then $m_A$ has distinct roots.

**B5.** (Companion matrices and $t^2e^t$.) Let $C = \begin{pmatrix}0&1&0\\0&0&1\\1&-3&3\end{pmatrix}$, the companion matrix of $y''' - 3y'' + 3y' - y = 0$ (with $\mathbf{x} = (y, y', y'')^\top$, $\mathbf{x}' = C\mathbf{x}$).
(a) Show $p_C = (x-1)^3$, that $C$ has a one-dimensional eigenspace, and hence that its Jordan form is $J_3(1)$ and $m_C = p_C$.
(b) Use Cayley–Hamilton to write $C^{-1}$ as a polynomial in $C$ and compute it. Then show $C^{10} = 45C^2 - 80C + 36I$ in two ways: by reducing $x^{10}$ modulo $(x-1)^3$, and by expanding $(I + N)^{10}$ with $N = C - I$. Compute the matrix.
(c) Find a Jordan chain for $C$ (top first) and use (3.9) to write the general solution of the ODE; identify $e^t, te^t, t^2e^t$.
(d) Prove that for *any* companion matrix (last row $(-c_0, -c_1, \dots, -c_{n-1})$, ones on the superdiagonal, zeros elsewhere) every eigenvalue has geometric multiplicity $1$, so there is one Jordan block per eigenvalue and $m_C = p_C$. (Such matrices are called *nonderogatory*; this is why repeated roots in ODEs and AR processes always produce $t^k$ terms.)

**B6.** (Symmetric matrices and the three-mass ring.)
(a) Prove directly, for real symmetric $A$ and real $\lambda$: $\ker(A-\lambda I)^2 = \ker(A-\lambda I)$. Conclude $A$ has no Jordan block of size $\ge2$.
(b) Three equal masses joined in a ring by three equal springs have stiffness matrix $K = \begin{pmatrix}2&-1&-1\\-1&2&-1\\-1&-1&2\end{pmatrix}$ (units $k = m = 1$), and $\ddot{\mathbf{x}} = -K\mathbf{x}$. Find the eigenvalues of $K$ and a basis of each eigenspace. The frequency $\sqrt3$ is repeated: explain in one sentence why there is nevertheless no $t\cos\sqrt3t$ solution.
(c) Write $\ddot{\mathbf{x}} = -K\mathbf{x}$ as a first-order system $\dot{\mathbf{z}} = A\mathbf{z}$ with $A = \begin{pmatrix}0&I\\-K&0\end{pmatrix}$. Show that the eigenvalue $0$ of $A$ has algebraic multiplicity $2$ and geometric multiplicity $1$, exhibit the chain, and interpret the resulting solution physically. Give the full Jordan form of $A$.

**B7.** (Code.) In `project.py`, Part 6 defines `minimal_polynomial_by_dependence(A)`. Write your own version before reading it: build the $n^2\times(d+1)$ matrix whose columns are the flattened $I, A, \dots, A^d$, and return the monic polynomial from the first nontrivial null vector. Test that it agrees with `minimal_polynomial_from_blocks` on every matrix in the script. Then add `is_diagonalizable(A)` that uses *only* $m_A$ (hint: $m_A$ has a repeated root iff $\gcd(m_A, m_A') \ne 1$; `sympy.gcd` does polynomial gcds) and confirm it on `A_ode` (no), `A5` (yes), `A45` (no), `K_ring` (yes). On which of the script's matrices is $\deg m_A < \deg p_A$?

## C. Challenge (do at least one)

**C1.** (Triangularization by induction.) Prove without inner products that every $A\in\mathbb{C}^{n\times n}$ is similar to an upper triangular matrix. (Pick an eigenvector $\mathbf{v}_1$ with eigenvalue $\lambda$, extend to a basis, and look at the matrix of $A$ in that basis; induct on the lower-right block.) Show that the diagonal entries are the eigenvalues with algebraic multiplicity and that, by choosing which eigenvector to take at each step, you can arrange them in any order. Confirm that this is all Theorems 3.9(c) and 3.16 used.

**C2.** (The $\varepsilon$-demo, exactly.) Let $A_\varepsilon = \begin{pmatrix}1&1\\ \varepsilon&1\end{pmatrix}$ with $\varepsilon>0$.
(a) Compute its eigenvalues and the eigenvector matrix $V = \begin{pmatrix}1&1\\ \sqrt\varepsilon&-\sqrt\varepsilon\end{pmatrix}$.
(b) Compute $V^\top V$ and its eigenvalues, and deduce that the singular values of $V$ are $\sqrt2$ and $\sqrt{2\varepsilon}$, so $\kappa_2(V) = 1/\sqrt\varepsilon$ exactly.
(c) Show the angle between the eigenvectors is $2\arctan\sqrt\varepsilon$.
(d) For the other perturbation $B_\varepsilon = \begin{pmatrix}1&1\\0&1+\varepsilon\end{pmatrix}$, find the eigenvectors and show $\kappa_2$ of its eigenvector matrix grows like $2/\varepsilon$, which is worse.
(e) At $\varepsilon = 10^{-16}$, how far apart are the eigenvalues of $A_\varepsilon$, and how does that compare with the size of the perturbation?

**C3.** (The adjugate proof of Cayley–Hamilton.) Let $B(x) = \operatorname{adj}(xI - A)$, so that $B(x)(xI - A) = p_A(x)I$ (LA-6).
(a) Explain why each entry of $B(x)$ is a polynomial in $x$ of degree $\le n-1$, so $B(x) = \sum_{k=0}^{n-1}B_kx^k$ with constant matrices $B_k$.
(b) Write $p_A(x) = \sum_{k=0}^n c_kx^k$ ($c_n = 1$). Expand $B(x)(xI - A)$ and match coefficients of $x^k$ on both sides to get $n+1$ matrix identities.
(c) Multiply the $k$-th identity on the right by $A^k$ and add them all up. Show the right side telescopes to $0$ and the left side is $p_A(A)$.
(d) State precisely why this is not the fake proof of §3.8.

**C4.** (Every matrix function is a polynomial in $A$.)
(a) Using (3.7) and (3.10), prove: if $q$ is a polynomial with $q^{(k)}(\lambda) = f^{(k)}(\lambda)$ for every eigenvalue $\lambda$ and every $0\le k<s_{\max}(\lambda)$, then $f(A) = q(A)$.
(b) Show that such a $q$ exists with $\deg q < \deg m_A$ and is unique. (Count conditions versus unknowns, and show the only polynomial of degree $<\deg m_A$ with all those derivatives zero is $0$.)
(c) For the matrix $A$ of Example 4.1, find the degree-$1$ polynomial $q$ with $q(1) = e^t$, $q'(1) = te^t$, and show $e^{tA} = e^t\big((1-t)I + tA\big)$. Compare with the $e^{At}$ computed in Example 4.1.
(d) When $A$ has $n$ distinct eigenvalues, what is $q$? Which classical formula is it?
