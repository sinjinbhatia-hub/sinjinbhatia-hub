# LA-12 Solutions

## A. Warm-up

**A1.** $E = \begin{pmatrix}2&1\\1&1\end{pmatrix}$, $\det E = 2 - 1 = 1$, $E^{-1} = \begin{pmatrix}1&-1\\-1&2\end{pmatrix}$. By (3.3) the dual basis vectors are the rows: $\mathbf{e}^1 = (1,-1)$, $\mathbf{e}^2 = (-1,2)$.

Checks: $\mathbf{e}^1(\mathbf{e}_1) = 2 - 1 = 1$, $\mathbf{e}^1(\mathbf{e}_2) = 1 - 1 = 0$, $\mathbf{e}^2(\mathbf{e}_1) = -2 + 2 = 0$, $\mathbf{e}^2(\mathbf{e}_2) = -1 + 2 = 1$. ✓

Rows of $E^\top$: $(2,1)$ and $(1,1)$. Then $(2,1)\cdot\mathbf{e}_2 = 2 + 1 = 3\ne 0$ and $(2,1)\cdot\mathbf{e}_1 = 5\ne 1$. They fail both conditions.

Coordinates of $\mathbf{v} = (3,2)^\top$: $v'^1 = \mathbf{e}^1(\mathbf{v}) = 3 - 2 = 1$, $v'^2 = \mathbf{e}^2(\mathbf{v}) = -3 + 4 = 1$. Check: $1\cdot(2,1)^\top + 1\cdot(1,1)^\top = (3,2)^\top$. ✓

**A2.** (a) $\varphi_i = \varphi(\mathbf{e}_i)$: $\varphi_1 = 2$, $\varphi_2 = -1$, $\varphi_3 = 3$. Over $\mathbb{R}$ there is no conjugate, so by Theorem 3.9 $\mathbf{v} = \sum_i\varphi(\mathbf{e}_i)\mathbf{e}_i = (2,-1,3)^\top$, and indeed $\langle\mathbf{v},\mathbf{u}\rangle = 2u^1 - u^2 + 3u^3$. ✓

(b) $\varphi(\mathbf{e}_1) = 1 + i$, $\varphi(\mathbf{e}_2) = -2i$. Theorem 3.9: $\mathbf{v} = \sum_i\overline{\varphi(\mathbf{e}_i)}\mathbf{e}_i = (1 - i,\ 2i)^\top$. Check with the first-slot conjugation: $\langle\mathbf{v},\mathbf{u}\rangle = \overline{(1-i)}u^1 + \overline{(2i)}u^2 = (1+i)u^1 - 2i\,u^2 = \varphi(\mathbf{u})$. ✓

Forgetting the conjugate gives $\mathbf{w} = (1+i, -2i)^\top$ and $\langle\mathbf{w},\mathbf{u}\rangle = (1-i)u^1 + 2i\,u^2$: the functional whose coefficients are the conjugates of $\varphi$'s, not $\varphi$ (test on $\mathbf{u} = \mathbf{e}_1$: it returns $1 - i$ where $\varphi$ returns $1 + i$). The Riesz map $\mathbf{v}\mapsto\langle\mathbf{v},\cdot\rangle$ is conjugate-linear, so the vector representing $\varphi$ carries the conjugated components.

**A3.** $P^{-1} = \operatorname{diag}(\tfrac12, 1)$. $MP = \begin{pmatrix}1&2\\0&1\end{pmatrix}\begin{pmatrix}2&0\\0&1\end{pmatrix} = \begin{pmatrix}2&2\\0&1\end{pmatrix}$.

Operator: $P^{-1}MP = \begin{pmatrix}\tfrac12&0\\0&1\end{pmatrix}\begin{pmatrix}2&2\\0&1\end{pmatrix} = \begin{pmatrix}1&1\\0&1\end{pmatrix}$.

Form: $P^\top MP = \begin{pmatrix}2&0\\0&1\end{pmatrix}\begin{pmatrix}2&2\\0&1\end{pmatrix} = \begin{pmatrix}4&4\\0&1\end{pmatrix}$.

The $(1,2)$ entry: as an operator, $M\mathbf{e}_2 = 2\mathbf{e}_1 + \mathbf{e}_2$, and in the new basis $2\mathbf{e}_1 = \mathbf{e}'_1$, so $M\mathbf{e}'_2 = \mathbf{e}'_1 + \mathbf{e}'_2$: the output is a vector, measured against a basis vector twice as long, so the coefficient *halves* ($2\to1$; contravariant output slot, $P^{-1}$). As a form, $B(\mathbf{e}'_1,\mathbf{e}'_2) = B(2\mathbf{e}_1,\mathbf{e}_2) = 2B(\mathbf{e}_1,\mathbf{e}_2) = 4$: the first slot *ate* a vector twice as long, so the reading *doubles* (covariant input slot, $P$). The $(1,1)$ entry: operator $M\mathbf{e}'_1 = 2M\mathbf{e}_1 = 2\mathbf{e}_1 = \mathbf{e}'_1$, coefficient $1$ (one $P$ from the input and one $P^{-1}$ from the output cancel); form $B(2\mathbf{e}_1, 2\mathbf{e}_1) = 4$ (two input slots, two factors of 2). The $(2,2)$ entry involves only $\mathbf{e}_2 = \mathbf{e}'_2$ and is unchanged in both.

**A4.** (a) By (3.12), $\sigma_x\otimes\sigma_z = \begin{pmatrix}0\cdot\sigma_z & 1\cdot\sigma_z\\ 1\cdot\sigma_z & 0\cdot\sigma_z\end{pmatrix} = \begin{pmatrix}0&0&1&0\\0&0&0&-1\\1&0&0&0\\0&-1&0&0\end{pmatrix}$.

$\sigma_x$ is the left factor, so it acts on qubit 0. Check: $(\sigma_x\otimes\sigma_z)|00\rangle = \sigma_x|0\rangle\otimes\sigma_z|0\rangle = |1\rangle\otimes|0\rangle = |10\rangle$, which has index $1\cdot2 + 0 = 2$, i.e. $(0,0,1,0)^\top$; and the first column of the matrix is $(0,0,1,0)^\top$. ✓ Qubit 0 flipped, qubit 1 picked up the sign $+1$ from $\sigma_z|0\rangle = |0\rangle$.

(b) In the basis $|00\rangle,|01\rangle,|10\rangle,|11\rangle$ (indices $0,1,2,3$), CNOT fixes indices $0,1$ and swaps $2\leftrightarrow3$:
$$\mathrm{CNOT} = \begin{pmatrix}1&0&0&0\\0&1&0&0\\0&0&0&1\\0&0&1&0\end{pmatrix} = \begin{pmatrix}I&0\\0&\sigma_x\end{pmatrix}.$$
If $\mathrm{CNOT} = A\otimes B$, then by (3.12) its $2\times2$ blocks are $A_{ij}B$. The $(1,1)$ block is $A_{11}B = I$ and the $(2,2)$ block is $A_{22}B = \sigma_x$. The first forces $B = I/A_{11}$ (so $A_{11}\ne0$), and then the second says $\sigma_x = (A_{22}/A_{11})I$, a multiple of the identity, which $\sigma_x$ is not. Contradiction. CNOT is not a Kronecker product, which is exactly why it can turn the product state $|+\rangle|0\rangle$ into the entangled Bell state: a Kronecker product $A\otimes B$ sends products to products by (3.10).

**A5.** (a) Let $U = \begin{pmatrix}a&b\\c&d\end{pmatrix}$ with $\det U = ad - bc = 1$. The $2\times2$ inverse formula gives $U^{-1} = \tfrac1{\det U}\begin{pmatrix}d&-b\\-c&a\end{pmatrix} = \begin{pmatrix}d&-b\\-c&a\end{pmatrix}$. Unitarity says $U^{-1} = U^\dagger = \begin{pmatrix}\bar a&\bar c\\ \bar b&\bar d\end{pmatrix}$. Comparing entries: $d = \bar a$ and $-c = \bar b$, i.e. $c = -\bar b$. So $U = \begin{pmatrix}a&b\\-\bar b&\bar a\end{pmatrix}$, and $\det U = a\bar a + b\bar b = |a|^2 + |b|^2 = 1$. Conversely any such matrix is unitary with determinant 1: column 1 is $(a, -\bar b)^\top$, column 2 is $(b, \bar a)^\top$, each has squared length $|a|^2 + |b|^2 = 1$, their inner product is $\bar a b + \overline{(-\bar b)}\,\bar a = \bar ab - b\bar a = 0$, and the determinant is $a\bar a + b\bar b = 1$. ∎

(b) With $\mathbf{n} = \mathbf{e}_z$, (A.7) gives $U = \cos\tfrac\theta2 I - i\sin\tfrac\theta2\sigma_z = \begin{pmatrix}\cos\tfrac\theta2 - i\sin\tfrac\theta2 & 0\\ 0 & \cos\tfrac\theta2 + i\sin\tfrac\theta2\end{pmatrix} = \begin{pmatrix}e^{-i\theta/2}&0\\0&e^{i\theta/2}\end{pmatrix}$. This is the form of (a) with $a = e^{-i\theta/2}$, $b = 0$: $\bar a = e^{i\theta/2}$ ✓, $|a|^2 = 1$ ✓.

## B. Core

**B1.** Define $\psi^i = \sum_k(P^{-1})_{ik}\mathbf{e}^k$. Evaluate on the new basis:
$$\psi^i(\mathbf{e}'_j) = \sum_k(P^{-1})_{ik}\,\mathbf{e}^k\!\Big(\sum_lP_{lj}\mathbf{e}_l\Big) = \sum_{k,l}(P^{-1})_{ik}P_{lj}\delta^k_l = \sum_k(P^{-1})_{ik}P_{kj} = (P^{-1}P)_{ij} = \delta_{ij}.$$
So the $\psi^i$ satisfy (3.2) for the primed basis. Theorem 3.2 says the dual basis is the *unique* family with this property, hence $\mathbf{e}'^i = \psi^i = \sum_k(P^{-1})_{ik}\mathbf{e}^k$. ∎

Components of a fixed $\varphi$: by definition $\varphi'_j = \varphi(\mathbf{e}'_j) = \varphi\!\big(\sum_iP_{ij}\mathbf{e}_i\big) = \sum_iP_{ij}\varphi(\mathbf{e}_i) = \sum_iP_{ij}\varphi_i = \sum_i(P^\top)_{ji}\varphi_i = (P^\top\boldsymbol\varphi)_j$. So $\boldsymbol\varphi' = P^\top\boldsymbol\varphi$. (Equivalently: $\varphi = \sum_i\varphi_i\mathbf{e}^i = \sum_i\varphi_i\sum_j P_{ij}\mathbf{e}'^j$ by inverting the first result, so $\varphi'_j = \sum_iP_{ij}\varphi_i$. Same answer.)

Orthonormal case: if the columns of $E$ are orthonormal, then $E^\top E = I$ (entry $(i,j)$ is $\mathbf{e}_i\cdot\mathbf{e}_j = \delta_{ij}$), so $E^{-1} = E^\top$. By (3.3) the dual basis vectors are the rows of $E^{-1} = E^\top$, which are the columns of $E$ written as rows, i.e. the basis vectors themselves: $\mathbf{e}^i(\mathbf{v}) = \mathbf{e}_i\cdot\mathbf{v}$. Conversely if the rows of $E^\top$ are the dual basis then $E^\top E = I$ and the basis is orthonormal. So "dual basis $=$ transpose" holds exactly for orthonormal bases. ∎

**B2.** (a) For $\psi\in U^*$: $(ST)^*\psi = \psi\circ(S\circ T) = (\psi\circ S)\circ T = T^*(\psi\circ S) = T^*(S^*\psi)$. So $(ST)^* = T^*\circ S^*$, with the order reversed. ∎

Matrices: let $A = [T]$, $B = [S]$ in some bases, so $[ST] = BA$. Theorem 3.5: $[(ST)^*] = [ST]^\top = (BA)^\top$. Also $[T^*S^*] = [T^*][S^*] = A^\top B^\top$ (matrix of a composition is the product; each factor by Theorem 3.5). Equating: $(BA)^\top = A^\top B^\top$. This is LA-1 Theorem 3.8, now with a reason: composition reverses under pulling back.

(b) By Definition 3.4 applied to $T^*$: for $\Lambda\in V^{**}$ and $\psi\in W^*$, $(T^{**}\Lambda)(\psi) = \Lambda(T^*\psi)$. Take $\Lambda = \Lambda_{\mathbf{v}}$: $(T^{**}\Lambda_{\mathbf{v}})(\psi) = \Lambda_{\mathbf{v}}(T^*\psi) = (T^*\psi)(\mathbf{v}) = \psi(T\mathbf{v}) = \Lambda_{T\mathbf{v}}(\psi)$. Since this holds for every $\psi$, $T^{**}\Lambda_{\mathbf{v}} = \Lambda_{T\mathbf{v}}$, i.e. $T^{**}\circ\Lambda_V = \Lambda_W\circ T$. After identifying $V$ with $V^{**}$ and $W$ with $W^{**}$ through $\Lambda$, this says $T^{**} = T$. ∎ The matrix identity: the matrix of $T^{**}$ in the double-dual bases (which Theorem 3.2 twice makes the original bases under $\Lambda$) is $(A^\top)^\top = A$.

**B3.** (a) Row index = qubit 0, column = qubit 1: $C = \tfrac1{\sqrt{10}}\begin{pmatrix}2&1\\1&2\end{pmatrix}$. $\det C = \tfrac1{10}(4 - 1) = \tfrac3{10}\ne0$, so rank 2: **entangled** by Theorem 3.15.

(b) $C^\dagger C = \tfrac1{10}\begin{pmatrix}2&1\\1&2\end{pmatrix}^2 = \tfrac1{10}\begin{pmatrix}5&4\\4&5\end{pmatrix}$. Eigenvalues: $\tfrac1{10}(5\pm4) = \tfrac9{10}, \tfrac1{10}$. Singular values are their square roots: $\sigma_1 = \tfrac3{\sqrt{10}}\approx0.9487$, $\sigma_2 = \tfrac1{\sqrt{10}}\approx0.3162$. These are the Schmidt coefficients; $\sigma_1^2 + \sigma_2^2 = 1$ ✓.

(c) $C$ is real symmetric with eigenvectors $\tfrac1{\sqrt2}(1,1)^\top$ (eigenvalue $\tfrac3{\sqrt{10}}$) and $\tfrac1{\sqrt2}(1,-1)^\top$ (eigenvalue $\tfrac1{\sqrt{10}}$), both positive, so the SVD is the eigendecomposition: $U = V = \tfrac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}$. By Theorem 3.16, $\mathbf{a}_1 = \mathbf{b}_1 = |+\rangle$, $\mathbf{a}_2 = \mathbf{b}_2 = |-\rangle$:
$$\psi = \tfrac3{\sqrt{10}}\,|+\rangle\otimes|+\rangle + \tfrac1{\sqrt{10}}\,|-\rangle\otimes|-\rangle.$$
Verify: $|+\rangle|+\rangle = \tfrac12(|00\rangle + |01\rangle + |10\rangle + |11\rangle)$, $|-\rangle|-\rangle = \tfrac12(|00\rangle - |01\rangle - |10\rangle + |11\rangle)$. Then $\tfrac1{2\sqrt{10}}\big[(3+1)|00\rangle + (3-1)|01\rangle + (3-1)|10\rangle + (3+1)|11\rangle\big] = \tfrac1{\sqrt{10}}(2|00\rangle + |01\rangle + |10\rangle + 2|11\rangle)$. ✓

(d) $S = -0.9\ln0.9 - 0.1\ln0.1 = 0.0948 + 0.2303 = 0.3251$ nats (about $0.469$ bits). Less than the maximum $\ln2 = 0.6931$: entangled, but not maximally.

(e) In the Schmidt form the state is a superposition of $|+\rangle|+\rangle$ and $|-\rangle|-\rangle$ only. Measuring both qubits in the $\pm$ basis: $P(+,+) = \sigma_1^2 = 0.9$, $P(-,-) = \sigma_2^2 = 0.1$, $P(+,-) = P(-,+) = 0$. The outcomes are perfectly correlated; this is what Schmidt form means physically, and QM-7 builds Bell's inequality out of such correlations.

**B4.** (a) $q = \mathbf{x}^\top B\mathbf{x}$ with the off-diagonal coefficients halved: $B = \begin{pmatrix}1&2&1\\2&2&0\\1&0&1\end{pmatrix}$.

(b) Complete the square in $x$: $(x + 2y + z)^2 = x^2 + 4y^2 + z^2 + 4xy + 2xz + 4yz$, so
$$q - (x+2y+z)^2 = 2y^2 + z^2 - 4y^2 - z^2 - 4yz = -2y^2 - 4yz = -2(y+z)^2 + 2z^2.$$
Hence $q = u^2 - v^2 + w^2$ with $u = x + 2y + z$, $v = \sqrt2(y+z)$, $w = \sqrt2z$. Signature $(p,q,z) = (2,1,0)$.

Explicit $P$: the new coordinates are $(u,v,w)^\top = Q(x,y,z)^\top$ with $Q = \begin{pmatrix}1&2&1\\0&\sqrt2&\sqrt2\\0&0&\sqrt2\end{pmatrix}$, and $q = \mathbf{x}^\top Q^\top DQ\,\mathbf{x}$ with $D = \operatorname{diag}(1,-1,1)$; so $B = Q^\top DQ$. Check: $DQ = \begin{pmatrix}1&2&1\\0&-\sqrt2&-\sqrt2\\0&0&\sqrt2\end{pmatrix}$ and $Q^\top(DQ)$ has rows $(1,2,1)$; $2(1,2,1) + \sqrt2(0,-\sqrt2,-\sqrt2) = (2,2,0)$; $(1,2,1) + \sqrt2(0,-\sqrt2,-\sqrt2) + \sqrt2(0,0,\sqrt2) = (1,0,1)$. That is $B$. ✓ Then $P = Q^{-1}$ gives $P^\top BP = Q^{-\top}Q^\top DQQ^{-1} = D$. Inverting $Q$ (back-substitute $z = w/\sqrt2$, $y = v/\sqrt2 - z$, $x = u - 2y - z$):
$$P = Q^{-1} = \begin{pmatrix}1 & -\sqrt2 & \tfrac1{\sqrt2}\\ 0 & \tfrac1{\sqrt2} & -\tfrac1{\sqrt2}\\ 0&0&\tfrac1{\sqrt2}\end{pmatrix}, \qquad P^\top BP = \operatorname{diag}(1,-1,1).$$
(The project checks this in Part 2: `change_basis_form(B4, P4)` with these entries returns $\operatorname{diag}(1,-1,1)$ to machine precision.)

(c) $\det B = 1(2\cdot1 - 0) - 2(2\cdot1 - 0) + 1(0 - 2) = 2 - 4 - 2 = -4 < 0$. The product of the three eigenvalues is negative, so an odd number are negative; $\operatorname{tr}B = 4 > 0$ rules out all three being negative. So exactly one negative eigenvalue: signature $(2,1,0)$, matching (b). (Numerically the eigenvalues are $-0.903, 1.194, 3.709$.)

(d) Theorem 3.12: the signature is invariant under every congruence $P^\top BP$. A basis in which $q = \alpha u^2 + \beta v^2 + \gamma w^2$ with all coefficients positive would have signature $(3,0,0)\ne(2,1,0)$. Impossible. A negative vector: take $v = 1$, $u = w = 0$, i.e. $z = 0$, $y = \tfrac1{\sqrt2}$, $x = -\sqrt2$: $q(-\sqrt2, \tfrac1{\sqrt2}, 0) = 2 + 4(-\sqrt2)(\tfrac1{\sqrt2}) + 0 + 2\cdot\tfrac12 + 0 = 2 - 4 + 1 = -1 < 0$. ✓

**B5.** (a) One upper index gets $P^{-1}$, two lower indices get $P$:
$$T'^i{}_{jk} = (P^{-1})^i{}_a\,P^b{}_j\,P^c{}_k\,T^a{}_{bc}.$$

(b) Substitute $T^a{}_{bc} = A^a{}_b\varphi_c$:
$$T'^i{}_{jk} = \big[(P^{-1})^i{}_aA^a{}_bP^b{}_j\big]\big[P^c{}_k\varphi_c\big] = (P^{-1}AP)^i{}_j\,(P^\top\boldsymbol\varphi)_k = A'^i{}_j\,\varphi'_k.$$
The transformed tensor is the tensor product of the transformed operator and the transformed functional. ∎ (This is a special case of a general fact: (3.18) is compatible with tensor products, because each index transforms independently.)

(c) Code (also in `project.py`, Part 5):
```python
rng = np.random.default_rng(0)
A = rng.standard_normal((3, 3)); phi = rng.standard_normal(3)
P = rng.standard_normal((3, 3)); Pinv = np.linalg.inv(P)
T = np.einsum('ab,c->abc', A, phi)                      # T^a_{bc} = A^a_b phi_c
T_law = np.einsum('ia,bj,ck,abc->ijk', Pinv, P, P, T)   # (3.18) for type (1,2)
T_parts = np.einsum('ij,k->ijk', Pinv @ A @ P, P.T @ phi)
assert np.allclose(T_law, T_parts)
```
The string `'ia,bj,ck,abc->ijk'` reads: $P^{-1}$ contracts its *second* index with the upper index $a$ of $T$ (output index $i$ first); each $P$ contracts its *first* index with a lower index of $T$ (output index $j$ or $k$ second). The usual mistake is writing `'ia,jb,kc,abc->ijk'`, which uses $P^\top$ in place of $P$ for the lower indices, giving $P^{-1}AP^\top$ for the operator part: wrong unless $P$ is orthogonal, which is why the bug is invisible on rotations and shows up only when $P$ is a shear or a stretch. Test with $P = \operatorname{diag}(2,1,1)$ if in doubt.

**B6.** (a) Both sides are antisymmetric under $j\leftrightarrow k$ (left side because $\varepsilon_{ikj} = -\varepsilon_{ijk}$; right side becomes $\delta_{km}\delta_{jn} - \delta_{kn}\delta_{jm}$, the negative) and under $m\leftrightarrow n$. So both vanish if $j = k$ or $m = n$, and it suffices to check $j\ne k$, $m\ne n$. In $\mathbb{R}^3$ the sets $\{j,k\}$ and $\{m,n\}$ are 2-element subsets of $\{1,2,3\}$.

*Case $\{j,k\} = \{m,n\}$.* On the left, the only $i$ with $\varepsilon_{ijk}\ne0$ is the single element $c\notin\{j,k\}$, and then $\varepsilon_{cmn}\ne0$ too. If $(m,n) = (j,k)$: left $= \varepsilon_{cjk}^2 = 1$; right $= 1\cdot1 - 0 = 1$. If $(m,n) = (k,j)$: left $= \varepsilon_{cjk}\varepsilon_{ckj} = -1$; right $= 0 - 1\cdot1 = -1$. ✓

*Case $\{j,k\}\ne\{m,n\}$.* Right side: $\delta_{jm}\delta_{kn}$ needs $j = m, k = n$ and $\delta_{jn}\delta_{km}$ needs $j = n, k = m$; either would make the sets equal, so the right side is $0$. Left side: a nonzero term needs $i\notin\{j,k\}$, so $i = c$, the complement of $\{j,k\}$. Two distinct 2-subsets of a 3-set share exactly one element, so the one element *not* in $\{j,k\}$ must be in $\{m,n\}$: $c\in\{m,n\}$, and $\varepsilon_{cmn} = 0$. Left side $0$. ✓ ∎

(b) $[\mathbf{a}\times(\mathbf{b}\times\mathbf{c})]_i = \varepsilon_{ijk}a_j(\mathbf{b}\times\mathbf{c})_k = \varepsilon_{ijk}a_j\varepsilon_{kmn}b_mc_n$. Cycle the first symbol, $\varepsilon_{ijk} = \varepsilon_{kij}$, and apply (a) with the summed index first: $\varepsilon_{kij}\varepsilon_{kmn} = \delta_{im}\delta_{jn} - \delta_{in}\delta_{jm}$. Then $(\delta_{im}\delta_{jn} - \delta_{in}\delta_{jm})a_jb_mc_n = b_i(a_jc_j) - c_i(a_jb_j) = b_i(\mathbf{a}\cdot\mathbf{c}) - c_i(\mathbf{a}\cdot\mathbf{b})$. ∎

(c) $|\mathbf{a}\times\mathbf{b}|^2 = (\mathbf{a}\times\mathbf{b})_i(\mathbf{a}\times\mathbf{b})_i = \varepsilon_{ijk}a_jb_k\,\varepsilon_{imn}a_mb_n = (\delta_{jm}\delta_{kn} - \delta_{jn}\delta_{km})a_jb_ka_mb_n = (a_ja_j)(b_kb_k) - (a_jb_j)(b_ka_k) = |\mathbf{a}|^2|\mathbf{b}|^2 - (\mathbf{a}\cdot\mathbf{b})^2$. ∎ (With $\mathbf{a}\cdot\mathbf{b} = |\mathbf{a}||\mathbf{b}|\cos\theta$ this is $|\mathbf{a}\times\mathbf{b}| = |\mathbf{a}||\mathbf{b}|\sin\theta$.)

**B7.** (a) $AX = AXI_m$, so by (3.13) $\operatorname{vec}(AX) = (I_m^\top\otimes A)\operatorname{vec}X = (I_m\otimes A)\operatorname{vec}X$. And $XB = I_nXB$, so $\operatorname{vec}(XB) = (B^\top\otimes I_n)\operatorname{vec}X$. Adding, $AX + XB = C$ becomes $(I_m\otimes A + B^\top\otimes I_n)\operatorname{vec}X = \operatorname{vec}C$: an $nm\times nm$ linear system. ∎

(b) Reference implementation (in `project.py`, Part 3):
```python
def solve_sylvester_vec(A, B, C):
    """Solve AX + XB = C via vec(AXB) = (B^T kron A) vec(X), eq. (3.13), column-major."""
    n, m = C.shape
    K = np.kron(np.eye(m), A) + np.kron(B.T, np.eye(n))
    x = np.linalg.solve(K, C.flatten(order="F"))
    return x.reshape(n, m, order="F")

A = rng.standard_normal((3, 3)); B = rng.standard_normal((2, 2)); C = rng.standard_normal((3, 2))
X = solve_sylvester_vec(A, B, C)
assert np.abs(A @ X + X @ B - C).max() < 1e-12
from scipy.linalg import solve_sylvester
assert np.allclose(X, solve_sylvester(A, B, C))
```
The `order="F"` is essential: (3.13) is a column-major identity, and `C.flatten()` (row-major) silently gives the wrong system.

(c) Let $A\mathbf{x} = \lambda\mathbf{x}$ and $B^\top\mathbf{y} = \mu\mathbf{y}$ ($B^\top$ has the same eigenvalues as $B$). By (3.10), $(I_m\otimes A)(\mathbf{y}\otimes\mathbf{x}) = \mathbf{y}\otimes\lambda\mathbf{x} = \lambda(\mathbf{y}\otimes\mathbf{x})$ and $(B^\top\otimes I_n)(\mathbf{y}\otimes\mathbf{x}) = \mu\mathbf{y}\otimes\mathbf{x} = \mu(\mathbf{y}\otimes\mathbf{x})$. So $\mathbf{y}\otimes\mathbf{x}$ is an eigenvector of the sum with eigenvalue $\lambda + \mu$. When $A$ and $B$ are diagonalizable these $nm$ vectors form a basis (Theorem 3.14), so the eigenvalues of $I\otimes A + B^\top\otimes I$ are exactly the $nm$ sums $\lambda_i + \mu_j$ (for the general case, use Schur as in Theorem 3.17(e): both summands become upper triangular simultaneously). The system is uniquely solvable iff no $\lambda_i + \mu_j = 0$, i.e. iff $A$ and $-B$ have no common eigenvalue. (This is why `solve_sylvester` fails for $B = -A$: then $X = I$ and $X = 2I$ both solve $AX + XB = 0$.)

(d) Row-major: $\operatorname{vec}_r(X) = \operatorname{vec}(X^\top)$. So $\operatorname{vec}_r(AXB) = \operatorname{vec}(B^\top X^\top A^\top) = ((A^\top)^\top\otimes B^\top)\operatorname{vec}(X^\top) = (A\otimes B^\top)\operatorname{vec}_r(X)$. The Sylvester system becomes $(A\otimes I_m + I_n\otimes B^\top)\operatorname{vec}_rX = \operatorname{vec}_rC$, with `C.ravel()` and `reshape(n, m)`. The project checks both forms on random matrices.

## C. Challenge

**C1.** (a) The map $\beta: V\times W\to W\otimes V$, $\beta(\mathbf{v},\mathbf{w}) = \mathbf{w}\otimes\mathbf{v}$, is bilinear (because $\otimes$ is bilinear in $W\times V$). By the universal property there is a unique linear $\tau: V\otimes W\to W\otimes V$ with $\tau(\mathbf{v}\otimes\mathbf{w}) = \mathbf{w}\otimes\mathbf{v}$. Symmetrically there is a linear $\tau': W\otimes V\to V\otimes W$ with $\tau'(\mathbf{w}\otimes\mathbf{v}) = \mathbf{v}\otimes\mathbf{w}$. Then $\tau'\circ\tau$ is a linear map $V\otimes W\to V\otimes W$ fixing every simple tensor; simple tensors span, and a linear map that is the identity on a spanning set is the identity. Likewise $\tau\circ\tau' = \mathrm{id}$. So $\tau$ is an isomorphism. ∎ (Note what was never needed: a basis, or a check that different expressions of the same tensor give the same image. The universal property did that.)

(b) For $\varphi\in V^*$, $\psi\in W^*$, the map $(\mathbf{v},\mathbf{w})\mapsto\varphi(\mathbf{v})\psi(\mathbf{w})$ is bilinear $V\times W\to\mathbb{F}$, so by the universal property it induces a unique linear functional $\Phi(\varphi,\psi)\in(V\otimes W)^*$ with $\Phi(\varphi,\psi)(\mathbf{v}\otimes\mathbf{w}) = \varphi(\mathbf{v})\psi(\mathbf{w})$. The assignment $(\varphi,\psi)\mapsto\Phi(\varphi,\psi)$ is bilinear: e.g. $\Phi(a\varphi + b\varphi',\psi)$ and $a\Phi(\varphi,\psi) + b\Phi(\varphi',\psi)$ agree on all simple tensors, hence are equal. So $\Phi$ induces a linear $\tilde\Phi: V^*\otimes W^*\to(V\otimes W)^*$ with $\tilde\Phi(\varphi\otimes\psi) = \Phi(\varphi,\psi)$. On the basis: $\tilde\Phi(\mathbf{e}^i\otimes\mathbf{f}^j)(\mathbf{e}_a\otimes\mathbf{f}_b) = \mathbf{e}^i(\mathbf{e}_a)\mathbf{f}^j(\mathbf{f}_b) = \delta^i_a\delta^j_b$. So $\tilde\Phi(\mathbf{e}^i\otimes\mathbf{f}^j)$ is exactly the $(i,j)$-th dual basis vector of the basis $\{\mathbf{e}_a\otimes\mathbf{f}_b\}$ of $V\otimes W$ (Theorem 3.14 plus Theorem 3.2). A linear map sending a basis of $V^*\otimes W^*$ (dimension $nm$, Theorem 3.14) to a basis of $(V\otimes W)^*$ (dimension $nm$, Theorem 3.2) is an isomorphism. ∎

**C2.** (a) A normalized 3-qubit state has $8$ complex amplitudes $= 16$ real numbers, minus $1$ for normalization and $1$ for the irrelevant global phase: $14$ real parameters. The proposed form: an orthonormal basis $\{\mathbf{a}_1,\mathbf{a}_2\}$ of $\mathbb{C}^2$ is an element of $U(2)$, $4$ real parameters, so the three bases give $12$; plus $2$ real $\sigma_k$ with $\sigma_1^2 + \sigma_2^2 = 1$, giving $1$ more: $13$ at most. But there is redundancy: multiplying $\mathbf{a}_k$ by a phase $e^{i\alpha_k}$, $\mathbf{b}_k$ by $e^{i\beta_k}$ and $\mathbf{c}_k$ by $e^{-i(\alpha_k + \beta_k)}$ leaves each term unchanged (2 parameters per $k$, 4 total), and a common phase on both $\mathbf{a}_k$ only changes the global phase of the state, which the count of $14$ has already discarded (1 more). So the family has at most $13 - 4 - 1 = 8$ dimensions, fewer than $14$. A generic state lies outside it. (The exact count is not the point; the point is that for two qubits the analogous count, $8 - 2 = 6$ versus $4 + 4 + 1 - 2 - 1 = 6$, comes out equal, and for three it does not.)

(b) Suppose $W = \sum_{k=1}^2\sigma_k\,\mathbf{a}_k\otimes\mathbf{b}_k\otimes\mathbf{c}_k$ with orthonormal $\{\mathbf{a}_k\}$, $\{\mathbf{b}_k\}$, $\{\mathbf{c}_k\}$. Group as qubit 0 versus qubits 12: $W = \sum_k\sigma_k\,\mathbf{a}_k\otimes(\mathbf{b}_k\otimes\mathbf{c}_k)$. The vectors $\mathbf{b}_k\otimes\mathbf{c}_k$ are orthonormal in $\mathbb{C}^4$ (inner product $\langle\mathbf{b}_1,\mathbf{b}_2\rangle\langle\mathbf{c}_1,\mathbf{c}_2\rangle = 0$, norms $1$), so this *is* a Schmidt decomposition of $W$ across that cut, with right Schmidt vectors that are product states.

Now compute the actual Schmidt decomposition across that cut. The coefficient matrix is $2\times4$ (row = qubit 0, column = qubits 12 in order $00,01,10,11$): $W = \tfrac1{\sqrt3}\big(|0\rangle\otimes(|01\rangle + |10\rangle) + |1\rangle\otimes|00\rangle\big)$, so
$$C = \tfrac1{\sqrt3}\begin{pmatrix}0&1&1&0\\1&0&0&0\end{pmatrix}.$$
The two rows are orthogonal with norms $\sqrt{2/3}$ and $\sqrt{1/3}$, so $CC^\dagger = \operatorname{diag}(\tfrac23,\tfrac13)$ and the singular values are $\sigma_1 = \sqrt{2/3}$, $\sigma_2 = \sqrt{1/3}$, **distinct**. Right singular vectors: the normalized rows, $\mathbf{r}_1 = \tfrac1{\sqrt2}(|01\rangle + |10\rangle)$ and $\mathbf{r}_2 = |00\rangle$. When the singular values are distinct the singular vectors are unique up to a phase (LA-9: they are the eigenvectors of $C^\dagger C$ for distinct nonzero eigenvalues). So any Schmidt decomposition of $W$ across this cut has right vectors $e^{i\alpha}\mathbf{r}_1$, $e^{i\beta}\mathbf{r}_2$. But $\mathbf{r}_1 = \tfrac1{\sqrt2}(|01\rangle + |10\rangle)$ has $2\times2$ coefficient matrix $\tfrac1{\sqrt2}\begin{pmatrix}0&1\\1&0\end{pmatrix}$ of rank 2: it is *not* a product $\mathbf{b}\otimes\mathbf{c}$ (Theorem 3.15). Contradiction. ∎

(c) GHZ $= \tfrac1{\sqrt2}|0\rangle\otimes|0\rangle\otimes|0\rangle + \tfrac1{\sqrt2}|1\rangle\otimes|1\rangle\otimes|1\rangle$ is already of the required form with $\sigma_1 = \sigma_2 = \tfrac1{\sqrt2}$ and $\mathbf{a}_k = \mathbf{b}_k = \mathbf{c}_k = |k\rangle$. So some three-party states have the form and some do not; there is no single "Schmidt rank" for three parties, and three-party entanglement comes in inequivalent kinds (GHZ-type and W-type are the two for three qubits).

What Theorem 3.16 really used: the identification of $V\otimes W$ with *matrices*, which have an SVD. Three factors give a $2\times2\times2$ array, and arrays with three indices have no SVD: the best analogue (the CP decomposition) does not have orthogonal factors in general and its rank is not even continuous. The two-particle case was using the fact that a bipartite state is a linear map $V^*\to W$ (§3.12), and linear maps are the one thing linear algebra fully understands.

**C3.** (a) $\{L_x, L_y, L_z\}$ is a basis of $\mathfrak{so}(3)$ (A.2) and $\{-\tfrac i2\sigma_x, -\tfrac i2\sigma_y, -\tfrac i2\sigma_z\}$ is a basis of $\mathfrak{su}(2)$ (A.5: anti-Hermitian traceless $2\times2$ matrices form a real 3-dimensional space, and these three are independent). A linear map sending a basis to a basis is a bijection. Brackets: for basis elements, (A.3) gives $[L_i, L_j] = \varepsilon_{ijk}L_k$ and (A.6) gives $[\phi(L_i),\phi(L_j)] = \varepsilon_{ijk}\phi(L_k) = \phi([L_i,L_j])$. Both sides of $\phi([X,Y]) = [\phi(X),\phi(Y)]$ are bilinear in $(X,Y)$, so agreement on basis pairs gives agreement everywhere. ∎

(b) $U(\theta) = \cos\tfrac\theta2I - i\sin\tfrac\theta2\,\mathbf{n}\cdot\boldsymbol\sigma$, so $\dot U(0) = -\tfrac i2\mathbf{n}\cdot\boldsymbol\sigma$ and $\dot U(0)^\dagger = \tfrac i2\mathbf{n}\cdot\boldsymbol\sigma$ (Pauli matrices are Hermitian). Product rule at $\theta = 0$ (where $U = I$):
$$\frac{d}{d\theta}\Big[U(\mathbf{v}\cdot\boldsymbol\sigma)U^\dagger\Big]_{\theta=0} = \dot U(0)(\mathbf{v}\cdot\boldsymbol\sigma) + (\mathbf{v}\cdot\boldsymbol\sigma)\dot U(0)^\dagger = -\tfrac i2\big[(\mathbf{n}\cdot\boldsymbol\sigma)(\mathbf{v}\cdot\boldsymbol\sigma) - (\mathbf{v}\cdot\boldsymbol\sigma)(\mathbf{n}\cdot\boldsymbol\sigma)\big] = -\tfrac i2[\mathbf{n}\cdot\boldsymbol\sigma, \mathbf{v}\cdot\boldsymbol\sigma].$$
By (A.5), $[\sigma_j,\sigma_k] = 2i\varepsilon_{jkl}\sigma_l$, so $[\mathbf{n}\cdot\boldsymbol\sigma,\mathbf{v}\cdot\boldsymbol\sigma] = n_jv_k\cdot2i\varepsilon_{jkl}\sigma_l = 2i(\mathbf{n}\times\mathbf{v})\cdot\boldsymbol\sigma$ by (3.23). Hence the derivative is $-\tfrac i2\cdot2i\,(\mathbf{n}\times\mathbf{v})\cdot\boldsymbol\sigma = (\mathbf{n}\times\mathbf{v})\cdot\boldsymbol\sigma$. ∎ Since $U(\mathbf{v}\cdot\boldsymbol\sigma)U^\dagger = (R(U)\mathbf{v})\cdot\boldsymbol\sigma$ by definition of $R$, and $\mathbf{v}\mapsto\mathbf{v}\cdot\boldsymbol\sigma$ is linear and injective, this says $\frac{d}{d\theta}R(U(\theta))\big|_0\mathbf{v} = \mathbf{n}\times\mathbf{v}$ for all $\mathbf{v}$, i.e. $\frac{d}{d\theta}R(U(\theta))\big|_0 = [\mathbf{n}]_\times$.

(c) Since $-\tfrac{i\theta}2\mathbf{n}\cdot\boldsymbol\sigma$ and $-\tfrac{is}2\mathbf{n}\cdot\boldsymbol\sigma$ commute, $U(\theta+s) = U(\theta)U(s)$ (LA-7's $e^Xe^Y = e^{X+Y}$ for commuting $X, Y$). $R$ is a homomorphism (A.6), so $R(U(\theta+s)) = R(U(\theta))R(U(s))$. Write $\mathcal{R}(\theta) = R(U(\theta))$. Differentiate in $s$ at $s = 0$: $\dot{\mathcal R}(\theta) = \mathcal{R}(\theta)\dot{\mathcal R}(0) = \mathcal{R}(\theta)[\mathbf{n}]_\times$ by (b), with $\mathcal{R}(0) = R(I) = I$. This is a linear ODE with constant coefficients; its unique solution (LA-7) is $\mathcal{R}(\theta) = \exp(\theta[\mathbf{n}]_\times)$, which by (A.4) is the rotation by $\theta$ about $\mathbf{n}$. ∎ So $R(U(\mathbf{n},\theta)) = R(\mathbf{n},\theta)$ for every axis, which is the claim of A.6; the demonstration there for $\mathbf{n} = \mathbf{e}_z$ was the special case. And $\theta\mapsto\theta + 2\pi$ changes $U$ to $-U$ (A.7) but not $\mathcal{R}$: the 2-to-1 property, now with a proof that nothing else is identified (if $R(U) = I$ then $U$ commutes with every $\mathbf{v}\cdot\boldsymbol\sigma$, hence with all $2\times2$ matrices, hence $U = cI$ with $c^2 = \det U = 1$, so $U = \pm I$).
