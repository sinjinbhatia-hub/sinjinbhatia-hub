# LA-8: Inner Product Spaces and the Spectral Theorem

**Time:** ~14 hours · **Prerequisites:** LA-5 (orthogonality, projections, Gram–Schmidt, QR), LA-7 (eigenvalues, diagonalization, complex eigenvalues) · **Unlocks:** LA-9 (SVD), LA-12 (bilinear forms, tensors), F2 (Fourier series), F5 (multivariate Gaussian), CM-6 (normal modes), CM-7 (rigid bodies), QM-1, QM-4, QC-1

> **Conventions in this module.** Vectors are columns in bold lowercase, $\mathbf{v} \in \mathbb{R}^n$ or $\mathbb{C}^n$; matrices are uppercase. $A^\top$ is the transpose, $\bar{A}$ the entrywise complex conjugate, and $A^\dagger = \bar{A}^\top$ the conjugate transpose ("dagger"). $\langle u, v\rangle$ is an abstract inner product. **Complex inner products are conjugate-linear in the first slot and linear in the second:** $\langle c\,u, v\rangle = \bar c\,\langle u, v\rangle$ and $\langle u, c\,v\rangle = c\,\langle u, v\rangle$. On $\mathbb{C}^n$ the standard one is therefore $\langle \mathbf{u}, \mathbf{v}\rangle = \mathbf{u}^\dagger \mathbf{v} = \sum_i \bar u_i v_i$. This is the physics convention, and it is chosen so that Dirac's $\langle u | v\rangle$ means exactly $\langle u, v\rangle$. Axler and most pure-mathematics texts put the linearity in the *first* slot; when you read them, swap the slots. $A^*$ denotes the adjoint of an operator with respect to whatever inner product is in force; for matrices with the standard inner product, $A^* = A^\top$ (real) or $A^\dagger$ (complex). *Symmetric* means $A^\top = A$; *Hermitian* means $A^\dagger = A$; *orthogonal* means $Q^\top Q = I$; *unitary* means $U^\dagger U = I$. Eigenvalues of a symmetric or Hermitian matrix are real and we order them $\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_n$, so $\lambda_1 = \lambda_{\max}$; NumPy's `eigh` returns them in *increasing* order. $\operatorname{tr}$ is the trace (sum of diagonal entries), $\delta_{ij}$ is the Kronecker delta (1 if $i = j$, else 0). Pauli matrices: $\sigma_x = \begin{pmatrix}0&1\\1&0\end{pmatrix}$, $\sigma_y = \begin{pmatrix}0&-i\\i&0\end{pmatrix}$, $\sigma_z = \begin{pmatrix}1&0\\0&-1\end{pmatrix}$. Stresses are in MPa with tension positive; the shear stress on a plane with unit normal $\mathbf{n}$ is $\tau_n = \mathbf{t}^\top\sigma\mathbf{n}$ with $\mathbf{t} = (-n_y, n_x)^\top$ ($\mathbf{n}$ rotated $+90°$), plotted upward on Mohr's circle (many engineering texts plot the opposite sign). Qubit basis: $|0\rangle = \mathbf{e}_1$, $|1\rangle = \mathbf{e}_2$.

## 0. Where this sits

You know what an inner product is in $\mathbb{R}^n$ (LA-1, LA-5): the dot product, and the three properties it has. You know how to project, how to run Gram–Schmidt, and what an orthonormal basis buys you (LA-5). You know eigenvalues, diagonalization $A = X\Lambda X^{-1}$, and that some matrices have complex eigenvalues or cannot be diagonalized at all (LA-7).

This module does two things. First, it takes the three dot-product properties as *axioms*, so that functions, polynomials, and matrices get lengths and angles, and Fourier series (F2) become orthogonal projection. Second, it proves the single most used theorem in applied linear algebra: a **symmetric matrix has an orthonormal basis of eigenvectors with real eigenvalues**, $A = Q\Lambda Q^\top$. Everything downstream leans on it: the SVD in LA-9 is this theorem applied to $A^\top A$; covariance matrices in F5 and ML-6 are symmetric positive semidefinite; the inertia tensor in CM-7 and the stress tensor in solid mechanics have principal axes because of it; quantum observables in QM-1 and QM-4 are Hermitian operators and their measurable values are the eigenvalues in this theorem.

## 1. The problem that forced this idea

**Heat, 1807.** Joseph Fourier submitted a memoir on the propagation of heat in solids to the Institut de France. To solve the heat equation on a rod he claimed that an essentially arbitrary function $f$ on $[-\pi, \pi]$ could be written as a sum of sines and cosines,
$$f(x) = \frac{a_0}{2} + \sum_{n \ge 1} \left(a_n \cos nx + b_n \sin nx\right),$$
and he gave formulas for the coefficients: $b_n = \frac{1}{\pi}\int_{-\pi}^{\pi} f(x)\sin nx\,dx$. He first got them by truncating an infinite system of linear equations for the coefficients and letting the size grow. The argument that survived is shorter, and Euler had already used it in 1777: multiply both sides by $\sin mx$, integrate, and use the fact that $\int_{-\pi}^{\pi}\sin mx \sin nx\,dx = 0$ when $m \neq n$ and $= \pi$ when $m = n$, so every term dies except one. The committee (Lagrange, Laplace, Monge, Lacroix) was unconvinced by the convergence claims, and the full *Théorie analytique de la chaleur* did not appear until 1822. But look at what the coefficient formula *is*. It is "how much of $f$ points along $\sin nx$," divided by "the squared length of $\sin nx$" ($\int\sin^2 nx\,dx = \pi$): the projection formula $t^* = (\mathbf{u}\cdot\mathbf{v})/\|\mathbf{v}\|^2$ from LA-1 Example 4.5, with the dot product replaced by an integral and the vectors replaced by functions. Fourier was projecting onto an orthogonal basis of a function space. He did not have the words.

**Planets, 1829.** The words were slow in coming, but the matrix theorem was not. Cauchy, studying the "secular inequalities" of planetary orbits (a problem that reduces to a symmetric matrix), proved in 1829 that the roots of the characteristic equation of a symmetric matrix are all real, and that the eigenvectors for distinct roots are perpendicular. In the language of the time this was the *principal axes theorem*: a quadratic form $\sum a_{ij}x_i x_j$ can be rotated into a sum of squares $\sum \lambda_i y_i^2$. Hermite extended the reality of the roots to complex matrices with $A^\dagger = A$ in 1855, which is why they carry his name.

**Integral equations, 1904–1910.** Hilbert, working on Fredholm's integral equations, treated the kernel $K(s, t)$ of $\int K(s,t)\phi(t)\,dt = \lambda\phi(s)$ as a symmetric "matrix with infinitely many rows," and proved a version of the principal axes theorem for it. He adopted the word *spectrum* for the set of eigenvalues (Wirtinger had used it for a similar set in 1897, in work on Hill's equation), and it stuck. His student Erhard Schmidt, in papers of 1907 and 1908, supplied the geometry: he treated square-summable sequences as vectors, defined their lengths and inner products, spoke of orthogonality and projection, and wrote down the orthogonalization procedure that LA-5 calls Gram–Schmidt. The space Schmidt described is what everyone now calls a *Hilbert space*. The move was: forget that the vectors are arrows, keep only the three properties of the dot product, and see how much of Euclidean geometry survives. The answer is: all of it that matters.

**Quantum mechanics, 1925–1927.** Heisenberg's 1925 paper represented position and momentum by arrays of numbers that Born recognized as matrices. Born and Jordan, then Born, Heisenberg, and Jordan (1925–26), required the arrays for physical quantities to satisfy $A^\dagger = A$, and showed that the energy levels are the eigenvalues of the energy matrix, found by rotating it to principal axes. The general rule, that the possible measured values of *any* quantity are the eigenvalues of its matrix, came with Dirac's and Jordan's transformation theory (1926–27) and von Neumann's Hilbert-space formulation (1927). Hermite's theorem is why those values are real numbers, as measurements must be. Born's statistical interpretation (1926) made the probability of observing the eigenvalue $\lambda_i$ equal to $|\langle q_i, \psi\rangle|^2$, where $q_i$ is the unit eigenvector. For these probabilities to sum to 1, the eigenvectors have to be orthonormal and complete. That is the spectral theorem. Von Neumann (1927–1932) built out that framework and proved the theorem for the infinite-dimensional operators QM actually uses. The point you should take from the history: the spectral theorem is not decorative. It is the reason measurement outcomes are well-defined numbers with well-defined probabilities.

## 2. Intuition

### Angles between functions

Sample a function $f$ on $[a, b]$ at $n$ equally spaced points $x_1, \dots, x_n$. You get a vector $\mathbf{f} = (f(x_1), \dots, f(x_n))^\top \in \mathbb{R}^n$. Do the same for $g$. Their dot product, scaled by the spacing $\Delta x = (b-a)/n$, is
$$\Delta x\;\mathbf{f}\cdot\mathbf{g} = \sum_{k=1}^n f(x_k)\,g(x_k)\,\Delta x \;\xrightarrow{\;n\to\infty\;}\; \int_a^b f(x)g(x)\,dx.$$
A Riemann sum. So $\int f g\,dx$ is the dot product of two functions, in the limit where you sample them everywhere. Lengths, angles, orthogonality, projections: all of LA-5 should go through for functions, with $\int$ in place of $\sum$.

Try it: on $[-\pi, \pi]$, $\int \sin x\,\sin 2x\,dx = 0$. The functions $\sin x$ and $\sin 2x$ are *perpendicular*. On $[-1, 1]$, $\int 1\cdot x\,dx = 0$: the constant function and the identity function are perpendicular. But $\int 1\cdot x^2\,dx = 2/3 \ne 0$, so $x^2$ has a component along $1$. Subtract it off and you get $x^2 - \tfrac13$, which *is* perpendicular to $1$. That is Gram–Schmidt on polynomials, and §4 does it in full.

### Principal axes of an ellipse

The curve $5x^2 + 4xy + 2y^2 = 1$ is an ellipse, tilted. In the right rotated coordinates $(y_1, y_2)$ it has no cross term: $6y_1^2 + 1\,y_2^2 = 1$. Its axes are perpendicular lines through the origin, and along them the curve is at distance $1/\sqrt{6}$ and $1/\sqrt{1}$ from the center.

```
                y
                |
        \...    |
       ... .....|
      ..  \    ...
      .    \    |..
      ..    \   | ...
       .     \  |   ..
       .      \ |    ..===     minor axis (===) along (2, 1),
        .      \| ====..         semi-axis 1/sqrt(6) ≈ 0.41
    ----..----=====----..---- x
         ..==== |\      ..
       ===..    | \      .     major axis (\\\) along (1, -2),
           ..   |  \     .       semi-axis 1
            ... |   \    ..
              ..|    \    .
               ...    \  ..
                |..... ...
                |    ...\
```
(Drawn to scale, with each character cell twice as tall as it is wide. The project's `quadratic_forms.png` draws the same ellipse properly.)

Now write the quadratic as a matrix product: $5x^2 + 4xy + 2y^2 = \begin{pmatrix}x & y\end{pmatrix}\begin{pmatrix}5 & 2\\2 & 2\end{pmatrix}\begin{pmatrix}x\\y\end{pmatrix}$, splitting the cross coefficient $4$ as $2 + 2$ so the matrix is symmetric. The numbers $6$ and $1$ are its eigenvalues (trace $7$, determinant $6$). The axes are its eigenvectors. Changing to the axis coordinates is exactly diagonalizing the matrix, and the fact that the axes are *perpendicular* is the fact that the eigenvectors of a symmetric matrix are orthogonal.

### A symmetric matrix stretches without twisting

Apply any invertible $2\times 2$ matrix $A$ to the unit circle. You get an ellipse (LA-9 proves this for every invertible matrix; a singular one flattens the circle to a segment or a point). For a *symmetric* $A$, something extra is true: the two perpendicular directions that get mapped onto the ellipse's axes **are** the axes. $A$ stretches the plane along two perpendicular lines and does nothing else. No rotation, no shear.

```
   Symmetric A = [[2,1],[1,2]]:            Shear A = [[1,1],[0,1]]:

       ^                 ^                      ^                ^
     . | .             .   .                  . | .            .   .  .
   .   |   .    -->  .   |   .              .   |   .   -->  .    |     .
   .   +   .        .    +    .             .   +   .       .     +      .
     . | .           .   |   .                . | .           .   |    .
       v               .   .                    v               .  .
   circle           ellipse, axes along         circle          ellipse, but the vectors
                    (1,1) and (1,-1),                           that land on its axes were
                    which were ALREADY                          NOT the axes beforehand
                    the directions stretched
```

For the shear, the unit circle still becomes an ellipse, but the preimages of the ellipse's axes are some other pair of perpendicular directions, and the eigenvectors (there is only one direction, $(1,0)^\top$) have nothing to do with the axes. For a symmetric matrix the three notions coincide: eigenvectors, directions of pure stretch, axes of the image ellipse.

### Guess the theorem

Put the three pictures together and guess. A symmetric $n\times n$ matrix should have $n$ **mutually perpendicular eigenvectors**, with **real eigenvalues** (stretch factors, possibly negative). Writing those unit eigenvectors as the columns of $Q$, you'd expect $Q^\top Q = I$ and $A = Q\Lambda Q^\top$: rotate to the eigen-axes, stretch, rotate back. The quadratic form $\mathbf{x}^\top A\mathbf{x}$ becomes $\sum\lambda_i y_i^2$ in those axes, so its level sets are ellipses when all $\lambda_i > 0$ and hyperbolas when the signs mix. And for a *rotation* (an orthogonal matrix) nothing gets stretched, so its eigenvalues should all have absolute value $1$.

Every one of these guesses is a theorem in §3.

## 3. Formal development

### 3.1 Inner products

**Definition 3.1 (Real inner product).** An inner product on a real vector space $V$ is a function $\langle\cdot,\cdot\rangle : V \times V \to \mathbb{R}$ such that for all $u, v, w \in V$ and $c \in \mathbb{R}$:
1. (symmetric) $\langle u, v\rangle = \langle v, u\rangle$;
2. (linear in each slot) $\langle u, cv + w\rangle = c\langle u, v\rangle + \langle u, w\rangle$, and by symmetry the same in the first slot;
3. (positive definite) $\langle v, v\rangle \ge 0$, with equality only when $v = 0$.

These are exactly the three dot-product properties listed after LA-1 Definition 3.9, now taken as axioms.

**Definition 3.2 (Complex inner product).** On a complex vector space, $\langle\cdot,\cdot\rangle : V\times V\to\mathbb{C}$ with
1. (conjugate symmetric) $\langle u, v\rangle = \overline{\langle v, u\rangle}$;
2. (linear in the **second** slot) $\langle u, cv + w\rangle = c\langle u, v\rangle + \langle u, w\rangle$;
3. (positive definite) $\langle v, v\rangle \ge 0$ (it is real by axiom 1), with equality only when $v = 0$.

Axioms 1 and 2 force the first slot to be **conjugate-linear**: $\langle cu, v\rangle = \overline{\langle v, cu\rangle} = \overline{c\langle v, u\rangle} = \bar c\,\langle u, v\rangle$. This is the convention flagged in the box: conjugate the first slot. The reason a conjugate is needed at all: with $\langle \mathbf{u}, \mathbf{v}\rangle = \sum u_i v_i$ and no conjugate, $\mathbf{v} = (1, i)^\top$ would have $\langle\mathbf{v},\mathbf{v}\rangle = 1 + i^2 = 0$, so axiom 3 would fail. The reason to put the conjugate on the *first* slot rather than the second is purely conventional, and we follow physics so that Dirac's $\langle\phi|\psi\rangle$ (§3.12) means $\langle\phi, \psi\rangle$ with no translation. Axler conjugates the second slot. Both are inner products; do not mix them in one computation.

**Examples.** Each is an inner product; checking the axioms is a short exercise in each case.

- *Standard.* On $\mathbb{R}^n$: $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^\top\mathbf{v}$. On $\mathbb{C}^n$:
$$\langle \mathbf{u}, \mathbf{v}\rangle = \mathbf{u}^\dagger\mathbf{v} = \sum_{i=1}^n \bar u_i v_i. \tag{3.1}$$
Conjugate-linear in $\mathbf{u}$, linear in $\mathbf{v}$, and $\langle\mathbf{v},\mathbf{v}\rangle = \sum|v_i|^2 > 0$ for $\mathbf{v}\ne\mathbf{0}$.
- *Weighted.* For positive weights $w_i$, $\langle\mathbf{u},\mathbf{v}\rangle_w = \sum_i w_i \bar u_i v_i$ (LA-1 Exercise C1). More generally $\langle\mathbf{u},\mathbf{v}\rangle_W = \mathbf{u}^\dagger W\mathbf{v}$ for any Hermitian positive definite $W$ (defined after Theorem 3.14; the tests are in §3.8); Exercise B4 shows *every* inner product on $\mathbb{C}^n$ has this form.
- *Functions.* On continuous complex-valued functions on $[a, b]$:
$$\langle f, g\rangle = \int_a^b \overline{f(x)}\,g(x)\,dx. \tag{3.2}$$
Axiom 3 needs continuity: a continuous $f$ with $\int|f|^2 = 0$ is identically zero, because if $|f(x_0)| > 0$ at one point then, by continuity, $|f(x)|^2 > \tfrac12|f(x_0)|^2$ on some interval around $x_0$, making the integral positive. F2 uses (3.2) on $[-\pi,\pi]$; Example 4.1 uses the real version on $[-1,1]$.
- *Matrices (Frobenius).* On $m\times n$ complex matrices:
$$\langle A, B\rangle_F = \operatorname{tr}(A^\dagger B) = \sum_{i,j} \bar a_{ij} b_{ij}. \tag{3.3}$$
The second equality: $(A^\dagger B)_{jj} = \sum_i (A^\dagger)_{ji} b_{ij} = \sum_i \bar a_{ij} b_{ij}$, then sum over $j$. So (3.3) is (3.1) applied to the matrix flattened into a vector of length $mn$, and the axioms follow. The norm it induces is $\|A\|_F = \sqrt{\sum|a_{ij}|^2}$, the Frobenius norm of LA-1 Exercise C3.

**Definition 3.3 (Induced norm, orthogonality).** $\|v\| = \sqrt{\langle v, v\rangle}$. Vectors $u, v$ are *orthogonal*, written $u \perp v$, if $\langle u, v\rangle = 0$. A set is *orthonormal* if its members are pairwise orthogonal and each has norm 1.

**Theorem 3.4 (Cauchy–Schwarz).** In any inner product space,
$$|\langle u, v\rangle| \le \|u\|\,\|v\|, \tag{3.4}$$
with equality if and only if one of $u, v$ is a scalar multiple of the other.

▶ **Try it first.** Reread the LA-1 proof of Theorem 3.11. Replace every dot with $\langle\cdot,\cdot\rangle$. Does any step use a coordinate, a geometric fact, or anything beyond the three axioms?

**Derivation (real case).** This is the LA-1 proof, word for word, with the dot product renamed.

1. If $v = 0$ both sides are $0$ and we're done. Assume $v \ne 0$.
2. For every $t\in\mathbb{R}$, by axiom 3 and then axioms 1–2:
$$0 \le \|u - tv\|^2 = \langle u - tv, u - tv\rangle = \|u\|^2 - 2t\langle u, v\rangle + t^2\|v\|^2.$$
3. The right side is a quadratic $q(t) = \|v\|^2 t^2 - 2\langle u, v\rangle t + \|u\|^2$ with positive leading coefficient that is never negative, so it has at most one real root, so its discriminant is $\le 0$:
$$4\langle u, v\rangle^2 - 4\|v\|^2\|u\|^2 \le 0.$$
4. Rearranging: $\langle u, v\rangle^2 \le \|u\|^2\|v\|^2$. Take square roots.
5. Equality holds iff the discriminant is $0$ iff $q$ has a real root $t_0$ iff $\|u - t_0 v\| = 0$ iff $u = t_0 v$. ∎

**Complex case.** A few extra lines. With $t$ real, expanding $\|u - tv\|^2$ gives $\|u\|^2 - 2t\operatorname{Re}\langle u, v\rangle + t^2\|v\|^2$ (the cross terms are $-t\langle u,v\rangle - t\langle v, u\rangle = -2t\operatorname{Re}\langle u,v\rangle$ by axiom 1), and steps 3–4 yield $|\operatorname{Re}\langle u, v\rangle| \le \|u\|\|v\|$. Now write $\langle u, v\rangle = |\langle u,v\rangle|e^{i\theta}$ and apply this to $u' = e^{i\theta}u$: by conjugate-linearity in the first slot, $\langle u', v\rangle = e^{-i\theta}\langle u, v\rangle = |\langle u, v\rangle|$, which is real, and $\|u'\| = \|u\|$. So $|\langle u,v\rangle| = \operatorname{Re}\langle u', v\rangle \le \|u'\|\|v\| = \|u\|\|v\|$. *Equality:* if $|\langle u,v\rangle| = \|u\|\|v\|$ (and $v \ne 0$), then $\operatorname{Re}\langle u', v\rangle = \|u'\|\|v\|$, so the quadratic $\|u' - tv\|^2 = \|u'\|^2 - 2t\operatorname{Re}\langle u', v\rangle + t^2\|v\|^2$ has discriminant $0$ and step 5 gives a real $t_0$ with $u' = t_0 v$, i.e. $u = t_0e^{-i\theta}v$. Conversely, if $u = cv$ then $|\langle cv, v\rangle| = |\bar c|\,\|v\|^2 = \|u\|\,\|v\|$. ∎

■ **What just happened.** Nothing changed. The LA-1 proof never used coordinates or geometry, only the three axioms, so it holds for functions, matrices, and anything else with an inner product. In particular $\left|\int_a^b \bar f g\right| \le \sqrt{\int|f|^2}\sqrt{\int|g|^2}$, an inequality you will use in F2 and in QM-4's uncertainty principle.

**Theorem 3.5 (Triangle inequality).** $\|u + v\| \le \|u\| + \|v\|$.

*Proof.* $\|u+v\|^2 = \|u\|^2 + 2\operatorname{Re}\langle u, v\rangle + \|v\|^2 \le \|u\|^2 + 2|\langle u,v\rangle| + \|v\|^2 \le (\|u\| + \|v\|)^2$ by Cauchy–Schwarz. ∎

So the induced norm is a norm in the sense of LA-1 §3.6. It also satisfies the parallelogram law $\|u+v\|^2 + \|u-v\|^2 = 2\|u\|^2 + 2\|v\|^2$ (LA-1 Exercise B6: expand both sides). The converse, that a norm satisfying the parallelogram law comes from an inner product, is the Jordan–von Neumann theorem (1935); Exercise C2 does the easy half, recovering the inner product from the norm.

### 3.2 Orthonormal bases, Parseval, Gram–Schmidt

**Theorem 3.6 (Coordinates are inner products; Parseval).** Let $q_1, \dots, q_n$ be an orthonormal basis of $V$. For every $v \in V$,
$$v = \sum_{i=1}^n \langle q_i, v\rangle\, q_i, \tag{3.5}$$
$$\|v\|^2 = \sum_{i=1}^n |\langle q_i, v\rangle|^2. \tag{3.6}$$

▶ **Try it first.** Write $v = \sum_j c_j q_j$ (a basis, so the $c_j$ exist). Take the inner product of both sides with $q_i$, keeping $q_i$ in the first slot.

**Derivation.**
1. $v = \sum_j c_j q_j$ for some scalars $c_j$ since the $q_j$ form a basis.
2. $\langle q_i, v\rangle = \langle q_i, \sum_j c_j q_j\rangle = \sum_j c_j\langle q_i, q_j\rangle = \sum_j c_j\delta_{ij} = c_i$, using linearity in the second slot. So $c_i = \langle q_i, v\rangle$, which is (3.5).
3. $\|v\|^2 = \langle\sum_i c_i q_i, \sum_j c_j q_j\rangle = \sum_{i,j}\bar c_i c_j\langle q_i, q_j\rangle = \sum_i |c_i|^2$, using conjugate-linearity in the first slot. That is (3.6). ∎

■ **What just happened.** In an orthonormal basis, a coordinate is a single inner product, and Pythagoras holds in $n$ dimensions. **Convention alert:** the coefficient is $\langle q_i, v\rangle$ with the basis vector in the *first* slot. In Axler's convention it would be written $\langle v, q_i\rangle$; the two differ by a complex conjugate. In Dirac notation (3.5) reads $|v\rangle = \sum_i |q_i\rangle\langle q_i | v\rangle$, which is the reason for our choice. Fourier's coefficient formula from §1 is (3.5) for the orthonormal basis $\{1/\sqrt{2\pi},\ \cos nx/\sqrt\pi,\ \sin nx/\sqrt\pi\}$ of (3.2), or equivalently $\langle q, f\rangle/\|q\|^2$ for the merely orthogonal $q = \sin nx$ (that is where the $1/\pi$ comes from), with the infinite-dimensional subtleties handled in F2; (3.6) there is Parseval's identity for Fourier series.

**Gram–Schmidt** works verbatim in any inner product space (LA-5 §3.8): given independent $v_1, v_2, \dots$, set
$$u_k = v_k - \sum_{j<k}\frac{\langle u_j, v_k\rangle}{\langle u_j, u_j\rangle}\,u_j, \qquad q_k = \frac{u_k}{\|u_k\|}. \tag{3.7}$$
Each subtraction removes the component of $v_k$ along an earlier $u_j$, by the projection formula of LA-1 Example 4.5 with the dot replaced by $\langle\cdot,\cdot\rangle$. Example 4.1 runs (3.7) on $1, x, x^2, x^3$ with the inner product $\int_{-1}^1 fg\,dx$ and produces the Legendre polynomials.

### 3.3 The adjoint

**Definition 3.7 (Adjoint).** Let $A : V\to V$ be linear on a finite-dimensional inner product space. The adjoint $A^*$ is the linear map with
$$\langle Au, v\rangle = \langle u, A^* v\rangle \qquad\text{for all } u, v \in V. \tag{3.8}$$

Two things need checking: that such a map exists, and that it is unique. Uniqueness: if $\langle u, Bv\rangle = \langle u, Cv\rangle$ for all $u, v$, then $\langle u, (B - C)v\rangle = 0$ for all $u$; take $u = (B-C)v$ to get $\|(B-C)v\|^2 = 0$, so $Bv = Cv$ for every $v$. Existence we get by writing it down:

**Theorem 3.8.** With the standard inner product on $\mathbb{C}^n$, the adjoint of a matrix $A$ is $A^\dagger$. On $\mathbb{R}^n$ it is $A^\top$.

▶ **Try it first.** Write $\langle A\mathbf{u}, \mathbf{v}\rangle$ as a matrix product using (3.1) and move $A$ to the other side using LA-1 Theorem 3.8, $(AB)^\top = B^\top A^\top$, in its conjugate-transpose form $(AB)^\dagger = B^\dagger A^\dagger$.

**Derivation.**
1. $(A\mathbf{u})^\dagger = \mathbf{u}^\dagger A^\dagger$: conjugate both sides of $(A\mathbf{u})^\top = \mathbf{u}^\top A^\top$.
2. $\langle A\mathbf{u}, \mathbf{v}\rangle = (A\mathbf{u})^\dagger\mathbf{v} = \mathbf{u}^\dagger A^\dagger\mathbf{v} = \mathbf{u}^\dagger(A^\dagger\mathbf{v}) = \langle\mathbf{u}, A^\dagger\mathbf{v}\rangle$.
3. So $A^\dagger$ satisfies (3.8); by uniqueness $A^* = A^\dagger$. For real matrices $A^\dagger = A^\top$. ∎

*Any finite-dimensional $V$.* Pick an orthonormal basis $q_1, \dots, q_n$ (Gram–Schmidt, §3.2) and write $\mathbf{c}(v) = (\langle q_i, v\rangle)_i$ for the coordinates of $v$. Step 3 of Theorem 3.6's derivation, run with two different vectors, gives $\langle u, v\rangle = \mathbf{c}(u)^\dagger\mathbf{c}(v)$: in these coordinates every inner product *is* the standard one (3.1). So if $[A]$ is the matrix of $A$ in this basis ($\mathbf{c}(Au) = [A]\,\mathbf{c}(u)$), then $\langle Au, v\rangle = ([A]\mathbf{c}(u))^\dagger\mathbf{c}(v) = \mathbf{c}(u)^\dagger\big([A]^\dagger\mathbf{c}(v)\big)$, so the map whose matrix is $[A]^\dagger$ satisfies (3.8), and the adjoint exists on every finite-dimensional inner product space. Theorems 3.9 and 3.10 below, and Exercise B1, use it for inner products other than (3.1).

■ **What just happened.** The transpose, which LA-1 introduced as "rows become columns," is really the operator that lets you move a matrix from one side of an inner product to the other. That is the right way to think about it from now on. Properties that follow directly from (3.8): $(A^*)^* = A$, $(AB)^* = B^*A^*$, $(cA)^* = \bar c A^*$. Exercise B1 shows that with a *different* inner product $\langle\mathbf{u},\mathbf{v}\rangle_W = \mathbf{u}^\top W\mathbf{v}$ the adjoint is $W^{-1}A^\top W$, not $A^\top$; the adjoint depends on the inner product.

**Definition (Self-adjoint).** $A$ is self-adjoint if $A^* = A$. For matrices with the standard inner product this means $A^\top = A$ (real: **symmetric**) or $A^\dagger = A$ (complex: **Hermitian**). A Hermitian matrix has real diagonal entries ($a_{ii} = \bar a_{ii}$) and conjugate-mirrored off-diagonal entries ($a_{ji} = \bar a_{ij}$). Note: a complex matrix with $A^\top = A$ is **not** self-adjoint in general; see §8 and Exercise A5.

### 3.4 Hermitian operators: real eigenvalues, orthogonal eigenvectors

**Theorem 3.9.** Every eigenvalue of a self-adjoint operator is real.

▶ **Try it first.** Let $Av = \lambda v$ with $v\ne 0$. Compute $\langle v, Av\rangle$ two ways: once by pulling $\lambda$ out of the second slot, once by moving $A$ to the first slot and pulling $\lambda$ out of there.

**Derivation.**
1. $\langle v, Av\rangle = \langle v, \lambda v\rangle = \lambda\langle v, v\rangle = \lambda\|v\|^2$ (linear in the second slot).
2. $\langle v, Av\rangle = \langle A^* v, v\rangle = \langle Av, v\rangle = \langle\lambda v, v\rangle = \bar\lambda\|v\|^2$ (self-adjoint, then **conjugate**-linear in the first slot; this is the convention from the box doing real work).
3. So $\lambda\|v\|^2 = \bar\lambda\|v\|^2$ and $\|v\|^2 \ne 0$, hence $\lambda = \bar\lambda$: $\lambda$ is real. ∎

**Theorem 3.10.** Eigenvectors of a self-adjoint operator belonging to distinct eigenvalues are orthogonal.

**Derivation.** Let $Av = \lambda v$, $Aw = \mu w$, $\lambda\ne\mu$.
1. $\lambda\langle v, w\rangle = \bar\lambda\langle v, w\rangle = \langle\lambda v, w\rangle = \langle Av, w\rangle$, using that $\lambda$ is real (Theorem 3.9) and conjugate-linearity in the first slot.
2. $\langle Av, w\rangle = \langle v, Aw\rangle = \langle v, \mu w\rangle = \mu\langle v, w\rangle$.
3. So $(\lambda - \mu)\langle v, w\rangle = 0$ and $\lambda\ne\mu$, hence $\langle v, w\rangle = 0$. ∎

■ **What just happened.** Both proofs are two lines because the adjoint lets $A$ hop across the inner product, and self-adjointness means it lands as itself. One more fact for the real case: if $A$ is real symmetric and $\lambda$ is an eigenvalue, then $\lambda$ is real, so $A - \lambda I$ is a *real* singular matrix, so its null space contains a *real* nonzero vector. Real symmetric matrices have real eigenvectors. LA-7's complex eigenvectors for rotation matrices cannot happen here.

### 3.5 The spectral theorem

**Theorem 3.11 (Spectral theorem, real symmetric).** Let $A$ be a real symmetric $n\times n$ matrix. Then there is an orthogonal matrix $Q$ ($Q^\top Q = I$) and a real diagonal matrix $\Lambda = \operatorname{diag}(\lambda_1, \dots, \lambda_n)$ with
$$A = Q\Lambda Q^\top. \tag{3.9}$$
Equivalently: the columns $\mathbf{q}_1, \dots, \mathbf{q}_n$ of $Q$ form an orthonormal basis of $\mathbb{R}^n$ consisting of eigenvectors of $A$, $A\mathbf{q}_i = \lambda_i\mathbf{q}_i$.

The proof has two parts. Part I finds *one* eigenvector by maximizing a function on the sphere. Part II shows that the directions perpendicular to it form a smaller problem of the same type, and inducts. Part II is where people stop understanding; take it slowly.

**Part I: a largest eigenvalue exists.** Define the **Rayleigh quotient**
$$R(\mathbf{x}) = \frac{\mathbf{x}^\top A\mathbf{x}}{\mathbf{x}^\top\mathbf{x}}, \qquad \mathbf{x}\ne\mathbf{0}$$
(it gets its equation number, (3.13), when §3.9 studies it in its own right). Since $R(c\mathbf{x}) = R(\mathbf{x})$ for $c\ne 0$, it is enough to study $g(\mathbf{x}) = \mathbf{x}^\top A\mathbf{x}$ on the unit sphere $S = \{\mathbf{x} : \|\mathbf{x}\| = 1\}$.

▶ **Try it first.** $g$ is continuous and $S$ is closed and bounded, so $g$ has a maximum on $S$ at some $\mathbf{q}$. Set $\lambda = g(\mathbf{q})$. You want to show $A\mathbf{q} = \lambda\mathbf{q}$. Take any $\mathbf{w}\perp\mathbf{q}$ and look at $R(\mathbf{q} + t\mathbf{w})$ as a function of $t$ near $0$. What must its derivative at $t = 0$ be?

**Derivation.**
1. (The one fact imported from analysis.) A continuous function on a closed bounded subset of $\mathbb{R}^n$ attains its maximum (extreme value theorem). So there is $\mathbf{q}\in S$ with $g(\mathbf{q}) \ge g(\mathbf{x})$ for all $\mathbf{x}\in S$; equivalently $R(\mathbf{q}) \ge R(\mathbf{x})$ for all $\mathbf{x}\ne\mathbf{0}$. Let $\lambda = R(\mathbf{q}) = \mathbf{q}^\top A\mathbf{q}$.
2. Fix any $\mathbf{w}$ with $\mathbf{q}^\top\mathbf{w} = 0$ and define $h(t) = R(\mathbf{q} + t\mathbf{w})$ for real $t$. Since $h(t)\le h(0)$ for all $t$, $h$ has a maximum at $t = 0$, so $h'(0) = 0$.
3. Compute $h$. Denominator: $(\mathbf{q}+t\mathbf{w})^\top(\mathbf{q}+t\mathbf{w}) = 1 + 2t\,\mathbf{q}^\top\mathbf{w} + t^2\|\mathbf{w}\|^2 = 1 + t^2\|\mathbf{w}\|^2$. Numerator: $(\mathbf{q}+t\mathbf{w})^\top A(\mathbf{q}+t\mathbf{w}) = \mathbf{q}^\top A\mathbf{q} + t(\mathbf{q}^\top A\mathbf{w} + \mathbf{w}^\top A\mathbf{q}) + t^2\mathbf{w}^\top A\mathbf{w}$. Now use symmetry: $\mathbf{q}^\top A\mathbf{w}$ is a number, equal to its own transpose $\mathbf{w}^\top A^\top\mathbf{q} = \mathbf{w}^\top A\mathbf{q}$. So the numerator is $\lambda + 2t\,\mathbf{w}^\top A\mathbf{q} + t^2\mathbf{w}^\top A\mathbf{w}$.
4. Differentiate the quotient at $t = 0$. The denominator equals $1$ there and its derivative $2t\|\mathbf{w}\|^2$ vanishes there, so $h'(0) = $ (derivative of numerator at $0$) $= 2\,\mathbf{w}^\top A\mathbf{q}$.
5. Hence $\mathbf{w}^\top(A\mathbf{q}) = 0$ for every $\mathbf{w}\perp\mathbf{q}$. The vector $A\mathbf{q}$ is orthogonal to the entire hyperplane $\mathbf{q}^\perp$, so it lies in $(\mathbf{q}^\perp)^\perp = \operatorname{span}\{\mathbf{q}\}$ (LA-5: $(W^\perp)^\perp = W$ in finite dimensions). So $A\mathbf{q} = \mu\mathbf{q}$ for some $\mu$, and $\mu = \mathbf{q}^\top A\mathbf{q} = \lambda$.
6. Moreover $\lambda$ is the *largest* eigenvalue: if $A\mathbf{v} = \mu\mathbf{v}$ with $\|\mathbf{v}\| = 1$ then $\mu = R(\mathbf{v}) \le R(\mathbf{q}) = \lambda$. ∎ (Part I)

■ **What just happened.** The top eigenvector is the direction that maximizes $\mathbf{x}^\top A\mathbf{x}$ on the sphere. Step 5 is the whole idea: at a maximum on the sphere, the gradient $2A\mathbf{q}$ of $g$ has no component tangent to the sphere, so it points along the normal $\mathbf{q}$. If you did Lagrange multipliers in multivariable calculus, steps 2–5 are $\nabla(\mathbf{x}^\top A\mathbf{x}) = \mu\nabla(\mathbf{x}^\top\mathbf{x})$, i.e. $2A\mathbf{x} = 2\mu\mathbf{x}$, done by hand. This "eigenvector = constrained optimum" picture is the one LA-9 uses to derive the SVD and ML-6 uses to derive PCA. (An alternative existence proof: the characteristic polynomial has a complex root by the fundamental theorem of algebra (LA-7), that root is real by Theorem 3.9, and its eigenvector is real by the remark after Theorem 3.10. Shorter, but it tells you nothing about *which* eigenvalue you found.)

**Part II: the orthogonal complement is invariant, and induction.** Let $\mathbf{q}_1 = \mathbf{q}$ and $\lambda_1 = \lambda$ from Part I.

▶ **Try it first.** Let $W = \mathbf{q}_1^\perp$. Show that if $\mathbf{x}\in W$ then $A\mathbf{x}\in W$. You need one line, and it must use $A^\top = A$. Then ask: where would this fail for the shear $\begin{pmatrix}1&1\\0&1\end{pmatrix}$?

**Derivation.**
7. (Invariance.) Let $\mathbf{x}\in W$, i.e. $\mathbf{q}_1^\top\mathbf{x} = 0$. Then
$$\mathbf{q}_1^\top(A\mathbf{x}) = (A^\top\mathbf{q}_1)^\top\mathbf{x} = (A\mathbf{q}_1)^\top\mathbf{x} = \lambda_1\mathbf{q}_1^\top\mathbf{x} = 0.$$
So $A\mathbf{x}\in W$. **This is the step.** The first equality is LA-1 Theorem 3.8; the second is symmetry; the third is that $\mathbf{q}_1$ is an eigenvector. Without symmetry, $A^\top\mathbf{q}_1$ is some unrelated vector and the argument dies. For the shear, $\mathbf{q}_1 = \mathbf{e}_1$, $W = \operatorname{span}\{\mathbf{e}_2\}$, and $A\mathbf{e}_2 = (1, 1)^\top\notin W$. The complement leaks.
8. (Write it as a block matrix.) By Gram–Schmidt (LA-5) choose an orthonormal basis $\mathbf{v}_2, \dots, \mathbf{v}_n$ of $W$ and put them in the $n\times(n-1)$ matrix $V$. Then $Q_1 = \begin{pmatrix}\mathbf{q}_1 & V\end{pmatrix}$ is $n\times n$ orthogonal, $V^\top V = I_{n-1}$, and $V^\top\mathbf{q}_1 = \mathbf{0}$. Compute
$$Q_1^\top A Q_1 = \begin{pmatrix}\mathbf{q}_1^\top A\mathbf{q}_1 & \mathbf{q}_1^\top A V\\ V^\top A\mathbf{q}_1 & V^\top A V\end{pmatrix} = \begin{pmatrix}\lambda_1 & \mathbf{0}^\top\\ \mathbf{0} & B\end{pmatrix}, \qquad B = V^\top A V.$$
The entries: $\mathbf{q}_1^\top A\mathbf{q}_1 = \lambda_1\mathbf{q}_1^\top\mathbf{q}_1 = \lambda_1$. Bottom-left: $V^\top A\mathbf{q}_1 = \lambda_1 V^\top\mathbf{q}_1 = \mathbf{0}$. Top-right: $\mathbf{q}_1^\top A V = (A^\top\mathbf{q}_1)^\top V = (A\mathbf{q}_1)^\top V = \lambda_1\mathbf{q}_1^\top V = \mathbf{0}^\top$; this is step 7 again, in matrix form. The zero blocks *are* the statement that $A$ maps $W$ into $W$ and $\operatorname{span}\{\mathbf{q}_1\}$ into itself.
9. ($B$ is a smaller instance of the problem.) $B$ is $(n-1)\times(n-1)$ and symmetric: $B^\top = V^\top A^\top V = V^\top A V = B$. Concretely, $B$ is "$A$ restricted to $W$," written in the coordinates $V$: if $\mathbf{x} = V\mathbf{z}\in W$ then $A\mathbf{x}\in W$ by step 7, and its coordinates in the basis $V$ are $V^\top A V\mathbf{z} = B\mathbf{z}$.
10. (Induction.) Induct on $n$. For $n = 1$, $A = (a)$, $Q = (1)$, $\Lambda = (a)$. For $n > 1$, apply the theorem to $B$: $B = Q_2\Lambda_2 Q_2^\top$ with $Q_2$ orthogonal, $(n-1)\times(n-1)$. Set
$$Q = Q_1\begin{pmatrix}1 & \mathbf{0}^\top\\ \mathbf{0} & Q_2\end{pmatrix}.$$
A product of orthogonal matrices is orthogonal ($(Q_1 D)^\top Q_1 D = D^\top Q_1^\top Q_1 D = D^\top D = I$), so $Q$ is orthogonal, and
$$Q^\top A Q = \begin{pmatrix}1 & \mathbf{0}^\top\\ \mathbf{0} & Q_2^\top\end{pmatrix}\begin{pmatrix}\lambda_1 & \mathbf{0}^\top\\ \mathbf{0} & B\end{pmatrix}\begin{pmatrix}1 & \mathbf{0}^\top\\ \mathbf{0} & Q_2\end{pmatrix} = \begin{pmatrix}\lambda_1 & \mathbf{0}^\top\\ \mathbf{0} & Q_2^\top B Q_2\end{pmatrix} = \begin{pmatrix}\lambda_1 & \mathbf{0}^\top\\ \mathbf{0} & \Lambda_2\end{pmatrix} = \Lambda.$$
Multiply on the left by $Q$ and on the right by $Q^\top$: $A = Q\Lambda Q^\top$. ∎

■ **What just happened.** Peel off one eigenvector. What is left (the action of $A$ on the perpendicular directions) is again symmetric, one dimension smaller, so peel again. Symmetry was used three times: in Part I (step 3, to get $\mathbf{w}^\top A\mathbf{q}$ as the derivative), in steps 7–8 (to make the off-diagonal blocks vanish), and in step 9 ($B$ is again symmetric, which is the only reason the induction hypothesis applies to it). The columns of $Q$ other than $\mathbf{q}_1$ are $V\mathbf{z}$ for the eigenvectors $\mathbf{z}$ of $B$, and Exercise B8 asks you to verify directly that each $V\mathbf{z}$ is an eigenvector of $A$ lying in $W$. The project's Jacobi algorithm (§7) is this proof turned into arithmetic: it applies rotations until the off-diagonal blocks are zero.

**Theorem 3.12 (Spectral theorem, Hermitian).** If $A$ is Hermitian ($A^\dagger = A$), there is a unitary $U$ ($U^\dagger U = I$) and a *real* diagonal $\Lambda$ with $A = U\Lambda U^\dagger$.

*Proof.* Identical to Theorem 3.11 with $\dagger$ in place of $\top$ throughout. For existence of the first eigenvalue use the fundamental theorem of algebra (the characteristic polynomial has a complex root), then Theorem 3.9 to see it is real. Step 7 becomes $\mathbf{q}_1^\dagger(A\mathbf{x}) = (A^\dagger\mathbf{q}_1)^\dagger\mathbf{x} = (A\mathbf{q}_1)^\dagger\mathbf{x} = \bar\lambda_1\mathbf{q}_1^\dagger\mathbf{x} = 0$, and $\lambda_1$ is real anyway. The block computation in step 8 goes through with $V^\dagger$. ∎

### 3.6 The spectral decomposition: a sum of projections

**Theorem 3.13.** Let $A = Q\Lambda Q^\top$ as in (3.9), with columns $\mathbf{q}_i$. Then
$$A = \sum_{i=1}^n \lambda_i\,\mathbf{q}_i\mathbf{q}_i^\top = \sum_{i=1}^n\lambda_i P_i, \qquad P_i = \mathbf{q}_i\mathbf{q}_i^\top, \tag{3.10}$$
where each $P_i$ is the orthogonal projection onto the line through $\mathbf{q}_i$: $P_i^2 = P_i$, $P_i^\top = P_i$, $P_iP_j = 0$ for $i\ne j$, and $\sum_i P_i = I$. Consequently for any polynomial $f$, and by definition for any function $f$ defined on the eigenvalues,
$$f(A) = \sum_{i=1}^n f(\lambda_i)\,P_i = Q\,f(\Lambda)\,Q^\top. \tag{3.11}$$

▶ **Try it first.** Write $Q\Lambda Q^\top$ as $(Q\Lambda)(Q^\top)$ and use LA-1's column picture: the $k$-th column of $Q\Lambda$ is $\lambda_k\mathbf{q}_k$. A product $CD$ equals $\sum_k(\text{column } k \text{ of } C)(\text{row } k \text{ of } D)$. Then compute $P_iP_j$.

**Derivation.**
1. $(CD)_{ab} = \sum_k c_{ak}d_{kb}$, so $CD = \sum_k \mathbf{c}_k\mathbf{d}_k^\top$ where $\mathbf{c}_k$ is column $k$ of $C$ and $\mathbf{d}_k^\top$ is row $k$ of $D$ (check the $(a,b)$ entry of each outer product: $c_{ak}d_{kb}$). Apply with $C = Q\Lambda$ (columns $\lambda_k\mathbf{q}_k$) and $D = Q^\top$ (rows $\mathbf{q}_k^\top$): $A = \sum_k\lambda_k\mathbf{q}_k\mathbf{q}_k^\top$.
2. $P_iP_j = \mathbf{q}_i(\mathbf{q}_i^\top\mathbf{q}_j)\mathbf{q}_j^\top = \delta_{ij}\mathbf{q}_i\mathbf{q}_j^\top$, so $P_i^2 = P_i$ and $P_iP_j = 0$ for $i\ne j$. $P_i^\top = (\mathbf{q}_i^\top)^\top\mathbf{q}_i^\top = P_i$. And $\sum_i P_i = QQ^\top = I$ (step 1 with $\Lambda = I$; $QQ^\top = I$ because a square matrix with $Q^\top Q = I$ has $Q^{-1} = Q^\top$).
3. $A^2 = \sum_{i,j}\lambda_i\lambda_jP_iP_j = \sum_i\lambda_i^2P_i$ by step 2, and inductively $A^k = \sum_i\lambda_i^kP_i$ (with $A^0 = I = \sum P_i$). For a polynomial $f(t) = \sum_k c_kt^k$, $f(A) = \sum_k c_kA^k = \sum_i\left(\sum_kc_k\lambda_i^k\right)P_i = \sum_i f(\lambda_i)P_i$. In matrix form that is $Q f(\Lambda) Q^\top$ with $f(\Lambda) = \operatorname{diag}(f(\lambda_i))$. ∎

■ **What just happened.** A symmetric matrix is a weighted sum of projections onto perpendicular lines: project $\mathbf{x}$ onto each eigen-axis, scale by $\lambda_i$, add up. Any function you can apply to numbers you can now apply to $A$ by applying it to the eigenvalues: $A^{-1} = \sum\lambda_i^{-1}P_i$ (if no $\lambda_i = 0$), $\sqrt{A} = \sum\sqrt{\lambda_i}P_i$ (if all $\lambda_i\ge 0$), $e^{tA} = \sum e^{t\lambda_i}P_i$, which agrees with LA-7's power series since (3.11) holds for every partial sum. Example 4.2 computes $\sqrt A$ and $A^{-1}$ this way. In quantum mechanics (3.10) is written $A = \sum\lambda_i|q_i\rangle\langle q_i|$ and $\sum_i|q_i\rangle\langle q_i| = I$ is called *completeness*; see §3.12.

**Remark ($f(A)$ does not depend on $Q$; the square root is unique).** When an eigenvalue repeats, $Q$ is not unique (Example 4.2 picks one orthonormal basis of a plane out of infinitely many), so (3.11) needs a check.
1. Fix a distinct eigenvalue $\mu$ of $A$ and set $E_\mu = \sum_{i:\,\lambda_i = \mu}P_i$. The $\mathbf{q}_i$ with $\lambda_i = \mu$ are orthonormal and lie in the eigenspace $N(A - \mu I)$. They span it: if $A\mathbf{x} = \mu\mathbf{x}$ and $\mathbf{y} = Q^\top\mathbf{x}$, then $\Lambda\mathbf{y} = Q^\top A\mathbf{x} = \mu\mathbf{y}$, so $y_i = 0$ whenever $\lambda_i\ne\mu$ and $\mathbf{x} = Q\mathbf{y} = \sum_{\lambda_i = \mu}y_i\mathbf{q}_i$. So $E_\mu$ is the orthogonal projection onto $N(A - \mu I)$ (LA-5), and that projection is determined by the subspace, not by the basis chosen inside it.
2. Grouping equal eigenvalues in (3.11) gives $f(A) = \sum_\mu f(\mu)E_\mu$, a formula that mentions only $A$. So $f(A)$ is well defined. For $\mathbf{x}\in N(A - \mu I)$, $E_\mu\mathbf{x} = \mathbf{x}$ and $E_\nu\mathbf{x} = \mathbf{0}$ for $\nu\ne\mu$ (Theorem 3.10), so $f(A)\mathbf{x} = f(\mu)\mathbf{x}$.
3. Let $S$ be symmetric with $S^2 = A$ and all eigenvalues of $S$ nonnegative. (Then $A = S^\top S$ has $\mathbf{x}^\top A\mathbf{x} = \|S\mathbf{x}\|^2\ge 0$, so every $\lambda_i = \mathbf{q}_i^\top A\mathbf{q}_i\ge 0$ and $\sqrt A$ is defined.) By Theorem 3.11, $S$ has an orthonormal eigenbasis $\mathbf{r}_1, \dots, \mathbf{r}_n$ with $S\mathbf{r}_j = s_j\mathbf{r}_j$, $s_j\ge 0$. Then $A\mathbf{r}_j = S^2\mathbf{r}_j = s_j^2\mathbf{r}_j$, so by step 2, $\sqrt A\,\mathbf{r}_j = \sqrt{s_j^2}\,\mathbf{r}_j = s_j\mathbf{r}_j = S\mathbf{r}_j$. Two matrices that agree on a basis are equal: $S = \sqrt A = Q\Lambda^{1/2}Q^\top$. ∎

In the language of §3.7: a positive semidefinite matrix has exactly one positive semidefinite square root. LA-9's polar decomposition relies on this.

### 3.7 Quadratic forms and principal axes

**Definition.** A quadratic form on $\mathbb{R}^n$ is a function $q(\mathbf{x}) = \mathbf{x}^\top A\mathbf{x} = \sum_{i,j}a_{ij}x_ix_j$ for some $n\times n$ matrix $A$.

**Theorem 3.14 (Principal axes).** (a) $\mathbf{x}^\top A\mathbf{x} = \mathbf{x}^\top S\mathbf{x}$ where $S = \tfrac12(A + A^\top)$ is the symmetric part of $A$; so assume $A$ symmetric. (b) With $A = Q\Lambda Q^\top$ and $\mathbf{y} = Q^\top\mathbf{x}$ (the coordinates of $\mathbf{x}$ in the eigenbasis, $y_i = \mathbf{q}_i^\top\mathbf{x}$),
$$\mathbf{x}^\top A\mathbf{x} = \sum_{i=1}^n\lambda_i y_i^2. \tag{3.12}$$
(c) The level set $L = \{\mathbf{x} : \mathbf{x}^\top A\mathbf{x} = 1\}$ is:
- **empty** if every $\lambda_i\le 0$;
- an **ellipsoid** if every $\lambda_i > 0$, with axes along the lines through the $\mathbf{q}_i$ and semi-axis length $1/\sqrt{\lambda_i}$ along $\mathbf{q}_i$;
- a **hyperboloid** (a hyperbola in $\mathbb{R}^2$) if no $\lambda_i$ is $0$ and both signs occur;
- a **cylinder** if some $\lambda_i = 0$ and some $\lambda_j > 0$: $L$ is unbounded along every $\mathbf{q}_i$ with $\lambda_i = 0$, and is the level set of the remaining coordinates swept along those directions (a pair of parallel lines in $\mathbb{R}^2$; an elliptic or hyperbolic cylinder in $\mathbb{R}^3$).

▶ **Try it first.** For (a), note that $\mathbf{x}^\top A\mathbf{x}$ is a $1\times 1$ matrix and so equals its own transpose. For (c), set all $y_j = 0$ except $y_i$ in (3.12).

**Derivation.**
1. $\mathbf{x}^\top A\mathbf{x} = (\mathbf{x}^\top A\mathbf{x})^\top = \mathbf{x}^\top A^\top\mathbf{x}$. Average the two expressions: $\mathbf{x}^\top A\mathbf{x} = \mathbf{x}^\top\tfrac12(A+A^\top)\mathbf{x}$. The antisymmetric part $\tfrac12(A - A^\top)$ contributes nothing. This is why $5x^2 + 4xy + 2y^2$ in §2 is encoded with $2$ and $2$ off the diagonal, not $4$ and $0$.
2. $\mathbf{x}^\top A\mathbf{x} = \mathbf{x}^\top Q\Lambda Q^\top\mathbf{x} = (Q^\top\mathbf{x})^\top\Lambda(Q^\top\mathbf{x}) = \mathbf{y}^\top\Lambda\mathbf{y} = \sum_i\lambda_iy_i^2$. The map $\mathbf{x}\mapsto\mathbf{y} = Q^\top\mathbf{x}$ preserves lengths and angles ($\|\mathbf{y}\|^2 = \mathbf{x}^\top QQ^\top\mathbf{x} = \|\mathbf{x}\|^2$; §3.10 says more), so the level set in $\mathbf{y}$-coordinates is the level set in $\mathbf{x}$-coordinates, rigidly rotated. The $\mathbf{y}$-axes correspond to the directions $\mathbf{q}_i$ in $\mathbf{x}$-space.
3. On the $i$-th axis ($y_j = 0$ for $j\ne i$), the equation $\sum_j\lambda_jy_j^2 = 1$ reads $\lambda_iy_i^2 = 1$, so $y_i = \pm1/\sqrt{\lambda_i}$ if $\lambda_i > 0$. Those are the points where the surface crosses the axis: semi-axis length $1/\sqrt{\lambda_i}$. **Large eigenvalue, short axis.** If $\lambda_i < 0$ the axis never meets the surface (hyperbolic direction).
4. The other cases. If every $\lambda_i\le 0$, then $\sum_j\lambda_jy_j^2\le 0 < 1$ for every $\mathbf{y}$, so $L$ is empty (Exercise A4(d)). If $\lambda_i = 0$, the coordinate $y_i$ does not appear in the equation, so whenever $\mathbf{y}$ solves it so does $\mathbf{y} + t\mathbf{e}_i$ for every $t$: $L$ is a union of lines parallel to $\mathbf{q}_i$, a cylinder over the level set of the other coordinates. That smaller level set is nonempty exactly when some $\lambda_j > 0$ (take $y_j = 1/\sqrt{\lambda_j}$ and the rest $0$). ∎

■ **What just happened.** Rotating to the eigenbasis removes every cross term, and then the shape of the level set is read off from the signs of the $\lambda_i$: positive ones bound it, negative ones open it up, zero ones make it a cylinder, and with no positive one there is nothing to draw. The cross term $4xy$ in §2 was only a sign that the axes were tilted.

**Definition (Definiteness).** A symmetric $A$ is **positive definite** (PD) if $\mathbf{x}^\top A\mathbf{x} > 0$ for all $\mathbf{x}\ne\mathbf{0}$; **positive semidefinite** (PSD) if $\mathbf{x}^\top A\mathbf{x}\ge 0$ for all $\mathbf{x}$; **negative (semi)definite** if $-A$ is positive (semi)definite; **indefinite** if $\mathbf{x}^\top A\mathbf{x}$ takes both signs. By (3.12) these are statements about the signs of the eigenvalues; §3.8 makes that precise. Written "$A\succ 0$" (PD) and "$A\succeq 0$" (PSD). For a complex **Hermitian** $A$ the same words mean the same inequalities for $\mathbf{x}^\dagger A\mathbf{x}$ with $\mathbf{x}\in\mathbb{C}^n$: Hermitian PD means $\mathbf{x}^\dagger A\mathbf{x} > 0$ for $\mathbf{x}\ne\mathbf{0}$, Hermitian PSD means $\mathbf{x}^\dagger A\mathbf{x}\ge 0$. (The number $\mathbf{x}^\dagger A\mathbf{x}$ is real because it equals its own conjugate: $\overline{\mathbf{x}^\dagger A\mathbf{x}} = (\mathbf{x}^\dagger A\mathbf{x})^\dagger = \mathbf{x}^\dagger A^\dagger\mathbf{x} = \mathbf{x}^\dagger A\mathbf{x}$.)

### 3.8 Positive definite matrices

**Theorem 3.15 (Tests for positive definiteness).** For a real symmetric $A$, the following are equivalent:
(i) every eigenvalue $\lambda_i > 0$;
(ii) $\mathbf{x}^\top A\mathbf{x} > 0$ for every $\mathbf{x}\ne\mathbf{0}$;
(iii) $A = B^\top B$ for some invertible matrix $B$;
(iv) every pivot in Gaussian elimination without row exchanges is $> 0$.

We prove (i)⇔(ii), (i)⇒(iii), (iii)⇒(ii). The equivalence with (iv) is proved in LA-10 via the $LDL^\top$ factorization ($D$ holds the pivots, and $\mathbf{x}^\top A\mathbf{x} = (L^\top\mathbf{x})^\top D(L^\top\mathbf{x})$ is a sum of pivots times squares); for $A = \begin{pmatrix}a&b\\b&c\end{pmatrix}$ the pivots are $a$ and $(ac - b^2)/a$, so (iv) says $a > 0$ and $\det A > 0$, which Exercise B3(d) and the two-variable rule after Theorem 3.19 use.

▶ **Try it first.** For (i)⇒(ii), use (3.12). For (ii)⇒(i), plug in $\mathbf{x} = \mathbf{q}_i$. For (i)⇒(iii), guess $B$ from $A = Q\Lambda Q^\top$ and $\Lambda = \Lambda^{1/2}\Lambda^{1/2}$.

**Derivation.**
1. (i)⇒(ii). Let $\mathbf{x}\ne\mathbf{0}$ and $\mathbf{y} = Q^\top\mathbf{x}$. Then $\mathbf{y}\ne\mathbf{0}$ ($Q^\top$ is invertible), so some $y_i\ne 0$, and $\mathbf{x}^\top A\mathbf{x} = \sum_i\lambda_iy_i^2 > 0$ since every term is $\ge 0$ and the $i$-th is $> 0$.
2. (ii)⇒(i). $\mathbf{q}_i^\top A\mathbf{q}_i = \lambda_i\mathbf{q}_i^\top\mathbf{q}_i = \lambda_i$. If (ii) holds, this is $> 0$.
3. (i)⇒(iii). Let $\Lambda^{1/2} = \operatorname{diag}(\sqrt{\lambda_i})$, real because $\lambda_i > 0$, and $B = \Lambda^{1/2}Q^\top$. Then $B^\top B = Q\Lambda^{1/2}\Lambda^{1/2}Q^\top = Q\Lambda Q^\top = A$, and $B$ is invertible as a product of invertible matrices ($B^{-1} = Q\Lambda^{-1/2}$). Another valid choice is $B = \sqrt{A} = Q\Lambda^{1/2}Q^\top$ from (3.11), which is symmetric. The choice that makes $B$ upper triangular is the **Cholesky factor** $A = R^\top R$; LA-10 constructs it from the pivots, and Example 4.6 shows one.
4. (iii)⇒(ii). $\mathbf{x}^\top A\mathbf{x} = \mathbf{x}^\top B^\top B\mathbf{x} = (B\mathbf{x})^\top(B\mathbf{x}) = \|B\mathbf{x}\|^2 \ge 0$, with equality iff $B\mathbf{x} = \mathbf{0}$ iff $\mathbf{x} = \mathbf{0}$ because $B$ is invertible. ∎

■ **What just happened.** Positive definite means "a sum of squares with positive weights, in rotated coordinates" (by (3.12)), or "$\|B\mathbf{x}\|^2$ for an invertible $B$" (by (iii)). Both say the same thing: $\mathbf{x}^\top A\mathbf{x}$ is a *squared length* in a distorted geometry. That is why a PD matrix defines an inner product $\langle\mathbf{u},\mathbf{v}\rangle_A = \mathbf{u}^\top A\mathbf{v}$ (Exercise B4), and why covariance matrices (F5) and kernel matrices (ML-4) are PD or PSD.

**Hermitian case.** For Hermitian $A$, (i)⇔(ii)⇔(iii) hold verbatim with $\dagger$ in place of $\top$, $\mathbf{x}\in\mathbb{C}^n$, and Theorem 3.12's unitary $U$ in place of $Q$: with $\mathbf{y} = U^\dagger\mathbf{x}$, $\mathbf{x}^\dagger A\mathbf{x} = \sum_i\lambda_i|y_i|^2$; $\mathbf{u}_i^\dagger A\mathbf{u}_i = \lambda_i$; $B = \Lambda^{1/2}U^\dagger$ gives $B^\dagger B = A$; and $\mathbf{x}^\dagger B^\dagger B\mathbf{x} = \|B\mathbf{x}\|^2$.

**Theorem 3.16 (Gram matrices are PSD).** For any real $m\times n$ matrix $A$, the matrix $A^\top A$ is symmetric positive semidefinite, and it is positive definite if and only if the columns of $A$ are linearly independent. More generally, for vectors $v_1, \dots, v_n$ in any inner product space, the **Gram matrix** $G_{ij} = \langle v_i, v_j\rangle$ is PSD (real case) or Hermitian PSD (complex case). A symmetric or Hermitian matrix is PSD if and only if all its eigenvalues are $\ge 0$.

▶ **Try it first.** Write $\mathbf{x}^\top A^\top A\mathbf{x}$ as the squared length of one vector. For the Gram matrix, expand $\|\sum_ix_iv_i\|^2$ with the axioms and compare.

**Derivation.**
1. $(A^\top A)^\top = A^\top A$. For any $\mathbf{x}$, $\mathbf{x}^\top A^\top A\mathbf{x} = \|A\mathbf{x}\|^2\ge 0$: PSD. Equality holds iff $A\mathbf{x} = \mathbf{0}$, and that has a nonzero solution iff the columns of $A$ are dependent (LA-3). So PD iff the null space of $A$ is $\{\mathbf{0}\}$ iff the columns are independent.
2. Gram matrix. $G$ is Hermitian: $G_{ji} = \langle v_j, v_i\rangle = \overline{\langle v_i, v_j\rangle} = \overline{G_{ij}}$ (symmetric in the real case). For $\mathbf{x}\in\mathbb{C}^n$, conjugate-linearity in the first slot and linearity in the second give
$$\mathbf{x}^\dagger G\mathbf{x} = \sum_{i,j}\bar x_ix_j\langle v_i, v_j\rangle = \Big\langle\sum_ix_iv_i, \sum_jx_jv_j\Big\rangle = \Big\|\sum_ix_iv_i\Big\|^2\ge 0.$$
In a real space $\bar x_i = x_i$ and this reads $\mathbf{x}^\top G\mathbf{x}\ge 0$. (For $\mathbb{R}^m$ with the standard inner product and $v_i$ the columns of $A$, $G = A^\top A$.)
3. PSD ⇔ eigenvalues $\ge 0$: repeat steps 1–2 of Theorem 3.15 (or their Hermitian versions) with $\ge$ in place of $>$. ∎

■ **What just happened.** "PSD" is the matrix shadow of "squared lengths are never negative." Any matrix built as inner products of vectors with each other inherits it, and it fails to be PD exactly when some nontrivial combination $\sum_ix_iv_i$ of the vectors is zero.

The Gram matrix is where the SVD starts: LA-9 applies Theorem 3.11 to $A^\top A$, and the singular values of $A$ are $\sqrt{\lambda_i(A^\top A)}$, real because of Theorem 3.16. Exercise B5 shows $A^\top A$ and $AA^\top$ have the same nonzero eigenvalues, which LA-9 needs.

### 3.9 The Rayleigh quotient, Courant–Fischer, and the second-derivative test

**Theorem 3.17 (Rayleigh quotient bounds).** For symmetric $A$ with eigenvalues $\lambda_1\ge\cdots\ge\lambda_n$ and the Rayleigh quotient of §3.5,
$$R(\mathbf{x}) = \frac{\mathbf{x}^\top A\mathbf{x}}{\mathbf{x}^\top\mathbf{x}} \qquad(\mathbf{x}\ne\mathbf{0};\ \mathbf{x}^\dagger \text{ in the complex case}), \tag{3.13}$$
we have
$$\lambda_n \le R(\mathbf{x}) \le \lambda_1 \qquad\text{for all } \mathbf{x}\ne\mathbf{0},$$
with $R(\mathbf{x}) = \lambda_1$ exactly when $\mathbf{x}$ is an eigenvector for $\lambda_1$, and $R(\mathbf{x}) = \lambda_n$ exactly when $\mathbf{x}$ is an eigenvector for $\lambda_n$.

▶ **Try it first.** Put (3.12) over $\|\mathbf{y}\|^2$ and read the result as a weighted average.

**Derivation.**
1. With $\mathbf{y} = Q^\top\mathbf{x}$: $R(\mathbf{x}) = \dfrac{\sum_i\lambda_iy_i^2}{\sum_iy_i^2} = \sum_i\lambda_iw_i$, where $w_i = y_i^2/\|\mathbf{y}\|^2$ satisfy $w_i\ge 0$ and $\sum_iw_i = 1$.
2. A weighted average of numbers lies between their minimum and maximum: $\lambda_n = \sum_i\lambda_nw_i \le \sum_i\lambda_iw_i \le \sum_i\lambda_1w_i = \lambda_1$.
3. Equality on the right: $\lambda_1 - R(\mathbf{x}) = \sum_i(\lambda_1 - \lambda_i)w_i$ is a sum of nonnegative terms, so it is $0$ iff $w_i = 0$ for every $i$ with $\lambda_i < \lambda_1$, iff $\mathbf{y}$ is supported on the indices with $\lambda_i = \lambda_1$, iff $\mathbf{x} = Q\mathbf{y}$ lies in the $\lambda_1$-eigenspace. Same on the left. ∎

■ **What just happened.** $R(\mathbf{x})$ is the average of the eigenvalues weighted by how much of $\mathbf{x}$ lies along each eigenvector. Two reasons to care now, before LA-9 and ML-6 give more: (1) In quantum mechanics $R(\psi) = \langle\psi|A|\psi\rangle/\langle\psi|\psi\rangle$ is the *expected value* of a measurement of $A$ in state $\psi$, and step 1 is literally "sum of outcomes times probabilities," with $w_i = |\langle q_i|\psi\rangle|^2$ the Born probabilities of §3.12. (2) In mechanics, the normal stress on a plane with unit normal $\mathbf{n}$ is $\mathbf{n}^\top\sigma\mathbf{n} = R(\mathbf{n})$, so the theorem says the normal stress on any plane is bounded by the principal stresses (Example 4.5).

**Every eigenvector is a critical point of $R$.** Theorem 3.17 is about the maximum and minimum, but the other eigenvectors are critical points too. By the quotient rule, with $\nabla(\mathbf{x}^\top A\mathbf{x}) = 2A\mathbf{x}$ (this uses $A^\top = A$, as in Part I step 3) and $\nabla(\mathbf{x}^\top\mathbf{x}) = 2\mathbf{x}$,
$$\nabla R(\mathbf{x}) = \frac{2A\mathbf{x}\,(\mathbf{x}^\top\mathbf{x}) - (\mathbf{x}^\top A\mathbf{x})\,2\mathbf{x}}{(\mathbf{x}^\top\mathbf{x})^2} = \frac{2\big(A\mathbf{x} - R(\mathbf{x})\,\mathbf{x}\big)}{\mathbf{x}^\top\mathbf{x}}.$$
This vanishes exactly when $A\mathbf{x} = R(\mathbf{x})\mathbf{x}$: the critical points of $R$ are exactly the eigenvectors, and the critical values are the eigenvalues. Eigenvectors for eigenvalues strictly between $\lambda_n$ and $\lambda_1$ are saddle points of $R$ (moving toward $\mathbf{q}_1$ raises $R$, moving toward $\mathbf{q}_n$ lowers it). LA-13 uses this to explain why the Rayleigh quotient is so accurate.

**Theorem 3.18 (Courant–Fischer min-max).** For symmetric $A$ with $\lambda_1\ge\cdots\ge\lambda_n$ and $1\le k\le n$,
$$\lambda_k = \max_{\substack{S\subseteq\mathbb{R}^n\\ \dim S = k}}\;\min_{\substack{\mathbf{x}\in S\\ \mathbf{x}\ne\mathbf{0}}} R(\mathbf{x}). \tag{3.14}$$

▶ **Try it first.** Take $k = 2$ in $\mathbb{R}^3$. Any plane $S$ through the origin meets the plane $T = \operatorname{span}\{\mathbf{q}_2, \mathbf{q}_3\}$ in at least a line, and on $T$ the weighted average $R$ is at most $\lambda_2$. What does that say about $\min_SR$? Which plane $S$ makes $\min_SR$ as large as possible?

*The first case, $k = 1$.* A one-dimensional $S$ is $\operatorname{span}\{\mathbf{x}\}$ and $R$ is constant on it, so the inner minimum is $R(\mathbf{x})$, and (3.14) says $\lambda_1 = \max_{\mathbf{x}\ne\mathbf{0}}R(\mathbf{x})$. That is Theorem 3.17. ∎

*General $k$* (one dimension count more). ($\ge$): take $S = \operatorname{span}\{\mathbf{q}_1, \dots, \mathbf{q}_k\}$. For $\mathbf{x}\in S$, $\mathbf{y}$ is supported on $i\le k$, so $R(\mathbf{x}) = \sum_{i\le k}\lambda_iw_i\ge\lambda_k$. Thus $\min_SR\ge\lambda_k$ and the max is at least $\lambda_k$. ($\le$): let $S$ be any $k$-dimensional subspace and $T = \operatorname{span}\{\mathbf{q}_k, \dots, \mathbf{q}_n\}$, of dimension $n - k + 1$. By the dimension formula (LA-3 Theorem 3.35, consequence (ii)), $\dim(S\cap T) = \dim S + \dim T - \dim(S + T) \ge k + (n-k+1) - n = 1$, so there is $\mathbf{x}\ne\mathbf{0}$ in $S\cap T$. Its $\mathbf{y}$ is supported on $i\ge k$, so $R(\mathbf{x}) = \sum_{i\ge k}\lambda_iw_i\le\lambda_k$. Hence $\min_SR\le\lambda_k$ for every $S$, and the max is at most $\lambda_k$. ∎

■ **What just happened.** $\lambda_k$ is described without mentioning any eigenvector: it is the best guaranteed value of $R$ on a $k$-dimensional subspace. The only tool was a dimension count forcing two subspaces to intersect. Because the description is eigenvector-free, it survives perturbations and restrictions, which is how it gets used:

**Corollary (Weyl's inequality).** For symmetric $X, Y$ of the same size, with eigenvalues in decreasing order and $\varepsilon = \max_i|\lambda_i(X - Y)|$,
$$|\lambda_k(X) - \lambda_k(Y)| \le \varepsilon \le \|X - Y\|_F \qquad\text{for every } k.$$

*Proof.* 1. Let $E = X - Y$. For $\mathbf{x}\ne\mathbf{0}$, $R_X(\mathbf{x}) = R_Y(\mathbf{x}) + R_E(\mathbf{x})$, and Theorem 3.17 for $E$ gives $-\varepsilon\le R_E(\mathbf{x})\le\varepsilon$. So $R_Y(\mathbf{x}) - \varepsilon\le R_X(\mathbf{x})\le R_Y(\mathbf{x}) + \varepsilon$ for every $\mathbf{x}\ne\mathbf{0}$.
2. Take the minimum over $\mathbf{x}$ in a $k$-dimensional $S$, then the maximum over $S$, in (3.14): $\lambda_k(Y) - \varepsilon\le\lambda_k(X)\le\lambda_k(Y) + \varepsilon$.
3. With $E = Q_E\Lambda_EQ_E^\top$, $\|E\|_F^2 = \operatorname{tr}(E^\top E) = \operatorname{tr}(Q_E\Lambda_E^2Q_E^\top) = \operatorname{tr}(\Lambda_E^2Q_E^\top Q_E) = \sum_i\lambda_i(E)^2\ge\varepsilon^2$. ∎

So a symmetric perturbation of size $\varepsilon$ moves every eigenvalue by at most $\varepsilon$: the eigenvalues of a symmetric matrix are perfectly conditioned. Exercise B7(c) uses this with $Y = \operatorname{diag}(X)$; LA-13 uses it to bound the error of computed eigenvalues; its singular-value cousin is LA-9 Exercise C2(c).

Exercise C1 uses (3.14) to prove Cauchy's interlacing theorem; LA-9 uses it to characterize singular values; ML-6 uses the $k = 2$ case to find the second principal component.

**Theorem 3.19 (Second-derivative test is the spectral theorem for the Hessian).** Let $f : \mathbb{R}^n\to\mathbb{R}$ be twice continuously differentiable, let $\mathbf{a}$ be a critical point ($\nabla f(\mathbf{a}) = \mathbf{0}$), and let $H$ be the Hessian at $\mathbf{a}$, $H_{ij} = \partial^2f/\partial x_i\partial x_j$. Then $H$ is symmetric, and:
- $H\succ 0$ (all $\lambda_i > 0$) ⇒ strict local minimum;
- $H\prec 0$ (all $\lambda_i < 0$) ⇒ strict local maximum;
- $H$ has eigenvalues of both signs, **singular or not** ⇒ saddle;
- $H$ semidefinite and singular (all $\lambda_i\ge 0$ or all $\lambda_i\le 0$, with some $\lambda_i = 0$) ⇒ no conclusion.

A singular Hessian is inconclusive only in the last case: $f = x^2 - y^2 + z^4$ has $H = \operatorname{diag}(2, -2, 0)$ at the origin, singular, and the origin is still a saddle.

▶ **Try it first.** Write the second-order Taylor expansion at $\mathbf{a}$ in the eigen-coordinates $\mathbf{y} = Q^\top\mathbf{h}$. What does $\tfrac12\sum_i\lambda_iy_i^2$ do along $\mathbf{q}_1$? Along $\mathbf{q}_n$? What can the remainder do when it is not dominated?

**Derivation.**
1. $H$ is symmetric because mixed partials of a $C^2$ function commute (Clairaut/Schwarz theorem, multivariable calculus). So Theorem 3.11 applies: $H = Q\Lambda Q^\top$.
2. Taylor's theorem to second order at a critical point: $f(\mathbf{a} + \mathbf{h}) = f(\mathbf{a}) + \tfrac12\mathbf{h}^\top H\mathbf{h} + r(\mathbf{h})$ with $r(\mathbf{h})/\|\mathbf{h}\|^2\to 0$ as $\mathbf{h}\to\mathbf{0}$.
3. If $H\succ 0$, Theorem 3.17 gives $\mathbf{h}^\top H\mathbf{h}\ge\lambda_n\|\mathbf{h}\|^2$ with $\lambda_n > 0$. For $\|\mathbf{h}\|$ small enough that $|r(\mathbf{h})| < \tfrac14\lambda_n\|\mathbf{h}\|^2$, we get $f(\mathbf{a}+\mathbf{h}) - f(\mathbf{a}) \ge \tfrac14\lambda_n\|\mathbf{h}\|^2 > 0$: strict local minimum. Negative definite: apply this to $-f$.
4. If $H$ has eigenvalues of both signs, $\lambda_1 > 0 > \lambda_n$ (whether or not some other $\lambda_i$ is $0$; this step never needs $H$ invertible). Along $\mathbf{h} = t\mathbf{q}_1$, $f(\mathbf{a}+\mathbf{h}) - f(\mathbf{a}) = \tfrac12\lambda_1t^2 + r \ge \tfrac14\lambda_1 t^2 > 0$ for small $t$; along $\mathbf{h} = t\mathbf{q}_n$ it is $< 0$. $f$ goes up in one direction and down in a perpendicular one: a saddle, and the eigenvectors tell you the up and down directions.
5. If $H$ is semidefinite with some $\lambda_i = 0$, the quadratic term vanishes along $\mathbf{q}_i$ and $r$ decides; $f = x^4$, $-x^4$, $x^3$ all have $H = 0$ at $0$ with different behavior (minimum, maximum, neither). ∎

■ **What just happened.** At a critical point, $f$ looks like the quadratic form $\tfrac12\mathbf{h}^\top H\mathbf{h}$ to second order, and Theorem 3.14 sorts quadratic forms by the signs of their eigenvalues. Positive and negative eigenvalues are robust: the remainder is too small to overturn them. Zero eigenvalues are directions where the quadratic form has nothing to say. If the nonzero eigenvalues include both signs, the saddle is already settled; only when they all share one sign does the verdict pass to the higher-order terms.

In two variables this is the textbook rule. $\det H = \lambda_1\lambda_2$, so $\det H > 0$ means the eigenvalues have the same sign. And $f_{xx} = \mathbf{e}_1^\top H\mathbf{e}_1 = R(\mathbf{e}_1)$ lies between $\lambda_2$ and $\lambda_1$ (Theorem 3.17), so when both have the same sign, $f_{xx}$ has it too. "$f_{xx}f_{yy} - f_{xy}^2 > 0$ and $f_{xx} > 0$ ⇒ minimum" is Theorem 3.15 (iv): the two pivots $f_{xx}$ and $\det H/f_{xx}$ are positive. Exercise B3 works a case where the test fails and a case where it doesn't; LA-14 and ML-2 build on this.

### 3.10 Unitary, orthogonal, and normal matrices

**Theorem 3.20 (Unitary matrices).** For $U\in\mathbb{C}^{n\times n}$ the following are equivalent: (a) $U^\dagger U = I$; (b) $\langle U\mathbf{u}, U\mathbf{v}\rangle = \langle\mathbf{u},\mathbf{v}\rangle$ for all $\mathbf{u},\mathbf{v}$; (c) the columns of $U$ are orthonormal. Each implies (d) $\|U\mathbf{x}\| = \|\mathbf{x}\|$ for all $\mathbf{x}$, and (d) implies them in turn (Exercise C2). A unitary matrix has $U^{-1} = U^\dagger$ (so also $UU^\dagger = I$), every eigenvalue satisfies $|\lambda| = 1$, eigenvectors for distinct eigenvalues are orthogonal, and $|\det U| = 1$. A real unitary matrix is called **orthogonal**, $Q^\top Q = I$.

▶ **Try it first.** Write $\langle U\mathbf{u}, U\mathbf{v}\rangle$ as a matrix product and use $(U\mathbf{u})^\dagger = \mathbf{u}^\dagger U^\dagger$. For the eigenvalues, take the length of both sides of $U\mathbf{v} = \lambda\mathbf{v}$.

**Derivation.**
1. (a)⇒(b): $\langle U\mathbf{u}, U\mathbf{v}\rangle = (U\mathbf{u})^\dagger U\mathbf{v} = \mathbf{u}^\dagger U^\dagger U\mathbf{v} = \mathbf{u}^\dagger\mathbf{v}$.
2. (b)⇒(a): $\langle\mathbf{e}_i, (U^\dagger U - I)\mathbf{e}_j\rangle = \langle U\mathbf{e}_i, U\mathbf{e}_j\rangle - \langle\mathbf{e}_i,\mathbf{e}_j\rangle = 0$ for all $i, j$; these are the entries of $U^\dagger U - I$.
3. (a)⇔(c): $(U^\dagger U)_{ij} = \mathbf{u}_i^\dagger\mathbf{u}_j = \langle\mathbf{u}_i,\mathbf{u}_j\rangle$ where $\mathbf{u}_j$ are the columns; $U^\dagger U = I$ says exactly $\langle\mathbf{u}_i,\mathbf{u}_j\rangle = \delta_{ij}$.
4. (b)⇒(d): take $\mathbf{u} = \mathbf{v}$.
5. Inverse: a square matrix with a left inverse has that as its two-sided inverse (LA-2), so $U^{-1} = U^\dagger$.
6. Eigenvalues: if $U\mathbf{v} = \lambda\mathbf{v}$, $\mathbf{v}\ne\mathbf{0}$, then by (d) $\|\mathbf{v}\| = \|U\mathbf{v}\| = |\lambda|\,\|\mathbf{v}\|$, so $|\lambda| = 1$.
7. Orthogonality: if also $U\mathbf{w} = \mu\mathbf{w}$ with $\mu\ne\lambda$, then by (b) $\langle\mathbf{v},\mathbf{w}\rangle = \langle U\mathbf{v}, U\mathbf{w}\rangle = \bar\lambda\mu\langle\mathbf{v},\mathbf{w}\rangle$. Since $|\lambda| = 1$, $\bar\lambda = 1/\lambda$, so $\bar\lambda\mu = \mu/\lambda\ne 1$, forcing $\langle\mathbf{v},\mathbf{w}\rangle = 0$.
8. $|\det U|^2 = \overline{\det U}\,\det U = \det U^\dagger\det U = \det(U^\dagger U) = 1$. The middle step is $\det A^\dagger = \overline{\det A}$: conjugating every entry conjugates each term of the Leibniz formula, and transposing does not change the determinant (LA-6, Theorem 3.11 and the complex-matrices remark after it). ∎

■ **What just happened.** Unitary matrices are the rigid motions of $\mathbb{C}^n$: they preserve every length and angle. Nothing gets stretched, so every eigenvalue lies on the unit circle. For a real orthogonal matrix the eigenvalues are $\pm1$ or complex-conjugate pairs $e^{\pm i\theta}$: a $2\times 2$ rotation by $\theta$ has eigenvalues $e^{\pm i\theta}$ (LA-7), and a reflection has $\pm1$. In QC-1 every single-qubit gate is a $2\times 2$ unitary. The change-of-basis matrix $Q$ in the spectral theorem is orthogonal, which is what made $Q^{-1} = Q^\top$ and made step 2 of Theorem 3.14 a rigid rotation.

**Definition (Normal).** $A$ is **normal** if $AA^* = A^*A$. Hermitian ($A^* = A$) and unitary ($A^*A = I = AA^*$) matrices are both normal.

**Theorem 3.21.** If $A = U\Lambda U^\dagger$ with $U$ unitary and $\Lambda$ diagonal (any complex diagonal), then $A$ is normal.

**Derivation.**
1. $A^\dagger = (U\Lambda U^\dagger)^\dagger = U\Lambda^\dagger U^\dagger = U\bar\Lambda U^\dagger$.
2. $AA^\dagger = U\Lambda U^\dagger U\bar\Lambda U^\dagger = U\Lambda\bar\Lambda U^\dagger$ and $A^\dagger A = U\bar\Lambda\Lambda U^\dagger$.
3. Diagonal matrices commute: $\Lambda\bar\Lambda = \bar\Lambda\Lambda = \operatorname{diag}(|\lambda_i|^2)$. So $AA^\dagger = A^\dagger A$. ∎

The converse is also true: **every normal matrix is unitarily diagonalizable**, and normality is exactly the condition for an orthonormal eigenbasis to exist. The proof goes through the Schur factorization $A = UTU^\dagger$ with $T$ upper triangular, then shows a normal triangular matrix is diagonal; LA-10 does it. Exercise C4 does the $2\times 2$ case by hand, and it contains the whole idea. A matrix that is diagonalizable but not normal, like $\begin{pmatrix}1&1\\0&2\end{pmatrix}$ with eigenvectors $(1,0)^\top$ and $(1,1)^\top$, has an eigenbasis that is not orthogonal, so you cannot write it as $U\Lambda U^\dagger$.

### 3.11 The Pauli matrices

**Proposition 3.22.** Each Pauli matrix is Hermitian and unitary, with $\sigma_k^2 = I$ and eigenvalues $+1, -1$. The eigenvectors are
$$\sigma_z: \ \mathbf{e}_1\ (+1),\ \mathbf{e}_2\ (-1); \qquad \sigma_x: \ \tfrac{1}{\sqrt2}(1, 1)^\top\ (+1),\ \tfrac1{\sqrt2}(1, -1)^\top\ (-1); \qquad \sigma_y: \ \tfrac1{\sqrt2}(1, i)^\top\ (+1),\ \tfrac1{\sqrt2}(1, -i)^\top\ (-1).$$
For a real unit vector $\mathbf{n} = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$, the matrix
$$\mathbf{n}\cdot\boldsymbol\sigma = n_x\sigma_x + n_y\sigma_y + n_z\sigma_z = \begin{pmatrix}\cos\theta & \sin\theta\,e^{-i\phi}\\ \sin\theta\,e^{i\phi} & -\cos\theta\end{pmatrix} \tag{3.15}$$
is Hermitian and unitary with eigenvalues $\pm1$ and orthonormal eigenvectors
$$|\mathbf{n},+\rangle = \begin{pmatrix}\cos(\theta/2)\\ e^{i\phi}\sin(\theta/2)\end{pmatrix}, \qquad |\mathbf{n},-\rangle = \begin{pmatrix}\sin(\theta/2)\\ -e^{i\phi}\cos(\theta/2)\end{pmatrix}.$$

▶ **Try it first.** Square $\sigma_y$ and $\mathbf{n}\cdot\boldsymbol\sigma$. A matrix that is Hermitian and squares to $I$ is also unitary; what do Theorems 3.9 and 3.20 then allow its eigenvalues to be? Use the trace to decide how many of each.

**Derivation.**
1. Hermitian: $\sigma_x$ and $\sigma_z$ are real symmetric. $\sigma_y^\dagger = \overline{\sigma_y}^\top = \begin{pmatrix}0&i\\-i&0\end{pmatrix}^\top = \begin{pmatrix}0&-i\\i&0\end{pmatrix} = \sigma_y$. For (3.15): the diagonal is real and the off-diagonal entries are conjugates of each other.
2. $\sigma_k^2 = I$ by direct multiplication; e.g. $\sigma_y^2 = \begin{pmatrix}(-i)(i) & 0\\ 0 & (i)(-i)\end{pmatrix} = I$. For (3.15), the top-left entry of $(\mathbf{n}\cdot\boldsymbol\sigma)^2$ is $\cos^2\theta + \sin^2\theta\,e^{-i\phi}e^{i\phi} = 1$, the top-right is $\cos\theta\sin\theta e^{-i\phi} - \sin\theta e^{-i\phi}\cos\theta = 0$, and the other two entries follow the same way. So $(\mathbf{n}\cdot\boldsymbol\sigma)^2 = I$.
3. Hermitian with $A^2 = I$ means $A^\dagger A = A^2 = I$: unitary. So the eigenvalues are real (Theorem 3.9) and on the unit circle (Theorem 3.20): $\lambda = \pm1$. The trace is $0$ and equals the eigenvalue sum (LA-7), so one eigenvalue is $+1$ and the other $-1$.
4. Eigenvectors for $\sigma_y$: $\sigma_y(1, i)^\top = (-i\cdot i, i\cdot 1)^\top = (1, i)^\top$ and $\sigma_y(1, -i)^\top = ((-i)(-i), i)^\top = (-1, i)^\top = -(1, -i)^\top$. Orthogonal: $\langle(1,i), (1,-i)\rangle = \bar1\cdot1 + \bar i\cdot(-i) = 1 + (-i)(-i) = 1 - 1 = 0$, remembering to conjugate the first slot.
5. For $|\mathbf{n},+\rangle$: top entry of $(\mathbf{n}\cdot\boldsymbol\sigma)|\mathbf{n},+\rangle$ is $\cos\theta\cos\tfrac\theta2 + \sin\theta e^{-i\phi}e^{i\phi}\sin\tfrac\theta2 = \cos(\theta - \tfrac\theta2) = \cos\tfrac\theta2$; bottom entry is $e^{i\phi}(\sin\theta\cos\tfrac\theta2 - \cos\theta\sin\tfrac\theta2) = e^{i\phi}\sin\tfrac\theta2$. Eigenvalue $+1$. For $|\mathbf{n},-\rangle$: top is $\cos\theta\sin\tfrac\theta2 - \sin\theta\cos\tfrac\theta2 = -\sin\tfrac\theta2$, bottom is $e^{i\phi}(\sin\theta\sin\tfrac\theta2 + \cos\theta\cos\tfrac\theta2) = e^{i\phi}\cos\tfrac\theta2$, which is $-1$ times the vector. Orthogonality: $\langle|\mathbf{n},+\rangle, |\mathbf{n},-\rangle\rangle = \cos\tfrac\theta2\sin\tfrac\theta2 + \overline{e^{i\phi}}\sin\tfrac\theta2\cdot(-e^{i\phi}\cos\tfrac\theta2) = 0$. Both have norm $\cos^2\tfrac\theta2 + \sin^2\tfrac\theta2 = 1$. ∎

■ **What just happened.** The half-angles are the first sign of something deep: the eigenvector of $\mathbf{n}\cdot\boldsymbol\sigma$ depends on $\theta/2$, so rotating $\mathbf{n}$ by $2\pi$ (a full turn, back to the same $\mathbf{n}$) sends $|\mathbf{n},+\rangle$ to $-|\mathbf{n},+\rangle$. This is the double cover $SU(2)\to SO(3)$, worked out in LA-12's appendix. In QC-1, the unit vector $\mathbf{n}$ is a point on the **Bloch sphere** and $|\mathbf{n},+\rangle$ is the qubit state it labels; the project plots a few.

### 3.12 Dirac notation: a preview of what QM needs

Physicists write a vector $\boldsymbol\psi\in\mathbb{C}^n$ as a **ket** $|\psi\rangle$ and its conjugate transpose $\boldsymbol\psi^\dagger$ as a **bra** $\langle\psi|$. Then:

- $\langle\phi|\psi\rangle = \boldsymbol\phi^\dagger\boldsymbol\psi = \langle\boldsymbol\phi,\boldsymbol\psi\rangle$: the inner product (3.1), conjugate-linear in $\phi$, linear in $\psi$. **This is the whole reason for our slot convention.**
- $|q\rangle\langle q| = \mathbf{q}\mathbf{q}^\dagger$: the projector $P$ of Theorem 3.13 onto the line through $\mathbf{q}$. Acting on $|\psi\rangle$ it gives $|q\rangle\langle q|\psi\rangle$, "the component along $q$ times $q$."
- (3.5) becomes $|\psi\rangle = \sum_i|q_i\rangle\langle q_i|\psi\rangle$, and $\sum_i|q_i\rangle\langle q_i| = I$ is **completeness**.
- (3.10) becomes $A = \sum_i\lambda_i|q_i\rangle\langle q_i|$, and (3.11) becomes $f(A) = \sum_if(\lambda_i)|q_i\rangle\langle q_i|$.
- The Rayleigh quotient is $\langle\psi|A|\psi\rangle/\langle\psi|\psi\rangle$.

Now the sentence that QM-1 opens with and QM-4 derives: **an observable is a Hermitian operator $A$; the possible outcomes of measuring it are the eigenvalues $\lambda_i$; the probability of outcome $\lambda_i$ in the normalized state $|\psi\rangle$ is**
$$p_i = |\langle q_i|\psi\rangle|^2. \tag{3.16}$$
Theorem 3.9 makes the outcomes real numbers. Theorem 3.12 makes the $|q_i\rangle$ an orthonormal basis, so Parseval (3.6) gives $\sum_ip_i = \|\psi\|^2 = 1$: the probabilities are a probability distribution. The expected outcome is $\sum_i\lambda_ip_i = \langle\psi|A|\psi\rangle$ by (3.10). For a qubit, measuring "spin along $\mathbf{n}$" is measuring $\mathbf{n}\cdot\boldsymbol\sigma$; the outcomes are $\pm1$ and the probabilities are $|\langle\mathbf{n},\pm|\psi\rangle|^2$. None of this is physics yet. It is the statement that the spectral theorem is what makes the recipe consistent, and that is why QM needs it.

## 4. Worked examples

### Example 4.1: Legendre polynomials by Gram–Schmidt

Inner product $\langle f, g\rangle = \int_{-1}^1 f(x)g(x)\,dx$ on polynomials. Run (3.7) on $v_0 = 1$, $v_1 = x$, $v_2 = x^2$, $v_3 = x^3$. The integrals you need: $\int_{-1}^1x^k\,dx = 2/(k+1)$ for even $k$ and $0$ for odd $k$.

*Step 0.* $u_0 = 1$. $\langle u_0, u_0\rangle = \int_{-1}^11\,dx = 2$.

*Step 1.* $\langle u_0, x\rangle = \int x\,dx = 0$, so $u_1 = x - 0\cdot u_0 = x$. $\langle u_1, u_1\rangle = \int x^2\,dx = 2/3$.

*Step 2.* $\langle u_0, x^2\rangle = 2/3$, $\langle u_1, x^2\rangle = \int x^3\,dx = 0$. So
$$u_2 = x^2 - \frac{2/3}{2}\cdot 1 - \frac{0}{2/3}\cdot x = x^2 - \tfrac13.$$
$\langle u_2, u_2\rangle = \int(x^4 - \tfrac23x^2 + \tfrac19)\,dx = \tfrac25 - \tfrac23\cdot\tfrac23 + \tfrac19\cdot2 = \tfrac25 - \tfrac49 + \tfrac29 = \tfrac25 - \tfrac29 = \tfrac{8}{45}$.

*Step 3.* $\langle u_0, x^3\rangle = 0$, $\langle u_1, x^3\rangle = \int x^4\,dx = 2/5$, $\langle u_2, x^3\rangle = \int(x^5 - \tfrac13x^3)\,dx = 0$. So
$$u_3 = x^3 - \frac{2/5}{2/3}\,x = x^3 - \tfrac35x.$$
$\langle u_3, u_3\rangle = \int(x^6 - \tfrac65x^4 + \tfrac9{25}x^2)\,dx = \tfrac27 - \tfrac65\cdot\tfrac25 + \tfrac9{25}\cdot\tfrac23 = \tfrac27 - \tfrac{12}{25} + \tfrac6{25} = \tfrac27 - \tfrac6{25} = \tfrac{8}{175}$.

*Normalization.* The Legendre polynomials $P_n$ are these $u_n$ rescaled so that $P_n(1) = 1$ (that is the classical convention, not unit norm). $u_2(1) = 2/3$, so $P_2 = \tfrac32u_2 = \tfrac12(3x^2 - 1)$. $u_3(1) = 2/5$, so $P_3 = \tfrac52u_3 = \tfrac12(5x^3 - 3x)$. Hence
$$P_0 = 1,\quad P_1 = x,\quad P_2 = \tfrac12(3x^2 - 1),\quad P_3 = \tfrac12(5x^3 - 3x),$$
and $\|P_n\|^2 = 2/(2n+1)$: check $\|P_2\|^2 = \tfrac94\cdot\tfrac8{45} = \tfrac25$ and $\|P_3\|^2 = \tfrac{25}4\cdot\tfrac8{175} = \tfrac27$. The orthonormal versions are $q_n = \sqrt{(2n+1)/2}\,P_n$. The project repeats this with numerical integration and compares to `scipy.special.legendre`.

Why you care: the $P_n$ are the orthogonal basis for expanding functions on $[-1,1]$ exactly as $\sin nx, \cos nx$ are on $[-\pi,\pi]$. Evaluated at $\cos\theta\in[-1,1]$, they are the $\theta$-dependence of the $m = 0$ spherical harmonics (QM-5), hence of hydrogen's $s$, $p_z$, $d_{z^2}$, ... orbitals (QM-6); associated Legendre functions handle $m\ne 0$. And the nodes of $P_n$ are the Gauss–Legendre quadrature points the project uses to integrate.

### Example 4.2: Full orthogonal diagonalization of a $3\times3$ symmetric matrix

$$A = \begin{pmatrix}2&1&1\\1&2&1\\1&1&2\end{pmatrix}.$$

*Eigenvalues.* Let $\mu = 2 - \lambda$. Then $\det(A - \lambda I) = \det\begin{pmatrix}\mu&1&1\\1&\mu&1\\1&1&\mu\end{pmatrix} = \mu^3 - 3\mu + 2$ (expand: $\mu(\mu^2 - 1) - (\mu - 1) + (1 - \mu)$). Factor: $\mu^3 - 3\mu + 2 = (\mu - 1)^2(\mu + 2)$. So $\mu = 1$ (double) or $\mu = -2$, i.e. $\lambda = 1, 1, 4$. All real, as Theorem 3.9 demands. Ordered: $\lambda_1 = 4$, $\lambda_2 = \lambda_3 = 1$.

*Eigenvector for $\lambda = 4$.* $(A - 4I)\mathbf{x} = \mathbf{0}$: $-2x_1 + x_2 + x_3 = 0$ and its two permutations. Solution $x_1 = x_2 = x_3$: $\mathbf{q}_1 = \tfrac1{\sqrt3}(1, 1, 1)^\top$.

*Eigenspace for $\lambda = 1$.* $(A - I)\mathbf{x} = \mathbf{0}$: all three rows say $x_1 + x_2 + x_3 = 0$. A plane, orthogonal to $\mathbf{q}_1$, as Theorem 3.10 promised. **A repeated eigenvalue hands you a subspace, not a basis; you must orthonormalize inside it.** Pick $\mathbf{v}_2 = (1, -1, 0)^\top$ and $\mathbf{v}_3 = (1, 0, -1)^\top$, both in the plane. Gram–Schmidt: $\mathbf{u}_3 = \mathbf{v}_3 - \dfrac{\mathbf{v}_2\cdot\mathbf{v}_3}{\mathbf{v}_2\cdot\mathbf{v}_2}\mathbf{v}_2 = (1, 0, -1)^\top - \tfrac12(1, -1, 0)^\top = (\tfrac12, \tfrac12, -1)^\top$. Normalize:
$$\mathbf{q}_2 = \tfrac1{\sqrt2}(1, -1, 0)^\top, \qquad \mathbf{q}_3 = \tfrac1{\sqrt6}(1, 1, -2)^\top.$$

*Verify $Q$ is orthogonal.* $\mathbf{q}_1\cdot\mathbf{q}_2 = (1 - 1 + 0)/\sqrt6 = 0$; $\mathbf{q}_1\cdot\mathbf{q}_3 = (1 + 1 - 2)/\sqrt{18} = 0$; $\mathbf{q}_2\cdot\mathbf{q}_3 = (1 - 1 + 0)/\sqrt{12} = 0$; each has unit length. So
$$Q = \begin{pmatrix}1/\sqrt3 & 1/\sqrt2 & 1/\sqrt6\\ 1/\sqrt3 & -1/\sqrt2 & 1/\sqrt6\\ 1/\sqrt3 & 0 & -2/\sqrt6\end{pmatrix}, \qquad \Lambda = \operatorname{diag}(4, 1, 1).$$

*Verify $A = Q\Lambda Q^\top$ via (3.10).* Let $J$ be the all-ones matrix. $P_1 = \mathbf{q}_1\mathbf{q}_1^\top = \tfrac13J$. Since $P_1 + P_2 + P_3 = I$, $P_2 + P_3 = I - \tfrac13J$. So
$$Q\Lambda Q^\top = 4P_1 + 1\cdot(P_2 + P_3) = \tfrac43J + I - \tfrac13J = I + J = A. \ ✓$$
(Indeed $A = I + J$ by inspection: diagonal $1 + 1 = 2$, off-diagonal $0 + 1 = 1$.)

*Functions of $A$ by (3.11).* Using $J^2 = 3J$ to check:
$A^{-1} = \tfrac14P_1 + (P_2 + P_3) = \tfrac1{12}J + I - \tfrac13J = I - \tfrac14J$. Check: $(I + J)(I - \tfrac14J) = I + J - \tfrac14J - \tfrac14J^2 = I + \tfrac34J - \tfrac34J = I$. ✓
$\sqrt A = 2P_1 + (P_2 + P_3) = \tfrac23J + I - \tfrac13J = I + \tfrac13J$. Check: $(I + \tfrac13J)^2 = I + \tfrac23J + \tfrac19\cdot3J = I + J$. ✓

The project runs Jacobi's algorithm on this $A$ and recovers $\Lambda$ and (up to signs and the rotation freedom inside the $\lambda = 1$ plane) $Q$.

### Example 4.3: Classify the conic $5x^2 + 4xy + 2y^2 = 1$

Matrix (symmetric part, splitting $4xy$ as $2xy + 2yx$): $A = \begin{pmatrix}5&2\\2&2\end{pmatrix}$. Trace $7$, determinant $10 - 4 = 6$, so $\lambda^2 - 7\lambda + 6 = 0$: $\lambda_1 = 6$, $\lambda_2 = 1$. Both positive: by Theorem 3.14 the curve is an **ellipse**.

*Eigenvectors.* $\lambda = 6$: $(A - 6I) = \begin{pmatrix}-1&2\\2&-4\end{pmatrix}$, so $x = 2y$: $\mathbf{q}_1 = \tfrac1{\sqrt5}(2, 1)^\top$. $\lambda = 1$: $(A - I) = \begin{pmatrix}4&2\\2&1\end{pmatrix}$, so $y = -2x$: $\mathbf{q}_2 = \tfrac1{\sqrt5}(1, -2)^\top$. Check $\mathbf{q}_1\cdot\mathbf{q}_2 = (2 - 2)/5 = 0$. ✓

*Axes and semi-axes.* Along $\mathbf{q}_1$ (direction $(2,1)$, angle $\arctan\tfrac12\approx26.57°$ above the $x$-axis) the semi-axis is $1/\sqrt{\lambda_1} = 1/\sqrt6\approx0.408$. Along $\mathbf{q}_2$ (direction $(1,-2)$, angle $\approx-63.43°$, equivalently $116.57°$) the semi-axis is $1/\sqrt{\lambda_2} = 1$. The **major axis** is along $\mathbf{q}_2$, the eigenvector of the *smaller* eigenvalue.

*Check a point.* The tip of the major axis is $\mathbf{q}_2\cdot1 = (1/\sqrt5, -2/\sqrt5)$: $5\cdot\tfrac15 + 4\cdot\tfrac{-2}{5} + 2\cdot\tfrac45 = 1 - \tfrac85 + \tfrac85 = 1$. ✓ The tip of the minor axis is $\mathbf{q}_1/\sqrt6 = (2, 1)/\sqrt{30}$: $(5\cdot4 + 4\cdot2 + 2\cdot1)/30 = 30/30 = 1$. ✓

*In words.* An ellipse centered at the origin, about 2.45 times longer than it is wide, with its long axis pointing from the upper-left to the lower-right at about $63°$ below the $x$-axis, and its short axis perpendicular to that at about $27°$ above the $x$-axis. In the rotated coordinates $y_1 = \mathbf{q}_1^\top\mathbf{x}$, $y_2 = \mathbf{q}_2^\top\mathbf{x}$ the equation is $6y_1^2 + y_2^2 = 1$, which is (3.12). The project draws it with the axes.

### Example 4.4: Eigen decomposition of $\sigma_y$ and of $\mathbf{n}\cdot\boldsymbol\sigma$ (Bloch sphere preview)

*$\sigma_y$.* From Proposition 3.22, eigenvalues $\pm1$ with eigenvectors $|{+y}\rangle = \tfrac1{\sqrt2}(1, i)^\top$ and $|{-y}\rangle = \tfrac1{\sqrt2}(1, -i)^\top$. The projectors:
$$|{+y}\rangle\langle{+y}| = \tfrac12\begin{pmatrix}1\\i\end{pmatrix}\begin{pmatrix}1 & -i\end{pmatrix} = \tfrac12\begin{pmatrix}1&-i\\i&1\end{pmatrix}, \qquad |{-y}\rangle\langle{-y}| = \tfrac12\begin{pmatrix}1&i\\-i&1\end{pmatrix}.$$
(Remember the bra conjugates: $\langle{+y}| = \tfrac1{\sqrt2}(1, -i)$.) Spectral decomposition (3.10): $(+1)\tfrac12\begin{pmatrix}1&-i\\i&1\end{pmatrix} + (-1)\tfrac12\begin{pmatrix}1&i\\-i&1\end{pmatrix} = \begin{pmatrix}0&-i\\i&0\end{pmatrix} = \sigma_y$. ✓ Completeness: the two projectors sum to $I$. ✓

*A tilted axis.* Take $\theta = \pi/3$, $\phi = \pi/2$, so $\mathbf{n} = (0, \tfrac{\sqrt3}2, \tfrac12)$ and by (3.15)
$$\mathbf{n}\cdot\boldsymbol\sigma = \begin{pmatrix}\tfrac12 & -i\tfrac{\sqrt3}2\\ i\tfrac{\sqrt3}2 & -\tfrac12\end{pmatrix}.$$
Proposition 3.22 says $|\mathbf{n},+\rangle = (\cos\tfrac\pi6, e^{i\pi/2}\sin\tfrac\pi6)^\top = (\tfrac{\sqrt3}2, \tfrac i2)^\top$. Check: top entry of the product is $\tfrac12\cdot\tfrac{\sqrt3}2 + (-i\tfrac{\sqrt3}2)(\tfrac i2) = \tfrac{\sqrt3}4 + \tfrac{\sqrt3}4 = \tfrac{\sqrt3}2$; bottom is $(i\tfrac{\sqrt3}2)(\tfrac{\sqrt3}2) + (-\tfrac12)(\tfrac i2) = \tfrac{3i}4 - \tfrac i4 = \tfrac i2$. Eigenvalue $+1$. ✓ And $|\mathbf{n},-\rangle = (\sin\tfrac\pi6, -e^{i\pi/2}\cos\tfrac\pi6)^\top = (\tfrac12, -i\tfrac{\sqrt3}2)^\top$.

*Born probabilities.* If the qubit is in state $|\psi\rangle = |0\rangle = \mathbf{e}_1$ and you measure spin along this $\mathbf{n}$, (3.16) gives $p_+ = |\langle\mathbf{n},+|0\rangle|^2 = \cos^2\tfrac\pi6 = \tfrac34$ and $p_- = \sin^2\tfrac\pi6 = \tfrac14$. They sum to $1$ by Parseval. In general $p_\pm = \tfrac12(1\pm\cos\theta)$, which depends only on the angle between $\mathbf{n}$ and the $z$-axis, the Bloch-sphere picture that QC-1 develops.

### Example 4.5: The 2D stress tensor, principal stresses, and Mohr's circle

A plane-stress state at a point is a symmetric $2\times2$ matrix (symmetric because of moment balance on an infinitesimal element):
$$\sigma = \begin{pmatrix}\sigma_{xx} & \tau_{xy}\\ \tau_{xy} & \sigma_{yy}\end{pmatrix} = \begin{pmatrix}80 & 30\\ 30 & 20\end{pmatrix}\ \text{MPa}.$$
The traction on a plane with unit normal $\mathbf{n}$ is $\sigma\mathbf{n}$; its normal component is $\sigma_n = \mathbf{n}^\top\sigma\mathbf{n}$ (the Rayleigh quotient!) and its shear component is $\tau_n = \mathbf{t}^\top\sigma\mathbf{n}$ with $\mathbf{t}\perp\mathbf{n}$.

*Principal stresses = eigenvalues.* Mean $m = \tfrac12(80 + 20) = 50$, half-difference $d = \tfrac12(80 - 20) = 30$. Characteristic equation $\lambda^2 - 100\lambda + (1600 - 900) = 0$, so $\lambda = 50\pm\sqrt{2500 - 700} = 50\pm\sqrt{1800} = 50\pm30\sqrt2$:
$$\sigma_1 = 50 + 30\sqrt2 \approx 92.43\ \text{MPa}, \qquad \sigma_2 = 50 - 30\sqrt2\approx7.57\ \text{MPa}.$$
Both tensile. Radius $r = \sqrt{d^2 + \tau_{xy}^2} = \sqrt{900 + 900} = 30\sqrt2\approx42.43$.

*Principal directions = eigenvectors.* For $\sigma_1$: $(80 - \sigma_1)x + 30y = 0$ gives $y/x = (\sigma_1 - 80)/30 = (30\sqrt2 - 30)/30 = \sqrt2 - 1 = \tan22.5°$. So $\mathbf{q}_1 = (\cos22.5°, \sin22.5°)^\top\approx(0.924, 0.383)^\top$, and $\mathbf{q}_2$ is perpendicular at $112.5°$. The textbook formula $\tan2\theta_p = 2\tau_{xy}/(\sigma_{xx} - \sigma_{yy}) = 60/60 = 1$ gives $\theta_p = 22.5°$. Same thing. Check: $\sigma\mathbf{q}_1 = (80\cdot0.924 + 30\cdot0.383,\ 30\cdot0.924 + 20\cdot0.383)^\top = (85.4, 35.4)^\top$ and $\sigma_1\mathbf{q}_1 = 92.43\cdot(0.924, 0.383)^\top = (85.4, 35.4)^\top$. ✓ On the principal planes the shear stress is zero: that is what "eigenvector" means here ($\sigma\mathbf{n}\parallel\mathbf{n}$, no tangential part).

*Mohr's circle is the spectral theorem.* With $\mathbf{n} = (\cos\theta, \sin\theta)^\top$ and $\mathbf{t} = (-\sin\theta, \cos\theta)^\top$, expand:
$$\sigma_n = \sigma_{xx}\cos^2\theta + \sigma_{yy}\sin^2\theta + 2\tau_{xy}\sin\theta\cos\theta = m + d\cos2\theta + \tau_{xy}\sin2\theta,$$
$$\tau_n = -d\sin2\theta + \tau_{xy}\cos2\theta,$$
using $\cos^2\theta = \tfrac12(1 + \cos2\theta)$, $\sin^2\theta = \tfrac12(1 - \cos2\theta)$, $2\sin\theta\cos\theta = \sin2\theta$. Then $(\sigma_n - m)^2 + \tau_n^2 = d^2 + \tau_{xy}^2 = r^2$. As $\theta$ runs over all planes, $(\sigma_n, \tau_n)$ traces a circle centered at $(m, 0) = (50, 0)$ with radius $r = 42.43$: **Mohr's circle**. Its rightmost and leftmost points are $m\pm r = \sigma_1, \sigma_2$, the eigenvalues; those are exactly the bounds $\lambda_{\min}\le R(\mathbf{n})\le\lambda_{\max}$ of Theorem 3.17. The top of the circle is the maximum shear $\tau_{\max} = r = 42.43$ MPa, on planes at $45°$ to the principal directions, i.e. $\theta = 22.5°\pm45°$. In the eigenbasis the circle is $\sigma_n = m + r\cos2\theta'$, $\tau_n = -r\sin2\theta'$ with $\theta'$ measured from $\mathbf{q}_1$: (3.12) in polar form. The project draws it by sweeping $\theta$ through all planes, asserts the circle equation, and marks the eigenvalues at its two ends. In 3D the stress tensor is a symmetric $3\times3$ matrix, Theorem 3.11 gives three perpendicular principal directions, and the same construction gives three Mohr circles.

### Example 4.6: Four tests for positive definiteness on one matrix

$$A = \begin{pmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{pmatrix},$$
the second-difference matrix (it discretizes $-d^2/dx^2$ with zero boundary conditions). You meet it again in F2, as the stiffness matrix of the beads-on-a-string preview of CM-6, and in CM-6 itself, where it is $K$ for $N$ equal masses on a string and its continuum limit is the wave equation.

(i) *Eigenvalues.* The eigenvectors are sampled sines: for $\mathbf{v}_k = (\sin\tfrac{k\pi}4, \sin\tfrac{2k\pi}4, \sin\tfrac{3k\pi}4)^\top$, $A\mathbf{v}_k = (2 - 2\cos\tfrac{k\pi}4)\mathbf{v}_k$ (check the middle row: $-\sin\tfrac{k\pi}4 + 2\sin\tfrac{2k\pi}4 - \sin\tfrac{3k\pi}4 = 2\sin\tfrac{2k\pi}4 - 2\sin\tfrac{2k\pi}4\cos\tfrac{k\pi}4$ by the sum-to-product identity). So $\lambda = 2 - \sqrt2,\ 2,\ 2 + \sqrt2 \approx 0.586, 2, 3.414$. All positive. ✓

(ii) *Quadratic form.* $\mathbf{x}^\top A\mathbf{x} = 2x_1^2 - 2x_1x_2 + 2x_2^2 - 2x_2x_3 + 2x_3^2 = x_1^2 + (x_1 - x_2)^2 + (x_2 - x_3)^2 + x_3^2 > 0$ unless $\mathbf{x} = \mathbf{0}$. ✓

(iv) *Pivots.* Eliminate: pivot $2$; row 2 becomes $(0, 2 - \tfrac12, -1) = (0, \tfrac32, -1)$, pivot $\tfrac32$; row 3 becomes $(0, 0, 2 - \tfrac1{3/2}) = (0, 0, \tfrac43)$, pivot $\tfrac43$. All positive. ✓ Their product $2\cdot\tfrac32\cdot\tfrac43 = 4 = \det A = (2-\sqrt2)\cdot2\cdot(2+\sqrt2) = 2\cdot2 = 4$. ✓ (Pivots and eigenvalues differ but have the same product and the same signs. The second fact is LA-10's proof of (iv).)

(iii) *Cholesky.* $A = R^\top R$ with $R$ upper triangular whose diagonal is the square roots of the pivots:
$$R = \begin{pmatrix}\sqrt2 & -1/\sqrt2 & 0\\ 0 & \sqrt{3/2} & -\sqrt{2/3}\\ 0 & 0 & \sqrt{4/3}\end{pmatrix}.$$
Check $(R^\top R)_{12} = \sqrt2\cdot(-1/\sqrt2) = -1$, $(R^\top R)_{22} = \tfrac12 + \tfrac32 = 2$, $(R^\top R)_{23} = \sqrt{3/2}\cdot(-\sqrt{2/3}) = -1$, $(R^\top R)_{33} = \tfrac23 + \tfrac43 = 2$. ✓ The project checks all four on this matrix, and runs the eigenvalue, Cholesky, and pivot tests on random SPD and non-SPD matrices, confirming they always agree.

## 5. Why it matters and where it connects

- **LA-9 (SVD).** Apply Theorem 3.11 to $A^\top A$ (symmetric PSD by Theorem 3.16). Its orthonormal eigenvectors are the right singular vectors, $\sqrt{\lambda_i}$ are the singular values, and Theorem 3.17 for $A^\top A$ says $\max\|A\mathbf{x}\|/\|\mathbf{x}\| = \sigma_1$: the operator norm asked about in LA-1 Exercise C3. Exercise B5 here is the lemma LA-9 needs about $A^\top A$ versus $AA^\top$.
- **LA-10 (decompositions).** Cholesky $A = R^\top R$ is Theorem 3.15 (iii) with a triangular $B$; Schur $A = UTU^\dagger$ proves the converse of Theorem 3.21; the polar decomposition is $\sqrt{A^\top A}$ from (3.11).
- **LA-12.** Quadratic forms are the symmetric bilinear forms of LA-12; Sylvester's law of inertia says the *signs* in (3.12) are basis-independent. The half-angles of Proposition 3.22 become the $SU(2)\to SO(3)$ double cover.
- **LA-13.** The symmetric QR algorithm is how (3.9) is computed in practice (Jacobi's algorithm, §7, is the older alternative). Weyl's inequality, the corollary to Courant–Fischer (Theorem 3.18), is why eigenvalues of symmetric matrices are perfectly conditioned, and LA-13 uses it to bound the error of computed eigenvalues; LA-13 also uses the critical-point property of $R$ after Theorem 3.17.
- **LA-14 / ML-2.** Theorem 3.19 reappears as LA-14 Theorem 3.14 (the second-derivative test in matrix-calculus form), and ML-2 upgrades it from one critical point to every point: a $C^2$ function is convex iff its Hessian is PSD everywhere; the condition number $\lambda_1/\lambda_n$ of the Hessian sets the speed of gradient descent, and the principal axes of Example 4.3 are the directions it zigzags between.
- **F2 (Fourier series).** The functions $e^{inx}/\sqrt{2\pi}$ are orthonormal for (3.2) on $[-\pi,\pi]$; Fourier's coefficients are (3.5); Parseval is (3.6); Legendre (Example 4.1) is the same story on $[-1,1]$.
- **F5 (multivariate Gaussian).** A covariance matrix $\Sigma = \mathbb{E}[(\mathbf{x}-\boldsymbol\mu)(\mathbf{x}-\boldsymbol\mu)^\top]$ is PSD by Theorem 3.16 (it is a Gram matrix). The density's level sets $(\mathbf{x}-\boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{x}-\boldsymbol\mu) = c$ are the ellipsoids of Theorem 3.14 with semi-axes $\sqrt{c\lambda_i(\Sigma)}$, and $\Sigma^{1/2}$ from (3.11) is how you sample from it.
- **ML-4, ML-6.** A kernel matrix $K_{ij} = k(\mathbf{x}_i, \mathbf{x}_j)$ is a Gram matrix and hence PSD by Theorem 3.16; that is what "valid kernel" means. PCA is Courant–Fischer (3.14) applied to the sample covariance: the first principal component maximizes the Rayleigh quotient (variance), the second maximizes it over the orthogonal complement.
- **ST-3.** $X^\top X$ in least squares is PD iff the columns of $X$ are independent (Theorem 3.16); ridge regression adds $\lambda I$ to push every eigenvalue up by $\lambda$ via (3.11).
- **CM-6 (normal modes).** The equations $M\ddot{\mathbf{x}} = -K\mathbf{x}$ lead to $K\mathbf{v} = \omega^2M\mathbf{v}$ with $M$ PD and $K$ symmetric; Exercises B1 and C3 show this is a self-adjoint problem in the $M$-inner product, so $\omega^2$ is real and the modes are $M$-orthogonal.
- **CM-7 (rigid bodies).** The inertia tensor is symmetric; Theorem 3.11 is why principal axes of rotation exist, and the moments of inertia are its eigenvalues.
- **QM-1, QM-4.** §3.12 verbatim, then extended to infinite dimensions. The uncertainty principle is Cauchy–Schwarz (Theorem 3.4) for two particular vectors.
- **QC-1.** Qubit states are unit vectors in $\mathbb{C}^2$; gates are unitary (Theorem 3.20); the Bloch sphere is Proposition 3.22; measurement is (3.16).

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 drill the Frobenius and function inner products, a $2\times2$ spectral decomposition, classifying quadratic forms, the symmetric-versus-Hermitian distinction, and Rayleigh quotients. B1–B8 are the module: B1 (adjoint under a weighted inner product), B3 (Hessian classification), B5 ($A^\top A$ versus $AA^\top$), B7 (Jacobi convergence, code), and B8 (redo the invariance step of the spectral theorem, and watch it fail for a non-symmetric matrix) are essential. C1–C4: Cauchy interlacing from Courant–Fischer, polarization and isometries, the generalized eigenproblem of CM-6, and the $2\times2$ case of "normal ⇒ unitarily diagonalizable."

## 7. Coding project

`project.py` builds and checks everything numerical in this module:

1. **Inner products.** Frobenius and function inner products, Cauchy–Schwarz on random data, Parseval (3.6) with a random orthonormal basis, and the orthogonality of $\sin mx, \sin nx$.
2. **Legendre by Gram–Schmidt.** Runs (3.7) on $1, x, \dots, x^5$ with the inner product $\int_{-1}^1$ computed by Gauss–Legendre quadrature, confirms Example 4.1's coefficients exactly, and compares to `scipy.special.legendre`.
3. **Jacobi's eigenvalue algorithm.** Repeatedly finds the largest off-diagonal entry $a_{pq}$ and applies a Givens rotation in the $(p, q)$ plane that zeroes it. The rotation angle is derived in a markdown cell: $\tan2\theta = 2a_{pq}/(a_{pp} - a_{qq})$. Each rotation is an orthogonal similarity, so the eigenvalues are preserved and the accumulated product of rotations converges to $Q$. Checked against `np.linalg.eigh` on Example 4.2 and random symmetric matrices; the squared off-diagonal norm $\operatorname{off}(A)^2$ decreases by exactly $2a_{pq}^2$ per step (Exercise B7), which is asserted after every rotation and plotted.
4. **Spectral decomposition.** $A = \sum\lambda_iP_i$, $\sqrt A$, $A^{-1}$, and $e^A$ via (3.11), checked against direct computation and `scipy.linalg.expm` on Example 4.2.
5. **Quadratic-form visualizer.** Contour plots of $\mathbf{x}^\top A\mathbf{x}$ for a PD (Example 4.3), an indefinite, and a PSD matrix, with eigenvector axes drawn and the semi-axis tips $1/\sqrt{\lambda_i}$ marked.
6. **Rayleigh quotient on the unit circle** for a $2\times2$: $R(\theta)$ plotted, maximum and minimum confirmed to be the eigenvalues at the eigenvector angles.
7. **Positive-definiteness tests.** Eigenvalues, Cholesky success/failure, and pivots (hand-written elimination) on Example 4.6, on random SPD matrices $B^\top B + 0.05I$ (SPD by Theorem 3.15 (iii), shifted to stay away from singular), and on random symmetric matrices that are usually indefinite; all three verdicts agree every time.
8. **Pauli matrices and the Bloch sphere.** Hermitian and unitary checks, eigenvalues $\pm1$, the eigenvectors of Proposition 3.22 verified for random $(\theta, \phi)$, the Born probabilities of Example 4.4, random unitaries having eigenvalues on the unit circle, and a 3D plot of the sphere with eigenvector directions.
9. **Mohr's circle** for Example 4.5: $(\sigma_n, \tau_n)$ computed on every plane, the circle equation asserted, and its two ends asserted to be the eigenvalues (the Rayleigh bounds of Theorem 3.17) at $\theta_p = 22.5°$.

**Build it yourself first:** write `jacobi_rotate(A, p, q)` (return the rotated matrix and the rotation used) and `gram_schmidt_poly(k)` (return the first $k$ orthogonal polynomials on $[-1,1]$) before reading the versions in the file. Then the "Your turn" section asks you to show numerically that symmetric and antisymmetric matrices are orthogonal complements under (3.3), to compute the Hessian of $f(x,y) = x^3 - 3xy^2$ with SymPy and classify its critical point by Theorem 3.19 (Exercise B3), to confirm that $A^\top A$ and $AA^\top$ share nonzero eigenvalues (Exercise B5), to write cyclic Jacobi (Exercise B7(d)), and to watch Weyl's inequality pin the diagonal of each Jacobi iterate to the eigenvalues (Exercise B7(c)). Reference implementations follow the prompts.

## 8. Where you'll struggle

- **The conjugate-linear slot.** You will get it backwards at least once, and the symptom is a sign error in an imaginary part or a probability that isn't real. Our rule: *conjugate the first slot*, $\langle\mathbf{u},\mathbf{v}\rangle = \mathbf{u}^\dagger\mathbf{v}$, so that the coefficient in (3.5) is $\langle q_i, v\rangle$ and Dirac's $\langle q_i|v\rangle$ is the same symbol. We chose this because QM-1 through QM-9 and QC-1 through QC-9 all use it, and translating in your head for eighteen modules is worse than translating once when you open Axler. When you read Axler, his $\langle u, v\rangle$ is our $\langle v, u\rangle$; his "linear in the first slot" is our "linear in the second." Theorem 3.9's proof is the place to test yourself: in step 2, the $\bar\lambda$ comes out of the *first* slot. If you had put $\lambda$ there, you would have "proved" nothing.
- **The induction step of the spectral theorem.** Part I (one eigenvector exists) is a calculus argument; people follow it. Part II is where they nod and don't actually see it. The claim is that $A$ maps $\mathbf{q}_1^\perp$ into itself, and the proof is the one-line computation in step 7. Do Exercise B8: redo step 7, then try it on the non-symmetric matrix of B8(c) (or the shear of step 7) and watch $A\mathbf{e}_2$ leave $\mathbf{e}_1^\perp = \operatorname{span}\{\mathbf{e}_2\}$. Then look at the block matrix in step 8 and say out loud what each zero block means. If the zero in the top-right still looks like algebra rather than like "$A$ doesn't leak out of $W$," you are not done.
- **"Symmetric" is not enough over $\mathbb{C}$.** The matrix $\begin{pmatrix}1&i\\i&-1\end{pmatrix}$ is symmetric and squares to zero; it has eigenvalue $0$ twice and only one eigenvector. The conclusion of the spectral theorem (Theorem 3.12) fails for it: no orthonormal eigenbasis, not even a non-orthogonal one. (Its eigenvalues happen to be real, but nothing guaranteed that; Exercise A5's $S$ has eigenvalues $\pm i$.) The condition that works is $A^\dagger = A$. Conversely $\sigma_y$ is Hermitian but *not* symmetric ($\sigma_y^\top = -\sigma_y$), and it behaves perfectly. Exercise A5.
- **PSD versus PD.** "Semidefinite" allows $\mathbf{x}^\top A\mathbf{x} = 0$ for nonzero $\mathbf{x}$; equivalently a zero eigenvalue, a nontrivial null space, a Gram matrix of *dependent* vectors (Theorem 3.16), a level set that is a cylinder rather than an ellipsoid. Covariance matrices of data with fewer samples than dimensions are PSD, not PD, and that is why ML-6 needs the pseudoinverse. Also: definiteness is about eigenvalues, not entries. $\begin{pmatrix}1&2\\2&1\end{pmatrix}$ has all positive entries and is indefinite; $\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ has negative entries and is PSD.
- **The Rayleigh quotient feels unmotivated.** It is introduced here because it proves the spectral theorem (Part I) and bounds quadratic forms (Theorem 3.17), and both of those are useful but not gripping. It becomes gripping in LA-9, where $\max R$ for $A^\top A$ is the squared operator norm $\sigma_1^2$; in ML-6, where it is the variance captured by a principal component; in QM-4, where it is the expected value of a measurement; and in CM-6 and LA-13 as the basis of every iterative eigenvalue method. For now, hold onto the mechanics reading: $R(\mathbf{n})$ is the normal stress on the plane with normal $\mathbf{n}$, and Mohr's circle is its graph.
- **Repeated eigenvalues.** When an eigenvalue has multiplicity $k$, the eigenspace is $k$-dimensional and you must run Gram–Schmidt *inside it* to get orthonormal columns of $Q$ (Example 4.2). Students take any two eigenvectors and then wonder why $Q^\top Q\ne I$. Theorem 3.10 guarantees orthogonality only *across* distinct eigenvalues.
- **Semi-axes are $1/\sqrt{\lambda_i}$, not $\lambda_i$ or $\sqrt{\lambda_i}$.** Large eigenvalue, short axis. For a covariance ellipse in F5 it is the other way ($\Sigma^{-1}$ appears in the exponent), so there large eigenvalue means long axis. Derive it each time rather than memorizing.
- **Eigenvalue ordering.** These notes order $\lambda_1\ge\cdots\ge\lambda_n$. `np.linalg.eigh` returns *ascending* order. The project reverses the output; if you forget, "$\mathbf{q}_1$" is the *bottom* eigenvector and Courant–Fischer looks false.

## 9. Check yourself

1. State the three axioms of a complex inner product. Which slot is conjugate-linear in this curriculum, and why was that chosen?
2. Write the Frobenius inner product of $A$ and $B$ in two ways and prove they agree.
3. Prove Cauchy–Schwarz for $\int_a^b\bar fg$ in five lines, from the axioms alone.
4. If $\{q_i\}$ is an orthonormal basis, what is the coefficient of $q_i$ in $v$? Write Parseval's identity, and use it to say why the Born probabilities (3.16) sum to $1$.
5. Define the adjoint by its one defining equation, and show that it is $A^\dagger$ for the standard inner product.
6. Prove that a Hermitian matrix has real eigenvalues and that eigenvectors for distinct eigenvalues are orthogonal. Mark where the conjugate-linear slot is used.
7. In the spectral theorem's proof, why is $\mathbf{q}_1^\perp$ invariant under $A$? Give a $2\times2$ matrix for which it is not, and say what goes wrong.
8. Write $A = Q\Lambda Q^\top$ as a sum of projections. What is $|q\rangle\langle q|$ as a matrix, and what is $\sum_i|q_i\rangle\langle q_i|$? What is $\sqrt A$? What is $e^A$?
9. For $\mathbf{x}^\top A\mathbf{x} = 1$ with $A$ PD, what are the directions and lengths of the axes?
10. Give four equivalent conditions for positive definiteness. Prove any two implications.
11. Why is $\lambda_{\min}\le\mathbf{x}^\top A\mathbf{x}/\mathbf{x}^\top\mathbf{x}\le\lambda_{\max}$? What does this say about normal stresses?
12. Why are the eigenvalues of a unitary matrix on the unit circle? Why are the eigenvalues of $\sigma_y$ exactly $\pm1$?

## 10. Pair with

- **Axler, *Linear Algebra Done Right*, ch. 6–7.** The cleanest coordinate-free treatment of inner product spaces, adjoints, and the spectral theorem (his real-case proof stays over $\mathbb{R}$, using the lemma that $T^2 + bT + cI$ is invertible for self-adjoint $T$ when $b^2 < 4c$ to find an eigenvalue: a different route from our Rayleigh-quotient argument, worth seeing). **Warning:** Axler is linear in the *first* slot. Translate.
- **Strang, *Introduction to Linear Algebra*, §6.4–6.5.** Symmetric matrices and positive definite matrices, with the pivots–eigenvalues–energy triple that Example 4.6 follows. Strang's "energy" $\mathbf{x}^\top A\mathbf{x}$ language is the best intuition for positive definiteness in print.
- **MIT 18.06, lectures 25–28 (Strang, OCW).** Lecture 25 (symmetric matrices and positive definiteness) and 27 (positive definite matrices and minima, including the second-derivative test and the ellipse $\mathbf{x}^\top A\mathbf{x} = 1$) are the ones to watch; 26 does complex matrices and the FFT; 28 starts similar matrices and Jordan form (LA-11). Watch 27 after §3.7–3.9.
- **Trefethen & Bau, *Numerical Linear Algebra*, lectures 24, 27, and 30** (optional). Lecture 24 sets up eigenvalue problems from the numerical side; lecture 27 covers the Rayleigh quotient and its quadratic accuracy (the critical-point remark after Theorem 3.17); lecture 30 covers the Jacobi algorithm the project implements, whose asymptotically quadratic convergence the project's Jacobi plot shows.
