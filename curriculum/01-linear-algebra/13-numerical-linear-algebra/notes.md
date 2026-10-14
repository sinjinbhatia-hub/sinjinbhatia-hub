# LA-13: Numerical Linear Algebra: Floating Point, Conditioning, Stability, Householder QR, Eigenvalue Algorithms, Krylov Methods, Randomized LA

**Time:** ~16 hours (the longest module in the course; budget three of those hours for §3.10 alone) · **Prerequisites:** LA-10 (and through it LA-5, LA-7, LA-8, LA-9) · **Unlocks:** every coding project in every later course; ML-2 (optimization) leans on §3.10 directly

> **Conventions in this module.** $\varepsilon_{\text{mach}} = 2^{-53} \approx 1.11\times 10^{-16}$ is the *unit roundoff* of IEEE double precision: the largest relative error of a single rounding. NumPy's `np.finfo(float).eps` $= 2^{-52}$ is the *gap* from $1$ to the next double, which is $2\varepsilon_{\text{mach}}$; Trefethen & Bau and Higham use our convention, many websites use NumPy's. $\mathrm{fl}(x)$ is the double nearest to the real number $x$. $\|\cdot\|$ with no subscript is the Euclidean norm on vectors and the induced $2$-norm $\|A\|_2 = \sigma_1(A)$ on matrices (LA-9). $\kappa(A) = \|A\|\,\|A^{-1}\| = \sigma_1/\sigma_n$ is the $2$-norm condition number. Householder vectors are chosen with the sign $v = x + \operatorname{sign}(x_1)\|x\|\mathbf{e}_1$, so the reflection sends $x$ to $-\operatorname{sign}(x_1)\|x\|\mathbf{e}_1$; the opposite sign is correct in exact arithmetic and wrong in floating point (§3.5). "SPD" means symmetric positive definite (LA-8). Flop counts are leading-order terms, with a multiply-add counted as two flops (Trefethen & Bau's convention). The Poisson grid has $N$ interior points per side and spacing $h = 1/(N+1)$.

## 0. Where this sits

You can now factor any matrix (LA-10) and you know what each factorization means. What you do not yet know is what happens when the factoring is done by a machine with 16 digits. This module is the bridge between the mathematics of LA-1 to LA-12 and the code you will write for the rest of the curriculum: every `np.linalg.solve`, every `scipy.sparse.linalg.cg`, every PCA in ML-6, every Hamiltonian diagonalization in QM-4 and QC-9 goes through the ideas here.

You need from earlier modules: QR and Gram–Schmidt (LA-5), eigenvalues and power iteration (LA-7), the spectral theorem and the Rayleigh quotient (LA-8), the SVD and $\kappa = \sigma_1/\sigma_n$ (LA-9), and the Schur form (LA-10). The module makes possible: the convergence proofs of gradient descent in ML-2 (the conjugate gradient derivation in §3.10 is the template), the "which solver" decisions in every physics project, and the honest reading of any numerical result you ever produce.

## 1. The problem that forced this idea

In 1943 Harold Hotelling estimated that solving $n$ linear equations by elimination could amplify rounding errors by a factor of order $4^n$. For $n = 100$ that is $10^{60}$. The first electronic computers were about to be built, and the people building them concluded that Gaussian elimination, the oldest algorithm in this course, would be useless on them. John von Neumann and Herman Goldstine (1947) showed it was fine for symmetric positive definite matrices. Alan Turing (1948), in a paper that introduced the word *condition number* and the LU factorization in its modern form, argued it was fine more generally. But nobody had a clean way to say *why*.

James Wilkinson supplied it in the early 1960s, and the idea is the heart of this module. He asked a different question. Not "how far is the computed answer from the true answer?" but "**is the computed answer the exact answer to a nearby problem?**" For Gaussian elimination with pivoting the answer is yes: the computed $\hat{x}$ satisfies $(A + \delta A)\hat{x} = b$ exactly, with $\delta A$ tiny (as long as the numbers in the elimination do not grow). That property is called **backward stability**. Whether $\hat{x}$ is then *close* to $x$ is a separate question about the *problem*: if $A$ is such that a tiny change $\delta A$ moves the solution a lot, no algorithm can help you. That separate sensitivity is the **condition number**. Hotelling had conflated the two. Separating them is the single most important thing you will learn here, and the roadmap warns you that it takes a few projects to internalize.

Wilkinson also discovered the opposite case. On the Pilot ACE computer he tested a root-finder on the polynomial $(x-1)(x-2)\cdots(x-20)$, expecting an easy success, and watched it fail completely. He later called this "the most traumatic experience in my career as a numerical analyst." The polynomial is perfectly well behaved as a *function*. But its *roots as a function of its coefficients* are catastrophically ill-conditioned: changing the coefficient of $x^{19}$ from $-210$ to $-210 - 2^{-23}$, a relative change of $6\times 10^{-10}$, turns ten of the twenty real roots into five complex-conjugate pairs with imaginary parts up to $\pm 2.81$ (§3.7, with the numbers). Every textbook of the time computed eigenvalues by forming the characteristic polynomial. Wilkinson's example killed that method and forced the field to find algorithms that never leave the matrix: the QR algorithm of John Francis and Vera Kublanovskaya (1961, independently) is the result, and it is still what `np.linalg.eig` runs.

Meanwhile a different line had already produced the best method for large systems and then forgotten it. Magnus Hestenes and Eduard Stiefel published the conjugate gradient method in 1952, as did Cornelius Lanczos in a closely related form. In exact arithmetic it solves an $n\times n$ symmetric positive definite system in at most $n$ steps, so it was filed as a *direct* method, and as a direct method it looked bad: in floating point the "exact after $n$ steps" property fails, and the bookkeeping is worse than elimination. It sat unused for nearly twenty years. In 1971 John Reid pointed out that for large sparse matrices you should stop after far fewer than $n$ steps, treating it as an *iterative* method, and that in that role nothing else came close. With preconditioning (Meijerink and van der Vorst, 1977) it became the standard solver for discretized PDEs, which is what it was for all along.

The most recent chapter is randomness. For matrices with thousands of columns, even a single pass over the data is expensive and a full SVD is out of the question. Nathan Halko, Per-Gunnar Martinsson and Joel Tropp (2011), building on work of Frieze, Kannan, Vempala, Rokhlin, Tygert and others, showed that multiplying $A$ by a *random* matrix with a few more columns than the rank you want captures the dominant singular subspace with high probability, with an error bound you can prove. It feels like cheating until you see the bound (§3.13).

So the module has four stories: finite digits and their arithmetic, the problem/algorithm split that Wilkinson made precise, the algorithms that survive it (Householder, QR iteration, conjugate gradients), and the modern trade of exactness for probability.

## 2. Intuition

### Digits are finite, and errors are relative

A double has 53 binary digits of precision, about 16 decimal digits. Every operation rounds its result to the nearest representable number. Rounding does not add a fixed absolute error; it multiplies by $(1 + \delta)$ with $|\delta| \le 1.1\times 10^{-16}$. So a number near $10^{20}$ is stored with absolute error up to $10^{4}$, and a number near $10^{-20}$ with absolute error $10^{-36}$. **Relative** error is the natural currency, and every bound in this module is a relative bound.

The one operation that breaks this picture is subtraction of nearly equal numbers. Take $x = 9000$ and $y = 8999.9998$, each known to 8 digits. Their difference, $0.0002$, is known to *one* digit: the seven leading digits that agreed carried no information about the difference, and they cancelled. Nothing went wrong in the subtraction itself (it was exact). What went wrong is that the inputs' tiny relative errors became a huge relative error of the output. This is **catastrophic cancellation**, and you should be able to guess the rule now: the amplification factor is $(|x| + |y|)/|x - y|$.

### Problem versus algorithm

Picture a function $f$ you want to evaluate, and a point $x$ where you want it. Two different things can go wrong.

```
   input x  ----(problem f)---->  output f(x)
     |                                 |
   tiny nudge δx                   how big is δf ?     <- CONDITIONING: a property of f at x
```

If a relative nudge of $10^{-16}$ in $x$ moves $f(x)$ by a relative $10^{-6}$, the problem is **ill-conditioned** at $x$, with condition number $\kappa \approx 10^{10}$. No algorithm, however clever, can give you more than 6 correct digits of $f(x)$ from a 16-digit $x$, because you do not *know* $x$ to better than 16 digits and the problem does not *determine* $f(x)$ to better than 6.

```
   input x  ----(algorithm f~)---->  computed f~(x)

   BACKWARD STABLE means:  f~(x) = f(x~)  for some x~ with |x~ - x| / |x| ≈ ε_mach
                           "the exact answer to a nearby question"
```

A good algorithm does not promise a small forward error $|\tilde f(x) - f(x)|$. It promises something humbler and achievable: its output is the *exact* value of $f$ at some input $\tilde x$ that is within rounding distance of $x$. If the problem is well-conditioned, that forces a small forward error. If it is ill-conditioned, the forward error is large, and that is the problem's fault, not the algorithm's.

Put the two together and you can guess the master formula of the module:

$$\text{forward error} \;\lesssim\; \kappa \times \text{backward error} \;\approx\; \kappa\, \varepsilon_{\text{mach}}.$$

The number of correct digits you get is $16$ minus $\log_{10}\kappa$. Memorize that sentence.

### Two ways to compute the same QR

Gram–Schmidt (LA-5) builds $q_j$ by subtracting from $a_j$ its components along $q_1,\dots,q_{j-1}$. When $a_j$ is nearly in the span of the earlier columns, that subtraction is a difference of nearly equal vectors: catastrophic cancellation, and the resulting $q_j$ is not orthogonal to the others to anything like 16 digits. Householder's method never subtracts. It *reflects*: it applies an orthogonal matrix to the whole column and reads off the result. Orthogonal matrices preserve lengths, so errors are never amplified, step after step. Guess: Householder QR is backward stable and Gram–Schmidt is not. §3.5 confirms it and the project shows Gram–Schmidt losing all orthogonality on a $7\times 7$ matrix while Householder sits at $10^{-15}$ on a $14\times 14$.

### Eigenvalues: you cannot finish, so iterate well

The eigenvalues of a $5\times 5$ matrix are the roots of a degree-5 polynomial, and Abel and Ruffini proved there is no formula for those in radicals. So every eigenvalue algorithm is iterative, and the question is only how fast it converges. You know one iteration from LA-7: multiply by $A$ repeatedly and the vector turns toward the dominant eigenvector, at the rate $|\lambda_2/\lambda_1|$ per step. Two ideas make this fast enough to be the industry standard. First, run power iteration on *all* $n$ columns of the identity at once, re-orthonormalizing after every step so they do not all collapse onto the same eigenvector; the orthonormalization *is* a QR factorization, and the whole thing turns out to be the innocent-looking iteration "factor $A = QR$, form $RQ$, repeat." Second, shift: the eigenvalues of $(A - \mu I)^{-1}$ are $1/(\lambda_i - \mu)$, so a shift $\mu$ near an eigenvalue makes that eigenvalue enormously dominant, and if you take $\mu$ to be the current Rayleigh quotient the convergence is *cubic*: the number of correct digits triples every step.

### Large sparse systems: never form $A^{-1}$, barely even form $A$

The 2D Poisson equation on a $200\times 200$ grid is a $40{,}000\times 40{,}000$ system with five nonzeros per row. Elimination fills in the zeros and costs $\tfrac23 n^3 \approx 4\times 10^{13}$ flops. But multiplying by $A$ costs $5n = 2\times 10^5$ flops. So the natural question is: how much can you learn about $x = A^{-1}b$ from the vectors $b, Ab, A^2b, \dots, A^{k-1}b$ for $k \ll n$? The answer is "the best approximation in that subspace," and for symmetric positive definite $A$ the conjugate gradient method finds it with a two-term recurrence, converging in about $\sqrt{\kappa}$ steps rather than $\kappa$. For the grid above that is about $600$ matrix–vector products, a third of a second, versus hours and a terabyte for elimination. The project measures exactly this.

### Randomness

If you multiply a $400\times 300$ matrix of rank $10$ by a random $300\times 10$ Gaussian matrix, the ten output columns span the column space of $A$ with probability $1$: a random vector has no reason to be orthogonal to anything in particular. If the rank is only *approximately* $10$, add a few extra random columns as insurance, and the span of the output captures the dominant subspace up to an error of a small constant times $\sigma_{11}$. That is the whole randomized SVD. The bound in §3.13 says precisely how small the constant is.

## 3. Formal development

### 3.1 IEEE double precision and the rounding model

**Definition 3.1 (IEEE 754 double precision).** A double is $x = \pm (1.b_1 b_2 \cdots b_{52})_2 \times 2^{e}$ with $52$ stored fraction bits (plus the implicit leading $1$, so a $53$-bit significand) and an exponent $e \in [-1022, 1023]$ stored in $11$ bits. The largest finite double is $\approx 1.80\times 10^{308}$, the smallest positive normal double is $\approx 2.23\times 10^{-308}$ (below that, "subnormal" numbers lose precision gradually down to $4.9\times 10^{-324}$), and the special values $\pm\infty$ and NaN absorb overflow and $0/0$. The doubles between $1$ and $2$ are spaced $2^{-52}$ apart; between $2^e$ and $2^{e+1}$ they are spaced $2^{e-52}$ apart. Rounding is to the nearest double, ties to the one with an even last bit.

Because the spacing scales with the magnitude, rounding error is *relative*:

$$\mathrm{fl}(x) = x(1 + \delta), \qquad |\delta| \le \varepsilon_{\text{mach}} = 2^{-53} \approx 1.11\times 10^{-16} \tag{3.1}$$

for every real $x$ in the normal range. (Half the spacing $2^{-52}$ relative to a number in $[1,2)$ is at most $2^{-53}$.) The IEEE standard guarantees the same for every arithmetic operation $\circ \in \{+, -, \times, \div\}$ and for $\sqrt{\cdot}$: the result is the *correctly rounded* exact result, so

$$\mathrm{fl}(x \circ y) = (x \circ y)(1 + \delta), \qquad |\delta| \le \varepsilon_{\text{mach}}. \tag{3.2}$$

Equations (3.1) and (3.2) are the only facts about the hardware this module uses. Everything else is consequences. The project verifies (3.2) against exact rational arithmetic on $8{,}000$ random operations; the worst $|\delta|$ it finds is $1.109\times 10^{-16}$, just under the bound.

Two consequences that trip people up. Floating-point addition is commutative but *not associative*: $(10^{16} + 1) - 10^{16} = 0$ while $10^{16} + (1 - 10^{16}) = 1$. And the test `1.0 + eps == 1.0` searches for the *gap* $2^{-52}$, not $\varepsilon_{\text{mach}}$: $1 + 2^{-53}$ is a tie that rounds down to $1$, anything larger rounds up.

### 3.2 Catastrophic cancellation

Suppose $x$ and $y$ are not exact but already carry rounding errors: we hold $\tilde x = x(1 + \delta_1)$ and $\tilde y = y(1+\delta_2)$ with $|\delta_i| \le \varepsilon$ (write $\varepsilon$ for $\varepsilon_{\text{mach}}$ from here on). How accurate is $\tilde x - \tilde y$ as an estimate of $x - y$?

▶ **Try it first.** Expand $\tilde x - \tilde y - (x - y)$ and bound it by the triangle inequality. Then divide by $|x - y|$.

**Derivation.**

1. $\tilde x - \tilde y = (x - y) + x\delta_1 - y\delta_2$. The subtraction itself we take as exact (by (3.2) it contributes one more factor $(1+\delta_3)$, which changes nothing below).
2. The absolute error is $|x\delta_1 - y\delta_2| \le (|x| + |y|)\varepsilon$.
3. The relative error is therefore bounded by
$$\frac{|(\tilde x - \tilde y) - (x - y)|}{|x - y|} \le \frac{|x| + |y|}{|x - y|}\,\varepsilon. \tag{3.3}$$
4. The factor in front of $\varepsilon$ is the **condition number of subtraction**,
$$\kappa_{\text{sub}}(x, y) = \frac{|x| + |y|}{|x - y|}, \tag{3.4}$$
which is about $1$ when $x$ and $y$ have opposite signs and unbounded as $y \to x$. ∎

■ **What just happened.** Subtraction is a perfectly stable *operation* (it is even exact when $x/2 \le y \le 2x$, a fact called Sterbenz's lemma). It is an ill-conditioned *problem* when the inputs are close. The errors it exposes were already present in the inputs; it merely amplifies them by (3.4). This is the conditioning/stability distinction in miniature, and you will meet it again as Theorem 3.8.

**The quadratic formula.** The roots of $ax^2 + bx + c = 0$ with $b^2 \gg 4ac$ are one large and one small. In $x_\pm = (-b \pm \sqrt{b^2 - 4ac})/2a$, one of the two signs makes $-b \pm \sqrt{\cdot}$ a sum of same-sign numbers (safe) and the other makes it a difference of nearly equal numbers (catastrophic). The fix is to compute the safe root $x_{\text{big}}$ first and get the other from the product of the roots, $x_+ x_- = c/a$:

$$x_{\text{big}} = \frac{-b - \operatorname{sign}(b)\sqrt{b^2 - 4ac}}{2a}, \qquad x_{\text{small}} = \frac{c}{a\, x_{\text{big}}}. \tag{3.5}$$

Example 4.1 runs this in 8-digit arithmetic: the naive formula returns $-0.0001$ for a root whose true value is $-0.00011111111$, one correct digit out of eight, and (3.5) returns all eight.

### 3.3 Errors, conditioning, backward stability

**Definition 3.2 (Errors).** Let $f: X \to Y$ be the problem (a map between normed spaces) and $\tilde f$ the algorithm (what the computer actually returns). For an input $x$:
- The **absolute error** is $\|\tilde f(x) - f(x)\|$; the **relative error**, which is what we always mean, is $\|\tilde f(x) - f(x)\|/\|f(x)\|$. This is also called the **forward error**.
- The **backward error** is the smallest $\|\tilde x - x\|/\|x\|$ over all $\tilde x$ with $f(\tilde x) = \tilde f(x)$: how far you must move the *input* for the computed output to become exactly right.
- $\tilde f$ is **backward stable** if for every $x$ the backward error is $O(\varepsilon)$: $\tilde f(x) = f(\tilde x)$ for some $\tilde x$ with $\|\tilde x - x\|/\|x\| \le C\varepsilon$, where $C$ may depend on the dimension but not on $x$.

Read the last line again, slowly. It does not say $\tilde f(x) \approx f(x)$. It says $\tilde f(x)$ is *exactly* $f$ of something *near* $x$. For $Ax = b$, a backward stable solver returns $\hat x$ with $(A + \delta A)\hat x = b$ and $\|\delta A\|/\|A\| = O(\varepsilon)$. You can *measure* this without knowing $x$: compute the residual $r = b - A\hat x$; then $\hat x$ solves $(A + \delta A)\hat x = b$ with $\delta A = r\hat x^\top/\|\hat x\|^2$, which has $\|\delta A\| = \|r\|/\|\hat x\|$. (Check: $\delta A\,\hat x = r\,(\hat x^\top \hat x)/\|\hat x\|^2 = r$.) So a small *relative residual* $\|r\|/(\|A\|\|\hat x\|)$ is exactly a small backward error. This is the one diagnostic you can always run.

**Definition 3.3 (Condition number).** The (relative) condition number of $f$ at $x$ is

$$\kappa(f; x) = \lim_{\delta \to 0}\ \sup_{\|\delta x\| \le \delta}\ \frac{\|f(x + \delta x) - f(x)\| / \|f(x)\|}{\|\delta x\| / \|x\|}: \tag{3.6}$$

the worst-case ratio of relative output change to relative input change, for infinitesimal changes. If $f$ is differentiable with Jacobian $J$ (LA-14), then $\kappa = \|J\|\,\|x\|/\|f(x)\|$.

**Proposition 3.4 (Matrix–vector product).** For $f(x) = Ax$ with fixed $A$ and $Ax \ne 0$,

$$\kappa(f; x) = \frac{\|A\|\,\|x\|}{\|Ax\|}, \qquad\text{and if $A$ is invertible,}\qquad \kappa(f; x) \le \|A\|\,\|A^{-1}\|.$$

*Proof.* The map is linear, so $f(x + \delta x) - f(x) = A\,\delta x$ exactly, and the ratio in (3.6) is $\dfrac{\|A\,\delta x\|/\|Ax\|}{\|\delta x\|/\|x\|} = \dfrac{\|A\,\delta x\|}{\|\delta x\|}\cdot\dfrac{\|x\|}{\|Ax\|}$. The supremum of the first factor over $\delta x$ is the definition of the operator norm $\|A\|$ (LA-9 Theorem 3.5), which gives the equality. If $A$ is invertible then $\|x\| = \|A^{-1}Ax\| \le \|A^{-1}\|\,\|Ax\|$, so $\|x\|/\|Ax\| \le \|A^{-1}\|$ and the bound follows. The bound is attained when $x$ is a right singular vector for $\sigma_n$: then $\|Ax\| = \sigma_n\|x\|$ and $\kappa = \sigma_1/\sigma_n$. ∎

That product $\|A\|\|A^{-1}\|$ is about to become the most important number in the module.

**Theorem 3.5 (Conditioning of $Ax = b$, perturbed right-hand side).** Let $A$ be invertible, $Ax = b$ and $A(x + \delta x) = b + \delta b$. Then

$$\frac{\|\delta x\|}{\|x\|} \le \kappa(A)\, \frac{\|\delta b\|}{\|b\|}, \qquad \kappa(A) = \|A\|\,\|A^{-1}\|. \tag{3.7}$$

▶ **Try it first.** Subtract the two equations to get $\delta x$ in terms of $\delta b$. Bound $\|\delta x\|$ from above and $\|b\|$ from above; combine.

**Derivation.**

1. Subtracting, $A\,\delta x = \delta b$, so $\delta x = A^{-1}\delta b$ and $\|\delta x\| \le \|A^{-1}\|\,\|\delta b\|$.
2. From $b = Ax$, $\|b\| \le \|A\|\,\|x\|$, i.e. $1/\|x\| \le \|A\|/\|b\|$.
3. Multiply: $\dfrac{\|\delta x\|}{\|x\|} \le \|A^{-1}\|\,\|\delta b\| \cdot \dfrac{\|A\|}{\|b\|} = \kappa(A)\dfrac{\|\delta b\|}{\|b\|}$. ∎

■ **What just happened.** The two inequalities are each attainable: step 1 is tight when $\delta b$ is along the direction $A^{-1}$ stretches most (the smallest singular direction, $v_n$), and step 2 is tight when $x$ is along the direction $A$ stretches most. Example 4.2 hits both at once and shows the bound is exact, then shows a $\delta b$ for which it is loose by $10^4$. The condition number is a *worst case*.

**Theorem 3.6 (Perturbed matrix).** Let $Ax = b$ and $(A + \delta A)\tilde x = b$ with $A$ invertible. Then

$$\frac{\|\tilde x - x\|}{\|\tilde x\|} \le \kappa(A)\,\frac{\|\delta A\|}{\|A\|}. \tag{3.8}$$

*Proof.* Subtract: $A(\tilde x - x) + \delta A\,\tilde x = 0$, so $\tilde x - x = -A^{-1}\delta A\,\tilde x$ and $\|\tilde x - x\| \le \|A^{-1}\|\|\delta A\|\|\tilde x\|$. Divide by $\|\tilde x\|$ and write $\|A^{-1}\|\|\delta A\| = \kappa(A)\|\delta A\|/\|A\|$. ∎

(The left side is relative to $\tilde x$ rather than $x$; to first order in $\delta A$ the two agree. There is also a combined version with both $\delta A$ and $\delta b$, giving $\kappa(A)(\|\delta A\|/\|A\| + \|\delta b\|/\|b\|)$ to first order; Exercise B2.)

**Theorem 3.7 ($\kappa$ in the 2-norm).** $\kappa_2(A) = \sigma_1/\sigma_n$, the ratio of the largest to the smallest singular value, and $\kappa_2(A) \ge 1$ with equality if and only if $A$ is a multiple of an orthogonal matrix.

*Proof.* From LA-9, $\|A\|_2 = \sigma_1$ and, since $A^{-1} = V\Sigma^{-1}U^\top$ has singular values $1/\sigma_i$, $\|A^{-1}\|_2 = 1/\sigma_n$. Equality $\sigma_1 = \sigma_n$ means $\Sigma = \sigma_1 I$, so $A = \sigma_1 UV^\top$. ∎

For symmetric $A$, $\sigma_i = |\lambda_i|$ and $\kappa_2 = |\lambda|_{\max}/|\lambda|_{\min}$. Orthogonal matrices have $\kappa = 1$: they are perfectly conditioned, which is why the whole of §3.5 is built from them.

**Theorem 3.8 (The master principle: forward error $\lesssim$ condition number $\times$ backward error).** If $\tilde f$ is a backward stable algorithm for a problem $f$ with condition number $\kappa$ at $x$, then

$$\frac{\|\tilde f(x) - f(x)\|}{\|f(x)\|} = O(\kappa(f; x)\,\varepsilon). \tag{3.9}$$

▶ **Try it first.** Write $\tilde f(x) = f(\tilde x)$ (that is what backward stable means) and apply the definition (3.6) of $\kappa$ to the pair $x, \tilde x$.

**Derivation.**

1. By backward stability there is $\tilde x$ with $\tilde f(x) = f(\tilde x)$ and $\|\tilde x - x\|/\|x\| \le C\varepsilon$.
2. By (3.6), for $\|\tilde x - x\|$ small enough, $\dfrac{\|f(\tilde x) - f(x)\|}{\|f(x)\|} \le (\kappa + o(1))\dfrac{\|\tilde x - x\|}{\|x\|}$.
3. Substitute step 1 into step 2: $\dfrac{\|\tilde f(x) - f(x)\|}{\|f(x)\|} \le (\kappa + o(1))\,C\varepsilon = O(\kappa\varepsilon)$. ∎

■ **What just happened.** Three lines, and they organize everything that follows. To know how many digits you will get, you need two numbers: the backward error of your algorithm (a property of the code; for the good algorithms in this module it is $\approx \varepsilon$) and the condition number of your problem (a property of the data; you must estimate it). The project runs `np.linalg.solve` on matrices with $\kappa$ from $1$ to $10^{12}$: the backward error stays between $5\times 10^{-17}$ and $2\times 10^{-15}$ throughout, while the forward error climbs from $10^{-15}$ to $5\times 10^{-6}$, tracking $\kappa\varepsilon$ exactly as (3.9) says.

### 3.4 Gaussian elimination: pivoting and the growth factor

LA-2 gave you $A = LU$ by elimination, and LA-10 told you to pivot. Here is why, quantitatively. Without pivoting, elimination on $\begin{pmatrix} 10^{-20} & 1 \\ 1 & 1\end{pmatrix}$ computes the multiplier $10^{20}$ and $u_{22} = 1 - 10^{20} \cdot 1 = \mathrm{fl}(1 - 10^{20}) = -10^{20}$: the $1$ in $a_{22}$ has been rounded away entirely, and the factors satisfy $LU = \begin{pmatrix} 10^{-20} & 1 \\ 1 & 0 \end{pmatrix}$, a matrix that differs from $A$ by $1$ in the $(2,2)$ entry: backward error $O(1)$, not $O(\varepsilon)$. Swapping the rows first (partial pivoting: always eliminate with the largest entry in the column) gives multiplier $10^{-20}$ and $u_{22} = 1 - 10^{-20} = 1$, and the backward error is $O(\varepsilon)$.

What went wrong in the first case is that the entries of $U$ became enormous compared with those of $A$. Define the **growth factor**

$$\rho = \frac{\max_{i,j} |u_{ij}|}{\max_{i,j} |a_{ij}|}. \tag{3.10}$$

**Theorem 3.9 (Wilkinson, 1961; stated, not proved).** Gaussian elimination with partial pivoting computes $\tilde L, \tilde U$ with $\tilde L\tilde U = PA + \delta A$ and $\|\delta A\|/\|A\| = O(\rho\,\varepsilon)$ (the constant depends polynomially on $n$). With partial pivoting every multiplier satisfies $|l_{ij}| \le 1$, which implies $\rho \le 2^{n-1}$.

The proof is a careful accounting of (3.2) through the elimination and lives in Higham, *Accuracy and Stability*, §9.3, and Trefethen & Bau lecture 22. The content: pivoting is backward stable *if the growth factor is modest*. And here you must be careful about what is true. The bound $2^{n-1}$ is attained, by the matrix with $1$ on the diagonal, $-1$ below it and $1$ in the last column (the project checks $\rho = 512$ for $n = 10$ and $524{,}288$ for $n = 20$). So Gaussian elimination with partial pivoting is **not** backward stable in the worst case. What is true is that for matrices arising in practice the growth is small, typically of order $n^{1/2}$ to $n^{2/3}$ (the project measures $\rho \approx 4$ on a random $100\times 100$ Gaussian matrix, against a worst case of $6\times 10^{29}$), and that in sixty years of use no naturally occurring matrix with exponential growth has turned up. Why the worst case is so rare is only partly understood (Trefethen and Schreiber, 1990, give the random-matrix statistics). Say "backward stable in practice" and mean exactly that.

### 3.5 Householder reflections and Householder QR

**Theorem 3.10 (Householder reflector).** For any nonzero $v \in \mathbb{R}^m$, the matrix

$$H = I - 2\frac{vv^\top}{v^\top v} \tag{3.11}$$

is symmetric, orthogonal, and involutory ($H^2 = I$). It is the reflection across the hyperplane $v^\perp$: $Hv = -v$ and $Hw = w$ for every $w \perp v$.

▶ **Try it first.** Compute $H^\top$ and $H^2$ directly from (3.11); $(vv^\top)(vv^\top) = v(v^\top v)v^\top$ is the only identity you need.

**Derivation.**

1. Symmetric: $(vv^\top)^\top = vv^\top$, so $H^\top = H$.
2. Involutory: $H^2 = I - 4\dfrac{vv^\top}{v^\top v} + 4\dfrac{v(v^\top v)v^\top}{(v^\top v)^2} = I - 4\dfrac{vv^\top}{v^\top v} + 4\dfrac{vv^\top}{v^\top v} = I$.
3. Orthogonal: $H^\top H = H^2 = I$.
4. Action: $Hv = v - 2v\dfrac{v^\top v}{v^\top v} = -v$; if $v^\top w = 0$ then $Hw = w - 2v\cdot 0 = w$. So $H$ fixes the hyperplane $v^\perp$ pointwise and flips its normal: a reflection. ∎

The use of $H$ is to zero out all but one entry of a vector with a single orthogonal transformation. We want $Hx = \alpha\mathbf{e}_1$; since $H$ preserves length, $|\alpha| = \|x\|$. The reflection that maps $x$ to $\|x\|\mathbf{e}_1$ must reflect across the hyperplane that bisects them, whose normal is their difference:

$$v = x - \|x\|\mathbf{e}_1 \quad\Longrightarrow\quad Hx = \|x\|\mathbf{e}_1. \tag{3.12}$$

▶ **Try it first.** Verify (3.12): compute $v^\top x$ and $v^\top v$ and show $v^\top v = 2v^\top x$.

**Derivation.**

1. $v^\top x = \|x\|^2 - \|x\|x_1$.
2. $v^\top v = \|x\|^2 - 2\|x\|x_1 + \|x\|^2 = 2(\|x\|^2 - \|x\|x_1) = 2\,v^\top x$.
3. $Hx = x - 2v\dfrac{v^\top x}{v^\top v} = x - 2v\cdot\dfrac12 = x - v = \|x\|\mathbf{e}_1$. ∎

**The sign.** Replacing $\|x\|\mathbf{e}_1$ by $-\|x\|\mathbf{e}_1$ gives an equally valid reflection with $v = x + \|x\|\mathbf{e}_1$. In exact arithmetic either works. In floating point only one is safe. The first component of $v$ is $x_1 \mp \|x\|$. If $x$ is nearly parallel to $\mathbf{e}_1$, then $\|x\| \approx |x_1|$ and the choice with a minus sign is a catastrophic cancellation: by (3.4) the computed $v_1$ can have relative error of order $\varepsilon\|x\|/|x_1 - \|x\||$, which can be $O(1)$. A wrong $v_1$ means a wrong *direction* for $v$, and $H$ depends only on the direction of $v$: the computed $H$ is still an exact reflection (Theorem 3.10 holds for whatever $v$ you have), but across the wrong hyperplane, so $Hx$ is not a multiple of $\mathbf{e}_1$. The project builds this failure: for $x = (1, 10^{-8}, 10^{-8})$ the bad sign leaves a residual of $1.4\times 10^{-8}$ in the entries that should be zero; the good sign leaves $5\times 10^{-25}$. The rule, which we adopt throughout:

$$v = x + \operatorname{sign}(x_1)\|x\|\mathbf{e}_1, \qquad Hx = -\operatorname{sign}(x_1)\|x\|\mathbf{e}_1, \tag{3.13}$$

with $\operatorname{sign}(0) = +1$. Then $v_1 = x_1 + \operatorname{sign}(x_1)\|x\|$ is a sum of two numbers of the same sign and cannot cancel.

**Algorithm 3.1 (Householder QR).** For $A$ of size $m\times n$ with $m \ge n$:

```
R = A
for k = 1, …, n:
    x = R[k:m, k]                         # the part of column k on and below the diagonal
    v = x + sign(x_1) ‖x‖ e_1 ;  v = v/‖v‖  # (3.13)
    R[k:m, k:n] = R[k:m, k:n] − 2 v (vᵀ R[k:m, k:n])   # apply H_k to the active block
    store v
Q = I_m ;  for k = n, …, 1:  Q[k:m, :] = Q[k:m, :] − 2 v_k (v_kᵀ Q[k:m, :])
```

Step $k$ applies $Q_k = \operatorname{diag}(I_{k-1}, H_k)$, which fixes the first $k-1$ rows and columns (already in triangular form) and zeros column $k$ below the diagonal. After $n$ steps $Q_n\cdots Q_1 A = R$ is upper triangular, so $A = QR$ with $Q = Q_1 Q_2\cdots Q_n$ (each $Q_k$ is its own inverse and transpose). The second loop forms $Q$ by applying $Q_1\cdots Q_n$ to $I$ from the inside out; in practice you often skip it and apply the stored reflectors to whatever you need $Q^\top b$ for.

**Cost.** At step $k$ the active block is $(m-k+1)\times(n-k+1)$. Forming $v^\top R$ costs $2(m-k)(n-k)$ flops and the rank-one update another $2(m-k)(n-k)$, so about $4(m-k)(n-k)$ per step. Summing, with $j = n - k$,
$$\sum_{k=1}^{n} 4(m-k)(n-k) \approx 4\int_0^n (m - n + j)\,j\,dj = 4\Big[(m-n)\frac{n^2}{2} + \frac{n^3}{3}\Big] = 2mn^2 - \frac{2}{3}n^3 \tag{3.14}$$
flops. For square matrices that is $\tfrac43 n^3$, twice the cost of LU; the price of orthogonality. Forming $Q$ explicitly costs about the same again.

**Theorem 3.11 (Backward stability of Householder QR; Trefethen & Bau Theorem 16.1, stated).** Let $\tilde Q, \tilde R$ be the factors computed by Algorithm 3.1 in floating point, where $\tilde Q$ denotes the *exactly* orthogonal matrix $Q_1\cdots Q_n$ built from the computed reflector vectors. Then $\tilde Q\tilde R = A + \delta A$ with $\|\delta A\|/\|A\| = O(\varepsilon)$.

The proof (T&B lecture 16; Higham ch. 19) is an induction using (3.2) and one lemma: applying a computed Householder reflector to a vector is backward stable, $\mathrm{fl}(\tilde H x) = \tilde H(x + \delta x)$ with $\|\delta x\| = O(\varepsilon\|x\|)$. The idea you should retain is *why* this is true and Gram–Schmidt is not. Each Householder step applies an orthogonal matrix, which neither grows nor shrinks any error already present; new rounding errors enter at relative size $\varepsilon$ of the *current* block, and they stay that size. Gram–Schmidt instead computes $q_j$ by subtracting from $a_j$ its projections onto $q_1, \dots, q_{j-1}$; when $a_j$ is nearly in their span this is a difference of nearly equal vectors, and by (3.4) the computed $q_j$ has relative error of order $\kappa\varepsilon$ in its *direction*, which is exactly a loss of orthogonality. Quantitatively (Björck 1967; Higham §19.8): classical Gram–Schmidt loses orthogonality like $\|Q^\top Q - I\| \approx \kappa(A)^2\varepsilon$, modified Gram–Schmidt like $\kappa(A)\varepsilon$, Householder like $\varepsilon$ regardless of $\kappa$. The project measures all three on Hilbert matrices (LA-10's standard ill-conditioned example): at $n = 10$, $\kappa = 1.6\times 10^{13}$, the three numbers are $3.5$, $2.5\times 10^{-4}$ and $2.2\times 10^{-15}$. Classical Gram–Schmidt has lost *all* orthogonality by $n = 7$. Note also that both Gram–Schmidt variants still satisfy $QR \approx A$ to $\varepsilon$; what they lose is the orthogonality of $Q$, which is the property you were paying for.

### 3.6 Least squares: normal equations versus QR

For $A$ of size $m\times n$ with $m > n$ and full column rank, the least-squares solution $\hat x = \arg\min\|Ax - b\|$ satisfies the normal equations $A^\top A\hat x = A^\top b$ (LA-5). The matrix $A^\top A$ is SPD, and solving by Cholesky costs $mn^2 + n^3/3$, about half of QR. It is nonetheless the wrong method.

**Theorem 3.12.** $\kappa_2(A^\top A) = \kappa_2(A)^2$.

▶ **Try it first.** Write $A = U\Sigma V^\top$ (LA-9, thin form) and compute $A^\top A$.

**Derivation.**

1. $A^\top A = V\Sigma^\top U^\top U\Sigma V^\top = V\Sigma^\top\Sigma V^\top = V\operatorname{diag}(\sigma_1^2, \dots, \sigma_n^2)V^\top$, using $U^\top U = I$.
2. This is a spectral (and singular value) decomposition of $A^\top A$: its singular values are $\sigma_i^2$.
3. By Theorem 3.7, $\kappa_2(A^\top A) = \sigma_1^2/\sigma_n^2 = \kappa_2(A)^2$. ∎

■ **What just happened.** Forming $A^\top A$ squares the condition number *before* any algorithm runs. By (3.9), a backward stable solve of the normal equations gives forward error $\kappa(A)^2\varepsilon$. QR instead solves $R\hat x = Q^\top b$: $Q$ is orthogonal so $\kappa(R) = \kappa(A)$ (same singular values), and Householder QR followed by back substitution is backward stable for least squares (T&B Theorem 19.3), so the error is of order $\kappa(A)\varepsilon$ when the residual is small. With $\kappa(A) = 10^8$ the normal equations are at $\kappa^2\varepsilon = 1$: no correct digits. The project: error $6\times 10^{-3}$ for the normal equations, $7.5\times 10^{-10}$ for Householder QR, same problem. (A caveat for §8: when the residual $\|b - A\hat x\|$ is *large*, the least-squares problem itself has a $\kappa(A)^2$ term in its conditioning, T&B Theorem 18.1, and no algorithm beats it. QR is the right method; it is not a miracle.)

### 3.7 Eigenvalues: never the characteristic polynomial, so iterate well

LA-7 taught you to find eigenvalues as roots of $p_A(\lambda) = \det(A - \lambda I)$. Every textbook before 1960 computed them that way. Here is why nobody does now.

**The roots of a polynomial are an ill-conditioned function of its coefficients.** Let $p(x) = \sum_{j=0}^{n} a_j x^j$ have a simple root $r$ (so $p(r) = 0$, $p'(r) \ne 0$) and perturb one coefficient, $a_j \to a_j + \delta a_j$.

▶ **Try it first.** The perturbed root $r + \delta r$ satisfies $\sum_j a_j (r+\delta r)^j + \delta a_j\,r^j = 0$ to first order. Expand and solve for $\delta r$.

**Derivation.**

1. To first order in small quantities, $p(r + \delta r) = p(r) + p'(r)\,\delta r = p'(r)\,\delta r$, since $p(r) = 0$.
2. The perturbed polynomial at the perturbed root is $p(r + \delta r) + \delta a_j\,(r + \delta r)^j = 0$. To first order, $p'(r)\,\delta r + \delta a_j\,r^j = 0$.
3. Hence
$$\delta r \approx -\frac{r^j}{p'(r)}\,\delta a_j, \qquad \frac{|\delta r|/|r|}{|\delta a_j|/|a_j|} \approx \frac{|a_j|\,|r|^{j-1}}{|p'(r)|}. \tag{3.15}$$
The right-hand quantity is the relative condition number of the root $r$ with respect to the coefficient $a_j$, in the sense of (3.6). ∎

■ **What just happened.** The sensitivity is $|r|^j$ divided by $|p'(r)| = \prod_{i \ne r}|r - r_i|$. For roots that are large and close to their neighbours, the numerator is huge and the denominator modest. For **Wilkinson's polynomial** $w(x) = \prod_{k=1}^{20}(x - k)$, the coefficient of $x^{19}$ is $a_{19} = -210$, and at the root $r = 16$ we have $p'(16) = 15!\cdot 4! \approx 3.1\times 10^{13}$ while $16^{19} \approx 7.6\times 10^{22}$. So $|r^{19}/p'(r)| \approx 2.4\times 10^{9}$ and the relative condition number (3.15) is about $3.2\times 10^{10}$: a relative change of $10^{-10}$ in one coefficient moves the root by a relative amount of order $1$. The perturbation Wilkinson used, $\delta a_{19} = -2^{-23} \approx -1.2\times 10^{-7}$, is so far outside the linear regime that (3.15) predicts a move of $287$, which is meaningless; what actually happens is that the roots $10$ through $19$ collide in pairs and leave the real axis, the pair near $16.73$ acquiring imaginary parts $\pm 2.81$. The root $r = 1$, by contrast, has condition number $1.7\times 10^{-15}$ and does not move. The project computes all twenty perturbed roots to 30 digits and plots them. Worse: the coefficients of $w$ run up to $20! \approx 2.4\times 10^{18} > 2^{53}$, so they cannot even be *stored* exactly in double precision, and `np.roots` on the unperturbed coefficients already returns roots off by $0.085$. Forming $p_A$ throws away the matrix and keeps twenty numbers that do not determine its eigenvalues to any useful accuracy. (For a symmetric matrix the eigenvalues themselves are perfectly conditioned, $|\delta\lambda_i| \le \|\delta A\|$ by Weyl's inequality, LA-8 Exercise C1; it is the *polynomial* that is the problem, not the eigenvalues.)

**So every eigenvalue algorithm is iterative.** Abel (1824) and Ruffini proved there is no formula in radicals for the roots of a general polynomial of degree $5$ or more. Since the eigenvalues of any $n\times n$ matrix are the roots of its characteristic polynomial, and every monic polynomial is the characteristic polynomial of its companion matrix, a finite exact algorithm for eigenvalues would be a finite exact algorithm for polynomial roots, which does not exist. This is the opposite of the situation for $Ax = b$ or $A = QR$, which finite algorithms solve exactly in exact arithmetic. Eigenvalue algorithms must converge; the only questions are "to what" and "how fast."

**Theorem 3.13 (Power iteration; LA-7 §3.6).** Let $A$ be diagonalizable with eigenvalues $|\lambda_1| > |\lambda_2| \ge \cdots \ge |\lambda_n|$ and eigenvectors $x_1, \dots, x_n$. If $v_0 = \sum c_i x_i$ with $c_1 \ne 0$, then $v_k = A^k v_0/\|A^k v_0\|$ satisfies $\|v_k - (\pm x_1/\|x_1\|)\| = O(|\lambda_2/\lambda_1|^k)$. If $A$ is also symmetric, the Rayleigh quotient $\mu_k = v_k^\top A v_k$ satisfies $|\mu_k - \lambda_1| = O(|\lambda_2/\lambda_1|^{2k})$.

*Derivation.* $A^k v_0 = \sum_i c_i\lambda_i^k x_i = c_1\lambda_1^k\big(x_1 + \sum_{i\ge 2}(c_i/c_1)(\lambda_i/\lambda_1)^k x_i\big)$, and every term in the sum is $O(|\lambda_2/\lambda_1|^k)$. For the Rayleigh quotient, write $v_k = (x_1 + \varepsilon w)/\|x_1 + \varepsilon w\|$ with $w \perp x_1$ a unit vector (possible when $A$ is symmetric, since the eigenvectors are orthogonal) and $\varepsilon = O(|\lambda_2/\lambda_1|^k)$. Then $v_k^\top A v_k = (\lambda_1 + \varepsilon^2 w^\top A w)/(1 + \varepsilon^2)$, which differs from $\lambda_1$ by $O(\varepsilon^2)$: the first-order term vanishes because $x_1^\top A w = \lambda_1 x_1^\top w = 0$. ∎

The Rayleigh quotient is *quadratically* accurate in the eigenvector error, for the same reason a function is flat at a critical point: the eigenvectors are the critical points of $v \mapsto v^\top A v/v^\top v$ (LA-8 Theorem 3.17). Keep this; it is half of the cubic convergence below.

**Inverse iteration with a shift.** Power iteration finds only $\lambda_1$ and is slow when $|\lambda_2/\lambda_1|$ is near $1$. Both problems have one fix.

▶ **Try it first.** If $Ax = \lambda x$, what does $(A - \mu I)^{-1}$ do to $x$?

**Derivation.** From $Ax = \lambda x$, $(A - \mu I)x = (\lambda - \mu)x$, so for $\mu$ not an eigenvalue,
$$(A - \mu I)^{-1}x = \frac{1}{\lambda - \mu}\,x. \tag{3.16}$$
The eigenvectors are unchanged and the eigenvalues become $1/(\lambda_i - \mu)$. Power iteration on $(A - \mu I)^{-1}$ therefore converges to the eigenvector whose eigenvalue is *closest to $\mu$*, at the rate $|\lambda_{\text{closest}} - \mu|/|\lambda_{\text{next}} - \mu|$ per step, which is tiny when $\mu$ is a good estimate. Each step is a linear solve with $A - \mu I$; factor it once ($\tfrac23 n^3$) and each step costs $2n^2$. ∎

■ **What just happened.** With eigenvalues $1, 2, 3$ and the shift $\mu = 2.9$, the transformed eigenvalues are $-0.53, -1.11, 10$: the ratio drops from $2/3$ (power iteration) to $0.11$, and reaching $10^{-10}$ takes $11$ steps instead of $57$. Note the apparent paradox: $A - \mu I$ is nearly singular when $\mu$ is nearly an eigenvalue, so the solve is ill-conditioned (Theorem 3.5), and yet the method works. The error of the solve is large but lies almost entirely along the eigenvector you want, so after normalizing it does no harm (Trefethen & Bau lecture 27, Exercise 27.5). Ill-conditioning that is aligned with what you are computing is harmless; this is a useful thing to know.

**Algorithm 3.2 (Rayleigh quotient iteration).** Start with a unit vector $v_0$. For $k = 0, 1, 2, \dots$: $\mu_k = v_k^\top A v_k$; solve $(A - \mu_k I)w = v_k$; $v_{k+1} = w/\|w\|$.

**Theorem 3.14 (Cubic convergence of RQI for symmetric $A$; T&B Theorem 27.3, sketched).** If $v_k$ converges to an eigenvector $x_J$ of a symmetric $A$ with simple eigenvalue $\lambda_J$, then eventually $\|v_{k+1} - (\pm x_J)\| = O(\|v_k - (\pm x_J)\|^3)$ and $|\mu_{k+1} - \lambda_J| = O(|\mu_k - \lambda_J|^3)$.

*Sketch (this is the whole idea; the constants are in T&B).* Let the eigenvector error at step $k$ be $\varepsilon$. By Theorem 3.13 the Rayleigh quotient shift is accurate to $\varepsilon^2$: $|\mu_k - \lambda_J| = O(\varepsilon^2)$. One step of inverse iteration with that shift multiplies the eigenvector error by $|\lambda_J - \mu_k|/|\lambda_{\text{next}} - \mu_k| = O(\varepsilon^2)$. So the new error is $O(\varepsilon\cdot\varepsilon^2) = O(\varepsilon^3)$. ∎

Cubic means the number of correct digits roughly triples each step. The project's run on a $6\times 6$ matrix gives eigenvalue errors $1.7\times 10^{-1}$, $2.3\times 10^{-3}$, $4.9\times 10^{-9}$, $8.9\times 10^{-16}$: $0, 2, 8, 15$ digits. Three steps from a rough guess to machine precision. The catch is the same as for inverse iteration: each step is a fresh factorization, since the shift changes, so RQI is for polishing one eigenpair, not for finding all $n$. For that we need §3.8.

### 3.8 The QR algorithm

The algorithm that `np.linalg.eig` and `eigh` run was found by John Francis (1961) and Vera Kublanovskaya (1961) independently. In its purest form it is one line:

$$A_0 = A; \qquad A_k = Q_k R_k \ \text{(QR factorization)}, \qquad A_{k+1} = R_k Q_k. \tag{3.17}$$

Factor, multiply back in the opposite order, repeat. Two theorems explain why this innocent-looking iteration computes eigenvalues.

**Theorem 3.15 (Every iterate is orthogonally similar to $A$).** With $Q^{(k)} = Q_1 Q_2\cdots Q_k$,
$$A_{k+1} = Q_k^\top A_k Q_k = (Q^{(k)})^\top A\, Q^{(k)}.$$
In particular every $A_k$ has the same eigenvalues as $A$, and if $A$ is symmetric so is every $A_k$.

▶ **Try it first.** From $A_k = Q_k R_k$, solve for $R_k$ and substitute into $A_{k+1} = R_k Q_k$.

**Derivation.**

1. $Q_k$ is orthogonal, so $R_k = Q_k^\top A_k$.
2. $A_{k+1} = R_k Q_k = Q_k^\top A_k Q_k$.
3. Induct: $A_{k+1} = Q_k^\top Q_{k-1}^\top \cdots Q_1^\top A\, Q_1\cdots Q_k = (Q^{(k)})^\top A Q^{(k)}$.
4. Orthogonal similarity preserves eigenvalues (LA-7 Theorem 3.18) and symmetry: $(Q^\top A Q)^\top = Q^\top A^\top Q$. ∎

**Theorem 3.16 (The QR algorithm is simultaneous power iteration).** With $R^{(k)} = R_k R_{k-1}\cdots R_1$ (upper triangular, as a product of upper triangular matrices),
$$A^k = Q^{(k)} R^{(k)}.$$
That is: $Q^{(k)}$ is the orthogonal factor in the QR factorization of $A^k$.

▶ **Try it first.** Check $k = 1$ (trivial) and $k = 2$: write $A^2 = A\cdot A = Q_1 R_1 Q_1 R_1$ and use $R_1 Q_1 = A_1 = Q_2 R_2$.

**Derivation (induction on $k$).**

1. $k = 1$: $A = Q_1 R_1$ by definition.
2. Assume $A^{k-1} = Q^{(k-1)}R^{(k-1)}$. Then $A^k = A\,A^{k-1} = A\,Q^{(k-1)} R^{(k-1)}$.
3. By Theorem 3.15, $A\,Q^{(k-1)} = Q^{(k-1)} A_{k}$ (move $(Q^{(k-1)})^\top$ to the other side of $A_k = (Q^{(k-1)})^\top A Q^{(k-1)}$). And $A_k = Q_k R_k$. So $A\,Q^{(k-1)} = Q^{(k-1)} Q_k R_k = Q^{(k)} R_k$.
4. Substitute: $A^k = Q^{(k)} R_k R^{(k-1)} = Q^{(k)} R^{(k)}$. ∎

■ **What just happened.** Read the first column of $A^k = Q^{(k)}R^{(k)}$: $A^k e_1 = r^{(k)}_{11}\, q^{(k)}_1$. So the first column of $Q^{(k)}$ is the normalized $A^k e_1$: **power iteration**, converging to the dominant eigenvector at rate $|\lambda_2/\lambda_1|$. Read the first $j$ columns: $A^k[e_1, \dots, e_j] = Q^{(k)}_{:,1:j}\, R^{(k)}_{1:j,1:j}$, so the first $j$ columns of $Q^{(k)}$ are an orthonormal basis for $A^k\operatorname{span}\{e_1,\dots,e_j\}$. That is power iteration on $j$ vectors at once, re-orthonormalized at every step so they do not all collapse onto $x_1$: **simultaneous iteration**, and the $j$-dimensional space converges to the span of the top $j$ eigenvectors at rate $|\lambda_{j+1}/\lambda_j|$ (T&B Theorem 28.4 and the argument of Theorem 3.13 applied blockwise). Now use Theorem 3.15: $A_{k+1} = (Q^{(k)})^\top A Q^{(k)}$ is $A$ written in the basis $Q^{(k)}$. When the first $j$ columns of $Q^{(k)}$ span an invariant subspace of $A$, the block of $A_{k+1}$ below the first $j$ columns is zero. So $A_k$ converges to upper triangular form, the **Schur form** of LA-10 Theorem 3.8, with the eigenvalues on the diagonal in decreasing order of modulus, the $(j+1, j)$ entry decaying like $|\lambda_{j+1}/\lambda_j|^k$. For symmetric $A$, triangular and symmetric means diagonal, and $Q^{(k)}$ converges to the eigenvector matrix. Example 4.4 does two steps by hand on a $2\times 2$ and watches the off-diagonal shrink by $\lambda_2/\lambda_1 = 3 - 2\sqrt2 \approx 0.172$ per step; the project runs a $6\times 6$ with eigenvalues $1, \dots, 6$ and confirms that the five subdiagonal entries decay at rates $5/6, 4/5, 3/4, 2/3, 1/2$, each at its own speed.

(What the argument assumes: $|\lambda_1| > |\lambda_2| > \cdots > |\lambda_n|$ and a genericity condition on the leading principal minors of the eigenvector matrix. Equal-modulus pairs, including complex conjugate pairs of a real matrix, do not separate: the real QR algorithm converges to the *real* Schur form of LA-10 Theorem 3.11, with $2\times 2$ blocks for complex pairs. We stay with symmetric matrices for the practical algorithm.)

**Making it practical.** As written, (3.17) is hopeless: each step is a full QR at $\tfrac43 n^3$ flops, and convergence at rate $|\lambda_{j+1}/\lambda_j|$ is slow when eigenvalues are close. Three ideas fix this, and together they are the algorithm in LAPACK.

**(a) Reduce to Hessenberg (tridiagonal) form first.** An upper Hessenberg matrix has $h_{ij} = 0$ for $i > j+1$: upper triangular plus one subdiagonal. Any $A$ can be brought to Hessenberg form by $n-2$ Householder reflections applied from *both* sides, $H = Q^\top A Q$: at step $k$ choose the reflector that zeros column $k$ below row $k+1$ (not below row $k$: a reflector that acts on rows $k+1\ldots n$ leaves row $k$ alone, so when you also apply it from the right, which mixes columns $k+1\ldots n$, the zeros already created in column $k$ survive). If you tried to zero column $k$ below row $k$, the right-multiplication would mix columns $k\ldots n$ and destroy them; you cannot reach triangular form this way, which is just as well, since triangular form *is* the Schur form and would mean a finite eigenvalue algorithm. For symmetric $A$ the Hessenberg form is symmetric, hence **tridiagonal**. Cost: $\tfrac{10}{3}n^3$ in general, $\tfrac43 n^3$ for the symmetric case (T&B lecture 26).

**Theorem 3.17 (A QR step preserves Hessenberg form).** If $A_k$ is upper Hessenberg and invertible, so is $A_{k+1} = R_k Q_k$.

*Proof.* $Q_k = A_k R_k^{-1}$. The inverse of an upper triangular matrix is upper triangular (LA-2), and (Hessenberg)$\times$(upper triangular) is Hessenberg: entry $(i,j)$ of the product is $\sum_l h_{il}u_{lj}$, and $h_{il} \ne 0$ needs $l \ge i-1$ while $u_{lj} \ne 0$ needs $l \le j$, so the entry vanishes when $i - 1 > j$. So $Q_k$ is Hessenberg, and $R_k Q_k$ = (upper triangular)$\times$(Hessenberg) is Hessenberg by the same count. ∎

Why this matters: the QR factorization of a Hessenberg matrix needs only $n-1$ **Givens rotations**, one per subdiagonal entry. A Givens rotation in the $(k, k+1)$ plane is

$$G_k = \begin{pmatrix} c & s \\ -s & c\end{pmatrix} \text{ in rows/columns } k, k+1, \qquad c = \frac{a}{\sqrt{a^2 + b^2}},\ s = \frac{b}{\sqrt{a^2+b^2}} \ \text{ zeros } b \text{ below } a. \tag{3.18}$$

Each rotation touches two rows, so one QR step on a Hessenberg matrix costs $O(n^2)$ instead of $O(n^3)$, and on a tridiagonal matrix $O(n)$. The project's `qr_step_tridiagonal` is exactly this: $n-1$ rotations from the left to get $R$, then the same rotations from the right to form $RQ$.

**(b) Shift.** Replace (3.17) by $A_k - \mu_k I = Q_k R_k$, $A_{k+1} = R_k Q_k + \mu_k I$. Then $A_{k+1} = Q_k^\top(A_k - \mu_k I)Q_k + \mu_k I = Q_k^\top A_k Q_k$: still an orthogonal similarity, same eigenvalues. By Theorem 3.16 applied to $A - \mu I$, the *last* column of $Q^{(k)}$ is now inverse iteration with shift $\mu$ (the last column of $Q$ in $B = QR$ is orthogonal to the first $n-1$ columns of $B$, which makes it proportional to $B^{-\top}e_n$; T&B lecture 29 does the bookkeeping). Choosing $\mu_k = (A_k)_{nn}$, the Rayleigh quotient of $e_n$, makes the last row and column of $A_k$ undergo Rayleigh quotient iteration, cubically convergent by Theorem 3.14. The $(n, n-1)$ entry is cubed every step.

The Rayleigh shift can fail: for $A = \begin{pmatrix}0 & 1\\1 & 0\end{pmatrix}$, $\mu = a_{22} = 0$ is exactly midway between the eigenvalues $\pm 1$ and the iteration never moves. The repair is the **Wilkinson shift**: take $\mu$ to be the eigenvalue of the trailing $2\times 2$ block $\begin{pmatrix} a & b\\ b & c\end{pmatrix}$ closer to $c$, computed as

$$\mu = c - \frac{\operatorname{sign}(\delta)\, b^2}{|\delta| + \sqrt{\delta^2 + b^2}}, \qquad \delta = \frac{a - c}{2}, \tag{3.19}$$

with $\operatorname{sign}(0) = 1$. (Check: the eigenvalues of the block are $\frac{a+c}{2} \pm \sqrt{\delta^2 + b^2}$; the one closer to $c$ is $c + \delta - \operatorname{sign}(\delta)\sqrt{\delta^2+b^2}$, and multiplying numerator and denominator by $|\delta| + \sqrt{\delta^2+b^2}$ gives (3.19). This is the stable quadratic formula (3.5) again: the subtraction $|\delta| - \sqrt{\delta^2 + b^2}$ would cancel catastrophically when $b \ll |\delta|$, which is exactly when the algorithm is converging.) For symmetric tridiagonal matrices, QR with the Wilkinson shift always converges, at least linearly and in practice cubically (Wilkinson 1968; T&B lecture 29). For $A = \begin{pmatrix}3&1\\1&1\end{pmatrix}$ of Example 4.4 the Wilkinson shift is $2 - \sqrt2 = 0.5858$, the exact smaller eigenvalue: one shifted step would finish what the unshifted iteration takes a dozen steps to do.

**(c) Deflate.** When a subdiagonal entry $t_{j+1,j}$ becomes negligible, $|t_{j+1,j}| \le \varepsilon(|t_{jj}| + |t_{j+1,j+1}|)$, set it to zero. The matrix is now block diagonal and the two blocks have separate eigenvalue problems; recurse on each. Since the shifted iteration drives the *last* subdiagonal entry to zero cubically, the usual pattern is: a few steps, deflate one eigenvalue off the bottom, repeat on the $(n-1)\times(n-1)$ block. Empirically about two QR steps per eigenvalue; the project's `qr_algorithm_symmetric` takes $20$, $105$ and $304$ steps for $n = 10, 50, 150$ ($2.0$, $2.1$, $2.0$ per eigenvalue) and matches `eigvalsh` to $6\times 10^{-15}$ relative.

**Total cost (symmetric).** Tridiagonalization $\tfrac43 n^3$, then $O(n)$ per QR step times $O(n)$ steps: the eigenvalues cost $\tfrac43 n^3 + O(n^2)$, essentially the price of the reduction. Accumulating the eigenvectors (applying every rotation to $Q$) raises this to about $\tfrac{9}{2}n^3$ or so; LAPACK's `dsyev` (which `eigh` calls) is this algorithm with a divide-and-conquer or relatively-robust-representation variant for the tridiagonal phase. For nonsymmetric $A$ the same structure holds with the Francis double shift (two complex conjugate shifts combined into one real step, so complex arithmetic is never needed) and the real Schur form as the target: roughly $10n^3$ flops for eigenvalues, $25n^3$ with the Schur vectors (Golub & Van Loan §7.5). All of it is backward stable: the computed Schur form is the exact Schur form of $A + \delta A$ with $\|\delta A\| = O(\varepsilon\|A\|)$, so by Weyl's inequality the computed eigenvalues of a symmetric matrix have *absolute* error $O(\varepsilon\|A\|)$; small eigenvalues of a matrix with large ones carry large *relative* error, which is a property of the problem, not the algorithm.

**The SVD: Golub–Kahan.** The temptation is to compute the SVD as the eigendecomposition of $A^\top A$. Theorem 3.12 says this squares the condition number: singular values below $\sqrt{\varepsilon}\,\sigma_1 \approx 10^{-8}\sigma_1$ are lost entirely (LA-9's project shows it). Golub and Kahan (1965) found the way around: reduce $A$ to upper **bidiagonal** form $A = U B V^\top$ by Householder reflections alternating from the left (zero column $k$ below the diagonal) and the right (zero row $k$ to the right of the superdiagonal), at $4mn^2 - \tfrac43 n^3$ flops; $B$ has the same singular values as $A$ (unitary equivalence, LA-10 §3.1). Then observe that $B^\top B$ is tridiagonal, so the symmetric QR algorithm on $B^\top B$ would be cheap, and arrange the rotations so that they act on $B$ itself without ever forming $B^\top B$ (the "implicit" shifted QR step of Golub and Reinsch 1970; Demmel and Kahan 1990 refined it to compute even the tiny singular values to high relative accuracy). That is `np.linalg.svd`. The project implements the bidiagonalization and checks that $B$'s singular values match $A$'s to machine precision; the implicit QR sweep is Trefethen & Bau lecture 31 and Golub & Van Loan §8.6.

### 3.9 Large sparse systems: the model problem and stationary iterations

**The running example.** Poisson's equation $-\nabla^2 u = f$ on the unit square with $u = 0$ on the boundary is the electrostatic potential of a charge density (EM-3), the steady temperature of a heated plate, and the deflection of a membrane. Put $N$ interior grid points per side with spacing $h = 1/(N+1)$, write $u_{ij} \approx u(ih, jh)$, and replace each second derivative by its centered difference $(u_{i-1,j} - 2u_{ij} + u_{i+1,j})/h^2$. The result is the **five-point stencil**

$$\frac{4u_{ij} - u_{i-1,j} - u_{i+1,j} - u_{i,j-1} - u_{i,j+1}}{h^2} = f_{ij}, \tag{3.20}$$

one equation per interior point: $n = N^2$ equations in $N^2$ unknowns. Ordering the unknowns row by row and using the Kronecker product of LA-12,

$$A = \frac{1}{h^2}\left(I_N\otimes L_N + L_N\otimes I_N\right), \qquad L_N = \begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & \ddots & \ddots & \ddots \\ & & -1 & 2\end{pmatrix}, \tag{3.21}$$

an $N^2\times N^2$ matrix with at most five nonzeros per row. For $N = 200$ that is $40{,}000$ unknowns, $200{,}000$ nonzeros, and $1.6\times 10^9$ entries in the dense matrix you must never form (12.8 GB).

**Theorem 3.18 (Spectrum of the model problem).** The eigenvalues of $A$ in (3.21) are

$$\lambda_{jk} = \frac{4}{h^2}\left(\sin^2\frac{j\pi h}{2} + \sin^2\frac{k\pi h}{2}\right), \qquad j, k = 1, \dots, N, \tag{3.22}$$

so $A$ is SPD with $\lambda_{\min} \approx 2\pi^2$, $\lambda_{\max} \approx 8/h^2$, and $\kappa(A) = \cot^2(\pi h/2) \approx \dfrac{4}{\pi^2 h^2} \approx 0.4\,N^2$.

▶ **Try it first.** Guess that the eigenvectors of $L_N$ are the discrete sines $v^{(j)}_i = \sin(ij\pi h)$, and apply $L_N$ to one. You will need $\sin(a - b) + \sin(a+b) = 2\sin a\cos b$.

**Derivation.**

1. Row $i$ of $L_N v^{(j)}$ is $-\sin((i-1)j\pi h) + 2\sin(ij\pi h) - \sin((i+1)j\pi h) = 2\sin(ij\pi h) - 2\sin(ij\pi h)\cos(j\pi h) = 2(1 - \cos j\pi h)\sin(ij\pi h)$. (At $i = 1$ and $i = N$ the missing neighbour is $\sin 0 = 0$ and $\sin((N+1)j\pi h) = \sin(j\pi) = 0$, so the formula holds there too; this is why the boundary condition $u = 0$ and the spacing $h = 1/(N+1)$ fit together.)
2. So $L_N v^{(j)} = \mu_j v^{(j)}$ with $\mu_j = 2 - 2\cos(j\pi h) = 4\sin^2(j\pi h/2)$, $j = 1, \dots, N$: $N$ distinct eigenvalues, hence all of them.
3. For the Kronecker sum, $(I\otimes L + L\otimes I)(v^{(j)}\otimes v^{(k)}) = v^{(j)}\otimes L v^{(k)} + L v^{(j)}\otimes v^{(k)} = (\mu_k + \mu_j)(v^{(j)}\otimes v^{(k)})$, using $(B\otimes C)(x\otimes y) = Bx\otimes Cy$ (LA-12). Divide by $h^2$ to get (3.22).
4. $\lambda_{\min} = \lambda_{11} = (8/h^2)\sin^2(\pi h/2) \to 2\pi^2$ as $h\to 0$ (since $\sin x \approx x$); $\lambda_{\max} = \lambda_{NN} = (8/h^2)\cos^2(\pi h/2)$ (since $\sin(N\pi h/2) = \cos(\pi h/2)$); their ratio is $\cot^2(\pi h/2) \approx 4/(\pi^2h^2)$. ∎

■ **What just happened.** The condition number grows like $N^2$: refining the grid to get a better answer makes the linear system harder, and the digits-lost rule of §3.3 says a $1000\times 1000$ grid ($\kappa \approx 4\times 10^5$) still leaves ten digits, which is plenty. The project checks (3.22) against `eigvalsh` for $N = 8$. The same calculation gives the vibration frequencies of a square drum (CM-6) and the eigenstates of a particle in a 2D box (QM-2); the matrix $L_N$ is the discrete Laplacian whose eigenvectors are the discrete sine transform (F2). The $\sin^2$ structure, small eigenvalues for smooth modes and large ones for oscillatory modes, is what makes multigrid methods work, and what makes Jacobi slow.

**Sparse storage.** Only the nonzeros are stored. Compressed sparse row (CSR) format keeps three arrays: the nonzero values row by row, their column indices, and a pointer to where each row starts; $A x$ costs $2\,\mathrm{nnz}(A)$ flops, here $10n$. For $N = 200$, one matrix–vector product is $4\times 10^5$ flops against $\tfrac23 n^3 = 4\times 10^{13}$ for dense elimination: a factor of $10^8$. So every method below is built from matrix–vector products only.

**Stationary iterations from a splitting.** Write $A = M - N$ with $M$ easy to invert (diagonal, triangular). Then $Ax = b$ becomes $Mx = Nx + b$, and the fixed-point iteration is

$$x_{k+1} = M^{-1}(N x_k + b) = x_k + M^{-1} r_k, \qquad r_k = b - Ax_k. \tag{3.23}$$

(The second form: $M^{-1}(Nx_k + b) = M^{-1}((M - A)x_k + b) = x_k + M^{-1}(b - Ax_k)$. It says: correct $x_k$ by an approximate solve of the residual equation $A\,e = r_k$ with $M$ in place of $A$.) **Jacobi** takes $M = D$, the diagonal of $A$; **Gauss–Seidel** takes $M = D + L_{\text{low}}$, the lower triangle including the diagonal, so that each updated component is used immediately for the next one. Both cost about one matrix–vector product per step.

**Theorem 3.19 (Convergence of stationary iterations).** Let $B = M^{-1}N = I - M^{-1}A$ be the iteration matrix. The iteration (3.23) converges to $x = A^{-1}b$ for every starting vector $x_0$ if and only if the spectral radius $\rho(B) = \max_i|\lambda_i(B)| < 1$.

▶ **Try it first.** Subtract $x = M^{-1}(Nx + b)$ (true for the exact solution) from (3.23) to get a recurrence for the error $e_k = x_k - x$. Then ask when $B^k \to 0$.

**Derivation.**

1. $e_{k+1} = x_{k+1} - x = M^{-1}N(x_k - x) = B e_k$, so $e_k = B^k e_0$. Convergence for every $e_0$ is the statement $B^k \to 0$.
2. ($\Leftarrow$) Suppose $\rho(B) < 1$. Take the Jordan form $B = XJX^{-1}$ (LA-11 Theorem 3.10); then $B^k = XJ^kX^{-1}$ and it suffices that $J^k \to 0$ block by block. For a Jordan block $J_s(\lambda) = \lambda I + N_s$ with $N_s^s = 0$, LA-11 Theorem 3.17 gives $J_s(\lambda)^k = \sum_{i=0}^{s-1}\binom{k}{i}\lambda^{k-i}N_s^i$. Each term is bounded by $k^i|\lambda|^{k-i}$, and $k^i|\lambda|^k \to 0$ for $|\lambda| < 1$ because the exponential beats the polynomial. So every block, hence $J^k$, hence $B^k$, tends to zero.
3. ($\Rightarrow$) Suppose $|\lambda| \ge 1$ for some eigenvalue with eigenvector $v$ (possibly complex). If $B^k e_0 \to 0$ for every real $e_0$ then, by linearity, also for every complex $e_0 = u + iw$. But $B^k v = \lambda^k v$ has $\|B^kv\| = |\lambda|^k\|v\| \ge \|v\|$, which does not tend to zero. Contradiction. ∎

■ **What just happened.** The asymptotic error reduction per step is $\rho(B)$: the number of iterations to gain one decimal digit is $1/(-\log_{10}\rho)$. Note that $\|B\| < 1$ in some norm is sufficient but $\rho(B) < 1$ is the exact condition; a matrix can have $\|B\|_2 > 1$ and still $B^k \to 0$ (any Jordan block with $|\lambda| < 1$ and a large off-diagonal), in which case the error first grows and then decays.

**Jacobi on the model problem.** Here $D = (4/h^2)I$, so $B_J = I - D^{-1}A = I - (h^2/4)A$ has the *same eigenvectors* as $A$ with eigenvalues $1 - (h^2/4)\lambda_{jk} = 1 - \sin^2\frac{j\pi h}{2} - \sin^2\frac{k\pi h}{2} = \frac{\cos j\pi h + \cos k\pi h}{2}$. The largest in modulus is at $j = k = 1$ (and $j = k = N$, with the opposite sign):

$$\rho(B_J) = \cos\pi h \approx 1 - \frac{\pi^2h^2}{2}. \tag{3.24}$$

For Gauss–Seidel on this matrix, $\rho(B_{GS}) = \rho(B_J)^2 = \cos^2\pi h$ (David Young's 1950 thesis, for matrices with "property A," which the five-point stencil has; proof in Demmel, *Applied Numerical Linear Algebra* §6.5.3). Gauss–Seidel is exactly twice as fast as Jacobi, and both need $O(1/h^2) = O(N^2)$ iterations to gain a fixed number of digits, since $-\ln\cos\pi h \approx \pi^2h^2/2$. Each iteration is $O(N^2)$ work, so the total is $O(N^4)$: no better than banded elimination. For $N = 20$ ($\kappa = 178$, $\rho_J = 0.98883$) the project needs $1416$ Jacobi and $703$ Gauss–Seidel iterations to reduce the residual by $10^{-8}$, with the observed late-stage rates $0.98883$ and $0.97779$ matching $\cos\pi h$ and $\cos^2\pi h$ to five digits. Conjugate gradients, next, does it in $63$.

The slowness has a picture. The eigenvalue of $B_J$ for the mode $(j, k)$ is $(\cos j\pi h + \cos k\pi h)/2$: near $1$ for the smoothest modes ($j, k$ small) and near $-1$ for the most oscillatory ones. Jacobi barely touches either. Weighted Jacobi ($M = D/\omega$, $\omega = 2/3$) kills the oscillatory half of the spectrum in a few sweeps and leaves the smooth half, which looks oscillatory on a coarser grid. That observation is multigrid, the only method that solves the model problem in $O(n)$ work; it is beyond this module (Briggs, Henson & McCormick, *A Multigrid Tutorial*).

### 3.10 Conjugate gradients, derived

Budget three hours for this section and do not skip the derivation steps; ML-2 reuses this argument for gradient descent on quadratics, and every Krylov method in §3.11 is a variation on it. Throughout, $A$ is $n\times n$ SPD, and $\|v\|_A = \sqrt{v^\top A v}$ is the **$A$-norm** (a genuine norm because $A$ is positive definite, LA-8 §3.8), with the associated **$A$-inner product** $\langle u, v\rangle_A = u^\top A v$.

**(i) Solving $Ax = b$ is minimizing a quadratic.**

**Theorem 3.20.** For SPD $A$, the function $\varphi(x) = \tfrac12 x^\top A x - b^\top x$ has a unique minimizer, namely $x_* = A^{-1}b$, and for every $x$,
$$\varphi(x) - \varphi(x_*) = \tfrac12\|x - x_*\|_A^2. \tag{3.25}$$

▶ **Try it first.** Compute the gradient of $\varphi$ (LA-14 if you want the formal rule; or expand in components). Then expand $\tfrac12(x - x_*)^\top A(x - x_*)$.

**Derivation.**

1. $\nabla\varphi(x) = Ax - b = -r(x)$: the negative residual. ($\partial/\partial x_i$ of $\tfrac12\sum_{jk}a_{jk}x_jx_k$ is $\tfrac12\sum_k a_{ik}x_k + \tfrac12\sum_j a_{ji}x_j = (Ax)_i$ by symmetry.) Critical points satisfy $Ax = b$; there is exactly one, $x_*$.
2. $\tfrac12(x - x_*)^\top A(x - x_*) = \tfrac12x^\top Ax - x^\top Ax_* + \tfrac12x_*^\top Ax_* = \tfrac12x^\top Ax - x^\top b + \tfrac12x_*^\top b$, using $Ax_* = b$. And $\varphi(x_*) = \tfrac12x_*^\top b - x_*^\top b = -\tfrac12x_*^\top b$. So the right side equals $\varphi(x) - \varphi(x_*)$, which is (3.25).
3. Since $\|x - x_*\|_A^2 > 0$ for $x \ne x_*$, the critical point is the strict global minimizer. ∎

■ **What just happened.** The level sets of $\varphi$ are ellipsoids centered at $x_*$ with axes along the eigenvectors of $A$ and semi-axis lengths proportional to $1/\sqrt{\lambda_i}$ (LA-8 Theorem 3.14): long in the direction of small eigenvalues, short along large ones, with aspect ratio $\sqrt\kappa$. Solving the system means finding the bottom of a long thin valley.

**(ii) Steepest descent zigzags.** The obvious method moves downhill: $x_{k+1} = x_k + \alpha_k r_k$ along the negative gradient $r_k = b - Ax_k$, with $\alpha_k$ chosen to minimize $\varphi$ along the line.

▶ **Try it first.** Minimize $\varphi(x_k + \alpha r_k)$ over $\alpha$; then compute $r_{k+1}^\top r_k$.

**Derivation.**

1. $\frac{d}{d\alpha}\varphi(x_k + \alpha r_k) = \nabla\varphi(x_k + \alpha r_k)^\top r_k = (A(x_k + \alpha r_k) - b)^\top r_k = -r_k^\top r_k + \alpha\, r_k^\top A r_k$. Setting it to zero: $\alpha_k = \dfrac{r_k^\top r_k}{r_k^\top A r_k}$.
2. The new residual is $r_{k+1} = b - A(x_k + \alpha_k r_k) = r_k - \alpha_k A r_k$, and $r_{k+1}^\top r_k = r_k^\top r_k - \alpha_k r_k^\top A r_k = 0$ by the choice of $\alpha_k$. ∎

Consecutive search directions are *orthogonal*: the path turns $90°$ at every step, zigzagging down the valley. In the eigenbasis the error component along the $i$th axis is multiplied by $(1 - \alpha_k\lambda_i)$ at each step, and no single $\alpha$ can be good for both $\lambda_{\min}$ and $\lambda_{\max}$. The resulting bound (Kantorovich; Luenberger, *Linear and Nonlinear Programming* §8.6) is $\|e_{k}\|_A \le \left(\frac{\kappa - 1}{\kappa + 1}\right)^k\|e_0\|_A$: about $\kappa/2$ steps per factor of $e$, i.e. $O(\kappa)$ iterations, no better than Jacobi on the model problem. Exercise B7 runs two steps by hand on the $2\times 2$ system of Example 4.5: steepest descent reaches $(1/12, 7/12)$ while the answer is $(1/11, 7/11)$, and it never terminates; CG, which starts with the *same* first step, is exact at step two.

**(iii) Conjugate directions finish in $n$ steps.** The repair is to choose directions that do not undo each other's work: directions $p_0, \dots, p_{n-1}$ that are **$A$-conjugate**, $p_i^\top A p_j = 0$ for $i \ne j$ (orthogonal in the $A$-inner product).

**Theorem 3.21 (Conjugate directions).** Let $p_0, \dots, p_{n-1}$ be nonzero and $A$-conjugate. Then they are linearly independent, and the iteration $x_{k+1} = x_k + \alpha_k p_k$ with the exact line-search step
$$\alpha_k = \frac{p_k^\top r_k}{p_k^\top A p_k} \tag{3.26}$$
reaches $x_n = x_*$ exactly, for any $x_0$.

▶ **Try it first.** Expand the initial error $e_0 = x_* - x_0$ in the basis $\{p_i\}$ and use conjugacy to read off the coefficients. Then show those coefficients are the $\alpha_k$.

**Derivation.**

1. Independence: if $\sum_i c_ip_i = 0$, multiply by $p_j^\top A$: $c_j\,p_j^\top Ap_j = 0$, and $p_j^\top Ap_j > 0$, so $c_j = 0$.
2. Since $\{p_i\}$ is a basis, $e_0 = x_* - x_0 = \sum_{i=0}^{n-1}c_ip_i$. Multiply by $p_k^\top A$: $p_k^\top Ae_0 = c_k\,p_k^\top Ap_k$, so $c_k = \dfrac{p_k^\top Ae_0}{p_k^\top Ap_k} = \dfrac{p_k^\top r_0}{p_k^\top Ap_k}$, using $Ae_0 = b - Ax_0 = r_0$.
3. The iteration gives $x_k = x_0 + \sum_{i<k}\alpha_ip_i$, so $r_k = r_0 - \sum_{i<k}\alpha_iAp_i$ and $p_k^\top r_k = p_k^\top r_0 - \sum_{i<k}\alpha_i\,p_k^\top Ap_i = p_k^\top r_0$ by conjugacy. Therefore $\alpha_k$ in (3.26) equals $c_k$.
4. Hence $x_n = x_0 + \sum_{k<n}c_kp_k = x_0 + e_0 = x_*$.
5. (That (3.26) is the exact line search: as in (ii), $\frac{d}{d\alpha}\varphi(x_k + \alpha p_k) = -p_k^\top r_k + \alpha\,p_k^\top Ap_k$.) ∎

■ **What just happened.** In the $A$-inner product the directions are orthogonal, and the method is nothing but "expand the error in an orthogonal basis, one coefficient per step." Each step kills one component of the error for good, because later steps are $A$-orthogonal to it. The question is where to get $A$-conjugate directions cheaply. Any basis can be $A$-orthogonalized by Gram–Schmidt (LA-5 §3.8 with $\langle\cdot,\cdot\rangle_A$), but that costs $O(n^2)$ per direction plus all previous directions in storage. The miracle of CG is that for one particular choice, Gram–Schmidt collapses to a two-term recurrence.

**(iv) Gram–Schmidt on the residuals collapses.** Take $p_0 = r_0$ and, at step $k$, build $p_k$ by $A$-orthogonalizing the current residual $r_k$ against $p_0, \dots, p_{k-1}$:
$$p_k = r_k - \sum_{j<k}\frac{p_j^\top Ar_k}{p_j^\top Ap_j}\,p_j. \tag{3.27}$$

**Theorem 3.22 (The CG orthogonality relations).** With $p_0 = r_0$, (3.26), (3.27) and $x_{k+1} = x_k + \alpha_kp_k$, for as long as $r_k \ne 0$:
(a) $\operatorname{span}\{p_0,\dots,p_k\} = \operatorname{span}\{r_0,\dots,r_k\} = \operatorname{span}\{r_0, Ar_0, \dots, A^kr_0\} =: \mathcal{K}_{k+1}$ (the **Krylov subspace**);
(b) $r_k \perp \mathcal{K}_k$, i.e. $r_k^\top p_j = 0$ and $r_k^\top r_j = 0$ for all $j < k$;
(c) in (3.27) every term with $j < k-1$ vanishes, so
$$p_k = r_k + \beta_{k-1}p_{k-1}, \qquad \beta_{k-1} = \frac{r_k^\top r_k}{r_{k-1}^\top r_{k-1}}; \tag{3.28}$$
(d) $\alpha_k = \dfrac{r_k^\top r_k}{p_k^\top Ap_k}$.

▶ **Try it first.** Prove (b) from the fact that $\alpha_k$ is an exact line search (so $r_{k+1} \perp p_k$) plus conjugacy. Everything else follows from (b).

**Derivation (induction on $k$; all three statements in (a) hold trivially for $k = 0$).**

1. *(a), the spans.* By (3.27), $p_k \in \operatorname{span}\{r_k, p_0,\dots,p_{k-1}\}$ and conversely $r_k \in \operatorname{span}\{p_0,\dots,p_k\}$, so the first two spans agree. For the third: $r_{k} = r_{k-1} - \alpha_{k-1}Ap_{k-1}$, and $Ap_{k-1} \in A\mathcal{K}_{k} \subseteq \mathcal{K}_{k+1}$, so by induction $r_k \in \mathcal{K}_{k+1}$ and $\operatorname{span}\{r_0..r_k\} \subseteq \mathcal{K}_{k+1}$. The dimensions match (the $p_j$ are independent by Theorem 3.21, as long as $r_k \ne 0$ makes $p_k \ne 0$: if $p_k = 0$ then $r_k \in \operatorname{span}\{p_0..p_{k-1}\}$, contradicting (b) below unless $r_k = 0$), so the spans are equal.
2. *(b), residuals are orthogonal to the past.* We show $r_{k+1} \perp p_j$ for all $j \le k$. For $j = k$: $r_{k+1}^\top p_k = r_k^\top p_k - \alpha_k\,p_k^\top Ap_k = 0$ by (3.26). For $j < k$: $r_{k+1}^\top p_j = r_k^\top p_j - \alpha_k\,p_k^\top Ap_j = 0 - 0$, by the induction hypothesis and conjugacy. Since the $p_j$ span $\mathcal{K}_{k+1}$ by (a), $r_{k+1} \perp \mathcal{K}_{k+1}$; in particular $r_{k+1} \perp r_j$ for $j \le k$. **The residuals are mutually orthogonal.**
3. *(c), the collapse.* The coefficient of $p_j$ in (3.27) has numerator $p_j^\top Ar_k = r_k^\top(Ap_j)$. From $r_{j+1} = r_j - \alpha_jAp_j$, $Ap_j = (r_j - r_{j+1})/\alpha_j \in \operatorname{span}\{r_j, r_{j+1}\}$. For $j < k-1$ both $r_j$ and $r_{j+1}$ are orthogonal to $r_k$ by (b), so the numerator is zero. Only $j = k-1$ survives, with numerator $r_k^\top(r_{k-1} - r_k)/\alpha_{k-1} = -r_k^\top r_k/\alpha_{k-1}$ and denominator $p_{k-1}^\top Ap_{k-1} = r_{k-1}^\top r_{k-1}/\alpha_{k-1}$ (by (d) at step $k-1$). The ratio, with the minus sign from (3.27), is $\beta_{k-1} = r_k^\top r_k/r_{k-1}^\top r_{k-1}$.
4. *(d).* From (3.28), $p_k^\top r_k = r_k^\top r_k + \beta_{k-1}p_{k-1}^\top r_k = r_k^\top r_k$ by (b). Substitute into (3.26). ∎

■ **What just happened.** This is the center of the module. Gram–Schmidt against all previous directions should cost $O(nk)$ at step $k$ and require storing everything. Because the search directions are built from residuals, and residuals are orthogonal (step 2), and $Ap_j$ lives in the span of two consecutive residuals (step 3), all but one term vanish *exactly*. The algorithm needs one matrix–vector product, two inner products, three vector updates, and four vectors of storage per step, regardless of $k$. Hestenes and Stiefel found this in 1952; the fact that it is "Gram–Schmidt in disguise" is why it was misfiled as a direct method for twenty years.

**(v) The algorithm.**

**Algorithm 3.3 (Conjugate gradients; Hestenes–Stiefel 1952).** Given SPD $A$, $b$, and $x_0$ (often $0$): $r_0 = b - Ax_0$, $p_0 = r_0$. For $k = 0, 1, 2, \dots$ until $\|r_k\| \le \tau\|b\|$:

```
α_k = (r_kᵀ r_k) / (p_kᵀ A p_k)
x_{k+1} = x_k + α_k p_k
r_{k+1} = r_k − α_k A p_k
β_k = (r_{k+1}ᵀ r_{k+1}) / (r_kᵀ r_k)
p_{k+1} = r_{k+1} + β_k p_k
```

Cost per step: one product $Ap_k$ ($2\,\mathrm{nnz}$ flops), two inner products and three axpys ($10n$ flops). Example 4.5 runs it by hand on a $2\times 2$ and terminates at step $2$ with $r_2 = 0$, as Theorem 3.21 requires. In floating point the orthogonality of step 2 degrades and exact termination at step $n$ is lost; that does not matter, because the method is used with $k \ll n$, and the next two theorems say why that works.

**(vi) CG is optimal in the $A$-norm.**

**Theorem 3.23 (Optimality).** The CG iterate $x_k$ minimizes $\|x - x_*\|_A$ over all $x \in x_0 + \mathcal{K}_k$. Equivalently, by (3.25), it minimizes $\varphi$ over that affine space.

▶ **Try it first.** $e_k = x_* - x_k$ satisfies $Ae_k = r_k \perp \mathcal{K}_k$. Translate that into a statement about the $A$-inner product, then use Pythagoras.

**Derivation.**

1. For any $y \in x_0 + \mathcal{K}_k$, write $y = x_k + d$ with $d \in \mathcal{K}_k$ (since $x_k \in x_0 + \mathcal{K}_k$ by Theorem 3.22(a)).
2. $\langle e_k, d\rangle_A = e_k^\top Ad = (Ae_k)^\top d = r_k^\top d = 0$ by Theorem 3.22(b). So $e_k \perp_A \mathcal{K}_k$.
3. $\|x_* - y\|_A^2 = \|e_k - d\|_A^2 = \|e_k\|_A^2 - 2\langle e_k, d\rangle_A + \|d\|_A^2 = \|e_k\|_A^2 + \|d\|_A^2 \ge \|e_k\|_A^2$, with equality only for $d = 0$. ∎

The polynomial reading, which is how convergence is analyzed: every $x \in x_0 + \mathcal{K}_k$ has $x - x_0 = q(A)r_0$ for some polynomial $q$ of degree $< k$, hence $e = e_0 - q(A)Ae_0 = p(A)e_0$ with $p(t) = 1 - tq(t)$, a polynomial of degree $\le k$ with $p(0) = 1$. Expanding $e_0 = \sum c_ix_i$ in the orthonormal eigenbasis, $\|p(A)e_0\|_A^2 = \sum_i\lambda_ic_i^2p(\lambda_i)^2 \le \max_i p(\lambda_i)^2\,\|e_0\|_A^2$. Theorem 3.23 then gives

$$\frac{\|e_k\|_A}{\|e_0\|_A} \le \min_{\deg p \le k,\ p(0) = 1}\ \max_{\lambda\in\sigma(A)}|p(\lambda)|. \tag{3.29}$$

Two consequences. If $A$ has only $m$ distinct eigenvalues, the polynomial $\prod_{i=1}^m(1 - t/\lambda_i)$ vanishes on all of them, so CG terminates in $m$ steps. And if the eigenvalues are clustered, a low-degree polynomial can be small on all the clusters: CG responds to the *distribution* of the spectrum, not only to $\kappa$. (The project's jump-coefficient example has $\kappa$ barely changed by preconditioning yet converges ten times faster, because the preconditioned spectrum is clustered.)

**(vii) The $\sqrt\kappa$ bound.**

**Theorem 3.24 (CG convergence; Chebyshev bound).** With $\kappa = \lambda_{\max}/\lambda_{\min}$,
$$\frac{\|e_k\|_A}{\|e_0\|_A} \le 2\left(\frac{\sqrt\kappa - 1}{\sqrt\kappa + 1}\right)^k. \tag{3.30}$$

*Sketch (Exercise C2 fills in the Chebyshev facts).* In (3.29) we may replace $\sigma(A)$ by the whole interval $[\lambda_{\min}, \lambda_{\max}]$, which can only increase the maximum. The polynomial of degree $k$ with $p(0) = 1$ that is smallest on an interval not containing $0$ is a scaled and shifted Chebyshev polynomial: with $T_k(\cos\theta) = \cos k\theta$ (so $|T_k| \le 1$ on $[-1, 1]$) and $\ell(\lambda) = \dfrac{\lambda_{\max} + \lambda_{\min} - 2\lambda}{\lambda_{\max} - \lambda_{\min}}$, which maps $[\lambda_{\min}, \lambda_{\max}]$ onto $[-1, 1]$ and $0$ to $\gamma = \dfrac{\kappa + 1}{\kappa - 1} > 1$, take $p(\lambda) = T_k(\ell(\lambda))/T_k(\gamma)$. Then $\max_{[\lambda_{\min},\lambda_{\max}]}|p| = 1/|T_k(\gamma)|$. Outside $[-1,1]$, $T_k(\gamma) = \tfrac12\big[(\gamma + \sqrt{\gamma^2 - 1})^k + (\gamma - \sqrt{\gamma^2-1})^k\big] \ge \tfrac12(\gamma + \sqrt{\gamma^2-1})^k$, and the algebra $\gamma^2 - 1 = \dfrac{4\kappa}{(\kappa-1)^2}$ gives $\gamma + \sqrt{\gamma^2 - 1} = \dfrac{\kappa + 1 + 2\sqrt\kappa}{\kappa - 1} = \dfrac{(\sqrt\kappa+1)^2}{(\sqrt\kappa - 1)(\sqrt\kappa+1)} = \dfrac{\sqrt\kappa + 1}{\sqrt\kappa - 1}$. Hence $1/T_k(\gamma) \le 2\left(\frac{\sqrt\kappa-1}{\sqrt\kappa+1}\right)^k$. ∎

■ **What just happened.** Since $\frac{\sqrt\kappa - 1}{\sqrt\kappa+1} \approx 1 - 2/\sqrt\kappa$, the error drops by a factor $e$ roughly every $\sqrt\kappa/2$ steps: **$O(\sqrt\kappa)$ iterations** against $O(\kappa)$ for steepest descent and Jacobi. For the model problem $\sqrt\kappa \approx 0.64N$, so CG needs $O(N)$ iterations of $O(N^2)$ work each: $O(N^3) = O(n^{3/2})$ total, against $O(N^4)$ for Gauss–Seidel or banded elimination and $O(N^6)$ for dense elimination. The project confirms the scaling: from $N = 20$ to $N = 200$ the iteration count to reach $\|r\|/\|b\| < 10^{-8}$ is $4.7\sqrt\kappa \pm 1.5\%$, i.e. $600$ iterations for $40{,}000$ unknowns, and the bound (3.30) on the $A$-norm error holds at every step (it is loose by a factor of about $400$ at $k = 50$ on the $N = 40$ grid, because the actual spectrum is not the whole interval). Note what is bounded: the $A$-norm of the *error*, which you cannot compute without $x_*$; the residual you *can* compute is not monotone in CG and can jump around, which confuses people the first time (§8).

**Why it was forgotten for twenty years.** In exact arithmetic CG terminates in $n$ steps, so it was judged as a direct method and compared with elimination: for dense $A$ its $n$ steps cost $2n^3$ flops, three times LU, and in floating point the "exact at step $n$" property fails because the residuals lose orthogonality. Reid (1971) pointed out that for sparse, well-conditioned or preconditioned $A$ one should stop after $k \ll n$ steps and judge it by (3.30). Everything in this section was in Hestenes and Stiefel's paper; the use was not.
