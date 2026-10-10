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

**Proposition 3.4 (Matrix–vector product).** For $f(x) = Ax$ with fixed $A$, $\kappa \le \|A\|\,\|x\|/\|Ax\|$, with equality for $A$ square and invertible if and only if... no: with the bound $\kappa \le \|A\|\|A^{-1}\|$ when $A$ is invertible.

*Proof.* The Jacobian of $x \mapsto Ax$ is $A$, so by the remark after (3.6), $\kappa = \|A\|\|x\|/\|Ax\|$ exactly; the inequality in the statement is this with "$\le$" because some authors define $\kappa$ with the supremum taken differently. If $A$ is invertible then $\|x\| = \|A^{-1}Ax\| \le \|A^{-1}\|\|Ax\|$, so $\kappa \le \|A\|\|A^{-1}\|$. ∎

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
