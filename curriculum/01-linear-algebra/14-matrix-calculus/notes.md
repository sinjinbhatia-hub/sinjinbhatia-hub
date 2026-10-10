# LA-14: Matrix Calculus: Gradients, Jacobians, Hessians, the Chain Rule in Matrix Form, and Backpropagation

**Time:** ~10 hours · **Prerequisites:** LA-9 (and through it LA-1, LA-5, LA-6, LA-8; LA-12 is referenced for the Riesz map and `vec`, not required) · **Unlocks:** ML-2, ML-7, ST-1

> **Conventions in this module.** Vectors are columns, $\mathbf{x} \in \mathbb{R}^n$. For a scalar function $f : \mathbb{R}^n \to \mathbb{R}$ the **gradient** $\nabla f(\mathbf{x})$ is a **column** vector in $\mathbb{R}^n$. For $\mathbf{f} : \mathbb{R}^n \to \mathbb{R}^m$ the **Jacobian** $J$ is the $m \times n$ matrix with $J_{ij} = \partial f_i / \partial x_j$, so that $\mathbf{f}(\mathbf{x} + \mathbf{h}) \approx \mathbf{f}(\mathbf{x}) + J\mathbf{h}$; for scalar $f$ the Jacobian is the $1 \times n$ **row** $\nabla f^\top$, i.e. $\nabla f = J^\top$. The derivative of a scalar with respect to a matrix $X$ is a matrix $\nabla_X f$ of the **same shape as $X$**. We compute every derivative by **differentials**: write $df$ as a linear expression in $dX$, put it in the form $df = \operatorname{tr}(G^\top dX)$, and read off $\nabla_X f = G$; for vectors, $df = \mathbf{g}^\top d\mathbf{x}$ gives $\nabla f = \mathbf{g}$. We never use "numerator layout" or "denominator layout"; the differential method makes that choice unnecessary, and §3.3 says why. Other notation: $\odot$ is the elementwise (Hadamard) product; $\mathbf{1}$ is the all-ones column; $\operatorname{tr}$ is the trace; $\|A\|_F^2 = \operatorname{tr}(A^\top A)$ is the Frobenius norm; $\sigma$ is a scalar function applied elementwise and $\sigma'$ its derivative applied elementwise; $\operatorname{diag}(\mathbf{v})$ is the diagonal matrix with $\mathbf{v}$ on the diagonal; $\|\cdot\|$ is the Euclidean norm; $o(\|\mathbf{h}\|)$ means a quantity $r(\mathbf{h})$ with $\|r(\mathbf{h})\|/\|\mathbf{h}\| \to 0$ as $\mathbf{h} \to \mathbf{0}$. All functions are assumed continuously differentiable (twice, where Hessians appear).

## 0. Where this sits

You know multivariable calculus (partial derivatives, the one-variable chain rule, Taylor's theorem) and, from LA-1 to LA-9, how to think of a matrix as a linear map, of the dot product as an inner product, of symmetric matrices through the spectral theorem, and of least squares through projections. This module puts the two together: the derivative of a map between vector spaces *is* a linear map, so everything you know about matrices applies to derivatives. ML-2 (optimization: gradient descent, Newton, convergence proofs) needs every identity here. ML-7 (neural networks) is backpropagation, which §3.7 derives in full. ST-1 (maximum likelihood, Fisher information) needs the gradient of a log-likelihood with respect to a mean vector and a covariance matrix, which is §3.9.

## 1. The problem that forced this idea

Calculus of one variable gives you a number, the slope. Already with two variables you get two slopes, and the moment you ask how a *vector* of outputs responds to a *vector* of inputs you need a rectangular array of them. Nobody minded this while the arrays were small. Three episodes made it unavoidable.

**Fitting curves by iteration.** Gauss, fitting orbital elements to astronomical observations (*Theoria Motus*, 1809), minimized a sum of squared residuals in which the parameters entered nonlinearly. His method was to linearize the residuals about a current guess, solve the resulting *linear* least-squares problem (LA-5), move to the solution, and repeat. The linearization is a Jacobian; the linear solve is the normal equations. The method is now called Gauss–Newton and is still the workhorse of nonlinear least squares (§3.5, Example 4.5).

**Likelihoods with matrix parameters.** Fisher's theory of maximum likelihood (1922, 1925) asks you to maximize $\ell(\theta) = \log p(\text{data} \mid \theta)$ by setting its derivative to zero. For the multivariate normal, $\theta$ contains a covariance *matrix* $\Sigma$, and $\ell$ contains $\log\det\Sigma$ and $(\mathbf{x}-\boldsymbol{\mu})^\top\Sigma^{-1}(\mathbf{x}-\boldsymbol{\mu})$. "Differentiate with respect to $\Sigma$" is not something one-variable calculus prepares you for, and doing it entry by entry is a swamp. Statisticians developed a cleaner tool, the **differential method**, systematized by Magnus and Neudecker (1988) and taught compactly by Minka (2000): differentiate without ever choosing how to lay out a derivative array, then read off the answer at the end. We adopt it wholesale (§3.3).

**Compositions of hundreds of maps.** A neural network is a composition $f_k \circ \cdots \circ f_1$ of simple vector maps, and training it means computing the gradient of a scalar loss with respect to every parameter in every layer. The chain rule says the derivative of a composition is the product of the derivatives, so in principle this is "multiply $k$ Jacobians." Done naively that costs a factor of the layer width more than evaluating the network. Werbos (1974 thesis) derived the efficient way to do it through a general ordered system; Linnainmaa (1970 thesis) had found the same reverse accumulation for tracking rounding errors; Rumelhart, Hinton and Williams (1986) stated it for neural networks, named it **backpropagation**, and showed it trains hidden layers. Griewank and others in automatic differentiation proved the general fact that makes all this work: the gradient of a scalar function can always be computed for at most a small constant multiple of the cost of the function itself (the "cheap gradient principle," building on Baur and Strassen, 1983). Griewank's 2012 essay "Who invented the reverse mode of differentiation?" is the honest history.

Three episodes, one module. The thesis in one line: **the derivative is a linear map; the chain rule is matrix multiplication; backpropagation is the chain rule evaluated in the cheap order.**

## 2. Intuition

### The derivative is the best linear approximation, so it is a matrix

In one variable, $f(x + h) \approx f(x) + f'(x)\,h$: near $x$ the function is a line with slope $f'(x)$. The slope is the thing that multiplies $h$. In many variables, $\mathbf{f}(\mathbf{x} + \mathbf{h}) \approx \mathbf{f}(\mathbf{x}) + (\text{something})\,\mathbf{h}$, and "something" must turn a vector $\mathbf{h} \in \mathbb{R}^n$ into a vector in $\mathbb{R}^m$ linearly. By LA-4, that something is an $m \times n$ matrix. That matrix is the derivative. There is nothing more to the definition than this.

Picture $\mathbf{f} : \mathbb{R}^2 \to \mathbb{R}^2$ as a deformation of the plane: draw a square grid, apply $\mathbf{f}$ to every grid point, look at the bent grid. Zoom in on one point $\mathbf{p}$. A tiny square around $\mathbf{p}$ maps to a tiny *parallelogram*, because up close every smooth map looks linear, and linear maps send squares to parallelograms (LA-1, LA-6). The matrix that sends the square to the parallelogram is the Jacobian $J(\mathbf{p})$. The project draws this for $\mathbf{f}(x, y) = (x^2 - y^2,\; 2xy)$ (figure `jacobian_grid.png`): the image of a small square at $(1, 0.5)$ and the parallelogram $\mathbf{f}(\mathbf{p}) + J(\mathbf{p})(\mathbf{q} - \mathbf{p})$ are visually indistinguishable.

```
   domain                            image under f
   +----+----+                       the square near p becomes
   |    |    |      f                a parallelogram: f(p) + J (q - p)
   +----#----+   ------->                 ___
   |    |    |                           /  /   <- edges are J e_1, J e_2
   +----+----+                          /__/       (the columns of J)
       p
```

The columns of $J$ are the images of the edges of the square: column $j$ is $\partial\mathbf{f}/\partial x_j$, "how fast the output moves when input $j$ moves." That is Theorem 3.5 of LA-1 ($A\mathbf{e}_j$ is column $j$) applied to the derivative.

### The gradient points uphill, and it is a column for a reason

For a scalar $f : \mathbb{R}^n \to \mathbb{R}$, the Jacobian is a single row $(\partial f/\partial x_1, \dots, \partial f/\partial x_n)$. Multiply it by a direction $\mathbf{u}$ and you get the rate of change of $f$ in that direction. By Cauchy–Schwarz (LA-1, Theorem 3.11), the rate is largest when $\mathbf{u}$ is parallel to the vector of partials. So the vector of partials, stood up as a column, is the direction of steepest ascent. That column is the gradient. Row and column hold the same numbers; the row *acts on* directions, the column *is* a direction. Keeping them straight is most of what "layout conventions" are about, and in §3.1 the distinction gets a name: the row is a linear functional, the column is its Riesz representative (LA-12).

### The Hessian is the local shape

Go one order further. Near a critical point, $f(\mathbf{x} + \mathbf{h}) \approx f(\mathbf{x}) + \tfrac12 \mathbf{h}^\top H \mathbf{h}$, a quadratic form (LA-8). $H$ is symmetric, so by the spectral theorem it has orthogonal eigenvectors: those are the principal directions of curvature, and the eigenvalues are the curvatures. All positive: a bowl, a minimum. All negative: a dome. Mixed: a saddle. You already know this picture from LA-8's quadratic forms; the only new fact is that the matrix is made of second partial derivatives.

### The chain rule is matrix multiplication, and there is a cheap order

If $\mathbf{f}$ is locally the matrix $J_f$ and $\mathbf{g}$ is locally $J_g$, then $\mathbf{g} \circ \mathbf{f}$ is locally $J_g J_f$: do one linear map, then the other, and LA-1 §3.3 says the composition is the product. So a network of $k$ layers, each of width $d$, has derivative $J_k J_{k-1} \cdots J_1$, a product of $k$ matrices of size $d \times d$ (the last one $1 \times d$ because the loss is a scalar).

Here is the thing to guess before §3.6 proves it. Multiplying these from the right, $J_2 J_1$ first, costs about $d^3$ per product. Multiplying from the *left*, starting with the $1 \times d$ row, every product is a row times a matrix: $d^2$. Same answer, a factor of $d$ cheaper. For $d = 1000$ that is the difference between a training run and no training run. Backpropagation is the left-to-right order. Nothing else.

**Guess the main results now.** (i) The derivative at a point is a matrix, the Jacobian, and $\mathbf{f}(\mathbf{x} + \mathbf{h}) = \mathbf{f}(\mathbf{x}) + J\mathbf{h} + (\text{smaller than } \mathbf{h})$. (ii) Jacobians of compositions multiply, outermost on the left. (iii) For a scalar output, multiplying from the left is a factor $d$ cheaper. (iv) There should be a way to differentiate with respect to a matrix that never asks you to index anything. All four are right.

## 3. Formal development

### 3.1 The derivative is a linear map

**Definition 3.1 (Fréchet derivative).** $\mathbf{f} : \mathbb{R}^n \to \mathbb{R}^m$ is differentiable at $\mathbf{x}$ if there is a linear map $D\mathbf{f}(\mathbf{x}) : \mathbb{R}^n \to \mathbb{R}^m$ such that
$$
\mathbf{f}(\mathbf{x} + \mathbf{h}) = \mathbf{f}(\mathbf{x}) + D\mathbf{f}(\mathbf{x})[\mathbf{h}] + r(\mathbf{h}), \qquad \frac{\|r(\mathbf{h})\|}{\|\mathbf{h}\|} \to 0 \text{ as } \mathbf{h} \to \mathbf{0}. \tag{3.1}
$$
We write $r(\mathbf{h}) = o(\|\mathbf{h}\|)$. The square brackets mean "apply the linear map to." $D\mathbf{f}(\mathbf{x})$ is the derivative; it is a map, not a number.

**Theorem 3.2 (The matrix of the derivative is the Jacobian).** If $D\mathbf{f}(\mathbf{x})$ exists, it is unique, and its matrix in the standard bases is
$$
J(\mathbf{x}) = \begin{pmatrix} \partial f_1/\partial x_1 & \cdots & \partial f_1/\partial x_n \\ \vdots & & \vdots \\ \partial f_m/\partial x_1 & \cdots & \partial f_m/\partial x_n \end{pmatrix}, \qquad J_{ij} = \frac{\partial f_i}{\partial x_j}, \qquad \mathbf{f}(\mathbf{x} + \mathbf{h}) = \mathbf{f}(\mathbf{x}) + J\mathbf{h} + o(\|\mathbf{h}\|). \tag{3.2}
$$

▶ **Try it first.** Put $\mathbf{h} = t\mathbf{e}_j$ into (3.1), divide by $t$, and let $t \to 0$. What do you get on each side?

**Derivation.**

1. With $\mathbf{h} = t\mathbf{e}_j$, (3.1) reads $\mathbf{f}(\mathbf{x} + t\mathbf{e}_j) - \mathbf{f}(\mathbf{x}) = t\, D\mathbf{f}(\mathbf{x})[\mathbf{e}_j] + r(t\mathbf{e}_j)$, using linearity to pull out $t$.
2. Divide by $t$: the left side becomes $\big(\mathbf{f}(\mathbf{x} + t\mathbf{e}_j) - \mathbf{f}(\mathbf{x})\big)/t$, whose $i$-th entry tends to $\partial f_i/\partial x_j$ by the definition of a partial derivative. The right side is $D\mathbf{f}(\mathbf{x})[\mathbf{e}_j] + r(t\mathbf{e}_j)/t$, and $\|r(t\mathbf{e}_j)\|/|t| \to 0$.
3. So $D\mathbf{f}(\mathbf{x})[\mathbf{e}_j]$ is the column of partials $(\partial f_i/\partial x_j)_{i}$. By LA-1 Theorem 3.5, $D\mathbf{f}(\mathbf{x})[\mathbf{e}_j]$ is column $j$ of the matrix of the map. That fixes every column, hence the matrix, hence the map. ∎

■ **What just happened.** Partial derivatives are not the definition of the derivative; they are the *entries* of it. The definition is (3.1). The distinction matters: a function can have all partials and fail (3.1) (standard example: $f(x, y) = xy/(x^2 + y^2)$ at the origin, with $f(0,0) = 0$). When the partials are continuous, (3.1) holds; that is the standing assumption in this module.

**Definition 3.3 (Directional derivative).** $D_{\mathbf{u}}\mathbf{f}(\mathbf{x}) = \lim_{t \to 0} \big(\mathbf{f}(\mathbf{x} + t\mathbf{u}) - \mathbf{f}(\mathbf{x})\big)/t$. Putting $\mathbf{h} = t\mathbf{u}$ in (3.1) and repeating steps 1–2 gives $D_{\mathbf{u}}\mathbf{f}(\mathbf{x}) = D\mathbf{f}(\mathbf{x})[\mathbf{u}] = J\mathbf{u}$. Every directional derivative is one matrix–vector product.

**Scalar functions: the gradient.** Let $f : \mathbb{R}^n \to \mathbb{R}$. Then $Df(\mathbf{x})$ is a linear map $\mathbb{R}^n \to \mathbb{R}$, a **linear functional**, and its matrix $J$ is $1 \times n$: a row. A linear functional on $\mathbb{R}^n$ has the form $\mathbf{h} \mapsto \mathbf{g}^\top\mathbf{h}$ for a unique $\mathbf{g} \in \mathbb{R}^n$ (take $g_i = Df(\mathbf{x})[\mathbf{e}_i]$ and use linearity; this is the Riesz representation theorem of LA-12 in its simplest case, with the dot product of LA-1 as the inner product).

**Definition 3.4 (Gradient).** $\nabla f(\mathbf{x})$ is the unique vector with
$$
Df(\mathbf{x})[\mathbf{h}] = \nabla f(\mathbf{x})^\top \mathbf{h} \quad \text{for all } \mathbf{h}, \qquad \text{so} \qquad \nabla f = J^\top = \begin{pmatrix} \partial f/\partial x_1 \\ \vdots \\ \partial f/\partial x_n \end{pmatrix}. \tag{3.3}
$$

This is why the gradient is a column and the Jacobian of a scalar function is a row: $J$ *is* the functional (an element of the dual space), $\nabla f$ is the *vector that represents it* through the inner product (an element of $\mathbb{R}^n$). They contain the same numbers only because we used the dot product. With a different inner product $\langle \mathbf{u}, \mathbf{v}\rangle_M = \mathbf{u}^\top M \mathbf{v}$, the representative becomes $M^{-1}J^\top$ (check: $\langle M^{-1}J^\top, \mathbf{h}\rangle_M = J\mathbf{h}$). The Jacobian never changes; the gradient does. That is the "natural gradient" of ML-2 and the Fisher metric of ST-1, and it is the honest reason the row/column distinction exists.

**Theorem 3.5 (Steepest ascent).** Let $\nabla f(\mathbf{x}) \ne \mathbf{0}$. Over unit vectors $\mathbf{u}$, the directional derivative $D_{\mathbf{u}}f(\mathbf{x}) = \nabla f(\mathbf{x})^\top\mathbf{u}$ is maximized by $\mathbf{u} = \nabla f(\mathbf{x})/\|\nabla f(\mathbf{x})\|$, with maximum value $\|\nabla f(\mathbf{x})\|$, and minimized by the opposite direction.

*Proof.* By Cauchy–Schwarz (LA-1, Theorem 3.11), $\nabla f^\top \mathbf{u} \le \|\nabla f\|\|\mathbf{u}\| = \|\nabla f\|$, with equality iff $\mathbf{u}$ is a positive multiple of $\nabla f$; the unit such vector is $\nabla f/\|\nabla f\|$. The minimum is the same argument with $-\mathbf{u}$. ∎

Also: $\nabla f$ is perpendicular to the level set through $\mathbf{x}$. If $\boldsymbol{\gamma}(t)$ is a curve in the level set, $f(\boldsymbol{\gamma}(t))$ is constant, so by the chain rule below $0 = \nabla f^\top \boldsymbol{\gamma}'(t)$.

### 3.2 The chain rule: Jacobians multiply

**Theorem 3.6 (Chain rule).** Let $\mathbf{f} : \mathbb{R}^n \to \mathbb{R}^m$ be differentiable at $\mathbf{x}$ and $\mathbf{g} : \mathbb{R}^m \to \mathbb{R}^p$ be differentiable at $\mathbf{y} = \mathbf{f}(\mathbf{x})$. Then $\mathbf{g} \circ \mathbf{f}$ is differentiable at $\mathbf{x}$ and
$$
D(\mathbf{g} \circ \mathbf{f})(\mathbf{x}) = D\mathbf{g}(\mathbf{f}(\mathbf{x})) \circ D\mathbf{f}(\mathbf{x}), \qquad J_{\mathbf{g} \circ \mathbf{f}}(\mathbf{x}) = J_{\mathbf{g}}(\mathbf{f}(\mathbf{x}))\, J_{\mathbf{f}}(\mathbf{x}). \tag{3.4}
$$

▶ **Try it first.** Write $\mathbf{f}(\mathbf{x} + \mathbf{h}) = \mathbf{y} + \mathbf{k}$ with $\mathbf{k} = J_f\mathbf{h} + r_f(\mathbf{h})$, substitute into the expansion of $\mathbf{g}$ at $\mathbf{y}$, and collect the terms that are not $J_g J_f \mathbf{h}$. You have to show they are $o(\|\mathbf{h}\|)$. There are two of them.

**Derivation.** Write $J_f = J_{\mathbf{f}}(\mathbf{x})$, $J_g = J_{\mathbf{g}}(\mathbf{y})$.

1. By (3.1) for $\mathbf{f}$: $\mathbf{k}(\mathbf{h}) := \mathbf{f}(\mathbf{x} + \mathbf{h}) - \mathbf{y} = J_f\mathbf{h} + r_f(\mathbf{h})$ with $\|r_f(\mathbf{h})\| = o(\|\mathbf{h}\|)$.
2. By (3.1) for $\mathbf{g}$ at $\mathbf{y}$: $\mathbf{g}(\mathbf{y} + \mathbf{k}) = \mathbf{g}(\mathbf{y}) + J_g\mathbf{k} + r_g(\mathbf{k})$ with $\|r_g(\mathbf{k})\| = o(\|\mathbf{k}\|)$ and $r_g(\mathbf{0}) = \mathbf{0}$.
3. Substitute 1 into 2 and use linearity of $J_g$:
   $$\mathbf{g}(\mathbf{f}(\mathbf{x} + \mathbf{h})) = \mathbf{g}(\mathbf{y}) + J_g J_f\mathbf{h} + \underbrace{J_g\, r_f(\mathbf{h}) + r_g(\mathbf{k}(\mathbf{h}))}_{\text{must be } o(\|\mathbf{h}\|)}.$$
4. First term: $\|J_g r_f(\mathbf{h})\| \le \|J_g\|_F\,\|r_f(\mathbf{h})\|$ (LA-1, exercise C3(b)), and the right side is $o(\|\mathbf{h}\|)$.
5. Second term, in two steps. (a) $\mathbf{k}$ is at most linear in $\mathbf{h}$: for $\|\mathbf{h}\|$ small enough that $\|r_f(\mathbf{h})\| \le \|\mathbf{h}\|$, $\|\mathbf{k}(\mathbf{h})\| \le \|J_f\|_F\|\mathbf{h}\| + \|\mathbf{h}\| = C\|\mathbf{h}\|$ with $C = \|J_f\|_F + 1$. (b) Given $\varepsilon > 0$, choose $\delta$ so that $\|\mathbf{k}\| < \delta$ implies $\|r_g(\mathbf{k})\| \le (\varepsilon/C)\|\mathbf{k}\|$ (possible because $r_g(\mathbf{k}) = o(\|\mathbf{k}\|)$; trivially true if $\mathbf{k} = \mathbf{0}$). Then for $\|\mathbf{h}\| < \delta/C$ we have $\|\mathbf{k}(\mathbf{h})\| < \delta$, so $\|r_g(\mathbf{k}(\mathbf{h}))\| \le (\varepsilon/C)\, C\|\mathbf{h}\| = \varepsilon\|\mathbf{h}\|$. That is the definition of $r_g(\mathbf{k}(\mathbf{h})) = o(\|\mathbf{h}\|)$.
6. So $\mathbf{g} \circ \mathbf{f}$ satisfies (3.1) at $\mathbf{x}$ with linear part $\mathbf{h} \mapsto J_g J_f \mathbf{h}$. By uniqueness (Theorem 3.2) that linear part is $D(\mathbf{g} \circ \mathbf{f})(\mathbf{x})$, and its matrix is $J_g J_f$. ∎

■ **What just happened.** "The derivative of a composition is the composition of the derivatives." Since composition of linear maps is matrix multiplication (LA-1 §3.3, LA-4), Jacobians multiply. **Order:** the outermost function's Jacobian is on the left, evaluated at the inner function's output. You *compute* $\mathbf{f}$ first and $\mathbf{g}$ second; you *write* $J_g J_f$. This reversal is a persistent source of error (§8). For scalar $g \circ \mathbf{f}$, transpose (3.4): $\nabla(g \circ \mathbf{f})(\mathbf{x}) = J_f^\top \nabla g(\mathbf{f}(\mathbf{x}))$. The transpose of the inner Jacobian carries the outer gradient back. You will see $J^\top$ carrying things backwards in every backprop formula.

### 3.3 The differential method

Writing out Jacobians entry by entry is fine for $\mathbb{R}^2 \to \mathbb{R}^2$ and hopeless for $f(X) = \log\det X$. The differential method keeps everything in matrix form.

**Definition 3.7 (Differential).** For a differentiable $F$ (scalar-, vector-, or matrix-valued, of a vector or matrix argument $X$), the **differential** at $X$ is the linear part of the change,
$$dF = DF(X)[dX],$$
where $dX$ is a free increment of the same shape as $X$. So (3.1) reads $F(X + dX) = F(X) + dF + o(\|dX\|)$. The rule of the game: expand to first order and discard any product of two differentials, because $\|dA\|\|dB\| = o(\|(dA, dB)\|)$.

**Theorem 3.8 (Rules for differentials).** For matrices $A, B$ of compatible sizes, scalars $c$, and an invertible $A$ where $A^{-1}$ appears:

(a) $d(A + B) = dA + dB$ and $d(cA) = c\,dA$ (linearity).

(b) Product rule:
$$
d(AB) = (dA)\,B + A\,(dB). \tag{3.5}
$$

(c) $d(A^\top) = (dA)^\top$ and $d(\operatorname{tr}A) = \operatorname{tr}(dA)$ (transpose and trace are linear).

(d) Inverse:
$$
d(A^{-1}) = -A^{-1}(dA)A^{-1}. \tag{3.6}
$$

(e) Jacobi's formula:
$$
d(\det A) = \det(A)\,\operatorname{tr}(A^{-1}dA). \tag{3.7}
$$

(f) $d(\log\det A) = \operatorname{tr}(A^{-1}dA)$ (from (e)).

(g) $d(\sigma(X)) = \sigma'(X) \odot dX$ for elementwise $\sigma$ (entrywise one-variable chain rule).

(h) Chain rule in differential form: if $Y = F(X)$ and $Z = G(Y)$, then $dZ = DG(Y)[dY]$ with $dY = DF(X)[dX]$. In practice: compute $dY$, substitute it wherever $dY$ appears in $dZ$.

▶ **Try it first.** Derive (b) by expanding $(A + dA)(B + dB)$. Derive (d) from $d(AA^{-1}) = d(I) = 0$ and (b). For (e), recall from LA-6 that $\det A = \sum_j a_{ij}C_{ij}$ (cofactor expansion along row $i$) and that $A^{-1} = C^\top/\det A$ where $C$ is the cofactor matrix.

**Derivation.**

1. (b): $(A + dA)(B + dB) - AB = (dA)B + A(dB) + (dA)(dB)$. The last term is a product of two differentials, $\|(dA)(dB)\|_F \le \|dA\|_F\|dB\|_F$, which is $o$ of the size of the increment. Drop it. ∎
2. (c): $X \mapsto X^\top$ and $X \mapsto \operatorname{tr}X$ are linear maps, and a linear map is its own derivative: $L(X + dX) = L(X) + L(dX)$ exactly, with zero remainder. ∎
3. (d): $AA^{-1} = I$ for every invertible $A$, so differentiating both sides, by (b): $(dA)A^{-1} + A\,d(A^{-1}) = 0$. Multiply on the left by $A^{-1}$ and solve: $d(A^{-1}) = -A^{-1}(dA)A^{-1}$. ∎ (Sanity check in one dimension: $d(1/a) = -da/a^2$.)
4. (e): In the cofactor expansion $\det A = \sum_j a_{ij}C_{ij}$ along row $i$, the cofactor $C_{ij}$ is a determinant of a submatrix that omits row $i$, so it does not contain $a_{ij}$. Therefore $\partial(\det A)/\partial a_{ij} = C_{ij}$. Summing over all entries, $d(\det A) = \sum_{i,j} C_{ij}\,da_{ij} = \operatorname{tr}(C^\top dA)$ (the identity $\sum_{ij}G_{ij}dX_{ij} = \operatorname{tr}(G^\top dX)$ is checked in Theorem 3.9 below). By the adjugate formula of LA-6, $C^\top = \det(A)A^{-1}$. So $d(\det A) = \det(A)\operatorname{tr}(A^{-1}dA)$. ∎
5. (f): $d\log u = du/u$ with $u = \det A$. ∎
6. (g): entry by entry, $d(\sigma(X_{ij})) = \sigma'(X_{ij})\,dX_{ij}$. ∎

■ **What just happened.** Every rule is the product rule or linearity, applied without ever writing an index except in Jacobi's formula, where LA-6 did the index work for us. Note what (b) does **not** say: $d(A^2) = (dA)A + A(dA)$, which is *not* $2A\,dA$ unless $A$ and $dA$ commute. Order is preserved in every product.

**Trace identities you will use constantly** (each is a one-line index check, exercise A3): $\operatorname{tr}(AB) = \operatorname{tr}(BA)$ (cyclic); $\operatorname{tr}(A^\top) = \operatorname{tr}A$; a scalar is its own trace, so $\mathbf{u}^\top\mathbf{v} = \operatorname{tr}(\mathbf{u}^\top\mathbf{v}) = \operatorname{tr}(\mathbf{v}\mathbf{u}^\top)$; and for elementwise products
$$\operatorname{tr}\big(A^\top(B \odot C)\big) = \sum_{ij}A_{ij}B_{ij}C_{ij} = \operatorname{tr}\big((A \odot B)^\top C\big).$$

**Theorem 3.9 (Identification rule).** Let $f$ be a differentiable scalar function of a matrix $X \in \mathbb{R}^{m \times n}$. If
$$
df = \operatorname{tr}(G^\top dX) \quad \text{for all } dX, \qquad \text{then} \qquad \nabla_X f = G, \quad \text{i.e. } \frac{\partial f}{\partial X_{ij}} = G_{ij}. \tag{3.8}
$$
Such a $G$ always exists and is unique. For a vector argument, $df = \mathbf{g}^\top d\mathbf{x}$ gives $\nabla f = \mathbf{g}$; for a vector-valued function, $d\mathbf{f} = M\,d\mathbf{x}$ gives $J = M$.

*Proof.* Compute the trace: $(G^\top dX)_{jj} = \sum_i G_{ij}dX_{ij}$, so $\operatorname{tr}(G^\top dX) = \sum_{ij}G_{ij}\,dX_{ij}$. Take $dX = tE_{ij}$ (the matrix with a single $1$ in position $(i,j)$): then $f(X + tE_{ij}) - f(X) = tG_{ij} + o(t)$, so $\partial f/\partial X_{ij} = G_{ij}$. Uniqueness: if $\operatorname{tr}((G - G')^\top dX) = 0$ for all $dX$, take $dX = G - G'$ to get $\|G - G'\|_F^2 = 0$. Existence: $df = Df(X)[dX]$ is a linear functional on $\mathbb{R}^{m \times n}$, which with the Frobenius inner product $\langle G, X\rangle = \operatorname{tr}(G^\top X)$ (LA-1, C3) is an inner-product space, so Riesz gives $G$. The vector cases are the same with $\operatorname{tr}(\mathbf{g}^\top d\mathbf{x}) = \mathbf{g}^\top d\mathbf{x}$ and with $m$ scalar functions stacked. ∎

**Why this avoids the layout wars.** Textbooks disagree about whether $\partial\mathbf{y}/\partial\mathbf{x}$ should be the $m \times n$ Jacobian ("numerator layout") or its transpose ("denominator layout"), and the *Matrix Cookbook* uses one while Wikipedia's tables use both. The fight is about how to arrange the numbers of a linear map in a rectangle. The differential $d\mathbf{y} = M\,d\mathbf{x}$ does not arrange anything: it states the linear map itself, and $M$ multiplying $d\mathbf{x}$ from the left is the Jacobian in our convention by (3.2). For scalars, $df = \mathbf{g}^\top d\mathbf{x}$ names the vector $\mathbf{g}$ with no room for a transpose to hide. The one convention we do fix, that $\nabla_X f$ has the shape of $X$, is forced by (3.8). Everything downstream of this point uses (3.5)–(3.8) and nothing else.

### 3.4 The catalogue, each derived

Each line below is: expand $df$ with Theorem 3.8, rearrange into $\mathbf{g}^\top d\mathbf{x}$ or $\operatorname{tr}(G^\top dX)$ using the trace identities, read off with (3.8). Do each one yourself before reading the right-hand column; they are two or three lines each, and the project checks every one against finite differences.

**(3.9)**
$$
\nabla_{\mathbf{x}}(\mathbf{a}^\top\mathbf{x}) = \mathbf{a}. \tag{3.9}
$$
*Derivation.* $d(\mathbf{a}^\top\mathbf{x}) = \mathbf{a}^\top d\mathbf{x}$ (linearity; $\mathbf{a}$ is constant). Read off $\mathbf{g} = \mathbf{a}$. ∎

**(3.10)**
$$
\nabla_{\mathbf{x}}(\mathbf{x}^\top A\mathbf{x}) = (A + A^\top)\mathbf{x}, \qquad = 2A\mathbf{x} \text{ when } A = A^\top. \tag{3.10}
$$
*Derivation.* $d(\mathbf{x}^\top A\mathbf{x}) = (d\mathbf{x})^\top A\mathbf{x} + \mathbf{x}^\top A\,d\mathbf{x}$ by (3.5). The first term is a scalar, so it equals its transpose: $(d\mathbf{x})^\top A\mathbf{x} = \mathbf{x}^\top A^\top d\mathbf{x}$. Sum: $\mathbf{x}^\top(A^\top + A)\,d\mathbf{x} = \big((A + A^\top)\mathbf{x}\big)^\top d\mathbf{x}$. ∎

**(3.11)**
$$
\nabla_{\mathbf{x}}\|\mathbf{x}\|^2 = 2\mathbf{x}. \tag{3.11}
$$
*Derivation.* (3.10) with $A = I$. ∎

**(3.12)**
$$
\nabla_{\mathbf{x}}\|A\mathbf{x} - \mathbf{b}\|^2 = 2A^\top(A\mathbf{x} - \mathbf{b}). \tag{3.12}
$$
*Derivation.* Let $\mathbf{r} = A\mathbf{x} - \mathbf{b}$, so $d\mathbf{r} = A\,d\mathbf{x}$ (linearity). Then $d(\mathbf{r}^\top\mathbf{r}) = 2\mathbf{r}^\top d\mathbf{r}$ by (3.11) and rule (h), $= 2\mathbf{r}^\top A\,d\mathbf{x} = (2A^\top\mathbf{r})^\top d\mathbf{x}$. Setting this to zero gives $A^\top A\mathbf{x} = A^\top\mathbf{b}$: **the normal equations of LA-5 are $\nabla = \mathbf{0}$.** Least squares is the projection picture (LA-5) and the calculus picture at once. ∎

**(3.13)**
$$
J_{\mathbf{x} \mapsto A\mathbf{x}} = A. \tag{3.13}
$$
*Derivation.* $d(A\mathbf{x}) = A\,d\mathbf{x}$; read off $M = A$ from Theorem 3.9. (A linear map is its own derivative at every point.) ∎

**(3.14)**
$$
\nabla_A(\mathbf{x}^\top A\mathbf{y}) = \mathbf{x}\mathbf{y}^\top. \tag{3.14}
$$
*Derivation.* $d(\mathbf{x}^\top A\mathbf{y}) = \mathbf{x}^\top(dA)\mathbf{y} = \operatorname{tr}(\mathbf{x}^\top dA\,\mathbf{y}) = \operatorname{tr}(\mathbf{y}\mathbf{x}^\top dA)$ (cyclic) $= \operatorname{tr}\big((\mathbf{x}\mathbf{y}^\top)^\top dA\big)$. Read off $G = \mathbf{x}\mathbf{y}^\top$. Shape check: $A$ is $m \times n$, $\mathbf{x} \in \mathbb{R}^m$, $\mathbf{y} \in \mathbb{R}^n$, $\mathbf{x}\mathbf{y}^\top$ is $m \times n$. ✓ ∎

**(3.15)**
$$
\nabla_A\operatorname{tr}(AB) = B^\top. \tag{3.15}
$$
*Derivation.* $d\operatorname{tr}(AB) = \operatorname{tr}(dA\,B) = \operatorname{tr}(B\,dA) = \operatorname{tr}\big((B^\top)^\top dA\big)$. ∎

**(3.16)**
$$
\nabla_A\operatorname{tr}(A^\top BA) = (B + B^\top)A. \tag{3.16}
$$
*Derivation.* $d\operatorname{tr}(A^\top BA) = \operatorname{tr}\big((dA)^\top BA\big) + \operatorname{tr}(A^\top B\,dA)$. First term: $\operatorname{tr}\big((dA)^\top BA\big) = \operatorname{tr}\big((BA)^\top dA\big)$ (transpose inside the trace). Second: $\operatorname{tr}(A^\top B\,dA) = \operatorname{tr}\big((B^\top A)^\top dA\big)$. Sum: $\operatorname{tr}\big((BA + B^\top A)^\top dA\big)$. ∎

**(3.17)**
$$
\nabla_A\log\det A = A^{-\top} := (A^{-1})^\top. \tag{3.17}
$$
*Derivation.* $d\log\det A = \operatorname{tr}(A^{-1}dA) = \operatorname{tr}\big((A^{-\top})^\top dA\big)$ by rule (f). ∎

**(3.18)**
$$
\nabla_A\|A\|_F^2 = 2A. \tag{3.18}
$$
*Derivation.* $\|A\|_F^2 = \operatorname{tr}(A^\top A)$; $d = \operatorname{tr}\big((dA)^\top A\big) + \operatorname{tr}(A^\top dA) = 2\operatorname{tr}(A^\top dA)$. ∎

**(3.19)**
$$
J_{\sigma}(\mathbf{x}) = \operatorname{diag}(\sigma'(\mathbf{x})), \qquad \nabla_{\mathbf{x}}\big(\mathbf{a}^\top\sigma(\mathbf{x})\big) = \operatorname{diag}(\sigma'(\mathbf{x}))\,\mathbf{a} = \sigma'(\mathbf{x}) \odot \mathbf{a}. \tag{3.19}
$$
*Derivation.* $d\sigma(\mathbf{x}) = \sigma'(\mathbf{x}) \odot d\mathbf{x} = \operatorname{diag}(\sigma'(\mathbf{x}))\,d\mathbf{x}$ by rule (g), which is the Jacobian claim by Theorem 3.9. Then $d(\mathbf{a}^\top\sigma(\mathbf{x})) = \mathbf{a}^\top\operatorname{diag}(\sigma'(\mathbf{x}))\,d\mathbf{x} = \big(\operatorname{diag}(\sigma'(\mathbf{x}))\mathbf{a}\big)^\top d\mathbf{x}$, using that a diagonal matrix is symmetric. ∎

**(3.20) Jacobian of softmax.** $\mathbf{s} = \operatorname{softmax}(\mathbf{z})$, $s_i = e^{z_i}/\sum_k e^{z_k}$. Then
$$
J_{\operatorname{softmax}}(\mathbf{z}) = \operatorname{diag}(\mathbf{s}) - \mathbf{s}\mathbf{s}^\top. \tag{3.20}
$$

▶ **Try it first.** Differentiate $\log s_i = z_i - \log\sum_k e^{z_k}$ instead of $s_i$. It is shorter.

**Derivation.**
1. $d(\log s_i) = dz_i - \dfrac{\sum_k e^{z_k}dz_k}{\sum_k e^{z_k}} = dz_i - \sum_k s_k\,dz_k = dz_i - \mathbf{s}^\top d\mathbf{z}$.
2. $ds_i = s_i\,d(\log s_i) = s_i\,dz_i - s_i\,(\mathbf{s}^\top d\mathbf{z})$.
3. Stack over $i$: $d\mathbf{s} = \mathbf{s} \odot d\mathbf{z} - \mathbf{s}(\mathbf{s}^\top d\mathbf{z}) = \operatorname{diag}(\mathbf{s})\,d\mathbf{z} - \mathbf{s}\mathbf{s}^\top d\mathbf{z}$. Read off the Jacobian. ∎

■ **What just happened.** The Jacobian is symmetric, its rows sum to zero ($\mathbf{1}^\top J = \mathbf{s}^\top - (\mathbf{1}^\top\mathbf{s})\mathbf{s}^\top = \mathbf{0}$ because $\sum s_k = 1$), and $J\mathbf{1} = \mathbf{s} - \mathbf{s} = \mathbf{0}$, so it is singular with null vector $\mathbf{1}$: adding a constant to every $z_i$ does not change the softmax. Its rank is $n - 1$. Example 4.2 writes it out.

**(3.21) Cross-entropy with softmax.** For a target $\mathbf{y}$ with $\mathbf{1}^\top\mathbf{y} = 1$ (one-hot or any probability vector) and $L = -\mathbf{y}^\top\log\mathbf{s}$,
$$
\nabla_{\mathbf{z}}L = \mathbf{s} - \mathbf{y}. \tag{3.21}
$$
*Derivation.* $dL = -\sum_i y_i\,d(\log s_i) = -\sum_i y_i\big(dz_i - \mathbf{s}^\top d\mathbf{z}\big) = -\mathbf{y}^\top d\mathbf{z} + (\mathbf{1}^\top\mathbf{y})\,\mathbf{s}^\top d\mathbf{z} = (\mathbf{s} - \mathbf{y})^\top d\mathbf{z}$. ∎ This is the famous simplification: the log cancels the exponential, and the softmax Jacobian never has to be formed. It is also why softmax and cross-entropy are always implemented as one fused operation (ML-3, ML-7).

### 3.5 Hessians, Taylor, Newton, Gauss–Newton

**Definition 3.10 (Hessian).** For $f : \mathbb{R}^n \to \mathbb{R}$ twice differentiable, the Hessian $H = \nabla^2 f(\mathbf{x})$ is the Jacobian of the gradient: $d(\nabla f) = H\,d\mathbf{x}$, so $H_{ij} = \partial(\nabla f)_i/\partial x_j = \partial^2 f/\partial x_j\partial x_i$.

**Theorem 3.11 (Schwarz/Clairaut).** If the second partials of $f$ are continuous, then $H$ is symmetric.

*Proof.* Fix $i \ne j$ and $\mathbf{x}$, and for small $h, k$ let $\Delta = f(\mathbf{x} + h\mathbf{e}_i + k\mathbf{e}_j) - f(\mathbf{x} + h\mathbf{e}_i) - f(\mathbf{x} + k\mathbf{e}_j) + f(\mathbf{x})$. Put $\phi(s) = f(\mathbf{x} + s\mathbf{e}_i + k\mathbf{e}_j) - f(\mathbf{x} + s\mathbf{e}_i)$, so $\Delta = \phi(h) - \phi(0)$. By the mean value theorem, $\Delta = h\,\phi'(\xi) = h\big[\partial_i f(\mathbf{x} + \xi\mathbf{e}_i + k\mathbf{e}_j) - \partial_i f(\mathbf{x} + \xi\mathbf{e}_i)\big]$ for some $\xi$ between $0$ and $h$; by the mean value theorem again in the $j$ direction, $\Delta = hk\,\partial_j\partial_i f(\mathbf{x} + \xi\mathbf{e}_i + \eta\mathbf{e}_j)$ for some $\eta$ between $0$ and $k$. The same $\Delta$ is also $\psi(k) - \psi(0)$ with $\psi(s) = f(\mathbf{x} + h\mathbf{e}_i + s\mathbf{e}_j) - f(\mathbf{x} + s\mathbf{e}_j)$, and the same two steps in the other order give $\Delta = hk\,\partial_i\partial_j f(\mathbf{x} + \xi'\mathbf{e}_i + \eta'\mathbf{e}_j)$. Divide both expressions by $hk$ and let $h, k \to 0$: by continuity both sides tend to the second partials at $\mathbf{x}$, so $\partial_j\partial_i f(\mathbf{x}) = \partial_i\partial_j f(\mathbf{x})$. ∎

**Theorem 3.12 (Hessian of a quadratic form).** For $f(\mathbf{x}) = \tfrac12\mathbf{x}^\top A\mathbf{x}$,
$$
\nabla f = \tfrac12(A + A^\top)\mathbf{x}, \qquad H = \tfrac12(A + A^\top), \tag{3.22}
$$
so for symmetric $A$, $\nabla f = A\mathbf{x}$ and $H = A$: the Hessian of a quadratic form is its (symmetrized) matrix.

*Proof.* $\nabla f$ is (3.10) halved. Then $d(\nabla f) = \tfrac12(A + A^\top)\,d\mathbf{x}$ by linearity, and Theorem 3.9 identifies $H$. ∎

**Theorem 3.13 (Second-order Taylor).** If $f$ is $C^2$,
$$
f(\mathbf{x} + \mathbf{h}) = f(\mathbf{x}) + \nabla f(\mathbf{x})^\top\mathbf{h} + \tfrac12\mathbf{h}^\top H(\mathbf{x})\mathbf{h} + o(\|\mathbf{h}\|^2). \tag{3.22'}
$$

▶ **Try it first.** Reduce to one variable: $\varphi(t) = f(\mathbf{x} + t\mathbf{h})$. Compute $\varphi'$ and $\varphi''$ with the chain rule (3.4), then use the one-variable Taylor theorem you already know.

**Derivation.**
1. The inner map $t \mapsto \mathbf{x} + t\mathbf{h}$ has Jacobian $\mathbf{h}$ (an $n \times 1$ column). By (3.4), $\varphi'(t) = J_f(\mathbf{x} + t\mathbf{h})\,\mathbf{h} = \nabla f(\mathbf{x} + t\mathbf{h})^\top\mathbf{h}$.
2. Apply (3.4) again to $t \mapsto \nabla f(\mathbf{x} + t\mathbf{h})$, whose Jacobian is $H(\mathbf{x} + t\mathbf{h})\,\mathbf{h}$ by Definition 3.10: $\varphi''(t) = \mathbf{h}^\top H(\mathbf{x} + t\mathbf{h})\,\mathbf{h}$.
3. One-variable Taylor with Lagrange remainder: $\varphi(1) = \varphi(0) + \varphi'(0) + \tfrac12\varphi''(\theta)$ for some $\theta \in (0, 1)$. That is $f(\mathbf{x} + \mathbf{h}) = f(\mathbf{x}) + \nabla f(\mathbf{x})^\top\mathbf{h} + \tfrac12\mathbf{h}^\top H(\mathbf{x} + \theta\mathbf{h})\mathbf{h}$.
4. By continuity of $H$, $H(\mathbf{x} + \theta\mathbf{h}) = H(\mathbf{x}) + E(\mathbf{h})$ with $\|E\| \to 0$, and $|\mathbf{h}^\top E\mathbf{h}| \le \|E\|_F\|\mathbf{h}\|^2 = o(\|\mathbf{h}\|^2)$. ∎

**Theorem 3.14 (Second-derivative test).** Let $\nabla f(\mathbf{x}) = \mathbf{0}$ and $H = H(\mathbf{x})$ with eigenvalues $\lambda_1, \dots, \lambda_n$ (real, by Theorem 3.11 and LA-8). If all $\lambda_i > 0$, $\mathbf{x}$ is a strict local minimum; if all $\lambda_i < 0$, a strict local maximum; if some are positive and some negative, a saddle; if some are zero and the rest have one sign, the test is inconclusive.

*Proof.* By the spectral theorem (LA-8), $H = Q\Lambda Q^\top$ with $Q$ orthogonal, so with $\mathbf{c} = Q^\top\mathbf{h}$ (and $\|\mathbf{c}\| = \|\mathbf{h}\|$), $\mathbf{h}^\top H\mathbf{h} = \sum_i\lambda_i c_i^2$. If all $\lambda_i \ge \lambda_{\min} > 0$, then by (3.22′) $f(\mathbf{x} + \mathbf{h}) - f(\mathbf{x}) \ge \tfrac12\lambda_{\min}\|\mathbf{h}\|^2 + o(\|\mathbf{h}\|^2) > 0$ for all small $\mathbf{h} \ne \mathbf{0}$. The maximum case is the same with signs flipped. If $\lambda_i > 0 > \lambda_j$, moving along $\mathbf{q}_i$ increases $f$ and along $\mathbf{q}_j$ decreases it. With a zero eigenvalue the quadratic term vanishes along $\mathbf{q}_i$ and the $o(\|\mathbf{h}\|^2)$ term decides: $x^4$, $-x^4$ and $x^3$ all have $H = 0$ at the origin with different behavior. ∎

The eigenvectors of $H$ are the principal directions of curvature and the eigenvalues are the curvatures; the level sets near a minimum are ellipses with axes along the $\mathbf{q}_i$ and lengths $\propto 1/\sqrt{\lambda_i}$. This is the LA-8 picture of a quadratic form, now attached to every smooth function at every critical point.

**Newton's method.** Replace $f$ near $\mathbf{x}$ by its quadratic model $m(\mathbf{h}) = f + \nabla f^\top\mathbf{h} + \tfrac12\mathbf{h}^\top H\mathbf{h}$ and go to the model's critical point. By (3.9) and (3.22) (with $H$ symmetric), $\nabla_{\mathbf{h}}m = \nabla f + H\mathbf{h}$, which vanishes at
$$
\mathbf{h}_{\text{Newton}} = -H^{-1}\nabla f(\mathbf{x}), \qquad \mathbf{x} \leftarrow \mathbf{x} + \mathbf{h}_{\text{Newton}}. \tag{3.23}
$$
It is the minimizer of $m$ iff $H \succ 0$. For a quadratic $f$, $m = f$ exactly, so one step lands on the exact minimizer (Example 4.1). Compare gradient descent $\mathbf{x} \leftarrow \mathbf{x} - \eta\nabla f$: Newton replaces the scalar step $\eta$ by the matrix $H^{-1}$, which rescales each principal direction by its own curvature. ML-2 proves the convergence rates; the cost is forming and solving with $H$, which is $O(n^3)$ per step.

**Gauss–Newton.** Nonlinear least squares: $f(\mathbf{x}) = \tfrac12\|\mathbf{r}(\mathbf{x})\|^2$ with residual $\mathbf{r} : \mathbb{R}^n \to \mathbb{R}^m$ and Jacobian $J_r$.

▶ **Try it first.** Compute $df$ with (3.11) and rule (h), identify $\nabla f$, then differentiate $\nabla f$ once more with the product rule. There are two terms; which one needs second derivatives of $\mathbf{r}$?

**Derivation.**
1. $df = \tfrac12\,d(\mathbf{r}^\top\mathbf{r}) = \mathbf{r}^\top d\mathbf{r} = \mathbf{r}^\top J_r\,d\mathbf{x} = (J_r^\top\mathbf{r})^\top d\mathbf{x}$. So $\nabla f = J_r^\top\mathbf{r}$.
2. $d(\nabla f) = (dJ_r)^\top\mathbf{r} + J_r^\top d\mathbf{r} = (dJ_r)^\top\mathbf{r} + J_r^\top J_r\,d\mathbf{x}$.
3. Row $i$ of $J_r$ is $\nabla r_i^\top$, so $(dJ_r)^\top\mathbf{r} = \sum_i r_i\,d(\nabla r_i) = \sum_i r_i\,\nabla^2 r_i\,d\mathbf{x}$. Therefore
   $$H = J_r^\top J_r + \sum_{i=1}^m r_i\,\nabla^2 r_i(\mathbf{x}).$$
4. Gauss–Newton **drops the sum**: $H \approx J_r^\top J_r$. The dropped term is small when the residuals $r_i$ are small at the solution (a good fit) or when $\mathbf{r}$ is nearly linear ($\nabla^2 r_i \approx 0$). Newton's step (3.23) becomes
   $$
   (J_r^\top J_r)\,\mathbf{h} = -J_r^\top\mathbf{r}. \tag{3.25}
   $$
   ∎

■ **What just happened.** (3.25) is the normal equations (LA-5, and (3.12) above) for the *linear* least-squares problem $\min_{\mathbf{h}}\|\mathbf{r} + J_r\mathbf{h}\|^2$. Gauss–Newton is "linearize the residual, solve linear least squares, move, repeat," exactly Gauss's 1809 procedure. The approximation buys two things: no second derivatives, and $J_r^\top J_r \succeq 0$ always, so the step is a descent direction whenever $J_r$ has full column rank. Written with the pseudoinverse of LA-9, $\mathbf{h} = -J_r^+\mathbf{r}$. The failure mode is overshooting when the linearization is poor; adding $\lambda I$ to $J_r^\top J_r$ (Levenberg–Marquardt, ML-2) fixes that. Example 4.5 runs the undamped version on a decay-rate fit.

### 3.6 Compositions: forward mode, reverse mode, and why backprop exists

Let $f = f_k \circ f_{k-1} \circ \cdots \circ f_1$ with $\mathbf{x}_0 = \mathbf{x} \in \mathbb{R}^d$, $\mathbf{x}_\ell = f_\ell(\mathbf{x}_{\ell-1}) \in \mathbb{R}^d$ for $\ell < k$, and a scalar final layer $f_k : \mathbb{R}^d \to \mathbb{R}$. By the chain rule (3.4) applied $k - 1$ times,
$$J_f(\mathbf{x}) = J_k\,J_{k-1}\cdots J_2\,J_1, \qquad J_k \in \mathbb{R}^{1 \times d},\quad J_\ell \in \mathbb{R}^{d \times d}\ (\ell < k),$$
each $J_\ell$ evaluated at $\mathbf{x}_{\ell-1}$. Multiplying a $p \times q$ matrix by a $q \times r$ matrix costs $2pqr$ flops.

▶ **Try it first.** Count the flops of evaluating this product right-to-left (start with $J_2J_1$) and left-to-right (start with $J_kJ_{k-1}$).

**Derivation (the cost comparison).**
1. **Right to left.** $J_2J_1$ costs $2d^3$; so does each of the next $k - 3$ products, since every partial product is $d \times d$. The final product with the $1 \times d$ row costs $2d^2$. Total $\approx 2(k - 2)d^3 + 2d^2 \sim 2kd^3$.
2. **Left to right.** $J_kJ_{k-1}$ is a $1 \times d$ row times a $d \times d$ matrix: $2d^2$. The result is again a $1 \times d$ row, so every subsequent product is also $2d^2$. Total $2(k - 1)d^2 \sim 2kd^2$. **A factor of $d$ cheaper.** This is **reverse mode**: starting from the output, propagate a row vector backwards through the Jacobians. Each step is a **vector–Jacobian product** (VJP) $\mathbf{v}^\top J_\ell$, equivalently $J_\ell^\top\mathbf{v}$.
3. **Forward mode.** The other natural primitive is a **Jacobian–vector product** (JVP) $J_f\mathbf{v} = J_k(J_{k-1}(\cdots(J_1\mathbf{v})))$ for a chosen direction $\mathbf{v}$: right-to-left with a column, each product $d \times d$ times $d \times 1$, total $2kd^2$. That gives one directional derivative. To get the whole gradient you need $J_f\mathbf{e}_j$ for $j = 1, \dots, d$: $d$ passes, $2kd^3$ in total, the same as step 1.
4. **General rule.** For $f : \mathbb{R}^n \to \mathbb{R}^m$, the full Jacobian costs $n$ JVPs (forward) or $m$ VJPs (reverse). Pick the smaller side. For a scalar loss, $m = 1$: one reverse pass, and reverse wins by a factor $n$. For one parameter and many outputs (sensitivity of an ODE solution to one constant, $n = 1$), forward wins. ∎

■ **What just happened.** Backpropagation exists because the loss is a scalar. In practice nobody forms a $J_\ell$: the VJP through a layer is computed from the layer's formula, and for $\mathbf{z} = W\mathbf{a}$ the VJP is $\boldsymbol{\delta} \mapsto W^\top\boldsymbol{\delta}$, which costs the same $2hd$ flops as the forward product $W\mathbf{a}$. This is the content of the **cheap gradient principle** (Griewank and Walther, *Evaluating Derivatives*): the gradient of a scalar function costs at most a small constant (their bound is between 3 and 4) times the cost of the function, regardless of the number of inputs. The price is memory: the VJP at layer $\ell$ needs the forward values $\mathbf{x}_{\ell-1}$ (to evaluate $J_\ell$ there), so reverse mode stores the entire forward pass. Forward mode stores nothing extra. Part 10 of the project times both orders.

Two VJP rules make an autodiff engine. For the loss $L$ and an intermediate $C$, write $G = \nabla_C L$ so that $dL = \operatorname{tr}(G^\top dC)$ (this is (3.8) read backwards: $G$ is what the layers above have already propagated to $C$).

**Matrix product.** If $C = AB$ then, by (3.5), $dC = (dA)B + A(dB)$ and
$$dL = \operatorname{tr}(G^\top dA\,B) + \operatorname{tr}(G^\top A\,dB) = \operatorname{tr}\big((GB^\top)^\top dA\big) + \operatorname{tr}\big((A^\top G)^\top dB\big),$$
using the cyclic rule on the first term, $\operatorname{tr}(G^\top dA\,B) = \operatorname{tr}(BG^\top dA)$ and $BG^\top = (GB^\top)^\top$. By (3.8),
$$
\nabla_A L = G\,B^\top, \qquad \nabla_B L = A^\top G. \tag{3.24}
$$
Shape check: if $A$ is $m \times p$ and $B$ is $p \times n$, then $G$ is $m \times n$, $GB^\top$ is $m \times p$ and $A^\top G$ is $p \times n$. ✓ These two lines *are* the `__matmul__` backward of every autodiff library.

**Broadcast addition.** If $Z = Y + \mathbf{1}\mathbf{b}^\top$ (add the row $\mathbf{b}^\top$ to every row of $Y$, which is what `Y + b` does in NumPy), then $dZ = dY + \mathbf{1}\,d\mathbf{b}^\top$ and $\operatorname{tr}(G^\top\mathbf{1}\,d\mathbf{b}^\top) = \operatorname{tr}(d\mathbf{b}^\top G^\top\mathbf{1}) = (G^\top\mathbf{1})^\top d\mathbf{b}$, so
$$
\nabla_Y L = G, \qquad \nabla_{\mathbf{b}}L = G^\top\mathbf{1} = \text{column sums of } G. \tag{3.26}
$$
Broadcasting copies $\mathbf{b}$ into every row; its adjoint sums over the rows. This is the single most common shape bug in hand-written backprop (§8).

### 3.7 Backpropagation for a two-layer network, completely

The network, one input $\mathbf{x} \in \mathbb{R}^d$, hidden width $h$, output dimension $m$, target $\mathbf{y} \in \mathbb{R}^m$:
$$
\mathbf{z}_1 = W_1\mathbf{x} + \mathbf{b}_1,\quad \mathbf{a}_1 = \sigma(\mathbf{z}_1),\quad \mathbf{z}_2 = W_2\mathbf{a}_1 + \mathbf{b}_2,\quad L = \tfrac12\|\mathbf{z}_2 - \mathbf{y}\|^2,
$$
with shapes $W_1 \in \mathbb{R}^{h \times d}$, $\mathbf{b}_1, \mathbf{z}_1, \mathbf{a}_1 \in \mathbb{R}^h$, $W_2 \in \mathbb{R}^{m \times h}$, $\mathbf{b}_2, \mathbf{z}_2 \in \mathbb{R}^m$, $L \in \mathbb{R}$. We want $\nabla_{W_1}L, \nabla_{\mathbf{b}_1}L, \nabla_{W_2}L, \nabla_{\mathbf{b}_2}L$.

▶ **Try it first.** Compute $dL$ in terms of $d\mathbf{z}_2$. Then expand $d\mathbf{z}_2$ in terms of $dW_2$, $d\mathbf{a}_1$, $d\mathbf{b}_2$ and identify two of the four gradients. Then push $d\mathbf{a}_1$ down to $d\mathbf{z}_1$ and $d\mathbf{z}_1$ down to $dW_1$, $d\mathbf{b}_1$. Write the shape next to every object.

**Derivation.** Every step is one of (3.5), (3.8), (3.14), (3.19), (3.24).

1. **Top.** $dL = (\mathbf{z}_2 - \mathbf{y})^\top d\mathbf{z}_2$ by (3.11) and rule (h). Define
   $$\boldsymbol{\delta}_2 := \nabla_{\mathbf{z}_2}L = \mathbf{z}_2 - \mathbf{y} \in \mathbb{R}^m. \tag{3.28}$$
2. **Layer 2.** $d\mathbf{z}_2 = (dW_2)\mathbf{a}_1 + W_2\,d\mathbf{a}_1 + d\mathbf{b}_2$ by (3.5) (hold $\mathbf{x}$ fixed; the data are not variables). So
   $$dL = \boldsymbol{\delta}_2^\top(dW_2)\mathbf{a}_1 + \boldsymbol{\delta}_2^\top W_2\,d\mathbf{a}_1 + \boldsymbol{\delta}_2^\top d\mathbf{b}_2.$$
   - First term: by (3.14) with $\mathbf{x} \to \boldsymbol{\delta}_2$, $\mathbf{y} \to \mathbf{a}_1$, it is $\operatorname{tr}\big((\boldsymbol{\delta}_2\mathbf{a}_1^\top)^\top dW_2\big)$:
     $$\nabla_{W_2}L = \boldsymbol{\delta}_2\,\mathbf{a}_1^\top \in \mathbb{R}^{m \times h}. \tag{3.29}$$
   - Third term: $$\nabla_{\mathbf{b}_2}L = \boldsymbol{\delta}_2 \in \mathbb{R}^m. \tag{3.30}$$
   - Second term: $\boldsymbol{\delta}_2^\top W_2\,d\mathbf{a}_1 = (W_2^\top\boldsymbol{\delta}_2)^\top d\mathbf{a}_1$, so $\nabla_{\mathbf{a}_1}L = W_2^\top\boldsymbol{\delta}_2 \in \mathbb{R}^h$. The transpose carries the error signal backwards through the layer, as in §3.2.
3. **Nonlinearity.** $d\mathbf{a}_1 = \sigma'(\mathbf{z}_1) \odot d\mathbf{z}_1$ by rule (g), so $(W_2^\top\boldsymbol{\delta}_2)^\top\big(\sigma'(\mathbf{z}_1) \odot d\mathbf{z}_1\big) = \big((W_2^\top\boldsymbol{\delta}_2) \odot \sigma'(\mathbf{z}_1)\big)^\top d\mathbf{z}_1$ by the $\odot$ trace identity (or (3.19)). Define
   $$\boldsymbol{\delta}_1 := \nabla_{\mathbf{z}_1}L = (W_2^\top\boldsymbol{\delta}_2) \odot \sigma'(\mathbf{z}_1) \in \mathbb{R}^h. \tag{3.31}$$
4. **Layer 1.** $d\mathbf{z}_1 = (dW_1)\mathbf{x} + d\mathbf{b}_1$ (plus $W_1\,d\mathbf{x}$ if you also want $\nabla_{\mathbf{x}}L = W_1^\top\boldsymbol{\delta}_1$). Exactly as in step 2:
   $$\nabla_{W_1}L = \boldsymbol{\delta}_1\,\mathbf{x}^\top \in \mathbb{R}^{h \times d}, \tag{3.32}$$
   $$\nabla_{\mathbf{b}_1}L = \boldsymbol{\delta}_1 \in \mathbb{R}^h. \tag{3.33}$$
   ∎

■ **What just happened.** Only two kinds of object appear: **error signals** $\boldsymbol{\delta}_\ell = \nabla_{\mathbf{z}_\ell}L$, each the shape of the layer's pre-activation, and **gradients** which are outer products of an error signal with the layer's input. The recursion for any depth is
$$\boldsymbol{\delta}_{\ell-1} = (W_\ell^\top\boldsymbol{\delta}_\ell) \odot \sigma'(\mathbf{z}_{\ell-1}), \qquad \nabla_{W_\ell}L = \boldsymbol{\delta}_\ell\,\mathbf{a}_{\ell-1}^\top, \qquad \nabla_{\mathbf{b}_\ell}L = \boldsymbol{\delta}_\ell,$$
with $\mathbf{a}_0 = \mathbf{x}$. Each step is one VJP, costing the same as the forward product; the whole backward pass costs about twice the forward pass, which is the cheap gradient principle in the flesh. Notice that $\sigma'(\mathbf{z}_1)$ multiplies the signal: a saturated $\tanh$ unit ($|z| \gg 1$, $\sigma' \approx 0$) stops the signal dead. Example 4.4 has one such unit; that is the "vanishing gradient" of ML-7, visible in a $2 \times 2$ computation.

**The batch version.** Stack $n$ inputs as the rows of $X \in \mathbb{R}^{n \times d}$ and targets as $Y \in \mathbb{R}^{n \times m}$, keep the same $W_1, W_2$, and average the loss:
$$Z_1 = XW_1^\top + \mathbf{1}\mathbf{b}_1^\top \in \mathbb{R}^{n \times h},\quad A_1 = \sigma(Z_1),\quad Z_2 = A_1W_2^\top + \mathbf{1}\mathbf{b}_2^\top \in \mathbb{R}^{n \times m},\quad L = \tfrac{1}{2n}\|Z_2 - Y\|_F^2.$$
Row $i$ of $Z_1$ is $(W_1\mathbf{x}_i + \mathbf{b}_1)^\top$, so this is the single-sample network applied to every row, transposed. The transposes on $W$ are forced by the row convention: $(W\mathbf{x})^\top = \mathbf{x}^\top W^\top$.

▶ **Try it first.** Repeat steps 1–4 with (3.18), (3.24) and (3.26) in place of (3.11), (3.14) and "$\nabla_{\mathbf{b}} = \boldsymbol{\delta}$." Every $\boldsymbol{\delta}$ becomes an $n$-row matrix $\Delta$.

**Derivation.**
1. $dL = \frac1n\operatorname{tr}\big((Z_2 - Y)^\top dZ_2\big)$ by (3.18), so
   $$\Delta_2 := \nabla_{Z_2}L = (Z_2 - Y)/n \in \mathbb{R}^{n \times m}. \tag{3.34}$$
2. $Z_2 = A_1W_2^\top + \mathbf{1}\mathbf{b}_2^\top$. By (3.24) with $A \to A_1$, $B \to W_2^\top$, $G \to \Delta_2$: $\nabla_{A_1}L = \Delta_2W_2 \in \mathbb{R}^{n \times h}$ and $\nabla_{W_2^\top}L = A_1^\top\Delta_2$, i.e. $\nabla_{W_2}L = \Delta_2^\top A_1 \in \mathbb{R}^{m \times h}$ (the gradient with respect to $W^\top$ is the transpose of the gradient with respect to $W$, since $\operatorname{tr}(G^\top dW^\top) = \operatorname{tr}(G\,dW)$). By (3.26), $\nabla_{\mathbf{b}_2}L = \Delta_2^\top\mathbf{1}$. Collected:
   $$\nabla_{W_2}L = \Delta_2^\top A_1 \in \mathbb{R}^{m \times h}, \qquad \nabla_{\mathbf{b}_2}L = \Delta_2^\top\mathbf{1} \in \mathbb{R}^m, \qquad \nabla_{A_1}L = \Delta_2 W_2 \in \mathbb{R}^{n \times h}. \tag{3.35}$$
3. $dA_1 = \sigma'(Z_1) \odot dZ_1$, so
   $$\Delta_1 := \nabla_{Z_1}L = (\Delta_2W_2) \odot \sigma'(Z_1) \in \mathbb{R}^{n \times h}. \tag{3.36}$$
4. $Z_1 = XW_1^\top + \mathbf{1}\mathbf{b}_1^\top$, so by the same two rules
   $$\nabla_{W_1}L = \Delta_1^\top X \in \mathbb{R}^{h \times d}, \qquad \nabla_{\mathbf{b}_1}L = \Delta_1^\top\mathbf{1} \in \mathbb{R}^h. \tag{3.37}$$
   ∎

Row $i$ of $\Delta_\ell$ is $\boldsymbol{\delta}_\ell^{(i)\top}/n$, and $\Delta_2^\top A_1 = \frac1n\sum_i\boldsymbol{\delta}_2^{(i)}\mathbf{a}_1^{(i)\top}$: the batch gradient is the average of the single-sample outer products, computed as one matrix product. The project's Part 6 implements (3.34)–(3.37) and trains a network with them.

### 3.8 Gradient checking by finite differences

Every derivation above can be wrong by a transpose, and the way to know is to compare with a derivative-free estimate.

**Central difference.** For scalar $f$ and the $i$-th coordinate,
$$\frac{f(\mathbf{x} + h\mathbf{e}_i) - f(\mathbf{x} - h\mathbf{e}_i)}{2h} = \frac{\partial f}{\partial x_i} + \frac{h^2}{6}\frac{\partial^3 f}{\partial x_i^3} + O(h^4),$$
because the one-variable Taylor expansions $f(x \pm h) = f \pm hf' + \tfrac12h^2f'' \pm \tfrac16h^3f''' + \cdots$ cancel the even terms when subtracted. Truncation error $O(h^2)$. The forward difference $(f(\mathbf{x} + h\mathbf{e}_i) - f(\mathbf{x}))/h$ has truncation error $\tfrac12hf'' = O(h)$.

**Choosing $h$.** Each evaluation of $f$ carries rounding error about $\epsilon|f|$ with $\epsilon \approx 2.2 \times 10^{-16}$ in double precision (LA-13). The difference of two such numbers has error $\sim 2\epsilon|f|$, and dividing by $2h$ gives $\epsilon|f|/h$. Total error of the central difference:
$$E(h) \approx \frac{|f'''|}{6}h^2 + \frac{\epsilon|f|}{h}.$$
Minimize: $E'(h) = \tfrac13|f'''|h - \epsilon|f|/h^2 = 0$ gives $h^3 = 3\epsilon|f|/|f'''|$, so for $f$ and $f'''$ of order one, $h \approx \epsilon^{1/3} \approx 6 \times 10^{-6}$, with best achievable error $\sim\epsilon^{2/3} \approx 4 \times 10^{-11}$. The same computation for the forward difference gives $h \approx \sqrt\epsilon \approx 1.5 \times 10^{-8}$ and error $\sim\sqrt\epsilon$. Too small an $h$ and rounding wins (cancellation); too large and truncation wins. The project plots $E(h)$ over twelve decades (`fd_error_vs_h.png`): a V with its bottom at $h \approx 3 \times 10^{-6}$ for central and $10^{-8}$ for forward. Scale $h$ by $\max(1, |x_i|)$ when the coordinates are large.

**The metric.** Compare the analytic gradient $\mathbf{g}$ with the numerical $\hat{\mathbf{g}}$ by the **relative error**
$$
\frac{\|\mathbf{g} - \hat{\mathbf{g}}\|}{\|\mathbf{g}\| + \|\hat{\mathbf{g}}\|}, \tag{3.27}
$$
which is scale-free and does not blow up when one gradient is zero. For smooth functions with central differences, below $10^{-7}$ is a pass, $10^{-5}$ is suspicious, above $10^{-3}$ is a bug. (A kink, as in ReLU, can legitimately produce a large error if $\mathbf{x} \pm h\mathbf{e}_i$ straddle it.) The project's `grad_check` is this; every identity in §3.4 passes it at around $10^{-11}$.

### 3.9 A statistics application: the Gaussian log-likelihood

The density of $\mathbf{x} \in \mathbb{R}^p$ under $\mathcal{N}(\boldsymbol{\mu}, \Sigma)$ with $\Sigma \succ 0$ (F5) gives the log-likelihood
$$\ell(\boldsymbol{\mu}, \Sigma) = -\tfrac{p}{2}\log 2\pi - \tfrac12\log\det\Sigma - \tfrac12(\mathbf{x} - \boldsymbol{\mu})^\top\Sigma^{-1}(\mathbf{x} - \boldsymbol{\mu}).$$

▶ **Try it first.** Let $\mathbf{u} = \mathbf{x} - \boldsymbol{\mu}$. For $\nabla_{\boldsymbol{\mu}}$, note $d\mathbf{u} = -d\boldsymbol{\mu}$ and use (3.10). For $\nabla_\Sigma$, use rule (f) for the log-determinant and (3.6) for the inverse, then the cyclic rule to isolate $d\Sigma$.

**Derivation.**
1. **With respect to $\boldsymbol{\mu}$.** $\Sigma^{-1}$ is symmetric, so by (3.10) $d(\mathbf{u}^\top\Sigma^{-1}\mathbf{u}) = 2\mathbf{u}^\top\Sigma^{-1}d\mathbf{u} = -2\mathbf{u}^\top\Sigma^{-1}d\boldsymbol{\mu}$. Hence $d\ell = \mathbf{u}^\top\Sigma^{-1}d\boldsymbol{\mu} = (\Sigma^{-1}\mathbf{u})^\top d\boldsymbol{\mu}$ and
   $$\nabla_{\boldsymbol{\mu}}\ell = \Sigma^{-1}(\mathbf{x} - \boldsymbol{\mu}). \tag{3.38}$$
2. **With respect to $\Sigma$**, holding $\boldsymbol{\mu}$ fixed. $d\log\det\Sigma = \operatorname{tr}(\Sigma^{-1}d\Sigma)$ by rule (f). And $d(\mathbf{u}^\top\Sigma^{-1}\mathbf{u}) = \mathbf{u}^\top d(\Sigma^{-1})\mathbf{u} = -\mathbf{u}^\top\Sigma^{-1}(d\Sigma)\Sigma^{-1}\mathbf{u}$ by (3.6), $= -\operatorname{tr}\big(\Sigma^{-1}\mathbf{u}\mathbf{u}^\top\Sigma^{-1}d\Sigma\big)$ by the cyclic rule (move $\mathbf{u}^\top\Sigma^{-1}$ from the front to the back).
3. Assemble: $d\ell = -\tfrac12\operatorname{tr}(\Sigma^{-1}d\Sigma) + \tfrac12\operatorname{tr}(\Sigma^{-1}\mathbf{u}\mathbf{u}^\top\Sigma^{-1}d\Sigma) = \operatorname{tr}\big(G^\top d\Sigma\big)$ with $G^\top = \tfrac12(\Sigma^{-1}\mathbf{u}\mathbf{u}^\top\Sigma^{-1} - \Sigma^{-1})$. This matrix is symmetric, so $G = G^\top$:
   $$\nabla_\Sigma\ell = \tfrac12\big(\Sigma^{-1}(\mathbf{x} - \boldsymbol{\mu})(\mathbf{x} - \boldsymbol{\mu})^\top\Sigma^{-1} - \Sigma^{-1}\big). \tag{3.39}$$
   ∎

**The MLE (preview of ST-1 and ML-6).** For an i.i.d. sample $\mathbf{x}_1, \dots, \mathbf{x}_n$ the log-likelihood is the sum, and so are the gradients. Setting $\sum_i\Sigma^{-1}(\mathbf{x}_i - \boldsymbol{\mu}) = \mathbf{0}$ and multiplying by $\Sigma$ gives $\hat{\boldsymbol{\mu}} = \bar{\mathbf{x}} = \frac1n\sum_i\mathbf{x}_i$. With $S = \sum_i(\mathbf{x}_i - \bar{\mathbf{x}})(\mathbf{x}_i - \bar{\mathbf{x}})^\top$, setting $\tfrac12(\Sigma^{-1}S\Sigma^{-1} - n\Sigma^{-1}) = 0$ and multiplying by $\Sigma$ on both sides gives $S = n\Sigma$, so $\hat\Sigma = S/n$, the (biased) sample covariance. ST-1 shows the score $\nabla\ell$ has mean zero and that $-\mathbb{E}[\nabla^2\ell]$ is the Fisher information; ML-6's EM algorithm is this computation with weights. The project recovers $\hat{\boldsymbol{\mu}}, \hat\Sigma$ by gradient ascent using (3.38)–(3.39).

*A subtlety you should know exists.* (3.39) treats all $p^2$ entries of $\Sigma$ as independent. If you instead parametrize by the $p(p+1)/2$ entries on and above the diagonal, the derivative with respect to an off-diagonal entry doubles, because that entry appears twice in $\Sigma$. Stationary points are the same either way, which is why the MLE is unaffected, but if you hand a gradient to an optimizer the two parametrizations differ. Minka's note discusses this; exercise C3 touches it.

### 3.10 The `vec` trick, one example

For matrix-valued functions of matrices, the Jacobian is a linear map on matrices. To write it as a matrix, flatten with $\operatorname{vec}$ (stack the columns; LA-12) and use the Kronecker identity $\operatorname{vec}(AXB) = (B^\top \otimes A)\operatorname{vec}(X)$ from LA-12. Then for $F(X) = AXB$ with $A \in \mathbb{R}^{m \times n}$, $X \in \mathbb{R}^{n \times q}$, $B \in \mathbb{R}^{q \times r}$:
$$
dF = A\,(dX)\,B \quad\Longrightarrow\quad \operatorname{vec}(dF) = (B^\top \otimes A)\operatorname{vec}(dX), \qquad \frac{\partial\operatorname{vec}F}{\partial\operatorname{vec}X} = B^\top \otimes A \in \mathbb{R}^{mr \times nq}. \tag{3.40}
$$
The same move turns (3.6) into the Jacobian of inversion, $-(X^{-\top} \otimes X^{-1})$ (exercise C1), and $d(X^2) = (dX)X + X(dX)$ into $X^\top \otimes I + I \otimes X$. You rarely want these matrices (they are huge), but they exist, and they are what "the derivative is a linear map" means for matrix-to-matrix functions. The project checks (3.40) numerically.

## 4. Worked examples

### Example 4.1: Gradient, Hessian, and one Newton step on a quadratic

$f(\mathbf{x}) = \tfrac12\mathbf{x}^\top A\mathbf{x} - \mathbf{b}^\top\mathbf{x}$ with $A = \begin{pmatrix}2 & 1\\1 & 3\end{pmatrix}$, $\mathbf{b} = \begin{pmatrix}1\\1\end{pmatrix}$.

*Gradient and Hessian.* $A$ is symmetric, so by (3.22) and (3.9), $\nabla f = A\mathbf{x} - \mathbf{b}$ and $H = A$. Eigenvalues of $A$: $\lambda = \tfrac12(5 \pm \sqrt5) \approx 1.382, 3.618$, both positive, so $f$ is a bowl (Theorem 3.14) and its unique minimizer is where $A\mathbf{x} = \mathbf{b}$.

*Newton from $\mathbf{x}_0 = \mathbf{0}$.* $\nabla f(\mathbf{0}) = -\mathbf{b} = (-1, -1)^\top$. By (3.23), $\mathbf{x}_1 = \mathbf{0} - A^{-1}(-\mathbf{b}) = A^{-1}\mathbf{b}$. With $A^{-1} = \tfrac15\begin{pmatrix}3 & -1\\-1 & 2\end{pmatrix}$ (LA-2), $\mathbf{x}_1 = \tfrac15(2, 1)^\top = (0.4, 0.2)^\top$. Check: $A\mathbf{x}_1 = (0.8 + 0.2,\; 0.4 + 0.6)^\top = (1, 1)^\top = \mathbf{b}$. ✓ One step, exact minimizer, as §3.5 promised for a quadratic. $f(\mathbf{x}_1) = \tfrac12\mathbf{x}_1^\top\mathbf{b} - \mathbf{b}^\top\mathbf{x}_1 = -\tfrac12(0.6) = -0.3$.

*Gradient descent for comparison.* $\mathbf{x}_{k+1} = \mathbf{x}_k - \eta(A\mathbf{x}_k - \mathbf{b})$. With $\eta = 1/\lambda_{\max}$ it takes 36 iterations to get within $10^{-8}$ of $\mathbf{x}^*$; with the optimal fixed step $\eta = 2/(\lambda_{\min} + \lambda_{\max})$, 22 iterations, and the error shrinks by the factor $(\lambda_{\max} - \lambda_{\min})/(\lambda_{\max} + \lambda_{\min}) = 1/\sqrt5 \approx 0.447$ per step (ML-2 proves this rate; the project measures it). Figure `newton_vs_gd.png`.

### Example 4.2: The softmax Jacobian, written out

At $\mathbf{s} = (0.2, 0.3, 0.5)^\top$ (any $\mathbf{z}$ with this softmax, e.g. $\mathbf{z} = \log\mathbf{s}$), (3.20) gives
$$J = \operatorname{diag}(\mathbf{s}) - \mathbf{s}\mathbf{s}^\top = \begin{pmatrix}0.2 & 0 & 0\\0 & 0.3 & 0\\0 & 0 & 0.5\end{pmatrix} - \begin{pmatrix}0.04 & 0.06 & 0.10\\0.06 & 0.09 & 0.15\\0.10 & 0.15 & 0.25\end{pmatrix} = \begin{pmatrix}0.16 & -0.06 & -0.10\\-0.06 & 0.21 & -0.15\\-0.10 & -0.15 & 0.25\end{pmatrix}.$$
Every row sums to zero; $J$ is symmetric; $J\mathbf{1} = \mathbf{0}$; rank 2. The diagonal entries $s_i(1 - s_i)$ are the derivative of a single logistic unit, and the off-diagonal $-s_is_j$ says pushing up $z_j$ pulls down $s_i$. With cross-entropy and target $\mathbf{y} = \mathbf{e}_3$, (3.21) gives $\nabla_{\mathbf{z}}L = \mathbf{s} - \mathbf{y} = (0.2, 0.3, -0.5)^\top$: the gradient pushes down the wrong logits in proportion to their probability and pushes up the right one. The project verifies $J$ against the numerical Jacobian.

### Example 4.3: $\nabla\log\det$ for a $2 \times 2$, checked by direct differentiation

For $A = \begin{pmatrix}a & b\\c & d\end{pmatrix}$, $\det A = ad - bc$ and $\log\det A = \log(ad - bc)$. Direct partials:
$$\frac{\partial}{\partial a} = \frac{d}{ad - bc},\quad \frac{\partial}{\partial b} = \frac{-c}{ad - bc},\quad \frac{\partial}{\partial c} = \frac{-b}{ad - bc},\quad \frac{\partial}{\partial d} = \frac{a}{ad - bc},$$
so, arranged in the shape of $A$, $\nabla_A\log\det A = \dfrac{1}{ad - bc}\begin{pmatrix}d & -c\\-b & a\end{pmatrix}$. And (3.17) says $A^{-\top}$: $A^{-1} = \dfrac{1}{ad - bc}\begin{pmatrix}d & -b\\-c & a\end{pmatrix}$, transposed is $\dfrac{1}{ad - bc}\begin{pmatrix}d & -c\\-b & a\end{pmatrix}$. ✓ Note the transpose matters: the $(1,2)$ entry of the gradient is $-c$, not $-b$. Numerically, for $A = \begin{pmatrix}2 & 1\\1 & 3\end{pmatrix}$, $\det = 5$ and $\nabla\log\det A = \begin{pmatrix}0.6 & -0.2\\-0.2 & 0.4\end{pmatrix}$ (symmetric because $A$ is). The project does this with SymPy and with finite differences.

### Example 4.4: Backprop on a two-layer network, every number

Take $d = 2$, $h = 2$, $m = 1$, $\sigma = \tanh$ (so $\sigma'(z) = 1 - \tanh^2 z = 1 - a^2$), and
$$\mathbf{x} = \begin{pmatrix}1\\2\end{pmatrix},\quad W_1 = \begin{pmatrix}1 & 1\\2 & -1\end{pmatrix},\quad \mathbf{b}_1 = \begin{pmatrix}-1\\1\end{pmatrix},\quad W_2 = \begin{pmatrix}2 & -1\end{pmatrix},\quad b_2 = 1,\quad y = 1.$$

*Forward.*
- $\mathbf{z}_1 = W_1\mathbf{x} + \mathbf{b}_1 = (1 + 2 - 1,\; 2 - 2 + 1)^\top = (2, 1)^\top$.
- $\mathbf{a}_1 = \tanh(\mathbf{z}_1) = (0.9640, 0.7616)^\top$.
- $z_2 = W_2\mathbf{a}_1 + b_2 = 2(0.9640) - 0.7616 + 1 = 2.1665$.
- $L = \tfrac12(2.1665 - 1)^2 = \tfrac12(1.1665)^2 = 0.6803$.

*Backward*, with shapes.
- (3.28): $\delta_2 = z_2 - y = 1.1665$ $\;(1 \times 1)$.
- (3.29): $\nabla_{W_2}L = \delta_2\mathbf{a}_1^\top = 1.1665\,(0.9640,\; 0.7616) = (1.1245,\; 0.8884)$ $\;(1 \times 2)$.
- (3.30): $\nabla_{b_2}L = 1.1665$.
- $W_2^\top\delta_2 = (2, -1)^\top(1.1665) = (2.3329, -1.1665)^\top$ $\;(2 \times 1)$.
- $\sigma'(\mathbf{z}_1) = 1 - \mathbf{a}_1^2 = (1 - 0.9293,\; 1 - 0.5800)^\top = (0.0707, 0.4200)^\top$.
- (3.31): $\boldsymbol{\delta}_1 = (2.3329 \times 0.0707,\; -1.1665 \times 0.4200)^\top = (0.1648, -0.4899)^\top$ $\;(2 \times 1)$.
- (3.32): $\nabla_{W_1}L = \boldsymbol{\delta}_1\mathbf{x}^\top = \begin{pmatrix}0.1648\\-0.4899\end{pmatrix}(1\;\;2) = \begin{pmatrix}0.1648 & 0.3296\\-0.4899 & -0.9798\end{pmatrix}$ $\;(2 \times 2)$.
- (3.33): $\nabla_{\mathbf{b}_1}L = (0.1648, -0.4899)^\top$.

Look at the first hidden unit: $z = 2$ puts $\tanh$ near saturation, $\sigma' = 0.07$, and the error signal $2.33$ arriving from above is cut to $0.16$. That is a vanishing gradient in miniature. The project computes these numbers four ways (hand formulas, the `Tensor` engine, the scalar `Value` engine, and central differences) and asserts they agree to $10^{-7}$.

### Example 4.5: Gauss–Newton on an exponential decay fit (a mechanics example)

The peak amplitudes of a damped oscillation (or the temperature excess of a cooling part, or the voltage on a discharging capacitor) follow $y = a e^{-bt}$. Measurements at $t = 0, 1, 2$: $y = 2.0, 1.2, 0.8$. Fit $(a, b)$ by least squares: $r_i(a, b) = ae^{-bt_i} - y_i$, $f = \tfrac12\|\mathbf{r}\|^2$.

*Jacobian.* $\partial r_i/\partial a = e^{-bt_i}$, $\partial r_i/\partial b = -at_ie^{-bt_i}$, so $J_r$ is $3 \times 2$ with rows $(e^{-bt_i},\; -at_ie^{-bt_i})$.

*Start at $(a_0, b_0) = (1.5, 0.5)$.*
$$\mathbf{r} = \begin{pmatrix}-0.5\\-0.2902\\-0.2482\end{pmatrix},\qquad J_r = \begin{pmatrix}1 & 0\\0.6065 & -0.9098\\0.3679 & -1.1036\end{pmatrix},$$
$$J_r^\top J_r = \begin{pmatrix}1.5032 & -0.9578\\-0.9578 & 2.0457\end{pmatrix},\qquad J_r^\top\mathbf{r} = \begin{pmatrix}-0.7673\\0.5379\end{pmatrix} = \nabla f.$$
Solve (3.25): $\mathbf{h} = -(J_r^\top J_r)^{-1}J_r^\top\mathbf{r} = (0.4887, -0.0341)^\top$, so $(a_1, b_1) = (1.9887, 0.4659)$. The residual norm drops from $f = 0.198$ to $0.00136$ in one step. Iterating, $\|\mathbf{h}\|$ goes $0.49,\; 8\times10^{-3},\; 2\times10^{-4},\; 5\times10^{-6},\; 8\times10^{-8}, \dots$ (faster than linear; near-quadratic because the residual is small), converging to $(a, b) = (1.9884, 0.4743)$ with $f^* = 1.215 \times 10^{-3}$. SciPy's `least_squares` agrees to $10^{-8}$. At the solution, the dropped Hessian term $\sum_ir_i\nabla^2r_i = \begin{pmatrix}0 & 0\\0 & -0.046\end{pmatrix}$ against $J_r^\top J_r = \begin{pmatrix}1.537 & -1.366\\-1.366 & 3.903\end{pmatrix}$: about 1%, which is why Gauss–Newton behaved like Newton here. (Its off-diagonal is exactly zero at the solution: $\sum_i r_i\,\partial^2 r_i/\partial a\partial b = -\sum_ir_it_ie^{-bt_i} = (J_r^\top\mathbf{r})_2/a = 0$. Exercise B7.) Figure `gauss_newton_fit.png`.

## 5. Why it matters and where it connects

- **ML-2 (optimization).** Gradient descent is Theorem 3.5 iterated; its convergence proofs bound the Hessian's eigenvalues ($\lambda_{\min}$-strong convexity, $\lambda_{\max}$-smoothness), which is Theorem 3.14 applied globally. Newton (3.23), Gauss–Newton (3.25) and Levenberg–Marquardt are the second-order methods there, and the natural gradient is Definition 3.4 with the Fisher metric as the inner product.
- **ML-7 (neural networks).** §3.7 is backpropagation. ML-7 adds depth, different losses (via (3.21)), initialization (keeping $\sigma'(\mathbf{z}_\ell)$ away from the saturation seen in Example 4.4), and the autodiff viewpoint of §3.6, which the project's `Value` and `Tensor` classes implement.
- **ML-3, ML-9.** Logistic and softmax regression use (3.21) directly. The attention mechanism of ML-9 differentiates a softmax of a matrix product, which is (3.20) and (3.24) composed.
- **ST-1 (inference).** The score is $\nabla_\theta\ell$, the Fisher information is $-\mathbb{E}[\nabla^2_\theta\ell]$, and the Cramér–Rao bound is a statement about that Hessian's inverse. §3.9 is the first MLE computed there.
- **ST-3 (regression).** (3.12) is the normal equations; the ridge version is exercise B2; the Hessian $2A^\top A$ is what makes least squares convex and its condition number (LA-9, LA-13) what makes it numerically delicate.
- **ML-6 (unsupervised).** The M-step of EM for Gaussian mixtures is (3.38)–(3.39) with responsibilities as weights; PCA's variance maximization is a Rayleigh-quotient gradient (LA-8).
- **CM-2, CM-6.** The Euler–Lagrange equations are "$\nabla = 0$" for a functional (F6), and the stiffness matrix of small oscillations in CM-6 is the Hessian of the potential energy at equilibrium; Theorem 3.14 is why stable equilibria have positive-definite stiffness.
- **LA-13.** The conditioning of $H$ governs how fast gradient descent converges (Example 4.1's factor $(\kappa - 1)/(\kappa + 1)$) and how precisely you can solve (3.23).

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills: a gradient and Hessian by hand, a Jacobian and its linear approximation, trace identities, a softmax Jacobian, and the finite-difference U-curve. B1–B8 are the module: B1 (product rule and $d(A^{-1})$ from the definition), B3 (derive and code-check $\nabla_A\|AX - B\|_F^2$), B5 (the forward/reverse flop count), B6 (three-layer backprop), B8 (the Gaussian MLE from the gradients) are essential. C1–C3: the Jacobian of matrix inversion via `vec`, forward-mode autodiff timed against reverse, and the concavity of $\log\det$ from its second differential.

## 7. Coding project

`project.py` builds a reverse-mode autodiff engine and uses it to check everything in these notes:

1. **Finite differences.** `num_grad`, `num_jac`, `rel_error` (3.27), `grad_check`, and the error-vs-$h$ plot of §3.8.
2. **The catalogue.** Every identity (3.9)–(3.21) on random inputs; Jacobi's formula and $d(A^{-1})$ as first-order approximations; the Hessian (3.22); the `vec` trick (3.40); Examples 4.2 and 4.3 (the latter also symbolically with SymPy).
3. **`Value`.** A scalar autodiff class with `+`, `*`, `**`, `exp`, `log`, `tanh`, `relu` and a `backward()` that topologically sorts the graph and runs one VJP per node from the output backwards.
4. **`Tensor`.** The same for arrays: `@` with the rule (3.24), broadcasting `+` and `*` with the adjoint (3.26), elementwise maps with the diagonal-Jacobian rule (3.19), `sum`, transpose. Enough to train a two-layer net.
5. **Example 4.4 four ways:** the hand formulas (3.28)–(3.33), `Tensor`, `Value`, and finite differences, all asserted equal, with every intermediate printed.
6. **Batch backprop** (3.34)–(3.37), asserted equal to the engine, then used to train a 32-unit $\tanh$ network on $y = \sin x_1 + x_2^2$ by plain gradient descent; loss curve and fitted surface saved.
7. **Example 4.1:** Newton in one step vs gradient descent in 36 (or 22) steps, with the linear rate measured.
8. **Example 4.5:** Gauss–Newton with the full trace, compared with `scipy.optimize.least_squares`, plus the exact Hessian vs $J_r^\top J_r$ at the solution.
9. **Gaussian MLE:** (3.38)–(3.39) checked numerically, then gradient ascent from $(\mathbf{0}, I)$ recovers $\bar{\mathbf{x}}$ and $S/n$ on a sample of 400.
10. **Forward vs reverse:** the two multiplication orders of §3.6 timed on a chain of eight $400 \times 400$ Jacobians.

**Build it yourself first:** before reading Part 3, write `Value.__add__`, `Value.__mul__`, `Value.tanh` and `Value.backward` (the topological sort is the only non-obvious part). Before reading Part 5, write `backprop_single` from (3.28)–(3.33) with a shape comment on every line, then run it against `num_grad`. If a relative error comes out around $0.5$, you have a transpose wrong; §8 says where.

## 8. Where you'll struggle

- **Layout conventions are the swamp.** Half the matrix-calculus material online is a fight about whether $\partial\mathbf{y}/\partial\mathbf{x}$ is $m \times n$ or $n \times m$, and formulas copied between sources with different conventions are silently wrong by a transpose. Do not enter the swamp. Commit to the differential method: write $df$ as a linear expression in $d\mathbf{x}$ or $dX$, massage it into $\mathbf{g}^\top d\mathbf{x}$ or $\operatorname{tr}(G^\top dX)$, read off the answer. The only conventions you then need are the two in the box (gradient is a column, $\nabla_X f$ has the shape of $X$), and they are forced, not chosen. When you look something up in the *Matrix Cookbook*, re-derive it with differentials before using it; the lookup tells you what to expect, the derivation tells you the transpose.
- **The universal first bug: $\boldsymbol{\delta}_1 = (W_2\boldsymbol{\delta}_2) \odot \sigma'$ instead of $(W_2^\top\boldsymbol{\delta}_2) \odot \sigma'$.** When $W_2$ is square the wrong version runs and produces a wrong number; when it is not square it crashes, which is the lucky case. The fix is to write shapes next to everything: $W_2$ is $m \times h$, $\boldsymbol{\delta}_2$ is $m \times 1$, and the only product that gives an $h \times 1$ is $W_2^\top\boldsymbol{\delta}_2$. The deeper reason is the transposed chain rule at the end of §3.2: gradients flow backwards through $J^\top$, forward values flow through $J$.
- **"Jacobians multiply" is in the opposite order from how you compute.** You evaluate $f_1$ first, then $f_2$, but the derivative is $J_2J_1$ with $J_1$ on the *right*. Students write $J_1J_2$ because that is the order they computed in, and with square matrices it again fails silently. Say it as "outermost on the left," and check with (3.4) on a scalar-valued outer function: the product must end up $1 \times n$, which forces the $1 \times d$ Jacobian of the loss to be leftmost.
- **Gradient (column) versus Jacobian row for scalar functions.** The numbers are the same and the shapes differ, and the difference is not pedantry: $\nabla f$ depends on the inner product (Definition 3.4), $J$ does not. If you ever change the metric (natural gradient, Riemannian optimization), the gradient changes and the row does not. Until then, the practical rule is: a gradient is something you *add to $\mathbf{x}$*, so it must have the shape of $\mathbf{x}$.
- **Finite differences fail at both ends.** Students pick $h = 10^{-12}$ "to be safe" and get pure rounding noise, or $h = 10^{-2}$ and get truncation error, then conclude the analytic gradient is wrong. Look at `fd_error_vs_h.png`: the usable window is roughly $10^{-7}$ to $10^{-4}$ for central differences, and the best error is about $10^{-11}$, not $10^{-16}$. Use (3.27), use $h \approx 6 \times 10^{-6}$ scaled by $|x_i|$, and expect $10^{-9}$ to $10^{-11}$ for a correct smooth gradient. For ReLU networks, check at points where no pre-activation is within $h$ of zero.
- **Symmetric-matrix arguments.** (3.17) and (3.39) are gradients over *all* entries. If your optimizer only stores the upper triangle, off-diagonal derivatives double (§3.9). Know which one you are computing before comparing with a finite-difference check that perturbs one entry at a time.
- **Exercise B5 is harder than it looks** because the counting has to be done for both orders *and* for forward mode, and the temptation is to stop after reverse mode. Do all three columns of the table; the "why backprop exists" claim is the comparison, not any one count.

## 9. Check yourself

1. State the definition of $D\mathbf{f}(\mathbf{x})$ with the little-$o$ term. What is its matrix?
2. Why is the gradient a column and the Jacobian of a scalar function a row? Use the word "functional."
3. Write the chain rule for $\mathbf{g} \circ \mathbf{f}$ as a matrix product. Which Jacobian is on the left, and at what point is it evaluated?
4. Derive $d(A^{-1})$ from $AA^{-1} = I$ in two lines.
5. State the identification rule (3.8). Why does it make layout conventions unnecessary?
6. Derive $\nabla_{\mathbf{x}}\|A\mathbf{x} - \mathbf{b}\|^2$ with differentials and say what "$\nabla = \mathbf{0}$" is called in LA-5.
7. What is $\nabla_A\log\det A$? What is $\nabla_A\operatorname{tr}(AB)$? Which one has a transpose you could forget?
8. Write the Jacobian of softmax. What is its null vector, and why?
9. Why is $\nabla_{\mathbf{z}}$ of softmax-cross-entropy equal to $\mathbf{s} - \mathbf{y}$ without ever forming the softmax Jacobian?
10. For a chain of $k$ layers of width $d$ ending in a scalar, what do reverse mode and forward mode cost, and why?
11. Write the backprop recursion $\boldsymbol{\delta}_{\ell-1}$ in terms of $\boldsymbol{\delta}_\ell$. Where is the transpose, and what are the shapes?
12. In the batch version with $X$ as $n \times d$, what is $\nabla_{W_1}L$ and what is $\nabla_{\mathbf{b}_1}L$?
13. Why is the best finite-difference step about $\epsilon^{1/3}$ for central differences? What goes wrong on each side?
14. What does Gauss–Newton drop from the Hessian of $\tfrac12\|\mathbf{r}\|^2$, and when is that justified?

## 10. Pair with

- **Parr & Howard, *The Matrix Calculus You Need for Deep Learning* (arXiv 1802.01528, free).** Written for exactly the reader who has calculus and NumPy and no matrix calculus; builds the Jacobian from scalar partials slowly, with the neuron example throughout. Better than these notes at the entry-by-entry picture and at the vector chain rule as a sum over paths. Uses numerator layout; translate with (3.2).
- **Petersen & Pedersen, *The Matrix Cookbook* (free).** A lookup table, not a textbook. Use §2 (derivatives) to find what the answer should be, then derive it with §3.3. Its gradient-of-a-matrix formulas use the same shape convention as ours.
- **Magnus & Neudecker, *Matrix Differential Calculus with Applications in Statistics and Econometrics*, ch. 5–9.** The source of the differential method, with full rigor on identification and on symmetric and patterned matrices (the §3.9 subtlety is their ch. 9 duplication matrix). Dense; read ch. 5 and 8 first.
- **Minka, *Old and New Matrix Algebra Useful for Statistics* (free notes).** The differential method applied to exactly the Gaussian computations of §3.9 in ten pages, including the symmetric-parametrization issue. Read after §3.9.
- **Griewank & Walther, *Evaluating Derivatives*, ch. 3–4**, for the cheap gradient principle and forward versus reverse mode done properly; and **Griewank, "Who invented the reverse mode of differentiation?"** (*Optimization Stories*, 2012) for the history in §1.
- **Karpathy's micrograd lecture (video, "The spelled-out intro to neural networks and backpropagation").** The `Value` class of the project is this, built live over two hours; watch it after writing your own.
