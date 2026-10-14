# F5: Probability II: Joint Distributions, the Multivariate Gaussian, Limit Theorems, Conditional Expectation, Markov Chains, Entropy and KL

**Time:** ~14 hours · **Prerequisites:** F4 (random variables, named distributions, expectation, MGFs), LA-8 (spectral theorem, positive semidefinite matrices) · **Unlocks:** ST-1, SM-4, ML-1, ML-6

> **Conventions in this module.** Random variables are capitals ($X, Y$), their values lowercase ($x, y$). A random *vector* is bold: $\mathbf{X} = (X_1, \dots, X_k)^\top$, a column. PMFs and PDFs are both written $p(x)$ or $f(x)$ with the variable named by the argument; $p_{X,Y}(x,y)$ is the joint, $p_X(x)$ the marginal, $p_{Y|X}(y\mid x)$ the conditional. $\mathbb{E}$, $\mathrm{Var}$, $\mathrm{Cov}$ as in F4. The covariance matrix is $\Sigma$; it is always symmetric PSD. $\mathcal{N}(\boldsymbol\mu, \Sigma)$ is the multivariate normal with mean vector $\boldsymbol\mu$ and covariance $\Sigma$; $\varphi$ and $\Phi$ are the standard normal density and CDF. **Markov chains use row vectors:** $\boldsymbol\pi_{t+1} = \boldsymbol\pi_t P$, with $P_{ij} = P(\text{next} = j \mid \text{now} = i)$ and rows of $P$ summing to 1. (LA-7 used columns and $P^\top$; see §3.12.) Logs in information theory are natural (nats) unless written $\log_2$ (bits). $\Gamma(a) = \int_0^\infty x^{a-1}e^{-x}\,dx$; the Gamma distribution is parametrized by shape $a$ and **rate** $\lambda$, so its mean is $a/\lambda$. The Jacobian in a change of variables is the matrix of partial derivatives of the *old* variables with respect to the *new* ones, and its determinant always appears in absolute value.

## 0. Where this sits

F4 gave you one random variable at a time: its distribution, its expectation, its MGF. Almost nothing you will do later involves one random variable. Regression (ST-3) is a statement about the joint distribution of $(X, Y)$. The Kalman filter (ST-8) is a statement about the conditional distribution of a Gaussian vector. Statistical mechanics (SM-4) is the maximum-entropy distribution over microstates. Every loss function in ML-3 is a cross-entropy. This module is where "several random variables at once" gets built, and the linear algebra you rebuilt in LA-1 through LA-8 is what makes it tractable: a covariance matrix is a PSD matrix, a Gaussian's level sets are the ellipsoids of a quadratic form, a Markov chain is a matrix acting on a vector, and conditional expectation is an orthogonal projection.

You need F4 fluently (LOTUS, MGFs, the named distributions, Bayes' rule) and LA-8 (spectral theorem, $\mathbf{a}^\top \Sigma \mathbf{a} \ge 0$, Cholesky from LA-10 helps). This module is the direct prerequisite of ST-1 (inference), SM-4 (canonical ensemble as max-entropy), ML-1 (the learning problem) and ML-6 (Gaussian mixtures, PCA).

## 1. The problem that forced this idea

**Galton and the heights (1886).** Francis Galton collected the heights of 928 adult children and their parents (he averaged the two parents into a "mid-parent" height). He wanted to know how heredity worked numerically, and the question cannot be asked about one variable. He needed the *joint* distribution of (parent height, child height). What he found surprised him: tall parents had tall children, but less tall than themselves; short parents had short children, but less short. He called it "regression towards mediocrity" and at first thought it was a biological force pulling the population toward the mean. It is not. It is a property of any pair of imperfectly correlated variables, and the direction is symmetric: tall children also have, on average, less tall parents. Explaining this correctly requires the conditional expectation $\mathbb{E}[Y \mid X]$ and the correlation coefficient, both of which Galton (and then Karl Pearson) invented for this purpose. Example 4.5 redoes his computation.

**Gauss and the errors (1809).** Gauss needed the distribution of the error in a *combination* of many measurements. If each measurement has an error, the error of the least-squares estimate is a linear combination of them. He argued (not rigorously) that such errors should be normal, and that a vector of jointly normal errors is described by a matrix. Laplace proved the first version of the central limit theorem in 1810, which is the actual justification: sums of many small independent errors are approximately normal whatever the individual errors look like. The multivariate normal, with its covariance matrix, is the object that makes "a vector of correlated normal errors" precise.

**Markov and Nekrasov (1906).** Pavel Nekrasov, a Moscow mathematician with theological commitments, claimed that the law of large numbers *requires* independence, and drew from this the conclusion that observed statistical regularities in human affairs proved free will. Andrei Markov, who disliked Nekrasov on every available axis, constructed sequences of *dependent* random variables, each depending on the one before, for which the law of large numbers still holds. These are Markov chains. The dependence structure is the simplest possible beyond independence, and it turned out to describe diffusion, queues, genetics, Google's ranking of web pages, and the Monte Carlo methods that make Bayesian statistics computable (ST-6).

**Shannon and Boltzmann (1948, 1877).** Claude Shannon, designing communication systems at Bell Labs, needed a number that measured how much information a source produces per symbol. He wrote down a short list of properties such a measure should have (continuity, monotone in the number of equally likely outcomes, additive over independent choices) and proved that the only function with these properties is $-\sum_i p_i \log p_i$. He then noticed that this was the formula Boltzmann had written in 1877 for the entropy of a gas, where the $p_i$ are probabilities of microstates. The coincidence is not a coincidence: both are measuring the same thing, how spread out a probability distribution is, and SM-4 will derive the Boltzmann distribution from Shannon's entropy. The same quantity, as *relative* entropy or Kullback–Leibler divergence, is the loss function behind most of modern machine learning.

Four problems, one theme: **one random variable is never enough.** You need the joint distribution, and you need linear algebra to handle it.

## 2. Intuition

### A joint distribution is a landscape

Picture the joint density $f(x, y)$ of two continuous random variables as a landscape over the $xy$-plane: a hill whose height at $(x, y)$ says how likely that pair is. The total volume under the hill is 1.

```
   y
   ^      . . : : : . .
   |    . : : # # # : : .          height = f(x, y)
   |   . : # # # # # # : .         volume under the surface = 1
   |    . : : # # # : : .
   |      . . : : : . .
   +---------------------------> x
```

Three operations on the landscape account for half the module:

- **Marginal = shadow.** The marginal density $f_X(x)$ is what you get by collapsing the hill onto the $x$-axis: for each $x$, add up (integrate) all the height along the vertical line through $x$. It is the shadow the hill casts on the $x$-axis when lit from the $y$-direction. You lose all information about how $y$ relates to $x$. Two very different hills can cast the same two shadows.
- **Conditional = slice, renormalized.** Fix $Y = y_0$ and cut the hill with the vertical plane $y = y_0$. The cross-section is a curve in $x$; it does not have area 1, so divide by its area, which is $f_Y(y_0)$. That renormalized slice is $f_{X|Y}(x \mid y_0)$.
- **Independence = the hill is a product.** If the landscape factors as $f(x, y) = g(x)h(y)$, every slice at fixed $y_0$ has the *same shape* (only the scale $h(y_0)$ changes), so knowing $y$ tells you nothing about $x$. Any hill that is not separable in this way encodes dependence.

### Covariance: do they lean together?

Covariance asks: when $X$ is above its mean, is $Y$ usually above its mean too? Positive if the hill is stretched along the diagonal, negative along the anti-diagonal, zero if it has no diagonal lean. Note the word *lean*: a hill shaped like a ring, or like a U, can have zero covariance while $X$ and $Y$ are completely dependent. Covariance detects linear leaning and nothing else.

### The Gaussian hill is an ellipse

The multivariate normal is the hill whose level curves are ellipses, all with the same center and the same axes. You already know what determines those axes. The exponent is $-\tfrac12 (\mathbf{x} - \boldsymbol\mu)^\top \Sigma^{-1} (\mathbf{x} - \boldsymbol\mu)$, a quadratic form in LA-8's sense, and LA-8 Theorem (spectral theorem) says its level sets are ellipsoids whose axes are the eigenvectors of $\Sigma$. A long axis means a direction of large variance. The guess you should make now: the semi-axis lengths are proportional to $\sqrt{\lambda_i}$, since variance is the square of a length.

Slice this hill at $x_2 = c$ and you get a 1D Gaussian whose center has moved (linearly in $c$) and whose width has shrunk. The moved center *is* the regression line. The shrunken width *is* the residual variance. That one slice contains linear regression, the Kalman filter, and Gaussian process regression.

### A Markov chain is a matrix acting on a vector

Weather tomorrow depends on weather today and nothing earlier. Write today's probabilities as a row vector $\boldsymbol\pi_t = (P(\text{sun}), P(\text{cloud}), P(\text{rain}))$ and the transition probabilities as a matrix $P$. Then tomorrow is $\boldsymbol\pi_t P$, the day after is $\boldsymbol\pi_t P^2$, and the long-run weather is whatever $\boldsymbol\pi_t P^t$ converges to. From LA-7 you know what $P^t$ does: it is dominated by the eigenvalue of largest modulus. For a stochastic matrix that eigenvalue is 1, and its eigenvector is the long-run climate. The second eigenvalue tells you how fast you forget today's weather.

### Entropy is average surprise

If an event has probability $p$, its surprise is $-\log p$: certain events surprise you not at all, rare ones a lot, and the surprises of independent events add. The entropy of a distribution is the expected surprise. A fair coin has entropy 1 bit; a coin that always lands heads has entropy 0; a six-sided die has $\log_2 6 \approx 2.58$ bits. The guess to make: the uniform distribution has the largest entropy for a given number of outcomes, because it maximizes the average surprise. KL divergence will be the *extra* surprise you suffer when you believe $q$ but the world runs on $p$. It can never be negative, and that single fact is behind maximum likelihood, maximum entropy, and the ELBO (ML-10).

## 3. Formal development

### 3.1 Joint, marginal, conditional

**Definition 3.1 (Joint distribution).** For discrete $X, Y$, the **joint PMF** is $p_{X,Y}(x, y) = P(X = x, Y = y)$. For continuous $X, Y$, the **joint PDF** is the function $f_{X,Y} \ge 0$ with
$$P\big((X, Y) \in A\big) = \iint_A f_{X,Y}(x, y)\,dx\,dy \quad \text{for every region } A \subseteq \mathbb{R}^2,$$
so the joint CDF is $F_{X,Y}(x, y) = P(X \le x, Y \le y) = \int_{-\infty}^x\int_{-\infty}^y f_{X,Y}(s, t)\,dt\,ds$ and $f_{X,Y} = \partial^2 F_{X,Y}/\partial x\,\partial y$ wherever that derivative exists. Both sum/integrate to 1 over the plane.

▶ **Try it first.** Write the event $\{X = x\}$ as a disjoint union of events involving both $X$ and $Y$. Then use additivity.

**Derivation (marginals).**

1. Discrete case. The events $\{X = x, Y = y\}$ for different $y$ are disjoint, and their union over all $y$ is $\{X = x\}$. By additivity of $P$,
$$p_X(x) = \sum_y p_{X,Y}(x, y). \tag{3.1}$$
2. Continuous case. $F_X(x) = P(X \le x) = P\big((X, Y) \in (-\infty, x] \times \mathbb{R}\big) = \int_{-\infty}^x \left[\int_{-\infty}^\infty f_{X,Y}(s, y)\,dy\right] ds$. Differentiate in $x$:
$$f_X(x) = \int_{-\infty}^\infty f_{X,Y}(x, y)\,dy. \tag{3.2}$$ ∎

■ **What just happened.** To get the distribution of $X$ alone, sum or integrate $Y$ *out*. This is the shadow of §2. The joint determines the marginals; the marginals do not determine the joint.

**Definition 3.2 (Conditional distribution).** For discrete $X, Y$ and $p_X(x) > 0$, the **conditional PMF** of $Y$ given $X = x$ is
$$p_{Y|X}(y \mid x) = \frac{p_{X,Y}(x, y)}{p_X(x)}, \tag{3.3}$$
which is the definition of conditional probability from F4 applied to the events $\{Y = y\}$ and $\{X = x\}$. For continuous $X, Y$ and $f_X(x) > 0$, the **conditional PDF** is defined by the same formula,
$$f_{Y|X}(y \mid x) = \frac{f_{X,Y}(x, y)}{f_X(x)}. \tag{3.4}$$

Why (3.4) is the right definition when $P(X = x) = 0$: condition on the event $X \in [x, x + \delta]$, which has probability $\approx f_X(x)\delta$. Then $P(Y \in [y, y+\epsilon] \mid X \in [x, x+\delta]) \approx f_{X,Y}(x,y)\,\delta\epsilon / (f_X(x)\delta) = \frac{f_{X,Y}(x,y)}{f_X(x)}\epsilon$. Let $\delta \to 0$. The $\delta$ cancels and what is left is a density in $y$. For each fixed $x$, (3.4) is a genuine density in $y$: it is nonnegative and $\int f_{Y|X}(y\mid x)\,dy = f_X(x)/f_X(x) = 1$ by (3.2). This is the renormalized slice.

Two immediate consequences, both read off (3.3)–(3.4):

- **Chain rule.** $f_{X,Y}(x, y) = f_X(x)\, f_{Y|X}(y \mid x) = f_Y(y)\, f_{X|Y}(x \mid y)$.
- **Bayes' rule for densities.** $f_{X|Y}(x \mid y) = \dfrac{f_{Y|X}(y \mid x)\, f_X(x)}{f_Y(y)}$, with $f_Y(y) = \int f_{Y|X}(y \mid x) f_X(x)\,dx$ (the law of total probability, continuous form). This is the whole of Bayesian inference (ST-5) in one line: $x$ is the parameter, $y$ the data.

Everything extends to $k$ variables: $f(x_1, \dots, x_k)$, marginals by integrating out whichever variables you don't want, and the chain rule $f(x_1, \dots, x_k) = f(x_1) f(x_2 \mid x_1) f(x_3 \mid x_1, x_2) \cdots$.

### 3.2 Independence as factorization

**Definition 3.3.** $X$ and $Y$ are **independent** (written $X \perp Y$) if $P(X \in A, Y \in B) = P(X \in A)\,P(Y \in B)$ for all sets $A, B$.

**Theorem 3.4 (Three characterizations of independence).** Let $(X, Y)$ have joint PMF or PDF $p(x, y)$ with marginals $p_X, p_Y$. The following are equivalent:

(a) $p(x, y) = p_X(x)\,p_Y(y)$ for all $x, y$.
(b) $p_{Y|X}(y \mid x) = p_Y(y)$ for all $y$ and all $x$ with $p_X(x) > 0$.
(c) $p(x, y) = g(x)\,h(y)$ for some functions $g, h$ (not necessarily densities).

▶ **Try it first.** (a)⇔(b) is one line each way from (3.3). For (c)⇒(a), sum (3.1) and its mirror image, and figure out what $\sum_x g(x)$ times $\sum_y h(y)$ must be.

**Derivation.** We write sums; replace by integrals for the continuous case.

1. (a)⇒(b): $p_{Y|X}(y\mid x) = p(x,y)/p_X(x) = p_X(x)p_Y(y)/p_X(x) = p_Y(y)$.
2. (b)⇒(a): $p(x, y) = p_X(x)\,p_{Y|X}(y \mid x) = p_X(x)\,p_Y(y)$ when $p_X(x) > 0$; when $p_X(x) = 0$, both sides are $0$ because $0 \le p(x,y) \le p_X(x)$.
3. (a)⇒(c): take $g = p_X$, $h = p_Y$.
4. (c)⇒(a): Let $G = \sum_x g(x)$ and $H = \sum_y h(y)$. By (3.1), $p_X(x) = \sum_y g(x)h(y) = g(x)H$, and likewise $p_Y(y) = h(y)G$. Summing the joint over everything, $1 = \sum_{x,y} g(x)h(y) = GH$. Therefore $p_X(x)p_Y(y) = g(x)h(y)\,GH = g(x)h(y) = p(x,y)$.
5. Finally (a) is equivalent to Definition 3.3: summing (a) over $x \in A$, $y \in B$ gives the product rule for events; conversely taking $A = \{x\}$, $B = \{y\}$ (or shrinking intervals in the continuous case) recovers (a). ∎

■ **What just happened.** Characterization (c) is the one you will use most. To check independence you do not need to compute the marginals: if the joint density *factors* into something depending only on $x$ times something depending only on $y$, done. This is how we will prove that Box–Muller produces independent normals and that uncorrelated jointly-Gaussian variables are independent. One caution: the *domain* must factor too. A density that is $2$ on the triangle $0 < x < y < 1$ looks like "$g(x)h(y)$ with both constant," but the indicator of the triangle does not factor, and $X, Y$ are dependent.

**Corollary 3.5.** If $X \perp Y$ then $\mathbb{E}[g(X)h(Y)] = \mathbb{E}[g(X)]\,\mathbb{E}[h(Y)]$ for any functions $g, h$ for which the expectations exist.

*Proof.* By LOTUS in two variables (F4), $\mathbb{E}[g(X)h(Y)] = \sum_{x,y} g(x)h(y)p_X(x)p_Y(y) = \left(\sum_x g(x)p_X(x)\right)\left(\sum_y h(y)p_Y(y)\right)$. ∎

### 3.3 Change of variables in two dimensions

From F4 you know the one-dimensional rule: if $Y = g(X)$ with $g$ strictly monotone, $f_Y(y) = f_X(g^{-1}(y))\,\big|\tfrac{d}{dy}g^{-1}(y)\big|$. The absolute value of the derivative is the factor by which $g^{-1}$ stretches a small interval. In two dimensions, "stretch factor of a small region" is a determinant, by LA-6.

**Theorem 3.6 (Change of variables).** Let $(X, Y)$ have joint density $f_{X,Y}$, and let $T : (x, y) \mapsto (u, v)$ be a bijection between open regions of the plane, continuously differentiable with continuously differentiable inverse. Put $(U, V) = T(X, Y)$. Then
$$f_{U,V}(u, v) = f_{X,Y}\big(T^{-1}(u, v)\big)\,\left|\det \frac{\partial(x, y)}{\partial(u, v)}\right|, \tag{3.5}$$
where $\dfrac{\partial(x,y)}{\partial(u,v)} = \begin{pmatrix} \partial x/\partial u & \partial x/\partial v \\ \partial y/\partial u & \partial y/\partial v\end{pmatrix}$ is the Jacobian matrix of the *inverse* map (old variables as functions of the new ones), evaluated at $(u, v)$.

▶ **Try it first.** Start from $P((U,V) \in B) = P((X, Y) \in T^{-1}(B))$, write the right side as an integral, and substitute. What is $dx\,dy$ in terms of $du\,dv$?

**Derivation.**

1. For any region $B$ in the $(u, v)$-plane, $P((U, V) \in B) = P\big((X, Y) \in T^{-1}(B)\big) = \iint_{T^{-1}(B)} f_{X,Y}(x, y)\,dx\,dy$.
2. Substitute $(x, y) = T^{-1}(u, v)$. By the multivariable substitution rule from calculus, the area element transforms as $dx\,dy = \left|\det \tfrac{\partial(x,y)}{\partial(u,v)}\right| du\,dv$. (LA-6: the linear approximation of $T^{-1}$ near a point is multiplication by the Jacobian matrix, and $|\det|$ of a matrix is the factor by which it scales area. The sign of the determinant records orientation, which probability does not care about: hence the absolute value.) The region $T^{-1}(B)$ becomes $B$:
$$P((U, V) \in B) = \iint_B f_{X,Y}\big(T^{-1}(u, v)\big)\left|\det \frac{\partial(x,y)}{\partial(u,v)}\right| du\,dv.$$
3. This holds for every $B$. A function whose integral over every region equals $P((U,V) \in B)$ is, by Definition 3.1, the joint density of $(U, V)$. ∎

■ **What just happened.** Same rule as 1D, with $|dx/dy|$ replaced by $|\det|$ of the Jacobian. Two things people get wrong (see §8): it is the Jacobian of the map from *new* to *old*, and it carries an absolute value. If it is easier to differentiate $T$ than $T^{-1}$, use the inverse-function fact $\det \tfrac{\partial(x,y)}{\partial(u,v)} = 1 \big/ \det \tfrac{\partial(u,v)}{\partial(x,y)}$ (LA-2/LA-6: the Jacobian matrices of $T$ and $T^{-1}$ are inverses, so their determinants are reciprocals).

**Example (Box–Muller).** Let $U_1, U_2$ be independent $\mathrm{Uniform}(0, 1)$. Define
$$R = \sqrt{-2 \ln U_1}, \qquad \Theta = 2\pi U_2, \qquad Z_1 = R\cos\Theta, \qquad Z_2 = R\sin\Theta. \tag{3.6}$$
**Claim:** $Z_1, Z_2$ are independent $\mathcal{N}(0, 1)$.

▶ **Try it first.** First find the distribution of $S = R^2$ (one line from F4). Then change variables from $(S, \Theta)$ to $(Z_1, Z_2)$; the inverse map is $s = z_1^2 + z_2^2$, $\theta = \operatorname{atan2}(z_2, z_1)$, and it is easier to differentiate than the forward map.

**Derivation.**

1. $S = -2\ln U_1$. For $s > 0$, $P(S > s) = P(U_1 < e^{-s/2}) = e^{-s/2}$. So $S \sim \mathrm{Exponential}(\tfrac12)$ with density $f_S(s) = \tfrac12 e^{-s/2}$ on $s > 0$. And $\Theta \sim \mathrm{Uniform}(0, 2\pi)$ with density $\tfrac{1}{2\pi}$. They are independent (functions of independent $U_1, U_2$), so $f_{S,\Theta}(s, \theta) = \tfrac{1}{4\pi} e^{-s/2}$ on $(0, \infty) \times (0, 2\pi)$.
2. The map $(s, \theta) \mapsto (z_1, z_2) = (\sqrt{s}\cos\theta, \sqrt{s}\sin\theta)$ is a bijection from $(0,\infty)\times(0,2\pi)$ onto the plane minus a ray (a set of probability zero, which we ignore). Its inverse is $s = z_1^2 + z_2^2$, $\theta = \operatorname{atan2}(z_2, z_1)$.
3. Jacobian of the inverse (old $(s,\theta)$ with respect to new $(z_1, z_2)$):
$$\frac{\partial(s, \theta)}{\partial(z_1, z_2)} = \begin{pmatrix} 2z_1 & 2z_2 \\ -\dfrac{z_2}{z_1^2 + z_2^2} & \dfrac{z_1}{z_1^2 + z_2^2}\end{pmatrix}, \qquad \det = \frac{2z_1^2 + 2z_2^2}{z_1^2 + z_2^2} = 2.$$
(Using $\partial_{z_1}\operatorname{atan2}(z_2, z_1) = -z_2/(z_1^2+z_2^2)$ and $\partial_{z_2} = z_1/(z_1^2+z_2^2)$.)
4. By (3.5), with the roles "old = $(s,\theta)$, new = $(z_1, z_2)$":
$$f_{Z_1, Z_2}(z_1, z_2) = \frac{1}{4\pi} e^{-(z_1^2 + z_2^2)/2}\cdot 2 = \frac{1}{\sqrt{2\pi}}e^{-z_1^2/2}\cdot\frac{1}{\sqrt{2\pi}}e^{-z_2^2/2} = \varphi(z_1)\varphi(z_2).$$
5. The joint density factors (Theorem 3.4(c)) into two standard normal densities. So $Z_1 \perp Z_2$ and each is $\mathcal{N}(0,1)$. ∎

■ **What just happened.** Two uniforms in, two independent normals out, exactly. Reading it backwards is also useful: a pair of independent standard normals has a squared radius that is $\mathrm{Exponential}(\tfrac12)$ (which is $\chi^2_2$, §3.5) and an angle that is uniform and independent of the radius. That is why the 2D standard normal is rotationally symmetric. Every normal random number your computer has ever produced came from (3.6) or a close relative. The project implements it and compares to `scipy.stats.norm` with a QQ plot.

### 3.4 Sums of independent random variables: convolution

**Theorem 3.7 (Convolution).** Let $X \perp Y$ and $S = X + Y$. Then
$$p_S(s) = \sum_x p_X(x)\,p_Y(s - x) \quad \text{(discrete)}, \qquad f_S(s) = \int_{-\infty}^\infty f_X(x)\,f_Y(s - x)\,dx \quad \text{(continuous)}. \tag{3.7}$$

▶ **Try it first.** Discrete: condition on the value of $X$ and use the law of total probability. Continuous: compute the CDF $P(X + Y \le s)$ by conditioning on $X = x$, then differentiate.

**Derivation.**

1. Discrete. $P(S = s) = \sum_x P(X = x)\,P(X + Y = s \mid X = x) = \sum_x P(X = x)\,P(Y = s - x \mid X = x) = \sum_x p_X(x)\,p_Y(s - x)$, the last step by independence.
2. Continuous. $F_S(s) = P(X + Y \le s) = \int f_X(x)\,P(Y \le s - x \mid X = x)\,dx = \int f_X(x)\,F_Y(s - x)\,dx$. Differentiate under the integral: $f_S(s) = \int f_X(x) f_Y(s - x)\,dx$. ∎

■ **What just happened.** The density of a sum is the convolution of the densities. In F3 you met convolution as the operation that multiplication of Fourier transforms turns into. The same is true here: the MGF (or characteristic function) of a sum of independent variables is the *product* of the MGFs, $M_{X+Y}(t) = \mathbb{E}[e^{tX}e^{tY}] = M_X(t)M_Y(t)$ by Corollary 3.5. When the MGFs are nice, that is the faster route.

Three sums you must know:

**Normal + Normal.** If $X \sim \mathcal{N}(\mu_1, \sigma_1^2) \perp Y \sim \mathcal{N}(\mu_2, \sigma_2^2)$, then $X + Y \sim \mathcal{N}(\mu_1 + \mu_2, \sigma_1^2 + \sigma_2^2)$. *Proof.* From F4, $M_X(t) = \exp(\mu_1 t + \tfrac12\sigma_1^2 t^2)$. So $M_{X+Y}(t) = \exp\big((\mu_1+\mu_2)t + \tfrac12(\sigma_1^2+\sigma_2^2)t^2\big)$, which is the MGF of $\mathcal{N}(\mu_1+\mu_2, \sigma_1^2+\sigma_2^2)$. MGFs determine distributions (F4, cited there). ∎

**Poisson + Poisson.** If $X \sim \mathrm{Pois}(\lambda) \perp Y \sim \mathrm{Pois}(\mu)$ then $X + Y \sim \mathrm{Pois}(\lambda + \mu)$. *Proof by (3.7).*
$$p_S(s) = \sum_{k=0}^{s} \frac{e^{-\lambda}\lambda^k}{k!}\cdot\frac{e^{-\mu}\mu^{s-k}}{(s-k)!} = \frac{e^{-(\lambda+\mu)}}{s!}\sum_{k=0}^s \binom{s}{k}\lambda^k\mu^{s-k} = \frac{e^{-(\lambda+\mu)}(\lambda+\mu)^s}{s!}$$
by the binomial theorem. ∎

**Exponential + ⋯ + Exponential = Gamma.** Define the **Gamma distribution** $\mathrm{Gamma}(a, \lambda)$ (shape $a > 0$, rate $\lambda > 0$) by the density
$$f(x) = \frac{\lambda^a}{\Gamma(a)}x^{a-1}e^{-\lambda x}, \qquad x > 0. \tag{3.8}$$
It integrates to 1 because substituting $u = \lambda x$ gives $\int_0^\infty x^{a-1}e^{-\lambda x}dx = \lambda^{-a}\int_0^\infty u^{a-1}e^{-u}du = \Gamma(a)/\lambda^a$. Note $\mathrm{Gamma}(1, \lambda) = \mathrm{Exponential}(\lambda)$.

**Theorem 3.8.** If $X_1, \dots, X_n$ are iid $\mathrm{Exponential}(\lambda)$, then $S_n = X_1 + \cdots + X_n \sim \mathrm{Gamma}(n, \lambda)$.

▶ **Try it first.** Induction on $n$ using (3.7). The integrand simplifies to a power of $x$ because the exponentials combine.

**Derivation.**

1. Base case $n = 1$: $\mathrm{Gamma}(1,\lambda)$ is the exponential density, since $\Gamma(1) = 1$.
2. Inductive step. Assume $S_{n-1} \sim \mathrm{Gamma}(n-1, \lambda)$, with $\Gamma(n-1) = (n-2)!$. Then $S_n = S_{n-1} + X_n$ with $S_{n-1} \perp X_n$, and for $s > 0$,
$$f_{S_n}(s) = \int_0^s \frac{\lambda^{n-1}x^{n-2}e^{-\lambda x}}{(n-2)!}\cdot\lambda e^{-\lambda(s-x)}\,dx = \frac{\lambda^n e^{-\lambda s}}{(n-2)!}\int_0^s x^{n-2}\,dx = \frac{\lambda^n e^{-\lambda s}}{(n-2)!}\cdot\frac{s^{n-1}}{n-1} = \frac{\lambda^n s^{n-1}e^{-\lambda s}}{(n-1)!}.$$
(The limits are $0$ to $s$ because both densities vanish for negative arguments.) This is (3.8) with $a = n$. ∎

■ **What just happened.** In a Poisson process of rate $\lambda$ (F4), the gaps between arrivals are iid $\mathrm{Exponential}(\lambda)$. So the time of the $n$-th arrival is $\mathrm{Gamma}(n, \lambda)$. That is the "Poisson process link": counting arrivals in a fixed time gives Poisson; waiting for a fixed number of arrivals gives Gamma. The two are the same process looked at from two directions.

### 3.5 The Gamma distribution, $\chi^2$, and Gamma–Poisson

**MGF of the Gamma.** For $t < \lambda$,
$$M(t) = \int_0^\infty e^{tx}\frac{\lambda^a}{\Gamma(a)}x^{a-1}e^{-\lambda x}dx = \frac{\lambda^a}{\Gamma(a)}\int_0^\infty x^{a-1}e^{-(\lambda - t)x}dx = \frac{\lambda^a}{\Gamma(a)}\cdot\frac{\Gamma(a)}{(\lambda-t)^a} = \left(\frac{\lambda}{\lambda - t}\right)^a. \tag{3.9}$$
Consequences: $M'(0) = a/\lambda$ is the mean and $M''(0) - M'(0)^2 = a/\lambda^2$ is the variance (differentiate (3.9) twice; F4). And if $X \sim \mathrm{Gamma}(a, \lambda) \perp Y \sim \mathrm{Gamma}(b, \lambda)$ with the *same rate*, then $M_{X+Y}(t) = (\lambda/(\lambda-t))^{a+b}$, so $X + Y \sim \mathrm{Gamma}(a+b, \lambda)$: shapes add. Theorem 3.8 is the case $a = b = 1$ repeated.

**$\chi^2$ is a Gamma.** Let $Z \sim \mathcal{N}(0,1)$ and $W = Z^2$. For $w > 0$, $P(W \le w) = P(-\sqrt w \le Z \le \sqrt w) = 2\Phi(\sqrt w) - 1$. Differentiating, $f_W(w) = 2\varphi(\sqrt w)\cdot\frac{1}{2\sqrt w} = \frac{1}{\sqrt{2\pi}}w^{-1/2}e^{-w/2}$. Compare with (3.8) at $a = \tfrac12$, $\lambda = \tfrac12$: $\frac{(1/2)^{1/2}}{\Gamma(1/2)}w^{-1/2}e^{-w/2} = \frac{1}{\sqrt2\sqrt\pi}w^{-1/2}e^{-w/2}$, using $\Gamma(\tfrac12) = \sqrt\pi$. Same. So $Z^2 \sim \mathrm{Gamma}(\tfrac12, \tfrac12) =: \chi^2_1$, and by "shapes add," $Z_1^2 + \cdots + Z_n^2 \sim \mathrm{Gamma}(\tfrac n2, \tfrac12) =: \chi^2_n$, with mean $n$ and variance $2n$. In particular $\chi^2_2 = \mathrm{Gamma}(1, \tfrac12) = \mathrm{Exponential}(\tfrac12)$, which is the squared radius in Box–Muller.

**Gamma–Poisson conjugacy.** If the rate $\Lambda \sim \mathrm{Gamma}(a, b)$ and, given $\Lambda = \lambda$, $X \sim \mathrm{Pois}(\lambda)$, then by Bayes' rule for densities the posterior is
$$f(\lambda \mid X = x) \propto \lambda^{a-1}e^{-b\lambda}\cdot e^{-\lambda}\lambda^x \propto \lambda^{a + x - 1}e^{-(b+1)\lambda},$$
so $\Lambda \mid X = x \sim \mathrm{Gamma}(a + x, b + 1)$. The prior shape is "pseudo-counts" and the rate is "pseudo-exposure." ST-5 generalizes this to exponential families.

### 3.6 Order statistics and the Beta distribution

Let $X_1, \dots, X_n$ be iid with CDF $F$ and density $f$. Sort them: $X_{(1)} \le X_{(2)} \le \cdots \le X_{(n)}$. These are the **order statistics**; $X_{(n)} = \max$, $X_{(1)} = \min$.

▶ **Try it first.** The max is $\le m$ exactly when *every* $X_i$ is $\le m$. Write that as a product, then differentiate.

**Derivation (max and min).**

1. $P(X_{(n)} \le m) = P(X_1 \le m, \dots, X_n \le m) = F(m)^n$ by independence. So $f_{X_{(n)}}(m) = nF(m)^{n-1}f(m)$.
2. $P(X_{(1)} > m) = P(\text{all } X_i > m) = (1 - F(m))^n$. So $f_{X_{(1)}}(m) = n(1 - F(m))^{n-1}f(m)$. ∎

**Theorem 3.9 (Density of the $k$-th order statistic).**
$$f_{X_{(k)}}(x) = \frac{n!}{(k-1)!\,(n-k)!}\,F(x)^{k-1}\big(1 - F(x)\big)^{n-k}f(x). \tag{3.10}$$

*Derivation.* $f_{X_{(k)}}(x)\,dx \approx P(X_{(k)} \in [x, x+dx])$. For this to happen, exactly one of the $n$ variables must land in $[x, x+dx]$ (probability $f(x)dx$ each), exactly $k-1$ must be below $x$ (probability $F(x)$ each), and the remaining $n-k$ above (probability $1 - F(x)$ each). The number of ways to assign the $n$ variables to the three groups is the multinomial coefficient $\frac{n!}{(k-1)!\,1!\,(n-k)!}$. Multiply. (Two variables landing in $[x, x+dx]$ has probability of order $dx^2$ and drops out.) ∎

For $X_i \sim \mathrm{Uniform}(0,1)$: $F(x) = x$, $f(x) = 1$, and (3.10) becomes $\propto x^{k-1}(1-x)^{n-k}$ on $[0,1]$. That shape gets a name.

**Definition 3.10 (Beta distribution).** $\mathrm{Beta}(a, b)$ for $a, b > 0$ has density
$$f(x) = \frac{1}{B(a,b)}x^{a-1}(1-x)^{b-1}, \quad 0 < x < 1, \qquad B(a,b) = \int_0^1 x^{a-1}(1-x)^{b-1}dx. \tag{3.11}$$
$\mathrm{Beta}(1,1)$ is $\mathrm{Uniform}(0,1)$. Theorem 3.9 says: **the $k$-th smallest of $n$ iid Uniforms is $\mathrm{Beta}(k, n-k+1)$**, and comparing (3.10) with (3.11) also gives $B(k, n-k+1) = \frac{(k-1)!(n-k)!}{n!}$ for integers. The general formula is next.

**Theorem 3.11 (The Beta function).** $B(a, b) = \dfrac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}$.

▶ **Try it first.** Take $X \sim \mathrm{Gamma}(a, 1) \perp Y \sim \mathrm{Gamma}(b, 1)$ and change variables to the total $T = X + Y$ and the fraction $W = X/(X+Y)$. Use Theorem 3.6. You know the distribution of $T$ already from §3.5.

**Derivation.**

1. The inverse map is $x = tw$, $y = t(1 - w)$, from $(0,\infty)\times(0,1)$ onto the positive quadrant. Jacobian of the inverse:
$$\frac{\partial(x,y)}{\partial(t,w)} = \begin{pmatrix} w & t \\ 1 - w & -t\end{pmatrix}, \qquad |\det| = |{-tw} - t(1-w)| = t.$$
2. By (3.5), for $t > 0$, $0 < w < 1$:
$$f_{T,W}(t, w) = \frac{(tw)^{a-1}e^{-tw}}{\Gamma(a)}\cdot\frac{(t(1-w))^{b-1}e^{-t(1-w)}}{\Gamma(b)}\cdot t = \underbrace{\frac{t^{a+b-1}e^{-t}}{\Gamma(a)\Gamma(b)}}_{\text{depends on } t}\;\underbrace{w^{a-1}(1-w)^{b-1}}_{\text{depends on } w}.$$
3. The density factors, so $T \perp W$ (Theorem 3.4(c)). Integrating out $w$: $f_T(t) = \frac{t^{a+b-1}e^{-t}}{\Gamma(a)\Gamma(b)}B(a,b)$.
4. But from §3.5, $T \sim \mathrm{Gamma}(a+b, 1)$, so $f_T(t) = \frac{t^{a+b-1}e^{-t}}{\Gamma(a+b)}$. Matching the two expressions: $\frac{B(a,b)}{\Gamma(a)\Gamma(b)} = \frac{1}{\Gamma(a+b)}$. ∎

■ **What just happened.** Three results for the price of one: the Beta normalizing constant; $W = X/(X+Y) \sim \mathrm{Beta}(a, b)$ whenever $X, Y$ are independent Gammas with the same rate; and the fraction $W$ is *independent* of the total $T$. (Blitzstein calls this the bank–post office story.) The derivation was Theorem 3.6 plus Theorem 3.4(c), which is the pattern for most of this module. Moments: $\mathbb{E}[\mathrm{Beta}(a,b)] = \frac{B(a+1,b)}{B(a,b)} = \frac{a}{a+b}$ using $\Gamma(a+1) = a\Gamma(a)$.

**Theorem 3.12 (Beta–Binomial conjugacy).** Let $p \sim \mathrm{Beta}(a, b)$ be an unknown success probability, and given $p$ let $X \sim \mathrm{Binomial}(n, p)$. Then
$$p \mid X = k \;\sim\; \mathrm{Beta}(a + k,\; b + n - k), \qquad \mathbb{E}[p \mid X = k] = \frac{a + k}{a + b + n}. \tag{3.12}$$

*Derivation.* Bayes' rule for densities:
$$f(p \mid X = k) = \frac{P(X = k \mid p)\,f(p)}{P(X = k)} \propto \binom nk p^k(1-p)^{n-k}\cdot\frac{p^{a-1}(1-p)^{b-1}}{B(a,b)} \propto p^{a+k-1}(1-p)^{b+n-k-1}.$$
Everything not depending on $p$ was absorbed into the proportionality constant. The surviving shape is the kernel of $\mathrm{Beta}(a+k, b+n-k)$, and since the posterior is a density it must be exactly that. The mean follows from the Beta mean formula. ∎

■ **What just happened.** Your first Bayesian computation. The prior parameters $a, b$ act as counts of imaginary prior successes and failures; the data adds real ones. With $a = b = 1$ (uniform prior) the posterior mean is $(k+1)/(n+2)$, Laplace's rule of succession. The posterior mean is a weighted average of the prior mean $a/(a+b)$ (weight $a+b$) and the sample proportion $k/n$ (weight $n$); exercise B4 has you verify this. ST-5 is this computation repeated for every exponential family.

### 3.7 Covariance and correlation

**Definition 3.13.** For random variables $X, Y$ with means $\mu_X, \mu_Y$ and finite variances,
$$\mathrm{Cov}(X, Y) = \mathbb{E}\big[(X - \mu_X)(Y - \mu_Y)\big] = \mathbb{E}[XY] - \mathbb{E}[X]\mathbb{E}[Y]. \tag{3.13}$$
The second form follows by expanding the product and using linearity: $\mathbb{E}[XY - \mu_Y X - \mu_X Y + \mu_X\mu_Y] = \mathbb{E}[XY] - 2\mu_X\mu_Y + \mu_X\mu_Y$. Note $\mathrm{Cov}(X, X) = \mathrm{Var}(X)$, and by Corollary 3.5, $X \perp Y \Rightarrow \mathrm{Cov}(X, Y) = 0$.

**Theorem 3.14 (Bilinearity).** $\mathrm{Cov}(aX + bY, Z) = a\,\mathrm{Cov}(X, Z) + b\,\mathrm{Cov}(Y, Z)$, and symmetrically in the second slot. Also $\mathrm{Cov}(X + c, Y) = \mathrm{Cov}(X, Y)$ for constants $c$.

*Proof.* Using (3.13) and linearity of $\mathbb{E}$: $\mathbb{E}[(aX + bY)Z] - \mathbb{E}[aX + bY]\mathbb{E}[Z] = a(\mathbb{E}[XZ] - \mathbb{E}[X]\mathbb{E}[Z]) + b(\mathbb{E}[YZ] - \mathbb{E}[Y]\mathbb{E}[Z])$. Symmetry of (3.13) in $X, Y$ gives the second slot. Shifting by $c$ does not change $X - \mu_X$. ∎

**Corollary 3.15.** $\mathrm{Var}(X + Y) = \mathrm{Var}(X) + \mathrm{Var}(Y) + 2\,\mathrm{Cov}(X, Y)$, and more generally
$$\mathrm{Var}\Big(\sum_i a_i X_i\Big) = \sum_{i,j} a_i a_j\,\mathrm{Cov}(X_i, X_j) = \mathbf{a}^\top\Sigma\,\mathbf{a}, \tag{3.14}$$
where $\Sigma_{ij} = \mathrm{Cov}(X_i, X_j)$.

*Proof.* $\mathrm{Var}(\sum a_i X_i) = \mathrm{Cov}(\sum_i a_i X_i, \sum_j a_j X_j)$, expand both slots by bilinearity. ∎ For independent (or merely uncorrelated) variables the cross terms vanish and variances add. This is why $\mathrm{Var}(\bar X_n) = \sigma^2/n$, the fact behind every error bar.

**Definition 3.16 (Correlation).** $\rho = \mathrm{Corr}(X, Y) = \dfrac{\mathrm{Cov}(X,Y)}{\sigma_X\sigma_Y}$. It is the covariance of the standardized variables $(X-\mu_X)/\sigma_X$ and $(Y - \mu_Y)/\sigma_Y$, so it is dimensionless and unchanged by $X \mapsto aX + b$ with $a > 0$.

**Theorem 3.17.** $-1 \le \rho \le 1$, with $|\rho| = 1$ if and only if $Y = \alpha X + \beta$ for some constants (with probability 1).

▶ **Try it first.** Define $\langle U, V\rangle = \mathbb{E}[UV]$ on random variables with finite variance. Check the three inner-product properties from LA-1 Definition 3.9. Then apply Cauchy–Schwarz (LA-1 Theorem 3.11), whose proof used *only* those three properties.

**Derivation.**

1. $\langle U, V\rangle = \mathbb{E}[UV]$ is symmetric, linear in each slot (linearity of $\mathbb{E}$), and $\langle U, U\rangle = \mathbb{E}[U^2] \ge 0$ with equality iff $U = 0$ with probability 1. These are the three properties of an inner product. (LA-8 makes this a genuine inner-product space, called $L^2$; we do not need the name.)
2. The LA-1 proof of Cauchy–Schwarz used nothing except those three properties: for all $t$, $0 \le \langle U - tV, U - tV\rangle = \langle U,U\rangle - 2t\langle U,V\rangle + t^2\langle V,V\rangle$, so the discriminant is $\le 0$, so $\langle U, V\rangle^2 \le \langle U,U\rangle\langle V,V\rangle$.
3. Apply it to the centered variables $U = X - \mu_X$, $V = Y - \mu_Y$: $\mathrm{Cov}(X,Y)^2 = \langle U, V\rangle^2 \le \mathbb{E}[U^2]\mathbb{E}[V^2] = \sigma_X^2\sigma_Y^2$. Divide and take square roots: $|\rho| \le 1$.
4. Equality iff the discriminant is zero iff $\langle U - t_0 V, U - t_0 V\rangle = 0$ for some $t_0$ iff $U = t_0 V$ with probability 1 iff $Y - \mu_Y = (X - \mu_X)/t_0$, a linear relation. ∎

■ **What just happened.** Correlation is the cosine of the angle between two centered random variables, in the inner product $\mathbb{E}[UV]$. This is not an analogy; it is the same theorem as LA-1 Theorem 3.11, applied in a vector space whose "vectors" are random variables. Keep this picture: §3.9 will show that conditional expectation is orthogonal projection in the same space.

**Uncorrelated does not mean independent.** Let $X$ be uniform on $\{-1, 0, 1\}$ and $Y = X^2$. Then $\mathbb{E}[X] = 0$ and $\mathbb{E}[XY] = \mathbb{E}[X^3] = 0$, so $\mathrm{Cov}(X, Y) = 0$. But $P(X = 0, Y = 0) = \tfrac13$ while $P(X = 0)P(Y = 0) = \tfrac13\cdot\tfrac13 = \tfrac19$. Not independent: $Y$ is a *function* of $X$. The same works with $X \sim \mathcal{N}(0,1)$, $Y = X^2$: $\mathbb{E}[X^3] = 0$ by symmetry, so $\rho = 0$, yet $Y$ is determined by $X$. Covariance sees only linear dependence. The one important setting where uncorrelated *does* imply independent is the multivariate normal (Theorem 3.23).

**Definition 3.18 (Covariance matrix).** For a random vector $\mathbf{X} \in \mathbb{R}^k$ with mean $\boldsymbol\mu = \mathbb{E}[\mathbf{X}]$ (componentwise),
$$\Sigma = \mathrm{Cov}(\mathbf{X}) = \mathbb{E}\big[(\mathbf{X} - \boldsymbol\mu)(\mathbf{X} - \boldsymbol\mu)^\top\big], \qquad \Sigma_{ij} = \mathrm{Cov}(X_i, X_j). \tag{3.15}$$
The expectation of a matrix is taken entrywise. The diagonal holds the variances; $\Sigma$ is symmetric because $\mathrm{Cov}$ is.

**Theorem 3.19.** $\Sigma$ is positive semidefinite. It is positive definite unless some nontrivial linear combination $\mathbf{a}^\top\mathbf{X}$ is constant with probability 1. Moreover, for any matrix $A$ and vector $\mathbf{b}$, $\mathrm{Cov}(A\mathbf{X} + \mathbf{b}) = A\Sigma A^\top$.

*Proof.* For any $\mathbf{a} \in \mathbb{R}^k$, pull the constant vector through the expectation:
$$\mathbf{a}^\top\Sigma\,\mathbf{a} = \mathbb{E}\big[\mathbf{a}^\top(\mathbf{X} - \boldsymbol\mu)(\mathbf{X} - \boldsymbol\mu)^\top\mathbf{a}\big] = \mathbb{E}\big[(\mathbf{a}^\top(\mathbf{X} - \boldsymbol\mu))^2\big] = \mathrm{Var}(\mathbf{a}^\top\mathbf{X}) \ge 0.$$
This is (3.14) again. It is zero iff $\mathbf{a}^\top\mathbf{X}$ has zero variance, i.e., is constant. For the transformation rule: $A\mathbf{X} + \mathbf{b}$ has mean $A\boldsymbol\mu + \mathbf{b}$, so its covariance is $\mathbb{E}[A(\mathbf{X}-\boldsymbol\mu)(\mathbf{X}-\boldsymbol\mu)^\top A^\top] = A\Sigma A^\top$. ∎

■ **What just happened.** Every covariance matrix is symmetric PSD, so LA-8's spectral theorem applies: $\Sigma = Q\Lambda Q^\top$ with orthonormal eigenvectors and eigenvalues $\lambda_i \ge 0$. The eigenvector $\mathbf{q}_1$ is the direction in which $\mathbf{X}$ varies most (the Rayleigh quotient from LA-8: $\max_{\|\mathbf{a}\|=1}\mathrm{Var}(\mathbf{a}^\top\mathbf{X}) = \lambda_1$). That is PCA (LA-9, ML-6) in one sentence. Conversely every symmetric PSD matrix is a covariance matrix: take $\mathbf{X} = \Sigma^{1/2}\mathbf{Z}$ for any $\mathbf{Z}$ with $\mathrm{Cov}(\mathbf{Z}) = I$ and apply the transformation rule.

### 3.8 The multivariate normal

**Definition 3.20.** Let $\mathbf{Z} = (Z_1, \dots, Z_k)^\top$ with $Z_i$ iid $\mathcal{N}(0, 1)$. A random vector $\mathbf{X} \in \mathbb{R}^k$ is **multivariate normal** (MVN, or jointly Gaussian) if $\mathbf{X} = \boldsymbol\mu + A\mathbf{Z}$ for some vector $\boldsymbol\mu$ and some $k \times m$ matrix $A$. By Theorem 3.19, $\mathbb{E}[\mathbf{X}] = \boldsymbol\mu$ and $\mathrm{Cov}(\mathbf{X}) = AA^\top =: \Sigma$. We write $\mathbf{X} \sim \mathcal{N}(\boldsymbol\mu, \Sigma)$.

The density of $\mathbf{Z}$ is the product of $k$ standard normal densities (Theorem 3.4(a)):
$$f_{\mathbf{Z}}(\mathbf{z}) = \prod_{i=1}^k\frac{1}{\sqrt{2\pi}}e^{-z_i^2/2} = (2\pi)^{-k/2}\exp\!\left(-\tfrac12\|\mathbf{z}\|^2\right). \tag{3.16}$$

**Theorem 3.21 (Density of the MVN).** If $A$ is $k\times k$ and invertible (equivalently, $\Sigma$ is positive definite), then
$$f_{\mathbf{X}}(\mathbf{x}) = \frac{1}{\sqrt{(2\pi)^k\det\Sigma}}\exp\!\left(-\tfrac12(\mathbf{x} - \boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{x} - \boldsymbol\mu)\right). \tag{3.17}$$

▶ **Try it first.** The map $\mathbf{z} \mapsto \mathbf{x} = \boldsymbol\mu + A\mathbf{z}$ is linear plus a shift. What is the Jacobian of its inverse? How is $\det A$ related to $\det\Sigma$?

**Derivation.** This is Theorem 3.6 in $k$ dimensions (the proof there did not use $k = 2$).

1. The inverse map is $\mathbf{z} = A^{-1}(\mathbf{x} - \boldsymbol\mu)$. Its Jacobian matrix $\partial\mathbf{z}/\partial\mathbf{x}$ is the constant matrix $A^{-1}$, so $|\det\partial\mathbf{z}/\partial\mathbf{x}| = |\det A|^{-1}$.
2. $\det\Sigma = \det(AA^\top) = \det A\cdot\det A^\top = (\det A)^2$ (LA-6), so $|\det A| = \sqrt{\det\Sigma}$.
3. By (3.5), $f_{\mathbf{X}}(\mathbf{x}) = f_{\mathbf{Z}}(A^{-1}(\mathbf{x} - \boldsymbol\mu))\,|\det A|^{-1}$. Substitute (3.16):
$$\|A^{-1}(\mathbf{x} - \boldsymbol\mu)\|^2 = (\mathbf{x} - \boldsymbol\mu)^\top (A^{-1})^\top A^{-1}(\mathbf{x} - \boldsymbol\mu) = (\mathbf{x} - \boldsymbol\mu)^\top(AA^\top)^{-1}(\mathbf{x} - \boldsymbol\mu) = (\mathbf{x} - \boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{x} - \boldsymbol\mu),$$
using $(A^{-1})^\top A^{-1} = (A^\top)^{-1}A^{-1} = (AA^\top)^{-1}$ (LA-2: inverse of a product reverses the order). Assemble with the constant $(2\pi)^{-k/2}(\det\Sigma)^{-1/2}$. ∎

■ **What just happened.** The density depends on $A$ only through $\Sigma = AA^\top$. So two different square roots of $\Sigma$ (say the Cholesky factor from LA-10 and the symmetric square root $Q\Lambda^{1/2}Q^\top$ from LA-8) produce the *same distribution*. The distribution is parametrized by $(\boldsymbol\mu, \Sigma)$, which is why we write $\mathcal{N}(\boldsymbol\mu, \Sigma)$ and not $\mathcal{N}(\boldsymbol\mu, A)$. For $k = 1$, (3.17) is the familiar $\frac{1}{\sqrt{2\pi\sigma^2}}e^{-(x-\mu)^2/2\sigma^2}$.

**Geometry.** The level sets of (3.17) are the sets where the quadratic form $(\mathbf{x} - \boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{x} - \boldsymbol\mu) = c$ is constant. Write $\Sigma = Q\Lambda Q^\top$ (LA-8 spectral theorem; $Q$ orthogonal, $\Lambda = \mathrm{diag}(\lambda_1, \dots, \lambda_k)$, all $\lambda_i > 0$). Then $\Sigma^{-1} = Q\Lambda^{-1}Q^\top$. Put $\mathbf{y} = Q^\top(\mathbf{x} - \boldsymbol\mu)$ (coordinates along the eigenvectors). The level set becomes
$$\sum_{i=1}^k\frac{y_i^2}{\lambda_i} = c, \qquad\text{i.e.}\qquad \sum_i\frac{y_i^2}{(\sqrt{c\lambda_i})^2} = 1,$$
an ellipsoid centered at $\boldsymbol\mu$ with axes along the eigenvectors $\mathbf{q}_i$ and semi-axis lengths $\sqrt{c\lambda_i}$. **The axes of the Gaussian are the eigenvectors of $\Sigma$; the semi-axes are proportional to $\sqrt{\lambda_i}$.** This is the guess from §2 confirmed, and it is PCA's picture (ML-6): the first principal component is the long axis of the ellipse.

**The MGF of a random vector.** Define $M_{\mathbf{X}}(\mathbf{t}) = \mathbb{E}[e^{\mathbf{t}^\top\mathbf{X}}]$ for $\mathbf{t} \in \mathbb{R}^k$. For $\mathbf{Z}$: $M_{\mathbf{Z}}(\mathbf{t}) = \mathbb{E}[\prod_i e^{t_iZ_i}] = \prod_i e^{t_i^2/2} = e^{\|\mathbf{t}\|^2/2}$ (independence and the F4 normal MGF). For $\mathbf{X} = \boldsymbol\mu + A\mathbf{Z}$:
$$M_{\mathbf{X}}(\mathbf{t}) = \mathbb{E}\big[e^{\mathbf{t}^\top\boldsymbol\mu + (A^\top\mathbf{t})^\top\mathbf{Z}}\big] = e^{\mathbf{t}^\top\boldsymbol\mu}M_{\mathbf{Z}}(A^\top\mathbf{t}) = \exp\!\left(\mathbf{t}^\top\boldsymbol\mu + \tfrac12\mathbf{t}^\top AA^\top\mathbf{t}\right) = \exp\!\left(\mathbf{t}^\top\boldsymbol\mu + \tfrac12\mathbf{t}^\top\Sigma\,\mathbf{t}\right). \tag{3.18}$$
This holds for *any* $A$, square or not, invertible or not, and depends only on $(\boldsymbol\mu, \Sigma)$. Since the MGF determines the distribution (the multivariate version of the F4 uniqueness theorem; we cite it as we did there), Definition 3.20 defines a single distribution for each $(\boldsymbol\mu, \Sigma)$ even when $\Sigma$ is singular and no density exists.

**Theorem 3.22 (Linear transforms of MVN are MVN).** If $\mathbf{X} \sim \mathcal{N}(\boldsymbol\mu, \Sigma)$, $B$ is $m \times k$ and $\mathbf{c} \in \mathbb{R}^m$, then $B\mathbf{X} + \mathbf{c} \sim \mathcal{N}(B\boldsymbol\mu + \mathbf{c},\; B\Sigma B^\top)$.

*Proof.* $M_{B\mathbf{X}+\mathbf{c}}(\mathbf{s}) = \mathbb{E}[e^{\mathbf{s}^\top(B\mathbf{X} + \mathbf{c})}] = e^{\mathbf{s}^\top\mathbf{c}}M_{\mathbf{X}}(B^\top\mathbf{s}) = \exp\big(\mathbf{s}^\top(B\boldsymbol\mu + \mathbf{c}) + \tfrac12\mathbf{s}^\top B\Sigma B^\top\mathbf{s}\big)$ by (3.18), which is (3.18) for $\mathcal{N}(B\boldsymbol\mu + \mathbf{c}, B\Sigma B^\top)$. ∎ (Or directly: $B\mathbf{X} + \mathbf{c} = (B\boldsymbol\mu + \mathbf{c}) + (BA)\mathbf{Z}$ is of the form in Definition 3.20.)

**Corollary (marginals are Gaussian).** Partition $\mathbf{X} = (\mathbf{X}_1, \mathbf{X}_2)$ into the first $p$ and last $q$ coordinates, with $\boldsymbol\mu = (\boldsymbol\mu_1, \boldsymbol\mu_2)$ and $\Sigma = \begin{pmatrix}\Sigma_{11} & \Sigma_{12}\\ \Sigma_{21} & \Sigma_{22}\end{pmatrix}$. Take $B = (I_p\;\; 0)$ in Theorem 3.22: $\mathbf{X}_1 \sim \mathcal{N}(\boldsymbol\mu_1, \Sigma_{11})$. In particular every coordinate $X_i \sim \mathcal{N}(\mu_i, \Sigma_{ii})$ and every linear combination $\mathbf{a}^\top\mathbf{X} \sim \mathcal{N}(\mathbf{a}^\top\boldsymbol\mu, \mathbf{a}^\top\Sigma\mathbf{a})$. (The converse is the standard alternative definition: $\mathbf{X}$ is MVN iff every linear combination is univariate normal.) **Warning:** Gaussian marginals do *not* imply a Gaussian joint; exercise B5 gives a counterexample.

**Theorem 3.23 (For MVN, uncorrelated ⇒ independent).** If $\mathbf{X} \sim \mathcal{N}(\boldsymbol\mu, \Sigma)$ with $\Sigma_{12} = 0$ (so $\Sigma$ is block diagonal), then $\mathbf{X}_1 \perp \mathbf{X}_2$.

*Proof.* Assume first $\Sigma$ is positive definite so (3.17) applies. Block diagonal $\Sigma$ has $\Sigma^{-1} = \mathrm{diag}(\Sigma_{11}^{-1}, \Sigma_{22}^{-1})$ and $\det\Sigma = \det\Sigma_{11}\det\Sigma_{22}$ (LA-6). So the quadratic form splits, $(\mathbf{x}-\boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{x}-\boldsymbol\mu) = (\mathbf{x}_1-\boldsymbol\mu_1)^\top\Sigma_{11}^{-1}(\mathbf{x}_1-\boldsymbol\mu_1) + (\mathbf{x}_2-\boldsymbol\mu_2)^\top\Sigma_{22}^{-1}(\mathbf{x}_2-\boldsymbol\mu_2)$, and the exponential of a sum is a product: (3.17) factors as $f(\mathbf{x}_1, \mathbf{x}_2) = g(\mathbf{x}_1)h(\mathbf{x}_2)$, with $g, h$ the $\mathcal{N}(\boldsymbol\mu_1, \Sigma_{11})$ and $\mathcal{N}(\boldsymbol\mu_2, \Sigma_{22})$ densities. Theorem 3.4(c) gives independence. In the singular case the same factorization happens in the MGF (3.18): $M_{\mathbf{X}}(\mathbf{t}_1, \mathbf{t}_2) = M_{\mathbf{X}_1}(\mathbf{t}_1)M_{\mathbf{X}_2}(\mathbf{t}_2)$, which also characterizes independence. ∎

■ **What just happened.** The Gaussian's dependence structure lives *entirely* in $\Sigma$. Zero covariance kills the cross terms in the exponent, and the density factors. For the $X, X^2$ counterexample of §3.7 this fails because $(X, X^2)$ is not jointly Gaussian, even though $X$ is.

**Theorem 3.24 (Conditional distribution of a Gaussian).** Let $\mathbf{X} = (\mathbf{X}_1, \mathbf{X}_2) \sim \mathcal{N}(\boldsymbol\mu, \Sigma)$ with $\Sigma$ positive definite, partitioned as above. Then
$$\mathbf{X}_1 \mid \mathbf{X}_2 = \mathbf{x}_2 \;\sim\; \mathcal{N}\!\Big(\boldsymbol\mu_1 + \Sigma_{12}\Sigma_{22}^{-1}(\mathbf{x}_2 - \boldsymbol\mu_2),\;\; \Sigma_{11} - \Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21}\Big). \tag{3.19}$$
The matrix $S := \Sigma_{11} - \Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21}$ is the **Schur complement** of $\Sigma_{22}$ in $\Sigma$ (LA-10).

▶ **Try it first.** You need $f(\mathbf{x}_1 \mid \mathbf{x}_2) = f(\mathbf{x}_1, \mathbf{x}_2)/f(\mathbf{x}_2)$, where the numerator is (3.17) and the denominator is the marginal $\mathcal{N}(\boldsymbol\mu_2, \Sigma_{22})$. Everything hinges on writing the big quadratic form $(\mathbf{x}-\boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{x}-\boldsymbol\mu)$ as [something in $\mathbf{x}_1$ given $\mathbf{x}_2$] + [the marginal's quadratic form in $\mathbf{x}_2$]. That is "completing the square" in block form. Start by doing block elimination (LA-2) on $\Sigma$: clear the $\Sigma_{21}$ block.

**Derivation.** Write $\mathbf{a} = \mathbf{x}_1 - \boldsymbol\mu_1$, $\mathbf{b} = \mathbf{x}_2 - \boldsymbol\mu_2$, and $C = \Sigma_{12}\Sigma_{22}^{-1}$ (a $p\times q$ matrix; note $C^\top = \Sigma_{22}^{-1}\Sigma_{21}$ since $\Sigma$ is symmetric).

1. **Block LDU factorization of $\Sigma$.** Multiply out to verify:
$$\Sigma = \begin{pmatrix} I & C \\ 0 & I\end{pmatrix}\begin{pmatrix} S & 0 \\ 0 & \Sigma_{22}\end{pmatrix}\begin{pmatrix} I & 0 \\ C^\top & I\end{pmatrix}.$$
The last two factors give $\begin{pmatrix} S & 0 \\ \Sigma_{22}C^\top & \Sigma_{22}\end{pmatrix} = \begin{pmatrix} S & 0 \\ \Sigma_{21} & \Sigma_{22}\end{pmatrix}$ since $\Sigma_{22}\Sigma_{22}^{-1}\Sigma_{21} = \Sigma_{21}$. Then the first factor times that is $\begin{pmatrix} S + C\Sigma_{21} & C\Sigma_{22} \\ \Sigma_{21} & \Sigma_{22}\end{pmatrix} = \begin{pmatrix}\Sigma_{11} & \Sigma_{12} \\ \Sigma_{21} & \Sigma_{22}\end{pmatrix}$, because $S + \Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21} = \Sigma_{11}$ and $\Sigma_{12}\Sigma_{22}^{-1}\Sigma_{22} = \Sigma_{12}$. ✓ (This is exactly block Gaussian elimination from LA-2: the Schur complement is what remains in the $(1,1)$ block after eliminating with the $(2,2)$ pivot.)
2. **Determinant.** The two triangular factors have determinant 1 (block triangular with identity diagonal blocks), so $\det\Sigma = \det S\cdot\det\Sigma_{22}$. In particular $S$ is positive definite (both other determinants are positive and $S$ is a principal-block Schur complement of a PD matrix; or directly: $S = \mathrm{Cov}$ of something, as step 6 shows).
3. **Inverse.** Inverting a product reverses the order (LA-2), and $\begin{pmatrix} I & C \\ 0 & I\end{pmatrix}^{-1} = \begin{pmatrix} I & -C \\ 0 & I\end{pmatrix}$ (check: the product is $\begin{pmatrix} I & -C + C \\ 0 & I\end{pmatrix} = I$). Hence
$$\Sigma^{-1} = \begin{pmatrix} I & 0 \\ -C^\top & I\end{pmatrix}\begin{pmatrix} S^{-1} & 0 \\ 0 & \Sigma_{22}^{-1}\end{pmatrix}\begin{pmatrix} I & -C \\ 0 & I\end{pmatrix}.$$
4. **The quadratic form splits.** Apply the rightmost factor to $(\mathbf{a}; \mathbf{b})$: $\begin{pmatrix} I & -C \\ 0 & I\end{pmatrix}\begin{pmatrix}\mathbf{a}\\ \mathbf{b}\end{pmatrix} = \begin{pmatrix}\mathbf{a} - C\mathbf{b} \\ \mathbf{b}\end{pmatrix}$. The leftmost factor is the transpose of the rightmost, so
$$(\mathbf{x} - \boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{x} - \boldsymbol\mu) = (\mathbf{a} - C\mathbf{b})^\top S^{-1}(\mathbf{a} - C\mathbf{b}) + \mathbf{b}^\top\Sigma_{22}^{-1}\mathbf{b}. \tag{3.20}$$
This is the completed square. The second term is exactly the quadratic form in the marginal density of $\mathbf{X}_2 \sim \mathcal{N}(\boldsymbol\mu_2, \Sigma_{22})$.
5. **Divide.** Using (3.17) for the joint ($k = p + q$) and for the marginal ($k = q$), with $\det\Sigma = \det S\det\Sigma_{22}$ from step 2:
$$f(\mathbf{x}_1 \mid \mathbf{x}_2) = \frac{(2\pi)^{-(p+q)/2}(\det S\det\Sigma_{22})^{-1/2}\exp\!\big(-\tfrac12[(\mathbf{a} - C\mathbf{b})^\top S^{-1}(\mathbf{a} - C\mathbf{b}) + \mathbf{b}^\top\Sigma_{22}^{-1}\mathbf{b}]\big)}{(2\pi)^{-q/2}(\det\Sigma_{22})^{-1/2}\exp\!\big(-\tfrac12\mathbf{b}^\top\Sigma_{22}^{-1}\mathbf{b}\big)} = \frac{\exp\!\big(-\tfrac12(\mathbf{a} - C\mathbf{b})^\top S^{-1}(\mathbf{a} - C\mathbf{b})\big)}{\sqrt{(2\pi)^p\det S}}.$$
6. **Read it off.** As a function of $\mathbf{x}_1$ this is the $\mathcal{N}(\boldsymbol\mu_1 + C\mathbf{b}, S)$ density, since $\mathbf{a} - C\mathbf{b} = \mathbf{x}_1 - (\boldsymbol\mu_1 + \Sigma_{12}\Sigma_{22}^{-1}(\mathbf{x}_2 - \boldsymbol\mu_2))$. ∎

■ **What just happened.** Three facts to memorize from (3.19). (i) The conditional mean is **linear in $\mathbf{x}_2$**, with coefficient matrix $\Sigma_{12}\Sigma_{22}^{-1}$. In one dimension this is $\mu_1 + \rho\frac{\sigma_1}{\sigma_2}(x_2 - \mu_2)$: **this formula is linear regression** (ST-3), and §3.9 shows it is the best linear predictor for *any* joint distribution, Gaussian or not. (ii) The conditional covariance $S$ **does not depend on $\mathbf{x}_2$**: observing $\mathbf{X}_2$ shrinks your uncertainty about $\mathbf{X}_1$ by the same amount whatever value you observe. In 1D, $S = \sigma_1^2(1 - \rho^2)$. (iii) The conditional covariance is $\Sigma_{11}$ minus something PSD: conditioning never increases Gaussian uncertainty. Equation (3.19) **is the Kalman filter update** (ST-8: $\mathbf{X}_1$ = state, $\mathbf{X}_2$ = measurement, $\Sigma_{12}\Sigma_{22}^{-1}$ = Kalman gain) and **is Gaussian-process regression** (ML-4: $\mathbf{X}_1$ = function values at test points, $\mathbf{X}_2$ = observed values, $\Sigma$ = kernel matrix). You will do this derivation again in ST-8; it is worth doing twice.

**Sampling.** To draw $\mathbf{X} \sim \mathcal{N}(\boldsymbol\mu, \Sigma)$: compute the Cholesky factor $L$ with $LL^\top = \Sigma$ (LA-10; exists since $\Sigma$ is PD), draw $\mathbf{Z}$ with iid $\mathcal{N}(0,1)$ entries (Box–Muller), return $\boldsymbol\mu + L\mathbf{Z}$. By Definition 3.20 with $A = L$, this has the right distribution. Any $A$ with $AA^\top = \Sigma$ works; Cholesky is the cheapest.

**Mahalanobis distance.** $d(\mathbf{x}) = \sqrt{(\mathbf{x} - \boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{x} - \boldsymbol\mu)}$. Since $\Sigma^{-1} = L^{-\top}L^{-1}$, $d(\mathbf{x}) = \|L^{-1}(\mathbf{x} - \boldsymbol\mu)\| = \|\mathbf{z}\|$: it is the ordinary Euclidean distance *after* undoing the transformation that made the Gaussian elliptical. It measures "how many standard deviations away, in the right direction." If $\mathbf{X} \sim \mathcal{N}(\boldsymbol\mu, \Sigma)$ then $d(\mathbf{X})^2 = \|\mathbf{Z}\|^2 \sim \chi^2_k$ (§3.5), which is how you set a confidence ellipse: the ellipse $d(\mathbf{x})^2 \le \chi^2_{k, 0.95}$ contains 95% of the mass. For $k = 2$, $\chi^2_{2,0.95} = -2\ln 0.05 \approx 5.99$ (since $\chi^2_2$ is $\mathrm{Exponential}(\tfrac12)$). The weighted inner product of LA-1 exercise C1 with weight matrix $\Sigma^{-1}$ is exactly this.

### 3.9 Conditional expectation

**Definition 3.25.** For each $x$ with $p_X(x) > 0$, the **conditional expectation of $Y$ given $X = x$** is the mean of the conditional distribution (3.3)/(3.4):
$$g(x) := \mathbb{E}[Y \mid X = x] = \sum_y y\,p_{Y|X}(y \mid x) \quad\text{or}\quad \int y\,f_{Y|X}(y \mid x)\,dy.$$
This is an ordinary function of $x$. The **conditional expectation $\mathbb{E}[Y \mid X]$** is the random variable $g(X)$: plug the random $X$ into the function $g$.

The second sentence is the conceptual leap (§8). $\mathbb{E}[Y \mid X = 3]$ is a number. $\mathbb{E}[Y \mid X]$ is a random variable, a function of $X$, whose value is $\mathbb{E}[Y\mid X = 3]$ on the event $X = 3$. It is "the best guess of $Y$ once you know $X$," considered before you know $X$.

Three properties, each immediate from the definition since conditioning on $X = x$ means computing with the conditional distribution:

- **Linearity:** $\mathbb{E}[aY + bW \mid X] = a\,\mathbb{E}[Y \mid X] + b\,\mathbb{E}[W \mid X]$.
- **Taking out what is known:** $\mathbb{E}[h(X)\,Y \mid X] = h(X)\,\mathbb{E}[Y \mid X]$. Given $X = x$, $h(X)$ is the constant $h(x)$, and constants come out of expectations. In particular $\mathbb{E}[h(X) \mid X] = h(X)$.
- **Independence drops the conditioning:** if $X \perp Y$ then $p_{Y|X}(y\mid x) = p_Y(y)$ (Theorem 3.4(b)), so $\mathbb{E}[Y \mid X] = \mathbb{E}[Y]$, a constant.

**Theorem 3.26 (Adam's law, or the tower property).** $\mathbb{E}\big[\mathbb{E}[Y \mid X]\big] = \mathbb{E}[Y]$.

▶ **Try it first.** $\mathbb{E}[Y\mid X]$ is $g(X)$; compute $\mathbb{E}[g(X)]$ by LOTUS and swap the order of summation.

**Derivation.** $\mathbb{E}[g(X)] = \sum_x g(x)\,p_X(x) = \sum_x\Big(\sum_y y\,p_{Y|X}(y\mid x)\Big)p_X(x) = \sum_x\sum_y y\,p_{X,Y}(x,y) = \sum_y y\sum_x p_{X,Y}(x,y) = \sum_y y\,p_Y(y) = \mathbb{E}[Y]$, using the chain rule and then (3.1). Integrals for the continuous case. ∎ The same argument with extra conditioning gives $\mathbb{E}\big[\mathbb{E}[Y\mid X, W]\,\big|\,W\big] = \mathbb{E}[Y\mid W]$: averaging out the finer information leaves the coarser.

■ **What just happened.** To compute $\mathbb{E}[Y]$, you may first compute it on each slice $X = x$ and then average the slice answers with the weights $p_X(x)$. This is the law of total probability for expectations, and it is the single most useful computational trick in probability: condition on the thing that makes the problem easy, then average.

**Definition 3.27.** $\mathrm{Var}(Y \mid X) := \mathbb{E}[Y^2 \mid X] - (\mathbb{E}[Y\mid X])^2$, the variance of the conditional distribution, as a random variable (function of $X$).

**Theorem 3.28 (Eve's law, or the law of total variance).**
$$\mathrm{Var}(Y) = \mathbb{E}\big[\mathrm{Var}(Y \mid X)\big] + \mathrm{Var}\big(\mathbb{E}[Y\mid X]\big). \tag{3.21}$$

▶ **Try it first.** Take $\mathbb{E}$ of Definition 3.27 using Adam's law, write down $\mathrm{Var}(g(X))$ by the usual formula, and add.

**Derivation.** Write $g(X) = \mathbb{E}[Y\mid X]$.

1. $\mathbb{E}[\mathrm{Var}(Y\mid X)] = \mathbb{E}\big[\mathbb{E}[Y^2\mid X]\big] - \mathbb{E}[g(X)^2] = \mathbb{E}[Y^2] - \mathbb{E}[g(X)^2]$, by Adam's law applied to $Y^2$.
2. $\mathrm{Var}(g(X)) = \mathbb{E}[g(X)^2] - (\mathbb{E}[g(X)])^2 = \mathbb{E}[g(X)^2] - (\mathbb{E}[Y])^2$, by Adam's law applied to $Y$.
3. Add: $\mathbb{E}[Y^2] - (\mathbb{E}[Y])^2 = \mathrm{Var}(Y)$. ∎

■ **What just happened.** Total variance = average *within*-slice variance (how spread $Y$ is once $X$ is known; the "unexplained" part) + variance of the slice *means* (how much knowing $X$ moves your guess; the "explained" part). The ratio $\mathrm{Var}(\mathbb{E}[Y\mid X])/\mathrm{Var}(Y)$ is the fraction of variance explained by $X$; for a jointly Gaussian pair it equals $\rho^2$ (check: $\mathrm{Var}(\mathbb{E}[Y\mid X]) = \mathrm{Var}(\rho\frac{\sigma_Y}{\sigma_X}X) = \rho^2\sigma_Y^2$). This is $R^2$ in ST-3 and the ANOVA decomposition in ST-2. Students swap the two terms constantly (§8); say the words "average of the conditional variances, plus variance of the conditional means."

**Theorem 3.29 (Conditional expectation is the best predictor in mean square; projection).** Among all functions $h$, the mean squared error $\mathbb{E}[(Y - h(X))^2]$ is minimized by $h(X) = \mathbb{E}[Y\mid X]$. Moreover, the residual $Y - \mathbb{E}[Y\mid X]$ is orthogonal to every function of $X$:
$$\mathbb{E}\big[(Y - \mathbb{E}[Y\mid X])\,h(X)\big] = 0 \quad\text{for all } h. \tag{3.22}$$

▶ **Try it first.** Prove (3.22) with Adam's law and "taking out what is known." Then write $Y - h(X) = (Y - \mathbb{E}[Y\mid X]) + (\mathbb{E}[Y\mid X] - h(X))$ and expand the square. Which cross term vanishes, and why?

**Derivation.**

1. Orthogonality. $\mathbb{E}\big[(Y - g(X))h(X)\big] = \mathbb{E}\Big[\mathbb{E}\big[(Y - g(X))h(X)\,\big|\,X\big]\Big]$ by Adam's law $= \mathbb{E}\big[h(X)\,\mathbb{E}[Y - g(X)\mid X]\big]$ by taking out what is known $= \mathbb{E}\big[h(X)(g(X) - g(X))\big] = 0$, since $\mathbb{E}[Y\mid X] = g(X)$ and $\mathbb{E}[g(X)\mid X] = g(X)$.
2. Pythagoras. Let $R = Y - g(X)$ (the residual) and $D = g(X) - h(X)$ (a function of $X$). Then $Y - h(X) = R + D$ and
$$\mathbb{E}[(Y - h(X))^2] = \mathbb{E}[R^2] + 2\,\mathbb{E}[RD] + \mathbb{E}[D^2] = \mathbb{E}[R^2] + \mathbb{E}[D^2] \ge \mathbb{E}[R^2],$$
because $\mathbb{E}[RD] = 0$ by step 1 with $h = D$. Equality iff $\mathbb{E}[D^2] = 0$ iff $h(X) = g(X)$ with probability 1. ∎

■ **What just happened.** In the inner product $\langle U, V\rangle = \mathbb{E}[UV]$ of Theorem 3.17, the set of all (square-integrable) functions of $X$ is a subspace. Equation (3.22) says the residual $Y - \mathbb{E}[Y\mid X]$ is orthogonal to that whole subspace. So $\mathbb{E}[Y\mid X]$ is the **orthogonal projection of $Y$ onto the functions of $X$** (LA-5: the projection is the unique vector in the subspace whose residual is orthogonal to the subspace), and the mean-square optimality is LA-5's "projection is the closest point," proved by the same Pythagoras argument as LA-1 exercise B7(c). This is the geometric heart of regression: *regression is projection*. ST-3 does it with a finite data vector and the column space of a design matrix; here it is the infinite-dimensional original.

**The best linear predictor.** Suppose you restrict $h$ to lines, $h(X) = \alpha + \beta X$. Minimize $J(\alpha, \beta) = \mathbb{E}[(Y - \alpha - \beta X)^2]$.

▶ **Try it first.** Set both partial derivatives to zero. They are two linear equations in $(\alpha, \beta)$; the $\alpha$ equation is "the residual has mean zero."

**Derivation.**

1. $\partial J/\partial\alpha = -2\,\mathbb{E}[Y - \alpha - \beta X] = 0 \Rightarrow \alpha = \mu_Y - \beta\mu_X$.
2. $\partial J/\partial\beta = -2\,\mathbb{E}[X(Y - \alpha - \beta X)] = 0 \Rightarrow \mathbb{E}[XY] - \alpha\mu_X - \beta\,\mathbb{E}[X^2] = 0$. Substitute step 1: $\mathbb{E}[XY] - \mu_X\mu_Y + \beta\mu_X^2 - \beta\,\mathbb{E}[X^2] = 0$, i.e., $\mathrm{Cov}(X,Y) - \beta\,\mathrm{Var}(X) = 0$.
3. Therefore
$$\beta^* = \frac{\mathrm{Cov}(X, Y)}{\mathrm{Var}(X)} = \rho\,\frac{\sigma_Y}{\sigma_X}, \qquad \hat Y = \mu_Y + \rho\frac{\sigma_Y}{\sigma_X}(X - \mu_X). \tag{3.23}$$
($J$ is a convex quadratic in $(\alpha, \beta)$, so the stationary point is the minimum.) The residual variance is $\mathrm{Var}(Y - \beta^*X) = \sigma_Y^2 - 2\beta^*\mathrm{Cov}(X,Y) + \beta^{*2}\sigma_X^2 = \sigma_Y^2 - 2\rho^2\sigma_Y^2 + \rho^2\sigma_Y^2 = \sigma_Y^2(1 - \rho^2)$. ∎

■ **What just happened.** (3.23) is (3.19) in one dimension. For a jointly Gaussian pair, the best linear predictor *is* the conditional expectation (because the conditional mean happens to be linear). For any other joint distribution, (3.23) is still the best line, but $\mathbb{E}[Y\mid X]$ may be a curve that does better. Now **regression to the mean**: in standardized units, (3.23) reads $\hat z_Y = \rho\,z_X$. A parent two standard deviations tall ($z_X = 2$) has a child predicted at $2\rho$ standard deviations, which is *less* than 2 whenever $\rho < 1$. Nothing is pulling anyone toward the mean: it is the statement that the line of best fit has slope $\rho < 1$ in standardized units, and it holds equally for predicting parents from children. Galton's "force" was an artifact of imperfect correlation. Example 4.5 does the numbers.

### 3.10 Inequalities

**Theorem 3.30 (Markov).** If $X \ge 0$ and $a > 0$, then $P(X \ge a) \le \mathbb{E}[X]/a$.

*Derivation.* For every outcome, $a\cdot\mathbf{1}\{X \ge a\} \le X$: if $X \ge a$ the left side is $a \le X$; otherwise the left side is $0 \le X$. Take expectations: $a\,P(X \ge a) \le \mathbb{E}[X]$. ∎

**Theorem 3.31 (Chebyshev).** If $X$ has mean $\mu$ and variance $\sigma^2$, then for any $k > 0$, $P(|X - \mu| \ge k) \le \sigma^2/k^2$.

*Derivation.* $|X - \mu| \ge k$ iff $(X - \mu)^2 \ge k^2$. Apply Markov to the nonnegative variable $(X - \mu)^2$ with $a = k^2$: $P((X-\mu)^2 \ge k^2) \le \mathbb{E}[(X-\mu)^2]/k^2 = \sigma^2/k^2$. ∎

Both are crude (for a Gaussian, $P(|X - \mu| \ge 2\sigma) \approx 0.046$ while Chebyshev gives $0.25$) but they assume *nothing* about the shape of the distribution, and that is what the law of large numbers needs.

**Theorem 3.32 (Jensen).** If $g$ is convex and $\mathbb{E}[X]$ exists, then $\mathbb{E}[g(X)] \ge g(\mathbb{E}[X])$. If $g$ is concave the inequality reverses. For strictly convex $g$, equality holds iff $X$ is constant.

*Two-point case.* If $X = x_1$ with probability $p$ and $x_2$ with probability $1 - p$, then $\mathbb{E}[g(X)] = p\,g(x_1) + (1-p)\,g(x_2) \ge g(p x_1 + (1-p)x_2) = g(\mathbb{E}[X])$. That inequality *is* the definition of convexity: the chord lies above the graph.

*General case.* Let $\mu = \mathbb{E}[X]$. A convex function has a supporting line at every point: there is a constant $c$ with $g(x) \ge g(\mu) + c(x - \mu)$ for all $x$ (the tangent line if $g$ is differentiable; otherwise any slope between the one-sided derivatives). Take expectations of both sides: $\mathbb{E}[g(X)] \ge g(\mu) + c\,\mathbb{E}[X - \mu] = g(\mu)$. For strictly convex $g$ the supporting line touches only at $\mu$, so equality forces $X = \mu$ with probability 1. ∎

Uses: $g(x) = x^2$ gives $\mathbb{E}[X^2] \ge (\mathbb{E}X)^2$, i.e., $\mathrm{Var} \ge 0$. $g = -\log$ (convex) gives $\mathbb{E}[\log X] \le \log\mathbb{E}[X]$, which is the proof of Gibbs' inequality in §3.13. Cauchy–Schwarz (Theorem 3.17) completes the set; together these four are most of the inequalities you will use in ST and ML.

### 3.11 Limit theorems

**Definition 3.33.** $X_n \to X$ **in probability** ($X_n \xrightarrow{p} X$) if for every $\epsilon > 0$, $P(|X_n - X| > \epsilon) \to 0$. $X_n \to X$ **in distribution** ($X_n \xrightarrow{d} X$) if $F_{X_n}(x) \to F_X(x)$ at every $x$ where $F_X$ is continuous.

The difference in one sentence each: convergence in probability says the random variables themselves become close (realization by realization, with high probability), while convergence in distribution says only that their *histograms* become close ($X_n$ and $X$ need not even live on the same sample space). In probability implies in distribution; the converse fails: if $Z \sim \mathcal{N}(0,1)$ and $X_n = (-1)^n Z$, every $X_n$ is exactly $\mathcal{N}(0,1)$, so $X_n \xrightarrow{d} Z$, but $|X_n - Z| = 2|Z|$ for odd $n$, so there is no convergence in probability. Almost-sure convergence ($P(\lim X_n = X) = 1$) is stronger than both.

**Theorem 3.34 (Weak law of large numbers).** Let $X_1, X_2, \dots$ be iid with mean $\mu$ and finite variance $\sigma^2$, and $\bar X_n = \frac1n\sum_{i=1}^n X_i$. Then $\bar X_n \xrightarrow{p} \mu$.

▶ **Try it first.** Compute $\mathbb{E}[\bar X_n]$ and $\mathrm{Var}(\bar X_n)$ using (3.14); then Chebyshev.

**Derivation.** $\mathbb{E}[\bar X_n] = \mu$ by linearity. By (3.14) with all $a_i = 1/n$ and all covariances zero (independence), $\mathrm{Var}(\bar X_n) = \frac{1}{n^2}\cdot n\sigma^2 = \sigma^2/n$. Chebyshev: $P(|\bar X_n - \mu| \ge \epsilon) \le \dfrac{\sigma^2}{n\epsilon^2} \to 0$. ∎

■ **What just happened.** Averaging $n$ independent copies divides the variance by $n$; Chebyshev converts small variance into high probability of being close. Independence was used once, to kill the covariances; the proof goes through verbatim for *uncorrelated* variables (exercise B7), which is Markov's point against Nekrasov. The **strong law** says more: $P(\bar X_n \to \mu) = 1$ (almost-sure convergence), and needs only $\mathbb{E}|X_1| < \infty$, no variance. Its proof is beyond this module (Blitzstein §10.2 sketches it; Durrett, *Probability: Theory and Examples*, §2.4 proves it). Operationally: the weak law says at each large $n$ the average is probably near $\mu$; the strong law says the *sequence* of running averages converges, for almost every realization.

**Theorem 3.35 (Central limit theorem).** Let $X_1, X_2, \dots$ be iid with mean $\mu$ and variance $\sigma^2 \in (0, \infty)$. Then
$$\frac{\bar X_n - \mu}{\sigma/\sqrt n} = \frac{1}{\sqrt n}\sum_{i=1}^n\frac{X_i - \mu}{\sigma} \;\xrightarrow{d}\; \mathcal{N}(0, 1). \tag{3.24}$$

We prove it under one extra assumption: the MGF of $X_1$ exists in a neighborhood of 0. (The general finite-variance case replaces MGFs with characteristic functions $\mathbb{E}[e^{itX}]$, which always exist; the structure of the proof is identical.)

▶ **Try it first.** Let $Y_i = (X_i - \mu)/\sigma$, with MGF $M$. Write the MGF of the standardized sum as a power of $M(t/\sqrt n)$. Taylor-expand $M$ to second order around 0 (you know $M(0), M'(0), M''(0)$) and take $n\to\infty$. What limit do you recognize?

**Derivation.**

1. $Y_i$ has mean 0 and variance 1, so its MGF satisfies $M(0) = 1$, $M'(0) = \mathbb{E}[Y] = 0$, $M''(0) = \mathbb{E}[Y^2] = 1$. Since $M$ exists near 0 it is infinitely differentiable there, and Taylor's theorem gives $M(s) = 1 + 0\cdot s + \tfrac12 s^2 + r(s)$ with $r(s)/s^2 \to 0$ as $s \to 0$.
2. Let $S_n = \frac{1}{\sqrt n}\sum Y_i$. By independence and Corollary 3.5, $M_{S_n}(t) = \prod_i\mathbb{E}[e^{tY_i/\sqrt n}] = M(t/\sqrt n)^n = \left(1 + \frac{t^2}{2n} + r\!\left(\tfrac{t}{\sqrt n}\right)\right)^n$.
3. Take logs and use $\log(1 + u) = u + O(u^2)$ for small $u$, with $u_n = \frac{t^2}{2n} + r(t/\sqrt n)$:
$$\log M_{S_n}(t) = n\,u_n + n\,O(u_n^2) = \frac{t^2}{2} + n\,r\!\left(\tfrac{t}{\sqrt n}\right) + O\!\left(\tfrac{1}{n}\right).$$
The middle term: $n\,r(t/\sqrt n) = t^2\cdot\dfrac{r(t/\sqrt n)}{(t/\sqrt n)^2} \to t^2\cdot 0 = 0$ by step 1. So $\log M_{S_n}(t) \to t^2/2$, i.e., $M_{S_n}(t) \to e^{t^2/2}$ for every $t$ near 0.
4. $e^{t^2/2}$ is the MGF of $\mathcal{N}(0,1)$ (F4). The **continuity theorem** (Lévy; for MGFs, Curtiss 1942) states: if $M_{S_n}(t) \to M_S(t)$ for all $t$ in a neighborhood of 0, where $M_S$ is the MGF of some $S$, then $S_n \xrightarrow{d} S$. We cite this; it is not proved here. (Billingsley, *Probability and Measure*, §30; Blitzstein §10.3 states it the same way.) Hence $S_n \xrightarrow{d} \mathcal{N}(0,1)$. ∎

■ **What just happened.** Every distribution with mean 0 and variance 1 has an MGF that looks like $1 + t^2/2$ to second order; the third-order and higher terms are what make distributions different from each other, and raising $M(t/\sqrt n)$ to the $n$-th power wipes them out, because they are $O(n^{-3/2})$ each and there are only $n$ of them. The CLT is the statement that only the first two moments survive aggregation. The honest gap is step 4: we proved the MGFs converge and cited that this forces the distributions to converge.

**How fast?** **Berry–Esseen (1941–42):** if additionally $\rho_3 = \mathbb{E}|X_1 - \mu|^3 < \infty$, then for all $n$ and all $x$,
$$\Big|P\!\left(\tfrac{\bar X_n - \mu}{\sigma/\sqrt n} \le x\right) - \Phi(x)\Big| \le \frac{C\,\rho_3}{\sigma^3\sqrt n}, \qquad C \le 0.4748 \text{ (Shevtsova 2011)}. \tag{3.25}$$
The error decays like $1/\sqrt n$, and the constant depends on the skewness-like ratio $\rho_3/\sigma^3$. For a $\mathrm{Bernoulli}(0.1)$ this ratio is large, which is why that case converges visibly slower in the project.

**What the CLT does not say.**
1. *Heavy tails.* If the variance is infinite, (3.24) is false. The Cauchy distribution has no mean at all, and the average of $n$ iid Cauchys is again Cauchy with the same scale (provable via characteristic functions: $\varphi(t) = e^{-|t|}$, so $\varphi(t/n)^n = e^{-|t|}$): averaging does nothing. A Pareto with tail exponent $\alpha \le 2$ has infinite variance, and its standardized sums converge to a non-Gaussian *stable* law. The project shows the histograms refusing to become bell-shaped.
2. *Dependence.* If $X_i = X_1$ for all $i$ then $\bar X_n = X_1$ forever. Mild dependence (e.g., a Markov chain that mixes fast, §3.12) still gives a CLT, but with a different variance in the limit; strong dependence can destroy it.
3. *The tails at fixed $n$.* (3.25) bounds the *absolute* error uniformly in $x$. A bound of $0.01$ is excellent near the center and worthless at $x = 5$, where $1 - \Phi(5) \approx 3\times 10^{-7}$. The relative error in the tails can be enormous; large-deviation theory (not in this curriculum) handles that regime.

**Normal approximation to the Binomial.** $X \sim \mathrm{Bin}(n, p)$ is a sum of $n$ iid Bernoullis, so $X \approx \mathcal{N}(np, np(1-p))$ for $np$ and $n(1-p)$ both at least about 10. Because $X$ is an integer and the normal is continuous, approximate $P(X \le k)$ by the normal mass up to $k + \tfrac12$ (**continuity correction**):
$$P(X \le k) \approx \Phi\!\left(\frac{k + \tfrac12 - np}{\sqrt{np(1-p)}}\right).$$
For $n = 100$, $p = \tfrac12$, $k = 45$: exact $P(X \le 45) = 0.1841$; without correction $\Phi(-1) = 0.1587$; with correction $\Phi(-0.9) = 0.1841$. The project prints all three.

### 3.12 Markov chains

**Definition 3.36.** A sequence $X_0, X_1, X_2, \dots$ taking values in a finite **state space** $\{1, \dots, m\}$ is a (time-homogeneous) **Markov chain** with **transition matrix** $P$ if
$$P(X_{t+1} = j \mid X_t = i, X_{t-1} = i_{t-1}, \dots, X_0 = i_0) = P(X_{t+1} = j \mid X_t = i) = P_{ij}$$
for all $t$ and all histories. The first equality is the **Markov property**: given the present, the past is irrelevant. $P$ is $m\times m$, has nonnegative entries, and each row sums to 1 ($\sum_j P_{ij} = 1$: from state $i$ you must go somewhere). Such a matrix is called **row-stochastic**.

**Convention (row vectors).** Write the distribution of $X_t$ as a row vector $\boldsymbol\pi_t = (P(X_t = 1), \dots, P(X_t = m))$. LA-7 wrote probability vectors as columns and used a *column*-stochastic matrix acting on the left, $\mathbf{p}_{t+1} = M\mathbf{p}_t$. The two conventions are transposes of each other: $M = P^\top$, $\mathbf{p}_t = \boldsymbol\pi_t^\top$. Probability texts (Blitzstein, Norris, Levin–Peres) use rows; we follow them here. Whenever you read $P_{ij}$, say "from $i$ to $j$."

**Theorem 3.37 (Evolution).** $\boldsymbol\pi_{t+1} = \boldsymbol\pi_t P$, hence $\boldsymbol\pi_t = \boldsymbol\pi_0 P^t$.

*Derivation.* Law of total probability, conditioning on $X_t$: $P(X_{t+1} = j) = \sum_i P(X_t = i)P(X_{t+1} = j \mid X_t = i) = \sum_i (\boldsymbol\pi_t)_i P_{ij} = (\boldsymbol\pi_t P)_j$. Iterate. ∎

**Chapman–Kolmogorov.** Let $P^{(n)}_{ij} = P(X_{t+n} = j \mid X_t = i)$ be the $n$-step transition probability. Then $P^{(m+n)}_{ij} = \sum_k P^{(m)}_{ik}P^{(n)}_{kj}$, i.e., $P^{(m+n)} = P^{(m)}P^{(n)}$, and so $P^{(n)} = P^n$. *Derivation:* condition on the state at the intermediate time $t + m$, use the Markov property to drop the earlier history. ∎ The $n$-step probabilities are the entries of the matrix power: LA-7's $P^n$ via diagonalization is exactly what you need.

**Classification (brief).** State $j$ is **accessible** from $i$ if $P^n_{ij} > 0$ for some $n$. The chain is **irreducible** if every state is accessible from every other. The **period** of state $i$ is $\gcd\{n \ge 1 : P^n_{ii} > 0\}$; the chain is **aperiodic** if every period is 1 (a single $P_{ii} > 0$ anywhere in an irreducible chain suffices). A state is **recurrent** if, starting there, the chain returns with probability 1, **transient** otherwise; in a finite irreducible chain every state is recurrent (stated; Blitzstein §11.2). A bipartite random walk, which alternates between two sets of states, has period 2 (exercise B8).

**Definition 3.38.** A **stationary distribution** is a row vector $\boldsymbol\pi$ with $\pi_i \ge 0$, $\sum_i\pi_i = 1$, and
$$\boldsymbol\pi P = \boldsymbol\pi. \tag{3.26}$$
If $\boldsymbol\pi_0 = \boldsymbol\pi$ then $\boldsymbol\pi_t = \boldsymbol\pi$ forever. Equation (3.26) says $\boldsymbol\pi$ is a **left eigenvector of $P$ with eigenvalue 1**; transposing, $P^\top\boldsymbol\pi^\top = \boldsymbol\pi^\top$, an ordinary eigenvector of $P^\top$, which is how you compute it (LA-7 and the project).

**Theorem 3.39 (Existence for finite chains).** Every finite row-stochastic $P$ has a stationary distribution.

*Derivation (with one citation).* (i) $P\mathbf{1} = \mathbf{1}$ because rows sum to 1, so 1 is an eigenvalue of $P$ with right eigenvector $\mathbf{1}$. $P$ and $P^\top$ have the same eigenvalues (same characteristic polynomial, LA-6), so there is a nonzero row vector $\mathbf{v}$ with $\mathbf{v}P = \mathbf{v}$. (ii) No eigenvalue exceeds 1 in modulus: if $P\mathbf{x} = \lambda\mathbf{x}$, pick the index $i$ with $|x_i|$ largest; then $|\lambda||x_i| = |\sum_j P_{ij}x_j| \le \sum_j P_{ij}|x_j| \le |x_i|$, so $|\lambda| \le 1$. (iii) What remains is that $\mathbf{v}$ can be chosen with nonnegative entries. This is the **Perron–Frobenius theorem** (stated, not proved; Horn & Johnson ch. 8): a nonnegative matrix has its spectral radius as an eigenvalue with a nonnegative eigenvector, and if the matrix is irreducible that eigenvalue is simple and the eigenvector is strictly positive. Normalize $\mathbf{v}$ to sum to 1. ∎

**Theorem 3.40 (Convergence).** If $P$ is irreducible and aperiodic, the stationary distribution $\boldsymbol\pi$ is unique, $\pi_i > 0$ for all $i$, and $\boldsymbol\pi_0P^t \to \boldsymbol\pi$ for every initial distribution $\boldsymbol\pi_0$. Moreover the long-run fraction of time spent in state $i$ converges to $\pi_i$ (the ergodic theorem for Markov chains, which is the law of large numbers that Markov built).

*The two-state case, proved completely.* Let $P = \begin{pmatrix}1-a & a\\ b & 1-b\end{pmatrix}$ with $0 < a, b < 1$ (irreducible since $a, b > 0$; aperiodic since the diagonal is positive). Write $\boldsymbol\pi_t = (p_t, 1 - p_t)$.

1. Evolution: $p_{t+1} = p_t(1 - a) + (1 - p_t)b = b + (1 - a - b)p_t$.
2. Fixed point: $p^* = b + (1 - a - b)p^*$ gives $p^* = \dfrac{b}{a + b}$, so $\boldsymbol\pi = \left(\dfrac{b}{a+b}, \dfrac{a}{a+b}\right)$. Check (3.26): $\frac{b(1-a) + ab}{a+b} = \frac{b}{a+b}$. ✓
3. Subtract: $p_{t+1} - p^* = (1 - a - b)(p_t - p^*)$, so $p_t - p^* = (1 - a - b)^t(p_0 - p^*)$.
4. Since $0 < a, b < 1$, we have $-1 < 1 - a - b < 1$, so $(1-a-b)^t \to 0$ and $p_t \to p^*$ from any $p_0$. The rate is geometric with ratio $|\lambda_2|$, where $\lambda_2 = 1 - a - b$ is the second eigenvalue of $P$ (trace $= 2 - a - b = 1 + \lambda_2$). ∎

*General case, sketched.* If $P$ is diagonalizable, $\boldsymbol\pi_0P^t = \sum_i c_i\lambda_i^t\mathbf{v}_i$ in the basis of left eigenvectors (LA-7). Irreducible + aperiodic is exactly the hypothesis under which Perron–Frobenius guarantees that $\lambda_1 = 1$ is simple and every other eigenvalue has $|\lambda_i| < 1$, so all terms but the first die and $\boldsymbol\pi_0P^t \to c_1\mathbf{v}_1 = \boldsymbol\pi$ (the constant $c_1 = 1$ because probabilities sum to 1 at every $t$ and $\mathbf{v}_1 = \boldsymbol\pi$ sums to 1). The error decays like $|\lambda_2|^t$; the **mixing time** is of order $1/\log(1/|\lambda_2|)$. The project measures this slope. For the non-diagonalizable case and the full proof, Levin & Peres, *Markov Chains and Mixing Times*, ch. 4.

**Theorem 3.41 (Detailed balance ⇒ stationary).** If a distribution $\boldsymbol\pi$ satisfies
$$\pi_iP_{ij} = \pi_jP_{ji} \quad\text{for all } i, j, \tag{3.27}$$
then $\boldsymbol\pi$ is stationary. A chain with such a $\boldsymbol\pi$ is called **reversible**.

*Derivation.* $(\boldsymbol\pi P)_j = \sum_i\pi_iP_{ij} = \sum_i\pi_jP_{ji} = \pi_j\sum_iP_{ji} = \pi_j$, using (3.27) and then that row $j$ of $P$ sums to 1. ∎

■ **What just happened.** Detailed balance says the probability flow $i \to j$ equals the flow $j \to i$ for *every pair*. Stationarity only says total inflow equals total outflow at every state. Detailed balance is stronger, and false for some stationary chains (a chain that cycles $1\to2\to3\to1$ with high probability has uniform stationary distribution but net circulation; the project checks). Its power: (3.27) involves only *ratios* $\pi_j/\pi_i$. If you want to sample from a $\boldsymbol\pi$ known only up to a normalizing constant (a Bayesian posterior in ST-6, a Boltzmann distribution in SM-6), you can design $P$ to satisfy (3.27) without ever computing the constant. That is the Metropolis algorithm: propose a move, accept it with probability $\min(1, \pi_j/\pi_i)$, and (3.27) holds by construction (exercise B10 verifies). This theorem is the engine of all of Markov chain Monte Carlo.

**Random walk on a graph.** Let $G$ be a connected undirected graph with vertex degrees $d_i$ and $|E|$ edges. The random walk moves from $i$ to a uniformly random neighbor: $P_{ij} = 1/d_i$ if $i \sim j$, else 0.

▶ **Try it first.** Guess $\boldsymbol\pi$ (where does a random walker spend more time: a vertex with many roads or few?) and check detailed balance.

**Derivation.** Try $\pi_i = d_i/(2|E|)$ (this sums to 1 since $\sum_i d_i = 2|E|$, each edge counted from both ends). For $i \sim j$: $\pi_iP_{ij} = \frac{d_i}{2|E|}\cdot\frac{1}{d_i} = \frac{1}{2|E|}$, and by the same computation $\pi_jP_{ji} = \frac{1}{2|E|}$. Equal. For $i \not\sim j$ both sides are 0. So (3.27) holds and by Theorem 3.41, $\boldsymbol\pi$ is stationary: **the stationary distribution is proportional to degree.** ∎ (If $G$ is bipartite the walk has period 2 and $\boldsymbol\pi_t$ does not converge, though $\boldsymbol\pi$ is still stationary; the lazy walk that stays put with probability $\tfrac12$ fixes this. Exercise B8.)

**PageRank, in one paragraph.** The web is a directed graph; the random surfer follows a random out-link with probability $\alpha \approx 0.85$ and jumps to a uniformly random page with probability $1 - \alpha$. The transition matrix is $G = \alpha P + (1 - \alpha)\tfrac{1}{m}\mathbf{1}\mathbf{1}^\top$ (with dangling pages patched to jump always). The teleportation makes $G$ strictly positive, hence irreducible and aperiodic, so by Theorem 3.40 there is a unique stationary $\boldsymbol\pi$ and $\boldsymbol\pi_0G^t \to \boldsymbol\pi$ from any start. That $\boldsymbol\pi$ is the PageRank vector; it is computed by power iteration (LA-7), and the second eigenvalue of $G$ is at most $\alpha$, so about $\log(10^{-8})/\log(0.85) \approx 113$ iterations suffice for eight digits. A page's rank is the long-run fraction of time the surfer spends there: high if many high-rank pages link to it.

### 3.13 Information theory: entropy and KL divergence

**Surprise.** The **surprise** (or information content) of an outcome with probability $p$ is $-\log p$. Shannon's requirements: surprise should be a decreasing function of $p$; a certain event ($p = 1$) should carry none; and the surprise of two *independent* events should add. Since the probability of two independent events multiplies, additivity forces a logarithm: $s(pq) = s(p) + s(q)$ with $s$ continuous and decreasing has only the solutions $s(p) = -c\log p$ with $c > 0$ (stated; it is the Cauchy functional equation). The constant $c$ is the choice of base: $\log_2$ gives **bits**, $\ln$ gives **nats**; $1\text{ nat} = 1/\ln 2 \approx 1.443$ bits.

**Definition 3.42 (Entropy).** For a discrete $X$ with PMF $p$,
$$H(X) = \mathbb{E}[-\log p(X)] = -\sum_x p(x)\log p(x), \tag{3.28}$$
with the convention $0\log 0 = 0$ (justified by $\lim_{p\to0^+}p\log p = 0$). $H$ depends only on the multiset of probabilities, not on the labels or values of $X$.

**Theorem 3.43.** (a) $H(X) \ge 0$, with equality iff $X$ is constant. (b) If $X$ takes at most $m$ values, $H(X) \le \log m$, with equality iff $X$ is uniform on $m$ values.

*Derivation.* (a) Each term $-p(x)\log p(x) \ge 0$ since $0 < p \le 1$; the sum is 0 iff every $p(x) \in \{0, 1\}$. (b) Jensen (Theorem 3.32) with the concave function $\log$ applied to the random variable $1/p(X)$:
$$H(X) = \mathbb{E}\!\left[\log\frac{1}{p(X)}\right] \le \log\mathbb{E}\!\left[\frac{1}{p(X)}\right] = \log\sum_{x : p(x) > 0}p(x)\cdot\frac{1}{p(x)} = \log(\#\text{support}) \le \log m.$$
Equality in Jensen (strictly concave $\log$) iff $1/p(X)$ is constant iff $p$ is uniform on its support; equality in the last step iff the support has all $m$ values. ∎ (Exercise C3 gives the Lagrange-multiplier proof instead.)

**Entropy of a coin.** $h(p) := H(\mathrm{Bernoulli}(p)) = -p\log_2p - (1-p)\log_2(1-p)$ bits. In words: zero at $p = 0$ and $p = 1$ (no surprise), symmetric about $p = \tfrac12$, maximum of exactly 1 bit at $p = \tfrac12$, concave, and with *vertical* tangents at the endpoints ($h'(p) = \log_2\frac{1-p}{p} \to \pm\infty$). At $p = 0.1$, $h \approx 0.469$ bits: a 90/10 coin carries less than half the information of a fair one. The project plots it.

**Joint and conditional entropy.** $H(X, Y) = -\sum_{x,y}p(x,y)\log p(x,y)$ is the entropy of the pair. The **conditional entropy** is the average, over $x$, of the entropy of the slice:
$$H(Y \mid X) = \sum_xp_X(x)H(Y \mid X = x) = -\sum_{x,y}p(x,y)\log p(y\mid x).$$
**Chain rule:** $H(X, Y) = H(X) + H(Y\mid X)$. *Derivation:* $\log p(x,y) = \log p_X(x) + \log p(y\mid x)$ by the chain rule for probabilities; multiply by $-p(x,y)$ and sum; the first term gives $-\sum_xp_X(x)\log p_X(x)$ by (3.1). ∎

**Definition 3.44 (KL divergence / relative entropy).** For PMFs $p, q$ on the same set,
$$D(p\,\|\,q) = \sum_xp(x)\log\frac{p(x)}{q(x)} = \mathbb{E}_p\!\left[\log\frac{p(X)}{q(X)}\right], \tag{3.29}$$
with $0\log(0/q) = 0$ and $D = +\infty$ if some $p(x) > 0$ where $q(x) = 0$. In words: the expected *extra* surprise when you use $q$ to code outcomes that actually follow $p$.

**Theorem 3.45 (Gibbs' inequality).** $D(p\,\|\,q) \ge 0$, with equality if and only if $p = q$.

▶ **Try it first.** Write $-D$ as $\mathbb{E}_p[\log(q/p)]$ and apply Jensen with the concave $\log$. Watch what $\mathbb{E}_p[q(X)/p(X)]$ sums to.

**Derivation.** Restrict to the support $A = \{x : p(x) > 0\}$ (terms outside contribute 0 to $D$; if $q$ vanishes somewhere on $A$ then $D = \infty > 0$ and we are done). Jensen for concave $\log$:
$$-D(p\|q) = \sum_{x\in A}p(x)\log\frac{q(x)}{p(x)} = \mathbb{E}_p\!\left[\log\frac{q(X)}{p(X)}\right] \le \log\mathbb{E}_p\!\left[\frac{q(X)}{p(X)}\right] = \log\sum_{x\in A}q(x) \le \log 1 = 0.$$
Equality requires equality in Jensen (so $q/p$ is constant on $A$, since $\log$ is strictly concave) *and* $\sum_{x\in A}q(x) = 1$ (so $q$ has no mass outside $A$). Constant ratio with both summing to 1 forces $q = p$. ∎

■ **What just happened.** One application of Jensen gives the inequality on which maximum likelihood (ST-1), maximum entropy (SM-4), the EM algorithm (ML-6) and variational inference (ST-6, ML-10) all rest. Note what KL is **not**: it is not symmetric ($D(p\|q) \ne D(q\|p)$ in general; Example 4.7 has $0.511$ nats one way and $0.368$ the other) and it does not satisfy the triangle inequality, so it is not a distance. It is a *divergence*: zero iff equal, positive otherwise, and that is all.

**Cross-entropy.** $H(p, q) := -\sum_xp(x)\log q(x) = \mathbb{E}_p[-\log q(X)]$, the average surprise when the world is $p$ and you believe $q$. Adding and subtracting $\log p$:
$$H(p, q) = -\sum_xp\log p + \sum_xp\log\frac pq = H(p) + D(p\,\|\,q). \tag{3.30}$$
So $H(p, q) \ge H(p)$: believing the wrong distribution always costs. In ML-3, $p$ is the data distribution (fixed) and $q_\theta$ is the model; minimizing cross-entropy over $\theta$ is minimizing $D(p\|q_\theta)$ because $H(p)$ is a constant. And with $p$ the empirical distribution of a training set, $H(p, q_\theta) = -\frac1n\sum_i\log q_\theta(x_i)$ is the negative average log-likelihood: **cross-entropy minimization is maximum likelihood** (ST-1), and the gap to the minimum possible loss is exactly $D(p\|q_\theta)$.

**Mutual information.** $I(X; Y) := D(p_{X,Y}\,\|\,p_Xp_Y) = \sum_{x,y}p(x,y)\log\dfrac{p(x,y)}{p_X(x)p_Y(y)}$. By Gibbs, $I(X;Y) \ge 0$ with equality iff $p_{X,Y} = p_Xp_Y$, i.e., iff $X \perp Y$ (Theorem 3.4). Expanding the log: $I(X;Y) = H(X) + H(Y) - H(X,Y) = H(Y) - H(Y\mid X)$ (chain rule). The second form says: mutual information is the reduction in uncertainty about $Y$ from learning $X$, and it is never negative: **conditioning reduces entropy on average** (exercise B9). Unlike covariance, $I$ detects every kind of dependence, including $Y = X^2$.

**Differential entropy.** For a density $f$, $h(X) := -\int f(x)\log f(x)\,dx$. Two warnings. It can be **negative**: $\mathrm{Uniform}(0, \tfrac12)$ has $f = 2$ and $h = -\log 2$. And it is not scale-invariant: $h(aX) = h(X) + \log|a|$ (change of variables). The reason is that $h$ is really "entropy relative to the Lebesgue measure $dx$," and $dx$ is not a probability distribution; differences and KL divergences between densities behave perfectly, and KL is scale-invariant. Do not read a negative $h$ as a bug (§8).

**Entropy of a Gaussian.** For $X \sim \mathcal{N}(\mu, \sigma^2)$, $\log f(x) = -\tfrac12\log(2\pi\sigma^2) - \frac{(x-\mu)^2}{2\sigma^2}$, so
$$h(X) = -\mathbb{E}[\log f(X)] = \tfrac12\log(2\pi\sigma^2) + \frac{\mathbb{E}[(X-\mu)^2]}{2\sigma^2} = \tfrac12\log(2\pi\sigma^2) + \tfrac12 = \tfrac12\log(2\pi e\,\sigma^2). \tag{3.31}$$
For $\sigma = 1$: $1.4189$ nats $= 2.047$ bits. For $\mathbf{X}\sim\mathcal{N}(\boldsymbol\mu, \Sigma)$ in $\mathbb{R}^k$, the same computation with $\mathbb{E}[(\mathbf{X}-\boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{X}-\boldsymbol\mu)] = \mathbb{E}[\|\mathbf{Z}\|^2] = k$ gives $h(\mathbf{X}) = \tfrac12\log\big((2\pi e)^k\det\Sigma\big)$.

**Theorem 3.46 (The Gaussian maximizes entropy for fixed variance).** Among all densities $f$ on $\mathbb{R}$ with mean $\mu$ and variance $\sigma^2$, $h(f) \le \tfrac12\log(2\pi e\sigma^2)$, with equality iff $f$ is $\mathcal{N}(\mu, \sigma^2)$.

▶ **Try it first.** Let $g$ be the $\mathcal{N}(\mu,\sigma^2)$ density. Expand $0 \le D(f\|g) = -h(f) - \int f\log g$. Compute $-\int f\log g$: notice it involves only the mean and variance of $f$.

**Derivation.** Gibbs' inequality holds for densities with the same proof (replace sums by integrals). So
$$0 \le D(f\,\|\,g) = \int f\log f - \int f\log g = -h(f) - \int f(x)\Big[-\tfrac12\log(2\pi\sigma^2) - \frac{(x-\mu)^2}{2\sigma^2}\Big]dx = -h(f) + \tfrac12\log(2\pi\sigma^2) + \frac{\sigma^2}{2\sigma^2},$$
because $\int f(x)(x-\mu)^2dx = \sigma^2$ by the assumption on $f$. The right side is $-h(f) + h(g)$ by (3.31). Hence $h(f) \le h(g)$, with equality iff $D(f\|g) = 0$ iff $f = g$. ∎

■ **What just happened.** $-\int f\log g$ depends on $f$ only through the constraints (mean and variance) because $\log g$ is a *quadratic*. That is the general principle: the maximum-entropy density under constraints $\mathbb{E}[T_j(X)] = c_j$ has $\log f$ linear in the $T_j$. Variance constraint ⇒ Gaussian. Mean constraint on $[0,\infty)$ ⇒ exponential (next). Energy constraint ⇒ Boltzmann.

**The exponential maximizes entropy on $[0, \infty)$ for fixed mean.** Maximize $h(f) = -\int_0^\infty f\log f$ subject to $\int_0^\infty f = 1$ and $\int_0^\infty xf = \mu$. Lagrangian (F6, calculus of variations, with multipliers $\alpha, \beta$):
$$\mathcal{L}[f] = -\int f\log f\,dx - \alpha\Big(\int f\,dx - 1\Big) - \beta\Big(\int xf\,dx - \mu\Big).$$
Varying $f$ pointwise, the integrand's derivative with respect to $f(x)$ must vanish: $-\log f(x) - 1 - \alpha - \beta x = 0$, so $f(x) = e^{-1-\alpha}e^{-\beta x} \propto e^{-\beta x}$. Normalizing on $[0,\infty)$ and matching the mean gives $f(x) = \frac1\mu e^{-x/\mu}$, the $\mathrm{Exponential}(1/\mu)$ density, with $h = 1 + \log\mu$. (The KL argument of Theorem 3.46 with $g$ exponential confirms it is the maximum, not a saddle: $-\int f\log g = \log\mu + \mathbb{E}_f[X]/\mu = \log\mu + 1 = h(g)$.) **SM-4 is this calculation** with $x$ replaced by the energy $E_i$ of microstate $i$ and the sum over microstates: the maximum-entropy distribution with fixed mean energy is $p_i \propto e^{-\beta E_i}$, the Boltzmann distribution, and $\beta = 1/k_BT$.

**KL between two Gaussians (1D).** Let $f_1 = \mathcal{N}(\mu_1, \sigma_1^2)$, $f_2 = \mathcal{N}(\mu_2, \sigma_2^2)$. Then $D(f_1\|f_2) = \mathbb{E}_1[\log f_1] - \mathbb{E}_1[\log f_2] = -h(f_1) - \mathbb{E}_1[\log f_2]$. Using $\mathbb{E}_1[(X - \mu_2)^2] = \sigma_1^2 + (\mu_1 - \mu_2)^2$ (expand around $\mu_1$):
$$-\mathbb{E}_1[\log f_2] = \tfrac12\log(2\pi\sigma_2^2) + \frac{\sigma_1^2 + (\mu_1-\mu_2)^2}{2\sigma_2^2}, \qquad D(f_1\,\|\,f_2) = \log\frac{\sigma_2}{\sigma_1} + \frac{\sigma_1^2 + (\mu_1 - \mu_2)^2}{2\sigma_2^2} - \frac12. \tag{3.32}$$
Check: $\mu_1 = \mu_2$, $\sigma_1 = \sigma_2$ gives $0 + \tfrac12 - \tfrac12 = 0$. ✓ The multivariate version, by the same route, is $D = \tfrac12\big[\mathrm{tr}(\Sigma_2^{-1}\Sigma_1) + (\boldsymbol\mu_2 - \boldsymbol\mu_1)^\top\Sigma_2^{-1}(\boldsymbol\mu_2 - \boldsymbol\mu_1) - k + \log\frac{\det\Sigma_2}{\det\Sigma_1}\big]$; it is the KL term in every VAE loss (ML-10).

## 4. Worked examples

### Example 4.1: A joint PMF table, end to end

$X \in \{0, 1, 2\}$, $Y \in \{0, 1\}$, joint PMF:

| $p(x,y)$ | $y = 0$ | $y = 1$ | $p_X(x)$ |
|---|---|---|---|
| $x = 0$ | 0.10 | 0.20 | **0.30** |
| $x = 1$ | 0.20 | 0.20 | **0.40** |
| $x = 2$ | 0.20 | 0.10 | **0.30** |
| $p_Y(y)$ | **0.50** | **0.50** | 1 |

*Marginals* (3.1): row sums and column sums, shown in bold.

*Conditionals* (3.3): $p_{Y|X}(\cdot\mid 0) = (0.10, 0.20)/0.30 = (\tfrac13, \tfrac23)$; $p_{Y|X}(\cdot\mid 1) = (\tfrac12, \tfrac12)$; $p_{Y|X}(\cdot\mid 2) = (\tfrac23, \tfrac13)$. The slices have different shapes, so $X$ and $Y$ are dependent. Directly: $p(0,0) = 0.10 \ne p_X(0)p_Y(0) = 0.15$. One failing cell is enough.

*Moments.* $\mathbb{E}[X] = 0(0.3) + 1(0.4) + 2(0.3) = 1.0$. $\mathbb{E}[X^2] = 0.4 + 4(0.3) = 1.6$, so $\mathrm{Var}(X) = 0.6$. $\mathbb{E}[Y] = 0.5$, $\mathrm{Var}(Y) = 0.25$. $\mathbb{E}[XY] = \sum xy\,p(x,y) = 1\cdot1\cdot0.20 + 2\cdot1\cdot0.10 = 0.40$ (only cells with $x, y$ both nonzero contribute). $\mathrm{Cov}(X,Y) = 0.40 - (1.0)(0.5) = -0.10$. $\rho = -0.10/\sqrt{0.6\cdot0.25} = -0.10/0.3873 = -0.258$. Negative: larger $X$ makes $Y = 1$ less likely, as the conditionals showed.

*Conditional expectation.* $\mathbb{E}[Y\mid X = x] = p_{Y|X}(1\mid x)$: $(\tfrac23, \tfrac12, \tfrac13)$ for $x = 0, 1, 2$. Adam's law check: $0.3(\tfrac23) + 0.4(\tfrac12) + 0.3(\tfrac13) = 0.2 + 0.2 + 0.1 = 0.5 = \mathbb{E}[Y]$. ✓ Eve's law: $\mathrm{Var}(Y\mid X = x) = p(1-p)$ for $p = \tfrac23, \tfrac12, \tfrac13$: $(\tfrac29, \tfrac14, \tfrac29)$, average $0.3(\tfrac29) + 0.4(\tfrac14) + 0.3(\tfrac29) = 0.1333 + 0.1 = 0.2333$. $\mathrm{Var}(\mathbb{E}[Y\mid X]) = \mathbb{E}[g(X)^2] - 0.25 = 0.3(\tfrac49) + 0.4(\tfrac14) + 0.3(\tfrac19) - 0.25 = 0.1333 + 0.1 + 0.0333 - 0.25 = 0.0167$. Sum: $0.2333 + 0.0167 = 0.25 = \mathrm{Var}(Y)$. ✓ Only 6.7% of the variance of $Y$ is "explained" by $X$; the best linear predictor (3.23) has slope $-0.10/0.6 = -0.167$. The project reproduces every number.

### Example 4.2: Box–Muller, with numbers

The derivation is §3.3. Here is what the computer does. Given $u_1 = 0.5$, $u_2 = 0.25$: $R = \sqrt{-2\ln 0.5} = \sqrt{1.3863} = 1.1774$, $\Theta = 2\pi(0.25) = \pi/2$, so $(Z_1, Z_2) = (1.1774\cos\tfrac\pi2, 1.1774\sin\tfrac\pi2) = (0, 1.1774)$. Given $u_1 = 0.1$, $u_2 = 0.125$: $R = \sqrt{-2\ln0.1} = \sqrt{4.6052} = 2.1460$, $\Theta = \pi/4$, $(Z_1, Z_2) = (1.5174, 1.5174)$. Check the inverse map: $z_1^2 + z_2^2 = 2(1.5174)^2 = 4.605 = -2\ln 0.1$. ✓

Why this and not "add twelve uniforms and subtract 6" (which has mean 0 and variance 1)? Because the sum of twelve uniforms is *approximately* normal by the CLT, with a hard cutoff at $\pm6$ and the wrong kurtosis; Box–Muller is *exactly* normal. The project draws $10^5$ pairs, compares the quantiles against `scipy.stats.norm.ppf` in a QQ plot (a straight line through the origin with slope 1), checks that $Z_1$ and $Z_2$ have sample correlation within $\pm 0.01$ of zero, and runs a Kolmogorov–Smirnov test.

### Example 4.3: Sum of two independent Uniforms, by convolution

$X, Y \sim \mathrm{Uniform}(0,1)$ independent, $S = X + Y$. By (3.7), $f_S(s) = \int f_X(x)f_Y(s - x)\,dx$. The integrand is 1 when both $0 \le x \le 1$ and $0 \le s - x \le 1$ (i.e., $s - 1 \le x \le s$), else 0. So $f_S(s)$ is the length of the interval $[\max(0, s-1), \min(1, s)]$:
$$f_S(s) = \begin{cases} s & 0 \le s \le 1 \\ 2 - s & 1 \le s \le 2 \\ 0 & \text{otherwise.}\end{cases}$$
A triangle. Checks: area $= 1$; $\mathbb{E}[S] = 1 = \tfrac12 + \tfrac12$; $\mathrm{Var}(S) = \tfrac1{12} + \tfrac1{12} = \tfrac16$ (Corollary 3.15 with zero covariance); $P(S \le \tfrac12) = \int_0^{1/2}s\,ds = \tfrac18 = 0.125$. The sum of three uniforms is piecewise quadratic (the Irwin–Hall distribution) and already visibly bell-shaped: the CLT begins at $n = 3$ for this distribution, which is why the $n = 2, 5$ panels in the project look so different for Uniform versus Exponential.

### Example 4.4: The bivariate normal with $\Sigma = \begin{pmatrix}2 & 1\\1 & 2\end{pmatrix}$

Take $\boldsymbol\mu = \mathbf{0}$. *Eigenstructure* (LA-8): $\det(\Sigma - \lambda I) = (2-\lambda)^2 - 1 = 0$ gives $\lambda_1 = 3$, $\lambda_2 = 1$, with eigenvectors $\mathbf{q}_1 = (1, 1)/\sqrt2$ and $\mathbf{q}_2 = (1, -1)/\sqrt2$. So the density's level curves are ellipses with long axis along the diagonal (semi-axis $\propto\sqrt3 = 1.732$) and short axis along the anti-diagonal (semi-axis $\propto 1$). $\det\Sigma = 3$, $\Sigma^{-1} = \tfrac13\begin{pmatrix}2 & -1\\ -1 & 2\end{pmatrix}$, and the peak density is $f(\mathbf 0) = 1/(2\pi\sqrt3) = 0.0919$.

*Correlation:* $\rho = 1/\sqrt{2\cdot2} = 0.5$.

*The 95% ellipse:* $(\mathbf{x}^\top\Sigma^{-1}\mathbf{x} \le \chi^2_{2,0.95} = 5.991)$ has semi-axes $\sqrt{5.991\cdot 3} = 4.24$ along $\mathbf{q}_1$ and $\sqrt{5.991} = 2.45$ along $\mathbf{q}_2$.

*Mahalanobis distance:* for $\mathbf{x} = (1, 1)$: $d^2 = \tfrac13(2 - 1 - 1 + 2) = \tfrac23$. For $\mathbf{x} = (1, -1)$: $d^2 = \tfrac13(2 + 1 + 1 + 2) = 2$. Same Euclidean length $\sqrt2$, but $(1,-1)$ is three times as "far" in Mahalanobis terms, because it points along the short axis. $P(d^2 > 2) = e^{-1} = 0.368$ while $P(d^2 > \tfrac23) = e^{-1/3} = 0.717$ (since $d^2 \sim \chi^2_2 = \mathrm{Exponential}(\tfrac12)$).

*Conditional of $X_1$ given $X_2 = 1$* by (3.19), with $\Sigma_{11} = 2$, $\Sigma_{12} = 1$, $\Sigma_{22} = 2$: mean $= 0 + 1\cdot\tfrac12\cdot(1 - 0) = 0.5$; variance $= 2 - 1\cdot\tfrac12\cdot1 = 1.5$. So $X_1\mid X_2 = 1 \sim \mathcal{N}(0.5, 1.5)$. Sanity: the conditional variance $1.5 = \sigma_1^2(1 - \rho^2) = 2(0.75)$ ✓, and the regression slope $\Sigma_{12}/\Sigma_{22} = 0.5 = \rho\sigma_1/\sigma_2$ ✓. The project slices $4\times10^5$ Cholesky samples at $|x_2 - 1| < 0.05$ and asserts the sample mean and variance of the slice match $0.5$ and $1.5$.

*Cholesky factor:* $L = \begin{pmatrix}\sqrt2 & 0\\ 1/\sqrt2 & \sqrt{3/2}\end{pmatrix} = \begin{pmatrix}1.4142 & 0 \\ 0.7071 & 1.2247\end{pmatrix}$. Check $LL^\top = \begin{pmatrix}2 & 1\\1 & \tfrac12 + \tfrac32\end{pmatrix} = \Sigma$. ✓ To sample: $\mathbf{x} = L\mathbf{z}$, i.e., $x_1 = \sqrt2 z_1$, $x_2 = z_1/\sqrt2 + \sqrt{3/2}\,z_2$.

### Example 4.5: Galton's heights (a statistics example)

Model mid-parent height $P$ and adult child height $C$ as jointly Gaussian, both with mean $\mu = 68$ inches and standard deviation $\sigma = 2.5$ inches, correlation $\rho = 0.5$. (Galton's 1886 table of 928 children had the mid-parent SD smaller than the children's, because averaging two parents reduces variance, and his fitted slope was about $\tfrac23$. With equal SDs the slope *is* $\rho$; we use round numbers.)

*The regression line* (3.19)/(3.23): $\mathbb{E}[C\mid P = p] = 68 + 0.5\,(p - 68)$. A mid-parent of $73$ inches (5 above the mean, $z = 2$) has children expected at $68 + 2.5 = 70.5$ inches ($z = 1$). Half the excess is gone: regression to the mean.

*Conditional spread:* $\mathrm{Var}(C\mid P) = \sigma^2(1 - \rho^2) = 6.25\cdot0.75 = 4.6875$, SD $= 2.165$ inches, the same for every $p$. So $P(C > 73 \mid P = 73) = 1 - \Phi\big(\tfrac{73 - 70.5}{2.165}\big) = 1 - \Phi(1.155) = 0.124$: only about one child in eight out-grows a 73-inch mid-parent.

*It is symmetric:* $\mathbb{E}[P\mid C = 73] = 68 + 0.5(73 - 68) = 70.5$ as well. Exceptionally tall children have, on average, less exceptional parents. If "regression to mediocrity" were a biological force acting forward in time, it could not also act backward. It is not a force; it is $\rho < 1$.

*Eve's law:* $\mathrm{Var}(C) = \mathbb{E}[\mathrm{Var}(C\mid P)] + \mathrm{Var}(\mathbb{E}[C\mid P]) = 4.6875 + 0.25\cdot6.25 = 4.6875 + 1.5625 = 6.25$. ✓ The explained fraction is $\rho^2 = 0.25$: knowing the parents' height removes a quarter of the variance in the child's. The project simulates $10^5$ families and recovers the slope, the intercept, and the conditional mean of children of 73-inch parents.

### Example 4.6: A three-state weather chain

States (Sunny, Cloudy, Rainy), transition matrix
$$P = \begin{pmatrix} 0.7 & 0.3 & 0 \\ 0.3 & 0.4 & 0.3 \\ 0 & 0.6 & 0.4\end{pmatrix}.$$
Rows sum to 1. Irreducible ($S\leftrightarrow C\leftrightarrow R$) and aperiodic ($P_{SS} > 0$), so Theorem 3.40 applies.

*Via the eigenvector.* Solve $\boldsymbol\pi P = \boldsymbol\pi$, i.e., $P^\top\boldsymbol\pi^\top = \boldsymbol\pi^\top$, for $\boldsymbol\pi = (s, c, r)$:
- column S: $0.7s + 0.3c = s \Rightarrow c = s$;
- column R: $0.3c + 0.4r = r \Rightarrow r = c/2$;
- column C is then automatic (the three equations are dependent, as they must be since $\det(P - I) = 0$): $0.3s + 0.4c + 0.6r = 0.3s + 0.4s + 0.3s = s = c$. ✓
Normalize $s + s + s/2 = 1$: $\boldsymbol\pi = (0.4, 0.4, 0.2)$.

*Via detailed balance.* $S$ and $R$ never communicate directly, so the chain is a path, and (3.27) along the two edges gives $\pi_S(0.3) = \pi_C(0.3)$ and $\pi_C(0.3) = \pi_R(0.6)$: $\pi_S = \pi_C = 2\pi_R$. Same answer, with no linear algebra. Any chain whose transition graph is a tree is reversible for this reason.

*Mixing.* $\mathrm{tr}P = 1.5 = 1 + \lambda_2 + \lambda_3$ and $\det P = 0.7(0.16 - 0.18) - 0.3(0.12) = -0.05 = \lambda_2\lambda_3$. So $\lambda_{2,3}$ solve $\lambda^2 - 0.5\lambda - 0.05 = 0$: $\lambda_2 = 0.5854$, $\lambda_3 = -0.0854$. Starting from a sunny day, $\boldsymbol\pi_1 = (0.7, 0.3, 0)$ and $\boldsymbol\pi_2 = (0.58, 0.33, 0.09)$; the $\ell^1$ distances to $\boldsymbol\pi$ are $0.6$ and $0.36$, ratio $0.6 \approx \lambda_2$. The project plots $\|\boldsymbol\pi_t - \boldsymbol\pi\|_1$ on a log scale and fits the slope to $\log\lambda_2 = -0.535$.

### Example 4.7: KL divergences, both ways

*Bernoulli.* $p = \mathrm{Bern}(0.5)$, $q = \mathrm{Bern}(0.9)$.
$$D(p\|q) = 0.5\ln\frac{0.5}{0.9} + 0.5\ln\frac{0.5}{0.1} = 0.5\ln\frac{25}{9} = 0.5108\text{ nats} = 0.7370\text{ bits},$$
$$D(q\|p) = 0.9\ln\frac{0.9}{0.5} + 0.1\ln\frac{0.1}{0.5} = 0.9(0.5878) + 0.1(-1.6094) = 0.3681\text{ nats} = 0.5310\text{ bits}.$$
Not symmetric. Interpretation: believing the coin is 90/10 when it is fair costs more than the reverse, because under $p$ the outcome tails (probability 0.5) is one you assigned only 0.1 to (surprise $\ln 5$), while under $q$ the rarely-seen tails is the only outcome where $p$ is badly wrong.

*Gaussians,* by (3.32). $D(\mathcal{N}(0,1)\,\|\,\mathcal{N}(1, 2^2)) = \ln\tfrac21 + \frac{1 + 1}{2\cdot4} - \tfrac12 = 0.6931 + 0.25 - 0.5 = 0.4431$ nats. $D(\mathcal{N}(1, 4)\,\|\,\mathcal{N}(0, 1)) = \ln\tfrac12 + \frac{4 + 1}{2\cdot1} - \tfrac12 = -0.6931 + 2.5 - 0.5 = 1.3069$ nats. The wide distribution sees the narrow one as very wrong (it puts almost no mass where the wide one lives); the narrow one sees the wide one as merely diffuse. The project verifies both closed forms by Monte Carlo: sample $X\sim f_1$, average $\log f_1(X) - \log f_2(X)$.

## 5. Why it matters and where it connects

- **ST-1** is built on four results from here. The likelihood of iid data is a product because the joint density of independent variables factors (Theorem 3.4). The MLE is consistent because the average log-likelihood obeys the law of large numbers (Theorem 3.34), and asymptotically normal because of the CLT (3.24). Maximizing likelihood is minimizing $D(p_{\text{data}}\,\|\,q_\theta)$ by (3.30). The Cramér–Rao bound is Cauchy–Schwarz in the inner product $\mathbb{E}[UV]$ of Theorem 3.17.
- **ST-3** takes Theorem 3.29 and (3.23) and replaces the population inner product $\mathbb{E}[UV]$ by the sample dot product. The normal equations are the orthogonality condition (3.22) restricted to linear functions. $R^2$ is the explained-variance fraction from Eve's law (3.21).
- **ST-5** generalizes Theorem 3.12 and the Gamma–Poisson posterior of §3.5 to every exponential family.
- **ST-6 and SM-6.** Metropolis–Hastings is Theorem 3.41 applied to a proposal-and-accept chain (exercise B10). MCMC averages are justified by the ergodic statement in Theorem 3.40. The factor $(1+\lambda_2)/(1-\lambda_2)$ in exercise B7 is the integrated autocorrelation time, which ST-6 turns into the effective sample size. SM-6 runs the same algorithm on the Ising model.
- **ST-8.** Equation (3.19) is the Kalman measurement update: $\mathbf{X}_1$ is the state, $\mathbf{X}_2$ the measurement, $\Sigma_{12}\Sigma_{22}^{-1}$ the gain, and the Schur complement the posterior covariance.
- **ML-1.** The target of squared-error learning is the regression function $\mathbb{E}[Y\mid X]$ (Theorem 3.29). The Pythagoras step in that proof is the first half of the bias–variance decomposition. The cross-entropy loss is (3.30).
- **ML-4.** In Gaussian-process regression, the posterior mean and covariance at test points are (3.19), with $\Sigma$ the kernel Gram matrix.
- **ML-6.** PCA finds the eigenvectors of $\Sigma$ (the remark after Theorem 3.19, and the long axis in Example 4.4). The mixture in exercise B3 is a two-component Gaussian mixture, and $P(Z = 1\mid Y = 3)$ in B3(d) is an E-step responsibility. EM never decreases the likelihood because of Gibbs' inequality (Theorem 3.45).
- **ML-10 and SM-7.** The KL term in a VAE loss is the formula of exercise C1 with $\Sigma_2 = I$. The gap between the log-likelihood and the ELBO is a KL divergence, so it is $\ge 0$ by Theorem 3.45.
- **SM-4.** The exponential calculation at the end of §3.13, with $x$ replaced by microstate energy, is the derivation of the Boltzmann distribution. Exercise C3 is the discrete version: the multiplier $\beta$ there plays the role of $1/k_BT$. The partition function is the normalizer $\sum_i e^{-\beta E_i}$.
- **LA links used here.** The spectral theorem from LA-8 gives the Gaussian ellipses. Cholesky and the Schur complement are from LA-10. Perron–Frobenius and power iteration come from LA-7 (PageRank). Orthogonal projection from LA-5 is conditional expectation.

## 6. Exercises

See [exercises.md](exercises.md). **A1–A5** are drills: a table that is uncorrelated but dependent, covariance algebra, Exponential order statistics and Poisson sums, a bivariate normal conditional, and entropies and KLs to compute by hand. **B1–B10** are the module: Box–Muller run backwards and the Cauchy ratio (B1), the vector best linear predictor (B2), Adam and Eve on a Gaussian mixture (B3), Beta–Binomial updating (B4), Gaussian marginals without a Gaussian joint (B5), Chebyshev against the CLT (B6), a law of large numbers for a Markov chain (B7), periodicity and the lazy walk (B8), conditioning and entropy (B9), and Metropolis (B10, which needs code). **B2, B3, B7 and B10 are essential.** **C1–C3** are the multivariate Gaussian KL, the mutual information of a Gaussian pair, and maximum entropy with a mean constraint (Jaynes' die).

## 7. Coding project

`project.py` builds every object in the module and asserts every number quoted in §4 and in the solutions. It runs in about 20 seconds.

1. **Discrete joint tables.** The script has `marginal`, `conditional`, `joint_moments`, and the conditional mean and variance by row. It reproduces Example 4.1 in floating point and Eve's law exactly in SymPy rationals ($7/30 + 1/60 = 1/4$).
2. **Box–Muller.** It implements (3.6), checks the Jacobian determinant 2 symbolically, and draws a QQ plot against `scipy.stats.norm.ppf` (figure `box_muller_qq.png`). It also runs B1: the radius is Rayleigh and $Z_1/Z_2$ is Cauchy.
3. **Sums and order statistics.** It convolves two Uniform densities numerically to get the triangle, checks that Exponentials sum to a Gamma (watch SciPy's *scale* $= 1/\text{rate}$ convention), that Uniform order statistics are Beta, and the Beta–Binomial posterior of B4.
4. **Covariance algebra and the two counterexamples** ($X, X^2$, and $(X, SX)$).
5. **The multivariate normal.** Cholesky sampling with the sample covariance asserted against $\Sigma$, then contours with the eigenvectors drawn at length $\sqrt{5.991\,\lambda_i}$ (figure `mvn_contours.png`). Mahalanobis $d^2$ is checked against $\chi^2_2$, and a slice at $|x_2 - 1| < 0.05$ is compared with (3.19).
6. **Conditional expectation.** It estimates $\mathbb{E}[Y\mid X]$ by binning. The residual is uncorrelated with $X$, $X^2$, $\sin 5X$ and $e^X$, while the best line's residual is correlated with $X^2$ (figure `conditional_expectation.png`). Also here: B2's vector predictor, B3's mixture, and Galton's heights (figure `galton.png`).
7. **Limit theorems.** CLT histograms for Uniform, Exponential, Bernoulli(0.1) and Pareto(1.5) at $n = 1, 2, 5, 30$ (figure `clt_histograms.png`). The exact sup error $\sup_x|F_n - \Phi|$ is plotted against $n$ on log axes (figure `clt_error_vs_n.png`). Its slope is $-\tfrac12$ for the skewed cases, as Berry–Esseen predicts, and $-1$ for the symmetric Uniform; Pareto(1.5) has no slope. The continuity correction, B6 and B7 are checked here too.
8. **Markov chains.** `stationary(P)` returns the eigenvector of $P^\top$, and `simulate_chain` runs a chain. For the weather chain, the mixing curve has slope $\ln\lambda_2$ (figure `markov_mixing.png`). The script then checks detailed balance on the weather chain and shows it failing on a circulating chain. It also covers the graph walk with $\boldsymbol\pi \propto$ degree, B8's bipartite walk and its lazy fix, B10's Metropolis chain, and PageRank by power iteration.
9. **Information theory.** The Bernoulli entropy curve (figure `bernoulli_entropy.png`), and Gibbs and $H(p,q) = H(p) + D(p\|q)$ on 1000 random pairs. KLs between Gaussians in closed form are compared with Monte Carlo and with SymPy integration. Entropies are compared at equal variance (Gaussian > Laplace > Uniform > Exponential), and A5, B9, C1–C3 are checked.

**Build it yourself first:** `marginal` and `conditional` (Part 1), `box_muller` (Part 2), `sample_mvn` and `mvn_conditional` (Part 5), and `stationary` (Part 8). Each is under ten lines. Test your `mvn_conditional` on a $\Sigma$ with unequal diagonal entries too: in Example 4.4 both diagonal blocks equal 2, so inverting the wrong block goes unnoticed there.

**Your turn** (end of the script) has four tasks. First, the return probability of a 2D random walk, which is recurrent but slow; then 3D, which returns with probability about 0.34. Second, the entropy of B3's mixture by Monte Carlo, which must come out below $2.249$ nats by Theorem 3.46. Third, $I(X;Y) = -\tfrac12\ln(1-\rho^2)$ across $\rho$, then $I(X; X^2)$ by histogram, to see mutual information detect what covariance misses. Fourth, Metropolis with a continuous target.

## 8. Where you'll struggle

- **Jacobians.** You will use the Jacobian of the wrong map, or drop the absolute value. The rule: it is $\partial(\text{old})/\partial(\text{new})$, because you substitute the old variables *out*. If $\partial(\text{new})/\partial(\text{old})$ is easier to compute, compute it and take the reciprocal of its determinant. Then check two things every time: that the result integrates to 1, and that the new *domain* is right. Most wrong answers in exercise B1 come from the domain, not the determinant. Theorem 3.11 is worth redoing from a blank page.
- **The MVN conditional derivation is long.** It is the longest derivation in the module and you will not hold it in your head on first reading. Do it once in 2D with scalars, completing the square in $x_1$ by hand, and then once with blocks. Memorize the three facts in the "What just happened" after (3.19), not the algebra. ST-8 makes you do it again, and the second time is when it sticks.
- **$\mathbb{E}[Y\mid X]$ as a random variable** takes about a week to settle. The fix is a type discipline. First write the function $g(x) = \mathbb{E}[Y\mid X = x]$ as an ordinary function, then substitute $X$. Ask of every expression "is this a number or a random variable?" Adam's law is LOTUS applied to $g$; once you see that, it stops being mysterious.
- **Eve's law: you will swap the terms.** Say "expected conditional variance plus variance of the conditional expectation." Check with two extremes. If $Y \perp X$, the second term is 0. If $Y = X$, the first term is 0.
- **Row versus column vectors for chains.** Here $\boldsymbol\pi P$ uses rows and rows of $P$ sum to 1. LA-7 used columns. `np.linalg.eig(P)` returns *right* eigenvectors, so for a stationary distribution you must call it on `P.T`. If your "stationary distribution" has negative entries or does not satisfy `pi @ P == pi`, this is why.
- **"KL is a distance."** It is not: it is asymmetric (Example 4.7) and fails the triangle inequality. It is infinite when $q = 0$ somewhere $p > 0$. Always say which argument holds the truth. ML minimizes $D(p_{\text{data}}\,\|\,q_\theta)$, while variational inference minimizes $D(q\,\|\,p)$, and the two behave differently (ML-10).
- **Differential entropy can be negative**, and it changes under rescaling. Neither is a bug. If a quantity must be invariant, use KL or mutual information, which are.
- **The CLT is about standardized sums**, not about your data becoming normal. "$n = 30$ is enough" is folklore. The project's Bernoulli(0.1) panel is still visibly skewed at $n = 30$, and the Pareto(1.5) panel never converges.
- **Exercises harder than they look.** B2: do not differentiate with respect to a matrix; use the orthogonality trick of Theorem 3.29. B7(c): the double sum over $|s - t|$ needs care at the diagonal. B8: you need to say *why* the distance stays at 1, not only observe it. C3: the multiplier has no closed form, so you need a root-finder.

## 9. Check yourself

1. Given a joint PMF table, how do you get $p_X$, $p_{Y|X}(\cdot\mid x)$, and test independence? Which single cell suffices to refute it?
2. Why does $f(x,y) = 2$ on $0 < x < y < 1$ fail characterization (c) of Theorem 3.4?
3. Write the 2D change-of-variables formula. Which Jacobian goes in, and why the absolute value?
4. What is the distribution of the sum of $n$ iid $\mathrm{Exponential}(\lambda)$? Of the $k$-th smallest of $n$ iid Uniforms?
5. Prove $|\rho| \le 1$ in three lines. Which inner product did you use?
6. Give a pair that is uncorrelated but dependent. Why can't it be jointly Gaussian?
7. Write the density of $\mathcal{N}(\boldsymbol\mu, \Sigma)$. What are the axes and semi-axes of its level sets?
8. State (3.19) from memory. Which part does not depend on $\mathbf{x}_2$?
9. State Adam's and Eve's laws. What is $\mathbb{E}[Y\mid X]$ the orthogonal projection of, onto what?
10. Prove the weak law of large numbers in two lines. Where is independence used?
11. State the CLT precisely. What step of the MGF proof is cited, not proved? Name a distribution for which the CLT fails.
12. What is a stationary distribution, as an eigenvector statement? Prove detailed balance implies stationarity.
13. Show $D(p\|q) \ge 0$ using Jensen. Why is cross-entropy minimization the same as maximum likelihood?
14. Which density maximizes entropy at fixed variance? At fixed mean on $[0, \infty)$?

## 10. Pair with

- **Blitzstein & Hwang, *Introduction to Probability*, 2nd ed. (free at probabilitybook.net), ch. 7–11.** The chapters map to this module as follows: ch. 7 is joint distributions and covariance, ch. 8 transformations, convolutions, Beta, Gamma and order statistics, ch. 9 conditional expectation, ch. 10 inequalities and limit theorems, ch. 11 Markov chains. Their story proofs (bank–post office, Beta–Gamma) are better than any algebra. Their ch. 9 is the best treatment of $\mathbb{E}[Y\mid X]$ as a random variable at this level. Do the end-of-chapter problems.
- **Harvard Stat 110, lectures 16–30 on YouTube**, plus 31–33 for Markov chains. Lectures 19–22 cover joint distributions, covariance and transformations, 23–25 Beta, Gamma and order statistics, 25–27 conditional expectation, and 28–30 inequalities, LLN, CLT and the multivariate normal. Watch the lectures on conditional expectation even if you skip the rest.
- **MacKay, *Information Theory, Inference, and Learning Algorithms* (free online), ch. 1–2 and 8.** Ch. 2 defines entropy, KL and Gibbs' inequality with many small examples. Ch. 8 covers joint and conditional entropy and mutual information with diagrams. MacKay explains *why* entropy measures information (compression) far better than §3.13 has room for.
- **Bishop, *Pattern Recognition and Machine Learning*, §2.3.** This is the Gaussian conditional and marginal derived through the precision matrix $\Sigma^{-1}$ rather than the block LDU of §3.8. Reading both derivations is the best way to make (3.19) stick. §2.3.3 (Bayes' theorem for linear-Gaussian models) is the bridge to ST-8 and ML-4.
- **Optional:** Levin & Peres, *Markov Chains and Mixing Times* (free online), ch. 1–4, for the full convergence theorem behind Theorem 3.40.
