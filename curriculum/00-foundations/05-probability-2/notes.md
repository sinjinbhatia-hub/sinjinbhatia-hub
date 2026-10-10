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

