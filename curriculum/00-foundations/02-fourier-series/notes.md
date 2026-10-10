# F2: Fourier Series

**Time:** ~10 hours · **Prerequisites:** LA-5 (projections, orthonormal bases), LA-8 (inner product spaces, self-adjoint operators) · **Unlocks:** F3, EM-3, QM-3, ST-8

> **Conventions in this module.** Functions live on $[-\pi, \pi]$ and are extended $2\pi$-periodically. The inner product is $\langle f, g\rangle = \int_{-\pi}^{\pi}\overline{f(x)}\,g(x)\,dx$, conjugate-linear in the *first* slot (the physics convention of LA-8), with norm $\|f\| = \sqrt{\langle f,f\rangle}$. **Real form:** $f \sim \dfrac{a_0}{2} + \sum_{n=1}^{\infty}(a_n\cos nx + b_n\sin nx)$ with $a_n = \dfrac{1}{\pi}\int_{-\pi}^{\pi} f(x)\cos nx\,dx$ for $n \ge 0$ and $b_n = \dfrac{1}{\pi}\int_{-\pi}^{\pi} f(x)\sin nx\,dx$ for $n\ge1$. **Complex form:** $f \sim \sum_{n=-\infty}^{\infty} c_n e^{inx}$ with $c_n = \dfrac{1}{2\pi}\int_{-\pi}^{\pi} f(x)e^{-inx}\,dx$; the doubly infinite sum always means the symmetric limit $\lim_{N\to\infty}\sum_{|n|\le N}$. **Relationship:** $c_0 = a_0/2$, and for $n \ge 1$, $c_n = (a_n - ib_n)/2$, $c_{-n} = (a_n + ib_n)/2$; conversely $a_n = c_n + c_{-n}$ and $b_n = i(c_n - c_{-n})$. For real $f$, $c_{-n} = \overline{c_n}$. **General period $T$:** replace $e^{inx}$ by $e^{2\pi i n t/T}$ and $\frac{1}{2\pi}\int_{-\pi}^{\pi}$ by $\frac{1}{T}\int$ over any interval of length $T$. The symbol $\sim$ means "has these Fourier coefficients"; it becomes $=$ only after a convergence theorem (§3.8). **Other books differ:** Osgood uses period $1$ and $e^{2\pi i n t}$; some texts write the constant term as $a_0$ (not $a_0/2$) and then need $a_0 = \frac{1}{2\pi}\int f$; some put $\frac{1}{\sqrt{2\pi}}$ on both basis functions and coefficients so the basis is orthonormal. Every one of these is the same projection with the constants moved around. When in doubt, recompute the constant from the projection formula (3.5).

## 0. Where this sits

You have met Fourier series before, probably as a recipe: write down the integrals for $a_n$ and $b_n$, compute them for a square wave, watch a sum of sines approximate a step. This module replaces the recipe with a picture you now own from LA-5 and LA-8. A periodic function is a vector in an inner product space; $\{1, \cos nx, \sin nx\}$ is an orthogonal basis of that space; the Fourier coefficients are the components of $f$ in that basis, which is orthogonal projection with the dot product replaced by an integral. The reason the basis is orthogonal at all is LA-8's spectral theorem applied to the operator $d^2/dx^2$.

What it makes possible: F3 (the Fourier transform is the $T\to\infty$ limit of the period-$T$ series; the DFT is its sampled version), EM-3 (separation of variables in boundary-value problems, where Legendre polynomials are the next orthogonal basis), QM-3 (the infinite well's energy eigenstates are $\sin(n\pi x/L)$, and expanding a wavefunction in them is this module), ST-8 (the periodogram of a time series is $|c_n|^2$), and CM-6 (normal modes of a string are the continuum limit of the coupled-oscillator eigenproblem).

## 1. The problem that forced this idea

**The vibrating string, 1747–1753.** A string of length $L$ fixed at both ends has displacement $u(x,t)$ obeying the wave equation $u_{tt} = c^2 u_{xx}$. Jean le Rond d'Alembert solved it in 1747: $u = \phi(x+ct) + \psi(x - ct)$, two traveling waves, with $\phi, \psi$ fixed by the initial shape and velocity. Euler, in 1748, insisted that the initial shape could be *any* curve you could draw, including one with a corner, like a plucked string; d'Alembert replied that his solution only made sense for shapes with two derivatives.

Daniel Bernoulli entered in 1753 with a different kind of answer. The string has standing-wave solutions $\sin(n\pi x/L)\cos(n\pi c t/L)$ for every $n = 1, 2, 3, \dots$ (the fundamental and its overtones, which musicians had known about for a century). The wave equation is linear, so any sum $\sum_n b_n \sin(n\pi x/L)\cos(n\pi ct/L)$ is also a solution. Bernoulli claimed that *this* is the general solution: every motion of the string is a superposition of its harmonics. Set $t = 0$ and the claim becomes: every initial shape is $\sum_n b_n \sin(n\pi x/L)$.

Euler and d'Alembert both rejected this, and the objection is worth understanding because it is almost right. A sum of sines is odd and $2L$-periodic; how can it equal an arbitrary function that is neither? Each sine is infinitely differentiable; how can a sum of them have a corner? The resolution is that the equality is only claimed on $[0, L]$, and that an *infinite* sum of smooth functions need not be smooth. Neither point was clear in 1753, and nobody had a formula for the $b_n$. (Euler did derive coefficient integrals in 1777 by the orthogonality argument of §3.3, for a related problem, without connecting it to the general claim.) The dispute went unresolved for fifty years.

**Heat, 1807.** Joseph Fourier, who had gone to Egypt with Napoleon and was by then prefect of the Isère, took up the problem of how heat spreads in a solid. He derived the heat equation $u_t = k\,u_{xx}$ and, for a bar with its ends held at zero temperature, found the structure Bernoulli had found: special solutions $\sin(n\pi x/L)e^{-k(n\pi/L)^2 t}$, linearity, and so the need to write an arbitrary initial temperature as $\sum b_n \sin(n\pi x/L)$. Unlike Bernoulli, Fourier gave the formula for the coefficients, computed many examples including functions with jumps, and claimed flatly that *any* function on an interval has such an expansion.

He submitted this in December 1807 to the Institut de France as a memoir on the propagation of heat in solid bodies. The referees were Lagrange, Laplace, Monge, and Lacroix. Lagrange, who had worked on the string problem himself in 1759 and did not believe a trigonometric series could represent a function with a corner, objected, and the memoir was not published. The Institut set the propagation of heat as the subject of its prize for 1811; Fourier's revised essay won, but the committee's report said the analysis left something to be desired in generality and even in rigor, and again it was not printed. Fourier published the theory himself in 1822 as the *Théorie analytique de la chaleur*.

**Dirichlet, 1829.** The rigor Lagrange wanted arrived in a paper by Peter Gustav Lejeune Dirichlet, who proved: if $f$ is bounded on $[-\pi,\pi]$ with finitely many maxima and minima and finitely many jumps, then its Fourier series converges at every point to the average of the left and right limits of $f$. It was the first convergence theorem for Fourier series, and the first time anyone had been precise about what "function" means. Lagrange's instinct was half right: the series does not always converge, and for merely continuous functions it can fail at a point (du Bois-Reymond, 1873). Fourier's instinct was more right: for anything you will meet in physics, it works.

**The problem in one sentence.** The heat equation and the wave equation on a finite interval have easy solutions for sinusoidal initial data; linearity turns sums of easy solutions into solutions; the question is whether *every* initial condition is such a sum and, if so, how to find the coefficients. Linear algebra answers both halves: the sines are an orthogonal basis, and the coefficients are projections.

## 2. Intuition

### A function is a vector with infinitely many components

In LA-1 a vector in $\mathbb{R}^n$ was a list of $n$ numbers. A function $f$ on $[-\pi, \pi]$ is a list of numbers too, one for each $x$: uncountably many components, indexed by a real number instead of an integer. You add functions pointwise and scale them pointwise, so they form a vector space (LA-3). What LA-8 added is that this space can carry an inner product. For vectors, $\mathbf{u}\cdot\mathbf{v} = \sum_i u_i v_i$. For functions the index $i$ becomes $x$ and the sum becomes an integral:

$$\langle f, g \rangle = \int_{-\pi}^{\pi} f(x)\,g(x)\,dx \qquad\text{(real functions; conjugate $f$ for complex ones).}$$

Everything LA-5 did with the dot product survives: lengths $\|f\|^2 = \int f^2$, angles, perpendicularity, projection.

### Sampling makes this concrete

If you don't trust the integral, sample. Pick $N$ equally spaced points $x_j$ in $[-\pi, \pi)$ and replace $f$ by the vector $\mathbf{f} = (f(x_0), \dots, f(x_{N-1}))^\top \in \mathbb{R}^N$. Then $\mathbf{f}\cdot\mathbf{g} = \sum_j f(x_j)g(x_j) \approx \frac{N}{2\pi}\int f g\,dx$, a Riemann sum. The integral inner product is the $N\to\infty$ limit of the ordinary dot product of sample vectors, up to the constant $2\pi/N$. The project does this and finds that sampled cosines are orthogonal *as vectors in $\mathbb{R}^N$*, exactly, not approximately.

### The sines and cosines are perpendicular

Draw $\cos x$ and $\cos 2x$ on $[-\pi, \pi]$. Their product is positive where they have the same sign and negative where they don't, and the regions balance exactly: $\int \cos x\cos 2x\,dx = 0$. The same happens for any two distinct members of

$$\{1,\ \cos x,\ \sin x,\ \cos 2x,\ \sin 2x,\ \cos 3x,\ \sin 3x,\ \dots\}.$$

This is a set of mutually perpendicular vectors in function space, infinitely many of them. In $\mathbb{R}^n$ you can have at most $n$ mutually perpendicular nonzero vectors (LA-3), so this space is infinite-dimensional. Their lengths are not $1$: $\|1\|^2 = \int 1\,dx = 2\pi$, and $\|\cos nx\|^2 = \int\cos^2 nx\,dx = \pi$ (the average of $\cos^2$ over whole periods is $\frac12$, times the interval length $2\pi$). The set is orthogonal but not orthonormal.

```
 cos x:      ____                     cos 2x:   __          __
          __/    \__                           /  \        /  \
   -pi  _/          \_  pi              -pi  _/    \__  __/    \_  pi
                                                        \/
 product cos x * cos 2x: (+)(-)(+)(-) in four lobes that cancel; integral = 0
```

### Guess the coefficient formula before reading §3

LA-5 gave you this: if $\mathbf{q}_1, \dots, \mathbf{q}_k$ are mutually orthogonal (not necessarily unit) vectors and $\mathbf{v}$ is in their span, then

$$\mathbf{v} = \sum_k \frac{\mathbf{q}_k\cdot\mathbf{v}}{\mathbf{q}_k\cdot\mathbf{q}_k}\,\mathbf{q}_k,$$

and if $\mathbf{v}$ is *not* in their span, the same formula gives the closest point in the span, the projection. The coefficient of $\mathbf{q}_k$ is "the shadow of $\mathbf{v}$ on $\mathbf{q}_k$, divided by the squared length of $\mathbf{q}_k$." Dot products between different $\mathbf{q}$'s vanish, which is why each coefficient can be read off independently.

Now replace $\mathbf{q}_k$ by $\cos nx$ and the dot product by the integral. The coefficient of $\cos nx$ in $f$ should be

$$\frac{\langle \cos nx, f\rangle}{\langle \cos nx, \cos nx\rangle} = \frac{\int_{-\pi}^\pi f(x)\cos nx\,dx}{\pi}.$$

That is $a_n$. The coefficient of $\sin nx$ should be $\frac1\pi\int f\sin nx\,dx = b_n$. The coefficient of the constant function $1$ should be $\frac{\int f\,dx}{2\pi}$, which is $\frac{a_0}{2}$ with $a_0 = \frac1\pi\int f\,dx$: the odd-looking $\frac12$ on $a_0$ is nothing but $\|1\|^2 = 2\pi$ versus $\|\cos nx\|^2 = \pi$. Every normalization constant in this subject is a $\|\mathbf{q}_k\|^2$. Remember only that and you will never get a factor of $2$ or $\pi$ wrong.

Things you should be able to guess by the end of this section, and that §3 proves:

1. The partial sum $S_N f$ built from the first $N$ harmonics is the *closest* trigonometric polynomial of degree $N$ to $f$ in mean square (LA-5's projection theorem).
2. The coefficients can't be too large: $\sum (a_n^2 + b_n^2) \le \frac1\pi\|f\|^2$ (a shadow is no longer than the vector; LA-5's Bessel inequality), with equality when the basis is "complete."
3. Smooth functions have small high-frequency coefficients, because high harmonics are what wiggles are made of.
4. The basis is orthogonal for a *reason*: $\cos nx$ and $\sin nx$ are eigenfunctions of $d^2/dx^2$ with eigenvalue $-n^2$, and that operator is self-adjoint (LA-8).

## 3. Formal development

### 3.1 The space and its inner product

**Definition 3.1 (Piecewise continuous, piecewise smooth).** $PC$ is the set of complex-valued functions on $[-\pi,\pi]$ that are continuous except at finitely many points, at each of which the one-sided limits $f(x^-)$ and $f(x^+)$ exist and are finite. Every $f\in PC$ is extended to $\mathbb{R}$ by $f(x + 2\pi) = f(x)$. A function is **piecewise smooth** if both $f$ and $f'$ are in $PC$ (the derivative taken on each open subinterval between breaks, with one-sided limits at the ends).

On $PC$ define

$$\langle f, g\rangle = \int_{-\pi}^{\pi} \overline{f(x)}\,g(x)\,dx, \qquad \|f\| = \langle f, f\rangle^{1/2} = \left(\int_{-\pi}^{\pi}|f(x)|^2\,dx\right)^{1/2}. \tag{3.1}$$

Check the LA-8 axioms. *Conjugate symmetry:* $\langle g, f\rangle = \int\bar g f = \overline{\int \bar f g} = \overline{\langle f, g\rangle}$. *Linearity in the second slot:* $\langle f, \alpha g + \beta h\rangle = \alpha\langle f, g\rangle + \beta\langle f, h\rangle$ by linearity of the integral; conjugate-linearity in the first slot follows from the first two. *Positivity:* $\langle f, f\rangle = \int|f|^2 \ge 0$. *Definiteness* needs one sentence of care: $\int|f|^2 = 0$ forces $f = 0$ at every point where $f$ is continuous but says nothing about the finitely many exceptional points, so we agree to identify two functions that differ only at finitely many points (in general, on a set of measure zero). With that identification (3.1) is an inner product and $PC$ is an inner product space in LA-8's sense, so Cauchy–Schwarz $|\langle f,g\rangle| \le \|f\|\|g\|$ and the triangle inequality hold, because LA-8 proved them from the axioms alone.

$\|f - g\|$ is the **mean-square distance** between $f$ and $g$, the natural notion of "close" in this module. It is not the only one: $\max_x|f(x) - g(x)|$ (uniform distance) is stronger, and §3.9–3.10 are about when the two agree and when they do not.

### 3.2 Orthogonality

**Theorem 3.2 (Orthogonality of the trigonometric system).** For integers $m, n \ge 1$,

$$\int_{-\pi}^{\pi}\cos mx\cos nx\,dx = \pi\,\delta_{mn}, \qquad \int_{-\pi}^{\pi}\sin mx\sin nx\,dx = \pi\,\delta_{mn}, \qquad \int_{-\pi}^{\pi}\cos mx\sin nx\,dx = 0, \tag{3.2}$$

where $\delta_{mn} = 1$ if $m = n$ and $0$ otherwise. Also $\int_{-\pi}^\pi \cos nx\,dx = \int_{-\pi}^\pi \sin nx\,dx = 0$ for $n\ge1$, and $\int_{-\pi}^\pi 1\,dx = 2\pi$.

▶ **Try it first.** Everything rests on two facts: $\int_{-\pi}^\pi\cos kx\,dx = 0$ for every nonzero integer $k$, and $\int_{-\pi}^\pi \sin kx\,dx = 0$ for every integer $k$. Prove those, then turn each product in (3.2) into a sum using the addition formulas.

**Derivation.**

1. For a nonzero integer $k$, $\int_{-\pi}^\pi\cos kx\,dx = \left[\frac{\sin kx}{k}\right]_{-\pi}^{\pi} = \frac{\sin k\pi - \sin(-k\pi)}{k} = 0$, because $\sin k\pi = 0$. For $k = 0$ the integral is $2\pi$. For any integer $k \ne 0$, $\int_{-\pi}^{\pi}\sin kx\,dx = \left[-\frac{\cos kx}{k}\right]_{-\pi}^{\pi} = \frac{-\cos k\pi + \cos(-k\pi)}{k} = 0$ because cosine is even; for $k = 0$ the integrand is $0$.
2. Product-to-sum formulas. Add $\cos(A - B) = \cos A\cos B + \sin A\sin B$ and $\cos(A+B) = \cos A\cos B - \sin A\sin B$ to get $2\cos A\cos B = \cos(A-B) + \cos(A+B)$; subtract them to get $2\sin A\sin B = \cos(A-B) - \cos(A+B)$. Add $\sin(A+B) = \sin A\cos B + \cos A\sin B$ and $\sin(A - B) = \sin A\cos B - \cos A\sin B$ to get $2\sin A\cos B = \sin(A-B) + \sin(A+B)$.
3. Cosine–cosine, with $A = mx$, $B = nx$:
   $$\int_{-\pi}^{\pi}\cos mx\cos nx\,dx = \tfrac12\int_{-\pi}^{\pi}\cos((m-n)x)\,dx + \tfrac12\int_{-\pi}^{\pi}\cos((m+n)x)\,dx.$$
   Since $m, n\ge1$, $m + n \ne 0$, so the second integral is $0$ by step 1. The first is $0$ if $m \ne n$ and $\frac12\cdot 2\pi = \pi$ if $m = n$.
4. Sine–sine: the same two integrals with a minus sign on the second, which is $0$ anyway. Result $\pi\delta_{mn}$.
5. Sine–cosine: $\int\sin nx\cos mx\,dx = \frac12\int\sin((n-m)x)\,dx + \frac12\int\sin((n+m)x)\,dx = 0 + 0$ by step 1, for all $m, n$ (including $m = n$, where the first integrand is $\sin 0 = 0$).
6. The constant: $\int 1\cdot\cos nx = 0$ and $\int 1\cdot\sin nx = 0$ are step 1 with $k = n \ge 1$; $\int 1\cdot 1 = 2\pi$. ∎

■ **What just happened.** Every product of two basis functions is a sum of pure sinusoids of frequency $m \pm n$, and a sinusoid of nonzero integer frequency integrates to zero over a full period. The only survivor is the frequency-zero term that appears when $m = n$, and it contributes $\frac12\cdot 2\pi = \pi$. So the norms are
$$\|1\|^2 = 2\pi, \qquad \|\cos nx\|^2 = \|\sin nx\|^2 = \pi \quad (n \ge 1).$$

**Theorem 3.3 (Orthogonality of the complex exponentials).** For integers $m, n$,
$$\langle e^{imx}, e^{inx}\rangle = \int_{-\pi}^{\pi} e^{-imx}e^{inx}\,dx = \int_{-\pi}^{\pi} e^{i(n-m)x}\,dx = 2\pi\,\delta_{mn}. \tag{3.3}$$

*Proof.* If $m = n$ the integrand is $1$ and the integral is $2\pi$. If $k = n - m \ne 0$, then $\int_{-\pi}^{\pi} e^{ikx}dx = \left[\frac{e^{ikx}}{ik}\right]_{-\pi}^{\pi} = \frac{e^{ik\pi} - e^{-ik\pi}}{ik} = \frac{(-1)^k - (-1)^k}{ik} = 0$, using $e^{\pm ik\pi} = \cos k\pi \pm i\sin k\pi = (-1)^k$. ∎

Note the conjugate on the first slot: $\langle e^{imx}, e^{inx}\rangle$ has $e^{-imx}$ in the integrand. That minus sign is where the $e^{-inx}$ in the coefficient formula comes from. One line, versus six steps for the real case: the first of many reasons to prefer the complex form.

### 3.3 The coefficients are projections

▶ **Try it first.** Suppose $f(x) = \frac{a_0}{2} + \sum_{k=1}^{N}(a_k\cos kx + b_k\sin kx)$ with unknown coefficients. Take the inner product of both sides with $\cos nx$ and use Theorem 3.2 to isolate $a_n$. Then do the same with $\sin nx$ and with $1$.

**Derivation.**

1. Pair both sides with $\cos nx$, $1\le n\le N$, using linearity in the second slot:
   $$\langle\cos nx, f\rangle = \frac{a_0}{2}\langle\cos nx, 1\rangle + \sum_{k=1}^N\Big(a_k\langle\cos nx,\cos kx\rangle + b_k\langle\cos nx, \sin kx\rangle\Big).$$
2. By Theorem 3.2 every inner product on the right is zero except $\langle\cos nx,\cos nx\rangle = \pi$. So $\langle\cos nx, f\rangle = \pi a_n$.
3. Pair with $\sin nx$: only $b_n\langle\sin nx, \sin nx\rangle = \pi b_n$ survives.
4. Pair with $1$: only $\frac{a_0}{2}\langle 1, 1\rangle = \frac{a_0}{2}\cdot 2\pi = \pi a_0$ survives, so $\langle 1, f\rangle = \pi a_0$, which is the $n = 0$ case of the cosine formula since $\cos 0x = 1$. ∎

This gives the real coefficient formulas (repeating the conventions box):

$$a_n = \frac{1}{\pi}\int_{-\pi}^{\pi} f(x)\cos nx\,dx \ \ (n\ge0), \qquad b_n = \frac{1}{\pi}\int_{-\pi}^{\pi} f(x)\sin nx\,dx\ \ (n\ge1). \tag{3.4}$$

Each is an instance of one rule, the LA-5 projection rule:

$$\boxed{\ \text{coefficient of } \varphi_n \text{ in } f \;=\; \frac{\langle\varphi_n, f\rangle}{\langle\varphi_n,\varphi_n\rangle}\ } \tag{3.5}$$

for any orthogonal family $\{\varphi_n\}$. The $\frac1\pi$ is $1/\|\cos nx\|^2$; the $\frac{a_0}{2}$ is $\langle 1, f\rangle/\|1\|^2 = \frac{1}{2\pi}\int f$. (3.5) is the formula to memorize. Everything in the conventions box is (3.5) with a particular $\|\varphi_n\|^2$ plugged in.

■ **What just happened.** The derivation assumed $f$ is a trigonometric polynomial. But the *formulas* (3.4) make sense for any $f \in PC$, so we turn them into a definition.

**Definition 3.4 (Fourier coefficients, series, partial sums).** For $f \in PC$, the numbers $a_n, b_n$ of (3.4) are its **Fourier coefficients**; the formal series $\frac{a_0}{2} + \sum_{n\ge1}(a_n\cos nx + b_n\sin nx)$ is its **Fourier series**; $f \sim \frac{a_0}{2} + \sum(\cdots)$ means "these are the coefficients of $f$," with no claim yet about convergence. The **$N$-th partial sum** is

$$S_N f(x) = \frac{a_0}{2} + \sum_{n=1}^{N}\big(a_n\cos nx + b_n\sin nx\big). \tag{3.6}$$

$S_N f$ lies in $\mathcal{T}_N = \operatorname{span}\{1, \cos x, \sin x, \dots, \cos Nx, \sin Nx\}$, the **trigonometric polynomials of degree $\le N$**, a subspace of dimension $2N+1$.

**Theorem 3.5 (Best mean-square approximation).** Let $f \in PC$. Then $f - S_N f$ is orthogonal to every element of $\mathcal{T}_N$, and for every $g \in \mathcal{T}_N$,
$$\|f - S_N f\| \le \|f - g\|,$$
with equality if and only if $g = S_N f$. In words: $S_N f$ is the orthogonal projection of $f$ onto $\mathcal{T}_N$, the unique trigonometric polynomial of degree $\le N$ closest to $f$ in mean square.

▶ **Try it first.** This is LA-5's projection theorem word for word. Show the residual $f - S_Nf$ is orthogonal to each basis function $\varphi \in \{1, \cos kx, \sin kx\}_{k\le N}$ by computing $\langle \varphi, f - S_N f\rangle$, then use Pythagoras on $f - g = (f - S_N f) + (S_N f - g)$.

**Derivation.**

1. Let $\varphi$ be one of the $2N+1$ basis functions of $\mathcal{T}_N$ and $\gamma_\varphi = \langle\varphi, f\rangle/\|\varphi\|^2$ its coefficient in $S_Nf$, which by (3.5) is what (3.6) says. Then
   $$\langle\varphi, f - S_Nf\rangle = \langle\varphi, f\rangle - \sum_{\psi}\gamma_\psi\langle\varphi,\psi\rangle = \langle\varphi, f\rangle - \gamma_\varphi\|\varphi\|^2 = 0,$$
   because by Theorem 3.2 the only nonzero term in the sum over basis functions $\psi$ is $\psi = \varphi$.
2. Every $g \in \mathcal{T}_N$ is a linear combination of the $\varphi$'s, so by linearity $\langle g, f - S_N f\rangle = 0$: the residual is orthogonal to all of $\mathcal{T}_N$.
3. For $g \in \mathcal{T}_N$ write $f - g = (f - S_Nf) + (S_N f - g)$. The second piece is in $\mathcal{T}_N$, so it is orthogonal to the first by step 2. Expanding $\|u + v\|^2 = \|u\|^2 + 2\operatorname{Re}\langle u, v\rangle + \|v\|^2$ with $\langle u,v\rangle = 0$ (Pythagoras):
   $$\|f - g\|^2 = \|f - S_Nf\|^2 + \|S_Nf - g\|^2 \ \ge\ \|f - S_Nf\|^2.$$
4. Equality holds iff $\|S_Nf - g\| = 0$ iff $g = S_Nf$. ∎

■ **What just happened.** Nothing new: this is Example 4.5 of LA-1 and the projection theorem of LA-5 with the dot product replaced by (3.1). The ambient space is now infinite-dimensional, and that changes nothing in the proof, because every step happened inside the finite-dimensional $\mathcal{T}_N$ plus the one vector $f$. Taking $g = 0$ in step 3 gives the identity
$$\|f\|^2 = \|S_Nf\|^2 + \|f - S_Nf\|^2,$$
which runs §3.6.

### 3.4 The complex form

Euler's formula gives $\cos nx = \frac{e^{inx} + e^{-inx}}{2}$ and $\sin nx = \frac{e^{inx} - e^{-inx}}{2i}$. Substitute into one term of the real series, using $\frac1i = -i$:

$$a_n\cos nx + b_n\sin nx = \frac{a_n}{2}(e^{inx} + e^{-inx}) - \frac{ib_n}{2}(e^{inx} - e^{-inx}) = \frac{a_n - ib_n}{2}e^{inx} + \frac{a_n + ib_n}{2}e^{-inx}.$$

So $S_N f = \sum_{|n|\le N} c_n e^{inx}$ with

$$c_0 = \frac{a_0}{2}, \qquad c_n = \frac{a_n - ib_n}{2}, \qquad c_{-n} = \frac{a_n + ib_n}{2} \quad (n\ge1). \tag{3.7}$$

Check that this agrees with the projection rule (3.5) applied to the orthogonal family $\{e^{inx}\}$ of Theorem 3.3, whose norms are $\|e^{inx}\|^2 = 2\pi$:

$$c_n = \frac{\langle e^{inx}, f\rangle}{\langle e^{inx}, e^{inx}\rangle} = \frac{1}{2\pi}\int_{-\pi}^{\pi}f(x)e^{-inx}\,dx = \frac{1}{2\pi}\int f(x)(\cos nx - i\sin nx)\,dx = \frac{a_n - ib_n}{2}. \tag{3.8}$$

✓. Inverting (3.7): $a_n = c_n + c_{-n}$ and $b_n = i(c_n - c_{-n})$. If $f$ is real-valued, $\overline{c_n} = \frac{1}{2\pi}\int f e^{+inx}dx = c_{-n}$, so the pair $c_n, c_{-n}$ carries the same information as $a_n, b_n$: $|c_n| = \frac12\sqrt{a_n^2 + b_n^2}$ is half the **amplitude** of harmonic $n$ and $\arg c_n$ is its **phase**.

**General period.** If $f$ has period $T$, the substitution $x = 2\pi t/T$ maps $t\in[-T/2, T/2]$ to $[-\pi,\pi]$, turns $e^{inx}$ into $e^{2\pi i n t/T}$, and gives $dx = \frac{2\pi}{T}dt$:

$$f(t) \sim \sum_{n=-\infty}^{\infty} c_n e^{2\pi i n t/T}, \qquad c_n = \frac{1}{T}\int_{-T/2}^{T/2} f(t)\,e^{-2\pi i n t/T}\,dt \quad\text{(any interval of length $T$ will do).} \tag{3.9}$$

The frequencies are $n/T$: the fundamental $1/T$ and its integer multiples. Osgood's convention is (3.9) with $T = 1$. F3 starts from (3.9): as $T \to \infty$ the frequencies $n/T$ crowd into a continuum, the sum becomes an integral, and $c_n$ becomes the Fourier transform.

### 3.5 Even and odd functions; half-range series

If $f$ is even, $f(x)\sin nx$ is odd, so $b_n = 0$; and $f(x)\cos nx$ is even, so $a_n = \frac{2}{\pi}\int_0^\pi f\cos nx\,dx$. If $f$ is odd, $a_n = 0$ and $b_n = \frac2\pi\int_0^\pi f\sin nx\,dx$. An even function has a **cosine series**; an odd function has a **sine series**.

Now suppose $f$ is given only on $[0, L]$. You may extend it to $[-L, 0]$ any way you like; each extension gives a different series that agrees with $f$ on $[0,L]$. Two choices matter. The **odd extension** $f(-x) = -f(x)$ gives a pure sine series, the **even extension** a pure cosine series. With $x$ rescaled by $\pi/L$ and the projection rule (3.5), using $\|\sin(n\pi x/L)\|^2 = \|\cos(n\pi x/L)\|^2 = L/2$ on $[0, L]$ (the $\pi$ of Theorem 3.2 scaled by $L/\pi$ and halved for the half interval):

$$f \sim \sum_{n\ge1} b_n\sin\frac{n\pi x}{L},\ \ b_n = \frac{2}{L}\int_0^L f(x)\sin\frac{n\pi x}{L}\,dx; \qquad f \sim \frac{a_0}{2} + \sum_{n\ge1} a_n\cos\frac{n\pi x}{L},\ \ a_n = \frac2L\int_0^L f(x)\cos\frac{n\pi x}{L}\,dx.$$

The sine series is what a problem with $u = 0$ at both ends wants (heat with cold ends, string with fixed ends, §3.12), because every $\sin(n\pi x/L)$ already satisfies the boundary conditions. The cosine series is what zero-slope (insulated, free) ends want. Exercise B7 does the insulated case.

### 3.6 Bessel's inequality, Parseval's identity, Riemann–Lebesgue

**Theorem 3.6 (Bessel's inequality).** For $f \in PC$ and every $N$,
$$\frac{|a_0|^2}{2} + \sum_{n=1}^N\big(|a_n|^2 + |b_n|^2\big) = 2\sum_{|n|\le N}|c_n|^2 \;\le\; \frac{1}{\pi}\|f\|^2 = \frac1\pi\int_{-\pi}^{\pi}|f|^2\,dx. \tag{3.10}$$
Hence the series $\sum(|a_n|^2 + |b_n|^2)$ and $\sum_{n\in\mathbb{Z}}|c_n|^2$ converge.

▶ **Try it first.** Compute $\|S_Nf\|^2$ by expanding (3.6) in both slots and using Theorem 3.2 to kill the cross terms. Then use the identity at the end of §3.3.

**Derivation.**

1. $\|S_Nf\|^2 = \langle S_Nf, S_Nf\rangle$. Expand both slots using (3.6). Cross terms between different basis functions vanish by Theorem 3.2. The diagonal terms give $\left|\frac{a_0}{2}\right|^2\|1\|^2 + \sum_{n=1}^N\big(|a_n|^2\|\cos nx\|^2 + |b_n|^2\|\sin nx\|^2\big) = \frac{|a_0|^2}{4}\cdot 2\pi + \pi\sum_{n=1}^N(|a_n|^2 + |b_n|^2)$.
2. So $\frac1\pi\|S_Nf\|^2 = \frac{|a_0|^2}{2} + \sum_{n=1}^N(|a_n|^2 + |b_n|^2)$. The same computation in complex form with Theorem 3.3 gives $\|S_Nf\|^2 = 2\pi\sum_{|n|\le N}|c_n|^2$; the two agree because, by (3.7), $|c_n|^2 + |c_{-n}|^2 = \frac{|a_n|^2 + |b_n|^2}{2}$ for real $a_n, b_n$ (and in general, after a short expansion).
3. From §3.3, $\|f\|^2 = \|S_Nf\|^2 + \|f - S_Nf\|^2 \ge \|S_Nf\|^2$. Divide by $\pi$.
4. The partial sums of the nonnegative series $\sum(|a_n|^2 + |b_n|^2)$ are bounded above by $\frac1\pi\|f\|^2$, so the series converges. ∎

**Corollary 3.7 (Riemann–Lebesgue lemma, $PC$ version).** For $f \in PC$, $a_n \to 0$, $b_n \to 0$, and $c_n\to0$ as $|n|\to\infty$. Equivalently: for every $h \in PC$, $\int_{-\pi}^\pi h(x)\cos nx\,dx \to 0$ and $\int_{-\pi}^{\pi}h(x)\sin nx\,dx\to0$.

*Proof.* The terms of a convergent series tend to zero; apply this to $\sum(|a_n|^2 + |b_n|^2)$. ∎

This one-line corollary of projection is the engine of the convergence proof in §3.8. Its content: a rapidly oscillating sinusoid, multiplied by any fixed piecewise continuous function and integrated, averages out to nothing.

**Theorem 3.8 (Parseval's identity).** For $f\in PC$, equality holds in (3.10) as $N\to\infty$:
$$\frac1\pi\int_{-\pi}^\pi|f(x)|^2\,dx = \frac{|a_0|^2}{2} + \sum_{n=1}^{\infty}\big(|a_n|^2 + |b_n|^2\big), \qquad \frac{1}{2\pi}\int_{-\pi}^{\pi}|f(x)|^2\,dx = \sum_{n=-\infty}^{\infty}|c_n|^2. \tag{3.11}$$

*Proof, assuming mean-square convergence.* From step 3 of the Bessel derivation, $\frac1\pi\|f\|^2 - \frac1\pi\|S_Nf\|^2 = \frac1\pi\|f - S_Nf\|^2$. So Parseval holds if and only if $\|f - S_Nf\|\to0$: **if and only if $S_N f \to f$ in mean square.** That mean-square convergence holds for every $f\in PC$ (indeed for every $f$ with $\int|f|^2 < \infty$) is true, but its proof needs something this module does not establish: that trigonometric polynomials can approximate any such $f$ arbitrarily well in mean square (the basis is **complete**). The clean route is through Fejér's theorem (Exercise C3 builds the kernel) or the Weierstrass approximation theorem; Stein & Shakarchi, *Fourier Analysis*, ch. 3 §1, does it in four pages. We take completeness as given and use Parseval freely. ∎

Parseval is Pythagoras in infinitely many dimensions: the squared length of $f$ is the sum of the squared lengths of its components. Physically, $\frac{1}{2\pi}\int|f|^2$ is the mean power of a signal and $|c_n|^2$ is the power in harmonic $n$; ST-8's spectral density is exactly this. Example 4.5 uses it to compute $\sum 1/n^2$.

### 3.7 Smoothness means decay

▶ **Try it first.** Let $g$ be piecewise smooth, with jumps allowed. Integrate $c_n(g) = \frac{1}{2\pi}\int g\,e^{-inx}dx$ by parts on each smooth piece, differentiating $g$ and integrating $e^{-inx}$. Collect the boundary terms at each break. What do they add up to?

**Theorem 3.9 (A jump costs one power of $n$).** Let $g$ be $2\pi$-periodic and piecewise smooth, with jumps $J_j = g(x_j^+) - g(x_j^-)$ at the points $x_1, \dots, x_m$ of $(-\pi, \pi]$ (counting the seam $x = \pi$, where $g(\pi^+)$ means $g(-\pi^+)$). Then for $n \ne 0$,
$$c_n(g) = \frac{1}{in}\left[\frac{1}{2\pi}\sum_{j=1}^m J_j\,e^{-inx_j} \;+\; c_n(g')\right], \qquad\text{so}\qquad |c_n(g)| \le \frac{1}{|n|}\cdot\frac{\sum_j|J_j| + \int_{-\pi}^\pi|g'|\,dx}{2\pi}. \tag{3.12}$$
If $g$ is continuous (all $J_j = 0$), this reads $c_n(g') = in\,c_n(g)$.

**Derivation.**

1. On a smooth piece $(\alpha, \beta)$ between consecutive breaks, integrate by parts with $u = g$, $dv = e^{-inx}dx$, $v = \frac{e^{-inx}}{-in}$:
   $$\int_{\alpha}^{\beta} g\,e^{-inx}dx = \left[g(x)\frac{e^{-inx}}{-in}\right]_{\alpha^+}^{\beta^-} + \frac{1}{in}\int_\alpha^\beta g'(x)e^{-inx}dx.$$
2. Sum over the pieces. At an interior break $x_j$, the piece ending there contributes $g(x_j^-)\frac{e^{-inx_j}}{-in}$ and the piece starting there contributes $-g(x_j^+)\frac{e^{-inx_j}}{-in}$. Together: $\frac{e^{-inx_j}}{-in}\big(g(x_j^-) - g(x_j^+)\big) = \frac{J_j\,e^{-inx_j}}{in}$.
3. At the seam, the piece ending at $\pi$ gives $g(\pi^-)\frac{e^{-in\pi}}{-in}$ and the piece starting at $-\pi$ gives $-g(-\pi^+)\frac{e^{in\pi}}{-in}$. Since $e^{\pm in\pi} = (-1)^n$, together they give $\frac{(-1)^n}{in}\big(g(-\pi^+) - g(\pi^-)\big) = \frac{J_{\text{seam}}e^{-in\pi}}{in}$: the same form as step 2 with $x_j = \pi$. If $g$ is continuous and periodic this term is zero.
4. Divide by $2\pi$: $c_n(g) = \frac{1}{in}\Big[\frac{1}{2\pi}\sum_j J_je^{-inx_j} + \frac{1}{2\pi}\int g'e^{-inx}dx\Big]$, which is (3.12). The bound follows from $|e^{-inx_j}| = 1$ and $\left|\int g'e^{-inx}\right| \le \int|g'|$. ∎

■ **What just happened.** Integration by parts trades one power of $n$ for one derivative. Where $g$ is continuous the boundary terms cancel in pairs; where it jumps they don't, and each jump leaves behind a term of size $|J_j|/(2\pi|n|)$ that no amount of smoothness elsewhere can remove. Try it on the square wave: jumps $J = +2$ at $x = 0$ and $J = -2$ at the seam, $g' = 0$, so $c_n = \frac{1}{2\pi i n}\big(2 - 2(-1)^n\big) = -\frac{i(1 - (-1)^n)}{\pi n}$, which is $-\frac{2i}{\pi n}$ for odd $n$: Example 4.1 without computing an integral. On the sawtooth: one jump $J = -\pi - \pi = -2\pi$ at the seam and $g' = 1$ (whose $c_n$ vanishes for $n\ne0$), so $c_n = \frac{1}{in}\cdot\frac{-2\pi(-1)^n}{2\pi} = \frac{i(-1)^n}{n}$: Example 4.2.

**Theorem 3.10 (Decay rate).** Suppose $f, f', \dots, f^{(k)}$ are continuous and $2\pi$-periodic, and $f^{(k+1)}$ is piecewise smooth. Then
$$|c_n(f)| \le \frac{C}{|n|^{k+2}}, \qquad C = \frac{1}{2\pi}\Big(\sum_{\text{jumps of } f^{(k+1)}}|J_j| + \int_{-\pi}^{\pi}|f^{(k+2)}|\,dx\Big),$$
with the same bound (doubled) for $|a_n|$ and $|b_n|$. If you only know that $f^{(k+1)}$ is bounded and piecewise continuous, you get the weaker $|c_n| \le \max|f^{(k+1)}|/|n|^{k+1}$.

*Proof.* Apply the continuous case of (3.12) $k+1$ times: $c_n(f) = \frac{c_n(f')}{in} = \frac{c_n(f'')}{(in)^2} = \cdots = \frac{c_n(f^{(k+1)})}{(in)^{k+1}}$; each step is legitimate because the function being differentiated is continuous, periodic, and piecewise smooth. Then apply the jump case of (3.12) to $g = f^{(k+1)}$ for one more factor of $1/|n|$ with the stated constant. For the weaker bound stop one step earlier and use $|c_n(f^{(k+1)})| \le \frac{1}{2\pi}\int|f^{(k+1)}| \le \max|f^{(k+1)}|$. ∎

■ **What just happened.** Each derivative you can take continuously *across the periodic seam* buys a power of $1/n$, and the first derivative that jumps buys one more and then stops the process. Read backwards: the coefficients tell you how smooth the periodic extension is. The examples of §4 sort themselves out:

| Function on $[-\pi,\pi]$ | Periodic extension | Coefficients | Decay |
|---|---|---|---|
| square wave $\operatorname{sgn} x$ | jumps at $0$ and $\pm\pi$ | $b_n = \frac{4}{\pi n}$, $n$ odd | $1/n$ |
| sawtooth $x$ | jumps at $\pm\pi$ | $b_n = \frac{2(-1)^{n+1}}{n}$ | $1/n$ |
| triangle $\lvert x\rvert$ | continuous; $f'$ jumps at $0, \pm\pi$ | $a_n = -\frac{4}{\pi n^2}$, $n$ odd | $1/n^2$ |
| parabola $x^2$ | continuous; $f'$ jumps at $\pm\pi$ | $a_n = \frac{4(-1)^n}{n^2}$ | $1/n^2$ |
| cubic $x(\pi^2 - x^2)$ | $f, f'$ continuous; $f''$ jumps at $\pm\pi$ | $b_n = \frac{12(-1)^{n+1}}{n^3}$ | $1/n^3$ |

$x^2$ is a polynomial, as smooth as you like on the open interval, yet its *periodic extension* has a corner at $x = \pm\pi$ (slope $-2\pi$ meets slope $+2\pi$), and that corner is what the coefficients see. The cubic $x(\pi^2 - x^2)$ was chosen so that both it and its derivative match at the seam ($f(\pm\pi) = 0$, $f'(\pm\pi) = \pi^2 - 3\pi^2 = -2\pi^2$); its coefficients follow from (3.12) in two lines (since $f' = \pi^2 - 3x^2$, $c_n(f') = -3c_n(x^2) = -6(-1)^n/n^2$, so $c_n(f) = \frac{-6(-1)^n}{in\cdot n^2} = \frac{6i(-1)^n}{n^3}$ and $b_n = i(c_n - c_{-n}) = \frac{12(-1)^{n+1}}{n^3}$), and its $1/n^3$ is visible on the project's log-log plot. If $f$ is infinitely differentiable and periodic, $|c_n|$ decays faster than any power of $n$.

### 3.8 Convergence: the Dirichlet kernel and Dirichlet's theorem

So far we know $S_Nf$ is the best approximation in mean square and (taking completeness on faith) converges to $f$ in mean square. That says nothing about convergence at a *point*. For that we need $S_Nf(x)$ as an integral against $f$.

**Definition 3.11 (Dirichlet kernel).** $D_N(x) = \sum_{n=-N}^{N} e^{inx} = 1 + 2\sum_{n=1}^N\cos nx$.

**Theorem 3.12 (Closed form).** For $x$ not a multiple of $2\pi$,
$$D_N(x) = \frac{\sin\big((N+\tfrac12)x\big)}{\sin(x/2)}, \tag{3.13}$$
and $D_N(0) = 2N+1$. $D_N$ is even and $2\pi$-periodic, and $\frac{1}{2\pi}\int_{-\pi}^{\pi}D_N(x)\,dx = 1$.

▶ **Try it first.** $D_N$ is a geometric series in $r = e^{ix}$ with $2N+1$ terms starting at $r^{-N}$. Sum it, then multiply top and bottom by $e^{-ix/2}$ to make the result visibly real.

**Derivation.**

1. Factor out the first term: $D_N(x) = e^{-iNx}\sum_{k=0}^{2N}e^{ikx} = e^{-iNx}\,\dfrac{e^{i(2N+1)x} - 1}{e^{ix} - 1}$ (geometric sum, valid for $e^{ix}\ne1$).
2. Multiply through: $\dfrac{e^{i(N+1)x} - e^{-iNx}}{e^{ix} - 1}$. Multiply numerator and denominator by $e^{-ix/2}$: $\dfrac{e^{i(N+\frac12)x} - e^{-i(N+\frac12)x}}{e^{ix/2} - e^{-ix/2}}$.
3. Each difference is $2i\sin(\cdot)$: $\dfrac{2i\sin((N+\frac12)x)}{2i\sin(x/2)}$, which is (3.13).
4. At $x = 0$ every term of the sum is $1$, so $D_N(0) = 2N+1$ (also the limit of (3.13)). Evenness: $e^{inx} + e^{-inx} = 2\cos nx$. Integral: $\frac{1}{2\pi}\int D_N = \frac{1}{2\pi}\sum_{|n|\le N}\int e^{inx}dx = \frac{1}{2\pi}\cdot 2\pi = 1$, since only $n = 0$ survives (Theorem 3.3). ∎

**Theorem 3.13 (Partial sums are a convolution).** For $f \in PC$,
$$S_Nf(x) = \frac{1}{2\pi}\int_{-\pi}^{\pi} f(y)\,D_N(x - y)\,dy = \frac{1}{2\pi}\int_{-\pi}^{\pi} f(x - t)\,D_N(t)\,dt. \tag{3.14}$$

*Derivation.* Insert (3.8) into $S_Nf(x) = \sum_{|n|\le N}c_ne^{inx}$:
$$S_Nf(x) = \sum_{|n|\le N}\left(\frac{1}{2\pi}\int_{-\pi}^{\pi}f(y)e^{-iny}dy\right)e^{inx} = \frac{1}{2\pi}\int_{-\pi}^{\pi}f(y)\sum_{|n|\le N}e^{in(x-y)}\,dy = \frac{1}{2\pi}\int_{-\pi}^{\pi}f(y)D_N(x-y)\,dy.$$
Swapping a finite sum with an integral is always allowed. For the second form substitute $t = x - y$ and use the $2\pi$-periodicity of $f$ and $D_N$ to shift the interval of integration back to $[-\pi,\pi]$. ∎

■ **What just happened.** $S_Nf$ is the **convolution** of $f$ with $D_N$, normalized by $\frac{1}{2\pi}$. Truncating a Fourier series to $N$ terms is the same as smearing $f$ with the kernel $D_N$: a spike of height $2N+1$ and width about $2\pi/(2N+1)$ at the origin, surrounded by side lobes that change sign and decay only like $1/|x|$. Convergence $S_Nf(x)\to f(x)$ is the statement that $\frac{1}{2\pi}D_N$ acts more and more like a unit spike, a delta function, as $N\to\infty$. F3 makes convolution and the delta function central; this is where they first appear. The sign-changing side lobes are what make Fourier series delicate; a kernel that is positive (Fejér's, Exercise C3) is far better behaved.

**Theorem 3.14 (Dirichlet's theorem, piecewise smooth version).** Let $f$ be $2\pi$-periodic and piecewise smooth. Then for every $x$,
$$\lim_{N\to\infty}S_Nf(x) = \frac{f(x^+) + f(x^-)}{2}.$$
Where $f$ is continuous this is $f(x)$. At a jump, the series converges to the midpoint of the jump.

**Proof sketch.** Each step below is correct, but steps 3 and 4 lean on facts whose full verification is left to the references, and that is why this is a sketch.

1. **Split the kernel integral into two halves.** $D_N$ is even and $\frac{1}{2\pi}\int_{-\pi}^\pi D_N = 1$, so $\frac{1}{2\pi}\int_0^\pi D_N = \frac12$. From the second form of (3.14), splitting at $t = 0$ and substituting $t\to-t$ on the negative half,
   $$S_Nf(x) = \frac{1}{2\pi}\int_0^\pi\big[f(x+t) + f(x-t)\big]D_N(t)\,dt.$$
   Subtract $\frac{f(x^+) + f(x^-)}{2} = \frac{1}{2\pi}\int_0^\pi\big[f(x^+) + f(x^-)\big]D_N(t)\,dt$ and insert (3.13):
   $$S_Nf(x) - \frac{f(x^+)+f(x^-)}{2} = \frac{1}{2\pi}\int_0^\pi g_x(t)\,\sin\!\big((N+\tfrac12)t\big)\,dt, \qquad g_x(t) = \frac{\big[f(x+t) - f(x^+)\big] + \big[f(x-t) - f(x^-)\big]}{\sin(t/2)}.$$
2. **The goal.** Show the right side tends to $0$. If $g_x$ were in $PC[0,\pi]$, Riemann–Lebesgue would do it. The only danger is $t\to0$, where the denominator $\sin(t/2)\to0$.
3. **$g_x$ is bounded near $t = 0$.** This is where piecewise smoothness is used. Write
   $$g_x(t) = \left[\frac{f(x+t)-f(x^+)}{t} + \frac{f(x-t) - f(x^-)}{t}\right]\cdot\frac{t}{\sin(t/2)}.$$
   As $t\to0^+$ the first fraction tends to the right-hand derivative $f'(x^+)$ and the second to $-f'(x^-)$, both of which exist because $f'\in PC$; and $\frac{t}{\sin(t/2)}\to2$. So $g_x$ has a finite limit at $0^+$, and away from $0$ it is a quotient of $PC$ functions with nonvanishing denominator. Hence $g_x\in PC[0,\pi]$. *What this hides:* for $f$ merely continuous there need be no one-sided derivatives, $g_x$ need not be integrable, and the theorem is false (du Bois-Reymond's example).
4. **Riemann–Lebesgue at a half-integer frequency.** Corollary 3.7 was stated for $\int h\cos nx$ and $\int h\sin nx$ with integer $n$ over $[-\pi,\pi]$. Expand $\sin((N+\frac12)t) = \sin Nt\cos\frac{t}{2} + \cos Nt\sin\frac t2$; then
   $$\int_0^\pi g_x(t)\sin((N+\tfrac12)t)\,dt = \int_{-\pi}^{\pi}h_1(t)\sin Nt\,dt + \int_{-\pi}^{\pi}h_2(t)\cos Nt\,dt,$$
   with $h_1 = g_x\cos\frac t2$ and $h_2 = g_x\sin\frac t2$ on $(0,\pi]$, both extended by $0$ on $[-\pi,0]$. Both are in $PC[-\pi,\pi]$, so both integrals are $\pi$ times Fourier coefficients of $PC$ functions and tend to $0$ by Corollary 3.7. *What this hides:* nothing, for piecewise smooth $f$. But Dirichlet's original hypothesis (finitely many extrema) and the modern bounded-variation version need Riemann–Lebesgue for merely *integrable* $h$, which requires approximating $h$ by step functions; our Bessel-based proof only covers square-integrable $h$. Stein & Shakarchi ch. 3 §2 has the general local result. ∎

■ **What just happened.** The convolution formula localizes the question: whether $S_Nf(x)$ converges depends only on how $f$ behaves near $x$ (the **Riemann localization principle**), because away from $t = 0$ the kernel's oscillation kills everything by Riemann–Lebesgue. Near $t = 0$ the kernel is a spike of total mass $1$ and picks out the average of the two one-sided limits. The whole difficulty of the full theorem is controlling the quotient $g_x$ near $0$; piecewise smoothness makes it a non-issue.

### 3.9 The Gibbs phenomenon

Dirichlet's theorem is pointwise. It does not say the convergence is uniform, and near a jump it cannot be: each $S_Nf$ is continuous, so if $S_Nf\to f$ uniformly then $f$ would be continuous. What actually happens near a jump is specific and permanent. (It was first noticed by Henry Wilbraham in 1848, forgotten, and rediscovered in an 1899 letter to *Nature* by J. Willard Gibbs, whose name it carries.)

**Theorem 3.15 (Gibbs overshoot for the square wave).** Let $f(x) = \operatorname{sgn}(x)$ on $(-\pi,\pi)$, so $f \sim \frac4\pi\sum_{k\text{ odd}}\frac{\sin kx}{k}$ (Example 4.1). For $N = 2M - 1$, so that $S_N$ contains the odd harmonics $1, 3, \dots, 2M-1$, the first maximum of $S_Nf$ to the right of $0$ is at $x_M = \frac{\pi}{2M}$, and
$$\lim_{M\to\infty}S_{2M-1}f(x_M) = \frac{2}{\pi}\int_0^\pi\frac{\sin u}{u}\,du = \frac{2}{\pi}\operatorname{Si}(\pi) \approx 1.17898. \tag{3.15}$$
The limiting overshoot above the value $1$ is $\approx 0.17898$, which is $\frac{\operatorname{Si}(\pi)}{\pi} - \frac12 \approx 0.08949$, about $9\%$, of the jump $2$. It does not shrink as $N\to\infty$: the location $x_M\to0$, the height does not.

▶ **Try it first.** Differentiate $S_Nf$ termwise (a finite sum, so this is legal), sum the resulting finite cosine series in closed form, find its first positive zero, and write $S_Nf(x_M) = \int_0^{x_M}S_Nf'(x)\,dx$. Then substitute $u = 2Mx$.

**Derivation.**

1. $S_Nf(x) = \frac4\pi\sum_{m=1}^{M}\frac{\sin((2m-1)x)}{2m-1}$, so $S_Nf'(x) = \frac4\pi\sum_{m=1}^M\cos((2m-1)x)$.
2. Sum the odd cosines: multiply by $2\sin x$ and use $2\sin x\cos((2m-1)x) = \sin(2mx) - \sin((2m-2)x)$ (product-to-sum, §3.2). The sum telescopes to $\sin(2Mx) - \sin 0 = \sin 2Mx$. So $\sum_{m=1}^M\cos((2m-1)x) = \dfrac{\sin 2Mx}{2\sin x}$ and
   $$S_Nf'(x) = \frac{2}{\pi}\,\frac{\sin 2Mx}{\sin x}.$$
3. For $0 < x < \pi$ the denominator is positive, so the sign of $S_Nf'$ is the sign of $\sin 2Mx$: positive on $(0, \frac{\pi}{2M})$, negative just after. The first maximum is at $x_M = \frac{\pi}{2M}$.
4. Since $S_Nf(0) = 0$, $S_Nf(x_M) = \int_0^{x_M}\frac{2}{\pi}\frac{\sin 2Mx}{\sin x}\,dx$. Substitute $u = 2Mx$, $dx = du/2M$:
   $$S_Nf(x_M) = \frac{2}{\pi}\int_0^{\pi}\frac{\sin u}{2M\sin\!\big(\frac{u}{2M}\big)}\,du.$$
5. As $M\to\infty$, $2M\sin(u/2M)\to u$ uniformly for $u\in[0,\pi]$ (since $\sin\theta = \theta - \theta^3/6 + \cdots$), so the integrand converges uniformly to $\frac{\sin u}{u}$, which is bounded, and the limit is $\frac2\pi\int_0^\pi\frac{\sin u}{u}du = \frac2\pi\operatorname{Si}(\pi)$. Numerically $\operatorname{Si}(\pi) = 1.8519370\ldots$ (`scipy.special.sici`), so the limit is $1.1789797\ldots$ ∎

■ **What just happened.** Near the jump, $S_Nf$ is a copy of one fixed function, $\frac2\pi\operatorname{Si}(u)$, compressed horizontally by $2M$ but not vertically. Horizontal compression toward the jump is why Dirichlet's theorem is still true (at any fixed $x > 0$, eventually $x \gg x_M$ and the ripples there have died down); the fixed vertical scale is why the overshoot is permanent. Because $2M\sin(u/2M) < u$, the finite-$M$ value is slightly *above* the limit, so the project measures overshoots a little over $0.1790$ that decrease toward it. By the same mechanism (Example 4.8 and Exercise C1), any piecewise smooth $f$ overshoots by $\approx 0.0895\times$(size of the jump) on each side of each jump; for the sawtooth, with its jump of $2\pi$, that is $\approx 0.562$. Gibbs is **not** a numerical artifact and not a defect of your code: it is a theorem about the partial sums. It is also why $\|f - S_Nf\|\to0$ can hold while $\max|f - S_Nf|$ stays at $0.18$: the region where $S_Nf$ is wrong has width $\sim1/N$, so its contribution to $\int|f-S_Nf|^2$ is $\sim(0.18)^2/N\to0$.

### 3.10 Uniform convergence

**Theorem 3.16 (Uniform convergence).** If $\sum_{n\in\mathbb{Z}}|c_n| < \infty$ (equivalently $\sum_n(|a_n| + |b_n|) < \infty$), then $S_Nf$ converges absolutely and uniformly on $\mathbb{R}$ to a continuous function, and for every $N$,
$$\max_x|f(x) - S_Nf(x)| \le \sum_{|n| > N}|c_n|.$$
In particular, if $f$ is continuous, $2\pi$-periodic, and piecewise smooth (so $|c_n|\le C/n^2$ by Theorem 3.10 with $k = 0$), then $S_Nf\to f$ uniformly, with error at most $\sum_{|n|>N}C/n^2 < 2C/N$.

*Proof.* $|c_ne^{inx}| = |c_n|$ for all $x$ and $\sum|c_n| < \infty$, so the Weierstrass M-test gives absolute and uniform convergence of $\sum c_ne^{inx}$ to some limit $g$, continuous as a uniform limit of continuous functions. The tail bound is $|g(x) - S_Nf(x)| = \big|\sum_{|n|>N}c_ne^{inx}\big| \le \sum_{|n|>N}|c_n|$. When $f$ is continuous and piecewise smooth, Dirichlet's theorem says $S_Nf(x)\to f(x)$ at every $x$, so $g = f$. Finally $\sum_{n > N}n^{-2} < \int_N^\infty u^{-2}du = 1/N$, doubled for $\pm n$. ∎

So for $x^2$ and $|x|$ the convergence is uniform, and the project finds how many terms it takes to get within $0.01$ everywhere ($N = 400$ and $N = 63$; the tail bounds $4/N$ and $2/(\pi N)$ predict $400$ and $64$). For the square wave and sawtooth $\sum|b_n| = \infty$ (harmonic series), and the convergence is not uniform (Gibbs).

### 3.11 Termwise integration and differentiation

**Theorem 3.17 (Integration termwise, always allowed).** Let $f\in PC$ with coefficients $a_n, b_n$. Then for every $x$,
$$\int_0^x f(t)\,dt = \frac{a_0 x}{2} + \sum_{n=1}^{\infty}\left[\frac{a_n\sin nx}{n} + \frac{b_n(1 - \cos nx)}{n}\right], \tag{3.16}$$
and the series converges uniformly. You may integrate a Fourier series term by term whether or not it converges, and the result converges to the integral of $f$.

▶ **Try it first.** The function $F(x) = \int_0^x f(t)dt - \frac{a_0x}{2}$ is continuous. Check that it is $2\pi$-periodic (compute $F(x + 2\pi) - F(x)$), then find its coefficients by integration by parts using $F' = f - a_0/2$ wherever $f$ is continuous.

**Derivation.**

1. $F(x+2\pi) - F(x) = \int_x^{x+2\pi}f(t)\,dt - \pi a_0 = \int_{-\pi}^{\pi}f\,dt - \pi a_0 = \pi a_0 - \pi a_0 = 0$, using periodicity of $f$ to shift the interval and $a_0 = \frac1\pi\int f$. So $F$ is continuous, periodic, and piecewise smooth, with $F' = f - \frac{a_0}{2}\in PC$.
2. By Theorem 3.14 applied to the continuous $F$, its Fourier series converges to $F(x)$ at every $x$; by Theorem 3.16 (since $|c_n(F)|\le C/n^2$ by Theorem 3.10) the convergence is uniform.
3. Coefficients of $F$ for $n\ge1$, by parts, boundary terms vanishing because $F$ is periodic:
   $$a_n(F) = \frac1\pi\int_{-\pi}^\pi F\cos nx\,dx = \frac1\pi\left[F\frac{\sin nx}{n}\right]_{-\pi}^{\pi} - \frac{1}{n\pi}\int_{-\pi}^\pi F'\sin nx\,dx = -\frac{b_n(F')}{n} = -\frac{b_n}{n},$$
   $$b_n(F) = \frac1\pi\int_{-\pi}^\pi F\sin nx\,dx = \frac1\pi\left[-F\frac{\cos nx}{n}\right]_{-\pi}^\pi + \frac{1}{n\pi}\int_{-\pi}^\pi F'\cos nx\,dx = \frac{a_n(F')}{n} = \frac{a_n}{n},$$
   where $a_n(F') = a_n$ and $b_n(F') = b_n$ for $n\ge1$ because the constant $\frac{a_0}{2}$ in $F' = f - \frac{a_0}{2}$ contributes only to the $n = 0$ coefficient.
4. The constant term $\frac{a_0(F)}{2}$: rather than integrate $F$, evaluate the convergent series at $x = 0$, where $F(0) = 0$: $0 = \frac{a_0(F)}{2} + \sum_n a_n(F)$, so $\frac{a_0(F)}{2} = \sum_n\frac{b_n}{n}$.
5. Assemble: $F(x) = \sum_n\frac{b_n}{n} + \sum_n\left(-\frac{b_n}{n}\cos nx + \frac{a_n}{n}\sin nx\right)$, and $\int_0^xf = F(x) + \frac{a_0x}{2}$ gives (3.16). ∎

■ **What just happened.** Integration divides $c_n$ by $in$, which improves convergence, so the integrated series always converges, and the proof needed only Dirichlet's theorem for the *continuous* function $F$. The term $\frac{a_0x}{2}$ is not periodic; the integral of a periodic function is periodic only if its mean is zero. Example: integrating the square-wave series $\frac4\pi\sum_{k\text{ odd}}\frac{\sin kx}{k}$ from $0$ to $x$ gives $\frac4\pi\sum_{k\text{ odd}}\frac{1 - \cos kx}{k^2} = \frac{4}{\pi}\cdot\frac{\pi^2}{8} - \frac4\pi\sum_{k \text{ odd}}\frac{\cos kx}{k^2}$, which is the triangle series of Example 4.3, as it must be since $\int_0^x\operatorname{sgn}(t)dt = |x|$.

**Theorem 3.18 (Differentiation termwise, when allowed).** If $f$ is continuous, $2\pi$-periodic, and piecewise smooth, and $f'$ is itself piecewise smooth, then the termwise derivative of the Fourier series of $f$ is the Fourier series of $f'$: $c_n(f') = in\,c_n(f)$, i.e. $a_n(f') = n\,b_n(f)$ and $b_n(f') = -n\,a_n(f)$. By Theorem 3.14 it converges to $\frac{f'(x^+) + f'(x^-)}{2}$ at every $x$.

*Proof.* The coefficient relation is the continuous case of (3.12), whose only hypothesis was that $f$ is continuous and periodic with $f'\in PC$. Convergence of the differentiated series is Dirichlet's theorem applied to the piecewise smooth $f'$. ∎

The hypothesis that $f$ is **continuous as a periodic function** is the one people forget. Differentiate the sawtooth series $x = 2\sum(-1)^{n+1}\frac{\sin nx}{n}$ termwise and you get $2\sum(-1)^{n+1}\cos nx$, whose terms do not even tend to zero: it diverges everywhere. Nothing went wrong with the theorem; the sawtooth has jumps at $\pm\pi$, so the jump term of (3.12) is present, and it is exactly that term, $\frac{1}{2\pi}J e^{-in\pi}$ with $J = -2\pi$, which when multiplied by $in$ refuses to decay. The divergent series is trying to tell you that the derivative of the periodic sawtooth is $1$ minus a train of infinitely tall spikes at the jumps. F3 makes this precise: $2\sum(-1)^{n+1}\cos nx$ is the Fourier series of $1 - 2\pi\sum_k\delta\big(x - (2k+1)\pi\big)$, a "delta comb," and the series is meaningful as a distribution.

### 3.12 Two PDEs: heat and the vibrating string

Both of §1's problems are now solvable. The method is **separation of variables**, and this module is what completes it.

**The heat equation.** A rod of length $L$, thermal diffusivity $k$ (units m²/s), ends held at temperature $0$, initial temperature $\phi(x)$:
$$u_t = k\,u_{xx}\quad (0 < x < L,\ t > 0), \qquad u(0,t) = u(L,t) = 0, \qquad u(x, 0) = \phi(x).$$

▶ **Try it first.** Look for solutions $u(x,t) = X(x)T(t)$. Substitute, divide by $kXT$, and argue that a function of $x$ alone equal to a function of $t$ alone must be constant. Solve the two ODEs with the boundary conditions. Then use linearity and §3.5.

**Derivation.**

1. Substitute $u = XT$: $XT' = kX''T$. Divide by $kXT$: $\dfrac{T'}{kT} = \dfrac{X''}{X}$. The left side depends only on $t$, the right only on $x$; they can be equal for all $x$ and $t$ only if both equal the same constant, call it $-\lambda$.
2. The $x$-problem: $X'' + \lambda X = 0$ with $X(0) = X(L) = 0$ (from $u(0,t) = X(0)T(t) = 0$ for all $t$ and $T\not\equiv0$). If $\lambda < 0$, $X = Ae^{\sqrt{-\lambda}x} + Be^{-\sqrt{-\lambda}x}$; $X(0) = 0$ gives $B = -A$, then $X(L) = 2A\sinh(\sqrt{-\lambda}L) = 0$ forces $A = 0$. If $\lambda = 0$, $X = Ax + B$ and the boundary conditions force $A = B = 0$. If $\lambda > 0$, $X = A\cos\sqrt\lambda x + B\sin\sqrt\lambda x$; $X(0) = 0$ gives $A = 0$, and $X(L) = B\sin(\sqrt\lambda L) = 0$ with $B\ne0$ requires $\sqrt\lambda L = n\pi$. So
   $$\lambda_n = \left(\frac{n\pi}{L}\right)^2, \qquad X_n(x) = \sin\frac{n\pi x}{L}, \qquad n = 1, 2, 3, \dots$$
3. The $t$-problem: $T' = -k\lambda_nT$, so $T_n(t) = e^{-k(n\pi/L)^2t}$.
4. Each $u_n = X_nT_n$ solves the PDE and the boundary conditions. By linearity so does any finite sum, and (with the convergence theorems of this module to justify the limit) the infinite sum
   $$u(x,t) = \sum_{n=1}^\infty b_n\sin\frac{n\pi x}{L}\,e^{-k(n\pi/L)^2 t}. \tag{3.17}$$
5. The initial condition requires $\phi(x) = \sum b_n\sin\frac{n\pi x}{L}$ on $(0, L)$: the **sine series** of $\phi$ (§3.5), so $b_n = \frac2L\int_0^L\phi(x)\sin\frac{n\pi x}{L}dx$. This is the step that was missing in 1753 and that Fourier supplied. ∎

■ **What just happened.** Mode $n$ decays with time constant $\tau_n = \dfrac{L^2}{k\,n^2\pi^2}$: the second harmonic dies four times faster than the first, the tenth a hundred times faster. After even a short time the high modes are gone and the profile is dominated by $\sin(\pi x/L)$. This is **smoothing**: by Theorem 3.10 read backwards, the factor $e^{-kn^2\pi^2t/L^2}$ makes the coefficients decay faster than any power of $n$ for any $t > 0$, so $u(\cdot, t)$ is infinitely differentiable even if $\phi$ had jumps. Heat flow erases detail, highest frequencies first. Example 4.6 does this for a step.

**The wave equation.** A string of length $L$, wave speed $c = \sqrt{\text{tension}/\text{linear density}}$, fixed ends, initial shape $\phi$ and initial velocity $\psi$:
$$u_{tt} = c^2u_{xx}, \qquad u(0,t) = u(L,t) = 0, \qquad u(x,0) = \phi(x),\quad u_t(x,0) = \psi(x).$$
Separation gives the same $X_n = \sin\frac{n\pi x}{L}$ (the $x$-problem is identical) and now $T_n'' = -c^2\lambda_nT_n$, so $T_n = A_n\cos\omega_nt + B_n\sin\omega_nt$ with
$$\omega_n = \frac{n\pi c}{L}, \qquad u(x,t) = \sum_{n=1}^\infty\sin\frac{n\pi x}{L}\big(A_n\cos\omega_nt + B_n\sin\omega_nt\big), \tag{3.18}$$
where the $A_n$ are the sine coefficients of $\phi$ and the $\omega_nB_n$ those of $\psi$. Each term is a **standing wave** or **normal mode**: a fixed shape $\sin(n\pi x/L)$ with $n-1$ nodes, oscillating at $\omega_n$. The frequencies are integer multiples of the fundamental $\omega_1 = \pi c/L$, which is why a string sounds harmonious and a drum (whose modes are Bessel functions with incommensurate zeros, EM-3) does not. Unlike heat, nothing decays: energy put into mode $n$ stays there (no damping in this model). Example 4.7 plucks the string.

*CM-6 preview.* Replace the string by $N$ beads of mass $m$ on a massless string with spacing $h$. Newton's law gives $m\ddot y_j = \frac{\tau}{h}(y_{j+1} - 2y_j + y_{j-1})$, a system $\ddot{\mathbf{y}} = -\frac{\tau}{mh}K\mathbf{y}$ with $K$ the tridiagonal second-difference matrix. Its eigenvectors are the sampled sines $\sin\frac{n\pi j}{N+1}$ (CM-6), with eigenvalues $2 - 2\cos\frac{n\pi}{N+1} \approx \left(\frac{n\pi}{N+1}\right)^2$ for small $n$. As $N\to\infty$ with $Nh = L$ fixed, the eigenvectors become $\sin(n\pi x/L)$ and the eigenfrequencies become $\omega_n$: the continuous normal modes are the limit of the coupled-oscillator eigenproblem, and orthogonality of the modes is orthogonality of eigenvectors of a symmetric matrix. §5 makes the same statement for the continuous operator directly.

## 4. Worked examples

All integrals use (3.4) with the even/odd shortcuts of §3.5. Every number below is reproduced in `project.py`, symbolically (SymPy) and numerically.

### Example 4.1: The square wave

$f(x) = \operatorname{sgn}(x)$: $-1$ on $(-\pi, 0)$, $+1$ on $(0, \pi)$, extended periodically. Jumps of size $2$ at $x = 0$ and at $x = \pm\pi$.

$f$ is odd, so $a_n = 0$ for all $n$ (including $a_0$: the mean is zero). For $b_n$ the integrand $f(x)\sin nx$ is even, and we **split at the break** $x = 0$:
$$b_n = \frac1\pi\int_{-\pi}^\pi\operatorname{sgn}(x)\sin nx\,dx = \frac{2}{\pi}\int_0^\pi\sin nx\,dx = \frac2\pi\left[-\frac{\cos nx}{n}\right]_0^\pi = \frac{2}{\pi}\cdot\frac{1 - \cos n\pi}{n} = \frac{2\big(1 - (-1)^n\big)}{\pi n}.$$
For even $n$ this is $0$; for odd $n$ it is $\frac{4}{\pi n}$:
$$\operatorname{sgn}(x) \sim \frac4\pi\left(\sin x + \frac{\sin 3x}{3} + \frac{\sin 5x}{5} + \cdots\right), \qquad b_1 = \frac4\pi \approx 1.2732,\ \ b_3 \approx 0.4244,\ \ b_5\approx 0.2546.$$
Complex form: $c_n = -\frac{ib_n}{2} = -\frac{2i}{\pi n}$ for odd $n$, $0$ for even $n$. Decay $1/n$, as Theorem 3.9 predicts for a jump (and as its jump formula gives directly, §3.7). Dirichlet check at $x = \pi/2$: every term is $\frac{4}{\pi}\frac{\sin(k\pi/2)}{k} = \frac4\pi\frac{(-1)^{(k-1)/2}}{k}$, so $1 = \frac4\pi\left(1 - \frac13 + \frac15 - \cdots\right)$, Leibniz's series for $\pi/4$. ✓ At $x = 0$ every term is $0$ and the series gives $0 = \frac{f(0^+) + f(0^-)}{2}$. ✓

### Example 4.2: The sawtooth

$f(x) = x$ on $(-\pi, \pi)$, extended periodically: jumps of size $2\pi$ at $x = \pm\pi$. Odd, so $a_n = 0$.
$$b_n = \frac2\pi\int_0^\pi x\sin nx\,dx.$$
By parts with $u = x$, $dv = \sin nx\,dx$: $\int_0^\pi x\sin nx\,dx = \left[-\frac{x\cos nx}{n}\right]_0^\pi + \frac1n\int_0^\pi\cos nx\,dx = -\frac{\pi\cos n\pi}{n} + \frac{1}{n}\left[\frac{\sin nx}{n}\right]_0^\pi = -\frac{\pi(-1)^n}{n} + 0$. So
$$b_n = \frac2\pi\cdot\left(-\frac{\pi(-1)^n}{n}\right) = \frac{2(-1)^{n+1}}{n}, \qquad x \sim 2\left(\sin x - \frac{\sin 2x}{2} + \frac{\sin 3x}{3} - \frac{\sin 4x}{4} + \cdots\right).$$
$b_1 = 2$, $b_2 = -1$, $b_3 = \frac23$, $b_4 = -\frac12$. Complex: $c_n = -\frac{ib_n}{2} = \frac{i(-1)^n}{n}$ for $n\ne0$, $c_0 = 0$. Decay $1/n$. At $x = \pi$ the series gives $0$, the average of $\pi$ and $-\pi$. ✓ At $x = \pi/2$: $\frac\pi2 = 2\left(1 - 0 - \frac13 + 0 + \frac15 - \cdots\right)$, Leibniz again. ✓

### Example 4.3: The triangle wave

$f(x) = |x|$ on $[-\pi,\pi]$: continuous as a periodic function ($f(\pi) = f(-\pi) = \pi$), with corners at $0$ and $\pm\pi$. Even, so $b_n = 0$.
$$a_0 = \frac2\pi\int_0^\pi x\,dx = \frac2\pi\cdot\frac{\pi^2}{2} = \pi, \qquad\text{so the constant term is } \frac{a_0}{2} = \frac\pi2, \text{ the mean of } |x|.$$
For $n\ge1$, by parts with $u = x$, $dv = \cos nx\,dx$:
$$\int_0^\pi x\cos nx\,dx = \left[\frac{x\sin nx}{n}\right]_0^\pi - \frac1n\int_0^\pi\sin nx\,dx = 0 - \frac1n\cdot\frac{1 - (-1)^n}{n} = -\frac{1 - (-1)^n}{n^2},$$
$$a_n = \frac2\pi\cdot\left(-\frac{1-(-1)^n}{n^2}\right) = -\frac{2\big(1 - (-1)^n\big)}{\pi n^2} = \begin{cases}-\dfrac{4}{\pi n^2} & n\text{ odd}\\[4pt] 0 & n \text{ even}\end{cases}$$
$$|x| \sim \frac\pi2 - \frac4\pi\left(\cos x + \frac{\cos 3x}{9} + \frac{\cos 5x}{25} + \cdots\right), \qquad a_1 \approx -1.2732,\ \ a_3\approx-0.1415,\ \ a_5\approx -0.0509.$$
Decay $1/n^2$, one order better than the square wave, because $|x|$ is continuous. Indeed $\frac{d}{dx}|x| = \operatorname{sgn}x$, and Theorem 3.18 gives $b_n(f') = -na_n(f) = \frac{4}{\pi n}$: differentiating Example 4.3 termwise gives Example 4.1. ✓ At $x = 0$: $0 = \frac\pi2 - \frac4\pi\sum_{k\text{ odd}}\frac1{k^2}$, so $\sum_{k\text{ odd}}\frac{1}{k^2} = \frac{\pi^2}{8}$. ✓ The convergence is uniform by Theorem 3.16, so there is no Gibbs overshoot; the project shows the partial sums hugging the corners.

### Example 4.4: The parabola $x^2$

$f(x) = x^2$ on $[-\pi,\pi]$: continuous as a periodic function, corner at $\pm\pi$ (slope jumps from $2\pi$ to $-2\pi$). Even, so $b_n = 0$.
$$a_0 = \frac2\pi\int_0^\pi x^2dx = \frac2\pi\cdot\frac{\pi^3}{3} = \frac{2\pi^2}{3}, \qquad \frac{a_0}{2} = \frac{\pi^2}{3}.$$
For $n\ge1$, by parts with $u = x^2$, $dv = \cos nx\,dx$, then reuse $\int_0^\pi x\sin nx\,dx = -\frac{\pi(-1)^n}{n}$ from Example 4.2:
$$\int_0^\pi x^2\cos nx\,dx = \left[\frac{x^2\sin nx}{n}\right]_0^\pi - \frac2n\int_0^\pi x\sin nx\,dx = 0 - \frac2n\left(-\frac{\pi(-1)^n}{n}\right) = \frac{2\pi(-1)^n}{n^2},$$
$$a_n = \frac2\pi\cdot\frac{2\pi(-1)^n}{n^2} = \frac{4(-1)^n}{n^2}, \qquad x^2 \sim \frac{\pi^2}{3} + 4\sum_{n=1}^\infty\frac{(-1)^n}{n^2}\cos nx = \frac{\pi^2}{3} - 4\cos x + \cos 2x - \frac49\cos 3x + \cdots$$
Complex: $c_0 = \frac{\pi^2}{3}$, $c_n = \frac{2(-1)^n}{n^2}$. Decay $1/n^2$. Two free checks, both legitimate because $x^2$ is continuous and piecewise smooth so Theorem 3.14 gives equality everywhere: at $x = \pi$, $\pi^2 = \frac{\pi^2}{3} + 4\sum\frac{1}{n^2}$, so $\sum_{n\ge1}\frac1{n^2} = \frac{\pi^2}{6}$; at $x = 0$, $0 = \frac{\pi^2}{3} + 4\sum\frac{(-1)^n}{n^2}$, so $\sum\frac{(-1)^{n+1}}{n^2} = \frac{\pi^2}{12}$.

### Example 4.5: Basel via Parseval, and $\sum 1/n^4$

Parseval (3.11) for the sawtooth of Example 4.2. Left side: $\frac1\pi\int_{-\pi}^\pi x^2\,dx = \frac1\pi\cdot\frac{2\pi^3}{3} = \frac{2\pi^2}{3}$. Right side: $a_n = 0$ and $b_n^2 = \frac{4}{n^2}$, so $\sum b_n^2 = 4\sum\frac1{n^2}$. Therefore
$$4\sum_{n=1}^{\infty}\frac1{n^2} = \frac{2\pi^2}{3} \quad\Longrightarrow\quad \sum_{n=1}^\infty\frac1{n^2} = \frac{\pi^2}{6} \approx 1.644934.$$
This is Euler's 1735 answer to the Basel problem, here a one-line consequence of Pythagoras in function space. Partial sums: $\sum_{n\le10} = 1.549768$, $\sum_{n\le100} = 1.634984$, $\sum_{n\le1000} = 1.643935$; the error is about $1/N$.

Now the parabola of Example 4.4. Left side: $\frac1\pi\int_{-\pi}^\pi x^4\,dx = \frac1\pi\cdot\frac{2\pi^5}{5} = \frac{2\pi^4}{5}$. Right side: $\frac{a_0^2}{2} + \sum a_n^2 = \frac12\left(\frac{2\pi^2}{3}\right)^2 + \sum\frac{16}{n^4} = \frac{2\pi^4}{9} + 16\sum\frac1{n^4}$. So
$$16\sum_{n=1}^\infty\frac{1}{n^4} = \frac{2\pi^4}{5} - \frac{2\pi^4}{9} = 2\pi^4\cdot\frac{9 - 5}{45} = \frac{8\pi^4}{45} \quad\Longrightarrow\quad \sum_{n=1}^\infty\frac1{n^4} = \frac{\pi^4}{90}\approx 1.082323.$$
Each new smooth function hands you a new zeta value; Exercise B3 gets $\sum_{\text{odd}}1/n^4 = \pi^4/96$ from the triangle and $\sum_{\text{odd}} 1/n^6 = \pi^6/960$ from $x(\pi-x)$.

### Example 4.6: Heat in a rod with a step initial temperature

A rod of length $L$ has its left half at temperature $u_0$ and its right half at $0$ at $t = 0$; both ends are held at $0$ forever. So $\phi(x) = u_0$ for $0 < x < L/2$ and $\phi(x) = 0$ for $L/2 < x < L$. By (3.17) we need the sine coefficients of $\phi$, splitting the integral at the break $x = L/2$:
$$b_n = \frac2L\int_0^{L/2}u_0\sin\frac{n\pi x}{L}\,dx = \frac{2u_0}{L}\left[-\frac{L}{n\pi}\cos\frac{n\pi x}{L}\right]_0^{L/2} = \frac{2u_0}{n\pi}\left(1 - \cos\frac{n\pi}{2}\right).$$
$\cos\frac{n\pi}{2}$ cycles through $0, -1, 0, 1$ for $n = 1, 2, 3, 4$, so
$$b_1 = \frac{2u_0}{\pi} \approx 0.6366u_0,\quad b_2 = \frac{2u_0}{\pi}\approx0.6366u_0,\quad b_3 = \frac{2u_0}{3\pi}\approx0.2122u_0,\quad b_4 = 0,\quad b_5 = \frac{2u_0}{5\pi}\approx 0.1273u_0,\quad b_6 = \frac{2u_0}{3\pi}\approx0.2122u_0.$$
Unlike Example 4.1, even modes appear: $\phi$ is not symmetric about the center, and $\sin(2\pi x/L)$ is exactly the mode that is positive on the left half and negative on the right. The first three terms of the solution:
$$u(x,t) \approx \frac{2u_0}{\pi}\left[\sin\frac{\pi x}{L}\,e^{-\pi^2kt/L^2} + \sin\frac{2\pi x}{L}\,e^{-4\pi^2kt/L^2} + \frac13\sin\frac{3\pi x}{L}\,e^{-9\pi^2kt/L^2}\right].$$
Time scales: $\tau_n = \dfrac{L^2}{n^2\pi^2k}$. For a copper rod with $L = 10\,\text{cm}$ and $k = 1.1\times10^{-4}\,\text{m}^2/\text{s}$: $\tau_1 = \dfrac{0.01}{\pi^2\cdot1.1\times10^{-4}} = 9.21\,\text{s}$, $\tau_2 = 2.30\,\text{s}$, $\tau_3 = 1.02\,\text{s}$. At $t = \tau_1$, mode 1 has decayed to $e^{-1} = 37\%$ of its initial amplitude, mode 2 to $e^{-4} = 1.8\%$, mode 3 to $e^{-9} = 0.012\%$. After about ten seconds the step has become a single sine arch of amplitude $\approx 0.6366u_0e^{-1} \approx 0.23u_0$, and the memory of where the step was survives only in the dying trace of mode 2. The project overlays the truncated series on a finite-difference solution at several times; they agree to about $10^{-3}u_0$ at the earliest time and better later.

### Example 4.7: A string plucked at its center

A string of length $L$ is pulled up a height $h$ at its midpoint and released from rest: $\phi(x) = \frac{2h}{L}x$ for $0\le x\le L/2$, $\phi(x) = \frac{2h}{L}(L - x)$ for $L/2\le x\le L$, and $\psi = 0$, so $B_n = 0$ in (3.18). We need $A_n = \frac2L\int_0^L\phi\sin\frac{n\pi x}{L}dx$.

*Symmetry first.* $\phi$ is symmetric about $x = L/2$. The mode $\sin\frac{n\pi x}{L}$ is symmetric about $L/2$ for odd $n$ and antisymmetric for even $n$: $\sin\frac{n\pi(L - x)}{L} = \sin\big(n\pi - \frac{n\pi x}{L}\big) = (-1)^{n+1}\sin\frac{n\pi x}{L}$. So $A_n = 0$ for even $n$: a symmetric pluck cannot excite an antisymmetric mode. For odd $n$ the two halves contribute equally:
$$A_n = \frac4L\int_0^{L/2}\frac{2h}{L}x\sin\frac{n\pi x}{L}\,dx = \frac{8h}{L^2}\int_0^{L/2}x\sin\frac{n\pi x}{L}\,dx.$$
With $\alpha = n\pi/L$, by parts: $\int_0^{L/2}x\sin\alpha x\,dx = \left[-\frac{x\cos\alpha x}{\alpha}\right]_0^{L/2} + \frac1\alpha\int_0^{L/2}\cos\alpha x\,dx = -\frac{L}{2\alpha}\cos\frac{n\pi}{2} + \frac{1}{\alpha^2}\sin\frac{n\pi}{2}$. For odd $n$, $\cos\frac{n\pi}{2} = 0$ and $\sin\frac{n\pi}2 = (-1)^{(n-1)/2}$, so the integral is $\frac{L^2}{n^2\pi^2}(-1)^{(n-1)/2}$ and
$$A_n = \frac{8h}{n^2\pi^2}(-1)^{(n-1)/2}\ \ (n\text{ odd}): \qquad A_1 = \frac{8h}{\pi^2}\approx0.8106h,\quad A_3 = -\frac{8h}{9\pi^2}\approx-0.0901h,\quad A_5 = \frac{8h}{25\pi^2}\approx 0.0324h.$$
$$u(x,t) = \frac{8h}{\pi^2}\sum_{n\text{ odd}}\frac{(-1)^{(n-1)/2}}{n^2}\sin\frac{n\pi x}{L}\cos\frac{n\pi ct}{L}.$$
Amplitudes fall off as $1/n^2$, Theorem 3.10 for a continuous shape with a corner. Check at the midpoint at $t = 0$: $u(L/2, 0) = \frac{8h}{\pi^2}\sum_{n\text{ odd}}\frac{(-1)^{(n-1)/2}\sin(n\pi/2)}{n^2} = \frac{8h}{\pi^2}\sum_{n\text{ odd}}\frac1{n^2} = \frac{8h}{\pi^2}\cdot\frac{\pi^2}{8} = h$. ✓ The energy in mode $n$ is proportional to $\omega_n^2A_n^2 \propto n^2\cdot n^{-4} = n^{-2}$, so the fundamental carries $\frac{1}{\sum_{\text{odd}}n^{-2}} = \frac{8}{\pi^2}\approx 81\%$ of the energy, the third harmonic $9.0\%$, the fifth $3.2\%$. A string plucked at the center sounds hollow because every even harmonic is missing; pluck it at $L/4$ and they return (project, "Your turn").

### Example 4.8: The Gibbs overshoot, numerically

Theorem 3.15 says the limiting peak of $S_Nf$ for the square wave is $\frac2\pi\operatorname{Si}(\pi)$. Evaluate: $\operatorname{Si}(\pi) = \int_0^\pi\frac{\sin u}{u}du = 1.8519370\ldots$, so the peak is $\frac{2}{\pi}(1.8519370) = 1.1789797$, the overshoot is $0.1789797$, and its fraction of the jump is $\frac{0.1789797}{2} = 0.0894899$. The project computes $S_Nf$ for $N = 1, 3, 5, 11, 51$ on a fine grid and measures $\max S_Nf - 1$: for $N = 11$ it is $0.1813$ and for $N = 51$ it is $0.1791$, both slightly above the limit as §3.9 predicts and converging downward to $0.1790$. It also evaluates the finite-$M$ integral of step 4 by quadrature and finds it matches the measured maxima to four decimals, which checks the derivation and not just the limit.

For the sawtooth, write $x = \pi - t$ in $S_Nf(x) = 2\sum_{n=1}^N\frac{(-1)^{n+1}\sin nx}{n}$. Since $\sin(n\pi - nt) = -(-1)^n\sin nt$, $S_Nf(\pi - t) = 2\sum_{n=1}^N\frac{\sin nt}{n}$. The same telescoping argument (Exercise C1) puts the first maximum at $t_N = \frac{\pi}{N+1}$ with limiting value $2\operatorname{Si}(\pi) = 3.7039$. The function itself is not constant there: $f(\pi - t_N) = \pi - t_N$, so the right quantity is the overshoot *above $f$*, $S_Nf(x) - f(x)$, whose limit is $2\operatorname{Si}(\pi) - \pi = 0.5623 = 0.0895\times2\pi$. ✓ The project measures $0.5643$ at $N = 11$ and $0.5624$ at $N = 51$. Nine percent of whatever the jump is.

## 5. Why it matters and where it connects

### The centerpiece: a Fourier series is an eigenfunction expansion

Everything in §3 rested on Theorem 3.2: the trigonometric functions are orthogonal. We proved it by computing integrals. That proof is correct and tells you nothing about *why*. Here is why.

Let $V$ be the space of twice continuously differentiable functions on $[-\pi,\pi]$ that are $2\pi$-periodic together with their first derivative: $f(-\pi) = f(\pi)$ and $f'(-\pi) = f'(\pi)$. Let $L = \dfrac{d^2}{dx^2}$ act on $V$. Recall from LA-8 that an operator $L$ on an inner product space is **self-adjoint** if $\langle f, Lg\rangle = \langle Lf, g\rangle$ for all $f, g$.

**Theorem 5.1.** $L = d^2/dx^2$ with periodic boundary conditions is self-adjoint with respect to the inner product (3.1).

▶ **Try it first.** Integrate $\int\bar f g''$ by parts twice. Watch what happens to the boundary terms.

**Derivation.**

1. $\langle f, Lg\rangle = \displaystyle\int_{-\pi}^\pi\overline{f}\,g''\,dx = \Big[\overline{f}\,g'\Big]_{-\pi}^{\pi} - \int_{-\pi}^\pi\overline{f'}\,g'\,dx$.
2. The boundary term is $\overline{f(\pi)}g'(\pi) - \overline{f(-\pi)}g'(-\pi) = 0$, because $f$ and $g'$ each take the same value at $\pm\pi$. So $\langle f, Lg\rangle = -\int\overline{f'}g'\,dx = -\langle f', g'\rangle$.
3. Integrate by parts once more, now moving the derivative off $g'$: $-\int\overline{f'}g'\,dx = -\Big[\overline{f'}g\Big]_{-\pi}^\pi + \int\overline{f''}g\,dx = 0 + \langle Lf, g\rangle$, the boundary term vanishing because $f'$ and $g$ are periodic. ∎

■ **What just happened.** Two integrations by parts move $d^2/dx^2$ from one slot to the other, and periodicity kills both boundary terms. The intermediate form says $\langle f, Lf\rangle = -\|f'\|^2 \le 0$: $L$ is **negative semidefinite** in LA-8's sense, with $\langle f, Lf\rangle = 0$ only for constants. Note which boundary conditions make the argument work: periodic ($f(-\pi) = f(\pi)$, $f'(-\pi) = f'(\pi)$), Dirichlet ($f = 0$ at both ends), and Neumann ($f' = 0$ at both ends) all kill the boundary terms. Those are exactly the three boundary conditions of §3.12 and Exercise B7.

**Eigenfunctions.** $Lg = \lambda g$ with periodic boundary conditions means $g'' = \lambda g$. As in step 2 of the heat derivation: $\lambda > 0$ gives real exponentials, which are not periodic; $\lambda = 0$ gives $g = $ const (the linear solution $Ax + B$ is periodic only if $A = 0$); $\lambda = -\mu^2 < 0$ gives $g = Ae^{i\mu x} + Be^{-i\mu x}$, which is $2\pi$-periodic iff $\mu$ is an integer. So the eigenvalues are $\lambda_n = -n^2$, $n = 0, 1, 2, \dots$, with eigenfunctions $e^{inx}$ and $e^{-inx}$ (equivalently $\cos nx$ and $\sin nx$) for $n \ge 1$, and the constant for $n = 0$. Every eigenvalue is real and $\le 0$, as self-adjointness and negative semidefiniteness require (LA-8).

**Theorem 5.2 (Eigenfunctions of a self-adjoint operator with distinct eigenvalues are orthogonal).** If $L$ is self-adjoint, $Lf = \lambda f$, $Lg = \mu g$, and $\lambda\ne\mu$, then $\langle f, g\rangle = 0$.

*Proof.* Eigenvalues of a self-adjoint operator are real (LA-8: $\lambda\|f\|^2 = \langle f, Lf\rangle = \langle Lf, f\rangle = \bar\lambda\|f\|^2$). So $\lambda\langle f, g\rangle = \langle\lambda f, g\rangle = \langle Lf, g\rangle = \langle f, Lg\rangle = \mu\langle f, g\rangle$, where the first equality uses $\bar\lambda = \lambda$ when $\lambda$ comes out of the conjugate-linear slot. Hence $(\lambda - \mu)\langle f, g\rangle = 0$, and $\lambda\ne\mu$. ∎

This is LA-8's spectral theorem argument verbatim, and it proves at once that $\cos mx\perp\cos nx$, $\sin mx\perp\sin nx$, $\cos mx\perp\sin nx$, $1\perp\cos nx$, $1\perp\sin nx$ whenever $m\ne n$: different eigenvalues $-m^2\ne-n^2$. That is five of the six orthogonality relations of Theorem 3.2 with no integral computed. The sixth, $\cos nx\perp\sin nx$ for the *same* $n$, is not forced by $L$ alone: both have eigenvalue $-n^2$, and inside a two-dimensional eigenspace any basis is a basis. We chose an orthogonal one (parity: one is even, the other odd, so their product integrates to zero over a symmetric interval). Equivalently, $\{e^{inx}, e^{-inx}\}$ is the orthogonal basis of that eigenspace singled out by a second self-adjoint operator, $-i\,d/dx$, which commutes with $L$ and has eigenvalues $n$ and $-n$ on them. This is the pattern of simultaneous diagonalization of commuting self-adjoint operators from LA-8, and it is the pattern of QM-5, where $L^2$ and $L_z$ play the roles of $d^2/dx^2$ and $-i\,d/dx$.

So the picture is: **a Fourier series is the expansion of $f$ in the eigenbasis of a self-adjoint operator.** Orthogonality of the basis is the spectral theorem. Completeness of the basis, which is what Parseval needs, is the infinite-dimensional part of the spectral theorem and the one piece this module takes on faith. And the reason the PDEs of §3.12 separate is that $u_t = ku_{xx}$ and $u_{tt} = c^2u_{xx}$ are both of the form "$\dot{\mathbf{u}} = $ (self-adjoint operator)$\,\mathbf{u}$," and in the eigenbasis of that operator they decouple into independent scalar ODEs, one per mode, exactly as $\dot{\mathbf{x}} = A\mathbf{x}$ decouples in the eigenbasis of a symmetric $A$ (LA-7).

### Sturm–Liouville: the same story with a different operator

Replace $d^2/dx^2$ by any operator $Lf = \dfrac{1}{w(x)}\left[\big(p(x)f'\big)' + q(x)f\right]$ on an interval, with inner product $\langle f, g\rangle_w = \int\bar f g\,w\,dx$ and boundary conditions that kill the boundary terms in the same two integrations by parts. The same proof shows $L$ is self-adjoint, so its eigenfunctions are orthogonal with respect to the weight $w$, and (the hard theorem again) complete. Legendre polynomials are the eigenfunctions for $p = 1 - x^2$, $q = 0$, $w = 1$ on $[-1,1]$ (EM-3, from separating Laplace's equation in spherical coordinates). Bessel functions come from $p = x$, $q = -\nu^2/x$, $w = x$ (EM-3 and QM-6, cylindrical and spherical problems). Hermite functions are the eigenfunctions of the harmonic-oscillator Hamiltonian $-\frac{d^2}{dx^2} + x^2$ (QM-3). In every case: separate variables, get a self-adjoint eigenproblem, expand the data in its orthogonal eigenfunctions by (3.5), solve mode by mode. Fourier series is the case you can do with trigonometry; the others need special functions, but the logic is this module's.

### The DFT: sample, and the integrals become sums

Sample $f$ at $N$ equally spaced points $x_j = -\pi + 2\pi j/N$, $j = 0, \dots, N-1$. The trapezoid rule for (3.8) with these nodes (the two endpoint half-weights merge into one weight because $f$ is periodic) is
$$\tilde c_n = \frac1N\sum_{j=0}^{N-1}f(x_j)\,e^{-inx_j}. \tag{5.1}$$
Up to the phase $e^{in\pi} = (-1)^n$ from starting the grid at $-\pi$ instead of $0$, this is the **discrete Fourier transform** of the sample vector $\mathbf{f}\in\mathbb{C}^N$. The sampled exponentials $\mathbf{e}_n = (e^{inx_j})_j$ satisfy $\mathbf{e}_m^\dagger\mathbf{e}_n = \sum_j e^{i(n-m)x_j} = N\delta_{mn}$ for $|m - n| < N$ (a finite geometric sum, the same calculation as Theorem 3.12), so (5.1) is again the projection formula (3.5), now in $\mathbb{C}^N$ with the ordinary dot product. Two things are new. Only $N$ of the $\mathbf{e}_n$ are distinct ($\mathbf{e}_{n+N} = \mathbf{e}_n$: **aliasing**). And the trapezoid rule is *exactly* $\tilde c_n = \sum_{k\in\mathbb{Z}}c_{n+kN}$ (substitute $f = \sum c_m e^{imx}$ into (5.1) and use the discrete orthogonality), so it is accurate precisely when the coefficients beyond $|n|\approx N/2$ are negligible: for a smooth periodic $f$ that is spectacular (the project gets machine precision for $e^{\cos x}$ with 16 points), for the square wave the error is $O(1/N^2)$. F3 builds the DFT and FFT from here, proves the sampling theorem, and shows that the second-difference matrix of §3.12's beaded string is a *circulant*, diagonalized by the DFT with eigenvalues $-\frac{4}{h^2}\sin^2\frac{nh}{2}\approx -n^2$: the discrete Theorem 5.1, which the project checks numerically.

### Where it goes next

- **F3.** Period $T\to\infty$ in (3.9) turns $c_n$ into $\hat f(\xi)$ and the series into the inverse transform. The Dirichlet kernel becomes the sinc function; Theorem 3.13's convolution becomes the convolution theorem; §3.11's delta comb becomes the Dirac comb and the sampling theorem; (5.1) becomes the DFT and FFT.
- **EM-3.** Separation of variables for Laplace's equation in a box gives double Fourier series; in a sphere, Legendre polynomials replace cosines. Boundary data is expanded by (3.5).
- **QM-3.** The infinite square well's energy eigenstates are $\sqrt{2/L}\sin(n\pi x/L)$ with $E_n\propto n^2$: (3.17) with $e^{-iE_nt/\hbar}$ in place of $e^{-k\lambda_nt}$. Expanding $\psi(x,0)$ in them is a sine series; $|b_n|^2$ is the probability of measuring $E_n$, and Parseval is the statement that probabilities sum to one.
- **ST-8.** The periodogram of a time series is $|\tilde c_n|^2$ from (5.1); Parseval is the decomposition of variance by frequency; the spectral density is its expectation.
- **CM-6.** Normal modes of coupled oscillators are eigenvectors of a symmetric matrix; §3.12's beaded string is the example whose limit is a Fourier sine series.
- **ML-8.** A convolutional layer computes Theorem 3.13's convolution with a learned kernel; the DFT diagonalizes it, which is why convolution is cheap in the frequency domain.
- **QC-5.** The quantum Fourier transform is (5.1) as a unitary matrix on $N = 2^n$ amplitudes.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are drills: orthogonality integrals, reading coefficients off a trigonometric polynomial, converting between real and complex forms, a half-range sine series, and translating Osgood's period-1 convention. B1–B8 are the module: B1 (sine series of $x(\pi - x)$), B2 (Bessel's inequality from the projection identity, in complex form), B3 (two more zeta values via Parseval), B4 (decay rate of $(\pi^2 - x^2)^2$ via the derivative trick), B5 (Dirichlet kernel properties), B6 (termwise integration from the sawtooth to $x^2$), B7 (heat equation with insulated ends), and B8 (extend `project.py`). B2, B4, B6, and B7 are essential. C1 (Gibbs for the sawtooth), C2 (the isoperimetric inequality via Parseval), C3 (Fejér's kernel and the disappearance of Gibbs) go one step further.

## 7. Coding project

`project.py` builds and checks everything in this module:

1. **`fourier_coeffs(f, N, T)`** computes $c_{-N},\dots,c_N$ by the composite trapezoid rule, and `real_coeffs` converts to $a_n, b_n$ via (3.7). For the four standard waves the result is compared with SymPy's exact integrals and must agree to $10^{-8}$. The periodic-trapezoid accuracy is demonstrated: for $e^{\cos x}$ the error reaches machine precision with 16 sample points; for the square wave it decays like $1/M^2$.
2. **Partial sums** $S_Nf$ for $N = 1, 3, 5, 11, 51$ are plotted for the square wave and sawtooth. The **Gibbs overshoot** is measured and compared with $0.0895\times$(jump), $\operatorname{Si}(\pi)$ is checked against (3.15), and the finite-$M$ integral of the derivation is checked against the measured peaks.
3. **Parseval** is verified for all four waves; the Basel partial sums converge to $\pi^2/6$ and the $x^2$ sums give $\pi^4/90$.
4. **Decay rates:** $\log|c_n|$ versus $\log n$ for the square wave, triangle, and cubic $x(\pi^2 - x^2)$ from *numerical* coefficients, with fitted slopes $-1, -2, -3$.
5. **Dirichlet kernel** plotted for several $N$; its closed form (3.13) checked against the sum, its mean checked to be $1$; the convolution formula (3.14) checked against the direct partial sum.
6. **Heat equation** on a rod with the step of Example 4.6: truncated series versus an explicit finite-difference solver, overlaid at several times; the copper time scales $\tau_n$ printed.
7. **Plucked string:** mode amplitudes $A_n$ checked against Example 4.7 and the energy fractions; a figure grid of the first six normal modes and of snapshots of the plucked string.
8. **Self-adjointness numerically:** the periodic second-difference matrix is symmetric; in the DFT basis it is diagonal and real with entries $\approx -n^2$ (Theorem 5.1 and the F3 circulant preview); the spectral second derivative is exactly $-n^2$; sampled cosines are exactly orthogonal in $\mathbb{R}^N$ (§2).

**Build it yourself first:** write `fourier_coeffs` (sample, weight, sum) and `partial_sum(n, c, x)` before reading them. Then write the finite-difference heat solver (one line of NumPy per time step) and compare it with the series. The "Your turn" section asks for the full-wave rectified sine and its DC content, the sawtooth on $[0, 2\pi)$ and the phase shift of its coefficients, the number of terms for uniform error $< 0.01$ on $x^2$ and $|x|$, a pluck at $L/4$, and Fejér averaging.

## 8. Where you'll struggle

- **Factors of $2$ and $\pi$.** Everyone gets them wrong, repeatedly, and every book has a different convention. The fix is not to memorize. Write the projection formula $\frac{\langle\varphi_n, f\rangle}{\langle\varphi_n,\varphi_n\rangle}$ (3.5) and compute $\langle\varphi_n,\varphi_n\rangle$ for the basis in front of you: $2\pi$ for $1$ and for $e^{inx}$ on $[-\pi,\pi]$; $\pi$ for $\cos nx$ and $\sin nx$; $L/2$ for $\sin(n\pi x/L)$ on $[0, L]$; $T$ for $e^{2\pi int/T}$ over a period. The constants fall out. If you ever write $a_0$ as the constant term instead of $a_0/2$, you have mixed conventions; go back to (3.5).
- **Splitting the integral at the breaks.** For a piecewise function, $\int_{-\pi}^{\pi}f\cos nx\,dx$ is a *sum* of integrals, one per piece, each with its own formula for $f$. Writing a single integral with "$f$" in it and integrating by parts as if $f$ were one formula is the most common computational error in this module, and it silently drops the boundary terms that Theorem 3.9 shows are the entire $1/n$ content of a jump. Examples 4.1 and 4.6 show the split.
- **Dirichlet's theorem is genuinely hard, and we only sketched it.** What the sketch in §3.8 glosses: the boundedness of $g_x(t)$ near $t = 0$ (fine for piecewise smooth $f$; fatal for merely continuous $f$, where the theorem is false), and Riemann–Lebesgue for functions that are integrable but not square-integrable. Mean-square convergence and Parseval for general $f$ (completeness of the basis) were taken on faith entirely; Stein & Shakarchi ch. 3 fills that in. If you want to know what you are standing on, read it. If you want to use Fourier series, piecewise smooth is all you will ever need.
- **Gibbs is not a bug.** Every student who first plots $S_{51}f$ for a square wave assumes the $9\%$ ears near the jump are a numerical artifact and tries to fix the code. They are a theorem (3.15). Adding terms moves the ears toward the jump and does not shrink them. If you need uniform convergence near a jump, change the summation method (Fejér, Exercise C3), not the number of terms.
- **Mean-square versus pointwise versus uniform.** Three different statements. $\|f - S_Nf\|\to0$ (always, for $f\in PC$). $S_Nf(x)\to\frac{f(x^+)+f(x^-)}{2}$ at each $x$ (piecewise smooth). $\max|f - S_Nf|\to0$ (needs $\sum|c_n|<\infty$, e.g. continuous and piecewise smooth). The square wave satisfies the first two and not the third. Keep them separate.
- **Differentiating termwise when you may not.** The sawtooth's termwise derivative $2\sum(-1)^{n+1}\cos nx$ diverges, and this is correct behavior: the periodic sawtooth has jumps, so its derivative contains delta functions, which have non-decaying Fourier coefficients. Theorem 3.18 requires $f$ continuous *as a periodic function*. Check $f(-\pi) = f(\pi)$ before differentiating any series. Integration (Theorem 3.17) never has this problem.
- **The periodic extension is what the series sees.** $x^2$ is smooth; its periodic extension has a corner at $\pm\pi$, and the $1/n^2$ decay reflects the corner, not the smoothness. When predicting decay rates, draw two periods of the extension and count derivatives across the seam.
- **Example 4.7's symmetry argument.** Students compute $A_n$ for even $n$ by brute force, get $0$, and don't know why. It is parity about the midpoint. Learn to look for it first: it halves the work and explains the result.

## 9. Check yourself

1. State the projection formula (3.5) and use it to derive the $\frac{a_0}{2}$ in the real Fourier series.
2. Compute $\int_{-\pi}^\pi\cos 3x\cos 5x\,dx$ and $\int_{-\pi}^\pi\cos^2 3x\,dx$ from the product-to-sum formulas, without looking anything up.
3. Why does $\langle e^{imx}, e^{inx}\rangle$ have $e^{-imx}$ in the integrand? Which slot of the inner product is conjugate-linear in this curriculum?
4. Write $c_n$ in terms of $a_n, b_n$, and $a_n, b_n$ in terms of $c_{\pm n}$. For real $f$, what is $c_{-n}$?
5. State Theorem 3.5 in one sentence using the word "projection," and write the identity $\|f\|^2 = \|S_Nf\|^2 + \|f - S_Nf\|^2$ from it.
6. Write Bessel's inequality and Parseval's identity in complex form. Which one is proved in this module without any further theorem, and what does the other one require?
7. If $|c_n|\sim 1/n^3$, how smooth is the periodic extension of $f$? If $f$ has a jump, how fast can $|c_n|$ decay? Use Theorem 3.9 to write the square wave's $c_n$ without an integral.
8. Derive $D_N(x) = \frac{\sin((N+\frac12)x)}{\sin(x/2)}$ in three lines. What is $\frac{1}{2\pi}\int_{-\pi}^\pi D_N$?
9. State Dirichlet's theorem with its hypothesis. What does the series converge to at a jump? Which step of the proof uses piecewise smoothness?
10. What is the Gibbs overshoot as a fraction of the jump, and why doesn't it go away as $N\to\infty$?
11. Why may you always integrate a Fourier series termwise but not always differentiate it? Give the example that fails.
12. Write the solution of $u_t = ku_{xx}$ on $[0,L]$ with zero ends and initial data $\phi$. What is the time constant of mode $n$?
13. Prove in three lines that $d^2/dx^2$ with periodic boundary conditions is self-adjoint. Which five of the six orthogonality relations in Theorem 3.2 follow from this and Theorem 5.2? Which one does not, and why?

## 10. Pair with

- **Osgood, *Lectures on the Fourier Transform and Its Applications*, ch. 1** (the Stanford EE261 notes). Osgood uses period $1$ and $e^{2\pi int}$ throughout, the convention F3 adopts; his chapter 1 is the best informal account of why the complex form is the right one, and his discussion of what the Fourier series of a square wave *means* is worth reading after §3.9.
- **Stein & Shakarchi, *Fourier Analysis: An Introduction*, ch. 1–3.** Ch. 1 is the history and the heat/wave motivation of §1 with full derivations; ch. 2 has the convergence theory of §3.8–3.9 with complete proofs, and Fejér's theorem (Exercise C3) in §5; ch. 3 proves mean-square convergence and Parseval (the completeness we took on faith) and applies them to the isoperimetric inequality (Exercise C2). The rigor missing here is all there.
- **MIT 18.06SC (Strang), the Fourier series unit**, and Strang, *Introduction to Linear Algebra*, §10.5 "Fourier Series: Linear Algebra for Functions" in the editions that have it. Strang does exactly what §2 does, functions as infinite vectors and coefficients as projections, in twenty minutes. Watch it before §3 if §2 felt abstract.
- **Körner, *Fourier Analysis*, Part I.** Short chapters with historical asides, a complete proof of Dirichlet's theorem, and du Bois-Reymond's continuous function with a divergent Fourier series. Read the counterexample once to understand why the hypothesis of Theorem 3.14 matters.
