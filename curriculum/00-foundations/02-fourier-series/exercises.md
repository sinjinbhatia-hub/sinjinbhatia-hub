# F2 Exercises

Conventions as in the notes: $\langle f, g\rangle = \int_{-\pi}^{\pi}\overline{f}\,g\,dx$ (conjugate on the first slot), $a_n = \frac1\pi\int f\cos nx$, $b_n = \frac1\pi\int f\sin nx$, $c_n = \frac{1}{2\pi}\int f e^{-inx}$, constant term $a_0/2 = c_0$. Equation numbers refer to `notes.md`. Every number you are asked to compute is also computed in `project.py` (Part 9, "Exercise checks"), so you can check yourself.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Orthogonality by hand.) Using only the product-to-sum formulas of §3.2 step 2 and the fact that $\int_{-\pi}^{\pi}\cos kx\,dx = 0$ for every nonzero integer $k$, compute
(a) $\int_{-\pi}^{\pi}\cos 2x\cos 5x\,dx$, (b) $\int_{-\pi}^{\pi}\sin^2 3x\,dx$, (c) $\int_{-\pi}^{\pi}\cos 4x\sin 4x\,dx$, (d) $\int_{-\pi}^{\pi}e^{-2ix}e^{3ix}\,dx$, (e) $\langle e^{3ix}, e^{3ix}\rangle$.
Then write down $\|1\|^2$, $\|\cos 7x\|^2$, $\|\sin 7x\|^2$, $\|e^{7ix}\|^2$ without integrating.

**A2.** (Reading coefficients off a trigonometric polynomial.) Let $f(x) = 2 + 3\cos x - \sin x + 4\sin 2x - 5\cos 3x$.
(a) Write down $a_0$, every $a_n$ and $b_n$, and $c_n$ for $|n|\le3$, without integrating. Say in one sentence why no integral is needed (Theorem 3.5).
(b) Compute $\|f\|^2$ from the coefficients (step 1 of the Bessel derivation), in both real and complex form.
(c) What is the closest trigonometric polynomial of degree $\le1$ to $f$ in mean square, and what is $\|f - S_1f\|^2$? Same question for degree $\le2$.

**A3.** (Real $\leftrightarrow$ complex.)
(a) A real function has $c_0 = 3$, $c_{\pm1} = 1\mp2i$, $c_{\pm2} = \pm i/2$, all other $c_n = 0$. Write its real Fourier series.
(b) Without integrating, find the Fourier series of $\cos^2x$ and of $\sin^3x$ in both real and complex form (trig identities, or write $\cos x = (e^{ix}+e^{-ix})/2$ and expand).
(c) Verify your $a_2$ for $\cos^2x$ by computing (3.4) directly.

**A4.** (A half-range sine series.) Expand $f(x) = 1$ on $(0, L)$ in a sine series (§3.5). Evaluate the series at $x = L/2$ and at $x = 0$. Which function on $(-L, L)$ does the series represent, and how is it related to Example 4.1? Compute $S_5f(L/2)$ numerically.

**A5.** (Translating a convention.) Osgood writes a $1$-periodic function as $f(t)\sim\sum_n\hat c_n e^{2\pi int}$ with $\hat c_n = \int_{-1/2}^{1/2}f(t)e^{-2\pi int}\,dt$. Compute $\hat c_n$ for $f(t) = t$ on $[-\tfrac12,\tfrac12)$ directly by parts. Then get the same answer from Example 4.2 by the substitution $x = 2\pi t$ in (3.9). Convert to real form and check the result at $t = 1/4$.

## B. Core (do all; these are the module)

**B1.** (Sine series of $x(\pi - x)$.) Expand $f(x) = x(\pi - x)$ on $[0,\pi]$ in a sine series.
(a) Compute $b_n$ and show $b_n = \dfrac{8}{\pi n^3}$ for odd $n$ and $0$ for even $n$. Give $b_1, b_3, b_5$ to four decimals.
(b) Evaluate at $x = \pi/2$ to get $\sum_{n\text{ odd}}\dfrac{(-1)^{(n-1)/2}}{n^3} = \dfrac{\pi^3}{32}$.
(c) Sketch two periods of the odd $2\pi$-periodic extension and explain the $1/n^3$ decay with Theorem 3.10: how many derivatives are continuous across the seam, which derivative jumps, and where? Then recover $b_n$ from Theorem 3.9 applied to that jumping derivative, with no integral.
(d) Why do the even coefficients vanish? (Parity about $x = \pi/2$.)

**B2.** (Bessel in complex form.) Let $f\in PC$ with coefficients $c_n$ as in (3.8) and $S_Nf = \sum_{|n|\le N}c_ne^{inx}$.
(a) Show $\langle e^{imx}, f - S_Nf\rangle = 0$ for $|m|\le N$, hence $\langle S_Nf, f - S_Nf\rangle = 0$.
(b) Show $\|S_Nf\|^2 = 2\pi\sum_{|n|\le N}|c_n|^2$.
(c) Conclude $\|f - S_Nf\|^2 = \|f\|^2 - 2\pi\sum_{|n|\le N}|c_n|^2$ and Bessel's inequality in the form $\sum_{|n|\le N}|c_n|^2\le\frac{1}{2\pi}\|f\|^2$.
(d) For the square wave of Example 4.1 compute the relative mean-square error $\frac{1}{2\pi}\|f - S_Nf\|^2$ for $N = 1, 5, 11, 51$, and compare with the overshoot $\max_x S_Nf(x) - 1$ (the Gibbs ears, about $0.18$ for every large $N$) and with the worst-case error $\sup_x|f(x) - S_Nf(x)|$ (what is it, and why?). Which of these does Theorem 3.5 control?

**B3.** (Two more zeta values.)
(a) Apply Parseval (3.11) to the triangle wave of Example 4.3 and deduce $\sum_{n\text{ odd}}1/n^4 = \pi^4/96$. Combine with Example 4.5's $\sum1/n^4 = \pi^4/90$ to find $\sum_{n\text{ even}}1/n^4$, and check it against $\frac{1}{16}\sum1/n^4$.
(b) Show that Parseval for a half-range sine series on $[0,\pi]$ reads $\frac2\pi\int_0^\pi f^2\,dx = \sum b_n^2$, and apply it to B1 to get $\sum_{n\text{ odd}}1/n^6 = \pi^6/960$. Deduce $\sum_{n\ge1}1/n^6$.
(c) Explain why Parseval applied to piecewise-polynomial functions only ever produces *even* zeta values, and what B1(b) produced instead.

**B4.** (Decay rate by the derivative trick.) Let $f(x) = (\pi^2 - x^2)^2$ on $[-\pi,\pi]$, extended periodically.
(a) Show $f, f', f''$ are continuous as periodic functions and $f'''$ jumps at the seam $x = \pm\pi$; find the jump.
(b) Without computing $\int f\cos nx$, use $c_n(g') = in\,c_n(g)$ (Theorem 3.18) three times together with the sawtooth's $c_n(x) = i(-1)^n/n$ to show $c_n(f) = -\dfrac{24(-1)^n}{n^4}$ for $n\ne0$, so $a_n = -\dfrac{48(-1)^n}{n^4}$. Check it against Theorem 3.9 applied to $g = f'''$.
(c) Compute $a_0/2 = 8\pi^4/15$ directly.
(d) Evaluate the series at $x = \pi$ and at $x = 0$ to recover $\sum1/n^4 = \pi^4/90$ and to find $\sum(-1)^{n+1}/n^4$.
(e) Which $k$ in Theorem 3.10 applies, what decay does it predict, and what constant $C$ does it give? Compare $C$ with the actual $|c_n|n^4$.

**B5.** (The Dirichlet kernel.)
(a) Show $D_N(x) = 1 + 2\sum_{n=1}^N\cos nx$, that $D_N$ is even, and that $\frac{1}{2\pi}\int_{-\pi}^\pi D_N = 1$.
(b) Where are the zeros of $D_N$ in $(0, 2\pi)$? How wide is the central lobe?
(c) Show $\int_{-\pi}^{\pi}D_N(x)^2\,dx = 2\pi(2N+1)$. Evaluate for $N = 10$.
(d) Prove that the Lebesgue constant $L_N = \frac{1}{2\pi}\int_{-\pi}^{\pi}|D_N(x)|\,dx$ satisfies $L_N\ge\frac{4}{\pi^2}\ln(N+1)$, so $L_N\to\infty$. Hint: $0 < \sin(x/2)\le x/2$ for $0 < x\le\pi$; substitute $u = (N+\tfrac12)x$; on $[(k-1)\pi, k\pi]$ use $1/u\ge1/(k\pi)$.
(e) Show $|S_Nf(x)|\le L_N\max|f|$ and explain in two or three sentences why (d) is the reason $S_Nf$ can fail to converge for a merely continuous $f$, while a *positive* kernel of unit mass (C3) can never misbehave this way.

**B6.** (Termwise integration: sawtooth $\to x^2\to$ cubic.)
(a) Integrate the sawtooth series $x\sim2\sum(-1)^{n+1}\frac{\sin nx}{n}$ from $0$ to $x$ term by term (Theorem 3.17 permits it). Using $\sum(-1)^{n+1}/n^2 = \pi^2/12$ from Example 4.4, recover the series of $x^2$ in Example 4.4.
(b) Theorem 3.17 requires removing the mean before the integrated series is periodic. Subtract $\pi^2/3$ from $x^2$, integrate from $0$ to $x$, and show $x(\pi^2 - x^2) = 12\sum(-1)^{n+1}\frac{\sin nx}{n^3}$, the cubic of the §3.7 table.
(c) Evaluate (b) at $x = \pi/2$ and reconcile with B1(b).
(d) Differentiate the $x^2$ series termwise. What do you get, is it legitimate, and what does it converge to at $x = \pi$? Why is the same operation on the sawtooth series not legitimate?

**B7.** (Heat with insulated ends.) Solve $u_t = k\,u_{xx}$ on $[0,L]$ with $u_x(0,t) = u_x(L,t) = 0$ and $u(x,0) = \phi(x)$.
(a) Separate variables as in §3.12 and show the eigenfunctions are $\cos(n\pi x/L)$ for $n = 0, 1, 2,\dots$, with $\lambda_0 = 0$ included. Treat the three sign cases of $\lambda$.
(b) Write the solution as a cosine series with the $n = 0$ term explicit. What is $\lim_{t\to\infty}u(x,t)$? Show directly from the PDE that $\frac{d}{dt}\int_0^Lu\,dx = 0$ and explain what that has to do with the limit.
(c) For the step $\phi = u_0$ on $(0, L/2)$, $0$ on $(L/2, L)$ (the data of Example 4.6), compute $a_0/2$ and $a_1,\dots,a_6$ to four decimals and write the first three nonzero terms of $u$. Why do the even modes vanish here when they did not in Example 4.6?
(d) Why does the Neumann problem have the eigenvalue $0$ while the Dirichlet problem does not? Answer with the intermediate form $\langle f, Lf\rangle = -\|f'\|^2$ of Theorem 5.1.

**B8.** (Code.) Extend `project.py`.
(a) Add `quartic(x) = (π² − x²)²` to the decay-rate study of Part 4: compute its $c_n$ numerically with `fourier_coeffs`, fit the log–log slope (expect $-4$), and assert agreement with B4's closed form $c_n = -24(-1)^n/n^4$ and $c_0 = 8\pi^4/15$ to $10^{-8}$.
(b) Write `terms_for_uniform_error(name, tol)` returning the smallest $N$ with $\max_x|f - S_Nf| < \mathrm{tol}$ on a fine grid, and run it for $x^2$ and $|x|$ with tol $= 0.01$. Compare with the tail bounds $4/N$ and $2/(\pi N)$ from Theorem 3.16. Recomputing $S_Nf$ from scratch for each $N$ costs $O(N^2\times\text{grid})$; find the $O(N\times\text{grid})$ way.

## C. Challenge (do at least one)

**C1.** (Gibbs for the sawtooth.) Let $f(x) = x$ on $(-\pi,\pi)$ and $S_Nf(x) = 2\sum_{n=1}^N(-1)^{n+1}\frac{\sin nx}{n}$.
(a) Show $S_Nf(\pi - t) = 2\sum_{n=1}^N\frac{\sin nt}{n} =: T_N(t)$.
(b) Show $T_N'(t) = D_N(t) - 1 = \dfrac{2\cos\frac{(N+1)t}{2}\sin\frac{Nt}{2}}{\sin\frac t2}$ and that the first maximum of $T_N$ on $(0,\pi)$ is at $t_N = \dfrac{\pi}{N+1}$.
(c) Show $T_N(t_N)$ is a Riemann sum for $2\int_0^\pi\frac{\sin u}{u}du$ and conclude $\lim_N T_N(t_N) = 2\operatorname{Si}(\pi)\approx3.7039$.
(d) The function itself is $f(\pi - t_N) = \pi - t_N$ there, so the overshoot *above $f$* tends to $2\operatorname{Si}(\pi) - \pi\approx0.5623 = 0.0895\times2\pi$. Compute $T_N(t_N) - (\pi - t_N)$ for $N = 11$ and $N = 51$, compare with the project's grid measurement of $\max_x(S_Nf - f)$ ($0.5643$ and $0.5624$), and explain why the two differ at finite $N$ and by how much.
(e) Generalize: let $f$ be $2\pi$-periodic and piecewise smooth with a jump $J = f(x_0^+) - f(x_0^-) > 0$ at $x_0$. Show that the peak of $S_Nf$ just to the right of $x_0$ overshoots $f(x_0^+)$ by $(\frac{\operatorname{Si}(\pi)}{\pi} - \frac12)J\approx0.0895\,J$ in the limit: precisely, $\lim_{\delta\to0^+}\limsup_{N\to\infty}\max_{x_0 < x < x_0+\delta}\big(S_Nf(x) - f(x_0^+)\big) = 0.0895\,J$, and the same with $\liminf$. (By the same argument $S_Nf$ undershoots $f(x_0^-)$ by $0.0895\,J$ just to the left; for a downward jump apply the result to $-f$.) Hint: let $\sigma$ be the $2\pi$-periodic function equal to $\frac{\pi - x}{2\pi}$ on $(0, 2\pi)$, a sawtooth with a single jump of $+1$ per period, at $0$, and $S_N\sigma = T_N/2\pi$ by (a). Subtract $J_j\,\sigma(x - x_j)$ at every jump $x_j$ of $f$ and apply Theorem 3.16 to what is left. Then show, using $T_N' = D_N - 1$ and one integration by parts, that $T_N(t)\to\pi - t$ uniformly on $[\delta, 2\pi - \delta]$, and that $T_N(t) - (\pi - t)$ is largest on $(0, \pi]$ at the first zero of $D_N$.

**C2.** (The isoperimetric inequality, Hurwitz 1902.) Let $\gamma(s) = (x(s), y(s))$ be a simple closed piecewise-smooth curve of length $2\pi$, parametrized by arc length, so $x'^2 + y'^2 = 1$ and $x, y$ are continuous and $2\pi$-periodic. Let $\alpha_n, \beta_n$ be the complex Fourier coefficients of $x$ and $y$.
(a) Derive the polarized Parseval identity $\langle f, g\rangle = 2\pi\sum_n\overline{c_n(f)}\,c_n(g)$ from (3.11) applied to $f + g$ and $f + ig$.
(b) Show $\sum_n n^2\big(|\alpha_n|^2 + |\beta_n|^2\big) = 1$.
(c) Green's theorem gives the enclosed area $A = \int_0^{2\pi}x\,y'\,ds$. Show $A = 2\pi\sum_n in\,\overline{\alpha_n}\beta_n$.
(d) Prove $|A|\le\pi$ using $|n|\le n^2$ and $2|\alpha_n||\beta_n|\le|\alpha_n|^2 + |\beta_n|^2$. Show equality forces $\alpha_n = \beta_n = 0$ for $|n|\ge2$ and $|\alpha_{\pm1}| = |\beta_{\pm1}|$, and that the curve is then a circle of radius $1$.
(e) Rescale to length $L$: $A\le L^2/(4\pi)$. Check numerically on the ellipse $x = 2\cos\theta$, $y = \sin\theta$ (the project computes $L$, $L^2/4\pi$ and $A$), and verify (c) on that ellipse with $\theta$ as the parameter. Why is (c) valid for any parametrization while (b) needs arc length?

**C3.** (Fejér's kernel and the end of Gibbs.) Define $\sigma_Nf = \dfrac{S_0f + S_1f + \cdots + S_Nf}{N+1}$.
(a) Show $\sigma_Nf = \sum_{|n|\le N}\Big(1 - \dfrac{|n|}{N+1}\Big)c_ne^{inx}$.
(b) Show $\sigma_Nf(x) = \frac{1}{2\pi}\int_{-\pi}^\pi f(x-t)F_N(t)\,dt$ with $F_N = \frac{D_0 + \cdots + D_N}{N+1}$, and derive $F_N(t) = \dfrac{1}{N+1}\cdot\dfrac{\sin^2\frac{(N+1)t}{2}}{\sin^2\frac t2}$. Hint: multiply $\sum_{k=0}^N\sin\big((k+\tfrac12)t\big)$ by $2\sin\frac t2$ and telescope.
(c) Show $F_N\ge0$ and $\frac{1}{2\pi}\int_{-\pi}^\pi F_N = 1$.
(d) Conclude $\min f\le\sigma_Nf(x)\le\max f$ for every $x$ and $N$: no overshoot, ever. For the square wave, what are $\max_x\sigma_{11}f$ and $\max_x\sigma_{51}f$? (The project computes them.)
(e) Prove Fejér's theorem: if $f$ is continuous and $2\pi$-periodic then $\sigma_Nf\to f$ uniformly. Hint: for $\delta\le|t|\le\pi$, $F_N(t)\le\frac{1}{(N+1)\sin^2(\delta/2)}$. Then explain how this supplies the completeness that Parseval needed in §3.6.
(f) Add $\sigma_N$ to `project.py` and plot it next to $S_N$ for the square wave at $N = 51$.
