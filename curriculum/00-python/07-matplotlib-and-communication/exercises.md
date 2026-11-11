# PY-7 Exercises

Conventions as in the notes: `import matplotlib.pyplot as plt`, the OO API, `ph = plot_helpers` (the file in this folder), `MPLBACKEND=Agg`, `rng = np.random.default_rng(seed)`. **Write your prediction down before you run anything**, test with `assert`s on *data* (limits, tick text, array values, artist properties) and look at the saved file once. "Write a function" means: write it, then test it, including an edge case. Put answers in a file next to `project.py` (importing `project` re-runs its checks, about half a minute, so copy what you need).

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Object model.) For `fig, axs = plt.subplots(1, 3)` and `fig2, ax2 = plt.subplots(2, 2, sharey="row")`, predict `type(axs)`, `axs.shape`, `ax2.shape`, `len(fig.axes)`, and whether `ax2[0, 0]`'s $y$-axis is shared with `ax2[0, 1]` and with `ax2[1, 0]`. Then check each with code (`get_shared_y_axes().joined(a, b)`).

**A2.** (The state machine.) Predict where each call draws, then rewrite the function with `fig, ax = plt.subplots()` so it returns the two figures, and assert that the label ended up on the first figure's axes.

```python
def plot_pair(x, y1, y2):
    plt.figure()
    plt.plot(x, y1)
    plt.figure()
    plt.plot(x, y2)
    plt.xlabel("x (m)")          # intended for the FIRST figure
    plt.title("both")
```

**A3.** (Which scale?) For each model, say which axes make it a straight line (linear, `semilogy`, `semilogx`, `loglog`) and the slope you would read off: (a) $y = 3x^2$; (b) $y = 2e^{-x/5}$; (c) $y = 5/\sqrt x$; (d) $y = 2 + 3\log_{10}x$. Confirm each with `np.polyfit` in the right coordinates and an assertion that $R^2 > 0.9999$.

**A4.** (`imshow`.) `A = np.zeros((4, 6)); A[0, 3] = 1`. Where is the bright cell, in data coordinates and on the page (top or bottom), for (i) `imshow(A)` and (ii) `imshow(A, origin="lower", extent=(0, 6, 0, 4))`? Check by reading the rendered pixel colour with `ph.rgba_at` at your predicted cell centre and at the mirror-image cell.

**A5.** (Tick formatter.) Write a `FuncFormatter` function `deg(v, pos)` that labels a radian axis in degrees (`"90°"`), and put ticks every $\pi/6$ on $[0, 2\pi]$ with `MultipleLocator`. Test it on `np.pi / 6 * np.arange(13)` and on the value `-1e-16`, which a locator can produce for zero. What goes wrong with the naive `f"{np.degrees(v):.0f}°"`, and what is the fix?

**A6.** (Norms and maps.) Predict `Normalize(2, 6)(5)`, `LogNorm(1, 100)(10)`, `Normalize(-3, 3)(0)` and `TwoSlopeNorm(vmin=-1, vcenter=0, vmax=5)(2.5)`. Then, for each dataset, name a colour map family and a norm: (i) a temperature anomaly in K, symmetric about 0; (ii) wind direction in degrees; (iii) population density from $1$ to $10^5$ per km²; (iv) cluster labels 0 to 5.

## B. Core (do all; these are the module)

**B1.** (Three integrators, one picture.) Solve $y' = -y$, $y(0) = 1$ on $[0, 2]$ with forward Euler, Heun's method (RK2) and classical RK4 for $n = 8, 12, 16, 24, 32, 48, 64, 96, 128, 192$ steps; the error is $|y_n - e^{-2}|$. Plot all three on one log-log axes with `ph.plot_convergence` (orders 1, 2, 4), label each series directly with `ph.label_end` instead of a legend, and assert each fitted slope is within 0.1 of its order. Which of the three slopes is furthest from its order, and why?

**B2.** (Two slopes in one error curve.) The forward difference $D_h = (\sin(1+h) - \sin 1)/h$ approximates $\cos 1$. Plot $|D_h - \cos 1|$ for $h = 10^{-1}, \dots, 10^{-16}$ on log-log axes. Truncation error is $\approx \tfrac h2|f''|$ (slope $+1$), roundoff is $\approx \varepsilon|f|/h$ (slope $-1$), with $\varepsilon \approx 2.2\times10^{-16}$. Split the data at the minimum, fit each side separately with `ph.loglog_fit`, draw a slope guide on each branch, and assert: the right branch ($h \ge 10^{-6}$) has slope $1\pm0.05$, the left branch ($h \le 10^{-10}$) slope $-1\pm0.3$, and the minimiser is within a factor of 10 of $\sqrt\varepsilon$. Why is the left branch noisy?

**B3.** (A limit cycle, CM-1 preview.) The Van der Pol oscillator $\ddot x - \mu(1-x^2)\dot x + x = 0$, $\mu = 1$, as the system $\dot x = y$, $\dot y = \mu(1-x^2)y - x$. Draw `quiver` (unit-length arrows with `angles="xy"`, as in Example 4.7), `streamplot`, and two `solve_ivp` trajectories, one starting inside the cycle at $(0.1, 0)$ and one outside at $(3, 3)$, up to $t = 60$. Assert that both trajectories settle on the same cycle: after $t = 40$, $\max|x|$ is within 0.05 of 2.0 for both. Label the axes with their quantities.

**B4.** (Intervals and bandwidths.) (a) Draw $K = 4000$ samples of size $n = 5$ from $N(0, 3^2)$. For each, form $\bar x \pm 1.96\,\mathrm{SE}$ and the $t$-interval $\bar x \pm t_{0.975,4}\,\mathrm{SE}$ (3.9). Estimate each interval's coverage of the true mean 0 and assert that the $t$-interval is within 0.02 of 0.95 and the $1.96$ interval is below 0.92. (b) For the Gamma sample of Example 4.3, vary the bandwidth factor of `gaussian_kde` over 25 log-spaced values from 0.02 to 1.5; plot the total-variation distance to the true density against the factor on log-$x$ axes with Scott's factor marked, and assert that the best factor is within a factor of 3 of Scott's and that both extremes are worse than the best.

**B5.** (Colour maps by lightness.) Using $L^*$ (3.7), write `classify(name)` returning `"sequential"` (monotone in one direction and a range of more than 50), `"diverging"` (one ridge within 25 steps of the centre and end lightnesses within 5 of each other), or `"other"`. Run it on `viridis, cividis, magma, plasma, inferno, gray, coolwarm, RdBu_r, jet, turbo` and assert your table: the first six sequential, the next two diverging, the last two other. Both `"other"` maps are rainbows: where along each map is the lightest colour, and what does that do to a single smooth bump?

**B6.** (Slide and paper.) Wrap Example 4.3's histogram figure in `make_fig()`. Produce (a) a slide version: `figsize=(10, 5.6)`, `font.size` 16, thicker lines, at most five tick intervals per axis (`MaxNLocator(5)`, so at most six ticks), PNG; (b) a paper version: `figsize=(3.4, 2.6)`, `font.size` 8, `pdf.fonttype` 42, `fmt="pdf"`. Use `plt.rc_context` so nothing leaks, and assert: the figure sizes, the tick counts in view, that the PDF starts with `%PDF`, and that `plt.rcParams["font.size"]` is unchanged afterwards. In two sentences: what does the paper version cut that the slide keeps?

**B7.** (Extend `project.py`: small multiples with a shared norm.) Add Part 7b: six snapshots of 2D diffusion, $u(x,y,t) = \exp\!\big(-r^2/2s^2\big)/(2\pi s^2)$ with $s^2 = \sigma_0^2 + 2Dt$, $\sigma_0 = 0.3$, $D = 0.1$, at $t = 0, 0.5, 1, 2, 4, 8$ on $[-6, 6]^2$ ($161\times161$). Draw a $2\times3$ grid of `imshow` panels with `origin="lower"` and a correct `extent`, shared axes, **one** `Normalize(0, u_max(t=0))` and **one** colorbar. Assert: all six images have the same colour limits, `len(fig.axes) == 7`, and the total mass $\sum u\,\Delta x\,\Delta y$ is 1 within $10^{-3}$ in every panel. Then redo it with per-panel autoscaling and print the six colour limits. What does the autoscaled version hide?

## C. Challenge (do at least one)

**C1.** (Mislead, then fix.) Two *independent* random walks of 200 steps (PY-6: cumulative sums of $\pm1$ steps). (a) Draw them on a `twinx` figure with limits you choose so they appear to move together, and write the false claim ("$X$ tracks $Y$") as the title. (b) Over 2,000 independent pairs, compute the Pearson correlation of the *levels* and of the *increments*; assert that the median $|r|$ of levels exceeds 0.3 and that of increments is below 0.08. (c) Draw the honest figure: a histogram of $r$ for levels and for increments, or the increments scatter with its $r$. State in one sentence why (a) misleads.

**C2.** (Direct labels that do not collide.) Write `label_lines(ax, lines, names, min_gap=4)` that labels each line at its right end (like `ph.label_end`) but nudges the labels vertically so that no two label boxes overlap, keeping the top-to-bottom order of the line ends. Test with eight lines whose right ends lie within 2% of each other. Assert with `Text.get_window_extent` that no two label boxes overlap and that the order is preserved. What does your method do when there is no room (more labels than the axes can hold)?

**C3.** (File size and time.) For $N = 10^3, 10^4, 10^5$ random points, save (i) `scatter` as PDF, (ii) `scatter(rasterized=True)` as PDF, (iii) `plot(".")` as PDF, (iv) `hexbin` as PDF and (v) `scatter` as PNG at 150 dpi, recording the file size and the time of `savefig`. Print a table. Assert: the vector scatter's PDF grows by at least $5\times$ from $10^3$ to $10^5$ points, the rasterised one by less than $3\times$, and the `hexbin` PDF by less than $2\times$. Then write three sentences: which you would use for a paper, for a slide, and for a figure with $10^7$ points, and what you would change about the data first.
