# PY-7 Solutions

Every code block was run, and the ones with `assert`s passed them; the outputs shown are the real ones (timings and file sizes come from one machine and will differ on yours). The blocks share one session, top to bottom, started in this folder: start your own file with the imports below. Figures are saved to `figures/` as `sol_*.png`; open them.

```python
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from matplotlib.ticker import FuncFormatter, MaxNLocator, MultipleLocator
from scipy import stats
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

import plot_helpers as ph

ph.use_house_style()
```

## A. Warm-up

**A1.** `subplots(1, 3)` returns `(fig, ndarray of shape (3,))`; `subplots(2, 2, sharey="row")` an array of shape `(2, 2)`. `fig.axes` lists the figure's Axes: 3 and 4. `sharey="row"` joins the $y$-axes of each *row* into one group, so `[0, 0]` shares with `[0, 1]` but not with `[1, 0]`.

```python
fig, axs = plt.subplots(1, 3)
fig2, ax2 = plt.subplots(2, 2, sharey="row")
print(type(axs).__name__, axs.shape, ax2.shape, len(fig.axes), len(fig2.axes))
same_row = ax2[0, 0].get_shared_y_axes().joined(ax2[0, 0], ax2[0, 1])
same_col = ax2[0, 0].get_shared_y_axes().joined(ax2[0, 0], ax2[1, 0])
print("shared with [0, 1]:", same_row, "| shared with [1, 0]:", same_col)
assert (same_row, same_col) == (True, False)
plt.close("all")
```
```text
ndarray (3,) (2, 2) 3 4
shared with [0, 1]: True | shared with [1, 0]: False
```

**A2.** `plt.xlabel` and `plt.title` act on the *current* axes, the one created last. So both land on the **second** figure and the first is bare. The fix holds each figure and its axes in variables:

```python
def plot_pair_buggy(x, y1, y2):
    plt.figure()
    plt.plot(x, y1)
    plt.figure()
    plt.plot(x, y2)
    plt.xlabel("x (m)")          # intended for the FIRST figure
    plt.title("both")


x = np.linspace(0, 1, 5)
plot_pair_buggy(x, x, x**2)
f1, f2 = (plt.figure(n) for n in plt.get_fignums())
print("first figure:", repr(f1.axes[0].get_xlabel()), repr(f1.axes[0].get_title()),
      "| second figure:", repr(f2.axes[0].get_xlabel()), repr(f2.axes[0].get_title()))
plt.close("all")


def plot_pair(x, y1, y2):
    fig1, ax1 = plt.subplots()
    ax1.plot(x, y1)
    ax1.set(xlabel="x (m)", title="both")
    fig2, ax2 = plt.subplots()
    ax2.plot(x, y2)
    ax2.set(title="both")
    return fig1, fig2


f1, f2 = plot_pair(x, x, x**2)
assert f1.axes[0].get_xlabel() == "x (m)" and f2.axes[0].get_xlabel() == ""
assert f1.axes[0].get_title() == f2.axes[0].get_title() == "both"
plt.close("all")
```
```text
first figure: '' '' | second figure: 'x (m)' 'both'
```

**A3.** Take logs of the model. (a) $\log_{10}y = \log_{10}3 + 2\log_{10}x$: `loglog`, slope 2. (b) $\log_{10}y = \log_{10}2 - (\log_{10}e/5)\,x$: `semilogy`, slope $-\log_{10}e/5 = -0.0869$ per unit of $x$. (c) $\log_{10}y = \log_{10}5 - \tfrac12\log_{10}x$: `loglog`, slope $-1/2$. (d) already linear in $\log_{10}x$: `semilogx`, slope 3 *per decade*.

```python
xx = np.linspace(1, 100, 300)
ident = lambda v: v
cases = {
    "(a) y = 3 x^2": (3 * xx**2, np.log10, np.log10, "loglog", 2.0),
    "(b) y = 2 exp(-x/5)": (2 * np.exp(-xx / 5), ident, np.log10, "semilogy", -np.log10(np.e) / 5),
    "(c) y = 5 / sqrt(x)": (5 / np.sqrt(xx), np.log10, np.log10, "loglog", -0.5),
    "(d) y = 2 + 3 log10(x)": (2 + 3 * np.log10(xx), np.log10, ident, "semilogx", 3.0),
}
for name, (yy, fx, fy, kind, slope) in cases.items():
    X_, Y_ = fx(xx), fy(yy)
    a, b = np.polyfit(X_, Y_, 1)
    r2 = np.corrcoef(X_, Y_)[0, 1] ** 2
    assert r2 > 0.9999 and np.isclose(a, slope)
    print(f"{name:24s} straight on {kind:8s}, slope {a:8.4f}; R^2 on linear axes would be {np.corrcoef(xx, yy)[0, 1] ** 2:.3f}")
```
```text
(a) y = 3 x^2            straight on loglog  , slope   2.0000; R^2 on linear axes would be 0.939
(b) y = 2 exp(-x/5)      straight on semilogy, slope  -0.0869; R^2 on linear axes would be 0.272
(c) y = 5 / sqrt(x)      straight on loglog  , slope  -0.5000; R^2 on linear axes would be 0.530
(d) y = 2 + 3 log10(x)   straight on semilogx, slope   3.0000; R^2 on linear axes would be 0.815
```

**A4.** (i) The default `imshow(A)` has extent $(-0.5, 5.5, 3.5, -0.5)$ and row 0 at the top, so `A[0, 3]` is the cell centred at $(x, y) = (3, 0)$, at the **top** of the page; its mirror image is $(3, 3)$. (ii) With `origin="lower"` and `extent=(0, 6, 0, 4)`, cell `[i, j]` fills $[j, j+1]\times[i, i+1]$, so `A[0, 3]` is centred at $(3.5, 0.5)$ at the **bottom**; the mirror is $(3.5, 3.5)$.

```python
A = np.zeros((4, 6))
A[0, 3] = 1.0
bright, dark = plt.get_cmap("viridis")(1.0), plt.get_cmap("viridis")(0.0)
for name, kw, cell, mirror in [("default", {}, (3.0, 0.0), (3.0, 3.0)),
                               ("origin='lower', extent=(0, 6, 0, 4)", {"origin": "lower", "extent": (0, 6, 0, 4)}, (3.5, 0.5), (3.5, 3.5))]:
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(A, **kw)
    ax.grid(False)
    c_cell, c_mirror = ph.rgba_at(ax, *cell), ph.rgba_at(ax, *mirror)
    (_, y_cell), (_, y_mirror) = ax.transData.transform([cell, mirror])
    assert np.allclose(c_cell, bright, atol=0.02) and np.allclose(c_mirror, dark, atol=0.02)
    print(f"{name}: extent {tuple(float(v) for v in im.get_extent())}; bright cell at {cell}, {'top' if y_cell > y_mirror else 'bottom'} of the page")
    plt.close(fig)
```
```text
default: extent (-0.5, 5.5, 3.5, -0.5); bright cell at (3.0, 0.0), top of the page
origin='lower', extent=(0, 6, 0, 4): extent (0.0, 6.0, 0.0, 4.0); bright cell at (3.5, 0.5), bottom of the page
```

**A5.** `MultipleLocator` computes ticks by floating-point arithmetic, so "zero" can arrive as $-10^{-16}$; `np.degrees` turns that into $-5.7\times10^{-15}$, and `format(..., ".0f")` of a tiny negative number is `"-0"`. Rounding to an `int` first removes the sign.

```python
def deg(v, pos=None):
    return f"{int(round(np.degrees(v)))}°"


naive = lambda v, pos=None: f"{np.degrees(v):.0f}°"
ticks = np.pi / 6 * np.arange(13)
labels = [deg(v) for v in ticks]
assert labels[3] == "90°" and labels[-1] == "360°" and naive(-1e-16) == "-0°" and deg(-1e-16) == "0°"
print(labels)
print("naive(-1e-16) =", naive(-1e-16), "| deg(-1e-16) =", deg(-1e-16))
fig, ax = plt.subplots()
tt = np.linspace(0, 2 * np.pi, 300)
ax.plot(tt, np.sin(tt))
ax.xaxis.set_major_locator(MultipleLocator(np.pi / 6))
ax.xaxis.set_major_formatter(FuncFormatter(deg))
ax.set(xlabel="angle (degrees)", ylabel="$\\sin\\theta$ (dimensionless)", xlim=(0, 2 * np.pi))
fmt = ax.xaxis.get_major_formatter()
print("ticks in view:", [fmt(v) for v in ax.xaxis.get_majorticklocs() if -1e-9 <= v <= 2 * np.pi + 1e-9][:5], "...")
ph.savefig(fig, "sol_a5")
```
```text
['0°', '30°', '60°', '90°', '120°', '150°', '180°', '210°', '240°', '270°', '300°', '330°', '360°']
naive(-1e-16) = -0° | deg(-1e-16) = 0°
ticks in view: ['0°', '30°', '60°', '90°', '120°'] ...
```

**A6.** $\text{Normalize}(2, 6)(5) = (5-2)/4 = 0.75$; $\text{LogNorm}(1, 100)(10) = \ln 10/\ln 100 = 0.5$; $\text{Normalize}(-3, 3)(0) = 0.5$; `TwoSlopeNorm` maps $[v_{\min}, v_{\text{center}}]$ to $[0, 0.5]$ and $[v_{\text{center}}, v_{\max}]$ to $[0.5, 1]$ linearly, so $2.5 \mapsto 0.5 + 0.5\cdot2.5/5 = 0.75$.

```python
vals = [float(Normalize(2, 6)(5)), float(LogNorm(1, 100)(10)), float(Normalize(-3, 3)(0)), float(TwoSlopeNorm(vmin=-1, vcenter=0, vmax=5)(2.5))]
print(vals)
assert np.allclose(vals, [0.75, 0.5, 0.5, 0.75])
```
```text
[0.75, 0.5, 0.5, 0.75]
```

(i) The anomaly is signed with a meaningful zero: a **diverging** map (`RdBu_r`) with `Normalize(-m, m)`, so white means zero. (ii) Wind direction is an angle: a **cyclic** map (`twilight`), so 359° and 1° get nearly the same colour. (iii) Density over five decades: a **sequential** map (`viridis`) with `LogNorm(1, 1e5)`, since equal steps are ratios. (iv) Cluster labels are unordered: a **categorical** palette (Okabe–Ito), never a sequential map, which would suggest cluster 5 is "more" than cluster 0.

## B. Core

**B1.** One step of each method for $y' = f(y) = -y$ with step $h$: Euler $y + hf(y)$; Heun $y + \tfrac h2[f(y) + f(y + hf(y))]$; RK4 $y + \tfrac h6(k_1 + 2k_2 + 2k_3 + k_4)$ with $k_1 = f(y)$, $k_2 = f(y + \tfrac h2k_1)$, $k_3 = f(y + \tfrac h2k_2)$, $k_4 = f(y + hk_3)$. With $n$ listed *descending*, $h$ ascends, so the last point of each series is its largest $h$, at the right edge, where `ph.label_end` puts the name. `plot_convergence` draws the data first, so the data line is `ax.lines[before]`.

```python
f = lambda y: -y
euler = lambda y, h: y + h * f(y)
heun = lambda y, h: y + h / 2 * (f(y) + f(y + h * f(y)))


def rk4(y, h):
    k1 = f(y)
    k2 = f(y + h / 2 * k1)
    k3 = f(y + h / 2 * k2)
    k4 = f(y + h * k3)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def global_error(step, n, T=2.0):
    y, h = 1.0, T / n
    for _ in range(n):
        y = step(y, h)
    return abs(y - np.exp(-T))


ns = np.array([192, 128, 96, 64, 48, 32, 24, 16, 12, 8])
h = 2.0 / ns
fig, ax = plt.subplots()
slopes = {}
for name, step, order in [("Euler", euler, 1), ("Heun", heun, 2), ("RK4", rk4, 4)]:
    err = np.array([global_error(step, n) for n in ns])
    before = len(ax.lines)
    slopes[name] = ph.plot_convergence(ax, h, err, order)
    ph.label_end(ax, ax.lines[before], name)
    assert abs(slopes[name] - order) < 0.1
ax.set(ylabel="global error at $t = 2$ (dimensionless)")
ph.savefig(fig, "sol_b1")
print({k: round(v, 3) for k, v in slopes.items()})
```
```text
{'Euler': 1.011, 'Heun': 2.054, 'RK4': 4.057}
```

RK4's slope is furthest from its order (4.057 against 4). The fit uses $h$ from 0.01 to 0.25, where the next term of the error expansion ($h^5$) is still a visible fraction of the leading $h^4$ term, the same effect as Simpson's 4.018 in the notes. A fitted slope gives the order *only* where the leading term dominates: report it with its $h$ range.

**B2.** The fits use the branches the exercise names; the minimum is where the two error sources cross: $\tfrac h2|f''| \approx \varepsilon|f|/h$ gives $h \approx \sqrt{2\varepsilon|f|/|f''|} = \sqrt{2\varepsilon} \approx 2.1\times10^{-8}$ here ($|f| = |f''| = \sin 1$), within a factor of 2 of the grid point $10^{-8}$ where the minimum lands.

```python
hs = 10.0 ** -np.arange(1, 17)
err = np.abs((np.sin(1 + hs) - np.sin(1)) / hs - np.cos(1))
k = np.argmin(err)
right, left = hs >= 1e-6, hs <= 1e-10
s_r, b_r, _ = ph.loglog_fit(hs[right], err[right])
s_l, b_l, _ = ph.loglog_fit(hs[left], err[left])
eps = np.finfo(float).eps
print(f"minimum at h = {hs[k]:.0e} (sqrt(eps) = {np.sqrt(eps):.1e}), error {err[k]:.1e}; slopes: right {s_r:.3f}, left {s_l:.3f}")
assert abs(s_r - 1) < 0.05 and abs(s_l + 1) < 0.3 and abs(np.log10(hs[k] / np.sqrt(eps))) < 1

fig, ax = plt.subplots()
ax.loglog(hs, err, "o")
ax.loglog(hs[right], 10**b_r * hs[right] ** s_r, "-", lw=0.8, color="0.4")
ax.loglog(hs[left], 10**b_l * hs[left] ** s_l, "-", lw=0.8, color="0.4")
ph.add_slope_guide(ax, 1, 1e-5, 1e-8)
ph.add_slope_guide(ax, -1, 1e-15, 1e-5)
ax.set(xlabel="step size $h$", ylabel="absolute error of $D_h$", title="forward difference of $\\sin$ at $x=1$")
ph.savefig(fig, "sol_b2")
print("errors at h = 1e-10 ... 1e-14:", np.array2string(err[9:14], precision=1))
```
```text
minimum at h = 1e-08 (sqrt(eps) = 1.5e-08), error 3.0e-09; slopes: right 1.001, left -1.108
errors at h = 1e-10 ... 1e-14: [5.8e-08 1.2e-06 4.3e-05 7.3e-04 3.7e-03]
```

The left branch is noisy because $\sin(1+h) - \sin 1$ is a difference of two doubles near 0.84, which are spaced $1.1\times10^{-16}$ apart. For $h \sim 10^{-10}$ the numerator is a multiple of that spacing with only a few significant digits, so the error is a staircase in $h$, not a smooth $\varepsilon/h$; at $h = 10^{-16}$, $1+h = 1$ exactly and $D_h = 0$, error $= \cos 1$. The fitted $-1.1$ is the trend through the staircase. This is why the optimal step for first derivatives is about $\sqrt\varepsilon$, not "as small as possible" (LA-13, PY-8).

**B3.** The field is $(\dot x, \dot y) = (y, \mu(1-x^2)y - x)$. Both trajectories are integrated to $t = 60$ and their late amplitude is measured over $t > 40$, when transients have decayed (the cycle's period is about 6.7).

```python
mu = 1.0
vdp = lambda t, s: [s[1], mu * (1 - s[0] ** 2) * s[1] - s[0]]
tt = np.linspace(0, 60, 3000)
trajs = {}
for name, s0 in {"inside (0.1, 0)": (0.1, 0.0), "outside (3, 3)": (3.0, 3.0)}.items():
    sol = solve_ivp(vdp, (0, 60), s0, rtol=1e-9, atol=1e-11, dense_output=True)
    trajs[name] = sol.sol(tt)
    amp = np.abs(trajs[name][0][tt > 40]).max()
    print(f"{name:16s} late max|x| = {amp:.4f}")
    assert abs(amp - 2.0) < 0.05

L = (3.5, 4.0)
xg, yg = np.meshgrid(np.linspace(-L[0], L[0], 25), np.linspace(-L[1], L[1], 15))
u, v = yg, mu * (1 - xg**2) * yg - xg
nrm = np.maximum(np.hypot(u, v), 1e-9)                        # unit length in data units
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
ax1.quiver(xg, yg, u / nrm, v / nrm, angles="xy", pivot="mid", color="0.3")   # angles="xy": tangent on the page (Example 4.7)
xf, yf = np.meshgrid(np.linspace(-L[0], L[0], 200), np.linspace(-L[1], L[1], 200))
ax2.streamplot(xf, yf, yf, mu * (1 - xf**2) * yf - xf, color="0.65", density=1.2, linewidth=0.8)
for (name, (xt, yt)), c in zip(trajs.items(), ph.OKABE_ITO[1:]):
    ax2.plot(xt, yt, color=c, lw=1.6, label=name)
ax2.legend(loc="lower right")
for a in (ax1, ax2):
    a.set(xlim=(-L[0], L[0]), ylim=(-L[1], L[1]), xlabel="$x$ (dimensionless)")
ax1.set(ylabel="$\\dot x$ (dimensionless)", title="quiver: direction only")
ax2.set_title("two trajectories, one limit cycle")
ph.savefig(fig, "sol_b3")
```
```text
inside (0.1, 0)  late max|x| = 2.0086
outside (3, 3)   late max|x| = 2.0086
```

Open `figures/sol_b3.png`: the trajectory from inside spirals out and the one from outside spirals in, onto the same closed curve, whose amplitude in $x$ is about 2.0 (for $\mu = 1$ it is 2.0086). Unlike the pendulum there is no conserved energy, which is why a closed curve here is an *attractor* and not one of a family.

**B4.** (a) For $n = 5$ the SE is itself random and often too small, so $\pm1.96\,\mathrm{SE}$ undercovers: the exact coverage is $P(|t_4| < 1.96) = 0.878$. The $t$ interval has the 95% it claims.

```python
rng = np.random.default_rng(4)
K, n = 4000, 5
x = rng.normal(0, 3, (K, n))
xbar, se = x.mean(axis=1), x.std(axis=1, ddof=1) / np.sqrt(n)
cov_196 = np.mean(np.abs(xbar) < 1.96 * se)
cov_t = np.mean(np.abs(xbar) < stats.t.ppf(0.975, n - 1) * se)
print(f"coverage of the true mean 0 at n = 5: +-1.96 SE {cov_196:.3f}, t interval {cov_t:.3f} (exact: {2 * stats.t.cdf(1.96, 4) - 1:.3f}, 0.950)")
assert abs(cov_t - 0.95) < 0.02 and cov_196 < 0.92
```
```text
coverage of the true mean 0 at n = 5: +-1.96 SE 0.866, t interval 0.948 (exact: 0.878, 0.950)
```

(b) A bandwidth factor $f$ in `gaussian_kde(sample, bw_method=f)` sets the kernel SD to $f\cdot s$ (4.1). Too small a factor traces the sample's noise, too large smears the shape; the total-variation distance to the truth has a minimum between.

```python
true = stats.gamma(4, scale=1.5)
sample = true.rvs(1500, random_state=np.random.default_rng(0))
xs = np.linspace(-6, 24, 1500)
factors = np.logspace(np.log10(0.02), np.log10(1.5), 25)
tv = np.array([0.5 * np.trapezoid(np.abs(stats.gaussian_kde(sample, bw_method=f)(xs) - true.pdf(xs)), xs) for f in factors])
scott = len(sample) ** -0.2
best = factors[np.argmin(tv)]
print(f"best factor {best:.3f} (TV {tv.min():.3f}); Scott {scott:.3f}; TV at 0.02: {tv[0]:.3f}, at 1.5: {tv[-1]:.3f}")
assert abs(np.log(best / scott)) < np.log(3) and tv[0] > tv.min() and tv[-1] > tv.min()
fig, ax = plt.subplots()
ax.semilogx(factors, tv, "o-")
ax.axvline(scott, color="k", ls=":", lw=0.8)
ax.annotate("Scott: $n^{-1/5}$", xy=(scott, tv.max() * 0.6), xytext=(6, 0), textcoords="offset points")
ax.set(xlabel="bandwidth factor $f$ (kernel SD = $f\\,s$)", ylabel="total-variation distance to the truth")
ph.savefig(fig, "sol_b4")
```
```text
best factor 0.207 (TV 0.018); Scott 0.232; TV at 0.02: 0.084, at 1.5: 0.308
```

Scott's rule is within 12% of the best factor here, but it is tuned for roughly normal data; for a sample with a sharp feature the best factor is smaller.

**B5.** Sequential: monotone with a range over 50. Diverging: one ridge near the middle and matching ends. Anything else is `"other"`.

```python
def lstar(rgb):
    """CIE L* of sRGB colours (rows of `rgb`, floats in [0, 1]), eq. (3.7)."""
    rgb = np.atleast_2d(rgb)[:, :3]
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    Y = lin @ np.array([0.2126, 0.7152, 0.0722])
    f = np.where(Y > 216 / 24389, np.cbrt(Y), (24389 / 27 * Y + 16) / 116)
    return 116 * f - 16


def lightness(name):
    return lstar(plt.get_cmap(name)(np.linspace(0, 1, 256)))


def classify(name):
    L = lightness(name)
    d = np.diff(L)
    if (np.all(d > 0) or np.all(d < 0)) and L.max() - L.min() > 50:
        return "sequential"
    peaks = find_peaks(L, prominence=2.0)[0]
    if len(peaks) == 1 and abs(peaks[0] - 127.5) <= 25 and abs(L[0] - L[-1]) < 5:
        return "diverging"
    return "other"


expected = {**dict.fromkeys(["viridis", "cividis", "magma", "plasma", "inferno", "gray"], "sequential"),
            **dict.fromkeys(["coolwarm", "RdBu_r"], "diverging"), **dict.fromkeys(["jet", "turbo"], "other")}
for name, want in expected.items():
    L = lightness(name)
    print(f"{name:9s} {classify(name):10s} L* ends {L[0]:5.1f}, {L[-1]:5.1f}; lightest colour at {L.argmax() / 255:.2f} of the map")
    assert classify(name) == want
```
```text
viridis   sequential L* ends  14.9,  90.9; lightest colour at 1.00 of the map
cividis   sequential L* ends  13.9,  91.2; lightest colour at 1.00 of the map
magma     sequential L* ends   0.1,  97.8; lightest colour at 1.00 of the map
plasma    sequential L* ends  15.4,  94.4; lightest colour at 1.00 of the map
inferno   sequential L* ends   0.1,  97.9; lightest colour at 1.00 of the map
gray      sequential L* ends   0.0, 100.0; lightest colour at 1.00 of the map
coolwarm  diverging  L* ends  37.7,  37.7; lightest colour at 0.50 of the map
RdBu_r    diverging  L* ends  20.1,  20.0; lightest colour at 0.50 of the map
jet       other      L* ends  12.9,  25.4; lightest colour at 0.64 of the map
turbo     other      L* ends  12.0,  24.5; lightest colour at 0.49 of the map
```

`jet` puts its lightest colour (yellow, $L^*\approx96$) at 64% of the way along the map and `turbo` puts a peak ($L^*\approx91$) in the middle, at 49%. For one smooth bump, in either map the *brightest* ring is at an intermediate data value and the maximum is drawn darker (`jet` ends at $L^* = 25$): the eye reads a ridge or a hole that the data do not have. `turbo` fails the diverging test only because its two ends differ in lightness (12 and 24.5), so it is not a symmetric map either; it was designed as a better rainbow, not as a diverging map.

**B6.** `make_fig` takes the size; sizes of text and lines come from an `rc_context`, so nothing leaks. The embedded font type is visible in the PDF bytes: `pdf.fonttype` 42 stores the fonts as TrueType (a `CIDFontType2` object), the default 3 as `Type3`.

```python
def make_fig(figsize, max_ticks=None):
    true = stats.gamma(4, scale=1.5)
    sample = true.rvs(1500, random_state=np.random.default_rng(0))
    mu, sd = stats.norm.fit(sample)
    kde = stats.gaussian_kde(sample)
    xs = np.linspace(-6, 24, 600)
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(sample, bins=30, range=(0, 20), density=True, color="0.8", edgecolor="white")
    ax.plot(xs, stats.norm.pdf(xs, mu, sd), label="Gaussian fit")
    ax.plot(xs, kde(xs), label="KDE")
    ax.set(xlim=(-3, 22), xlabel="waiting time (s)", ylabel="probability density (1/s)")
    ax.legend()
    if max_ticks:
        ax.xaxis.set_major_locator(MaxNLocator(max_ticks))
        ax.yaxis.set_major_locator(MaxNLocator(max_ticks))
    return fig, ax


def ticks_in_view(axis, lim):
    return len([t for t in axis.get_majorticklocs() if lim[0] <= t <= lim[1]])


before = plt.rcParams["font.size"]
with plt.rc_context({"font.size": 16, "lines.linewidth": 3.5}):
    fig, ax = make_fig((10, 5.6), max_ticks=5)
    n_x, n_y = ticks_in_view(ax.xaxis, ax.get_xlim()), ticks_in_view(ax.yaxis, ax.get_ylim())
    assert tuple(fig.get_size_inches()) == (10.0, 5.6) and n_x <= 6 and n_y <= 6
    slide = ph.savefig(fig, "sol_b6_slide")
with plt.rc_context({"font.size": 8, "pdf.fonttype": 42}):
    fig, ax = make_fig((3.4, 2.6))
    assert tuple(fig.get_size_inches()) == (3.4, 2.6)
    paper = ph.savefig(fig, "sol_b6_paper", fmt="pdf")
data = paper.read_bytes()
assert data[:4] == b"%PDF" and b"CIDFontType2" in data and b"/Type3" not in data and plt.rcParams["font.size"] == before
print(f"slide: {slide.name}, ticks per axis {n_x} and {n_y}; paper: {paper.name}, {len(data) / 1e3:.0f} kB, TrueType fonts (no Type 3); font.size after = {plt.rcParams['font.size']}")
```
```text
slide: sol_b6_slide.png, ticks per axis 5 and 5; paper: sol_b6_paper.pdf, 13 kB, TrueType fonts (no Type 3); font.size after = 10.0
```

The paper version cuts the size, the line weight and the tick density to fit a 3.4-inch column, where the caption and the surrounding text carry the explanation, and it is vector so it prints sharply at any size. The slide keeps large type, thick lines and few ticks because the audience is far away and has seconds, not minutes, so it also has to carry its own message in the title and label.

**B7.** The mass check is a Riemann sum $\sum u\,\Delta x\,\Delta y$; the Gaussian is well resolved ($\sigma_0 = 0.3$ against $\Delta x = 0.075$) and the box is wide enough that only $10^{-5}$ is lost at $t = 8$ ($\sigma = 1.3$, edge at $4.6\sigma$). One `Normalize` object is passed to every panel and one colorbar is built from the last image; `len(fig.axes) == 7` is six panels plus the colorbar's own Axes.

```python
D, s0 = 0.1, 0.3
g = np.linspace(-6, 6, 161)
dxy = g[1] - g[0]
GX, GY = np.meshgrid(g, g)
times = [0, 0.5, 1, 2, 4, 8]


def snapshot(t):
    s2 = s0**2 + 2 * D * t
    return np.exp(-(GX**2 + GY**2) / (2 * s2)) / (2 * np.pi * s2)


U = [snapshot(t) for t in times]
masses = [u.sum() * dxy**2 for u in U]
assert all(abs(m - 1) < 1e-3 for m in masses)
extent = (g[0] - dxy / 2, g[-1] + dxy / 2, g[0] - dxy / 2, g[-1] + dxy / 2)         # pixel EDGES
norm = Normalize(0, U[0].max())                                                      # ONE norm for all six panels
fig, axs = plt.subplots(2, 3, figsize=(10, 6.5), sharex=True, sharey=True)
images = []
for ax, u, t in zip(axs.flat, U, times):
    images.append(ax.imshow(u, origin="lower", extent=extent, norm=norm))
    ax.grid(False)
    ax.set(title=f"$t = {t:g}$")
for ax in axs[1]:
    ax.set_xlabel("$x$ (length)")
for ax in axs[:, 0]:
    ax.set_ylabel("$y$ (length)")
fig.colorbar(images[-1], ax=axs, shrink=0.8, label="concentration $u$ (1/length$^2$)")
ph.savefig(fig, "sol_b7")
assert len(fig.axes) == 7 and len({im.get_clim() for im in images}) == 1
print("masses:", np.round(masses, 5), "| shared clim:", tuple(round(float(v), 3) for v in images[0].get_clim()))

fig, axs = plt.subplots(2, 3, figsize=(10, 6.5), sharex=True, sharey=True)
auto = [ax.imshow(u, origin="lower", extent=extent) for ax, u in zip(axs.flat, U)]       # no shared norm: each autoscales
print("autoscaled clims:", [tuple(round(float(v), 3) for v in im.get_clim()) for im in auto])
plt.close(fig)
```
```text
masses: [1.      1.      1.      1.      1.      0.99999] | shared clim: (0.0, 1.768)
autoscaled clims: [(0.0, 1.768), (0.0, 0.838), (0.0, 0.549), (0.0, 0.325), (0.0, 0.179), (0.0, 0.094)]
```

With six different colour limits, every panel's brightest pixel is the same yellow, so the six pictures look like one blob of *constant* height getting wider. The data say the peak falls from 1.77 to about 0.05 while the total stays 1; the shared norm shows that (the late panels are nearly flat), and the shared colorbar is only valid because the norm is shared.

## C. Challenge

**C1.** (a) With `twinx`, each axis autoscales to its own series, so *any* two series fill the frame and wiggle together visually; the pair below is the most positively correlated of 2,000. (b) Each random walk is dominated by a few slow, low-frequency wiggles, so the effective number of independent points is far below 200 and the correlation of two independent *levels* is widely spread; increments are independent noise, so their $r \approx N(0, 1/\sqrt{199})$.

```python
rng = np.random.default_rng(7)
steps = rng.choice([-1, 1], size=(2000, 2, 200))
walks = steps.cumsum(axis=2)


def corr(a, b):
    a, b = a - a.mean(axis=1, keepdims=True), b - b.mean(axis=1, keepdims=True)
    return (a * b).sum(axis=1) / np.sqrt((a * a).sum(axis=1) * (b * b).sum(axis=1))


r_lev, r_inc = corr(walks[:, 0], walks[:, 1]), corr(steps[:, 0], steps[:, 1])
print(f"median |r|: levels {np.median(np.abs(r_lev)):.3f}, increments {np.median(np.abs(r_inc)):.3f}; largest |r| of levels {np.abs(r_lev).max():.3f}")
assert np.median(np.abs(r_lev)) > 0.3 and np.median(np.abs(r_inc)) < 0.08

k = np.argmax(r_lev)                                  # the most positively correlated of the 2,000 independent pairs
print(f"the pair drawn below has r = {r_lev[k]:.3f}")
t = np.arange(1, 201)
fig, ax = plt.subplots()
(l1,) = ax.plot(t, walks[k, 0], color=ph.OKABE_ITO[0])
ax2 = ax.twinx()
(l2,) = ax2.plot(t, walks[k, 1], color=ph.OKABE_ITO[1])
ax.set(xlabel="step $t$", ylabel="walk X (steps)", title="MISLEADING: \"X tracks Y\"")
ax2.set_ylabel("walk Y (steps)")
ph.savefig(fig, "sol_c1_misleading")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6), sharex=True, sharey=True)
bins = np.linspace(-1, 1, 41)
a1.hist(r_lev, bins=bins, color=ph.OKABE_ITO[0])
a2.hist(r_inc, bins=bins, color=ph.OKABE_ITO[1])
a1.set(title="levels of two independent walks", xlabel="Pearson $r$ (dimensionless)", ylabel="pairs (of 2,000)")
a2.set(title="their increments", xlabel="Pearson $r$ (dimensionless)")
ph.savefig(fig, "sol_c1_honest")
```
```text
median |r|: levels 0.391, increments 0.047; largest |r| of levels 0.948
the pair drawn below has r = 0.940
```

The misleading figure claims "X tracks Y", but the pair is *independent* by construction (`corr` of increments is 0.05 in the median) and was picked as the most positively correlated of 2,000. Two defects compound: autoscaled twin axes erase the scale of each series, and a correlation of levels of two trending series carries almost no evidence of a relationship. The honest figure shows the whole distribution of $r$ under independence: a typical pair of random walks has $|r| \approx 0.4$ by chance, so such an $r$ is not evidence.

**C2.** For $n$ labels of equal height $d$ sorted by desired position $p_1 \le \dots \le p_n$, the non-overlap constraint $q_{i+1} \ge q_i + d$ becomes $z_{i+1} \ge z_i$ for $z_i = q_i - id$. Minimising $\sum(q_i - p_i)^2$ is then *isotonic regression* of $p_i - id$, solved by pool-adjacent-violators; the answer is added back as $q_i = z_i + id$. That moves labels as little as possible and keeps their order.

```python
def isotonic(y):
    """Least-squares non-decreasing fit to `y` by pool-adjacent-violators."""
    blocks = []                                    # each block: [mean, count]
    for v in y:
        blocks.append([float(v), 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            (m2, c2), (m1, c1) = blocks.pop(), blocks.pop()
            blocks.append([(m1 * c1 + m2 * c2) / (c1 + c2), c1 + c2])
    return np.concatenate([[m] * c for m, c in blocks])


def label_lines(ax, lines, names, min_gap=4.0, fontsize=9):
    """Label each line at its right end; nudge labels vertically so no two overlap, keeping order."""
    fig = ax.figure
    fig.canvas.draw()
    d = fontsize * 1.2 * fig.dpi / 72 + min_gap                      # label pitch in pixels
    ends = np.array([ax.transData.transform(ln.get_xydata()[-1]) for ln in lines])
    order = np.argsort(ends[:, 1])
    if d * len(lines) > ax.bbox.height:
        raise ValueError("no room: more labels than fit in the axes height")
    q = isotonic(ends[order, 1] - d * np.arange(len(lines))) + d * np.arange(len(lines))
    texts = [None] * len(lines)
    for rank, i in enumerate(order):
        dy_points = (q[rank] - ends[i, 1]) * 72 / fig.dpi
        texts[i] = ax.annotate(names[i], xy=lines[i].get_xydata()[-1], xytext=(5, dy_points), textcoords="offset points",
                               color=lines[i].get_color(), fontsize=fontsize, va="center")
    return texts


fig, ax = plt.subplots()
t = np.linspace(0, 1, 50)
lines = [ax.plot(t, 1 + (0.0025 * i - 0.01) * t + 0.1 * (7 - i) * (1 - t), color=f"C{i}")[0] for i in range(8)]
names = [f"series {i}" for i in range(8)]
ax.set(xlabel="time $t$ (s)", ylabel="signal (V)", xlim=(0, 1.15))
texts = label_lines(ax, lines, names)
fig.canvas.draw()
boxes = [tx.get_window_extent() for tx in texts]
overlaps = sum(boxes[i].overlaps(boxes[j]) for i in range(8) for j in range(i + 1, 8))
y_end = [ln.get_xydata()[-1][1] for ln in lines]
order_ok = list(np.argsort([b.y0 for b in boxes])) == list(np.argsort(y_end))
print("ends span", round(max(y_end) - min(y_end), 4), "of", round(np.ptp(ax.get_ylim()), 2), "| overlapping pairs:", overlaps, "| order kept:", order_ok)
assert overlaps == 0 and order_ok
ph.savefig(fig, "sol_c2")
try:
    label_lines(ax, lines * 12, names * 12)
except ValueError as e:
    print("ValueError:", e)
```
```text
ends span 0.0175 of 0.78 | overlapping pairs: 0 | order kept: True
ValueError: no room: more labels than fit in the axes height
```

Labels are placed by an offset in points relative to the line end, so they stay attached when the figure is saved at another dpi. With more labels than the axes can hold, no placement exists; the function raises instead of drawing a pile. The real remedy is fewer series, or a legend, or a small-multiples layout (§3.7).

**C3.** Sizes are exact for one run; times depend on the machine.

```python
rng = np.random.default_rng(11)
variants = {
    "scatter, vector PDF": ("pdf", lambda ax, x, y: ax.scatter(x, y, s=2)),
    "scatter, rasterized PDF": ("pdf", lambda ax, x, y: ax.scatter(x, y, s=2, rasterized=True)),
    "plot('.'), vector PDF": ("pdf", lambda ax, x, y: ax.plot(x, y, ".", ms=2)),
    "hexbin, vector PDF": ("pdf", lambda ax, x, y: ax.hexbin(x, y, gridsize=40)),
    "scatter, PNG 150 dpi": ("png", lambda ax, x, y: ax.scatter(x, y, s=2)),
}
sizes, times = {}, {}
for N in (10**3, 10**4, 10**5):
    x, y = rng.normal(size=(2, N))
    for name, (fmt, draw) in variants.items():
        fig, ax = plt.subplots(figsize=(4, 3))
        draw(ax, x, y)
        ax.set(xlabel="$x$ (arb.)", ylabel="$y$ (arb.)")
        t0 = time.perf_counter()
        path = ph.savefig(fig, f"sol_c3_{fmt}", fmt=fmt, tight=False)
        times[name, N], sizes[name, N] = time.perf_counter() - t0, path.stat().st_size / 1e3
print(f"{'variant':26s}" + "".join(f"{'N=' + format(N, ','):>22s}" for N in (10**3, 10**4, 10**5)))
for name in variants:
    print(f"{name:26s}" + "".join(f"{sizes[name, N]:12.0f} kB {times[name, N]:5.2f} s" for N in (10**3, 10**4, 10**5)))
grow = lambda name: sizes[name, 10**5] / sizes[name, 10**3]
assert grow("scatter, vector PDF") >= 5 and grow("scatter, rasterized PDF") < 3 and grow("hexbin, vector PDF") < 2
print("growth 1e3 -> 1e5:", {n: round(grow(n), 2) for n in variants})
```
```text
variant                                  N=1,000              N=10,000             N=100,000
scatter, vector PDF                 25 kB  0.11 s         162 kB  0.24 s        1518 kB  1.43 s
scatter, rasterized PDF             25 kB  0.14 s          31 kB  0.26 s          24 kB  0.21 s
plot('.'), vector PDF               25 kB  0.11 s         162 kB  0.24 s        1518 kB  1.31 s
hexbin, vector PDF                  11 kB  0.20 s          14 kB  0.22 s          16 kB  0.19 s
scatter, PNG 150 dpi                43 kB  0.09 s          63 kB  0.11 s          44 kB  0.13 s
growth 1e3 -> 1e5: {'scatter, vector PDF': 59.8, 'scatter, rasterized PDF': 0.95, "plot('.'), vector PDF": 59.93, 'hexbin, vector PDF': 1.48, 'scatter, PNG 150 dpi': 1.04}
```

For a **paper**: a vector PDF with `rasterized=True` on the dense layer (axes, text and ticks stay vector and crisp, the cloud is one image), or `hexbin`/a 2D histogram if density is the message. For a **slide**: PNG at 150 dpi, since a slide is a raster anyway and the PNG's size is nearly independent of $N$. For **$10^7$ points**: no vector layer can hold them and even a raster scatter is slow and saturates to a solid blob; bin the data first (`np.histogram2d` or `hexbin`) and draw the counts, with a log colour norm if the range is wide. Change the *data* (aggregate) before changing the file format.
