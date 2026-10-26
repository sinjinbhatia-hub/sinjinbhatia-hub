# PY-7: Matplotlib and Scientific Communication: Figures That Say Something True

**Time:** ~7 hours · **Prerequisites:** PY-6 (PY-5's integrators supply the data for §1) · **Unlocks:** every figure in the curriculum; take it before CM-1's project

> **Conventions in this module.** Python 3.11+; every output below was pasted from a real run (Python 3.13, NumPy 2.5.3, SciPy 1.18, Matplotlib 3.11.2). `import matplotlib.pyplot as plt` and the **object-oriented API only**: `fig, ax = plt.subplots()`, then `ax.plot(...)`; pyplot is for creating and closing figures and for global settings. Headless: run with `MPLBACKEND=Agg`; scripts call `fig.savefig`, never `plt.show()`. Every axis label carries a quantity and a unit. Default colour map: viridis; default cycle: Okabe–Ito. "Slope" on a log-log plot means slope in $\log_{10}$–$\log_{10}$ coordinates. Code blocks run top to bottom in one session started in this folder, where `plot_helpers.py` sits beside `notes.md` (§3.7); the `text` block after a block is its real output. Images are gitignored (`figures/`), so figures are described in words and in numbers; open `figures/*.png` as you go.

## 0. Where this sits

PY-6 gave you arrays; PY-5 §3.10 gave you `trapezoid` and `simpson` and their orders, 2 and 4. This module shows that to a skeptic in one picture. After it, every project ends in a figure you can defend: CM-1's phase portraits (§4.7), EM-1's vector fields, EM-2's potentials (§4.5), F3's spectra on log axes, F4's histograms (§4.3), F5's Monte Carlo error bars (§4.4), ST-1's intervals. PY-9's `DataFrame.plot` is Matplotlib underneath.

## 1. The problem that forced this idea

A figure can hide a result as easily as show it. The forcing case: PY-5's `trapezoid` and `simpson` integrate $\sin$ over $[0,\pi]$ (exact value 2), and the claim is that their errors shrink like $h^2$ and $h^4$. You have the errors for ten step sizes (vectorized here; `project.py` asserts they match PY-5's `integrators.py`) and must convince a skeptic with one picture.

```python
import numpy as np
import matplotlib.pyplot as plt
import plot_helpers as ph          # this module's helpers (built in §3.7); copy the file to reuse them

ph.use_house_style()


def trapezoid(f, a, b, n):
    x = np.linspace(a, b, n + 1)
    y = f(x)
    return (b - a) / n * (y.sum() - 0.5 * (y[0] + y[-1]))


def simpson(f, a, b, n):
    x = np.linspace(a, b, n + 1)
    y = f(x)
    return (b - a) / n / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


ns = np.array([4, 6, 8, 12, 16, 24, 32, 48, 64, 96])      # panels (Simpson needs n even)
h = np.pi / ns                                             # step sizes
err_t = np.array([abs(trapezoid(np.sin, 0, np.pi, n) - 2) for n in ns])
err_s = np.array([abs(simpson(np.sin, 0, np.pi, n) - 2) for n in ns])
print("n = 4 and 96:  trapezoid", err_t[[0, -1]].round(6), " Simpson", err_s[[0, -1]])
```
```text
n = 4 and 96:  trapezoid [0.103881 0.000178]  Simpson [4.55975498e-03 1.27446356e-08]
```

The obvious first picture:

```python
fig, ax = plt.subplots()
ax.plot(h, err_t, "o-", label="trapezoid")
ax.plot(h, err_s, "s-", label="Simpson")
ax.set(xlabel="step size $h$", ylabel="absolute error")
ax.legend()
ph.savefig(fig, "s1_linear")
px = err_s / err_t.max() * 400                    # each Simpson error as a height in pixels on a 400-pixel-tall axis
print(f"{(px < 0.25).sum()} of Simpson's 10 points are less than a quarter pixel above zero; the largest is {px.max():.0f} pixels")
```
```text
7 of Simpson's 10 points are less than a quarter pixel above zero; the largest is 18 pixels
```

Open `figures/s1_linear.png`: the trapezoid curve bends and most of Simpson's points are on the floor. No reader can read an *order* off a bending curve. What works is a straight line whose slope is the order, with a reference triangle beside it. §4.2 draws it; getting there needs the object model (§3.1), scales and ticks (§3.4) and the habits of §3.9.

## 2. Intuition

A figure is a tree of objects, and you talk to the nodes directly.

```text
Figure                          the canvas: size in inches, dpi, background
 ├── Axes                       ONE plot: a rectangle with its own data coordinates
 │    ├── XAxis, YAxis          scale, tick locator, tick formatter, tick marks, label
 │    ├── 4 Spines              the frame lines
 │    └── Artists               Line2D (plot), PathCollection (scatter), Rectangle (bars),
 │                              AxesImage (imshow), Text, Legend, ...
 ├── Axes                       another plot (a colorbar is an Axes too)
 └── Text, Legend, ...          artists that belong to the figure itself
```

Everything you can see is an **Artist**. `ax.plot(...)` builds a `Line2D`, adds it to `ax.lines` and *returns* it; `fig, ax = plt.subplots()` hands you the two objects, and every call you need is a method of one of them. Four coordinate systems matter: **data**, **axes** (0 to 1 across one Axes), **figure** (0 to 1 across the canvas) and **display** (pixels). Each artist carries a transform from its coordinates to pixels; a log scale is the first stage of the data transform.

You can already guess the rule that fixes §1. On a log axis, distance is proportional to $\log_{10}$ of the value. A power law $e = Ch^p$ has $\log_{10} e = \log_{10} C + p\log_{10} h$: a straight line of slope $p$ when *both* axes are logarithmic; an exponential is straight with a log $y$ axis and a linear $x$ axis. **The scale that makes your model a straight line is the scale that lets a reader check it.** Every choice (range, scale, aspect, colour map, what an error bar means) is a claim the reader takes at face value; §3.9 makes that into rules.

## 3. Formal development

Seven **Build** sequences are marked; do them in order, and write each ▶ Try it first prediction down before running.

### 3.1 The object model, and why only the OO API

A **Figure** owns the canvas and a list of Axes. An **Axes** owns a patch, two **Axis** objects (`ax.xaxis`, `ax.yaxis`: each owns a scale, a *locator* that chooses tick positions, a *formatter* that writes tick text, tick marks and the label), four spines, and the artists you add. Artist properties (`color`, `linewidth`, `alpha`, `zorder`, ...) are read with `get_*` and written with `set_*` or `artist.set(...)`.

▶ **Try it first.** For `plt.subplots(2, 3, sharex=True)`, predict: the type and shape of `axs`; what `plt.subplots()` returns in place of `axs`; and `axs[1, 2].get_xlim()` after `axs[0, 0].set_xlim(0, 10)`.

```python
fig, axs = plt.subplots(2, 3, figsize=(6, 3.5), sharex=True)
print(type(fig).__name__, type(axs).__name__, axs.shape, type(axs[0, 0]).__name__, type(axs[0, 0].xaxis).__name__)
fig1, ax1 = plt.subplots()
print("one panel returns a bare", type(ax1).__name__, "| squeeze=False gives a 2-D array:", plt.subplots(2, 1, squeeze=False)[1].shape)
axs[0, 0].set_xlim(0, 10)
print("after set_xlim on [0, 0], panel [1, 2] reads", tuple(float(v) for v in axs[1, 2].get_xlim()),
      "| same x-axis group:", axs[0, 0].get_shared_x_axes().joined(axs[0, 0], axs[1, 2]))
plt.close("all")
```
```text
Figure ndarray (2, 3) Axes XAxis
one panel returns a bare Axes | squeeze=False gives a 2-D array: (2, 1)
after set_xlim on [0, 0], panel [1, 2] reads (0.0, 10.0) | same x-axis group: True
```

■ **What just happened.** `subplots` returns the Figure and a NumPy array of Axes (a bare Axes for one panel). `sharex=True` does not copy limits once; it *joins the x-axes into one group*, so a later change to one changes all.

Now the reason for the rule. pyplot keeps a hidden pointer to the "current" axes, the last one created, and `plt.plot` draws there. ▶ **Try it first.** Two figures are made, then one `plt.plot` call. Which axes gets the line?

```python
fig_a, ax_a = plt.subplots()
fig_b, ax_b = plt.subplots()
plt.plot([1, 2, 3])                      # THE TRAP: do not write this in a script
print("lines on ax_a:", len(ax_a.lines), "| lines on ax_b:", len(ax_b.lines))
plt.close("all")
```
```text
lines on ax_a: 0 | lines on ax_b: 1
```

■ **What just happened.** The line went to the figure created last, silently. With one figure that is a convenience; with two, or a function that makes its own, it is a bug (§8). The rule: **every drawing call is a method of an object you hold.** A function that takes an `ax` works on one panel, one cell of a grid, or an inset; Build 1 is written that way.

### 3.2 Coordinates, scales and normalisation

For a linear axis with limits $[x_{\min}, x_{\max}]$, a value $x$ sits a fraction $(x - x_{\min})/(x_{\max} - x_{\min})$ of the way along the axes. A log axis applies $\log_{10}$ first:

$$u(x) = \frac{\log_{10}x - \log_{10}x_{\min}}{\log_{10}x_{\max} - \log_{10}x_{\min}}.\tag{3.1}$$

A colour map is indexed by a number in $[0,1]$, and a *norm* produces it. `Normalize(vmin, vmax)` is linear and `LogNorm` is (3.1) again:

$$n(v) = \frac{v - v_{\min}}{v_{\max} - v_{\min}}, \qquad n_{\log}(v) = \frac{\ln v - \ln v_{\min}}{\ln v_{\max} - \ln v_{\min}}.\tag{3.2}$$

With $y \in [1, 1000]$ on a log axis, $y = 10$ should sit one third up and $y = 100$ two thirds:

```python
from matplotlib.colors import LogNorm, Normalize

fig, ax = plt.subplots()
ax.set(xlim=(0, 10), ylim=(1, 1000), yscale="log")
fig.canvas.draw()                                      # transforms are final only after a draw
to_axes = ax.transAxes.inverted().transform            # display pixels -> axes fraction
for y in (10, 100):
    print(f"data (5, {y}) is at axes fraction", np.round(to_axes(ax.transData.transform((5, y))), 4))
print("Normalize(0, 10)(2.5) =", float(Normalize(0, 10)(2.5)), "| LogNorm(1, 1000)(10) =", round(float(LogNorm(1, 1000)(10)), 4))
plt.close(fig)
```
```text
data (5, 10) is at axes fraction [0.5    0.3333]
data (5, 100) is at axes fraction [0.5    0.6667]
Normalize(0, 10)(2.5) = 0.25 | LogNorm(1, 1000)(10) = 0.3333
```

A log axis and a `LogNorm` are the same map, along a ruler and along a colour bar; both need positive data (§8).

### 3.3 The plot types, and the arguments that matter

| call | use it for | arguments that matter |
|---|---|---|
| `plot`, `step`, `stem` | ordered data $y(x)$; piecewise-constant; discrete sequences | `marker`, `ls`, `lw`, `color`, `label`; `where="mid"` |
| `scatter` | unordered points; a 3rd and 4th variable | `s` (area, pt²), `c` (array) + `cmap`/`norm`, `alpha`; slow beyond $10^5$ points |
| `bar`, `barh` | amounts or counts, with a **zero baseline** | `width`, `align`, `yerr`, `capsize` |
| `hist` | the distribution of a sample | `bins`, `range`, `density=True` |
| `errorbar` | an estimate with an uncertainty | `yerr`/`xerr`, `fmt`, `capsize` |
| `fill_between` | a band over $x$ | `alpha`, `where` |
| `imshow` | values on a **uniform grid** as an image | `origin`, `extent`, `aspect`, `interpolation`, `vmin`/`vmax` |
| `pcolormesh` | values on a grid with coordinates (non-uniform too) | `X, Y`, `shading` |
| `contour`, `contourf` | level sets of a field | `levels`, `colors`, `clabel` |
| `quiver`, `streamplot` | vector fields: arrows; field lines | `angles`, `pivot`; `density`, `color` |

Bars, histograms, error bars, images, contours and vector fields appear in §4, and `scatter` with `c` and `norm` in §3.7. The three left over:

```python
x = np.linspace(0, 4, 9)
fig, axs = plt.subplots(1, 3, figsize=(10, 3))
axs[0].plot(x, x**2, "o--", label="plot")
axs[0].step(x, x**2, where="mid", label="step")                     # piecewise constant; where= says where the jump sits
axs[0].legend()
axs[1].stem(x, np.cos(x))                                           # a discrete sequence
axs[2].fill_between(x, np.sin(x) - 0.3, np.sin(x) + 0.3, alpha=0.3)  # a band of half-width 0.3 around sin
for ax, name in zip(axs, ["plot, step", "stem", "fill_between"]):
    ax.set(title=name, xlabel="$x$ (arb.)", ylabel="$y$ (arb.)")
ph.savefig(fig, "gallery")
```

**`density=True` is a claim about *area*.** If $n_i$ of $N$ samples fall in a bin of width $w_i$, the bar height that makes the total area 1 is

$$f_i = \frac{n_i}{N\,w_i}, \qquad \sum_i f_i w_i = \frac{1}{N}\sum_i n_i = 1.\tag{3.3}$$

For equal bins the heights sum to $1/w$, not 1 (§4.3 prints both). "Density" is mass per unit of $x$: the y-label is "probability density (1/unit of $x$)", never "frequency", which is what lets you lay a PDF over a histogram (§4.3).

**`imshow` and its pixel rule.** `imshow(A)` draws `A` of shape `(ny, nx)` as colored rectangles. With `origin="lower"` and `extent=(x0, x1, y0, y1)`, element `A[i, j]` fills

$$x\in[x_0 + j\,\Delta x,\; x_0 + (j+1)\Delta x],\quad y\in[y_0 + i\,\Delta y,\; y_0 + (i+1)\Delta y],\qquad \Delta x = \tfrac{x_1-x_0}{n_x},\;\Delta y = \tfrac{y_1-y_0}{n_y}.\tag{3.4}$$

The default is `origin="upper"`: row 0 at the *top*, the convention of images and matrices, not of physics. And `extent` is the outer edge of the whole image, so for samples at pixel centres, $x_0 = x_{\text{first}} - \Delta x/2$. §4.5 tests both rules on a rendered figure. `pcolormesh(X, Y, A)` takes coordinates instead (cell centres with `shading="nearest"`, the `(ny+1, nx+1)` corners with `"flat"`) and handles non-uniform grids, which `imshow` cannot.

### 3.4 Scales, slopes and ticks: the convergence recipe (Builds 1 and 2)

If $e = C h^p$ then

$$\log_{10} e = \log_{10} C + p\,\log_{10} h,\tag{3.5}$$

a straight line of slope $p$ and intercept $\log_{10}C$ in log-log coordinates. The order of convergence is a slope: fit it with `np.polyfit(np.log10(h), np.log10(e), 1)`. Fit in log space; least squares in linear space is dominated by the largest errors, log space weights every decade equally.

**Build 2: the same data on three scales.** ▶ **Try it first.** Predict which scale gives the straightest line for each series. "Straightness" is the $R^2$ of a line fitted *in the plotted coordinates*.

```python
def r_squared(x, y):
    res = y - np.polyval(np.polyfit(x, y, 1), x)
    return 1 - res.var() / y.var()


scales = {"linear": (lambda v: v, lambda v: v), "semilogy": (lambda v: v, np.log10), "log-log": (np.log10, np.log10)}
print("R^2 of a straight line in the plotted coordinates")
for name, (fx, fy) in scales.items():
    print(f"  {name:9s} trapezoid {r_squared(fx(h), fy(err_t)):.4f}   Simpson {r_squared(fx(h), fy(err_s)):.4f}")

fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
for ax, title in zip(axs, scales):
    ax.plot(h, err_t, "o-", label="trapezoid")
    ax.plot(h, err_s, "s-", label="Simpson")
    ax.set(title=title, xlabel="step size $h$", ylabel="absolute error")
axs[1].set_yscale("log")
axs[2].set(xscale="log", yscale="log")
axs[0].legend()
ph.savefig(fig, "three_scales")

for name, e in [("trapezoid", err_t), ("Simpson", err_s)]:
    p, b, se = ph.loglog_fit(h, e)
    print(f"{name:9s} fitted order {p:.3f}; the 'standard error' from scatter alone is {se:.3f}")
```
```text
R^2 of a straight line in the plotted coordinates
  linear    trapezoid 0.9241   Simpson 0.7301
  semilogy  trapezoid 0.8448   Simpson 0.8466
  log-log   trapezoid 1.0000   Simpson 1.0000
trapezoid fitted order 2.002; the 'standard error' from scatter alone is 0.001
Simpson   fitted order 4.018; the 'standard error' from scatter alone is 0.005
```

■ **What just happened.** Only log-log makes both series straight ($R^2 > 0.9999$); semilog bends both, and on linear axes Simpson is the *worst* line because it hugs zero. One caution: the fitted orders are 2.002 and 4.018, and 4.018 is several "standard errors" from 4. But that standard error measures scatter about a line, and these points do not scatter; they lie on a slightly curved line because the next term of the error expansion (order $h^6$) still contributes. Report "4.0 ± 0.1 for $h \in [0.03, 0.8]$", range included.

**Build 1: from `plt.plot` to `plot_convergence(ax, h, err, order)`.** (1) The §1 picture. (2) A function of `ax` using `ax.loglog` that returns the fitted slope. (3) The fitted line $10^b h^p$. (4) A **slope triangle** for the *expected* order, corner $(x_0, y_0)$, width $d$ decades:

$$(x_0,\,y_0),\qquad (x_0 10^{d},\,y_0),\qquad (x_0 10^{d},\,y_0 10^{pd}),\tag{3.6}$$

whose hypotenuse has $\log_{10}$-slope $pd/d = p$; label the legs "1" and "$p$". (5) Raise `ValueError` for $h \le 0$ or $e \le 0$, which a log axis would drop silently. Write these yourself; here is step 4, checked against the helper of §3.7:

```python
def my_guide(ax, order, x0, y0, decades=0.5):
    x1, y1 = x0 * 10**decades, y0 * 10 ** (order * decades)               # (3.6)
    ax.plot([x0, x1, x1, x0], [y0, y0, y1, y0], color="0.3", lw=0.9)        # the labels "1" and p are two ax.annotate calls


fig, ax = plt.subplots()
ax.loglog(h, err_t, "o")
my_guide(ax, 2, 0.15, 3e-4)
ph.add_slope_guide(ax, 2, 0.15, 3e-4)
print("my triangle equals the helper's:", np.allclose(ax.lines[1].get_xydata(), ax.lines[2].get_xydata()))
plt.close(fig)
```
```text
my triangle equals the helper's: True
```

**Ticks.** An Axis has a *locator* (where ticks go) and a *formatter* (what they say): `MultipleLocator(base)` ticks every `base`; `FuncFormatter(f)` calls your `f(value, position)` per label; `NullFormatter` silences them (on a short log axis it stops minor labels colliding). For trigonometric axes, multiples of $\pi/2$ read better than 1.57:

```python
fig, ax = plt.subplots()
xx = np.linspace(-2 * np.pi, 2 * np.pi, 400)
ax.plot(xx, np.sin(xx))
ph.pi_ticks(ax, "x", np.pi / 2)                       # MultipleLocator(pi/2) + a FuncFormatter built on fractions.Fraction
ax.set(xlabel="angle $\\theta$ (rad)", ylabel="$\\sin\\theta$ (dimensionless)")
fmt = ax.xaxis.get_major_formatter()
print([fmt(v) for v in (-np.pi, np.pi / 2, 3 * np.pi / 2)])
ph.savefig(fig, "pi_ticks")
```
```text
['$-\\pi$', '$\\frac{\\pi}{2}$', '$\\frac{3\\pi}{2}$']
```

**Zeros and negatives on a log axis vanish.** ▶ **Try it first.** Five points, two $\le 0$, on a log $y$-axis: does Matplotlib warn, and where do they go?

```python
y = np.array([1e-1, 1e-2, 0.0, -1e-3, 1e-4])
fig, ax = plt.subplots()
(ln,) = ax.semilogy(np.arange(5), y, "o-")
fig.canvas.draw()
print("display y of the five points (pixels):", np.round(ax.transData.transform(ln.get_xydata())[:, 1]).astype(int))
print(f"canvas is {fig.bbox.height:.0f} px tall; autoscaled y-limits ignore the two bad points:", tuple(float(f"{v:.3g}") for v in ax.get_ylim()))
plt.close(fig)
```
```text
display y of the five points (pixels): [    379     269 -110158 -110158      47]
canvas is 400 px tall; autoscaled y-limits ignore the two bad points: (7.08e-05, 0.141)
```

■ **What just happened.** No warning. The two points map to about $-1.1\times10^{5}$ pixels, far below the canvas, and the line plunges off the bottom edge; a reader sees a curve that stops. Remedies: plot what you meant (`np.abs`, noted), use `symlog` (linear within $\pm$`linthresh`, log outside) if the sign matters, or fix the data. `plot_convergence` raises instead.

### 3.5 Limits, aspect, twin axes, text

- **Limits and aspect.** `ax.set_ylim(100, 0)` (larger first) *inverts* the axis, right for depth; setting limits turns autoscaling off. `ax.set_aspect("equal")` makes one data unit the same length on both axes, as geometry needs; it changes the box, not the data. `ax.set_box_aspect(1)` makes a square panel whatever the data.
- **Twin axes.** `ax.twinx()` adds a second $y$ scale over the same $x$; nothing ties the scales together, so the curves cross where you choose (§3.9).
- **Legends.** `loc`, `ncol`, `bbox_to_anchor` move it off the data; with constrained layout a figure legend can sit outside: `fig.legend(loc="outside lower center")`. Better: label the line (`ph.label_end`).
- **Text.** `ax.annotate(text, xy, xytext, arrowprops=..., textcoords="offset points")` points at a data location (§4.1 uses it with an artist as the coordinate system); `ax.text(x, y, s, transform=ax.transAxes)` puts text at a fixed place *in the panel* (0 to 1), surviving rescaling. Labels use **mathtext**, a built-in LaTeX subset needing no LaTeX install: `r"$\alpha$ (rad)"`, `r"$\mathrm{m\,s^{-1}}$"`; the labels of every figure here use it. Unsupported commands (`\begin{align}`) fail at *draw* time, not when the string is set.

### 3.6 Colour (Build 4)

Pick the map by the data. **Sequential** (ordered magnitudes): `viridis`, `magma`, `cividis`. **Diverging** (signed, meaningful zero): `RdBu_r`, `coolwarm`, centred at zero (`Normalize(-m, m)` or `TwoSlopeNorm(vcenter=0)`). **Cyclic** (phase, angle): `twilight`. **Categorical** (unordered classes): at most about eight hues, such as Okabe and Ito's palette.

A sequential map is honest when equal data steps look like equal steps; a necessary condition is that perceived *lightness* rises monotonically along it. CIE 1976 lightness $L^*$ (0 black, 100 white) of an sRGB colour is

$$L^* = 116\,f(Y) - 16,\quad f(t)=\begin{cases}t^{1/3}, & t > 216/24389\\ (24389\,t/27 + 16)/116, & \text{otherwise,}\end{cases}\quad Y = 0.2126\,R_\ell + 0.7152\,G_\ell + 0.0722\,B_\ell,\tag{3.7}$$

where $C_\ell = C/12.92$ for $C \le 0.04045$ and $((C + 0.055)/1.055)^{2.4}$ otherwise (undoing the sRGB gamma). Monotone $L^*$ does not guarantee uniformity, but a map where it fails is not uniform.

▶ **Try it first.** Predict for `viridis` and `jet`: is $L^*$ monotone? If not, how many ridges (peaks rising at least 2 $L^*$ units above their surroundings), and is the top end of `jet` dark or light?

```python
from scipy.signal import find_peaks


def lstar(rgb):
    """CIE L* of sRGB colours: rows of `rgb`, floats in [0, 1], eq. (3.7)."""
    rgb = np.atleast_2d(rgb)[:, :3]
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    Y = lin @ np.array([0.2126, 0.7152, 0.0722])
    f = np.where(Y > 216 / 24389, np.cbrt(Y), (24389 / 27 * Y + 16) / 116)
    return 116 * f - 16


def lightness(name):
    return lstar(plt.get_cmap(name)(np.linspace(0, 1, 256)))


def ridges(a, prominence=2.0):
    """Peaks rising at least `prominence` above their surroundings (flat tops are handled)."""
    return len(find_peaks(a, prominence=prominence)[0])


print("map        L* at 0, 1/4, 1/2, 3/4, 1   monotone  ridges")
for name in ["viridis", "magma", "jet", "hsv"]:
    L = lightness(name)
    print(f"{name:9s}  {np.round(L[[0, 64, 128, 192, 255]]).astype(int)!s:24s} {np.all(np.diff(L) > 0)!s:8s}  {ridges(L)}")
```
```text
map        L* at 0, 1/4, 1/2, 3/4, 1   monotone  ridges
viridis    [15 36 54 73 91]         True      0
magma      [ 0 22 44 69 98]         True      0
jet        [13 55 90 71 25]         False     1
hsv        [53 90 91 39 53]         False     3
```

■ **What just happened.** `viridis` and `magma` rise monotonically. `jet` climbs to $L^*\approx96$ (yellow) then *falls* to about 25 at dark red: its top is darker than its middle. (The 2-unit threshold is this module's convention, not a perceptual measurement.) **Build 4** shows what that does to a picture: one smooth bump $f = e^{-r^2/2}$ with a single maximum, drawn with both maps, and the $L^*$ the eye receives along $y=0$:

```python
g = np.linspace(-3, 3, 241)
GX, GY = np.meshgrid(g, g)
bump = np.exp(-(GX**2 + GY**2) / 2)
profile = bump[120]                                       # row 120 is y = 0
fig, axs = plt.subplots(2, 2, figsize=(9, 7))
for ax, name in zip(axs[0], ["jet", "viridis"]):
    im = ax.imshow(bump, origin="lower", extent=(-3.025, 3.025, -3.025, 3.025), cmap=name, vmin=0, vmax=1)
    ax.grid(False)
    ax.set(title=name, xlabel="$x$ (units of $\\sigma$)", ylabel="$y$ (units of $\\sigma$)")
    fig.colorbar(im, ax=ax, label="$f$")
axs[1, 0].plot(g, profile, "k")
axs[1, 0].set(title="the data along $y=0$", xlabel="$x$ (units of $\\sigma$)", ylabel="$f(x, 0)$")
seen = {}
for name in ("jet", "viridis"):
    seen[name] = np.interp(profile, np.linspace(0, 1, 256), lightness(name))   # L* of the colour drawn at each x
    axs[1, 1].plot(g, seen[name], label=name)
axs[1, 1].set(title="perceived lightness along $y=0$", xlabel="$x$ (units of $\\sigma$)", ylabel="CIE $L^*$ (0 to 100)")
axs[1, 1].legend()
ph.savefig(fig, "jet_vs_viridis")
print(f"peaks along y = 0: data {ridges(profile, 0.01)}, viridis {ridges(seen['viridis'])}, jet {ridges(seen['jet'])}")
print(f"jet L* at the centre {seen['jet'][120]:.0f}, on its two bright rings {seen['jet'].max():.0f}; viridis centre {seen['viridis'][120]:.0f}")
```
```text
peaks along y = 0: data 1, viridis 1, jet 2
jet L* at the centre 25, on its two bright rings 96; viridis centre 91
```

The data have one peak and `viridis` shows one. `jet` shows a dark centre between two bright rings that exist in the colour map, not in the data. Crameri, Shephard and Heron (2020) make this case at length against rainbow maps; the numbers above are a small piece of it. Course rule for `jet`: never, unless the figure is *about* `jet`.

**Colour-blind safety.** Red–green colour-vision deficiency is common (on the order of one man in twelve), and many figures are read in greyscale. Okabe and Ito's eight colours were chosen to stay distinguishable under the common deficiencies; this module's cycle is that palette, first four colours most separated. They are separated mainly by hue; some pairs have almost the same lightness (printed below), so a figure that must survive greyscale also needs distinct markers or line styles.

```python
from matplotlib.colors import to_rgb

print("L* of the palette:", " ".join(f"{c} {lstar(to_rgb(c))[0]:.0f}" for c in ph.OKABE_ITO))
```
```text
L* of the palette: #0072B2 46 #D55E00 54 #009E73 58 #CC79A7 61 #E69F00 71 #56B4E9 70 #F0E442 89 #000000 0
```

### 3.7 Layout, styles, and the house helpers (Builds 5 and 6)

- **Constrained layout.** `plt.subplots(..., layout="constrained")` (also `constrained_layout=True`, or the rc key `figure.constrained_layout.use`) solves for margins at draw time so labels, titles and colorbars do not collide. Leave it on; the older `tight_layout()` fights with colorbars and `suptitle` (§8). Hand-placed axes (`fig.add_axes`) need `layout="none"`.
- **Unequal panels and insets.** `fig.add_gridspec(2, 3, width_ratios=[2, 1, 1])`, or `plt.subplot_mosaic([["a", "a", "b"], ["c", "d", "b"]])` (a dict of Axes). `ax.inset_axes((x0, y0, w, h))` takes axes fractions; `ax.indicate_inset_zoom(inset)` draws the zoom box. `fig.suptitle` titles the whole figure (a caption is usually better).
- **Small multiples** repeat the same axes with one thing varied: `sharex=True, sharey=True`, and for images one shared `Normalize` and one `fig.colorbar(im, ax=axs)`. Letting each panel autoscale its own colours is the commonest way they lie.
- **Styles.** `plt.rcParams` holds defaults, `plt.rc_context({...})` scopes a change to a `with` block, `plt.style.use("name")` loads a bundle. The **house style**, `ph.HOUSE_STYLE`: constrained layout, `savefig.dpi` 150, no top/right spines, a faint grid, the Okabe–Ito cycle (`axes.prop_cycle` set to `cycler(color=[...])`), `image.interpolation: nearest` (never blur a computed grid), no legend frame, viridis. One trap: the grid is drawn *over* images, so image panels need `ax.grid(False)`.

**Build 5: `plot_helpers.py`.** Write a module with (a) `HOUSE_STYLE` and `use_house_style()`, (b) `savefig(fig, name, ...)` that creates `figures/`, trims with `bbox_inches="tight"`, closes the figure and returns the path, (c) `add_slope_guide(ax, order)` from (3.6), (d) `plot_convergence`. Then compare with the file beside these notes, which adds `pi_ticks`, `label_end`, `loglog_fit` and `rgba_at` (§4.5). Later modules may import it by copying it next to their own project; `project.py` exercises every function in it.

**Build 6: small multiples.** Six damped oscillators, $\ddot x + 2\zeta\dot x + x = 0$, $x(0)=1$, $\dot x(0)=0$ (time in units of $1/\omega_0$): one phase plane per damping ratio $\zeta$, shared axes, points coloured by time, one colour scale, one colorbar. ▶ **Try it first.** For which of $\zeta = 0.05, 0.2, 0.5, 1, 1.5, 3$ does $x$ change sign?

```python
from scipy.integrate import solve_ivp

zetas = [0.05, 0.2, 0.5, 1.0, 1.5, 3.0]
tt = np.linspace(0, 30, 600)
tnorm = Normalize(0, 30)                                   # ONE normalisation for all six panels
fig, axs = plt.subplots(2, 3, figsize=(10, 6), sharex=True, sharey=True)
crosses = []
for ax, z in zip(axs.flat, zetas):
    sol = solve_ivp(lambda t, s: [s[1], -2 * z * s[1] - s[0]], (0, 30), [1.0, 0.0], t_eval=tt, rtol=1e-10, atol=1e-12)
    x_, v_ = sol.y
    crosses.append(bool(np.any(x_ < 0)))
    sc = ax.scatter(x_, v_, c=tt, cmap="viridis", norm=tnorm, s=6)
    ax.set_title(f"$\\zeta = {z:g}$")
    ax.set_box_aspect(1)
for ax in axs[1]:
    ax.set_xlabel("$x / x(0)$ (dimensionless)")
for ax in axs[:, 0]:
    ax.set_ylabel("$\\dot x / (\\omega_0 x(0))$ (dimensionless)")
cb = fig.colorbar(sc, ax=axs, shrink=0.8, label="time $\\omega_0 t$ (rad)")
ph.savefig(fig, "small_multiples")
print("x changes sign:", dict(zip(zetas, crosses)), "| colorbar range:", (float(cb.vmin), float(cb.vmax)))
```
```text
x changes sign: {0.05: True, 0.2: True, 0.5: True, 1.0: False, 1.5: False, 3.0: False} | colorbar range: (0.0, 30.0)
```

■ **What just happened.** Only $\zeta < 1$ crosses $x=0$: critical and over-damped oscillators released from rest creep to zero. The figure says it at a glance because the panels share limits and the colour means the same everywhere. `project.py` also checks each panel against the closed-form $x(t)$.

### 3.8 Saving, 3D, animation, and reading data (Build 7)

**Saving.** PNG is raster: pixel size is `figsize` $\times$ `dpi`. PDF and SVG are vector, so text stays text and lines stay curves, but the file grows with the number of artists; `rasterized=True` turns one dense artist into an image inside the vector file. `bbox_inches="tight"` changes the pixel size (`tight=False` keeps it); `transparent=True` drops the white background; `rcParams["pdf.fonttype"] = 42` embeds editable fonts, which journals often require. **Design at final size:** font sizes are in points, so a 6.4-inch figure shrunk to a 3.4-inch column shrinks its text too.

```python
import matplotlib.image as mimg

big = np.random.default_rng(2).normal(size=(20_000, 2))
fig, ax = plt.subplots(figsize=(4, 3))
dots = ax.scatter(*big.T, s=2)
ax.set(xlabel="$x$ (arb.)", ylabel="$y$ (arb.)")
kb = {fmt: ph.savefig(fig, f"scatter20k_{fmt}", fmt=fmt, tight=False, close=False).stat().st_size / 1e3 for fmt in ("png", "pdf", "svg")}
dots.set_rasterized(True)
kb["pdf, rasterized"] = ph.savefig(fig, "scatter20k_raster", fmt="pdf", tight=False, close=False).stat().st_size / 1e3
print("png shape at dpi 150, 4x3 in:", mimg.imread(ph.FIGDIR / "scatter20k_png.png").shape)
print("file sizes (kB), 20,000 points:", {k: round(v) for k, v in kb.items()})
plt.close(fig)
```
```text
png shape at dpi 150, 4x3 in: (450, 600, 4)
file sizes (kB), 20,000 points: {'png': 52, 'pdf': 313, 'svg': 2136, 'pdf, rasterized': 28}
```

**3D in one paragraph.** `fig.add_subplot(projection="3d")` gives an `Axes3D`; `ax.plot_surface(X, Y, Z)` and `ax.view_init(elev, azim)` draw and aim it. It suits *exploring* a shape and figures whose message is the shape (a saddle), but depth is ambiguous and the surface hides itself, so it is a poor way to show *values*: use a 2D colour map with a colorbar (§4.5) or line plots of slices.

**Build 7: a pendulum animation, saved headless.** `FuncAnimation(fig, update, frames=N)` calls `update(i)` per frame; `update` should change existing artists (`set_data`), not create new ones. Saving needs a *writer*: `PillowWriter` makes a GIF from Pillow, which Matplotlib already requires. Test `update` as data code *before* animating: call it for every frame and assert on the artists' data. Then write the GIF (fallback: a strip of frames).

```python
import math

from matplotlib import animation

th0, n_frames, t_end = 2.5, 60, 10.0
sol = solve_ivp(lambda t, s: [s[1], -np.sin(s[0])], (0, t_end), [th0, 0.0], dense_output=True, method="DOP853", rtol=1e-11, atol=1e-12)
t_frames = np.linspace(0, t_end, n_frames, endpoint=False)
theta_f = sol.sol(t_frames)[0]

fig, ax = plt.subplots(figsize=(3.2, 3.2))
ax.set(xlim=(-1.2, 1.2), ylim=(-1.2, 0.3), xlabel="$x$ (units of $L$)", ylabel="$y$ (units of $L$)")
ax.set_aspect("equal")
(rod,) = ax.plot([], [], "-o", color="0.2", ms=9, mfc=ph.OKABE_ITO[1])


def update(i):
    """Move the bob to frame i; return the changed artists."""
    rod.set_data([0, math.sin(theta_f[i])], [0, -math.cos(theta_f[i])])
    return (rod,)


for i in range(n_frames):                                  # data checks, not pixels
    update(i)
    xs, ys = rod.get_data()
    assert math.isclose(math.hypot(xs[1], ys[1]), 1.0, rel_tol=1e-12)               # rod length stays 1
    assert math.isclose(math.atan2(xs[1], -ys[1]), theta_f[i], abs_tol=1e-12)       # bob is at angle theta_i
gif = ph.FIGDIR / "pendulum_demo.gif"
if animation.writers.is_available("pillow"):
    anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=50, blit=False)
    anim.save(gif, writer=animation.PillowWriter(fps=20), dpi=60)
    print(f"{gif.name}: header {gif.read_bytes()[:6]!r}, {gif.stat().st_size / 1e3:.0f} kB, {n_frames} frames")
else:                                                       # fallback: save a strip of frames instead
    print("no Pillow writer")
plt.close(fig)
```
```text
pendulum_demo.gif: header b'GIF89a', 66 kB, 60 frames
```

**Reading data.** `np.loadtxt` for a numeric table, the `csv` module (PY-2) for named or mixed columns, `json` for nested records. Missing values need `np.genfromtxt` or pandas (PY-9). `project.py` Part 2 round-trips the convergence table through CSV and JSON exactly.

### 3.9 What makes a figure honest

Rules of thumb from practice, stated as guidance, not as measured laws.

1. **Position may be truncated; length may not.** A bar's *length* is its value, so its baseline is zero. Tufte's **lie factor** puts a number on a violation:

$$\text{lie factor} = \frac{\text{size of the effect shown in the graphic}}{\text{size of the effect in the data}},\qquad \text{effect} = \frac{\text{second} - \text{first}}{\text{first}}.\tag{3.8}$$

   A dot or line is read by *position*, so a visible truncated range can be fine when the claim concerns differences. §4.6 computes one.
2. **Aspect ratio sets the apparent slope** (measured below). Pick a shape in which the slopes that carry the message are neither flat nor vertical, and put the numbers on the figure.
3. **Dual axes are two arguments sharing a frame.** With `twinx` you choose both limits, so you choose where the curves cross. If both series matter, stack two panels with a shared $x$-axis.
4. **An error bar must say what it is.** Standard deviation (spread of the data), standard error (precision of a mean) and a confidence interval are different numbers:

$$\mathrm{SE} = \frac{s}{\sqrt n},\qquad \text{95\% CI half-width} = t_{0.975,\,n-1}\,\mathrm{SE},\tag{3.9}$$

   with $s$ the sample SD (`ddof=1`); the CI is for the mean of roughly normal data (F4, ST-2). Say which in the caption or legend.
5. **Overplotting hides density:** $10^5$ scatter points are a blob; use `alpha`, `hexbin`, a 2D histogram, or rasterise.
6. **One message per figure.** A figure that needs a paragraph to say where to look has two messages: split it.
7. **Label directly.** A legend sends the reader back and forth; label the line's end (`ph.label_end`) unless labels collide.

```python
from scipy import stats

rng = np.random.default_rng(3)
print("   n   mean    SD    SE   95% CI half-width")
for n in (5, 50, 500):
    s = rng.normal(10, 2, n)
    sd = s.std(ddof=1)
    se = sd / np.sqrt(n)                                    # (3.9)
    half = stats.t.ppf(0.975, n - 1) * se
    print(f"{n:4d}  {s.mean():5.2f}  {sd:5.2f}  {se:5.2f}   {half:5.2f}")

for w, hgt in [(4, 4), (8, 2), (2, 8)]:                     # axes of w x h inches, both data ranges 0 to 10
    fig = plt.figure(figsize=(w + 1, hgt + 1), layout="none")
    ax = fig.add_axes((0.5 / (w + 1), 0.5 / (hgt + 1), w / (w + 1), hgt / (hgt + 1)))
    ax.plot([0, 10], [0, 10])
    ax.set(xlim=(0, 10), ylim=(0, 10))
    (x0, y0), (x1, y1) = ax.transData.transform([(0, 0), (10, 10)])
    print(f"axes {w} x {hgt} in: the line y = x is drawn at {np.degrees(np.arctan2(y1 - y0, x1 - x0)):.1f} degrees")
    plt.close(fig)
```
```text
   n   mean    SD    SE   95% CI half-width
   5   9.55   3.34   1.49    4.15
  50  10.09   2.02   0.29    0.57
 500  10.14   1.95   0.09    0.17
axes 4 x 4 in: the line y = x is drawn at 45.0 degrees
axes 8 x 2 in: the line y = x is drawn at 14.0 degrees
axes 2 x 8 in: the line y = x is drawn at 76.0 degrees
```

The SD settles near the true value 2 as $n$ grows (the 3.34 at $n=5$ is the noise of a five-point sample) because it describes the data; the SE and CI shrink like $1/\sqrt n$ because they describe the mean. At $n=500$ an error bar drawn as SD is about 20 times longer than one drawn as SE. And one line, three axes shapes: 45.0, 14.0 and 76.0 degrees.

## 4. Worked examples

Seven complete programs, increasing in difficulty; each also appears, fuller, in `project.py`. Every number quoted is printed in one of them.

### Example 4.1: Anatomy of a figure, labelled by itself

Every Artist knows its own bounding box, and `annotate` accepts an *artist* as its coordinate system (`xycoords=artist`: `xy` is then a fraction of that artist's box). This figure marks seven artists with numbers and keys them beside the axes.

```python
fig, ax = plt.subplots(figsize=(8.5, 4.5))
t = np.linspace(0, 2 * np.pi, 100)
(curve,) = ax.plot(t, np.sin(t), label="$x = \\sin t$")
ax.set(xlabel="time $t$ (rad)", ylabel="displacement $x$ (m)", title="title")
ax.legend(loc="lower left")
fig.canvas.draw()                                            # artists report their boxes only after a draw
named = {"Axes: one plot, with its own patch": (ax.patch, (0.95, 0.9)), "Line2D: the curve": (curve, (0.38, 1.0)),
         "Legend": (ax.get_legend(), (1.0, 0.5)), "title (a Text)": (ax.title, (1.0, 0.5)),
         "xlabel (a Text)": (ax.xaxis.label, (1.0, 0.5)), "tick label (a Text)": (ax.get_xticklabels()[3], (1.0, 0.5)),
         "Spine: a frame line": (ax.spines["bottom"], (0.9, 0.5))}
circle = {"boxstyle": "circle,pad=0.25", "fc": "0.2", "ec": "none"}
for k, (text, (artist, frac)) in enumerate(named.items(), 1):
    ax.annotate(str(k), xy=frac, xycoords=artist, xytext=(12, 0), textcoords="offset points", ha="center", va="center",
                color="white", fontsize=8, bbox=circle, annotation_clip=False)
    ax.text(1.03, 1.0 - 0.1 * k, f"{k}  {text}", transform=ax.transAxes, fontsize=8, va="center")
ph.savefig(fig, "anatomy_notes")

from collections import Counter

for name, parent in [("Figure", fig), ("Axes", ax), ("XAxis", ax.xaxis)]:
    print(f"{name:6s} children:", dict(Counter(type(c).__name__ for c in parent.get_children())))
```
```text
Figure children: {'Rectangle': 1, 'Axes': 1}
Axes   children: {'Line2D': 1, 'Annotation': 7, 'Text': 10, 'Spine': 4, 'XAxis': 1, 'YAxis': 1, 'Legend': 1, 'Rectangle': 1}
XAxis  children: {'Text': 2, 'XTick': 9}
```

The tree of §2 is real: the Axes holds one `Line2D`, four `Spine`s, the two Axis objects and some texts; each Axis holds one tick object per tick, owning its mark, grid line and labels. When tick text is wrong, the fix is a method of the Axis (`ax.xaxis.set_major_formatter`), not a search through the figure. `project.py` Part 1 marks ten artists (adding the ylabel, an Axis and the Figure) and asserts their types.

### Example 4.2: The convergence plot, with slopes 2 and 4 annotated

The claim: the trapezoid error is $\approx h^2/6$ and Simpson's is $\approx h^4/90$ for $\int_0^\pi\sin$. PY-5 §3.10 (its equation (3.2), a Taylor argument) gives $T - I \approx \tfrac{(b-a)h^2}{12}f''(\xi)$ and $S - I \approx \tfrac{(b-a)h^4}{180}f^{(4)}(\xi)$. That $(b-a)f''(\xi)$ is a sum of per-panel terms $\sum_i h\,f''(\xi_i)$, a Riemann sum for $\int_a^b f'' = f'(b) - f'(a)$; the same for $f^{(4)}$. So the leading terms are $T - I \approx \tfrac{h^2}{12}[f'(b) - f'(a)]$ and $S - I \approx \tfrac{h^4}{180}[f'''(b) - f'''(a)]$, the first terms of the Euler–Maclaurin expansion. SymPy evaluates them; the figure checks them.

```python
import sympy as sp

x_, h_ = sp.symbols("x h", positive=True)
f = sp.sin(x_)
lead_T = h_**2 / 12 * (sp.diff(f, x_).subs(x_, sp.pi) - sp.diff(f, x_).subs(x_, 0))
lead_S = h_**4 / 180 * (sp.diff(f, x_, 3).subs(x_, sp.pi) - sp.diff(f, x_, 3).subs(x_, 0))
print("leading terms of T - I and S - I:", sp.simplify(lead_T), "and", sp.simplify(lead_S))

fig, ax = plt.subplots()
p_t = ph.plot_convergence(ax, h, err_t, 2, label="trapezoid")
p_s = ph.plot_convergence(ax, h, err_s, 4, label="Simpson")
ax.legend(loc="lower right")
ph.savefig(fig, "convergence")
print(f"fitted orders: trapezoid {p_t:.3f}, Simpson {p_s:.3f}")
print(f"err / (h^2/6) at n = 96: {err_t[-1] / (h[-1] ** 2 / 6):.5f}; err / (h^4/90) at n = 96: {err_s[-1] / (h[-1] ** 4 / 90):.5f}")
```
```text
leading terms of T - I and S - I: -h**2/6 and h**4/90
fitted orders: trapezoid 2.002, Simpson 4.018
err / (h^2/6) at n = 96: 1.00002; err / (h^4/90) at n = 96: 1.00013
```

Open `figures/convergence.png`: two straight lines with fitted slopes 2.00 and 4.02, each with a triangle of its *expected* slope half a decade below, so the eye compares parallel edges. At $n=96$ the errors are within $2\times10^{-5}$ and $1.3\times10^{-4}$ (relative) of $h^2/6$ and $h^4/90$. Markers are measurements, the thin line is what the data say (the *fit*), the triangle is what theory says (the *reference*); a skeptic can compare the last two without trusting the first. The figure does not show that the orders hold below $h=0.03$ or for non-smooth integrands: state the range.

### Example 4.3: A histogram, a fitted Gaussian and a KDE (F4)

1,500 waiting times from a skewed Gamma(4, scale 1.5). The fitted Gaussian has the sample's mean and SD; a **kernel density estimate** puts a Gaussian bump of width $b$ on each point,

$$\hat f(x) = \frac{1}{nb}\sum_{i=1}^n \varphi\!\Big(\frac{x - x_i}{b}\Big),\qquad b = n^{-1/5}\, s\ \ \text{(Scott's rule in 1D)},\tag{4.1}$$

with $\varphi$ the standard normal density and $s$ the sample SD; SciPy's `gaussian_kde` calls $n^{-1/5}$ its `factor`.

```python
from scipy import stats

true = stats.gamma(4, scale=1.5)
sample = true.rvs(1500, random_state=np.random.default_rng(0))
mu, sd = stats.norm.fit(sample)                                       # maximum likelihood: sample mean and ddof=0 SD
kde = stats.gaussian_kde(sample)
xs = np.linspace(-6, 24, 1500)
fig, ax = plt.subplots()
counts, edges, _ = ax.hist(sample, bins=30, range=(0, 20), density=True, color="0.8", edgecolor="white")
ax.plot(xs, stats.norm.pdf(xs, mu, sd), label=f"Gaussian fit ($\\mu$={mu:.2f}, $\\sigma$={sd:.2f})")
ax.plot(xs, kde(xs), label="Gaussian KDE (Scott)")
ax.plot(xs, true.pdf(xs), "k:", label="true Gamma(4, 1.5) density")
ax.set(xlim=(-3, 22), xlabel="waiting time (s)", ylabel="probability density (1/s)")
ax.legend()
ph.savefig(fig, "hist_fit_kde")

tv = lambda f: 0.5 * np.trapezoid(np.abs(f - true.pdf(xs)), xs)         # total-variation distance to the truth
print(f"histogram area {np.sum(counts * np.diff(edges)):.6f}, sum of heights {counts.sum():.2f} = 1/w = {1 / np.diff(edges)[0]:.2f}; Scott factor {kde.factor:.4f} = n^(-1/5) = {len(sample) ** -0.2:.4f}; kernel sd b = {kde.factor * sample.std(ddof=1):.3f} s")
print(f"total-variation distance to the true density: Gaussian fit {tv(stats.norm.pdf(xs, mu, sd)):.3f}, KDE {tv(kde(xs)):.3f}")
```
```text
histogram area 1.000000, sum of heights 1.50 = 1/w = 1.50; Scott factor 0.2316 = n^(-1/5) = 0.2316; kernel sd b = 0.690 s
total-variation distance to the true density: Gaussian fit 0.132, KDE 0.019
```

The message is one sentence: *a Gaussian with the right mean and SD is the wrong shape for these data* (it puts mass below zero for a time that cannot be negative, and misses the skew). The numbers back it: 0.132 from the truth in total variation for the Gaussian, 0.019 for the KDE. `density=True` made the histogram area 1 (3.3), so the curves share its axis. The KDE's width $b$ is a modelling choice: small gives spikes, large over-smooths (Exercise B4).

### Example 4.4: Monte Carlo error bars with a $1/\sqrt n$ guide (F5)

Estimate $\pi$ by throwing $n$ darts at the unit square and counting the fraction inside the quarter disk: $\hat\pi = 4\hat p$, $p = \pi/4$. Each dart contributes $4\cdot\mathbf 1[\text{inside}]$, of variance $16p(1-p) = \pi(4-\pi)$, so over repeated runs

$$\operatorname{SD}(\hat\pi_n) = \sqrt{\pi(4-\pi)/n} \approx 1.64/\sqrt n.\tag{4.2}$$

▶ **Try it first.** What is the log-log slope of SD against $n$, and by what factor must $n$ grow to halve the error bar?

```python
rng5 = np.random.default_rng(5)
R = 400                                                     # replicate runs per n
n_list = np.unique(np.logspace(2, 4.5, 8).astype(int))


def pi_hat(n):
    """4 * (fraction of n uniform darts that land inside the quarter disk)."""
    xy = rng5.random((n, 2))
    return 4.0 * np.mean(np.sum(xy**2, axis=1) < 1)


est = np.array([[pi_hat(n) for _ in range(R)] for n in n_list])      # shape (len(n_list), R)
mean_est, sd_est = est.mean(axis=1), est.std(axis=1, ddof=1)
sd_theory = np.sqrt(np.pi * (4 - np.pi) / n_list)                    # (4.2)
slope = ph.loglog_fit(n_list, sd_est)[0]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
ax1.errorbar(n_list, mean_est, yerr=sd_est, fmt="o-", capsize=3, label="mean $\\pm$ 1 SD of 400 runs")
ax1.axhline(np.pi, color="k", lw=0.8, ls=":")
ax1.set(xscale="log", xlabel="darts per run $n$", ylabel="estimate of $\\pi$ (dimensionless)")
ax1.legend()
ax2.errorbar(n_list, sd_est, yerr=sd_est / np.sqrt(2 * (R - 1)), fmt="o", capsize=3, label="SD across runs ($\\pm$ 1 SE)")
ax2.plot(n_list, sd_theory, "-", lw=0.8, label="$\\sqrt{\\pi(4-\\pi)/n}$")
ax2.set(xscale="log", yscale="log", xlabel="darts per run $n$", ylabel="standard deviation of the estimate")
ph.add_slope_guide(ax2, -0.5, x0=n_list[2], y0=sd_est[2] * 0.35)
ax2.legend()
ph.savefig(fig, "mc_errorbars")
print(f"fitted log-log slope {slope:.3f} (theory -0.5); SD*sqrt(n) at each n: {np.round(sd_est * np.sqrt(n_list), 2)}  (theory {np.sqrt(np.pi * (4 - np.pi)):.2f})")

dart = 4.0 * (np.sum(rng5.random((1000, 2)) ** 2, axis=1) < 1)         # ONE run of n = 1000 darts
sd1 = dart.std(ddof=1)
se1 = sd1 / np.sqrt(1000)
ci1 = stats.t.ppf(0.975, 999) * se1                                    # (3.9)
print(f"one run, n=1000: SD of a single dart {sd1:.3f}, SE of the mean {se1:.4f}, 95% CI half-width {ci1:.4f}")
```
```text
fitted log-log slope -0.501 (theory -0.5); SD*sqrt(n) at each n: [1.74 1.61 1.57 1.69 1.61 1.64 1.7  1.62]  (theory 1.64)
one run, n=1000: SD of a single dart 1.674, SE of the mean 0.0529, 95% CI half-width 0.1039
```

Left: the means sit on $\pi$ and the bars shrink; the legend says what a bar *means* (one SD of the estimator over 400 runs). Right: that SD on log-log axes, theory underneath, a triangle of slope $-1/2$. The fitted slope is $-0.501$ (asserted: $-0.5\pm0.05$) and $\mathrm{SD}\times\sqrt n$ is flat near 1.64, as (4.2) says. Halving the bar costs four times the darts. The last lines are §3.9's distinction: from *one* run, the SD of a single dart (1.67) and the SE of the mean (0.053) differ by $\sqrt{1000}$, and the SE is what the SD across repeated runs estimates. The right panel's bars are $\mathrm{SD}/\sqrt{2(R-1)}$, the standard error of an SD from $R$ replicates of a roughly normal quantity.

### Example 4.5: A 2D potential: `imshow` versus `pcolormesh`, contours, and the pixel test (EM-2)

The softened potential of charges $+q$ at $(0.5, 0.3)$ and $-q$ at $(-0.5, -0.3)$, in units of $q/4\pi\epsilon_0\ell$ for a length unit $\ell$ (formula in the code). The tilt matters: a dipole on the $x$ axis is mirror-symmetric in $y$ and would hide a vertical flip. `np.meshgrid(x, y)` gives `X, Y` of shape `(len(y), len(x))`, so `V[i, j]` belongs to `y[i], x[j]`: **row index is $y$**. Grid points are pixel centres, so by (3.4) the extent reaches half a cell beyond them.

▶ **Try it first.** What are the extent and origin of a bare `imshow(V)` for `V` of shape `(61, 81)`? For physical data, which is right: `extent=...` with the default origin, or with `origin="lower"`?

```python
x = np.linspace(-2, 2, 81)
y = np.linspace(-1.5, 1.5, 61)
X, Y = np.meshgrid(x, y)
V = 1 / np.sqrt((X - 0.5) ** 2 + (Y - 0.3) ** 2 + 0.05**2) - 1 / np.sqrt((X + 0.5) ** 2 + (Y + 0.3) ** 2 + 0.05**2)
dx, dy = x[1] - x[0], y[1] - y[0]
extent = (x[0] - dx / 2, x[-1] + dx / 2, y[0] - dy / 2, y[-1] + dy / 2)       # outer EDGES of the pixels
norm, cmap = Normalize(-4, 4), plt.get_cmap("RdBu_r")
fig, ax = plt.subplots()
print("bare imshow(V): extent", ax.imshow(V).get_extent())
plt.close(fig)


def draw_potential(ax, how):
    """Draw V three ways on `ax`: imshow, pcolormesh, or imshow without origin='lower'."""
    ax.grid(False)                                          # the house grid would paint over the pixels
    if how == "imshow":
        return ax.imshow(V, origin="lower", extent=extent, cmap=cmap, norm=norm, aspect="equal")
    if how == "pcolormesh":
        return ax.pcolormesh(X, Y, V, shading="nearest", cmap=cmap, norm=norm)
    return ax.imshow(V, extent=extent, cmap=cmap, norm=norm, aspect="equal")


rng6 = np.random.default_rng(6)
rows, cols = rng6.integers(5, 56, 25), rng6.integers(5, 76, 25)          # 25 random grid points
worst = {}
for how in ("imshow", "pcolormesh", "flipped"):
    fig, ax = plt.subplots(figsize=(7, 5))
    draw_potential(ax, how)
    ax.set(xlim=(x[0], x[-1]), ylim=(y[0], y[-1]))
    worst[how] = np.abs(ph.rgba_at(ax, x[cols], y[rows]) - cmap(norm(V[rows, cols]))).max()       # colour under each grid point vs cmap(norm(V))
    plt.close(fig)
print("worst colour error over 25 grid points:", {k: float(f"{v:.2g}") for k, v in worst.items()})
```
```text
bare imshow(V): extent [-0.5, 80.5, 60.5, -0.5]
worst colour error over 25 grid points: {'imshow': 0.0038, 'pcolormesh': 0.0018, 'flipped': 0.93}
```

`rgba_at` (in `plot_helpers.py`) draws the figure, maps data points to pixels with `ax.transData`, and reads the rendered RGBA buffer there (display $y$ starts at the bottom, the buffer at the top, so the row is flipped). `imshow` and `pcolormesh` match the colour map to within 0.004, one of its 256 steps; the image without `origin="lower"` is off by 0.93 in a channel: upside down. The bare `imshow(V)` extent printed above is in pixel *indices*, $y$ increasing downward: nothing to do with your data. Now the half-pixel rule on a coarse $9\times7$ grid, then the finished figure:

```python
xc, yc = np.linspace(-2, 2, 9), np.linspace(-1.5, 1.5, 7)
Xc, Yc = np.meshgrid(xc, yc)
Vc = np.sign(Xc + 0.1) * (Yc + 2.0)
for label, ext in [("edges", (-2.25, 2.25, -1.75, 1.75)), ("centres", (-2, 2, -1.5, 1.5))]:
    fig, ax = plt.subplots()
    ax.imshow(Vc, origin="lower", extent=ext, cmap=cmap, norm=norm)
    ax.grid(False)
    ax.set(xlim=(-2.25, 2.25), ylim=(-1.75, 1.75))
    bad = np.abs(ph.rgba_at(ax, Xc.ravel(), Yc.ravel()) - cmap(norm(Vc.ravel()))).max(axis=1) > 1e-2
    print(f"extent through the pixel {label}: {bad.sum()} of 63 centre tests fail")
    plt.close(fig)

fig, ax = plt.subplots(figsize=(6.5, 4.4))
im = draw_potential(ax, "imshow")
cs = ax.contour(X, Y, V, levels=[-3, -2, -1, -0.5, 0.5, 1, 2, 3], colors="k", linewidths=0.6)
ax.clabel(cs, fmt="%g", fontsize=7)
ax.set(xlabel="$x$ (length unit $\\ell$)", ylabel="$y$ (length unit $\\ell$)")
cb = fig.colorbar(im, ax=ax, extend="both", label="potential $V$ (units of $q / 4\\pi\\epsilon_0 \\ell$)")
print("colorbar axis limits:", tuple(float(v) for v in cb.ax.get_ylim()), "| inverted:", cb.ax.yaxis_inverted())
ph.savefig(fig, "potential_notes")
```
```text
extent through the pixel edges: 0 of 63 centre tests fail
extent through the pixel centres: 9 of 63 centre tests fail
colorbar axis limits: (-4.0, 4.0) | inverted: False
```

An extent through the pixel *centres* shifts every cell by half a cell: 9 of 63 tests fail on the coarse grid, none with edges. On the fine grid the shift is 0.025, easy to miss, and the kind of error that puts a field line half a cell from a charge. The colorbar runs $-4\to4$ upward, with `extend="both"` arrows because $|V|$ exceeds 4 near the charges (such values take the end colours). Use `imshow` for fast uniform grids, `pcolormesh` for non-uniform grids or explicit coordinates; contour labels give the values colours cannot. `project.py` draws all three variants side by side and checks the zero contour (the charges' perpendicular bisector) and the colorbar's top colour.

### Example 4.6: A truncated axis and its lie factor

Two alloys, tensile strength $412 \pm 3$ and $421 \pm 3$ MPa (±1 standard error). Draw the *same* bars twice.

```python
labels, means, se = ["alloy A", "alloy B"], np.array([412.0, 421.0]), np.array([3.0, 3.0])
fig, axs = plt.subplots(2, 1, figsize=(5, 6), sharex=True)
for ax, ylim in zip(axs, [(400, 425), (0, 450)]):
    ax.bar(labels, means, yerr=se, capsize=5, color=[ph.OKABE_ITO[0], ph.OKABE_ITO[1]], width=0.5)
    ax.set(ylim=ylim, ylabel="tensile strength (MPa)")
axs[0].set_title("truncated axis (starts at 400 MPa)")
axs[1].set_title("zero baseline")
ph.savefig(fig, "truncated_vs_zero")

shown = (means[1] - 400) / (means[0] - 400) - 1               # relative bar-length difference as drawn on the truncated axis
actual = means[1] / means[0] - 1
print(f"data: B is {actual:.1%} above A; truncated bars show {shown:.0%}; lie factor (3.8) = {shown / actual:.0f}")
print(f"difference / its standard error = {(means[1] - means[0]) / np.hypot(*se):.2f}")
```
```text
data: B is 2.2% above A; truncated bars show 75%; lie factor (3.8) = 34
difference / its standard error = 2.12
```

The top panel draws bar B 75% longer than bar A for a 2.2% difference in the data: lie factor 34. Each panel supports a different claim. **Top:** "B's mean is higher than A's", an *ordering* claim; the difference is 2.1 standard errors of the difference ($9/\sqrt{3^2+3^2}$), so even this wants modest wording. **Bottom:** "the alloys differ by about 2%", a *magnitude* claim. If the reader needs both, replace the top panel by a dot plot (position encodes the value, so a stated truncated range is legitimate) or by the difference with its interval.

### Example 4.7: The pendulum phase portrait (Build 3; CM-1 preview)

$\ddot\theta = -\sin\theta$ in time units of $\sqrt{L/g}$, state $(\theta, \omega)$, $\omega = \dot\theta$. The energy $E = \tfrac12\omega^2 + (1-\cos\theta)$ is conserved: $\dot E = \omega\dot\omega + \sin\theta\,\dot\theta = -\omega\sin\theta + \omega\sin\theta = 0$, so trajectories are level sets of $E$. The **separatrix** is the level through the unstable rest point $(\pm\pi, 0)$, $E = 1-\cos\pi = 2$; there $\tfrac12\omega^2 = 1+\cos\theta = 2\cos^2(\theta/2)$, so $\omega = \pm2\cos(\theta/2)$. Below it the pendulum swings (libration); above it, it goes round (rotation).

▶ **Try it first.** For release states $(\theta_0, \omega_0) = (0.5, 0),\ (2.5, 0),\ (0, 2.3)$, predict each energy and which are librations. Build 3 in three steps: (1) `quiver` on a coarse grid; (2) `streamplot` on a fine one; (3) trajectories from `solve_ivp` on top. One subtlety in step 1: $\theta$ and $\omega$ have different units and different lengths per unit on the page, so an arrow drawn at the data-space angle of $(\dot\theta, \dot\omega)$ (`angles="uv"`) is not tangent to the curves as drawn. `angles="xy"` makes Matplotlib point each arrow from $(\theta, \omega)$ to $(\theta + \dot\theta,\ \omega + \dot\omega)$ *through the data-to-display transform*, which is tangent whatever the panel's shape; give $(\dot\theta, \dot\omega)$ unit length first so every arrow is the same size. (Dividing the components by the axis ranges instead, "axes fractions", is tangent only on a square panel; the block measures how far off that would be here.)

```python
TH_LIM, W_LIM = 2 * np.pi, 3.0
field = lambda th, w: (w, -np.sin(th))                          # (d theta/dt, d omega/dt)
energy = lambda th, w: 0.5 * w**2 + 1 - np.cos(th)

starts = {"small swing": (0.5, 0.0), "large swing": (2.5, 0.0), "rotation": (0.0, 2.3)}
traj = {}
for name, (th0, w0) in starts.items():
    sol = solve_ivp(lambda t, s: field(s[0], s[1]), (0, 14), [th0, w0], method="DOP853", rtol=1e-12, atol=1e-13, dense_output=True)
    traj[name] = sol.sol(np.linspace(0, 14, 1400))
    drift = np.max(np.abs(energy(*traj[name]) - energy(th0, w0))) / energy(th0, w0)
    print(f"{name:11s} E0 = {energy(th0, w0):.4f}, max relative energy drift {drift:.1e}")
    assert drift < 1e-6

th_g, w_g = np.meshgrid(np.linspace(-TH_LIM, TH_LIM, 25), np.linspace(-W_LIM, W_LIM, 13))
u, v = field(th_g, w_g)
n_ = np.maximum(np.hypot(u, v), 1e-9)                           # the field vanishes at equilibria: do not divide by 0
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
ax1.quiver(th_g, w_g, u / n_, v / n_, angles="xy", pivot="mid", color="0.3")   # unit length; direction computed on the page
th_f, w_f = np.meshgrid(np.linspace(-TH_LIM, TH_LIM, 240), np.linspace(-W_LIM, W_LIM, 120))
uf, vf = field(th_f, w_f)
stream = ax2.streamplot(th_f, w_f, uf, vf, color=np.hypot(uf, vf), cmap="cividis", density=1.3, linewidth=0.9)
fig.colorbar(stream.lines, ax=ax2, label="phase-plane speed $\\sqrt{\\dot\\theta^2 + \\dot\\omega^2}$ (time in units of $\\sqrt{L/g}$)")
tg = np.linspace(-TH_LIM, TH_LIM, 400)
for sign in (1, -1):
    ax2.plot(tg, sign * 2 * np.cos(tg / 2), "k--", lw=0.8)                  # the separatrix E = 2
for (name, (th_t, w_t)), c in zip(traj.items(), ph.OKABE_ITO[1:]):
    ax2.plot(th_t, w_t, color=c, lw=2, label=name)
for a in (ax1, ax2):
    a.set(xlim=(-TH_LIM, TH_LIM), ylim=(-W_LIM, W_LIM), xlabel="angle $\\theta$ (rad)")
    ph.pi_ticks(a, "x", np.pi)
ax1.set_ylabel("angular velocity $\\omega$ (rad per $\\sqrt{L/g}$)")
fig.canvas.draw()
box = ax1.get_window_extent()
pts = np.column_stack([th_g.ravel(), w_g.ravel()])
page = ax1.transData.transform(pts + 1e-6 * np.column_stack([u.ravel(), v.ravel()])) - ax1.transData.transform(pts)
on_page = np.degrees(np.arctan2(page[:, 1], page[:, 0]))                                  # what angles="xy" draws
axes_frac = np.degrees(np.arctan2(v.ravel() / (2 * W_LIM), u.ravel() / (2 * TH_LIM)))     # the shortcut: right only on a square panel
off = np.abs((on_page - axes_frac + 180) % 360 - 180)[n_.ravel() > 1e-6]
print(f"panel is {box.width / box.height:.2f} times wider than tall; axes-fraction arrows would be off by up to {off.max():.1f} degrees")
ph.savefig(fig, "phase_portrait_notes")
```
```text
small swing E0 = 0.1224, max relative energy drift 5.3e-12
large swing E0 = 1.8011, max relative energy drift 2.9e-12
rotation    E0 = 2.6450, max relative energy drift 3.6e-12
panel is 1.33 times wider than tall; axes-fraction arrows would be off by up to 8.1 degrees
```

The last line is the promised measurement: this panel is 1.33 times wider than tall, so arrows built from axes fractions would miss the tangent by up to about 8°; with `angles="xy"` the error is zero to rounding. (`project.py` also asserts the large swing never exceeds its release angle and the rotation never reverses.) Closed loops around $(0,0)$: the small swing nearly an ellipse (the harmonic limit), the large swing squared off; saddles at $(\pm\pi, 0)$ where the separatrix crosses; the rotation *above* the separatrix, wavy but never turning back. The drift of $E$ is below $6\times10^{-12}$, far inside the $10^{-6}$ asserted, so the curves are level sets, not integrator artifacts. Streamlines come from the field alone, trajectories are integrated; where they overlay, the agreement is the check.

## 5. Why it matters and where it connects

- **CM-1, CM-3 to CM-8.** Phase portraits (§4.7) and energy-drift asserts are the standard check on an integrator; CM-5's orbits need `set_aspect("equal")`.
- **EM-1, EM-2.** `quiver` and `streamplot` draw $\mathbf E = -\nabla V$; §4.5 is the EM-2 figure, and its `origin`/`extent` test is the check before believing any field map.
- **F3, F4, F5.** Power laws are straight on log axes (§3.4); histograms with densities and KDEs (§4.3); the $1/\sqrt n$ plot (§4.4) is the figure of F5's limit theorems.
- **ST-1 onward, ML-1 onward.** Intervals, bootstraps and learning curves depend on whether a bar says SD, SE or CI (§3.9); weight-matrix heat maps need §3.6's colour rules.
- **LA-13, PY-8 to PY-10.** Convergence and timing plots reuse `plot_convergence`; `DataFrame.plot` returns an Axes, so §3.1 is how you finish a pandas figure.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 drill the object model, scales, `imshow` placement, tick formatters, norms and colour-map choice. B1–B7 are the module and all are essential: multi-series `plot_convergence`, a V-shaped error curve, a limit cycle, SD/SE/CI and KDE bandwidth, colour maps by lightness, slide and paper versions, small multiples with a shared norm. C1–C3: a misleading figure and its fix, a direct-label collision resolver, a file-size and timing study.

## 7. Coding project

`plot_helpers.py` plus `project.py`. The script builds and asserts: (1) the anatomy figure and artist types; (2) the convergence plot, slopes within 0.1 of 2 and 4, leading constants, the $R^2$ table, CSV/JSON round trips; (3) the pendulum phase portrait, energy conserved to $10^{-6}$ relative; (4) the histogram, Gaussian fit and KDE; (5) Monte Carlo error bars, slope $-0.5\pm0.05$; (6) the potential map with the pixel test of `origin`/`extent` and of the colorbar; (7) six damped-oscillator small multiples against closed forms; (8) `jet` against `viridis` with $L^*$; (9) the pendulum animation (GIF via Pillow, every frame's data checked); (10) the truncated axis and its lie factor; (11) the numbers quoted in §3. Asserts check data, not pixels, except the `extent` check. Run `MPLBACKEND=Agg python3 project.py` (about half a minute).

**Build it yourself:** before reading `plot_helpers.py`, write `savefig`, `add_slope_guide` and `plot_convergence` (Builds 1 and 5) and run them on your PY-5 errors. Then the "Your turn" cell: slide and paper versions of one figure, an inset zoom, a deliberately misleading figure (claim stated, then fixed), and the PY-6 random-walk MSD with its guide line.

## 8. Where you'll struggle

- **The pyplot state machine works until you have two figures.** `plt.plot`, `plt.xlabel` draw on whatever figure is "current" (§3.1 printed the proof). It feels faster for a week. Switch to `ax` now; keep only `plt.subplots`, `plt.close(fig)` and `plt.rcParams`.
- **`imshow` puts row 0 at the top and your physics puts $y=0$ at the bottom.** Give every `imshow` of physical data `origin="lower"` *and* an `extent` from pixel edges (§4.5), or you plot a mirror image half a cell off. A symmetric function hides a flip, so test with an asymmetric one.
- **Log axes drop zero and negative data without a warning** (§3.4). Check `np.all(y > 0)` before `loglog`; use `symlog` if the sign matters.
- **`tight_layout` and `constrained_layout` fight.** Keep constrained layout on (the house style does) and never call `tight_layout()` or `subplots_adjust` after it. Hand-placed axes need `layout="none"` (§3.9).
- **The grid is drawn over images**, and `bbox_inches="tight"` changes the pixel size. Both are one-line fixes, found only by looking at the bad figure.
- **A legend is a failure of direct labelling that we tolerate.** Label the line unless lines collide; keep legend order equal to line order and the legend off the data.
- **`jet` is still the default in people's heads.** Colleagues will ask for "the rainbow one"; show them the $L^*$ plot (§3.6).
- **A figure that needs a paragraph to explain has two messages.** Split it.
- **Headless saving is what this curriculum requires.** `plt.show()` in a script hangs on an interactive backend or does nothing under `Agg`. Call `fig.savefig`, close the figure (`plt.close(fig)`; Matplotlib warns when more than 20 stay open), and open the saved file.
- **Global state leaks.** `plt.rcParams[...] = ...` in one function changes every later figure: use `with plt.rc_context(...)`.

## 9. Check yourself

1. What do `plt.subplots(2, 3)` and `plt.subplots()` return, and what are the type and shape of the second item in each?
2. Two figures exist and a script calls `plt.plot(x, y)`. Where does the line go, and what is the fix?
3. If $e = Ch^p$, which axes make a straight line, with what slope? How do you fit it, and why in log space?
4. What does `LogNorm(1, 1000)` map 10 to? Why must the data be positive?
5. For `V` of shape `(ny, nx)` sampled at pixel centres `x`, `y`, which `origin` and `extent` put `V[i, j]` at `(x[j], y[i])`?
6. What does `density=True` normalise, and what do $n$ equal bars of width $w$ sum to?
7. Name two ways `jet` misleads, and the $L^*$ property of `viridis` that prevents them.
8. When is `twinx` misleading, and what do you do instead?
9. What do SD, SE and a 95% CI half-width describe, and which shrink like $1/\sqrt n$?
10. Why may a bar chart not be truncated? What is the lie factor of bars drawn from 400 for 412 and 421?

## 10. Pair with

- **Matplotlib "Quick start guide" and "Artist tutorial" (matplotlib.org).** The Quick start is the OO API in §3.1's order; the Artist tutorial is the authoritative account of the artist tree and of transforms.
- **Rougier, *Scientific Visualization: Python + Matplotlib* (free), ch. 1–4.** How the rendering pipeline works and how to make figures look designed.
- **VanderPlas, *Python Data Science Handbook* 2e, ch. 4.** A catalogue of recipes (legends, colour bars, subplots, annotations, 3D) to adapt.
- **Tufte, *The Visual Display of Quantitative Information*, ch. 1–2.** Graphical integrity and the lie factor of (3.8); read it for judgment about what a graphic may claim.
- **Crameri, Shephard & Heron, "The misuse of colour in science communication," *Nature Communications* 11 (2020).** The argument against rainbow maps, with real examples; §3.6's $L^*$ plot is a small piece of it.
- **Wilke, *Fundamentals of Data Visualization* (free online), chapters on axes, colour and uncertainty.** What an error bar should say and how to choose a map by data type.
