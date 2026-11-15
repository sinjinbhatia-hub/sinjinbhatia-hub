# %% [markdown]
# # PY-7 Project: Canonical Scientific Figures
#
# Builds the figures the curriculum keeps drawing and checks the *data* behind each (notes
# sections 3 and 4): (1) anatomy of a figure, (2) a convergence plot with fitted slopes, (3) the
# pendulum phase portrait with conserved energy, (4) histogram, Gaussian fit and KDE, (5) Monte
# Carlo error bars with a 1/sqrt(n) guide, (6) a 2D potential with the pixel test of `origin` and
# `extent`, (7) small multiples of damped oscillators, (8) jet against viridis, (9) a pendulum
# animation, (10) a truncated axis and its lie factor. Helpers: `plot_helpers.py` (notes 3.7).
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.image as mimg
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import NullFormatter
from scipy import stats
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

import plot_helpers as ph

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
ph.use_house_style()

# %% [markdown]
# ## Part 1: Anatomy of a figure (notes Example 4.1)
#
# Every artist knows its box; `annotate(xycoords=artist)` puts a numbered marker on it.

# %%
fig, ax = plt.subplots(figsize=(8.5, 4.5))
t = np.linspace(0, 2 * np.pi, 100)
(curve,) = ax.plot(t, np.sin(t), label="$x = \\sin t$")
ax.set(xlabel="time $t$ (rad)", ylabel="displacement $x$ (m)", title="title")
ax.legend(loc="lower left")
fig.canvas.draw()  # artists report their boxes only after a draw
named = {"Axes: one plot, with its own patch": (ax.patch, (0.95, 0.9)), "Line2D: the curve": (curve, (0.38, 1.0)),
         "Legend": (ax.get_legend(), (1.0, 0.5)), "title (a Text)": (ax.title, (1.0, 0.5)),
         "xlabel (a Text)": (ax.xaxis.label, (1.0, 0.5)), "ylabel (a Text)": (ax.yaxis.label, (0.5, 1.0)),
         "tick label (a Text)": (ax.get_xticklabels()[3], (1.0, 0.5)), "Spine: a frame line": (ax.spines["bottom"], (0.9, 0.5)),
         "XAxis: ticks, tick labels, label": (ax.xaxis, (0.5, 0.0)), "Figure: the canvas": (fig.patch, (0.02, 0.97))}
for k, (text, (artist, frac)) in enumerate(named.items(), 1):
    ax.annotate(str(k), xy=frac, xycoords=artist, xytext=(12, 0), textcoords="offset points", ha="center", va="center", color="white",
                fontsize=8, bbox={"boxstyle": "circle,pad=0.25", "fc": "0.2", "ec": "none"}, annotation_clip=False)
    ax.text(1.03, 1.02 - 0.1 * k, f"{k}  {text}", transform=ax.transAxes, fontsize=8, va="center")
ph.savefig(fig, "anatomy")
kinds = {type(a).__name__ for a in ax.get_children()}
assert {"Line2D", "XAxis", "YAxis", "Spine", "Legend", "Text", "Annotation"} <= kinds, kinds
assert ax.xaxis.label.get_text() == "time $t$ (rad)" and len(ax.lines) == 1 and ax.lines[0] is curve and fig.axes == [ax]
print("Part 1: ax children are", sorted(kinds))

# %% [markdown]
# ## Part 2: The convergence plot (notes Example 4.2, Builds 1 and 2)
#
# Trapezoid and Simpson for the integral of sin over [0, pi] (exact value 2), ten step sizes.
# Leading errors (PY-5 section 3.10 eq. (3.2), with (b - a) f''(xi) read as the Riemann sum for the
# integral of f''): T - I ~ (h^2/12)(f'(b) - f'(a)) = -h^2/6 and S - I ~ (h^4/180)(f'''(b) - f'''(a)) = h^4/90 for f = sin.


# %%
def trapezoid(f, a, b, n):
    """Composite trapezoid rule, h (y0/2 + y1 + ... + y_{n-1} + yn/2); vectorized over x."""
    x = np.linspace(a, b, n + 1)
    y = f(x)
    return (b - a) / n * (y.sum() - 0.5 * (y[0] + y[-1]))


def simpson(f, a, b, n):
    """Composite Simpson rule, (h/3)(y0 + 4 y1 + 2 y2 + ... + 4 y_{n-1} + yn); n even."""
    x = np.linspace(a, b, n + 1)
    y = f(x)
    return (b - a) / n / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


ns = np.array([4, 6, 8, 12, 16, 24, 32, 48, 64, 96])
h = np.pi / ns
err_t = np.array([abs(trapezoid(np.sin, 0, np.pi, n) - 2) for n in ns])
err_s = np.array([abs(simpson(np.sin, 0, np.pi, n) - 2) for n in ns])

# Cross-check against the PY-5 scalar implementations when the sibling folder is present.
sib = HERE.parent / "05-errors-testing-debugging"
if (sib / "integrators.py").exists():
    sys.path.insert(0, str(sib))
    import integrators as ig

    for n in ns:
        for mine, theirs in ((trapezoid, ig.trapezoid), (simpson, ig.simpson)):
            assert math.isclose(theirs(math.sin, 0.0, math.pi, int(n)), mine(np.sin, 0, np.pi, n), rel_tol=1e-12)
    sys.path.remove(str(sib))
    print("Part 2: matches PY-5 integrators.py at all ten n")

fig, ax = plt.subplots()
p_t = ph.plot_convergence(ax, h, err_t, 2, label="trapezoid")
p_s = ph.plot_convergence(ax, h, err_s, 4, label="Simpson")
ax.legend(loc="lower right")
ph.savefig(fig, "convergence")
fig, ax = plt.subplots()                                       # add_slope_guide(ax, order) with the default corner
(pts,) = ax.loglog(h, err_t, "o")
ph.add_slope_guide(ax, 2)
xg, yg = ax.lines[1].get_xydata()[[0, 2]].T
assert math.isclose(math.log10(yg[1] / yg[0]) / math.log10(xg[1] / xg[0]), 2.0)   # the hypotenuse has slope 2
ph.label_end(ax, pts, "trapezoid")                             # direct label: an Annotation at the line's last point, in its colour
tag = ax.texts[-1]
assert tag.get_text() == "trapezoid" and tag.get_color() == pts.get_color() and tuple(tag.xy) == tuple(pts.get_xydata()[-1])
try:                                                           # notes 3.4: an error of exactly 0 must raise, not vanish from a log axis
    ph.plot_convergence(ax, h, np.append(err_t[:-1], 0.0), 2)
except ValueError as exc:
    assert "err > 0" in str(exc)
else:
    raise AssertionError("plot_convergence accepted an error of exactly 0")
plt.close(fig)
assert abs(p_t - 2) < 0.1 and abs(p_s - 4) < 0.1
c_t, c_s = err_t[-1] / (h[-1] ** 2 / 6), err_s[-1] / (h[-1] ** 4 / 90)           # leading constants: both tend to 1 as h -> 0
assert abs(c_t - 1) < 2e-5 and abs(c_s - 1) < 1.3e-4                           # notes Example 4.2: within 2e-5 and 1.3e-4
print(f"Part 2: fitted slopes: trapezoid {p_t:.3f}, Simpson {p_s:.3f}; err/(h^2/6) at n=96: {c_t:.5f}, err/(h^4/90): {c_s:.5f}")


# %% [markdown]
# Same data on linear, semilog and log-log axes (Build 2). "How straight" is the coefficient
# of determination $R^2$ of a line fitted *in the plotted coordinates*.

# %%
def r_squared(x, y):
    """R^2 of the least-squares line y ~ a x + b."""
    res = y - np.polyval(np.polyfit(x, y, 1), x)
    return 1 - res.var() / y.var()


scales = {"linear": (lambda v: v, lambda v: v), "semilogy": (lambda v: v, np.log10), "log-log": (np.log10, np.log10)}
r2 = {k: (r_squared(fx(h), fy(err_t)), r_squared(fx(h), fy(err_s))) for k, (fx, fy) in scales.items()}
for k, (a, b) in r2.items():
    print(f"Part 2: R^2 on {k:8s} axes: trapezoid {a:.4f}, Simpson {b:.4f}")
assert r2["log-log"][0] > 0.9999 and r2["log-log"][1] > 0.9999 and r2["linear"][1] < 0.8 and r2["semilogy"][0] < 0.9
frac = err_s[-1] / err_t.max()
print(f"Part 2: Simpson at n=96 is {frac:.1e} of the trapezoid axis height: {frac * 400:.1e} pixels on a 400 px axis")

fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
for ax, title in zip(axs, scales):
    ax.plot(h, err_t, "o-", label="trapezoid")
    ax.plot(h, err_s, "s-", label="Simpson")
    ax.set(title=title, xlabel="step size $h$", ylabel="absolute error")
axs[1].set_yscale("log")
axs[2].set(xscale="log", yscale="log")
axs[2].xaxis.set_minor_formatter(NullFormatter())  # a short log axis labels its minor ticks and they collide
axs[0].legend()
ph.savefig(fig, "three_scales")

# Reading the data back: CSV (csv module and loadtxt) and JSON (PY-2) agree to the last digit.
table = FIG / "convergence.csv"
with table.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["n", "h", "err_trapezoid", "err_simpson"])
    w.writerows(zip(ns.tolist(), h.tolist(), err_t.tolist(), err_s.tolist()))
back = np.loadtxt(table, delimiter=",", skiprows=1)
with table.open(newline="") as fh:
    rows = list(csv.DictReader(fh))
(FIG / "convergence.json").write_text(json.dumps({"h": h.tolist(), "err_t": err_t.tolist()}))
js = json.loads((FIG / "convergence.json").read_text())
assert back.shape == (10, 4) and np.array_equal(back[:, 2], err_t) and np.array_equal(back[:, 3], err_s)
assert float(rows[3]["err_simpson"]) == err_s[3] and np.array_equal(np.array(js["err_t"]), err_t)
print("Part 2: CSV (csv + loadtxt) and JSON round trips are exact")

# %% [markdown]
# ## Part 3: Pendulum phase portrait (notes Example 4.7, Build 3; CM-1 preview)
#
# theta'' = -sin(theta), state (theta, omega); E = omega^2/2 + 1 - cos(theta) is conserved and the
# separatrix is E = 2, i.e. omega = +-2 cos(theta/2). Quiver arrows use angles="xy": the direction is
# computed through the data-to-display transform, so they are tangent on the page although theta and
# omega have unlike units and the panel is not square (dividing by the axis ranges is right only if it is).


# %%
def pendulum_field(th, w):
    """(d theta/dt, d omega/dt) = (omega, -sin theta)."""
    return w, -np.sin(th)


def energy(th, w):
    """E = omega^2/2 + 1 - cos(theta), zero at rest at the bottom."""
    return 0.5 * w**2 + 1 - np.cos(th)


TH_LIM, W_LIM = 2 * np.pi, 3.0
starts = {"small swing": (0.5, 0.0), "large swing": (2.5, 0.0), "rotation": (0.0, 2.3)}
traj = {}
for name, (th0, w0) in starts.items():
    sol = solve_ivp(lambda t, s: pendulum_field(s[0], s[1]), (0, 14), [th0, w0], method="DOP853",
                    rtol=1e-12, atol=1e-13, dense_output=True)
    tt = np.linspace(0, 14, 1400)
    traj[name] = sol.sol(tt)
    drift = np.max(np.abs(energy(*traj[name]) - energy(th0, w0))) / energy(th0, w0)
    print(f"Part 3: {name:11s} E0 = {energy(th0, w0):.4f}, max relative energy drift {drift:.1e}")
    assert drift < 1e-6
assert np.max(np.abs(traj["large swing"][0])) < 2.5 + 1e-6          # a libration never exceeds its release angle
assert np.all(traj["rotation"][1] > 0) and energy(0.0, 2.3) > 2     # above the separatrix: always turning one way

th_g, w_g = np.meshgrid(np.linspace(-TH_LIM, TH_LIM, 25), np.linspace(-W_LIM, W_LIM, 13))
u, v = pendulum_field(th_g, w_g)
norm_ = np.maximum(np.hypot(u, v), 1e-9)  # the field vanishes at the equilibria; do not divide by 0
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
ax1.quiver(th_g, w_g, u / norm_, v / norm_, angles="xy", pivot="mid", color="0.3")  # unit length; page direction
ax1.set_title("quiver: direction only")
th_f, w_f = np.meshgrid(np.linspace(-TH_LIM, TH_LIM, 240), np.linspace(-W_LIM, W_LIM, 120))
uf, vf = pendulum_field(th_f, w_f)
stream = ax2.streamplot(th_f, w_f, uf, vf, color=np.hypot(uf, vf), cmap="cividis", density=1.3, linewidth=0.9)
fig.colorbar(stream.lines, ax=ax2, label="phase-plane speed $\\sqrt{\\dot\\theta^2 + \\dot\\omega^2}$ (time in units of $\\sqrt{L/g}$)")
tgrid = np.linspace(-TH_LIM, TH_LIM, 400)
for sign in (1, -1):
    ax2.plot(tgrid, sign * 2 * np.cos(tgrid / 2), "k--", lw=0.8)
ax2.text(0.02, 0.93, "separatrix $E=2$", transform=ax2.transAxes, fontsize=8)
offsets = {"small swing": (0, 14), "large swing": (6, -16), "rotation": (-4, 8)}
for (name, (th_t, w_t)), c in zip(traj.items(), ph.OKABE_ITO[1:]):
    ax2.plot(th_t, w_t, color=c, lw=2)
    xy = (0.0, w_t.max()) if name == "small swing" else (th_t[0], w_t[0])       # label each trajectory in place
    ax2.annotate(name, xy=xy, xytext=offsets[name], textcoords="offset points", color=c, fontsize=8,
                 ha="right" if name == "rotation" else "left" if name == "large swing" else "center",
                 bbox={"fc": "white", "ec": "none", "pad": 0.5, "alpha": 0.8})
for a in (ax1, ax2):
    a.set(xlim=(-TH_LIM, TH_LIM), ylim=(-W_LIM, W_LIM), xlabel="angle $\\theta$ (rad)")
    ph.pi_ticks(a, "x", np.pi)
ax1.set_ylabel("angular velocity $\\omega$ (rad per $\\sqrt{L/g}$)")
ax2.set_title("streamplot and three trajectories")
fig.canvas.draw()                                             # the axes-fraction shortcut is off unless the panel is square: measure it
box = ax1.get_window_extent()
pts = np.column_stack([th_g.ravel(), w_g.ravel()])
page = ax1.transData.transform(pts + 1e-6 * np.column_stack([u.ravel(), v.ravel()])) - ax1.transData.transform(pts)
on_page = np.degrees(np.arctan2(page[:, 1], page[:, 0]))                                  # what angles="xy" draws
axes_frac = np.degrees(np.arctan2(v.ravel() / (2 * W_LIM), u.ravel() / (2 * TH_LIM)))     # the shortcut
off = np.abs((on_page - axes_frac + 180) % 360 - 180)[norm_.ravel() > 1e-6]
k = box.height / box.width                                    # tan(page angle) = k tan(shortcut angle); the gap peaks at tan = 1/sqrt(k)
bound = np.degrees(np.arctan(1 / np.sqrt(k)) - np.arctan(np.sqrt(k)))
assert k < 0.85 and 5 < off.max() <= bound + 1e-9
print(f"Part 3: panel W/H = {1 / k:.2f}; axes-fraction arrows would be off by up to {off.max():.1f} deg (bound {bound:.1f})")
ph.savefig(fig, "phase_portrait")

# %% [markdown]
# ## Part 4: Histogram, fitted Gaussian and KDE (notes Example 4.3; F4)

# %%
true = stats.gamma(4, scale=1.5)
sample = true.rvs(1500, random_state=np.random.default_rng(0))
mu, sd = stats.norm.fit(sample)
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

assert math.isclose(np.sum(counts * np.diff(edges)), 1.0, rel_tol=1e-12)     # density=True: total AREA is 1
assert math.isclose(kde.factor, len(sample) ** (-1 / 5)) and math.isclose(np.trapezoid(kde(xs), xs), 1.0, abs_tol=2e-3)
tv = lambda f: 0.5 * np.trapezoid(np.abs(f - true.pdf(xs)), xs)
tv_gauss, tv_kde = tv(stats.norm.pdf(xs, mu, sd)), tv(kde(xs))
assert tv_kde < tv_gauss and 0.12 < tv_gauss < 0.14 and 0.015 < tv_kde < 0.025        # notes Example 4.3: 0.132 and 0.019
print(f"Part 4: hist area {np.sum(counts * np.diff(edges)):.6f}; Scott factor n^(-1/5) = {kde.factor:.4f}, "
      f"kernel sd = {kde.factor * sample.std(ddof=1):.3f} s")
print(f"        total-variation distance to the true density: Gaussian fit {tv_gauss:.3f}, KDE {tv_kde:.3f}")

# %% [markdown]
# ## Part 5: Monte Carlo error bars (notes Example 4.4; F5)
#
# Hit-or-miss estimate of pi from n darts. Over R replicates the estimates scatter with
# standard deviation sqrt(pi (4 - pi) / n): slope -1/2 on log-log axes.

# %%
rng5 = np.random.default_rng(5)
R = 400
n_list = np.unique(np.logspace(2, 4.5, 8).astype(int))


def pi_hat(n: int) -> float:
    """4 * (fraction of n uniform darts in the unit square that land inside the quarter disk)."""
    xy = rng5.random((n, 2))
    return 4.0 * np.mean(np.sum(xy**2, axis=1) < 1)


est = np.array([[pi_hat(n) for _ in range(R)] for n in n_list])      # shape (len(n_list), R)
mean_est, sd_est = est.mean(axis=1), est.std(axis=1, ddof=1)
sd_theory = np.sqrt(np.pi * (4 - np.pi) / n_list)
slope_mc = ph.loglog_fit(n_list, sd_est)[0]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
ax1.errorbar(n_list, mean_est, yerr=sd_est, fmt="o-", capsize=3, label="mean $\\pm$ 1 SD of 400 runs")
ax1.axhline(np.pi, color="k", lw=0.8, ls=":")
ax1.set(xscale="log", xlabel="darts per run $n$", ylabel="estimate of $\\pi$ (dimensionless)")
ax1.legend(loc="upper right")
ax2.errorbar(n_list, sd_est, yerr=sd_est / np.sqrt(2 * (R - 1)), fmt="o", capsize=3, label="SD across runs ($\\pm$ 1 SE)")
ax2.plot(n_list, sd_theory, "-", lw=0.8, label="$\\sqrt{\\pi(4-\\pi)/n}$")
ax2.set(xscale="log", yscale="log", xlabel="darts per run $n$", ylabel="standard deviation of the estimate")
ph.add_slope_guide(ax2, -0.5, x0=n_list[2], y0=sd_est[2] * 0.35)
ax2.legend(loc="upper right")
ph.savefig(fig, "mc_errorbars")
print(f"Part 5: fitted log-log slope {slope_mc:.3f} (expected -0.5), SD*sqrt(n) at each n: {np.round(sd_est * np.sqrt(n_list), 3)}")
assert abs(slope_mc + 0.5) < 0.05
assert np.all(np.abs(mean_est - np.pi) < 4 * sd_est / np.sqrt(R))                # means within 4 SE of pi

# What the bar means: one run of n = 1000 darts, per-dart values v_i = 4 * 1[inside].
dart = 4.0 * (np.sum(rng5.random((1000, 2)) ** 2, axis=1) < 1)
sd1, se1 = dart.std(ddof=1), dart.std(ddof=1) / np.sqrt(1000)
ci1 = stats.t.ppf(0.975, 999) * se1
print(f"Part 5: one run, n=1000: SD of a dart {sd1:.3f}, SE of the mean {se1:.4f}, 95% CI half-width {ci1:.4f}; "
      f"SD across 400 runs at n={n_list[3]}: {sd_est[3]:.4f} (SE predicted for that n: {sd_theory[3]:.4f})")
assert math.isclose(se1, sd1 / math.sqrt(1000)) and ci1 > 1.95 * se1

# %% [markdown]
# ## Part 6: A 2D potential: imshow, pcolormesh, contours, and the pixel test (notes Example 4.5; EM-2)
#
# Tilted-dipole potential (units of q/(4 pi eps0) per length) from charges +q at (+0.5, +0.3)
# and -q at (-0.5, -0.3), softened by 0.05. `V[i, j]` belongs to `y[i], x[j]`. Grid points are
# pixel *centres*, so imshow's `extent` must reach half a cell beyond them. The tilt matters:
# a dipole on the x axis is mirror-symmetric in y and would hide a vertically flipped image.

# %%
x = np.linspace(-2, 2, 81)
y = np.linspace(-1.5, 1.5, 61)
X, Y = np.meshgrid(x, y)  # default indexing="xy": X.shape == (len(y), len(x))
V = (1 / np.sqrt((X - 0.5) ** 2 + (Y - 0.3) ** 2 + 0.05**2)
     - 1 / np.sqrt((X + 0.5) ** 2 + (Y + 0.3) ** 2 + 0.05**2))
assert V.shape == (61, 81) and np.allclose(V, -V[::-1, ::-1])                   # V(-x, -y) = -V(x, y)
dx, dy = x[1] - x[0], y[1] - y[0]
extent = (x[0] - dx / 2, x[-1] + dx / 2, y[0] - dy / 2, y[-1] + dy / 2)         # outer EDGES of the pixels
norm = Normalize(-4, 4)
cmap = plt.get_cmap("RdBu_r")

def draw_potential(ax, how: str):
    """Draw V on `ax` three ways; returns the mappable. `how` is imshow, pcolormesh or flipped."""
    ax.grid(False)  # the house-style grid would be drawn over the image and change the pixels
    if how == "imshow":
        return ax.imshow(V, origin="lower", extent=extent, cmap=cmap, norm=norm, aspect="equal")
    if how == "pcolormesh":
        return ax.pcolormesh(X, Y, V, shading="nearest", cmap=cmap, norm=norm)
    return ax.imshow(V, extent=extent, cmap=cmap, norm=norm, aspect="equal")  # origin left at 'upper': the classic mistake


# The one pixel test, on clean figures (contour lines would paint over the pixels we sample):
# the rendered colour at a grid point must be the colour map applied to V at that point.
rng6 = np.random.default_rng(6)
rows, cols = rng6.integers(5, 56, 25), rng6.integers(5, 76, 25)          # 25 random grid points (row i is y[i], column j is x[j])
worst = {}
for how in ("imshow", "pcolormesh", "flipped"):
    fig, ax = plt.subplots(figsize=(7, 5))
    draw_potential(ax, how)
    ax.set(xlim=(x[0], x[-1]), ylim=(y[0], y[-1]))
    worst[how] = np.abs(ph.rgba_at(ax, x[cols], y[rows]) - cmap(norm(V[rows, cols]))).max()
    plt.close(fig)
assert worst["imshow"] < 5e-3 and worst["pcolormesh"] < 5e-3, worst                    # within about one colour-map step, 1/255
assert worst["flipped"] > 0.1                                                  # the flipped image fails the test
print(f"Part 6: pixel test, worst colour error over 25 grid points: imshow {worst['imshow']:.1e}, "
      f"pcolormesh {worst['pcolormesh']:.1e}, flipped imshow {worst['flipped']:.2f}")

fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), sharex=True, sharey=True)
titles = ["imshow(origin='lower', extent=edges)", "pcolormesh(X, Y, V, shading='nearest')", "imshow(extent=...) without origin='lower'"]
for ax, how, title in zip(axs, ("imshow", "pcolormesh", "flipped"), titles):
    im = draw_potential(ax, how)
    ax.set(title=title, xlabel="$x$ (length unit $\\ell$)")
    if how != "flipped":
        cs = ax.contour(X, Y, V, levels=[-3, -2, -1, -0.5, 0.5, 1, 2, 3], colors="k", linewidths=0.6)
        ax.clabel(cs, fmt="%g", fontsize=7)
        zero = ax.contour(X, Y, V, levels=[0], colors="k", linewidths=1.4, linestyles="--")
    ax.set_aspect("equal")
axs[0].set_ylabel("$y$ (length unit $\\ell$)")
cb = fig.colorbar(im, ax=axs, extend="both", label="potential $V$ (units of $q / 4\\pi\\epsilon_0 \\ell$)")
ph.savefig(fig, "potential_heatmap", close=False)
# colorbar orientation: values increase upward, ends are the Normalize limits, top colour is cmap(1)
assert cb.ax.get_ylim() == (-4.0, 4.0) and not cb.ax.yaxis_inverted()
assert np.abs(ph.rgba_at(cb.ax, 0.5, 3.9) - cmap(norm(3.9))).max() < 1e-2
assert np.abs(ph.rgba_at(cb.ax, 0.5, -3.9) - cmap(norm(-3.9))).max() < 1e-2
# the zero contour is the perpendicular bisector of the charges, 0.5 x + 0.3 y = 0 (marching squares is O(dx^2) accurate)
seg = np.vstack([s_ for s_ in zero.allsegs[0] if len(s_)])
dev = np.abs(0.5 * seg[:, 0] + 0.3 * seg[:, 1]).max()
assert dev < 2e-3, dev
print(f"Part 6: zero contour deviates from 0.5 x + 0.3 y = 0 by at most {dev:.1e}")
plt.close(fig)

# half-pixel error: on a coarse grid, an extent through the pixel centres shifts every cell
xc, yc = np.linspace(-2, 2, 9), np.linspace(-1.5, 1.5, 7)
Xc, Yc = np.meshgrid(xc, yc)
Vc = np.sign(Xc + 0.1) * (Yc + 2.0)
for label, ext in [("edges", (-2.25, 2.25, -1.75, 1.75)), ("centres", (-2, 2, -1.5, 1.5))]:
    fig, ax = plt.subplots()
    ax.imshow(Vc, origin="lower", extent=ext, cmap=cmap, norm=Normalize(-4, 4))
    ax.grid(False)
    ax.set(xlim=(-2.25, 2.25), ylim=(-1.75, 1.75))
    ii, jj = np.divmod(np.arange(63), 9)                                  # every (row, column) pair
    bad = int(np.sum(np.abs(ph.rgba_at(ax, xc[jj], yc[ii]) - cmap(Normalize(-4, 4)(Vc[ii, jj]))).max(axis=1) > 1e-2))
    plt.close(fig)
    print(f"Part 6: coarse 9x7 grid, extent through {label}: {bad} of 63 pixel-centre tests fail")
    assert (bad == 0) == (label == "edges")

# %% [markdown]
# ## Part 7: Small multiples: damped oscillators (notes Build 6)
#
# $\ddot x + 2\zeta\dot x + x = 0$, $x(0) = 1$, $\dot x(0) = 0$, time in units of $1/\omega_0$.
# Six phase planes, shared axes, one colour scale (time) and so one colorbar.

# %%
zetas = [0.05, 0.2, 0.5, 1.0, 1.5, 3.0]


def analytic_x(z, t):
    """x(t) for x'' + 2 z x' + x = 0, x(0)=1, x'(0)=0 (under-, critically and over-damped)."""
    if z < 1:
        wd = math.sqrt(1 - z * z)
        return np.exp(-z * t) * (np.cos(wd * t) + z / wd * np.sin(wd * t))
    if z == 1:
        return (1 + t) * np.exp(-t)
    r1, r2 = -z + math.sqrt(z * z - 1), -z - math.sqrt(z * z - 1)
    return (-r2 * np.exp(r1 * t) + r1 * np.exp(r2 * t)) / (r1 - r2)


T_END = 30.0
tt = np.linspace(0, T_END, 600)
tnorm = Normalize(0, T_END)
fig, axs = plt.subplots(2, 3, figsize=(10, 6), sharex=True, sharey=True)
for ax, z in zip(axs.flat, zetas):
    sol = solve_ivp(lambda t, s: [s[1], -2 * z * s[1] - s[0]], (0, T_END), [1.0, 0.0], t_eval=tt, rtol=1e-10, atol=1e-12)
    xx, vv = sol.y
    assert np.max(np.abs(xx - analytic_x(z, tt))) < 1e-7                      # numerics agree with the closed form
    assert np.all(np.diff(0.5 * (xx**2 + vv**2)) <= 1e-12)                     # E decays: dE/dt = -2 z v^2 <= 0
    assert (np.any(xx < 0)) == (z < 1)                                         # only under-damped oscillators cross x = 0
    sc = ax.scatter(xx, vv, c=tt, cmap="viridis", norm=tnorm, s=6)
    ax.set_title(f"$\\zeta = {z:g}$")
    ax.set_box_aspect(1)
for ax in axs[1]:
    ax.set_xlabel("$x / x(0)$ (dimensionless)")
for ax in axs[:, 0]:
    ax.set_ylabel("$\\dot x\\,/\\,(\\omega_0 x(0))$ (dimensionless)")
fig.colorbar(sc, ax=axs, shrink=0.8, label="time $\\omega_0 t$ (rad)")
assert axs[0, 0].get_shared_x_axes().joined(axs[0, 0], axs[1, 2]) and axs[0, 0].get_shared_y_axes().joined(axs[0, 0], axs[1, 2])
ph.savefig(fig, "small_multiples")
print("Part 7: all six panels agree with the closed forms; energy decays; sign changes only for zeta < 1")

# %% [markdown]
# ## Part 8: jet against viridis (notes Build 4)
#
# Lightness L* (CIE 1976) of each colour map along its length. A map whose L* is not monotone
# draws ridges and edges into smooth data. L* = 116 f(Y) - 16 with Y the relative luminance.


# %%
def lightness(cmap_name: str, n: int = 256) -> np.ndarray:
    """CIE L* of `n` evenly spaced colours of a Matplotlib colour map (sRGB, D65 white)."""
    rgb = plt.get_cmap(cmap_name)(np.linspace(0, 1, n))[:, :3]
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    Y = lin @ np.array([0.2126, 0.7152, 0.0722])
    f = np.where(Y > 216 / 24389, np.cbrt(Y), (24389 / 27 * Y + 16) / 116)
    return 116 * f - 16


def n_ridges(a: np.ndarray, prominence: float = 2.0) -> int:
    """Number of peaks of a sampled curve that rise at least `prominence` above their surroundings.

    `find_peaks` copes with the flat tops a 256-entry colour map produces; the 2 L* threshold
    is a convention here (a small step), not a measured perceptual limit.
    """
    return len(find_peaks(a, prominence=prominence)[0])


L_vir, L_jet = lightness("viridis"), lightness("jet")
assert np.all(np.diff(L_vir) > 0)                                       # viridis lightness rises monotonically
assert n_ridges(L_jet) == 1 and not np.all(np.diff(L_jet) > 0) and L_jet[-1] < L_jet.max() - 60
print(f"Part 8: L* at 0, 1/4, 1/2, 3/4, 1 along the map: viridis {np.round(L_vir[[0, 64, 128, 192, 255]], 0)}, "
      f"jet {np.round(L_jet[[0, 64, 128, 192, 255]], 0)}; peaks of L*: viridis {n_ridges(L_vir)}, jet {n_ridges(L_jet)}")

g = np.linspace(-3, 3, 241)
GX, GY = np.meshgrid(g, g)
bump = np.exp(-(GX**2 + GY**2) / 2)                                    # smooth, one maximum at the centre
profile = bump[120]                                                    # the line y = 0
fig, axs = plt.subplots(2, 2, figsize=(9, 7))
for ax, name in zip(axs[0], ["jet", "viridis"]):
    im = ax.imshow(bump, origin="lower", extent=(-3.025, 3.025, -3.025, 3.025), cmap=name, vmin=0, vmax=1)
    ax.axhline(0, color="w", lw=0.8, ls="--")
    ax.grid(False)
    ax.set(title=f"{name}", xlabel="$x$ (units of $\\sigma$)", ylabel="$y$ (units of $\\sigma$)")
    fig.colorbar(im, ax=ax, label="$f$")
axs[1, 0].plot(g, profile, "k")
axs[1, 0].set(title="the data along y = 0", xlabel="$x$ (units of $\\sigma$)", ylabel="$f(x, 0)$")
for name, c in [("jet", ph.OKABE_ITO[1]), ("viridis", ph.OKABE_ITO[0])]:
    axs[1, 1].plot(g, np.interp(profile, np.linspace(0, 1, 256), lightness(name)), color=c, label=name)
axs[1, 1].set(title="perceived lightness along y = 0", xlabel="$x$ (units of $\\sigma$)", ylabel="CIE $L^*$ (0 to 100)")
axs[1, 1].legend()
ph.savefig(fig, "jet_vs_viridis")
Lp_jet = np.interp(profile, np.linspace(0, 1, 256), L_jet)
Lp_vir = np.interp(profile, np.linspace(0, 1, 256), L_vir)
# f has ONE maximum along the line; viridis shows one, jet shows two ridges and a dark centre the data do not have
assert n_ridges(profile, 0.01) == 1 and n_ridges(Lp_vir) == 1 and n_ridges(Lp_jet) == 2
print(f"Part 8: peaks along y = 0: data {n_ridges(profile, 0.01)}, viridis L* {n_ridges(Lp_vir)}, jet L* {n_ridges(Lp_jet)}; "
      f"jet L* at the centre {Lp_jet[120]:.0f} vs {Lp_jet.max():.0f} on the ridges")

# %% [markdown]
# ## Part 9: Pendulum animation (notes Build 7)
#
# `update(i)` sets the artists for frame i; we call it ourselves to check the *data*, then
# let `FuncAnimation` + `PillowWriter` write a GIF. Without a Pillow writer we save a strip of frames.

# %%
th0 = 2.5
sol = solve_ivp(lambda t, s: [s[1], -np.sin(s[0])], (0, 10), [th0, 0.0], dense_output=True, method="DOP853", rtol=1e-11, atol=1e-12)
n_frames = 60
t_frames = np.linspace(0, 10, n_frames, endpoint=False)
theta_f = sol.sol(t_frames)[0]

fig, (axp, axt) = plt.subplots(1, 2, figsize=(7, 3.2), gridspec_kw={"width_ratios": [1, 1.6]})
axp.set(xlim=(-1.2, 1.2), ylim=(-1.2, 0.3), xlabel="$x$ (units of $L$)", ylabel="$y$ (units of $L$)")
axp.set_aspect("equal")
(rod,) = axp.plot([], [], "-o", color="0.2", ms=9, mfc=ph.OKABE_ITO[1])
axt.plot(np.linspace(0, 10, 400), sol.sol(np.linspace(0, 10, 400))[0], color="0.8")
axt.set(xlim=(0, 10), ylim=(-3, 3), xlabel="time $t$ ($\\sqrt{L/g}$)", ylabel="angle $\\theta$ (rad)")
(dot,) = axt.plot([], [], "o", color=ph.OKABE_ITO[1])


def update(i: int):
    """Move the bob and the marker to frame i; return the changed artists."""
    th = theta_f[i]
    rod.set_data([0, math.sin(th)], [0, -math.cos(th)])
    dot.set_data([t_frames[i]], [th])
    return rod, dot


for i in range(n_frames):                                    # data checks, not pixels
    update(i)
    xs_, ys_ = rod.get_data()
    assert math.isclose(math.hypot(xs_[1], ys_[1]), 1.0, rel_tol=1e-12)            # rod length constant
    assert math.isclose(math.atan2(xs_[1], -ys_[1]), theta_f[i], abs_tol=1e-12)    # bob at angle theta_i
gif = FIG / "pendulum.gif"
if animation.writers.is_available("pillow"):
    anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=50, blit=False)
    anim.save(gif, writer=animation.PillowWriter(fps=20), dpi=70)
    assert gif.read_bytes()[:6] in (b"GIF87a", b"GIF89a") and gif.stat().st_size > 5_000
    print(f"Part 9: wrote {gif.name} ({gif.stat().st_size / 1e3:.0f} kB, {n_frames} frames)")
else:  # no writer: a strip of six frames, drawn from the same update() logic
    strip, axs_s = plt.subplots(1, 6, figsize=(12, 2.4), sharey=True)
    for a, i in zip(axs_s, np.linspace(0, n_frames - 1, 6).astype(int)):
        a.plot([0, math.sin(theta_f[i])], [0, -math.cos(theta_f[i])], "-o", color="0.2")
        a.set(xlim=(-1.2, 1.2), ylim=(-1.2, 0.3), aspect="equal", title=f"t = {t_frames[i]:.1f}", xlabel="$x$ (units of $L$)")
    axs_s[0].set_ylabel("$y$ (units of $L$)")
    print("Part 9: no Pillow writer; saved", ph.savefig(strip, "pendulum_strip").name)
plt.close(fig)

# %% [markdown]
# ## Part 10: A truncated axis and its lie factor (notes Example 4.6)
#
# Tufte's *lie factor* is (size of the effect shown in the graphic) / (size of the effect in
# the data). Bars encode value by *length*, so their baseline must be zero.

# %%
labels, means, se = ["alloy A", "alloy B"], np.array([412.0, 421.0]), np.array([3.0, 3.0])
fig, axs = plt.subplots(2, 1, figsize=(5, 6), sharex=True)
for ax, ylim in zip(axs, [(400, 425), (0, 450)]):
    ax.bar(labels, means, yerr=se, capsize=5, color=[ph.OKABE_ITO[0], ph.OKABE_ITO[1]], width=0.5)
    ax.set(ylim=ylim, ylabel="tensile strength (MPa)")
axs[0].set_title("truncated axis (starts at 400 MPa)")
axs[1].set_title("zero baseline")
ph.savefig(fig, "truncated_vs_zero")
shown = (means[1] - 400) / (means[0] - 400) - 1          # relative bar-length difference on the truncated axis
actual = means[1] / means[0] - 1
lie = shown / actual
assert math.isclose(shown, 0.75) and math.isclose(actual, 9 / 412) and 30 < lie < 40
print(f"Part 10: B is {actual:.1%} higher than A; the truncated bars show {shown:.0%}; lie factor {lie:.0f}")
z = (means[1] - means[0]) / math.hypot(*se)
print(f"         difference / its standard error = {z:.2f}")

# %% [markdown]
# ## Part 11: Numbers quoted in section 3 of the notes
#
# A log axis and a `LogNorm` are one map (3.1, 3.2); zeros vanish from log axes; SD, SE and the
# CI differ (3.9); a line's slope on the page depends on the axes' shape; PNG size is figsize x dpi.

# %%
fig, ax = plt.subplots()
ax.set(xlim=(0, 10), ylim=(1, 1000), yscale="log")
fig.canvas.draw()
frac = ax.transAxes.inverted().transform(ax.transData.transform([(5, 10), (5, 100)]))[:, 1]
assert np.allclose(frac, [1 / 3, 2 / 3]) and math.isclose(float(LogNorm(1, 1000)(10)), 1 / 3)
plt.close(fig)

yb = np.array([1e-1, 1e-2, 0.0, -1e-3, 1e-4])                       # two points <= 0 on a log axis
fig, ax = plt.subplots()
(ln,) = ax.semilogy(np.arange(5), yb, "o-")
fig.canvas.draw()
disp_y = ax.transData.transform(ln.get_xydata())[:, 1]
assert np.all(disp_y[[2, 3]] < -1e4) and np.all(disp_y[[0, 1, 4]] > 0)         # off the canvas, no warning
plt.close(fig)

rng3 = np.random.default_rng(3)
for n in (5, 50, 500):                                                # SD describes the data; SE and CI shrink like 1/sqrt(n)
    sample_n = rng3.normal(10, 2, n)
    sd_n = sample_n.std(ddof=1)
    se_n = sd_n / np.sqrt(n)
    half = stats.t.ppf(0.975, n - 1) * se_n
    assert math.isclose(se_n, stats.sem(sample_n)) and half > 1.96 * se_n
    print(f"Part 11: n={n:3d}  mean {sample_n.mean():5.2f}  SD {sd_n:5.2f}  SE {se_n:5.2f}  95% CI half-width {half:5.2f}")

angles = []
for w, hgt in [(4, 4), (8, 2), (2, 8)]:                               # axes of w x h inches, both data ranges 0 to 10
    fig = plt.figure(figsize=(w + 1, hgt + 1), layout="none")
    ax = fig.add_axes((0.5 / (w + 1), 0.5 / (hgt + 1), w / (w + 1), hgt / (hgt + 1)))
    ax.set(xlim=(0, 10), ylim=(0, 10))
    (x0, y0), (x1, y1) = ax.transData.transform([(0, 0), (10, 10)])
    angles.append(np.degrees(np.arctan2(y1 - y0, x1 - x0)))
    plt.close(fig)
assert np.allclose(angles, [45.0, np.degrees(np.arctan(2 / 8)), np.degrees(np.arctan(8 / 2))])
print("Part 11: the line y = x is drawn at", np.round(angles, 1), "degrees in 4x4, 8x2, 2x8 inch axes")

fig, ax = plt.subplots(figsize=(4, 3))
ax.scatter(*np.random.default_rng(2).normal(size=(2, 2000)), s=2)
ax.set(xlabel="$x$ (arb.)", ylabel="$y$ (arb.)")
path = ph.savefig(fig, "size_check", tight=False)
assert mimg.imread(path).shape == (450, 600, 4)                        # figsize x dpi = (3 x 150, 4 x 150)

# %% [markdown]
# ## Your turn
#
# 1. **Slide and paper versions.** Redraw the convergence plot for a slide (`figsize=(10, 5.6)`,
#    `font.size` 16, at most five tick intervals per axis via `MaxNLocator(5)`) and for a paper (`figsize=(3.4, 2.6)`,
#    saved with `fmt="pdf"` and `pdf.fonttype` 42). Assert the sizes and that the PDF starts with `%PDF`.
# 2. **Inset.** Add `ax.inset_axes((0.55, 0.1, 0.4, 0.4))` to the Part 6 potential, zoomed on the
#    region within 0.3 of the +q charge, and check with `rgba_at` that its pixels match the main panel's.
# 3. **Mislead, then fix.** Make a figure that supports a false claim using only axis limits, aspect
#    ratio or a second y-axis; state the claim in a comment; redraw it so the data cannot support it.
# 4. **PY-6 random walk.** Take the PY-6 ensemble (10^4 walkers, 10^3 steps), plot the MSD on
#    log-log axes with `ph.add_slope_guide(ax, 1, ...)` and fit the slope.
# 5. **Direct labels.** Replace each legend here by `ph.label_end`. Which figures need the legend after all?

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
