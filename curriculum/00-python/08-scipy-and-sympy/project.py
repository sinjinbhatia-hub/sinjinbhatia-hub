# %% [markdown]
# # PY-8 Project: The Pendulum Pipeline, Controlled Error, and the SciPy Cabinet
#
# This script runs the whole SymPy -> lambdify -> solve_ivp pipeline of the notes
# (section 3.3) and checks every number the notes quote. Parts: (1) the pendulum:
# Euler-Lagrange, energy, period from events against (3.6); (2) the tolerance sweep
# and your RK4 (3.4, Example 4.3); (3) stiffness (3.5); (4) the double pendulum and
# chaos (3.6); (5) curve_fit, the covariance (3.10) and its coverage (3.7, Example 4.5);
# (6) the sparse Poisson solve and Theorem 3.1 (3.10); (7) expm against solve_ivp (3.9);
# (8) special, stats, fft, interpolate, quad, lambdify (3.2, 3.8, 3.11, 3.12).
#
# Run headless: `MPLBACKEND=Agg python3 project.py` (about a minute). Figures go to ./figures/.
# The dense n = 10_000 solve in Part 6 needs about 1.6 GB of RAM: set DENSE_BIG = False if
# you have less, and the script extrapolates from n = 4000 instead.

# %%
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from scipy import constants, integrate, interpolate, linalg, optimize, sparse, special, stats, fft
from scipy.sparse import linalg as spla

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T_START = time.perf_counter()
DENSE_BIG = True
G, L = 9.81, 1.0


def require(res, what=""):
    """Return a SciPy result object, or raise if it reports failure (.success or .converged)."""
    ok = getattr(res, "success", getattr(res, "converged", None))
    if not ok:
        raise RuntimeError(f"{what}: {getattr(res, 'message', getattr(res, 'flag', res))}")
    return res


# %% [markdown]
# ## Part 1: the pendulum pipeline (notes 3.3)
#
# Stages: Lagrangian -> Euler-Lagrange (3.3) -> solve for theta'' -> substitute derivatives
# first -> lambdify -> solve_ivp -> energy (3.4) -> period from events against (3.6).

# %%
def pendulum_pipeline():
    """Compile the pendulum from its Lagrangian: returns (EL, hand result, theta'' expression, rhs, energy).

    EL is (3.3) for L = m L^2 thetadot^2 / 2 + m g L cos(theta); rhs(th, om, g, L) is the
    first-order system [omega, theta'']; energy(th, om, m, g, L) is (3.4).
    """
    t = sp.Symbol("t", real=True)
    m, Ls, gs = sp.symbols("m L g", positive=True)
    theta = sp.Function("theta")(t)
    lag = sp.Rational(1, 2) * m * Ls**2 * theta.diff(t) ** 2 + m * gs * Ls * sp.cos(theta)
    el = sp.diff(lag.diff(theta.diff(t)), t) - lag.diff(theta)
    acc = sp.solve(el, theta.diff(t, 2))[0]
    th, om = sp.symbols("theta omega", real=True)
    state = lambda ex: ex.subs(theta.diff(t), om).subs(theta, th)          # derivative first!
    rhs = sp.lambdify((th, om, gs, Ls), [om, state(acc)], modules="numpy")
    energy = sp.lambdify((th, om, m, gs, Ls), state(theta.diff(t) * lag.diff(theta.diff(t)) - lag), modules="numpy")
    hand = m * Ls**2 * theta.diff(t, 2) + m * gs * Ls * sp.sin(theta)
    return el, hand, acc, rhs, energy, (theta, Ls, gs)


EL, HAND, ACC, RHS_NUM, ENERGY_NUM, (THETA, LS, GS) = pendulum_pipeline()
assert sp.simplify(EL - HAND) == 0                                          # Example 4.1
assert sp.simplify(ACC + GS * sp.sin(THETA) / LS) == 0
print("Euler-Lagrange:", EL, "= 0   ->   theta'' =", ACC)

mm = sp.Symbol("m")
k_series = sp.series(sp.elliptic_k(mm), mm, 0, 3).removeO()                 # step 6 of the period derivation
assert sp.simplify(k_series - (sp.pi / 2 + sp.pi * mm / 8 + 9 * sp.pi * mm**2 / 128)) == 0


def pend_rhs(t, y, g=G, length=L):
    """rhs(t, y) for y = (theta, omega), compiled from the Euler-Lagrange equation (3.3)."""
    return RHS_NUM(y[0], y[1], g, length)


def spec_energy(y):
    """Specific energy (3.4) divided by m L^2 for L = 1: omega^2/2 - g cos(theta); y has shape (2, n)."""
    return 0.5 * y[1] ** 2 - G * np.cos(y[0])


def period(th0, length=L, g=G):
    """Exact pendulum period (3.6): 4 sqrt(L/g) K(m) with the PARAMETER m = sin^2(theta0/2)."""
    return 4 * np.sqrt(length / g) * special.ellipk(np.sin(th0 / 2) ** 2)


ys = rng.uniform(-3, 3, size=(5, 2))
assert np.allclose(spec_energy(ys.T), ENERGY_NUM(ys[:, 0], ys[:, 1], 1.0, G, 1.0))
assert np.allclose([pend_rhs(0, y) for y in ys], [[y[1], -G * np.sin(y[0])] for y in ys])

TH0 = 2.0
T_EX = period(TH0)
wrong = 4 * np.sqrt(L / G) * special.ellipk(np.sin(TH0 / 2))                 # modulus instead of parameter
assert abs(wrong / T_EX - 1) > 0.05, "passing k instead of m must be visibly wrong"

N_PER = 20
te = np.linspace(0, N_PER * T_EX, 2001)
sol = require(integrate.solve_ivp(pend_rhs, (0, N_PER * T_EX), [TH0, 0.0], method="DOP853", rtol=1e-10, atol=1e-12, t_eval=te))
e_t = spec_energy(sol.y)
drift = np.abs(e_t / e_t[0] - 1)
print(f"energy drift over {N_PER} periods at rtol=1e-10: {drift.max():.2e}")
assert drift.max() < 1e-8


def crossing(t, y):
    """Event function: theta passes through 0 ..."""
    return y[0]


crossing.direction = -1                                                     # ... moving downward


ev = require(integrate.solve_ivp(pend_rhs, (0, N_PER * T_EX), [TH0, 0.0], method="DOP853", rtol=1e-10, atol=1e-12, events=crossing))
T_ev = np.diff(ev.t_events[0])
print(f"period from {len(T_ev)} events: {T_ev.mean():.9f}, exact (3.6): {T_EX:.9f}, relative error {abs(T_ev.mean() / T_EX - 1):.1e}")
assert len(T_ev) >= N_PER - 2 and abs(T_ev.mean() / T_EX - 1) < 1e-6

T_SMALL = 2 * np.pi * np.sqrt(L / G)
th_small = 0.01                                                             # small-angle limit (3.7)
ev_s = integrate.solve_ivp(pend_rhs, (0, 10 * T_SMALL), [th_small, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, events=crossing)
ratio_s = np.diff(ev_s.t_events[0]).mean() / T_SMALL - 1
print(f"theta0=0.01: T/T_small - 1 = {ratio_s:.3e}, theta0^2/16 = {th_small**2 / 16:.3e}")
assert abs(ratio_s - th_small**2 / 16) < 1e-8

fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
axes[0].plot(sol.t, sol.y[0]); axes[0].set(xlabel="t (s)", ylabel="theta (rad)", title="theta(t), theta0 = 2 rad")
axes[1].plot(sol.y[0], sol.y[1]); axes[1].set(xlabel="theta", ylabel="omega", title="phase portrait")
axes[2].semilogy(sol.t, drift + 1e-18); axes[2].set(xlabel="t (s)", ylabel="relative energy drift", title="DOP853, rtol = 1e-10")
for ax in axes:
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG / "pendulum_pipeline.png", dpi=120)
plt.close(fig)


# %% [markdown]
# ## Part 2: tolerance sweep and RK4 (notes 3.4, Example 4.3)
#
# Energy drift over 100 periods against `rtol` (with atol = rtol), log-log. The notes predict
# slope one. `rk4` is the fixed-step method you were asked to write.

# %%
def rk4(f, y0, t_end, n_steps):
    """Classical RK4 with a fixed step; returns the states at the n_steps+1 grid times, shape (len(y0), n_steps+1)."""
    h = t_end / n_steps
    y = np.array(y0, dtype=float)
    out = np.empty((y.size, n_steps + 1))
    out[:, 0] = y
    for k in range(n_steps):
        tk = k * h
        k1 = np.asarray(f(tk, y))
        k2 = np.asarray(f(tk + h / 2, y + h / 2 * k1))
        k3 = np.asarray(f(tk + h / 2, y + h / 2 * k2))
        k4 = np.asarray(f(tk + h, y + h * k3))
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[:, k + 1] = y
    return out


N100 = 100
t100 = np.linspace(0, N100 * T_EX, 2001)
RTOLS = np.array([1e-3, 1e-5, 1e-7, 1e-9, 1e-11])
sweep = {}
for method in ("RK45", "DOP853"):
    for rtol in RTOLS:
        s = require(integrate.solve_ivp(pend_rhs, (0, N100 * T_EX), [TH0, 0.0], method=method, rtol=rtol, atol=rtol, t_eval=t100))
        e = spec_energy(s.y)
        sweep[method, rtol] = (s.nfev, np.max(np.abs(e / e[0] - 1)))
        print(f"{method:7s} rtol={rtol:.0e} nfev={s.nfev:7d} drift={sweep[method, rtol][1]:.2e}")
slopes = {}
for method in ("RK45", "DOP853"):
    d = np.array([sweep[method, r][1] for r in RTOLS])
    slopes[method] = np.polyfit(np.log10(RTOLS[1:]), np.log10(d[1:]), 1)[0]    # skip rtol=1e-3, outside the linear regime
    print(f"{method}: log-log slope of drift against rtol = {slopes[method]:.2f}")
    assert 0.8 < slopes[method] < 1.2
assert sweep["DOP853", 1e-11][0] < sweep["RK45", 1e-11][0]                   # DOP853 is cheaper at tight tolerance

yr = rk4(pend_rhs, [TH0, 0.0], N100 * T_EX, N100 * 50)
e_rk = spec_energy(yr)
rk4_drift = np.max(np.abs(e_rk / e_rk[0] - 1))
print(f"RK4, h = T/50, {4 * N100 * 50} calls: drift {rk4_drift:.2e}")
assert 1e-3 < rk4_drift < 1e-2

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
for method, mk in (("RK45", "o-"), ("DOP853", "s-")):
    ax1.loglog(RTOLS, [sweep[method, r][1] for r in RTOLS], mk, label=f"{method} (slope {slopes[method]:.2f})")
    ax2.loglog([sweep[method, r][0] for r in RTOLS], [sweep[method, r][1] for r in RTOLS], mk, label=method)
ax1.loglog(RTOLS, RTOLS, "k--", lw=0.8, label="slope 1")
ax2.loglog([4 * N100 * 50], [rk4_drift], "k^", label="your RK4")
ax1.set(xlabel="rtol", ylabel="max relative energy drift", title="100 periods"); ax2.set(xlabel="function calls", title="cost against accuracy")
for ax in (ax1, ax2):
    ax.grid(alpha=0.3, which="both"); ax.legend()
fig.tight_layout()
fig.savefig(FIG / "tolerance_sweep.png", dpi=120)
plt.close(fig)


# %% [markdown]
# ## Part 3: stiffness (notes 3.5, equation 3.9)
#
# RK45's stability limit read off the solver itself with brentq, then nfev on a linear
# stiff problem and on Robertson's chemical kinetics.

# %%
def rk45_multiplier(z):
    """RK45's one-step multiplier R(z): one forced-accept step of y' = z y with h = 1 from y = 1."""
    s = integrate.RK45(lambda t, y: z * y, 0.0, [1.0], 1.0, first_step=1.0, rtol=1e10, atol=1e10)
    s.step()
    return s.y[0]


x_lim = optimize.brentq(lambda x: abs(rk45_multiplier(-x)) - 1, 2.0, 5.0)
R6 = lambda z: 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24 + z**5 / 120 + z**6 / 600     # Dormand-Prince 5 stability polynomial
assert abs(rk45_multiplier(-1.0) - R6(-1.0)) < 1e-12 and abs(x_lim - 3.3066) < 1e-3

LAM = 1000.0
stiff = lambda t, y: -LAM * (y - np.cos(t))
exact_stiff = lambda t: (LAM**2 * np.cos(t) + LAM * np.sin(t) - LAM**2 * np.exp(-LAM * t)) / (1 + LAM**2)
runs = {m_: require(integrate.solve_ivp(stiff, (0, 10), [0.0], method=m_, rtol=1e-6, atol=1e-9)) for m_ in ("RK45", "Radau", "LSODA")}
for m_, s in runs.items():
    assert abs(s.y[0, -1] - exact_stiff(10.0)) < 1e-5
ratio1 = runs["RK45"].nfev / runs["Radau"].nfev
print(f"stability limit {x_lim:.4f}; linear stiff problem: nfev RK45 {runs['RK45'].nfev}, Radau {runs['Radau'].nfev}, ratio {ratio1:.0f}")
assert ratio1 > 10 and np.diff(runs["RK45"].t).mean() < x_lim / LAM


def robertson(t, y):
    """Robertson's kinetics: three species, rate constants 0.04, 1e4, 3e7 (stiff)."""
    a, b, c = y
    return [-0.04 * a + 1e4 * b * c, 0.04 * a - 1e4 * b * c - 3e7 * b * b, 3e7 * b * b]


rob = {m_: require(integrate.solve_ivp(robertson, (0, 10), [1.0, 0.0, 0.0], method=m_, rtol=1e-8, atol=1e-10)) for m_ in ("RK45", "Radau", "LSODA")}
ratio2 = rob["RK45"].nfev / rob["Radau"].nfev
print(f"Robertson on [0, 10]: nfev RK45 {rob['RK45'].nfev}, Radau {rob['Radau'].nfev}, LSODA {rob['LSODA'].nfev}; ratio {ratio2:.0f}")
assert ratio2 > 10 and np.allclose(rob["RK45"].y[:, -1], rob["Radau"].y[:, -1], rtol=1e-4, atol=1e-8)
assert abs(rob["Radau"].y[:, -1].sum() - 1) < 1e-8                           # mass is conserved


# %% [markdown]
# ## Part 4: the double pendulum and chaos (notes 3.6, Builds 2 and 8)

# %%
def double_pendulum():
    """Euler-Lagrange for the double pendulum by SymPy: returns rhs(a, b, wa, wb, m1, m2, l1, l2, g) and E(...) via (3.4)."""
    t = sp.Symbol("t", real=True)
    gs = sp.Symbol("g", positive=True)
    M1, M2, l1, l2 = sp.symbols("m1 m2 l1 l2", positive=True)
    q1, q2 = sp.Function("theta1")(t), sp.Function("theta2")(t)
    x1, y1 = l1 * sp.sin(q1), -l1 * sp.cos(q1)
    x2, y2 = x1 + l2 * sp.sin(q2), y1 - l2 * sp.cos(q2)
    sp2 = lambda xx, yy: xx.diff(t) ** 2 + yy.diff(t) ** 2
    lag = sp.simplify(sp.Rational(1, 2) * M1 * sp2(x1, y1) + sp.Rational(1, 2) * M2 * sp2(x2, y2) - gs * (M1 * y1 + M2 * y2))
    qs = [q1, q2]
    el = [sp.diff(lag.diff(q.diff(t)), t) - lag.diff(q) for q in qs]
    accs = [q.diff(t, 2) for q in qs]
    solved = sp.solve(el, accs, dict=True)[0]
    a1, a2, w1, w2 = sp.symbols("a1 a2 w1 w2", real=True)
    state = lambda ex: ex.subs({q1.diff(t): w1, q2.diff(t): w2}).subs({q1: a1, q2: a2})
    args = (a1, a2, w1, w2, M1, M2, l1, l2, gs)
    rhs = sp.lambdify(args, [state(solved[a]) for a in accs], modules="numpy", cse=True)
    energy = sp.lambdify(args, state(sum(q.diff(t) * lag.diff(q.diff(t)) for q in qs) - lag), modules="numpy", cse=True)
    return rhs, energy


def textbook(a, b, wa, wb, g=G):
    """Published hand-derived accelerations of the double pendulum for m1 = m2 = l1 = l2 = 1."""
    d = a - b
    den = 3 - np.cos(2 * d)
    A1 = (-3 * g * np.sin(a) - g * np.sin(a - 2 * b) - 2 * np.sin(d) * (wb**2 + wa**2 * np.cos(d))) / den
    A2 = 2 * np.sin(d) * (2 * wa**2 + 2 * g * np.cos(a) + wb**2 * np.cos(d)) / den
    return np.array([A1, A2])


t0 = time.perf_counter()
DP_RHS, DP_ENERGY = double_pendulum()
print(f"double pendulum derived and compiled in {time.perf_counter() - t0:.1f} s")
pts = rng.uniform(-3, 3, size=(4, 8))
mismatch = max(np.max(np.abs(np.array(DP_RHS(*p, 1.0, 1.0, 1.0, 1.0, G)) - textbook(*p))) for p in pts.T)
print(f"pipeline vs published formula, max |difference| = {mismatch:.1e}")
assert mismatch < 1e-12


def dp_rhs(t, y):
    return [y[2], y[3], *DP_RHS(y[0], y[1], y[2], y[3], 1.0, 1.0, 1.0, 1.0, G)]


y0 = np.array([2.0, 2.0, 0.0, 0.0])
te2 = np.linspace(0, 30, 3001)
runs2 = [require(integrate.solve_ivp(dp_rhs, (0, 30), y0 + d, method="DOP853", rtol=1e-12, atol=1e-12, t_eval=te2))
         for d in (0.0, np.array([1e-8, 0.0, 0.0, 0.0]))]
E_dp = DP_ENERGY(*runs2[0].y, 1.0, 1.0, 1.0, 1.0, G)
print(f"double pendulum energy drift over 30 s: {np.max(np.abs(E_dp / E_dp[0] - 1)):.1e}")
assert np.max(np.abs(E_dp / E_dp[0] - 1)) < 1e-9
sep = np.linalg.norm(runs2[0].y - runs2[1].y, axis=0)
win = (sep > 1e-7) & (sep < 1e-1)
slope, intercept = np.polyfit(te2[win], np.log(sep[win]), 1)
print(f"separation grows like exp({slope:.2f} t) between 1e-7 and 1e-1 (a finite-time estimate, not a Lyapunov exponent)")
assert 0.3 < slope < 3.0 and sep[-1] > 1.0

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.8))
for r_, lab in zip(runs2, ("start A", "start B (+1e-8 in theta1)")):
    ax1.plot(r_.t, r_.y[1], lw=0.9, label=lab)
ax1.set(xlabel="t (s)", ylabel="theta2 (rad)", title="two nearby starts"); ax1.legend()
ax2.semilogy(te2, sep, label="separation in (theta1, theta2, omega1, omega2)")
ax2.semilogy(te2[win], np.exp(intercept + slope * te2[win]), "k--", label=f"exp({slope:.2f} t)")
ax2.set(xlabel="t (s)", title="exponential divergence, then saturation"); ax2.legend()
for ax in (ax1, ax2):
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG / "double_pendulum.png", dpi=120)
plt.close(fig)

# Build 8: lambdify against subs + evalf, same expression, same value.
x_ = sp.Symbol("x")
expr = sp.series(sp.sin(x_) * sp.cos(x_), x_, 0, 12).removeO() + sp.exp(-x_**2) * sp.sin(3 * x_)
f_expr = sp.lambdify(x_, expr, modules="numpy")
t0 = time.perf_counter()
for _ in range(100):
    v_sym = float(expr.subs(x_, 0.37).evalf())
t_sub = (time.perf_counter() - t0) / 100
t0 = time.perf_counter()
for _ in range(10000):
    v_num = f_expr(0.37)
t_lam = (time.perf_counter() - t0) / 10000
print(f"subs+evalf {t_sub * 1e6:.0f} us, lambdify {t_lam * 1e6:.1f} us: speedup {t_sub / t_lam:.0f}x")
assert abs(v_sym - v_num) < 1e-12 and t_sub / t_lam > 10


# %% [markdown]
# ## Part 5: curve_fit, the covariance (3.10), and its coverage (notes 3.7, Example 4.5)

# %%
LS_M, TH_AMP, SIG_T, G_TRUE = np.linspace(0.2, 2.0, 10), 0.5, 0.01, 9.81
MODELS = {"small-angle": lambda Lm, gg: 2 * np.pi * np.sqrt(Lm / gg), "exact (3.6)": lambda Lm, gg: period(TH_AMP, Lm, gg)}


def fit_g(model, T_meas):
    """Weighted least squares for g with known sigma: returns (g_hat, sigma_g) from popt and sqrt(pcov), absolute_sigma=True."""
    popt, pcov = optimize.curve_fit(model, LS_M, T_meas, p0=[8.0], sigma=np.full(LS_M.size, SIG_T), absolute_sigma=True)
    return popt[0], np.sqrt(pcov[0, 0])


T_meas = period(TH_AMP, LS_M, G_TRUE) + SIG_T * np.random.default_rng(0).standard_normal(LS_M.size)
z = {}
for name, model in MODELS.items():
    g_hat, sg = fit_g(model, T_meas)
    z[name] = (g_hat - G_TRUE) / sg
    print(f"{name:12s} g = {g_hat:.3f} +/- {sg:.3f} ({z[name]:+.1f} sigma)")
assert abs(z["exact (3.6)"]) < 2 and abs(z["small-angle"]) > 5               # the covariance does not know the model is wrong

Lq, gq = sp.symbols("Lq gq", positive=True)                                  # (3.10) with a SymPy Jacobian
dTdg = sp.lambdify((Lq, gq), sp.diff(2 * sp.pi * sp.sqrt(Lq / gq), gq), modules="numpy")
g_small, sg_small = fit_g(MODELS["small-angle"], T_meas)
Jm = (-dTdg(LS_M, g_small) / SIG_T).reshape(-1, 1)
sg_310 = np.sqrt(np.linalg.inv(Jm.T @ Jm)[0, 0])
print(f"(3.10) from the SymPy Jacobian: {sg_310:.6f}, curve_fit: {sg_small:.6f}")
assert abs(sg_310 / sg_small - 1) < 1e-5


def coverage(model, reps=500, seed=1):
    """Fraction of repeated experiments whose fitted g is within 1 and 2 reported sigmas of the truth, plus all z-scores."""
    r = np.random.default_rng(seed)
    zs = []
    for _ in range(reps):
        Tm = period(TH_AMP, LS_M, G_TRUE) + SIG_T * r.standard_normal(LS_M.size)
        g_hat, sg = fit_g(model, Tm)
        zs.append((g_hat - G_TRUE) / sg)
    zs = np.array(zs)
    return (np.abs(zs) < 1).mean(), (np.abs(zs) < 2).mean(), zs


c1, c2, z_exact = coverage(MODELS["exact (3.6)"])
s1, s2, _ = coverage(MODELS["small-angle"])
print(f"coverage, exact model: {c1:.3f} within 1 sigma, {c2:.3f} within 2 (Gaussian: 0.683, 0.954); small-angle model: {s1:.3f}, {s2:.3f}")
assert 0.62 < c1 < 0.74 and 0.92 < c2 < 0.98 and s2 == 0.0



# %% [markdown]
# ## Part 6: sparse Poisson and Theorem 3.1 (notes 3.10, Build 6)

# %%
def lap1d(n):
    """CSC form of (3.12): (1/h^2) tridiag(-1, 2, -1), h = 1/(n+1)."""
    h = 1.0 / (n + 1)
    return sparse.diags([-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1], format="csc") / h**2


def best(fn, repeat=3):
    """Smallest wall time of `repeat` calls of fn(), and the last result."""
    times = []
    for _ in range(repeat):
        t0_ = time.perf_counter()
        out = fn()
        times.append(time.perf_counter() - t0_)
    return min(times), out


T6 = sparse.diags([-np.ones(5), 2 * np.ones(6), -np.ones(5)], [-1, 0, 1], format="csr")
an6 = 2 - 2 * np.cos(np.arange(1, 7) * np.pi / 7)                           # Theorem 3.1 for n = 6
assert np.allclose(np.linalg.eigvalsh(T6.toarray()), an6, atol=1e-13) and T6.nnz == 16
assert list(T6.indptr) == [0, 2, 5, 8, 11, 14, 16]                          # Example 4.6

spla.spsolve(lap1d(10), np.ones(10))                                        # warm-up
poisson = {}
for n in (1000, 4000, 10000):
    h = 1 / (n + 1)
    xg = np.arange(1, n + 1) * h
    f_src = np.pi**2 * np.sin(np.pi * xg)
    A_sp = lap1d(n)
    t_sp, u_sp = best(lambda: spla.spsolve(A_sp, f_src))
    err = np.abs(u_sp - np.sin(np.pi * xg)).max()
    pred = np.pi**2 * h**2 / (4 * np.sin(np.pi * h / 2) ** 2) - 1             # exact discrete error from Theorem 3.1
    t_dn = None
    if n < 10000 or DENSE_BIG:
        A_dn = A_sp.toarray()
        t_dn, u_dn = best(lambda: np.linalg.solve(A_dn, f_src), repeat=1 if n > 3000 else 3)
        assert np.abs(u_sp - u_dn).max() < 1e-9
        del A_dn
    poisson[n] = (t_sp, t_dn, err)
    print(f"n={n:6d}: sparse {t_sp * 1e3:8.2f} ms, dense {'n/a' if t_dn is None else f'{t_dn * 1e3:9.1f} ms'}, error {err:.3e}, predicted {pred:.3e}")
    assert err < 1e-6 and abs(err / pred - 1) < 0.02
if not DENSE_BIG:                                                           # extrapolate O(n^3)
    poisson[10000] = (poisson[10000][0], poisson[4000][1] * (10000 / 4000) ** 3, poisson[10000][2])
ratio_1e3 = poisson[1000][1] / poisson[1000][0]
ratio_1e4 = poisson[10000][1] / poisson[10000][0]
print(f"dense/sparse time: {ratio_1e3:.0f}x at n=1000, {ratio_1e4:.0f}x at n=10000 (reported, machine dependent)")
assert ratio_1e3 > 5 and ratio_1e4 > 50
assert poisson[1000][2] / poisson[10000][2] > 90                              # O(h^2): n x 10 -> error / 100

n = 10000
w4 = np.sort(spla.eigsh(lap1d(n), k=4, sigma=0, which="LM", return_eigenvectors=False))
lam_exact = 4 * np.sin(np.arange(1, 5) * np.pi / (2 * (n + 1))) ** 2 * (n + 1) ** 2
print("eigsh smallest four:", np.round(w4, 4), "max relative error vs Theorem 3.1:", np.max(np.abs(w4 / lam_exact - 1)))
assert np.max(np.abs(w4 / lam_exact - 1)) < 1e-9

# sparse * is not @ for the array family (notes 3.10)
Bm = sparse.csr_matrix(sparse.diags([np.ones(2), 2 * np.ones(3), np.ones(2)], [-1, 0, 1]))
Ba = sparse.csr_array(Bm)
assert np.allclose((Bm * Bm).toarray(), (Bm @ Bm).toarray()) and not np.allclose((Ba * Ba).toarray(), (Ba @ Ba).toarray())

ns = np.array(sorted(poisson))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
ax1.loglog(ns, [poisson[k][0] for k in ns], "o-", label="sparse (spsolve)")
ax1.loglog(ns, [poisson[k][1] for k in ns], "s-", label="dense (np.linalg.solve)")
ax1.set(xlabel="n", ylabel="seconds", title="1D Poisson solve"); ax1.legend()
ax2.loglog(1 / (ns + 1), [poisson[k][2] for k in ns], "o-", label="max error")
ax2.loglog(1 / (ns + 1), np.pi**2 / 12 * (1 / (ns + 1)) ** 2, "k--", label="pi^2 h^2 / 12")
ax2.set(xlabel="h", title="second-order convergence"); ax2.legend()
for ax in (ax1, ax2):
    ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(FIG / "sparse_poisson.png", dpi=120)
plt.close(fig)


# %% [markdown]
# ## Part 7: expm against solve_ivp (notes 3.9, Build 7, equation 3.11)

# %%
w0, zeta = 3.0, 0.1
Aosc = np.array([[0.0, 1.0], [-(w0**2), -2 * zeta * w0]])
x0v = np.array([1.0, 0.0])
tt = np.linspace(0, 10, 201)
x_exp = np.array([linalg.expm(Aosc * ti) @ x0v for ti in tt]).T
x_ivp = require(integrate.solve_ivp(lambda ti, y: Aosc @ y, (0, 10), x0v, method="DOP853", rtol=1e-12, atol=1e-14, t_eval=tt)).y
print(f"expm vs solve_ivp: max |difference| = {np.abs(x_exp - x_ivp).max():.1e}")
assert np.abs(x_exp - x_ivp).max() < 1e-8

Mx = rng.standard_normal((5, 5))
Pm, Lo, Up = linalg.lu(Mx)
Qm, Rm = linalg.qr(Mx)
Um, sv, Vh = linalg.svd(Mx)
assert np.allclose(Pm @ Lo @ Up, Mx) and np.allclose(Qm @ Rm, Mx) and np.allclose((Um * sv) @ Vh, Mx)
rank1 = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])
assert linalg.null_space(rank1).shape == (2, 1) and linalg.orth(rank1).shape == (3, 1)


# %% [markdown]
# ## Part 8: quad, interpolate, special, stats, fft (notes 3.8, 3.11, 3.12, Example 4.4)

# %%
kink = lambda x: abs(x - 0.3)
v0, _, i0 = integrate.quad(kink, 0, 1, full_output=True)[:3]
v1, _, i1 = integrate.quad(kink, 0, 1, points=[0.3], full_output=True)[:3]
assert abs(v0 - 0.29) < 1e-12 and abs(v1 - 0.29) < 1e-12 and i1["neval"] < i0["neval"]
spike = lambda x: np.exp(-(((x - 0.3) / 1e-3) ** 2))
exact_spike = 1e-3 * np.sqrt(np.pi)
lie, lie_err = integrate.quad(spike, 0, 10)
pieces = sum(integrate.quad(spike, a, b)[0] for a, b in ((0, 0.29), (0.29, 0.31), (0.31, 10)))
print(f"quad on a narrow spike: {lie} (abserr {lie_err}), pieces {pieces:.6e}, exact {exact_spike:.6e}")
assert lie == 0.0 and abs(pieces / exact_spike - 1) < 1e-8                   # abserr is an estimate, not a bound
assert abs(integrate.dblquad(lambda y, x: x * y**2, 0, 2, 0, 1)[0] - 2 / 3) < 1e-12   # inner variable first
assert abs(integrate.quad(lambda x: np.exp(-(x**2)), -np.inf, np.inf)[0] - np.sqrt(np.pi)) < 1e-7

runge = lambda x: 1 / (1 + 25 * x**2)
xx = np.linspace(-1, 1, 2001)
xn = np.linspace(-1, 1, 13)
xc = np.cos(np.pi * (2 * np.arange(13) + 1) / 26)
poly_eq = np.polynomial.Polynomial.fit(xn, runge(xn), 12)
poly_ch = np.polynomial.Polynomial.fit(xc, runge(xc), 12)
cs = interpolate.CubicSpline(xn, runge(xn))
errs = [np.abs(f_(xx) - runge(xx)).max() for f_ in (poly_eq, poly_ch, cs)]
print(f"Runge, n=13: equispaced polynomial {errs[0]:.3f}, Chebyshev nodes {errs[1]:.4f}, CubicSpline {errs[2]:.4f}")
assert errs[0] > 1 and errs[1] < 0.1 and errs[2] < 0.01
fig, ax = plt.subplots(figsize=(6, 3.8))
ax.plot(xx, runge(xx), "k", label="f"); ax.plot(xx, poly_eq(xx), label="degree 12, equispaced"); ax.plot(xx, cs(xx), label="CubicSpline"); ax.plot(xn, runge(xn), "ko", ms=3)
ax.set(ylim=(-1, 2), title="Runge's phenomenon, 13 nodes"); ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG / "runge.png", dpi=120)
plt.close(fig)

xg_ = np.linspace(-1, 1, 41)
xs_ = sp.Symbol("x")
for deg in range(7):
    leg = sp.lambdify(xs_, sp.legendre(deg, xs_), modules="numpy")
    assert np.allclose(special.eval_legendre(deg, xg_), leg(xg_ * np.ones_like(xg_)), atol=1e-12)
assert abs(special.sph_harm_y(2, 1, 0.7, 0.3) - complex(sp.N(sp.Ynm(2, 1, 0.7, 0.3).expand(func=True), 15))) < 1e-12
assert abs(special.erfc(5.0) - 1.5374597944280348e-12) < 1e-26 and abs((1 - special.erf(5.0)) / special.erfc(5.0) - 1) > 1e-6
assert abs(special.jn_zeros(0, 1)[0] - 2.404825557695773) < 1e-12 and special.gamma(5) == 24.0

assert abs(stats.norm.ppf(0.975) - 1.959964) < 1e-6 and abs(stats.t.ppf(0.975, 9) - 2.262157) < 1e-6
sample = stats.norm.rvs(loc=1.0, scale=2.0, size=500, random_state=rng)
mu_hat, sd_hat = stats.norm.fit(sample)
assert abs(mu_hat - 1.0) < 3 * 2.0 / np.sqrt(500) and abs(sd_hat**2 * 500 / 499 - stats.describe(sample).variance) < 1e-9
assert abs(constants.g - 9.80665) < 1e-12

N, fs = 64, 64.0
sig = np.sin(2 * np.pi * 5 * np.arange(N) / fs)
F, fr = fft.fft(sig), fft.fftfreq(N, d=1 / fs)
kpk = int(np.argmax(np.abs(F[: N // 2])))
assert kpk == 5 and fr[kpk] == 5.0 and abs(abs(F[5]) - N / 2) < 1e-9 and fr[-1] == -1.0 and fr[N // 2] == -32.0
assert fft.fftshift(fr)[0] == -32.0 and fft.rfftfreq(N, 1 / fs).size == N // 2 + 1
sig2 = np.sin(2 * np.pi * 5.5 * np.arange(N) / fs)                           # between bins: leakage
assert np.abs(fft.fft(sig2)).max() < 0.8 * N / 2

print(f"Total run time so far: {time.perf_counter() - T_START:.1f} s.")


# %% [markdown]
# ## Your turn
#
# 1. **Resonance.** Add damping and driving to the pendulum, theta'' + gamma theta' + (g/L) sin(theta) = F cos(Omega t),
#    by adding a generalized force to the Euler-Lagrange equation in `pendulum_pipeline`. Sweep Omega, use events (omega = 0)
#    to read the steady-state amplitude, and plot the resonance curve. Compare with the linear amplitude
#    F / sqrt((w0^2 - Omega^2)^2 + gamma^2 Omega^2) and watch the peak bend as F grows.
# 2. **Rosenbrock.** `minimize` f(x, y) = (1 - x)^2 + 100 (y - x^2)^2 from (-1.2, 1) with BFGS with and without the analytic
#    gradient (from SymPy, lambdified). Compare `nfev`, `njev`, and the iterations; then try `Nelder-Mead` and `trust-exact`
#    with the SymPy Hessian.
# 3. **A boundary-value problem.** Use `integrate.solve_bvp` for u'' + exp(u) = 0, u(0) = u(1) = 0, from two different initial
#    guesses, and find the two solutions. Check one with `solve_ivp` shooting from the returned u'(0).
# 4. **Closed form against numerics.** `sp.dsolve` the damped oscillator x'' + 2 zeta w0 x' + w0^2 x = 0, x(0) = 1, x'(0) = 0
#    (use `ics=`), lambdify the result, and compare with `solve_ivp` and with Part 7's `expm`.

# %%
if __name__ == "__main__":
    print("All checks passed.")
