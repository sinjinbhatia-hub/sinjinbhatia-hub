# %% [markdown]
# # F6 Project: A Discretized Functional Minimizer
#
# A function y(x) is stored as its values at N interior nodes; the functional
# J[y] is evaluated on the polygon through the nodes (Euler's 1744 picture,
# notes §2 and §3.5); scipy.optimize.minimize does the rest. Parts: (1) arc
# length -> straight line, with the Euler–Lagrange residual tracked along the
# iterations; (2) the brachistochrone -> cycloid, plus a discretization that
# cheats and the case that dips below the endpoint; (3) shooting on the E–L ODE
# and the Beltrami constant as a test; (4) the catenary with a length
# constraint and its Lagrange multiplier; (5) Snell's law from Fermat;
# (6) geodesics on a sphere (two dependent variables); (7) J[y + eps*eta] vs eps;
# (8) gradient of J_h versus h times the E–L expression, eq. (3.9);
# (9) SymPy checks of the identities used in the notes and exercise numbers.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import minimize, minimize_scalar, brentq, approx_fprime, NonlinearConstraint
from scipy.integrate import quad, solve_ivp
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
g = 9.81  # m/s^2


def with_ends(y_int, ya, yb):
    """Interior values -> full node vector including the fixed endpoints."""
    return np.concatenate([[ya], np.asarray(y_int, float), [yb]])


# %% [markdown]
# ## Part 1: Arc length -> straight line (Example 4.1)
#
# J[y] = ∫ sqrt(1 + y'^2) dx. On the polygon through the nodes, the trapezoid
# rule is exact: J_h = Σ sqrt(h^2 + Δy_i^2). Its gradient is the discrete E–L
# expression (3.8) with L_y = 0 and L_y' = p/sqrt(1+p^2) (this is the N-point
# computation of §2). **Build it yourself first:** `polygon_length` and
# `el_residual_arclength`.


# %%
def polygon_length(y_int, x, ya, yb):
    """Trapezoid-rule arc length of the polygon through (x_i, y_i): Σ sqrt(dx^2 + dy^2)."""
    y = with_ends(y_int, ya, yb)
    return float(np.sum(np.hypot(np.diff(x), np.diff(y))))


def polygon_length_grad(y_int, x, ya, yb):
    """∂J_h/∂y_i = p_{i-1}/sqrt(1+p_{i-1}^2) - p_i/sqrt(1+p_i^2): the §2 vertex condition."""
    y = with_ends(y_int, ya, yb)
    dx, dy = np.diff(x), np.diff(y)
    s = dy / np.hypot(dx, dy)          # sin of each segment's angle = p/sqrt(1+p^2)
    return s[:-1] - s[1:]


def el_residual_arclength(y_full, h):
    """E–L expression for L = sqrt(1+y'^2): -y''/(1+y'^2)^{3/2}, by central differences at interior nodes."""
    yp = (y_full[2:] - y_full[:-2]) / (2 * h)
    ypp = (y_full[2:] - 2 * y_full[1:-1] + y_full[:-2]) / h**2
    return -ypp / (1 + yp**2) ** 1.5


N = 60
a_, b_, A_, B_ = 0.0, 1.0, 0.0, 1.0
x = np.linspace(a_, b_, N + 2); h = x[1] - x[0]
y0 = x[1:-1] + 0.3 * np.sin(3 * np.pi * x[1:-1]) + 0.1 * rng.standard_normal(N)   # wiggly start
history = [np.max(np.abs(el_residual_arclength(with_ends(y0, A_, B_), h)))]
res_line = minimize(polygon_length, y0, args=(x, A_, B_), jac=polygon_length_grad, method="L-BFGS-B",
                    callback=lambda yk: history.append(np.max(np.abs(el_residual_arclength(with_ends(yk, A_, B_), h)))),
                    options=dict(maxiter=5000, ftol=1e-15, gtol=1e-12))
y_line_exact = A_ + (B_ - A_) * (x[1:-1] - a_) / (b_ - a_)
dev_line = np.max(np.abs(res_line.x - y_line_exact))
assert dev_line < 1e-3, dev_line
assert np.isclose(res_line.fun, np.hypot(b_ - a_, B_ - A_), rtol=1e-10)
assert history[-1] < 1e-3 * history[0], "E–L residual should collapse along the iterations"
print(f"Part 1: arc-length minimizer -> line, max deviation {dev_line:.1e}; "
      f"E–L residual {history[0]:.2f} -> {history[-1]:.1e} in {len(history) - 1} iterations")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(x, with_ends(y0, A_, B_), ".-", alpha=.5, label="initial guess")
ax1.plot(x, with_ends(res_line.x, A_, B_), "k-", label="minimizer")
ax1.set_title("Arc length: minimizer is the straight line"); ax1.legend(); ax1.grid(alpha=.3)
ax2.semilogy(history, ".-"); ax2.set_xlabel("L-BFGS iteration"); ax2.set_ylabel("max |E–L residual|")
ax2.set_title("Euler–Lagrange residual -> 0 along the iterations"); ax2.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "arclength_minimizer.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 2: The brachistochrone (Example 4.2)
#
# y is measured DOWNWARD. Cycloid (4.1): x = a(θ − sin θ), y = a(1 − cos θ);
# endpoint condition (4.2); time (4.3): T = θ_1 sqrt(a/g).
# Numerics: the exact sliding time along a straight segment with constant
# acceleration is Δs / (mean speed) = 2Δs/(v_i + v_{i+1}), v = sqrt(2 g y).
# **Build it yourself first:** `polygon_descent_time`.


# %%
def cycloid_params(x1, y1):
    """Solve (4.2) for θ_1, then a = y1/(1 − cos θ_1). Returns (θ_1, a)."""
    f = lambda th: (th - np.sin(th)) * y1 - x1 * (1 - np.cos(th))
    th1 = brentq(f, 1e-3, 2 * np.pi - 1e-3)   # (4.2) has exactly one root in (0, 2π); avoid the trivial θ = 0
    return th1, y1 / (1 - np.cos(th1))


def cycloid_y_of_x(xs, a, th_max):
    """Depth y(x) on the cycloid arc θ ∈ [0, th_max] (θ is monotone in x since dx/dθ = a(1−cos θ) ≥ 0)."""
    out = np.empty_like(xs, dtype=float)
    for k, xv in enumerate(xs):
        th = 0.0 if xv <= 0 else brentq(lambda t: a * (t - np.sin(t)) - xv, 0.0, th_max + 1e-9)
        out[k] = a * (1 - np.cos(th))
    return out


def polygon_descent_time(y_int, x, ya, yb):
    """Exact frictionless sliding time along the polygon, from rest at depth 0: Σ 2Δs/(v_i + v_{i+1}), v = sqrt(2 g y)."""
    y = with_ends(y_int, ya, yb)
    ds = np.hypot(np.diff(x), np.diff(y)); v = np.sqrt(2 * g * np.maximum(y, 0.0))
    return float(np.sum(2 * ds / (v[1:] + v[:-1])))


def midpoint_descent_time(y_int, x, ya, yb):
    """Midpoint-rule quadrature of ∫ ds / sqrt(2 g y). Looks reasonable; Part 2b shows it lies."""
    y = with_ends(y_int, ya, yb)
    ds = np.hypot(np.diff(x), np.diff(y)); ym = 0.5 * (y[1:] + y[:-1])
    return float(np.sum(ds / np.sqrt(2 * g * np.maximum(ym, 1e-300))))


# Analytic numbers quoted in Example 4.2.
x1, y1 = 1.0, 1.0
th1, a_cyc = cycloid_params(x1, y1)
T_cyc = th1 * np.sqrt(a_cyc / g)
T_line = 2 / np.sqrt(g)
assert np.isclose(th1, 2.4120, atol=5e-5) and np.isclose(a_cyc, 0.57292, atol=5e-6)
assert np.isclose(T_cyc, 0.5829, atol=5e-5) and np.isclose(T_line, 0.6386, atol=5e-5)
# T_line two ways: incline kinematics, and the functional with y' = 1.
assert np.isclose(quad(lambda s: np.sqrt(2) / np.sqrt(2 * g * s), 0, 1)[0], T_line)
# T_cyc also by direct quadrature of ∫ ds/v along the cycloid (dt = sqrt(a/g) dθ).
assert np.isclose(quad(lambda th: 2 * a_cyc * np.sin(th / 2) / (2 * np.sqrt(g * a_cyc) * np.sin(th / 2)), 0, th1)[0], T_cyc)
end_slope = np.sin(th1) / (1 - np.cos(th1))
assert np.isclose(end_slope, 1 / np.tan(th1 / 2)) and np.isclose(np.degrees(np.arctan(end_slope)), 20.9, atol=0.05)
print(f"Part 2: θ1 = {th1:.4f} rad ({np.degrees(th1):.1f}°), a = {a_cyc:.5f} m, "
      f"T_cyc = {T_cyc:.4f} s, T_line = {T_line:.4f} s, cycloid faster by {100 * (1 - T_cyc / T_line):.1f}%")

# (b) Numerical brachistochrone with the exact polygon time.
N = 100
xb = np.linspace(0, x1, N + 2)
yb0 = xb[1:-1].copy()  # start from the straight line
res_br = minimize(polygon_descent_time, yb0, args=(xb, 0.0, y1), method="L-BFGS-B",
                  bounds=[(1e-9, None)] * N, options=dict(maxiter=5000, maxfun=10**6, ftol=1e-14, gtol=1e-10))
y_cyc_nodes = cycloid_y_of_x(xb[1:-1], a_cyc, th1)
rel_T = (res_br.fun - T_cyc) / T_cyc
assert abs(rel_T) < 5e-3, rel_T
assert rel_T > 0, "a faithful discretization cannot beat the true minimum"
dev_br = np.max(np.abs(res_br.x - y_cyc_nodes))
assert dev_br < 2e-2, dev_br
assert polygon_descent_time(yb0, xb, 0.0, y1) > res_br.fun  # the line is slower
print(f"        numerical T_h = {res_br.fun:.5f} s (analytic {T_cyc:.5f}, rel. error {rel_T:+.2e}); "
      f"max shape deviation {dev_br:.1e} (slow convergence near the vertical start)")

# A discretization that lies: the midpoint rule near the singular start.
Nm = 50
xm = np.linspace(0, x1, Nm + 2)
res_mid = minimize(midpoint_descent_time, xm[1:-1].copy(), args=(xm, 0.0, y1), method="L-BFGS-B",
                   bounds=[(1e-9, None)] * Nm, options=dict(maxiter=5000, maxfun=10**6, ftol=1e-14))
assert res_mid.fun < 0.9 * T_cyc, "midpoint rule should 'beat' the true optimum, which is impossible"
assert res_mid.x[0] > 0.3, "the cheat is a near-vertical dive at the first node"
print(f"        midpoint-rule 'minimum' {res_mid.fun:.4f} s < true minimum {T_cyc:.4f} s: "
      f"the optimizer dives to depth {res_mid.x[0]:.2f} at the first node (x = {xm[1]:.3f}). Quadrature lied.")

# The dip below the endpoint: (x1, y1) = (3, 1) has x1/y1 > π/2.
th_dip, a_dip = cycloid_params(3.0, 1.0)
assert th_dip > np.pi and 2 * a_dip > 1.0
assert np.isclose(th_dip, 4.052, atol=5e-4) and np.isclose(a_dip, 0.6197, atol=5e-5)
# Threshold: θ1 = π exactly when x1/y1 = π/2.
assert np.isclose(cycloid_params(np.pi / 2, 1.0)[0], np.pi)
print(f"        dip case (3, 1): θ1 = {th_dip:.3f} > π, a = {a_dip:.4f}, lowest point 2a = {2 * a_dip:.3f} > 1")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
th = np.linspace(0, th1, 300)
ax1.plot(a_cyc * (th - np.sin(th)), a_cyc * (1 - np.cos(th)), "k-", lw=2, label=f"cycloid, T = {T_cyc:.4f} s")
ax1.plot(xb, with_ends(res_br.x, 0, y1), "C1.", ms=3, label=f"L-BFGS-B polygon, T_h = {res_br.fun:.4f} s")
ax1.plot([0, 1], [0, 1], "C0--", label=f"straight line, T = {T_line:.4f} s")
ax1.plot(xm, with_ends(res_mid.x, 0, y1), "C3:", label=f"midpoint-rule cheat, 'T' = {res_mid.fun:.4f} s")
ax1.invert_yaxis(); ax1.set_aspect("equal"); ax1.legend(fontsize=8); ax1.grid(alpha=.3)
ax1.set_title("Brachistochrone (0,0) -> (1,1), y downward")
th = np.linspace(0, th_dip, 400)
ax2.plot(a_dip * (th - np.sin(th)), a_dip * (1 - np.cos(th)), "k-", lw=2)
ax2.plot([0, 3], [0, 1], "ko"); ax2.axhline(1, color="C0", ls="--", lw=1, label="endpoint depth")
ax2.invert_yaxis(); ax2.set_aspect("equal"); ax2.grid(alpha=.3); ax2.legend()
ax2.set_title(f"(0,0) -> (3,1): the fastest wire dips to depth {2 * a_dip:.2f}")
fig.tight_layout(); fig.savefig(FIG / "brachistochrone.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 3: Shooting on the E–L ODE, and the Beltrami constant as a test
#
# Beltrami (Example 4.2): y (1 + y'^2) = 2a, i.e. dx/dy = sqrt(y/(2a − y)) with
# x(0) = 0. The initial slope dy/dx is infinite, so the shooting parameter is
# the first-integral constant 2a; brentq adjusts it until x(y1) = x1.
# Then H = L − y' L_y' = 1/sqrt(2 g y (1+y'^2)) is evaluated along the cycloid
# (constant) and along the straight line (not constant).


# %%
def x_at_end(two_a, y_end=y1):
    """Integrate dx/dy = sqrt(y/(2a − y)) from y = 0 to y_end; return x(y_end)."""
    sol = solve_ivp(lambda y, xx: [np.sqrt(y / (two_a - y))], (0.0, y_end), [0.0], rtol=1e-11, atol=1e-13)
    return sol.y[0, -1]


two_a_shoot = brentq(lambda s: x_at_end(s) - x1, y1 * (1 + 1e-6), 20.0)
assert abs(two_a_shoot - 2 * a_cyc) < 1e-9, (two_a_shoot, 2 * a_cyc)
print(f"Part 3: shooting gives 2a = {two_a_shoot:.9f}; cycloid 2a = {2 * a_cyc:.9f}")


def beltrami_H(y_nodes, yp_nodes):
    """H = L − y' L_y' for L = sqrt(1+y'^2)/sqrt(2 g y): H = 1/sqrt(2 g y (1+y'^2)), constant on solutions (Example 4.2)."""
    return 1 / np.sqrt(2 * g * y_nodes * (1 + yp_nodes**2))


def central_slope(x_nodes, y_nodes):
    """y' at interior nodes by central differences."""
    return (y_nodes[2:] - y_nodes[:-2]) / (x_nodes[2:] - x_nodes[:-2])


th = np.linspace(0.3, th1, 400)
H_cyc = beltrami_H(a_cyc * (1 - np.cos(th)), np.sin(th) / (1 - np.cos(th)))       # exact slope cot(θ/2)
H_exact = 1 / np.sqrt(2 * g * 2 * a_cyc)
assert np.ptp(H_cyc) / H_cyc.mean() < 1e-12 and np.isclose(H_cyc.mean(), H_exact, rtol=1e-12)
y_num = with_ends(res_br.x, 0.0, y1)                                               # the L-BFGS-B polygon of Part 2
H_num = beltrami_H(y_num[1:-1], central_slope(xb, y_num))[10:]                     # skip the first nodes near the vertical start
assert np.ptp(H_num) / H_num.mean() < 2e-3 and abs(H_num.mean() - H_exact) / H_exact < 2e-3
xs_line = np.linspace(0.01, 1, 400)
H_line = beltrami_H(xs_line[1:-1], central_slope(xs_line, xs_line))                # the straight wire y = x
assert H_line.max() / H_line.min() > 5
print(f"        Beltrami H = 1/sqrt(4ga) = {H_exact:.5f}: along the cycloid spread/mean = {np.ptp(H_cyc) / H_cyc.mean():.1e}; "
      f"along the numerical polygon {np.ptp(H_num) / H_num.mean():.1e} (mean {H_num.mean():.5f}); "
      f"along the line max/min = {H_line.max() / H_line.min():.1f}")

# %% [markdown]
# ## Part 4: The catenary with a length constraint (Example 4.4)
#
# Minimize V = ∫ y ds subject to ∫ ds = ℓ. On the polygon both are exact
# (trapezoid rule): V_h = Σ Δs_i (y_i + y_{i+1})/2, ℓ_h = Σ Δs_i. SLSQP handles
# the equality constraint; trust-constr also returns the multiplier, which
# should equal λ in (4.7).

# %%
ell = 3.0
c_cat = brentq(lambda c: 2 * c * np.sinh(1 / c) - ell, 0.2, 5.0)
lam_cat = c_cat * np.cosh(1 / c_cat)
y_cat = lambda xx: c_cat * np.cosh(xx / c_cat) - lam_cat          # (4.7) with x0 = 0, y(±1) = 0
assert np.isclose(c_cat, 0.61647, atol=5e-6) and np.isclose(lam_cat, 1.6217, atol=5e-5)
assert np.isclose(y_cat(0.0), -1.0053, atol=5e-5)
assert np.isclose(quad(lambda xx: np.sqrt(1 + np.sinh(xx / c_cat) ** 2), -1, 1)[0], ell)   # length check
print(f"Part 4: catenary c = {c_cat:.5f}, λ = {lam_cat:.4f}, lowest point y(0) = {y_cat(0.0):.4f}")


def polygon_potential(y_int, x, ya, yb):
    """∫ y ds on the polygon (trapezoid rule, exact for a straight segment): Σ Δs (y_i + y_{i+1})/2."""
    y = with_ends(y_int, ya, yb)
    return float(np.sum(np.hypot(np.diff(x), np.diff(y)) * 0.5 * (y[1:] + y[:-1])))


Nc = 40
xc = np.linspace(-1, 1, Nc + 2)
yc0 = -0.9 * (1 - xc[1:-1] ** 2)
res_cat = minimize(polygon_potential, yc0, args=(xc, 0.0, 0.0), method="SLSQP",
                   constraints=[dict(type="eq", fun=lambda yy: polygon_length(yy, xc, 0.0, 0.0) - ell)],
                   options=dict(maxiter=1000, ftol=1e-12))
assert res_cat.success and np.isclose(polygon_length(res_cat.x, xc, 0, 0), ell, atol=1e-8)
dev_cat = np.max(np.abs(res_cat.x - y_cat(xc[1:-1])))
assert dev_cat < 1e-3, dev_cat
res_tc = minimize(polygon_potential, yc0, args=(xc, 0.0, 0.0), method="trust-constr",
                  constraints=[NonlinearConstraint(lambda yy: polygon_length(yy, xc, 0.0, 0.0), ell, ell)],
                  options=dict(maxiter=3000, xtol=1e-10, gtol=1e-10))
mult = float(np.ravel(res_tc.v[0])[0])
assert abs(mult - lam_cat) < 2e-3, (mult, lam_cat)
print(f"        SLSQP: max deviation from cosh {dev_cat:.1e}; trust-constr multiplier {mult:.4f} vs λ = {lam_cat:.4f}")

fig, ax = plt.subplots(figsize=(6, 4))
xx = np.linspace(-1, 1, 300)
ax.plot(xx, y_cat(xx), "k-", lw=2, label=r"$c\cosh(x/c)-\lambda$")
ax.plot(xc, with_ends(res_cat.x, 0, 0), "C1.", label="SLSQP with length constraint")
ax.plot(xc, with_ends(yc0, 0, 0), "C0--", alpha=.6, label="initial guess")
ax.axhline(-lam_cat, color="gray", ls=":", label=r"directrix $y=-\lambda$")
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_title("Catenary: 3 m chain, supports 2 m apart")
fig.tight_layout(); fig.savefig(FIG / "catenary.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5: Snell's law from Fermat's principle (Example 4.3)
# Two media; the path is straight in each, so the only unknown is the crossing point x_c.

# %%
n1, n2 = 1.0, 1.5
A_pt, B_pt = (0.0, 1.0), (2.0, -1.0)
optical_path = lambda xc_: n1 * np.hypot(xc_ - A_pt[0], A_pt[1]) + n2 * np.hypot(B_pt[0] - xc_, B_pt[1])
r = minimize_scalar(optical_path, bounds=(A_pt[0], B_pt[0]), method="bounded", options=dict(xatol=1e-12))
xcross = r.x
sin1 = (xcross - A_pt[0]) / np.hypot(xcross - A_pt[0], A_pt[1])
sin2 = (B_pt[0] - xcross) / np.hypot(B_pt[0] - xcross, B_pt[1])
assert np.isclose(n1 * sin1, n2 * sin2, atol=1e-7)
assert np.isclose(xcross, 1.3627, atol=5e-5) and np.isclose(np.degrees(np.arcsin(sin1)), 53.73, atol=0.01)
assert np.isclose(np.degrees(np.arcsin(sin2)), 32.51, atol=0.01) and np.isclose(n1 * sin1, 0.8062, atol=5e-5)
print(f"Part 5: Snell: crossing x_c = {xcross:.4f}, θ1 = {np.degrees(np.arcsin(sin1)):.2f}°, "
      f"θ2 = {np.degrees(np.arcsin(sin2)):.2f}°, n1 sin θ1 = {n1 * sin1:.4f} = n2 sin θ2 = {n2 * sin2:.4f}")

# %% [markdown]
# ## Part 6: Geodesics on the unit sphere, two dependent variables (§3.7)
#
# Length is invariant under reparametrization, so minimizing it over (θ(t), φ(t))
# is degenerate. Minimize the energy E = ∫ (θ'^2 + sin^2θ φ'^2) dt instead: its
# minimizer is the great circle traversed at constant speed (so E_min = Ω^2 with
# Ω the central angle), which is exactly spherical linear interpolation (slerp).

# %%
def sph(th_, ph_):
    return np.stack([np.sin(th_) * np.cos(ph_), np.sin(th_) * np.sin(ph_), np.cos(th_)], axis=-1)


th_a, ph_a, th_b, ph_b = np.pi / 3, 0.0, np.pi / 4, np.pi / 2
P0, P1 = sph(th_a, ph_a), sph(th_b, ph_b)
Omega = np.arccos(P0 @ P1)


def energy_and_grad(z, N, h):
    """E_h = (1/h) Σ [Δθ_i^2 + sin^2(θ̄_i) Δφ_i^2] (midpoint θ̄) and its exact gradient w.r.t. interior (θ, φ)."""
    th_ = np.concatenate([[th_a], z[:N], [th_b]]); ph_ = np.concatenate([[ph_a], z[N:], [ph_b]])
    dth, dph = np.diff(th_), np.diff(ph_); thm = 0.5 * (th_[1:] + th_[:-1]); s2 = np.sin(thm) ** 2
    E = np.sum(dth**2 + s2 * dph**2) / h
    sc = np.sin(thm) * np.cos(thm) * dph**2
    g_th = (-2 * dth[1:] + 2 * dth[:-1] + sc[1:] + sc[:-1]) / h
    g_ph = (-2 * s2[1:] * dph[1:] + 2 * s2[:-1] * dph[:-1]) / h
    return E, np.concatenate([g_th, g_ph])


Ns = 60; hs = 1 / (Ns + 1); ts = np.linspace(0, 1, Ns + 2)
z0 = np.concatenate([np.linspace(th_a, th_b, Ns + 2)[1:-1], np.linspace(ph_a, ph_b, Ns + 2)[1:-1]])
assert np.max(np.abs(approx_fprime(z0, lambda z: energy_and_grad(z, Ns, hs)[0], 1e-7) - energy_and_grad(z0, Ns, hs)[1])) < 1e-4
res_geo = minimize(energy_and_grad, z0, args=(Ns, hs), jac=True, method="L-BFGS-B",
                   options=dict(maxiter=10000, ftol=1e-15, gtol=1e-12, maxcor=30))
th_n = np.concatenate([[th_a], res_geo.x[:Ns], [th_b]]); ph_n = np.concatenate([[ph_a], res_geo.x[Ns:], [ph_b]])
P_num = sph(th_n, ph_n)
P_slerp = (np.sin((1 - ts)[:, None] * Omega) * P0 + np.sin(ts[:, None] * Omega) * P1) / np.sin(Omega)
ang_dev = np.arccos(np.clip(np.sum(P_num * P_slerp, axis=1), -1, 1))
len_num = np.sum(np.linalg.norm(np.diff(P_num, axis=0), axis=1))
assert ang_dev.max() < 1e-4 and abs(len_num - Omega) < 1e-4 and abs(res_geo.fun - Omega**2) < 1e-3
print(f"Part 6: sphere geodesic: max angular deviation from the great circle {ang_dev.max():.1e} rad, "
      f"length {len_num:.6f} vs central angle {Omega:.6f}, E_h = {res_geo.fun:.6f} vs Ω^2 = {Omega**2:.6f}")

fig = plt.figure(figsize=(5, 5)); ax = fig.add_subplot(111, projection="3d")
u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 30), np.linspace(0, np.pi, 15))
ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color="lightgray", lw=.5)
ax.plot(*P_slerp.T, "k-", lw=2, label="great circle (slerp)")
ax.plot(*P_num.T, "C1.", ms=4, label="energy minimizer")
ax.plot(*sph(np.linspace(th_a, th_b, 50), np.linspace(ph_a, ph_b, 50)).T, "C0--", label="initial guess")
ax.legend(fontsize=8); ax.set_title("Geodesic on the sphere"); ax.view_init(25, 35)
fig.tight_layout(); fig.savefig(FIG / "sphere_geodesic.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: J[y + ε η] versus ε — what "stationary" looks like (§2)
# Arc length on [0, 1], η = sin(πx). True solution y = x: flat at ε = 0. Non-solution y = x^2: sloped.

# %%
eta_p = lambda xx: np.pi * np.cos(np.pi * xx)
J_eps = lambda yp, eps: quad(lambda xx: np.sqrt(1 + (yp(xx) + eps * eta_p(xx)) ** 2), 0, 1)[0]
slope_true = quad(lambda xx: 1.0 * eta_p(xx) / np.sqrt(2.0), 0, 1)[0]                 # (3.3) with L_y = 0
slope_wrong = quad(lambda xx: 2 * xx * eta_p(xx) / np.sqrt(1 + 4 * xx**2), 0, 1)[0]
assert abs(slope_true) < 1e-12 and np.isclose(slope_wrong, -0.5277, atol=5e-5)
eps = np.linspace(-0.6, 0.6, 61)
Jt = np.array([J_eps(lambda xx: 1.0 + 0 * xx, e) for e in eps]); Jw = np.array([J_eps(lambda xx: 2 * xx, e) for e in eps])
# Finite-difference slopes at ε = 0 agree with the first variation.
assert abs((Jt[31] - Jt[29]) / (eps[31] - eps[29]) - slope_true) < 1e-6
assert abs((Jw[31] - Jw[29]) / (eps[31] - eps[29]) - slope_wrong) < 1e-3
print(f"Part 7: dJ/dε at 0: true solution {slope_true:+.1e} (flat), y = x^2 {slope_wrong:+.4f} (sloped)")
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(eps, Jt, label=r"$y = x$ (solution): flat at $\varepsilon=0$")
ax.plot(eps, Jw, label=r"$y = x^2$ (not a solution): sloped")
ax.axvline(0, color="gray", lw=.8); ax.set_xlabel(r"$\varepsilon$"); ax.set_ylabel(r"$J[y+\varepsilon\eta]$")
ax.legend(); ax.grid(alpha=.3); ax.set_title(r"Stationary means $dJ/d\varepsilon = 0$ for every $\eta$")
fig.tight_layout(); fig.savefig(FIG / "J_vs_eps.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 8: Gradient of J_h versus h × (Euler–Lagrange expression), eq. (3.9)
#
# Generic L = y sqrt(1+y'^2) + x y^2 on a smooth NON-solution y(x). The §3.5
# formula (3.8) is checked against finite differences of J_h, and then
# max_i |∂J_h/∂y_i − h·EL(x_i)| is shown to drop by 4 when h is halved.

# %%
X, Y, P = sp.symbols("x y p")
L_sym = Y * sp.sqrt(1 + P**2) + X * Y**2
L_f, Ly_f, Lp_f = (sp.lambdify((X, Y, P), e) for e in (L_sym, sp.diff(L_sym, Y), sp.diff(L_sym, P)))
y_sym = 1 + X**2 / 2 + sp.sin(2 * X) / 3
EL_sym = sp.diff(L_sym, Y).subs({Y: y_sym, P: sp.diff(y_sym, X)}) - sp.diff(sp.diff(L_sym, P).subs({Y: y_sym, P: sp.diff(y_sym, X)}), X)
y_f, EL_f = sp.lambdify(X, y_sym), sp.lambdify(X, EL_sym)


def J_h(y_int, x):
    """Left Riemann sum with forward-difference slopes, eq. (3.7)."""
    y = with_ends(y_int, y_f(x[0]), y_f(x[-1])); hh = x[1] - x[0]
    return hh * np.sum(L_f(x[:-1], y[:-1], np.diff(y) / hh))


def grad_J_h(y_int, x):
    """Eq. (3.8): h [ L_y(x_i, y_i, p_i) − (L_p(i) − L_p(i−1))/h ]."""
    y = with_ends(y_int, y_f(x[0]), y_f(x[-1])); hh = x[1] - x[0]; p = np.diff(y) / hh
    return hh * (Ly_f(x[1:-1], y[1:-1], p[1:]) - (Lp_f(x[1:-1], y[1:-1], p[1:]) - Lp_f(x[:-2], y[:-2], p[:-1])) / hh)


hs_list, errs = [], []
for Ng in (50, 100, 200, 400):
    xg = np.linspace(0, 1, Ng + 2); hg = xg[1] - xg[0]; yg = y_f(xg[1:-1])
    assert np.max(np.abs(approx_fprime(yg, J_h, 1e-7, xg) - grad_J_h(yg, xg))) < 1e-4   # (3.8) is the gradient
    errs.append(np.max(np.abs(grad_J_h(yg, xg) - hg * EL_f(xg[1:-1])))); hs_list.append(hg)
ratios = np.array(errs[:-1]) / np.array(errs[1:])
assert np.all(ratios > 3.5) and np.all(ratios < 4.5), ratios
print(f"Part 8: max|∂J_h/∂y_i − h·EL(x_i)| = {errs[0]:.2e}, {errs[1]:.2e}, {errs[2]:.2e}, {errs[3]:.2e}; "
      f"ratios on halving h: {np.round(ratios, 2)} -> O(h^2)")
fig, ax = plt.subplots(figsize=(5, 4))
ax.loglog(hs_list, errs, "o-", label=r"$\max_i|\partial J_h/\partial y_i - h\,\mathrm{EL}(x_i)|$")
ax.loglog(hs_list, errs[0] * (np.array(hs_list) / hs_list[0]) ** 2, "k--", label=r"$\propto h^2$")
ax.set_xlabel("h"); ax.legend(fontsize=8); ax.grid(alpha=.3, which="both"); ax.set_title("Discrete gradient vs Euler–Lagrange, eq. (3.9)")
fig.tight_layout(); fig.savefig(FIG / "gradient_vs_EL.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 9: SymPy checks of the identities, and the numbers in the exercises

# %%
xs_ = sp.symbols("x"); yf = sp.Function("y")(xs_); yp_ = sp.diff(yf, xs_)
EL_of = lambda L: sp.diff(L, yf) - sp.diff(sp.diff(L, yp_), xs_)
beltrami_defect = lambda L: sp.simplify(sp.diff(L - yp_ * sp.diff(L, yp_), xs_) - yp_ * EL_of(L))   # (3.12) with L_x = 0: dH/dx = y' · EL
for Lc in (sp.sqrt(1 + yp_**2), sp.sqrt((1 + yp_**2) / yf), yf * sp.sqrt(1 + yp_**2), yp_**2 / 2 - yf):
    assert beltrami_defect(Lc) == 0, Lc
# Example 4.5: falling particle. H = −(y'^2/2 + y); dH/dx = −y'(y'' + 1); constant on the solution.
L45 = yp_**2 / 2 - yf
H45 = L45 - yp_ * sp.diff(L45, yp_)
assert sp.simplify(H45 + (yp_**2 / 2 + yf)) == 0 and sp.simplify(EL_of(L45) + sp.diff(yf, xs_, 2) + 1) == 0
v0, yy0 = sp.symbols("v0 y0")
sol45 = -xs_**2 / 2 + v0 * xs_ + yy0
assert sp.simplify((yp_**2 / 2 + yf).subs(yf, sol45).doit() - (v0**2 / 2 + yy0)) == 0
assert sp.diff(H45.subs(yf, xs_**3).doit(), xs_) != 0   # not constant on a non-solution
# Example 4.1 / 4.2 / 4.4: the named curves satisfy their E–L equations.
assert sp.simplify(EL_of(sp.sqrt(1 + yp_**2)) * (1 + yp_**2) ** sp.Rational(3, 2) + sp.diff(yf, xs_, 2)) == 0
tht, aa, cc, x0s = sp.symbols("theta a c x0", positive=True)
Xc, Yc = aa * (tht - sp.sin(tht)), aa * (1 - sp.cos(tht))
dydx = sp.diff(Yc, tht) / sp.diff(Xc, tht)
assert sp.simplify(Yc * (1 + dydx**2) - 2 * aa) == 0                                      # Beltrami first integral
assert sp.simplify(2 * Yc * sp.diff(dydx, tht) / sp.diff(Xc, tht) + 1 + dydx**2) == 0      # second-order E–L: 2yy'' + 1 + y'^2 = 0
assert sp.simplify(EL_of(yf * sp.sqrt(1 + yp_**2)).subs(yf, cc * sp.cosh((xs_ - x0s) / cc)).doit()) == 0   # catenary / catenoid
print("Part 9: SymPy: Beltrami identity on four Lagrangians, falling-particle energy, line, cycloid and cosh all check.")

# Example 4.6: catenoid exists iff d/R <= 1 / min_u (cosh u / u); the minimizer solves u tanh u = 1.
u_star = brentq(lambda u: u * np.tanh(u) - 1, 0.5, 2)
d_crit = u_star / np.cosh(u_star)
assert np.isclose(d_crit, 0.6627, atol=5e-5)
# Exercise B6: y = 6x(1−x), λ = −24, J = 12.   B7: free end on x = 1 -> θ1 = π, T = sqrt(π/g).   C2: same height, T = sqrt(2π/g).
assert np.isclose(quad(lambda t: 6 * t * (1 - t), 0, 1)[0], 1) and np.isclose(quad(lambda t: (6 - 12 * t) ** 2, 0, 1)[0], 12)
th_B7, a_B7 = cycloid_params(1.0, 2 / np.pi)
assert np.isclose(th_B7, np.pi) and np.isclose(th_B7 * np.sqrt(a_B7 / g), np.sqrt(np.pi / g))
assert np.isclose(2 * np.pi * np.sqrt((1 / (2 * np.pi)) / g), np.sqrt(2 * np.pi / g))
print(f"        catenoid critical half-separation d/R = {d_crit:.4f}; B7: T = sqrt(π/g) = {np.sqrt(np.pi / g):.4f} s; "
      f"C2: T = sqrt(2π/g) = {np.sqrt(2 * np.pi / g):.4f} s")

# %% [markdown]
# ## Your turn
#
# 1. **Catenoid.** Minimize the surface area 2π Σ Δs_i (y_i + y_{i+1})/2 of the
#    polygon between rings y(±d) = 1 (no length constraint: Example 4.6 with
#    λ = 0). For d = 0.5 compare with c cosh(x/c), where c solves c cosh(d/c) = 1
#    (two roots; the smaller area belongs to the larger c). Then try d = 0.8 > 0.6627
#    and watch the minimizer pinch toward the axis: no catenoid exists, and the
#    polygon heads for the two-disk (Goldschmidt) configuration. Exercise B8.
# 2. **Brachistochrone with friction.** With Coulomb friction μ on a wire whose
#    slope stays below 45°, energy gives v^2 = 2g(y − μx) (y down). Replace the
#    speed law in `polygon_descent_time`, minimize for μ = 0.1 from (0,0) to (1,1),
#    and compare the shape and time with the frictionless cycloid. Does the start
#    still look vertical?
# 3. **Newton's law from the action.** Store q(t) at N interior times on [0, T]
#    with q(0) = 0, q(T) = 0, T = 1, and minimize the discretized action
#    Σ h [ (Δq_i/h)^2 / 2 − g (q_i + q_{i+1})/2 ]. Compare the minimizer with the
#    parabola q = (g/2) t (T − t) (that is q̈ = −g with these endpoints), and check
#    that the second difference of the numerical q equals −g at every node. This
#    is §3.13 done by the optimizer; CM-2 does it by hand.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
