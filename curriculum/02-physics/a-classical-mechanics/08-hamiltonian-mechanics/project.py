# %% [markdown]
# # CM-8 Project: Phase Space, Liouville, Poisson Brackets, Canonical Maps, Hamilton–Jacobi
#
# This script builds the geometry of Hamiltonian mechanics and checks every
# numerical claim in the notes. Parts: (1) the Legendre transform as a
# maximization and its tangent-line picture, (2) phase portraits for the
# oscillator, the pendulum (with separatrix) and a double well, (3) the
# spherical pendulum reduced by p_phi (Example 4.2), (4) Liouville's theorem:
# a blob of 5000 points transported by three integrators, with its area
# tracked, (5) the module's signature figure: energy drift over 10^4 periods,
# RK4 vs Störmer–Verlet, (6) a SymPy Poisson bracket that verifies the
# angular-momentum algebra, Jacobi's identity and conservation laws,
# (7) a canonical-transformation checker (M^T J M = J), (8) action-angle
# variables for the pendulum and nu = dE/dJ, (9) adiabatic invariance,
# (10) the Hamilton–Jacobi function of the oscillator checked against the
# PDE and against p along a trajectory, (11) WKB: |psi|^2 ~ 1/p vs the exact
# n = 20 oscillator state, (12) symplectic Euler as two shears and its exactly
# conserved modified Hamiltonian, (13) checks for every number quoted in the
# exercise solutions.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import math
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from pathlib import Path
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize_scalar
from scipy.special import ellipk, ellipe, eval_hermite, gammaln

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T_START = time.time()

# %% [markdown]
# ## Part 1: The Legendre transform as a maximization (§3.1, Example 4.1)
#
# For convex f, g(p) = max_v [p v - f(v)]  (3.1). The maximizer satisfies
# p = f'(v) (3.2) and g'(p) = v (3.3); applying the transform twice returns f.
# Geometrically, -g(p) is the y-intercept of the tangent line to f of slope p.
#
# **Build it yourself first:** write `legendre(f, p)` with `minimize_scalar`
# before reading this one.


# %%
def legendre(f, p, bracket=(-20, 20)):
    """g(p) = max_v [p v - f(v)], eq. (3.1), by numerical maximization over v."""
    res = minimize_scalar(lambda v: -(p * v - f(v)), bounds=bracket, method="bounded",
                          options={"xatol": 1e-10})
    return -res.fun, res.x


m_ex = 2.0
f_kin = lambda v: 0.5 * m_ex * v**2
f_exp = lambda v: np.exp(v)
for p in (0.5, 1.0, 2.0, 3.0):
    g, v_star = legendre(f_kin, p)
    assert np.isclose(g, p**2 / (2 * m_ex), atol=1e-7), (p, g)      # g = p^2/2m
    assert np.isclose(v_star, p / m_ex, atol=1e-5)                   # v = p/m
    g2, v2 = legendre(f_exp, p)
    assert np.isclose(g2, p * np.log(p) - p, atol=1e-7), (p, g2)     # g = p ln p - p
    assert np.isclose(v2, np.log(p), atol=1e-5)                      # v = ln p
g_at_2 = 2 * np.log(2) - 2
print(f"Example 4.1: Legendre of e^v at p = 2 is 2 ln 2 - 2 = {g_at_2:.5f}; of m v^2/2 (m=2) at p=3 is {9/4:.4f}")

# Involution: the transform of g is f again.
g_fun = lambda p: p * np.log(p) - p if p > 0 else np.inf
for v in (-1.0, 0.0, 0.7, 1.5):
    f_back, p_star = legendre(g_fun, v, bracket=(1e-9, 50))
    assert np.isclose(f_back, np.exp(v), atol=1e-6), (v, f_back)
print("Legendre transform is an involution on e^v: checked at four points.")

# Tangent-line picture: f = v^2/2 (m = 1); the tangent of slope p has intercept -g(p) = -p^2/2.
v = np.linspace(-2.5, 2.5, 300)
fig, ax = plt.subplots(figsize=(5.5, 4.5))
ax.plot(v, 0.5 * v**2, "k", lw=2, label=r"$f(v)=v^2/2$")
for p in (-2, -1, 0.5, 1.5, 2.2):
    ax.plot(v, p * v - 0.5 * p**2, lw=1, alpha=.7)
    ax.plot(0, -0.5 * p**2, "o", ms=4, color="C3")
ax.set_ylim(-3, 3.5); ax.set_xlabel("v"); ax.set_ylabel("y")
ax.set_title("Tangent lines of slope p hit the y-axis at -g(p) = -p²/2")
ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "legendre_tangent_lines.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 2: Phase portraits (§3.4, Example 4.3)
#
# For one degree of freedom, trajectories are level sets of H. We draw the
# vector field (dq/dt, dp/dt) = (dH/dp, -dH/dq) and contours of H for the
# oscillator, the pendulum (separatrix E = m g l, drawn thick), and the double
# well V = -x^2/2 + x^4/4 (separatrix E = 0).


# %%
def H_osc(q, p, m=1.0, w=1.0):
    return p**2 / (2 * m) + 0.5 * m * w**2 * q**2


def H_pend(th, p, m=1.0, l=1.0, g=1.0):
    """Pendulum Hamiltonian H = p^2 / (2 m l^2) - m g l cos(theta), theta from the downward vertical."""
    return p**2 / (2 * m * l**2) - m * g * l * np.cos(th)


def H_dwell(x, p):
    return p**2 / 2 - x**2 / 2 + x**4 / 4


def portrait(ax, H, qlim, plim, levels, title, sep=None):
    q = np.linspace(*qlim, 300); p = np.linspace(*plim, 300)
    Q, P = np.meshgrid(q, p)
    ax.contour(Q, P, H(Q, P), levels=levels, colors="C0", linewidths=.8)
    if sep is not None:
        ax.contour(Q, P, H(Q, P), levels=[sep], colors="C3", linewidths=2.2)
    qq = np.linspace(*qlim, 18); pp = np.linspace(*plim, 14)
    QQ, PP = np.meshgrid(qq, pp)
    eps = 1e-6
    dq = (H(QQ, PP + eps) - H(QQ, PP - eps)) / (2 * eps)       # dH/dp
    dp = -(H(QQ + eps, PP) - H(QQ - eps, PP)) / (2 * eps)      # -dH/dq
    ax.quiver(QQ, PP, dq, dp, color="gray", alpha=.6, width=.003)
    ax.set_title(title); ax.set_xlabel("q"); ax.set_ylabel("p")


fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
portrait(axs[0], H_osc, (-2.5, 2.5), (-2.5, 2.5), np.linspace(.1, 3, 10), "Harmonic oscillator: ellipses")
portrait(axs[1], H_pend, (-2.2 * np.pi, 2.2 * np.pi), (-3.2, 3.2), np.linspace(-.9, 4, 14),
         "Pendulum: libration inside, rotation outside", sep=1.0)
portrait(axs[2], H_dwell, (-1.8, 1.8), (-1.3, 1.3), np.linspace(-.24, .8, 12),
         "Double well: figure-eight separatrix", sep=0.0)
fig.tight_layout(); fig.savefig(FIG / "phase_portraits.png", dpi=120); plt.close(fig)

# Example 4.3: separatrix energy and the area it encloses, 16 m l sqrt(g l).
m, l, g = 1.0, 1.0, 9.81
E_sep = m * g * l                                       # H at the hyperbolic point (pi, 0)
assert np.isclose(H_pend(np.pi, 0, m, l, g), E_sep)
p_sep = lambda th: 2 * m * l * np.sqrt(g * l) * np.cos(th / 2)     # separatrix branch p >= 0
# check it is a level set of H at E_sep
th_test = np.linspace(-np.pi, np.pi, 50)
assert np.allclose(H_pend(th_test, p_sep(th_test), m, l, g), E_sep)
area_sep, _ = quad(lambda th: 2 * p_sep(th), -np.pi, np.pi)           # upper + lower branch
assert np.isclose(area_sep, 16 * m * l * np.sqrt(g * l), rtol=1e-10)
print(f"Example 4.3: separatrix energy E_s = m g l = {E_sep:.2f} J; enclosed area = 16 m l sqrt(g l) = {area_sep:.3f}")
# Motion on the separatrix: sin(theta/2) = tanh(w0 t) solves theta_dot = 2 w0 cos(theta/2), reaching pi only as t -> inf.
w0 = np.sqrt(g / l); ts_ = np.linspace(-2, 2, 4001)
th_sx = 2 * np.arcsin(np.tanh(w0 * ts_))
dth_num = np.gradient(th_sx, ts_, edge_order=2)
assert np.allclose(dth_num, 2 * w0 * np.cos(th_sx / 2), rtol=2e-3, atol=1e-6)
t_sym = sp.symbols("t", real=True)
assert sp.simplify(sp.diff(2 * sp.asin(sp.tanh(t_sym)), t_sym) - 2 * sp.cos(sp.asin(sp.tanh(t_sym)))) == 0

# %% [markdown]
# ## Part 3: Spherical pendulum reduced by p_phi (Example 4.2)
#
# H = p_th^2/(2 m l^2) + p_phi^2/(2 m l^2 sin^2 th) - m g l cos th. phi is cyclic,
# so p_phi is conserved and theta moves in V_eff. The conical orbit sits at
# V_eff'(theta_0) = 0, which gives p_phi^2 = m^2 g l^3 sin^4(th0)/cos(th0) and
# phidot^2 = g/(l cos th0). Small oscillations about it have
# Omega^2 = (g/l)(1 + 3 cos^2 th0)/cos th0.

# %%
m, l, g = 1.0, 1.0, 9.81
th0 = np.pi / 6
p_phi = np.sqrt(m**2 * g * l**3 * np.sin(th0)**4 / np.cos(th0))
phidot_pred = np.sqrt(g / (l * np.cos(th0)))
Omega_pred = np.sqrt((g / l) * (1 + 3 * np.cos(th0)**2) / np.cos(th0))
print(f"Example 4.2: p_phi = {p_phi:.6f}, phidot = {phidot_pred:.5f} rad/s, "
      f"T_orbit = {2*np.pi/phidot_pred:.5f} s, Omega_theta = {Omega_pred:.4f} rad/s, ratio = {Omega_pred/phidot_pred:.4f}")
assert np.isclose(p_phi, 0.841413, atol=1e-6) and np.isclose(phidot_pred, 3.36565, atol=1e-5)
assert np.isclose(Omega_pred / phidot_pred, np.sqrt(1 + 3 * np.cos(th0)**2))

# Symbolic check that V_eff'(th0) = 0 with that p_phi.
th_s, pphi_s, m_s, l_s, g_s = sp.symbols("theta p_phi m l g", positive=True)
V_eff = pphi_s**2 / (2 * m_s * l_s**2 * sp.sin(th_s)**2) - m_s * g_s * l_s * sp.cos(th_s)
dV = sp.diff(V_eff, th_s).subs(pphi_s**2, m_s**2 * g_s * l_s**3 * sp.sin(th_s)**4 / sp.cos(th_s))
assert sp.simplify(dV) == 0
d2V = sp.diff(V_eff, th_s, 2).subs(pphi_s**2, m_s**2 * g_s * l_s**3 * sp.sin(th_s)**4 / sp.cos(th_s))
assert sp.simplify(d2V - m_s * g_s * l_s * (1 + 3 * sp.cos(th_s)**2) / sp.cos(th_s)) == 0


def spherical_rhs(t, y):
    """Hamilton's equations for (theta, p_theta, phi) at fixed p_phi (Example 4.2)."""
    th, pth, phi = y
    dth = pth / (m * l**2)
    dpth = p_phi**2 * np.cos(th) / (m * l**2 * np.sin(th)**3) - m * g * l * np.sin(th)
    dphi = p_phi / (m * l**2 * np.sin(th)**2)
    return [dth, dpth, dphi]


sol = solve_ivp(spherical_rhs, (0, 5), [th0, 0.0, 0.0], rtol=1e-10, atol=1e-12, dense_output=True)
assert np.max(np.abs(sol.y[0] - th0)) < 1e-10, "conical orbit should keep theta fixed (notes: to about 1e-12)"
assert np.isclose(sol.y[2, -1] / sol.t[-1], phidot_pred, rtol=1e-8)
# Perturb theta slightly and measure the oscillation frequency about the conical orbit.
sol2 = solve_ivp(spherical_rhs, (0, 20), [th0 + 0.01, 0.0, 0.0], rtol=1e-10, atol=1e-12,
                 t_eval=np.linspace(0, 20, 20001))
th_dev = sol2.y[0] - th0
crossings = np.where(np.diff(np.sign(th_dev)) != 0)[0]
t_cross = sol2.t[crossings]
Omega_meas = np.pi / np.mean(np.diff(t_cross))
assert np.isclose(Omega_meas, Omega_pred, rtol=2e-3), (Omega_meas, Omega_pred)
assert np.isclose(Omega_pred, 6.0675, atol=5e-5) and np.isclose(Omega_meas, 6.0673, atol=1e-4)   # Example 4.2 numbers
print(f"  measured Omega_theta = {Omega_meas:.4f} (predicted {Omega_pred:.4f}); theta stays at th0 to {np.max(np.abs(sol.y[0]-th0)):.1e}")

# %% [markdown]
# ## Part 4: Liouville's theorem with a blob of points (§3.5)
#
# A disk of 5000 initial conditions in the pendulum's phase space is carried
# by the flow with three integrators: scipy's RK45 (adaptive, not symplectic),
# symplectic Euler (an exact area-preserving map), and explicit Euler (area
# grows by the factor 1 + h^2 cos(theta) per step). The area is measured by the
# shoelace formula on the transported *boundary* polygon (4000 points), which is
# the exact area of the image region up to polygon resolution.

# %%
m, l, g = 1.0, 1.0, 1.0           # H = p^2/2 - cos(theta)
N_IN, N_BD = 5000, 4000
r_disk, c_disk = 0.5, np.array([1.5, 0.0])
rr = r_disk * np.sqrt(rng.uniform(size=N_IN)); aa = rng.uniform(0, 2 * np.pi, N_IN)
interior0 = np.c_[c_disk[0] + rr * np.cos(aa), c_disk[1] + rr * np.sin(aa)]
ang = np.linspace(0, 2 * np.pi, N_BD, endpoint=False)
boundary0 = np.c_[c_disk[0] + r_disk * np.cos(ang), c_disk[1] + r_disk * np.sin(ang)]
pts0 = np.vstack([interior0, boundary0])        # shape (9000, 2)
assert np.all(H_pend(pts0[:, 0], pts0[:, 1]) < 1.0), "all points librate"


def shoelace(P):
    """Area enclosed by the closed polygon with vertices P (n x 2)."""
    x, y = P[:, 0], P[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def pend_rhs_batched(t, y):
    th, p = y[:len(y) // 2], y[len(y) // 2:]
    return np.concatenate([p, -np.sin(th)])


def run_rk45(pts, t_eval):
    y0 = np.concatenate([pts[:, 0], pts[:, 1]])
    s = solve_ivp(pend_rhs_batched, (0, t_eval[-1]), y0, t_eval=t_eval, rtol=1e-3, atol=1e-6)
    return [np.c_[s.y[:len(pts), k], s.y[len(pts):, k]] for k in range(len(t_eval))]


def run_map(pts, t_eval, h, step):
    out, cur, t = [], pts.copy(), 0.0
    for tk in t_eval:
        n_steps = int(round((tk - t) / h))
        for _ in range(n_steps):
            cur = step(cur, h)
        t = tk; out.append(cur.copy())
    return out


def symplectic_euler_step(P, h):
    """Kick then drift: p' = p - h sin(theta); theta' = theta + h p'. Each is a shear (det 1)."""
    th, p = P[:, 0], P[:, 1]
    p1 = p - h * np.sin(th)
    return np.c_[th + h * p1, p1]


def explicit_euler_step(P, h):
    th, p = P[:, 0], P[:, 1]
    return np.c_[th + h * p, p - h * np.sin(th)]


t_eval = np.array([0.0, 8.0, 16.0, 24.0])
h_map = 0.02
results = {"RK45 (rtol 1e-3)": run_rk45(pts0, t_eval),
           "symplectic Euler (h=0.02)": run_map(pts0, t_eval, h_map, symplectic_euler_step),
           "explicit Euler (h=0.02)": run_map(pts0, t_eval, h_map, explicit_euler_step)}
A0 = shoelace(boundary0)
assert np.isclose(A0, np.pi * r_disk**2, rtol=1e-5)
fig, axs = plt.subplots(3, 4, figsize=(16, 10.5))
print("Liouville blob areas (relative to pi r^2):")
areas = {}
for i, (name, frames) in enumerate(results.items()):
    areas[name] = []
    for k, (tk, P) in enumerate(zip(t_eval, frames)):
        A = shoelace(P[N_IN:])
        areas[name].append(A / A0)
        ax = axs[i, k]
        ax.scatter(P[:N_IN, 0], P[:N_IN, 1], s=1, alpha=.4, color="C0")
        ax.plot(np.r_[P[N_IN:, 0], P[N_IN, 0]], np.r_[P[N_IN:, 1], P[N_IN, 1]], "C3", lw=.8)
        ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.7, 1.7)
        ax.set_title(f"{name}\nt = {tk:g}, area/area0 = {A/A0:.5f}", fontsize=9)
        ax.set_xlabel(r"$\theta$"); ax.set_ylabel("p")
    print(f"  {name:28s}: " + "  ".join(f"{a:.6f}" for a in areas[name]))
fig.suptitle("Liouville: the blob deforms but its area is invariant under a canonical map", y=1.0)
fig.tight_layout(); fig.savefig(FIG / "liouville_blob.png", dpi=110); plt.close(fig)
assert max(abs(a - 1) for a in areas["symplectic Euler (h=0.02)"]) < 1e-5, "symplectic map must preserve area"
assert np.isclose(areas["symplectic Euler (h=0.02)"][-1], 0.999998, atol=1e-6)        # §7: 0.999998 at t = 24
assert np.isclose(areas["RK45 (rtol 1e-3)"][-1], 0.960, atol=5e-4)                   # §7: 0.960
assert np.isclose(areas["explicit Euler (h=0.02)"][-1], 1.202, atol=5e-4), "explicit Euler should inflate the blob"   # §7: 1.202
# Predicted inflation for explicit Euler: product of (1 + h^2 cos theta) along the path ~ exp(h^2 * sum cos theta).
# Check this on the centroid trajectory as a sanity estimate.
P = pts0[:N_IN].mean(axis=0)[None, :]; log_growth = 0.0
for _ in range(int(round(24.0 / h_map))):
    log_growth += np.log(1 + h_map**2 * np.cos(P[0, 0])); P = explicit_euler_step(P, h_map)
print(f"  explicit Euler predicted area factor (centroid estimate) exp(sum ln(1+h^2 cos th)) = {np.exp(log_growth):.4f}")

# %% [markdown]
# ## Part 5: Signature figure — energy drift over 10^4 periods, RK4 vs Verlet (§3.11)
#
# Pendulum with theta0 = 2 rad. Both methods use 50 steps per period. RK4 is
# 4th order but not symplectic: its energy error grows secularly. Verlet is
# only 2nd order but symplectic: its energy error oscillates but does not grow
# (for exponentially long times, §3.11).

# %%
th_amp = 2.0
k_mod = np.sin(th_amp / 2)**2
T_pend = 4 * ellipk(k_mod)                   # period, l = g = 1 (Exercise C1)
N_PER, STEPS = 10_000, 50
h = T_pend / STEPS
n_total = N_PER * STEPS
E0 = H_pend(th_amp, 0.0)
sin, cos = math.sin, math.cos


def rk4_pendulum(th, p, h, n, every):
    out = []
    for i in range(n):
        k1t, k1p = p, -sin(th)
        k2t, k2p = p + 0.5 * h * k1p, -sin(th + 0.5 * h * k1t)
        k3t, k3p = p + 0.5 * h * k2p, -sin(th + 0.5 * h * k2t)
        k4t, k4p = p + h * k3p, -sin(th + h * k3t)
        th += h * (k1t + 2 * k2t + 2 * k3t + k4t) / 6
        p += h * (k1p + 2 * k2p + 2 * k3p + k4p) / 6
        if i % every == every - 1:
            out.append(0.5 * p * p - cos(th))
    return np.array(out)


def verlet_pendulum(th, p, h, n, every):
    """Störmer–Verlet (3.32): half kick, drift, half kick."""
    out = []
    for i in range(n):
        p -= 0.5 * h * sin(th)
        th += h * p
        p -= 0.5 * h * sin(th)
        if i % every == every - 1:
            out.append(0.5 * p * p - cos(th))
    return np.array(out)


every = STEPS                        # record once per period
t0 = time.time()
E_rk4 = rk4_pendulum(th_amp, 0.0, h, n_total, every)
E_ver = verlet_pendulum(th_amp, 0.0, h, n_total, every)
t_periods = np.arange(1, N_PER + 1)
err_rk4 = np.abs(E_rk4 - E0) / abs(E0); err_ver = np.abs(E_ver - E0) / abs(E0)
print(f"Energy drift after {N_PER} periods (h = T/{STEPS}): RK4 {err_rk4[-1]:.3e}, Verlet {err_ver[-1]:.3e}  "
      f"(Verlet max {err_ver.max():.3e}); loop time {time.time()-t0:.1f}s")
assert err_ver.max() < 2e-2, "Verlet energy error must stay bounded"
# §3.11 numbers: RK4 0.257 at 10^4 periods; Verlet max 0.0181, 0.0019 at the end
assert np.isclose(err_rk4[-1], 0.257, atol=1e-3) and np.isclose(err_ver.max(), 0.0181, atol=1e-4) and np.isclose(err_ver[-1], 0.0019, atol=1e-4)
assert err_rk4[-1] > 10 * err_ver.max(), "RK4 should have drifted far past Verlet's bounded error"
# RK4 drift should be roughly linear in time: compare error at 10^3 and 10^4 periods.
ratio = err_rk4[-1] / err_rk4[N_PER // 10 - 1]
assert 5 < ratio < 20, ratio
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.loglog(t_periods, err_rk4 + 1e-17, label="RK4 (4th order, not symplectic)")
ax.loglog(t_periods, err_ver + 1e-17, label="Störmer–Verlet (2nd order, symplectic)", alpha=.8)
ax.set_xlabel("time (periods)"); ax.set_ylabel(r"$|E - E_0| / |E_0|$")
ax.set_title(f"Pendulum, $\\theta_0 = 2$ rad, {STEPS} steps per period, $10^4$ periods")
ax.grid(alpha=.3, which="both"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "energy_drift.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 6: Poisson brackets in SymPy (§3.6, Example 4.4)
#
# {F, G} = sum_i (dF/dq_i dG/dp_i - dF/dp_i dG/dq_i), eq. (3.12). We verify the
# angular-momentum algebra {L_i, L_j} = eps_ijk L_k (3.17), {L^2, L_i} = 0,
# Jacobi's identity (3.13) on random polynomials, {F, p} = dF/dq (3.18),
# conservation {L_z, H} = {L^2, H} = 0 in a central potential, and the
# so(3)-like algebra of the 2D isotropic oscillator (Exercise B3).

# %%
def poisson_bracket(F, G, qs, ps):
    """Poisson bracket {F, G} with respect to canonical pairs (qs[i], ps[i]), eq. (3.12)."""
    return sp.expand(sum(sp.diff(F, q) * sp.diff(G, p) - sp.diff(F, p) * sp.diff(G, q)
                         for q, p in zip(qs, ps)))


x, y, z, px, py, pz = sp.symbols("x y z p_x p_y p_z", real=True)
qs, ps = [x, y, z], [px, py, pz]
Lx, Ly, Lz = y * pz - z * py, z * px - x * pz, x * py - y * px
L = [Lx, Ly, Lz]
# fundamental brackets
for i in range(3):
    for j in range(3):
        assert poisson_bracket(qs[i], ps[j], qs, ps) == (1 if i == j else 0)
        assert poisson_bracket(qs[i], qs[j], qs, ps) == 0 and poisson_bracket(ps[i], ps[j], qs, ps) == 0
# angular momentum algebra (Example 4.4 is {Lx, Ly} = Lz)
assert poisson_bracket(Lx, Ly, qs, ps) == sp.expand(Lz)
eps = sp.LeviCivita
for i in range(3):
    for j in range(3):
        rhs = sp.expand(sum(eps(i, j, k) * L[k] for k in range(3)))
        assert poisson_bracket(L[i], L[j], qs, ps) == rhs, (i, j)
L2 = sp.expand(Lx**2 + Ly**2 + Lz**2)
assert all(poisson_bracket(L2, Li, qs, ps) == 0 for Li in L)
print("Poisson brackets: {L_i, L_j} = eps_ijk L_k and {L^2, L_i} = 0 verified symbolically.")

# Jacobi identity on random polynomials of degree <= 2 in (x, y, p_x, p_y)
q2, p2 = [x, y], [px, py]
monos = [1, x, y, px, py, x*y, x*px, x*py, y*px, y*py, x**2, y**2, px**2, py**2, px*py]
def rand_poly():
    return sum(int(c) * mnm for c, mnm in zip(rng.integers(-3, 4, len(monos)), monos))
for _ in range(5):
    F, G, K = rand_poly(), rand_poly(), rand_poly()
    jac = (poisson_bracket(F, poisson_bracket(G, K, q2, p2), q2, p2)
           + poisson_bracket(G, poisson_bracket(K, F, q2, p2), q2, p2)
           + poisson_bracket(K, poisson_bracket(F, G, q2, p2), q2, p2))
    assert sp.expand(jac) == 0
    # Leibniz rule and antisymmetry
    assert sp.expand(poisson_bracket(F * G, K, q2, p2) - F * poisson_bracket(G, K, q2, p2) - G * poisson_bracket(F, K, q2, p2)) == 0
    assert sp.expand(poisson_bracket(F, G, q2, p2) + poisson_bracket(G, F, q2, p2)) == 0
    # {F, p_x} = dF/dx (3.18): p generates translations
    assert sp.expand(poisson_bracket(F, px, q2, p2) - sp.diff(F, x)) == 0
print("Jacobi identity, Leibniz rule, antisymmetry, and {F, p} = dF/dq verified on 5 random polynomial triples.")

# Conservation in a central potential: {L_z, H} = {L^2, H} = 0 for H = |p|^2/2m + V(r).
r = sp.sqrt(x**2 + y**2 + z**2); Vf = sp.Function("V")
m_s = sp.symbols("m", positive=True)
H_cent = (px**2 + py**2 + pz**2) / (2 * m_s) + Vf(r)
assert sp.simplify(poisson_bracket(Lz, H_cent, qs, ps)) == 0
assert sp.simplify(poisson_bracket(L2, H_cent, qs, ps)) == 0
# Infinitesimal rotation generated by L_z: {x, L_z} = -y, {y, L_z} = x (§3.6)
assert poisson_bracket(x, Lz, qs, ps) == -y and poisson_bracket(y, Lz, qs, ps) == x
# Exercise B3: 2D isotropic oscillator (m = omega = 1): A, B, L_z all conserved, and {A,B} = 2L, {B,L} = 2A, {L,A} = 2B
H2 = (px**2 + py**2) / 2 + (x**2 + y**2) / 2
A2 = (px**2 - py**2) / 2 + (x**2 - y**2) / 2; B2 = px * py + x * y; L2d = x * py - y * px
for Fc in (A2, B2, L2d):
    assert poisson_bracket(Fc, H2, q2, p2) == 0
assert poisson_bracket(A2, B2, q2, p2) == sp.expand(2 * L2d)
assert poisson_bracket(B2, L2d, q2, p2) == sp.expand(2 * A2)
assert poisson_bracket(L2d, A2, q2, p2) == sp.expand(2 * B2)
print("Central potential: {L_z, H} = {L^2, H} = 0.  2D oscillator: {A,B} = 2L_z, {B,L_z} = 2A, {L_z,A} = 2B (B3).")

# %% [markdown]
# ## Part 7: Canonical-transformation checker (§3.7, Example 4.5)
#
# A map z -> Z is canonical iff its Jacobian M = dZ/dz satisfies M^T J M = J
# (3.20), with J = [[0, I], [-I, 0]]. For one degree of freedom this is
# det M = 1. We test the oscillator's action-angle map, the exchange map, a
# point transformation, and a non-canonical map.

# %%
def symplectic_matrix(n):
    Jm = sp.zeros(2 * n, 2 * n)
    Jm[:n, n:] = sp.eye(n); Jm[n:, :n] = -sp.eye(n)
    return Jm


def canonical_check(Z, z):
    """Return (M, M^T J M - J) for the map z -> Z given as sympy expressions; zero matrix means canonical."""
    n = len(z) // 2
    M = sp.Matrix([[sp.diff(Zi, zj) for zj in z] for Zi in Z])
    Jm = symplectic_matrix(n)
    return M, sp.simplify(M.T * Jm * M - Jm)


q, p, Q, P = sp.symbols("q p Q P", real=True)
w_s = sp.symbols("omega", positive=True)
# Example 4.5: (Q, P) -> (q, p); the map is canonical iff its inverse is, so checking this direction suffices.
q_of = sp.sqrt(2 * P / (m_s * w_s)) * sp.sin(Q); p_of = sp.sqrt(2 * P * m_s * w_s) * sp.cos(Q)
M_aa, resid = canonical_check([q_of, p_of], [Q, P])
assert resid == sp.zeros(2, 2), resid
assert sp.simplify(M_aa.det() - 1) == 0
H_new = sp.simplify(p_of**2 / (2 * m_s) + m_s * w_s**2 * q_of**2 / 2)
assert sp.simplify(H_new - w_s * P) == 0
print("Example 4.5: action-angle map has M^T J M = J (det M = 1) and H(q(Q,P), p(Q,P)) =", H_new)
# ... and the inverse map (q, p) -> (Q, P), checked directly: {Q, P}_{q,p} = 1
Q_of = sp.atan(m_s * w_s * q / p); P_of = p**2 / (2 * m_s * w_s) + m_s * w_s * q**2 / 2
assert sp.simplify(poisson_bracket(Q_of, P_of, [q], [p])) == 1
# Exchange, point transformation, and a non-canonical map
_, r1 = canonical_check([p, -q], [q, p]); assert r1 == sp.zeros(2, 2)
f = sp.Function("f")
_, r2 = canonical_check([f(q), p / sp.diff(f(q), q)], [q, p]); assert r2 == sp.zeros(2, 2)
_, r3 = canonical_check([q, p**2], [q, p]); assert r3 != sp.zeros(2, 2)
print("  exchange (Q,P) = (p,-q): canonical; point map (f(q), p/f'(q)): canonical; (q, p^2): NOT canonical, residual =", r3.tolist())
# Numerical Jacobian check at random points for the arctan map (float), and for a 2-DOF map.
Jn = np.array(symplectic_matrix(1), dtype=float)
def num_jac(fun, z, eps=1e-6):
    z = np.asarray(z, float); n = len(z); M = np.zeros((n, n))
    for j in range(n):
        dz = np.zeros(n); dz[j] = eps
        M[:, j] = (np.asarray(fun(z + dz)) - np.asarray(fun(z - dz))) / (2 * eps)
    return M
aa_map = lambda z: (np.arctan(z[0] / z[1]), 0.5 * (z[0]**2 + z[1]**2))        # m = omega = 1
for _ in range(20):
    z0 = rng.normal(size=2)
    Mn = num_jac(aa_map, z0)
    assert np.allclose(Mn.T @ Jn @ Mn, Jn, atol=1e-6)
# 2 DOF: polar coordinates as a point transformation (r, phi, p_r, p_phi) from (x, y, p_x, p_y)
def polar_map(zv):
    xx, yy, pxx, pyy = zv
    rr_ = np.hypot(xx, yy); ph = np.arctan2(yy, xx)
    return (rr_, ph, (xx * pxx + yy * pyy) / rr_, xx * pyy - yy * pxx)
J2 = np.array(symplectic_matrix(2), dtype=float)
for _ in range(20):
    z0 = rng.normal(size=4); Mn = num_jac(polar_map, z0)
    assert np.allclose(Mn.T @ J2 @ Mn, J2, atol=1e-5)
print("  numerical Jacobians: arctan map (1 DOF) and polar point transformation (2 DOF) satisfy M^T J M = J.")

# %% [markdown]
# ## Part 8: Action-angle variables for the pendulum (§3.9, Exercise C1)
#
# J(E) = closed integral of p dtheta = area enclosed by the orbit of energy E.
# Closed form (C1): J = 16 m l sqrt(g l) [E(k) - (1 - k^2) K(k)], k = sin(theta0/2),
# and the period is T = 4 sqrt(l/g) K(k). The frequency is nu = dE/dJ = 1/T.

# %%
m, l, g = 1.0, 1.0, 1.0


def action_numeric(E):
    """J(E) = 2 * integral_{-th0}^{th0} sqrt(2 m l^2 (E + m g l cos th)) dth, for libration (-mgl < E < mgl)."""
    th0 = np.arccos(-E / (m * g * l))
    pfun = lambda th: np.sqrt(np.maximum(2 * m * l**2 * (E + m * g * l * np.cos(th)), 0.0))
    val, _ = quad(pfun, -th0, th0, limit=200)
    return 2 * val


def action_closed(E):
    k2 = (1 + E / (m * g * l)) / 2                                  # k^2 = sin^2(th0/2) = (1 + E/mgl)/2
    return 16 * m * l * np.sqrt(g * l) * (ellipe(k2) - (1 - k2) * ellipk(k2))


def period_closed(E):
    k2 = (1 + E / (m * g * l)) / 2
    return 4 * np.sqrt(l / g) * ellipk(k2)


def period_measured(E):
    th0 = np.arccos(-E / (m * g * l))
    ev = lambda t, y: y[1]; ev.direction = 1          # p crosses zero upward: one full period after start? start has p=0 at max theta
    s = solve_ivp(lambda t, y: [y[1] / (m * l**2), -m * g * l * np.sin(y[0])], (0, 60), [th0, 0.0],
                  rtol=1e-11, atol=1e-13, events=ev)
    te = s.t_events[0]
    # p starts at 0 and goes negative; it crosses zero upward at each theta = -th0 turning point, i.e. at
    # T/2, 3T/2, ...; consecutive upward crossings are one period apart.
    return te[1] - te[0]


Es = np.array([-0.9, -0.5, 0.0, 0.5, 0.9])
print("Pendulum action-angle (m = l = g = 1):   E     J_numeric   J_closed    T=4K(k)   T_measured   dE/dJ*T")
for E in Es:
    Jn_, Jc = action_numeric(E), action_closed(E)
    assert np.isclose(Jn_, Jc, rtol=1e-7), (E, Jn_, Jc)
    T_c, T_m = period_closed(E), period_measured(E)
    assert np.isclose(T_c, T_m, rtol=1e-7), (E, T_c, T_m)
    dE = 1e-5
    dJdE = (action_closed(E + dE) - action_closed(E - dE)) / (2 * dE)
    assert np.isclose(dJdE, T_c, rtol=1e-6), (E, dJdE, T_c)          # nu = dE/dJ = 1/T  (3.25)
    print(f"   {E:5.2f}   {Jn_:9.5f}   {Jc:9.5f}   {T_c:8.5f}   {T_m:8.5f}    {T_c/dJdE:.7f}")
# §3.9 numbers at E = 0: J = 6.77770, T = 7.41630
assert np.isclose(action_closed(0.0), 6.77770, atol=5e-6) and np.isclose(period_closed(0.0), 7.41630, atol=5e-6)
# Small-amplitude limit: J -> 2 pi E_osc / omega with E_osc = E + mgl, omega = sqrt(g/l)
E_small = -m * g * l + 1e-4
assert np.isclose(action_closed(E_small), 2 * np.pi * (E_small + m * g * l) / np.sqrt(g / l), rtol=1e-4)
# Separatrix limit: J -> 16 m l sqrt(g l) (Example 4.3)
assert np.isclose(action_closed(m * g * l - 1e-12), 16 * m * l * np.sqrt(g * l), rtol=1e-6)
E_grid = np.linspace(-0.999, 0.999, 300)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
axs[0].plot(E_grid, [action_closed(E) for E in E_grid], label="J(E) closed form")
axs[0].plot(Es, [action_numeric(E) for E in Es], "o", label="numerical ∮p dθ")
axs[0].axhline(16, ls="--", color="gray", label="separatrix area 16")
axs[0].set_xlabel("E"); axs[0].set_ylabel("J"); axs[0].legend(); axs[0].grid(alpha=.3)
axs[0].set_title("Action of the pendulum")
axs[1].plot(E_grid, [period_closed(E) for E in E_grid], label="T = 4K(k) = dJ/dE")
axs[1].axhline(2 * np.pi, ls="--", color="gray", label="small-amplitude 2π")
axs[1].set_xlabel("E"); axs[1].set_ylabel("period"); axs[1].set_ylim(5, 25); axs[1].legend(); axs[1].grid(alpha=.3)
axs[1].set_title("Period diverges logarithmically at the separatrix")
fig.tight_layout(); fig.savefig(FIG / "action_angle_pendulum.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 9: Adiabatic invariance (§3.9, Example 4.7)
#
# Pendulum whose length shrinks from l0 to l0/2. Slowly (100 periods): E/omega
# is nearly constant, so the amplitude grows by 2^(3/4) = 1.682. Fast (half a
# period): E/omega is not conserved.

# %%
m, g, l0, th_init = 1.0, 1.0, 1.0, 0.10


def run_shrinking(tau):
    """Integrate the pendulum with l(t) going linearly from l0 to l0/2 over time tau, then holding."""
    l_of = lambda t: l0 - 0.5 * l0 * min(t / tau, 1.0)
    rhs = lambda t, yv: [yv[1] / (m * l_of(t)**2), -m * g * l_of(t) * np.sin(yv[0])]
    T_end = tau + 4 * 2 * np.pi * np.sqrt(l0 / 2 / g)             # a few periods after the shrink finishes
    s = solve_ivp(rhs, (0, T_end), [th_init, 0.0], rtol=1e-10, atol=1e-12, t_eval=np.linspace(0, T_end, 20000))
    lt = np.array([l_of(t) for t in s.t])
    E = s.y[1]**2 / (2 * m * lt**2) + m * g * lt * (1 - np.cos(s.y[0]))
    om = np.sqrt(g / lt)
    return s.t, s.y[0], E / om, lt


tau_slow = 100 * 2 * np.pi * np.sqrt(l0 / g)
tau_fast = 0.5 * 2 * np.pi * np.sqrt(l0 / g)
t_s, th_slow, I_slow, l_slow = run_shrinking(tau_slow)
t_f, th_fast, I_fast, l_fast = run_shrinking(tau_fast)
amp_final_slow = np.max(np.abs(th_slow[t_s > tau_slow]))
ratio_pred = 2**0.75
print(f"Example 4.7: slow shrink  -> amplitude ratio {amp_final_slow/th_init:.4f} (adiabatic prediction 2^(3/4) = {ratio_pred:.4f}); "
      f"max |E/omega - 1| deviation = {np.max(np.abs(I_slow/I_slow[0]-1)):.4f}")
amp_final_fast = np.max(np.abs(th_fast[t_f > tau_fast]))
print(f"             fast shrink  -> amplitude ratio {amp_final_fast/th_init:.4f}; E/omega changes by {abs(I_fast[-1]/I_fast[0]-1):.3f}")
assert abs(amp_final_slow / th_init - ratio_pred) / ratio_pred < 0.01
assert np.max(np.abs(I_slow / I_slow[0] - 1)) < 0.002, "E/omega should be adiabatically invariant (0.14%)"
assert abs(I_fast[-1] / I_fast[0] - 1) > 0.03, "fast variation should break the invariant"
# Example 4.7 numbers: slow ratio 1.6823, fast ratio 1.5669, fast E/omega off by 13.3%
assert np.isclose(amp_final_slow / th_init, 1.6823, atol=1e-4) and np.isclose(amp_final_fast / th_init, 1.5669, atol=1e-4)
assert np.isclose(abs(I_fast[-1] / I_fast[0] - 1), 0.133, atol=5e-4)
fig, axs = plt.subplots(1, 2, figsize=(11, 4))
for ax, (t_, I_, l_, lab) in zip(axs, [(t_s, I_slow, l_slow, "slow: 100 periods"), (t_f, I_fast, l_fast, "fast: 0.5 period")]):
    ax.plot(t_, I_ / I_[0], label=r"$E/\omega$ (normalized)")
    ax.plot(t_, l_ / l0, "--", label=r"$l(t)/l_0$")
    ax.set_title(f"Pendulum length halves, {lab}"); ax.set_xlabel("t"); ax.grid(alpha=.3); ax.legend()
axs[0].set_ylim(0.4, 1.1); axs[1].set_ylim(0.4, 1.6)
fig.tight_layout(); fig.savefig(FIG / "adiabatic_invariant.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 10: Hamilton–Jacobi for the oscillator (§3.10, Example 4.6)
#
# S(q, E, t) = W(q, E) - E t with W = ∫ sqrt(2mE - m^2 w^2 q^2) dq
#   = (q/2) sqrt(2mE - m^2 w^2 q^2) + (E/w) arcsin(q sqrt(m w^2 / 2E)).
# Checks: (i) S satisfies H(q, dS/dq) + dS/dt = 0 symbolically; (ii) along the
# true trajectory dS/dq = p(t) and Q = dS/dE is constant.

# %%
t_s, E_s = sp.symbols("t E", positive=True)
W = q / 2 * sp.sqrt(2 * m_s * E_s - m_s**2 * w_s**2 * q**2) + E_s / w_s * sp.asin(q * sp.sqrt(m_s * w_s**2 / (2 * E_s)))
S = W - E_s * t_s
assert sp.simplify(sp.diff(W, q) - sp.sqrt(2 * m_s * E_s - m_s**2 * w_s**2 * q**2)) == 0     # dW/dq = p(q, E)
HJ_residual = sp.simplify(sp.diff(S, q)**2 / (2 * m_s) + m_s * w_s**2 * q**2 / 2 + sp.diff(S, t_s))
assert HJ_residual == 0, HJ_residual
print("Hamilton–Jacobi: S = W - Et satisfies (1/2m)(dS/dq)^2 + (1/2) m w^2 q^2 + dS/dt = 0 symbolically.")
# numerical check along a trajectory: m = 1.3, omega = 0.7, E = 2.0, on the branch p > 0
mv, wv, Ev = 1.3, 0.7, 2.0
S_num = sp.lambdify((q, t_s), S.subs({m_s: mv, w_s: wv, E_s: Ev}), "numpy")
dSdq = sp.lambdify((q, t_s), sp.diff(S, q).subs({m_s: mv, w_s: wv, E_s: Ev}), "numpy")
dSdE = sp.lambdify((q, t_s), sp.diff(S, E_s).subs({m_s: mv, w_s: wv, E_s: Ev}), "numpy")
A_amp = np.sqrt(2 * Ev / (mv * wv**2)); beta = 0.3
tt = np.linspace(-beta + 0.02, -beta + np.pi / (2 * wv) - 0.02, 2000)      # quarter period, q: 0 -> A, p > 0
q_t = A_amp * np.sin(wv * (tt + beta)); p_t = np.sqrt(2 * mv * Ev) * np.cos(wv * (tt + beta))
assert np.allclose(dSdq(q_t, tt), p_t, rtol=1e-9), "dS/dq must equal p along the motion"
Qc = dSdE(q_t, tt)
assert np.allclose(Qc, Qc[0], atol=1e-9) and np.isclose(Qc[0], beta, atol=1e-9), "Q = dS/dE is the constant beta"
# S equals the action integral of L = p qdot - H along the path (up to a constant)
Lag = p_t**2 / (2 * mv) - 0.5 * mv * wv**2 * q_t**2
S_path = np.concatenate([[0], np.cumsum(0.5 * (Lag[1:] + Lag[:-1]) * np.diff(tt))])
assert np.allclose(S_num(q_t, tt) - S_num(q_t[0], tt[0]), S_path, atol=1e-5), np.max(np.abs(S_num(q_t, tt) - S_num(q_t[0], tt[0]) - S_path))
print(f"  along the trajectory: dS/dq = p to 1e-9, Q = dS/dE = beta = {Qc[0]:.6f} constant, and S - S(0) = ∫L dt to 1e-5.")
fig, ax = plt.subplots(figsize=(6, 4.5))
qq = np.linspace(-A_amp * .999, A_amp * .999, 200)
for tk in (0.0, 0.5, 1.0, 1.5):
    ax.plot(qq, S_num(qq, tk), label=f"t = {tk}")
ax.set_xlabel("q"); ax.set_ylabel("S(q, t)"); ax.set_title("Hamilton's principal function for the oscillator (p > 0 branch)")
ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "hamilton_jacobi_S.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 11: WKB — the classical density 1/p vs the exact n = 20 oscillator state (§3.10)
#
# With m = omega = hbar = 1 and E = n + 1/2, the WKB amplitude gives
# |psi|^2 ∝ 1/p(q) = 1/sqrt(A^2 - q^2), normalized to rho_cl = 1/(pi sqrt(A^2 - q^2)).
# The exact |psi_n|^2 oscillates about rho_cl. Probability in |q| < A/2:
# classical 1/3 exactly; quantum 0.3436 (integrated with quad between the exact
# limits +-A/2, not on the grid, whose points miss them).

# %%
n_q = 20
E_q = n_q + 0.5; A_q = np.sqrt(2 * E_q)
xq = np.linspace(-9, 9, 6001)
log_norm = -0.5 * (n_q * np.log(2) + gammaln(n_q + 1) + 0.5 * np.log(np.pi))
psi = np.exp(log_norm) * eval_hermite(n_q, xq) * np.exp(-xq**2 / 2)
rho_q = psi**2
assert np.isclose(np.trapezoid(rho_q, xq), 1.0, atol=1e-8)
rho_cl = np.where(np.abs(xq) < A_q, 1 / (np.pi * np.sqrt(np.maximum(A_q**2 - xq**2, 1e-300))), 0.0)
P_cl = 2 / np.pi * np.arcsin(0.5)                      # = 1/3 exactly
psi_sq = lambda x_: (np.exp(log_norm) * eval_hermite(n_q, x_) * np.exp(-x_**2 / 2))**2
P_q = quad(psi_sq, -A_q / 2, A_q / 2, limit=400, epsabs=1e-13)[0]
assert np.isclose(P_cl, 1 / 3)
assert abs(P_q - P_cl) < 0.02, (P_q, P_cl)
assert np.isclose(P_q, 0.343582, atol=2e-6), P_q                       # §3.10: 0.3436
print(f"§3.10 (WKB): P(|q| < A/2): classical 1/p density = {P_cl:.4f}, exact n = 20 state = {P_q:.4f}")
fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(xq, rho_q, lw=.9, label=rf"$|\psi_{{{n_q}}}(q)|^2$ exact")
ax.plot(xq, rho_cl, "C3", lw=2, label=r"classical $\rho \propto 1/p(q)$")
ax.axvline(A_q, color="gray", ls="--"); ax.axvline(-A_q, color="gray", ls="--")
ax.set_ylim(0, 0.35); ax.set_xlabel("q"); ax.set_ylabel("probability density")
ax.set_title("WKB: the quantum density oscillates about the classical dwell-time density")
ax.legend(); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "wkb_oscillator.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 12: Symplectic Euler as two shears; its exactly conserved modified Hamiltonian (§3.11)
#
# Kick K: (q, p) -> (q, p - h m w^2 q), drift D: (q, p) -> (q + h p/m, p). Both have
# det 1. Symplectic Euler is D∘K with matrix [[1 - h^2 w^2, h/m], [-h m w^2, 1]].
# It conserves H~ = p^2/2m + m w^2 q^2/2 - (h/2) w^2 q p exactly (to machine
# precision numerically), so its orbits are the ellipses H~ = const. (The map is
# the exact time-h flow of c(h) H~ with c = 1 + O(h^2), not of H~ itself.) Explicit Euler has det = 1 + h^2 w^2 and
# inflates area by that factor every step (Exercise B2: (1.01)^1000 = 20959).

# %%
h_s = sp.symbols("h", positive=True)
K_shear = sp.Matrix([[1, 0], [-h_s * m_s * w_s**2, 1]])
D_shear = sp.Matrix([[1, h_s / m_s], [0, 1]])
M_SE = sp.simplify(D_shear * K_shear)
assert M_SE == sp.Matrix([[1 - h_s**2 * w_s**2, h_s / m_s], [-h_s * m_s * w_s**2, 1]])
assert K_shear.det() == 1 and D_shear.det() == 1 and sp.simplify(M_SE.det()) == 1
# modified Hamiltonian, exactly conserved by the map
H_tilde = p**2 / (2 * m_s) + m_s * w_s**2 * q**2 / 2 - h_s / 2 * w_s**2 * q * p
q1, p1 = M_SE * sp.Matrix([q, p])
assert sp.simplify(H_tilde.subs({q: q1, p: p1}, simultaneous=True) - H_tilde) == 0
print("Symplectic Euler = drift ∘ kick (two shears, det 1 each); H~ = H - (h/2) w^2 q p is exactly conserved (symbolic).")
# numerical: m = 1.3, omega = 0.7, h = 0.1, 10^4 steps
mv, wv, hv = 1.3, 0.7, 0.1
zq, zp = 1.0, 0.3
Ht = lambda qq_, pp_: pp_**2 / (2 * mv) + mv * wv**2 * qq_**2 / 2 - hv / 2 * wv**2 * qq_ * pp_
Hs = []; Hts = []
for _ in range(10_000):
    zp = zp - hv * mv * wv**2 * zq; zq = zq + hv * zp / mv
    Hs.append(H_osc(zq, zp, mv, wv)); Hts.append(Ht(zq, zp))
Hs, Hts = np.array(Hs), np.array(Hts)
assert np.max(np.abs(Hts - Hts[0])) < 1e-12, "modified Hamiltonian must be exactly conserved"
assert np.max(np.abs(Hs - Hs[0])) / Hs[0] < hv * wv * 1.0 and np.max(np.abs(Hs - Hs[0])) > 1e-4   # bounded, O(h), not zero
print(f"  numerically: H~ varies by {np.max(np.abs(Hts-Hts[0])):.1e}; H oscillates within {np.max(np.abs(Hs-Hs[0]))/Hs[0]:.3%} (bounded, never drifts)")
# explicit Euler area growth: det M^N = (1 + h^2 w^2)^N
M_EE = np.array([[1, 0.1], [-0.1, 1]])          # h*omega = 0.1 with m = omega = 1
det_growth = np.linalg.det(np.linalg.matrix_power(M_EE, 1000))
assert np.isclose(det_growth, 1.01**1000, rtol=1e-8)
print(f"  explicit Euler, h w = 0.1: area factor after 1000 steps = (1.01)^1000 = {det_growth:.1f}  (Exercise B2)")
# Exercise B6: a Gaussian blob under free flow — covariance determinant is constant
sq, sp_, mm, tt_ = 0.7, 0.4, 1.0, 3.0
Sig0 = np.diag([sq**2, sp_**2]); Mfree = np.array([[1, tt_ / mm], [0, 1]])
Sig_t = Mfree @ Sig0 @ Mfree.T
assert np.isclose(Sig_t[0, 0], sq**2 + sp_**2 * tt_**2 / mm**2) and np.isclose(np.linalg.det(Sig_t), np.linalg.det(Sig0))
assert np.isclose(Sig_t[0, 0], 1.93) and np.isclose(Sig_t[0, 1], 0.48) and np.isclose(np.linalg.det(Sig_t), 0.0784)   # B6 numbers
print(f"  B6: sigma_q^2(t=3) = {Sig_t[0, 0]:.4f} (from {sq**2:.2f}), cov(q,p) = {Sig_t[0, 1]:.4f}, "
      f"det Sigma = {np.linalg.det(Sig_t):.4f} = {np.linalg.det(Sig0):.4f} at t = 0")

# %% [markdown]
# ## Part 13: Checks for the exercise solutions
#
# Every number quoted in solutions.md is computed here (B2, B3, B6 are in Parts 6 and 12).

# %%
v_, p_, c_, m_ = sp.symbols("v p c m", positive=True)
# A1(a): f = v^4/4 -> g = (3/4) p^(4/3);  A1(b): f = -m c^2 sqrt(1 - v^2/c^2) -> H = sqrt(p^2 c^2 + m^2 c^4)
assert sp.simplify(p_ * p_**sp.Rational(1, 3) - (p_**sp.Rational(1, 3))**4 / 4 - sp.Rational(3, 4) * p_**sp.Rational(4, 3)) == 0
v_of_p = p_ / (m_ * sp.sqrt(1 + p_**2 / (m_**2 * c_**2)))                     # inverts p = m v / sqrt(1 - v^2/c^2)
assert sp.simplify(m_ * v_of_p / sp.sqrt(1 - v_of_p**2 / c_**2) - p_) == 0
H_rel = sp.simplify(p_ * v_of_p + m_ * c_**2 * sp.sqrt(1 - v_of_p**2 / c_**2))
assert sp.simplify(H_rel - sp.sqrt(p_**2 * c_**2 + m_**2 * c_**4)) == 0
# A3: {q^2, p^2} = 4qp and {qp, H_osc} = p^2/m - m w^2 q^2 = 2(T - V)
assert poisson_bracket(q**2, p**2, [q], [p]) == 4 * q * p
assert sp.simplify(poisson_bracket(q * p, p**2 / (2 * m_s) + m_s * w_s**2 * q**2 / 2, [q], [p]) - (p**2 / m_s - m_s * w_s**2 * q**2)) == 0
# A4: determinants of five maps
al = sp.symbols("alpha", real=True)
maps_A4 = {"rotation": [q * sp.cos(al) - p * sp.sin(al), q * sp.sin(al) + p * sp.cos(al)], "scaling": [2 * q, p / 2],
           "Q=q^2,P=p/2q": [q**2, p / (2 * q)], "Q=q,P=p+q^2": [q, p + q**2], "Q=p,P=q": [p, q]}
dets_A4 = {k_: sp.simplify(sp.Matrix([[sp.diff(Z, zz) for zz in (q, p)] for Z in v]).det()) for k_, v in maps_A4.items()}
assert dets_A4 == {"rotation": 1, "scaling": 1, "Q=q^2,P=p/2q": 1, "Q=q,P=p+q^2": 1, "Q=p,P=q": -1}
# A5: V = x^3/3 - x: fixed points x = +-1, separatrix energy 2/3, loop spans [-1, 2]
xs = sp.symbols("x", real=True); V5 = xs**3 / 3 - xs
assert set(sp.solve(sp.diff(V5, xs), xs)) == {-1, 1} and V5.subs(xs, -1) == sp.Rational(2, 3) and V5.subs(xs, 1) == -sp.Rational(2, 3)
assert sp.factor(V5 - sp.Rational(2, 3)) == (xs - 2) * (xs + 1)**2 / 3
# A6: oscillator m = 2 kg, w = 3 rad/s, amplitude 0.5 m
E6 = 0.5 * 2 * 3**2 * 0.5**2; J6 = 2 * np.pi * E6 / 3; h_planck = 6.62607015e-34
assert np.isclose(E6, 2.25) and np.isclose(J6, 4.7124, atol=5e-5) and np.isclose(3 / (2 * np.pi), 0.4775, atol=5e-5)
assert np.isclose(J6 / h_planck, 7.112e33, rtol=1e-4)
print(f"A6: E = {E6:.2f} J, J = 2 pi E / w = {J6:.4f} J s, nu = {3/(2*np.pi):.4f} Hz, n = J/h = {J6/h_planck:.3e}")
# B1: Hamilton's equations for |p - qA|^2/2m + q phi give m dv/dt = q(E + v x B) (general A, phi)
tq, mq, qe = sp.symbols("t m q_e", positive=True)
X3, P3 = sp.symbols("x1:4", real=True), sp.symbols("P1:4", real=True)
A3f = [sp.Function(f"A{i}")(*X3, tq) for i in range(3)]; phif = sp.Function("phi")(*X3, tq)
Hq = sum((P3[i] - qe * A3f[i])**2 for i in range(3)) / (2 * mq) + qe * phif
vel = [sp.diff(Hq, P3[i]) for i in range(3)]
Ef = [-sp.diff(phif, X3[i]) - sp.diff(A3f[i], tq) for i in range(3)]
Bf = [sp.diff(A3f[2], X3[1]) - sp.diff(A3f[1], X3[2]), sp.diff(A3f[0], X3[2]) - sp.diff(A3f[2], X3[0]),
      sp.diff(A3f[1], X3[0]) - sp.diff(A3f[0], X3[1])]
for i in range(3):
    j, k_ = (i + 1) % 3, (i + 2) % 3
    dmv = -sp.diff(Hq, X3[i]) - qe * (sum(sp.diff(A3f[i], X3[jj]) * vel[jj] for jj in range(3)) + sp.diff(A3f[i], tq))
    assert sp.simplify(dmv - qe * (Ef[i] + vel[j] * Bf[k_] - vel[k_] * Bf[j])) == 0
print("B1: Hamilton's equations for (3.9) reproduce the Lorentz force (symbolic, arbitrary A and phi).")
# B3: A^2 + B^2 + L^2 = H^2 for the 2D isotropic oscillator
assert sp.expand(A2**2 + B2**2 + L2d**2 - H2**2) == 0
# B4: bead on a hoop rotating at Omega: H = p^2/(2mR^2) - m R^2 Om^2 sin^2/2 - m g R cos; H - (T + V) = -m R^2 Om^2 sin^2
R_, Om, gq, thq, thd = sp.symbols("R Omega g theta thetadot", positive=True)
Lh = mq * R_**2 * (thd**2 + Om**2 * sp.sin(thq)**2) / 2 + mq * gq * R_ * sp.cos(thq)
ph = sp.diff(Lh, thd); Hh = sp.simplify(ph * thd - Lh)
TV = mq * R_**2 * (thd**2 + Om**2 * sp.sin(thq)**2) / 2 - mq * gq * R_ * sp.cos(thq)
assert sp.simplify(Hh - TV + mq * R_**2 * Om**2 * sp.sin(thq)**2) == 0
Vh = -mq * R_**2 * Om**2 * sp.sin(thq)**2 / 2 - mq * gq * R_ * sp.cos(thq)
assert sp.simplify(sp.diff(Vh, thq) - mq * R_ * sp.sin(thq) * (gq - R_ * Om**2 * sp.cos(thq))) == 0
assert sp.simplify(sp.diff(Vh, thq, 2).subs(thq, 0) - mq * R_ * (gq - R_ * Om**2)) == 0
# B5: F1 = (1/2) m w q^2 cot Q generates (3.21)
F1 = m_s * w_s * q**2 / 2 * sp.cot(Q)
p_F1, P_F1 = sp.diff(F1, q), -sp.diff(F1, Q)
assert sp.simplify(P_F1.subs(q, q_of) - P) == 0 and sp.simplify(p_F1.subs(q, q_of) - p_of) == 0
# B7: free particle. S = m (q - q0)^2 / 2t solves H-J and equals the action of the straight path.
q0s, tt7 = sp.symbols("q0 t", positive=True)
S7 = m_s * (q - q0s)**2 / (2 * tt7)
assert sp.simplify(sp.diff(S7, q)**2 / (2 * m_s) + sp.diff(S7, tt7)) == 0
tau_ = sp.symbols("tau", positive=True)
assert sp.simplify(sp.integrate(m_s / 2 * ((q - q0s) / tt7)**2, (tau_, 0, tt7)) - S7) == 0   # L is constant on the path
print("A1, A3, A4, A5, B3, B4, B5, B7 symbolic checks passed.")


# %%
def energy_errors(steps, periods):
    """Return (RK4 final, Verlet max) relative energy errors for the theta0 = 2 pendulum (Exercise B8)."""
    hh = T_pend / steps; n = periods * steps
    e_r = rk4_pendulum(th_amp, 0.0, hh, n, n)[-1]
    e_v = verlet_pendulum(th_amp, 0.0, hh, n, 1)
    return abs(e_r - E0) / abs(E0), np.max(np.abs(e_v - E0)) / abs(E0)


(r50, v50), (r100, v100) = energy_errors(50, 1000), energy_errors(100, 1000)
print(f"B8: 10^3 periods. RK4 final error {r50:.3e} -> {r100:.3e} (ratio {r50/r100:.1f}); "
      f"Verlet max error {v50:.3e} -> {v100:.3e} (ratio {v50/v100:.2f})")
assert 3.8 < v50 / v100 < 4.2 and 25 < r50 / r100 < 40
# solutions B8 numbers, and Verlet's max over 10^3 periods equals Part 5's over 10^4 (no growth)
assert np.isclose(r50, 2.639e-2, rtol=5e-4) and np.isclose(r100, 8.347e-4, rtol=5e-4)
assert np.isclose(v50, 1.815e-2, rtol=5e-4) and np.isclose(v100, 4.533e-3, rtol=5e-4)
assert np.isclose(v50, err_ver.max(), rtol=1e-3)
yv = sp.symbols("y", real=True)
R4 = 1 + yv * sp.I + (yv * sp.I)**2 / 2 + (yv * sp.I)**3 / 6 + (yv * sp.I)**4 / 24
assert sp.expand(R4 * sp.conjugate(R4)) == sp.expand(1 - yv**6 / 72 + yv**8 / 576)      # |R(iy)|^2 for RK4

# C1: log divergence of the period, T ~ 2 sqrt(l/g) ln(32 m g l / (E_s - E))  (m = g = l = 1)
for Ec in (0.99, 0.999):
    T_ex = period_closed(Ec); T_as = 2 * np.log(32 / (1 - Ec))
    print(f"C1: E = {Ec}: T = 4K(k) = {T_ex:.4f}, asymptotic 2 ln(32/(1-E)) = {T_as:.4f}")
    assert abs(T_ex - T_as) < 0.02
    assert np.isclose(T_ex, {0.99: 16.1570, 0.999: 20.7491}[Ec], atol=5e-5)
    assert np.isclose(T_as, {0.99: 16.1418, 0.999: 20.7470}[Ec], atol=5e-5)
k2s = sp.symbols("k", positive=True)
Jk = sp.elliptic_e(k2s**2) - (1 - k2s**2) * sp.elliptic_k(k2s**2)              # sympy uses the parameter m = k^2
assert abs(float((sp.diff(Jk, k2s) - k2s * sp.elliptic_k(k2s**2)).subs(k2s, 0.6))) < 1e-12     # dJ/dk = 16 k K(k)

# C2: Kepler radial action and the degenerate frequencies (m = k = 1, E = -0.3, l = 0.8)
mk, kk, Ek, lk = 1.0, 1.0, -0.3, 0.8
r1, r2 = np.sort(np.roots([2 * mk * Ek, 2 * mk * kk, -lk**2]).real)
Jr_num = 2 * quad(lambda r: np.sqrt(max(2 * mk * (Ek + kk / r) - lk**2 / r**2, 0.0)), r1, r2, limit=200)[0]
Jr_cl = np.pi * kk * np.sqrt(2 * mk / abs(Ek)) - 2 * np.pi * lk
assert np.isclose(Jr_num, Jr_cl, rtol=1e-9)
T_kep = np.pi * kk * np.sqrt(2 * mk) / (2 * abs(Ek)**1.5)
ev_peri = lambda t, yy: yy[1]; ev_peri.direction = 1                           # p_r crosses 0 upward: perihelion
sk = solve_ivp(lambda t, yy: [yy[1] / mk, lk**2 / (mk * yy[0]**3) - kk / yy[0]**2, lk / (mk * yy[0]**2)],
               (0, 3 * T_kep), [r2, 0.0, 0.0], events=ev_peri, dense_output=True, rtol=1e-11, atol=1e-12)
t_peri = sk.t_events[0]; dphi = sk.sol(t_peri[1])[2] - sk.sol(t_peri[0])[2]
assert np.isclose(t_peri[1] - t_peri[0], T_kep, rtol=1e-7) and np.isclose(dphi, 2 * np.pi, rtol=1e-7)
assert np.allclose([r1, r2, Jr_cl, T_kep], [0.35857, 2.97476, 3.08501, 13.5193], atol=5e-5)   # solutions C2(c)
print(f"C2: r in [{r1:.5f}, {r2:.5f}], J_r numeric {Jr_num:.5f} = closed form {Jr_cl:.5f}; "
      f"radial period {t_peri[1]-t_peri[0]:.4f} = 1/nu = {T_kep:.4f}; phi advances {dphi/np.pi:.6f} pi per radial period")

# C3: standard map p' = p + K sin(theta), theta' = theta + p'. Jacobian det = 1; KAM confinement vs diffusion.
thm, pm, Km = sp.symbols("theta p K", real=True)
Msm = sp.Matrix([[sp.diff(e, z) for z in (thm, pm)] for e in (thm + pm + Km * sp.sin(thm), pm + Km * sp.sin(thm))])
assert sp.simplify(Msm.det()) == 1
rng_sm = np.random.default_rng(1)
spread = {}
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, Kv in zip(axs, (0.5, 0.97, 1.5)):
    th_m = rng_sm.uniform(0, 2 * np.pi, 50); p_m = rng_sm.uniform(0, 2 * np.pi, 50); p_start = p_m.copy()
    pts_th, pts_p, mx = [], [], 0.0
    for _ in range(2000):
        p_m = p_m + Kv * np.sin(th_m); th_m = (th_m + p_m) % (2 * np.pi)
        mx = max(mx, np.max(np.abs(p_m - p_start))); pts_th.append(th_m); pts_p.append(p_m % (2 * np.pi))
    spread[Kv] = mx
    ax.scatter(np.concatenate(pts_th), np.concatenate(pts_p), s=.05, color="k")
    ax.set_title(f"standard map, K = {Kv}: max |p - p0| = {mx:.2f}"); ax.set_xlabel(r"$\theta$"); ax.set_ylabel("p mod 2π")
fig.tight_layout(); fig.savefig(FIG / "standard_map.png", dpi=110); plt.close(fig)
assert spread[0.5] < 2 * np.pi and spread[0.97] < 2 * np.pi and spread[1.5] > 4 * np.pi
assert np.allclose([spread[0.5], spread[0.97], spread[1.5]], [2.45, 4.13, 28.71], atol=5e-3)   # solutions C3(c)
print("C3: max |p_n - p_0| over 50 orbits x 2000 kicks: " + ", ".join(f"K = {k_}: {v:.2f}" for k_, v in spread.items())
      + f"  (2 pi = {2*np.pi:.2f}); Greene's literature threshold K_c = 0.9716")

# %% [markdown]
# ## Your turn
#
# 1. **Fourth-order symplectic integrator.** Implement Forest–Ruth / Yoshida:
#    with c1 = 1/(2(2 - 2^(1/3))), c2 = (1 - 2^(1/3)) c1, do kick(c1 h), drift(d1 h),
#    kick(c2 h), drift(d2 h), kick(c2 h), drift(d1 h), kick(c1 h) with
#    d1 = 1/(2 - 2^(1/3)), d2 = -2^(1/3) d1 (so the drift coefficients sum to 1).
#    Add it to Part 5 and compare the drift plot with RK4 (same order, same cost
#    per step). Then halve h and confirm the energy error falls by 16.
# 2. **Kepler action-angle variables.** For H = p_r^2/2m + p_phi^2/(2 m r^2) - k/r
#    compute J_r = ∮ p_r dr numerically for an orbit with E < 0 and compare with
#    the closed form J_r = pi k sqrt(2m/|E|) - 2 pi p_phi (Exercise C2). Verify that
#    E depends only on J_r + J_phi, so nu_r = nu_phi: the orbit closes.
# 3. **The standard map.** Iterate p' = p + K sin(theta), theta' = theta + p' (mod 2 pi)
#    for 50 random initial points and 2000 iterations each, for K = 0.5, 0.97, 1.5.
#    Scatter (theta, p mod 2 pi). Show analytically that the map has Jacobian
#    determinant 1 (it is symplectic Euler for a kicked rotor), and watch the KAM
#    curves break into islands and a chaotic sea as K grows (Exercise C3).
# 4. **Liouville with a density estimate.** Replace the boundary-polygon area in
#    Part 4 with a 2D histogram of the 5000 interior points; verify that the
#    maximum density is (approximately) invariant under the symplectic map and
#    falls under explicit Euler. Say why it falls rather than rises.

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T_START:.1f} s. Figures in {FIG}")
