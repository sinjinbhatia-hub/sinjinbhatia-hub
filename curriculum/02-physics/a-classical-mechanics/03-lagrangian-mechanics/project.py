# %% [markdown]
# # CM-3 Project: From a Lagrangian to Equations of Motion, Multipliers, and the Worked Systems
#
# This script builds one tool and uses it on every system in the notes. The tool is
# `lagrangian_to_eom`: give it a SymPy Lagrangian in generalized coordinates and it
# returns the accelerations (by extracting the mass matrix M of (3.14) and solving
# M qdd = F) and a `rhs(t, y)` function for `scipy.integrate.solve_ivp`. With it the
# script checks every number quoted in notes §4: the Atwood tension by the multiplier
# DAE (3.21), the rotating hoop's conserved energy function and its pitchfork, the
# double pendulum's hand-expanded equations (4.10), its normal modes and its Lyapunov
# exponent, the sliding wedge, the spherical pendulum, the cyclotron orbit, and the
# Coriolis/centrifugal forces from the rotating-frame Lagrangian (4.19). The last part
# holds reference solutions to the coded exercises (B2d, B4, B7c, C2d, C3d).
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
T_START = time.perf_counter()

import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
G = 9.81                       # m/s^2, the value used in every numerical check in the notes
t = sp.symbols("t")
half = sp.Rational(1, 2)
TIGHT = dict(method="DOP853", rtol=1e-11, atol=1e-11)   # integrator settings for conservation checks


def period_from_crossings(tt, x):
    """Mean period of an oscillation from the upward zero crossings of x(t) (linear interpolation)."""
    s = np.sign(x)
    idx = np.where((s[:-1] < 0) & (s[1:] > 0))[0]
    tc = tt[idx] - x[idx] * (tt[idx + 1] - tt[idx]) / (x[idx + 1] - x[idx])
    assert len(tc) >= 3, "need at least three upward crossings"
    return np.mean(np.diff(tc))


# %% [markdown]
# ## Part 1: The symbolic-to-numeric pipeline (notes §3.6, eq. (3.14))
#
# For L(q, qd, t), each Euler–Lagrange expression E_i = d/dt(∂L/∂qd_i) − ∂L/∂q_i − Q_i is
# *linear* in the accelerations qdd (notes §3.6): E = M(q, t) qdd − F(q, qd, t). So
# M_ij = ∂E_i/∂qdd_j is the mass matrix, F = −E|_{qdd=0}, and qdd = M⁻¹F.
#
# **Build it yourself first.** Write `lagrangian_to_eom` before reading this one. The
# only subtle step is replacing q(t), q'(t), q''(t) by plain symbols in the right
# order (highest derivative first) so that differentiation with respect to qdd works.
# Then write `dae_accelerations` as a single `np.linalg.solve` on the block matrix (3.21).


# %%
class EOM:
    """Result of lagrangian_to_eom. Fields: q, qd, qdd (plain symbols), M, F (eq. 3.14),
    acc (accelerations in the original q(t) notation), acc_s (in plain symbols), rhs(t, y)."""
    def __init__(self, **kw):
        self.__dict__.update(kw)


def _plain_symbols(qs):
    names = [q.func.__name__ for q in qs]
    return ([sp.Symbol(n) for n in names],
            [sp.Symbol(n + "_d") for n in names],
            [sp.Symbol(n + "_dd") for n in names])


def lagrangian_to_eom(L, qs, t, Q=None, params=None, build_rhs=True):
    """Euler–Lagrange equations (3.11) of L in the sp.Function coordinates qs, solved for
    the accelerations via the mass matrix M of (3.14). Q: list of generalized forces
    Q_i^nc (same order as qs); params: dict of symbol -> number substituted before the
    numeric rhs is built. Returns an EOM object; rhs(t, y) takes y = [q, qd]."""
    n = len(qs)
    Q = [0] * n if Q is None else list(Q)
    EL = [sp.diff(sp.diff(L, q.diff(t)), t) - sp.diff(L, q) - Q[i] for i, q in enumerate(qs)]
    q_s, qd_s, qdd_s = _plain_symbols(qs)
    sub1 = [(q.diff(t, 2), s) for q, s in zip(qs, qdd_s)]
    sub2 = [(q.diff(t), s) for q, s in zip(qs, qd_s)]
    sub3 = [(q, s) for q, s in zip(qs, q_s)]
    EL_s = [sp.expand(e.subs(sub1).subs(sub2).subs(sub3)) for e in EL]
    M = sp.Matrix(n, n, lambda i, j: sp.diff(EL_s[i], qdd_s[j]))
    assert not any(mij.has(*qdd_s) for mij in M), "Euler–Lagrange expressions must be linear in qdd"
    F = -sp.Matrix([e.subs({s: 0 for s in qdd_s}) for e in EL_s])
    acc_s = M.LUsolve(F)
    back = [(s, q) for s, q in zip(q_s, qs)] + [(s, q.diff(t)) for s, q in zip(qd_s, qs)]
    acc = [a.subs(back) for a in acc_s]
    rhs = None
    if build_rhs:
        acc_num = acc_s.subs(params) if params else acc_s
        f = sp.lambdify([t, *q_s, *qd_s], list(acc_num), modules="math", cse=True)

        def rhs(tt, y):
            return np.concatenate([y[n:], f(tt, *y[:n], *y[n:])])
    return EOM(qs=qs, q=q_s, qd=qd_s, qdd=qdd_s, EL=EL, EL_s=EL_s, M=M, F=F, acc_s=acc_s, acc=acc, rhs=rhs)


def energy_function(L, qs, t, params=None):
    """The energy function h = Σ qd_i ∂L/∂qd_i − L of (3.15). Returns (h symbolic, h_num(t, y))
    where h_num accepts y of shape (2n,) or (2n, N)."""
    n = len(qs)
    h = sum(q.diff(t) * sp.diff(L, q.diff(t)) for q in qs) - L
    q_s, qd_s, _ = _plain_symbols(qs)
    h_s = h.subs([(q.diff(t), s) for q, s in zip(qs, qd_s)]).subs([(q, s) for q, s in zip(qs, q_s)])
    if params:
        h_s = h_s.subs(params)
    f = sp.lambdify([t, *q_s, *qd_s], h_s, modules="numpy", cse=True)
    return h, (lambda tt, y: f(tt, *y[:n], *y[n:]))


def dae_accelerations(M, F, J, gamma):
    """Solve the KKT system (3.21): [[M, -J^T], [J, 0]] [qdd; lam] = [F; gamma], where
    gamma_k = -qd^T (∇²f_k) qd (plus any stabilization terms). Returns (qdd, lam)."""
    M = np.atleast_2d(np.asarray(M, float)); J = np.atleast_2d(np.asarray(J, float))
    F = np.atleast_1d(np.asarray(F, float)); gamma = np.atleast_1d(np.asarray(gamma, float))
    n, m = M.shape[0], J.shape[0]
    K = np.block([[M, -J.T], [J, np.zeros((m, m))]])
    sol = np.linalg.solve(K, np.concatenate([F, gamma]))
    return sol[:n], sol[n:]


# %% [markdown]
# ## Part 2: Simple pendulum (Example 4.1) — calibration

# %%
m, ell, g = sp.symbols("m ell g", positive=True)
th = sp.Function("theta")(t)
L_pend = half * m * ell**2 * th.diff(t)**2 + m * g * ell * sp.cos(th)
pend = lagrangian_to_eom(L_pend, [th], t, params={m: 1, ell: 1, g: G})
assert sp.simplify(pend.acc[0] + g / ell * sp.sin(th)) == 0          # (4.1)
assert pend.M[0, 0] == m * ell**2                                      # M = moment of inertia
print("Example 4.1: pipeline gives theta_dd =", pend.acc[0])

# %% [markdown]
# ## Part 3: Atwood machine (Example 4.2) — reduced Lagrangian, then the multiplier DAE
#
# Reduced: L = ½(m1+m2) xd² + (m1−m2) g x gives (4.2). Multipliers: keep x1, x2 with
# f = x1 + x2 − ℓ, so J = (1, 1), ∇²f = 0, and (3.21) is solved at every step. The
# recorded λ must equal −2 m1 m2 g/(m1+m2) = −23.544 N for m1 = 3 kg, m2 = 2 kg.

# %%
m1, m2 = sp.symbols("m1 m2", positive=True)
x = sp.Function("x")(t)
L_atw = half * (m1 + m2) * x.diff(t)**2 + (m1 - m2) * g * x
atw = lagrangian_to_eom(L_atw, [x], t, build_rhs=False)
assert sp.simplify(atw.acc[0] - (m1 - m2) / (m1 + m2) * g) == 0      # (4.2)

M1, M2, ELL = 3.0, 2.0, 2.0
M_atw = np.diag([M1, M2]); F_atw = np.array([M1 * G, M2 * G]); J_atw = np.array([[1.0, 1.0]])


def atwood_dae_rhs(tt, y):
    """y = (x1, x2, x1d, x2d); the KKT solve of (3.21) gives the accelerations."""
    qdd, lam = dae_accelerations(M_atw, F_atw, J_atw, [0.0])    # ∇²f = 0, so gamma = 0
    return np.concatenate([y[2:], qdd])


y0 = np.array([1.0, ELL - 1.0, 0.0, 0.0])
sol = solve_ivp(atwood_dae_rhs, (0, 1.0), y0, t_eval=np.linspace(0, 1, 201), **TIGHT)
lam_t = np.array([dae_accelerations(M_atw, F_atw, J_atw, [0.0])[1][0] for _ in sol.t])
lam_exact = -2 * M1 * M2 * G / (M1 + M2)
acc_exact = (M1 - M2) / (M1 + M2) * G
assert np.isclose(lam_exact, -23.544) and np.isclose(acc_exact, 1.962)
assert np.allclose(lam_t, lam_exact, atol=1e-9), "lambda must equal -T at every step"
assert np.allclose(sol.y[0], 1.0 + 0.5 * acc_exact * sol.t**2, atol=1e-9)
drift_atwood = np.abs(sol.y[0] + sol.y[1] - ELL).max()
print(f"Example 4.2: x1_dd = {acc_exact:.4f} m/s^2, lambda = {lam_t[0]:.4f} N = -T, "
      f"T = {-lam_exact:.3f} N; constraint drift over 1 s = {drift_atwood:.1e} m (linear constraint)")

fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(sol.t, sol.y[0], label="$x_1$"); ax[0].plot(sol.t, sol.y[1], label="$x_2$")
ax[0].set_xlabel("t [s]"); ax[0].set_ylabel("m"); ax[0].legend(); ax[0].set_title("Atwood by the DAE (3.21)")
ax[1].plot(sol.t, lam_t, label=r"$\lambda(t)$ from the KKT solve")
ax[1].axhline(lam_exact, color="k", ls="--", label=r"$-2m_1m_2g/(m_1+m_2)$")
ax[1].set_xlabel("t [s]"); ax[1].set_ylabel("N"); ax[1].legend(); ax[1].set_title("multiplier = minus the tension")
fig.tight_layout(); fig.savefig(FIG / "atwood_dae.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 4: Bead on a rotating hoop (Example 4.3) — the showpiece
#
# (4.6) symbolically; h = T2 − T0 + V conserved while T + V is not; the small-oscillation
# period about θ* at ω = 2ω0 (notes: 0.7325 s); the pitchfork diagram with integrated
# trajectories confirming stability on each branch.

# %%
R, w = sp.symbols("R omega", positive=True)
L_hoop = half * m * R**2 * th.diff(t)**2 + half * m * R**2 * w**2 * sp.sin(th)**2 + m * g * R * sp.cos(th)
Rv = 0.5
w0 = np.sqrt(G / Rv)
hoop_params = {m: 1, R: Rv, g: G}
hoop = lagrangian_to_eom(L_hoop, [th], t, params={**hoop_params, w: 2 * w0})
assert sp.simplify(hoop.acc[0] - sp.sin(th) * (w**2 * sp.cos(th) - g / R)) == 0   # (4.6)

h_sym, h_num = energy_function(L_hoop, [th], t, params={**hoop_params, w: 2 * w0})
assert sp.simplify(h_sym - (half * m * R**2 * th.diff(t)**2 - half * m * R**2 * w**2 * sp.sin(th)**2
                            - m * g * R * sp.cos(th))) == 0                      # h = T2 - T0 + V

# numbers quoted in the notes
W = 2 * w0
th_star = np.arccos((w0 / W)**2)
Omega = W * np.sin(th_star)
assert np.isclose(w0, 4.4294, atol=5e-5) and np.isclose(W, 8.8589, atol=5e-5)
assert np.isclose(th_star, 1.3181, atol=5e-5) and np.isclose(np.degrees(th_star), 75.52, atol=5e-3)
assert np.isclose(np.sin(th_star), np.sqrt(15) / 4) and np.isclose(np.sin(th_star), 0.96825, atol=5e-6)
assert np.isclose(Omega, np.sqrt(15) / 2 * w0) and np.isclose(Omega, 8.5776, atol=5e-5)
assert np.isclose(2 * np.pi / Omega, 0.7325, atol=5e-5) and np.isclose(2 * np.pi / w0, 1.4185, atol=5e-5)
assert np.isclose(Omega**2, W**2 - w0**4 / W**2)                                   # (4.8)
print(f"Example 4.3: omega0 = {w0:.4f}, omega = {W:.4f}, theta* = {th_star:.4f} rad = {np.degrees(th_star):.2f} deg, "
      f"Omega = {Omega:.4f} rad/s, period = {2*np.pi/Omega:.4f} s (pendulum period {2*np.pi/w0:.4f} s)")

# conservation: h yes, T + V no
y0 = np.array([th_star + 0.02, 0.0])
tt = np.linspace(0, 10, 4001)
sol = solve_ivp(hoop.rhs, (0, 10), y0, t_eval=tt, **TIGHT)
h_t = h_num(sol.t, sol.y)
T0_t = 0.5 * Rv**2 * W**2 * np.sin(sol.y[0])**2
E_t = h_t + 2 * T0_t                      # T + V = h + 2 T0
assert np.abs(h_t - h_t[0]).max() < 1e-9, np.abs(h_t - h_t[0]).max()
assert np.abs(E_t - E_t[0]).max() > 1e-3
print(f"   h conserved to {np.abs(h_t - h_t[0]).max():.1e} J; T+V varies by {np.ptp(E_t):.3e} J")

# period about theta*
P_meas = period_from_crossings(sol.t, sol.y[0] - th_star)
assert abs(P_meas / (2 * np.pi / Omega) - 1) < 5e-3, P_meas
print(f"   measured period about theta* = {P_meas:.4f} s vs 2pi/Omega = {2*np.pi/Omega:.4f} s "
      f"({100*abs(P_meas/(2*np.pi/Omega)-1):.3f}% off)")

fig, ax = plt.subplots(figsize=(6, 3.4))
ax.plot(sol.t, h_t - h_t[0], label=r"$h - h(0)$  (conserved)")
ax.plot(sol.t, E_t - E_t[0], label=r"$(T+V) - (T+V)(0)$  (not conserved)")
ax.set_xlabel("t [s]"); ax.set_ylabel("J"); ax.legend(); ax.set_title(r"Rotating hoop: $h = T_2 - T_0 + V$ vs $T + V$")
fig.tight_layout(); fig.savefig(FIG / "hoop_energy.png", dpi=120); plt.close(fig)


# bifurcation diagram + stability by integration
def hoop_rhs_factory(wv):
    return lagrangian_to_eom(L_hoop, [th], t, params={**hoop_params, w: wv}).rhs


def settles(rhs, th_eq, eps=0.02, T=8.0, tol=0.06):
    """Start eps away from th_eq; stable if the trajectory stays within tol of th_eq."""
    s = solve_ivp(rhs, (0, T), [th_eq + eps, 0.0], rtol=1e-9, atol=1e-11, t_eval=np.linspace(0, T, 801))
    dev = np.abs(s.y[0] - th_eq).max()
    return bool(dev < tol), s


ws = np.linspace(0.2, 3.0, 400) * w0
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(ws[ws < w0] / w0, 0 * ws[ws < w0], "C0", lw=2, label="stable")
ax.plot(ws[ws >= w0] / w0, 0 * ws[ws >= w0], "C3--", lw=2, label="unstable")
branch = np.arccos((w0 / ws[ws >= w0])**2)
ax.plot(ws[ws >= w0] / w0, branch, "C0", lw=2); ax.plot(ws[ws >= w0] / w0, -branch, "C0", lw=2)
ax.plot(ws / w0, np.pi + 0 * ws, "C3--", lw=2)
stab_results = {}
for ratio in (0.5, 1.5, 2.0, 3.0):
    wv = ratio * w0
    rhs = hoop_rhs_factory(wv)
    ok_bottom, s_b = settles(rhs, 0.0)
    stab_results[(ratio, "bottom")] = ok_bottom
    ax.plot(ratio, s_b.y[0].mean(), "ko" if ok_bottom else "kx", ms=7)
    ax.vlines(ratio, s_b.y[0].min(), s_b.y[0].max(), color="k", lw=1)
    if ratio > 1:
        ts = np.arccos((w0 / wv)**2)
        ok_off, s_o = settles(rhs, ts)
        stab_results[(ratio, "off-axis")] = ok_off
        ax.plot(ratio, s_o.y[0].mean(), "ko" if ok_off else "kx", ms=7)
        ax.vlines(ratio, s_o.y[0].min(), s_o.y[0].max(), color="k", lw=1)
    ok_top, _ = settles(rhs, np.pi)
    stab_results[(ratio, "top")] = ok_top
assert stab_results[(0.5, "bottom")] is True
assert all(stab_results[(r, "bottom")] is False for r in (1.5, 2.0, 3.0))
assert all(stab_results[(r, "off-axis")] is True for r in (1.5, 2.0, 3.0))
assert all(stab_results[(r, "top")] is False for r in (0.5, 1.5, 2.0, 3.0))
ax.plot([], [], "ko", label="integrated: stays near equilibrium"); ax.plot([], [], "kx", label="integrated: departs")
ax.set_xlabel(r"$\omega/\omega_0$"); ax.set_ylabel(r"equilibrium $\theta$ [rad]")
ax.set_title(r"Pitchfork at $\omega = \omega_0 = \sqrt{g/R}$ (Example 4.3)"); ax.legend(fontsize=8, loc="upper left")
ax.set_ylim(-1.8, 3.5); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "hoop_bifurcation.png", dpi=120); plt.close(fig)
print("   bifurcation diagram: bottom stable only below omega0, off-axis branches stable, top never:", stab_results)

# %% [markdown]
# ## Part 5: The double pendulum (Example 4.4)
#
# (a) SymPy's accelerations satisfy the hand-expanded (4.10a), (4.10b) identically.
# (b) Small oscillations: M0 and K as Hessians, `eigh(K, M0)` gives Ω² = (2 ∓ √2) g/ℓ,
#     mode ratios ±√2. (c) Two trajectories 1e-8 apart: log-divergence and Lyapunov fit.
# (d) Snapshot grid. (e) Energy conservation.

# %%
l1, l2 = sp.symbols("ell1 ell2", positive=True)
th1, th2 = sp.Function("theta1")(t), sp.Function("theta2")(t)
Delta = th1 - th2
T_dp = (half * (m1 + m2) * l1**2 * th1.diff(t)**2 + half * m2 * l2**2 * th2.diff(t)**2
        + m2 * l1 * l2 * th1.diff(t) * th2.diff(t) * sp.cos(Delta))
V_dp = -(m1 + m2) * g * l1 * sp.cos(th1) - m2 * g * l2 * sp.cos(th2)
L_dp = T_dp - V_dp                                                                   # (4.9)
dp_params = {m1: 1, m2: 1, l1: 1, l2: 1, g: G}
dp = lagrangian_to_eom(L_dp, [th1, th2], t, params=dp_params)

# mass matrix and its determinant (4.9 paragraph)
M_hand = sp.Matrix([[(m1 + m2) * l1**2, m2 * l1 * l2 * sp.cos(Delta)], [m2 * l1 * l2 * sp.cos(Delta), m2 * l2**2]])
M_back = dp.M.subs(list(zip(dp.q, dp.qs)))
assert sp.simplify(M_back - M_hand) == sp.zeros(2, 2)
assert sp.simplify(dp.M.det() - m2 * l1**2 * l2**2 * (m1 + m2 * sp.sin(dp.q[0] - dp.q[1])**2)) == 0

# (a) hand-expanded (4.10) holds for SymPy's accelerations
a1, a2 = dp.acc
eq_a = ((m1 + m2) * l1 * th1.diff(t, 2) + m2 * l2 * th2.diff(t, 2) * sp.cos(Delta)
        + m2 * l2 * th2.diff(t)**2 * sp.sin(Delta) + (m1 + m2) * g * sp.sin(th1))                   # (4.10a)
eq_b = (l2 * th2.diff(t, 2) + l1 * th1.diff(t, 2) * sp.cos(Delta)
        - l1 * th1.diff(t)**2 * sp.sin(Delta) + g * sp.sin(th2))                                    # (4.10b)
sub_acc = [(th1.diff(t, 2), a1), (th2.diff(t, 2), a2)]
assert sp.simplify(eq_a.subs(sub_acc)) == 0 and sp.simplify(eq_b.subs(sub_acc)) == 0
print("Example 4.4: SymPy accelerations satisfy (4.10a) and (4.10b) identically")
# held straight (theta1 = theta2 = vt): (4.10b) forces vt_dd = -g sin(vt)/(l1+l2); (4.10a) keeps a residual m1 g l2 sin(vt)/(l1+l2)
vt = sp.Function("vartheta")(t)
straight = [(th1.diff(t, 2), -g * sp.sin(vt) / (l1 + l2)), (th2.diff(t, 2), -g * sp.sin(vt) / (l1 + l2)), (th1, vt), (th2, vt)]
assert sp.simplify(eq_b.subs(straight)) == 0
assert sp.simplify(eq_a.subs(straight) - m1 * g * l2 * sp.sin(vt) / (l1 + l2)) == 0

# (b) normal modes: Hessians at the hanging equilibrium
qd_syms = [sp.Symbol("theta1_d"), sp.Symbol("theta2_d")]
T_s = T_dp.subs([(th1.diff(t), qd_syms[0]), (th2.diff(t), qd_syms[1])]).subs([(th1, 0), (th2, 0)])
M0_sym = sp.hessian(T_s, qd_syms)
K_sym = sp.hessian(V_dp, [th1, th2]).subs([(th1, 0), (th2, 0)])
unit = {m1: 1, m2: 1, l1: 1, l2: 1}
assert M0_sym.subs(unit) == sp.Matrix([[2, 1], [1, 1]]) and K_sym.subs(unit) == g * sp.Matrix([[2, 0], [0, 1]])
M0 = np.array(M0_sym.subs(dp_params), float); K = np.array(K_sym.subs(dp_params), float)
Om2, A = eigh(K, M0)
assert np.allclose(Om2, [(2 - np.sqrt(2)) * G, (2 + np.sqrt(2)) * G])
assert np.allclose(Om2, [5.7466, 33.4934], atol=5e-5) and np.allclose(np.sqrt(Om2), [2.3972, 5.7874], atol=5e-5)
ratios = A[1] / A[0]
assert np.allclose(ratios, [np.sqrt(2), -np.sqrt(2)])
print(f"   normal modes: Omega^2 = {Om2[0]:.4f}, {Om2[1]:.4f}  (Omega = {np.sqrt(Om2[0]):.4f}, {np.sqrt(Om2[1]):.4f} rad/s); "
      f"a2/a1 = {ratios[0]:.4f}, {ratios[1]:.4f}")

# (c) chaos: two trajectories 1e-8 apart
y0 = np.array([2.0, 2.0, 0.0, 0.0])
tt = np.arange(0, 30.0001, 0.01)
solA = solve_ivp(dp.rhs, (0, 30), y0, t_eval=tt, **TIGHT)
solB = solve_ivp(dp.rhs, (0, 30), y0 + np.array([1e-8, 0, 0, 0]), t_eval=tt, **TIGHT)
delta = np.linalg.norm(solA.y - solB.y, axis=0)
window = (delta > 1e-7) & (delta < 1e-2)
lyap, intercept = np.polyfit(tt[window], np.log(delta[window]), 1)
assert 0.9 < lyap < 1.4, lyap                      # notes: "about 1.1 s^-1"
assert delta[1500:].max() > 1.0                    # by 15 s the separation is order 1 or larger
t_sat = tt[np.argmax(delta > 1.0)]
print(f"   Lyapunov exponent estimate = {lyap:.3f} 1/s (fit over {tt[window][0]:.2f}-{tt[window][-1]:.2f} s); "
      f"one decimal digit lost every {np.log(10)/lyap:.2f} s; separation reaches 1 at t = {t_sat:.1f} s")

# (e) energy conservation on trajectory A (scleronomic: h = T + V)
_, h_dp = energy_function(L_dp, [th1, th2], t, params=dp_params)
E_dp = h_dp(solA.t, solA.y)
assert np.abs(E_dp - E_dp[0]).max() < 1e-7, np.abs(E_dp - E_dp[0]).max()
print(f"   energy conserved to {np.abs(E_dp - E_dp[0]).max():.1e} J over 30 s")

fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.semilogy(tt, delta, lw=1, label=r"$\delta(t)$ (state-space distance)")
ax.semilogy(tt[window], np.exp(intercept + lyap * tt[window]), "k--", label=rf"fit: slope $\lambda_L$ = {lyap:.2f} s$^{{-1}}$")
ax.set_xlabel("t [s]"); ax.set_ylabel(r"$\delta$"); ax.legend(); ax.set_title(r"Double pendulum: two starts $10^{-8}$ rad apart")
ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(FIG / "double_pendulum_lyapunov.png", dpi=120); plt.close(fig)


# (d) snapshot grid
def dp_positions(th1v, th2v, L1=1.0, L2=1.0):
    x1, y1 = L1 * np.sin(th1v), -L1 * np.cos(th1v)
    return (x1, y1), (x1 + L2 * np.sin(th2v), y1 - L2 * np.cos(th2v))


snap_times = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24]
fig, axes = plt.subplots(3, 4, figsize=(11, 8.2))
for ax, ts in zip(axes.flat, snap_times):
    k = int(round(ts / 0.01))
    for s, col, lab in ((solA, "C0", "A"), (solB, "C3", "B")):
        (x1, y1), (x2, y2) = dp_positions(s.y[0, k], s.y[1, k])
        ax.plot([0, x1, x2], [0, y1, y2], "-o", color=col, alpha=.8, label=lab)
    ax.set_xlim(-2.1, 2.1); ax.set_ylim(-2.1, 2.1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"t = {ts} s,  δ = {delta[k]:.1e}", fontsize=9)
axes[0, 0].legend(fontsize=8)
fig.suptitle("Double pendulum: trajectories A and B agree, then do not", y=0.995)
fig.tight_layout(); fig.savefig(FIG / "double_pendulum_snapshots.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 6: Pendulum with a horizontally oscillating support (Example 4.5)
# A rheonomic Lagrangian with explicit t. Only the symbolic check of (4.11).

# %%
a, gam = sp.symbols("a gamma", positive=True)
x_b = a * sp.cos(gam * t) + ell * sp.sin(th); y_b = -ell * sp.cos(th)
L_osc = half * m * (x_b.diff(t)**2 + y_b.diff(t)**2) + m * g * ell * sp.cos(th)
osc = lagrangian_to_eom(L_osc, [th], t, build_rhs=False)
assert sp.simplify(osc.acc[0] - (-(g / ell) * sp.sin(th) + a * gam**2 / ell * sp.cos(gam * t) * sp.cos(th))) == 0  # (4.11)
print("Example 4.5: (4.11) confirmed symbolically; T1 term =",
      sp.simplify(sp.expand(L_osc).coeff(th.diff(t), 1)), "* theta_d")

# %% [markdown]
# ## Part 7: Block on a sliding wedge (Example 4.6)

# %%
Mw, alpha = sp.symbols("M alpha", positive=True)
X, s_ = sp.Function("X")(t), sp.Function("s")(t)
xb = X + s_ * sp.cos(alpha); yb = -s_ * sp.sin(alpha)
L_wedge = half * Mw * X.diff(t)**2 + half * m * (xb.diff(t)**2 + yb.diff(t)**2) - m * g * yb
wedge_params = {Mw: 4, m: 1, alpha: sp.pi / 6, g: G}
wedge = lagrangian_to_eom(L_wedge, [X, s_], t, params=wedge_params)
s_dd_hand = (Mw + m) * g * sp.sin(alpha) / (Mw + m * sp.sin(alpha)**2)
X_dd_hand = -m * g * sp.sin(alpha) * sp.cos(alpha) / (Mw + m * sp.sin(alpha)**2)
assert sp.simplify(wedge.acc[1] - s_dd_hand) == 0 and sp.simplify(wedge.acc[0] - X_dd_hand) == 0   # (4.13)
M_wedge_hand = sp.Matrix([[Mw + m, m * sp.cos(alpha)], [m * sp.cos(alpha), m]])
assert sp.simplify(wedge.M - M_wedge_hand) == sp.zeros(2, 2)
s_dd = float(s_dd_hand.subs(wedge_params)); X_dd = float(X_dd_hand.subs(wedge_params))
assert np.isclose(s_dd, 5.7706, atol=5e-5) and np.isclose(X_dd, -0.9995, atol=5e-5)
sol = solve_ivp(wedge.rhs, (0, 1.0), [0, 0, 0, 0], t_eval=np.linspace(0, 1, 101), **TIGHT)
assert np.allclose(sol.y[1], 0.5 * s_dd * sol.t**2, atol=1e-9)          # s(t) exactly quadratic
assert np.allclose(sol.y[0], 0.5 * X_dd * sol.t**2, atol=1e-9)
pX = (4 + 1) * sol.y[2] + 1 * np.cos(np.pi / 6) * sol.y[3]                # (4.12)
assert np.abs(pX).max() < 1e-10
_, h_wedge = energy_function(L_wedge, [X, s_], t, params=wedge_params)
E_w = h_wedge(sol.t, sol.y)
assert np.abs(E_w - E_w[0]).max() < 1e-9
print(f"Example 4.6: s_dd = {s_dd:.4f} m/s^2, X_dd = {X_dd:.4f} m/s^2; s(t) quadratic, p_X = 0 conserved, "
      f"energy conserved to {np.abs(E_w - E_w[0]).max():.1e} J")

# %% [markdown]
# ## Part 8: Spherical pendulum (Example 4.7, Exercise B6)
# p_φ conserved (4.15); the θ equation is (4.16); the conical solution stays at θ0; the
# small oscillation of θ about θ0 has Ω = φ̇0 √(1 + 3cos²θ0) (Exercise B6).

# %%
ph = sp.Function("phi")(t)
L_sph = half * m * ell**2 * (th.diff(t)**2 + sp.sin(th)**2 * ph.diff(t)**2) + m * g * ell * sp.cos(th)   # (4.14)
sph_params = {m: 1, ell: 1, g: G}
sph = lagrangian_to_eom(L_sph, [th, ph], t, params=sph_params)
p_phi = sp.diff(L_sph, ph.diff(t))
assert sp.simplify(p_phi - m * ell**2 * sp.sin(th)**2 * ph.diff(t)) == 0            # (4.15)
# phi cyclic: d p_phi/dt = 0 along the motion
assert sp.simplify(p_phi.diff(t).subs([(th.diff(t, 2), sph.acc[0]), (ph.diff(t, 2), sph.acc[1])])) == 0
pphi = sp.symbols("p_phi")
th_dd_eff = pphi**2 * sp.cos(th) / (m**2 * ell**4 * sp.sin(th)**3) - g / ell * sp.sin(th)   # (4.16)/(m ell^2)
assert sp.simplify(sph.acc[0].subs(ph.diff(t), pphi / (m * ell**2 * sp.sin(th)**2)) - th_dd_eff) == 0

th0 = np.radians(40.0)
phid0 = np.sqrt(G / np.cos(th0))                        # conical pendulum, ell = 1
Om_B6 = phid0 * np.sqrt(1 + 3 * np.cos(th0)**2)
print(f"Example 4.7 / B6: theta0 = 40 deg: phi_dot0 = {phid0:.4f} rad/s, Omega = {Om_B6:.4f} rad/s, "
      f"Omega/phi_dot0 = {Om_B6/phid0:.4f}")
assert 1 < Om_B6 / phid0 < 2
# conical motion: theta stays constant
sol_c = solve_ivp(sph.rhs, (0, 5), [th0, 0.0, 0.0, phid0], t_eval=np.linspace(0, 5, 501), **TIGHT)
assert np.abs(sol_c.y[0] - th0).max() < 1e-9
# perturbed conical motion: measure the theta oscillation period
sol_p = solve_ivp(sph.rhs, (0, 12), [th0 + 0.01, 0.0, 0.0, phid0], t_eval=np.linspace(0, 12, 6001), **TIGHT)
P_th = period_from_crossings(sol_p.t, sol_p.y[0] - sol_p.y[0].mean())
assert abs(P_th / (2 * np.pi / Om_B6) - 1) < 1e-2, P_th
pphi_t = np.sin(sol_p.y[0])**2 * sol_p.y[3]
_, h_sph = energy_function(L_sph, [th, ph], t, params=sph_params)
E_sph = h_sph(sol_p.t, sol_p.y)
assert np.abs(pphi_t - pphi_t[0]).max() < 1e-10 and np.abs(E_sph - E_sph[0]).max() < 1e-9
print(f"   conical motion holds theta0 to {np.abs(sol_c.y[0]-th0).max():.1e}; perturbed: theta period {P_th:.4f} s vs "
      f"2pi/Omega = {2*np.pi/Om_B6:.4f} s; p_phi and E conserved to 1e-10, 1e-9")

fig, ax = plt.subplots(1, 2, figsize=(9, 4))
xs = np.sin(sol_p.y[0]) * np.cos(sol_p.y[1]); ys = np.sin(sol_p.y[0]) * np.sin(sol_p.y[1])
ax[0].plot(xs, ys, lw=.8); ax[0].set_aspect("equal"); ax[0].set_title("spherical pendulum, top view: precessing orbit")
ax[1].plot(sol_p.t, sol_p.y[0]); ax[1].axhline(th0, color="k", ls="--"); ax[1].set_xlabel("t [s]"); ax[1].set_ylabel(r"$\theta$ [rad]")
ax[1].set_title(rf"$\theta$ oscillates about $\theta_0$ with period {P_th:.3f} s")
fig.tight_layout(); fig.savefig(FIG / "spherical_pendulum.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 9: Charged particle in a uniform B (Example 4.8)
# Landau gauge A = (0, Bx, 0): L = ½m|v|² + q_e B x ẏ. p_y = m ẏ + q_e B x is conserved,
# m ẏ is not; radius m v/(q_e B) = 0.5, period 2π m/(q_e B) = π for q_e = m = 1, B = 2.

# %%
qe, B = sp.symbols("q_e B", positive=True)
xc, yc, zc = sp.Function("x")(t), sp.Function("y")(t), sp.Function("z")(t)
L_em = half * m * (xc.diff(t)**2 + yc.diff(t)**2 + zc.diff(t)**2) + qe * B * xc * yc.diff(t)
em_params = {m: 1, qe: 1, B: 2}
em = lagrangian_to_eom(L_em, [xc, yc, zc], t, params=em_params)
assert sp.simplify(em.acc[0] - qe * B / m * yc.diff(t)) == 0 and sp.simplify(em.acc[1] + qe * B / m * xc.diff(t)) == 0
assert em.acc[2] == 0
p_sym = [sp.diff(L_em, c.diff(t)) for c in (xc, yc, zc)]
assert p_sym[1] == m * yc.diff(t) + qe * B * xc                                   # (4.18), y-component
sol = solve_ivp(em.rhs, (0, np.pi), [0, 0, 0, 1, 0, 0], t_eval=np.linspace(0, np.pi, 629), **TIGHT)
py_t = 1 * sol.y[4] + 1 * 2 * sol.y[0]
assert np.abs(py_t - py_t[0]).max() < 1e-10 and np.ptp(sol.y[4]) > 1.9       # p_y conserved, m*ydot swings by 2v
radius = np.sqrt(sol.y[0]**2 + (sol.y[1] + 0.5)**2)     # centre at (0, -v/omega_c) = (0, -0.5) for this start
assert np.allclose(radius, 0.5, atol=1e-9)
assert np.allclose(sol.y[:3, -1], [0, 0, 0], atol=1e-9)  # back at the start after one period pi
print("Example 4.8: p_y = m*ydot + q_e B x conserved to 1e-10 while m*ydot swings by", f"{np.ptp(sol.y[4]):.3f};",
      f"orbit radius {radius.mean():.4f} = m v/(q_e B); period pi confirmed")
fig, ax = plt.subplots(figsize=(4, 4))
ax.plot(sol.y[0], sol.y[1]); ax.set_aspect("equal"); ax.set_title(r"cyclotron orbit, $\omega_c = q_eB/m$"); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "cyclotron.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 10: Non-inertial frames from the Lagrangian (Example 4.9)
# (a) For ω = ω(t) ẑ, the Euler–Lagrange equations of (4.19) are (4.20) *including* the
# Euler term −m ω̇ × r'. (b) A projectile integrated in the rotating frame (ω constant)
# agrees with the inertial parabola rotated by R(t)ᵀ.

# %%
xr, yr, zr = sp.Function("xp")(t), sp.Function("yp")(t), sp.Function("zp")(t)
wt = sp.Function("omega")(t)
rp = sp.Matrix([xr, yr, zr]); rpd = rp.diff(t)
omega_vec = sp.Matrix([0, 0, wt])
v_in = rpd + omega_vec.cross(rp)
L_rot = half * m * (v_in.dot(v_in)) - m * g * zr                                   # (4.19)
rot = lagrangian_to_eom(L_rot, [xr, yr, zr], t, build_rhs=False)
rhs_420 = (sp.Matrix([0, 0, -g]) - 2 * omega_vec.cross(rpd) - omega_vec.cross(omega_vec.cross(rp))
           - omega_vec.diff(t).cross(rp))                                          # (4.20) / m
assert all(sp.simplify(rot.acc[i] - rhs_420[i]) == 0 for i in range(3))
# the Euler term is really there: with omega constant the omega' terms vanish, otherwise not
assert any(rot.acc[i].has(wt.diff(t)) for i in range(3))
print("Example 4.9: Euler–Lagrange of (4.19) equals (4.20) componentwise, Euler term included")

# numeric overlay: omega = 1 rad/s about z, projectile with V = m g z
Wr = 1.0
rot_num = lagrangian_to_eom(L_rot.subs(wt, w), [xr, yr, zr], t, params={m: 1, g: G, w: Wr})
r0 = np.array([1.0, 0.0, 0.0]); v0 = np.array([0.0, 2.0, 5.0])        # inertial initial data
omega_np = np.array([0, 0, Wr])
y0 = np.concatenate([r0, v0 - np.cross(omega_np, r0)])                 # r'(0) = r(0), r'd(0) = v - omega x r'
tt = np.linspace(0, 1.0, 201)
sol = solve_ivp(rot_num.rhs, (0, 1.0), y0, t_eval=tt, **TIGHT)
r_inertial = r0[:, None] + v0[:, None] * tt + 0.5 * np.array([0, 0, -G])[:, None] * tt**2
c, s_ = np.cos(Wr * tt), np.sin(Wr * tt)
r_rotated = np.array([c * r_inertial[0] + s_ * r_inertial[1], -s_ * r_inertial[0] + c * r_inertial[1], r_inertial[2]])
err_rot = np.abs(sol.y[:3] - r_rotated).max()
assert err_rot < 1e-8, err_rot
print(f"   projectile: rotating-frame Lagrangian vs rotated inertial parabola agree to {err_rot:.1e} m")
fig, ax = plt.subplots(figsize=(5, 4.2))
ax.plot(r_inertial[0], r_inertial[1], "k--", label="inertial frame (parabola, seen from above)")
ax.plot(r_rotated[0], r_rotated[1], "C0", lw=3, alpha=.5, label=r"inertial solution rotated by $\mathcal{R}^\top$")
ax.plot(sol.y[0], sol.y[1], "C3", lw=1, label="integrated from the rotating-frame L (4.19)")
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect("equal"); ax.legend(fontsize=8); ax.grid(alpha=.3)
ax.set_title(r"Coriolis + centrifugal from $L$: $\omega = 1$ rad/s")
fig.tight_layout(); fig.savefig(FIG / "rotating_frame.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 11: Rayleigh dissipation on a damped pendulum (notes §3.9, eq. 3.23)
# Drag −c ṙ on the bob gives 𝓕 = ½ c ℓ² θ̇², Q_θ = −∂𝓕/∂θ̇ = −c ℓ² θ̇. Check dE/dt = −2𝓕.

# %%
c = sp.symbols("c", positive=True)
F_ray = half * c * ell**2 * th.diff(t)**2
Q_ray = [-sp.diff(F_ray, th.diff(t))]
ray_params = {m: 1, ell: 1, g: G, c: 0.3}
ray = lagrangian_to_eom(L_pend, [th], t, Q=Q_ray, params=ray_params)
assert sp.simplify(ray.acc[0] - (-(g / ell) * sp.sin(th) - c / m * th.diff(t))) == 0
_, h_ray = energy_function(L_pend, [th], t, params=ray_params)
tt = np.linspace(0, 10, 2001)
sol = solve_ivp(ray.rhs, (0, 10), [1.5, 0.0], t_eval=tt, **TIGHT)
E_r = h_ray(sol.t, sol.y)
dEdt = np.gradient(E_r, tt)
minus2F = -2 * (0.5 * 0.3 * 1**2 * sol.y[1]**2)
assert np.abs(dEdt - minus2F)[5:-5].max() < 2e-3 * np.abs(minus2F).max()
print(f"Rayleigh: dE/dt vs -2F agree to {np.abs(dEdt - minus2F)[5:-5].max():.1e} W "
      f"(finite-difference derivative); energy lost over 10 s = {E_r[0]-E_r[-1]:.4f} J")
fig, ax = plt.subplots(figsize=(6, 3.4))
ax.plot(tt, dEdt, label="dE/dt (finite difference)"); ax.plot(tt, minus2F, "--", label=r"$-2\mathcal{F}$")
ax.set_xlabel("t [s]"); ax.set_ylabel("W"); ax.legend(); ax.set_title("Damped pendulum: (3.23)")
fig.tight_layout(); fig.savefig(FIG / "rayleigh.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Your turn
#
# Do these before reading the reference section below.
#
# 1. **Slowly varying length (Exercise C3, adiabatic invariant).** ℓ(t) = 1 + 0.02 t for
#    0 ≤ t ≤ 50 s, θ(0) = 0.1 rad at rest. Derive ℓθ̈ + 2ℓ̇θ̇ + g sinθ = 0 (it is *not*
#    θ̈ = −(g/ℓ)sinθ), integrate, and plot E/ω with E = ½mℓ²θ̇² + ½mgℓθ² and ω = √(g/ℓ)
#    at the turning points. It should be flat to well under 1% while E itself drops by about 30%.
#    Also plot θ_max ℓ^{3/4}.
# 2. **Bead on a rotating parabolic wire (Exercise B4).** L = ½m(1 + k²r²)ṙ² + ½m(ω² − gk)r².
#    Build the rhs with `lagrangian_to_eom`, then integrate for ω² < gk, = gk, > gk. There is
#    no pitchfork: at ω² = gk every point of the wire is an equilibrium. Confirm the
#    small-oscillation frequency √(gk − ω²) about r = 0 below threshold, and that above it
#    the bead leaves.
# 3. **Rayleigh damping on the double pendulum (Exercise B7).** Write 𝓕 for drag −cṙ_a on
#    both bobs by replacing (m1+m2, m2) → (2c, c) in T of (4.9), pass Q_i = −∂𝓕/∂θ̇_i to
#    `lagrangian_to_eom`, and verify dE/dt = −2𝓕 along the trajectory from θ1 = θ2 = 2 rad.
# 4. **Pendulum tension as a DAE (Exercise B2d).** Use `dae_accelerations` with M = mI,
#    J = (2x, 2y), gamma = −2(ẋ² + ẏ²). Compare λ with −T/(2ℓ), measure the drift of
#    f = x² + y² − ℓ² over 10 s, then add Baumgarte stabilization with β = 5.
# 5. **Kapitza (Exercise C2d).** Integrate θ̈ + (g − aγ²cosγt) sinθ/ℓ = 0 from θ = π + 0.1
#    for γ = 150 and γ = 60 rad/s (a = 0.05 m, ℓ = 1 m). Which one stays inverted?

# %% [markdown]
# ## Reference solutions to the coded exercises (B2d, B4, B7c, C2d, C3d)
# Every number quoted in `solutions.md` for these exercises is printed here. The final cell
# checks the symbolic results of the pen-and-paper exercises (A3, A4, B1, B3, B5, B6, B8, C1–C3).

# %%
# --- B2(d): pendulum tension by the multiplier DAE, drift, and Baumgarte stabilization
ELL_P, M_P = 1.0, 1.0
M_pend = M_P * np.eye(2); F_pend = np.array([0.0, -M_P * G])


def pendulum_dae_rhs_factory(beta=0.0):
    def rhs(tt, y):
        xx, yy, xd, yd = y
        J = np.array([[2 * xx, 2 * yy]])
        f = xx**2 + yy**2 - ELL_P**2
        fd = 2 * (xx * xd + yy * yd)
        gamma = -2 * (xd**2 + yd**2) - 2 * beta * fd - beta**2 * f
        qdd, lam = dae_accelerations(M_pend, F_pend, J, [gamma])
        return np.concatenate([[xd, yd], qdd])
    return rhs


def pendulum_lambda(y):
    xx, yy, xd, yd = y
    _, lam = dae_accelerations(M_pend, F_pend, np.array([[2 * xx, 2 * yy]]), [-2 * (xd**2 + yd**2)])
    return lam[0]


th_p0 = 1.0
y0 = np.array([ELL_P * np.sin(th_p0), -ELL_P * np.cos(th_p0), 0.0, 0.0])
tt = np.linspace(0, 10, 2001)
drifts = {}
for beta in (0.0, 5.0):
    sol = solve_ivp(pendulum_dae_rhs_factory(beta), (0, 10), y0, t_eval=tt, rtol=1e-6, atol=1e-8)   # RK45, modest tolerance
    f_t = sol.y[0]**2 + sol.y[1]**2 - ELL_P**2
    drifts[beta] = np.abs(f_t).max()
    if beta == 0.0:
        lam_t = np.array([pendulum_lambda(sol.y[:, k]) for k in range(len(tt))])
        theta_t = np.arctan2(sol.y[0], -sol.y[1])
        thd_t = (sol.y[0] * sol.y[3] - sol.y[1] * sol.y[2]) / ELL_P**2
        T_t = M_P * (G * np.cos(theta_t) + ELL_P * thd_t**2)       # tension from B2(b)
        assert np.allclose(lam_t, -T_t / (2 * ELL_P), rtol=1e-3)     # polar conversion inherits the drift
        sol_B2 = sol
T_B2c = 1.0 * (G * np.cos(0.0) + 1.0 * 2.0**2)
assert np.isclose(T_B2c, 13.81)
assert drifts[5.0] < drifts[0.0] / 50
print(f"B2(c): T(theta=0, theta_dot=2) = {T_B2c:.2f} N.  B2(d): lambda = -T/(2 ell) along the DAE solution; "
      f"max|f| over 10 s (RK45, rtol 1e-6): {drifts[0.0]:.2e} m^2 (plain), {drifts[5.0]:.2e} m^2 (Baumgarte beta=5), "
      f"ratio {drifts[0.0]/drifts[5.0]:.0f}")
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(tt, lam_t, label=r"$\lambda$ from KKT"); ax[0].plot(tt, -T_t / (2 * ELL_P), "--", label=r"$-T/(2\ell)$, $T = m(g\cos\theta + \ell\dot\theta^2)$")
ax[0].set_xlabel("t [s]"); ax[0].legend(fontsize=8); ax[0].set_title("B2: pendulum multiplier")
ax[1].semilogy(tt, np.abs(sol_B2.y[0]**2 + sol_B2.y[1]**2 - 1) + 1e-18, label=r"$|f|$, plain (3.21)")
sol_b = solve_ivp(pendulum_dae_rhs_factory(5.0), (0, 10), y0, t_eval=tt, rtol=1e-6, atol=1e-8)
ax[1].semilogy(tt, np.abs(sol_b.y[0]**2 + sol_b.y[1]**2 - 1) + 1e-18, label=r"$|f|$, Baumgarte $\beta = 5$")
ax[1].set_xlabel("t [s]"); ax[1].legend(fontsize=8); ax[1].set_title("constraint drift")
fig.tight_layout(); fig.savefig(FIG / "pendulum_dae_B2.png", dpi=120); plt.close(fig)

# --- B4: bead on a rotating parabolic wire
k = sp.symbols("k", positive=True)
r_ = sp.Function("r")(t)
rpos = sp.Matrix([r_ * sp.cos(w * t), r_ * sp.sin(w * t), half * k * r_**2])
vpos = rpos.diff(t)
L_par = half * m * vpos.dot(vpos) - m * g * rpos[2]
L_par_hand = half * m * (1 + k**2 * r_**2) * r_.diff(t)**2 + half * m * (w**2 - g * k) * r_**2
assert sp.simplify(sp.expand(L_par - L_par_hand)) == 0                        # B4(a)
par = lagrangian_to_eom(L_par_hand, [r_], t, build_rhs=False)
rdd_hand = ((w**2 - g * k) * r_ - k**2 * r_ * r_.diff(t)**2) / (1 + k**2 * r_**2)
assert sp.simplify(par.acc[0] - rdd_hand) == 0                                 # B4(b)
KV = 2.0; w_crit = np.sqrt(G * KV)
outcomes = {}
for label, wv, r_start in (("below", 0.5 * w_crit, 0.3), ("at", w_crit, 0.3), ("above", 1.5 * w_crit, 0.3),
                           ("below, small amplitude", 0.5 * w_crit, 0.02)):
    rhs = lagrangian_to_eom(L_par_hand, [r_], t, params={m: 1, g: G, k: KV, w: wv}).rhs
    s = solve_ivp(rhs, (0, 6), [r_start, 0.0], t_eval=np.linspace(0, 6, 1201), rtol=1e-10, atol=1e-12)
    outcomes[label] = s
assert np.abs(outcomes["at"].y[0] - 0.3).max() < 1e-9                         # every r is an equilibrium
assert np.abs(outcomes["below"].y[0]).max() < 0.31                             # bounded oscillation about 0
P_par = period_from_crossings(outcomes["below, small amplitude"].t, outcomes["below, small amplitude"].y[0])
Om_par = np.sqrt(G * KV - (0.5 * w_crit)**2)                                   # linearized frequency about r = 0
assert abs(P_par / (2 * np.pi / Om_par) - 1) < 2e-3, P_par
P_par_big = period_from_crossings(outcomes["below"].t, outcomes["below"].y[0])  # the (1 + k^2 r^2) mass term slows it
assert P_par_big > 1.02 * P_par
assert outcomes["above"].y[0][-1] > 10                                         # the bead leaves
outcomes.pop("below, small amplitude")
print(f"B4: omega_crit = sqrt(gk) = {w_crit:.4f} rad/s (k = 2/m). Below (omega = omega_crit/2): small-amplitude period {P_par:.4f} s "
      f"vs 2pi/sqrt(gk - omega^2) = {2*np.pi/Om_par:.4f} s, amplitude-0.3 period {P_par_big:.4f} s; "
      f"at threshold r stays at 0.3 exactly; above: r(6 s) = {outcomes['above'].y[0][-1]:.1f} m")
fig, ax = plt.subplots(figsize=(6, 3.4))
for label, s in outcomes.items():
    ax.plot(s.t, s.y[0], label=f"omega {label} sqrt(gk)")
ax.set_ylim(-0.5, 3); ax.set_xlabel("t [s]"); ax.set_ylabel("r [m]"); ax.legend(); ax.set_title("B4: parabolic wire: no pitchfork")
fig.tight_layout(); fig.savefig(FIG / "parabolic_wire_B4.png", dpi=120); plt.close(fig)

# --- B7(c): Rayleigh damping on the double pendulum
F_dp = T_dp.subs({m1: c, m2: c})               # (m1 + m2, m2) -> (2c, c): F = ½ c (|r1d|² + |r2d|²)
Q_dp = [-sp.diff(F_dp, th1.diff(t)), -sp.diff(F_dp, th2.diff(t))]
dp_damp_params = {**dp_params, c: 0.2}
dpd = lagrangian_to_eom(L_dp, [th1, th2], t, Q=Q_dp, params=dp_damp_params)
tt = np.linspace(0, 60, 12001)
sol = solve_ivp(dpd.rhs, (0, 60), [2.0, 2.0, 0.0, 0.0], t_eval=tt, rtol=1e-10, atol=1e-12)
E_d = h_dp(sol.t, sol.y)
F_num = sp.lambdify([sp.Symbol("theta1"), sp.Symbol("theta2"), sp.Symbol("theta1_d"), sp.Symbol("theta2_d")],
                    F_dp.subs(c, 0.2).subs([(th1.diff(t), sp.Symbol("theta1_d")), (th2.diff(t), sp.Symbol("theta2_d"))])
                    .subs([(th1, sp.Symbol("theta1")), (th2, sp.Symbol("theta2"))]).subs(dp_params), "numpy")
F_t = F_num(*sol.y)
dEdt = np.gradient(E_d, tt)
rel_err_B7 = np.abs(dEdt + 2 * F_t)[5:-5].max() / np.abs(2 * F_t).max()
assert rel_err_B7 < 5e-3, rel_err_B7
V_min_dp = float(V_dp.subs([(th1, 0), (th2, 0)]).subs(dp_params))      # -3 g = -29.43 J: both bobs hanging at rest
assert np.isclose(V_min_dp, -29.43) and E_d[-1] - V_min_dp < 1e-3 * (E_d[0] - V_min_dp)
print(f"B7(c): double pendulum with c = 0.2: E(0) = {E_d[0]:.3f} J -> E(60 s) = {E_d[-1]:.3f} J, i.e. the energy above the "
      f"rest value V_min = {V_min_dp:.2f} J goes from {E_d[0]-V_min_dp:.2f} J to {E_d[-1]-V_min_dp:.2e} J; "
      f"dE/dt = -2F holds to {rel_err_B7:.1e} (relative to max |2F|)")
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(tt, E_d); ax[0].set_xlabel("t [s]"); ax[0].set_ylabel("E [J]"); ax[0].set_title("B7: damped double pendulum energy")
ax[1].plot(tt[:2000], dEdt[:2000], label="dE/dt"); ax[1].plot(tt[:2000], -2 * F_t[:2000], "--", label=r"$-2\mathcal{F}$")
ax[1].set_xlabel("t [s]"); ax[1].legend(); ax[1].set_title("first 10 s")
fig.tight_layout(); fig.savefig(FIG / "double_pendulum_damped_B7.png", dpi=120); plt.close(fig)

# --- C2(d): Kapitza's pendulum
A_K, ELL_K = 0.05, 1.0


def kapitza_rhs(gamma_v):
    def rhs(tt, y):
        return [y[1], -(G - A_K * gamma_v**2 * np.cos(gamma_v * tt)) * np.sin(y[0]) / ELL_K]
    return rhs


kap = {}
for gamma_v in (150.0, 60.0):
    s = solve_ivp(kapitza_rhs(gamma_v), (0, 10), [np.pi + 0.1, 0.0], t_eval=np.linspace(0, 10, 4001),
                  rtol=1e-9, atol=1e-11, max_step=0.005)
    kap[gamma_v] = s
dev150 = np.abs(kap[150.0].y[0] - np.pi).max(); dev60 = np.abs(kap[60.0].y[0] - np.pi).max()
assert A_K**2 * 150**2 > 2 * G * ELL_K > A_K**2 * 60**2       # 56.25 > 19.62 > 9
assert dev150 < 0.3 and dev60 > 2.0
Om_eff = np.sqrt(A_K**2 * 150**2 / (2 * ELL_K**2) - G / ELL_K)
print(f"C2(d): a^2 gamma^2 = {A_K**2*150**2:.2f} (gamma=150) and {A_K**2*60**2:.2f} (gamma=60) vs 2 g ell = {2*G*ELL_K:.2f}. "
      f"gamma=150: max|theta - pi| = {dev150:.3f} rad (stays inverted; slow period {2*np.pi/Om_eff:.2f} s); "
      f"gamma=60: max|theta - pi| = {dev60:.2f} rad (falls)")
fig, ax = plt.subplots(figsize=(6.5, 3.4))
for gamma_v, s in kap.items():
    ax.plot(s.t, s.y[0] - np.pi, lw=.8, label=rf"$\gamma$ = {gamma_v:.0f} rad/s")
ax.set_xlabel("t [s]"); ax.set_ylabel(r"$\theta - \pi$ [rad]"); ax.legend(); ax.set_title("C2: Kapitza's inverted pendulum")
fig.tight_layout(); fig.savefig(FIG / "kapitza_C2.png", dpi=120); plt.close(fig)


# --- C3(d): slowly varying length, adiabatic invariant
def adiabatic_rhs(tt, y):
    ell_t, elld = 1 + 0.02 * tt, 0.02
    return [y[1], (-2 * elld * y[1] - G * np.sin(y[0])) / ell_t]


tt = np.linspace(0, 50, 50001)
sol = solve_ivp(adiabatic_rhs, (0, 50), [0.1, 0.0], t_eval=tt, rtol=1e-10, atol=1e-12)
ell_t = 1 + 0.02 * tt
E_ad = 0.5 * ell_t**2 * sol.y[1]**2 + 0.5 * G * ell_t * sol.y[0]**2
om_ad = np.sqrt(G / ell_t)
turning = np.where(np.diff(np.sign(sol.y[1])) != 0)[0]           # theta_dot changes sign: turning points
I_ad = (E_ad / om_ad)[turning]
amp_ad = (np.abs(sol.y[0]) * ell_t**0.75)[turning]
spread_I = np.ptp(I_ad) / I_ad.mean(); spread_amp = np.ptp(amp_ad) / amp_ad.mean()
assert spread_I < 0.012 and spread_amp < 0.012
assert E_ad[turning][-1] / E_ad[turning][0] < 0.75                 # E itself falls by ~ ell^{-1/2} = 1/sqrt(2)
print(f"C3(d): E/omega flat to {100*spread_I:.2f}% and theta_max ell^(3/4) flat to {100*spread_amp:.2f}% over {len(turning)} turning points, "
      f"while E falls to {E_ad[turning][-1]/E_ad[turning][0]:.3f} of its start (prediction ell^(-1/2) = {2**-0.5:.3f})")
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(tt[turning], I_ad / I_ad[0], ".-", ms=3, label=r"$E/\omega$ (normalized)")
ax[0].plot(tt[turning], amp_ad / amp_ad[0], ".-", ms=3, label=r"$\theta_{\max}\ell^{3/4}$ (normalized)")
ax[0].plot(tt[turning], E_ad[turning] / E_ad[turning][0], "--", label=r"$E$ (normalized)")
ax[0].set_xlabel("t [s]"); ax[0].legend(fontsize=8); ax[0].set_title("C3: adiabatic invariant")
ax[1].plot(tt, sol.y[0], lw=.5); ax[1].set_xlabel("t [s]"); ax[1].set_ylabel(r"$\theta$"); ax[1].set_title(r"$\ell = 1 + 0.02t$")
fig.tight_layout(); fig.savefig(FIG / "adiabatic_C3.png", dpi=120); plt.close(fig)

# --- A5 and B5 numbers quoted in solutions.md
dEdt_A5 = -0.5 * 1.0**2
assert dEdt_A5 == -2 * (0.5 * 0.5 * 1.0**2)
r_B5 = np.sqrt(1 / 2)
mu_B5 = np.array([1 / (1 + r_B5), 1 / (1 - r_B5)])
assert np.allclose(mu_B5, [2 - np.sqrt(2), 2 + np.sqrt(2)])
print(f"A5: dE/dt(0) = -c v^2 = {dEdt_A5:.2f} W (E(0) = {0.5*1*1**2 + 0.5*4*1**2:.2f} J).  "
      f"B5: mu_-+ for m1 = m2: {mu_B5[0]:.4f}, {mu_B5[1]:.4f} = 2 -+ sqrt2")

# %%
# --- Symbolic checks of the pen-and-paper solutions (each line mirrors a result in solutions.md)
Z = lambda e: sp.simplify(e) == 0                                   # noqa: E731
Zm = lambda A: sp.simplify(A) == sp.zeros(*A.shape)                 # noqa: E731
# A3: constant horizontal force: Q = F ell cos(theta); same EOM as V = -F x; equilibrium tan(theta) = F/(m g)
Fh = sp.symbols("F", positive=True)
Q_A3 = (sp.Matrix([Fh, 0]).T * sp.Matrix([ell * sp.sin(th), -ell * sp.cos(th)]).diff(th))[0]
a3 = lagrangian_to_eom(L_pend, [th], t, Q=[Q_A3], build_rhs=False)
assert Z(Q_A3 - Fh * ell * sp.cos(th)) and Z(a3.acc[0] - lagrangian_to_eom(L_pend + Fh * ell * sp.sin(th), [th], t, build_rhs=False).acc[0])
assert Z(a3.acc[0].subs(th, sp.atan(Fh / (m * g))))
# A4: J^T diag(M, m, m) J is the wedge mass matrix, det = m (M + m sin^2 alpha)
Xs, ss_, hh = sp.symbols("X s h")
J_A4 = sp.Matrix([Xs, Xs + ss_ * sp.cos(alpha), hh - ss_ * sp.sin(alpha)]).jacobian([Xs, ss_])
M_A4 = J_A4.T * sp.diag(Mw, m, m) * J_A4
assert Zm(M_A4 - M_wedge_hand) and Z(M_A4.det() - m * (Mw + m * sp.sin(alpha)**2)) and Z(M_A4.trace() - (Mw + 2 * m))
# B1(c): cancellation of dots (3.7) and interchange (3.8) for the rotating hoop, with the explicit vector
r_h = sp.Matrix([R * sp.sin(th) * sp.cos(w * t), R * sp.sin(th) * sp.sin(w * t), -R * sp.cos(th)])
assert Zm(r_h.diff(t).diff(th.diff(t)) - r_h.diff(th)) and Zm(r_h.diff(t).diff(th) - r_h.diff(th).diff(t))
assert Zm(r_h.diff(th).diff(t) - (R * th.diff(t) * sp.Matrix([-sp.sin(th) * sp.cos(w * t), -sp.sin(th) * sp.sin(w * t), sp.cos(th)])
                                   + R * w * sp.cos(th) * sp.Matrix([-sp.sin(w * t), sp.cos(w * t), 0])))
# B3: massive pulley. Unknowns xdd1 = a1 (xdd2 = -a1, phidd = a1/Rp) and the two multipliers
Ip, Rp, la1, la2, a1s = sp.symbols("I R_p lambda1 lambda2 a1")
s3 = sp.solve([m1 * a1s - m1 * g - la1 - la2, -m2 * a1s - m2 * g - la1, Ip * a1s / Rp + Rp * la2], [a1s, la1, la2], dict=True)[0]
D3 = m1 + m2 + Ip / Rp**2
assert Z(s3[a1s] - (m1 - m2) * g / D3) and Z(s3[la1] + m2 * g * (2 * m1 + Ip / Rp**2) / D3) and Z(s3[la2] + Ip * (m1 - m2) * g / (Rp**2 * D3))
assert Z(-(s3[la1] + s3[la2]) - m1 * (g - s3[a1s])) and Z(-s3[la1] - m2 * (g + s3[a1s]))          # T1, T2 from Newton
assert Z(s3[la1].subs(Ip, 0) + 2 * m1 * m2 * g / (m1 + m2))                                         # I -> 0 gives (4.4)
# B5: general masses, ell1 = ell2: det(K - mu (g/ell) M0) = 0 at mu = 1/(1 -+ sqrt(m2/(m1+m2)))
mu = sp.symbols("mu")
M0g = ell**2 * sp.Matrix([[m1 + m2, m2], [m2, m2]]); Kg = g * ell * sp.diag(m1 + m2, m2)
char = (Kg - mu * g / ell * M0g).det()
assert all(Z(char.subs(mu, 1 / (1 - sgn * sp.sqrt(m2 / (m1 + m2))))) for sgn in (1, -1))
# B6: V_eff'' at the conical equilibrium = m ell^2 phid0^2 (1 + 3 cos^2 theta0)
t0 = sp.symbols("theta0", positive=True)
Veff6 = pphi**2 / (2 * m * ell**2 * sp.sin(t0)**2) - m * g * ell * sp.cos(t0)
pphi0 = m * ell**2 * sp.sin(t0)**2 * sp.sqrt(g / (ell * sp.cos(t0)))
assert Z(Veff6.diff(t0).subs(pphi, pphi0)) and Z(Veff6.diff(t0, 2).subs(pphi, pphi0) - m * g * ell * (1 + 3 * sp.cos(t0)**2) / sp.cos(t0))
# B8(a): inertially fixed point seen from the rotating frame: Coriolis = -2 w^2 r', centrifugal = +w^2 r' (per unit mass)
r0s = sp.symbols("r0", positive=True)
rp8 = sp.Matrix([r0s * sp.cos(w * t), -r0s * sp.sin(w * t), 0]); om8 = sp.Matrix([0, 0, w])
cor8 = -2 * om8.cross(rp8.diff(t)); cen8 = -om8.cross(om8.cross(rp8))
assert Zm(cor8 + 2 * w**2 * rp8) and Zm(cen8 - w**2 * rp8) and Zm(rp8.diff(t, 2) - cor8 - cen8)
# C1(c): rolling disk, Lagrange–d'Alembert with a1 = (1,0,0,-R cos th), a2 = (0,1,0,-R sin th) and the constraints differentiated
Ith, Iph, thv, thdv, phdv = sp.symbols("I_theta I_phi theta theta_d phi_d")
xdd, ydd, thdd, phdd, l1s, l2s = sp.symbols("xdd ydd thdd phdd l1 l2")
eqC1 = [m * xdd - l1s, m * ydd - l2s, Ith * thdd, Iph * phdd + R * sp.cos(thv) * l1s + R * sp.sin(thv) * l2s,
        xdd - R * phdd * sp.cos(thv) + R * phdv * thdv * sp.sin(thv), ydd - R * phdd * sp.sin(thv) - R * phdv * thdv * sp.cos(thv)]
sC1 = sp.solve(eqC1, [xdd, ydd, thdd, phdd, l1s, l2s], dict=True)[0]
assert sC1[thdd] == 0 and Z(sC1[phdd]) and Z(sC1[l1s] + m * R * phdv * thdv * sp.sin(thv)) and Z(sC1[l2s] - m * R * phdv * thdv * sp.cos(thv))
# C1(d): skate on an incline. Substituting ydot = xdot tan(theta) into L and varying gives I thdd = m xd^2 sec^2 tan  (wrong)
Isk, beta = sp.symbols("I beta", positive=True)
xk, thk = sp.Function("x")(t), sp.Function("theta")(t)
L_vak = half * m * (xk.diff(t)**2 + (xk.diff(t) * sp.tan(thk))**2) + half * Isk * thk.diff(t)**2 + m * g * xk * sp.sin(beta)
vak = lagrangian_to_eom(L_vak, [xk, thk], t, build_rhs=False)
assert Z(vak.acc[1] - m * xk.diff(t)**2 * sp.tan(thk) / (Isk * sp.cos(thk)**2))
# C2(a): vertically driven pivot; (c) Kapitza's effective potential and the averaged fast force
y_v = a * sp.cos(gam * t) - ell * sp.cos(th)
L_vert = half * m * ((ell * sp.sin(th)).diff(t)**2 + y_v.diff(t)**2) - m * g * y_v
assert Z(lagrangian_to_eom(L_vert, [th], t, build_rhs=False).acc[0] + (g - a * gam**2 * sp.cos(gam * t)) * sp.sin(th) / ell)
Th, u = sp.symbols("Theta u")
xi = -(a / ell) * sp.cos(u) * sp.sin(Th)                                   # fast response, u = gamma t
avg_force = a * gam**2 / ell * sp.integrate(sp.cos(u) * xi * sp.cos(Th), (u, 0, 2 * sp.pi)) / (2 * sp.pi)
VK = -m * g * ell * sp.cos(Th) + m * a**2 * gam**2 / 4 * sp.sin(Th)**2
assert Z(-(g / ell) * sp.sin(Th) + avg_force + VK.diff(Th) / (m * ell**2))
assert Z(VK.diff(Th, 2).subs(Th, sp.pi) - (m * a**2 * gam**2 / 2 - m * g * ell))
# C3(a): ell(t) pendulum keeps the full T (with ½ m elld^2): same EOM; (c) exact dE/dt for the linearized motion
lf = sp.Function("ell")(t)
rv3 = lf * sp.Matrix([sp.sin(th), -sp.cos(th)])
L_var = half * m * rv3.diff(t).dot(rv3.diff(t)) - m * g * rv3[1]
c3 = lagrangian_to_eom(L_var, [th], t, build_rhs=False)
assert Z(c3.acc[0] + (2 * lf.diff(t) * th.diff(t) + g * sp.sin(th)) / lf)
E3 = half * m * lf**2 * th.diff(t)**2 + half * m * g * lf * th**2
dE3 = E3.diff(t).subs(th.diff(t, 2), -(2 * lf.diff(t) * th.diff(t) + g * th) / lf)
assert Z(dE3 - lf.diff(t) / lf * (half * m * g * lf * th**2 - m * lf**2 * th.diff(t)**2))
print("Symbolic checks of solutions A3, A4, B1(c), B3, B5, B6, B8(a), C1(c,d), C2(a,c), C3(a,c) passed")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.perf_counter() - T_START:.1f} s. Figures in {FIG}")
