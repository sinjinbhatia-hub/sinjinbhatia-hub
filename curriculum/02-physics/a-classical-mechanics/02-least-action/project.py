# %% [markdown]
# # CM-2 Project: Least action, discretized, minimized, and checked
#
# Hamilton's principle says the true path makes S = ∫(T − V) dt stationary.
# This script makes that literal: it chops a path into N points, writes the
# action as an ordinary function of N variables (notes eq. 2.1), and hands it
# to L-BFGS-B. The minimizer lands on the ODE solution. Then it checks the
# claims of the notes one by one: S(ε) is flat at ε = 0 (§3.1), the oscillator's
# path stops being a minimum past half a period (Theorem 3.6), Euler–Lagrange
# in polar coordinates and for a charged particle (SymPy, §3.5, §3.9),
# cyclotron motion, and coordinate covariance (Theorem 3.5) on the pendulum.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import warnings
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import minimize, brentq
from scipy.integrate import solve_ivp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
G = 9.81
XG, WG = np.polynomial.legendre.leggauss(200)  # Gauss–Legendre nodes on [-1, 1]


def gl_integrate(f, a, b):
    """∫_a^b f(t) dt by 200-point Gauss–Legendre (exact to ~1e-14 for smooth f)."""
    t = 0.5 * (b - a) * XG + 0.5 * (b + a)
    return 0.5 * (b - a) * np.sum(WG * f(t))


def action(q, qd, a, b, M, V):
    """S = ∫ [½ M q̇² − V(q)] dt along an explicit path q(t), eq. (3.1)."""
    return gl_integrate(lambda t: 0.5 * M * qd(t) ** 2 - V(q(t)), a, b)


# %% [markdown]
# ## Part 1: The discretized action, minimized (§2, Example 4.1)
#
# Nodes q_0..q_{N+1} with the ends fixed; eq. (2.1):
#   S_h = Σ_k h [ ½ M ((q_{k+1} − q_k)/h)² − ½ (V(q_k) + V(q_{k+1})) ].
# Its gradient, eq. (2.2), is −h [M (q_{k+1} − 2q_k + q_{k−1})/h² + V'(q_k)]:
# minus h times the Newton residual at node k. **Build it yourself first**:
# write `discrete_action` returning (S_h, ∇S_h) before reading it.


# %%
def discrete_action(q_int, qa, qb, tau, M, V, dV):
    """S_h of eq. (2.1) and its gradient (2.2) with respect to the interior nodes."""
    q = np.concatenate([[qa], q_int, [qb]])
    h = tau / (len(q) - 1)
    dq = np.diff(q)
    S = np.sum(0.5 * M * dq**2 / h) - h * np.sum(0.5 * (V(q[:-1]) + V(q[1:])))
    grad = M * (q[1:-1] - q[:-2]) / h - M * (q[2:] - q[1:-1]) / h - h * dV(q[1:-1])
    return S, grad


def minimize_action(qa, qb, tau, M, V, dV, N=199):
    """Minimize S_h over the N interior nodes with L-BFGS-B, starting from a straight line."""
    t = np.linspace(0, tau, N + 2)
    res = minimize(discrete_action, np.interp(t[1:-1], [0, tau], [qa, qb]),
                   args=(qa, qb, tau, M, V, dV), jac=True, method="L-BFGS-B",
                   options=dict(maxiter=20000, gtol=1e-14, ftol=1e-16, maxcor=30))
    return t, np.concatenate([[qa], res.x, [qb]]), res.fun


def shoot(rhs, qa, qb, tau, lo=-20, hi=20):
    """Solve the two-point problem q(0)=qa, q(tau)=qb by shooting on q̇(0); returns dense solution."""
    end = lambda v0: solve_ivp(rhs, (0, tau), [qa, v0], rtol=1e-12, atol=1e-12).y[0, -1] - qb
    v0 = brentq(end, lo, hi, xtol=1e-14)
    return solve_ivp(rhs, (0, tau), [qa, v0], rtol=1e-12, atol=1e-12, dense_output=True).sol


systems = {  # name: (M, V, dV, rhs, qa, qb, tau)
    "projectile": (1.0, lambda y: G * y, lambda y: G * np.ones_like(y),
                   lambda t, s: [s[1], -G], 0.0, 0.0, 1.0),
    "oscillator": (1.0, lambda q: 0.5 * q**2, lambda q: q,
                   lambda t, s: [s[1], -s[0]], 0.0, 1.0, 2.0),
    "pendulum": (1.0, lambda th: -G * np.cos(th), lambda th: G * np.sin(th),
                 lambda t, s: [s[1], -G * np.sin(s[0])], 0.3, -0.2, 0.8),
}
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
for ax, (name, (M, V, dV, rhs, qa, qb, tau)) in zip(axes, systems.items()):
    t, q_min, S_min = minimize_action(qa, qb, tau, M, V, dV)
    sol = shoot(rhs, qa, qb, tau)
    err = np.abs(q_min - sol(t)[0]).max()
    assert err < (1e-6 if name == "projectile" else 1e-4), (name, err)  # Verlet is exact for parabolas
    print(f"{name:10s}: max |L-BFGS-B minimizer − ODE solution| = {err:.1e}, S_h = {S_min:.5f}")
    tt = np.linspace(0, tau, 300)
    ax.plot(tt, sol(tt)[0], "k-", lw=3, alpha=.3, label="ODE (shooting)")
    ax.plot(t[::10], q_min[::10], "o", ms=4, label="argmin of discrete action")
    ax.set_title(name); ax.set_xlabel("t"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "discrete_action_vs_ode.png", dpi=120); plt.close(fig)

# Example 4.1 / §2: thrown ball, family y_h = 4h t(τ−t)/τ², S(h) = 8mh²/(3τ) − (2/3) m g h τ, eq. (2.3).
m, tau = 1.0, 1.0
S_of_h = lambda h: action(lambda t: 4 * h * t * (tau - t) / tau**2, lambda t: 4 * h * (tau - 2 * t) / tau**2,
                          0, tau, m, lambda y: m * G * y)
for h in (0.3, 1.0, 2.0):
    assert np.isclose(S_of_h(h), 8 * m * h**2 / (3 * tau) - 2 / 3 * m * G * h * tau)
h_star = G * tau**2 / 8
S_cl_ball = -m * G**2 * tau**3 / 24
assert np.isclose(h_star, 1.22625) and np.isclose(S_cl_ball, -4.00984, atol=5e-6)   # numbers quoted in Example 4.1
assert np.isclose(S_of_h(h_star), S_cl_ball) and S_of_h(h_star - .1) > S_cl_ball < S_of_h(h_star + .1)
S_pert = action(lambda t: 0.5 * G * t * (tau - t) + np.sin(np.pi * t),   # η = sin(πt), η(0) = η(1) = 0
                lambda t: 0.5 * G * (tau - 2 * t) + np.pi * np.cos(np.pi * t), 0, tau, m, lambda y: m * G * y)
ke_eta = gl_integrate(lambda t: 0.5 * m * (np.pi * np.cos(np.pi * t))**2, 0, tau)
assert np.isclose(ke_eta, np.pi**2 / 4)
assert np.isclose(S_pert, S_cl_ball + ke_eta)  # S[q_cl + η] = S[q_cl] + ½m∫η̇² exactly (Example 4.1)
print(f"Example 4.1: h* = gτ²/8 = {h_star:.5f} m, S_cl = −mg²τ³/24 = {S_cl_ball:.5f} J·s; "
      f"S[q_cl+η] − S_cl = ½m∫η̇² = {ke_eta:.4f}")
print(f"§3.10: |S_cl|/ħ for this ball = {abs(S_cl_ball) / 1.054571817e-34:.1e} radians of phase")

# %% [markdown]
# ## Part 2: S(ε) is flat at ε = 0 on the true path, and only there (§3.1, Example 4.4)

# %%
def S_eps(q, qd, eta, etad, a, b, M, V, eps):
    """S[q + εη] for a perturbation η with η(a) = η(b) = 0 (§3.1)."""
    return action(lambda t: q(t) + eps * eta(t), lambda t: qd(t) + eps * etad(t), a, b, M, V)


def dS_deps(*args, e=1e-3):
    """Central difference for dS[q + εη]/dε at ε = 0: the first variation (3.2)."""
    return (S_eps(*args, e) - S_eps(*args, -e)) / (2 * e)


etas = [(lambda tau: (lambda t: np.sin(np.pi * t / tau), lambda t: np.pi / tau * np.cos(np.pi * t / tau))),
        (lambda tau: (lambda t: t * (tau - t)**2, lambda t: (tau - t)**2 - 2 * t * (tau - t)))]
fig, ax = plt.subplots(figsize=(6, 4))
for name in ("oscillator", "pendulum"):
    M, V, dV, rhs, qa, qb, tau = systems[name]
    sol = shoot(rhs, qa, qb, tau)
    q_cl, qd_cl = (lambda t: sol(t)[0]), (lambda t: sol(t)[1])
    q_ln, qd_ln = (lambda t: qa + (qb - qa) * t / tau), (lambda t: (qb - qa) / tau + 0 * t)
    for k, mk in enumerate(etas):
        eta, etad = mk(tau)
        slope_true = dS_deps(q_cl, qd_cl, eta, etad, 0, tau, M, V)
        slope_line = dS_deps(q_ln, qd_ln, eta, etad, 0, tau, M, V)
        assert abs(slope_true) < 1e-6 and abs(slope_line) > 1e-2, (name, slope_true, slope_line)
        print(f"{name:10s} η{k+1}: dS/dε at ε=0 on true path = {slope_true:+.1e}; on straight line = {slope_line:+.3f}")
        eps = np.linspace(-.3, .3, 61)
        S0 = S_eps(q_cl, qd_cl, eta, etad, 0, tau, M, V, 0)
        ax.plot(eps, [S_eps(q_cl, qd_cl, eta, etad, 0, tau, M, V, e) - S0 for e in eps], label=f"{name}, η{k+1}")
ax.axvline(0, color="k", lw=.5); ax.set_xlabel("ε"); ax.set_ylabel("S[q_cl + εη] − S[q_cl]")
ax.legend(fontsize=8); ax.set_title("Flat at ε = 0: the first variation vanishes")
fig.tight_layout(); fig.savefig(FIG / "S_of_epsilon.png", dpi=120); plt.close(fig)

# Example 4.4: free particle, m = 1, x(0) = 0, x(1) = 1. Straight line vs the bent path via x = 2 at t = ½.
seg_action = lambda xs, ts, m=1.0: sum(0.5 * m * ((x1 - x0) / (t1 - t0))**2 * (t1 - t0)
                                       for x0, x1, t0, t1 in zip(xs[:-1], xs[1:], ts[:-1], ts[1:]))
S_straight, S_bent = seg_action([0, 1], [0, 1]), seg_action([0, 2, 1], [0, .5, 1])
eta_pen = seg_action([0, 1.5, 0], [0, .5, 1])  # ½m∫η̇² for η = bent − straight
assert np.isclose(S_straight, 0.5) and np.isclose(S_bent, 5.0) and np.isclose(S_bent, S_straight + eta_pen)
print(f"Example 4.4: S_straight = {S_straight}, S_bent = {S_bent} = {S_straight} + ½m∫η̇² ({eta_pen})")

# %% [markdown]
# ## Part 3: Is it a minimum? The oscillator and the T/2 threshold (Theorem 3.6, Example 4.5)
#
# m = k = 1, so ω = 1, T = 2π, T/2 = π. For q_cl from q(0) = 0 to q(τ) = 1 and η = sin(πt/τ),
# S[q_cl + εη] − S[q_cl] = c(τ) ε² exactly, with c(τ) = (τ/4) ((π/τ)² − 1), eq. (4.2).

# %%
def c_coeff(tau):
    return 0.5 * (tau / 2) * ((np.pi / tau)**2 - 1)


osc_V = lambda q: 0.5 * q**2
fig, ax = plt.subplots(figsize=(6, 4))
for tau in (2.0, 4.0):
    q, qd = (lambda t, T=tau: np.sin(t) / np.sin(T)), (lambda t, T=tau: np.cos(t) / np.sin(T))
    eta, etad = (lambda t, T=tau: np.sin(np.pi * t / T)), (lambda t, T=tau: np.pi / T * np.cos(np.pi * t / T))
    S0 = action(q, qd, 0, tau, 1.0, osc_V)
    assert np.isclose(S0, 0.5 / np.tan(tau))                      # S_cl = ½ cot τ for these endpoints
    for e in (0.1, 0.7, 2.0):                                      # exact parabola in ε, any size
        assert np.isclose(S_eps(q, qd, eta, etad, 0, tau, 1.0, osc_V, e) - S0, c_coeff(tau) * e**2)
    eps = np.linspace(-1, 1, 81)
    ax.plot(eps, [S_eps(q, qd, eta, etad, 0, tau, 1.0, osc_V, e) for e in eps],
            label=f"τ = {tau:g} ({'<' if tau < np.pi else '>'} T/2 = π)")
    print(f"Example 4.5: τ = {tau}: S_cl = ½cot τ = {S0:+.5f}, c(τ) = {c_coeff(tau):+.5f}")
assert np.isclose(c_coeff(2.0), 0.73370, atol=5e-6) and np.isclose(c_coeff(4.0), -0.38315, atol=5e-6)
assert np.isclose(0.5 / np.tan(2.0), -0.22883, atol=5e-6) and np.isclose(0.5 / np.tan(4.0), 0.43185, atol=5e-6)
print(f"Example 4.5: π²/4 = {np.pi**2 / 4:.5f}, π²/16 = {np.pi**2 / 16:.5f}, threshold τ = π = {np.pi:.4f}")
ax.set_xlabel("ε"); ax.set_ylabel("S[q_cl + ε sin(πt/τ)]"); ax.legend(); ax.set_title("Minimum, then saddle")
fig.tight_layout(); fig.savefig(FIG / "oscillator_minimum_vs_saddle.png", dpi=120); plt.close(fig)


def discrete_hessian(tau, N=400, M=1.0, k=1.0):
    """Hessian of S_h for V = ½kq²: (M/h) tridiag(−1, 2, −1) − h k I."""
    h = tau / (N + 1)
    return (M / h) * (2 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)) - h * k * np.eye(N)


n_neg = {}
for tau in (2.0, 4.0, 5.0, 7.0):
    w, U = np.linalg.eigh(discrete_hessian(tau))
    n_neg[tau] = int(np.sum(w < 0))
    if tau == 4.0:  # the descending direction: lowest eigenvector should be sin(πt/τ)
        s = np.linspace(0, tau, 402)[1:-1]
        ref = np.sin(np.pi * s / tau); ref /= np.linalg.norm(ref)
        assert abs(U[:, 0] @ ref) > 0.9999
assert n_neg == {2.0: 0, 4.0: 1, 5.0: 1, 7.0: 2}, n_neg
print("Negative Hessian eigenvalues (descending directions) by τ:", n_neg, "(B4)")
# The discrete threshold is not exactly π: the lowest eigenvalue (4M/h) sin²(π/(2(N+1))) − hk vanishes at
# τ_h = 2(N+1) sin(π/(2(N+1))) = π − π³/(24(N+1)²) + ...  (Example 4.5), here π − 8.0e-6.
N_h = 400
tau_h = 2 * (N_h + 1) * np.sin(np.pi / (2 * (N_h + 1)))
assert np.isclose(np.pi - tau_h, 8.0e-6, rtol=1e-2)
assert np.linalg.eigvalsh(discrete_hessian(tau_h - 1e-6, N_h))[0] > 0 > np.linalg.eigvalsh(discrete_hessian(tau_h + 1e-6, N_h))[0]  # τ_h + 1e-6 < π
print(f"Discrete convexity threshold for N = {N_h}: τ_h = {tau_h:.7f} = π − {np.pi - tau_h:.1e}")
with warnings.catch_warnings(), np.errstate(all="ignore"):  # L-BFGS-B runs away past T/2
    warnings.simplefilter("ignore")
    t_ = np.linspace(0, 4.0, 101)[1:-1]
    S_run = minimize(discrete_action, t_ / 4.0, args=(0.0, 1.0, 4.0, 1.0, osc_V, lambda q: q),
                     jac=True, method="L-BFGS-B", options=dict(maxiter=30)).fun
assert S_run < -1e3, S_run
print(f"τ = 4 > T/2: L-BFGS-B on the discrete action runs off to S_h = {S_run:.2e} (unbounded below)")

# §3.8 and B5: S as a function of the endpoints. Oscillator S_cl(q1, q2, τ) = [(q1²+q2²)cos τ − 2q1q2]/(2 sin τ).
S_osc = lambda q1, q2, tau: ((q1**2 + q2**2) * np.cos(tau) - 2 * q1 * q2) / (2 * np.sin(tau))
d = 1e-6
dS_dq2 = (S_osc(0, 1 + d, 2) - S_osc(0, 1 - d, 2)) / (2 * d)
dS_dq1 = (S_osc(d, 1, 2) - S_osc(-d, 1, 2)) / (2 * d)
dS_dtau = (S_osc(0, 1, 2 + d) - S_osc(0, 1, 2 - d)) / (2 * d)
p0, p2 = 1 / np.sin(2), np.cos(2) / np.sin(2)                   # m q̇(0), m q̇(τ) on q = sin t / sin 2
assert np.isclose(dS_dq2, p2) and np.isclose(dS_dq1, -p0) and np.isclose(dS_dtau, -0.5 * p0**2)
print(f"B5: S_cl = {S_osc(0, 1, 2):.5f}; ∂S/∂q2 = {dS_dq2:.5f} = p(τ); ∂S/∂q1 = {dS_dq1:.5f} = −p(0); "
      f"∂S/∂τ = {dS_dtau:.5f} = −E")
S_free = lambda q2, t2, m=1.0: m * q2**2 / (2 * t2)            # free particle from (0, 0)
assert np.isclose(S_free(1, 1), .5)
assert np.isclose((S_free(1 + d, 1) - S_free(1 - d, 1)) / (2 * d), 1.0)        # = p2 = mv
assert np.isclose((S_free(1, 1 + d) - S_free(1, 1 - d)) / (2 * d), -0.5)       # = −E
print("§3.8 free particle: S = 0.5, ∂S/∂q2 = 1 = p2, ∂S/∂t2 = −0.5 = −E")

# %% [markdown]
# ## Part 4: Euler–Lagrange in SymPy (§3.2–3.5, §3.9)
#
# Slots are plain symbols (q, v, a); d/dt is the total derivative along a path:
# D f = ∂f/∂t + Σ v_i ∂f/∂q_i + Σ a_i ∂f/∂v_i. We return E_i = D(∂L/∂v_i) − ∂L/∂q_i, eq. (3.3).

# %%
t = sp.symbols("t", real=True)


def total_d(f, qs, vs, accs):
    """Total time derivative along a path: ∂f/∂t + Σ q̇_i ∂f/∂q_i + Σ q̈_i ∂f/∂q̇_i (the d/dt of (3.3))."""
    return sp.diff(f, t) + sum(sp.diff(f, q) * v for q, v in zip(qs, vs)) + \
        sum(sp.diff(f, v) * a for v, a in zip(vs, accs))


def euler_lagrange(L, qs, vs, accs):
    """E_i = d/dt(∂L/∂q̇_i) − ∂L/∂q_i, eq. (3.3); the equations of motion are E_i = 0."""
    return [sp.expand(total_d(sp.diff(L, v), qs, vs, accs) - sp.diff(L, q)) for q, v in zip(qs, vs)]


m_, k_, g_, l_, qe, B, al, be, c_ = sp.symbols("m k g ell q_e B alpha beta c", positive=True)
x, y, z, vx, vy, vz, ax_, ay_, az_ = sp.symbols("x y z v_x v_y v_z a_x a_y a_z", real=True)
r, th, vr, vth, ar, ath = sp.symbols("r theta v_r v_theta a_r a_theta", real=True)
Vf, Ff, Phi = sp.Function("V"), sp.Function("F"), sp.Function("Phi")
A = [sp.Function(f"A_{i}")(x, y, z, t) for i in "xyz"]
X3, V3, A3 = [x, y, z], [vx, vy, vz], [ax_, ay_, az_]

# Theorem 3.2: Cartesian, time-dependent V: E_x = m a_x + ∂V/∂x, and likewise for y, z.
E = euler_lagrange(m_ / 2 * (vx**2 + vy**2 + vz**2) - Vf(x, y, z, t), X3, V3, A3)
assert all(sp.simplify(Ei - (m_ * ai + sp.diff(Vf(x, y, z, t), xi))) == 0 for Ei, ai, xi in zip(E, A3, X3))
# Theorem 3.3: L = αT + βV gives α m ẍ − β V'(x).
E_ab = euler_lagrange(al * m_ / 2 * vx**2 + be * Vf(x), [x], [vx], [ax_])[0]
assert sp.simplify(E_ab - (al * m_ * ax_ - be * sp.diff(Vf(x), x))) == 0
# Theorem 3.4: the E–L expression of dF/dt vanishes identically, for any F(q1, q2, t); and E[cL] = c E[L].
q1, q2, v1, v2, a1, a2 = sp.symbols("q1 q2 v1 v2 a1 a2", real=True)
dFdt = total_d(Ff(q1, q2, t), [q1, q2], [v1, v2], [a1, a2])
assert all(sp.simplify(Ei) == 0 for Ei in euler_lagrange(dFdt, [q1, q2], [v1, v2], [a1, a2]))
L12 = m_ / 2 * (v1**2 + v2**2) - Vf(q1, q2, t)
assert all(sp.simplify(Ec - c_ * Ei) == 0 for Ec, Ei in zip(euler_lagrange(c_ * L12, [q1, q2], [v1, v2], [a1, a2]),
                                                           euler_lagrange(L12, [q1, q2], [v1, v2], [a1, a2])))
# Example 4.2: polar coordinates.
Lpol = m_ / 2 * (vr**2 + r**2 * vth**2) - Vf(r)
Er, Eth = euler_lagrange(Lpol, [r, th], [vr, vth], [ar, ath])
assert sp.simplify(Er - (m_ * ar - m_ * r * vth**2 + sp.diff(Vf(r), r))) == 0
assert sp.simplify(Eth - sp.expand(total_d(m_ * r**2 * vth, [r, th], [vr, vth], [ar, ath]))) == 0
assert sp.simplify(sp.solve(Er.subs({ar: 0, Vf(r): -k_ / r}).doit(), vth**2)[0] - k_ / (m_ * r**3)) == 0
# Theorem 3.5 (covariance) for Cartesian -> polar with a generic V(x, y): Ẽ_j = Σ_i (∂q_i/∂Q_j) E_i.
Ex, Ey = euler_lagrange(m_ / 2 * (vx**2 + vy**2) - Vf(x, y), [x, y], [vx, vy], [ax_, ay_])
X, Y = r * sp.cos(th), r * sp.sin(th)
D = lambda f: total_d(f, [r, th], [vr, vth], [ar, ath])
to_polar = {x: X, y: Y, vx: D(X), vy: D(Y), ax_: D(D(X)), ay_: D(D(Y))}
Er_c, Eth_c = euler_lagrange((m_ / 2 * (vx**2 + vy**2) - Vf(x, y)).subs(to_polar), [r, th], [vr, vth], [ar, ath])
Exp, Eyp = Ex.subs(to_polar), Ey.subs(to_polar)
assert sp.simplify((Er_c - (sp.diff(X, r) * Exp + sp.diff(Y, r) * Eyp)).doit()) == 0
assert sp.simplify((Eth_c - (sp.diff(X, th) * Exp + sp.diff(Y, th) * Eyp)).doit()) == 0
# Theorem 3.5 with a time-dependent map (B3): rotating axes x = X cos Ωt − Y sin Ωt, y = X sin Ωt + Y cos Ωt.
Om = sp.symbols("Omega", positive=True)
Xr, Yr, VXr, VYr, AXr, AYr = sp.symbols("X_r Y_r V_X V_Y A_X A_Y", real=True)
xr, yr = Xr * sp.cos(Om * t) - Yr * sp.sin(Om * t), Xr * sp.sin(Om * t) + Yr * sp.cos(Om * t)
Dr = lambda f: total_d(f, [Xr, Yr], [VXr, VYr], [AXr, AYr])
to_rot = {x: xr, y: yr, vx: Dr(xr), vy: Dr(yr), ax_: Dr(Dr(xr)), ay_: Dr(Dr(yr))}
EX_r, EY_r = euler_lagrange((m_ / 2 * (vx**2 + vy**2) - Vf(x, y)).subs(to_rot), [Xr, Yr], [VXr, VYr], [AXr, AYr])
Exr, Eyr = Ex.subs(to_rot), Ey.subs(to_rot)
assert sp.simplify((EX_r - (sp.diff(xr, Xr) * Exr + sp.diff(yr, Xr) * Eyr)).doit()) == 0
assert sp.simplify((EY_r - (sp.diff(xr, Yr) * Exr + sp.diff(yr, Yr) * Eyr)).doit()) == 0
print("SymPy: Newton (Thm 3.2), αT+βV (Thm 3.3), dF/dt and cL (Thm 3.4), polar E–L, and covariance (Thm 3.5) "
      "for polar and rotating coordinates verified.")

# Theorem 3.7: L = ½mv² − qΦ + q v·A gives m a = q(E + v × B) for arbitrary Φ(r, t), A(r, t).
Lem = m_ / 2 * sum(v**2 for v in V3) - qe * Phi(x, y, z, t) + qe * sum(v * Ai for v, Ai in zip(V3, A))
E_em = euler_lagrange(Lem, X3, V3, A3)
Efield = [-sp.diff(Phi(x, y, z, t), xi) - sp.diff(Ai, t) for xi, Ai in zip(X3, A)]
Bfield = [sp.diff(A[2], y) - sp.diff(A[1], z), sp.diff(A[0], z) - sp.diff(A[2], x), sp.diff(A[1], x) - sp.diff(A[0], y)]
vxB = [vy * Bfield[2] - vz * Bfield[1], vz * Bfield[0] - vx * Bfield[2], vx * Bfield[1] - vy * Bfield[0]]
assert all(sp.simplify(E_em[i] - (m_ * A3[i] - qe * (Efield[i] + vxB[i]))) == 0 for i in range(3))
# Example 4.6 / B6: uniform B ẑ, symmetric gauge vs Landau gauge differ by q dχ/dt with χ = Bxy/2.
L_sym = m_ / 2 * sum(v**2 for v in V3) + qe * B / 2 * (x * vy - y * vx)
L_lan = m_ / 2 * sum(v**2 for v in V3) + qe * B * x * vy
assert sp.simplify(L_lan - L_sym - qe * total_d(B * x * y / 2, X3, V3, A3)) == 0
E_sym, E_lan = euler_lagrange(L_sym, X3, V3, A3), euler_lagrange(L_lan, X3, V3, A3)
assert [sp.simplify(e) for e in E_sym] == [sp.simplify(e) for e in E_lan] == \
    [m_ * ax_ - qe * B * vy, m_ * ay_ + qe * B * vx, m_ * az_]
assert sp.simplify(sp.diff(L_sym, vx) - (m_ * vx - qe * B * y / 2)) == 0  # p_x ≠ m v_x
print("SymPy: Lorentz force from L (Thm 3.7); symmetric and Landau gauges give identical E–L (B6).")

# %% [markdown]
# ## Part 5: Cyclotron motion (Example 4.6): m = 2, q_e = 1, B = 1, start at origin with v = (1, 0, 0)

# %%
mc, qc, Bc = 2.0, 1.0, 1.0
wc = qc * Bc / mc
assert np.isclose(wc, 0.5)
# Integrate the SymPy E–L equations of Part 4 themselves: solve E_sym = 0 for the accelerations, then lambdify.
acc_sym = sp.solve(E_sym, A3, dict=True)[0]
acc_num = sp.lambdify((x, y, z, vx, vy, vz), [acc_sym[ai].subs({m_: mc, qe: qc, B: Bc}) for ai in A3], "numpy")
rhs_B = lambda t, s: [s[3], s[4], s[5], *acc_num(*s)]
Tc = 2 * np.pi / wc
sol = solve_ivp(rhs_B, (0, Tc), [0, 0, 0, 1, 0, 0], rtol=1e-11, atol=1e-12, dense_output=True)
tt = np.linspace(0, Tc, 400)
P = sol.sol(tt)
radius = np.hypot(P[0] - 0.0, P[1] - (-2.0))
assert np.allclose(radius, mc * 1.0 / (qc * Bc), atol=1e-8)          # R = mv/(qB) = 2, center (0, −2)
assert np.allclose(sol.y[:3, -1], 0, atol=1e-8)                        # back home after T_c = 4π
assert np.allclose(np.hypot(P[3], P[4]), 1.0)                           # speed constant: B does no work
quarter = sol.sol(Tc / 4)
assert np.allclose(quarter[:2], [2, -2], atol=1e-8)                     # clockwise seen from +z
print(f"Example 4.6: ω_c = {wc}, T_c = 4π = {Tc:.4f}, R = 2, center (0, −2); at T_c/4 the particle is at "
      f"({quarter[0]:.3f}, {quarter[1]:.3f})")
fig, ax = plt.subplots(figsize=(4.5, 4.5))
ax.plot(P[0], P[1]); ax.plot(0, -2, "k+"); ax.set_aspect("equal"); ax.grid(alpha=.3)
ax.set_title("Cyclotron orbit from L = ½mv² + q v·A"); fig.tight_layout()
fig.savefig(FIG / "cyclotron.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 6: Covariance, numerically (Theorem 3.5, Example 4.3)
#
# CM-1 §3.10(a): Cartesian pendulum with the tension as an unknown (a multiplier),
# T = m v²/ℓ − m g y/ℓ (CM-1 eq. 3.12). Against: the one-line θ equation from
# L = ½ m ℓ² θ̇² + m g ℓ cos θ. Same motion.

# %%
ell, th0 = 1.0, 1.0


def rhs_cart(t, s):
    x_, y_, vx_, vy_ = s
    T_over_m = (vx_**2 + vy_**2 - G * y_) / ell
    return [vx_, vy_, -T_over_m * x_ / ell, -T_over_m * y_ / ell - G]


tt = np.linspace(0, 5, 501)
cart = solve_ivp(rhs_cart, (0, 5), [ell * np.sin(th0), -ell * np.cos(th0), 0, 0], t_eval=tt, rtol=1e-12, atol=1e-12)
lag = solve_ivp(lambda t, s: [s[1], -G / ell * np.sin(s[0])], (0, 5), [th0, 0], t_eval=tt, rtol=1e-12, atol=1e-12)
th_cart = np.arctan2(cart.y[0], -cart.y[1])
dev = np.abs(th_cart - lag.y[0]).max()
assert dev < 1e-9, dev   # Example 4.3 quotes "about 1e-11"
print(f"Example 4.3: Cartesian-with-tension vs θ-Lagrangian pendulum, max |Δθ| over 5 s = {dev:.1e}")
assert np.isclose(G * np.sin(np.pi / 6), 4.905)
print(f"Example 4.3: inclined plane at 30°: s̈ = g sin α = {G * np.sin(np.pi / 6):.3f} m/s²")
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(tt, lag.y[0], lw=4, alpha=.35, label="θ from L(θ, θ̇): one line")
ax.plot(tt, th_cart, "k--", lw=1, label="atan2(x, −y) from CM-1's Cartesian + tension")
ax.set_xlabel("t (s)"); ax.set_ylabel("θ (rad)"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "pendulum_covariance.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: Numbers quoted in solutions.md (A1, A2, A4, B7, C2, C3)

# %%
# A1: m = 2, x(0) = 0, x(2) = 3. Straight line vs x = 3t²/4.
S_a1 = action(lambda t: 1.5 * t, lambda t: 1.5 + 0 * t, 0, 2, 2.0, lambda q: 0 * q)
S_a1c = action(lambda t: .75 * t**2, lambda t: 1.5 * t, 0, 2, 2.0, lambda q: 0 * q)
assert np.isclose(S_a1, 4.5) and np.isclose(S_a1c, 6.0) and np.isclose(S_a1c - S_a1, 1.5)
# A2: Atwood, L = ½(m1+m2)ẏ² − (m1 − m2) g y  =>  ÿ = (m2 − m1) g/(m1 + m2), CM-1 (3.16).
m1, m2 = sp.symbols("m1 m2", positive=True)
E_at = euler_lagrange((m1 + m2) / 2 * vy**2 - (m1 - m2) * g_ * y, [y], [vy], [ay_])[0]
assert sp.simplify(sp.solve(E_at, ay_)[0] - (m2 - m1) * g_ / (m1 + m2)) == 0
# A4: which added terms change the E–L equation of L = ½ m q̇²?
q, v, a = sp.symbols("q v a", real=True)
extra = {"3 q q̇": 3 * q * v, "t q̇": t * v, "q̇ cos q": v * sp.cos(q), "q + t q̇": q + t * v}
change = {k: sp.simplify(euler_lagrange(G_, [q], [v], [a])[0]) for k, G_ in extra.items()}
assert change == {"3 q q̇": 0, "t q̇": 1, "q̇ cos q": 0, "q + t q̇": 0}, change
# B7 reference: L' = L + dF/dt with F = q² sin t changes S by F(q2,t2) − F(q1,t1) = sin 2 on every path.
Lp = lambda q, qd, t: 0.5 * qd**2 - 0.5 * q**2 + 2 * q * qd * np.sin(t) + q**2 * np.cos(t)
paths = [(lambda t: np.sin(t) / np.sin(2), lambda t: np.cos(t) / np.sin(2)),
         (lambda t: t / 2, lambda t: 0.5 + 0 * t), (lambda t: t**2 / 4, lambda t: t / 2)]
for qf, qdf in paths:
    S, Sp = action(qf, qdf, 0, 2, 1.0, osc_V), gl_integrate(lambda t: Lp(qf(t), qdf(t), t), 0, 2)
    assert np.isclose(Sp - S, np.sin(2))
print(f"A1: S = 4.5, 6.0 (diff 1.5); A2, A4 by SymPy; B7: S' − S = sin 2 = {np.sin(2):.6f} on all three paths")

# C2: Kepler, m = k = 1. Circle r = 1 over one period 2π vs the a = 1 ellipse launched from (1, 0) at angle γ.
kep = lambda t, s: [s[2], s[3], -s[0] / np.hypot(s[0], s[1])**3, -s[1] / np.hypot(s[0], s[1])**3,
                    0.5 * (s[2]**2 + s[3]**2) + 1 / np.hypot(s[0], s[1])]   # last slot accumulates S
for gam in (0.0, 0.3):
    kep_sol = solve_ivp(kep, (0, 2 * np.pi), [1, 0, np.sin(gam), np.cos(gam), 0], rtol=1e-12, atol=1e-12,
                        dense_output=True)
    s = kep_sol.y[:, -1]
    assert np.allclose(s[:2], [1, 0], atol=1e-7) and np.isclose(s[4], 3 * np.pi)
    rr = np.hypot(*kep_sol.sol(np.linspace(0, 2 * np.pi, 4001))[:2])
    e_orbit = (rr.max() - rr.min()) / (rr.max() + rr.min())          # e = (r_a − r_p)/(r_a + r_p), CM-5 (3.13)
    assert np.isclose(e_orbit, np.sin(gam), atol=1e-6) and np.isclose(0.5 * (rr.max() + rr.min()), 1.0, atol=1e-6)
    print(f"C2: γ = {gam}: measured e = {e_orbit:.4f} = sin γ, a = 1, returns to (1, 0) at t = 2π, S = {s[4]:.6f} = 3π")
S_rad = lambda e: gl_integrate(lambda t: 0.5 * ((0.5 * e * np.cos(t / 2))**2 + (1 + e * np.sin(t / 2))**2)
                               + 1 / (1 + e * np.sin(t / 2)), 0, 2 * np.pi)
e = 1e-3
curv = (S_rad(e) + S_rad(-e) - 2 * S_rad(0)) / (2 * e**2)
assert np.isclose(S_rad(0), 3 * np.pi) and np.isclose(curv, 13 * np.pi / 8, rtol=1e-5)
print(f"C2: radial perturbation r = 1 + ε sin(t/2): S = 3π + {curv:.4f} ε² + O(ε³)  (13π/8 = {13 * np.pi / 8:.4f})")
# C3: L = e^{γt}(½ m q̇² − ½ k q²) gives e^{γt}(m q̈ + m γ q̇ + k q).
gm = sp.symbols("gamma", positive=True)
E_damp = euler_lagrange(sp.exp(gm * t) * (m_ * v**2 / 2 - k_ * q**2 / 2), [q], [v], [a])[0]
assert sp.simplify(E_damp - sp.exp(gm * t) * (m_ * a + m_ * gm * v + k_ * q)) == 0
print("C3: E–L of e^{γt}(T − V) is e^{γt}(m q̈ + mγ q̇ + k q)  (SymPy)")

# %% [markdown]
# ## Your turn
#
# 1. (B7) Add F(q, t) = q² sin t to the oscillator Lagrangian and check that
#    S' − S is the same number for every path with the same endpoints. Part 7
#    has a reference; write yours first. Then discretize L' and minimize it.
# 2. Push the pendulum interval in Part 1 from τ = 0.8 s to τ = 1.5 s (past half
#    a small-oscillation period, π√(ℓ/g) ≈ 1.003 s). The minimizer and `shoot`
#    now disagree, and neither is "wrong": scan q̇(0) over [−20, 20] and you find
#    three true paths (q̇(0) ≈ −5.62, 0.57, 5.68). Unlike the oscillator, the
#    pendulum's V is bounded, so S_h is bounded below and L-BFGS-B converges, to
#    the large swing q̇(0) ≈ 5.68. Build the Hessian of S_h at each path, as Part 3
#    does: the small-amplitude path (q̇(0) ≈ 0.57) has one negative eigenvalue (a
#    saddle), the two large swings have none. Then continue Part 1's τ = 0.8
#    solution to τ = 1.5 in small steps of τ: it grows into the q̇(0) ≈ 5.68 swing.
# 3. (C1) The brachistochrone as mechanics: treat x as time and y as the
#    coordinate, with "Lagrangian" √((1 + y'²)/(2 g y)). Minimize its discrete
#    "action" between (0, small) and (1, 0.65) and compare with the cycloid of F6
#    Example 4.2. What is the conserved "energy"?
# 4. (C2) Replace the circle in Part 7 by a slightly non-Keplerian path with
#    the same endpoints (perturb θ(t) too, not just r) and confirm the action
#    goes up. Then find a perturbation direction for which it does not.
# 5. Add a uniform E = E₀ x̂ to Part 5 (Φ = −E₀ x). Predict the drift velocity
#    E × B / B² and measure it from the orbit.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
