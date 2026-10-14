# %% [markdown]
# # CM-1 Project: Conservation Laws, Constraints, and Three Integrators
#
# This script checks every number in the CM-1 notes and makes §3.11 visible.
# Parts: (1) Euler vs RK4 vs symplectic Euler on the nonlinear pendulum for
# 100 periods, plus the Euler energy-growth laws (3.17), (3.18); (2) the pendulum
# phase portrait with its separatrix; (3) period vs amplitude: numerics vs the
# elliptic integral (4.1) vs the series (4.2); (4) projectiles with linear and
# quadratic drag; (5) the driven damped oscillator's resonance curve (3.10);
# (6) the bead on the parabolic wire, eq. (3.14), with the normal force;
# (7) the Atwood machine solved the Newtonian way at every step; (8) the
# exploding shell and the rocket; (9) SymPy checks of the derivations;
# (10) the remaining numbers quoted in solutions.md.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import ellipk, gamma as Gamma
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
g = 9.81

# %% [markdown]
# ## Part 1: Three integrators on the nonlinear pendulum (§3.11)
#
# x'' = a(x) as the system (x, v)' = (v, a(x)). Euler inflates phase-space area
# by 1 + h^2 ω^2 per step on the harmonic oscillator, eq. (3.17), and gains
# ½h²(a² + V''v²) per step in general, eq. (3.18); symplectic
# Euler has Jacobian determinant exactly 1 (exercise C2); RK4 is accurate but
# not area-preserving. Pendulum: θ'' = -(g/l) sin θ, energy per unit mass
# E = ½ l² θ'² + g l (1 - cos θ).
#
# **Build it yourself first:** write `euler_step`, `rk4_step`,
# `symplectic_euler_step` for a generic a(x) before reading these.

# %%
def euler_step(a, x, v, h):
    """Explicit Euler for x'' = a(x): x_{n+1} = x_n + h v_n, v_{n+1} = v_n + h a(x_n). §3.11."""
    return x + h * v, v + h * a(x)

def symplectic_euler_step(a, x, v, h):
    """Symplectic Euler: update v, then x with the NEW v. Area-preserving (C2). §3.11."""
    v_new = v + h * a(x)
    return x + h * v_new, v_new

def rk4_step(a, x, v, h):
    """Classical fourth-order Runge–Kutta on (x, v)' = (v, a(x))."""
    k1x, k1v = v, a(x)
    k2x, k2v = v + 0.5 * h * k1v, a(x + 0.5 * h * k1x)
    k3x, k3v = v + 0.5 * h * k2v, a(x + 0.5 * h * k2x)
    k4x, k4v = v + h * k3v, a(x + h * k3x)
    return (x + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x),
            v + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v))

def integrate(step, a, x0, v0, h, n):
    """Run `step` n times from (x0, v0); return arrays x, v of length n+1."""
    xs = np.empty(n + 1); vs = np.empty(n + 1)
    xs[0], vs[0] = x, v = x0, v0
    for i in range(1, n + 1):
        x, v = step(a, x, v, h)
        xs[i], vs[i] = x, v
    return xs, vs

# (3.17) exactly, on the harmonic oscillator: E_N = (1 + h^2 ω^2)^N E_0.
w_sho, h_sho, N_sho = 3.0, 0.01, 500
xs, vs = integrate(euler_step, lambda x: -w_sho**2 * x, 1.0, 0.0, h_sho, N_sho)
E_sho = 0.5 * vs**2 + 0.5 * w_sho**2 * xs**2
assert np.allclose(E_sho, E_sho[0] * (1 + h_sho**2 * w_sho**2) ** np.arange(N_sho + 1), rtol=1e-10)
print(f"(3.17) verified: Euler multiplies SHO energy by exactly {1 + h_sho**2 * w_sho**2:.6f} per step.")
# The step used in §3.11 and below: h = T0/1000, so hω = 2π/1000. SHO with ω = 1 for 100 periods.
hw = 2 * np.pi / 1000
xs, vs = integrate(euler_step, lambda x: -x, 1.0, 0.0, hw, 100_000)
E_r = vs**2 + xs**2                                                    # E/E0, since E0 = 1/2
assert np.isclose(hw, 0.006283, atol=1e-6) and np.isclose(hw**2, 3.948e-5, rtol=1e-3)
assert np.isclose(E_r[10_000], 1.484, atol=1e-3) and np.isclose(E_r[-1], 51.8, atol=0.05)
print(f"§3.11: hω = {hw:.6f}, factor 1 + {hw**2:.3e}/step; SHO E/E0 = {E_r[10_000]:.3f} after 10 periods, {E_r[-1]:.1f} after 100")

# Nonlinear pendulum, 100 periods, two step sizes.
l_p = 1.0
theta0 = 1.0
a_pend = lambda th: -(g / l_p) * math.sin(th)
E_pend = lambda th, om: 0.5 * l_p**2 * om**2 + g * l_p * (1 - np.cos(th))
T0 = 2 * np.pi * np.sqrt(l_p / g)
T_exact = 4 * np.sqrt(l_p / g) * ellipk(np.sin(theta0 / 2) ** 2)       # eq. (4.1)
E0 = E_pend(theta0, 0.0)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
results = {}
for ax, steps_per_period in zip(axes, (1000, 100)):
    h = T0 / steps_per_period
    n = int(round(100 * T_exact / h))
    t = h * np.arange(n + 1)
    for name, step in (("Euler", euler_step), ("RK4", rk4_step), ("symplectic Euler", symplectic_euler_step)):
        xs, vs = integrate(step, a_pend, theta0, 0.0, h, n)
        relerr = (E_pend(xs, vs) - E0) / E0
        results[(steps_per_period, name)] = relerr
        if (steps_per_period, name) == (1000, "Euler"):
            th_eu = xs
        ax.semilogy(t / T_exact, np.abs(relerr) + 1e-17, lw=0.8, label=name)
    ax.set_title(f"h = T0/{steps_per_period}: |E - E0| / E0 over 100 periods")
    ax.set_xlabel("time / period"); ax.grid(alpha=.3); ax.legend()
axes[0].set_ylabel("relative energy error")
fig.tight_layout(); fig.savefig(FIG / "pendulum_energy_three_integrators.png", dpi=120); plt.close(fig)

# Fine step (1000/period). Euler: slower than geometric growth, (3.18), since V'' ∝ cos θ softens at
# large swings, but it still crosses the separatrix E = 2gl and turns oscillation into rotation.
eu, rk, sy = (results[(1000, k)] for k in ("Euler", "RK4", "symplectic Euler"))
n10 = len(eu) // 10
per_sep = np.argmax(E0 * (1 + eu) > 2 * g * l_p) * (T0 / 1000) / T_exact     # first step past the separatrix
print(f"Pendulum Euler: E/E0 = {1 + eu[n10]:.3f} after 10 periods; crosses the separatrix (E/E0 = {2*g*l_p/E0:.2f}) "
      f"at period {per_sep:.1f}; after 100 periods E/E0 = {1 + eu[-1]:.2f} and θ = {th_eu[-1]:.0f} rad (spinning)")
print(f"RK4 max |dE|/E0 = {np.abs(rk).max():.2e}; symplectic max |dE|/E0 = {np.abs(sy).max():.2e}")
assert np.isclose(1 + eu[n10], 1.41, atol=0.01), "Euler looks almost fine at 10 periods"
assert 60 < per_sep < 62.5 and np.isclose(1 + eu[-1], 5.9, atol=0.05), "Euler pushes the bob over the top"
assert np.abs(th_eu).max() > 100 * np.pi, "by period 100 the Euler 'pendulum' has made dozens of full turns"
assert np.abs(rk).max() < 1e-6
assert np.abs(sy).max() < 2e-2
drift = abs(sy[-n10:].mean() - sy[:n10].mean())
assert drift < 0.1 * np.abs(sy).max(), "symplectic energy error must not drift"
# Coarse step (100/period): RK4 drift is now visible but still small; symplectic still bounded.
rk_c, sy_c = results[(100, "RK4")], results[(100, "symplectic Euler")]
print(f"h = T0/100: RK4 max |dE|/E0 = {np.abs(rk_c).max():.2e} (final {rk_c[-1]:+.2e}); symplectic max {np.abs(sy_c).max():.2e}")
assert 1e-6 < np.abs(rk_c).max() < 1e-2, np.abs(rk_c).max()
assert np.abs(sy_c).max() < 0.2
assert abs(sy_c[-n10 // 10:].mean() - sy_c[: n10 // 10].mean()) < 0.1 * np.abs(sy_c).max()
# Hand RK4 agrees with solve_ivp over 10 periods.
h = T0 / 1000; n = int(round(10 * T_exact / h))
xs, vs = integrate(rk4_step, a_pend, theta0, 0.0, h, n)
ref = solve_ivp(lambda t, y: [y[1], a_pend(y[0])], [0, n * h], [theta0, 0.0], method="DOP853", rtol=1e-12, atol=1e-12)
assert abs(xs[-1] - ref.y[0, -1]) < 1e-6
print("Part 1 passed: Euler turns oscillation into rotation, RK4 drifts slightly, symplectic Euler stays bounded.")

# %% [markdown]
# ## Part 2: Phase portrait of the pendulum (§2, §3.6)
#
# Level curves of E(θ, θ') = ½ θ'² + (g/l)(1 - cos θ) (per unit m l²). The
# separatrix is the level E = 2g/l through the saddles at θ = ±π.

# %%
TH, OM = np.meshgrid(np.linspace(-2 * np.pi, 2 * np.pi, 601), np.linspace(-9, 9, 401))
EE = 0.5 * OM**2 + (g / l_p) * (1 - np.cos(TH))
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.contour(TH, OM, EE, levels=np.linspace(0.5, 60, 24), colors="C0", linewidths=0.7)
ax.contour(TH, OM, EE, levels=[2 * g / l_p], colors="r", linewidths=2)
ax.plot([-np.pi, np.pi], [0, 0], "rx", ms=8); ax.plot([-2 * np.pi, 0, 2 * np.pi], [0, 0, 0], "ko", ms=5)
ax.set_xlabel("θ"); ax.set_ylabel("θ-dot"); ax.set_title("Pendulum phase portrait; red = separatrix E = 2g/l")
ax.set_xticks([-2 * np.pi, -np.pi, 0, np.pi, 2 * np.pi]); ax.set_xticklabels(["-2π", "-π", "0", "π", "2π"]); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "pendulum_phase_portrait.png", dpi=120); plt.close(fig)
print("Part 2: phase portrait saved.")

# %% [markdown]
# ## Part 3: Period vs amplitude (Example 4.2, exercise C1)
#
# Exact: T = 4 sqrt(l/g) K(k), k = sin(θ0/2), eq. (4.1). NOTE: `ellipk(m)`
# takes the parameter m = k², not k. Series: T/T0 = 1 + θ0²/16 + 11 θ0⁴/3072.

# %%
def pendulum_period_elliptic(theta0, l=1.0):
    """Eq. (4.1): T = 4 sqrt(l/g) K(sin²(θ0/2)); ellipk takes m = k²."""
    return 4 * np.sqrt(l / g) * ellipk(np.sin(theta0 / 2) ** 2)

def pendulum_period_numeric(theta0, l=1.0):
    """Period from solve_ivp: θ = 0 crossings occur at T/4 and 3T/4, so T = 2 (t2 - t1)."""
    T_guess = pendulum_period_elliptic(theta0, l)
    sol = solve_ivp(lambda t, y: [y[1], -(g / l) * np.sin(y[0])], [0, 1.1 * T_guess], [theta0, 0.0],
                    events=lambda t, y: y[0], rtol=1e-11, atol=1e-13)
    te = sol.t_events[0]
    return 2 * (te[1] - te[0])

deg = np.arange(5, 180, 5)
th0 = np.radians(deg)
T_ell = pendulum_period_elliptic(th0, l_p)
T_num = np.array([pendulum_period_numeric(t, l_p) for t in th0])
T_1 = T0 * (1 + th0**2 / 16)
T_2 = T0 * (1 + th0**2 / 16 + 11 * th0**4 / 3072)
assert np.allclose(T_num, T_ell, rtol=1e-7)
r90 = pendulum_period_elliptic(np.pi / 2) / T0; r30 = pendulum_period_elliptic(np.pi / 6) / T0
assert np.isclose(r90, 1.18034, atol=1e-5) and np.isclose(r30, 1.01741, atol=1e-5)
assert np.isclose(1 + (np.pi / 2) ** 2 / 16, 1.15421, atol=1e-5)
assert np.isclose(1 + (np.pi / 2) ** 2 / 16 + 11 * (np.pi / 2) ** 4 / 3072, 1.17601, atol=1e-5)
assert np.isclose(ellipk(0.5), 1.85407, atol=1e-5) and np.isclose(np.sin(np.pi / 12) ** 2, 0.06699, atol=1e-5)
assert np.isclose(1 + (np.pi / 6) ** 2 / 16, 1.01713, atol=1e-5)                                 # 30°, one term
s1 = 1 + (np.pi / 2) ** 2 / 16; s2 = s1 + 11 * (np.pi / 2) ** 4 / 3072                          # 90°, Ex. 4.2 and C1
assert np.isclose((np.pi / 2) ** 4, 6.0881, atol=1e-4) and np.isclose(s2 - s1, 0.02180, atol=1e-5)
assert np.isclose(1 - s1 / r90, 0.0221, atol=1e-4) and np.isclose(1 - s2 / r90, 0.0037, atol=1e-4)
print(f"Example 4.2: T/T0 at 90° = {r90:.5f} (exact), {1 + (np.pi/2)**2/16:.5f} (1 term), "
      f"{1 + (np.pi/2)**2/16 + 11*(np.pi/2)**4/3072:.5f} (2 terms); at 30° = {r30:.5f}")

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(deg, T_ell / T0, "k-", label="elliptic integral (4.1)")
ax.plot(deg, T_num / T0, "o", ms=4, label="solve_ivp")
ax.plot(deg, T_1 / T0, "--", label="1 + θ0²/16")
ax.plot(deg, T_2 / T0, ":", label="+ 11θ0⁴/3072 (C1)")
ax.set_xlabel("θ0 (degrees)"); ax.set_ylabel("T / T0"); ax.set_ylim(0.95, 2.7); ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "pendulum_period_vs_amplitude.png", dpi=120); plt.close(fig)
print("Part 3 passed: numerical period matches the elliptic integral at every amplitude.")

# %% [markdown]
# ## Part 4: Projectiles with drag (§3.7, Example 4.3, exercise A3)
#
# Linear drag closed form (3.9) vs solve_ivp; quadratic drag vs no drag (no
# closed form); the tanh law for a falling body with quadratic drag.

# %%
m_lin, b_lin = 0.1, 0.05
tau = m_lin / b_lin; v_inf = g * tau
v0, ang = 30.0, np.radians(45)
vx0, vy0 = v0 * np.cos(ang), v0 * np.sin(ang)
x_lin = lambda t: vx0 * tau * (1 - np.exp(-t / tau))                       # (3.9)
y_lin = lambda t: (vy0 + v_inf) * tau * (1 - np.exp(-t / tau)) - v_inf * t   # (3.9)
t_f = brentq(y_lin, 1.0, 20.0)
t_top = tau * np.log((vy0 + v_inf) / v_inf)
assert np.isclose(tau, 2.0) and np.isclose(v_inf, 19.62)
assert np.isclose(vx0 * tau, 42.426, atol=1e-3)
assert np.isclose(t_f, 3.4032, atol=1e-4) and np.isclose(x_lin(t_f), 34.688, atol=1e-3)
assert np.isclose(t_top, 1.4659, atol=1e-4) and np.isclose(y_lin(t_top), 13.666, atol=1e-3)
assert np.isclose(v0**2 / g, 91.743, atol=1e-3) and np.isclose(2 * vy0 / g, 4.3248, atol=1e-4)
assert np.isclose((vy0 + v_inf) * tau, 81.666, atol=1e-3) and np.isclose((vy0 + v_inf) / v_inf, 2.0812, atol=1e-4)
assert np.isclose(vy0**2 / (2 * g), 22.94, atol=1e-2) and np.isclose(1 - x_lin(t_f) / (v0**2 / g), 0.62, atol=5e-3)
sol = solve_ivp(lambda t, y: [y[2], y[3], -y[2] / tau, -g - y[3] / tau], [0, t_f], [0, 0, vx0, vy0],
                rtol=1e-11, atol=1e-13, dense_output=True)
tt = np.linspace(0, t_f, 200)
assert np.allclose(sol.sol(tt)[0], x_lin(tt), atol=1e-6) and np.allclose(sol.sol(tt)[1], y_lin(tt), atol=1e-6)
print(f"Example 4.3: flight time {t_f:.4f} s, range {x_lin(t_f):.3f} m (no drag {v0**2/g:.2f} m), "
      f"max height {y_lin(t_top):.3f} m at t = {t_top:.4f} s, asymptote x = {vx0*tau:.3f} m")

# Quadratic drag (baseball, Taylor's c ≈ 1.225e-3 kg/m) vs no drag.
m_q, c_q = 0.145, 1.225e-3
def quad_drag(t, y):
    vx, vy = y[2], y[3]; s = np.hypot(vx, vy)
    return [vx, vy, -(c_q / m_q) * s * vx, -g - (c_q / m_q) * s * vy]
hit = lambda t, y: y[1]; hit.terminal = True; hit.direction = -1
sq = solve_ivp(quad_drag, [0, 20], [0, 1e-9, vx0, vy0], events=hit, rtol=1e-10, atol=1e-12, dense_output=True)
sn = solve_ivp(lambda t, y: [y[2], y[3], 0, -g], [0, 20], [0, 1e-9, vx0, vy0], events=hit, rtol=1e-10, atol=1e-12, dense_output=True)
range_q, range_n = sq.y[0, -1], sn.y[0, -1]
assert np.isclose(range_n, v0**2 / g, rtol=1e-6) and range_q < range_n
print(f"Quadratic drag baseball at 30 m/s, 45°: range {range_q:.2f} m vs {range_n:.2f} m without drag")

# A3: falling from rest with quadratic drag, s(t) = s_inf tanh(t/τ).
s_inf = np.sqrt(m_q * g / c_q); tau_q = s_inf / g
sf = solve_ivp(lambda t, y: [g - (c_q / m_q) * y[0] ** 2], [0, 12], [0.0], rtol=1e-11, atol=1e-13, dense_output=True)
tt = np.linspace(0, 12, 100)
assert np.allclose(sf.sol(tt)[0], s_inf * np.tanh(tt / tau_q), atol=1e-7)
assert np.isclose(s_inf, 34.08, atol=0.01) and np.isclose(tau_q * np.arctanh(0.9), 5.114, atol=1e-3)
assert np.isclose(tau_q, 3.474, atol=1e-3) and np.isclose(np.arctanh(0.9), 1.4722, atol=1e-4) and np.isclose(5 * g, 49.05)
print(f"A3: terminal speed {s_inf:.2f} m/s, τ = {tau_q:.3f} s, 90% of terminal at t = {tau_q*np.arctanh(0.9):.3f} s")

fig, ax = plt.subplots(figsize=(8, 4.5))
tn = np.linspace(0, sn.t[-1], 300); tq = np.linspace(0, sq.t[-1], 300); tl = np.linspace(0, t_f, 300)
ax.plot(sn.sol(tn)[0], sn.sol(tn)[1], label="no drag")
ax.plot(x_lin(tl), y_lin(tl), label="linear drag, τ = 2 s (Example 4.3)")
ax.plot(sq.sol(tq)[0], sq.sol(tq)[1], label="quadratic drag (baseball)")
ax.axvline(vx0 * tau, ls=":", color="C1", label="linear-drag asymptote x = v_x0 τ")
ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_ylim(0, 25); ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "projectiles_with_drag.png", dpi=120); plt.close(fig)
print("Part 4 passed.")

# %% [markdown]
# ## Part 5: The driven damped oscillator (§3.8, Example 4.4, exercise B5)
#
# m x'' + b x' + k x = F0 cos ωt. Integrate from rest for 20 s (transient
# ~ e^{-t}), fit the last whole periods to A cos(ωt - δ), compare with (3.10).

# %%
m_o, k_o, b_o, F0 = 1.0, 100.0, 2.0, 1.0
w0 = np.sqrt(k_o / m_o); gam = b_o / m_o; Q = w0 / gam
A_an = lambda w: (F0 / m_o) / np.sqrt((w0**2 - w**2) ** 2 + gam**2 * w**2)       # (3.10)
d_an = lambda w: np.arctan2(gam * w, w0**2 - w**2)                               # (3.10), δ in [0, π]

def steady_state(w, t_end=20.0):
    """Numerical steady-state amplitude and phase lag of the driven oscillator at drive frequency w."""
    rhs = lambda t, y: [y[1], (F0 * np.cos(w * t) - b_o * y[1] - k_o * y[0]) / m_o]
    sol = solve_ivp(rhs, [0, t_end], [0.0, 0.0], rtol=1e-10, atol=1e-13, dense_output=True)
    Tw = 2 * np.pi / w; nper = max(1, int(2.0 // Tw)); t = np.linspace(t_end - nper * Tw, t_end, 4001)
    x = sol.sol(t)[0]
    c = 2 / (nper * Tw) * np.trapezoid(x * np.cos(w * t), t)      # x = A cos δ cos ωt + A sin δ sin ωt
    s = 2 / (nper * Tw) * np.trapezoid(x * np.sin(w * t), t)
    return np.hypot(c, s), np.arctan2(s, c)

assert np.isclose(w0, 10.0) and np.isclose(gam, 2.0) and np.isclose(Q, 5.0)
assert np.isclose(A_an(w0), 0.05) and np.isclose(d_an(w0), np.pi / 2) and np.isclose(Q * F0 / k_o, 0.05)
assert np.isclose(A_an(8.0), 0.02538, atol=1e-5) and np.isclose(d_an(8.0), 0.4182, atol=1e-4)
assert np.isclose(A_an(12.0), 0.01995, atol=1e-5) and np.isclose(d_an(12.0), 2.6422, atol=1e-4)
w_max = np.sqrt(w0**2 - gam**2 / 2)
assert np.isclose(w_max, 9.8995, atol=1e-4) and np.isclose(A_an(w_max), 0.050252, atol=1e-6)
assert np.isclose(A_an(w_max), (F0 / m_o) / (gam * np.sqrt(w0**2 - gam**2 / 4)))
assert np.isclose(np.sqrt(w0**2 - gam**2 / 4), 9.950, atol=1e-3)                                # ω1
assert np.isclose(np.exp(-5), 0.0067, atol=1e-4) and np.isclose(np.exp(-15), 3e-7, atol=1e-8)   # transient e^{-γt/2}
P_an = lambda w: 0.5 * b_o * w**2 * A_an(w) ** 2
w_lo = brentq(lambda w: P_an(w) - 0.5 * P_an(w0), 1, w0); w_hi = brentq(lambda w: P_an(w) - 0.5 * P_an(w0), w0, 30)
assert np.isclose(w_hi - w_lo, gam, atol=1e-8), "power FWHM should equal γ (B5c)"
print(f"Example 4.4: A(ω0) = {A_an(w0):.4f} m = Q × F0/k; ω_max = {w_max:.4f}, A_max = {A_an(w_max):.6f}; "
      f"A(8) = {A_an(8):.5f}, δ(8) = {np.degrees(d_an(8)):.2f}°; A(12) = {A_an(12):.5f}, δ(12) = {np.degrees(d_an(12)):.2f}°")
print(f"Half-power frequencies {w_lo:.4f}, {w_hi:.4f} rad/s: FWHM = {w_hi - w_lo:.4f} = γ")

ws = np.linspace(2, 20, 37)
num = np.array([steady_state(w) for w in ws])
assert np.allclose(num[:, 0], A_an(ws), rtol=2e-3), np.abs(num[:, 0] / A_an(ws) - 1).max()
assert np.allclose(num[:, 1], d_an(ws), atol=5e-3)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
wf = np.linspace(2, 20, 400)
ax1.plot(wf, A_an(wf), "k-", label="(3.10)"); ax1.plot(ws, num[:, 0], "o", ms=4, label="solve_ivp steady state")
ax1.axhline(F0 / k_o, ls=":", color="gray", label="static F0/k"); ax1.set_xlabel("ω (rad/s)"); ax1.set_ylabel("amplitude (m)")
ax1.set_title(f"Resonance: Q = {Q:.0f}, peak at ω = {w_max:.2f}"); ax1.grid(alpha=.3); ax1.legend()
ax2.plot(wf, np.degrees(d_an(wf)), "k-"); ax2.plot(ws, np.degrees(num[:, 1]), "o", ms=4)
ax2.set_xlabel("ω (rad/s)"); ax2.set_ylabel("phase lag δ (deg)"); ax2.set_title("Phase lag: 0 → 90° at ω0 → 180°"); ax2.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "driven_oscillator_resonance.png", dpi=120); plt.close(fig)
print("Part 5 passed: numerical steady state matches (3.10) in amplitude and phase.")

# %% [markdown]
# ## Part 6: Bead on the parabolic wire (§3.10(b), Example 4.5)
#
# x'' = -2kx (g + 2k x'²) / (1 + 4k²x²), eq. (3.14).
# N = m [g + 2k (x'² + x x'')] sqrt(1 + 4k²x²); at the bottom N = m (g + 2k x'²), (3.15).

# %%
k_w, m_b, x0_b = 0.5, 0.2, 1.0
xdd_bead = lambda x, v: -2 * k_w * x * (g + 2 * k_w * v**2) / (1 + 4 * k_w**2 * x**2)       # (3.14)

def normal_force(x, v):
    """N along the trajectory, from step 6 of §3.10(b)."""
    return m_b * (g + 2 * k_w * (v**2 + x * xdd_bead(x, v))) * np.sqrt(1 + 4 * k_w**2 * x**2)

bottom = lambda t, y: y[0]
sb = solve_ivp(lambda t, y: [y[1], xdd_bead(y[0], y[1])], [0, 3.0], [x0_b, 0.0], events=bottom,
               rtol=1e-12, atol=1e-14, dense_output=True)
t_bot = sb.t_events[0][0]; xb, vb = sb.sol(t_bot)
E_bead = lambda x, v: 0.5 * m_b * (1 + 4 * k_w**2 * x**2) * v**2 + m_b * g * k_w * x**2
tt = np.linspace(0, 3, 600); xx, vv = sb.sol(tt)
assert np.allclose(E_bead(xx, vv), E_bead(x0_b, 0.0), rtol=1e-9), "normal force must do no work"
assert np.isclose(vb**2, 2 * g * k_w * x0_b**2, rtol=1e-8) and np.isclose(vb**2, 9.81, atol=1e-6)
assert np.isclose(normal_force(xb, vb), 3.924, atol=1e-6) and np.isclose(normal_force(xb, vb), 2 * m_b * g)
assert np.isclose(normal_force(x0_b, 0.0), 1.3873, atol=1e-4) and np.isclose(normal_force(x0_b, 0.0), m_b * g / np.sqrt(2))
assert np.isclose(abs(vb), 3.132, atol=1e-3) and np.isclose(xdd_bead(x0_b, 0.0), -4.905)
# Independent check: N - mg cos α = m v_tot² / ρ with ρ = (1 + f'²)^{3/2} / f'' and cos α = 1/sqrt(1 + f'²).
fp = 2 * k_w * xx; v_tot2 = vv**2 * (1 + fp**2); rho = (1 + fp**2) ** 1.5 / (2 * k_w)
assert np.allclose(normal_force(xx, vv) - m_b * g / np.sqrt(1 + fp**2), m_b * v_tot2 / rho, rtol=1e-8)
print(f"Example 4.5: bottom reached at t = {t_bot:.4f} s with v = {abs(vb):.4f} m/s; N_bottom = {normal_force(xb, vb):.4f} N = 2mg;"
      f" N at release = {normal_force(x0_b, 0.0):.4f} N = mg/√2")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
ax1.plot(tt, xx); ax1.set_xlabel("t (s)"); ax1.set_ylabel("x (m)"); ax1.set_title("Bead on y = 0.5 x²"); ax1.grid(alpha=.3)
ax2.plot(tt, normal_force(xx, vv)); ax2.axhline(m_b * g, ls=":", color="gray", label="mg")
ax2.set_xlabel("t (s)"); ax2.set_ylabel("N (N)"); ax2.set_title("Normal force along the trajectory"); ax2.grid(alpha=.3); ax2.legend()
fig.tight_layout(); fig.savefig(FIG / "bead_on_parabola.png", dpi=120); plt.close(fig)
print("Part 6 passed.")

# %% [markdown]
# ## Part 7: The Atwood machine, the Newtonian way (§3.10(c), Example 4.6)
#
# At every step solve the 3×3 linear system for (y1'', y2'', T):
# m1 y1'' - T = -m1 g,  m2 y2'' - T = -m2 g,  y1'' + y2'' = 0.

# %%
m1, m2 = 1.0, 1.5

def atwood_unknowns(m1, m2):
    """Solve for (a1, a2, T) from Newton's equations plus the string constraint."""
    A = np.array([[m1, 0.0, -1.0], [0.0, m2, -1.0], [1.0, 1.0, 0.0]])
    return np.linalg.solve(A, [-m1 * g, -m2 * g, 0.0])

def rk4_vec(f, t, y, h):
    k1 = f(t, y); k2 = f(t + h / 2, y + h / 2 * k1); k3 = f(t + h / 2, y + h / 2 * k2); k4 = f(t + h, y + h * k3)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

def atwood_rhs(t, y):
    a1, a2, T = atwood_unknowns(m1, m2)
    return np.array([y[1], a1, y[3], a2])

a1, a2, T_at = atwood_unknowns(m1, m2)
assert np.isclose(a1, (m2 - m1) * g / (m1 + m2)) and np.isclose(T_at, 2 * m1 * m2 * g / (m1 + m2))   # (3.16)
assert np.isclose(a1, 1.962, atol=1e-6) and np.isclose(T_at, 11.772, atol=1e-6) and np.isclose(a2, -a1)
assert m1 * g < T_at < m2 * g and np.isclose(m2 * g, 14.715)
y = np.array([0.0, 0.0, 0.0, 0.0]); h = 0.01
for i in range(200):
    y = rk4_vec(atwood_rhs, i * h, y, h)
assert np.isclose(y[0], 0.5 * a1 * 2.0**2, rtol=1e-10) and np.isclose(y[0] + y[2], 0.0, atol=1e-12)
print(f"Example 4.6 (Atwood): a = {a1:.4f} m/s², T = {T_at:.4f} N; y1(2 s) = {y[0]:.4f} m = ½ a t²")

# %% [markdown]
# ## Part 8: The exploding shell (Example 4.1) and the rocket (§3.9, Example 4.6, B6)

# %%
v0_s, ang_s = 50.0, np.radians(60)
vx_s, vy_s = v0_s * np.cos(ang_s), v0_s * np.sin(ang_s)
t_apex = vy_s / g; x_apex = vx_s * t_apex; h_apex = vy_s**2 / (2 * g); R = v0_s**2 * np.sin(2 * ang_s) / g
assert np.isclose(vx_s, 25.0) and np.isclose(vy_s, 43.301, atol=1e-3) and np.isclose(t_apex, 4.4140, atol=1e-4)
assert np.isclose(h_apex, 95.57, atol=1e-2) and np.isclose(R, 220.70, atol=1e-2) and np.isclose(x_apex, 110.35, atol=1e-2)
# Fragments after the explosion (equal masses; fragment 1 at rest): momentum gives v2 = (50, 0).
v2 = 2 * np.array([vx_s, 0.0]) - np.array([0.0, 0.0])
def land(x0, y0, vx, vy):
    s = solve_ivp(lambda t, y: [y[2], y[3], 0, -g], [0, 30], [x0, y0, vx, vy], events=hit, rtol=1e-11, atol=1e-13)
    return s.y[0, -1], s.t[-1]
x1_land, t1 = land(x_apex, h_apex, 0.0, 0.0); x2_land, t2 = land(x_apex, h_apex, v2[0], v2[1])
assert np.isclose(x1_land, 110.35, atol=1e-2) and np.isclose(x2_land, 331.05, atol=1e-2) and np.isclose(t1, t2)
assert np.isclose(0.5 * (x1_land + x2_land), R, rtol=1e-8) and np.isclose(x2_land, 1.5 * R, rtol=1e-8)
# The CM of the fragments follows the unexploded parabola at every instant after the explosion.
tt = np.linspace(0, t1, 50)
cm_x = 0.5 * (x_apex + (x_apex + v2[0] * tt)); cm_y = h_apex - 0.5 * g * tt**2
assert np.allclose(cm_x, vx_s * (t_apex + tt)) and np.allclose(cm_y, vy_s * (t_apex + tt) - 0.5 * g * (t_apex + tt) ** 2)
print(f"Example 4.1: fragments land at {x1_land:.2f} m and {x2_land:.2f} m; CM lands at {0.5*(x1_land+x2_land):.2f} m = R")

# Rocket: integrate m dv/dt = -u dm/dt numerically and compare with u ln(m0/m), eq. (3.11).
m0_r, mf_r, u_r, t_b = 1000.0, 300.0, 2500.0, 60.0
beta = (m0_r - mf_r) / t_b
sr = solve_ivp(lambda t, y: [u_r * beta / (m0_r - beta * t)], [0, t_b], [0.0], rtol=1e-12, atol=1e-12)
dv_closed = u_r * np.log(m0_r / mf_r)
assert np.isclose(sr.y[0, -1], dv_closed, rtol=1e-9) and np.isclose(dv_closed, 3009.9, atol=0.1)
sg = solve_ivp(lambda t, y: [u_r * beta / (m0_r - beta * t) - g, y[0]], [0, t_b], [0.0, 0.0], rtol=1e-12, atol=1e-12)
h_B6 = u_r * t_b - (u_r * mf_r / beta) * np.log(m0_r / mf_r) - 0.5 * g * t_b**2
assert np.isclose(sg.y[0, -1], dv_closed - g * t_b, rtol=1e-9) and np.isclose(sg.y[0, -1], 2421.3, atol=0.1)
assert np.isclose(sg.y[1, -1], h_B6, rtol=1e-9) and np.isclose(h_B6, 54944, atol=1)
assert u_r * beta > m0_r * g and np.isclose(u_r * beta, 29167, atol=1)
assert np.isclose(g * t_b, 588.6) and np.isclose(beta, 11.67, atol=5e-3) and np.isclose(np.log(m0_r / mf_r), 1.2040, atol=1e-4)
assert np.isclose(u_r * mf_r / beta, 64285.7, atol=0.1) and np.isclose(0.5 * g * t_b**2, 17658, atol=0.5)       # B6
assert np.isclose(u_r * mf_r / beta * np.log(m0_r / mf_r), 77398.3, atol=0.1)
assert np.isclose(u_r * (m0_r - mf_r) / 300.0, 5833, atol=1) and u_r * (m0_r - mf_r) / 300.0 < m0_r * g   # 5-minute burn
print(f"Example 4.6 (rocket): Δv = {dv_closed:.1f} m/s in free space, {sg.y[0, -1]:.1f} m/s after a 60 s vertical burn; "
      f"burnout height {sg.y[1, -1]:.0f} m (B6); thrust {u_r*beta:.0f} N > weight {m0_r*g:.0f} N")
print("Parts 7 and 8 passed.")

# %% [markdown]
# ## Part 9: SymPy checks of the derivations
#
# (3.17) Euler energy factor; C2 symplectic determinant and modified energy;
# the C1 series from `elliptic_k`; the resonance peak ω² = ω0² - γ²/2;
# elimination of λ in the bead problem, eq. (3.14); curl tests for A2.

# %%
x, v, h, w = sp.symbols("x v h omega", positive=True)
E_of = lambda X, V: V**2 / 2 + w**2 * X**2 / 2
assert sp.simplify(E_of(x + h * v, v - h * w**2 * x) / E_of(x, v)) == 1 + h**2 * w**2          # (3.17)
v1 = v - h * w**2 * x; x1 = x + h * v1
J = sp.Matrix([[sp.diff(x1, x), sp.diff(x1, v)], [sp.diff(v1, x), sp.diff(v1, v)]])
assert sp.simplify(J.det()) == 1                                                                # C2(b)
E_mod = lambda X, V: V**2 / 2 + w**2 * X**2 / 2 - h * w**2 * X * V / 2
assert sp.expand(E_mod(x1, v1) - E_mod(x, v)) == 0                                             # C2(c)

mm, th = sp.symbols("m theta", positive=True)
K_ser = sp.series(sp.elliptic_k(mm), mm, 0, 3).removeO()
ratio = sp.series((2 / sp.pi) * K_ser.subs(mm, sp.sin(th / 2) ** 2), th, 0, 6).removeO()
assert sp.expand(ratio - (1 + th**2 / 16 + 11 * th**4 / 3072)) == 0                           # C1
print("SymPy: T/T0 =", ratio)

gm, W0, Om = sp.symbols("gamma omega_0 Omega", positive=True)
D = (W0**2 - Om**2) ** 2 + gm**2 * Om**2
peak = sp.solve(sp.diff(D, Om), Om)
assert len(peak) == 1 and sp.simplify(peak[0] ** 2 - (W0**2 - gm**2 / 2)) == 0                 # §3.8 resonance peak

kk, M, G, xd, xdd, lam = sp.symbols("k m g xdot xddot lambda", real=True)
eqs = [sp.Eq(M * xdd, -2 * kk * x * lam), sp.Eq(M * (2 * kk * xd**2 + 2 * kk * x * xdd), lam - M * G)]
sol_bead = sp.solve(eqs, [xdd, lam], dict=True)[0]
assert sp.simplify(sol_bead[xdd] - (-2 * kk * x * (G + 2 * kk * xd**2) / (1 + 4 * kk**2 * x**2))) == 0     # (3.14)
assert sp.simplify(sol_bead[lam] * sp.sqrt(1 + 4 * kk**2 * x**2) - M * (G + 2 * kk * (xd**2 + x * sol_bead[xdd])) * sp.sqrt(1 + 4 * kk**2 * x**2)) == 0

X, Y, Z = sp.symbols("x y z", real=True)
curl = lambda F: sp.Matrix([sp.diff(F[2], Y) - sp.diff(F[1], Z), sp.diff(F[0], Z) - sp.diff(F[2], X), sp.diff(F[1], X) - sp.diff(F[0], Y)])
assert curl([Y, X, 0]) == sp.zeros(3, 1) and curl([-Y, X, 0]) == sp.Matrix([0, 0, 2])            # A2(a), (b)
assert curl([2 * X * Y, X**2 + Z**2, 2 * Y * Z]) == sp.zeros(3, 1)                              # A2(c)
Vc = -(X**2 * Y + Y * Z**2)
assert [-sp.diff(Vc, s) for s in (X, Y, Z)] == [2 * X * Y, X**2 + Z**2, 2 * Y * Z]
Fv = [-Y / (X**2 + Y**2), X / (X**2 + Y**2), sp.Integer(0)]
assert sp.simplify(curl(Fv)) == sp.zeros(3, 1)                                                  # B2(a): curl-free...
tt_ = sp.symbols("t", real=True)
loop = sp.integrate(sp.simplify(sum(f.subs({X: sp.cos(tt_), Y: sp.sin(tt_)}) * d for f, d in zip(Fv, [-sp.sin(tt_), sp.cos(tt_), 0]))), (tt_, 0, 2 * sp.pi))
assert loop == 2 * sp.pi                                                                        # ...but ∮ = 2π
I_quartic = quad(lambda u: 1 / np.sqrt(1 - u**4), 0, 1)[0]
assert np.isclose(I_quartic, Gamma(0.25) * Gamma(0.5) / (4 * Gamma(0.75))) and np.isclose(I_quartic, 1.31103, atol=1e-5)  # B4
print(f"B4: ∫₀¹ du/√(1-u⁴) = {I_quartic:.5f}.  Part 9 passed: all symbolic identities hold.")

# %% [markdown]
# ## Part 10: Numbers quoted in solutions.md (A1, A5, B3, B4, B7)
#
# B7 is checked against its exact piecewise solution (each half-cycle is SHM about
# x_c = ±μmg/k lasting π/ω0), not simulated: simulating it is "Your turn" item 1.

# %%
ms_ = np.array([2.0, 1.0, 3.0]); rs_ = np.array([[1, 0, 0], [0, 2, 0], [0, 0, -1.0]]); vv_ = np.array([[0, 3, 0], [-1, 0, 0], [2, 1, 0.0]])
P_, R_, L_ = ms_ @ vv_, ms_ @ rs_ / ms_.sum(), np.cross(rs_, ms_[:, None] * vv_).sum(0)
assert np.allclose(P_, [5, 9, 0]) and np.allclose(R_, [1 / 3, 1 / 3, -0.5]) and np.allclose(L_, [3, -6, 8])   # A1
a5, _, T5 = atwood_unknowns(0.8, 1.2)                                                                           # A5
assert np.isclose(a5, 1.962) and np.isclose(T5, 9.4176) and np.isclose(2 * T5, 18.835, atol=1e-3)
assert np.isclose(2.0 * g - 2 * T5, 0.785, atol=1e-3) and np.isclose((0.8 - 1.2) * a5 / 2.0, -0.392, atol=1e-3)
Vdw = lambda x: -x**2 / 2 + x**4 / 4                                                                            # B3(d)
assert np.isclose(Vdw(2.0), 2.0) and np.isclose(np.sqrt(2 * (Vdw(2.0) - Vdw(0.0))), 2.0) and np.isclose(Vdw(-2.0), Vdw(2.0))
T_q = lambda A: 4 / A * np.sqrt(1 / 2) * I_quartic                                                              # B4, m = α = 1
assert np.allclose([T_q(1), T_q(2), T_q(4)], [3.708, 1.854, 0.927], atol=1e-3)
assert np.allclose([Gamma(0.25), Gamma(0.5), Gamma(0.75)], [3.6256, 1.7725, 1.2254], atol=1e-4)
dA = 2 * 0.1 * g / 100.0; A_n = 0.5 - dA * np.arange(30); n_stick = int(np.argmax(100.0 * A_n <= 0.1 * g))     # B7
assert np.isclose(dA, 0.01962) and np.isclose(2 * dA, 0.03924) and np.allclose(A_n[1:4], [0.48038, 0.46076, 0.44114])
assert n_stick == 25 and np.isclose(A_n[24], 0.0291, atol=1e-4) and np.isclose(A_n[25], 0.0095) and np.isclose(25 * np.pi / 10, 7.854, atol=1e-3)
print(f"Part 10 passed: B7 sticks after {n_stick} half-cycles at |x| = {A_n[n_stick]:.4f} m, t = {n_stick * np.pi / 10:.3f} s")

# %% [markdown]
# ## Your turn
#
# 1. **Coulomb friction** (B7). Add `a_coulomb(x, v) = (-k x - μ m g sgn(v)) / m`
#    with m = 1, k = 100, μ = 0.1, x0 = 0.5, and integrate with a fixed-step RK4
#    (`rk4_step` takes a(x) only; write a two-argument variant). Record |x| at
#    each turning point (sign change of v) and stop when k|x| ≤ μ m g there.
#    Assert the amplitude drops by 4μmg/k = 0.03924 per period. How many
#    half-cycles before it sticks, and where? (Solutions B7 has the numbers.)
# 2. **Quartic well** (B4). For V = α x⁴ with α = m = 1, measure the period at
#    A = 0.5, 1, 2, 4 with `solve_ivp` events (as in `pendulum_period_numeric`)
#    and check T = (4/A) sqrt(m/2α) × 1.31103. Doubling A should halve T.
# 3. **Rocket with gravity** (B6). Integrate m v' = u β - m g with m = m0 - βt
#    using your own `rk4_vec`, compare with u ln(m0/m) - g t at every step, and
#    find the burn rate β below which the rocket never leaves the pad.
# 4. **Symplectic vs RK4 at equal cost.** RK4 evaluates a(x) four times per step;
#    give symplectic Euler a step four times smaller and rerun Part 1. Which
#    wins at 100 periods? At 10,000? (Reduce `theta0` to keep the run under a
#    minute, or vectorize over several initial conditions at once.)
# 5. **Phase lag by hand.** In Part 5, replace the projection fit with a direct
#    read-off: time of the first maximum of x minus time of the nearest earlier
#    maximum of the drive, times ω. Check it agrees with (3.10) at ω = 8, 10, 12.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
