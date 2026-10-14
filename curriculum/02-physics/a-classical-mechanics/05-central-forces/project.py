# %% [markdown]
# # CM-5 Project: Central Forces and the Two-Body Problem
#
# This script integrates the two-body problem *without* using the reduction to
# one body, then checks that everything the notes derive about the reduced
# problem is true of the raw numerical motion: the centre of mass moves in a
# straight line, the relative orbit is a conic with the eccentricity and
# semi-major axis predicted by $E$ and $\ell$, Kepler's third law has slope
# 3/2, circular orbits of power-law forces are stable iff $n > -3$, a $1/r^3$
# force makes the ellipse precess by exactly the predicted angle, the
# Laplace–Runge–Lenz vector points at perihelion with length $\mu k e$, a Monte
# Carlo of scattered alpha particles reproduces Rutherford's $\sin^{-4}(\theta/2)$,
# and a Hohmann transfer with two impulsive burns lands on a circular orbit.
# Every number quoted in the worked examples of the notes is printed here.
#
# Units: Parts 1–7 use "Kepler units" $\mu = 1$, $k = 1$ unless stated. Parts 8
# and 9 use SI. Run headless: `MPLBACKEND=Agg python3 project.py`.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T_START = time.perf_counter()

# %% [markdown]
# ## Part 0: Tools shared by every part
#
# The relative problem is a single particle of mass $\mu$ acted on by a central
# force $\mathbf F = F(r)\,\hat{\mathbf r}$ (eq. (3.4) of the notes). We integrate
# $\ddot{\mathbf r} = F(r)\,\mathbf r/(\mu r)$ in Cartesian coordinates with
# `solve_ivp` (DOP853, tight tolerances) and detect apsides with events on
# $\mathbf r\cdot\dot{\mathbf r} = 0$: $r\dot r$ passes from negative to positive at
# perihelion and from positive to negative at aphelion.


# %%
def kepler_accel(k, mu):
    """Acceleration of the reduced particle for V = -k/r: a = -(k/mu) r_vec / r^3."""
    def accel(x, y):
        r3 = (x * x + y * y) ** 1.5
        return -(k / mu) * x / r3, -(k / mu) * y / r3
    return accel


def power_law_accel(k, mu, n):
    """Acceleration for the attractive power law F(r) = -k r^n (notes §3.9)."""
    def accel(x, y):
        r = np.hypot(x, y)
        f_over_r = -(k / mu) * r ** (n - 1)
        return f_over_r * x, f_over_r * y
    return accel


def perihelion_event(t, y):
    """r . v, increasing through zero at a minimum of r."""
    return y[0] * y[2] + y[1] * y[3]


perihelion_event.direction = +1


def aphelion_event(t, y):
    return y[0] * y[2] + y[1] * y[3]


aphelion_event.direction = -1


def integrate_orbit(r0, v0, accel, t_end, events=None, n_out=2000, rtol=1e-11, atol=1e-13):
    """Integrate x'' = accel(x, y) from (r0, v0) to t_end. Returns the solve_ivp result."""
    def rhs(t, y):
        ax, ay = accel(y[0], y[1])
        return [y[2], y[3], ax, ay]
    t_eval = np.linspace(0.0, t_end, n_out)
    return solve_ivp(rhs, (0.0, t_end), [r0[0], r0[1], v0[0], v0[1]], method="DOP853",
                     rtol=rtol, atol=atol, t_eval=t_eval, events=events)


def energy_and_ell(x, y, vx, vy, mu, k):
    """E = mu v^2/2 - k/r and ell = mu (x vy - y vx), eqs. (3.8), (3.5)."""
    r = np.hypot(x, y)
    E = 0.5 * mu * (vx**2 + vy**2) - k / r
    ell = mu * (x * vy - y * vx)
    return E, ell


def ecc_from_E_ell(E, ell, mu, k):
    """Energy–eccentricity relation (3.14): e = sqrt(1 + 2 E ell^2 / (mu k^2)).

    For an exactly circular orbit the radicand is 0 up to roundoff; clip it.
    """
    return np.sqrt(np.maximum(0.0, 1.0 + 2.0 * E * ell**2 / (mu * k**2)))


def fit_conic(x, y):
    """Least-squares fit of 1/r = A + B cos(phi) + C sin(phi) to orbit points.

    Returns (c, e, phi_peri): r = c / (1 + e cos(phi - phi_peri)), eq. (3.12).
    """
    r = np.hypot(x, y)
    phi = np.arctan2(y, x)
    M = np.c_[np.ones_like(phi), np.cos(phi), np.sin(phi)]
    (A, B, C), *_ = np.linalg.lstsq(M, 1.0 / r, rcond=None)
    return 1.0 / A, np.hypot(B, C) / A, np.arctan2(C, B)


def unwrapped_angle_at_events(sol_, t_ev, y_ev):
    """Polar angle at the event points, continued through full turns.

    The wrapped angle of successive perihelia of a *closed* orbit is the same
    number, so the winding count has to come from the densely sampled path:
    phi_event = wrapped angle + 2 pi * (turns completed by that time).
    """
    phi_dense = np.unwrap(np.arctan2(sol_.y[1], sol_.y[0]))
    phi_wrapped = np.arctan2(y_ev[:, 1], y_ev[:, 0])
    phi_near = np.interp(t_ev, sol_.t, phi_dense)
    turns = np.round((phi_near - phi_wrapped) / (2 * np.pi))
    return phi_wrapped + 2 * np.pi * turns


def kepler_period(a, mu, k):
    """Kepler's third law (3.16): T = 2 pi sqrt(mu a^3 / k)."""
    return 2.0 * np.pi * np.sqrt(mu * a**3 / k)


# %% [markdown]
# ## Part 1: Two bodies, integrated directly, are one body in disguise
#
# Two masses $m_1 = 1$, $m_2 = 0.3$ with $G = 1$ (so $k = Gm_1m_2 = 0.3$,
# $\mu = m_1m_2/(m_1+m_2) = 3/13$). We integrate the eight-dimensional Newtonian
# system, then check §3.1–3.2: $\mathbf R(t)$ is a straight line and
# $\mathbf r = \mathbf r_1 - \mathbf r_2$ traces an ellipse whose $e$ and $a$ are
# given by (3.14) and (3.15) in terms of the relative $E$ and $\ell$.
#
# **Build it yourself first:** write `two_body_rhs` and the conic fit before
# reading on.

# %%
m1, m2, G = 1.0, 0.3, 1.0
M_tot, mu2, k2 = m1 + m2, m1 * m2 / (m1 + m2), G * m1 * m2


def two_body_rhs(t, y):
    """Newton's equations for two bodies: y = (x1, y1, x2, y2, vx1, vy1, vx2, vy2)."""
    x1, y1, x2, y2, vx1, vy1, vx2, vy2 = y
    dx, dy = x1 - x2, y1 - y2
    r3 = (dx * dx + dy * dy) ** 1.5
    # force on 1 is -G m1 m2 (r1 - r2)/r^3
    ax1, ay1 = -G * m2 * dx / r3, -G * m2 * dy / r3
    ax2, ay2 = +G * m1 * dx / r3, +G * m1 * dy / r3
    return [vx1, vy1, vx2, vy2, ax1, ay1, ax2, ay2]


R0, V0 = np.array([0.0, 0.0]), np.array([0.10, 0.05])       # CM position and velocity
rel0, vrel0 = np.array([1.0, 0.0]), np.array([0.0, 1.0])     # relative position and velocity
r1_0, r2_0 = R0 + (m2 / M_tot) * rel0, R0 - (m1 / M_tot) * rel0
v1_0, v2_0 = V0 + (m2 / M_tot) * vrel0, V0 - (m1 / M_tot) * vrel0

E_rel, ell_rel = energy_and_ell(*rel0, *vrel0, mu2, k2)
e_pred = ecc_from_E_ell(E_rel, ell_rel, mu2, k2)
a_pred = -k2 / (2 * E_rel)
T_pred = kepler_period(a_pred, mu2, k2)
assert E_rel < 0, "chosen initial conditions must be bound"

t_eval = np.sort(np.append(np.linspace(0, 3 * T_pred, 3000), T_pred))   # include t = T exactly
sol = solve_ivp(two_body_rhs, (0, 3 * T_pred), [*r1_0, *r2_0, *v1_0, *v2_0],
                method="DOP853", rtol=1e-11, atol=1e-13, t_eval=t_eval)
x1, y1, x2, y2 = sol.y[:4]
Rx = (m1 * x1 + m2 * x2) / M_tot
Ry = (m1 * y1 + m2 * y2) / M_tot
# (3.4): the centre of mass moves uniformly.
assert np.max(np.abs(Rx - (R0[0] + V0[0] * sol.t))) < 1e-9
assert np.max(np.abs(Ry - (R0[1] + V0[1] * sol.t))) < 1e-9

xr, yr = x1 - x2, y1 - y2
c_fit, e_fit, phi_peri_fit = fit_conic(xr, yr)
a_fit = c_fit / (1 - e_fit**2)
assert abs(e_fit - e_pred) < 1e-7, (e_fit, e_pred)
assert abs(a_fit - a_pred) < 1e-7, (a_fit, a_pred)
# c = ell^2/(mu k), eq. (3.12)
assert abs(c_fit - ell_rel**2 / (mu2 * k2)) < 1e-7
# ell is conserved along the raw two-body motion
ell_t = mu2 * (xr * (sol.y[5] - sol.y[7]) - yr * (sol.y[4] - sol.y[6]))
assert np.max(np.abs(ell_t - ell_rel)) < 1e-9
# the orbit closes after one period
i_T = np.argmin(np.abs(sol.t - T_pred))
assert sol.t[i_T] == T_pred
assert np.hypot(xr[i_T] - xr[0], yr[i_T] - yr[0]) < 1e-7
print(f"Part 1: m1={m1}, m2={m2}: mu={mu2:.4f}, k={k2}, E={E_rel:.4f}, ell={ell_rel:.4f}")
print(f"        e: formula {e_pred:.6f}, fitted {e_fit:.6f};  a: formula {a_pred:.6f}, fitted {a_fit:.6f}")
print(f"        period T = {T_pred:.4f}; CM straight-line residual < 1e-9. ")

fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 5))
axL.plot(x1, y1, lw=0.8, label="$m_1$")
axL.plot(x2, y2, lw=0.8, label="$m_2$")
axL.plot(Rx, Ry, "k--", lw=1, label="centre of mass")
axL.set_aspect("equal"); axL.legend(); axL.grid(alpha=.3)
axL.set_title("Lab frame: two bodies, CM drifts uniformly")
axR.plot(xr, yr, lw=1, label="relative orbit $\\mathbf{r}_1-\\mathbf{r}_2$")
phi_grid = np.linspace(0, 2 * np.pi, 400)
r_grid = c_fit / (1 + e_fit * np.cos(phi_grid - phi_peri_fit))
axR.plot(r_grid * np.cos(phi_grid), r_grid * np.sin(phi_grid), "r:", lw=2, label="fitted conic")
axR.plot(0, 0, "ko", label="focus = position of $m_2$")
axR.set_aspect("equal"); axR.legend(fontsize=8); axR.grid(alpha=.3)
axR.set_title(f"Relative orbit: e = {e_fit:.4f}, a = {a_fit:.4f}")
fig.tight_layout(); fig.savefig(FIG / "two_body_direct.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 2: Kepler's third law across ten orbits
#
# In Kepler units ($\mu = k = 1$) we launch ten orbits from perihelion with
# random $a \in [0.3, 5]$ and random $e$, measure the period as the time between
# successive perihelion events, and fit $\log T$ against $\log a$. Theorem 3.7
# says the slope is exactly $3/2$ and the intercept is $\log 2\pi$.

# %%
mu, k = 1.0, 1.0
a_list = np.logspace(np.log10(0.3), np.log10(5.0), 10)
e_list = rng.uniform(0.0, 0.8, size=10)
T_meas = []
for a, e in zip(a_list, e_list):
    rp = a * (1 - e)
    vp = np.sqrt((k / mu) * (2 / rp - 1 / a))           # vis-viva (3.18) at perihelion
    T_th = kepler_period(a, mu, k)
    s = integrate_orbit((rp, 0), (0, vp), kepler_accel(k, mu), 2.6 * T_th,
                        events=perihelion_event, n_out=50)
    tp = s.t_events[0]
    tp = tp[tp > 1e-6]                                   # drop the t = 0 apsis if reported
    T_meas.append(tp[1] - tp[0])
    assert abs((tp[1] - tp[0]) - T_th) / T_th < 1e-7, (a, e, tp[1] - tp[0], T_th)
T_meas = np.array(T_meas)
slope, intercept = np.polyfit(np.log(a_list), np.log(T_meas), 1)
assert abs(slope - 1.5) < 1e-6, slope
assert abs(np.exp(intercept) - 2 * np.pi) < 1e-5
print(f"Part 2: Kepler III log-log slope = {slope:.8f} (theory 1.5); intercept exp = {np.exp(intercept):.6f} (2pi = {2*np.pi:.6f})")

fig, ax = plt.subplots(figsize=(5.5, 4.5))
ax.loglog(a_list, T_meas, "o", label="measured period")
ax.loglog(a_list, kepler_period(a_list, mu, k), "-", label="$T = 2\\pi\\sqrt{\\mu a^3/k}$")
ax.set_xlabel("semi-major axis $a$"); ax.set_ylabel("period $T$"); ax.grid(alpha=.3, which="both"); ax.legend()
ax.set_title(f"Kepler's third law: slope {slope:.4f}")
fig.tight_layout(); fig.savefig(FIG / "kepler_third_law.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 3: The effective potential and the four conics
#
# $V_{\rm eff}(r) = -k/r + \ell^2/(2\mu r^2)$ with $\ell = \mu = k = 1$, so the
# minimum is at $r_0 = \ell^2/(\mu k) = 1$ with $V_{\rm eff}(r_0) = -\mu k^2/(2\ell^2) = -1/2$.
# Four energies give the circle, an ellipse, the parabola and a hyperbola
# (§3.4). Turning points solve $E = V_{\rm eff}(r)$, a quadratic in $r$.

# %%
ell = 1.0
V_eff = lambda r: -k / r + ell**2 / (2 * mu * r**2)
energies = [-0.5, -0.3, 0.0, 0.3]
names = ["circle", "ellipse", "parabola", "hyperbola"]
fig = plt.figure(figsize=(13, 6.5))
axV = fig.add_subplot(1, 2, 1)
rr = np.linspace(0.3, 6, 600)
axV.plot(rr, V_eff(rr), "k", lw=2, label="$V_{\\rm eff}$")
axV.plot(rr, -k / rr, "k:", lw=1, label="$-k/r$")
axV.plot(rr, ell**2 / (2 * mu * rr**2), "k--", lw=1, label="$\\ell^2/2\\mu r^2$")
for E, name, col in zip(energies, names, ["C0", "C1", "C2", "C3"]):
    e = ecc_from_E_ell(E, ell, mu, k)
    c = ell**2 / (mu * k)
    rp = c / (1 + e)
    axV.axhline(E, color=col, lw=1, alpha=.8, label=f"E = {E:+.1f} ({name}, e = {e:.3f})")
    # turning points: E r^2 + k r - ell^2/(2 mu) = 0
    if abs(E) > 1e-12:
        roots = np.roots([E, k, -ell**2 / (2 * mu)])
        roots = np.sort(roots[roots > 0].real)
    else:
        roots = np.array([ell**2 / (2 * mu * k)])
    assert abs(roots[0] - rp) < 1e-12, (roots, rp)
    if E < 0:
        assert abs(roots[-1] - c / (1 - e)) < 1e-12
    axV.plot(roots, [E] * len(roots), "o", color=col)
axV.set_ylim(-0.7, 0.6); axV.set_xlim(0.3, 6); axV.grid(alpha=.3); axV.legend(fontsize=8)
axV.set_xlabel("r"); axV.set_ylabel("energy"); axV.set_title("Effective potential, energies, turning points")

for j, (E, name, col) in enumerate(zip(energies, names, ["C0", "C1", "C2", "C3"])):
    e = ecc_from_E_ell(E, ell, mu, k)
    c = ell**2 / (mu * k)
    rp = c / (1 + e)
    vp = ell / (mu * rp)
    if E < 0:
        t_end = kepler_period(c / (1 - e**2), mu, k)
    else:
        t_end = 14.0
    s = integrate_orbit((rp, 0), (0, vp), kepler_accel(k, mu), t_end, n_out=1500)
    ax = fig.add_subplot(2, 4, 3 + j + (4 if j >= 2 else 0) - (2 if j >= 2 else 0))
    ax.plot(s.y[0], s.y[1], color=col, lw=1.2)
    ax.plot(0, 0, "k*", ms=9)
    ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_title(f"{name}: E={E:+.1f}, e={e:.3f}", fontsize=10)
    ax.set_xlim(-4.5, 2.5); ax.set_ylim(-3.5, 3.5)
    # Energy is conserved along the numerical orbit
    Et, _ = energy_and_ell(s.y[0], s.y[1], s.y[2], s.y[3], mu, k)
    assert np.max(np.abs(Et - E)) < 1e-8
fig.tight_layout(); fig.savefig(FIG / "effective_potential_orbits.png", dpi=120); plt.close(fig)
# Example 4.2: the E = -0.3, ell = 1 orbit in Kepler units, constants fixed by E and ell
E42 = -0.3
e42 = ecc_from_E_ell(E42, ell, mu, k); c42 = ell**2 / (mu * k)
a42 = -k / (2 * E42); b42 = a42 * np.sqrt(1 - e42**2)
assert abs(a42 - c42 / (1 - e42**2)) < 1e-12
# Launch at r = 1, phi = 0 with transverse speed ell/(mu r) = 1 and outward radial speed from E:
vr42 = np.sqrt(2 * (E42 + k) / mu - ell**2 / mu**2)
s42 = integrate_orbit((1.0, 0.0), (vr42, ell / mu), kepler_accel(k, mu), kepler_period(a42, mu, k), n_out=800)
c_f, e_f, phi_p42 = fit_conic(s42.y[0], s42.y[1])
assert abs(c_f - c42) < 1e-8 and abs(e_f - e42) < 1e-8 and abs(phi_p42 + np.pi / 2) < 1e-8
print(f"Part 3, Ex 4.2: launch r=1, phi=0, v_r = {vr42:.4f}, v_phi = 1 -> fitted c = {c_f:.6f}, e = {e_f:.6f}, "
      f"perihelion at phi = {np.degrees(phi_p42):.4f} deg, so r = 1/(1 - {e42:.4f} sin phi)")
print(f"Part 3, Ex 4.2: E={E42}, ell=1: e = {e42:.4f}, c = {c42:.4f}, r_p = {c42/(1+e42):.4f}, r_a = {c42/(1-e42):.4f}, "
      f"a = {a42:.4f}, b = {b42:.4f}, T = {kepler_period(a42, mu, k):.3f}, v_p = {ell/(mu*c42/(1+e42)):.4f}, v_a = {ell/(mu*c42/(1-e42)):.4f}")
print("Part 3: turning points from the quadratic match c/(1±e); energy conserved to 1e-8 on all four conics.")

# %% [markdown]
# ## Part 4: Power-law forces $F = -kr^n$: closed, precessing, unstable (Example 4.5)
#
# Theorem 3.9: a circular orbit is stable iff $n > -3$, and small radial
# oscillations have $\omega_r/\omega_\phi = \sqrt{n+3}$, so the angle from one
# perihelion to the next is $\Delta\phi = 2\pi/\sqrt{n+3}$. We launch from $r_0 = 1$
# with the circular speed multiplied by $1 + \epsilon$ and measure $\Delta\phi$.
# $n = -2$ and $n = 1$ close ($\Delta\phi = 2\pi$ and $\pi$); $n = -1$ and
# $n = -2.1$ precess ($\sqrt{2}$ and $\sqrt{0.9}$ are irrational); $n = -3.5$
# falls in or escapes depending on the sign of the perturbation.

# %%
def apsidal_data(n, eps, n_orbits=6):
    """Return (times, angles) of successive perihelia for F = -k r^n launched at
    r0 = 1 with speed (1+eps) * circular speed."""
    r0 = 1.0
    v_circ = np.sqrt(k * r0 ** (n + 1) / mu)
    omega_phi = v_circ / r0
    t_end = n_orbits * 2 * np.pi / omega_phi * 1.3
    s = integrate_orbit((r0, 0), (0, (1 + eps) * v_circ), power_law_accel(k, mu, n), t_end,
                        events=perihelion_event, n_out=3000)
    te = s.t_events[0]; ye = s.y_events[0]
    keep = te > 1e-6
    te, ye = te[keep], ye[keep]
    phi = unwrapped_angle_at_events(s, te, ye)
    return s, te, phi


fig, axes = plt.subplots(1, 5, figsize=(18, 4))
cases = [(-2, "n = -2: closed (ellipse, focus at origin)"),
         (-1, "n = -1: precesses, apsidal angle $\\pi/\\sqrt{2}$"),
         (1, "n = 1: closed (ellipse centred at origin)"),
         (-2.1, "n = -2.1: precesses slowly"),
         (-3.5, "n = -3.5: circular orbit unstable")]
ratios = {}
for ax, (n, title) in zip(axes, cases):
    if n > -3:
        s, te, phi = apsidal_data(n, 0.15)
        ax.plot(s.y[0], s.y[1], lw=0.8)
        # small-amplitude measurement for the frequency ratio
        s2, te2, phi2 = apsidal_data(n, 0.03)
        dphi = np.diff(phi2)                      # perihelion-to-perihelion angle
        ratios[n] = 2 * np.pi / dphi.mean()       # omega_r / omega_phi = 2 pi / Delta phi
    else:
        for eps, col, lab in [(-0.02, "C3", "slightly slow: falls in"), (+0.02, "C0", "slightly fast: escapes")]:
            v_circ = np.sqrt(k / mu)
            def stop(t, y):
                return (np.hypot(y[0], y[1]) - 0.03) * (12.0 - np.hypot(y[0], y[1]))
            stop.terminal = True
            s = integrate_orbit((1, 0), (0, (1 + eps) * v_circ), power_law_accel(k, mu, n), 60.0,
                                events=stop, n_out=4000, rtol=1e-10, atol=1e-12)
            ax.plot(s.y[0], s.y[1], color=col, lw=0.8, label=lab)
        ax.legend(fontsize=7, loc="upper right")
        ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
    ax.plot(0, 0, "k*", ms=8); ax.set_aspect("equal"); ax.grid(alpha=.3); ax.set_title(title, fontsize=9)
fig.tight_layout(); fig.savefig(FIG / "power_law_orbits.png", dpi=120); plt.close(fig)

# Theorem 3.9: omega_r/omega_phi = sqrt(n+3). Exact for n = -2 and n = 1 (closed orbits).
assert abs(ratios[-2] - 1.0) < 1e-6, ratios[-2]
assert abs(ratios[1] - 2.0) < 1e-6, ratios[1]
# For other n the small-amplitude formula holds to O(eps^2).
assert abs(ratios[-1] - np.sqrt(2.0)) < 3e-3, ratios[-1]
assert abs(ratios[-2.1] - np.sqrt(0.9)) < 3e-3, ratios[-2.1]
for n in (-2, -1, 1, -2.1):
    print(f"Part 4: n = {n:+.1f}: omega_r/omega_phi measured {ratios[n]:.6f}, theory sqrt(n+3) = {np.sqrt(n+3):.6f}; "
          f"apsidal angle {np.pi/ratios[n]*180/np.pi:.2f} deg")

# %% [markdown]
# ## Part 5: Precession under a $1/r^3$ force, and the GR number for Mercury
#
# For $F = -k/r^2 - \lambda/r^3$ Binet's equation stays linear:
# $u'' + \beta^2 u = \mu k/\ell^2$ with $\beta^2 = 1 - \mu\lambda/\ell^2$, so the
# perihelion advances by exactly $\Delta = 2\pi(1/\beta - 1)$ per orbit, which is
# $\pi\mu\lambda/\ell^2$ to first order (§3.10, eq. (3.22)). We measure it from successive
# perihelion events. Then the general-relativistic formula
# $\delta\phi = 6\pi GM/(c^2 a(1-e^2))$ is evaluated for Mercury.

# %%
lam = 0.02
a0, e0 = 1.0, 0.4
ell0 = np.sqrt(mu * k * a0 * (1 - e0**2))
rp0 = a0 * (1 - e0)
vp0 = ell0 / (mu * rp0)


def perturbed_accel(x, y):
    r = np.hypot(x, y)
    f_over_r = -(k / mu) / r**3 - (lam / mu) / r**4      # F(r)/(mu r) with F = -k/r^2 - lam/r^3
    return f_over_r * x, f_over_r * y


T0 = kepler_period(a0, mu, k)
s = integrate_orbit((rp0, 0), (0, vp0), perturbed_accel, 12.5 * T0, events=perihelion_event, n_out=6000)
te = s.t_events[0]; ye = s.y_events[0]
keep = te > 1e-6
phi_peri = unwrapped_angle_at_events(s, te[keep], ye[keep])
prec_meas = np.diff(phi_peri) - 2 * np.pi
beta = np.sqrt(1 - mu * lam / ell0**2)
prec_exact = 2 * np.pi * (1 / beta - 1)
prec_first = np.pi * mu * lam / ell0**2
s_prec, ye_prec = s, ye[keep]                         # kept for exercise B8 in Part 10
assert np.max(np.abs(prec_meas - prec_exact)) < 1e-6, (prec_meas, prec_exact)
assert abs(prec_first - prec_exact) / prec_exact < 0.03
print(f"Part 5: 1/r^3 perturbation lam={lam}, ell^2={ell0**2:.3f}: precession/orbit measured {prec_meas.mean():.6f} rad, "
      f"exact 2pi(1/beta-1) = {prec_exact:.6f}, first order pi*mu*lam/ell^2 = {prec_first:.6f} "
      f"({100*(prec_first-prec_exact)/prec_exact:+.2f}%)")

fig, ax = plt.subplots(figsize=(5.5, 5.5))
ax.plot(s.y[0], s.y[1], lw=0.5)
ax.plot(ye[keep, 0], ye[keep, 1], "ro", ms=4, label="successive perihelia")
ax.plot(0, 0, "k*", ms=9); ax.set_aspect("equal"); ax.grid(alpha=.3); ax.legend()
ax.set_title(f"$F=-k/r^2-\\lambda/r^3$: perihelion advances {np.degrees(prec_exact):.2f}° per orbit")
fig.tight_layout(); fig.savefig(FIG / "precession.png", dpi=120); plt.close(fig)

# GR perihelion shift of Mercury (formula (3.23), stated in §3.10, not derived here)
GM_sun = 1.32712440018e20           # m^3/s^2
c_light = 2.99792458e8
a_merc, e_merc, T_merc_days = 5.7909e10, 0.20563, 87.969
dphi_orbit = 6 * np.pi * GM_sun / (c_light**2 * a_merc * (1 - e_merc**2))
orbits_per_century = 36525.0 / T_merc_days
arcsec = dphi_orbit * orbits_per_century * (180 / np.pi) * 3600
assert 42.8 < arcsec < 43.2, arcsec
print(f"        GR: Mercury perihelion shift {dphi_orbit:.3e} rad/orbit = {arcsec:.2f} arcsec/century (observed anomalous shift about 43 arcsec/century)")

# %% [markdown]
# ## Part 6: Rutherford scattering by Monte Carlo
#
# Repulsive Coulomb force $F = +k/r^2$ in units $\mu = k = 1$, $E = 1/2$
# ($v_\infty = 1$), so $k/2E = 1$ and the deflection law (3.26) is
# $\cot(\theta/2) = b$. We launch $N$ particles with impact parameters drawn
# uniformly over a disk $b_{\min} \le b \le b_{\max}$ (density $\propto b$), integrate
# every trajectory (all at once, as one big ODE), histogram the deflection angles,
# convert to $d\sigma/d\Omega$ by (3.24), and compare with the Rutherford formula
# (3.27), which in these units is $1/(4\sin^4(\theta/2))$.

# %%
N_mc = 12000
b_min, b_max = 0.1, 6.0
b = np.sqrt(b_min**2 + rng.uniform(size=N_mc) * (b_max**2 - b_min**2))
R_start = 300.0
E_mc = 0.5
v_inf = np.sqrt(2 * E_mc / mu)
# Launch at x = -R with the speed that makes the total energy exactly E (the particle
# has already climbed part of the Coulomb hill) and with the transverse offset that
# makes ell = mu v_inf b exactly. Then E and ell, hence theta, are exact; the only
# error is that the velocity at R is not quite parallel to the asymptote, ~ b/(2R^2).
v_R = np.sqrt(2 * (E_mc - k / np.hypot(R_start, b)) / mu)
b_R = b * v_inf / v_R
Y0 = np.concatenate([-R_start * np.ones(N_mc), b_R, v_R, np.zeros(N_mc)])


def coulomb_rhs(t, Y):
    x, y, vx, vy = Y.reshape(4, -1)
    r3 = (x * x + y * y) ** 1.5
    return np.concatenate([vx, vy, (k / mu) * x / r3, (k / mu) * y / r3])


t_mc = time.perf_counter()
sol_mc = solve_ivp(coulomb_rhs, (0, 2 * R_start + 60), Y0, method="DOP853", rtol=1e-9, atol=1e-12)
xf, yf, vxf, vyf = sol_mc.y[:, -1].reshape(4, -1)
assert np.all(np.hypot(xf, yf) > 0.8 * R_start), "all particles must be far away again"
theta_mc = np.arccos(vxf / np.hypot(vxf, vyf))                       # deflection from the +x direction
theta_th = 2 * np.arctan(k / (2 * E_mc * b))                         # eq. (3.26)
err_theta = np.max(np.abs(theta_mc - theta_th))
assert err_theta < 5e-4, err_theta                                   # finite start distance ~ b/(2R^2)
# energy conservation per particle
E_final = 0.5 * mu * (vxf**2 + vyf**2) + k / np.hypot(xf, yf)
assert np.max(np.abs(E_final - E_mc)) < 1e-8

sigma_disk = np.pi * (b_max**2 - b_min**2)
edges = np.linspace(theta_th.min(), np.radians(160), 15)
counts, _ = np.histogram(theta_mc, bins=edges)
mid = 0.5 * (edges[1:] + edges[:-1]); dth = np.diff(edges)
dsdo_mc = (counts / N_mc) * sigma_disk / (2 * np.pi * np.sin(mid) * dth)
dsdo_err = (np.sqrt(np.maximum(counts, 1)) / N_mc) * sigma_disk / (2 * np.pi * np.sin(mid) * dth)
rutherford = (k / (4 * E_mc))**2 / np.sin(mid / 2)**4
good = counts >= 50
rel = np.abs(dsdo_mc[good] - rutherford[good]) / rutherford[good]
assert np.all(rel < 0.35), rel
assert np.median(rel) < 0.12, np.median(rel)
print(f"Part 6: {N_mc} trajectories integrated in {time.perf_counter()-t_mc:.1f} s; max |theta_num - theta(b)| = {err_theta:.1e} rad")
print(f"        d sigma/d Omega vs Rutherford on {good.sum()} well-populated bins: median rel. error {np.median(rel):.3f}, max {rel.max():.3f}")

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5))
axA.plot(b, np.degrees(theta_mc), ".", ms=2, label="integrated trajectories")
bb = np.linspace(b_min, b_max, 300)
axA.plot(bb, np.degrees(2 * np.arctan(k / (2 * E_mc * bb))), "r-", lw=1, label="$\\theta = 2\\,\\mathrm{arccot}(2Eb/k)$")
axA.set_xlabel("impact parameter $b$ (units of $k/2E$)"); axA.set_ylabel("deflection $\\theta$ (deg)"); axA.grid(alpha=.3); axA.legend()
axB.errorbar(np.degrees(mid), dsdo_mc, yerr=dsdo_err, fmt="o", ms=4, label="Monte Carlo")
tt = np.linspace(mid[0], mid[-1], 300)
axB.semilogy(np.degrees(tt), (k / (4 * E_mc))**2 / np.sin(tt / 2)**4, "r-", label="Rutherford $(k/4E)^2\\sin^{-4}(\\theta/2)$")
axB.set_xlabel("scattering angle $\\theta$ (deg)"); axB.set_ylabel("$d\\sigma/d\\Omega$ (units of $(k/2E)^2$)")
axB.grid(alpha=.3, which="both"); axB.legend()
fig.tight_layout(); fig.savefig(FIG / "rutherford_mc.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: The Laplace–Runge–Lenz vector points at perihelion
#
# $\mathbf A = \mathbf p\times\mathbf L - \mu k\,\hat{\mathbf r}$ (§3.12, eq. (3.29)). Along the
# raw two-body orbit of Part 1 we check that $\mathbf A$ is constant, that
# $|\mathbf A| = \mu k e$, and that its direction is the direction of the
# fitted perihelion.

# %%
vxr, vyr = sol.y[4] - sol.y[6], sol.y[5] - sol.y[7]
r_rel = np.hypot(xr, yr)
Ax = mu2 * ell_rel * vyr - mu2 * k2 * xr / r_rel
Ay = -mu2 * ell_rel * vxr - mu2 * k2 * yr / r_rel
assert np.max(np.abs(Ax - Ax[0])) < 1e-8 and np.max(np.abs(Ay - Ay[0])) < 1e-8
A_mag = np.hypot(Ax[0], Ay[0])
assert abs(A_mag / (mu2 * k2) - e_fit) < 1e-7
ang_A = np.arctan2(Ay[0], Ax[0])
assert abs(np.angle(np.exp(1j * (ang_A - phi_peri_fit)))) < 1e-6
i_min = np.argmin(r_rel)
assert abs(np.angle(np.exp(1j * (ang_A - np.arctan2(yr[i_min], xr[i_min]))))) < 1e-2
# |A|^2 = mu^2 k^2 + 2 mu E ell^2  (reproduces the energy–eccentricity relation)
assert abs(A_mag**2 - (mu2**2 * k2**2 + 2 * mu2 * E_rel * ell_rel**2)) < 1e-10
print(f"Part 7: LRL vector constant to 1e-8; |A|/(mu k) = {A_mag/(mu2*k2):.6f} = e; direction {np.degrees(ang_A):.2f} deg = perihelion direction {np.degrees(phi_peri_fit):.2f} deg")
# |A| = mu k e across circle, ellipses, parabola, hyperbola (mu = k = 1, launched from perihelion r_p = 1)
for e_t in (0.0, 0.3, 0.6, 0.9, 1.0, 1.5):
    s_t = integrate_orbit((1.0, 0.0), (0.0, np.sqrt(1 + e_t)), kepler_accel(k, mu), 6.0, n_out=400)
    x_, y_, vx_, vy_ = s_t.y
    l_t = mu * (x_ * vy_ - y_ * vx_)
    A_t = np.hypot(mu * l_t * vy_ - mu * k * x_ / np.hypot(x_, y_), -mu * l_t * vx_ - mu * k * y_ / np.hypot(x_, y_))
    assert np.max(np.abs(A_t / (mu * k) - e_t)) < 1e-8, (e_t, A_t[:3])
print("        |A|/(mu k) = e to 1e-8 along orbits with e = 0, 0.3, 0.6, 0.9, 1 (parabola), 1.5 (hyperbola)")

# %% [markdown]
# ## Part 8: Hohmann transfer, LEO to GEO, with two impulsive burns (Example 4.4)
#
# SI units, Earth $GM = 3.986004418\times10^{14}$ m$^3$/s$^2$. Start on a circular
# orbit at $r_1 = 6771$ km (400 km altitude). Burn 1 raises the speed to the
# perigee speed of the transfer ellipse ($a_t = (r_1+r_2)/2$, vis-viva (3.18));
# coast half a period to apogee; burn 2 raises the speed to the circular speed
# at $r_2 = 42164$ km. We simulate it and check the final orbit is circular.

# %%
GM_E = 3.986004418e14
r1, r2 = 6.771e6, 4.2164e7
a_t = 0.5 * (r1 + r2)
v_c1, v_c2 = np.sqrt(GM_E / r1), np.sqrt(GM_E / r2)
v_p_t = np.sqrt(GM_E * (2 / r1 - 1 / a_t))
v_a_t = np.sqrt(GM_E * (2 / r2 - 1 / a_t))
dv1, dv2 = v_p_t - v_c1, v_c2 - v_a_t
t_transfer = np.pi * np.sqrt(a_t**3 / GM_E)
assert abs(v_a_t - v_p_t * r1 / r2) < 1e-9          # angular momentum on the transfer ellipse

acc_E = kepler_accel(GM_E, 1.0)                       # per unit mass: k/mu = GM
# leg 0: a quarter of the parking orbit
s0 = integrate_orbit((r1, 0), (0, v_c1), acc_E, 0.25 * kepler_period(r1, 1.0, GM_E), n_out=300)
x, y, vx, vy = s0.y[:, -1]
# burn 1: tangential, add dv1 along the velocity
sp = np.hypot(vx, vy); vx, vy = vx * (sp + dv1) / sp, vy * (sp + dv1) / sp
s1 = integrate_orbit((x, y), (vx, vy), acc_E, 1.05 * t_transfer, events=aphelion_event, n_out=800)
t_ap = s1.t_events[0][s1.t_events[0] > 1.0][0]
assert abs(t_ap - t_transfer) / t_transfer < 1e-8
x, y, vx, vy = s1.y_events[0][s1.t_events[0] > 1.0][0]
assert abs(np.hypot(x, y) - r2) / r2 < 1e-8
# burn 2
sp = np.hypot(vx, vy); vx, vy = vx * (sp + dv2) / sp, vy * (sp + dv2) / sp
s2 = integrate_orbit((x, y), (vx, vy), acc_E, kepler_period(r2, 1.0, GM_E), n_out=600)
r_final = np.hypot(s2.y[0], s2.y[1])
assert np.max(np.abs(r_final - r2)) / r2 < 1e-7, "final orbit is not circular at GEO"
E_f, l_f = energy_and_ell(s2.y[0], s2.y[1], s2.y[2], s2.y[3], 1.0, GM_E)
assert np.max(ecc_from_E_ell(E_f, l_f, 1.0, GM_E)) < 1e-6
print(f"Part 8: Hohmann LEO({r1/1e3:.0f} km) -> GEO({r2/1e3:.0f} km), a_t = {a_t/1e3:.1f} km: v_c1 = {v_c1:.1f} m/s, v_p(transfer) = {v_p_t:.1f} m/s, "
      f"dv1 = {dv1:.1f} m/s; v_a(transfer) = {v_a_t:.1f} m/s, v_c2 = {v_c2:.1f} m/s, dv2 = {dv2:.1f} m/s; "
      f"total {dv1+dv2:.1f} m/s; transfer time {t_transfer/3600:.3f} h")
fig, ax = plt.subplots(figsize=(6, 6))
for s_, col, lab in [(s0, "C0", "parking orbit (LEO)"), (s1, "C1", "transfer ellipse"), (s2, "C2", "final orbit (GEO)")]:
    mask = s_.t <= (t_transfer if s_ is s1 else np.inf)
    ax.plot(s_.y[0][mask] / 1e6, s_.y[1][mask] / 1e6, color=col, lw=1.2, label=lab)
earth = plt.Circle((0, 0), 6.371, color="0.7"); ax.add_patch(earth)
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_xlabel("x (Mm)"); ax.set_ylabel("y (Mm)")
ax.set_title(f"Hohmann transfer: Δv = {dv1:.0f} + {dv2:.0f} m/s, {t_transfer/3600:.2f} h")
fig.tight_layout(); fig.savefig(FIG / "hohmann.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 9: The numbers in the worked examples
#
# Example 4.1 (Earth–Sun), 4.3 (an interstellar visitor), 4.5 (a radial kick in
# low Earth orbit), 4.6 (alpha particles on gold), and the reference values for
# "Your turn" 3 (Halley's comet). Solar-system numbers use the measured products
# $GM$ (known to about ten digits), not $G\times M$ ($G$ is known to about 2e-5).

# %%
GM_sun, GM_earth = 1.32712440018e20, 3.986004418e14      # m^3/s^2
AU, yr = 1.495978707e11, 365.25 * 86400
G_SI, M_sun_kg = 6.67430e-11, 1.9885e30                   # only to show why G*M is the wrong input

# Example 4.1: Earth–Sun. mu/m_E = M/(M+m); k/mu = G(M+m).
ratio = GM_sun / (GM_sun + GM_earth)
a_E, e_E = AU, 0.0167
GMtot = GM_sun + GM_earth
v_mean = np.sqrt(GMtot / a_E)
v_peri = v_mean * np.sqrt((1 + e_E) / (1 - e_E))          # (3.17)
v_apo = v_mean * np.sqrt((1 - e_E) / (1 + e_E))
T_E, T_sun_only = kepler_period(a_E, 1.0, GMtot), kepler_period(a_E, 1.0, GM_sun)
T_GxM = kepler_period(a_E, 1.0, G_SI * M_sun_kg)
wobble = a_E * GM_earth / GMtot                           # Sun's distance from the barycentre
assert abs(v_peri / v_apo - (1 + e_E) / (1 - e_E)) < 1e-12
assert abs(v_peri**2 / (GMtot * (2 / (a_E * (1 - e_E)) - 1 / a_E)) - 1) < 1e-12   # vis-viva (3.18) agrees
print(f"Part 9, Ex 4.1 Earth-Sun: mu/m_E = {ratio:.9f} (1 - {1-ratio:.3e}); Sun-barycentre distance {wobble/1e3:.0f} km; "
      f"r_p = {a_E*(1-e_E)/1e9:.2f} Gm, r_a = {a_E*(1+e_E)/1e9:.2f} Gm")
print(f"        sqrt(G(M+m)/a) = {v_mean/1e3:.3f} km/s, factor sqrt((1+e)/(1-e)) = {np.sqrt((1+e_E)/(1-e_E)):.5f}, "
      f"v_peri = {v_peri/1e3:.3f} km/s, v_apo = {v_apo/1e3:.3f} km/s, v_p/v_a = {v_peri/v_apo:.4f}")
print(f"        T = {T_E/86400:.5f} d with M+m, {T_sun_only/86400:.5f} d with M alone (difference {T_sun_only-T_E:.1f} s); "
      f"with G*M: {T_GxM/86400:.5f} d (off by {(T_sun_only-T_GxM)/86400:.4f} d = {(T_sun_only-T_GxM)/(T_sun_only-T_E):.0f}x the mass effect). "
      f"Sidereal year: 365.25636 d")

# Example 4.3: an interstellar visitor (numbers close to 1I/'Oumuamua): v_inf = 26.4 km/s, q = 0.255 AU
v_inf_c, q_c = 26.4e3, 0.255 * AU
a_hyp = GM_sun / v_inf_c**2                    # |a| = k/2E, so v_inf^2 = (k/mu)/|a|
e_c = 1 + q_c / a_hyp                          # q = |a|(e - 1)
theta_c = 2 * np.arcsin(1 / e_c)               # (3.26): sin(theta/2) = 1/e
v_p_c = np.sqrt(GM_sun * (2 / q_c + 1 / a_hyp))   # vis-viva (3.18) with a = -|a|
b_c = a_hyp * np.sqrt(e_c**2 - 1)              # (3.26): cot(theta/2) = 2Eb/k = b/|a|
assert abs(b_c / (v_p_c * q_c / v_inf_c) - 1) < 1e-12      # ell = mu v_inf b = mu v_p q
assert abs(1 / np.tan(theta_c / 2) - b_c / a_hyp) < 1e-12
print(f"Part 9, Ex 4.3 visitor: |a| = {a_hyp:.3e} m = {a_hyp/AU:.3f} AU, e = {e_c:.3f}, sin(theta/2) = {1/e_c:.3f}, b = {b_c/AU:.3f} AU, "
      f"deflection = {np.degrees(theta_c):.1f} deg, v_peri = {v_p_c/1e3:.1f} km/s, E/mu = {0.5*v_inf_c**2/1e6:.1f} MJ/kg, "
      f"v_inf b = v_p q = {v_inf_c*b_c:.3e} m^2/s")

# Example 4.5: a 10 m/s radial kick to a circular orbit at r0 = 6771 km (n = -2, so omega_r = omega_phi)
r0_k, dv_r = 6.771e6, 10.0
v_c0 = np.sqrt(GM_E / r0_k)
omega0 = v_c0 / r0_k                           # omega_phi = omega_r, (3.21) with n = -2
amp_lin = dv_r / omega0                        # linearized radial oscillation amplitude
e_k = dv_r / v_c0                              # exact, from (3.14)
E_k = 0.5 * (v_c0**2 + dv_r**2) - GM_E / r0_k
assert abs(ecc_from_E_ell(E_k, r0_k * v_c0, 1.0, GM_E) - e_k) < 1e-12
rp_k, ra_k = r0_k / (1 + e_k), r0_k / (1 - e_k)   # c = ell^2/(mu k) = r0 because the kick leaves ell alone
s_k = integrate_orbit((r0_k, 0), (dv_r, v_c0), acc_E, kepler_period(r0_k, 1.0, GM_E), n_out=4000)
r_k = np.hypot(s_k.y[0], s_k.y[1])
assert abs(r_k.max() - ra_k) < 1.0 and abs(r_k.min() - rp_k) < 1.0       # metres
print(f"Part 9, Ex 4.5 radial kick: omega = {omega0:.4e} rad/s, period {2*np.pi/omega0/60:.1f} min; linear amplitude "
      f"{amp_lin/1e3:.2f} km; exact e = {e_k:.4e}, r_a - r0 = {(ra_k-r0_k)/1e3:.3f} km, r0 - r_p = {(r0_k-rp_k)/1e3:.3f} km")

# Example 4.6: alpha particles (Z=2) on gold (Z=79) at E = 5 MeV
e2_4pi_eps0 = 1.439964548          # MeV fm
Z1, Z2, E_alpha = 2, 79, 5.0
k_coul = Z1 * Z2 * e2_4pi_eps0     # MeV fm
r_min_headon = k_coul / E_alpha    # fm
k_over_4E = k_coul / (4 * E_alpha)
dsdo_90 = k_over_4E**2 / np.sin(np.pi / 4)**4           # fm^2/sr at 90 degrees
sigma_gt90 = np.pi * (k_coul / (2 * E_alpha))**2        # fm^2: every b < b(90 deg) = k/2E, (3.28)
rho_Au, A_Au, N_A = 19.32, 196.967, 6.02214076e23       # g/cm^3, g/mol
n_Au = rho_Au / A_Au * N_A                              # atoms / cm^3
t_foil = 0.4e-4                                         # 0.4 micrometre in cm
frac_gt90 = n_Au * t_foil * sigma_gt90 * 1e-26          # 1 fm^2 = 1e-26 cm^2
R_Au = 1.2 * A_Au ** (1 / 3)
ratio_15 = (np.sin(np.radians(15 / 2)) / np.sin(np.radians(150 / 2)))**4   # dsigma(150)/dsigma(15)
ratio_5_150 = (np.sin(np.radians(75)) / np.sin(np.radians(2.5)))**4         # dsigma(5)/dsigma(150)
E_cm_factor = A_Au / (A_Au + 4.0026)                    # E_cm/E_lab
sigma_hard = np.pi * R_Au**2                            # a hard sphere of the nuclear radius, (3.25)
theta_thomson = k_coul / (E_alpha * 1e5)                # Thomson atom: deflection ~ k/(E R_atom), R_atom ~ 1e5 fm
layers = t_foil / n_Au ** (-1 / 3)                      # atomic layers in the foil
theta_rms = np.sqrt(layers) * theta_thomson
log10_P90 = -(np.pi / 2 / theta_rms)**2 / np.log(10)
print(f"Part 9, Ex 4.6 Rutherford: k = {k_coul:.1f} MeV fm, head-on r_min = k/E = {r_min_headon:.1f} fm (gold nuclear radius ~ {R_Au:.1f} fm); E_cm/E_lab = {E_cm_factor:.4f}")
print(f"        k/2E = {2*k_over_4E:.2f} fm, k/4E = {k_over_4E:.2f} fm, (k/4E)^2 = {k_over_4E**2:.1f} fm^2/sr = {k_over_4E**2/100:.2f} b/sr; d sigma/d Omega(90 deg) = {dsdo_90:.0f} fm^2/sr = {dsdo_90/100:.2f} b/sr")
print(f"        sigma(theta > 90 deg) = pi (k/2E)^2 = {sigma_gt90:.0f} fm^2 = {sigma_gt90/100:.1f} b; gold n = {n_Au:.3e} /cm^3; "
      f"fraction for a 0.4 um foil = {frac_gt90:.2e} (about 1 in {1/frac_gt90:,.0f})")
print(f"        dsigma(150)/dsigma(15) = 1/{1/ratio_15:,.0f}; dsigma(5)/dsigma(150) = {ratio_5_150:,.0f}; hard sphere R = {R_Au:.1f} fm: sigma = {sigma_hard:.0f} fm^2")
print(f"        Thomson atom: theta per atom ~ {theta_thomson:.2e} rad, {layers:.0f} layers, rms {theta_rms:.4f} rad = {np.degrees(theta_rms):.2f} deg, "
      f"P(>90 deg) ~ 10^{log10_P90:.0f}")

# Halley's comet (reference values for "Your turn" 3)
a_H, e_H = 17.8 * AU, 0.967
q_H = a_H * (1 - e_H)
v_pH = np.sqrt(GM_sun * (2 / q_H - 1 / a_H))
T_H = kepler_period(a_H, 1.0, GM_sun) / yr
print(f"Part 9, Halley (Your turn 3): q = {q_H/AU:.3f} AU, Q = {a_H*(1+e_H)/AU:.1f} AU, v_peri = {v_pH/1e3:.1f} km/s, "
      f"v_apo = {v_pH*q_H/(a_H*(1+e_H))/1e3:.2f} km/s, T = {T_H:.1f} yr")

# %% [markdown]
# ## Part 10: Symbolic checks of §3, and reference values for the exercises
#
# SymPy re-derives the chain-rule steps of Binet's equation (3.10), the Kepler
# solution (3.11), the energy–eccentricity relation (3.14), Kepler III (3.16),
# the stability ratio (3.21), the hard-sphere and Rutherford cross-sections
# (3.25), (3.27), and the LRL magnitude (3.29). Then every number quoted in
# `solutions.md` is printed. Do the exercises before reading this cell.

# %%
import sympy as sp
ph, mu_s, l_s, k_s, E_s, e_s, a_s, r_s, R_s, r0_s = sp.symbols("phi mu ell k E e a r R r0", positive=True)
th, n_s = sp.symbols("theta n", real=True)
u = sp.Function("u")(ph)
dt_dphi = mu_s / (l_s * u**2)                                        # phi-dot = ell u^2 / mu
rdot = sp.diff(1 / u, ph) / dt_dphi
rddot = sp.diff(rdot, ph) / dt_dphi
assert sp.simplify(rdot + l_s / mu_s * u.diff(ph)) == 0                     # Binet step 1
assert sp.simplify(rddot + l_s**2 * u**2 / mu_s**2 * u.diff(ph, 2)) == 0    # Binet step 2
uK = mu_s * k_s / l_s**2 * (1 + e_s * sp.cos(ph))
assert sp.simplify(uK.diff(ph, 2) + uK - mu_s * k_s / l_s**2) == 0          # (3.11) solved
up = uK.subs(ph, 0)
assert sp.simplify(l_s**2 * up**2 / (2 * mu_s) - k_s * up - mu_s * k_s**2 * (e_s**2 - 1) / (2 * l_s**2)) == 0   # (3.14)
T_area = 2 * mu_s * sp.pi * a_s * a_s * sp.sqrt(1 - e_s**2) / sp.sqrt(mu_s * k_s * a_s * (1 - e_s**2))
assert sp.simplify(T_area - 2 * sp.pi * sp.sqrt(mu_s * a_s**3 / k_s)) == 0  # (3.16)
l2 = mu_s * k_s * r0_s**(n_s + 3)                                           # circular orbit of F = -k r^n
Veff = k_s * r_s**(n_s + 1) / (n_s + 1) + l2 / (2 * mu_s * r_s**2)
assert sp.simplify(sp.diff(Veff, r_s).subs(r_s, r0_s)) == 0
assert sp.simplify(sp.diff(Veff, r_s, 2).subs(r_s, r0_s) / mu_s / (l2 / (mu_s**2 * r0_s**4)) - (n_s + 3)) == 0   # (3.21)
x = sp.symbols("x", positive=True)
assert sp.series(2 * sp.pi * (1 / sp.sqrt(1 - x) - 1), x, 0, 2).removeO() == sp.pi * x                        # (3.22)
b_ruth, b_hard = k_s / (2 * E_s) * sp.cot(th / 2), R_s * sp.cos(th / 2)
ds_ruth = sp.lambdify((th, k_s, E_s), b_ruth / sp.sin(th) * sp.Abs(sp.diff(b_ruth, th)))
ds_hard = sp.lambdify((th, R_s), b_hard / sp.sin(th) * sp.Abs(sp.diff(b_hard, th)))
tt = np.linspace(0.1, 3.0, 30)
assert np.allclose(ds_ruth(tt, 1.3, 0.7), (1.3 / 2.8)**2 / np.sin(tt / 2)**4, rtol=1e-12)                    # (3.27)
assert np.allclose(ds_hard(tt, 1.7), 1.7**2 / 4, rtol=1e-12)                                                   # (3.25)
rp_s = l_s**2 / (mu_s * k_s * (1 + e_s))
assert sp.simplify(mu_s * (l_s / (mu_s * rp_s)) * l_s - mu_s * k_s - mu_s * k_s * e_s) == 0                    # (3.29)
print("Part 10: SymPy confirms Binet's chain rule, (3.11), (3.14), (3.16), (3.21), (3.22), (3.25), (3.27), (3.29)")

# A1: reduced masses
mp_me, md_me, Ry = 1836.15267343, 3670.48296788, 13.605693
red_H, red_D = mp_me / (1 + mp_me), md_me / (1 + md_me)
m_E, m_M, d_EM = 5.9722e24, 7.342e22, 3.844e8
print(f"A1: mu_H/m_e = {red_H:.8f}, mu_D/m_e = {red_D:.8f}; H ground state {-Ry*red_H:.4f} eV, positronium {-Ry/2:.4f} eV; "
      f"H-alpha isotope shift {656.28*(1-red_H/red_D):.4f} nm")
print(f"    Earth-Moon: mu/m_Moon = {m_E/(m_E+m_M):.5f}; barycentre {d_EM*m_M/(m_E+m_M)/1e3:.0f} km from Earth's centre")
# A2: r = 2/(1 + 0.5 cos phi), mu = k = 1
cA, eA = 2.0, 0.5
aA = cA / (1 - eA**2); lA = np.sqrt(cA); EA = -1 / (2 * aA)
assert abs(ecc_from_E_ell(EA, lA, 1.0, 1.0) - eA) < 1e-12
print(f"A2: r_p = {cA/(1+eA):.4f}, r_a = {cA/(1-eA):.4f}, a = {aA:.4f}, b = {aA*np.sqrt(1-eA**2):.4f}, ell = {lA:.4f}, E = {EA:.4f}, "
      f"T = {kepler_period(aA, 1, 1):.2f}, v_p = {lA/(cA/(1+eA)):.4f}, v_a = {lA/(cA/(1-eA)):.4f} "
      f"(sqrt(1-e^2) = {np.sqrt(1-eA**2):.5f}, a^1.5 = {aA**1.5:.4f})")
# A3: Mercury
vm = np.sqrt(GM_sun / a_merc)
print(f"A3: Mercury sqrt(GM/a) = {vm/1e3:.2f} km/s, sqrt((1+e)/(1-e)) = {np.sqrt((1+e_merc)/(1-e_merc)):.5f}, v_p = {vm*np.sqrt((1+e_merc)/(1-e_merc))/1e3:.2f} km/s, v_a = {vm*np.sqrt((1-e_merc)/(1+e_merc))/1e3:.2f} km/s, "
      f"ratio {(1+e_merc)/(1-e_merc):.4f}")
# A4: Jupiter
q_J, T_J = 1.26686534e17 / GM_sun, 4332.59
print(f"A4: m_J/M = {q_J:.4e}, sqrt(1+q) - 1 = {np.sqrt(1+q_J)-1:.3e}; Sun-only formula overestimates Jupiter's period by {T_J*(np.sqrt(1+q_J)-1):.2f} d")
# A5
print(f"A5: sigma(theta > 1 deg)/sigma(theta > 90 deg) = cot^2(0.5 deg) = {1/np.tan(np.radians(0.5))**2:,.0f}")
# B2: the circle r = 2R cos(phi) through the centre needs F = -8 ell^2 R^2/(mu r^5), and has E = 0
uc = 1 / (2 * R_s * sp.cos(ph))
assert sp.simplify(uc.diff(ph, 2) + uc - 8 * R_s**2 * uc**3) == 0
assert sp.simplify(uc.diff(ph)**2 + uc**2 - 4 * R_s**2 * uc**4) == 0          # v^2 = (ell/mu)^2 (u'^2 + u^2)
print("B2: SymPy: u'' + u = 8 R^2 u^3 and u'^2 + u^2 = 4 R^2 u^4 on r = 2R cos(phi)")
# B3: Molniya
rpM, raM = (6371 + 600) * 1e3, (6371 + 39700) * 1e3
aM, eM = 0.5 * (rpM + raM), (raM - rpM) / (raM + rpM)
vpM, vaM = np.sqrt(GM_E * (2 / rpM - 1 / aM)), np.sqrt(GM_E * (2 / raM - 1 / aM))
assert abs(vpM * rpM - vaM * raM) / (vpM * rpM) < 1e-12
print(f"B3: Molniya a = {aM/1e3:.0f} km, e = {eM:.4f}, T = {kepler_period(aM, 1, GM_E)/3600:.3f} h, v_p = {vpM/1e3:.3f} km/s, v_a = {vaM/1e3:.3f} km/s, "
      f"r_p v_p = r_a v_a = {rpM*vpM/1e6:.4e} km^2/s")
# B4: Yukawa, units k = lambda = 1: stability condition d ln f/d ln r > -3
f_y = sp.exp(-r_s) * (1 / r_s**2 + 1 / r_s)
cond = sp.factor(sp.simplify(r_s * sp.diff(sp.log(f_y), r_s) + 3))
root = [s_ for s_ in sp.solve(sp.numer(cond), r_s) if s_.is_positive][0]
assert sp.simplify(root - (1 + sp.sqrt(5)) / 2) == 0
print(f"B4: Yukawa: d ln f/d ln r + 3 = {cond}; stable iff r0/lambda < {float(root):.6f}")
# B6: closest approach
rmin = lambda th_deg, kk, EE: kk / (2 * EE) * (1 + 1 / np.sin(np.radians(th_deg) / 2))
k_Al, d_Al = 2 * 13 * e2_4pi_eps0, 1.2 * (27 ** (1 / 3) + 4 ** (1 / 3))
assert abs(rmin(180, k_coul, E_alpha) - r_min_headon) < 1e-12
print(f"B6: gold, 5 MeV: r_min(60 deg) = {rmin(60, k_coul, E_alpha):.1f} fm, r_min(180 deg) = {rmin(180, k_coul, E_alpha):.1f} fm; "
      f"Al: k = {k_Al:.2f} MeV fm, contact {d_Al:.2f} fm, threshold E = {k_Al/d_Al:.2f} MeV (at 5 MeV head-on r_min = {k_Al/E_alpha:.2f} fm)")
# B7: Earth -> Mars Hohmann (circular, coplanar)
r1M, r2M = AU, 1.524 * AU
atM = 0.5 * (r1M + r2M)
dv1M = np.sqrt(GM_sun * (2 / r1M - 1 / atM)) - np.sqrt(GM_sun / r1M)
dv2M = np.sqrt(GM_sun / r2M) - np.sqrt(GM_sun * (2 / r2M - 1 / atM))
t_tM = 0.5 * kepler_period(atM, 1, GM_sun)
T_Mars = kepler_period(r2M, 1, GM_sun)
lead = 180 - 360 * t_tM / T_Mars
syn = 1 / (1 / kepler_period(r1M, 1, GM_sun) - 1 / T_Mars)
print(f"B7: v_E = {np.sqrt(GM_sun/r1M)/1e3:.2f}, v_p = {np.sqrt(GM_sun*(2/r1M-1/atM))/1e3:.2f}, v_a = {np.sqrt(GM_sun*(2/r2M-1/atM))/1e3:.2f}, "
      f"v_Mars = {np.sqrt(GM_sun/r2M)/1e3:.2f} km/s; a_t = {atM/AU:.3f} AU; T_Mars = {T_Mars/86400:.1f} d, Mars moves {360*t_tM/T_Mars:.1f} deg in transit")
print(f"B7: Earth->Mars: dv1 = {dv1M/1e3:.3f} km/s, dv2 = {dv2M/1e3:.3f} km/s, transfer {t_tM/86400:.1f} d, Mars leads by {lead:.1f} deg, "
      f"windows every {syn/86400:.0f} d")
# B8: LRL along the precessing orbit of Part 5 (mu = k = 1)
xs, ys, vxs, vys = s_prec.y
ls = xs * vys - ys * vxs
Axs, Ays = ls * vys - xs / np.hypot(xs, ys), -ls * vxs - ys / np.hypot(xs, ys)
xe, ye_, vxe, vye = ye_prec.T
le = xe * vye - ye_ * vxe
angA = np.unwrap(np.arctan2(-le * vxe - ye_ / np.hypot(xe, ye_), le * vye - xe / np.hypot(xe, ye_)))
assert np.max(np.abs(np.diff(angA) - prec_exact)) < 1e-6
assert np.max(np.abs(np.angle(np.exp(1j * (angA - np.arctan2(ye_, xe)))))) < 1e-6   # A points at each perihelion
u_c = mu * k / (ell0**2 * beta**2)                     # centre of the u oscillation, Theorem 3.11
e_pert = (1 / rp0) / u_c - 1
A_ap = mu * k - ell0**2 * u_c * (1 - e_pert)           # |A| at aphelion, where A = (l^2/r - mu k) r_hat
assert abs(A_ap - np.hypot(Axs, Ays).min()) < 1e-6
print(f"B8: mu k/(l^2 beta^2) = {u_c:.5f}, e' = {e_pert:.5f}, u_a = {u_c*(1-e_pert):.5f}, r_a = {1/(u_c*(1-e_pert)):.5f}, "
      f"l^2 u_a = {ell0**2*u_c*(1-e_pert):.4f}, |A|/(mu k) at aphelion = {A_ap:.4f}")
print(f"B8: LRL at successive perihelia advances {np.mean(np.diff(angA)):.6f} rad/orbit (precession {prec_exact:.6f}); "
      f"|A|/(mu k) ranges over [{np.hypot(Axs, Ays).min():.4f}, {np.hypot(Axs, Ays).max():.4f}] within an orbit")
# C1: modified Binet equation w'' + w = 1 + delta w^2 (w = u h^2/GM), e = 0.2: advance per orbit vs 2 pi delta
delta = 1e-3
def peri(p, w): return w[1]
peri.direction = -1
sw = solve_ivp(lambda p, w: [w[1], 1 + delta * w[0]**2 - w[0]], (0, 12 * np.pi), [1.2, 0.0], events=peri,
               rtol=1e-12, atol=1e-12, method="DOP853")
adv = np.diff(sw.t_events[0]) - 2 * np.pi
assert abs(adv.mean() / (2 * np.pi * delta) - 1) < 0.01
print(f"C1: modified Binet, delta = {delta}: advance per orbit {adv.mean():.6f} rad vs first-order 2 pi delta = {2*np.pi*delta:.6f} "
      f"(ratio {adv.mean()/(2*np.pi*delta):.5f}; near-circular second order 1 + 2.5 delta = {1+2.5*delta:.5f})")
# C2: Kepler's equation, time averages, Molniya dwell time
Mgrid = np.linspace(0, 2 * np.pi, 200001)[:-1]
psi = Mgrid.copy()
for _ in range(50):
    psi -= (psi - eM * np.sin(psi) - Mgrid) / (1 - eM * np.cos(psi))
rr_t = 1 - eM * np.cos(psi)                       # r/a, sampled uniformly in time
assert abs(rr_t.mean() - (1 + eM**2 / 2)) < 1e-9 and abs((1 / rr_t).mean() - 1) < 1e-9
print(f"C2: Molniya <r>/a = {rr_t.mean():.6f} = 1 + e^2/2; fraction of time with r > a: {np.mean(rr_t > 1):.4f} (1/2 + e/pi = {0.5+eM/np.pi:.4f})")
# C3: repulsive alpha/r^2, alpha = 1, E = 1/2 (mu = 1): theta(b) = pi (1 - b/sqrt(b^2 + alpha/E))
th_c3 = lambda bb: np.pi * (1 - bb / np.sqrt(bb**2 + 2.0))
ds_c3 = lambda T: 2.0 * np.pi**2 * (np.pi - T) / (T**2 * (2 * np.pi - T)**2 * np.sin(T))
bb = np.linspace(0.2, 5, 50); h_ = 1e-6
num = bb / np.sin(th_c3(bb)) / np.abs((th_c3(bb + h_) - th_c3(bb - h_)) / (2 * h_))
assert np.allclose(num, ds_c3(th_c3(bb)), rtol=1e-6)
R0 = 200.0; vR = np.sqrt(1 - 2 / R0**2); bR = 1.0 / vR       # E = 1/2 and ell = b = 1 exactly at launch
sc3 = solve_ivp(lambda t, Y: [Y[2], Y[3], 2 * Y[0] / np.hypot(Y[0], Y[1])**4, 2 * Y[1] / np.hypot(Y[0], Y[1])**4],
                (0, 2 * R0 + 50), [-R0, bR, vR, 0], rtol=1e-11, atol=1e-12, method="DOP853")
th_num = np.arctan2(sc3.y[3, -1], sc3.y[2, -1])
assert abs(th_num - th_c3(1.0)) < 1e-3
print(f"C3: b = 1: theta formula {np.degrees(th_c3(1.0)):.2f} deg, integrated {np.degrees(th_num):.2f} deg; dsigma/dOmega(90 deg) = {ds_c3(np.pi/2):.4f}")
# C4: L4 is an equilibrium of the rotating-frame potential; Routh's stability bound
m1c, m2c = 1.0, 9.546e-4; Mc = m1c + m2c; Om2 = Mc                        # G = d = 1
r1c, r2c = np.array([-m2c / Mc, 0]), np.array([m1c / Mc, 0])
L4 = 0.5 * (r1c + r2c) + np.array([0, np.sqrt(3) / 2])
grad = m1c * (L4 - r1c) / np.linalg.norm(L4 - r1c)**3 + m2c * (L4 - r2c) / np.linalg.norm(L4 - r2c)**3 - Om2 * L4
assert np.max(np.abs(grad)) < 1e-14
print(f"C4: |grad Phi_eff| at L4 = {np.max(np.abs(grad)):.1e}; Routh bound m2/M < {(1-np.sqrt(23/27))/2:.5f}; Sun-Jupiter {m2c/Mc:.3e}, "
      f"Earth-Moon {m_M/(m_E+m_M):.5f}")

# %% [markdown]
# ## Your turn
#
# 1. **Restricted three-body problem.** Put the Sun ($m_1 = 1$) and Jupiter
#    ($m_2 = 1/1047$) on a circular orbit of radius 1 (units $G = 1$), and add a
#    massless test particle. Integrate in the *rotating* frame (Coriolis and
#    centrifugal terms: $\ddot{\mathbf r} = -\nabla\Phi + 2\Omega\,\dot{\mathbf r}\times\hat{\mathbf z}
#    + \Omega^2\mathbf r$, with $\Omega^2 = G(m_1+m_2)/d^3$). Start the particle at
#    the equilateral point L4 (60° ahead of Jupiter) with a small displacement
#    and confirm it stays nearby (librates) for 100 Jupiter years. Then do the
#    same at L1 (between Sun and Jupiter, find it by solving
#    $\partial\Phi_{\rm eff}/\partial x = 0$ numerically) and watch it leave.
#    Exercise C4 shows why L4 is an equilibrium; stability needs
#    $m_2/(m_1+m_2) < 0.0385$, which Jupiter satisfies.
# 2. **Orbital decay under drag.** Add a drag acceleration $-\gamma |\mathbf v|\,\mathbf v$
#    to the Kepler integrator in Part 2 (small $\gamma$). Record $a$ and $e$
#    from $E$ and $\ell$ at each perihelion. You should see $a$ shrink and $e$
#    fall (drag circularizes), and $T^2 \propto a^3$ still hold at each instant.
#    Why does the speed *increase* as energy is lost? (Vis-viva.)
# 3. **Halley's comet.** From $a = 17.8$ AU and $e = 0.967$ compute $q$, $Q$, the
#    period, and the perihelion and aphelion speeds with vis-viva; then integrate
#    one full orbit and confirm the speeds and the period. Reference values are
#    printed in Part 9.

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.perf_counter() - T_START:.1f} s. Figures in {FIG}")
