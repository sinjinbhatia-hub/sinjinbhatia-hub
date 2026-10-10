# %% [markdown]
# # CM-4 Project: Symmetry and Conservation
#
# Two halves. **Symbolic (SymPy):** Noether's theorem as a function.
# `noether_charge(L, q, v, K)` checks the invariance condition (3.3) and returns
# the charge Q = sum_i p_i K_i of Theorem 3.3; `noether_charge_quasi` handles the
# total-derivative case of Theorem 3.4 by *constructing* F; `energy_function`
# builds h of (3.8); `on_shell_dt` differentiates any quantity along solutions of
# Euler–Lagrange, so every "is conserved" claim in the notes is checked exactly.
# **Numerical (SciPy):** integrate the Kepler problem and watch E, L and the
# Laplace–Runge–Lenz vector stay flat; add a 1/r^3 term and watch the orbit
# precess while E and L stay put; integrate the bead on the rotating hoop and
# see h flat while T + V is not; verify mechanical similarity (Kepler's third
# law, amplitude-independence of the harmonic period) and the virial theorem.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp, quad
from pathlib import Path

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
# Fixed categorical colors (assigned in order, never cycled).
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.25, "lines.linewidth": 1.6})

# %% [markdown]
# ## Part 1: Noether's theorem in SymPy
#
# Coordinates are plain symbols: `q_i`, velocities `v_i` (for q̇_i), accelerations
# `a_i` (for q̈_i), and `t`. A point transformation is a list K of expressions
# K_i(q, t). The induced velocity change is K̇_i = sum_j (∂K_i/∂q_j) v_j + ∂K_i/∂t,
# and the invariance condition (3.3) is
#     δL = sum_i (∂L/∂q_i K_i + ∂L/∂v_i K̇_i) = 0   for all (q, v, t).
# If it holds, Q = sum_i p_i K_i is conserved (Theorem 3.3). If instead δL = dF/dt
# for some F(q, t), then Q = sum_i p_i K_i - F is conserved (Theorem 3.4).
#
# **Build it yourself first:** write `delta_L` and `noether_charge` before reading on.

# %%
t = sp.Symbol("t")


def coords(names):
    """Symbols (q, v, a) for generalized coordinates, velocities, accelerations."""
    q = [sp.Symbol(n) for n in names]
    v = [sp.Symbol(n + "dot") for n in names]
    a = [sp.Symbol(n + "ddot") for n in names]
    return q, v, a


def total_dt(expr, q, v, a):
    """Total time derivative of expr(q, v, t) along a path:
    d expr/dt = sum_i (∂expr/∂q_i v_i + ∂expr/∂v_i a_i) + ∂expr/∂t."""
    return sum(sp.diff(expr, qi) * vi + sp.diff(expr, vi) * ai
               for qi, vi, ai in zip(q, v, a)) + sp.diff(expr, t)


def delta_L(L, q, v, K, a=None):
    """First-order change of L under q_i -> q_i + eps K_i, eq. (3.3):
    δL = sum_i (∂L/∂q_i K_i + ∂L/∂v_i K̇_i).  K may depend on (q, t); if it also
    depends on v (velocity-dependent generator, §3.4 and C1) pass the acceleration
    symbols a so that K̇ can be formed."""
    if a is None:
        a = [sp.Integer(0)] * len(q)
    Kdot = [total_dt(Ki, q, v, a) for Ki in K]
    return sp.simplify(sum(sp.diff(L, qi) * Ki + sp.diff(L, vi) * Kdi
                           for qi, vi, Ki, Kdi in zip(q, v, K, Kdot)))


def noether_charge(L, q, v, K):
    """Noether's theorem for a point symmetry (Theorem 3.3).
    Checks δL == 0 and returns Q = sum_i p_i K_i, p_i = ∂L/∂v_i.
    Raises ValueError (carrying δL) if the transformation is not a symmetry."""
    dL = delta_L(L, q, v, K)
    if dL != 0:
        raise ValueError(f"not a symmetry of L: δL = {dL}")
    return sp.simplify(sum(sp.diff(L, vi) * Ki for vi, Ki in zip(v, K)))


def total_derivative_potential(dL, q, v):
    """Given δL(q, v, t) linear in v, construct F(q, t) with dF/dt = δL (Theorem 3.4).
    Write δL = sum_i a_i(q, t) v_i + b(q, t). dF/dt = δL means ∂F/∂q_i = a_i and
    ∂F/∂t = b: (a, b) must be a gradient field on (q, t)-space. Poincaré's lemma
    (star-shaped domain) gives F(q, t) = ∫_0^1 [sum_i a_i(sq, st) q_i + b(sq, st) t] ds.
    The result is verified; ValueError if δL is not a total derivative."""
    coef = [sp.diff(dL, vi) for vi in v]
    if any(c.has(vj) for c in coef for vj in v):
        raise ValueError("δL is not linear in the velocities; no F(q, t) exists")
    b = sp.simplify(dL - sum(c * vi for c, vi in zip(coef, v)))
    s = sp.Symbol("s", positive=True)
    scale = {qi: s * qi for qi in q}
    scale[t] = s * t
    integrand = sum(c.subs(scale, simultaneous=True) * qi for c, qi in zip(coef, q)) \
        + b.subs(scale, simultaneous=True) * t
    F = sp.simplify(sp.integrate(sp.expand(integrand), (s, 0, 1)))
    check = sp.simplify(sum(sp.diff(F, qi) * vi for qi, vi in zip(q, v)) + sp.diff(F, t) - dL)
    if check != 0:
        raise ValueError(f"δL = {dL} is not a total time derivative of any F(q, t)")
    return F


def noether_charge_quasi(L, q, v, K):
    """Theorem 3.4: if δL = dF/dt, return (Q, F) with Q = sum_i p_i K_i - F."""
    dL = delta_L(L, q, v, K)
    F = sp.Integer(0) if dL == 0 else total_derivative_potential(dL, q, v)
    Q = sp.simplify(sum(sp.diff(L, vi) * Ki for vi, Ki in zip(v, K)) - F)
    return Q, F


def el_accelerations(L, q, v, a):
    """Solve Euler–Lagrange d/dt(∂L/∂v_i) = ∂L/∂q_i for the accelerations.
    d/dt(∂L/∂v_i) = sum_j W_ij a_j + sum_j (∂²L/∂v_i∂q_j) v_j + ∂²L/∂v_i∂t with the
    Hessian W_ij = ∂²L/∂v_i∂v_j. L is *regular* when W is invertible (used in
    Theorem 3.5); then a = W^{-1}(∂L/∂q - mixed terms)."""
    n = len(q)
    W = sp.Matrix(n, n, lambda i, j: sp.diff(L, v[i], v[j]))
    rhs = sp.Matrix([sp.diff(L, q[i])
                     - sum(sp.diff(L, v[i], q[j]) * v[j] for j in range(n))
                     - sp.diff(L, v[i], t) for i in range(n)])
    acc = W.LUsolve(rhs)
    return {ai: sp.simplify(acc[i]) for i, ai in enumerate(a)}


def on_shell_dt(expr, L, q, v, a):
    """d expr/dt evaluated along solutions of the Euler–Lagrange equations of L."""
    return sp.simplify(total_dt(expr, q, v, a).subs(el_accelerations(L, q, v, a)))


def energy_function(L, q, v):
    """Energy function h = sum_i v_i ∂L/∂v_i - L, eq. (3.8)."""
    return sp.simplify(sum(vi * sp.diff(L, vi) for vi in v) - L)


# %% [markdown]
# ### 1a. Cyclic coordinate: central force in polar coordinates (Example in §3.1)

# %%
m, kappa, g, R, omega = sp.symbols("m kappa g R omega", positive=True)
Vf = sp.Function("V")  # an arbitrary potential

(r, th), (rdot, thdot), (rddot, thddot) = coords(["r", "theta"])
L_polar = sp.Rational(1, 2) * m * (rdot**2 + r**2 * thdot**2) - Vf(r)
Q_rot = noether_charge(L_polar, [r, th], [rdot, thdot], K=[0, 1])      # theta -> theta + eps
assert sp.simplify(Q_rot - m * r**2 * thdot) == 0
assert on_shell_dt(Q_rot, L_polar, [r, th], [rdot, thdot], [rddot, thddot]) == 0
print("1a  polar central force: Q =", Q_rot, " (p_theta; dQ/dt = 0 on shell)")

# %% [markdown]
# ### 1b. Two-body translation and the Galilean boost (Example 4.2, 1D version)
# Translation leaves L invariant: Q = total momentum. The boost x_a -> x_a + eps t
# changes L by a total derivative; `noether_charge` refuses, `noether_charge_quasi`
# constructs F = m1 x1 + m2 x2 and returns Q = P t - M X (the center-of-mass theorem).

# %%
m1, m2 = sp.symbols("m1 m2", positive=True)
(x1, x2), (v1, v2), (a1, a2) = coords(["x1", "x2"])
L_2body = sp.Rational(1, 2) * m1 * v1**2 + sp.Rational(1, 2) * m2 * v2**2 - Vf(x1 - x2)
P = noether_charge(L_2body, [x1, x2], [v1, v2], K=[1, 1])
assert sp.simplify(P - (m1 * v1 + m2 * v2)) == 0
assert on_shell_dt(P, L_2body, [x1, x2], [v1, v2], [a1, a2]) == 0
print("1b  two-body translation: Q =", P)

try:
    noether_charge(L_2body, [x1, x2], [v1, v2], K=[t, t])
    raise AssertionError("boost should not be a strict symmetry")
except ValueError as err:
    print("1b  boost, strict check:", err)
G_boost, F_boost = noether_charge_quasi(L_2body, [x1, x2], [v1, v2], K=[t, t])
assert sp.simplify(F_boost - (m1 * x1 + m2 * x2)) == 0
assert sp.simplify(G_boost - (t * (m1 * v1 + m2 * v2) - (m1 * x1 + m2 * x2))) == 0
assert on_shell_dt(G_boost, L_2body, [x1, x2], [v1, v2], [a1, a2]) == 0
print("1b  boost (quasi-symmetry): F =", F_boost, "   Q =", G_boost)

# %% [markdown]
# ### 1c. Rotation about an axis for a 3D particle in a central potential (Example 4.1)
# K = n̂ × r. The charge should be n̂ · (r × p): the scalar triple product identity
# p · (n × r) = n · (r × p) done by SymPy.

# %%
(x, y, z), (vx, vy, vz), (ax, ay, az) = coords(["x", "y", "z"])
rvec = sp.Matrix([x, y, z]); vvec = sp.Matrix([vx, vy, vz])
L_3d = sp.Rational(1, 2) * m * vvec.dot(vvec) - Vf(sp.sqrt(rvec.dot(rvec)))

Lz = noether_charge(L_3d, [x, y, z], [vx, vy, vz], K=list(sp.Matrix([0, 0, 1]).cross(rvec)))
assert sp.simplify(Lz - m * (x * vy - y * vx)) == 0
assert on_shell_dt(Lz, L_3d, [x, y, z], [vx, vy, vz], [ax, ay, az]) == 0
print("1c  rotation about z: Q =", Lz)

n1, n2, n3 = sp.symbols("n1 n2 n3", real=True)
nhat = sp.Matrix([n1, n2, n3])
Ln = noether_charge(L_3d, [x, y, z], [vx, vy, vz], K=list(nhat.cross(rvec)))
assert sp.simplify(Ln - nhat.dot(rvec.cross(m * vvec))) == 0
print("1c  rotation about n: Q = n·(r × p) =", sp.factor(Ln))

try:  # a translation is NOT a symmetry of a central potential
    noether_charge(L_3d, [x, y, z], [vx, vy, vz], K=[1, 0, 0])
    raise AssertionError
except ValueError as err:
    print("1c  translation of a central-force particle:", str(err)[:60], "...")

# %% [markdown]
# ### 1d. Time translation and the energy function (§3.4)
# Bead on the rotating hoop: ∂L/∂t = 0, so h is conserved (3.9), but h ≠ T + V.
# Driven oscillator: dh/dt = -∂L/∂t exactly. Time translation as a Noether symmetry
# with the velocity-dependent generator K = q̇: δL = dL/dt, so F = L and Q = h.

# %%
(theta,), (thetadot,), (thetaddot,) = coords(["theta"])
T_hoop = sp.Rational(1, 2) * m * R**2 * (thetadot**2 + omega**2 * sp.sin(theta)**2)
V_hoop = -m * g * R * sp.cos(theta)
L_hoop = T_hoop - V_hoop
h_hoop = energy_function(L_hoop, [theta], [thetadot])
E_hoop = T_hoop + V_hoop
assert sp.simplify(h_hoop - (sp.Rational(1, 2) * m * R**2 * thetadot**2
                             - sp.Rational(1, 2) * m * R**2 * omega**2 * sp.sin(theta)**2
                             - m * g * R * sp.cos(theta))) == 0
assert on_shell_dt(h_hoop, L_hoop, [theta], [thetadot], [thetaddot]) == 0
dE_hoop = on_shell_dt(E_hoop, L_hoop, [theta], [thetadot], [thetaddot])
assert sp.simplify(dE_hoop - 2 * m * R**2 * omega**2 * sp.sin(theta) * sp.cos(theta) * thetadot) == 0
assert sp.simplify(E_hoop - h_hoop - m * R**2 * omega**2 * sp.sin(theta)**2) == 0
# time translation as a Noether quasi-symmetry with K = thetadot (needs accelerations):
dL_time = delta_L(L_hoop, [theta], [thetadot], K=[thetadot], a=[thetaddot])
assert sp.simplify(dL_time - total_dt(L_hoop, [theta], [thetadot], [thetaddot])) == 0  # F = L
assert sp.simplify(sp.diff(L_hoop, thetadot) * thetadot - L_hoop - h_hoop) == 0        # Q = h
print("1d  rotating hoop: h =", h_hoop)
print("    dh/dt = 0 on shell;  dE/dt =", dE_hoop, " (motor power)")

(xq,), (xv,), (xa,) = coords(["x"])
k_s, F0, Om = sp.symbols("k F_0 Omega", positive=True)
L_drive = sp.Rational(1, 2) * m * xv**2 - sp.Rational(1, 2) * k_s * xq**2 + xq * F0 * sp.cos(Om * t)
h_drive = energy_function(L_drive, [xq], [xv])
assert sp.simplify(on_shell_dt(h_drive, L_drive, [xq], [xv], [xa]) + sp.diff(L_drive, t)) == 0
print("1d  driven oscillator: dh/dt = -∂L/∂t =", sp.simplify(-sp.diff(L_drive, t)))

# %% [markdown]
# ## Part 2: Numerical conservation monitoring for the Kepler problem
#
# Units m = κ = 1. Orbit with a = 1, e = 0.5, tilted so every component of L is
# nonzero. E = -κ/(2a), |L| = sqrt(mκa(1-e²)), |A| = mκe (Example 4.6).


# %%
def kepler_rhs(t_, s, kappa=1.0, m=1.0, eps=0.0):
    """r'' = -(κ/r³ + 3ε/r⁵) r / m, i.e. V = -κ/r - ε/r³ (ε = 0: pure Kepler)."""
    rr, vv = s[:3], s[3:]
    rn = np.linalg.norm(rr)
    return np.concatenate([vv, -(kappa / rn**3 + 3 * eps / rn**5) * rr / m])


def invariants(S, kappa=1.0, m=1.0):
    """E, L (3-vector) and the Laplace–Runge–Lenz vector A = p × L - mκ r̂ along a trajectory.
    S has shape (N, 6): rows (x, y, z, vx, vy, vz)."""
    rr, vv = S[:, :3], S[:, 3:]
    rn = np.linalg.norm(rr, axis=1)
    E = 0.5 * m * np.sum(vv**2, axis=1) - kappa / rn
    Lv = np.cross(rr, m * vv)
    A = np.cross(m * vv, Lv) - m * kappa * rr / rn[:, None]
    return E, Lv, A


def rot_x(ang):
    c, s_ = np.cos(ang), np.sin(ang)
    return np.array([[1, 0, 0], [0, c, -s_], [0, s_, c]])


def rot_z(ang):
    c, s_ = np.cos(ang), np.sin(ang)
    return np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])


def kepler_ic(a=1.0, e=0.5, kappa=1.0, m=1.0, tilt=np.radians(35), node=np.radians(20)):
    """State at pericenter of an orbit with semi-major axis a, eccentricity e, rotated out of the xy-plane."""
    rp = a * (1 - e)
    vp = np.sqrt(kappa * (1 + e) / (m * rp))
    Rm = rot_z(node) @ rot_x(tilt)
    return np.concatenate([Rm @ np.array([rp, 0, 0]), Rm @ np.array([0, vp, 0])])


def integrate(rhs, s0, t_end, n_out=4000, events=None, **kw):
    ts = np.linspace(0, t_end, n_out)
    sol = solve_ivp(rhs, (0, t_end), s0, t_eval=ts, method="DOP853",
                    rtol=1e-11, atol=1e-13, events=events, **kw)
    assert sol.success, sol.message
    return sol


a0, e0 = 1.0, 0.5
T0 = 2 * np.pi * np.sqrt(a0**3)          # period for m = κ = 1
s0 = kepler_ic(a0, e0)
sol = integrate(kepler_rhs, s0, 10 * T0)
E, Lv, A = invariants(sol.y.T)

assert np.isclose(E[0], -0.5 / a0), E[0]
assert np.isclose(np.linalg.norm(Lv[0]), np.sqrt(a0 * (1 - e0**2)))
assert np.isclose(np.linalg.norm(A[0]), e0)
assert np.all(np.abs(Lv[0]) > 0.1), "all three components of L should be nonzero for this test"
assert np.allclose(np.einsum("ij,ij->i", A, Lv), 0, atol=1e-9), "A ⟂ L"
dE = np.abs(E - E[0]).max() / abs(E[0])
dL = np.abs(Lv - Lv[0]).max()
dA = np.abs(A - A[0]).max()
assert dE < 1e-8 and dL < 1e-8 and dA < 1e-8, (dE, dL, dA)
print(f"Part 2  Kepler, 10 periods: max |ΔE|/|E| = {dE:.1e}, max |ΔL| = {dL:.1e}, max |ΔA| = {dA:.1e}")

fig, axs = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
axs[0].plot(sol.t / T0, E - E[0], color=C[0])
axs[0].set_ylabel("E − E(0)")
for i, lab in enumerate("xyz"):
    axs[1].plot(sol.t / T0, Lv[:, i] - Lv[0, i], color=C[i], label=f"$L_{lab}$")
    axs[2].plot(sol.t / T0, A[:, i] - A[0, i], color=C[i], label=f"$A_{lab}$")
axs[1].set_ylabel("L − L(0)"); axs[1].legend(ncol=3, fontsize=8)
axs[2].set_ylabel("A − A(0)"); axs[2].legend(ncol=3, fontsize=8)
axs[2].set_xlabel("time / period")
axs[0].set_title("Kepler problem: energy, angular momentum, Laplace–Runge–Lenz vector (deviations)")
fig.tight_layout(); fig.savefig(FIG / "kepler_conservation.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 3: Break the hidden symmetry
#
# V = -κ/r - ε/r³ is still central, so E and L stay conserved, but the orbit is no
# longer closed: the pericenter advances each revolution and A rotates with it.
# First-order perturbation theory (Landau §15, problem 3, with γ = -ε) predicts a
# shift per revolution of Δφ = 6π κ ε m² / L⁴.


# %%
def pericenter(t_, s):
    """r · v changes sign from - to + at pericenter."""
    return np.dot(s[:3], s[3:])


pericenter.direction = 1

eps = 0.002
n_orb = 12
sol_p = integrate(lambda t_, s: kepler_rhs(t_, s, eps=eps), s0, n_orb * T0, n_out=12000, events=pericenter)
E_p, L_p, A_p = invariants(sol_p.y.T)     # Kepler E and A evaluated along the perturbed orbit
E_true = E_p - eps / np.linalg.norm(sol_p.y[:3].T, axis=1)**3   # the real energy includes the -ε/r³ term
assert np.abs(E_true - E_true[0]).max() < 1e-8
assert np.abs(L_p - L_p[0]).max() < 1e-8
assert np.abs(A_p - A_p[0]).max() > 0.1, "the LRL vector should visibly move"

# in-plane basis: e1 = r̂(0), e2 = L̂ × e1 (L is conserved, so the plane is fixed)
Lhat = L_p[0] / np.linalg.norm(L_p[0])
e1 = sol_p.y[:3, 0] / np.linalg.norm(sol_p.y[:3, 0])
e2 = np.cross(Lhat, e1)
peri_states = sol_p.y_events[0]
peri_ang = np.unwrap(np.arctan2(peri_states[:, :3] @ e2, peri_states[:, :3] @ e1))
shift_per_orbit = np.mean(np.diff(peri_ang))
Lmag = np.linalg.norm(L_p[0])
shift_pred = 6 * np.pi * eps / Lmag**4
rel_err = abs(shift_per_orbit - shift_pred) / shift_pred
assert rel_err < 0.1, (shift_per_orbit, shift_pred)
print(f"Part 3  ε = {eps}: pericenter advance per orbit = {np.degrees(shift_per_orbit):.3f}°, "
      f"first-order prediction 6πκεm²/L⁴ = {np.degrees(shift_pred):.3f}° (rel. err {rel_err:.1%})")

X, Y = sol_p.y[:3].T @ e1, sol_p.y[:3].T @ e2
fig, (axo, axa) = plt.subplots(1, 2, figsize=(11, 4.6))
axo.plot(X, Y, color=C[0], lw=0.7)
axo.plot(peri_states[:, :3] @ e1, peri_states[:, :3] @ e2, "o", color=C[1], ms=4, label="pericenters")
axo.plot(0, 0, "k*", ms=9)
axo.set_aspect("equal"); axo.legend(fontsize=8)
axo.set_title(f"V = −κ/r − ε/r³, ε = {eps}: {n_orb} orbits precess")
axa.plot(sol_p.t / T0, A_p @ e1, color=C[0], label="$A_1$ (in plane)")
axa.plot(sol_p.t / T0, A_p @ e2, color=C[1], label="$A_2$ (in plane)")
axa.plot(sol_p.t / T0, E_true - E_true[0], color=C[2], label="E − E(0)")
axa.plot(sol_p.t / T0, np.linalg.norm(L_p, axis=1) - Lmag, color=C[3], label="|L| − |L|(0)")
axa.set_xlabel("time / period"); axa.legend(fontsize=8)
axa.set_title("A rotates; E and L stay flat")
fig.tight_layout(); fig.savefig(FIG / "lrl_precession.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 4: Bead on a rotating hoop — h is conserved, T + V is not (Example 4.3)
#
# θ'' = ω² sin θ cos θ − (g/R) sin θ. h = ½mR²θ̇² − ½mR²ω² sin²θ − mgR cos θ (flat),
# E = T + V = ½mR²θ̇² + ½mR²ω² sin²θ − mgR cos θ (not flat), E − h = mR²ω² sin²θ.

# %%
m_h, R_h, g_h, w_h = 1.0, 1.0, 9.81, 2.0
th0, thd0 = 1.0, 0.0


def hoop_rhs(t_, s):
    th_, thd_ = s
    return [thd_, w_h**2 * np.sin(th_) * np.cos(th_) - (g_h / R_h) * np.sin(th_)]


def hoop_h(th_, thd_):
    return 0.5 * m_h * R_h**2 * thd_**2 - 0.5 * m_h * R_h**2 * w_h**2 * np.sin(th_)**2 - m_h * g_h * R_h * np.cos(th_)


def hoop_E(th_, thd_):
    return 0.5 * m_h * R_h**2 * thd_**2 + 0.5 * m_h * R_h**2 * w_h**2 * np.sin(th_)**2 - m_h * g_h * R_h * np.cos(th_)


sol_h = integrate(hoop_rhs, [th0, thd0], 10.0, n_out=4000)
th_t, thd_t = sol_h.y
h_t, E_t = hoop_h(th_t, thd_t), hoop_E(th_t, thd_t)
assert np.abs(h_t - h_t[0]).max() < 1e-9
assert E_t.max() - E_t.min() > 2.0, "E should swing by E(0) - h ≈ 2.83"
assert np.allclose(E_t - h_t, m_h * R_h**2 * w_h**2 * np.sin(th_t)**2)
assert np.isclose(h_t[0], -6.7165, atol=5e-4) and np.isclose(E_t[0], -3.8842, atol=5e-4)
# motor power = dE/dt = 2 m R² ω² sinθ cosθ θ̇ (checked symbolically in 1d; here numerically)
dEdt_num = np.gradient(E_t, sol_h.t)
dEdt_formula = 2 * m_h * R_h**2 * w_h**2 * np.sin(th_t) * np.cos(th_t) * thd_t
assert np.abs(dEdt_num - dEdt_formula)[5:-5].max() < 0.05 * np.abs(dEdt_formula).max()
print(f"Part 4  hoop: h = {h_t[0]:.4f} (flat to {np.abs(h_t - h_t[0]).max():.1e}); "
      f"E(0) = {E_t[0]:.4f}, E ranges over [{E_t.min():.4f}, {E_t.max():.4f}]")

fig, ax_ = plt.subplots(figsize=(8, 4))
ax_.plot(sol_h.t, h_t, color=C[0], label="h = Σ q̇ ∂L/∂q̇ − L  (conserved)")
ax_.plot(sol_h.t, E_t, color=C[1], label="E = T + V  (not conserved)")
ax_.set_xlabel("t [s]"); ax_.set_ylabel("energy [J]"); ax_.legend(fontsize=9)
ax_.set_title("Bead on a hoop rotating at ω = 2 rad/s: energy function vs mechanical energy")
fig.tight_layout(); fig.savefig(FIG / "rotating_hoop_h_vs_E.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5: Mechanical similarity (§3.6)
#
# V homogeneous of degree k: scaling lengths by α and times by α^{1−k/2} maps
# solutions to solutions, so periods scale as T ∝ (size)^{1−k/2}.
# k = −1: T² ∝ a³ (Kepler III). k = 2: T independent of amplitude. k = 4: T ∝ 1/A.

# %%
alphas = np.array([0.5, 1.0, 2.0, 4.0])
periods = []
for al in alphas:
    s_al = np.concatenate([al * s0[:3], al**-0.5 * s0[3:]])    # r -> α r, v -> α^{k/2} v with k = -1
    sol_al = integrate(kepler_rhs, s_al, 2.6 * T0 * al**1.5, n_out=200, events=pericenter)
    te = sol_al.t_events[0]
    assert len(te) >= 2, te
    periods.append(np.mean(np.diff(te)))
periods = np.array(periods)
ratio = periods**2 / alphas**3
assert np.allclose(ratio, 4 * np.pi**2, rtol=1e-6), ratio
print("Part 5  Kepler: a =", alphas, " T²/a³ =", np.round(ratio, 6), " (4π² =", round(4 * np.pi**2, 6), ")")
print(f"        Mars check: (1.5237 AU)^{{3/2}} = {1.5237**1.5:.4f} yr  (observed 1.8808 yr)")


def period_quad(k, A, m=1.0):
    """Period of m x'' = -|x|^{k-1} sign(x), V = |x|^k / k, amplitude A:
    T = 4 ∫_0^A dx / sqrt(2 (V(A) − V(x)) / m)."""
    V = lambda xx: abs(xx)**k / k
    val, _ = quad(lambda xx: 1.0 / np.sqrt(2 * (V(A) - V(xx)) / m), 0, A, limit=200)
    return 4 * val


amps = np.logspace(-1, 1, 7)
slopes = {}
for k_hom, col in [(2, C[0]), (4, C[1])]:
    Tk = np.array([period_quad(k_hom, A) for A in amps])
    slope = np.polyfit(np.log(amps), np.log(Tk), 1)[0]
    slopes[k_hom] = (Tk, slope)
    assert np.isclose(slope, 1 - k_hom / 2, atol=1e-6), (k_hom, slope)
assert np.allclose(slopes[2][0], 2 * np.pi, rtol=1e-6)          # harmonic: T = 2π for m = k = 1
# ODE cross-check for the quartic oscillator at A = 1


def quartic_rhs(t_, s):
    return [s[1], -s[0]**3]


def x_zero(t_, s):
    return s[0]


x_zero.direction = -1
sol_q = solve_ivp(quartic_rhs, (0, 20), [1.0, 0.0], method="DOP853", rtol=1e-11, atol=1e-13, events=x_zero, dense_output=True)
T_quartic_ode = 4 * sol_q.t_events[0][0]
T_quartic_quad = period_quad(4, 1.0)
assert np.isclose(T_quartic_ode, T_quartic_quad, rtol=1e-7)
print(f"        quartic oscillator, A = 1: T = {T_quartic_quad:.6f} (quadrature) vs {T_quartic_ode:.6f} (ODE); "
      f"log-log slopes: k=2 → {slopes[2][1]:+.4f}, k=4 → {slopes[4][1]:+.4f}")

fig, ax_ = plt.subplots(figsize=(7, 4.5))
ax_.loglog(amps, slopes[2][0], "o-", color=C[0], label="k = 2 (harmonic): slope 0")
ax_.loglog(amps, slopes[4][0], "s-", color=C[1], label="k = 4 (quartic): slope −1")
ax_.loglog(alphas, periods, "^-", color=C[2], label="k = −1 (Kepler): slope 3/2")
ax_.set_xlabel("size (amplitude A or semi-major axis a)"); ax_.set_ylabel("period")
ax_.set_title("Mechanical similarity: period $\\propto$ size$^{1-k/2}$"); ax_.legend(fontsize=9)
fig.tight_layout(); fig.savefig(FIG / "mechanical_similarity.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 6: The virial theorem, 2⟨T⟩ = k⟨V⟩ (§3.7)
#
# Time averages over one period: Kepler (k = −1) should give 2⟨T⟩ = −⟨V⟩ and
# E = −⟨T⟩; the quartic oscillator 2⟨T⟩ = 4⟨V⟩; the harmonic ⟨T⟩ = ⟨V⟩.


# %%
def time_average(f, ts):
    return np.trapezoid(f, ts) / (ts[-1] - ts[0])


sol_1 = integrate(kepler_rhs, s0, T0, n_out=20001)
rr = sol_1.y[:3].T; vv = sol_1.y[3:].T
T_kin = 0.5 * np.sum(vv**2, axis=1); V_pot = -1.0 / np.linalg.norm(rr, axis=1)
Tavg, Vavg = time_average(T_kin, sol_1.t), time_average(V_pot, sol_1.t)
assert np.isclose(2 * Tavg, -Vavg, rtol=1e-6), (Tavg, Vavg)
assert np.isclose(E[0], -Tavg, rtol=1e-6)
print(f"Part 6  Kepler: ⟨T⟩ = {Tavg:.6f}, ⟨V⟩ = {Vavg:.6f}, 2⟨T⟩ + ⟨V⟩ = {2*Tavg + Vavg:.1e}, E = {E[0]:.6f} = −⟨T⟩")

tq = np.linspace(0, T_quartic_ode, 20001)
xq_, vq_ = sol_q.sol(tq)
Tq, Vq = 0.5 * vq_**2, xq_**4 / 4
assert np.isclose(2 * time_average(Tq, tq), 4 * time_average(Vq, tq), rtol=1e-6)
print(f"        quartic: 2⟨T⟩ = {2*time_average(Tq, tq):.6f}, 4⟨V⟩ = {4*time_average(Vq, tq):.6f}")
tho = np.linspace(0, 2 * np.pi, 20001); xo = np.cos(tho); vo = -np.sin(tho)
assert np.isclose(time_average(0.5 * vo**2, tho), time_average(0.5 * xo**2, tho))
print(f"        harmonic: ⟨T⟩ = {time_average(0.5*vo**2, tho):.6f} = ⟨V⟩ = {time_average(0.5*xo**2, tho):.6f}")

# Example 4.5(b): the Sun. 2⟨T⟩ = -⟨V⟩ with V = -(3/5) G M²/R (uniform sphere) and
# T = (3/2) N k_B T_mean, N = M / (μ m_p), μ ≈ 0.6 for ionized solar composition.
G_N, M_sun, R_sun, m_p, k_B, mu = 6.674e-11, 1.989e30, 6.957e8, 1.6726e-27, 1.3807e-23, 0.6
T_sun = G_N * M_sun * mu * m_p / (5 * k_B * R_sun)
assert 2.5e6 < T_sun < 3.0e6
print(f"        Sun (virial, uniform sphere): k_B T = G M μ m_p / (5R)  →  T ≈ {T_sun:.2e} K")

# running averages figure: 2<T> and -<V> converge to each other over a Kepler period
run = lambda f: np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(sol_1.t)) / sol_1.t[1:]
fig, ax_ = plt.subplots(figsize=(7, 4))
ax_.plot(sol_1.t[1:] / T0, 2 * run(T_kin), color=C[0], label="running 2⟨T⟩")
ax_.plot(sol_1.t[1:] / T0, -run(V_pot), color=C[1], label="running −⟨V⟩")
ax_.axhline(-E[0], color="k", lw=0.8, ls="--", label="−E")
ax_.set_xlabel("time / period"); ax_.legend(fontsize=9)
ax_.set_title("Virial theorem for a Kepler orbit: 2⟨T⟩ = −⟨V⟩ after one full period")
fig.tight_layout(); fig.savefig(FIG / "virial_kepler.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Your turn
#
# 1. **Isotropic 2D oscillator** L = ½m(ẋ² + ẏ²) − ½mω²(x² + y²). Rotation gives
#    L_z = m(xẏ − yẋ). Show that Q₂ = m(ẋẏ + ω²xy) and Q₁ = ½m(ẋ² − ẏ²) + ½mω²(x² − y²)
#    are also conserved (use `on_shell_dt`), then find velocity-dependent generators
#    K for them (hint for Q₂: K = (ẏ, ẋ); δL is a total derivative of something that
#    depends on velocities, so Theorem 3.4 still applies but `total_derivative_potential`
#    does not; check by hand). Verify E² = Q₁² + Q₂² + ω²L_z² numerically. (C1.)
# 2. **Double pendulum** (CM-3): L = ½(m₁+m₂)l₁²θ̇₁² + ½m₂l₂²θ̇₂² + m₂l₁l₂θ̇₁θ̇₂cos(θ₁−θ₂)
#    + (m₁+m₂)gl₁cosθ₁ + m₂gl₂cosθ₂. Integrate it, and plot h, p_θ₁, p_θ₂ and
#    p_θ₁ + p_θ₂ against time. Only h is flat. Which symmetry is missing, and which
#    term in L breaks it?
# 3. **V = −c/r²** (C3). Work out α^{1−k/2} for k = −2, then check numerically that
#    two orbits related by the scaling have the *same* angular momentum. Integrate
#    with ℓ² slightly above and slightly below 2mc and describe what happens.
# 4. Change `eps` in Part 3 to 0.0005, 0.001, 0.004 and plot the measured advance
#    per orbit against ε. Where does the first-order formula start to fail? (B8.)

# %%
# Reference check for item 1 (do it yourself first).
(ox, oy), (ovx, ovy), (oax, oay) = coords(["x", "y"])
w = sp.Symbol("omega", positive=True)
L_iso = sp.Rational(1, 2) * m * (ovx**2 + ovy**2) - sp.Rational(1, 2) * m * w**2 * (ox**2 + oy**2)
Q2 = m * (ovx * ovy + w**2 * ox * oy)
assert on_shell_dt(Q2, L_iso, [ox, oy], [ovx, ovy], [oax, oay]) == 0
dL_Q2 = delta_L(L_iso, [ox, oy], [ovx, ovy], K=[ovy, ovx], a=[oax, oay])
F_Q2 = m * ovx * ovy - m * w**2 * ox * oy
assert sp.simplify(dL_Q2 - total_dt(F_Q2, [ox, oy], [ovx, ovy], [oax, oay])) == 0
assert sp.simplify(sp.diff(L_iso, ovx) * ovy + sp.diff(L_iso, ovy) * ovx - F_Q2 - Q2) == 0
print("Your turn 1 (reference): Q2 = m(ẋẏ + ω²xy) is conserved; generator K = (ẏ, ẋ), F =", F_Q2)

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
