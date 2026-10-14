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


_ACCEL = set()   # every acceleration symbol ever created, so a δL can be screened for them


def coords(names):
    """Symbols (q, v, a) for generalized coordinates, velocities, accelerations."""
    q = [sp.Symbol(n) for n in names]
    v = [sp.Symbol(n + "dot") for n in names]
    a = [sp.Symbol(n + "ddot") for n in names]
    _ACCEL.update(a)
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
    K = [sp.sympify(Ki) for Ki in K]
    if a is None:
        if any(Ki.has(*v) for Ki in K):
            raise ValueError("K depends on the velocities: pass the acceleration symbols a so that K̇ can be formed")
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
    The result is verified; ValueError if δL is not a total derivative, or if it
    contains accelerations (velocity-dependent generator: any F would depend on
    velocities, outside this construction; see C1)."""
    dL = sp.sympify(dL)
    if dL.free_symbols & _ACCEL:
        raise ValueError("δL contains accelerations (velocity-dependent generator); no F(q, t) exists. "
                         "Theorem 3.4 may still hold with an F that depends on velocities: check by hand (C1)")
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

# Example 4.2 numbers: m1 = 1, m2 = 3, r1(0) = 0, r2(0) = (4, 0, 0), v1(0) = (2, 1, 0), v2(0) = (0, -1, 0).
m1n, m2n = 1.0, 3.0
r1_0, r2_0 = np.array([0.0, 0.0, 0.0]), np.array([4.0, 0.0, 0.0])
v1_0, v2_0 = np.array([2.0, 1.0, 0.0]), np.array([0.0, -1.0, 0.0])
P_num = m1n * v1_0 + m2n * v2_0
R_0 = (m1n * r1_0 + m2n * r2_0) / (m1n + m2n)
G_num = P_num * 0.0 - (m1n + m2n) * R_0                 # G = P t − M R, evaluated at t = 0
R_2 = (P_num * 2.0 - G_num) / (m1n + m2n)               # R(t) = (P t − G)/M, whatever V is
assert np.allclose(P_num, [2, -2, 0]) and np.allclose(R_0, [3, 0, 0]) and np.allclose(G_num, [-12, 0, 0])
assert np.allclose(R_2, [4, -1, 0])
print("1b  Example 4.2 numbers: P =", P_num, " R(0) =", R_0, " G =", G_num, " R(2) =", R_2)

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
assert np.isclose(G_N * M_sun / R_sun, 1.908e11, rtol=1e-3)          # GM/R in Example 4.5(b)
assert np.isclose(T_sun, 2.77e6, rtol=1e-2) and 2.5e6 < T_sun < 3.0e6
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

# %% [markdown]
# ## Reference checks for the exercises
#
# Every number quoted in `solutions.md` is produced below, and every symbolic
# claim there is re-derived by `delta_L`, `noether_charge_quasi` and `on_shell_dt`.
# Do the exercise on paper first; then run this and compare.

# %%
# --- A1(c): helical potential V(r, z - a*theta). Nothing is cyclic, but the screw
#     motion K = (0, 1, a) is a symmetry with charge p_theta + a p_z.
a_h = sp.Symbol("a", positive=True)
Vh = sp.Function("V")
(rh, thh, zh), (rhd, thhd, zhd), _ = coords(["r", "theta", "z"])
L_helix = sp.Rational(1, 2) * m * (rhd**2 + rh**2 * thhd**2 + zhd**2) - Vh(rh, zh - a_h * thh)
Q_helix = noether_charge(L_helix, [rh, thh, zh], [rhd, thhd, zhd], K=[0, 1, a_h])
assert sp.simplify(Q_helix - (m * rh**2 * thhd + a_h * m * zhd)) == 0
print("A1  helical potential: screw generator K = (0, 1, a) gives Q =", Q_helix)

# --- A2(c): the magnetic term drops out of h.
qe, Bf = sp.symbols("q_e B", positive=True)
L_mag = sp.Rational(1, 2) * m * (vx**2 + vy**2) + qe * Bf * x * vy        # Landau gauge A = (0, Bx, 0)
assert sp.simplify(energy_function(L_mag, [x, y], [vx, vy]) - sp.Rational(1, 2) * m * (vx**2 + vy**2)) == 0

# --- A5: four candidate generators for L = ½m(ẋ² + ẏ²) − V(y).
L_Vy = sp.Rational(1, 2) * m * (vx**2 + vy**2) - Vf(y)
dVy = sp.diff(Vf(y), y)
assert delta_L(L_Vy, [x, y], [vx, vy], K=[1, 0]) == 0
assert sp.simplify(delta_L(L_Vy, [x, y], [vx, vy], K=[0, 1]) + dVy) == 0
assert sp.simplify(delta_L(L_Vy, [x, y], [vx, vy], K=[-y, x]) + x * dVy) == 0
dL_shear = delta_L(L_Vy, [x, y], [vx, vy], K=[y, 0])
assert sp.simplify(dL_shear - m * vx * vy) == 0
try:
    total_derivative_potential(dL_shear, [x, y], [vx, vy])
    raise AssertionError("m ẋ ẏ must not be a total derivative of any F(q, t)")
except ValueError:
    pass
print("A5  shear K = (y, 0): δL =", dL_shear, " (quadratic in velocities: not a quasi-symmetry)")

# --- A6: time translation K = ẋ for L = ½mẋ² − V(x): δL = dL/dt, F = L, Q = h.
L_1d = sp.Rational(1, 2) * m * xv**2 - Vf(xq)
dL_t = delta_L(L_1d, [xq], [xv], K=[xv], a=[xa])
assert sp.simplify(dL_t - total_dt(L_1d, [xq], [xv], [xa])) == 0
assert sp.simplify(sp.diff(L_1d, xv) * xv - L_1d - energy_function(L_1d, [xq], [xv])) == 0

# %%
# --- B1: uniform gravity L = ½m|v|² − mgz.
L_grav = sp.Rational(1, 2) * m * (vx**2 + vy**2 + vz**2) - m * g * z
assert noether_charge(L_grav, [x, y, z], [vx, vy, vz], K=[1, 0, 0]) == m * vx
Q_z, F_z = noether_charge_quasi(L_grav, [x, y, z], [vx, vy, vz], K=[0, 0, 1])
assert sp.simplify(F_z + m * g * t) == 0 and sp.simplify(Q_z - (m * vz + m * g * t)) == 0
Q_bz, F_bz = noether_charge_quasi(L_grav, [x, y, z], [vx, vy, vz], K=[0, 0, t])
assert sp.simplify(F_bz - (m * z - m * g * t**2 / 2)) == 0
assert sp.simplify(Q_bz - (m * vz * t - m * z + m * g * t**2 / 2)) == 0
for Qg in (Q_z, Q_bz):
    assert on_shell_dt(Qg, L_grav, [x, y, z], [vx, vy, vz], [ax, ay, az]) == 0
assert sp.simplify(noether_charge(L_grav, [x, y, z], [vx, vy, vz], K=[-y, x, 0]) - m * (x * vy - y * vx)) == 0
for Kbad in ([0, -z, y], [z, 0, -x]):            # rotations about x and y are not even quasi-symmetries
    try:
        noether_charge_quasi(L_grav, [x, y, z], [vx, vy, vz], K=Kbad)
        raise AssertionError
    except ValueError:
        pass
print("B1  gravity: z-translation F =", F_z, ", Q =", Q_z, ";  z-boost F =", F_bz, ", Q =", Q_bz)

# --- B2: two-body separation and L_tot = M R×V + μ s×ṡ.
Mtot = m1 + m2
mu = m1 * m2 / Mtot
VX, vs_ = sp.symbols("V_X v_s", real=True)
T_two = sp.Rational(1, 2) * m1 * v1**2 + sp.Rational(1, 2) * m2 * v2**2
T_sep = T_two.subs({v1: VX + m2 / Mtot * vs_, v2: VX - m1 / Mtot * vs_}, simultaneous=True)
assert sp.simplify(T_sep - (sp.Rational(1, 2) * Mtot * VX**2 + sp.Rational(1, 2) * mu * vs_**2)) == 0
Rv = sp.Matrix(sp.symbols("R1 R2 R3", real=True)); Vv = sp.Matrix(sp.symbols("V1 V2 V3", real=True))
sv = sp.Matrix(sp.symbols("s1 s2 s3", real=True)); sdv = sp.Matrix(sp.symbols("sd1 sd2 sd3", real=True))
Ltot = m1 * (Rv + m2 / Mtot * sv).cross(Vv + m2 / Mtot * sdv) + m2 * (Rv - m1 / Mtot * sv).cross(Vv - m1 / Mtot * sdv)
assert sp.simplify(Ltot - (Mtot * Rv.cross(Vv) + mu * sv.cross(sdv))) == sp.zeros(3, 1)
print("B2  two-body: T = ½M V² + ½μ v_s² and L_tot = M R×V + μ s×ṡ verified symbolically")

# --- B3: charged particle in B = B ẑ, Landau gauge (CM-3 Example 4.8).
py_mag = noether_charge(L_mag, [x, y], [vx, vy], K=[0, 1])
assert sp.simplify(py_mag - (m * vy + qe * Bf * x)) == 0
Qx_mag, Fx_mag = noether_charge_quasi(L_mag, [x, y], [vx, vy], K=[1, 0])
assert sp.simplify(Fx_mag - qe * Bf * y) == 0 and sp.simplify(Qx_mag - (m * vx - qe * Bf * y)) == 0
Qrot_mag, Frot_mag = noether_charge_quasi(L_mag, [x, y], [vx, vy], K=[-y, x])
assert sp.simplify(Frot_mag - qe * Bf * (x**2 - y**2) / 2) == 0
assert sp.simplify(Qrot_mag - (m * (x * vy - y * vx) + qe * Bf * (x**2 + y**2) / 2)) == 0
for Qm in (py_mag, Qx_mag, Qrot_mag):
    assert on_shell_dt(Qm, L_mag, [x, y], [vx, vy], [ax, ay]) == 0
Xc, Yc = py_mag / (qe * Bf), -Qx_mag / (qe * Bf)                       # guiding centre
rhoL2 = (x - Xc)**2 + (y - Yc)**2
assert sp.simplify(rhoL2 - m**2 * (vx**2 + vy**2) / (qe * Bf)**2) == 0   # Larmor radius² = (mv/q_eB)²
assert sp.simplify(Qrot_mag - qe * Bf * (Xc**2 + Yc**2 - rhoL2) / 2) == 0
# symmetric gauge: rotation is strict, both translations quasi, same charges
L_sym = sp.Rational(1, 2) * m * (vx**2 + vy**2) + qe * Bf * (x * vy - y * vx) / 2
assert sp.simplify(noether_charge(L_sym, [x, y], [vx, vy], K=[-y, x]) - Qrot_mag) == 0
Qxs, Fxs = noether_charge_quasi(L_sym, [x, y], [vx, vy], K=[1, 0])
Qys, Fys = noether_charge_quasi(L_sym, [x, y], [vx, vy], K=[0, 1])
assert sp.simplify(Fxs - qe * Bf * y / 2) == 0 and sp.simplify(Qxs - Qx_mag) == 0
assert sp.simplify(Fys + qe * Bf * x / 2) == 0 and sp.simplify(Qys - py_mag) == 0
print("B3  magnetic field: p_y =", py_mag, ";  x-translation F =", Fx_mag, ", Q =", Qx_mag)
print("    rotation: F =", Frot_mag, ", Q =", Qrot_mag, " = ½q_eB(R_c² − ρ_L²)")

# --- B4: anisotropic oscillator; δL for a rotation equals dL_z/dt on shell.
w1, w2 = sp.symbols("omega1 omega2", positive=True)
L_aniso = sp.Rational(1, 2) * m * (ovx**2 + ovy**2) - sp.Rational(1, 2) * m * (w1**2 * ox**2 + w2**2 * oy**2)
Lz_xy = m * (ox * ovy - oy * ovx)
dL_rot = delta_L(L_aniso, [ox, oy], [ovx, ovy], K=[-oy, ox])
assert sp.simplify(dL_rot - m * (w1**2 - w2**2) * ox * oy) == 0
assert sp.simplify(on_shell_dt(Lz_xy, L_aniso, [ox, oy], [ovx, ovy], [oax, oay]) - dL_rot) == 0
Ex_, Ey_ = (sp.Rational(1, 2) * m * ovx**2 + sp.Rational(1, 2) * m * w1**2 * ox**2,
            sp.Rational(1, 2) * m * ovy**2 + sp.Rational(1, 2) * m * w2**2 * oy**2)
assert on_shell_dt(Ex_, L_aniso, [ox, oy], [ovx, ovy], [oax, oay]) == 0
assert on_shell_dt(Ey_, L_aniso, [ox, oy], [ovx, ovy], [oax, oay]) == 0
Q3 = (ovx**2 - w1**2 * ox**2) * ovy + 4 * w1**2 * ox * ovx * oy       # cubic charge for ω2 = 2ω1
L_12 = L_aniso.subs(w2, 2 * w1)
assert on_shell_dt(Q3, L_12, [ox, oy], [ovx, ovy], [oax, oay]) == 0
assert on_shell_dt(Q3, L_aniso, [ox, oy], [ovx, ovy], [oax, oay]) != 0
print("B4  anisotropic oscillator: δL(rotation) = dL_z/dt =", dL_rot, ";  Q3 conserved iff ω2 = 2ω1")

# --- B5: bead on a horizontal rod rotating at ω: h flat, E = T grows.
(rb,), (rbd,), (rbdd,) = coords(["rho"])
L_rod = sp.Rational(1, 2) * m * (rbd**2 + omega**2 * rb**2)
h_rod = energy_function(L_rod, [rb], [rbd])
assert sp.simplify(h_rod - sp.Rational(1, 2) * m * (rbd**2 - omega**2 * rb**2)) == 0
assert on_shell_dt(h_rod, L_rod, [rb], [rbd], [rbdd]) == 0
assert sp.simplify(on_shell_dt(L_rod, L_rod, [rb], [rbd], [rbdd]) - 2 * m * omega**2 * rb * rbd) == 0  # dE/dt = τω
assert el_accelerations(L_rod, [rb], [rbd], [rbdd])[rbdd] == omega**2 * rb
m_r, w_r, r0_r, t1 = 0.5, 2.0, 0.1, 1.0
r1, rd1 = r0_r * np.cosh(w_r * t1), r0_r * w_r * np.sinh(w_r * t1)
h_rod_num = 0.5 * m_r * (rd1**2 - w_r**2 * r1**2)
E_rod_0, E_rod_1 = 0.5 * m_r * w_r**2 * r0_r**2, 0.5 * m_r * (rd1**2 + w_r**2 * r1**2)
P_motor = m_r * w_r**3 * r0_r**2 * np.sinh(2 * w_r * t1)
assert np.isclose(h_rod_num, -0.5 * m_r * w_r**2 * r0_r**2)
assert np.isclose(E_rod_1, 0.5 * m_r * w_r**2 * r0_r**2 * np.cosh(2 * w_r * t1))
assert np.isclose(E_rod_1 - h_rod_num, m_r * w_r**2 * r1**2)
sol_rod = solve_ivp(lambda t_, s: [s[1], w_r**2 * s[0]], (0, t1), [r0_r, 0.0], method="DOP853", rtol=1e-11, atol=1e-13)
assert np.allclose(sol_rod.y[:, -1], [r1, rd1], rtol=1e-8)
print(f"B5  rod (m = {m_r}, ω = {w_r}, r₀ = {r0_r}): at t = 1 s  r = {r1:.5f} m, ṙ = {rd1:.5f} m/s, "
      f"h = {h_rod_num:.4f} J, E(0) = {E_rod_0:.4f} J, E(1 s) = {E_rod_1:.5f} J, motor power = {P_motor:.4f} W")
print(f"    cosh 2 = {np.cosh(2):.5f}, sinh 2 = {np.sinh(2):.5f}, cosh 4 = {np.cosh(4):.4f}, sinh 4 = {np.sinh(4):.4f}")

# --- B6: quartic period in closed form; T(A) ∝ A^{-1}; Kepler ⟨1/r⟩ = 1/a.
from scipy.special import gamma as Gamma
T_quartic_exact = np.sqrt(2 * np.pi) * Gamma(0.25) / Gamma(0.75)
assert np.isclose(T_quartic_quad, T_quartic_exact, rtol=1e-9)
assert np.isclose(period_quad(4, 2.0), T_quartic_exact / 2, rtol=1e-8)
inv_r_avg = time_average(1 / np.linalg.norm(rr, axis=1), sol_1.t)
assert np.isclose(inv_r_avg, 1 / a0, rtol=1e-6)
r_avg = time_average(np.linalg.norm(rr, axis=1), sol_1.t)
assert np.isclose(r_avg, a0 * (1 + e0**2 / 2), rtol=1e-6)                 # ⟨r⟩_t = a(1 + e²/2), not a
Tq_avg, Vq_avg = time_average(Tq, tq), time_average(Vq, tq)
assert np.isclose(Tq_avg, 2 / 3 * 0.25, rtol=1e-6) and np.isclose(Vq_avg, 1 / 12, rtol=1e-6)
print(f"B6  quartic: T(1) = {T_quartic_exact:.6f} = √(2π)Γ(1/4)/Γ(3/4), T(2) = {period_quad(4, 2.0):.6f}; "
      f"⟨T⟩ = {Tq_avg:.6f} = 2E/3, ⟨V⟩ = {Vq_avg:.6f} = E/3;  Kepler ⟨1/r⟩ = {inv_r_avg:.6f} = 1/a, ⟨r⟩ = {r_avg:.6f} = a(1 + e²/2)")
print(f"    Γ(1/4) = {Gamma(0.25):.6f}, Γ(3/4) = {Gamma(0.75):.6f}, √(2π) = {np.sqrt(2 * np.pi):.6f}")

# --- B7: A² = m²κ² + 2mEL² along the orbit; e from E and L.
A2 = np.sum(A**2, axis=1); L2 = np.sum(Lv**2, axis=1)
assert np.allclose(A2, 1 + 2 * E * L2, atol=1e-9)
e_from_EL = np.sqrt(1 + 2 * E[0] * L2[0])
assert np.isclose(e_from_EL, e0)
print(f"B7  A² − m²κ² − 2mEL² = {np.abs(A2 - 1 - 2 * E * L2).max():.1e} along the orbit;  e = √(1 + 2EL²/mκ²) = {e_from_EL:.4f}")

# %%
# --- B8: pericenter advance vs ε (extends Part 3).
eps_list = np.array([0.0005, 0.001, 0.002, 0.004])
adv_meas = []
for ep in eps_list:
    sol_e = integrate(lambda t_, s, ep=ep: kepler_rhs(t_, s, eps=ep), s0, 8 * T0, n_out=800, events=pericenter)
    ps = sol_e.y_events[0]
    ang = np.unwrap(np.arctan2(ps[:, :3] @ e2, ps[:, :3] @ e1))
    adv_meas.append(np.mean(np.diff(ang)))
adv_meas = np.degrees(adv_meas)
adv_pred = np.degrees(6 * np.pi * eps_list / Lmag**4)
rel_err = (adv_meas - adv_pred) / adv_pred
c_fit = np.sum(eps_list * rel_err) / np.sum(eps_list**2)   # least-squares slope of a line through the origin
assert np.isclose(adv_meas[2], np.degrees(shift_per_orbit), rtol=1e-3)     # same as Part 3
assert np.all(rel_err > 0) and np.all(np.diff(rel_err) > 0), rel_err      # error grows with ε
assert np.allclose(rel_err, c_fit * eps_list, atol=0.004), (rel_err, c_fit * eps_list)   # ... linearly
print("B8  ε            =", eps_list)
print("    measured [°] =", np.round(adv_meas, 3))
print("    1st order [°]=", np.round(adv_pred, 3))
print("    rel. error   =", np.round(rel_err, 4), f"  ≈ {c_fit:.1f}·ε  (least-squares line through the origin)")
print(f"    first-order formula off by 1% at ε ≈ {0.01 / c_fit:.5f}, by 10% at ε ≈ {0.1 / c_fit:.4f} (linear extrapolation)")
fig, (axl, axr) = plt.subplots(1, 2, figsize=(10, 4))
axl.plot(eps_list, adv_meas, "o", color=C[0], label="measured")
axl.plot(eps_list, adv_pred, "-", color=C[1], label="6πκεm²/L⁴")
axl.set_xlabel("ε"); axl.set_ylabel("advance per orbit [°]"); axl.legend(fontsize=9)
axr.plot(eps_list, 100 * rel_err, "o", color=C[0], label="measured")
axr.plot(eps_list, 100 * c_fit * eps_list, "--", color=C[2], label=f"{c_fit:.1f}·ε")
axr.set_xlabel("ε"); axr.set_ylabel("relative error of first order [%]"); axr.legend(fontsize=9)
fig.suptitle("B8: first-order perturbation theory fails linearly in ε")
fig.tight_layout(); fig.savefig(FIG / "precession_vs_eps.png", dpi=120); plt.close(fig)

# %%
# --- C1: the isotropic oscillator's three charges and E² = Q1² + Q2² + ω²L_z².
Q1 = sp.Rational(1, 2) * m * (ovx**2 - ovy**2) + sp.Rational(1, 2) * m * w**2 * (ox**2 - oy**2)
Lz_w = m * (ox * ovy - oy * ovx)
E_iso = energy_function(L_iso, [ox, oy], [ovx, ovy])
for Qc in (Q1, Q2, Lz_w, E_iso):
    assert on_shell_dt(Qc, L_iso, [ox, oy], [ovx, ovy], [oax, oay]) == 0
assert sp.simplify(E_iso**2 - (Q1**2 + Q2**2 + w**2 * Lz_w**2)) == 0
dL_Q1 = delta_L(L_iso, [ox, oy], [ovx, ovy], K=[ovx, -ovy], a=[oax, oay])
Lx_part = sp.Rational(1, 2) * m * ovx**2 - sp.Rational(1, 2) * m * w**2 * ox**2
Ly_part = sp.Rational(1, 2) * m * ovy**2 - sp.Rational(1, 2) * m * w**2 * oy**2
assert sp.simplify(dL_Q1 - total_dt(Lx_part - Ly_part, [ox, oy], [ovx, ovy], [oax, oay])) == 0
assert sp.simplify(sp.diff(L_iso, ovx) * ovx - sp.diff(L_iso, ovy) * ovy - (Lx_part - Ly_part) - Q1) == 0
for dL_vel in (dL_Q2, dL_Q1):           # F depends on velocities here, so the F(q, t) constructor must refuse
    try:
        total_derivative_potential(dL_vel, [ox, oy], [ovx, ovy])
        raise AssertionError("a velocity-dependent δL must be refused")
    except ValueError:
        pass
try:                                    # and delta_L itself refuses a velocity-dependent K without accelerations
    delta_L(L_iso, [ox, oy], [ovx, ovy], K=[ovy, ovx])
    raise AssertionError
except ValueError:
    pass
# complex amplitudes: m a_x conj(a_y) = Q2 + i ω L_z
a_x, a_y_conj = ovx + sp.I * w * ox, ovy - sp.I * w * oy          # symbols are real, so conj(a_y) = ẏ − iωy
assert sp.simplify(sp.expand(m * a_x * a_y_conj) - (Q2 + sp.I * w * Lz_w)) == 0
# the trajectory x = cos t, y = 2 sin t (m = ω = 1)
tt = np.linspace(0, 2 * np.pi, 7)
xs_, ys_, vxs_, vys_ = np.cos(tt), 2 * np.sin(tt), -np.sin(tt), 2 * np.cos(tt)
E_n = 0.5 * (vxs_**2 + vys_**2) + 0.5 * (xs_**2 + ys_**2)
Q1_n = 0.5 * (vxs_**2 - vys_**2) + 0.5 * (xs_**2 - ys_**2)
Q2_n = vxs_ * vys_ + xs_ * ys_
Lz_n = xs_ * vys_ - ys_ * vxs_
assert np.allclose(E_n, 2.5) and np.allclose(Q1_n, -1.5) and np.allclose(Q2_n, 0) and np.allclose(Lz_n, 2)
assert np.allclose(E_n**2, Q1_n**2 + Q2_n**2 + Lz_n**2)
print(f"C1  x = cos t, y = 2 sin t: E = {E_n[0]:.2f}, Q1 = {Q1_n[0]:.2f}, Q2 = {Q2_n[0]:.2f}, L_z = {Lz_n[0]:.2f};  "
      f"E² = {E_n[0]**2:.4f} = Q1² + Q2² + ω²L_z² = {Q1_n[0]**2 + Q2_n[0]**2 + Lz_n[0]**2:.4f}")

# --- C2: hydrostatic virial for a uniform sphere; negative heat capacity numbers.
rho_u = 1.3;  R_u = 1.0;  G_u = 1.0
M_u = 4 / 3 * np.pi * R_u**3 * rho_u
P_u = lambda r_: 2 * np.pi * G_u / 3 * rho_u**2 * (R_u**2 - r_**2)       # dP/dr = -G m(r) ρ / r²
three_int_P = 3 * quad(lambda r_: P_u(r_) * 4 * np.pi * r_**2, 0, R_u)[0]
V_grav_u = -quad(lambda r_: G_u * (4 / 3 * np.pi * r_**3 * rho_u) * (4 * np.pi * r_**2 * rho_u) / r_, 0, R_u)[0]
assert np.isclose(three_int_P, -V_grav_u) and np.isclose(-V_grav_u, 0.6 * G_u * M_u**2 / R_u)
# hydrostatic equilibrium check: dP/dr + G m(r) ρ / r² = 0
r_chk = np.linspace(0.05, 0.95, 7)
assert np.allclose(np.gradient(P_u(r_chk), r_chk, edge_order=2), -G_u * (4 / 3 * np.pi * r_chk**3 * rho_u) * rho_u / r_chk**2, rtol=1e-2)
T_sun_real = 3 * T_sun                                                     # |V| ≈ 1.8 GM²/R instead of 0.6
GM_E, R_E = 3.986004e14, 6.371e6
r_sat = R_E + 400e3
v_sat = np.sqrt(GM_E / r_sat)
dE_per_km = GM_E / (2 * r_sat**2) * 1e3                                    # J/kg lost when r drops by 1 km
dv_per_km = 0.5 * v_sat / r_sat * 1e3                                      # m/s gained
assert np.isclose(np.sqrt(GM_E / (r_sat - 1e3)) - v_sat, dv_per_km, rtol=1e-3)
print(f"C2  uniform sphere: 3∫P dV = {three_int_P:.5f} = −V_grav = (3/5)GM²/R = {0.6 * G_u * M_u**2 / R_u:.5f}")
print(f"    Sun with |V| = 1.8 GM²/R: T ≈ {T_sun_real:.2e} K;  satellite at 400 km: v = {v_sat:.1f} m/s, "
      f"loses {dE_per_km:.0f} J/kg and gains {dv_per_km:.3f} m/s per km of descent")

# --- C3: V = −c/r²: fall to the centre (ℓ² < 2mc) or escape (ℓ² > 2mc), compared with Binet.
def inv_sq_rhs(t_, s, c=1.0, m_=1.0):
    """V = −c/r² in the plane: r̈ = −(2c/(m r⁴)) r (force −2c/r³ toward the centre)."""
    x_, y_, vx_, vy_ = s
    r2 = x_ * x_ + y_ * y_
    f = -2 * c / (m_ * r2 * r2)
    return [vx_, vy_, f * x_, f * y_]


def run_inv_sq(ell2, r_stop, t_max=200.0):
    """Start at r = 1, ṟ = 0, angular momentum ℓ (m = c = 1); stop when r = r_stop."""
    ev = lambda t_, s: np.hypot(s[0], s[1]) - r_stop
    ev.terminal = True
    sol_ = solve_ivp(inv_sq_rhs, (0, t_max), [1.0, 0.0, 0.0, np.sqrt(ell2)], method="DOP853",
                     rtol=1e-11, atol=1e-13, events=ev, dense_output=True)
    assert sol_.t_events[0].size == 1, (ell2, r_stop)
    t_hit = sol_.t_events[0][0]
    ts = np.linspace(0, t_hit, 40001)
    xs, ys = sol_.sol(ts)[:2]
    return t_hit, np.unwrap(np.arctan2(ys, xs))[-1], xs, ys


ell2_in, ell2_out = 1.9, 2.1
t_in, th_in, x_in, y_in = run_inv_sq(ell2_in, 0.05)
lam_in = np.sqrt(2 / ell2_in - 1)
t_in_pred, th_in_pred = np.sqrt(1 - 0.05**2) / np.sqrt(2 - ell2_in), np.arccosh(1 / 0.05) / lam_in
assert np.isclose(t_in, t_in_pred, rtol=1e-6) and np.isclose(th_in, th_in_pred, rtol=1e-5), (t_in, t_in_pred, th_in, th_in_pred)
t_out, th_out, x_out, y_out = run_inv_sq(ell2_out, 20.0)
lam_out = np.sqrt(1 - 2 / ell2_out)
t_out_pred, th_out_pred = np.sqrt(20.0**2 - 1) / np.sqrt(ell2_out - 2), np.arccos(1 / 20.0) / lam_out
th_inf = np.pi / (2 * lam_out)
assert np.isclose(t_out, t_out_pred, rtol=1e-6) and np.isclose(th_out, th_out_pred, rtol=1e-5)
sol_c = solve_ivp(inv_sq_rhs, (0, 20.0), [1.0, 0.0, 0.0, np.sqrt(2.0)], method="DOP853", rtol=1e-11, atol=1e-13, t_eval=np.linspace(0, 20, 500))
assert np.abs(np.hypot(sol_c.y[0], sol_c.y[1]) - 1).max() < 1e-8          # ℓ² = 2mc: circular at any radius
print(f"C3  ℓ² = {ell2_in}: reaches r = 0.05 at t = {t_in:.4f} (predicted {t_in_pred:.4f}) after θ = {th_in:.3f} rad "
      f"= {th_in / (2 * np.pi):.2f} turns (Binet: arccosh(20)/λ = {th_in_pred:.3f})")
print(f"    ℓ² = {ell2_out}: reaches r = 20 at t = {t_out:.3f} (predicted {t_out_pred:.3f}) after θ = {th_out:.4f} rad "
      f"(Binet {th_out_pred:.4f}); r → ∞ at θ∞ = π/(2λ') = {th_inf:.4f} rad")
print(f"    λ = √(2mc/ℓ² − 1) = {lam_in:.6f} (ℓ² = {ell2_in}),  λ' = √(1 − 2mc/ℓ²) = {lam_out:.6f} (ℓ² = {ell2_out}), "
      f"arccosh 20 = {np.arccosh(20):.6f}, arccos 0.05 = {np.arccos(0.05):.6f}")
fig, (axi, axo2) = plt.subplots(1, 2, figsize=(10, 4.6))
axi.plot(x_in, y_in, color=C[0], lw=0.8); axi.plot(0, 0, "k*"); axi.set_aspect("equal")
axi.set_title(f"V = −c/r², ℓ² = {ell2_in} < 2mc: spiral into the centre")
axo2.plot(x_out, y_out, color=C[1], lw=0.8); axo2.plot(0, 0, "k*"); axo2.set_aspect("equal")
axo2.set_title(f"ℓ² = {ell2_out} > 2mc: escape after {th_inf / (2 * np.pi):.2f} turns")
fig.tight_layout(); fig.savefig(FIG / "inverse_square_orbits.png", dpi=120); plt.close(fig)

# --- C4: damped oscillator L = e^{γt}(½mẋ² − ½kx²): h not conserved, but Q is.
gam = sp.Symbol("gamma", positive=True)
L_damp = sp.exp(gam * t) * (sp.Rational(1, 2) * m * xv**2 - sp.Rational(1, 2) * k_s * xq**2)
assert sp.simplify(el_accelerations(L_damp, [xq], [xv], [xa])[xa] + gam * xv + k_s * xq / m) == 0
h_damp = energy_function(L_damp, [xq], [xv])
assert sp.simplify(on_shell_dt(h_damp, L_damp, [xq], [xv], [xa]) + gam * L_damp) == 0          # dh/dt = −∂L/∂t = −γL
dL_damp = delta_L(L_damp, [xq], [xv], K=[xv + gam * xq / 2], a=[xa])
assert sp.simplify(dL_damp - total_dt(L_damp, [xq], [xv], [xa])) == 0                           # F = L
Q_damp = sp.simplify(sp.diff(L_damp, xv) * (xv + gam * xq / 2) - L_damp)
Q_damp_expected = sp.exp(gam * t) * (sp.Rational(1, 2) * m * xv**2 + sp.Rational(1, 2) * k_s * xq**2 + gam * m * xq * xv / 2)
assert sp.simplify(Q_damp - Q_damp_expected) == 0
assert on_shell_dt(Q_damp, L_damp, [xq], [xv], [xa]) == 0
gam_n, Om_n = 0.2, np.sqrt(1 - 0.2**2 / 4)
sol_d = solve_ivp(lambda t_, s: [s[1], -gam_n * s[1] - s[0]], (0, 20), [1.0, 0.0], method="DOP853",
                  rtol=1e-11, atol=1e-13, t_eval=np.linspace(0, 20, 401))
xd, vd = sol_d.y
Qd = np.exp(gam_n * sol_d.t) * (0.5 * vd**2 + 0.5 * xd**2 + 0.5 * gam_n * xd * vd)
hd = np.exp(gam_n * sol_d.t) * (0.5 * vd**2 + 0.5 * xd**2)
assert np.abs(Qd - 0.5).max() < 1e-8 and hd.max() - hd.min() > 0.05
assert np.isclose(hd.min(), 0.5 / 1.1, atol=2e-4) and np.isclose(hd.max(), 0.5 / 0.9, atol=2e-4)   # h ∈ [Qω₀/(ω₀+γ/2), Qω₀/(ω₀−γ/2)]
x20_pred = np.exp(-gam_n * 10) * (np.cos(Om_n * 20) + gam_n / (2 * Om_n) * np.sin(Om_n * 20))
assert np.isclose(xd[-1], x20_pred, atol=1e-9)
print(f"C4  damped oscillator (m = k = 1, γ = 0.2, x₀ = 1, ẋ₀ = 0): Q = {Qd[0]:.6f} at t = 0, {Qd[-1]:.6f} at t = 20 "
      f"(|ΔQ| ≤ {np.abs(Qd - 0.5).max():.1e}); h ranges over [{hd.min():.4f}, {hd.max():.4f}]; x(20) = {xd[-1]:.6f}")
print(f"    Ω = √(1 − γ²/4) = {Om_n:.6f}, γ/(2Ω) = {gam_n / (2 * Om_n):.6f}, e^(−2) = {np.exp(-2):.6f}")

# %%
# --- Further checks behind the written solutions (A1, A2, A3, A4, A6c, B1e, B2e, B3a/e, B4c/d, C1f, C3e, C4c/d).
# A1(a): spherical pendulum, phi cyclic.  A2(a): h = T + V for the polar central force.
ell_ = sp.Symbol("ell", positive=True)
(thp, php), (thpd, phpd), _ = coords(["theta_p", "phi_p"])
L_sph = sp.Rational(1, 2) * m * ell_**2 * (thpd**2 + sp.sin(thp)**2 * phpd**2) + m * g * ell_ * sp.cos(thp)
assert sp.simplify(noether_charge(L_sph, [thp, php], [thpd, phpd], K=[0, 1]) - m * ell_**2 * sp.sin(thp)**2 * phpd) == 0
h_polar = energy_function(L_polar, [r, th], [rdot, thdot])
assert sp.simplify(h_polar - (sp.Rational(1, 2) * m * (rdot**2 + r**2 * thdot**2) + Vf(r))) == 0
assert on_shell_dt(h_polar, L_polar, [r, th], [rdot, thdot], [rddot, thddot]) == 0

# A3: Euler's theorem (3.10) for three homogeneous functions; its failure for x² + y; the T = ½ q̇ᵀM(q) q̇ case.
u1, u2, u3 = sp.symbols("u1 u2 u3", positive=True)


def euler_defect(f, vars_, k):
    """sum_i x_i ∂f/∂x_i − k f, which (3.10) says vanishes when f is homogeneous of degree k."""
    return sp.simplify(sum(xi * sp.diff(f, xi) for xi in vars_) - k * f)


assert euler_defect(-kappa / sp.sqrt(u1**2 + u2**2 + u3**2), [u1, u2, u3], -1) == 0
assert euler_defect(sp.Rational(1, 2) * k_s * (u1**2 + u2**2), [u1, u2], 2) == 0
assert euler_defect(u1**2 * u2 / u3, [u1, u2, u3], 2) == 0
assert euler_defect(u1**2 + u2, [u1, u2], 1) == u1**2 and euler_defect(u1**2 + u2, [u1, u2], 2) == -u2
M11, M12, M22 = [sp.Function(n)(ox, oy) for n in ("M11", "M12", "M22")]
T_gen = sp.Rational(1, 2) * (M11 * ovx**2 + 2 * M12 * ovx * ovy + M22 * ovy**2)
assert sp.simplify(ovx * sp.diff(T_gen, ovx) + ovy * sp.diff(T_gen, ovy) - 2 * T_gen) == 0

# A4: similarity exponents (time, speed, energy, angular momentum) ∝ size^ν, and the ratios for a' = 4a.
def similarity_exponents(k):
    """ν in (time, speed, energy, angular momentum) ∝ size^ν for V homogeneous of degree k (Theorem 3.8)."""
    return {"time": 1 - k / 2, "speed": k / 2, "energy": k, "angmom": 1 + k / 2}


tab = {k: similarity_exponents(k) for k in (-2, -1, 1, 2, 4)}
assert tab[-1] == {"time": 1.5, "speed": -0.5, "energy": -1, "angmom": 0.5}
assert tab[2]["time"] == 0 and tab[4]["time"] == -1 and tab[-2]["angmom"] == 0 and tab[1]["time"] == 0.5
ratios_4 = {key: 4.0**nu for key, nu in tab[-1].items()}
assert np.allclose([ratios_4[k_] for k_ in ("time", "speed", "energy", "angmom")], [8, 0.5, 0.25, 2])
print("A4  k = -1, a' = 4a: period ×%.0f, speed ×%.2f, energy ×%.2f, angular momentum ×%.0f"
      % (ratios_4["time"], ratios_4["speed"], ratios_4["energy"], ratios_4["angmom"]))

# A6(c): V(x, t) = V0(x) + f(t) is the only time dependence that keeps K = ẋ a quasi-symmetry; Q = ½mẋ² + V0.
ff = sp.Function("f")
L_sep = sp.Rational(1, 2) * m * xv**2 - Vf(xq) - ff(t)
dL_sep = delta_L(L_sep, [xq], [xv], K=[xv], a=[xa])
assert sp.simplify(dL_sep - (total_dt(L_sep, [xq], [xv], [xa]) - sp.diff(L_sep, t))) == 0   # δL = dL/dt − ∂L/∂t
Q_sep = sp.simplify(sp.diff(L_sep, xv) * xv - (L_sep + ff(t)))                               # F = L + f
assert sp.simplify(Q_sep - (sp.Rational(1, 2) * m * xv**2 + Vf(xq))) == 0
assert on_shell_dt(Q_sep, L_sep, [xq], [xv], [xa]) == 0

# B1: h for uniform gravity; x-boost; the Galilean count; the deformed rotation that survives.
h_grav = energy_function(L_grav, [x, y, z], [vx, vy, vz])
assert sp.simplify(h_grav - (sp.Rational(1, 2) * m * (vx**2 + vy**2 + vz**2) + m * g * z)) == 0
assert on_shell_dt(h_grav, L_grav, [x, y, z], [vx, vy, vz], [ax, ay, az]) == 0
Q_bx, F_bx = noether_charge_quasi(L_grav, [x, y, z], [vx, vy, vz], K=[t, 0, 0])
assert sp.simplify(F_bx - m * x) == 0 and sp.simplify(Q_bx - (m * vx * t - m * x)) == 0
gens = {"Px": [1, 0, 0], "Py": [0, 1, 0], "Pz": [0, 0, 1], "Lx": [0, -z, y], "Ly": [z, 0, -x], "Lz": [-y, x, 0],
        "Gx": [t, 0, 0], "Gy": [0, t, 0], "Gz": [0, 0, t]}
surviving = []
for name, Kg in gens.items():
    try:
        noether_charge_quasi(L_grav, [x, y, z], [vx, vy, vz], Kg)
        surviving.append(name)
    except ValueError:
        pass
assert surviving == ["Px", "Py", "Pz", "Lz", "Gx", "Gy", "Gz"]        # plus h: eight of the ten
K_def = [0, -z - g * t**2 / 2, y]                                        # rotation about x̂ in the freely falling frame
Q_def, F_def = noether_charge_quasi(L_grav, [x, y, z], [vx, vy, vz], K=K_def)
assert sp.simplify(F_def + m * g * t * y) == 0
assert sp.simplify(Q_def - (m * (y * vz - z * vy) + m * g * (t * y - t**2 * vy / 2))) == 0
assert on_shell_dt(Q_def, L_grav, [x, y, z], [vx, vy, vz], [ax, ay, az]) == 0
Q_defy, F_defy = noether_charge_quasi(L_grav, [x, y, z], [vx, vy, vz], K=[z + g * t**2 / 2, 0, -x])   # about ŷ
assert sp.simplify(F_defy - m * g * t * x) == 0
assert sp.simplify(Q_defy - (m * (z * vx - x * vz) + m * g * (t**2 * vx / 2 - t * x))) == 0
assert on_shell_dt(Q_defy, L_grav, [x, y, z], [vx, vy, vz], [ax, ay, az]) == 0
print("B1  gravity: surviving Galilean generators", surviving, "+ h  (8 of 10);")
print("    deformed x-rotation K = (0, −z − gt²/2, y): F =", F_def, ", Q =", sp.factor(Q_def))

# B2(e): in (X, s) coordinates E_rel = ½μṡ² + V(s) and E_cm = ½MẊ² are conserved separately.
(Xc_, sc_), (Xcd, scd), (Xcdd, scdd) = coords(["X", "s"])
L_cm = sp.Rational(1, 2) * Mtot * Xcd**2 + sp.Rational(1, 2) * mu * scd**2 - Vf(sc_)
L_from_2body = L_2body.subs({x1: Xc_ + m2 / Mtot * sc_, x2: Xc_ - m1 / Mtot * sc_,
                             v1: Xcd + m2 / Mtot * scd, v2: Xcd - m1 / Mtot * scd}, simultaneous=True)
assert sp.simplify(L_from_2body - L_cm) == 0
E_rel = sp.Rational(1, 2) * mu * scd**2 + Vf(sc_)
assert on_shell_dt(E_rel, L_cm, [Xc_, sc_], [Xcd, scd], [Xcdd, scdd]) == 0
assert on_shell_dt(sp.Rational(1, 2) * Mtot * Xcd**2, L_cm, [Xc_, sc_], [Xcd, scd], [Xcdd, scdd]) == 0
assert sp.simplify(energy_function(L_cm, [Xc_, sc_], [Xcd, scd]) - (sp.Rational(1, 2) * Mtot * Xcd**2 + E_rel)) == 0

# B3(a): ẍ = ω_c ẏ, ÿ = −ω_c ẋ.  B3(e): Q_rot = (p_y² + Q_x²)/(2q_eB) − m h/(q_eB) is not independent.
acc_mag = el_accelerations(L_mag, [x, y], [vx, vy], [ax, ay])
wc = qe * Bf / m
assert sp.simplify(acc_mag[ax] - wc * vy) == 0 and sp.simplify(acc_mag[ay] + wc * vx) == 0
h_mag = energy_function(L_mag, [x, y], [vx, vy])
assert sp.simplify(Qrot_mag - ((py_mag**2 + Qx_mag**2) / (2 * qe * Bf) - m * h_mag / (qe * Bf))) == 0

# B4(c): the isotropic orbit x = A cos ωt, y = B cos(ωt − φ) is a centred ellipse with L_z = mωAB sin φ.
A_, B_amp, phi_ = sp.symbols("A B_amp phi", real=True)
xt, yt = A_ * sp.cos(w * t), B_amp * sp.cos(w * t - phi_)
assert sp.simplify((xt / A_)**2 + (yt / B_amp)**2 - 2 * xt * yt * sp.cos(phi_) / (A_ * B_amp) - sp.sin(phi_)**2) == 0
assert sp.simplify(m * (xt * sp.diff(yt, t) - yt * sp.diff(xt, t)) - m * w * A_ * B_amp * sp.sin(phi_)) == 0
# B4(d): a_x² conj(a_y) = Q3 + i Q4 is constant when ω2 = 2ω1.
a_x1, a_y2c = ovx + sp.I * w1 * ox, ovy - 2 * sp.I * w1 * oy
assert sp.simplify(on_shell_dt(a_x1, L_12, [ox, oy], [ovx, ovy], [oax, oay]) - sp.I * w1 * a_x1) == 0
Q4 = 2 * w1 * (ox * ovx * ovy - oy * (ovx**2 - w1**2 * ox**2))
assert sp.simplify(sp.expand(a_x1**2 * a_y2c) - (Q3 + sp.I * Q4)) == 0
assert on_shell_dt(Q4, L_12, [ox, oy], [ovx, ovy], [oax, oay]) == 0

# C1(f): a rotation of the plane by α rotates (Q1, Q2) by 2α and fixes L_z.
alpha_ = sp.Symbol("alpha", real=True)
rot = {ox: sp.cos(alpha_) * ox - sp.sin(alpha_) * oy, oy: sp.sin(alpha_) * ox + sp.cos(alpha_) * oy,
       ovx: sp.cos(alpha_) * ovx - sp.sin(alpha_) * ovy, ovy: sp.sin(alpha_) * ovx + sp.cos(alpha_) * ovy}
Q1r, Q2r, Lzr = [sp.expand(Qc.subs(rot, simultaneous=True)) for Qc in (Q1, Q2, Lz_w)]
assert sp.simplify(Q1r - (sp.cos(2 * alpha_) * Q1 - sp.sin(2 * alpha_) * Q2)) == 0
assert sp.simplify(Q2r - (sp.sin(2 * alpha_) * Q1 + sp.cos(2 * alpha_) * Q2)) == 0
assert sp.simplify(Lzr - Lz_w) == 0
print("C1  rotation by α: (Q1, Q2) → rotated by 2α, L_z fixed  (verified symbolically)")

# C3(e): no length and no time can be built from m, c, ℓ; ℓ²/(mc) is dimensionless.
Dm = np.array([[1, 1, 1], [0, 4, 2], [0, -2, -1]], float)       # rows M, L, T; columns m, c, ℓ
for target in (np.array([0, 1, 0.]), np.array([0, 0, 1.])):     # a length, a time
    assert np.linalg.matrix_rank(np.column_stack([Dm, target])) > np.linalg.matrix_rank(Dm)
assert np.allclose(Dm @ np.array([-1, -1, 2]), 0)
print("C3  dimensions: rank", np.linalg.matrix_rank(Dm), "of 3; no length or time from (m, c, ℓ); ℓ²/(mc) dimensionless")

# C4(c): time shift alone and rescaling alone each change L by a multiple of L; only their sum is a total derivative.
dL_shift = delta_L(L_damp, [xq], [xv], K=[xv], a=[xa])
assert sp.simplify(dL_shift - (total_dt(L_damp, [xq], [xv], [xa]) - gam * L_damp)) == 0
assert sp.simplify(delta_L(L_damp, [xq], [xv], K=[gam * xq / 2]) - gam * L_damp) == 0
# C4(d): x = e^{−γt/2} X turns Q into ½mẊ² + ½mΩ²X², Ω² = k/m − γ²/4.
XX, XXd = sp.symbols("X_ Xdot_", real=True)
Q_X = sp.simplify(Q_damp.subs({xq: sp.exp(-gam * t / 2) * XX, xv: sp.exp(-gam * t / 2) * (XXd - gam * XX / 2)},
                              simultaneous=True))
assert sp.simplify(Q_X - (sp.Rational(1, 2) * m * XXd**2 + sp.Rational(1, 2) * m * (k_s / m - gam**2 / 4) * XX**2)) == 0
# C4(e): the h range [Qω₀/(ω₀ + γ/2), Qω₀/(ω₀ − γ/2)] in closed form, compared with the integration above.
h_lo, h_hi = 0.5 * 1.0 / (1.0 + gam_n / 2), 0.5 * 1.0 / (1.0 - gam_n / 2)
assert np.isclose(hd.min(), h_lo, atol=2e-4) and np.isclose(hd.max(), h_hi, atol=2e-4)
print(f"C4  Q = ½mẊ² + ½mΩ²X² verified; h ∈ [{h_lo:.4f}, {h_hi:.4f}] in closed form")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
