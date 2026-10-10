# %% [markdown]
# # CM-7 Project: Rigid Bodies — the Inertia Tensor, Euler's Equations, Euler Angles, the Top
#
# This script builds every object in the notes and checks every number quoted in
# the worked examples and solutions. Parts: (1) the inertia tensor from a point
# cloud via `einsum` and its symmetry/positivity; (2) Monte-Carlo bodies (cube,
# plate, cone, sphere) against the closed forms, the parallel- and
# perpendicular-axis theorems; (3) the tensor law I' = R I R^T and principal axes,
# with the spinning-plate bearing torque; (4) angular velocity extracted from
# R-dot R^T; (5) Euler's equations integrated torque-free: stable/stable/unstable
# with the flip, Poinsot polhodes, the symmetric top's two precession rates, the
# frisbee and the Chandler wobble; (6) Euler angles: R(phi, theta, psi), the
# omega formula checked by finite differences, gimbal-lock conditioning; (7) the
# heavy symmetric top from its Lagrangian with SymPy: V_eff, steady precession,
# nutation, the sleeping-top threshold, the axis path on the unit sphere;
# (8) precession of the equinoxes.
#
# Conventions (same as the notes): body axes e_1, e_2, e_3 with body components
# subscripted 1,2,3; space components x,y,z; r_space = R r_body, columns of R are
# the body axes in space coordinates; omega is extracted from R-dot R^T (space)
# or R^T R-dot (body); Euler angles (phi, theta, psi) in the z-x'-z'' sequence,
# R = R_z(phi) R_x(theta) R_z(psi).
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp
from pathlib import Path

T0 = time.time()
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)


def hat(v):
    """The antisymmetric matrix [v]_x with [v]_x u = v x u, eq. (3.3) of the notes."""
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def unhat(K):
    """Inverse of hat: the vector v with K = [v]_x. Averages the two off-diagonal copies."""
    return 0.5 * np.array([K[2, 1] - K[1, 2], K[0, 2] - K[2, 0], K[1, 0] - K[0, 1]])


def rodrigues(n, theta):
    """Rotation by theta about unit axis n, R = I + sin(theta) K + (1 - cos(theta)) K^2 (LA-12 (A.4))."""
    n = np.asarray(n, float) / np.linalg.norm(n)
    K = hat(n)
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * K @ K


def R_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def R_x(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


# %% [markdown]
# ## Part 1: The inertia tensor from a point cloud (§3.3, eq. (3.8))
#
# I_ij = sum_a m_a (|r_a|^2 delta_ij - x_ai x_aj). The einsum string 'a,ai,aj->ij' is
# literally the sum over a of m_a x_ai x_aj (LA-12 Example 4.6). We check: I is
# symmetric, L = I omega agrees with sum m r x (omega x r), T = omega.I.omega / 2
# agrees with sum m |omega x r|^2 / 2, I is positive definite for a non-collinear
# cloud, and the triangle inequalities I_1 + I_2 >= I_3 hold.
#
# **Build it yourself first:** write `inertia_tensor` with one einsum call per term.


# %%
def inertia_tensor(masses, positions):
    """Inertia tensor I_ij = sum_a m_a (r_a^2 delta_ij - x_ai x_aj), notes eq. (3.8), about the origin."""
    m = np.asarray(masses, float)
    X = np.asarray(positions, float)
    r2 = np.einsum("ai,ai->a", X, X)
    return np.einsum("a,a->", m, r2) * np.eye(3) - np.einsum("a,ai,aj->ij", m, X, X)


def angular_momentum_direct(masses, positions, omega):
    """L = sum_a m_a r_a x (omega x r_a), the definition before BAC-CAB (eq. (3.7))."""
    L = np.zeros(3)
    for m, r in zip(masses, positions):
        L += m * np.cross(r, np.cross(omega, r))
    return L


def kinetic_energy_direct(masses, positions, omega):
    """T = (1/2) sum_a m_a |omega x r_a|^2, eq. (3.5)."""
    return 0.5 * sum(m * np.dot(np.cross(omega, r), np.cross(omega, r)) for m, r in zip(masses, positions))


# LA-12 Example 4.6: masses (1,2,1) at (1,1,0), (-1,0,1), (0,1,1): principal moments 3, 5, 8.
m3 = np.array([1.0, 2.0, 1.0])
X3 = np.array([[1, 1, 0], [-1, 0, 1], [0, 1, 1]], float)
I3pt = inertia_tensor(m3, X3)
assert np.allclose(I3pt, [[5, -1, 2], [-1, 6, -1], [2, -1, 5]])
assert np.allclose(np.linalg.eigvalsh(I3pt), [3, 5, 8])
print("Part 1: LA-12 Example 4.6 reproduced; principal moments", np.linalg.eigvalsh(I3pt))

for _ in range(20):
    n = rng.integers(3, 12)
    m = rng.uniform(0.5, 2.0, n)
    X = rng.standard_normal((n, 3))
    I = inertia_tensor(m, X)
    w = rng.standard_normal(3)
    assert np.allclose(I, I.T)                                              # symmetric
    assert np.allclose(I @ w, angular_momentum_direct(m, X, w))             # L = I omega
    assert np.isclose(0.5 * w @ I @ w, kinetic_energy_direct(m, X, w))      # T = omega.I.omega/2
    lam = np.linalg.eigvalsh(I)
    assert lam.min() > 0                                                    # positive definite
    assert lam[0] + lam[1] >= lam[2] - 1e-12                                # triangle inequality (3.12)

# Positive semidefinite but singular: all mass on a line (the rod has I_3 = 0).
Xline = np.outer(np.linspace(-1, 1, 7), [1.0, 2.0, 2.0])
lam_line = np.linalg.eigvalsh(inertia_tensor(np.ones(7), Xline))
assert np.isclose(lam_line[0], 0) and lam_line[1] > 0
print("Part 1 checks passed: symmetric, PD (singular only for a line of masses), L = I w, T = w.I.w/2.")

# %% [markdown]
# ## Part 2: Monte-Carlo bodies against the closed forms (§3.6–3.7)
#
# Sample N uniform points in each body, build I/M from the sample, compare with the
# table in §3.7. Parallel-axis and perpendicular-axis theorems are exact identities
# on a sample, so they are asserted tightly; the closed forms are statistical and are
# asserted to about 1.5 %.


# %%
N = 200_000


def sample_cube(a, n):
    return rng.uniform(-a / 2, a / 2, (n, 3))


def sample_plate(a, b, n):
    P = np.zeros((n, 3))
    P[:, 0] = rng.uniform(-a / 2, a / 2, n)
    P[:, 1] = rng.uniform(-b / 2, b / 2, n)
    return P


def sample_sphere(R, n):
    v = rng.standard_normal((n, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    r = R * rng.uniform(0, 1, n) ** (1 / 3)
    return v * r[:, None]


def sample_cone(R, h, n):
    """Uniform in a solid cone with apex at the origin, axis +z, base radius R at z = h."""
    z = h * rng.uniform(0, 1, n) ** (1 / 3)        # pdf proportional to z^2 (disc area)
    rho = (R * z / h) * np.sqrt(rng.uniform(0, 1, n))
    phi = rng.uniform(0, 2 * np.pi, n)
    return np.c_[rho * np.cos(phi), rho * np.sin(phi), z]


def I_per_mass(P):
    return inertia_tensor(np.full(len(P), 1.0 / len(P)), P)


def parallel_axis(I_cm, M, d):
    """I about a point displaced by d from the CM: I_cm + M(d^2 delta - d d^T), eq. (3.13)."""
    d = np.asarray(d, float)
    return I_cm + M * (d @ d * np.eye(3) - np.outer(d, d))


a = 1.0
cube = sample_cube(a, N)
I_cube = I_per_mass(cube)
assert np.allclose(np.diag(I_cube), a**2 / 6, rtol=0.015) and np.allclose(I_cube - np.diag(np.diag(I_cube)), 0, atol=5e-3)

# Cube about a corner (Example 4.1): closed form (M a^2)[(11/12) delta - (1/4) J], eigenvalues 1/6, 11/12, 11/12.
I_corner_closed = a**2 * ((11 / 12) * np.eye(3) - 0.25 * np.ones((3, 3)))
assert np.allclose(np.linalg.eigvalsh(I_corner_closed), [1 / 6, 11 / 12, 11 / 12])
lam_c, Q_c = np.linalg.eigh(I_corner_closed)
assert np.isclose(abs(Q_c[:, 0] @ (np.ones(3) / np.sqrt(3))), 1.0)           # smallest moment along the diagonal
# parallel-axis theorem applied to the SAME sample is exact, provided I_cm is about the sample's own CM.
c_cube = cube.mean(axis=0)
I_cube_cm = I_per_mass(cube - c_cube)
I_corner_sample = I_per_mass(cube + a / 2)                                   # corner at (-a/2,-a/2,-a/2) moved to origin
assert np.allclose(I_corner_sample, parallel_axis(I_cube_cm, 1.0, c_cube + a / 2), atol=1e-12)
assert np.allclose(I_corner_sample, I_corner_closed, atol=0.01)
print("Part 2: cube about corner eigenvalues (M a^2 units):", np.round(np.linalg.eigvalsh(I_corner_sample), 3), "closed form 1/6, 11/12, 11/12")

# Rectangular plate a x b: (b^2/12, a^2/12, (a^2+b^2)/12); perpendicular-axis exact on the sample.
pa, pb = 0.3, 0.1
plate = sample_plate(pa, pb, N)
I_plate = I_per_mass(plate)
assert np.allclose(np.diag(I_plate), [pb**2 / 12, pa**2 / 12, (pa**2 + pb**2) / 12], rtol=0.015)
assert np.isclose(I_plate[2, 2], I_plate[0, 0] + I_plate[1, 1])                                    # perpendicular-axis theorem (3.14)

# Sphere: 2/5 R^2 on the diagonal.
Rs = 1.0
I_sph = I_per_mass(sample_sphere(Rs, N))
assert np.allclose(np.diag(I_sph), 0.4 * Rs**2, rtol=0.015)

# Cone about its apex: I_3 = 3/10 R^2, I_1 = I_2 = 3/20 R^2 + 3/5 h^2; CM at 3h/4 so I_1^cm = 3/20 R^2 + 3/80 h^2.
Rc, hc = 1.0, 2.0
cone = sample_cone(Rc, hc, N)
I_cone_apex = I_per_mass(cone)
assert np.isclose(I_cone_apex[2, 2], 0.3 * Rc**2, rtol=0.015)
assert np.allclose([I_cone_apex[0, 0], I_cone_apex[1, 1]], 0.15 * Rc**2 + 0.6 * hc**2, rtol=0.015)
cm_cone = cone.mean(axis=0)
assert np.allclose(cm_cone, [0, 0, 0.75 * hc], atol=0.01)
I_cone_cm = I_per_mass(cone - cm_cone)
assert np.allclose(I_cone_apex, parallel_axis(I_cone_cm, 1.0, -cm_cone), atol=1e-12)       # exact on the sample
assert np.allclose([I_cone_cm[0, 0], I_cone_cm[1, 1]], 0.15 * Rc**2 + (3 / 80) * hc**2, rtol=0.02)
print(f"Part 2 checks passed: cube, plate, sphere, cone within 1.5% of closed forms (N = {N}); parallel/perpendicular-axis exact.")

# %% [markdown]
# ## Part 3: I is a (0,2) tensor: I' = R I R^T; principal axes; the spinning plate (§3.4–3.5, Example 4.2)
#
# Rotating every point of the cloud by R rotates I by R I R^T (eq. (3.9)); the
# eigenvalues do not move. `eigh` returns the principal axes; spinning about one gives
# L parallel to omega; spinning a plate about its diagonal does not, and the bearing
# must supply the torque omega x L.


# %%
for _ in range(10):
    n_ax = rng.standard_normal(3)
    R = rodrigues(n_ax, rng.uniform(0, np.pi))
    assert np.allclose(R.T @ R, np.eye(3)) and np.isclose(np.linalg.det(R), 1)
    mm = rng.uniform(0.5, 2, 9)
    XX = rng.standard_normal((9, 3))
    I_old = inertia_tensor(mm, XX)
    I_new = inertia_tensor(mm, XX @ R.T)                 # each row r -> R r
    assert np.allclose(I_new, R @ I_old @ R.T)
    assert np.allclose(np.linalg.eigvalsh(I_new), np.linalg.eigvalsh(I_old))
    lam, Q = np.linalg.eigh(I_old)
    assert np.allclose(Q.T @ I_old @ Q, np.diag(lam), atol=1e-10)
    for k in range(3):                                   # L parallel to omega along a principal axis
        Lk = I_old @ Q[:, k]
        assert np.allclose(Lk, lam[k] * Q[:, k])

# Example 4.2: rectangular plate M = 2 kg, a = 0.3 m, b = 0.1 m spun at 50 Hz about its diagonal.
Mp, ap, bp = 2.0, 0.3, 0.1
I_plate_closed = Mp / 12 * np.diag([bp**2, ap**2, ap**2 + bp**2])
w_mag = 2 * np.pi * 50
n_diag = np.array([ap, bp, 0.0]) / np.hypot(ap, bp)
w_plate = w_mag * n_diag
L_plate = I_plate_closed @ w_plate
cos_ang = L_plate @ w_plate / (np.linalg.norm(L_plate) * w_mag)
ang_LW = np.degrees(np.arccos(cos_ang))
tau_bearing = np.cross(w_plate, L_plate)                 # dL/dt = omega x L for L fixed in the body
assert np.isclose(cos_ang, 2 * ap * bp / (ap**2 + bp**2))                     # = 0.6
assert np.isclose(ang_LW, 53.130, atol=2e-3)
assert np.isclose(np.linalg.norm(L_plate), Mp * ap * bp * w_mag / 12)         # = 1.5708
assert np.isclose(np.linalg.norm(tau_bearing), w_mag**2 * Mp * ap * bp / 12 * (ap**2 - bp**2) / (ap**2 + bp**2))
assert np.isclose(np.linalg.norm(tau_bearing), 394.78, rtol=1e-4)
assert np.isclose(0.5 * w_plate @ I_plate_closed @ w_plate, 148.04, rtol=1e-4)  # kinetic energy in J
print(f"Example 4.2: angle(L, w) = {ang_LW:.2f} deg, |L| = {np.linalg.norm(L_plate):.4f} kg m^2/s, "
      f"bearing torque = {np.linalg.norm(tau_bearing):.1f} N m, T = {0.5 * w_plate @ I_plate_closed @ w_plate:.2f} J")
print("Part 3 checks passed: I' = R I R^T, principal axes from eigh, spinning plate numbers.")

# %% [markdown]
# ## Part 4: Angular velocity from R-dot R^T (§3.2, Theorem 3.1)
#
# For R(t) orthogonal, R-dot R^T is antisymmetric and equals [omega]_x. We check
# this by finite differences on a rotation about a fixed axis (omega = Omega n) and
# on a rotation whose axis itself moves. Finite rotations do not commute;
# infinitesimal ones do to first order.


# %%
def omega_from_R(R_func, t, h=1e-6):
    """omega (space components) from the central difference of R(t): [omega]_x = R-dot R^T."""
    Rd = (R_func(t + h) - R_func(t - h)) / (2 * h)
    Om = Rd @ R_func(t).T
    assert np.allclose(Om, -Om.T, atol=1e-6), "R-dot R^T must be antisymmetric"
    return unhat(Om)


n_fixed = np.array([1.0, 2.0, 2.0]) / 3
Om0 = 1.7
w_fd = omega_from_R(lambda t: rodrigues(n_fixed, Om0 * t), 0.83)
assert np.allclose(w_fd, Om0 * n_fixed, atol=1e-6)

# Moving axis: R(t) = R_z(0.5 t) R_x(0.3 t). Space omega = 0.5 z + 0.3 R_z(0.5t) x  (each factor contributes its own axis).
t_s = 0.6
w_fd2 = omega_from_R(lambda t: R_z(0.5 * t) @ R_x(0.3 * t), t_s)
assert np.allclose(w_fd2, 0.5 * np.array([0, 0, 1.0]) + 0.3 * R_z(0.5 * t_s) @ np.array([1.0, 0, 0]), atol=1e-6)

# Body-frame omega is R^T R-dot: for a body-fixed point r_b, r-dot_space = omega_space x r_space.
R_now = R_z(0.5 * t_s) @ R_x(0.3 * t_s)
rb = np.array([0.2, -0.4, 0.9])
h = 1e-6
rdot = ((R_z(0.5 * (t_s + h)) @ R_x(0.3 * (t_s + h))) @ rb - (R_z(0.5 * (t_s - h)) @ R_x(0.3 * (t_s - h))) @ rb) / (2 * h)
assert np.allclose(rdot, np.cross(w_fd2, R_now @ rb), atol=1e-6)
w_body = R_now.T @ w_fd2
assert np.allclose(hat(w_body), R_now.T @ ((R_z(0.5 * (t_s + h)) @ R_x(0.3 * (t_s + h)) - R_z(0.5 * (t_s - h)) @ R_x(0.3 * (t_s - h))) / (2 * h)), atol=1e-6)

# Finite rotations do not commute (LA-1 Example 4.3); infinitesimal ones do to first order.
A, B = rodrigues([1, 0, 0], np.pi / 2), rodrigues([0, 1, 0], np.pi / 2)
assert not np.allclose(A @ B, B @ A)
eps = 1e-4
Ae, Be = rodrigues([1, 0, 0], eps), rodrigues([0, 1, 0], eps)
assert np.allclose(Ae @ Be, Be @ Ae, atol=2 * eps**2)                         # commutator is O(eps^2)
assert np.allclose(Ae @ Be, np.eye(3) + eps * (hat([1, 0, 0]) + hat([0, 1, 0])), atol=2 * eps**2)
print("Part 4 checks passed: R-dot R^T antisymmetric, omega recovered, r-dot = omega x r, infinitesimal rotations add.")

# %% [markdown]
# ## Part 5: Euler's equations, torque-free (§3.8–3.10, Examples 4.3, 4.4)
#
# I_1 w1-dot = (I_2 - I_3) w2 w3 and cyclic. For I = (1, 2, 3): spin near axis 1 and
# axis 3 is stable (small oscillation), near axis 2 unstable with growth rate
# sigma = w0 sqrt((I_2 - I_1)(I_3 - I_2)/(I_1 I_3)) = w0/sqrt(3). T and |L| are conserved.
# The symmetric top (I_1 = I_2) precesses in the body frame at Omega_b = w3 (I_3 - I_1)/I_1
# and its axis precesses in space at Omega_s = L/I_1.


# %%
def euler_rhs(t, w, I):
    """Euler's equations (3.17) in the principal body frame, torque-free."""
    I1, I2, I3 = I
    w1, w2, w3 = w
    return [(I2 - I3) * w2 * w3 / I1, (I3 - I1) * w3 * w1 / I2, (I1 - I2) * w1 * w2 / I3]


def euler_with_attitude_rhs(t, y, I):
    """Euler's equations plus R-dot = R [omega_body]_x, so the attitude R (columns = body axes in space) is carried along."""
    w = y[:3]
    R = y[3:].reshape(3, 3)
    dw = euler_rhs(t, w, I)
    dR = R @ hat(w)
    return np.concatenate([dw, dR.ravel()])


I123 = np.array([1.0, 2.0, 3.0])
w0 = 1.0
delta = 0.01
runs = {}
for axis in range(3):
    w_init = np.full(3, delta)
    w_init[axis] = w0
    t_end = 60.0
    sol = solve_ivp(euler_rhs, (0, t_end), w_init, args=(I123,), method="DOP853", rtol=1e-10, atol=1e-12,
                    t_eval=np.linspace(0, t_end, 6001))
    W = sol.y
    T_t = 0.5 * np.einsum("i,it->t", I123, W**2)
    L2_t = np.einsum("i,it->t", I123**2, W**2)
    assert np.allclose(T_t, T_t[0], rtol=1e-8) and np.allclose(L2_t, L2_t[0], rtol=1e-8)
    runs[axis] = sol

# Stable axes: the deviation stays O(delta); unstable axis: w2 flips sign (passes through zero).
for axis in (0, 2):
    W = runs[axis].y
    others = [k for k in range(3) if k != axis]
    assert np.abs(W[others]).max() < 5 * delta
    assert W[axis].min() > 0.99 * w0
W_mid = runs[1].y
assert W_mid[1].min() < -0.99 * w0, "spin about the intermediate axis must flip"

# Growth rate about axis 2, Example 4.4: sigma = w0 sqrt((I2-I1)(I3-I2)/(I1 I3)) = w0/sqrt(3).
# The linearization (3.18) has a growing mode (dw1, dw3) ~ (1, -sigma I1/((I2-I3) w0)) ... for I = (1,2,3): dw3 = -sigma dw1.
# A generic perturbation also excites the decaying mode, so |w1| first dips through zero (visible in the figure);
# to measure sigma cleanly we start on the growing eigenvector with a tiny amplitude.
sigma_theory = w0 * np.sqrt((I123[1] - I123[0]) * (I123[2] - I123[1]) / (I123[0] * I123[2]))
assert np.isclose(sigma_theory, 1 / np.sqrt(3))
amp = 1e-4
w_grow = np.array([amp, w0, -sigma_theory * amp * I123[0] / ((I123[1] - I123[2]) * w0) * (-1)])   # dw3 = (I1 sigma/((I2-I3) w0)) dw1 ... sign worked out below
# From (3.18): dw1-dot = (I2-I3) w0 dw3 / I1  =>  for growth e^{sigma t}: dw3 = sigma I1 dw1 / ((I2-I3) w0) = -sigma dw1 here.
w_grow[2] = sigma_theory * I123[0] * amp / ((I123[1] - I123[2]) * w0)
assert np.isclose(w_grow[2], -sigma_theory * amp)
sol_g = solve_ivp(euler_rhs, (0, 12), w_grow, args=(I123,), method="DOP853", rtol=1e-11, atol=1e-14,
                  t_eval=np.linspace(0, 12, 4001))
mask = (np.abs(sol_g.y[0]) > 1e-3) & (np.abs(sol_g.y[0]) < 3e-2)
slope = np.polyfit(sol_g.t[mask], np.log(np.abs(sol_g.y[0][mask])), 1)[0]
assert np.isclose(slope, sigma_theory, rtol=0.01), (slope, sigma_theory)
# e-folding: amplitude grows by e^(sigma t); from 1e-4 to 1e-2 takes ln(100)/sigma = 7.98 time units.
t_100 = sol_g.t[np.argmax(np.abs(sol_g.y[0]) > 1e-2)]
assert np.isclose(t_100, np.log(100) / sigma_theory, rtol=0.02)
print(f"Example 4.4: fitted growth rate {slope:.4f} w0, theory 1/sqrt(3) = {sigma_theory:.4f} w0; "
      f"x100 growth in {t_100:.2f}/w0 (theory {np.log(100) / sigma_theory:.2f})")


def oscillation_frequency(t, x):
    """Angular frequency from the mean spacing of zero crossings of x(t) - mean(x)."""
    xc = x - x.mean()
    idx = np.where(np.diff(np.sign(xc)) != 0)[0]
    tc = t[idx] - xc[idx] * (t[idx + 1] - t[idx]) / (xc[idx + 1] - xc[idx])
    return np.pi / np.mean(np.diff(tc))


# Small-oscillation frequencies about the stable axes: axis 1: w0 sqrt((I2-I1)(I3-I1)/(I2 I3)) = w0/sqrt(3); axis 3: w0 sqrt((I3-I1)(I3-I2)/(I1 I2)) = w0.
f1 = oscillation_frequency(runs[0].t, runs[0].y[1])
f3 = oscillation_frequency(runs[2].t, runs[2].y[0])
assert np.isclose(f1, w0 * np.sqrt((I123[1] - I123[0]) * (I123[2] - I123[0]) / (I123[1] * I123[2])), rtol=1e-3)
assert np.isclose(f3, w0 * np.sqrt((I123[2] - I123[0]) * (I123[2] - I123[1]) / (I123[0] * I123[1])), rtol=1e-3)
print(f"Stable-axis wobble frequencies: about axis 1: {f1:.4f} w0 (theory 0.5774), about axis 3: {f3:.4f} w0 (theory 1)")

fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
for axis, ax in zip(range(3), axes):
    sol = runs[axis]
    for k in range(3):
        ax.plot(sol.t, sol.y[k], label=f"$\\omega_{k+1}$")
    ax.set_ylabel(f"spin near axis {axis + 1}\n{'stable' if axis != 1 else 'UNSTABLE (flips)'}")
    ax.grid(alpha=.3)
axes[0].legend(ncol=3)
axes[-1].set_xlabel("t  (units of 1/$\\omega_0$)")
axes[0].set_title("Torque-free Euler equations, $I = (1, 2, 3)$, start within 1% of each principal axis")
fig.tight_layout()
fig.savefig(FIG / "euler_torque_free.png", dpi=120)
plt.close(fig)

# Poinsot: L in body components lives on the sphere |L| = const and on the ellipsoid sum L_i^2/I_i = 2T (polhodes).
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection="3d")
u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color="0.8", linewidth=0.4)
starts = []
for ang in np.linspace(0.08, np.pi / 2 - 0.08, 7):
    starts.append([np.cos(ang), np.sin(ang), 0.0])
    starts.append([np.cos(ang), 0.0, np.sin(ang)])
    starts.append([0.0, np.cos(ang), np.sin(ang)])
for s in starts:
    for sign in (1, -1):
        w_init = sign * np.array(s)
        sol = solve_ivp(euler_rhs, (0, 40), w_init, args=(I123,), method="DOP853", rtol=1e-9, atol=1e-12,
                        t_eval=np.linspace(0, 40, 2000))
        Lb = (I123[:, None] * sol.y)
        Lb /= np.linalg.norm(Lb, axis=0)
        ax.plot(Lb[0], Lb[1], Lb[2], lw=0.8)
for k, c in zip(range(3), "rgb"):
    e = np.zeros(3); e[k] = 1.1
    ax.plot([-e[0], e[0]], [-e[1], e[1]], [-e[2], e[2]], c, lw=1.5)
ax.set_title("Poinsot: $\\mathbf{L}$ in the body frame on the unit sphere ($I = 1, 2, 3$)\nclosed loops around axes 1 and 3; separatrices through axis 2")
ax.set_box_aspect((1, 1, 1))
fig.savefig(FIG / "poinsot_sphere.png", dpi=120)
plt.close(fig)

# Symmetric top, torque-free: thin disk I_3 = 2 I_1 (Example 4.3, the frisbee).
I_disk = np.array([1.0, 1.0, 2.0])
spin_hz, tilt = 5.0, np.radians(5.0)
w_abs = 2 * np.pi * spin_hz
w_init = np.array([w_abs * np.sin(tilt), 0.0, w_abs * np.cos(tilt)])
y0 = np.concatenate([w_init, np.eye(3).ravel()])
t_end = 2.0
sol_sym = solve_ivp(euler_with_attitude_rhs, (0, t_end), y0, args=(I_disk,), method="DOP853", rtol=1e-10, atol=1e-12,
                    t_eval=np.linspace(0, t_end, 20001))
Wb = sol_sym.y[:3]
Rt = sol_sym.y[3:].reshape(3, 3, -1)
assert np.allclose(Wb[2], Wb[2, 0])                                              # w3 constant
Omega_b_theory = w_init[2] * (I_disk[2] - I_disk[0]) / I_disk[0]
Omega_b_meas = oscillation_frequency(sol_sym.t, Wb[0])
assert np.isclose(Omega_b_meas, Omega_b_theory, rtol=2e-3)
assert np.isclose(Omega_b_theory, w_init[2])                                      # thin disk: Omega_b = w3
# Space frame: L constant; e3 precesses about L at Omega_s = |L|/I_1.
L_space = np.einsum("ijt,j,jt->it", Rt, I_disk, Wb)
assert np.allclose(L_space, L_space[:, [0]], atol=1e-8)
L_mag = np.linalg.norm(L_space[:, 0])
Omega_s_theory = L_mag / I_disk[0]
e3_space = Rt[:, 2, :]
Lhat = L_space[:, 0] / L_mag
perp = e3_space - np.outer(Lhat, Lhat @ e3_space)                                   # component of e3 perpendicular to L
u1 = perp[:, 0] / np.linalg.norm(perp[:, 0])
u2 = np.cross(Lhat, u1)
phase = np.unwrap(np.arctan2(u2 @ perp, u1 @ perp))
Omega_s_meas = np.polyfit(sol_sym.t, phase, 1)[0]
assert np.isclose(abs(Omega_s_meas), Omega_s_theory, rtol=1e-3)
ratio = Omega_s_theory / w_init[2]
assert np.isclose(ratio, np.sqrt(np.tan(tilt)**2 + 4), rtol=1e-6)
assert np.isclose(ratio, 2.0019, atol=2e-4)
theta_L = np.degrees(np.arctan(I_disk[0] * w_init[0] / (I_disk[2] * w_init[2])))
assert np.isclose(theta_L, 2.502, atol=2e-3)
assert np.isclose(np.degrees(np.arccos(Lhat @ e3_space[:, 0])), theta_L, atol=1e-6)
print(f"Example 4.3 (frisbee): Omega_b = {Omega_b_meas:.3f} rad/s (= w3 = {w_init[2]:.3f}); "
      f"Omega_s / w3 = {ratio:.4f}; angle(e3, L) = {theta_L:.3f} deg; wobble/spin = 2 to 0.1%")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
ax1.plot(sol_sym.t, Wb[0], label="$\\omega_1$"); ax1.plot(sol_sym.t, Wb[1], label="$\\omega_2$")
ax1.set_xlim(0, 0.5); ax1.set_xlabel("t (s)"); ax1.set_title("Body frame: $\\omega_\\perp$ rotates at $\\Omega_b = \\omega_3$")
ax1.legend(); ax1.grid(alpha=.3)
ax2.plot(sol_sym.t, np.degrees(np.arccos(np.clip(Lhat @ e3_space, -1, 1))), label="angle($\\hat e_3$, L)")
w_space = np.einsum("ijt,jt->it", Rt, Wb)
ax2.plot(sol_sym.t, np.degrees(np.arccos(np.clip(Lhat @ (w_space / np.linalg.norm(w_space, axis=0)), -1, 1))), label="angle($\\omega$, L)")
ax2.set_xlabel("t (s)"); ax2.set_ylabel("degrees"); ax2.set_title("Space frame: both angles constant (cones)")
ax2.legend(); ax2.grid(alpha=.3)
fig.tight_layout()
fig.savefig(FIG / "symmetric_top_free.png", dpi=120)
plt.close(fig)

# Chandler wobble: I_3/I_1 - 1 = 1/305 gives a body-frame period of 305 sidereal days for a rigid Earth.
eps_earth = 1 / 305
Omega_b_earth = (2 * np.pi) * eps_earth          # per sidereal day, with w3 = 2 pi / day
period_days = 2 * np.pi / Omega_b_earth
assert np.isclose(period_days, 305.0)
print(f"Chandler wobble, rigid Earth: period = {period_days:.0f} d (observed about 433 d: the Earth is not rigid)")
print("Part 5 checks passed.")

# %% [markdown]
# ## Part 6: Euler angles z-x'-z'' (§3.11, Example 4.6)
#
# R(phi, theta, psi) = R_z(phi) R_x(theta) R_z(psi). Body components of omega:
# w1 = phi-dot sin(theta) sin(psi) + theta-dot cos(psi), w2 = phi-dot sin(theta) cos(psi) - theta-dot sin(psi),
# w3 = phi-dot cos(theta) + psi-dot. We check it by finite-differencing R along a
# trajectory, and show the Jacobian d(omega)/d(rates) has determinant -sin(theta):
# gimbal lock at theta = 0, pi.


# %%
def euler_R(phi, theta, psi):
    """Rotation matrix of the z-x'-z'' Euler angles, eq. (3.19): r_space = R r_body."""
    return R_z(phi) @ R_x(theta) @ R_z(psi)


def omega_body_from_rates(phi, theta, psi, dphi, dtheta, dpsi):
    """Body components of omega from Euler-angle rates, eq. (3.21)."""
    return np.array([dphi * np.sin(theta) * np.sin(psi) + dtheta * np.cos(psi),
                     dphi * np.sin(theta) * np.cos(psi) - dtheta * np.sin(psi),
                     dphi * np.cos(theta) + dpsi])


def omega_space_from_rates(phi, theta, psi, dphi, dtheta, dpsi):
    """Space components of omega from Euler-angle rates, eq. (3.22)."""
    return np.array([dtheta * np.cos(phi) + dpsi * np.sin(theta) * np.sin(phi),
                     dtheta * np.sin(phi) - dpsi * np.sin(theta) * np.cos(phi),
                     dphi + dpsi * np.cos(theta)])


def rate_matrix(theta, psi):
    """J with omega_body = J (phi-dot, theta-dot, psi-dot); det J = -sin(theta)."""
    return np.array([[np.sin(theta) * np.sin(psi), np.cos(psi), 0.0],
                     [np.sin(theta) * np.cos(psi), -np.sin(psi), 0.0],
                     [np.cos(theta), 0.0, 1.0]])


for _ in range(20):
    ang = rng.uniform(0.2, np.pi - 0.2, 3)
    R = euler_R(*ang)
    assert np.allclose(R.T @ R, np.eye(3)) and np.isclose(np.linalg.det(R), 1)
    # columns are the body axes in space coordinates; the third column is e3 = (sin th sin ph, -sin th cos ph, cos th)
    assert np.allclose(R[:, 2], [np.sin(ang[1]) * np.sin(ang[0]), -np.sin(ang[1]) * np.cos(ang[0]), np.cos(ang[1])])
    # omega by finite difference along a trajectory with random rates
    rates = rng.standard_normal(3)
    def R_traj(t, ang=ang, rates=rates):
        return euler_R(*(ang + rates * t))
    w_space_fd = omega_from_R(R_traj, 0.0)
    assert np.allclose(w_space_fd, omega_space_from_rates(*ang, *rates), atol=1e-6)
    assert np.allclose(R.T @ w_space_fd, omega_body_from_rates(*ang, *rates), atol=1e-6)
    assert np.allclose(rate_matrix(ang[1], ang[2]) @ rates, omega_body_from_rates(*ang, *rates))
    assert np.isclose(np.linalg.det(rate_matrix(ang[1], ang[2])), -np.sin(ang[1]))

# Example 4.6: phi = 30 deg, theta = 60 deg, psi = 45 deg; rates (1, 0.5, 2) rad/s.
ex_ang = np.radians([30.0, 60.0, 45.0])
ex_rates = np.array([1.0, 0.5, 2.0])
wb_ex = omega_body_from_rates(*ex_ang, *ex_rates)
ws_ex = omega_space_from_rates(*ex_ang, *ex_rates)
assert np.allclose(wb_ex, [0.96593, 0.25882, 2.5], atol=1e-5)
assert np.allclose(ws_ex, [1.29904, -1.25, 2.0], atol=1e-5)
assert np.isclose(np.linalg.norm(wb_ex), np.sqrt(7.25)) and np.isclose(np.linalg.norm(ws_ex), np.sqrt(7.25))
assert np.allclose(euler_R(*ex_ang) @ wb_ex, ws_ex)
# |omega|^2 = phi-dot^2 + theta-dot^2 + psi-dot^2 + 2 phi-dot psi-dot cos(theta)
assert np.isclose(np.linalg.norm(wb_ex)**2, ex_rates @ ex_rates + 2 * ex_rates[0] * ex_rates[2] * np.cos(ex_ang[1]))
print(f"Example 4.6: omega_body = {np.round(wb_ex, 5)}, omega_space = {np.round(ws_ex, 5)}, |omega| = {np.linalg.norm(wb_ex):.4f} rad/s")

# Gimbal lock: at theta = 0 the matrix R depends only on phi + psi, and J is singular.
assert np.allclose(euler_R(0.7, 0.0, 0.4), euler_R(0.2, 0.0, 0.9))
thetas = np.logspace(-4, np.log10(np.pi / 2), 200)
conds = np.array([np.linalg.cond(rate_matrix(th, 0.3)) for th in thetas])
assert conds[-1] < 3 and conds[0] > 5e3
assert np.isclose(np.polyfit(np.log(thetas[:50]), np.log(conds[:50]), 1)[0], -1.0, atol=0.02)   # cond ~ 1/theta
fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(thetas, conds)
ax.set_xlabel("$\\theta$ (rad)"); ax.set_ylabel("cond$(J)$")
ax.set_title("Gimbal lock: condition number of the Euler-rate map $\\sim 1/\\theta$")
ax.grid(alpha=.3, which="both")
fig.tight_layout()
fig.savefig(FIG / "gimbal_lock_conditioning.png", dpi=120)
plt.close(fig)
print("Part 6 checks passed: R orthogonal, omega formulas by finite difference, det J = -sin(theta), cond ~ 1/theta.")

# %% [markdown]
# ## Part 7: The heavy symmetric top from its Lagrangian (§3.12, Example 4.5)
#
# L = (I_1/2)(theta-dot^2 + phi-dot^2 sin^2 theta) + (I_3/2)(psi-dot + phi-dot cos theta)^2 - M g l cos theta.
# SymPy forms the Euler–Lagrange equations (as in CM-3's `lagrangian_to_eom`), checks
# that phi and psi are cyclic, verifies the V_eff reduction (3.27) symbolically, and
# hands `solve_ivp` the accelerations. Then: steady precession rates (slow and fast),
# nutation frequency, the sleeping-top threshold, and the axis path on the unit sphere.


# %%
t_sym = sp.symbols("t")
I1s, I3s, Ms, gs, ls = sp.symbols("I_1 I_3 M g ell", positive=True)
th_f, ph_f, ps_f = (sp.Function(n)(t_sym) for n in ("theta", "phi", "psi"))
thd, phd, psd = th_f.diff(t_sym), ph_f.diff(t_sym), ps_f.diff(t_sym)
Lag = I1s / 2 * (thd**2 + phd**2 * sp.sin(th_f)**2) + I3s / 2 * (psd + phd * sp.cos(th_f))**2 - Ms * gs * ls * sp.cos(th_f)

# Cyclic coordinates (CM-4 Theorem 3.1): dL/dphi = dL/dpsi = 0, so p_phi and p_psi are conserved.
assert sp.simplify(sp.diff(Lag, ph_f)) == 0 and sp.simplify(sp.diff(Lag, ps_f)) == 0
p_phi_expr = sp.diff(Lag, phd)
p_psi_expr = sp.diff(Lag, psd)
assert sp.simplify(p_psi_expr - I3s * (psd + phd * sp.cos(th_f))) == 0
assert sp.simplify(p_phi_expr - (I1s * phd * sp.sin(th_f)**2 + I3s * (psd + phd * sp.cos(th_f)) * sp.cos(th_f))) == 0

# Euler–Lagrange equations and the accelerations.
qs = [th_f, ph_f, ps_f]
EL = [sp.diff(sp.diff(Lag, q.diff(t_sym)), t_sym) - sp.diff(Lag, q) for q in qs]
acc_syms = [q.diff(t_sym, 2) for q in qs]
acc_sol = sp.solve(EL, acc_syms, dict=True)[0]

# The theta equation in the form used in the notes (3.26): I_1 theta-ddot = I_1 phi-dot^2 sin cos - p_psi phi-dot sin + M g l sin.
theta_eq_notes = I1s * phd**2 * sp.sin(th_f) * sp.cos(th_f) - p_psi_expr * phd * sp.sin(th_f) + Ms * gs * ls * sp.sin(th_f)
assert sp.simplify(I1s * acc_sol[th_f.diff(t_sym, 2)] - theta_eq_notes) == 0

# Effective potential (3.27): E - p_psi^2/(2 I_3) - I_1 theta-dot^2 / 2 = (p_phi - p_psi cos)^2 / (2 I_1 sin^2) + M g l cos.
pphi, ppsi = sp.symbols("p_phi p_psi")
E_expr = thd * sp.diff(Lag, thd) + phd * p_phi_expr + psd * p_psi_expr - Lag          # energy function h = E (CM-4 (3.8))
phd_of_p = (pphi - ppsi * sp.cos(th_f)) / (I1s * sp.sin(th_f)**2)
psd_of_p = ppsi / I3s - phd_of_p * sp.cos(th_f)
E_reduced = E_expr.subs({psd: psd_of_p}).subs({phd: phd_of_p})
V_eff_sym = (pphi - ppsi * sp.cos(th_f))**2 / (2 * I1s * sp.sin(th_f)**2) + Ms * gs * ls * sp.cos(th_f)
assert sp.simplify(E_reduced - I1s * thd**2 / 2 - ppsi**2 / (2 * I3s) - V_eff_sym) == 0
# Steady precession: dV_eff/dtheta = sin(theta) [p_psi phi-dot - I_1 phi-dot^2 cos(theta) - M g l] with phi-dot = phd_of_p.
dV = sp.diff(V_eff_sym, th_f)
assert sp.simplify(dV - sp.sin(th_f) * (ppsi * phd_of_p - I1s * phd_of_p**2 * sp.cos(th_f) - Ms * gs * ls)) == 0
print("Part 7: SymPy confirms cyclic momenta, the theta equation, the V_eff reduction and the steady-precession condition.")

# Numerical right-hand side for the toy gyroscope of Example 4.5.
M_g, R_g, l_g, g_g = 0.1, 0.03, 0.04, 9.81
I3_g = 0.5 * M_g * R_g**2                      # 4.5e-5 kg m^2
I1_g = 0.25 * M_g * R_g**2 + M_g * l_g**2      # about the pivot: 1.825e-4 kg m^2
params = {I1s: I1_g, I3s: I3_g, Ms: M_g, gs: g_g, ls: l_g}
th_s, dth_s, dph_s, dps_s = sp.symbols("th dth dph dps")
subs_state = {thd: dth_s, phd: dph_s, psd: dps_s, th_f: th_s}
acc_funcs = [sp.lambdify((th_s, dth_s, dph_s, dps_s), acc_sol[a].subs(params).subs(subs_state), "numpy") for a in acc_syms]


def top_rhs(t, y):
    """State y = (theta, phi, psi, theta-dot, phi-dot, psi-dot) for the heavy symmetric top."""
    th, ph, ps, dth, dph, dps = y
    return [dth, dph, dps, acc_funcs[0](th, dth, dph, dps), acc_funcs[1](th, dth, dph, dps), acc_funcs[2](th, dth, dph, dps)]


def momenta(y):
    th, ph, ps, dth, dph, dps = y
    p_psi = I3_g * (dps + dph * np.cos(th))
    p_phi = I1_g * dph * np.sin(th)**2 + p_psi * np.cos(th)
    E = 0.5 * I1_g * (dth**2 + dph**2 * np.sin(th)**2) + p_psi**2 / (2 * I3_g) + M_g * g_g * l_g * np.cos(th)
    return p_phi, p_psi, E


def V_eff(th, p_phi, p_psi):
    return (p_phi - p_psi * np.cos(th))**2 / (2 * I1_g * np.sin(th)**2) + M_g * g_g * l_g * np.cos(th)


def turning_points(p_phi, p_psi, E):
    """Roots in [-1, 1] of f(u) = (alpha - beta u)(1 - u^2) - (b - a u)^2, u = cos(theta), eq. (3.29)."""
    Ep = E - p_psi**2 / (2 * I3_g)
    alpha, beta = 2 * Ep / I1_g, 2 * M_g * g_g * l_g / I1_g
    a_, b_ = p_psi / I1_g, p_phi / I1_g
    # f(u) = beta u^3 - alpha u^2 - (beta + a^2) u + (alpha - b^2 + 2 a b) ... expand explicitly:
    coeffs = [beta, -alpha - a_**2, -beta + 2 * a_ * b_, alpha - b_**2]
    roots = np.roots(coeffs)
    roots = np.sort(roots[np.abs(roots.imag) < 1e-9].real)
    return roots[(roots >= -1 - 1e-9) & (roots <= 1 + 1e-9)]


# Example 4.5: 3000 rpm, axle at theta = 60 degrees from the vertical.
w3_g = 3000 * 2 * np.pi / 60                   # 314.16 rad/s
theta0 = np.radians(60.0)
Mgl = M_g * g_g * l_g
p_psi_g = I3_g * w3_g
# Steady precession: I_1 phi-dot^2 cos(theta) - p_psi phi-dot + M g l = 0.
disc = p_psi_g**2 - 4 * I1_g * np.cos(theta0) * Mgl
phi_slow = (p_psi_g - np.sqrt(disc)) / (2 * I1_g * np.cos(theta0))
phi_fast = (p_psi_g + np.sqrt(disc)) / (2 * I1_g * np.cos(theta0))
phi_slow_approx = Mgl / p_psi_g
phi_fast_approx = p_psi_g / (I1_g * np.cos(theta0))
nut_fast_top = p_psi_g / I1_g
w3_sleep = np.sqrt(4 * Mgl * I1_g) / I3_g
assert np.isclose(I3_g, 4.5e-5) and np.isclose(I1_g, 1.825e-4) and np.isclose(Mgl, 0.03924)
assert np.isclose(phi_slow, 2.827, atol=1e-3) and np.isclose(phi_slow_approx, 2.776, atol=1e-3)
assert np.isclose(phi_fast, 152.1, atol=0.1) and np.isclose(phi_fast_approx, 154.9, atol=0.1)
assert np.isclose(nut_fast_top, 77.46, atol=0.01)
assert np.isclose(w3_sleep, 118.9, atol=0.1) and np.isclose(w3_sleep * 60 / (2 * np.pi), 1136, atol=1)
print(f"Example 4.5: I_3 = {I3_g:.3e}, I_1 = {I1_g:.4e} kg m^2, Mgl = {Mgl:.5f} N m, p_psi = {p_psi_g:.6f}")
print(f"  steady precession at 60 deg: slow {phi_slow:.3f} rad/s (approx Mgl/p_psi = {phi_slow_approx:.3f}), "
      f"fast {phi_fast:.1f} rad/s (approx p_psi/(I_1 cos) = {phi_fast_approx:.1f})")
print(f"  slow precession period {2 * np.pi / phi_slow:.3f} s; fast-top nutation estimate p_psi/I_1 = {nut_fast_top:.2f} rad/s")
print(f"  sleeping-top threshold w3 = {w3_sleep:.1f} rad/s = {w3_sleep * 60 / (2 * np.pi):.0f} rpm")

# (a) Start exactly at slow steady precession: theta and phi-dot must stay constant.
y_steady = [theta0, 0.0, 0.0, 0.0, phi_slow, w3_g - phi_slow * np.cos(theta0)]
sol_st = solve_ivp(top_rhs, (0, 1.0), y_steady, method="DOP853", rtol=1e-10, atol=1e-12, t_eval=np.linspace(0, 1, 2001))
assert np.allclose(sol_st.y[0], theta0, atol=1e-7) and np.allclose(sol_st.y[4], phi_slow, rtol=1e-6)

# (b) Nutation: kick theta-dot. Measured frequency vs sqrt(V_eff''(theta0)/I_1), closed form (3.28).
y_nut = [theta0, 0.0, 0.0, 0.8, phi_slow, w3_g - phi_slow * np.cos(theta0)]
t_end = 2.0
sol_nut = solve_ivp(top_rhs, (0, t_end), y_nut, method="DOP853", rtol=1e-10, atol=1e-12, t_eval=np.linspace(0, t_end, 40001))
p_phi_n, p_psi_n, E_n = momenta(sol_nut.y[:, 0])
P = np.array([momenta(sol_nut.y[:, k]) for k in range(0, 40001, 2000)])
assert np.allclose(P[:, 0], p_phi_n, rtol=1e-8) and np.allclose(P[:, 1], p_psi_n, rtol=1e-8) and np.allclose(P[:, 2], E_n, rtol=1e-8)
w_nut_meas = oscillation_frequency(sol_nut.t, sol_nut.y[0])
w_nut_closed = np.sqrt(((p_psi_g - 2 * I1_g * phi_slow * np.cos(theta0))**2 + I1_g**2 * phi_slow**2 * np.sin(theta0)**2) / I1_g**2)
h = 1e-5
V2 = (V_eff(theta0 + h, p_phi_n, p_psi_n) - 2 * V_eff(theta0, p_phi_n, p_psi_n) + V_eff(theta0 - h, p_phi_n, p_psi_n)) / h**2
assert np.isclose(np.sqrt(V2 / I1_g), w_nut_closed, rtol=1e-4)                   # closed form (3.28) = V_eff''/I_1
assert np.isclose(w_nut_meas, w_nut_closed, rtol=5e-3), (w_nut_meas, w_nut_closed)
assert np.isclose(w_nut_closed, 74.76, atol=0.02)
# Turning points from the cubic agree with the simulated min/max theta.
tp = turning_points(p_phi_n, p_psi_n, E_n)
assert len(tp) == 2
th_lo, th_hi = np.arccos(tp[1]), np.arccos(tp[0])
assert np.isclose(sol_nut.y[0].min(), th_lo, atol=1e-4) and np.isclose(sol_nut.y[0].max(), th_hi, atol=1e-4)
print(f"  nutation: measured {w_nut_meas:.2f} rad/s, closed form sqrt(V_eff''/I_1) = {w_nut_closed:.2f} rad/s "
      f"(fast-top estimate {nut_fast_top:.2f}); theta swings between {np.degrees(th_lo):.2f} and {np.degrees(th_hi):.2f} deg")

# (c) Sleeping top: bisection on w3 for the largest angle reached from theta0 = 0.02 rad at rest.
def theta_max_from_rest(w3, th0=0.02):
    y = [th0, 0.0, 0.0, 0.0, 0.0, w3]
    p_phi, p_psi, E = momenta(y)
    tp = turning_points(p_phi, p_psi, E)
    return np.arccos(tp.min())


lo, hi = 50.0, 300.0
assert theta_max_from_rest(lo) > 0.5 and theta_max_from_rest(hi) < 0.1
for _ in range(40):
    mid = 0.5 * (lo + hi)
    if theta_max_from_rest(mid) > 0.3:
        lo = mid
    else:
        hi = mid
w3_sleep_num = 0.5 * (lo + hi)
assert np.isclose(w3_sleep_num, w3_sleep, rtol=0.01), (w3_sleep_num, w3_sleep)
print(f"  sleeping threshold found numerically: {w3_sleep_num:.1f} rad/s (theory 2 sqrt(M g l I_1)/I_3 = {w3_sleep:.1f})")
# Two integrations straddling the threshold.
sols_sleep = {}
for w3_try in (0.8 * w3_sleep, 1.25 * w3_sleep):
    y = [0.02, 0.0, 0.0, 0.0, 0.0, w3_try]
    sols_sleep[w3_try] = solve_ivp(top_rhs, (0, 3.0), y, method="DOP853", rtol=1e-9, atol=1e-12, t_eval=np.linspace(0, 3, 6001))
assert sols_sleep[0.8 * w3_sleep].y[0].max() > 0.5 and sols_sleep[1.25 * w3_sleep].y[0].max() < 0.1

# (d) A general motion: released from rest at theta = 30 degrees with spin 3000 rpm: cusped nutation.
y_cusp = [np.radians(30.0), 0.0, 0.0, 0.0, 0.0, w3_g]
sol_cusp = solve_ivp(top_rhs, (0, 1.5), y_cusp, method="DOP853", rtol=1e-10, atol=1e-12, t_eval=np.linspace(0, 1.5, 30001))
p_phi_c, p_psi_c, E_c = momenta(y_cusp)
tp_c = turning_points(p_phi_c, p_psi_c, E_c)
assert np.isclose(np.arccos(tp_c.max()), np.radians(30.0), atol=1e-9)             # upper turning point is the release angle
assert np.isclose(sol_cusp.y[0].max(), np.arccos(tp_c.min()), atol=2e-4)
assert np.all(sol_cusp.y[4] >= -1e-9)                                                 # phi-dot never reverses: cusps, not loops
mean_prec = (sol_cusp.y[1, -1] - sol_cusp.y[1, 0]) / 1.5
print(f"  released from rest at 30 deg: theta dips to {np.degrees(sol_cusp.y[0].max()):.3f} deg, "
      f"mean precession {mean_prec:.3f} rad/s (slow-precession formula {Mgl / p_psi_g:.3f})")
assert np.isclose(mean_prec, Mgl / p_psi_g, rtol=0.05)

# Figures: V_eff with turning points; theta(t), phi(t); axis path on the unit sphere.
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
th_grid = np.linspace(0.05, np.pi - 0.05, 2000)
axs[0].plot(np.degrees(th_grid), V_eff(th_grid, p_phi_n, p_psi_n), label="$V_{\\rm eff}(\\theta)$")
axs[0].axhline(E_n - p_psi_n**2 / (2 * I3_g), color="k", ls="--", label="$E'$")
for tpk in (th_lo, th_hi):
    axs[0].axvline(np.degrees(tpk), color="C3", alpha=.6)
axs[0].set_ylim(V_eff(th_grid, p_phi_n, p_psi_n).min() - 0.02, E_n - p_psi_n**2 / (2 * I3_g) + 0.15)
axs[0].set_xlim(40, 80); axs[0].set_xlabel("$\\theta$ (deg)"); axs[0].legend(); axs[0].grid(alpha=.3)
axs[0].set_title("Effective potential and turning points")
axs[1].plot(sol_nut.t, np.degrees(sol_nut.y[0]))
axs[1].set_xlim(0, 0.5); axs[1].set_xlabel("t (s)"); axs[1].set_ylabel("$\\theta$ (deg)"); axs[1].grid(alpha=.3)
axs[1].set_title(f"Nutation at {w_nut_meas:.1f} rad/s on top of precession at {phi_slow:.2f} rad/s")
for w3_try, lab in zip(sols_sleep, ("below threshold", "above threshold")):
    axs[2].plot(sols_sleep[w3_try].t, np.degrees(sols_sleep[w3_try].y[0]), label=f"{lab}: $\\omega_3$ = {w3_try:.0f} rad/s")
axs[2].set_xlabel("t (s)"); axs[2].set_ylabel("$\\theta$ (deg)"); axs[2].legend(); axs[2].grid(alpha=.3)
axs[2].set_title("Sleeping top: released at 0.02 rad")
fig.tight_layout()
fig.savefig(FIG / "heavy_top_veff_nutation.png", dpi=120)
plt.close(fig)

fig = plt.figure(figsize=(12, 6))
for k, (sol, title) in enumerate([(sol_nut, "kicked from steady precession"), (sol_cusp, "released from rest: cusps")]):
    ax = fig.add_subplot(1, 2, k + 1, projection="3d")
    th_t, ph_t = sol.y[0], sol.y[1]
    e3 = np.array([np.sin(th_t) * np.sin(ph_t), -np.sin(th_t) * np.cos(ph_t), np.cos(th_t)])
    ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color="0.85", linewidth=0.3)
    ax.plot(e3[0], e3[1], e3[2], lw=0.7)
    ax.set_title(f"Path of $\\hat e_3$ on the unit sphere: {title}")
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=35, azim=-60)
fig.savefig(FIG / "heavy_top_axis_path.png", dpi=120)
plt.close(fig)
print("Part 7 checks passed.")

# %% [markdown]
# ## Part 8: Precession of the equinoxes (§3.13, Example 4.7)
#
# Quadrupole torque on an oblate body: tau = (3 G M / r^3) r-hat x (I r-hat). Averaged over
# a circular orbit in the ecliptic and applied to L = I_3 w3 e3 it gives a retrograde
# precession of e3 about the ecliptic pole at |Omega| = (3/2)(n^2 / w3) H cos(eps), with
# n^2 = GM/r^3 and H = (I_3 - I_1)/I_3 = 1/305.45. Sun plus Moon give about 51 arcsec/yr,
# a period of about 25,500 yr; the observed value is 50.3 arcsec/yr (25,772 yr).


# %%
GM_sun, r_sun = 1.32712e20, 1.495979e11
GM_moon, r_moon = 4.9028e12, 3.844e8
w3_earth = 7.2921e-5
H_earth = 1 / 305.45
eps_obliq = np.radians(23.44)
year = 365.25 * 86400
arcsec = np.degrees(1) * 3600

# Check the torque formula on a random symmetric body against the direct sum over point masses in a 1/r^2 field.
I_sym = np.diag([1.0, 1.0, 1.3])
masses_q = np.full(4000, 1.0 / 4000)
pts_q = sample_sphere(1.0, 4000) * np.array([1.0, 1.0, 1.5])           # oblate-ish cloud (not exactly I_sym; we use its own I)
pts_q -= pts_q.mean(axis=0)
I_q = inertia_tensor(masses_q, pts_q)
r_vec = 200.0 * np.array([0.3, -0.5, 0.8]) / np.linalg.norm([0.3, -0.5, 0.8])   # attracting mass M sits at r (GM = 1)
tau_direct = np.zeros(3)
for mq, pq in zip(masses_q, pts_q):
    d = r_vec - pq                                                      # from the mass element at r_a to the attractor at r
    tau_direct += mq * np.cross(pq, d / np.linalg.norm(d)**3)           # F_a = GM m_a (r - r_a)/|r - r_a|^3, tau = sum r_a x F_a
rhat = r_vec / np.linalg.norm(r_vec)
tau_quad = 3 / np.linalg.norm(r_vec)**3 * np.cross(rhat, I_q @ rhat)
assert np.allclose(tau_direct, tau_quad, rtol=0.02, atol=1e-9 * np.linalg.norm(tau_quad))

Omega_sun = 1.5 * (GM_sun / r_sun**3) / w3_earth * H_earth * np.cos(eps_obliq)
Omega_moon = 1.5 * (GM_moon / r_moon**3) / w3_earth * H_earth * np.cos(eps_obliq)
rate_sun = Omega_sun * year * arcsec
rate_moon = Omega_moon * year * arcsec
rate_total = rate_sun + rate_moon
period_yr = 360 * 3600 / rate_total
assert np.isclose(rate_sun, 15.94, atol=0.05)
assert np.isclose(Omega_moon / Omega_sun, 2.18, atol=0.01)
assert 45 < rate_total < 55 and 23_000 < period_yr < 29_000
print(f"Example 4.7: solar {rate_sun:.2f} arcsec/yr, lunar {rate_moon:.2f} arcsec/yr (ratio {Omega_moon / Omega_sun:.2f}), "
      f"total {rate_total:.1f} arcsec/yr, period {period_yr:.0f} yr (observed 50.3 arcsec/yr, 25,772 yr)")
print("Part 8 checks passed.")

# %% [markdown]
# ## Your turn
#
# 1. **Satellite with a deployed boom.** A spacecraft spins about its axis of largest
#    moment (I_3). At t = 10 s it extends a boom along e_1, which raises I_2 and I_3 but
#    leaves I_1 alone. Modify `euler_rhs` to take I(t) and make I_3 the *intermediate*
#    moment after deployment. Does the spin flip? Angular momentum in space is still
#    conserved; what happens to T = (1/2) w.I.w while I changes? (Hint: Euler's
#    equations with time-dependent I pick up a term I-dot w; derive it from
#    d(I w)/dt in the body frame before coding.)
# 2. **Dzhanibekov T-handle.** Build a T-handle from two rods (use `inertia_tensor`
#    on a dense sample of two cylinders), find its principal axes with `eigh`, and
#    integrate Euler's equations from a spin about the intermediate axis with a 1%
#    perturbation. Carry the attitude R along as in `euler_with_attitude_rhs` and
#    plot the handle's axis in space: it flips by 180 degrees periodically while L
#    stays fixed. Measure the flip period and compare with 1/sigma from §3.10.
# 3. **Spinning coin.** A coin of radius a rolling in a circle while leaning at angle
#    alpha (Euler's disk, slow regime). Use the rolling constraint and the steady
#    precession condition for a disk in contact with the table to derive the relation
#    between the precession rate Omega and the lean angle, then check it by extending
#    the heavy-top code with the contact-point constraint forces as Lagrange
#    multipliers (CM-3 §3.8). Why can you not just substitute the rolling constraint
#    into L (CM-3 §3.2)?

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s. Figures in {FIG}")
