# %% [markdown]
# # LA-7 Project: Eigenvalues in Action
#
# Every object in LA-7 computed and checked: (0) a 2x2 matrix acting on the
# unit circle, with the two directions it only stretches; (1) power iteration
# and shifted inverse iteration against `np.linalg.eig`, with the convergence
# rate |lambda_2/lambda_1|^k of eq. (3.9) measured; (2) the weather Markov
# chain of Example 4.1, its steady state, and a simulation; (3) Fibonacci via
# A^k = X Lambda^k X^{-1} and Binet's formula; (4) du/dt = Au by
# eigendecomposition versus `solve_ivp`; (5) the matrix exponential by series
# versus `scipy.linalg.expm`, d/dt e^{At} = A e^{At}, and the e^A e^B != e^{A+B}
# counterexample; (6) phase portraits of a saddle, a stable node, a stable
# spiral and a center, and the trace-determinant plane; (7) the two-mass three-spring normal modes of
# Example 4.3, animated as frames; (8) trace/det identities, similar matrices,
# the defective shear, Example 4.5, and the symmetric 3x3 of Example 4.6.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
from fractions import Fraction
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.linalg import expm, eigh, lu_factor, lu_solve

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
SQ3 = np.sqrt(3.0)


def same_direction(v, w, tol=1e-8):
    """True if v and w are nonzero scalar multiples of each other (real or complex
    scalar). Uses the equality case of the complex Cauchy-Schwarz inequality
    (LA-8, Theorem 3.4; LA-1, Theorem 3.11 for real vectors): |<v, w>| = ||v|| ||w||.
    This is how to compare eigenvectors, which are only defined up to scale
    (Theorem 3.2 and section 8 of the notes)."""
    v = np.asarray(v, complex).ravel()
    w = np.asarray(w, complex).ravel()
    nv, nw = np.linalg.norm(v), np.linalg.norm(w)
    return abs(abs(np.vdot(v, w)) - nv * nw) <= tol * max(1.0, nv * nw)


def eig_sorted(A):
    """Eigenvalues and eigenvectors of A, sorted by decreasing |lambda|."""
    lam, X = np.linalg.eig(np.asarray(A, float))
    order = np.argsort(-np.abs(lam), kind="stable")
    return lam[order], X[:, order]


# %% [markdown]
# ## Part 0: A matrix acting on the unit circle (section 2)
#
# A = [[2, 1], [1, 2]] sends unit vectors to points on an ellipse. Most are
# rotated; (1,1) and (1,-1) are only stretched (by 3 and by 1). The table in
# section 2 is reproduced here.

# %%
A_sec2 = np.array([[2.0, 1.0], [1.0, 2.0]])
th = np.linspace(0, 2 * np.pi, 48, endpoint=False)
U = np.c_[np.cos(th), np.sin(th)]                 # rows: unit vectors
V = U @ A_sec2.T                                  # rows: their images A u
rot = np.degrees((np.arctan2(V[:, 1], V[:, 0]) - th + np.pi) % (2 * np.pi) - np.pi)
stretch = np.linalg.norm(V, axis=1)


def rotation_and_stretch(u):
    """Angle (degrees) by which A_sec2 rotates u, and the factor it stretches it."""
    u = np.asarray(u, float); v = A_sec2 @ u
    ang = np.degrees(np.arctan2(v[1], v[0]) - np.arctan2(u[1], u[0]))
    ang = (ang + 180) % 360 - 180
    return ang, np.linalg.norm(v) / np.linalg.norm(u)


# The table in section 2.
for u, exp_rot, exp_str in [([1, 0], 26.565, np.sqrt(5)), ([0, 1], -26.565, np.sqrt(5)),
                            ([2, 1], 12.095, np.sqrt(41 / 5)), ([1, 2], -12.095, np.sqrt(41 / 5)),
                            ([1, 1], 0.0, 3.0), ([1, -1], 0.0, 1.0)]:
    r, s = rotation_and_stretch(u)
    assert abs(r - exp_rot) < 1e-2 and abs(s - exp_str) < 1e-9, (u, r, s)
print("Section 2 table reproduced: (1,1) stretched x3 and (1,-1) x1 with no rotation.")

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4.8))
for u, v, r in zip(U, V, rot):
    ax0.annotate("", xy=u, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="0.6", lw=0.8))
    ax0.annotate("", xy=v, xytext=(0, 0),
                 arrowprops=dict(arrowstyle="->", color=plt.cm.coolwarm((r + 30) / 60), lw=1.0))
t = np.linspace(-3.3, 3.3, 2)
ax0.plot(t, t, "k--", lw=1, label="eigenline (1,1): stretched x3")
ax0.plot(t, -t, "k:", lw=1, label="eigenline (1,-1): stretched x1")
ax0.set_aspect("equal"); ax0.set_xlim(-3.3, 3.3); ax0.set_ylim(-3.3, 3.3); ax0.grid(alpha=.3)
ax0.legend(fontsize=8, loc="lower right"); ax0.set_title("unit vectors (gray) and A u (colored by rotation)")
ax1.plot(np.degrees(th), rot, "o-", ms=3)
ax1.axhline(0, color="k", lw=.8)
for ang in (45, 135, 225, 315):
    ax1.axvline(ang, color="k", ls="--", lw=.8)
ax1.set_xlabel("direction of u (degrees)"); ax1.set_ylabel("rotation of A u (degrees)")
ax1.set_title("zero rotation exactly at the eigen-directions"); ax1.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "unit_vectors_and_images.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 1: Power iteration and shifted inverse iteration (section 3.6)
#
# **Build it yourself first:** multiply, normalize, repeat. The Rayleigh
# quotient v^T A v / v^T v estimates the eigenvalue. Inverse iteration is the
# same loop with a solve against (A - sigma I) in place of the multiply.
# Equation (3.9): the direction error decays like |lambda_2/lambda_1|^k.

# %%
def power_iteration(A, k_max=500, tol=1e-13, v0=None):
    """Power iteration (section 3.6): v_{k+1} = A v_k / ||A v_k||.
    Returns (rayleigh_estimate, unit_vector, n_steps). Converges to the eigenvector
    of the eigenvalue of largest modulus at rate |lambda_2/lambda_1|^k, eq. (3.9)."""
    A = np.asarray(A, float)
    v = rng.standard_normal(A.shape[0]) if v0 is None else np.asarray(v0, float)
    v = v / np.linalg.norm(v)
    for k in range(1, k_max + 1):
        w = A @ v
        v_new = w / np.linalg.norm(w)
        # Direction converged? Compare after aligning signs (lambda_1 < 0 flips v each step).
        if np.linalg.norm(v_new - np.sign(v_new @ v) * v) < tol:
            v = v_new
            break
        v = v_new
    return float(v @ A @ v), v, k


def inverse_iteration(A, shift, k_max=200, tol=1e-13, v0=None):
    """Shifted inverse iteration (section 3.6): power iteration on (A - shift I)^{-1},
    whose eigenvalues are 1/(lambda_i - shift) with the same eigenvectors (Theorem 3.5b).
    Converges to the eigenpair whose lambda is closest to `shift`. One LU factorization,
    then a triangular solve per step."""
    A = np.asarray(A, float); n = A.shape[0]
    lu = lu_factor(A - shift * np.eye(n))
    v = rng.standard_normal(n) if v0 is None else np.asarray(v0, float)
    v = v / np.linalg.norm(v)
    for k in range(1, k_max + 1):
        w = lu_solve(lu, v)
        v_new = w / np.linalg.norm(w)
        if np.linalg.norm(v_new - np.sign(v_new @ v) * v) < tol:
            v = v_new
            break
        v = v_new
    return float(v @ A @ v), v, k


def convergence_rate(A, v0, k_max=60, floor=1e-10, shift=None):
    """Exercise B8: direction error ||v_k - x_1|| and Rayleigh error |lambda^(k) - lambda_1|
    along the iteration, with fitted slopes of log(error) vs k, and the predicted slope.
    shift=None: power iteration; x_1 belongs to the largest |lambda|, predicted slope
    log|lambda_2/lambda_1| (eq. 3.9).  shift=sigma: inverse iteration on (A - sigma I)^{-1};
    x_1 belongs to the lambda closest to sigma, predicted slope
    log(|lambda_closest - sigma| / |lambda_next - sigma|) (section 3.6)."""
    A = np.asarray(A, float); n = A.shape[0]
    lam, X = np.linalg.eig(A)
    if shift is None:
        order = np.argsort(-np.abs(lam), kind="stable")
        step = lambda v: A @ v
    else:
        order = np.argsort(np.abs(lam - shift), kind="stable")
        lu = lu_factor(A - shift * np.eye(n))
        step = lambda v: lu_solve(lu, v)
    lam, X = lam[order], X[:, order]
    x1 = X[:, 0].real; x1 = x1 / np.linalg.norm(x1); lam1 = lam[0].real
    v = np.asarray(v0, float); v = v / np.linalg.norm(v)
    dir_err, ray_err = [], []
    for _ in range(k_max):
        s = np.sign(v @ x1) or 1.0
        dir_err.append(np.linalg.norm(s * v - x1))
        ray_err.append(abs(v @ A @ v - lam1))
        w = step(v); v = w / np.linalg.norm(w)
    dir_err, ray_err = np.array(dir_err), np.array(ray_err)

    def slope(err):
        k = np.arange(len(err)); ok = (err > floor) & (err < 0.5)
        return np.polyfit(k[ok], np.log(err[ok]), 1)[0]

    if shift is None:
        pred = np.log(abs(lam[1] / lam[0]))
    else:
        pred = np.log(abs(lam[0] - shift) / abs(lam[1] - shift))
    return dir_err, ray_err, slope(dir_err), slope(ray_err), pred


# Exercise A1's matrix: eigenvalues 5 and 2, eigenvectors (1,1) and (1,-2).
A1 = np.array([[4.0, 1.0], [2.0, 3.0]])
lam_pi, v_pi, steps = power_iteration(A1, v0=[1, 0])
assert abs(lam_pi - 5) < 1e-10 and same_direction(v_pi, [1, 1]), (lam_pi, v_pi)
lam_ii, v_ii, steps_ii = inverse_iteration(A1, shift=2.2, v0=[1, 0])
assert abs(lam_ii - 2) < 1e-10 and same_direction(v_ii, [1, -2]), (lam_ii, v_ii)
print(f"Power iteration on A1: lambda = {lam_pi:.12f} in {steps} steps; "
      f"inverse iteration (shift 2.2): lambda = {lam_ii:.12f} in {steps_ii} steps.")

# Compare with np.linalg.eig on a random symmetric 6x6 with a clear gap.
S6 = rng.standard_normal((6, 6)); S6 = S6 + S6.T
lam_ref, X_ref = eig_sorted(S6)
lam_top, v_top, _ = power_iteration(S6)
assert abs(lam_top - lam_ref[0].real) < 1e-9 and same_direction(v_top, X_ref[:, 0], tol=1e-6)
sigma = lam_ref[3].real + 0.05 * (lam_ref[2].real - lam_ref[3].real)
lam_mid, v_mid, _ = inverse_iteration(S6, shift=sigma)
assert abs(lam_mid - lam_ref[3].real) < 1e-9 and same_direction(v_mid, X_ref[:, 3], tol=1e-6)
print("Random symmetric 6x6: power iteration matches eig's largest eigenvalue; "
      "shifted inverse iteration lands on the eigenvalue nearest the shift.")

# Convergence rate, eq. (3.9): slope of log(error) should be log|lambda_2/lambda_1|.
d1, r1, s_dir1, s_ray1, pred1 = convergence_rate(A1, [1, 0])            # ratio 2/5
d2, r2, s_dir2, s_ray2, pred2 = convergence_rate(A_sec2, [1, 0])        # ratio 1/3, symmetric
assert abs(s_dir1 - pred1) < 0.03 * abs(pred1), (s_dir1, pred1)
assert abs(s_dir2 - pred2) < 0.03 * abs(pred2), (s_dir2, pred2)
assert abs(s_ray2 - 2 * pred2) < 0.05 * abs(2 * pred2), (s_ray2, 2 * pred2)  # B8(b): doubled slope
assert abs(s_ray1 - pred1) < 0.05 * abs(pred1), (s_ray1, pred1)              # B8(b): nonsymmetric, first order
# B8(c): shifted inverse iteration on A1 with sigma = 2.2 targets lambda = 2;
# predicted ratio |2 - 2.2| / |5 - 2.2| = 1/14.
d3, r3, s_dir3, s_ray3, pred3 = convergence_rate(A1, [1, 0], k_max=30, shift=2.2)
assert np.isclose(pred3, np.log(1 / 14))
assert abs(s_dir3 - pred3) < 0.05 * abs(pred3), (s_dir3, pred3)
print(f"Direction-error slopes: {s_dir1:.3f} (predicted log(2/5) = {pred1:.3f}), "
      f"{s_dir2:.3f} (predicted log(1/3) = {pred2:.3f}).")
print(f"Rayleigh-error slope for the symmetric matrix: {s_ray2:.3f} (predicted 2 log(1/3) = {2*pred2:.3f}); "
      f"for the nonsymmetric one: {s_ray1:.3f} (first order, like the direction; predicted {pred1:.3f}).")
print(f"Inverse iteration on A1, shift 2.2: direction slope {s_dir3:.3f} (predicted log(1/14) = {pred3:.3f}), "
      f"{-s_dir3 / np.log(10):.2f} digits per step versus {-s_dir1 / np.log(10):.2f} for power iteration.")

fig, ax = plt.subplots(figsize=(6.5, 4.5))
k = np.arange(len(d1))
ax.semilogy(k, d1, "o-", ms=3, label="direction error, A = [[4,1],[2,3]]")
ax.semilogy(k, d1[0] * (2 / 5) ** k, "k--", lw=1, label="(2/5)^k")
ax.semilogy(k, d2, "s-", ms=3, label="direction error, A = [[2,1],[1,2]]")
ax.semilogy(k, d2[0] * (1 / 3) ** k, "k:", lw=1, label="(1/3)^k")
ax.semilogy(k, r2, "^-", ms=3, label="Rayleigh error, [[2,1],[1,2]]")
ax.semilogy(k, r2[0] * (1 / 9) ** k, "k-.", lw=1, label="(1/9)^k")
ax.set_ylim(1e-16, 10); ax.set_xlim(0, 40)
ax.set_xlabel("iteration k"); ax.set_ylabel("error"); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.set_title("Power iteration converges like |lambda_2/lambda_1|^k, eq. (3.9)")
fig.tight_layout(); fig.savefig(FIG / "power_iteration_convergence.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 2: The weather Markov chain (Example 4.1, section 3.7)
#
# A = [[0.9, 0.5], [0.1, 0.5]]: eigenvalues 1 and 0.4, steady state (5/6, 1/6).
# Starting rainy, u_k = (5/6, 1/6) - (5/6)(0.4)^k (1, -1). A 20,000-day
# simulation should spend about 5/6 of its days sunny.

# %%
A_w = np.array([[0.9, 0.5], [0.1, 0.5]])
assert np.all(A_w >= 0) and np.allclose(A_w.sum(axis=0), 1), "not a Markov matrix"
lam_w, X_w = eig_sorted(A_w)
assert np.allclose(lam_w, [1.0, 0.4]), lam_w
steady = X_w[:, 0].real / X_w[:, 0].real.sum()
assert np.allclose(steady, [5 / 6, 1 / 6]), steady
assert same_direction(X_w[:, 1], [1, -1])
assert abs(X_w[:, 1].real.sum()) < 1e-12, "transient eigenvector must sum to zero (B4b)"

u = np.array([0.0, 1.0])                                   # rainy today
traj = [u.copy()]
for k in range(30):
    u = A_w @ u; traj.append(u.copy())
traj = np.array(traj)
assert np.allclose(traj[1], [0.5, 0.5]) and np.allclose(traj[2], [0.7, 0.3])
ks = np.arange(31)[:, None]
closed = np.array([5 / 6, 1 / 6]) - (5 / 6) * 0.4 ** ks * np.array([1.0, -1.0])
assert np.allclose(traj, closed), "closed form (3.8) disagrees with iteration"
assert np.allclose(traj[-1], steady, atol=1e-11)
print("Weather chain: u_1 = (0.5, 0.5), u_2 = (0.7, 0.3), u_30 ~ (5/6, 1/6); closed form matches iteration.")


def simulate_markov(A, state0, n_steps):
    """Sample a trajectory of the chain with column-stochastic A: P(next = i | now = j) = A[i, j]."""
    n = A.shape[0]; s = state0; states = np.empty(n_steps, int)
    for k in range(n_steps):
        s = rng.choice(n, p=A[:, s]); states[k] = s
    return states


N_DAYS = 20_000
path = simulate_markov(A_w, 1, N_DAYS)
empirical = np.bincount(path, minlength=2) / N_DAYS
assert np.abs(empirical - steady).max() < 0.015, (empirical, steady)
print(f"Simulated {N_DAYS} days: empirical (sunny, rainy) = {empirical}, eigenvector says {steady}.")

# Theorem 3.11 on random Markov matrices: 1 is an eigenvalue, all |lambda| <= 1.
for _ in range(200):
    n = rng.integers(2, 7)
    Mk = rng.random((n, n)); Mk /= Mk.sum(axis=0)
    ev = np.linalg.eigvals(Mk)
    assert np.min(np.abs(ev - 1)) < 1e-8 and np.max(np.abs(ev)) <= 1 + 1e-10
print("Theorem 3.11 verified on 200 random Markov matrices.")

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4))
ax0.plot(traj[:12, 0], "o-", label="P(sunny)"); ax0.plot(traj[:12, 1], "s-", label="P(rainy)")
ax0.axhline(5 / 6, color="C0", ls="--", lw=.8); ax0.axhline(1 / 6, color="C1", ls="--", lw=.8)
ax0.set_xlabel("day k"); ax0.set_title("u_k from a rainy day; dashed = steady state"); ax0.legend(); ax0.grid(alpha=.3)
running = np.cumsum(path == 0) / np.arange(1, N_DAYS + 1)
ax1.plot(running, lw=1); ax1.axhline(5 / 6, color="k", ls="--", lw=.8)
ax1.set_xscale("log"); ax1.set_ylim(0.5, 1.0); ax1.set_xlabel("days simulated")
ax1.set_title("running fraction of sunny days -> 5/6"); ax1.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "markov_weather.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 3: Fibonacci (Example 4.2)
#
# u_{k+1} = [[1,1],[1,0]] u_k.  F_k = (phi^k - psi^k)/sqrt(5), Binet's formula.
# Exact integers from the recurrence; floats from the closed form; A^k from (3.7).

# %%
A_f = np.array([[1.0, 1.0], [1.0, 0.0]])
lam_f, X_f = eig_sorted(A_f)
phi, psi = (1 + np.sqrt(5)) / 2, (1 - np.sqrt(5)) / 2
assert np.allclose(lam_f, [phi, psi])
assert same_direction(X_f[:, 0], [phi, 1]) and same_direction(X_f[:, 1], [psi, 1])

F = [0, 1]
for k in range(2, 101):
    F.append(F[-1] + F[-2])
binet = lambda k: (phi ** k - psi ** k) / np.sqrt(5)
first_bad = next(k for k in range(101) if round(binet(k)) != F[k])
assert first_bad > 70, first_bad
# Binet in doubles is good to ~3e-15 relative (14-15 digits), not 16: rounding phi to a double
# (relative error delta ~ 3e-17) is amplified to ~100 delta in phi^100 (Example 4.2).
rel_binet = abs(binet(100) - F[100]) / F[100]
assert rel_binet < 1e-14, rel_binet
assert F[100] == 354224848179261915075
# A^k via (3.7) in floating point, compared with the exact integers.
for k in (5, 10, 30, 60):
    Ak = (X_f * lam_f ** k) @ np.linalg.inv(X_f)
    assert np.allclose(Ak, [[F[k + 1], F[k]], [F[k], F[k - 1]]], rtol=1e-12)
# Exact A^100 with Python integers.
A100 = np.linalg.matrix_power(np.array([[1, 1], [1, 0]], dtype=object), 100)
assert A100[0, 1] == F[100]
print(f"F_100 = {F[100]} (exact); Binet in floating point = {binet(100):.15e} "
      f"(relative error {rel_binet:.1e}); rounding Binet first fails at k = {first_bad}.")

# %% [markdown]
# ## Part 4: du/dt = Au by eigendecomposition and by solve_ivp (sections 3.8, 3.11)
#
# **Build it yourself first:** eigendecompose, c = X^{-1} u0, assemble
# sum_i c_i e^{lambda_i t} x_i (eq. 3.11), take the real part (eq. 3.18).

# %%
def solve_linear_ode(A, u0, t):
    """u(t) = sum_i c_i e^{lambda_i t} x_i, eq. (3.11), with c = X^{-1} u0.
    For real A and u0 the imaginary parts cancel (eq. 3.18); we assert that and return the real part.
    Returns an array of shape (n, len(t))."""
    lam, X = np.linalg.eig(np.asarray(A, float))
    c = np.linalg.solve(X, np.asarray(u0, float))
    Ut = (X * c) @ np.exp(np.outer(lam, np.asarray(t, float)))
    assert np.abs(Ut.imag).max() < 1e-9 * max(1.0, np.abs(Ut.real).max())
    return Ut.real


def solve_ivp_linear(A, u0, t):
    sol = solve_ivp(lambda _t, u: A @ u, (t[0], t[-1]), np.asarray(u0, float),
                    t_eval=t, rtol=1e-10, atol=1e-12, method="DOP853")
    return sol.y


# Example 4.4: the spiral. u(t) = e^{-t} (2 cos 2t - sin 2t, sin 2t).
A_sp = np.array([[-2.0, -5.0], [1.0, 0.0]])
u0_sp = np.array([2.0, 0.0])
lam_sp = np.linalg.eigvals(A_sp)
assert np.allclose(sorted(lam_sp, key=lambda z: z.imag), [-1 - 2j, -1 + 2j])
t_sp = np.linspace(0, 6, 601)
U_eig = solve_linear_ode(A_sp, u0_sp, t_sp)
U_ivp = solve_ivp_linear(A_sp, u0_sp, t_sp)
U_ana = np.exp(-t_sp) * np.array([2 * np.cos(2 * t_sp) - np.sin(2 * t_sp), np.sin(2 * t_sp)])
assert np.allclose(U_eig, U_ana, atol=1e-12) and np.allclose(U_ivp, U_ana, atol=1e-7)
u07 = expm(A_sp * 0.7) @ u0_sp
assert np.allclose(u07, [-0.320553, 0.489360], atol=1e-6), u07
print(f"Spiral: u(0.7) = {u07} by expm; eigen-solution, analytic formula and solve_ivp agree.")

# Exercise B2: the saddle. u(t) = 2 e^{3t} (1,1) + e^{-t} (1,-1).
A_sd = np.array([[1.0, 2.0], [2.0, 1.0]])
u0_sd = np.array([3.0, 1.0])
t_sd = np.linspace(0, 1.5, 151)
U_sd = solve_linear_ode(A_sd, u0_sd, t_sd)
U_sd_ana = 2 * np.exp(3 * t_sd) * np.array([[1.0], [1.0]]) + np.exp(-t_sd) * np.array([[1.0], [-1.0]])
assert np.allclose(U_sd, U_sd_ana) and np.allclose(solve_ivp_linear(A_sd, u0_sd, t_sd), U_sd_ana, rtol=1e-7)
print("Saddle (B2): eigen-solution matches 2e^{3t}(1,1) + e^{-t}(1,-1) and solve_ivp.")

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
axes[0].plot(t_sp, U_eig[0], label="u1, eigen (3.11)"); axes[0].plot(t_sp, U_eig[1], label="u2, eigen (3.11)")
axes[0].plot(t_sp[::20], U_ivp[0, ::20], "ko", ms=3, label="solve_ivp"); axes[0].plot(t_sp[::20], U_ivp[1, ::20], "ko", ms=3)
axes[0].set_title("Example 4.4: components"); axes[0].set_xlabel("t"); axes[0].legend(fontsize=8); axes[0].grid(alpha=.3)
axes[1].plot(U_eig[0], U_eig[1]); axes[1].plot(*u0_sp, "ko"); axes[1].set_aspect("equal")
axes[1].set_title("Example 4.4: stable spiral, lambda = -1 +/- 2i"); axes[1].grid(alpha=.3)
axes[2].plot(t_sd, U_sd[0], label="u1"); axes[2].plot(t_sd, U_sd[1], label="u2")
axes[2].plot(t_sd[::10], solve_ivp_linear(A_sd, u0_sd, t_sd)[:, ::10].T, "k.", ms=4)
axes[2].set_yscale("log"); axes[2].set_title("B2: saddle, growth along (1,1)"); axes[2].legend(fontsize=8); axes[2].grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "ode_eigen_vs_solve_ivp.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5: The matrix exponential (section 3.9)
#
# e^A = sum A^k/k!, eq. (3.12). Checks: truncated series vs `expm`;
# e^{At} = X e^{Lambda t} X^{-1}, eq. (3.13); d/dt e^{At} = A e^{At} by finite
# differences; e^A e^B = e^{A+B} only when AB = BA; e^{J theta} = R(theta), eq. (3.20).

# %%
def expm_series(A, n_terms=40):
    """e^A by the truncated power series (3.12): sum_{k<n_terms} A^k / k!."""
    A = np.asarray(A, float); n = A.shape[0]
    S = np.eye(n); term = np.eye(n)
    for k in range(1, n_terms):
        term = term @ A / k
        S = S + term
    return S


for _ in range(20):
    n = rng.integers(2, 6)
    A = rng.standard_normal((n, n))
    assert np.allclose(expm_series(A), expm(A), atol=1e-11), "series vs scipy expm"
    # (3.13): e^{At} = X e^{Lambda t} X^{-1}, complex arithmetic allowed, result real.
    lam, X = np.linalg.eig(A); tt = 0.8
    E = (X * np.exp(lam * tt)) @ np.linalg.inv(X)
    assert np.allclose(E.real, expm(A * tt), atol=1e-9) and np.abs(E.imag).max() < 1e-9
    # Theorem 3.15(b): d/dt e^{At} = A e^{At} = e^{At} A (central difference).
    h = 1e-5
    dE = (expm(A * (tt + h)) - expm(A * (tt - h))) / (2 * h)
    assert np.allclose(dE, A @ expm(A * tt), rtol=1e-6, atol=1e-8)
    assert np.allclose(A @ expm(A * tt), expm(A * tt) @ A)
    # Theorem 3.15(d): e^{At} e^{-At} = I.
    assert np.allclose(expm(A * tt) @ expm(-A * tt), np.eye(n))
print("Series == expm; (3.13) holds; d/dt e^{At} = A e^{At}; e^{At} e^{-At} = I on 20 random matrices.")

# Example 4.5: e^{At} = e^{-t}(I + tN) for the critically damped oscillator.
A_cd = np.array([[0.0, 1.0], [-1.0, -2.0]])
N_cd = A_cd + np.eye(2)
assert np.allclose(N_cd @ N_cd, 0)
assert np.allclose(np.linalg.eigvals(A_cd), [-1, -1]) and np.linalg.matrix_rank(A_cd + np.eye(2)) == 1
t13 = 1.3
assert np.allclose(expm(A_cd * t13), np.exp(-t13) * (np.eye(2) + t13 * N_cd))
x_cd = lambda t: (1 + t) * np.exp(-t)
assert np.allclose((expm(A_cd * t13) @ [1, 0])[0], x_cd(t13))
print("Example 4.5: e^{At} = e^{-t}(I + tN) matches expm at t = 1.3; x(t) = (1+t)e^{-t}.")

# e^A e^B != e^{A+B} when AB != BA (the pair from section 3.9), equality when they commute.
A_n = np.array([[0.0, 1.0], [0.0, 0.0]]); B_n = np.array([[0.0, 0.0], [1.0, 0.0]])
assert not np.allclose(A_n @ B_n, B_n @ A_n)
assert np.allclose(expm(A_n) @ expm(B_n), [[2, 1], [1, 1]], rtol=0, atol=1e-14)
assert np.allclose(expm(A_n + B_n), [[np.cosh(1), np.sinh(1)], [np.sinh(1), np.cosh(1)]], rtol=0, atol=1e-14)
assert not np.allclose(expm(A_n) @ expm(B_n), expm(A_n + B_n))
assert not np.allclose(expm(A_n) @ expm(B_n), expm(B_n) @ expm(A_n))
C1 = rng.standard_normal((3, 3)); C2 = 2 * C1 @ C1 - 0.5 * C1 + np.eye(3)      # a polynomial in C1 commutes with C1
assert np.allclose(C1 @ C2, C2 @ C1) and np.allclose(expm(C1 + C2), expm(C1) @ expm(C2))
print("e^A e^B = [[2,1],[1,1]] but e^{A+B} = [[cosh1, sinh1],[sinh1, cosh1]]: the counterexample holds; "
      "commuting matrices satisfy e^{A+B} = e^A e^B.")

# (3.20): e^{(aI + bJ)t} = e^{at} R(bt); and R(theta) has eigenvalues e^{+-i theta}, eq. (3.17).
J = np.array([[0.0, -1.0], [1.0, 0.0]])
R = lambda th: np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
for a, b, tt in [(-0.3, 2.0, 1.1), (0.5, -1.0, 0.4), (0.0, 1.0, np.pi / 3)]:
    assert np.allclose(expm((a * np.eye(2) + b * J) * tt), np.exp(a * tt) * R(b * tt))
theta = 0.6
assert np.allclose(sorted(np.linalg.eigvals(R(theta)), key=lambda z: z.imag), [np.exp(-1j * theta), np.exp(1j * theta)])
lamR, XR = np.linalg.eig(R(theta))
assert same_direction(XR[:, np.argmax(lamR.imag)], [1, -1j])
print("e^{J theta} = R(theta); R(theta) has eigenvalues e^{+-i theta} with eigenvector (1, -i) for e^{i theta}.")

# %% [markdown]
# ## Part 6: Phase portraits and the trace-determinant plane (section 3.15)

# %%
def classify_2x2(A):
    """Section 3.15: the type of the phase portrait of u' = Au from T = tr A, D = det A."""
    A = np.asarray(A, float); T = np.trace(A); D = np.linalg.det(A); disc = T * T - 4 * D
    if abs(D) < 1e-12:
        return "degenerate (zero eigenvalue)"
    if D < 0:
        return "saddle"
    if disc > 1e-12:
        return "stable node" if T < 0 else "unstable node"
    if disc < -1e-12:
        if abs(T) < 1e-12:
            return "center"
        return "stable spiral" if T < 0 else "unstable spiral"
    return "repeated eigenvalue (star or degenerate node)"


portraits = {
    "saddle": np.array([[1.0, 2.0], [2.0, 1.0]]),                 # lambda = 3, -1
    "stable node": np.array([[-3.0, 1.0], [1.0, -3.0]]),           # lambda = -2, -4
    "stable spiral": A_sp,                                         # lambda = -1 +/- 2i
    "center": np.array([[0.0, -2.0], [1.0, 0.0]]),                # lambda = +/- i sqrt2
}
for name, A in portraits.items():
    assert classify_2x2(A) == name, (name, classify_2x2(A))
    ev = np.linalg.eigvals(A)
    if name == "saddle":
        assert np.isreal(ev).all() and np.prod(ev.real) < 0
    elif name == "stable node":
        assert np.isreal(ev).all() and np.all(ev.real < 0)
    elif name == "stable spiral":
        assert np.all(ev.imag != 0) and np.all(ev.real < 0)
    else:
        assert np.allclose(ev.real, 0) and np.all(ev.imag != 0)
assert classify_2x2(A_cd) == "repeated eigenvalue (star or degenerate node)"
assert classify_2x2(np.array([[1.0, -1.0], [1.0, 1.0]])) == "unstable spiral"   # B5
print("Trace-determinant classification agrees with the eigenvalues for all test matrices.")

fig, axes = plt.subplots(2, 2, figsize=(10, 9.5))
g = np.linspace(-3, 3, 24); GX, GY = np.meshgrid(g, g)
for ax, (name, A) in zip(axes.ravel(), portraits.items()):
    DX = A[0, 0] * GX + A[0, 1] * GY; DY = A[1, 0] * GX + A[1, 1] * GY
    ax.streamplot(GX, GY, DX, DY, color="0.55", density=1.1, linewidth=0.8, arrowsize=0.9)
    lam, X = np.linalg.eig(A)
    if np.isreal(lam).all():
        for j in range(2):
            x = X[:, j].real; s = 4.0 / np.abs(x).max()
            ax.plot([-s * x[0], s * x[0]], [-s * x[1], s * x[1]], lw=2,
                    color="C3" if lam[j].real > 0 else "C0",
                    label=f"eigenvector, lambda = {lam[j].real:.3g}")
    else:
        Uu = solve_linear_ode(A, [2.0, 0.0], np.linspace(0, 8, 800))
        ax.plot(Uu[0], Uu[1], "C3", lw=2, label=f"one trajectory, lambda = {lam[0]:.3g}")
    ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_aspect("equal")
    ax.set_title(f"{name}:  T = {np.trace(A):g},  D = {np.linalg.det(A):g}"); ax.legend(fontsize=8, loc="lower left")
fig.tight_layout(); fig.savefig(FIG / "phase_portraits.png", dpi=120); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 5.5))
Tg = np.linspace(-6, 6, 400)
ax.plot(Tg, Tg ** 2 / 4, "k", lw=1.5, label="D = T^2/4: repeated eigenvalue")
ax.axhline(0, color="k", lw=0.8); ax.axvline(0, color="k", lw=0.8)
ax.fill_between(Tg, Tg ** 2 / 4, 10, where=Tg < 0, color="C0", alpha=.12)
ax.fill_between(Tg, Tg ** 2 / 4, 10, where=Tg > 0, color="C3", alpha=.12)
ax.fill_between(Tg, 0, Tg ** 2 / 4, where=Tg < 0, color="C0", alpha=.25)
ax.fill_between(Tg, 0, Tg ** 2 / 4, where=Tg > 0, color="C3", alpha=.25)
ax.fill_between(Tg, -4, 0, color="C2", alpha=.15)
for txt, xy in [("stable spiral", (-2.2, 6.5)), ("unstable spiral", (1.2, 6.5)), ("stable node", (-5.2, 2.0)),
                ("unstable node", (3.4, 2.0)), ("saddle", (-1.0, -2.5)), ("centers (T = 0)", (0.15, 9.0))]:
    ax.text(*xy, txt, fontsize=9)
for name, A in list(portraits.items()) + [("critically damped (Ex. 4.5)", A_cd)]:
    ax.plot(np.trace(A), np.linalg.det(A), "ko"); ax.annotate(name, (np.trace(A), np.linalg.det(A)),
                                                               textcoords="offset points", xytext=(6, -10), fontsize=8)
ax.set_xlim(-6, 6); ax.set_ylim(-4, 10); ax.set_xlabel("T = trace"); ax.set_ylabel("D = determinant")
ax.set_title("The trace-determinant plane (section 3.15)"); ax.legend(loc="upper left", fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "trace_det_plane.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: Normal modes of two masses and three springs (Example 4.3, section 3.10)
#
# M = diag(1, 2), K = [[2,-1],[-1,2]].  K v = omega^2 M v, eq. (3.15).
# omega^2 = (3 -/+ sqrt3)/2, modes (1, (1+sqrt3)/2) and (1, (1-sqrt3)/2),
# M-orthogonal but not orthogonal. Released from x = (1, 0) at rest:
# x(t) = a1 v1 cos(w1 t) + a2 v2 cos(w2 t) with a1 = (3-sqrt3)/6, a2 = (3+sqrt3)/6.

# %%
M_sp = np.diag([1.0, 2.0]); K_sp = np.array([[2.0, -1.0], [-1.0, 2.0]])
w2, Vm = eigh(K_sp, M_sp)                        # generalized symmetric-definite eigenproblem
assert np.allclose(w2, [(3 - SQ3) / 2, (3 + SQ3) / 2]), w2
omega = np.sqrt(w2)
assert np.allclose(omega, [0.796225, 1.538189], atol=1e-6)
v1 = np.array([1.0, (1 + SQ3) / 2]); v2 = np.array([1.0, (1 - SQ3) / 2])
assert same_direction(Vm[:, 0], v1) and same_direction(Vm[:, 1], v2)
assert np.allclose(K_sp @ v1, w2[0] * M_sp @ v1) and np.allclose(K_sp @ v2, w2[1] * M_sp @ v2)
assert abs(v1 @ M_sp @ v2) < 1e-12 and abs(v1 @ v2 - 0.5) < 1e-12        # M-orthogonal, not orthogonal
assert np.allclose(Vm.T @ M_sp @ Vm, np.eye(2))                           # eigh returns M-orthonormal modes
# Cross-check with the ordinary eigenproblem for M^{-1}K (not symmetric).
MK = np.linalg.solve(M_sp, K_sp)
assert np.allclose(np.trace(MK), 3) and np.allclose(np.linalg.det(MK), 1.5)
assert np.allclose(np.sort(np.linalg.eigvals(MK).real), w2)
# Initial condition x(0) = (1, 0), xdot(0) = 0.
a = np.linalg.solve(np.c_[v1, v2], [1.0, 0.0])
assert np.allclose(a, [(3 - SQ3) / 6, (3 + SQ3) / 6]), a
print(f"Two-mass system: omega = {omega}, modes {v1} and {v2}, v1^T M v2 = {v1 @ M_sp @ v2:.1e}, "
      f"a = {a}.")


def modes_solution(t):
    """x(t) from eq. (3.16) with b_j = 0: sum_j a_j v_j cos(omega_j t). Shape (2, len(t))."""
    return a[0] * np.outer(v1, np.cos(omega[0] * t)) + a[1] * np.outer(v2, np.cos(omega[1] * t))


t_m = np.linspace(0, 40, 2001)
X_modes = modes_solution(t_m)
Minv_K = np.linalg.solve(M_sp, K_sp)
sol = solve_ivp(lambda _t, y: np.r_[y[2:], -Minv_K @ y[:2]], (0, 40), [1, 0, 0, 0],
                t_eval=t_m, rtol=1e-10, atol=1e-12, method="DOP853")
assert np.allclose(sol.y[:2], X_modes, atol=1e-6), "normal-mode solution disagrees with solve_ivp"
# The frequency ratio sqrt(2 + sqrt3) is irrational: the motion is not periodic.
assert np.isclose(omega[1] / omega[0], np.sqrt(2 + SQ3))
print("Normal-mode superposition matches solve_ivp over 40 time units; omega2/omega1 = sqrt(2 + sqrt3).")

fig, ax = plt.subplots(figsize=(10, 3.6))
ax.plot(t_m, X_modes[0], label="x1 (mass 1), eq. (3.16)"); ax.plot(t_m, X_modes[1], label="x2 (mass 2), eq. (3.16)")
ax.plot(t_m[::50], sol.y[0, ::50], "k.", ms=4, label="solve_ivp"); ax.plot(t_m[::50], sol.y[1, ::50], "k.", ms=4)
ax.set_xlabel("t"); ax.set_title("Example 4.3: energy sloshes between the masses (two incommensurate frequencies)")
ax.legend(fontsize=8, ncol=3); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "two_mass_trajectory.png", dpi=120); plt.close(fig)


def draw_frame(ax, x, masses, amp=0.35):
    """Walls at 0 and 3, equilibrium positions 1 and 2, masses drawn at 1 + amp*x1, 2 + amp*x2."""
    pos = np.array([1.0, 2.0]) + amp * np.asarray(x)
    pts = np.r_[0.0, pos, 3.0]
    for p, q in zip(pts[:-1], pts[1:]):                    # springs as zigzags
        s = np.linspace(p, q, 30); z = 0.08 * np.sign(np.sin(np.linspace(0, 12 * np.pi, 30)))
        z[0] = z[-1] = 0; ax.plot(s, z, "0.4", lw=1)
    for p, m in zip(pos, masses):
        ax.plot(p, 0, "o", ms=10 + 6 * m, color="C0" if m == 1 else "C3")
    for wall in (0.0, 3.0):
        ax.plot([wall, wall], [-0.3, 0.3], "k", lw=3)
    ax.set_xlim(-0.1, 3.1); ax.set_ylim(-0.4, 0.4); ax.set_xticks([]); ax.set_yticks([])


n_frames = 8
fig, axes = plt.subplots(2, n_frames, figsize=(14, 3.2))
for row, (vj, wj, label) in enumerate([(v1, omega[0], "mode 1 (slow, in phase)"), (v2, omega[1], "mode 2 (fast, opposed)")]):
    for col in range(n_frames):
        tj = col * (2 * np.pi / wj) / n_frames
        draw_frame(axes[row, col], vj * np.cos(wj * tj), masses=[1, 2])
        axes[row, col].set_title(f"t = {tj:.2f}", fontsize=8)
    axes[row, 0].set_ylabel(f"{label}\nomega = {wj:.3f}", fontsize=8)
fig.suptitle("Example 4.3: one period of each normal mode (blue m = 1, red m = 2)", fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "normal_modes_frames.png", dpi=120); plt.close(fig)

# Three-mass chain (C3 and 'Your turn'): eigenvalues of tridiag(-1, 2, -1) are 2 - 2cos(j pi/4).
T3 = 2 * np.eye(3) - np.eye(3, k=1) - np.eye(3, k=-1)
mu3 = np.linalg.eigvalsh(T3)
assert np.allclose(mu3, [2 - np.sqrt(2), 2, 2 + np.sqrt(2)])
assert np.allclose(mu3, [2 - 2 * np.cos(j * np.pi / 4) for j in (1, 2, 3)])
print(f"Three-mass chain: eigenvalues of tridiag(-1,2,-1) = {mu3} = 2 - 2cos(j pi/4).")

# %% [markdown]
# ## Part 8: Identities, similarity, defectiveness, and the symmetric 3x3 (sections 3.3-3.5, 3.12-3.14)

# %%
for _ in range(50):
    n = rng.integers(2, 7)
    A = rng.standard_normal((n, n))
    lam = np.linalg.eigvals(A)
    # Theorem 3.4: sum = trace, product = det.
    assert np.isclose(lam.sum().real, np.trace(A)) and abs(lam.sum().imag) < 1e-9
    assert np.isclose(np.prod(lam).real, np.linalg.det(A)) and abs(np.prod(lam).imag) < 1e-8
    # Theorem 3.5(b),(c): A^2, A^{-1}, A + cI, A^T.
    key = lambda z: (round(z.real, 6), round(z.imag, 6))
    assert np.allclose(sorted(np.linalg.eigvals(A @ A), key=key), sorted(lam ** 2, key=key))
    assert np.allclose(sorted(np.linalg.eigvals(np.linalg.inv(A)), key=key), sorted(1 / lam, key=key))
    assert np.allclose(sorted(np.linalg.eigvals(A + 2.5 * np.eye(n)), key=key), sorted(lam + 2.5, key=key))
    assert np.allclose(sorted(np.linalg.eigvals(A.T), key=key), sorted(lam, key=key))
    # Theorem 3.5(a): triangular.
    Tri = np.triu(A)
    assert np.allclose(sorted(np.linalg.eigvals(Tri).real), sorted(np.diag(Tri)))
    # Theorem 3.18: similar matrices, same eigenvalues, eigenvectors M^{-1} x.
    Mm = rng.standard_normal((n, n)) + 3 * np.eye(n)
    B = np.linalg.solve(Mm, A @ Mm)
    assert np.allclose(sorted(np.linalg.eigvals(B), key=key), sorted(lam, key=key))
    lamA, XA = np.linalg.eig(A)
    for j in range(n):
        y = np.linalg.solve(Mm, XA[:, j])
        assert np.allclose(B @ y, lamA[j] * y)
    # Theorem 3.17 on a random symmetric matrix: real eigenvalues, orthonormal eigenvectors.
    S = A + A.T
    lamS, QS = np.linalg.eigh(S)
    assert np.allclose(QS.T @ QS, np.eye(n)) and np.allclose(QS @ np.diag(lamS) @ QS.T, S)
    assert np.abs(np.linalg.eigvals(S).imag).max() < 1e-10
print("Trace/det, powers/inverse/shift/transpose, triangular, similarity and symmetric checks passed (50 random matrices).")

# The shear: AM(1) = 2, GM(1) = 1, S^k = [[1,k],[0,1]], e^{St} = e^t [[1,t],[0,1]].
S_sh = np.array([[1.0, 1.0], [0.0, 1.0]])
assert np.allclose(np.linalg.eigvals(S_sh), [1, 1])
assert 2 - np.linalg.matrix_rank(S_sh - np.eye(2)) == 1
for k in (3, 10, 25):
    assert np.allclose(np.linalg.matrix_power(S_sh, k), [[1, k], [0, 1]])
assert np.allclose(expm(S_sh * 0.7), np.exp(0.7) * np.array([[1, 0.7], [0, 1]]))
# Section 3.4 warning: eigenvalues of A + B are not sums (A_n, B_n from Part 5).
assert np.allclose(np.linalg.eigvals(A_n), 0) and np.allclose(np.linalg.eigvals(B_n), 0)
assert np.allclose(sorted(np.linalg.eigvals(A_n + B_n).real), [-1, 1])
assert np.allclose(sorted(np.linalg.eigvals(A_n @ B_n).real), [0, 1])
# ...and a shared eigenbasis is sufficient (B1c) but not necessary: these two do not commute,
# yet, being upper triangular, their sum and product have eigenvalues 1+2, 2+1 and 1*2, 2*1.
A_t = np.array([[1.0, 1.0], [0.0, 2.0]]); B_t = np.diag([2.0, 1.0])
assert not np.allclose(A_t @ B_t, B_t @ A_t)
assert np.allclose(np.linalg.eigvals(A_t + B_t), [3, 3]) and np.allclose(np.linalg.eigvals(A_t @ B_t), [2, 2])
print("Shear: AM = 2, GM = 1, S^k = [[1,k],[0,1]]; eigenvalues of A + B are not lambda + mu in general.")

# B7(d): A_eps = [[3, 1], [0, 3 + eps]] is diagonalizable but X_eps is nearly singular. The top-right
# entry of X Lambda^k X^{-1} is ((3+eps)^k - 3^k)/eps, a difference quotient that loses about
# log10(1/eps) ~ 8 of the 16 digits for eps = 1e-8; matrix_power stays accurate to ~1e-16.
eps = 1e-8
A_eps = np.array([[3.0, 1.0], [0.0, 3.0 + eps]])
lam_e, X_e = np.linalg.eig(A_eps)
err_eig, err_pow = [], []
for kk in (5, 10, 30):
    e_ex = Fraction(eps)                                   # the double eps, exactly
    true = float(((3 + e_ex) ** kk - 3 ** kk) / e_ex)     # exact top-right entry, rounded once
    err_eig.append(abs(((X_e * lam_e ** kk) @ np.linalg.inv(X_e))[0, 1] - true) / true)
    err_pow.append(abs(np.linalg.matrix_power(A_eps, kk)[0, 1] - true) / true)
assert max(err_pow) < 1e-14, err_pow
assert 1e-11 < max(err_eig) < 1e-6, err_eig                  # about half the digits lost, not all
print(f"B7(d), eps = 1e-8: relative error of the top-right entry of A^k (k = 5, 10, 30) via X Lambda^k X^-1: "
      f"{', '.join(f'{e:.1e}' for e in err_eig)}; via matrix_power: {max(err_pow):.1e}.")

# Example 4.6: the symmetric 3x3.
A6 = np.array([[3.0, 1.0, -1.0], [1.0, 3.0, -1.0], [-1.0, -1.0, 5.0]])
lam6, Q6 = np.linalg.eigh(A6)
assert np.allclose(lam6, [2, 3, 6]), lam6
assert np.allclose(np.poly(A6), [1, -11, 36, -36])            # monic char poly: lambda^3 - 11 lambda^2 + 36 lambda - 36
assert np.isclose(np.trace(A6), 11) and np.isclose(np.linalg.det(A6), 36)
for j, x in enumerate([[1, -1, 0], [1, 1, 1], [1, 1, -2]]):
    assert same_direction(Q6[:, j], x), (j, Q6[:, j])
    assert np.allclose(A6 @ np.array(x, float), lam6[j] * np.array(x, float))
Xh = np.array([[1, 1, 1], [-1, 1, 1], [0, 1, -2]], float)      # hand eigenvectors as columns
assert np.allclose(Xh.T @ Xh, np.diag([2, 3, 6]))             # mutually orthogonal (Theorem 3.17b)
assert np.allclose(Q6 @ np.diag(lam6) @ Q6.T, A6)             # A = Q Lambda Q^T (LA-8 preview)
print(f"Example 4.6: eigenvalues {lam6}; NumPy's eigenvectors are the hand ones up to sign/scale; "
      f"last column from eigh = {Q6[:, 2]} vs (1,1,-2)/sqrt6 = {np.array([1, 1, -2]) / np.sqrt(6)}.")

# Numbers from the exercise solutions.
lamA1, XA1 = eig_sorted(A1)                                                     # A1
assert np.allclose(lamA1, [5, 2]) and same_direction(XA1[:, 0], [1, 1]) and same_direction(XA1[:, 1], [1, -2])
U_a2 = np.array([[2.0, 7.0, -3.0], [0.0, -1.0, 4.0], [0.0, 0.0, 5.0]])           # A2
assert np.allclose(U_a2 @ [7, -3, 0], -1 * np.array([7, -3, 0])) and np.allclose(U_a2 @ [5, 6, 9], 5 * np.array([5, 6, 9]))
A_a3 = np.array([[0.8, 0.3], [0.2, 0.7]])                                      # A3
assert np.allclose(np.linalg.matrix_power(A_a3, 2) @ [1, 0], [0.7, 0.3]) and np.allclose(A_a3 @ [0.6, 0.4], [0.6, 0.4])
k5 = 5                                                                          # B3
assert np.allclose(np.linalg.matrix_power(A_sec2, k5), 0.5 * np.array([[3**k5 + 1, 3**k5 - 1], [3**k5 - 1, 3**k5 + 1]]))
assert np.allclose(np.linalg.matrix_power(A_sec2, 40) / 3.0**40, 0.5 * np.ones((2, 2)), atol=1e-12)
A_b5 = np.array([[1.0, -1.0], [1.0, 1.0]])                                     # B5: e^{At}(1,0) = e^t (cos t, sin t)
tb = 0.9
assert np.allclose(expm(A_b5 * tb) @ [1, 0], np.exp(tb) * np.array([np.cos(tb), np.sin(tb)]))
A_b7 = np.array([[3.0, 1.0], [0.0, 3.0]])                                      # B7: A^k and e^{At}
for kk in (2, 7):
    assert np.allclose(np.linalg.matrix_power(A_b7, kk), [[3**kk, kk * 3**(kk - 1)], [0, 3**kk]])
assert np.allclose(expm(A_b7 * tb), np.exp(3 * tb) * np.array([[1, tb], [0, 1]]))
assert abs(phi**10 / np.sqrt(5) - 55.0036) < 1e-4 and round(binet(10)) == F[10] == 55   # Example 4.2 check value
print("Exercise numbers (A1, A2, A3, B3, B5, B7) verified.")

# A^100 two ways (a preview of 'Your turn' 3): repeated multiplication vs X Lambda^100 X^{-1}.
A_big = rng.standard_normal((60, 60)) / np.sqrt(60)
t0 = time.perf_counter(); P1 = np.linalg.matrix_power(A_big, 100); t1 = time.perf_counter()
lamB, XB = np.linalg.eig(A_big); P2 = ((XB * lamB ** 100) @ np.linalg.inv(XB)).real; t2 = time.perf_counter()
assert np.allclose(P1, P2, atol=1e-8 * max(1.0, np.abs(P1).max()))
print(f"A^100 for a random 60x60: matrix_power {1e3*(t1-t0):.2f} ms, eigendecomposition {1e3*(t2-t1):.2f} ms, "
      f"max difference {np.abs(P1-P2).max():.2e}.")

# %% [markdown]
# ## Your turn
#
# 1. **A three-state Markov chain.** Invent a 3x3 column-stochastic matrix for
#    something you care about (weather with "cloudy"; a web graph with three
#    pages; a machine with states working / degraded / broken). Compute its
#    eigenvalues, confirm lambda = 1 and |lambda_i| <= 1 (Theorem 3.11), find
#    the steady state from the eigenvector, and simulate 50,000 steps with
#    `simulate_markov`. Plot ||u_k - u_inf|| on a log scale and read off
#    |lambda_2| from the slope. Then break convergence. A single zero entry is
#    not enough: with the other entries positive, A^2 is strictly positive and
#    Theorem 3.12 still applies. Make the chain periodic (a cyclic permutation
#    of the three states) or reducible (two closed classes), check that a
#    second eigenvalue now has modulus 1, and watch u_k stop converging.
# 2. **A three-mass chain.** Build K = tridiag(-1, 2, -1) for n = 3 with
#    M = I, solve K v = omega^2 v, and compare omega_j with
#    2 sin(j pi / 8) from exercise C3. Reuse `draw_frame` (extend it to three
#    masses and four springs) to draw one period of each of the three modes.
#    Which mode leaves the middle mass still, and why?
# 3. **A^100 two ways, timed.** Repeat the last cell of Part 8 for n = 10, 100,
#    400. For a nearly defective matrix such as [[3, 1], [0, 3 + 1e-9]] (B7d),
#    compare the two answers for A^30 and explain the discrepancy.
# 4. **Inverse iteration as a root polisher.** Run `inverse_iteration` on A6
#    with shifts 2.9, 3.1 and 3.0 + 1e-12 and count steps: the closer the
#    shift, the smaller the ratio |3 - sigma| / |2 - sigma| and the fewer the
#    steps. Now try shift exactly 3.0. A6 has integer entries, so A6 - 3I is
#    *exactly* singular in floating point: `lu_factor` reports a zero pivot and
#    the solve returns inf/NaN. Move the shift by 1e-13 and it works in two
#    steps. Explain why a nearly singular solve is what inverse iteration
#    wants: the huge component of the solution lies along the eigenvector you
#    are after, and normalizing throws away its size.
# 5. (C1) Verify Cayley-Hamilton numerically: `np.poly(A)` gives the monic
#    characteristic polynomial; evaluate it at A with `matrix_power` and
#    check you get the zero matrix.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
