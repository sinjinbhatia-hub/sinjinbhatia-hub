# %% [markdown]
# # CM-6 Project: Normal Modes, the Generalized Eigenproblem, and Modal Superposition
#
# This script builds the module's one algorithm, `normal_modes(K, M)` (Cholesky
# reduction of Theorem 3.4, then a symmetric eigensolver), checks it against
# `scipy.linalg.eigh(K, M)`, and then runs it on every worked system of §4:
# the coupled pendulums and their beats (against the exact nonlinear motion),
# the two unequal masses, CO2, the double pendulum, the loaded string and its
# continuum limit, and the triangle with its degenerate pair. Part 7 does the
# damped, driven frequency response (3.16); Part 8 fits initial conditions by
# M-projection (3.13), zero mode included. Part 9 checks the numbers quoted in
# `solutions.md`. Every number printed in the notes is asserted here.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.linalg as sla
from scipy.integrate import solve_ivp
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
SQ2, SQ3 = np.sqrt(2.0), np.sqrt(3.0)

# %% [markdown]
# ## Part 1: `normal_modes(K, M)` by Cholesky reduction (Theorem 3.4)
#
# K a = w^2 M a (3.5). With M = L L^T, set Kt = L^{-1} K L^{-T} (3.9), symmetric;
# eigh(Kt) gives w^2 and orthonormal y; a = L^{-T} y is M-orthonormal (3.10).
# Collected: A^T M A = I, A^T K A = diag(w^2) (3.11).
#
# **Build it yourself first:** write `normal_modes` and `fit_initial_conditions`
# (Part 8) before reading them. The two classic mistakes: forgetting to
# M-normalize, and dividing by omega = 0 on a zero mode.


# %%
def normal_modes(K, M):
    """Generalized eigenpairs of (K, M), eq. (3.5), via the Cholesky reduction (3.9).

    Returns (w2, A): w2 ascending (omega_j^2, may include ~0 or negatives), A with
    columns a_j satisfying A^T M A = I and A^T K A = diag(w2), eq. (3.11)."""
    K = np.asarray(K, float); M = np.asarray(M, float)
    assert np.allclose(K, K.T) and np.allclose(M, M.T), "K and M must be symmetric"
    L = np.linalg.cholesky(M)                                  # M = L L^T (fails unless M > 0)
    LinvK = sla.solve_triangular(L, K, lower=True)             # L^{-1} K
    Kt = sla.solve_triangular(L, LinvK.T, lower=True).T        # (L^{-1} K) L^{-T}
    Kt = 0.5 * (Kt + Kt.T)                                     # kill round-off asymmetry
    w2, Y = np.linalg.eigh(Kt)                                 # spectral theorem on Kt
    A = sla.solve_triangular(L.T, Y, lower=False)              # a_j = L^{-T} y_j
    return w2, A


def frequencies(w2, rel_tol=1e-9):
    """omega_j = sqrt(omega_j^2) with zero modes thresholded (§3.5, §8): |omega_j^2| < rel_tol * max|omega^2|
    counts as zero (an all-zero spectrum is all zero modes). Returns (omega, is_zero)."""
    w2 = np.asarray(w2, float); scale = np.abs(w2).max()
    zero = np.abs(w2) <= rel_tol * scale
    assert np.all(w2[~zero] > 0), "negative omega^2: unstable equilibrium (Theorem 3.5)"
    return np.sqrt(np.where(zero, 0.0, w2)), zero


def sign_fix(A):
    """Flip each column so its largest-|entry| is positive (eigenvectors are defined up to sign).
    Ties (e.g. an antisymmetric mode, +x and -x) go to the first such entry, so the rule is deterministic."""
    A = np.array(A, float)
    idx = np.argmax(np.abs(A) >= np.abs(A).max(axis=0) * (1 - 1e-8), axis=0)
    s = np.sign(A[idx, np.arange(A.shape[1])]); s[s == 0] = 1
    return A * s


def check_modes(K, M, w2, A, tol=1e-9):
    """Assert (3.5) and (3.11) hold for the returned modes."""
    n = len(w2); scale = max(1.0, np.abs(w2).max())
    assert np.allclose(K @ A, (M @ A) * w2, atol=tol * scale), "(3.5) fails"
    assert np.allclose(A.T @ M @ A, np.eye(n), atol=tol), "A^T M A != I"
    assert np.allclose(A.T @ K @ A, np.diag(w2), atol=tol * scale), "A^T K A != diag(w2)"


def eigenspace_projectors(w2, A, rel_tol=1e-8):
    """Sum_j a_j a_j^T over each cluster of equal w2: unique even when eigenvectors are not."""
    groups, start = [], 0
    for j in range(1, len(w2) + 1):
        if j == len(w2) or abs(w2[j] - w2[start]) > rel_tol * max(1.0, abs(w2[start])):
            groups.append(A[:, start:j] @ A[:, start:j].T); start = j
    return groups


def random_spd(n):
    Q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    return Q @ np.diag(rng.uniform(0.5, 5, n)) @ Q.T


for n in (2, 3, 5, 8):
    for _ in range(5):
        M = random_spd(n); K = rng.standard_normal((n, n)); K = K + K.T   # K indefinite is allowed
        w2, A = normal_modes(K, M); check_modes(K, M, w2, A)
        w2s, As = sla.eigh(K, M)                                          # LAPACK sygv does the same reduction
        assert np.allclose(w2, w2s) and np.allclose(sign_fix(A), sign_fix(As))
        assert np.allclose(np.sort(np.linalg.eigvals(np.linalg.solve(M, K)).real), w2)  # M^{-1}K: same eigenvalues
print("Part 1: normal_modes agrees with scipy.linalg.eigh(K, M) on random pairs; (3.5), (3.11) hold.")

# %% [markdown]
# ## Part 2: Examples 4.1, 4.2, 4.4 — every printed number
#
# 4.1: coupled pendulums, M = m l^2 I, K = [[mgl + kl^2, -kl^2], [-kl^2, mgl + kl^2]].
# 4.2: masses 1, 2 with three unit springs: eigenvectors of M^{-1}K are NOT orthogonal, but are M-orthogonal.
# 4.4: double pendulum, M = ml^2 [[2,1],[1,1]], K = mgl diag(2,1): w^2 = (2 -+ sqrt2) g/l.
# The Hessians (3.1) are also computed symbolically from the exact T and V.

# %%
m_p, l_p, g, k_p, D_p = 1.0, 1.0, 9.81, 1.0, 2.0      # Example 4.1 (D = pivot spacing; drops out of K)
th1, th2 = sp.symbols("theta1 theta2")
dx = D_p + l_p * (sp.sin(th2) - sp.sin(th1)); dy = l_p * (sp.cos(th1) - sp.cos(th2))
V_exact = m_p * g * l_p * (2 - sp.cos(th1) - sp.cos(th2)) + sp.Rational(1, 2) * k_p * (sp.sqrt(dx**2 + dy**2) - D_p) ** 2
K41_sym = sp.hessian(V_exact, (th1, th2)).subs({th1: 0, th2: 0})
K41 = np.array(K41_sym.evalf(), float); M41 = m_p * l_p**2 * np.eye(2)
assert np.allclose(K41, [[10.81, -1], [-1, 10.81]])
w2, A41 = normal_modes(K41, M41); check_modes(K41, M41, w2, A41)
om41 = np.sqrt(w2)
assert np.allclose(w2, [9.81, 11.81]) and np.allclose(om41, [3.1321, 3.4366], atol=5e-5)
assert np.allclose(np.abs(A41), 1 / SQ2) and np.allclose(sign_fix(A41)[:, 0], [1, 1] / SQ2) and np.allclose(sign_fix(A41)[:, 1], [1, -1] / SQ2)
om_bar, eps = om41.mean(), (om41[1] - om41[0]) / 2
T_transfer, T_beat = np.pi / (om41[1] - om41[0]), 2 * np.pi / (om41[1] - om41[0])
assert np.isclose(om_bar, 3.2843, atol=5e-5) and np.isclose(eps, 0.15224, atol=5e-6)
assert np.isclose(om41[1] - om41[0], 0.30448, atol=5e-6)
assert np.isclose(T_transfer, 10.318, atol=5e-4) and np.isclose(T_beat, 20.636, atol=5e-4)
T_beat_weak = 2 * np.pi * m_p * om41[0] / k_p
assert np.isclose(T_beat_weak, 19.68, atol=5e-3) and np.isclose((T_beat - T_beat_weak) / T_beat, 0.046, atol=5e-4)   # the 4.6% gap
print(f"Example 4.1: omega = {om41.round(4)}, omega_bar = {om_bar:.4f}, eps = {eps:.5f}, "
      f"T_transfer = {T_transfer:.3f} s, T_beat = {T_beat:.3f} s (weak-coupling approx {T_beat_weak:.2f} s)")

M42 = np.diag([1.0, 2.0]); K42 = np.array([[2.0, -1], [-1, 2]])
w2, A42 = normal_modes(K42, M42); check_modes(K42, M42, w2, A42)
assert np.allclose(w2, [(3 - SQ3) / 2, (3 + SQ3) / 2]) and np.allclose(w2, [0.63397, 2.36603], atol=5e-6)
assert np.allclose(np.sqrt(w2), [0.79623, 1.53819], atol=5e-6)
S42 = A42 / A42[0]                                   # shapes with first component 1
assert np.allclose(S42[1], [(1 + SQ3) / 2, (1 - SQ3) / 2]) and np.allclose(S42[1], [1.36603, -0.36603], atol=5e-6)
assert np.isclose(S42[:, 0] @ S42[:, 1], 0.5) and np.isclose(S42[:, 0] @ M42 @ S42[:, 1], 0.0)   # dot 1/2, M-dot 0
assert np.allclose([S42[:, 0] @ M42 @ S42[:, 0], S42[:, 1] @ M42 @ S42[:, 1]], [4.73205, 1.26795], atol=5e-6)
assert np.allclose(sign_fix(A42), [[0.45970, 0.88807], [0.62796, -0.32506]], atol=5e-6)
MinvK = np.linalg.solve(M42, K42)
assert np.allclose(MinvK, [[2, -1], [-0.5, 1]]) and not np.allclose(MinvK, MinvK.T)
lam, V = np.linalg.eig(MinvK); V = V / V[0]          # eig of the nonsymmetric M^{-1}K: same eigenvalues, same shapes
assert np.allclose(np.sort(lam.real), w2) and np.allclose(np.sort(V[1]), np.sort(S42[1]))
print(f"Example 4.2: omega^2 = {w2.round(5)}; shapes {S42[1].round(5)}; dot = {S42[:, 0] @ S42[:, 1]:.3f}, "
      f"M-dot = {S42[:, 0] @ M42 @ S42[:, 1]:.1e}; A = {sign_fix(A42).round(5).tolist()}")

td1, td2 = sp.symbols("thetadot1 thetadot2"); m, l, gs = sp.symbols("m l g", positive=True)
T_dp = sp.Rational(1, 2) * m * l**2 * (2 * td1**2 + td2**2 + 2 * td1 * td2 * sp.cos(th1 - th2))
V_dp = -m * gs * l * (2 * sp.cos(th1) + sp.cos(th2))
M_dp = sp.hessian(T_dp, (td1, td2)).subs({th1: 0, th2: 0}); K_dp = sp.hessian(V_dp, (th1, th2)).subs({th1: 0, th2: 0})
assert M_dp == m * l**2 * sp.Matrix([[2, 1], [1, 1]]) and K_dp == m * gs * l * sp.diag(2, 1)
lam_s = sp.symbols("lambda")
char_poly = sp.expand((K_dp - lam_s * (gs / l) * M_dp).det() / (m**2 * l**2 * gs**2))  # lambda = omega^2 l/g
assert sp.simplify(char_poly - (lam_s**2 - 4 * lam_s + 2)) == 0                         # (3.6) for Example 4.4
assert set(sp.solve(char_poly, lam_s)) == {2 - sp.sqrt(2), 2 + sp.sqrt(2)}
M44 = np.array([[2.0, 1], [1, 1]]); K44 = g * np.diag([2.0, 1])            # m = l = 1
w2, A44 = normal_modes(K44, M44); check_modes(K44, M44, w2, A44)
assert np.allclose(w2 / g, [2 - SQ2, 2 + SQ2]) and np.allclose(w2 / g, [0.58579, 3.41421], atol=5e-6)
assert np.allclose(np.sqrt(w2), [2.3972, 5.7874], atol=5e-5) and np.isclose(np.sqrt(w2[1] / w2[0]), 1 + SQ2)
S44 = A44 / A44[0]
assert np.allclose(S44[1], [SQ2, -SQ2]) and np.isclose(S44[:, 0] @ S44[:, 1], -1.0)
assert np.allclose([S44[:, 0] @ M44 @ S44[:, 0], S44[:, 1] @ M44 @ S44[:, 1]], [4 + 2 * SQ2, 4 - 2 * SQ2])
assert np.allclose([4 + 2 * SQ2, 4 - 2 * SQ2], [6.8284, 1.1716], atol=5e-5)
print(f"Example 4.4: omega^2 l/g = {(w2 / g).round(5)}, omega = {np.sqrt(w2).round(4)}, shapes a2/a1 = {S44[1].round(5)}")

# %% [markdown]
# ## Part 3: Beats — exact nonlinear pendulums vs the linear modal solution (4.1)
#
# The exact equations m l^2 theta_i'' = -dV/dtheta_i (full sin, exact spring geometry)
# are derived symbolically and integrated; the overlay is theta_1 = theta_0 cos(eps t) cos(wbar t),
# theta_2 = theta_0 sin(eps t) sin(wbar t). Small amplitude: agreement to about 2% of theta_0
# over two beat periods. Large amplitude: the overlay is meaningless (§8).

# %%
dV_fun = sp.lambdify((th1, th2), [sp.diff(V_exact, th1), sp.diff(V_exact, th2)], "numpy")


def pendulum_rhs(t, y):
    """Exact coupled pendulums: y = (theta1, theta2, theta1dot, theta2dot)."""
    f1, f2 = dV_fun(y[0], y[1])
    return [y[2], y[3], -f1 / (m_p * l_p**2), -f2 / (m_p * l_p**2)]


def beats_linear(t, th0):
    """Eq. (4.1): the two-mode solution for q0 = (theta0, 0), v0 = 0."""
    return th0 * np.cos(eps * t) * np.cos(om_bar * t), th0 * np.sin(eps * t) * np.sin(om_bar * t)


t_b = np.linspace(0, 2 * T_beat, 4001)
fig, axes = plt.subplots(2, 2, figsize=(12, 6), sharex=True)
beat_err = {}
for row, th0 in enumerate([0.05, 1.0]):
    sol = solve_ivp(pendulum_rhs, (0, 2 * T_beat), [th0, 0, 0, 0], t_eval=t_b, method="DOP853", rtol=1e-10, atol=1e-12)
    lin1, lin2 = beats_linear(t_b, th0)
    beat_err[th0] = max(np.abs(sol.y[0] - lin1).max(), np.abs(sol.y[1] - lin2).max()) / th0
    for col, (num, lin, env, name) in enumerate([(sol.y[0], lin1, np.abs(np.cos(eps * t_b)), "theta_1"),
                                                 (sol.y[1], lin2, np.abs(np.sin(eps * t_b)), "theta_2")]):
        ax = axes[row, col]
        ax.plot(t_b, num / th0, "k", lw=0.8, label="exact (solve_ivp)")
        ax.plot(t_b, lin / th0, "C1", lw=0.8, alpha=0.8, label="linear modes (4.1)")
        ax.plot(t_b, env, "C0--", lw=1); ax.plot(t_b, -env, "C0--", lw=1, label="envelope")
        ax.set_title(f"{name}/theta_0, theta_0 = {th0} rad; max deviation {beat_err[th0]:.3f} theta_0", fontsize=9)
        ax.grid(alpha=.3)
axes[0, 0].legend(fontsize=7, loc="upper right"); axes[1, 0].set_xlabel("t [s]"); axes[1, 1].set_xlabel("t [s]")
fig.tight_layout(); fig.savefig(FIG / "beats.png", dpi=120); plt.close(fig)
assert beat_err[0.05] < 0.0205, beat_err         # notes: "agree to 2% of theta_0" (measured 0.0198)
assert beat_err[1.0] > 0.5, beat_err             # notes: "meaningless"
sol = solve_ivp(pendulum_rhs, (0, 2 * T_beat), [0.3, 0, 0, 0], t_eval=t_b, method="DOP853", rtol=1e-10, atol=1e-12)
lin03 = beats_linear(t_b, 0.3)[0]
err03 = np.abs(sol.y[0] - lin03).max() / 0.3
win = int(round(len(t_b) * (2 * np.pi / om_bar) / (2 * T_beat)))   # one carrier period of samples
env = lambda x: np.array([np.abs(x[i:i + win]).max() for i in range(0, len(x) - win, win)])
env_err03 = np.abs(env(sol.y[0]) - env(lin03)).max() / 0.3
shift03 = om41[0] * 0.3**2 / 16                 # §8: anharmonic frequency shift omega_0 theta_0^2 / 16
assert 0.055 < shift03 / (om41[1] - om41[0]) < 0.065 and 0.6 < shift03 * 2 * T_beat < 0.8
assert err03 > 0.5 and env_err03 < 0.15         # §8: at 0.3 rad the carriers are out of step, the envelopes still roughly agree
print(f"Part 3: max |exact - linear|/theta_0 = {beat_err[0.05]:.4f} at 0.05 rad, {err03:.3f} at 0.3 rad "
      f"(envelopes differ by {env_err03:.3f}), {beat_err[1.0]:.3f} at 1.0 rad; "
      f"at 0.3 rad the anharmonic shift is {100 * shift03 / (om41[1] - om41[0]):.1f}% of the splitting, {shift03 * 2 * T_beat:.2f} rad of drift over 2 T_beat")

# %% [markdown]
# ## Part 4: CO2 (Example 4.3) — zero mode, symmetric and antisymmetric stretch
#
# K = B^T C B with B the two-spring difference matrix; M = diag(16, 12, 16); k = 1.
# w_s^2 = k/m = 0.0625, w_a^2 = k/m + 2k/M_c = 11/48; ratio sqrt(11/3) = 1.9149.

# %%
def spring_K(B, c):
    """Stiffness K = B^T C B of a spring network (§3.1): rows of B are first-order elongations."""
    return B.T @ (np.asarray(c, float)[:, None] * B)


m_O, M_C, k_co2 = 16.0, 12.0, 1.0
B_co2 = np.array([[-1.0, 1, 0], [0, -1, 1]]); K_co2 = spring_K(B_co2, [k_co2, k_co2]); M_co2 = np.diag([m_O, M_C, m_O])
assert np.allclose(K_co2, [[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
w2_co2, A_co2 = normal_modes(K_co2, M_co2); check_modes(K_co2, M_co2, w2_co2, A_co2)
om_co2, zero_co2 = frequencies(w2_co2)
assert zero_co2.tolist() == [True, False, False] and np.allclose(w2_co2[1:], [1 / 16, 11 / 48])
assert np.allclose(w2_co2[1:], [0.0625, 0.229167], atol=5e-7) and np.allclose(om_co2[1:], [0.25, 0.47871], atol=5e-6)
assert np.isclose(om_co2[2] / om_co2[1], np.sqrt(11 / 3)) and np.isclose(np.sqrt(11 / 3), 1.9149, atol=5e-5)
shapes_co2 = np.array([[1, 1, 1], [1, 0, -1], [1, -2 * m_O / M_C, 1]], float).T        # a_0, a_s, a_a (unnormalized)
norms_co2 = np.sqrt(np.einsum("ij,ij->j", shapes_co2, M_co2 @ shapes_co2))
assert np.allclose(norms_co2**2, [44, 32, 117.33], atol=5e-3) and np.isclose(norms_co2[2], 10.832, atol=5e-4)
assert np.allclose(sign_fix(A_co2), sign_fix(shapes_co2 / norms_co2))
assert np.allclose(M_co2.diagonal() @ A_co2[:, 1:], 0)                               # vibrations do not move the centre of mass
S_co2 = -np.eye(3)[::-1]                                                              # reflection through C: (x1,x2,x3) -> (-x3,-x2,-x1)
assert np.allclose(S_co2.T @ K_co2 @ S_co2, K_co2) and np.allclose(S_co2.T @ M_co2 @ S_co2, M_co2)
assert np.allclose(S_co2 @ shapes_co2, shapes_co2 * [-1, 1, -1])                       # zero mode odd, symmetric stretch even, antisymmetric stretch odd
assert np.isclose(1e4 / 2349, 4.3, atol=0.05) and np.isclose(1e4 / 667, 15, atol=0.05)  # band positions in micrometres
assert np.isclose(2349 / 1388, 1.69, atol=5e-3) and np.isclose(1.9149 / (2349 / 1388), 1.13, atol=5e-3)   # measured ratio 1.69; model 13% high
print(f"CO2: omega^2 = {w2_co2.round(6)}, omega_a/omega_s = {om_co2[2] / om_co2[1]:.4f}, a_a^T M a_a = {norms_co2[2]**2:.2f}")

fig, axes = plt.subplots(1, 3, figsize=(12, 2.8))
for j, (ax, name) in enumerate(zip(axes, ["zero mode (translation)", "symmetric stretch", "antisymmetric stretch"])):
    a = sign_fix(A_co2)[:, j]; a = a / np.abs(a).max() * 0.45
    ax.scatter([-1, 0, 1], [0, 0, 0], s=[600, 400, 600], c=["C3", "0.3", "C3"], zorder=3)
    ax.quiver([-1, 0, 1], [0, 0, 0], a, [0, 0, 0], angles="xy", scale_units="xy", scale=1, width=0.012, zorder=4)
    ax.set_xlim(-1.8, 1.8); ax.set_ylim(-0.6, 0.6); ax.set_yticks([]); ax.set_title(f"{name}: omega^2 = {w2_co2[j]:.4f} k", fontsize=9)
fig.tight_layout(); fig.savefig(FIG / "co2_modes.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5: N masses on a string (Example 4.5) and the continuum limit
#
# K = (tau/d) tridiag(-1, 2, -1), M = m I. Closed form (4.5): a_k^(j) = sin(jk pi/(N+1)),
# omega_j = 2 omega_0 sin(j pi / 2(N+1)), omega_0^2 = tau/(m d). Dispersion (4.7):
# omega(kappa) = (2c/d) sin(kappa d/2) -> c kappa for kappa d << 1, with c = sqrt(tau/mu).

# %%
def string_KM(N, tau=1.0, L=1.0, mu=1.0):
    """Loaded string: N masses m = mu d at spacing d = L/(N+1). Returns K, M, d, m (Derivation 5)."""
    d = L / (N + 1); m = mu * d
    B = np.zeros((N + 1, N))                                    # (N+1) string segments, N masses
    for s in range(N + 1):
        if s < N: B[s, s] += 1.0
        if s > 0: B[s, s - 1] -= 1.0
    K = spring_K(B, np.full(N + 1, tau / d))
    assert np.allclose(K, (tau / d) * (2 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)))
    return K, m * np.eye(N), d, m


def string_closed_form(N, tau=1.0, L=1.0, mu=1.0):
    """Eq. (4.5): omega_j and M-normalized shapes sqrt(2/(m(N+1))) sin(jk pi/(N+1))."""
    d = L / (N + 1); m = mu * d; om0 = np.sqrt(tau / (m * d)); j = np.arange(1, N + 1); k = np.arange(1, N + 1)
    om = 2 * om0 * np.sin(j * np.pi / (2 * (N + 1)))
    A = np.sqrt(2 / (m * (N + 1))) * np.sin(np.outer(k, j) * np.pi / (N + 1))
    return om, A, om0


for N in (3, 5, 20, 100):
    K, M, d, m = string_KM(N)
    w2, A = normal_modes(K, M); check_modes(K, M, w2, A, tol=1e-7)
    om_cf, A_cf, om0 = string_closed_form(N)
    assert np.allclose(np.sqrt(w2), om_cf, rtol=1e-10) and np.allclose(sign_fix(A), sign_fix(A_cf), atol=1e-7)
    if N == 3: assert np.allclose(w2 / om0**2, [2 - SQ2, 2, 2 + SQ2])
    if N == 5:
        assert np.allclose(om_cf / om0, [0.5176, 1, 1.4142, 1.7321, 1.9319], atol=5e-5)
        fig, axes = plt.subplots(1, 5, figsize=(14, 2.6))
        xs = np.linspace(0, 1, 200)
        for j, ax in enumerate(axes, 1):
            a = sign_fix(A)[:, j - 1]; scale = np.sqrt(2 / (m * (N + 1)))
            ax.plot(xs, np.sin(j * np.pi * xs), "C0", lw=1, label="sin(j pi x/L)")
            ax.plot(np.arange(1, N + 1) * d, a / scale, "ko", label="eigh shape")
            ax.axhline(0, color="0.7", lw=0.5); ax.set_title(f"j={j}: omega/omega_0 = {om_cf[j - 1] / om0:.4f}", fontsize=9)
        axes[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(FIG / "string_modes_N5.png", dpi=120); plt.close(fig)
    if N == 20:
        om_snap, A_snap, om0_snap = om_cf, sign_fix(A), om0; d_snap = d
    if N == 100:
        c = 1.0                                                  # sqrt(tau/mu)
        kappa = np.arange(1, N + 1) * np.pi / 1.0
        assert np.allclose(np.sqrt(w2), 2 * c / d * np.sin(kappa * d / 2))                  # (4.7)
        assert np.isclose(np.sqrt(w2[0]) / (c * kappa[0]), np.sinc(kappa[0] * d / (2 * np.pi)))  # omega_1 = c kappa_1 (1 - (kappa d)^2/24 ...)
        assert abs(np.sqrt(w2[0]) / (c * kappa[0]) - 1) < 1e-4 and np.isclose(np.sqrt(w2[-1]) / (2 * om0), np.sin(50 * np.pi / 101))
        fig, ax = plt.subplots(figsize=(6, 4))
        kk = np.linspace(0, np.pi / d, 400)
        ax.plot(kk * d, 2 * c / d * np.sin(kk * d / 2) * d / c, "C0", label="chain: (2c/d) sin(kappa d/2)")
        ax.plot(kk * d, kk * d, "C1--", label="continuum: omega = c kappa")
        ax.plot(kappa * d, np.sqrt(w2) * d / c, "k.", ms=4, label=f"eigh, N = {N}")
        ax.set_xlabel("kappa d"); ax.set_ylabel("omega d / c"); ax.set_title("Dispersion relation (4.7): linear at long wavelength, saturates at kappa = pi/d")
        ax.legend(fontsize=8); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(FIG / "string_dispersion.png", dpi=120); plt.close(fig)
print("Part 5: string frequencies and shapes match (4.5) for N = 3, 5, 20, 100; dispersion (4.7) checked.")

phases = [0, 1 / 8, 1 / 4, 3 / 8, 1 / 2]                         # fractions of each mode's own period
fig, axes = plt.subplots(3, 5, figsize=(14, 6), sharex=True, sharey=True)
x_nodes = np.arange(0, 22) * d_snap
for j in range(3):
    a = np.r_[0, A_snap[:, j], 0] / np.abs(A_snap[:, j]).max()
    for p, frac in enumerate(phases):
        ax = axes[j, p]; ax.plot(x_nodes, a * np.cos(2 * np.pi * frac), "k.-", lw=0.8, ms=4)
        ax.set_ylim(-1.2, 1.2); ax.grid(alpha=.3)
        if j == 0: ax.set_title(f"t = {frac:g} T_j", fontsize=9)
    axes[j, 0].set_ylabel(f"mode {j + 1}\nomega/omega_0 = {om_snap[j] / om0_snap:.3f}", fontsize=8)
fig.suptitle("Loaded string, N = 20: snapshots of the first three modes over half a period"); fig.tight_layout()
fig.savefig(FIG / "string_snapshots.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 6: Three masses on an equilateral triangle (Example 4.6) — degeneracy
#
# K = k B^T B assembled from geometry; M = m I_6. Spectrum (0,0,0, 3/2, 3/2, 3) k/m.
# The computed shear pair is some rotation of the "textbook" pair; the computed zero-mode
# triple is some rotation of (translation x, translation y, rigid rotation).

# %%
k_t, m_t = 1.0, 1.0
ang = np.deg2rad([90.0, 210.0, 330.0]); R_t = np.c_[np.cos(ang), np.sin(ang)]
B_t = np.zeros((3, 6))
for r, (i, j) in enumerate([(0, 1), (1, 2), (2, 0)]):
    nhat = (R_t[j] - R_t[i]) / np.linalg.norm(R_t[j] - R_t[i])
    B_t[r, 2 * i:2 * i + 2] = -nhat; B_t[r, 2 * j:2 * j + 2] = nhat
K_t = spring_K(B_t, [k_t] * 3); M_t = m_t * np.eye(6)
assert np.linalg.matrix_rank(B_t) == 3 and np.isclose(np.trace(K_t), 6 * k_t)
w2_t, A_t = normal_modes(K_t, M_t); check_modes(K_t, M_t, w2_t, A_t)
om_t, zero_t = frequencies(w2_t)
assert zero_t.sum() == 3 and np.allclose(w2_t, [0, 0, 0, 1.5, 1.5, 3.0], atol=1e-9)
rigid = np.array([[1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1], np.c_[-np.sin(ang), np.cos(ang)].ravel()], float).T
assert np.allclose(K_t @ rigid, 0)                                            # translations and rotation stretch no spring
breathing = R_t.ravel(); assert np.allclose(K_t @ breathing, 3 * k_t * breathing)   # omega_b^2 = 3k/m
# textbook shear pair: push vertex 1 outward, project out the four known modes (M-Gram-Schmidt), rotate by 120 deg
known = np.c_[rigid, breathing]; Q_known, _ = np.linalg.qr(known)
proj_out = lambda v: v - Q_known @ (Q_known.T @ v)
s1 = proj_out(np.array([0, 1, 0, 0, 0, 0.0])); s1 /= np.sqrt(s1 @ M_t @ s1)
Rot = np.array([[np.cos(2 * np.pi / 3), -np.sin(2 * np.pi / 3)], [np.sin(2 * np.pi / 3), np.cos(2 * np.pi / 3)]])
S_t = np.zeros((6, 6))                                                        # rotate triangle by 120 deg: mass i -> i+1, displacement rotated
for i in range(3):
    S_t[2 * ((i + 1) % 3):2 * ((i + 1) % 3) + 2, 2 * i:2 * i + 2] = Rot
assert np.allclose(S_t.T @ K_t @ S_t, K_t) and np.allclose(S_t.T @ M_t @ S_t, M_t)   # a symmetry in the sense of §3.6
s2 = S_t @ s1; s2 -= (s1 @ M_t @ s2) * s1; s2 /= np.sqrt(s2 @ M_t @ s2)
shear_text = np.c_[s1, s2]
assert np.allclose(K_t @ shear_text, 1.5 * k_t * shear_text)                  # both are modes at 3k/2m
P_text, P_comp = shear_text @ shear_text.T, A_t[:, 3:5] @ A_t[:, 3:5].T
assert np.allclose(P_text, P_comp)                                            # same eigenspace ...
R_pair = shear_text.T @ M_t @ A_t[:, 3:5]
assert np.allclose(R_pair.T @ R_pair, np.eye(2))                              # ... computed pair = textbook pair times an orthogonal 2x2
assert np.allclose(eigenspace_projectors(w2_t, A_t)[0], rigid @ np.linalg.inv(rigid.T @ M_t @ rigid) @ rigid.T)  # zero-mode space = rigid motions
print(f"Triangle: omega^2 m/k = {w2_t.round(6)}; computed shear pair = textbook pair times R with R^T R = I, "
      f"R = {R_pair.round(3).tolist()}; rotation angle {np.degrees(np.arctan2(R_pair[1, 0], R_pair[0, 0])):.1f} deg")

fig, axes = plt.subplots(2, 3, figsize=(11, 7))
labels = ["zero mode", "zero mode", "zero mode", "shear (degenerate)", "shear (degenerate)", "breathing"]
tri = np.r_[R_t, R_t[:1]]
for j, ax in enumerate(axes.ravel()):
    a = sign_fix(A_t)[:, j].reshape(3, 2); a = 0.45 * a / np.abs(a).max()
    ax.plot(tri[:, 0], tri[:, 1], "0.6", lw=1); ax.scatter(R_t[:, 0], R_t[:, 1], s=120, c="k", zorder=3)
    ax.quiver(R_t[:, 0], R_t[:, 1], a[:, 0], a[:, 1], angles="xy", scale_units="xy", scale=1, color="C3", width=0.012, zorder=4)
    ax.set_aspect("equal"); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.3, 1.6); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{labels[j]}: omega^2 = {abs(w2_t[j]):.3f} k/m", fontsize=9)
fig.suptitle("Triangle modes as eigh returned them (degenerate subspaces in an arbitrary basis)")
fig.tight_layout(); fig.savefig(FIG / "triangle_modes.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: Frequency response of a damped, driven 3-DOF chain (3.16)
#
# Wall-k-m-k-m-k-m-k-wall, Rayleigh damping C = alpha M + beta K (Theorem 3.6).
# (3.16) summed over modes is compared with the direct solve of (K - w^2 M + i w C) q = F.
# A force at mass 1 shows three peaks; a force at mass 2 misses mode 2, whose shape has a node there.

# %%
def receptance(w2, A, zeta, F, omega):
    """Eq. (3.16): q_hat(omega) = sum_j a_j (a_j^T F) / (omega_j^2 - omega^2 + 2 i zeta_j omega_j omega)."""
    om_j = np.sqrt(w2); f = A.T @ F
    denom = w2[None, :] - omega[:, None]**2 + 2j * zeta[None, :] * om_j[None, :] * omega[:, None]
    return (f[None, :] / denom) @ A.T                              # shape (n_omega, n)


K3 = np.array([[2.0, -1, 0], [-1, 2, -1], [0, -1, 2]]); M3 = np.eye(3)
w2_3, A3 = normal_modes(K3, M3); check_modes(K3, M3, w2_3, A3); om3 = np.sqrt(w2_3)
assert np.allclose(w2_3, [2 - SQ2, 2, 2 + SQ2])
alpha, beta = 0.02, 0.01
C3 = alpha * M3 + beta * K3
zeta3 = (alpha + beta * w2_3) / (2 * om3)                          # (3.15): 2 zeta_j omega_j = alpha + beta omega_j^2
assert np.allclose(A3.T @ C3 @ A3, np.diag(2 * zeta3 * om3))      # Rayleigh damping decouples (Theorem 3.6)
omega_sweep = np.linspace(0.05, 2.5, 6000)
fig, ax = plt.subplots(figsize=(7, 4)); peaks_at = {}
for F, color, name in [(np.array([1.0, 0, 0]), "C0", "force on mass 1"), (np.array([0, 1.0, 0]), "C1", "force on mass 2")]:
    Q_modal = receptance(w2_3, A3, zeta3, F, omega_sweep)
    Q_direct = np.array([np.linalg.solve(K3 - w**2 * M3 + 1j * w * C3, F) for w in omega_sweep])
    assert np.allclose(Q_modal, Q_direct, rtol=1e-9, atol=1e-12)
    amp = np.abs(Q_modal[:, 0])                                     # response of mass 1
    peaks = np.where((amp[1:-1] > amp[:-2]) & (amp[1:-1] > amp[2:]))[0] + 1
    excited = np.where(np.abs(A3.T @ F) > 1e-12)[0]
    assert len(peaks) == len(excited) and np.allclose(omega_sweep[peaks], om3[excited], rtol=0.01)
    peaks_at[name] = omega_sweep[peaks]
    ax.semilogy(omega_sweep, amp, color, lw=1, label=f"|q_1(omega)|, {name}: {len(peaks)} peaks")
    print(f"Part 7 ({name}): peaks at omega = {omega_sweep[peaks].round(3)}, modes {om3.round(3)}; "
          f"a_j^T F = {(A3.T @ F).round(4)}")
for w in om3: ax.axvline(w, color="0.6", ls="--", lw=0.8)
ax.set_xlabel("drive frequency omega"); ax.set_ylabel("|q_1|"); ax.set_title("Receptance (3.16) vs direct solve; dashed = omega_j"); ax.legend(fontsize=8); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "frequency_response.png", dpi=120); plt.close(fig)
assert np.isclose(np.abs(A3[1, 1]), 0.0, atol=1e-12)              # mode 2 has a node at the middle mass
assert np.allclose(peaks_at["force on mass 2"], [0.765, 1.847], atol=5e-4)   # B8(c): the two peaks quoted in solutions.md

# %% [markdown]
# ## Part 8: `fit_initial_conditions` by M-projection (3.13), zero mode included
#
# xi_j(0) = <a_j, q0>_M, xidot_j(0) = <a_j, v0>_M; q(t) = sum_j a_j [xi_j(0) cos w_j t + xidot_j(0) sin(w_j t)/w_j],
# with sin(w t)/w -> t on a zero mode. Checked against direct integration of M q'' = -K q for CO2 given a push.

# %%
def fit_initial_conditions(q0, v0, modes):
    """Modal solution (3.13) for M q'' = -K q with q(0) = q0, q'(0) = v0.

    modes = (w2, A, M) from normal_modes. Returns (q, xi0, xidot0) where q(t) evaluates (3.13)
    at an array of times (shape (len(t), n)); zero modes drift linearly (§3.5)."""
    w2, A, M = modes; om, zero = frequencies(w2)
    xi0 = A.T @ M @ np.asarray(q0, float); xidot0 = A.T @ M @ np.asarray(v0, float)

    def q(t):
        t = np.atleast_1d(np.asarray(t, float))[:, None]
        wt = t * om[None, :]
        sinc = np.where(zero[None, :], t, np.sin(wt) / np.where(zero, 1.0, om)[None, :])
        return (xi0[None, :] * np.cos(wt) + xidot0[None, :] * sinc) @ A.T

    return q, xi0, xidot0


q0 = np.array([0.1, 0.0, -0.05]); v0 = np.array([0.3, 0.1, 0.0])        # total momentum 16*0.3 + 12*0.1 = 6: it drifts
q_modal, xi0, xidot0 = fit_initial_conditions(q0, v0, (w2_co2, A_co2, M_co2))
t_c = np.linspace(0, 60, 1201)
sol = solve_ivp(lambda t, y: np.r_[y[3:], -np.linalg.solve(M_co2, K_co2 @ y[:3])], (0, 60), np.r_[q0, v0],
                t_eval=t_c, method="DOP853", rtol=1e-11, atol=1e-13)
Qm = q_modal(t_c)
assert np.abs(Qm - sol.y[:3].T).max() < 1e-7
assert np.allclose(Qm[0], q0) and np.allclose((q_modal(1e-6)[0] - q_modal(-1e-6)[0]) / 2e-6, v0, atol=1e-6)
P_tot, M_tot = M_co2.diagonal() @ v0, M_co2.trace()
x_cm = (M_co2.diagonal() @ Qm.T) / M_tot
assert np.allclose(x_cm, x_cm[0] + P_tot / M_tot * t_c)                    # the zero mode carries the centre of mass
assert np.isclose(xidot0[0] * A_co2[0, 0], P_tot / M_tot)                 # a_0 <a_0, v0>_M = v_cm (1,1,1)
print(f"Part 8: modal solution matches solve_ivp to {np.abs(Qm - sol.y[:3].T).max():.1e}; v_cm = {P_tot / M_tot:.5f}")
fig, ax = plt.subplots(figsize=(7, 3.5))
for i in range(3):
    ax.plot(t_c, sol.y[i], "k", lw=0.8); ax.plot(t_c, Qm[:, i], f"C{i}", lw=1.5, alpha=.6, label=f"x_{i + 1}: modal (color), solve_ivp (black)")
ax.plot(t_c, x_cm, "0.5", ls="--", label="centre of mass (zero mode)"); ax.legend(fontsize=7); ax.set_xlabel("t"); ax.grid(alpha=.3)
ax.set_title("CO2 given a push: vibration on top of drift, (3.13) with the zero mode"); fig.tight_layout(); fig.savefig(FIG / "co2_drift.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 9: Numbers quoted in `solutions.md`
#
# Each block is labelled with its exercise. Nothing new is derived here; this makes the
# solutions file checkable.

# %%
# A1: m = (1, 2), springs (1, 2, 2): K = [[3,-2],[-2,4]], omega^2 = 1, 4, shapes (1,1)/sqrt3, (2,-1)/sqrt6.
w2, A = normal_modes([[3.0, -2], [-2, 4]], np.diag([1.0, 2]))
assert np.allclose(w2, [1, 4]) and np.allclose(sign_fix(A), np.c_[[1, 1] / SQ3, [2, -1] / np.sqrt(6)])
assert np.isclose(abs(A[:, 0] @ A[:, 1]), 1 / (3 * SQ2)) and np.isclose(1 / (3 * SQ2), 0.2357, atol=5e-5)   # ordinary dot != 0
# A2: double-pendulum shapes M-normalized; A^T K A = g diag(2 -+ sqrt2).
assert np.allclose(A44 * np.sign(A44[0]), [[0.38268, 0.92388], [0.54120, -1.30656]], atol=5e-6)   # first entries positive
assert np.allclose(A44.T @ K44 @ A44, np.diag([5.7466, 33.4934]), atol=5e-4)
assert np.allclose(A44 @ A44.T, np.linalg.inv(M44)) and np.allclose(np.linalg.inv(M44), [[1, -1], [-1, 2]])   # A A^T = M^{-1}
# A3: M = [[1,2],[2,1]] is symmetric but not PD: v = (1,-1) gives -2; eigenvalues 3, -1.
assert np.array([1, -1]) @ np.array([[1, 2], [2, 1]]) @ np.array([1, -1]) == -2 and np.allclose(np.linalg.eigvalsh([[1, 2], [2, 1]]), [-1, 3])
Mb = np.array([[1.0, 1], [0, 1]]); assert np.allclose(np.linalg.eigvalsh((Mb + Mb.T) / 2), [0.5, 1.5]) and Mb[0, 1] != Mb[1, 0]   # PD form, not symmetric
# A4: V = (x^2-1)^2/4 + (y-x)^2: equilibria (0,0) saddle, (+-1,+-1) stable.
x, y = sp.symbols("x y"); V_a4 = (x**2 - 1)**2 / 4 + (y - x)**2
eqs = sp.solve([sp.diff(V_a4, x), sp.diff(V_a4, y)], [x, y], dict=True)
assert sorted((e[x], e[y]) for e in eqs) == [(-1, -1), (0, 0), (1, 1)]
H0 = np.array(sp.hessian(V_a4, (x, y)).subs({x: 0, y: 0}), float); H1 = np.array(sp.hessian(V_a4, (x, y)).subs({x: 1, y: 1}), float)
assert np.allclose(H0, [[1, -2], [-2, 2]]) and np.allclose(H1, [[4, -2], [-2, 2]])
assert np.allclose(np.linalg.eigvalsh(H0), [(3 - np.sqrt(17)) / 2, (3 + np.sqrt(17)) / 2]) and np.allclose(np.linalg.eigvalsh(H0), [-0.5616, 3.5616], atol=5e-5)
assert np.isclose(np.sqrt(0.5616), 0.7494, atol=5e-5) and np.isclose(np.sqrt((3 + np.sqrt(17)) / 2), 1.8872, atol=5e-5)
assert np.allclose(np.linalg.eigvalsh(H1), [3 - np.sqrt(5), 3 + np.sqrt(5)]) and np.allclose(np.linalg.eigvalsh(H1), [0.7639, 5.2361], atol=5e-5)
assert np.allclose(np.sqrt(np.linalg.eigvalsh(H1)), [0.8740, 2.2882], atol=5e-5)
# A5: N = 3 string shapes (1/sqrt2, 1, 1/sqrt2), (1, 0, -1), (1/sqrt2, -1, 1/sqrt2) are eigenvectors with 2 -+ sqrt2, 2.
K, M, d, m = string_KM(3); om0sq = 1 / d**2
for a, lam in [([1 / SQ2, 1, 1 / SQ2], 2 - SQ2), ([1, 0, -1], 2), ([1 / SQ2, -1, 1 / SQ2], 2 + SQ2)]:
    assert np.allclose(K @ a, lam * om0sq * (M @ a))
# A6: l = 0.6, m = 0.2, T_beat = 12 s -> k = m (omega2^2 - omega1^2)/2 = 0.4508; weak-coupling k ~ m omega1 dOmega = 0.4234.
om1 = np.sqrt(9.81 / 0.6); dom = 2 * np.pi / 12; k_exact = 0.2 * ((om1 + dom)**2 - om1**2) / 2; k_weak = 0.2 * om1 * dom
assert np.isclose(k_exact, 0.4508, atol=5e-5) and np.isclose(k_weak, 0.4234, atol=5e-5) and np.isclose(om1, 4.0435, atol=5e-5)
assert np.isclose(dom, 0.5236, atol=5e-5) and np.isclose((k_exact - k_weak) / k_exact, 0.061, atol=5e-4)
# B1: Example 4.2 with a damper on mass 1 only: A^T C A is full. Rayleigh alpha=0.1, beta=0.05: zeta = (0.0827, 0.0710). Fit zeta = 0.02 twice.
ACA = A42.T @ np.diag([1.0, 0]) @ A42
assert np.allclose(ACA, [[(3 - SQ3) / 6, 1 / np.sqrt(6)], [1 / np.sqrt(6), (3 + SQ3) / 6]]) and np.allclose(ACA, [[0.21132, 0.40825], [0.40825, 0.78868]], atol=5e-6)
om42 = np.sqrt(np.diag(A42.T @ K42 @ A42))
zeta_b1 = (0.1 + 0.05 * om42**2) / (2 * om42); assert np.allclose(zeta_b1, [0.08270, 0.07096], atol=5e-6)
beta_fit = 0.04 / (om42[0] + om42[1]); alpha_fit = 0.04 * om42[0] - beta_fit * om42[0]**2
assert np.allclose((alpha_fit + beta_fit * om42**2) / (2 * om42), 0.02) and np.allclose([alpha_fit, beta_fit], [0.020986, 0.017135], atol=5e-7)
# B2: double pendulum m1 = 3, m2 = 1: lambda = 1/(1 +- sqrt(mu)), mu = m2/(m1+m2) = 1/4 -> 2/3, 2; shapes (1, 2), (1, -2).
m1, m2 = 3.0, 1.0; Mb2 = np.array([[m1 + m2, m2], [m2, m2]]); Kb2 = np.diag([m1 + m2, m2])    # l = g = 1
w2, A = normal_modes(Kb2, Mb2); assert np.allclose(w2, [2 / 3, 2]) and np.allclose((A / A[0])[1], [2, -2])
for mu in rng.uniform(0.05, 0.95, 5):
    w2, _ = normal_modes(np.diag([1.0, mu]), np.array([[1.0, mu], [mu, mu]]))
    assert np.allclose(w2, [1 / (1 + np.sqrt(mu)), 1 / (1 - np.sqrt(mu))])
# B3: Kt for Example 4.2, its orthonormal eigenvectors y = L^T a, and the similarity M^{-1}K = L^{-T} Kt L^T.
L42 = np.linalg.cholesky(M42); Kt42 = np.linalg.solve(L42, np.linalg.solve(L42, K42).T).T
assert np.allclose(Kt42, [[2, -1 / SQ2], [-1 / SQ2, 1]]) and np.allclose(np.linalg.eigvalsh(Kt42), [(3 - SQ3) / 2, (3 + SQ3) / 2])
Y42 = L42.T @ A42; assert np.allclose(Y42.T @ Y42, np.eye(2)) and np.allclose(np.abs(Y42), [[0.45970, 0.88807], [0.88807, 0.45970]], atol=5e-6)
assert np.allclose(np.linalg.solve(M42, K42), np.linalg.solve(L42.T, Kt42) @ L42.T)
# B4: modal energies for the beat initial condition of Example 4.1: E_1 = 2.4525 th0^2, E_2 = 2.9525 th0^2, sum 5.405 th0^2; each constant.
q_b, xi_b, xid_b = fit_initial_conditions([1.0, 0.0], [0.0, 0.0], (np.array([9.81, 11.81]), A41, M41))
E_modes = 0.5 * np.array([9.81, 11.81]) * xi_b**2
assert np.allclose(E_modes, [2.4525, 2.9525]) and np.isclose(E_modes.sum(), 0.5 * np.array([1, 0]) @ K41 @ np.array([1, 0]))
tt = np.linspace(0, 30, 7); Qb = q_b(tt); Vb = (q_b(tt + 1e-6) - q_b(tt - 1e-6)) / 2e-6
Xi_t, Xid_t = Qb @ M41 @ A41, Vb @ M41 @ A41
assert np.allclose(0.5 * (Xid_t**2 + np.array([9.81, 11.81]) * Xi_t**2), E_modes, atol=1e-6)
assert np.allclose(q_b(T_transfer)[0, 0], 0, atol=1e-12)                 # pendulum 1 passes through 0 at T_transfer ...
thd1 = (q_b(T_transfer + 1e-6) - q_b(T_transfer - 1e-6))[0, 0] / 2e-6      # ... with |theta1'| <= eps theta0 (theta0 = 1 here)
E_pend1 = 0.5 * thd1**2 * M41[0, 0]; assert abs(thd1) <= eps + 1e-9 and E_pend1 / E_modes.sum() < 0.0025
assert np.isclose(abs(thd1), 0.119, atol=5e-4) and np.isclose(E_pend1 / E_modes.sum(), 0.0013, atol=5e-5)   # values quoted in B4(c)
print(f"B4: E_modes = {E_modes} theta0^2; at T_transfer theta1' = {thd1:.4f} theta0, pendulum-1 share of E = {E_pend1 / E_modes.sum():.5f}")
# B5: N = 5, q0 = (1,2,3,2,1): sine coefficients c_j = (2/(N+1)) sum_k q0_k sin(jk pi/6).
N = 5; kk = np.arange(1, N + 1); S5 = np.sin(np.outer(kk, kk) * np.pi / (N + 1))
assert np.allclose(S5.T @ S5, (N + 1) / 2 * np.eye(N))
c5 = (2 / (N + 1)) * S5.T @ np.array([1, 2, 3, 2, 1.0])
assert np.allclose(c5, [2.48803, 0, -1 / 3, 0, 0.17863], atol=5e-6) and np.allclose(S5 @ c5, [1, 2, 3, 2, 1])
om_b5 = 2 * np.sin(np.arange(1, 6) * np.pi / 12)                                     # omega_j / omega_0 for N = 5
assert np.isclose(om_b5[2] / om_b5[0], 1 + SQ3) and np.isclose(om_b5[4] / om_b5[0], 2 + SQ3)   # B5(d): irrational ratios
# B6: CO2 kicked at atom 1 with v = 1: v_cm = 4/11; modal velocities 16/sqrt44, 16/sqrt32, 16/10.832; x_1 coefficients 4/11, 1/2, 3/22.
q_k, xi_k, xid_k = fit_initial_conditions([0, 0, 0.0], [1.0, 0, 0], (w2_co2, A_co2, M_co2))
assert np.allclose(np.abs(xid_k), [16 / np.sqrt(44), 16 / np.sqrt(32), 16 / norms_co2[2]]) and np.allclose(np.abs(xid_k), [2.41209, 2.82843, 1.47710], atol=5e-6)
coef_x1 = A_co2[0] * xid_k; assert np.allclose(coef_x1, [4 / 11, 0.5, 3 / 22]) and np.allclose(coef_x1, [0.36364, 0.5, 0.13636], atol=5e-6)
assert np.allclose(A_co2[1] * xid_k, [4 / 11, 0, -4 / 11]) and np.allclose(A_co2[2] * xid_k, [4 / 11, -0.5, 3 / 22])
assert np.allclose(np.abs(xi_k), 0) and np.allclose([A_co2[i] @ xid_k for i in range(3)], [1, 0, 0])   # velocities reproduce v0
assert np.allclose([0.5 / om_co2[1], (3 / 22) / om_co2[2], (4 / 11) / om_co2[2]], [2, 0.28485, 0.75961], atol=5e-6)  # sine amplitudes in x_i(t)
# B7 (reference solution): double pendulum (m = l = 1) released from (theta0, 0) at rest.
# Linear: xi = (0.13066, 0.054120) for theta0 = 0.1, E_1 = E_2 = g theta0^2 / 2 exactly. Exact (CM-3) equations:
# [[2, c], [c, 1]] theta'' = (-w2^2 s - 2 g sin th1, w1^2 s - g sin th2), c = cos(th1 - th2), s = sin(th1 - th2).
def double_pendulum_rhs(t, y):
    t1, t2, w1, w2 = y; c, s = np.cos(t1 - t2), np.sin(t1 - t2)
    acc = np.linalg.solve([[2.0, c], [c, 1.0]], [-w2**2 * s - 2 * g * np.sin(t1), w1**2 * s - g * np.sin(t2)])
    return [w1, w2, acc[0], acc[1]]


modes44 = (np.array([2 - SQ2, 2 + SQ2]) * g, A44, M44); t_d = np.linspace(0, 10, 2001); b7 = {}
for th0 in (0.1, 1.0):
    q_d, xi_d, _ = fit_initial_conditions([th0, 0.0], [0, 0.0], modes44)
    E_d = 0.5 * modes44[0] * xi_d**2
    assert np.allclose(E_d, g * th0**2 / 2) and np.isclose(E_d.sum(), 0.5 * np.array([th0, 0]) @ K44 @ np.array([th0, 0]))
    sol = solve_ivp(double_pendulum_rhs, (0, 10), [th0, 0, 0, 0], t_eval=t_d, method="DOP853", rtol=1e-10, atol=1e-12)
    Xi, Xid = sol.y[:2].T @ M44 @ A44, sol.y[2:].T @ M44 @ A44              # project the exact motion onto the modes
    E_t = 0.5 * (Xid**2 + modes44[0] * Xi**2)
    b7[th0] = (np.abs(q_d(t_d) - sol.y[:2].T).max() / th0, np.ptp(E_t[:, 0]) / E_t[0, 0], xi_d)
assert np.allclose(np.abs(b7[0.1][2]), [0.13066, 0.054120], atol=5e-6) and np.isclose(g * 0.1**2 / 2, 0.04905)
assert 0.15 < b7[0.1][0] < 0.2 and b7[0.1][1] < 0.01 and b7[1.0][0] > 1 and b7[1.0][1] > 0.5
assert np.allclose([b7[0.1][0], b7[1.0][0], b7[1.0][1]], [0.171, 1.434, 0.852], atol=5e-4) and np.isclose(b7[0.1][1], 0.0085, atol=5e-5)  # B7(b), (c)
print(f"B7: max|exact - modal|/theta0 over 10 s = {b7[0.1][0]:.3f} (0.1 rad), {b7[1.0][0]:.3f} (1.0 rad); "
      f"relative swing of E_1 = {b7[0.1][1]:.4f}, {b7[1.0][1]:.3f}")
# B8: K^{-1} = A Omega^{-2} A^T; static response to F = e1 is (3/4, 1/2, 1/4); mode 2 has a node at mass 2.
assert np.allclose(np.linalg.inv(K3), A3 @ np.diag(1 / w2_3) @ A3.T) and np.allclose(np.linalg.solve(K3, [1, 0, 0]), [0.75, 0.5, 0.25])
assert np.allclose(np.linalg.inv(K3), np.array([[3, 2, 1], [2, 4, 2], [1, 2, 3]]) / 4)
assert np.allclose(A3[0]**2 / w2_3, [0.42678, 0.25, 0.07322], atol=5e-6)            # modal pieces of the static 3/4
pieces = A3 * (A3[0] / w2_3)                                                           # column j: a_j a_1j / omega_j^2
assert np.allclose(np.abs(pieces.T), [[0.42678, 0.60355, 0.42678], [0.25, 0, 0.25], [0.07322, 0.10355, 0.07322]], atol=5e-6)
assert np.allclose(np.abs(A3[:, 1]), [1 / SQ2, 0, 1 / SQ2]) and np.allclose(zeta3, [0.016892, 0.014142, 0.014651], atol=5e-7)
peak_exact = np.abs(receptance(w2_3, A3, zeta3, np.array([1.0, 0, 0]), om3)[:, 0])  # |q_1| driven at mass 1, at omega = omega_j
peak_single = A3[0]**2 / (2 * zeta3 * w2_3)                                          # one-mode estimate
assert np.allclose(peak_exact, [12.649, 8.853, 2.557], atol=5e-4) and np.allclose(peak_single, [12.632, 8.839, 2.499], atol=5e-4)
print(f"B8: peaks |q_1(omega_j)| = {peak_exact.round(3)}, single-mode estimates {peak_single.round(3)}")
# C1: 5x5 lattice: K = K1 (x) I + I (x) K1, omega^2_{jl} = s_j + s_l; 13 distinct values, multiplicity 5 at 4 omega_0^2, lowest 2(2 - sqrt3).
K1 = 2 * np.eye(5) - np.eye(5, k=1) - np.eye(5, k=-1); K_lat = np.kron(K1, np.eye(5)) + np.kron(np.eye(5), K1)
s = 4 * np.sin(np.arange(1, 6) * np.pi / 12)**2; lam_lat = np.sort(np.add.outer(s, s).ravel())
assert np.allclose(np.linalg.eigvalsh(K_lat), lam_lat) and np.isclose(lam_lat[0], 2 * (2 - SQ3))
distinct = np.unique(np.round(lam_lat, 8)); assert len(distinct) == 13 and np.sum(np.isclose(lam_lat, 4)) == 5
assert np.allclose(lam_lat[:6], [0.5359, 1.2679, 1.2679, 2, 2.2679, 2.2679], atol=5e-5)                 # lowest six (Your turn 1)
# C2: four-DOF car: body heave z, pitch theta (nose up), front/rear axles z_f, z_r.
mb, Ib, mf, mr, a_f, b_r = 1200.0, 1800.0, 90.0, 70.0, 1.1, 1.5
ksf, ksr, ktf, ktr = 45e3, 35e3, 200e3, 200e3
B_car = np.array([[1, a_f, -1, 0], [1, -b_r, 0, -1], [0, 0, 1, 0], [0, 0, 0, 1]], float)   # elongations: front susp, rear susp, front tyre, rear tyre
K_car = spring_K(B_car, [ksf, ksr, ktf, ktr]); M_car = np.diag([mb, Ib, mf, mr])
w2_car, A_car = normal_modes(K_car, M_car); check_modes(K_car, M_car, w2_car, A_car, tol=1e-6)
f_car = np.sqrt(w2_car) / (2 * np.pi); shapes_car = A_car / A_car[np.abs(A_car).argmax(axis=0), np.arange(4)]
k_eff_f, k_eff_r = ksf * ktf / (ksf + ktf), ksr * ktr / (ksr + ktr)                      # suspension and tyre in series
f_bounce, f_pitch = np.sqrt((k_eff_f + k_eff_r) / mb) / (2 * np.pi), np.sqrt((k_eff_f * a_f**2 + k_eff_r * b_r**2) / Ib) / (2 * np.pi)
f_hop_f, f_hop_r = np.sqrt((ksf + ktf) / mf) / (2 * np.pi), np.sqrt((ksr + ktr) / mr) / (2 * np.pi)
x_node = -shapes_car[0, :2] / shapes_car[1, :2]                                           # body point that stays still, m ahead of CG
print(f"C2 car: f = {f_car.round(3)} Hz; one-DOF estimates bounce {f_bounce:.3f}, pitch {f_pitch:.3f}, hop {f_hop_f:.3f}/{f_hop_r:.3f} Hz; "
      f"body-mode nodes at x = {x_node.round(3)} m")
print("C2 car shapes (columns = modes; rows z [m], theta [rad], z_f, z_r [m]; largest entry 1):\n", shapes_car.round(4))
assert np.allclose(K_car[0], [ksf + ksr, ksf * a_f - ksr * b_r, -ksf, -ksr]) and np.isclose(K_car[0, 1], -3000)
assert np.allclose(f_car, [1.170, 1.261, 8.323, 9.237], atol=5e-4) and np.allclose([f_bounce, f_pitch, f_hop_f, f_hop_r], [1.185, 1.252, 8.304, 9.222], atol=5e-4)
assert np.allclose(x_node, [-3.170, 0.474], atol=5e-4) and np.allclose(shapes_car[:2, :2], [[1, -0.4741], [0.3155, 1]], atol=5e-4)
w2_bounce, w2_pitch = (k_eff_f + k_eff_r) / mb, (k_eff_f * a_f**2 + k_eff_r * b_r**2) / Ib
coup = (k_eff_f * a_f - k_eff_r * b_r) / np.sqrt(mb * Ib)                                # mass-normalized bounce-pitch coupling
assert np.allclose([k_eff_f, k_eff_r], [36734.7, 29787.2], atol=0.05) and np.allclose([w2_bounce, w2_pitch, coup], [55.435, 61.928, -2.907], atol=5e-4)
# C3: Caughey series C = M sum_l c_l (M^{-1}K)^l decouples: A^T C A = sum_l c_l Omega^{2l}; fit zeta = (0.02, 0.03, 0.05) on the 3-DOF chain.
Vand = np.vander(w2_3, 3, increasing=True); c_l = np.linalg.solve(Vand, 2 * np.array([0.02, 0.03, 0.05]) * om3)
MinvK3 = np.linalg.solve(M3, K3); C_caughey = M3 @ sum(c * np.linalg.matrix_power(MinvK3, l) for l, c in enumerate(c_l))
assert np.allclose(C_caughey, C_caughey.T) and np.allclose(A3.T @ C_caughey @ A3, np.diag(2 * np.array([0.02, 0.03, 0.05]) * om3))
y_c = 2 * np.array([0.02, 0.03, 0.05]) * om3; d_c = [y_c[1], (y_c[2] - y_c[0]) / (2 * SQ2), (y_c[0] + y_c[2] - 2 * y_c[1]) / 4]
assert np.allclose(d_c, [0.084853, 0.054504, 0.011421], atol=5e-7) and np.allclose([d_c[0] - 2 * d_c[1] + 4 * d_c[2], d_c[1] - 4 * d_c[2], d_c[2]], c_l)
assert np.allclose(c_l, [0.021529, 0.008819, 0.011421], atol=5e-7) and np.allclose(2 * np.array([0.02, 0.03, 0.05]) * om3, [0.030615, 0.084853, 0.184776], atol=5e-7)
print(f"C3 Caughey coefficients (c0, c1, c2) = {c_l.round(6)}")
# Intermediate values quoted in solutions.md.
assert np.allclose(np.sqrt([4 + 2 * SQ2, 4 - 2 * SQ2]), [2.61313, 1.08239], atol=5e-6)                                   # A2
assert np.allclose(M44 @ (A44 * np.sign(A44[0])), [[1.30656, 0.54120], [0.92388, -0.38268]], atol=5e-6) and np.isclose(2 - SQ2, 0.58579, atol=5e-6)
assert np.isclose(om1 + dom / 2, 4.3053, atol=5e-5) and np.isclose(k_exact / k_weak, 1.065, atol=5e-4) and np.isclose(2 * k_exact / 0.2, 4.508, atol=5e-4)  # A6
assert np.allclose([0.1 + 0.05 * 0.63397, 2 * om42[0], 0.1 + 0.05 * 2.36603, 2 * om42[1]], [0.13170, 1.59246, 0.21830, 3.07638], atol=5e-5)  # B1(c)
assert np.isclose(om42.sum(), 2.33442, atol=5e-6) and np.isclose(om42.prod(), np.sqrt(1.5)) and np.isclose(np.sqrt(1.5), 1.22474, atol=5e-6)  # B1(d)
assert np.allclose(Y42[1] / Y42[0], [SQ2 * (2 - w) for w in [(3 - SQ3) / 2, (3 + SQ3) / 2]]) and np.allclose(Y42[1] / Y42[0], [1.93185, -0.51764], atol=5e-6)  # B3(b)
assert np.allclose([0.5 / (2 - SQ2), 0.5 / (2 + SQ2)], [0.85355, 0.14645], atol=5e-6) and np.isclose(0.42678 / 0.75, 0.569, atol=5e-4)  # B8(b)
assert np.isclose(peak_exact[2] / peak_single[2] - 1, 0.023, atol=5e-4) and np.all(np.abs(peak_exact[:2] / peak_single[:2] - 1) < 0.002)  # B8(d)
front_rear = [shapes_car[0, 0] + a_f * shapes_car[1, 0], shapes_car[0, 0] - b_r * shapes_car[1, 0]]
assert np.allclose(front_rear, [1.347, 0.527], atol=5e-4)                                                                   # C2(d)
assert np.all(np.abs(np.array([f_bounce, f_pitch, f_hop_f, f_hop_r]) / f_car - 1) < 0.015) and 7 < f_car[2] / f_car[0] and f_car[3] / f_car[1] < 8
assert np.isclose(w2_pitch - w2_bounce, 6.493, atol=5e-4) and np.isclose(2 * abs(coup) / (w2_pitch - w2_bounce), 0.90, atol=5e-3)
assert np.isclose(np.degrees(np.arctan(2 * abs(coup) / (w2_pitch - w2_bounce)) / 2), 21, atol=0.5)
assert np.allclose(om3, [0.76537, 1.41421, 1.84776], atol=5e-6)                                                             # C3(c)
assert np.allclose([0.38268 * 0.54120, 0.92388 * 1.30656, 0.54120**2, 1.30656**2, 2 * 0.38268**2], [0.20711, 1.20711, 0.29289, 1.70711, 0.29289], atol=5e-5)  # A2 products
assert np.isclose(9.81 / 0.6, 16.35) and np.isclose(E_modes.sum(), 5.405)                                                   # A6, B4(b)
assert 9 < 10 / (2 * np.pi / np.sqrt((2 + SQ2) * g)) < 9.5                                                # B7(e): ~nine fast periods in 10 s
assert np.allclose(s, [2 - SQ3, 1, 2, 3, 2 + SQ3])                                                                           # C1(b)
D_alt = np.diag((-1.0) ** np.arange(1, 6)); assert np.allclose(D_alt @ K1 @ D_alt, 4 * np.eye(5) - K1) and np.allclose(s + s[::-1], 4)  # C1(d)
mult = np.unique(np.unique(np.round(lam_lat, 8), return_counts=True)[1], return_counts=True)
assert dict(zip(mult[0].tolist(), mult[1].tolist())) == {1: 4, 2: 8, 5: 1}                                                   # C1(d) tally
assert np.allclose(K_car, [[80000, -3000, -45000, -35000], [-3000, 133200, -49500, 52500],
                           [-45000, -49500, 245000, 0], [-35000, 52500, 0, 235000]])                                         # C2(a)
print("Part 9: all exercise numbers reproduced.")

# %% [markdown]
# ## Your turn
#
# 1. **5x5 lattice.** Build the 25-DOF transverse lattice (masses on a 5x5 grid, unit
#    springs to the four neighbours, frame fixed): K = K1 (x) I + I (x) K1 with K1 the N = 5
#    string matrix (C1 derives this). Run `normal_modes`, plot the lowest six shapes as
#    5x5 heat maps, and compare the eigenvalues with s_j + s_l. Which of the six are a
#    degenerate pair, and what does `eigh` return for them? Then break the symmetry by
#    making one mass 1% heavier and watch the pair split.
# 2. **Rayleigh damping.** For the triangle (Part 6) with C = 0.05 M + 0.02 K, verify
#    A^T C A is diagonal, then write `damped_free_response(q0, v0, modes, alpha, beta)`
#    using xi_j(t) = e^{-zeta_j w_j t}[...] per mode (solve the damped (3.15) with f = 0,
#    zero modes included: they are damped only by alpha) and check it against `solve_ivp`
#    on M q'' + C q' + K q = 0. Then try C = diag(1, 0, 0, 0, 0, 0) and show it does not decouple.
# 3. **Car model.** Take the C2 parameters, add a road input: the tyre springs attach to
#    a profile y_f(t) = h sin(2 pi v t / lambda) at the front and y_r(t) = y_f(t - (a+b)/v) at
#    the rear. Write F(t), use (3.15) per mode with Rayleigh damping zeta ~ 0.3 in the body
#    modes, and plot body heave and pitch amplitude against speed v for lambda = 10 m.
#    Where are the peaks, and which mode does the wheelbase delay suppress?

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
