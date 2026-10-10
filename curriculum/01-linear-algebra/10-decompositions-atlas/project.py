# %% [markdown]
# # LA-10 Project: The Decompositions Atlas in Code
#
# This script turns the atlas into working code and checks every number in the
# notes. Parts: (1) Cholesky from entry matching (3.3), with a loud failure that
# names the bad pivot; LDL^T and unpivoted LU on SPD matrices. (2) Schur form by
# hand and by SciPy, eigenvectors from T, the real Schur form, and the
# normal-matrix theorem checked numerically. (3) Polar decomposition and the
# orthogonal Procrustes problem, including registering a scanned 3D part to its
# CAD model. (4) Timing LU / Cholesky / QR / SVD solves of one SPD system against
# the flop-count predictions. (5) A `recommend(A, task)` helper that inspects A
# and names the decomposition. (6) The Example 4.4 numbers, Exercise B8, C1.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
import numpy as np
import scipy.linalg as sla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)


def random_spd(n, cond=10.0):
    """A random SPD matrix with eigenvalues spread between 1 and `cond`."""
    Q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    lam = np.exp(rng.uniform(0, np.log(cond), n))
    return Q @ np.diag(lam) @ Q.T


def same_multiset(a, b, tol=1e-8):
    """True if the complex lists a and b agree as multisets (greedy nearest matching)."""
    b = list(b)
    for z in a:
        j = int(np.argmin(np.abs(np.array(b) - z)))
        if abs(b[j] - z) > tol:
            return False
        b.pop(j)
    return not b


# %% [markdown]
# ## Part 1: Cholesky from entry matching
#
# A = L L^T with L lower triangular, positive diagonal. Matching entries for
# j <= i gives (3.3):
#     l_jj = sqrt(a_jj - sum_{k<j} l_jk^2),   l_ij = (a_ij - sum_{k<j} l_ik l_jk) / l_jj.
# The radicand is the j-th pivot d_j of A = L D L^T (Theorem 3.4), so a radicand
# <= 0 proves that A is not SPD, and j is the size of the smallest leading block
# that already fails. Theorem 3.5: for SPD A this never happens.
#
# **Build it yourself first:** write `cholesky(A)` from (3.3) before reading on.


# %%
class NotSPDError(ValueError):
    """Raised when a Cholesky pivot is not positive: A is not SPD."""


def cholesky(A):
    """Cholesky factor L of an SPD matrix A, A = L L^T, eq. (3.3).

    L is lower triangular with positive diagonal. Raises NotSPDError naming the
    first pivot whose radicand is <= 0.
    """
    A = np.asarray(A, float)
    n = A.shape[0]
    if A.shape != (n, n) or not np.allclose(A, A.T):
        raise NotSPDError("A is not square and symmetric")
    L = np.zeros((n, n))
    for j in range(n):
        radicand = A[j, j] - L[j, :j] @ L[j, :j]        # a_jj - sum_{k<j} l_jk^2  (= pivot d_j)
        if radicand <= 0:
            raise NotSPDError(
                f"pivot {j + 1} has radicand {radicand:.4g} <= 0: "
                f"the leading {j + 1}x{j + 1} block of A is not positive definite")
        L[j, j] = np.sqrt(radicand)
        L[j + 1:, j] = (A[j + 1:, j] - L[j + 1:, :j] @ L[j, :j]) / L[j, j]
    return L


def ldlt(A):
    """A = L diag(d) L^T with L unit lower triangular (Theorem 3.3), no pivoting.

    Same recurrence as (3.3) without square roots:
        d_j = a_jj - sum_{k<j} l_jk^2 d_k,    l_ij = (a_ij - sum_{k<j} l_ik d_k l_jk) / d_j.
    Works for indefinite A as long as no pivot is exactly zero. (Reference for Exercise B7.)
    """
    A = np.asarray(A, float)
    n = A.shape[0]
    L = np.eye(n)
    d = np.zeros(n)
    for j in range(n):
        d[j] = A[j, j] - (L[j, :j] ** 2) @ d[:j]
        if d[j] == 0:
            raise ZeroDivisionError(f"pivot {j + 1} is exactly zero; symmetric pivoting needed")
        L[j + 1:, j] = (A[j + 1:, j] - (L[j + 1:, :j] * d[:j]) @ L[j, :j]) / d[j]
    return L, d


def lu_nopivot(A):
    """Plain Gaussian elimination with no row swaps: A = L U. Returns L, U and the pivots diag(U)."""
    A = np.asarray(A, float).copy()
    n = A.shape[0]
    L = np.eye(n)
    for k in range(n - 1):
        if A[k, k] == 0:
            raise ZeroDivisionError(f"zero pivot at step {k + 1}")
        L[k + 1:, k] = A[k + 1:, k] / A[k, k]
        A[k + 1:, k:] -= np.outer(L[k + 1:, k], A[k, k:])
    return L, np.triu(A), np.diag(A).copy()


def cholesky_solve(L, b):
    """Solve (L L^T) x = b: forward substitution L y = b, then back substitution L^T x = y."""
    y = sla.solve_triangular(L, b, lower=True)
    return sla.solve_triangular(L.T, y, lower=False)


# Example 4.1: L = [[2,0,0],[1,2,0],[1,1,2]], pivots (4,4,4), det 64, x = (1,1,1) for b = (8,10,11).
A41 = np.array([[4, 2, 2], [2, 5, 3], [2, 3, 6.]])
L41 = cholesky(A41)
assert np.allclose(L41, [[2, 0, 0], [1, 2, 0], [1, 1, 2]]), L41
assert np.allclose(L41 @ L41.T, A41)
assert np.allclose(L41, sla.cholesky(A41, lower=True))
Lu41, d41 = ldlt(A41)
assert np.allclose(d41, [4, 4, 4]) and np.allclose(Lu41, L41 / 2)       # L_chol = L_unit D^{1/2}
assert np.isclose(np.prod(np.diag(L41)) ** 2, 64) and np.isclose(np.linalg.det(A41), 64)
x41 = cholesky_solve(L41, [8, 10, 11])
assert np.allclose(x41, [1, 1, 1]), x41
print("Example 4.1: L =", L41.tolist(), " pivots =", d41.tolist(), " x =", x41.round(12).tolist())

# Failure diagnostics. [[1,2],[2,1]] fails at pivot 2 (radicand -3); modified 4.1 fails at pivot 3 (radicand -1).
B_mod = np.array([[4, 2, 2], [2, 5, 3], [2, 3, 1.]])
for name, B, bad_pivot in [("[[1,2],[2,1]]", np.array([[1, 2], [2, 1.]]), 2), ("Example 4.1 with a33 = 1", B_mod, 3)]:
    try:
        cholesky(B)
        raise AssertionError("should have failed")
    except NotSPDError as e:
        assert f"pivot {bad_pivot}" in str(e)
        print(f"cholesky({name}) failed as it should: {e}")
    _, dB = ldlt(B)
    n_neg_piv = int(np.sum(dB < 0))
    n_neg_eig = int(np.sum(np.linalg.eigvalsh(B) < 0))
    assert n_neg_piv == n_neg_eig, (dB, np.linalg.eigvalsh(B))              # Sylvester's law of inertia
    print(f"   pivots {np.round(dB, 4).tolist()} -> {n_neg_piv} negative; eigenvalues {np.round(np.linalg.eigvalsh(B), 3).tolist()}")
assert np.allclose(ldlt(B_mod)[1], [4, 4, -1])
assert np.allclose(np.sort(np.linalg.eigvalsh(B_mod)), [-0.7828, 2.4543, 8.3285], atol=1e-4)

# Random SPD matrices: agree with SciPy; unpivoted LU never fails and its pivots are the d_j = l_jj^2.
for _ in range(30):
    n = int(rng.integers(2, 9))
    A = random_spd(n, cond=100)
    L = cholesky(A)
    assert np.allclose(L, sla.cholesky(A, lower=True)) and np.allclose(L @ L.T, A)
    _, U, piv = lu_nopivot(A)
    _, dA = ldlt(A)
    assert np.all(piv > 0) and np.allclose(piv, dA) and np.allclose(piv, np.diag(L) ** 2)
    assert np.linalg.det(A) <= np.prod(np.diag(A)) * (1 + 1e-10)           # Hadamard's inequality (C3)
print("Part 1 checks passed: cholesky == scipy on 30 random SPD matrices; unpivoted LU never failed.")

# %% [markdown]
# ## Part 2: Schur form, and the normal-matrix theorem
#
# Example 4.2: A = [[1,2],[3,0]], eigenvalues 3 and -2. q1 = (1,1)/sqrt2 is the
# eigenvector for 3; q2 = (-1,1)/sqrt2 completes Q; T = Q^T A Q = [[3,-1],[0,-2]].
# Corollary 3.8a: |t12|^2 = ||A||_F^2 - sum |lambda_i|^2 = 14 - 13 = 1 in *every* Schur form.
# Corollary 3.8b: eigenvectors of A are Q y with (T - lambda I) y = 0.
# Theorem 3.10: A normal  <=>  T diagonal.


# %%
A42 = np.array([[1, 2], [3, 0.]])
Q_hand = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
T_hand = Q_hand.T @ A42 @ Q_hand
assert np.allclose(T_hand, [[3, -1], [0, -2]]), T_hand
assert np.allclose(Q_hand.T @ Q_hand, np.eye(2))
T_sp, Z_sp = sla.schur(A42)                       # real Schur form; triangular here since the eigenvalues are real
assert np.allclose(Z_sp @ T_sp @ Z_sp.T, A42)
assert np.isclose(T_sp[1, 0], 0) and np.allclose(sorted(np.diag(T_sp)), [-2, 3])
departure = np.sqrt(np.sum(A42 ** 2) - np.sum(np.diag(T_sp) ** 2))        # Corollary 3.8a
assert np.isclose(abs(T_sp[0, 1]), 1) and np.isclose(departure, 1)
assert not np.allclose(A42 @ A42.T, A42.T @ A42)                            # not normal
print("Example 4.2: hand T =", T_hand.round(12).tolist(), " SciPy T =", T_sp.round(12).tolist(), " |t12| =", round(abs(T_sp[0, 1]), 12))


def eigvec_from_schur(T, Q, k):
    """Unit eigenvector of A = Q T Q^* for the eigenvalue T[k,k], by back substitution on T (Cor. 3.8b).

    Assumes T[k,k] differs from T[i,i] for i < k (distinct eigenvalues).
    """
    n = T.shape[0]
    y = np.zeros(n, dtype=complex)
    y[k] = 1.0
    for i in range(k - 1, -1, -1):                 # row i of (T - lam I) y = 0:  (t_ii - lam) y_i + sum_{j>i} t_ij y_j = 0
        y[i] = -(T[i, i + 1:k + 1] @ y[i + 1:k + 1]) / (T[i, i] - T[k, k])
    x = Q @ y
    return x / np.linalg.norm(x)


# Eigenvector for -2 from the hand Schur form: y = (1/5, 1) -> x ∝ (-2, 3).
x_m2 = eigvec_from_schur(T_hand.astype(complex), Q_hand.astype(complex), 1).real
assert np.allclose(np.abs(x_m2), np.abs([-2, 3]) / np.sqrt(13))
assert np.allclose(A42 @ x_m2, -2 * x_m2)
print("   eigenvector for -2 from T:", x_m2.round(6).tolist(), " (∝ (-2, 3))")

# Random real matrix: the real Schur form is quasi-triangular; the complex Schur form is triangular.
A6 = rng.standard_normal((6, 6))
T_r, Z_r = sla.schur(A6, output="real")
assert np.allclose(Z_r @ T_r @ Z_r.T, A6)
assert np.allclose(np.tril(T_r, -2), 0)                                     # nothing below the first subdiagonal
T_c, Z_c = sla.schur(A6, output="complex")
assert np.allclose(Z_c @ T_c @ Z_c.conj().T, A6) and np.allclose(np.tril(T_c, -1), 0)
assert same_multiset(np.diag(T_c), np.linalg.eigvals(A6))
sub = np.diag(T_r, -1)
for i in np.flatnonzero(np.abs(sub) > 1e-12):                               # each 2x2 block carries a conjugate pair
    assert np.all(np.abs(np.linalg.eigvals(T_r[i:i + 2, i:i + 2]).imag) > 1e-12)
print("   random 6x6: real Schur form has", int(np.sum(np.abs(sub) > 1e-12)), "2x2 blocks; complex Schur form is triangular.")


def random_unitary(n):
    """Random unitary matrix: QR of a complex Gaussian matrix, with the phases of diag(R) normalized away."""
    Z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def strict_upper_norm(T):
    """Frobenius norm of the strictly upper triangular part of T: the departure from normality."""
    return np.linalg.norm(np.triu(T, 1))


# Theorem 3.10 numerically: A = U Λ U^* with complex Λ is normal and its Schur T is diagonal.
for n in (3, 5, 8):
    U = random_unitary(n)
    assert np.allclose(U.conj().T @ U, np.eye(n))
    lam = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    A_norm = U @ np.diag(lam) @ U.conj().T
    assert np.allclose(A_norm @ A_norm.conj().T, A_norm.conj().T @ A_norm)   # normal
    T, Z = sla.schur(A_norm, output="complex")
    assert strict_upper_norm(T) < 1e-10 * np.linalg.norm(A_norm), strict_upper_norm(T)
    assert same_multiset(np.diag(T), lam)
    # Hermitian: real diagonal. Unitary: diagonal on the unit circle. (Corollaries of Theorem 3.10.)
    H = U @ np.diag(lam.real) @ U.conj().T
    assert np.allclose(np.diag(sla.schur(H, output="complex")[0]).imag, 0)
    V = U @ np.diag(np.exp(1j * lam.real)) @ U.conj().T
    assert np.allclose(np.abs(np.diag(sla.schur(V, output="complex")[0])), 1)
    # Non-normal: the strictly upper part has size sqrt(||A||_F^2 - sum|lam|^2), whatever Schur form is returned.
    A_gen = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    T_g, _ = sla.schur(A_gen, output="complex")
    predicted = np.sqrt(np.linalg.norm(A_gen) ** 2 - np.sum(np.abs(np.diag(T_g)) ** 2))
    assert np.isclose(strict_upper_norm(T_g), predicted) and strict_upper_norm(T_g) > 0.1
print("   normal-matrix theorem: Schur T is diagonal to 1e-10 for random U Λ U*; non-normal departure matches Cor. 3.8a.")

# Figure: |T| for a normal and a non-normal 8x8 matrix.
n = 8
U = random_unitary(n)
lam = rng.standard_normal(n) + 1j * rng.standard_normal(n)
T_n, _ = sla.schur(U @ np.diag(lam) @ U.conj().T, output="complex")
T_g, _ = sla.schur(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)), output="complex")
fig, axes = plt.subplots(1, 2, figsize=(9, 4))
for ax, T, title in zip(axes, (T_n, T_g), ("normal A = U Λ U*: T is diagonal", "non-normal A: T is triangular")):
    im = ax.imshow(np.abs(T), cmap="magma")
    ax.set_title(f"{title}\n||strict upper(T)||_F = {strict_upper_norm(T):.2e}")
    ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.046)
fig.suptitle("Schur form T = Q* A Q: entry magnitudes")
fig.tight_layout()
fig.savefig(FIG / "schur_normal_vs_not.png", dpi=120)
plt.close(fig)
print("Part 2 checks passed.")

# %% [markdown]
# ## Part 3: Polar decomposition and orthogonal Procrustes
#
# Theorem 3.13: A = W Σ V^T  =>  A = (W V^T)(V Σ V^T) = U P. For invertible A,
# P^2 = A^T A forces P, then U = A P^{-1}.
# Theorem 3.14: for point sets X (source) and Y (target) as columns, the orthogonal
# R minimizing ||R X - Y||_F is W V^T from the SVD of M = Y X^T, and the minimum
# is ||X||_F^2 + ||Y||_F^2 - 2 sum sigma_i  (3.6).
#
# **Build it yourself first:** write `procrustes(X, Y)`: center both, M = Yc Xc^T,
# SVD, R = W V^T, fix det, t = ybar - R xbar.


# %%
def polar(A):
    """Polar decomposition A = U P (Theorem 3.13): U orthogonal, P = (A^T A)^{1/2} PSD, from the SVD."""
    W, s, Vt = np.linalg.svd(np.asarray(A, float))
    return W @ Vt, Vt.T @ np.diag(s) @ Vt


# A5: [[0,-2],[3,0]] = (rotation by 90°) diag(3,2).  B6: [[3,0],[4,5]] = U P with sqrt5 U = [[2,-1],[1,2]], P = sqrt5 [[2,1],[1,2]].
U5, P5 = polar([[0, -2], [3, 0]])
assert np.allclose(U5, [[0, -1], [1, 0]]) and np.allclose(P5, np.diag([3, 2]))
A_B6 = np.array([[3, 0], [4, 5.]])
U6, P6 = polar(A_B6)
assert np.allclose(U6 * np.sqrt(5), [[2, -1], [1, 2]]) and np.allclose(P6 / np.sqrt(5), [[2, 1], [1, 2]])
assert np.isclose(np.sum((A_B6 - U6) ** 2), 52 - 8 * np.sqrt(5))                     # B6(c)
assert np.isclose(np.degrees(np.arctan2(U6[1, 0], U6[0, 0])), 26.565, atol=1e-3)      # B6(d)
print("A5: U =", U5.round(12).tolist(), " P =", P5.round(12).tolist())
print("B6: sqrt5*U =", (U6 * np.sqrt(5)).round(12).tolist(), " P/sqrt5 =", (P6 / np.sqrt(5)).round(12).tolist(),
      f" ||A-U||_F^2 = {np.sum((A_B6 - U6)**2):.4f} = 52 - 8 sqrt5")

# Uniqueness and the nearest-orthogonal property on random invertible matrices.
for _ in range(20):
    n = int(rng.integers(2, 7))
    A = rng.standard_normal((n, n))
    U, P = polar(A)
    assert np.allclose(U @ P, A) and np.allclose(U.T @ U, np.eye(n))
    assert np.allclose(P, P.T) and np.allclose(P @ P, A.T @ A)        # P is *a* PSD square root of A^T A ...
    cholesky(P)                                                        # ... and positive definite, since A is invertible
    dU = np.linalg.norm(A - U)
    for _ in range(20):                                                # ... and no orthogonal Q is closer to A
        Qr, _ = np.linalg.qr(rng.standard_normal((n, n)))
        assert np.linalg.norm(A - Qr) >= dU - 1e-12
print("   polar: U P = A, P^2 = A^T A, P SPD, and U is the nearest orthogonal matrix (20 random trials).")


def procrustes(X, Y, proper=True):
    """Best rigid alignment of source points X (d x k) onto target points Y (d x k), Theorem 3.14.

    Returns R (orthogonal; a proper rotation if `proper`), the translation t, and
    the residual ||R X + t - Y||_F^2, so that R @ X + t[:, None] ≈ Y.
    """
    X = np.asarray(X, float); Y = np.asarray(Y, float)
    xbar = X.mean(axis=1, keepdims=True); ybar = Y.mean(axis=1, keepdims=True)
    Xc, Yc = X - xbar, Y - ybar
    M = Yc @ Xc.T                                        # cross-covariance, d x d
    W, s, Vt = np.linalg.svd(M)
    D = np.eye(X.shape[0])
    if proper and np.linalg.det(W @ Vt) < 0:             # reflection fix (Kabsch 1976 / Umeyama 1991)
        D[-1, -1] = -1
    R = W @ D @ Vt
    t = (ybar - R @ xbar).ravel()
    resid = np.sum((R @ X + t[:, None] - Y) ** 2)
    return R, t, resid


# Example 4.3: two triangles. M = [[-3,-6],[6,3]], sigma = (9, 3), R = rotation by 90°, residual 16 (40 before).
X43 = np.array([[-1, 2, -1], [-1, -1, 2.]])
Y43 = np.array([[3, -1, -2], [-3, 2, 1.]])
assert np.allclose(X43.sum(1), 0) and np.allclose(Y43.sum(1), 0)
M43 = Y43 @ X43.T
assert np.allclose(M43, [[-3, -6], [6, 3]])
W43, s43, Vt43 = np.linalg.svd(M43)
assert np.allclose(s43, [9, 3])
R43, t43, res43 = procrustes(X43, Y43)
assert np.allclose(R43, [[0, -1], [1, 0]]) and np.allclose(t43, 0)
assert np.isclose(res43, 16) and np.isclose(np.sum(X43 ** 2) + np.sum(Y43 ** 2) - 2 * s43.sum(), 16)   # (3.6)
assert np.isclose(np.sum((X43 - Y43) ** 2), 40)
assert np.allclose(R43 @ (Vt43.T @ np.diag(s43) @ Vt43), M43)                    # M = R (polar stretch)
print(f"Example 4.3: sigma = {s43.round(12).tolist()}, R = {R43.round(12).tolist()}, residual^2 = {res43:.1f} (before: 40), RMS/point = {np.sqrt(res43/3):.3f}")

# Reflection trap: a mirrored target gives det(W V^T) = -1; the fix returns a proper rotation instead.
Y_mirror = np.diag([1, -1]) @ X43
R_refl, _, _ = procrustes(X43, Y_mirror, proper=False)
R_prop, _, res_prop = procrustes(X43, Y_mirror, proper=True)
assert np.isclose(np.linalg.det(R_refl), -1) and np.isclose(np.linalg.det(R_prop), 1)
print(f"   mirrored target: unconstrained answer is a reflection (det -1); best proper rotation leaves residual^2 = {res_prop:.2f}")


def rotation_axis_angle(axis, angle):
    """Rodrigues' formula: the rotation by `angle` (radians) about the unit vector along `axis`."""
    k = np.asarray(axis, float); k = k / np.linalg.norm(k)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


def rotation_angle_between(R1, R2):
    """Angle (degrees) of the rotation R1^T R2: how far apart two rotations are."""
    c = (np.trace(R1.T @ R2) - 1) / 2
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def cad_bracket(n_plate=200, n_wall=150, n_hole=50):
    """Synthetic 'CAD model' of an L-bracket with two holes, as a 3 x N point cloud."""
    plate = np.c_[rng.uniform(0, 4, n_plate), rng.uniform(0, 2, n_plate), np.zeros(n_plate)]
    wall = np.c_[np.zeros(n_wall), rng.uniform(0, 2, n_wall), rng.uniform(0, 3, n_wall)]
    th = rng.uniform(0, 2 * np.pi, n_hole)
    holes = np.vstack([np.c_[cx + 0.3 * np.cos(th), 1 + 0.3 * np.sin(th), np.zeros(n_hole)] for cx in (1.0, 3.0)])
    return np.vstack([plate, wall, holes]).T


# ME application: register a noisy scan of a part to its CAD model.
X_cad = cad_bracket()
R_true = rotation_axis_angle([1, 2, 0.5], np.radians(40))
t_true = np.array([1.5, -0.7, 2.0])
sigma_noise = 0.02
Y_scan = R_true @ X_cad + t_true[:, None] + sigma_noise * rng.standard_normal(X_cad.shape)
R_est, t_est, res_after = procrustes(X_cad, Y_scan)
rms_before = np.sqrt(np.mean(np.sum((X_cad - Y_scan) ** 2, axis=0)))
rms_after = np.sqrt(res_after / X_cad.shape[1])
assert rotation_angle_between(R_est, R_true) < 0.5, rotation_angle_between(R_est, R_true)
assert np.linalg.norm(t_est - t_true) < 0.05
assert rms_before > 1.0 and rms_after < 2.5 * sigma_noise * np.sqrt(3)
print(f"   registration: rotation error {rotation_angle_between(R_est, R_true):.3f} deg, "
      f"translation error {np.linalg.norm(t_est - t_true):.4f}, RMS before {rms_before:.3f}, after {rms_after:.4f} "
      f"(noise level sigma*sqrt3 = {sigma_noise*np.sqrt(3):.4f})")

fig = plt.figure(figsize=(11, 4.5))
X_aligned = R_est @ X_cad + t_est[:, None]
for k, (pts, title) in enumerate([(X_cad, f"before: RMS = {rms_before:.2f}"),
                                  (X_aligned, f"after Procrustes: RMS = {rms_after:.3f}")], 1):
    ax = fig.add_subplot(1, 2, k, projection="3d")
    ax.scatter(*Y_scan, s=3, c="C3", alpha=0.5, label="scan")
    ax.scatter(*pts, s=3, c="C0", alpha=0.5, label="CAD (transformed)" if k == 2 else "CAD")
    ax.set_title(title); ax.legend(fontsize=8)
    both = np.c_[Y_scan, pts]
    ax.set_xlim(both[0].min(), both[0].max()); ax.set_ylim(both[1].min(), both[1].max()); ax.set_zlim(both[2].min(), both[2].max())
fig.suptitle("Registering a scanned part to its CAD model: the polar factor of the cross-covariance")
fig.tight_layout()
fig.savefig(FIG / "procrustes_registration.png", dpi=120)
plt.close(fig)
print("Part 3 checks passed.")

# %% [markdown]
# ## Part 4: Timing one SPD solve by LU, Cholesky, QR, SVD
#
# Flop counts (leading terms, n x n): LU 2n^3/3; Cholesky n^3/3; Householder QR
# 4n^3/3 plus 4n^3/3 to form Q explicitly (8n^3/3 total); full SVD about 20n^3.
# The exponent 3 should be visible. The ratios will not match exactly: BLAS is
# multithreaded and blocked, and each factorization blocks differently (notes §8).
# The dashed curves are the flop predictions scaled to the flop rate LU achieves.


# %%
def solve_lu(A, b):
    lu, piv = sla.lu_factor(A)
    return sla.lu_solve((lu, piv), b)


def solve_chol(A, b):
    c, low = sla.cho_factor(A)
    return sla.cho_solve((c, low), b)


def solve_qr(A, b):
    Q, R = sla.qr(A, mode="economic")
    return sla.solve_triangular(R, Q.T @ b)


def solve_svd(A, b):
    W, s, Vt = np.linalg.svd(A)
    return Vt.T @ ((W.T @ b) / s)


methods = {"LU": (solve_lu, 2 / 3), "Cholesky": (solve_chol, 1 / 3), "QR": (solve_qr, 8 / 3), "SVD": (solve_svd, 20.0)}
ns = [200, 400, 800, 1200, 1600, 2000]
times = {name: [] for name in methods}
for n in ns:
    B = rng.standard_normal((n, n))
    A = B.T @ B + n * np.eye(n)                       # SPD, moderately conditioned
    b = rng.standard_normal(n)
    x_ref = np.linalg.solve(A, b)
    for name, (f, _) in methods.items():
        t0 = time.perf_counter()
        x = f(A, b)
        times[name].append(time.perf_counter() - t0)
        assert np.allclose(x, x_ref, rtol=1e-6, atol=1e-8), name
ns_arr = np.array(ns, float)
# Fit t = c n^3 on the three largest n; express every prediction at the flop rate measured for LU.
c_meas = {name: float(np.mean(np.array(times[name])[-3:] / ns_arr[-3:] ** 3)) for name in methods}
rate_LU = (2 / 3) / c_meas["LU"]                      # flops per second achieved by LU
print(f"{'method':9s} {'flops/n^3':>10s} {'meas t/n^3 [s]':>15s} {'ratio to LU':>12s} {'predicted':>10s}")
for name, (_, coef) in methods.items():
    print(f"{name:9s} {coef:10.3f} {c_meas[name]:15.3e} {c_meas[name]/c_meas['LU']:12.2f} {coef/(2/3):10.2f}")

fig, ax = plt.subplots(figsize=(7, 4.8))
for k, (name, (_, coef)) in enumerate(methods.items()):
    ax.loglog(ns, times[name], "o-", color=f"C{k}", label=f"{name} measured")
    ax.loglog(ns, coef * ns_arr ** 3 / rate_LU, "--", color=f"C{k}", alpha=0.6,
              label=f"{name} predicted: {coef:.2g} n³ at LU's flop rate")
ax.set_xlabel("n"); ax.set_ylabel("seconds"); ax.set_title("Solving one SPD system: time vs n")
ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(FIG / "timing_solve.png", dpi=120)
plt.close(fig)
print("Part 4: timings recorded; see figures/timing_solve.png.")

# %% [markdown]
# ## Part 5: `recommend(A, task)`: the decision table as code
#
# Inspect shape, symmetry, definiteness (by attempting Cholesky), rank and the
# condition number; then apply §3.9.


# %%
def inspect_matrix(A, tol=1e-10):
    """The facts about A that the decision table needs."""
    A = np.asarray(A, float)
    m, n = A.shape
    s = np.linalg.svd(A, compute_uv=False)
    rank = int(np.sum(s > s[0] * max(m, n) * np.finfo(float).eps * 10)) if s[0] > 0 else 0
    facts = dict(shape=(m, n), square=(m == n), rank=rank, full_rank=(rank == min(m, n)),
                 cond=(s[0] / s[-1] if rank == min(m, n) else np.inf),
                 symmetric=False, spd=False, normal=False, spd_failure=None)
    if m == n:
        scale = max(np.abs(A).max(), 1e-300)
        facts["symmetric"] = bool(np.allclose(A, A.T, atol=tol * scale))
        facts["normal"] = bool(np.allclose(A @ A.T, A.T @ A, atol=tol * scale ** 2))
        if facts["symmetric"]:
            try:
                cholesky(A)
                facts["spd"] = True
            except NotSPDError as e:
                facts["spd_failure"] = str(e)
    return facts


def recommend(A, task):
    """Which decomposition to use for `task` on A, and why (decision table, notes §3.9).

    task in {'solve', 'solve_many', 'lstsq', 'eig', 'rank', 'lowrank', 'nearest_orthogonal', 'matrix_function'}.
    Returns (decomposition, explanation) and prints them.
    """
    f = inspect_matrix(A)
    m, n = f["shape"]
    shape = f"{m}x{n}"
    if task in ("solve", "solve_many"):
        if not f["square"]:
            rec, why = "QR (least squares)", f"A is {shape}, not square: there is no 'solve', only least squares."
        elif not f["full_rank"]:
            rec, why = "SVD (pseudoinverse)", f"A is square but has rank {f['rank']} < {n}: the SVD shows the null space and gives the minimum-norm solution."
        elif f["spd"]:
            rec, why = "Cholesky", f"A is SPD (all {n} pivots positive): n^3/3 flops, no pivoting, and the factorization doubles as the definiteness proof."
        elif f["symmetric"]:
            rec, why = "LDL^T with symmetric pivoting", f"A is symmetric but not SPD ({f['spd_failure']}): keep the n^3/3 cost and get the inertia for free."
        else:
            rec, why = "PLU", "A is square, nonsingular and unstructured: 2n^3/3 flops with partial pivoting."
        if task == "solve_many" and f["square"] and f["full_rank"]:
            why += " Factor once, then 2n^2 flops per right-hand side."
        if np.isfinite(f["cond"]) and f["cond"] > 1e8:
            why += f" Warning: condition number {f['cond']:.1e}, expect to lose about {np.log10(f['cond']):.0f} digits (consider the SVD)."
    elif task == "lstsq":
        if f["full_rank"] and f["cond"] < 1e8:
            rec, why = "QR (Householder)", f"A is {shape} with full column rank and cond {f['cond']:.1e}: QR avoids squaring the condition number (Theorem 3.7)."
        else:
            rec, why = "SVD (truncated / pseudoinverse)", f"A is {shape} with rank {f['rank']} and cond {f['cond']:.1e}: only the SVD shows which directions are undetermined."
    elif task == "eig":
        if not f["square"]:
            rec, why = "none", "Eigenvalues need a square matrix; you probably want singular values (SVD)."
        elif f["symmetric"]:
            rec, why = "spectral (symmetric QR algorithm / eigh)", "A is symmetric: orthonormal eigenbasis, real eigenvalues, perfectly conditioned basis change."
        elif f["normal"]:
            rec, why = "Schur (will come out diagonal)", "A is normal but not symmetric: unitarily diagonalizable with complex eigenvalues (Theorem 3.10)."
        else:
            rec, why = "Schur, then eigenvectors from T", "A is not normal: Schur always exists and is computed stably; derive eigenvectors from T if you need them (Cor. 3.8b)."
    elif task in ("rank", "lowrank"):
        rec, why = "SVD", f"Rank, null space, conditioning and low-rank approximation all read off the singular values (here rank {f['rank']}, cond {f['cond']:.1e})."
    elif task == "nearest_orthogonal":
        rec, why = "polar via SVD", "The nearest orthogonal matrix is the orthogonal polar factor W V^T (Theorem 3.14)."
    elif task == "matrix_function":
        if f["symmetric"]:
            rec, why = "spectral", "f(A) = Q f(Λ) Q^T: apply f to the eigenvalues."
        else:
            rec, why = "Schur (Schur–Parlett)", "f(A) = Q f(T) Q^*, with f(T) by a triangular recurrence; never via the Jordan form."
    else:
        raise ValueError(f"unknown task {task!r}")
    print(f"recommend({shape}, {task}): {rec}\n    because {why}")
    return rec, why


spd = random_spd(6, cond=50)
gen = rng.standard_normal((6, 6))
Qo, _ = np.linalg.qr(rng.standard_normal((6, 6)))
sym_indef = Qo @ np.diag([3, 1, -2, 0.5, -1, 4.]) @ Qo.T
tall = rng.standard_normal((50, 4))
tall_def = tall @ np.array([[1, 0, 0, 2], [0, 1, 0, -1], [0, 0, 1, 0], [0, 0, 0, 0.]])    # rank 3
assert recommend(spd, "solve")[0] == "Cholesky"
assert recommend(spd, "solve_many")[0] == "Cholesky"
assert recommend(gen, "solve")[0] == "PLU"
assert recommend(sym_indef, "solve")[0].startswith("LDL^T")
assert recommend(tall, "lstsq")[0].startswith("QR")
assert recommend(tall_def, "lstsq")[0].startswith("SVD")
assert recommend(sym_indef, "eig")[0].startswith("spectral")
assert recommend(gen, "eig")[0].startswith("Schur")
assert recommend(tall_def, "rank")[0] == "SVD"
assert recommend(gen, "nearest_orthogonal")[0].startswith("polar")
assert recommend(spd, "matrix_function")[0] == "spectral"
assert recommend(gen, "matrix_function")[0].startswith("Schur")
print("Part 5 checks passed.")

# %% [markdown]
# ## Part 6: The Example 4.4 numbers, the chain of springs (B8), transient growth (C1)


# %%
# (a) Stiffness: the fixed-free 2-spring chain is SPD with pivots 2, 1/2; the free-free chain fails at pivot 2 with radicand 0.
K_ff = np.array([[2, -1], [-1, 1.]])
assert np.allclose(ldlt(K_ff)[1], [2, 0.5]) and np.allclose(np.diag(cholesky(K_ff)) ** 2, [2, 0.5])
try:
    cholesky([[1, -1], [-1, 1.]])
    raise AssertionError
except NotSPDError as e:
    assert "pivot 2" in str(e) and "radicand 0" in str(e)
    print("free-free chain:", e)

# (b) Modal analysis: M = diag(2,1), K = [[3,-1],[-1,1]] -> omega^2 = 1/2, 2; modes (1,2), (1,-1); M-orthogonal.
M = np.diag([2.0, 1.0])
K = np.array([[3, -1], [-1, 1.]])
L_M = cholesky(M)
L_inv = np.linalg.inv(L_M)
C = L_inv @ K @ L_inv.T                                                      # L^{-1} K L^{-T}: symmetric
assert np.allclose(C, [[1.5, -1 / np.sqrt(2)], [-1 / np.sqrt(2), 1]]) and np.allclose(C, C.T)
w2, Y = np.linalg.eigh(C)
assert np.allclose(w2, [0.5, 2])
Phi = L_inv.T @ Y                                                            # phi = L^{-T} y
assert np.allclose(Phi.T @ M @ Phi, np.eye(2))                               # M-orthonormal modes
assert np.allclose(Phi[:, 0] / Phi[0, 0], [1, 2]) and np.allclose(Phi[:, 1] / Phi[0, 1], [1, -1])
w2_sp, Phi_sp = sla.eigh(K, M)
assert np.allclose(w2_sp, [0.5, 2]) and np.allclose(np.abs(Phi_sp), np.abs(Phi))
MinvK = np.linalg.inv(M) @ K
assert not np.allclose(MinvK, MinvK.T)                                        # M^{-1} K is not symmetric
print("modal analysis: omega^2 =", w2.round(12).tolist(), " modes (phi_1 = 1):", (Phi / Phi[0]).T.round(6).tolist())

# (c) Calibration: x in [10, 10.1], columns 1, x, x^2.  kappa(A^T A) = kappa(A)^2.
x = np.linspace(10, 10.1, 200)
A_cal = np.c_[np.ones_like(x), x, x ** 2]
kA, kAtA = np.linalg.cond(A_cal), np.linalg.cond(A_cal.T @ A_cal)
assert np.isclose(kAtA, kA ** 2, rtol=1e-3)
assert 1e7 < kA < 2e7 and 1e14 < kAtA < 2e14
print(f"calibration matrix: kappa(A) = {kA:.3e}, kappa(A^T A) = {kAtA:.3e} = kappa(A)^2 -> normal equations lose ~{np.log10(kAtA):.0f} digits")
c_true = np.array([1.0, -0.2, 0.01])
y = A_cal @ c_true + 1e-6 * rng.standard_normal(x.size)
Qc, Rc = sla.qr(A_cal, mode="economic")
c_qr = sla.solve_triangular(Rc, Qc.T @ y)
c_ne = np.linalg.solve(A_cal.T @ A_cal, A_cal.T @ y)
print(f"   coefficient error: QR {np.linalg.norm(c_qr - c_true):.2e}, normal equations {np.linalg.norm(c_ne - c_true):.2e}")


# B8: the fixed-free chain K_n has pivots (k+1)/k and last pivot 1/n, det 1; the free-free chain fails at the last step.
def chain_stiffness(n, fixed_left=True):
    """Stiffness matrix of n unit springs in series, right end free, left end fixed or free."""
    K = 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    K[-1, -1] = 1
    if not fixed_left:
        K[0, 0] = 1
    return K


for n in (3, 5, 8):
    _, piv = ldlt(chain_stiffness(n))
    expected = np.r_[(np.arange(1, n) + 1) / np.arange(1, n), 1 / n]
    assert np.allclose(piv, expected) and np.isclose(np.prod(piv), 1)
    try:
        cholesky(chain_stiffness(n, fixed_left=False))
        raise AssertionError
    except NotSPDError as e:
        assert f"pivot {n}" in str(e)
assert np.allclose(ldlt(chain_stiffness(3))[1], [2, 1.5, 1 / 3])
print("B8: pivots of K_3 =", ldlt(chain_stiffness(3))[1].round(12).tolist(), " det = 1; the free-free chain fails at the last pivot.")

# C1(c): transient growth of a non-normal matrix with spectral radius 0.9.
A_tr = np.array([[0.9, 1000.0], [0, 0.9]])
norms = [np.linalg.norm(np.linalg.matrix_power(A_tr, k), 2) for k in range(1, 120)]
k_peak = int(np.argmax(norms)) + 1
assert 3800 < max(norms) < 3900 and k_peak in (9, 10) and norms[-1] < 1
print(f"C1: ||A^k||_2 for rho(A) = 0.9 peaks at k = {k_peak} with value {max(norms):.1f} before decaying (||A^119||_2 = {norms[-1]:.3f}).")
print("Part 6 checks passed.")

# %% [markdown]
# ## Your turn
#
# 1. **`sqrtm_spd(A)`.** Implement the SPD square root via the spectral
#    decomposition: A = Q Λ Q^T  =>  A^{1/2} = Q Λ^{1/2} Q^T (notes §3.9, "matrix
#    functions"). Check `sqrtm_spd(A) @ sqrtm_spd(A) ≈ A`, compare with
#    `scipy.linalg.sqrtm`, and use it to recompute the polar factor P of a random A
#    as `sqrtm_spd(A.T @ A)` (Theorem 3.13 step 3 says you must get the same P as
#    `polar`). Then feed it the modified Example 4.1 matrix `B_mod`: what goes
#    wrong, and what does `scipy.linalg.sqrtm` return instead?
# 2. **Unpivoted LU cannot fail on an SPD matrix.** Try to construct an SPD A for
#    which `lu_nopivot(A)` meets a zero pivot. You will not succeed. Write the proof
#    as a comment below: the k-th pivot is the (1,1) entry of the Schur complement
#    left after k-1 steps; Theorem 3.5 step 4 shows that Schur complement is SPD, so
#    its (1,1) entry e_1^T S e_1 is positive. (Equivalently d_k = det A_k / det A_{k-1}
#    with both determinants positive, because leading principal submatrices of an SPD
#    matrix are SPD.) Which of the two arguments generalizes to `ldlt` on a matrix
#    that is merely symmetric and invertible, and why does it fail there?
# 3. **Procrustes with scale (Exercise C2).** Extend `procrustes` to return a scale s:
#    s* = (sum of singular values of M, with the sign fix) / ||X_c||_F^2, R = W D V^T
#    unchanged. Generate `Y = 1.7 * R_true @ X_cad + t[:, None] + noise` and recover
#    s to three digits. Then explain why plain `procrustes` already returns the right
#    rotation when the scale is wrong (hint: what does scaling Y by c > 0 do to M, W, V?).
# 4. In Part 4, replace the SPD matrix by a general one and drop Cholesky; then add a
#    "50 right-hand sides" variant that factors once and solves 50 times. Check the
#    ratio against Exercise A3.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
