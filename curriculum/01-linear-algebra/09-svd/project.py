# %% [markdown]
# # LA-9 Project: The SVD, Low-Rank Approximation, PCA, Conditioning, Pseudoinverse, Polar
#
# This script implements the §3.2 derivation as code and then checks every
# numerical claim in the notes. Parts: (1) `svd_via_eig` and why it is not how
# the SVD is computed, (2) every number in Examples 4.1–4.5, (3) image
# compression with the Eckart–Young error curve, (4) PCA on a 3-D Gaussian
# cloud two ways, (5) condition number and error amplification, (6) the
# pseudoinverse as the minimum-norm least-squares solution, (7) polar
# decomposition, (8) the spectral norm by random sampling.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.linalg import sqrtm

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
EPS = np.finfo(float).eps


def random_orthogonal(n, rng=rng):
    """A Haar-ish random orthogonal n x n matrix: Q from the QR of a Gaussian matrix."""
    Q, R = np.linalg.qr(rng.standard_normal((n, n)))
    return Q * np.sign(np.diag(R))  # fix column signs so the distribution is not biased


def projector(Q):
    """Orthogonal projector Q Q^T onto the span of the (orthonormal) columns of Q."""
    return Q @ Q.T


# %% [markdown]
# ## Part 1: the SVD via `eigh(A^T A)` (the §3.2 derivation as code)
#
# Steps 1–6 of Theorem 3.2: diagonalize A^T A = V Λ V^T (LA-8), set σ_i = √λ_i,
# u_i = A v_i / σ_i for σ_i > 0 (eq. 3.4), complete U by Gram–Schmidt, return
# A = U Σ V^T in economy shape (U is m x n, s has n entries, Vt is n x n) so it
# can be compared with `np.linalg.svd(A, full_matrices=False)`.
#
# **Build it yourself first.** Write `svd_via_eig` before reading this one.


# %%
def complete_orthonormal(Q, m, rng=rng):
    """Extend the orthonormal columns of Q (m x k) to an orthonormal basis of R^m (m x m).

    Gram–Schmidt (via QR) of [Q | random] and keep the new directions; step 5 of §3.2.
    """
    k = Q.shape[1]
    if k == m:
        return Q
    M = np.hstack([Q, rng.standard_normal((m, m - k))])
    Qfull, _ = np.linalg.qr(M)
    # QR may flip signs of the first k columns; keep the originals exactly.
    return np.hstack([Q, Qfull[:, k:]])


def svd_via_eig(A, tol=None):
    """Economy SVD of A (m >= n) computed from the eigendecomposition of A^T A.

    Returns U (m x n), s (n,), Vt (n x n) with A = U diag(s) Vt and s decreasing.
    Columns of U for s_i <= tol are completed to be orthonormal to the rest.
    This is the EXISTENCE PROOF as an algorithm; Part 1b shows why it is a bad
    algorithm (small singular values are lost). Exercise B3 extends it to m < n.
    """
    A = np.asarray(A, float)
    m, n = A.shape
    if m < n:
        raise NotImplementedError("wide case: use A A^T instead (Exercise B3)")
    lam, V = np.linalg.eigh(A.T @ A)          # (3.3): A^T A = V Λ V^T, Λ ascending
    order = np.argsort(lam)[::-1]              # sort to σ_1 >= σ_2 >= ...
    lam, V = lam[order], V[:, order]
    s = np.sqrt(np.clip(lam, 0.0, None))       # σ_i = sqrt(λ_i); clip rounding negatives
    if tol is None:
        # The A^T A route cannot tell σ_i < sqrt(ε) σ_1 from zero (λ_i drowns in ε σ_1²),
        # so this is the honest rank cutoff for THIS method; LAPACK can use ε σ_1.
        tol = max(m, n) * np.sqrt(EPS) * (s[0] if s.size else 0.0)
    r = int(np.sum(s > tol))
    U_r = (A @ V[:, :r]) / s[:r]               # (3.4): u_i = A v_i / σ_i
    U = complete_orthonormal(U_r, m)[:, :n]    # step 5, trimmed to economy shape
    return U, s, V.T


def four_subspaces(A, tol=1e-10):
    """Orthonormal bases of C(A), N(A^T), C(A^T), N(A) read off the SVD (Theorem 3.3)."""
    U, s, Vt = np.linalg.svd(A)
    r = int(np.sum(s > tol * s[0]))
    return U[:, :r], U[:, r:], Vt[:r].T, Vt[r:].T


# Theorem 3.2 / 3.3 checks on random tall matrices (full rank and rank deficient).
for trial in range(30):
    m = int(rng.integers(2, 8)); n = int(rng.integers(1, m + 1))
    r = int(rng.integers(1, n + 1))
    A = rng.standard_normal((m, r)) @ rng.standard_normal((r, n))   # rank r
    U, s, Vt = svd_via_eig(A)
    U2, s2, Vt2 = np.linalg.svd(A, full_matrices=False)
    assert np.allclose(U @ np.diag(s) @ Vt, A), "reconstruction A = U Σ V^T"
    # zero singular values come back as sqrt(rounding noise) ~ sqrt(EPS)*σ_1 ~ 1e-8 σ_1 (Part 1b)
    assert np.allclose(s, s2, atol=1e-6 * s2[0]), "singular values agree"
    assert np.allclose(U.T @ U, np.eye(n)) and np.allclose(Vt @ Vt.T, np.eye(n)), "orthonormal"
    # Signs and null-space bases differ; subspaces (as projectors) must not.
    assert np.allclose(projector(U[:, :r]), projector(U2[:, :r])), "C(A) projector"
    assert np.allclose(projector(Vt[:r].T), projector(Vt2[:r].T)), "C(A^T) projector"
    # (3.1'): A v_i = σ_i u_i, and A kills v_{r+1..n}.
    assert np.allclose(A @ Vt.T, U * s, atol=1e-6 * s[0])
    # Theorem 3.3: the four subspaces, checked as "A maps / kills" statements.
    CA, NAt, CAt, NA = four_subspaces(A)
    assert CA.shape[1] == r and NA.shape[1] == n - r and NAt.shape[1] == m - r
    assert np.allclose(A @ NA, 0) and np.allclose(A.T @ NAt, 0)
    assert np.allclose(projector(CA) @ A, A), "A's columns lie in span(u_1..u_r)"
    assert np.allclose(A @ projector(CAt), A), "A's rows lie in span(v_1..v_r)"
    # §3.3: u_i are eigenvectors of A A^T with eigenvalue σ_i^2.
    assert np.allclose((A @ A.T) @ U[:, :r], U[:, :r] * s[:r] ** 2)
    # (3.7), (3.8): norms from singular values.
    assert np.isclose(np.linalg.norm(A, 2), s[0]) and np.isclose(np.linalg.norm(A, "fro") ** 2, np.sum(s ** 2))
print("Part 1a: svd_via_eig agrees with np.linalg.svd on 30 random matrices (up to signs).")

# %% [markdown]
# ### Part 1b: where `svd_via_eig` fails
# λ_i of A^T A is computed with absolute error ≈ ε‖A^T A‖ = ε σ_1². For σ_3 = 1e-9,
# λ_3 = 1e-18 is far below that noise (ε σ_1² ≈ 2e-16), so σ_3 comes back as
# garbage. `np.linalg.svd` (Golub–Kahan, LA-13) gets absolute error ≈ ε σ_1 in
# σ_i itself, so σ_3 is correct to about 7 digits.

# %%
true_s = np.array([1.0, 1e-5, 1e-9])
A_bad = random_orthogonal(3) @ np.diag(true_s) @ random_orthogonal(3).T
_, s_eig, _ = svd_via_eig(A_bad)
s_lapack = np.linalg.svd(A_bad, compute_uv=False)
rel_eig = np.abs(s_eig - true_s) / true_s
rel_lapack = np.abs(s_lapack - true_s) / true_s
print("Part 1b: true σ        =", true_s)
print("         via eigh(A^TA) =", s_eig, " rel. err =", rel_eig)
print("         np.linalg.svd  =", s_lapack, " rel. err =", rel_lapack)
assert rel_lapack[2] < 1e-5, "LAPACK gets the tiny singular value to ~7 digits"
assert rel_eig[2] > 1e-2, "the A^T A route loses the tiny singular value"
assert rel_eig[0] < 1e-12 and rel_lapack[0] < 1e-12, "both get σ_1"

# %% [markdown]
# ## Part 2: every number in the worked examples
#
# Example 4.1 (3,0;4,5), 4.2 (rank-one 3x2), 4.3 (PCA on four points),
# 4.4 (condition number), 4.5 (polar decomposition of simple shear).


# %%
def rank_k(U, s, Vt, k):
    """A_k = sum_{i<=k} σ_i u_i v_i^T, eq. (3.6) truncated."""
    return (U[:, :k] * s[:k]) @ Vt[:k]


def same_up_to_sign(a, b):
    return np.allclose(a, b) or np.allclose(a, -b)


# --- Example 4.1 ---
A41 = np.array([[3.0, 0.0], [4.0, 5.0]])
U, s, Vt = np.linalg.svd(A41)
assert np.allclose(A41.T @ A41, [[25, 20], [20, 25]])
assert np.allclose(s, [np.sqrt(45), np.sqrt(5)])
assert same_up_to_sign(Vt[0], [1, 1] / np.sqrt(2)) and same_up_to_sign(Vt[1], [1, -1] / np.sqrt(2))
assert same_up_to_sign(U[:, 0], [1, 3] / np.sqrt(10)) and same_up_to_sign(U[:, 1], [3, -1] / np.sqrt(10))
# u_i = A v_i / σ_i with the hand-chosen signs, and the pair flips together:
v1, v2 = np.array([1, 1]) / np.sqrt(2), np.array([1, -1]) / np.sqrt(2)
assert np.allclose(A41 @ v1 / np.sqrt(45), [1, 3] / np.sqrt(10)) and np.allclose(A41 @ v2 / np.sqrt(5), [3, -1] / np.sqrt(10))
assert np.isclose(s[0] * s[1], abs(np.linalg.det(A41))) and np.isclose(np.sum(s ** 2), 50)
A1 = rank_k(U, s, Vt, 1)
assert np.allclose(A1, [[1.5, 1.5], [4.5, 4.5]]) and np.isclose(np.linalg.norm(A41 - A1, "fro") ** 2, 5)
assert np.isclose(np.linalg.cond(A41), 3.0)
print("Example 4.1: σ =", s, "; NumPy's signs:", np.sign(U[0, 0]), np.sign(Vt[0, 0]), "(both flipped together)")

# Figure: the unit circle and its image ellipse, with the singular vectors.
th = np.linspace(0, 2 * np.pi, 400)
circle = np.c_[np.cos(th), np.sin(th)].T
ellipse = A41 @ circle
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4.8))
ax0.plot(circle[0], circle[1], "k", lw=1)
for v, c, lab in [(v1, "C0", "v1"), (v2, "C1", "v2")]:
    ax0.quiver(0, 0, *v, angles="xy", scale_units="xy", scale=1, color=c, label=lab)
ax0.set_title("input: unit circle, right singular vectors"); ax0.set_aspect("equal"); ax0.set_xlim(-1.5, 1.5); ax0.set_ylim(-1.5, 1.5)
ax0.grid(alpha=.3); ax0.legend()
ax1.plot(ellipse[0], ellipse[1], "k", lw=1)
for v, sig, c, lab in [(v1, s[0], "C0", "σ1 u1 = A v1"), (v2, s[1], "C1", "σ2 u2 = A v2")]:
    ax1.quiver(0, 0, *(A41 @ v), angles="xy", scale_units="xy", scale=1, color=c, label=lab)
ax1.set_title("output: the ellipse A(circle), semi-axes √45 and √5"); ax1.set_aspect("equal")
ax1.set_xlim(-8, 8); ax1.set_ylim(-8, 8); ax1.grid(alpha=.3); ax1.legend()
fig.tight_layout(); fig.savefig(FIG / "svd_ellipse.png", dpi=120); plt.close(fig)

# --- Example 4.2 ---
A42 = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])
U, s, Vt = np.linalg.svd(A42)
assert np.allclose(A42.T @ A42, [[14, 28], [28, 56]])
assert np.isclose(s[0], np.sqrt(70)) and np.isclose(s[1], 0)
assert same_up_to_sign(U[:, 0], [1, 2, 3] / np.sqrt(14)) and same_up_to_sign(Vt[0], [1, 2] / np.sqrt(5))
assert same_up_to_sign(Vt[1], [2, -1] / np.sqrt(5)) and np.allclose(A42 @ Vt[1], 0)
u2, u3 = np.array([2, -1, 0]) / np.sqrt(5), np.array([3, 6, -5]) / np.sqrt(70)
Uhand = np.c_[[1, 2, 3] / np.sqrt(14), u2, u3]
assert np.allclose(Uhand.T @ Uhand, np.eye(3)), "hand-completed U is orthogonal"
assert np.allclose(projector(Uhand[:, 1:]), projector(U[:, 1:])), "N(A^T) = plane x+2y+3z=0 either way"
assert np.allclose(np.linalg.pinv(A42), A42.T / 70)
print("Example 4.2: σ1 =", s[0], "= √70; N(A^T) spanned by (2,-1,0)/√5 and (3,6,-5)/√70; A+ = A^T/70.")

# --- Example 4.3 ---
raw = np.array([[5.0, 4.0], [4.0, 1.0], [1.0, 0.0], [2.0, 3.0]])
mu = raw.mean(axis=0); X43 = raw - mu
assert np.allclose(mu, [3, 2]) and np.allclose(X43, [[2, 2], [1, -1], [-2, -2], [-1, 1]])
S43 = X43.T @ X43 / 3
assert np.allclose(X43.T @ X43, [[10, 6], [6, 10]])
lam, Vs = np.linalg.eigh(S43)
assert np.allclose(lam[::-1], [16 / 3, 4 / 3])
U, s, Vt = np.linalg.svd(X43, full_matrices=False)
assert np.allclose(s, [4, 2]) and np.allclose(s ** 2 / 3, lam[::-1])          # Theorem 3.12(a)
assert same_up_to_sign(Vt[0], [1, 1] / np.sqrt(2)) and same_up_to_sign(Vt[1], [1, -1] / np.sqrt(2))
scores1 = X43 @ (np.array([1, 1]) / np.sqrt(2)); scores2 = X43 @ (np.array([1, -1]) / np.sqrt(2))
assert np.allclose(scores1, [2 * np.sqrt(2), 0, -2 * np.sqrt(2), 0]) and np.allclose(scores2, [0, np.sqrt(2), 0, -np.sqrt(2)])
assert np.isclose(np.sum(scores1 ** 2), 16) and np.isclose(np.sum(scores2 ** 2), 4)
assert np.isclose(16 / 20, s[0] ** 2 / np.sum(s ** 2))                          # 80% explained
X1 = rank_k(U, s, Vt, 1)
assert np.allclose(X1, [[2, 2], [0, 0], [-2, -2], [0, 0]]) and np.isclose(np.linalg.norm(X43 - X1) ** 2, 4)
lam_raw, V_raw = np.linalg.eigh(raw.T @ raw)
assert np.allclose(raw.T @ raw, [[46, 30], [30, 26]]) and np.allclose(lam_raw[::-1], [36 + np.sqrt(1000), 36 - np.sqrt(1000)])
ang_raw = np.degrees(np.arctan2(*V_raw[::-1, -1] * np.sign(V_raw[0, -1])))
assert abs(ang_raw - 35.8) < 0.1, ang_raw
print(f"Example 4.3: σ = {s}, eigenvalues of S = {lam[::-1]}, PC1 at 45°, uncentered 'PC1' at {ang_raw:.1f}°.")

# --- Example 4.4 ---
A44 = np.array([[1.0, 1.0], [1.0, 1.0001]]); b44 = np.array([2.0, 2.0001])
s44 = np.linalg.svd(A44, compute_uv=False)
kappa44 = s44[0] / s44[1]
assert np.allclose(s44, [2.00005000125, 0.00004999875], rtol=1e-9) and abs(kappa44 - 40002) < 0.01
assert np.isclose(kappa44, np.linalg.cond(A44))
x44 = np.linalg.solve(A44, b44); assert np.allclose(x44, [1, 1])
db = np.array([1e-4, -1e-4]); dx = np.linalg.solve(A44, b44 + db) - x44
assert np.allclose(dx, [2.0001, -2.0], atol=1e-8)
assert np.allclose(np.linalg.inv(A44), [[10001, -10000], [-10000, 10000]], rtol=1e-8)
rel_x = np.linalg.norm(dx) / np.linalg.norm(x44); rel_b = np.linalg.norm(db) / np.linalg.norm(b44)
assert rel_x <= kappa44 * rel_b * (1 + 1e-9), "bound (3.10)"
assert abs(rel_x / rel_b - kappa44) < 1e-3 * kappa44, "bound (3.10) essentially attained"
print(f"Example 4.4: κ2 = {kappa44:.5f}; rel. change in b = {rel_b:.3e}, in x = {rel_x:.5f}, ratio = {rel_x / rel_b:.2f}")

# --- Example 4.5 ---
F = np.array([[1.0, 1.0], [0.0, 1.0]]); phi = (1 + np.sqrt(5)) / 2
U, s, Vt = np.linalg.svd(F)
assert np.allclose(s, [phi, 1 / phi])
Q = U @ Vt; S = Vt.T @ np.diag(s) @ Vt                      # (3.13): Q = U V^T, S = V Σ V^T
assert np.allclose(Q, np.array([[2, 1], [-1, 2]]) / np.sqrt(5)) and np.allclose(Q @ S, F)
theta = np.degrees(np.arctan2(Q[1, 0], Q[0, 0]))
assert np.isclose(theta, -np.degrees(np.arctan(0.5))) and abs(theta + 26.565) < 1e-2
assert np.allclose(S * np.sqrt(5), [[2, 1], [1, 3]]) and np.allclose(S @ S, F.T @ F) and np.allclose(S, sqrtm(F.T @ F).real)
assert same_up_to_sign(Vt[0], np.array([1, phi]) / np.sqrt(1 + phi ** 2)) and same_up_to_sign(Vt[1], np.array([phi, -1]) / np.sqrt(1 + phi ** 2))
assert same_up_to_sign(U[:, 0], np.array([phi, 1]) / np.sqrt(phi + 2)) and same_up_to_sign(U[:, 1], np.array([1, -phi]) / np.sqrt(phi + 2))
assert np.isclose(np.degrees(np.arctan(phi)), 58.28, atol=0.01) and np.isclose(np.degrees(np.arctan(1 / phi)), 31.72, atol=0.01)
V_left = Q @ S @ Q.T
assert np.allclose(V_left * np.sqrt(5), [[3, 1], [1, 2]]) and np.allclose(V_left @ V_left, F @ F.T)
E = 0.5 * (F.T @ F - np.eye(2))
assert np.allclose(np.sort(np.linalg.eigvalsh(E)), [-1 / (2 * phi), phi / 2])
assert np.isclose(np.linalg.cond(F), phi ** 2) and np.isclose(np.linalg.cond(F, 1), 4) and np.isclose(np.linalg.cond(F, np.inf), 4)
print(f"Example 4.5: stretches {s[0]:.4f}, {s[1]:.4f} (φ, 1/φ); rotation {theta:.2f}°; κ2 = φ² = {phi**2:.4f}, κ1 = κ∞ = 4.")

# %% [markdown]
# ## Part 3: image compression and Eckart–Young
#
# A 256 x 256 synthetic image (no downloads): a linear gradient (rank 2), a disc,
# a rectangle (rank 1), a diagonal band (Toeplitz, not low rank), and a
# sinusoidal texture. We compress at ranks 1, 5, 20, 50 and check Theorem 3.7:
# ‖A − A_k‖_F = sqrt(sum_{i>k} σ_i²) exactly, for every k. Then we add noise and
# watch the singular values develop a floor (§3.13).


# %%
def synthetic_image(N=256):
    """Overlapping shapes and gradients on an N x N grid, values in [0, 1]."""
    y, x = np.mgrid[0:N, 0:N] / (N - 1)
    img = 0.3 * x + 0.2 * y                                               # gradient (rank 2)
    img += 0.5 * (((x - 0.35) ** 2 + (y - 0.4) ** 2) < 0.06)               # disc
    img += 0.4 * ((x > 0.55) & (x < 0.9) & (y > 0.55) & (y < 0.8))        # rectangle (rank 1)
    img += 0.3 * (np.abs(x - y) < 0.04)                                   # diagonal band
    img += 0.15 * np.sin(14 * np.pi * x) * np.cos(10 * np.pi * y)          # texture (rank 1)
    return np.clip(img / img.max(), 0, 1)


img = synthetic_image()
U, s, Vt = np.linalg.svd(img, full_matrices=False)
ranks = [1, 5, 20, 50]
fig, axes = plt.subplots(1, len(ranks) + 1, figsize=(16, 3.6))
axes[0].imshow(img, cmap="gray", vmin=0, vmax=1); axes[0].set_title("original (rank %d)" % np.linalg.matrix_rank(img))
for ax, k in zip(axes[1:], ranks):
    Ak = rank_k(U, s, Vt, k)
    err = np.linalg.norm(img - Ak, "fro")
    pred = np.sqrt(np.sum(s[k:] ** 2))                                   # Eckart–Young prediction
    assert np.isclose(err, pred, rtol=1e-10), (k, err, pred)
    ax.imshow(Ak, cmap="gray", vmin=0, vmax=1)
    ax.set_title(f"rank {k}: {k*(2*256+1)} numbers\n‖A−A_k‖_F/‖A‖_F = {err/np.linalg.norm(img):.3f}")
for ax in axes:
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(FIG / "image_compression.png", dpi=110); plt.close(fig)

# Error-vs-rank for every k, overlaid with the Eckart–Young prediction.
ks = np.arange(0, 256)
errs = np.array([np.linalg.norm(img - rank_k(U, s, Vt, k), "fro") for k in ks])
preds = np.sqrt(np.cumsum(s[::-1] ** 2)[::-1])                            # sqrt(sum_{i>k} σ_i²), k = 0..255
preds = np.append(preds, 0.0)[ks]
assert np.allclose(errs, preds, atol=1e-9), "Eckart–Young equality for all k"
# Any other rank-k matrix does worse: random rank-k guess, and a rank-k "crop" of A.
k = 20
B_rand = rng.standard_normal((256, k)) @ rng.standard_normal((k, 256)) * s[0] / 256
B_crop = img.copy(); B_crop[k:, :] = 0                                     # rank <= k (only k nonzero rows)
assert np.linalg.norm(img - B_rand) > errs[k] and np.linalg.norm(img - B_crop) > errs[k]
assert np.isclose(np.linalg.norm(img - rank_k(U, s, Vt, k), 2), s[k])      # Theorem 3.8 at k = 20

# Noise: a low-rank-plus-noise spectrum has a flat floor.
noisy = img + 0.05 * rng.standard_normal(img.shape)
s_noisy = np.linalg.svd(noisy, compute_uv=False)
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4))
ax0.semilogy(ks[1:], errs[1:], "k", lw=2, label="‖A − A_k‖_F (computed)")
ax0.semilogy(ks[1:], preds[1:], "C3--", lw=1.5, label="Eckart–Young  √(Σ_{i>k} σ_i²)")
ax0.set_xlabel("rank k"); ax0.set_ylabel("Frobenius error"); ax0.set_title("Theorem 3.7: the curves coincide"); ax0.legend(); ax0.grid(alpha=.3)
ax1.semilogy(np.arange(1, 257), s, "k", label="σ_i of the image")
ax1.semilogy(np.arange(1, 257), s_noisy, "C0", alpha=.8, label="σ_i of image + N(0, 0.05²) noise")
ax1.axhline(0.05 * np.sqrt(256) * 2, color="C3", ls=":", label="≈ 2·τ√N: top of the noise floor")
ax1.set_xlabel("i"); ax1.set_ylabel("σ_i"); ax1.set_title("scree plot: signal decays, noise is flat"); ax1.legend(); ax1.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "eckart_young_curve.png", dpi=120); plt.close(fig)
s_pure = np.linalg.svd(0.05 * rng.standard_normal(img.shape), compute_uv=False)
assert s_noisy[100] > 2 * s[100], "the noise floor lifts the tail"
assert s[0] / s[100] > 50 * s_pure[0] / s_pure[100], "image spectrum decays ~280x by i=100; pure noise only ~2x (flat floor)"
print(f"Part 3: Eckart–Young equality verified for k = 0..255; rank-20 relative error {errs[20]/np.linalg.norm(img):.4f}.")

# %% [markdown]
# ## Part 4: PCA on a correlated 3-D Gaussian cloud, two ways
#
# Theorem 3.12: eigenvectors of S = X^T X/(n−1) are the right singular vectors of
# the centered X, with eigenvalues σ_i²/(n−1). We generate data with a KNOWN
# covariance R diag(9, 1, 0.09) R^T (standard deviations 3, 1, 0.3 along a random
# rotation R), then recover it.
#
# **Build it yourself first:** write `pca(X_raw)` returning (components V, explained
# variances σ_i²/(n−1), scores X V) before reading this one.


# %%
def pca(X_raw):
    """PCA by the SVD of the centered data (Theorem 3.12).

    Returns components (columns = v_i), variances σ_i²/(n−1), scores U Σ = X V, and the mean.
    """
    X_raw = np.asarray(X_raw, float)
    mu = X_raw.mean(axis=0)
    X = X_raw - mu                                      # centering is mandatory (§8)
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    n = X.shape[0]
    return Vt.T, s ** 2 / (n - 1), U * s, mu


def pca_via_covariance(X_raw):
    """PCA by eigh of the sample covariance (3.12): the same thing, formed the slow way."""
    X = X_raw - X_raw.mean(axis=0)
    S = X.T @ X / (X.shape[0] - 1)
    lam, V = np.linalg.eigh(S)
    order = np.argsort(lam)[::-1]
    return V[:, order], lam[order]


n_pts = 2000
R_true = random_orthogonal(3)
sd_true = np.array([3.0, 1.0, 0.3])
cloud = np.array([1.0, -2.0, 0.5]) + (rng.standard_normal((n_pts, 3)) * sd_true) @ R_true.T
V_pca, var_pca, scores, mu_c = pca(cloud)
V_cov, var_cov = pca_via_covariance(cloud)
assert np.allclose(var_pca, var_cov), "variances agree: eigh(S) vs SVD"
for i in range(3):
    assert abs(abs(V_pca[:, i] @ V_cov[:, i]) - 1) < 1e-10, "components agree up to sign"
    assert abs(abs(V_pca[:, i] @ R_true[:, i]) - 1) < 0.02, "and recover the true axes"
assert np.allclose(np.sqrt(var_pca), sd_true, rtol=0.08), "recovered standard deviations ≈ 3, 1, 0.3"
assert np.allclose(scores.mean(axis=0), 0) and np.allclose(scores.T @ scores / (n_pts - 1), np.diag(var_pca)), "scores are uncorrelated"
assert np.isclose(var_pca.sum(), np.var(cloud, axis=0, ddof=1).sum()), "total variance = trace S"

fig = plt.figure(figsize=(15, 4.5))
ax = fig.add_subplot(1, 3, 1, projection="3d")
ax.scatter(*cloud.T, s=2, alpha=.25)
for i, c in enumerate("C3 C2 C1".split()):
    d = 2 * np.sqrt(var_pca[i]) * V_pca[:, i]          # 2 standard deviations along each axis
    ax.quiver(*mu_c, *d, color=c, lw=2.5, arrow_length_ratio=.15)
ax.set_title("cloud with principal axes (length 2σ_i/√(n−1))"); ax.set_box_aspect((1, 1, 1))
lim = 8; ax.set_xlim(mu_c[0] - lim, mu_c[0] + lim); ax.set_ylim(mu_c[1] - lim, mu_c[1] + lim); ax.set_zlim(mu_c[2] - lim, mu_c[2] + lim)
ax = fig.add_subplot(1, 3, 2)
ax.bar([1, 2, 3], var_pca, color=["C3", "C2", "C1"]); ax.plot([1, 2, 3], sd_true ** 2, "k_", ms=25, mew=2, label="true variances 9, 1, 0.09")
ax.set_xticks([1, 2, 3]); ax.set_xlabel("component"); ax.set_ylabel("variance σ_i²/(n−1)"); ax.set_title(f"scree plot: {100*var_pca[:2].sum()/var_pca.sum():.1f}% in two components"); ax.legend()
ax = fig.add_subplot(1, 3, 3)
ax.scatter(scores[:, 0], scores[:, 1], s=2, alpha=.3); ax.set_aspect("equal")
ax.set_xlabel("score 1 (σ_1 u_1)"); ax.set_ylabel("score 2 (σ_2 u_2)"); ax.set_title("projection onto the first two components"); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "pca_cloud.png", dpi=120); plt.close(fig)
print("Part 4: PCA via eigh(S) == PCA via SVD; recovered sd =", np.sqrt(var_pca).round(3), "(true 3, 1, 0.3).")

# %% [markdown]
# ## Part 5: condition number and error amplification
#
# Theorem 3.9: ‖δx‖/‖x‖ ≤ κ_2 ‖δb‖/‖b‖, attained for b ∝ u_1, δb ∝ u_n. We build
# 6 x 6 matrices with κ_2 from 1 to 1e12 and measure (i) the worst-case
# amplification, (ii) the amplification of a random δb, (iii) the actual rounding
# error of `np.linalg.solve`, which behaves like κ_2 · ε (the "lose log10 κ
# digits" rule).

# %%
kappas = np.logspace(0, 12, 13)
amp_worst, amp_rand, fwd_err = [], [], []
n = 6
for kappa in kappas:
    Uk, Vk = random_orthogonal(n), random_orthogonal(n)
    sig = np.logspace(0, -np.log10(kappa), n)                 # σ_1 = 1, σ_n = 1/κ
    A = Uk @ np.diag(sig) @ Vk.T
    assert np.isclose(np.linalg.cond(A), kappa, rtol=1e-6)
    # (i) worst case: b = u_1 (x = v_1/σ_1), δb = η u_n (δx = η v_n/σ_n)
    eta = 1e-6
    b = Uk[:, 0]; x = np.linalg.solve(A, b); dx = np.linalg.solve(A, b + eta * Uk[:, -1]) - x
    amp_worst.append((np.linalg.norm(dx) / np.linalg.norm(x)) / (eta / np.linalg.norm(b)))
    # (ii) random δb on a random b
    b = rng.standard_normal(n); db = 1e-6 * rng.standard_normal(n)
    x = np.linalg.solve(A, b); dx = np.linalg.solve(A, b + db) - x
    amp_rand.append((np.linalg.norm(dx) / np.linalg.norm(x)) / (np.linalg.norm(db) / np.linalg.norm(b)))
    # (iii) forward error of solve with a known x: rounding of (A, b) is the perturbation
    x_true = rng.standard_normal(n); b = A @ x_true
    fwd_err.append(np.linalg.norm(np.linalg.solve(A, b) - x_true) / np.linalg.norm(x_true))
amp_worst, amp_rand, fwd_err = map(np.array, (amp_worst, amp_rand, fwd_err))
assert np.allclose(amp_worst[:9], kappas[:9], rtol=1e-4), "worst case attains κ (while rounding is negligible)"
assert np.all(amp_rand <= kappas * (1 + 1e-6)), "random perturbations never exceed the bound (3.10)"
assert np.all(fwd_err <= 1e3 * kappas * EPS + 1e-15), "solve's error is a modest multiple of κ·ε"
fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.loglog(kappas, amp_worst, "C3o-", label="worst case: b ∝ u_1, δb ∝ u_n")
ax.loglog(kappas, amp_rand, "C0s-", label="random b, δb")
ax.loglog(kappas, kappas, "k--", lw=1, label="κ_2 (the bound 3.10)")
ax.loglog(kappas, fwd_err / EPS, "C2^-", label="rounding error of solve ÷ ε")
ax.set_xlabel("κ_2(A)"); ax.set_ylabel("relative error amplification"); ax.set_title("digits lost ≈ log10 κ")
ax.legend(fontsize=8); ax.grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig(FIG / "condition_number.png", dpi=120); plt.close(fig)
print("Part 5: worst-case amplification equals κ; random amplification ≤ κ; solve error ≈ κ·ε. Digits kept at κ=1e12:",
      f"{-np.log10(fwd_err[-1]):.1f}")

# %% [markdown]
# ## Part 6: the pseudoinverse on a rank-deficient least-squares problem
#
# Theorem 3.10: x+ = A+ b minimizes ‖Ax − b‖, lies in C(A^T), and has the smallest
# norm among all minimizers x+ + z, z ∈ N(A). Theorem 3.11(a): A+ = (A^T A)^{-1} A^T
# when the columns are independent.


# %%
def pinv_svd(A, tol=1e-12):
    """Moore–Penrose pseudoinverse A+ = V Σ+ U^T, eq. (3.11)."""
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    r = int(np.sum(s > tol * s[0]))
    return (Vt[:r].T / s[:r]) @ U[:, :r].T


A = rng.standard_normal((8, 3)) @ rng.standard_normal((3, 5))          # 8 x 5, rank 3
b = rng.standard_normal(8)
Ap = pinv_svd(A)
assert np.allclose(Ap, np.linalg.pinv(A))
x_plus = Ap @ b
x_lstsq = np.linalg.lstsq(A, b, rcond=None)[0]
assert np.allclose(x_plus, x_lstsq), "lstsq returns the minimum-norm solution too"
CA, NAt, CAt, NA = four_subspaces(A)
assert NA.shape[1] == 2
res_min = np.linalg.norm(A @ x_plus - b)
assert np.isclose(res_min, np.linalg.norm(projector(NAt) @ b)), "residual = part of b outside C(A)"
assert np.allclose(NA.T @ x_plus, 0), "x+ ⟂ N(A), i.e. x+ ∈ C(A^T)"
for _ in range(200):
    z = NA @ rng.standard_normal(2)
    assert np.isclose(np.linalg.norm(A @ (x_plus + z) - b), res_min), "same residual for every x+ + z"
    assert np.linalg.norm(x_plus + z) >= np.linalg.norm(x_plus) - 1e-12, "x+ has minimum norm"
    assert np.isclose(np.linalg.norm(x_plus + z) ** 2, np.linalg.norm(x_plus) ** 2 + np.linalg.norm(z) ** 2)  # Pythagoras
# A random x that is NOT a minimizer has a larger residual.
assert np.linalg.norm(A @ rng.standard_normal(5) - b) > res_min
# Penrose equations and the two projectors (Theorem 3.11(d), Exercise B5).
assert np.allclose(A @ Ap @ A, A) and np.allclose(Ap @ A @ Ap, Ap)
assert np.allclose(A @ Ap, projector(CA)) and np.allclose(Ap @ A, projector(CAt))
# Full column rank: A+ = (A^T A)^{-1} A^T.
A_full = rng.standard_normal((8, 5))
assert np.allclose(pinv_svd(A_full), np.linalg.inv(A_full.T @ A_full) @ A_full.T)
assert np.isclose(np.linalg.cond(A_full.T @ A_full), np.linalg.cond(A_full) ** 2), "κ(A^T A) = κ(A)^2 (B7)"
print(f"Part 6: pseudoinverse checks passed; min residual {res_min:.4f}, ‖x+‖ = {np.linalg.norm(x_plus):.4f}.")

# %% [markdown]
# ## Part 7: polar decomposition A = Q S of random matrices
#
# Theorem 3.13: Q = U V^T orthogonal, S = V Σ V^T symmetric PSD, S = sqrt(A^T A).
# The left version A = S' Q has S' = U Σ U^T = sqrt(A A^T).


# %%
def polar(A):
    """Polar decomposition A = Q S (3.13): Q = U V^T orthogonal, S = V Σ V^T PSD."""
    U, s, Vt = np.linalg.svd(A)
    return U @ Vt, Vt.T @ np.diag(s) @ Vt


for n in (2, 3, 4, 6):
    for _ in range(5):
        A = rng.standard_normal((n, n))
        Q, S = polar(A)
        assert np.allclose(Q.T @ Q, np.eye(n)) and np.allclose(Q @ S, A)
        assert np.allclose(S, S.T) and np.all(np.linalg.eigvalsh(S) >= -1e-12)
        assert np.allclose(S, sqrtm(A.T @ A).real), "S = sqrt(A^T A)"
        assert np.isclose(np.linalg.det(Q), np.sign(np.linalg.det(A))), "det Q = sign det A"
        S_left = Q @ S @ Q.T
        assert np.allclose(S_left @ Q, A) and np.allclose(S_left, sqrtm(A @ A.T).real)
        # Exercise C1: Q is the nearest orthogonal matrix to A (check against random orthogonals).
        dQ = np.linalg.norm(A - Q)
        assert all(np.linalg.norm(A - random_orthogonal(n)) >= dQ - 1e-12 for _ in range(50))
print("Part 7: polar decomposition verified on random 2x2 .. 6x6 matrices; Q is nearest orthogonal.")

# %% [markdown]
# ## Part 8: ‖A‖_2 = σ_1 by sampling random unit vectors (Theorem 3.5)

# %%
A = rng.standard_normal((3, 3))
sig = np.linalg.svd(A, compute_uv=False)
Xs = rng.standard_normal((20_000, 3)); Xs /= np.linalg.norm(Xs, axis=1, keepdims=True)
stretch = np.linalg.norm(Xs @ A.T, axis=1)
assert stretch.max() <= sig[0] * (1 + 1e-12) and stretch.min() >= sig[-1] * (1 - 1e-12), "σ_n ≤ ‖Ax‖ ≤ σ_1 on the sphere"
assert stretch.max() > 0.98 * sig[0], "random sampling gets within 2% of σ_1 in 3-D"
U, s, Vt = np.linalg.svd(A)
assert np.isclose(np.linalg.norm(A @ Vt[0]), s[0]) and np.isclose(np.linalg.norm(A @ Vt[-1]), s[-1]), "v_1 and v_n attain the extremes"
print(f"Part 8: σ_1 = {sig[0]:.5f}, max over 20000 random unit x of ‖Ax‖ = {stretch.max():.5f}; σ_3 = {sig[-1]:.5f}, min = {stretch.min():.5f}.")

# %% [markdown]
# ## Your turn
#
# 1. **Your own image.** Replace `synthetic_image` with a pattern of your own
#    (text rendered as a boolean mask, a checkerboard, concentric rings, a
#    Gaussian blob plus stripes). Before running, predict: which of your
#    ingredients are exactly low rank, and roughly how fast will the singular
#    values decay? Compress at ranks 1, 5, 20, 50 and compare with your prediction.
# 2. **Recover the rank.** Generate `X = L @ R.T + tau * noise` with L (n x k),
#    R (p x k) Gaussian, n = 500, p = 40, k = 4, and tau = 0.1, 1, 3. Run `pca`
#    and plot the scree. For which tau can you still read off k = 4? Compare the
#    noise floor to tau * (sqrt(n) + sqrt(p)), the rough size of the largest
#    singular value of a pure-noise matrix.
# 3. **Procrustes (Exercise C1).** Implement `procrustes(A, B)` returning the
#    orthogonal Q minimizing ‖A − BQ‖_F via the SVD of B^T A. Test it on
#    `A = B @ Q0 + 0.01 * noise` for a random rotation Q0, and in 2-D compare
#    against a brute-force search over rotation angles.
# 4. **(Exercise B3)** Extend `svd_via_eig` to m < n via A A^T, and plot the
#    relative error of each σ_i against σ_i for a matrix with singular values
#    1, 1e-2, ..., 1e-12.
# 5. **(Exercise B8)** Add ‖A − A_k‖_2 = σ_{k+1} to the error plot, then find the
#    rank k that minimizes ‖(noisy image)_k − (clean image)‖_F.

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
