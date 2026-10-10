# %% [markdown]
# # LA-4 Project: Linear Maps, Their Matrices, and Change of Basis
#
# This script builds the central objects of LA-4 and checks every number in
# the notes. Parts: (1) coordinate vectors and `matrix_of_map`, applied to the
# differentiation operator D on P_3 in two bases and checked against SymPy;
# (2) the two theorems that make matrices work, [Tv]_C = [T]_{B,C}[v]_B and
# composition = multiplication; (3) change of basis, the 90-degree rotation in
# the basis {(1,1),(1,-1)}, and the similarity invariants trace and determinant
# on random P; (4) the projection onto y = 2x with P^2 = P; (5) grid
# deformations, the image of the unit circle (an ellipse: SVD preview), and
# areas as |det|; (6) rank–nullity on random maps; (7) the 2D stress tensor in
# a rotated frame and Mohr's circle.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from scipy.linalg import null_space
from pathlib import Path

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# One fixed categorical palette for every figure (blue = input/before,
# orange = output/after, aqua = a third object), grays for grids.
BLUE, ORANGE, AQUA, INK, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#c8c7c2"

# %% [markdown]
# ## Part 1: Coordinates and the matrix of a map
#
# Definition 3.3: [v]_B is the unique coefficient vector with v = sum_j c_j b_j.
# Definition 3.5 / eq. (3.2): column j of [T]_{B,C} is [T(b_j)]_C.
#
# Everything here is written against *ambient* coordinates: a polynomial in
# P_3 is stored as its coefficient vector (c_0, c_1, c_2, c_3) in the monomial
# basis, a 2x2 matrix is stored flattened, a vector in R^n is itself. A "basis"
# is a list of such ambient vectors. `coords` solves for coordinates in that
# basis; `matrix_of_map` applies the definition column by column.
#
# **Build it yourself first:** write `coords(v, basis)` and
# `matrix_of_map(T, basis_in, basis_out)` before reading these.


# %%
def coords(v, basis):
    """[v]_B (Definition 3.3): solve sum_j c_j b_j = v for the coefficient vector c.

    `basis` is a list of ambient vectors; v must lie in their span (asserted).
    """
    Bmat = np.column_stack([np.asarray(b, float).ravel() for b in basis])
    v = np.asarray(v, float).ravel()
    c, *_ = np.linalg.lstsq(Bmat, v, rcond=None)
    assert np.allclose(Bmat @ c, v, atol=1e-10), "v is not in the span of the given basis"
    return c


def matrix_of_map(T, basis_in, basis_out):
    """[T]_{B,C}, eq. (3.2): column j is the C-coordinate vector of T(b_j).

    T is any Python function on ambient vectors; linearity is assumed, not checked.
    """
    return np.column_stack([coords(T(np.asarray(b, float)), basis_out) for b in basis_in])


# Polynomials in P_3 as coefficient vectors (c_0, c_1, c_2, c_3): p(x) = sum c_k x^k.
def D(p):
    """Differentiation on P_3: (c_0, c_1, c_2, c_3) -> (c_1, 2c_2, 3c_3, 0)."""
    p = np.asarray(p, float)
    return np.array([k * p[k] for k in range(1, len(p))] + [0.0])


monomial = [np.eye(4)[k] for k in range(4)]                        # 1, x, x^2, x^3
factorial = [np.eye(4)[k] / math.factorial(k) for k in range(4)]  # 1, x, x^2/2, x^3/6

D_mono = matrix_of_map(D, monomial, monomial)
D_fact = matrix_of_map(D, factorial, factorial)

D_mono_expected = np.array([[0, 1, 0, 0],
                            [0, 0, 2, 0],
                            [0, 0, 0, 3],
                            [0, 0, 0, 0]], float)
shift = np.array([[0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1],
                  [0, 0, 0, 0]], float)
assert np.allclose(D_mono, D_mono_expected), D_mono
assert np.allclose(D_fact, shift), D_fact
print("Example 4.1: [D] in the monomial basis =\n", D_mono.astype(int))
print("Example 4.1: [D] in the basis {1, x, x^2/2, x^3/6} =\n", D_fact.astype(int), "(a shift)")


# Independent check with SymPy: differentiate symbolically, then read off coordinates.
def sympy_matrix_of_map(Tsym, basis_in, basis_out, x):
    """Same construction as matrix_of_map, done exactly on SymPy expressions."""
    cs = sp.symbols(f"c0:{len(basis_out)}")
    cols = []
    for b in basis_in:
        residual = sp.expand(sum(c * q for c, q in zip(cs, basis_out)) - Tsym(b))
        sol = sp.solve(sp.Poly(residual, x).all_coeffs(), cs, dict=True)[0]
        cols.append([sol.get(c, sp.Integer(0)) for c in cs])
    return sp.Matrix(cols).T


x = sp.symbols("x")
mono_sym = [sp.Integer(1), x, x**2, x**3]
fact_sym = [sp.Integer(1), x, x**2 / 2, x**3 / 6]
D_sym = lambda p: sp.diff(p, x)
assert sympy_matrix_of_map(D_sym, mono_sym, mono_sym, x) == sp.Matrix(D_mono_expected.astype(int))
assert sympy_matrix_of_map(D_sym, fact_sym, fact_sym, x) == sp.Matrix(shift.astype(int))
print("SymPy agrees with matrix_of_map in both bases.")

# Nilpotency (section 3.9): D^3 != 0 but D^4 = 0, as operators and as matrices.
p_rand = rng.standard_normal(4)
q = p_rand.copy()
for _ in range(4):
    q = D(q)
assert np.allclose(q, 0) and not np.allclose(D(D(D(p_rand))), 0)
for M in (D_mono, D_fact):
    assert not np.allclose(np.linalg.matrix_power(M, 3), 0)
    assert np.allclose(np.linalg.matrix_power(M, 4), 0)
assert np.allclose(np.linalg.matrix_power(D_mono, 3), np.array([[0, 0, 0, 6], [0] * 4, [0] * 4, [0] * 4]))
print("D^4 = 0 on P_3 (nilpotent of index 4); D^3 has the single entry 6 in the monomial basis.")

# Heatmaps of the two matrices: same operator, two bases.
fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.4))
for ax, M, name in zip(axes, (D_mono, D_fact), ("monomial basis (1, x, x², x³)", "basis (1, x, x²/2, x³/6)")):
    ax.imshow(M, cmap="Blues", vmin=0, vmax=3)
    for (i, j), val in np.ndenumerate(M):
        ax.text(j, i, f"{val:g}", ha="center", va="center", color=INK if val < 2 else "white", fontsize=11)
    ax.set_title(f"[D] in the {name}", fontsize=10)
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xlabel("input basis vector j"); ax.set_ylabel("output coordinate i")
fig.suptitle("One operator, two matrices", fontsize=11)
fig.tight_layout()
fig.savefig(FIG / "differentiation_matrices.png", dpi=120)
plt.close(fig)

# %% [markdown]
# ## Part 2: The two theorems that make matrices work
#
# Theorem 3.6, eq. (3.3): [T v]_C = [T]_{B,C} [v]_B for every v.
# Theorem 3.7, eq. (3.4): [S o T]_{B,D} = [S]_{C,D} [T]_{B,C}.
#
# Checked on D (both bases), on D o D, and on the integration map
# J: P_2 -> P_3, J(p) = integral_0^x p, for which D o J = I on P_2 but
# J o D != I on P_3 (Exercise B2).

# %%
for basis, M in ((monomial, D_mono), (factorial, D_fact)):
    for _ in range(20):
        p = rng.standard_normal(4)
        assert np.allclose(coords(D(p), basis), M @ coords(p, basis))     # (3.3)
    assert np.allclose(matrix_of_map(lambda p: D(D(p)), basis, basis), M @ M)   # (3.4)

# P_2 stored as (c_0, c_1, c_2); P_3 as (c_0, ..., c_3).
mono2 = [np.eye(3)[k] for k in range(3)]


def J(p):
    """Integration P_2 -> P_3: (c_0, c_1, c_2) -> (0, c_0, c_1/2, c_2/3)."""
    p = np.asarray(p, float)
    return np.array([0.0] + [p[k] / (k + 1) for k in range(len(p))])


def D_32(p):
    """Differentiation P_3 -> P_2 (drop the vanishing top coefficient)."""
    return D(p)[:-1]


J_mat = matrix_of_map(J, mono2, monomial)          # 4 x 3
D_32_mat = matrix_of_map(D_32, monomial, mono2)    # 3 x 4
assert np.allclose(J_mat, [[0, 0, 0], [1, 0, 0], [0, 1 / 2, 0], [0, 0, 1 / 3]])
assert np.allclose(D_32_mat @ J_mat, np.eye(3))                        # D o J = identity on P_2
assert np.allclose(J_mat @ D_32_mat, np.diag([0, 1, 1, 1]))            # J o D kills constants
assert np.allclose(matrix_of_map(lambda p: D_32(J(p)), mono2, mono2), D_32_mat @ J_mat)
print("Part 2: (3.3) and (3.4) hold; D o J = I_3 but J o D = diag(0,1,1,1).")

# Rank–nullity for D on P_3 (Example 4.1): rank 3, nullity 1, 3 + 1 = 4 = dim P_3.
rank_D = np.linalg.matrix_rank(D_mono)
null_D = null_space(D_mono)
assert rank_D == 3 and null_D.shape[1] == 1
assert np.allclose(np.abs(null_D.ravel()), [1, 0, 0, 0])   # kernel = constants
print(f"Rank–nullity for D: rank {rank_D} + nullity {null_D.shape[1]} = {rank_D + null_D.shape[1]} = dim P_3.")

# %% [markdown]
# ## Part 3: Change of basis and similarity
#
# Convention (3.5): P has the NEW basis vectors as columns, written in OLD
# coordinates, so P [v]_new = [v]_old. Theorem 3.8, eq. (3.6):
# [T]_new = P^{-1} [T]_old P. Example 4.2: rotation by 90 degrees in the basis
# {(1,1), (1,-1)}.


# %%
def change_of_basis(P, A):
    """P^{-1} A P, eq. (3.6): the matrix of the operator with matrix A, rewritten in the
    basis whose vectors are the columns of P (expressed in the old basis)."""
    return np.linalg.solve(P, A @ P)


A_rot = np.array([[0, -1], [1, 0]], float)          # rotation by +90 degrees, standard basis
P_new = np.array([[1, 1], [1, -1]], float)          # columns = new basis vectors
v = np.array([3.0, 1.0])
v_new = np.linalg.solve(P_new, v)                   # [v]_{B'} = P^{-1} v
assert np.allclose(v_new, [2, 1]) and np.allclose(P_new @ v_new, v)   # direction check
assert np.allclose(np.linalg.inv(P_new), 0.5 * P_new)                 # P^2 = 2I here

A_rot_new = change_of_basis(P_new, A_rot)
assert np.allclose(A_rot @ P_new, [[-1, 1], [1, 1]])                  # AP, shown in the notes
assert np.allclose(A_rot_new, [[0, 1], [-1, 0]])
# Second route: apply the definition directly in the new basis. Must agree.
assert np.allclose(matrix_of_map(lambda u: A_rot @ u, list(P_new.T), list(P_new.T)), A_rot_new)
assert np.linalg.det(P_new) < 0   # the new basis is left-handed: that is why the sign flipped
print("Example 4.2: P^{-1} A P =\n", A_rot_new.astype(int), " (det P = %g < 0 explains the apparent reversal)" % np.linalg.det(P_new))

# Exercise B3: the symmetric matrix [[2,1],[1,2]] is diagonal in the same basis.
assert np.allclose(change_of_basis(P_new, np.array([[2, 1], [1, 2]], float)), np.diag([3, 1]))

# Similarity invariants (Theorem 3.10 and the LA-6 / LA-7 previews) on random P.
for n in (2, 3, 4, 5):
    for _ in range(20):
        A = rng.standard_normal((n, n))
        P = rng.standard_normal((n, n))
        while abs(np.linalg.det(P)) < 1e-2:
            P = rng.standard_normal((n, n))
        B = change_of_basis(P, A)
        assert np.isclose(np.trace(B), np.trace(A))
        assert np.isclose(np.linalg.det(B), np.linalg.det(A))
        ev_A = np.sort_complex(np.linalg.eigvals(A)); ev_B = np.sort_complex(np.linalg.eigvals(B))
        assert np.allclose(ev_A, ev_B, atol=1e-6)
        # tr(XY) = tr(YX), the lemma behind trace invariance
        Y = rng.standard_normal((n, n))
        assert np.isclose(np.trace(A @ Y), np.trace(Y @ A))
print("Part 3: trace, determinant and eigenvalues are unchanged under A -> P^{-1} A P (random P, n = 2..5).")

# %% [markdown]
# ## Part 4: A projection and P^2 = P
#
# Example 4.3: projection onto the line y = 2x. Unit direction u = (1,2)/sqrt5,
# P = u u^T. Kernel = the perpendicular line through (2,-1); image = the line
# itself; in the basis {(1,2), (2,-1)} the matrix is diag(1, 0).

# %%
u = np.array([1.0, 2.0]) / np.sqrt(5)
P_proj = np.outer(u, u)
assert np.allclose(P_proj, np.array([[1, 2], [2, 4]]) / 5)
assert np.allclose(P_proj @ P_proj, P_proj)                              # P^2 = P
assert np.allclose(P_proj @ np.array([2, -1.0]), 0)                       # kernel direction
assert np.allclose(P_proj @ np.array([1, 2.0]), [1, 2])                   # fixed on its image
assert np.allclose(change_of_basis(np.array([[1, 2], [2, -1.0]]), P_proj), np.diag([1.0, 0.0]))
assert np.linalg.matrix_rank(P_proj) == 1 and null_space(P_proj).shape[1] == 1
assert np.isclose(np.trace(P_proj), 1) and np.isclose(np.linalg.det(P_proj), 0)
print("Example 4.3: P = (1/5)[[1,2],[2,4]], P^2 = P, ker = span(2,-1), im = span(1,2); diag(1,0) in the adapted basis.")

# %% [markdown]
# ## Part 5: Grid deformations, the unit circle, and areas
#
# Section 2's picture: a linear map sends the square grid to a grid of
# parallelograms (parallel, evenly spaced lines stay parallel and evenly
# spaced; the origin stays put). Section 3.9: the image of the unit square has
# area |ad - bc|, eq. (3.8). The image of the unit circle is an ellipse whose
# semi-axes are the singular values of the matrix. That is the SVD, LA-9; here
# it is only drawn and compared, not derived.


# %%
def plot_grid_deformation(ax, M, title, extent=2.0, step=0.5):
    """Draw the grid lines of [-extent, extent]^2 before (gray) and after (colored) applying M,
    together with the images of e_1, e_2 and of the unit square."""
    M = np.asarray(M, float)
    lines = np.arange(-extent, extent + step / 2, step)
    t = np.linspace(-extent, extent, 2)
    for c in lines:
        for seg in (np.vstack([np.full_like(t, c), t]), np.vstack([t, np.full_like(t, c)])):
            ax.plot(*seg, color=GRID, lw=0.8, zorder=1)
            img = M @ seg
            ax.plot(*img, color=BLUE, lw=0.9, alpha=0.8, zorder=2)
    square = np.array([[0, 1, 1, 0, 0], [0, 0, 1, 1, 0]], float)
    ax.fill(*(M @ square), color=ORANGE, alpha=0.35, zorder=3)
    ax.plot(*(M @ square), color=ORANGE, lw=1.5, zorder=3)
    for j, lab in enumerate(("T(e₁)", "T(e₂)")):
        ax.annotate("", xy=M[:, j], xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=INK, lw=1.6), zorder=4)
        ax.text(*(M[:, j] * 1.08), lab, fontsize=9, color=INK, zorder=4)
    ax.set_aspect("equal"); ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.6, 2.6)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color(GRID)


def rotation(theta):
    """Rotation of the plane by theta (counterclockwise), eq. in section 3.9."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def reflection(phi):
    """Reflection across the line through the origin at angle phi."""
    c, s = np.cos(2 * phi), np.sin(2 * phi)
    return np.array([[c, s], [s, -c]])


def shear(k):
    """Horizontal shear: (x, y) -> (x + k y, y)."""
    return np.array([[1, k], [0, 1]], float)


M_random = rng.standard_normal((2, 2))
maps = [("shear, k = 1.5", shear(1.5)),
        ("rotation by 60°", rotation(np.pi / 3)),
        ("projection onto y = 2x (rank 1)", P_proj),
        (f"random map, det = {np.linalg.det(M_random):.2f}", M_random)]
fig, axes = plt.subplots(2, 2, figsize=(8, 8))
for ax, (name, M) in zip(axes.ravel(), maps):
    plot_grid_deformation(ax, M, name)
fig.suptitle("Linear maps as grid deformations (gray: before, blue: after, orange: unit square)", fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "grid_deformations.png", dpi=120)
plt.close(fig)


def polygon_area(pts):
    """Shoelace formula for a simple polygon given as an (N, 2) array of vertices in order."""
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


unit_square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
for name, M in [("rotation 60°", rotation(np.pi / 3)), ("reflection 30°", reflection(np.pi / 6)),
                ("shear 1.5", shear(1.5)), ("projection", P_proj), ("[[2,1],[1,3]]", np.array([[2, 1], [1, 3.]]))]:
    area = polygon_area(unit_square @ M.T)
    assert np.isclose(area, abs(np.linalg.det(M)), atol=1e-12), (name, area)
    # The identity (ad - bc)^2 + (ab + cd)^2 = (a^2 + c^2)(b^2 + d^2) behind (3.8):
    (a, b), (c, d) = M
    assert np.isclose((a * d - b * c) ** 2 + (a * b + c * d) ** 2, (a**2 + c**2) * (b**2 + d**2))
assert np.isclose(polygon_area(unit_square @ np.array([[2, 1], [1, 3.]]).T), 5.0)
assert np.isclose(np.linalg.det(reflection(np.pi / 6)), -1) and np.isclose(np.linalg.det(rotation(0.7)), 1)
assert np.isclose(np.linalg.det(shear(1.5)), 1)
print("Part 5: area of the image of the unit square = |det| for rotation (1), reflection (1, det -1), shear (1), projection (0), [[2,1],[1,3]] (5).")

# Image of the unit circle: an ellipse whose semi-axes are the singular values (LA-9 preview).
th = np.linspace(0, 2 * np.pi, 400)
circle = np.vstack([np.cos(th), np.sin(th)])
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
for ax, (name, M) in zip(axes, [("shear, k = 1.5", shear(1.5)), ("projection onto y = 2x", P_proj),
                                 (f"random map", M_random)]):
    img = M @ circle
    U, s, Vt = np.linalg.svd(M)
    ax.plot(*circle, color=GRID, lw=1.2, label="unit circle")
    ax.plot(*img, color=ORANGE, lw=1.8, label="its image")
    for k in range(2):
        ax.annotate("", xy=s[k] * U[:, k], xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color=AQUA, lw=1.4))
    ax.plot([], [], color=AQUA, label="semi-axes = singular values")
    # Numerical check: the farthest and nearest points of the image from the origin are s[0], s[1].
    r = np.linalg.norm(img, axis=0)
    assert np.isclose(r.max(), s[0], atol=1e-3) and np.isclose(r.min(), s[1], atol=1e-3)
    ax.set_aspect("equal"); ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.6, 2.6)
    ax.set_title(f"{name}: σ = ({s[0]:.2f}, {s[1]:.2f})", fontsize=10)
    ax.tick_params(labelsize=8)
axes[0].legend(fontsize=8, loc="lower left")
fig.suptitle("The image of the unit circle is an ellipse (a preview of the SVD, LA-9)", fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "unit_circle_images.png", dpi=120)
plt.close(fig)
print("Unit circle images saved; semi-axes match np.linalg.svd.")

# %% [markdown]
# ## Part 6: Rank–nullity on random maps
#
# Theorem 3.14, eq. (3.7): dim V = dim ker T + dim im T. For a matrix A of
# shape m x n: n = nullity(A) + rank(A). We build matrices of a chosen rank r
# (a product of an m x r and an r x n random matrix) and check.

# %%
for _ in range(300):
    m, n = rng.integers(1, 8, size=2)
    r = int(rng.integers(0, min(m, n) + 1))
    A = rng.standard_normal((m, r)) @ rng.standard_normal((r, n)) if r > 0 else np.zeros((m, n))
    rank = np.linalg.matrix_rank(A)
    nullity = null_space(A).shape[1]
    assert rank == r, (rank, r)
    assert rank + nullity == n, (m, n, r, rank, nullity)
    # Corollary 3.15 for square matrices: injective (nullity 0) iff surjective (rank m).
    if m == n:
        assert (nullity == 0) == (rank == m)
print("Part 6: rank + nullity = number of columns on 300 random maps of chosen rank (m, n <= 7).")

# %% [markdown]
# ## Part 7: The 2D stress tensor in a rotated frame (Example 4.4)
#
# The Cauchy stress tensor is a linear map n -> t (unit normal to traction).
# Its matrix in a frame rotated by theta is sigma' = R^T sigma R with
# R = rotation(theta), whose columns are the new axes in old coordinates:
# this is eq. (3.6) with P = R and P^{-1} = R^T. As theta varies, the pair
# (sigma'_x, tau') traces Mohr's circle.

# %%
sx, sy, txy = 80.0, 20.0, 40.0          # MPa
sigma = np.array([[sx, txy], [txy, sy]])
centre, radius = (sx + sy) / 2, np.hypot((sx - sy) / 2, txy)
assert np.isclose(centre, 50) and np.isclose(radius, 50)

theta_p = 0.5 * np.arctan2(2 * txy, sx - sy)
assert np.isclose(np.degrees(theta_p), 26.565, atol=1e-3)
assert np.isclose(np.cos(2 * theta_p), 3 / 5) and np.isclose(np.sin(2 * theta_p), 4 / 5)
sigma_p = change_of_basis(rotation(theta_p), sigma)
assert np.allclose(sigma_p, np.diag([100.0, 0.0]), atol=1e-9)
assert np.allclose(change_of_basis(rotation(np.pi / 4), sigma), [[90, -30], [-30, 10]])
assert np.isclose(np.linalg.det(sigma), 0)      # one principal stress is zero: det is invariant
print("Example 4.4: at theta_p = %.2f deg the stress matrix is diag(100, 0) MPa; at 45 deg it is [[90,-30],[-30,10]]." % np.degrees(theta_p))

thetas = np.linspace(0, np.pi, 361)
sig_x, tau = [], []
for th_ in thetas:
    S = change_of_basis(rotation(th_), sigma)
    assert np.isclose(np.trace(S), sx + sy) and np.isclose(np.linalg.det(S), np.linalg.det(sigma))
    assert np.allclose(S, S.T)                    # symmetry survives the change of frame
    sig_x.append(S[0, 0]); tau.append(S[0, 1])
sig_x, tau = np.array(sig_x), np.array(tau)
# Every rotated (sigma'_x, tau') lies on the circle: centre (50, 0), radius 50.
assert np.allclose(np.hypot(sig_x - centre, tau), radius)
# Closed forms used in the notes:
assert np.allclose(sig_x, centre + (sx - sy) / 2 * np.cos(2 * thetas) + txy * np.sin(2 * thetas))
assert np.allclose(tau, -(sx - sy) / 2 * np.sin(2 * thetas) + txy * np.cos(2 * thetas))

fig, ax = plt.subplots(figsize=(5.2, 5))
ax.plot(sig_x, tau, color=BLUE, lw=1.8, label="(σ′ₓ, τ′) as θ runs 0…180°")
ax.plot([sx], [txy], "o", color=ORANGE, ms=7, label=f"θ = 0: ({sx:g}, {txy:g})")
ax.plot([100, 0], [0, 0], "s", color=AQUA, ms=7, label="principal stresses 100, 0")
ax.plot(sigma_p[0, 0], sigma_p[0, 1], "x", color=INK, ms=9, mew=2, label=f"θ = θₚ = {np.degrees(theta_p):.1f}°")
ax.axhline(0, color=GRID, lw=0.8); ax.axvline(0, color=GRID, lw=0.8)
ax.set_aspect("equal"); ax.grid(alpha=.25)
ax.set_xlabel("normal stress σ′ₓ (MPa)"); ax.set_ylabel("shear stress τ′ (MPa)")
ax.set_title("Mohr's circle is a change of basis", fontsize=10)
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout()
fig.savefig(FIG / "mohr_circle.png", dpi=120)
plt.close(fig)
print("Mohr's circle saved: centre (50, 0), radius 50.")

# %% [markdown]
# ## Your turn
#
# 1. **The transpose map.** Let T(A) = A^T on 2x2 matrices, stored flattened
#    row-major as (a11, a12, a21, a22). Use `matrix_of_map` with the basis
#    E11, E12, E21, E22 to get its 4x4 matrix. Check that it is a permutation
#    matrix and that T^2 = I. What are its kernel and image?
# 2. **Kernel and image of A -> A - A^T.** Same storage. Compute the 4x4
#    matrix, its rank and nullity with `np.linalg.matrix_rank` and
#    `null_space`, and identify the kernel (symmetric matrices) and image
#    (skew-symmetric matrices) by inspecting the basis vectors returned.
#    Confirm 3 + 1 = 4.
# 3. **Rotation about the z-axis in a tilted basis.** R_z(30°) in the standard
#    basis is block-diagonal. Take the tilted basis b1 = (1,0,0),
#    b2 = (0, cos45°, sin45°), b3 = (0, -sin45°, cos45°), build P from its
#    columns, and compute `change_of_basis(P, Rz)`. It is a full 3x3 matrix
#    with no zeros. Check the trace is still 1 + 2cos30° and that the result
#    is orthogonal (P is orthogonal here, so P^{-1} = P^T). Then redo it with
#    a non-orthogonal tilt such as b2 = (0, 1, 1): the trace survives but
#    orthogonality of the matrix does not. Why?
# 4. **A linearity check.** `matrix_of_map` assumes T is linear. Add a function
#    `is_linear(T, basis, trials=20)` that tests T(a u + b v) = a T(u) + b T(v)
#    on random a, b and random combinations of basis vectors. Run it on D, on
#    the transpose map, and on the non-linear map p -> p^2 (coefficient-wise
#    square). Explain why `matrix_of_map` on the non-linear map returns a
#    matrix anyway, and what that matrix predicts wrongly.

# %% [markdown]
# ### Reference checks for "Your turn" 1–3 (compare after you have tried)

# %%
E = [np.eye(4)[k] for k in range(4)]            # E11, E12, E21, E22 flattened row-major
transpose_map = lambda a: a.reshape(2, 2).T.ravel()
T_mat = matrix_of_map(transpose_map, E, E)
assert np.allclose(T_mat, [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
assert np.allclose(T_mat @ T_mat, np.eye(4)) and np.linalg.matrix_rank(T_mat) == 4

skew_part = lambda a: a - transpose_map(a)
K_mat = matrix_of_map(skew_part, E, E)
assert np.allclose(K_mat, np.eye(4) - T_mat)
assert np.linalg.matrix_rank(K_mat) == 1 and null_space(K_mat).shape[1] == 3
ker_basis = null_space(K_mat)
for k in range(3):   # every kernel vector is a symmetric matrix
    Kk = ker_basis[:, k].reshape(2, 2)
    assert np.allclose(Kk, Kk.T)
img = K_mat @ rng.standard_normal(4)
assert np.allclose(img.reshape(2, 2), -img.reshape(2, 2).T)   # the image is skew-symmetric
print("Your turn 1–2: transpose map is a 4x4 permutation with T^2 = I; A -> A - A^T has rank 1, nullity 3.")

th30, a45 = np.radians(30), np.radians(45)
Rz = np.array([[np.cos(th30), -np.sin(th30), 0], [np.sin(th30), np.cos(th30), 0], [0, 0, 1]])
P_tilt = np.array([[1, 0, 0], [0, np.cos(a45), -np.sin(a45)], [0, np.sin(a45), np.cos(a45)]])
Rz_tilt = change_of_basis(P_tilt, Rz)
assert np.isclose(np.trace(Rz_tilt), 1 + 2 * np.cos(th30))
assert np.allclose(Rz_tilt @ Rz_tilt.T, np.eye(3))              # orthogonal P keeps orthogonality
assert np.count_nonzero(np.abs(Rz_tilt) > 1e-12) == 9           # no zeros left
P_skew = np.array([[1, 0, 0], [0, 1, -1], [0, 1, 1.]])          # non-orthogonal tilt (b2 = (0,1,1) is not unit)
Rz_skew = change_of_basis(P_skew, Rz)
assert np.isclose(np.trace(Rz_skew), 1 + 2 * np.cos(th30))
print("Your turn 3: R_z(30°) in the tilted basis has trace %.4f = 1 + 2cos30° and no zero entries." % np.trace(Rz_tilt))

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
