# %% [markdown]
# # LA-1 Project: The Two Pictures, Span, and the Dot Product
#
# This script builds the geometric intuition of LA-1 and checks every
# numerical claim in the notes. Parts: (1) row vs column picture of Ax = b,
# (2) a span explorer in R^3, (3) Ax computed two ways, (3.3) checked column
# by column, and the non-commuting products of Example 4.3, (4) Cauchy–Schwarz,
# the angle formula, and the triangle inequality on random data, (5) unit
# spheres of three norms. A reference solution for exercise B8 is the last cell.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# %% [markdown]
# ## Part 1: Row picture vs column picture
#
# Example 4.1: 2x - y = 1, x + y = 5. Row picture: two lines meeting at (2, 3).
# Column picture: 2*(2,1) + 3*(-1,1) = (1,5). Both are drawn below.


# %%
def plot_two_pictures(A, b, name):
    """Draw the row picture (lines) and the column picture (vector recipe) of Ax = b."""
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    x = np.linalg.solve(A, b)  # unique solution assumed for the picture

    fig, (ax_row, ax_col) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Row picture: each row a_i1 x + a_i2 y = b_i is a line.
    t = np.linspace(x[0] - 6, x[0] + 6, 200)
    for i in range(2):
        a1, a2 = A[i]
        if abs(a2) > 1e-12:
            ax_row.plot(t, (b[i] - a1 * t) / a2, label=f"row {i+1}: {a1:g}x + {a2:g}y = {b[i]:g}")
        else:
            ax_row.axvline(b[i] / a1, label=f"row {i+1}")
    ax_row.plot(*x, "ko", ms=7, label=f"solution ({x[0]:g}, {x[1]:g})")
    ax_row.set_title("Row picture: intersect the lines")
    ax_row.set_aspect("equal"); ax_row.grid(alpha=.3); ax_row.legend(fontsize=8)
    ax_row.set_xlim(x[0] - 4, x[0] + 4); ax_row.set_ylim(x[1] - 4, x[1] + 4)  # centred on the solution

    # Column picture: x1 * a_1 + x2 * a_2 = b, drawn head to tail.
    a1, a2 = A[:, 0], A[:, 1]
    ax_col.quiver(0, 0, *a1, angles="xy", scale_units="xy", scale=1, color="C0", label="column 1")
    ax_col.quiver(0, 0, *a2, angles="xy", scale_units="xy", scale=1, color="C1", label="column 2")
    ax_col.quiver(0, 0, *(x[0] * a1), angles="xy", scale_units="xy", scale=1, color="C0", alpha=.4, width=.012,
                  label=f"{x[0]:g} × column 1")
    ax_col.quiver(*(x[0] * a1), *(x[1] * a2), angles="xy", scale_units="xy", scale=1, color="C1", alpha=.4, width=.012,
                  label=f"{x[1]:g} × column 2")
    ax_col.quiver(0, 0, *b, angles="xy", scale_units="xy", scale=1, color="k", width=.006, label="b")
    ax_col.set_title("Column picture: combine the columns to hit b")
    ax_col.set_aspect("equal"); ax_col.grid(alpha=.3); ax_col.legend(fontsize=8)
    lim = max(np.abs(b).max(), np.abs(A).max() * np.abs(x).max()) + 1
    ax_col.set_xlim(-lim, lim); ax_col.set_ylim(-lim, lim)

    fig.tight_layout()
    fig.savefig(FIG / f"two_pictures_{name}.png", dpi=120)
    plt.close(fig)
    return x


A1 = [[2, -1], [1, 1]]
b1 = [1, 5]
x1 = plot_two_pictures(A1, b1, "example_4_1")
assert np.allclose(x1, [2, 3]), x1
# The solution IS the recipe: x1[0] * col1 + x1[1] * col2 == b.
assert np.allclose(x1[0] * np.array(A1)[:, 0] + x1[1] * np.array(A1)[:, 1], b1)
print("Example 4.1 solution:", x1)

A_rand = rng.integers(-3, 4, size=(2, 2)).astype(float)
while abs(np.linalg.det(A_rand)) < 0.5:
    A_rand = rng.integers(-3, 4, size=(2, 2)).astype(float)
b_rand = rng.integers(-4, 5, size=2).astype(float)
x_rand = plot_two_pictures(A_rand, b_rand, "random")
print("Random system A =", A_rand.tolist(), "b =", b_rand.tolist(), "x =", x_rand)

# Example 4.2 (and Theorem 3.6): parallel columns span only a line.
A2 = np.array([[1, 2], [2, 4.]])
assert np.linalg.matrix_rank(A2) == 1                      # columns (1,2), (2,4) are parallel
assert np.linalg.matrix_rank(np.c_[A2, [1, 3]]) == 2       # b = (1,3) is off their line: no solution
assert np.linalg.matrix_rank(np.c_[A2, [1, 2]]) == 1       # b = (1,2) is on it: solutions exist
assert np.allclose(A2 @ [1, 0], [1, 2]) and np.allclose(A2 @ [0, 0.5], [1, 2])  # two different recipes
print("Example 4.2: rank A = 1; (1,3) not in the column span, (1,2) is, with recipes (1,0) and (0,1/2).")

# %% [markdown]
# ## Part 2: Span explorer
#
# Sample random coefficients c and plot sum_i c_i v_i. For two independent
# vectors in R^3 you get a plane; for three independent vectors you fill space;
# for the dependent triple of exercise B4 you still get only a plane.


# %%
def sample_span(vectors, n=3000, scale=2.0):
    """Return n random linear combinations of the given vectors (rows of V)."""
    V = np.asarray(vectors, float)
    C = rng.uniform(-scale, scale, size=(n, V.shape[0]))
    return C @ V  # row i is sum_k C[i,k] * V[k]


cases = {
    "two_vectors_plane": [[1, 0, 1], [0, 1, 1]],
    "three_independent_space": [[1, 0, 1], [0, 1, 1], [1, 1, 0]],
    "three_dependent_plane_B4": [[1, 0, 1], [0, 1, 1], [1, 1, 2]],
}
fig = plt.figure(figsize=(13, 4))
for k, (name, vecs) in enumerate(cases.items(), 1):
    P = sample_span(vecs)
    ax = fig.add_subplot(1, 3, k, projection="3d")
    ax.scatter(P[:, 0], P[:, 1], P[:, 2], s=2, alpha=.3)
    for v in vecs:
        ax.quiver(0, 0, 0, *v, color="r", arrow_length_ratio=.1)
    ax.set_title(name.replace("_", " "))
    ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_zlim(-4, 4)
fig.tight_layout()
fig.savefig(FIG / "span_explorer.png", dpi=120)
plt.close(fig)

# Numerical check of B4(b): every combination of the dependent triple satisfies x + y - z = 0.
P_dep = sample_span(cases["three_dependent_plane_B4"])
n_normal = np.array([1, 1, -1.0])
assert np.allclose(P_dep @ n_normal, 0), "dependent triple should stay in the plane x + y - z = 0"
# ... while the independent triple leaves that plane,
P_ind = sample_span(cases["three_independent_space"])
assert np.abs(P_ind @ n_normal).max() > 1.0
# and the actual claims: the independent triple spans all of R^3 (rank 3), the dependent one a plane (rank 2).
assert np.linalg.matrix_rank(np.array(cases["three_independent_space"], float)) == 3
assert np.linalg.matrix_rank(np.array(cases["three_dependent_plane_B4"], float)) == 2
assert np.linalg.matrix_rank(np.array(cases["two_vectors_plane"], float)) == 2
print("Span checks passed: dependent triple has rank 2 and lies in x + y - z = 0; independent triple has rank 3.")

# %% [markdown]
# ## Part 3: Ax two ways, and AB from its columns
#
# (3.1) column form: Ax = sum_j x_j a_j.   (3.2) row form: (Ax)_i = r_i · x.
# (3.3) AB has k-th column A b_k.  Theorem 3.5: A e_j = a_j.
#
# **Build it yourself first** (notes §7): write `matvec_columns`, `matvec_rows`,
# and Part 4's `project_onto_line` before reading these. Exercise B8 then asks
# you to add `matmat_columns(A, B)`; this part uses `A @ B` so as not to give it away.


# %%
def matvec_columns(A, x):
    """Ax as a linear combination of the columns of A, eq. (3.1)."""
    A = np.asarray(A, float); x = np.asarray(x, float)
    m, n = A.shape
    assert x.shape == (n,), "x must have one entry per column of A"
    out = np.zeros(m)
    for j in range(n):
        out += x[j] * A[:, j]
    return out


def matvec_rows(A, x):
    """Ax entry by entry as dot products of rows with x, eq. (3.2)."""
    A = np.asarray(A, float); x = np.asarray(x, float)
    m, n = A.shape
    out = np.zeros(m)
    for i in range(m):
        out[i] = sum(A[i, j] * x[j] for j in range(n))
    return out


for _ in range(50):
    m, n, p = rng.integers(1, 7, size=3)
    A = rng.standard_normal((m, n)); x = rng.standard_normal(n); B = rng.standard_normal((n, p))
    assert np.allclose(matvec_columns(A, x), matvec_rows(A, x))
    assert np.allclose(matvec_columns(A, x), A @ x)
    # (3.3): the k-th column of AB is A times the k-th column of B.
    for k in range(p):
        assert np.allclose((A @ B)[:, k], matvec_columns(A, B[:, k]))
    # Theorem 3.5: A e_j = column j.
    for j in range(n):
        e = np.zeros(n); e[j] = 1
        assert np.allclose(matvec_columns(A, e), A[:, j])
    # (AB)x = A(Bx), the reason (3.3) is the right definition.
    y = rng.standard_normal(p)
    assert np.allclose((A @ B) @ y, A @ (B @ y))
    # (AB)^T = B^T A^T, Theorem 3.8.
    assert np.allclose((A @ B).T, B.T @ A.T)

# Example 4.3: rotation and stretch do not commute.
R = np.array([[0, -1], [1, 0]], float)
S = np.array([[2, 0], [0, 1]], float)
RS = R @ S; SR = S @ R
for k in range(2):  # the column-by-column computation in the notes
    assert np.allclose(RS[:, k], matvec_columns(R, S[:, k])) and np.allclose(SR[:, k], matvec_columns(S, R[:, k]))
assert np.allclose(RS, [[0, -1], [2, 0]]) and np.allclose(SR, [[0, -2], [1, 0]])
assert not np.allclose(RS, SR)
assert np.allclose(RS @ [1, 0], [0, 2]) and np.allclose(SR @ [1, 0], [0, 1])
print("Example 4.3: RS =", RS.tolist(), " SR =", SR.tolist(), " (not equal)")

# Example 4.3, rigid body: 90° turns about the fixed x- and y-axes, applied to a phone's screen normal e_3.
Rx = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]], float)
Ry = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]], float)
c90, s90 = np.cos(np.pi / 2), np.sin(np.pi / 2)
assert np.allclose(Rx, [[1, 0, 0], [0, c90, -s90], [0, s90, c90]])   # standard rotation about x by 90°
assert np.allclose(Ry, [[c90, 0, s90], [0, 1, 0], [-s90, 0, c90]])   # standard rotation about y by 90°
e3 = np.array([0, 0, 1.0])
assert np.allclose(Ry @ Rx, [[0, 1, 0], [0, 0, -1], [-1, 0, 0]]) and np.allclose(Ry @ Rx @ e3, [0, -1, 0])
assert np.allclose(Rx @ Ry, [[0, 0, 1], [1, 0, 0], [0, 1, 0]]) and np.allclose(Rx @ Ry @ e3, [1, 0, 0])
print("Example 4.3 (3D): x-then-y sends e3 to", (Ry @ Rx @ e3).tolist(), "; y-then-x sends it to", (Rx @ Ry @ e3).tolist())
print("Part 3 checks passed: (3.1) == (3.2) == NumPy; (3.3) column by column; (AB)x = A(Bx).")

# %% [markdown]
# ## Part 4: Dot product, angle formula, Cauchy–Schwarz, triangle inequality
#
# Theorem 3.10: u·v = |u||v|cos θ in R^2 and R^3.  Theorem 3.11: |u·v| ≤ |u||v|.
# (3.5): cos θ = u·v / (|u||v|) defines θ in any R^n.  Theorem 3.12: |u+v| ≤ |u| + |v|.
# §3.6: the l1 and l-inf triangle inequalities.  Example 4.5 / B7: projection formula.


# %%
def dot(u, v):
    """Dot product u·v = sum_i u_i v_i, Definition 3.9."""
    return float(np.sum(np.asarray(u, float) * np.asarray(v, float)))


def norm(v):
    """Euclidean norm ||v|| = sqrt(v·v), Definition 3.9."""
    return np.sqrt(dot(v, v))


def angle_between(u, v):
    """θ in [0, π] (radians) with cos θ = u·v / (||u|| ||v||), eq. (3.5)."""
    c = dot(u, v) / (norm(u) * norm(v))
    return np.arccos(np.clip(c, -1, 1))


def project_onto_line(u, v):
    """Projection of u onto the line through v: t* v with t* = (u·v)/(v·v)."""
    t = dot(u, v) / dot(v, v)
    return t * np.asarray(v, float)


N = 10_000
for n in (2, 3, 4, 10):
    U = rng.standard_normal((N, n)); V = rng.standard_normal((N, n))
    uv = np.sum(U * V, axis=1); nu = np.linalg.norm(U, axis=1); nv = np.linalg.norm(V, axis=1)
    assert np.all(np.abs(uv) <= nu * nv * (1 + 1e-12)), "Cauchy–Schwarz violated"
    assert np.all(np.linalg.norm(U + V, axis=1) <= nu + nv + 1e-12), "triangle inequality violated"
    # §3.6: triangle inequality for l1 and l-inf (entrywise |a + b| <= |a| + |b|, then sum / take max).
    assert np.all(np.abs(U + V).sum(1) <= np.abs(U).sum(1) + np.abs(V).sum(1) + 1e-12)
    assert np.all(np.abs(U + V).max(1) <= np.abs(U).max(1) + np.abs(V).max(1) + 1e-12)
    # Parallelogram law (B6) holds for the Euclidean norm...
    lhs = np.linalg.norm(U + V, axis=1) ** 2 + np.linalg.norm(U - V, axis=1) ** 2
    assert np.allclose(lhs, 2 * nu**2 + 2 * nv**2)
    # ... and fails for the l1 norm (so l1 is not an inner-product norm).
    l1 = lambda M: np.abs(M).sum(axis=1)
    assert not np.allclose(l1(U + V) ** 2 + l1(U - V) ** 2, 2 * l1(U) ** 2 + 2 * l1(V) ** 2)

# Example 4.4 numbers.
u = [1, 2, 0, -1]; v = [3, -1, 2, 1]; w = [1, 1, 1, 1]
assert dot(u, v) == 0
assert np.isclose(norm(u), np.sqrt(6)) and np.isclose(norm(v), np.sqrt(15))
assert dot(u, w) == 2 and np.isclose(norm(w), 2)
theta_uw = np.degrees(angle_between(u, w))
assert np.isclose(theta_uw, 65.9, atol=0.05), theta_uw          # the value printed in Example 4.4
print(f"Example 4.4: u·v = 0 (orthogonal); angle(u, w) = {theta_uw:.1f} degrees")

# Theorem 3.10 checked independently in R^2: the geometric angle between two arrows is the
# difference of their polar angles (np.arctan2), folded into [0, π]. No dot products involved.
U2 = rng.standard_normal((N, 2)); V2 = rng.standard_normal((N, 2))
d = np.abs(np.arctan2(U2[:, 1], U2[:, 0]) - np.arctan2(V2[:, 1], V2[:, 0]))
theta_geom = np.minimum(d, 2 * np.pi - d)
theta_dot = np.array([angle_between(a, b) for a, b in zip(U2, V2)])
# arccos loses digits near 0 and π (its slope blows up), so allow 1e-9 rather than 1e-15.
assert np.allclose(theta_dot, theta_geom, rtol=0, atol=1e-9)
print(f"Theorem 3.10 in R^2: dot-product angle matches the arctan2 angle on {N} random pairs "
      f"(max diff {np.abs(theta_dot - theta_geom).max():.1e})")

# Example 4.5 and B7(d): projection and the orthogonal residual.
p = project_onto_line([3, 1], [1, 1]); assert np.allclose(p, [2, 2])
assert np.isclose(dot(np.array([3, 1]) - p, [1, 1]), 0)
u = np.array([2, 3, 6.]); v = np.array([1, 1, 1.])
p = project_onto_line(u, v); r = u - p
assert np.allclose(p, 11 / 3) and np.isclose(dot(r, v), 0)
assert np.isclose(norm(u) ** 2, norm(p) ** 2 + norm(r) ** 2)  # Pythagoras, B7(c)
print("Example 4.5 / B7(d): projection =", p, " residual·v =", round(dot(r, v), 12))

# Cauchy–Schwarz equality holds exactly when the vectors are parallel.
a = rng.standard_normal(5); b = -2.5 * a
assert np.isclose(abs(dot(a, b)), norm(a) * norm(b))
print("Part 4 checks passed on", N, "random pairs in several dimensions.")

# %% [markdown]
# ## Part 5: Unit balls of the l1, l2, l-infinity norms
#
# The figure draws each unit sphere {v : ||v|| = 1}, the boundary of the unit ball {v : ||v|| <= 1}.

# %%
th = np.linspace(0, 2 * np.pi, 400)
circle = np.c_[np.cos(th), np.sin(th)]
fig, ax = plt.subplots(figsize=(5, 5))
for q, (name, f) in enumerate([("l1 (diamond)", lambda P: np.abs(P).sum(1)),
                               ("l2 (disk)", lambda P: np.linalg.norm(P, axis=1)),
                               ("l-inf (square)", lambda P: np.abs(P).max(1))]):
    # Scale each direction so its norm is exactly 1: that traces the unit sphere of the norm.
    boundary = circle / f(circle)[:, None]
    ax.plot(boundary[:, 0], boundary[:, 1], label=name)
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.legend(); ax.set_title("Unit spheres {v : ||v|| = 1}\n(boundaries of the unit balls)")
fig.savefig(FIG / "unit_balls.png", dpi=120)
plt.close(fig)

# C2's chain holds on random vectors: ||v||_inf <= ||v||_2 <= ||v||_1 <= sqrt(n)||v||_2 <= n||v||_inf
V = rng.standard_normal((N, 6)); n = 6
linf = np.abs(V).max(1); l2 = np.linalg.norm(V, axis=1); l1 = np.abs(V).sum(1)
assert np.all(linf <= l2 + 1e-12) and np.all(l2 <= l1 + 1e-12)
assert np.all(l1 <= np.sqrt(n) * l2 + 1e-12) and np.all(np.sqrt(n) * l2 <= n * linf + 1e-12)
print("Part 5: norm chain (C2) verified; unit spheres saved.")

# %% [markdown]
# ## Your turn
#
# 1. In Part 1, pick an A whose columns are parallel and a b off their line.
#    `np.linalg.solve` will fail. Catch the error and instead draw the column
#    picture to *show* why there is no recipe. Then move b onto the line.
# 2. In Part 2, add a fourth case: four random vectors in R^3. What does the
#    span look like? Why can it never be more than all of R^3?
# 3. Write `matvec_rows` using a single NumPy expression (no Python loops).
#    Hint: `(A * x).sum(axis=1)`. Why does broadcasting make this correct?
# 4. (C3) Implement the Frobenius norm and check `||A x|| <= ||A||_F ||x||`
#    on random matrices. Then estimate the *smallest* constant K by maximizing
#    `||A x||` over 10,000 random unit vectors x. Compare K with ||A||_F.
#    (LA-9 will tell you K is the largest singular value.)
# 5. Exercise B8: add `matmat_columns(A, B)` here, built from `matvec_columns`.
#    A reference solution is the next cell. Write your own first, then delete
#    the reference cell (it would otherwise redefine your function) and rerun.


# %% [markdown]
# ## Reference solution for B8 (write your own first)


# %%
def matmat_columns(A, B):
    """AB whose k-th column is A @ B[:, k], eq. (3.3). Reference solution for exercise B8."""
    A = np.asarray(A, float); B = np.asarray(B, float)
    m, n = A.shape; n2, p = B.shape
    assert n == n2, "inner dimensions must match"
    C = np.zeros((m, p))
    for k in range(p):
        C[:, k] = matvec_columns(A, B[:, k])
    return C


A = rng.standard_normal((4, 3)); B = rng.standard_normal((3, 5))
assert np.allclose(matmat_columns(A, B), A @ B)
x = rng.standard_normal(5)
assert np.allclose(matmat_columns(A, B) @ x, A @ (B @ x))
# Why (AB)x = A(Bx) is automatic: by (3.1), (AB)x = sum_k x_k (column k of AB) = sum_k x_k (A b_k).
# By linearity of A (Theorem 3.4) that is A(sum_k x_k b_k) = A(Bx). Definition (3.3) was chosen to
# make this one-line argument work; "row times column" is the entry-level consequence.
print("B8 reference: matmat_columns agrees with A @ B on 4x3 times 3x5.")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
