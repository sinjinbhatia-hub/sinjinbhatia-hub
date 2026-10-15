# %% [markdown]
# # LA-6 Project: Determinants Three Ways, Volume Scaling, and the Jacobian
#
# This script implements the determinant three ways (elimination with swap
# tracking, the Leibniz sum over permutations, and recursive cofactor expansion),
# checks them against each other and against NumPy, and then tests every
# numerical claim in the notes: the identities of §3.4, the adjugate and Cramer
# (§3.6–3.7), the cross product (§3.8), the Vandermonde formula (§3.12), volume
# scaling by |det A| via Monte Carlo (§3.10), the Jacobians of polar, cylindrical
# and spherical coordinates (§3.11), and the scale-dependence that makes det
# useless as a measure of "how singular" (§3.13). It also times cofactor against
# elimination so you can watch n! happen.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import itertools
import math
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from scipy.spatial import ConvexHull

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# Fixed categorical colors (series 1-3) and neutral ink; never cycled.
BLUE, ORANGE, AQUA, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"

# %% [markdown]
# ## Part 1: Three algorithms for det A
#
# (3.7) elimination: det A = (-1)^s * product of pivots, s = number of row swaps.
# (3.9) Leibniz: det A = sum over permutations sigma of sgn(sigma) prod_i a_{i,sigma(i)},
#       with sgn(sigma) = (-1)^(number of inversions).
# (3.12) Laplace: det A = sum_j (-1)^(1+j) a_{1j} M_{1j}, recursively.
#
# **Build it yourself first:** write `det_elimination` and `det_leibniz` before
# reading these.


# %%
def inversions(perm):
    """N(sigma): the number of pairs i < j with perm[i] > perm[j] (§3.3)."""
    n = len(perm)
    return sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])


def sign(perm):
    """sgn(sigma) = (-1)^N(sigma)."""
    return -1 if inversions(perm) % 2 else 1


def leibniz_terms(A):
    """The n! signed products of the Leibniz formula (3.9), keyed by permutation (0-based)."""
    A = np.asarray(A, float)
    n = A.shape[0]
    terms = {}
    for perm in itertools.permutations(range(n)):
        term = float(sign(perm))
        for i, j in enumerate(perm):
            term *= A[i, j]
        terms[perm] = float(term)
    return terms


def det_leibniz(A):
    """det A by the Leibniz formula (3.9). n! terms, so refuse past n = 7."""
    A = np.asarray(A, float)
    assert A.shape[0] <= 7, "n! terms: use det_elimination for n > 7"
    return sum(leibniz_terms(A).values())


def det_cofactor(A, count=False):
    """det A by recursive Laplace expansion along the first row, eq. (3.12).

    Works on index tuples instead of copying submatrices, and stops the recursion
    at 2x2. The sign (-1)^(1+j) is applied by adding or subtracting, so the only
    multiplications are the products a_1j * M_1j: T(1) = 0, T(n) = n T(n-1) + n,
    i.e. T(n) = sum_{k=1}^{n-1} n!/k!  ~  (e - 1) n!  (§3.5). Zero entries are
    skipped, so the count is exact for matrices with no zeros.
    With count=True, returns (det, number of multiplications)."""
    M = [[float(x) for x in row] for row in np.asarray(A)]
    mults = [0]

    def rec(rows, cols):
        n = len(cols)
        if n == 1:
            return M[rows[0]][cols[0]]
        if n == 2:
            (r0, r1), (c0, c1) = rows, cols
            mults[0] += 2
            return M[r0][c0] * M[r1][c1] - M[r0][c1] * M[r1][c0]
        i, rest = rows[0], rows[1:]
        total = 0.0
        for k, j in enumerate(cols):
            a = M[i][j]
            if a != 0.0:
                mults[0] += 1
                term = a * rec(rest, cols[:k] + cols[k + 1:])
                total = total + term if k % 2 == 0 else total - term
        return total

    n = len(M)
    d = rec(tuple(range(n)), tuple(range(n)))
    return (d, mults[0]) if count else d


def det_elimination(A, pivoting="partial"):
    """det A = (-1)^s * prod(pivots), eq. (3.7), by Gaussian elimination.

    pivoting="partial" picks the largest |entry| in the column (what LAPACK does);
    pivoting="first_nonzero" picks the first nonzero entry (what you do by hand).
    Returns (det, pivots, swaps)."""
    U = np.array(A, float)
    n = U.shape[0]
    # "No pivot" must be judged relative to the scale of A, never by an absolute
    # threshold: det(1e-15 * I_2) = 1e-30 is a perfectly invertible matrix (§3.13, §8).
    tol = n * np.finfo(float).eps * (np.abs(U).max() if U.size else 0.0)
    swaps, pivots = 0, []
    for k in range(n):
        col = U[k:, k]
        if pivoting == "partial":
            p = k + int(np.argmax(np.abs(col)))
        else:
            nz = np.flatnonzero(np.abs(col) > tol)
            p = k + int(nz[0]) if nz.size else k
        if abs(U[p, k]) <= tol:
            return 0.0, [float(x) for x in pivots] + [0.0] * (n - k), swaps   # no pivot: singular
        if p != k:
            U[[k, p]] = U[[p, k]]
            swaps += 1
        pivots.append(U[k, k])
        for i in range(k + 1, n):
            U[i, k:] -= (U[i, k] / U[k, k]) * U[k, k:]
    return float((-1) ** swaps * np.prod(pivots)), [float(x) for x in pivots], swaps


def cofactor(A, i, j):
    """C_ij = (-1)^(i+j) * det(A with row i and column j deleted)."""
    minor = np.delete(np.delete(np.asarray(A, float), i, axis=0), j, axis=1)
    return (-1) ** (i + j) * det_elimination(minor)[0]


# All three agree with NumPy on random matrices (dense and integer).
for n in range(1, 8):
    for _ in range(10):
        A = rng.standard_normal((n, n))
        ref = np.linalg.det(A)
        assert np.isclose(det_elimination(A)[0], ref, rtol=1e-9, atol=1e-12)
        assert np.isclose(det_elimination(A, "first_nonzero")[0], ref, rtol=1e-9, atol=1e-12)
        assert np.isclose(det_cofactor(A), ref, rtol=1e-9, atol=1e-12)
        if n <= 6:
            assert np.isclose(det_leibniz(A), ref, rtol=1e-9, atol=1e-12)
for n in range(8, 13):
    A = rng.standard_normal((n, n))
    assert np.isclose(det_elimination(A)[0], np.linalg.det(A), rtol=1e-9)
# Scale: a tiny but perfectly invertible matrix is not "singular" (§3.13).
assert np.isclose(det_elimination(1e-15 * np.eye(2))[0], 1e-30, rtol=1e-12, atol=0)
assert det_elimination([[1.0, 2.0], [2.0, 4.0]])[0] == 0.0                 # genuinely singular
print("Part 1: det_elimination, det_cofactor, det_leibniz agree with np.linalg.det.")

# Example 4.1, all four ways.
A41 = np.array([[0, 2, 1], [1, 1, 3], [2, -1, 1]], float)
d_hand, piv_hand, s_hand = det_elimination(A41, "first_nonzero")
assert s_hand == 1 and np.allclose(piv_hand, [1, 2, -3.5]) and np.isclose(d_hand, 7)
d_pp, piv_pp, s_pp = det_elimination(A41, "partial")
assert s_pp == 2 and np.allclose(piv_pp, [2, 2, 1.75]) and np.isclose(d_pp, 7)
print(f"Example 4.1 by hand-style elimination: pivots {piv_hand}, swaps {s_hand}, det {d_hand:g}")
print(f"Example 4.1 by partial pivoting:        pivots {piv_pp}, swaps {s_pp}, det {d_pp:g}")

row1 = sum(A41[0, j] * cofactor(A41, 0, j) for j in range(3))
col1 = sum(A41[i, 0] * cofactor(A41, i, 0) for i in range(3))
assert np.isclose(cofactor(A41, 0, 1), 5) and np.isclose(cofactor(A41, 0, 2), -3)
assert np.isclose(cofactor(A41, 1, 0), -3) and np.isclose(cofactor(A41, 2, 0), 5)
assert np.isclose(row1, 7) and np.isclose(col1, 7)

terms = leibniz_terms(A41)
expected = {(0, 1, 2): 0, (0, 2, 1): 0, (1, 0, 2): -2, (1, 2, 0): 12, (2, 0, 1): -1, (2, 1, 0): -2}
for perm, val in expected.items():
    assert np.isclose(terms[perm], val), (perm, terms[perm])
assert np.isclose(sum(terms.values()), 7)
print("Example 4.1 Leibniz terms (1-based permutations):",
      {tuple(p + 1 for p in k): v for k, v in terms.items()})

# %% [markdown]
# ## Part 2: Cofactor cost (n!) against elimination cost (n^3)
#
# §3.5/§3.13: T_cofactor(n) = sum_{k=1}^{n-1} n!/k! ~ (e - 1) n! multiplications,
# T_elimination(n) ~ (2/3) n^3. First the exact multiplication count, then wall time:
# each is timed on one random matrix per size; elimination is repeated to get a
# measurable time.

# %%
for n in range(1, 10):
    _, mults = det_cofactor(rng.standard_normal((n, n)), count=True)
    assert mults == sum(math.factorial(n) // math.factorial(k) for k in range(1, n)), (n, mults)
_, mults9 = det_cofactor(rng.standard_normal((9, 9)), count=True)
assert abs(mults9 / math.factorial(9) - (math.e - 1)) < 1e-5
T20 = sum(math.factorial(20) // math.factorial(k) for k in range(1, 20))
assert 4.1e18 < T20 < 4.25e18                                   # "about 4.2 x 10^18" (§3.13)
print(f"Part 2: cofactor multiplications for n = 9: {mults9} = {mults9 / math.factorial(9):.6f} * 9!"
      f" (e - 1 = {math.e - 1:.6f}); for n = 20: {T20:.2e}")

ns = list(range(2, 10))
t_cof, t_elim = [], []
for n in ns:
    A = rng.standard_normal((n, n))
    best = np.inf
    for _ in range(5 if n <= 8 else 1):                   # best of 5 where it is cheap
        t0 = time.perf_counter()
        det_cofactor(A)
        best = min(best, time.perf_counter() - t0)
    t_cof.append(best)
    reps = 100
    t0 = time.perf_counter()
    for _ in range(reps):
        det_elimination(A)
    t_elim.append((time.perf_counter() - t0) / reps)
assert t_cof[-1] > 50 * t_cof[2], "cofactor time should explode between n = 4 and n = 9"
# Past the small sizes where overhead dominates, each step multiplies cofactor time by ~n
# (geometric mean of the ratios for n = 7, 8, 9 should be near 8), while elimination
# time grows by a small factor over the whole range.
step = (t_cof[-1] / t_cof[-4]) ** (1 / 3)
assert 4 < step < 16, step
assert (t_elim[-1] / t_elim[0]) < 0.01 * (t_cof[-1] / t_cof[0])

fig, ax = plt.subplots(figsize=(6.5, 4.2))
ax.semilogy(ns, t_cof, "o-", color=BLUE, lw=2, ms=6, label="cofactor expansion (3.12)")
ax.semilogy(ns, t_elim, "s-", color=ORANGE, lw=2, ms=6, label="elimination (3.7)")
fact = np.array([math.factorial(n) for n in ns], float)
ax.semilogy(ns, t_cof[-1] * fact / fact[-1], "--", color=MUTED, lw=1.2, label="c · n!  (reference)")
ax.set_xlabel("matrix size n"); ax.set_ylabel("seconds per determinant")
ax.set_title("Cofactor expansion is n!, elimination is n³")
ax.grid(alpha=.25); ax.legend(frameon=False)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(FIG / "cofactor_vs_elimination.png", dpi=120); plt.close(fig)
print(f"Part 2: cofactor n=9 took {t_cof[-1]:.2f} s (time ratio per step for n = 7..9: {step:.1f});"
      f" elimination n=9 took {t_elim[-1]*1e6:.0f} µs ({t_elim[-1] / t_elim[0]:.1f}x its n=2 time).")

# %% [markdown]
# ## Part 3: The identities of §3.2 and §3.4 on random matrices
#
# Row swap flips sign (D2); adding a multiple of a row changes nothing (Thm 3.3);
# scaling a row scales det (D3); equal rows give 0 (Thm 3.2); triangular = product
# of diagonal (Thm 3.6); det(AB) = det A det B (Thm 3.10); det A^T = det A (Thm 3.11);
# det A^{-1} = 1/det A, det(cA) = c^n det A (Cor 3.10'); det P_sigma = sgn(sigma).

# %%
for _ in range(50):
    n = int(rng.integers(2, 7))
    A, B = rng.standard_normal((n, n)), rng.standard_normal((n, n))
    dA = np.linalg.det(A)
    S = A.copy(); S[[0, 1]] = S[[1, 0]]
    assert np.isclose(np.linalg.det(S), -dA)                                  # (D2)
    R = A.copy(); R[1] += 3.7 * R[0]
    assert np.isclose(np.linalg.det(R), dA)                                   # Thm 3.3
    Sc = A.copy(); Sc[0] *= -2.5
    assert np.isclose(np.linalg.det(Sc), -2.5 * dA)                           # (D3)
    E = A.copy(); E[1] = E[0]
    assert np.isclose(np.linalg.det(E), 0, atol=1e-10)                        # Thm 3.2
    assert np.isclose(np.linalg.det(np.triu(A)), np.prod(np.diag(A)))         # Thm 3.6, upper
    assert np.isclose(np.linalg.det(np.tril(A)), np.prod(np.diag(A)))         # Thm 3.6, lower
    assert np.isclose(det_elimination(np.tril(A))[0], np.prod(np.diag(A)))
    assert np.isclose(np.linalg.det(A @ B), dA * np.linalg.det(B))            # Thm 3.10
    assert np.isclose(np.linalg.det(A.T), dA)                                 # Thm 3.11
    assert np.isclose(np.linalg.det(np.linalg.inv(A)), 1 / dA)                # Cor 3.10'(a)
    c = rng.uniform(0.5, 3)
    assert np.isclose(np.linalg.det(c * A), c ** n * dA)                      # Cor 3.10'(d)
    Sinv = np.linalg.inv(B)
    assert np.isclose(np.linalg.det(Sinv @ A @ B), dA)                        # Cor 3.10'(c)
    perm = tuple(rng.permutation(n))
    P = np.eye(n)[list(perm)]                                                 # rows e_{perm(i)}
    assert np.isclose(np.linalg.det(P), sign(perm))                           # det P_sigma = sgn
# Complex matrices (remark at the end of §3.4): same formulas, and det A^dagger = conj(det A).
for _ in range(10):
    Z = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
    assert np.isclose(sum(sign(p) * np.prod([Z[i, p[i]] for i in range(4)])
                          for p in itertools.permutations(range(4))), np.linalg.det(Z))
    assert np.isclose(np.linalg.det(Z.conj()), np.conj(np.linalg.det(Z)))
    assert np.isclose(np.linalg.det(Z.conj().T), np.conj(np.linalg.det(Z)))
print("Part 3: row-operation rules, product/transpose/inverse/scaling identities, det P = sgn,"
      " complex conjugate rules: all pass.")

# %% [markdown]
# ## Part 4: Adjugate (3.14) and Cramer (3.13)
#
# A adj A = det(A) I numerically and symbolically; Cramer on the circuit of
# Example 4.5 (symbolic and numeric) and on the interpolation of Example 4.3.


# %%
def adjugate(A):
    """adj A = C^T with C the cofactor matrix; A @ adj A = det(A) I, eq. (3.14)."""
    A = np.asarray(A, float)
    n = A.shape[0]
    C = np.array([[cofactor(A, i, j) for j in range(n)] for i in range(n)])
    return C.T


def cramer(A, b):
    """x_j = det(A_j) / det(A), A_j = A with column j replaced by b, eq. (3.13)."""
    A = np.asarray(A, float); b = np.asarray(b, float)
    d = det_elimination(A)[0]
    x = np.empty(len(b))
    for j in range(len(b)):
        Aj = A.copy(); Aj[:, j] = b
        x[j] = det_elimination(Aj)[0] / d
    return x


for _ in range(20):
    n = int(rng.integers(2, 6))
    A = rng.standard_normal((n, n)); b = rng.standard_normal(n)
    assert np.allclose(A @ adjugate(A), np.linalg.det(A) * np.eye(n))
    assert np.allclose(adjugate(A) @ A, np.linalg.det(A) * np.eye(n))        # B7(c)
    assert np.allclose(cramer(A, b), np.linalg.solve(A, b))
    assert np.isclose(np.linalg.det(adjugate(A)), np.linalg.det(A) ** (n - 1))  # B7(a)

# Symbolic 3x3: A adj A = det(A) I, exactly.
a = sp.symbols("a11 a12 a13 a21 a22 a23 a31 a32 a33")
As = sp.Matrix(3, 3, a)
assert sp.simplify(As * As.adjugate() - As.det() * sp.eye(3)) == sp.zeros(3, 3)
# Symbolic Leibniz formula for 3x3 equals SymPy's det.
leib3 = sum(sign(p) * sp.prod([As[i, p[i]] for i in range(3)]) for p in itertools.permutations(range(3)))
assert sp.expand(leib3 - As.det()) == 0

# Example 4.5: two-loop circuit, Cramer in symbols and numbers.
R1, R2, R3, V = sp.symbols("R1 R2 R3 V", positive=True)
K = sp.Matrix([[R1 + R2, -R2], [-R2, R2 + R3]])
detK = sp.expand(K.det())
assert detK == sp.expand(R1 * R2 + R1 * R3 + R2 * R3)
I1 = K.copy(); I1[:, 0] = sp.Matrix([V, 0]); I2 = K.copy(); I2[:, 1] = sp.Matrix([V, 0])
assert sp.simplify(I1.det() / detK - V * (R2 + R3) / detK) == 0
assert sp.simplify(I2.det() / detK - V * R2 / detK) == 0
I_num = cramer([[2, -1], [-1, 2]], [3, 0])
assert np.allclose(I_num, [2, 1])
print("Example 4.5: det =", detK, "; currents with R=1,1,1 ohm, V=3 V:", I_num)

# Example 4.3: parabola through (1,2), (2,3), (4,11) via Cramer.
V43 = np.vander([1, 2, 4], increasing=True).astype(float)
y43 = np.array([2, 3, 11], float)
nums = []
for j in range(3):
    Vj = V43.copy(); Vj[:, j] = y43
    nums.append(det_elimination(Vj)[0])
assert np.allclose(nums, [18, -12, 6]) and np.isclose(det_elimination(V43)[0], 6)
c43 = cramer(V43, y43)
assert np.allclose(c43, [3, -2, 1])
assert np.allclose(np.polyval(c43[::-1], [1, 2, 4]), y43)
print("Example 4.3: Cramer numerators", nums, "-> coefficients (c0, c1, c2) =", c43)

# %% [markdown]
# ## Part 5: The Vandermonde determinant (3.19)
#
# det V = prod_{i<j} (x_j - x_i), checked numerically on random nodes and
# symbolically for n = 3 (Example 4.3's factorization); nodes 0..3 give 1!2!3! = 12 (B4).


# %%
def vandermonde_product(x):
    x = np.asarray(x, float)
    return float(np.prod([x[j] - x[i] for i in range(len(x)) for j in range(i + 1, len(x))]))


for n in range(2, 7):
    x = rng.uniform(-2, 2, n)
    assert np.isclose(np.linalg.det(np.vander(x, increasing=True)), vandermonde_product(x), rtol=1e-8)
assert np.isclose(np.linalg.det(V43), 6) and np.isclose(vandermonde_product([1, 2, 4]), 6)
assert np.isclose(np.linalg.det(np.vander([0, 1, 2, 3], increasing=True)), 12)

x1, x2, x3 = sp.symbols("x1 x2 x3")
Vs = sp.Matrix([[1, x1, x1 ** 2], [1, x2, x2 ** 2], [1, x3, x3 ** 2]])
assert sp.expand(Vs.det() - (x2 - x1) * (x3 - x1) * (x3 - x2)) == 0
print("Part 5: Vandermonde formula verified numerically (n = 2..6) and symbolically (n = 3).")

# %% [markdown]
# ## Part 6: The cross product is a determinant (3.15), Theorem 3.15
#
# Cofactor expansion along a symbolic first row; triple product = det; Lagrange's
# identity |u x v|^2 = |u|^2 |v|^2 - (u . v)^2 symbolically; A5's numbers.


# %%
def cross_by_cofactors(u, v):
    """u x v = sum_k (-1)^(1+k) M_{1k} e_k, the expansion of (3.15) along row 1."""
    M = np.array([[0, 0, 0], u, v], float)
    return np.array([(-1) ** k * np.linalg.det(np.delete(M[1:], k, axis=1)) for k in range(3)])


u1, u2, u3, v1, v2, v3, w1, w2, w3 = sp.symbols("u1 u2 u3 v1 v2 v3 w1 w2 w3")
us, vs, ws = sp.Matrix([u1, u2, u3]), sp.Matrix([v1, v2, v3]), sp.Matrix([w1, w2, w3])
cross_s = us.cross(vs)
assert sp.expand(ws.dot(cross_s) - sp.Matrix.vstack(ws.T, us.T, vs.T).det()) == 0      # 3.15(a)
assert sp.expand(cross_s.dot(cross_s) - (us.dot(us) * vs.dot(vs) - us.dot(vs) ** 2)) == 0  # 3.15(c)
for _ in range(20):
    u, v = rng.standard_normal(3), rng.standard_normal(3)
    assert np.allclose(cross_by_cofactors(u, v), np.cross(u, v))
    assert np.isclose(np.dot(u, np.cross(u, v)), 0) and np.isclose(np.dot(v, np.cross(u, v)), 0)
uA5, vA5, wA5 = [1, 2, 0], [0, 1, 3], [1, 1, 1]
assert np.allclose(np.cross(uA5, vA5), [6, -3, 1])
assert np.isclose(np.linalg.norm(np.cross(uA5, vA5)) ** 2, 46)
assert np.isclose(np.dot(wA5, np.cross(uA5, vA5)), 4) and np.isclose(np.linalg.det(np.array([wA5, uA5, vA5], float)), 4)
print("Part 6: cross product via cofactors, triple product = det, Lagrange identity: all pass.")

# %% [markdown]
# ## Part 7: |det A| is the volume scaling factor (Theorem 3.17), by Monte Carlo
#
# Map 100k uniform points in the unit square (cube) through A. Estimate the image's
# area (volume) two ways: the convex hull of the mapped points, and bounding-box
# rejection (sample the bounding box of the image, count points whose preimage
# A^{-1} y lies in [0,1]^n). Both should match |det A|.


# %%
def image_volume_estimates(A, n_pts=100_000):
    """(hull volume of A @ uniform points, rejection estimate, |det A|) for the unit cube."""
    A = np.asarray(A, float)
    n = A.shape[0]
    pts = rng.uniform(0, 1, (n_pts, n))
    img = pts @ A.T                                      # row i is A @ pts[i]
    hull_vol = ConvexHull(img).volume                    # 'volume' is area in 2D
    corners = np.array(list(itertools.product([0, 1], repeat=n))) @ A.T
    lo, hi = corners.min(axis=0), corners.max(axis=0)
    box = rng.uniform(lo, hi, (n_pts, n))
    pre = box @ np.linalg.inv(A).T                       # A^{-1} applied to each sample
    inside = np.all((pre >= 0) & (pre <= 1), axis=1).mean()
    rejection_vol = inside * np.prod(hi - lo)
    return hull_vol, rejection_vol, abs(np.linalg.det(A)), img


M42 = np.array([[np.sqrt(3) / 2, -0.5], [0.5, np.sqrt(3) / 2]]) @ np.array([[1, 2], [0, 1]])   # Example 4.2
A2 = np.array([[2.0, 1.0], [0.5, 1.5]])                                                        # det = 2.5
A3 = np.array([[2, 1, 0], [0, 1, 1], [1, 0, 3]], float)                                       # det = 7


def random_well_conditioned(n):
    """A random n x n matrix with cond < 3, so Monte Carlo estimates are not starved."""
    while True:
        A = rng.standard_normal((n, n))
        if np.linalg.cond(A) < 3:
            return A


rand_mats = [random_well_conditioned(2) for _ in range(3)] + [random_well_conditioned(3) for _ in range(2)]
for A in [M42, A2, A3] + rand_mats:
    hull, rej, dv, img = image_volume_estimates(A)
    assert abs(hull - dv) / dv < 0.01, (hull, dv)
    assert abs(rej - dv) / dv < 0.02, (rej, dv)
    print(f"  n={A.shape[0]}: hull {hull:.4f}, rejection {rej:.4f}, |det A| = {dv:.4f}")
assert np.isclose(np.linalg.det(A3), 7)

# Theorem 3.17(c) for a region that is not a parallelepiped: the unit disk (area pi).
# Its image under A is an ellipse; estimate its area by rejection in a bounding box
# (y is in the image iff ||A^{-1} y|| <= 1) and compare with pi |det A|.
for A in [A2] + rand_mats[:3]:
    R_box = np.linalg.norm(A, 2)                         # the image lies in the disk of radius ||A||_2
    ys = rng.uniform(-R_box, R_box, (200_000, 2))
    inside = np.linalg.norm(ys @ np.linalg.inv(A).T, axis=1) <= 1
    disk_img = inside.mean() * (2 * R_box) ** 2
    assert abs(disk_img - np.pi * abs(np.linalg.det(A))) / (np.pi * abs(np.linalg.det(A))) < 0.02
    print(f"  unit disk -> area {disk_img:.4f}, pi |det A| = {np.pi * abs(np.linalg.det(A)):.4f}")
print("Part 7: Monte Carlo area/volume of images of the unit square, cube and disk matches |det A| * vol.")

hull, rej, dv, img = image_volume_estimates(A2, 20_000)
fig, ax = plt.subplots(figsize=(5.5, 5))
ax.scatter(img[:, 0], img[:, 1], s=2, color=BLUE, alpha=.25, label="A · (uniform points in unit square)")
para = np.array([[0, 0], A2[:, 0], A2[:, 0] + A2[:, 1], A2[:, 1], [0, 0]])
ax.plot(para[:, 0], para[:, 1], color=INK, lw=1.5, label=f"image of unit square, area |det A| = {dv:.2f}")
ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=ORANGE, lw=1.5, ls="--", label="unit square, area 1")
ax.set_aspect("equal"); ax.grid(alpha=.25); ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.set_title(f"Monte Carlo area: hull {hull:.3f}, rejection {rej:.3f}")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(FIG / "volume_monte_carlo.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 8: Example 4.2 drawn: shear, then rotate, then reflect; orientation
#
# Signed area of the image of the unit square = det. Shear (det 1), rotation (det 1),
# reflection (det -1), stretch diag(3, 1/2) (det 3/2). The arrows show where e1 and
# e2 go; the sign of det says whether (A e1, A e2) is still a counterclockwise pair.

# %%
S = np.array([[1, 2], [0, 1]], float)
theta = np.pi / 6
R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
F = np.diag([1.0, -1.0])
D = np.diag([3.0, 0.5])
RS = R @ S
assert np.allclose(RS, M42)
assert np.allclose(RS[:, 0], [0.8660254, 0.5], atol=1e-6) and np.allclose(RS[:, 1], [1.2320508, 1.8660254], atol=1e-6)
assert np.isclose(np.linalg.det(S), 1) and np.isclose(np.linalg.det(R), 1) and np.isclose(np.linalg.det(RS), 1)
assert np.isclose(np.linalg.det(F @ RS), -1) and np.isclose(np.linalg.det(RS @ D), 1.5)
shoelace = RS[0, 0] * RS[1, 1] - RS[0, 1] * RS[1, 0]
assert np.isclose(shoelace, 1.0)
print(f"Example 4.2: M e1 = {RS[:, 0].round(3)}, M e2 = {RS[:, 1].round(3)}, det M = {shoelace:.6f};"
      f" det(F M) = {np.linalg.det(F @ RS):.3f}; det(M D) = {np.linalg.det(RS @ D):.3f}")


def draw_image_of_square(ax, A, title):
    """Draw A applied to the unit square, with A e1 (blue) and A e2 (orange) as arrows."""
    sq = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]], float)
    img = sq @ A.T
    ax.fill(img[:, 0], img[:, 1], color=AQUA, alpha=.25)
    ax.plot(img[:, 0], img[:, 1], color=INK, lw=1.2)
    ax.plot(sq[:, 0], sq[:, 1], color=MUTED, lw=1, ls="--")
    ax.annotate("", xy=A[:, 0], xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2))
    ax.annotate("", xy=A[:, 1], xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2))
    ax.text(*(A[:, 0] * 1.08), "A e₁", color=BLUE, fontsize=9)
    ax.text(*(A[:, 1] * 1.08), "A e₂", color=ORANGE, fontsize=9)
    d = np.linalg.det(A)
    ax.set_title(f"{title}\ndet = {d:+.2f} ({'CCW kept' if d > 0 else 'orientation flipped'})", fontsize=9)
    ax.set_aspect("equal"); ax.grid(alpha=.25)
    ax.set_xlim(-0.5, 3.6); ax.set_ylim(-2.2, 2.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


fig, axes = plt.subplots(1, 5, figsize=(17.5, 4.8))
draw_image_of_square(axes[0], np.eye(2), "unit square")
draw_image_of_square(axes[1], S, "shear S")
draw_image_of_square(axes[2], RS, "rotate after shear: RS")
draw_image_of_square(axes[3], F @ RS, "then reflect: F R S")
draw_image_of_square(axes[4], RS @ D, "stretch D first: R S D")
assert np.isclose(abs(np.linalg.det(RS @ D)), 1.5)        # area 3/2 (Example 4.2)
fig.suptitle("Example 4.2: the determinant is signed area of the image of the unit square (dashed = original)")
fig.tight_layout(); fig.savefig(FIG / "example_4_2_shear_rotate_reflect.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 9: Jacobians (§3.11, Example 4.4, B6) and the quarter annulus
#
# Symbolic Jacobian determinants: polar r, cylindrical r, spherical r^2 sin(phi).
# Then a Monte Carlo estimate of the area of {1 <= r <= 2, 0 <= theta <= pi/2}:
# sampling (r, theta) uniformly in the parameter rectangle and weighting by |det J| = r
# gives 3 pi / 4; forgetting the weight gives pi / 2.

# %%
r, th, ph, z = sp.symbols("r theta phi z", positive=True)
polar = sp.Matrix([r * sp.cos(th), r * sp.sin(th)])
assert sp.simplify(polar.jacobian([r, th]).det()) == r
cyl = sp.Matrix([r * sp.cos(th), r * sp.sin(th), z])
assert sp.simplify(cyl.jacobian([r, th, z]).det()) == r
sph = sp.Matrix([r * sp.sin(ph) * sp.cos(th), r * sp.sin(ph) * sp.sin(th), r * sp.cos(ph)])
assert sp.simplify(sph.jacobian([r, ph, th]).det() - r ** 2 * sp.sin(ph)) == 0
assert sp.simplify(sph.jacobian([r, th, ph]).det() + r ** 2 * sp.sin(ph)) == 0      # column swap flips sign
Jp = polar.jacobian([r, th])
assert sp.simplify(Jp[:, 0].dot(Jp[:, 1])) == 0 and sp.simplify(Jp[:, 0].norm()) == 1 and sp.simplify(Jp[:, 1].norm()) == r
print("Part 9: det J = r (polar), r (cylindrical), r^2 sin(phi) (spherical), symbolically.")

N = 200_000
rs = rng.uniform(1, 2, N); ths = rng.uniform(0, np.pi / 2, N)
rect_area = (2 - 1) * (np.pi / 2)
area_weighted = rect_area * np.mean(rs)              # ∫∫ r dr dθ ≈ (rect area) × mean of r
area_naive = rect_area                               # ∫∫ 1 dr dθ
exact = 3 * np.pi / 4
assert abs(area_weighted - exact) / exact < 0.005, (area_weighted, exact)
assert np.isclose(exact / area_naive, 1.5)           # off by exactly the factor 3/2 = mean of r on [1, 2]
# Independent check: rejection sampling in the plane.
pts = rng.uniform([0, 0], [2, 2], (N, 2))
rad = np.hypot(pts[:, 0], pts[:, 1])
area_plane = 4 * np.mean((rad >= 1) & (rad <= 2))
assert abs(area_plane - exact) / exact < 0.01
print(f"Quarter annulus: exact 3π/4 = {exact:.4f}; with Jacobian weight r: {area_weighted:.4f}; "
      f"without: {area_naive:.4f}; direct rejection in the plane: {area_plane:.4f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.6))
rg = np.linspace(1, 2, 6); tg = np.linspace(0, np.pi / 2, 7)
for rv in rg:
    ax1.axvline(rv, color=MUTED, lw=.8); ax2.plot(rv * np.cos(np.linspace(0, np.pi / 2, 100)), rv * np.sin(np.linspace(0, np.pi / 2, 100)), color=MUTED, lw=.8)
for tv in tg:
    ax1.axhline(tv, color=MUTED, lw=.8); ax2.plot([np.cos(tv), 2 * np.cos(tv)], [np.sin(tv), 2 * np.sin(tv)], color=MUTED, lw=.8)
# Shade one cell near r = 1 and one near r = 2 to show the area ratio (r1^2 - r0^2)/2 * dθ.
cell_ratio = (rg[5] ** 2 - rg[4] ** 2) / (rg[1] ** 2 - rg[0] ** 2)
assert np.isclose(cell_ratio, 0.76 / 0.44)
for (i, col) in ((0, BLUE), (4, ORANGE)):
    r0, r1 = rg[i], rg[i + 1]; t0, t1 = tg[2], tg[3]
    ax1.fill([r0, r1, r1, r0], [t0, t0, t1, t1], color=col, alpha=.5)
    tt = np.linspace(t0, t1, 30)
    xs = np.r_[r0 * np.cos(tt), r1 * np.cos(tt[::-1])]; ys = np.r_[r0 * np.sin(tt), r1 * np.sin(tt[::-1])]
    ax2.fill(xs, ys, color=col, alpha=.5)
ax1.set_xlabel("r"); ax1.set_ylabel("θ"); ax1.set_title("(r, θ) rectangle: all cells have equal area dr dθ")
ax2.set_aspect("equal"); ax2.set_xlabel("x"); ax2.set_ylabel("y")
ax2.set_title(f"image: cell area ≈ r dr dθ (orange/blue = {cell_ratio:.2f})")
for ax in (ax1, ax2):
    ax.grid(alpha=.15)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(FIG / "polar_jacobian.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 10: Scale dependence, or why det does not measure "how singular" (§3.13)
#
# det(2A) = 2^n det A. diag(1e6, 1e-6) has det 1 and condition number 1e12.
# 0.1 * I_100 has det 1e-100 and condition number 1. Use slogdet for large n.

# %%
for n in (2, 5, 10):
    A = rng.standard_normal((n, n))
    assert np.isclose(np.linalg.det(2 * A), 2 ** n * np.linalg.det(A))
bad = np.diag([1e6, 1e-6])
good = 0.1 * np.eye(100)
sgn_g, logdet_g = np.linalg.slogdet(good)
assert np.isclose(np.linalg.det(bad), 1.0) and np.isclose(np.linalg.cond(bad), 1e12)
assert sgn_g == 1 and np.isclose(logdet_g, -100 * np.log(10)) and np.isclose(np.linalg.cond(good), 1.0)
print(f"diag(1e6, 1e-6): det = {np.linalg.det(bad):.1f}, cond = {np.linalg.cond(bad):.1e}  (nearly singular, det says nothing)")
print(f"0.1 * I_100:     det = 10^{logdet_g / np.log(10):.0f}, cond = {np.linalg.cond(good):.1f}  (perfectly conditioned, det is tiny)")
Ar = rng.standard_normal((5, 5))
print(f"random 5x5: det = {np.linalg.det(Ar):+.3f}, det(2A) = {np.linalg.det(2 * Ar):+.3f} = 32 * det(A), cond = {np.linalg.cond(Ar):.1f}")

# %% [markdown]
# ## Part 11: Numbers from the exercises (A1, A3, A4, B1(c), B3(c), C1(a), C2(c))

# %%
assert np.isclose(det_elimination([[2, 1, 0], [4, 3, 1], [-2, 0, 5]])[0], 8)                        # A1
assert [sign(p) for p in ((2, 0, 1), (1, 0, 3, 2), (3, 2, 0, 1), (2, 0, 3, 1))] == [1, 1, -1, -1]   # A3, 0-based
assert sum(sign(p) == 1 for p in itertools.permutations(range(4))) == 12                             # A3
assert np.allclose(cramer([[3, 1], [4, 2]], [5, 8]), [1, 2])                                        # A4
assert np.allclose(adjugate([[3, 1], [4, 2]]), [[2, -1], [-4, 3]])                                   # A4
assert np.isclose(-2 * np.linalg.det(np.array([[3, 0], [0, 2]], float)), -12)                       # B3(c)
M_b1c = np.block([[np.zeros((2, 2)), np.eye(2)], [np.eye(2), np.zeros((2, 2))]])                     # B1(c)
assert np.isclose(np.linalg.det(M_b1c), 1) and sign((2, 3, 0, 1)) == 1                              # det = +1, not -1
t, lam = sp.symbols("t lamda")
Bs = sp.Matrix(3, 3, sp.symbols("b11 b12 b13 b21 b22 b23 b31 b32 b33"))
assert sp.expand(sp.diff((sp.eye(3) + t * Bs).det(), t).subs(t, 0) - Bs.trace()) == 0              # C1(a)
p2 = (sp.Matrix([[2, 1], [1, 2]]) - lam * sp.eye(2)).det()
assert sp.expand(p2 - (lam ** 2 - 4 * lam + 3)) == 0 and sorted(sp.solve(p2, lam)) == [1, 3]        # C2(c)
print("Part 11: exercise numbers A1, A3, A4, B1(c), B3(c), C1(a), C2(c) verified.")

# %% [markdown]
# ## Your turn
#
# 1. **Adjugate inverse (B8).** `adjugate` is in Part 4. Compute the inverse of the
#    4x4 Hilbert matrix as `adjugate(H) / det(H)` and as `np.linalg.inv(H)`, compare
#    both with the exact rational inverse from SymPy, and report the errors and
#    `np.linalg.cond(H)`. Does the formula buy accuracy? What does it cost?
# 2. **Tetrahedron volume (B5).** Write `tetra_volume(P0, P1, P2, P3)` using
#    |det(P1 - P0, P2 - P0, P3 - P0)| / 6. Check it on the vertices of B5(c)
#    (volume 4), then rotate all four vertices by a random rotation and translate
#    them: the volume should not move.
# 3. **Repeated row.** Build a 3x3 matrix whose first two rows are equal and whose
#    third row is random. Perturb the third row 1000 times and record det. It should
#    stay at 0 (to rounding) every time: Theorem 3.2 does not care about row 3.
#    Then perturb row 2 instead and watch det become nonzero.
# 4. **Hadamard and Gram (C3).** Check |det A| <= prod of column norms on 1000
#    random 5x5 matrices, equality for orthogonal columns, and det(A^T A) = |a1 x a2|^2
#    for a random 3x2 matrix.
#
# Reference implementations follow; try yours first.

# %%
# 1. Adjugate inverse vs np.linalg.inv on the 4x4 Hilbert matrix.
H = np.array([[1 / (i + j + 1) for j in range(4)] for i in range(4)])
H_exact_inv = np.array(sp.Matrix(4, 4, lambda i, j: sp.Rational(1, i + j + 1)).inv(), float)
err_adj = np.abs(adjugate(H) / det_elimination(H)[0] - H_exact_inv).max()
err_inv = np.abs(np.linalg.inv(H) - H_exact_inv).max()
kappa = np.linalg.cond(H)
# B8's claim: both errors are ~6e-10 in the worst entry, i.e. a relative error ~1e-13 that sits
# inside the kappa * eps bound conditioning allows, and the two methods agree to within a few
# percent: the adjugate formula buys no accuracy.
assert err_adj < 1e-8 and err_inv < 1e-8            # both fine for a 4x4; neither is exact
assert max(err_adj, err_inv) / np.abs(H_exact_inv).max() < kappa * np.finfo(float).eps
assert 0.5 < err_adj / err_inv < 2.0, (err_adj, err_inv)
print(f"Hilbert 4x4: cond = {kappa:.3g}; max |error| adjugate/det = {err_adj:.2e}, np.linalg.inv = {err_inv:.2e};"
      f" largest exact entry = {np.abs(H_exact_inv).max():.0f}")


# 2. Tetrahedron volume.
def tetra_volume(P0, P1, P2, P3):
    """Volume = |det(P1-P0, P2-P0, P3-P0)| / 6 (B5): the image of the standard simplex, whose volume is 1/6."""
    P0, P1, P2, P3 = (np.asarray(P, float) for P in (P0, P1, P2, P3))
    return abs(np.linalg.det(np.column_stack([P1 - P0, P2 - P0, P3 - P0]))) / 6


verts = [[0, 0, 0], [2, 0, 0], [0, 3, 0], [1, 1, 4]]
assert np.isclose(tetra_volume(*verts), 4.0)
Qrot, _ = np.linalg.qr(rng.standard_normal((3, 3)))
if np.linalg.det(Qrot) < 0:
    Qrot[:, 0] *= -1                                  # make it a proper rotation, det +1
assert np.isclose(np.linalg.det(Qrot), 1.0)
moved = [Qrot @ np.array(P, float) + np.array([5, -2, 1.5]) for P in verts]
assert np.isclose(tetra_volume(*moved), 4.0)
# Same via the 4x4 determinant of B5(d).
M4 = np.vstack([np.ones(4), np.array(verts, float).T])
assert np.isclose(abs(np.linalg.det(M4)) / 6, 4.0)
print("Tetrahedron of B5(c): volume", tetra_volume(*verts), "(unchanged by rotation + translation)")

# 3. Repeated rows survive perturbation of the other row.
base = rng.standard_normal(3)
dets_row3 = []
for _ in range(1000):
    A = np.array([base, base, rng.standard_normal(3)])
    dets_row3.append(np.linalg.det(A))
assert np.max(np.abs(dets_row3)) < 1e-12
A = np.array([base, base + 1e-3 * rng.standard_normal(3), rng.standard_normal(3)])
assert abs(np.linalg.det(A)) > 1e-6
print(f"Rows 1 = 2, row 3 random (1000 trials): max |det| = {np.max(np.abs(dets_row3)):.1e}; "
      f"perturb row 2 by 1e-3 instead: det = {np.linalg.det(A):.2e}")

# 4. Hadamard's inequality and the Gram determinant.
for _ in range(1000):
    A = rng.standard_normal((5, 5))
    assert abs(np.linalg.det(A)) <= np.prod(np.linalg.norm(A, axis=0)) * (1 + 1e-12)
Qo, _ = np.linalg.qr(rng.standard_normal((5, 5)))
A_orth = Qo @ np.diag(rng.uniform(0.5, 2, 5))
assert np.isclose(abs(np.linalg.det(A_orth)), np.prod(np.linalg.norm(A_orth, axis=0)))
A32 = rng.standard_normal((3, 2))
assert np.isclose(np.linalg.det(A32.T @ A32), np.linalg.norm(np.cross(A32[:, 0], A32[:, 1])) ** 2)
print("Hadamard |det A| <= prod ||a_j|| on 1000 matrices, equality for orthogonal columns; Gram det = |a1 x a2|^2.")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
