# %% [markdown]
# # EM-1 Project: The Vector Calculus Toolkit, Checked
#
# Every operator, theorem, and worked-example number of EM-1 is checked here.
# Parts: (1) grad/div/curl, the second-derivative identities and product rules
# in SymPy; (2) curvilinear coordinates: Jacobians, scale factors, and every
# line of the §3.10 table checked against Cartesian; (3) finite-difference div
# and curl on a grid vs the analytic answers; (4) a divergence-theorem checker;
# (5) a Stokes checker; (6) the flux of r̂/r² through spheres; (7) the delta
# function three ways; (8) Helmholtz decomposition by FFT; (9) potentials;
# (10) numbers quoted in the exercise solutions. Reference solutions for B8,
# C1(c), and C4 sit after "Your turn".
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import itertools
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.special import erf
from scipy.integrate import quad

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
x, y, z = sp.symbols("x y z", real=True)
X = (x, y, z)
PI = np.pi

def gl(a, b, n):
    """Gauss–Legendre nodes and weights on [a, b]; exact for polynomials of degree < 2n."""
    t, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * t + 0.5 * (b + a), 0.5 * (b - a) * w

def numfield(F):
    """SymPy 3-vector in (x, y, z) -> NumPy function returning three broadcast arrays."""
    fs = [sp.lambdify(X, c, "numpy") for c in F]

    def Fnum(X_, Y_, Z_):
        shape = np.broadcast(X_, Y_, Z_).shape
        return [np.broadcast_to(np.asarray(fi(X_, Y_, Z_), float), shape) for fi in fs]
    return Fnum

# %% [markdown]
# ## Part 1: The operators and their identities (§3.1, §3.4–3.6, Example 4.1)
#
# (3.5) ∇·F = ∂_i F_i; (3.6) (∇×F)_i = ε_ijk ∂_j F_k. Theorem 3.5: curl grad = 0,
# div curl = 0. Theorem 3.6: ∇×(∇×A) = ∇(∇·A) − ∇²A. The six product rules (P1)–(P6).

# %%
def grad(f):
    """∇f, the vector of (3.4)."""
    return sp.Matrix([sp.diff(f, v) for v in X])

def div(F):
    """∇·F = ∂_i F_i, eq. (3.5)."""
    return sum(sp.diff(F[i], X[i]) for i in range(3))

def curl(F):
    """(∇×F)_i = ε_ijk ∂_j F_k, eq. (3.6)."""
    return sp.Matrix([sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x),
                      sp.diff(F[1], x) - sp.diff(F[0], y)])

lap = lambda f: div(grad(f))                                                   # ∇² = ∇·∇
adv = lambda A, B: sp.Matrix([sum(A[j] * sp.diff(B[i], X[j]) for j in range(3)) for i in range(3)])  # (A·∇)B
zero = lambda e: sp.simplify(e) == 0

Z3 = sp.zeros(3, 1)
rr = sp.sqrt(x**2 + y**2 + z**2)

# Example 4.1: the four fields of §2.
F_r, F_c = sp.Matrix([x, y, z]), sp.Matrix([x, y, z]) / rr**3
F_w, F_v = sp.Matrix([-y, x, 0]), sp.Matrix([-y, x, 0]) / (x**2 + y**2)
assert div(F_r) == 3 and curl(F_r) == Z3
assert zero(div(F_c)) and all(zero(c) for c in curl(F_c))
assert div(F_w) == 0 and curl(F_w) == sp.Matrix([0, 0, 2])
assert zero(div(F_v)) and all(zero(c) for c in curl(F_v))
assert zero(lap(1 / rr))                              # Example 4.5(a)
print("Example 4.1: div r = 3, div(r̂/r²) = 0 and curl = 0 off the origin; "
      "curl(-y,x,0) = (0,0,2); vortex curl = 0 off the axis.")

# Identities with arbitrary smooth functions.
f, g = sp.Function("f")(x, y, z), sp.Function("g")(x, y, z)
A = sp.Matrix([sp.Function(n)(x, y, z) for n in ("A1", "A2", "A3")])
B = sp.Matrix([sp.Function(n)(x, y, z) for n in ("B1", "B2", "B3")])
assert sp.simplify(curl(grad(f))) == Z3 and zero(div(curl(A)))             # (3.11), (3.12)
lapA = sp.Matrix([lap(c) for c in A])
assert sp.simplify(curl(curl(A)) - grad(div(A)) + lapA) == Z3              # (3.13)
P = {
    "P1": grad(f * g) - (f * grad(g) + g * grad(f)),
    "P2": grad(A.dot(B)) - (A.cross(curl(B)) + B.cross(curl(A)) + adv(A, B) + adv(B, A)),
    "P3": sp.Matrix([div(f * A) - (f * div(A) + A.dot(grad(f)))]),
    "P4": sp.Matrix([div(A.cross(B)) - (B.dot(curl(A)) - A.dot(curl(B)))]),
    "P5": curl(f * A) - (f * curl(A) - A.cross(grad(f))),
    "P6": curl(A.cross(B)) - (adv(B, A) - adv(A, B) + A * div(B) - B * div(A)),
}
assert all(zero(e) for expr in P.values() for e in expr)
assert sp.simplify(curl(f * grad(g)) - grad(f).cross(grad(g))) == Z3       # exercise B1
# (3.2) BAC-CAB and (3.1) the ε–δ identity.
a, b, c = (sp.Matrix(sp.symbols(f"{v}1:4")) for v in "abc")
assert sp.expand(a.cross(b.cross(c)) - (b * a.dot(c) - c * a.dot(b))) == Z3
eps = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    eps[p] = np.linalg.det(np.eye(3)[list(p)])
d = np.eye(3)
assert np.allclose(np.einsum("ijk,imn->jkmn", eps, eps),
                   np.einsum("jm,kn->jkmn", d, d) - np.einsum("jn,km->jkmn", d, d))
print("Part 1: (3.1), (3.2), (3.11)-(3.13) and product rules P1-P6 verified symbolically.")

# Figure: the four fields of §2 / Example 4.1 in the xy-plane.
gq = np.linspace(-2, 2, 16)          # even count: no grid point at the origin
QX, QY = np.meshgrid(gq, gq)
fields = {"r (div 3)": (QX, QY), "r̂/r² (div 4πδ³)": (QX / np.hypot(QX, QY)**3, QY / np.hypot(QX, QY)**3),
          "(−y, x) (curl 2ẑ)": (-QY, QX), "(−y, x)/s² (curl 0 off axis)": (-QY / (QX**2 + QY**2), QX / (QX**2 + QY**2))}
fig, axs = plt.subplots(1, 4, figsize=(16, 4))
for ax, (name, (U, V)) in zip(axs, fields.items()):
    M = np.hypot(U, V)
    ax.quiver(QX, QY, U / M, V / M, np.log(M), cmap="viridis")
    ax.set_title(name); ax.set_aspect("equal")
fig.tight_layout(); fig.savefig(FIG / "four_fields.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 2: Curvilinear coordinates (§3.10)
#
# Jacobians r² sinθ and s; scale factors (1, r, r sinθ) and (1, s, 1). Then every
# formula in the §3.10 table is checked against Cartesian at random points.

# %%
r_, th, ph, s_ = sp.symbols("r theta phi s", positive=True)
sph = sp.Matrix([r_ * sp.sin(th) * sp.cos(ph), r_ * sp.sin(th) * sp.sin(ph), r_ * sp.cos(th)])
cyl = sp.Matrix([s_ * sp.cos(ph), s_ * sp.sin(ph), z])
J_sph, J_cyl = sph.jacobian([r_, th, ph]), cyl.jacobian([s_, ph, z])
assert zero(J_sph.det() - r_**2 * sp.sin(th)) and zero(J_cyl.det() - s_)
G_sph, G_cyl = sp.simplify(J_sph.T * J_sph), sp.simplify(J_cyl.T * J_cyl)   # metric = diag(h_i²)
assert G_sph == sp.diag(1, r_**2, r_**2 * sp.sin(th)**2) and G_cyl == sp.diag(1, s_**2, 1)
print("Jacobians: det = r² sinθ (spherical), s (cylindrical); scale factors (1, r, r sinθ), (1, s, 1).")

e_sph = [J_sph[:, 0], J_sph[:, 1] / r_, J_sph[:, 2] / (r_ * sp.sin(th))]   # r̂, θ̂, φ̂
e_cyl = [J_cyl[:, 0], J_cyl[:, 1] / s_, J_cyl[:, 2]]                       # ŝ, φ̂, ẑ
S = sp.sin(th)
sph_ops = {
    "grad": lambda f: [sp.diff(f, r_), sp.diff(f, th) / r_, sp.diff(f, ph) / (r_ * S)],
    "div": lambda F: sp.diff(r_**2 * F[0], r_) / r_**2 + sp.diff(S * F[1], th) / (r_ * S) + sp.diff(F[2], ph) / (r_ * S),
    "curl": lambda F: [(sp.diff(S * F[2], th) - sp.diff(F[1], ph)) / (r_ * S),
                       (sp.diff(F[0], ph) / S - sp.diff(r_ * F[2], r_)) / r_,
                       (sp.diff(r_ * F[1], r_) - sp.diff(F[0], th)) / r_],
    "lap": lambda f: (sp.diff(r_**2 * sp.diff(f, r_), r_) / r_**2 + sp.diff(S * sp.diff(f, th), th) / (r_**2 * S)
                      + sp.diff(f, ph, 2) / (r_**2 * S**2)),
}
cyl_ops = {
    "grad": lambda f: [sp.diff(f, s_), sp.diff(f, ph) / s_, sp.diff(f, z)],
    "div": lambda F: sp.diff(s_ * F[0], s_) / s_ + sp.diff(F[1], ph) / s_ + sp.diff(F[2], z),
    "curl": lambda F: [sp.diff(F[2], ph) / s_ - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], s_),
                       (sp.diff(s_ * F[1], s_) - sp.diff(F[0], ph)) / s_],
    "lap": lambda f: sp.diff(s_ * sp.diff(f, s_), s_) / s_ + sp.diff(f, ph, 2) / s_**2 + sp.diff(f, z, 2),
}
f_test = x**2 * y + z**3 + x * z
F_test = sp.Matrix([x * y + z**2, y * z - x**2, x + y * z**2])

def check_table(ops, emb, evecs, syms, sample):
    """Compare the curvilinear formulas with Cartesian ones, evaluated at random points."""
    sub = dict(zip(X, emb)) if syms[2] != z else {x: emb[0], y: emb[1]}
    to = lambda e: e.subs(sub, simultaneous=True)
    fc, Fc = to(f_test), F_test.subs(sub, simultaneous=True)
    Fcomp = [Fc.dot(e) for e in evecs]
    pairs = [(ops["grad"](fc), [to(grad(f_test)).dot(e) for e in evecs]),
             ([ops["div"](Fcomp)], [to(div(F_test))]),
             (ops["curl"](Fcomp), [to(curl(F_test)).dot(e) for e in evecs]),
             ([ops["lap"](fc)], [to(lap(f_test))])]
    for lhs, rhs in pairs:
        fl, fr = sp.lambdify(syms, lhs, "numpy"), sp.lambdify(syms, rhs, "numpy")
        for p in sample:
            assert np.allclose(np.array(fl(*p), float), np.array(fr(*p), float), rtol=1e-10, atol=1e-10)

pts = np.c_[rng.uniform(0.5, 2, 10), rng.uniform(0.3, 2.8, 10), rng.uniform(0, 2 * PI, 10)]
check_table(sph_ops, sph, e_sph, (r_, th, ph), pts)
pts_c = np.c_[rng.uniform(0.5, 2, 10), rng.uniform(0, 2 * PI, 10), rng.uniform(-1, 1, 10)]
check_table(cyl_ops, cyl, e_cyl, (s_, ph, z), pts_c)
n_ = sp.symbols("n", integer=True)
assert zero(sph_ops["div"]([r_**n_, 0, 0]) - (n_ + 2) * r_**(n_ - 1))       # A4(a)
assert zero(sph_ops["lap"](r_**n_) - n_ * (n_ + 1) * r_**(n_ - 2))           # A4(c), Example 4.5
assert sph_ops["div"]([1 / r_**2, 0, 0]) == 0 and sph_ops["lap"](1 / r_) == 0  # Example 4.4(a)
assert cyl_ops["div"]([s_, 0, 0]) == 2                                        # A4(b)
assert cyl_ops["curl"]([0, s_, 0])[2] == 2 and cyl_ops["curl"]([0, 1 / s_, 0])[2] == 0   # B2
assert sp.simplify(sph_ops["div"]([1 / r_, 0, 0]) - 1 / r_**2) == 0           # B6(a)
mu = sp.symbols("mu", positive=True)
assert zero(sph_ops["lap"](sp.exp(-mu * r_) / r_) - mu**2 * sp.exp(-mu * r_) / r_)   # C3, r ≠ 0
print("Part 2: all eight formulas of the §3.10 table agree with Cartesian at random points.")

# %% [markdown]
# ## Part 3: Finite-difference div and curl vs the analytic answers
#
# Second-order central differences (np.gradient, edge_order=2). Exact for fields
# that are quadratic in each variable; O(h²) otherwise.
#
# **Build it yourself first:** write `fd_div_curl` before reading it.

# %%
def fd_div_curl(Fnum, n, L=1.0):
    """Grid div and curl of F on [-L, L]^3 with spacing h = 2L/(n-1), eqs. (3.5)-(3.6)."""
    g1 = np.linspace(-L, L, n); h = g1[1] - g1[0]
    G = np.meshgrid(g1, g1, g1, indexing="ij")
    Fx, Fy, Fz = Fnum(*G)
    D = lambda Fc, ax: np.gradient(Fc, h, axis=ax, edge_order=2)
    return G, D(Fx, 0) + D(Fy, 1) + D(Fz, 2), np.stack([D(Fz, 1) - D(Fy, 2), D(Fx, 2) - D(Fz, 0), D(Fy, 0) - D(Fx, 1)])

test_fields = {"Example 4.2": sp.Matrix([x**2 * y + z, 3 * y**2 * z, 2 * x * z**2]),
               "Example 4.3": sp.Matrix([z - y, x, x * y]),
               "whirlpool": F_w,
               "nonpolynomial": sp.Matrix([sp.sin(x * y) + z, sp.exp(z) * sp.cos(x), x * y * z])}
fd_err = {}
for name, F in test_fields.items():
    for n in (21, 41):
        G, dv, cu = fd_div_curl(numfield(F), n)
        dv_exact = np.broadcast_to(sp.lambdify(X, div(F), "numpy")(*G), dv.shape)
        cu_exact = np.stack(numfield(curl(F))(*G))
        fd_err[name, n] = max(np.abs(dv - dv_exact).max(), np.abs(cu - cu_exact).max())
    print(f"  FD vs analytic, {name:14s}: max error h=0.1: {fd_err[name, 21]:.1e}, h=0.05: {fd_err[name, 41]:.1e}")
    assert fd_err[name, 41] < 5e-3
for name in ("Example 4.2", "Example 4.3", "whirlpool"):
    assert fd_err[name, 41] < 1e-10                          # quadratic in each variable: exact
order = np.log2(fd_err["nonpolynomial", 21] / fd_err["nonpolynomial", 41])
assert 1.8 < order < 2.2
print(f"Part 3: finite differences match; observed order for the nonpolynomial field = {order:.2f}.")

# %% [markdown]
# ## Part 4: Divergence-theorem checker (Theorems 3.3 and 3.8, Example 4.2, B3)

# %%
def flux_box(Fnum, box, n=12):
    """Outward flux through each of the six faces of an axis-aligned box (surface side of 3.7)."""
    faces = {}
    for ax, (lo, hi) in enumerate(box):
        j, k = [i for i in range(3) if i != ax]
        (u, wu), (v, wv) = gl(*box[j], n), gl(*box[k], n)
        U, V = np.meshgrid(u, v, indexing="ij")
        for sign, val in ((-1, lo), (1, hi)):
            pts3 = [None] * 3
            pts3[ax], pts3[j], pts3[k] = np.full_like(U, val), U, V
            faces[("xyz"[ax], val)] = sign * np.sum(np.outer(wu, wv) * Fnum(*pts3)[ax])
    return faces

def volume_box(fnum, box, n=12):
    """∫ f dτ over a box by tensor Gauss–Legendre (volume side of 3.7)."""
    nodes = [gl(lo, hi, n) for lo, hi in box]
    G = np.meshgrid(*[p for p, _ in nodes], indexing="ij")
    W = np.einsum("i,j,k->ijk", *[w for _, w in nodes])
    return np.sum(W * np.broadcast_to(fnum(*G), W.shape))

F42, unit = test_fields["Example 4.2"], [(0, 1)] * 3
faces = flux_box(numfield(F42), unit)
expected = {("x", 0): -0.5, ("x", 1): 1.0, ("y", 0): 0.0, ("y", 1): 1.5, ("z", 0): 0.0, ("z", 1): 1.0}
assert all(np.isclose(faces[key], val, atol=1e-12) for key, val in expected.items()), faces
vol42 = volume_box(sp.lambdify(X, div(F42), "numpy"), unit)
assert np.isclose(sum(faces.values()), 3) and np.isclose(vol42, 3)
print("Example 4.2: face fluxes", {f"{k[0]}={k[1]}": round(float(v), 12) + 0.0 for k, v in faces.items()},
      f"total = {sum(faces.values()):.12f}, ∫div = {vol42:.12f}")

Fnl = test_fields["nonpolynomial"]
box = [(-0.3, 1.1), (0.2, 0.9), (-1.0, 0.5)]
lhs, rhs = sum(flux_box(numfield(Fnl), box, 20).values()), volume_box(sp.lambdify(X, div(Fnl), "numpy"), box, 20)
assert np.isclose(lhs, rhs, rtol=1e-12)
print(f"  random box, nonpolynomial field: flux = {lhs:.12f}, ∫div = {rhs:.12f}")

# B3: half-ball of radius R for F = (xz, yz, z² + 1). Expected: hemisphere πR⁴ + πR², base −πR².
R, F_B3 = 2.0, sp.Matrix([x * z, y * z, z**2 + 1])
FB3 = numfield(F_B3)
(t1, w1), (p1, w2) = gl(0, PI / 2, 40), gl(0, 2 * PI, 40)
T, Pp = np.meshgrid(t1, p1, indexing="ij")
nrm = np.array([np.sin(T) * np.cos(Pp), np.sin(T) * np.sin(Pp), np.cos(T)])
hemi = np.sum(np.outer(w1, w2) * np.einsum("i...,i...", np.array(FB3(*(R * nrm))), nrm) * R**2 * np.sin(T))
s1, ws = gl(0, R, 40)
base = -np.sum(np.outer(ws, w2) * np.outer(s1, np.ones_like(p1)))          # F_z = 1 on z = 0, normal −ẑ
(rq, wr) = gl(0, R, 30)
RR, TT, PP = np.meshgrid(rq, t1, p1, indexing="ij")
divB3 = sp.lambdify(X, div(F_B3), "numpy")
vol = np.sum(np.einsum("i,j,k->ijk", wr, w1, w2) * divB3(RR * np.sin(TT) * np.cos(PP), RR * np.sin(TT) * np.sin(PP),
                                                          RR * np.cos(TT)) * RR**2 * np.sin(TT))
assert np.isclose(hemi, PI * R**4 + PI * R**2) and np.isclose(base, -PI * R**2)
assert np.isclose(hemi + base, PI * R**4) and np.isclose(vol, PI * R**4)
print(f"B3 (R = 2): hemisphere {hemi:.6f} = 20π, base {base:.6f} = −4π, total {hemi + base:.6f}, ∫div = {vol:.6f} = 16π")

# %% [markdown]
# ## Part 5: Stokes checker (Theorems 3.4 and 3.9, Example 4.3, B4)

# %%
def circulation(Fnum, curve, dcurve, t0, t1, n=200):
    """∮ F·dl = ∫ F(r(t))·r'(t) dt by Gauss–Legendre."""
    t, w = gl(t0, t1, n)
    return np.sum(w * np.einsum("i...,i...", np.array(Fnum(*curve(t))), dcurve(t)))

def polyline(Fnum, corners, n=40):
    """∫ F·dl along straight segments corners[0] → corners[1] → … (each r(t) = a + (b − a)t)."""
    cs = [np.asarray(c, float) for c in corners]
    return sum(circulation(Fnum, lambda t, a=a, b=b: a[:, None] + (b - a)[:, None] * t,
                           lambda t, a=a, b=b: np.broadcast_to((b - a)[:, None], (3, t.size)), 0, 1, n)
               for a, b in zip(cs[:-1], cs[1:]))

def flux_param(Gnum, surf, ru, rv, urange, vrange, n=60):
    """∫ G·da with da = (r_u × r_v) du dv (the recipe of §3.7)."""
    (u, wu), (v, wv) = gl(*urange, n), gl(*vrange, n)
    U, V = np.meshgrid(u, v, indexing="ij")
    da = np.cross(ru(U, V), rv(U, V), axis=0)
    return np.sum(np.outer(wu, wv) * np.einsum("i...,i...", np.array(Gnum(*surf(U, V))), da))

F43 = test_fields["Example 4.3"]
assert curl(F43) == sp.Matrix([x, 1 - y, 2])
C43 = numfield(curl(F43))
hemi43 = flux_param(C43, lambda u, v: R * np.array([np.sin(u) * np.cos(v), np.sin(u) * np.sin(v), np.cos(u)]),
                    lambda u, v: R * np.array([np.cos(u) * np.cos(v), np.cos(u) * np.sin(v), -np.sin(u)]),
                    lambda u, v: R * np.array([-np.sin(u) * np.sin(v), np.sin(u) * np.cos(v), 0 * u]),
                    (0, PI / 2), (0, 2 * PI))
disk43 = flux_param(C43, lambda u, v: np.array([u * np.cos(v), u * np.sin(v), 0 * u]),
                    lambda u, v: np.array([np.cos(v), np.sin(v), 0 * u]),
                    lambda u, v: np.array([-u * np.sin(v), u * np.cos(v), 0 * u]), (0, R), (0, 2 * PI))
loop43 = circulation(numfield(F43), lambda t: R * np.array([np.cos(t), np.sin(t), 0 * t]),
                     lambda t: R * np.array([-np.sin(t), np.cos(t), 0 * t]), 0, 2 * PI)
assert all(np.isclose(val, 2 * PI * R**2) for val in (hemi43, disk43, loop43))
print(f"Example 4.3 (R = 2): hemisphere {hemi43:.10f}, disk {disk43:.10f}, equator {loop43:.10f}; 8π = {8 * PI:.10f}")

def parallelogram_check(F, p, al, be):
    """Circulation around p → p+α → p+α+β → p+β vs ∫(∇×F)·(α×β) du dv over the flat parallelogram."""
    Fn, Cn = numfield(F), numfield(curl(F))
    circ = polyline(Fn, [p, p + al, p + al + be, p + be, p])
    flux = flux_param(Cn, lambda u, v: p[:, None, None] + al[:, None, None] * u + be[:, None, None] * v,
                      lambda u, v: np.broadcast_to(al[:, None, None], (3,) + u.shape),
                      lambda u, v: np.broadcast_to(be[:, None, None], (3,) + u.shape), (0, 1), (0, 1), 30)
    return circ, flux

p0 = np.array([0.2, -0.1, 0.3])
al0, be0 = np.array([0.6, 0.2, -0.1]), np.array([-0.1, 0.5, 0.4])
circ, flux = parallelogram_check(Fnl, p0, al0, be0)
assert np.isclose(circ, flux, rtol=1e-11)
print(f"  tilted parallelogram, nonpolynomial field: circulation {circ:.12f}, flux of curl {flux:.12f}")
# Theorem 3.4: circulation / area -> n̂·(∇×F) at the centre as the loop shrinks.
Cnl = numfield(curl(Fnl))
nhat = np.cross(al0, be0) / np.linalg.norm(np.cross(al0, be0))
errs = []
for ell in (1.0, 0.1, 0.01):
    circ, _ = parallelogram_check(Fnl, p0, ell * al0, ell * be0)
    m = p0 + ell * (al0 + be0) / 2
    errs.append(abs(circ / np.linalg.norm(np.cross(ell * al0, ell * be0)) - nhat @ np.array(Cnl(*m)).ravel()))
ratios = [errs[0] / errs[1], errs[1] / errs[2]]
assert all(50 < q < 200 for q in ratios) and errs[2] < 1e-5, (errs, ratios)   # §3.4: error falls by ~100 per factor-10 shrink
print("  Theorem 3.4 limit, |circ/area − n̂·curl| for ℓ = 1, 0.1, 0.01:", [f"{e:.1e}" for e in errs])

# B4: triangle (1,0,0) → (0,1,0) → (0,0,1), F = (y, z, x).
F_B4 = sp.Matrix([y, z, x])
V3 = [np.array(v, float) for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0))]
circB4 = polyline(numfield(F_B4), V3, 10)
fluxB4 = np.array(curl(F_B4), float).ravel() @ np.cross(V3[1] - V3[0], V3[2] - V3[0]) / 2
assert np.isclose(circB4, -1.5) and np.isclose(fluxB4, -1.5)
print(f"B4: circulation {circB4:.6f}, flux of curl through triangle {fluxB4:.6f}")

# %% [markdown]
# ## Part 6: The flux of r̂/r² through spheres (§3.11, Example 4.4, C1)
#
# Concentric spheres: 4π for every radius. Off-centre spheres that enclose the
# origin: 4π. Spheres that do not: 0.

# %%
def coulomb(X_, Y_, Z_):
    R3 = (X_**2 + Y_**2 + Z_**2) ** 1.5
    return X_ / R3, Y_ / R3, Z_ / R3

def flux_sphere(Fnum, center, radius, n=80):
    """∮ F·da over a sphere, da = radius² sinθ dθ dφ n̂ (Gauss–Legendre in θ and φ)."""
    (t, wt), (p, wp) = gl(0, PI, n), gl(0, 2 * PI, 2 * n)
    T, Pp = np.meshgrid(t, p, indexing="ij")
    nv = np.array([np.sin(T) * np.cos(Pp), np.sin(T) * np.sin(Pp), np.cos(T)])
    Fv = np.array(Fnum(*(np.asarray(center, float)[:, None, None] + radius * nv)))
    return np.sum(np.outer(wt, wp) * np.einsum("i...,i...", Fv, nv) * radius**2 * np.sin(T))

assert all(np.isclose(flux_sphere(coulomb, (0, 0, 0), Rs), 4 * PI, rtol=1e-12) for Rs in (0.01, 1.0, 10.0, 1000.0))
cases = {"centre (0,0,0.5), radius 1": ((0, 0, 0.5), 1.0, 4 * PI), "centre (0,0,3), radius 1": ((0, 0, 3), 1.0, 0.0),
         "centre (2,2,0), radius 2": ((2, 2, 0), 2.0, 0.0)}
for name, (cen, rad, want) in cases.items():
    val = flux_sphere(coulomb, cen, rad, 120)
    assert abs(val - want) < 5e-9, (name, val)                 # Example 4.4: 4π or 0 to eight decimals
    print(f"  flux of r̂/r² through sphere {name}: {val:.8f}")
print(f"Example 4.4: concentric spheres R = 0.01 … 1000 all give 4π = {4 * PI:.10f}")

offsets = np.r_[np.linspace(0, 0.9, 10), np.linspace(1.1, 2, 10)]
fl = [flux_sphere(coulomb, (0, 0, d_), 1.0, 160) for d_ in offsets]

def disk_solid_angle(a, z0, n=80):
    """Flux of r̂/r² through the disk s ≤ a at height z0 (normal +ẑ) = solid angle (exercise C1)."""
    (s, ws), (p, wp) = gl(0, a, n), gl(0, 2 * PI, 8)
    return np.sum(np.outer(ws, wp) * (z0 / (s**2 + z0**2) ** 1.5 * s)[:, None])

omega = disk_solid_angle(1.0, 1.0)
assert np.isclose(omega, 2 * PI * (1 - 1 / np.sqrt(2)))
assert np.isclose(disk_solid_angle(1e4, 1.0, 400), 2 * PI, rtol=1e-3)
print(f"C1: solid angle of a unit disk seen from height 1 = {omega:.6f} = 2π(1 − 1/√2)")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
gx = np.linspace(-1.55, 4.05, 22); gz = gx
GX, GZ = np.meshgrid(gx, gz)
Ex, _, Ez = coulomb(GX, 0 * GX, GZ); Mg = np.hypot(Ex, Ez)
ax1.quiver(GX, GZ, Ex / Mg, Ez / Mg, np.log(Mg), cmap="viridis")
for (cen, rad, _), col in zip([((0, 0, 0), 1.0, 0)] + list(cases.values()), "rbgm"):
    ax1.add_patch(plt.Circle((cen[0], cen[2]), rad, fill=False, color=col, lw=2))
ax1.set_aspect("equal"); ax1.set_title("r̂/r² in the xz-plane and four test spheres"); ax1.set_xlabel("x"); ax1.set_ylabel("z")
ax2.plot(offsets, np.array(fl) / PI, "o"); ax2.axvline(1, color="k", ls=":")
ax2.set_xlabel("distance of unit-sphere centre from origin"); ax2.set_ylabel("flux / π"); ax2.set_title("4π if the origin is inside, 0 if not")
fig.tight_layout(); fig.savefig(FIG / "coulomb_flux.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 7: The delta function three ways (§3.11, Theorem 3.14, Example 4.5)
#
# (a) Gaussian cloud: u_σ = erf(r/√2σ)/r has ∇²u_σ = −4π g_σ exactly.
# (b) Plummer: ∇²(r² + ε²)^(−1/2) = −3ε²/(r² + ε²)^(5/2), total −4π.
# (c) A grid: the 7-point Laplacian of 1/r summed over a cube is ≈ −4π.

# %%
rs_, sg = sp.symbols("r sigma", positive=True)
u_sig = sp.erf(rs_ / (sp.sqrt(2) * sg)) / rs_
g_sig = (2 * sp.pi * sg**2) ** sp.Rational(-3, 2) * sp.exp(-rs_**2 / (2 * sg**2))
lap_rad = lambda e: sp.diff(rs_**2 * sp.diff(e, rs_), rs_) / rs_**2
assert zero(lap_rad(u_sig) + 4 * sp.pi * g_sig)                       # Example 4.5(b), exact

def pair_u_with_test(u):
    """⟨∇²u, ψ⟩ = ∫ u ∇²ψ dτ for the test function ψ = e^{-r²}, ∇²ψ = (4r² − 6)e^{-r²}."""
    return 4 * PI * quad(lambda r: r**2 * u(r) * (4 * r**2 - 6) * np.exp(-r**2), 0, np.inf, limit=200)[0]

print("  σ      ⟨∇²u_σ, e^{-r²}⟩     −4π(1+2σ²)^(-3/2)")
for s in (0.5, 0.2, 0.1, 0.05):
    val = pair_u_with_test(lambda r, s=s: erf(r / (np.sqrt(2) * s)) / r)
    assert np.isclose(val, -4 * PI * (1 + 2 * s**2) ** -1.5, rtol=1e-8)
    print(f"  {s:<5}  {val:.6f}          {-4 * PI * (1 + 2 * s**2) ** -1.5:.6f}")
val0 = pair_u_with_test(lambda r: 1 / r)
assert np.isclose(val0, -4 * PI)
print(f"  σ → 0 (u = 1/r): ⟨∇²(1/r), ψ⟩ = {val0:.10f} = −4πψ(0) = {-4 * PI:.10f}")
print(f"Example 4.5: σ = 0.1 gives {-4 * PI * 1.02 ** -1.5:.4f}")

eps_ = sp.symbols("epsilon", positive=True)
plummer = lap_rad(1 / sp.sqrt(rs_**2 + eps_**2))
assert zero(plummer + 3 * eps_**2 / (rs_**2 + eps_**2) ** sp.Rational(5, 2))
tot = 4 * PI * quad(lambda r: r**2 * (-3 * 0.3**2 / (r**2 + 0.3**2) ** 2.5), 0, np.inf)[0]
assert np.isclose(tot, -4 * PI)
print(f"  Plummer ε = 0.3: ∫∇²(r²+ε²)^(-1/2) dτ = {tot:.10f}")

# (c) Discrete: grid points at (i + 1/2)h so none sits on the origin.
Nh = 30; hgrid = 1.0 / Nh
g1 = (np.arange(-Nh, Nh) + 0.5) * hgrid
GX, GY, GZ = np.meshgrid(g1, g1, g1, indexing="ij")
fgrid = 1 / np.sqrt(GX**2 + GY**2 + GZ**2)
L7 = (sum(np.roll(fgrid, sft, ax) for ax in range(3) for sft in (1, -1)) - 6 * fgrid)[1:-1, 1:-1, 1:-1] / hgrid**2
rin = np.sqrt(GX**2 + GY**2 + GZ**2)[1:-1, 1:-1, 1:-1]
total = L7.sum() * hgrid**3
near = L7[rin < 4 * hgrid].sum() * hgrid**3
assert abs(total - (-12.568)) < 5e-4 and near / total > 0.998, (total, near / total)   # Example 4.5: −12.568, 99.9 % within 4h
print(f"  7-point Laplacian of 1/r summed over the cube: {total:.5f} (−4π = {-4 * PI:.5f}); "
      f"fraction from r < 4h: {near / total:.4f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
rp = np.linspace(1e-3, 2, 400)
for s in (0.4, 0.2, 0.1):
    ax1.plot(rp, erf(rp / (np.sqrt(2) * s)) / rp, label=f"u_σ, σ = {s}")
    ax2.plot(rp, (2 * PI * s**2) ** -1.5 * np.exp(-rp**2 / (2 * s**2)) * 4 * PI * rp**2, label=f"σ = {s}")
ax1.plot(rp, 1 / rp, "k--", label="1/r"); ax1.set_ylim(0, 10); ax1.legend(); ax1.set_xlabel("r")
ax1.set_title("Smoothed potentials approach 1/r")
ax2.set_xlabel("r"); ax2.set_title("−∇²u_σ · r² (area 4π for every σ)"); ax2.legend()
fig.tight_layout(); fig.savefig(FIG / "delta_regularized.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 8: Helmholtz decomposition by FFT (Theorem 3.15, eq. 3.35, Example 4.7)
#
# On a periodic box, ∇ → ik. The gradient part is k̂(k̂·F̂), the curl part is
# F̂ − k̂(k̂·F̂), the k = 0 mean is a third piece. Potentials: Û = i k·F̂/k²,
# Ŵ = i k×F̂/k² (this W is divergence-free: the Coulomb gauge).

# %%
n = 32
g1 = np.arange(n) * 2 * PI / n
GX, GY, GZ = np.meshgrid(g1, g1, g1, indexing="ij")
k1 = np.fft.fftfreq(n, d=1.0 / n)
K = np.array(np.meshgrid(k1, k1, k1, indexing="ij"))
K2 = np.sum(K**2, axis=0); K2s = np.where(K2 == 0, 1, K2)
fftv = lambda F: np.fft.fftn(F, axes=(1, 2, 3))
ifftv = lambda F: np.real(np.fft.ifftn(F, axes=(1, 2, 3)))

def helmholtz_fft(F):
    """Return gradient part, curl part, mean, U, W of a periodic field F (shape 3×n×n×n)."""
    Fh = fftv(F)
    kF = np.sum(K * Fh, axis=0)
    mean_h = np.zeros_like(Fh); mean_h[:, 0, 0, 0] = Fh[:, 0, 0, 0]
    par_h = K * kF / K2s
    perp_h = Fh - par_h - mean_h
    U = np.real(np.fft.ifftn(1j * kF / K2s))
    W = ifftv(1j * np.cross(K, Fh, axis=0) / K2s)
    return ifftv(par_h), ifftv(perp_h), ifftv(mean_h), U, W

sdiv = lambda F: np.real(np.fft.ifftn(1j * np.sum(K * fftv(F), axis=0)))
scurl = lambda F: ifftv(1j * np.cross(K, fftv(F), axis=0))
sgrad = lambda f: ifftv(1j * K * np.fft.fftn(f)[None])

# Example 4.7: F = (sin(x+y), 0, 0).
S1 = np.sin(GX + GY)
F47 = np.array([S1, 0 * S1, 0 * S1])
par, perp, mean, U, W = helmholtz_fft(F47)
assert np.allclose(par, [S1 / 2, S1 / 2, 0 * S1]) and np.allclose(perp, [S1 / 2, -S1 / 2, 0 * S1])
assert np.allclose(U, np.cos(GX + GY) / 2) and np.allclose(W, [0 * S1, 0 * S1, -np.cos(GX + GY) / 2])
assert np.isclose(np.mean(np.sum(par * perp, axis=0)), 0)
print("Example 4.7: F = (sin(x+y),0,0) splits as ½(sin, sin, 0) + ½(sin, −sin, 0); U = ½cos(x+y), W = −½cos(x+y) ẑ")

# A random smooth field plus a uniform piece.
noise = rng.standard_normal((3, n, n, n))
Fr_ = ifftv(fftv(noise) * np.exp(-K2 / 8)[None]) + np.array([0.3, 0.0, -0.1])[:, None, None, None]
par, perp, mean, U, W = helmholtz_fft(Fr_)
checks = {"curl(grad part)": np.abs(scurl(par)).max(), "div(curl part)": np.abs(sdiv(perp)).max(),
          "reassembly": np.abs(par + perp + mean - Fr_).max(), "−∇U − grad part": np.abs(-sgrad(U) - par).max(),
          "∇×W − curl part": np.abs(scurl(W) - perp).max(), "div W": np.abs(sdiv(W)).max(),
          "⟨grad part, curl part⟩": abs(np.mean(np.sum(par * perp, axis=0)))}
assert all(val < 1e-13 for val in checks.values()), checks            # Example 4.7: every property to better than 1e-13
print("Part 8 (random field): max residuals", {k: f"{v:.1e}" for k, v in checks.items()})
print(f"  mean (k = 0) piece = {mean[:, 0, 0, 0].round(12)}; it is both div- and curl-free")

fig, axs = plt.subplots(1, 3, figsize=(14, 4.5))
sl = (slice(None, None, 2), slice(None, None, 2), 0)
for ax, (name, Fv) in zip(axs, (("F", Fr_), ("gradient part (curl-free)", par), ("curl part (div-free)", perp))):
    ax.quiver(GX[sl], GY[sl], Fv[0][sl], Fv[1][sl]); ax.set_title(name + ", slice z = 0"); ax.set_aspect("equal")
fig.tight_layout(); fig.savefig(FIG / "helmholtz_fft.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 9: Potentials (Theorems 3.16–3.17, Example 4.6, B7)

# %%
F46 = sp.Matrix([2 * x * y - z, x**2 + z**2, 2 * y * z - x])
U46 = -x**2 * y - y * z**2 + x * z
assert curl(F46) == Z3 and sp.expand(-grad(U46) - F46) == Z3
F46n = numfield(F46)
paths = {"straight line": (lambda t: np.array([t, 2 * t, 3 * t]), lambda t: np.array([1 + 0 * t, 2 + 0 * t, 3 + 0 * t])),
         "curve (t², 2 sin(πt/2), 3t³)": (lambda t: np.array([t**2, 2 * np.sin(PI * t / 2), 3 * t**3]),
                                           lambda t: np.array([2 * t, PI * np.cos(PI * t / 2), 9 * t**2]))}
work = {"axis path": polyline(F46n, [(0, 0, 0), (1, 0, 0), (1, 2, 0), (1, 2, 3)], 10)}
work.update({name: circulation(F46n, *pc, 0, 1, 60) for name, pc in paths.items()})
assert all(np.isclose(wv, 17) for wv in work.values()), work
print("Example 4.6: U(1,2,3) =", U46.subs({x: 1, y: 2, z: 3}), "; work along three paths:", {k: round(float(v), 10) for k, v in work.items()})

# B7: F = (y², z², x²), W from (3.36), then the Coulomb-gauge W' = W + ∇(y³z/3).
F_B7 = sp.Matrix([y**2, z**2, x**2])
zp, xp = sp.symbols("zp xp")
W_B7 = sp.Matrix([sp.integrate(F_B7[1].subs(z, zp), (zp, 0, z)),
                  sp.integrate(F_B7[2].subs({x: xp, z: 0}), (xp, 0, x)) - sp.integrate(F_B7[0].subs(z, zp), (zp, 0, z)), 0])
assert div(F_B7) == 0 and sp.expand(curl(W_B7) - F_B7) == Z3
assert sp.expand(W_B7 - sp.Matrix([z**3 / 3, x**3 / 3 - y**2 * z, 0])) == Z3 and div(W_B7) == -2 * y * z
Wp = W_B7 + grad(y**3 * z / 3)
assert sp.expand(Wp - sp.Matrix([z**3, x**3, y**3]) / 3) == Z3 and div(Wp) == 0 and sp.expand(curl(Wp) - F_B7) == Z3
print("B7: W =", list(W_B7), "; Coulomb gauge W' =", list(Wp))

# %% [markdown]
# ## Part 10: Numbers quoted in the exercise solutions

# %%
fA1, FA1 = x**2 * y + y * z**3, sp.Matrix([x * y, y * z, z * x])
pt = {x: 1, y: 1, z: 1}
assert list(grad(fA1).subs(pt)) == [2, 2, 3] and lap(fA1).subs(pt) == 8
assert div(FA1).subs(pt) == 3 and list(curl(FA1)) == [-y, -z, -x]
assert curl(grad(fA1)) == Z3 and div(curl(FA1)) == 0
Av, Bv, Cv = np.array([1, 2, 0]), np.array([0, 1, 3]), np.array([2, 0, 1])
assert Av @ np.cross(Bv, Cv) == 13 and Bv @ np.cross(Cv, Av) == 13
assert np.array_equal(np.cross(Bv, Cv), [1, 6, -2]) and np.array_equal(np.cross(Av, Bv), [6, -3, 1])
assert np.array_equal(np.cross(Cv, Av), [-2, 1, 4]) and (Av @ Cv, Av @ Bv, Bv @ Cv) == (2, 2, 3)
assert np.array_equal(np.cross(Av, np.cross(Bv, Cv)), [-4, 2, 4]) and np.array_equal(np.cross(np.cross(Av, Bv), Cv), [-3, -4, 6])
T_A3 = x**2 + 2 * y**2 + 3 * z**2 + x * y
gT = np.array(grad(T_A3).subs({x: 1, y: -1, z: 2}), float).ravel()
assert np.array_equal(gT, [1, -3, 12])
print(f"A1: ∇f(1,1,1) = (2,2,3), ∇²f = 8, ∇·F = 3, ∇×F = (−1,−1,−1).  A2: volume 13.  "
      f"A3: |∇T| = √154 = {np.linalg.norm(gT):.4f}, D_u T = 10/√3 = {gT.sum() / np.sqrt(3):.4f}")
print("A5: (a) 1 + 4 + 1 = 6, (b) 0, (c) −8π =", f"{-8 * PI:.4f}")
# A5(c) through a Gaussian regularization: ∫(r²+1) g_σ(r − b) dτ = |b|² + 3σ² + 1.
bvec = np.array([0, 0, 1.0])
for s in (0.2, 0.05):
    nodes = [gl(bi - 7 * s, bi + 7 * s, 30) for bi in bvec]
    G = np.meshgrid(*[p for p, _ in nodes], indexing="ij")
    Wt = np.einsum("i,j,k->ijk", *[w for _, w in nodes])
    gs = (2 * PI * s**2) ** -1.5 * np.exp(-((G[0] - bvec[0])**2 + (G[1] - bvec[1])**2 + (G[2] - bvec[2])**2) / (2 * s**2))
    val = -4 * PI * np.sum(Wt * (G[0]**2 + G[1]**2 + G[2]**2 + 1) * gs)
    assert np.isclose(val, -4 * PI * (2 + 3 * s**2), rtol=1e-8)
    print(f"  A5(c) regularized, σ = {s}: {val:.6f} → −8π as σ → 0")
t_ = sp.symbols("t", real=True)
helix = sp.Matrix([sp.cos(t_), sp.sin(t_), t_ / (2 * sp.pi)])
F_A6 = sp.Matrix([y, -x, z])
I_A6 = sp.integrate(F_A6.subs({x: helix[0], y: helix[1], z: helix[2]}).dot(helix.diff(t_)), (t_, 0, 2 * sp.pi))
assert sp.simplify(I_A6 - (sp.Rational(1, 2) - 2 * sp.pi)) == 0 and curl(F_A6) == sp.Matrix([0, 0, -2])
seg_A6 = polyline(numfield(F_A6), [(1, 0, 0), (1, 0, 1)], 10)
assert np.isclose(seg_A6, 0.5) and np.isclose(float(I_A6) - seg_A6, -2 * PI)
proj_area = circulation(numfield(sp.Matrix([-y, x, 0]) / 2), lambda t: np.array([np.cos(t), np.sin(t), t / (2 * PI)]),
                        lambda t: np.array([-np.sin(t), np.cos(t), 0 * t + 1 / (2 * PI)]), 0, 2 * PI)
assert np.isclose(proj_area, PI)                      # ẑ-component of the vector area spanned by helix + segment
print(f"A6: helix {float(I_A6):.4f} = 1/2 − 2π, segment {seg_A6:.4f}, difference {float(I_A6) - seg_A6:.4f} = −2π")
Fv_n = numfield(F_v)
ring = lambda c, a: (lambda t: np.array([c + a * np.cos(t), a * np.sin(t), 0 * t]),
                     lambda t: np.array([-a * np.sin(t), a * np.cos(t), 0 * t]))
circ_B2 = [circulation(Fv_n, *ring(0, a), 0, 2 * PI) for a in (0.5, 3.0)] + [circulation(Fv_n, *ring(2, 1), 0, 2 * PI)]
assert np.allclose(circ_B2, [2 * PI, 2 * PI, 0], atol=1e-10)
print(f"B2: vortex circulation around the axis (a = 0.5, 3): {circ_B2[0]:.6f}, {circ_B2[1]:.6f}; "
      f"around the circle centred (2,0,0): {circ_B2[2]:.1e}")
edges_B4 = [polyline(numfield(F_B4), V3[i:i + 2], 10) for i in range(3)]
assert np.allclose(edges_B4, -0.5)
big = [disk_solid_angle(a, 1.0, 400) for a in (10, 100, 1000)]
assert np.allclose(big, [2 * PI * (1 - 1 / np.sqrt(1 + a * a)) for a in (10, 100, 1000)], rtol=1e-10)
print(f"B4: edge integrals {np.round(edges_B4, 12)}; C1(b): z0 = 1, a = 10, 100, 1000: "
      f"{np.round(big, 6)} → 2π = {2 * PI:.6f}")
yuk = 4 * PI * quad(lambda r: r * np.exp(-1.5 * r) * (4 * r**2 - 6 - 1.5**2) * np.exp(-r**2), 0, np.inf)[0]
assert np.isclose(yuk, -4 * PI)
print(f"C3: ⟨(∇² − μ²)(e^(-μr)/r), e^(-r²)⟩ with μ = 1.5: {yuk:.10f} = −4π")

# %% [markdown]
# ## Your turn
#
# 1. (Exercise B8.) Product rule (P3) on a grid: with `fd_div_curl`'s grid, take
#    f = e^{-(x²+y²+z²)} and A = (xy, yz, zx). Compute ∇·(fA) by finite
#    differences directly, and f∇·A + A·∇f from separately differenced pieces.
#    Assert they agree to O(h²); halve h and confirm the error drops by ≈ 4.
# 2. (Exercise C1.) Write `square_solid_angle(z0)`: the flux of r̂/r² through the
#    square [−½, ½]² in the plane z = z0 (2D Gauss–Legendre in x and y). For
#    z0 = ½ this is one face of the unit cube seen from its centre: check that
#    you get 4π/6 = 2π/3, so the six faces give 4π.
# 3. (Exercise B8.) For F = (2xy − z, x² + z², 2yz − x) add a fourth path of your
#    own to `paths` and confirm the work is still 17. Then do the same with the
#    vortex (−y, x, 0)/(x² + y²) from (1, 0, 0) to (−1, 0, 0) along the upper and
#    lower half circles: you get +π and −π. Which hypothesis of Theorem 3.16 fails?
# 4. In Part 8, replace the random field by (cos z, sin z, sin x) and predict the
#    split by hand before running it.
# 5. (Exercise C4.) Check that the curl of F' (x') = R F(Rᵀx') equals det(R) R curl F
#    for a rotation and for a reflection.
#
# Reference solutions for B8, C1(c), and C4 are in the next cell. Write your own
# first, then delete that cell and rerun.

# %% [markdown]
# ## Reference solutions: B8, C1(c), C4

# %%
# B8(a): (P3) on a grid. Difference the product, or apply the product rule to differenced pieces.
f_b8, A_b8 = sp.exp(-(x**2 + y**2 + z**2)), sp.Matrix([x * y, y * z, z * x])
p3_gap = {}
for n in (21, 41):
    G, dv_fA, _ = fd_div_curl(numfield(f_b8 * A_b8), n)
    _, dv_A, _ = fd_div_curl(numfield(A_b8), n)
    h = G[0][1, 0, 0] - G[0][0, 0, 0]
    fg = np.broadcast_to(sp.lambdify(X, f_b8, "numpy")(*G), dv_A.shape)
    rule = fg * dv_A + sum(Ac * np.gradient(fg, h, axis=ax, edge_order=2) for ax, Ac in enumerate(numfield(A_b8)(*G)))
    p3_gap[n] = np.abs(dv_fA - rule).max()
assert 3.5 < p3_gap[21] / p3_gap[41] < 4.5 and p3_gap[41] < 1e-2
print(f"B8(a): max |FD(∇·(fA)) − (f FD(∇·A) + A·FD(∇f))|: n=21 {p3_gap[21]:.2e}, n=41 {p3_gap[41]:.2e}, "
      f"ratio {p3_gap[21] / p3_gap[41]:.2f}")
# B8(b): a fourth path for Example 4.6, and the vortex along two half circles.
w4 = circulation(F46n, lambda t: np.array([np.sin(PI * t / 2), 2 * t, 3 * t**2]),
                 lambda t: np.array([PI / 2 * np.cos(PI * t / 2), 2 + 0 * t, 6 * t]), 0, 1, 60)
up, down = (circulation(Fv_n, *ring(0, 1), 0, end) for end in (PI, -PI))
assert np.isclose(w4, 17) and np.isclose(up, PI) and np.isclose(down, -PI)
print(f"B8(b): path (sin(πt/2), 2t, 3t²) gives {w4:.10f}; vortex upper half {up:.6f}, lower half {down:.6f}")

# C1(c): one face of the unit cube seen from its centre subtends 4π/6.
def square_solid_angle(z0, n=40):
    """Flux of r̂/r² through the square [-1/2, 1/2]² at height z0 (normal +ẑ)."""
    (u, wu), (v, wv) = gl(-0.5, 0.5, n), gl(-0.5, 0.5, n)
    U_, V_ = np.meshgrid(u, v, indexing="ij")
    return np.sum(np.outer(wu, wv) * z0 / (U_**2 + V_**2 + z0**2) ** 1.5)

assert np.isclose(square_solid_angle(0.5), 2 * PI / 3, rtol=1e-10)
print(f"C1(c): square face from the cube centre: {square_solid_angle(0.5):.10f}; 2π/3 = {2 * PI / 3:.10f}")

# C4: under x' = R x, F'(x') = R F(x): div is a scalar, curl picks up det R.
c45 = (sp.Rational(3, 5), sp.Rational(4, 5))
for Rm in (sp.Matrix([[c45[0], -c45[1], 0], [c45[1], c45[0], 0], [0, 0, 1]]), sp.diag(1, 1, -1)):
    old = dict(zip(X, Rm.T * sp.Matrix(X)))
    for F in (F_test, F_w):
        Fp = Rm * F.subs(old, simultaneous=True)
        assert sp.expand(div(Fp) - div(F).subs(old, simultaneous=True)) == 0
        assert sp.expand(curl(Fp) - Rm.det() * Rm * curl(F).subs(old, simultaneous=True)) == Z3
Fw_mirror = sp.diag(1, 1, -1) * F_w.subs(z, -z)
assert curl(Fw_mirror) == sp.Matrix([0, 0, 2]) and sp.diag(1, 1, -1) * curl(F_w) == sp.Matrix([0, 0, -2])
print("C4: div invariant and curl' = det(R) R curl under a rotation and the reflection z → −z; "
      "mirrored whirlpool still has curl (0, 0, 2)")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
