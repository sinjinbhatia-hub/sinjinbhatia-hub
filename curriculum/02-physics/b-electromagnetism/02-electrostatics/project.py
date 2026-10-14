# %% [markdown]
# # EM-2 Project: Electrostatics on a Computer
#
# This script builds electrostatics from Coulomb's law and checks the notes
# against it: fields and potentials of point charges (field lines of a dipole
# and a quadrupole), a numerical Gauss's-law checker on spheres and cubes, the
# uniformly charged sphere (field, potential, and energy two ways), the line
# charge's reference-point problem, a finite-difference Poisson solver for a
# parallel-plate capacitor with fringing fields, a relaxation solution for a
# grounded conductor in a uniform field, and the charge in an off-centre cavity.
# Every number quoted in the notes' worked examples and in solutions.md is
# asserted or printed here.
#
# SI units throughout. Run headless: `MPLBACKEND=Agg python3 project.py`.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sympy as sp
from scipy import integrate, sparse
from scipy.sparse.linalg import spsolve
from scipy.special import roots_legendre
from scipy.spatial.distance import pdist
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize_scalar

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.perf_counter()
EPS0 = 8.8541878128e-12            # F/m
K = 1 / (4 * np.pi * EPS0)          # Coulomb's constant, 8.98755e9 N m^2/C^2
nC = 1e-9

# %% [markdown]
# ## Part 0: Symbolic checks of the derivations (§3, §4)
#
# SymPy confirms: the Gauss's-law sphere potential (Example 4.1) solves Poisson's
# equation (3.9) inside and Laplace's outside with V and dV/dr continuous; the
# sphere's energy is 3/5 kQ^2/R by (3.15), by (3.16), and shell by shell (C1);
# the dipole field (4.3) is -grad of p.r/(4 pi eps0 r^3); the spherical capacitor
# tends to the parallel-plate and isolated-sphere limits (B2); the shell
# potential integral (B3) gives kQ/R inside and kQ/r outside.

# %%
r, R, Q, e0, z, a_, b_, d_ = sp.symbols("r R Q epsilon_0 z a b d", positive=True)
rho = Q / (sp.Rational(4, 3) * sp.pi * R**3)
V_in = Q / (8 * sp.pi * e0 * R) * (3 - r**2 / R**2)
V_out = Q / (4 * sp.pi * e0 * r)
lap = lambda f: sp.diff(r**2 * sp.diff(f, r), r) / r**2
assert sp.simplify(lap(V_in) + rho / e0) == 0 and sp.simplify(lap(V_out)) == 0
assert sp.simplify((V_in - V_out).subs(r, R)) == 0
assert sp.simplify(sp.diff(V_in - V_out, r).subs(r, R)) == 0
E_in, E_out = Q * r / (4 * sp.pi * e0 * R**3), Q / (4 * sp.pi * e0 * r**2)
W_rhoV = sp.Rational(1, 2) * sp.integrate(rho * V_in * 4 * sp.pi * r**2, (r, 0, R))
W_in = e0 / 2 * sp.integrate(E_in**2 * 4 * sp.pi * r**2, (r, 0, R))
W_out = e0 / 2 * sp.integrate(E_out**2 * 4 * sp.pi * r**2, (r, R, sp.oo))
q_r = Q * r**3 / R**3                                     # charge inside radius r
W_shells = sp.integrate(q_r / (4 * sp.pi * e0 * r) * sp.diff(q_r, r), (r, 0, R))
W_sph = 3 * Q**2 / (20 * sp.pi * e0 * R)
assert all(sp.simplify(w - W_sph) == 0 for w in (W_rhoV, W_in + W_out, W_shells))
assert sp.simplify(W_in / (W_in + W_out)) == sp.Rational(1, 6)
x, y, zz, px, py, pz = sp.symbols("x y z p_x p_y p_z", real=True)
rv, pv = sp.Matrix([x, y, zz]), sp.Matrix([px, py, pz]); rr = sp.sqrt(x**2 + y**2 + zz**2)
Vdip = pv.dot(rv) / (4 * sp.pi * e0 * rr**3)
Edip = -sp.Matrix([sp.diff(Vdip, s) for s in (x, y, zz)])
Eform = (3 * pv.dot(rv) * rv / rr**2 - pv) / (4 * sp.pi * e0 * rr**3)
assert sp.simplify(Edip - Eform) == sp.zeros(3, 1)
C_sph = 4 * sp.pi * e0 * a_ * b_ / (b_ - a_)
assert sp.simplify(sp.limit(C_sph.subs(b_, a_ + d_) / (e0 * 4 * sp.pi * a_**2 / d_), a_, sp.oo)) == 1
assert sp.limit(C_sph, b_, sp.oo) == 4 * sp.pi * e0 * a_
u = sp.symbols("u")                      # u = cos(theta'); shell potential, B3
shell = lambda rr_: sp.integrate(1 / sp.sqrt(R**2 + rr_**2 - 2 * R * rr_ * u), (u, -1, 1)) / 2
assert abs(float(shell(sp.Rational(3, 10)).subs(R, 1)) - 1.0) < 1e-12        # inside: 1/R
assert abs(float(shell(sp.Integer(3)).subs(R, 1)) - 1 / 3) < 1e-12           # outside: 1/r
print("Part 0: symbolic checks passed (Poisson for the sphere, W = 3Q^2/20 pi eps0 R three ways, "
      "interior share 1/6, dipole field, capacitor limits, shell potential).")

# %% [markdown]
# ## Part 1: Fields of point charges; dipole and quadrupole (§3.1, Example 4.3)
#
# E(r) = (1/4 pi eps0) sum_i q_i (r - r_i)/|r - r_i|^3 (3.2) and V = (1/4 pi eps0)
# sum_i q_i/|r - r_i|. **Build it yourself first:** write `field_potential` with
# broadcasting, no loops over charges.


# %%
def field_potential(q, pos, pts):
    """E (M,3) and V (M,) at points pts (M,3) from charges q (N,) at pos (N,3); eq. (3.2)."""
    sep = pts[:, None, :] - pos[None, :, :]                     # script-r = r - r'
    dist = np.linalg.norm(sep, axis=-1)
    E = K * np.einsum("n,mnk->mk", q, sep / dist[..., None] ** 3)
    return E, K * (q / dist).sum(axis=1)


def plane_plot(ax, q, pos, title, lim=2.0):
    """Streamlines of E and contours of V in the z = 0 plane."""
    g = np.linspace(-lim, lim, 240)              # even count: no grid node sits on a charge
    X, Y = np.meshgrid(g, g)
    E, V = field_potential(q, pos, np.c_[X.ravel(), Y.ravel(), 0 * X.ravel()])
    Ex, Ey, V = E[:, 0].reshape(X.shape), E[:, 1].reshape(X.shape), V.reshape(X.shape)
    lv = K * nC * np.array([-8, -4, -2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2, 4, 8])
    ax.contour(X, Y, V, levels=lv, colors="gray", linewidths=0.7)
    ax.streamplot(X, Y, Ex, Ey, density=1.4, color=np.log(np.hypot(Ex, Ey)), cmap="viridis", linewidth=0.8)
    ax.scatter(pos[:, 0], pos[:, 1], c=["r" if s > 0 else "b" for s in q], s=50, zorder=3)
    ax.set_aspect("equal"); ax.set_title(title)


dip_q, dip_pos = np.array([1, -1.0]) * nC, np.array([[0, 0.5, 0], [0, -0.5, 0.0]])
quad_q = np.array([1, -1, 1, -1.0]) * nC
quad_pos = np.array([[0.5, 0.5, 0], [-0.5, 0.5, 0], [-0.5, -0.5, 0], [0.5, -0.5, 0.0]])
fig, axs = plt.subplots(1, 2, figsize=(11, 5.2))
plane_plot(axs[0], dip_q, dip_pos, "Dipole: field lines (colour = log|E|), equipotentials (grey)")
plane_plot(axs[1], quad_q, quad_pos, "Quadrupole")
fig.tight_layout(); fig.savefig(FIG / "dipole_quadrupole.png", dpi=110); plt.close(fig)

# Example 4.3: dipole along z, charges +-q at +-d/2. On axis at r = 5d the
# leading-order potential is exactly 1 - d^2/(4 r^2) = 0.99 of the true one.
q0, d = nC, 1e-3
pz_pos = np.array([[0, 0, d / 2], [0, 0, -d / 2]]); pz_q = np.array([q0, -q0])
p = np.array([0, 0, q0 * d])
_, Vex = field_potential(pz_q, pz_pos, np.array([[0, 0, 5 * d]]))
Vdip_axis = K * q0 * d / (5 * d) ** 2
assert np.isclose(Vdip_axis / Vex[0], 0.99)
dirs = rng.standard_normal((50, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
pts = 20 * d * dirs
E_num, _ = field_potential(pz_q, pz_pos, pts)
rn = np.linalg.norm(pts, axis=1)[:, None]; rh = pts / rn
E_formula = K * (3 * (rh @ p)[:, None] * rh - p) / rn**3                     # (4.3)
rel = np.linalg.norm(E_num - E_formula, axis=1).max() / np.linalg.norm(E_formula, axis=1).min()
assert rel < 5e-3
print(f"Example 4.3: V_dipole/V_exact on axis at r = 5d: {Vdip_axis/Vex[0]:.4f}; "
      f"max relative error of the (3(p.r)r - p) field at r = 20d: {rel:.1e}")

# Curl-free (Theorem 3.2) and E = -grad V: loop integral of E around a circle
# threading between the quadrupole charges vanishes; central differences of V give E.
t = np.linspace(0, 2 * np.pi, 4001)
loop = np.c_[0.7 * np.cos(t), 0.7 * np.sin(t) + 0.1, 0.05 + 0 * t]
E_loop, _ = field_potential(quad_q, quad_pos, loop)
circ = integrate.trapezoid(np.einsum("ij,ij->i", E_loop, np.gradient(loop, t, axis=0)), t)
pt, hh = np.array([[0.3, -0.2, 0.4]]), 1e-5
Eg = [-(field_potential(quad_q, quad_pos, pt + hh * e)[1] - field_potential(quad_q, quad_pos, pt - hh * e)[1])[0] / (2 * hh)
      for e in np.eye(3)]
assert abs(circ) < 1e-6 * np.abs(E_loop).max() and np.allclose(Eg, field_potential(quad_q, quad_pos, pt)[0][0], rtol=1e-6)

# B4(c): a dipole field line obeys r = r0 sin^2(theta). Follow E from the equator.
def rhs(_, s):
    E, _v = field_potential(pz_q, pz_pos, s[None, :]); return E[0] / np.linalg.norm(E[0])
stop = lambda _, s: np.linalg.norm(s) - 0.1; stop.terminal = True
sol = integrate.solve_ivp(rhs, (0, 5), [1.0, 0, 0], events=stop, rtol=1e-10, atol=1e-12, max_step=0.01)
rl = np.linalg.norm(sol.y, axis=0); thl = np.arccos(sol.y[2] / rl)
assert np.abs(rl / np.sin(thl) ** 2 - 1).max() < 1e-3
print(f"B4: along a dipole field line r/sin^2(theta) stays within {np.abs(rl/np.sin(thl)**2-1).max():.1e} of 1")

# %% [markdown]
# ## Part 2: Gauss's-law checker (Theorem 3.1)
#
# Flux through a sphere (Gauss-Legendre in cos(theta), uniform in phi) and through
# the six faces of a cube (Gauss-Legendre squared), compared with Q_enc/eps0.


# %%
def flux_sphere(q, pos, c, Rs, n=160):
    """Phi = closed-surface integral of E . n da over the sphere |r - c| = Rs."""
    ct, w = roots_legendre(n); ph = 2 * np.pi * (np.arange(2 * n) + 0.5) / (2 * n)
    CT, PH = np.meshgrid(ct, ph, indexing="ij"); ST = np.sqrt(1 - CT**2)
    nhat = np.c_[(ST * np.cos(PH)).ravel(), (ST * np.sin(PH)).ravel(), CT.ravel()]
    E, _ = field_potential(q, pos, c + Rs * nhat)
    wts = (w[:, None] * np.full(2 * n, np.pi / n)[None, :]).ravel() * Rs**2
    return np.sum(np.einsum("ij,ij->i", E, nhat) * wts)


def flux_cube(q, pos, c, L, n=120):
    """Phi through the cube of side L centred at c: sum over six faces, n x n Gauss-Legendre."""
    s, w = roots_legendre(n); s, W2 = s * L / 2, np.outer(w, w).ravel() * (L / 2) ** 2
    U, Vv = [m.ravel() for m in np.meshgrid(s, s, indexing="ij")]
    tot = 0.0
    for ax in range(3):
        for sgn in (1, -1):
            P = np.zeros((U.size, 3)); others = [k for k in range(3) if k != ax]
            P[:, ax] = sgn * L / 2; P[:, others[0]] = U; P[:, others[1]] = Vv
            E, _ = field_potential(q, pos, c + P)
            tot += sgn * np.sum(E[:, ax] * W2)
    return tot


surfaces = [("sphere", np.zeros(3), 1.2), ("sphere", np.array([0.5, 0.3, -0.2]), 0.6),
            ("cube", np.array([0.2, 0.0, 0.0]), 1.6), ("cube", np.array([-0.5, 0.4, 0.3]), 1.1)]
def dist_to_surface(P, kind, c, size):
    if kind == "sphere": return np.abs(np.linalg.norm(P - c, axis=1) - size)
    dd = np.abs(P - c) - size / 2           # distance to the cube surface, inside or outside
    out = np.linalg.norm(np.maximum(dd, 0), axis=1); return np.where(dd.max(1) > 0, out, -dd.max(1))
gq = rng.uniform(-3, 3, 12) * nC; gpos = np.empty((12, 3)); k_ = 0
while k_ < 12:                              # keep every charge at least 0.12 m from every surface
    cand = rng.uniform(-1.0, 1.0, (1, 3)) if k_ < 8 else rng.uniform(-1.6, 1.6, (1, 3))
    if all(dist_to_surface(cand, *s_)[0] > 0.12 for s_ in surfaces): gpos[k_] = cand; k_ += 1
for kind, c, size in surfaces:
    inside = (np.linalg.norm(gpos - c, axis=1) < size) if kind == "sphere" else np.all(np.abs(gpos - c) < size / 2, axis=1)
    assert 0 < inside.sum() < len(gq)
    Phi = (flux_sphere if kind == "sphere" else flux_cube)(gq, gpos, c, size)
    assert abs(Phi - gq[inside].sum() / EPS0) < 1e-7 * np.abs(gq).sum() / EPS0, (kind, Phi)
    print(f"Gauss {kind:6s} size {size}: {inside.sum()} charges inside, flux*eps0 = {Phi*EPS0/nC:+.9f} nC, "
          f"Q_enc = {gq[inside].sum()/nC:+.9f} nC")
# A2: charge at a corner of a cube of side L; flux through one far face is q/(24 eps0).
s, w = roots_legendre(200); s, w = (s + 1) / 2, w / 2                   # nodes on [0, 1]
Yc, Zc = np.meshgrid(s, s, indexing="ij")
far = np.sum(np.outer(w, w) * K * nC / (1 + Yc**2 + Zc**2) ** 1.5)         # E_x on the face x = L = 1
assert np.isclose(far * EPS0 / nC, 1 / 24, rtol=1e-10)
print(f"A2: flux through a far face from a corner charge = q/eps0 x {far*EPS0/nC:.10f} (1/24 = {1/24:.10f})")

# %% [markdown]
# ## Part 3: The uniformly charged sphere: E, V, and energy (Examples 4.1, 4.4)

# %%
Qs, Rs = 1 * nC, 0.1
E_sph = lambda rr_: K * Qs * np.where(rr_ < Rs, rr_ / Rs**3, 1 / np.maximum(rr_, Rs) ** 2)
V_closed = lambda rr_: np.where(rr_ < Rs, K * Qs / (2 * Rs) * (3 - rr_**2 / Rs**2), K * Qs / np.maximum(rr_, Rs))
V_num = lambda r0: integrate.quad(E_sph, r0, Rs)[0] + integrate.quad(E_sph, max(r0, Rs), np.inf)[0] if r0 < Rs \
    else integrate.quad(E_sph, r0, np.inf)[0]
for r0 in (0.0, 0.05, 0.1, 0.2):
    assert np.isclose(V_num(r0), V_closed(r0), rtol=1e-9)
print(f"Example 4.1: V(0) = {V_num(0):.2f} V, V(R) = {V_num(Rs):.2f} V, E(R) = {E_sph(Rs):.1f} V/m, "
      f"V(0)/V(R) = {V_num(0)/V_num(Rs):.3f}")
assert np.isclose(V_num(0) / V_num(Rs), 1.5)
rg = np.linspace(0, 4 * Rs, 400)
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(rg / Rs, E_sph(rg) / E_sph(Rs), label="E / E(R)"); ax.plot(rg / Rs, V_closed(rg) / V_closed(Rs), label="V / V(R)")
ax.axvline(1, color="k", lw=0.5); ax.set_xlabel("r / R"); ax.legend(); ax.set_title("Uniformly charged sphere")
fig.tight_layout(); fig.savefig(FIG / "charged_sphere.png", dpi=110); plt.close(fig)

rho_s = Qs / (4 / 3 * np.pi * Rs**3)
W1 = 0.5 * integrate.quad(lambda rr_: rho_s * V_closed(rr_) * 4 * np.pi * rr_**2, 0, Rs)[0]
Win = EPS0 / 2 * integrate.quad(lambda rr_: E_sph(rr_) ** 2 * 4 * np.pi * rr_**2, 0, Rs)[0]
Wout = EPS0 / 2 * integrate.quad(lambda rr_: E_sph(rr_) ** 2 * 4 * np.pi * rr_**2, Rs, np.inf)[0]
W_exact = 3 * Qs**2 / (20 * np.pi * EPS0 * Rs)
assert np.isclose(W1, W_exact, rtol=1e-10) and np.isclose(Win + Wout, W_exact, rtol=1e-10)
assert np.isclose(Win / W_exact, 1 / 6)
print(f"Example 4.4: W = {W_exact:.4e} J; (1/2)int rho V = {W1:.4e}; (eps0/2)int E^2 = {Win+Wout:.4e}; "
      f"share of E^2 energy inside the ball = {Win/W_exact:.4f}")
# Discrete version: N equal point charges spread uniformly through the ball. The pair
# sum (3.14) has no self-energy terms; it approaches the continuum 3/5 kQ^2/R.
Nmc = 3000
uu = rng.standard_normal((Nmc, 3)); uu /= np.linalg.norm(uu, axis=1)[:, None]
P = Rs * uu * rng.random(Nmc)[:, None] ** (1 / 3)
W_pairs = K * (Qs / Nmc) ** 2 * np.sum(1 / pdist(P))
assert abs(W_pairs / W_exact - 1) < 0.02
print(f"  pair sum over {Nmc} point charges / (3/5 kQ^2/R) = {W_pairs/W_exact:.4f}")
# Point-charge self-energy (eps0/2) int_a^inf E^2 = kq^2/(2a), in units of the electron's rest energy.
ke2_MeVfm = K * 1.602176634e-19**2 / 1.602176634e-28             # e^2/(4 pi eps0) = 1.44 MeV fm
r_e = ke2_MeVfm / 0.51099895                                      # classical electron radius, fm
for a_cut in (1e5, 1.0, 1e-3):                                    # cutoff radius in fm
    print(f"  electron self-energy outside a = {a_cut:g} fm: {r_e/(2*a_cut):.3g} m_e c^2")
# Nuclear Coulomb energy (Example 4.4): a_C = (3/5) e^2/(4 pi eps0 r0) with r0 = 1.2 fm.
r0n = 1.2; a_C = 0.6 * ke2_MeVfm / r0n
EC = lambda Z, A: a_C * Z * (Z - 1) / A ** (1 / 3)
print(f"  e^2/4pi eps0 = {ke2_MeVfm:.4f} MeV fm; a_C = {a_C:.3f} MeV; U-238: {EC(92, 238):.0f} MeV; "
      f"two Pd-119: {2*EC(46, 119):.0f} MeV; released {EC(92, 238)-2*EC(46, 119):.0f} MeV")
assert np.isclose(a_C, 0.720, atol=1e-3) and abs(EC(92, 238) - 973) < 1 and abs(2 * EC(46, 119) - 606) < 1
print(f"C1: electrostatic electron radius: ball {0.6*r_e:.3f} fm, shell {0.5*r_e:.3f} fm (r_e = {r_e:.4f} fm)")

# %% [markdown]
# ## Part 4: The infinite line and its reference point (Example 4.2)
#
# A finite line of half-length L: V(s) = (lambda/2 pi eps0) asinh(L/s). V itself
# grows like ln L; differences converge to (lambda/2 pi eps0) ln(s2/s1).

# %%
lam = 1 * nC
for L in (1.0, 1e3, 1e6):
    Vs = lambda s_: integrate.quad(lambda zz_: K * lam / np.hypot(s_, zz_), -L, L, points=[0], limit=200)[0]
    V1, V2 = Vs(0.01), Vs(0.10)
    assert np.isclose(V1, lam / (2 * np.pi * EPS0) * np.arcsinh(L / 0.01), rtol=1e-8)
    print(f"Example 4.2: L = {L:g} m: V(1 cm) = {V1:.2f} V, V(10 cm) = {V2:.2f} V, difference = {V1-V2:.3f} V")
dV_inf = lam / (2 * np.pi * EPS0) * np.log(10)
print(f"  infinite line: lambda/(2 pi eps0) = {lam/(2*np.pi*EPS0):.3f} V; V(1 cm) - V(10 cm) = {dV_inf:.3f} V")
assert np.isclose(dV_inf, 41.39, atol=0.005)

# %% [markdown]
# ## Part 5: Poisson on a grid: the parallel-plate capacitor (Example 4.5, §3.11)
#
# Five-point stencil for (3.9): V_ij = average of its 4 neighbours + rho h^2/(4 eps0).
# Conductors and the grounded outer box are "fixed" nodes. We solve the sparse
# system directly (LA-13), then compute: the field at the centre (vs V0/d), the
# surface charge from the jump in dV/dn (3.12), and C two ways: from the energy
# (eps0/2) sum (dV)^2 over grid edges (3.16) and from Q/V0. Everything is per unit
# depth (2D), plate width a = 1 m, V0 = 1 V.


# %%
def solve_dirichlet(V, fixed, h, rho=None):
    """Solve the 5-point discretisation of lap V = -rho/eps0 (3.9) for the free nodes;
    fixed nodes (conductors, outer boundary) keep their values. Returns the full grid."""
    free = ~fixed; nf = free.sum()
    idx = -np.ones(V.shape, int); idx[free] = np.arange(nf)
    rhs = np.zeros(nf) if rho is None else -rho[free] * h * h / EPS0
    jj, ii = np.nonzero(free); me = idx[jj, ii]
    rows, cols, vals = [me], [me], [-4.0 * np.ones(nf)]
    for dj, di in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        nj, ni = jj + dj, ii + di; nb = free[nj, ni]
        rows.append(me[nb]); cols.append(idx[nj[nb], ni[nb]]); vals.append(np.ones(nb.sum()))
        np.subtract.at(rhs, me[~nb], V[nj[~nb], ni[~nb]])
    A = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nf, nf))
    out = V.copy(); out[free] = spsolve(A, rhs); return out


def node_charge(V):
    """Discrete Gauss's law: charge per unit depth at each node, Q_ij = -eps0 (sum of neighbours - 4 V_ij)."""
    Qn = np.zeros_like(V)
    Qn[1:-1, 1:-1] = -EPS0 * (V[2:, 1:-1] + V[:-2, 1:-1] + V[1:-1, 2:] + V[1:-1, :-2] - 4 * V[1:-1, 1:-1])
    return Qn


def edge_energy(V):
    """(eps0/2) int |grad V|^2 dA on the grid = (eps0/2) sum over edges of (Delta V)^2 (the h^2 cancels)."""
    return 0.5 * EPS0 * (np.sum(np.diff(V, axis=0) ** 2) + np.sum(np.diff(V, axis=1) ** 2))


def capacitor_grid(d, h=0.01, half=1.5, a=1.0, V0=1.0):
    n = int(round(2 * half / h)) + 1; xg = np.linspace(-half, half, n)
    V = np.zeros((n, n)); fixed = np.zeros((n, n), bool)
    fixed[[0, -1], :] = fixed[:, [0, -1]] = True
    jt, jb = np.argmin(abs(xg - d / 2)), np.argmin(abs(xg + d / 2))
    assert np.isclose(xg[jt] - xg[jb], d), "d must be a multiple of 2h"
    plate = np.abs(xg) <= a / 2 + 1e-9
    for j, v in ((jt, V0 / 2), (jb, -V0 / 2)): fixed[j, plate] = True; V[j, plate] = v
    return xg, solve_dirichlet(V, fixed, h), jt, jb, plate


V0 = 1.0
rows_out = []
for d_ in (0.32, 0.16, 0.08, 0.04):
    xg, Vg, jt, jb, plate = capacitor_grid(d_)
    h = xg[1] - xg[0]; mid = len(xg) // 2
    Qn = node_charge(Vg); Q_top = Qn[jt, plate].sum()
    C_E, C_Q, C_id = 2 * edge_energy(Vg) / V0**2, Q_top / V0, EPS0 * 1.0 / d_
    assert np.isclose(C_E, C_Q, rtol=1e-9)                 # discrete Green's identity: (3.15) = (3.16)
    E_c = (Vg[mid + 1, mid] - Vg[mid - 1, mid]) / (2 * h)  # = -E_y at the centre
    sig_in = EPS0 * (Vg[jt, plate] - Vg[jt - 1, plate]) / h   # (3.12) on the gap side of the top plate
    sig_out = EPS0 * (Vg[jt, plate] - Vg[jt + 1, plate]) / h
    assert np.isclose(E_c * d_ / V0, 1, atol=1e-3) and np.isclose(sig_in[len(sig_in) // 2] / (EPS0 * V0 / d_), 1, atol=1e-3)
    rows_out.append((d_, C_E / C_id, sig_out[len(sig_out) // 2] / sig_in[len(sig_in) // 2]))
    if d_ == 0.32: keep = (xg, Vg, sig_in, sig_out, plate, d_)
for d_, ratio, so in rows_out:
    print(f"Capacitor d/a = {d_:.2f}: C/C_ideal = {ratio:.3f}, centre field = V0/d, "
          f"outer/inner surface charge at the centre = {so:.3f}")
ratios = [rw[1] for rw in rows_out]
assert all(np.diff(ratios) < 0) and ratios[-1] < 1.1 and ratios[0] > 1.4
xg, Vg, sig_in, sig_out, plate, dk = keep
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
Xg, Yg = np.meshgrid(xg, xg); Ey, Ex = np.gradient(-Vg, xg, xg)
axs[0].contour(Xg, Yg, Vg, levels=np.linspace(-0.45, 0.45, 19), colors="gray", linewidths=0.6)
axs[0].streamplot(Xg, Yg, Ex, Ey, density=1.6, color="C0", linewidth=0.7)
axs[0].plot([-0.5, 0.5], [dk / 2] * 2, "r", lw=3); axs[0].plot([-0.5, 0.5], [-dk / 2] * 2, "b", lw=3)
axs[0].set_xlim(-1.2, 1.2); axs[0].set_ylim(-1, 1); axs[0].set_aspect("equal"); axs[0].set_title(f"Fringing field, d/a = {dk}")
axs[1].plot(xg[plate], sig_in / (EPS0 * V0 / dk), label="gap side"); axs[1].plot(xg[plate], sig_out / (EPS0 * V0 / dk), label="outer side")
axs[1].set_xlabel("x / a"); axs[1].set_ylabel("sigma / (eps0 V0 / d)"); axs[1].legend(); axs[1].set_title("Surface charge on the top plate")
fig.tight_layout(); fig.savefig(FIG / "capacitor_fringing.png", dpi=110); plt.close(fig)

# Relaxation is slow (§8, LA-13): Jacobi on the d/a = 0.2 capacitor at two resolutions.
def jacobi(V, fixed, tol, ref):
    """Jacobi sweeps V_ij <- mean of 4 neighbours until max|V - ref| < tol; returns the sweep count."""
    V = V.copy(); free = ~fixed; k = 0
    while np.abs(V - ref).max() > tol:
        avg = 0.25 * (V[2:, 1:-1] + V[:-2, 1:-1] + V[1:-1, 2:] + V[1:-1, :-2])
        V[1:-1, 1:-1] = np.where(free[1:-1, 1:-1], avg, V[1:-1, 1:-1]); k += 1
    return k
its = []
for hj in (0.05, 0.025):
    xg_, Vd, jt_, jb_, pl_ = capacitor_grid(0.2, h=hj)
    V_start = np.zeros_like(Vd); fx = np.zeros(Vd.shape, bool); fx[[0, -1], :] = fx[:, [0, -1]] = True
    for j, v in ((jt_, 0.5), (jb_, -0.5)): fx[j, pl_] = True; V_start[j, pl_] = v
    its.append(jacobi(V_start, fx, 1e-6, Vd))
print(f"Jacobi sweeps to 1e-6 V: h = 0.05: {its[0]}, h = 0.025: {its[1]} (ratio {its[1]/its[0]:.2f}; LA-13 predicts ~4)")
assert 3.2 < its[1] / its[0] < 4.5

# Example 4.5 numbers (3D plates, ideal formulas).
A_p, d_p, V_p = 0.01, 1e-3, 100.0
C_p = EPS0 * A_p / d_p; Q_p = C_p * V_p; W_p = 0.5 * C_p * V_p**2; F_p = Q_p**2 / (2 * EPS0 * A_p)
print(f"Example 4.5: C = {C_p*1e12:.2f} pF, Q = {Q_p*1e9:.3f} nC, sigma = {Q_p/A_p*1e6:.4f} uC/m^2, E = {V_p/d_p:.0f} V/m, "
      f"W = {W_p:.3e} J, F = {F_p:.3e} N = W/d = {W_p/d_p:.3e} N, pressure {F_p/A_p:.4f} Pa, weight of {F_p/9.81*1e6:.0f} mg")
assert np.isclose(F_p, W_p / d_p) and np.isclose(Q_p / A_p / EPS0, V_p / d_p)

# %% [markdown]
# ## Part 6: A grounded conductor in a uniform field, by relaxation (§3.10)
#
# A grounded conducting cylinder of radius a (a disk in the 2D cross-section) in
# a uniform field E0 x-hat. Red-black successive over-relaxation (Gauss-Seidel
# with over-correction omega; LA-13). The outer box carries the exact
# solution V = -E0 (s - a^2/s) cos(phi) (EM-3 derives it), so the only error is
# the staircase disk. Checks: V on the disk is 0 by construction; the field
# meets it perpendicularly; the induced charge is sigma = 2 eps0 E0 cos(phi).


# %%
def sor_red_black(V, fixed, omega, tol=1e-10, maxit=20000):
    """Red-black SOR for lap V = 0: each colour V <- V + omega (mean of neighbours - V)."""
    J, I = np.indices(V.shape); inner = np.zeros(V.shape, bool); inner[1:-1, 1:-1] = True
    colours = [((I + J) % 2 == c) & ~fixed & inner for c in (0, 1)]
    for k in range(maxit):
        dmax = 0.0
        for m in colours:
            avg = np.zeros_like(V)
            avg[1:-1, 1:-1] = 0.25 * (V[2:, 1:-1] + V[:-2, 1:-1] + V[1:-1, 2:] + V[1:-1, :-2])
            delta = omega * (avg - V) * m; V += delta; dmax = max(dmax, np.abs(delta).max())
        if dmax < tol: return V, k + 1
    raise RuntimeError("SOR did not converge")


a_d, E0, Lb, hd = 1.0, 1.0, 3.0, 1.0 / 40
xd = np.linspace(-Lb, Lb, int(round(2 * Lb / hd)) + 1); Xd, Yd = np.meshgrid(xd, xd)
Sd, PHd = np.hypot(Xd, Yd), np.arctan2(Yd, Xd)
fixed_d = Sd <= a_d + 1e-9                 # tolerance keeps the mask mirror-symmetric
V_ex = np.where(~fixed_d, -E0 * (Sd - a_d**2 / np.maximum(Sd, 1e-12)) * np.cos(PHd), 0.0)
cond = fixed_d.copy(); fixed_d[[0, -1], :] = fixed_d[:, [0, -1]] = True
Vd, n_sor = sor_red_black(np.where(fixed_d, V_ex, 0.0), fixed_d, 2 / (1 + np.sin(np.pi * hd / (2 * Lb))))
err_d = np.abs(Vd - V_ex)[~cond].max()
Qd = node_charge(Vd)
c1 = np.sum(Qd[cond] * np.cos(PHd[cond])) / (np.pi * a_d) / (2 * EPS0 * E0)   # cos(phi) Fourier coefficient of sigma
assert abs(Qd[cond].sum()) < 1e-12 * EPS0 * E0 * 1e3 and abs(c1 - 1) < 0.02 and err_d < 0.05
gy, gx = np.gradient(Vd, hd)
Ix, Iy = RegularGridInterpolator((xd, xd), -gx.T), RegularGridInterpolator((xd, xd), -gy.T)
phs = np.linspace(-np.pi, np.pi, 721)[:-1]; tang = []
for dl in (0.2, 0.1, 0.05):
    sR = a_d + dl; Pp = np.c_[sR * np.cos(phs), sR * np.sin(phs)]
    Exx, Eyy = Ix(Pp), Iy(Pp)
    Es, Ep = Exx * np.cos(phs) + Eyy * np.sin(phs), -Exx * np.sin(phs) + Eyy * np.cos(phs)
    num, ex = np.sqrt(np.mean(Ep**2) / np.mean(Es**2)), (1 - a_d**2 / sR**2) / (1 + a_d**2 / sR**2)
    tang.append(num); assert abs(num - ex) < 0.02
assert tang[0] > tang[1] > tang[2]
print(f"Disk in uniform field: SOR {n_sor} sweeps on {len(xd)}^2 grid, max |V - exact| = {err_d:.3f} E0 a; "
      f"sigma = {c1:.3f} x 2 eps0 E0 cos(phi); rms E_tan/E_norm at s = a + 0.2, 0.1, 0.05: "
      + ", ".join(f"{t_:.3f}" for t_ in tang))
fig, ax = plt.subplots(figsize=(6, 6))
ax.contour(Xd, Yd, Vd, levels=np.linspace(-2.6, 2.6, 27), colors="gray", linewidths=0.6)
ax.streamplot(Xd, Yd, -gx, -gy, density=1.5, linewidth=0.7)
ax.add_patch(plt.Circle((0, 0), a_d, color="k", alpha=0.3)); ax.set_aspect("equal")
ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5); ax.set_title("Grounded conductor in a uniform field (SOR)")
fig.tight_layout(); fig.savefig(FIG / "conductor_uniform_field.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 7: A charge in an off-centre cavity (Example 4.6)
#
# Spherical cavity of radius a, charge q at distance b from its centre. The wall
# charge is sigma(theta) = -q (a^2 - b^2) / (4 pi a (a^2 + b^2 - 2ab cos theta)^{3/2})
# (EM-3 gets this from an image charge -q a/b at a^2/b; here we check that the
# image field reproduces it, that it totals -q, and the max/min ratio).

# %%
qc, ac, bc, Rc = 1 * nC, 0.02, 0.01, 0.10
sig_wall = lambda th: -qc * (ac**2 - bc**2) / (4 * np.pi * ac * (ac**2 + bc**2 - 2 * ac * bc * np.cos(th)) ** 1.5)
Qwall = integrate.quad(lambda th: sig_wall(th) * 2 * np.pi * ac**2 * np.sin(th), 0, np.pi)[0]
th_ = np.linspace(0, np.pi, 7); wall = np.c_[ac * np.sin(th_), 0 * th_, ac * np.cos(th_)]
E_img, _ = field_potential(np.array([qc, -qc * ac / bc]), np.array([[0, 0, bc], [0, 0, ac**2 / bc]]), wall)
sig_img = -EPS0 * np.einsum("ij,ij->i", E_img, wall / ac)            # n-hat points into the cavity
assert np.isclose(Qwall, -qc, rtol=1e-10) and np.allclose(sig_img, sig_wall(th_), rtol=1e-9)
ratio_wall = sig_wall(0) / sig_wall(np.pi)
assert np.isclose(ratio_wall, 27)
print(f"Example 4.6: wall charge = {Qwall/nC:.6f} nC, sigma(near)/sigma(far) = {ratio_wall:.1f}; outer sigma = "
      f"{qc/(4*np.pi*Rc**2)*1e9:.2f} nC/m^2 uniform; conductor potential = {K*qc/Rc:.2f} V")

# %% [markdown]
# ## Part 8: Numbers quoted in solutions.md

# %%
E_A1, V_A1 = field_potential(np.array([2, -1.0]) * nC, np.array([[0, 0, 0], [0.3, 0, 0.0]]), np.array([[0, 0.4, 0.0]]))
print(f"A1: E = ({E_A1[0,0]:.2f}, {E_A1[0,1]:.2f}) V/m, |E| = {np.linalg.norm(E_A1):.2f} V/m, V = {V_A1[0]:.2f} V")
sig4 = 5e-6; print(f"A4: sigma/eps0 = {sig4/EPS0:.4e} V/m; E_above = (3, 0, {2+sig4/EPS0/1e5:.3f}) x 1e5 V/m")
Cs = 1 / (1 / 2e-6 + 1 / 3e-6); Ceq = Cs + 4e-6
print(f"A5: series {Cs*1e6:.2f} uF, C_eq = {Ceq*1e6:.2f} uF, Q_tot = {Ceq*12*1e6:.1f} uC, Q(2uF) = {Cs*12*1e6:.1f} uC, "
      f"W = {0.5*Ceq*144*1e6:.1f} uJ")
sq = np.array([[0, 0, 0], [0.1, 0, 0], [0.1, 0.1, 0], [0, 0.1, 0.0]]); qsq = np.full(4, nC)
W_sq = K * nC**2 * np.sum(1 / pdist(sq)); _, Vsq = field_potential(qsq[1:], sq[1:], sq[:1])
assert np.isclose(W_sq, 0.5 * 4 * nC * Vsq[0]) and np.isclose(W_sq, (4 + np.sqrt(2)) * K * nC**2 / 0.1)
print(f"A6: W = {W_sq:.4e} J = (4 + sqrt 2) kq^2/a")
GM, RE = 3.986004e14, 6.371e6
print(f"B1: tunnel period 2 pi sqrt(R^3/GM) = {2*np.pi*np.sqrt(RE**3/GM):.0f} s = {2*np.pi*np.sqrt(RE**3/GM)/60:.1f} min, "
      f"speed at centre {np.sqrt(GM/RE)/1e3:.2f} km/s")
sigE = -EPS0 * 100
print(f"B2: Earth C = {4*np.pi*EPS0*RE*1e6:.0f} uF; fair-weather sigma = {sigE:.3e} C/m^2, Q = {sigE*4*np.pi*RE**2:.2e} C")
Qb6, Rb6 = 1e-6, 0.1
F_hemi = integrate.quad(lambda th: (Qb6 / (4 * np.pi * Rb6**2)) ** 2 / (2 * EPS0) * np.cos(th) * 2 * np.pi * Rb6**2 * np.sin(th), 0, np.pi / 2)[0]
assert np.isclose(F_hemi, Qb6**2 / (32 * np.pi * EPS0 * Rb6**2))
print(f"B6: force on the northern hemisphere = {F_hemi:.4f} N")
print(f"C2: predicted Jacobi sweeps for a 1e-6 reduction, N = 100: {np.log(1e-6)/np.log(np.cos(np.pi/101)):.0f}; "
      f"N = 30: {np.log(1e-6)/np.log(np.cos(np.pi/31)):.0f}")
Nsq = 30; Vsqr = np.zeros((Nsq + 2, Nsq + 2)); Vsqr[-1, :] = 1.0; fsq = np.zeros_like(Vsqr, bool)
fsq[[0, -1], :] = fsq[:, [0, -1]] = True
ref_sq = solve_dirichlet(Vsqr, fsq, 1.0)
n_meas = jacobi(Vsqr, fsq, 1e-6 * np.abs(ref_sq).max(), ref_sq)
assert abs(n_meas / (np.log(1e-6) / np.log(np.cos(np.pi / 31))) - 1) < 0.15
print(f"C2: measured Jacobi sweeps, N = 30, error down by 1e-6: {n_meas}")

# %% [markdown]
# ## Your turn
#
# 1. **Ring and disk on the axis** (Exercise B8). Closed forms: ring E_z =
#    kQz/(z^2+R^2)^{3/2}; disk E_z = (sigma/2 eps0)(1 - z/sqrt(z^2+R^2)). Build
#    the disk out of point charges on a polar grid, sum with `field_potential`,
#    and compare on the axis. Where is the ring's field largest?
# 2. **Cylindrical capacitor** (Exercise B2(d)). Use `solve_dirichlet` on a grid
#    with an inner disk s <= a at V0 and everything s >= b grounded; compute C'
#    from `edge_energy` and compare with 2 pi eps0 / ln(b/a). Why is the error a
#    few percent, and which way does it go?
# 3. **The field inside a uniformly charged shell vanishes.** Put N point charges
#    on a Fibonacci lattice on a sphere and evaluate E at random interior points.
#    How does the residual depend on N and on the distance from the shell?
#
# Reference solutions follow. Write your own before reading them.

# %%
Qr, Rr = 1 * nC, 0.1
ring_E = lambda zz_: K * Qr * zz_ / (zz_**2 + Rr**2) ** 1.5
disk_E = lambda zz_: Qr / (np.pi * Rr**2) / (2 * EPS0) * (1 - zz_ / np.sqrt(zz_**2 + Rr**2))
nr, nphi = 400, 400                                 # midpoint rule in s and phi
sm = (np.arange(nr) + 0.5) * Rr / nr; pm = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
SM, PM = np.meshgrid(sm, pm, indexing="ij")
dq = (Qr / (np.pi * Rr**2)) * SM * (Rr / nr) * (2 * np.pi / nphi)
dpos = np.c_[(SM * np.cos(PM)).ravel(), (SM * np.sin(PM)).ravel(), 0 * SM.ravel()]
zs = np.array([0.01, 0.05, 0.1, 1.0])
Ed, _ = field_potential(dq.ravel(), dpos, np.c_[0 * zs, 0 * zs, zs])
print("  disk superposition rel. errors:", np.abs(Ed[:, 2] / disk_E(zs) - 1)); assert np.allclose(Ed[:, 2], disk_E(zs), rtol=1e-3)
zr = minimize_scalar(lambda zz_: -ring_E(zz_), bounds=(0, Rr), method="bounded").x
assert np.isclose(zr, Rr / np.sqrt(2), rtol=1e-4)
print(f"Your turn 1 / B8: ring E(0.1 m) = {ring_E(0.1):.1f} V/m, max {ring_E(zr):.1f} V/m at z = {zr:.4f} m; "
      f"disk E at z = 0.01, 0.05, 0.1, 1 m: " + ", ".join(f"{e:.2f}" for e in disk_E(zs))
      + f" V/m (sigma/2eps0 = {Qr/(np.pi*Rr**2)/(2*EPS0):.1f}; kQ/z^2 at 1 m = {K*Qr:.3f})")
ac_, bc_, hc = 0.25, 1.0, 0.01
xc = np.linspace(-1.05, 1.05, int(round(2.1 / hc)) + 1); Xc, Yc_ = np.meshgrid(xc, xc); Sc = np.hypot(Xc, Yc_)
fixed_c = (Sc <= ac_ + 1e-9) | (Sc >= bc_ - 1e-9); Vc = np.where(Sc <= ac_ + 1e-9, 1.0, 0.0)
Vc = solve_dirichlet(Vc, fixed_c, hc)
C_num, C_th = 2 * edge_energy(Vc), 2 * np.pi * EPS0 / np.log(bc_ / ac_)
assert abs(C_num / C_th - 1) < 0.05
print(f"Your turn 2 / B2(d): coax a = {ac_}, b = {bc_}: C' numerical = {C_num*1e12:.2f} pF/m, "
      f"2 pi eps0/ln(b/a) = {C_th*1e12:.2f} pF/m, ratio {C_num/C_th:.4f}")
Nf = 4000; ii_ = np.arange(Nf) + 0.5; thf, phf = np.arccos(1 - 2 * ii_ / Nf), np.pi * (1 + 5**0.5) * ii_
shell_pos = np.c_[np.sin(thf) * np.cos(phf), np.sin(thf) * np.sin(phf), np.cos(thf)]
uin = rng.standard_normal((300, 3)); uin /= np.linalg.norm(uin, axis=1)[:, None]
for rin in (0.3, 0.5, 0.9, 2.0):
    Esh, _ = field_potential(np.full(Nf, 1.0 / Nf), shell_pos, rin * uin)
    ratio_sh = np.linalg.norm(Esh, axis=1).max() / K          # in units of kQ/R^2, Q = R = 1
    print(f"Your turn 3: shell of {Nf} charges, max |E| at r = {rin} R: {ratio_sh:.2e} kQ/R^2")
    if rin <= 0.5: assert ratio_sh < 1e-4
    if rin == 2.0: assert np.isclose(ratio_sh, 0.25, rtol=1e-3)

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.perf_counter() - T0:.1f} s. Figures in {FIG}")
