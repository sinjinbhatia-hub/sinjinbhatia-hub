# %% [markdown]
# # EM-4 Project: Bound Charges, the Polarized Sphere, and Dielectrics on a Grid
#
# This script checks every number in the EM-4 notes and solutions and builds the
# module's three pictures. Parts: (1) bound charges from a given P, symbolically
# (SymPy divergence) and on a grid, with the total asserted zero; (2) the
# uniformly polarized sphere three ways: surface patches of sigma_b = P cos(theta),
# two offset uniformly charged balls, and a cubic lattice of point dipoles (where
# the microscopic field at a lattice site is zero and only the average is
# -P/3eps0); the averaging lemma of section 3.7 checked by quadrature;
# (3) a 2D finite-difference solver for div(eps grad V) = -rho_f with harmonic-mean
# eps on faces, used for the half-filled capacitor both ways; (4) a dielectric
# cylinder and a cavity in a uniform field against 2 E0/(eps_r + 1), with the error
# split into box truncation and staircasing, and the sphere against 3 E0/(eps_r + 2)
# (SymPy); (5) energy 1/2 int D.E vs 1/2 C V^2; (6) the force on a partially inserted
# slab from C(x), at fixed V and fixed Q, and where along the slab it is delivered;
# (7) Clausius-Mossotti for N2, argon, water (with Debye's term and Onsager's repair),
# and the numbers in every worked example and exercise; "Your turn" with reference
# solutions (coaxial cable analytically and with the solver, water, quartz).
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from pathlib import Path
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
EPS0 = 8.8541878128e-12          # F/m
K_E = 1 / (4 * np.pi * EPS0)     # N m^2 / C^2
E_CH = 1.602176634e-19           # C
K_B = 1.380649e-23               # J/K
A3 = 1e-30                       # one cubic angstrom in m^3

# %% [markdown]
# ## Part 1: Bound charges from P, symbolically and on a grid
#
# (3.8): sigma_b = P.n, rho_b = -div P. (3.10): the total bound charge of a bounded
# object is zero because int(-div P) = -closed_int P.da cancels the surface term.
# Also the crude-atom polarizability (3.2) and the identity grad'(1/scr) = scr_hat/scr^2 (3.9).

# %%
x, y, z, k, R, a, q, eps0, Eext, d = sp.symbols("x y z k R a q epsilon_0 E d", positive=True)


def div(F):
    """Cartesian divergence of a 3-component SymPy vector field."""
    return sum(sp.diff(Fi, v) for Fi, v in zip(F, (x, y, z)))


# (3.2): a uniform ball of charge -q, radius a, pulls the nucleus back with E = q d/(4 pi eps0 a^3)
# (EM-2 Example 4.1's interior field). Equilibrium with an external E fixes d and hence p = q d.
d_eq = sp.solve(sp.Eq(Eext, q * d / (4 * sp.pi * eps0 * a**3)), d)[0]
assert sp.simplify(q * d_eq / Eext - 4 * sp.pi * eps0 * a**3) == 0, "alpha = 4 pi eps0 a^3"

# (3.9): the source-point gradient of 1/scr is scr_hat/scr^2 = (r - r')/scr^3.
X, Y, Z, xp, yp, zp = sp.symbols("X Y Z x' y' z'", real=True)
scr = sp.sqrt((X - xp)**2 + (Y - yp)**2 + (Z - zp)**2)
for v, comp in zip((xp, yp, zp), (X - xp, Y - yp, Z - zp)):
    assert sp.simplify(sp.diff(1 / scr, v) - comp / scr**3) == 0

# Two radial polarizations in a sphere of radius R: rho_b, sigma_b, and the zero total (Exercises B1 and B3).
rr = sp.sqrt(x**2 + y**2 + z**2)
for P_r, name in [(k * rr, "k r"), (k * rr**2, "k r^2")]:
    P_vec = [P_r * x / rr, P_r * y / rr, P_r * z / rr]
    rho_b = sp.simplify(-div(P_vec))
    rho_b_r = rho_b.subs({x: 0, y: 0, z: sp.Symbol("r", positive=True)})
    Q_vol = sp.integrate(rho_b_r * 4 * sp.pi * sp.Symbol("r", positive=True)**2, (sp.Symbol("r", positive=True), 0, R))
    Q_surf = P_r.subs({x: 0, y: 0, z: R}) * 4 * sp.pi * R**2       # sigma_b = P.r_hat = P(R)
    assert sp.simplify(Q_vol + Q_surf) == 0
    print(f"P = {name} r_hat:  rho_b = {rho_b_r},  sigma_b = {sp.simplify(P_r.subs({x: 0, y: 0, z: R}))},  "
          f"Q_vol = {sp.simplify(Q_vol)}, Q_surf = {sp.simplify(Q_surf)}")
assert sp.simplify(-div([k * x, k * y, k * z]) + 3 * k) == 0                 # P = k r r_hat: rho_b = -3k
assert sp.simplify(-div([k * x, 0, 0]) + k) == 0                              # Exercise A2: P = k x x_hat


def bound_charge_on_grid(P_fn, L=1.5, n=121):
    """rho_b = -div P by central differences on an n^3 grid of half-width L, and its integral.

    P_fn(X, Y, Z) returns the three Cartesian components. Returns (Q_b, int |rho_b|, rho_b array, h)."""
    s = np.linspace(-L, L, n)
    h = s[1] - s[0]
    Xg, Yg, Zg = np.meshgrid(s, s, s, indexing="ij")
    Px, Py, Pz = P_fn(Xg, Yg, Zg)
    divP = np.gradient(Px, h, axis=0) + np.gradient(Py, h, axis=1) + np.gradient(Pz, h, axis=2)
    rho = -divP
    return rho.sum() * h**3, np.abs(rho).sum() * h**3, rho, h


def P_blob(Xg, Yg, Zg, P0=1.0, Rb=1.0):
    """A bounded, C^1 polarization: P = P0 (1 - r^2/R^2)^2 (x, 0.3 y, z) inside r < R, zero outside."""
    r2 = Xg**2 + Yg**2 + Zg**2
    f = np.where(r2 < Rb**2, P0 * (1 - r2 / Rb**2)**2, 0.0)
    return f * Xg, 0.3 * f * Yg, f * Zg


Qb, Qabs, rho_grid, h_grid = bound_charge_on_grid(P_blob)
assert Qabs > 1.0 and abs(Qb) < 1e-3 * Qabs, (Qb, Qabs)
# Compare the grid rho_b with SymPy at one interior point.
f_sym = (1 - rr**2)**2
rho_sym = sp.lambdify((x, y, z), -div([f_sym * x, sp.Rational(3, 10) * f_sym * y, f_sym * z]))
s = np.linspace(-1.5, 1.5, 121)
i0 = np.searchsorted(s, 0.3); j0 = np.searchsorted(s, -0.2); k0 = np.searchsorted(s, 0.4)
assert abs(rho_grid[i0, j0, k0] - rho_sym(s[i0], s[j0], s[k0])) < 2e-3 * np.abs(rho_grid).max()
print(f"Grid bound charge of a bounded blob: Q_b = {Qb:.2e} vs int|rho_b| = {Qabs:.3f}  (zero to {abs(Qb)/Qabs:.1e})")

# Exercise B2: the dipole moment of the bound charges is the integral of P, int r rho_b dtau = int P dtau (no surface term here: P -> 0 smoothly).
def P_blob_offset(Xg, Yg, Zg, offset=0.5):
    """P = f(r) (offset + x, 0.3 y, z) with f = (1 - r^2)^2 for r < 1: int P_x dtau = offset * int f = offset * 32 pi/105."""
    r2 = Xg**2 + Yg**2 + Zg**2
    f = np.where(r2 < 1.0, (1 - r2)**2, 0.0)
    return f * (offset + Xg), 0.3 * f * Yg, f * Zg


s_off = np.linspace(-1.5, 1.5, 121)
Xo, Yo, Zo = np.meshgrid(s_off, s_off, s_off, indexing="ij")
for off_ in (0.5, -1.0):                                                   # Exercise B2(e): predict p = offset * 32 pi/105
    _, _, rho_off, h_off = bound_charge_on_grid(lambda X, Y, Z: P_blob_offset(X, Y, Z, off_))
    p_bound = np.array([(c * rho_off).sum() * h_off**3 for c in (Xo, Yo, Zo)])
    assert np.isclose(p_bound[0], off_ * 32 * np.pi / 105, rtol=2e-3) and np.abs(p_bound[1:]).max() < 1e-12, p_bound
    print(f"B2: offset {off_:+.1f}: dipole moment of the grid bound charge = ({p_bound[0]:+.5f}, {p_bound[1]:.1e}, {p_bound[2]:.1e}); "
          f"int P dtau = offset * 32 pi/105 = {off_ * 32 * np.pi / 105:+.5f}")

# Exercise A2: P = k x x_hat in the cube 0 <= x, y, z <= L; k = 1e-4 C/m^3, L = 2 cm.
k_A2, L_A2 = 1e-4, 0.02
Qv_A2 = -k_A2 * L_A2**3; Qs_A2 = k_A2 * L_A2 * L_A2**2
assert np.isclose(Qv_A2 + Qs_A2, 0) and np.isclose(Qs_A2, 8e-10) and np.isclose(k_A2 * L_A2, 2e-6)
print(f"A2: rho_b = {-k_A2:.1e} C/m^3, sigma_b(x = L) = {k_A2 * L_A2:.1e} C/m^2, Q_vol = {Qv_A2:.1e} C, Q_surf = {Qs_A2:.1e} C")

# Exercises B1 and B3: P = k r^n r_hat in a sphere of radius R. Gauss for E with the bound charge, then D = 0, then V(0) two ways.
rs, us = sp.symbols("r u", positive=True)
for n_, name in [(1, "k r"), (2, "k r^2")]:
    rho_b_r = -sp.diff(us**2 * k * us**n_, us) / us**2                              # -(1/r^2) d(r^2 P_r)/dr
    Q_enc = sp.integrate(rho_b_r * 4 * sp.pi * us**2, (us, 0, rs))
    E_in_r = sp.simplify(Q_enc / (4 * sp.pi * eps0 * rs**2))
    assert sp.simplify(E_in_r + k * rs**n_ / eps0) == 0                                # E_in = -P/eps0, so D = 0
    sigma_R = k * R**n_
    assert sp.simplify(Q_enc.subs(rs, R) + sigma_R * 4 * sp.pi * R**2) == 0           # E_out = 0
    V0_E = sp.integrate(E_in_r.subs(rs, us), (us, 0, R))                              # V(0) = int_0^R E dr, E_out = 0
    V0_direct = (sigma_R * 4 * sp.pi * R**2 / R + sp.integrate(rho_b_r * 4 * sp.pi * us**2 / us, (us, 0, R))) / (4 * sp.pi * eps0)
    assert sp.simplify(V0_E - V0_direct) == 0
    print(f"P = {name} r_hat: E_in = {E_in_r}, V(0) = {sp.simplify(V0_E)} (from E) = {sp.simplify(V0_direct)} (from rho_b, sigma_b)")
k_B1, R_B, k_B3 = 1e-4, 0.01, 1e-2
assert np.isclose(-k_B1 * R_B / EPS0, -1.1294e5, rtol=1e-4) and np.isclose(-3 * k_B1, -3e-4)
assert np.isclose(-k_B3 * R_B**3 / (3 * EPS0), -376.47, rtol=1e-4) and np.isclose(-4 * k_B3 * R_B, -4e-4)
print(f"B1: E_in(R) = {-k_B1 * R_B / EPS0:.4e} V/m, sigma_b = {k_B1 * R_B:.1e} C/m^2, rho_b = {-3 * k_B1:.1e} C/m^3.  "
      f"B3: V(0) = {-k_B3 * R_B**3 / (3 * EPS0):.2f} V, E_in(R) = {-k_B3 * R_B**2 / EPS0:.4e} V/m, rho_b(R) = {-4 * k_B3 * R_B:.1e} C/m^3")

# %% [markdown]
# ## Part 2: The uniformly polarized sphere, and what "the field inside" means
#
# Example 4.1 / (3.12): sigma_b = P cos(theta) produces E = -P/(3 eps0) inside and a pure
# dipole p = (4/3) pi R^3 P outside. Route 1 sums Coulomb fields of surface patches.
# Route 2 superposes two uniformly charged balls offset by d with rho d = P (section 3.6).
# Route 3 is a cubic lattice of point dipoles: the microscopic field at a lattice site
# is zero (depolarizing field -P/3eps0 plus Lorentz local field +P/3eps0), while the
# field outside is the dipole field. Then the averaging lemma (3.13) of section 3.7.
# Units in this part: eps0 = 1, so -P/(3 eps0) = -1/3 for P = 1.


# %%
def coulomb_field(points, src_pos, src_q):
    """E(points) = (1/4pi) sum_i q_i (r - r_i)/|r - r_i|^3 with eps0 = 1 (EM-2 (3.2))."""
    dr = points[:, None, :] - src_pos[None, :, :]
    r3 = np.linalg.norm(dr, axis=2)**3
    return (src_q[None, :, None] * dr / r3[:, :, None]).sum(axis=1) / (4 * np.pi)


def dipole_field(points, dip_pos, p):
    """E of point dipoles p (3,) at dip_pos: (1/4pi)[3(p.r_hat)r_hat - p]/r^3, EM-2 (4.3), eps0 = 1."""
    dr = points[:, None, :] - dip_pos[None, :, :]
    r = np.linalg.norm(dr, axis=2)
    rhat = dr / r[:, :, None]
    pr = rhat @ p
    return ((3 * pr[:, :, None] * rhat - p) / r[:, :, None]**3).sum(axis=1) / (4 * np.pi)


def polarized_sphere_patches(P=1.0, R=1.0, n_th=48, n_ph=96):
    """Surface charge sigma_b = P cos(theta) on a sphere, as patches (positions, charges)."""
    mu, w = leggauss(n_th)                         # Gauss-Legendre in cos(theta)
    phi = (np.arange(n_ph) + 0.5) * 2 * np.pi / n_ph
    MU, PH = np.meshgrid(mu, phi, indexing="ij")
    W = np.broadcast_to(w[:, None], MU.shape)
    sin = np.sqrt(1 - MU**2)
    pos = R * np.stack([sin * np.cos(PH), sin * np.sin(PH), MU], axis=-1).reshape(-1, 3)
    chg = (P * MU * R**2 * W * 2 * np.pi / n_ph).reshape(-1)   # sigma_b da
    return pos, chg


pos_s, chg_s = polarized_sphere_patches()
assert abs(chg_s.sum()) < 1e-12                                     # total bound charge zero (3.10)
assert np.isclose((chg_s * pos_s[:, 2]).sum(), 4 * np.pi / 3)       # dipole moment = (4/3) pi R^3 P
# Interior points (r <= 0.8 R): uniform field -P/3.
u = rng.standard_normal((150, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
pts_in = u * (0.8 * rng.random(150)**(1 / 3))[:, None]
E_in = coulomb_field(pts_in, pos_s, chg_s)
err_in = np.abs(E_in - np.array([0, 0, -1 / 3])).max()
assert err_in < 2e-3, err_in
# Exterior points (1.5 R to 4 R): pure dipole with p = (4/3) pi R^3 P z_hat.
pts_out = u * (1.5 + 2.5 * rng.random(150))[:, None]
E_out = coulomb_field(pts_out, pos_s, chg_s)
E_dip = dipole_field(pts_out, np.zeros((1, 3)), np.array([0, 0, 4 * np.pi / 3]))
err_out = (np.linalg.norm(E_out - E_dip, axis=1) / np.linalg.norm(E_dip, axis=1)).max()
assert err_out < 2e-3, err_out
print(f"Patch sum: interior field uniform -P/3 to {err_in:.1e}; exterior is a dipole to {err_out:.1e} relative.")


def ball_field(points, centre, rho, R=1.0):
    """Field of a uniformly charged ball (EM-2 Example 4.1), eps0 = 1: rho r/3 inside, rho R^3 r_hat/(3 r^2) outside."""
    dr = points - centre
    r = np.linalg.norm(dr, axis=1)
    inside = r < R
    f = np.where(inside, rho / 3, rho * R**3 / (3 * r**3))
    return f[:, None] * dr


# Route 2: two offset balls, rho d = P = 1. Inside the overlap the field is exactly -P/3 (no discretization).
d_off = 1e-3
E_two = ball_field(pts_in, np.array([0, 0, d_off / 2]), 1 / d_off) + ball_field(pts_in, np.array([0, 0, -d_off / 2]), -1 / d_off)
assert np.allclose(E_two, [0, 0, -1 / 3], atol=1e-9)
E_two_out = ball_field(pts_out, np.array([0, 0, d_off / 2]), 1 / d_off) + ball_field(pts_out, np.array([0, 0, -d_off / 2]), -1 / d_off)
assert np.allclose(E_two_out, E_dip, rtol=1e-5, atol=1e-8)       # exact up to O(d^2/r^2)
print("Two offset balls: interior field -P/3 to 1e-9; exterior dipole to 1e-5.")

# Route 3: cubic lattice of dipoles (spacing 1, p = z_hat, so P = 1) filling a sphere of radius R_lat.
R_lat = 8.0
g = np.arange(-9, 10)
G = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3).astype(float)
sites = G[np.linalg.norm(G, axis=1) <= R_lat]
p_hat = np.array([0, 0, 1.0])
probe = sites[np.linalg.norm(sites, axis=1) <= 4.0]            # deep interior lattice sites
E_site = np.array([dipole_field(s[None, :], sites[~np.all(sites == s, axis=1)], p_hat)[0] for s in probe])
micro = np.abs(E_site).max()
assert micro < 0.03, micro                                       # microscopic field at a site: ~0, not -1/3
far = u * (2 * R_lat + R_lat * rng.random(150))[:, None]        # 2R to 3R
E_far = dipole_field(far, sites, p_hat)
E_far_dip = dipole_field(far, np.zeros((1, 3)), len(sites) * p_hat)
err_far = (np.linalg.norm(E_far - E_far_dip, axis=1) / np.linalg.norm(E_far_dip, axis=1)).max()
assert err_far < 0.02, err_far
print(f"Lattice of {len(sites)} dipoles: |E| at interior lattice sites <= {micro:.3f} (macroscopic value would be 0.333); "
      f"exterior field is the dipole of the total moment to {err_far:.1e}.")

# Averaging lemma (3.13): the average over a ball of radius R_s of a point charge's field is
# -q r'/(4 pi R_s^3) for r' inside, and the field at the centre for r' outside. Computed as
# -(1/Vol) closed_int V n_hat da, which has no singularity, by Gauss-Legendre on the sphere.
mu, w = leggauss(80); phi = (np.arange(160) + 0.5) * 2 * np.pi / 160
MU, PH = np.meshgrid(mu, phi, indexing="ij"); W = np.broadcast_to(w[:, None], MU.shape) * 2 * np.pi / 160
nhat = np.stack([np.sqrt(1 - MU**2) * np.cos(PH), np.sqrt(1 - MU**2) * np.sin(PH), MU], -1).reshape(-1, 3)
for rp, expect in [(np.array([0.2, -0.1, 0.5]), -np.array([0.2, -0.1, 0.5]) / (4 * np.pi)),
                   (np.array([0.0, 1.0, 2.0]), -np.array([0.0, 1.0, 2.0]) / (4 * np.pi * 5**1.5))]:
    Vsurf = 1 / (4 * np.pi * np.linalg.norm(nhat - rp, axis=1))        # q = 1, R_s = 1
    E_avg = -(Vsurf[:, None] * nhat * W.reshape(-1, 1)).sum(axis=0) / (4 * np.pi / 3)
    assert np.allclose(E_avg, expect, atol=1e-8), (E_avg, expect)
print("Averaging lemma (3.13) verified for a charge inside and a charge outside the sphere.")

# Figure: macroscopic vs microscopic field along a line through the sphere.
zz = np.linspace(-3, 3, 241)
line = np.stack([0.37 * np.ones_like(zz), 0.21 * np.ones_like(zz), zz], -1)
E_line = coulomb_field(line, pos_s, chg_s)[:, 2]
site_r = np.linalg.norm(probe, axis=1)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(zz, E_line, label="patch sum, $E_z$ along a line")
ax1.axhline(-1 / 3, color="k", ls="--", lw=.8, label="$-P/3\\varepsilon_0$")
ax1.axvspan(-np.sqrt(1 - 0.37**2 - 0.21**2), np.sqrt(1 - 0.37**2 - 0.21**2), alpha=.1, color="C1", label="inside the sphere")
ax1.set_xlabel("z / R"); ax1.set_ylabel("$E_z$ (units $P/\\varepsilon_0$)"); ax1.legend(fontsize=8)
ax1.set_title("Macroscopic field of $\\sigma_b = P\\cos\\theta$")
ax2.plot(site_r / R_lat, np.linalg.norm(E_site, axis=1), "o", ms=3, label="|E| at lattice sites (microscopic)")
ax2.axhline(1 / 3, color="k", ls="--", lw=.8, label="$|{-P/3\\varepsilon_0}|$ (macroscopic)")
ax2.set_xlabel("r / R"); ax2.set_ylabel("$|E|$ (units $P/\\varepsilon_0$)"); ax2.legend(fontsize=8)
ax2.set_title("Lattice of dipoles: the site field is not the average")
fig.tight_layout(); fig.savefig(FIG / "polarized_sphere.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 3: A finite-difference solver for div(eps grad V) = -rho_f, and the half-filled capacitor
#
# (3.14) in the form div(eps_r grad V) = -rho_f/eps0 with E = -grad V. Nodes carry V; faces
# between neighbouring nodes carry eps_r (harmonic mean is what a sharp interface between two
# node values would give; here eps is sampled at face midpoints, so each face lies in one medium
# and the normal component of D = eps_r E is continuous across every face by construction).
# Dirichlet nodes are fixed; other boundary nodes are zero-flux.
#
# **Build it yourself first** (notes section 7): write `solve_poisson2d` and `node_charges` before reading on.


# %%
def weighted_faces(eps_xf, eps_yf):
    """Multiply face eps by the dual-cell face length in units of h: 1 inside, 1/2 along the grid boundary.

    With these weights, sum_f eps_f (V_nb - V_i) is exactly closed_int D.da/eps0 over node i's dual cell,
    including the half cells of boundary nodes (so zero-flux walls are the mirror-image discretization)."""
    wx = np.ones_like(eps_xf); wx[0, :] = wx[-1, :] = 0.5
    wy = np.ones_like(eps_yf); wy[:, 0] = wy[:, -1] = 0.5
    return eps_xf * wx, eps_yf * wy


def solve_poisson2d(eps_xf, eps_yf, rho_e0, h, fixed, fixed_val):
    """Solve div(eps_r grad V) = -rho_f/eps0 on an (ny, nx) node grid with spacing h.

    eps_xf (ny, nx-1): eps_r on faces between x-neighbours; eps_yf (ny-1, nx): between y-neighbours.
    Discrete form at node i: sum_faces eps_f (V_nb - V_i) = -rho_i h^2 / eps0, the flux balance
    closed_int D.da = Q_f (3.15) on the dual cell. Returns V with shape (ny, nx)."""
    ny, nx = fixed.shape
    N = ny * nx
    eps_xf, eps_yf = weighted_faces(eps_xf, eps_yf)
    idx = np.arange(N).reshape(ny, nx)
    diag = np.zeros((ny, nx)); rhs = -rho_e0 * h**2
    rhs[0, :] *= 0.5; rhs[-1, :] *= 0.5; rhs[:, 0] *= 0.5; rhs[:, -1] *= 0.5     # half cells hold half the charge
    rows, cols, vals = [], [], []
    S = slice
    for sl_n, sl_b, ef in [((S(None), S(None, -1)), (S(None), S(1, None)), eps_xf),
                           ((S(None), S(1, None)), (S(None), S(None, -1)), eps_xf),
                           ((S(None, -1), S(None)), (S(1, None), S(None)), eps_yf),
                           ((S(1, None), S(None)), (S(None, -1), S(None)), eps_yf)]:
        diag[sl_n] -= ef
        nb_fixed = fixed[sl_b]
        rhs[sl_n] -= np.where(nb_fixed, ef * fixed_val[sl_b], 0.0)        # known neighbours to the right side
        m = (~nb_fixed) & (~fixed[sl_n])
        rows.append(idx[sl_n][m]); cols.append(idx[sl_b][m]); vals.append(ef[m])
    diag[fixed] = 1.0; rhs[fixed] = fixed_val[fixed]
    rows.append(idx.ravel()); cols.append(idx.ravel()); vals.append(diag.ravel())
    A = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return spsolve(A, rhs.ravel()).reshape(ny, nx)


def node_charges(V, eps_xf, eps_yf):
    """Discrete closed_int D.da around each node in units of eps0 (per unit depth): Q_i = -sum_f eps_f (V_nb - V_i)."""
    eps_xf, eps_yf = weighted_faces(eps_xf, eps_yf)
    Q = np.zeros_like(V)
    fx = eps_xf * (V[:, 1:] - V[:, :-1]); fy = eps_yf * (V[1:, :] - V[:-1, :])
    Q[:, :-1] -= fx; Q[:, 1:] += fx; Q[:-1, :] -= fy; Q[1:, :] += fy
    return Q


def field_energy(V, eps_xf, eps_yf):
    """W/eps0 per unit depth = 1/2 sum_faces eps_f (dV)^2: the discrete 1/2 int D.E d tau of (3.23)."""
    eps_xf, eps_yf = weighted_faces(eps_xf, eps_yf)
    return 0.5 * (np.sum(eps_xf * (V[:, 1:] - V[:, :-1])**2) + np.sum(eps_yf * (V[1:, :] - V[:-1, :])**2))


def face_eps(eps_fn, xs, ys):
    """Sample eps_r(x, y) at the midpoints of x-faces and y-faces: use when every interface lies on a grid line."""
    xm, ym = 0.5 * (xs[1:] + xs[:-1]), 0.5 * (ys[1:] + ys[:-1])
    XX, YY = np.meshgrid(xm, ys); e_x = eps_fn(XX, YY)
    XX, YY = np.meshgrid(xs, ym); e_y = eps_fn(XX, YY)
    return e_x, e_y


def harmonic_faces(eps_nodes):
    """Face eps_r from node eps_r by the harmonic mean 2 e1 e2/(e1 + e2): the two half-cells a face joins are in
    series for the normal flux, so this is exact for a planar interface midway between the nodes (section 7)."""
    en = eps_nodes
    e_x = 2 * en[:, :-1] * en[:, 1:] / (en[:, :-1] + en[:, 1:])
    e_y = 2 * en[:-1, :] * en[1:, :] / (en[:-1, :] + en[1:, :])
    return e_x, e_y


# Example 4.3: plates at y = 0 (V = V0) and y = d (V = 0), width w, zero-flux side walls (no fringing).
# 40 x 40 nodes, so the interfaces at d/2 and w/2 fall midway between nodes (19.5 h) and the harmonic mean is exact.
d_cap, w_cap, V0, er = 1.0, 1.0, 1.0, 4.0
nc = 40
xs = np.linspace(0, w_cap, nc); ys = np.linspace(0, d_cap, nc)
h = xs[1] - xs[0]
XX, YY = np.meshgrid(xs, ys)
fixed = np.zeros((nc, nc), bool); fixed[0, :] = fixed[-1, :] = True
fixed_val = np.zeros((nc, nc)); fixed_val[0, :] = V0                   # bottom plate at V0, top grounded: E points +y
C0 = w_cap / d_cap                                                     # eps0 w/d per unit depth, EM-2 (3.18a)

# (a) series: slab of thickness d/2 against the bottom plate.
ex, ey = harmonic_faces(np.where(YY < d_cap / 2, er, 1.0))
V_ser = solve_poisson2d(ex, ey, np.zeros((nc, nc)), h, fixed, fixed_val)
Ey_face = -(V_ser[1:, :] - V_ser[:-1, :]) / h                        # E_y on y-faces
Dy_face = ey * Ey_face
assert np.allclose(Dy_face, 2 * er / (er + 1) * V0 / d_cap)          # D = sigma_f = 1.6 V0/d on every face (3.16)
assert np.allclose(Ey_face[:19], 2 / (er + 1) * V0 / d_cap) and np.allclose(Ey_face[20:], 2 * er / (er + 1) * V0 / d_cap)
assert np.allclose(ey[19], 2 * er / (er + 1)) and np.allclose(Ey_face[19], 1.0 * V0 / d_cap)   # the interface face: harmonic eps, mean E
C_ser = node_charges(V_ser, ex, ey)[0, :].sum() / V0                 # charge on the positive (bottom) plate
assert np.isclose(C_ser, 2 * er / (er + 1) * C0)                     # C = 1.6 C0
assert np.isclose(field_energy(V_ser, ex, ey), 0.5 * C_ser * V0**2)  # 1/2 int D.E = 1/2 C V^2 (discrete Green identity)
W_E0 = field_energy(V_ser, np.ones_like(ex), np.ones_like(ey))       # (eps0/2) int E^2 on the same grid
assert np.isclose(field_energy(V_ser, ex, ey) - W_E0, field_energy(V_ser, ex - 1, ey - 1))   # difference is 1/2 int P.E (3.24)
print(f"Series half-fill: C/C0 = {C_ser / C0:.4f} (theory {2 * er / (er + 1):.4f}); E_slab/E_gap = {Ey_face[0, 0] / Ey_face[-1, 0]:.4f}; "
      f"W = {field_energy(V_ser, ex, ey):.4f}, (eps0/2) int E^2 = {W_E0:.4f}, difference 1/2 int P.E = {field_energy(V_ser, ex, ey) - W_E0:.4f}")
# Exercise B6: the analytic split. Per unit area, W = D^2 d (er + 1)/(4 eps0 er) and 1/2 int P.E = D^2 d (er - 1)/(4 eps0 er^2),
# so the "spring" fraction is (er - 1)/(er (er + 1)) = 0.15 for er = 4. The grid's 0.1246/0.8 = 0.156 differs only because the
# single interface face carries the harmonic-mean eps and the mean E, which is exact for D.E but not for E^2.
D_ser_u = 2 * er / (er + 1)                                                          # units eps0 V0/d
W_an = D_ser_u**2 * (er + 1) / (4 * er); W_PE_an = D_ser_u**2 * (er - 1) / (4 * er**2)
assert np.isclose(W_an, field_energy(V_ser, ex, ey)) and np.isclose(W_PE_an / W_an, (er - 1) / (er * (er + 1))) and np.isclose(W_PE_an / W_an, 0.15)
assert np.isclose(W_an - W_PE_an, 0.5 * (0.5 * (D_ser_u / er)**2 + 0.5 * D_ser_u**2))    # (eps0/2) int E^2 exactly: 0.68
print(f"B6: W = {W_an:.4f}, (eps0/2) int E^2 = {W_an - W_PE_an:.4f}, 1/2 int P.E = {W_PE_an:.4f}, fraction {W_PE_an / W_an:.4f}")
# In units of D^2 d/eps0 (Exercise B6): W = (er + 1)/(4 er), 1/2 int P.E = (er - 1)/(4 er^2), (eps0/2) int E^2 = W - 1/2 int P.E.
assert np.isclose((er + 1) / (4 * er), 0.3125) and np.isclose((er - 1) / (4 * er**2), 0.046875) and np.isclose(0.3125 - 0.046875, 0.265625)

# (b) parallel: slab filling the left half of the width, full gap. V = V0 (1 - y/d) exactly.
ex2, ey2 = harmonic_faces(np.where(XX < w_cap / 2, er, 1.0))
V_par = solve_poisson2d(ex2, ey2, np.zeros((nc, nc)), h, fixed, fixed_val)
assert np.allclose(V_par, V0 * (1 - ys[:, None] / d_cap) * np.ones((1, nc)), atol=1e-10)
C_par = node_charges(V_par, ex2, ey2)[0, :].sum() / V0
assert np.isclose(C_par, (er + 1) / 2 * C0)                          # C = 2.5 C0
print(f"Parallel half-fill: C/C0 = {C_par / C0:.4f} (theory {(er + 1) / 2:.4f}); V is exactly linear.")

fig, ax = plt.subplots(figsize=(6, 4))
ym = 0.5 * (ys[1:] + ys[:-1])
ax.step(ym, Ey_face[:, nc // 2], where="mid", label="$E_y$ (units $V_0/d$)")
ax.plot(ym, Dy_face[:, nc // 2], "k--", label="$D_y$ (units $\\varepsilon_0V_0/d$)")
ax.axvline(d_cap / 2, color="gray", lw=.8); ax.text(0.26, 1.0, "$\\varepsilon_r = 4$", ha="center"); ax.text(0.76, 1.0, "vacuum", ha="center")
ax.set_xlabel("y / d"); ax.legend(); ax.set_title("Half-filled capacitor: $D_\\perp$ continuous, $E$ jumps")
fig.tight_layout(); fig.savefig(FIG / "slab_capacitor_DE.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 4: A dielectric cylinder (and a cavity) in a uniform field; the sphere by SymPy
#
# 2D analogue of (3.22): a cylinder of eps_in in a medium eps_out, uniform E0, has the uniform interior
# field E_in = 2 eps_out E0/(eps_in + eps_out) (Exercise B7). The box is 8R wide with V = -E0 x on its
# boundary, so the dipole term of the exterior solution is truncated at the 1% level.


# %%
def cylinder_in_field(eps_in, eps_out=1.0, Rc=1.0, L=8.0, n=241, E0=1.0, exact_bc=False):
    """Solve the cylinder-in-uniform-field problem on [-L, L]^2; return (E_x at centre, V, xs, eps_xf).

    The box boundary carries V = -E0 x (a uniform field, which truncates the exterior dipole term) or, with
    exact_bc=True, the exact exterior potential -E0 x + B x/s^2, B = Rc^2 E0 (eps_in - eps_out)/(eps_in + eps_out)
    (no truncation error, so what is left is the staircasing of the circle; Exercise B7(d))."""
    xs = np.linspace(-L, L, n); h = xs[1] - xs[0]
    ex, ey = face_eps(lambda X, Y: np.where(X**2 + Y**2 < Rc**2, eps_in, eps_out), xs, xs)
    fixed = np.zeros((n, n), bool); fixed[0, :] = fixed[-1, :] = fixed[:, 0] = fixed[:, -1] = True
    XG, YG = np.meshgrid(xs, xs)
    fixed_val = -E0 * XG
    if exact_bc:
        B_dip = Rc**2 * E0 * (eps_in - eps_out) / (eps_in + eps_out)
        fixed_val = fixed_val + B_dip * XG / np.where(XG**2 + YG**2 > 0, XG**2 + YG**2, 1.0)
    V = solve_poisson2d(ex, ey, np.zeros((n, n)), h, fixed, fixed_val)
    c = n // 2
    return -(V[c, c + 1] - V[c, c - 1]) / (2 * h), V, xs, ex


results = {}
for eps_in, eps_out in [(4.0, 1.0), (2.0, 1.0), (8.0, 1.0), (80.0, 1.0), (1.0, 4.0)]:
    E_c, V_cyl, xs_c, ex_c = cylinder_in_field(eps_in, eps_out)
    E_th = 2 * eps_out / (eps_in + eps_out)
    # Two error sources, separated in the budget below: the box (uniform-field data truncate the dipole term)
    # and the staircased circle. Together about 1-2% at h = R/15, L = 8R.
    assert abs(E_c - E_th) < 0.02, (eps_in, eps_out, E_c, E_th)
    results[(eps_in, eps_out)] = (E_c, E_th)
    print(f"cylinder eps_in = {eps_in:>4}, eps_out = {eps_out}: E_in/E0 = {E_c:.4f}  (theory {E_th:.4f}, error {E_c - E_th:+.4f})")
    if (eps_in, eps_out) == (4.0, 1.0):
        V4, xs4 = V_cyl, xs_c
# For eps_r -> infinity the cylinder becomes a conductor: E_in -> 0 (EM-3's conducting cylinder in a uniform field).
assert results[(80.0, 1.0)][0] < results[(8.0, 1.0)][0] < results[(4.0, 1.0)][0] < results[(2.0, 1.0)][0]
# Convergence check on the worst case (the cavity, eps_in = 1 in eps_out = 4): the error shrinks when h does.
E_cav_fine = cylinder_in_field(1.0, 4.0, n=321)[0]
assert abs(E_cav_fine - 1.6) < abs(results[(1.0, 4.0)][0] - 1.6)
print(f"cavity error {results[(1.0, 4.0)][0] - 1.6:+.4f} at n = 241 -> {E_cav_fine - 1.6:+.4f} at n = 321")

# Error budget for eps_in = 4 (Exercise B7(d)). Same h = R/15 and three box sizes give the *same* staircase, so the
# difference between uniform-field and exact boundary data isolates the box truncation, which scales as (R/L)^2.
budget = {}
for L_, n_ in [(4.0, 121), (8.0, 241), (16.0, 481)]:
    budget[L_] = (cylinder_in_field(4.0, 1.0, L=L_, n=n_)[0], cylinder_in_field(4.0, 1.0, L=L_, n=n_, exact_bc=True)[0])
trunc = {L_: (u_ - x_) * L_**2 for L_, (u_, x_) in budget.items()}
assert max(trunc.values()) / min(trunc.values()) < 1.05                       # truncation error ~ 0.21 (R/L)^2 E0
assert max(x_ for _, x_ in budget.values()) - min(x_ for _, x_ in budget.values()) < 5e-4   # staircase error independent of L
print("Cylinder error budget (eps_in = 4, h = R/15): " + "; ".join(
    f"L = {L_:g}R: uniform-BC {u_:.5f}, exact-BC {x_:.5f}, truncation*L^2 = {trunc[L_]:.4f}" for L_, (u_, x_) in budget.items()))
# Staircase error alone (exact boundary data): the sign flips as the grid changes; only the envelope shrinks, roughly like h.
stair = {n_: cylinder_in_field(4.0, 1.0, n=n_, exact_bc=True)[0] - 0.4 for n_ in (61, 121, 241, 481)}
assert stair[61] > 0.01 > max(abs(stair[241]), abs(stair[481])) and stair[121] < 0 < stair[241]
assert abs(stair[481]) < abs(stair[61]) / 5
print("Staircase error E_in/E0 - 0.4 (exact BC, L = 8R): " + ", ".join(f"n = {n_}: {v:+.4f}" for n_, v in stair.items()))

fig, ax = plt.subplots(figsize=(5.5, 5))
XX, YY = np.meshgrid(xs4, xs4)
Ex4, Ey4 = -np.gradient(V4, xs4[1] - xs4[0], axis=1), -np.gradient(V4, xs4[1] - xs4[0], axis=0)
ax.contour(XX, YY, V4, levels=np.linspace(-3, 3, 25), colors="C0", linewidths=.7)
ax.streamplot(XX, YY, Ex4, Ey4, density=0.9, color="gray", linewidth=.6, arrowsize=.7)
ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, color="C3", lw=1.5))
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_aspect("equal")
ax.set_title(f"Dielectric cylinder $\\varepsilon_r=4$ in $E_0$: $E_{{in}}/E_0 = {results[(4.0, 1.0)][0]:.3f}$ (theory 0.4)")
fig.tight_layout(); fig.savefig(FIG / "cylinder_uniform_field.png", dpi=120); plt.close(fig)

# The sphere (3.22) by separation of variables (section 3.11, route 2), symbolically.
A_, B_, E0s, Rs, ers = sp.symbols("A B E_0 R epsilon_r", positive=True)
sol = sp.solve([sp.Eq(-A_ * Rs, -E0s * Rs + B_ / Rs**2),                 # V continuous at r = R
                sp.Eq(ers * A_, E0s + 2 * B_ / Rs**3)],                  # D_perp continuous: eps_r dV_in/dr = dV_out/dr
               [A_, B_])
assert sp.simplify(sol[A_] - 3 * E0s / (ers + 2)) == 0
assert sp.simplify(sol[B_] - Rs**3 * E0s * (ers - 1) / (ers + 2)) == 0
# 2D version used above: V_in = -A s cos(phi), V_out = -E0 s cos(phi) + B cos(phi)/s.
# Matching: V continuous, eps dV/ds continuous (Exercise B7(b)), general eps_in, eps_out.
ein, eout = sp.symbols("epsilon_in epsilon_out", positive=True)
sol2 = sp.solve([sp.Eq(-A_ * Rs, -E0s * Rs + B_ / Rs), sp.Eq(ein * A_, eout * (E0s + B_ / Rs**2))], [A_, B_])
assert sp.simplify(sol2[A_] - 2 * eout * E0s / (ein + eout)) == 0
assert sp.simplify(sol2[B_] - Rs**2 * E0s * (ein - eout) / (ein + eout)) == 0
# Example 4.4 numbers: E0 = 1e5 V/m, R = 2 cm, eps_r = 4.
E0n, Rn, ern = 1e5, 0.02, 4.0
E_in_n = 3 * E0n / (ern + 2); P_n = EPS0 * (ern - 1) * E_in_n; p_n = 4 / 3 * np.pi * Rn**3 * P_n
assert np.isclose(E_in_n, 5e4) and np.isclose(P_n, 1.5 * EPS0 * E0n) and np.isclose(p_n, 2 * np.pi * EPS0 * Rn**3 * E0n)
E_pole_out = E0n + 2 * p_n / (4 * np.pi * EPS0 * Rn**3)
assert np.isclose(E_pole_out, 2 * E0n) and np.isclose(E_pole_out - E_in_n, P_n / EPS0)      # jump = sigma_b/eps0 at the pole
assert np.isclose(ern * EPS0 * E_in_n, EPS0 * E_pole_out)                                     # D_perp continuous at the pole
print(f"Example 4.4: E_in = {E_in_n:.3e} V/m, P = {P_n:.4e} C/m^2, p = {p_n:.4e} C m, E just outside the pole = {E_pole_out:.3e} V/m")

# %% [markdown]
# ## Part 5: Energy. The metal sphere in a dielectric shell (Example 4.2), two ways
#
# (3.23) W = 1/2 int D.E d tau for linear media equals 1/2 Q V for a capacitor. Example 4.2:
# Q = 1 nC on a metal sphere a = 1 cm inside a shell of eps_r = 4 out to b = 2 cm.

# %%
Qn, an, bn, ern = 1e-9, 0.01, 0.02, 4.0
chi = ern - 1
E_a_in = K_E * Qn / (ern * an**2); E_b_in = K_E * Qn / (ern * bn**2); E_b_out = K_E * Qn / bn**2
sig_b_a = -chi * Qn / (4 * np.pi * ern * an**2); sig_b_b = chi * Qn / (4 * np.pi * ern * bn**2)
Qb_a = sig_b_a * 4 * np.pi * an**2; Qb_b = sig_b_b * 4 * np.pi * bn**2
assert np.isclose(Qb_a + Qb_b, 0) and np.isclose(Qb_b, 0.75e-9)
V_a = K_E * Qn * (1 / (ern * an) + (1 - 1 / ern) / bn)
C_shell = Qn / V_a
# Gauss check with total charge: field just outside the metal from (Q + Q_b(a))/(4 pi eps0 a^2)
assert np.isclose(E_a_in, K_E * (Qn + Qb_a) / an**2) and np.isclose(E_b_out, K_E * (Qn + Qb_a + Qb_b) / bn**2)
# Energy: 1/2 int D.E = int_a^b (1/2) eps E^2 4 pi r^2 dr + int_b^inf (1/2) eps0 E^2 4 pi r^2 dr, done by quadrature.
rr1 = np.linspace(an, bn, 20001); rr2 = np.linspace(bn, 50.0, 400001)
W_int = np.trapezoid(0.5 * ern * EPS0 * (K_E * Qn / (ern * rr1**2))**2 * 4 * np.pi * rr1**2, rr1) \
      + np.trapezoid(0.5 * EPS0 * (K_E * Qn / rr2**2)**2 * 4 * np.pi * rr2**2, rr2) + 0.5 * EPS0 * (K_E * Qn)**2 * 4 * np.pi / 50.0
assert np.isclose(W_int, 0.5 * Qn * V_a, rtol=1e-6)
print(f"Example 4.2: E(a+) = {E_a_in:.1f} V/m, E(b-) = {E_b_in:.1f}, E(b+) = {E_b_out:.1f} V/m; sigma_b(a) = {sig_b_a:.4e}, "
      f"sigma_b(b) = {sig_b_b:.4e} C/m^2; Q_b = {Qb_a:.3e} + {Qb_b:.3e}; V(a) = {V_a:.2f} V (vacuum {K_E * Qn / an:.2f}); "
      f"C = {C_shell * 1e12:.4f} pF; W = {0.5 * Qn * V_a:.4e} J = 1/2 int D.E")

# Example 4.1 numbers: R = 1 cm, P = 1 uC/m^2.
Pn, Rn1 = 1e-6, 0.01
assert np.isclose(-Pn / (3 * EPS0), -3.7647e4, rtol=1e-4) and np.isclose(4 / 3 * np.pi * Rn1**3 * Pn, 4.18879e-12, rtol=1e-5)
assert np.isclose(2 * Pn / (3 * EPS0), 7.5294e4, rtol=1e-4) and np.isclose(Pn * Rn1 / (3 * EPS0), 376.47, rtol=1e-4)
print(f"Example 4.1: E_in = {-Pn / (3 * EPS0):.4e} V/m, p = {4 / 3 * np.pi * Rn1**3 * Pn:.4e} C m, E(pole, outside) = {2 * Pn / (3 * EPS0):.4e} V/m, V(pole) = {Pn * Rn1 / (3 * EPS0):.2f} V")

# Example 4.3 numbers: A = 100 cm^2, d = 2 mm, eps_r = 4, V = 100 V.
C0n = EPS0 * 0.01 / 0.002
assert np.isclose(C0n, 44.27e-12, rtol=1e-3) and np.isclose(1.6 * C0n, 70.83e-12, rtol=1e-3) and np.isclose(2.5 * C0n, 110.68e-12, rtol=1e-3)
D_ser = 1.6 * C0n * 100 / 0.01
assert np.isclose(D_ser, 7.0834e-7, rtol=1e-4) and np.isclose(D_ser / EPS0, 8.0e4) and np.isclose(0.75 * D_ser, 5.3125e-7, rtol=1e-4)
print(f"Example 4.3: C0 = {C0n * 1e12:.2f} pF, series {1.6 * C0n * 1e12:.2f} pF, parallel {2.5 * C0n * 1e12:.2f} pF; "
      f"series at 100 V: D = {D_ser:.4e} C/m^2, E_gap = {D_ser / EPS0:.3e}, E_slab = {D_ser / (4 * EPS0):.3e} V/m, sigma_b = -/+{0.75 * D_ser:.4e}")

# %% [markdown]
# ## Part 6: Force on a partially inserted slab
#
# Section 3.13 / (3.25): F = 1/2 V^2 dC/dx at fixed V (battery included) and also at fixed Q, with
# dC/dx = eps0 chi_e w/d for a slab filling the gap. Here the plates are finite (length 8d) inside a
# large box, the slab enters from the left and fills the gap to depth x. C(x) from the plate charge
# should be linear in x with slope chi_e/d per unit depth (units eps0).


# %%
def slab_capacitor(x_ins, er=4.0, h=0.05, full=False):
    """Finite plates y = +-1/2 for 0 < x < 8 at V = +-1/2, slab eps_r filling |y| < 1/2 for x < x_ins.

    Returns (C, W)/eps0 per unit depth, plus (V, xs, ys) when full=True."""
    xs = np.round(np.arange(-5, 13 + h / 2, h), 10); ys = np.round(np.arange(-5, 5 + h / 2, h), 10)   # exact +-0.5 rows
    nx, ny = len(xs), len(ys)
    ex, ey = face_eps(lambda X, Y: np.where((np.abs(Y) < 0.5) & (X < x_ins), er, 1.0), xs, ys)
    fixed = np.zeros((ny, nx), bool); fixed[0, :] = fixed[-1, :] = fixed[:, 0] = fixed[:, -1] = True
    fixed_val = np.zeros((ny, nx))
    jt, jb = np.argmin(np.abs(ys - 0.5)), np.argmin(np.abs(ys + 0.5))
    on_plate = (xs >= -h / 4) & (xs <= 8 + h / 4)
    fixed[jt, on_plate] = fixed[jb, on_plate] = True
    fixed_val[jt, on_plate] = 0.5; fixed_val[jb, on_plate] = -0.5
    V = solve_poisson2d(ex, ey, np.zeros((ny, nx)), h, fixed, fixed_val)
    Q_top = node_charges(V, ex, ey)[jt, on_plate].sum()
    if full:
        return Q_top / 1.0, field_energy(V, ex, ey), V, xs, ys
    return Q_top / 1.0, field_energy(V, ex, ey)


x_list = np.array([2.0, 3.0, 4.0, 5.0, 6.0])
CW = np.array([slab_capacitor(xi) for xi in x_list])
C_x, W_x = CW[:, 0], CW[:, 1]
assert np.allclose(W_x, 0.5 * C_x * 1.0**2, rtol=1e-10)                      # W = 1/2 C V^2 on the grid
slope, icpt = np.polyfit(x_list, C_x, 1)
chi_e = 3.0
assert abs(slope - chi_e / 1.0) < 0.03 * chi_e, slope                         # dC/dx = eps0 chi_e w/d per unit w, d = 1
assert np.abs(C_x - (slope * x_list + icpt)).max() < 1e-3 * C_x.mean()        # linear in x
F_V = 0.5 * 1.0**2 * slope                                                    # fixed V, battery's work included
Q0 = C_x[2] * 1.0                                                             # fixed Q: hold the charge of the x = 4 state
W_Q = Q0**2 / (2 * C_x)
C_lo, C_hi = slab_capacitor(3.9)[0], slab_capacitor(4.1)[0]                  # central difference of W(x) at fixed Q, step 0.1 d
F_Q = -(Q0**2 / (2 * C_hi) - Q0**2 / (2 * C_lo)) / 0.2                       # -dW/dx at fixed Q
assert np.isclose(F_Q, F_V, rtol=3e-3)
print(f"Slab force: C(x) = {slope:.4f} x + {icpt:.3f} (units eps0 per unit depth; theory slope chi_e/d = {chi_e:.1f}); "
      f"F/eps0 V^2 = {F_V:.4f} (fixed V) vs {F_Q:.4f} (fixed Q, finite difference of W = Q^2/2C).")

# The fringing-field picture (section 3.13): force density on the induced dipoles f = (P.grad)E = (1/2) eps0 chi_e grad(E^2).
# Integrated over the slab, F_x = (1/2) chi_e [g(end) - g(start)] with g(x) = int E^2 dy across the gap, so the cumulative
# force up to x is (1/2) chi_e g(x): it is delivered where E^2 rises, at the plates' leading edge x = 0, not at the slab's end face.
_, _, V_k, xs_k, ys_k = slab_capacitor(4.0, full=True)
h_k = 0.05
Ey_c = -(V_k[1:, :-1] + V_k[1:, 1:] - V_k[:-1, :-1] - V_k[:-1, 1:]) / (2 * h_k)     # E at cell centres
Ex_c = -(V_k[:-1, 1:] + V_k[1:, 1:] - V_k[:-1, :-1] - V_k[1:, :-1]) / (2 * h_k)
xc_k, yc_k = xs_k[:-1] + h_k / 2, ys_k[:-1] + h_k / 2
g_k = h_k * ((Ex_c**2 + Ey_c**2)[np.abs(yc_k) < 0.5, :]).sum(axis=0)                # int E^2 dy across the gap, per column
col = xc_k < 4.0
F_kelvin = 0.5 * chi_e * (g_k[col][-1] - g_k[col][0])
cum = lambda xq: 0.5 * chi_e * g_k[np.searchsorted(xc_k, xq)] / F_kelvin
assert np.isclose(F_kelvin, F_V, rtol=0.02) and cum(-1.0) < 0.05 and cum(1.0) > 0.97
print(f"Fringing-field picture: (P.grad)E integrates to {F_kelvin:.4f} (energy method {F_V:.4f}); "
      f"cumulative fraction of the force at x = -1, 0, +1: {cum(-1.0):.3f}, {cum(0.0):.3f}, {cum(1.0):.3f} (plate edge at 0, slab end at 4).")
# Exercise B8 numbers: d = 1 cm, w = 10 cm, V = 1 kV, eps_r = 4.
F_B8 = EPS0 * 3 * 0.10 * 1e3**2 / (2 * 0.01)
assert np.isclose(F_B8, 1.3281e-4, rtol=1e-4)
print(f"Exercise B8: F = {F_B8:.4e} N (about the weight of {F_B8 / 9.81 * 1e6:.1f} mg)")
# B8 energy accounting, plates 10 cm long: slab pulled in from x1 = 2 cm to x2 = 6 cm, d = 1 cm, w = 10 cm, V = 1 kV, eps_r = 4.
w8, d8, V8, L8, x1, x2 = 0.10, 0.01, 1e3, 0.10, 0.02, 0.06
C8 = lambda xx: EPS0 * w8 / d8 * (L8 + 3.0 * xx)                          # C(x) = eps0 w (l + chi_e x)/d
dC8 = C8(x2) - C8(x1)
Wb, Wf, Wm = V8**2 * dC8, 0.5 * V8**2 * dC8, F_B8 * (x2 - x1)             # battery, field energy gain, mechanical work (fixed V)
assert np.isclose(dC8, 1.0625e-11, rtol=1e-4) and np.isclose(Wb, Wf + Wm) and np.isclose(Wm, 5.3125e-6, rtol=1e-4)
Q8 = C8(x1) * V8; W1, W2 = Q8**2 / (2 * C8(x1)), Q8**2 / (2 * C8(x2))   # fixed Q: disconnect the battery at x1
V2_8 = Q8 / C8(x2); F2_8 = 0.5 * V2_8**2 * EPS0 * 3 * w8 / d8
assert np.isclose(C8(x1), 14.167e-12, rtol=1e-3) and np.isclose(W1 - W2, 3.0357e-6, rtol=1e-3) and np.isclose(V2_8, 571.4, rtol=1e-3)
assert np.isclose(F2_8, 4.336e-5, rtol=1e-3) and np.isclose(0.5 * V8**2 * EPS0 * 3 * w8 / d8, F_B8)
print(f"B8: fixed V: dC = {dC8 * 1e12:.3f} pF, battery {Wb:.4e} J = field {Wf:.4e} + mechanical {Wm:.4e}.  "
      f"Fixed Q (disconnect at x = 2 cm): C1 = {C8(x1) * 1e12:.3f} pF, C2 = {C8(x2) * 1e12:.3f} pF, Q = {Q8:.4e} C, W1 = {W1:.4e} J, W2 = {W2:.4e} J, "
      f"V falls {V8:.0f} -> {V2_8:.1f} V, work done by the field = {W1 - W2:.4e} J, dC/dx = {EPS0 * 3 * w8 / d8:.4e} F/m, "
      f"F falls {F_B8:.4e} -> {F2_8:.4e} N")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(x_list, C_x, "o", label="C(x) from plate charge")
ax1.plot(x_list, slope * x_list + icpt, "k--", lw=.8, label=f"fit slope {slope:.3f} (theory 3)")
ax1.set_xlabel("insertion x / d"); ax1.set_ylabel("$C/\\varepsilon_0$ per unit width"); ax1.legend()
ax2.plot(x_list, W_x, "o-", label="fixed $V$: $W = \\frac{1}{2} CV^2$ rises")
ax2.plot(x_list, W_Q, "s-", label="fixed $Q$: $W = Q^2/2C$ falls")
ax2.set_xlabel("insertion x / d"); ax2.set_ylabel("$W/\\varepsilon_0V^2$"); ax2.legend()
ax2.set_title("Same force $\\frac{1}{2}V^2\\,dC/dx$ either way")
fig.tight_layout(); fig.savefig(FIG / "slab_force.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: Clausius-Mossotti, and the numbers in the exercises
#
# (3.27): (eps_r - 1)/(eps_r + 2) = N alpha/(3 eps0), equivalently chi_e = 3x/(1 - x) with x = N alpha/3eps0.
# Example 4.5: N2 at STP with alpha/4pi eps0 = 1.74 A^3 (tabulated; the crude model with a = 1.2 A gives 1.73).


# %%
def clausius_mossotti(N, alpha_over_4pieps0):
    """chi_e from (3.27); alpha_over_4pieps0 in m^3 (alpha/4pi eps0, the 'polarizability volume')."""
    xcm = 4 * np.pi * N * alpha_over_4pieps0 / 3
    return 3 * xcm / (1 - xcm)


xs_ = sp.symbols("x_cm", positive=True)
assert sp.simplify((1 + 3 * xs_ / (1 - xs_) - 1) / (1 + 3 * xs_ / (1 - xs_) + 2) - xs_) == 0   # the two forms of (3.27) agree
N_stp = 101325 / (K_B * 273.15)
assert np.isclose(N_stp, 2.6868e25, rtol=1e-4)
chi_N2 = clausius_mossotti(N_stp, 1.74 * A3)
chi_N2_crude = clausius_mossotti(N_stp, 1.2**3 * A3)
chi_tab = 0.00055 * 293.15 / 273.15                       # Griffiths Table 4.2: eps_r(N2, 20 C, 1 atm) = 1.00055, scaled to 0 C
assert np.isclose(chi_N2, 5.876e-4, rtol=1e-3) and np.isclose(chi_N2_crude, 5.834e-4, rtol=1e-3)
assert abs(chi_N2 - chi_tab) / chi_tab < 0.01 and abs(chi_N2 - 5.5e-4) / 5.5e-4 < 0.07
print(f"Example 4.5: N = {N_stp:.4e} /m^3; chi_e(N2) = {chi_N2:.4e} (alpha tabulated), {chi_N2_crude:.4e} (a = 1.2 A); "
      f"table 5.5e-4 at 20 C -> {chi_tab:.3e} at 0 C")
# C1: liquid argon (works) and water (fails for the static value, works for the optical one).
N_Ar = 1.40 / 39.948 * 6.02214e23 * 1e6; chi_Ar = clausius_mossotti(N_Ar, 1.64 * A3)
N_w = 1.00 / 18.015 * 6.02214e23 * 1e6; chi_w = clausius_mossotti(N_w, 1.45 * A3)
assert np.isclose(4 * np.pi / 3 * N_Ar * 1.64 * A3, 0.14495, rtol=1e-3) and np.isclose(4 * np.pi / 3 * N_w * 1.45 * A3, 0.2030, rtol=1e-3)
assert np.isclose(N_Ar, 2.11e28, rtol=5e-3) and np.isclose(1 + chi_Ar, 1.51, atol=0.01)
assert np.isclose(N_w, 3.34e28, rtol=5e-3) and np.isclose(1 + chi_w, 1.76, atol=0.01) and abs(1 + chi_w - 1.333**2) < 0.02
print(f"C1: liquid argon N = {N_Ar:.3e}, x = {4 * np.pi / 3 * N_Ar * 1.64 * A3:.4f}, eps_r(CM) = {1 + chi_Ar:.3f} (measured about 1.5); "
      f"water N = {N_w:.3e}, x = {4 * np.pi / 3 * N_w * 1.45 * A3:.4f}, eps_r(CM, electronic alpha) = {1 + chi_w:.3f} vs n^2 = {1.333**2:.3f} and static 80.1")
assert np.isclose(1 + clausius_mossotti(N_Ar, 1.64 * A3), 1.5086, atol=2e-3) and np.isclose(1 + 4 * np.pi * N_Ar * 1.64 * A3, 1.4349, atol=1e-3)
# C1(c)-(d): add Debye's orientational polarizability p^2/(3 k T) to water's alpha, then Onsager's repair (quoted, not derived):
#   (eps - n^2)(2 eps + n^2) / (eps (n^2 + 2)^2) = N p^2 / (9 eps0 k T).
T_w, p_w = 293.15, 6.17e-30
a_or = p_w**2 / (3 * K_B * T_w) / (4 * np.pi * EPS0)                          # orientational alpha/4pi eps0 in m^3
a_tot = 1.45 * A3 + a_or
x_deb = 4 * np.pi * N_w * a_tot / 3                                           # N alpha/3eps0 with the Debye term
eps_dilute_deb = 1 + 4 * np.pi * N_w * a_tot                                  # no local-field correction
eps_cm_deb = (1 + 2 * x_deb) / (1 - x_deb)                                    # Clausius-Mossotti with the Debye term
n2_w = 1.333**2
rhs_ons = N_w * p_w**2 / (9 * EPS0 * K_B * T_w)
eps_ons = brentq(lambda e: (e - n2_w) * (2 * e + n2_w) / (e * (n2_w + 2)**2) - rhs_ons, n2_w, 1e3)
assert np.isclose(a_or / A3, 28.18, rtol=1e-3) and np.isclose(x_deb, 4.149, rtol=1e-3) and x_deb > 1 and eps_cm_deb < 0
assert np.isclose(eps_dilute_deb, 13.4, atol=0.1) and np.isclose(rhs_ons, 3.946, rtol=1e-3) and np.isclose(eps_ons, 29.2, atol=0.3)
assert 2.5 < 80.1 / eps_ons < 3.0
lhs_80 = (80.1 - n2_w) * (2 * 80.1 + n2_w) / (80.1 * (n2_w + 2)**2); g_K = lhs_80 / rhs_ons       # p^2 -> g p^2 reaches 80.1
assert np.isclose(g_K, 2.81, atol=0.02) and np.isclose(lhs_80, 11.10, atol=0.01)
print(f"C1(c,d): alpha_or/4pi eps0 = {a_or / A3:.2f} A^3 (electronic 1.45); N alpha/3eps0 = {x_deb:.3f} > 1, so CM gives eps_r = {eps_cm_deb:.2f} (nonsense); "
      f"dilute+Debye {eps_dilute_deb:.1f}; Onsager (RHS {rhs_ons:.3f}) eps_r = {eps_ons:.1f}; 80.1/Onsager = {80.1 / eps_ons:.2f}; "
      f"left side at 80.1 = {lhs_80:.3f}, g needed = {g_K:.2f}")

Nalpha = np.linspace(0, 2.4, 300)                                           # N alpha / eps0
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(Nalpha, Nalpha, "k--", lw=.8, label="dilute: $\\chi_e = N\\alpha/\\varepsilon_0$")
ax.plot(Nalpha, 3 * (Nalpha / 3) / (1 - Nalpha / 3), label="Clausius-Mossotti (3.27)")
for lab, Nn, al in [("N$_2$ gas", N_stp, 1.74 * A3), ("liquid Ar", N_Ar, 1.64 * A3), ("water (electronic)", N_w, 1.45 * A3)]:
    v = 4 * np.pi * Nn * al; ax.plot(v, clausius_mossotti(Nn, al), "o"); ax.annotate(lab, (v, clausius_mossotti(Nn, al)), textcoords="offset points", xytext=(5, -10), fontsize=8)
ax.set_ylim(0, 8); ax.set_xlabel("$N\\alpha/\\varepsilon_0$"); ax.set_ylabel("$\\chi_e$"); ax.legend()
ax.set_title("Local field makes $\\chi_e$ superlinear in density")
fig.tight_layout(); fig.savefig(FIG / "clausius_mossotti.png", dpi=120); plt.close(fig)

# Exercise numbers (solutions.md).
a_H = (0.667) ** (1 / 3); p_H = 4 * np.pi * EPS0 * 0.667 * A3 * 1e6; d_H = p_H / E_CH         # A1
assert np.isclose(a_H, 0.8737, rtol=1e-3) and np.isclose(p_H, 7.421e-35, rtol=1e-3) and np.isclose(d_H, 4.632e-16, rtol=1e-3)
E_atom = K_E * E_CH / (a_H * 1e-10)**2                                                           # A1(c): restoring field at d = a
assert np.isclose(E_atom, 1.886e11, rtol=1e-3) and np.isclose(d_H / (a_H * 1e-10), 1e6 / E_atom) and np.isclose(1e6 / E_atom, 5.30e-6, rtol=2e-3)
N_w_tq = 6.17e-30 * 1e5; ratio_kT = 2 * N_w_tq / (K_B * 300); E_unit = K_B * 300 / (2 * 6.17e-30)  # A3 (water, p = 1.85 D)
assert np.isclose(N_w_tq, 6.17e-25) and np.isclose(ratio_kT, 2.979e-4, rtol=1e-3) and np.isclose(E_unit, 3.356e8, rtol=1e-3)
cos_avg = 6.17e-30 * 1e5 / (K_B * 300) / 3                                                       # Langevin <cos> = x/3, x = pE/kT
assert np.isclose(cos_avg, 4.97e-5, rtol=2e-3) and np.isclose(E_unit / 3e6, 112, rtol=5e-3)
C_A4 = 2.1 * EPS0 * 50e-4 / 1e-3; Q_A4 = C_A4 * 100; sig_A4 = Q_A4 / 50e-4; P_A4 = sig_A4 - EPS0 * 1e5     # A4
assert np.isclose(C_A4, 92.97e-12, rtol=1e-3) and np.isclose(sig_A4, 1.859e-6, rtol=1e-3) and np.isclose(P_A4, 9.74e-7, rtol=2e-3)
assert np.isclose(sig_A4 - P_A4, EPS0 * 1e5)                                                   # sigma_f + sigma_b = eps0 E
th2 = np.degrees(np.arctan(5 * np.tan(np.radians(30))))                                         # A5
assert np.isclose(th2, 70.89, rtol=1e-3) and np.isclose(np.sin(np.radians(30)) / np.sin(np.radians(th2)), 0.5292, rtol=1e-3)
sig_B4 = (2 / 3) * 1e-9 / (4 * np.pi * 0.05**2); E_B4 = K_E * 1e-9 / (3 * 0.02**2)               # B4
assert np.isclose(sig_B4, 2.122e-8, rtol=1e-3) and np.isclose(E_B4, 7490, rtol=1e-3)
Qc_B4 = -(3 - 1) / 3 * 1e-9                                                                      # bound charge at the point charge: -(chi_e/eps_r) q
assert np.isclose(Qc_B4 + sig_B4 * 4 * np.pi * 0.05**2, 0) and np.isclose(Qc_B4, -0.6667e-9, rtol=1e-3)
assert np.isclose(K_E * 1e-9 / 0.05**2, 3595.2, rtol=1e-3)                                       # E just outside (and unchanged): kq/R^2
E_Rin_B4 = K_E * 1e-9 / (3 * 0.05**2); jump_B4 = sig_B4 / EPS0                                    # E just inside R and the jump sigma_b/eps0
assert np.isclose(E_Rin_B4, 1198.4, rtol=1e-3) and np.isclose(jump_B4, 2397, rtol=1e-3) and np.isclose(K_E * 1e-9 / 0.05**2 - E_Rin_B4, jump_B4)
D_B5 = EPS0 * 200 / (1e-3 / 2 + 2e-3 / 5); E1, E2 = D_B5 / (2 * EPS0), D_B5 / (5 * EPS0)        # B5
assert np.isclose(D_B5, 1.9676e-6, rtol=1e-3) and np.isclose(E1 * 1e-3 + E2 * 2e-3, 200)
assert np.isclose(D_B5 * (1 / 5 - 1 / 2), EPS0 * (E2 - E1)) and np.isclose(D_B5 * 20e-4 / 200, 19.68e-12, rtol=1e-3)
P1_B5, P2_B5 = D_B5 * (1 - 1 / 2), D_B5 * (1 - 1 / 5)                                            # P in each layer; faces: -P1, +P1, -P2, +P2
assert np.isclose(P1_B5, 9.838e-7, rtol=1e-3) and np.isclose(P2_B5, 1.5741e-6, rtol=1e-3) and np.isclose(P1_B5 - P2_B5, D_B5 * (1 / 5 - 1 / 2))
assert np.isclose(D_B5 * 20e-4, 3.935e-9, rtol=1e-3)                                            # Q = sigma_f A = C V
C1_B5, C2_B5 = 2 * EPS0 * 20e-4 / 1e-3, 5 * EPS0 * 20e-4 / 2e-3                                 # series: 1/C = 1/C1 + 1/C2
assert np.isclose(C1_B5, 35.42e-12, rtol=1e-3) and np.isclose(C2_B5, 44.27e-12, rtol=1e-3) and np.isclose(C1_B5 * C2_B5 / (C1_B5 + C2_B5), 19.68e-12, rtol=1e-3)
E_A5 = 1.0                                                                                       # A5: |E2|/|E1| = sin(30)/sin(theta_2)
print(f"A1: a = {a_H:.4f} A, p = {p_H:.3e} C m, d = {d_H:.3e} m, E_atom = {E_atom:.3e} V/m, d/a = {1e6 / E_atom:.2e}.  "
      f"A3: N = {N_w_tq:.2e} N m, 2pE/kT = {ratio_kT:.2e}, <cos> = {cos_avg:.2e}, field for 2pE/kT = 1: {E_unit:.3e} V/m = {E_unit / 3e6:.0f} x air breakdown.  "
      f"A4: C = {C_A4 * 1e12:.2f} pF, Q = {Q_A4 * 1e9:.3f} nC, sigma_f = {sig_A4:.3e}, P = {P_A4:.3e} C/m^2, sigma_f + sigma_b = {(sig_A4 - P_A4):.3e}.  "
      f"A5: theta_2 = {th2:.2f} deg, |E2|/|E1| = {np.sin(np.radians(30)) / np.sin(np.radians(th2)):.4f}.")
print(f"B4: sigma_b(R) = {sig_B4:.3e} C/m^2, bound charge at the centre = {Qc_B4:.4e} C, E(2 cm) = {E_B4:.0f} V/m, E(R-) = {E_Rin_B4:.1f}, E(R+) = {K_E * 1e-9 / 0.05**2:.1f}, sigma_b/eps0 = {jump_B4:.0f} V/m.  B5: D = {D_B5:.4e} C/m^2, "
      f"E1 = {E1:.4e}, E2 = {E2:.4e} V/m, P1 = {P1_B5:.4e}, P2 = {P2_B5:.4e}, sigma_b(interface) = {D_B5 * (1 / 5 - 1 / 2):.3e} C/m^2, "
      f"C = {D_B5 * 20e-4 / 200 * 1e12:.2f} pF (C1 = {C1_B5 * 1e12:.2f}, C2 = {C2_B5 * 1e12:.2f} pF in series), Q = {D_B5 * 20e-4 * 1e9:.3f} nC")

# %% [markdown]
# ## Your turn
#
# 1. **Coaxial cable with two dielectrics.** Inner conductor radius a, dielectric eps_1 to radius b,
#    eps_2 to the outer conductor at c. Use Gauss's law for D (3.15) to get D(s), then E in each
#    layer, then V and C per unit length. Find sigma_b at s = b and check that the total bound charge
#    per length is zero. Then solve it with `solve_poisson2d` (annulus in a square box) and compare.
# 2. **Where Clausius-Mossotti fails.** Water has N = 3.34e28 /m^3 and electronic alpha/4pi eps0 = 1.45 A^3.
#    Compute eps_r from (3.27). It gives about 1.76, which is n^2 at optical frequencies, not the static 80.
#    Explain which assumption of section 3.14 breaks for a polar liquid, and what Debye's p^2/(3 k T) term adds
#    (Part 7 shows that adding it to (3.27) overshoots the polarization catastrophe; Exercise C1 goes on to Onsager).
# 3. **A piezoelectric bar.** Quartz has d_11 = 2.3 pC/N: a stress sigma along the x axis produces P_x = d_11 sigma.
#    For sigma = 1 MPa on a 1 cm x 1 cm face of a 1 mm thick plate with eps_r = 4.5, compute sigma_b, the surface
#    charge Q, and the open-circuit voltage across the plate (treat it as a capacitor with D = 0 inside).
#
# Reference solutions follow. Write yours first.

# %%
a_c, b_c, c_c, e1, e2 = 1e-3, 2e-3, 4e-3, 2.3, 4.0
C_coax = 2 * np.pi * EPS0 / (np.log(b_c / a_c) / e1 + np.log(c_c / b_c) / e2)       # F/m
assert np.isclose(C_coax, 1.172e-10, rtol=2e-3)
lam = 1e-9                                                                          # C/m on the inner conductor
P1_b = (1 - 1 / e1) * lam / (2 * np.pi * b_c); P2_b = (1 - 1 / e2) * lam / (2 * np.pi * b_c)   # P = (1 - 1/eps_r) D in each layer
sig_b_b = P1_b - P2_b                      # outer face of layer 1 (n = +s) plus inner face of layer 2 (n = -s): (3.8) twice
assert np.isclose(sig_b_b, (1 / e2 - 1 / e1) * lam / (2 * np.pi * b_c))            # = eps0 (E2 - E1) at b: Gauss, no free charge there
Qb_per_len = -(1 - 1 / e1) * lam + sig_b_b * 2 * np.pi * b_c + (1 - 1 / e2) * lam  # at a (n = -s), at b, at c (n = +s)
assert abs(Qb_per_len) < 1e-15 * lam
print(f"Your turn 1: C' = {C_coax * 1e12:.1f} pF/m; sigma_b(b) = {sig_b_b:.4e} C/m^2 for lambda = 1 nC/m (negative: eps_2 > eps_1); "
      f"total bound charge per length = {Qb_per_len:.1e}")
# 1, numerically with the module's own solver: a = 1, b = 2, c = 4 (units of a), eps_r = 2.3, 4.0, V = 1 on the inner conductor, 0 outside r > c.
# Conductors are fixed nodes (staircased circles); C'/eps0 = 2 W/V^2 from the discrete energy; E_s on the +x axis from central differences.
n_cx, L_cx = 361, 4.5
xs_cx = np.linspace(-L_cx, L_cx, n_cx); h_cx = xs_cx[1] - xs_cx[0]
ex_cx, ey_cx = face_eps(lambda X, Y: np.where(X**2 + Y**2 < 2.0**2, 2.3, 4.0), xs_cx, xs_cx)
r_cx = np.hypot(*np.meshgrid(xs_cx, xs_cx))
fixed_cx = (r_cx < 1.0) | (r_cx > 4.0)
V_cx = solve_poisson2d(ex_cx, ey_cx, np.zeros((n_cx, n_cx)), h_cx, fixed_cx, np.where(r_cx < 1.0, 1.0, 0.0))
C_cx_num = 2 * field_energy(V_cx, ex_cx, ey_cx)                                        # C'/eps0 per unit length, V = 1
C_cx_th = 2 * np.pi / (np.log(2) / 2.3 + np.log(2) / 4.0)
c0 = n_cx // 2
E_s_num = lambda s_: -(V_cx[c0, c0 + round(s_ / h_cx) + 1] - V_cx[c0, c0 + round(s_ / h_cx) - 1]) / (2 * h_cx)
E1_th, E2_th = C_cx_th / (2 * np.pi * 2.3 * 1.5), C_cx_th / (2 * np.pi * 4.0 * 3.0)
assert abs(C_cx_num / C_cx_th - 1) < 0.02 and abs(E_s_num(1.5) / E1_th - 1) < 0.03 and abs(E_s_num(3.0) / E2_th - 1) < 0.03
print(f"Your turn 1 (solver): C'/eps0 = {C_cx_num:.4f} vs analytic {C_cx_th:.4f} ({C_cx_num / C_cx_th - 1:+.2%}); "
      f"E_s(1.5) = {E_s_num(1.5):.4f} vs {E1_th:.4f}, E_s(3.0) = {E_s_num(3.0):.4f} vs {E2_th:.4f}")
# 2: water, computed in Part 7 (eps_r(CM) = 1.76 vs n^2 = 1.78 vs static 80.1; Debye + CM overshoots to eps_r < 0; Onsager 29).
# 3: piezoelectric quartz plate.
d11, stress, area, thick, er_q = 2.3e-12, 1e6, 1e-4, 1e-3, 4.5
P_pz = d11 * stress; Q_pz = P_pz * area
V_oc = P_pz * thick / (er_q * EPS0)          # D = 0 open circuit: eps0 eps_r E = -P ... |V| = P t / (eps0 eps_r)
assert np.isclose(P_pz, 2.3e-6) and np.isclose(Q_pz, 2.3e-10) and np.isclose(V_oc, 57.73, rtol=1e-3)
C_pz = er_q * EPS0 * area / thick                                         # plate capacitance, clamped
assert np.isclose(C_pz, 3.984e-12, rtol=1e-3) and np.isclose(0.5 * C_pz * V_oc**2, 0.5 * Q_pz * V_oc) and np.isclose(0.5 * Q_pz * V_oc, 6.64e-9, rtol=2e-3)
assert np.isclose(stress * area, 100.0) and np.isclose(d11 * (100.0 / thick) * thick, 2.3e-10)     # converse effect: S = d E, 100 V gives 0.23 nm
print(f"Your turn 3: sigma_b = {P_pz:.2e} C/m^2, Q = {Q_pz:.2e} C, open-circuit V = {V_oc:.1f} V, C = {C_pz * 1e12:.3f} pF, "
      f"stored energy {0.5 * Q_pz * V_oc:.3e} J, force {stress * area:.0f} N, converse: 100 V changes the thickness by {d11 * 100.0 * 1e9:.2f} nm")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
