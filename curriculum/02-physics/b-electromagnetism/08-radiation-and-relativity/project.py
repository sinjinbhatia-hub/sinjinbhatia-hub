# %% [markdown]
# # EM-8 Project: The Retarded-Time Solver, Dipole Radiation, and the Field Tensor
#
# Parts: (1) SymPy: the retarded point-source potential solves the wave equation (Example 4.1),
# the Lorentz transformation preserves eta, E^2 = p^2c^2 + m^2c^4, the velocity-addition gamma
# identity; (2) the retarded-time solver and Lienard-Wiechert fields, the head-on factor 1/(1 - beta)
# (Example 4.2), the radiation kink figure; (3) dipole radiation from the numerical fields: the
# sin^2 pattern, power independent of radius, the omega^4 law, Larmor; (4) the field tensor:
# random E, B boosted with Lambda F Lambda^T vs the boost rules, the invariants, Example 4.6, the
# squashed Coulomb field of a moving charge = the boosted static field; (5) SymPy: F = dA - dA
# reproduces E and B for the dipole potentials, Lorenz gauge, and finite-difference
# d_nu F^{mu nu} = 0, d_nu G^{mu nu} = 0 off the source; (6) numbers: the antenna (Example 4.3),
# the classical atom (Example 4.4), Purcell's wire (Example 4.5); Your turn with reference
# solutions; (7) every number in solutions.md.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from pathlib import Path
from scipy.integrate import solve_ivp, quad
from scipy.constants import mu_0 as MU0, epsilon_0 as EPS0, c as C, e as QE, m_e as ME
from scipy.constants import physical_constants as PC

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()
K = 1 / (4 * np.pi * EPS0)
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])          # metric (-,+,+,+), x^mu = (ct, x, y, z)


def close(a, b, rtol, label=""):
    """assert a ~ b to relative tolerance rtol; return the relative error"""
    err = np.max(np.abs(np.asarray(a) - np.asarray(b))) / np.max(np.abs(b))
    assert err < rtol, f"{label}: {a} vs {b}, rel err {err:.2e}"
    return err

# %% [markdown]
# ## Part 1: Symbolic checks (Example 4.1, notes 3.2, 3.12-3.13)
# V = f(t - r/c)/r solves the homogeneous wave equation for r > 0 for every f (Example 4.1).
# Lambda^T eta Lambda = eta for a boost (3.23). u.u = -c^2 and E^2 = p^2 c^2 + m^2 c^4 (3.25).
# gamma(w) = gamma(u) gamma(v) (1 + uv/c^2) for w = (u + v)/(1 + uv/c^2), used in 3.15 and 3.16.

# %%
t, x, y, z = sp.symbols("t x y z", real=True)
cs, bs, ms = sp.symbols("c beta m", positive=True)
f = sp.Function("f")
rr = sp.sqrt(x**2 + y**2 + z**2)
Vpt = f(t - rr / cs) / rr
box = sum(sp.diff(Vpt, s, 2) for s in (x, y, z)) - sp.diff(Vpt, t, 2) / cs**2
assert sp.simplify(box) == 0
g = 1 / sp.sqrt(1 - bs**2)
Lam = sp.Matrix([[g, -g * bs, 0, 0], [-g * bs, g, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
assert sp.simplify(Lam.T * sp.diag(-1, 1, 1, 1) * Lam - sp.diag(-1, 1, 1, 1)) == sp.zeros(4, 4)
vx = sp.symbols("v_x", real=True)
gv = 1 / sp.sqrt(1 - vx**2 / cs**2)
u4 = sp.Matrix([gv * cs, gv * vx, 0, 0])
assert sp.simplify((u4.T * sp.diag(-1, 1, 1, 1) * u4)[0] + cs**2) == 0
Erel, prel = gv * ms * cs**2, gv * ms * vx
assert sp.simplify(Erel**2 - prel**2 * cs**2 - ms**2 * cs**4) == 0
uu, vv = sp.symbols("u v", positive=True)
w_ = (uu + vv) / (1 + uu * vv / cs**2)
gam = lambda s: 1 / sp.sqrt(1 - s**2 / cs**2)
num = {cs: 1, uu: sp.Rational(3, 10), vv: sp.Rational(1, 2)}
assert sp.simplify((gam(w_) - gam(uu) * gam(vv) * (1 + uu * vv / cs**2)).subs(num)) == 0
print("Part 1: box[f(t-r/c)/r] = 0, Lambda^T eta Lambda = eta, E^2 = p^2c^2 + m^2c^4, gamma addition: OK")

# %% [markdown]
# ## Part 2: The retarded-time solver and the Lienard-Wiechert fields (notes 3.4-3.6)
# For a trajectory w(t), the retarded time solves g(t_r) = c(t - t_r) - |r - w(t_r)| = 0 (3.6).
# g is strictly decreasing when |v| < c (dg/dt_r = -c + rhat.v < 0), so there is exactly one
# root; we bracket it (g(t) <= 0, then step back until g > 0) and bisect, vectorized over every
# field point. Then (3.9): E = q/(4 pi eps0) * scr/(scr.u)^3 [(c^2 - v^2) u + scr x (u x a)],
# u = c scrhat - v, B = scrhat x E / c, and V = q c/(4 pi eps0 scr.u), A = v V / c^2 (3.8).

# %%
def retarded_time(traj, r, t, iters=64):
    """Solve c (t - t_r) = |r - w(t_r)| (3.6) by bracketing and bisection.
    r: (..., 3) field points; t: array broadcastable to r[..., 0]; traj(t) -> (w, v, a)."""
    t = np.broadcast_to(np.asarray(t, float), r.shape[:-1]).copy()
    gfun = lambda tr: C * (t - tr) - np.linalg.norm(r - traj(tr)[0], axis=-1)
    hi = t.copy()
    step = np.linalg.norm(r - traj(t)[0], axis=-1) / C + 1e-30
    lo = t - step
    while np.any(bad := gfun(lo) <= 0):
        step = np.where(bad, 2 * step, step)
        lo = np.where(bad, t - step, lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        pos = gfun(mid) > 0
        lo, hi = np.where(pos, mid, lo), np.where(pos, hi, mid)
    return 0.5 * (lo + hi)


def lw_fields(traj, q, r, t):
    """Lienard-Wiechert E, B (3.9) and potentials V, A (3.8) of a point charge q on traj."""
    tr = retarded_time(traj, r, t)
    w, v, a = traj(tr)
    R = r - w
    Rn = np.linalg.norm(R, axis=-1, keepdims=True)
    u = C * R / Rn - v
    Ru = np.sum(R * u, axis=-1, keepdims=True)
    v2 = np.sum(v * v, axis=-1, keepdims=True)
    E = q * K * Rn / Ru**3 * ((C**2 - v2) * u + np.cross(R, np.cross(u, a)))
    B = np.cross(R / Rn, E) / C
    V = q * K * C / Ru[..., 0]
    return E, B, V, v * V[..., None] / C**2, tr


def uniform(vel):
    """Trajectory w = vel t through the origin at t = 0."""
    vel = np.asarray(vel, float)
    def traj(tt):
        tt = np.asarray(tt, float)[..., None]
        return tt * vel, np.broadcast_to(vel, tt.shape[:-1] + (3,)), np.zeros(tt.shape[:-1] + (3,))
    return traj


# Example 4.2: charge moving straight at the observer, beta = 0.6, present distance 1 m.
q = 1e-9
beta = 0.6
E42, B42, V42, A42, tr42 = lw_fields(uniform([beta * C, 0, 0]), q, np.array([[1.0, 0, 0]]), 0.0)
scr42 = -C * tr42[0]
print(f"Example 4.2: t_r = {tr42[0]:.4e} s, retarded distance = {scr42:.4f} m, "
      f"1/(1 - beta) = {1/(1-beta):.4f}, V = {V42[0]:.4f} V, Coulomb at present position = {q*K:.4f} V, "
      f"Coulomb at retarded position = {q*K/scr42:.4f} V, E_x = {E42[0,0]:.4f} V/m = (1-beta^2) q/(4 pi eps0 R^2) = {(1-beta**2)*q*K:.4f}")
close(scr42, 1 / (1 - beta), 1e-12, "retarded distance")
close(V42[0], q * K, 1e-12, "head-on potential")
close(E42[0, 0], (1 - beta**2) * q * K, 1e-12, "head-on field")

# The kink: a charge at rest until t = 0, accelerated over tau to 0.6 c, then coasting.
VF, TAU = 0.6 * C, 0.05 / C


def kicked(tt):
    """x(t) = (v_f/2)(t + tau ln(2 cosh(t/tau))), v = (v_f/2)(1 + tanh), a = v_f sech^2/(2 tau)."""
    tt = np.asarray(tt, float)
    s = tt / TAU
    lncosh = np.logaddexp(s, -s) - np.log(2)
    w = np.zeros(tt.shape + (3,)); v = np.zeros_like(w); a = np.zeros_like(w)
    w[..., 0] = 0.5 * VF * (tt + TAU * lncosh + TAU * np.log(2))   # at rest at the origin as t -> -inf
    v[..., 0] = 0.5 * VF * (1 + np.tanh(s))
    a[..., 0] = 0.5 * VF / TAU / np.cosh(np.clip(s, -300, 300))**2
    return w, v, a


n = 240
xs = np.linspace(-3, 3, n)
XX, YY = np.meshgrid(xs, xs)
P = np.stack([XX, YY, np.zeros_like(XX)], axis=-1)
tsnap = 2.0 / C
Ek, *_ = lw_fields(kicked, 1.0, P, tsnap)
# outside the light sphere |r| > c t (to within the kick duration) the field is still the
# Coulomb field of the charge at rest at the origin: information has not arrived.
far = np.linalg.norm(P, axis=-1) > 2.0 + 12 * C * TAU
Ecoul = K * P / np.linalg.norm(P, axis=-1, keepdims=True)**3
close(Ek[far], Ecoul[far], 1e-6, "field outside the light sphere is the old Coulomb field")

# oscillating charge for the second panel: beta_max = 0.5, lambda = 1 m
LAM = 1.0; OMG = 2 * np.pi * C / LAM; AMP_OSC = 0.5 * C / OMG


def osc_traj(amp, omega):
    def traj(tt):
        tt = np.asarray(tt, float)
        w = np.zeros(tt.shape + (3,)); v = np.zeros_like(w); a = np.zeros_like(w)
        w[..., 2] = amp * np.cos(omega * tt)
        v[..., 2] = -amp * omega * np.sin(omega * tt)
        a[..., 2] = -amp * omega**2 * np.cos(omega * tt)
        return w, v, a
    return traj


Pxz = np.stack([XX, np.zeros_like(XX), YY], axis=-1)
Eo, *_ = lw_fields(osc_traj(AMP_OSC, OMG), 1.0, Pxz, 0.0)
fig, ax = plt.subplots(1, 2, figsize=(12, 5.6))
ax[0].streamplot(XX, YY, Ek[..., 0], Ek[..., 1], density=1.6, color="k", linewidth=0.6, arrowsize=0.5)
ax[0].add_patch(plt.Circle((0, 0), 2.0, fill=False, ls="--", color="tab:red"))
ax[0].plot(kicked(tsnap)[0][0], 0, "o", color="tab:blue")
ax[0].set(title="E lines 2 m/c after a kick to 0.6c: the kink shell", xlabel="x (m)", ylabel="y (m)", aspect="equal")
rE = np.linalg.norm(Pxz, axis=-1) * np.linalg.norm(Eo, axis=-1)
ax[1].imshow(np.log10(rE), extent=[-3, 3, -3, 3], origin="lower", cmap="magma", vmin=np.percentile(np.log10(rE), 5))
ax[1].set(title="log10(r|E|): oscillating charge, beta_max = 0.5, lambda = 1 m", xlabel="x (m)", ylabel="z (m)")
fig.tight_layout(); fig.savefig(FIG / "radiation_kink.png", dpi=110); plt.close(fig)
print("Part 2: retarded-time solver OK; outside the light sphere the field is the old Coulomb field")

# %% [markdown]
# ## Part 3: Dipole radiation from the numerical fields (notes 3.7-3.9)
# One charge q oscillating as z = A cos(wt) with beta_max = A w/c small is a dipole p0 = qA.
# Average S_r over a period on spheres of several radii. Theory (3.15)-(3.16):
# <S_r> = mu0 p0^2 w^4 sin^2(theta) / (32 pi^2 c r^2), P = mu0 p0^2 w^4 / (12 pi c), which is
# Larmor (3.19) with <a^2> = A^2 w^4 / 2. The time-averaged flux is the same through every
# sphere (energy conservation for a periodic source), so the radius test works even at r < lambda.

# %%
mu_nodes, mu_w = np.polynomial.legendre.leggauss(40)    # nodes in cos(theta)
NT = 64


def dipole_flux(amp, omega, R, qq=1e-9):
    """Period-averaged S_r(theta) and total power through a sphere of radius R."""
    th = np.arccos(mu_nodes)
    T = 2 * np.pi / omega
    tt = (np.arange(NT) * T / NT)[:, None]
    pts = R * np.stack([np.sin(th), 0 * th, np.cos(th)], axis=-1)[None, :, :].repeat(NT, 0)
    E, B, *_ = lw_fields(osc_traj(amp, omega), qq, pts, tt)
    Sr = np.sum(np.cross(E, B) * pts / R, axis=-1).mean(0) / MU0
    return th, Sr, 2 * np.pi * R**2 * np.sum(mu_w * Sr)


def p_dipole(p0, omega):
    return MU0 * p0**2 * omega**4 / (12 * np.pi * C)


qq, w0 = 1e-9, 2 * np.pi * C / LAM
amp = 0.005 * C / w0                                   # beta_max = 0.005
Pth = p_dipole(qq * amp, w0)
for R in [0.3, 1.0, 5.0, 50.0]:
    th, Sr, Pnum = dipole_flux(amp, w0, R)
    err = close(Pnum, Pth, 2e-3, f"power through r = {R} lambda")
    shape = Sr * R**2 / (MU0 * (qq * amp)**2 * w0**4 / (32 * np.pi**2 * C))
    close(shape, np.sin(th)**2, 2e-2, "sin^2 pattern")
    print(f"Part 3: r = {R:5.1f} lambda: P_num = {Pnum:.5e} W, P_dipole = {Pth:.5e} W, rel err {err:.1e}")
larmor = MU0 * qq**2 * (amp**2 * w0**4 / 2) / (6 * np.pi * C)
close(larmor, Pth, 1e-12, "Larmor = dipole formula")
oms = w0 * np.array([0.25, 0.5, 1.0])
amp_s = 0.004 * C / oms[-1]
Ps = [dipole_flux(amp_s, om, 3.0)[2] for om in oms]
slope = np.polyfit(np.log(oms), np.log(Ps), 1)[0]
assert abs(slope - 4) < 2e-3, slope
print(f"Part 3: fitted power-law exponent in omega = {slope:.5f} (theory 4)")
fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(5, 5))
th, Sr, _ = dipole_flux(amp, w0, 5.0)
ax.plot(th, Sr / Sr.max(), "o", ms=3, label="numerical <S_r>")
ax.plot(np.r_[th, -th[::-1]], np.r_[np.sin(th)**2, np.sin(th[::-1])**2], "-", label="sin^2 theta")
ax.set_theta_zero_location("N"); ax.legend(loc="lower right", fontsize=8)
fig.savefig(FIG / "dipole_pattern.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 4: The field tensor (notes 3.16-3.17) and the squashed Coulomb field (3.10)
# Griffiths' convention: F^{0i} = E_i/c, F^{ij} = eps_ijk B_k; G^{0i} = B_i, G^{ij} = -eps_ijk E_k/c.
# A boost with velocity v (beta = v/c) is Lambda^0_0 = gamma, Lambda^0_i = Lambda^i_0 = -gamma beta_i,
# Lambda^i_j = delta_ij + (gamma - 1) beta_i beta_j / beta^2. F' = Lambda F Lambda^T (rank-2 rule, 3.28).
# Boost rules (3.26) in vector form: E' = gamma(E + v x B) - gamma^2/(gamma+1) beta(beta.E),
# B' = gamma(B - v x E/c^2) - gamma^2/(gamma+1) beta(beta.B).

# %%
def F_tensor(E, B):
    Ex, Ey, Ez = np.asarray(E) / C
    Bx, By, Bz = B
    return np.array([[0, Ex, Ey, Ez], [-Ex, 0, Bz, -By], [-Ey, -Bz, 0, Bx], [-Ez, By, -Bx, 0]])


def EB_from_F(F):
    return C * F[0, 1:], np.array([F[2, 3], F[3, 1], F[1, 2]])


def G_tensor(E, B):
    Ex, Ey, Ez = np.asarray(E) / C
    Bx, By, Bz = B
    return np.array([[0, Bx, By, Bz], [-Bx, 0, -Ez, Ey], [-By, Ez, 0, -Ex], [-Bz, -Ey, Ex, 0]])


def boost(vel):
    b = np.asarray(vel, float) / C
    b2 = b @ b
    gm = 1 / np.sqrt(1 - b2)
    L = np.eye(4)
    L[0, 0] = gm; L[0, 1:] = L[1:, 0] = -gm * b
    L[1:, 1:] += (gm - 1) * np.outer(b, b) / b2
    return L


def boost_rules(E, B, vel):
    vel = np.asarray(vel, float); b = vel / C; gm = 1 / np.sqrt(1 - b @ b)
    Ep = gm * (E + np.cross(vel, B)) - gm**2 / (gm + 1) * b * (b @ E)
    Bp = gm * (B - np.cross(vel, E) / C**2) - gm**2 / (gm + 1) * b * (b @ B)
    return Ep, Bp


lower = lambda F: ETA @ F @ ETA
for trial in range(20):
    E = rng.normal(size=3) * 1e3; B = rng.normal(size=3) * 1e-5
    vel = rng.normal(size=3); vel *= rng.uniform(0.05, 0.95) * C / np.linalg.norm(vel)
    L = boost(vel)
    assert np.allclose(L.T @ ETA @ L, ETA)
    F = F_tensor(E, B); Fp = L @ F @ L.T
    Ep, Bp = EB_from_F(Fp)
    Er, Br = boost_rules(E, B, vel)
    close(Ep, Er, 1e-10, "E'"); close(Bp, Br, 1e-10, "B'")
    close(Ep @ Ep - C**2 * Bp @ Bp, E @ E - C**2 * B @ B, 1e-9, "E^2 - c^2 B^2")
    close(Ep @ Bp, E @ B, 1e-9, "E.B")
    close(np.sum(lower(F) * F), 2 * (B @ B - E @ E / C**2), 1e-12, "F_mn F^mn")
    close(np.sum(lower(F) * G_tensor(E, B)), -4 * E @ B / C, 1e-12, "F_mn G^mn")
print("Part 4: Lambda F Lambda^T = boost rules, both invariants, F F = 2(B^2 - E^2/c^2): OK (20 random trials)")

# Example 4.6: pure E = E0 x-hat, frame moving along y at 0.6 c
E0 = 1e6
Ep, Bp = EB_from_F(boost([0, 0.6 * C, 0]) @ F_tensor([E0, 0, 0], [0, 0, 0]) @ boost([0, 0.6 * C, 0]).T)
print(f"Example 4.6: E' = {Ep} V/m, B' = {Bp} T, E'^2 - c^2 B'^2 = {Ep@Ep - C**2*Bp@Bp:.6e} (E0^2 = {E0**2:.6e})")
close(Ep, [1.25e6, 0, 0], 1e-12); close(Bp, [0, 0, 1.25 * 0.6 * E0 / C], 1e-12)

# Squashed Coulomb field: Lienard-Wiechert for uniform motion vs (3.10) vs boosted static field.
bu = 0.8; vel = np.array([bu * C, 0, 0]); gu = 1 / np.sqrt(1 - bu**2)
ang = np.linspace(0, np.pi, 13)
Rp = np.stack([np.cos(ang), np.sin(ang), 0 * ang], axis=-1)       # present-position vectors, |R| = 1 m
Elw, Blw, *_ = lw_fields(uniform(vel), q, Rp, 0.0)
Eform = q * K * (1 - bu**2) / (1 - bu**2 * np.sin(ang)**2)**1.5
close(np.linalg.norm(Elw, axis=-1), Eform, 1e-10, "(3.10)")
L = boost(vel)                                    # lab -> rest frame
for i, R in enumerate(Rp):
    xr = L @ np.r_[0.0, R]                        # event (t = 0, R) in the charge's rest frame
    Erest = q * K * xr[1:] / np.linalg.norm(xr[1:])**3
    Elab, Blab = EB_from_F(np.linalg.inv(L) @ F_tensor(Erest, np.zeros(3)) @ np.linalg.inv(L).T)
    close(Elab, Elw[i], 1e-10, "boosted Coulomb = LW"); assert np.allclose(Blab, Blw[i], rtol=1e-10, atol=1e-10 * np.linalg.norm(Elw[i]) / C)
    assert np.allclose(Blw[i], np.cross(vel, Elw[i]) / C**2, rtol=1e-10, atol=1e-10 * np.linalg.norm(Elw[i]) / C)
print(f"Part 4: uniform motion beta = 0.8: E_perp/E_par at equal R = {Eform[6]/Eform[0]:.4f} = gamma^3 = {gu**3:.4f};"
      " LW field = boosted Coulomb field at 13 angles")

# %% [markdown]
# ## Part 5: Covariant Maxwell (notes 3.18, eqs. 3.30-3.31), symbolically and on a grid
# Exact potentials of a point dipole p(t) = p0 cos(wt) z-hat (the retarded integrals (3.2) without approximations (i)-(iii) of 3.8):
# V = (1/4 pi eps0)(z/r)[p/r^2 + pdot/(c r)], A = (mu0/4 pi r) pdot z-hat, all at t_r = t - r/c.
# Build A^mu = (V/c, A), F^{mu nu} = d^mu A^nu - d^nu A^mu with d^mu = eta^{mu nu} d_nu, and compare
# with E = -grad V - dA/dt, B = curl A. Then check d_nu F^{mu nu} = 0 and d_nu G^{mu nu} = 0 by
# fourth-order finite differences at points off the source, where J^mu = 0.

# %%
p0s, ws, e0s = sp.symbols("p_0 omega epsilon_0", positive=True)
mu0s = 1 / (e0s * cs**2)
R_ = sp.sqrt(x**2 + y**2 + z**2); tr_ = t - R_ / cs
pfun = p0s * sp.cos(ws * tr_); pdot = -p0s * ws * sp.sin(ws * tr_)
Vd = z / R_ * (pfun / R_**2 + pdot / (cs * R_)) / (4 * sp.pi * e0s)
Ad = [0, 0, mu0s * pdot / (4 * sp.pi * R_)]
X4 = [cs * t, x, y, z]
Aup = [Vd / cs] + Ad
d_low = lambda expr, mu: sp.diff(expr, t) / cs if mu == 0 else sp.diff(expr, X4[mu])
d_up = lambda expr, mu: -d_low(expr, 0) if mu == 0 else d_low(expr, mu)
Fs = [[d_up(Aup[nu], m) - d_up(Aup[m], nu) for nu in range(4)] for m in range(4)]
Evec = [-sp.diff(Vd, s) - sp.diff(Ad[i], t) for i, s in enumerate((x, y, z))]
Bvec = [sp.diff(Ad[2], y) - sp.diff(Ad[1], z), sp.diff(Ad[0], z) - sp.diff(Ad[2], x), sp.diff(Ad[1], x) - sp.diff(Ad[0], y)]
lorenz = sum(d_low(Aup[m], m) for m in range(4))
vals = {cs: 1.3, p0s: 0.7, ws: 2.1, e0s: 0.9}
checks = [Fs[0][i + 1] - Evec[i] / cs for i in range(3)] + [Fs[2][3] - Bvec[0], Fs[3][1] - Bvec[1], Fs[1][2] - Bvec[2], lorenz]
fchk = sp.lambdify((t, x, y, z), [ch.subs(vals) for ch in checks], "numpy")
for pt in rng.uniform(0.5, 2.0, size=(5, 4)):
    assert np.max(np.abs(fchk(*pt))) < 1e-12
print("Part 5: SymPy: F^{0i} = E_i/c, F^{ij} = eps_ijk B_k from d^mu A^nu - d^nu A^mu, Lorenz gauge holds")

Fnum = sp.lambdify((t, x, y, z), sp.Matrix(Fs).subs(vals), "numpy")


def Gnum(*X):
    Fm = np.array(Fnum(*X), float)
    return G_tensor(*EB_from_F(Fm))    # G only needs E/c = F^{0i} and B, so the units of c cancel


def divergence(Tfun, X, h=1e-3):
    """d_nu T^{mu nu} at event X = (t, x, y, z), fourth-order central differences; x^0 = c t."""
    out = np.zeros(4); scale = 0.0
    for nu in range(4):
        dX = np.zeros(4); dX[nu] = h / (1.3 if nu == 0 else 1.0)    # step in t is h/c
        D = (-np.array(Tfun(*(X + 2 * dX)), float) + 8 * np.array(Tfun(*(X + dX)), float)
             - 8 * np.array(Tfun(*(X - dX)), float) + np.array(Tfun(*(X - 2 * dX)), float)) / (12 * h)
        out += D[:, nu]; scale = max(scale, np.max(np.abs(D[:, nu])))
    return out, scale


for X in rng.uniform(0.6, 1.8, size=(4, 4)):
    for name, fn in (("F", Fnum), ("G", Gnum)):
        res, scale = divergence(fn, X)
        assert np.max(np.abs(res)) < 1e-8 * scale, (name, res, scale)
print("Part 5: finite-difference d_nu F^{mu nu} and d_nu G^{mu nu} vanish off the source (residual < 1e-8 of largest term)")

# %% [markdown]
# ## Part 6: The worked-example numbers (Examples 4.3-4.5)

# %%
# Example 4.3: Hertzian dipole, d = 1 m, f = 100 MHz, I0 = 1 A (uniform current assumed)
d, f0, I0 = 1.0, 100e6, 1.0
w = 2 * np.pi * f0; lam = C / f0
p0 = I0 * d / w
Pant = p_dipole(p0, w)
Rrad = 2 * Pant / I0**2
S1km = MU0 * p0**2 * w**4 / (32 * np.pi**2 * C) / 1e3**2
E1km = MU0 * p0 * w**2 / (4 * np.pi * 1e3)
print(f"Example 4.3: lambda = {lam:.4f} m, d/lambda = {d/lam:.4f}, p0 = {p0:.4e} C m, <P> = {Pant:.3f} W, "
      f"R_rad = {Rrad:.2f} ohm (= 2pi/3 Z0 (d/lam)^2 = {2*np.pi/3*MU0*C*(d/lam)**2:.2f}), "
      f"<S> broadside at 1 km = {S1km:.4e} W/m^2, E amplitude = {E1km:.4e} V/m")
close(Rrad, 2 * np.pi / 3 * MU0 * C * (d / lam)**2, 1e-10)
# half-wave dipole for comparison: R = (Z0/2pi) int_0^pi cos^2(pi/2 cos th)/sin th dth
Ihw = quad(lambda th: np.cos(np.pi / 2 * np.cos(th))**2 / np.sin(th), 0, np.pi)[0]
print(f"   2pi/3 Z0 = {2*np.pi/3*MU0*C:.1f} ohm, (d/lambda)^2 = {(d/lam)**2:.4f}, cyclotron period at 1 T = {2*np.pi*ME/QE:.3e} s")
print(f"   half-wave dipole: integral = {Ihw:.4f}, R_rad = {MU0*C/(2*np.pi)*Ihw:.2f} ohm")

# Example 4.4: classical hydrogen collapse
a0 = PC["Bohr radius"][0]; re = PC["classical electron radius"][0]
tf = a0**3 / (4 * re**2 * C)
sol = solve_ivp(lambda tt, r: -4 / 3 * re**2 * C / r**2, [0, tf], [a0], rtol=1e-10, atol=1e-25,
                events=lambda tt, r: r[0] - 1e-3 * a0)
v0 = np.sqrt(K * QE**2 / (ME * a0)); T_orb = 2 * np.pi * a0 / v0
acc0 = K * QE**2 / (ME * a0**2)
P0 = MU0 * QE**2 * acc0**2 / (6 * np.pi * C)
print(f"Example 4.4: a = {acc0:.4e} m/s^2, P = {P0:.4e} W, E = {-K*QE**2/(2*a0)/QE:.3f} eV, "
      f"t_collapse = {tf:.4e} s (ODE reaches 0.001 a0 at {sol.t_events[0][0]:.4e} s), "
      f"orbit period = {T_orb:.4e} s, revolutions = 2 t_f/T0 = {2*tf/T_orb:.4e}")
close(sol.t_events[0][0], tf * (1 - 1e-9), 1e-6)

# Example 4.5: Purcell's wire, I = 10 A, s = 1 cm, q = e moving at u = 1e6 m/s parallel to I
I, s, u = 10.0, 0.01, 1e6
gu = 1 / np.sqrt(1 - (u / C)**2)
F_lab = QE * u * MU0 * I / (2 * np.pi * s)
lam_bar = -gu * I * u / C**2
E_bar = lam_bar / (2 * np.pi * EPS0 * s)
F_bar = QE * abs(E_bar)
print(f"Example 4.5: B = {MU0*I/(2*np.pi*s):.4e} T, F_lab = {F_lab:.6e} N; in the charge's frame lambda' = {lam_bar:.4e} C/m, "
      f"E' = {E_bar:.4f} V/m, F' = {F_bar:.6e} N, F'/gamma_u = {F_bar/gu:.6e} N, gamma_u - 1 = {gu-1:.3e}")
close(F_bar / gu, F_lab, 1e-10)
# The density bookkeeping from first principles with a real drift speed (copper, 1 mm^2):
lam0 = 8.5e28 * QE * 1e-6; vd = I / lam0
gv = 1 / np.sqrt(1 - (vd / C)**2)
wv = (-vd - u) / (1 + u * vd / C**2); gw = 1 / np.sqrt(1 - (wv / C)**2)
close(gw, gu * gv * (1 + u * vd / C**2), 1e-14, "gamma addition")
print(f"   copper: lambda0 = {lam0:.4e} C/m, drift = {vd:.4e} m/s, fractional imbalance u v/c^2 = {u*vd/C**2:.3e}")

# %% [markdown]
# ## Your turn
# 1. **Synchrotron beaming.** Replace the oscillating charge by a circulating one,
#    w = R(cos wt, sin wt, 0), and plot r^2 |E|^2 around the orbit plane for beta = 0.1, 0.5, 0.9.
#    The reference below uses Lienard's angular formula (power per unit retarded time),
#    dP/dOmega = q^2/(16 pi^2 eps0 c) |n x ((n - beta) x betadot)|^2 / (1 - n.beta)^5,
#    checks that it integrates to Lienard's total gamma^4 mu0 q^2 a^2/(6 pi c) and that the
#    forward lobe narrows like 1/gamma.
# 2. **Beam steering (B9 extends this).** Two z-directed dipoles at x = +-d/2 with phases +-delta/2.
#    Show that in the xy-plane the power goes as cos^2((kd cos phi + delta)/2) and steer the beam.
# 3. **Abraham-Lorentz.** Integrate m a = m tau adot (no external force) forward: a runaway
#    a0 e^{t/tau}. Then impose the "no runaway" condition by integrating backward and see the
#    pre-acceleration a(t) = (F0/m) e^{t/tau} for t < 0 from a force switched on at t = 0.
# 4. **Exercise B9.** Run `dipole_flux` with the circulating trajectory and check P = mu0 q^2 w^4 R^2/(6 pi c)
#    and the pattern (1 + cos^2 theta)/2 (reference in Part 7).

# %%
def lienard_dOmega(nvec, beta_v, betadot, qq=1.0):
    nb = nvec - beta_v
    num_ = np.cross(nvec, np.cross(nb, betadot))
    return qq**2 / (16 * np.pi**2 * EPS0 * C) * np.sum(num_**2, -1) / (1 - nvec @ beta_v)**5


mu_n, mu_wt = np.polynomial.legendre.leggauss(80)
phis = np.linspace(0, 2 * np.pi, 7200, endpoint=False)
fig, ax = plt.subplots(figsize=(6, 4))
widths = []
for bb in [0.1, 0.5, 0.9, 0.99]:
    gb = 1 / np.sqrt(1 - bb**2)
    bv = np.array([bb, 0, 0]); acc = 1e15; bd = np.array([0, -acc / C, 0])   # a toward the centre
    TH, PH = np.meshgrid(np.arccos(mu_n), phis, indexing="ij")
    nv = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    dP = lienard_dOmega(nv, bv, bd)
    Ptot = np.sum(mu_wt[:, None] * dP) * (2 * np.pi / len(phis))
    close(Ptot, gb**4 * MU0 * acc**2 / (6 * np.pi * C), 1e-6, "Lienard total")
    plane = lienard_dOmega(np.stack([np.cos(phis), np.sin(phis), 0 * phis], -1), bv, bd)
    ang = np.angle(np.exp(1j * phis))
    half = ang[plane >= plane.max() / 2]
    widths.append((half.max() - half.min()) * gb)
    ax.plot(np.degrees(ang[np.argsort(ang)]), (plane / plane.max())[np.argsort(ang)], label=f"beta = {bb}")
ax.set(xlabel="angle from v in orbit plane (deg)", ylabel="dP/dOmega (normalized)", title="Synchrotron beaming")
ax.legend(); fig.tight_layout(); fig.savefig(FIG / "synchrotron.png", dpi=110); plt.close(fig)
print(f"[YT1] gamma x (full width at half max of forward lobe) for beta = 0.1, 0.5, 0.9, 0.99: {np.round(widths, 3)}")
assert abs(widths[3] / widths[2] - 1) < 0.15       # the lobe width scales like 1/gamma

kd = np.pi; deltas = [0, np.pi / 2, -np.pi / 2]
for dl in deltas:
    AF = np.cos((kd * np.cos(phis) + dl) / 2)**2
    pk = np.degrees(phis[np.argmax(AF)])
    expect = np.degrees(np.arccos(np.clip(-dl / kd, -1, 1)))
    assert abs(pk - expect) < 1.0 or abs(360 - pk - expect) < 1.0, (dl, pk, expect)
    print(f"[YT2] d = lambda/2, delta = {dl:.4f}: beam peak at phi = {pk:.1f} deg (theory {expect:.1f})")

tau_e = MU0 * QE**2 / (6 * np.pi * ME * C)
# time in units of tau: m a = F + m tau adot  =>  da/ds = a - F/m with s = t/tau
sol = solve_ivp(lambda s_, a: a, [0, 20], [1e-3], rtol=1e-11, atol=1e-15)
close(sol.y[0, -1], 1e-3 * np.exp(20), 1e-6, "runaway")
F0m = 1.0     # force switched on at t = 0; for t > 0 the non-runaway solution is a = F0/m
back = solve_ivp(lambda s_, a: a, [0, -5], [F0m], rtol=1e-11, atol=1e-15, dense_output=True)  # F = 0 for t < 0
close(back.sol(-1.0)[0], F0m * np.exp(-1), 1e-6, "pre-acceleration")
print(f"[YT3] tau_e = {tau_e:.4e} s; runaway grows by e^20 = {np.exp(20):.3e} in 20 tau; "
      f"pre-acceleration a(-tau)/(F0/m) = {back.sol(-1.0)[0]:.4f}")

# %% [markdown]
# ## Part 7: Numbers quoted in solutions.md

# %%
mp = PC["proton mass energy equivalent in MeV"][0]
taumu = 2.197e-6; bmu = 0.995; gmu = 1 / np.sqrt(1 - bmu**2)
print(f"[A1] gamma = {gmu:.4f}, lab lifetime = {gmu*taumu*1e6:.3f} us, range = {gmu*bmu*C*taumu/1e3:.4f} km, "
      f"naive = {bmu*C*taumu:.1f} m, 10 km atmosphere in muon frame = {10/gmu:.4f} km")
surv_lab = np.exp(-10e3 / (gmu * bmu * C * taumu)); surv_mu = np.exp(-(10e3 / gmu) / (bmu * C * taumu))
close(surv_lab, surv_mu, 1e-12)
print(f"[A1] surviving fraction: Earth frame {surv_lab:.4f}, muon frame {surv_mu:.4f}, without dilation {np.exp(-10e3/(bmu*C*taumu)):.2e}")
gp = 1 + 1000 / mp; bp = np.sqrt(1 - 1 / gp**2); pc_ = np.sqrt((mp + 1000)**2 - mp**2)
print(f"[A2] gamma = {gp:.4f}, beta = {bp:.4f}, E = {mp+1000:.2f} MeV, pc = {pc_:.1f} MeV, check gamma beta m c^2 = {gp*bp*mp:.1f} MeV")
close(gp * bp * mp, pc_, 1e-12)
print(f"[A3] (700/450)^4 = {(700/450)**4:.4f}")
aA4 = QE * 1e6 * 1.0 / ME; PA4 = MU0 * QE**2 * aA4**2 / (6 * np.pi * C); KE = 0.5 * ME * 1e12
print(f"[A4] a = {aA4:.4e} m/s^2, P = {PA4:.4e} W, KE = {KE:.4e} J, KE/P = {KE/PA4:.4f} s, "
      f"energy e-folding time = {3*np.pi*EPS0*ME**3*C**3/(QE**4*1.0**2):.4f} s")
ga5 = 1 / np.sqrt(1 - (1e3 / C)**2)
E5, B5 = boost_rules(np.zeros(3), np.array([0, 0, 1.0]), [1e3, 0, 0])
print(f"[A5] E' = {E5} V/m, B' = {B5} T, gamma - 1 = {ga5-1:.3e}, E^2 - c^2B^2 = {-C**2:.4e} V^2/m^2")
b3 = 0.9; g3 = 1 / np.sqrt(1 - b3**2)
print(f"[B3] gamma = {g3:.4f}, along-axis factor 1 - beta^2 = {1-b3**2:.4f}, broadside factor gamma = {g3:.4f}, ratio gamma^3 = {g3**3:.4f}")
bB4 = 0.01; fB4 = 100e6
print(f"[B4] pi b w / c = {np.pi*bB4*2*np.pi*fB4/C:.5f}, P_m/P_e for b = d = 1 cm at 100 MHz = (pi b w / c)^2 = {(np.pi*bB4*2*np.pi*fB4/C)**2:.4e}")
print(f"[B5] t_f = {tf:.4e} s, T0 = {T_orb:.4e} s, N = 2 t_f / T0 = {2*tf/T_orb:.4e}, "
      f"energy lost in first orbit / |E| = {P0*T_orb/(K*QE**2/(2*a0)):.3e}")
Ev7, Bv7 = np.array([0, 1e3, 0]), np.array([0, 0, 1.0])
vdrift = np.cross(Ev7, Bv7) / (Bv7 @ Bv7)
E7, B7 = boost_rules(Ev7, Bv7, vdrift)
assert np.allclose(E7, 0, atol=1e-9)
print(f"[B7] drift v = {vdrift} m/s, E' = {E7}, B' = {B7}")

# B9: circulating charge, R, w with beta = 0.005: P = mu0 q^2 w^4 R^2/(6 pi c), pattern (1 + cos^2)/2
Rc = 0.005 * C / w0


def circ(tt):
    tt = np.asarray(tt, float)
    w_ = np.zeros(tt.shape + (3,)); v_ = np.zeros_like(w_); a_ = np.zeros_like(w_)
    cth, sth = np.cos(w0 * tt), np.sin(w0 * tt)
    w_[..., 0], w_[..., 1] = Rc * cth, Rc * sth
    v_[..., 0], v_[..., 1] = -Rc * w0 * sth, Rc * w0 * cth
    a_[..., 0], a_[..., 1] = -Rc * w0**2 * cth, -Rc * w0**2 * sth
    return w_, v_, a_


th = np.arccos(mu_nodes); Rs = 5.0
tt = (np.arange(NT) * (2 * np.pi / w0) / NT)[:, None]
pts = Rs * np.stack([np.sin(th), 0 * th, np.cos(th)], -1)[None].repeat(NT, 0)
E, B, *_ = lw_fields(circ, qq, pts, tt)
Sr = np.sum(np.cross(E, B) * pts / Rs, -1).mean(0) / MU0
Pc = 2 * np.pi * Rs**2 * np.sum(mu_w * Sr)
Pc_th = MU0 * qq**2 * w0**4 * Rc**2 / (6 * np.pi * C)
close(Pc, Pc_th, 2e-3, "circulating power"); close(Sr / Sr.max(), (1 + np.cos(th)**2) / 2, 2e-2, "1 + cos^2 pattern")
print(f"[B9] circulating charge: P_num = {Pc:.5e} W, theory {Pc_th:.5e} W, ratio to one linear dipole of amplitude R: {Pc/p_dipole(qq*Rc, w0):.4f}")
print(f"[C1] tau_e = {tau_e:.4e} s, c tau_e = {C*tau_e:.4e} m (2/3 r_e = {2*re/3:.4e} m)")
for gC in [10, 1000]:
    print(f"[C2] gamma = {gC}: emission at theta' = 90 deg appears at sin(theta) = 1/gamma, theta = {np.degrees(np.arcsin(1/gC)):.4f} deg")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s. Figures in {FIG}")
