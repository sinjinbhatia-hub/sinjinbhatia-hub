# %% [markdown]
# # EM-6 Project: Induction, Maxwell's Equations, Energy and Momentum Flow, Gauge Freedom, and FDTD
#
# Parts: (1) Drude, the sliding bar and the falling loop (energy balances); (2) the induced E of a
# ramping solenoid as a loop integral on an off-centre circle; (3) inductance: Neumann's formula
# done numerically, coaxial loops against the exact elliptic result, a small loop in a solenoid
# computed both ways (M12 = M21), magnetic energy; (4) the charging capacitor; (5) the Poynting
# field of a resistor; (6) symbolic checks of Poynting's theorem and of the stress-tensor identity,
# and the force between two charges as a flux of stress through the midplane; (7) gauge invariance
# on a grid and the Lorenz-gauge dipole (SymPy); (8) a 1D Yee/FDTD solver: speed, PEC reflection,
# a dielectric interface, discrete energy conservation, open ends and the Poynting flux, the Courant
# limit; (9) a 2D TE FDTD point source; Your turn (lossy FDTD, Faraday's disk, two coils);
# (10) the numbers quoted in solutions.md.
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
from scipy.integrate import quad, solve_ivp
from scipy.special import ellipk, ellipe
from scipy.constants import mu_0 as MU0, epsilon_0 as EPS0, c as C, e as QE, m_e, g as G0

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()
Z0 = MU0 * C

def close(a, b, rtol, label=""):
    """assert a ~ b to relative tolerance rtol; return the relative error"""
    err = abs(a - b) / abs(b)
    assert err < rtol, f"{label}: {a} vs {b}, rel err {err:.2e}"
    return err

# %% [markdown]
# ## Part 1: Ohm, Drude, motional EMF (notes 3.1-3.2, Example 4.1)
# Drude: m dv/dt = qE - m v / tau gives sigma = n q^2 tau / m (3.2). Sliding bar: EMF = B l v,
# and the mechanical power F v equals I^2 R (Example 4.1). Coasting: all of (1/2) m v0^2 ends as heat.

# %%
n_cu, sig_cu = 8.49e28, 5.96e7
tau_cu = m_e * sig_cu / (n_cu * QE**2)
vd_cu = 1.0 / (n_cu * QE * 1e-6)
print(f"[1] copper: Drude tau = {tau_cu:.2e} s, drift speed at 1 A in 1 mm^2 = {vd_cu:.2e} m/s, eps0/sigma = {EPS0/sig_cu:.2e} s")
close(tau_cu, 2.5e-14, 0.01, "tau Cu")

Bb, lb, Rb, vb, mb = 0.8, 0.5, 2.0, 3.0, 0.1
emf_b = Bb * lb * vb; I_b = emf_b / Rb; F_b = I_b * lb * Bb
assert abs(F_b * vb - I_b**2 * Rb) < 1e-12
tau_b = mb * Rb / (Bb * lb)**2
sol = solve_ivp(lambda t, y: [-(Bb * lb)**2 * y[0] / (mb * Rb), (Bb * lb * y[0])**2 / Rb], (0, 20), [vb, 0], rtol=1e-10, atol=1e-12)
close(sol.y[1, -1], 0.5 * mb * vb**2, 1e-6, "coasting bar heat")
print(f"[1] Ex 4.1: EMF = {emf_b:.3f} V, I = {I_b:.3f} A, F = {F_b:.3f} N, Fv = I^2R = {F_b*vb:.3f} W; coasting tau = {tau_b:.3f} s, heat = {sol.y[1,-1]:.4f} J = mv0^2/2")

# %% [markdown]
# ## Part 2: Faraday's law: the induced E of a ramping solenoid (notes 3.3, Example 4.2)
# E_phi = -(mu0 n s / 2) dI/dt inside, -(mu0 n R^2 / 2 s) dI/dt outside (3.11). Check: the circulation
# of this E around an *off-centre* circle equals -dPhi/dt through it (lens-shaped overlap area).
# SymPy confirms curl E = -dB/dt and div E = 0 in Cartesian components.

# %%
n_s, R_s, dIdt = 2000.0, 0.02, 100.0
def E_sol(s, n=n_s, R=R_s, dI=dIdt):
    """|E_phi| of (3.11) (sign: opposes the change; returned as the -dI/dt coefficient)"""
    return np.where(s < R, MU0 * n * s / 2, MU0 * n * R**2 / (2 * np.maximum(s, 1e-30))) * dI
print(f"[2] Ex 4.2: mu0 n dI/dt = {MU0*n_s*dIdt:.4f} T/s; |E| at s = 1, 2, 4 cm: " + ", ".join(f"{E_sol(s)*1e3:.3f}" for s in (0.01, 0.02, 0.04)) + " mV/m")
xc, rc = 0.02, 0.03                                    # circle of radius 3 cm centred 2 cm off axis
ph = np.linspace(0, 2 * np.pi, 40001)[:-1]
X, Y = xc + rc * np.cos(ph), rc * np.sin(ph)
S = np.hypot(X, Y)
circ = np.sum(E_sol(S) * (-Y * (-rc * np.sin(ph)) + X * (rc * np.cos(ph))) / S) * (ph[1] - ph[0])
def lens(r1, r2, d):
    """overlap area of circles of radii r1, r2 with centres d apart"""
    a1 = r1**2 * np.arccos((d**2 + r1**2 - r2**2) / (2 * d * r1)); a2 = r2**2 * np.arccos((d**2 + r2**2 - r1**2) / (2 * d * r2))
    return a1 + a2 - 0.5 * np.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
close(circ, MU0 * n_s * dIdt * lens(R_s, rc, xc), 1e-6, "off-centre Faraday loop")
print(f"[2] off-centre circle: overlap area = {lens(R_s, rc, xc):.4e} m^2; circulation of E = {circ*1e6:.4f} uV, mu0 n (dI/dt) x overlap area = {MU0*n_s*dIdt*lens(R_s, rc, xc)*1e6:.4f} uV")
x, y, z, t = sp.symbols("x y z t", real=True)
k_, R_ = sp.symbols("k R", positive=True)              # k = mu0 n dI/dt
E_in = sp.Matrix([k_ * y / 2, -k_ * x / 2, 0])         # E_phi = -k s/2 phi-hat
E_out = sp.Matrix([k_ * R_**2 * y / (2 * (x**2 + y**2)), -k_ * R_**2 * x / (2 * (x**2 + y**2)), 0])
curl = lambda F: sp.Matrix([sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x), sp.diff(F[1], x) - sp.diff(F[0], y)])
div = lambda F: sp.diff(F[0], x) + sp.diff(F[1], y) + sp.diff(F[2], z)
assert sp.simplify(curl(E_in) - sp.Matrix([0, 0, -k_])) == sp.zeros(3, 1) and sp.simplify(curl(E_out)) == sp.zeros(3, 1)
assert sp.simplify(div(E_in)) == 0 and sp.simplify(div(E_out)) == 0
print("[2] SymPy: curl E = -dB/dt inside, 0 outside; div E = 0 both")

# %% [markdown]
# ## Part 3: Inductance and magnetic energy (notes 3.4-3.5, Example 4.3)
# Neumann (3.12): M = (mu0/4 pi) loop-loop integral of dl1.dl2/|r1 - r2|, done with the trapezoid rule
# on the exact parametrization (spectrally accurate for smooth closed curves). Coaxial loops are
# compared with the exact elliptic formula and the dipole approximation. A small tilted loop in a
# solenoid: M12 = flux of the solenoid's B through the loop (elliptic ring fields), M21 = sum over the
# solenoid's turns of the loop integral of the small loop's A. They must agree.

# %%
def circle_pts(R, n, c=(0, 0, 0), u=(1, 0, 0), v=(0, 1, 0)):
    """points r(phi) = c + R(cos phi u + sin phi v) and line elements dl = r'(phi) dphi"""
    p = 2 * np.pi * np.arange(n) / n
    u, v, c = map(np.asarray, (u, v, c))
    return c + R * (np.cos(p)[:, None] * u + np.sin(p)[:, None] * v), R * (-np.sin(p)[:, None] * u + np.cos(p)[:, None] * v) * 2 * np.pi / n

def neumann(P1, D1, P2, D2):
    """M = mu0/(4 pi) sum dl1.dl2 / |r1 - r2| (3.12)"""
    r = np.linalg.norm(P1[:, None, :] - P2[None, :, :], axis=2)
    return MU0 / (4 * np.pi) * np.sum((D1 @ D2.T) / r)

def M_coax(R1, R2, z):
    """exact mutual inductance of coaxial circles (Maxwell's elliptic-integral form)"""
    m = 4 * R1 * R2 / ((R1 + R2)**2 + z**2); k = np.sqrt(m)
    return MU0 * np.sqrt(R1 * R2) * ((2 / k - k) * ellipk(m) - 2 / k * ellipe(m))

R1, R2 = 0.05, 0.03
for zz in (0.01, 0.10, 0.50):
    Mn = neumann(*circle_pts(R1, 512), *circle_pts(R2, 512, c=(0, 0, zz)))
    Md = MU0 * np.pi * R1**2 * R2**2 / (2 * (R1**2 + zz**2)**1.5)          # small loop 2 in loop 1's axial field
    close(Mn, M_coax(R1, R2, zz), 1e-9, "Neumann vs elliptic")
    print(f"[3] Ex 4.3a: coaxial loops 5 cm, 3 cm, z = {zz*100:.0f} cm: Neumann {Mn:.5e} H, exact {M_coax(R1,R2,zz):.5e} H, on-axis-dipole {Md:.5e} H, far dipole mu0 pi R1^2 R2^2/2z^3 {MU0*np.pi*R1**2*R2**2/(2*zz**3):.5e} H")
zs = np.geomspace(0.01, 1, 60)
fig, ax = plt.subplots(figsize=(5, 3.6))
ax.loglog(zs, M_coax(R1, R2, zs), "C0", label="exact (elliptic)"); ax.loglog(zs, MU0 * np.pi * R1**2 * R2**2 / (2 * zs**3), "C1--", label=r"far dipole $\mu_0\pi R_1^2R_2^2/2z^3$")
ax.set_xlabel("separation z (m)"); ax.set_ylabel("M (H)"); ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(FIG / "mutual_coax.png", dpi=100); plt.close(fig)

def ring_B(a, I, rho, zr):
    """(B_rho, B_z) of a ring of radius a, current I, at cylindrical (rho, z) (elliptic integrals)"""
    m = 4 * a * rho / ((a + rho)**2 + zr**2); K, Ee = ellipk(m), ellipe(m)
    pre = MU0 * I / (2 * np.pi * np.sqrt((a + rho)**2 + zr**2)); den = (a - rho)**2 + zr**2
    return pre * zr / rho * (-K + (a**2 + rho**2 + zr**2) / den * Ee), pre * (K + (a**2 - rho**2 - zr**2) / den * Ee)

Lsol, rs, th = 0.5, 0.005, np.pi / 6
zk = (np.arange(int(n_s * Lsol)) + 0.5) / n_s - Lsol / 2                       # 1000 turns
nh, uh = np.array([np.sin(th), 0, np.cos(th)]), np.array([np.cos(th), 0, -np.sin(th)])
gq, gw = np.polynomial.legendre.leggauss(10)
rq, wq = rs * (gq + 1) / 2, rs / 2 * gw
pq = 2 * np.pi * np.arange(24) / 24
Pts = (rq[:, None, None] * (np.cos(pq)[None, :, None] * uh + np.sin(pq)[None, :, None] * np.array([0, 1, 0]))).reshape(-1, 3)
wts = (wq[:, None] * rq[:, None] * np.ones(24) * 2 * np.pi / 24).ravel()
rho_p = np.hypot(Pts[:, 0], Pts[:, 1])
Br, Bz = ring_B(R_s, 1.0, rho_p[None, :], Pts[None, :, 2] - zk[:, None])
Br, Bz = Br.sum(0), Bz.sum(0)
Bvec = np.stack([Br * Pts[:, 0] / rho_p, Br * Pts[:, 1] / rho_p, Bz], 1)
M12 = np.sum((Bvec @ nh) * wts)
Ps, Ds = circle_pts(rs, 96, u=uh, v=(0, 1, 0))
M21 = sum(neumann(*circle_pts(R_s, 96, c=(0, 0, zz)), Ps, Ds) for zz in zk)
M_inf = MU0 * n_s * np.pi * rs**2 * np.cos(th)
close(M12, M21, 1e-6, "M12 = M21"); close(M12, M_inf, 5e-3, "vs infinite solenoid")
fin = (Lsol / 2) / np.hypot(Lsol / 2, R_s); close(M12, M_inf * fin, 1e-4, "finite-solenoid centre factor")
print(f"[3] Ex 4.3b: small loop (5 mm, 30 deg) in 1000-turn solenoid: M12 (flux of B_sol) = {M12:.5e} H, M21 (sum of loop integrals of A_loop) = {M21:.5e} H, infinite-solenoid formula {M_inf:.5e} H, ratio {M12/M_inf:.5f}, centre-field factor (l/2)/sqrt((l/2)^2+R^2) = {fin:.5f}")

L_sol = MU0 * n_s**2 * np.pi * R_s**2 * Lsol
Nt, ht, at, bt = 500, 0.01, 0.04, 0.06
L_tor = MU0 * Nt**2 * ht / (2 * np.pi) * np.log(bt / at)
W_tor = quad(lambda s: (MU0 * Nt * 1.0 / (2 * np.pi * s))**2 / (2 * MU0) * ht * 2 * np.pi * s, at, bt)[0]
close(W_tor, 0.5 * L_tor, 1e-10, "toroid energy")
print(f"[3] L(solenoid n = 2000/m, R = 2 cm, 0.5 m) = {L_sol:.4e} H; L(toroid N = 500, h = 1 cm, 4-6 cm) = {L_tor:.4e} H; (1/2mu0) int B^2 at 1 A = {W_tor:.4e} J = L/2")
print(f"[3] energy densities: B^2/2mu0 at 1 T = {1/(2*MU0):.4e} J/m^3; eps0 E^2/2 at 3 MV/m = {EPS0*9e12/2:.2f} J/m^3")

# %% [markdown]
# ## Part 4: Maxwell's correction: the charging capacitor (notes 3.6, Example 4.4)

# %%
ac, dc, Ic, Vc = 0.05, 1e-3, 2.0, 100.0
dEdt = Ic / (EPS0 * np.pi * ac**2)
B_gap = lambda s: MU0 * Ic * s / (2 * np.pi * ac**2)
Cc = EPS0 * np.pi * ac**2 / dc
S_rim = (Vc / dc) * B_gap(ac) / MU0
close(S_rim * 2 * np.pi * ac * dc, Vc * Ic, 1e-12, "capacitor power")
print(f"[4] Ex 4.4: dE/dt = {dEdt:.4e} V/(m s); B(2 cm) = {B_gap(0.02)*1e6:.3f} uT, B(a) = {B_gap(ac)*1e6:.3f} uT, wire B at 10 cm = {MU0*Ic/(2*np.pi*0.1)*1e6:.3f} uT; C = {Cc*1e12:.2f} pF, dV/dt = {Ic/Cc:.3e} V/s; at V = 100 V: S(rim) = {S_rim:.4e} W/m^2, power in = {S_rim*2*np.pi*ac*dc:.1f} W")

# %% [markdown]
# ## Part 5: Poynting's theorem: energy enters a resistor through its sides (notes 3.9, Example 4.5)
# Uniform J along a cylinder; E = J/sigma (tangential E continuous at the surface); B from Ampere.
# S = E x B / mu0 on a grid in the cross-section. The inward flux through a *square* around the wire
# is I^2 R per length; through a square inside, it is the Joule heat inside; div S = -E.J on the grid.

# %%
rho_nc, a_w, L_w, I_w = 1.10e-6, 0.5e-3, 1.0, 2.0
R_w = rho_nc * L_w / (np.pi * a_w**2); E_w = I_w * R_w / L_w; J_w = I_w / (np.pi * a_w**2)
def S_field(X, Y):
    """Poynting vector (x, y components) of the resistor: S = -(E B_phi / mu0) s-hat"""
    s = np.hypot(X, Y); Bphi = np.where(s < a_w, MU0 * I_w * s / (2 * np.pi * a_w**2), MU0 * I_w / (2 * np.pi * np.maximum(s, 1e-30)))
    Sm = -E_w * Bphi / MU0
    return Sm * X / np.maximum(s, 1e-30), Sm * Y / np.maximum(s, 1e-30)
def inflow_square(h, n=4000):
    """-(closed line integral of S.n) around the square |x|,|y| <= h (per unit length of wire)"""
    u = -h + 2 * h * (np.arange(n) + 0.5) / n; du = 2 * h / n
    return -sum(np.sum(sgn * S_field(*(xy))[k]) * du for xy, k, sgn in (((h + 0 * u, u), 0, 1), ((-h + 0 * u, u), 0, -1), ((u, h + 0 * u), 1, 1), ((u, -h + 0 * u), 1, -1)))
P_out, P_in = inflow_square(1.5 * a_w), inflow_square(0.5 * a_w)
close(P_out * L_w, I_w**2 * R_w, 1e-6, "I^2R through a square"); close(P_in, E_w * J_w * (a_w)**2, 1e-6, "Joule heat inside a square")
g = np.linspace(-2 * a_w, 2 * a_w, 401); GX, GY = np.meshgrid(g, g, indexing="ij"); Sx, Sy = S_field(GX, GY)
divS = np.gradient(Sx, g, axis=0) + np.gradient(Sy, g, axis=1)
inside = np.hypot(GX, GY) < 0.9 * a_w
close(divS[inside].mean(), -E_w * J_w, 1e-3, "div S = -E.J")
assert np.abs(divS[np.hypot(GX, GY) > 1.1 * a_w]).max() < 1e-3 * E_w * J_w
print(f"[5] Ex 4.5: R = {R_w:.4f} Ohm, V = {I_w*R_w:.4f} V, E = {E_w:.4f} V/m, B(a) = {MU0*I_w/(2*np.pi*a_w)*1e3:.3f} mT, |S(a)| = {E_w*I_w/(2*np.pi*a_w):.1f} W/m^2; inflow through square 3a x 3a = {P_out*L_w:.4f} W = I^2R = {I_w**2*R_w:.4f} W; div S inside = {divS[inside].mean():.4e} vs -E.J = {-E_w*J_w:.4e} W/m^3")
fig, ax = plt.subplots(figsize=(4.5, 4.2)); sk = slice(None, None, 25)
ax.quiver(GX[sk, sk] * 1e3, GY[sk, sk] * 1e3, Sx[sk, sk], Sy[sk, sk], color="C3"); ax.add_patch(plt.Circle((0, 0), a_w * 1e3, fill=False))
ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_title("Poynting vector of a resistor (cross-section)"); fig.tight_layout(); fig.savefig(FIG / "poynting_resistor.png", dpi=100); plt.close(fig)

# %% [markdown]
# ## Part 6: Poynting's theorem and the stress tensor, symbolically and as a midplane integral (notes 3.9-3.10, Example 4.6)
# Fields from arbitrary polynomial potentials satisfy div B = 0 and Faraday automatically. Define rho and
# J by Gauss and Maxwell-Ampere. Then E.J + du/dt + div S = 0 (3.28) and rho E + J x B = div T - eps0 mu0 dS/dt
# (Theorem 3.8, step 6) must hold identically. Then F on charge 1 = -integral of T.z over the midplane (outward normal -z).

# %%
e0, m0 = sp.symbols("epsilon_0 mu_0", positive=True)
X4 = (x, y, z, t)
rp = lambda: sum(sp.Integer(int(c_)) * x**i * y**j * z**k * t**l for c_, (i, j, k, l) in zip(rng.integers(-3, 4, 12), rng.integers(0, 3, (12, 4))))
Vp, Ap = rp(), sp.Matrix([rp(), rp(), rp()])
grad = lambda f: sp.Matrix([sp.diff(f, v) for v in (x, y, z)])
Ef, Bf = -grad(Vp) - sp.diff(Ap, t), curl(Ap)
rho_f, Jf = e0 * div(Ef), curl(Bf) / m0 - e0 * sp.diff(Ef, t)
u_f, S_f = (e0 * Ef.dot(Ef) + Bf.dot(Bf) / m0) / 2, Ef.cross(Bf) / m0
assert sp.expand(Ef.dot(Jf) + sp.diff(u_f, t) + div(S_f)) == 0
T_f = sp.Matrix(3, 3, lambda i, j: e0 * (Ef[i] * Ef[j] - sp.KroneckerDelta(i, j) * Ef.dot(Ef) / 2) + (Bf[i] * Bf[j] - sp.KroneckerDelta(i, j) * Bf.dot(Bf) / 2) / m0)
divT = sp.Matrix([sum(sp.diff(T_f[i, j], X4[j]) for j in range(3)) for i in range(3)])
assert sp.expand(rho_f * Ef + Jf.cross(Bf) - divT + e0 * m0 * sp.diff(S_f, t)) == sp.zeros(3, 1)
print("[6] SymPy: Poynting's theorem and f = div T - eps0 mu0 dS/dt hold identically for random polynomial potentials")

def maxwell_T(E, B):
    """T_ij = eps0 (E_i E_j - d_ij E^2/2) + (B_i B_j - d_ij B^2/2)/mu0, for arrays (..., 3) (3.29)"""
    I3 = np.eye(3)
    return EPS0 * (np.einsum("...i,...j->...ij", E, E) - 0.5 * np.einsum("...k,...k", E, E)[..., None, None] * I3) + (np.einsum("...i,...j->...ij", B, B) - 0.5 * np.einsum("...k,...k", B, B)[..., None, None] * I3) / MU0

def E_point(q, r0, P):
    d = P - r0; return q / (4 * np.pi * EPS0) * d / np.linalg.norm(d, axis=-1, keepdims=True)**3

qq, dd = 1e-6, 0.01
uq, wu = np.polynomial.legendre.leggauss(200)
th_ = np.pi / 4 * (uq + 1); s_ = dd * np.tan(th_); ws_ = np.pi / 4 * wu * dd / np.cos(th_)**2       # s = d tan(theta) maps (0, inf)
phis = 2 * np.pi * np.arange(32) / 32
P = np.stack([s_[:, None] * np.cos(phis), s_[:, None] * np.sin(phis), 0 * s_[:, None] + 0 * phis], -1)
dA = (ws_ * s_)[:, None] * (2 * np.pi / 32) * np.ones(32)
F_coul = qq**2 / (4 * np.pi * EPS0 * (2 * dd)**2)
for q2, sign in ((qq, +1), (-qq, -1)):
    E_ = E_point(qq, np.array([0, 0, dd]), P) + E_point(q2, np.array([0, 0, -dd]), P)
    T_ = maxwell_T(E_, 0 * E_)
    F = -np.einsum("abiz,ab->i", T_[..., :, 2:3], dA)
    close(F[2], sign * F_coul, 1e-8, "midplane stress force"); assert abs(F[0]) + abs(F[1]) < 1e-10 * F_coul
    print(f"[6] Ex 4.6: charges {qq:.0e} C and {q2:+.0e} C, 2 cm apart: stress-tensor force on the upper charge F_z = {F[2]:+.4f} N; Coulomb {sign*F_coul:+.4f} N")
ev = np.linalg.eigvalsh(maxwell_T(np.array([1e6, 0, 0]), np.zeros(3)))
close(ev.max(), EPS0 * 1e12 / 2, 1e-12, "tension"); close(ev.min(), -EPS0 * 1e12 / 2, 1e-12, "pressure")
print(f"[6] T for E = 1 MV/m x-hat: eigenvalues {ev} Pa (tension along the line, pressure across); B = 1 T: B^2/2mu0 = {1/(2*MU0):.4e} Pa = {1/(2*MU0)/101325:.2f} atm")

# %% [markdown]
# ## Part 7: Gauge invariance on a grid, and the Lorenz-gauge oscillating dipole (notes 3.11-3.12)
# E = -grad V - dA/dt, B = curl A by central differences on an (x, y, z, t) grid. A gauge change with an
# arbitrary lambda changes V and A by O(1) but E and B only by round-off, because difference operators
# along different axes commute (the discrete version of the proof of Theorem 3.9).

# %%
h4 = 0.05; ax1 = np.arange(-0.6, 0.6 + 1e-9, h4); tt = np.array([-h4, 0, h4])
TT, XX, YY, ZZ = np.meshgrid(tt, ax1, ax1, ax1, indexing="ij")
V0 = np.sin(2 * XX - 3 * TT) * np.exp(-YY**2) + XX * ZZ
A0 = np.stack([np.cos(YY) * ZZ * TT, XX * TT**2 + YY, np.sin(XX * YY + TT)])
lam = 5 * np.exp(-(XX**2 + 2 * YY**2 + ZZ**2)) * np.cos(4 * TT) + XX * YY * ZZ * TT
D = lambda f, k: np.gradient(f, h4, axis=k)                    # axes: 0 = t, 1 = x, 2 = y, 3 = z
def fields(V, A):
    E = np.stack([-D(V, i + 1) - D(A[i], 0) for i in range(3)])
    B = np.stack([D(A[2], 2) - D(A[1], 3), D(A[0], 3) - D(A[2], 1), D(A[1], 1) - D(A[0], 2)])
    return E[:, 1], B[:, 1]
E1g, B1g = fields(V0, A0)
E2g, B2g = fields(V0 - D(lam, 0), A0 + np.stack([D(lam, k) for k in (1, 2, 3)]))
gerr = max(np.abs(E2g - E1g).max() / np.abs(E1g).max(), np.abs(B2g - B1g).max() / np.abs(B1g).max())
assert gerr < 1e-12
print(f"[7] gauge change: max |dV| = {np.abs(D(lam,0)).max():.2f}, max |dA| = {np.abs(np.stack([D(lam,k) for k in (1,2,3)])).max():.2f}; relative change in E and B = {gerr:.1e}")

p0, w_, c_ = sp.symbols("p_0 omega c", positive=True)
r = sp.sqrt(x**2 + y**2 + z**2); tr = t - r / c_
Vd = p0 * z / (4 * sp.pi * e0 * r**2) * (-(w_ / c_) * sp.sin(w_ * tr) + sp.cos(w_ * tr) / r)
Ad = sp.Matrix([0, 0, -p0 * w_ * sp.sin(w_ * tr) / (4 * sp.pi * e0 * c_**2 * r)])          # mu0 = 1/(eps0 c^2)
lorenz = div(Ad) + sp.diff(Vd, t) / c_**2
box = lambda f: sum(sp.diff(f, v, 2) for v in (x, y, z)) - sp.diff(f, t, 2) / c_**2
Ed, Bd = -grad(Vd) - sp.diff(Ad, t), curl(Ad)
checks = [lorenz, box(Vd), box(Ad[2]), div(Ed)] + list(curl(Bd) - sp.diff(Ed, t) / c_**2) + list(curl(Ed) + sp.diff(Bd, t))
for pt in rng.uniform(-2, 2, (3, 4)):
    sub = {v: sp.Rational(int(1000 * q), 1000) for v, q in zip((x, y, z, t), pt)} | {p0: 1, w_: sp.Rational(7, 3), c_: sp.Rational(3, 2), e0: sp.Rational(1, 2)}
    assert max(abs(sp.N(f.subs(sub), 40)) for f in checks) < 1e-30
print("[7] SymPy: oscillating-dipole potentials satisfy the Lorenz condition, box V = box A = 0 (r > 0), and E, B from them satisfy all four vacuum Maxwell equations")

# %% [markdown]
# ## Part 8: 1D FDTD (Yee) for E = E_y(x) y-hat, B = B_z(x) z-hat (notes section 7)
# Faraday: dB_z/dt = -dE_y/dx. Maxwell-Ampere with conductivity: eps dE_y/dt = -(1/mu0) dB_z/dx - sigma E_y.
# Put E_y on integer nodes x_i at integer times t_n, B_z on half nodes x_{i+1/2} at half times t_{n+1/2}.
# Centred differences in both space and time give the leapfrog
#   B^{n+1/2}_{i+1/2} = B^{n-1/2}_{i+1/2} - (dt/dx)(E^n_{i+1} - E^n_i)
#   E^{n+1}_i = c_a E^n_i - c_b (B^{n+1/2}_{i+1/2} - B^{n+1/2}_{i-1/2}),
#   c_a = (1 - sigma dt/2eps)/(1 + sigma dt/2eps),  c_b = (dt/(eps mu0 dx))/(1 + sigma dt/2eps)
# (sigma E averaged between n and n+1). Both updates are second-order accurate. Stability needs the
# Courant number S = c dt/dx <= 1: the scheme conserves W = sum[eps E^n.E^n + B^{n-1/2}.B^{n+1/2}/mu0] dx/2
# exactly, and W is positive definite only when S <= 1 (notes section 7).

# %%
def yee1d(E, B, S, steps, eps_r=1.0, sig=0.0, dx=1e-3, bc="pec", probe=None, snaps=()):
    """Run the leapfrog above. E: (N,) at t=0; B: (N-1,) at t=-dt/2. Returns E, B, energy invariant, probe flux, snapshots."""
    E, B = E.copy(), B.copy(); dt = S * dx / C
    eps = EPS0 * np.broadcast_to(eps_r, E.shape); sg = np.broadcast_to(sig, E.shape)
    ca = (1 - sg * dt / (2 * eps)) / (1 + sg * dt / (2 * eps)); cb = dt / (eps * MU0 * dx) / (1 + sg * dt / (2 * eps))
    kmur = (S - 1) / (S + 1); W, flux, out = [], 0.0, {}
    for n in range(steps):
        Bold = B.copy()
        B -= dt / dx * (E[1:] - E[:-1])
        W.append(0.5 * dx * (np.sum(eps * E**2) + np.sum(Bold * B) / MU0))
        Eold = E.copy()
        E[1:-1] = ca[1:-1] * E[1:-1] - cb[1:-1] * (B[1:] - B[:-1])
        if bc == "mur":                                    # first-order Mur absorbing ends
            E[0] = Eold[1] + kmur * (E[1] - Eold[0]); E[-1] = Eold[-2] + kmur * (E[-2] - Eold[-1])
        if probe is not None:                              # S_x = E_y B_z / mu0 at node `probe`, time n + 1/2
            flux += dt * 0.5 * (E[probe] + Eold[probe]) * 0.5 * (B[probe - 1] + B[probe]) / MU0
        if n + 1 in snaps: out[n + 1] = E.copy()
    return E, B, np.array(W), flux, out

def gauss_pulse(N, x0, w, S, dx=1e-3, n_ref=1.0):
    """right-moving pulse: E(x, 0) = exp(-(x - x0)^2/w^2), B(x, -dt/2) = E(x + v dt/2)/v, v = c/n"""
    xs = np.arange(N) * dx; xh = xs[:-1] + dx / 2; v = C / n_ref; dt = S * dx / C
    return xs, np.exp(-((xs - x0) / w)**2), np.exp(-((xh + v * dt / 2 - x0) / w)**2) / v

N1, S1 = 2001, 0.5
xs, E0, B0 = gauss_pulse(N1, 0.4, 0.02, S1)
snaps = (0, 1000, 2000, 4000, 6000)
cent = []
for k in range(1, 6):                                      # speed: energy centroid after 200 k steps
    Ek = yee1d(E0, B0, S1, 200 * k)[0]; cent.append(np.sum(xs * Ek**2) / np.sum(Ek**2))
v_num = np.polyfit(200 * np.arange(1, 6) * S1 * 1e-3 / C, cent, 1)[0]
close(v_num, C, 1e-3, "FDTD pulse speed")
Ep, Bp, Wp, _, snp = yee1d(E0, B0, S1, 6000, snaps=snaps)
assert np.ptp(Wp) / Wp[0] < 1e-12
assert -1.01 < snp[6000].min() < -0.98 and snp[6000].max() < 0.01        # PEC: reflected E inverted
print(f"[8] FDTD S = 0.5: pulse speed = {v_num:.6e} m/s (c = {C:.6e}, rel. err {abs(v_num/C-1):.1e}); discrete energy constant to {np.ptp(Wp)/Wp[0]:.1e}; after PEC reflection E ranges [{snp[6000].min():.4f}, {snp[6000].max():.4f}]")

i_int = 1000; epsr = np.where(np.arange(N1) > i_int, 4.0, 1.0); epsr[i_int] = 2.5
_, Ew0, Bw0 = gauss_pulse(N1, 0.4, 0.04, S1)                 # wider pulse: it is compressed by n = 2 in the glass
Ed_, Bd_, Wd, _, snd = yee1d(Ew0, Bw0, S1, 2600, eps_r=epsr, snaps=(2600,))
Ed_ = snd[2600]
r_num, t_num = Ed_[:i_int].min(), Ed_[i_int:].max()
U_E = 0.5 * EPS0 * epsr * Ed_**2
R_num = U_E[:i_int].sum() / U_E.sum(); T_num = U_E[i_int:].sum() / U_E.sum()
close(r_num, -1 / 3, 5e-3, "reflection amplitude"); close(t_num, 2 / 3, 5e-3, "transmission amplitude")
close(R_num, 1 / 9, 5e-3, "R"); close(T_num, 8 / 9, 1e-3, "T"); assert np.ptp(Wd) / Wd[0] < 1e-12
print(f"[8] dielectric n = 2: reflected amplitude {r_num:.4f} (theory -1/3), transmitted {t_num:.4f} (2/3); energy fractions R = {R_num:.4f} (1/9 = 0.1111), T = {T_num:.4f} (8/9 = 0.8889)")

Eo, Bo, Wo, flux, _ = yee1d(E0, B0, S1, 6000, bc="mur", probe=1800)
W_pulse = Wp[0]
close(flux, W_pulse, 2e-3, "Poynting flux through the probe plane")
assert Wo[-1] / W_pulse < 1e-4
print(f"[8] open ends (Mur): energy through x = 1.8 m = integral of S dt = {flux:.5e} J/m^2 vs pulse energy {W_pulse:.5e} J/m^2; left in domain {Wo[-1]/W_pulse:.1e}")

Eu = yee1d(E0, B0, 1.02, 600)[0]
assert np.abs(Eu).max() > 1e3
print(f"[8] Courant violated (S = 1.02): max |E| after 600 steps = {np.abs(Eu).max():.2e} (started at 1)")
fig, axs = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
for n_, Es in snp.items(): axs[0].plot(xs, Es, lw=1, label=f"step {n_}")
axs[0].plot(xs, E0, "k:", lw=1); axs[0].set_title("vacuum, PEC walls (S = 0.5)"); axs[0].legend(fontsize=7); axs[0].set_ylabel("E_y")
axs[1].plot(xs, Ed_, "C3"); axs[1].axvspan(xs[i_int], xs[-1], color="0.9"); axs[1].set_title("dielectric interface: n = 2 for x > 1 m (step 2600)")
axs[1].set_xlabel("x (m)"); axs[1].set_ylabel("E_y"); fig.tight_layout(); fig.savefig(FIG / "fdtd_1d.png", dpi=100); plt.close(fig)

# %% [markdown]
# ## Part 9: 2D FDTD, TE mode (E_x, E_y, H_z) with a point source
# dB_z/dt = -(dE_y/dx - dE_x/dy), eps0 dE_x/dt = dH_z/dy, eps0 dE_y/dt = -dH_z/dx; 2D Courant limit S <= 1/sqrt 2.

# %%
N2, dx2, S2 = 401, 1e-3, 0.5
dt2 = S2 * dx2 / C
Hz = np.zeros((N2, N2)); Ex = np.zeros((N2, N2 + 1)); Ey = np.zeros((N2 + 1, N2)); ic = N2 // 2
src = lambda n: (1 - 2 * (np.pi * (n - 60) / 30)**2) * np.exp(-(np.pi * (n - 60) / 30)**2)      # Ricker pulse
fig, axs = plt.subplots(2, 2, figsize=(7, 7)); shots, rp1, rp2, rec = (120, 220, 320, 420), 40, 120, []
for n in range(1, 521):
    Hz -= dt2 / (MU0 * dx2) * ((Ey[1:, :] - Ey[:-1, :]) - (Ex[:, 1:] - Ex[:, :-1]))
    Hz[ic, ic] += src(n)
    Ex[:, 1:-1] += dt2 / (EPS0 * dx2) * (Hz[:, 1:] - Hz[:, :-1])
    Ey[1:-1, :] -= dt2 / (EPS0 * dx2) * (Hz[1:, :] - Hz[:-1, :])
    rec.append((Hz[ic + rp1, ic], Hz[ic + rp2, ic]))
    if n in shots:
        a_ = axs.flat[shots.index(n)]; a_.imshow(Hz.T, origin="lower", cmap="RdBu", vmin=-0.03, vmax=0.03, extent=[-200, 200, -200, 200])
        a_.set_title(f"H_z, step {n}"); a_.set_xticks([]); a_.set_yticks([])
        assert np.abs(Hz - Hz.T).max() < 1e-12 * np.abs(Hz).max()
fig.tight_layout(); fig.savefig(FIG / "fdtd_2d.png", dpi=100); plt.close(fig)
rec = np.array(rec); xc_ = np.correlate(rec[:, 1], rec[:, 0], "full"); k0 = np.argmax(xc_)
lag = k0 - (len(rec) - 1) + 0.5 * (xc_[k0 - 1] - xc_[k0 + 1]) / (xc_[k0 - 1] - 2 * xc_[k0] + xc_[k0 + 1])
v2 = (rp2 - rp1) * dx2 / (lag * dt2); amp_ratio = np.abs(rec[:, 1]).max() / np.abs(rec[:, 0]).max()
close(v2, C, 0.02, "2D wave speed"); close(amp_ratio, np.sqrt(rp1 / rp2), 0.1, "1/sqrt(r) spreading")
print(f"[9] 2D TE: arrival lag between r = {rp1} and {rp2} cells = {lag:.2f} steps: speed = {v2:.4e} m/s; peak ratio {amp_ratio:.3f} vs sqrt(r1/r2) = {np.sqrt(rp1/rp2):.3f}; Hz symmetric under x <-> y")

# %% [markdown]
# ## Your turn
# 1. **Lossy FDTD.** Give the right half of the 1D domain a conductivity sigma (c_a, c_b above). Launch the
#    pulse into it and measure how its peak decays. Compare with exp(-alpha x), alpha = sigma Z0 / 2 (low loss).
# 2. **Faraday's disk.** A disk of radius a spins at omega in uniform B along its axis. Compute the EMF as a
#    line integral of (v x B).dl from axis to rim along a straight radius and along a spiral. Compare with omega B a^2/2.
# 3. **Two coils.** Model two coaxial coils as stacks of rings; compute M with `neumann` and with `M_coax`
#    summed over ring pairs; compare with the dipole estimate. Reference solutions follow.

# %%
sigL = 2 * 1.0 / Z0                                         # alpha = 1 per metre
Ns, Sl = 3001, 0.5
xl, El0, Bl0 = gauss_pulse(Ns, 0.3, 0.02, Sl)
sgv = np.where(np.arange(Ns) * 1e-3 > 0.6, sigL, 0.0)
_, _, _, _, snl = yee1d(El0, Bl0, Sl, 4400, sig=sgv, snaps=(1600, 4400))
pk = {n: (xl[np.argmax(snl[n])], snl[n].max()) for n in snl}
alpha_num = -np.log(pk[4400][1] / pk[1600][1]) / (pk[4400][0] - pk[1600][0])
close(alpha_num, 1.0, 0.02, "lossy attenuation")
print(f"[YT1] lossy medium sigma = {sigL:.3e} S/m (sigma/(eps0 c/w) = {sigL/(EPS0*C/0.02):.3f}): peak at x = {pk[1600][0]:.3f} m has {pk[1600][1]:.4f}, at x = {pk[4400][0]:.3f} m has {pk[4400][1]:.4f}: alpha = {alpha_num:.4f} /m (theory sigma Z0/2 = 1)")

a_d, om_d, B_d = 0.10, 3000 * 2 * np.pi / 60, 0.5
u_ = np.linspace(0, 1, 20001)
for name, path in (("radius", np.stack([a_d * u_, 0 * u_], 1)), ("spiral", np.stack([a_d * u_ * np.cos(3 * np.pi * u_), a_d * u_ * np.sin(3 * np.pi * u_)], 1))):
    v = om_d * np.stack([-path[:, 1], path[:, 0], 0 * u_], 1); f = np.cross(v, [0, 0, B_d])[:, :2]
    emf = np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(path, axis=0))
    close(emf, 0.5 * om_d * B_d * a_d**2, 1e-7, "Faraday disk")
    print(f"[YT2] Faraday disk a = 10 cm, 3000 rpm, 0.5 T, path {name}: EMF = {emf:.5f} V (omega B a^2/2 = {0.5*om_d*B_d*a_d**2:.5f} V)")

zc1 = (np.arange(20) + 0.5) * 0.04 / 20 - 0.02; zc2 = 0.10 + (np.arange(30) + 0.5) * 0.03 / 30 - 0.015
M_ell = sum(M_coax(0.03, 0.02, b - a) for a in zc1 for b in zc2)
P1 = np.concatenate([circle_pts(0.03, 128, c=(0, 0, a))[0] for a in zc1]); D1 = np.concatenate([circle_pts(0.03, 128, c=(0, 0, a))[1] for a in zc1])
P2 = np.concatenate([circle_pts(0.02, 128, c=(0, 0, b))[0] for b in zc2]); D2 = np.concatenate([circle_pts(0.02, 128, c=(0, 0, b))[1] for b in zc2])
M_neu = neumann(P1, D1, P2, D2); M_dip = 20 * 30 * MU0 * np.pi * 0.03**2 * 0.02**2 / (2 * 0.10**3)
close(M_neu, M_ell, 1e-9, "two coils")
print(f"[YT3] coils (20 turns, 3 cm, 4 cm long) and (30 turns, 2 cm, 3 cm long), 10 cm apart: M = {M_neu:.5e} H (ring sum {M_ell:.5e}); far-dipole estimate {M_dip:.5e} H")

# %% [markdown]
# ## Part 10: The numbers quoted in solutions.md

# %%
n_al, s_al, A_al = 1.81e29, 3.77e7, np.pi * 1e-6
print(f"[10] A1: R = {10/(s_al*A_al):.5f} Ohm, P = {25*10/(s_al*A_al):.4f} W, v_d = {5/(n_al*QE*A_al):.3e} m/s, tau = {m_e*s_al/(n_al*QE**2):.3e} s, E = {5/(s_al*A_al):.4e} V/m")
print(f"[10] A2: peak EMF = {0.05*0.02*2*np.pi*60:.4f} V, peak I = {0.05*0.02*2*np.pi*60/0.5:.4f} A")
L3 = MU0 * 400**2 * np.pi * 0.015**2 / 0.2; B3 = MU0 * 2000 * 3
close(B3**2 / (2 * MU0) * np.pi * 0.015**2 * 0.2, 0.5 * L3 * 9, 1e-12, "A3 energy")
print(f"[10] A3: volume = {np.pi*0.015**2*0.2:.4e} m^3, L = {L3:.4e} H, W = {0.5*L3*9:.4e} J, B = {B3*1e3:.4f} mT, u = {B3**2/(2*MU0):.3f} J/m^3, cB = {C*B3:.3e} V/m")
C4 = EPS0 * np.pi * 0.03**2 / 2e-3; Id4 = C4 * 50 * 2 * np.pi * 1e3
print(f"[10] A4: C = {C4:.4e} F, peak I_d = {Id4:.4e} A, peak J_d = {Id4/(np.pi*0.03**2):.4e} A/m^2, peak B(rim) = {MU0*Id4/(2*np.pi*0.03):.3e} T")
Trot = maxwell_T(1e6 * np.array([1, 1, 0]) / np.sqrt(2), np.zeros(3)); ev6, V6 = np.linalg.eigh(Trot)
assert np.allclose(Trot[:2, :2], EPS0 * 1e12 / 2 * np.array([[0, 1], [1, 0]])) and abs(abs(V6[:, 2] @ np.array([1, 1, 0])) / np.sqrt(2) - 1) < 1e-12
assert abs(np.trace(maxwell_T(np.array([1e6, 2e6, -3e6]), np.array([0.1, 0.2, 0.3]))) + 0.5 * (EPS0 * 14e12 + 0.14 / MU0)) < 1e-6
print(f"[10] A6: eps0 E^2/2 at 1 MV/m = {EPS0*1e12/2:.4f} Pa; B^2/2mu0 at 1 T = {1/(2*MU0):.4e} Pa; rotated E: eigenvalues {ev6} Pa, top eigenvector {V6[:,2]}")
tq = 10 * B_d * a_d**2 / 2
print(f"[10] B1: omega = {om_d:.3f} rad/s, EMF = {0.5*om_d*B_d*a_d**2:.4f} V, rim speed = {om_d*a_d:.2f} m/s; with I = 10 A: torque = I B a^2/2 = {tq:.4f} N m, power = {tq*om_d:.4f} W = EMF x I = {0.5*om_d*B_d*a_d**2*10:.4f} W")
m2, l2, B2_, R2_ = 5e-3, 0.05, 1.0, 0.01
vt = m2 * G0 * R2_ / (B2_ * l2)**2
s2 = solve_ivp(lambda t, y: [G0 - (B2_ * l2)**2 * y[0] / (m2 * R2_)], (0, 0.1), [0.0], rtol=1e-10, atol=1e-12, dense_output=True)
close(s2.sol(0.05)[0], vt * (1 - np.exp(-G0 * 0.05 / vt)), 1e-7, "falling loop")
print(f"[10] B2: v_t = {vt:.5f} m/s, time constant v_t/g = {vt/G0*1e3:.3f} ms, v(50 ms) = {s2.sol(0.05)[0]:.5f} m/s, I_t = {B2_*l2*vt/R2_:.4f} A, m g v_t = {m2*G0*vt:.5e} W = I^2 R = {(B2_*l2*vt/R2_)**2*R2_:.5e} W")
E3 = MU0 * 1000 * 0.01 * 500 / 2
print(f"[10] B3: mu0 n dI/dt = {MU0*1000*500:.4f} T/s, E(1 cm) = {E3:.4e} V/m, electron acceleration = {QE*E3/m_e:.3e} m/s^2, E(5 cm) = {MU0*1000*0.03**2*500/(2*0.05):.4e} V/m, ring EMF = {MU0*1000*500*np.pi*0.03**2:.4e} V, ring current = {MU0*1000*500*np.pi*0.03**2/0.2*1e3:.4f} mA")
M4 = MU0 * 0.20 / (2 * np.pi) * np.log(3)
close(quad(lambda s: MU0 * 10 / (2 * np.pi * s) * 0.2, 0.05, 0.15)[0] / 10, M4, 1e-10, "B4 flux")
print(f"[10] B4: M = {M4:.4e} H, peak EMF = {M4*2*np.pi*60*10*1e3:.4f} mV")
Lc = MU0 / (2 * np.pi) * np.log(4)
close(2 * quad(lambda s: (MU0 * 1 / (2 * np.pi * s))**2 / (2 * MU0) * 2 * np.pi * s, 0.5e-3, 2e-3)[0], Lc, 1e-10, "coax L via energy")
P5 = quad(lambda s: 100 / (s * np.log(4)) * MU0 * 2 / (2 * np.pi * s) / MU0 * 2 * np.pi * s, 0.5e-3, 2e-3)[0]
close(P5, 200, 1e-10, "coax power")
print(f"[10] B5: L/l = {Lc:.4e} H/m, power through annulus = {P5:.2f} W, S(a) = {100*2/(2*np.pi*(0.5e-3)**2*np.log(4)):.4e} W/m^2, S(b) = {100*2/(2*np.pi*(2e-3)**2*np.log(4)):.4e} W/m^2")
print(f"[10] B6: current sheet K0 = 1 A/m: |E0| = Z0 K0 / 2 = {Z0/2:.2f} V/m; source power <-E.K> = Z0 K0^2/4 = {Z0/4:.2f} W/m^2 = 2 x <S> = {2*(Z0/2)**2/(2*Z0):.2f} W/m^2")
print(f"[10] B7: solenoid of A3: outward pressure on windings B^2/2mu0 = {B3**2/(2*MU0):.3f} Pa")
print(f"[10] B7: plates Q = 1 nC, A = 1 cm^2: E = {1e-9/(EPS0*1e-4):.4e} V/m, force = eps0 E^2 A / 2 = {EPS0*(1e-9/(EPS0*1e-4))**2*1e-4/2:.4e} N = Q^2/(2 eps0 A)")
pc = QE * 0.5 * 0.4 * C / QE / 1e6; mec2 = m_e * C**2 / QE / 1e6; KE1 = np.hypot(pc, mec2) - mec2; Vturn = np.pi * 0.25 * 0.8 / 4e-3
Eg = lambda K: C * np.sqrt(1 - (mec2 / (K + mec2))**2)                                  # speed at kinetic energy K (MeV)
K09 = mec2 * (1 / np.sqrt(1 - 0.81) - 1); t09 = quad(lambda K: np.pi / (Eg(K) * Vturn * 1e-6), 0, K09)[0]
print(f"[10] C1: v = 0.9c at K = {K09:.3f} MeV, reached after {t09*1e6:.1f} us = {t09/4e-3*100:.2f}% of the ramp")
print(f"[10] C1: betatron R = 0.5 m, B(R) = 0.4 T: pc = {pc:.2f} MeV, kinetic energy = {KE1:.2f} MeV, mean B inside orbit = 0.8 T, flux = {np.pi*0.25*0.8:.4f} Wb; 4 ms ramp: {Vturn:.2f} eV per turn, {KE1*1e6/Vturn:.3e} turns, path {KE1*1e6/Vturn*np.pi:.3e} m vs c x 4 ms = {C*4e-3:.3e} m")

# B8: numerical group velocity of a modulated pulse, and the growth rate beyond the Courant limit
k0 = 0.5
for Sb in (0.25, 0.5, 0.9, 1.0):
    xs8 = np.arange(1501) * 1e-3; xh8 = xs8[:-1] + 5e-4; dt8 = Sb * 1e-3 / C
    env = lambda X: np.exp(-((X - 0.3) / 0.04)**2) * np.cos(k0 * (X - 0.3) / 1e-3)
    n1, n2 = round(0.1 / (Sb * 1e-3)), round(0.7 / (Sb * 1e-3))
    sn8 = yee1d(env(xs8), env(xh8 + C * dt8 / 2) / C, Sb, n2, snaps=(n1, n2))[4]
    cx = [np.sum(xs8 * sn8[n]**2) / np.sum(sn8[n]**2) for n in (n1, n2)]
    vg_meas = (cx[1] - cx[0]) / ((n2 - n1) * dt8); vg_th = C * np.cos(k0 / 2) / np.sqrt(1 - Sb**2 * np.sin(k0 / 2)**2)
    close(vg_meas, vg_th, 2e-3, "numerical group velocity")
    print(f"[10] B8: S = {Sb}: measured v/c = {vg_meas/C:.5f}, numerical group velocity v_g/c = {vg_th/C:.5f}")
g400, g600 = np.abs(yee1d(E0, B0, 1.02, 400)[0]).max(), np.abs(yee1d(E0, B0, 1.02, 600)[0]).max()
gf, gth = (g600 / g400)**(1 / 200), np.exp(2 * np.arccosh(1.02))
close(gf, gth, 1e-2, "Courant growth factor")
print(f"[10] B8: S = 1.02: eta = arccosh(S) = {np.arccosh(1.02):.4f}; growth per step {gf:.4f}, predicted exp(2 eta) = {gth:.4f}")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time()-T0:.1f} s. Figures in", FIG)
