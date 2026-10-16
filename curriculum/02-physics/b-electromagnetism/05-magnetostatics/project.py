# %% [markdown]
# # EM-5 Project: Biot-Savart, Ampere, the Vector Potential, and Magnetized Matter
#
# Everything in the EM-5 notes is checked here. Parts: (1) Lorentz force: a Boris pusher
# (the magnetic step is a pure rotation, so it does no work), cyclotron and mass-spectrometer
# numbers; (2) a Biot-Savart integrator for arbitrary polygonal wires using the *exact*
# straight-segment formula, checked on the loop, Helmholtz pair and solenoid; (3) Ampere's law
# as a numerical loop integral, div B = 0 and curl B = mu0 J on grids, and the surface-current
# jump; (4) the vector potential of a solenoid, B = curl A, and the loop integral of A equals
# the flux; (5) the uniformly magnetized sphere, the spinning shell, and the dipole built from
# current loops; (6) field-line pictures including a bar magnet's B and H; (7) a magnetic circuit
# against a 2D finite-difference solve of curl(nu curl A) = J for a C-core; (8) Your turn: the
# cycloid, a toroid, and a Preisach hysteresis loop (reference solutions at the bottom).
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
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.constants import mu_0 as MU0, e as QE, m_p, m_e, atomic_mass as AMU, k as KB, physical_constants

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()
MUB = physical_constants["Bohr magneton"][0]
C4 = MU0 / (4 * np.pi)

def close(a, b, rtol, label=""):
    """assert a ~ b to relative tolerance rtol and return the relative error"""
    err = abs(a - b) / abs(b)
    assert err < rtol, f"{label}: {a} vs {b}, rel err {err:.2e}"
    return err

# %% [markdown]
# ## Part 1: The Lorentz force (notes 3.1-3.2)
# Boris pusher: half electric kick, exact rotation about B, half electric kick. The rotation
# preserves |v| to round-off, which is "magnetic forces do no work" (Theorem 3.1) made numerical.

# %%
def boris(q, m, E, B, x0, v0, dt, n):
    """Integrate m dv/dt = q(E + v x B) (3.1). Returns positions and velocities, shape (n+1, 3)."""
    x, v = np.array(x0, float), np.array(v0, float)
    X, V = [x.copy()], [v.copy()]
    qm = q / m
    for _ in range(n):
        vm = v + 0.5 * dt * qm * E
        t = 0.5 * dt * qm * B
        vp = vm + np.cross(vm + np.cross(vm, t), 2 * t / (1 + t @ t))
        v = vp + 0.5 * dt * qm * E
        x = x + dt * v
        X.append(x.copy()); V.append(v.copy())
    return np.array(X), np.array(V)

Bz = np.array([0, 0, 1.0])
K_J = 1e6 * QE                                    # 1 MeV proton, nonrelativistic (gamma - 1 = 1e-3)
v_p = np.sqrt(2 * K_J / m_p)
r_cyc = m_p * v_p / (QE * 1.0)
w_cyc = QE * 1.0 / m_p
X, V = boris(QE, m_p, np.zeros(3), Bz, [0, 0, 0], [v_p, 0, 0], 2 * np.pi / w_cyc / 200, 400)
close((X[:, 1].max() - X[:, 1].min()) / 2, r_cyc, 1e-3, "cyclotron radius")
close(np.linalg.norm(V, axis=1).max(), np.linalg.norm(V, axis=1).min(), 1e-12, "|v| constant")
print(f"[1] proton 1 MeV in 1 T: v = {v_p:.3e} m/s, r = {r_cyc*100:.2f} cm, omega = {w_cyc:.3e} rad/s, f = {w_cyc/2/np.pi/1e6:.2f} MHz")

def spectrometer_r(m, V_acc, B, q=QE):
    """r = (1/B) sqrt(2 m V / q): accelerate through V (qV = mv^2/2), then circle in B (3.3)."""
    return np.sqrt(2 * m * V_acc / q) / B

r12, r13 = spectrometer_r(12 * AMU, 2000, 0.5), spectrometer_r(13.003355 * AMU, 2000, 0.5)
print(f"[1] spectrometer, 2 kV and 0.5 T: r(12C+) = {r12*100:.3f} cm, r(13C+) = {r13*100:.3f} cm, landing gap 2*dr = {2*(r13-r12)*1e3:.3f} mm")
close(2 * (r13 - r12), 3.655e-3, 1e-3, "spectrometer gap")
print(f"[1] crossed fields E = 1 kV/m, B = 10 mT: drift E/B = {1e3/1e-2:.1e} m/s for any charge or mass")

# %% [markdown]
# ## Part 2: A Biot-Savart integrator (notes 3.4-3.6)
# A straight segment from a to b carries current I. With l = unit vector along it, u = (r-a).l,
# d = perpendicular distance, L = length, (3.9) integrates exactly:
# B = (mu0 I / 4pi) (l x (r-a)) / d^2 * [ (L-u)/sqrt(d^2+(L-u)^2) + u/sqrt(d^2+u^2) ].
# Any wire is a polygon; a closed polygon is a real circuit, so div B = 0 holds to round-off.

# %%
def _seg(P, a, b):
    L = b - a
    Ln = np.linalg.norm(L, axis=1)
    lh = L / Ln[:, None]
    w = P[:, None, :] - a[None]
    u = np.einsum("psk,sk->ps", w, lh)
    d = np.linalg.norm(w - u[..., None] * lh[None], axis=2)
    return lh, Ln, w, u, d

def biot_savart(P, a, b, I, chunk=48):
    """B(P) = (mu0/4pi) sum_segments I * (exact straight-segment integral of dl x r_hat / r^2), eq. (3.9)."""
    P = np.atleast_2d(P).astype(float)
    out = np.empty((len(P), 3))
    for k in range(0, len(P), chunk):
        lh, Ln, w, u, d = _seg(P[k:k + chunk], a, b)
        f = ((Ln - u) / np.sqrt(d ** 2 + (Ln - u) ** 2) + u / np.sqrt(d ** 2 + u ** 2)) / d ** 2
        out[k:k + chunk] = C4 * np.einsum("ps,psk,s->pk", f, np.cross(lh[None], w), I)
    return out

def vector_potential(P, a, b, I, chunk=48):
    """A(P) = (mu0/4pi) sum I l asinh-terms, the exact segment integral of dl'/r, eq. (3.30)."""
    P = np.atleast_2d(P).astype(float)
    out = np.empty((len(P), 3))
    for k in range(0, len(P), chunk):
        lh, Ln, w, u, d = _seg(P[k:k + chunk], a, b)
        out[k:k + chunk] = C4 * np.einsum("ps,sk,s->pk", np.arcsinh((Ln - u) / d) + np.arcsinh(u / d), lh, I)
    return out

def ring(R, z0=0.0, n=96, I=1.0, c=(0.0, 0.0)):
    """Closed n-gon of radius R in the plane z = z0, counterclockwise seen from +z."""
    ph = np.linspace(0, 2 * np.pi, n + 1)
    p = np.stack([c[0] + R * np.cos(ph), c[1] + R * np.sin(ph), z0 + 0 * ph], 1)
    return p[:-1], p[1:], np.full(n, float(I))

def join(parts):
    a, b, I = zip(*parts)
    return np.vstack(a), np.vstack(b), np.concatenate(I)

def solenoid(R, Lsol, n_turns, I, nrings=200, nseg=48):
    """n_turns*I distributed over nrings equal rings (the ring-stack model of a winding)."""
    zs = (np.arange(nrings) + 0.5) * Lsol / nrings - Lsol / 2
    return join([ring(R, z, nseg, I * n_turns / nrings) for z in zs])

def B_loop_axis(I, R, z):
    """Closed form (3.10): B_z = mu0 I R^2 / (2 (R^2+z^2)^{3/2})."""
    return MU0 * I * R ** 2 / (2 * (R ** 2 + z ** 2) ** 1.5)

def B_sol_axis(n, I, R, Lsol, z):
    """Closed form (3.12): finite solenoid on axis."""
    zp, zm = z + Lsol / 2, z - Lsol / 2
    return MU0 * n * I / 2 * (zp / np.sqrt(R ** 2 + zp ** 2) - zm / np.sqrt(R ** 2 + zm ** 2))

# loop on axis (Example 4.1 uses R = 0.1, I = 2 A)
R0, I0 = 0.10, 2.0
a, b, I = ring(R0, 0, 360, I0)
zz = np.array([0, 0.5, 1, 2, 5]) * R0
Bn = biot_savart(np.c_[0 * zz, 0 * zz, zz], a, b, I)
assert np.allclose(Bn[:, :2], 0, atol=1e-12)
for z_, bz_ in zip(zz, Bn[:, 2]):
    close(bz_, B_loop_axis(I0, R0, z_), 2e-4, "loop on axis")
print(f"[2] loop R = 10 cm, I = 2 A: B(0) = {Bn[0,2]*1e6:.3f} uT (closed form {B_loop_axis(I0,R0,0)*1e6:.3f} uT)")

# dipole approximation of the loop (Example 4.4)
m0 = I0 * np.pi * R0 ** 2
for k_ in (5, 10):
    ratio = B_loop_axis(I0, R0, k_ * R0) / (MU0 * m0 / (2 * np.pi * (k_ * R0) ** 3))
    print(f"[2] loop / dipole on axis at z = {k_}R: {ratio:.4f}")
close(B_loop_axis(I0, R0, 5 * R0) / (MU0 * m0 / (2 * np.pi * (5 * R0) ** 3)), (1 + 1 / 25) ** -1.5, 1e-12)

# Helmholtz: second derivative of the pair's axial field vanishes iff d = R (symbolic)
ds, Rs, zs_ = sp.symbols("d R z", positive=True)
f_ = (Rs ** 2 + zs_ ** 2) ** sp.Rational(-3, 2)
pair = f_.subs(zs_, zs_ - ds / 2) + f_.subs(zs_, zs_ + ds / 2)
d2 = sp.simplify(sp.diff(pair, zs_, 2).subs(zs_, 0))
assert sp.simplify(d2 - 6 * (ds ** 2 - Rs ** 2) / (Rs ** 2 + ds ** 2 / 4) ** sp.Rational(7, 2)) == 0
assert sp.solve(d2, ds) == [Rs]
ser = sp.series((pair / pair.subs(zs_, 0)).subs(ds, Rs), zs_, 0, 7).removeO()
c4 = sp.nsimplify(sp.simplify(ser.coeff(zs_, 4) * Rs ** 4))
c2 = sp.simplify(ser.coeff(zs_, 2))
assert c2 == 0 and c4 == sp.Rational(-144, 125), (c2, c4)
Nt = 100
B_helm = (4 / 5) ** 1.5 * MU0 * Nt * I0 / R0
print(f"[2] Helmholtz d = R: B(0)/B(z) series 1 + ({c4})(z/R)^4; N = 100, I = 2 A, R = 10 cm gives B0 = {B_helm*1e3:.4f} mT; single loop of 100 turns {MU0*Nt*I0/2/R0*1e3:.4f} mT")

def helm_map(d, n=25, ext=0.4):
    a, b, I = join([ring(R0, -d / 2, 200, 1.0), ring(R0, d / 2, 200, 1.0)])
    xs, zs = np.linspace(0, ext * R0, n), np.linspace(-ext * R0, ext * R0, 2 * n - 1)
    Xg, Zg = np.meshgrid(xs, zs)
    pts = np.c_[Xg.ravel(), 0 * Xg.ravel(), Zg.ravel()]
    B0 = biot_savart([[0, 0, 0]], a, b, I)[0, 2]
    return Xg, Zg, (biot_savart(pts, a, b, I)[:, 2] / B0 - 1).reshape(Xg.shape), B0

Xg, Zg, dev_H, B0_H = helm_map(R0)
_, _, dev_W, _ = helm_map(1.5 * R0)
r_pts = np.hypot(Xg, Zg)
dH, dW = np.abs(dev_H[r_pts <= 0.2 * R0]).max(), np.abs(dev_W[r_pts <= 0.2 * R0]).max()
close(B0_H * Nt * I0, B_helm, 5e-4, "Helmholtz centre")
print(f"[2] max |B_z/B_0 - 1| within 0.2 R of the centre: d = R: {dH:.2e}; d = 1.5 R: {dW:.2e}")
assert dH < 3e-3 and dW > 10 * dH

# finite solenoid (Example 4.2): 2000 turns on 0.5 m, R = 2 cm, 3 A
Rs_, Ls_, Ns_, Is_ = 0.02, 0.5, 2000, 3.0
n_s = Ns_ / Ls_
SOL = solenoid(Rs_, Ls_, Ns_, Is_)
zax = np.array([-0.2, -0.1, 0.0, 0.1, 0.2])
Bax = biot_savart(np.c_[0 * zax, 0 * zax, zax], *SOL)[:, 2]
for z_, b_ in zip(zax, Bax):
    close(b_, B_sol_axis(n_s, Is_, Rs_, Ls_, z_), 1e-3, "solenoid on axis")
B_inf = MU0 * n_s * Is_
close(Bax[2], B_inf, 4e-3, "finite vs infinite")
B_out = biot_savart([[2 * Rs_, 0, 0]], *SOL)[0, 2]
print(f"[2] solenoid n = {n_s:.0f}/m, I = 3 A: mu0 n I = {B_inf*1e3:.3f} mT; centre {Bax[2]*1e3:.3f} mT; B_z at s = 2R, z = 0: {B_out*1e6:.1f} uT ({B_out/Bax[2]*100:.3f}% of inside)")
print(f"    runtime so far {time.time()-T0:.1f} s")

# %% [markdown]
# ## Part 3: Ampere's law, div B, curl B, and the jump at a surface current (notes 3.7-3.9, 3.14)
# The circulation of the *computed* B around closed paths is mu0 times the current threading
# the path (any path, any shape) and zero for paths that do not link the wire. Then the
# differential equations: div B = 0 on a 3D grid next to a loop, and curl B = mu0 J on a 2D grid
# of smeared (Gaussian) line currents, where J is smooth enough to differentiate.

# %%
def circulation(B_fun, centre, radius, e1, e2, n=720):
    """Closed integral of B . dl around the circle centre + radius (cos t e1 + sin t e2); trapezoid is spectral here."""
    t = np.arange(n) * 2 * np.pi / n
    e1, e2 = np.asarray(e1, float), np.asarray(e2, float)
    pts = np.asarray(centre, float) + radius * (np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2)
    dl = radius * (-np.sin(t)[:, None] * e1 + np.cos(t)[:, None] * e2) * (2 * np.pi / n)
    return np.einsum("pk,pk->", B_fun(pts), dl)

WIRE = (np.array([[0, 0, -5000.0]]), np.array([[0, 0, 5000.0]]), np.array([1.0]))     # I = 1 A, "infinite"
Bw = lambda p: biot_savart(p, *WIRE)
ex, ey, ez = np.eye(3)
c1 = circulation(Bw, [0, 0, 0], 0.5, ex, ey)
c2 = circulation(Bw, [0.2, 0.1, 0], 0.5, ex, ey)            # off-centre, still encloses
c3 = circulation(Bw, [2.0, 0, 0], 0.5, ex, ey)              # does not enclose
close(c1, MU0, 1e-7, "Ampere circle"); close(c2, MU0, 1e-7, "Ampere off-centre")
assert abs(c3) < 1e-6 * MU0
RING = ring(1.0, 0, 400, 1.0)
Br = lambda p: biot_savart(p, *RING)
c5 = circulation(Br, [1, 0, 0], 0.3, ez, ex)                # small circle linking the ring's wire
c6 = circulation(Br, [3, 0, 0], 0.3, ez, ex)                # same circle moved away
close(c5, MU0, 1e-6, "linking loop"); assert abs(c6) < 1e-6 * MU0
print(f"[3] circulation / (mu0 I): centred circle {c1/MU0:.9f}, offset circle {c2/MU0:.9f}, non-enclosing {c3/MU0:.1e}; linking a ring {c5/MU0:.9f}, unlinked {c6/MU0:.1e}")

# div B and curl B on a 3D grid near a loop, central differences
h = 1e-3
P = np.array([0.0, 0.0, 0.0]) + rng.uniform(-0.6, 0.6, (30, 3)) * np.array([1, 1, 1])
P = P[np.abs(np.hypot(np.hypot(P[:, 0], P[:, 1]) - 1, P[:, 2])) > 0.25]                 # stay 0.25 R from the wire
def jac(P):
    D = np.empty((len(P), 3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        D[:, :, j] = (biot_savart(P + e, *RING) - biot_savart(P - e, *RING)) / (2 * h)     # D[p, i, j] = d_j B_i
    return D
D = jac(P); Bs = np.linalg.norm(biot_savart(P, *RING), axis=1)
div = D[:, 0, 0] + D[:, 1, 1] + D[:, 2, 2]
curl = np.stack([D[:, 2, 1] - D[:, 1, 2], D[:, 0, 2] - D[:, 2, 0], D[:, 1, 0] - D[:, 0, 1]], 1)
rel_div = np.abs(div).max() / Bs.min(); rel_curl = np.linalg.norm(curl, axis=1).max() / (Bs).min()
assert rel_div < 1e-4 and rel_curl < 1e-4
print(f"[3] ring, {len(P)} random points away from the wire: max|div B| R/|B| = {rel_div:.1e}, max|curl B| R/|B| = {rel_curl:.1e}")

def smeared_B(X, Y, xs, ys, Is, w):
    """B of z-directed Gaussian current filaments, (mu0 I / 2 pi d^2)(1 - e^{-d^2/2w^2}) (-dy, dx); J_z = I e^{-d^2/2w^2}/(2 pi w^2)."""
    Bx, By, Jz = 0 * X, 0 * X, 0 * X
    for x0, y0, I_ in zip(xs, ys, Is):
        dx, dy = X - x0, Y - y0
        d2 = dx ** 2 + dy ** 2 + 1e-300
        g = -np.expm1(-d2 / (2 * w ** 2)) / d2
        Bx += -MU0 * I_ / (2 * np.pi) * g * dy; By += MU0 * I_ / (2 * np.pi) * g * dx
        Jz += I_ * np.exp(-d2 / (2 * w ** 2)) / (2 * np.pi * w ** 2)
    return Bx, By, Jz

hx = 0.01
xg = np.arange(-1.0, 1.0 + hx / 2, hx)
Xg2, Yg2 = np.meshgrid(xg, xg)
fx, fy, fI = rng.uniform(-0.5, 0.5, 6), rng.uniform(-0.5, 0.5, 6), rng.uniform(-3, 3, 6)
Bx, By, Jz = smeared_B(Xg2, Yg2, fx, fy, fI, 0.2)
curl_z = np.gradient(By, hx, axis=1) - np.gradient(Bx, hx, axis=0)
div2 = np.gradient(Bx, hx, axis=1) + np.gradient(By, hx, axis=0)
inn = (slice(2, -2), slice(2, -2))
err_curl = np.abs(curl_z - MU0 * Jz)[inn].max() / np.abs(MU0 * Jz).max()
err_div = np.abs(div2)[inn].max() / np.abs(MU0 * Jz).max()
assert err_curl < 2e-3 and err_div < 1e-3
print(f"[3] six Gaussian filaments (w = 0.2): max|curl_z B - mu0 J_z| / max|mu0 J| = {err_curl:.1e}; max|div B| / same = {err_div:.1e}")

# jump of the parallel field across the solenoid wall: B_in - B_out = mu0 K, K = n I (3.31)
Bin = biot_savart([[0.8 * Rs_, 0, 0]], *SOL)[0, 2]; Bout = biot_savart([[1.2 * Rs_, 0, 0]], *SOL)[0, 2]
close(Bin - Bout, MU0 * n_s * Is_, 5e-3, "surface current jump")
print(f"[3] solenoid wall: B_in - B_out = {(Bin-Bout)*1e3:.3f} mT vs mu0 K = {MU0*n_s*Is_*1e3:.3f} mT")

# %% [markdown]
# ## Part 4: The vector potential (notes 3.10-3.13)
# A of a polygon is the exact segment integral of dl'/r. Checks: curl A = B by finite differences;
# the loop integral of A around a circle *outside* the winding equals the flux through the disk
# (Stokes), although B = 0 there to 0.3%; and A -> A + grad(lambda) leaves that integral unchanged.

# %%
def curlA(P, segs, hh=2e-5):
    out = np.empty((len(P), 3))
    J = np.empty((len(P), 3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = hh
        J[:, :, j] = (vector_potential(P + e, *segs) - vector_potential(P - e, *segs)) / (2 * hh)
    return np.stack([J[:, 2, 1] - J[:, 1, 2], J[:, 0, 2] - J[:, 2, 0], J[:, 1, 0] - J[:, 0, 1]], 1)

Pt = np.array([[0.005, 0.002, 0.03], [0.015, 0.0, -0.1], [0.03, 0.01, 0.0], [0.0, 0.0, 0.0]])
Bfd, Bbs = curlA(Pt, SOL), biot_savart(Pt, *SOL)
assert np.abs(Bfd - Bbs).max() < 2e-4 * np.abs(Bbs).max()
print(f"[4] curl A (FD) vs Biot-Savart B at 4 points: max difference {np.abs(Bfd-Bbs).max()/np.abs(Bbs).max():.1e} of max|B|")

gl_x, gl_w = np.polynomial.legendre.leggauss(80)
def flux_disk(rho, segs):
    """Flux of B_z through the midplane disk of radius rho > R (exact B_z(s) from Biot-Savart, Gauss-Legendre on [0,R] and [R,rho])."""
    tot = 0.0
    for lo, hi in ((0, Rs_), (Rs_, rho)):
        s = 0.5 * (hi - lo) * gl_x + 0.5 * (hi + lo)
        tot += np.sum(gl_w * 0.5 * (hi - lo) * 2 * np.pi * s * biot_savart(np.c_[s, 0 * s, 0 * s], *segs)[:, 2])
    return tot

phi_inf = MU0 * n_s * Is_ * np.pi * Rs_ ** 2
for rho in (1.5 * Rs_, 3 * Rs_):
    A_circ = circulation(lambda p: vector_potential(p, *SOL), [0, 0, 0], rho, ex, ey, n=96)
    Phi = flux_disk(rho, SOL)
    close(A_circ, Phi, 2e-3, "oint A.dl = flux")
    print(f"[4] circle s = {rho/Rs_:.1f} R in the midplane: oint A.dl = {A_circ:.4e} Wb, flux of B = {Phi:.4e} Wb, ratio to mu0 n I pi R^2 = {A_circ/phi_inf:.4f}; |A_phi| = {A_circ/(2*np.pi*rho)*1e6:.2f} uT m")
lam_grad = lambda p: np.c_[p[:, 1] * p[:, 2], p[:, 0] * p[:, 2], p[:, 0] * p[:, 1]]       # grad(x y z)
assert abs(circulation(lam_grad, [0, 0.01, 0.02], 0.03, ex, ey)) < 1e-12
print("[4] gauge: oint grad(xyz) . dl = 0 around a circle, so oint A.dl is unchanged by A -> A + grad(lambda).")
# infinite solenoid: A_phi = mu0 n I s/2 inside, mu0 n I R^2/(2 s) outside; flux quantum for an electron, Phi = h/e
h_planck = physical_constants["Planck constant"][0]
print(f"[4] h/e = {h_planck/QE:.4e} Wb; flux of the example solenoid = {phi_inf:.3e} Wb = {phi_inf/(h_planck/QE):.3e} flux quanta")
print(f"    runtime so far {time.time()-T0:.1f} s")

# %% [markdown]
# ## Part 5: Magnetized matter built from currents (notes 3.15-3.19, Examples 4.3-4.5, 4.7)
# A uniformly magnetized sphere has K_b = M sin(theta) phi_hat and nothing else (3.40). Build it
# from 200 rings carrying dI = M R sin(theta) dtheta, and check B_in = (2/3) mu0 M, H_in = -M/3,
# the exterior dipole, continuity of B_r at the surface, and A of the spinning shell (3.28).

# %%
def sphere_rings(R, M, nth=200, nseg=64):
    """Rings at polar angles theta_j (midpoints) with current K R dtheta, K = M sin(theta): the surface current M x n_hat."""
    th = (np.arange(nth) + 0.5) * np.pi / nth
    return join([ring(R * np.sin(t), R * np.cos(t), nseg, M * R * np.sin(t) * np.pi / nth) for t in th])

Rm, Mm = 0.01, 1.0e6                                        # 1 cm sphere, M = 1e6 A/m (NdFeB-like)
SPH = sphere_rings(Rm, Mm)
pin = np.array([[0, 0, 0], [0.3 * Rm, 0, 0], [0, 0.2 * Rm, 0.4 * Rm], [0, 0, -0.5 * Rm]])
Bin_s = biot_savart(pin, *SPH)
assert np.allclose(Bin_s[:, :2], 0, atol=2e-3 * MU0 * Mm) and np.allclose(Bin_s[:, 2], 2 / 3 * MU0 * Mm, rtol=1e-3)
H_in = Bin_s[0, 2] / MU0 - Mm
close(H_in, -Mm / 3, 2e-3, "H inside sphere")
m_s = 4 / 3 * np.pi * Rm ** 3 * Mm
pout = np.array([[0, 0, 2 * Rm], [2 * Rm, 0, 0], [0, 0, 1.0001 * Rm], [0, 0, -1.0001 * Rm]])
Bout_s = biot_savart(pout, *SPH)
close(Bout_s[0, 2], MU0 * 2 * m_s / (4 * np.pi * (2 * Rm) ** 3), 2e-3, "dipole on axis")
close(Bout_s[1, 2], -MU0 * m_s / (4 * np.pi * (2 * Rm) ** 3), 2e-3, "dipole on equator")
close(Bout_s[2, 2], 2 / 3 * MU0 * Mm, 5e-3, "B_r continuous at the pole")
print(f"[5] sphere R = 1 cm, M = 1e6 A/m: B_in = {Bin_s[0,2]:.4f} T (2/3 mu0 M = {2/3*MU0*Mm:.4f} T), H_in = {H_in:.4e} A/m (-M/3 = {-Mm/3:.4e}); B(2R, axis) = {Bout_s[0,2]:.4f} T = mu0 M/12 = {MU0*Mm/12:.4f} T; just above the pole {Bout_s[2,2]:.4f} T")

# spinning charged shell: M = sigma R omega; A_in = (mu0 sigma R / 3) omega x r  (3.28)
sig, om = 1e-6, 100.0
Msh = sig * Rm * om
SHELL = sphere_rings(Rm, Msh)
pa = np.array([[0.3 * Rm, 0, 0.0], [0.0, 0.5 * Rm, 0.2 * Rm]])
A_num = vector_potential(pa, *SHELL)
A_th = MU0 * sig * Rm / 3 * np.cross([0, 0, om], pa)
assert np.allclose(A_num, A_th, rtol=2e-3, atol=1e-6 * np.abs(A_th).max())
A_out = vector_potential([[0, 2 * Rm, 0]], *SHELL)[0]
A_dip = C4 * np.cross([0, 0, 4 / 3 * np.pi * Rm ** 4 * sig * om], [0, 2 * Rm, 0]) / (2 * Rm) ** 3
assert np.allclose(A_out, A_dip, rtol=2e-3, atol=1e-6 * np.abs(A_dip).max())
E_surf = sig / 8.8541878128e-12
print(f"[5] spinning shell sigma = 1e-6 C/m^2, omega = 100 rad/s: B_in = {2/3*MU0*Msh*1e12:.3f} pT (E at the surface sigma/eps0 = {E_surf:.3e} V/m); A_y(0.3R) = {A_num[0,1]:.4e} vs {A_th[0,1]:.4e} T m")

# uniformly magnetized cylinder = solenoid with nI = M (Example 4.7c): R = 5 mm, L = 20 mm
Rc, Lc = 5e-3, 20e-3
CYL = solenoid(Rc, Lc, 1, Mm * Lc, nrings=100)                # total current M L spread over L: K = M
Bc = biot_savart([[0, 0, 0]], *CYL)[0, 2]
Bc_th = MU0 * Mm * (Lc / 2) / np.sqrt(Rc ** 2 + (Lc / 2) ** 2)
close(Bc, Bc_th, 1e-3, "magnetized cylinder centre")
print(f"[5] magnetized cylinder (R = 5 mm, L = 20 mm, M = 1e6): B(centre) = {Bc:.4f} T = mu0 M L/2/sqrt(R^2+L^2/4) = {Bc_th:.4f} T (long limit mu0 M = {MU0*Mm:.4f} T)")

# materials numbers (Example 4.7)
nI, mur_iron, chi_al = 1000 * 0.1, 2000.0, 2.1e-5                 # n = 1000/m, I = 0.1 A
B_fill = MU0 * mur_iron * nI
close(B_fill, MU0 * nI + MU0 * (mur_iron - 1) * nI, 1e-12)          # B = mu0 (H + M)
print(f"[5] long solenoid nI = {nI:.0f} A/m filled with mu_r = 2000: B = {B_fill:.4f} T (vacuum {MU0*nI*1e3:.4f} mT), M = {(mur_iron-1)*nI:.4e} A/m, K_b = {(mur_iron-1)*nI:.4e} A/m = {mur_iron-1:.0f} x the free K = {nI:.0f} A/m; aluminium (chi_m = 2.1e-5): K_b = {chi_al*nI*1e3:.3f} mA/m")
Nt_, It_, rm_ = 500, 0.5, 0.10
Ht = Nt_ * It_ / (2 * np.pi * rm_); Bt = MU0 * mur_iron * Ht
print(f"[5] toroid N = 500, I = 0.5 A, mean radius 10 cm, mu_r = 2000: H = {Ht:.1f} A/m, B = {Bt:.3f} T, M = {(mur_iron-1)*Ht:.3e} A/m")
print(f"[5] mu_B B / k T at 1 T, 300 K = {MUB*1.0/(KB*300):.2e}; Curie chi = mu0 n mu_B^2/(k T) for n = 6e28 /m^3 at 300 K: {MU0*6e28*MUB**2/(KB*300):.2e}")
mE = 8e22                                                  # Earth's dipole moment, A m^2
print(f"[5] Earth dipole m = 8e22 A m^2, R_E = 6.371e6 m: B_equator = {C4*mE/6.371e6**3*1e6:.1f} uT, B_pole = {2*C4*mE/6.371e6**3*1e6:.1f} uT")

# %% [markdown]
# ## Part 6: Field-line pictures (notes 2, 3.19)
# Streamplots in the x-z plane (y = 0) for a loop, a finite solenoid, and a bar magnet modeled
# as a stack of rings with nI = M. For the bar, B closes on itself, while H = B/mu0 - M inside
# points *against* M: H lines start on the north face and end on the south face.

# %%
def grid_field(segs, xlim, zlim, n=60):
    xs, zs = np.linspace(*xlim, n), np.linspace(*zlim, n)
    Xp, Zp = np.meshgrid(xs, zs)
    Bf = biot_savart(np.c_[Xp.ravel(), 0 * Xp.ravel(), Zp.ravel()], *segs)
    return Xp, Zp, Bf[:, 0].reshape(Xp.shape), Bf[:, 2].reshape(Xp.shape)

fig, axs = plt.subplots(2, 2, figsize=(10, 9))
Rb, Lb = 1.0, 4.0
BAR = solenoid(Rb, Lb, 1, 1.0 * Lb, nrings=60, nseg=48)       # M = 1 (arbitrary units: B in units of mu0 M)
panels = [("single loop (B)", ring(1.0, 0, 96, 1.0), (-3, 3), (-3, 3), None),
          ("solenoid R = 1, L = 6 (B)", solenoid(1.0, 6.0, 1, 6.0, 60, 48), (-3, 3), (-5, 5), None),
          ("bar magnet R = 1, L = 4: B", BAR, (-3, 3), (-5, 5), "B"),
          ("bar magnet: H = B/mu0 - M", BAR, (-3, 3), (-5, 5), "H")]
for ax, (ttl, segs, xl, zl, kind) in zip(axs.ravel(), panels):
    Xp, Zp, Bx_, Bz_ = grid_field(segs, xl, zl)
    if kind == "H":                                         # H = B/mu0 - M, with M = 1 A/m inside the bar
        inside = (np.abs(Xp) < Rb) & (np.abs(Zp) < Lb / 2)
        Bx_, Bz_ = Bx_ / MU0, Bz_ / MU0 - inside * 1.0
    sp_ = np.hypot(Bx_, Bz_)
    ax.streamplot(Xp, Zp, Bx_, Bz_, color=np.log10(sp_ + 1e-30), cmap="viridis", density=1.4, linewidth=1)
    if kind:
        ax.add_patch(plt.Rectangle((-Rb, -Lb / 2), 2 * Rb, Lb, fill=False, ec="r", lw=1.5))
    ax.set_title(ttl); ax.set_aspect("equal"); ax.set_xlim(xl); ax.set_ylim(zl)
fig.tight_layout(); fig.savefig(FIG / "field_lines.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 7: A magnetic circuit versus a 2D finite-difference solve (notes 3.20, Example 4.6)
# 2D, z-independent: A = A_z(x,y) z_hat, B = (dA/dy, -dA/dx), and curl(nu curl A) = J_f with
# nu = 1/mu becomes -div(nu grad A_z) = J_z (3.47). Cell-centred finite volumes, harmonic-mean
# nu on faces (so the normal flux nu dA/dn is continuous across an iron-air interface), A = 0 on the
# box. The C-core: 20 cm square, 4 cm legs, a 4 mm gap across the right leg, coil on the left leg
# carrying NI = 1000 A-turns (current out of the page on one side of the leg, into it on the other).

# %%
def fd_solve(nu, Jz, h):
    """Solve -div(nu grad A) = J on a cell-centred grid; A = 0 at the box edge (ghost face at h/2)."""
    ny, nx = nu.shape
    idx = np.arange(nx * ny).reshape(ny, nx)
    hm = lambda p, q: 2 * p * q / (p + q)
    cx, cy = hm(nu[:, :-1], nu[:, 1:]) / h ** 2, hm(nu[:-1, :], nu[1:, :]) / h ** 2
    diag = np.zeros((ny, nx))
    diag[:, :-1] += cx; diag[:, 1:] += cx; diag[:-1, :] += cy; diag[1:, :] += cy
    for sl in ((slice(None), 0), (slice(None), -1), (0, slice(None)), (-1, slice(None))):
        diag[sl] += 2 * nu[sl] / h ** 2
    rows, cols, vals = [idx.ravel()], [idx.ravel()], [diag.ravel()]
    for i1, i2, c in ((idx[:, :-1], idx[:, 1:], cx), (idx[:-1, :], idx[1:, :], cy)):
        rows += [i1.ravel(), i2.ravel()]; cols += [i2.ravel(), i1.ravel()]; vals += [-c.ravel(), -c.ravel()]
    K = sparse.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nx * ny,) * 2)
    return spsolve(K, Jz.ravel()).reshape(ny, nx)

def ccore(mur, NI, g=4e-3, h=1e-3, box=0.32):
    n = int(round(box / h)); xc = (np.arange(n) + 0.5) * h
    X, Y = np.meshgrid(xc, xc)
    core = (X > .06) & (X < .26) & (Y > .06) & (Y < .26) & ~((X > .10) & (X < .22) & (Y > .10) & (Y < .22))
    core &= ~((X > .22) & (np.abs(Y - .16) < g / 2))
    nu = np.where(core, 1 / (MU0 * mur), 1 / MU0)
    cout = (X > .03) & (X < .06) & (Y > .12) & (Y < .20); cin = (X > .10) & (X < .13) & (Y > .12) & (Y < .20)
    J = np.where(cout, NI / (cout.sum() * h ** 2), 0.0) - np.where(cin, NI / (cin.sum() * h ** 2), 0.0)
    return X, Y, core, nu, J, fd_solve(nu, J, h)

t7 = time.time()
mur_c, NI_c, g_c, h_c = 1000.0, 1000.0, 4e-3, 1e-3
X7, Y7, core7, nu7, J7, A7 = ccore(mur_c, NI_c, g_c, h_c)
By7 = -(A7[:, 2:] - A7[:, :-2]) / (2 * h_c)                      # centred, columns 1..n-2
Bx7 = (A7[2:, :] - A7[:-2, :]) / (2 * h_c)                       # rows 1..n-2
jm, ig = 240, 159                                                # column x = 0.2405 m (centre of right leg), rows y ~ 0.16 m
B_gap_centre = abs(By7[ig:ig + 2, jm - 1]).mean()                # By7 column k holds A column k+1
Phi_gap = abs(A7[ig:ig + 2, 259].mean() - A7[ig:ig + 2, 220].mean())      # flux per unit depth between x = 0.2205 and 0.2595
B_gap_avg = Phi_gap / (0.039)
l_core = 4 * 0.16 - g_c
B_model = MU0 * NI_c / (g_c + l_core / mur_c)
frac_gap = g_c / (g_c + l_core / mur_c)
# Ampere's law for H around the centre-line of the core (columns 80 and 240, rows 80 and 240), 1/mu0 outside iron
H_x = (A7[2:, :] - A7[:-2, :]) / (2 * h_c) * nu7[1:-1, :]        # H_x = nu B_x, rows 1..n-2
H_y = -(A7[:, 2:] - A7[:, :-2]) / (2 * h_c) * nu7[:, 1:-1]       # H_y = nu B_y, columns 1..n-2
left = H_y[80:240, 79].sum() * h_c                               # up the left leg  (column 80 -> index 79 in H_y)
top = H_x[238, 80:240].sum() * h_c                               # along the top    (row 240 -> index 238 in H_x)
right = -H_y[80:240, 239].sum() * h_c                            # down the right leg through the gap
bottom = -H_x[78, 80:240].sum() * h_c                            # back along the bottom
circ = left + top + right + bottom
close(circ, NI_c, 2e-2, "Ampere for H around the core")
print(f"[7] C-core mu_r = 1000, g = 4 mm, NI = 1000 A, h = 1 mm ({core7.shape[0]**2} cells, {time.time()-t7:.1f} s): oint H.dl = {circ:.1f} A (NI = 1000)")
print(f"[7] reluctance model: B = mu0 NI / (g + l/mu_r) = {B_model:.4f} T with l = {l_core:.3f} m ({frac_gap*100:.1f}% of the mmf drops across the gap)")
print(f"    2D solve: B_gap(centre) = {B_gap_centre:.4f} T, B_gap(flux / 39 mm) = {B_gap_avg:.4f} T; ratios to model {B_gap_centre/B_model:.3f}, {B_gap_avg/B_model:.3f}")
assert abs(B_gap_avg / B_model - 1) < 0.04 and abs(B_gap_centre / B_model - 1) < 0.04
# without the gap the same NI would push linear iron to mu0 mu_r NI / l, already at saturation; the gap makes B ~ mu0 NI / g
print(f"[7] same core with g = 0 (linear iron): B = mu0 mu_r NI / l = {MU0*mur_c*NI_c/0.64:.2f} T; upper bound mu0 NI / g = {MU0*NI_c/g_c:.4f} T; mu_r = 5000 gives {MU0*NI_c/(g_c+l_core/5000):.4f} T")
close(MU0 * NI_c / (g_c + l_core / 5000), MU0 * NI_c / g_c, 0.04, "gap dominates")

fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
Bm = np.hypot(By7[1:-1, :], Bx7[:, 1:-1])
im = axs[0].imshow(Bm, origin="lower", extent=[0, 0.32, 0, 0.32], cmap="magma", vmax=0.5)
axs[0].contour(X7, Y7, core7, [0.5], colors="c", linewidths=.6); fig.colorbar(im, ax=axs[0], label="|B| (T)")
axs[0].set_title("C-core |B|: gap field vs leakage"); axs[0].set_xlabel("x (m)"); axs[0].set_ylabel("y (m)")
axs[1].contour(X7, Y7, A7, 30, cmap="coolwarm"); axs[1].contour(X7, Y7, core7, [0.5], colors="k", linewidths=.6)
axs[1].set_title("flux lines = contours of A_z"); axs[1].set_aspect("equal")
fig.tight_layout(); fig.savefig(FIG / "c_core.png", dpi=100); plt.close(fig)
print(f"    runtime so far {time.time()-T0:.1f} s")

# %% [markdown]
# ## Your turn
# 1. **The cycloid.** Charge q in E = E y_hat, B = B z_hat, starting from rest: derive
#    x = (E / omega B)(omega t - sin omega t), y = (E / omega B)(1 - cos omega t) (notes 3.2), and
#    check the Boris integrator against it. Then start with a nonzero velocity and watch the
#    trochoid. 2. **A toroid.** B_phi = mu0 N I / (2 pi s) varies across the cross-section; compute
#    the flux for a rectangular cross-section and the error of the "mean radius" estimate used in
#    Example 4.7. 3. **Hysteresis.** A Preisach model: each hysteron flips up at H > up_i and
#    down at H < dn_i; the loop is the ensemble average. Change the distribution of (hu, hc) and
#    describe what makes a "soft" versus a "hard" material. Write your own first; references follow.

# %%
# 1. cycloid
Eyy, Bzz = 1e3, 1e-2
wc = QE * Bzz / m_p
Tc = 2 * np.pi / wc
Xc, Vc = boris(QE, m_p, np.array([0, Eyy, 0.0]), np.array([0, 0, Bzz]), [0, 0, 0], [0, 0, 0], Tc / 400, 1200)
tt = np.arange(1201) * Tc / 400
x_an = Eyy / (wc * Bzz) * (wc * tt - np.sin(wc * tt)); y_an = Eyy / (wc * Bzz) * (1 - np.cos(wc * tt))
arch = 2 * Eyy / (wc * Bzz)
cyc_err = max(np.abs(Xc[:, 0] - x_an).max(), np.abs(Xc[:, 1] - y_an).max()) / arch      # O(omega dt) start-up phase of leapfrog; halve dt to halve it
assert cyc_err < 1e-2
close(Xc[-1, 0] / tt[-1], Eyy / Bzz, 1e-2, "drift speed E/B")
print(f"[8] cycloid (proton, E = 1 kV/m, B = 10 mT): drift {Eyy/Bzz:.2e} m/s, arch height 2E/(omega B) = {arch*1e3:.3f} mm, period {Tc*1e6:.3f} us")

# 2. toroid with a rectangular cross-section a < s < b, height hz
a_t, b_t, hz_t, N_t, I_t = 0.08, 0.12, 0.04, 500, 0.5
Phi_t = MU0 * N_t * I_t * hz_t * np.log(b_t / a_t) / (2 * np.pi)
Phi_mean = MU0 * N_t * I_t / (2 * np.pi * (a_t + b_t) / 2) * hz_t * (b_t - a_t)
close(Phi_t, MU0 * N_t * I_t * hz_t / (2 * np.pi) * np.log(1.5), 1e-12)
print(f"[8] toroid a = 8 cm, b = 12 cm: B(a)/B(b) = {b_t/a_t:.2f}, B(a) = {MU0*N_t*I_t/(2*np.pi*a_t)*1e3:.3f} mT, B(b) = {MU0*N_t*I_t/(2*np.pi*b_t)*1e3:.3f} mT; "
      f"flux {Phi_t:.4e} Wb vs mean-radius estimate {Phi_mean:.4e} Wb ({(Phi_mean/Phi_t-1)*100:.2f}% high)")

# 3. Preisach hysteresis
Nh = 20000
hu, hc = rng.normal(0, 0.8, Nh), np.abs(rng.normal(0.5, 0.3, Nh))
up, dn = hu + hc, hu - hc
sgn = rng.choice([-1.0, 1.0], Nh)
Hs = np.concatenate([np.linspace(0, 3, 150), np.linspace(3, -3, 300), np.linspace(-3, 3, 300)])
Ms = []
for Hh in Hs:
    sgn[Hh > up] = 1.0; sgn[Hh < dn] = -1.0
    Ms.append(sgn.mean())
Ms = np.array(Ms)
Hd, Md = Hs[150:450], Ms[150:450]                       # descending branch
Hu, Mu = Hs[450:], Ms[450:]                              # ascending branch
M_rem = np.interp(0, Hd[::-1], Md[::-1])
H_coer = -np.interp(0, Md[::-1], Hd[::-1])
area = abs(np.sum(0.5 * (Ms[151:] + Ms[150:-1]) * np.diff(Hs[150:])))     # |oint M dH| over the closed loop
assert 0 < M_rem < 1 and H_coer > 0 and area > 0
assert abs(np.interp(1.0, Hu, Mu) + np.interp(-1.0, Hd[::-1], Md[::-1])) < 0.03       # symmetry M_up(H) = -M_down(-H)
print(f"[8] Preisach loop: remanence {M_rem:.3f} M_s, coercivity {H_coer:.3f}, loop area {area:.3f} (energy lost per cycle per volume = mu0 M_s H_s x area)")

fig, axs = plt.subplots(1, 3, figsize=(13, 4))
axs[0].plot(Xc[:, 0] * 1e3, Xc[:, 1] * 1e3, "C0", label="Boris"); axs[0].plot(x_an * 1e3, y_an * 1e3, "k--", lw=.8, label="analytic")
axs[0].set_title("cycloid in crossed E, B"); axs[0].set_xlabel("x (mm)"); axs[0].set_ylabel("y (mm)"); axs[0].legend()
ss = np.linspace(a_t, b_t, 100); axs[1].plot(ss * 100, MU0 * N_t * I_t / (2 * np.pi * ss) * 1e3)
axs[1].set_title("toroid B(s) over the cross-section"); axs[1].set_xlabel("s (cm)"); axs[1].set_ylabel("B (mT)")
axs[2].plot(Hs[:150], Ms[:150], "C1", label="virgin curve"); axs[2].plot(Hs[150:], Ms[150:], "C0", label="loop")
axs[2].set_title("Preisach hysteresis"); axs[2].set_xlabel("H (arb.)"); axs[2].set_ylabel("M / M_s"); axs[2].legend(); axs[2].grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "your_turn.png", dpi=100); plt.close(fig)

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time()-T0:.1f} s. Figures in", FIG)
