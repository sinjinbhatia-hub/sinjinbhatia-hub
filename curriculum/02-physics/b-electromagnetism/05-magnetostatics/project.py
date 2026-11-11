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
# cycloid, a toroid, and a Preisach hysteresis loop (reference solutions follow); (9) every
# numerical answer quoted in solutions.md, including the extension asked for in Exercise B8.
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
from scipy.optimize import brentq
from scipy.constants import mu_0 as MU0, epsilon_0 as EPS0, c as C_LIGHT, e as QE, m_p, m_e, atomic_mass as AMU, k as KB, physical_constants

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
# ## Part 2: A Biot-Savart integrator (notes 3.3-3.4)
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

def biot_savart_midpoint(P, a, b, I):
    """Naive (3.9): one term mu0 I dl x r_hat / (4 pi r^2) per segment, evaluated at the segment midpoint."""
    P = np.atleast_2d(P).astype(float)
    dl, mid = b - a, 0.5 * (a + b)
    w = P[:, None, :] - mid[None]
    return C4 * np.einsum("psk,s,ps->pk", np.cross(dl[None], w), I, np.linalg.norm(w, axis=2) ** -3)

def vector_potential(P, a, b, I, chunk=48):
    """A(P) = (mu0/4pi) sum I l (asinh terms): the exact segment integral of dl'/r in eq. (3.24)."""
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
    """Closed form (3.13): B_z = mu0 I R^2 / (2 (R^2+z^2)^{3/2})."""
    return MU0 * I * R ** 2 / (2 * (R ** 2 + z ** 2) ** 1.5)

def B_sol_axis(n, I, R, Lsol, z):
    """Closed form (3.14): finite solenoid on axis."""
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

# "build it yourself": the midpoint rule converges as 1/N^2 for a closed polygon
Pq = np.array([[0.04, 0.02, 0.05]])
ref_q = biot_savart(Pq, *ring(R0, 0, 4000, I0))[0]
err_mid = [np.linalg.norm(biot_savart_midpoint(Pq, *ring(R0, 0, N_, I0))[0] - ref_q) / np.linalg.norm(ref_q) for N_ in (16, 32, 64)]
assert all(3.8 < err_mid[i] / err_mid[i + 1] < 4.2 for i in range(2))
print(f"[2] midpoint rule, N = 16, 32, 64 segments: errors {err_mid[0]:.2e}, {err_mid[1]:.2e}, {err_mid[2]:.2e} (ratios {err_mid[0]/err_mid[1]:.2f}, {err_mid[1]/err_mid[2]:.2f}: 1/N^2)")

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
print(f"[2] series (3.16) at z = 0.2 R on the axis: -{144/125*0.2**4:.2e} B_0")
fig, axs = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
for ax, dv, ttl in zip(axs, (dev_H, dev_W), ("d = R", "d = 1.5 R")):
    cs = ax.contourf(Xg / R0, Zg / R0, dv, levels=np.linspace(-0.03, 0.03, 13), cmap="coolwarm", extend="both")
    ax.set_title(f"Helmholtz pair, {ttl}: B_z/B_0 - 1"); ax.set_xlabel("x / R"); ax.set_ylabel("z / R")
fig.colorbar(cs, ax=axs); fig.savefig(FIG / "helmholtz.png", dpi=100); plt.close(fig)

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
# ## Part 3: Ampere's law, div B, curl B, and the jump at a surface current (notes 3.5, 3.6, 3.9)
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

# jump of the parallel field across the solenoid wall: B_in - B_out = mu0 K, K = n I (3.28)
Bin = biot_savart([[0.8 * Rs_, 0, 0]], *SOL)[0, 2]; Bout = biot_savart([[1.2 * Rs_, 0, 0]], *SOL)[0, 2]
close(Bin - Bout, MU0 * n_s * Is_, 5e-3, "surface current jump")
print(f"[3] solenoid wall: B_in - B_out = {(Bin-Bout)*1e3:.3f} mT vs mu0 K = {MU0*n_s*Is_*1e3:.3f} mT")

# %% [markdown]
# ## Part 4: The vector potential (notes 3.7-3.8)
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
    A_inf = MU0 * n_s * Is_ * Rs_ ** 2 / (2 * rho)                  # infinite solenoid, outside (3.26)
    close(A_circ / (2 * np.pi * rho), A_inf, 4e-2, "A_phi vs infinite solenoid")
    print(f"[4] circle s = {rho/Rs_:.1f} R in the midplane: oint A.dl = {A_circ:.4e} Wb, flux of B = {Phi:.4e} Wb, ratio to mu0 n I pi R^2 = {A_circ/phi_inf:.4f}; |A_phi| = {A_circ/(2*np.pi*rho)*1e6:.2f} uT m (infinite solenoid (3.26): {A_inf*1e6:.2f})")
lam_grad = lambda p: np.c_[p[:, 1] * p[:, 2], p[:, 0] * p[:, 2], p[:, 0] * p[:, 1]]       # grad(x y z)
assert abs(circulation(lam_grad, [0, 0.01, 0.02], 0.03, ex, ey)) < 1e-12
print("[4] gauge: oint grad(xyz) . dl = 0 around a circle, so oint A.dl is unchanged by A -> A + grad(lambda).")
# infinite solenoid: A_phi = mu0 n I s/2 inside, mu0 n I R^2/(2 s) outside; flux quantum for an electron, Phi = h/e
h_planck = physical_constants["Planck constant"][0]
print(f"[4] h/e = {h_planck/QE:.4e} Wb; flux of the example solenoid = {phi_inf:.3e} Wb = {phi_inf/(h_planck/QE):.3e} flux quanta")
print(f"    runtime so far {time.time()-T0:.1f} s")

# %% [markdown]
# ## Part 5: Magnetized matter built from currents (notes 3.8, 3.12-3.15, Examples 4.3-4.5, 4.7)
# A uniformly magnetized sphere has K_b = M sin(theta) phi_hat and nothing else (3.40). Build it
# from 200 rings carrying dI = M R sin(theta) dtheta, and check B_in = (2/3) mu0 M, H_in = -M/3,
# the exterior dipole, continuity of B_r at the surface, and A of the spinning shell (3.27).

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

H_above = Bout_s[2, 2] / MU0
close(H_above - H_in, Mm, 5e-3, "H_perp jumps by M")
print(f"[5] H_perp at the pole: inside {H_in:.3e}, just above {H_above:.3e} A/m, jump {H_above-H_in:.3e} = M")

# spinning charged shell: M = sigma R omega; A_in = (mu0 sigma R / 3) omega x r  (3.27)
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
E_surf = sig / EPS0
print(f"[5] spinning shell sigma = 1e-6 C/m^2, omega = 100 rad/s: B_in = {2/3*MU0*Msh*1e12:.3f} pT (E at the surface sigma/eps0 = {E_surf:.3e} V/m); A_y(0.3R) = {A_num[0,1]:.4e} vs {A_th[0,1]:.4e} T m")

# uniformly magnetized cylinder = solenoid with nI = M (Example 4.7c): R = 5 mm, L = 20 mm
Rc, Lc = 5e-3, 20e-3
CYL = solenoid(Rc, Lc, 1, Mm * Lc, nrings=100)                # total current M L spread over L: K = M
Bc = biot_savart([[0, 0, 0]], *CYL)[0, 2]
Bc_th = MU0 * Mm * (Lc / 2) / np.sqrt(Rc ** 2 + (Lc / 2) ** 2)
close(Bc, Bc_th, 1e-3, "magnetized cylinder centre")
print(f"[5] magnetized cylinder (R = 5 mm, L = 20 mm, M = 1e6): B(centre) = {Bc:.4f} T = mu0 M L/2/sqrt(R^2+L^2/4) = {Bc_th:.4f} T (long limit mu0 M = {MU0*Mm:.4f} T)")
Be = biot_savart([[0, 0, Lc / 2]], *CYL)[0, 2]
Be_th = 0.5 * MU0 * Mm * Lc / np.sqrt(Rc ** 2 + Lc ** 2)
close(Be, Be_th, 1e-3, "end-face centre")
print(f"[5] same magnet, centre of an end face: B = {Be:.4f} T (formula {Be_th:.4f}); H just inside {Be/MU0-Mm:.4e}, just outside {Be/MU0:.4e} A/m")

# materials numbers (Example 4.7)
nI, mur_iron, chi_al = 1000 * 0.1, 2000.0, 2.1e-5                 # n = 1000/m, I = 0.1 A
B_fill = MU0 * mur_iron * nI
close(B_fill, MU0 * nI + MU0 * (mur_iron - 1) * nI, 1e-12)          # B = mu0 (H + M)
print(f"[5] long solenoid nI = {nI:.0f} A/m filled with mu_r = 2000: B = {B_fill:.4f} T (vacuum {MU0*nI*1e3:.4f} mT), M = {(mur_iron-1)*nI:.4e} A/m, K_b = {(mur_iron-1)*nI:.4e} A/m = {mur_iron-1:.0f} x the free K = {nI:.0f} A/m; aluminium (chi_m = 2.1e-5): K_b = {chi_al*nI*1e3:.3f} mA/m")
Nt_, It_, rm_ = 500, 0.5, 0.10
Ht = Nt_ * It_ / (2 * np.pi * rm_); Bt = MU0 * mur_iron * Ht
print(f"[5] toroid N = 500, I = 0.5 A, mean radius 10 cm, mu_r = 2000: H = {Ht:.1f} A/m, B = {Bt:.3f} T, M = {(mur_iron-1)*Ht:.3e} A/m")
chi_dia = -MU0 * 6e28 * 10 * QE ** 2 * (1e-10) ** 2 / (6 * m_e)
print(f"[5] orbital diamagnetism (3.36): chi_m = -mu0 n Z e^2 <r^2> / 6 m_e = {chi_dia:.2e} (n = 6e28, Z = 10, <r^2>^(1/2) = 1 A)")
print(f"[5] mu_B B / k T at 1 T, 300 K = {MUB*1.0/(KB*300):.2e}; Curie chi = mu0 n mu_B^2/(k T) for n = 6e28 /m^3 at 300 K: {MU0*6e28*MUB**2/(KB*300):.2e}")
mE = 8e22                                                  # Earth's dipole moment, A m^2
print(f"[5] Earth dipole m = 8e22 A m^2, R_E = 6.371e6 m: B_equator = {C4*mE/6.371e6**3*1e6:.1f} uT, B_pole = {2*C4*mE/6.371e6**3*1e6:.1f} uT")

# %% [markdown]
# ## Part 6: Field-line pictures (notes 2, 3.15)
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
# ## Part 7: A magnetic circuit versus a 2D finite-difference solve (notes 3.17, Example 4.6)
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

def magnetic_circuit(NI, segs):
    """Series circuit (3.45): segs = [(length, area, mu_r), ...]. Returns flux Phi = NI / sum(R_k), the B in each segment, and R_k = l / (mu0 mu_r A)."""
    R_k = np.array([l / (MU0 * mr * A) for l, A, mr in segs])
    Phi = NI / R_k.sum()
    return Phi, np.array([Phi / A for _, A, _ in segs]), R_k

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
A_leg = 0.04                                                     # leg width x unit depth
_, (B_model, _), R_seg = magnetic_circuit(NI_c, [(l_core, A_leg, mur_c), (g_c, A_leg, 1.0)])
close(B_model, MU0 * NI_c / (g_c + l_core / mur_c), 1e-12, "(3.45) vs (3.46)")
frac_gap = R_seg[1] / R_seg.sum()
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
print(f"[7] mmf bookkeeping: H_core = {B_model/(MU0*mur_c):.1f} A/m, H_gap = {B_model/MU0:.4e} A/m; H_core l + H_gap g = {B_model/(MU0*mur_c)*l_core:.1f} + {B_model/MU0*g_c:.1f} = {B_model/(MU0*mur_c)*l_core + B_model/MU0*g_c:.1f} A")
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
#    x = (E / omega B)(omega t - sin omega t), y = (E / omega B)(1 - cos omega t) (notes 3.1), and
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
      f"flux {Phi_t:.4e} Wb vs mean-radius estimate {Phi_mean:.4e} Wb ({(Phi_mean/Phi_t-1)*100:+.2f}% relative to the exact value)")

# 3. Preisach hysteresis
def preisach_ensemble(rng_, Nh=20000):
    """Hysterons flip up at H > hu + hc and down at H < hu - hc (hc > 0); returns up, down thresholds and random initial states."""
    hu, hc = rng_.normal(0, 0.8, Nh), np.abs(rng_.normal(0.5, 0.3, Nh))
    return hu + hc, hu - hc, rng_.choice([-1.0, 1.0], Nh)

def preisach_run(Hseq, up, dn, state):
    """Apply the field sequence Hseq; returns the mean state M/M_s after each step (state is updated in place)."""
    out = []
    for Hh in Hseq:
        state[Hh > up] = 1.0; state[Hh < dn] = -1.0
        out.append(state.mean())
    return np.array(out)

up, dn, sgn = preisach_ensemble(rng)
sgn_init = sgn.copy()
Hs = np.concatenate([np.linspace(0, 3, 150), np.linspace(3, -3, 300), np.linspace(-3, 3, 300)])
Ms = preisach_run(Hs, up, dn, sgn)
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

# %% [markdown]
# ## Part 9: The numbers quoted in solutions.md
# One block per exercise group. Theory checks are asserted; the rest are printed so that
# solutions.md can quote them. B8 is the extension of this file that Exercise B8 asks you to write.

# %%
def small_loop_force_torque(src, c, normal, side, I_s):
    """F = I_s sum dl x B(mid) and N = I_s sum (mid - c) x (dl x B(mid)) on a small square loop (3.6); also returns m = I_s side^2 n_hat."""
    n_ = np.asarray(normal, float); n_ = n_ / np.linalg.norm(n_)
    e1 = np.cross(n_, [0, 1.0, 0]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(n_, e1)                                          # e1 x e2 = n_: counterclockwise about n_
    c = np.asarray(c, float); hs = side / 2
    P4 = np.array([c - hs * e1 - hs * e2, c + hs * e1 - hs * e2, c + hs * e1 + hs * e2, c - hs * e1 + hs * e2])
    dl = np.roll(P4, -1, axis=0) - P4
    mid = P4 + dl / 2
    f = np.cross(dl, biot_savart(mid, *src))
    return I_s * f.sum(0), I_s * np.cross(mid - c, f).sum(0), I_s * side ** 2 * n_

# A1: a 5 MeV alpha particle in 0.8 T (nonrelativistic, gamma - 1 ~ 1e-3)
m_a = physical_constants["alpha particle mass"][0]
K_a, B_a = 5e6 * QE, 0.8
v_a = np.sqrt(2 * K_a / m_a); r_a = m_a * v_a / (2 * QE * B_a); f_a = 2 * QE * B_a / (2 * np.pi * m_a)
K_p = K_a * m_a / (4 * m_p)
print(f"[9] A1: alpha, 5 MeV, 0.8 T: v = {v_a:.4e} m/s, r = {r_a:.4f} m, f = {f_a/1e6:.3f} MHz, T = {1e6/f_a:.4f} us, gamma-1 = {K_a/(m_a*C_LIGHT**2):.1e}; a proton on the same circle has K = {K_p/QE/1e6:.3f} MeV")

# A2: wires
Iw, sw = 30.0, 0.05
B_inf_w = MU0 * Iw / (2 * np.pi * sw)
B_fin_w = biot_savart([[sw, 0, 0]], np.array([[0, 0, -0.1]]), np.array([[0, 0, 0.1]]), np.array([Iw]))[0, 1]
sin_th = 0.1 / np.hypot(0.1, sw)
close(B_fin_w, MU0 * Iw / (4 * np.pi * sw) * 2 * sin_th, 1e-12, "finite wire (3.10)")
F_pw = MU0 * 10 * 25 / (2 * np.pi * 0.02)
print(f"[9] A2: infinite wire {B_inf_w*1e6:.1f} uT; 20 cm segment (sin theta = {sin_th:.4f}) {B_fin_w*1e6:.2f} uT = {B_fin_w/B_inf_w:.4f} x infinite; 10 A and 25 A antiparallel at 2 cm: {F_pw*1e3:.2f} mN/m")

# A3: a coil on axis
Nc, Ic, R3c, z3c = 50, 1.5, 0.04, 0.03
B_c0, B_cz = MU0 * Nc * Ic / (2 * R3c), B_loop_axis(Nc * Ic, R3c, z3c)
z_half = R3c * np.sqrt(2 ** (2 / 3) - 1)
close(B_loop_axis(1.0, R3c, z_half), 0.5 * B_loop_axis(1.0, R3c, 0.0), 1e-12, "half-field distance")
close(biot_savart([[0, 0, z3c]], *ring(R3c, 0, 360, Nc * Ic))[0, 2], B_cz, 2e-4, "coil on axis")
print(f"[9] A3: N = 50, R = 4 cm, I = 1.5 A: B(0) = {B_c0*1e3:.4f} mT, B(3 cm) = {B_cz*1e3:.4f} mT (factor {B_cz/B_c0:.3f}); B falls to half at z = {z_half*100:.3f} cm = {z_half/R3c:.4f} R")

# A4: a toroid
N4, I4 = 1200, 0.8
print("[9] A4: toroid N = 1200, I = 0.8 A: B(6, 8, 10 cm) = " + ", ".join(f"{MU0*N4*I4/(2*np.pi*s_)*1e3:.2f}" for s_ in (0.06, 0.08, 0.10)) + " mT")

# A5: right-hand rule, checked with the integrator; units
xw = (np.array([[-50.0, 0, 0]]), np.array([[50.0, 0, 0]]), np.array([1.0]))             # current along +x
B_up, B_side = biot_savart([[0, 0, 1.0]], *xw)[0], biot_savart([[0, 1.0, 0]], *xw)[0]
assert B_up[1] < 0 and abs(B_up[0]) + abs(B_up[2]) < 1e-12 and B_side[2] > 0                # -y above, +z beside
cw = biot_savart([[0, 0, 0]], *ring(0.1, 0, 96, -1.0))[0, 2]
assert cw < 0                                                                              # clockwise from +z: m and B along -z
F_rep, _, _ = small_loop_force_torque(ring(0.1, 0, 200, 2.0), [0, 0, 0.1], [0, 0, -1.0], 2e-3, 1.0)
assert F_rep[2] > 0                                                                        # opposite circulation: repelled
close(1 / np.sqrt(MU0 * EPS0), C_LIGHT, 1e-9, "c = 1/sqrt(mu0 eps0)")
print(f"[9] A5: wire along +x: B at +z is ({B_up[0]:.0e}, {B_up[1]:.2e}, {B_up[2]:.0e}) (-y), at +y ({B_side[0]:.0e}, {B_side[1]:.0e}, {B_side[2]:.2e}) (+z); clockwise loop: B_z(centre) = {cw:.2e} < 0; opposite loop above is pushed up: F_z = {F_rep[2]:.2e} N > 0; 1/sqrt(mu0 eps0) = {1/np.sqrt(MU0*EPS0):.6e} m/s")

# B2: J = J0 s/a inside a wire of radius a
mu_s, s_s, a_s, J0_s = sp.symbols("mu0 s a J0", positive=True)
B2_in = mu_s * J0_s * s_s ** 2 / (3 * a_s)
assert sp.simplify(sp.diff(s_s * B2_in, s_s) / s_s - mu_s * J0_s * s_s / a_s) == 0                     # (curl B)_z = mu0 J
a2, I2 = 0.01, 10.0
J02 = 3 * I2 / (2 * np.pi * a2 ** 2)
B2a = MU0 * I2 / (2 * np.pi * a2)
print(f"[9] B2: a = 1 cm, I = 10 A: J0 = {J02:.4e} A/m^2, B(a) = {B2a*1e6:.1f} uT, B(a/2) = {MU0*J02*(a2/2)**2/(3*a2)*1e6:.1f} uT (uniform J would give {B2a/2*1e6:.1f}), half the current inside s = {a2*2**(-1/3)*1e3:.3f} mm")
close(MU0 * J02 * (a2 / 2) ** 2 / (3 * a2), B2a / 4, 1e-12, "B(a/2) = B(a)/4")

# B3: anti-Helmholtz pair, d = sqrt(3) R, ring at +d/2 carries +NI, ring at -d/2 carries -NI
zz_, dd_ = zs_, ds
Bs_sym = (f_.subs(zz_, zz_ - dd_ / 2) - f_.subs(zz_, zz_ + dd_ / 2)) / 2                                  # B_z / (mu0 I R^2)
assert sp.solve(sp.simplify(sp.diff(Bs_sym, zz_, 3).subs(zz_, 0)), dd_) == [sp.sqrt(3) * Rs]
assert sp.simplify(sp.diff(Bs_sym, zz_).subs(zz_, 0).subs(dd_, sp.sqrt(3) * Rs) - 48 * sp.sqrt(3) / (sp.Integer(7) ** sp.Rational(5, 2) * Rs ** 4)) == 0
R3, NI3 = 0.05, 1000.0
d3 = np.sqrt(3) * R3
AH = join([ring(R3, d3 / 2, 200, NI3), ring(R3, -d3 / 2, 200, -NI3)])
Bz3 = lambda z_: biot_savart([[0, 0, z_]], *AH)[0, 2]
grad3 = (Bz3(1e-4) - Bz3(-1e-4)) / 2e-4
grad3_th = 48 * np.sqrt(3) / 7 ** 2.5 * MU0 * NI3 / R3 ** 2
close(grad3, grad3_th, 1e-3, "anti-Helmholtz gradient")
dev3 = [Bz3(f_ * R3) / (grad3 * f_ * R3) - 1 for f_ in (0.25, 0.5)]
Bs3 = biot_savart([[0.1 * R3, 0, 0]], *AH)[0, 0]
close(Bs3, -0.1 * R3 / 2 * grad3, 1e-2, "radial field -(s/2) B_z'")
AH2 = join([ring(R3, R3 / 2, 200, NI3), ring(R3, -R3 / 2, 200, -NI3)])
g_R = (biot_savart([[0, 0, 1e-4]], *AH2)[0, 2] - biot_savart([[0, 0, -1e-4]], *AH2)[0, 2]) / 2e-4
dev_R = biot_savart([[0, 0, 0.25 * R3]], *AH2)[0, 2] / (g_R * 0.25 * R3) - 1
print(f"[9] B3: R = 5 cm, NI = 1000 A each: B_z'(0) = {grad3:.4f} T/m (formula 48 sqrt3 / 7^(5/2) mu0 NI/R^2 = {grad3_th:.4f}); departure from linear at z = R/4, R/2: {dev3[0]*100:.2f}%, {dev3[1]*100:.2f}% (d = R: {dev_R*100:.1f}% at R/4); B_s(0.1 R) = {Bs3*1e6:.1f} uT vs -(s/2)B_z' = {-0.1*R3/2*grad3*1e6:.1f} uT")

# B4: two gauges for a uniform field B0 z_hat
x_, y_, z_ = sp.symbols("x y z"); B0s = sp.Symbol("B0")
A1s, A2s = sp.Matrix([-B0s * y_ / 2, B0s * x_ / 2, 0]), sp.Matrix([0, B0s * x_, 0])
curl_ = lambda A: sp.Matrix([sp.diff(A[2], y_) - sp.diff(A[1], z_), sp.diff(A[0], z_) - sp.diff(A[2], x_), sp.diff(A[1], x_) - sp.diff(A[0], y_)])
div_ = lambda A: sp.diff(A[0], x_) + sp.diff(A[1], y_) + sp.diff(A[2], z_)
assert curl_(A1s) == sp.Matrix([0, 0, B0s]) and curl_(A2s) == sp.Matrix([0, 0, B0s]) and div_(A1s) == 0 and div_(A2s) == 0
lam = B0s * x_ * y_ / 2
assert sp.simplify(A2s - A1s - sp.Matrix([sp.diff(lam, x_), sp.diff(lam, y_), sp.diff(lam, z_)])) == sp.zeros(3, 1)
assert sp.diff(lam, x_, 2) + sp.diff(lam, y_, 2) + sp.diff(lam, z_, 2) == 0
B04, a4 = 0.5, 0.1
sq = np.array([[0, 0, 0], [a4, 0, 0], [a4, a4, 0], [0, a4, 0], [0, 0, 0]], float)
def loop_integral(Afun, corners):
    dl = np.diff(corners, axis=0)
    return np.einsum("pk,pk->", Afun(corners[:-1] + dl / 2), dl)
I_A1 = loop_integral(lambda p: 0.5 * B04 * np.c_[-p[:, 1], p[:, 0], 0 * p[:, 0]], sq)
I_A2 = loop_integral(lambda p: B04 * np.c_[0 * p[:, 0], p[:, 0], 0 * p[:, 0]], sq)
close(I_A1, B04 * a4 ** 2, 1e-12, "oint A1"); close(I_A2, B04 * a4 ** 2, 1e-12, "oint A2")
print(f"[9] B4: B0 = 0.5 T, square of side 10 cm: oint A1.dl = {I_A1:.4e} Wb, oint A2.dl = {I_A2:.4e} Wb, flux = {B04*a4**2:.4e} Wb")

# B5: sheet K = 100 A/m
K5 = 100.0
print(f"[9] B5: K = 100 A/m: |B| = mu0 K / 2 = {MU0*K5/2*1e6:.2f} uT on either side; jump of dA_x/dz = -mu0 K = {-MU0*K5*1e6:.2f} uT")

# B6: magnetized cylinder, B and H (reuses CYL, Bc, Be, Rc, Lc, Mm from Part 5)
Hc_th = Bc_th / MU0 - Mm
close(Bc / MU0 - Mm, Hc_th, 1e-2, "H at the centre")
print(f"[9] B6: R = 5 mm, L = 20 mm, M = 1e6 A/m: centre B = {Bc_th:.4f} T (formula; ring stack {Bc:.4f}), H = {Hc_th:.4e} A/m = {Hc_th/Mm:.4f} M (long rod 0, sphere -M/3); end face B = {Be_th:.4f} T, H inside {Be_th/MU0-Mm:.4e}, outside {Be_th/MU0:.4e} A/m")

# B7: a permeable wire, mu_r = 200, a = 5 mm, I = 20 A
a7, I7, mur7 = 5e-3, 20.0, 200.0
chi7 = mur7 - 1
H7 = I7 / (2 * np.pi * a7)
Jf7 = I7 / (np.pi * a7 ** 2)
Jb7, Kb7 = chi7 * Jf7, -chi7 * H7                                   # K_b = M x s_hat = -M z_hat on the surface
assert abs(Jb7 * np.pi * a7 ** 2 + Kb7 * 2 * np.pi * a7) < 1e-9 * chi7 * I7
print(f"[9] B7: H(a) = {H7:.2f} A/m, B(a) = {MU0*mur7*H7:.4f} T, M(a) = {chi7*H7:.4e} A/m, J_b = {Jb7:.4e} A/m^2 (= {chi7:.0f} J_f), K_b = {Kb7:.4e} A/m; volume bound current {Jb7*np.pi*a7**2:.1f} A, surface {Kb7*2*np.pi*a7:.1f} A; B(2a) = {MU0*I7/(2*np.pi*2*a7)*1e3:.3f} mT")

# B8: force and torque on a small square loop (side 2 mm, 1 A) in the field of the 10 cm, 2 A ring, versus grad(m.B) and m x B
SRC8 = ring(0.10, 0, 400, 2.0)
c8, n8 = np.array([0.03, 0, 0.10]), np.array([np.sin(np.pi / 6), 0, np.cos(np.pi / 6)])
F8, N8, m8 = small_loop_force_torque(SRC8, c8, n8, 2e-3, 1.0)
B8c = biot_savart([c8], *SRC8)[0]
gradmB = np.array([(m8 @ biot_savart([c8 + e], *SRC8)[0] - m8 @ biot_savart([c8 - e], *SRC8)[0]) / 2e-4 for e in 1e-4 * np.eye(3)])
eF, eN = np.linalg.norm(F8 - gradmB) / np.linalg.norm(gradmB), np.linalg.norm(N8 - np.cross(m8, B8c)) / np.linalg.norm(np.cross(m8, B8c))
assert eF < 1e-4 and eN < 1e-3
print(f"[9] B8: m = {m8} A m^2, B = {B8c} T; F = {F8} N vs grad(m.B) = {gradmB} (rel. diff {eF:.1e}); N = {N8} N m vs m x B = {np.cross(m8, B8c)} (rel. diff {eN:.1e}); U = -m.B = {-m8@B8c:.4e} J")

errs8 = []
for w_ in (4e-3, 2e-3, 1e-3):                                           # errors fall as w^2: the midpoint rule on a small loop
    F_, N_, m_ = small_loop_force_torque(SRC8, c8, n8, w_, 1.0)
    g_ = np.array([(m_ @ biot_savart([c8 + e], *SRC8)[0] - m_ @ biot_savart([c8 - e], *SRC8)[0]) / 2e-4 for e in 1e-4 * np.eye(3)])
    errs8.append((np.linalg.norm(F_ - g_) / np.linalg.norm(g_), np.linalg.norm(N_ - np.cross(m_, B8c)) / np.linalg.norm(np.cross(m_, B8c))))
assert all(3.5 < errs8[i][j] / errs8[i + 1][j] < 4.5 for i in range(2) for j in range(2))
print("[9] B8: side 4, 2, 1 mm: force errors " + ", ".join(f"{e[0]:.1e}" for e in errs8) + "; torque errors " + ", ".join(f"{e[1]:.1e}" for e in errs8) + " (x4 per halving)")

# C1: the integral of B over a sphere is (2/3) mu0 m: the uniformly magnetized sphere has uniform B_in = (2/3) mu0 M
V_s = 4 / 3 * np.pi * Rm ** 3
close(Bin_s[0, 2] * V_s, 2 / 3 * MU0 * m_s, 2e-3, "volume integral of B over the sphere")
print(f"[9] C1: sphere of Part 5: B_in V = {Bin_s[0,2]*V_s:.4e} T m^3 vs (2/3) mu0 m = {2/3*MU0*m_s:.4e}")

# C2: a saturating core, B = Bs tanh(mu0 mu_r H / Bs), same C-core as Example 4.6
def gap_field(NI, Bs=1.6, mur=1000.0, l=0.636, g=4e-3):
    """Solve H(B) l + B g / mu0 = NI for B (same flux in core and gap), H(B) = (Bs / mu0 mur) artanh(B / Bs)."""
    return brentq(lambda B: Bs / (MU0 * mur) * np.arctanh(B / Bs) * l + B * g / MU0 - NI, 0, Bs * (1 - 1e-12))
def NI_needed(B, Bs=1.6, mur=1000.0, l=0.636, g=4e-3):
    return Bs / (MU0 * mur) * np.arctanh(B / Bs) * l + B * g / MU0
lin = lambda NI: MU0 * NI / (4e-3 + 0.636 / 1000)
close(gap_field(1000.0), lin(1000.0), 3e-3, "small NI is linear")
assert gap_field(1e4) < 1.6 and gap_field(1e4) > 1.599 and lin(1e4) > 2.7
close(NI_needed(gap_field(3000.0)), 3000.0, 1e-9, "solver round trip")
print("[9] C2: B_gap(NI) nonlinear vs linear: " + "; ".join(f"NI = {n_:.0f}: {gap_field(n_):.4f} vs {lin(n_):.4f} T" for n_ in (1000, 3000, 5000, 10000)))
print("[9] C2: NI needed for B = 0.5, 1.0, 1.5 T: " + ", ".join(f"{NI_needed(b_):.1f} (linear {b_*(4e-3+0.636/1000)/MU0:.1f})" for b_ in (0.5, 1.0, 1.5)))

# C3: ac demagnetization of the Preisach ensemble of Your turn 3
st = sgn_init.copy()
preisach_run(np.linspace(0, 3, 150), up, dn, st)                    # saturate at H = +3
M_sat = st.mean()
M_direct = preisach_run([0.0], up, dn, st.copy())[0]                 # remove the field at once
K_cyc = 80
Hdem = np.concatenate([3.0 * 0.9 ** np.arange(1, K_cyc + 1) * (-1.0) ** np.arange(1, K_cyc + 1), [0.0]])
M_dem = preisach_run(Hdem, up, dn, st)
assert abs(M_dem[-1]) < 0.05 < M_direct
st5 = sgn_init.copy(); preisach_run(np.linspace(0, 3, 150), up, dn, st5)
Hcoarse = np.concatenate([3.0 * 0.5 ** np.arange(1, 25) * (-1.0) ** np.arange(1, 25), [0.0]])
M_coarse = preisach_run(Hcoarse, up, dn, st5)
# the final state of a hysteron is set by whichever of its two thresholds is smaller in magnitude, i.e. by the sign of its bias hu = (up + dn)/2
M_pred = -np.mean(np.sign((up + dn) / 2))
st6 = sgn_init.copy(); preisach_run(np.linspace(0, 3, 150), up, dn, st6)
M_fine = preisach_run(np.concatenate([3.0 * 0.99 ** np.arange(1, 801) * (-1.0) ** np.arange(1, 801), [0.0]]), up, dn, st6)[-1]
assert abs(M_fine - M_pred) < 1e-3
H_up = np.linspace(0, 3, 150)
fig, ax = plt.subplots(figsize=(4.5, 4))
ax.plot(np.concatenate([H_up, Hdem]), np.concatenate([preisach_run(H_up, up, dn, sgn_init.copy()), M_dem]), "C0", lw=.8); ax.set_xlabel("H (arb.)"); ax.set_ylabel("M / M_s"); ax.set_title("ac demagnetization (factor 0.9)"); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "demag.png", dpi=100); plt.close(fig)
print(f"[9] C3: after saturating M = {M_sat:.3f}; field removed directly: M = {M_direct:.3f}; alternating field 3(-0.9)^k: M at H_k (k = 1, 10, 20, 40, 80) = " + ", ".join(f"{M_dem[k-1]:+.3f}" for k in (1, 10, 20, 40, 80)) + f"; final M(H = 0) = {M_dem[-1]:+.4f}; with factor 0.5 instead of 0.9 (24 steps): M = {M_coarse[-1]:+.4f}; factor 0.99 (800 steps): {M_fine:+.4f}; bias-sign prediction -mean(sign(hu)) = {M_pred:+.4f}")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time()-T0:.1f} s. Figures in", FIG)
