# %% [markdown]
# # EM-7 Project: Plane Waves, Polarization, Fresnel, Dispersion, Wave Packets, Waveguides
#
# Parts: (1) plane waves in vacuum checked symbolically, the time-average rule, Example 4.1;
# (2) Jones calculus: polarizers, wave plates, Malus, Example 4.2, the Poincare sphere;
# (3) Fresnel coefficients vs angle with the boundary conditions checked directly, Brewster and
# critical angles, the evanescent decay length (Examples 4.3-4.4); (4) 1D FDTD with a Lorentz
# auxiliary differential equation: r and t at a glass half-space; (5) a chirped pulse in a dispersive
# medium: FDTD vs FFT propagation vs the Gaussian spreading formula (3.30); (6) the Lorentz model:
# n and absorption for a two-resonance toy, Cauchy fit, plasma frequency; (7) conductors and skin
# depth (Example 4.5); (8) waveguides: TE fields checked against Maxwell with SymPy, dispersion,
# X-band numbers (Example 4.6), coax; Your turn (with reference solutions); (9) solutions.md numbers.
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
from scipy.constants import mu_0 as MU0, epsilon_0 as EPS0, c as C, e as QE, m_e as ME, hbar as HBAR

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()
Z0 = MU0 * C

def close(a, b, rtol, label=""):
    """assert a ~ b to relative tolerance rtol"""
    err = abs(a - b) / abs(b)
    assert err < rtol, f"{label}: {a} vs {b}, rel err {err:.2e}"
    return err

# %% [markdown]
# ## Part 1: Plane waves in vacuum (notes 3.2, 3.4-3.5, Example 4.1)
# E = E0 cos(kz - wt) x-hat, B = (E0/c) cos(kz - wt) y-hat must satisfy all four vacuum Maxwell
# equations when w = ck (3.6). The time average of a product of two phasor fields is (1/2)Re(A B*)
# (3.3). Two strings joined at z = 0 (3.4): check the amplitudes satisfy both junction conditions.

# %%
x, y, z, t = sp.symbols("x y z t", real=True)
E0s, ks, cs = sp.symbols("E_0 k c", positive=True)
ws = cs * ks
Ev = sp.Matrix([E0s * sp.cos(ks * z - ws * t), 0, 0])
Bv = sp.Matrix([0, E0s / cs * sp.cos(ks * z - ws * t), 0])
X = (x, y, z)
div = lambda F: sum(sp.diff(F[i], X[i]) for i in range(3))
curl = lambda F: sp.Matrix([sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x),
                            sp.diff(F[1], x) - sp.diff(F[0], y)])
assert div(Ev) == 0 and div(Bv) == 0 and sp.simplify(curl(Ev) + sp.diff(Bv, t)) == sp.zeros(3, 1)
assert sp.simplify(curl(Bv) - sp.diff(Ev, t) / cs**2) == sp.zeros(3, 1)
a1, a2, b1, b2 = sp.symbols("a1 a2 b1 b2", real=True); w_ = sp.symbols("omega", positive=True)
A_, B_ = a1 + sp.I * a2, b1 + sp.I * b2
prod = sp.re(A_ * sp.exp(-sp.I * w_ * t)) * sp.re(B_ * sp.exp(-sp.I * w_ * t))
avg = sp.integrate(sp.expand(sp.expand_complex(prod)), (t, 0, 2 * sp.pi / w_)) * w_ / (2 * sp.pi)
assert sp.simplify(avg - sp.re(A_ * sp.conjugate(B_)) / 2) == 0
k1s, k2s = sp.symbols("k1 k2", positive=True)                       # two strings, (3.4)
AR, AT = (k1s - k2s) / (k1s + k2s), 2 * k1s / (k1s + k2s)
assert sp.simplify(1 + AR - AT) == 0 and sp.simplify(k1s * (1 - AR) - k2s * AT) == 0

def plane_wave_numbers(I):
    """E0 = sqrt(2I/(c eps0)) from (3.7), B0 = E0/c, radiation pressure I/c absorbed, 2I/c reflected (3.8)"""
    E0 = np.sqrt(2 * I / (C * EPS0))
    return E0, E0 / C, I / C, 2 * I / C

E01, B01, Pa1, Pr1 = plane_wave_numbers(1000.0)
close(E01, 868.0, 1e-3, "E0"); close(B01, 2.895e-6, 1e-3, "B0"); close(Pr1, 6.671e-6, 1e-3, "P refl")
Es, Bs, Pas, Prs = plane_wave_numbers(1361.0)
print(f"[1] Maxwell's equations, the time-average rule and the string junction verified symbolically")
print(f"[1] Ex 4.1  I = 1 kW/m^2: E0 = {E01:.1f} V/m, B0 = {B01:.4e} T, <u> = {1000/C:.3e} J/m^3, "
      f"P_abs = {Pa1:.3e} Pa, P_mirror = {Pr1:.3e} Pa")
print(f"[1] sunlight I = 1361 W/m^2: E0 = {Es:.1f} V/m, P_abs = {Pas:.3e} Pa, P_mirror = {Prs:.3e} Pa; "
      f"force on a 100 m x 100 m mirror sail {Prs*1e4*1e3:.1f} mN")

# %% [markdown]
# ## Part 2: Jones calculus and the Poincare sphere (notes 3.6, Example 4.2)
# Jones vector J = (Ex, Ey); physical field Re[J e^{i(kz - wt)}]. Polarizer along u = (cos th, sin th):
# P = u u^T (a projection, LA-5). Wave plate, fast axis along x, retardance G: diag(1, e^{iG}) (the slow
# component falls behind, which with e^{-iwt} means an extra e^{+iG}). Rotated element: R M R^T.
# Stokes: S1 = |Ex|^2 - |Ey|^2, S2 = 2 Re(Ex* Ey), S3 = 2 Im(Ex* Ey), divided by S0 = |Ex|^2 + |Ey|^2.

# %%
rot = lambda th: np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])

def polarizer(th):
    """projection onto the transmission axis at angle th from x: P = u u^T (3.10)"""
    return np.outer([np.cos(th), np.sin(th)], [np.cos(th), np.sin(th)])

def waveplate(G, th=0.0):
    """retarder with fast axis at angle th, retardance G: R(th) diag(1, e^{iG}) R(-th) (3.10)"""
    return rot(th) @ np.diag([1, np.exp(1j * G)]) @ rot(-th)

def stokes(J):
    """normalized Stokes vector (S1, S2, S3)/S0, a point on the Poincare sphere (3.11)"""
    Ex, Ey = J; S0 = abs(Ex)**2 + abs(Ey)**2
    return np.array([abs(Ex)**2 - abs(Ey)**2, 2 * (np.conj(Ex) * Ey).real, 2 * (np.conj(Ex) * Ey).imag]) / S0

H, D45 = np.array([1, 0], complex), np.array([1, 1], complex) / np.sqrt(2)
Jc = waveplate(np.pi / 2) @ D45                                      # Example 4.2
assert np.allclose(Jc, np.array([1, 1j]) / np.sqrt(2)) and np.allclose(stokes(Jc), [0, 0, 1])
tt = np.linspace(0, 2 * np.pi, 400)
Ex_t, Ey_t = (Jc[0] * np.exp(-1j * tt)).real, (Jc[1] * np.exp(-1j * tt)).real   # field at z = 0
assert np.all(Ex_t * np.gradient(Ey_t) - Ey_t * np.gradient(Ex_t) > 0)          # turns from x toward y
ths = np.linspace(0, np.pi, 181)
I_after = [np.linalg.norm(polarizer(th) @ Jc)**2 for th in ths]
assert np.allclose(I_after, 0.5)                                     # circular: half through any polarizer
J2 = waveplate(np.pi / 2, np.pi / 4) @ H
assert np.allclose(J2, np.exp(1j * np.pi / 4) * np.array([1, -1j]) / np.sqrt(2))
assert np.allclose([np.linalg.norm(polarizer(th) @ H)**2 for th in ths], np.cos(ths)**2)   # Malus
I3 = 0.5 * np.linalg.norm(polarizer(np.pi / 2) @ polarizer(np.pi / 4) @ H)**2              # A3
assert abs(I3 - 1 / 8) < 1e-12
Jhwp = waveplate(np.pi, np.pi / 8) @ H                               # half-wave plate at 22.5 deg
assert np.allclose(stokes(Jhwp), stokes(D45))
assert np.allclose(waveplate(np.pi / 2), np.diag([1, 1j]))            # QWP = the phase gate S (QC-1)
assert np.allclose(waveplate(np.pi, np.pi / 8), np.array([[1, 1], [1, -1]]) / np.sqrt(2))   # HWP at 22.5 deg = Hadamard
traj = np.array([stokes(waveplate(np.pi / 2, th) @ H) for th in ths])          # rotating QWP
assert np.allclose(np.linalg.norm(traj, axis=1), 1)                  # unitary plate: stays on sphere
assert np.allclose(traj, np.c_[np.cos(2 * ths)**2, np.sin(2 * ths) * np.cos(2 * ths), -np.sin(2 * ths)])
print(f"[2] Ex 4.2  QWP(fast x) on 45-deg linear -> {np.round(Jc, 4)}, Stokes {stokes(Jc).round(3)}; "
      f"through any polarizer: {min(I_after):.3f}..{max(I_after):.3f}; three polarizers: I/I0 = {I3:.4f}")
fig = plt.figure(figsize=(5, 5)); ax = fig.add_subplot(projection="3d")
u_, v_ = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
ax.plot_wireframe(np.cos(u_) * np.sin(v_), np.sin(u_) * np.sin(v_), np.cos(v_), color="0.85", lw=0.5)
ax.plot(*traj.T, "C3", lw=2, label="H through a rotating QWP")
hw = np.array([stokes(waveplate(np.pi, th) @ H) for th in ths]); ax.plot(*hw.T, "C0", lw=1.5, label="rotating HWP")
for lab, p in {"H": (1, 0, 0), "V": (-1, 0, 0), "D": (0, 1, 0), "+circ": (0, 0, 1)}.items(): ax.text(*np.array(p) * 1.15, lab)
ax.set_xlabel("S1"); ax.set_ylabel("S2"); ax.set_zlabel("S3"); ax.legend(fontsize=7)
ax.set_title("Poincare sphere (= Bloch sphere, QC-1)"); fig.savefig(FIG / "poincare.png", dpi=100); plt.close(fig)

# %% [markdown]
# ## Part 3: Fresnel equations, Brewster, total internal reflection (notes 3.8-3.11, Examples 4.3-4.4)
# s: r = (n1 ci - n2 ct)/(n1 ci + n2 ct), t = 2 n1 ci/(n1 ci + n2 ct).  p (Griffiths' sign, all three
# p-vectors equal at normal incidence): r = (n1 ct - n2 ci)/(n1 ct + n2 ci), t = 2 n1 ci/(n1 ct + n2 ci).
# Beyond the critical angle ct = +i sqrt(sin^2 th_t - 1): the + sign makes e^{i k_z z} decay under e^{-iwt}.

# %%
def fresnel(n1, n2, th):
    """(r_s, t_s, r_p, t_p, cos th_i, cos th_t) from (3.16)-(3.17); complex beyond the critical angle"""
    ci = np.cos(th) + 0j
    ct = np.sqrt(1 - (n1 * np.sin(th) / n2)**2 + 0j)                  # principal root: Im(ct) >= 0
    rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct); ts = 2 * n1 * ci / (n1 * ci + n2 * ct)
    rp = (n1 * ct - n2 * ci) / (n1 * ct + n2 * ci); tp = 2 * n1 * ci / (n1 * ct + n2 * ci)
    return rs, ts, rp, tp, ci, ct

def RT(n1, n2, th):
    """reflectances and transmittances (R_s, T_s, R_p, T_p); T = Re(n2 ct)/(n1 ci) |t|^2"""
    rs, ts, rp, tp, ci, ct = fresnel(n1, n2, th); f = (n2 * ct).real / (n1 * ci.real)
    return abs(rs)**2, f * abs(ts)**2, abs(rp)**2, f * abs(tp)**2

rs0, ts0, rp0, tp0, _, _ = fresnel(1.0, 1.5, 0.0)                     # Example 4.3
Rs0, Ts0, _, _ = RT(1.0, 1.5, 0.0)
assert np.isclose(rs0, -0.2) and np.isclose(ts0, 0.8) and np.isclose(rp0, -0.2) and np.isclose(Rs0, 0.04) and np.isclose(Ts0, 0.96)
win2, win_inc = (1 - Rs0)**2, (1 - Rs0)**2 / (1 - Rs0**2)
th45 = np.pi / 4                                                      # Example 4.4
rs, ts, rp, tp, ci, ct = fresnel(1.0, 1.5, th45); Rs, Ts, Rp, Tp = RT(1.0, 1.5, th45)
for val, ref, lab in [(rs, -0.30334, "rs"), (ts, 0.69666, "ts"), (rp, -0.09201, "rp"), (tp, 0.72801, "tp"),
                      (Rs, 0.09202, "Rs"), (Rp, 0.008466, "Rp")]:
    close(val.real, ref, 2e-4, lab)
thB = np.arctan(1.5); thc = np.arcsin(1 / 1.5)
assert abs(fresnel(1.0, 1.5, thB)[2]) < 1e-14 and abs(fresnel(1.5, 1.0, np.arctan(1 / 1.5))[2]) < 1e-14
for n1, n2 in [(1.0, 1.5), (1.5, 1.0)]:                               # energy: R + T = 1, every angle
    for th in np.linspace(0, np.pi / 2 - 1e-6, 500):
        R_s, T_s, R_p, T_p = RT(n1, n2, th)
        assert abs(R_s + T_s - 1) < 1e-12 and abs(R_p + T_p - 1) < 1e-12

def bc_residuals(n1, n2, th):
    """Build the three waves at z = 0 and return the jumps in E_par, eps E_perp, H_par, B_perp (all should be 0)."""
    rs, ts, rp, tp, ci, ct = fresnel(n1, n2, th); si, st = np.sin(th), n1 * np.sin(th) / n2
    kI, kR, kT = np.array([si, 0, ci]), np.array([si, 0, -ci]), np.array([st, 0, ct])
    res = []
    for EI, ER, ET in [(np.array([0, 1, 0]), rs * np.array([0, 1, 0]), ts * np.array([0, 1, 0])),        # s
                       (np.array([ci, 0, -si]), rp * np.array([ci, 0, si]), tp * np.array([ct, 0, -st]))]:  # p
        B1, B2, E1, E2 = n1 * (np.cross(kI, EI) + np.cross(kR, ER)), n2 * np.cross(kT, ET), EI + ER, ET   # cB = n k x E
        res += [abs(E1[:2] - E2[:2]).max(), abs(n1**2 * E1[2] - n2**2 * E2[2]), abs(B1 - B2).max()]
    return max(res)

assert max(bc_residuals(1.0, 1.5, th) for th in np.linspace(0, 1.5, 40)) < 1e-12
assert bc_residuals(1.5, 1.0, np.radians(60)) < 1e-12                 # also beyond the critical angle
assert np.isclose(abs(fresnel(1.5, 1.0, np.radians(60))[0]), 1) and np.isclose(abs(fresnel(1.5, 1.0, np.radians(60))[2]), 1)

def decay_length(n1, n2, th, lam0=1.0):
    """evanescent 1/e field decay length d = lam0 / (2 pi sqrt(n1^2 sin^2 th - n2^2))  (3.19)"""
    return lam0 / (2 * np.pi * np.sqrt((n1 * np.sin(th))**2 - n2**2))

d45 = decay_length(1.5, 1.0, th45)
assert np.isclose(d45, 1 / (2 * np.pi * 1.5 * fresnel(1.5, 1.0, th45)[5].imag / 1.5 * 1.0))    # = 1/Im(k_z)
print(f"[3] Ex 4.3  n = 1.5, normal: r = {rs0.real:.3f}, t = {ts0.real:.3f}, R = {Rs0:.3f}, T = {Ts0:.3f}; "
      f"window (two faces) T = {win2:.4f}, with incoherent multiple reflections {win_inc:.4f}; 12 faces: {Ts0**12:.3f}")
print(f"[3] Ex 4.4  45 deg air->glass: theta_t = {np.degrees(np.arccos(ct.real)):.3f} deg, cos = {ct.real:.5f}; "
      f"r_s = {rs.real:.5f}, t_s = {ts.real:.5f}, r_p = {rp.real:.5f}, t_p = {tp.real:.5f}")
print(f"[3]          R_s = {Rs:.5f}, T_s = {Ts:.5f}, R_p = {Rp:.6f}, T_p = {Tp:.5f}, unpolarized R = {(Rs+Rp)/2:.5f}, "
      f"polarization of reflected light (R_s-R_p)/(R_s+R_p) = {(Rs-Rp)/(Rs+Rp):.4f}")
print(f"[3] Brewster {np.degrees(thB):.2f} deg (glass->air {np.degrees(np.arctan(1/1.5)):.2f}), critical {np.degrees(thc):.2f} deg; "
      f"decay length glass->air at 45 deg = {d45:.4f} lambda0 = {d45*633:.0f} nm at 633 nm; boundary conditions hold to 1e-12")
angs = np.linspace(0, np.pi / 2 - 1e-4, 600)
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
for ax, (n1, n2) in zip(axs[:2], [(1.0, 1.5), (1.5, 1.0)]):
    r_s, _, r_p, _, _, _ = fresnel(n1, n2, angs); deg = np.degrees(angs)
    ax.plot(deg, r_s.real, "C0", label="r_s (real part)"); ax.plot(deg, r_p.real, "C1", label="r_p (real part)")
    ax.plot(deg, abs(r_s)**2, "C0--", label="R_s"); ax.plot(deg, abs(r_p)**2, "C1--", label="R_p")
    ax.axvline(np.degrees(np.arctan(n2 / n1)), color="k", ls=":", label="Brewster")
    if n1 > n2: ax.axvline(np.degrees(np.arcsin(n2 / n1)), color="C3", ls=":", label="critical")
    ax.set_title(f"n1 = {n1}, n2 = {n2}"); ax.set_xlabel("angle of incidence (deg)"); ax.legend(fontsize=7); ax.grid(alpha=.3)
ab = np.linspace(thc + 1e-3, np.pi / 2, 300)
axs[2].plot(np.degrees(ab), decay_length(1.5, 1.0, ab)); axs[2].set_ylim(0, 2); axs[2].grid(alpha=.3)
axs[2].set_title("evanescent decay length / lambda0 (glass->air)"); axs[2].set_xlabel("angle (deg)")
fig.tight_layout(); fig.savefig(FIG / "fresnel.png", dpi=100); plt.close(fig)

# %% [markdown]
# ## Part 4: 1D FDTD with a Lorentz/Drude medium (notes 3.8, 3.13; extends EM-6's yee1d)
# Grid units: dx = c = eps0 = mu0 = 1, dt = S. E_y on integer nodes, B_z on half nodes (EM-6 section 7).
#   B^{n+1/2}_{i+1/2} = B^{n-1/2}_{i+1/2} - S (E^n_{i+1} - E^n_i)
#   eps_inf (E^{n+1} - E^n) = -S (B_{i+1/2} - B_{i-1/2}) - (P^{n+1} - P^n)        [from dD/dt = -dB/dx]
#   P'' + gam P' + w0^2 P = wp2 E  (the driven-electron equation (3.23); Drude/plasma when w0 = 0),
#   centred differences: P^{n+1} = [(2 - w0^2 dt^2) P^n - (1 - gam dt/2) P^{n-1} + wp2 dt^2 E^n]/(1 + gam dt/2).
# Static limit P = (wp2/w0^2) E, so eps(0) = eps_inf + wp2/w0^2. First-order Mur absorbing ends (in vacuum).

# %%
def fdtd(N, steps, f_launch, x_ref, eps_inf=1.0, wp2=0.0, w0=0.0, gam=0.0, S=0.5, probes=(), snaps=()):
    """Right-moving pulse with time signal f(t) at x_ref: E(x, t) = f(t + x_ref - x) in vacuum at t = 0.
    Returns E snapshots {step: E} and probe records (steps, len(probes))."""
    dt = S; xs = np.arange(N, dtype=float); xh = xs[:-1] + 0.5
    E = f_launch(x_ref - xs); B = f_launch(-dt / 2 + x_ref - xh)
    ei, wp2, w0, g = [np.broadcast_to(np.asarray(v, float), (N,)) for v in (eps_inf, wp2, w0, gam)]
    a1 = (2 - w0**2 * dt**2) / (1 + g * dt / 2); a2 = (1 - g * dt / 2) / (1 + g * dt / 2); a3 = wp2 * dt**2 / (1 + g * dt / 2)
    P, Pm = np.zeros(N), np.zeros(N); kmur = (S - 1) / (S + 1)
    rec, out = np.zeros((steps, len(probes))), {}             # probe records cost nothing when probes = ()
    for n in range(steps):
        B -= S * (E[1:] - E[:-1])
        Pn = a1 * P - a2 * Pm + a3 * E; Eold = E.copy()
        E[1:-1] += (-S * (B[1:] - B[:-1]) - (Pn[1:-1] - P[1:-1])) / ei[1:-1]
        E[0] = Eold[1] + kmur * (E[1] - Eold[0]); E[-1] = Eold[-2] + kmur * (E[-2] - Eold[-1])
        Pm, P = P, Pn; rec[n] = E[list(probes)]
        if n + 1 in snaps: out[n + 1] = E.copy()
    return out, rec

N4, i4 = 2000, 1000
eps4 = np.where(np.arange(N4) > i4, 2.25, 1.0); eps4[i4] = (1 + 2.25) / 2
gauss = lambda s: np.exp(-((s - 500.0) / 40.0)**2)                    # pulse centred 500 cells left of x_ref
E4 = fdtd(N4, 1800, gauss, x_ref=1000.0, eps_inf=eps4, snaps=(1800,))[0][1800]
r_num, t_num = E4[:i4].min(), E4[i4:].max()
R_num = np.sum(E4[:i4]**2) / np.sum(eps4 * E4**2); T_num = np.sum(eps4[i4:] * E4[i4:]**2) / np.sum(eps4 * E4**2)
close(r_num, -0.2, 5e-3, "FDTD r"); close(t_num, 0.8, 5e-3, "FDTD t"); close(R_num, 0.04, 1e-2, "FDTD R"); close(T_num, 0.96, 1e-3, "FDTD T")
print(f"[4] FDTD glass half-space (n = 1.5): r = {r_num:.4f} (theory -0.2), t = {t_num:.4f} (0.8), "
      f"R = {R_num:.4f} (0.04), T = {T_num:.4f} (0.96)")

# %% [markdown]
# ## Part 5: A chirped pulse in a dispersive (Lorentz) medium: FDTD vs FFT propagation vs (3.30)
# Medium: eps(w) = 1 + Deps w0^2/(w0^2 - w^2), Deps = 1.25, w0 = 0.3 (grid units), lossless. A Gaussian whose
# spectrum is exp(-tau0^2 q^2) e^{+i b0 q^2/2} (q = w - wc) is a transform-limited pulse that has already
# "propagated" a distance -b0/k''. After depth z its intensity width is tau0 sqrt(1 + ((b0 + k'' z)/(2 tau0^2))^2) (3.30).
# We launch b0 < 0 (opposite to the medium's k'' > 0): the pulse first compresses, then GVD spreads it.
# FFT propagation multiplies each frequency by e^{i k(w) z} (3.28). NumPy's FFT synthesizes with e^{+iwt}, the
# engineering sign, so in code the factor is exp(-1j * conj(k(w)) * z): the sec. 8 sign trap, in one line.

# %%
DEPS, W0L, WC, TAU0 = 1.25, 0.3, 0.12, 60.0
eps_lor = lambda w: 1 + DEPS * W0L**2 / (W0L**2 - w**2 + 0j)
k_lor = lambda w: w * np.sqrt(eps_lor(w))
hh = 1e-4; k1L = ((k_lor(WC + hh) - k_lor(WC - hh)) / (2 * hh)).real; k2L = ((k_lor(WC + hh) - 2 * k_lor(WC) + k_lor(WC - hh)) / hh**2).real
Z1, Z2 = 1000, 3000                                                   # probe depths in the medium
B0C = -k2L * Z1                                                      # chirp that compresses at depth Z1
a0 = TAU0**2 - 1j * B0C / 2                                          # envelope exp(-t^2/(4 a)) with a = tau0^2 - i b/2
launch = lambda s: (np.exp(-s**2 / (4 * a0)) * np.exp(-1j * WC * s)).real   # f(t), centred at t = 0
width = lambda b: TAU0 * np.sqrt(1 + (b / (2 * TAU0**2))**2)

def rms_width(sig, dt):
    """standard deviation in time of |analytic signal|^2 (the intensity envelope)"""
    F = np.fft.fft(sig); F[len(sig) // 2:] = 0; env2 = abs(np.fft.ifft(2 * F))**2; tt = np.arange(len(sig)) * dt
    m = np.sum(tt * env2) / np.sum(env2); return np.sqrt(np.sum((tt - m)**2 * env2) / np.sum(env2))

def apply_filter(sig, dt, Hfun):
    """Multiply each physical frequency component e^{-iwt} of a real signal by H(w). NumPy's rfft bins are
    e^{+iwt} components, i.e. physical frequency -w, where the response is H(-w) = conj(H(w))."""
    w = 2 * np.pi * np.fft.rfftfreq(len(sig), dt)
    Hw = np.where(w > 0, Hfun(np.maximum(w, 1e-9)), 1.0)
    return np.fft.irfft(np.fft.rfft(sig) * np.conj(Hw), len(sig))

def propagate(sig, dt, kfun, zz):
    """E(z, t) from E(0, t): multiply each component by e^{i k(w) z} (3.28)"""
    return apply_filter(sig, dt, lambda w: np.exp(1j * kfun(w) * zz))

N5, x5, S5 = 4600, 1000, 0.9                                         # interface at node 1000, pulse starts 600 cells left
steps5 = int(6800 / S5)                                              # S = 0.9: less grid dispersion
lor = np.arange(N5) >= x5
_, rec5 = fdtd(N5, steps5, lambda s: launch(s - 600.0), x_ref=float(x5), wp2=np.where(lor, DEPS * W0L**2, 0.0),
               w0=np.where(lor, W0L, 0.0), S=S5, probes=(x5 + 5, x5 + Z1, x5 + Z2))
dt5 = S5; ts5 = np.arange(steps5) * dt5
w_fd = [rms_width(rec5[:, j], dt5) for j in range(3)]
# continuous-model reference: the transmitted signal just inside the interface, propagated by FFT
inc = launch(ts5 - 600.0)                                            # incident signal at the interface
sig_in = apply_filter(inc, dt5, lambda w: 2 / (1 + np.sqrt(eps_lor(w))))  # times t(w) from (3.13)
sig_in = propagate(sig_in, dt5, k_lor, 5.0)
w_fft = [rms_width(propagate(sig_in, dt5, k_lor, zz - 5.0), dt5) for zz in (5.0, Z1, Z2)]
w_th = [width(B0C + k2L * zz) for zz in (5.0, Z1, Z2)]
for j in range(3):
    close(w_fd[j], w_fft[j], 0.03, f"FDTD vs FFT width at probe {j}"); close(w_fft[j], w_th[j], 0.03, f"FFT vs (3.30) at {j}")
assert w_fd[1] < w_fd[0] and w_fd[2] > w_fd[1]                       # compress, then spread
tpk = [ts5[np.argmax(abs(np.fft.ifft(np.where(np.arange(steps5) < steps5 // 2, 2, 0) * np.fft.fft(rec5[:, j]))))] for j in range(3)]
vg_fd = (Z2 - Z1) / (tpk[2] - tpk[1])
close(vg_fd, 1 / k1L, 0.01, "group velocity")
print(f"[5] Lorentz medium at w = {WC}: n = {np.sqrt(eps_lor(WC)).real:.4f}, v_p = {WC/k_lor(WC).real:.4f} c, "
      f"v_g = 1/k' = {1/k1L:.4f} c (FDTD peak arrival: {vg_fd:.4f} c), k'' = {k2L:.3f}")
print(f"[5] chirped pulse widths at depths 5, {Z1}, {Z2}: FDTD {np.round(w_fd, 1)}, FFT {np.round(w_fft, 1)}, "
      f"Gaussian formula {np.round(w_th, 1)} (tau0 = {TAU0})")
fig, axs = plt.subplots(3, 1, figsize=(8, 6))
for j, (ax, zz) in enumerate(zip(axs, (5, Z1, Z2))):
    ax.plot(ts5, rec5[:, j], lw=0.8, label=f"FDTD, depth {zz}")
    ax.plot(ts5, propagate(sig_in, dt5, k_lor, zz - 5.0), "k:", lw=0.8, label="FFT propagation")
    c0 = tpk[j]; ax.set_xlim(c0 - 600, c0 + 600); ax.legend(fontsize=7, loc="upper right")
axs[-1].set_xlabel("time (grid units)"); axs[0].set_title("chirped pulse: compressed by GVD at depth 1000, then spread")
fig.tight_layout(); fig.savefig(FIG / "dispersive_fdtd.png", dpi=100); plt.close(fig)
# Pure quadratic k(w): the Gaussian formula is exact, so FFT propagation must match it to round-off.
kq = lambda w: 0.5 * w + 0.5 * k2L * (w - WC)**2
ref = launch(np.arange(16384) * 1.0 - 4000.0)
for zz in (300.0, 1000.0):
    close(rms_width(propagate(ref, 1.0, kq, zz), 1.0), width(B0C + k2L * zz), 1e-4, "quadratic k: exact spreading")

# %% [markdown]
# ## Part 6: The Lorentz model, Sellmeier/Cauchy, Rayleigh, the plasma frequency (notes 3.13)
# eps_r(w) = 1 + sum_j F_j/(w_j^2 - w^2 - i g_j w) (3.24). With g_j = 0 and lambda = 2 pi c/w this is the Sellmeier
# form n^2 = 1 + sum_j B_j lambda^2/(lambda^2 - lambda_j^2), B_j = F_j/w_j^2. Toy "glass": a UV resonance at
# 0.09 um (B = 1.1) and an IR one at 9.9 um (B = 0.9). Below the UV resonance n ~ A + B/lambda^2 (Cauchy, (3.26)).

# %%
BJ, LJ = np.array([1.1, 0.9]), np.array([0.09, 9.9])                 # um
GJ = np.array([0.03, 0.05])                                          # damping / w_j, for the plot only
def n_toy(lam, damp=False):
    """complex index from (3.24) written with lambda in um; damping g_j = GJ w_j if damp"""
    wr = LJ[:, None] / lam[None, :]                                  # w / w_j
    g = GJ[:, None] if damp else 0.0
    return np.sqrt(1 + np.sum(BJ[:, None] / (1 - wr**2 - 1j * g * wr), axis=0) + 0j)

lam_v = np.linspace(0.4, 0.7, 61); n_v = n_toy(lam_v).real
Ac, Bc = np.polyfit(1 / lam_v**2, n_v, 1)[::-1]
n0 = np.sqrt(1 + BJ[0]); B_th = BJ[0] * LJ[0]**2 / (2 * n0)
assert np.abs(Ac + Bc / lam_v**2 - n_v).max() < 2e-4
close(Bc, B_th, 0.15, "Cauchy B")
nD = n_toy(np.array([0.5893]))[0].real; assert np.all(np.diff(n_v) < 0)   # normal dispersion: n falls with lambda
lam_all = np.logspace(-1.6, 1.6, 3000); nn = n_toy(lam_all, damp=True)
alpha = 2 * (2 * np.pi / lam_all) * nn.imag                          # absorption coefficient 2 kappa (1/um)
ray, wp_ion, wp_cu = (700 / 450)**4, np.sqrt(1e12 * QE**2 / (ME * EPS0)), np.sqrt(8.47e28 * QE**2 / (ME * EPS0))
close(wp_ion / (2 * np.pi), 8.98e6, 1e-3, "ionosphere f_p")
print(f"[6] toy glass: n(589.3 nm) = {nD:.4f}; Cauchy fit 0.4-0.7 um: A = {Ac:.5f}, B = {Bc:.5f} um^2 "
      f"(one-resonance theory n0 = {n0:.5f}, B = {B_th:.5f}); max fit error {np.abs(Ac + Bc/lam_v**2 - n_v).max():.1e}")
print(f"[6] Rayleigh (700/450)^4 = {ray:.2f}; ionosphere N = 1e12 m^-3: w_p = {wp_ion:.4e} rad/s, f_p = {wp_ion/2/np.pi/1e6:.2f} MHz; "
      f"copper's free electrons (N = 8.47e28 m^-3): w_p = {wp_cu:.3e} rad/s, lambda_p = 2 pi c/w_p = {2*np.pi*C/wp_cu*1e9:.0f} nm")
fig, axs = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
axs[0].semilogx(lam_all, nn.real); axs[0].axvspan(0.4, 0.7, color="0.9"); axs[0].set_ylabel("n"); axs[0].set_ylim(0, 3.5)
axs[0].set_title("two-resonance Lorentz toy: anomalous dispersion inside each absorption band")
axs[1].loglog(lam_all, alpha); axs[1].set_ylabel("alpha (1/um)"); axs[1].set_xlabel("wavelength (um)"); axs[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "lorentz_model.png", dpi=100); plt.close(fig)

# %% [markdown]
# ## Part 7: Conductors and skin depth (notes 3.12, Example 4.5)
# k~^2 = mu eps w^2 + i mu sigma w (3.20). Exact kappa = Im k~; good-conductor limit delta = sqrt(2/(mu sigma w)) (3.21).

# %%
SIG_CU = 5.96e7
def k_cond(f, sig, eps=EPS0, mu=MU0):
    """complex wave number in an ohmic conductor, root with Im k > 0 (decays along +z)"""
    w = 2 * np.pi * f; return np.sqrt(mu * eps * w**2 + 1j * mu * sig * w)

skin = {}
for f in (60.0, 1e9):
    kc = k_cond(f, SIG_CU); d_ex, d_ap = 1 / kc.imag, np.sqrt(2 / (MU0 * SIG_CU * 2 * np.pi * f))
    close(d_ex, d_ap, 1e-6, "skin depth"); skin[f] = (d_ap, np.degrees(np.angle(kc)), C * abs(kc) / (2 * np.pi * f))
close(skin[60.0][0], 8.42e-3, 2e-3, "Cu 60 Hz"); close(skin[1e9][0], 2.062e-6, 2e-3, "Cu 1 GHz")
bt = C * k_cond(1e9, SIG_CU) / (2 * np.pi * 1e9)                     # beta~ = c k~/w at a vacuum-copper face
R_cu = abs((1 - bt) / (1 + bt))**2
print(f"[7] sigma/(eps0 w) at 1 GHz = {SIG_CU/(EPS0*2*np.pi*1e9):.3e}")
print(f"[7] Ex 4.5 copper: delta(60 Hz) = {skin[60.0][0]*1e3:.2f} mm, delta(1 GHz) = {skin[1e9][0]*1e6:.3f} um; "
      f"phase of k~ = {skin[1e9][1]:.4f} deg; cB/E = {skin[1e9][2]:.3e} at 1 GHz; R = {R_cu:.6f} (1 - R = {1-R_cu:.2e})")

# %% [markdown]
# ## Part 8: Rectangular waveguide (notes 3.15, Example 4.6)
# TE_mn: B_z = B0 cos(m pi x/a) cos(n pi y/b) e^{i(kz - wt)}, transverse fields (3.32), k^2 = (w^2 - w_mn^2)/c^2,
# w_mn = c pi sqrt(m^2/a^2 + n^2/b^2). SymPy checks all of Maxwell's equations and E_par = 0 on the walls.

# %%
a_, b_, wv, B0_ = sp.symbols("a b omega B_0", positive=True)
m_, n_ = sp.symbols("m", integer=True, positive=True), sp.symbols("n", integer=True, nonnegative=True)
g2 = sp.pi**2 * (m_**2 / a_**2 + n_**2 / b_**2); kz = sp.sqrt(wv**2 / cs**2 - g2)
Bz = B0_ * sp.cos(m_ * sp.pi * x / a_) * sp.cos(n_ * sp.pi * y / b_) * sp.exp(sp.I * (kz * z - wv * t))
Eg = sp.Matrix([sp.I * wv / g2 * sp.diff(Bz, y), -sp.I * wv / g2 * sp.diff(Bz, x), 0])
Bg = sp.Matrix([sp.I * kz / g2 * sp.diff(Bz, x), sp.I * kz / g2 * sp.diff(Bz, y), Bz])
for expr in [div(Eg), div(Bg)] + list(curl(Eg) + sp.diff(Bg, t)) + list(curl(Bg) - sp.diff(Eg, t) / cs**2):
    assert sp.simplify(expr) == 0
assert all(sp.simplify(e) == 0 for e in (Eg[1].subs(x, 0), Eg[1].subs(x, a_), Eg[0].subs(y, 0), Eg[0].subs(y, b_)))
A_WG, B_WG = 0.023, 0.010
fc = lambda m, n: C / 2 * np.sqrt((m / A_WG)**2 + (n / B_WG)**2)
modes = sorted([(fc(m, n), m, n) for m in range(4) for n in range(3) if (m, n) != (0, 0)])[:7]
f10, f = fc(1, 0), 10e9
root = np.sqrt(1 - (f10 / f)**2); vp, vg = C / root, C * root
close(f10, 6.517e9, 1e-3, "TE10 cutoff"); assert np.isclose(vp * vg, C**2) and vp > C > vg
close(fc(1, 1), 16.35e9, 1e-3, "TE11 cutoff"); close(fc(2, 0), 13.03e9, 1e-3, "TE20 cutoff")
print("[8] TE fields satisfy Maxwell's equations and E_par = 0 on all four walls (SymPy)")
print(f"[8] Ex 4.6  a = 2.3 cm, b = 1.0 cm: f_c(TE10) = {f10/1e9:.3f} GHz; at 10 GHz sqrt(1 - (fc/f)^2) = {root:.5f}, "
      f"v_p = {vp:.4e} m/s = {vp/C:.4f} c, v_g = {vg:.4e} m/s = {vg/C:.4f} c, lambda_g = {C/f/root*100:.3f} cm, "
      f"k = {2*np.pi*f/C*root:.2f} rad/m")
print("[8] lowest cutoffs (GHz): " + ", ".join(f"TE{m}{n} {fv/1e9:.2f}" for fv, m, n in modes))
Zcoax = Z0 / (2 * np.pi * 1.5) * np.log(3.5)
print(f"[8] coax, polyethylene (eps_r = 2.25), b/a = 3.5: Z = {Zcoax:.2f} ohm; v = c/1.5 = {C/1.5:.4e} m/s")
fig, axs = plt.subplots(1, 5, figsize=(19, 3.2))
xg, yg = np.meshgrid(np.linspace(0, A_WG, 41), np.linspace(0, B_WG, 19))
for ax, (m, n) in zip(axs[:3], [(1, 0), (2, 0), (1, 1)]):
    Ex_ = -n * np.pi / B_WG * np.cos(m * np.pi * xg / A_WG) * np.sin(n * np.pi * yg / B_WG)   # E ~ (dBz/dy, -dBz/dx)
    Ey_ = m * np.pi / A_WG * np.sin(m * np.pi * xg / A_WG) * np.cos(n * np.pi * yg / B_WG)
    ax.pcolormesh(xg * 100, yg * 100, np.hypot(Ex_, Ey_), shading="auto", cmap="viridis")
    ax.quiver(xg[::3, ::3] * 100, yg[::3, ::3] * 100, Ex_[::3, ::3], Ey_[::3, ::3], color="w")
    ax.set_title(f"TE{m}{n}: transverse E"); ax.set_aspect("equal"); ax.set_xlabel("x (cm)")
kk_ = np.linspace(0, 600, 300)
for fv, m, n in modes[:5]:
    axs[3].plot(kk_, np.sqrt((2 * np.pi * fv)**2 + (C * kk_)**2) / 2 / np.pi / 1e9, label=f"TE{m}{n}")
axs[3].plot(kk_, C * kk_ / 2 / np.pi / 1e9, "k--", label="light line w = ck"); axs[3].axhline(10, color="0.5", ls=":")
axs[3].set_xlabel("k (rad/m)"); axs[3].set_ylabel("f (GHz)"); axs[3].legend(fontsize=6); axs[3].set_title("dispersion")
ff = np.linspace(1.001 * f10, 3 * f10, 300); rt = np.sqrt(1 - (f10 / ff)**2)
axs[4].plot(ff / 1e9, 1 / rt, label="v_p / c"); axs[4].plot(ff / 1e9, rt, label="v_g / c"); axs[4].axhline(1, color="k", lw=0.5)
axs[4].set_ylim(0, 4); axs[4].set_xlabel("f (GHz)"); axs[4].legend(fontsize=7); axs[4].set_title("TE10: v_p v_g = c^2")
fig.tight_layout(); fig.savefig(FIG / "waveguide.png", dpi=100); plt.close(fig)

# %% [markdown]
# ## Your turn
# 1. **Anti-reflection coating (exercise B8).** Derive that a layer of index n_c = sqrt(n_s) and thickness
#    lambda0/(4 n_c) kills normal-incidence reflection at lambda0. Then verify in `fdtd`: substrate n_s = 2.25,
#    coating n_c = 1.5 of 20 cells, carrier wavelength 120 cells, long envelope. Compare reflected energy bare vs coated.
# 2. **Evanescent decay vs angle.** Part 3 plots d(theta). Add frustrated TIR: a gap of thickness g between two glass
#    blocks; compute the transmitted fraction vs g/lambda0 at 45 deg with a 3-layer transfer matrix (s polarization).
# 3. **Plasma mirror (the FDTD side of exercise B5).** Make the right half a Drude plasma (w0 = 0, wp2 = w_p^2) and send
#    pulses with carrier below and above w_p. Measure the reflected energy fraction and compare with |r|^2 from (3.27).
# 4. Flip the sign of B0C in Part 5. Predict the widths at the two probes from (3.30) before you run it.
# Reference solutions for 1 and 3 follow. Write your own first.

# %%
def reflected_fraction(eps_profile, carrier_w, tau, N, x_int, extra=None):
    """energy reflected into x < x_int - 100 from a long pulse, after it has fully interacted"""
    pulse = lambda s: np.exp(-s**2 / (4 * tau**2)) * np.cos(carrier_w * s)
    off = 6 * tau; steps = int((2 * off + 300) / 0.5); inc = np.sum(pulse(np.arange(-4000, 4000, 1.0))**2)
    Ef = fdtd(N, steps, lambda s: pulse(s - off), x_ref=float(x_int), eps_inf=eps_profile, snaps=(steps,), **(extra or {}))[0][steps]
    return np.sum(Ef[: x_int - 100]**2) / inc

N9, x9 = 3200, 1600; idx = np.arange(N9)
eps_bare = np.where(idx >= x9, 2.25**2, 1.0)
eps_ar = np.where(idx >= x9 + 20, 2.25**2, np.where(idx >= x9, 1.5**2, 1.0))
R_bare = reflected_fraction(eps_bare, 2 * np.pi / 120, 100.0, N9, x9)
R_ar = reflected_fraction(eps_ar, 2 * np.pi / 120, 100.0, N9, x9)
close(R_bare, ((1 - 2.25) / 3.25)**2, 0.02, "bare substrate R"); assert R_ar < 0.01
WP = 0.1; plasma = dict(wp2=np.where(idx >= x9, WP**2, 0.0)); n_hi = np.sqrt(1 - (WP / 0.2)**2)
R_lo, R_hi = [reflected_fraction(np.ones(N9), wc_, 80.0, N9, x9, plasma) for wc_ in (0.05, 0.2)]
assert R_lo > 0.97 and abs(R_hi - ((1 - n_hi) / (1 + n_hi))**2) < 3e-3
print(f"[YT] AR coating: R bare = {R_bare:.4f} (theory {((1-2.25)/3.25)**2:.4f}), R with quarter-wave n = 1.5 layer = {R_ar:.5f}")
print(f"[YT] plasma w_p = {WP}: R(w = 0.05) = {R_lo:.4f} (below w_p: mirror), R(w = 0.2) = {R_hi:.5f} "
      f"(theory {((1-n_hi)/(1+n_hi))**2:.5f})")

# %% [markdown]
# ## Part 9: The numbers quoted in solutions.md

# %%
lam = 500e-9; E0a = 300.0                                             # A1
print(f"[A1] f = {C/lam:.4e} Hz, w = {2*np.pi*C/lam:.4e} rad/s, k = {2*np.pi/lam:.4e} rad/m, B0 = {E0a/C:.4e} T, I = {0.5*C*EPS0*E0a**2:.1f} W/m^2")
Ilas = 5e-3 / (np.pi * 0.5e-3**2)                                     # A2
print(f"[A2] I = {Ilas:.0f} W/m^2, E0 = {np.sqrt(2*Ilas/(C*EPS0)):.0f} V/m, E_rms = {np.sqrt(Ilas/(C*EPS0)):.0f} V/m, F = 2P/c = {2*5e-3/C:.3e} N")
nw = 1.33; print(f"[A4] water: R = {((nw-1)/(nw+1))**2:.5f}, Brewster {np.degrees(np.arctan(nw)):.2f} deg, critical {np.degrees(np.arcsin(1/nw)):.2f} deg")
for fq in (1e4, 1e6):                                                 # A5
    print(f"[A5] seawater sigma = 4 S/m, f = {fq:.0e} Hz: delta = {np.sqrt(2/(MU0*4*2*np.pi*fq)):.4f} m, "
          f"sigma/(eps w) = {4/(81*EPS0*2*np.pi*fq):.3e}")
Rs3, Ts3, Rp3, Tp3 = RT(1.0, 1.5, th45)                               # B3 (same numbers as Example 4.4)
print(f"[B3] 45 deg: degree of polarization of reflected light = {(Rs3-Rp3)/(Rs3+Rp3):.4f}; reflected fraction of unpolarized = {(Rs3+Rp3)/2:.5f}")
th60 = np.radians(60); rs60 = fresnel(1.5, 1.0, th60)[0]              # B4
print(f"[B4] glass->air 60 deg: d = {decay_length(1.5, 1.0, th60, 633):.1f} nm at 633 nm; r_s = {rs60:.4f}, |r_s| = {abs(rs60):.6f}, "
      f"phase = {np.degrees(np.angle(rs60)):.2f} deg")
fpi = wp_ion / 2 / np.pi; ng = np.sqrt(1 - (fpi / 20e6)**2)            # B5
print(f"[B5] 20 MHz in the N = 1e12 ionosphere: n = {ng:.4f}, v_g = {ng:.4f} c, v_p = {1/ng:.4f} c; "
      f"extra delay over 100 km = {1e5/C*(1/ng - 1)*1e6:.2f} us")
print(f"[B6] group delay of 1 m of X-band guide at 10 GHz: {1/vg*1e9:.3f} ns (vacuum {1/C*1e9:.3f} ns); single-mode band "
      f"{f10/1e9:.2f}-{fc(2,0)/1e9:.2f} GHz")
print(f"[B7] Z(vacuum, b/a = 3.5) = {Z0/(2*np.pi)*np.log(3.5):.2f} ohm; Z(eps_r = 2.25) = {Zcoax:.2f} ohm")
r1, r2 = (1 - 1.38) / 2.38, (1.38 - 1.5) / 2.88; rAR = (r1 - r2) / (1 - r1 * r2)          # B8, quarter-wave: e^{2i delta} = -1
assert np.isclose(rAR, (1.5 - 1.38**2) / (1.5 + 1.38**2))
print(f"[B8] MgF2 (n = 1.38) quarter-wave on n = 1.5 glass: r1 = {r1:.5f}, r2 = {r2:.5f}, r = {rAR:.5f}, R = {rAR**2:.4f} vs bare 0.0400; "
      f"thickness at 550 nm = {550/(4*1.38):.1f} nm")
for Fd in (5e-5, 2e-4):                                               # C1: dilute Lorentz gas, w0 = 1, gamma = 0.01
    kd = lambda w: w * np.sqrt(1 + Fd / (1 - w**2 - 0.01j * w))
    ng = ((kd(1 + 1e-6) - kd(1 - 1e-6)) / 2e-6).real                # c/v_g = d(Re k)/dw at w0
    close(ng, 1 - Fd / 0.01**2, 1e-2, "group index at resonance")
    print(f"[C1] F = {Fd:.0e}, gamma = 0.01: n_g(w0) = {ng:.4f} (dilute formula 1 - F/gamma^2 = {1 - Fd/1e-4:.2f}), "
          f"v_g = {1/ng:.3f} c; intensity absorption length 1/alpha = {1/(2*kd(1.0).imag):.1f} (theory c gamma/F = {0.01/Fd:.1f})")
thr = np.radians(np.linspace(41.6, 89.9, 200001)); nr = 1.51                       # C2: Fresnel rhomb
rsr, _, rpr, _, _, _ = fresnel(nr, 1.0, thr)
dphi = np.degrees(np.abs(np.angle(-rpr / rsr)))                      # -r_p: right-handed (p, s, k) basis per beam
sols = thr[np.where(np.diff(np.sign(dphi - 45)))[0]]
dachro = [np.degrees(abs(np.angle(-fresnel(nn1, 1.0, np.radians(54.62))[2] / fresnel(nn1, 1.0, np.radians(54.62))[0]))) for nn1 in (1.50, 1.52)]
print(f"[C2] n = 1.51 (critical {np.degrees(np.arcsin(1/nr)):.2f} deg): |phase(-r_p) - phase(r_s)| = 45 deg at {np.degrees(sols).round(2)} deg; "
      f"max {dphi.max():.2f} deg at {np.degrees(thr[np.argmax(dphi)]):.2f} deg; at 54.62 deg: n = 1.50 -> {dachro[0]:.2f}, n = 1.52 -> {dachro[1]:.2f} deg; "
      f"toy glass n(400 nm) = {n_toy(np.array([0.4]))[0].real:.4f}, n(700 nm) = {n_toy(np.array([0.7]))[0].real:.4f}")
td = lambda m, s0: 2 * np.sqrt(3) * m * s0**2 / HBAR                  # C3: sigma doubles
print(f"[C3] width doubling time: electron, sigma0 = 1 A: {td(ME, 1e-10):.3e} s; 1 g marble, sigma0 = 1 um: "
      f"{td(1e-3, 1e-6):.3e} s = {td(1e-3, 1e-6)/3.156e7:.2e} yr")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s.")
