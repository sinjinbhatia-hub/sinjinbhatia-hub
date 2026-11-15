# %% [markdown]
# # SM-1 Project: A Hard-Disk Gas, Pistons, Cycles, and a Hot Rod
#
# This script checks every number in the SM-1 notes and builds a small
# molecular-dynamics gas from scratch. Parts: (1) the worked examples of §4;
# (2) a 2D hard-disk gas: equipartition, relaxation to the Maxwell distribution
# (3.4), pressure from wall impulses vs N k_B T / A; (3) piston compressions:
# slow adiabatic (gamma = 2 in 2D), slow isothermal, fast (irreversible);
# (4) the Otto cycle on the P-V plane; (5) path dependence of Q and W;
# (6) the heat equation on a rod with Fourier-law fluxes; (7) numbers used in
# the solutions. "Your turn" and a reference solution for B8 come last.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from scipy import constants as C
from scipy.integrate import quad
from scipy.spatial import cKDTree

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()
kB, NA, R = C.k, C.N_A, C.R
ATM = 101325.0

# %% [markdown]
# ## Part 1: Worked examples (§4)

# %%
# Example 4.1: air at STP.
T, P, M = 273.15, ATM, 28.97e-3
m = M / NA
nden = P / (kB * T)
rho = nden * m
v2 = 3 * P / rho                                 # eq. (3.2)
vrms, vp, vbar = np.sqrt(v2), np.sqrt(2 * kB * T / m), np.sqrt(8 * kB * T / (np.pi * m))
print(f"Ex 4.1: m={m:.4g} kg, n={nden:.4g} m^-3, rho={rho:.4g} kg/m^3, <v^2>={v2:.4g}, "
      f"v_rms={vrms:.1f}, v_p={vp:.1f}, v_bar={vbar:.1f} m/s, KE={1.5*kB*T:.4g} J "
      f"= {1.5*kB*T/C.e:.4f} eV, molar volume={R*T/P*1e3:.2f} L, sound={vrms*np.sqrt(1.4/3):.0f} m/s")
assert abs(nden - 2.687e25) / 2.687e25 < 1e-3 and abs(rho - 1.292) < 1e-3
assert round(vrms) == 485 and round(vp) == 396 and round(vbar) == 447
assert np.isclose(vrms, np.sqrt(3 * kB * T / m))

# Example 4.2: C_P - C_V.
CV_air, CP_air = 2.5 * R, 3.5 * R
print(f"Ex 4.2: air C_V={CV_air:.2f}, C_P={CP_air:.2f}, gamma={CP_air/CV_air:.3f}; tabulated "
      f"1.005*28.97={1.005*28.97:.2f}; He C_V={1.5*R:.2f}, C_P={2.5*R:.2f}; dV={R/ATM:.3g} m^3")
assert np.isclose(CP_air - CV_air, R) and abs(1.005 * 28.97 - CP_air) < 0.05

# Example 4.3: adiabatic compression by 10.
g = 1.4
Tf, Pf = 300 * 10 ** (g - 1), 10 ** g
W_ad = CV_air * (Tf - 300)
W_iso10 = R * 300 * np.log(10)
print(f"Ex 4.3: T_f={Tf:.1f} K, P_f={Pf:.2f} atm, W={W_ad:.0f} J/mol, isothermal W={W_iso10:.0f} J/mol")
assert round(Tf, 1) == 753.6 and round(Pf, 2) == 25.12 and round(W_ad) == 9428 and round(W_iso10) == 5743
assert np.isclose(Pf / 10, Tf / 300)


# Example 4.5: Otto cycle, leg by leg.
def otto(r=8.0, V1=5e-4, P1=1e5, T1=300.0, T3=1800.0, gam=1.4):
    """Air-standard Otto cycle (notes Example 4.5). Returns states and per-leg (W_on, Q) in J."""
    n = P1 * V1 / (R * T1)
    nCV = n * R / (gam - 1)
    T2 = T1 * r ** (gam - 1); P2 = P1 * r ** gam
    P3 = P2 * T3 / T2
    T4 = T3 / r ** (gam - 1); P4 = P3 / r ** gam
    legs = {"12": (nCV * (T2 - T1), 0.0), "23": (0.0, nCV * (T3 - T2)),
            "34": (nCV * (T4 - T3), 0.0), "41": (0.0, nCV * (T1 - T4))}
    states = dict(V=[V1, V1 / r, V1 / r, V1], P=[P1, P2, P3, P4], T=[T1, T2, T3, T4], n=n, nCV=nCV)
    return states, legs


st, legs = otto()
Wnet = sum(w for w, _ in legs.values()); Qnet = sum(q for _, q in legs.values())
print(f"Ex 4.5: n={st['n']:.5f} mol, nC_V={st['nCV']:.4f} J/K, T={np.round(st['T'], 1)}, "
      f"P(kPa)={np.round(np.array(st['P'])/1e3, 1)}")
print("        legs (W_on, Q) J:", {k: (round(w, 1), round(q, 1)) for k, (w, q) in legs.items()},
      f"net W={Wnet:.1f}, net Q={Qnet:.1f}, ratio={-Wnet/legs['23'][1]:.4f}, power at 3000 rpm="
      f"{-Wnet*25/1e3:.2f} kW")
assert np.isclose(Wnet + Qnet, 0) and round(-Wnet, 1) == 261.4 and round(legs["23"][1], 1) == 462.8
assert np.isclose(-Wnet / legs["23"][1], 1 - 8 ** -0.4) and np.isclose(st["P"][2], 4.8e6)

# Example 4.6: cyclic rule, symbolically, for the ideal gas and for van der Waals.
Ps, Vs, Ts, Nk, a_, b_, n_ = sp.symbols("P V T Nk a b n", positive=True)
for name, Pexpr in [("ideal", Nk * Ts / Vs), ("vdW", n_ * sp.Symbol("R") * Ts / (Vs - n_ * b_) - a_ * n_**2 / Vs**2)]:
    dPdT_V = sp.diff(Pexpr, Ts)
    dPdV_T = sp.diff(Pexpr, Vs)
    dTdV_P = -dPdV_T / dPdT_V             # implicit differentiation of P(T,V) = const
    dVdP_T = 1 / dPdV_T                   # reciprocity (3.17)
    prod = sp.simplify(dPdT_V * dTdV_P * dVdP_T)
    print(f"Ex 4.6: cyclic product for {name} gas = {prod}")
    assert prod == -1
beta_w, kap_w = 2.07e-4, 4.59e-10
print(f"Ex 4.6: water (dP/dT)_V = {beta_w/kap_w:.3g} Pa/K = {beta_w/kap_w/ATM:.2f} atm/K")
assert round(beta_w / kap_w / ATM, 2) == 4.45

# §3.7 latent heat and §3.10 conduction, §3.11 CO2 critical point.
Wvap = 1 / 18.015e-3 * R * 373.15
print(f"§3.7: P dV = {Wvap/1e6:.3f} MJ, Delta U = {(2.257e6 - Wvap)/1e6:.3f} MJ; §3.10 brick wall q = {0.7*20/0.2:.0f} W/m^2")
assert round(Wvap / 1e6, 3) == 0.172 and round((2.257e6 - Wvap) / 1e6, 3) == 2.085
aC, bC = 0.364, 4.267e-5
Tc, Pc, Vc = 8 * aC / (27 * R * bC), aC / (27 * bC**2), 3 * bC
print(f"§3.11: CO2 T_c={Tc:.1f} K, P_c={Pc/1e6:.2f} MPa, V_c={Vc:.3g} m^3/mol, Z_c={Pc*Vc/(R*Tc):.3f}")
assert round(Tc) == 304 and round(Pc / 1e6, 2) == 7.40

# %% [markdown]
# ## Part 2: A 2D hard-disk gas
#
# Units: disk mass m = 1, k_B = 1, initial k_B T = 1 (so each disk has energy 1
# on average in 2D, eq. (3.3) with d = 2). Disks of diameter SIG move in straight
# lines for a time step DT; overlapping pairs that are approaching exchange the
# normal component of their relative velocity (an elastic collision of equal
# masses). Walls reflect. The impulse delivered to the walls gives the pressure
# (force per unit length), eq. (3.2) in 2D: P A = N k_B T.

# %%
N, L, SIG, DT = 400, 20.0, 0.1, 0.01
RAD = SIG / 2


def init_gas(seed, mode):
    """Disks on a jittered grid. mode='x': all speeds sqrt(2) along +-x (far from equilibrium);
    mode='iso': speed sqrt(2), random directions. Both have total energy N (k_B T = 1)."""
    r = np.random.default_rng(seed)
    k = int(np.ceil(np.sqrt(N)))
    gx = (np.arange(k) + 0.5) * L / k
    X, Y = np.meshgrid(gx, gx)
    pos = np.c_[X.ravel(), Y.ravel()][:N] + r.uniform(-0.3, 0.3, (N, 2))
    if mode == "x":
        vel = np.zeros((N, 2)); vel[:, 0] = r.choice([-1.0, 1.0], N) * np.sqrt(2)
    else:
        th = r.uniform(0, 2 * np.pi, N); vel = np.sqrt(2) * np.c_[np.cos(th), np.sin(th)]
    return pos, vel


def collide(pos, vel):
    """Elastic equal-mass collisions for every overlapping, approaching pair (processed in turn
    so that energy and momentum are conserved exactly even when a disk is in two pairs)."""
    for a, b in cKDTree(pos).query_pairs(SIG, output_type="ndarray"):
        d = pos[a] - pos[b]; dv = vel[a] - vel[b]; p = d @ dv
        if p < 0:
            dvn = (p / (d @ d)) * d
            vel[a] -= dvn; vel[b] += dvn


def thermal_kick(vel, mask, axis, sign, Tb, r):
    """Thermalizing wall: re-draw velocities of disks hitting it from the wall-flux distribution
    at temperature Tb (normal component Rayleigh, tangential Gaussian). Returns impulse delivered."""
    k = mask.sum()
    if k == 0:
        return 0.0
    old = vel[mask, axis].copy()
    vn = np.sqrt(-2 * Tb * np.log(r.uniform(size=k)))
    vel[mask, axis] = sign * vn
    vel[mask, 1 - axis] = r.normal(0, np.sqrt(Tb), k)
    return np.abs(old).sum() + vn.sum()


def step(pos, vel, Lx, uw, thermal=None, r=None):
    """Advance one DT. Right wall (piston) at x = Lx moving with velocity uw (uw < 0 compresses).
    A disk hitting it reflects in the piston frame: v_x -> 2 uw - v_x. Returns (wall impulse on
    the fixed walls, impulse on piston, work done on the gas by the piston)."""
    pos += vel * DT
    imp = 0.0
    for axis, lo, hi in ((0, True, False), (1, True, True)):
        walls = [(pos[:, axis] < RAD) & (vel[:, axis] < 0, +1)] if lo else []
        if hi:
            walls.append(((pos[:, axis] > L - RAD) & (vel[:, axis] > 0), -1))
        for m_, sgn in walls:
            if thermal is not None:
                imp += thermal_kick(vel, m_, axis, sgn, thermal, r)
            else:
                imp += 2 * np.abs(vel[m_, axis]).sum(); vel[m_, axis] *= -1
    hit = (pos[:, 0] > Lx - RAD) & (vel[:, 0] > uw)
    e0 = 0.5 * (vel[hit] ** 2).sum()
    imp_p = 2 * (vel[hit, 0] - uw).sum()
    vel[hit, 0] = 2 * uw - vel[hit, 0]
    work = 0.5 * (vel[hit] ** 2).sum() - e0
    collide(pos, vel)
    return imp, imp_p, work


# Run 1: start with all motion along x; watch equipartition and the speed distribution.
pos, vel = init_gas(1, "x")
D2 = lambda v, T: 1 - np.exp(-v**2 / (2 * T))          # 2D Maxwell CDF, from (3.4)
speeds0 = np.linalg.norm(vel, axis=1).copy()
ts, vx2, vy2, snaps, imp_tot, t_meas = [], [], [], [], 0.0, 0.0
T_RELAX, T_END = 60.0, 160.0
for s in range(int(T_END / DT)):
    imp, imp_p, _ = step(pos, vel, L, 0.0)
    t = (s + 1) * DT
    if t > T_RELAX:
        imp_tot += imp + imp_p; t_meas += DT
        if s % 100 == 0:
            snaps.append(np.linalg.norm(vel, axis=1))
    if s % 20 == 0:
        ts.append(t); vx2.append((vel[:, 0] ** 2).mean()); vy2.append((vel[:, 1] ** 2).mean())
ts, vx2, vy2 = map(np.array, (ts, vx2, vy2))
kT = 0.5 * (vel ** 2).sum(1).mean()
assert abs(kT - 1) < 1e-9, "energy must be conserved exactly by elastic collisions"
late = ts > T_RELAX
ratio = vx2[late].mean() / vy2[late].mean()
spd = np.concatenate(snaps)
grid = np.sort(spd)
ks = np.max(np.abs(np.arange(1, grid.size + 1) / grid.size - D2(grid, kT)))
ks0 = np.max(np.abs(np.arange(1, N + 1) / N - D2(np.sort(speeds0), kT)))
print(f"Part 2: <v_x^2>/<v_y^2> after relaxation = {ratio:.3f}; KS distance to 2D Maxwell: "
      f"start {ks0:.3f}, relaxed {ks:.4f}")
assert abs(ratio - 1) < 0.05 and ks0 > 0.3 and ks < 0.03

# Pressure: total impulse / (time * perimeter seen by disk centres). Disk centres live in
# (L - SIG)^2. Finite disk size adds the second virial correction B2 n = (pi SIG^2 / 2) n.
Aeff = (L - SIG) ** 2
P_meas = imp_tot / (t_meas * 4 * (L - SIG))
P_ideal = N * kT / Aeff
B2 = np.pi * SIG**2 / 2
P_vir = P_ideal * (1 + B2 * N / Aeff)
# Statistical error from 8 independent short runs at equilibrium.
reps = []
for seed in range(8):
    p2, v2_ = init_gas(100 + seed, "iso")
    it, tt = 0.0, 0.0
    for s in range(3000):
        a1, a2_, _ = step(p2, v2_, L, 0.0)
        if s >= 1000:
            it += a1 + a2_; tt += DT
    reps.append(it / (tt * 4 * (L - SIG)))
P_err = np.std(reps, ddof=1) / np.sqrt(len(reps))
P_rep = np.mean(reps)
print(f"Part 2: P_meas={P_meas:.5f}, replicate mean={P_rep:.5f} +- {P_err:.5f}, ideal NkT/A={P_ideal:.5f}, "
      f"with virial correction={P_vir:.5f}")
assert abs(P_rep - P_vir) < 3 * P_err + 0.005 * P_vir
assert abs(P_meas / P_ideal - 1) < 0.05

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(ts, vx2, label=r"$\langle v_x^2\rangle$"); ax[0].plot(ts, vy2, label=r"$\langle v_y^2\rangle$")
ax[0].axhline(1, color="k", lw=0.8, ls="--"); ax[0].set_xlabel("time"); ax[0].legend()
ax[0].set_title("Equipartition from an all-x start")
vv = np.linspace(0, 5, 200)
ax[1].hist(spd, bins=60, density=True, alpha=0.6, label="simulation (t > 60)")
ax[1].plot(vv, vv / kT * np.exp(-vv**2 / (2 * kT)), "k", label="2D Maxwell (3.4)")
ax[1].axvline(np.sqrt(2), color="r", ls=":", label="initial speed")
ax[1].set_xlabel("speed"); ax[1].legend(); ax[1].set_title("Speed distribution")
fig.tight_layout(); fig.savefig(FIG / "hard_disk_relaxation.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 3: Piston compressions
#
# The right wall moves inward at speed u from Lx = L to Lx = L/2. Slow and insulated:
# eq. (3.15) with gamma = 2 says T A = const. Slow with thermalizing fixed walls at
# T = 1: W = N k_B T ln(A_i/A_f), eq. (3.14). Fast: the gas piles up on the piston
# and ends hotter than the adiabat predicts (quasistatic assumption of (3.8) broken).


# %%
def compress(u, thermal=None, seed=7, Lx_f=L / 2):
    """Move the piston from L to Lx_f at speed u. Returns final k_B T (averaged over the last
    stretch), total work on the gas by the piston, and the (A_eff, T) trace."""
    r = np.random.default_rng(seed)
    p_, v_ = init_gas(seed, "iso")
    for _ in range(1000):                         # equilibrate in the full box
        step(p_, v_, L, 0.0, thermal, r)
    Lx, W, trace = L, 0.0, []
    nsteps = int(round((L - Lx_f) / (u * DT)))
    for s in range(nsteps):
        Lx = L - u * DT * (s + 1)
        _, _, w = step(p_, v_, Lx, -u, thermal, r); W += w
        if s % max(1, nsteps // 200) == 0:
            trace.append(((Lx - SIG) * (L - SIG), 0.5 * (v_**2).sum(1).mean()))
    Ts_ = []
    for _ in range(1500):                         # settle with the piston stopped
        step(p_, v_, Lx_f, 0.0, thermal, r); Ts_.append(0.5 * (v_**2).sum(1).mean())
    return np.mean(Ts_[500:]), W, np.array(trace)


Ai, Af = (L - SIG) ** 2, (L / 2 - SIG) * (L - SIG)
T_ad, W_ad_sim, tr_ad = compress(0.02)
T_fast, W_fast, tr_fast = compress(1.0)
T_iso, W_iso_sim, tr_iso = compress(0.02, thermal=1.0)
T_pred = (Ai / Af) * np.exp(B2 * N * (1 / Af - 1 / Ai))     # adiabat including B2 (solutions B1(c))
gam_eff = 1 + np.log(T_ad) / np.log(Ai / Af)
W_iso_th = N * 1.0 * (np.log(Ai / Af) + B2 * N * (1 / Af - 1 / Ai))
print(f"Part 3: slow adiabatic T_f={T_ad:.3f} (ideal gamma=2: {Ai/Af:.3f}, with B2: {T_pred:.3f}), "
      f"gamma_eff={gam_eff:.3f}, W/N={W_ad_sim/N:.3f}")
print(f"        fast piston T_f={T_fast:.3f}, W/N={W_fast/N:.3f}  (extra heating {T_fast-T_ad:.3f})")
print(f"        slow isothermal T={T_iso:.3f}, W/(N kT)={W_iso_sim/N:.3f} vs ln(Ai/Af)={np.log(Ai/Af):.3f} "
      f"(with B2: {W_iso_th/N:.3f})")
assert abs(gam_eff - 2) < 0.1 and abs(T_ad / T_pred - 1) < 0.03
assert np.isclose(W_ad_sim / N, T_ad - 1, atol=0.03)        # first law, Q = 0: W = Delta U
assert T_fast > 1.1 * T_ad
assert abs(T_iso - 1) < 0.06 and abs(W_iso_sim / W_iso_th - 1) < 0.08

fig, ax = plt.subplots(figsize=(6, 4.5))
Agrid = np.linspace(Af, Ai, 100)
for tr, lab in [(tr_ad, "slow, insulated"), (tr_fast, "fast, insulated"), (tr_iso, "slow, thermal walls")]:
    ax.loglog(tr[:, 0], N * tr[:, 1] / tr[:, 0], ".", ms=3, label=lab)
ax.loglog(Agrid, N / Ai * (Ai / Agrid) ** 2, "k--", label=r"$PA^2$ = const ($\gamma=2$)")
ax.loglog(Agrid, N / Agrid, "k:", label=r"$PA$ = const (isotherm)")
ax.set_xlabel("area A"); ax.set_ylabel(r"$P = Nk_BT/A$"); ax.legend(fontsize=8)
ax.set_title("2D hard-disk gas under a piston")
fig.tight_layout(); fig.savefig(FIG / "piston_compressions.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 4: The Otto cycle on the P-V plane (Example 4.5)


# %%
def loop_work(segments):
    """W_on = -sum over segments of the integral of P dV, eq. (3.8). Each segment is (P(V), V_a, V_b)."""
    return -sum(quad(Pf_, Va, Vb)[0] for Pf_, Va, Vb in segments)


V1, V2 = st["V"][0], st["V"][1]
P1, P2, P3, P4 = st["P"]
otto_segs = [(lambda V: P1 * (V1 / V) ** 1.4, V1, V2), (lambda V: 0.0, V2, V2),
             (lambda V: P3 * (V2 / V) ** 1.4, V2, V1), (lambda V: 0.0, V1, V1)]
W_num = loop_work(otto_segs)
print(f"Part 4: Otto -oint P dV numerically = {W_num:.2f} J (analytic {Wnet:.2f} J)")
assert np.isclose(W_num, Wnet, rtol=1e-8)
fig, ax = plt.subplots(figsize=(5.5, 4))
Vc_ = np.linspace(V2, V1, 200)
ax.plot(Vc_ * 1e6, P1 * (V1 / Vc_) ** 1.4 / 1e3, label="1-2 adiabatic compression")
ax.plot([V2 * 1e6] * 2, [P2 / 1e3, P3 / 1e3], label="2-3 heat in")
ax.plot(Vc_ * 1e6, P3 * (V2 / Vc_) ** 1.4 / 1e3, label="3-4 adiabatic expansion")
ax.plot([V1 * 1e6] * 2, [P4 / 1e3, P1 / 1e3], label="4-1 heat out")
ax.set_xlabel("V (cm$^3$)"); ax.set_ylabel("P (kPa)"); ax.legend(fontsize=8); ax.set_title("Otto cycle, r = 8")
fig.tight_layout(); fig.savefig(FIG / "otto_cycle.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5: Path dependence (Example 4.4)
#
# Monatomic ideal gas, U = (3/2) P V. Units kPa and L, so products are joules.
# Each path is parametrized by s in [0, 1]; we integrate dU, dW = -P dV, and dQ = dU - dW.


# %%
def path_integrals(Pfun, Vfun, n=20001):
    """Return (Delta U, W, Q, integral of dQ/T * T1) along (P(s), V(s)), trapezoid rule."""
    s = np.linspace(0, 1, n)
    Pp, Vv = Pfun(s), Vfun(s)
    U = 1.5 * Pp * Vv
    dV, dU = np.diff(Vv), np.diff(U)
    Pm = 0.5 * (Pp[1:] + Pp[:-1])
    dW = -Pm * dV
    dQ = dU - dW
    Tm = 0.5 * (Pp[1:] * Vv[1:] + Pp[:-1] * Vv[:-1]) / 1000.0   # T/T1, since P1 V1 = 1000 J
    return U[-1] - U[0], dW.sum(), dQ.sum(), (dQ / Tm).sum()


pw = lambda a, b: (lambda s: np.where(s < 0.5, a(2 * s), b(2 * s - 1)))
paths = {
    "a": (pw(lambda x: 100 + 0 * x, lambda x: 100 + 100 * x), pw(lambda x: 10 + 10 * x, lambda x: 20 + 0 * x)),
    "b": (pw(lambda x: 100 + 100 * x, lambda x: 200 + 0 * x), pw(lambda x: 10 + 0 * x, lambda x: 10 + 10 * x)),
    "c": (lambda s: 100 + 100 * s, lambda s: 10 + 10 * s),
}
res = {k: path_integrals(*v) for k, v in paths.items()}
for k, (dU, W, Q, S) in res.items():
    print(f"Part 5: path {k}: dU={dU:.1f} J, W={W:.1f} J, Q={Q:.1f} J, T1*int(dQ/T)={S:.1f} J")
for k, (Wexp, Qexp) in {"a": (-1000, 5500), "b": (-2000, 6500), "c": (-1500, 6000)}.items():
    assert np.isclose(res[k][0], 4500) and np.isclose(res[k][1], Wexp, atol=0.5) and np.isclose(res[k][2], Qexp, atol=0.5)
S_exact = 1000 * (1.5 * np.log(4) + np.log(2))
assert all(np.isclose(r_[3], S_exact, rtol=1e-4) for r_ in res.values())   # dQ/T is exact (§3.4)
print(f"Part 5: cycle b then a reversed: W_net={res['b'][1]-res['a'][1]:.1f} J, Q_net={res['b'][2]-res['a'][2]:.1f} J; "
      f"T1*Delta S={S_exact:.1f} J on every path")

# %% [markdown]
# ## Part 6: Conduction in an aluminium rod (§3.10, F2 §3.12)
#
# Explicit finite differences for (3.23). Left end held at 100 C, right end at 0 C,
# initially 0 C. Check: (i) energy stored = time integral of Fourier flux in minus out;
# (ii) the profile matches the sine-series solution; (iii) the steady flux k*100/L.

# %%
kth, rho_al, c_al, Lr = 237.0, 2700.0, 897.0, 0.5
alpha = kth / (rho_al * c_al)
nx = 51
x = np.linspace(0, Lr, nx); dx = x[1] - x[0]
dt = 0.4 * dx**2 / alpha
Tr = np.zeros(nx); Tr[0] = 100.0
E_in, t_r, saved = 0.0, 0.0, None
t_snap = 500.0
for s in range(int(8000 / dt)):
    q_in = -kth * (-3 * Tr[0] + 4 * Tr[1] - Tr[2]) / (2 * dx)        # Fourier's law (3.22), 2nd order
    q_out = -kth * (3 * Tr[-1] - 4 * Tr[-2] + Tr[-3]) / (2 * dx)
    Tr[1:-1] += alpha * dt / dx**2 * (Tr[2:] - 2 * Tr[1:-1] + Tr[:-2])
    q_in2 = -kth * (-3 * Tr[0] + 4 * Tr[1] - Tr[2]) / (2 * dx)
    q_out2 = -kth * (3 * Tr[-1] - 4 * Tr[-2] + Tr[-3]) / (2 * dx)
    E_in += 0.5 * dt * ((q_in - q_out) + (q_in2 - q_out2))
    t_r += dt
    if saved is None and t_r >= t_snap:
        saved, t_saved = Tr.copy(), t_r
stored = rho_al * c_al * np.trapezoid(Tr, x) if hasattr(np, "trapezoid") else rho_al * c_al * np.trapz(Tr, x)
nn = np.arange(1, 400)[:, None]
series = 100 * (1 - x / Lr) - (200 / (nn * np.pi) * np.sin(nn * np.pi * x / Lr)
                               * np.exp(-alpha * (nn * np.pi / Lr) ** 2 * t_saved)).sum(0)
q_ss = kth * 100 / Lr
print(f"Part 6: alpha={alpha:.3g} m^2/s, stored={stored/1e6:.4f} MJ/m^2 vs flux integral={E_in/1e6:.4f} MJ/m^2, "
      f"max |FD - series| at t={t_saved:.0f} s = {np.max(np.abs(saved-series)[1:-1]):.3f} K, "
      f"end flux={q_in2:.0f} vs k*100/L={q_ss:.0f} W/m^2")
assert abs(stored / E_in - 1) < 0.01 and np.max(np.abs(saved - series)[1:-1]) < 0.5
assert abs(q_in2 / q_ss - 1) < 0.02
fig, ax = plt.subplots(figsize=(5.5, 4))
ax.plot(x, saved, "o", ms=3, label=f"finite differences, t={t_saved:.0f} s")
ax.plot(x, series, label="sine series (F2)"); ax.plot(x, Tr, label="t = 8000 s (near steady)")
ax.set_xlabel("x (m)"); ax.set_ylabel("T (C)"); ax.legend(); ax.set_title("Aluminium rod")
fig.tight_layout(); fig.savefig(FIG / "rod_heat_equation.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: Numbers used in the solutions

# %%
print(f"A1: v_rms He={np.sqrt(3*R*300/4.0026e-3):.0f}, N2={np.sqrt(3*R*300/28.014e-3):.0f} m/s; "
      f"A3: W=+{R*300*np.log(2):.0f} J")
# B4: Diesel cycle per mole, r = 18, cutoff ratio 2.
r, rc, gm = 18.0, 2.0, 1.4
CV, CP = R / (gm - 1), gm * R / (gm - 1)
D1 = 300.0; D2_ = D1 * r ** (gm - 1); D3 = 2 * D2_; D4 = D3 * (rc / r) ** (gm - 1)
dW = {"12": CV * (D2_ - D1), "23": -R * (D3 - D2_), "34": CV * (D4 - D3), "41": 0.0}
dQ = {"12": 0.0, "23": CP * (D3 - D2_), "34": 0.0, "41": CV * (D1 - D4)}
eta_D = -sum(dW.values()) / dQ["23"]
eta_formula = 1 - r ** (1 - gm) * (rc**gm - 1) / (gm * (rc - 1))
print(f"B4: T2={D2_:.1f}, P2={100*r**gm:.0f} kPa, T3={D3:.1f}, T4={D4:.1f} K; W={ {k: round(v) for k, v in dW.items()} }, "
      f"Q={ {k: round(v) for k, v in dQ.items()} }; net W={sum(dW.values()):.0f}, ratio={eta_D:.4f}")
assert np.isclose(eta_D, eta_formula) and np.isclose(sum(dW.values()) + sum(dQ.values()), 0)
# B5: rectangle cycle, diatomic, U = 5/2 PV, corners (V in L, P in kPa).
cor = [(1, 100), (1, 300), (4, 300), (4, 100), (1, 100)]
tot_W = tot_Q = Qin = 0.0
for (Va, Pa), (Vb, Pb) in zip(cor[:-1], cor[1:]):
    W_ = -Pa * (Vb - Va) if Pa == Pb else 0.0
    Q_ = 2.5 * (Pb * Vb - Pa * Va) - W_
    tot_W += W_; tot_Q += Q_; Qin += max(Q_, 0)
    print(f"B5: leg ({Va},{Pa})->({Vb},{Pb}): dU={2.5*(Pb*Vb-Pa*Va):.0f}, W={W_:.0f}, Q={Q_:.0f} J")
print(f"B5: net W={tot_W:.0f}, net Q={tot_Q:.0f}, Q_in={Qin:.0f}, ratio={-tot_W/Qin:.4f}")
assert tot_W == -600 and tot_Q == 600 and Qin == 3650
print(f"B6: wrong adiabat T_f={300*2**(-2/3):.1f} K, wrong dU={1.5*R*(300*2**(-2/3)-300):.0f} J; "
      f"B7: Delta S(Ex 4.4, n=1)={S_exact/(1000/R):.2f} J/K with T1={1000/R:.2f} K")
print(f"C1: CO2 V_c={Vc*1e6:.1f} cm^3/mol, Z_c=3/8={3/8}; C2: g/c_p={9.81/1005*1e3:.2f} K/km")
mb, Ab, Vb_ = 0.0165, np.pi * 0.008**2, 0.010
Pb = ATM + mb * 9.81 / Ab
w2 = 1.4 * Pb * Ab**2 / (mb * Vb_)
print(f"C3: A={Ab:.3g} m^2, P={Pb:.0f} Pa, omega={np.sqrt(w2):.2f} rad/s, period={2*np.pi/np.sqrt(w2):.3f} s")

# %% [markdown]
# ## Your turn
#
# 1. **Van der Waals isotherms for CO2** (SM-6 preview). Plot P(V) from (3.25) for
#    T = 280, 304, 330 K, V from 6e-5 to 6e-4 m^3/mol. Find the critical point
#    numerically: solve dP/dV = d2P/dV2 = 0 with `scipy.optimize.fsolve` and compare
#    with T_c = 8a/(27Rb). Where does dP/dV > 0, and why is that unphysical?
# 2. **The Diesel cycle** (exercise B8): build the four segments for `loop_work` and
#    check the B4 numbers. A reference solution follows; write yours first.
# 3. **3D hard spheres.** Copy Part 2 to 3D (positions and velocities with three
#    columns, six walls). Check P V = N k_B T with <(1/2) m v^2> = (3/2) k_B T, and
#    repeat the slow insulated compression: T V^(2/3) should be constant (gamma = 5/3).
# 4. Make the disks 4x larger (SIG = 0.4). How much does P/P_ideal grow? Compare with
#    1 + B2 N/A. When does the second-virial estimate start to fail?

# %% [markdown]
# ## Reference solution for B8 (write your own first)

# %%
V1d = 1.0; V2d = V1d / r; V3d = rc * V2d
P1d = R * D1 / V1d; P2d = P1d * r ** gm; P4d = R * D4 / V1d
diesel = [(lambda V: P1d * (V1d / V) ** gm, V1d, V2d), (lambda V: P2d, V2d, V3d),
          (lambda V: P2d * (V3d / V) ** gm, V3d, V1d), (lambda V: 0.0, V1d, V1d)]
W_d = loop_work(diesel)
print(f"B8: Diesel -oint P dV = {W_d:.1f} J/mol (B4: {sum(dW.values()):.1f})")
assert np.isclose(W_d, sum(dW.values()), rtol=1e-8) and np.isclose(P4d, P2d * (V3d / V1d) ** gm)

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s. Figures in {FIG}")
