# %% [markdown]
# # SM-2 Project: A Carnot engine you can run, entropy bookkeeping, potentials, and phase coexistence
#
# Parts: (1) a particle-level Carnot engine (1D ideal gas, thermal walls, a slow
# piston) whose measured efficiency is compared with 1 - Tc/Th; (2) the numbers
# of the worked examples; (3) entropy production between two bodies; (4) Maxwell
# relations, C_P - C_V and the Joule-Thomson coefficient checked with SymPy;
# (5) the Legendre transform U(S) -> F(T) done numerically; (6) Maxwell's
# equal-area construction for van der Waals CO2; (7) Clausius-Clapeyron for
# water; (8) the van der Waals Joule-Thomson inversion curve; (9) exercise numbers.
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from pathlib import Path
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
R = 8.314462618  # J/(mol K)

# %% [markdown]
# ## Part 1: a Carnot engine made of particles (notes §3.1, §7)
#
# N non-interacting particles (m = k_B = 1) move on [0, L]. The wall at x = 0 is
# either specular (adiabatic legs) or a *thermal wall* at temperature T: a particle
# that hits it leaves with a fresh speed drawn from the wall-flux distribution
# p(v) ∝ v exp(-v²/2T). That energy change is heat. The wall at x = L is a piston
# moving at speed ±u; a bounce gives v -> 2u_p - v, and that energy change is work.
# A 1D monatomic gas has f = 1, γ = 3, so an adiabat is T L² = const (SM-1's
# T V^(γ-1) = const). The cycle: hot isotherm L1->L2, adiabat L2->L3, cold
# isotherm L3->L4, adiabat L4->L1.


# %%
def run_leg(x, v, L0, L1, T_wall, u=0.001, dt=0.1, record=None):
    """Move the piston from L0 to L1. Returns (x, v, Q_in, W_on): heat into the gas
    at the thermal wall and work done on the gas by the piston."""
    n = int(np.ceil(abs(L1 - L0) / (u * dt)))
    up = (L1 - L0) / (n * dt)
    L, Q, W = L0, 0.0, 0.0
    for i in range(n):
        L += up * dt
        x += v * dt
        hit = x < 0
        if hit.any():
            if T_wall is None:
                x[hit] = -x[hit]; v[hit] = -v[hit]
            else:
                t_rem = x[hit] / v[hit]                       # time since the bounce (> 0)
                s = np.sqrt(-2 * T_wall * np.log(1 - rng.random(hit.sum())))
                Q += 0.5 * np.sum(s**2 - v[hit]**2)
                v[hit] = s; x[hit] = np.minimum(s * t_rem, L)
        hit = x > L
        if hit.any():
            vn = 2 * up - v[hit]
            W += 0.5 * np.sum(vn**2 - v[hit]**2)
            v[hit] = vn; x[hit] = 2 * L - x[hit]
        if record is not None and i % 50 == 0:
            record.append((L, np.sum(v**2) / L))              # 1D pressure P = 2E/L
    return x, v, Q, W


def carnot_engine(Th, Tc, N=20000, L1=1.0, L2=1.5, u=0.001):
    """Run one Carnot cycle; return measured efficiency, Q_h, Q_c, W_by, and the P-L path."""
    r = np.sqrt(Th / Tc)                                      # adiabat: T L^2 = const
    L3, L4 = L2 * r, L1 * r
    x = rng.uniform(0, L1, N); v = rng.normal(0, np.sqrt(Th), N)
    E0, path = 0.5 * np.sum(v**2), []
    x, v, Qh, W1 = run_leg(x, v, L1, L2, Th, u=u, record=path)
    x, v, Q2, W2 = run_leg(x, v, L2, L3, None, u=u, record=path)
    T_after_adiabat = np.mean(v**2)
    x, v, Qc, W3 = run_leg(x, v, L3, L4, Tc, u=u, record=path)
    x, v, Q4, W4 = run_leg(x, v, L4, L1, None, u=u, record=path)
    E1 = 0.5 * np.sum(v**2)
    W_on = W1 + W2 + W3 + W4
    assert Q2 == 0 and Q4 == 0
    assert abs((E1 - E0) - (Qh + Qc + W_on)) < 1e-6 * E0      # first law, exactly
    return -W_on / Qh, Qh, Qc, -W_on, T_after_adiabat, np.array(path)


fig, ax = plt.subplots(1, 2, figsize=(11, 4))
etas = []
for Tc in (0.5, 0.7):
    eta, Qh, Qc, Wby, Tad, path = carnot_engine(1.0, Tc)
    etas.append((Tc, eta))
    print(f"Carnot engine Th=1, Tc={Tc}: eta measured {eta:.4f}, 1-Tc/Th = {1 - Tc:.4f}, "
          f"T after hot->cold adiabat {Tad:.4f}, Q_c/Q_h = {Qc / Qh:.4f}")
    assert abs(Tad - Tc) < 0.02 * Tc                          # T L^2 = const along the adiabat
    assert abs(eta - (1 - Tc)) < 0.03
    ax[0].plot(path[:, 0], path[:, 1], label=f"$T_c$={Tc}")
    # entropy bookkeeping per notes (3.4): Q_h/T_h + Q_c/T_c ≈ 0
    print(f"   Q_h/T_h + Q_c/T_c = {Qh / 1.0 + Qc / Tc:.2f} (vs Q_h/T_h = {Qh:.1f})")
    assert abs(Qh + Qc / Tc) < 0.05 * Qh
eta_fast = carnot_engine(1.0, 0.5, u=0.01)[0]
print(f"Same cycle with the piston 10x faster (u = 0.01): eta = {eta_fast:.4f} (finite-speed loss)")
assert eta_fast < etas[0][1]
ax[0].set_xlabel("piston position L"); ax[0].set_ylabel("pressure 2E/L")
ax[0].set_title("Particle Carnot cycle (N = 20000)"); ax[0].legend()
Tcs = np.linspace(0.2, 0.95, 50)
ax[1].plot(Tcs, 1 - Tcs, "k-", label="$1-T_c/T_h$")
ax[1].plot(*zip(*etas), "ro", label="measured")
ax[1].set_xlabel("$T_c/T_h$"); ax[1].set_ylabel("efficiency"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "carnot_particles.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 2: worked-example numbers (notes §3.7, §4)

# %%
# Example 4.1: Carnot 500 K -> 300 K, 1 mol ideal gas, hot isotherm doubles the volume
Th, Tc, n = 500.0, 300.0, 1.0
Qh = n * R * Th * np.log(2); Qc = n * R * Tc * np.log(2); W = Qh - Qc
print(f"Ex 4.1: eta = {1 - Tc / Th:.3f}, Q_h = {Qh:.1f} J, Q_c = {Qc:.1f} J, W = {W:.1f} J, "
      f"dS_gas(hot leg) = {Qh / Th:.3f} J/K")
assert np.isclose(W / Qh, 0.4) and np.isclose(Qh / Th, Qc / Tc)
# Example 4.2: free expansion vs reversible isothermal, factor 2, 1 mol, 300 K
dS = R * np.log(2)
print(f"Ex 4.2: dS_gas = R ln2 = {dS:.3f} J/K; reversible: Q = {300 * dS:.1f} J, dS_res = {-dS:.3f}")
# Example 4.3: two bricks, C = 2 kg * 840 J/(kg K)
C, T1, T2 = 2 * 840.0, 400.0, 200.0
Tf = (T1 + T2) / 2
dS_tot = C * np.log(Tf / T1) + C * np.log(Tf / T2)
print(f"Ex 4.3: Tf = {Tf:.0f} K, dS_hot = {C * np.log(Tf / T1):.1f}, dS_cold = {C * np.log(Tf / T2):.1f}, "
      f"dS_total = {dS_tot:.1f} J/K (= C ln(9/8), ln(9/8) = {np.log(9 / 8):.4f})")
assert dS_tot > 0
# §3.5 melting and mixing; §3.7 Otto, Diesel, heat pump
print(f"Melting 1 kg ice: dS = {334e3 / 273.15:.1f} J/K; mixing 1 mol + 1 mol: dS = {2 * R * np.log(2):.2f} J/K")
g = 1.4
eta_otto = 1 - 10 ** (1 - g)
rD, rc = 18.0, 2.0
eta_diesel = 1 - rD ** (1 - g) * (rc**g - 1) / (g * (rc - 1))
print(f"Otto r=10: {eta_otto:.3f}; Diesel r=18, rc=2: {eta_diesel:.3f}; Otto r=18: {1 - 18**(1 - g):.3f}")
cop_hp = 293.15 / 20.0
print(f"Heat pump 0 C -> 20 C: Carnot COP = {cop_hp:.2f}")

# %% [markdown]
# ## Part 3: entropy production between two bodies (notes §3.5, eq. 3.9)
# S_i(U_i) = C_i ln U_i (constant heat capacity, k_B = 1). With U_A + U_B fixed,
# S_tot(U_A) is maximal where 1/T_A = dS_A/dU_A equals 1/T_B.

# %%
CA, CB, Utot = 3.0, 5.0, 10.0
UA = np.linspace(0.05, Utot - 0.05, 2001)
Stot = CA * np.log(UA) + CB * np.log(Utot - UA)
UA_star = UA[np.argmax(Stot)]
assert abs(UA_star - Utot * CA / (CA + CB)) < 0.01          # equal temperatures U_i/C_i
for _ in range(1000):                                        # random initial splits relax
    u0 = rng.uniform(0.1, Utot - 0.1)
    TA, TB = u0 / CA, (Utot - u0) / CB
    dU = 1e-4 * np.sign(TA - TB)                             # heat flows hot -> cold
    dS_small = CA * np.log((u0 - dU) / u0) + CB * np.log((Utot - u0 + dU) / (Utot - u0))
    dS_relax = (CA * np.log(Utot * CA / (CA + CB) / u0)
                + CB * np.log(Utot * CB / (CA + CB) / (Utot - u0)))
    assert dS_small >= -1e-15 and dS_relax >= -1e-12
fig, ax = plt.subplots(figsize=(5.5, 3.8))
ax.plot(UA, Stot); ax.axvline(UA_star, ls="--", c="k")
ax.set_xlabel("$U_A$"); ax.set_ylabel("$S_A+S_B$"); ax.set_title("Total entropy is maximal at $T_A=T_B$")
fig.tight_layout(); fig.savefig(FIG / "entropy_two_bodies.png", dpi=120); plt.close(fig)
print("Two bodies: S_tot max at U_A =", round(UA_star, 3), "; every hot->cold flow raises S_tot.")

# %% [markdown]
# ## Part 4: Maxwell relations and their consequences, symbolically (notes §3.9)

# %%
T, V, nn, a, b, cv, S0 = sp.symbols("T V n a b c_v S_0", positive=True)
for name, F in [("ideal", -nn * R * T * sp.log(V) + nn * cv * T * (1 - sp.log(T))),
                ("vdW", -nn * R * T * sp.log(V - nn * b) - a * nn**2 / V + nn * cv * T * (1 - sp.log(T)))]:
    P = -sp.diff(F, V); S = -sp.diff(F, T)
    assert sp.simplify(sp.diff(S, V) - sp.diff(P, T)) == 0          # (∂S/∂V)_T = (∂P/∂T)_V
    CV = sp.simplify(T * sp.diff(S, T))
    assert sp.simplify(sp.diff(CV, V)) == 0                          # C_V independent of V
    # C_P - C_V = T V α²/κ_T = -T (∂P/∂T)_V² / (∂P/∂V)_T   (eq. 3.27)
    diff_CP_CV = sp.simplify(-T * sp.diff(P, T)**2 / sp.diff(P, V))
    print(f"{name}: S = {sp.simplify(S)},  C_P - C_V = {diff_CP_CV}")
    if name == "ideal":
        assert sp.simplify(diff_CP_CV - nn * R) == 0
    else:
        assert sp.simplify(diff_CP_CV.subs({a: 0, b: 0}) - nn * R) == 0
        vdW_P, vdW_CPCV = P, diff_CP_CV
# Maxwell relation from U(S, V) of an ideal monatomic gas: (∂T/∂V)_S = -(∂P/∂S)_V
Ss, N_k = sp.symbols("S N_k", positive=True)
U_SV = V**sp.Rational(-2, 3) * sp.exp(2 * Ss / (3 * N_k))
assert sp.simplify(sp.diff(sp.diff(U_SV, Ss), V) - sp.diff(-sp.diff(U_SV, V), Ss) * -1) == 0
# Joule-Thomson: μ_JT = (V/C_P)(Tα - 1); zero for an ideal gas
Pp = sp.symbols("P", positive=True)
Vid = nn * R * T / Pp
assert sp.simplify(T * sp.diff(Vid, T) - Vid) == 0
# Euler: G = μ N for the ideal gas  G = N kT ln(P/P0) + N g(T)  is extensive in N at fixed T, P
print("SymPy: Maxwell relations, C_V(V)=const, C_P-C_V, ideal-gas μ_JT = 0 all verified.")

# %% [markdown]
# ## Part 5: the Legendre transform U(S) -> F(T) numerically (Example 4.4)
# Ideal monatomic gas, 1 mol, fixed V = V0: S(U) = S0 + (3/2) R ln(U/U0), so
# U(S) = U0 exp(2(S-S0)/3R). F(T) = min_S [U(S) - T S] should equal
# (3/2)RT - T[S0 + (3/2)R ln(T/T0)].

# %%
T0 = 300.0; U0 = 1.5 * R * T0; S0v = 0.0
S_grid = np.linspace(-60, 60, 200001)
U_tab = U0 * np.exp(2 * (S_grid - S0v) / (3 * R))
assert np.all(np.diff(U_tab, 2) > 0)                          # U(S) convex <=> C_V > 0
F_closed = lambda T_: 1.5 * R * T_ - T_ * (S0v + 1.5 * R * np.log(T_ / T0))
for T_ in (150.0, 300.0, 600.0):
    F_num = np.min(U_tab - T_ * S_grid)
    assert abs(F_num - F_closed(T_)) < 1e-3 * abs(F_closed(T_)) + 1e-2
    print(f"Ex 4.4: T = {T_:.0f} K: F numeric {F_num:.2f} J, closed form {F_closed(T_):.2f} J")
dT = 1e-3
S_from_F = -(F_closed(600 + dT) - F_closed(600 - dT)) / (2 * dT)
assert np.isclose(S_from_F, 1.5 * R * np.log(2), rtol=1e-6)    # -∂F/∂T = S
print(f"   S(600 K) = -dF/dT = {S_from_F:.3f} J/K = (3/2) R ln 2")

# %% [markdown]
# ## Part 6: C_P - C_V for liquid water at 25 C (Example 4.5)

# %%
Tw, vw, alpha, kappa = 298.15, 18.015e-3 / 997.05, 2.572e-4, 4.524e-10
dC_mol = Tw * vw * alpha**2 / kappa
print(f"Ex 4.5: molar volume {vw * 1e6:.3f} cm^3, C_P - C_V = {dC_mol:.3f} J/(mol K) "
      f"= {dC_mol / 18.015e-3:.1f} J/(kg K); C_P/C_V = {4181.3 / (4181.3 - dC_mol / 18.015e-3):.4f}")
assert 0.7 < dC_mol < 0.85

# %% [markdown]
# ## Part 7: Maxwell's equal-area construction for van der Waals CO2 (notes §3.11)

# %%
aC, bC = 0.3640, 4.267e-5                                     # Pa m^6/mol^2, m^3/mol
Tcr, Pcr = 8 * aC / (27 * R * bC), aC / (27 * bC**2)
P_vdw = lambda v, T_: R * T_ / (v - bC) - aC / v**2


def coexistence(T_):
    """Return (P*, v_l, v_g) with ∫_{v_l}^{v_g} (P - P*) dv = 0, eq. (3.40)."""
    sp_roots = np.roots([R * T_, -2 * aC, 4 * aC * bC, -2 * aC * bC**2])   # dP/dv = 0
    vs = np.sort(sp_roots[np.abs(sp_roots.imag) < 1e-12].real)
    vs = vs[vs > bC]
    Plo, Phi = max(P_vdw(vs[0], T_), 1.0), P_vdw(vs[1], T_)

    def outer_roots(Ps):
        r = np.roots([Ps, -(Ps * bC + R * T_), aC, -aC * bC])
        r = np.sort(r[np.abs(r.imag) < 1e-9].real)
        return r[0], r[-1]

    def area(Ps):
        vl, vg = outer_roots(Ps)
        return R * T_ * np.log((vg - bC) / (vl - bC)) + aC * (1 / vg - 1 / vl) - Ps * (vg - vl)

    Ps = brentq(area, Plo * (1 + 1e-9), Phi * (1 - 1e-9), xtol=1e-6)
    return (Ps, *outer_roots(Ps))


mu_vdw = lambda v, T_: -R * T_ * np.log(v - bC) - 2 * aC / v + R * T_ * v / (v - bC)  # μ minus f0(T)
print(f"CO2 vdW: Tc = {Tcr:.1f} K, Pc = {Pcr / 1e6:.2f} MPa")
Ts_coex = np.linspace(230, 303, 30)
coex = np.array([coexistence(T_) for T_ in Ts_coex])
for T_, (Ps, vl, vg) in zip(Ts_coex, coex):
    assert abs(mu_vdw(vl, T_) - mu_vdw(vg, T_)) < 1e-6 * R * T_     # equal μ <=> equal area
P280 = coexistence(280.0)
print(f"vdW CO2 coexistence at 280 K: P* = {P280[0] / 1e6:.2f} MPa, v_l = {P280[1] * 1e6:.1f}, "
      f"v_g = {P280[2] * 1e6:.1f} cm^3/mol")
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
v = np.linspace(6e-5, 1.2e-3, 2000)
for T_ in (260.0, 280.0, 295.0):
    Ps, vl, vg = coexistence(T_)
    ln, = ax[0].plot(v * 1e6, P_vdw(v, T_) / 1e6, label=f"{T_:.0f} K")
    ax[0].hlines(Ps / 1e6, vl * 1e6, vg * 1e6, colors=ln.get_color(), linestyles="--")
ax[0].set_ylim(0, 12); ax[0].set_xlabel("v (cm³/mol)"); ax[0].set_ylabel("P (MPa)")
ax[0].set_title("vdW CO₂ isotherms and equal-area tie lines"); ax[0].legend()
ax[1].plot(coex[:, 1] * 1e6, Ts_coex, "b-", coex[:, 2] * 1e6, Ts_coex, "r-")
ax[1].set_xscale("log"); ax[1].set_xlabel("v (cm³/mol)"); ax[1].set_ylabel("T (K)")
ax[1].set_title("Coexistence curve (liquid blue, gas red)")
fig.tight_layout(); fig.savefig(FIG / "vdw_equal_area.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 8: Clausius-Clapeyron for water (Example 4.6)
# Ideal vapour, liquid volume neglected: d ln P/dT = L(T)/(R T²), with
# L(T) = L0 + ΔC_P (T - T0). Compared with IAPWS saturation pressures.

# %%
L0, T0w, P0w, dCp = 40660.0, 373.15, 101325.0, -42.0
P_const = lambda T_: P0w * np.exp(-L0 / R * (1 / T_ - 1 / T0w))
sol = solve_ivp(lambda T_, y: [(L0 + dCp * (T_ - T0w)) / (R * T_**2)], (T0w, 273.16), [np.log(P0w)],
                dense_output=True, rtol=1e-10, atol=1e-12)
P_var = lambda T_: np.exp(sol.sol(T_)[0])
ref = {273.16: 611.66, 298.15: 3169.9, 323.15: 12352.0, 348.15: 38595.0}
for T_, Pref in ref.items():
    ec, ev = P_const(T_) / Pref - 1, P_var(T_) / Pref - 1
    print(f"T = {T_:.2f} K: IAPWS {Pref:8.1f} Pa, constant L {P_const(T_):8.1f} ({100 * ec:+.1f}%), "
          f"L(T) {P_var(T_):8.1f} ({100 * ev:+.1f}%)")
    assert abs(ev) < abs(ec) and abs(ev) < 0.03
P3000 = 101325 * (1 - 2.25577e-5 * 3000) ** 5.25588          # standard atmosphere
Tb_const = 1 / (1 / T0w - R / L0 * np.log(P3000 / P0w))
Tb_var = brentq(lambda T_: P_var(T_) - P3000, 340, 373.15)
print(f"Ex 4.6: P(3000 m) = {P3000 / 1e3:.2f} kPa; boiling point {Tb_const - 273.15:.1f} C (constant L), "
      f"{Tb_var - 273.15:.1f} C (L(T))")
assert abs(Tb_const - 273.15 - 89.8) < 0.1
Ts = np.linspace(273.16, 373.15, 200)
fig, ax = plt.subplots(figsize=(5.5, 3.8))
ax.semilogy(Ts - 273.15, P_const(Ts), label="constant L"); ax.semilogy(Ts - 273.15, P_var(Ts), label="L(T)")
ax.semilogy([t - 273.15 for t in ref], list(ref.values()), "ko", label="IAPWS")
ax.set_xlabel("T (°C)"); ax.set_ylabel("P_sat (Pa)"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "clausius_clapeyron_water.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 9: Joule-Thomson inversion curve of a van der Waals gas (notes §3.9)
# Inversion: T(∂V/∂T)_P = V, which for vdW gives RT = (2a/b)(1 - b/v)², P = RT/(v-b) - a/v².

# %%
vv = np.linspace(1.0001, 40, 4000) * 3 * bC                  # in units of v_c = 3b
T_inv = (2 * aC / (R * bC)) * (1 - bC / vv) ** 2
P_inv = R * T_inv / (vv - bC) - aC / vv**2
ok = P_inv > 0
Tr, Pr = T_inv[ok] / Tcr, P_inv[ok] / Pcr
assert np.allclose(Pr, 24 * np.sqrt(3 * Tr) - 12 * Tr - 27, atol=1e-8)
assert abs(Pr.max() - 9) < 1e-3 and abs(Tr[np.argmax(Pr)] - 3) < 0.01
# direct check of the sign of μ_JT ∝ T(∂V/∂T)_P - V on either side of the curve at T_r = 3
def mu_sign(T_, P_):
    vroot = brentq(lambda v_: P_vdw(v_, T_) - P_, 1.01 * bC, 1.0)
    dPdT, dPdv = R / (vroot - bC), -R * T_ / (vroot - bC)**2 + 2 * aC / vroot**3
    return np.sign(T_ * (-dPdT / dPdv) - vroot)
assert mu_sign(3 * Tcr, 4 * Pcr) > 0 and mu_sign(3 * Tcr, 12 * Pcr) < 0
fig, ax = plt.subplots(figsize=(5.5, 3.8))
ax.plot(Tr, Pr); ax.set_xlabel("$T/T_c$"); ax.set_ylabel("$P/P_c$")
ax.set_title("vdW Joule–Thomson inversion curve (cools inside)")
fig.tight_layout(); fig.savefig(FIG / "jt_inversion.png", dpi=120); plt.close(fig)
print(f"JT inversion (vdW): T_r from {Tr.min():.3f} to {Tr.max():.3f}, max P_r = {Pr.max():.3f} at T_r = 3")

# %% [markdown]
# ## Part 10: numbers used in the exercise solutions

# %%
print(f"A1: fridge Carnot COP = {275 / 20:.2f}; power for 300 W at COP 3 = {300 / 3:.0f} W")
print(f"A2: dS water 20->80 C = {4186 * np.log(353.15 / 293.15):.1f} J/K")
print(f"A3: dS = {2 * R * (1.5 * np.log(2) + np.log(3)):.2f} J/K")
print(f"A5: Otto r=8: {1 - 8**(1 - g):.3f}")
print(f"A6: ice {0.05 * 334e3 / 273.15:.2f}, room {-0.05 * 334e3 / 293.15:.2f}, "
      f"total {0.05 * 334e3 * (1 / 273.15 - 1 / 293.15):.2f} J/K")
print(f"B2: dS_gas {-R * np.log(2):.3f}, dS_res {R:.3f}, total {R * (1 - np.log(2)):.3f} J/K")
cvC = 2.5 * R; V1, V2, Ti = 1e-3, 2e-3, 300.0
Tf4 = Ti - aC / cvC * (1 / V1 - 1 / V2)
dS4 = cvC * np.log(Tf4 / Ti) + R * np.log((V2 - bC) / (V1 - bC))
print(f"B4: vdW CO2 free expansion: Tf = {Tf4:.2f} K, dS = {dS4:.3f} J/K (ideal: {R * np.log(2):.3f})")
v5 = 1e-3
cpcv5 = R / (1 - 2 * aC * (v5 - bC)**2 / (R * 300 * v5**3))
assert np.isclose(cpcv5, float(vdW_CPCV.subs({T: 300, V: v5, nn: 1, a: aC, b: bC})))
print(f"B5: vdW CO2 at 300 K, 1 L/mol: C_P - C_V = {cpcv5:.2f} J/(mol K)")
for gas, aa, bb in (("N2", 0.1370, 3.87e-5), ("H2", 0.02476, 2.661e-5)):
    print(f"B6: {gas}: low-density inversion T = 2a/Rb = {2 * aa / (R * bb):.0f} K")
Tpc = 1 / (1 / T0w - R / L0 * np.log(2))
dPdT_ice = 333.6e3 / (273.15 * (1.0002e-3 - 1.0908e-3))
print(f"B7: boiling at 2 atm {Tpc - 273.15:.1f} C; ice line slope {dPdT_ice / 1e5:.1f} bar/K")
print(f"C1: W_max bricks = {C * (T1 + T2 - 2 * np.sqrt(T1 * T2)):.0f} J, T_final = {np.sqrt(T1 * T2):.1f} K")
print(f"C3: plant 0.40 x heat pump COP 3 -> {0.40 * 3:.2f} J heat per J fuel")

# %% [markdown]
# ## Your turn
# 1. **Stirling cycle.** Replace the two adiabats in `carnot_engine` with isochores
#    (piston fixed, wall thermal, temperature stepped through Th -> Tc). Measure the
#    efficiency without a regenerator and compare with 1 - Tc/Th. Then count the
#    isochoric heat as perfectly regenerated and show the efficiency returns to Carnot.
# 2. **Gibbs paradox numerically.** Two halves, N particles each, partition removed.
#    Compute ΔS = 2Nk ln 2 for distinguishable gases from S(U,V,N) of eq. (3.17) with the
#    N inside the log treated two ways (ln V vs ln(V/N)); show only ln(V/N) gives 0 for
#    identical gases. SM-3 explains why ln(V/N) is right.
# 3. **Heat pump COP vs outdoor temperature.** Plot T_in/(T_in - T_out) for T_in = 293 K,
#    T_out from -20 C to 15 C, and a "real" COP = 0.45 x Carnot. At what outdoor temperature
#    does the real COP fall below 1/0.40 (break-even with a gas furnace via a 40% power plant)?
# 4. **Exercise B9** (reference solution below): latent heat of vdW CO2 two ways.

# %% [markdown]
# ## Reference solution for exercise B9 (write your own first)

# %%
Tb9 = np.linspace(240, 300, 13)
Pb9 = np.array([coexistence(T_) for T_ in Tb9])
for i in range(1, len(Tb9) - 1):
    dPdT = (Pb9[i + 1, 0] - Pb9[i - 1, 0]) / (Tb9[i + 1] - Tb9[i - 1])
    L_cc = Tb9[i] * (Pb9[i, 2] - Pb9[i, 1]) * dPdT            # Clausius-Clapeyron
    L_s = R * Tb9[i] * np.log((Pb9[i, 2] - bC) / (Pb9[i, 1] - bC))  # T Δs, from (3.28)
    assert abs(L_cc / L_s - 1) < 2e-3
print(f"B9: vdW CO2 latent heat at {Tb9[6]:.0f} K: T Δs = "
      f"{R * Tb9[6] * np.log((Pb9[6, 2] - bC) / (Pb9[6, 1] - bC)) / 1e3:.2f} kJ/mol, matches Clausius-Clapeyron")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
