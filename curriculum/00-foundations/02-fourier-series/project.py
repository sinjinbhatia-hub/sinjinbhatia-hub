# %% [markdown]
# # F2 Project: Fourier Series as Orthogonal Projection
#
# This script builds the Fourier machinery of F2 from the projection formula (3.5)
# and checks every number in the notes. Parts: (1) `fourier_coeffs` by the
# composite trapezoid rule, checked against SymPy's exact integrals for the four
# standard waves, plus the periodic-trapezoid accuracy; (2) partial sums and the
# Gibbs overshoot measured against Si(pi)/pi - 1/2; (3) Parseval and the Basel
# sum; (4) decay rates of |c_n| for three regularity classes; (5) the Dirichlet
# kernel and the convolution formula (3.14); (6) the heat equation on a rod,
# truncated series vs explicit finite differences; (7) the plucked string's
# normal modes; (8) d^2/dx^2 is self-adjoint: the periodic second-difference
# matrix is diagonal and real in the DFT basis (F3's circulant preview); (9)
# every number in solutions.md, exercise by exercise.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sympy as sp
from scipy.special import sici, iv
from scipy.integrate import quad

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
PI = np.pi

# %% [markdown]
# ## Part 0: The test functions and their exact coefficients
#
# The four standard waves on [-pi, pi] (Examples 4.1-4.4) and the C^1 cubic of
# §3.7. `np.sign(0) = 0` is deliberate: at a jump the trapezoid rule needs the
# midpoint value, which is also what Dirichlet's theorem says the series gives.
# Exact complex coefficients c_n come from (3.7) applied to the closed forms.


# %%
def square(x):
    """Square wave sgn(x): -1 on (-pi, 0), +1 on (0, pi). Example 4.1."""
    return np.sign(x)


def sawtooth(x):
    """Sawtooth f(x) = x on (-pi, pi). Example 4.2."""
    return np.asarray(x, float)


def triangle(x):
    """Triangle wave |x|. Example 4.3."""
    return np.abs(x)


def parabola(x):
    """f(x) = x^2. Example 4.4."""
    return np.asarray(x, float) ** 2


def cubic(x):
    """f(x) = x(pi^2 - x^2): C^1 as a periodic function, f'' jumps at +-pi. §3.7 table."""
    x = np.asarray(x, float)
    return x * (PI**2 - x**2)


def exact_c(name, n):
    """Exact c_n (complex form) for the named wave, from (3.7): c_n = (a_n - i b_n)/2, c_0 = a_0/2."""
    n = np.asarray(n, int)
    c = np.zeros(n.shape, complex)
    nz = n != 0
    m = n[nz].astype(float)
    odd = (n[nz] % 2 != 0)
    if name == "square":            # b_n = 4/(pi n) for odd n          -> c_n = -2i/(pi n)
        c[nz] = np.where(odd, -2j / (PI * m), 0)
    elif name == "sawtooth":        # b_n = 2(-1)^{n+1}/n              -> c_n = i(-1)^n / n
        c[nz] = 1j * (-1.0) ** m / m
    elif name == "triangle":        # a_0/2 = pi/2, a_n = -4/(pi n^2) odd -> c_n = -2/(pi n^2)
        c[~nz] = PI / 2
        c[nz] = np.where(odd, -2 / (PI * m**2), 0)
    elif name == "parabola":        # a_0/2 = pi^2/3, a_n = 4(-1)^n/n^2  -> c_n = 2(-1)^n/n^2
        c[~nz] = PI**2 / 3
        c[nz] = 2 * (-1.0) ** m / m**2
    elif name == "cubic":           # b_n = 12(-1)^{n+1}/n^3           -> c_n = 6i(-1)^n/n^3
        c[nz] = 6j * (-1.0) ** m / m**3
    else:
        raise ValueError(name)
    return c


WAVES = {"square": square, "sawtooth": sawtooth, "triangle": triangle, "parabola": parabola}

# %% [markdown]
# ## Part 1: Fourier coefficients by the trapezoid rule, checked against SymPy
#
# (3.9): c_n = (1/T) ∫ f(t) e^{-2 pi i n t/T} dt over one period. The composite
# trapezoid rule on M subintervals gives the two endpoints half weight; since f
# is periodic they are the same point, so this is the equal-weight sum (5.1),
# i.e. the DFT. For a periodic integrand the trapezoid rule is *spectrally*
# accurate: its error is the aliasing sum Σ_{k≠0} c_{n+kM}. For a smooth f that
# is exponentially small; for a jump it is O(1/M^2) (with midpoint values at the
# jumps, which the trapezoid rule supplies automatically at the seam ±pi).
#
# **Build it yourself first:** write `fourier_coeffs` and `partial_sum` before
# reading them.


# %%
def fourier_coeffs(f, N, T=2 * PI, M=2**17):
    """Complex Fourier coefficients c_n, |n| <= N, of a T-periodic f, eq. (3.9),
    by the composite trapezoid rule with M subintervals on [-T/2, T/2].
    Returns (n, c) with n = -N..N. f must be vectorized."""
    t = np.linspace(-T / 2, T / 2, M + 1)
    w = np.full(M + 1, T / M)
    w[0] = w[-1] = T / (2 * M)              # trapezoid end weights; they merge at the periodic seam
    ft = f(t) * w
    n = np.arange(-N, N + 1)
    c = np.empty(len(n), complex)
    for k, nk in enumerate(n):              # loop keeps memory small for large M
        c[k] = np.sum(ft * np.exp(-2j * PI * nk * t / T)) / T
    return n, c


def real_coeffs(n, c):
    """Convert c_n (n = -N..N) to real coefficients via (3.7) inverted:
    a_0 = 2 c_0, a_n = c_n + c_{-n}, b_n = i (c_n - c_{-n}). Returns (a, b) indexed 0..N (b_0 = 0)."""
    N = n.max()
    idx = {int(m): k for k, m in enumerate(n)}
    a = np.zeros(N + 1); b = np.zeros(N + 1)
    a[0] = 2 * c[idx[0]].real
    for m in range(1, N + 1):
        a[m] = (c[idx[m]] + c[idx[-m]]).real
        b[m] = (1j * (c[idx[m]] - c[idx[-m]])).real
    return a, b


def partial_sum(n, c, x):
    """S_N f(x) = Σ_{|n|<=N} c_n e^{inx}, eq. (3.6) in complex form (real part returned)."""
    return np.real(np.exp(1j * np.outer(x, n)) @ c)


# --- SymPy: exact a_n, b_n by (3.4) for n = 0..6, splitting at the breaks via sign/Abs.
xs = sp.symbols("x", real=True)
SYM = {"square": sp.sign(xs), "sawtooth": xs, "triangle": sp.Abs(xs), "parabola": xs**2}
NSYM = 6
for name, f in WAVES.items():
    n, c = fourier_coeffs(f, NSYM)
    a_num, b_num = real_coeffs(n, c)
    for m in range(NSYM + 1):
        a_ex = sp.integrate(SYM[name] * sp.cos(m * xs), (xs, -sp.pi, sp.pi)) / sp.pi
        b_ex = sp.integrate(SYM[name] * sp.sin(m * xs), (xs, -sp.pi, sp.pi)) / sp.pi if m else 0
        assert abs(a_num[m] - float(a_ex)) < 1e-8, (name, m, a_num[m], a_ex)
        assert abs(b_num[m] - float(b_ex)) < 1e-8, (name, m, b_num[m], b_ex)
    # and against the closed forms of §4
    assert np.allclose(c, exact_c(name, n), atol=1e-8), name
    print(f"{name:9s} trapezoid == SymPy == closed form to 1e-8 for n <= {NSYM}")

# Closed forms with a symbolic integer n, so you can see where (-1)^n comes from.
nn = sp.symbols("n", integer=True, positive=True)
print("  sawtooth b_n =", sp.simplify(sp.integrate(xs * sp.sin(nn * xs), (xs, -sp.pi, sp.pi)) / sp.pi))
print("  parabola a_n =", sp.simplify(sp.integrate(xs**2 * sp.cos(nn * xs), (xs, -sp.pi, sp.pi)) / sp.pi))
print("  triangle a_n =", sp.simplify(2 * sp.integrate(xs * sp.cos(nn * xs), (xs, 0, sp.pi)) / sp.pi))
print("  square   b_n =", sp.simplify(2 * sp.integrate(sp.sin(nn * xs), (xs, 0, sp.pi)) / sp.pi))

# Example 4.1-4.4 headline numbers.
assert np.isclose(4 / PI, 1.2732, atol=5e-5) and np.isclose(4 / (3 * PI), 0.4244, atol=5e-5)
assert np.isclose(-4 / (9 * PI), -0.1415, atol=5e-5) and np.isclose(-4 / (25 * PI), -0.0509, atol=5e-5)

# --- Periodic-trapezoid accuracy. Smooth periodic f: c_n(e^{cos x}) = I_n(1) (modified Bessel).
for M in (4, 8, 16, 32):
    n, c = fourier_coeffs(lambda t: np.exp(np.cos(t)), 2, M=M)
    err = abs(c[n == 1][0] - iv(1, 1))
    print(f"  e^cos x, M={M:3d}: |c_1 - I_1(1)| = {err:.1e}")
    if M >= 16:
        assert err < 1e-12, "trapezoid should be spectrally accurate for smooth periodic f"
# Square wave: error decays like 1/M^2.
Ms = 2 ** np.arange(6, 14)
errs = []
for M in Ms:
    n, c = fourier_coeffs(square, 5, M=M)
    errs.append(np.abs(c - exact_c("square", n)).max())
slope_trap = np.polyfit(np.log(Ms), np.log(errs), 1)[0]
print(f"  square wave trapezoid error ~ M^{slope_trap:.2f}  (expect -2)")
assert abs(slope_trap + 2) < 0.1
print("Part 1 passed.")

# %% [markdown]
# ## Part 2: Partial sums and the Gibbs phenomenon
#
# Theorem 3.15: for the square wave the first peak of S_{2M-1} is at pi/(2M) and
# tends to (2/pi) Si(pi) = 1.17898; the overshoot tends to Si(pi)/pi - 1/2 = 0.0895
# of the jump, and is slightly larger than that at every finite N. For the sawtooth the peak sits where f itself is
# pi - t, so measure S_N f(x) - f(x); the limit is 2 Si(pi) - pi = 0.5623 =
# 0.0895 * 2pi.

# %%
Si_pi = sici(PI)[0]
Si_quad = quad(lambda u: np.sin(u) / u, 0, PI)[0]
assert np.isclose(Si_pi, Si_quad) and np.isclose(Si_pi, 1.8519370, atol=1e-7)
gibbs_frac = Si_pi / PI - 0.5
print(f"Si(pi) = {Si_pi:.7f};  (2/pi) Si(pi) = {2 / PI * Si_pi:.7f};  overshoot/jump = {gibbs_frac:.7f}")
assert np.isclose(gibbs_frac, 0.08949, atol=1e-5)

NS = [1, 3, 5, 11, 51]
xg = np.linspace(-PI, PI, 20001)
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
over = {"square": {}, "sawtooth": {}}
for ax, name in zip(axes, ("square", "sawtooth")):
    f = WAVES[name]
    ax.plot(xg, f(xg), "k", lw=1, label="f")
    for k, N in enumerate(NS):
        n = np.arange(-N, N + 1)
        S = partial_sum(n, exact_c(name, n), xg)
        ax.plot(xg, S, color=f"C{k}", lw=1, label=f"N={N}")
        pos = (xg > 0) & (xg < PI)
        over[name][N] = (S - f(xg))[pos].max()      # overshoot above f, to the right of the jump at 0 / left of pi
    ax.set_title(f"{name}: partial sums S_N f"); ax.grid(alpha=.3); ax.legend(fontsize=8)
    ax.set_xlabel("x")
fig.tight_layout(); fig.savefig(FIG / "partial_sums_gibbs.png", dpi=120); plt.close(fig)

limit_sq, limit_saw = 2 / PI * Si_pi - 1, 2 * Si_pi - PI
for N in NS:
    print(f"  N={N:2d}: square overshoot {over['square'][N]:.4f} (limit {limit_sq:.4f});  "
          f"sawtooth overshoot {over['sawtooth'][N]:.4f} (limit {limit_saw:.4f} = 0.0895*2pi)")
assert abs(over["square"][51] - limit_sq) < 2e-3 and abs(over["sawtooth"][51] - limit_saw) < 2e-3
assert abs(over["square"][11] - 0.1813) < 1e-3 and abs(over["square"][51] - 0.1791) < 1e-3   # Example 4.8
# Finite-M formula from step 4 of the derivation: S_N f(x_M) = (2/pi) ∫_0^pi sin u / (2M sin(u/2M)) du.
for N in (11, 51):
    M = (N + 1) // 2
    val = 2 / PI * quad(lambda u: np.sin(u) / (2 * M * np.sin(u / (2 * M))), 0, PI)[0]
    assert abs(val - 1 - over["square"][N]) < 1e-4, (N, val, over["square"][N])
# The overshoot is above the limit at finite M and decreases toward it: Gibbs does not go away.
assert over["square"][11] > over["square"][51] > limit_sq
print("Part 2 passed: Gibbs overshoot ~ 9% of the jump and does not shrink with N.")

# %% [markdown]
# ## Part 3: Parseval and the Basel problem
#
# (3.11): (1/pi) ∫|f|^2 = a_0^2/2 + Σ (a_n^2 + b_n^2) = 2 Σ_{n in Z} |c_n|^2.
# Example 4.5: sawtooth gives 4 Σ 1/n^2 = 2 pi^2/3, parabola gives Σ 1/n^4 = pi^4/90.

# %%
norm2 = {"square": 2.0, "sawtooth": 2 * PI**2 / 3, "triangle": 2 * PI**2 / 3, "parabola": 2 * PI**4 / 5}
Nbig = 10_000
n = np.arange(-Nbig, Nbig + 1)
for name in WAVES:
    lhs = norm2[name]                       # (1/pi) ∫ |f|^2 dx, by hand
    rhs = 2 * np.sum(np.abs(exact_c(name, n)) ** 2)
    print(f"  Parseval {name:9s}: (1/pi)||f||^2 = {lhs:.6f}, 2 Σ|c_n|^2 (|n|<={Nbig}) = {rhs:.6f}")
    assert abs(lhs - rhs) < 1e-3            # tail ~ 1/N for the discontinuous waves
    # the trapezoid-rule norm agrees too (to O(h): the single sample at a jump carries the midpoint value)
    Mt = 2**14; t = np.linspace(-PI, PI, Mt + 1); w = np.full(Mt + 1, 2 * PI / Mt); w[0] = w[-1] = PI / Mt
    assert abs(np.sum(WAVES[name](t) ** 2 * w) / PI - lhs) < 1e-3

k = np.arange(1, 1001)
basel = np.cumsum(1 / k**2)
print(f"  Basel: Σ_{{n<=10}} = {basel[9]:.6f}, Σ_{{n<=100}} = {basel[99]:.6f}, Σ_{{n<=1000}} = {basel[999]:.6f}, pi^2/6 = {PI**2 / 6:.6f}")
assert np.isclose(basel[9], 1.549768, atol=1e-6) and np.isclose(basel[999], PI**2 / 6, atol=1.1e-3)
# The Parseval arithmetic of Example 4.5, literally.
assert np.isclose((2 * PI**2 / 3) / 4, PI**2 / 6)
assert np.isclose((2 * PI**4 / 5 - 2 * PI**4 / 9) / 16, PI**4 / 90)
assert np.isclose(np.sum(1 / k**4), PI**4 / 90, atol=1e-9)
fig, ax = plt.subplots(figsize=(5.5, 3.5))
ax.semilogx(k, basel, label="partial sums of Σ 1/n²"); ax.axhline(PI**2 / 6, color="k", ls="--", lw=1, label="π²/6")
ax.set_xlabel("N"); ax.grid(alpha=.3); ax.legend(); ax.set_title("Basel via Parseval of the sawtooth")
fig.tight_layout(); fig.savefig(FIG / "basel.png", dpi=120); plt.close(fig)
print("Part 3 passed.")

# %% [markdown]
# ## Part 4: Decay rates (Theorems 3.9 and 3.10)
#
# A jump gives |c_n| ~ 1/n (square wave); a continuous function with a corner
# gives 1/n^2 (triangle); C^1 with a jump in f'' gives 1/n^3 (cubic x(pi^2-x^2)).
# Coefficients here are *numerical* (trapezoid), to show the slopes come out of
# the computation, not the formula.

# %%
Ndec = 61
fig, ax = plt.subplots(figsize=(6, 4.5))
slopes = {}
for k, (name, f, expect) in enumerate([("square", square, -1), ("triangle", triangle, -2), ("cubic", cubic, -3)]):
    n, c = fourier_coeffs(f, Ndec, M=2**16)
    mask = (n >= 1) & (np.abs(c) > 1e-12)        # skip the exact zeros (even n of square/triangle)
    nn_, cc_ = n[mask], np.abs(c[mask])
    slope = np.polyfit(np.log(nn_), np.log(cc_), 1)[0]
    slopes[name] = slope
    ax.loglog(nn_, cc_, "o", ms=4, color=f"C{k}", label=f"{name}: fitted slope {slope:.2f}")
    ax.loglog(nn_, cc_[0] * (nn_ / nn_[0]) ** expect, "-", lw=1, color=f"C{k}", alpha=.6)
    assert abs(slope - expect) < 0.05, (name, slope)
    if name == "cubic":                            # check the closed form 6/n^3 from §3.7
        assert np.allclose(cc_, 6 / nn_**3, rtol=1e-6)
ax.set_xlabel("n"); ax.set_ylabel("|c_n|"); ax.set_title("Smoothness ⇒ decay: |c_n| vs n (lines: n⁻¹, n⁻², n⁻³)")
ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "decay_rates.png", dpi=120); plt.close(fig)
print("Part 4 passed: slopes", {k: round(v, 3) for k, v in slopes.items()})

# %% [markdown]
# ## Part 5: The Dirichlet kernel and S_N f as a convolution
#
# (3.13): D_N(x) = Σ_{|n|<=N} e^{inx} = sin((N+1/2)x)/sin(x/2), with
# (1/2pi) ∫ D_N = 1 and D_N(0) = 2N+1. (3.14): S_N f(x) = (1/2pi) ∫ f(y) D_N(x-y) dy.

# %%
def dirichlet_kernel(N, x):
    """D_N(x) by the closed form (3.13), with the removable singularity at multiples of 2pi filled in."""
    x = np.asarray(x, float)
    s = np.sin(x / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        D = np.where(np.abs(s) > 1e-12, np.sin((N + 0.5) * x) / s, 2 * N + 1.0)
    return D


xk = np.linspace(-PI, PI, 4001)
fig, ax = plt.subplots(figsize=(7, 4))
for k, N in enumerate((2, 5, 10, 20)):
    D_closed = dirichlet_kernel(N, xk)
    D_sum = np.real(np.exp(1j * np.outer(xk, np.arange(-N, N + 1))).sum(axis=1))
    assert np.allclose(D_closed, D_sum, atol=1e-9), "geometric-sum derivation of (3.13)"
    assert np.isclose(np.trapezoid(D_closed, xk) / (2 * PI), 1.0, atol=1e-6)   # (1/2pi) ∫ D_N = 1
    assert np.isclose(dirichlet_kernel(N, 0.0), 2 * N + 1)
    ax.plot(xk, D_closed, color=f"C{k}", lw=1, label=f"N={N}")
ax.set_title("Dirichlet kernel D_N(x) = sin((N+½)x)/sin(x/2)"); ax.set_xlabel("x"); ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "dirichlet_kernel.png", dpi=120); plt.close(fig)

# Convolution formula (3.14) vs the direct partial sum, for the parabola at N = 8.
N = 8
y = np.linspace(-PI, PI, 2**12 + 1); wy = np.full(y.size, 2 * PI / (y.size - 1)); wy[0] = wy[-1] = wy[0] / 2
for x0 in (0.0, 1.0, 2.5):
    conv = np.sum(parabola(y) * dirichlet_kernel(N, x0 - y) * wy) / (2 * PI)
    direct = partial_sum(np.arange(-N, N + 1), exact_c("parabola", np.arange(-N, N + 1)), np.array([x0]))[0]
    assert abs(conv - direct) < 1e-6, (x0, conv, direct)
print("Part 5 passed: D_N closed form, unit mass, and S_N f = f * D_N / 2pi.")

# %% [markdown]
# ## Part 6: The heat equation on a rod (Example 4.6)
#
# u_t = k u_xx on [0, L], u(0,t) = u(L,t) = 0, u(x,0) = step (u_0 on the left
# half). Series (3.17) with b_n = (2u_0/(n pi))(1 - cos(n pi/2)), each mode
# decaying as exp(-k (n pi/L)^2 t). Compared with an explicit finite-difference
# (FTCS) solver with r = k dt/dx^2 = 0.4 < 1/2.
#
# **Build it yourself first:** the FTCS update is one line of NumPy.

# %%
def heat_series(x, t, L=1.0, k=1.0, nmax=400):
    """u(x,t) from (3.17) for the half-step initial condition, truncated at nmax modes."""
    n = np.arange(1, nmax + 1)
    b = 2 / (n * PI) * (1 - np.cos(n * PI / 2))
    return (b * np.exp(-(n * PI / L) ** 2 * k * t)) @ np.sin(np.outer(n, PI * x / L))


def heat_fd(x, times, L=1.0, k=1.0, r=0.4):
    """Explicit FTCS: u_j^{m+1} = u_j^m + r (u_{j+1} - 2u_j + u_{j-1}), r = k dt/dx^2.
    Returns the solution at each requested time (times must be increasing)."""
    dx = x[1] - x[0]; dt = r * dx**2 / k
    u = np.where(x < L / 2, 1.0, 0.0); u[0] = u[-1] = 0.0; u[np.isclose(x, L / 2)] = 0.5
    t, out = 0.0, []
    for T in times:
        while t < T - 1e-12:
            u[1:-1] += r * (u[2:] - 2 * u[1:-1] + u[:-2]); t += dt
        out.append(u.copy())
    return out


# Example 4.6 coefficients and time scales.
nb = np.arange(1, 7)
b_step = 2 / (nb * PI) * (1 - np.cos(nb * PI / 2))
print("  step b_n/u_0, n=1..6:", np.round(b_step, 4))
assert np.allclose(b_step, [0.6366, 0.6366, 0.2122, 0.0, 0.1273, 0.2122], atol=5e-5)
L_cu, k_cu = 0.10, 1.1e-4                                  # 10 cm copper rod
tau = L_cu**2 / (nb**2 * PI**2 * k_cu)
print("  copper tau_n (s), n=1..3:", np.round(tau[:3], 2))
assert np.allclose(tau[:3], [9.21, 2.30, 1.02], atol=6e-3)
assert np.isclose(np.exp(-4), 0.0183, atol=5e-5) and np.isclose(np.exp(-9), 1.234e-4, atol=1e-6)

xr = np.linspace(0, 1, 201)
times = [0.001, 0.005, 0.02, 0.1]
fd = heat_fd(xr, times)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(xr, np.where(xr < 0.5, 1, 0), "k:", lw=1, label="t=0 (step)")
for k, (T, u_fd) in enumerate(zip(times, fd)):
    u_ser = heat_series(xr, T)
    err = np.abs(u_fd - u_ser).max()
    print(f"  t={T}: max |FD - series| = {err:.1e}")
    assert err < 2e-3, "finite differences and the truncated series must agree"
    ax.plot(xr, u_ser, color=f"C{k}", lw=1.5, label=f"series t={T}")
    ax.plot(xr[::8], u_fd[::8], "o", color=f"C{k}", ms=3, label="FD" if k == 0 else None)
ax.set_xlabel("x/L"); ax.set_ylabel("u/u₀"); ax.set_title("Heat equation: step initial data (L=k=1)")
ax.grid(alpha=.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "heat_rod.png", dpi=120); plt.close(fig)
# Smoothing: at t = tau_1 the mode-2 amplitude is down by e^{-4} relative to its start.
assert np.isclose(np.exp(-(2 * PI) ** 2 * (1 / PI**2)), np.exp(-4))
print("Part 6 passed.")

# %% [markdown]
# ## Part 7: The plucked string (Example 4.7)
#
# Triangle initial shape of height h at L/2, released from rest. A_n are the sine
# coefficients: A_n = 8h/(n^2 pi^2) (-1)^{(n-1)/2} for odd n, 0 for even n.
# Each mode oscillates at omega_n = n pi c/L (3.18).

# %%
L, h, c_wave = 1.0, 1.0, 1.0


def pluck(x):
    return np.where(x < L / 2, 2 * h * x / L, 2 * h * (L - x) / L)


def sine_coeffs(phi, nmax, L=1.0, M=2**14):
    """Half-range sine coefficients b_n = (2/L) ∫_0^L phi sin(n pi x/L) dx (§3.5), trapezoid rule."""
    x = np.linspace(0, L, M + 1); w = np.full(M + 1, L / M); w[0] = w[-1] = L / (2 * M)
    n = np.arange(1, nmax + 1)
    return n, (2 / L) * (np.sin(np.outer(n, PI * x / L)) @ (phi(x) * w))


nA, A_num = sine_coeffs(pluck, 12)
A_exact = np.where(nA % 2 == 1, 8 * h / (nA**2 * PI**2) * (-1.0) ** ((nA - 1) // 2), 0.0)
assert np.allclose(A_num, A_exact, atol=1e-8)
print("  A_n/h, n=1,3,5:", np.round(A_exact[[0, 2, 4]], 4), "  even n:", np.round(A_exact[[1, 3, 5]], 10))
assert np.allclose(A_exact[[0, 2, 4]], [0.8106, -0.0901, 0.0324], atol=5e-5)
# Midpoint check: Σ A_n sin(n pi/2) = h, i.e. (8/pi^2) Σ_odd 1/n^2 = 1.
nodd = np.arange(1, 20001, 2)
assert np.isclose(8 / PI**2 * np.sum(1 / nodd**2), 1.0, atol=1e-4)
# Energy fractions ∝ omega_n^2 A_n^2 ∝ 1/n^2 over odd n: 81%, 9.0%, 3.2%.
E = (nodd**2) * (8 / (nodd**2 * PI**2)) ** 2
frac = E / E.sum()
print("  energy fractions n=1,3,5:", np.round(frac[:3], 4))
assert np.allclose(frac[:3], [0.8106, 0.0901, 0.0324], atol=5e-4)

# Figure grid: (top) the first six normal modes; (bottom) the plucked string at several times.
xs_ = np.linspace(0, L, 400)
nmodes = np.arange(1, 201)
A_all = np.where(nmodes % 2 == 1, 8 * h / (nmodes**2 * PI**2) * (-1.0) ** ((nmodes - 1) // 2), 0.0)
T1 = 2 * L / c_wave                                        # period of the fundamental
fig, axes = plt.subplots(2, 6, figsize=(14, 5))
for j, n in enumerate(range(1, 7)):
    for ph in np.linspace(0, PI, 5):
        axes[0, j].plot(xs_, np.sin(n * PI * xs_ / L) * np.cos(ph), color="C0", alpha=.25 + .75 * (ph == 0))
    axes[0, j].set_title(f"mode n={n}, ω={n}π c/L"); axes[0, j].set_ylim(-1.1, 1.1); axes[0, j].set_xticks([])
for j, frac_t in enumerate([0, 1 / 8, 1 / 4, 3 / 8, 1 / 2, 3 / 4]):
    t = frac_t * T1
    u = (A_all * np.cos(nmodes * PI * c_wave * t / L)) @ np.sin(np.outer(nmodes, PI * xs_ / L))
    axes[1, j].plot(xs_, u, color="C1"); axes[1, j].set_ylim(-1.1, 1.1)
    axes[1, j].set_title(f"pluck, t = {frac_t:.3g} T₁"); axes[1, j].set_xlabel("x/L")
    if j == 0:
        assert np.allclose(u, pluck(xs_), atol=5e-3)        # 200 modes reproduce the triangle
fig.suptitle("Top: normal modes sin(nπx/L). Bottom: center pluck = Σ A_n sin(nπx/L) cos(ω_n t)")
fig.tight_layout(); fig.savefig(FIG / "string_modes.png", dpi=120); plt.close(fig)
print("Part 7 passed.")

# %% [markdown]
# ## Part 8: d²/dx² is self-adjoint, numerically (Theorem 5.1, F3 preview)
#
# Sample [-pi, pi) at N points. The periodic second-difference matrix D2 is
# symmetric (self-adjoint in R^N). In the DFT basis it is diagonal and real
# with entries -(4/h²) sin²(nh/2) ≈ -n²: this is the circulant diagonalization
# of F3. The *spectral* second derivative (FFT, multiply by -n², inverse FFT)
# has eigenvalues exactly -n². Also: the sampled cosines/sines of §2 are
# exactly orthogonal in R^N.

# %%
Ng = 64
hg = 2 * PI / Ng
xj = -PI + hg * np.arange(Ng)
D2 = (-2 * np.eye(Ng) + np.eye(Ng, k=1) + np.eye(Ng, k=-1) + np.eye(Ng, k=Ng - 1) + np.eye(Ng, k=-(Ng - 1))) / hg**2
assert np.allclose(D2, D2.T), "second-difference matrix with periodic BC is symmetric"
for _ in range(20):                                        # <f, D2 g> = <D2 f, g> on random vectors
    f_, g_ = rng.standard_normal(Ng) + 1j * rng.standard_normal(Ng), rng.standard_normal(Ng) + 1j * rng.standard_normal(Ng)
    assert np.isclose(np.vdot(f_, D2 @ g_), np.vdot(D2 @ f_, g_))

freqs = np.arange(-Ng // 2, Ng // 2)
F = np.exp(1j * np.outer(xj, freqs)) / np.sqrt(Ng)        # columns: sampled e^{inx}/sqrt(N), orthonormal in C^N
assert np.allclose(F.conj().T @ F, np.eye(Ng), atol=1e-12), "sampled exponentials are orthonormal"
Lam = F.conj().T @ D2 @ F
off = Lam - np.diag(np.diag(Lam))
assert np.abs(off).max() < 1e-10, "D2 is diagonal in the Fourier basis"
lam = np.diag(Lam)
assert np.abs(lam.imag).max() < 1e-10, "eigenvalues are real (self-adjoint)"
assert np.allclose(lam.real, -(4 / hg**2) * np.sin(freqs * hg / 2) ** 2)
small = np.abs(freqs) <= 3
print("  D2 eigenvalues in Fourier basis for n=-3..3:", np.round(lam.real[small], 3), " vs -n²:", -(freqs[small] ** 2))
assert np.allclose(lam.real[small], -(freqs[small] ** 2), rtol=0.01)
# Spectral second derivative: exact -n² on trigonometric polynomials.
f_trig = np.cos(3 * xj) + 0.5 * np.sin(5 * xj)
fpp_exact = -9 * np.cos(3 * xj) - 12.5 * np.sin(5 * xj)
kf = np.fft.fftfreq(Ng, d=1 / Ng)                          # integer frequencies
fpp_spec = np.real(np.fft.ifft(-(kf**2) * np.fft.fft(f_trig)))
assert np.allclose(fpp_spec, fpp_exact, atol=1e-10)
# The DFT computes the trapezoid-rule coefficients (5.1): compare with fourier_coeffs on the parabola.
c_dft = np.fft.fft(parabola(xj)) / Ng * np.exp(-1j * kf * PI)   # phase from starting the grid at -pi
n_, c_tr = fourier_coeffs(parabola, 5, M=Ng)
for m in range(-5, 6):
    assert np.isclose(c_dft[kf == m][0], c_tr[n_ == m][0], atol=1e-12)
# §2: sampled cosines are orthogonal vectors in R^N (exactly), with ||cos||² = N/2.
C = np.cos(np.outer(np.arange(1, Ng // 2), xj))
G = C @ C.T
assert np.allclose(G, (Ng / 2) * np.eye(Ng // 2 - 1), atol=1e-9)
print("Part 8 passed: D2 symmetric, diagonal & real in the DFT basis, ≈ -n² for small n.")

# %% [markdown]
# ## Your turn
#
# 1. **Full-wave rectified sine (EE/ME).** Compute the coefficients of |sin x|
#    with `fourier_coeffs`. Confirm a_0/2 = 2/pi (the DC level a rectifier
#    delivers), that all odd harmonics vanish, and that a_{2m} = -4/(pi(4m²-1)).
#    Why is the ripple at *twice* the line frequency?
# 2. **Non-centered interval.** Let g(x) = x on [0, 2pi) (periodic). Compute its
#    c_n with `fourier_coeffs(g, N)` after writing g for inputs in [-pi, pi)
#    (hint: `np.mod(x, 2*pi)`). Compare |c_n| and arg c_n with the centered
#    sawtooth: the magnitudes agree and the phases differ by n·pi, because
#    g(x) = pi + f(x - pi). This is the shift theorem of F3.
# 3. **How many terms for uniform error < 0.01?** For x² and for |x|, compute
#    max_x |f - S_N f| on a fine grid and find the smallest N that gets below
#    0.01. Compare with the tail bounds 4/N and 2/(pi N) from Theorem 3.16.
# 4. **Pluck at L/4 instead of L/2.** Which modes return? Which is now largest
#    after the fundamental?
# 5. **Fejér (C3).** Average the partial sums: sigma_N = (S_0 + ... + S_N)/(N+1).
#    Plot for the square wave at N = 51 and measure the overshoot. It should be
#    gone (the Fejér kernel is nonnegative).

# %%
# Reference for B8 (the project exercise), so you can compare after trying it yourself.
def quartic(x):
    """(pi² - x²)²: C² as a periodic function, f''' jumps at ±pi -> |c_n| ~ 1/n^4 (Exercise B4)."""
    x = np.asarray(x, float)
    return (PI**2 - x**2) ** 2


n, c = fourier_coeffs(quartic, 40, M=2**16)
mask = n >= 1
slope4 = np.polyfit(np.log(n[mask]), np.log(np.abs(c[mask])), 1)[0]
assert abs(slope4 + 4) < 0.05, slope4
assert np.allclose(c[mask], -24 * (-1.0) ** n[mask] / n[mask] ** 4, atol=1e-8)   # B4's closed form, c_n = a_n/2
assert np.isclose(c[n == 0][0].real, 8 * PI**4 / 15, atol=1e-8)
print(f"  B8: quartic slope {slope4:.3f} (expect -4); a_0/2 = 8pi^4/15 = {8 * PI**4 / 15:.4f}")


def terms_for_uniform_error(name, tol=0.01, Nmax=600):
    """Smallest N with max_x |f - S_N f| < tol. Builds every term a_n cos(nx) once (both waves are
    even, so b_n = 0) and accumulates with cumsum, so the cost is one (Nmax x grid) array, not Nmax partial sums."""
    xf = np.linspace(-PI, PI, 8001)
    fx = WAVES[name](xf)
    nn_ = np.arange(1, Nmax + 1)
    a_n = 2 * exact_c(name, nn_).real                      # a_n = c_n + c_{-n} = 2 Re c_n for real even f
    terms = a_n[:, None] * np.cos(np.outer(nn_, xf))       # row n-1: a_n cos(n x)
    S = exact_c(name, [0])[0].real + np.cumsum(terms, axis=0)   # row N-1: S_N f on the grid
    err = np.abs(S - fx).max(axis=1)
    hits = np.nonzero(err < tol)[0]
    return int(nn_[hits[0]]) if hits.size else None


N_par, N_tri = terms_for_uniform_error("parabola"), terms_for_uniform_error("triangle")
print(f"  B8: uniform error < 0.01 needs N = {N_par} for x² (bound 4/N -> 400) and N = {N_tri} for |x| (bound 2/(pi N) -> 64)")
assert N_par == 400 and N_tri == 63

# %% [markdown]
# ## Part 9: Exercise checks
#
# Every number quoted in `solutions.md` is recomputed and asserted here, in
# exercise order, so you can check your own work against the script rather than
# against the text. (B8's reference is the block just above.) Reuses
# `fourier_coeffs`, `real_coeffs`, `exact_c`, `partial_sum`, `dirichlet_kernel`,
# `sine_coeffs`, `cubic`, Part 2's grid `xg` and measured overshoots `over`.

# %%
def trap_weights(a, b, M):
    """Composite-trapezoid nodes and weights on [a, b] with M subintervals."""
    t_ = np.linspace(a, b, M + 1)
    w_ = np.full(M + 1, (b - a) / M); w_[0] = w_[-1] = (b - a) / (2 * M)
    return t_, w_


# --- A1: orthogonality integrals (Theorems 3.2, 3.3).
assert abs(quad(lambda x: np.cos(2 * x) * np.cos(5 * x), -PI, PI)[0]) < 1e-12
assert np.isclose(quad(lambda x: np.sin(3 * x) ** 2, -PI, PI)[0], PI)
assert abs(quad(lambda x: np.cos(4 * x) * np.sin(4 * x), -PI, PI)[0]) < 1e-12
assert abs(quad(np.cos, -PI, PI)[0]) < 1e-12 and abs(quad(np.sin, -PI, PI)[0]) < 1e-12      # (d): ∫ e^{ix} = 0

# --- A2: coefficients read off a trigonometric polynomial; norms by Pythagoras.
f_A2 = lambda x: 2 + 3 * np.cos(x) - np.sin(x) + 4 * np.sin(2 * x) - 5 * np.cos(3 * x)
n, c = fourier_coeffs(f_A2, 4)
a, b = real_coeffs(n, c)
assert np.allclose(a, [4, 3, 0, -5, 0], atol=1e-10) and np.allclose(b, [0, -1, 4, 0, 0], atol=1e-10)
assert np.allclose(c[(n >= -3) & (n <= 3)], [-2.5, 2j, (3 - 1j) / 2, 2, (3 + 1j) / 2, -2j, -2.5], atol=1e-10)
t9, w9 = trap_weights(-PI, PI, 2**12)
assert np.isclose(np.sum(f_A2(t9) ** 2 * w9), 59 * PI) and np.isclose(2 * PI * np.sum(np.abs(c) ** 2), 59 * PI)
r1 = f_A2(t9) - (2 + 3 * np.cos(t9) - np.sin(t9)); r2 = r1 - 4 * np.sin(2 * t9)              # f - S_1 f, f - S_2 f
assert np.isclose(np.sum(r1**2 * w9), 41 * PI) and np.isclose(np.sum(r2**2 * w9), 25 * PI)
print(f"  A2: ||f||² = 59π = {59*PI:.3f};  ||f-S_1f||² = 41π = {41*PI:.3f};  ||f-S_2f||² = 25π = {25*PI:.3f}")

# --- A3: cos²x and sin³x (n = -3..3); a_2 of cos²x by (3.4).
n, c = fourier_coeffs(lambda x: np.cos(x) ** 2, 3)
assert np.allclose(c, [0, 0.25, 0, 0.5, 0, 0.25, 0], atol=1e-12)
n, c = fourier_coeffs(lambda x: np.sin(x) ** 3, 3)
assert np.allclose(c, [-1j / 8, 0, 3j / 8, 0, -3j / 8, 0, 1j / 8], atol=1e-12)
assert np.isclose(quad(lambda x: np.cos(x) ** 2 * np.cos(2 * x), -PI, PI)[0] / PI, 0.5)

# --- A4: sine series of 1 on (0, L): b_n = 4/(nπ), n odd; S_5 f(L/2) = 52/(15π).
nA4, bA4 = sine_coeffs(lambda x: np.ones_like(x), 6, L=1.0, M=2**16)
assert np.allclose(bA4, np.where(nA4 % 2 == 1, 4 / (nA4 * PI), 0), atol=1e-8)
S5_mid = np.sum(bA4[:5] * np.sin(nA4[:5] * PI / 2))
assert np.isclose(S5_mid, 52 / (15 * PI), atol=1e-8) and np.isclose(S5_mid, 1.1035, atol=5e-5)

# --- A5: Osgood's period-1 convention, f(t) = t: ĉ_n = i(-1)^n/(2πn), b_n = (-1)^{n+1}/(πn).
n, c = fourier_coeffs(lambda t_: np.asarray(t_, float), 3, T=1.0)
nz = n != 0
assert np.allclose(c[nz], 1j * (-1.0) ** n[nz] / (2 * PI * n[nz]), atol=1e-8) and abs(c[~nz][0]) < 1e-12
assert np.isclose(c[n == 1][0].imag, -0.1592, atol=5e-5)
kk = np.arange(1, 200001, dtype=float)
assert np.isclose(np.sum((-1.0) ** (kk + 1) / (PI * kk) * np.sin(PI * kk / 2)), 0.25, atol=1e-5)   # series at t = 1/4
print(f"  A4: S_5 f(L/2) = 52/(15π) = {S5_mid:.4f};  A5: ĉ_1 = {c[n == 1][0]:.4f}")

# --- B1: x(π - x) on [0, π]: b_n = 8/(πn³), n odd; Σ_odd (-1)^{(n-1)/2}/n³ = π³/32; Theorem 3.9 route.
nB1, bB1 = sine_coeffs(lambda x: x * (PI - x), 7, L=PI, M=2**16)
assert np.allclose(bB1, np.where(nB1 % 2 == 1, 8 / (PI * nB1**3), 0), atol=1e-8)
assert np.allclose(bB1[[0, 2, 4]], [2.5465, 0.0943, 0.0204], atol=5e-5)
kodd = np.arange(1, 100001, 2, dtype=float)
assert np.isclose(np.sum((-1.0) ** ((kodd - 1) // 2) / kodd**3), PI**3 / 32, atol=1e-10) and np.isclose(PI**3 / 32, 0.96895, atol=5e-6)
n, c = fourier_coeffs(lambda x: x * (PI - np.abs(x)), 7, M=2**16)       # odd extension; c_n = -4i/(πn³), n odd
odd = n % 2 != 0
assert np.allclose(c[odd], -4j / (PI * n[odd] ** 3), atol=1e-8) and np.allclose(c[~odd], 0, atol=1e-8)
print(f"  B1: b_1, b_3, b_5 = {bB1[0]:.4f}, {bB1[2]:.4f}, {bB1[4]:.4f};  π³/32 = {PI**3/32:.5f}")

# --- B2(d): square wave. Mean-square error -> 0; overshoot max S_N f - 1 -> 0.179; sup|f - S_N f| = 1 always.
mse_B2 = {}
for N in (1, 5, 11, 51):
    mse_B2[N] = 1 - 8 / PI**2 * np.sum(1 / np.arange(1, N + 1, 2) ** 2)
    nN = np.arange(-N, N + 1)
    resid = square(xg) - partial_sum(nN, exact_c("square", nN), xg)
    assert abs(np.trapezoid(resid**2, xg) / (2 * PI) - mse_B2[N]) < 1e-3          # (1/2π)∫|f - S_N f|² directly
    assert np.abs(resid).max() > 0.98                                              # sup error ≈ 1 next to the jump
assert np.allclose([mse_B2[N] for N in (1, 5, 11, 51)], [0.1894, 0.0669, 0.0337, 0.0078], atol=5e-5)
assert np.allclose([over["square"][N] for N in (1, 5, 11, 51)], [0.2732, 0.1884, 0.1813, 0.1791], atol=5e-5)
print("  B2(d): mean-square error", {N: round(float(v), 4) for N, v in mse_B2.items()},
      " overshoot max S_N f - 1", {N: round(float(over["square"][N]), 4) for N in (1, 5, 11, 51)}, " sup|f - S_N f| -> 1")

# --- B3: zeta values from Parseval.
kk = np.arange(1, 100001, dtype=float)                                                  # float: k**6 overflows int64
assert np.isclose(np.sum(1 / kk[::2] ** 4), PI**4 / 96, atol=1e-9) and np.isclose(PI**4 / 96, 1.014678, atol=5e-7)
assert np.isclose(np.sum(1 / kk[1::2] ** 4), PI**4 / 1440, atol=1e-9) and np.isclose(PI**4 / 1440, 0.067645, atol=5e-7)
assert np.isclose(2 / PI * quad(lambda x: (x * (PI - x)) ** 2, 0, PI)[0], PI**4 / 15)       # half-range Parseval, B1
assert np.isclose(np.sum(1 / kk[::2] ** 6), PI**6 / 960, atol=1e-12) and np.isclose(PI**6 / 960, 1.001447, atol=5e-7)
assert np.isclose(np.sum(1 / kk**6), PI**6 / 945, atol=1e-12) and np.isclose(PI**6 / 945, 1.017343, atol=5e-7)

# --- B4: (π² - x²)²: f''' = 24x jumps by -48π; c_n = -24(-1)^n/n⁴ (B8 above checks it by trapezoid).
assert np.isclose(quad(lambda x: (PI**2 - x**2) ** 2, -PI, PI)[0] / (2 * PI), 8 * PI**4 / 15) and np.isclose(8 * PI**4 / 15, 51.9515, atol=5e-5)
nn9 = np.arange(1, 2001, dtype=float)
assert abs(8 * PI**4 / 15 - 48 * np.sum(1 / nn9**4)) < 1e-8                                  # series at x = π: f(π) = 0
assert np.isclose(np.sum((-1.0) ** (nn9 + 1) / nn9**4), 7 * PI**4 / 720, atol=1e-12) and np.isclose(7 * PI**4 / 720, 0.947033, atol=5e-7)
assert np.isclose((48 * PI + 24 * 2 * PI) / (2 * PI), 48)                                   # Theorem 3.10's C, vs actual 24

# --- B5: Dirichlet kernel: zeros, ∫D_N² = 2π(2N+1), Lebesgue constants and their log growth.
tk, wk = trap_weights(-PI, PI, 2**18)
for N in (3, 10):
    assert np.allclose(dirichlet_kernel(N, 2 * np.arange(1, 2 * N + 1) * PI / (2 * N + 1)), 0, atol=1e-9)
    assert np.isclose(np.sum(dirichlet_kernel(N, tk) ** 2 * wk), 2 * PI * (2 * N + 1), rtol=1e-10)
assert np.isclose(42 * PI, 131.95, atol=5e-3)
LN = {N: np.sum(np.abs(dirichlet_kernel(N, tk)) * wk) / (2 * PI) for N in (1, 10, 100)}   # (1/2π)∫|D_N|
bound = {N: 4 / PI**2 * np.log(N + 1) for N in LN}
assert np.allclose([LN[N] for N in (1, 10, 100)], [1.436, 2.223, 3.139], atol=5e-4) and all(LN[N] > bound[N] for N in LN)
assert np.allclose([bound[N] for N in (1, 10, 100)], [0.281, 0.972, 1.870], atol=5e-4)
assert np.allclose([LN[N] - 4 / PI**2 * np.log(N) for N in (10, 100)], [1.29, 1.27], atol=5e-3)
f_ext = np.sign(dirichlet_kernel(10, tk))                                                 # (e): sgn D_N attains L_N
assert np.isclose(np.sum(f_ext * dirichlet_kernel(10, -tk) * wk) / (2 * PI), LN[10], atol=1e-6)
print(f"  B5: L_N = {LN[1]:.3f}, {LN[10]:.3f}, {LN[100]:.3f} (N = 1, 10, 100) > (4/π²)ln(N+1) = "
      f"{bound[1]:.3f}, {bound[10]:.3f}, {bound[100]:.3f}")

# --- B6: termwise integration: sawtooth -> x²/2, then x² - π²/3 -> x(π² - x²); termwise derivative of x².
xa = np.linspace(-PI, PI, 2001)
assert np.abs(2 * ((-1.0) ** (nn9 + 1) / nn9**2) @ (1 - np.cos(np.outer(nn9, xa))) - xa**2 / 2).max() < 5e-3
cubic_series = 12 * ((-1.0) ** (nn9 + 1) / nn9**3) @ np.sin(np.outer(nn9, xa))
assert np.abs(cubic_series - cubic(xa)).max() < 1e-5                                      # (b): uniform
assert np.isclose(cubic(PI / 2), 3 * PI**3 / 8) and np.isclose(12 * PI**3 / 32, 3 * PI**3 / 8)   # (c)
deriv_x2 = -4 * ((-1.0) ** nn9 / nn9) @ np.sin(np.outer(nn9, xa))                         # (d)
inner = np.abs(xa) < PI - 0.3
assert np.abs(deriv_x2[inner] - 2 * xa[inner]).max() < 0.05 and abs(deriv_x2[-1]) < 1e-6

# --- B7: insulated ends. Cosine coefficients of the step; mean conserved; u -> u_0/2.
xb = np.linspace(0, 1, 2**14 + 1); wb = np.full(xb.size, 1 / 2**14); wb[0] = wb[-1] = 0.5 / 2**14
stepB7 = np.where(xb < 0.5, 1.0, np.where(xb > 0.5, 0.0, 0.5))
aB7 = 2 * np.cos(np.outer(np.arange(7), PI * xb)) @ (stepB7 * wb)                         # a_0..a_6 (L = 1, u_0 = 1)
assert np.isclose(aB7[0] / 2, 0.5, atol=1e-8)
assert np.allclose(aB7[1:], [0.6366, 0.0, -0.2122, 0.0, 0.1273, 0.0], atol=5e-5)
nB7 = np.arange(1, 401)
u_neu = lambda x, t_: 0.5 + (2 * np.sin(nB7 * PI / 2) / (nB7 * PI) * np.exp(-(nB7 * PI) ** 2 * t_)) @ np.cos(np.outer(nB7, PI * x))
for t_ in (0.01, 0.05, 0.2):
    assert np.isclose(np.sum(u_neu(xb, t_) * wb), 0.5, atol=1e-6)                         # ∫u dx conserved
assert np.abs(u_neu(xb, 5.0) - 0.5).max() < 1e-12                                         # u -> mean
print("  B7: Neumann a_1..a_6 =", np.round(aB7[1:], 4))

# --- C1: Gibbs for the sawtooth. (b) T_N peaks at π/(N+1); (c) Riemann sum; (d) overshoot above f;
#     (e) E_N = T_N - (π - t) peaks at 2π/(2N+1) and obeys |E_N| < 2/((N+1/2) sin(δ/2)) on [δ, 2π-δ].
T_N = lambda N, t_: 2 * np.sin(np.outer(np.atleast_1d(t_), np.arange(1, N + 1))) @ (1 / np.arange(1, N + 1))
C1v = {}
for N in (11, 51):
    tN = PI / (N + 1)
    assert abs(dirichlet_kernel(N, tN) - 1) < 1e-9                                         # T_N'(t_N) = D_N(t_N) - 1 = 0
    tt = np.linspace(1e-4, PI, 40001)
    assert abs(tt[np.argmax(T_N(N, tt))] - tN) < 1e-4
    E = T_N(N, tt) - (PI - tt)
    assert abs(tt[np.argmax(E)] - 2 * PI / (2 * N + 1)) < 1e-4                             # (e) step 4
    for d in (0.3, 1.0):
        ts = np.linspace(d, 2 * PI - d, 4001)
        assert np.abs(T_N(N, ts) - (PI - ts)).max() < 2 / ((N + 0.5) * np.sin(d / 2))     # (e) step 3
    val = T_N(N, tN)[0]
    u_ = np.arange(1, N + 1) * tN
    assert np.isclose(2 * np.sum(np.sin(u_) / u_) * tN, val)                               # (c): Riemann sum
    C1v[N] = (val, PI - tN, val - (PI - tN))
assert np.allclose(C1v[11], [3.4384, 2.8798, 0.5586], atol=5e-5) and np.allclose(C1v[51], [3.6433, 3.0812, 0.5621], atol=5e-5)
assert np.isclose(2 * Si_pi, 3.7038741, atol=1e-7) and np.isclose(2 * Si_pi - PI, 0.5623, atol=5e-5)
gap = {N: over["sawtooth"][N] - C1v[N][2] for N in (11, 51)}                               # grid max minus value at t_N
assert np.isclose(gap[11], 0.0057, atol=3e-4) and np.isclose(gap[51], 0.0003, atol=3e-4) and 10 < gap[11] / gap[51] < 40
print(f"  C1: T_11(π/12) - (π - π/12) = {C1v[11][2]:.4f}, T_51(π/52) - (π - π/52) = {C1v[51][2]:.4f};  "
      f"grid maxima {over['sawtooth'][11]:.4f}, {over['sawtooth'][51]:.4f};  limit 2Si(π) - π = {2*Si_pi - PI:.4f}")

# --- C2(e): isoperimetric check on the ellipse x = 2cos θ, y = sin θ, and the area formula (c) by FFT.
L_ell = quad(lambda th: np.sqrt(4 * np.sin(th) ** 2 + np.cos(th) ** 2), 0, 2 * PI)[0]
assert np.isclose(L_ell, 9.6884, atol=5e-5) and np.isclose(L_ell**2 / (4 * PI), 7.4696, atol=5e-5) and 2 * PI < L_ell**2 / (4 * PI)
Mc = 256; th = 2 * PI * np.arange(Mc) / Mc; nf = np.fft.fftfreq(Mc, d=1 / Mc)
alpha = np.fft.fft(2 * np.cos(th)) / Mc; beta = np.fft.fft(np.sin(th)) / Mc               # α_n, β_n by (5.1)
assert np.isclose(2 * PI * np.sum(1j * nf * np.conj(alpha) * beta), 2 * PI, atol=1e-10)   # (c): A = 2π
assert np.isclose(np.sum(nf**2 * (np.abs(alpha) ** 2 + np.abs(beta) ** 2)), 2.5, atol=1e-10)   # θ is not arc length
print(f"  C2: ellipse L = {L_ell:.4f}, L²/4π = {L_ell**2/(4*PI):.4f} > A = 2π = {2*PI:.4f}")

# --- C3: Fejér. (a) taper formula = average of S_0..S_N; (b) closed-form kernel; (c) F_N >= 0, unit mass; (d) no overshoot.
def fejer_sum(n, c, x):
    """σ_N f(x) = Σ_{|n|<=N} (1 - |n|/(N+1)) c_n e^{inx}, the average of S_0 f, ..., S_N f (C3(a))."""
    N = n.max()
    return np.real(np.exp(1j * np.outer(x, n)) @ ((1 - np.abs(n) / (N + 1)) * c))


def fejer_kernel(N, t_):
    """F_N(t) = sin²((N+1)t/2) / ((N+1) sin²(t/2)), F_N(0) = N+1 (C3(b))."""
    s = np.sin(np.asarray(t_, float) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.abs(s) > 1e-12, np.sin((N + 1) * t_ / 2) ** 2 / s**2 / (N + 1), N + 1.0)


N = 11; nN = np.arange(-N, N + 1)
sig_avg = np.mean([partial_sum(np.arange(-j, j + 1), exact_c("square", np.arange(-j, j + 1)), xg) for j in range(N + 1)], axis=0)
assert np.allclose(fejer_sum(nN, exact_c("square", nN), xg), sig_avg, atol=1e-10)                          # (a)
assert np.allclose(fejer_kernel(N, tk), np.mean([dirichlet_kernel(j, tk) for j in range(N + 1)], axis=0), atol=1e-9)   # (b)
assert fejer_kernel(N, tk).min() >= 0 and np.isclose(np.sum(fejer_kernel(N, tk) * wk) / (2 * PI), 1.0, atol=1e-10)    # (c)
sig_max = {}
for N in (11, 51):
    nN = np.arange(-N, N + 1)
    sig = fejer_sum(nN, exact_c("square", nN), xg)
    assert sig.max() <= 1 + 1e-12 and sig.min() >= -1 - 1e-12                                                # (d)
    sig_max[N] = sig.max()
assert np.allclose([sig_max[11], sig_max[51]], [0.9473, 0.9878], atol=5e-5)
N = 51; nN = np.arange(-N, N + 1); cN = exact_c("square", nN)                                                # (f)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(xg, square(xg), "k", lw=1, label="f")
ax.plot(xg, partial_sum(nN, cN, xg), lw=1, label=f"S_{N} f (Dirichlet): max {1 + over['square'][N]:.4f}")
ax.plot(xg, fejer_sum(nN, cN, xg), lw=1, label=f"σ_{N} f (Fejér): max {sig_max[N]:.4f}")
ax.set_xlim(-1, 1); ax.set_ylim(-1.3, 1.3); ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_xlabel("x")
ax.set_title("Square wave: the partial sum overshoots, the Fejér mean never does")
fig.tight_layout(); fig.savefig(FIG / "fejer_vs_dirichlet.png", dpi=120); plt.close(fig)
print(f"  C3: max σ_11 f = {sig_max[11]:.4f}, max σ_51 f = {sig_max[51]:.4f} (no overshoot); S_N peaks near 1.18")
print("Part 9 passed: every number in solutions.md reproduced.")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
