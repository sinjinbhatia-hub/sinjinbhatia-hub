# %% [markdown]
# # F3 Project: Transform, DFT/FFT, Circulants, Oscillator, Sampling, Spectrogram
#
# This script builds the whole module and checks every number quoted in the
# notes. Parts: (1) the DFT matrix, a matrix DFT and a recursive radix-2 FFT,
# unitarity of F/sqrt(N), the hand-computed DFTs of Examples 4.5-4.6;
# (2) timing N^2 vs N log N; (3) circular convolution and circulant
# diagonalization; (4) continuous transforms of rect, Gaussians, exponentials,
# triangle by the Riemann-sum relation (3.29), overlaid on the analytic formulas,
# the Example 4.1-4.3 numbers by quadrature and by FFT, and the convolution
# theorem; (5) the uncertainty product; (6) the damped
# oscillator: impulse response, |H|, resonance, chirp response; (7) aliasing and
# sinc interpolation; (8) Poisson summation and the comb; (9) a spectrogram.
# "Your turn" reference solutions (B7, C3) are at the end.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.perf_counter()

# %% [markdown]
# ## Part 1: DFT matrix, matrix DFT, recursive FFT, unitarity
#
# (3.26): F[k] = sum_n f[n] w^{kn}, w = exp(-2 pi i / N). Theorem 3.19: F^* F = N I.
# Theorem 3.26 (butterfly): F[k] = E[k] + w^k O[k], F[k+N/2] = E[k] - w^k O[k].
#
# **Build it yourself first:** write `fft_recursive` from (3.31) before reading it.


# %%
def dft_matrix(N):
    """The N x N Fourier matrix F_{kn} = w^{kn}, w = exp(-2 pi i/N) (Definition 3.18)."""
    k = np.arange(N)
    return np.exp(-2j * np.pi * np.outer(k, k) / N)


def dft(x):
    """DFT by matrix multiplication, eq. (3.26). O(N^2)."""
    x = np.asarray(x, complex)
    return dft_matrix(len(x)) @ x


def dft_blockwise(x, block=256):
    """Same O(N^2) arithmetic as `dft`, but building the matrix a few rows at a time
    (so N = 2^14 fits in memory). Used only for timing."""
    x = np.asarray(x, complex)
    N = len(x)
    roots = np.exp(-2j * np.pi * np.arange(N) / N)  # w^j for j = 0..N-1
    n = np.arange(N)
    out = np.empty(N, complex)
    for k0 in range(0, N, block):
        kk = np.arange(k0, min(k0 + block, N))
        out[k0:k0 + len(kk)] = roots[(kk[:, None] * n[None, :]) % N] @ x
    return out


def fft_recursive(x):
    """Radix-2 Cooley-Tukey FFT, Theorem 3.26 / eq. (3.31). Requires N = 2^m."""
    x = np.asarray(x, complex)
    N = len(x)
    if N == 1:
        return x.copy()
    assert N % 2 == 0, "radix-2 FFT needs N = 2^m (zero-pad otherwise)"
    E = fft_recursive(x[0::2])          # DFT of even samples, length N/2
    O = fft_recursive(x[1::2])          # DFT of odd samples
    twiddle = np.exp(-2j * np.pi * np.arange(N // 2) / N) * O
    return np.concatenate([E + twiddle, E - twiddle])   # the butterflies


for N in (1, 2, 4, 8, 16, 64, 256, 1024):
    x = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    ref = np.fft.fft(x)
    assert np.max(np.abs(dft(x) - ref)) < 1e-10, N
    assert np.max(np.abs(fft_recursive(x) - ref)) < 1e-10, N
    assert np.max(np.abs(dft_blockwise(x) - ref)) < 1e-10, N
    # Theorem 3.19: F/sqrt(N) is unitary; F^{-1} = conj(F)/N.
    F = dft_matrix(N)
    U = F / np.sqrt(N)
    assert np.allclose(U.conj().T @ U, np.eye(N), atol=1e-10)
    assert np.allclose(np.linalg.inv(F), F.conj() / N)
    # Parseval for the DFT.
    assert np.isclose(np.sum(np.abs(x) ** 2), np.sum(np.abs(ref) ** 2) / N)
print("Part 1: dft, dft_blockwise, fft_recursive agree with np.fft.fft to 1e-10; F/sqrt(N) unitary.")

# F^2 = N P (reversal) and F^4 = N^2 I; eigenvalues of F/sqrt(N) in {1,-1,i,-i} (C2).
N = 8
F = dft_matrix(N)
P = np.zeros((N, N))
for k in range(N):
    P[k, (-k) % N] = 1
assert np.allclose(F @ F, N * P) and np.allclose(F @ F @ F @ F, N**2 * np.eye(N))
lam = np.linalg.eigvals(F / np.sqrt(N))
assert np.allclose(np.sort_complex(np.round(lam, 8)), np.sort_complex(np.round(np.array([1, 1, 1, -1, -1, 1j, -1j, -1j]), 8)))
print("  F^2 = N P, F^4 = N^2 I; eigenvalues of F/sqrt(8):", np.round(lam, 6))

# Example 4.5: small DFTs by hand, integer vs non-integer periods (leakage).
assert np.allclose(dft([1, 0, 0, 0]), [1, 1, 1, 1])
assert np.allclose(dft([1, 1, 1, 1]), [4, 0, 0, 0])
n8 = np.arange(8)
F_int = dft(np.cos(2 * np.pi * 2 * n8 / 8))
assert np.allclose(F_int, [0, 0, 4, 0, 0, 0, 4, 0], atol=1e-12)
F_leak = dft(np.cos(2 * np.pi * 2.5 * n8 / 8))
leak_expected = [1.000, 1.192, 2.798, 2.398, 1.000, 2.398, 2.798, 1.192]
assert np.allclose(np.abs(F_leak), leak_expected, atol=5e-4)
print("  Example 4.5: |DFT(cos 2.5 cycles)| =", np.round(np.abs(F_leak), 3))

# Example 4.6: the length-8 FFT by hand.
s2 = np.sqrt(2)
hand = np.array([10, (1 - s2) - 1j * (3 + 3 * s2), -2 + 2j, (1 + s2) + 1j * (3 - 3 * s2),
                 -2, (1 + s2) - 1j * (3 - 3 * s2), -2 - 2j, (1 - s2) + 1j * (3 + 3 * s2)])
assert np.allclose(fft_recursive([1, 2, 3, 4, 0, 0, 0, 0]), hand)
assert np.allclose(np.fft.fft([1, 2, 3, 4, 0, 0, 0, 0]), hand)
print("  Example 4.6: FFT([1,2,3,4,0,0,0,0]) matches the hand computation.")

# Exercise A4.
assert np.allclose(dft([1, -1, 1, -1]), [0, 0, 4, 0]) and np.allclose(dft([0, 1, 0, -1]), [0, -2j, 0, 2j])

# %% [markdown]
# ## Part 2: Timing N^2 versus N log N
#
# Eq. (3.32): the FFT needs (N/2) log2 N complex multiplications; the matrix method N^2.

# %%
Ns = [2**p for p in range(6, 15)]
t_mat, t_fft, t_np = [], [], []
for N in Ns:
    x = rng.standard_normal(N) + 0j
    reps = 3 if N <= 2**10 else 1
    t = time.perf_counter()
    for _ in range(reps):
        y1 = dft_blockwise(x)
    t_mat.append((time.perf_counter() - t) / reps)
    t = time.perf_counter()
    for _ in range(reps):
        y2 = fft_recursive(x)
    t_fft.append((time.perf_counter() - t) / reps)
    t = time.perf_counter()
    for _ in range(20):
        y3 = np.fft.fft(x)
    t_np.append((time.perf_counter() - t) / 20)
    assert np.max(np.abs(y1 - y3)) < 1e-8 and np.max(np.abs(y2 - y3)) < 1e-8
assert t_mat[-1] > 3 * t_fft[-1], "the O(N^2) method should be much slower at N = 2^14"
# Crossover (§7 item 2): the vectorized matrix method wins at small N because the recursion pays
# Python call overhead at every level; the FFT overtakes it around N = 2^9-2^10.
N_cross = next(N for N, tm, tf in zip(Ns, t_mat, t_fft) if tf < tm)
assert N_cross <= 2**11 and t_mat[Ns.index(2**12)] > 2 * t_fft[Ns.index(2**12)], (N_cross, t_mat, t_fft)
# Operation counts quoted in §1 and §3.14 from (3.32).
assert 2**10 // 2 * 10 == 5120 and 2**20 == 1_048_576 and 200 < 2**20 / 5120 < 210          # N = 2^10: "a factor of 200"
assert 2**40 / (2**20 * 20) == 52428.8                                                     # §1: N log2 N vs N^2 at 2^20, "50,000"

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.loglog(Ns, t_mat, "o-", label="matrix DFT, O(N^2)")
ax.loglog(Ns, t_fft, "s-", label="recursive FFT (pure Python), O(N log N)")
ax.loglog(Ns, t_np, "^-", label="np.fft.fft")
Nn = np.array(Ns, float)
ax.loglog(Ns, t_mat[-1] * (Nn / Nn[-1]) ** 2, "k:", lw=1, label="slope N^2")
ax.loglog(Ns, t_fft[-1] * (Nn * np.log2(Nn)) / (Nn[-1] * np.log2(Nn[-1])), "k--", lw=1, label="slope N log N")
ax.set_xlabel("N"); ax.set_ylabel("seconds"); ax.set_title("DFT cost: N^2 vs N log N")
ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "fft_timing.png", dpi=120); plt.close(fig)
print(f"Part 2: N=2^14: matrix {t_mat[-1]:.3f} s, recursive FFT {t_fft[-1]:.4f} s, np.fft {t_np[-1]*1e3:.3f} ms; "
      f"recursive FFT first beats the matrix method at N = {N_cross}.")

# %% [markdown]
# ## Part 3: Circular convolution, circulant matrices, C = F^{-1} Lambda F
#
# Definition 3.21, Theorem 3.22 (DFT(f (*) g) = F G), Theorem 3.24 (eigenvalues of a
# circulant are the DFT of its first column), Theorem 3.25 (CS = SC iff circulant).


# %%
def circular_convolution(f, g):
    """(f (*) g)[n] = sum_m f[m] g[(n-m) mod N] via the DFT, Theorem 3.22."""
    f = np.asarray(f, complex); g = np.asarray(g, complex)
    assert f.shape == g.shape
    out = np.fft.ifft(np.fft.fft(f) * np.fft.fft(g))
    return out.real if np.isrealobj(f) and np.isrealobj(g) else out


def circular_convolution_direct(f, g):
    """The same sum, done directly, O(N^2)."""
    N = len(f)
    return np.array([sum(f[m] * g[(n - m) % N] for m in range(N)) for n in range(N)])


def circulant(c):
    """Circulant matrix with first column c: C_{jk} = c[(j-k) mod N] (Definition 3.23)."""
    c = np.asarray(c)
    N = len(c)
    j, k = np.indices((N, N))
    return c[(j - k) % N]


def cyclic_shift(N):
    """(S x)[n] = x[n-1 mod N]: ones on the subdiagonal and in the top-right corner."""
    S = np.zeros((N, N))
    for n in range(N):
        S[n, (n - 1) % N] = 1
    return S


N = 16
c = rng.standard_normal(N); x = rng.standard_normal(N)
assert np.allclose(circular_convolution(c, x), circular_convolution_direct(c, x))
C = circulant(c)
assert np.allclose(C @ x, circular_convolution(c, x))                 # convolution = circulant times vector
assert np.allclose(C, sum(c[m] * np.linalg.matrix_power(cyclic_shift(N), m) for m in range(N)))  # C = sum c[m] S^m
S = cyclic_shift(N)
assert np.allclose(C @ S, S @ C)                                      # Theorem 3.25
lam_dft = np.fft.fft(c)
lam_eig = np.linalg.eigvals(C)
assert np.allclose(np.sort_complex(np.round(lam_dft, 9)), np.sort_complex(np.round(lam_eig, 9)))
F = dft_matrix(N)
assert np.allclose(np.linalg.inv(F) @ np.diag(lam_dft) @ F, C)         # (3.30)
# The discrete tones are eigenvectors of S (eigenvalue w^k) and of C (eigenvalue DFT(c)[k]).
w = np.exp(-2j * np.pi / N)
for k in range(N):
    e_k = np.exp(2j * np.pi * k * np.arange(N) / N)
    assert np.allclose(S @ e_k, w**k * e_k) and np.allclose(C @ e_k, lam_dft[k] * e_k)
# Exercise B3: the ring Laplacian.
C_lap = circulant([2, -1, 0, -1])
assert np.allclose(np.fft.fft([2, -1, 0, -1]).real, [0, 2, 4, 2])
assert np.allclose(C_lap @ np.array([1, 1j, -1, -1j]), 2 * np.array([1, 1j, -1, -1j]))
# Linear convolution by zero padding.
a = rng.standard_normal(5); b = rng.standard_normal(7)
pad = len(a) + len(b) - 1
lin = circular_convolution(np.r_[a, np.zeros(pad - 5)], np.r_[b, np.zeros(pad - 7)])
assert np.allclose(lin, np.convolve(a, b))
print("Part 3: circulant eigenvalues = DFT(first column); C = F^-1 Lambda F; circular conv via FFT = direct; zero-padded = np.convolve.")

# %% [markdown]
# ## Part 4: Continuous transforms by the DFT (Theorem 3.20), overlaid on the formulas
#
# Samples f(x0 + n dx), n = 0..N-1. Then (3.29) with the start-point phase correction:
# f^(k/(N dx)) ~ dx * exp(-2 pi i x0 xi_k) * F[k].


# %%
def ft_numeric(samples, x0, dx):
    """Approximate f^(xi) on the DFT frequency grid from samples f(x0 + n dx), eq. (3.29)
    plus the shift rule (3.6) for the record start x0. Returns (xi, f^) with xi increasing."""
    N = len(samples)
    xi = np.fft.fftfreq(N, d=dx)                       # k/(N dx), negative frequencies in the second half
    Fh = dx * np.fft.fft(samples) * np.exp(-2j * np.pi * x0 * xi)
    return np.fft.fftshift(xi), np.fft.fftshift(Fh)


def rect(x):
    """Pi(x): 1 on |x|<1/2, 1/2 at the jumps (so the Riemann sum is a trapezoid rule there)."""
    return np.where(np.abs(x) < 0.5, 1.0, np.where(np.abs(x) == 0.5, 0.5, 0.0))


def sinc(x):
    return np.sinc(x)   # numpy's sinc is sin(pi x)/(pi x), our convention


L, dx = 16.0, 1 / 128
x = np.arange(-L, L, dx)
xi, _ = ft_numeric(np.zeros_like(x), -L, dx)
band = np.abs(xi) <= 6
cases = {
    "rect -> sinc (3.3)": (rect(x), sinc(xi), 2e-3),
    "exp(-pi x^2) -> exp(-pi xi^2) (3.11)": (np.exp(-np.pi * x**2), np.exp(-np.pi * xi**2), 1e-10),
    "exp(-x^2/2), sigma=1 (Ex. 4.2)": (np.exp(-x**2 / 2), np.sqrt(2 * np.pi) * np.exp(-2 * np.pi**2 * xi**2), 1e-10),
    "exp(-2x)u(x) -> 1/(2+2 pi i xi) (3.4)": (np.where(x > 0, np.exp(-2 * x), np.where(x == 0, 0.5, 0.0)), 1 / (2 + 2j * np.pi * xi), 3e-3),
    "exp(-2|x|) -> 4/(4+4pi^2 xi^2) (3.5)": (np.exp(-2 * np.abs(x)), 4 / (4 + 4 * np.pi**2 * xi**2), 1e-4),
    "triangle -> sinc^2": (np.clip(1 - np.abs(x), 0, None), sinc(xi) ** 2, 1e-4),
}
fig, axes = plt.subplots(2, 3, figsize=(13, 7))
for ax, (name, (fx, analytic, tol)) in zip(axes.ravel(), cases.items()):
    _, Fh = ft_numeric(fx, -L, dx)
    err = np.max(np.abs(Fh[band] - analytic[band]))
    assert err < tol, (name, err)
    ax.plot(xi[band], analytic[band].real, "k-", lw=2, label="analytic (Re)")
    ax.plot(xi[band], Fh[band].real, "C1--", label="DFT via (3.29)")
    if np.iscomplexobj(analytic):
        ax.plot(xi[band], analytic[band].imag, "k-", lw=1, alpha=.5)
        ax.plot(xi[band], Fh[band].imag, "C0--", label="imag")
    ax.set_title(f"{name}\nmax err {err:.1e}", fontsize=9); ax.grid(alpha=.3); ax.legend(fontsize=7)
    ax.set_xlabel("xi")
fig.tight_layout(); fig.savefig(FIG / "transforms_overlay.png", dpi=120); plt.close(fig)
print("Part 4: rect, Gaussians, exponentials, triangle: DFT via (3.29) matches the analytic transforms.")


def ft_quadrature(samples, x0, dx, xi0):
    """f^(xi0) at one arbitrary frequency by the Riemann sum dx * sum_n f(x_n) exp(-2 pi i x_n xi0),
    x_n = x0 + n dx: the sum behind (3.29), evaluated OFF the DFT grid. O(N) per frequency."""
    xn = x0 + dx * np.arange(len(samples))
    return dx * np.sum(samples * np.exp(-2j * np.pi * xn * xi0))


# Example 4.2 numbers (sigma = 1). The DFT grid here is k/(N dx) = k/32, spacing 0.03125, so
# xi = 0.1 is NOT a grid point: the nearest bin (0.09375) holds f^(0.09375) = 2.1074, not f^(0.1).
# Read f^ at exactly 0.1 by quadrature (the Riemann sum at that xi), which is what the notes claim.
gauss1 = np.exp(-x**2 / 2)
_, Fg = ft_numeric(gauss1, -L, dx)
for xi0, val in ((0.0, 2.50663), (0.1, 2.05761), (0.25, 0.72996)):
    num = ft_quadrature(gauss1, -L, dx, xi0).real
    assert abs(num - val) < 1e-5, (xi0, num, val)
    assert abs(num - np.sqrt(2 * np.pi) * np.exp(-2 * np.pi**2 * xi0**2)) < 1e-10   # Poisson summation makes it exact
i_near = np.argmin(np.abs(xi - 0.1))
print(f"  Example 4.2 (quadrature): f^(0) = {ft_quadrature(gauss1, -L, dx, 0.0).real:.5f}, "
      f"f^(0.1) = {ft_quadrature(gauss1, -L, dx, 0.1).real:.5f}, f^(0.25) = {ft_quadrature(gauss1, -L, dx, 0.25).real:.5f}; "
      f"nearest DFT bin to 0.1 is xi = {xi[i_near]:.5f} with value {Fg[i_near].real:.4f} (grid spacing {xi[1] - xi[0]:.5f})")

# Example 4.1: the triangle at xi = 0.3, and a 2 ms pulse: f^(0) = 0.002, nulls at 500 and 1000 Hz.
assert np.isclose(sinc(0.3) ** 2, 0.73684, atol=5e-6)
assert np.isclose((1 - np.cos(0.6 * np.pi)) / (2 * np.pi**2 * 0.09), 0.73684, atol=5e-6)
tau_p, t0_p, dt_p = 2e-3, 3e-3, 1e-5                                   # width 2 ms, centre 3 ms, sampled at 10 us
tp = np.arange(2000) * dt_p                                            # 20 ms record starting at 0: frequency grid 50 Hz
d = np.abs(tp - t0_p)
pulse = np.where(d < tau_p / 2 - 1e-9, 1.0, np.where(d < tau_p / 2 + 1e-9, 0.5, 0.0))   # 1/2 at the two jumps
Gp = dt_p * np.fft.rfft(pulse)                                         # (3.29); record starts at 0, no phase correction
xi_p = np.fft.rfftfreq(len(tp), d=dt_p)
kp = lambda nu: int(np.argmin(np.abs(xi_p - nu)))
assert np.isclose(abs(Gp[0]), 0.002) and abs(Gp[kp(500)]) < 1e-9 and abs(Gp[kp(1000)]) < 1e-9
# Example 4.3: the delay leaves |G| = tau |sinc(tau xi)| alone and makes the phase -2 pi t0 xi.
sel_p = xi_p <= 2000
G_an = np.exp(-2j * np.pi * t0_p * xi_p[sel_p]) * tau_p * sinc(tau_p * xi_p[sel_p])
err_p = np.max(np.abs(Gp[sel_p] - G_an))
assert err_p < 5e-7, err_p            # measured 1.8e-7 up to 2 kHz: trapezoid error of a 10 us grid, 1e-4 of the 0.002 peak
assert np.isclose(abs(Gp[kp(100)]), 1.871e-3, atol=1e-6) and np.isclose(abs(Gp[kp(100)]), 0.002 * sinc(0.2), atol=1e-8)
assert np.isclose(np.degrees(np.angle(Gp[kp(100)])), -108.0, atol=1e-6)
print(f"  Example 4.1/4.3: sinc^2(0.3) = {sinc(0.3)**2:.5f}; 2 ms pulse delayed 3 ms, from the FFT at 100 Hz: "
      f"|G| = {abs(Gp[kp(100)]):.3e}, phase = {np.degrees(np.angle(Gp[kp(100)])):.1f} deg; nulls at 500, 1000 Hz.")

# Convolution theorem: Pi * Pi = triangle, transform = sinc^2 (Theorem 3.7).
# np.convolve's full output lives on the grid x0 + x0 + m dx; the samples at our x[j] = x0 + j dx are
# m = j + N/2. (mode="same" on an even-length grid is off by one sample, which shifts the phase by
# exp(-2 pi i xi dx) and breaks the frequency-domain check below; slice explicitly instead.)
Nx = len(x)
tri_num = np.convolve(rect(x), rect(x))[Nx // 2:Nx // 2 + Nx] * dx
assert np.max(np.abs(tri_num - np.clip(1 - np.abs(x), 0, None))) < dx     # Riemann-sum error at the peak, dx/2
_, Ftri = ft_numeric(tri_num, -L, dx)
_, Frect = ft_numeric(rect(x), -L, dx)
assert np.max(np.abs(Ftri[band] - sinc(xi[band]) ** 2)) < 1e-4
assert np.max(np.abs(Ftri[band] - Frect[band] ** 2)) < 1e-12         # (f*g)^ = f^ g^ EXACTLY on the grid (Theorem 3.22)
# The mirror (3.14): (f g)^ = f^ * g^ for two Gaussians, where everything is explicit.
g1, g2 = np.exp(-np.pi * x**2), np.exp(-2 * np.pi * x**2)
_, Fprod = ft_numeric(g1 * g2, -L, dx)
dxi = xi[1] - xi[0]
conv_hat = np.convolve(np.exp(-np.pi * xi**2), np.sqrt(0.5) * np.exp(-np.pi * xi**2 / 2))[Nx // 2:Nx // 2 + Nx] * dxi
assert np.max(np.abs(Fprod[band].real - conv_hat[band])) < 1e-12     # exact: Gaussians sum under Poisson summation
assert np.max(np.abs(Fprod[band].real - np.exp(-np.pi * xi[band]**2 / 3) / np.sqrt(3))) < 1e-12   # (3.12) with a = 3 pi
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(x, tri_num, label="Pi * Pi (np.convolve)"); ax[0].plot(x, np.clip(1 - np.abs(x), 0, None), "k--", label="triangle")
ax[0].set_xlim(-2, 2); ax[0].legend(); ax[0].grid(alpha=.3); ax[0].set_title("Convolution: rect * rect = triangle")
ax[1].plot(xi[band], Ftri[band].real, label="transform of Pi*Pi"); ax[1].plot(xi[band], sinc(xi[band])**2, "k--", label="sinc^2")
ax[1].plot(xi[band], Frect[band].real, ":", label="sinc = transform of Pi"); ax[1].legend(); ax[1].grid(alpha=.3)
ax[1].set_title("Convolution theorem: (Pi*Pi)^ = sinc^2")
fig.tight_layout(); fig.savefig(FIG / "convolution_theorem.png", dpi=120); plt.close(fig)
print("  Convolution theorem checked: (Pi*Pi)^ = sinc^2 = (Pi^)^2; (f g)^ = f^ * g^ for Gaussians.")

# %% [markdown]
# ## Part 5: The uncertainty product (Theorem 3.11)
#
# sigma_x^2 = int x^2 |f|^2, sigma_xi^2 = int xi^2 |f^|^2, for ||f|| = 1. Gaussian: 1/(4 pi).
# Two-sided exponential (B2): sqrt(2)/(4 pi). Rect: infinite, and grows like L/pi^2 on |xi| < L.


# %%
def spreads(fx, x, Fh, xi):
    fx = fx / np.sqrt(np.trapezoid(np.abs(fx) ** 2, x))
    Fh = Fh / np.sqrt(np.trapezoid(np.abs(Fh) ** 2, xi))
    return np.sqrt(np.trapezoid(x**2 * np.abs(fx) ** 2, x)), np.sqrt(np.trapezoid(xi**2 * np.abs(Fh) ** 2, xi))


xg = np.linspace(-12, 12, 60001); xig = xg.copy()
sx, sxi = spreads(np.exp(-np.pi * xg**2), xg, np.exp(-np.pi * xig**2), xig)
assert abs(sx * sxi - 1 / (4 * np.pi)) < 1e-6
print(f"Part 5: Gaussian sigma_x*sigma_xi = {sx*sxi:.7f}  (1/(4 pi) = {1/(4*np.pi):.7f})")
a = 1.0
xe = np.linspace(-60, 60, 600001)
sx, sxi = spreads(np.exp(-a * np.abs(xe)), xe, 2 * a / (a**2 + 4 * np.pi**2 * xe**2), xe)
assert abs(sx - 1 / (a * np.sqrt(2))) < 1e-6 and abs(sxi - a / (2 * np.pi)) < 1e-3
assert abs(sx * sxi - np.sqrt(2) / (4 * np.pi)) < 2e-3
print(f"  two-sided exponential: product = {sx*sxi:.5f}  (sqrt2/(4 pi) = {np.sqrt(2)/(4*np.pi):.5f}), 41% above the minimum")
for Lxi in (20, 80):
    xr = np.linspace(-1, 1, 20001); xir = np.linspace(-Lxi, Lxi, int(4000 * Lxi) + 1)
    sx, sxi = spreads(rect(xr), xr, sinc(xir), xir)
    print(f"  rect on |xi|<{Lxi}: sigma_x = {sx:.4f} (1/sqrt12 = {1/np.sqrt(12):.4f}), sigma_xi^2 = {sxi**2:.3f} (L/pi^2 = {Lxi/np.pi**2:.3f}), product {sx*sxi:.3f}")
    assert sx * sxi > 1 / (4 * np.pi) and abs(sxi**2 - Lxi / np.pi**2) < 0.05 * Lxi / np.pi**2

# %% [markdown]
# ## Part 6: The damped driven oscillator m x'' + c x' + k x = f(t)
#
# (3.23): H(w) = 1/(k - m w^2 + i c w). (3.24): peak at w^2 = k/m - c^2/(2m^2),
# height 1/(c sqrt(k/m - c^2/(4m^2))). (3.25): h(t) = exp(-ct/2m) sin(w_d t)/(m w_d).
# Example 4.4: m=1, c=0.2, k=4: w_peak = 1.99499, |H|max = 2.50313.

# %%
m, c, k = 1.0, 0.2, 4.0
H = lambda w: 1 / (k - m * w**2 + 1j * c * w)
w0 = np.sqrt(k / m); wd = np.sqrt(k / m - c**2 / (4 * m**2))
w_peak = np.sqrt(k / m - c**2 / (2 * m**2)); H_max = 1 / (c * np.sqrt(k / m - c**2 / (4 * m**2)))
assert np.isclose(w_peak, 1.99499, atol=1e-5) and np.isclose(H_max, 2.50313, atol=1e-5)
assert np.isclose(abs(H(w0)), 2.5) and np.isclose(abs(H(0)), 0.25)
assert np.isclose(abs(H(1)), 0.33260, atol=1e-5) and np.isclose(np.degrees(np.angle(H(1))), -3.814, atol=1e-3)
assert np.isclose(abs(H(3)), 0.19858, atol=1e-5) and np.isclose(np.degrees(np.angle(H(3))), -173.157, atol=1e-3)
assert np.isclose(np.degrees(np.angle(H(2))), -90.0)
ws = np.linspace(1.5, 2.5, 200001); hp = ws[np.abs(H(ws)) >= H_max / np.sqrt(2)]
assert np.isclose(hp[-1] - hp[0], 0.2005, atol=5e-4)
assert np.isclose(hp[0], 1.892, atol=1e-3) and np.isclose(hp[-1], 2.093, atol=1e-3)   # half-power edges, Example 4.4
assert np.isclose(abs(H(10)), 0.0104, atol=5e-5) and np.isclose(abs(H(10)), 1 / 96, rtol=3e-4)   # 1/(m w^2) roll-off
assert np.isclose(wd, 1.99750, atol=1e-5) and np.isclose(1 / (m * wd), 0.50063, atol=1e-5)       # (3.25) amplitude
assert np.isclose(H_max / abs(H(0)), 10.0, atol=0.02)                                           # Q ~ w0 m / c = 10
assert np.isclose(np.angle(H(1)), -0.0666, atol=5e-5) and np.isclose(np.angle(H(3)), -3.0222, atol=5e-5)   # B6(b), radians
print(f"Part 6: oscillator w_peak = {w_peak:.5f}, |H|max = {H_max:.5f}, half-power edges {hp[0]:.3f}, {hp[-1]:.3f} "
      f"(width {hp[-1]-hp[0]:.4f} rad/s); |H(10)| = {abs(H(10)):.4f}; w_d = {wd:.5f}, 1/(m w_d) = {1/(m*wd):.5f}")

# Impulse response by solve_ivp: x(0) = 0, x'(0+) = 1/m, then free motion.
dt, T_rec = 0.02, 80.0
t = np.arange(0, T_rec, dt)
sol = solve_ivp(lambda tt, y: [y[1], (-c * y[1] - k * y[0]) / m], (0, T_rec), [0.0, 1.0 / m],
                t_eval=t, rtol=1e-10, atol=1e-12, method="DOP853")
h_num = sol.y[0]
h_formula = np.exp(-c * t / (2 * m)) * np.sin(wd * t) / (m * wd)
assert np.max(np.abs(h_num - h_formula)) < 1e-7                      # (3.25)
# FFT of h, zero-padded to refine the frequency grid (3.28); compare with |H| (3.23).
Npad = 2**16
xi_h = np.fft.rfftfreq(Npad, d=dt)
Hh = dt * np.fft.rfft(h_num, Npad)                                     # (3.29); record starts at 0
w_grid = 2 * np.pi * xi_h
sel = w_grid < 6
assert np.max(np.abs(np.abs(Hh[sel]) - np.abs(H(w_grid[sel])))) < 2e-3   # truncation at T=80 costs e^{-8}
i_pk = np.argmax(np.abs(Hh))
assert abs(w_grid[i_pk] - w_peak) < 0.01 and abs(np.abs(Hh[i_pk]) - H_max) / H_max < 1e-2
print(f"  FFT of impulse response: peak at w = {w_grid[i_pk]:.4f} (theory {w_peak:.4f}), height {np.abs(Hh[i_pk]):.4f} (theory {H_max:.4f})")

# Steady-state response to a slow chirp: envelope traces |H(w(t))|.
T_c = 600.0; wa, wb = 0.5, 3.5
rate = (wb - wa) / T_c
phase = lambda tt: wa * tt + 0.5 * rate * tt**2
w_inst = lambda tt: wa + rate * tt
tc = np.arange(0, T_c, 0.05)
solc = solve_ivp(lambda tt, y: [y[1], (np.cos(phase(tt)) - c * y[1] - k * y[0]) / m], (0, T_c), [0.0, 0.0],
                 t_eval=tc, rtol=1e-8, atol=1e-10, max_step=0.05, method="DOP853")
xc = solc.y[0]
# Envelope via the analytic signal: zero the negative frequencies, double the positive (a Hilbert transform in the DFT domain).
Xc = np.fft.fft(xc); fr = np.fft.fftfreq(len(xc), d=0.05)
Xc[fr < 0] = 0; Xc[fr > 0] *= 2
env = np.abs(np.fft.ifft(Xc))
i_env = np.argmax(env[(tc > 50) & (tc < T_c - 50)]) + np.sum(tc <= 50)
assert abs(w_inst(tc[i_env]) - w_peak) < 0.1, w_inst(tc[i_env])
assert abs(env[i_env] - H_max) / H_max < 0.1
print(f"  chirp: envelope peaks at w(t) = {w_inst(tc[i_env]):.3f} with amplitude {env[i_env]:.3f} (|H|max = {H_max:.3f})")

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes[0, 0].plot(t, h_num, label="solve_ivp"); axes[0, 0].plot(t, h_formula, "k--", lw=1, label="(3.25)")
axes[0, 0].set_xlim(0, 40); axes[0, 0].set_title("Impulse response h(t)"); axes[0, 0].legend(); axes[0, 0].grid(alpha=.3)
axes[0, 1].loglog(w_grid[1:][sel[1:]], np.abs(Hh[1:][sel[1:]]), label="|FFT of h| (zero-padded)")
axes[0, 1].loglog(w_grid[1:][sel[1:]], np.abs(H(w_grid[1:][sel[1:]])), "k--", label="|H(w)| (3.23)")
axes[0, 1].axvline(w_peak, color="r", lw=.8, ls=":"); axes[0, 1].set_title("Bode magnitude: flat 1/k, peak, 1/(m w^2) roll-off")
axes[0, 1].set_xlabel("w (rad/s)"); axes[0, 1].legend(); axes[0, 1].grid(alpha=.3, which="both")
axes[1, 0].plot(w_grid[sel], np.degrees(np.angle(H(w_grid[sel])))); axes[1, 0].set_title("Phase of H(w)")
axes[1, 0].set_xlabel("w (rad/s)"); axes[1, 0].set_ylabel("degrees"); axes[1, 0].grid(alpha=.3)
axes[1, 1].plot(w_inst(tc), xc, lw=.4, alpha=.6, label="x(t) driven by chirp")
axes[1, 1].plot(w_inst(tc), env, "C1", label="envelope")
axes[1, 1].plot(w_inst(tc), np.abs(H(w_inst(tc))), "k--", label="|H(w(t))|")
axes[1, 1].set_xlabel("instantaneous drive frequency w(t)"); axes[1, 1].set_title("Steady-state response to a slow chirp"); axes[1, 1].legend(); axes[1, 1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "oscillator.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: Aliasing (Example 4.7) and sinc interpolation (Theorem 3.28)

# %%
fs = 10.0; Ns_ = 50
tn = np.arange(Ns_) / fs
s3 = np.sin(2 * np.pi * 3 * tn); s7 = np.sin(2 * np.pi * 7 * tn); s12 = np.sin(2 * np.pi * 12 * tn)
assert np.allclose(s7, -s3)                                           # the samples are identical up to sign
k_peak = lambda s: np.argmax(np.abs(np.fft.rfft(s)))
assert k_peak(s3) == 15 and k_peak(s7) == 15 and k_peak(s12) == 10   # 3 Hz * 5 s = bin 15; 12 Hz aliases to 2 Hz = bin 10
fold = lambda nu, fs_: abs(nu - fs_ * np.round(nu / fs_))
assert fold(7, 10) == 3 and fold(12, 10) == 2 and fold(60, 50) == 10 and fold(60, 100) == 40 and fold(60, 130) == 60 and fold(9, 8) == 1
print("Part 7: 7 Hz sampled at 10 Hz == -(3 Hz samples); DFT peaks at bin 15 for both; 12 Hz -> bin 10 (2 Hz).")

fig, axes = plt.subplots(2, 2, figsize=(12, 7))
tf = np.linspace(0, 1, 2000)
axes[0, 0].plot(tf, np.sin(2 * np.pi * 7 * tf), "C0", lw=.8, label="7 Hz"); axes[0, 0].plot(tf, -np.sin(2 * np.pi * 3 * tf), "C1", lw=.8, label="-(3 Hz)")
axes[0, 0].plot(tn[tn < 1], s7[tn < 1], "ko", label="samples at 10 Hz"); axes[0, 0].legend(); axes[0, 0].set_title("Aliasing: the samples cannot tell 7 Hz from 3 Hz")
freqs = np.fft.rfftfreq(Ns_, d=1 / fs)
for s, lab in ((s3, "3 Hz"), (s7, "7 Hz"), (s12, "12 Hz")):
    axes[0, 1].plot(freqs, np.abs(np.fft.rfft(s)), "o-", ms=3, label=f"input {lab}")
axes[0, 1].set_xlabel("frequency (Hz) up to Nyquist = 5"); axes[0, 1].set_title("|DFT| of the samples"); axes[0, 1].legend(); axes[0, 1].grid(alpha=.3)


def sinc_interpolate(samples, dt_, t_eval, t_start=0.0):
    """Shannon reconstruction (3.34): f(t) = sum_n f(n dt) sinc((t - n dt)/dt)."""
    n = np.arange(len(samples))
    return np.array([np.sum(samples * np.sinc((tt - (t_start + n * dt_)) / dt_)) for tt in t_eval])


# A band-limited test signal: sinc^2(3t) has spectrum (1/3) Lambda(xi/3), support |xi| < 3;
# sinc(2(t-1.5)) has spectrum (1/2) Pi(xi/2), support |xi| < 1. Bandwidth B = 3 Hz.
f_bl = lambda tt: np.sinc(3 * tt) ** 2 + 0.6 * np.sinc(2 * (tt - 1.5))
t_fine = np.linspace(-5, 5, 2001)
for fs_, ok in ((10.0, True), (4.0, False)):
    dt_ = 1 / fs_; t_s = np.arange(-30, 30, dt_)
    rec = sinc_interpolate(f_bl(t_s), dt_, t_fine, t_start=t_s[0])
    err = np.max(np.abs(rec - f_bl(t_fine)))
    assert (err < 1e-2) == ok, (fs_, err)
    axes[1, 0 if ok else 1].plot(t_fine, f_bl(t_fine), "k", lw=2, label="band-limited f (B = 3 Hz)")
    axes[1, 0 if ok else 1].plot(t_fine, rec, "C1--", label=f"sinc interpolation, fs = {fs_:g} Hz")
    axes[1, 0 if ok else 1].plot(t_s[(t_s > -5) & (t_s < 5)], f_bl(t_s[(t_s > -5) & (t_s < 5)]), "ko", ms=3, label="samples")
    axes[1, 0 if ok else 1].set_title(f"fs = {fs_:g} Hz {'> 2B: exact' if ok else '< 2B: fails'}  (max err {err:.1e})"); axes[1, 0 if ok else 1].legend(fontsize=8); axes[1, 0 if ok else 1].grid(alpha=.3)
    print(f"  sinc interpolation at fs = {fs_:g} Hz (2B = 6): max error {err:.2e}")
fig.tight_layout(); fig.savefig(FIG / "sampling.png", dpi=120); plt.close(fig)
# C1(a): orthogonality of the shifted sincs, <s_n, s_m> = Delta delta_nm.
D = 0.5; xs = np.linspace(-400, 400, 1600001)
assert abs(np.trapezoid(np.sinc(xs / D) * np.sinc(xs / D - 1), xs)) < 1e-3 and abs(np.trapezoid(np.sinc(xs / D) ** 2, xs) - D) < 1e-3

# %% [markdown]
# ## Part 8: Poisson summation (3.20) and the comb (3.33)

# %%
nn = np.arange(-40, 41)
for tt in (0.5, 0.3, 2.0):
    lhs = np.sum(np.exp(-np.pi * nn**2 * tt)); rhs = tt**-0.5 * np.sum(np.exp(-np.pi * nn**2 / tt))
    assert abs(lhs - rhs) < 1e-13, (tt, lhs, rhs)
assert np.isclose(np.sum(np.exp(-np.pi * nn**2 * 0.5)), 1.419495, atol=1e-6)
print("Part 8: theta identity sum exp(-pi n^2 t) = t^(-1/2) sum exp(-pi k^2/t) holds to 1e-13; at t=1/2 both = 1.419495")
# Sampled Gaussian: Delta * DFT equals the periodized spectrum (3.33) exactly, and equals f^ itself when the copies do not overlap.
for Delta in (0.5, 0.125):
    xs_ = np.arange(-16, 16, Delta)
    xi_s, Fs = ft_numeric(np.exp(-np.pi * xs_**2), -16, Delta)
    periodized = sum(np.exp(-np.pi * (xi_s - kk / Delta) ** 2) for kk in range(-6, 7))
    assert np.max(np.abs(Fs - periodized)) < 1e-12, Delta
    alias_err = np.max(np.abs(Fs - np.exp(-np.pi * xi_s**2)))
    print(f"  Delta = {Delta}: Delta*DFT = periodized f^ to 1e-12; deviation from f^ itself (aliasing) = {alias_err:.1e}")
    if Delta == 0.125:
        assert alias_err < 1e-12   # copies 8 apart: e^{-pi 16} ~ 1e-22 overlap

# %% [markdown]
# ## Part 9: Spectrogram of a chirp with a click (STFT, §3.18), two window lengths

# %%
fs_s = 400.0; T_s = 8.0
ts = np.arange(0, T_s, 1 / fs_s)
f_start, f_end = 2.0, 40.0
sig = np.cos(2 * np.pi * (f_start * ts + 0.5 * (f_end - f_start) / T_s * ts**2))
sig[int(5.0 * fs_s)] += 20.0                                           # a click at t = 5 s


def hann(N):
    """Hann window w[n] = sin^2(pi n / N)."""
    return np.sin(np.pi * np.arange(N) / N) ** 2


def stft(x, Nw, hop):
    """Short-time Fourier transform: Hann-windowed DFTs of length Nw every `hop` samples.
    Returns (frame centers in samples, frequencies in cycles/sample, |S|)."""
    w = hann(Nw)
    starts = np.arange(0, len(x) - Nw, hop)
    S = np.array([np.abs(np.fft.rfft(x[s:s + Nw] * w)) for s in starts])
    return starts + Nw / 2, np.fft.rfftfreq(Nw), S.T


def fwhm_bins(col):
    """Full width at half maximum (in bins) of the main peak of a spectrum column."""
    i = np.argmax(col); half = col[i] / 2
    lo = i
    while lo > 0 and col[lo] > half: lo -= 1
    hi = i
    while hi < len(col) - 1 and col[hi] > half: hi += 1
    return hi - lo


hop = 8
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
widths_hz, click_frames = [], []
for ax, Nw in zip(axes, (32, 256)):
    centers, fr, S = stft(sig, Nw, hop)
    # frequency resolution at t = 2 s (chirp at 2 + 38/8*2 = 11.5 Hz)
    j = np.argmin(np.abs(centers / fs_s - 2.0))
    widths_hz.append(fwhm_bins(S[:, j]) * fs_s / Nw)
    # time smearing of the click: frames where the top-frequency band is above threshold
    hi_band = S[fr * fs_s > 150].max(axis=0)
    click_frames.append(np.sum(hi_band > 0.2 * hi_band.max()))
    ax.imshow(20 * np.log10(S + 1e-6), origin="lower", aspect="auto", cmap="magma", vmin=-20, vmax=40,
              extent=[centers[0] / fs_s, centers[-1] / fs_s, 0, fs_s / 2])
    ax.set_ylim(0, 60); ax.set_xlabel("time (s)"); ax.set_ylabel("frequency (Hz)")
    ax.set_title(f"Hann window {Nw} samples = {Nw/fs_s*1e3:.0f} ms: ridge FWHM {widths_hz[-1]:.1f} Hz, click smeared over {click_frames[-1]} frames")
fig.tight_layout(); fig.savefig(FIG / "spectrogram.png", dpi=120); plt.close(fig)
assert widths_hz[0] > 3 * widths_hz[1], widths_hz     # short window: poor frequency resolution
assert click_frames[1] > 3 * click_frames[0], click_frames  # long window: poor time resolution
print(f"Part 9: spectrogram ridge FWHM {widths_hz[0]:.1f} Hz (32-pt) vs {widths_hz[1]:.1f} Hz (256-pt); click smeared over {click_frames[0]} vs {click_frames[1]} frames.")

# %% [markdown]
# ## Your turn
#
# 1. **Inverse FFT via conjugation (B7a).** From F^{-1} = conj(F)/N: write
#    `ifft_via_conjugation(X)` that calls `fft_recursive` once. Assert against np.fft.ifft.
# 2. **Low-pass filter (B7c).** Build sin(2 pi 2 t) + 0.5 sin(2 pi 5 t) + noise at 100 Hz for
#    4 s; zero the DFT bins with |xi| > 3 Hz; inverse transform. Confirm the 5 Hz tone is gone.
#    Look at the two ends of the record and explain the ringing (Theorem 3.22: circular!).
# 3. **2D convolution theorem.** Make a 64x64 image (a bright square on a dark background),
#    blur it with a 2D Gaussian kernel two ways: direct 2D circular convolution (nested
#    loops or np.roll) and `ifft2(fft2(img) * fft2(kernel))`. Assert they agree. Plot the
#    image, its |fft2| (fftshifted, log scale), and the blurred result. Which spatial
#    frequencies did the blur remove?
# 4. **Windows (C3).** Implement `first_sidelobe_db` and compare rect vs Hann on the
#    2.5-cycle cosine of Example 4.5 at N = 64.
#
# Reference solutions for 1, 2 and 4 follow (B7 and C3); compare after you have tried.


# %%
def ifft_via_conjugation(X):
    """Inverse DFT from F^{-1} = conj(F)/N: f = conj(F conj(X))/N, eq. (3.27). (B7a reference.)"""
    X = np.asarray(X, complex)
    return np.conj(fft_recursive(np.conj(X))) / len(X)


Xr = np.fft.fft(rng.standard_normal(256) + 1j * rng.standard_normal(256))
assert np.max(np.abs(ifft_via_conjugation(Xr) - np.fft.ifft(Xr))) < 1e-10

# B7c: frequency-domain low-pass.
fs_lp, T_lp = 100.0, 4.0
t_lp = np.arange(0, T_lp, 1 / fs_lp)
x_lp = np.sin(2 * np.pi * 2 * t_lp) + 0.5 * np.sin(2 * np.pi * 5 * t_lp) + 0.3 * rng.standard_normal(t_lp.size)
X_lp = np.fft.fft(x_lp); xi_lp = np.fft.fftfreq(t_lp.size, d=1 / fs_lp)
Y_lp = X_lp.copy(); Y_lp[np.abs(xi_lp) > 3.0] = 0.0
y_lp = np.fft.ifft(Y_lp).real
k2, k5 = int(2 * T_lp), int(5 * T_lp)                                  # bins of 2 Hz and 5 Hz: k = nu * record length
assert abs(Y_lp[k5]) == 0 and np.isclose(abs(Y_lp[k2]), abs(X_lp[k2]))
# the 5 Hz tone is gone from y: project y onto sin(2 pi 5 t)
assert abs(np.sum(y_lp * np.sin(2 * np.pi * 5 * t_lp)) * 2 / t_lp.size) < 1e-12
assert abs(np.sum(y_lp * np.sin(2 * np.pi * 2 * t_lp)) * 2 / t_lp.size - 1.0) < 0.1
print("Your turn (B7): ifft via conjugation OK; low-pass removed the 5 Hz tone, kept 2 Hz.")


# C3: windows and leakage.
def first_sidelobe_db(A):
    """Level (dB) of the first sidelobe of a magnitude spectrum A with A[0] the main peak."""
    i = 1
    while A[i] < A[i - 1]: i += 1
    while A[i + 1] >= A[i]: i += 1
    return 20 * np.log10(A[i] / A[0]), i


Nw_, pad_ = 4096, 64
W_ = np.abs(np.fft.fft(hann(Nw_), pad_ * Nw_)); W_ /= W_.max()
R_ = np.abs(np.fft.fft(np.ones(Nw_), pad_ * Nw_)); R_ /= R_.max()
sl_rect, _ = first_sidelobe_db(R_); sl_hann, _ = first_sidelobe_db(W_)
assert abs(sl_rect - (-13.26)) < 0.05 and abs(sl_hann - (-31.47)) < 0.05
n64 = np.arange(64); tone = np.cos(2 * np.pi * 2.5 * n64 / 64)
leak_rect = np.abs(np.fft.fft(tone)); leak_hann = np.abs(np.fft.fft(tone * hann(64)))
far_rect = leak_rect[10:32].max() / leak_rect.max(); far_hann = leak_hann[10:32].max() / leak_hann.max()
bin30_rect = leak_rect[30] / leak_rect.max(); bin30_hann = leak_hann[30] / leak_hann.max()
assert far_rect > 0.05 and far_hann < 2e-3            # past bin 10: rect still leaks ~10% of the peak, Hann ~0.1%
assert bin30_rect > 0.04 and bin30_hann < 1e-4        # at bin 30: rect ~4.5% (1/|k| tails from +-2.5 add), Hann < 1e-5
print(f"Your turn (C3): first sidelobe rect {sl_rect:.2f} dB, Hann {sl_hann:.2f} dB; leakage beyond bin 10: "
      f"rect {far_rect:.1e}, Hann {far_hann:.1e}; at bin 30: rect {bin30_rect:.1e}, Hann {bin30_hann:.1e}")

# Your turn 3 (reference): the 2D convolution theorem, Section 3.17. A Gaussian blur is a low-pass filter.
img = np.zeros((64, 64)); img[24:40, 24:40] = 1.0                            # a bright square on a dark background
yy, xx = np.meshgrid(np.arange(64) - 32, np.arange(64) - 32, indexing="ij")
ker = np.exp(-(xx**2 + yy**2) / (2 * 2.0**2)); ker /= ker.sum()
ker = np.roll(ker, (-32, -32), axis=(0, 1))                                   # kernel centred at index (0, 0): no shift
blur_fft = np.fft.ifft2(np.fft.fft2(img) * np.fft.fft2(ker)).real            # (3.13) in 2D, circular
blur_direct = sum(ker[a, b] * np.roll(img, (a, b), axis=(0, 1))             # (ker (*) img)[i,j] = sum ker[a,b] img[i-a,j-b]
                  for a in range(64) for b in range(64))                     # all 4096 shifts: exact circular convolution
assert np.allclose(blur_fft, blur_direct, atol=1e-12)
K2 = np.fft.fft2(ker)
assert np.isclose(K2[0, 0].real, 1.0) and abs(K2[32, 32]) < 1e-3            # passes DC untouched, kills the Nyquist corner
assert np.isclose(blur_fft.sum(), img.sum())                                 # total brightness preserved (K2[0,0] = 1)
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
axes[0].imshow(img, cmap="gray"); axes[0].set_title("image")
axes[1].imshow(np.log10(np.abs(np.fft.fftshift(np.fft.fft2(img))) + 1e-3), cmap="magma"); axes[1].set_title("log|fft2| (fftshifted)")
axes[2].imshow(blur_fft, cmap="gray"); axes[2].set_title("Gaussian blur via fft2")
fig.tight_layout(); fig.savefig(FIG / "convolution_2d.png", dpi=120); plt.close(fig)
print("Your turn (2D): ifft2(fft2(img) fft2(ker)) == direct circular convolution; the blur's spectrum is 1 at DC and < 1e-3 at Nyquist.")

# %% [markdown]
# ## Numbers quoted in the notes and solutions that are not checked above
# (Examples 4.2, 4.3, 4.8; §3.9; solutions A1, A3, A5d, B1c, B4b, B5b, B6b, C3)

# %%
assert np.isclose(1 / (2 * np.pi), 0.159, atol=5e-4)                                                    # Example 4.2: sigma_xi = 0.159 Hz for sigma = 1 s
assert np.isclose(sinc(0.2), 0.93549, atol=5e-6)                                                        # Example 4.3
assert np.isclose(np.sqrt(2) - 1, 0.41, atol=5e-3)                                                      # §3.9 / B2(d): exponential is 41% above the minimum
assert np.isclose(np.exp(-9 * np.pi / 2), 7e-7, atol=5e-8)                                              # Example 4.8: third term 0.0000007
assert np.isclose(np.degrees(np.arctan(0.6 / 5)), 6.843, atol=5e-4)                                     # B6(b)
xi_sl = np.linspace(1.0, 2.0, 200001); sl_sinc = np.abs(sinc(xi_sl)).max()                              # B1(c): first sidelobe of sinc
assert np.isclose(sl_sinc, 0.2172, atol=5e-5) and np.isclose(20 * np.log10(sl_sinc), -13.26, atol=0.01)  # "about 22%", -13 dB
first_zero = lambda A: int(np.argmax(np.diff(A) > 0))                                                   # fine-grid index where the main lobe stops descending
assert abs(first_zero(R_) - pad_) <= 1 and abs(first_zero(W_) - 2 * pad_) <= 1                          # C3: main lobe 1 bin (rect) vs 2 bins (Hann)
assert np.isclose(1 / np.sqrt(10), 0.3162, atol=5e-5)                                                   # A1
assert np.isclose(2 / np.sqrt(2 * np.pi), 0.7979, atol=5e-5)                                            # A3(c)
assert fold(13, 12) == 1 and 13 - 12 * np.round(13 / 12) == 1                                           # A5(d): +1 rev/s, forward
assert 2**19 * 20 == 10_485_760 and 2**40 == 1_099_511_627_776 and 2**40 / (2**19 * 20) == 104857.6    # B4(b)
assert np.isclose(2 * np.exp(-np.pi / 2), 0.415759, atol=5e-7) and np.isclose(2 * np.exp(-2 * np.pi), 0.003735, atol=5e-7)   # B5(b)
assert np.isclose(np.exp(-np.pi / 2), 0.2078796, atol=5e-8) and np.isclose(np.exp(-2 * np.pi), 0.0018674, atol=5e-8)        # Example 4.8
assert np.isclose(np.sqrt(2) * (1 + 2 * np.exp(-2 * np.pi)), 1.419495, atol=5e-7)
print("Remaining quoted numbers (Examples 4.2/4.3/4.8, §3.9, A1, A3, A5d, B1c, B4b, B5b, B6b, C3 lobe widths) verified.")

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.perf_counter() - T0:.1f} s. Figures in", FIG)
