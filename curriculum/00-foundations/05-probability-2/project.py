# %% [markdown]
# # F5 Project: Joint Distributions, the Gaussian, Limit Theorems, Chains, Entropy
#
# This script builds the objects of F5 from scratch and checks every number quoted
# in the notes, the worked examples and the solutions. Parts: (1) discrete joint
# tables: marginalize, condition, moments, Adam and Eve; (2) Box–Muller, a QQ plot
# against scipy, polar coordinates and the Cauchy ratio; (3) sums and order
# statistics: convolution, Gamma from Exponentials, Beta from Uniforms, the
# Beta–Binomial posterior; (4) covariance algebra and two counterexamples;
# (5) the multivariate normal: Cholesky sampling, contours with scaled
# eigenvectors, Mahalanobis distance, the conditional slice against the Schur
# complement; (6) conditional expectation by binning, the orthogonal residual,
# Adam/Eve on a mixture, Galton's heights; (7) Chebyshev, CLT histograms, the
# error-versus-n plot, continuity correction, a law of large numbers for a Markov
# chain; (8) Markov chains: stationary vector, mixing slope, detailed balance,
# random walk on a graph, lazy walk, Metropolis, PageRank; (9) entropy, KL,
# cross-entropy, mutual information, maximum entropy.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import math
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gamma as Gamma_fn, comb
import sympy as sp

T_START = time.time()
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
LN2 = np.log(2.0)

# %% [markdown]
# ## Part 1: Discrete joint distributions (Example 4.1, exercises A1 and B9)
#
# A joint PMF is a table p[i, j] = P(X = x_i, Y = y_j). Marginals (3.1) are row and
# column sums; conditionals (3.3) are renormalized rows; independence is
# Theorem 3.4(a), joint == outer product of marginals. E[Y | X = x_i] is the mean of
# row i after renormalizing (Definition 3.25).
#
# **Build it yourself first:** write `marginal`, `conditional`, `joint_moments`
# before reading these.


# %%
def marginal(p, axis):
    """Marginal PMF, eq. (3.1): sum the joint table over the OTHER variable.
    axis=0 returns p_X (sum over columns j); axis=1 returns p_Y (sum over rows i)."""
    p = np.asarray(p, float)
    return p.sum(axis=1) if axis == 0 else p.sum(axis=0)


def conditional(p, given, index):
    """Conditional PMF, eq. (3.3). given='X', index=i: p_{Y|X}(. | x_i) = p[i, :] / p_X(x_i).
    given='Y', index=j: p_{X|Y}(. | y_j) = p[:, j] / p_Y(y_j)."""
    p = np.asarray(p, float)
    row = p[index, :] if given == "X" else p[:, index]
    return row / row.sum()


def is_independent(p, tol=1e-12):
    """Theorem 3.4(a): X and Y independent iff p == outer(p_X, p_Y)."""
    p = np.asarray(p, float)
    return bool(np.allclose(p, np.outer(marginal(p, 0), marginal(p, 1)), atol=tol))


def joint_moments(p, xs, ys):
    """Means, variances, E[XY], Cov (3.13) and rho (Definition 3.16) of a joint PMF table."""
    p = np.asarray(p, float); xs = np.asarray(xs, float); ys = np.asarray(ys, float)
    pX, pY = marginal(p, 0), marginal(p, 1)
    EX, EY = pX @ xs, pY @ ys
    VX, VY = pX @ xs**2 - EX**2, pY @ ys**2 - EY**2
    EXY = float(np.sum(p * np.outer(xs, ys)))
    cov = EXY - EX * EY
    return dict(EX=EX, EY=EY, VarX=VX, VarY=VY, EXY=EXY, Cov=cov, rho=cov / np.sqrt(VX * VY))


def cond_expectation_rows(p, ys):
    """g(x_i) = E[Y | X = x_i] for every row i (Definition 3.25), as an array."""
    p = np.asarray(p, float); ys = np.asarray(ys, float)
    return (p @ ys) / marginal(p, 0)


def cond_variance_rows(p, ys):
    """Var(Y | X = x_i) for every row i (Definition 3.27)."""
    p = np.asarray(p, float); ys = np.asarray(ys, float)
    pX = marginal(p, 0)
    return (p @ ys**2) / pX - cond_expectation_rows(p, ys) ** 2


# Example 4.1
P41 = np.array([[0.10, 0.20], [0.20, 0.20], [0.20, 0.10]])
xs41, ys41 = np.array([0, 1, 2]), np.array([0, 1])
assert np.allclose(marginal(P41, 0), [0.3, 0.4, 0.3]) and np.allclose(marginal(P41, 1), [0.5, 0.5])
assert np.allclose(conditional(P41, "X", 0), [1 / 3, 2 / 3])
assert np.allclose(conditional(P41, "X", 1), [0.5, 0.5])
assert np.allclose(conditional(P41, "X", 2), [2 / 3, 1 / 3])
assert not is_independent(P41) and not np.isclose(P41[0, 0], 0.3 * 0.5)
m = joint_moments(P41, xs41, ys41)
assert np.isclose(m["EX"], 1.0) and np.isclose(m["VarX"], 0.6) and np.isclose(m["EY"], 0.5) and np.isclose(m["VarY"], 0.25)
assert np.isclose(m["EXY"], 0.40) and np.isclose(m["Cov"], -0.10) and np.isclose(m["rho"], -0.258, atol=5e-4)
g41 = cond_expectation_rows(P41, ys41)
assert np.allclose(g41, [2 / 3, 1 / 2, 1 / 3])
pX41 = marginal(P41, 0)
assert np.isclose(pX41 @ g41, 0.5)                                   # Adam's law, Theorem 3.26
EcondVar = pX41 @ cond_variance_rows(P41, ys41)
VarCondMean = pX41 @ g41**2 - (pX41 @ g41) ** 2
assert np.isclose(EcondVar, 0.2333, atol=5e-5) and np.isclose(VarCondMean, 0.0167, atol=5e-5)
assert np.isclose(EcondVar + VarCondMean, m["VarY"])                 # Eve's law, eq. (3.21)
assert np.isclose(m["Cov"] / m["VarX"], -1 / 6)                      # best linear slope, eq. (3.23)
assert np.isclose(VarCondMean / m["VarY"], 0.0667, atol=5e-4)        # explained fraction 6.7%
print(f"Example 4.1: Cov = {m['Cov']:.3f}, rho = {m['rho']:.4f}, E[Y|X=x] = {np.round(g41, 4)}, "
      f"E[Var(Y|X)] = {EcondVar:.4f}, Var(E[Y|X]) = {VarCondMean:.4f}, sum = {EcondVar + VarCondMean:.4f}")

# Exercise A1: uncorrelated but dependent table
PA1 = np.array([[0.10, 0.20, 0.10], [0.20, 0.20, 0.20]])
xsA1, ysA1 = np.array([1, 2]), np.array([1, 2, 3])
assert np.allclose(marginal(PA1, 0), [0.4, 0.6]) and np.allclose(marginal(PA1, 1), [0.3, 0.4, 0.3])
assert np.allclose(conditional(PA1, "X", 0), [0.25, 0.5, 0.25]) and np.allclose(conditional(PA1, "X", 1), [1 / 3] * 3)
assert not is_independent(PA1)
mA1 = joint_moments(PA1, xsA1, ysA1)
assert np.isclose(mA1["EX"], 1.6) and np.isclose(mA1["EY"], 2.0) and np.isclose(mA1["EXY"], 3.2) and np.isclose(mA1["Cov"], 0.0)
assert np.allclose(cond_expectation_rows(PA1, ysA1), [2.0, 2.0])
assert np.allclose(cond_variance_rows(PA1, ysA1), [0.5, 2 / 3])
assert np.isclose(marginal(PA1, 0) @ cond_variance_rows(PA1, ysA1), 0.6) and np.isclose(mA1["VarY"], 0.6)
print("Exercise A1: Cov = 0 exactly, E[Y|X] constant = 2, yet the conditionals differ (dependent).")

# Exact Eve's law on Example 4.1 with SymPy rationals (no floating point at all).
R = sp.Rational
Pq = sp.Matrix([[R(1, 10), R(2, 10)], [R(2, 10), R(2, 10)], [R(2, 10), R(1, 10)]])
pXq = [sum(Pq.row(i)) for i in range(3)]
gq = [Pq[i, 1] / pXq[i] for i in range(3)]
E_cond_var_q = sum(pXq[i] * gq[i] * (1 - gq[i]) for i in range(3))
var_cond_mean_q = sum(pXq[i] * gq[i] ** 2 for i in range(3)) - R(1, 4)
assert E_cond_var_q == R(7, 30) and var_cond_mean_q == R(1, 60) and E_cond_var_q + var_cond_mean_q == R(1, 4)
print("SymPy: E[Var(Y|X)] = 7/30, Var(E[Y|X]) = 1/60, sum = 1/4 = Var(Y) exactly.")

# %% [markdown]
# ## Part 2: Box–Muller, QQ plot, polar coordinates, the Cauchy ratio (Example 4.2, B1)
#
# Eq. (3.6): R = sqrt(-2 ln U1), Theta = 2 pi U2, Z = (R cos Theta, R sin Theta).
# The derivation in section 3.3 used the Jacobian of the inverse map (s, theta) <-
# (z1, z2), whose determinant is 2; SymPy confirms it. Reading the map backwards
# (exercise B1): for iid standard normals, R is Rayleigh, Theta is Uniform, they
# are independent, and Z1/Z2 is Cauchy.
#
# **Build it yourself first:** `box_muller(u1, u2)`.


# %%
def box_muller(u1, u2):
    """Eq. (3.6): two independent Uniform(0,1) arrays -> two independent N(0,1) arrays."""
    r = np.sqrt(-2.0 * np.log(u1))
    th = 2.0 * np.pi * u2
    return r * np.cos(th), r * np.sin(th)


# SymPy: Jacobian of the inverse map (old (s, theta) with respect to new (z1, z2)) has det 2.
z1s, z2s = sp.symbols("z1 z2", real=True)
S_expr, TH_expr = z1s**2 + z2s**2, sp.atan2(z2s, z1s)
Jac = sp.Matrix([[sp.diff(S_expr, z1s), sp.diff(S_expr, z2s)], [sp.diff(TH_expr, z1s), sp.diff(TH_expr, z2s)]])
assert sp.simplify(Jac.det()) == 2
print("SymPy: det of the Box–Muller inverse Jacobian =", sp.simplify(Jac.det()))

# Example 4.2 numbers
a1, a2 = box_muller(np.array(0.5), np.array(0.25))
assert np.isclose(a1, 0.0, atol=1e-12) and np.isclose(a2, 1.1774, atol=5e-5)
b1, b2 = box_muller(np.array(0.1), np.array(0.125))
assert np.isclose(b1, 1.5174, atol=5e-5) and np.isclose(b2, 1.5174, atol=5e-5)
assert np.isclose(b1**2 + b2**2, -2 * np.log(0.1)) and np.isclose(b1**2 + b2**2, 4.605, atol=5e-4)
print(f"Example 4.2: (u1,u2)=(0.5,0.25) -> ({a1:.4f}, {a2:.4f});  (0.1,0.125) -> ({b1:.4f}, {b2:.4f}); "
      f"z1^2+z2^2 = {b1**2 + b2**2:.4f} = -2 ln 0.1")

N = 100_000
Z1, Z2 = box_muller(rng.random(N), rng.random(N))
Z = np.concatenate([Z1, Z2])
assert abs(Z.mean()) < 0.01 and abs(Z.var() - 1) < 0.02
assert abs(np.corrcoef(Z1, Z2)[0, 1]) < 0.01, "Z1 and Z2 should be uncorrelated"
ks = stats.kstest(Z, "norm")
assert ks.pvalue > 1e-3, ks
print(f"Box–Muller: {N} pairs; mean {Z.mean():+.4f}, var {Z.var():.4f}, corr(Z1,Z2) = {np.corrcoef(Z1, Z2)[0, 1]:+.4f}, "
      f"KS p = {ks.pvalue:.3f}")

# QQ plot against scipy's Normal quantiles
Zs = np.sort(Z1)
q_theory = stats.norm.ppf((np.arange(1, N + 1) - 0.5) / N)
fig, (axq, axh) = plt.subplots(1, 2, figsize=(10, 4.3))
axq.plot(q_theory, Zs, ".", ms=2, alpha=.5, label="Box–Muller sample quantiles")
axq.plot([-4.5, 4.5], [-4.5, 4.5], "k--", lw=1, label="slope 1 through origin")
axq.set_xlabel("scipy.stats.norm.ppf"); axq.set_ylabel("sample quantile"); axq.set_title("QQ plot: Box–Muller vs N(0,1)")
axq.legend(fontsize=8); axq.grid(alpha=.3)
axh.hist(Z, bins=80, density=True, alpha=.6, label="Box–Muller")
xx = np.linspace(-4.5, 4.5, 400); axh.plot(xx, stats.norm.pdf(xx), "k", lw=1.5, label="N(0,1) density")
axh.set_title("Histogram"); axh.legend(fontsize=8); axh.grid(alpha=.3)
fig.tight_layout(); fig.savefig(FIG / "box_muller_qq.png", dpi=120); plt.close(fig)
# the sample quantiles should lie on the line: slope and intercept of a fit through the central 99%
mask = np.abs(q_theory) < 2.6
slope_qq, icpt_qq = np.polyfit(q_theory[mask], Zs[mask], 1)
assert abs(slope_qq - 1) < 0.01 and abs(icpt_qq) < 0.01, (slope_qq, icpt_qq)

# Exercise B1: polar coordinates of the 2D standard normal and the Cauchy ratio
Rr = np.hypot(Z1, Z2); Th = np.arctan2(Z2, Z1)
assert stats.kstest(Rr, "rayleigh").pvalue > 1e-3                      # R has density r e^{-r^2/2}
assert abs(Rr.mean() - np.sqrt(np.pi / 2)) < 0.01                       # E[R] = sqrt(pi/2) = 1.2533
assert stats.kstest(Rr**2, "expon", args=(0, 2)).pvalue > 1e-3          # R^2 ~ Exponential(1/2) = chi^2_2
assert stats.kstest(Th, "uniform", args=(-np.pi, 2 * np.pi)).pvalue > 1e-3
assert abs(np.corrcoef(Rr, Th)[0, 1]) < 0.01
Tc = Z1 / Z2
q25, q50, q75 = np.quantile(Tc, [0.25, 0.5, 0.75])
assert abs(q50) < 0.02 and abs(q25 + 1) < 0.03 and abs(q75 - 1) < 0.03   # Cauchy quartiles are -1, 0, +1
assert stats.kstest(Tc, "cauchy").pvalue > 1e-3
print(f"B1: E[R] = {Rr.mean():.4f} (sqrt(pi/2) = {np.sqrt(np.pi/2):.4f}); Z1/Z2 quartiles {q25:+.3f}, {q50:+.3f}, {q75:+.3f} "
      f"(Cauchy: -1, 0, 1); KS p vs cauchy = {stats.kstest(Tc, 'cauchy').pvalue:.3f}")

# %% [markdown]
# ## Part 3: Sums and order statistics (Example 4.3, A3, Theorems 3.8, 3.9, 3.12, B4)
#
# Convolution (3.7) of two Uniform(0,1) densities is the triangle of Example 4.3.
# Poissons add rates. n Exponentials add to Gamma(n, lambda) (Theorem 3.8); note
# scipy parametrizes Gamma by *scale* = 1/rate. The k-th of n Uniforms is
# Beta(k, n-k+1) (Theorem 3.9), and the Beta–Binomial posterior (3.12) is checked
# with the numbers of exercise B4.

# %%
dx = 1e-3
xg = np.arange(0.0, 1.0 + dx / 2, dx)
f_unif = np.ones_like(xg)
f_S = np.convolve(f_unif, f_unif) * dx                 # density of U1 + U2 on the grid s = k dx
s_grid = np.arange(len(f_S)) * dx
tri = np.where(s_grid <= 1, s_grid, 2 - s_grid)
assert np.abs(f_S - tri).max() < 2e-3, "numerical convolution should be the triangle"
P_half = np.trapezoid(f_S[s_grid <= 0.5], s_grid[s_grid <= 0.5])
assert np.isclose(P_half, 0.125, atol=2e-3)
S_mc = rng.random(200_000) + rng.random(200_000)
assert abs(S_mc.mean() - 1) < 5e-3 and abs(S_mc.var() - 1 / 6) < 3e-3 and abs(np.mean(S_mc <= 0.5) - 0.125) < 3e-3
print(f"Example 4.3: f_S(0.5) = {f_S[np.argmin(np.abs(s_grid - 0.5))]:.3f}, f_S(1) = {f_S[np.argmin(np.abs(s_grid - 1))]:.3f}, "
      f"P(S <= 1/2) = {P_half:.4f} (exact 0.125); MC var = {S_mc.var():.4f} (1/6 = {1/6:.4f})")

# A3(c): Poisson(2) + Poisson(3) = Poisson(5) by discrete convolution
k = np.arange(0, 40)
pmf_sum = np.convolve(stats.poisson.pmf(k, 2), stats.poisson.pmf(k, 3))[: len(k)]
assert np.allclose(pmf_sum, stats.poisson.pmf(k, 5), atol=1e-12)
assert np.isclose(pmf_sum[4], 0.1755, atol=5e-5)
print(f"A3(c): P(Pois(2) + Pois(3) = 4) = {pmf_sum[4]:.4f} = e^-5 5^4/4! = {np.exp(-5) * 625 / 24:.4f}")

# Theorem 3.8: five Exponential(rate 2) add to Gamma(shape 5, rate 2); scipy wants scale = 1/2.
E5 = rng.exponential(scale=0.5, size=(50_000, 5)).sum(axis=1)
assert stats.kstest(E5, "gamma", args=(5, 0, 0.5)).pvalue > 1e-3
assert abs(E5.mean() - 2.5) < 0.02 and abs(E5.var() - 1.25) < 0.05  # mean a/lambda, var a/lambda^2
# chi^2_3 = Gamma(3/2, 1/2) is the sum of three squared normals (section 3.5)
W3 = (rng.standard_normal((50_000, 3)) ** 2).sum(axis=1)
assert stats.kstest(W3, "chi2", args=(3,)).pvalue > 1e-3 and abs(W3.mean() - 3) < 0.03
print(f"Theorem 3.8: sum of 5 Exp(2): mean {E5.mean():.3f} (2.5), var {E5.var():.3f} (1.25), KS p = "
      f"{stats.kstest(E5, 'gamma', args=(5, 0, 0.5)).pvalue:.3f}")

# A3(a,b): min of n Exponentials is Exponential(n lambda); E[max of two Exp(1)] = 1.5
Emin = rng.exponential(scale=0.5, size=(100_000, 5)).min(axis=1)
assert stats.kstest(Emin, "expon", args=(0, 0.1)).pvalue > 1e-3 and abs(Emin.mean() - 0.1) < 2e-3
Emax2 = rng.exponential(size=(200_000, 2)).max(axis=1)
assert abs(Emax2.mean() - 1.5) < 0.01
print(f"A3: min of 5 Exp(2) has mean {Emin.mean():.4f} (1/10); max of 2 Exp(1) has mean {Emax2.mean():.4f} (3/2)")

# Theorem 3.9: 2nd smallest of 5 Uniforms is Beta(2, 4); max of 5 is Beta(5, 1)
U5 = np.sort(rng.random((100_000, 5)), axis=1)
assert stats.kstest(U5[:, 1], "beta", args=(2, 4)).pvalue > 1e-3
assert stats.kstest(U5[:, 4], "beta", args=(5, 1)).pvalue > 1e-3
assert abs(U5[:, 1].mean() - 2 / 6) < 3e-3 and abs(U5[:, 4].mean() - 5 / 6) < 3e-3
# Theorem 3.11: B(a, b) = Gamma(a)Gamma(b)/Gamma(a+b), checked by quadrature at (2.5, 1.5)
B_quad = quad(lambda t: t**1.5 * (1 - t) ** 0.5, 0, 1)[0]
assert np.isclose(B_quad, Gamma_fn(2.5) * Gamma_fn(1.5) / Gamma_fn(4.0))
print(f"Theorem 3.9/3.11: U_(2) of 5 mean {U5[:, 1].mean():.4f} (1/3); B(2.5,1.5) = {B_quad:.6f} by quad = Gamma formula")

# Theorem 3.12 and exercise B4: Beta(1,1) prior, 7 successes in 10 -> Beta(8, 4)
a0, b0, n_tr, k_s = 1, 1, 10, 7
post = stats.beta(a0 + k_s, b0 + n_tr - k_s)
w_prior = (a0 + b0) / (a0 + b0 + n_tr)
assert np.isclose(post.mean(), 2 / 3) and np.isclose(post.mean(), w_prior * a0 / (a0 + b0) + (1 - w_prior) * k_s / n_tr)
assert np.isclose(post.std(), 0.1307, atol=5e-5) and np.isclose(post.sf(0.5), 0.8867, atol=5e-5)
# Beta(8, 4) is the 8th smallest of 11 Uniforms (Theorem 3.9): P(p > 1/2) = P(Bin(11, 1/2) <= 7) = 227/256
assert np.isclose(post.sf(0.5), stats.binom.cdf(7, 11, 0.5)) and np.isclose(post.sf(0.5), 227 / 256)
assert sum(comb(11, j, exact=True) for j in range(8, 12)) == 232
print(f"B4: posterior Beta(8,4): mean {post.mean():.4f}, sd {post.std():.4f}, P(p > 1/2 | data) = {post.sf(0.5):.4f}; "
      f"weights {w_prior:.4f} on prior mean, {1 - w_prior:.4f} on 7/10")

# %% [markdown]
# ## Part 4: Covariance algebra and two counterexamples (A2, section 3.7, B5)
#
# Var(a^T X) = a^T Sigma a and Cov(a^T X, b^T X) = a^T Sigma b, eq. (3.14). Then the
# two things covariance cannot see: X and X^2 are uncorrelated but dependent, and
# (X, SX) with a random sign S has Gaussian marginals, zero correlation, and is
# neither independent nor jointly Gaussian (exercise B5).

# %%
SigA2 = np.array([[4.0, -3.0], [-3.0, 9.0]])
cov_lin = lambda a, b: np.asarray(a) @ SigA2 @ np.asarray(b)
assert np.isclose(cov_lin([1, 1], [1, 1]), 7) and np.isclose(cov_lin([1, -1], [1, -1]), 19)
assert np.isclose(cov_lin([2, -3], [2, -3]), 133) and np.isclose(cov_lin([2, 0], [-1, 3]), -26)
assert np.isclose(-3 / np.sqrt(4 * 9), -0.5)
evA2 = np.linalg.eigvalsh(SigA2)
assert np.allclose(evA2, [2.595, 10.405], atol=5e-4) and evA2.min() > 0
print(f"A2: Var(X+Y) = 7, Var(X-Y) = 19, Var(2X-3Y) = 133, Cov(2X+1, 3Y-X) = -26, rho = -0.5, eigenvalues {np.round(evA2, 3)}")

# Uncorrelated but dependent, section 3.7: X uniform on {-1, 0, 1}, Y = X^2.
PXX2 = np.zeros((3, 2)); xv, yv = np.array([-1, 0, 1]), np.array([0, 1])
PXX2[0, 1] = PXX2[2, 1] = 1 / 3; PXX2[1, 0] = 1 / 3
mm = joint_moments(PXX2, xv, yv)
assert np.isclose(mm["Cov"], 0) and not is_independent(PXX2)
assert np.isclose(PXX2[1, 0], 1 / 3) and np.isclose(marginal(PXX2, 0)[1] * marginal(PXX2, 1)[0], 1 / 9)
print("Section 3.7: X uniform on {-1,0,1}, Y = X^2: Cov = 0 but P(0,0) = 1/3 != 1/9. Dependent.")

# Exercise B5: Gaussian marginals, zero correlation, not independent, not jointly Gaussian.
Xb = rng.standard_normal(200_000); Sb = rng.choice([-1.0, 1.0], size=200_000); Yb = Sb * Xb
assert stats.kstest(Yb, "norm").pvalue > 1e-3                 # Y ~ N(0,1)
assert abs(np.corrcoef(Xb, Yb)[0, 1]) < 0.01                   # uncorrelated
frac_zero = np.mean(Xb + Yb == 0.0)
assert abs(frac_zero - 0.5) < 0.01                             # X + Y has an atom at 0: not Gaussian
assert np.mean((np.abs(Xb) < 1) & (np.abs(Yb) > 1)) == 0.0     # |Y| = |X|: not independent
print(f"B5: Y = S X: KS-normal p = {stats.kstest(Yb, 'norm').pvalue:.3f}, corr = {np.corrcoef(Xb, Yb)[0, 1]:+.4f}, "
      f"P(X + Y = 0) = {frac_zero:.4f} (1/2), P(|X|<1, |Y|>1) = 0.")

# %% [markdown]
# ## Part 5: The multivariate normal (Example 4.4, A4; Theorems 3.21–3.24)
#
# Sampling: X = mu + L Z with L L^T = Sigma (Cholesky). Level sets are ellipses with
# axes along the eigenvectors of Sigma and semi-axes sqrt(c lambda_i). The
# Mahalanobis distance squared is chi^2_2 = Exponential(1/2). The conditional
# X1 | X2 = x2 is N(mu1 + S12 S22^-1 (x2 - mu2), S11 - S12 S22^-1 S21), eq. (3.19),
# checked against a thin slice of samples at x2 ~ 1.
#
# **Build it yourself first:** `sample_mvn` and `mvn_conditional`.


# %%
def sample_mvn(mu, Sigma, n, rng):
    """n draws of N(mu, Sigma) as rows: mu + L Z with L the Cholesky factor (Definition 3.20, 'Sampling' in 3.8)."""
    L = np.linalg.cholesky(np.asarray(Sigma, float))
    Zn = rng.standard_normal((n, len(mu)))
    return np.asarray(mu, float) + Zn @ L.T, L


def mvn_conditional(mu, Sigma, idx1, idx2, x2):
    """Theorem 3.24, eq. (3.19): mean and covariance of X[idx1] given X[idx2] = x2."""
    mu = np.asarray(mu, float); Sigma = np.asarray(Sigma, float); x2 = np.atleast_1d(np.asarray(x2, float))
    S11 = Sigma[np.ix_(idx1, idx1)]; S12 = Sigma[np.ix_(idx1, idx2)]; S22 = Sigma[np.ix_(idx2, idx2)]
    C = S12 @ np.linalg.inv(S22)                     # regression coefficient / Kalman gain
    return mu[idx1] + C @ (x2 - mu[idx2]), S11 - C @ S12.T


def mahalanobis2(X, mu, Sigma):
    """d(x)^2 = (x - mu)^T Sigma^-1 (x - mu) = ||L^-1 (x - mu)||^2 for each row of X."""
    L = np.linalg.cholesky(np.asarray(Sigma, float))
    Zt = np.linalg.solve(L, (np.asarray(X, float) - mu).T)
    return (Zt**2).sum(axis=0)


mu44 = np.zeros(2); Sig44 = np.array([[2.0, 1.0], [1.0, 2.0]])
lam, Q = np.linalg.eigh(Sig44)                       # ascending: lam = (1, 3)
assert np.allclose(lam, [1, 3])
q_long, q_short = Q[:, 1] * np.sign(Q[0, 1]), Q[:, 0] * np.sign(Q[0, 0])
assert np.allclose(q_long, [1, 1] / np.sqrt(2)) and np.allclose(np.abs(q_short), [1, 1] / np.sqrt(2)) and np.isclose(q_short @ [1, -1], np.sqrt(2) * np.sign(q_short[0]))
assert np.isclose(np.linalg.det(Sig44), 3) and np.allclose(np.linalg.inv(Sig44), np.array([[2, -1], [-1, 2]]) / 3)
peak = stats.multivariate_normal(mu44, Sig44).pdf(mu44)
assert np.isclose(peak, 1 / (2 * np.pi * np.sqrt(3))) and np.isclose(peak, 0.0919, atol=5e-5)
rho44 = Sig44[0, 1] / np.sqrt(Sig44[0, 0] * Sig44[1, 1]); assert np.isclose(rho44, 0.5)
c95 = stats.chi2.ppf(0.95, 2)
assert np.isclose(c95, -2 * np.log(0.05)) and np.isclose(c95, 5.991, atol=5e-4)
assert np.isclose(np.sqrt(c95 * 3), 4.24, atol=5e-3) and np.isclose(np.sqrt(c95), 2.45, atol=5e-3)
d2_a, d2_b = mahalanobis2(np.array([[1, 1], [1, -1.0]]), mu44, Sig44)
assert np.isclose(d2_a, 2 / 3) and np.isclose(d2_b, 2.0)
assert np.isclose(np.exp(-1), 0.368, atol=5e-4) and np.isclose(np.exp(-1 / 3), 0.717, atol=5e-4)
cm, cv = mvn_conditional(mu44, Sig44, [0], [1], 1.0)
assert np.isclose(cm[0], 0.5) and np.isclose(cv[0, 0], 1.5) and np.isclose(cv[0, 0], 2 * (1 - rho44**2))
X44, L44 = sample_mvn(mu44, Sig44, 400_000, rng)
assert np.allclose(L44, [[np.sqrt(2), 0], [1 / np.sqrt(2), np.sqrt(1.5)]]) and np.allclose(L44 @ L44.T, Sig44)
assert np.allclose(L44, [[1.4142, 0], [0.7071, 1.2247]], atol=5e-5)
Shat = np.cov(X44.T)
assert np.allclose(Shat, Sig44, atol=0.02), Shat                         # sample covariance ~ Sigma
d2 = mahalanobis2(X44, mu44, Sig44)
assert abs(np.mean(d2 <= c95) - 0.95) < 3e-3                            # 95% ellipse holds 95%
assert abs(np.mean(d2 > 2) - np.exp(-1)) < 3e-3 and abs(np.mean(d2 > 2 / 3) - np.exp(-1 / 3)) < 3e-3
assert stats.kstest(d2, "chi2", args=(2,)).pvalue > 1e-3                # d^2 ~ chi^2_2
# the conditional slice at x2 ~ 1 versus the Schur-complement formula
slab = X44[np.abs(X44[:, 1] - 1.0) < 0.05, 0]
assert len(slab) > 3000
assert abs(slab.mean() - 0.5) < 0.05 and abs(slab.var() - 1.5) < 0.1, (slab.mean(), slab.var())
print(f"Example 4.4: eigenvalues {lam[::-1]}, L = {np.round(L44, 4).tolist()}, peak density {peak:.4f}, "
      f"95% semi-axes {np.sqrt(c95*3):.3f}/{np.sqrt(c95):.3f}; d^2(1,1) = {d2_a:.4f}, d^2(1,-1) = {d2_b:.4f}")
print(f"   sample covariance {np.round(Shat, 3).tolist()}; inside 95% ellipse: {np.mean(d2 <= c95):.4f}; "
      f"slice |x2 - 1| < 0.05 ({len(slab)} pts): mean {slab.mean():.3f} (0.5), var {slab.var():.3f} (1.5)")
# binned conditional means lie on the regression line x1 = 0.5 x2 (slope Sigma12/Sigma22)
edges = np.linspace(-3, 3, 25); centers = 0.5 * (edges[1:] + edges[:-1])
idx = np.digitize(X44[:, 1], edges) - 1
bin_means = np.array([X44[idx == j, 0].mean() for j in range(len(centers))])
slope_bins = np.polyfit(centers, bin_means, 1)[0]
assert abs(slope_bins - 0.5) < 0.02

# Contours with scaled eigenvectors and the 95% ellipse
gx = np.linspace(-5, 5, 300); GX, GY = np.meshgrid(gx, gx)
dens = stats.multivariate_normal(mu44, Sig44).pdf(np.dstack([GX, GY]))
fig, ax = plt.subplots(figsize=(5.8, 5.8))
ax.scatter(X44[:3000, 0], X44[:3000, 1], s=3, alpha=.25, color="gray", label="3000 Cholesky samples")
ax.contour(GX, GY, dens, levels=8, cmap="viridis")
tt = np.linspace(0, 2 * np.pi, 300)
ell = (Q * np.sqrt(c95 * lam)) @ np.vstack([np.cos(tt), np.sin(tt)])
ax.plot(ell[0], ell[1], "r-", lw=1.5, label="95% ellipse: d^2 = 5.991")
for lam_i, q_i, c in zip(lam, Q.T, ["C1", "C3"]):
    v = np.sqrt(c95 * lam_i) * q_i
    ax.annotate("", xy=v, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=c, lw=2))
    ax.text(*(v * 1.1), f"sqrt(5.991 x {lam_i:g}) q", color=c, fontsize=8)
# horizontal axis is x1, vertical is x2: plot the points (E[X1 | X2 = c], c), which lie on x1 = 0.5 x2
ax.plot(bin_means, centers, "ks", ms=3, label="binned E[X1 | X2]: the line x1 = 0.5 x2")
ax.axhline(1.0, color="k", lw=.6, ls=":"); ax.text(-4.8, 1.1, "slice x2 = 1", fontsize=8)
ax.plot(cm[0], 1.0, "r*", ms=12, label="E[X1 | X2 = 1] = 0.5, eq. (3.19)")
ax.set_aspect("equal"); ax.set_xlim(-5, 5); ax.set_ylim(-5, 5); ax.grid(alpha=.3); ax.legend(fontsize=7, loc="lower right")
ax.set_title("N(0, [[2,1],[1,2]]): axes = eigenvectors, semi-axes ∝ sqrt(λ)")
fig.tight_layout(); fig.savefig(FIG / "mvn_contours.png", dpi=120); plt.close(fig)

# Exercise A4: mu = (1, 2), Sigma = [[4, -2], [-2, 4]]
muA4 = np.array([1.0, 2.0]); SigA4 = np.array([[4.0, -2.0], [-2.0, 4.0]])
assert np.allclose(np.linalg.eigvalsh(SigA4), [2, 6]) and np.isclose(-2 / 4, -0.5)
cmA4, cvA4 = mvn_conditional(muA4, SigA4, [0], [1], 4.0)
assert np.isclose(cmA4[0], 0.0) and np.isclose(cvA4[0, 0], 3.0)
B_pm = np.array([[1, 1], [1, -1.0]])
cov_pm = B_pm @ SigA4 @ B_pm.T                        # Theorem 3.22
assert np.allclose(B_pm @ muA4, [3, -1]) and np.allclose(cov_pm, [[4, 0], [0, 12]])
print(f"A4: X1 | X2 = 4 ~ N({cmA4[0]:.1f}, {cvA4[0,0]:.1f}); (X1+X2, X1-X2) ~ N(({3},{-1}), diag(4, 12)): independent by Theorem 3.23")

# (3.19) on a partition where the blocks differ in size and scale (Example 4.4 and A4 have equal
# diagonal blocks, so inverting the wrong block would go unnoticed there). Independent route: the
# precision matrix Lam = Sigma^-1 gives X1 | X2 = x2 ~ N(mu1 - Lam11^-1 Lam12 (x2 - mu2), Lam11^-1)
# (Bishop PRML 2.3.1), which is (3.20) read through the block inverse of step 3.
mu3 = np.array([1.0, -2.0, 0.5])
Sig3 = np.array([[4.0, 1.2, -0.8], [1.2, 2.0, 0.3], [-0.8, 0.3, 0.5]])
x2_3 = np.array([-1.0, 1.0])
cm3, cv3 = mvn_conditional(mu3, Sig3, [0], [1, 2], x2_3)
Lam3 = np.linalg.inv(Sig3)
cv3_prec = np.linalg.inv(Lam3[:1, :1])
cm3_prec = mu3[:1] - cv3_prec @ Lam3[:1, 1:] @ (x2_3 - mu3[1:])
assert np.allclose(cm3, cm3_prec) and np.allclose(cv3, cv3_prec)
assert np.isclose(np.linalg.det(Sig3), np.linalg.det(cv3) * np.linalg.det(Sig3[1:, 1:]))    # det Sigma = det S det Sigma22
assert np.all(np.linalg.eigvalsh(Sig3[:1, :1] - cv3) >= -1e-12)                             # conditioning shrinks: Sigma11 - S is PSD
cm3b, cv3b = mvn_conditional(mu3, Sig3, [1, 2], [0], 3.0)                                    # the other direction, 2 | 1
cv3b_prec = np.linalg.inv(Lam3[1:, 1:])
assert np.allclose(cv3b, cv3b_prec) and np.allclose(cm3b, mu3[1:] - cv3b_prec @ Lam3[1:, :1] @ (np.array([3.0]) - mu3[:1]))
print(f"(3.19) on a 1+2 partition of a 3x3 Sigma: mean {cm3[0]:.4f}, var {cv3[0, 0]:.4f}; matches the precision-matrix form")

# %% [markdown]
# ## Part 6: Conditional expectation (section 3.9; Example 4.5; exercises B2, B3)
#
# E[Y | X] estimated by binning X. Theorem 3.29: the residual Y - E[Y|X] is
# orthogonal to every function of X, while the residual of the best *line* (3.23)
# is orthogonal only to linear functions. Adam's and Eve's laws on a Gaussian
# mixture (B3), the vector best-linear-predictor formula (B2), and Galton's heights
# (Example 4.5) are simulated.


# %%
def binned_conditional_mean(x, y, edges):
    """Nonparametric estimate of g(x) = E[Y | X = x]: average y inside each bin of x.
    Returns bin centers, bin means, and g_hat evaluated at every sample point."""
    idx = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
    sums = np.bincount(idx, weights=y, minlength=len(edges) - 1)
    counts = np.bincount(idx, minlength=len(edges) - 1)
    means = sums / np.maximum(counts, 1)
    return 0.5 * (edges[1:] + edges[:-1]), means, means[idx]


Nce = 200_000
Xc = rng.random(Nce)                                    # Uniform(0, 1)
Yc = Xc**2 + 0.1 * rng.standard_normal(Nce)              # E[Y | X] = X^2, Var(Y | X) = 0.01
edges_c = np.linspace(0, 1, 41)
cen, gbar, g_hat = binned_conditional_mean(Xc, Yc, edges_c)
assert np.abs(gbar - cen**2).max() < 0.01                # binned means recover x^2
resid_ce = Yc - g_hat
# best linear predictor: beta = Cov/Var = (1/4 - 1/6)/(1/12) = 1, alpha = 1/3 - 1/2 = -1/6
beta_hat = np.cov(Xc, Yc)[0, 1] / Xc.var(ddof=1); alpha_hat = Yc.mean() - beta_hat * Xc.mean()
assert abs(beta_hat - 1.0) < 0.01 and abs(alpha_hat + 1 / 6) < 0.01
resid_lin = Yc - (alpha_hat + beta_hat * Xc)
corr = lambda u, v: np.corrcoef(u, v)[0, 1]
for h_name, h in [("X", Xc), ("X^2", Xc**2), ("sin(5X)", np.sin(5 * Xc)), ("exp(X)", np.exp(Xc))]:
    assert abs(corr(resid_ce, h)) < 0.02, f"residual of E[Y|X] correlated with {h_name}"
assert abs(corr(resid_lin, Xc)) < 0.01                   # the line's residual is uncorrelated with X ...
assert corr(resid_lin, Xc**2) > 0.1                      # ... but not with X^2 (0.149 in theory)
mse_ce, mse_lin = np.mean(resid_ce**2), np.mean(resid_lin**2)
assert mse_ce < mse_lin and abs(mse_ce - 0.01) < 1e-3 and abs(mse_lin - 0.0156) < 1e-3
assert abs(Yc.var() - (0.01 + 4 / 45)) < 2e-3            # Eve: 0.01 + Var(X^2) = 0.01 + (1/5 - 1/9)
print(f"Part 6: beta = {beta_hat:.4f} (1), alpha = {alpha_hat:.4f} (-1/6); corr(resid_E[Y|X], X^2) = {corr(resid_ce, Xc**2):+.4f}, "
      f"corr(resid_line, X^2) = {corr(resid_lin, Xc**2):+.4f}; MSE {mse_ce:.4f} vs {mse_lin:.4f}")

fig, ax = plt.subplots(figsize=(6, 4.2))
ax.scatter(Xc[:4000], Yc[:4000], s=3, alpha=.3, color="gray", label="samples")
ax.plot(cen, gbar, "C3o-", ms=4, label="binned E[Y|X]")
ax.plot(cen, cen**2, "k--", lw=1, label="true x^2")
ax.plot(cen, alpha_hat + beta_hat * cen, "C0-", lw=1.5, label="best line y = x - 1/6")
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.legend(fontsize=8); ax.grid(alpha=.3)
ax.set_title("E[Y|X] is the projection; the line is the projection onto lines")
fig.tight_layout(); fig.savefig(FIG / "conditional_expectation.png", dpi=120); plt.close(fig)

# Exercise B2: vector best linear predictor equals Sigma12 Sigma22^-1 (sample version), non-Gaussian data
X2v = rng.random((Nce, 2)) * 2 - 1                                 # Uniform cube, not Gaussian
X1v = np.c_[X2v @ [1.0, -2.0] + 0.5 * X2v[:, 0] ** 2, X2v @ [0.5, 0.5]] + 0.2 * rng.standard_normal((Nce, 2))
Sfull = np.cov(np.c_[X1v, X2v].T)
S12, S22 = Sfull[:2, 2:], Sfull[2:, 2:]
B_hat = S12 @ np.linalg.inv(S22)
Xd = np.c_[np.ones(Nce), X2v]
B_ls = np.linalg.lstsq(Xd, X1v, rcond=None)[0][1:].T              # least squares gives the same B
assert np.allclose(B_hat, B_ls, atol=1e-6)
res_v = X1v - X1v.mean(0) - (X2v - X2v.mean(0)) @ B_hat.T
assert np.abs(np.cov(res_v.T, X2v.T)[:2, 2:]).max() < 1e-3           # residual uncorrelated with X2
assert np.allclose(np.cov(res_v.T), Sfull[:2, :2] - S12 @ np.linalg.inv(S22) @ S12.T, atol=1e-6)  # Schur complement
print(f"B2: B = Sigma12 Sigma22^-1 = {np.round(B_hat, 3).tolist()} equals least squares; residual cov = Schur complement.")

# Exercise B3: Adam and Eve on a Gaussian mixture Z ~ Bern(0.3), Y|Z=0 ~ N(0,1), Y|Z=1 ~ N(4, 2^2)
pz, m0, s0, m1, s1 = 0.3, 0.0, 1.0, 4.0, 2.0
Zm = rng.random(400_000) < pz
Ym = np.where(Zm, m1 + s1 * rng.standard_normal(400_000), m0 + s0 * rng.standard_normal(400_000))
EY_th = (1 - pz) * m0 + pz * m1; EV_th = (1 - pz) * s0**2 + pz * s1**2; VE_th = pz * (1 - pz) * (m1 - m0) ** 2
assert np.isclose(EY_th, 1.2) and np.isclose(EV_th, 1.9) and np.isclose(VE_th, 3.36) and np.isclose(EV_th + VE_th, 5.26)
assert np.isclose(0.5 * np.log(2 * np.pi * np.e * 5.26), 2.249, atol=5e-4)     # Gaussian entropy at the mixture's variance
assert abs(Ym.mean() - 1.2) < 0.01 and abs(Ym.var() - 5.26) < 0.05
within = (1 - pz) * Ym[~Zm].var() + pz * Ym[Zm].var()
between = (1 - pz) * (Ym[~Zm].mean() - Ym.mean()) ** 2 + pz * (Ym[Zm].mean() - Ym.mean()) ** 2
assert abs(within - 1.9) < 0.03 and abs(between - 3.36) < 0.05
P_gt3 = (1 - pz) * stats.norm.sf(3) + pz * stats.norm.sf((3 - m1) / s1)
post_z = pz * stats.norm.pdf(3, m1, s1) / ((1 - pz) * stats.norm.pdf(3, m0, s0) + pz * stats.norm.pdf(3, m1, s1))
assert np.isclose(P_gt3, 0.2084, atol=5e-5) and np.isclose(post_z, 0.9445, atol=5e-5)
assert abs(np.mean(Ym > 3) - P_gt3) < 3e-3
# B3(e): third central moment sum_z P(z) [3 d s^2 + d^3], d = m_z - 1.2; a normal has 0
mu3 = (1 - pz) * (3 * (m0 - 1.2) * s0**2 + (m0 - 1.2) ** 3) + pz * (3 * (m1 - 1.2) * s1**2 + (m1 - 1.2) ** 3)
assert np.isclose(mu3, 12.936) and abs(np.mean((Ym - Ym.mean()) ** 3) - 12.936) < 0.4
print(f"B3: E[Y] = {Ym.mean():.3f} (1.2); within {within:.3f} (1.9) + between {between:.3f} (3.36) = {within + between:.3f} (5.26); "
      f"P(Y > 3) = {P_gt3:.4f}, P(Z = 1 | Y = 3) = {post_z:.4f}, third central moment {mu3:.3f} (MC {np.mean((Ym - Ym.mean()) ** 3):.3f})")

# Example 4.5: Galton's heights, rho = 0.5, mu = 68, sigma = 2.5
mu_h, sd_h, rho_h = 68.0, 2.5, 0.5
Nf = 100_000
Pp = mu_h + sd_h * rng.standard_normal(Nf)
Cc = mu_h + rho_h * (Pp - mu_h) + sd_h * np.sqrt(1 - rho_h**2) * rng.standard_normal(Nf)
slope_g = np.cov(Pp, Cc)[0, 1] / Pp.var(ddof=1); icpt_g = Cc.mean() - slope_g * Pp.mean()
assert abs(slope_g - 0.5) < 0.01 and abs(icpt_g - 34.0) < 0.7
kids73 = Cc[np.abs(Pp - 73) < 0.25]
assert len(kids73) > 500 and abs(kids73.mean() - 70.5) < 0.25 and abs(np.mean(kids73 > 73) - 0.124) < 0.04
assert np.isclose(sd_h**2 * (1 - rho_h**2), 4.6875) and np.isclose(np.sqrt(4.6875), 2.165, atol=5e-4)
assert np.isclose(stats.norm.sf((73 - 70.5) / np.sqrt(4.6875)), 0.124, atol=5e-4)
assert abs(np.var(Cc - (icpt_g + slope_g * Pp)) - 4.6875) < 0.05                # E[Var(C|P)]
assert abs(np.var(icpt_g + slope_g * Pp) - 1.5625) < 0.03                       # Var(E[C|P]) = rho^2 sigma^2
assert abs(Cc.var() - 6.25) < 0.06
slope_rev = np.cov(Pp, Cc)[0, 1] / Cc.var(ddof=1)
assert abs(slope_rev - 0.5) < 0.01                                               # regression is symmetric
print(f"Example 4.5: fitted C = {icpt_g:.2f} + {slope_g:.3f} P; children of 73-inch mid-parents: mean {kids73.mean():.2f} (70.5), "
      f"P(C > 73) = {np.mean(kids73 > 73):.3f} (0.124); reverse slope {slope_rev:.3f}")
fig, ax = plt.subplots(figsize=(5.5, 5.5))
ax.scatter(Pp[:5000], Cc[:5000], s=3, alpha=.3, color="gray")
hh = np.linspace(60, 76, 2)
ax.plot(hh, icpt_g + slope_g * hh, "C3", lw=2, label="E[C | P]: slope 0.5")
ax.plot(Pp.mean() + slope_rev * (hh - Cc.mean()), hh, "C0", lw=2, label="E[P | C]: slope 0.5 the other way")
ax.plot(hh, hh, "k:", lw=1, label="C = P")
ax.plot([73, 73], [60, 76], "k--", lw=.8); ax.plot(73, 70.5, "ko", ms=6)
ax.set_xlabel("mid-parent height (in)"); ax.set_ylabel("child height (in)"); ax.set_xlim(60, 76); ax.set_ylim(60, 76)
ax.set_aspect("equal"); ax.grid(alpha=.3); ax.legend(fontsize=8); ax.set_title("Regression to the mean is symmetric")
fig.tight_layout(); fig.savefig(FIG / "galton.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 7: Inequalities and limit theorems (sections 3.10–3.11; exercises B6, B7)
#
# Chebyshev is crude but assumption-free. The CLT histograms for Uniform,
# Exponential, Bernoulli(0.1) and an infinite-variance Pareto at n = 1, 2, 5, 30
# show what converges and what does not. The sup-norm error |F_n - Phi| is
# computed exactly (Binomial, Gamma, Irwin–Hall) and decays like n^(-1/2) for the
# skewed distributions (Berry–Esseen, (3.25)) but like n^(-1) for the symmetric
# Uniform, whose skewness term in the Edgeworth expansion vanishes.

# %%
assert np.isclose(2 * stats.norm.sf(2), 0.046, atol=1e-3) and 2 * stats.norm.sf(2) < 0.25      # Gaussian vs Chebyshev at 2 sigma
Ej = rng.exponential(size=200_000)
assert np.mean(np.log(Ej)) < np.log(Ej.mean()) and abs(np.mean(np.log(Ej)) + np.euler_gamma) < 0.01  # Jensen: E log X <= log E X
print(f"Chebyshev at 2 sigma gives 0.25; the Gaussian truth is {2 * stats.norm.sf(2):.4f}. Jensen: E[ln X] = {np.mean(np.log(Ej)):.4f} "
      f"(-gamma = {-np.euler_gamma:.4f}) <= ln E[X] = 0 for X ~ Exp(1).")

samplers = {
    "Uniform(0,1)": (lambda size: rng.random(size), 0.5, np.sqrt(1 / 12)),
    "Exponential(1)": (lambda size: rng.exponential(size=size), 1.0, 1.0),
    "Bernoulli(0.1)": (lambda size: (rng.random(size) < 0.1).astype(float), 0.1, 0.3),
    "Pareto(alpha=1.5): infinite variance": (lambda size: rng.random(size) ** (-1 / 1.5), 3.0, None),
}
ns = [1, 2, 5, 30]
fig, axes = plt.subplots(4, 4, figsize=(14, 11))
zz = np.linspace(-4, 4, 300)
for r_, (name, (samp, mean_, sd_)) in enumerate(samplers.items()):
    for c_, n in enumerate(ns):
        Xn = samp((100_000, n)).sum(axis=1)
        if sd_ is None:                                   # infinite variance: standardize by the sample
            Sn = (Xn - Xn.mean()) / Xn.std()
        else:
            Sn = (Xn - n * mean_) / (sd_ * np.sqrt(n))
        ax = axes[r_, c_]
        if name.startswith("Bernoulli"):                  # lattice: one bin per integer value of the sum,
            edges_b = (np.arange(-0.5, n + 1.5) - n * mean_) / (sd_ * np.sqrt(n))   # so bar heights compare with phi
            ax.hist(Sn, bins=edges_b, density=True, alpha=.7)
        else:
            ax.hist(np.clip(Sn, -4, 4), bins=80, range=(-4, 4), density=True, alpha=.7)
        ax.plot(zz, stats.norm.pdf(zz), "k", lw=1.2)
        ks_n = stats.kstest(Sn, "norm").statistic
        ax.set_title(f"{name}, n = {n}, KS dist {ks_n:.3f}", fontsize=9); ax.set_xlim(-4, 4); ax.set_yticks([])
        if name.startswith("Pareto") and n == 30:
            assert ks_n > 0.08, "the Pareto sum must NOT look normal"
        if name.startswith("Uniform") and n == 30:
            assert ks_n < 0.01
fig.suptitle("Standardized sums (3.24): bell curve appears for finite variance, never for Pareto(1.5)")
fig.tight_layout(); fig.savefig(FIG / "clt_histograms.png", dpi=110); plt.close(fig)


def ks_binomial_vs_normal(n, p):
    """sup_x |P(S_n <= x) - Phi((x - np)/sqrt(np(1-p)))|, exact; the sup sits at the jumps."""
    kk = np.arange(0, n + 1)
    F = stats.binom.cdf(kk, n, p); Fprev = np.r_[0.0, F[:-1]]
    Ph = stats.norm.cdf((kk - n * p) / np.sqrt(n * p * (1 - p)))
    return max(np.abs(F - Ph).max(), np.abs(Fprev - Ph).max())


def ks_gamma_vs_normal(n):
    """sup_x |Gamma(n,1) CDF - Phi((x - n)/sqrt n)| on a fine grid (Exponential(1) sums)."""
    x = np.linspace(max(0, n - 10 * np.sqrt(n)), n + 10 * np.sqrt(n), 40_001)
    return np.abs(stats.gamma.cdf(x, n) - stats.norm.cdf((x - n) / np.sqrt(n))).max()


def irwin_hall_cdf(x, n):
    """CDF of the sum of n iid Uniform(0,1): (1/n!) sum_{k<=floor(x)} (-1)^k C(n,k) (x-k)^n."""
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    for kk in range(n + 1):
        term = (-1) ** kk * comb(n, kk, exact=True) * np.clip(x - kk, 0, None) ** n
        out += term
    return np.clip(out / math.factorial(n), 0, 1)


def ks_irwin_hall_vs_normal(n):
    x = np.linspace(0, n, 20_001)
    return np.abs(irwin_hall_cdf(x, n) - stats.norm.cdf((x - n / 2) / np.sqrt(n / 12))).max()


n_grid = np.array([1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000])
err_bern = np.array([ks_binomial_vs_normal(n, 0.1) for n in n_grid])
err_exp = np.array([ks_gamma_vs_normal(n) for n in n_grid])
n_unif = np.array([1, 2, 3, 5, 8, 12, 16, 20])
err_unif = np.array([ks_irwin_hall_vs_normal(n) for n in n_unif])
n_par = np.array([1, 2, 5, 10, 30, 100, 300, 1000])
err_par = []
for n in n_par:
    Xn = (rng.random((20_000, n)) ** (-1 / 1.5)).sum(axis=1)
    err_par.append(stats.kstest((Xn - Xn.mean()) / Xn.std(), "norm").statistic)
err_par = np.array(err_par)
fit = lambda nn, ee, lo: np.polyfit(np.log(nn[nn >= lo]), np.log(ee[nn >= lo]), 1)[0]
s_bern, s_exp, s_unif = fit(n_grid, err_bern, 10), fit(n_grid, err_exp, 10), fit(n_unif, err_unif, 3)
assert -0.7 < s_bern < -0.3 and -0.6 < s_exp < -0.4, (s_bern, s_exp)      # Berry–Esseen rate n^-1/2
assert -1.25 < s_unif < -0.8, s_unif                                       # symmetric: n^-1
assert err_par.min() > 0.05, "Pareto(1.5): the error never goes to zero"
assert np.isclose(0.4748 * 2.7333 / np.sqrt(100), 0.1298, atol=1e-3)         # Berry–Esseen bound for Bern(0.1), n = 100
print(f"Error vs n (log-log slopes): Bernoulli(0.1) {s_bern:.2f}, Exponential {s_exp:.2f} (Berry–Esseen: -0.5); "
      f"Uniform {s_unif:.2f} (symmetric: -1); Pareto(1.5) min error {err_par.min():.3f} (no convergence)")
fig, ax = plt.subplots(figsize=(6.2, 4.4))
ax.loglog(n_grid, err_bern, "o-", label=f"Bernoulli(0.1), exact; slope {s_bern:.2f}")
ax.loglog(n_grid, err_exp, "s-", label=f"Exponential(1), exact; slope {s_exp:.2f}")
ax.loglog(n_unif, err_unif, "^-", label=f"Uniform(0,1), exact; slope {s_unif:.2f}")
ax.loglog(n_par, err_par, "x--", label="Pareto(1.5), Monte Carlo")
ax.loglog(n_grid, 0.4748 * 2.7333 / np.sqrt(n_grid), "k:", label="Berry–Esseen bound, Bernoulli(0.1)")
ax.set_xlabel("n"); ax.set_ylabel("sup |F_n - Φ|"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7)
ax.set_title("CLT error: n^-1/2 if skewed, n^-1 if symmetric, never for infinite variance")
fig.tight_layout(); fig.savefig(FIG / "clt_error_vs_n.png", dpi=120); plt.close(fig)

# Normal approximation to the Binomial with continuity correction (section 3.11)
exact45 = stats.binom.cdf(45, 100, 0.5); no_cc = stats.norm.cdf((45 - 50) / 5); cc = stats.norm.cdf((45.5 - 50) / 5)
assert np.isclose(exact45, 0.1841, atol=5e-5) and np.isclose(no_cc, 0.1587, atol=5e-5) and np.isclose(cc, 0.1841, atol=5e-5)
print(f"Bin(100, 1/2): P(X <= 45) exact {exact45:.4f}; normal without correction {no_cc:.4f}; with +1/2 correction {cc:.4f}")

# Exercise B6: Chebyshev vs CLT vs exact for 100 fair coins, |Xbar - 1/2| >= 0.1
cheb = 0.0025 / 0.01; clt_b6 = 2 * stats.norm.sf(2); exact_b6 = stats.binom.cdf(40, 100, 0.5) + stats.binom.sf(59, 100, 0.5)
cc_b6 = 2 * stats.norm.cdf((40.5 - 50) / 5)
assert np.isclose(cheb, 0.25) and np.isclose(clt_b6, 0.0455, atol=5e-5) and np.isclose(exact_b6, 0.0569, atol=5e-5) and np.isclose(cc_b6, 0.0574, atol=5e-5)
assert 0.25 / (0.01 * 0.05) == 500
n_clt = int(np.ceil((stats.norm.ppf(0.975) / 0.2) ** 2))                    # 2(1 - Phi(0.2 sqrt n)) <= 0.05
assert n_clt == 97 and np.isclose((1.96 / 0.2) ** 2, 96.04)
assert exact_b6 < 0.25 / 4                                                  # Chebyshev is off by a factor > 4
print(f"B6: Chebyshev {cheb:.4f}; CLT {clt_b6:.4f}; CLT with correction {cc_b6:.4f}; exact {exact_b6:.4f}; "
      f"for 0.05 Chebyshev needs n >= 500, the CLT n >= {n_clt}")
# B6(d): Berry–Esseen for Bernoulli(0.1): rho_3 / sigma^3 = (0.9 * 0.1^3 + 0.1 * 0.9^3) / 0.3^3
ratio_be = (0.9 * 0.1**3 + 0.1 * 0.9**3) / 0.3**3
err100 = ks_binomial_vs_normal(100, 0.1)
assert np.isclose(ratio_be, 2.7333, atol=5e-5) and np.isclose(0.4748 * ratio_be / 10, 0.1298, atol=5e-5)
assert np.isclose(err100, 0.0832, atol=5e-5) and err100 < 0.4748 * ratio_be / 10
print(f"B6(d): Berry–Esseen bound for Bern(0.1), n = 100: {0.4748 * ratio_be / 10:.4f}; actual sup error {err100:.4f}")

# Exercise B7: law of large numbers for a two-state Markov chain with a = 0.2, b = 0.3 (lambda_2 = 0.5)
a_c, b_c = 0.2, 0.3
P2 = np.array([[1 - a_c, a_c], [b_c, 1 - b_c]])
pi2 = np.array([b_c, a_c]) / (a_c + b_c)                   # (0.6, 0.4)
lam2 = 1 - a_c - b_c
Tn = 200_000
cum2 = P2.cumsum(axis=1); x_ch = np.empty(Tn + 1, int); x_ch[0] = int(rng.random() < pi2[1]); u_ch = rng.random(Tn)
for t in range(Tn):
    x_ch[t + 1] = int(u_ch[t] > cum2[x_ch[t], 0])
xc = x_ch.astype(float)
var_x = pi2[0] * pi2[1]
for kk in range(1, 6):
    acov = np.mean((xc[:-kk] - pi2[1]) * (xc[kk:] - pi2[1]))
    assert abs(acov - lam2**kk * var_x) < 0.01, (kk, acov, lam2**kk * var_x)   # Cov(X_0, X_k) = lambda_2^k Var
running = np.cumsum(xc) / np.arange(1, Tn + 2)
assert abs(running[-1] - pi2[1]) < 0.01
var_bound = var_x / Tn * (1 + lam2) / (1 - lam2)
assert np.isclose(var_x, 0.24) and np.isclose(var_x * (1 + lam2) / (1 - lam2), 0.72)          # Var(Xbar_n) <= 0.72 / n
print(f"B7: two-state chain, lambda_2 = {lam2}: Cov(X_0, X_k)/Var = {[round(float(np.mean((xc[:-kk] - .4) * (xc[kk:] - .4)) / var_x), 3) for kk in range(1, 4)]} "
      f"(0.5^k); running mean {running[-1]:.4f} -> pi_1 = 0.4; Var(Xbar) <= {var_bound:.2e}")

# %% [markdown]
# ## Part 8: Markov chains (section 3.12; Example 4.6; exercises B8, B10)
#
# Row-vector convention: pi_{t+1} = pi_t P. The stationary distribution is the left
# eigenvector for eigenvalue 1, found as an eigenvector of P^T (3.26). Mixing is
# geometric with ratio |lambda_2|. Detailed balance (3.27) holds for the weather
# chain and the random walk on a graph, fails for a circulating chain, and holds
# by construction for Metropolis.
#
# **Build it yourself first:** `stationary(P)` and `simulate_chain`.


# %%
def stationary(P):
    """Stationary row vector pi with pi P = pi, eq. (3.26): eigenvector of P^T for eigenvalue 1, normalized."""
    w, V = np.linalg.eig(np.asarray(P, float).T)
    i = np.argmin(np.abs(w - 1))
    pi = np.real(V[:, i])
    return pi / pi.sum()


def simulate_chain(P, x0, T, rng):
    """Sample X_0 = x0, X_1, ..., X_T from the chain with transition matrix P (rows sum to 1)."""
    cum = np.asarray(P, float).cumsum(axis=1)
    x = np.empty(T + 1, int); x[0] = x0
    u = rng.random(T)
    for t in range(T):
        x[t + 1] = int(np.searchsorted(cum[x[t]], u[t]))
    return x


def detailed_balance_gap(P, pi):
    """max_ij |pi_i P_ij - pi_j P_ji|: zero iff (3.27) holds."""
    F = pi[:, None] * np.asarray(P, float)
    return float(np.abs(F - F.T).max())


# Example 4.6: the weather chain
Pw = np.array([[0.7, 0.3, 0.0], [0.3, 0.4, 0.3], [0.0, 0.6, 0.4]])
assert np.allclose(Pw.sum(axis=1), 1)
pi_w = stationary(Pw)
assert np.allclose(pi_w, [0.4, 0.4, 0.2]) and np.allclose(pi_w @ Pw, pi_w)
assert detailed_balance_gap(Pw, pi_w) < 1e-12                                # a path graph: reversible
ev_w = np.sort(np.linalg.eigvals(Pw).real)[::-1]
assert np.allclose(ev_w, [1, 0.5854, -0.0854], atol=5e-5) and np.isclose(np.trace(Pw), 1.5) and np.isclose(np.linalg.det(Pw), -0.05)
pi1 = np.array([1.0, 0, 0]) @ Pw; pi2_ = pi1 @ Pw
assert np.allclose(pi1, [0.7, 0.3, 0]) and np.allclose(pi2_, [0.58, 0.33, 0.09])
d1, d2_ = np.abs(pi1 - pi_w).sum(), np.abs(pi2_ - pi_w).sum()
assert np.isclose(d1, 0.6) and np.isclose(d2_, 0.36) and np.isclose(d2_ / d1, 0.6)
dists = []; pt = np.array([1.0, 0, 0])
for t in range(1, 31):
    pt = pt @ Pw; dists.append(np.abs(pt - pi_w).sum())
dists = np.array(dists); tt_ = np.arange(1, 31)
slope_mix = np.polyfit(tt_[4:25], np.log(dists[4:25]), 1)[0]
assert abs(slope_mix - np.log(ev_w[1])) < 0.01 and np.isclose(np.log(ev_w[1]), -0.535, atol=5e-4)
xw = simulate_chain(Pw, 0, 100_000, rng)
freq_w = np.bincount(xw, minlength=3) / len(xw)
assert np.abs(freq_w - pi_w).max() < 0.01                                    # ergodic theorem
print(f"Example 4.6: pi = {np.round(pi_w, 4)}, eigenvalues {np.round(ev_w, 4)}, ||pi_1 - pi|| = {d1:.2f}, ||pi_2 - pi|| = {d2_:.2f}; "
      f"mixing slope {slope_mix:.4f} vs ln(lambda_2) = {np.log(ev_w[1]):.4f}; simulated frequencies {np.round(freq_w, 4)}")

# A circulating chain: stationary (uniform) but NOT in detailed balance
Pcirc = np.array([[0.1, 0.8, 0.1], [0.1, 0.1, 0.8], [0.8, 0.1, 0.1]])
pi_c = stationary(Pcirc)
assert np.allclose(pi_c, 1 / 3) and detailed_balance_gap(Pcirc, pi_c) > 0.2
assert np.isclose(pi_c[0] * Pcirc[0, 1], 0.8 / 3) and np.isclose(pi_c[1] * Pcirc[1, 0], 0.1 / 3)
print(f"Circulating chain: pi uniform, but pi_1 P_12 = {pi_c[0]*Pcirc[0,1]:.4f} != pi_2 P_21 = {pi_c[1]*Pcirc[1,0]:.4f}: not reversible")

# Random walk on a graph: stationary proportional to degree
edges_g = [(0, 1), (0, 2), (1, 2), (2, 3), (3, 4), (2, 4)]      # contains triangles: aperiodic
Ag = np.zeros((5, 5))
for i, j in edges_g:
    Ag[i, j] = Ag[j, i] = 1
deg = Ag.sum(axis=1); Pg = Ag / deg[:, None]
pi_g = stationary(Pg)
assert np.allclose(pi_g, deg / deg.sum()) and np.allclose(deg / deg.sum(), [1 / 6, 1 / 6, 1 / 3, 1 / 6, 1 / 6])
assert detailed_balance_gap(Pg, pi_g) < 1e-12
xg_ = simulate_chain(Pg, 0, 100_000, rng)
assert np.abs(np.bincount(xg_, minlength=5) / len(xg_) - pi_g).max() < 0.01
print(f"Random walk on a graph: degrees {deg.astype(int)}, pi = d/2|E| = {np.round(pi_g, 4)}; detailed balance holds")

# Exercise B8: the path graph 1-2-3-4 is bipartite (period 2); the lazy walk fixes it
Ap = np.array([[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0.0]])
Pp4 = Ap / Ap.sum(axis=1)[:, None]
pi_p = stationary(Pp4)
assert np.allclose(pi_p, [1 / 6, 1 / 3, 1 / 3, 1 / 6])
ev_p = np.sort(np.linalg.eigvals(Pp4).real)
assert np.allclose(ev_p, [-1, -0.5, 0.5, 1])                                  # -1 is the bipartite signature
pt = np.array([1.0, 0, 0, 0])
for t in range(1, 21):
    pt = pt @ Pp4
    assert np.abs(pt - pi_p).sum() >= 1 - 1e-12                               # never converges: support alternates
    if t == 2:
        assert np.allclose(pt, [0.5, 0, 0.5, 0]) and np.isclose(np.abs(pt - pi_p).sum(), 1.0)
assert np.isclose(np.abs(pt - pi_p).sum(), 1.0, atol=1e-9)
Plazy = 0.5 * (np.eye(4) + Pp4)
assert np.allclose(stationary(Plazy), pi_p) and np.allclose(np.sort(np.linalg.eigvals(Plazy).real), [0, 0.25, 0.75, 1])
pl = np.array([1.0, 0, 0, 0])
for _ in range(20):
    pl = pl @ Plazy
assert np.abs(pl - pi_p).sum() < 0.01
print(f"B8: path graph eigenvalues {ev_p}; after 20 steps from vertex 1, ||pi_t - pi||_1 = {np.abs(pt - pi_p).sum():.3f} (walk) "
      f"vs {np.abs(pl - pi_p).sum():.4f} (lazy walk, eigenvalues {np.round(np.sort(np.linalg.eigvals(Plazy).real), 4) + 0.0})")

# Exercise B10: Metropolis on {1,...,5} with target pi_k proportional to k, nearest-neighbour proposals
target_unnorm = np.arange(1, 6, dtype=float)                     # only ratios are ever used
pi_true = target_unnorm / target_unnorm.sum()
P_met = np.zeros((5, 5))
for i in range(5):
    for j in (i - 1, i + 1):
        if 0 <= j < 5:
            P_met[i, j] = 0.5 * min(1.0, target_unnorm[j] / target_unnorm[i])
    P_met[i, i] = 1 - P_met[i].sum()
assert np.allclose(P_met.sum(axis=1), 1) and detailed_balance_gap(P_met, pi_true) < 1e-12
assert np.allclose(pi_true @ P_met, pi_true) and np.allclose(stationary(P_met), pi_true)
assert np.allclose(P_met, [[1 / 2, 1 / 2, 0, 0, 0], [1 / 4, 1 / 4, 1 / 2, 0, 0], [0, 1 / 3, 1 / 6, 1 / 2, 0],
                           [0, 0, 3 / 8, 1 / 8, 1 / 2], [0, 0, 0, 2 / 5, 3 / 5]])


def metropolis_discrete(log_target, n_states, steps, rng, x0=0):
    """Metropolis chain on {0..n_states-1}: propose x +- 1 (symmetric), accept with prob min(1, pi_y/pi_x)."""
    x = x0; out = np.empty(steps, int)
    u_prop = rng.random(steps); u_acc = rng.random(steps)
    for t in range(steps):
        y = x + (1 if u_prop[t] < 0.5 else -1)
        if 0 <= y < n_states and np.log(u_acc[t]) < log_target[y] - log_target[x]:
            x = y
        out[t] = x
    return out


chain_m = metropolis_discrete(np.log(target_unnorm), 5, 200_000, rng)
freq_m = np.bincount(chain_m, minlength=5) / len(chain_m)
assert np.abs(freq_m - pi_true).max() < 0.01, freq_m
print(f"B10: Metropolis frequencies {np.round(freq_m, 4)} vs target {np.round(pi_true, 4)}; detailed-balance gap {detailed_balance_gap(P_met, pi_true):.1e}")
# B10(c): target k^2 on {1..10}; the code only ever sees log-ratios, the constant 385 is used only to check
k10 = np.arange(1, 11, dtype=float)
assert (k10**2).sum() == 385
chain_10 = metropolis_discrete(2 * np.log(k10), 10, 400_000, rng)
freq_10 = np.bincount(chain_10, minlength=10) / len(chain_10)
assert np.abs(freq_10 - k10**2 / 385).max() < 0.01, freq_10
print(f"B10(c): target k^2/385 on 1..10: max |freq - target| = {np.abs(freq_10 - k10**2 / 385).max():.4f}")

# PageRank on a four-page web, alpha = 0.85
links = {0: [1, 2], 1: [2], 2: [0], 3: [2]}
Hm = np.zeros((4, 4))
for i, outs in links.items():
    Hm[i, outs] = 1 / len(outs)
alpha_pr = 0.85
Gm = alpha_pr * Hm + (1 - alpha_pr) / 4
pr = np.ones(4) / 4
for it in range(200):
    pr_new = pr @ Gm
    if np.abs(pr_new - pr).sum() < 1e-12:
        break
    pr = pr_new
assert np.allclose(pr, stationary(Gm)) and np.argmax(pr) == 2
assert np.sort(np.abs(np.linalg.eigvals(Gm)))[-2] <= alpha_pr + 1e-12 and np.isclose(np.log(1e-8) / np.log(alpha_pr), 113.3, atol=0.1)
print(f"PageRank: {np.round(pr, 4)} after {it} power iterations (page 2 has the most in-links); |lambda_2| <= alpha = 0.85")

fig, ax = plt.subplots(figsize=(6, 4.2))
ax.semilogy(tt_, dists, "o-", label="weather chain ||pi_t - pi||_1")
ax.semilogy(tt_, dists[4] * ev_w[1] ** (tt_ - 5), "k--", lw=1, label=f"slope ln(lambda_2) = {np.log(ev_w[1]):.3f}")
dl = []; pl = np.array([1.0, 0, 0, 0])
for t in range(1, 31):
    pl = pl @ Plazy; dl.append(np.abs(pl - pi_p).sum())
ax.semilogy(tt_, dl, "s-", ms=4, label="lazy walk on path, ratio 0.75")
ax.set_xlabel("t"); ax.set_ylabel("l1 distance to stationary"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)
ax.set_title("Mixing is geometric with ratio |lambda_2|")
fig.tight_layout(); fig.savefig(FIG / "markov_mixing.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 9: Entropy, KL, cross-entropy, mutual information, maximum entropy
# (section 3.13; Example 4.7; exercises A5, B9, C1, C2, C3)
#
# Entropy (3.28) in bits or nats, KL (3.29) >= 0 (Gibbs), cross-entropy = H + KL
# (3.30), the Gaussian entropy (3.31) and the KL between Gaussians (3.32) checked
# by Monte Carlo and by SymPy, and the maximum-entropy comparisons at equal
# variance.


# %%
def entropy(p, base=np.e):
    """H(p) = -sum p log p, eq. (3.28), with 0 log 0 = 0."""
    p = np.asarray(p, float); p = p[p > 0]
    return float(-(p * np.log(p)).sum() / np.log(base))


def kl(p, q, base=np.e):
    """D(p || q) = sum p log(p/q), eq. (3.29); +inf if q = 0 where p > 0."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    if np.any((q == 0) & (p > 0)):
        return np.inf
    m_ = p > 0
    return float((p[m_] * np.log(p[m_] / q[m_])).sum() / np.log(base))


def cross_entropy(p, q, base=np.e):
    """H(p, q) = -sum p log q."""
    p = np.asarray(p, float); q = np.asarray(q, float); m_ = p > 0
    return float(-(p[m_] * np.log(q[m_])).sum() / np.log(base))


h_bern = lambda p: entropy([p, 1 - p], base=2)
assert np.isclose(h_bern(0.5), 1.0) and np.isclose(h_bern(0.1), 0.469, atol=5e-4) and h_bern(1e-12) < 1e-9
assert np.isclose(1 / LN2, 1.443, atol=5e-4)                                  # 1 nat = 1.443 bits
pp = np.linspace(1e-6, 1 - 1e-6, 500)
fig, ax = plt.subplots(figsize=(5.2, 3.8))
ax.plot(pp, [h_bern(t) for t in pp], lw=2); ax.plot(0.1, h_bern(0.1), "ko"); ax.text(0.12, h_bern(0.1) - 0.08, "h(0.1) = 0.469 bits")
ax.set_xlabel("p"); ax.set_ylabel("H(Bernoulli(p)) [bits]"); ax.grid(alpha=.3); ax.set_title("Binary entropy: concave, max 1 bit at 1/2, vertical at the ends")
fig.tight_layout(); fig.savefig(FIG / "bernoulli_entropy.png", dpi=120); plt.close(fig)

# Theorem 3.43(b) and Gibbs (3.45) on random distributions
for _ in range(1000):
    m_states = int(rng.integers(2, 8))
    p_r = rng.dirichlet(np.ones(m_states)); q_r = rng.dirichlet(np.ones(m_states))
    assert entropy(p_r) <= np.log(m_states) + 1e-12
    assert kl(p_r, q_r) >= 0 and kl(p_r, p_r) < 1e-12
    assert np.isclose(cross_entropy(p_r, q_r), entropy(p_r) + kl(p_r, q_r))              # (3.30)
    assert kl(p_r, np.ones(m_states) / m_states) >= 0 and np.isclose(kl(p_r, np.ones(m_states) / m_states), np.log(m_states) - entropy(p_r))
print("Theorem 3.43/3.45/(3.30) hold on 1000 random distribution pairs.")

# Example 4.7: Bernoulli KL both ways
pb, qb = np.array([0.5, 0.5]), np.array([0.9, 0.1])
D_pq, D_qp = kl(pb, qb), kl(qb, pb)
assert np.isclose(D_pq, 0.5108, atol=5e-5) and np.isclose(D_pq / LN2, 0.7370, atol=5e-5)
assert np.isclose(D_qp, 0.3681, atol=5e-5) and np.isclose(D_qp / LN2, 0.5310, atol=5e-5)
assert np.isclose(D_pq, 0.5 * np.log(25 / 9))
print(f"Example 4.7: D(Bern(.5)||Bern(.9)) = {D_pq:.4f} nats = {D_pq/LN2:.4f} bits; D(Bern(.9)||Bern(.5)) = {D_qp:.4f} nats = {D_qp/LN2:.4f} bits")


def kl_gauss_1d(m1, s1, m2, s2):
    """Eq. (3.32): D(N(m1, s1^2) || N(m2, s2^2))."""
    return np.log(s2 / s1) + (s1**2 + (m1 - m2) ** 2) / (2 * s2**2) - 0.5


def kl_mc(sample, logf1, logf2):
    """Monte Carlo KL: average of log f1(X) - log f2(X) over X ~ f1."""
    return float(np.mean(logf1(sample) - logf2(sample)))


Dg1, Dg2 = kl_gauss_1d(0, 1, 1, 2), kl_gauss_1d(1, 2, 0, 1)
assert np.isclose(Dg1, 0.4431, atol=5e-5) and np.isclose(Dg2, 1.3069, atol=5e-5) and kl_gauss_1d(0.3, 2, 0.3, 2) == 0
Xs1 = rng.standard_normal(400_000); Xs2 = 1 + 2 * rng.standard_normal(400_000)
mc1 = kl_mc(Xs1, stats.norm(0, 1).logpdf, stats.norm(1, 2).logpdf)
mc2 = kl_mc(Xs2, stats.norm(1, 2).logpdf, stats.norm(0, 1).logpdf)
assert abs(mc1 - Dg1) < 0.01 and abs(mc2 - Dg2) < 0.03, (mc1, mc2)
print(f"Example 4.7: D(N(0,1)||N(1,4)) = {Dg1:.4f} (MC {mc1:.4f}); D(N(1,4)||N(0,1)) = {Dg2:.4f} (MC {mc2:.4f})")

# SymPy: (3.31) and (3.32) by exact integration
xsym = sp.symbols("x", real=True)
mu1s, mu2s, s1s, s2s = sp.symbols("mu1 mu2 sigma1 sigma2", real=True, positive=True)
f1s = sp.exp(-(xsym - mu1s) ** 2 / (2 * s1s**2)) / sp.sqrt(2 * sp.pi * s1s**2)
f2s = sp.exp(-(xsym - mu2s) ** 2 / (2 * s2s**2)) / sp.sqrt(2 * sp.pi * s2s**2)
h_gauss = sp.simplify(sp.integrate(-f1s * sp.log(f1s), (xsym, -sp.oo, sp.oo)))
assert sp.simplify(h_gauss - sp.log(2 * sp.pi * sp.E * s1s**2) / 2) == 0
kl_sym = sp.simplify(sp.integrate(f1s * (sp.log(f1s) - sp.log(f2s)), (xsym, -sp.oo, sp.oo)))
target_sym = sp.log(s2s / s1s) + (s1s**2 + (mu1s - mu2s) ** 2) / (2 * s2s**2) - sp.Rational(1, 2)
assert sp.simplify(kl_sym - target_sym) == 0
assert np.isclose(float(h_gauss.subs(s1s, 1)), 1.4189, atol=5e-5) and np.isclose(float(h_gauss.subs(s1s, 1)) / LN2, 2.047, atol=5e-4)
print(f"SymPy: h(N(mu, sigma^2)) = {h_gauss}  -> {float(h_gauss.subs(s1s, 1)):.4f} nats = {float(h_gauss.subs(s1s, 1))/LN2:.3f} bits at sigma = 1; "
      f"(3.32) verified symbolically")


def kl_gauss_mv(m1, S1, m2, S2):
    """Multivariate KL (exercise C1): 1/2 [tr(S2^-1 S1) + (m2-m1)^T S2^-1 (m2-m1) - k + ln(det S2/det S1)]."""
    k_ = len(m1); S2i = np.linalg.inv(S2); dm = m2 - m1
    return 0.5 * (np.trace(S2i @ S1) + dm @ S2i @ dm - k_ + np.log(np.linalg.det(S2) / np.linalg.det(S1)))


mA, SA = np.zeros(2), np.array([[2.0, 1.0], [1.0, 2.0]]); mB, SB = np.array([1.0, 0.0]), np.eye(2)
D_mv = kl_gauss_mv(mA, SA, mB, SB)
assert np.isclose(D_mv, 0.9507, atol=5e-5) and np.isclose(D_mv, 0.5 * (4 + 1 - 2 + np.log(1 / 3)))
XA, _ = sample_mvn(mA, SA, 400_000, rng)
mc_mv = kl_mc(XA, stats.multivariate_normal(mA, SA).logpdf, stats.multivariate_normal(mB, SB).logpdf)
assert abs(mc_mv - D_mv) < 0.03
assert np.isclose(kl_gauss_mv(np.zeros(1), np.eye(1), np.ones(1), 4 * np.eye(1)), Dg1)     # reduces to (3.32) for k = 1
print(f"C1: D(N(0, [[2,1],[1,2]]) || N((1,0), I)) = {D_mv:.4f} (MC {mc_mv:.4f})")

# C2: mutual information of a jointly Gaussian pair, I = -1/2 ln(1 - rho^2)
rho_mi = 0.5
Sig_mi = np.array([[1, rho_mi], [rho_mi, 1.0]])
I_th = -0.5 * np.log(1 - rho_mi**2)
assert np.isclose(I_th, 0.1438, atol=5e-5) and np.isclose(I_th / LN2, 0.2075, atol=5e-5)
XY, _ = sample_mvn(np.zeros(2), Sig_mi, 400_000, rng)
I_mc = np.mean(stats.multivariate_normal(np.zeros(2), Sig_mi).logpdf(XY) - stats.norm.logpdf(XY[:, 0]) - stats.norm.logpdf(XY[:, 1]))
assert abs(I_mc - I_th) < 0.01
h_joint = 0.5 * np.log((2 * np.pi * np.e) ** 2 * np.linalg.det(Sig_mi))
assert np.isclose(2 * 0.5 * np.log(2 * np.pi * np.e) - h_joint, I_th)
print(f"C2: I(X;Y) for rho = 0.5: {I_th:.4f} nats = {I_th/LN2:.4f} bits (MC {I_mc:.4f})")

# Maximum entropy at equal variance sigma^2 = 1 (Theorem 3.46), by numerical integration of -f ln f
h_num = lambda pdf, lo, hi: quad(lambda t: -pdf(t) * np.log(pdf(t)) if pdf(t) > 0 else 0.0, lo, hi, limit=200)[0]
h_gaussian = h_num(stats.norm(0, 1).pdf, -12, 12)
h_uniform = h_num(stats.uniform(-np.sqrt(3), 2 * np.sqrt(3)).pdf, -np.sqrt(3), np.sqrt(3))
h_laplace = h_num(stats.laplace(0, 1 / np.sqrt(2)).pdf, -25, 25)
h_expo = h_num(stats.expon(0, 1).pdf, 0, 40)
assert np.isclose(h_gaussian, 1.4189, atol=1e-4) and np.isclose(h_uniform, np.log(2 * np.sqrt(3)), atol=1e-6)
assert np.isclose(h_laplace, 1 + np.log(np.sqrt(2)), atol=1e-4) and np.isclose(h_expo, 1.0, atol=1e-4)
assert h_gaussian > h_laplace > h_uniform > h_expo
assert np.isclose(h_num(stats.uniform(0, 0.5).pdf, 0, 0.5), -np.log(2))                     # differential entropy can be negative
# the exponential maximizes on [0, inf) at fixed mean 1: compare with a Gamma(2, 2) (same mean), h = 1 + ln(mu) = 1
h_gamma22 = h_num(stats.gamma(2, scale=0.5).pdf, 0, 40)
assert h_gamma22 < h_expo
print(f"Max entropy at var 1: Gaussian {h_gaussian:.4f} > Laplace {h_laplace:.4f} > Uniform {h_uniform:.4f} > Exponential {h_expo:.4f}; "
      f"Uniform(0,1/2) has h = {h_num(stats.uniform(0, 0.5).pdf, 0, 0.5):.4f} < 0; on [0,inf) with mean 1: Exp {h_expo:.4f} > Gamma(2,2) {h_gamma22:.4f}")

# Exercise A5
assert np.isclose(entropy(np.ones(6) / 6, 2), 2.585, atol=5e-4) and np.isclose(entropy([0.5, 0.25, 0.25], 2), 1.5)
assert np.isclose(kl([0.5, 0.25, 0.25], np.ones(3) / 3, 2), 0.0850, atol=5e-5) and np.isclose(kl(np.ones(3) / 3, [0.5, 0.25, 0.25], 2), 0.0817, atol=5e-5)
assert np.isclose(cross_entropy([0.5, 0.25, 0.25], np.ones(3) / 3, 2), np.log2(3))
print(f"A5: H(die) = {entropy(np.ones(6)/6, 2):.4f} bits; H(1/2,1/4,1/4) = 1.5; D(p||u) = {kl([.5,.25,.25], np.ones(3)/3, 2):.4f}, "
      f"D(u||p) = {kl(np.ones(3)/3, [.5,.25,.25], 2):.4f} bits; H(p, u) = log2 3 = {np.log2(3):.4f}")

# Exercise B9: conditioning reduces entropy on average, but not slice by slice
PB9 = np.array([[0.9, 0.0], [0.05, 0.05]])          # rows X in {0,1}, cols Y in {0,1}
HY = entropy(marginal(PB9, 1), 2); HX = entropy(marginal(PB9, 0), 2); HXY = entropy(PB9.ravel(), 2)
HY_given_X1 = entropy(conditional(PB9, "X", 1), 2); HY_given_X0 = entropy(conditional(PB9, "X", 0), 2)
HY_given_X = marginal(PB9, 0)[0] * HY_given_X0 + marginal(PB9, 0)[1] * HY_given_X1
I_B9 = HY - HY_given_X
assert np.isclose(HY, 0.2864, atol=5e-5) and np.isclose(HY_given_X1, 1.0) and np.isclose(HY_given_X, 0.1)
assert HY_given_X1 > HY and HY_given_X <= HY and np.isclose(I_B9, 0.1864, atol=5e-5)
assert np.isclose(HXY, HX + HY_given_X) and np.isclose(HXY, 0.5690, atol=5e-5)
assert np.isclose(I_B9, kl(PB9.ravel(), np.outer(marginal(PB9, 0), marginal(PB9, 1)).ravel(), 2))
print(f"B9: H(Y) = {HY:.4f}, H(Y|X=1) = {HY_given_X1:.1f} > H(Y), but H(Y|X) = {HY_given_X:.2f} <= H(Y); I = {I_B9:.4f}; H(X,Y) = {HXY:.4f} bits")

# Exercise C3: maximum entropy on {1..6} with mean 4.5 -> p_i proportional to e^{-beta i} (Jaynes' die)
faces = np.arange(1, 7, dtype=float)
mean_given_beta = lambda b: (faces * np.exp(-b * faces)).sum() / np.exp(-b * faces).sum()
beta_die = brentq(lambda b: mean_given_beta(b) - 4.5, -5, 5)
p_die = np.exp(-beta_die * faces); p_die /= p_die.sum()
assert np.isclose(beta_die, -0.3710, atol=5e-5) and np.isclose(faces @ p_die, 4.5)
assert np.allclose(p_die, [0.0544, 0.0788, 0.1142, 0.1654, 0.2398, 0.3475], atol=5e-5)
assert np.isclose(entropy(p_die, 2), 2.328, atol=5e-4) and entropy(p_die, 2) < np.log2(6)
for _ in range(500):                                              # no other distribution with mean 4.5 beats it
    q_ = rng.dirichlet(np.ones(6))
    q_ = 0.5 * q_ + 0.5 * p_die                                   # mix toward p to keep the mean near 4.5...
    q_ = q_ * np.exp(-brentq(lambda b: (faces * q_ * np.exp(-b * faces)).sum() / (q_ * np.exp(-b * faces)).sum() - 4.5, -20, 20) * faces)
    q_ /= q_.sum()                                                # ... then tilt it exactly onto mean 4.5
    assert np.isclose(faces @ q_, 4.5) and entropy(q_, 2) <= entropy(p_die, 2) + 1e-9
print(f"C3: loaded die with mean 4.5: beta = {beta_die:.4f}, p = {np.round(p_die, 4)}, H = {entropy(p_die, 2):.3f} bits < log2 6 = {np.log2(6):.3f}")

# Intermediate values quoted in solutions.md (A1, A5, B3, B4, B9, C3), so every printed digit is checked.
assert np.isclose(marginal(PA1, 0)[0] * marginal(PA1, 1)[0], 0.12)                                   # A1
assert np.isclose(np.log2(2 / 3), -0.5850, atol=5e-5) and np.isclose(2 * np.log2(4 / 3), 0.8301, atol=5e-5)  # A5
assert np.isclose(stats.norm.sf(3), 0.00135, atol=5e-6) and np.isclose(stats.norm.sf(-0.5), 0.69146, atol=5e-6)  # B3(c)
assert np.isclose(stats.norm.pdf(3, 4, 2), 0.17603, atol=5e-6) and np.isclose(stats.norm.pdf(3), 0.004432, atol=5e-7)  # B3(d)
assert np.isclose(0.3 * stats.norm.pdf(3, 4, 2), 0.05281, atol=5e-6) and np.isclose(0.7 * stats.norm.pdf(3), 0.00310, atol=5e-6)
assert np.isclose(-5.328 * 0.7 + 55.552 * 0.3, 12.936) and np.isclose(3 * 2.8 * 4 + 2.8**3, 55.552)   # B3(e)
assert np.isclose(post.var(), 0.01709, atol=5e-6) and np.isclose(post.var(), 32 / (144 * 13))         # B4(a)
assert np.isclose(-0.95 * np.log2(0.95), 0.0703, atol=5e-5) and np.isclose(-0.05 * np.log2(0.05), 0.2161, atol=5e-5)  # B9
assert np.isclose(-np.log2(0.9), 0.1520, atol=5e-5) and np.isclose(-np.log2(0.05), 4.3219, atol=5e-5)
assert np.isclose(np.exp(-beta_die), 1.449, atol=5e-4) and np.isclose(np.log2(6) - entropy(p_die, 2), 0.257, atol=5e-4)  # C3
print("Solutions: intermediate values for A1, A5, B3, B4, B9, C3 check out.")

# Intermediate values quoted in the notes (sections 3.11 and 4), so those digits are checked too.
assert np.isclose(np.sqrt(0.6 * 0.25), 0.3873, atol=5e-5)                                          # Ex 4.1
assert np.isclose(-2 * np.log(0.5), 1.3863, atol=5e-5) and np.isclose(-2 * np.log(0.1), 4.6052, atol=5e-5)  # Ex 4.2
assert np.isclose(np.sqrt(-2 * np.log(0.1)), 2.1460, atol=5e-5) and np.isclose(np.sqrt(3), 1.732, atol=5e-4)  # Ex 4.2, 4.4
assert np.isclose(2.5 / np.sqrt(4.6875), 1.155, atol=5e-4)                                          # Ex 4.5: z = 1.155
assert np.isclose(np.log(5), 1.6094, atol=5e-5) and np.isclose(np.log(1.8), 0.5878, atol=5e-5)      # Ex 4.7
assert np.isclose(stats.norm.sf(5), 3e-7, rtol=0.05)                                                # 3.11: 1 - Phi(5)
print(f"Notes: sqrt(0.15) = {np.sqrt(0.15):.4f}; ln 5 = {np.log(5):.4f}, ln 1.8 = {np.log(1.8):.4f}; "
      f"z = {2.5 / np.sqrt(4.6875):.3f}; 1 - Phi(5) = {stats.norm.sf(5):.2e}")

# %% [markdown]
# ## Your turn
#
# 1. **2D random walk return probability.** Simulate a simple random walk on Z^2
#    (steps to the four neighbours with probability 1/4) for 10^4 walks of 10^4
#    steps each, and estimate the probability of having returned to the origin by
#    step T for T = 10, 100, 1000, 10000. The walk is recurrent (the probability
#    tends to 1 as T grows) but the convergence is logarithmically slow. Then do
#    Z^3: the return probability converges to about 0.34 (Pólya's theorem). Use the
#    bipartite warning of B8: returns can only happen at even steps.
# 2. **Gaussian-mixture entropy by Monte Carlo.** The mixture of exercise B3 has no
#    closed-form entropy. Estimate h(Y) = -E[log f_Y(Y)] by sampling Y and
#    averaging -log f_Y, with f_Y the two-component density. Compare with the
#    entropy of a single Gaussian with the same variance 5.26 (that is
#    (1/2) ln(2 pi e 5.26) = 2.249 nats); Theorem 3.46 says the mixture must be
#    lower. Then shrink the separation mu_1 - mu_0 toward 0 and watch the gap close.
# 3. **Mutual information of a Gaussian pair.** Repeat the C2 Monte Carlo for
#    rho in {0, 0.3, 0.6, 0.9, 0.99} and plot against -1/2 ln(1 - rho^2). Then
#    replace Y by X^2 with X ~ N(0,1): the correlation is zero, but estimate I(X; X^2)
#    by binning both variables into 30 bins and computing the discrete mutual
#    information of the histogram. It is far from zero. (The exact I(X; X^2) is in
#    fact infinite, since X^2 is a deterministic function of a continuous X, so
#    the binned estimate keeps growing as you refine the bins.) Covariance misses
#    this; mutual information does not.
# 4. **Metropolis on a continuous target.** Change `metropolis_discrete` to a
#    Gaussian random-walk proposal on R and sample from a target known only up to a
#    constant, e.g. exp(-x^4/4 - x^2/2). Check that detailed balance still holds
#    (it is the same proof) and compare the acceptance rate for proposal scales
#    0.1, 1, 10. ST-6 does this properly.

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T_START:.1f} s. Figures in {FIG}")
