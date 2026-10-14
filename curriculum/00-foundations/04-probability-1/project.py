# %% [markdown]
# # F4 Project: Probability by Simulation and by Formula
#
# Every claim in the F4 notes that can be checked numerically is checked here,
# twice where possible: once from the exact formula derived in the notes and
# once by Monte Carlo with a seeded generator. Parts: (1) the traps, (2) exact
# PMFs vs scipy.stats and Examples 4.9-4.10, (3) the Poisson limit, (4)
# inverse-CDF sampling, (5) the Gaussian integral, (6) 68-95-99.7, (7) MGFs,
# (8) the law of large numbers and the Cauchy, (9) linearity on dependent
# variables, (10) every number in the exercise solutions. Run headless:
# `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import itertools
import math
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from scipy import integrate, special, stats

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
E = math.e


def check(est, target, se, k=3.0, label=""):
    """Assert a Monte Carlo estimate lies within k standard errors of the exact value."""
    z = abs(est - target) / se
    assert z < k, f"{label}: {est:.5f} vs {target:.5f} is {z:.2f} SE away"
    print(f"  {label}: simulated {est:.5f}, exact {target:.5f} ({z:.2f} SE)")


def prop_se(p, n):
    """Standard deviation of a sample proportion, sqrt(pq/n) (Bin(n, p)/n)."""
    return math.sqrt(p * (1 - p) / n)


# %% [markdown]
# ## Part 1: The traps, exactly and by Monte Carlo (§3.3-3.7, Examples 4.1-4.5)

# %%
def birthday_exact(k, days=365):
    """P(at least one shared birthday) = 1 - prod_{j<k}(1 - j/days), §3.3 step 3."""
    return 1 - np.prod(1 - np.arange(k) / days)


def birthday_approx(k, days=365):
    """Approximation (3.3): 1 - exp(-k(k-1)/(2 days))."""
    return 1 - math.exp(-k * (k - 1) / (2 * days))


print("Part 1: traps")
assert round(birthday_exact(23), 4) == 0.5073 and round(birthday_approx(23), 4) == 0.5000
assert round(float(np.log(1 - np.arange(23) / 365).sum()), 4) == -0.7078
assert round(birthday_exact(50), 4) == 0.9704 and round(birthday_approx(50), 4) == 0.9651
assert all(birthday_approx(k) <= birthday_exact(k) + 1e-15 for k in range(1, 100))  # 1 - x <= e^{-x}
assert round(math.comb(23, 2) / 365, 3) == 0.693  # expected matching pairs, indicator trick (2)
print(f"  birthday k=23: exact {birthday_exact(23):.5f}, approx {birthday_approx(23):.5f}")

ks, T = np.arange(1, 61), 20_000
sim = np.array([(np.diff(np.sort(rng.integers(0, 365, (T, k)), axis=1), axis=1) == 0).any(axis=1).mean()
                for k in ks])
exact = np.array([birthday_exact(k) for k in ks])
se = np.sqrt(exact * (1 - exact) / T) + 1e-12
assert np.max(np.abs(sim - exact) / np.maximum(se, 1 / T)) < 4.5
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(ks, exact, "k-", label="exact product")
ax.plot(ks, [birthday_approx(k) for k in ks], "C1--", label=r"$1-e^{-k(k-1)/730}$")
ax.plot(ks, sim, "C0.", ms=4, label=f"simulated ({T} rooms each)")
ax.axhline(0.5, color="gray", lw=0.5); ax.axvline(23, color="gray", lw=0.5)
ax.set(xlabel="people k", ylabel="P(shared birthday)", title="Birthday problem"); ax.legend()
fig.savefig(FIG / "birthday.png", dpi=120); plt.close(fig)


def posterior(prior, sens=0.99, spec=0.99):
    """Bayes (3.10) with the LOTP denominator (3.11): P(D | +)."""
    return sens * prior / (sens * prior + (1 - spec) * (1 - prior))


p1 = posterior(1e-4); p2 = posterior(p1); p3 = posterior(p2)
assert round(p1, 6) == 0.009804 and round(p2, 4) == 0.4950 and round(p3, 3) == 0.990
assert round(0.99e-4 + 0.01 * 0.9999, 6) == 0.010098
assert [round(9999 / 99), round(99**2 / 9999, 2), round(99**3 / 9999)] == [101, 0.98, 97]  # odds after 1, 2, 3 positives
M = 2_000_000
sick = rng.random(M) < 1e-4
pos = np.where(sick, rng.random(M) < 0.99, rng.random(M) < 0.01)
check(sick[pos].mean(), p1, prop_se(p1, pos.sum()), label="P(D|+) among 2e6 simulated patients")

# Monty Hall: you pick door 0. A knowing host opens a goat door among {1, 2}.
N = 100_000
car = rng.integers(0, 3, N)
host = np.where(car == 0, rng.integers(1, 3, N), 3 - car)
sig = prop_se(2 / 3, N)
assert round(sig, 5) == 0.00149
check(((3 - host) == car).mean(), 2 / 3, sig, label="Monty Hall, switch wins")
check((car == 0).mean(), 1 / 3, sig, label="Monty Hall, stay wins")
opened = host == 2  # "host opens door 3" in the notes' 1-based labels
check((car[opened] == 1).mean(), 2 / 3, prop_se(2 / 3, opened.sum()), label="P(car behind 2 | opens 3)")
host_r = rng.integers(1, 3, N); goat = host_r != car  # random host who happened to show a goat
check(((3 - host_r)[goat] == car[goat]).mean(), 0.5, prop_se(0.5, goat.sum()), label="random host, switch wins")


def p_some_match(n):
    """(3.8): P(at least one fixed point of a uniform permutation of n) = 1 - sum_{k<=n} (-1)^k/k!."""
    return 1 - sum((-1) ** k / math.factorial(k) for k in range(n + 1))


assert math.isclose(p_some_match(4), 0.625)
assert sum(all(s[i] != i for i in range(4)) for s in itertools.permutations(range(4))) == 9  # derangements
assert round(p_some_match(7), 4) == round(1 - 1 / E, 4) == 0.6321
assert abs(p_some_match(10) - (1 - 1 / E)) < 1 / math.factorial(11) < 2.6e-8
print("  P(some match), n = 1..10:", [round(p_some_match(n), 6) for n in range(1, 11)])
perms = rng.permuted(np.tile(np.arange(10), (N, 1)), axis=1)
fixed = (perms == np.arange(10)).sum(axis=1)
check(fixed.mean(), 1.0, 1 / math.sqrt(N), label="mean fixed points, n=10 (Var = 1 by C4)")
check((fixed > 0).mean(), p_some_match(10), prop_se(0.632, N), label="P(at least one match), n=10")

# Two coins (§3.6), two children (Trap 3), Simpson (Trap 4), prosecutor (Trap 5), St Petersburg.
h, b = Fr(1, 2), Fr(3, 4)
PH1, PHH = h * h + h * b, h * h ** 2 + h * b ** 2
assert PH1 == Fr(5, 8) and PHH == Fr(13, 32) and PH1 ** 2 == Fr(25, 64) and PHH != PH1 ** 2
kids = list(itertools.product("BG", repeat=2))           # (elder, younger), each 1/4
met = [(c, i) for c in kids for i in (0, 1)]             # also which child you meet, each 1/8
cond = lambda A, B, space: Fr(sum(A(w) and B(w) for w in space), sum(B(w) for w in space))
assert cond(lambda c: c == ("B", "B"), lambda c: "B" in c, kids) == Fr(1, 3)
assert cond(lambda c: c == ("B", "B"), lambda c: c[0] == "B", kids) == Fr(1, 2)
assert cond(lambda w: w[0] == ("B", "B"), lambda w: w[0][w[1]] == "B", met) == Fr(1, 2)
simpson = {"A": [(81, 87), (192, 263)], "B": [(234, 270), (55, 80)]}
for arm, strata in simpson.items():
    s, n_ = sum(x for x, _ in strata), sum(y for _, y in strata)
    print(f"  Simpson arm {arm}: small {strata[0][0]/strata[0][1]:.0%}, large {strata[1][0]/strata[1][1]:.0%},"
          f" all {s/n_:.0%}, P(large | arm) = {strata[1][1]/n_:.2f}")
assert [round(263 / 350, 2), round(80 / 350, 2)] == [0.75, 0.23]
assert abs(1e-7 / (1e-7 + 1e-6 * (1 - 1e-7)) - 1 / 11) < 1e-6
assert 40 + sum(2 ** 40 / 2 ** n for n in range(41, 200)) == 41.0  # capped St Petersburg

# %% [markdown]
# ## Part 2: Exact PMFs from (3.20)-(3.24) vs scipy.stats; Examples 4.9 and 4.10
#
# **Build it yourself first:** write the five PMFs below with `math.comb`.
# Convention traps: our Geom counts failures (SciPy's `geom` counts trials, so
# shift with loc=-1; NumPy's `rng.geometric` also counts trials). SciPy's
# `hypergeom(M, n, N)` means (population, white, draws).

# %%
def binom_pmf(k, n, p):
    """(3.20) C(n,k) p^k q^(n-k)."""
    return math.comb(n, k) * p ** k * (1 - p) ** (n - k)


def geom_pmf(k, p):
    """(3.21) q^k p: k failures before the first success."""
    return (1 - p) ** k * p


def nbin_pmf(k, r, p):
    """(3.22) C(k+r-1, r-1) p^r q^k: k failures before the r-th success."""
    return math.comb(k + r - 1, r - 1) * p ** r * (1 - p) ** k


def hgeom_pmf(k, w, b, n):
    """(3.23) C(w,k) C(b,n-k) / C(w+b,n)."""
    return math.comb(w, k) * math.comb(b, n - k) / math.comb(w + b, n) if 0 <= n - k else 0.0


def pois_pmf(k, lam):
    """(3.24) e^-lam lam^k / k!."""
    return math.exp(-lam) * lam ** k / math.factorial(k)


print("Part 2: named distributions")
p, q, n, r, lam, (w, bk, m) = 0.3, 0.7, 12, 3, 4.0, (7, 13, 8); Nw, pw = w + bk, w / (w + bk)
cases = {  # name: (our pmf, scipy dist, support, mean from §3.10 table, var from table, sampler)
    "Bin(12,0.3)": (lambda k: binom_pmf(k, n, p), stats.binom(n, p), range(n + 1), n * p, n * p * q,
                    lambda s: rng.binomial(n, p, s)),
    "Geom(0.3)": (lambda k: geom_pmf(k, p), stats.geom(p, loc=-1), range(200), q / p, q / p ** 2,
                  lambda s: rng.geometric(p, s) - 1),
    "NBin(3,0.3)": (lambda k: nbin_pmf(k, r, p), stats.nbinom(r, p), range(300), r * q / p, r * q / p ** 2,
                    lambda s: rng.negative_binomial(r, p, s)),
    "HGeom(7,13,8)": (lambda k: hgeom_pmf(k, w, bk, m), stats.hypergeom(Nw, w, m), range(m + 1), m * pw,
                      m * pw * (1 - pw) * (Nw - m) / (Nw - 1), lambda s: rng.hypergeometric(w, bk, m, s)),
    "Pois(4)": (lambda k: pois_pmf(k, lam), stats.poisson(lam), range(120), lam, lam,
                lambda s: rng.poisson(lam, s)),
    "DUnif{1..6}": (lambda k: 1 / 6, stats.randint(1, 7), range(1, 7), 3.5, 35 / 12,
                    lambda s: rng.integers(1, 7, s)),
}
fig, axes = plt.subplots(2, 3, figsize=(12, 6.5))
for ax, (name, (pmf, dist, supp, mu, var, sampler)) in zip(axes.flat, cases.items()):
    k = np.array(supp); pk = np.array([pmf(int(i)) for i in k])
    assert np.allclose(pk, dist.pmf(k), atol=1e-12), name
    assert abs(pk.sum() - 1) < 1e-10
    mu_pmf = (k * pk).sum(); var_pmf = (k ** 2 * pk).sum() - mu_pmf ** 2   # (3.16), LOTUS, (3.17)
    assert math.isclose(mu_pmf, mu, rel_tol=1e-9) and math.isclose(var_pmf, var, rel_tol=1e-8), name
    assert math.isclose(dist.mean(), mu) and math.isclose(dist.var(), var)
    x = sampler(100_000)
    check(x.mean(), mu, math.sqrt(var / x.size), k=4, label=f"{name} sample mean")
    cut = int(np.searchsorted(np.cumsum(pk), 0.9995)) + 2
    kk = k[:cut]
    ax.bar(kk, np.bincount(x, minlength=kk.max() + 1)[kk] / x.size, color="C0", alpha=0.5, label="sample")
    ax.plot(kk, pk[:cut], "ko", ms=3, label="formula")
    ax.set_title(f"{name}: mean {mu:.3g}, var {var:.3g}"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(FIG / "named_pmfs.png", dpi=120); plt.close(fig)
assert round(math.sqrt(35 / 12), 3) == 1.708

# Example 4.9: 5 defective in a batch of 50, inspect 10.
P0 = Fr(math.comb(45, 10), math.comb(50, 10))
assert P0 == Fr(40 * 39 * 38 * 37 * 36, 50 * 49 * 48 * 47 * 46) == Fr(78_960_960, 254_251_200)
var_h = 10 * 0.1 * 0.9 * 40 / 49
assert round(float(P0), 5) == 0.31056 and round(1 - float(P0), 4) == 0.6894 and round(var_h, 4) == 0.7347
assert round(0.9 ** 10, 4) == 0.3487 and round(0.9 ** 10 - float(P0), 3) == 0.038
assert math.isclose(0.9 / var_h, 1.225) and round(40 / 49, 3) == 0.816
accept = [math.comb(45, s) / math.comb(50, s) for s in range(51)]
n_star = next(s for s in range(51) if accept[s] < 0.05)
assert n_star == 22 and round(accept[22], 4) == 0.0464 and round(accept[21], 4) == 0.0560
print(f"  Example 4.9: P(X=0) = {float(P0):.5f} (Bin approx {0.9**10:.4f}); need n = {n_star} for P(accept) < 0.05")
fig, ax = plt.subplots(figsize=(6, 3.5)); kk = np.arange(6)
ax.bar(kk - 0.2, [hgeom_pmf(k, 5, 45, 10) for k in kk], 0.4, label="HGeom(5, 45, 10)")
ax.bar(kk + 0.2, [binom_pmf(k, 10, 0.1) for k in kk], 0.4, label="Bin(10, 0.1)")
ax.set(xlabel="defective parts in sample", title="Example 4.9"); ax.legend()
fig.savefig(FIG / "defects_hgeom_vs_binom.png", dpi=120); plt.close(fig)

# Example 4.10: bearing lifetime T ~ Expo(1/2000 per hour).
lam_b = 1 / 2000
assert round(math.exp(-3000 * lam_b), 4) == 0.2231 and round(1 - math.exp(-500 * lam_b), 4) == 0.2212
assert round(math.log(2) / lam_b) == 1386
life = rng.exponential(1 / lam_b, 100_000); alive = life > 3000
check((life[alive] > 6000).mean(), math.exp(-1.5), prop_se(0.2231, alive.sum()), label="P(T>6000 | T>3000)")

# %% [markdown]
# ## Part 3: The Poisson limit (§3.10, Example 4.6)

# %%
print("Part 3: Poisson limit")
table = {0: (0.049563, 0.049787), 1: (0.149137, 0.149361), 2: (0.224154, 0.224042),
         3: (0.224379, 0.224042), 4: (0.168284, 0.168031), 5: (0.100869, 0.100819)}
for k, (bv, pv) in table.items():
    bb, pp = binom_pmf(k, 1000, 0.003), pois_pmf(k, 3)
    assert round(bb, 6) == bv and round(pp, 6) == pv
    print(f"  k={k}: Bin {bb:.6f}  Pois {pp:.6f}  diff {bb - pp:+.6f}")
assert round(1000 * math.log(0.997), 4) == -3.0045


def tv_bin_pois(n, lam=3.0):
    """Total variation distance (1/2) sum_k |Bin(n, lam/n)(k) - Pois(lam)(k)|."""
    k = np.arange(n + 1)
    return 0.5 * (np.abs(stats.binom.pmf(k, n, lam / n) - stats.poisson.pmf(k, lam)).sum()
                  + stats.poisson.sf(n, lam))


ns = np.unique(np.logspace(1, 4, 25).astype(int)); tv = np.array([tv_bin_pois(v) for v in ns])
assert np.all(np.diff(tv) < 0) and np.all((ns * tv)[ns >= 100] > 0.7) and np.all((ns * tv)[ns >= 100] < 0.8)
fig, ax = plt.subplots(figsize=(5.5, 4))
ax.loglog(ns, tv, "o-", label="TV(Bin(n, 3/n), Pois(3))"); ax.loglog(ns, 0.75 / ns, "k--", label="0.75 / n")
ax.set(xlabel="n", title="Poisson approximation error"); ax.legend()
fig.savefig(FIG / "poisson_tv.png", dpi=120); plt.close(fig)
print(f"  n * TV at n = 10^4: {1e4 * tv[-1]:.4f} (error decays like 1/n)")

# %% [markdown]
# ## Part 4: Universality of the Uniform (Theorem 3.15, Example 4.7)
#
# **Build it yourself first:** `expo_sample(lam, size, rng)` from X = -ln(1-U)/lam.

# %%
def expo_sample(lam, size, rng):
    """Inverse-CDF sampling, Theorem 3.15(a): F^{-1}(u) = -ln(1-u)/lam."""
    return -np.log1p(-rng.random(size)) / lam


def piece_cdf(x):
    """F(x) = x/2 on [0,1], 1/2 + (x-1)/4 on [1,3] (Example 4.7)."""
    x = np.clip(x, 0, 3)
    return np.where(x < 1, x / 2, 0.5 + (x - 1) / 4)


def piece_inv(u):
    """F^{-1}(u) = 2u for u < 1/2, 1 + 4(u - 1/2) for u >= 1/2."""
    return np.where(u < 0.5, 2 * u, 1 + 4 * (u - 0.5))


print("Part 4: inverse-CDF sampling")
xe = expo_sample(2.0, 100_000, rng); xp = piece_inv(rng.random(100_000))
check(xe.mean(), 0.5, 0.5 / math.sqrt(xe.size), label="Expo(2) via -ln(1-U)/2, mean")
pe = stats.kstest(xe, stats.expon(scale=0.5).cdf).pvalue; pp_ = stats.kstest(xp, piece_cdf).pvalue
z = rng.standard_normal(100_000); pu = stats.kstest(stats.norm.cdf(z), "uniform").pvalue  # Theorem 3.15(b)
assert min(pe, pp_, pu) > 1e-3, (pe, pp_, pu)
print(f"  KS p-values: Expo {pe:.3f}, piecewise {pp_:.3f}, Phi(Z) vs Unif {pu:.3f}")
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for ax, x, F, name in [(axes[0], xe, stats.expon(scale=0.5).cdf, "Expo(2)"), (axes[1], xp, piece_cdf, "piecewise")]:
    xs = np.sort(x)[::500]; grid = np.linspace(0, xs.max(), 300)
    ax.step(xs, (np.arange(xs.size) + 1) / xs.size, where="post", label="empirical CDF")
    ax.plot(grid, F(grid), "k--", label="target F"); ax.set_title(f"Inverse-CDF sampling: {name}"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "inverse_cdf.png", dpi=120); plt.close(fig)

# %% [markdown]
# ## Part 5: The Gaussian integral (Theorem 3.18, Example 4.8) and Gamma (Theorem 3.19)

# %%
print("Part 5: Gaussian integral")
I1 = integrate.quad(lambda x: np.exp(-x * x / 2), -np.inf, np.inf)[0]
I2 = integrate.dblquad(lambda y, x: np.exp(-(x * x + y * y) / 2), -np.inf, np.inf, -np.inf, np.inf)[0]
Ipolar = integrate.dblquad(lambda rr, th: np.exp(-rr * rr / 2) * rr, 0, 2 * np.pi, 0, np.inf)[0]
assert abs(I1 - math.sqrt(2 * math.pi)) < 1e-10 and round(I1, 6) == 2.506628
assert abs(I2 - 2 * math.pi) < 1e-6 and abs(Ipolar - 2 * math.pi) < 1e-6 and round(I2, 6) == 6.283185
phi = stats.norm.pdf
moments = [integrate.quad(lambda x, j=j: x ** j * phi(x), -np.inf, np.inf)[0] for j in (0, 2, 4)]
assert np.allclose(moments, [1, 1, 3])
assert math.isclose(special.gamma(0.5), math.sqrt(math.pi)) and round(math.sqrt(math.pi), 6) == 1.772454
assert math.isclose(integrate.quad(lambda x: x ** -0.5 * np.exp(-x), 0, np.inf)[0], math.sqrt(math.pi), rel_tol=1e-8)
assert all(math.isclose(special.gamma(k), math.factorial(k - 1)) for k in range(1, 12))
assert math.isclose(special.gamma(2.5), 0.75 * math.sqrt(math.pi))
print(f"  int e^(-z^2/2) = {I1:.10f}, sqrt(2 pi) = {math.sqrt(2*math.pi):.10f}; double integral {I2:.6f}")

# %% [markdown]
# ## Part 6: 68-95-99.7 (§3.12)

# %%
print("Part 6: 68-95-99.7")
zs = rng.standard_normal(1_000_000)
for kk_, target in [(1, 0.6827), (2, 0.9545), (3, 0.9973)]:
    exact_k = 2 * stats.norm.cdf(kk_) - 1
    assert round(exact_k, 4) == target
    check((np.abs(zs) <= kk_).mean(), exact_k, prop_se(exact_k, zs.size), label=f"P(|Z| <= {kk_})")
assert round(stats.norm.ppf(0.99), 3) == 2.326

# %% [markdown]
# ## Part 7: MGFs (§3.14)
#
# The sample MGF M^(t) = mean(exp(t X_i)) is a smooth function of t whose
# derivatives at 0 are the sample moments (Theorem 3.20 applied to the
# empirical distribution). Finite differences recover them; the sample moments
# in turn estimate the exact moments 2, 6, 22 of Pois(2).

# %%
print("Part 7: MGFs")
t, lam_s, nn, mm, pp_s = sp.symbols("t lambda n m p", positive=True)
Mpois = sp.exp(lam_s * (sp.exp(t) - 1))
exact_m = [sp.simplify(sp.diff(Mpois, t, j).subs(t, 0)) for j in (1, 2, 3)]
assert [sp.expand(e - f) for e, f in zip(exact_m, [lam_s, lam_s**2 + lam_s, lam_s**3 + 3*lam_s**2 + lam_s])] == [0, 0, 0]
assert [int(e.subs(lam_s, 2)) for e in exact_m] == [2, 6, 22]
Mb = lambda k: (1 - pp_s + pp_s * sp.exp(t)) ** k
assert sp.simplify(Mb(nn) * Mb(mm) - Mb(nn + mm)) == 0                         # Binomials add
assert sp.simplify(sp.limit((1 + lam_s * (sp.exp(t) - 1) / nn) ** nn, nn, sp.oo) - Mpois) == 0  # Poisson limit
zz, tr = sp.symbols("z t_r", real=True)
assert sp.simplify(sp.integrate(sp.exp(tr * zz - zz**2 / 2), (zz, -sp.oo, sp.oo)) / sp.sqrt(2 * sp.pi)
                   - sp.exp(tr**2 / 2)) == 0                                      # Normal MGF
assert [sp.series(sp.exp(t**2 / 2), t, 0, 7).removeO().coeff(t, j) * sp.factorial(j) for j in (2, 4, 6)] == [1, 3, 15]
assert [sp.series(1 / (1 - t / lam_s), t, 0, 4).removeO().coeff(t, j) * sp.factorial(j) for j in (1, 2, 3)] \
    == [1 / lam_s, 2 / lam_s**2, 6 / lam_s**3]                                   # Expo: n!/lam^n
X = rng.poisson(2.0, 1_000_000).astype(float)
Mhat = lambda s: np.exp(s * X).mean()
hh = 1e-2
d1 = (Mhat(hh) - Mhat(-hh)) / (2 * hh)
d2 = (Mhat(hh) - 2 * Mhat(0) + Mhat(-hh)) / hh ** 2
d3 = (Mhat(2 * hh) - 2 * Mhat(hh) + 2 * Mhat(-hh) - Mhat(-2 * hh)) / (2 * hh ** 3)
for j, d, ex in [(1, d1, 2), (2, d2, 6), (3, d3, 22)]:
    smp = (X ** j).mean()
    assert abs(d - smp) / smp < 2e-3                       # M^(j)(0) = j-th sample moment
    check(smp, ex, (X ** j).std() / math.sqrt(X.size), k=4, label=f"E[X^{j}] for Pois(2), FD gives {d:.4f}")

# %% [markdown]
# ## Part 8: Law of large numbers works for Bernoulli, fails for Cauchy

# %%
print("Part 8: running means")
nsteps = np.arange(1, 100_001)
bern_run = np.cumsum(rng.random(nsteps.size) < 0.3) / nsteps
check(bern_run[-1], 0.3, prop_se(0.3, nsteps.size), label="Bern(0.3) running mean at 1e5")
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
axes[0].semilogx(nsteps, bern_run); axes[0].axhline(0.3, color="k", ls="--")
axes[0].set(title="Bern(0.3): running mean converges", ylim=(0, 0.6))
finals = []
for seed in (1, 2, 3):
    c = np.random.default_rng(seed).standard_cauchy(nsteps.size)
    run = np.cumsum(c) / nsteps; finals.append(run[-1])
    axes[1].semilogx(nsteps, run, lw=0.8, label=f"seed {seed}")
axes[1].set(title="Cauchy: running mean never settles"); axes[1].legend()
fig.tight_layout(); fig.savefig(FIG / "lln_bernoulli_cauchy.png", dpi=120); plt.close(fig)
print("  Cauchy final running means (seeds 1,2,3):", np.round(finals, 3))

# %% [markdown]
# ## Part 9: Linearity needs no independence; variance additivity does

# %%
print("Part 9: linearity on dependent variables")
vals = np.array([-1, 0, 1]); pX = np.full(3, Fr(1, 3))
EX = sum(pX * vals); EY = sum(pX * vals ** 2); EXY = sum(pX * vals ** 3)
assert EX == 0 and EY == Fr(2, 3) and sum(pX * (vals + vals ** 2)) == EX + EY and EXY - EX * EY == 0
assert pX[1] != pX[1] * pX[1]   # P(X=0, Y=0) = P(X=0) = 1/3, but P(X=0)P(Y=0) = 1/9: dependent
xs_ = rng.choice(vals, 100_000); check((xs_ + xs_ ** 2).mean(), 2 / 3, (xs_ + xs_ ** 2).std() / math.sqrt(1e5),
                                       label="E[X + X^2]")
heads = rng.binomial(10, 0.5, 100_000); tails = 10 - heads
assert np.var(heads + tails) == 0 and abs(np.var(heads) - 2.5) < 0.05 and abs(np.cov(heads, tails)[0, 1] + 2.5) < 0.05

# %% [markdown]
# ## Part 10: Numbers in the exercise solutions

# %%
print("Part 10: exercise numbers")
assert [round(1 - (5 / 6) ** 4, 4), round(1 - (35 / 36) ** 24, 4), round(1 - (35 / 36) ** 25, 4)] == [0.5177, 0.4914, 0.5055]
assert min(n for n in range(1, 60) if 1 - (35 / 36) ** n > 0.5) == 25 and round(1 - math.exp(-2 / 3), 4) == 0.4866
assert round(math.log(2) / math.log(36 / 35), 1) == 24.6 and (5**4, 6**4) == (625, 1296)
assert sum(math.comb(5, k) ** 2 for k in range(6)) == math.comb(10, 5) == 252                          # A2
assert [Fr(150, 200), Fr(72, 80), Fr(78, 120), Fr(72, 150), Fr(8, 50)] == [Fr(3, 4), Fr(9, 10), Fr(13, 20), Fr(12, 25), Fr(4, 25)]
assert Fr(9, 10) * Fr(80, 200) / Fr(3, 4) == Fr(12, 25) and (Fr(72, 200), Fr(2, 5) * Fr(3, 4)) == (Fr(9, 25), Fr(3, 10))  # A3
pmax = [Fr(2 * k - 1, 36) for k in range(1, 7)]                                                         # A4
EXm = sum(k * pk for k, pk in zip(range(1, 7), pmax)); EX2m = sum(k * k * pk for k, pk in zip(range(1, 7), pmax))
assert (EXm, EX2m, EX2m - EXm ** 2) == (Fr(161, 36), Fr(791, 36), Fr(2555, 1296))
assert [round(1 - 0.97725, 5), round(100 + 15 * 2.3263, 2)] == [0.02275, 134.89]                        # A5(a)
assert [round(math.exp(-2), 4), round(1 - math.exp(-0.5), 4), round(10 * math.log(2), 3)] == [0.1353, 0.3935, 6.931]
assert Fr(3, 6) - Fr(3, 36) + Fr(1, 216) == 1 - Fr(5, 6) ** 3 == Fr(91, 216)                             # B1(d)
assert (h * h / PH1, h * h / (1 - PH1), PHH / PH1) == (Fr(2, 5), Fr(2, 3), Fr(13, 20))                   # B2
prior = Fr(2, 5); assert h * prior / (h * prior + b * (1 - prior)) == h ** 3 / PHH == Fr(4, 13)
assert round(1 / (1 + 1.5 ** 5), 4) == 0.1164
assert Fr(1, 6) == Fr(1, 2) * Fr(1, 2) * Fr(2, 3) and Fr(1, 6) != Fr(1, 4)                               # B3(b)
pg = 0.2                                                                                                  # B4(d)
assert [round((1 - pg) / pg, 6), round((1 - pg) / pg ** 2, 6), round(math.sqrt(20), 4), round(0.8 ** 10, 4)] == [4, 20, 4.4721, 0.1074]
cov_I = Fr(4, 52) * Fr(3, 51) - Fr(4, 52) ** 2                                                            # B5
assert cov_I == -Fr(1, 13) * Fr(12, 13) / 51 == Fr(-12, 8619)
var_aces = 5 * Fr(1, 13) * Fr(12, 13) * Fr(47, 51)
assert var_aces == Fr(2820, 8619) and round(float(var_aces), 4) == 0.3272 and round(60 / 169, 4) == 0.3550
p_ace = 1 - Fr(math.comb(48, 5), math.comb(52, 5)); assert (math.comb(48, 5), math.comb(52, 5)) == (1_712_304, 2_598_960)
assert round(float(p_ace), 4) == 0.3412 and round(1 - (12 / 13) ** 5, 4) == 0.3298 and round(5 / 13, 4) == 0.3846
assert np.allclose(stats.hypergeom(52, 4, 5).pmf(range(6)), [hgeom_pmf(k, 4, 48, 5) for k in range(6)])
g32 = stats.gamma(3, scale=1 / 2)                                                                         # B6(e)
assert math.isclose(g32.mean(), 1.5) and math.isclose(g32.var(), 0.75) and math.isclose(g32.moment(3), 7.5)


def cauchy_sample(size, rng):
    """B7(a): F^{-1}(u) = tan(pi (u - 1/2)) applied to U ~ Unif(0, 1)."""
    return np.tan(np.pi * (rng.random(size) - 0.5))


xc = cauchy_sample(100_000, rng)
assert stats.kstest(xc, stats.cauchy.cdf).pvalue > 1e-3
assert abs(np.median(xc)) < 0.03 and abs(np.quantile(xc, 0.75) - 1) < 0.04
print(f"  Cauchy sampler: median {np.median(xc):.4f} (target 0), upper quartile {np.quantile(xc, .75):.4f} (target 1)")
mu_, s_ = 1, 2                                                                                            # B8(b)
assert [mu_**2 + s_**2, mu_**3 + 3 * mu_ * s_**2, mu_**4 + 6 * mu_**2 * s_**2 + 3 * s_**4] == [5, 13, 73]
assert np.allclose([stats.norm(1, 2).moment(j) for j in (2, 3, 4)], [5, 13, 73])
days = [(s, d) for s in "BG" for d in range(7)]                                                           # C1
fam = list(itertools.product(days, repeat=2)); tb = lambda c: c == ("B", 2)
assert len(fam) == 196 and cond(lambda f: f[0][0] == f[1][0] == "B", lambda f: tb(f[0]) or tb(f[1]), fam) == Fr(13, 27)
assert all((2 - rr) / (4 - rr) == ((2 * rr - rr**2) / 4) / (rr - rr**2 / 4) for rr in (Fr(1, 7), Fr(1, 2), Fr(1, 365)))
sx = rng.integers(0, 2, (200_000, 2)); dy = rng.integers(0, 7, (200_000, 2)); tbs = (sx == 0) & (dy == 2)
sel = tbs.any(axis=1); check((sx[sel] == 0).all(axis=1).mean(), 13 / 27, prop_se(13 / 27, sel.sum()), label="C1 Tuesday boy")


def wait_for(pattern, p, runs, rng, L=200):
    """C2: tosses until `pattern` (e.g. 'HH') first completes; H with probability p."""
    tosses = rng.random((runs, L)) < p
    a, b2 = (pattern[0] == "H"), (pattern[1] == "H")
    hit = (tosses[:, :-1] == a) & (tosses[:, 1:] == b2)
    assert hit.any(axis=1).all()
    return hit.argmax(axis=1) + 2


for pat, pc, ex in [("HT", 0.5, 4), ("HH", 0.5, 6), ("HT", 0.6, 1 / 0.6 + 1 / 0.4), ("HH", 0.6, 1 / 0.6 + 1 / 0.36)]:
    wv = wait_for(pat, pc, 100_000, rng)
    check(wv.mean(), ex, wv.std() / math.sqrt(wv.size), k=4, label=f"C2 E[W_{pat}], p={pc}")
assert [round(1 / 0.6 + 1 / 0.4, 4), round(1 / 0.6 + 1 / 0.36, 4)] == [4.1667, 4.4444]
fl = np.floor(rng.exponential(1 / math.log(2), 100_000))                                                 # C3
check(fl.mean(), 1.0, fl.std() / math.sqrt(fl.size), k=4, label="C3 E[floor(T)], lambda = ln 2")
check((fl >= 3).mean(), 1 / 8, prop_se(1 / 8, fl.size), k=4, label="C3 P(floor(T) >= 3)")
for nperm in (4, 5, 6):                                                                                   # C4
    cnt = np.array([sum(s[i] == i for i in range(nperm)) for s in itertools.permutations(range(nperm))])
    assert Fr(int(cnt.sum()), cnt.size) == 1 and Fr(int((cnt ** 2).sum()), cnt.size) == 2
    if nperm == 4:
        assert np.bincount(cnt, minlength=5).tolist() == [9, 8, 6, 0, 1]
pk10 = [sum((-1) ** i / math.factorial(i) for i in range(10 - k + 1)) / math.factorial(k) for k in range(3)]
print("  C4 n=10: P(X=k), k=0,1,2:", np.round(pk10, 7), " Pois(1):", np.round([1 / E, 1 / E, 0.5 / E], 7))
assert np.allclose(pk10, [1 / E, 1 / E, 0.5 / E], atol=3e-6)

# %% [markdown]
# ## Your turn
#
# 1. **Two children by simulation.** Simulate 10^5 families. Estimate
#    P(BB | at least one boy) and P(BB | elder is a boy); then simulate
#    "meet a random child, it is a boy". Recover 1/3, 1/2, 1/2 (Trap 3).
# 2. **Near birthdays.** For 23 people, estimate P(some pair has birthdays
#    within one day of each other), circularly (Dec 31 and Jan 1 count).
#    Before simulating, guess using the indicator trick: each pair "nearly
#    matches" with probability 3/365, so the expected count is 253 * 3/365.
# 3. **Coupon collector.** Rolls of a die until all six faces appear is a sum
#    of independent First Success waits with p = 6/6, 5/6, ..., 1/6, so the
#    mean is 6(1 + 1/2 + ... + 1/6). Simulate 10^5 runs and compare with the
#    value asserted below. Then do n = 365 faces and compare with n ln n.
# 4. (Exercise B7) Replace `cauchy_sample` with the ray story of B7(d):
#    X = tan(Theta), Theta ~ Unif(-pi/2, pi/2). Why is this the same code?

# %%
assert round(6 * sum(1 / j for j in range(1, 7)), 1) == 14.7

if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
