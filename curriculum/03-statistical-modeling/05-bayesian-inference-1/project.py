# %% [markdown]
# # ST-5 Project: Bayesian Inference I
#
# Conjugate updating, predictive distributions, priors, credible intervals,
# Bayesian regression, the Laplace approximation, and model comparison. Every
# number quoted in notes.md Section 4 is asserted or printed below, and each
# derived formula is checked against Monte Carlo, quadrature, or a grid.
# Parts: 1 Beta-Binomial and Normal-Normal updating; 2 predictive vs plug-in;
# 3 credible-interval coverage; 4 priors and reparametrization (SymPy);
# 5 Normal-Inverse-Gamma and the t; 6 exponential-family conjugacy;
# 7 decisions and MAP = ridge/lasso; 8 Bayesian regression; 9 Laplace for
# logistic regression; 10 evidence, Occam, BIC; 11 Jeffreys-Lindley;
# 12 posterior predictive checks; 13 stopping rules. Then "Your turn".
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats, special, optimize, integrate
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)


def close(a, b, tol=1e-9):
    return bool(np.allclose(a, b, rtol=tol, atol=tol))


# %% [markdown]
# ## Part 1: Conjugate updating as data arrives (notes 3.2, 3.3; Example 4.1, 4.2)

# %%
def beta_update(a, b, y, n):
    """Beta(a,b) prior, y successes in n trials -> Beta(a+y, b+n-y), eq. (3.4)."""
    return a + y, b + n - y


def normal_update(mu0, tau2, ybar, n, s2):
    """N(mu0,tau2) prior, n obs with mean ybar and known variance s2: precisions add, eq. (3.10)."""
    prec = 1 / tau2 + n / s2
    return (mu0 / tau2 + n * ybar / s2) / prec, 1 / prec


def hpd_beta(a, b, level=0.95):
    """Highest-posterior-density interval of a unimodal Beta: the shortest interval with mass `level`."""
    d = stats.beta(a, b)
    r = optimize.minimize_scalar(lambda p: d.ppf(p + level) - d.ppf(p), bounds=(0, 1 - level),
                                 method="bounded", options={"xatol": 1e-12})
    return d.ppf(r.x), d.ppf(r.x + level)


# Example 4.1: Beta(2,2) prior, 7 successes in 10 trials.
a0 = b0 = 2.0
y, n = 7, 10
a1, b1 = beta_update(a0, b0, y, n)
post = stats.beta(a1, b1)
assert (a1, b1) == (9, 5)
w = n / (a0 + b0 + n)                                            # data weight in (3.6)
assert close(post.mean(), w * y / n + (1 - w) * a0 / (a0 + b0))  # shrinkage form (3.6)
assert close(post.var(), a1 * b1 / ((a1 + b1) ** 2 * (a1 + b1 + 1)))
eq_ci, hpd = post.interval(0.95), hpd_beta(a1, b1)
assert close(post.pdf(hpd[0]), post.pdf(hpd[1]), 1e-5)           # HPD: equal density at both ends
log_ev = (special.gammaln(n + 1) - special.gammaln(y + 1) - special.gammaln(n - y + 1)
          + special.betaln(a1, b1) - special.betaln(a0, b0))     # evidence (3.5)
ev_quad = integrate.quad(lambda p: stats.binom.pmf(y, n, p) * stats.beta.pdf(p, a0, b0), 0, 1)[0]
assert close(np.exp(log_ev), ev_quad, 1e-8)
assert all(close(stats.betabinom.pmf(k, n, 1, 1), 1 / (n + 1)) for k in range(n + 1))  # uniform prior: 1/(n+1)
assert close(integrate.quad(lambda t: t * post.pdf(t), 0, 1)[0], a1 / (a1 + b1))      # P(next success)=posterior mean
assert close(beta_update(1, 1, 7, 10)[0] / sum(beta_update(1, 1, 7, 10)), 8 / 12)    # Laplace's rule (y+1)/(n+2)
print(f"Ex 4.1: posterior Beta({a1:.0f},{b1:.0f}), mean {post.mean():.4f}, sd {post.std():.4f}, w={w:.4f}")
print(f"  95% equal-tailed ({eq_ci[0]:.4f}, {eq_ci[1]:.4f}); HPD ({hpd[0]:.4f}, {hpd[1]:.4f}); "
      f"evidence p(y=7)={np.exp(log_ev):.5f} (uniform prior: {1 / 11:.5f}); P(next success)={a1 / (a1 + b1):.4f}")

# Example 4.2: Normal mean, known sigma^2 = 4, prior N(0, 1), eight observations with mean 1.5.
yN = np.array([0.5, 3.1, 1.2, 2.4, -0.3, 2.9, 0.8, 1.4])
s2, mu0, tau02 = 4.0, 0.0, 1.0
assert close(yN.mean(), 1.5)
mun, taun2 = normal_update(mu0, tau02, yN.mean(), len(yN), s2)
assert close(1 / taun2, 1 / tau02 + len(yN) / s2) and close(mun, 1.0) and close(taun2, 1 / 3)  # (3.10)
print(f"Ex 4.2: posterior N({mun:.4f}, {taun2:.4f}); predictive variance {s2 + taun2:.4f}; "
      f"95% interval ({mun - 1.96 * taun2 ** .5:.4f}, {mun + 1.96 * taun2 ** .5:.4f})")

# Figure grid: posteriors as data arrive; sequential updating equals batch updating.
p_true, mu_true = 0.7, 1.5
bern = (rng.random(100) < p_true).astype(int)
norm = rng.normal(mu_true, np.sqrt(s2), 100)
marks = [0, 1, 3, 10, 30, 100]
fig, axs = plt.subplots(2, len(marks), figsize=(16, 5.5), sharey="row")
pg, mg = np.linspace(0, 1, 400), np.linspace(-3, 4, 400)
sa, sb, sm, st = 2.0, 2.0, 0.0, 1.0
for i in range(101):
    if i in marks:
        k = marks.index(i)
        assert close((sa, sb), beta_update(2, 2, bern[:i].sum(), i))                      # sequential = batch
        assert close((sm, st), normal_update(0, 1, norm[:i].mean() if i else 0.0, i, s2))
        axs[0, k].plot(pg, stats.beta.pdf(pg, sa, sb)); axs[0, k].axvline(p_true, c="k", ls=":")
        axs[0, k].set_title(f"n={i}: Beta({sa:.0f},{sb:.0f})", fontsize=9)
        axs[1, k].plot(mg, stats.norm.pdf(mg, sm, np.sqrt(st))); axs[1, k].axvline(mu_true, c="k", ls=":")
        axs[1, k].set_title(f"n={i}: N({sm:.2f},{st:.3f})", fontsize=9)
    if i < 100:
        sa, sb = beta_update(sa, sb, bern[i], 1)
        sm, st = normal_update(sm, st, norm[i], 1, s2)
fig.suptitle("Posteriors sharpen as data arrive (dotted = truth). Sequential and batch updates agree.")
fig.tight_layout(); fig.savefig(FIG / "p1_updating.png", dpi=110); plt.close(fig)
print("Part 1: sequential updating equals batch updating at every checkpoint.")

# %% [markdown]
# ## Part 2: Posterior predictive vs plug-in (notes 3.2, 3.3, 3.10)

# %%
m = 10                                                           # ten future trials
pred = stats.betabinom(m, a1, b1)                                # Beta-Binomial predictive (3.8)
theta = post.rvs(400_000, random_state=rng)
ytil = rng.binomial(m, theta)                                    # sample theta, then y~
emp = np.bincount(ytil, minlength=m + 1) / ytil.size
assert np.abs(emp - pred.pmf(np.arange(m + 1))).max() < 0.004
assert close(pred.var(), m * a1 * b1 * (a1 + b1 + m) / ((a1 + b1) ** 2 * (a1 + b1 + 1)))
plug_mle, plug_mean = stats.binom(m, y / n), stats.binom(m, a1 / (a1 + b1))
assert pred.var() > plug_mean.var() > plug_mle.var()
print(f"Part 2: predictive var of 10 future trials {pred.var():.4f} (MC {ytil.var():.4f}) vs plug-in "
      f"{plug_mle.var():.3f} (MLE), {plug_mean.var():.4f} (posterior mean)")
muT = rng.normal(mun, np.sqrt(taun2), 400_000)                   # Normal: theta then y~ (3.11)
ynew = rng.normal(muT, np.sqrt(s2))
assert abs(ynew.var() - (s2 + taun2)) < 0.03 and abs(ynew.mean() - mun) < 0.01
fig, ax = plt.subplots(figsize=(6, 3.5))
ks = np.arange(m + 1)
ax.bar(ks - .2, pred.pmf(ks), .4, label="posterior predictive"); ax.bar(ks + .2, plug_mle.pmf(ks), .4, label="plug-in MLE")
ax.legend(); ax.set_title("10 future trials: plug-in is too confident")
fig.tight_layout(); fig.savefig(FIG / "p2_predictive.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 3: Credible intervals and coverage (notes 3.8)

# %%
def eq_tail(a, b, level=0.95):
    lo = (1 - level) / 2
    return stats.beta.ppf(lo, a, b), stats.beta.ppf(1 - lo, a, b)


# (a) Bayesian calibration (3.25): theta from the prior, data from theta -> 95% interval covers 95%.
N, nn = 400_000, 20
th = rng.beta(2, 2, N); yy = rng.binomial(nn, th)
lo, hi = eq_tail(2 + yy, 2 + nn - yy)
cov = np.mean((lo <= th) & (th <= hi))
assert abs(cov - 0.95) < 4 * np.sqrt(0.95 * 0.05 / N)
print(f"Part 3a: theta~Beta(2,2), n=20: coverage over prior and data {cov:.4f}")


# (b) Frequentist coverage at a fixed theta: exact sum over y (no Monte Carlo).
def freq_coverage(theta, n, a, b):
    ys = np.arange(n + 1)
    lo, hi = eq_tail(a + ys, b + n - ys)
    return float(np.sum(stats.binom.pmf(ys, n, theta) * ((lo <= theta) & (theta <= hi))))


ns = [5, 10, 20, 50, 100, 200, 400, 1000]
priors = {"Jeffreys Beta(.5,.5)": (.5, .5), "uniform Beta(1,1)": (1, 1), "strong Beta(10,10)": (10, 10)}
fig, axs = plt.subplots(1, 2, figsize=(10, 3.5), sharey=True)
cov_tab = {}
for ax, t0 in zip(axs, [0.5, 0.1]):
    for name, (a, b) in priors.items():
        c = [freq_coverage(t0, k, a, b) for k in ns]
        cov_tab[(name, t0)] = c
        ax.semilogx(ns, c, "o-", label=name)
    ax.axhline(.95, c="k", ls=":"); ax.set_title(f"true theta = {t0}"); ax.set_xlabel("n")
axs[0].set_ylabel("exact frequentist coverage"); axs[0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "p3_coverage.png", dpi=110); plt.close(fig)
for t0 in [0.5, 0.1]:
    print(f"Part 3b: theta={t0} coverage at n=10,100,1000:",
          {k[0][:7]: [round(v[ns.index(j)], 3) for j in (10, 100, 1000)] for k, v in cov_tab.items() if k[1] == t0})
assert cov_tab[("strong Beta(10,10)", 0.1)][1] < 0.6                    # strong prior, wrong place, n=10
assert cov_tab[("strong Beta(10,10)", 0.5)][1] > 0.97                   # strong prior, right place: over-covers
assert all(abs(cov_tab[(k, 0.1)][-1] - .95) < .03 for k in ("Jeffreys Beta(.5,.5)", "uniform Beta(1,1)"))   # weak priors: agree at n=1000
sc_ = cov_tab[("strong Beta(10,10)", 0.1)]
assert sc_[1] < sc_[3] < sc_[5] < sc_[-1] < 0.93                        # the strong prior's damage fades slowly with n
tg = np.linspace(0.0005, 0.9995, 4000)                                  # average coverage over a uniform theta grid
avg_cov = np.mean([freq_coverage(t, 20, 1, 1) for t in tg])
assert abs(avg_cov - 0.95) < 0.003                                      # calibration = coverage averaged over the prior
print(f"  uniform prior, n=20: coverage averaged over theta~U(0,1) = {avg_cov:.4f}")

# %% [markdown]
# ## Part 4: Priors, Jeffreys, and reparametrization (notes 3.7; Example 4.5)

# %%
p_, lam, mu_, sg = sp.symbols("p lambda mu sigma", positive=True)
yv = sp.symbols("y", real=True)


def fisher(logpmf, theta, repl):
    """I(theta) = -E[d^2 log p / d theta^2] with y-dependence replaced by its expectation (ST-1 (3.11))."""
    return sp.simplify(-sp.diff(logpmf, theta, 2).subs(repl))


rr = sp.symbols("r", positive=True)                                    # r stands for (y - mu)^2, whose mean is sigma^2
I_bern = fisher(yv * sp.log(p_) + (1 - yv) * sp.log(1 - p_), p_, {yv: p_})
I_pois = fisher(yv * sp.log(lam) - lam, lam, {yv: lam})
I_mu = fisher(-(yv - mu_) ** 2 / (2 * sg ** 2), mu_, {yv: mu_})
I_sig = fisher(-rr / (2 * sg ** 2) - sp.log(sg), sg, {rr: sg ** 2})
assert sp.simplify(I_bern - 1 / (p_ * (1 - p_))) == 0 and sp.simplify(I_pois - 1 / lam) == 0
assert sp.simplify(I_mu - 1 / sg ** 2) == 0 and sp.simplify(I_sig - 2 / sg ** 2) == 0
print("Part 4: Fisher information  Bernoulli", I_bern, "| Poisson", I_pois, "| Normal mean", I_mu, "| Normal sigma", I_sig)

eta = np.linspace(-8, 8, 801)
pp = special.expit(eta)
flat_eta = stats.beta.pdf(pp, 1, 1) * pp * (1 - pp)                     # flat in p, pushed to logit: (3.22)
assert close(flat_eta, stats.logistic.pdf(eta))
jeff_eta = stats.beta.pdf(pp, .5, .5) * pp * (1 - pp)                   # Jeffreys in p, pushed to logit
assert close(jeff_eta, 1 / (2 * np.pi * np.cosh(eta / 2)))
assert close(jeff_eta, np.sqrt(pp * (1 - pp)) / np.pi)                  # = sqrt(I_eta)/pi with I_eta = p(1-p): invariance (3.24)
ew = np.linspace(-80, 80, 160001)                                     # tails decay like e^(-|eta|/2): integrate wide
assert close(np.trapezoid(1 / (2 * np.pi * np.cosh(ew / 2)), ew), 1, 1e-6)
# The mode is not invariant: mode of Beta(9,5) in p vs the mode of its density on the logit scale.
mode_p = (a1 - 1) / (a1 + b1 - 2)
r = optimize.minimize_scalar(lambda e: -(stats.beta.pdf(special.expit(e), a1, b1) * special.expit(e) * (1 - special.expit(e))),
                             bounds=(-5, 5), method="bounded", options={"xatol": 1e-10})
mode_eta_as_p = special.expit(r.x)
assert close(mode_p, 2 / 3) and close(mode_eta_as_p, 9 / 14, 1e-6)
print(f"  Beta(9,5): mode in p {mode_p:.4f}; mode on logit scale maps back to p={mode_eta_as_p:.4f}; mean {a1 / (a1 + b1):.4f}")
# Example 4.5 numbers: no events in 10 trials; zero counts in ten Poisson intervals.
fl, jf = stats.beta(1, 11), stats.beta(.5, 10.5)
print(f"Ex 4.5: Bernoulli y=0,n=10: uniform 95% HPD upper {fl.ppf(.95):.4f} (one-sided), mean {fl.mean():.4f}; "
      f"Jeffreys mean {jf.mean():.4f}, 95% eq-tailed ({jf.ppf(.025):.5f}, {jf.ppf(.975):.4f})")
# Prior predictive check: a "vague" N(0, 10^2) on the log-odds says p is almost surely near 0 or 1.
for sd in (10, 1.5):
    pr = special.expit(rng.normal(0, sd, 400_000))
    extreme = np.mean((pr < .01) | (pr > .99))
    assert abs(extreme - 2 * stats.norm.cdf(-special.logit(.99) / sd)) < 0.004
    print(f"  prior N(0,{sd}^2) on log-odds: P(p<.01 or p>.99) = {extreme:.4f}")
fig, axs = plt.subplots(1, 3, figsize=(13, 3.5))
pl = np.linspace(0.001, .999, 500)
axs[0].plot(pl, stats.beta.pdf(pl, 1, 1), label="flat"); axs[0].plot(pl, stats.beta.pdf(pl, .5, .5), label="Jeffreys"); axs[0].set_ylim(0, 4)
axs[0].set_title("priors on p"); axs[0].legend()
axs[1].plot(eta, flat_eta, label="induced by flat p (logistic)"); axs[1].plot(eta, jeff_eta, label="induced by Jeffreys p")
axs[1].set_title("the same priors on logit(p)"); axs[1].legend(fontsize=8)
for sd in (10, 1.5):
    axs[2].hist(special.expit(rng.normal(0, sd, 20000)), 40, alpha=.5, density=True, label=f"N(0,{sd}^2) on log-odds")
axs[2].set_title("prior predictive: implied p"); axs[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "p4_priors.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 5: Normal-Inverse-Gamma, the t marginal, and ST-2's interval (notes 3.4; Example 4.3)

# %%
g = np.array([9.79, 9.83, 9.76, 9.85, 9.81, 9.78, 9.88, 9.80, 9.74, 9.82])   # ST-2 Example 4.2


def nig_update(mu0, k0, a0, b0, y):
    """NIG hyperparameter update, eq. (3.13): returns (mu_n, kappa_n, alpha_n, beta_n)."""
    y = np.asarray(y, float); n = y.size; yb = y.mean(); S = ((y - yb) ** 2).sum()
    kn = k0 + n
    return (k0 * mu0 + n * yb) / kn, kn, a0 + n / 2, b0 + 0.5 * (S + k0 * n * (yb - mu0) ** 2 / kn)


def nig_t(mun, kn, an, bn, predictive=False):
    """Marginal t for mu (3.14) or for a new observation (3.15)."""
    sc2 = bn / (an * kn) * ((kn + 1) if predictive else 1)
    return stats.t(df=2 * an, loc=mun, scale=np.sqrt(sc2))


n_g, S_g = g.size, ((g - g.mean()) ** 2).sum()
# Reference prior p(mu,sigma^2) ~ 1/sigma^2 is the formal limit k0->0, a0->-1/2, b0->0 of (3.13).
ref = nig_update(0.0, 0.0, -0.5, 0.0, g)
assert close(ref[:2], (g.mean(), n_g)) and close(ref[2], (n_g - 1) / 2) and close(ref[3], S_g / 2)
ref_ci = nig_t(*ref).interval(0.95)
t_ci = g.mean() + np.array([-1, 1]) * stats.t.ppf(.975, n_g - 1) * g.std(ddof=1) / np.sqrt(n_g)   # ST-2 interval
assert close(ref_ci, t_ci, 1e-10) and close(np.round(t_ci, 4), [9.7762, 9.8358])
sig2_ci = (stats.invgamma(ref[2], scale=ref[3]).ppf(.025), stats.invgamma(ref[2], scale=ref[3]).ppf(.975))
assert close(sig2_ci, (S_g / stats.chi2.ppf(.975, 9), S_g / stats.chi2.ppf(.025, 9)))               # = ST-2 (3.4)
# A proper weakly informative prior.
prior = (9.5, 1.0, 2.0, 0.02)
post_n = nig_update(*prior, g)
seq = prior
for gi in g:
    seq = nig_update(*seq, [gi])
assert close(seq, post_n)                                                                        # sequential = batch
tm, tp = nig_t(*post_n), nig_t(*post_n, predictive=True)
# Check the t marginal against brute-force integration of likelihood x prior over a (mu, sigma^2) grid.
mus, s2s = np.linspace(9.45, 10.15, 1401), np.geomspace(1e-4, 0.05, 4000)
M, V = np.meshgrid(mus, s2s, indexing="ij")
logp = (-0.5 * n_g * np.log(2 * np.pi * V) - (S_g + n_g * (g.mean() - M) ** 2) / (2 * V)
        + stats.norm.logpdf(M, prior[0], np.sqrt(V / prior[1])) + stats.invgamma.logpdf(V, prior[2], scale=prior[3]))
dens = np.exp(logp - logp.max())
marg = np.trapezoid(dens, s2s, axis=1)
marg /= np.trapezoid(marg, mus)
assert np.abs(marg - tm.pdf(mus)).max() < 0.003 * tm.pdf(mus).max()
# Joint Jeffreys prior ~ sigma^-2 in (mu, sigma): t with n (not n-1) df.
jt = stats.t(n_g, g.mean(), np.sqrt(S_g) / n_g)
num = [integrate.quad(lambda s: s ** (-n_g - 2) * np.exp(-(S_g + n_g * (mu - g.mean()) ** 2) / (2 * s * s)), 1e-4, 5)[0]
       for mu in (g.mean(), g.mean() + .03, g.mean() - .05)]
assert close(np.array(num) / num[0], jt.pdf([g.mean(), g.mean() + .03, g.mean() - .05]) / jt.pdf(g.mean()), 1e-5)
print(f"Ex 4.3: reference-prior 95% interval ({ref_ci[0]:.4f}, {ref_ci[1]:.4f}) = ST-2 t interval; sigma^2 interval "
      f"({sig2_ci[0]:.5f}, {sig2_ci[1]:.5f}), sigma ({sig2_ci[0] ** .5:.4f}, {sig2_ci[1] ** .5:.4f})")
print(f"  proper NIG{prior}: posterior (mu_n,k_n,a_n,b_n)=({post_n[0]:.4f},{post_n[1]:.0f},{post_n[2]:.0f},{post_n[3]:.6f}); "
      f"mu 95% ({tm.interval(.95)[0]:.4f}, {tm.interval(.95)[1]:.4f}); predictive 95% ({tp.interval(.95)[0]:.4f}, {tp.interval(.95)[1]:.4f})")
print(f"  joint-Jeffreys (n df) half-width {stats.t.ppf(.975, n_g) * np.sqrt(S_g) / n_g:.4f} vs reference {(ref_ci[1] - ref_ci[0]) / 2:.4f}")
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(mus, marg, label="grid (likelihood x prior)"); ax.plot(mus, tm.pdf(mus), "--", label="closed-form t (3.14)")
ax.plot(mus, nig_t(*ref).pdf(mus), ":", label="reference prior = ST-2 t"); ax.legend(fontsize=8); ax.set_title("posterior of mu")
fig.tight_layout(); fig.savefig(FIG / "p5_nig.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 6: Exponential-family conjugacy reproduces the three cases (notes 3.5, 3.6; Example 4.4)

# %%
def ef_update(chi, nu, T_sum, n):
    """Conjugate update for prior exp(eta.chi - nu A(eta)): chi -> chi + sum T, nu -> nu + n, eq. (3.20)."""
    return chi + T_sum, nu + n


# Poisson-Gamma: (chi, nu) = (alpha0, beta0). Example 4.4.
yP = np.array([3, 5, 2, 4, 6])
chi, nu = ef_update(2.0, 1.0, yP.sum(), yP.size)
assert (chi, nu) == (22.0, 6.0)
lamg = np.linspace(1e-3, 12, 4000)
lp = np.log(stats.gamma.pdf(lamg, 2, scale=1)) + stats.poisson.logpmf(yP[:, None], lamg).sum(0)
gp = np.exp(lp - lp.max()); gp /= np.trapezoid(gp, lamg)
assert np.abs(gp - stats.gamma.pdf(lamg, chi, scale=1 / nu)).max() < 1e-4        # grid posterior = Gamma(22, 6)
# Prior in eta = log(lambda): exp(eta chi - nu e^eta) has normalizer Gamma(chi)/nu^chi; E[A'(eta)] = chi/nu (3.21).
Zq = integrate.quad(lambda e: np.exp(e * 2 - 1 * np.exp(e)), -40, 6)[0]
assert close(Zq, special.gamma(2) / 1 ** 2, 1e-6)
mean_param = integrate.quad(lambda e: np.exp(e) * np.exp(e * chi - nu * np.exp(e)), -40, 6)[0] / \
    integrate.quad(lambda e: np.exp(e * chi - nu * np.exp(e)), -40, 6)[0]
assert close(mean_param, chi / nu, 1e-6)
nb = stats.nbinom(chi, nu / (nu + 1))                                               # predictive (3.18)
assert close(nb.mean(), chi / nu) and close(nb.var(), chi * (nu + 1) / nu ** 2)
mc = rng.poisson(rng.gamma(chi, 1 / nu, 400_000))
assert abs(mc.mean() - nb.mean()) < 0.02 and abs(mc.var() - nb.var()) < 0.06
assert np.abs(np.bincount(mc, minlength=16)[:16] / mc.size - nb.pmf(np.arange(16))).max() < 0.002
print(f"Ex 4.4: (chi,nu)=({chi:.0f},{nu:.0f}) -> Gamma({chi:.0f},{nu:.0f}), mean {chi / nu:.4f} = (1*2+5*4)/6; "
      f"NB predictive mean {nb.mean():.4f}, var {nb.var():.4f} (Poisson plug-in var {chi / nu:.4f}); P(y~=0)={nb.pmf(0):.4f}")
# Bernoulli: (chi, nu) = (alpha0, alpha0+beta0). Check the Jacobian computation by quadrature.
Zb = integrate.quad(lambda e: np.exp(e * 2.0 - 4.0 * np.log1p(np.exp(e))), -40, 40)[0]
assert close(Zb, special.beta(2, 2), 1e-6)                                          # chi=2, nu=4 <-> Beta(2, 2)
# Normal, known sigma: nu = sigma^2/tau0^2, chi = nu*mu0 -> the same posterior as (3.10).
chiN, nuN = ef_update(s2 / tau02 * mu0, s2 / tau02, yN.sum(), yN.size)
assert close(chiN / nuN, mun) and close(s2 / nuN, taun2)
print("Part 6: Beta, Gamma and Normal updates agree with (chi+sum T, nu+n).")

# %% [markdown]
# ## Part 7: Decisions, MAP, ridge and lasso (notes 3.9; Example 4.6)

# %%
d = stats.beta(a1, b1)
risk = lambda act, loss: integrate.quad(lambda t: loss(act, t) * d.pdf(t), 0, 1)[0]
best_sq = optimize.minimize_scalar(lambda a: risk(a, lambda a, t: (a - t) ** 2), bounds=(0, 1), method="bounded", options={"xatol": 1e-9}).x
best_ab = optimize.minimize_scalar(lambda a: risk(a, lambda a, t: abs(a - t)), bounds=(0, 1), method="bounded", options={"xatol": 1e-9}).x
assert abs(best_sq - d.mean()) < 1e-5 and abs(best_ab - d.median()) < 1e-4
print(f"Part 7: Beta(9,5): squared-loss action {best_sq:.4f} = mean {d.mean():.4f}; abs-loss action {best_ab:.4f} = median {d.median():.4f}; "
      f"MAP {mode_p:.4f}")
# Example 4.6: slope-only regression (ST-3 five points, centered), sigma^2 = 0.8, prior N(0, tau^2 = 0.5): lambda = sigma^2/tau^2.
xc, yc = np.array([-2., -1, 0, 1, 2]), np.array([-2., 0, 1, 0, 1])
Sxx, Sxy, sig2, tau2 = xc @ xc, xc @ yc, 0.8, 0.5
lam_r = sig2 / tau2
map_slope = optimize.minimize_scalar(lambda b: ((yc - b * xc) ** 2).sum() / (2 * sig2) + b * b / (2 * tau2)).x
assert close(map_slope, Sxy / (Sxx + lam_r), 1e-6) and close(Sxx, 10) and close(Sxy, 6)
print(f"Ex 4.6: OLS slope {Sxy / Sxx:.4f}, lambda={lam_r:.2f}, MAP = ridge slope {Sxy / (Sxx + lam_r):.4f} (factor {Sxx / (Sxx + lam_r):.4f}); "
      f"posterior var {sig2 / (Sxx + lam_r):.5f} vs OLS var {sig2 / Sxx:.4f}")
# Lasso = MAP under a Laplace prior: orthonormal case y|b ~ N(b, 1), prior exp(-|b|/s), lambda = 1/s.
soft = lambda z, l: np.sign(z) * max(abs(z) - l, 0)
for z, sc in [(1.2, 2.0), (0.3, 2.0)]:
    nlp = lambda b: 0.5 * (z - b) ** 2 + abs(b) / sc
    b_map = optimize.minimize_scalar(nlp, bounds=(-3, 3), method="bounded", options={"xatol": 1e-10}).x
    pm = integrate.quad(lambda b: b * np.exp(-nlp(b)), -12, 12, points=[0])[0] / integrate.quad(lambda b: np.exp(-nlp(b)), -12, 12, points=[0])[0]
    assert abs(b_map - soft(z, 1 / sc)) < 1e-5
    print(f"  Laplace prior scale {sc}, z={z}: MAP {b_map:.4f} = soft-threshold {soft(z, 1 / sc):.4f}; posterior mean {pm:.4f}")

# %% [markdown]
# ## Part 8: Bayesian linear regression and the predictive band (notes 3.11; Example 4.7)

# %%
def blr_posterior(X, y, s2, mu0, S0):
    """beta ~ N(mu0,S0), y|beta ~ N(X beta, s2 I) -> posterior N(mu_n, S_n), eq. (3.31)."""
    P0 = np.linalg.inv(S0)
    Sn = np.linalg.inv(P0 + X.T @ X / s2)
    return Sn @ (P0 @ mu0 + X.T @ y / s2), Sn


def blr_predict(Xs, mun, Sn, s2):
    """Predictive mean x*'mu_n and variance s2 + x*' S_n x*, eq. (3.33)."""
    return Xs @ mun, s2 + np.einsum("ij,jk,ik->i", Xs, Sn, Xs)


xr, yr = np.arange(1., 6), np.array([2., 4, 5, 4, 5])                  # ST-3 Example 4.1
X = np.column_stack([np.ones(5), xr])
s2r, tau2r = 0.8, 10.0
mun_r, Sn_r = blr_posterior(X, yr, s2r, np.zeros(2), tau2r * np.eye(2))
lam2 = s2r / tau2r
ridge = np.linalg.solve(X.T @ X + lam2 * np.eye(2), X.T @ yr)           # ST-3 (3.39), intercept penalized here too
assert close(mun_r, ridge) and close(Sn_r, s2r * np.linalg.inv(X.T @ X + lam2 * np.eye(2)))      # (3.32)
nlpost = lambda b: ((yr - X @ b) ** 2).sum() / (2 * s2r) + b @ b / (2 * tau2r)
assert close(optimize.minimize(nlpost, np.zeros(2), tol=1e-12).x, mun_r, 1e-5)                   # MAP = posterior mean
P0 = np.eye(2) * tau2r                                                  # F5 (3.19) applied to the joint of (beta, y)
cm = np.zeros(2) + P0 @ X.T @ np.linalg.solve(X @ P0 @ X.T + s2r * np.eye(5), yr)
assert close(cm, mun_r)
xs = np.array([[1., 3.], [1., 8.]])
mp, vp = blr_predict(xs, mun_r, Sn_r, s2r)
big_mu, big_S = blr_posterior(X, yr, s2r, np.zeros(2), 1e8 * np.eye(2))                           # flat-prior limit
_, vflat = blr_predict(xs, big_mu, big_S, s2r)
h = 1 / 5 + (xs[:, 1] - 3) ** 2 / 10
assert close(vflat, s2r * (1 + h), 1e-6) and np.allclose(big_mu, [2.2, 0.6], atol=1e-5)          # = ST-3 prediction variance
print(f"Ex 4.7: posterior mean {mun_r.round(4)}, sd {np.sqrt(np.diag(Sn_r)).round(4)}; predictive mean/var at x=3: {mp[0]:.4f}/{vp[0]:.4f}, "
      f"at x=8: {mp[1]:.4f}/{vp[1]:.4f}; flat-prior var {vflat.round(4)} = 0.8(1+h), h={h.round(2)}")
# Unknown sigma^2: beta|s2 ~ N(mu0, s2 V0), s2 ~ IG(a0,b0); predictive is t_{2 a_n}  (3.34).
def nig_reg(X, y, mu0, V0, a0, b0):
    Vn = np.linalg.inv(np.linalg.inv(V0) + X.T @ X)
    mn = Vn @ (np.linalg.solve(V0, mu0) + X.T @ y)
    bn = b0 + 0.5 * (y @ y + mu0 @ np.linalg.solve(V0, mu0) - mn @ np.linalg.solve(Vn, mn))
    return mn, Vn, a0 + len(y) / 2, bn


def nig_reg_logev(X, y, mu0, V0, a0, b0):
    """log p(y) for the NIG regression model, eq. (3.39)."""
    mn, Vn, an, bn = nig_reg(X, y, mu0, V0, a0, b0)
    return (special.gammaln(an) - special.gammaln(a0) + 0.5 * (np.linalg.slogdet(Vn)[1] - np.linalg.slogdet(V0)[1])
            + a0 * np.log(b0) - an * np.log(bn) - 0.5 * len(y) * np.log(2 * np.pi))


V0, a0r, b0r = 10 * np.eye(2), 2.0, 1.0
chain = 0.0                                                             # chain rule: p(y) = prod p(y_i | y_<i)
for i in range(5):
    if i == 0:
        mn_, Vn_, an_, bn_ = np.zeros(2), V0, a0r, b0r
    else:
        mn_, Vn_, an_, bn_ = nig_reg(X[:i], yr[:i], np.zeros(2), V0, a0r, b0r)
    chain += stats.t(2 * an_, X[i] @ mn_, np.sqrt(bn_ / an_ * (1 + X[i] @ Vn_ @ X[i]))).logpdf(yr[i])
assert close(chain, nig_reg_logev(X, yr, np.zeros(2), V0, a0r, b0r), 1e-9)                      # (3.39) = product of t predictives
mn_, Vn_, an_, bn_ = nig_reg(X, yr, np.zeros(2), V0, a0r, b0r)
s2draw = stats.invgamma(an_, scale=bn_).rvs(200_000, random_state=rng)
bdraw = mn_[:, None] + np.linalg.cholesky(Vn_) @ rng.standard_normal((2, 200_000)) * np.sqrt(s2draw)
ydraw = xs[1] @ bdraw + np.sqrt(s2draw) * rng.standard_normal(200_000)
tq = stats.t(2 * an_, xs[1] @ mn_, np.sqrt(bn_ / an_ * (1 + xs[1] @ Vn_ @ xs[1])))
assert abs(np.mean(ydraw <= tq.ppf(.9)) - .9) < .003 and abs(np.mean(ydraw <= tq.ppf(.1)) - .1) < .003
fig, axs = plt.subplots(1, 3, figsize=(15, 4))
b0g, b1g = np.meshgrid(np.linspace(-1, 5, 200), np.linspace(-1, 2, 200))
Z = np.stack([b0g, b1g], -1) - mun_r
Q = np.einsum("...i,ij,...j", Z, np.linalg.inv(Sn_r), Z)
axs[0].contour(b0g, b1g, Q, [2.30, 6.18, 11.8]); axs[0].plot(*mun_r, "k*"); axs[0].set(xlabel="intercept", ylabel="slope", title="posterior contours (68/95/99.7%)")
xg = np.linspace(-2, 12, 200); Xg = np.column_stack([np.ones_like(xg), xg])
for bs in rng.multivariate_normal(mun_r, Sn_r, 50):
    axs[1].plot(xg, Xg @ bs, c="C0", alpha=.15)
axs[1].plot(xr, yr, "ko"); axs[1].set_title("50 posterior-sampled lines")
mg_, vg_ = blr_predict(Xg, mun_r, Sn_r, s2r)
axs[2].fill_between(xg, mg_ - 2 * np.sqrt(vg_), mg_ + 2 * np.sqrt(vg_), alpha=.3); axs[2].plot(xg, mg_); axs[2].plot(xr, yr, "ko")
axs[2].set_title("predictive band widens away from data")
fig.tight_layout(); fig.savefig(FIG / "p8_regression.png", dpi=110); plt.close(fig)
assert vg_[np.argmin(abs(xg - 3))] < vg_[0] and vg_[-1] > vg_[np.argmin(abs(xg - 8))]             # variance grows outward

# %% [markdown]
# ## Part 9: The Laplace approximation for Bayesian logistic regression (notes 3.11; Example 4.8)

# %%
def laplace_logistic(X, y, s0sq, iters=25):
    """MAP by Newton (ST-4 IRLS with a ridge term), Hessian at the MAP, and the Laplace posterior N(beta_hat, (-H)^-1), eq. (3.35)."""
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = special.expit(X @ b)
        g_ = X.T @ (y - p) - b / s0sq
        H = -(X.T * (p * (1 - p))) @ X - np.eye(len(b)) / s0sq
        b = b - np.linalg.solve(H, g_)
    p = special.expit(X @ b)
    H = -(X.T * (p * (1 - p))) @ X - np.eye(len(b)) / s0sq
    return b, -H


def log_joint(X, y, B, s0sq):
    """log p(y|beta) + log N(beta; 0, s0sq I) for the rows of B."""
    eta_ = B @ X.T
    return (y * eta_ - np.logaddexp(0, eta_)).sum(1) + stats.norm.logpdf(B, 0, np.sqrt(s0sq)).sum(1)


# Example 4.8: one slope parameter, six points, N(0, 2^2) prior.
x1, y1 = np.array([-2., -1, .5, 1, 2, 3]), np.array([0., 0, 1, 0, 1, 1])
bh1, A1 = laplace_logistic(x1[:, None], y1, 4.0)
sd_l = 1 / np.sqrt(A1[0, 0])
bg = np.linspace(-3, 8, 20001)
lj = log_joint(x1[:, None], y1, bg[:, None], 4.0)
ev1 = np.trapezoid(np.exp(lj), bg)
pm1 = np.trapezoid(bg * np.exp(lj), bg) / ev1
sd1 = np.sqrt(np.trapezoid((bg - pm1) ** 2 * np.exp(lj), bg) / ev1)
logev_lap = log_joint(x1[:, None], y1, bh1[None], 4.0)[0] + 0.5 * np.log(2 * np.pi) - 0.5 * np.log(A1[0, 0])
print(f"Ex 4.8: MAP {bh1[0]:.4f}, Laplace sd {sd_l:.4f}; exact mean {pm1:.4f}, sd {sd1:.4f}; log evidence Laplace {logev_lap:.4f} vs exact {np.log(ev1):.4f}")
assert abs(logev_lap - np.log(ev1)) < 0.1 and 0.2 < pm1 - bh1[0] < 0.5 and abs(sd_l / sd1 - 1) < 0.15   # evidence good, skew missed
# Two parameters (intercept, slope), n = 60: brute-force grid vs Laplace.
n2 = 60
x2 = rng.normal(size=n2); y2 = (rng.random(n2) < special.expit(-0.5 + 1.5 * x2)).astype(float)
X2 = np.column_stack([np.ones(n2), x2]); s0 = 9.0
bh, A = laplace_logistic(X2, y2, s0)
ga, gb = np.linspace(-2.5, 1.5, 300), np.linspace(-0.5, 4.5, 300)
GA, GB = np.meshgrid(ga, gb, indexing="ij")
LJ = log_joint(X2, y2, np.column_stack([GA.ravel(), GB.ravel()]), s0).reshape(GA.shape)
Wg = np.exp(LJ - LJ.max()); cell = (ga[1] - ga[0]) * (gb[1] - gb[0])
Zg = Wg.sum() * cell; Wn = Wg / Wg.sum()
gm = np.array([(Wn * GA).sum(), (Wn * GB).sum()])
gcov = np.array([[(Wn * (GA - gm[0]) ** 2).sum(), (Wn * (GA - gm[0]) * (GB - gm[1])).sum()],
                 [0, (Wn * (GB - gm[1]) ** 2).sum()]]); gcov[1, 0] = gcov[0, 1]
logev_grid = LJ.max() + np.log(Zg)
logev_lap2 = log_joint(X2, y2, bh[None], s0)[0] + np.log(2 * np.pi) - 0.5 * np.linalg.slogdet(A)[1]
print(f"Part 9: 2-D logistic n=60: MAP {bh.round(3)} vs grid mean {gm.round(3)}; sd Laplace {np.sqrt(np.diag(np.linalg.inv(A))).round(3)} "
      f"vs grid {np.sqrt(np.diag(gcov)).round(3)}; log evidence {logev_lap2:.3f} vs {logev_grid:.3f}")
assert np.all(np.abs(bh - gm) < 0.3 * np.sqrt(np.diag(gcov))) and abs(logev_lap2 - logev_grid) < 0.15
assert np.all(np.abs(np.sqrt(np.diag(np.linalg.inv(A))) / np.sqrt(np.diag(gcov)) - 1) < 0.08)
fig, ax = plt.subplots(figsize=(5, 4.5))
ax.contour(GA, GB, Wg, [np.exp(-.5 * c) for c in (6.2, 2.3)], colors="C0")
Zl = np.stack([GA - bh[0], GB - bh[1]], -1)
ax.contour(GA, GB, np.exp(-.5 * np.einsum("...i,ij,...j", Zl, A, Zl)), [np.exp(-.5 * c) for c in (6.2, 2.3)], colors="C3", linestyles="--")
ax.set(xlabel="intercept", ylabel="slope", title="grid posterior (blue) vs Laplace (red dashed)")
fig.tight_layout(); fig.savefig(FIG / "p9_laplace.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 10: Evidence, Occam factors, BIC (notes 3.12; Example 4.9)

# %%
def gauss_evidence(Phi, y, s2, alpha):
    """log p(y | degree-d model), prior beta ~ N(0, I/alpha), noise variance s2, closed form (3.37).
    Returns (log evidence, max-posterior log likelihood term, log Occam factor)."""
    n_, d_ = Phi.shape
    A_ = alpha * np.eye(d_) + Phi.T @ Phi / s2
    m_ = np.linalg.solve(A_, Phi.T @ y / s2)
    fit = -0.5 * n_ * np.log(2 * np.pi * s2) - ((y - Phi @ m_) ** 2).sum() / (2 * s2)
    occam = -0.5 * alpha * m_ @ m_ + 0.5 * d_ * np.log(alpha) - 0.5 * np.linalg.slogdet(A_)[1]
    return fit + occam, fit, occam


xp = np.linspace(-1, 1, 30); sig_n, alpha_p = 0.3, 0.5
yp = 0.5 + 1.0 * xp - 2.0 * xp ** 2 + rng.normal(0, sig_n, xp.size)
degs = list(range(8)); rows = []
for dg in degs:
    Phi = np.vander(xp, dg + 1, increasing=True)
    le, fit, occam = gauss_evidence(Phi, yp, sig_n ** 2, alpha_p)
    ref_lp = stats.multivariate_normal(np.zeros(xp.size), sig_n ** 2 * np.eye(xp.size) + Phi @ Phi.T / alpha_p).logpdf(yp)
    assert close(le, ref_lp, 1e-7)                                      # closed form = N(y; 0, s2 I + Phi Phi'/alpha)
    bhat = np.linalg.lstsq(Phi, yp, rcond=None)[0]
    ll_hat = stats.norm.logpdf(yp, Phi @ bhat, sig_n).sum()
    bic = -2 * ll_hat + (dg + 1) * np.log(xp.size)
    rows.append((dg, le, ll_hat, -bic / 2, occam, fit))
R = np.array(rows)
assert int(R[:, 1].argmax()) == 2 and int((R[:, 3]).argmax()) == 2      # evidence and BIC both peak at the true degree
assert np.all(np.diff(R[:, 2]) >= -1e-9)                                # the maximized likelihood never decreases
assert np.all(np.diff(R[:, 4]) < 0) and np.allclose(R[:, 1], R[:, 4] + R[:, 5])   # Occam term only worsens; evidence = fit + Occam
print("Ex 4.9: degree | log evidence | max log-lik | -BIC/2 | log Occam | fit at posterior mean")
for r_ in R:
    print(f"  d={int(r_[0])}: " + " ".join(f"{v:8.2f}" for v in r_[1:]))
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(degs, R[:, 2], "o-", label="max log-likelihood"); ax.plot(degs, R[:, 1], "s-", label="log evidence"); ax.plot(degs, R[:, 3], "^--", label="-BIC/2")
ax.set(xlabel="polynomial degree", title="evidence peaks at the true degree (2)"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "p10_evidence.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 11: The Jeffreys-Lindley paradox (notes 3.12; Example 4.10)

# %%
def bf10(xbar, n, sigma, tau):
    """BF10 for H0: mu=0 vs H1: mu~N(0,tau^2), xbar~N(mu, sigma^2/n); ratio of marginal densities of xbar, eq. (3.40)."""
    v0 = sigma ** 2 / n
    return np.exp(stats.norm.logpdf(xbar, 0, np.sqrt(v0 + tau ** 2)) - stats.norm.logpdf(xbar, 0, np.sqrt(v0)))


n_j, xb, sg_j, tau_j = 10 ** 4, 0.02, 1.0, 1.0
z = xb * np.sqrt(n_j) / sg_j
pval = 2 * stats.norm.sf(z)
u = n_j * tau_j ** 2 / sg_j ** 2
assert close(bf10(xb, n_j, sg_j, tau_j), (1 + u) ** -0.5 * np.exp(0.5 * z * z * u / (1 + u)), 1e-10)   # closed form
num_marg = integrate.quad(lambda m_: stats.norm.pdf(xb, m_, sg_j / np.sqrt(n_j)) * stats.norm.pdf(m_, 0, tau_j), -8, 8, points=[0, xb])[0]
assert close(num_marg / stats.norm.pdf(xb, 0, sg_j / np.sqrt(n_j)), bf10(xb, n_j, sg_j, tau_j), 1e-6)  # by quadrature
f_u = lambda uu: (1 + uu) ** -0.5 * np.exp(0.5 * z * z * uu / (1 + uu))
best = optimize.minimize_scalar(lambda lu: -f_u(np.exp(lu)), bounds=(-5, 15), method="bounded")
post_sd = np.sqrt(1 / (n_j / sg_j ** 2 + 1 / tau_j ** 2)); post_m = xb * (n_j / sg_j ** 2) * post_sd ** 2
print(f"Ex 4.10: z={z:.2f}, p={pval:.4f}; BF10={bf10(xb, n_j, sg_j, tau_j):.4f} (BF01={1 / bf10(xb, n_j, sg_j, tau_j):.2f}); "
      f"best-tau BF10={-best.fun:.3f} at tau={np.sqrt(np.exp(best.x) * sg_j ** 2 / n_j):.4f}; e^(z^2/2)={np.exp(z * z / 2):.3f}")
print(f"  under H1 the posterior for mu is N({post_m:.4f}, {post_sd:.4f}^2): 95% interval ({post_m - 1.96 * post_sd:.4f}, {post_m + 1.96 * post_sd:.4f})")
assert abs(pval - 0.0455) < 1e-4 and bf10(xb, n_j, sg_j, tau_j) < 1 / 13 and -best.fun < np.exp(z * z / 2)
assert post_m - 1.96 * post_sd > 0                                      # estimation says mu != 0, testing says H0
ns_ = np.logspace(1, 8, 100); xs_ = 2 / np.sqrt(ns_)
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.loglog(ns_, 1 / bf10(xs_, ns_, 1, 1), label="BF01, tau=1"); ax.loglog(ns_, 1 / bf10(xs_, ns_, 1, 0.1), label="BF01, tau=0.1")
ax.axhline(1, c="k", lw=.5); ax.set(xlabel="n (z fixed at 2, so p = 0.0455 throughout)", ylabel="BF01 (>1 favors H0)", title="Jeffreys-Lindley"); ax.legend()
fig.tight_layout(); fig.savefig(FIG / "p11_jl.png", dpi=110); plt.close(fig)
assert (1 / bf10(xs_, ns_, 1, 1))[-1] > 1e3 and (1 / bf10(xs_, ns_, 1, 1))[0] < 1

# %% [markdown]
# ## Part 12: Posterior predictive checks (notes 3.10; Example 4.11)

# %%
def ppc_pvalue(data, stat, draws=4000):
    """Fit a Normal (reference prior), simulate replicated datasets from p(y_rep|y), return P(T(y_rep) >= T(y) | y), eq. (3.30)."""
    nn_, s_ = data.size, ((data - data.mean()) ** 2).sum()
    s2d = stats.invgamma((nn_ - 1) / 2, scale=s_ / 2).rvs(draws, random_state=rng)
    mud = data.mean() + np.sqrt(s2d / nn_) * rng.standard_normal(draws)
    rep = mud[:, None] + np.sqrt(s2d)[:, None] * rng.standard_normal((draws, nn_))
    T_rep = np.apply_along_axis(stat, 1, rep)
    return np.mean(T_rep >= stat(data)), T_rep


skewed = rng.lognormal(0, 0.8, 60)
normal_data = rng.normal(0, 1, 60)
p_sk, T_sk = ppc_pvalue(skewed, stats.skew)
p_no, _ = ppc_pvalue(normal_data, stats.skew)
print(f"Ex 4.11: skewness check, lognormal data: observed {stats.skew(skewed):.3f}, p_B={p_sk:.4f}; Normal data: p_B={p_no:.3f}")
assert p_sk == 0.0 and 0.05 < p_no < 0.95        # none of 4000 replicates is as skewed as the data
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.hist(T_sk, 60, alpha=.6); ax.axvline(stats.skew(skewed), c="r"); ax.set_title("replicated skewness vs observed (red)")
fig.tight_layout(); fig.savefig(FIG / "p12_ppc.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 13: The likelihood principle and stopping rules (notes 3.13; Example 4.12)

# %%
p_bin = stats.binom.sf(8, 12, .5)           # P(Y >= 9 | n = 12 fixed)
p_nb = stats.nbinom.sf(8, 3, .5)            # P(>= 9 heads before the 3rd tail)
pg_ = np.linspace(1e-4, 1 - 1e-4, 4001)
Lb, Ln = stats.binom.pmf(9, 12, pg_), stats.nbinom.pmf(9, 3, 1 - pg_)
pb_, pn_ = Lb / np.trapezoid(Lb, pg_), Ln / np.trapezoid(Ln, pg_)
assert np.allclose(pb_, pn_) and np.allclose(pb_, stats.beta.pdf(pg_, 10, 4), atol=1e-3)   # flat prior: same Beta(10,4) posterior
print(f"Ex 4.12: 9 heads, 3 tails: p-value {p_bin:.4f} (n fixed) vs {p_nb:.4f} (stop at third tail); both posteriors Beta(10,4)")

# %% [markdown]
# ## Your turn
#
# 1. A/B test: with `prob_B_beats_A` below, add a third arm and report the posterior probability that each arm is best.
#    How does it change if the prior is Beta(20, 180) instead of uniform? (Exercise B7.)
# 2. Gamma-Poisson with an exposure offset: y_i ~ Poisson(lambda t_i) gives Gamma(alpha0 + sum y, beta0 + sum t) (Exercise B5).
#    Fit 3 events in 2.5 hours and 9 events in 4 hours; compare the predictive for the next hour with a Poisson plug-in.
# 3. Sensitivity: plot BF01 for Example 4.10 against tau on a log grid. Where is it largest? (Exercise C2.)
# 4. Extend Part 3 to the Jeffreys prior at theta = 0.01, n = 20: what is the exact coverage? Why?
# 5. Part 10 with sigma unknown: use `nig_reg_logev` from Part 8 and repeat the degree comparison.
#
# ## Reference solution for B7 (write your own first)

# %%
def prob_B_beats_A(yA, nA, yB, nB, prior=(1, 1)):
    """P(theta_B > theta_A | data) for independent Beta posteriors: integral of f_A(t)(1 - F_B(t)) dt."""
    A_ = stats.beta(prior[0] + yA, prior[1] + nA - yA); B_ = stats.beta(prior[0] + yB, prior[1] + nB - yB)
    return integrate.quad(lambda t: A_.pdf(t) * B_.sf(t), 0, 1)[0]


pq = prob_B_beats_A(120, 1000, 145, 1000)
mcA = rng.beta(121, 881, 400_000); mcB = rng.beta(146, 856, 400_000)
assert abs(np.mean(mcB > mcA) - pq) < 0.003
print(f"B7 reference: P(B > A) = {pq:.4f} (uniform priors); with Beta(20,180) priors {prob_B_beats_A(120, 1000, 145, 1000, (20, 180)):.4f}")

# %%
if __name__ == "__main__":
    print("All checks passed. Figures in", FIG)
