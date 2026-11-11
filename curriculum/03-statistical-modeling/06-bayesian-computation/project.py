# %% [markdown]
# # ST-6 Project: Bayesian computation from scratch
#
# Every sampler here is written with numpy/scipy only: random-walk Metropolis, Gibbs (Normal model, linear
# regression, Albert-Chib probit), HMC with leapfrog, R-hat and n_eff, the eight schools (centered vs
# non-centered), hierarchical logistic regression, WAIC and PSIS-LOO, CAVI and reparametrized Gaussian VI.
# Every number in notes.md Section 4 and solutions.md is asserted or printed. Parts: 1 Monte Carlo and the
# curse of dimension; 2 Metropolis; 3 Gibbs; 4 HMC; 5 diagnostics; 6 eight schools; 7 hierarchical logistic;
# 8 variational inference; 9 numbers for the exercises. Run: `MPLBACKEND=Agg python3 project.py`.

# %%
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats, optimize, special
import sympy as sp

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(0)
T0 = time.time()


def save(fig, name):
    fig.tight_layout(); fig.savefig(FIG / name, dpi=110); plt.close(fig)


# %% [markdown]
# ## Shared tools: samplers and diagnostics (notes 3.2, 3.4, 3.5)
# %%
def rwm(logp, x0, scale, S, rng, chol=None):
    """Random-walk Metropolis, proposal x' = x + scale * chol @ z; accept with min(1, pi(x')/pi(x)), eq. (3.5)."""
    d = len(x0); x = np.array(x0, float); lp = logp(x); out = np.empty((S, d)); acc = 0
    C = np.eye(d) if chol is None else chol
    Z = rng.standard_normal((S, d)) @ C.T * scale; U = np.log(rng.random(S))
    for s in range(S):
        y = x + Z[s]; ly = logp(y)
        if U[s] < ly - lp: x, lp = y, ly; acc += 1
        out[s] = x
    return out, acc / S


def leapfrog(th, r, grad, eps, L, g=None):
    """L leapfrog steps for H = -log pi(th) + |r|^2/2 (M = I), eq. (3.14). Returns th, r, grad at th."""
    g = grad(th) if g is None else g
    for _ in range(L):
        r = r + 0.5 * eps * g; th = th + eps * r; g = grad(th); r = r + 0.5 * eps * g
    return th, r, g


def hmc(logp, grad, x0, eps, L, S, rng, jitter=0.2, div=1000.0):
    """HMC (box 3.16): fresh momentum, leapfrog with jittered eps, accept with min(1, exp(-dH)).
    Returns draws, acceptance rate, number of divergences (dH > div), gradient evaluations."""
    x = np.array(x0, float); lp = logp(x); g = grad(x); out = np.empty((S, len(x))); acc = nd = 0
    for s in range(S):
        r0 = rng.standard_normal(len(x)); e = eps * (1 + jitter * (2 * rng.random() - 1))
        y, r, gy = leapfrog(x, r0, grad, e, L, g); ly = logp(y)
        dH = (-ly + r @ r / 2) - (-lp + r0 @ r0 / 2)
        if not np.isfinite(dH) or dH > div: nd += 1
        if np.isfinite(dH) and np.log(rng.random()) < -dH: x, lp, g = y, ly, gy; acc += 1
        out[s] = x
    return out, acc / S, nd, S * L


def acf(x, K=None):
    """Sample autocorrelation rho_k via FFT."""
    x = np.asarray(x, float) - np.mean(x); n = len(x); f = np.fft.rfft(x, 2 * n)
    a = np.fft.irfft(f * np.conj(f))[:n] / n
    return a[: (K or n)] / a[0]


def ess(x):
    """n_eff = S / (1 + 2 sum rho_k), eq. (3.18); sum truncated by Geyer's initial positive sequence."""
    S = len(x)
    if np.var(x) == 0: return 1.0
    r = acf(x); tau = -1.0
    for k in range(0, S - 1, 2):
        p = r[k] + r[k + 1]
        if p < 0: break
        tau += 2 * p
    return S / tau


def rhat(ch):
    """Split R-hat, eq. (3.20): chains (m, n) are halved; sqrt of ((n-1)/n W + B/n) / W."""
    ch = np.asarray(ch); n = ch.shape[1] // 2; c = np.vstack([ch[:, :n], ch[:, n:2 * n]])
    W = c.var(axis=1, ddof=1).mean(); B = n * c.mean(axis=1).var(ddof=1)
    return np.sqrt(((n - 1) / n * W + B / n) / W)


def check_grad(f, g, x, h=1e-6):
    num = np.array([(f(x + h * e) - f(x - h * e)) / (2 * h) for e in np.eye(len(x))])
    assert np.allclose(num, g(x), rtol=1e-4, atol=1e-5)


# %% [markdown]
# ## Part 1: Monte Carlo and the curse of dimension (notes 3.1)
# Target N(0, I_d), proposal N(0, s^2 I_d): rejection acceptance s^{-d}; importance-sampling ESS/N
# = (sqrt(2 - 1/s^2)/s)^d, eq. (3.3). Checked by simulation in low d.
# %%
s_ = 1.5
for d in (1, 10, 100):
    print(f"d={d:3d}: rejection acceptance {s_**-d:.3g}, IS ESS/N {(np.sqrt(2 - 1/s_**2)/s_)**d:.3g}")
for d in (1, 5):
    X = rng.standard_normal((200000, d)) * s_
    lw = stats.norm.logpdf(X).sum(1) - stats.norm.logpdf(X, scale=s_).sum(1); w = np.exp(lw)
    assert abs(w.mean()**2 / (w**2).mean() - (np.sqrt(2 - 1/s_**2)/s_)**d) < 0.02
assert np.isclose(s_**-100, 2.46e-18, rtol=1e-2, atol=0) and np.isclose((np.sqrt(2 - 1/s_**2)/s_)**100, 9.66e-9, rtol=1e-2, atol=0)

# %% [markdown]
# ## Part 2: Metropolis-Hastings (notes 3.2; Example 4.1)
# %%
a_rw = min(1, np.exp(-(1.8**2 - 1.0**2) / 2))                       # N(0,1) target, 1 -> 1.8
a_noh = (3**2 * np.exp(-3)) / (2**2 * np.exp(-2)); a_h = min(1, a_noh * 3 / 2)   # Gamma(3,1), log-normal proposal
print(f"Ex 4.1: RW alpha={a_rw:.4f}; Gamma(3,1) 2->3: ratio {a_noh:.4f}, with Hastings factor 1.5 -> {a_noh*1.5:.4f}")
assert np.isclose(a_rw, 0.3263, atol=1e-4) and np.isclose(a_noh, 0.8277, atol=1e-4) and a_h == 1
# The Hastings factor matters: sample Gamma(3,1) with the multiplicative proposal, with and without it.
def mult_mh(hastings, S=60000):
    x, out, Z, U = 2.0, np.empty(S), rng.normal(0, 0.8, S), np.log(rng.random(S))
    for s in range(S):
        y = x * np.exp(Z[s]); la = 2 * np.log(y / x) - (y - x) + (np.log(y / x) if hastings else 0)
        if U[s] < la: x = y
        out[s] = x
    return out
m_h, m_no = mult_mh(True).mean(), mult_mh(False).mean()
print(f"Gamma(3,1) mean 3: with Hastings {m_h:.3f}, without {m_no:.3f} (wrong target Gamma(2,1), mean 2)")
assert abs(m_h - 3) < 0.15 and abs(m_no - 2) < 0.15

rho = 0.9; Sig = np.array([[1, rho], [rho, 1]]); P2 = np.linalg.inv(Sig)
lp2 = lambda x: -0.5 * x @ P2 @ x; gr2 = lambda x: -P2 @ x
scales = np.geomspace(0.05, 6, 12); res = []
for sc in scales:
    dr, a = rwm(lp2, [0, 0], sc, 20000, rng); res.append((a, np.mean([ess(dr[:, i]) for i in range(2)]) / 20000))
res = np.array(res); best2 = res[np.argmax(res[:, 1]), 0]
d20 = 20; lp20 = lambda x: -0.5 * x @ x; res20 = []
for sc in np.geomspace(0.3, 3, 10) * 2.38 / np.sqrt(d20):
    dr, a = rwm(lp20, np.zeros(d20), sc, 20000, rng); res20.append((a, np.mean([ess(dr[:, i]) for i in range(d20)]) / 20000))
res20 = np.array(res20); best20 = res20[np.argmax(res20[:, 1]), 0]
print(f"RWM d=20 iid: most efficient acceptance {best20:.2f} (theory 0.234 as d -> inf); d=2 rho=0.9 isotropic proposal: {best2:.2f}")
print(np.c_[res20[:, 0], res20[:, 1] * d20].round(3).T)
assert 0.12 < best20 < 0.35
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
for sc, c in ((0.05, "C0"), (2.4, "C1"), (6.0, "C3")):
    ax[0].plot(rwm(lp2, [0, 0], sc, 1000, rng)[0][:, 0], c, lw=.7, label=f"scale {sc}")
ax[0].set(title="RWM traces, rho=0.9", xlabel="iteration"); ax[0].legend(fontsize=7)
ax[1].semilogx(scales, res[:, 0], "o-"); ax[1].set(title="acceptance vs scale", xlabel="proposal scale")
ax[2].plot(res[:, 0], res[:, 1], "o-", label="d=2"); ax[2].plot(res20[:, 0], res20[:, 1] * d20, "s-", label="d=20 (x d)")
ax[2].axvline(0.234, ls=":", c="k"); ax[2].set(xlabel="acceptance rate", ylabel="n_eff / S", title="efficiency"); ax[2].legend()
save(fig, "rwm_scale.png")

# %% [markdown]
# ## Part 3: Gibbs sampling (notes 3.3; Examples 4.2, 4.3)
# Normal model with NIG(mu0=0, kappa0=1, alpha0=2, beta0=2) prior: full conditionals (3.9), (3.10).
# %%
yN = np.array([4.2, 5.1, 3.8, 6.0, 5.4, 4.7]); n = len(yN); ybar = yN.mean(); S_ = ((yN - ybar)**2).sum()
mu0, k0, a0, b0 = 0.0, 1.0, 2.0, 2.0
kn = k0 + n; mun = (k0 * mu0 + n * ybar) / kn; an = a0 + n / 2; bn = b0 + 0.5 * (S_ + k0 * n / kn * (ybar - mu0)**2)
print(f"Ex 4.2: ybar={ybar:.4f}, S={S_:.4f}, kappa_n={kn}, mu_n={mun:.4f}, alpha_n={an}, beta_n={bn:.4f}")
print(f"  mu | s2=1: N({mun:.4f}, {1/kn:.4f});  s2 | mu=mu_n: IG({a0 + (n+1)/2}, {b0 + 0.5*(((yN-mun)**2).sum() + k0*(mun-mu0)**2):.4f})")
assert np.isclose(mun, 4.1714, atol=1e-4) and np.isclose(bn, 13.7671, atol=1e-4)
assert np.isclose(b0 + 0.5 * (((yN - mun)**2).sum() + k0 * (mun - mu0)**2), bn)   # Example 4.2: conditional scale at mu_n equals beta_n
def gibbs_normal(S, rng):
    mu, s2, out = 0.0, 1.0, np.empty((S, 2))
    for s in range(S):
        mu = rng.normal(mun, np.sqrt(s2 / kn))
        s2 = (b0 + 0.5 * (((yN - mu)**2).sum() + k0 * (mu - mu0)**2)) / rng.gamma(a0 + (n + 1) / 2)
        out[s] = mu, s2
    return out
gN = gibbs_normal(40000, rng)[1000::10]
p_mu = stats.kstest(gN[:, 0], stats.t(2 * an, mun, np.sqrt(bn / (an * kn))).cdf).pvalue
p_s2 = stats.kstest(gN[:, 1], stats.invgamma(an, scale=bn).cdf).pvalue
print(f"Gibbs vs exact NIG marginals: KS p-values mu {p_mu:.3f}, sigma^2 {p_s2:.3f}")
assert p_mu > 0.01 and p_s2 > 0.01
# Gibbs on a correlated Gaussian: x | y ~ N(rho y, 1 - rho^2); the x-chain is AR(1) with coefficient rho^2.
def gibbs_biv(rho, S, rng):
    x = y = 0.0; out = np.empty((S, 2)); sd = np.sqrt(1 - rho**2)
    for s in range(S):
        x = rho * y + sd * rng.standard_normal(); y = rho * x + sd * rng.standard_normal(); out[s] = x, y
    return out
gb = gibbs_biv(0.99, 100000, rng); r1 = acf(gb[:, 0], 2)[1]
print(f"Gibbs rho=0.99: lag-1 autocorrelation {r1:.4f} (theory rho^2 = 0.9801), n_eff/S {ess(gb[:,0])/1e5:.4f} (theory 1/99.5)")
assert abs(r1 - 0.9801) < 0.01

# Bayesian linear regression, conjugate NIG prior beta|s2 ~ N(0, s2 V0): Gibbs vs exact (ST-5 3.30).
nL = 50; XL = np.column_stack([np.ones(nL), rng.uniform(-2, 2, nL)]); yL = XL @ [1.0, -0.7] + rng.normal(0, 0.5, nL)
V0i = np.eye(2) / 10; Vn = np.linalg.inv(XL.T @ XL + V0i); mn = Vn @ XL.T @ yL; aL, bL = 2.0, 1.0
def gibbs_lr(S, rng):
    s2, out = 1.0, np.empty((S, 3)); Lc = np.linalg.cholesky(Vn)
    for s in range(S):
        b = mn + np.sqrt(s2) * Lc @ rng.standard_normal(2); res_ = yL - XL @ b
        s2 = (bL + 0.5 * (res_ @ res_ + b @ V0i @ b)) / rng.gamma(aL + (nL + 2) / 2); out[s] = *b, s2
    return out
gl = gibbs_lr(20000, rng)[500:]
s2_exact = (bL + 0.5 * (yL @ yL - mn @ (XL.T @ XL + V0i) @ mn)) / (aL + nL / 2 - 1)
print(f"Gibbs regression: beta {gl[:, :2].mean(0).round(4)} vs exact {mn.round(4)}; E s2 {gl[:, 2].mean():.4f} vs {s2_exact:.4f}")
assert np.allclose(gl[:, :2].mean(0), mn, atol=0.01) and abs(gl[:, 2].mean() / s2_exact - 1) < 0.03

# Albert-Chib probit: z_i | beta, y ~ N(x_i beta, 1) truncated by y_i; beta | z ~ N(B X^T z, B), (3.11)-(3.12).
m_ = 0.5; e1 = m_ + stats.norm.pdf(m_) / stats.norm.cdf(m_); e0 = m_ - stats.norm.pdf(m_) / stats.norm.sf(m_)
print(f"Ex 4.3: E[z | y=1, eta=0.5] = {e1:.4f}, E[z | y=0, eta=0.5] = {e0:.4f}")
assert np.isclose(e1, 1.0092, atol=1e-4) and np.isclose(e0, -0.6411, atol=1e-4)
nP = 400; XP = np.column_stack([np.ones(nP), rng.standard_normal(nP)]); yP = (XP @ [-0.3, 0.8] + rng.standard_normal(nP) > 0).astype(float)
def albert_chib(S, rng):
    B = np.linalg.inv(XP.T @ XP + np.eye(2) / 100); Lc = np.linalg.cholesky(B); b = np.zeros(2); out = np.empty((S, 2))
    for s in range(S):
        m = XP @ b; F0 = stats.norm.cdf(-m); u = rng.random(nP)
        u = np.where(yP == 1, F0 + u * (1 - F0), u * F0)
        z = m + stats.norm.ppf(np.clip(u, 1e-12, 1 - 1e-12))
        b = B @ XP.T @ z + Lc @ rng.standard_normal(2); out[s] = b
    return out
ac = albert_chib(4000, rng)[500:]
nll = lambda b: -(stats.norm.logcdf(XP @ b) * yP + stats.norm.logcdf(-XP @ b) * (1 - yP)).sum()
mle = optimize.minimize(nll, np.zeros(2), method="BFGS").x
print(f"Albert-Chib posterior mean {ac.mean(0).round(3)}, sd {ac.std(0).round(3)}; probit MLE (ST-4) {mle.round(3)}")
assert np.allclose(ac.mean(0), mle, atol=0.04)

# %% [markdown]
# ## Part 4: Hamiltonian Monte Carlo (notes 3.4; Example 4.4)
# %%
th, r, _ = leapfrog(np.array([1.0]), np.array([0.5]), lambda t: -t, 0.1, 1)
dH = (th[0]**2 + r[0]**2) / 2 - 0.625
print(f"Ex 4.4: one leapfrog step -> theta={th[0]:.4f}, rho={r[0]:.5f}, dH={dH:.3e}, accept prob {np.exp(-dH):.5f}")
assert np.isclose(th[0], 1.045) and np.isclose(r[0], 0.39775) and np.isclose(dH, 1.1503e-4, rtol=1e-3)
# Shadow Hamiltonian (3.15): for U = theta^2/2 leapfrog conserves H~ = rho^2/2 + (1 - eps^2/4) theta^2/2 exactly (SymPy),
# so dH = (eps^2/8)(theta'^2 - theta^2): for Ex 4.4 that is 0.00125 * 0.092025.
t_, r_, e_ = sp.symbols("theta rho epsilon"); r1_ = r_ - e_ / 2 * t_; t1_ = t_ + e_ * r1_; r2_ = r1_ - e_ / 2 * t1_
Ht_ = lambda a, b: b**2 / 2 + (1 - e_**2 / 4) * a**2 / 2
assert sp.simplify(Ht_(t1_, r2_) - Ht_(t_, r_)) == 0 and np.isclose(dH, 0.01 / 8 * (1.045**2 - 1))
# Energy error over a trajectory of fixed length T = 2 on the rho=0.9 Gaussian: O(eps^2).
z0, r00 = np.array([1.0, -0.5]), np.array([0.3, 0.8]); Hf = lambda t, p: -lp2(t) + p @ p / 2; errs = []
eps_list = np.array([0.2, 0.1, 0.05, 0.025])
for e in eps_list:
    t, p, g, mx = z0.copy(), r00.copy(), None, 0
    for _ in range(int(round(2 / e))):
        t, p, g = leapfrog(t, p, gr2, e, 1, g); mx = max(mx, abs(Hf(t, p) - Hf(z0, r00)))
    errs.append(mx)
slope = np.polyfit(np.log(eps_list), np.log(errs), 1)[0]
print(f"max |dH| over T=2: {np.array(errs).round(6)}, log-log slope {slope:.3f} (theory 2)")
assert 1.8 < slope < 2.2
acc_eps = [hmc(lp2, gr2, [0, 0], e, max(1, int(1.5 / e)), 1500, rng)[1] for e in (0.05, 0.1, 0.2, 0.3, 0.4, 0.45)]
print("HMC acceptance vs eps (0.05..0.45):", np.round(acc_eps, 3), "(leapfrog unstable above eps = 2 sqrt(1-rho) = 0.632)")
# Volume preservation and reversibility on the banana density log pi = -t1^2/2 - (t2 - t1^2)^2/(2 s^2), s = 0.5.
sb = 0.5
lpb = lambda t: -t[0]**2 / 2 - (t[1] - t[0]**2)**2 / (2 * sb**2)
grb = lambda t: np.array([-t[0] + 2 * t[0] * (t[1] - t[0]**2) / sb**2, -(t[1] - t[0]**2) / sb**2])
check_grad(lpb, grb, np.array([0.3, 0.7]))
def lf_map(z, eps=0.05, L=40):
    t, p, _ = leapfrog(z[:2], z[2:], grb, eps, L); return np.concatenate([t, p])
def euler_map(z, eps=0.05, L=40):
    t, p = z[:2].copy(), z[2:].copy()
    for _ in range(L): t, p = t + eps * p, p + eps * grb(t)
    return np.concatenate([t, p])
zz = np.array([0.5, 0.2, 1.0, -0.4]); h = 1e-6
J = lambda f: np.column_stack([(f(zz + h * e) - f(zz - h * e)) / (2 * h) for e in np.eye(4)])
dl, de = np.linalg.det(J(lf_map)), np.linalg.det(J(euler_map))
back = lf_map(np.concatenate([lf_map(zz)[:2], -lf_map(zz)[2:]]))
print(f"Jacobian det: leapfrog {dl:.8f}, explicit Euler {de:.4f}; reversibility error {np.abs(back[:2]-zz[:2]).max():.1e}")
assert abs(dl - 1) < 1e-6 and abs(de - 1) > 0.05 and np.allclose(back[:2], zz[:2]) and np.allclose(-back[2:], zz[2:])
hb, ab, _, _ = hmc(lpb, grb, [0, 0], 0.08, 25, 6000, rng)
print(f"HMC banana: acc {ab:.3f}, E t1 {hb[:,0].mean():.3f} (0), E t2 {hb[:,1].mean():.3f} (1), Var t1 {hb[:,0].var():.3f} (1), "
      f"Var t2 {hb[:,1].var():.3f} (2.25); n_eff t2 {ess(hb[:,1]):.0f}")
assert abs(hb[:, 1].mean() - 1) < 0.12 and abs(hb[:, 0].var() - 1) < 0.15
g1, g2 = np.meshgrid(np.linspace(-3.5, 3.5, 200), np.linspace(-2, 9, 200))
dens = np.exp(-g1**2 / 2 - (g2 - g1**2)**2 / (2 * sb**2))
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].contour(g1, g2, dens, 8, colors="0.6")
for k in range(4):
    t, p, g = np.array([0.0, 0.0]) + k * np.array([0.4, 0.3]), rng.standard_normal(2), None; path = [t]
    for _ in range(40): t, p, g = leapfrog(t, p, grb, 0.06, 1, g); path.append(t)
    ax[0].plot(*np.array(path).T, ".-", ms=2, lw=1)
ax[0].set(title="leapfrog trajectories on the banana", xlabel="theta1", ylabel="theta2")
ax[1].loglog(eps_list, errs, "o-", label="max |dH|"); ax[1].loglog(eps_list, errs[0] * (eps_list / 0.2)**2, "k:", label="slope 2")
ax[1].set(xlabel="eps", title="energy error, fixed T"); ax[1].legend(); save(fig, "hmc_banana.png")

# Scaling with dimension: n_eff per density/gradient evaluation, iid N(0, I_d).
scal = {}
for d in (2, 10, 50):
    lpd = lambda x: -0.5 * x @ x; grd = lambda x: -x
    dr, a_r = rwm(lpd, np.zeros(d), 2.38 / np.sqrt(d), 30000, rng)
    eps = 1.0 * d**-0.25; L = int(np.ceil(1.2 / eps))
    dh, a_h, _, ng = hmc(lpd, grd, np.zeros(d), eps, L, 3000, rng)
    er = np.mean([ess(dr[:, i]) for i in range(d)]) / 30000; eh = np.mean([ess(dh[:, i]) for i in range(d)]) / ng
    scal[d] = (er, eh, a_r, a_h)
    print(f"d={d:2d}: RWM acc {a_r:.2f}, n_eff/eval {er:.4f} | HMC eps {eps:.3f} L {L} acc {a_h:.2f}, n_eff/grad {eh:.4f} | ratio {eh/er:.1f}")
rat = {d: scal[d][1] / scal[d][0] for d in (2, 10, 50)}
assert rat[2] < rat[10] < rat[50] and rat[50] > 4 * rat[2]                       # HMC's advantage grows with d (notes: 3.8x -> 28x)
assert scal[2][0] > scal[10][0] > scal[50][0] and scal[2][1] > scal[10][1] > scal[50][1]
assert scal[50][0] / scal[2][0] < 0.1 and scal[50][1] / scal[2][1] > 0.25         # RWM falls ~1/d (0.04 of d=2); HMC far more slowly
fig, ax = plt.subplots(figsize=(5, 4)); ds = np.array([2, 10, 50])
ax.loglog(ds, [scal[d][0] for d in ds], "o-", label="RWM (per density eval)")
ax.loglog(ds, [scal[d][1] for d in ds], "s-", label="HMC (per gradient eval)")
ax.loglog(ds, 0.3 / ds, "k:", label="~1/d"); ax.loglog(ds, 0.3 * ds**-0.25, "k--", label="~d^(-1/4)")
ax.set(xlabel="dimension d", ylabel="n_eff per evaluation", title="the scaling plot"); ax.legend(fontsize=8); save(fig, "scaling.png")

# %% [markdown]
# ## Part 5: Diagnostics (notes 3.5; Example 4.5)
# %%
phi, Sa = 0.9, 200000; xa = np.empty(Sa); xa[0] = 0; eps_a = rng.standard_normal(Sa) * np.sqrt(1 - phi**2)
for s in range(1, Sa): xa[s] = phi * xa[s - 1] + eps_a[s]
print(f"Ex 4.5: AR(1) phi=0.9, S={Sa}: n_eff estimate {ess(xa):.0f}, theory S/19 = {Sa/19:.0f}")
assert abs(ess(xa) / (Sa / 19) - 1) < 0.15
mixed = np.array([rwm(lp2, st, 1.2, 4000, rng)[0][1000:, 0] for st in ([-8, -8], [8, 8], [-8, 8], [8, -8])])
lpbi = lambda x: np.logaddexp(-0.5 * (x[0] - 4)**2, -0.5 * (x[0] + 4)**2)
unmixed = np.array([rwm(lpbi, [st], 0.5, 4000, rng)[0][1000:, 0] for st in (-4, 4)])
print(f"R-hat: four dispersed chains on the Gaussian {rhat(mixed):.4f}; two chains stuck in separate modes {rhat(unmixed):.2f}")
assert rhat(mixed) < 1.01 and rhat(unmixed) > 1.1

# %% [markdown]
# ## Part 6: The eight schools (notes 3.6, 3.7; Example 4.6)
# y_j ~ N(theta_j, sigma_j^2), theta_j ~ N(mu, tau^2), mu ~ N(0, 5^2), tau ~ half-Cauchy(0, 5).
# %%
y8 = np.array([28., 8, -3, 7, -1, 1, 18, 12]); s8 = np.array([15., 10, 16, 11, 9, 11, 10, 18]); J8 = 8
def tau_post(tg):
    """Exact p(tau | y) on a grid, eq. (3.25): y ~ N(0, diag(s^2 + tau^2) + 25 11^T) after integrating out mu and theta."""
    lp = np.array([stats.multivariate_normal(np.zeros(J8), np.diag(s8**2 + t**2) + 25).logpdf(y8) for t in tg])
    lp += -np.log1p((tg / 5)**2); w = np.exp(lp - lp.max()); return w / np.trapezoid(w, tg)
def cond_mean(t):
    """E[theta_j | tau, y] = (1 - B_j) y_j + B_j mu_hat(tau), B_j = s_j^2 / (s_j^2 + tau^2), eq. (3.24)."""
    w = 1 / (s8**2 + t**2); muh = (w * y8).sum() / (w.sum() + 1 / 25); B = s8**2 * w; return (1 - B) * y8 + B * muh, B, muh
tg = np.linspace(1e-3, 60, 3000); pt = tau_post(tg); Etau = np.trapezoid(tg * pt, tg)
Eth = np.trapezoid(np.array([cond_mean(t)[0] for t in tg]).T * pt, tg, axis=1)
th10, B10, muh10 = cond_mean(10.0)
print(f"exact: E[tau|y]={Etau:.3f}; E[theta|y] = {Eth.round(2)}")
print(f"shrinkage at tau = E[tau|y]: B_j = {cond_mean(Etau)[1].round(3)}")
print(f"Ex 4.6 (tau=10): mu_hat={muh10:.3f}, B_j = {B10.round(3)}, E[theta_j|tau,y] = {th10.round(2)}")
assert np.isclose(B10[0], 225 / 325) and np.isclose(B10[4], 81 / 181)
assert np.isclose(muh10, 3.663, atol=1e-3) and np.allclose(th10[[0, 4, 2, 6]], [11.15, 1.09, 1.79, 10.83], atol=6e-3)   # Example 4.6
assert np.isclose(Etau, 3.60, atol=6e-3) and np.allclose(Eth, [6.21, 4.94, 3.93, 4.76, 3.62, 4.04, 6.30, 4.85], atol=6e-3)   # notes 3.6
ybw = (y8 / s8**2).sum() / (1 / s8**2).sum(); vbw = 1 / (1 / s8**2).sum()
print(f"complete pooling: theta ~ N({ybw:.3f}, {np.sqrt(vbw):.3f}^2)")

def lp_c(q):
    mu, lt, th = q[0], q[1], q[2:]; t = np.exp(lt)
    return -mu**2 / 50 - np.log1p(t**2 / 25) + lt - J8 * lt - ((th - mu)**2).sum() / (2 * t**2) - ((y8 - th)**2 / (2 * s8**2)).sum()
def gr_c(q):
    mu, lt, th = q[0], q[1], q[2:]; t = np.exp(lt)
    return np.concatenate([[-mu / 25 + (th - mu).sum() / t**2, -2 * t**2 / (25 + t**2) + 1 - J8 + ((th - mu)**2).sum() / t**2],
                           -(th - mu) / t**2 + (y8 - th) / s8**2])
def lp_nc(q):
    mu, lt, et = q[0], q[1], q[2:]; t = np.exp(lt)
    return -mu**2 / 50 - np.log1p(t**2 / 25) + lt - (et**2).sum() / 2 - ((y8 - mu - t * et)**2 / (2 * s8**2)).sum()
def gr_nc(q):
    mu, lt, et = q[0], q[1], q[2:]; t = np.exp(lt); rr = (y8 - mu - t * et) / s8**2
    return np.concatenate([[-mu / 25 + rr.sum(), -2 * t**2 / (25 + t**2) + 1 + t * (rr * et).sum()], -et + t * rr])
qq = rng.standard_normal(10); check_grad(lp_c, gr_c, qq); check_grad(lp_nc, gr_nc, qq)
runs = {}; np.seterr(over="ignore", invalid="ignore", divide="ignore")
for name, f, g in (("centered", lp_c, gr_c), ("non-centered", lp_nc, gr_nc)):
    ch = [hmc(f, g, rng.normal(0, 1, 10), 0.15, 25, 2500, rng) for _ in range(4)]
    dr = np.array([c[0][500:] for c in ch]); runs[name] = dr
    print(f"  per-chain divergences {[c[2] for c in ch]}, per-chain mean log tau {dr[:, :, 1].mean(1).round(2)}")
    ndv = sum(c[2] for c in ch); ne = ess(dr[:, :, 1].ravel()); rh = rhat(dr[:, :, 1])
    print(f"{name:12s}: acc {np.mean([c[1] for c in ch]):.2f}, divergences {ndv}, n_eff(log tau) {ne:.0f}, R-hat(log tau) {rh:.3f}, "
          f"E tau {np.exp(dr[:, :, 1]).mean():.2f}, min tau {np.exp(dr[:, :, 1]).min():.3f}")
    runs[name + "_stats"] = (ndv, ne)
assert runs["non-centered_stats"][1] > 2 * runs["centered_stats"][1] and runs["centered_stats"][0] > runs["non-centered_stats"][0]
nc = runs["non-centered"].reshape(-1, 10); tau_nc = np.exp(nc[:, 1]); th_nc = nc[:, :1] + tau_nc[:, None] * nc[:, 2:]
assert abs(tau_nc.mean() - Etau) < 0.5 and np.allclose(th_nc.mean(0), Eth, atol=1.0)
cc = runs["centered"].reshape(-1, 10)
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
ax[0].scatter(cc[:, 2], cc[:, 1], s=2, alpha=.3, label="centered"); ax[0].scatter(th_nc[:, 0], nc[:, 1], s=2, alpha=.2, label="non-centered")
ax[0].set(xlabel="theta_1", ylabel="log tau", title="the funnel"); ax[0].legend(markerscale=5)
ax[1].scatter(nc[:, 2], nc[:, 1], s=2, alpha=.3); ax[1].set(xlabel="eta_1", ylabel="log tau", title="non-centered coordinates")
for j in range(J8): ax[2].plot([0, 1, 2], [y8[j], th_nc[:, j].mean(), ybw], "o-", lw=1)
ax[2].set_xticks([0, 1, 2], ["no pooling", "hierarchical", "complete"]); ax[2].set(title="shrinkage"); save(fig, "eight_schools.png")
# Posterior predictive check: T = max y_j - min y_j, replicated with the same schools.
yrep = th_nc + rng.standard_normal(th_nc.shape) * s8; Trep = np.ptp(yrep, axis=1)
ppp = (Trep >= np.ptp(y8)).mean(); print(f"PPC: T(y) = range = {np.ptp(y8):.0f}, P(T(yrep) >= T(y)) = {ppp:.2f}")
assert 0.05 < ppp < 0.95

def psis_lw(lw):
    """Pareto-smoothed log importance weights (simple form): fit a generalized Pareto to the M largest ratios."""
    S = len(lw); M = int(min(S / 5, 3 * np.sqrt(S))); o = np.argsort(lw); w = np.exp(lw - lw.max())
    u = w[o[-M - 1]]; exc = w[o[-M:]] - u
    k, _, sc = stats.genpareto.fit(exc, floc=0)
    w[o[-M:]] = np.minimum(u + stats.genpareto.ppf((np.arange(1, M + 1) - 0.5) / M, k, scale=sc), w.max())
    return np.log(w), k
def waic_loo(ll):
    """ll: (S, n) pointwise log-likelihoods. Returns elpd_waic, p_waic, elpd_loo, max k-hat, eqs. (3.29)-(3.31)."""
    lppd = special.logsumexp(ll, axis=0) - np.log(len(ll)); pw = ll.var(axis=0, ddof=1)
    loo, ks = [], []
    for i in range(ll.shape[1]):
        lw, k = psis_lw(-ll[:, i]); loo.append(special.logsumexp(lw + ll[:, i]) - special.logsumexp(lw)); ks.append(k)
    return (lppd - pw).sum(), pw.sum(), np.sum(loo), max(ks)
Sd = len(th_nc)
draws = {"no pooling": y8 + rng.standard_normal((Sd, J8)) * s8,
         "complete": np.repeat(ybw + np.sqrt(vbw) * rng.standard_normal((Sd, 1)), J8, axis=1), "hierarchical": th_nc}
ic = {}
for k_, th_ in draws.items():
    ic[k_] = waic_loo(stats.norm.logpdf(y8, th_, s8))
    print(f"{k_:12s}: elpd_waic {ic[k_][0]:7.2f} (p_waic {ic[k_][1]:.2f}), elpd_loo {ic[k_][2]:7.2f}, max k-hat {ic[k_][3]:.2f}")
assert ic["no pooling"][2] < ic["hierarchical"][2] and ic["no pooling"][3] > 0.7

# %% [markdown]
# ## Part 7: Hierarchical logistic regression, random intercepts (notes 3.6)
# %%
Jg = 8; ng = rng.integers(5, 60, Jg); grp = np.repeat(np.arange(Jg), ng); xg = rng.standard_normal(len(grp))
alpha_true = -0.5 + 0.8 * rng.standard_normal(Jg); yg = (rng.random(len(grp)) < special.expit(alpha_true[grp] + 1.2 * xg)).astype(float)
def lp_h(q):
    mu, lt, b, et = q[0], q[1], q[2], q[3:]; t = np.exp(lt); eta = mu + t * et[grp] + b * xg
    return (yg * eta - np.logaddexp(0, eta)).sum() - mu**2 / 8 - b**2 / 8 - t**2 / 2 + lt - (et**2).sum() / 2
def gr_h(q):
    mu, lt, b, et = q[0], q[1], q[2], q[3:]; t = np.exp(lt); eta = mu + t * et[grp] + b * xg
    rr = yg - special.expit(eta); da = np.bincount(grp, rr, Jg)
    return np.concatenate([[da.sum() - mu / 4, t * (et * da).sum() - t**2 + 1, (rr * xg).sum() - b / 4], t * da - et])
check_grad(lp_h, gr_h, rng.standard_normal(Jg + 3) * 0.5)
hl = np.array([hmc(lp_h, gr_h, np.zeros(Jg + 3), 0.1, 20, 1500, rng)[0][300:] for _ in range(2)])
hf = hl.reshape(-1, Jg + 3); a_post = hf[:, :1] + np.exp(hf[:, 1:2]) * hf[:, 3:]
raw = np.array([special.logit(np.clip(yg[grp == j].mean(), 0.02, 0.98)) for j in range(Jg)])
print(f"hier. logistic: beta {hf[:,2].mean():.2f} (true 1.2), mu {hf[:,0].mean():.2f}, tau {np.exp(hf[:,1]).mean():.2f}; "
      f"R-hat max {max(rhat(hl[:,:,i]) for i in range(Jg+3)):.3f}")
print("  n_j:", ng, "\n  raw logit:", raw.round(2), "\n  posterior alpha_j:", a_post.mean(0).round(2), "\n  true alpha_j:", alpha_true.round(2))
lo, hi = np.percentile(hf[:, 2], [2.5, 97.5]); assert lo < 1.2 < hi

# %% [markdown]
# ## Part 8: Variational inference (notes 3.8; Example 4.7)
# CAVI for the Normal model with the NIG prior of Part 3: q(mu) = N(mu_n, 1/(kappa_n E[1/s2])), q(s2) = IG(A, B).
# %%
A = a0 + (n + 1) / 2; Eprec = 1.0; hist = []
for it in range(50):
    v = 1 / (kn * Eprec)
    B = b0 + 0.5 * (((yN - mun)**2).sum() + n * v + k0 * ((mun - mu0)**2 + v)); Eprec = A / B; hist.append((v, B))
assert np.allclose(hist[0], (1 / 7, 14.2671), atol=1e-4) and np.allclose(hist[1], (0.37058, 15.0642), atol=1e-4) and np.allclose(hist[2], (0.39128, 15.1366), atol=1e-4)
print("Ex 4.7: CAVI (Var_q mu, B) by iteration:", ", ".join(f"({a:.5f}, {b:.4f})" for a, b in hist[:3]), f"... -> ({v:.5f}, {B:.4f})")
var_exact = bn / (kn * (an - 1))
print(f"  fixed point: A={A}, B={B:.4f} (= beta_n (2a_n+1)/(2a_n) = {bn*(2*an+1)/(2*an):.4f}), E_q[1/s2] {Eprec:.4f} = a_n/b_n {an/bn:.4f}")
print(f"  Var_q(mu) {v:.5f} vs exact {var_exact:.5f}; ratio {v/var_exact:.4f} = (a_n-1)/a_n = {(an-1)/an:.4f}")
assert np.isclose(B, bn * (2 * an + 1) / (2 * an)) and np.isclose(v / var_exact, (an - 1) / an)
# Gaussian VI on the banana: q = N(m, L L^T), reparametrization theta = m + L z, ELBO = E[log pi] + log|det L| + const.
m, Lv = np.zeros(2), np.eye(2) * 0.5; lr = 0.02
gb_ = lambda T: np.column_stack([-T[:, 0] + 2 * T[:, 0] * (T[:, 1] - T[:, 0]**2) / sb**2, -(T[:, 1] - T[:, 0]**2) / sb**2])
for it in range(4000):
    Z = rng.standard_normal((100, 2)); T = m + Z @ Lv.T; G = gb_(T)
    gm = G.mean(0); gL = np.tril(G.T @ Z / len(Z) + np.diag(1 / np.diag(Lv)))
    step = lr / (1 + it / 1000); m = m + step * gm; Lv = Lv + step * gL
Cq = Lv @ Lv.T
print(f"Gaussian VI on banana: mean {m.round(3)}, cov diag {np.diag(Cq).round(3)} vs HMC mean {hb.mean(0).round(3)}, var {hb.var(0).round(3)}")
assert Cq[1, 1] < 0.5 * hb[:, 1].var() and Cq[0, 0] < hb[:, 0].var()
fig, ax = plt.subplots(figsize=(5, 4.5)); ax.contour(g1, g2, dens, 8, colors="0.6")
ax.scatter(hb[::3, 0], hb[::3, 1], s=2, alpha=.3, label="HMC")
qd = rng.multivariate_normal(m, Cq, 2000); ax.scatter(qd[:, 0], qd[:, 1], s=2, alpha=.3, c="C3", label="Gaussian VI")
ax.set(xlim=(-3.5, 3.5), ylim=(-2, 9), title="VI is mode-seeking"); ax.legend(markerscale=5); save(fig, "vi_banana.png")

# %% [markdown]
# ## Part 9: Numbers for the exercises
# %%
print(f"A1: rejection acceptance 1.2^-d for d=10, 50: {1.2**-10:.4f}, {1.2**-50:.3g}")
print(f"A2: N(0,1) target, RW from 0.5 to -1.2: alpha = {np.exp(-(1.44 - 0.25)/2):.4f}")
tA, rA, _ = leapfrog(np.array([0.]), np.array([1.]), lambda t: -t, 0.2, 2); dHA = (tA[0]**2 + rA[0]**2 - 1) / 2
print(f"A3: two leapfrog steps, eps=0.2 from (0,1): theta={tA[0]:.5f}, rho={rA[0]:.5f}, dH={dHA:.4e} (eps^2/8 theta'^2 = {0.04/8*tA[0]**2:.4e})")
assert np.isclose(dHA, 0.04 / 8 * tA[0]**2)
for ph in (0.5, 0.99, -0.5): print(f"A4: AR(1) phi={ph}: n_eff/S = {(1-ph)/(1+ph):.4f}")
print(f"B2: Gibbs rho=0.9 lag-1 = {0.81:.2f}, n_eff/S = {(1-0.81)/(1+0.81):.4f}; rho=0.99: {(1-0.9801)/(1.9801):.5f}")
Bs = 4**2 / (4**2 + np.array([1, 4, 16])**2); print("B5: shrinkage sigma=4, tau=1,4,16:", Bs.round(4), "; E[theta|tau] with mu_hat=10, y=30:", ((1 - Bs) * 30 + Bs * 10).round(2))
ch_b = np.array([[1.0, 1.2, 0.8, 1.1], [3.0, 3.2, 2.8, 3.1]]); W = ch_b.var(1, ddof=1).mean(); Bb = 4 * ch_b.mean(1).var(ddof=1)
print(f"B6: W={W:.4f}, B={Bb:.4f}, var+={3/4*W + Bb/4:.4f}, R-hat (no split) = {np.sqrt((3/4*W + Bb/4)/W):.3f}")
print(f"B7: AR(1) phi=0.9: n_eff/S unthinned {1/19:.5f}; thinned by 10: {(1-0.9**10)/(1+0.9**10)/10:.5f}")
print(f"B9(d): Var_q(s2) = {B**2/((A-1)**2*(A-2)):.4f} vs exact {bn**2/((an-1)**2*(an-2)):.4f}; E_q s2 {B/(A-1):.4f} vs {bn/(an-1):.4f}")
print(f"B9: CAVI variance ratio (a_n-1)/a_n for a_n = 2,5,50: {[(a-1)/a for a in (2,5,50)]}")
# B8 reference: the centered model at other step sizes (same trajectory length eps*L = 3.75).
for e in (0.05, 0.3):
    ch = [hmc(lp_c, gr_c, rng.normal(0, 1, 10), e, int(round(3.75 / e)), 1000, rng) for _ in range(4)]
    dr = np.array([c[0][200:, 1] for c in ch])
    print(f"B8: centered eps={e}: divergences {[c[2] for c in ch]}, n_eff(log tau) {ess(dr.ravel()):.0f}, "
          f"R-hat {rhat(dr):.2f}, min tau {np.exp(dr.min()):.3f}, E tau {np.exp(dr).mean():.2f} (exact {Etau:.2f})")
# C2 reference: parallel tempering on 0.5 N(-4,1) + 0.5 N(4,1), temperatures 1, 2, 4, 8.
def ptemper(x0, S, rng, Ts=(1, 2, 4, 8)):
    x = np.full(len(Ts), x0, float); lpx = np.array([lpbi([v]) for v in x]); out = np.empty(S); sw = 0
    for s in range(S):
        for i, T in enumerate(Ts):
            y = x[i] + rng.normal(0, 1.5 * np.sqrt(T)); ly = lpbi([y])
            if np.log(rng.random()) < (ly - lpx[i]) / T: x[i], lpx[i] = y, ly
        i = rng.integers(len(Ts) - 1); la = (lpx[i + 1] - lpx[i]) * (1 / Ts[i] - 1 / Ts[i + 1])
        if np.log(rng.random()) < la: x[[i, i + 1]] = x[[i + 1, i]]; lpx[[i, i + 1]] = lpx[[i + 1, i]]; sw += 1
        out[s] = x[0]
    return out, sw / S
pt = [ptemper(st, 8000, rng) for st in (-4, 4)]; ptc = np.array([c[0][1000:] for c in pt])
print(f"C2: parallel tempering R-hat {rhat(ptc):.3f} (vs {rhat(unmixed):.2f} without), P(x>0) {(ptc > 0).mean():.3f}, "
      f"swap rate {np.mean([c[1] for c in pt]):.2f}")
assert rhat(ptc) < 1.02
print(f"\nAll checks passed in {time.time() - T0:.1f} s.")

# %% [markdown]
# ## Your turn
# 1. **Ising denoising (SM-6).** Corrupt a 40x40 binary image by flipping 10% of pixels. Posterior over the clean image
#    x in {-1,1}^N: log pi(x) = J sum_<ij> x_i x_j + h sum_i x_i y_i. Write single-site Metropolis (or Gibbs: the full
#    conditional of x_i is logistic in its neighbours' sum). Report the posterior-mean image and the error rate vs J.
# 2. **NUTS in its simplest form.** Replace fixed L by "integrate until (theta - theta0) . rho < 0" (the U-turn) and
#    accept the endpoint with the usual correction. Show on the banana that it adapts trajectory length, then read
#    Hoffman & Gelman (2014) for why the naive version is not reversible and how the tree-doubling fixes it.
# 3. **Parallel tempering.** Target 0.5 N(-4,1) + 0.5 N(4,1) (Part 5's unmixed example). Run chains at
#    pi^(1/T) for T in {1, 2, 4, 8} with swap moves accepted with min(1, (pi_i(x_j) pi_j(x_i))/(pi_i(x_i) pi_j(x_j))).
#    Show that R-hat of the T=1 chain drops below 1.01. (Exercise C2.)
# 4. Change Part 6's step size to 0.1 and 0.4 for the centered model; watch the divergence count and n_eff(log tau).
#    Does a smaller step "fix" the funnel, or only move where the sampler gets stuck? (Exercise B8.)
