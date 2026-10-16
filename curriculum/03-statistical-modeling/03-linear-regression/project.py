# %% [markdown]
# # ST-3 Project: Least Squares, Inference, Diagnostics, Ridge and Lasso
#
# Everything in the notes is implemented from scratch (no statsmodels, no sklearn) and checked
# against closed forms, brute force, or Monte Carlo. Parts: (1) OLS three ways and the failure of
# the normal equations, (2) the inference table, Example 4.1 and FWL (Example 4.2), (3) Gauss-Markov,
# coverage of intervals over 10^4 replications, (4) leverage, Cook's distance, the LOO identity,
# (5) heteroscedasticity: sandwich and WLS, (6) multicollinearity: VIF and the SVD, (7) model
# selection: C_p, AIC, BIC, LOO, post-selection, (8) ridge, (9) lasso by coordinate descent,
# (10) numbers used in exercises and solutions, (11) "Your turn" with two reference solutions.
# Every number quoted in notes.md, exercises.md or solutions.md is asserted or printed here.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures go to ./figures/.

# %%
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg, optimize, stats

T0 = time.time()
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
np.set_printoptions(precision=4, suppress=True)
close = lambda a, b, tol=5e-4: bool(np.all(np.abs(np.asarray(a) - np.asarray(b)) < tol))


def hat_diag(X, lam=0.0):
    """h_ii of H = X (X'X + lam I)^{-1} X' (notes (3.5); lam > 0 is ridge, (3.40))."""
    return np.einsum("ij,ji->i", X, np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T))


# %% [markdown]
# ## Part 1: OLS three ways (notes 3.2; LA-5 (3.6), (3.16); LA-9 (3.11); LA-13 Theorem 3.12)

# %%
def ols_normal(X, y):
    """Normal equations X'X b = X'y, (3.4). Squares the condition number."""
    return np.linalg.solve(X.T @ X, X.T @ y)


def ols_qr(X, y):
    """X = QR, then R b = Q'y by back substitution (LA-5 (3.16))."""
    Q, R = np.linalg.qr(X)
    return linalg.solve_triangular(R, Q.T @ y)


def ols_svd(X, y):
    """b = X^+ y = V diag(1/s) U'y (LA-9 (3.11))."""
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    return Vt.T @ ((U.T @ y) / s)


def make_design(n, s, rng):
    """X = U diag(s) V' with random orthonormal U (n x p), V (p x p): singular values exactly s."""
    U = np.linalg.qr(rng.standard_normal((n, len(s))))[0]
    V = np.linalg.qr(rng.standard_normal((len(s), len(s))))[0]
    return U @ np.diag(s) @ V.T


rng = np.random.default_rng(1)
Xw = rng.standard_normal((60, 5))
yw = Xw @ np.arange(1.0, 6) + 0.1 * rng.standard_normal(60)
ref = np.linalg.lstsq(Xw, yw, rcond=None)[0]
for f in (ols_normal, ols_qr, ols_svd):
    assert np.allclose(f(Xw, yw), ref, atol=1e-10)
Xi = make_design(60, np.logspace(0, -8, 5), rng)  # kappa(X) = 1e8, kappa(X'X) = 1e16
yi = Xi @ np.ones(5)  # zero residual: the error below is pure conditioning
err = {f.__name__: np.linalg.norm(f(Xi, yi) - np.ones(5)) for f in (ols_normal, ols_qr, ols_svd)}
print("Part 1  kappa(X) = %.1e;  error of beta-hat:" % np.linalg.cond(Xi), {k: "%.1e" % v for k, v in err.items()})
assert err["ols_qr"] < 1e-6 and err["ols_svd"] < 1e-6 and err["ols_normal"] > 1e3 * max(err["ols_qr"], err["ols_svd"])

# %% [markdown]
# ## Part 2: the inference table (notes 3.5), Example 4.1 and Example 4.2 (FWL)

# %%
def ols_fit(X, y, alpha=0.05):
    """OLS with the full table: b, se, t, p, CI, RSS, TSS, R2, adjR2, s^2, F (overall), leverage h, residuals e.
    Assumes X has an intercept column. se_j = s sqrt([(X'X)^-1]_jj) from (3.20); F is (3.25) with q = p-1."""
    n, p = X.shape
    A = np.linalg.inv(X.T @ X)
    b = A @ X.T @ y
    e = y - X @ b
    rss, tss, df = e @ e, ((y - y.mean()) ** 2).sum(), n - p
    s2 = rss / df
    se = np.sqrt(s2 * np.diag(A))
    tc = stats.t.ppf(1 - alpha / 2, df)
    F = ((tss - rss) / (p - 1)) / s2
    return dict(b=b, se=se, t=b / se, p=2 * stats.t.sf(abs(b / se), df), lo=b - tc * se, hi=b + tc * se, rss=rss, tss=tss,
                R2=1 - rss / tss, adjR2=1 - (rss / df) / (tss / (n - 1)), s2=s2, F=F, pF=stats.f.sf(F, p - 1, df),
                h=np.einsum("ij,jk,ik->i", X, A, X), e=e, A=A, df=df, tc=tc)


x41, y41 = np.arange(1.0, 6), np.array([2.0, 4, 5, 4, 5])
X41 = np.c_[np.ones(5), x41]
f41 = ols_fit(X41, y41)
Sxx, Sxy = ((x41 - 3) ** 2).sum(), ((x41 - 3) * (y41 - 4)).sum()
assert (Sxx, Sxy) == (10.0, 6.0) and close(f41["b"], [2.2, 0.6])
assert close(f41["e"], [-0.8, 0.6, 1.0, -0.6, -0.2]) and close([f41["rss"], f41["tss"]], [2.4, 6.0])
assert close([f41["R2"], f41["adjR2"], f41["s2"]], [0.6, 0.46667, 0.8]) and close(f41["se"], [0.93808, 0.28284], 5e-5)
assert close(f41["t"][1], 2.1213, 5e-5) and close(f41["tc"], 3.1824, 5e-5) and close(f41["p"][1], 0.1240, 5e-4)
assert close([f41["lo"][1], f41["hi"][1]], [-0.3001, 1.5001], 5e-4) and close(f41["F"], f41["t"][1] ** 2, 1e-9) and close(f41["F"], 4.5)
assert close(f41["h"], [0.6, 0.3, 0.2, 0.3, 0.6]) and abs(f41["h"].sum() - 2) < 1e-12
for x0, ch, ph in [(3.0, 1.2730, 3.1181), (5.0, 2.2049, 3.6005)]:  # CI for the mean vs PI for a new y, (3.22)-(3.23)
    v = np.array([1.0, x0]); h0 = v @ f41["A"] @ v; w = f41["tc"] * np.sqrt(f41["s2"])
    assert close([w * np.sqrt(h0), w * np.sqrt(1 + h0)], [ch, ph], 5e-4)
    print("Example 4.1: x0=%g  yhat=%.1f  h0=%.1f  CI half-width %.4f  PI half-width %.4f" % (x0, v @ f41["b"], h0, w * np.sqrt(h0), w * np.sqrt(1 + h0)))
print("Example 4.1: b =", f41["b"], " RSS, TSS =", f41["rss"], f41["tss"], " t =", f41["t"][1].round(4), " p =", f41["p"][1].round(4), " CI =", f41["lo"][1].round(3), f41["hi"][1].round(3))

# Example 4.2: Frisch-Waugh-Lovell, three predictors plus an intercept, both routes.
x1 = np.arange(1.0, 9); x2 = np.array([2, 1, 4, 3, 6, 5, 8, 6.0]); x3 = np.array([1, 0, 2, 1, 3, 1, 2, 4.0])
y42 = np.array([3.1, 2.2, 5.9, 5.0, 8.1, 6.8, 9.9, 9.3])
X42 = np.c_[np.ones(8), x1, x2, x3]
b42 = ols_qr(X42, y42)
X_other = X42[:, :3]
M1 = np.eye(8) - X_other @ np.linalg.solve(X_other.T @ X_other, X_other.T)  # I - P_1 (annihilator of the other columns)
xt, yt = M1 @ x3, M1 @ y42
b_fwl = (xt @ yt) / (xt @ xt)
res_full = y42 - X42 @ b42
assert close(b42, [0.9948, 0.2529, 0.7944, 0.3883], 5e-5) and close(b_fwl, b42[3], 1e-12)
assert close(xt, [0.3727, -0.5206, 0.5674, -0.3258, 0.7622, -1.1311, -1.0431, 1.3184], 5e-5)
assert close([xt @ yt, xt @ xt], [2.1449, 5.5243], 5e-5) and close(yt - b_fwl * xt, res_full, 1e-12)  # residuals coincide
assert close(res_full @ res_full, 0.12579, 5e-5)
print("Example 4.2: full-model b =", b42, " FWL slope =", round(b_fwl, 5))
# Var(b_j) = sigma^2 / ||x~_j||^2 (3.18): the FWL route also gives the standard error (with df n-p, not n-1-...)
assert close(np.linalg.inv(X42.T @ X42)[3, 3], 1 / (xt @ xt), 1e-12)

# %% [markdown]
# ## Part 3: Gauss-Markov, E and Var of b, unbiased s^2, coverage over 10^4 replications (notes 3.4-3.5)

# %%
rng = np.random.default_rng(3)
# Example 4.3: sample-mean weights, n = 4.
c = np.array([0.4, 0.3, 0.2, 0.1]); d = c - 0.25
assert close([c @ c, 0.25, d @ d, 0.25 + d @ d], [0.30, 0.25, 0.05, 0.30], 1e-12) and abs(d.sum()) < 1e-15

n, p, sig = 30, 4, 2.0
X = np.c_[np.ones(n), rng.uniform(-1, 1, (n, 3))]
beta = np.array([1.0, 2.0, -1.0, 0.5])
x0 = np.array([1.0, 0.3, -0.2, 0.5])
R = 10_000
A = np.linalg.solve(X.T @ X, X.T)
H = X @ A
Y = (X @ beta)[:, None] + sig * rng.standard_normal((n, R))
B = A @ Y
RSS = ((Y - H @ Y) ** 2).sum(0)
s2 = RSS / (n - p)
V = np.linalg.inv(X.T @ X)
assert np.allclose(np.cov(B), sig**2 * V, rtol=0.1, atol=0.01) and abs(B.mean(1) - beta).max() < 0.05  # E b = beta, Var b = sigma^2 (X'X)^-1
assert abs(s2.mean() / sig**2 - 1) < 0.02 and abs((RSS / n).mean() / (sig**2 * (n - p) / n) - 1) < 0.02  # (3.12); MLE biased by (n-p)/n
tq = stats.t.ppf(0.975, n - p)
half = tq * np.sqrt(s2[None, :] * np.diag(V)[:, None])
cov_beta = (np.abs(B - beta[:, None]) <= half).mean(1)
h0 = x0 @ V @ x0
ynew = x0 @ beta + sig * rng.standard_normal(R)
cov_pred = (np.abs(ynew - x0 @ B) <= tq * np.sqrt(s2 * (1 + h0))).mean()
cov_mean = (np.abs(x0 @ beta - x0 @ B) <= tq * np.sqrt(s2 * h0)).mean()
print("Part 3  coverage of beta_j intervals:", cov_beta, " mean-response CI: %.4f  prediction interval: %.4f" % (cov_mean, cov_pred))
assert np.all(np.abs(cov_beta - 0.95) < 0.01) and abs(cov_pred - 0.95) < 0.01 and abs(cov_mean - 0.95) < 0.01
assert (np.abs(ynew - x0 @ B) <= tq * np.sqrt(s2 * h0)).mean() < 0.5  # using the CI half-width as if it were a PI covers far too little
# Under H0 (all slopes zero) R^2 ~ Beta((p-1)/2, (n-p)/2): mean (p-1)/(n-1), the reason for adjusted R^2 (3.9).
Y0 = sig * rng.standard_normal((n, R))
Yc = Y0 - Y0.mean(0)
R2 = 1 - ((Y0 - H @ Y0) ** 2).sum(0) / (Yc**2).sum(0)
assert abs(R2.mean() - (p - 1) / (n - 1)) < 0.003
adj = 1 - (1 - R2) * (n - 1) / (n - p)
print("Part 3  null R^2 mean %.4f (theory %.4f), null adjusted R^2 mean %.4f" % (R2.mean(), (p - 1) / (n - 1), adj.mean()))
assert abs(adj.mean()) < 0.005
# F test = LRT: n log(1 + qF/(n-p)) = n log(RSS0/RSS); and the nested F with q = 2 rejects 5% under H0.
X0 = X[:, :2]
RSS0 = ((Y - X0 @ np.linalg.solve(X0.T @ X0, X0.T @ Y)) ** 2).sum(0)  # here H0 (beta_3 = beta_4 = 0) is FALSE for 2 of the 4
Fs = ((RSS0 - RSS) / 2) / s2
assert np.allclose(n * np.log(1 + 2 * Fs / (n - p)), n * np.log(RSS0 / RSS))
Yn = (X0 @ beta[:2])[:, None] + sig * rng.standard_normal((n, R))  # data from the reduced model
r1 = ((Yn - H @ Yn) ** 2).sum(0); r0 = ((Yn - X0 @ np.linalg.solve(X0.T @ X0, X0.T @ Yn)) ** 2).sum(0)
rej = (((r0 - r1) / 2) / (r1 / (n - p)) > stats.f.ppf(0.95, 2, n - p)).mean()
assert abs(rej - 0.05) < 0.01
print("Part 3  nested F test size under H0: %.4f" % rej)

# %% [markdown]
# ## Part 4: leverage, studentized residuals, Cook's distance, the LOO identity (notes 3.6, Example 4.4)

# %%
def influence(X, y):
    """h_ii, r_i = e_i/(s sqrt(1-h_ii)), t_i (external), D_i = r_i^2 h_ii /(p (1-h_ii)), via closed forms (3.29)-(3.33)."""
    n, p = X.shape
    f = ols_fit(X, y)
    e, h = f["e"], f["h"]
    r = e / np.sqrt(f["s2"] * (1 - h))
    s2_i = (f["rss"] - e**2 / (1 - h)) / (n - p - 1)  # leave-one-out variance, from RSS_(i) = RSS - e_i^2/(1-h_ii)
    return h, r, e / np.sqrt(s2_i * (1 - h)), r**2 * h / (p * (1 - h)), f


def brute_loo(X, y):
    """Refit n times: returns (cook D by definition, LOO prediction errors, RSS_(i), s2_(i))."""
    n, p = X.shape
    b = ols_qr(X, y); s2 = ((y - X @ b) ** 2).sum() / (n - p)
    D, loo, rss_i = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        m = np.arange(n) != i
        bi = ols_qr(X[m], y[m])
        D[i] = np.sum((X @ (b - bi)) ** 2) / (p * s2)
        loo[i] = y[i] - X[i] @ bi
        rss_i[i] = ((y[m] - X[m] @ bi) ** 2).sum()
    return D, loo, rss_i


xl = np.array([1, 2, 3, 4, 10.0])
for y5, D5, r5, hasbig in [(9.9, 0.2163, 0.194, False), (3.0, 16.961, -1.7175, True)]:
    yl = np.array([1.1, 1.9, 3.2, 3.9, y5]); Xl = np.c_[np.ones(5), xl]
    h, r, tt, D, fl = influence(Xl, yl)
    Db, loo, _ = brute_loo(Xl, yl)
    assert close(h, [0.38, 0.28, 0.22, 0.2, 0.92]) and close(D5, D[4], 1e-3) and close(r5, r[4], 1e-3)
    assert close(D, Db, 1e-10) and close(loo, fl["e"] / (1 - h), 1e-10)
    assert (D[4] > 1) == hasbig and (hasbig or D[4] < D[:4].max())   # on-trend leverage point: D5 below its neighbours' maximum
    print("Example 4.4  y5=%.1f  slope=%.3f  e5=%.3f  h5=%.2f  r5=%.3f  D5=%.3f  max D others=%.3f" % (y5, fl["b"][1], fl["e"][4], h[4], r[4], D[4], D[:4].max()))

rng = np.random.default_rng(4)
xs = np.r_[rng.uniform(0, 10, 27), 4.0, 5.5, 25.0]
ys = 1 + 0.8 * xs + 0.8 * rng.standard_normal(30)
ys[27] += 6.0; ys[28] -= 5.0                       # two planted outliers in y at ordinary x; index 29: high leverage, on trend
Xs = np.c_[np.ones(30), xs]
h, r, tt, D, fs = influence(Xs, ys)
Db, loo, rss_i = brute_loo(Xs, ys)
assert np.isclose(h.sum(), 2) and np.all((h >= 1 / 30 - 1e-12) & (h <= 1)) and np.allclose(h, 1 / 30 + (xs - xs.mean()) ** 2 / ((xs - xs.mean()) ** 2).sum())
assert np.allclose(D, Db) and np.allclose(loo, fs["e"] / (1 - h)) and np.allclose(rss_i, fs["rss"] - fs["e"] ** 2 / (1 - h))
assert set(np.argsort(-np.abs(tt))[:2]) == {27, 28} and h[29] > 2 * 2 / 30 and h[29] == h.max()
print("Part 4  planted outliers: |t_ext| =", np.abs(tt[[27, 28]]).round(2), " leverage point h = %.3f (2p/n = %.3f), D = %.3f" % (h[29], 4 / 30, D[29]))
fig, ax = plt.subplots(2, 2, figsize=(10, 7))
yh = Xs @ fs["b"]
ax[0, 0].scatter(yh, fs["e"]); ax[0, 0].axhline(0, c="k", lw=.5); ax[0, 0].set(title="residuals vs fitted", xlabel="fitted", ylabel="e")
(q, o), _ = stats.probplot(tt, dist=stats.t, sparams=(30 - 3,)); ax[0, 1].plot(q, o, "o"); ax[0, 1].plot(q, q, "k-", lw=.5); ax[0, 1].set(title="QQ plot of studentized residuals", xlabel="t quantile")
ax[1, 0].scatter(h, r**2); ax[1, 0].axvline(4 / 30, ls="--", c="gray"); ax[1, 0].set(title="leverage vs squared std. residual", xlabel="h_ii", ylabel="r_i^2")
for i in (27, 28, 29): ax[1, 0].annotate(str(i), (h[i], r[i] ** 2))
ax[1, 1].stem(D); ax[1, 1].set(title="Cook's distance", xlabel="observation")
fig.tight_layout(); fig.savefig(FIG / "diagnostics.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 5: heteroscedasticity: naive vs sandwich standard errors, and WLS (notes 3.6)

# %%
rng = np.random.default_rng(5)
n, R = 60, 4000
xh = np.sort(rng.uniform(0, 10, n)); Xh = np.c_[np.ones(n), xh]
sd = 0.15 * xh**2 + 0.2                    # planted pattern: sd grows like x^2
Ah = np.linalg.solve(Xh.T @ Xh, Xh.T); hh = np.einsum("ij,ji->i", Xh, Ah)
Yh = (Xh @ [1.0, 0.5])[:, None] + sd[:, None] * rng.standard_normal((n, R))
Bh = Ah @ Yh; Eh = Yh - Xh @ Bh
w1 = Ah[1]                                  # slope b_1 - beta_1 = sum_i w1_i eps_i
tqh = stats.t.ppf(0.975, n - 2)
v_naive = (Eh**2).sum(0) / (n - 2) * (w1**2).sum()
v_hc0 = (Eh**2).T @ w1**2                   # sandwich (3.37): sum_i w_i^2 e_i^2
v_hc3 = ((Eh / (1 - hh)[:, None]) ** 2).T @ w1**2
cov = {k: (np.abs(Bh[1] - 0.5) <= tqh * np.sqrt(v)).mean() for k, v in [("naive", v_naive), ("HC0", v_hc0), ("HC3", v_hc3)]}
true_var = (sd**2 * w1**2).sum()
Wt = 1 / sd**2
Bw = np.linalg.solve(Xh.T @ (Wt[:, None] * Xh), Xh.T * Wt) @ Yh         # WLS (3.38)
var_wls = np.linalg.inv(Xh.T @ (Wt[:, None] * Xh))[1, 1]
print("Part 5  slope coverage:", cov, " mean naive var %.5f vs true %.5f; Var OLS %.5f, Var WLS %.5f" % (v_naive.mean(), true_var, true_var, var_wls))
assert cov["naive"] < 0.92 and abs(cov["HC3"] - 0.95) < 0.015 and cov["HC0"] > cov["naive"] and var_wls < true_var
assert abs(Bw[1].var() / var_wls - 1) < 0.08 and abs(Bh[1].mean() - 0.5) < 0.01 and abs(Bh[1].var() / true_var - 1) < 0.08  # both unbiased; sandwich target is right
fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
ax[0].scatter(Xh @ Bh[:, 0], Eh[:, 0], s=10); ax[0].set(title="fan-shaped residuals", xlabel="fitted", ylabel="e")
ax[1].bar(list(cov), list(cov.values())); ax[1].axhline(.95, c="k", ls="--"); ax[1].set(ylim=(.8, 1), title="coverage of slope interval")
fig.tight_layout(); fig.savefig(FIG / "heteroscedasticity.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 6: multicollinearity: VIF and the SVD variance decomposition (notes 3.7, Example 4.5)

# %%
rng = np.random.default_rng(6)
rho = 0.99
L = np.linalg.cholesky(np.array([[1, rho], [rho, 1]]))
Xc = np.linalg.qr(rng.standard_normal((100, 2)))[0] @ L.T      # unit-norm columns, X'X = [[1, rho],[rho, 1]]
assert np.allclose(Xc.T @ Xc, [[1, rho], [rho, 1]])
U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
assert close(s**2, [1.99, 0.01], 1e-12) and close(np.linalg.cond(Xc), np.sqrt(199), 1e-9)
Vcov = np.linalg.inv(Xc.T @ Xc)                                # Var(b)/sigma^2
by_svd = sum(np.outer(Vt[i], Vt[i]) / s[i] ** 2 for i in range(2))  # (3.43)
vif = 1 / (1 - rho**2)
assert np.allclose(Vcov, by_svd) and close(Vcov[0, 0], vif, 1e-9) and close(vif, 50.2513, 5e-5) and close(np.sqrt(vif), 7.0888, 5e-5)
assert close(Vcov[0, 1] / Vcov[0, 0], -rho, 1e-12)
vsum, vdiff = np.array([1, 1]) @ Vcov @ [1, 1], np.array([1, -1]) @ Vcov @ [1, -1]
assert close([vsum, vdiff], [2 / (1 + rho), 2 / (1 - rho)], 1e-9) and close([vsum, vdiff], [1.0050, 200.0], 5e-4)
assert close([0.5 / s[0] ** 2, 0.5 / s[1] ** 2], [0.2513, 50.0], 5e-4)  # the two SVD terms of Var(b_1)/sigma^2
Rr = 20_000
Yc = (Xc @ [1.0, 1.0])[:, None] + rng.standard_normal((100, Rr))
Bc = np.linalg.solve(Xc.T @ Xc, Xc.T @ Yc)
assert abs(Bc[0].var() / vif - 1) < 0.03 and abs(np.corrcoef(Bc)[0, 1] + rho) < 0.005 and abs(Bc.mean(1) - 1).max() < 0.1  # unbiased, huge variance
print("Example 4.5  s =", s.round(4), " VIF = %.4f  se inflation = %.4f  Var(b1+b2)/s2 = %.4f  Var(b1-b2)/s2 = %.1f" % (vif, np.sqrt(vif), vsum, vdiff))
# VIF_j = 1/(1-R_j^2) = diag of the inverse correlation matrix, on a 4-predictor design.
Z = rng.standard_normal((200, 4)); Z[:, 3] = Z[:, 0] + Z[:, 1] + 0.3 * Z[:, 2] + 0.2 * rng.standard_normal(200)
Zc = Z - Z.mean(0)
vifs = []
for j in range(4):
    o = np.delete(Zc, j, 1); res = Zc[:, j] - o @ np.linalg.lstsq(o, Zc[:, j], rcond=None)[0]
    vifs.append((Zc[:, j] @ Zc[:, j]) / (res @ res))                    # 1/(1-R_j^2), via residual of x_j on the others
assert np.allclose(vifs, np.diag(np.linalg.inv(np.corrcoef(Z.T))))
print("Part 6  VIFs:", np.array(vifs).round(2))
# Variance of b_1 versus kappa(X)^2 as rho grows.
rhos = np.array([0, .5, .9, .99, .999]); kap, var1 = [], []
for r_ in rhos:
    Xr = np.linalg.qr(rng.standard_normal((100, 2)))[0] @ np.linalg.cholesky(np.array([[1, r_], [r_, 1]])).T
    kap.append(np.linalg.cond(Xr)); var1.append(np.linalg.inv(Xr.T @ Xr)[0, 0])
assert np.allclose(var1, 1 / (1 - rhos**2)) and np.allclose(np.array(kap) ** 2, (1 + rhos) / (1 - rhos))
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
ax[0].loglog(np.array(kap) ** 2, var1, "o-"); ax[0].set(xlabel="kappa(X)^2", ylabel="Var(b_1)/sigma^2", title="variance vs condition number")
th = np.linspace(0, 2 * np.pi, 200); circ = np.c_[np.cos(th), np.sin(th)]
ell = circ @ (Vt.T / s).T                                              # {b : (b)'X'X(b) = 1}: semi-axes 1/s_i along v_i
ax[1].plot(*ell.T); ax[1].plot(*(circ @ np.eye(2)).T, "k:", lw=.5); ax[1].set(aspect="equal", title="confidence ellipse, rho = 0.99", xlabel="b_1", ylabel="b_2")
fig.tight_layout(); fig.savefig(FIG / "collinearity.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 7: model selection: C_p, AIC, BIC, LOO shortcut, post-selection (notes 3.8)

# %%
rng = np.random.default_rng(7)
f_true = lambda x: np.sin(5 * x + 1)
dmax, sig7 = 16, 0.5
leg = lambda x, m: np.polynomial.legendre.legvander(x, m - 1)         # columns P_0..P_{m-1}: orthogonal-ish, well conditioned


def nested_fits(x, y, dmax):
    """For degrees 0..dmax: RSS_d, leverages h_ii(d), LOO errors, from ONE QR of the Legendre design (columns nested)."""
    Q = np.linalg.qr(leg(x, dmax + 1))[0]
    C = Q * (Q.T @ y)                                  # n x (dmax+1): column k = (q_k'y) q_k
    fits = np.cumsum(C, 1); E = y[:, None] - fits; Hd = np.cumsum(Q**2, 1)
    return (E**2).sum(0), Hd, E / (1 - Hd), fits


n7 = 30; x7 = np.linspace(-1, 1, n7); mu7 = f_true(x7)
Qf = np.linalg.qr(leg(x7, dmax + 1))[0]
bias2 = np.array([((mu7 - (Qf[:, :d + 1] * (Qf[:, :d + 1].T @ mu7)).sum(1)) ** 2).sum() for d in range(dmax + 1)])  # ||(I-H)mu||^2 by degree
pp = np.arange(1, dmax + 2)
Rm = 4000
Y7 = mu7[:, None] + sig7 * rng.standard_normal((n7, Rm)); Ynew = mu7[:, None] + sig7 * rng.standard_normal((n7, Rm))
rss_mc, epe_mc = np.zeros(dmax + 1), np.zeros(dmax + 1)
for d in range(dmax + 1):
    Hq = Qf[:, :d + 1] @ Qf[:, :d + 1].T; Fit = Hq @ Y7
    rss_mc[d] = ((Y7 - Fit) ** 2).sum(0).mean(); epe_mc[d] = ((Ynew - Fit) ** 2).sum(0).mean()
epe = n7 * sig7**2 + bias2 + pp * sig7**2                  # (3.46)
assert np.allclose(rss_mc, bias2 + (n7 - pp) * sig7**2, rtol=0.03, atol=0.02) and np.allclose(epe_mc, epe, rtol=0.03)
cp_true = rss_mc + 2 * pp * sig7**2                         # E[RSS + 2 p sigma^2] = EPE - n sigma^2 + n sigma^2 ... see (3.47)
assert np.allclose(cp_true, epe, rtol=0.03)                 # optimism = 2 p sigma^2
print("Part 7  expected prediction error by degree:", epe.round(2), " best degree", epe.argmin())
# LOO shortcut = brute force, for one nested fit.
yq = Y7[:, 0]; rs, Hd, loo_e, _ = nested_fits(x7, yq, dmax)
d_ = 4; Xd = leg(x7, d_ + 1)
brute = np.array([yq[i] - leg(x7[i:i + 1], d_ + 1) @ ols_qr(np.delete(Xd, i, 0), np.delete(yq, i)) for i in range(n7)]).ravel()
assert np.allclose(brute, loo_e[:, d_])
# Selected degree versus n.
ns = [20, 40, 80, 160, 320, 640, 1280]; sel = {k: [] for k in ("AIC", "BIC", "LOO")}
for nn in ns:
    xx_ = np.linspace(-1, 1, nn); mm = f_true(xx_); picks = {k: [] for k in sel}
    for _ in range(100):
        yy_ = mm + sig7 * rng.standard_normal(nn)
        rs, Hd, le, _ = nested_fits(xx_, yy_, dmax)
        picks["AIC"].append(np.argmin(nn * np.log(rs / nn) + 2 * pp)); picks["BIC"].append(np.argmin(nn * np.log(rs / nn) + pp * np.log(nn)))
        picks["LOO"].append(np.argmin((le**2).mean(0)))
    for k in sel: sel[k].append(np.mean(picks[k]))
print("Part 7  mean selected degree vs n =", ns, "\n        ", {k: np.round(v, 2) for k, v in sel.items()})
assert sel["LOO"][-1] > sel["LOO"][1] + 0.5 and sel["BIC"][-1] > sel["BIC"][1] + 0.5 and all(b <= a + 1e-9 for a, b in zip(sel["AIC"], sel["BIC"]))
assert sel["AIC"][0] > 12 > sel["LOO"][0]   # p close to n: AIC's penalty 2p is too light and it overfits wildly; LOO does not
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot(range(dmax + 1), bias2 / n7, label="bias^2 / n"); ax[0].plot(range(dmax + 1), pp * sig7**2 / n7, label="p sigma^2 / n (variance)")
ax[0].plot(range(dmax + 1), epe / n7 - sig7**2, "k", label="reducible error"); ax[0].plot(range(dmax + 1), rss_mc / n7, "--", label="training MSE"); ax[0].set(yscale="log", xlabel="degree", title="bias-variance, n = 30"); ax[0].legend(fontsize=7)
for k, v in sel.items(): ax[1].semilogx(ns, v, "o-", label=k)
ax[1].set(xlabel="n", ylabel="mean selected degree", title="selected degree vs n"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "model_selection.png", dpi=110); plt.close(fig)

# Post-selection: pure noise, pick the 5 best of 50 predictors by |corr|, then refit as if chosen in advance.
n, P, k, Rps = 100, 50, 5, 500
sig_frac, f_rej = [], 0
for _ in range(Rps):
    Xn = rng.standard_normal((n, P)); yn = rng.standard_normal(n)
    top = np.argsort(-np.abs(Xn.T @ yn) / np.linalg.norm(Xn, axis=0))[:k]
    ft = ols_fit(np.c_[np.ones(n), Xn[:, top]], yn)
    sig_frac.append((ft["p"][1:] < 0.05).mean()); f_rej += ft["pF"] < 0.05
print("Part 7  post-selection (pure noise): %.2f of selected coefficients 'significant' at 5%%; overall F rejects in %.2f of runs" % (np.mean(sig_frac), f_rej / Rps))
assert np.mean(sig_frac) > 0.3 and f_rej / Rps > 0.9   # nominal: 0.05 and 0.05

# %% [markdown]
# ## Part 8: ridge (notes 3.9, Example 4.6)

# %%
rng = np.random.default_rng(8)
# Example 4.6: singular values (10, 1, 0.1), lambda = 1.
Xr = make_design(30, np.array([10, 1, 0.1]), rng)
U, s, Vt = np.linalg.svd(Xr, full_matrices=False); assert close(s, [10, 1, 0.1], 1e-12)
lam = 1.0
fac = s**2 / (s**2 + lam)
assert close(fac, [0.990099, 0.5, 0.009901], 5e-7) and close(fac.sum(), 1.5, 1e-12)           # shrinkage factors; effective df
assert close(s / (s**2 + lam), [0.0990, 0.5, 0.0990], 5e-5) and close(1 / s, [0.1, 1, 10])     # ridge vs OLS multiplier on u_i'y
yr = Xr @ [1, -1, 2.0] + rng.standard_normal(30)
b_dir = np.linalg.solve(Xr.T @ Xr + lam * np.eye(3), Xr.T @ yr)
b_svd = Vt.T @ (s / (s**2 + lam) * (U.T @ yr))                                                  # (3.41)
b_aug = np.linalg.lstsq(np.r_[Xr, np.sqrt(lam) * np.eye(3)], np.r_[yr, np.zeros(3)], rcond=None)[0]  # ridge = OLS on augmented data
assert np.allclose(b_dir, b_svd) and np.allclose(b_dir, b_aug) and np.isclose(np.trace(Xr @ np.linalg.solve(Xr.T @ Xr + lam * np.eye(3), Xr.T)), 1.5)
alpha, sg = np.array([1.0, 1, 1]), 1.0                                                          # coordinates a = V'beta, sigma^2 = 1
mse_ridge = lambda l: np.sum((sg**2 * s**2 + l**2 * alpha**2) / (s**2 + l) ** 2)               # (3.44)
assert close([mse_ridge(0.0), mse_ridge(1.0)], [101.01, 1.5], 5e-3)
assert mse_ridge(1e-3) < mse_ridge(0.0) and np.isclose((mse_ridge(1e-9) - mse_ridge(0.0)) / 1e-9, -2 * np.sum(1 / s**4), rtol=1e-3)  # Hoerl-Kennard slope at 0
Bt = np.linalg.solve(Xr.T @ Xr + lam * np.eye(3), Xr.T @ ((Xr @ (Vt.T @ alpha))[:, None] + rng.standard_normal((30, 4000))))
mse_mc = np.mean(np.sum((Bt - (Vt.T @ alpha)[:, None]) ** 2, 0))
assert abs(mse_mc / 1.5 - 1) < 0.05
print("Example 4.6  factors", fac.round(6), " df =", fac.sum(), " MSE(OLS)=%.2f  MSE(ridge, lam=1)=%.4f (MC %.3f)" % (mse_ridge(0.0), mse_ridge(1.0), mse_mc))
# Orthonormal case: lambda* = p sigma^2/||beta||^2 and MSE* = p sigma^2/(1 + p sigma^2/||beta||^2).
pb, bb = 5, np.array([1.0, 0.5, 0.5, 0.2, 0.1]); l_star = pb * 0.25 / (bb @ bb)
mse_o = lambda l: (pb * 0.25 + l**2 * (bb @ bb)) / (1 + l) ** 2
assert close(mse_o(l_star), pb * 0.25 / (1 + l_star), 1e-12) and mse_o(l_star) < mse_o(0.0) and optimize.minimize_scalar(mse_o, bounds=(0, 10), method="bounded").x - l_star < 1e-4

# Correlated design, p = 10, n = 40: coefficient path, LOO-CV via the shortcut, MSE vs lambda.
n, P = 40, 10
Sig = 0.8 ** np.abs(np.subtract.outer(range(P), range(P)))
Xp = rng.standard_normal((n, P)) @ np.linalg.cholesky(Sig).T
Xp = (Xp - Xp.mean(0)) / Xp.std(0)                              # standardize (notes 3.9): ridge penalizes in these units
bt = np.r_[3, -2, 1.5, np.zeros(7)] * 1.0
yp = Xp @ bt + 3 * rng.standard_normal(n); yp -= yp.mean()      # center y: intercept unpenalized
U, s, Vt = np.linalg.svd(Xp, full_matrices=False)
lams = np.logspace(-2, 3, 60)
path = np.array([Vt.T @ (s / (s**2 + l) * (U.T @ yp)) for l in lams])
loo = np.array([np.mean(((yp - Xp @ pth) / (1 - hat_diag(Xp, l))) ** 2) for pth, l in zip(path, lams)])
i0 = 20; l0 = lams[i0]
brute = np.array([yp[i] - Xp[i] @ np.linalg.solve(np.delete(Xp, i, 0).T @ np.delete(Xp, i, 0) + l0 * np.eye(P), np.delete(Xp, i, 0).T @ np.delete(yp, i)) for i in range(n)])
assert np.allclose(np.mean(brute**2), loo[i0])                   # LOO shortcut is exact for ridge
l_cv = lams[loo.argmin()]
al = Vt @ bt; mse_curve = np.array([np.sum((9 * s**2 + l**2 * al**2) / (s**2 + l) ** 2) for l in lams])  # sigma^2 = 9
assert mse_curve.min() < mse_curve[0] * 0.7 and 0 < mse_curve.argmin() < len(lams) - 1
print("Part 8  OLS MSE %.2f; ridge min MSE %.2f at lambda %.2f; LOO-CV picks lambda %.2f (MSE there %.2f)" % (
    np.sum(9 / s**2), mse_curve.min(), lams[mse_curve.argmin()], l_cv, np.interp(np.log(l_cv), np.log(lams), mse_curve)))
# p > n: lambda -> 0 gives the minimum-norm solution X^+ y (LA-9 Theorem 3.10).
Xw = rng.standard_normal((10, 25)); yw = rng.standard_normal(10)
b_mn = np.linalg.pinv(Xw) @ yw
assert np.allclose(Xw.T @ np.linalg.solve(Xw @ Xw.T + 1e-9 * np.eye(10), yw), b_mn, atol=1e-6) and np.allclose(Xw @ b_mn, yw)  # interpolates; dual form of ridge
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].semilogx(lams, path); ax[0].axvline(l_cv, c="k", ls=":"); ax[0].set(xlabel="lambda", ylabel="coefficient", title="ridge path (dotted: LOO-CV)")
ax[1].loglog(lams, mse_curve); ax[1].axhline(np.sum(9 / s**2), c="r", ls="--", label="OLS"); ax[1].set(xlabel="lambda", ylabel="E||b - beta||^2", title="ridge MSE: minimum below OLS"); ax[1].legend()
fig.tight_layout(); fig.savefig(FIG / "ridge.png", dpi=110); plt.close(fig)

# %% [markdown]
# ## Part 9: the lasso by coordinate descent (notes 3.10, Example 4.7)

# %%
def soft(z, t):
    """S_t(z) = sign(z) (|z| - t)_+, the prox of t|.| (3.50)."""
    return np.sign(z) * np.maximum(np.abs(z) - t, 0.0)


def lasso_cd(X, y, lam, b=None, tol=1e-9, max_sweeps=5000):
    """Minimize 0.5||y - Xb||^2 + lam ||b||_1 by cycling b_j <- S_lam(x_j'r_j)/||x_j||^2, (3.52), r kept up to date."""
    n, p = X.shape
    b = np.zeros(p) if b is None else b.copy()
    r = y - X @ b; cn = (X**2).sum(0)
    for _ in range(max_sweeps):
        delta = 0.0
        for j in range(p):
            bj = soft(X[:, j] @ r + cn[j] * b[j], lam) / cn[j]
            if bj != b[j]:
                r += X[:, j] * (b[j] - bj); delta = max(delta, abs(bj - b[j])); b[j] = bj
        if delta < tol: break
    return b


rng = np.random.default_rng(9)
# Example 4.7: orthonormal design, so the answer is soft-thresholding of z = X'y.
Qo = np.linalg.qr(rng.standard_normal((20, 5)))[0]
z = np.array([3.0, -1.5, 0.4, -0.2, 2.2])
w = rng.standard_normal(20); w -= Qo @ (Qo.T @ w)                  # component orthogonal to C(X)
yo = Qo @ z + w; assert np.allclose(Qo.T @ yo, z)
bl = lasso_cd(Qo, yo, 1.0)
assert close(bl, [2.0, -0.5, 0, 0, 1.2], 1e-7) and np.allclose(bl, soft(z, 1.0)) and close(z / 2, [1.5, -0.75, 0.2, -0.1, 1.1])
print("Example 4.7  soft(z, 1) =", soft(z, 1.0), " ridge(lambda=1) =", z / 2)
# Slow reference: split b = b+ - b-, smooth bound-constrained problem solved by L-BFGS-B; plus KKT conditions.
n, P = 30, 6
Xs6 = rng.standard_normal((n, P)); Xs6[:, 1] = Xs6[:, 0] + 0.3 * rng.standard_normal(n); ys6 = Xs6 @ [2, 0, -1, 0, 0, 0.5] + rng.standard_normal(n)
lam6 = 8.0
def obj(v): bv = v[:P] - v[P:]; r_ = ys6 - Xs6 @ bv; return 0.5 * r_ @ r_ + lam6 * v.sum(), np.r_[-Xs6.T @ r_ + lam6, Xs6.T @ r_ + lam6]
ref = optimize.minimize(obj, np.zeros(2 * P), jac=True, bounds=[(0, None)] * (2 * P), method="L-BFGS-B", options=dict(ftol=1e-15, gtol=1e-11, maxiter=10_000))
b_ref = ref.x[:P] - ref.x[P:]; b_cd = lasso_cd(Xs6, ys6, lam6)
g = Xs6.T @ (ys6 - Xs6 @ b_cd)
assert np.allclose(b_cd, b_ref, atol=1e-5) and np.allclose(g[b_cd != 0], lam6 * np.sign(b_cd[b_cd != 0]), atol=1e-6) and np.all(np.abs(g[b_cd == 0]) <= lam6 + 1e-9)
lmax = np.abs(Xs6.T @ ys6).max()
assert not lasso_cd(Xs6, ys6, lmax * 1.0001).any() and (lasso_cd(Xs6, ys6, lmax * 0.99) != 0).sum() == 1
print("Part 9  CD matches L-BFGS-B reference (max diff %.1e), KKT holds, lambda_max = %.3f; nonzeros at lam=8: %d" % (np.abs(b_cd - b_ref).max(), lmax, (b_cd != 0).sum()))

def lasso_path(X, y, lams):
    b, out = None, []
    for l in lams: b = lasso_cd(X, y, l, b, tol=1e-7); out.append(b)
    return np.array(out)


def cv_lasso(X, y, lams, K=5):
    idx = np.arange(len(y)) % K; err = np.zeros(len(lams))
    for k in range(K):
        pth = lasso_path(X[idx != k], y[idx != k], lams); err += ((y[idx == k][:, None] - X[idx == k] @ pth.T) ** 2).sum(0)
    return err / len(y)


# Sparse recovery: n = 100, p = 50, five nonzero coefficients.
n, P = 100, 50
Xl = rng.standard_normal((n, P)); Xl = (Xl - Xl.mean(0)) / Xl.std(0)
bs = np.zeros(P); bs[[3, 11, 20, 31, 44]] = [3, -2, 2.5, -1.5, 2]
yl = Xl @ bs + rng.standard_normal(n); yl -= yl.mean()
lams9 = np.geomspace(np.abs(Xl.T @ yl).max(), 1e-2 * np.abs(Xl.T @ yl).max(), 40)
pth = lasso_path(Xl, yl, lams9); cvm = cv_lasso(Xl, yl, lams9); ib = cvm.argmin()
b_l = pth[ib]; b_ols = ols_qr(Xl, yl)
supp = set(np.flatnonzero(b_l)); true = set(np.flatnonzero(bs))
print("Part 9  CV-lasso: lambda = %.2f, support size %d, true support recovered: %s, false positives %d; ||b-beta||: lasso %.2f, OLS %.2f" % (
    lams9[ib], len(supp), true <= supp, len(supp - true), np.linalg.norm(b_l - bs), np.linalg.norm(b_ols - bs)))
assert true <= supp and np.linalg.norm(b_l - bs) < 0.5 * np.linalg.norm(b_ols - bs)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].semilogx(lams9, pth, lw=1); ax[0].axvline(lams9[ib], c="k", ls=":"); ax[0].set(xlabel="lambda", ylabel="coefficient", title="lasso path, p = 50 (dotted: CV)")
zz = np.linspace(-3, 3, 200); ax[1].plot(zz, zz, "k:", label="OLS"); ax[1].plot(zz, soft(zz, 1), label="lasso (soft)"); ax[1].plot(zz, zz / 2, label="ridge"); ax[1].legend(); ax[1].set(xlabel="z", title="orthonormal design, lambda = 1")
fig.tight_layout(); fig.savefig(FIG / "lasso.png", dpi=110); plt.close(fig)

# Ridge vs lasso: sparse truth / dense truth / correlated predictors; lambda chosen on a validation set.
def compare(kind, reps=12):
    out = {"ridge": [], "lasso": [], "missed": [], "false+": []}
    for _ in range(reps):
        Sg = 0.9 ** np.abs(np.subtract.outer(range(P), range(P))) if "correlated" in kind else np.eye(P)
        Xa = rng.standard_normal((300, P)) @ np.linalg.cholesky(Sg).T; Xa = (Xa - Xa.mean(0)) / Xa.std(0)
        bt = 0.45 * rng.standard_normal(P) if kind.startswith("dense") else bs.copy()
        ya = Xa @ bt + 2 * rng.standard_normal(300); ya -= ya.mean()
        Xt, yt_, Xv, yv = Xa[:100], ya[:100], Xa[100:200], ya[100:200]
        U_, s_, Vt_ = np.linalg.svd(Xt, full_matrices=False)
        br = [Vt_.T @ (s_ / (s_**2 + l) * (U_.T @ yt_)) for l in np.geomspace(1e-1, 1e4, 40)]
        out["ridge"].append(np.linalg.norm(br[int(np.argmin([np.sum((yv - Xv @ b_) ** 2) for b_ in br]))] - bt))
        ll = np.geomspace(np.abs(Xt.T @ yt_).max(), 1e-2 * np.abs(Xt.T @ yt_).max(), 20); bl_ = lasso_path(Xt, yt_, ll)
        bhat = bl_[int(np.argmin(((yv[:, None] - Xv @ bl_.T) ** 2).sum(0)))]
        out["lasso"].append(np.linalg.norm(bhat - bt)); out["missed"].append(np.sum((bhat == 0) & (bs != 0))); out["false+"].append(np.sum((bhat != 0) & (bs == 0)))
    return {k: np.mean(v) for k, v in out.items()}
res = {k: compare(k) for k in ("sparse", "dense", "sparse-correlated", "dense-correlated")}
print("Part 9  ||b - beta|| (mean of 12 runs) and lasso support errors:")
for k, d in res.items(): print("        %-18s" % k, {m: round(float(v), 2) for m, v in d.items()})
assert res["sparse"]["lasso"] < 0.7 * res["sparse"]["ridge"] and res["dense"]["ridge"] < res["dense"]["lasso"]

# %% [markdown]
# ## Part 10: numbers used in exercises and solutions (all exercises.md / solutions.md figures are asserted here)

# %%
# A1/A2: four points.
xa, ya = np.array([0.0, 1, 2, 3]), np.array([1.0, 2, 4, 5]); Xa = np.c_[np.ones(4), xa]; fa = ols_fit(Xa, ya)
assert close(fa["b"], [0.9, 1.4]) and close(fa["e"], [0.1, -0.3, 0.3, -0.1]) and close([fa["rss"], fa["tss"], fa["R2"]], [0.2, 10, 0.98])
assert close(fa["h"], [0.7, 0.3, 0.3, 0.7]) and close(fa["e"] / (1 - fa["h"]), [0.3333, -0.4286, 0.4286, -0.3333], 5e-5)
press = ((fa["e"] / (1 - fa["h"])) ** 2).sum(); assert close([press, press / 4], [0.5896, 0.1474], 5e-5)
assert close(Xa.T @ fa["e"], [0, 0], 1e-12)
# A3.
s2_ = 45 / 17; t17 = stats.t.ppf(0.975, 17)
assert close([s2_, np.sqrt(s2_), t17, t17 * np.sqrt(s2_ * 0.1), t17 * np.sqrt(s2_ * 1.1)], [2.6471, 1.6270, 2.1098, 1.0855, 3.6002], 5e-5)
# A4/A5/A6.
za = np.array([2.5, -0.7, 1.1, -3.2]); assert close(soft(za, 1.0), [1.5, 0, 0.1, -2.2]) and close(za / 2, [1.25, -0.35, 0.55, -1.6])
sa = np.array([4, 2, 0.5]); assert close(sa**2 / (sa**2 + 1), [0.9412, 0.8, 0.2], 5e-5) and close((sa**2 / (sa**2 + 1)).sum(), 1.9412, 5e-5)
Xd = np.array([[1, 1, 0], [1, 1, 0], [1, 0, 1], [1, 0, 1.0]]); assert np.linalg.matrix_rank(Xd) == 2 and np.allclose(Xd @ [1, -1, -1], 0)
yd = np.array([3.0, 5, 6, 8]); bd = np.linalg.pinv(Xd) @ yd
assert close(bd, [3.6667, 0.3333, 3.3333], 5e-5) and close(Xd @ bd, [4, 4, 7, 7]) and close(np.linalg.norm(bd), 4.9666, 5e-5)  # min-norm: beta0 = 11/3
# B3: short coefficient = long coefficient + gamma * (omitted coefficient), Example 4.2 data, x2 omitted.
Xlong, Xshort = X42[:, :3], X42[:, :2]
bl_, bs_ = ols_qr(Xlong, y42), ols_qr(Xshort, y42)
gam = ols_qr(Xshort, x2)[1]
assert close(bs_[1], bl_[1] + gam * bl_[2], 1e-12) and close([bs_[1], bl_[1], gam, bl_[2]], [1.0440, 0.3103, 0.8214, 0.8933], 5e-5)
print("B3 numbers:", bs_[1], bl_[1], gam, bl_[2])
# B5: -2 log lambda for n=30, p=4, q=2, F=3.5.
assert close(30 * np.log(1 + 2 * 3.5 / 26), 7.1523, 5e-4) and close(stats.f.sf(3.5, 2, 26), 0.0451, 5e-4) and close(stats.chi2.sf(7.1523, 2), 0.0280, 5e-4)
# B6: intercept-only, y = (1,2,3,6).
yb = np.array([1.0, 2, 3, 6]); eb = yb - yb.mean()
assert close([eb @ eb / 3 / 4, eb @ eb / 16, (eb**2 / 0.5625).sum() / 16], [1.1667, 0.875, 1.5556], 5e-5)
# C2: two identical columns: lasso solution not unique, ridge splits evenly.
xc = np.array([1.0, 2, 3, 4]); Xc2 = np.c_[xc, xc]; yc2 = np.array([1.5, 2.5, 3.0, 4.5])
b_a = lasso_cd(Xc2, yc2, 1.0, np.array([5.0, 0.0])); b_b = lasso_cd(Xc2, yc2, 1.0, np.array([0.0, 5.0]))
lo = lambda b: 0.5 * np.sum((yc2 - Xc2 @ b) ** 2) + np.abs(b).sum()
assert np.isclose(b_a.sum(), b_b.sum()) and np.isclose(lo(b_a), lo(b_b)) and not np.allclose(b_a, b_b)
rb = np.linalg.solve(Xc2.T @ Xc2 + np.eye(2), Xc2.T @ yc2); assert np.isclose(rb[0], rb[1])
print("C2 numbers: lasso starts (5,0)->", b_a.round(4), " (0,5)->", b_b.round(4), " ridge", rb.round(4), " sum", b_a.sum().round(4))
# C3: mean-shift model on Example 4.4 data (y5 = 3.0): coefficient on the indicator = LOO error, t = external studentized residual.
yl4 = np.array([1.1, 1.9, 3.2, 3.9, 3.0]); Xl4 = np.c_[np.ones(5), xl]; di = np.eye(5)[:, 4]
Xm = np.c_[Xl4, di]; bm = np.linalg.solve(Xm.T @ Xm, Xm.T @ yl4); em = yl4 - Xm @ bm; s2m = em @ em / (5 - 3)
tm = bm[2] / np.sqrt(s2m * np.linalg.inv(Xm.T @ Xm)[2, 2])
h4, r4, t4, D4, f4 = influence(Xl4, yl4)
assert close(bm[2], f4["e"][4] / (1 - h4[4]), 1e-10) and close(tm, t4[4], 1e-10) and close(bm[2], -6.8, 1e-9)
print("C3 numbers: indicator coefficient %.4f, t = %.4f (external studentized residual %.4f), 2 df" % (bm[2], tm, t4[4]))

# %% [markdown]
# ## Your turn
#
# 1. **WLS** (notes 3.6): Part 5 already compares OLS and WLS variances. Extend it: if the weights are estimated
#    (fit OLS, regress log e_i^2 on log x_i, use the fitted variances), does coverage survive?
# 2. **Interaction and FWL** (exercise B8): reference solution below.
# 3. **Elastic net** (exercise C4): add `l2` to `lasso_cd` (the update becomes S_lam(x_j'r_j)/(||x_j||^2 + l2)); rerun
#    the duplicated-column example and the correlated comparison of Part 9.
# 4. **Omitted-variable bias** (notes 3.8): reference solution below.

# %%
# Reference solution 1 (B8): y ~ 1 + x + g + x*g. The interaction coefficient = FWL slope on the residualized product column.
rng = np.random.default_rng(10)
n = 80; xg = rng.uniform(0, 5, n); gg = (rng.uniform(size=n) < 0.5).astype(float)
yg = 1 + 0.5 * xg + 2 * gg + 0.8 * xg * gg + rng.standard_normal(n)
Xg = np.c_[np.ones(n), xg, gg, xg * gg]; bg = ols_qr(Xg, yg)
M = np.eye(n) - Xg[:, :3] @ np.linalg.solve(Xg[:, :3].T @ Xg[:, :3], Xg[:, :3].T)
assert np.isclose((M @ Xg[:, 3]) @ (M @ yg) / ((M @ Xg[:, 3]) @ (M @ Xg[:, 3])), bg[3])
sl0 = ols_qr(Xg[gg == 0][:, :2], yg[gg == 0])[1]; sl1 = ols_qr(Xg[gg == 1][:, :2], yg[gg == 1])[1]
assert np.isclose(bg[1], sl0) and np.isclose(bg[1] + bg[3], sl1)             # interaction = difference of the two within-group slopes
print("B8  interaction coefficient %.4f = slope(g=1) - slope(g=0) = %.4f" % (bg[3], sl1 - sl0))
# Reference solution 2: omitted-variable bias, E[b_short] = b_1 + gamma b_2 (3.45), by simulation.
n, Rr = 50, 20_000
xa_ = rng.standard_normal(n); xb_ = 0.6 * xa_ + 0.8 * rng.standard_normal(n); xb_ -= xb_.mean(); xa_ -= xa_.mean()
Xfull = np.c_[np.ones(n), xa_, xb_]; bt = np.array([1.0, 2.0, 1.5])
Yo = (Xfull @ bt)[:, None] + rng.standard_normal((n, Rr))
b_short = np.linalg.solve(Xfull[:, :2].T @ Xfull[:, :2], Xfull[:, :2].T @ Yo)[1]
gamma = (xa_ @ xb_) / (xa_ @ xa_)
print("OVB  mean of short-regression slope %.4f; theory b_1 + gamma b_2 = %.4f (gamma = %.4f)" % (b_short.mean(), bt[1] + gamma * bt[2], gamma))
assert abs(b_short.mean() - (bt[1] + gamma * bt[2])) < 0.01 and abs(gamma * bt[2]) > 0.5

# %%
if __name__ == "__main__":
    print("All checks passed in %.1f s." % (time.time() - T0))
