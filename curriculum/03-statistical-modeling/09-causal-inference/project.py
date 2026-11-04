# %% [markdown]
# # ST-9 Project: a causal simulation library
#
# A small linear-Gaussian structural causal model (SCM) class with interventions,
# the trek rule (Wright's path rules) and d-separation, then every trap and design in
# the notes: confounding vs randomization, standardization and IPW, colliders,
# mediators, M-bias, the backdoor criterion, Fisher's randomization test, Neyman's
# variance, blocking, factorial designs, IV and weak instruments, regression
# discontinuity, difference-in-differences, sensitivity analysis, Granger
# "causality" from a common cause, and a predictive model that gives bad policy
# advice. The last part prints every number quoted in the exercises' solutions.
# Run: `MPLBACKEND=Agg python3 project.py`.

# %%
import itertools
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
rng = np.random.default_rng(9)
T0 = time.time()


def ols(y, *cols):
    """Least-squares coefficients of y on [1, cols...] (ST-3); returns (beta, residuals)."""
    X = np.column_stack([np.ones(len(y)), *cols])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return b, y - X @ b


# %% [markdown]
# ## Part 1: the linear-Gaussian SCM
# Each node is X_j = sum_{parents k} B[j,k] X_k + eps_j with eps_j ~ N(0, omega_j).
# In vector form X = B X + eps, so X = (I - B)^{-1} eps and
# Sigma = (I - B)^{-1} Omega (I - B)^{-T}  (notes eq. 3.18). For a DAG, B is nilpotent,
# so (I - B)^{-1} = I + B + B^2 + ...: entry (i,t) is the sum over directed paths t -> i
# of products of coefficients. That is the trek rule (eq. 3.19).

# %%
class SCM:
    """Linear-Gaussian SCM on named nodes given in topological order."""

    def __init__(self, nodes, edges, omega=None):
        self.nodes = list(nodes)
        self.ix = {v: i for i, v in enumerate(self.nodes)}
        p = len(self.nodes)
        self.B = np.zeros((p, p))
        for (par, ch), c in edges.items():
            assert self.ix[par] < self.ix[ch], "nodes must be in topological order"
            self.B[self.ix[ch], self.ix[par]] = c
        self.omega = np.ones(p) if omega is None else np.array([omega.get(v, 1.0) for v in self.nodes])

    def parents(self, v):
        return [u for u in self.nodes if self.B[self.ix[v], self.ix[u]] != 0]

    def children(self, v):
        return [u for u in self.nodes if self.B[self.ix[u], self.ix[v]] != 0]

    def descendants(self, v):
        out, stack = set(), [v]
        while stack:
            for c in self.children(stack.pop()):
                if c not in out:
                    out.add(c); stack.append(c)
        return out

    def sample(self, n, do=None):
        """n draws; do={node: value} replaces that node's equation (the do-operator, eq. 3.13)."""
        do = do or {}
        X = np.zeros((n, len(self.nodes)))
        for j, v in enumerate(self.nodes):
            if v in do:
                X[:, j] = do[v]
            else:
                X[:, j] = X @ self.B[j] + np.sqrt(self.omega[j]) * rng.standard_normal(n)
        return {v: X[:, self.ix[v]] for v in self.nodes}

    def cov(self):
        """Sigma = (I - B)^{-1} Omega (I - B)^{-T}, eq. (3.18)."""
        M = np.linalg.inv(np.eye(len(self.nodes)) - self.B)
        return M @ np.diag(self.omega) @ M.T

    def directed_paths(self, s, t):
        if s == t:
            return [[s]]
        return [[s] + p for c in self.children(s) for p in self.directed_paths(c, t)]

    def path_weight(self, path):
        return np.prod([self.B[self.ix[b], self.ix[a]] for a, b in zip(path, path[1:])])

    def trek_cov(self, i, j):
        """Trek rule: Cov(X_i, X_j) = sum over tops t and directed paths t~>i, t~>j of
        (product of coefficients on both paths) * omega_t   (eq. 3.19)."""
        return sum(self.path_weight(p) * self.path_weight(q) * self.omega[self.ix[t]]
                   for t in self.nodes for p in self.directed_paths(t, i) for q in self.directed_paths(t, j))

    def pop_coef(self, y, xs):
        """Population regression coefficients of y on xs (no intercept needed: all means 0)."""
        S = self.cov(); a = [self.ix[x] for x in xs]
        return np.linalg.solve(S[np.ix_(a, a)], S[a, self.ix[y]])

    def skeleton_paths(self, s, t, path=None):
        path = path or [s]
        if path[-1] == t:
            yield path; return
        u = path[-1]
        for w in self.parents(u) + self.children(u):
            if w not in path:
                yield from self.skeleton_paths(s, t, path + [w])

    def blocked(self, path, Z):
        """d-separation rule (Def. 3.6): a path is blocked by Z if some interior node is a
        non-collider in Z, or a collider with neither it nor any descendant in Z."""
        Z = set(Z)
        for a, m, b in zip(path, path[1:], path[2:]):
            collider = m in self.children(a) and m in self.children(b)
            if collider and not ({m} | self.descendants(m)) & Z:
                return True
            if not collider and m in Z:
                return True
        return False

    def d_separated(self, x, y, Z):
        return all(self.blocked(p, Z) for p in self.skeleton_paths(x, y))

    def backdoor_ok(self, a, y, Z):
        """Backdoor criterion (Thm 3.7): no descendant of a in Z, and every path a <- ... y blocked."""
        if set(Z) & self.descendants(a):
            return False
        return all(self.blocked(p, Z) for p in self.skeleton_paths(a, y) if p[1] in self.parents(a))


# Wright's rules vs matrix formula vs sample covariance on the collider DAG of Example 4.5.
col = SCM(["A", "Y", "C"], {("A", "Y"): 0.5, ("A", "C"): 1.0, ("Y", "C"): 1.0})
S = col.cov()
for u, v in itertools.combinations_with_replacement(col.nodes, 2):
    assert np.isclose(col.trek_cov(u, v), S[col.ix[u], col.ix[v]])
d = col.sample(200_000)
Ssamp = np.cov(np.vstack([d[v] for v in col.nodes]))
assert np.allclose(Ssamp, S, atol=0.05)
print("Ex 4.5 covariances (A,Y,C):\n", np.round(S, 4))
bAC = col.pop_coef("Y", ["A", "C"])
print(f"Ex 4.5: Y on A alone = {col.pop_coef('Y', ['A'])[0]:.3f};  Y on (A, C) = {np.round(bAC, 4)}")
assert np.allclose(bAC, [-0.25, 0.5])
b_s, _ = ols(d["Y"], d["A"], d["C"])
assert abs(b_s[1] + 0.25) < 0.02
berk = SCM(["A", "Y", "C"], {("A", "C"): 1.0, ("Y", "C"): 1.0})
assert np.isclose(berk.pop_coef("Y", ["A", "C"])[0], -0.5)
print("Pure collider (no A->Y): coefficient on A after adjusting C =", berk.pop_coef("Y", ["A", "C"])[0])

# %% [markdown]
# ## Part 2: Simpson's paradox (Example 4.1, kidney stones)

# %%
rec = {("open", "small"): (81, 87), ("open", "large"): (192, 263),
       ("pcnl", "small"): (234, 270), ("pcnl", "large"): (55, 80)}
nsize = {s: sum(rec[(t, s)][1] for t in ("open", "pcnl")) for s in ("small", "large")}
N = sum(nsize.values())
for t in ("open", "pcnl"):
    crude = sum(rec[(t, s)][0] for s in nsize) / sum(rec[(t, s)][1] for s in nsize)
    strata = [rec[(t, s)][0] / rec[(t, s)][1] for s in nsize]
    std = sum(rec[(t, s)][0] / rec[(t, s)][1] * nsize[s] / N for s in nsize)
    print(f"Ex 4.1 {t}: small {strata[0]:.4f} large {strata[1]:.4f} crude {crude:.4f} standardized {std:.4f}")
assert np.isclose(sum(rec[("open", s)][0] for s in nsize) / 350, 0.78)
assert round(sum(rec[("open", s)][0] / rec[("open", s)][1] * nsize[s] / N for s in nsize), 4) == 0.8325
assert round(sum(rec[("pcnl", s)][0] / rec[("pcnl", s)][1] * nsize[s] / N for s in nsize), 4) == 0.7789

# %% [markdown]
# ## Part 3: potential outcomes, Neyman's variance, blocking (Example 4.2, eqs. 3.2, 3.6, 3.7)

# %%
Y0 = np.array([2, 3, 4, 5, 6, 7, 8, 9.]); Y1 = np.array([5, 6, 6, 7, 7, 8, 10, 11.])
A8 = np.array([1, 1, 1, 1, 0, 0, 0, 0])
ate = (Y1 - Y0).mean(); att = (Y1 - Y0)[A8 == 1].mean()
assoc = Y1[A8 == 1].mean() - Y0[A8 == 0].mean()
sel = Y0[A8 == 1].mean() - Y0[A8 == 0].mean()
print(f"Ex 4.2: ATE {ate}, ATT {att}, association {assoc}, selection bias {sel}")
assert (ate, att, assoc, sel) == (2.0, 2.5, -1.5, -4.0) and np.isclose(assoc, att + sel)


def neyman_exact(Y0, Y1, n1, blocks=None):
    """Exact randomization distribution of the difference in means over all assignments
    (complete randomization, or within-block randomization with n1 per block if blocks given)."""
    N = len(Y0)
    if blocks is None:
        assigns = [np.isin(np.arange(N), c) for c in itertools.combinations(range(N), n1)]
    else:
        per = [list(itertools.combinations(np.flatnonzero(blocks == b), n1)) for b in np.unique(blocks)]
        assigns = [np.isin(np.arange(N), sum(c, ())) for c in itertools.product(*per)]
    est = np.array([Y1[w].mean() - Y0[~w].mean() for w in assigns])
    return est


def neyman_var(Y0, Y1, n1):
    """Var = S1^2/n1 + S0^2/n0 - S_tau^2/N, eq. (3.6)."""
    N = len(Y0); n0 = N - n1
    return Y1.var(ddof=1) / n1 + Y0.var(ddof=1) / n0 - (Y1 - Y0).var(ddof=1) / N


est = neyman_exact(Y0, Y1, 4)
assert len(est) == 70 and np.isclose(est.mean(), ate)          # unbiased over the design
assert np.isclose(est.var(), neyman_var(Y0, Y1, 4)) and np.isclose(est.var(), 2.5)
assert np.isclose(Y1.var(ddof=1), 30 / 7) and np.isclose(Y0.var(ddof=1), 6)
blk = np.array([0, 0, 0, 0, 1, 1, 1, 1])
estb = neyman_exact(Y0, Y1, 2, blocks=blk)
assert len(estb) == 36 and np.isclose(estb.mean(), ate)
print(f"Neyman: exact Var {est.var():.4f} = formula {neyman_var(Y0, Y1, 4):.4f}; S_tau^2 = {(Y1-Y0).var(ddof=1):.4f}; "
      f"S1^2 = {Y1.var(ddof=1):.4f}; conservative bound {Y1.var(ddof=1)/4 + Y0.var(ddof=1)/4:.4f}; blocked Var {estb.var():.4f}")
assert np.isclose(estb.var(), 0.875) and estb.var() < est.var() / 2

# %% [markdown]
# ## Part 4: Fisher's randomization test (Example 4.3) and a small experiment vs the t test

# %%
def fisher_exact_test(y, a):
    """Sharp null Y^1 = Y^0: re-randomize labels over all C(N, n1) assignments (eq. 3.5)."""
    y = np.asarray(y, float); a = np.asarray(a, bool); N, n1 = len(y), a.sum()
    obs = y[a].mean() - y[~a].mean()
    null = np.array([y[list(c)].mean() - np.delete(y, list(c)).mean()
                     for c in itertools.combinations(range(N), n1)])
    return obs, null, np.mean(null >= obs - 1e-12), np.mean(np.abs(null) >= abs(obs) - 1e-12)


obs, null, p1, p2 = fisher_exact_test([7, 9, 5, 4, 6, 3], [1, 1, 1, 0, 0, 0])
print(f"Ex 4.3: diff {obs:.4f}, {len(null)} assignments, one-sided p {p1:.2f}, two-sided p {p2:.2f}")
assert len(null) == 20 and np.isclose(p1, 0.10) and np.isclose(p2, 0.20)
ysm = np.r_[rng.normal(1.0, 1, 5), rng.normal(0, 1, 5)]
_, _, _, pf = fisher_exact_test(ysm, [1] * 5 + [0] * 5)
pt = stats.ttest_ind(ysm[:5], ysm[5:]).pvalue
print(f"10-unit experiment: randomization p = {pf:.4f} (252 assignments), t-test p = {pt:.4f}")

# %% [markdown]
# ## Part 5: standardization and IPW (Example 4.4) and in simulation
# g-formula (3.10) and IPW (3.11); both recover the ATE when L is measured, both fail when it is not.

# %%
cnt = {(0, 0): (90, 450), (0, 1): (45, 150), (1, 0): (50, 100), (1, 1): (210, 300)}   # (L, A): (Y=1, n)
n_tot = sum(v[1] for v in cnt.values())
PL = {l: sum(cnt[(l, a)][1] for a in (0, 1)) / n_tot for l in (0, 1)}
pA = {l: cnt[(l, 1)][1] / (cnt[(l, 0)][1] + cnt[(l, 1)][1]) for l in (0, 1)}
g = {a: sum(cnt[(l, a)][0] / cnt[(l, a)][1] * PL[l] for l in (0, 1)) for a in (0, 1)}
ipw = {a: sum(cnt[(l, a)][0] / (pA[l] if a else 1 - pA[l]) for l in (0, 1)) / n_tot for a in (0, 1)}
crude4 = {a: sum(cnt[(l, a)][0] for l in (0, 1)) / sum(cnt[(l, a)][1] for l in (0, 1)) for a in (0, 1)}
print(f"Ex 4.4: g-formula {g[1]:.2f}/{g[0]:.2f}, IPW {ipw[1]:.2f}/{ipw[0]:.2f}, crude {crude4[1]:.4f}/{crude4[0]:.4f}")
assert np.allclose([g[1], g[0], ipw[1], ipw[0]], [0.46, 0.32, 0.46, 0.32])
assert round(crude4[1] - crude4[0], 4) == 0.3121


def sim_binary(n, randomize=False):
    """L ~ N(0,1); A ~ Bern(expit(1.5 L)) or Bern(1/2); Y = 1 + 2A + 1.5L + 0.5 A L + noise. ATE = 2."""
    L = rng.standard_normal(n)
    e = np.full(n, 0.5) if randomize else 1 / (1 + np.exp(-1.5 * L))
    A = (rng.random(n) < e).astype(float)
    Y = 1 + 2 * A + 1.5 * L + 0.5 * A * L + rng.standard_normal(n)
    return L, A, Y, e


L, A, Y, e = sim_binary(100_000)
assoc_c = Y[A == 1].mean() - Y[A == 0].mean()
Lr, Ar, Yr, _ = sim_binary(100_000, randomize=True)
assoc_r = Yr[Ar == 1].mean() - Yr[Ar == 0].mean()
assert abs(assoc_c - 2) > 1 and abs(assoc_r - 2) < 0.05            # project item (a)
b1, _ = ols(Y[A == 1], L[A == 1]); b0, _ = ols(Y[A == 0], L[A == 0])
std_est = np.mean(b1[0] + b1[1] * L) - np.mean(b0[0] + b0[1] * L)    # g-formula with outcome regressions
# Fit the propensity by Newton's method for logistic regression (ST-4's IRLS).
th = np.zeros(2); X = np.column_stack([np.ones_like(L), L])
for _ in range(25):
    pr = 1 / (1 + np.exp(-X @ th)); th += np.linalg.solve(X.T @ (X * (pr * (1 - pr))[:, None]), X.T @ (A - pr))
ehat = 1 / (1 + np.exp(-X @ th))
ipw_est = np.mean(A * Y / ehat) - np.mean((1 - A) * Y / (1 - ehat))
omit = Y[A == 1].mean() - Y[A == 0].mean()                         # "standardize" over nothing
print(f"Sim: crude {assoc_c:.3f}, randomized {assoc_r:.3f}, g-formula {std_est:.3f}, IPW {ipw_est:.3f}, propensity {th.round(3)}")
assert abs(std_est - 2) < 0.05 and abs(ipw_est - 2) < 0.1 and abs(omit - 2) > 1   # project item (b)

# %% [markdown]
# ## Part 6: the collider, the mediator, M-bias, and the backdoor criterion (Examples 4.5, 4.6)

# %%
# (c) collider: two independent causes become correlated within a selected sample.
talent, looks = rng.standard_normal((2, 50_000))
famous = talent + looks > 1.5
r_all = np.corrcoef(talent, looks)[0, 1]; r_sel = np.corrcoef(talent[famous], looks[famous])[0, 1]
print(f"Berkson: corr all {r_all:.3f}, corr among selected {r_sel:.3f}")
assert abs(r_all) < 0.02 and r_sel < -0.4
fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(talent[~famous][:3000], looks[~famous][:3000], s=2, c="0.7", label="not selected")
ax.scatter(talent[famous][:1500], looks[famous][:1500], s=2, c="C3", label="selected (collider = 1)")
ax.set(xlabel="talent", ylabel="looks", title="Conditioning on a collider creates correlation"); ax.legend()
fig.savefig(FIG / "berkson.png", dpi=120, bbox_inches="tight"); plt.close(fig)

# (d) mediator: A -> M -> Y plus A -> Y; adjusting M removes the indirect effect a*b (eq. 3.16).
med = SCM(["A", "M", "Y"], {("A", "M"): 0.8, ("M", "Y"): 0.5, ("A", "Y"): 0.3})
assert np.isclose(med.pop_coef("Y", ["A"])[0], 0.7) and np.isclose(med.pop_coef("Y", ["A", "M"])[0], 0.3)
# Mediator-outcome confounder U: adjusting M gives c - a g h / (1 + h^2)  (eq. 3.17).
medU = SCM(["U", "A", "M", "Y"], {("A", "M"): 0.8, ("M", "Y"): 0.5, ("A", "Y"): 0.3, ("U", "M"): 1.0, ("U", "Y"): 1.0})
cdir = medU.pop_coef("Y", ["A", "M"])[0]
assert np.isclose(cdir, 0.3 - 0.8 * 1 * 1 / 2)
print(f"Mediator: total 0.7, adjusted-for-M 0.3; with U->M,U->Y adjusted coefficient {cdir:.2f}")

# (e) M-bias: A <- U1 -> C <- U2 -> Y; C is pre-treatment but adjusting it is wrong.
mb = SCM(["U1", "U2", "C", "A", "Y"], {("U1", "C"): 1, ("U2", "C"): 1, ("U1", "A"): 1, ("U2", "Y"): 1, ("A", "Y"): 0.5})
mb0, mbC = mb.pop_coef("Y", ["A"])[0], mb.pop_coef("Y", ["A", "C"])[0]
print(f"M-bias: Y on A = {mb0:.4f}; Y on (A, C) = {mbC:.4f}")
assert np.isclose(mb0, 0.5) and not np.isclose(mbC, 0.5)
assert mb.backdoor_ok("A", "Y", []) and not mb.backdoor_ok("A", "Y", ["C"])

# Example 4.6: the five-node DAG. Enumerate all adjustment sets.
g5 = SCM(["Z1", "Z2", "Z3", "A", "Y"], {("Z1", "Z3"): 0.8, ("Z2", "Z3"): 0.8, ("Z1", "A"): 0.7,
                                         ("Z3", "A"): 0.6, ("Z3", "Y"): 0.5, ("Z2", "Y"): 0.7, ("A", "Y"): 1.0})
print("Ex 4.6 backdoor paths:", [p for p in g5.skeleton_paths("A", "Y") if p[1] in g5.parents("A")])
valid = []
for r in range(4):
    for Zs in itertools.combinations(["Z1", "Z2", "Z3"], r):
        ok = g5.backdoor_ok("A", "Y", Zs); coef = g5.pop_coef("Y", ["A", *Zs])[0]
        print(f"  adjust {set(Zs) or '{}'}: backdoor {'valid  ' if ok else 'invalid'}  coefficient {coef:.4f}")
        assert ok == np.isclose(coef, 1.0)
        valid += [set(Zs)] * ok
assert valid == [{"Z1", "Z3"}, {"Z2", "Z3"}, {"Z1", "Z2", "Z3"}]
d5 = g5.sample(100_000); bb, _ = ols(d5["Y"], d5["A"], d5["Z1"], d5["Z3"])
assert abs(bb[1] - 1) < 0.02

# %% [markdown]
# ## Part 7: factorial designs (Example 4.8)

# %%
def factorial_effects(y, k):
    """Main effects and interactions of a 2^k design in standard order: contrast/2^(k-1) (eq. 3.9)."""
    X = np.array(list(itertools.product([-1, 1], repeat=k)))[:, ::-1]
    out = {}
    for r in range(1, k + 1):
        for S in itertools.combinations(range(k), r):
            out["".join("ABC"[i] for i in S)] = np.prod(X[:, S], axis=1) @ y / 2 ** (k - 1)
    return out, X


eff2, _ = factorial_effects(np.array([60, 72, 54, 80.]), 2)
print("Ex 4.8 2^2 effects:", eff2)
assert eff2 == {"A": 19.0, "B": 1.0, "AB": 7.0}
truth = {"A": 4.0, "B": -2.0, "C": 0.0, "AB": 3.0, "AC": 0.0, "BC": 0.0, "ABC": 0.0}
_, X3 = factorial_effects(np.zeros(8), 3)
mean3 = 50 + 0.5 * (truth["A"] * X3[:, 0] + truth["B"] * X3[:, 1] + truth["AB"] * X3[:, 0] * X3[:, 1])
reps = mean3 + rng.normal(0, 1, (4, 8))
eff3, _ = factorial_effects(reps.mean(0), 3)
assert all(abs(eff3[k] - truth[k]) < 1.2 for k in truth)
fig, axs = plt.subplots(1, 2, figsize=(9, 3.5))
axs[0].plot([-1, 1], [60, 72], "o-", label="B low"); axs[0].plot([-1, 1], [54, 80], "s-", label="B high")
axs[0].set(xlabel="factor A", ylabel="mean response", title="Example 4.8: non-parallel = interaction"); axs[0].legend()
for b, m in ((-1, "o-"), (1, "s-")):
    axs[1].plot([-1, 1], [reps.mean(0)[(X3[:, 0] == a) & (X3[:, 1] == b)].mean() for a in (-1, 1)], m, label=f"B={b}")
axs[1].set(xlabel="factor A", title="Simulated 2^3, 4 replicates"); axs[1].legend()
fig.savefig(FIG / "interaction.png", dpi=120, bbox_inches="tight"); plt.close(fig)

# %% [markdown]
# ## Part 8: instrumental variables (Example 4.7) and weak instruments

# %%
wald = (5.00 - 4.10) / (0.70 - 0.20)
print(f"Ex 4.7: Wald = {wald:.2f}, complier share {0.70 - 0.20:.2f}")
assert np.isclose(wald, 1.8)


def tsls(y, a, Z):
    """2SLS: regress a on [1, Z], then y on [1, a_hat] (eq. 3.23); returns (beta, first-stage F)."""
    Zc = np.column_stack([np.ones(len(a)), Z])
    g, *_ = np.linalg.lstsq(Zc, a, rcond=None); ahat = Zc @ g
    b, _ = ols(y, ahat)
    r_full = a - ahat; r_null = a - a.mean(); k = Zc.shape[1] - 1
    F = ((r_null @ r_null - r_full @ r_full) / k) / (r_full @ r_full / (len(a) - k - 1))
    return b[1], F


def iv_sim(n, pi, beta=1.0, reps=1):
    U = rng.standard_normal((reps, n)); Z = rng.standard_normal((reps, n))
    Aiv = pi * Z + U + rng.standard_normal((reps, n))
    Yiv = beta * Aiv + 2 * U + rng.standard_normal((reps, n))
    return Z, Aiv, Yiv


Z, Aiv, Yiv = iv_sim(100_000, 0.5)
b_ols = ols(Yiv[0], Aiv[0])[0][1]; b_iv, F_iv = tsls(Yiv[0], Aiv[0], Z[0])
print(f"IV sim: OLS {b_ols:.3f} (plim {1 + 2/2.25:.3f}), 2SLS {b_iv:.3f}, first-stage F {F_iv:.0f}")
assert abs(b_ols - (1 + 2 / 2.25)) < 0.02 and abs(b_iv - 1) < 0.05
pis = np.array([0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0]); med_bias, meanF = [], []
for pi in pis:
    Zs, As, Ys = iv_sim(200, pi, reps=1000)
    res = [tsls(Ys[r], As[r], Zs[r]) for r in range(1000)]
    med_bias.append(np.median([b for b, _ in res]) - 1); meanF.append(np.mean([F for _, F in res]))
med_bias, meanF = np.array(med_bias), np.array(meanF)
ols_bias = 2 / (pis ** 2 + 2)
print("Weak IV: mean F", meanF.round(1), "\n  median 2SLS bias", med_bias.round(3), "\n  OLS bias", ols_bias.round(3))
assert med_bias[0] > 0.5 * ols_bias[0] and abs(med_bias[-1]) < 0.05
fig, ax = plt.subplots(figsize=(5.5, 4))
ax.semilogx(meanF, med_bias, "o-", label="2SLS median bias"); ax.semilogx(meanF, ols_bias, "s--", label="OLS bias")
ax.axvline(10, c="0.5", ls=":", label="F = 10")
ax.set(xlabel="mean first-stage F", ylabel="bias", title="Weak instruments: 2SLS drifts to OLS"); ax.legend()
fig.savefig(FIG / "weak_iv.png", dpi=120, bbox_inches="tight"); plt.close(fig)

# %% [markdown]
# ## Part 9: regression discontinuity and difference-in-differences (§3.9)

# %%
R = rng.uniform(-1, 1, 4000); D = (R >= 0).astype(float)
Yrd = 1 + 0.8 * R + 1.5 * R ** 2 + 2.0 * D + 0.5 * rng.standard_normal(4000)


def rd_local_linear(R, Y, h):
    """Local linear fits with triangular weights on each side; jump at 0."""
    ends = []
    for side in (R < 0, R >= 0):
        m = side & (np.abs(R) < h); w = np.sqrt(1 - np.abs(R[m]) / h)
        X = np.column_stack([np.ones(m.sum()), R[m]]) * w[:, None]
        ends.append(np.linalg.lstsq(X, Y[m] * w, rcond=None)[0][0])
    return ends[1] - ends[0]


hs = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
rds = [rd_local_linear(R, Yrd, h) for h in hs]
naive = Yrd[D == 1].mean() - Yrd[D == 0].mean()
print("RD estimates by bandwidth:", dict(zip(hs, np.round(rds, 3))), f"naive {naive:.3f}")
assert abs(rds[3] - 2) < 0.15 and abs(naive - 2) > 0.5
fig, ax = plt.subplots(figsize=(5.5, 4))
ax.scatter(R[::4], Yrd[::4], s=3, c="0.6")
for side in (R < 0, R >= 0):
    m = side & (np.abs(R) < 0.3); bfit, _ = ols(Yrd[m], R[m]); xs = np.sort(R[m]); ax.plot(xs, bfit[0] + bfit[1] * xs, "C3")
ax.axvline(0, c="k", lw=0.8); ax.set(xlabel="running variable", ylabel="Y", title="RD: local linear fits, h = 0.3")
fig.savefig(FIG / "rd.png", dpi=120, bbox_inches="tight"); plt.close(fig)


def did(parallel, n=20_000, tau=2.0):
    G = rng.integers(0, 2, n); T = rng.integers(0, 2, n)
    trend = 1.0 + (0 if parallel else 0.7) * G
    y = 3 + 1.5 * G + trend * T + tau * G * T + rng.standard_normal(n)
    b, _ = ols(y, G, T, G * T)
    return b[3]


did_p, did_np = did(True), did(False)
print(f"DiD: parallel trends {did_p:.3f}; diverging trends {did_np:.3f} (bias {did_np - 2:.2f})")
assert abs(did_p - 2) < 0.1 and abs(did_np - 2.7) < 0.1

# %% [markdown]
# ## Part 10: sensitivity analysis (eq. 3.24)
# An unmeasured U with effect gamma on Y and imbalance delta (coefficient of U on A given L)
# shifts the adjusted coefficient by gamma*delta.

# %%
beta_hat = 0.30
gam = np.linspace(0.01, 1.5, 200)
fig, ax = plt.subplots(figsize=(5.5, 4))
for dl in (0.1, 0.2, 0.4, 0.6):
    ax.plot(gam, beta_hat - gam * dl, label=f"delta = {dl}")
ax.axhline(0, c="k", lw=0.8); ax.set(xlabel="gamma (effect of U on Y)", ylabel="bias-corrected effect",
                                     title="How strong must U be to explain 0.30 away?"); ax.legend()
fig.savefig(FIG / "sensitivity.png", dpi=120, bbox_inches="tight"); plt.close(fig)
sens = SCM(["U", "A", "Y"], {("U", "A"): 0.6, ("U", "Y"): 0.5, ("A", "Y"): 0.0}, omega={"A": 0.64})
assert np.isclose(sens.pop_coef("Y", ["A"])[0], 0.5 * 0.6 / (0.36 + 0.64))   # gamma*delta with delta=Cov(A,U)/Var A
evalue = lambda rr: rr + np.sqrt(rr * (rr - 1))
print(f"E-value(RR=2) = {evalue(2):.3f}")

# %% [markdown]
# ## Part 11: Granger "causality" from a common cause, and prediction vs intervention (§3.11)

# %%
n = 3000; Zt = np.zeros(n)
for t in range(1, n):
    Zt[t] = 0.8 * Zt[t - 1] + rng.standard_normal()
Xt = np.r_[0, Zt[:-1]] + 0.5 * rng.standard_normal(n)        # X_t = Z_{t-1} + noise
Yt = np.r_[0, 0, Zt[:-2]] + 0.5 * rng.standard_normal(n)     # Y_t = Z_{t-2} + noise: no X -> Y arrow


def granger_F(y, x, p=2):
    """F test of y_t on its own p lags vs adding p lags of x (ST-7)."""
    T = len(y); Ylag = np.column_stack([y[p - k - 1:T - k - 1] for k in range(p)])
    Xlag = np.column_stack([x[p - k - 1:T - k - 1] for k in range(p)]); yy = y[p:]
    _, r0 = ols(yy, *Ylag.T); _, r1 = ols(yy, *Ylag.T, *Xlag.T)
    F = ((r0 @ r0 - r1 @ r1) / p) / (r1 @ r1 / (len(yy) - 2 * p - 1))
    return F, stats.f.sf(F, p, len(yy) - 2 * p - 1)


Fg, pg = granger_F(Yt, Xt)
print(f"Granger X -> Y: F = {Fg:.1f}, p = {pg:.2e}; intervening on X changes Y by exactly 0 (no arrow).")
assert pg < 1e-10

S_ = rng.standard_normal(50_000)                               # severity
A_ = (rng.random(50_000) < 1 / (1 + np.exp(-3 * S_))).astype(float)   # sicker patients get treated
Ybad = 2 * S_ - 1.0 * A_ + 0.5 * rng.standard_normal(50_000)  # treatment lowers the bad outcome by 1
proxy = S_ + 0.6 * rng.standard_normal(50_000)                     # what the model sees instead of S
bp, rp = ols(Ybad, A_, proxy)
r2 = 1 - rp.var() / Ybad.var()
print(f"Predictive model: R^2 = {r2:.2f}, coefficient on treatment = {bp[1]:+.2f} (true causal effect -1.00)")
assert r2 > 0.6 and bp[1] > 0
pol = {a: np.mean(2 * S_ - 1.0 * a) for a in (0, 1)}
print(f"  model says 'never treat'; true mean outcome under do(A=0) {pol[0]:+.2f} vs do(A=1) {pol[1]:+.2f}")

# %% [markdown]
# ## Part 12: A/B test sizing and peeking (§3.3)

# %%
za, zb = stats.norm.ppf(0.975), stats.norm.ppf(0.8)
n_ab = (za + zb) ** 2 * (0.1 * 0.9 + 0.11 * 0.89) / 0.01 ** 2
print(f"A/B: n per arm to detect 10% -> 11% at 80% power: {np.ceil(n_ab):.0f}")
assert 14_700 < n_ab < 14_800
looks = np.cumsum(rng.standard_normal((4000, 10)), axis=1) / np.sqrt(np.arange(1, 11))
fp = np.mean((np.abs(looks) > za).any(axis=1))
print(f"Peeking: A/A test, 10 looks, stop at first |z| > 1.96: false-positive rate {fp:.3f}")
assert fp > 0.15

# %% [markdown]
# ## Part 13: numbers for exercises.md / solutions.md

# %%
y0 = np.array([4, 6, 5, 8, 7, 9.]); y1 = np.array([6, 7, 8, 8, 9, 10.]); a6 = np.array([1, 1, 1, 0, 0, 0]) == 1
print(f"A1: ATE {(y1-y0).mean()}, ATT {(y1-y0)[a6].mean()}, assoc {y1[a6].mean()-y0[~a6].mean()}, "
      f"SB {y0[a6].mean()-y0[~a6].mean()}")
gA2 = SCM(["L", "A", "M", "Y", "C"], {("L", "A"): 1, ("L", "Y"): 1, ("A", "M"): 1, ("M", "Y"): 1, ("A", "Y"): 1,
                                     ("A", "C"): 1, ("Y", "C"): 1})
okA2 = [set(Zs) for r in range(4) for Zs in itertools.combinations(["L", "M", "C"], r) if gA2.backdoor_ok("A", "Y", Zs)]
print("A2 valid sets:", okA2); assert okA2 == [{"L"}]
o3, n3, p3, _ = fisher_exact_test([5, 7, 2, 3], [1, 1, 0, 0]); print(f"A3: diff {o3}, one-sided p {p3:.4f}")
print(f"A4: Wald {(12.6 - 12.0) / (0.45 - 0.15):.2f}")
print("A5:", factorial_effects(np.array([20, 30, 24, 26.]), 2)[0])
print(f"A6: DiD {(15 - 10) - (11 - 8)}")
cB = {(0, 0): (48, 160), (0, 1): (20, 40), (1, 0): (30, 60), (1, 1): (168, 240)}
nB = 500; PLB = {l: (cB[(l, 0)][1] + cB[(l, 1)][1]) / nB for l in (0, 1)}
pAB = {l: cB[(l, 1)][1] / (cB[(l, 0)][1] + cB[(l, 1)][1]) for l in (0, 1)}
gB = {a: sum(cB[(l, a)][0] / cB[(l, a)][1] * PLB[l] for l in (0, 1)) for a in (0, 1)}
iB = {a: sum(cB[(l, a)][0] / (pAB[l] if a else 1 - pAB[l]) for l in (0, 1)) / nB for a in (0, 1)}
crB = {a: sum(cB[(l, a)][0] for l in (0, 1)) / sum(cB[(l, a)][1] for l in (0, 1)) for a in (0, 1)}
print(f"B1: g {gB}, IPW {iB}, crude {crB[1]:.4f} - {crB[0]:.4f} = {crB[1]-crB[0]:.4f}")
assert np.allclose([gB[1], gB[0]], [0.62, 0.42]) and np.allclose([iB[1], iB[0]], [0.62, 0.42])
pseudo = {(l, a): cB[(l, a)][1] / (pAB[l] if a else 1 - pAB[l]) for l in (0, 1) for a in (0, 1)}
pseudoY = {(l, a): cB[(l, a)][0] / (pAB[l] if a else 1 - pAB[l]) for l in (0, 1) for a in (0, 1)}
print("B1(d) pseudo-population n:", {k: round(v) for k, v in pseudo.items()}, " Y=1:", {k: round(v) for k, v in pseudoY.items()})
assert np.allclose([pseudo[(1, a)] / (pseudo[(0, a)] + pseudo[(1, a)]) for a in (0, 1)], 0.6)
eB2 = neyman_exact(y0, y1, 3)
print(f"B2: 20 assignments, exact Var {eB2.var():.4f}, formula {neyman_var(y0, y1, 3):.4f}, "
      f"S1^2 {y1.var(ddof=1):.4f}, S0^2 {y0.var(ddof=1):.4f}, Stau^2 {(y1-y0).var(ddof=1):.4f}, "
      f"conservative {y1.var(ddof=1)/3 + y0.var(ddof=1)/3:.4f}")
assert np.isclose(eB2.var(), neyman_var(y0, y1, 3))
fB3 = SCM(["L", "A", "Y"], {("L", "A"): 1, ("L", "Y"): 1, ("A", "Y"): 0.5})
fcB3 = SCM(["L", "A", "Y", "C"], {("L", "A"): 1, ("L", "Y"): 1, ("A", "Y"): 0.5, ("A", "C"): 1, ("Y", "C"): 1})
print(f"B3: Cov(A,Y) {fB3.trek_cov('A','Y'):.2f}, Y on A {fB3.pop_coef('Y',['A'])[0]:.3f}, "
      f"Y on (A,L) {fB3.pop_coef('Y',['A','L'])[0]:.3f}, Y on (A,L,C) {fcB3.pop_coef('Y',['A','L','C'])[0]:.4f}")
pY = 0.2; pAsel = 0.3; s = lambda a, y: 0.1 + 0.4 * a + 0.4 * y
q = {a: pY * s(a, 1) / (pY * s(a, 1) + (1 - pY) * s(a, 0)) for a in (0, 1)}
print(f"B4: P(Y=1|A=1,S=1) {q[1]:.4f}, P(Y=1|A=0,S=1) {q[0]:.4f}, ratio {q[1]/q[0]:.4f}")
g6 = SCM(["L1", "A", "L2", "M", "Y", "C"], {("L1", "A"): 1.0, ("L1", "L2"): 1.0, ("L2", "Y"): 1.0, ("A", "M"): 1.0,
                                             ("M", "Y"): 0.8, ("A", "Y"): 0.5, ("A", "C"): 1.0, ("Y", "C"): 1.0})
S6 = g6.cov()
for Zs in ([], ["L1"], ["L2"], ["L1", "L2"], ["M"], ["L1", "C"]):
    xs = ["A", *Zs]; b = g6.pop_coef("Y", xs); ia = [g6.ix[x] for x in xs]
    resid = S6[g6.ix["Y"], g6.ix["Y"]] - S6[g6.ix["Y"], ia] @ b
    vA = 1 / np.linalg.inv(S6[np.ix_(ia, ia)])[0, 0]                 # Var(A | Zs)
    print(f"B5: adjust {Zs}: backdoor {g6.backdoor_ok('A','Y',Zs)}, coef {b[0]:.4f}, n*Var(coef) {resid / vA:.4f}")
F6 = []
for r in range(1000):
    Zr, Ar, Yr = iv_sim(400, 0.5); F6.append(tsls(Yr[0], Ar[0], Zr[0])[1])
print(f"B6: OLS plim {1 + 2/2.25:.4f}, concentration n*pi^2/Var(A|Z) = {400 * 0.25 / 2:.0f}, simulated mean F {np.mean(F6):.1f}")
assert 45 < np.mean(F6) < 57
print(f"B7: delta needed at gamma=0.5: {0.30/0.5:.2f}; E-value(1.5) = {evalue(1.5):.3f}")
pc, tc, pd_, td = 0.4, 2.0, 0.1, 5.0
print(f"C1: Wald with defiers = {(pc*tc - pd_*td)/(pc - pd_):.2f}")
for k, pi in ((1, 0.15), (5, 0.15)):
    bs, Fs = [], []
    for r in range(1000):
        Zk = rng.standard_normal((200, k)); U = rng.standard_normal(200)
        Ak = Zk @ np.full(k, pi) + U + rng.standard_normal(200); Yk = Ak + 2 * U + rng.standard_normal(200)
        b, F = tsls(Yk, Ak, Zk); bs.append(b); Fs.append(F)
    print(f"C3: k={k}, pi={pi}: mean F {np.mean(Fs):.2f}, median bias {np.median(bs)-1:.3f}, "
          f"mean bias {np.mean(bs)-1:.3f}, OLS bias {2/(k*pi**2+2):.3f}")

# %% [markdown]
# ## Your turn
# 1. Propensity-score matching: for each treated unit in Part 5, find the control with the
#    nearest ehat and average the differences. That estimates the ATT. Compute the true ATT
#    from the simulation equation (it is 2 + 0.5 E[L | A=1]) and compare. Then match on L itself.
# 2. Draw a DAG of your own with at least six nodes and three distinct valid adjustment sets.
#    Build it with SCM, confirm with backdoor_ok, and check by simulation that all three give
#    the same coefficient while their standard errors differ (exercise B5 is a template).
# 3. Monty Hall (F4): nodes Car, Pick, Host with Car -> Host <- Pick. Simulate 100,000 games and
#    check Car and Pick are independent but dependent given Host. Explain switching with the collider.
# 4. Replace the logistic propensity in sim_binary by expit(4 L). Watch IPW's variance explode
#    as ehat approaches 0 and 1 (positivity, §3.4). Try trimming weights at the 1st/99th percentiles.

# %%
if __name__ == "__main__":
    print(f"All checks passed in {time.time() - T0:.1f} s. Figures in {FIG}")
