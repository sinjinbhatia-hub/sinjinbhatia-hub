# PY-8 Solutions

Every block below was run, in order, and its output pasted from the run; each solution starts from a fresh session, and its asserts passed.

## A. Warm-up

**A1.**

```python
import sympy as sp

x, n = sp.Symbol("x", positive=True), sp.Symbol("n", integer=True)
print(sp.sqrt(18), sp.sqrt(18.0), sp.Rational(3, 4) * 0.5, sp.S(3) / 4 * 2)
print(sp.integrate(1 / x, x), sp.diff(sp.Abs(x), x), sp.solve(x**2 + 1, x), sp.sin(n * sp.pi))


def kinetic_bad(m, v):
    return 1 / 2 * m * v**2


def kinetic(m, v):
    return sp.Rational(1, 2) * m * v**2


def kinetic2(m, v):
    return sp.S(1) / 2 * m * v**2


m_, v_ = sp.symbols("m v")
print(kinetic_bad(m_, v_), "|", kinetic(m_, v_), "|", kinetic2(m_, v_))
assert kinetic_bad(m_, v_).has(sp.Float) and not kinetic(m_, v_).has(sp.Float) and kinetic(m_, v_) == kinetic2(m_, v_)
```
```text
3*sqrt(2) 4.24264068711928 0.375000000000000 3/2
log(x) 1 [] 0
0.5*m*v**2 | m*v**2/2 | m*v**2/2
```

`sqrt(18)` stays exact as $3\sqrt2$; one float (`18.0`, `0.5`) makes the result a float; `S(3) / 4 * 2` is the exact $3/2$. For a positive `x`, $|x|=x$ so the derivative is 1, and $x^2+1=0$ has no positive root, so `solve` returns `[]`; $\sin(n\pi)=0$ because `n` is an integer. In `kinetic_bad`, Python evaluates `1 / 2` to the float `0.5` before SymPy sees it. The fixes are `sp.Rational(1, 2)` or `sp.S(1) / 2`.

**A2.**

```python
import numpy as np
from scipy import integrate, optimize


def require(res):
    """Return a SciPy result object, or raise if it reports failure."""
    ok = getattr(res, "success", getattr(res, "converged", None))
    if not ok:
        raise RuntimeError(getattr(res, "message", getattr(res, "flag", "no message")))
    return res


rosen = lambda p: (1 - p[0]) ** 2 + 100 * (p[1] - p[0] ** 2) ** 2
f = lambda x: np.cos(x) - x
cases = {
    "(i) minimize, maxiter=3": lambda: optimize.minimize(rosen, [-1.2, 1.0], options={"maxiter": 3}),
    "(ii) solve_ivp, blow-up": lambda: integrate.solve_ivp(lambda t, y: [y[0] ** 2], (0, 2), [1.0]),
    "(iii) root_scalar, maxiter=3": lambda: optimize.root_scalar(f, bracket=[0, 1], method="brentq", maxiter=3),
    "(iv) root_scalar, no sign change": lambda: optimize.root_scalar(f, bracket=[1, 2], method="brentq"),
}
for name, run in cases.items():
    try:
        res = run()
    except ValueError as exc:
        print(f"{name}: the routine itself raised ValueError: {exc}")
        continue
    extra = f"status={res.status}, t_end={res.t[-1]:.5f}" if hasattr(res, "t") else f"status={getattr(res, 'status', '-')}"
    try:
        require(res)
    except RuntimeError as exc:
        print(f"{name}: returned a failure flag ({extra}); require raised: {exc}")
```
```text
(i) minimize, maxiter=3: returned a failure flag (status=1); require raised: Maximum number of iterations has been exceeded.
(ii) solve_ivp, blow-up: returned a failure flag (status=-1, t_end=0.99993); require raised: Required step size is less than spacing between numbers.
(iii) root_scalar, maxiter=3: returned a failure flag (status=-); require raised: convergence error
(iv) root_scalar, no sign change: the routine itself raised ValueError: f(a) and f(b) must have different signs
```

Only (iv) raises by itself, because the bracket is invalid before any work is done. (i) returns `success=False` with `status=1` and "Maximum number of iterations has been exceeded"; (ii) returns `success=False` with `status=-1` because the step size collapses at the blow-up $t=1$; (iii) returns `converged=False` with the flag `convergence error`. Three of four failures would have passed silently through a script that never looks, and which routines raise is not something to memorize: call `require` on every result.

**A3.**

```python
import sympy as sp

t = sp.Symbol("t", real=True)
q = sp.Function("q")(t)
m, k = sp.symbols("m k", positive=True)


def to_state(expr, funcs, t):
    """Replace each q_i(t) and its derivatives by plain symbols, highest derivative first."""
    derivs = sorted(expr.atoms(sp.Derivative), key=lambda d: -d.derivative_count)
    syms = {d: sp.Symbol(f"{d.expr.func.__name__}_d{d.derivative_count}", real=True) for d in derivs}
    syms.update({f: sp.Symbol(f.func.__name__, real=True) for f in funcs})
    out = expr
    for d in derivs:
        out = out.subs(d, syms[d])
    for f in funcs:
        out = out.subs(f, syms[f])
    return out, syms


Lag = sp.Rational(1, 2) * m * q.diff(t) ** 2 - sp.Rational(1, 2) * k * q**2
E = q.diff(t) * Lag.diff(q.diff(t)) - Lag                                  # (3.4)
good, syms = to_state(E, [q], t)
x = sp.Symbol("x")
bad = E.subs(q, x)
print("E:", E, "| right order:", good, "| wrong order:", bad, "->", bad.doit())
assert sp.simplify(good - (m * syms[q.diff(t)] ** 2 / 2 + k * syms[q] ** 2 / 2)) == 0
assert sp.simplify(bad.doit() - k * x**2 / 2) == 0                           # the kinetic term is gone
print(to_state(m * q.diff(t, 2) + 3 * q.diff(t) + k * q, [q], t)[0])
```
```text
E: k*q(t)**2/2 + m*Derivative(q(t), t)**2/2 | right order: k*q**2/2 + m*q_d1**2/2 | wrong order: k*x**2/2 + m*Derivative(x, t)**2/2 -> k*x**2/2
k*q + m*q_d2 + 3*q_d1
```

Sorting by `derivative_count` descending replaces $\ddot q$ before $\dot q$ before $q$, whatever the expression contains. Substituting $q\to x$ first turns $\dot q$ into `Derivative(x, t)`, the derivative of a symbol, which `doit()` evaluates to $0$: the energy loses its kinetic term and becomes $kx^2/2$ without any error message.

**A4.**

```python
import numpy as np
from scipy import integrate

c0, exact = 1e-8, 1e-8 * np.exp(-5)
for label, kw in [("default atol (1e-6)", {}), ("atol = rtol * 1e-10 = 1e-16", dict(atol=1e-16))]:
    s = integrate.solve_ivp(lambda t, c: -c, (0, 5), [c0], rtol=1e-6, **kw)
    rel = abs(s.y[0, -1] / exact - 1)
    print(f"{label:30s} c(5) = {s.y[0, -1]:.3e} (exact {exact:.3e}), relative error {rel:.1e}, nfev {s.nfev}")
    assert (rel > 1) == (not kw) and (rel < 1e-5) == bool(kw)
```
```text
default atol (1e-6)            c(5) = 1.213e-08 (exact 6.738e-11), relative error 1.8e+02, nfev 14
atol = rtol * 1e-10 = 1e-16    c(5) = 6.738e-11 (exact 6.738e-11), relative error 1.6e-06, nfev 134
```

With the default `atol=1e-6`, the scale in (3.8) is $\text{atol}+\text{rtol}\,|c|\approx10^{-6}$, a hundred times the initial value, so the solver accepts any step whose error is below $10^{-6}$ *absolute*: after three steps it reports $c(5)\approx10^{-8}$ against the true $6.7\times10^{-11}$, an error of more than 100 times the answer, with `success=True`. The smallest value we care about is $c(5)\sim10^{-10}$, so the rule `atol` $=$ `rtol` $\times$ that gives $10^{-16}$, and then the scale is dominated by `rtol` $|c|$ as intended.

**A5.**

```python
import numpy as np
from scipy import sparse

D = sparse.diags([1, 2, 1], [-1, 0, 1], shape=(4, 4))
v = np.ones(4)


def apply_twice(A, v):
    return A * (A * v)


def apply_twice_ok(A, v):
    return A @ (A @ v)


for name, A in [("csr_matrix", sparse.csr_matrix(D)), ("csr_array", sparse.csr_array(D)), ("ndarray", D.toarray())]:
    out = apply_twice(A, v)
    shown = out.toarray().tolist() if sparse.issparse(out) else out.tolist()
    print(f"{name:10s} A*(A*v) is a {type(out).__name__}: {shown}")
    assert np.allclose(apply_twice_ok(A, v), [10, 15, 15, 10])
print("with @, all three give", apply_twice_ok(D.toarray(), v))
```
```text
csr_matrix A*(A*v) is a ndarray: [10.0, 15.0, 15.0, 10.0]
csr_array  A*(A*v) is a csr_array: [[4.0, 1.0, 0.0, 0.0], [1.0, 4.0, 1.0, 0.0], [0.0, 1.0, 4.0, 1.0], [0.0, 0.0, 1.0, 4.0]]
ndarray    A*(A*v) is a ndarray: [[4.0, 1.0, 0.0, 0.0], [1.0, 4.0, 1.0, 0.0], [0.0, 1.0, 4.0, 1.0], [0.0, 0.0, 1.0, 4.0]]
with @, all three give [10. 15. 15. 10.]
```

Only `csr_matrix` treats `*` as matrix multiplication and returns the vector $(10,15,15,10)$. For `csr_array` and for the dense array `*` is elementwise with broadcasting: `A * v` scales the columns of `A` by `v`, and the outer `*` multiplies elementwise again, giving the $4\times4$ matrix $A\circ A$ (sparse in one case, dense in the other), with no error. Writing `@` gives $(10,15,15,10)$ for all three; use it always.

**A6.**

```python
import numpy as np
from scipy import fft, integrate, special

right = integrate.dblquad(lambda y, x: x * y**2, 0, 2, 0, 1)[0]            # func(y, x): inner variable first
swapped = integrate.dblquad(lambda x, y: x * y**2, 0, 2, 0, 1)[0]
print("dblquad with func(y, x):", right, "| with func(x, y):", swapped)
print("ellipk(0.5) =", special.ellipk(0.5), "| ellipk(sin(pi/4)) =", special.ellipk(np.sin(np.pi / 4)))
print("fftfreq(8, d=0.01) =", fft.fftfreq(8, d=0.01))
assert abs(right - 2 / 3) < 1e-12 and abs(swapped - 4 / 3) < 1e-12
assert abs(special.ellipk(0.5) - 1.8540746773013719) < 1e-12 and np.allclose(fft.fftfreq(8, d=0.01), [0, 12.5, 25, 37.5, -50, -37.5, -25, -12.5])
```
```text
dblquad with func(y, x): 0.6666666666666667 | with func(x, y): 1.3333333333333335
ellipk(0.5) = 1.8540746773013719 | ellipk(sin(pi/4)) = 2.0859736981949366
fftfreq(8, d=0.01) = [  0.   12.5  25.   37.5 -50.  -37.5 -25.  -12.5]
```

(a) `dblquad(func, a, b, gfun, hfun)` integrates `func(y, x)` with $x\in[a,b]$ and $y\in[g(x),h(x)]$, so the lambda `lambda x, y: x * y**2` is really $\int_0^2\!\int_0^1 y\,x^2\,dy\,dx=\tfrac12\cdot\tfrac83=\tfrac43$, while `lambda y, x: x * y**2` gives $\int_0^2 x\,dx\int_0^1y^2\,dy=2\cdot\tfrac13=\tfrac23$. (b) The period factor is $K(m)$ with $m=\sin^2(\pi/4)=0.5$, so `ellipk(0.5)` is right; passing $k=\sin(\pi/4)$ treats it as a parameter and is wrong by 12.5%. (c) `fftfreq(8, d=0.01)` is $k/(Nd)=k\cdot12.5$ Hz for $k=0,\dots,3$, then $-4,\dots,-1$: the negative frequencies sit in the second half.


## B. Core

**B1.**

```python
import numpy as np
import sympy as sp
from scipy import integrate

t = sp.Symbol("t", real=True)
m, R, g, W = sp.symbols("m R g Omega", positive=True)
th = sp.Function("theta")(t)
Lag = sp.Rational(1, 2) * m * R**2 * (th.diff(t) ** 2 + W**2 * sp.sin(th) ** 2) + m * g * R * sp.cos(th)
EL = sp.diff(Lag.diff(th.diff(t)), t) - Lag.diff(th)                           # (3.3)
acc = sp.solve(EL, th.diff(t, 2))[0]
hand = sp.sin(th) * (W**2 * sp.cos(th) - g / R)
assert sp.simplify(acc - hand) == 0                                            # (a)
a, w = sp.symbols("a w", real=True)
state = lambda ex: ex.subs(th.diff(t), w).subs(th, a)                          # derivative first
acc_s = state(acc)
print("(a) theta'' =", sp.factor(acc_s))
print("(b) equilibria:", sp.solve(sp.Eq(acc_s, 0), a))
a_star = sp.acos(g / (R * W**2))
slope = sp.simplify(sp.diff(acc_s, a).subs(a, a_star))                          # d(theta'')/d(theta) at theta* = -omega*^2
print("(c) d theta''/d theta at theta* =", slope)
assert sp.simplify(slope + W**2 * (1 - (g / (R * W**2)) ** 2)) == 0             # = -Omega^2 sin^2(theta*)

vals = {m: 1.0, R: 0.5, g: 9.81, W: 8.0}
rhs = sp.lambdify((a, w), [w, acc_s.subs(vals)], modules="numpy")
E_f = sp.lambdify((a, w), state(th.diff(t) * Lag.diff(th.diff(t)) - Lag).subs(vals), modules="numpy")      # (3.4)
a_s = float(a_star.subs(vals))
w_s = 8.0 * np.sin(a_s)


def hit(t, y):
    return y[1]


hit.direction = -1                                                             # maxima of theta
sol = integrate.solve_ivp(lambda t, y: rhs(*y), (0, 5), [a_s + 0.01, 0.0], method="DOP853", rtol=1e-11, atol=1e-13, events=hit, t_eval=np.linspace(0, 5, 501))
T_ev = np.diff(sol.t_events[0]).mean()
print(f"(d) theta* = {a_s:.4f}, omega* = {w_s:.4f}; period from {len(sol.t_events[0]) - 1} events {T_ev:.6f} vs 2 pi/omega* {2 * np.pi / w_s:.6f}, relative difference {T_ev / (2 * np.pi / w_s) - 1:.1e}")
e = E_f(sol.y[0], sol.y[1])
lab = 0.5 * 0.5**2 * (sol.y[1] ** 2 + 64 * np.sin(sol.y[0]) ** 2) - 9.81 * 0.5 * np.cos(sol.y[0])      # lab-frame T + V for m = 1
print(f"(e) (3.4) drifts by {np.max(np.abs(e / e[0] - 1)):.1e}; lab-frame T+V varies by {np.max(np.abs(lab / lab[0] - 1)):.1e}")
assert abs(T_ev / (2 * np.pi / w_s) - 1) < 1e-3 and np.max(np.abs(e / e[0] - 1)) < 1e-10 and np.max(np.abs(lab / lab[0] - 1)) > 1e-3
```
```text
(a) theta'' = -(-Omega**2*R*cos(a) + g)*sin(a)/R
(b) equilibria: [0, pi, -acos(g/(Omega**2*R)) + 2*pi, acos(g/(Omega**2*R))]
(c) d theta''/d theta at theta* = -Omega**2 + g**2/(Omega**2*R**2)
(d) theta* = 1.2592, omega* = 7.6148; period from 6 events 0.825155 vs 2 pi/omega* 0.825128, relative difference 3.3e-05
(e) (3.4) drifts by 4.1e-14; lab-frame T+V varies by 3.2e-02
```

(a) The equation matches the hand result. (b) The equilibria are $\theta=0$, $\theta=\pi$ and $\theta=\pm\arccos(g/R\Omega^2)$, the last pair existing only if $g/(R\Omega^2)<1$ (here $0.31$). (c) $\partial\ddot\theta/\partial\theta$ at $\theta_*$ is $-\Omega^2+g^2/(\Omega^2R^2)=-\Omega^2\sin^2\theta_*$, since $\sin^2\theta_*=1-(g/R\Omega^2)^2$; so $\delta\ddot\theta=-\Omega^2\sin^2\theta_*\,\delta\theta$ and $\omega_*^2=\Omega^2\sin^2\theta_*$. (d) At amplitude $0.01$ rad the measured period agrees with $2\pi/\omega_*$ to 3e-05, the size of the amplitude correction. (e) (3.4) is conserved to rounding but is **not** the lab-frame $T+V$, which varies by percent: here (3.4) equals $\tfrac12mR^2(\dot\theta^2-\Omega^2\sin^2\theta)-mgR\cos\theta$, kinetic energy minus a centrifugal potential. The motor holding $\Omega$ fixed does work on the bead, so the lab energy changes, while (3.4) is conserved because $\mathcal L$ has no explicit time dependence (CM-4).

**B2.**

```python
import time

import numpy as np
from scipy import integrate, sparse


def lap1d(n):
    """(1/h^2) tridiag(-1, 2, -1), h = 1/(n+1), as in (3.12)."""
    h = 1.0 / (n + 1)
    return sparse.diags([-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1], format="csc") / h**2


n, t_end = 200, 0.1
x = np.arange(1, n + 1) / (n + 1)
A = lap1d(n)
u0, exact = np.sin(np.pi * x), np.exp(-np.pi**2 * t_end) * np.sin(np.pi * x)
lam_max = 4 * np.sin(n * np.pi / (2 * (n + 1))) ** 2 * (n + 1) ** 2             # Theorem 3.1, k = n
pred = t_end * lam_max / 3.3066
print(f"lambda_max = {lam_max:.4e}; RK45 needs at least {pred:.0f} steps to reach t = {t_end}")
runs = {}
for label, method, kw in [("RK45", "RK45", {}), ("Radau", "Radau", {}), ("BDF", "BDF", {}), ("LSODA", "LSODA", {}), ("BDF, jac=-A", "BDF", dict(jac=-A))]:
    calls = [0]

    def rhs(t, u):
        calls[0] += 1
        return -(A @ u)

    t0 = time.perf_counter()
    s = integrate.solve_ivp(rhs, (0, t_end), u0, method=method, rtol=1e-8, atol=1e-10, **kw)
    runs[label] = (s, calls[0])
    print(f"{label:12s} steps={len(s.t) - 1:5d} nfev={s.nfev:6d} counted calls={calls[0]:6d} njev={s.njev:2d} nlu={s.nlu:3d} error={np.abs(s.y[:, -1] - exact).max():.1e} time={time.perf_counter() - t0:.3f}s")
assert all(np.abs(s.y[:, -1] - exact).max() < 1e-5 for s, _ in runs.values())
assert abs((len(runs["RK45"][0].t) - 1) / pred - 1) < 0.02 and runs["RK45"][0].nfev / runs["BDF"][0].nfev > 100
assert runs["BDF, jac=-A"][1] == runs["BDF, jac=-A"][0].nfev and runs["BDF"][1] > runs["BDF"][0].nfev + 100
```
```text
lambda_max = 1.6159e+05; RK45 needs at least 4887 steps to reach t = 0.1
RK45         steps= 4871 nfev= 34130 counted calls= 34130 njev= 0 nlu=  0 error=7.5e-06 time=0.596s
Radau        steps=   32 nfev=   232 counted calls=   632 njev= 2 nlu=  6 error=7.5e-06 time=0.055s
BDF          steps=   44 nfev=    92 counted calls=   293 njev= 1 nlu= 12 error=7.5e-06 time=0.025s
LSODA        steps=   75 nfev=  2125 counted calls=  2125 njev=10 nlu= 10 error=7.5e-06 time=0.026s
BDF, jac=-A  steps=   44 nfev=    92 counted calls=    92 njev= 0 nlu= 12 error=7.5e-06 time=0.012s
```

(a) Theorem 3.1 with $k=n$ gives $\lambda_{\max}\approx4(n+1)^2\approx1.6\times10^5$, so the stability limit $3.31/\lambda_{\max}=2\times10^{-5}$ forces at least about 4900 steps; `RK45` took within 0.4% of that. The solution itself is a smooth, slowly decaying sine, and the step size is set entirely by the fastest mode of the discretization, which is stiffness. (b) The implicit methods need 30 to 75 steps and 16 to 370 times fewer right-hand-side calls (`LSODA` starts in its non-stiff mode and pays for the switch), with the same error: that error is the spatial discretization error, $\sim10^{-5}$, which no time integrator can remove. (c) Without `jac`, `BDF` reports `nfev=92`, but our wrapper counted 293 calls: `nfev` leaves out the calls spent on finite-difference Jacobians (201 here, one per column plus one). With `jac=-A` the counts agree and `njev` is 0. For implicit solvers `nfev` therefore under-reports the cost unless you supply the Jacobian, which for a linear problem is free.

**B3.**

```python
import numpy as np
import sympy as sp
from scipy import integrate

tt = sp.Symbol("t", real=True)
F = sp.Function("F")
ode = sp.Eq(F(tt).diff(tt, 2) + sp.Rational(2, 5) * F(tt).diff(tt) + sp.Rational(981, 100) * F(tt), 0)
cf = sp.dsolve(ode, F(tt), ics={F(0): sp.Rational(1, 20), F(tt).diff(tt).subs(tt, 0): 0})
print("(a)", cf)
closed = sp.lambdify(tt, cf.rhs, modules="numpy")

gam, gl = 0.4, 9.81


def rhs(t, y):
    return [y[1], -gam * y[1] - gl * np.sin(y[0])]


def maxima(t, y):
    return y[1]


maxima.direction = -1
te = np.linspace(0, 10, 1001)
s = integrate.solve_ivp(rhs, (0, 10), [0.05, 0.0], method="DOP853", rtol=1e-11, atol=1e-13, t_eval=te, events=maxima)
diff = np.max(np.abs(s.y[0] - closed(te)))
print(f"(b) largest |nonlinear - linear closed form| on [0, 10]: {diff:.1e}")
tm, ym = s.t_events[0], s.y_events[0][:, 0]
slope = np.polyfit(tm, np.log(ym), 1)[0]
w_d = np.sqrt(gl - gam**2 / 4)
print(f"(c) {len(tm)} maxima; decay slope {slope:.5f} (theory {-gam / 2}); mean spacing {np.diff(tm).mean():.5f} vs 2 pi/omega_d = {2 * np.pi / w_d:.5f}")
assert diff < 1e-4 and abs(slope + gam / 2) < 1e-3 and abs(np.diff(tm).mean() / (2 * np.pi / w_d) - 1) < 1e-4
```
```text
(a) Eq(F(t), (sqrt(977)*sin(sqrt(977)*t/10)/9770 + cos(sqrt(977)*t/10)/20)*exp(-t/5))
(b) largest |nonlinear - linear closed form| on [0, 10]: 2.4e-05
(c) 5 maxima; decay slope -0.20001 (theory -0.2); mean spacing 2.01027 vs 2 pi/omega_d = 2.01017
```

(a) `dsolve` returns $\theta(t)=e^{-t/5}\left(\tfrac{\sqrt{977}}{9770}\sin\tfrac{\sqrt{977}\,t}{10}+\tfrac1{20}\cos\tfrac{\sqrt{977}\,t}{10}\right)$, with $\omega_d=\sqrt{977}/10=\sqrt{g/L-\gamma^2/4}$. (b) The neglected force is $\theta^3/6\approx2\times10^{-5}$ at $\theta=0.05$ (and decays); the oscillator's static gain is $1/\omega_0^2\approx0.1$ and resonance can raise it by up to $Q=\omega_0/\gamma\approx8$, so a difference of order $10^{-5}$ is expected, and the run shows 2.4e-05. (c) The maxima decay as $e^{-\gamma t/2}$ and are spaced by the damped period, both to better than $10^{-3}$ relative: the linear theory is excellent at this amplitude. The event on $\omega=0$ with `direction=-1` catches exactly the maxima.

**B4.**

```python
import numpy as np
from scipy import optimize, special


def period(th0, L, g):
    return 4 * np.sqrt(L / g) * special.ellipk(np.sin(th0 / 2) ** 2)


Lm, sig = np.linspace(0.2, 2.0, 10), 0.01
T = period(0.1, Lm, 9.81) + sig * np.random.default_rng(0).standard_normal(10)
model = lambda L, a, b: a * L**b
sigma = np.full(10, sig)
p, cv = optimize.curve_fit(model, Lm, T, p0=[2.0, 0.5], sigma=sigma, absolute_sigma=True)
res = optimize.least_squares(lambda q: (model(Lm, *q) - T) / sig, [2.0, 0.5])
cv_310 = np.linalg.inv(res.jac.T @ res.jac)                                    # (3.10)
print("(a) a, b =", np.round(p, 5), "| sigma_a, sigma_b =", np.round(np.sqrt(np.diag(cv)), 5), "| max relative difference of pcov:", np.abs(cv_310 / cv - 1).max())
corr = cv[0, 1] / np.sqrt(cv[0, 0] * cv[1, 1])
wgt = T**2 / sig**2                                                              # weights for (ln a, b): (dT/d ln a)^2 / sigma^2
x_ = np.log(Lm)
corr_formula = -np.sum(wgt * x_) / np.sqrt(np.sum(wgt) * np.sum(wgt * x_**2))
print(f"(b) corr(a, b) = {corr:.3f}; weighted mean of ln L = {np.sum(wgt * x_) / np.sum(wgt):.3f}; formula {corr_formula:.3f}")
zb = (p[1] - 0.5) / np.sqrt(cv[1, 1])
g2, sg2 = 4 * np.pi**2 / p[0] ** 2, 8 * np.pi**2 * np.sqrt(cv[0, 0]) / p[0] ** 3
p1, cv1 = optimize.curve_fit(lambda L, gg: 2 * np.pi * np.sqrt(L / gg), Lm, T, p0=[9.0], sigma=sigma, absolute_sigma=True)
print(f"(c) b is {zb:+.2f} sigma from 1/2 | (d) g = {g2:.3f} +/- {sg2:.3f} (two parameters) vs {p1[0]:.3f} +/- {np.sqrt(cv1[0, 0]):.3f} (b fixed at 1/2)")
p2, cv2 = optimize.curve_fit(model, Lm, T, p0=[2.0, 0.5], sigma=sigma, absolute_sigma=False)
chi2 = np.sum(((model(Lm, *p) - T) / sig) ** 2)
print(f"(e) uncertainty ratio {np.sqrt(cv2[0, 0] / cv[0, 0]):.4f} = sqrt(chi2/dof) = sqrt({chi2:.2f}/8) = {np.sqrt(chi2 / 8):.4f}")
assert np.abs(cv_310 / cv - 1).max() < 1e-5 and abs(zb) < 2 and corr < 0 and abs(corr - corr_formula) < 0.01 and abs(np.sqrt(cv2[0, 0] / cv[0, 0]) - np.sqrt(chi2 / 8)) < 1e-6
```
```text
(a) a, b = [2.00872 0.49796] | sigma_a, sigma_b = [0.00339 0.00314] | max relative difference of pcov: 8.373118676185953e-08
(b) corr(a, b) = -0.456; weighted mean of ln L = 0.246; formula -0.457
(c) b is -0.65 sigma from 1/2 | (d) g = 9.784 +/- 0.033 (two parameters) vs 9.794 +/- 0.029 (b fixed at 1/2)
(e) uncertainty ratio 0.8003 = sqrt(chi2/dof) = sqrt(5.12/8) = 0.8003
```

(a) `curve_fit`'s covariance equals $(J^\top J)^{-1}$ from `least_squares` to $10^{-6}$ or better. (b) The correlation is negative. Linearize $\ln T=\ln a+b\ln L$: the fit "pivots" about the weighted mean of $\ln L$, where the weights are $T_i^2/\sigma^2$ (an *absolute* timing error makes long pendulums count more). That mean is positive, so $\operatorname{Cov}(\ln a,b)=-\bar x_w/S_{xx}<0$: raising $b$ must be compensated by lowering $a$. The closed form $-\sum w x/\sqrt{\sum w\sum wx^2}$ reproduces the number. (c) $b$ is well within $2\sigma$ of $\tfrac12$. (d) Freeing $b$ raises the uncertainty of $g$ from 0.029 to 0.033, because $a$ and $b$ are correlated: more freedom, larger error bars. (e) With `absolute_sigma=False` the covariance is scaled by $\chi^2_{\min}/(N-P)$, so the uncertainties scale by $\sqrt{\chi^2_{\min}/8}$; here $\chi^2_{\min}=5.1$ against an expected $8\pm4$, an ordinary fluctuation, so the bars shrink by 20% for no physical reason: with ten points the estimated scale is itself uncertain.

**B5.**

```python
import time

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla


def ham(n, box=10.0):
    """H = -1/2 d^2/dx^2 + x^2/2 on [-box, box], n interior points, zero ends; returns (H as CSC, h)."""
    x = np.linspace(-box, box, n + 2)[1:-1]
    h = x[1] - x[0]
    T = sparse.diags([-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1]) / (2 * h**2)
    return (T + sparse.diags(0.5 * x**2)).tocsc(), h


exact = np.arange(5) + 0.5
hs, errs = [], []
for n in (100, 200, 400, 800, 1600):
    H, h = ham(n)
    w = np.sort(spla.eigsh(H, k=5, sigma=0, which="LM", return_eigenvectors=False))
    hs.append(h)
    errs.append(np.abs(w - exact).max())
    print(f"n={n:5d} h={h:.4f} eigenvalues {np.round(w, 5)} worst error {errs[-1]:.2e}")
order = np.polyfit(np.log(hs), np.log(errs), 1)[0]
print(f"fitted order {order:.3f}")
H, h = ham(1000)
t0 = time.perf_counter()
w_s = np.sort(spla.eigsh(H, k=5, sigma=0, which="LM", return_eigenvectors=False))
t_s = time.perf_counter() - t0
t0 = time.perf_counter()
w_d = np.linalg.eigvalsh(H.toarray())[:5]
t_d = time.perf_counter() - t0
print(f"n=1000: eigsh {t_s * 1e3:.1f} ms, dense eigvalsh {t_d * 1e3:.0f} ms, eigenvalues differ by {np.abs(w_s - w_d).max():.1e}")
assert abs(order - 2) < 0.05 and np.abs(w_s - w_d).max() < 1e-8
```
```text
n=  100 h=0.1980 eigenvalues [0.49877 1.49385 2.48396 3.46909 4.44918] worst error 5.08e-02
n=  200 h=0.0995 eigenvalues [0.49969 1.49845 2.49597 3.49225 4.48728] worst error 1.27e-02
n=  400 h=0.0499 eigenvalues [0.49992 1.49961 2.49899 3.49806 4.49681] worst error 3.19e-03
n=  800 h=0.0250 eigenvalues [0.49998 1.4999  2.49975 3.49951 4.4992 ] worst error 7.99e-04
n= 1600 h=0.0125 eigenvalues [0.5     1.49998 2.49994 3.49988 4.4998 ] worst error 2.00e-04
fitted order 2.004
n=1000: eigsh 3.0 ms, dense eigvalsh 65 ms, eigenvalues differ by 1.5e-12
```

(a, b) The worst error falls by a factor of four each time $h$ halves: the fitted order is 2.00. This is step 2 of the Laplacian derivation: the stencil error is $-\tfrac{h^2}{12}u''''$, so the eigenvalues are second order too. The box truncation at $\pm10$ is negligible, since the wave functions there are $\sim e^{-50}$. (c) The shift-invert iteration on the sparse matrix finds the five lowest eigenvalues in 3.0 ms; the dense solver, which finds all 1000 eigenvalues of a matrix that is almost entirely zeros, takes 22 times longer for the same five numbers (to $10^{-8}$), and its cost grows as $n^3$ while the sparse cost grows roughly linearly. This is the QM-3 computation.

**B6.**

```python
import numpy as np
from scipy import fft, signal

fs, N = 1000.0, 500
t = np.arange(N) / fs
x = np.sin(2 * np.pi * 50 * t) + 0.5 * np.sin(2 * np.pi * 120.3 * t) + 0.2 * np.random.default_rng(0).standard_normal(N)
f = fft.rfftfreq(N, 1 / fs)
amp = 2 * np.abs(fft.rfft(x)) / N
pk, _ = signal.find_peaks(amp, height=0.2)
print("(a) peaks at", f[pk], "Hz with amplitudes", np.round(amp[pk], 3))

tone = 0.5 * np.sin(2 * np.pi * 120.3 * t)
win = signal.windows.hann(N, sym=False)
a_rect = 2 * np.abs(fft.rfft(tone)) / N
a_hann = 2 * np.abs(fft.rfft(tone * win)) / win.sum()
i200 = int(np.argmin(np.abs(f - 200)))
print(f"(b) 120.3 Hz tone alone: peak rectangular {a_rect.max():.4f}, Hann {a_hann.max():.4f} (true 0.5); leakage at {f[i200]:.0f} Hz: rectangular {a_rect[i200]:.1e}, Hann {a_hann[i200]:.1e}")
print(f"(c) bin spacing {fs / N} Hz; 120.3 Hz is {(120.3 - 120) / (fs / N):.2f} bins above bin {int(120 / (fs / N))}")
assert list(f[pk]) == [50.0, 120.0] and abs(amp[pk][0] - 1) < 0.05 and a_rect.max() < 0.49 and a_hann.max() > 0.49 and a_rect[i200] > 100 * a_hann[i200]
```
```text
(a) peaks at [ 50. 120.] Hz with amplitudes [1.004 0.482]
(b) 120.3 Hz tone alone: peak rectangular 0.4810, Hann 0.4928 (true 0.5); leakage at 200 Hz: rectangular 2.4e-03, Hann 1.2e-06
(c) bin spacing 2.0 Hz; 120.3 Hz is 0.15 bins above bin 60
```

(a) The peaks are at the bins 50 Hz and 120 Hz (bin spacing $f_s/N=2$ Hz) with amplitudes near 1.0 and 0.48. The 50 Hz tone completes exactly 25 cycles in the record, so it falls on a bin and its amplitude is right. The 120.3 Hz tone does not complete a whole number of cycles (60.15), so its energy spreads over neighbouring bins (leakage), and the bin it mostly occupies reads $|\operatorname{sinc}(0.15)|=0.96$ of the true amplitude (scalloping loss); noise adds $\sim\pm0.01$. (b) The Hann window brings the peak to within 1.5% of the true 0.5 and lowers the leakage at 200 Hz by a factor of 2051, at the price of a wider main lobe. (c) 120.3 Hz lies 0.15 bins above bin 60 (120 Hz). Windowing is how F3's convolution theorem shows up in practice: the spectrum you see is the true one convolved with the window's transform.


**B7.**

```python
import numpy as np
import sympy as sp
from scipy import integrate

t = sp.Symbol("t", real=True)
m, Ls, gs, gam = sp.symbols("m L g gamma", positive=True)
theta = sp.Function("theta")(t)
Lag = sp.Rational(1, 2) * m * Ls**2 * theta.diff(t) ** 2 + m * gs * Ls * sp.cos(theta)
Q = -gam * m * Ls**2 * theta.diff(t)                                           # generalized force of the damping
EL = sp.diff(Lag.diff(theta.diff(t)), t) - Lag.diff(theta) - Q                 # (3.3) with a force on the right
acc = sp.solve(EL, theta.diff(t, 2))[0]
print("(a) theta'' =", acc)
th, om = sp.symbols("th om", real=True)
state = lambda ex: ex.subs(theta.diff(t), om).subs(theta, th)                  # derivative first
rhs_num = sp.lambdify((th, om, gs, Ls, gam), [om, state(acc)], modules="numpy")
E_num = sp.lambdify((th, om, m, gs, Ls), state(theta.diff(t) * Lag.diff(theta.diff(t)) - Lag), modules="numpy")

g, L, gamma = 9.81, 1.0, 0.3


def rhs(t, y):
    """y = (theta, omega, W) with W' = omega^2, the dissipated-work integral."""
    return [*rhs_num(y[0], y[1], g, L, gamma), y[1] ** 2]


te = np.linspace(0, 20, 401)
s = integrate.solve_ivp(rhs, (0, 20), [2.0, 0.0, 0.0], method="DOP853", rtol=1e-11, atol=1e-13, t_eval=te)
e = E_num(s.y[0], s.y[1], 1.0, g, L)
resid = np.abs(e - e[0] + gamma * s.y[2]).max()
print(f"(b) energy falls from {e[0]:.4f} to {e[-1]:.4f}; max |e(t) - e(0) + gamma W(t)| = {resid:.1e}")
W_trap = integrate.cumulative_trapezoid(s.y[1] ** 2, te, initial=0)             # (c): quadrature on the 401 samples instead
resid_trap = np.abs(e - e[0] + gamma * W_trap).max()
print(f"(c) with trapezoid on the 401 samples the residual is {resid_trap:.1e}")
assert resid < 1e-8 and e[-1] < e[0] and resid_trap > 1e-5
```
```text
(a) theta'' = -gamma*Derivative(theta(t), t) - g*sin(theta(t))/L
(b) energy falls from 4.0824 to -9.7693; max |e(t) - e(0) + gamma W(t)| = 2.8e-10
(c) with trapezoid on the 401 samples the residual is 4.0e-03
```

(b) By hand, with $e=\tfrac12\omega^2-(g/L)\cos\theta$ and $\dot\omega=-\gamma\omega-(g/L)\sin\theta$: $\dot e=\omega\dot\omega+(g/L)\omega\sin\theta=-\gamma\omega^2-(g/L)\omega\sin\theta+(g/L)\omega\sin\theta=-\gamma\omega^2$. Integrating, $e(t)-e(0)=-\gamma\int_0^t\omega^2\,dt'=-\gamma W(t)$, which the run confirms to 3e-10. This checks the sign and size of the generalized-force term independently of the trajectory. (c) The extra state is integrated with the same adaptive error control as $\theta$ and $\omega$; `trapezoid` on the 401 `t_eval` samples has an $O(h^2)$ error set by the sampling and unrelated to `rtol`, and its residual is 4e-03, five orders of magnitude worse.

## C. Challenge

**C1.**

```python
import time

import numpy as np
from scipy import integrate

G, GAM = 9.81, 0.5


def steady_amplitude(Om, F, t_end=60.0, keep=6):
    """Mean and standard deviation of the last `keep` maxima of theta (events on omega = 0, downward)."""
    rhs = lambda t, y: [y[1], -GAM * y[1] - G * np.sin(y[0]) + F * np.cos(Om * t)]
    ev = lambda t, y: y[1]
    ev.direction = -1
    s = integrate.solve_ivp(rhs, (0, t_end), [0.0, 0.0], method="DOP853", rtol=1e-9, atol=1e-11, events=ev)
    m = s.y_events[0][:, 0][-keep:]
    return m.mean(), m.std()


Oms = np.arange(2.0, 4.0001, 0.05)
peaks = {}
for F in (0.3, 1.0, 1.5):
    t0 = time.perf_counter()
    res = np.array([steady_amplitude(O, F) for O in Oms])
    lin = F / np.sqrt((G - Oms**2) ** 2 + GAM**2 * Oms**2)
    i, j = res[:, 0].argmax(), lin.argmax()
    peaks[F] = (Oms[i], res[i, 0], Oms[j], lin[j])
    print(f"F={F}: peak {res[i, 0]:.3f} at Omega={Oms[i]:.2f} | linear peak {lin[j]:.3f} at {Oms[j]:.2f} | worst spread of the last maxima {res[:, 1].max():.1e} ({time.perf_counter() - t0:.1f} s)")
assert abs(peaks[0.3][1] / peaks[0.3][3] - 1) < 0.02 and peaks[0.3][0] == peaks[0.3][2]
assert peaks[1.0][0] < peaks[1.0][2] + 1e-9 and peaks[1.5][0] < peaks[1.0][0] and peaks[1.5][1] > peaks[1.5][3]
```
```text
F=0.3: peak 0.193 at Omega=3.10 | linear peak 0.192 at 3.10 | worst spread of the last maxima 3.8e-07 (2.3 s)
F=1.0: peak 0.654 at Omega=3.05 | linear peak 0.640 at 3.10 | worst spread of the last maxima 1.1e-05 (2.9 s)
F=1.5: peak 1.021 at Omega=2.90 | linear peak 0.960 at 3.10 | worst spread of the last maxima 6.3e-04 (2.9 s)
```

For $F=0.3$ the peak position agrees with the linear amplitude and its height to 0.3%. As $F$ grows, the peak moves to lower $\Omega$ and rises above the linear value: the pendulum is a *softening* oscillator, since $\sin\theta<\theta$ lowers the restoring force and the natural frequency falls with amplitude ($T/T_{\text{small}}\approx1+\theta_0^2/16$ from (3.7)), so the resonance bends toward lower frequency, and large amplitudes are reached because the drive stays resonant as the amplitude grows. The last maxima agree to $10^{-3}$ or better, so the amplitudes are steady; at larger $F$ the response has a jump (the curve is multivalued, so starting from rest samples only one branch) and eventually chaos, where "the amplitude" stops being a single number and events become a Poincaré-section tool (CM-8).

**C2.**

```python
import numpy as np
import sympy as sp
from scipy import integrate, optimize

fun = lambda x, y: np.vstack([y[1], -np.exp(y[0])])                             # u' = v, v' = -e^u
bc = lambda ya, yb: np.array([ya[0], yb[0]])
x = np.linspace(0, 1, 11)
guesses = {"flat": np.zeros((2, x.size)), "tall": np.vstack([12 * x * (1 - x), 12 * (1 - 2 * x)])}
sols = {}
for name, y_guess in guesses.items():
    for max_nodes in (1000, 100000):
        s = integrate.solve_bvp(fun, bc, x, y_guess, tol=1e-8, max_nodes=max_nodes)
        print(f"{name} guess, max_nodes={max_nodes}: success={s.success}, '{s.message}', u(1/2) = {s.sol(0.5)[0]:.6f}, u'(0) = {s.sol(0)[1]:.6f}, {s.x.size} nodes")
    sols[name] = s
for name, s in sols.items():
    up0 = s.sol(0)[1]
    shot = integrate.solve_ivp(lambda t, y: [y[1], -np.exp(y[0])], (0, 1), [0.0, up0], method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True)
    print(f"{name}: shooting from u'(0) = {up0:.6f} gives u(1) = {shot.y[0, -1]:.1e}, u(1/2) = {shot.sol(0.5)[0]:.6f}")
    assert abs(shot.y[0, -1]) < 1e-8 and abs(shot.sol(0.5)[0] - s.sol(0.5)[0]) < 1e-6 and s.success
assert sols["flat"].sol(0.5)[0] < 0.2 < 4 < sols["tall"].sol(0.5)[0]

X, Y = sp.symbols("X Y", real=True)
rosen = (1 - X) ** 2 + 100 * (Y - X**2) ** 2
f_r = sp.lambdify((X, Y), rosen, modules="numpy")
g_r = sp.lambdify((X, Y), [rosen.diff(X), rosen.diff(Y)], modules="numpy")
h_r = sp.lambdify((X, Y), sp.hessian(rosen, (X, Y)), modules="numpy")
fun_r, jac_r, hess_r = (lambda p: f_r(*p)), (lambda p: np.array(g_r(*p))), (lambda p: np.array(h_r(*p), dtype=float))
rows = {}
for label, kw in [("Nelder-Mead", dict(method="Nelder-Mead")), ("BFGS, finite differences", dict(method="BFGS")), ("BFGS, SymPy gradient", dict(method="BFGS", jac=jac_r)),
                  ("trust-exact, SymPy Hessian", dict(method="trust-exact", jac=jac_r, hess=hess_r))]:
    r = rows[label] = optimize.minimize(fun_r, [-1.2, 1.0], **kw)
    print(f"{label:27s} success={r.success} x=({r.x[0]:.6f}, {r.x[1]:.6f}) nfev={r.nfev:3d} njev={getattr(r, 'njev', 0):3d} nhev={getattr(r, 'nhev', 0):3d}")
    assert r.success and np.allclose(r.x, 1, atol=1e-4)
assert rows["Nelder-Mead"].nfev > rows["BFGS, finite differences"].nfev > rows["BFGS, SymPy gradient"].nfev > rows["trust-exact, SymPy Hessian"].nfev
```
```text
flat guess, max_nodes=1000: success=True, 'The algorithm converged to the desired accuracy.', u(1/2) = 0.140539, u'(0) = 0.549353, 121 nodes
flat guess, max_nodes=100000: success=True, 'The algorithm converged to the desired accuracy.', u(1/2) = 0.140539, u'(0) = 0.549353, 121 nodes
tall guess, max_nodes=1000: success=False, 'The maximum number of mesh nodes is exceeded.', u(1/2) = 4.091467, u'(0) = 10.846899, 537 nodes
tall guess, max_nodes=100000: success=True, 'The algorithm converged to the desired accuracy.', u(1/2) = 4.091467, u'(0) = 10.846899, 1047 nodes
flat: shooting from u'(0) = 0.549353 gives u(1) = -2.9e-12, u(1/2) = 0.140539
tall: shooting from u'(0) = 10.846899 gives u(1) = 1.1e-12, u(1/2) = 4.091467
Nelder-Mead                 success=True x=(1.000022, 1.000042) nfev=159 njev=  0 nhev=  0
BFGS, finite differences    success=True x=(0.999993, 0.999987) nfev=114 njev= 38 nhev=  0
BFGS, SymPy gradient        success=True x=(1.000000, 1.000000) nfev= 39 njev= 39 nhev=  0
trust-exact, SymPy Hessian  success=True x=(1.000000, 1.000000) nfev= 26 njev= 23 nhev= 26
```

(a) Two solutions exist (Bratu's problem with $\lambda=1$): $u(\tfrac12)=0.1405$ from the flat guess and $u(\tfrac12)=4.0915$ from the tall one. With `tol=1e-8` and the default `max_nodes=1000` the *tall* solution fails: `success=False` and the message says the maximum number of mesh nodes is exceeded, so the returned arrays are a half-refined mesh, not a solution; raising `max_nodes` gives `success=True`, and shooting from the returned $u'(0)$ reproduces $u(1)=0$ and $u(\tfrac12)$ to $10^{-6}$ or better. The lesson is the same as for `solve_ivp`: BVP solvers return what they have and tell you through `success`. Which solution you get is decided by the initial guess; a continuation in $\lambda$ is the systematic way to follow both branches. (b) Each extra piece of exact derivative information removes calls: Nelder–Mead uses only values; BFGS estimates gradients by finite differences, which costs $n$ extra evaluations each and shows in `nfev`; the SymPy gradient removes them and tightens the line searches; `trust-exact` with the Hessian takes Newton-like steps with a trust region and converges in the fewest iterations, though each iteration costs more. For a two-variable function the SymPy derivatives are free; for large problems you would write the gradient by hand or use automatic differentiation (PY-4's dual numbers are the forward-mode version).

**C3.**

```python
import time

import numpy as np
import sympy as sp


def derive(N):
    """N-link pendulum, unit masses and lengths: returns (derive time, count_ops, M and f as compiled functions, and the EL data for solve)."""
    t = sp.Symbol("t", real=True)
    g = sp.Symbol("g", positive=True)
    t0 = time.perf_counter()
    q = [sp.Function(f"q{i}")(t) for i in range(N)]
    x = [sum(sp.sin(q[i]) for i in range(k + 1)) for k in range(N)]
    y = [-sum(sp.cos(q[i]) for i in range(k + 1)) for k in range(N)]
    Lag = sum((x[k].diff(t) ** 2 + y[k].diff(t) ** 2) / 2 for k in range(N)) - g * sum(y)
    EL = [sp.diff(Lag.diff(qi.diff(t)), t) - Lag.diff(qi) for qi in q]
    acc = [qi.diff(t, 2) for qi in q]
    M, f = sp.linear_eq_to_matrix(EL, acc)
    t_der = time.perf_counter() - t0
    a, w = sp.symbols(f"a0:{N}", real=True), sp.symbols(f"w0:{N}", real=True)
    sub = lambda ex: ex.subs({qi.diff(t): wi for qi, wi in zip(q, w)}).subs({qi: ai for qi, ai in zip(q, a)})
    t0 = time.perf_counter()
    Mf = sp.lambdify((*a, *w, g), sub(M), modules="numpy", cse=True)
    ff = sp.lambdify((*a, *w, g), sub(f), modules="numpy", cse=True)
    return t_der, sp.count_ops(M) + sp.count_ops(f), time.perf_counter() - t0, Mf, ff, EL, acc


def numeric_rhs(a, w, g=9.81):
    """theta'' = M^{-1} f for unit masses and lengths: M_ij = (N - max(i,j)) cos(a_i - a_j), f_i = -sum_j (N - max(i,j)) sin(a_i - a_j) w_j^2 - g (N - i) sin(a_i)."""
    N = len(a)
    idx = np.arange(N)
    above = N - np.maximum.outer(idx, idx)                                       # total mass beyond joint max(i, j)
    d = a[:, None] - a[None, :]
    return np.linalg.solve(above * np.cos(d), -(above * np.sin(d)) @ w**2 - g * (N - idx) * np.sin(a))


print(" N  derive(s)  ops  lambdify(s)  call(us)  numpy call(us)  full solve(s), ops")
for N in range(1, 6):
    t_der, ops, t_lam, Mf, ff, EL, acc = derive(N)
    st = np.random.default_rng(N).uniform(-1, 1, 2 * N)
    t0 = time.perf_counter()
    for _ in range(1000):
        qdd = np.linalg.solve(np.array(Mf(*st, 9.81), float), np.array(ff(*st, 9.81), float).ravel())
    t_call = (time.perf_counter() - t0) / 1000
    t0 = time.perf_counter()
    for _ in range(1000):
        qdd_np = numeric_rhs(st[:N], st[N:])
    t_np = (time.perf_counter() - t0) / 1000
    assert np.allclose(qdd, qdd_np, atol=1e-12)                                  # the hand-written formula equals the SymPy pipeline
    full = ""
    if N <= 3:
        t0 = time.perf_counter()
        solved = sp.solve(EL, acc, dict=True)[0]
        full = f"{time.perf_counter() - t0:6.2f}, {sum(sp.count_ops(v) for v in solved.values())}"
    print(f"{N:2d} {t_der:9.2f} {ops:5d} {t_lam:11.2f} {t_call * 1e6:9.0f} {t_np * 1e6:14.0f}  {full}")
```
```text
 N  derive(s)  ops  lambdify(s)  call(us)  numpy call(us)  full solve(s), ops
 1      0.03    11        0.00        11             21    0.03, 12
 2      0.12   117        0.02        15             22    0.41, 1248
 3      0.18   326        0.05        19             24    7.15, 18921
 4      0.30   666        0.11        26             27  
 5      0.51  1165        0.18        26             25  
```

The table says: deriving and compiling $(M,\mathbf f)$ stays cheap (under a second at $N=5$; the expression grows polynomially, 117 to 1165 operations from $N=2$ to $5$), while the full `sp.solve` for the accelerations explodes: 0.4 s at $N=2$, seconds at $N=3$ with tens of thousands of operations against a few hundred for $(M,\mathbf f)$. We did not time $N=4$; the trend says minutes. Per call, the compiled $(M,\mathbf f)$ plus `np.linalg.solve` costs tens of microseconds. The structure is visible in the output: $M_{ij}=(N-\max(i,j))\cos(\theta_i-\theta_j)$ and $f_i=-\sum_j(N-\max(i,j))\sin(\theta_i-\theta_j)\omega_j^2-g(N-i)\sin\theta_i$, and the NumPy function built from it (a loop-free six lines) matches the SymPy pipeline to $10^{-12}$ at every $N$ tested and works for any $N$.

**Where to stop (about 100 words).** Use SymPy to *derive and check* at $N=1,2,3$: it removed every algebra error from the formulas above. Never ask `solve` to invert the mass matrix symbolically beyond $N=2$ or $3$: keep the equations as $M\ddot\theta=\mathbf f$ and let `np.linalg.solve` do the inversion numerically, which is also better conditioned. Once the structure is understood, write it directly in NumPy so it scales to any $N$ with no compile step, and keep the SymPy version as the test that the NumPy one is right, exactly as the asserts above do.
