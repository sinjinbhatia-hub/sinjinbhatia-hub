# PY-8: SciPy and SymPy — The Numerical and Symbolic Toolkits, and the Pipeline Between Them

**Time:** ~9 hours · **Prerequisites:** PY-7 (so PY-6's arrays and PY-3's functions-as-values); calculus and ODEs; LA-2 helps · **Unlocks:** CM-3 onward (the Lagrangian pipeline), EM-3, QM-3, ST-1, PY-10, and every project that calls `solve_ivp`, `quad`, `minimize` or `curve_fit`

> **Conventions in this module.** Every output below came from a real run of the block above it (Python 3.13, NumPy 2.5.3, SciPy 1.18.1, SymPy 1.14.0; the first block prints them). Blocks run top to bottom in one session, so names persist; timings are from one run, and only their ratios will transfer. Imports: `from scipy import optimize, integrate, linalg, sparse, interpolate, special, stats, fft`, `import sympy as sp`. `solve_ivp` always gets explicit `rtol` and `atol` and `t_eval` or `dense_output`: the defaults are for pictures. **Pendulum:** $\theta$ is measured from the downward vertical (so $\theta=0$ is stable equilibrium), the state is $\mathbf{y}=(\theta,\omega)$ with $\omega=\dot\theta$, $g=9.81\ \mathrm{m/s^2}$, and the *specific energy* is $e=\tfrac12\omega^2-(g/L)\cos\theta$ (energy divided by $mL^2$). **Elliptic integral:** $K(m)=\int_0^{\pi/2}d\varphi/\sqrt{1-m\sin^2\varphi}$ takes the *parameter* $m=k^2$, as does `special.ellipk(m)`. **DFT:** $F[k]=\sum_n f[n]e^{-2\pi i kn/N}$, no $1/N$, as in F3. Equations you refer back to are tagged $(3.n)$.

## 0. Where this sits

PY-3 taught you that a function is a value you can hand to another function; PY-6 gave you arrays; PY-7 plots them. SciPy is a cabinet of algorithms whose interface is "hand me a function, get a result object back"; SymPy is exact algebra on expression trees, where the equations come from. This module makes possible CM-3 (Lagrange's equations for any system, by the pipeline of §3.3), EM-3 (`special`'s Legendre and Bessel functions), QM-3 (finite-difference Schrödinger is §3.10's sparse eigenproblem) and ST-1 (`stats` distributions, `minimize` for maximum likelihood).

## 1. The problem that forced this idea

A pendulum is released from rest at 2.0 rad (115°). How long is its period, how do you know the integrator got it right, and what does a stopwatch experiment say about $g$? Then the same for a *double* pendulum, whose Euler–Lagrange equations are two coupled ODEs with a mass matrix full of $\cos(\theta_1-\theta_2)$: half a page by hand, with a sign error waiting. Three tools must cooperate: a symbolic engine for the algebra, an ODE solver for the dynamics, and PY-6/7's arrays and plots for the numbers. The skill is the pipeline between them, and knowing when a stage has silently failed.

Here is what failure looks like. This is the pendulum integrated with every default, then checked against the one thing we know must be constant:

```python
import sys
import time

import numpy as np
import sympy as sp
from scipy import optimize, integrate, linalg, sparse, interpolate, special, stats, fft
from scipy.sparse import linalg as spla

print(sys.version.split()[0], np.__version__, __import__("scipy").__version__, sp.__version__)

rng = np.random.default_rng(0)
g, L = 9.81, 1.0


def pend_rhs(t, y):
    """y = (theta, omega); theta'' = -(g/L) sin(theta)."""
    return [y[1], -(g / L) * np.sin(y[0])]


def spec_energy(y):
    """Specific energy 0.5*omega^2 - (g/L)*cos(theta); y has shape (2, n)."""
    return 0.5 * y[1] ** 2 - (g / L) * np.cos(y[0])


t_end = 20 * 2.67                                  # about 20 swings
sol = integrate.solve_ivp(pend_rhs, (0, t_end), [2.0, 0.0], t_eval=np.linspace(0, t_end, 401))
e = spec_energy(sol.y)
print(f"success={sol.success}, energy {e[0]:.4f} -> {e[-1]:.4f}, max relative change {np.max(np.abs(e / e[0] - 1)):.1%}")
```
```text
3.13.16 2.5.3 1.18.1 1.14.0
success=True, energy 4.0824 -> 3.4651, max relative change 15.8%
```

`success=True`, and the energy of a frictionless pendulum has changed by 16%. Nothing raised, nothing warned, and the plot of $\theta(t)$ would look fine. The rest of the module is how not to.

## 2. Intuition

**SciPy is a cabinet of drawers with one handle.** Every drawer takes a function you wrote and returns an object you must read:

| You call | Your function has the form | You get back |
|---|---|---|
| `optimize.root_scalar(f, bracket=...)` | `f(x)` | `RootResults`: `.root`, `.converged`, `.iterations`, `.function_calls` |
| `optimize.minimize(f, x0)` | `f(x)` returning a number | `OptimizeResult`: `.x`, `.fun`, `.success`, `.message`, `.nfev` |
| `integrate.solve_ivp(rhs, span, y0)` | `rhs(t, y)` returning $d\mathbf y/dt$ | `OdeResult`: `.t`, `.y`, `.success`, `.nfev`, `.t_events` |

Extra parameters travel by `args=(...)`, a closure, or `functools.partial` (PY-3). Result classes are dict subclasses with attribute access (PY-4). The routine calls your function dozens to millions of times, so *the call count is the cost*, and the result object reports it. A failed run usually does **not** raise; it returns `success=False` and a `message` you must read.

**SymPy manipulates trees.** `x**2 + sin(x)` is a node `Add` whose children are `Pow(x, 2)` and `sin(x)`; nothing is evaluated numerically, so `sp.sin(1)` stays `sin(1)` and `sp.Rational(1, 3)` is exactly one third. Every node is a Python object (PY-1), which is why SymPy is exact and slow. `lambdify` is a compiler: it prints the tree as Python source over NumPy functions and `exec`s it, giving an ordinary fast, approximate function. Symbolic cost is paid once; numeric cost a million times.

```text
 SYMBOLIC (exact, slow, once)                          NUMERIC (approximate, fast, many times)
 Lagrangian L(q, q')  --diff-->  Euler-Lagrange eqs  --solve-->  q'' = F(q, q')
                                                                    |  lambdify
                                                                    v
                                  rhs(t, y)  numpy function  --> solve_ivp --> OdeResult
                                                                                 |
                                                       energy check, events, fits, plots
```

**Three guesses to make before §3.** (i) Dividing `rtol` by 1000 should shrink the energy error by about that factor, and the error should grow with the number of periods: `solve_ivp` bounds the error *per step*, not a conserved quantity's drift. (ii) An explicit method can be forced to tiny steps by *stability*, not accuracy, long after a fast mode has died. (iii) A fit's uncertainty shrinks like $1/\sqrt N$, and speaks about noise *given that the model is right*.


## 3. Formal development

Each section states a tool's semantics, then hands you a build: **Build 1** the pendulum pipeline (§3.3), **2** the double pendulum (§3.6), **3** RK4 and the tolerance sweep (§3.4, Example 4.3), **4** stiffness (§3.5), **5** `curve_fit` for $g$ (§3.7), **6** the sparse Laplacian (§3.10), **7** `expm` against `solve_ivp` (§3.9), **8** SymPy sanity checks (§3.2, §3.6).

### 3.1 Function in, result object out

**Contract (3.1).** SciPy treats your function as a black box it may call anywhere in its domain, in any order, any number of times. So it must be pure (no hidden state) and return exactly the type and shape expected. The routine reports in a result object, and **a failed run usually returns instead of raising**: read `.success` (or `.converged`, or `.status`) and the call count every time.

▶ **Try it first.** PY-1 solved $\cos x = x$ by bisection (20 iterations to $10^{-6}$) and by Newton (4 iterations). Predict how many iterations `root_scalar` needs at `xtol=1e-12` for `bisect`, `brentq` and `newton` (with the derivative). Then run.

```python
f = lambda x: np.cos(x) - x                       # PY-1's fixed-point problem
df = lambda x: -np.sin(x) - 1
runs = {"bisect": dict(bracket=[0, 1]), "brentq": dict(bracket=[0, 1]), "newton": dict(x0=1.0, fprime=df)}
for method, kw in runs.items():
    r = optimize.root_scalar(f, method=method, xtol=1e-12, rtol=1e-12, **kw)
    print(f"{method:8s} converged={r.converged}  iterations={r.iterations:2d}  function_calls={r.function_calls:2d}  root={r.root:.15f}")
```
```text
bisect   converged=True  iterations=40  function_calls=42  root=0.739085133215667
brentq   converged=True  iterations= 7  function_calls= 8  root=0.739085133215156
newton   converged=True  iterations= 4  function_calls= 9  root=0.739085133215161
```

**■ What just happened.** Bisection gains one bit per call: $\lceil\log_2 10^{12}\rceil=40$ iterations. Newton squares its error each step (4 iterations) but its 9 calls include the derivative. `brentq` keeps bisection's guarantee and adds secant and inverse-quadratic steps: 7 iterations, 8 calls, no derivative. Compare `function_calls`, because in physics `f` may be an ODE solve costing seconds. All roots agree with PY-1's $0.7390851332151607$ to $\sim10^{-12}$.

### 3.2 SymPy in one sitting

**Semantics (3.2).** A SymPy expression is an immutable tree. `expr.func` is the node type, `expr.args` the children, and `expr == expr.func(*expr.args)`. Every operation returns a new tree. `Symbol("x", real=True)` attaches assumptions, and assumptions license simplifications. Python evaluates `1/3` *before* SymPy sees it.

▶ **Try it first.** Predict each printed line: `1/3 + Rational(1, 6)` and `Rational(1, 3) + Rational(1, 6)`; `sqrt(x**2)` for a plain, a real and a positive `x`; `solve(sin(x) - 2, x)` for a plain and a real `x`.

```python
x, y = sp.symbols("x y")
xr, xp = sp.Symbol("x", real=True), sp.Symbol("x", positive=True)
print(1 / 3 + sp.Rational(1, 6), "|", sp.Rational(1, 3) + sp.Rational(1, 6))
print(sp.sqrt(x**2), sp.sqrt(xr**2), sp.sqrt(xp**2))
print(sp.solve(sp.sin(x) - 2, x), sp.solve(sp.sin(xr) - 2, xr))
e = (x + y) ** 2
print(e.func, e.args, "|", sp.srepr(x * (x + 1)))
```
```text
0.500000000000000 | 1/2
sqrt(x**2) Abs(x) x
[pi - asin(2), asin(2)] []
<class 'sympy.core.power.Pow'> (x + y, 2) | Mul(Symbol('x'), Add(Symbol('x'), Integer(1)))
```

**■ What just happened.** `1/3` became the float `0.333...`, and one float makes the whole expression a float: the number-one SymPy gotcha. Write `sp.Rational(1, 3)` or `sp.S(1) / 3`. A bare symbol may be complex, so $\sqrt{x^2}$ stays and $\sin x=2$ has complex solutions; declare `real=True` or `positive=True` whenever you know it. `srepr` spells the tree as the constructor calls that rebuild it.

The targeted simplifiers do one predictable thing; `simplify` tries many and returns the shortest. Prefer the targeted ones in code you will rerun.

```python
for name, ex in [("expand", sp.expand((x + y) ** 3)), ("factor", sp.factor(x**3 - x)), ("cancel", sp.cancel((x**2 - 1) / (x - 1))),
                 ("together", sp.together(1 / x + 1 / y)), ("apart", sp.apart(1 / (x * (x + 1)))),
                 ("collect", sp.collect(x * y + x**2 * y + 3 * x, x)), ("trigsimp", sp.trigsimp(sp.sin(x) ** 4 - sp.cos(x) ** 4))]:
    print(f"{name:9s}{ex}")
```
```text
expand   x**3 + 3*x**2*y + 3*x*y**2 + y**3
factor   x*(x - 1)*(x + 1)
cancel   x + 1
together (x + y)/(x*y)
apart    -1/(x + 1) + 1/x
collect  x**2*y + x*(y + 3)
trigsimp -cos(2*x)
```

Solving comes in three flavors with different return types:

```python
print(sp.solve(x**2 - 2, x), "|", sp.solve(sp.sin(xr) - sp.Rational(1, 2), xr))
print(sp.solveset(sp.sin(xr) - sp.Rational(1, 2), xr, sp.S.Reals))
print(sp.solveset(sp.cos(xr) - xr, xr, sp.S.Reals), "|", sp.nsolve(sp.cos(x) - x, x, 1.0))
try:
    sp.solve(sp.cos(xr) - xr, xr)
except NotImplementedError as exc:
    print("solve raised NotImplementedError:", str(exc).splitlines()[0])
```
```text
[-sqrt(2), sqrt(2)] | [pi/6, 5*pi/6]
Union(ImageSet(Lambda(_n, 2*_n*pi + pi/6), Integers), ImageSet(Lambda(_n, 2*_n*pi + 5*pi/6), Integers))
ConditionSet(x, Eq(-x + cos(x), 0), Reals) | 0.739085133215161
solve raised NotImplementedError: multiple generators [x, cos(x)]
```

`solve` returns a *list* (principal solutions only), `solveset` a *set* (the full periodic family, or a `ConditionSet` meaning "I could not"), `nsolve` a number from Newton iteration at your guess. $\cos x=x$ has a root, yet `solve` raises and `solveset` shrugs: transcendental equations are `nsolve`'s or `root_scalar`'s job.

```python
t, f, k = sp.Symbol("t", real=True), sp.Function("f"), sp.Symbol("k", integer=True, positive=True)
print(sp.diff(x**3 * y, x, 2, y), "|", sp.integrate(x * sp.exp(-x), x), "|", sp.integrate(sp.exp(-(x**2)), (x, -sp.oo, sp.oo)))
print(sp.limit(sp.sin(x) / x, x, 0), "|", sp.series(sp.cos(x), x, 0, 6), "|", sp.summation(1 / k**2, (k, 1, sp.oo)))
print(sp.integrate(sp.sin(sp.sin(x)), x), "|", sp.Integral(sp.sin(sp.sin(x)), (x, 0, 1)).evalf())
A = sp.Matrix([[2, 1], [1, 2]])
print(A.eigenvals(), A.det(), A.inv().tolist(), sp.Matrix([[1, 2, 3], [2, 4, 6], [1, 0, 1]]).nullspace()[0].T)
print(sp.dsolve(f(t).diff(t, 2) + 2 * f(t).diff(t) + 5 * f(t), f(t)))
```
```text
6*x | (-x - 1)*exp(-x) | sqrt(pi)
1 | 1 - x**2/2 + x**4/24 + O(x**6) | pi**2/6
Integral(sin(sin(x)), x) | 0.430606103120691
{3: 1, 1: 1} 3 [[2/3, -1/3], [-1/3, 2/3]] Matrix([[-1, -1, 1]])
Eq(f(t), (C1*sin(2*t) + C2*cos(2*t))*exp(-t))
```

When `integrate` finds no closed form it returns the unevaluated `Integral`, and `.evalf()` then calls a numerical rule (§3.8). `Matrix` gives exact `eigenvals`, `det`, `inv`, `nullspace` (and `eigenvects`, `rref`), a check on LA-2 and LA-7 hand work; `dsolve` returns closed-form ODE solutions as an `Eq`.

The bridge is `lambdify`, which prints a tree as NumPy source and compiles it; `cse` factors out common subexpressions (`lambdify(..., cse=True)`, §3.6).

```python
expr = sp.sin(x) ** 2 + sp.sin(x) * y
fn = sp.lambdify((x, y), expr, modules="numpy")
print(fn(np.array([0.0, np.pi / 2]), 2.0), "|", sp.cse([sp.sin(x) ** 2 + sp.sin(x), sp.sin(x) * y]))
try:
    sp.lambdify(x, sp.sin(x), modules="math")(np.array([0.0, 1.0]))
except TypeError as exc:
    print("modules='math':", exc)
print(sp.lambdify(x, sp.Integer(5), modules="numpy")(np.array([1.0, 2.0])), "<- a constant expression gives a scalar, not an array")
sp.pprint(sp.sqrt(x) / (1 + x**2), use_unicode=False)
print(sp.latex(sp.sqrt(x) / (1 + x**2)))
```
```text
[0. 3.] | ([(x0, sin(x))], [x0**2 + x0, x0*y])
modules='math': only 0-dimensional arrays can be converted to Python scalars
5 <- a constant expression gives a scalar, not an array
  ___ 
\/ x  
------
 2    
x  + 1
\frac{\sqrt{x}}{x^{2} + 1}
```

With NumPy installed, SymPy's default `modules` targets NumPy; name it anyway, because `"math"` gives `math.sin`, which takes only scalars. A lambdified *constant* returns a bare scalar, which breaks a right-hand side that must return an array. `init_printing()` turns on `pprint`-style output in a notebook.


### 3.3 Build 1: the pendulum pipeline (the spine of the module)

Seven stages: the Lagrangian, Euler–Lagrange by `diff`, solve for $\ddot\theta$, trade functions of $t$ for plain symbols, `lambdify`, `solve_ivp`, then the energy and period checks. For a coordinate $q$ and $\mathcal L=T-V$, the Euler–Lagrange equation (derived in CM-3; here a recipe) and the energy function (CM-4) are

$$\frac{d}{dt}\frac{\partial\mathcal L}{\partial\dot q}-\frac{\partial\mathcal L}{\partial q}=0,\tag{3.3}$$

$$E=\sum_q \dot q\,\frac{\partial\mathcal L}{\partial\dot q}-\mathcal L.\tag{3.4}$$

▶ **Try it first.** Stage 4 swaps `theta(t)` for a plain symbol `th` and `theta(t).diff(t)` for `om`. If you substitute `th` first, what happens to the term $\dot\theta^2$ in the Lagrangian? Predict, then run.

```python
t = sp.Symbol("t", real=True)
m, Ls, gs = sp.symbols("m L g", positive=True)
theta = sp.Function("theta")(t)
Lag = sp.Rational(1, 2) * m * Ls**2 * theta.diff(t) ** 2 + m * gs * Ls * sp.cos(theta)          # T - V
EL = sp.diff(Lag.diff(theta.diff(t)), t) - Lag.diff(theta)                                     # (3.3)
acc = sp.solve(EL, theta.diff(t, 2))[0]
print("EL:", EL, "= 0   ->   theta'' =", acc)
th, om = sp.symbols("theta omega", real=True)
print("wrong order:", Lag.subs(theta, th), "-> doit():", Lag.subs(theta, th).doit(), "| right order:", Lag.subs(theta.diff(t), om).subs(theta, th))
```
```text
EL: L**2*m*Derivative(theta(t), (t, 2)) + L*g*m*sin(theta(t)) = 0   ->   theta'' = -g*sin(theta(t))/L
wrong order: L**2*m*Derivative(theta, t)**2/2 + L*g*m*cos(theta) -> doit(): L*g*m*cos(theta) | right order: L**2*m*omega**2/2 + L*g*m*cos(theta)
```

**■ What just happened.** The equation is $mL^2\ddot\theta+mgL\sin\theta=0$. Substituting the symbol first turned $\dot\theta$ into `Derivative(theta, t)` of a symbol, which is **zero**: `lambdify` refuses it, but `doit()` or `simplify` quietly deletes the kinetic energy. Substitute the highest derivative first: $\ddot q$, then $\dot q$, then $q$.

Stages 4 and 5 finish the bridge. The energy comes from (3.4), and the compiled right-hand side is checked against §1's `pend_rhs`:

```python
state = lambda ex: ex.subs(theta.diff(t), om).subs(theta, th)            # derivative first
rhs_num = sp.lambdify((th, om, gs, Ls), [om, state(acc)], modules="numpy")
E_sym = state(theta.diff(t) * Lag.diff(theta.diff(t)) - Lag)             # (3.4)
E_num = sp.lambdify((th, om, m, gs, Ls), E_sym, modules="numpy")


def pend_rhs_sym(t, y, g_=9.81, L_=1.0):
    """rhs(t, y) for y = (theta, omega), compiled from the Euler-Lagrange equation."""
    return rhs_num(y[0], y[1], g_, L_)


ys = rng.uniform(-3, 3, size=(5, 2))
print("E =", E_sym, "| rhs mismatch vs hand-written:", max(np.max(np.abs(np.array(pend_rhs_sym(0, y)) - np.array(pend_rhs(0, y)))) for y in ys),
      "| E equals spec_energy:", np.allclose(E_num(ys[:, 0], ys[:, 1], 1.0, g, 1.0), spec_energy(ys.T)))
```
```text
E = L**2*m*omega**2/2 - L*g*m*cos(theta) | rhs mismatch vs hand-written: 0.0 | E equals spec_energy: True
```

▶ **Try it first (the period).** The energy gives $\omega(\theta)$, and the period is four times the time from $\theta_0$ to $0$. Find a substitution that turns the integral into the elliptic integral $K$ (hint: $\sin(\theta/2)=k\sin\varphi$, $k=\sin(\theta_0/2)$). Then predict: for $\theta_0=2$, is $T/T_{\text{small}}$ above or below $1+\theta_0^2/16=1.25$?

**Derivation of the exact period.**

1. Released from rest, $e=-(g/L)\cos\theta_0$ for all time, so $\omega^2=\dfrac{2g}{L}(\cos\theta-\cos\theta_0)$. $\quad(3.5)$
2. By symmetry $\theta_0\to 0$ takes $T/4$, so $T=4\displaystyle\int_0^{\theta_0}\frac{d\theta}{|\omega|}=4\sqrt{\frac{L}{2g}}\int_0^{\theta_0}\frac{d\theta}{\sqrt{\cos\theta-\cos\theta_0}}$.
3. With $\cos\theta=1-2\sin^2(\theta/2)$: $\cos\theta-\cos\theta_0=2\left(k^2-\sin^2(\theta/2)\right)$.
4. Put $\sin(\theta/2)=k\sin\varphi$, $\varphi:0\to\pi/2$. Then $k^2-\sin^2(\theta/2)=k^2\cos^2\varphi$, and differentiating, $\tfrac12\cos(\theta/2)\,d\theta=k\cos\varphi\,d\varphi$ with $\cos(\theta/2)=\sqrt{1-k^2\sin^2\varphi}$, so $d\theta=\dfrac{2k\cos\varphi\,d\varphi}{\sqrt{1-k^2\sin^2\varphi}}$.
5. The integrand becomes $\dfrac{2k\cos\varphi}{\sqrt{1-k^2\sin^2\varphi}\cdot\sqrt2\,k\cos\varphi}\,d\varphi=\dfrac{\sqrt2\,d\varphi}{\sqrt{1-k^2\sin^2\varphi}}$, hence
$$T=4\sqrt{\frac{L}{g}}\,K\!\left(\sin^2\tfrac{\theta_0}{2}\right).\tag{3.6}$$ ∎
6. Expanding $K(m)=\tfrac{\pi}{2}\left(1+\tfrac m4+\tfrac{9m^2}{64}+\cdots\right)$ (SymPy check below) with $m=\sin^2(\theta_0/2)\approx\theta_0^2/4$: $T\approx2\pi\sqrt{L/g}\,\bigl(1+\theta_0^2/16\bigr)$. $\quad(3.7)$

```python
mm = sp.Symbol("m")
print(sp.series(sp.elliptic_k(mm), mm, 0, 3))
period = lambda th0_, L_=1.0, g_=9.81: 4 * np.sqrt(L_ / g_) * special.ellipk(np.sin(th0_ / 2) ** 2)
T_small = 2 * np.pi * np.sqrt(L / g)
for th0_ in [0.01, 0.5, 1.0, 2.0, 3.0]:
    print(f"theta0={th0_:5.2f}  T/T_small={period(th0_) / T_small:.5f}   1+theta0^2/16={1 + th0_**2 / 16:.5f}")
print("passing the modulus k instead of m:", 4 * np.sqrt(L / g) * special.ellipk(np.sin(1.0)), "vs correct", period(2.0))
```
```text
pi/2 + pi*m/8 + 9*pi*m**2/128 + O(m**3)
theta0= 0.01  T/T_small=1.00001   1+theta0^2/16=1.00001
theta0= 0.50  T/T_small=1.01585   1+theta0^2/16=1.01562
theta0= 1.00  T/T_small=1.06633   1+theta0^2/16=1.06250
theta0= 2.00  T/T_small=1.32890   1+theta0^2/16=1.25000
theta0= 3.00  T/T_small=2.57123   1+theta0^2/16=1.56250
passing the modulus k instead of m: 3.018449440931173 vs correct 2.665870942834926
```

**■ What just happened.** The period has a closed form, so the integrator can be tested to many digits. Every correction term is positive, so (3.7) undershoots more as $\theta_0$ grows, and $K\to\infty$ as $\theta_0\to\pi$. The last line is the classic error: `ellipk` wants $m=k^2$.

Stages 6 and 7: integrate with `DOP853` at `rtol=1e-10`, check the energy, and measure the period *from the solution* with an **event**. `solve_ivp` monitors an event function $h(t,\mathbf y)$; when it changes sign it locates the crossing by root-finding on the interpolant and records it in `t_events`. `direction=-1` keeps downward crossings only, so successive events are one period apart.

```python
th0, n_per = 2.0, 100
T_ex = period(th0)
t_eval = np.linspace(0, n_per * T_ex, 2001)
sol = integrate.solve_ivp(pend_rhs_sym, (0, n_per * T_ex), [th0, 0.0], method="DOP853", rtol=1e-10, atol=1e-12, t_eval=t_eval)
e_t = spec_energy(sol.y)
drift = np.abs(e_t / e_t[0] - 1)
print(f"success={sol.success}, nfev={sol.nfev}, max relative energy drift over {n_per} periods = {drift.max():.2e}")


def crossing(t, y):
    return y[0]                       # zero when theta passes through 0 ...


crossing.direction = -1               # ... moving downward


ev = integrate.solve_ivp(pend_rhs_sym, (0, 20 * T_ex), [th0, 0.0], method="DOP853", rtol=1e-10, atol=1e-12, events=crossing)
T_ev = np.diff(ev.t_events[0])
print(f"{len(T_ev)} periods from events: mean {T_ev.mean():.9f} s, spread {T_ev.std():.1e}, exact (3.6) {T_ex:.9f}, relative error {abs(T_ev.mean() / T_ex - 1):.1e}")
```
```text
success=True, nfev=58277, max relative energy drift over 100 periods = 8.55e-09
19 periods from events: mean 2.665870943 s, spread 1.8e-10, exact (3.6) 2.665870943, relative error 1.1e-10
```


### 3.4 What `rtol` and `atol` mean, and Build 3 (a hand-written RK4)

(The older `integrate.odeint` takes `f(y, t)` and has no events; use `solve_ivp`.) At each trial step an explicit *embedded* Runge–Kutta pair (`RK45`, the default, or `DOP853`) produces two answers of different order; their difference is an error estimate $\mathbf{err}$, and the step is accepted only if

$$\Bigl\|\,\mathrm{err}_i\Big/\bigl(\mathrm{atol}+\mathrm{rtol}\cdot\max(|y_i^{\text{old}}|,|y_i^{\text{new}}|)\bigr)\Bigr\|_{\text{RMS}}<1,\tag{3.8}$$

otherwise $h$ shrinks and the step is retried (it is the line `scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol` in SciPy's source). `rtol` asks for that many correct digits *per step*; `atol` is the size below which a component counts as zero, and matters whenever a component hovers near zero ($\omega$ at every turning point; a concentration of $10^{-6}$). Rule of thumb: `atol` $=$ `rtol` $\times$ the smallest magnitude you care about. Consequences: (a) the defaults `rtol=1e-3, atol=1e-6` ask for three digits per step, so §1's drift is not a bug; (b) (3.8) bounds the *local* error, and nothing bounds the accumulated error or the drift of a conserved quantity, so you measure it (guess (i) of §2). `max_step` caps the step for features the error estimate cannot see because no step samples them, such as a short pulse in a forcing term.

▶ **Try it first (Build 3).** Write `rk4(f, y0, t_end, n_steps)`, the classical fixed-step fourth-order method: $\mathbf y_{k+1}=\mathbf y_k+\tfrac h6(\mathbf k_1+2\mathbf k_2+2\mathbf k_3+\mathbf k_4)$, slopes taken at the start, twice at the midpoint, and at the end. Run 100 periods with $h=T/50$ (20 000 calls). Predict the worst relative energy drift: nearer $10^{-3}$, $10^{-6}$ or $10^{-9}$? Example 4.3 sets it against `solve_ivp`.

```python
def rk4(f, y0, t_end, n_steps):
    """Classical RK4 with a fixed step; returns the states at the n_steps+1 grid times, shape (len(y0), n_steps+1)."""
    h = t_end / n_steps
    y = np.array(y0, dtype=float)
    out = np.empty((y.size, n_steps + 1))
    out[:, 0] = y
    for k in range(n_steps):
        tk = k * h
        k1 = np.asarray(f(tk, y))
        k2 = np.asarray(f(tk + h / 2, y + h / 2 * k1))
        k3 = np.asarray(f(tk + h / 2, y + h / 2 * k2))
        k4 = np.asarray(f(tk + h, y + h * k3))
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[:, k + 1] = y
    return out


yr = rk4(pend_rhs, [th0, 0.0], n_per * T_ex, n_per * 50)
e_rk4 = spec_energy(yr)
rk4_drift = np.max(np.abs(e_rk4 / e_rk4[0] - 1))
print(f"RK4, h=T/50, {4 * n_per * 50} calls: max relative energy drift {rk4_drift:.2e}")
```
```text
RK4, h=T/50, 20000 calls: max relative energy drift 2.66e-03
```

**■ What just happened.** Nearer $10^{-3}$: nothing adapts $h$, so the drift is whatever $T/50$ gives; Example 4.3 compares.

### 3.5 Build 4: stiffness, where accuracy is not what limits the step

▶ **Try it first.** Take $y'=-\lambda\,(y-\cos t)$, $y(0)=0$, $\lambda=1000$. The solution is smooth after a millisecond. Predict how many steps `RK45` takes on $[0,10]$ at `rtol=1e-6`: about 100, 3000 or 100 000?

A transient $e^{-\lambda t}$ dies in milliseconds, then the solution follows $\cos t$ smoothly, so a step of $0.1$ would be accurate. Explicit methods cannot take it. For $y'=-\lambda y$, explicit Euler gives $y_{n+1}=(1-h\lambda)y_n$, which decays only if $h<2/\lambda$; backward Euler gives $y_{n+1}=y_n/(1+h\lambda)$, which decays for **every** $h>0$. $\quad(3.9)$ Any explicit method has a bounded stability region; for `RK45` the real-axis limit is where its one-step multiplier $R(z)$ reaches modulus 1. We read $R$ off the solver (one accepted step of $y'=zy$ from $y=1$) and find the limit with `brentq`:

```python
def R(z):
    """RK45's one-step multiplier: one forced-accept step of y' = z*y, h = 1, from y = 1."""
    s = integrate.RK45(lambda t, y: z * y, 0.0, [1.0], 1.0, first_step=1.0, rtol=1e10, atol=1e10)
    s.step()
    return s.y[0]


x_lim = optimize.brentq(lambda x: abs(R(-x)) - 1, 2.0, 5.0)
lam = 1000.0
stiff = lambda t, y: -lam * (y - np.cos(t))
exact = lambda t: (lam**2 * np.cos(t) + lam * np.sin(t) - lam**2 * np.exp(-lam * t)) / (1 + lam**2)
print(f"RK45 stability limit h*lambda < {x_lim:.4f}  ->  at least {10 * lam / x_lim:.0f} steps to reach t=10")
stiff_runs = {}
for method in ["RK45", "Radau", "LSODA"]:
    s = stiff_runs[method] = integrate.solve_ivp(stiff, (0, 10), [0.0], method=method, rtol=1e-6, atol=1e-9)
    print(f"{method:7s} success={s.success} steps={len(s.t) - 1:5d} nfev={s.nfev:6d} njev={s.njev:3d} mean h={np.diff(s.t).mean():.4f} error at t=10: {abs(s.y[0, -1] - exact(10.0)):.1e}")
```
```text
RK45 stability limit h*lambda < 3.3066  ->  at least 3024 steps to reach t=10
RK45    success=True steps= 3587 nfev= 21548 njev=  0 mean h=0.0028 error at t=10: 4.5e-07
Radau   success=True steps=   81 nfev=   624 njev=  3 mean h=0.1235 error at t=10: 2.1e-08
LSODA   success=True steps=  275 nfev=   552 njev= 42 mean h=0.0364 error at t=10: 6.3e-11
```

**■ What just happened.** `RK45` took 3587 steps of mean size 0.0028, a little under its stability limit $3.31/\lambda=0.0033$, although the answer is smooth enough for steps 44 times larger (`Radau`'s mean step); that is 21548 right-hand-side calls. `Radau` (implicit: each step solves a nonlinear equation using the Jacobian, counted by `njev`) took 81 steps and 624 calls, 35 times fewer; `LSODA` switches between a non-stiff and a stiff method by itself. All errors are small: the explicit solver pays for *stability*, not accuracy. Diagnosis: a smooth answer with a huge `nfev` or tiny mean step means stiffness; switch to `Radau`, `BDF` or `LSODA`, passing `jac=` if you have it (SymPy can differentiate your right-hand side). The project repeats the test on Robertson's kinetics.


### 3.6 Build 2: the double pendulum, and Build 8 (why `lambdify` exists)

Two point masses $m_1,m_2$ on massless rods $l_1,l_2$, angles $\theta_1,\theta_2$ from the downward vertical. The Lagrangian is four lines of kinematics; the equations of motion are what you do not want to derive twice. The pipeline of §3.3 only gains a coordinate.

```python
q1, q2 = sp.Function("theta1")(t), sp.Function("theta2")(t)
M1, M2, l1, l2 = sp.symbols("m1 m2 l1 l2", positive=True)
x1, y1 = l1 * sp.sin(q1), -l1 * sp.cos(q1)
x2, y2 = x1 + l2 * sp.sin(q2), y1 - l2 * sp.cos(q2)
speed2 = lambda xx, yy: xx.diff(t) ** 2 + yy.diff(t) ** 2
Lag2 = sp.simplify(sp.Rational(1, 2) * M1 * speed2(x1, y1) + sp.Rational(1, 2) * M2 * speed2(x2, y2) - gs * (M1 * y1 + M2 * y2))
qs = [q1, q2]
EL2 = [sp.diff(Lag2.diff(q.diff(t)), t) - Lag2.diff(q) for q in qs]                       # (3.3) per coordinate
accs = [q.diff(t, 2) for q in qs]
sol2 = sp.solve(EL2, accs, dict=True)[0]
Mmat, _ = sp.linear_eq_to_matrix(EL2, accs)
print("mass matrix M(theta):", sp.simplify(Mmat).tolist())
print("operations in the two accelerations:", [sp.count_ops(sol2[a]) for a in accs])

a1, a2, w1, w2 = sp.symbols("a1 a2 w1 w2", real=True)
state2 = lambda ex: ex.subs({q1.diff(t): w1, q2.diff(t): w2}).subs({q1: a1, q2: a2})       # derivatives first
args2 = (a1, a2, w1, w2, M1, M2, l1, l2, gs)
rhs2 = sp.lambdify(args2, [state2(sol2[a]) for a in accs], modules="numpy", cse=True)
E2 = sp.lambdify(args2, state2(sum(q.diff(t) * Lag2.diff(q.diff(t)) for q in qs) - Lag2), modules="numpy", cse=True)


def dp_rhs(t_, y):
    return [y[2], y[3], *rhs2(y[0], y[1], y[2], y[3], 1.0, 1.0, 1.0, 1.0, g)]


y0 = np.array([2.0, 2.0, 0.0, 0.0])
te = np.linspace(0, 30, 3001)
runs2 = [integrate.solve_ivp(dp_rhs, (0, 30), y0 + d, method="DOP853", rtol=1e-12, atol=1e-12, t_eval=te) for d in (0.0, 1e-8 * np.array([1, 0, 0, 0]))]
E_a = E2(*runs2[0].y, 1.0, 1.0, 1.0, 1.0, g)
sep = np.linalg.norm(runs2[0].y - runs2[1].y, axis=0)
win = (sep > 1e-7) & (sep < 1e-1)
slope = np.polyfit(te[win], np.log(sep[win]), 1)[0]
print(f"energy drift {np.max(np.abs(E_a / E_a[0] - 1)):.1e}; separation of two starts 1e-8 apart:",
      ", ".join(f"t={tt}: {sep[np.searchsorted(te, tt)]:.1e}" for tt in (0, 5, 10, 15, 20)))
print(f"log-separation grows at {slope:.2f} per second between 1e-7 and 1e-1")
```
```text
mass matrix M(theta): [[l1**2*(m1 + m2), l1*l2*m2*cos(theta1(t) - theta2(t))], [l1*l2*m2*cos(theta1(t) - theta2(t)), l2**2*m2]]
operations in the two accelerations: [106, 151]
energy drift 6.1e-11; separation of two starts 1e-8 apart: t=0: 1.0e-08, t=5: 3.2e-06, t=10: 3.6e-04, t=15: 8.6e-01, t=20: 1.5e+01
log-separation grows at 1.20 per second between 1e-7 and 1e-1
```

The accelerations are 106 and 151 operations long, and the system has the form $M(\boldsymbol\theta)\ddot{\boldsymbol\theta}=\mathbf f(\boldsymbol\theta,\dot{\boldsymbol\theta})$ that CM-3 derives in general. Two checks that the machine-derived equations are right: the project compares them with the published hand-derived form, and the energy from (3.4) holds to a relative 6e-11 over 30 s. Two starts $10^{-8}$ apart separate exponentially, at about 1.2 per second, until the separation saturates: chaos in one table (the project plots it). The slope is a finite-time, single-direction estimate, not a Lyapunov exponent, and it depends on the initial condition.

**Build 8: `lambdify` against `subs`.** Evaluate the first acceleration at one state, symbolically and compiled:

```python
acc1 = state2(sol2[accs[0]]).subs({M1: 1, M2: 1, l1: 1, l2: 1, gs: g})
f_acc1 = sp.lambdify((a1, a2, w1, w2), acc1, modules="numpy", cse=True)
pt = {a1: 1.0, a2: 0.5, w1: 0.3, w2: -0.2}
t0 = time.perf_counter()
for _ in range(200):
    v_sym = acc1.subs(pt).evalf()
t_sub = (time.perf_counter() - t0) / 200
t0 = time.perf_counter()
for _ in range(20000):
    v_num = f_acc1(1.0, 0.5, 0.3, -0.2)
t_lam = (time.perf_counter() - t0) / 20000
big = rng.uniform(-1, 1, size=(4, 10**6))
t0 = time.perf_counter()
f_acc1(*big)
t_vec = (time.perf_counter() - t0) / 10**6
print(f"subs+evalf {t_sub * 1e3:.2f} ms/call | lambdify {t_lam * 1e6:.1f} us/call, {t_vec * 1e9:.0f} ns/point on 10^6 points | "
      f"same value: {abs(float(v_sym) - v_num) < 1e-12} | speedups {t_sub / t_lam:.0f}x, {t_sub / t_vec:.0f}x")
```
```text
subs+evalf 0.13 ms/call | lambdify 1.6 us/call, 333 ns/point on 10^6 points | same value: True | speedups 86x, 403x
```



### 3.7 `optimize`, and Build 5 (`curve_fit` for $g$)

`minimize(fun, x0, method=..., jac=..., bounds=..., constraints=...)` is local: it walks downhill from `x0`, so *which* minimum you reach depends on `x0`. Nelder–Mead uses values only (a shrinking simplex: slow, noise-tolerant). BFGS builds a Hessian approximation from gradients; without `jac=` it uses finite differences, $n$ extra calls per gradient, counted in `nfev`. The pipeline supplies the exact gradient (SymPy `diff`, `lambdify`, `jac=`). Himmelblau's function has four minima:

```python
X, Y = sp.symbols("X Y", real=True)
himmel = (X**2 + Y - 11) ** 2 + (X + Y**2 - 7) ** 2
f_h = sp.lambdify((X, Y), himmel, modules="numpy")
grad_h = sp.lambdify((X, Y), [himmel.diff(X), himmel.diff(Y)], modules="numpy")
fun, jac = (lambda p: f_h(*p)), (lambda p: np.array(grad_h(*p)))
for label, kw in [("Nelder-Mead", dict(method="Nelder-Mead")), ("BFGS, finite differences", dict(method="BFGS")),
                  ("BFGS, SymPy gradient", dict(method="BFGS", jac=jac))]:
    r = optimize.minimize(fun, [0.0, 0.0], **kw)
    print(f"{label:25s} success={r.success} x=({r.x[0]:.5f}, {r.x[1]:.5f}) nfev={r.nfev:3d} njev={getattr(r, 'njev', 0):2d}")
print("four starts, four minima:", [np.round(optimize.minimize(fun, s0, jac=jac, method="BFGS").x, 2).tolist() for s0 in [(-4, -4), (-4, 4), (4, -4), (4, 4)]])
de = optimize.differential_evolution(fun, [(-5, 5), (-5, 5)], seed=0)
c = optimize.minimize(lambda p: p @ p, [2.0, 0.0], method="SLSQP", constraints=[{"type": "ineq", "fun": lambda p: p[0] + p[1] - 1}], bounds=[(0, 2), (0, 2)])
rt = optimize.root(lambda p: [p[0] ** 2 + p[1] ** 2 - 4, p[1] - p[0] ** 2], [1.0, 1.0])
print("global search found", np.round(de.x, 3), "| constrained", np.round(c.x, 4), c.success, "| root of a 2x2 system", np.round(rt.x, 4), rt.success)
xs = np.linspace(0, 1, 20)
ys = 2 * xs + 1 + 0.05 * rng.standard_normal(20)
ys[[3, 15]] += 2.0                                                       # two outliers
for loss in ["linear", "soft_l1"]:
    print(f"least_squares loss={loss:8s} slope, intercept =", np.round(optimize.least_squares(lambda p: p[0] * xs + p[1] - ys, [0.0, 0.0], loss=loss, f_scale=0.1).x, 3))
```
```text
Nelder-Mead               success=True x=(3.00001, 1.99997) nfev=157 njev= 0
BFGS, finite differences  success=True x=(3.00000, 2.00000) nfev= 48 njev=16
BFGS, SymPy gradient      success=True x=(3.00000, 2.00000) nfev= 16 njev=16
four starts, four minima: [[-3.78, -3.28], [-2.81, 3.13], [3.58, -1.85], [3.0, 2.0]]
global search found [ 3.584 -1.848] | constrained [0.5 0.5] True | root of a 2x2 system [1.2496 1.5616] True
least_squares loss=linear   slope, intercept = [1.926 1.232]
least_squares loss=soft_l1  slope, intercept = [1.967 1.025]
```

Constraints go through `SLSQP` or `trust-constr`; `differential_evolution` is a *global* search that needs bounds and many more calls; `root` solves $F(\mathbf x)=\mathbf 0$; `minimize_scalar` is the one-variable version (Brent, or `method="bounded"`); `loss="soft_l1"` down-weights the outliers a squared loss chases.

**`curve_fit`.** `popt, pcov = curve_fit(model, x, y, p0=..., sigma=..., absolute_sigma=...)` minimizes $\chi^2(\mathbf p)=\sum_i r_i^2$, $r_i=(y_i-m(x_i;\mathbf p))/\sigma_i$. What does `pcov` mean?

▶ **Try it first.** Let $y_i=m(x_i;\mathbf p^*)+\varepsilon_i$ with independent zero-mean noise of standard deviation $\sigma_i$. Linearize $\mathbf r(\mathbf p)$ about $\mathbf p^*$, minimize the quadratic, and find the covariance of the minimizer.

**Derivation.** (1) At the truth $\mathbf r(\mathbf p^*)=\boldsymbol\eta$ with $\eta_i=\varepsilon_i/\sigma_i$, so $E[\boldsymbol\eta\boldsymbol\eta^\top]=I$. (2) Near $\mathbf p^*$, $\mathbf r(\mathbf p)\approx\boldsymbol\eta+J\,\delta$, $\delta=\mathbf p-\mathbf p^*$, with $J_{ij}=\partial r_i/\partial p_j$. (3) Minimizing $\|\boldsymbol\eta+J\delta\|^2$ is linear least squares (LA-5): $\hat\delta=-(J^\top J)^{-1}J^\top\boldsymbol\eta$. (4) Hence
$$\operatorname{Cov}(\hat{\mathbf p})=(J^\top J)^{-1}J^\top E[\boldsymbol\eta\boldsymbol\eta^\top]J(J^\top J)^{-1}=(J^\top J)^{-1}.\tag{3.10}$$ ∎

This is ST-3's $\sigma^2(X^\top X)^{-1}$ with $J$ in the role of $X$. With `absolute_sigma=False` (the **default**) `sigma` only sets relative weights and `pcov` is multiplied by the reduced $\chi^2=\chi^2_{\min}/(N-P)$ from the residuals. **Assumptions behind $\sqrt{\mathrm{pcov}_{ii}}$:** the model family contains the truth (else the bar describes scatter about the wrong answer), the noise is independent with the stated $\sigma_i$, the model is near-linear over $\pm\sigma_{\mathbf p}$, and Gaussian noise is needed only to read $\pm1\sigma$ as 68%.

Build 5: ten lengths from 0.2 to 2.0 m, true $g=9.81$, amplitude 0.5 rad, timing noise 0.01 s. Fit the small-angle model and the exact model (3.6):

```python
Ls_m, th_amp, sig_T, g_true = np.linspace(0.2, 2.0, 10), 0.5, 0.01, 9.81
rng_fit = np.random.default_rng(0)                                   # a fresh stream for this experiment
T_meas = period(th_amp, Ls_m, g_true) + sig_T * rng_fit.standard_normal(Ls_m.size)
models = {"small-angle": lambda Lm, gg: 2 * np.pi * np.sqrt(Lm / gg), "exact (3.6)": lambda Lm, gg: period(th_amp, Lm, gg)}
fits = {}
for name, model in models.items():
    fits[name] = optimize.curve_fit(model, Ls_m, T_meas, p0=[8.0], sigma=np.full(Ls_m.size, sig_T), absolute_sigma=True)
    p, cv = fits[name]
    print(f"{name:12s} g = {p[0]:.3f} +/- {np.sqrt(cv[0, 0]):.3f}  ({(p[0] - g_true) / np.sqrt(cv[0, 0]):+.1f} sigma from the truth)")
Lq, gq = sp.symbols("Lq gq", positive=True)
dTdg = sp.lambdify((Lq, gq), sp.diff(2 * sp.pi * sp.sqrt(Lq / gq), gq), modules="numpy")      # exact Jacobian from SymPy
Jm = (-dTdg(Ls_m, fits["small-angle"][0][0]) / sig_T).reshape(-1, 1)                           # r = (y - m)/sigma
print("sqrt of (3.10) from the SymPy Jacobian:", np.sqrt(np.linalg.inv(Jm.T @ Jm)[0, 0]), "| curve_fit:", np.sqrt(fits["small-angle"][1][0, 0]))
p_rel, cv_rel = optimize.curve_fit(models["exact (3.6)"], Ls_m, T_meas, p0=[8.0])
print(f"default absolute_sigma=False, no sigma given: g = {p_rel[0]:.3f} +/- {np.sqrt(cv_rel[0, 0]):.3f} (scaled by the residuals' own scatter)")
```
```text
small-angle  g = 9.503 +/- 0.028  (-10.9 sigma from the truth)
exact (3.6)  g = 9.806 +/- 0.029  (-0.1 sigma from the truth)
sqrt of (3.10) from the SymPy Jacobian: 0.0281135249709341 | curve_fit: 0.028113522667159954
default absolute_sigma=False, no sigma given: g = 9.806 +/- 0.023 (scaled by the residuals' own scatter)
```

**■ What just happened.** The small-angle model is wrong by $(1+\theta_0^2/16)^{-2}\approx-3\%$ in $g$: it lands many reported $\sigma$ away while printing a tiny error bar, because the covariance knows about noise, not about a wrong model. The exact model recovers $g$ within its uncertainty. Formula (3.10) with our SymPy Jacobian reproduces `curve_fit`'s number. Example 4.5 tests the claim by repetition.


### 3.8 `integrate.quad` and friends

`quad(f, a, b)` is QUADPACK's adaptive Gauss–Kronrod rule: a 21-point rule is compared with its embedded 10-point Gauss rule, and subintervals where they disagree are bisected (at most `limit=50`, else an `IntegrationWarning`). It returns `(value, abserr)`, where `abserr` is an *estimate* from that disagreement. Infinite limits are mapped to finite ones by a change of variables. `dblquad(func, a, b, gfun, hfun)` integrates `func(y, x)`, **inner variable first**; `nquad` does $n$ dimensions. For samples, use `integrate.trapezoid`, `simpson` or `cumulative_trapezoid` (`np.trapz` is gone; use `np.trapezoid`). PY-5's `adaptive_simpson` (`integrators.py`, used in its Example 4.6 and §8) is the same idea with a 3-point rule. Example 4.4 tests `abserr`.

### 3.9 `linalg` beyond NumPy, and Build 7 (`expm` against `solve_ivp`)

`scipy.linalg` adds the factorizations (`lu`, `qr`, `svd`, `eigh`), structured solvers (`solve_triangular`, `solve_banded`), subspaces (`null_space`, `orth`) and `expm`.

```python
Mx = rng.standard_normal((5, 5))
S = Mx + Mx.T
bv = rng.standard_normal(5)
P, Lo, Up = linalg.lu(Mx)                              # Mx = P @ Lo @ Up
Q, Rr = linalg.qr(Mx)
U_, sv, Vh = linalg.svd(Mx)
w_, V_ = linalg.eigh(S)
tri = np.diag(2.0 * np.ones(5)) - np.diag(np.ones(4), 1) - np.diag(np.ones(4), -1)
ab = np.zeros((3, 5)); ab[0, 1:], ab[1], ab[2, :-1] = -1.0, 2.0, -1.0           # banded storage: super, main, sub
rank1 = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])
res = [P @ Lo @ Up - Mx, Q @ Rr - Mx, (U_ * sv) @ Vh - Mx, S @ V_ - V_ * w_, Up @ linalg.solve_triangular(Up, bv) - bv,
       linalg.solve_banded((1, 1), ab, bv) - np.linalg.solve(tri, bv)]
print("lu, qr, svd, eigh, solve_triangular, solve_banded: worst residual", max(np.abs(r).max() for r in res))
print("rank-1 3x2 matrix: null_space shape", linalg.null_space(rank1).shape, "| orth shape", linalg.orth(rank1).shape)
```
```text
lu, qr, svd, eigh, solve_triangular, solve_banded: worst residual 1.1102230246251565e-14
rank-1 3x2 matrix: null_space shape (2, 1) | orth shape (3, 1)
```

▶ **Try it first.** For $\dot{\mathbf x}=A\mathbf x$, differentiate $\sum_k(At)^k/k!$ term by term and decide which ODE it solves.

**The matrix exponential.** For $\dot{\mathbf x}=A\mathbf x$, $\mathbf x(0)=\mathbf x_0$, define $e^{At}=\sum_{k\ge0}(At)^k/k!$ (converges for every square $A$). Differentiating term by term, $\frac{d}{dt}e^{At}=\sum_{k\ge1}A\,\frac{(At)^{k-1}}{(k-1)!}=Ae^{At}$ and $e^{A\cdot0}=I$, so $\mathbf x(t)=e^{At}\mathbf x_0$ solves the problem, and by uniqueness for linear ODEs it is *the* solution. $\quad(3.11)$ LA-7 computes it as $Ve^{\Lambda t}V^{-1}$ for diagonalizable $A$; LA-11 handles the rest. Build 7: a damped oscillator $\ddot x+2\zeta\omega_0\dot x+\omega_0^2x=0$ as a first-order system, exact by (3.11) and numerical by `solve_ivp`:

```python
w0, zeta = 3.0, 0.1
Aosc = np.array([[0.0, 1.0], [-(w0**2), -2 * zeta * w0]])
x0v = np.array([1.0, 0.0])
tt = np.linspace(0, 10, 201)
x_exp = np.array([linalg.expm(Aosc * ti) @ x0v for ti in tt]).T
x_ivp = integrate.solve_ivp(lambda ti, y: Aosc @ y, (0, 10), x0v, method="DOP853", rtol=1e-12, atol=1e-14, t_eval=tt).y
print(f"expm vs solve_ivp (rtol 1e-12): max |difference| = {np.abs(x_exp - x_ivp).max():.1e}")
```
```text
expm vs solve_ivp (rtol 1e-12): max |difference| = 2.4e-12
```

**■ What just happened.** The series solves $\dot{\mathbf x}=A\mathbf x$, so `expm` is the exact solution operator, and a tight `solve_ivp` matches it to $2\times10^{-12}$.

### 3.10 `sparse`, and Build 6 (the Laplacian)

A finite-difference operator has a few nonzeros per row; storing $n^2$ numbers wastes memory and time. The formats differ in what they make cheap:

| Format | Stores | Use it for |
|---|---|---|
| COO | `(row, col, value)` triples; duplicates are **summed** | assembly; convert after |
| CSR / CSC | rows / columns contiguous | CSR: row slices and `A @ v`; CSC: what `spsolve` and `splu` factor |
| DIA | whole diagonals | banded operators (`sparse.diags`) |
| LIL / DOK | cheap element assignment | slow incremental building, then `.tocsr()` |

Sparse **matrices are not arrays**, and `*` differs between the two families:

```python
Bm = sparse.csr_matrix(sparse.diags([np.ones(2), 2 * np.ones(3), np.ones(2)], [-1, 0, 1]))
Ba = sparse.csr_array(Bm)
v3 = np.array([1.0, 2.0, 3.0])
print("matrix * vec:", Bm * v3, "| array * vec gives a", type(Ba * v3).__name__, "(elementwise)")
print("matrix: * is @:", np.allclose((Bm * Bm).toarray(), (Bm @ Bm).toarray()), "| array: * is @:", np.allclose((Ba * Ba).toarray(), (Ba @ Ba).toarray()))
print("sparse.diags returns", type(sparse.diags([np.ones(2)], [1])).__name__, "| sparse.diags_array returns", type(sparse.diags_array([np.ones(2)], offsets=[1])).__name__)
```
```text
matrix * vec: [4. 8. 8.] | array * vec gives a coo_array (elementwise)
matrix: * is @: True | array: * is @: False
sparse.diags returns dia_matrix | sparse.diags_array returns dia_array
```

Use `@` always. Now the build: $-u''=f$ on $(0,1)$ with $u(0)=u(1)=0$, on the grid $x_j=jh$, $h=1/(n+1)$, $j=1,\dots,n$.

1. Taylor: $u(x\pm h)=u\pm hu'+\tfrac{h^2}{2}u''\pm\tfrac{h^3}{6}u'''+\tfrac{h^4}{24}u''''+\cdots$; add the two: $u(x+h)+u(x-h)=2u+h^2u''+\tfrac{h^4}{12}u''''+O(h^6)$.
2. So $\dfrac{-u_{j-1}+2u_j-u_{j+1}}{h^2}=-u''(x_j)-\dfrac{h^2}{12}u''''(x_j)+\cdots$, and with the zero boundary values the discrete problem is $A_n\mathbf u=\mathbf f$, $\;A_n=\dfrac1{h^2}\,\mathrm{tridiag}(-1,\,2,\,-1)$. $\quad(3.12)$

▶ **Try it first.** Apply $T=\mathrm{tridiag}(-1,2,-1)$ to $v_j=\sin(j\theta)$ and guess the eigenvalue.

**Theorem 3.1 (eigenpairs of $T=\mathrm{tridiag}(-1,2,-1)$).** $T$ ($n\times n$) has eigenvalues $\lambda_k=2-2\cos\frac{k\pi}{n+1}$ with eigenvectors $v_k[j]=\sin\frac{jk\pi}{n+1}$, $k=1,\dots,n$. $\quad(3.13)$

*Proof.* Put $\theta=k\pi/(n+1)$. Row $j$ of $Tv_k$ is $-\sin((j-1)\theta)+2\sin j\theta-\sin((j+1)\theta)$. The addition formulas give $\sin((j\pm1)\theta)=\sin j\theta\cos\theta\pm\cos j\theta\sin\theta$, so the row equals $(2-2\cos\theta)\sin j\theta$. For $j=1$ and $j=n$ the missing neighbours would be $\sin 0=0$ and $\sin((n+1)\theta)=\sin k\pi=0$, so the boundary rows obey the same formula. ∎

3. For $f=\pi^2\sin\pi x$, $\mathbf f=\pi^2v_1$ is an eigenvector, so the discrete solution is exactly $\mathbf u=\dfrac{\pi^2h^2}{2-2\cos\pi h}\,v_1$. With $2-2\cos z=z^2-z^4/12+\cdots$, its error against $u=\sin\pi x$ is the factor $\pi^2h^2/(2-2\cos\pi h)-1\approx\pi^2h^2/12$: second order in $h$, predicted to the digit.

```python
def lap1d(n):
    """CSC form of (3.12): (1/h^2) tridiag(-1, 2, -1), h = 1/(n+1)."""
    h = 1.0 / (n + 1)
    return sparse.diags([-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1], format="csc") / h**2


def best(fn, repeat=3):
    """Smallest wall time of repeat calls, and the last result."""
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        out = fn()
        times.append(time.perf_counter() - t0)
    return min(times), out


spla.spsolve(lap1d(10), np.ones(10))                       # warm-up: first call pays import and setup costs
for n in (1000, 10000):                                    # the dense n=10000 case needs about 1.6 GB of RAM
    h = 1 / (n + 1)
    xg = np.arange(1, n + 1) * h
    f_src = np.pi**2 * np.sin(np.pi * xg)
    A_sp = lap1d(n)
    t_sp, u_sp = best(lambda: spla.spsolve(A_sp, f_src))
    A_dn = A_sp.toarray()
    t_dn, u_dn = best(lambda: np.linalg.solve(A_dn, f_src), repeat=1 if n > 5000 else 3)
    pred = np.pi**2 * h**2 / (4 * np.sin(np.pi * h / 2) ** 2) - 1      # 4 sin^2(z/2) = 2 - 2 cos z without cancellation
    mem_sp = (A_sp.data.nbytes + A_sp.indices.nbytes + A_sp.indptr.nbytes) / 1e6
    print(f"n={n:5d}: sparse {t_sp * 1e3:7.2f} ms {mem_sp:5.2f} MB | dense {t_dn * 1e3:8.1f} ms {A_dn.nbytes / 1e6:6.1f} MB | "
          f"ratio {t_dn / t_sp:5.0f}x | error {np.abs(u_sp - np.sin(np.pi * xg)).max():.3e} vs predicted {pred:.3e}")
n = 10000
w4 = spla.eigsh(lap1d(n), k=4, sigma=0, which="LM", return_eigenvectors=False)           # shift-invert: smallest eigenvalues
lam_exact = 4 * np.sin(np.arange(1, 5) * np.pi / (2 * (n + 1))) ** 2 * (n + 1) ** 2
print("eigsh smallest 4:", np.round(np.sort(w4), 4), "| max relative error vs (3.13):", np.max(np.abs(np.sort(w4) / lam_exact - 1)))
```
```text
n= 1000: sparse    0.28 ms  0.04 MB | dense     14.1 ms    8.0 MB | ratio    51x | error 8.208e-07 vs predicted 8.208e-07
n=10000: sparse    5.16 ms  0.40 MB | dense   3697.9 ms  800.0 MB | ratio   717x | error 8.209e-09 vs predicted 8.223e-09
eigsh smallest 4: [  9.8696  39.4784  88.8264 157.9136] | max relative error vs (3.13): 1.3592904579695642e-11
```

**■ What just happened.** $Tv_k=(2-2\cos\theta)v_k$, and (3.13) then predicts the discrete error $\pi^2h^2/12$ to four digits and the lowest eigenvalues to $10^{-11}$. Add a diagonal $V(x_j)$ and this is QM-3's finite-difference Hamiltonian; `eigsh` with `sigma` near the energy you want returns those states without a dense matrix.


### 3.11 `interpolate`

Interpolation passes through every data point; fitting (§3.7) uses few parameters to average noise away, and interpolating noisy data reproduces the noise. `CubicSpline(x, y)` is a piecewise cubic with continuous first and second derivatives, with `.derivative()` and `.integrate()`. `interp1d` is legacy (no longer developed): prefer `CubicSpline`, the monotone `PchipInterpolator`, or `make_interp_spline`. `splrep`/`splev` are the older FITPACK interface; its smoothing factor `s` makes it a fit. **Runge's phenomenon:** the polynomial through $n$ equispaced points has error $\frac{f^{(n)}(\xi)}{n!}\prod_i(x-x_i)$ (a standard theorem, proved by repeated use of Rolle's theorem in any numerical analysis text; not derived here), and for $f=1/(1+25x^2)$ the derivatives grow faster than $n!$, so more points make it worse:

```python
runge = lambda x: 1 / (1 + 25 * x**2)
xx = np.linspace(-1, 1, 2001)
for n in (5, 13, 21):
    xn = np.linspace(-1, 1, n)
    xc = np.cos(np.pi * (2 * np.arange(n) + 1) / (2 * n))                          # Chebyshev nodes
    errs = [np.abs(p(xx) - runge(xx)).max() for p in (np.polynomial.Polynomial.fit(xn, runge(xn), n - 1),
                                                       np.polynomial.Polynomial.fit(xc, runge(xc), n - 1), interpolate.CubicSpline(xn, runge(xn)))]
    print(f"n={n:2d}  equispaced polynomial {errs[0]:8.3f} | Chebyshev nodes {errs[1]:7.4f} | CubicSpline {errs[2]:.4f}")
```
```text
n= 5  equispaced polynomial    0.438 | Chebyshev nodes  0.4020 | CubicSpline 0.3171
n=13  equispaced polynomial    3.663 | Chebyshev nodes  0.0692 | CubicSpline 0.0069
n=21  equispaced polynomial   59.822 | Chebyshev nodes  0.0153 | CubicSpline 0.0032
```

The cures are splines or Chebyshev nodes; the project draws the $n=13$ figure.

### 3.12 `special`, `stats`, `fft`, `signal`, `constants`

```python
from scipy import constants

x_ = sp.Symbol("x")
leg_sym = sp.lambdify(x_, sp.legendre(3, x_), modules="numpy")
xg_ = np.linspace(-1, 1, 7)
print("erfc(5) =", special.erfc(5.0), "but 1 - erf(5) =", 1 - special.erf(5.0), "| gamma(5) =", special.gamma(5), "| gamma(1/2)^2 =", special.gamma(0.5) ** 2)
print("P3: special vs SymPy-lambdified, max difference", np.abs(special.eval_legendre(3, xg_) - leg_sym(xg_)).max(), "| J0 at its first zero:", special.jv(0, special.jn_zeros(0, 1)[0]))
print("Y_2^1 at polar 0.7, azimuth 0.3:", np.round(special.sph_harm_y(2, 1, 0.7, 0.3), 6), "| SymPy:", sp.N(sp.Ynm(2, 1, 0.7, 0.3).expand(func=True), 6))
print("norm.ppf(0.975) =", stats.norm.ppf(0.975), "| t(9):", stats.t.ppf(0.975, 9), "| Poisson(3) pmf at 2:", stats.poisson.pmf(2, 3.0))
sample = stats.norm.rvs(loc=1.0, scale=2.0, size=500, random_state=rng)
d = stats.describe(sample)
tt_ = stats.ttest_ind(rng.normal(0, 1, 30), rng.normal(0.8, 1, 30))
print(f"norm.fit (MLE): {np.round(stats.norm.fit(sample), 3)} | describe: mean {d.mean:.3f} variance {d.variance:.3f} | two-sample t = {tt_.statistic:.3f}, p = {tt_.pvalue:.1e}, df = {tt_.df:.0f}")
N, fs = 64, 64.0
sig = np.sin(2 * np.pi * 5 * np.arange(N) / fs)                                  # a 5 Hz tone, one second of samples
F, fr = fft.fft(sig), fft.fftfreq(N, d=1 / fs)
kpk = np.argmax(np.abs(F[: N // 2]))
print(f"fft peak at bin {kpk} = {fr[kpk]} Hz, |F[5]| = {np.abs(F[5]):.1f} (= N/2) | fftfreq tail {fr[-3:]} | rfftfreq has {fft.rfftfreq(N, 1 / fs).size} bins")
print("constants.g =", constants.g, "| electron mass:", constants.physical_constants["electron mass"])
```
```text
erfc(5) = 1.5374597944280347e-12 but 1 - erf(5) = 1.5374368445009168e-12 | gamma(5) = 24.0 | gamma(1/2)^2 = 3.1415926535897927
P3: special vs SymPy-lambdified, max difference 2.7755575615628914e-16 | J0 at its first zero: 1.0239231817845155e-16
Y_2^1 at polar 0.7, azimuth 0.3: (-0.363652-0.112491j) | SymPy: -0.363652 - 0.112491*I
norm.ppf(0.975) = 1.959963984540054 | t(9): 2.262157162798205 | Poisson(3) pmf at 2: 0.22404180765538775
norm.fit (MLE): [1.104 1.912] | describe: mean 1.104 variance 3.664 | two-sample t = -5.385, p = 1.4e-06, df = 58
fft peak at bin 5 = 5.0 Hz, |F[5]| = 32.0 (= N/2) | fftfreq tail [-3. -2. -1.] | rfftfreq has 33 bins
constants.g = 9.80665 | electron mass: (9.1093837139e-31, 'kg', 2.8e-40)
```

`special` holds the solutions of the ODEs in EM-3 (`eval_legendre`, `jv`, `jn_zeros`), QM-5 (`sph_harm_y`) and QM-6 (`genlaguerre`). Two traps are visible: `erfc` exists because $1-\mathrm{erf}(5)$ cancels away most of its digits, and `sph_harm_y(n, m, polar, azimuth)` replaced `sph_harm(m, n, azimuth, polar)`, reversed on both counts. `stats` distributions share one API: `pdf` or `pmf`, `cdf`, `ppf` (inverse cdf), `rvs(random_state=rng)`, `fit` (maximum likelihood, ST-1; `norm.fit` divides by $N$, `describe` by $N-1$). `ttest_ind` previews ST-2. `fft.fft` follows F3; `fftfreq(N, d)` gives $k/(Nd)$ with negative frequencies in the second half; use `rfft` for real input. `scipy.signal` (`butter` with `sosfiltfilt`, `find_peaks`, `welch`) is where F3's sampling theory becomes tools. `constants.g` is the standard 9.80665, not our 9.81.


## 4. Worked examples

### Example 4.1: the pendulum's Euler–Lagrange equation, by hand and by SymPy

By hand: $\mathcal L=\tfrac12mL^2\dot\theta^2+mgL\cos\theta$ gives $\partial\mathcal L/\partial\dot\theta=mL^2\dot\theta$, with time derivative $mL^2\ddot\theta$, and $\partial\mathcal L/\partial\theta=-mgL\sin\theta$, so (3.3) reads $mL^2\ddot\theta+mgL\sin\theta=0$. SymPy subtracts the hand result from the `EL` of §3.3 and prints the solved equation:

```python
hand = m * Ls**2 * theta.diff(t, 2) + m * gs * Ls * sp.sin(theta)
print("SymPy EL minus hand result:", sp.simplify(EL - hand))
sp.pprint(sp.Eq(theta.diff(t, 2), acc), use_unicode=False)
```
```text
SymPy EL minus hand result: 0
 2                               
d               -g*sin(theta(t)) 
---(theta(t)) = -----------------
  2                     L        
dt                               
```

### Example 4.2: reading an `OdeResult`

One period with `dense_output=True` and the event of §3.3:

```python
one = integrate.solve_ivp(pend_rhs_sym, (0, T_ex), [th0, 0.0], method="DOP853", rtol=1e-10, atol=1e-12, dense_output=True, events=crossing)
print(f"fields {sorted(one.keys())}\nstatus={one.status} success={one.success} message={one.message!r}")
print(f"t {one.t.shape}, y {one.y.shape}, nfev={one.nfev}, njev={one.njev}, nlu={one.nlu}")
print("t_events:", [np.round(te_, 6).tolist() for te_ in one.t_events], "| y_events:", [np.round(ye, 6).tolist() for ye in one.y_events])
print("state at T/2 from the interpolant one.sol:", np.round(one.sol(T_ex / 2), 8))
```
```text
fields ['message', 'nfev', 'njev', 'nlu', 'sol', 'status', 'success', 't', 't_events', 'y', 'y_events']
status=0 success=True message='The solver successfully reached the end of the integration interval.'
t (39,), y (2, 39), nfev=644, njev=0, nlu=0
t_events: [[0.666468]] | y_events: [[[0.0, -5.271129]]]
state at T/2 from the interpolant one.sol: [-2. -0.]
```

`status` is $0$ if the span was completed, $1$ if a terminal event stopped the run, $-1$ if the step size collapsed (then `success` is `False`). `t` holds the *accepted step times*, not a uniform grid, unless you pass `t_eval`; `y` has shape `(n_states, len(t))`, so a trajectory is a column. `nfev` counts right-hand-side calls; `njev` and `nlu` are used only by implicit methods. `t_events` has one array per event function. `sol` is a callable interpolant over the span: at $T/2$ the pendulum is at $-\theta_0$ with $\omega=0$, as symmetry requires.

### Example 4.3: energy drift against tolerance (Build 3)

Each row is a 100-period run with `atol = rtol`; the last lines fit the slope, show growth in time, and bring in your `rk4`.

```python
rtols = (1e-4, 1e-8, 1e-12)
drifts, drift_nfev = {}, {}
for method in ("RK45", "DOP853"):
    for rtol in rtols:
        s3 = integrate.solve_ivp(pend_rhs_sym, (0, n_per * T_ex), [th0, 0.0], method=method, rtol=rtol, atol=rtol, t_eval=t_eval)
        e3 = spec_energy(s3.y)
        drifts[method, rtol], drift_nfev[method, rtol] = np.max(np.abs(e3 / e3[0] - 1)), s3.nfev
        print(f"{method:7s} rtol={rtol:.0e}  nfev={s3.nfev:7d}  max relative energy drift={drifts[method, rtol]:.2e}")
for method in ("RK45", "DOP853"):
    print(f"{method}: log-log slope of drift against rtol = {np.polyfit(np.log10(rtols), np.log10([drifts[method, r] for r in rtols]), 1)[0]:.2f}")
run_max = np.maximum.accumulate(drift)                              # DOP853, rtol 1e-10, from section 3.3
print("running maximum drift at 25, 50, 75, 100 periods:", [f"{run_max[k]:.1e}" for k in (500, 1000, 1500, 2000)], "| your RK4 at 20000 calls:", f"{rk4_drift:.1e}")
```
```text
RK45    rtol=1e-04  nfev=   8480  max relative energy drift=1.32e-01
RK45    rtol=1e-08  nfev=  43070  max relative energy drift=7.85e-06
RK45    rtol=1e-12  nfev= 259358  max relative energy drift=7.47e-10
DOP853  rtol=1e-04  nfev=  11900  max relative energy drift=3.38e-03
DOP853  rtol=1e-08  nfev=  37436  max relative energy drift=1.92e-07
DOP853  rtol=1e-12  nfev=  88277  max relative energy drift=1.27e-10
RK45: log-log slope of drift against rtol = 1.03
DOP853: log-log slope of drift against rtol = 0.93
running maximum drift at 25, 50, 75, 100 periods: ['2.5e-09', '4.4e-09', '6.5e-09', '8.5e-09'] | your RK4 at 20000 calls: 2.7e-03
```

Three readings. (a) The drift scales like `rtol` to the power one (slopes 1.03 and 0.93): a factor $10^3$ in tolerance buys three digits of energy, as guess (i) of §2 said. (b) The drift grows roughly linearly in time, since a non-symplectic integrator leaks energy steadily: for a long-time claim, choose the tolerance for the *last* period. (c) `DOP853` reaches 1e-10 in 88277 calls where `RK45` needs 259358 for 7e-10; your RK4 (3e-03 at 20 000 calls) has no error control at all, so its drift is whatever $h$ gives.

### Example 4.4: `quad` on a kink, and the error estimate that lies

```python
kink = lambda x: abs(x - 0.3)                                            # exact integral over [0, 1]: 0.29
for pts in (None, [0.3]):
    val, err, info = integrate.quad(kink, 0, 1, points=pts, full_output=True)[:3]
    print(f"|x-0.3|, points={pts}: value={val!r}, abserr={err:.1e}, neval={info['neval']}")
spike = lambda x: np.exp(-(((x - 0.3) / 1e-3) ** 2))                     # width 1e-3 at x = 0.3
print(f"spike, exact integral {1e-3 * np.sqrt(np.pi):.6e}")
for label, call in [("whole range", lambda: integrate.quad(spike, 0, 10)),
                    ("points=[0.3]", lambda: integrate.quad(spike, 0, 10, points=[0.3])),
                    ("three pieces split at .29, .31", lambda: tuple(map(sum, zip(integrate.quad(spike, 0, 0.29), integrate.quad(spike, 0.29, 0.31), integrate.quad(spike, 0.31, 10)))))]:
    val, err = call()
    print(f"{label}: value={val:.6e}, abserr={err:.1e}")
```
```text
|x-0.3|, points=None: value=0.2899999999999998, abserr=3.2e-16, neval=357
|x-0.3|, points=[0.3]: value=0.29, abserr=3.2e-15, neval=42
spike, exact integral 1.772454e-03
whole range: value=0.000000e+00, abserr=0.0e+00
points=[0.3]: value=8.862269e-04, abserr=7.8e-10
three pieces split at .29, .31: value=1.772454e-03, abserr=1.8e-09
```

A kink is easy: bisection zooms in either way, and `points=[0.3]` only saves calls. A feature narrower than the spacing of the first 21 nodes is invisible: both rules see zeros, agree, and the estimate says `0.0`. Declaring a breakpoint is not enough (one half of the spike is found, the other is not, and the estimate is still tiny). Integrate a window around features you know about and sum the pieces, and treat `abserr` as an estimate that assumes the function is typical.

### Example 4.5: `curve_fit` into $g\pm\sigma$, and whether to believe it

The exact-model fit of Build 5 licenses one sentence: "$g=$ 9.81 $\pm$ 0.03 m/s² (1σ), assuming the model (3.6), an exactly known amplitude, and independent timing noise of 0.01 s." Test the claim by repeating the experiment 500 times:

```python
def coverage(model, reps=500, seed=1):
    """Fraction of repeated experiments whose fitted g lies within 1 and 2 reported sigmas of the truth."""
    r = np.random.default_rng(seed)
    z = []
    for _ in range(reps):
        Tm = period(th_amp, Ls_m, g_true) + sig_T * r.standard_normal(Ls_m.size)
        p, cv = optimize.curve_fit(model, Ls_m, Tm, p0=[8.0], sigma=np.full(Ls_m.size, sig_T), absolute_sigma=True)
        z.append(abs(p[0] - g_true) / np.sqrt(cv[0, 0]))
    z = np.array(z)
    return (z < 1).mean(), (z < 2).mean()


for name, model in models.items():
    c1, c2 = coverage(model)
    print(f"{name:12s} within 1 sigma: {c1:.3f}   within 2 sigma: {c2:.3f}")
```
```text
small-angle  within 1 sigma: 0.000   within 2 sigma: 0.000
exact (3.6)  within 1 sigma: 0.676   within 2 sigma: 0.956
```

With the right model the fractions are the Gaussian 0.683 and 0.954 to sampling error; with the small-angle model they are zero. The covariance is only as honest as the model.

### Example 4.6: the 6×6 Laplacian and its eigenvalues

```python
T6 = sparse.diags([-np.ones(5), 2 * np.ones(6), -np.ones(5)], [-1, 0, 1], format="csr")
print(T6.toarray().astype(int))
print("stored", T6.nnz, "of 36 | CSR data", T6.data.astype(int).tolist(), "| indices", T6.indices.tolist(), "| indptr", T6.indptr.tolist())
ev = np.linalg.eigvalsh(T6.toarray())
an = 2 - 2 * np.cos(np.arange(1, 7) * np.pi / 7)
v1 = np.sin(np.arange(1, 7) * np.pi / 7)
print("eigenvalues :", np.round(ev, 6), "\nanalytic (3.13):", np.round(an, 6), "\nmax difference", np.abs(ev - an).max(), "| T v1 - lambda1 v1:", np.abs(T6 @ v1 - an[0] * v1).max())
```
```text
[[ 2 -1  0  0  0  0]
 [-1  2 -1  0  0  0]
 [ 0 -1  2 -1  0  0]
 [ 0  0 -1  2 -1  0]
 [ 0  0  0 -1  2 -1]
 [ 0  0  0  0 -1  2]]
stored 16 of 36 | CSR data [2, -1, -1, 2, -1, -1, 2, -1, -1, 2, -1, -1, 2, -1, -1, 2] | indices [0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5] | indptr [0, 2, 5, 8, 11, 14, 16]
eigenvalues : [0.198062 0.75302  1.554958 2.445042 3.24698  3.801938] 
analytic (3.13): [0.198062 0.75302  1.554958 2.445042 3.24698  3.801938] 
max difference 4.440892098500626e-16 | T v1 - lambda1 v1: 1.6653345369377348e-16
```

`indptr[i]:indptr[i+1]` slices row $i$'s values and column numbers out of `data` and `indices`: that is all CSR is. The dense eigenvalues and (3.13) agree to rounding.


## 5. Why it matters and where it connects

- **CM-3:** the §3.3 pipeline is that module's workflow for any system (Exercise B1: a bead on a hoop). **CM-4:** (3.4) is the energy function whose conservation Noether's theorem explains; Example 4.3 is its numerical test. **CM-5:** events find apsides. **CM-6:** `eigh` (§3.9) solves the normal-mode eigenproblem. **CM-8:** Hamilton's equations go through `solve_ivp` too, where a symplectic method would beat Example 4.3's drift.
- **EM-3:** `eval_legendre` and `jv` are its separated solutions; §3.10's Laplacian in two dimensions is the grid solver for $\nabla^2\phi=-\rho/\epsilon_0$. **QM-3:** §3.10 plus a diagonal $V$ and `eigsh` (Exercise B5). **QM-5, QM-6:** `sph_harm_y`, `genlaguerre`.
- **ST-1:** `minimize` of a negative log-likelihood is maximum likelihood, and (3.10) is its large-sample covariance (inverse Fisher information for Gaussian noise). **ST-2:** `ttest_ind`. **ST-3:** (3.10) is the OLS covariance. **LA-7, LA-11:** `expm` is the exponential they derive. **LA-13:** `eigsh` is a Krylov method; conditioning explains §3.10's rounding. **F3:** `fftfreq` follows its conventions box. **PY-10:** profiling finds the `subs` or loop hiding in a pipeline like this.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are predict-and-fix drills (SymPy gotchas, result objects, substitution order, `atol`, sparse `*`, argument orders). B1–B7 are the module; B1 (bead on a rotating hoop through the whole pipeline), B3 (damped pendulum: events, `dsolve`), B4 (covariance by hand) and B7 (extend `project.py`) are essential. C1–C3: a resonance curve from events, a boundary-value problem with two solutions plus Rosenbrock, and the measured cost of the symbolic pipeline for an $N$-link pendulum.

## 7. Coding project

`project.py` runs the pipeline with asserts, headless, in about 30 s (set `DENSE_BIG = False` if you have under 4 GB free): (1) the pendulum: Euler–Lagrange equal to the hand result, energy conserved to $10^{-8}$ at `rtol=1e-10`, period from events against (3.6) to $10^{-6}$, the small-angle limit; (2) the tolerance sweep and RK4 as a log-log figure; (3) stiffness (`nfev` ratio above 10 on two problems); (4) the double pendulum against the published form, with its divergence plot; (5) `curve_fit` recovering a planted $g$ within $2\sigma$, (3.10) from a SymPy Jacobian, a coverage experiment; (6) the sparse Poisson solve against the analytic solution, timing reported; (7) `expm` against `solve_ivp`; (8) `special`, `stats`, `fft`, `interpolate`, `lambdify` checks.

**Build it yourself:** before opening the file, write Part 1's seven stages for the bead on a hoop. Then "Your turn": damping and driving with a resonance curve; Rosenbrock with and without the analytic gradient; a boundary-value problem with `solve_bvp`; `dsolve` of the damped oscillator against `solve_ivp`.

## 8. Where you'll struggle

- **Default tolerances are for pictures.** §1's 16% drift is the default; every student plots a spiralling "conserved" energy once. Set `rtol` and `atol` and measure the invariant; tolerance buys energy digits at slope one, no faster (Example 4.3).
- **Failure is quiet.** `success=False`, a `message`, an `IntegrationWarning`, a `ConditionSet`, an `abserr` of `0.0` on a spike: none stops your script. Read `.success` and `nfev`, and test one case with a known answer (a period from (3.6), an eigenvalue from (3.13)) before trusting a pipeline.
- **SymPy traps.** `1/3` is a float; substitute the highest derivative first; `solve` returns a list, `solveset` a set, and either may return nothing or raise for a solvable equation; a symbol without `real=True` may be complex; SymPy is slow, so never in an inner loop: `lambdify` first.
- **Argument orders.** `ellipk(m)` takes $m=k^2$; `dblquad` wants `func(y, x)`; `sph_harm_y(n, m, polar, azimuth)`; `fftfreq(N, d)` is in cycles per unit of $d$. Check each on a case you know.
- **Sparse is not dense.** `*` is matrix multiply for `csr_matrix` and elementwise for `csr_array`: use `@`. `toarray()` at $n=10^4$ costs 800 MB. Build in COO or `diags`, solve in CSC.
- **`pcov` is a promise about the model.** It assumes the model is right and the noise independent with the stated $\sigma_i$; the small-angle fit printed $\pm0.03$ and was wrong by 0.3. Say so when you quote an uncertainty.
- **The hard part is not the API.** It is deciding which stage failed when the plot looks wrong: the equation, the compile, the integration or the model. Keep each stage checkable alone, as §3.3 does.

## 9. Check yourself

1. Why does default `solve_ivp` lose 16% of a pendulum's energy in 20 swings yet report `success=True`? What does (3.8) bound, and what not?
2. State the rule for `atol` and give a component for which the default is wrong.
3. Derive the exact pendulum period (3.6) in outline. What does `special.ellipk` take as argument?
4. In what order do you substitute `theta(t)` and its derivatives by plain symbols, and what goes wrong otherwise?
5. Why must an explicit method take steps below $3.3/\lambda$ on $y'=-\lambda(y-\cos t)$? What does backward Euler do instead?
6. State the assumptions behind $\sqrt{\mathrm{pcov}_{ii}}$ and derive (3.10) in three lines.
7. What does `quad` return, and how can its second value mislead?
8. Which sparse format for assembly, for `A @ v`, for `spsolve`? Why is `A * B` risky?
9. State the eigenvalues of $\mathrm{tridiag}(-1,2,-1)$ and the error they predict for $-u''=\pi^2\sin\pi x$.
10. Why does a lambdified function beat `subs` by one to two orders of magnitude per call, and what can go wrong with `modules` and with constants?

## 10. Pair with

- **SciPy user guide tutorials for `optimize`, `integrate`, `linalg`, `sparse`, `interpolate`, `special`, `stats`, `fft` (docs.scipy.org).** The authoritative list of options; read `integrate` and `optimize` after §3.4 and §3.7.
- **Johansson, *Numerical Python* 3e, ch. 3, 5, 6, 8, 9.** More worked SymPy and SciPy examples, with plots.
- **The SymPy tutorial ("Introductory Tutorial" and "Gotchas", docs.sympy.org).** Read "Gotchas" after §3.2.
- **SciPy Lecture Notes §1.5, §1.6.** A compact tour.
- **Virtanen et al., "SciPy 1.0", *Nature Methods* 17 (2020); Meurer et al., "SymPy: symbolic computing in Python", *PeerJ CS* 3 (2017).** Reference only.
