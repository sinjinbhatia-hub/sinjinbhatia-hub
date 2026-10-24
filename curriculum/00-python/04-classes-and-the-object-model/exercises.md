# PY-4 Exercises

Conventions as in the notes: Python 3.11+, `Vec`, `Rational`, `Dual`, `dmath`, `Matrix` as in `project.py` (put your answers in a file next to it and `from project import ...`; importing runs its checks, about two seconds). Write your prediction down before you run anything. "Write a class" means: write it, then test it with `assert`s you chose, including an edge case.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Class and instance attributes; bound methods.) Predict the four printed lines, then repair the class so that each account has its own history.

```python
class Account:
    rate = 0.02
    history = []

    def __init__(self, balance):
        self.balance = balance

    def deposit(self, x):
        self.balance += x
        self.history.append(x)


a, b = Account(100), Account(50)
a.rate = 0.05
a.deposit(10)
b.deposit(5)
print(a.rate, b.rate, Account.rate)
print(a.history, b.history)
f = a.deposit
f(1)
print(Account.deposit(b, 2), b.balance, a.balance)
```

**A2.** (`__eq__` and `__hash__`.) Which of the four lines raise, with what error, and why? Then fix `Money` so that `m == 5` is `False`, `{Money(5), Money(5)}` has one element, and `Money` objects work as dict keys.

```python
class Money:
    def __init__(self, cents):
        self.cents = cents

    def __eq__(self, other):
        return self.cents == other.cents


m = Money(5)
print(m == Money(5))
print(m != Money(6))
print(m == 5)
print({m})
```

**A3.** (Reflected operators.) With the `Vec` of the notes, `sum([Vec(1, 2), Vec(3, 4)])` raises `TypeError`. Say why, using the first step of `sum`. Fix it two ways: with an argument to `sum`, and with a `__radd__` on a subclass. Which fix would you rather ship? Then: a class defines `__iadd__` that mutates `self` but has no `return`. What is `x` after `x += y`?

**A4.** (Write a class, 5–15 lines.) `Temperature` stores kelvin. Properties: `kelvin` (setter rejects negative values with `ValueError`), `celsius` (get and set, in terms of `kelvin`). A classmethod `from_fahrenheit`. A `__repr__` that round-trips. Ordering by kelvin with `functools.total_ordering`, and `__eq__`. Then decide whether `__hash__` should exist, given the setters, and defend the decision with a failing example. Assert all of it.

**A5.** (`Dual` by hand.) Using only (3.1), (3.2) and (3.3), evaluate with `Dual(2, 1)` by hand: (a) $f(x) = (x^2 + 1)/x$, giving $f(2)$ and $f'(2)$; (b) $g(x) = x^3\cdot(1/x)$ at $x = 3$, and say what calculus fact the $\varepsilon$ part confirms. Check both with `derivative`. Why does `Dual(2.0, 1.0) ** 0.5` raise `TypeError`, and what does the message say?

**A6.** (Dataclasses.) Predict the output, say why `hash(v1)` raises, and change the class (two changes) so that versions are hashable and `tags` do not take part in comparisons.

```python
from dataclasses import dataclass, field


@dataclass(order=True)
class Version:
    major: int
    minor: int = 0
    tags: list = field(default_factory=list)


v1, v2 = Version(1, 2), Version(1, 10)
print(v1, v1 < v2, v1 == Version(1, 2, []))
print(hash(v1))
```

## B. Core (do all; these are the module)

**B1.** (`Vec` indexing; extends `project.py`.) Add `__getitem__` to `Vec`: a float for an `int` (negative indices allowed), a `Vec` for a `slice`, `IndexError` out of range, `TypeError` for other keys. Assert `v[1:]`, `v[::-1]`, `list(reversed(v))`, and that `np.asarray(v).shape` becomes `(3,)`. Then explain why `v[0] = 1.0` still fails and what it would take to allow it.

**B2.** (`Complex`.) Write `Complex(re, im)` with `+ - * /`, `abs`, `conjugate`, `==`, `hash`, `repr`, and the reflected operators so that `2 * z`, `1 + z` and `1 / z` work. Compare with the built-in `complex` on 1000 random pairs (`math.isclose`). Which comparison operators does the built-in *not* have, and does yours?

**B3.** (Real and dual exponents; more `dmath`.) Extend `Dual`: `__pow__` for a float exponent (base $a > 0$), for a `Dual` exponent (so `x ** x` works), and `__rpow__` (`2 ** x`). Add `dmath.tan` and `dmath.sinh`. Verify six functions against SymPy at five points each. What should `(-8.0) ** (1/3)` do for a `Dual` base, and what does your code do?

**B4.** (Second derivatives by nesting.) Write `second_derivative(f, x)` by evaluating $f$ at `Dual(Dual(x, 1), Dual(1, 0))`. (a) Expand the seed as $x + \varepsilon_1 + \varepsilon_2$ and say what $\varepsilon_1\varepsilon_2$ is, so that you can say which entry of the result holds $f''$. (b) Check $x e^x$, $\sin(x^2)$ and $1/(1 + x^2)$ against SymPy. (c) Seed `x` with $\varepsilon_1$ (the inner `Dual`) and `y` with $\varepsilon_2$ (the outer) to get the mixed partial $\partial^2 f/\partial x\,\partial y$ of $f(x, y) = x^2y + \sin(xy)$, and check it. This is a Hessian entry (LA-14).

**B5.** (`Matrix` extensions.) Add `det()` from the LU factorization (the product of the diagonal of $U$ times the sign of the permutation), `inverse()` by solving for the columns of the identity, scalar `*` and `__rmul__`, and `__rmatmul__` so that `Vec @ Matrix` is a row vector times a matrix. Compare all four with NumPy on random $5\times5$ matrices, and test the singular case.

**B6.** (`Polynomial` arithmetic.) Add `+`, `-`, `*`, scalar versions with reflected twins, and `from_roots(*roots)` (a classmethod). Compare with `numpy.polynomial.Polynomial` on 50 random pairs. Then evaluate `p(q)` for a *polynomial* `q` and compare with NumPy's composition: explain why Horner's rule composes polynomials without any new code.

**B7.** (Extend `project.py`.) Add an eleventh function to `FUNCS` in Part 3, built from two different `dmath` functions, and verify it like the others. Then draw its error-versus-$h$ curve in the style of `error_vs_h.png` and print the measured best step for forward and central differences. Predict those steps first. For $f(x+h) - f(x)$ the rounding error is about `eps` $\cdot|f|$ and the truncation error $\tfrac12|f''|h$ (forward) or $\tfrac16|f'''|h^2$ (central); minimize the sum over $h$ and compare your formula with the measured minima.

## C. Challenge (do at least one)

**C1.** (Reverse mode.) Write a class `Var` holding a value, a gradient slot, and its parents with their local partial derivatives, with `+ - * /`, `sin`, `exp`, and a `backward()` that propagates in reverse topological order. Check the gradient of $f(x, y) = xy + \sin x$ against two forward `Dual` passes and SymPy. Then take $f(\mathbf x) = \sum_{i=1}^{n} x_i^2\sin x_i$ with $n = 1000$ and time the full gradient by reverse mode (one backward pass) and by forward mode ($n$ passes). Report the ratio and relate it to the cost model in §3.7 and to ML-7.

**C2.** (Design with measurements.) The `Integrator` of §3.10 and the `rule(ys, h)` function both assume a uniform grid. Implement Gauss–Legendre quadrature (nodes from `numpy.polynomial.legendre.leggauss`) in the class design and in the function design and say what had to change in each. Redesign so that a rule is `rule(f, a, b, n)`. Count evaluations of $f$ with a callable counting object (`__call__` with state, §3.4). For $\int_0^\pi \sin x\,dx$ and $\int_0^1 \sqrt x\,dx$ tabulate error against evaluations for trapezoid, Simpson and Gauss–Legendre, and explain what the second integrand (not smooth at 0) does to each rule's rate of convergence. Write half a page on which design you would keep and why.

**C3.** (What does pure Python cost?) Time `Matrix @ Matrix` for $n = 25, 50, 100, 200$ and fit the log-log slope. Compare with NumPy at each $n$. Then speed up the inner loop three ways (for example `sum(map(operator.mul, row, col))`, `math.sumprod` on Python 3.12+, and a row-update (axpy) form that needs no transpose) and report each speedup honestly. Explain, in terms of what the interpreter does per scalar operation, why the remaining gap to NumPy is so large. This is the argument for PY-6 and PY-10.
