# PY-4 Solutions

Every code block was run, and the ones with `assert`s passed them; the outputs shown are the real ones (the timings in C1 and C3 vary by machine: trust the ratios and slopes, not the milliseconds). Solutions assume the classes of the notes: put your file next to `project.py` and begin it with `import math, numbers, operator, functools, time`, `import numpy as np`, `import sympy as sp` and `from project import Vec, Rational, Dual, dmath, Matrix, Polynomial, derivative, _elementary` (importing re-runs the project's checks, about two seconds). Where a solution changes a class from `project.py` it does so on a subclass or temporarily, so that earlier solutions are unaffected; in your own file, put the new methods in the class.

## A. Warm-up

**A1.** `a.rate = 0.05` writes an *instance* attribute (rule 4, §3.1), so `a` sees `0.05`, `b` still sees the class's `0.02`. `history` is a list on the class, and `self.history.append` mutates it: both accounts share it. `a.deposit` is a bound method; `f(1)` is `Account.deposit(a, 1)`. `Account.deposit(b, 2)` returns `None` (no `return`).

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


class FixedAccount:
    rate = 0.02                        # a constant shared by design: fine on the class

    def __init__(self, balance):
        self.balance = balance
        self.history = []              # per-instance state is created in __init__

    def deposit(self, x):
        self.balance += x
        self.history.append(x)


c, d = FixedAccount(100), FixedAccount(50)
c.deposit(10)
d.deposit(5)
assert c.history == [10] and d.history == [5]
```
```text
0.05 0.02 0.02
[10, 5] [10, 5]
None 57 111
```

**A2.** Line 1 prints `True`. Line 2 prints `True`: `!=` is derived from `__eq__`. Line 3 raises `AttributeError`: `m == 5` calls `Money.__eq__(m, 5)`, which reads `5.cents`; the method should have declined with `NotImplemented`. Line 4 raises `TypeError: unhashable type: 'Money'`, because defining `__eq__` set `__hash__` to `None`.

```python
class Money:
    def __init__(self, cents):
        self.cents = cents

    def __eq__(self, other):
        return self.cents == other.cents


m = Money(5)
print(m == Money(5), m != Money(6))
for action in (lambda: m == 5, lambda: {m}):
    try:
        action()
    except (AttributeError, TypeError) as err:
        print(type(err).__name__ + ":", err)


class FixedMoney:
    def __init__(self, cents):
        self._cents = cents

    cents = property(lambda self: self._cents)               # read-only, so a value hash is safe

    def __repr__(self):
        return f"FixedMoney({self._cents!r})"

    def __eq__(self, other):
        if not isinstance(other, FixedMoney):
            return NotImplemented
        return self._cents == other._cents

    def __hash__(self):
        return hash(self._cents)


assert (FixedMoney(5) == 5) is False and len({FixedMoney(5), FixedMoney(5)}) == 1
assert {FixedMoney(5): "x"}[FixedMoney(5)] == "x"
```
```text
True True
AttributeError: 'int' object has no attribute 'cents'
TypeError: unhashable type: 'Money'
```

`m == 5` is now `False`: both sides declined, so identity decided. The hash is safe only because `cents` cannot change after construction (a frozen dataclass would do the same job in four lines).

**A3.** `sum` starts from the integer `0` and computes `0 + Vec(1, 2)`. `int.__add__` declines, and `Vec` has no `__radd__`, so: `unsupported operand type(s) for +: 'int' and 'Vec'`.

```python
v1, v2, v3 = Vec(1, 2), Vec(3, 4), Vec(5, 6)
try:
    sum([v1, v2])
except TypeError as err:
    print("TypeError:", err)
print(sum([v1, v2, v3], start=Vec(0, 0)))


class SummableVec(Vec):
    def __radd__(self, other):
        print("  __radd__ asked; left operand is", repr(other))
        return self if other == 0 else NotImplemented       # only for sum()'s starting 0


total = sum([SummableVec(1, 2), SummableVec(3, 4), SummableVec(5, 6)])
assert total == Vec(9, 12)


class Broken:
    def __init__(self, n):
        self.n = n

    def __iadd__(self, other):
        self.n += other.n                                     # mutates, but forgets: return self


x, y = Broken(1), Broken(2)
x += y
print(x)
```
```text
TypeError: unsupported operand type(s) for +: 'int' and 'Vec'
(9, 12)
  __radd__ asked; left operand is 0
  __radd__ asked; left operand is Vec(4.0, 6.0)
None
```

The second fix works, and the trace shows both rules of §3.3. The first call is `0 + v`: `int.__add__` declines, so the reflected `__radd__` runs (rule 1). After `v1 + v2` the left operand is a plain `Vec` and the right is a `SummableVec`, a subclass that defines `__radd__` where `Vec` has none, so the subclass is asked *first* (rule 2); it declines (the left operand is not `0`) and Python falls back to `Vec.__add__`. I would ship the first fix: `start=` is explicit, whereas `__radd__` teaches `Vec` that `0 + v` is legal, which is false mathematically and only exists to please `sum`. `x += y` is `x = x.__iadd__(y)`, so a missing `return self` leaves `x` bound to `None`.

**A4.**

```python
import functools


@functools.total_ordering
class Temperature:
    def __init__(self, kelvin):
        self.kelvin = kelvin                                  # goes through the setter

    @property
    def kelvin(self):
        return self._k

    @kelvin.setter
    def kelvin(self, value):
        if value < 0:
            raise ValueError(f"below absolute zero: {value} K")
        self._k = float(value)

    @property
    def celsius(self):
        return self._k - 273.15

    @celsius.setter
    def celsius(self, value):
        self.kelvin = value + 273.15

    @classmethod
    def from_fahrenheit(cls, f):
        return cls((f - 32) * 5 / 9 + 273.15)

    def __repr__(self):
        return f"Temperature({self._k!r})"

    def __eq__(self, other):
        if not isinstance(other, Temperature):
            return NotImplemented
        return self._k == other._k

    def __lt__(self, other):
        if not isinstance(other, Temperature):
            return NotImplemented
        return self._k < other._k


t = Temperature(300)
t.celsius = 0
assert abs(t.kelvin - 273.15) < 1e-12 and eval(repr(t)) == t
assert Temperature.from_fahrenheit(32) == Temperature(273.15)
assert Temperature(250) < Temperature(251) <= Temperature(251)
assert sorted([Temperature(3), Temperature(1), Temperature(2)])[0] == Temperature(1)
try:
    t.kelvin = -1
except ValueError as err:
    print("ValueError:", err)
assert Temperature.__hash__ is None
```
```text
ValueError: below absolute zero: -1 K
```

`__eq__` without `__hash__` made the class unhashable, and that is the right outcome here: the setters make a `Temperature` mutable, and a value-based hash of a mutable object breaks sets (a stored key whose hash changes is never found again; compare PY-2 §3.4). To allow hashing you would remove the setters (or use a frozen dataclass); the decision follows from the mutability, not from taste.

**A5.** (a) $x = (2, 1)$. $x^2 = (4,\ 2\cdot2\cdot1) = (4, 4)$ by (3.1); adding $1$ gives $(5, 4)$; dividing by $x$ with (3.2): value $5/2$, slope $(4\cdot2 - 5\cdot1)/2^2 = 3/4$. Calculus: $f = x + 1/x$, $f' = 1 - 1/x^2 = 3/4$ at $2$. (b) $x^3 = (27, 27)$, $1/x = (1/3, -1/9)$, and their product is $(27\cdot\tfrac13,\ 27\cdot(-\tfrac19) + 27\cdot\tfrac13) = (9, 6)$: that is $x^2$ and $2x$ at $3$, the power rule recovered from the product and quotient rules.

```python
print(derivative(lambda x: (x * x + 1) / x, 2.0), derivative(lambda x: x ** 3 * (1 / x), 3.0))
assert derivative(lambda x: (x * x + 1) / x, 2.0) == (2.5, 0.75)
assert derivative(lambda x: x ** 3 * (1 / x), 3.0) == (9.0, 6.0)
try:
    Dual(2.0, 1.0) ** 0.5
except TypeError as err:
    print("TypeError:", err)
```
```text
(2.5, 0.75) (9.0, 6.0)
TypeError: unsupported operand type(s) for ** or pow(): 'Dual' and 'float'
```

`Dual.__pow__` returns `NotImplemented` for a non-`int` exponent, and `float` has no `__rpow__` that accepts a `Dual`, so both operands decline. Exercise B3 supplies the missing case.

**A6.** `order=True` and `eq=True` (the default) generate the methods; `v1 < v2` compares `(1, 2, [])` with `(1, 10, [])`: $2 < 10$, so `True`. Equality compares all three fields. `hash(v1)` raises because `eq=True` without `frozen=True` sets `__hash__` to `None`. The two changes: `frozen=True`, and `compare=False` on `tags`. The second removes the field from `__eq__`, from the order methods *and* from the generated `__hash__` (which hashes the compared fields), so even a list-valued `tags` no longer blocks hashing. A tuple default is still the better record: a frozen dataclass holding a list is frozen only at the surface (`f1.tags.append(...)` still works), and `default=()` needs no factory. Both versions are below.

```python
from dataclasses import dataclass, field


@dataclass(order=True)
class Version:
    major: int
    minor: int = 0
    tags: list = field(default_factory=list)


v1, v2 = Version(1, 2), Version(1, 10)
print(v1, v1 < v2, v1 == Version(1, 2, []))
try:
    hash(v1)
except TypeError as err:
    print("TypeError:", err)


@dataclass(order=True, frozen=True)
class FrozenVersion:                                            # the two changes
    major: int
    minor: int = 0
    tags: list = field(default_factory=list, compare=False)


@dataclass(order=True, frozen=True)
class TupleVersion:                                             # the same, with an immutable default
    major: int
    minor: int = 0
    tags: tuple = field(default=(), compare=False)


f1, f2 = FrozenVersion(1, 2, ["beta"]), FrozenVersion(1, 10)
assert f1 < f2 and f1 == FrozenVersion(1, 2) and hash(f1) == hash(FrozenVersion(1, 2))
assert sorted([f2, f1]) == [f1, f2] and len({f1, FrozenVersion(1, 2)}) == 1
f1.tags.append("rc")                                            # frozen at the surface only
t1 = TupleVersion(1, 2, ("beta",))
assert t1 == TupleVersion(1, 2) and hash(t1) == hash(TupleVersion(1, 2)) and len({t1, TupleVersion(1, 2)}) == 1
```
```text
Version(major=1, minor=2, tags=[]) True True
TypeError: unhashable type: 'Version'
```

## B. Core

**B1.** The method belongs in `Vec`; here it lives on a subclass so that the project's `Vec` stays as the notes define it. `operator.index` accepts ints (and NumPy ints) and raises `TypeError` for floats and strings; indexing the tuple raises `IndexError` out of range and handles negatives; a slice returns the same kind of vector.

```python
class SliceVec(Vec):
    def __getitem__(self, key):
        if isinstance(key, slice):
            return type(self)(*self._c[key])
        return self._c[operator.index(key)]


v = SliceVec(1, 2, 3)
assert v[0] == 1.0 and v[-1] == 3.0 and v[1:] == Vec(2, 3) and v[::-1] == Vec(3, 2, 1)
assert isinstance(v[1:], SliceVec) and list(reversed(v)) == [3.0, 2.0, 1.0]
for bad in (lambda: v[3], lambda: v["x"], lambda: v[1.0]):
    try:
        bad()
    except (IndexError, TypeError) as err:
        print(type(err).__name__ + ":", err)
print(np.asarray(Vec(1, 2, 3)).shape, np.asarray(v).shape, np.asarray(v).dtype)
try:
    v[0] = 5.0
except TypeError as err:
    print("TypeError:", err)
```
```text
IndexError: tuple index out of range
TypeError: 'str' object cannot be interpreted as an integer
TypeError: 'float' object cannot be interpreted as an integer
() (3,) float64
TypeError: 'SliceVec' object does not support item assignment
```

`np.asarray` now sees a sequence (`__len__` and `__getitem__`), so the shape is `(3,)` and the dtype `float64`, where plain `Vec` gave a 0-dimensional object array (§3.10). `reversed(v)` works from `__len__` and `__getitem__` alone. `v[0] = 1.0` fails because there is no `__setitem__`, and the storage is a tuple. To allow it: keep the components in a list, define `__setitem__`, and set `__hash__ = None`. A mutable `Vec` must not be hashable (the contract of §3.2), and the immutability that made `Vec` safe in sets would be gone.

**B2.** Division follows from multiplying numerator and denominator by the conjugate: $\dfrac{a + bi}{c + di} = \dfrac{(a + bi)(c - di)}{c^2 + d^2} = \dfrac{(ac + bd) + (bc - ad)i}{c^2 + d^2}$. The hash is delegated to the built-in `complex`, so that `Complex(2, 0) == 2` and the hashes agree (the contract again).

```python
import random


def _lift(x):
    if isinstance(x, Complex):
        return x
    return Complex(x, 0.0) if isinstance(x, numbers.Real) else None


def _binary(method):
    @functools.wraps(method)
    def wrapper(self, other):
        other = _lift(other)
        return NotImplemented if other is None else method(self, other)
    return wrapper


class Complex:
    def __init__(self, re, im=0.0):
        self.re, self.im = float(re), float(im)

    def __repr__(self):
        return f"Complex({self.re!r}, {self.im!r})"

    @_binary
    def __eq__(self, o):
        return (self.re, self.im) == (o.re, o.im)

    def __hash__(self):
        return hash(complex(self.re, self.im))

    def __neg__(self):
        return Complex(-self.re, -self.im)

    def __abs__(self):
        return math.hypot(self.re, self.im)

    def conjugate(self):
        return Complex(self.re, -self.im)

    @_binary
    def __add__(self, o):
        return Complex(self.re + o.re, self.im + o.im)

    __radd__ = __add__

    @_binary
    def __sub__(self, o):
        return self + (-o)

    @_binary
    def __rsub__(self, o):
        return o + (-self)

    @_binary
    def __mul__(self, o):
        return Complex(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    __rmul__ = __mul__

    @_binary
    def __truediv__(self, o):
        n = o.re ** 2 + o.im ** 2
        if n == 0.0:
            raise ZeroDivisionError("division by zero")
        num = self * o.conjugate()
        return Complex(num.re / n, num.im / n)

    @_binary
    def __rtruediv__(self, o):
        return o / self


rng = random.Random(0)


def close(z, c):
    return all(math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12) for a, b in ((z.re, c.real), (z.im, c.imag)))


for _ in range(1000):
    c1, c2 = (complex(rng.uniform(-5, 5), rng.uniform(-5, 5)) for _ in range(2))
    z1, z2 = Complex(c1.real, c1.imag), Complex(c2.real, c2.imag)
    assert close(z1 + z2, c1 + c2) and close(z1 - z2, c1 - c2) and close(z1 * z2, c1 * c2) and close(z1 / z2, c1 / c2)
    assert close(2 * z1, 2 * c1) and close(1 + z1, 1 + c1) and close(1 / z1, 1 / c1)
    assert close(z1.conjugate(), c1.conjugate()) and math.isclose(abs(z1), abs(c1)) and (z1 == z2) == (c1 == c2)
assert Complex(2, 0) == 2 and hash(Complex(2, 0)) == hash(2) and Complex(1, 2) != 1.5
for ordered in (lambda: sorted([Complex(1, 2), Complex(0, 1)]), lambda: sorted([1 + 2j, 1j])):
    try:
        ordered()
    except TypeError as err:
        print("TypeError:", err)
```
```text
TypeError: '<' not supported between instances of 'Complex' and 'Complex'
TypeError: '<' not supported between instances of 'complex' and 'complex'
```

Neither has `<`, `<=`, `>` or `>=`: the complex numbers are not an ordered field, so a total order compatible with arithmetic does not exist, and `sorted` fails for both with `TypeError`. `Complex` gets this for free by not defining `__lt__`.

**B3.** For a real exponent, (3.3) with $f(t) = t^p$ gives $(a + b\varepsilon)^p = a^p + p\,a^{p-1}b\,\varepsilon$, valid for $a > 0$. For a dual exponent use $x^y = e^{y\ln x}$ and let `dmath.exp` and `dmath.log` apply the chain rule; then $\frac{d}{dx}x^x = x^x(\ln x + 1)$ comes out without a new rule. `__rpow__` handles `base ** dual` by the same identity, $b^x = e^{x\ln b}$. The two methods are installed temporarily on `Dual` (put them in the class in your own file).

```python
dmath.tan = _elementary(math.tan, lambda a: 1 + dmath.tan(a) ** 2)
dmath.sinh = _elementary(math.sinh, lambda a: dmath.cosh(a))
dmath.cosh = _elementary(math.cosh, lambda a: dmath.sinh(a))

_int_pow = Dual.__pow__


def general_pow(self, p):
    if isinstance(p, int):
        return _int_pow(self, p)
    if isinstance(p, numbers.Real):
        if self.a <= 0:
            raise ValueError("a real power needs a positive base")
        return Dual(self.a ** p, p * self.a ** (p - 1) * self.b)
    if isinstance(p, Dual):
        return dmath.exp(p * dmath.log(self))
    return NotImplemented


def general_rpow(self, base):
    if isinstance(base, numbers.Real) and base > 0:
        return dmath.exp(self * math.log(base))
    return NotImplemented


Dual.__pow__, Dual.__rpow__ = general_pow, general_rpow
try:
    FUNCS = {
        "x**1.5": lambda m, x: x ** 1.5,
        "x**x": lambda m, x: x ** x,
        "2**x": lambda m, x: 2 ** x,
        "tan(x) sinh(x)": lambda m, x: m.tan(x) * m.sinh(x),
        "x**-0.5 + tan x": lambda m, x: x ** -0.5 + m.tan(x),
        "(1 + x**2)**0.5": lambda m, x: (1 + x ** 2) ** 0.5,
    }
    X = sp.symbols("x")
    for name, F in FUNCS.items():
        dexpr = sp.diff(F(sp, X), X)
        errors = (abs(derivative(lambda t: F(dmath, t), x0)[1] - float(dexpr.subs(X, x0))) for x0 in (0.3, 0.6, 0.9, 1.2, 1.4))
        worst = max(errors)
        assert worst < 1e-12, (name, worst)
        print(f"{name:18s} max |Dual - SymPy| = {worst:.1e}")
    try:
        Dual(-8.0, 1.0) ** (1 / 3)
    except ValueError as err:
        print("ValueError:", err)
    print(type((-8.0) ** (1 / 3)).__name__, (-8.0) ** (1 / 3))
finally:
    Dual.__pow__ = _int_pow
    del Dual.__rpow__
```
```text
x**1.5             max |Dual - SymPy| = 0.0e+00
x**x               max |Dual - SymPy| = 0.0e+00
2**x               max |Dual - SymPy| = 2.2e-16
tan(x) sinh(x)     max |Dual - SymPy| = 8.9e-16
x**-0.5 + tan x    max |Dual - SymPy| = 0.0e+00
(1 + x**2)**0.5    max |Dual - SymPy| = 0.0e+00
ValueError: a real power needs a positive base
complex (1.0000000000000002+1.7320508075688772j)
```

For a negative base, the built-in float power does not fail: `(-8.0) ** (1/3)` silently returns a *complex* number, and there is no real derivative to carry. Raising `ValueError`, as above, is the safer choice. The tangent has poles, so the test points stay below $\pi/2$.

**B4.** (a) With $\varepsilon_1^2 = \varepsilon_2^2 = 0$ and $\varepsilon_1\varepsilon_2 \ne 0$ (a *hyper-dual* number), Taylor gives $f(x + \varepsilon_1 + \varepsilon_2) = f + f'(\varepsilon_1 + \varepsilon_2) + \tfrac12 f''(\varepsilon_1 + \varepsilon_2)^2 = f + f'\varepsilon_1 + f'\varepsilon_2 + f''\varepsilon_1\varepsilon_2$, since $(\varepsilon_1 + \varepsilon_2)^2 = 2\varepsilon_1\varepsilon_2$. The seed `Dual(Dual(x, 1), Dual(1, 0))` has inner $\varepsilon_1$ and outer $\varepsilon_2$, so the result is `Dual(Dual(f, f'), Dual(f', f''))` and `y.b.b` holds $f''$. (c) For two variables, seed $x = (x_0 + \varepsilon_1) + 0\varepsilon_2$ and $y = y_0 + 1\varepsilon_2$; the $\varepsilon_1\varepsilon_2$ coefficient `out.b.b` is $\partial^2 f/\partial x\,\partial y$.

```python
def second_derivative(f, x):
    y = f(Dual(Dual(x, 1.0), Dual(1.0, 0.0)))
    return y.b.b


X = sp.symbols("x")
for name, F in (("x exp(x)", lambda m, x: x * m.exp(x)), ("sin(x^2)", lambda m, x: m.sin(x ** 2)),
                ("1/(1+x^2)", lambda m, x: 1 / (1 + x ** 2))):
    exact = float(sp.diff(F(sp, X), X, 2).subs(X, 0.7))
    got = second_derivative(lambda t: F(dmath, t), 0.7)
    print(f"{name:10s} f'' = {got:.15f}  SymPy {exact:.15f}")
    assert abs(got - exact) < 1e-12

x0, y0 = 0.8, 1.3
f2 = lambda x, y: x ** 2 * y + dmath.sin(x * y)
out = f2(Dual(Dual(x0, 1.0), Dual(0.0, 0.0)), Dual(Dual(y0, 0.0), Dual(1.0, 0.0)))
exact_xy = 2 * x0 + math.cos(x0 * y0) - x0 * y0 * math.sin(x0 * y0)
print(out.b.b, exact_xy, out.a.b, out.b.a)
assert abs(out.b.b - exact_xy) < 1e-12 and abs(out.a.b - (2 * x0 * y0 + y0 * math.cos(x0 * y0))) < 1e-12
```
```text
x exp(x)   f'' = 5.437132310170287  SymPy 5.437132310170286
sin(x^2)   f'' = 0.842238976404774  SymPy 0.842238976404774
1/(1+x^2)  f'' = 0.284163994063995  SymPy 0.284163994063995
1.2093198608997064 1.2093198608997067 2.7380863344026123 1.0449762057862229
```

The other entries come free: `out.a.b` is $\partial f/\partial x$ and `out.b.a` is $\partial f/\partial y$. Each second derivative costs one pass with a hyper-dual seed; the whole Hessian of an $n$-variable function costs about $n^2/2$ passes this way, which is why large models use other schemes (ML-2).

**B5.** The determinant is $\det A = \pm\prod_i u_{ii}$, with the sign of the row permutation $P$ ($PA = LU$ gives $\det P\cdot\det A = \det L\cdot\det U = \prod u_{ii}$, and $\det P = \pm1$ is $(-1)^{\text{inversions}}$). The inverse solves $A\mathbf x_j = \mathbf e_j$ for each column, factoring once. For `Vec @ Matrix`, note what must hold: `Vec.__matmul__` returns `NotImplemented` for a non-`Vec`, so Python reaches `Matrix.__rmatmul__`. (Had `__matmul__` been a bare alias of `dot`, it would have raised `AttributeError` on `other._c` and the reflected method would never run: the silent half-working type of §8.)

```python
class Mat(Matrix):
    def _substitute(self, perm, L, U, b):
        n = len(perm)
        y = [0.0] * n
        for i in range(n):
            y[i] = float(b[perm[i]]) - sum(L._rows[i][j] * y[j] for j in range(i))
        x = [0.0] * n
        for i in reversed(range(n)):
            x[i] = (y[i] - sum(U._rows[i][j] * x[j] for j in range(i + 1, n))) / U._rows[i][i]
        return x

    def det(self):
        try:
            perm, L, U = self.lu()
        except ValueError:                                    # exactly singular
            return 0.0
        n = len(perm)
        inversions = sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        return (-1) ** inversions * math.prod(U._rows[i][i] for i in range(n))

    def inverse(self):
        perm, L, U = self.lu()                                # raises ValueError if singular
        n = len(perm)
        cols = [self._substitute(perm, L, U, [float(i == j) for i in range(n)]) for j in range(n)]
        return Mat(zip(*cols))

    def __mul__(self, k):
        if not isinstance(k, numbers.Real):
            return NotImplemented
        return Mat([[k * a for a in row] for row in self._rows])

    __rmul__ = __mul__

    def __rmatmul__(self, other):                             # row vector times matrix
        if not isinstance(other, Vec):
            return NotImplemented
        return Vec(*(sum(a * b for a, b in zip(other, col)) for col in zip(*self._rows)))


rng = np.random.default_rng(1)
for _ in range(100):
    M, row = rng.normal(size=(5, 5)), rng.normal(size=5)
    A = Mat(M)
    assert math.isclose(A.det(), np.linalg.det(M), rel_tol=1e-9, abs_tol=1e-12)
    assert np.allclose(np.array(A.inverse()._rows), np.linalg.inv(M), atol=1e-8)
    assert np.allclose(list(Vec(*row) @ A), row @ M) and np.allclose(np.array((3 * A)._rows), 3 * M)
    assert np.allclose(np.array((A * 0.5)._rows), 0.5 * M)
S = Mat([[1, 2], [2, 4]])
print(S.det(), Mat([[0, 1], [1, 0]]).det(), Mat([[2, 0], [0, 3]]).det())
try:
    S.inverse()
except ValueError as err:
    print("ValueError:", err)
```
```text
0.0 -1.0 6.0
ValueError: matrix is singular
```

The swap matrix has $\det = -1$: one inversion, and $U$ has unit diagonal. Computing the inverse costs one factorization ($\tfrac23 n^3$ flops, LA-2 (3.7)) plus $n$ pairs of substitutions; solving $A\mathbf x = \mathbf b$ for a single $\mathbf b$ rarely needs the inverse at all.

**B6.** Addition pads the shorter coefficient list; multiplication is the convolution $c_k = \sum_{i+j=k}a_ib_j$; `from_roots` multiplies factors $(x - r)$. Subclassing a dataclass is fine because the generated `__init__` and `__post_init__` are inherited. Horner's loop only uses `*` and `+`, so with the reflected scalar operators in place it also works when `x` is a polynomial, and `p(q)` is composition.

```python
from itertools import zip_longest


class Poly(Polynomial):
    @staticmethod
    def _lift(x):
        if isinstance(x, Poly):
            return x
        return Poly((x,)) if isinstance(x, numbers.Real) else None

    def __add__(self, other):
        o = self._lift(other)
        if o is None:
            return NotImplemented
        return Poly(tuple(a + b for a, b in zip_longest(self.coeffs, o.coeffs, fillvalue=0.0)))

    __radd__ = __add__

    def __neg__(self):
        return Poly(tuple(-a for a in self.coeffs))

    def __sub__(self, other):
        o = self._lift(other)
        return NotImplemented if o is None else self + (-o)

    def __rsub__(self, other):
        o = self._lift(other)
        return NotImplemented if o is None else o + (-self)

    def __mul__(self, other):
        o = self._lift(other)
        if o is None:
            return NotImplemented
        out = [0.0] * (len(self.coeffs) + len(o.coeffs) - 1)
        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(o.coeffs):
                out[i + j] += a * b
        return Poly(tuple(out))

    __rmul__ = __mul__

    @classmethod
    def from_roots(cls, *roots):
        result = cls((1.0,))
        for r in roots:
            result = result * cls((-r, 1.0))
        return result


from numpy.polynomial import Polynomial as NP

rng = np.random.default_rng(2)
for _ in range(50):
    a, b = rng.normal(size=int(rng.integers(2, 6))), rng.normal(size=int(rng.integers(2, 6)))
    p, q = Poly(tuple(a)), Poly(tuple(b))
    assert np.allclose((p + q).coeffs, (NP(a) + NP(b)).coef) and np.allclose((p * q).coeffs, (NP(a) * NP(b)).coef)
    assert np.allclose((p - q).coeffs, (NP(a) - NP(b)).coef) and np.allclose((2 * p + 1).coeffs, (2 * NP(a) + 1).coef)
    assert np.allclose(p(q).coeffs, NP(a)(NP(b)).coef)                     # composition by duck typing
r = Poly.from_roots(1.0, 2.0, 3.0)
print(r.coeffs, [r(x) for x in (1.0, 2.0, 3.0)])
assert np.allclose(r.coeffs, NP.fromroots([1, 2, 3]).coef)
```
```text
(-6.0, 11.0, -6.0, 1.0) [0.0, 0.0, 0.0]
```

`p(q)` needed no code of its own: Horner evaluates `acc * x + c`, and with `x = q` the products and sums are polynomial products and sums. The same unedited `__call__` evaluated `Dual`s in §3.9. That is duck typing at work: the algorithm is written against an interface (`*` and `+`), not a type.

**B7.** *Prediction.* With $\epsilon$ the machine epsilon, the rounding error of $f(x+h) - f(x)$ is about $\epsilon|f|$, so the forward error is $E(h) \approx \tfrac12|f''|h + \epsilon|f|/h$. Setting $E'(h) = 0$: $\tfrac12|f''| = \epsilon|f|/h^2$, so $h_{\text{fwd}} = \sqrt{2\epsilon|f|/|f''|}$ with $E \approx \sqrt{2\epsilon|f||f''|}$. For the central difference, $E(h) \approx \tfrac16|f'''|h^2 + \epsilon|f|/h$, so $\tfrac13|f'''|h = \epsilon|f|/h^2$ and $h_{\text{cen}} = (3\epsilon|f|/|f'''|)^{1/3}$. (The central quotient divides by $2h$, so its rounding term is nearer $\epsilon|f|/(2h)$; keeping the cruder $\epsilon|f|/h$ moves $h_{\text{cen}}$ by $2^{1/3} \approx 1.26$, below what this model can resolve.) Take $F(x) = e^{-x}\arctan x$ at $x_0 = 1$.

```python
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

F11 = lambda m, x: m.exp(-x) * m.atan(x)
X = sp.symbols("x")
expr = F11(sp, X)
x0 = 1.0
f0, f1, f2, f3 = (float(sp.diff(expr, X, k).subs(X, x0)) for k in range(4))
f_num = lambda t: F11(math, t)
val, der = derivative(lambda t: F11(dmath, t), x0)
exact = float(sp.diff(expr, X).subs(X, sp.Rational(x0)).evalf(40))
assert abs(val - f0) < 1e-14 and abs(der - exact) < 1e-15
eps = np.finfo(float).eps
h_fwd, h_cen = math.sqrt(2 * eps * abs(f0) / abs(f2)), (3 * eps * abs(f0) / abs(f3)) ** (1 / 3)

hs = np.logspace(-1, -14, 300)
fwd = np.array([abs((f_num(x0 + h) - f_num(x0)) / h - exact) for h in hs])
cen = np.array([abs((f_num(x0 + h) - f_num(x0 - h)) / (2 * h) - exact) for h in hs])
E_fwd, E_cen = math.sqrt(2 * eps * abs(f0) * abs(f2)), 1.5 * eps * abs(f0) / h_cen      # E(h*) in both cases
typ = lambda e, h: float(np.median(e[(hs > h / 3) & (hs < 3 * h)]))          # typical error within a factor 3 of h*
print(f"predicted h*:     forward {h_fwd:.1e}, central {h_cen:.1e}")
print(f"measured argmin:  forward {hs[fwd.argmin()]:.1e}, central {hs[cen.argmin()]:.1e}")
print(f"predicted floor:  forward {E_fwd:.1e}, central {E_cen:.1e}")
print(f"typical error near h*: forward {typ(fwd, h_fwd):.1e}, central {typ(cen, h_cen):.1e}")
print(f"best grid point:  forward {fwd.min():.1e}, central {cen.min():.1e}; Dual {abs(der - exact):.1e}")
assert 0.1 < typ(fwd, h_fwd) / E_fwd < 10 and 0.1 < typ(cen, h_cen) / E_cen < 10 and abs(der - exact) < 1e-15
assert 0.05 < hs[fwd.argmin()] / h_fwd < 20 and 0.05 < hs[cen.argmin()] / h_cen < 20

fig, ax = plt.subplots(figsize=(6.4, 4))
ax.loglog(hs, fwd, label="forward")
ax.loglog(hs, cen, label="central")
for h, c in ((h_fwd, "C0"), (h_cen, "C1")):
    ax.axvline(h, color=c, ls=":")
ax.axhline(max(abs(der - exact), 1e-18), color="k", ls="--", label="Dual")
ax.set(xlabel="h", ylabel="|error|", title="exp(-x) atan(x) at x = 1: predicted minima dotted")
ax.legend()
fig.tight_layout()
fig.savefig(Path("figures") / "b7_error_vs_h.png", dpi=120)
plt.close(fig)
```
```text
predicted h*:     forward 2.2e-08, central 5.8e-06
measured argmin:  forward 5.0e-09, central 2.5e-06
predicted floor:  forward 5.8e-09, central 1.7e-11
typical error near h*: forward 3.9e-09, central 7.2e-12
best grid point:  forward 1.3e-10, central 1.8e-13; Dual 0.0e+00
```

The typical error near the predicted step agrees with the predicted floor to within a factor of a few, and the measured argmin lands within an order of magnitude of $h^*$. The single best grid point can undershoot the prediction by a large factor: the rounding error is a different pseudo-random number at every $h$, and with 300 tries some cancel by luck (the same effect gave the $3\times10^{-14}$ central dip in `error_vs_h.png`). Do not read a floor off the luckiest point. The scaling is the lesson: the forward floor grows like $\sqrt\epsilon$, the central floor like $\epsilon^{2/3}$, and the `Dual` error does not depend on any $h$.

## C. Challenge

**C1.** Reverse mode records, for every operation, its parents and the local partial derivatives $\partial\,\text{out}/\partial\,\text{parent}$ (the same numbers forward mode multiplies by). `backward()` then walks the graph from the output to the inputs in reverse topological order, adding `output.grad * local` into each parent's `grad` (the chain rule, summed over every path). The topological order is built with an explicit stack: a recursive version would overflow Python's recursion limit on the long addition chain of the $n = 1000$ case.

```python
import time


class Var:
    def __init__(self, value, parents=()):
        self.value, self.parents, self.grad = value, parents, 0.0

    def __repr__(self):
        return f"Var({self.value!r}, grad={self.grad!r})"

    @staticmethod
    def lift(x):
        return x if isinstance(x, Var) else Var(float(x))

    def __add__(self, o):
        o = Var.lift(o)
        return Var(self.value + o.value, ((self, 1.0), (o, 1.0)))

    __radd__ = __add__

    def __neg__(self):
        return Var(-self.value, ((self, -1.0),))

    def __sub__(self, o):
        return self + (-Var.lift(o))

    def __mul__(self, o):
        o = Var.lift(o)
        return Var(self.value * o.value, ((self, o.value), (o, self.value)))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = Var.lift(o)
        return Var(self.value / o.value, ((self, 1 / o.value), (o, -self.value / o.value ** 2)))

    def backward(self):
        order, seen, stack = [], set(), [(self, False)]
        while stack:                                          # iterative post-order DFS
            node, done = stack.pop()
            if done:
                order.append(node)
            elif id(node) not in seen:
                seen.add(id(node))
                stack.append((node, True))
                stack.extend((p, False) for p, _ in node.parents)
        self.grad = 1.0
        for node in reversed(order):                          # every node after all of its consumers
            for parent, local in node.parents:
                parent.grad += node.grad * local


rmath = type("rmath", (), {})()
rmath.sin = lambda v: Var(math.sin(v.value), ((v, math.cos(v.value)),)) if isinstance(v, Var) else dmath.sin(v)
rmath.exp = lambda v: Var(math.exp(v.value), ((v, math.exp(v.value)),)) if isinstance(v, Var) else dmath.exp(v)

f = lambda x, y, m: x * y + m.sin(x)
x0, y0 = 1.5, 2.0
vx, vy = Var(x0), Var(y0)
out = f(vx, vy, rmath)
out.backward()
fx = derivative(lambda t: f(t, y0, dmath), x0)[1]                  # forward mode: one pass per input
fy = derivative(lambda t: f(x0, t, dmath), y0)[1]
sx, sy = sp.symbols("x y")
gsym = [float(sp.diff(sx * sy + sp.sin(sx), s_).subs({sx: x0, sy: y0})) for s_ in (sx, sy)]     # SymPy, third opinion
print((vx.grad, vy.grad), (fx, fy), (y0 + math.cos(x0), x0))
assert abs(vx.grad - fx) < 1e-14 and abs(vy.grad - fy) < 1e-14 and abs(vx.grad - (y0 + math.cos(x0))) < 1e-14
assert abs(vx.grad - gsym[0]) < 1e-14 and abs(vy.grad - gsym[1]) < 1e-14

n = 1000
xs = np.linspace(0.1, 3.0, n)
t0 = time.perf_counter()
vs = [Var(float(v)) for v in xs]
total = Var(0.0)
for v in vs:
    total = total + v * v * rmath.sin(v)
total.backward()
t_rev = time.perf_counter() - t0
t0 = time.perf_counter()
grad_fwd = np.empty(n)
for i in range(n):                                                # n forward passes
    seeded = [Dual(float(v), 1.0 if j == i else 0.0) for j, v in enumerate(xs)]
    acc = Dual(0.0)
    for d in seeded:
        acc = acc + d * d * dmath.sin(d)
    grad_fwd[i] = acc.b
t_fwd = time.perf_counter() - t0
exact = 2 * xs * np.sin(xs) + xs ** 2 * np.cos(xs)
assert np.allclose([v.grad for v in vs], exact) and np.allclose(grad_fwd, exact)
print(f"n = {n}: reverse {t_rev * 1e3:.0f} ms, forward {t_fwd * 1e3:.0f} ms, ratio {t_fwd / t_rev:.0f}")
assert t_fwd / t_rev > 20
```
```text
(2.070737201667703, 1.5) (2.070737201667703, 1.5) (2.070737201667703, 1.5)
n = 1000: reverse 11 ms, forward 2959 ms, ratio 279
```

Reverse mode recorded the computation once and produced all 1000 partial derivatives in a single backward sweep, at a cost proportional to the cost of $f$; forward mode repeated the whole evaluation $n$ times, and the measured ratio is of the order of $n$ (below $n$ here because building and sweeping the graph costs more per operation than a `Dual` step). That is the cost model of §3.7, and it is why ML-7 backpropagates: a loss has one output and millions of inputs. Reverse mode's price is memory, since the whole graph must be kept until `backward()`.

**C2.** Both designs assume the integrator evaluates $f$ on a uniform grid and then combines the samples. Gauss–Legendre evaluates at the roots of a Legendre polynomial, in the interior and unevenly spaced, so it does not fit `_combine(ys, h)` or `rule(ys, h)`: in the class design the template method `integrate` must be overridden entirely (the base class stops being a template and `Integrator` only supplies `n`), and in the function design the rule needs `f`, not samples. The redesign makes the rule own the sampling: `rule(f, a, b, n)` returns the integral. Counting evaluations uses a callable object, the use of `__call__` from §3.4.

```python
from numpy.polynomial.legendre import leggauss


class Counted:
    def __init__(self, f):
        self.f, self.calls = f, 0

    def __repr__(self):
        return f"Counted({self.f!r}, calls={self.calls})"

    def __call__(self, x):
        self.calls += 1
        return self.f(x)


def trapezoid(f, a, b, n):
    h = (b - a) / n
    return h * (0.5 * f(a) + sum(f(a + i * h) for i in range(1, n)) + 0.5 * f(b))


def simpson(f, a, b, n):
    h = (b - a) / n
    odd = sum(f(a + i * h) for i in range(1, n, 2))
    even = sum(f(a + i * h) for i in range(2, n, 2))
    return h / 3 * (f(a) + f(b) + 4 * odd + 2 * even)


def gauss_legendre(f, a, b, n):
    t, w = leggauss(n)
    return (b - a) / 2 * sum(wi * f((b - a) / 2 * ti + (a + b) / 2) for ti, wi in zip(t, w))


class Integrator:                       # composition: an Integrator HAS a rule
    def __init__(self, rule, n):
        self.rule, self.n = rule, n

    def integrate(self, f, a, b):
        g = Counted(f)
        return self.rule(g, a, b, self.n), g.calls


problems = {"sin on [0, pi]": (math.sin, 0.0, math.pi, 2.0), "sqrt on [0, 1]": (math.sqrt, 0.0, 1.0, 2 / 3)}
for label, (f, a, b, exact) in problems.items():
    print(label)
    for name, rule, ns in (("trapezoid", trapezoid, (8, 16, 32, 64)), ("simpson", simpson, (8, 16, 32, 64)),
                           ("gauss-legendre", gauss_legendre, (2, 4, 8, 16))):
        rows = [(Integrator(rule, n).integrate(f, a, b), n) for n in ns]
        errs = [abs(v - exact) for (v, _), _ in rows]
        evals = [calls for (_, calls), _ in rows]
        slope = np.polyfit(np.log(evals), np.log(errs), 1)[0] if errs[-1] > 1e-13 else float("nan")
        print(f"  {name:15s} evals {evals} errors {[f'{e:.1e}' for e in errs]}  slope vs evals {slope:.1f}")
        if label.startswith("sqrt") and name != "gauss-legendre":
            assert -1.9 < slope < -1.1                                  # about -1.5: both degrade to O(h^1.5)
    assert errs[-1] < (1e-8 if label.startswith("sin") else 1e-4)               # Gauss-Legendre with 16 nodes
```
```text
sin on [0, pi]
  trapezoid       evals [9, 17, 33, 65] errors ['2.6e-02', '6.4e-03', '1.6e-03', '4.0e-04']  slope vs evals -2.1
  simpson         evals [9, 17, 33, 65] errors ['2.7e-04', '1.7e-05', '1.0e-06', '6.5e-08']  slope vs evals -4.2
  gauss-legendre  evals [2, 4, 8, 16] errors ['6.4e-02', '1.6e-05', '6.7e-15', '8.9e-16']  slope vs evals nan
sqrt on [0, 1]
  trapezoid       evals [9, 17, 33, 65] errors ['8.5e-03', '3.1e-03', '1.1e-03', '4.0e-04']  slope vs evals -1.6
  simpson         evals [9, 17, 33, 65] errors ['3.6e-03', '1.3e-03', '4.5e-04', '1.6e-04']  slope vs evals -1.6
  gauss-legendre  evals [2, 4, 8, 16] errors ['7.2e-03', '1.2e-03', '1.7e-04', '2.3e-05']  slope vs evals -2.8
```

On the smooth integrand the three rules separate by orders of magnitude: error falls like $h^2$ (trapezoid), $h^4$ (Simpson) and, for Gauss–Legendre, faster than any power (sixteen evaluations already reach rounding level; the slope is printed as `nan` once the error is at rounding level). On $\int_0^1\sqrt x\,dx$ the integrand has an infinite derivative at $0$, which breaks the Taylor expansion the error formulas rest on: both uniform-grid rules fall to the same $O(h^{3/2})$ (Simpson loses its advantage), and Gauss–Legendre drops from exponential to algebraic convergence, though it stays ahead per evaluation.

*Which design to keep.* I would keep the composition: an `Integrator` that holds a rule function. The first design's base class encoded an assumption (a uniform grid) that the third rule violated, so adding Gauss–Legendre meant rewriting the base class and every subclass; with `rule(f, a, b, n)` a new rule is one function and existing rules are untouched. What the class still earns is the state around the rule (here the call counter and `n`; in an adaptive integrator, a tolerance and a maximum depth). An abstract base class would add a type to check and nothing else, and a `typing.Protocol` for "a callable `(f, a, b, n) -> float`" documents the interface without any inheritance.

**C3.**

```python
import time

def best_of(fn, repeat=3):
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def matmul_map(A, B):                                        # sum(map(mul, row, col))
    cols = list(zip(*B))
    return [[sum(map(operator.mul, row, col)) for col in cols] for row in A]


def matmul_sumprod(A, B):                                    # math.sumprod runs the loop in C (Python >= 3.12)
    cols = list(zip(*B))
    return [[math.sumprod(row, col) for col in cols] for row in A]


def matmul_axpy(A, B):                                       # row update: no transpose needed
    out = []
    for row in A:
        acc = [0.0] * len(B[0])
        for a, brow in zip(row, B):
            acc = [x + a * y for x, y in zip(acc, brow)]
        out.append(acc)
    return out


rng = np.random.default_rng(3)
print(f"{'n':>4} {'Matrix @':>10} {'map':>8} {'sumprod':>8} {'axpy':>8} {'NumPy':>10}   ms")
sizes, base = (25, 50, 100, 200), []
for n in sizes:
    P, R = rng.normal(size=(n, n)), rng.normal(size=(n, n))
    MA, MB, LA, LB = Matrix(P), Matrix(R), P.tolist(), R.tolist()
    t = [best_of(lambda: MA @ MB), best_of(lambda: matmul_map(LA, LB)), best_of(lambda: matmul_sumprod(LA, LB)),
         best_of(lambda: matmul_axpy(LA, LB)), best_of(lambda: P @ R, repeat=20)]
    assert all(np.allclose(np.array(f(LA, LB)), P @ R) for f in (matmul_map, matmul_sumprod, matmul_axpy))
    base.append(t[0])
    print(f"{n:>4} " + " ".join(f"{x * 1e3:>9.2f}" for x in t))
slope = np.polyfit(np.log(sizes), np.log(base), 1)[0]
print(f"log-log slope of Matrix @: {slope:.2f}")
print(f"at n = 200, relative to Matrix @: map {base[-1] / t[1]:.1f}x, sumprod {base[-1] / t[2]:.1f}x, "
      f"axpy {base[-1] / t[3]:.1f}x, NumPy {base[-1] / t[4]:.0f}x faster")
assert 2.2 < slope < 3.6 and t[4] < t[2] < t[0]
```
```text
   n   Matrix @      map  sumprod     axpy      NumPy   ms
  25      1.20      0.66      0.28      0.94      0.00
  50      8.37      4.87      1.96      6.89      0.01
 100     65.72     35.74     15.26     53.98      0.04
 200    559.29    299.93    167.01    458.11      0.17
log-log slope of Matrix @: 2.96
at n = 200, relative to Matrix @: map 1.9x, sumprod 3.3x, axpy 1.2x, NumPy 3293x faster
```

The exponent is about 3: the algorithm does $2n^3$ flops and the constant per flop is large and constant. The three rewrites trim the interpreter's per-element work: `map(operator.mul, ...)` drops the generator and tuple unpacking of `zip`; `math.sumprod` moves the whole dot product into C but still reads boxed Python floats from lists; the axpy form avoids the transpose but allocates a list per update. Each is a small constant factor; none changes what the interpreter must do for every scalar: look up two objects, check their types, dispatch to `float.__mul__`, allocate a result, and update a reference count. NumPy stores unboxed 8-byte doubles contiguously and hands the whole product to a blocked BLAS routine that uses SIMD instructions and caches, so the gap printed above is more than three orders of magnitude at $n = 200$, and it widens with $n$ because BLAS reuses data in cache. PY-6 builds the array model that makes this possible, and PY-10 measures where the remaining time goes.
